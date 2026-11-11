"""
A small fully-convolutional emitter localizer for single-molecule localization frames, in the style of Deep-STORM and DECODE
==========================================================================================================================

A learned baseline for the 64 x 64 EM-CCD frames of `experiments/smlm_common.py`: a CPU-scale network (<= 300k parameters, minutes of training
on two threads) that turns one frame, in photons, into a per-pixel emitter probability, a sub-pixel offset and a photon count, and a decoder that
turns those maps into detections in the harness's convention. It is the counterpart of the sliding inverse scale space for the question "what does
a learned localizer do on these frames, in distribution and under a mismatch of the assumed PSF and camera" (`experiments/exp_smlm_learned.py`).
Torch is needed only for the network, the loss, the training loop and the weights (`require_torch`); the normalisation, the target encoding, the
decoder and the simulator are numpy.

Conventions (the places where a comparison can silently go wrong)
-----------------------------------------------------------------
* Positions inside this module are in PIXEL-INDEX coordinates (the centre of the pixel with 0-based row i, column j is at (x, y) = (j, i); x is the
  column axis), the convention of `src/core/imaging.py`. The decoder returns NANOMETRES in the corner convention of the ground truth,
  x_nm = (x_px + 0.5) * 100 (`smlm_common.pixels_to_nm`): (0, 0) is the outer corner of the top-left pixel, the field is [0, 6400] nm.
* The output grid is the pixel grid (native resolution, as DECODE): the emitter with centre (x, y) belongs to the pixel (round(y), round(x));
  channel 0 is the logit of "an emitter centre lies in this pixel", channels 1, 2 the offsets (dx, dy) = (x - j, y - i) in pixels, channel 3 the
  natural log of the emitter's photons / 1000.
* Input: `normalize` = (frame - background surface - median of that difference) / s, with s the noise sd of a background pixel (`noise_scale`:
  sqrt(a * mean background + b0), a and the surface from the label-free calibration of `exp_smlm_iss.calibrate`). The per-frame median removes the
  pedestal that the emitters' halos add (the sliding ISS does the same before its high-pass), so a wrong background LEVEL is not seen by the
  network; a wrong noise level (s) is. The network feeds [x / 8, asinh(x / 2)] to its first layer (a linear and a compressed channel: peaks reach
  100 sd, the faintest detectable emitters are a few sd).
* The target of training is the pixel that contains the emitter (logit) and, on the 3 x 3 neighbourhood of that pixel, the offset of the emitter from
  each neighbour's centre ("voting": the offsets are in [-1.5, 1.5]) and the log photons; the decoder reads a detection at a local maximum of the probability map: its score is the pixel's probability and its position the
  vote of that pixel (the options of summing the probability over the 3 x 3 window and of pooling the votes are available; neither helped the efficiency on development frames,
  see exp_smlm_learned.py).

Architecture (`LocalizerNet`): space-to-depth U-Net. The 64 x 64 input is unshuffled to 32 x 32 x 8 (lossless), two stride-2 stages and a
two-level decoder with skips (3 x 3 convolutions, ReLU, no normalisation layers), a 1 x 1 head with 16 channels shuffled back to 64 x 64 x 4. Fully
convolutional: any frame size divisible by 8 works. The default widths (32, 64, 64) give 289k parameters.

Loss (`loss_terms`): per-frame sum of the binary cross-entropy of the probability map + w_off * Huber(offsets, 0.3 px) on the voting mask + w_ph * Huber(log photons)
on the same mask, averaged over the frames of the batch. Decoder (`decode`): see there.
"""
import io
import json
import math
import time

import numpy as np

try:                                                              # the numpy part works without torch
    import torch
    from torch import nn
    import torch.nn.functional as F
except Exception:                                                 # pragma: no cover
    torch = nn = F = None

PIXEL_NM = 100.0
B0 = 0.078                           # photon^2: nominal read-out term of the variance (as exp_smlm_iss.B0)
DEFAULT_ARCH = dict(c1=32, c2=64, c3=64)
LOGPH_SCALE = 1000.0                 # photons are regressed as log(photons / LOGPH_SCALE)
PRIOR_LOGIT = -5.5                   # initial bias of the probability logit (about 0.4 % per pixel)


def require_torch():
    if torch is None:
        raise ImportError('torch is needed for the network, the loss and the training loop of learned_localizer (not for encode / decode / simulate)')


# ----------------------------------------------------------------------------------------------------------------------------------
# input normalisation
# ----------------------------------------------------------------------------------------------------------------------------------
def noise_scale(bg_mean, var_slope, b0=B0):
    """The noise sd of a background pixel in photons: sqrt(a * bg + b0) for the variance law var = a x + b0 (x in photons)."""
    return float(math.sqrt(var_slope * bg_mean + b0))


def normalize(frames, bg, scale):
    """(frames - bg - per-frame median of that difference) / scale as float32; frames (n, H, W) or (H, W) in photons, bg a scalar or an (H, W) surface."""
    f = np.asarray(frames, dtype=np.float64)
    single = f.ndim == 2
    f = f[None] if single else f
    f = f - np.asarray(bg, dtype=np.float64)
    f = f - np.median(f.reshape(len(f), -1), axis=1)[:, None, None]
    out = (f / float(scale)).astype(np.float32)
    return out[0] if single else out


# ----------------------------------------------------------------------------------------------------------------------------------
# targets and augmentation (numpy)
# ----------------------------------------------------------------------------------------------------------------------------------
def encode_targets(emitters, shape=(64, 64), vote=1):
    """Training targets for a batch from the emitters of each frame.

    emitters: list over the B frames of arrays (k, 3) = x_px, y_px (pixel-index coordinates), photons. Emitters whose centre lies outside the field
    (round(x) or round(y) outside the grid) are skipped. Returns dict of float32 arrays: prob (B, H, W) (1 at the pixel containing a centre), off (B, 2, H, W)
    (dx, dy of the nearest emitter relative to each pixel centre of its (2 vote + 1)^2 neighbourhood), logph (B, H, W) (log(photons / 1000) of that emitter) and
    mask (B, H, W) (1 on the voting neighbourhoods; where neighbourhoods of two emitters overlap the nearer emitter wins). n_emit (B,) counts the encoded emitters."""
    H, W = shape
    B = len(emitters)
    prob = np.zeros((B, H, W), np.float32)
    off = np.zeros((B, 2, H, W), np.float32)
    logph = np.zeros((B, H, W), np.float32)
    mask = np.zeros((B, H, W), np.float32)
    dist = np.full((B, H, W), np.inf, np.float32)
    n_emit = np.zeros(B, np.int64)
    for b, e in enumerate(emitters):
        e = np.asarray(e, dtype=np.float64).reshape(-1, 3)
        for x, y, ph in e:
            ix, iy = int(math.floor(x + 0.5)), int(math.floor(y + 0.5))
            if not (0 <= ix < W and 0 <= iy < H):
                continue
            n_emit[b] += 1
            prob[b, iy, ix] = 1.0
            lp = math.log(max(ph, 1.0) / LOGPH_SCALE)
            for di in range(-vote, vote + 1):
                for dj in range(-vote, vote + 1):
                    i, j = iy + di, ix + dj
                    if 0 <= i < H and 0 <= j < W:
                        d = math.hypot(x - j, y - i)
                        if d < dist[b, i, j]:
                            dist[b, i, j] = d
                            off[b, 0, i, j], off[b, 1, i, j] = x - j, y - i
                            logph[b, i, j] = lp
                            mask[b, i, j] = 1.0
    return dict(prob=prob, off=off, logph=logph, mask=mask, n_emit=n_emit)


def augment_d4(frame, emitters, k):
    """One element k = 0..7 of the symmetry group of the square (rotations by multiples of 90 degrees and flips) applied to a frame (H, W) and its emitters
    (m, 3) in pixel-index coordinates: bit 0 transposes (swaps x and y), bit 1 flips the columns (x -> W - 1 - x), bit 2 flips the rows (y -> H - 1 - y).
    The frame must be square when bit 0 is set. Returns (frame, emitters) as new arrays."""
    H, W = frame.shape
    e = np.array(emitters, dtype=np.float64, copy=True).reshape(-1, 3)
    f = frame
    if k & 1:
        f = f.T
        e[:, [0, 1]] = e[:, [1, 0]]
    if k & 2:
        f = f[:, ::-1]
        e[:, 0] = (f.shape[1] - 1) - e[:, 0]
    if k & 4:
        f = f[::-1, :]
        e[:, 1] = (f.shape[0] - 1) - e[:, 1]
    return np.ascontiguousarray(f), e


def d4_maps(arr, k):
    """The element k of the square's symmetry group (``augment_d4``) applied to the last two axes (rows, columns) of an array of maps (..., H, W)."""
    if k & 1:
        arr = np.swapaxes(arr, -1, -2)
    if k & 2:
        arr = arr[..., ::-1]
    if k & 4:
        arr = arr[..., ::-1, :]
    return arr


def d4_maps_inverse(arr, k):
    """Inverse of ``d4_maps`` (the operations undone in the opposite order)."""
    if k & 4:
        arr = arr[..., ::-1, :]
    if k & 2:
        arr = arr[..., ::-1]
    if k & 1:
        arr = np.swapaxes(arr, -1, -2)
    return arr


def d4_offsets_inverse(off, k):
    """Maps of offset VECTORS (..., 2, H, W) = (dx, dy) predicted in a frame transformed by element k, expressed in the original frame: the spatial inverse of
    ``d4_maps`` and the inverse of the linear part (a row flip changes the sign of dy, a column flip that of dx, the transposition swaps them)."""
    o = d4_maps_inverse(off, k)
    dx, dy = o[..., 0, :, :], o[..., 1, :, :]
    if k & 4:
        dy = -dy
    if k & 2:
        dx = -dx
    if k & 1:
        dx, dy = dy, dx
    return np.stack([dx, dy], axis=-3)


# ----------------------------------------------------------------------------------------------------------------------------------
# decoder (numpy)
# ----------------------------------------------------------------------------------------------------------------------------------
def _box3(a):
    """Sum over the 3 x 3 window of a (H, W) array, zero outside."""
    p = np.pad(a, 1)
    H, W = a.shape
    return sum(p[i:i + H, j:j + W] for i in range(3) for j in range(3))


def decode_frame(prob, off, logph, threshold=0.5, window=0, nms_px=1.5, min_score=0.01, pos_window=0):
    """Detections of one frame from the maps. prob (H, W) in [0, 1], off (2, H, W) pixel offsets (dx, dy), logph (H, W).

    A candidate is a local maximum of the probability map over its 3 x 3 neighbourhood (an exact tie between neighbours goes to the first in raster order). Its score is
    the summed probability of its (2 window + 1)^2 window, capped at 1 (window = 0, the default: the pixel's own probability, as DECODE's per-pixel threshold; window = 1: the
    probability that an emitter lies in the 3 x 3 window, which collects the probability that the network splits over neighbouring pixels for an emitter on a pixel boundary). Its position is the probability-weighted mean of the pixel centres plus offsets over the (2 pos_window + 1)^2 window
    (pos_window = 0: the vote of the pixel itself), its photons exp of the probability-weighted mean log photons times 1000 over the same window. The candidates are
    suppressed greedily in order of decreasing score (a candidate closer than ``nms_px`` pixels to a kept one is dropped) and those with a score below ``threshold`` are
    discarded; as a suppressor always has the higher score, the result at a threshold equals the result at a lower one restricted to the scores above it (a stored list of
    candidates can be thresholded later). Candidates scoring below ``min_score`` are never formed. Returns an (n, 4) array x_nm, y_nm, photons, score (nm in the corner
    convention), by decreasing score."""
    H, W = prob.shape
    pw = np.pad(prob, 1)
    mx = np.max([pw[i:i + H, j:j + W] for i in range(3) for j in range(3)], axis=0)
    pwin = np.pad(prob, window)
    S = prob.copy() if window == 0 else sum(pwin[i:i + H, j:j + W] for i in range(2 * window + 1) for j in range(2 * window + 1))
    idx = np.argwhere((prob >= mx) & (prob > 0.0) & (np.minimum(S, 1.0) >= min_score))
    if len(idx) == 0:
        return np.zeros((0, 4))
    jj, ii = np.meshgrid(np.arange(W), np.arange(H))
    X, Y = jj + off[0], ii + off[1]
    pp = pos_window
    PP, PX, PY, PL = np.pad(prob, pp), np.pad(prob * X, pp), np.pad(prob * Y, pp), np.pad(prob * logph, pp)
    rows = []
    claimed = np.zeros((H, W), bool)
    for i, j in idx:                                                    # raster order: equal neighbouring maxima (ties) are resolved to the first
        if claimed[max(i - 1, 0):i + 2, max(j - 1, 0):j + 2].any():
            continue
        claimed[i, j] = True
        sl = (slice(i, i + 2 * pp + 1), slice(j, j + 2 * pp + 1))
        w = PP[sl].sum()
        rows.append(((PX[sl].sum() / w + 0.5) * PIXEL_NM, (PY[sl].sum() / w + 0.5) * PIXEL_NM, LOGPH_SCALE * math.exp(PL[sl].sum() / w), min(float(S[i, j]), 1.0)))
    out = np.array(rows).reshape(-1, 4)
    out = out[np.argsort(-out[:, 3], kind='stable')]
    keep = []
    for r in range(len(out)):
        if all(math.hypot(out[r, 0] - out[q, 0], out[r, 1] - out[q, 1]) >= nms_px * PIXEL_NM for q in keep):
            keep.append(r)
    out = out[keep]
    return out[out[:, 3] >= threshold]


def decode(prob, off, logph, threshold=0.5, window=0, nms_px=1.5, min_score=0.01, pos_window=0):
    """``decode_frame`` over a batch: prob (B, H, W), off (B, 2, H, W), logph (B, H, W); returns a list of B arrays (n, 4) x_nm, y_nm, photons, score."""
    return [decode_frame(prob[b], off[b], logph[b], threshold, window, nms_px, min_score, pos_window) for b in range(len(prob))]


# ----------------------------------------------------------------------------------------------------------------------------------
# network, loss, training, inference (torch)
# ----------------------------------------------------------------------------------------------------------------------------------
if torch is not None:
    def _block(i, o):
        return nn.Sequential(nn.Conv2d(i, o, 3, padding=1), nn.ReLU(inplace=True), nn.Conv2d(o, o, 3, padding=1), nn.ReLU(inplace=True))

    class LocalizerNet(nn.Module):
        """Space-to-depth U-Net: (B, 1, H, W) normalised frames -> (B, 4, H, W) maps (probability logit, dx, dy, log photons / 1000); H, W divisible by 8."""

        def __init__(self, c1=32, c2=64, c3=64):
            super().__init__()
            self.arch = dict(c1=int(c1), c2=int(c2), c3=int(c3))
            self.stem = _block(8, c1)
            self.down1 = _block(c1, c2)
            self.down2 = _block(c2, c3)
            self.up2 = _block(c3 + c2, c2)
            self.up1 = _block(c2 + c1, c1)
            self.head = nn.Conv2d(c1, 16, 1)
            with torch.no_grad():
                self.head.bias.zero_()
                self.head.bias[:4] = PRIOR_LOGIT                            # channels 0..3 of the unshuffled head are the four phases of the logit
                self.head.weight.mul_(0.1)

        def forward(self, x):
            z = F.pixel_unshuffle(torch.cat([x / 8.0, torch.asinh(x / 2.0)], dim=1), 2)
            a = self.stem(z)
            b = self.down1(F.max_pool2d(a, 2))
            c = self.down2(F.max_pool2d(b, 2))
            d = self.up2(torch.cat([F.interpolate(c, scale_factor=2, mode='nearest'), b], dim=1))
            e = self.up1(torch.cat([F.interpolate(d, scale_factor=2, mode='nearest'), a], dim=1))
            return F.pixel_shuffle(self.head(e), 2)

    def n_params(model):
        return int(sum(p.numel() for p in model.parameters()))

    def make_model(seed=0, **arch):
        """A freshly initialised network (PyTorch default initialisation, seeded) with the architecture ``DEFAULT_ARCH`` updated by ``arch``."""
        torch.manual_seed(int(seed))
        cfg = dict(DEFAULT_ARCH)
        cfg.update(arch)
        return LocalizerNet(**cfg).to(memory_format=torch.channels_last)

    def loss_terms(out, tgt, w_off=2.0, w_ph=0.1, huber=0.3):
        """Loss of a batch: out (B, 4, H, W) network output, tgt dict of tensors prob (B, H, W), off (B, 2, H, W), logph (B, H, W), mask (B, H, W).
        Returns (total, dict of the three terms), each a per-frame sum averaged over the batch."""
        B = out.shape[0]
        bce = F.binary_cross_entropy_with_logits(out[:, 0], tgt['prob'], reduction='sum') / B
        m = tgt['mask']
        lo = (F.huber_loss(out[:, 1:3], tgt['off'], reduction='none', delta=huber) * m[:, None]).sum() / B
        lp = (F.huber_loss(out[:, 3], tgt['logph'], reduction='none', delta=1.0) * m).sum() / B
        return bce + w_off * lo + w_ph * lp, dict(bce=bce.detach(), off=lo.detach(), ph=lp.detach())

    def lr_at(step, steps, lr, warmup=100, end=0.02):
        """Linear warm-up then cosine decay to end * lr (a function of the step only, so a resumed run continues the same schedule)."""
        w = min(1.0, (step + 1) / max(warmup, 1))
        return lr * w * (end + (1.0 - end) * 0.5 * (1.0 + math.cos(math.pi * min(step, steps) / max(steps, 1))))

    def to_tensors(x, tgt=None):
        xt = torch.from_numpy(np.ascontiguousarray(x, dtype=np.float32))[:, None].contiguous(memory_format=torch.channels_last)
        if tgt is None:
            return xt
        return xt, {k: torch.from_numpy(np.ascontiguousarray(tgt[k])) for k in ('prob', 'off', 'logph', 'mask')}

    def train(model, batch_fn, steps, *, start_step=0, opt=None, lr=2e-3, weight_decay=1e-4, warmup=100, clip=2.0, w_off=2.0, w_ph=0.1,
              callback=None, callback_every=0, log=None):
        """Adam(W) training for steps start_step .. steps - 1. ``batch_fn(step)`` returns (x, targets): x (B, H, W) float32 normalised frames and the dict of
        ``encode_targets``; it must be a function of the step only (its own seeded generator), which makes a run reproducible and resumable. ``callback(step_done, model,
        opt, history)`` is called every ``callback_every`` steps and at the end. Returns (opt, history) with history a list of (step, total, bce, off, ph) per step."""
        require_torch()
        opt = opt or torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=weight_decay)
        hist = []
        model.train()
        for step in range(start_step, steps):
            for g in opt.param_groups:
                g['lr'] = lr_at(step, steps, lr, warmup)
            x, tgt = batch_fn(step)
            xt, tt = to_tensors(x, tgt)
            opt.zero_grad(set_to_none=True)
            out = model(xt)
            loss, parts = loss_terms(out, tt, w_off, w_ph)
            loss.backward()
            if clip:
                torch.nn.utils.clip_grad_norm_(model.parameters(), clip)
            opt.step()
            hist.append((step, float(loss.detach()), float(parts['bce']), float(parts['off']), float(parts['ph'])))
            if log and (step + 1) % max(callback_every or 100, 1) == 0:
                log('step %d/%d loss %.4f (bce %.4f, off %.4f, ph %.4f)' % ((step + 1, steps) + hist[-1][1:]))
            if callback and callback_every and ((step + 1) % callback_every == 0 or step + 1 == steps):
                callback(step + 1, model, opt, hist)
                model.train()
        return opt, hist

    def _predict_maps(model, x, batch):
        outs = []
        with torch.no_grad():
            for i in range(0, len(x), batch):
                outs.append(model(to_tensors(x[i:i + batch])).float().numpy())
        o = np.concatenate(outs) if outs else np.zeros((0, 4) + tuple(x.shape[1:]), np.float32)
        return 1.0 / (1.0 + np.exp(-o[:, 0].astype(np.float64))), o[:, 1:3].astype(np.float64), o[:, 3].astype(np.float64)

    def predict(model, x, batch=64, tta=False):
        """Maps of normalised frames x (n, H, W) (square frames when tta): (prob (n, H, W) in [0, 1], off (n, 2, H, W), logph (n, H, W)). With tta the frames are also
        presented in all 8 orientations of the square (``d4_maps``) and the maps, brought back to the original orientation, are averaged (for a network trained with
        the 8-fold augmentation of ``augment_d4``; a network that has learned an orientation-specific PSF offset would be hurt by it)."""
        require_torch()
        model.eval()
        if not tta:
            return _predict_maps(model, x, batch)
        acc = None
        for k in range(8):
            p, o, l = _predict_maps(model, np.ascontiguousarray(d4_maps(x, k)), batch)
            m = (d4_maps_inverse(p, k), d4_offsets_inverse(o, k), d4_maps_inverse(l, k))
            acc = m if acc is None else tuple(a + b for a, b in zip(acc, m))
        return tuple(a / 8.0 for a in acc)

    def detect(model, x, threshold=0.5, window=0, nms_px=1.5, min_score=0.01, batch=64, tta=False, pos_window=0):
        """Normalised frames (n, H, W) -> list of n arrays (m, 4) x_nm, y_nm, photons, score (``predict`` then ``decode``)."""
        p, o, l = predict(model, x, batch, tta)
        return decode(p, o, l, threshold, window, nms_px, min_score, pos_window)

    def save(model, path, meta=None):
        """Weights and the architecture as an .npz (arrays only, no pickle) with ``meta`` (a JSON-serialisable dict) in the entry '_meta'."""
        arrays = {'w|' + k: v.detach().cpu().numpy() for k, v in model.state_dict().items()}
        arrays['_arch'] = np.array(json.dumps(model.arch))
        arrays['_meta'] = np.array(json.dumps(meta or {}, default=lambda o: o.item() if hasattr(o, 'item') else str(o)))
        np.savez(path, **arrays)

    def load(path):
        """Inverse of ``save``: returns (model in eval mode, meta dict)."""
        require_torch()
        with np.load(path, allow_pickle=False) as z:
            arch = json.loads(str(z['_arch']))
            meta = json.loads(str(z['_meta']))
            model = LocalizerNet(**arch)
            model.load_state_dict({k[2:]: torch.from_numpy(z[k]) for k in z.files if k.startswith('w|')})
        return model.to(memory_format=torch.channels_last).eval(), meta
else:                                                                  # pragma: no cover
    LocalizerNet = None


# ----------------------------------------------------------------------------------------------------------------------------------
# simulation of frames (numpy; imaging.py supplies the PSF and the camera)
# ----------------------------------------------------------------------------------------------------------------------------------
def simulate_frames(n, rng, *, bg, density, photons, sigma, parts=((0.75, 1.0), (0.25, 3.0)), psf_scale=1.0, photon_scale=1.0, cam=None, shape=(64, 64),
                    camera='emccd'):
    """n frames of the EM-CCD of the challenge (`imaging.simulate_emccd`) with a random number of emitters each.

    Per frame k ~ Poisson(density) emitters at positions uniform over the field (x, y in [-0.5, W - 0.5) x [-0.5, H - 0.5) pixel-index coordinates), photons from
    ``photons`` (an array: resampled with replacement; or a callable (rng, k) -> (k,) array), widths from ``sigma`` (a float, or a callable (rng, k) -> (k,) core widths in
    px), rendered as the Gaussian mixture ``parts`` = ((fraction, width multiplier), ...) (the PSF of `exp_smlm_iss.CONFIGS['hp5h']`: a core of the emitter's width
    and a 3 x wider halo of a quarter of the photons), each emitter with a width of its own (`PixelPSF` with c = (x, y, sigma)). ``psf_scale`` multiplies every width
    and ``photon_scale`` every photon count (the deliberate mismatches of the sweep), ``bg`` is the background in photons (scalar or (H, W) surface), ``cam`` the
    camera dict of `imaging.EMCCD`. camera 'poisson' replaces the EM-CCD by Poisson shot noise only (a light test setting).
    Returns (frames (n, H, W) float64 in photons including the background, emitters: list of n arrays (k, 3) = x_px, y_px, photons as simulated)."""
    from src.core import imaging as IM
    H, W = shape
    cam = IM.EMCCD if cam is None else cam
    psf = IM.PixelPSF(H, W, 1.0, parts=parts)
    frames = np.empty((n, H, W))
    ems = []
    bgv = np.asarray(bg, dtype=np.float64) * np.ones((H, W))
    for f in range(n):
        k = int(rng.poisson(density))
        x, y = rng.uniform(-0.5, W - 0.5, k), rng.uniform(-0.5, H - 0.5, k)
        ph = (np.asarray(photons(rng, k), float) if callable(photons) else rng.choice(np.asarray(photons, float), size=k)) * photon_scale
        sg = (np.asarray(sigma(rng, k), float) if callable(sigma) else np.full(k, float(sigma))) * psf_scale
        c = np.stack([x, y, sg], axis=1)
        if camera == 'poisson':
            frames[f] = rng.poisson(bgv + psf.synth(c, ph)).astype(float)
        else:
            frames[f] = IM.simulate_emccd(psf, c, ph, bgv, rng, cam)
        ems.append(np.stack([x, y, ph], axis=1))
    return frames, ems
