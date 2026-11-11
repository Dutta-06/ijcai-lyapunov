"""
A CPU-scale learned localizer on the simulated 2D single-molecule benchmarks, against the sliding inverse scale space and Picasso, in distribution and under a mismatch
=====================================================================================================================================================================

The question of a reviewer: "a deep-learning localizer (DECODE, Deep-STORM) is better". DECODE needs a GPU and was never published for the 2D MT0 sets (its challenge numbers are
3D only; the authors' GPU runs are specified in `results/smlm/decode_run_spec.md`). This script trains the small fully-convolutional localizer of `src/core/learned_localizer.py`
(space-to-depth U-Net, <= 300k parameters, per-pixel probability + sub-pixel offsets + photons, decoder with non-maximum suppression) on a CPU, scores it with the harness of
`experiments/smlm_common.py` on the SAME test frames as the sliding ISS (`exp_smlm_iss.py`) and the Picasso ladder (`exp_smlm_picasso.py`), and measures how it degrades when the assumed
PSF and camera are wrong, against the ISS and Picasso under the same mismatch. Be honest: if the CNN wins in distribution it is reported; the questions are the size of the gap, the
training data it needs and the degradation under mismatch.

DESIGN, DECLARED BEFORE ANY TEST-FRAME EVALUATION (nothing below was tuned on test frames; what was chosen where is listed at the end)
----------------------------------------------------------------------------------------------------------------------------------
Protocol: as `exp_smlm_iss.PROTOCOL`. Development frames: the first 1000 (LD) or 100 (HD) frames of a set, with labels; test frames: LD 1000-1499, HD 100-499 (0-based stack indices).
Frames enter the network as (frame - background surface - frame median) / s in units of the background noise sd (`learned_localizer.normalize`), the background surface and the
variance slope a from the label-free calibration of `exp_smlm_iss.calibrate` (the dense sets use the camera of the sparse set of their noise class, as the ISS does).

(a) IN DISTRIBUTION, SUPERVISED (tag `sup`, one network per set). Trained on the labelled development frames of THAT set only: the first 80 % of them (LD 800 frames, HD 80), the last
    20 % (LD 200, HD 20) held out. Augmentation: the 8 flips/rotations of the square (positions transformed with the frames) and a photon rescaling of exp(U(ln 0.9, ln 1.1)) of the
    normalised frame and its labels. The held-out frames choose the checkpoint (every CKPT_EVERY steps), the probability threshold (THRESHOLDS) and a constant position offset (the mean
    signed error of the matched held-out pairs, removed from every detection: the "+c" of the Picasso ladder and of the E-tuned ISS), jointly by the efficiency E of the offset-corrected held-out
    detections. The test frames are evaluated once, at the chosen checkpoint and threshold.

(b) SIMULATION-TRAINED, NO CHALLENGE LABELS ON THE FRAMES (tag `sim`, one network per set). Trained on frames simulated with `imaging.simulate_emccd` (EM gain, read-out noise, quantisation):
    * camera: the calibrated background surface and variance slope of the set (label-free, the ISS's inputs); the physical EM-CCD noise law of the simulator;
    * PSF: the Gaussian mixture of the ISS (`exp_smlm_iss.CONFIGS['hp5h']`: 75 % of the photons in a core, 25 % in a halo of three times the width) at the ISS's entry width
      SIGMA0 = 1.4 px, with a width of its own per emitter (the ISS lets every emitter fit a width in [1.0, 3.2] px): the core width is 1.4 px times a ratio drawn from the
      empirical spread of `results/smlm/calibration.md` (four defocus classes with probabilities 0.34 / 0.38 / 0.155 / 0.124, median ratios 1.277 / 1.393 / 1.656 / 2.358 over 1.369,
      log-normal scatter 0.04 / 0.07 / 0.13 / 0.25), clipped to [1.0, 3.2] px. The real PSF's long halo (44 % of the photons beyond 5 px, cut at 25 px) and its orientation-dependent
      centre offset are NOT in the simulation: this is the realistic simulation-to-real gap;
    * emitters per frame: Poisson with the mean of the set's DEVELOPMENT labels (LD about 1.25, HD about 12.9), positions uniform over the field;
    * photons per emitter: resampled with replacement from the DEVELOPMENT labels' photon counts of the set (the illumination profile of the simulator is not modelled).
    The checkpoint and the threshold are chosen on simulated validation frames only (fresh seed, same simulator, same assumptions), by E; no offset correction (nothing on the real frames
    is used to choose anything). Labels of the development frames enter only through the emitter density, the photon counts and (via `calibration.md`) the width spread.

(c) MISMATCH SWEEP (MT0.N1.LD and MT0.N1.HD). The simulation of (b) is made deliberately wrong and the network retrained with the SAME recipe, seeds and steps (tags): `psf0.8`, `psf1.25`
    (every core width times the factor), `bg0.5`, `bg2` (the simulated background surface times the factor; the network's input is normalised with the believed background, and
    the noise scale s follows the believed background through the physical law, a unchanged: a wrong background level is removed by the frame-median step, the wrong noise level is not),
    `ph0.5` (every simulated photon count halved: lower SNR). The checkpoint and threshold are chosen on the validation frames of the WRONG simulation. They are tested on the true test
    frames (the first SWEEP_FRAMES of the test range: LD 200 frames, 1000-1199; HD 100 frames, 100-199, about 1,260 emitters: the ISS needs ~5 s per HD frame on a free core, several times more when it underestimates the
    noise, and the machine was shared). The sliding ISS gets the same wrong inputs with its a priori default (config `hp5h`, kappa 5, no operating
    point chosen anywhere), run by `exp_smlm_iss.run_frames` on the same frames: its entry width sigma0 times 0.8 / 1.25 (the width bounds [1.0, 3.2] stay), its background surface times
    0.5 / 2 with the variance slope a unchanged (the camera law applied to the believed background: the same belief as the network's simulation; the conditions `bg0.5/a`, `bg2/a` repeat it
    with the slope recomputed for the wrong surface, a / factor, i.e. the noise at a background pixel kept right: the milder reading of "slope consistent"). `ph0.5` has no counterpart: the
    ISS has no photon prior (its entry is the nominal run). Picasso has no PSF, background or photon input to mismatch: its default rung R0 is unchanged; the options of its tuned rung R2E
    include a width window on the fitted sigma, its only quasi-input: stage `picasso` re-runs Picasso with the options of R2E (offset correction included) and the window times 0.8 and 1.25
    ($SMLM_PICASSO_PYTHON must name an interpreter that can import picasso; the factor 1 re-runs the rung itself and must reproduce the stored detections).

(a') DATA SCALING (LD sets, tags `sup100`, `sup200`, `sup400`; `--scaling`). The network of (a) with only the first 100, 200 or 400 of the 800 training frames, the same steps; scored on the same
    held-out development frames only (never on test frames): how much labelled data the supervised network needs.

WHAT WAS CHOSEN WHERE (pilot runs of a few thousand steps on MT0.N1.LD, before any test-frame evaluation; no pilot read a test frame)
    * Network: the space-to-depth U-Net with widths (32, 64, 64), 289,040 parameters (the budget is 300k), fixed a priori by the budget. A pilot compared it with widths (16, 32, 32) (73k parameters):
      after 2000 steps of batch 16 the held-out development E of the supervised network was 57.9 (large) against 56.5 (small), the simulation-trained network's E on the 1000 LD development frames
      51.8 against 54.7 (best threshold, no offset removed): no difference beyond noise; the large one was kept so that capacity cannot be the objection.
    * Steps, batch, learning rate: batch 16, AdamW at 2e-3 with warm-up and cosine decay; the simulated validation E of the pilots had saturated by ~1000 steps (64.8, 66.8, 64.6, 67.6 at 500,
      1000, 1500, 2000 steps), so the simulation-trained networks get 1500 steps; the supervised network 2000 (its held-out E was still rising at 1600-2000).
    * (a) augmentation: without the flips/rotations the supervised network memorised the 800 training frames (training BCE 0.82 against 2.65 with them at 1200 steps; held-out Jaccard index
      falling to 0.63 at the best threshold); the augmentation was kept.
    * Decoder: the pixel's own probability as score and its own vote as position (`DECODE`), suppression within 1.5 px. Summing the probability over the 3 x 3 window and pooling the votes were
      compared on the 1000 LD development frames (simulation-trained pilot network; E at the best threshold with / without the offset removed: 54.4 / 56.9 for the chosen decoder, 51.8 / 54.3 pooled
      votes, 52.3 / 53.6 window score) and on the 200 held-out frames (supervised pilot, offset removed: 59.4 against 57.9): the plain decoder was kept. Test-time augmentation over the 8 orientations
      raised J (0.816 to 0.848) and RMSE (37.8 to 41.1 nm) and lowered E (57.9 to 56.2) on the held-out frames: not used.
    * Threads: one per process (`LEARNED_THREADS`): on the oversubscribed machine two OpenMP threads were up to 30 times slower on the unit tests' kernels and 2-4 times slower in the pilots.
    These choices used labelled development frames (the decoder and augmentation comparisons used real development frames, also for the simulation-trained pilot network); they are choices of the
    pipeline, like the ISS's model configurations. What the simulation-trained networks of (b), (c) choose per model (checkpoint, threshold) uses simulated frames only.

Comparisons and fairness (also in the report): the supervised network is compared with the ISS and the Picasso rows tuned on the same development labels (E-tuned ISS with its offset
removed, Picasso R2E+c); the simulation-trained network, which sees no challenge labels on the frames, with the label-free rows (the ISS's a priori default, Picasso R0). The CNN's tuning
effort is: one architecture, one learning-rate schedule, 11 thresholds and 6-8 checkpoints per network.

Usage (the learned part needs torch: a separate environment, e.g. `python -m venv /path/venv_learn; pip install torch numpy scipy pillow pytest`):
    python experiments/exp_smlm_learned.py --run [--stage prep|train|infer|iss|picasso|collect] [--keys 'MT0.N1.LD|sup' ...] [--sets ...] [--scaling] [--iss-conds psf0.8 ...] [--hd-extra]
        (stages in this order by default; every stage is resumable and skips what is done; models are claimed by a pid file, so several workers can share one queue: at most two processes were
        used in total. LEARNED_SMOKE=1 with LEARNED_WEIGHTS / LEARNED_OUT in a temporary directory runs the whole pipeline with 20 steps and 20 test frames.)
    python experiments/exp_smlm_learned.py --from-raw        (regenerates results/smlm/learned.md from results/smlm/learned_raw.npz alone; numpy and scipy only)
Environment: LEARNED_WEIGHTS (weights, checkpoints and caches; OUTSIDE the repository), LEARNED_OUT (directory of the raw and the md, default results/smlm), LEARNED_THREADS (torch
threads, default 1: on an oversubscribed machine two OpenMP threads were up to 30 times slower on kernels this small), SMLM_DATA, SMLM_PICASSO_PYTHON, LEARNED_SMOKE.
"""
import argparse
import io
import json
import math
import os
import re
import sys
import time
import zlib
from collections import namedtuple

THREADS = int(os.environ.get('LEARNED_THREADS', '1'))
os.environ.setdefault('OMP_NUM_THREADS', str(THREADS))
os.environ.setdefault('MKL_NUM_THREADS', str(THREADS))
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import numpy as np

HERE = os.path.abspath(os.path.dirname(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
sys.path.insert(0, REPO)
import smlm_common as C                      # noqa: E402
import exp_smlm_iss as X                     # noqa: E402
from src.core import learned_localizer as LL     # noqa: E402

OUT = os.environ.get('LEARNED_OUT', os.path.join(REPO, 'results', 'smlm'))
RAW = os.path.join(OUT, 'learned_raw.npz')
MD = os.path.join(OUT, 'learned.md')
WEIGHTS = os.environ.get('LEARNED_WEIGHTS', os.path.join(os.path.expanduser('~'), 'smlm_learned_weights'))      # weights, checkpoints and caches: outside the repository

SETS = X.SETS
PROTOCOL = X.PROTOCOL
SEED = 20261107
REP = int(os.environ.get('LEARNED_REP', '0'))     # training replicate: 0 = the networks of the report; 1 = a second seed of the nominal networks (robustness check, separate weights directory)
ARCH = dict(c1=32, c2=64, c3=64)             # 289k parameters
BATCH = 16
LR = 2e-3
STEPS = dict(sup=2000, sim=1500)             # optimisation steps of (a) and of every simulation-trained network (b), (c)
CKPT_EVERY = 250
THRESHOLDS = (0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95)
HOLDOUT = 0.2                                # share of the development frames held out in (a)
VAL_SIM = {'LD': 600, 'HD': 150}             # simulated validation frames in (b), (c)
SWEEP_FRAMES = {'LD': 200, 'HD': 100}        # test frames of the mismatch sweep (the first of the test range): the ISS costs ~0.5 s per LD frame and ~5 s per HD frame on a free core, several times more when its
                                             # noise level is underestimated (more entries accepted), and the machine was shared (a third of a core per process)
SUP_AUG = True                               # (a): flips / rotations and small photon rescaling
TTA = {'sup': False, 'sim': False}           # test-time augmentation (average over the 8 orientations)
SCALING = (100, 200, 400)                    # data-scaling models of (a) on LD sets: fewer labelled training frames, scored on the held-out development frames only (never on test frames)
MM_SETS = ('MT0.N1.LD', 'MT0.N1.HD')
MM_TAGS = ('psf0.8', 'psf1.25', 'bg0.5', 'bg2', 'ph0.5')
ISS_CONDS = ('psf0.8', 'psf1.25', 'bg0.5', 'bg2', 'bg0.5/a', 'bg2/a')       # the sliding ISS under the same wrong inputs
PIC_FACTORS = (0.8, 1.25)                    # Picasso R2E+c with its tuned width window (min_sigma_nm, max_sigma_nm) times the factor: the only input of Picasso that can be wrong
BOOT, BLOCK = 1000, C.DEFAULT_BLOCK
SimCfg = namedtuple('SimCfg', 'psf bg ph')
if os.environ.get('LEARNED_SMOKE') == '1':              # a pipeline check: a few steps, few frames (use with LEARNED_WEIGHTS and LEARNED_OUT in a temporary directory)
    STEPS, CKPT_EVERY, VAL_SIM, SWEEP_FRAMES = dict(sup=20, sim=20), 10, {'LD': 30, 'HD': 10}, {'LD': 6, 'HD': 6}
    PROTOCOL = {s_: dict(dev=PROTOCOL[s_]['dev'], test=(PROTOCOL[s_]['test'][0], PROTOCOL[s_]['test'][0] + 20)) for s_ in SETS}

# width spread of the simulated emitters (results/smlm/calibration.md, 664 isolated emitters of MT0.N1.LD's development frames)
WIDTH_P = np.array([226, 253, 103, 82]) / 664.0
WIDTH_RATIO = np.array([1.277, 1.393, 1.656, 2.358]) / 1.369
WIDTH_LOGSD = np.array([0.04, 0.07, 0.13, 0.25])


def width_sampler(sigma0=X.SIGMA0, lo=X.WIDTHS[0], hi=X.WIDTHS[1]):
    """(rng, k) -> k core widths in px: sigma0 times a ratio of the empirical defocus classes of calibration.md, clipped to the ISS's width bounds."""
    def f(rng, k):
        b = rng.choice(4, size=k, p=WIDTH_P)
        return np.clip(sigma0 * WIDTH_RATIO[b] * np.exp(WIDTH_LOGSD[b] * rng.standard_normal(k)), lo, hi)
    return f


def tag_cfg(tag):
    """The simulation assumptions of a model tag: 'sim' nominal; 'psf0.8', 'bg2', 'ph0.5' one factor off."""
    if tag in ('sim', 'sup') or tag.startswith('sup'):
        return SimCfg(1.0, 1.0, 1.0)
    m = re.match(r'^(psf|bg|ph)([0-9.]+)$', tag)
    if not m:
        raise ValueError('unknown tag %r' % tag)
    v = float(m.group(2))
    return SimCfg(*(v if k == m.group(1) else 1.0 for k in ('psf', 'bg', 'ph')))


def model_keys(sets=SETS):
    """All models, in the order they are trained: the nominal pair of every set, then the mismatch models of the sweep sets."""
    keys = [(s, t) for s in sets for t in ('sup', 'sim')]
    keys += [(s, t) for s in sets if s in MM_SETS for t in MM_TAGS]
    return keys


def scaling_keys(sets=SETS):
    return [(s_, 'sup%d' % n) for s_ in sets if s_.endswith('LD') for n in SCALING]


def kstr(name, tag):
    return '%s|%s' % (name, tag)


def seed_of(name):
    return zlib.crc32((name if REP == 0 else '%s#%d' % (name, REP)).encode()) & 0x7fffffff


def steps_of(tag):
    return STEPS['sup' if tag.startswith('sup') else 'sim']


# ----------------------------------------------------------------------------------------------------------------------------------
# data and calibration of a set
# ----------------------------------------------------------------------------------------------------------------------------------
_FR = {}


def frames_photons(name, a, b):
    key = (name, a, b)
    if key not in _FR:
        _FR[key] = C.to_photons(C.load_frames(C.data_path(name, 'stack'), a, b))
    return _FR[key]


def gt_emitters(gt, a, b):
    """Per frame a..b-1 the labelled emitters as (k, 3) = x_px, y_px (pixel-index coordinates), photons."""
    lo, hi = np.searchsorted(gt['frame'], np.arange(a, b), 'left'), np.searchsorted(gt['frame'], np.arange(a, b), 'right')
    return [np.stack([gt['x_nm'][l:h] / C.PIXEL_NM - 0.5, gt['y_nm'][l:h] / C.PIXEL_NM - 0.5, gt['photons'][l:h]], axis=1) for l, h in zip(lo, hi)]


def gt_from_emitters(ems):
    """A ground-truth array (C.GT_DTYPE) from per-frame emitter arrays (x_px, y_px, photons)."""
    g = np.zeros(sum(len(e) for e in ems), dtype=C.GT_DTYPE)
    k = 0
    for f, e in enumerate(ems):
        for x, y, p in e:
            g[k] = (k + 1, f, (x + 0.5) * C.PIXEL_NM, (y + 0.5) * C.PIXEL_NM, 0.0, p)
            k += 1
    return g


class Ctx:
    """Everything about a set that the models need: protocol, labels of the development frames, the label-free calibration, the simulation statistics."""

    def __init__(self, name, cal=None):
        self.name = name
        self.dev, self.test = PROTOCOL[name]['dev'], PROTOCOL[name]['test']
        self.gt = C.load_gt(C.data_path(name, 'gt'))
        cal = cal or load_cal(name)
        self.bg, self.a = cal['bg'], cal['a']
        self.bg_mean = float(self.bg.mean())
        self.scale = LL.noise_scale(self.bg_mean, self.a)
        d0, d1 = self.dev
        dev_gt = self.gt[(self.gt['frame'] >= d0) & (self.gt['frame'] < d1)]
        self.pool = np.asarray(dev_gt['photons'], float)                      # photon counts of the development labels
        self.lam = len(dev_gt) / float(d1 - d0)                                # emitters per development frame
        self.kind = 'LD' if name.endswith('LD') else 'HD'
        self.n_tr = int(round((1.0 - HOLDOUT) * (d1 - d0)))                    # training frames of (a); the rest is held out
        self.pmin = C.photon_threshold(self.gt, 0.75)

    def norm(self, frames, cfg=None):
        """Network input of frames in photons under the (believed) background factor of ``cfg``."""
        f = cfg.bg if cfg else 1.0
        return LL.normalize(frames, self.bg * f, LL.noise_scale(self.bg_mean * f, self.a))


def cal_cache(name):
    return os.path.join(WEIGHTS, name.replace('|', '_'), 'cal.npz')


def load_cal(name):
    """The label-free calibration the ISS used for the set (`exp_smlm_iss.get_cal` on iss_raw.npz), cached in the weights directory so that later stages do not read the raw again."""
    p = cal_cache(name)
    if os.path.exists(p):
        with np.load(p) as z:
            return dict(bg=z['bg'], a=float(z['a']))
    for attempt in range(5):
        try:
            cal = X.get_cal(X.load_raw(), name)
            break
        except Exception:                                        # the raw may be rewritten by another job at this very moment
            time.sleep(2.0)
    else:
        raise RuntimeError('cannot read the ISS calibration')
    os.makedirs(os.path.dirname(p), exist_ok=True)
    np.savez(p, bg=cal['bg'], a=cal['a'])
    return cal


# ----------------------------------------------------------------------------------------------------------------------------------
# validation: thresholds, offsets, efficiency
# ----------------------------------------------------------------------------------------------------------------------------------
def score_rows(cands, gt, idx, offset):
    """Over THRESHOLDS: the harness metrics of the candidates (score >= threshold). With ``offset`` the constant offset (the mean signed error of the matched pairs of these very
    frames) is removed first (as the ISS operating points and the Picasso '+c' rungs). Returns a list of dicts."""
    rows = []
    for thr in THRESHOLDS:
        dets = [c[c[:, 3] >= thr] if len(c) else c for c in cands]
        off = np.zeros(2)
        if offset:
            off = np.nan_to_num(np.array([C.score(dets, gt, idx)[k] for k in ('bias_x_nm', 'bias_y_nm')]))
            dets = X.shift(dets, off)
        r = C.score(dets, gt, idx)
        rows.append(dict(thr=thr, E=r['efficiency'], J=r['jaccard'], P=r['precision'], R=r['recall'], rmse=r['rmse_nm'], dx=float(off[0]), dy=float(off[1]), n=r['n_det']))
    return rows


def best_row(rows):
    """The row of the largest E (a NaN counts as -inf; ties: the larger threshold), as in the other methods' operating points."""
    return max(rows, key=lambda r: (-np.inf if math.isnan(r['E']) else round(r['E'], 6), r['thr']))


DECODE = dict(window=0, pos_window=0, nms_px=1.5)          # the decoder: the pixel's own probability as score, the pixel's own vote as position, suppression within 1.5 px


def predict_cands(model, x, tta, batch=50):
    p, o, l = LL.predict(model, x, batch, tta)
    return LL.decode(p, o, l, threshold=THRESHOLDS[0], **DECODE)


# ----------------------------------------------------------------------------------------------------------------------------------
# training (resumable) of one model
# ----------------------------------------------------------------------------------------------------------------------------------
def key_dir(name, tag):
    return os.path.join(WEIGHTS, name.replace('|', '_'), tag.replace('/', '_'))


def sup_pool(ctx):
    d0, d1 = ctx.dev
    x = ctx.norm(frames_photons(ctx.name, d0, d1))
    return x, gt_emitters(ctx.gt, d0, d1)


def make_batch_fn(ctx, tag):
    """batch_fn(step) -> (x, targets), a function of the step only (seeded), for the model ``tag`` of the set."""
    s = seed_of(ctx.name)
    if tag.startswith('sup'):
        x, ems = sup_pool(ctx)
        ntr = ctx.n_tr if tag == 'sup' else int(tag[3:])
        if tag != 'sup':
            assert 0 < ntr <= ctx.n_tr

        def fn(step):
            rng = np.random.default_rng([SEED, s, 1, step])
            idx = rng.integers(0, ntr, BATCH)
            xs, es = [], []
            for i in idx:
                f, e = LL.augment_d4(x[i], ems[i], int(rng.integers(0, 8)) if SUP_AUG else 0)
                if SUP_AUG:
                    r = float(np.exp(rng.uniform(np.log(0.9), np.log(1.1))))
                    f, e = f * r, np.array(e, copy=True)
                    e[:, 2] *= r
                xs.append(f)
                es.append(e)
            return np.stack(xs).astype(np.float32), LL.encode_targets(es)
        return fn
    cfg = tag_cfg(tag)
    ws = width_sampler()

    def fn(step):
        rng = np.random.default_rng([SEED, s, 2, step])
        fr, ems = LL.simulate_frames(BATCH, rng, bg=ctx.bg * cfg.bg, density=ctx.lam, photons=ctx.pool, sigma=ws, psf_scale=cfg.psf, photon_scale=cfg.ph)
        return ctx.norm(fr, cfg), LL.encode_targets(ems)
    return fn


def validation_data(ctx, tag):
    """(x, gt, frame indices, offset flag, n frames, n emitters) of the validation set of a model: the held-out real development frames (sup) or fresh simulated frames."""
    d0, d1 = ctx.dev
    if tag.startswith('sup'):
        x, ems = sup_pool(ctx)
        a = ctx.n_tr
        return x[a:], ctx.gt, np.arange(d0 + a, d1), True
    cfg = tag_cfg(tag)
    rng = np.random.default_rng([SEED, seed_of(ctx.name), 3])
    n = VAL_SIM[ctx.kind]
    fr, ems = LL.simulate_frames(n, rng, bg=ctx.bg * cfg.bg, density=ctx.lam, photons=ctx.pool, sigma=width_sampler(), psf_scale=cfg.psf, photon_scale=cfg.ph)
    return ctx.norm(fr, cfg), gt_from_emitters(ems), np.arange(n), False


def n_train_stats(ctx, tag):
    """(frames, emitters) the model has been trained on / has seen: the supervised pool's size, or the simulated frames and emitters generated (STEPS x BATCH frames)."""
    if tag.startswith('sup'):
        ntr = ctx.n_tr if tag == 'sup' else int(tag[3:])
        _, ems = sup_pool(ctx)
        return ntr, int(sum(len(e) for e in ems[:ntr]))
    return steps_of(tag) * BATCH, None


def validate_ckpt(model, vdata, tta):
    x, gt, idx, offset = vdata
    return score_rows(predict_cands(model, x, tta), gt, idx, offset)


def claim(d, fname='claim.pid'):
    """Take the right to work on a model directory: True if no live process holds it. The claim is a file with the pid of the holder (created exclusively); a claim whose process is gone
    is taken over, so an interrupted run resumes and two workers started on the same queue never train the same model."""
    p = os.path.join(d, fname)
    for _ in range(2):
        try:
            fd = os.open(p, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            os.write(fd, str(os.getpid()).encode())
            os.close(fd)
            return True
        except FileExistsError:
            try:
                pid = int(open(p).read().strip() or 0)
                if pid and pid != os.getpid():
                    os.kill(pid, 0)
                    return False                                           # the holder is alive
                if pid == os.getpid():
                    return True
            except (ProcessLookupError, ValueError):
                pass
            except PermissionError:
                return False
            os.remove(p)                                                    # stale claim
    return False


_TORCH_READY = []


def torch_threads():
    """Set torch's thread counts once per process (the inter-op count can be set only before any parallel work has started)."""
    import torch
    if not _TORCH_READY:
        torch.set_num_threads(THREADS)
        try:
            torch.set_num_interop_threads(1)
        except RuntimeError:
            pass
        _TORCH_READY.append(True)
    return torch


def train_model(name, tag, log=print, max_steps=None):
    """Train (or resume) the model ``tag`` of the set; checkpoints, their validation tables and the final selection go to the weights directory. Returns the selection dict."""
    torch = torch_threads()
    ctx = Ctx(name)
    d = key_dir(name, tag)
    os.makedirs(d, exist_ok=True)
    steps = steps_of(tag) if max_steps is None else max_steps
    done_p = os.path.join(d, 'selection.json')
    if os.path.exists(done_p):
        log('%s %s: done' % (name, tag))
        return json.load(open(done_p))
    if not claim(d):
        log('%s %s: being trained by another worker, skipped' % (name, tag))
        return None
    tta = TTA['sup' if tag.startswith('sup') else 'sim']
    have = sorted(int(m.group(1)) for f in os.listdir(d) for m in [re.match(r'^opt_(\d+)\.pt$', f)]
                  if m and os.path.exists(os.path.join(d, 'ckpt_%d.npz' % int(m.group(1)))) and os.path.exists(os.path.join(d, 'log_%d.json' % int(m.group(1)))))
    model = LL.make_model(seed_of(name), **ARCH)
    opt, start, secs0 = None, 0, 0.0
    if have:
        start = have[-1]
        model, _ = LL.load(os.path.join(d, 'ckpt_%d.npz' % start))
        model = model.to(memory_format=torch.channels_last)
        opt = torch.optim.AdamW(model.parameters(), lr=LR, weight_decay=1e-4)
        opt.load_state_dict(torch.load(os.path.join(d, 'opt_%d.pt' % start), weights_only=True))
        secs0 = float(json.load(open(os.path.join(d, 'log_%d.json' % start)))['seconds'])
        log('%s %s: resuming at step %d (%.0f s so far)' % (name, tag, start, secs0))
    vdata = validation_data(ctx, tag)
    fn0 = make_batch_fn(ctx, tag)
    t0 = time.time()
    state = dict(losses=[], val_time=0.0, n_emit=0)
    if have:
        state['losses'] = json.load(open(os.path.join(d, 'loss.json')))
        state['n_emit'] = int(json.load(open(os.path.join(d, 'log_%d.json' % start)))['n_emit'])

    def fn(step):
        x, tgt = fn0(step)
        state['n_emit'] += int(tgt['n_emit'].sum())
        return x, tgt

    def cb(step, model_, opt_, hist):
        tcb = time.time()
        seg = np.array([h[1:] for h in hist[-CKPT_EVERY:]])
        secs = secs0 + (tcb - t0) - state['val_time']                                    # training seconds: validation and checkpoint writing excluded
        LL.save(model_, os.path.join(d, 'ckpt_%d.npz' % step), dict(step=step, name=name, tag=tag, seed=seed_of(name), arch=ARCH))
        torch.save(opt_.state_dict(), os.path.join(d, 'opt_%d.pt' % step))
        state['losses'].append([step] + [float(v) for v in seg.mean(axis=0)])
        json.dump(state['losses'], open(os.path.join(d, 'loss.json'), 'w'))
        json.dump(dict(step=step, seconds=secs, n_emit=state['n_emit']), open(os.path.join(d, 'log_%d.json' % step), 'w'))        # written last: marks the checkpoint complete
        for old in [f for f in os.listdir(d) if re.match(r'^opt_(\d+)\.pt$', f) and int(re.match(r'^opt_(\d+)\.pt$', f).group(1)) != step]:
            os.remove(os.path.join(d, old))
        rows = validate_ckpt(model_, vdata, tta)
        json.dump(dict(step=step, rows=rows), open(os.path.join(d, 'val_%d.json' % step), 'w'))
        b = best_row(rows)
        log('%s %s step %d/%d  %.0f s  loss %.3f  val E %.1f (thr %.2f, J %.3f, RMSE %.1f)' % (name, tag, step, steps, secs, state['losses'][-1][1], b['E'], b['thr'], b['J'], b['rmse']))
        state['val_time'] += time.time() - tcb
    if start < steps:
        LL.train(model, fn, steps, start_step=start, opt=opt, lr=LR, callback=cb, callback_every=CKPT_EVERY)
    # make sure every checkpoint has its validation table (a run interrupted between the weights and the table)
    for f in sorted(os.listdir(d)):
        m = re.match(r'^ckpt_(\d+)\.npz$', f)
        if m and not os.path.exists(os.path.join(d, 'val_%s.json' % m.group(1))):
            mdl, _ = LL.load(os.path.join(d, f))
            json.dump(dict(step=int(m.group(1)), rows=validate_ckpt(mdl, vdata, tta)), open(os.path.join(d, 'val_%s.json' % m.group(1)), 'w'))
    best = None
    for f in sorted(os.listdir(d)):
        m = re.match(r'^val_(\d+)\.json$', f)
        if m:
            v = json.load(open(os.path.join(d, f)))
            r = best_row(v['rows'])
            cand = (-np.inf if math.isnan(r['E']) else round(r['E'], 6), v['step'], r)
            if best is None or cand[:2] > best[:2]:
                best = cand
    r = best[2]
    last = json.load(open(os.path.join(d, 'log_%d.json' % max(int(re.match(r'^log_(\d+)\.json$', f).group(1)) for f in os.listdir(d) if re.match(r'^log_(\d+)\.json$', f)))))
    ntr, nem = n_train_stats(ctx, tag)
    if nem is None:
        nem = int(last['n_emit'])                                          # the emitters of all simulated training frames (counted as the batches were drawn)
    sel = dict(name=name, tag=tag, step=int(best[1]), thr=r['thr'], offset=[r['dx'], r['dy']], val=r, steps=steps, seconds=last['seconds'], n_frames=ntr, n_emitters=nem,
               seed=seed_of(name), tta=bool(tta), n_params=int(LL.n_params(model)), threads=THREADS)
    json.dump(sel, open(done_p, 'w'))
    log('%s %s: selected step %d, threshold %.2f, validation E %.1f' % (name, tag, sel['step'], sel['thr'], r['E']))
    return sel


# ----------------------------------------------------------------------------------------------------------------------------------
# inference on the test frames
# ----------------------------------------------------------------------------------------------------------------------------------
def pack(dets):
    """list of (k, 4) arrays -> (offsets (F + 1,), flat (M, 4)): the packing of exp_smlm_iss.pack."""
    return X.pack([np.asarray(d, float).reshape(-1, 4) for d in dets])


def infer_model(name, tag, log=print):
    """Candidates (score >= THRESHOLDS[0]) of the selected checkpoint on the test frames, and the inference seconds per frame (CPU, THREADS threads, whole batch)."""
    torch_threads()
    d = key_dir(name, tag)
    out_p = os.path.join(d, 'test_dets.npz')
    if os.path.exists(out_p) or re.match(r'^sup\d+$', tag):                    # (the data-scaling models are never run on test frames)
        return
    if not os.path.exists(os.path.join(d, 'selection.json')):
        log('%s %s: not trained yet, skipped' % (name, tag))
        return
    sel = json.load(open(os.path.join(d, 'selection.json')))
    ctx = Ctx(name)
    cfg = tag_cfg(tag)
    model, _ = LL.load(os.path.join(d, 'ckpt_%d.npz' % sel['step']))
    a, b = ctx.test
    x = ctx.norm(frames_photons(name, a, b), cfg)
    LL.predict(model, x[:4], 50, sel['tta'])                                    # warm-up (primitive creation), not timed
    t0 = time.perf_counter()
    cands = predict_cands(model, x, sel['tta'])
    secs = time.perf_counter() - t0
    off, flat = pack(cands)
    np.savez(out_p, off=off, flat=flat, secs=np.array([secs, b - a]), frames=np.array([a, b]))
    log('%s %s: test frames %d-%d, %d candidates, %.4f s / frame' % (name, tag, a, b - 1, len(flat), secs / (b - a)))


# ----------------------------------------------------------------------------------------------------------------------------------
# the sliding ISS under the same wrong inputs
# ----------------------------------------------------------------------------------------------------------------------------------
def iss_cache(name, cond):
    return os.path.join(WEIGHTS, name.replace('|', '_'), 'iss', cond.replace('/', '_') + '.npz')


def iss_inputs(cal, cond):
    """(calibration dict, ISS options) of a mismatch condition for the a priori default (config DEFAULT['cfg'], kappa DEFAULT['kappa'])."""
    opts = dict(rules=[('snr', X.DEFAULT['kappa'])])
    cal = dict(bg=cal['bg'], a=cal['a'])
    m = re.match(r'^(psf|bg)([0-9.]+)(/a)?$', cond)
    f = float(m.group(2))
    if m.group(1) == 'psf':
        opts['sigma0'] = X.SIGMA0 * f
    else:
        cal['bg'] = cal['bg'] * f
        if m.group(3):
            cal['a'] = cal['a'] / f
    return cal, opts


def sweep_n(name):
    return SWEEP_FRAMES['LD' if name.endswith('LD') else 'HD']


def iss_condition(name, cond, n_frames=None, log=print):
    """Run the ISS default on the first ``n_frames`` test frames with the wrong inputs of ``cond`` (one process), cached; a cache that holds fewer frames is extended by running only the
    missing ones (the frames are independent: the path of one frame does not depend on the others)."""
    p = iss_cache(name, cond)
    ctx = Ctx(name)
    a, b = ctx.test
    b = min(b, a + (n_frames or sweep_n(name)))
    old, start = None, a
    if os.path.exists(p):
        with np.load(p) as z:
            old = {k: z[k] for k in z.files}
        start = int(old['frames'][1])
        if start >= b:
            return
    os.makedirs(os.path.dirname(p), exist_ok=True)
    if not claim(os.path.dirname(p), cond.replace('/', '_') + '.claim'):
        log('ISS %s %s: being run by another worker, skipped' % (name, cond))
        return
    cal, opts = iss_inputs(dict(bg=ctx.bg, a=ctx.a), cond)
    t0 = time.time()
    dets, _, _, secs = X.run_frames(frames_photons(name, start, b), X.DEFAULT['cfg'], cal, opts, workers=1)
    off, flat = pack(dets[0])
    if old is not None:
        off, flat, secs = np.r_[old['off'], off[1:] + old['off'][-1]], np.concatenate([old['flat'], flat]), np.r_[old['secs'], secs]
    np.savez(p, off=off, flat=flat, secs=secs, frames=np.array([a, b]))
    log('ISS %s %s: frames %d-%d (%d frames in the file), %.0f s, %.2f s / frame' % (name, cond, start, b - 1, b - a, time.time() - t0, np.mean(secs)))


def iss_nominal_check(name, n=8, log=print):
    """The stored a priori default of the ISS (T0) against a fresh run of its first n test frames: the slice of the 500-frame run equals a 200-frame run (frames are independent)."""
    p = os.path.join(WEIGHTS, name.replace('|', '_'), 'iss', 'nominal_check.json')
    if os.path.exists(p):
        return
    ctx = Ctx(name)
    a, _ = ctx.test
    raw = X.load_raw()
    k0 = 'T0|%s|run|snr%g|' % (name, X.DEFAULT['kappa'])
    stored = X.unpack(raw[k0 + 'off'], raw[k0 + 'flat'])[:n]
    dets, _, _, _ = X.run_frames(frames_photons(name, a, a + n), X.DEFAULT['cfg'], dict(bg=ctx.bg, a=ctx.a), dict(rules=[('snr', X.DEFAULT['kappa'])]), workers=1)
    diff = max([float(np.abs(s - f).max()) if s.shape == f.shape and s.size else (0.0 if s.shape == f.shape else np.inf) for s, f in zip(stored, dets[0])])
    os.makedirs(os.path.dirname(p), exist_ok=True)
    json.dump(dict(frames=n, max_abs_diff=diff), open(p, 'w'))
    log('ISS nominal check %s: first %d frames, max |stored - fresh| = %g' % (name, n, diff))


def picasso_cache(name, f):
    return os.path.join(WEIGHTS, name.replace('|', '_'), 'picasso', 'window%g.npz' % f)


def picasso_window(name, f, n_frames=None, log=print):
    """Picasso with the options of its rung R2E (tuned on the development labels, `picasso_ladder_raw.npz`) but the width window times ``f``, on the first ``n_frames`` test frames,
    with the R2E offset correction; needs an interpreter that can import picasso ($SMLM_PICASSO_PYTHON). f = 1 re-runs the rung itself and is compared with the stored detections."""
    p = picasso_cache(name, f)
    if os.path.exists(p):
        return
    py = os.environ.get('SMLM_PICASSO_PYTHON')
    if not py:
        raise RuntimeError('set SMLM_PICASSO_PYTHON to an interpreter that can import picasso')
    z = np.load(X.LADDER, allow_pickle=False)
    cfg = json.loads(str(z['%s|R2E|config' % name]))
    off = z['%s|R2E|offset_nm' % name]
    for k in ('min_sigma_nm', 'max_sigma_nm'):
        if cfg[k] is not None:
            cfg[k] = cfg[k] * f
    ctx = Ctx(name)
    a, b = ctx.test
    b = min(b, a + (n_frames or sweep_n(name)))
    dets = C.run_picasso(C.load_frames(C.data_path(name, 'stack'), a, b), cfg, python=py)
    dets = C.shift_detections(dets, off[0], off[1])
    o, fl = C.pack_detections(dets)
    extra = {}
    if f == 1.0:
        stored, _ = C.load_ladder_detections(X.LADDER, name, 'R2E+c', 'test')
        extra['max_abs_diff'] = np.array(max([float(np.abs(np.asarray(s_) - np.asarray(d_)).max()) if np.shape(s_) == np.shape(d_) and np.size(s_) else (0.0 if np.shape(s_) == np.shape(d_) else np.inf)
                                              for s_, d_ in zip(stored[:b - a], [x_[:, :4] for x_ in unpack4(o, fl)])]))
    os.makedirs(os.path.dirname(p), exist_ok=True)
    np.savez(p, off=o, flat=fl, frames=np.array([a, b]), **extra)
    log('Picasso R2E+c, width window x%g, %s: %d frames, %d detections%s' % (f, name, b - a, len(fl), ', max |stored - rerun| = %g' % extra['max_abs_diff'] if extra else ''))


def unpack4(o, fl):
    return [fl[o[i]:o[i + 1]] for i in range(len(o) - 1)]


# ----------------------------------------------------------------------------------------------------------------------------------
# collect everything into the raw file
# ----------------------------------------------------------------------------------------------------------------------------------
def comparator_dets(name):
    """Test-frame detections of the rows the CNN is compared with: ISS-0 (a priori default), ISS-E (E-tuned, offset removed), Picasso R0 and R2E+c (4 columns each)."""
    raw = X.load_raw()
    gt = C.load_gt(C.data_path(name, 'gt'))
    out = {}
    k0 = 'T0|%s|run|snr%g|' % (name, X.DEFAULT['kappa'])
    out['ISS-0'] = X.unpack(raw[k0 + 'off'], raw[k0 + 'flat'])
    out['ISS-0|secs'] = np.asarray(raw['T0|%s|secs' % name], float)                  # seconds per frame of the ISS default (one core)
    out['R0|secs'] = np.asarray(raw['P|%s|secs' % name], float)                       # Picasso defaults: [seconds, frames] of the first 100 test frames
    ops = X.choose_op(X.op_table(raw, name, gt), 'E')
    out['ISS-E'] = X.op_test_dets(raw, name, ops)
    out['ISS-E|config'] = json.dumps(dict(cfg=ops['cfg'], kappa=ops['kappa'], smax=ops['smax'], off=[float(v) for v in ops['off']]))
    for row in ('R0', 'R2E+c'):
        dets, fr = C.load_ladder_detections(X.LADDER, name, row, 'test')
        out[row] = [np.asarray(d, float).reshape(-1, 4) for d in dets]
    n = PROTOCOL[name]['test'][1] - PROTOCOL[name]['test'][0]                     # (a smoke run has fewer test frames)
    for k in ('ISS-0', 'ISS-E', 'R0', 'R2E+c'):
        out[k] = out[k][:n]
    return out


def collect(sets, out_path=RAW, log=print):
    import platform
    arrays, models = {}, []
    meta = dict(seed=SEED, arch=ARCH, batch=BATCH, lr=LR, steps=STEPS, ckpt_every=CKPT_EVERY, thresholds=list(THRESHOLDS), holdout=HOLDOUT, val_sim=VAL_SIM, sweep_frames=SWEEP_FRAMES,
                sup_aug=SUP_AUG, tta=TTA, threads=THREADS, sets=list(sets), mm_sets=[s for s in MM_SETS if s in sets], mm_tags=list(MM_TAGS), iss_conds=list(ISS_CONDS), pic_factors=list(PIC_FACTORS),
                platform=platform.platform(), numpy=np.__version__)
    try:
        import torch
        meta['torch'] = torch.__version__
    except Exception:
        meta['torch'] = None
    try:
        meta['cpu'] = [l.split(':', 1)[1].strip() for l in open('/proc/cpuinfo') if l.startswith('model name')][0]
    except Exception:
        meta['cpu'] = 'unknown'
    for name in sets:
        ctx = Ctx(name)
        a, b = ctx.test
        arrays['%s|test_frames' % name] = np.array([a, b])
        arrays['%s|dev_frames' % name] = np.array(ctx.dev)
        g = ctx.gt[(ctx.gt['frame'] >= a) & (ctx.gt['frame'] < b)]
        for k in ('frame', 'x_nm', 'y_nm', 'z_nm', 'photons'):
            arrays['%s|gt|%s' % (name, k)] = np.asarray(g[k])
        arrays['%s|pmin' % name] = np.array(ctx.pmin)
        arrays['%s|calib' % name] = np.array([ctx.a, ctx.bg_mean, ctx.scale, ctx.lam, len(ctx.pool), float(ctx.pool.mean()), ctx.n_tr])
        cmp_ = comparator_dets(name)
        for k in ('ISS-0', 'ISS-E', 'R0', 'R2E+c'):
            off, flat = X.pack(cmp_[k])
            arrays['%s|cmp|%s|off' % (name, k)], arrays['%s|cmp|%s|flat' % (name, k)] = off, flat
        arrays['%s|cmp|ISS-E|config' % name] = np.array(cmp_['ISS-E|config'])
        arrays['%s|cmp|ISS-0|secs' % name], arrays['%s|cmp|R0|secs' % name] = cmp_['ISS-0|secs'], cmp_['R0|secs']
        for n_, t in model_keys([name]) + scaling_keys([name]):
            d = key_dir(n_, t)
            if not os.path.exists(os.path.join(d, 'selection.json')) or (not re.match(r'^sup\d+$', t) and not os.path.exists(os.path.join(d, 'test_dets.npz'))):
                continue
            sel = json.load(open(os.path.join(d, 'selection.json')))
            k = kstr(n_, t)
            if os.path.exists(os.path.join(d, 'test_dets.npz')):
                z = np.load(os.path.join(d, 'test_dets.npz'))
                arrays['%s|off' % k], arrays['%s|flat' % k], arrays['%s|secs' % k] = z['off'], z['flat'], z['secs']
            arrays['%s|sel' % k] = np.array(json.dumps(sel))
            arrays['%s|loss' % k] = np.array(json.load(open(os.path.join(d, 'loss.json'))))
            vals = []
            for f in sorted(os.listdir(d), key=lambda s: (len(s), s)):
                m = re.match(r'^val_(\d+)\.json$', f)
                if m:
                    v = json.load(open(os.path.join(d, f)))
                    vals.append(dict(step=v['step'], best=best_row(v['rows'])))
            arrays['%s|val' % k] = np.array(json.dumps(sorted(vals, key=lambda r: r['step'])))
            models.append(k)
        for cond in ISS_CONDS:
            p = iss_cache(name, cond)
            if os.path.exists(p):
                z = np.load(p)
                k = '%s|iss|%s' % (name, cond)
                arrays[k + '|off'], arrays[k + '|flat'], arrays[k + '|secs'], arrays[k + '|frames'] = z['off'], z['flat'], z['secs'], z['frames']
        for f in (1.0,) + tuple(PIC_FACTORS):
            p = picasso_cache(name, f)
            if os.path.exists(p):
                z = np.load(p)
                k = '%s|pic|window%g' % (name, f)
                arrays[k + '|off'], arrays[k + '|flat'] = z['off'], z['flat']
                if 'max_abs_diff' in z.files:
                    arrays[k + '|check'] = np.array(float(z['max_abs_diff']))
        pchk = os.path.join(WEIGHTS, name.replace('|', '_'), 'iss', 'nominal_check.json')
        if os.path.exists(pchk):
            arrays['%s|iss|nominal_check' % name] = np.array(open(pchk).read())
    rep_dir = os.environ.get('LEARNED_REP_WEIGHTS')
    rep_models = []
    if rep_dir:
        for name in sets:
            for t in ('sup', 'sim'):
                d = os.path.join(rep_dir, name.replace('|', '_'), t)
                if os.path.exists(os.path.join(d, 'selection.json')) and os.path.exists(os.path.join(d, 'test_dets.npz')):
                    z = np.load(os.path.join(d, 'test_dets.npz'))
                    k = kstr(name, t) + '#rep1'
                    arrays['%s|off' % k], arrays['%s|flat' % k] = z['off'], z['flat']
                    arrays['%s|sel' % k] = np.array(open(os.path.join(d, 'selection.json')).read())
                    rep_models.append(k)
    meta['rep_models'] = rep_models
    meta['models'] = [m_ for m_ in models if not re.match(r'^.*\|sup\d+$', m_)]
    meta['scaling_models'] = [m_ for m_ in models if re.match(r'^.*\|sup\d+$', m_)]
    arrays['meta'] = np.array(json.dumps(meta, sort_keys=True))
    C.save_npz_deterministic(out_path, arrays)
    log('wrote %s (%.1f MB, %d models)' % (out_path, os.path.getsize(out_path) / 1e6, len(models)))


# ----------------------------------------------------------------------------------------------------------------------------------
# report (from the raw alone)
# ----------------------------------------------------------------------------------------------------------------------------------
def load_learned_raw(path=RAW):
    with np.load(path, allow_pickle=False) as z:
        return {k: z[k] for k in z.files}


def gt_of(raw, name):
    g = np.zeros(len(raw['%s|gt|frame' % name]), dtype=C.GT_DTYPE)
    g['id'] = np.arange(1, len(g) + 1)
    g['frame'], g['x_nm'], g['y_nm'], g['photons'] = (raw['%s|gt|%s' % (name, k)] for k in ('frame', 'x_nm', 'y_nm', 'photons'))
    if '%s|gt|z_nm' % name in raw:
        g['z_nm'] = raw['%s|gt|z_nm' % name]
    return g


def dets_of(raw, key):
    return X.unpack(raw[key + '|off'], raw[key + '|flat'])


def apply_sel(raw, k, n_frames=None):
    """The test detections of model k at its selected threshold with its selected offset removed."""
    sel = json.loads(str(raw[k + '|sel']))
    d = [c[c[:, 3] >= sel['thr']] if len(c) else c for c in dets_of(raw, k)]
    d = X.shift(d, np.array(sel['offset'])) if sel['tag'].startswith('sup') else d
    return d[:n_frames] if n_frames else d


def pct(v, ci, nd=1):
    return X.fmt_ci(v, ci, 100.0, nd)


def f1(v, ci):
    return X.fmt_ci(v, ci, 1.0, 1)


def score_both(dets, gt, idx, pmin):
    return C.score(dets, gt, idx), C.score(dets, gt, idx, border_nm=450.0, photon_min=pmin)


def label_of(k):
    name, tag = k.split('|')
    return {'sup': 'CNN, supervised on the development labels (a)', 'sim': 'CNN, trained on simulation only (b)'}.get(tag, 'CNN, simulation with %s (c)' % tag)


def report(raw, boot=BOOT):
    meta = json.loads(str(raw['meta']))
    sets = meta['sets']
    L = ['# A CPU-scale learned localizer against the sliding inverse scale space and Picasso, in distribution and under a mismatch', '',
         'Generated by `experiments/exp_smlm_learned.py --from-raw` from `learned_raw.npz`; every number below is printed by the script from that file (the ground truth of the test frames and '
         'the detections of the comparison rows are stored in it). The design, declared before any test-frame evaluation, is in the docstring of the script and summarised in '
         '"Design and fairness" at the end. Frame indices are 0-based stack indices, positions in nm (corner convention), metrics pooled over the frames with 95 %% intervals of a bootstrap in '
         'blocks of %d consecutive frames (%d resamples, seed 0); differences are paired (the same resampled blocks). Network: space-to-depth U-Net, %s parameters, trained with torch %s on %s '
         '(%d thread(s)).' % (BLOCK, boot, '{:,}'.format(json.loads(str(raw[meta['models'][0] + '|sel']))['n_params']) if meta['models'] else '?', meta.get('torch'), meta.get('cpu'), meta['threads']), '']
    res, dets_by = {}, {}
    for name in sets:
        gt = gt_of(raw, name)
        a, b = (int(v) for v in raw['%s|test_frames' % name])
        idx = np.arange(a, b)
        pmin = float(raw['%s|pmin' % name])
        meths = [('ISS-0', 'sliding ISS, a priori default (hp5h, kappa 5, no width bound, no offset correction)', dets_of(raw, '%s|cmp|ISS-0' % name)),
                 ('ISS-E', 'sliding ISS, E-tuned on the development labels (%s), offset removed' % json.loads(str(raw['%s|cmp|ISS-E|config' % name]))['cfg'], dets_of(raw, '%s|cmp|ISS-E' % name)),
                 ('R0', 'Picasso, defaults (R0)', dets_of(raw, '%s|cmp|R0' % name)),
                 ('R2E+c', 'Picasso, every option tuned for E on the development labels + offset correction (R2E+c)', dets_of(raw, '%s|cmp|R2E+c' % name))]
        meths = [(tag, label_of(kstr(name, tag)), apply_sel(raw, kstr(name, tag))) for tag in ('sup', 'sim') if kstr(name, tag) + '|sel' in raw] + meths
        res[name], dets_by[name] = {}, {}
        for key, m, dets in meths:
            assert len(dets) == len(idx), (name, key, len(dets), len(idx))
            dets_by[name][key] = dets
            r, r2 = score_both(dets, gt, idx, pmin)
            res[name][key] = (m, r, r2, C.bootstrap_ci(r, B=boot, seed=0, block=BLOCK), C.bootstrap_ci(r2, B=boot, seed=0, block=BLOCK))
    L += summary_lines(raw, res, sets, boot)
    L += replicate_lines(raw, res, sets, boot)
    for name in sets:
        L += set_lines(raw, name, res[name], boot, dets_by[name])
    L += sweep_lines(raw, res, boot)
    L += design_lines(raw)
    return '\n'.join(L) + '\n'


CNN_KEYS = ('sup', 'sim')


def summary_lines(raw, res, sets, boot):
    L = ['## Summary: efficiency E and Jaccard index J on the test frames (point estimates; intervals in the per-set tables)', '',
         'E = 100 - sqrt((100 - J)^2 + RMSE^2) (J in %, RMSE in nm), all ground-truth emitters counted. Rows marked (a) use the development LABELS to choose the checkpoint, the threshold and a constant offset; '
         'row (b) uses no label on the frames (a simulation with the density and photon counts of the development labels); the ISS a priori default and Picasso R0 are label-free, the ISS E-tuned and '
         'Picasso R2E+c are tuned on the development labels.', '']
    order = [('sup', 'CNN (a) supervised'), ('sim', 'CNN (b) simulation-trained'), ('ISS-0', 'ISS a priori default'), ('ISS-E', 'ISS E-tuned + offset'), ('R0', 'Picasso defaults R0'), ('R2E+c', 'Picasso R2E+c')]
    hdr = ['method'] + [s + ' E (J %)' for s in sets]
    rows = []
    for key, lab in order:
        row = [lab]
        for s in sets:
            if key in res[s]:
                r = res[s][key][1]
                row.append('%.1f (%.1f)' % (r['efficiency'], 100 * r['jaccard']))
            else:
                row.append('-')
        rows.append(row)
    L += [X.table(hdr, rows), '']
    full = [s_ for s_ in sets if all(k in res[s_] for k in ('sup', 'sim'))]
    if full:
        mean = lambda key: np.mean([res[s_][key][1]['efficiency'] for s_ in full])
        L += ['Mean E over the %d sets: %s.' % (len(full), '; '.join('%s %.1f' % (lab, mean(key)) for key, lab in order)), '']
    nm = {'sup': 'CNN (a)', 'sim': 'CNN (b)', 'ISS-0': 'ISS default', 'ISS-E': 'ISS E-tuned+c', 'R0': 'Picasso R0', 'R2E+c': 'Picasso R2E+c'}
    L += ['Paired difference of E (CNN minus comparator, all emitters) with its 95 % interval and the verdict (win: the interval excludes 0 in the CNN\'s favour; loss: against it; tie: it contains 0). '
          'Like-for-like pairs first (label use matches), then the unequal ones.', '']
    drows = []
    for ka, kb, note in (('sup', 'ISS-E', 'like-for-like'), ('sup', 'R2E+c', 'like-for-like'), ('sim', 'ISS-0', 'like-for-like'), ('sim', 'R0', 'like-for-like'), ('sup', 'ISS-0', 'CNN has the labels'),
                         ('sim', 'ISS-E', 'comparator has the labels')):
        row = ['%s - %s' % (nm[ka], nm[kb]), note]
        for s_ in sets:
            if ka in res[s_] and kb in res[s_]:
                ra, rb = res[s_][ka][1], res[s_][kb][1]
                d = C.bootstrap_diff_ci(ra['per_frame'], rb['per_frame'], B=boot, seed=0, block=BLOCK)['efficiency']
                row.append('%s %s' % (f1(ra['efficiency'] - rb['efficiency'], d), X.wtl(ra['efficiency'] - rb['efficiency'], d)))
            else:
                row.append('-')
        drows.append(row)
    L += [X.table(['comparison', 'fairness'] + list(sets), drows), '']
    return L


def replicate_lines(raw, res, sets, boot):
    """Robustness: the nominal networks trained again with another seed (initial weights, simulated or augmented batches, validation frames) and the same recipe; the test frames evaluated once."""
    meta = json.loads(str(raw['meta']))
    reps = meta.get('rep_models', [])
    if not reps:
        return []
    L = ['## Robustness: a second training seed of the nominal networks', '',
         'The networks of (a) and (b) trained again from another seed (the initial weights, the augmentation or simulation draws and the simulated validation frames change; the recipe, steps, and the choice of '
         'checkpoint and threshold on the validation frames do not), scored on the same test frames. This is the seed-to-seed variability of one cell of the tables above; the paired interval is a bootstrap over frames '
         'for the SAME two networks, so it does not contain the seed variability (the difference of the two point estimates is the measure of that).', '']
    rows = []
    for name in sets:
        gt = gt_of(raw, name)
        a, b = (int(v) for v in raw['%s|test_frames' % name])
        idx = np.arange(a, b)
        pmin = float(raw['%s|pmin' % name])
        for tag in ('sup', 'sim'):
            k = kstr(name, tag) + '#rep1'
            if k not in reps or tag not in res[name]:
                continue
            sel = json.loads(str(raw[k + '|sel']))
            r1 = score_both(apply_sel(raw, k), gt, idx, pmin)[0]
            r0 = res[name][tag][1]
            d = C.bootstrap_diff_ci(r1['per_frame'], r0['per_frame'], B=boot, seed=0, block=BLOCK)
            rows.append([name, label_of(kstr(name, tag)), '%.1f (%.1f)' % (r0['efficiency'], 100 * r0['jaccard']), '%.1f (%.1f)' % (r1['efficiency'], 100 * r1['jaccard']), f1(r1['efficiency'] - r0['efficiency'], d['efficiency']),
                         '%d / %.2f' % (sel['step'], sel['thr'])])
    L += [X.table(['set', 'network', 'seed 0: E (J %)', 'seed 1: E (J %)', 'seed 1 - seed 0: dE (paired over frames)', 'seed 1: selected step / threshold'], rows), '']
    if rows:
        diffs = [float(r[4].split()[0]) for r in rows]
        L += ['Largest |dE| between the two seeds: %.1f; mean |dE|: %.1f (%d networks).' % (max(abs(v) for v in diffs), np.mean([abs(v) for v in diffs]), len(diffs)), '']
    return L


def set_lines(raw, name, R, boot, dets_by):
    gt = gt_of(raw, name)
    a, b = (int(v) for v in raw['%s|test_frames' % name])
    idx = np.arange(a, b)
    cal = raw['%s|calib' % name]
    L = ['## %s' % name, '',
         'Test frames %d-%d (%d frames, %.3f ground-truth emitters per frame). Calibration used (label-free, the ISS\'s): variance slope a = %.3f, mean background %.2f photons, noise scale s = %.2f photons. Simulation statistics '
         '(development labels): %.3f emitters per frame, %d photon counts (mean %.0f). Held-out development frames of (a): the last %d%% of the development frames.\n'
         % (a, b - 1, b - a, len(gt) / float(b - a), cal[0], cal[1], cal[2], cal[3], int(cal[4]), cal[5], int(round(100 * HOLDOUT)))]
    L.append('### Training and inference cost\n')
    rows = []
    for tag in ('sup', 'sim'):
        k = kstr(name, tag)
        if k + '|sel' not in raw:
            continue
        sel = json.loads(str(raw[k + '|sel']))
        secs = raw[k + '|secs']
        rows.append([label_of(k), sel['n_frames'], sel['n_emitters'], '%d x %d' % (sel['steps'], BATCH), '%.0f' % sel['seconds'], sel['step'], '%.2f' % sel['thr'],
                     '(%.1f, %.1f)' % tuple(sel['offset']) if sel['tag'].startswith('sup') else 'none', '%.1f' % sel['val']['E'], '%.4f' % (secs[0] / secs[1])])
    if '%s|cmp|ISS-0|secs' % name in raw:
        rows.append(['sliding ISS, a priori default (no training; one core)', '-', '-', '-', '-', '-', '-', '-', '-', '%.4f' % float(np.mean(raw['%s|cmp|ISS-0|secs' % name]))])
        rows.append(['Picasso, defaults (no training; one core, identification and fit)', '-', '-', '-', '-', '-', '-', '-', '-', '%.4f' % (raw['%s|cmp|R0|secs' % name][0] / raw['%s|cmp|R0|secs' % name][1])])
    L.append(X.table(['model', 'training frames', 'training emitters', 'steps x batch', 'training s (%d thread)' % json.loads(str(raw['meta']))['threads'], 'selected step', 'threshold',
                      'offset removed (dx, dy) nm', 'validation E', 'inference s / frame'], rows))
    L.append('\nTraining frames / emitters: (a) the labelled development frames used for training (the held-out 20 %% are not in them); (b) the simulated frames and emitters generated (nothing is stored; every step '
             'draws fresh frames). Validation E: on the held-out development frames with the offset removed (a) or on simulated validation frames (b), at the selected step and threshold. Inference '
             'seconds per frame: network and decoder, CPU, all test frames in batches of 50, %d thread(s), on a machine shared with other jobs (a wall-clock measurement, noisy).\n' % json.loads(str(raw['meta']))['threads'])
    L += scaling_lines(raw, name)
    L.append('### Test frames: all ground-truth emitters (primary setting)\n')
    L.append(X.table(['method', 'J %', 'precision %', 'recall %', 'RMSE nm', 'efficiency', 'bias (dx, dy) nm', 'detected / frame'],
                     [[m, pct(r['jaccard'], ci['jaccard']), pct(r['precision'], ci['precision']), pct(r['recall'], ci['recall']), f1(r['rmse_nm'], ci['rmse_nm']), f1(r['efficiency'], ci['efficiency']),
                       '(%.1f, %.1f)' % (r['bias_x_nm'], r['bias_y_nm']), '%.3f' % (r['n_det'] / len(idx))] for m, r, r2, ci, ci2 in R.values()]))
    L.append('\nChallenge assessment setting (emitters within 450 nm of the border and the dimmest quarter of the activations dropped after the pairing):\n')
    L.append(X.table(['method', 'J %', 'RMSE nm', 'efficiency'], [[m, pct(r2['jaccard'], ci2['jaccard']), f1(r2['rmse_nm'], ci2['rmse_nm']), f1(r2['efficiency'], ci2['efficiency'])] for m, r, r2, ci, ci2 in R.values()]))
    L.append('\n### Paired differences (CNN minus comparator, same frames, same resampled blocks), 95 % intervals\n')
    rows = []
    pairs = [('sup', 'ISS-E', 'both tuned on the development labels'), ('sup', 'R2E+c', 'both tuned on the development labels'), ('sup', 'ISS-0', 'unequal: the CNN uses labels'), ('sup', 'R0', 'unequal: the CNN uses labels'),
             ('sim', 'ISS-0', 'both label-free on the frames'), ('sim', 'R0', 'both label-free on the frames'), ('sim', 'ISS-E', 'unequal: the comparator uses labels'), ('sim', 'R2E+c', 'unequal: the comparator uses labels'),
             ('sim', 'sup', 'simulation-trained minus supervised')]
    for ka, kb, note in pairs:
        if ka not in R or kb not in R:
            continue
        ra, rb = R[ka], R[kb]
        d = C.bootstrap_diff_ci(ra[1]['per_frame'], rb[1]['per_frame'], B=boot, seed=0, block=BLOCK)
        d2 = C.bootstrap_diff_ci(ra[2]['per_frame'], rb[2]['per_frame'], B=boot, seed=0, block=BLOCK)
        nm = {'sup': 'CNN (a)', 'sim': 'CNN (b)', 'ISS-0': 'ISS default', 'ISS-E': 'ISS E-tuned+c', 'R0': 'Picasso R0', 'R2E+c': 'Picasso R2E+c'}
        rows.append(['%s - %s' % (nm[ka], nm[kb]), note, pct(ra[1]['jaccard'] - rb[1]['jaccard'], d['jaccard']), f1(ra[1]['rmse_nm'] - rb[1]['rmse_nm'], d['rmse_nm']),
                     f1(ra[1]['efficiency'] - rb[1]['efficiency'], d['efficiency']), f1(ra[2]['efficiency'] - rb[2]['efficiency'], d2['efficiency'])])
    L.append(X.table(['comparison', 'fairness', 'dJ % (all)', 'dRMSE nm (all)', 'dE (all)', 'dE (challenge)'], rows))
    L.append('')
    L += recall_lines(R, dets_by, gt, idx)
    L += z_lines(R, dets_by, gt, idx)
    L += threshold_lines(raw, name, gt, idx)
    L += curve_lines(raw, name)
    return L


def threshold_lines(raw, name, gt, idx):
    """Descriptive only: the test-frame efficiency of the networks at every threshold of the grid (the selected offset removed for (a)). The threshold of each network was chosen on its own
    validation frames; this table shows how far that choice is from the best threshold in hindsight, which is NOT an operating point anyone could have chosen."""
    thr = json.loads(str(raw['meta']))['thresholds']
    rows = []
    for tag in ('sup', 'sim'):
        k = kstr(name, tag)
        if k + '|sel' not in raw:
            continue
        sel = json.loads(str(raw[k + '|sel']))
        cands = dets_of(raw, k)
        es = []
        for t in thr:
            d = [c[c[:, 3] >= t] if len(c) else c for c in cands]
            if tag == 'sup':
                d = X.shift(d, np.array(sel['offset']))
            es.append(C.score(d, gt, idx)['efficiency'])
        best = max(range(len(thr)), key=lambda i: -np.inf if math.isnan(es[i]) else es[i])
        rows.append([label_of(k)] + [('**%.1f**' if t == sel['thr'] else '%.1f') % e for t, e in zip(thr, es)] + ['%.2f (E %.1f)' % (thr[best], es[best])])
    if not rows:
        return []
    L = ['### Threshold sensitivity on the test frames (descriptive only; the thresholds were chosen on validation frames, bold)\n',
         X.table(['network'] + ['thr %.2f' % t for t in thr] + ['best in hindsight'], rows),
         '\nThe last column uses the test labels and is not an operating point; it shows how much of the gap to a tuned method is a mis-set threshold and how much is the quality of the detections.\n']
    return L


def z_lines(R, dets_by, gt, idx):
    """Recall and RMSE of the matched pairs by the axial distance |z| of the emitter from the focal plane (the real PSF broadens and shifts with |z|; the simulation of (b) has no z): the
    defocus classes of calibration.md. Descriptive."""
    if not np.any(gt['z_nm'] != 0):
        return []
    z = np.abs(np.asarray(gt['z_nm']))
    edges = [0.0, 100.0, 200.0, 350.0, np.inf]
    cols = ['|z| %s-%s nm (n = %d)' % ('%.0f' % lo, '%.0f' % hi if np.isfinite(hi) else 'max', int(((z >= lo) & (z < hi)).sum())) for lo, hi in zip(edges[:-1], edges[1:])]
    rows = []
    for key, (m, r, r2, ci, ci2) in R.items():
        mm = C.matched_errors(dets_by[key], gt, idx)
        zz = z[mm['gt_row']]
        cells = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            sel = (zz >= lo) & (zz < hi)
            ok = sel & mm['matched']
            cells.append('%.1f / %.1f' % (100.0 * mm['matched'][sel].mean(), math.sqrt(float(np.mean(mm['dx'][ok] ** 2 + mm['dy'][ok] ** 2)))) if ok.any() else '- / -')
        rows.append([m] + cells)
    return ['### Recall % / RMSE nm of the matched pairs by the axial position of the emitter (test frames; descriptive)\n', X.table(['method'] + cols, rows), '']


def recall_lines(R, dets_by, gt, idx):
    """Where the methods lose emitters: recall by the photon count of the ground-truth emitter, in bins of the test frames' own photon quantiles (descriptive, no choice depends on it)."""
    ph = np.asarray(gt['photons'])
    edges = np.quantile(ph, [0.0, 0.125, 0.25, 0.5, 0.75, 1.0])
    edges[-1] = np.inf
    L = ['### Recall by the photon count of the emitter (test frames; bins are quantiles of the test frames\' ground-truth photons)\n']
    cols, rows = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        cols.append('%.0f-%s photons (n = %d)' % (lo, '%.0f' % hi if np.isfinite(hi) else 'max', int(((ph >= lo) & (ph < hi)).sum())))
    for key, (m, r, r2, ci, ci2) in R.items():
        mm = C.matched_errors(dets_by[key], gt, idx)
        p_ = ph[mm['gt_row']]
        rows.append([m] + ['%.1f' % (100.0 * mm['matched'][(p_ >= lo) & (p_ < hi)].mean()) for lo, hi in zip(edges[:-1], edges[1:])])
    L += [X.table(['method (recall %)'] + cols, rows), '']
    return L


def scaling_lines(raw, name):
    """How much labelled data the supervised network needs: the held-out development E for 100, 200, 400 and all (800) training frames (never evaluated on test frames)."""
    ks = [kstr(name, 'sup%d' % n) for n in SCALING] + [kstr(name, 'sup')]
    ks = [k for k in ks if k + '|sel' in raw]
    if len(ks) < 2:
        return []
    L = ['### How much labelled data (a) needs: the held-out development frames only (200 frames; the test frames are not used here)\n']
    rows = []
    for k in ks:
        sel = json.loads(str(raw[k + '|sel']))
        v = sel['val']
        rows.append([sel['n_frames'], sel['n_emitters'], sel['step'], '%.2f' % sel['thr'], '%.1f' % (100 * v['J']), '%.1f' % v['rmse'], '%.1f' % v['E']])
    L.append(X.table(['labelled training frames', 'labelled emitters', 'selected step', 'threshold', 'held-out J %', 'held-out RMSE nm', 'held-out E (offset removed)'], rows))
    L.append('\nThe held-out E is the quantity the checkpoint and threshold were chosen on (so it is optimistic for every row alike); the number of steps is the same for all rows. The ISS and Picasso need no such data.\n')
    return L


def curve_lines(raw, name):
    L = ['### Training curves (validation efficiency E of the best threshold at each checkpoint)\n']
    cur = []
    for tag in ('sup', 'sim'):
        k = kstr(name, tag)
        if k + '|val' in raw:
            cur.append((label_of(k), json.loads(str(raw[k + '|val']))))
    if cur:
        steps = sorted({v['step'] for _, vals in cur for v in vals})
        rows = [[lab] + ['%.1f' % {v['step']: v['best']['E'] for v in vals}[st] if st in {v['step'] for v in vals} else '-' for st in steps] for lab, vals in cur]
        L += [X.table(['model'] + ['step %d' % st for st in steps], rows), '']
    return L


def sweep_lines(raw, res, boot):
    meta = json.loads(str(raw['meta']))
    L = []
    for name in meta['mm_sets']:
        if kstr(name, 'sim') + '|sel' not in raw:
            continue
        gt = gt_of(raw, name)
        a, b = (int(v) for v in raw['%s|test_frames' % name])
        nf = min(meta['sweep_frames']['LD' if name.endswith('LD') else 'HD'], b - a)
        idx = np.arange(a, a + nf)
        pmin = float(raw['%s|pmin' % name])
        def sub(k):
            return apply_sel(raw, k, nf)
        entries = {}                                              # label -> dets over the sweep frames
        entries['nominal'] = {'CNN': sub(kstr(name, 'sim')), 'ISS': dets_of(raw, '%s|cmp|ISS-0' % name)[:nf]}
        for tag in meta['mm_tags']:
            k = kstr(name, tag)
            if k + '|sel' in raw:
                entries.setdefault(tag, {})['CNN'] = sub(k)
        for cond in meta['iss_conds']:
            k = '%s|iss|%s' % (name, cond)
            if k + '|off' in raw:
                got = dets_of(raw, k)
                entries.setdefault(cond, {})['ISS'] = got[:nf]
        entries.setdefault('ph0.5', {})['ISS'] = entries['nominal']['ISS']                       # no photon input: the nominal run
        rp = {key: dets_of(raw, '%s|cmp|%s' % (name, key))[:nf] for key in ('R0', 'R2E+c')}
        sc = {}
        for lab, ent in entries.items():
            for meth, dets in ent.items():
                r, r2 = score_both(dets, gt, idx, pmin)
                sc[(lab, meth)] = (r, r2)
        for key, dets in rp.items():
            sc[('-', key)] = score_both(dets, gt, idx, pmin)
        L += ['## Mismatch sweep (C9): %s, test frames %d-%d (%d frames)' % (name, a, a + nf - 1, nf), '',
              'The network is retrained with the SAME recipe on a simulation whose assumption is wrong by the factor shown (PSF width, background level, photon counts; the checkpoint and threshold are chosen on the wrong '
              'simulation\'s own validation frames); the sliding ISS (a priori default, kappa 5) gets the same wrong inputs (its entry width times the factor; its background surface times the factor with the slope a unchanged, and, marked, with the slope divided by the factor; it has no photon prior, so its row for the photon counts is the nominal run); Picasso has no such inputs (its default rung R0 is constant; the width window of R2E+c is re-run below). All rows are scored on the SAME %d test frames '
              '(the first %d of the test range; the ISS is slow, ~0.5 s per frame for LD and several seconds for HD); every interval is a bootstrap in blocks of %d frames, differences are paired.\n' % (nf, nf, BLOCK)]
        order = [('psf0.8', 'PSF width', '0.8'), ('nominal', 'PSF width', '1 (nominal)'), ('psf1.25', 'PSF width', '1.25'),
                 ('bg0.5', 'background', '0.5'), ('bg0.5/a', 'background (ISS slope a / 0.5)', '0.5'), ('nominal', 'background', '1 (nominal)'), ('bg2', 'background', '2'), ('bg2/a', 'background (ISS slope a / 2)', '2'),
                 ('ph0.5', 'photon counts', '0.5'), ('nominal', 'photon counts', '1 (nominal)')]
        rows, drows = [], []
        for lab, kind, fac in order:
            cells = [kind, fac]
            for meth in ('CNN', 'ISS'):
                if (lab, meth) in sc:
                    r, r2 = sc[(lab, meth)]
                    ci = C.bootstrap_ci(r, B=boot, seed=0, block=BLOCK)
                    cells += [f1(r['efficiency'], ci['efficiency']), pct(r['jaccard'], ci['jaccard']), f1(r['rmse_nm'], ci['rmse_nm'])]
                else:
                    cells += ['-', '-', '-']
            rows.append(cells)
            if lab != 'nominal':
                drow = [kind, fac]
                for meth in ('CNN', 'ISS'):
                    if (lab, meth) in sc and ('nominal', meth) in sc:
                        ra, rb = sc[(lab, meth)][0], sc[('nominal', meth)][0]
                        d = C.bootstrap_diff_ci(ra['per_frame'], rb['per_frame'], B=boot, seed=0, block=BLOCK)
                        drow += [f1(ra['efficiency'] - rb['efficiency'], d['efficiency']), pct(ra['jaccard'] - rb['jaccard'], d['jaccard'])]
                    else:
                        drow += ['-', '-']
                if (lab, 'CNN') in sc and (lab, 'ISS') in sc:
                    ra, rb = sc[(lab, 'CNN')][0], sc[(lab, 'ISS')][0]
                    d = C.bootstrap_diff_ci(ra['per_frame'], rb['per_frame'], B=boot, seed=0, block=BLOCK)
                    drow += [f1(ra['efficiency'] - rb['efficiency'], d['efficiency']), pct(ra['jaccard'] - rb['jaccard'], d['jaccard'])]
                else:
                    drow += ['-', '-']
                drows.append(drow)
        L.append('Efficiency, Jaccard index and RMSE by condition:\n')
        L.append(X.table(['assumption (or input) that is wrong', 'factor', 'CNN E', 'CNN J %', 'CNN RMSE nm', 'ISS E', 'ISS J %', 'ISS RMSE nm'], rows))
        L.append('\nPicasso (no mismatch possible), same frames: ' + '; '.join(
            '%s E %s, J %s, RMSE %s' % (key, f1(sc[('-', key)][0]['efficiency'], C.bootstrap_ci(sc[('-', key)][0], B=boot, seed=0, block=BLOCK)['efficiency']),
                                        pct(sc[('-', key)][0]['jaccard'], C.bootstrap_ci(sc[('-', key)][0], B=boot, seed=0, block=BLOCK)['jaccard']),
                                        f1(sc[('-', key)][0]['rmse_nm'], C.bootstrap_ci(sc[('-', key)][0], B=boot, seed=0, block=BLOCK)['rmse_nm'])) for key in ('R0', 'R2E+c')) + '.\n')
        L.append('Degradation curve: change from the nominal run on the same frames (paired), and the CNN minus the ISS under the same mismatch (paired):\n')
        L.append(X.table(['assumption (or input) that is wrong', 'factor', 'CNN dE vs nominal', 'CNN dJ %', 'ISS dE vs nominal', 'ISS dJ %', 'CNN - ISS dE (same condition)', 'CNN - ISS dJ %'], drows))
        if '%s|pic|window1' % name + '|off' in raw:
            prow = []
            base = score_both(dets_of(raw, '%s|pic|window1' % name)[:nf], gt, idx, pmin)[0]
            for f in [1.0] + list(meta['pic_factors']):
                k = '%s|pic|window%g' % (name, f)
                if k + '|off' not in raw:
                    continue
                r = score_both(dets_of(raw, k)[:nf], gt, idx, pmin)[0]
                ci = C.bootstrap_ci(r, B=boot, seed=0, block=BLOCK)
                d = C.bootstrap_diff_ci(r['per_frame'], base['per_frame'], B=boot, seed=0, block=BLOCK) if f != 1.0 else None
                prow.append(['x%g%s' % (f, ' (the rung itself)' if f == 1.0 else ''), f1(r['efficiency'], ci['efficiency']), pct(r['jaccard'], ci['jaccard']), f1(r['rmse_nm'], ci['rmse_nm']),
                             f1(r['efficiency'] - base['efficiency'], d['efficiency']) if d else '-', pct(r['jaccard'] - base['jaccard'], d['jaccard']) if d else '-'])
            L.append('\nPicasso under a PSF-width mismatch: the only input it has is the width window of its tuned post-filter (R2E+c: min / max fitted sigma, tuned on the development labels); the window times a factor '
                     '(re-run of Picasso on these frames; the factor 1 re-runs the rung itself and reproduces the stored detections, largest absolute difference %g):\n' % float(raw['%s|pic|window1|check' % name]))
            L.append(X.table(['window', 'E', 'J %', 'RMSE nm', 'dE vs the rung (paired)', 'dJ %'], prow))
        chk = '%s|iss|nominal_check' % name
        if chk in raw:
            c = json.loads(str(raw[chk]))
            L.append('\nCheck: a fresh ISS run of the first %d test frames equals the stored a priori default (T0) there, largest absolute difference %g (frames are independent, so the first %d frames of the stored run are the '
                     'nominal ISS of this sweep).\n' % (c['frames'], c['max_abs_diff'], nf))
    return L


def design_lines(raw):
    meta = json.loads(str(raw['meta']))
    st = meta['steps']
    return ['## Design and fairness', '',
            '* Protocol, simulation, mismatch definitions and the list of what was chosen where are in the docstring of `experiments/exp_smlm_learned.py` (declared before any test-frame evaluation).',
            '* Network: space-to-depth U-Net with widths (%d, %d, %d); AdamW (weight decay 1e-4, gradient-norm clip 2; loss = BCE + 2 x Huber offsets + 0.1 x Huber log photons), learning rate %g (warm-up 100 steps, cosine decay to 2 %%), batch %d, %d steps for (a) and %d for (b), (c), a checkpoint every %d steps; seeds fixed (`SEED` %d, the initial '
            'weights and the simulated batches per set from the set\'s name, so the networks of one set that differ only in the simulation\'s assumption share their random numbers).' % (
                meta['arch']['c1'], meta['arch']['c2'], meta['arch']['c3'], meta['lr'], meta['batch'], st['sup'], st['sim'], meta['ckpt_every'], meta['seed']),
            '* Fair comparisons: the supervised network (a) uses the development labels exactly as the E-tuned ISS and Picasso R2E+c do (checkpoint / threshold / offset on held-out or development frames, the test frames '
            'evaluated once), so (a) against ISS-E and R2E+c is the like-for-like comparison. The simulation-trained network (b) uses no label on the frames, so it is compared with the label-free rows (ISS a priori '
            'default, Picasso R0); against the tuned rows it is the comparator that has the advantage. The tuning effort differs in kind: the network has one architecture, one schedule, %d thresholds and %d / %d '
            'checkpoints per model, the ISS 168 operating points, Picasso 12,000-15,000 configurations; none of them is tuned on test frames.' % (
                len(meta['thresholds']), st['sup'] // meta['ckpt_every'], st['sim'] // meta['ckpt_every']),
            '* Picasso under mismatch: it has no PSF, background or photon input, so its default rung R0 is constant in the sweeps; the only input of its tuned rung R2E+c that can be wrong is the width window of its '
            'post-filter, which was re-run with the window times 0.8 and 1.25 (a table in each sweep).',
            '* The sweeps use the first 200 (LD) and first 100 (HD) test frames because the ISS is slow (and several times slower when it underestimates the noise). The ISS conditions that repeat a background mismatch with the '
            'slope recomputed (`a / factor`) were run on the LD set only.',
            '* Not comparable with the CNN\'s numbers: DECODE\'s published leaderboard results (3D only, hidden test sets).',
            '* The MT0 sets are simulated and share one structure (three microtubules): the development and test frames are not independent datasets. A network trained with flips and rotations on frames cannot '
            'memorise the absolute position of the microtubules, but it may learn their local shape; the simulation-trained network has seen no microtubule at all.',
            '* Limits of the evidence: one network per cell (a single seed; the jumps of the validation E between neighbouring checkpoints in the training-curve tables, up to about 6 points for the supervised '
            'networks on 200 held-out frames, measure the noise of the choice of checkpoint and threshold); the supervised network trains on 80 % of the development frames (the ISS and Picasso are tuned on all '
            'of them; the labelled-data table shows what fewer frames cost); one architecture and a handful of pilot comparisons, no hyper-parameter search; the mismatch sweep changes one assumption at a time; the simulation of (b) has '
            'no axial dependence beyond the width spread, no long halo and uniform emitter positions (the real emitters lie on three microtubules); the frames of one set are not independent samples of the '
            'structure; the Picasso seconds per frame are those stored by `exp_smlm_iss.py --only P` (a wall-clock on a shared machine).',
            '* Reproduce: `python experiments/exp_smlm_learned.py --run` (needs torch, `SMLM_PICASSO_PYTHON` for the Picasso window stage, `LEARNED_WEIGHTS` outside the repository; stages train / infer / iss / picasso / '
            'collect, resumable, workers share a queue; the weights and caches are not part of the repository) writes `learned_raw.npz`; `python experiments/exp_smlm_learned.py --from-raw` regenerates this file '
            'from it byte for byte (numpy and scipy only). Tests: `tests/test_learned_localizer.py` (the torch tests need torch).',
            '']


# ----------------------------------------------------------------------------------------------------------------------------------
# command line
# ----------------------------------------------------------------------------------------------------------------------------------
def select_models(sets, keys):
    """The models of ``sets`` (in training order), restricted to the 'SET|TAG' strings in ``keys`` if given."""
    allk = model_keys(sets)
    if not keys:
        return allk
    want = {tuple(k.split('|')) for k in keys}
    return [(n, t) for n, t in allk if (n, t) in want]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', action='store_true')
    ap.add_argument('--from-raw', action='store_true')
    ap.add_argument('--stage', nargs='*', default=['prep', 'train', 'infer', 'iss', 'picasso', 'collect'], choices=['prep', 'train', 'infer', 'iss', 'picasso', 'collect'])
    ap.add_argument('--sets', nargs='*', default=list(SETS))
    ap.add_argument('--keys', nargs='*', default=None, help="models 'SET|TAG' (default: all of the sets)")
    ap.add_argument('--iss-conds', nargs='*', default=list(ISS_CONDS))
    ap.add_argument('--hd-extra', action='store_true', help="also run the '/a' conditions of the ISS on the HD set (expensive; skipped by default)")
    ap.add_argument('--scaling', action='store_true', help='also train the data-scaling models of (a) (LD sets)')
    ap.add_argument('--boot', type=int, default=BOOT)
    A = ap.parse_args()
    if A.from_raw:
        raw = load_learned_raw()
        with open(MD, 'w') as f:
            f.write(report(raw, A.boot))
        print('wrote', MD)
        return
    if not A.run:
        ap.error('give --run or --from-raw')
    models = select_models(A.sets, A.keys)
    if A.scaling:
        models += [m_ for m_ in scaling_keys(A.sets) if not A.keys or '|'.join(m_) in A.keys]
    if 'prep' in A.stage:
        for name in A.sets:
            ctx = Ctx(name)
            d0, d1 = ctx.dev
            src = X.CAL_FROM.get(name, name)
            dev_src = frames_photons(src, *PROTOCOL[src]['dev'])
            cal = X.calibrate(dev_src)
            diff = (float(np.abs(cal['bg'] - ctx.bg).max()), abs(cal['a'] - ctx.a))
            print('prep %s: calibration from iss_raw vs recomputed on the development frames of %s: max |bg diff| %.2e photons, |a diff| %.2e; %d dev labels, %.3f emitters per frame, photons mean %.0f; n_train %d' % (
                name, src, diff[0], diff[1], len(ctx.pool), ctx.lam, ctx.pool.mean(), ctx.n_tr), flush=True)
    if 'train' in A.stage:
        for name, tag in models:
            train_model(name, tag, log=lambda s: print(s, flush=True))
    if 'infer' in A.stage:
        for name, tag in models:
            infer_model(name, tag, log=lambda s: print(s, flush=True))
    if 'iss' in A.stage:
        for name in A.sets:
            if name in MM_SETS:
                iss_nominal_check(name, log=lambda s: print(s, flush=True))
                for cond in sorted(A.iss_conds, key=lambda c: (c != 'bg0.5', ISS_CONDS.index(c))):          # the expensive one (noise underestimated) first
                    if cond.endswith('/a') and not name.endswith('LD') and not A.hd_extra:
                        continue
                    iss_condition(name, cond, log=lambda s: print(s, flush=True))
    if 'picasso' in A.stage:
        for name in A.sets:
            if name in MM_SETS:
                for f in (1.0,) + tuple(PIC_FACTORS):
                    picasso_window(name, f, log=lambda s: print(s, flush=True))
    if 'collect' in A.stage:
        collect(A.sets, log=lambda s: print(s, flush=True))


if __name__ == '__main__':
    main()
