"""
Fourier ring correlation (FRC) resolution of a single-molecule localization set
==============================================================================

The resolution of a reconstruction without a ground truth, in the form of Nieuwenhuizen et al. (Nature Methods 10:557, 2013) and Banterle
et al. (J. Struct. Biol. 183:363, 2013): the localizations are split into two independent halves, each half is rendered on a pixel grid, the
two images are correlated in the Fourier domain in rings of equal spatial frequency, and the resolution is the inverse of the frequency at
which the correlation first drops below a threshold. Everything below is a plain, deterministic numpy implementation (no plotting).

Conventions (the places where a number can silently differ)
-----------------------------------------------------------
* Positions are in nm, x along the columns and y along the rows of the image, (0, 0) the outer corner of the first pixel (the convention of
  ``experiments/smlm_common.py``). The field is a SQUARE of ``field_nm`` per side, rendered on ``n = round(field_nm / pixel_nm)`` pixels per
  side (default pixel 10 nm); localizations outside it are ignored (``frc_localizations`` reports how many).
* Rendering is a plain histogram (each localization adds its weight, default 1, to the pixel it falls in), no Gaussian blob: a blob would
  impose a resolution limit of its own.
* Edge effects: both images are multiplied by the same smooth separable Tukey window (``alpha``: the tapered fraction of the width, default
  0.25, cosine tapers) after subtracting the window-weighted mean of each image. A structure cut by the field border would otherwise put the
  same broadband ridge along the two frequency axes into both halves and inflate the correlation at high frequency; the mean subtraction
  removes the leakage of the (strongly correlated) mean image into the lowest rings.
* FRC(q) = sum_ring Re(F1 conj(F2)) / sqrt(sum_ring |F1|^2 sum_ring |F2|^2), the rings being the integer radii k = round(sqrt(kx^2 + ky^2)),
  k = 1 ... n/2, in units of the frequency step 1 / (n pixel_nm) (ring 0, the mean, is dropped; the corners beyond n/2 are not used). The
  frequency of ring k is k / (n pixel_nm) cycles per nm; the highest, k = n/2, is the Nyquist frequency 1 / (2 pixel_nm).
* The curve is smoothed LIGHTLY before the threshold is read: a centred moving average over ``smooth`` rings (default 5, edges repeated).
* Resolution = 1 / f at the FIRST crossing of the smoothed curve below the threshold, f found by linear interpolation in frequency between
  the last ring at or above and the first ring below the threshold. The threshold is 1/7 = 0.142857 (Nieuwenhuizen et al.; the often quoted
  "0.143" is this number rounded: ``frc_resolution(..., threshold=0.143)`` gives the other reading; the two thresholds differ by 1.5e-4 in correlation, which moves the
  crossing by a fraction of a ring, and the experiment prints the largest difference it found). A curve that never drops below the
  threshold (identical halves; structure resolved down to the pixel limit) has ``crossed = False`` and ``resolution_nm = NaN``; the pixel
  limit 2 x pixel_nm is then reported as ``limit_nm``. If the very first ring is already below the threshold the resolution is the field
  size (``at_first_ring = True``): no structure correlates between the halves at all. ``last_nm`` is the same reading at the LAST crossing
  (the highest frequency at which the smoothed curve is still at or above the threshold); a first and a last crossing that differ by much
  flag a curve that wanders around the threshold.
* The resolution of the HALF sets is what is measured (each image holds half of the localizations); that of the full set is somewhat
  better, and no correction for it is made here (Nieuwenhuizen et al. discuss it).

Splitting. ``split`` is ``'random'`` (localizations dealt at random into two halves of equal size, seed), ``'oddeven'`` (odd against even
frame numbers), ``'blocks'`` (blocks of ``block`` consecutive frames dealt at random into the two halves) or ``'alternating'`` (the same blocks
alternately to the two halves: no randomness). A fluorophore that stays on for
several consecutive frames is localized several times, with independent errors but the same true position: under the first two splits
the repeated localizations of one fluorophore fall into BOTH halves and correlate the images at every frequency, which can make the
resolution look better than it is; only a split by blocks of frames keeps (nearly) all repeats in one half. ``bootstrap_blocks`` resamples
the blocks with replacement and splits the UNIQUE resampled blocks between the halves, so that no localization (a block drawn twice
counts twice, with weight 2, within its half) is ever in both halves.
"""
import numpy as np
from scipy.signal.windows import tukey

THRESHOLD = 1.0 / 7.0                  # Nieuwenhuizen et al. 2013 / Banterle et al. 2013 ("0.143")
PIXEL_NM = 10.0                        # default grid of the rendering
ALPHA = 0.25                           # tapered fraction of the Tukey window
SMOOTH = 5                             # rings of the moving average
SPLITS = ('random', 'oddeven', 'blocks', 'alternating')

_CACHE = {}


def grid_size(field_nm, pixel_nm=PIXEL_NM):
    """Pixels per side of the rendering grid: round(field_nm / pixel_nm) (at least 8)."""
    n = int(round(float(field_nm) / float(pixel_nm)))
    if n < 8:
        raise ValueError('the grid has %d pixels per side: field %r nm, pixel %r nm' % (n, field_nm, pixel_nm))
    return n


def _xy(a):
    """(n, 2) float64 positions from an (n, 2+) array (extra columns ignored); an empty input of any shape gives (0, 2)."""
    a = np.asarray(a, dtype=np.float64)
    if a.size == 0:
        return np.empty((0, 2))
    if a.ndim != 2 or a.shape[1] < 2:
        raise ValueError('expected an (n, 2+) array of positions, got shape %r' % (a.shape,))
    return a[:, :2]


def inside_mask(xy_nm, field_nm):
    """Boolean mask of the positions (n, 2+) inside the square field [0, field_nm) x [0, field_nm) (NaN is outside)."""
    xy = _xy(xy_nm)
    with np.errstate(invalid='ignore'):
        return (xy[:, 0] >= 0) & (xy[:, 0] < field_nm) & (xy[:, 1] >= 0) & (xy[:, 1] < field_nm)


def render(xy_nm, field_nm, pixel_nm=PIXEL_NM, weights=None):
    """Histogram image (n, n) of the positions (n_loc, 2+; x along the columns, y along the rows, nm) on the grid of the square field; positions
    outside the field are ignored. ``weights``: one weight per position (default 1). The pixel (i, j) covers x in [j, j + 1) pixel_nm,
    y in [i, i + 1) pixel_nm."""
    n = grid_size(field_nm, pixel_nm)
    xy = _xy(xy_nm)
    if len(xy) == 0:
        return np.zeros((n, n))
    ok = inside_mask(xy, n * pixel_nm)
    ix = np.floor(xy[ok, 0] / pixel_nm).astype(np.int64)
    iy = np.floor(xy[ok, 1] / pixel_nm).astype(np.int64)
    w = None if weights is None else np.asarray(weights, dtype=np.float64)[ok]
    return np.bincount(iy * n + ix, weights=w, minlength=n * n).reshape(n, n).astype(np.float64)


def window2d(n, alpha=ALPHA):
    """The separable Tukey window (n, n): cosine tapers over the fraction ``alpha`` of the width (alpha = 0: rectangular, 1: Hann)."""
    key = ('win', n, float(alpha))
    if key not in _CACHE:
        w = tukey(n, alpha, sym=True) if alpha > 0 else np.ones(n)
        _CACHE[key] = np.outer(w, w)
    return _CACHE[key]


def _rings(n):
    """Ring index, Hermitian weight and mask of the rfft2 half plane of an n x n image (cached)."""
    key = ('rings', n)
    if key not in _CACHE:
        ky = np.fft.fftfreq(n) * n
        kx = np.arange(n // 2 + 1, dtype=np.float64)
        ring = np.rint(np.sqrt(ky[:, None] ** 2 + kx[None, :] ** 2)).astype(np.int64)
        wt = np.full(kx.shape, 2.0)
        wt[0] = 1.0
        if n % 2 == 0:
            wt[-1] = 1.0
        keep = ring <= n // 2
        _CACHE[key] = (ring[keep], np.broadcast_to(wt[None, :], ring.shape)[keep], keep)
    return _CACHE[key]


def frc_curve(img_a, img_b, pixel_nm=PIXEL_NM, alpha=ALPHA, subtract_mean=True):
    """The raw FRC curve of two images of the same square shape: dict(ring (k = 1 .. n/2), freq (cycles per nm), frc, npix (Fourier pixels
    in the ring, counting both halves of the plane)). The images are windowed (Tukey, ``alpha``) after subtracting their window-weighted
    mean (``subtract_mean``); rings with no power in one of the images have NaN."""
    a, b = np.asarray(img_a, dtype=np.float64), np.asarray(img_b, dtype=np.float64)
    if a.ndim != 2 or a.shape != b.shape or a.shape[0] != a.shape[1]:
        raise ValueError('need two square images of the same shape, got %r and %r' % (a.shape, b.shape))
    n = a.shape[0]
    w = window2d(n, alpha)
    if subtract_mean:
        sw = w.sum()
        a = a - (w * a).sum() / sw
        b = b - (w * b).sum() / sw
    fa, fb = np.fft.rfft2(a * w), np.fft.rfft2(b * w)
    ring, wt, keep = _rings(n)
    m = n // 2 + 1
    cross = np.bincount(ring, weights=wt * (fa * np.conj(fb)).real[keep], minlength=m)
    pa = np.bincount(ring, weights=wt * (fa.real ** 2 + fa.imag ** 2)[keep], minlength=m)
    pb = np.bincount(ring, weights=wt * (fb.real ** 2 + fb.imag ** 2)[keep], minlength=m)
    npix = np.bincount(ring, weights=wt, minlength=m)
    with np.errstate(invalid='ignore', divide='ignore'):
        frc = cross / np.sqrt(pa * pb)
    k = np.arange(1, m)
    return dict(ring=k, freq=k / (n * float(pixel_nm)), frc=frc[1:], npix=npix[1:])


def smooth_curve(y, width=SMOOTH):
    """Centred moving average of ``width`` points (odd; 1 or 0: none), the first and last values repeated at the ends."""
    y = np.asarray(y, dtype=np.float64)
    width = int(width)
    if width <= 1:
        return y.copy()
    if width % 2 == 0:
        width += 1
    h = width // 2
    yp = np.concatenate([np.full(h, y[0]), y, np.full(h, y[-1])])
    return np.convolve(yp, np.ones(width) / width, mode='valid')


def frc_resolution(curve, pixel_nm=PIXEL_NM, smooth=SMOOTH, threshold=THRESHOLD):
    """The resolution read from a ``frc_curve`` (see the module docstring for the convention): dict(resolution_nm (first crossing; NaN if none),
    crossed, at_first_ring, last_nm (last crossing; NaN if none), limit_nm (2 x pixel, the best a crossing could show), smoothed (the smoothed
    curve), threshold). A NaN in the smoothed curve is treated as below the threshold."""
    f = np.asarray(curve['freq'], dtype=np.float64)
    ys = smooth_curve(np.nan_to_num(curve['frc'], nan=-1.0), smooth)
    below = ys < threshold
    out = dict(resolution_nm=float('nan'), crossed=False, at_first_ring=False, last_nm=float('nan'), limit_nm=2.0 * pixel_nm, smoothed=ys, threshold=threshold)
    if not np.isfinite(np.asarray(curve['frc'], dtype=np.float64)).any():
        return out                                      # no power in one of the images (an empty half): there is no curve to read
    if below.any():
        k = int(np.argmax(below))
        out['crossed'] = True
        out['resolution_nm'] = 1.0 / _interp(f, ys, k, threshold)
        out['at_first_ring'] = k == 0
        above = np.flatnonzero(~below)
        if above.size:
            j = int(above[-1])                                           # the last ring still at or above the threshold
            out['last_nm'] = 1.0 / _interp(f, ys, j + 1, threshold) if j + 1 < len(f) else float('nan')
        else:
            out['last_nm'] = out['resolution_nm']
    return out


def _interp(f, y, k, thr):
    """The frequency at which the polyline (f, y) falls through ``thr`` between index k - 1 (at or above) and k (below); the first frequency if k = 0."""
    if k == 0:
        return float(f[0])
    y0, y1 = y[k - 1], y[k]
    return float(f[k - 1] + (y0 - thr) / (y0 - y1) * (f[k] - f[k - 1]))


def frc_images(img_a, img_b, pixel_nm=PIXEL_NM, alpha=ALPHA, smooth=SMOOTH, threshold=THRESHOLD):
    """Curve and resolution of two images: the dict of ``frc_curve`` extended with the entries of ``frc_resolution``."""
    c = frc_curve(img_a, img_b, pixel_nm, alpha)
    c.update(frc_resolution(c, pixel_nm, smooth, threshold))
    return c


# ----------------------------------------------------------------------------------------------------------------------------------
# splitting the localizations into two halves
# ----------------------------------------------------------------------------------------------------------------------------------
def split_mask(n, mode='random', seed=0, frames=None, block=10):
    """Boolean array (n,): False = half A, True = half B.

    'random': a random permutation (``numpy.random.default_rng(seed)``) deals floor(n / 2) localizations to B and the rest to A;
    'oddeven': odd frame numbers (``frames``, integers) to B; 'blocks': the blocks frame // block, dealt at random into two halves of
    (nearly) equal numbers of BLOCKS (seed), all frames of a block to the same half; 'alternating': the blocks frame // block alternately to A and B
    (odd block numbers to B), a split without any randomness that balances any slow drift of the sample between the halves."""
    if mode not in SPLITS:
        raise ValueError('unknown split %r (one of %s)' % (mode, ', '.join(SPLITS)))
    if mode == 'random':
        mask = np.zeros(n, bool)
        mask[np.random.default_rng(seed).permutation(n)[: n // 2]] = True
        return mask
    if frames is None or len(frames) != n:
        raise ValueError("split %r needs one frame number per localization" % mode)
    frames = np.asarray(frames, dtype=np.int64)
    if mode == 'oddeven':
        return (frames % 2) == 1
    ids = frames // int(block)
    if mode == 'alternating':
        return (ids % 2) == 1
    uniq = np.unique(ids)
    deal = np.zeros(len(uniq), bool)
    deal[np.random.default_rng(seed).permutation(len(uniq))[: len(uniq) // 2]] = True
    return deal[np.searchsorted(uniq, ids)]


def frc_localizations(xy_nm, frames=None, field_nm=12800.0, pixel_nm=PIXEL_NM, split='random', seed=0, block=10, alpha=ALPHA, smooth=SMOOTH,
                      threshold=THRESHOLD):
    """FRC of a localization set: split it (``split_mask``), render the halves (``render``), correlate (``frc_images``).

    xy_nm: (n, 2+) positions in nm; frames: (n,) integer frame numbers (needed by 'oddeven' and 'blocks'); field_nm: side of the square
    field (default 12800 nm = 128 pixels of 100 nm). Returns the dict of ``frc_images`` plus n_a, n_b (localizations in the halves), n_outside
    (positions ignored because they are outside the field) and the split parameters."""
    xy = _xy(xy_nm)
    fr = None if frames is None else np.asarray(frames)
    mask = split_mask(len(xy), split, seed, fr, block)
    ok = inside_mask(xy, grid_size(field_nm, pixel_nm) * pixel_nm)
    a, b = render(xy[~mask], field_nm, pixel_nm), render(xy[mask], field_nm, pixel_nm)
    out = frc_images(a, b, pixel_nm, alpha, smooth, threshold)
    out.update(n_a=int(((~mask) & ok).sum()), n_b=int((mask & ok).sum()), n_outside=int((~ok).sum()), split=split, seed=seed)
    return out


# ----------------------------------------------------------------------------------------------------------------------------------
# bootstrap over blocks of frames
# ----------------------------------------------------------------------------------------------------------------------------------
def block_ids(frames, block=10, n_frames=None):
    """(ids, uniq, inv): the block number frame // block of every localization, the block numbers (the sorted unique ones present, or with ``n_frames`` all
    blocks 0 .. ceil(n_frames / block) - 1, so that two localization sets of the same frames have the same blocks), and the index of each localization's block in ``uniq``."""
    ids = np.asarray(frames, dtype=np.int64) // int(block)
    if n_frames is not None:
        uniq = np.arange(-(-int(n_frames) // int(block)))
        if len(ids) and (ids.min() < 0 or ids.max() >= len(uniq)):
            raise ValueError('frame numbers outside 0 .. %d' % (n_frames - 1))
        return ids, uniq, ids
    uniq, inv = np.unique(ids, return_inverse=True)
    return ids, uniq, inv


def bootstrap_draw(n_blocks, seed, n_boot):
    """(n_boot, n_blocks) multiplicities: the number of times each block is drawn in each replicate (n_blocks draws with replacement per replicate,
    ``numpy.random.default_rng(seed)``)."""
    rng = np.random.default_rng(seed)
    draws = rng.integers(0, n_blocks, size=(n_boot, n_blocks))
    return np.stack([np.bincount(d, minlength=n_blocks) for d in draws])


def bootstrap_halves(mult, rng):
    """Half assignment (n_blocks,) bool of one replicate: the blocks drawn at least once are dealt at random (the generator ``rng``) into two halves of
    (nearly) equal numbers of blocks; blocks not drawn are False (and have multiplicity 0, so they never enter)."""
    used = np.flatnonzero(mult > 0)
    half = np.zeros(len(mult), bool)
    half[used[rng.permutation(len(used))[: len(used) // 2]]] = True
    return half


def bootstrap_blocks(xy_nm, frames, n_boot=100, seed=0, block=10, field_nm=12800.0, pixel_nm=PIXEL_NM, alpha=ALPHA, smooth=SMOOTH, threshold=THRESHOLD,
                     n_frames=None):
    """Block bootstrap of the FRC resolution: in each of ``n_boot`` replicates the blocks of ``block`` consecutive frames are resampled with
    replacement (as many draws as there are blocks), the UNIQUE drawn blocks are dealt at random into two halves (so that no localization
    is in both; a block drawn k times enters its half with weight k), the halves are rendered with these weights and correlated.

    Returns dict(resolution_nm (n_boot,) (NaN: no crossing), crossed (n_boot,) bool, n_blocks, seed, block). The replicates are a function of
    ``seed`` alone (the multiplicities come from ``default_rng(seed)``, the half assignment of replicate r from ``default_rng([seed, r])``); with ``n_frames``
    the blocks are those of the whole sequence, so two localization sets of the same frames get IDENTICAL resamples and halves (paired comparisons).

    A replicate holds only ~63 % of the distinct blocks, and the resolution of halves with fewer localizations is coarser: the replicates are systematically
    coarser than the estimate on the full sample and their percentile interval does not bracket it. Their RELATIVE spread is what carries over (and, in a pair of
    methods resampled together, the ratio of their resolutions)."""
    xy = _xy(xy_nm)
    _, uniq, inv = block_ids(frames, block, n_frames)
    mult = bootstrap_draw(len(uniq), seed, n_boot)
    res = np.full(n_boot, np.nan)
    crossed = np.zeros(n_boot, bool)
    for r in range(n_boot):
        half = bootstrap_halves(mult[r], np.random.default_rng([seed, r]))
        w = mult[r][inv].astype(np.float64)
        in_b = half[inv]
        a, b = render(xy[~in_b], field_nm, pixel_nm, w[~in_b]), render(xy[in_b], field_nm, pixel_nm, w[in_b])
        c = frc_images(a, b, pixel_nm, alpha, smooth, threshold)
        res[r], crossed[r] = c['resolution_nm'], c['crossed']
    return dict(resolution_nm=res, crossed=crossed, n_blocks=len(uniq), seed=seed, block=block)
