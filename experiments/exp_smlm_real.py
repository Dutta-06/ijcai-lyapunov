"""
A REAL experimental 2D SMLM dataset (EPFL 2013 challenge, Tubulin Alexa647, "Real High Density"): resolution by Fourier ring correlation,
sliding inverse scale space against Picasso
==================================================================================================================================

What this experiment measures, and what it cannot
-------------------------------------------------
The 500 frames (128 x 128 px, 16 bit) of the real High Density sequence of the 2013 challenge (a fixed cell, anti-alpha-tubulin / Alexa647; data
from the hub https://srm.epfl.ch/srm/dataset/challenge-2D-real/Real_High_Density/index.html, fetched by ``data/smlm/fetch_real.py``) are localized
by (1) Picasso with its DEFAULT pipeline (the R0 settings of ``experiments/exp_smlm_picasso.py``) and (2) the sliding inverse scale space (ISS) of
``src/core/imaging.py`` with an A PRIORI configuration, and the two localization sets are compared by
  * the Fourier ring correlation (FRC, ``experiments/frc.py``) resolution of the pooled reconstruction, with a block-bootstrap interval over the frames
    and for several splits of the localizations,
  * one ground-truth-free precision proxy: the repeat-localization precision, the scale of the same-molecule component of the nearest-neighbour
    distances between localizations in consecutive frames (NeNA-like), see ``repeat_precision``,
  * the numbers of localizations and the seconds per frame.
THERE IS NO GROUND TRUTH for these frames (the page says "Ground-truth localization unknown"). The experiment therefore measures resolution and
precision PROXIES only. It does NOT measure detection accuracy: how many of the localizations of a method lie on the same structure as the other
method's, or on a real filament at all, is not knowable here. An FRC resolution rewards localizations that are mutually consistent between the
halves of the data; spurious localizations that are consistent (the same bright blob split in the same way in two halves) are rewarded as well, and
a method that localizes fewer, better molecules can look worse. The repeat precision describes the molecules that are localized in two consecutive
frames (brighter than average), not all.

What the pages give and what is assumed (see also the docstring of ``data/smlm/fetch_real.py``)
-----------------------------------------------------------------------------------------------
READ ON THE PAGES: the sample, the number of frames (500 for the High Density sequence), that the frames are 16 bit TIFF, the courtesy, "ground truth
unknown". READ FROM THE FILES: 500 frames of 128 x 128 px, one 16-bit TIFF per frame. NOT GIVEN ANYWHERE and therefore ASSUMED: the PIXEL SIZE, taken
as 100 nm (``PIXEL_NM``) so that positions can be expressed in nm and the FRC rendered on a 10 nm grid; every length in nm below scales with it. Two
second-hand facts are consistent with it and do not prove it: the density arithmetic of Nehme et al., arXiv:1801.09631 (their "experimentally measured microtubules"
dataset of 500 frames with about 520,000 emitters at 6.35 emitters per square micrometre is 163.8 square micrometres per frame = (12.8 um)^2 = 128 px of 100 nm; the arXiv
text does not say that this is the very sequence used here, and the arithmetic holds only if it is), and the fitted PSF width of 1.4-1.6 px = 140-160 nm of Alexa647
imaged with a high-NA objective (the ideal 0.21 lambda / NA is ~95 nm; the simulated sets of this repository, at 100 nm pixels, fit 140 nm). The camera (model, gain, e-/count, baseline) is not given either: Picasso gets the
harness defaults of ``smlm_common`` (baseline 100, EM gain 300, 45 e-/ADU, QE 0.9, which are those of the simulated sets, NOT of this camera; they
change only Picasso's photon column, not its identification or its positions), and the ISS works in raw ADU with the noise law estimated from the data.

A priori configuration of the ISS (DECLARED BEFORE any FRC of a real localization set was computed; nothing below was tuned on these data)
------------------------------------------------------------------------------------------------------------------------------------
Order of events, stated for the record: the frames were looked at (mean image, background level and its drift, pixel statistics), Picasso R0 and
a low-threshold Picasso run were counted and their fitted widths looked at, and the ISS was timed on two frames, all BEFORE the rules below were
written down and before any FRC was computed on real localizations. FRC was first computed on synthetic data (``tests/test_frc.py``). ONE REVISION
of the noise-law rule was made after that, on label-free diagnostics of two frames (never the FRC): the first version took the scale of the noise law from
the single-frame high-frequency statistics alone; the temporal variance of the darkest 10 % of the pixels came out 1.36 times larger (1.23 times for the darkest 5 %
that the final rule uses), and on frame 250 the ISS with the single-frame scale put 180 units (115 with a variance 1.36 times larger; Picasso R0 79) with a
weighted residual chi2 of 3.2, splitting bright blobs. The scale of the law is therefore anchored to the temporal variance of the dark pixels (the quantity of
the simulated sets' calibration), as below; the final scale is 1.23. Nothing else was changed in the ISS configuration.
 * Model and algorithm: ``imaging.fit_frame`` with the configuration 'hp5h' of ``exp_smlm_iss.CONFIGS`` (Gaussian core 0.75 + halo 0.25 with 3x the
   width, high-pass 5 px, pixel-integrated kernels, widths sliding per kernel), start lattice h = 0.5 px, stop rule ('snr', kappa = 5) (``exp_smlm_iss.DEFAULT``),
   no post-filter on the fitted width, no position offset correction. The same constants as the simulated sets except where a length is tied to the PSF
   width sigma0: the entry width sigma0, the bounds of the sliding widths (1.0, 3.2) px x sigma0 / 1.4 and the resolution limit r_ex = 2.0 px x sigma0 / 1.4
   (all three equal the simulated sets' values for sigma0 = 1.4 px).
 * sigma0 (label free): the median of sqrt(sx sy) (px) of the fitted widths of Picasso-R0 localizations that are BRIGHT (photons at or above the median of
   all R0 localizations) and ISOLATED (no other R0 localization of the same frame within 8 px), see ``sigma0_from_picasso``. The frame is not circular
   (the median sx and sy differ by ~25 %): the ISS kernel is circular and has no ellipticity parameter, so this is a model mismatch the experiment lives with.
 * Units and noise: raw ADU. The frames are shifted by a reference dark level x_ref (the mean of the darkest 1 % of the smoothed pixel values, so that
   the shifted frame is non-negative; it is a reference of the regression, NOT the camera baseline). Background: the robust quadratic surface fitted to the
   per-pixel temporal median (``exp_smlm_iss.fit_surface``, unchanged) plus, per frame, the median of the frame minus that surface (the frame's own offset: the
   background level falls from ~440 to ~290 ADU over the sequence). Noise law var = a x + b0 (x in ADU above x_ref), ``calibrate_noise``: its SHAPE from the
   heteroscedastic high-frequency noise of the SINGLE frames (``noise_law``: the high-pass residual, spatial frequency above 0.45 cycles per pixel where the
   PSF-blurred signal has no power, of every 5th frame, mirror-extended so that no frame-border ridge enters; its robust variance (MAD) in 16 quantile bins of
   the local mean (the frame smoothed by a Gaussian of 2 px); a straight line fitted to the bin variances by least squares in RELATIVE terms); its absolute SCALE from the TEMPORAL
   variance of the darkest 5 % of the pixels (per-pixel robust variance of the time series after removing each frame's median, median over the pixels;
   ``temporal_dark_variance``), by which the whole single-frame line is multiplied. The temporal variance of a pixel at higher levels, which the simulated sets
   used for the slope, is NOT usable here: in this dense sequence it is dominated by blinking (see the table of temporal variance by level in the results).
 * Crowded-field recipe (a computational necessity, the one of ``exp_astro_m2.py``, not a tuned value): refits of the neighbours of an entrant only, within
   4 px (``fit_radius``, instead of the transitive cluster, which in a 128 x 128 px frame with 100-250 kernels is the whole frame), in a window of 10 px beyond
   them; no final global polish of all kernels against the whole frame (``refine_box = 0``: its cost grows with the cube of the number of kernels); a cap of
   400 entries per frame (the 40 of the 64 x 64 px simulated frames, scaled by the area and then by ~2.5; it is a safety bound: the experiment reports in how
   many frames it binds). Everything else as ``fit_frame``'s defaults (tol 1e-8, three entries of tail after the first entrant below significance 3).
 * Picasso: ``smlm_common.run_picasso`` with ``PICASSO_DEFAULTS`` (min net gradient 5000, box 7, ``gausslq``, no padding, no post-filter) = rung R0.
 * FRC: 10 nm pixels over the 12.8 um field, Tukey window (0.25), mean subtracted, 5-ring moving average, threshold 1/7 (``frc.py``). Headline split: BLOCKS
   of 10 consecutive frames dealt alternately to the two halves (the only kind of split that does not put the repeated localizations of one fluorophore into both
   halves; no randomness); also reported: 10 random assignments of the blocks (seeds 0-9), random localization splits (seeds 0-9) and odd / even frames.
   Uncertainty: the block bootstrap as specified (``frc.bootstrap_blocks``: blocks of 10 frames resampled with replacement, the unique blocks dealt into the two halves,
   100 replicates) whose replicates turned out to be systematically coarser than the full-sample estimate (a replicate holds ~63 % of the distinct blocks), so
   that its percentile interval does not bracket the estimate; its relative spread (a "recentred" interval) and the PAIRED ratio of two methods resampled together are used
   instead. This reporting was changed after the first full run (the first headline was the median over the random block assignments and the percentile
   interval of the bootstrap was shown as it was); the localization sets, every calibration and the splits were not touched, and all the rows of the first report are still there.
 * Rows named "variation" are additional, clearly labelled and not selected by their result: the ISS units whose fitted width is not at a bound of the
   sliding widths (a post-filter on the fitted width, the kind Picasso's width window is; about a quarter of the units sit at the lower bound), and Picasso with its minimum
   net gradient lowered to 1000 (the smallest value of the baseline ladder's grid).

Usage:  python experiments/exp_smlm_real.py --run [--frames N] [--workers 2] [--cache-dir DIR]     (localizes, writes results/smlm/real_raw.npz and real_frc.md)
        python experiments/exp_smlm_real.py --bench 50                                             (adds the single-process CPU timing of the ISS to the raw, regenerates the md)
        python experiments/exp_smlm_real.py --calib-check                                          (adds the stability of the noise-law scale along the sequence to the raw)
        python experiments/exp_smlm_real.py --from-raw [--boot 100]                                (regenerates real_frc.md from the raw alone, byte for byte)
Environment: SMLM_DATA_REAL (default /home/user/smlm_data_real), SMLM_PICASSO_PYTHON (an interpreter with picasso), REAL_OUT (output directory).
"""
import argparse
import hashlib
import io
import json
import os
import platform
import sys
import time
import zipfile

os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import numpy as np                                    # noqa: E402
from scipy import ndimage as ndi                      # noqa: E402
from scipy.optimize import minimize                   # noqa: E402
from scipy.spatial import cKDTree                     # noqa: E402

HERE = os.path.abspath(os.path.dirname(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
sys.path.insert(0, REPO)
import smlm_common as C                               # noqa: E402
import frc as FRC                                     # noqa: E402
import exp_smlm_iss as ISS                            # noqa: E402
from src.core import imaging as IM                    # noqa: E402

DATA = os.environ.get('SMLM_DATA_REAL', '/home/user/smlm_data_real')
OUT = os.environ.get('REAL_OUT', os.path.join(REPO, 'results', 'smlm'))
RAW = os.path.join(OUT, 'real_raw.npz')
MD = os.path.join(OUT, 'real_frc.md')
SET = 'Real_High_Density'
PIXEL_NM = 100.0                  # ASSUMED (not on the pages, see the docstring)

# ---- the a priori ISS configuration (declared in the docstring) -------------------------------------------------------------------------
CAP = 400                         # entries per frame at most
FIT_RADIUS = 4.0                  # px: refit only the kernels within this distance of an entrant (exp_astro_m2.FIT_RADIUS)
WINDOW = 10.0                     # px: window of a local refit beyond the refitted kernels (exp_astro_m2: rtrunc = 10 / reach)
REFINE_BOX = 0.0                  # no final global polish
TOL = 1e-8
MARGIN = 1.0
ISO_PX = 8.0                      # isolation radius of the sigma0 rule (px), as in results/smlm/calibration.md
MIN_SIGMA0_SAMPLE = 50            # fewer qualifying localizations: fall back to ISS.SIGMA0 (flagged)
# ---- the noise law -----------------------------------------------------------------------------------------------------------------------
NL_STEP, NL_F_LO, NL_F_HI, NL_SIGMA_MU, NL_BORDER = 5, 0.35, 0.45, 2.0, 4
TD_QUANTILE = 0.05                # the darkest 5 % of the pixels (by temporal median) anchor the scale of the noise law
TD_MIN_FRAMES = 50                # fewer frames: the temporal variance is not reliable and the single-frame law is used as it is (scale 1, flagged)
NL_QUANTILES = (0.0, 0.01, 0.025, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 0.98, 0.995, 1.0)
# ---- the analysis ------------------------------------------------------------------------------------------------------------------------
FRC_PIXEL_NM = 10.0
BLOCK = 10                        # frames per block (bootstrap and block split)
N_SPLITS = 10                     # split seeds 0 .. 9
N_BOOT = 100
BOOT_SEED = 0
D_MAX_NM = 150.0                  # repeat-localization window
WIDTH_TOL = 0.01                  # px: a fitted width within this of a bound of the sliding widths counts as at the bound
PIC_LOW_NG = 1000.0               # variation: Picasso's minimum net gradient
THIN_SEEDS = (0, 1, 2, 3, 4)      # random thinnings of a set to the count of another (count-controlled comparison)
ZOOM_NM, ZOOM_STEP_NM = 3000.0, 100.0
AGREE_NM = 150.0
FREQS_UM = (5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0)     # frequencies (cycles per micrometre) of the printed FRC values

ROWS = [                          # (key, label, kind, main)
    ('pic', 'Picasso, default pipeline (R0)', 'pic', True),
    ('iss', 'sliding ISS, a priori configuration', 'iss', True),
    ('pic_ng1000', 'variation: Picasso, minimum net gradient 1000', 'pic', False),
    ('iss_interior', 'variation: ISS, units whose fitted width is not at a bound of the sliding widths', 'iss', False),
]


# ----------------------------------------------------------------------------------------------------------------------------------
# data
# ----------------------------------------------------------------------------------------------------------------------------------
def data_path(name=SET, root=None):
    return os.path.join(root or DATA, name, 'sequence.zip')


def load_real_frames(path, start=0, stop=None):
    """Frames ``start`` .. ``stop - 1`` (0-based positions in the sequence, stop exclusive; None: to the end) of a zip of per-frame 16-bit TIFFs
    (sequence/00001.tif, ...; read in memory with Pillow, nothing is extracted) as float64 (n, ny, nx) in ADU. The member numbers must be
    consecutive; frame 0 is the first member. Raises IndexError when the sequence has fewer than ``stop`` frames (no silent truncation)."""
    from PIL import Image
    with zipfile.ZipFile(path) as z:
        members = []
        for info in z.infolist():
            stem, ext = os.path.splitext(os.path.basename(info.filename))
            if info.file_size > 0 and ext.lower() in ('.tif', '.tiff') and stem.isdigit():
                members.append((int(stem), info))
        if not members:
            raise ValueError('no per-frame TIFF members in %s' % path)
        members.sort(key=lambda t: t[0])
        nums = [m[0] for m in members]
        if nums != list(range(nums[0], nums[0] + len(nums))):
            raise ValueError('the member numbers of %s are not consecutive' % path)
        if start < 0 or (stop is not None and stop < start):
            raise ValueError('need 0 <= start <= stop, got start=%r stop=%r' % (start, stop))
        if stop is not None and stop > len(members):
            raise IndexError('the sequence has only %d frames, requested stop=%d' % (len(members), stop))
        out = []
        for _, info in members[start:stop]:
            a = np.asarray(Image.open(io.BytesIO(z.read(info))))
            if a.ndim != 2:
                raise ValueError('expected single-channel frames, got shape %r' % (a.shape,))
            out.append(a.astype(np.float64))
    if not out:
        return np.empty((0, 0, 0))
    if len({a.shape for a in out}) != 1:
        raise ValueError('frames of different sizes in %s' % path)
    return np.stack(out)


def file_sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


# ----------------------------------------------------------------------------------------------------------------------------------
# calibration (no labels exist; nothing is tuned)
# ----------------------------------------------------------------------------------------------------------------------------------
def _highpass_mask(ny, nx, f_lo, f_hi):
    """Isotropic raised-cosine high-pass in the Fourier domain: 0 below f_lo, 1 above f_hi cycles per pixel."""
    f = np.hypot(np.fft.fftfreq(nx)[None, :], np.fft.fftfreq(ny)[:, None])
    t = np.clip((f - f_lo) / (f_hi - f_lo), 0.0, 1.0)
    return 0.5 - 0.5 * np.cos(np.pi * t)


def _mirror(a):
    """Even extension of a frame to (2 ny, 2 nx): continuous and periodic, so a Fourier filter has no wrap-around edge."""
    return np.block([[a, a[:, ::-1]], [a[::-1, :], a[::-1, ::-1]]])


def noise_law(frames, step=NL_STEP, f_lo=NL_F_LO, f_hi=NL_F_HI, sigma_mu=NL_SIGMA_MU, border=NL_BORDER, quantiles=NL_QUANTILES, min_bin=200):
    """The heteroscedastic noise law var = a x + b0 of the pixel values from the high-frequency content of SINGLE frames (frames (n, ny, nx) in ADU).

    Every ``step``-th frame is mirror-extended and filtered with a raised-cosine high-pass (0 below ``f_lo``, 1 above ``f_hi`` cycles per pixel, where
    the signal of a PSF of ~1.5 px width has no power), scaled so that white noise keeps its variance (divided by the root of the mean squared mask),
    the border of ``border`` px is dropped, and every pixel is paired with the local mean (the frame smoothed by a Gaussian of ``sigma_mu`` px). The pairs
    of all frames are cut into bins at the ``quantiles`` of the local mean; per bin the variance is (1.4826 MAD)^2 of the residuals (robust against
    the part of the signal that leaks through the filter), its abscissa the mean of the local mean. The line v = a x + c is fitted to the (x, v) of the bins
    by least squares in RELATIVE terms (minimising the sum of ((v - a x - c) / v)^2). x_ref is the abscissa of the darkest bin and b0 = a x_ref + c the
    variance there, so that var = a (x - x_ref) + b0 for the frame shifted by x_ref. Returns dict(a, c, x_ref, b0, x_bins, v_bins, n_bins, n_frames, n_pixels)."""
    frames = np.asarray(frames, dtype=np.float64)
    sel = frames[::max(1, int(step))]
    ny, nx = sel.shape[1:]
    if min(ny, nx) <= 2 * border + 8:
        raise ValueError('frames of %d x %d px are too small for a border of %d px' % (ny, nx, border))
    mask = _highpass_mask(2 * ny, 2 * nx, f_lo, f_hi)
    gain = float((mask ** 2).mean())
    mus, rs = [], []
    for f in sel:
        r = np.fft.ifft2(np.fft.fft2(_mirror(f)) * mask).real[:ny, :nx] / np.sqrt(gain)
        mu = ndi.gaussian_filter(f, sigma_mu, mode='reflect')
        sl = (slice(border, ny - border), slice(border, nx - border))
        mus.append(mu[sl].ravel())
        rs.append(r[sl].ravel())
    mu, r = np.concatenate(mus), np.concatenate(rs)
    edges = np.quantile(mu, quantiles)
    xs, vs, ns = [], [], []
    for i in range(len(edges) - 1):
        last = i == len(edges) - 2
        m = (mu >= edges[i]) & ((mu <= edges[i + 1]) if last else (mu < edges[i + 1]))
        if m.sum() < min_bin:
            continue
        rr = r[m]
        v = (1.4826 * np.median(np.abs(rr - np.median(rr)))) ** 2
        xs.append(float(mu[m].mean()))
        vs.append(float(v))
        ns.append(int(m.sum()))
    xs, vs = np.array(xs), np.array(vs)
    if len(xs) < 4 or not (vs > 0).all():
        raise RuntimeError('the noise law cannot be estimated: %d usable bins' % len(xs))
    wgt = 1.0 / vs
    a, c = np.linalg.lstsq(np.c_[xs, np.ones_like(xs)] * wgt[:, None], vs * wgt, rcond=None)[0]
    x_ref = float(xs[0])
    b0 = float(a * x_ref + c)
    if not (a > 0 and b0 > 0):
        raise RuntimeError('the noise law is not a positive line: a = %g, variance at the dark level = %g' % (a, b0))
    return dict(a=float(a), c=float(c), x_ref=x_ref, b0=b0, x_bins=xs, v_bins=vs, n_bins=np.array(ns), n_frames=len(sel), n_pixels=int(len(mu)))


def temporal_dark_variance(frames, quantile=TD_QUANTILE):
    """The temporal variance of the DARKEST pixels, where blinking does not contaminate the temporal statistics (the quantity ``exp_smlm_iss.calibrate`` measures, restricted
    to the pixels that carry no structure): the pixels whose temporal median is in the lowest ``quantile``, each frame's median removed, the robust variance
    (1.4826 MAD)^2 of every pixel's time series, its median over the pixels. Returns dict(v (ADU^2), x (the mean temporal median of those pixels, ADU), n (pixels))."""
    frames = np.asarray(frames, dtype=np.float64)
    n = len(frames)
    med = np.median(frames, axis=0)
    dark = med <= np.quantile(med, quantile)
    d = frames - np.median(frames.reshape(n, -1), axis=1)[:, None, None]
    sd = d[:, dark]
    mad = 1.4826 * np.median(np.abs(sd - np.median(sd, axis=0)), axis=0)
    return dict(v=float(np.median(mad ** 2)), x=float(med[dark].mean()), n=int(dark.sum()))


def calibrate_noise(frames):
    """The noise law used by the ISS on the frames (n, ny, nx) in ADU: the SHAPE from the single-frame high-frequency statistics (``noise_law``: independent pixel
    noise, a line in the local mean), the absolute SCALE from the temporal variance of the darkest pixels (``temporal_dark_variance``): the whole line is multiplied by
    s = (temporal variance of the dark pixels) / (the high-frequency law at their level). The single-frame law sees the white part of the noise only; the temporal
    variance of a pixel is its whole fluctuation in time (read-out and shot noise plus whatever else fluctuates there), the variance the fit has to live with. For purely white
    noise s = 1 (and with fewer than ``TD_MIN_FRAMES`` frames s = 1 by rule). Returns the dict of ``noise_law`` with the final a, b0, c (var = a (x - x_ref) + b0 = a x + c) and a_hf, b0_hf, scale, and the dark variance (v_dark, x_dark, n_dark)."""
    law = noise_law(frames)
    td = temporal_dark_variance(frames)
    v_hf = law['a'] * (td['x'] - law['x_ref']) + law['b0']
    scale = td['v'] / v_hf
    if len(frames) < TD_MIN_FRAMES:
        scale = 1.0                                     # too few frames for a temporal variance: the single-frame law as it is (flagged by n_frames_temporal)
    if not (np.isfinite(scale) and scale > 0):
        raise RuntimeError('the temporal dark variance (%g) and the single-frame law (%g) cannot be reconciled' % (td['v'], v_hf))
    out = dict(law)
    out.update(a_hf=law['a'], b0_hf=law['b0'], c_hf=law['c'], scale=float(scale), a=float(scale * law['a']), b0=float(scale * law['b0']), v_dark=td['v'], x_dark=td['x'], n_dark=td['n'])
    out['c'] = out['b0'] - out['a'] * out['x_ref']
    return out


def scale_by_segment(frames, seg=100):
    """The scale s of ``calibrate_noise`` (and its ingredients) on consecutive segments of ``seg`` frames: how stable the absolute scale of the noise law is along the sequence.
    Returns rows [first frame, stop frame, scale, temporal dark variance, level of the dark pixels]; the trailing frames that do not fill a segment are not used."""
    frames = np.asarray(frames, dtype=np.float64)
    rows = []
    for a in range(0, len(frames) - seg + 1, seg):
        law = calibrate_noise(frames[a:a + seg])
        rows.append([a, a + seg, law['scale'], law['v_dark'], law['x_dark']])
    return np.array(rows, dtype=np.float64).reshape(-1, 5)


def temporal_variance_by_level(frames, edges):
    """The temporal variance of the pixels in bins of their temporal median, each frame's median removed: the median over the pixels of the bin of (1.4826 MAD)^2
    of their time series. Shows what the simulated sets' calibration would have seen (blinking dominates above the dark level). Returns [(level, variance, n)]."""
    frames = np.asarray(frames, dtype=np.float64)
    n = len(frames)
    med = np.median(frames, axis=0)
    d = frames - np.median(frames.reshape(n, -1), axis=1)[:, None, None]
    mad = 1.4826 * np.median(np.abs(d - np.median(d, axis=0)), axis=0)
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (med >= lo) & (med < hi)
        if m.sum() >= 5:
            out.append((float(np.median(med[m])), float(np.median(mad[m] ** 2)), int(m.sum())))
    return out


def sigma0_from_picasso(dets, pixel_nm=PIXEL_NM, iso_px=ISO_PX):
    """sigma0 (px) of the label-free rule: the median of sqrt(sx sy) over the BRIGHT (photons >= the median of all localizations) and ISOLATED (no other
    localization of the same frame within ``iso_px`` px) Picasso localizations. dets: per-frame arrays with the columns of ``smlm_common.DET_COLUMNS``
    (x_nm, y_nm, photons, sx_nm, sy_nm, net_gradient). Returns dict(sigma0, n_sel, n_all, n_iso, n_bright, photon_threshold, median_sx, median_sy, median_all, fallback)."""
    xs, ph, sx, sy, iso = [], [], [], [], []
    for d in dets:
        d = np.asarray(d, dtype=np.float64)
        if d.size == 0:
            continue
        ok = np.isfinite(d[:, :5]).all(axis=1)
        d = d[ok]
        if len(d) == 0:
            continue
        xy = d[:, :2] / pixel_nm
        near = np.full(len(d), np.inf)
        if len(d) > 1:
            near = cKDTree(xy).query(xy, k=2)[0][:, 1]
        iso.append(near > iso_px)
        ph.append(d[:, C.PH_COL])
        sx.append(np.abs(d[:, C.SX_COL]) / pixel_nm)
        sy.append(np.abs(d[:, C.SY_COL]) / pixel_nm)
    if not ph:
        return dict(sigma0=ISS.SIGMA0, n_sel=0, n_all=0, n_iso=0, n_bright=0, photon_threshold=float('nan'), median_sx=float('nan'), median_sy=float('nan'),
                    median_all=float('nan'), fallback=True)
    iso, ph, sx, sy = np.concatenate(iso), np.concatenate(ph), np.concatenate(sx), np.concatenate(sy)
    thr = float(np.median(ph))
    bright = ph >= thr
    sel = iso & bright
    s = np.sqrt(sx * sy)
    fallback = int(sel.sum()) < MIN_SIGMA0_SAMPLE
    return dict(sigma0=float(ISS.SIGMA0 if fallback else np.median(s[sel])), n_sel=int(sel.sum()), n_all=int(len(ph)), n_iso=int(iso.sum()), n_bright=int(bright.sum()),
                photon_threshold=thr, median_sx=float(np.median(sx[sel])) if sel.any() else float('nan'), median_sy=float(np.median(sy[sel])) if sel.any() else float('nan'),
                median_all=float(np.median(s)), fallback=bool(fallback))


def calibrate(frames, pic_dets):
    """The label-free calibration of the ISS on the frames (n, ny, nx) in ADU and Picasso-R0 detections: background surface, noise law, sigma0."""
    bg = ISS.fit_surface(np.median(frames, axis=0))
    law = calibrate_noise(frames)
    sig = sigma0_from_picasso(pic_dets)
    return dict(bg=bg, law=law, sig=sig)


def iss_config(sigma0):
    """The a priori configuration of the ISS for a PSF width sigma0 (px): the constants of exp_smlm_iss, the widths and the resolution limit scaled by sigma0 / 1.4."""
    base = ISS.CONFIGS[ISS.DEFAULT['cfg']]
    s = float(sigma0) / ISS.SIGMA0
    return dict(cfg=ISS.DEFAULT['cfg'], hp=float(base['hp']), parts=[list(p) for p in base['parts']], sigma0=float(sigma0), widths=[ISS.WIDTHS[0] * s, ISS.WIDTHS[1] * s],
                r_ex=ISS.R_EX * s, lat_h=ISS.LAT_H, margin=MARGIN, kappa=float(ISS.DEFAULT['kappa']), cap=CAP, snr_extra=ISS.SNR_EXTRA, fit_radius=FIT_RADIUS,
                window=WINDOW, refine_box=REFINE_BOX, tol=TOL)


# ----------------------------------------------------------------------------------------------------------------------------------
# the ISS on frames
# ----------------------------------------------------------------------------------------------------------------------------------
def make_context(bg, a, b0, x_ref, cfg, shape):
    """The objects one worker needs: the PSF (with high-pass), the lattice and the calibration in the shifted units (frame - x_ref)."""
    ny, nx = shape
    P = IM.PixelPSF(ny, nx, cfg['sigma0'], parts=[tuple(p) for p in cfg['parts']], hp=cfg['hp'])
    L = IM.Lattice(P, h=cfg['lat_h'], margin=cfg['margin'])
    reach = P.reach(cfg['widths'][1])
    return dict(P=P, L=L, bg=np.asarray(bg, dtype=np.float64) - x_ref, a=float(a), b0=float(b0), x_ref=float(x_ref), cfg=cfg, rtrunc=cfg['window'] / reach)


def process_frame(frame, ctx):
    """One frame (ADU, float): the localizations of the sliding ISS (``imaging.fit_frame``) in the a priori configuration. Returns dict(det (k, 5): x_px, y_px
    (pixel-index coordinates, the pixel centre an integer), sigma_px (the fitted width), w_adu (the fitted integrated counts above the background), own (the
    own significance of the unit, ``imaging.unit_significance``); n_path (entries of the path), hit_cap (the path ended at the cap with every entrant
    significant), secs, delta (the frame's offset against the background surface, ADU), chi2 (the weighted mean squared residual of the high-passed frame after
    the chosen block: about 1 for a complete model and the right noise law))"""
    cfg, P = ctx['cfg'], ctx['P']
    y = np.asarray(frame, dtype=np.float64) - ctx['x_ref']
    delta = float(np.median(y - ctx['bg']))
    bgv = ctx['bg'] + delta
    t0 = time.perf_counter()
    out = IM.fit_frame(y, P, ctx['L'], bgv, ctx['a'], ctx['b0'], widths=tuple(cfg['widths']), stop=('snr', cfg['kappa']), max_entries=cfg['cap'], r_ex=cfg['r_ex'],
                       fit_radius=cfg['fit_radius'], rtrunc=ctx['rtrunc'], refine_box=cfg['refine_box'], snr_extra=cfg['snr_extra'], tol=cfg['tol'])
    c, w = out['c'], out['w']
    if len(w):
        own = IM.unit_significance(P, 1.0 / np.sqrt(ctx['a'] * bgv + ctx['b0']), c, w)
        det = np.c_[c[:, 0], c[:, 1], c[:, 2], w, own]
    else:
        det = np.zeros((0, 5))
    snr = out['snr']
    hit = bool(len(snr) >= cfg['cap'] and np.all(snr >= cfg['kappa']))
    chi2 = float(out['hist'][out['index']]['chi2'])
    return dict(det=det, n_path=int(out['n_entries']), hit_cap=hit, secs=time.perf_counter() - t0, delta=delta, chi2=chi2)


_CTX = {}


def _pool_init(args):
    _CTX.clear()
    _CTX.update(make_context(*args))


def _pool_work(job):
    k, frame = job
    return k, process_frame(frame, _CTX)


def _chunk_signature(frames, cal_args, cfg):
    h = hashlib.sha256()
    h.update(json.dumps(cfg, sort_keys=True).encode())
    for a in cal_args:
        h.update(np.ascontiguousarray(a, dtype=np.float64).tobytes())
    h.update(np.ascontiguousarray(frames, dtype=np.float64).tobytes())
    return h.hexdigest()[:16]


def _pack_results(res):
    """list of process_frame results -> (offsets, flat (M, 5), stats (F, 5: n_path, hit_cap, secs, delta, chi2))."""
    off = np.concatenate([[0], np.cumsum([len(r['det']) for r in res])]).astype(np.int64)
    flat = np.concatenate([r['det'] for r in res]) if off[-1] else np.zeros((0, 5))
    st = np.array([[r['n_path'], float(r['hit_cap']), r['secs'], r['delta'], r['chi2']] for r in res], dtype=np.float64).reshape(-1, 5)
    return off, flat, st


def _unpack_results(off, flat, st):
    return [dict(det=flat[off[i]:off[i + 1]], n_path=int(st[i, 0]), hit_cap=bool(st[i, 1]), secs=float(st[i, 2]), delta=float(st[i, 3]), chi2=float(st[i, 4])) for i in range(len(st))]


def run_iss(frames, cal, cfg, workers=2, cache_dir=None, chunk=25, log=print):
    """The ISS on every frame (n, ny, nx) in ADU. cal: dict(bg, a, b0, x_ref). Frames are processed in chunks of ``chunk`` in a pool of ``workers`` processes
    (1: in this process); with ``cache_dir`` a finished chunk is stored there under a signature of its input and the configuration and is read back by a later run
    (a crash does not lose the finished chunks). Returns the list of ``process_frame`` results in frame order."""
    args = (cal['bg'], cal['a'], cal['b0'], cal['x_ref'], cfg, frames.shape[1:])
    cal_args = (cal['bg'], [cal['a'], cal['b0'], cal['x_ref']])
    out = [None] * len(frames)
    pool = None
    if workers > 1:
        import multiprocessing
        pool = multiprocessing.get_context('fork').Pool(workers, initializer=_pool_init, initargs=(args,))
    else:
        _pool_init(args)
    t_start = time.time()
    try:
        for a in range(0, len(frames), chunk):
            b = min(a + chunk, len(frames))
            path = None
            if cache_dir:
                os.makedirs(cache_dir, exist_ok=True)
                path = os.path.join(cache_dir, 'chunk_%05d_%05d_%s.npz' % (a, b, _chunk_signature(frames[a:b], cal_args, cfg)))
                if os.path.exists(path):
                    with np.load(path, allow_pickle=False) as z:
                        out[a:b] = _unpack_results(z['off'], z['flat'], z['st'])
                    log('  frames %d-%d: read from the cache' % (a, b - 1))
                    continue
            jobs = [(k, frames[k]) for k in range(a, b)]
            if pool is not None:
                got = dict(pool.imap_unordered(_pool_work, jobs, chunksize=1))
            else:
                got = dict(_pool_work(j) for j in jobs)
            res = [got[k] for k in range(a, b)]
            out[a:b] = res
            if path:
                off, flat, st = _pack_results(res)
                tmp = path + '.tmp.npz'
                np.savez(tmp, off=off, flat=flat, st=st)
                os.replace(tmp, path)
            log('  frames %d-%d done: %.1f s per frame in the workers, %d localizations, %.0f s since the start' % (
                a, b - 1, float(np.mean([r['secs'] for r in res])), sum(len(r['det']) for r in res), time.time() - t_start))
    finally:
        if pool is not None:
            pool.close()
            pool.join()
    return out


# ----------------------------------------------------------------------------------------------------------------------------------
# packing localizations (ragged per-frame arrays <-> offsets + flat arrays)
# ----------------------------------------------------------------------------------------------------------------------------------
def pack_rows(per_frame, ncol):
    off = np.concatenate([[0], np.cumsum([len(d) for d in per_frame])]).astype(np.int64)
    flat = np.concatenate([np.asarray(d, dtype=np.float64).reshape(-1, ncol) for d in per_frame if len(d)]) if off[-1] else np.zeros((0, ncol))
    return off, flat


def frames_of(off):
    """Frame index of every row of a packed array."""
    return np.repeat(np.arange(len(off) - 1), np.diff(off)).astype(np.int64)


ISS_COLS = ('x_nm', 'y_nm', 'intensity_adu', 'sigma_px', 'own_significance')
PIC_COLS = ('x_nm', 'y_nm', 'intensity_adu', 'sx_px', 'sy_px', 'net_gradient')


def iss_table(res, pixel_nm=PIXEL_NM):
    """Packed ISS localizations (off, flat (M, 5): x_nm, y_nm, intensity_adu, sigma_px, own_significance), the corner convention nm = (p + 0.5) pixel_nm."""
    per = [np.c_[(r['det'][:, 0] + 0.5) * pixel_nm, (r['det'][:, 1] + 0.5) * pixel_nm, r['det'][:, 3], r['det'][:, 2], r['det'][:, 4]] if len(r['det']) else np.zeros((0, 5))
           for r in res]
    return pack_rows(per, 5)


def pic_table(dets, pixel_nm=PIXEL_NM):
    """Packed Picasso localizations (off, flat (M, 6): x_nm, y_nm, intensity_adu (the harness photons x 6 DN per photon: the fitted integrated counts), sx_px, sy_px,
    net_gradient) from the per-frame arrays of ``smlm_common.run_picasso``."""
    per = [np.c_[d[:, 0], d[:, 1], d[:, C.PH_COL] * C.DN_PER_PHOTON, np.abs(d[:, C.SX_COL]) / pixel_nm, np.abs(d[:, C.SY_COL]) / pixel_nm, d[:, C.NG_COL]] if len(d)
           else np.zeros((0, 6)) for d in dets]
    return pack_rows(per, 6)


def densest_window(xy_list, size_nm=ZOOM_NM, field_nm=12800.0, step_nm=ZOOM_STEP_NM):
    """The size x size nm window (top-left corner on a grid of step_nm inside the field) that contains the most localizations of all the sets in ``xy_list``
    together (ties: the first in row-major order, y then x). Returns (x0, y0, x1, y1) in nm."""
    n = int(round(field_nm / step_nm))
    k = int(round(size_nm / step_nm))
    if k > n:
        raise ValueError('window larger than the field')
    h = np.zeros((n, n))
    for xy in xy_list:
        xy = np.asarray(xy, dtype=np.float64).reshape(-1, 2) if np.size(xy) else np.zeros((0, 2))
        ok = (xy[:, 0] >= 0) & (xy[:, 0] < field_nm) & (xy[:, 1] >= 0) & (xy[:, 1] < field_nm)
        h += np.bincount((xy[ok, 1] // step_nm).astype(np.int64) * n + (xy[ok, 0] // step_nm).astype(np.int64), minlength=n * n).reshape(n, n)
    s = np.zeros((n + 1, n + 1))
    s[1:, 1:] = h.cumsum(0).cumsum(1)
    win = s[k:, k:] - s[:-k, k:] - s[k:, :-k] + s[:-k, :-k]                  # (n - k + 1, n - k + 1): [iy, ix]
    iy, ix = np.unravel_index(int(np.argmax(win)), win.shape)
    return (ix * step_nm, iy * step_nm, ix * step_nm + size_nm, iy * step_nm + size_nm)


# ----------------------------------------------------------------------------------------------------------------------------------
# the precision proxy: repeat-localization precision from consecutive-frame nearest neighbours
# ----------------------------------------------------------------------------------------------------------------------------------
def next_frame_distances(xy_nm, frames, n_frames):
    """For every localization the distance (nm) to the NEAREST localization of the next frame (NaN when the next frame has none or this is the last frame).
    xy_nm (n, 2+), frames (n,) integer frame indices in 0 .. n_frames - 1. Returns (n,)."""
    xy = np.asarray(xy_nm, dtype=np.float64).reshape(len(xy_nm), -1)[:, :2] if np.size(xy_nm) else np.zeros((0, 2))
    fr = np.asarray(frames, dtype=np.int64)
    d = np.full(len(xy), np.nan)
    order = np.argsort(fr, kind='stable')
    lo = np.searchsorted(fr[order], np.arange(n_frames + 1), side='left')
    for t in range(n_frames - 1):
        ia, ib = order[lo[t]:lo[t + 1]], order[lo[t + 1]:lo[t + 2]]
        if len(ia) and len(ib):
            d[ia] = cKDTree(xy[ib]).query(xy[ia], k=1)[0]
    return d


def _mixture_nll(theta, d, wts, dmax):
    w = 1.0 / (1.0 + np.exp(-theta[0]))
    s2 = np.exp(2.0 * theta[1])
    lam = np.exp(theta[2])
    pr = d / s2 * np.exp(-d * d / (2.0 * s2)) / (-np.expm1(-dmax * dmax / (2.0 * s2)))
    pu = 2.0 * lam * d * np.exp(-lam * d * d) / (-np.expm1(-lam * dmax * dmax))
    return -float(np.sum(wts * np.log(np.maximum(w * pr + (1.0 - w) * pu, 1e-300))))


def fit_repeat_mixture(d, dmax=D_MAX_NM, wts=None, starts=None):
    """Maximum-likelihood fit of the nearest-neighbour distances d (nm, all in (0, dmax]) to the mixture  w Rayleigh(s) + (1 - w) NN-of-unrelated-points(lambda),
    both truncated at dmax: p_r(d) = (d / s^2) exp(-d^2 / 2 s^2), p_u(d) = 2 lambda d exp(-lambda d^2) (the distance to the nearest of a uniform random field of
    density lambda / pi). The same molecule seen in two frames with independent per-axis errors sigma has an offset of per-axis sd sqrt(2) sigma, so its
    distance is Rayleigh with s = sqrt(2) sigma. ``wts``: weights of the observations (bootstrap multiplicities). ``starts``: list of (w, s, lambda) starting
    points (default three). Returns dict(w, s_nm, sigma_nm = s / sqrt(2), lam, nll, n)."""
    d = np.asarray(d, dtype=np.float64)
    wts = np.ones(len(d)) if wts is None else np.asarray(wts, dtype=np.float64)
    keep = (d > 0) & (d <= dmax) & (wts > 0)
    d, wts = d[keep], wts[keep]
    if len(d) < 30:
        return dict(w=float('nan'), s_nm=float('nan'), sigma_nm=float('nan'), lam=float('nan'), nll=float('nan'), n=int(wts.sum()))
    starts = starts or [(0.5, 20.0, 2e-5), (0.5, 40.0, 2e-5), (0.3, 70.0, 1e-4)]
    bounds = [(-8.0, 8.0), (np.log(2.0), np.log(dmax)), (np.log(1e-9), np.log(1e-3))]
    best = None
    for w0, s0, l0 in starts:
        th0 = [np.log(w0 / (1 - w0)), np.log(s0), np.log(l0)]
        r = minimize(_mixture_nll, th0, args=(d, wts, dmax), method='L-BFGS-B', bounds=bounds)
        if best is None or r.fun < best.fun:
            best = r
    th = best.x
    s = float(np.exp(th[1]))
    return dict(w=float(1.0 / (1.0 + np.exp(-th[0]))), s_nm=s, sigma_nm=s / np.sqrt(2.0), lam=float(np.exp(th[2])), nll=float(best.fun), n=int(round(wts.sum())))


def repeat_precision(xy_nm, frames, n_frames, dmax=D_MAX_NM, n_boot=0, seed=BOOT_SEED, block=BLOCK):
    """The repeat-localization precision of a localization set (the precision proxy of this experiment): the per-axis sd sigma of the same-molecule component
    of the nearest-neighbour distances between consecutive frames (``fit_repeat_mixture``) and, with ``n_boot`` > 0, a block-bootstrap sample of it (blocks of
    ``block`` consecutive frames, the multiplicities of ``frc.bootstrap_draw``; the source frame of a distance is the earlier one). A precision of the
    molecules that are localized in two consecutive frames (they are brighter than average), not of every localization. Returns dict(fit (the dict of the
    mixture fit), n_pairs (distances <= dmax), median_nm (their median), boot (n_boot,) sigma_nm)."""
    d = next_frame_distances(xy_nm, frames, n_frames)
    ok = np.isfinite(d) & (d <= dmax)
    fit = fit_repeat_mixture(d[ok], dmax)
    out = dict(fit=fit, n_pairs=int(ok.sum()), median_nm=float(np.median(d[ok])) if ok.any() else float('nan'), boot=np.zeros(0))
    if n_boot and np.isfinite(fit['sigma_nm']):
        fr = np.asarray(frames, dtype=np.int64)[ok]
        db = d[ok]
        nb = -(-n_frames // block)
        mult = FRC.bootstrap_draw(nb, seed, n_boot)
        start = [(0.5, fit['s_nm'], fit['lam'])] if np.isfinite(fit['lam']) else None
        boot = np.empty(n_boot)
        for r in range(n_boot):
            boot[r] = fit_repeat_mixture(db, dmax, wts=mult[r][fr // block].astype(np.float64), starts=start)['sigma_nm']
        out['boot'] = boot
    return out


# ----------------------------------------------------------------------------------------------------------------------------------
# agreement of the two sets (a consistency check, NOT an accuracy)
# ----------------------------------------------------------------------------------------------------------------------------------
def agreement(xy_a, fr_a, xy_b, fr_b, radius=AGREE_NM):
    """Fraction of the localizations of set A that have a localization of set B within ``radius`` nm in the SAME frame (and the number of A's localizations)."""
    xa, xb = np.asarray(xy_a, dtype=np.float64).reshape(-1, 2), np.asarray(xy_b, dtype=np.float64).reshape(-1, 2)
    fa, fb = np.asarray(fr_a, dtype=np.int64), np.asarray(fr_b, dtype=np.int64)
    hit = 0
    for t in np.unique(fa):
        ia, ib = np.flatnonzero(fa == t), np.flatnonzero(fb == t)
        if len(ib):
            hit += int((cKDTree(xb[ib]).query(xa[ia], k=1)[0] <= radius).sum())
    return hit / max(len(xa), 1), len(xa)


# ----------------------------------------------------------------------------------------------------------------------------------
# timing: CPU seconds per frame in one process
# ----------------------------------------------------------------------------------------------------------------------------------
def bench_iss(frames, cal, cfg, idx):
    """Wall-clock and CPU seconds of ``process_frame`` on the frames ``idx`` (one process, one warm-up frame first and not counted). The wall clock of the main run
    is inflated by the other jobs of a shared machine; the CPU time (user + system of this process) much less. Returns (len(idx), 4): frame, wall s, cpu s, localizations."""
    ctx = make_context(cal['bg'], cal['a'], cal['b0'], cal['x_ref'], cfg, frames.shape[1:])
    process_frame(frames[idx[0]], ctx)
    rows = []
    for k in idx:
        w0, c0 = time.perf_counter(), time.process_time()
        out = process_frame(frames[k], ctx)
        rows.append([k, time.perf_counter() - w0, time.process_time() - c0, len(out['det'])])
    return np.array(rows, dtype=np.float64)


def bench(raw_path=RAW, n=50, data_root=None, log=print):
    """Add the single-process timing of the ISS on ``n`` evenly spaced frames to an existing raw file (keys iss|bench and meta|bench), re-using its calibration and configuration."""
    with np.load(raw_path, allow_pickle=False) as z:
        raw = {k: z[k] for k in z.files}
    nf = int(raw['meta|shape'][0])
    frames = load_real_frames(data_path(str(raw['meta|set']), data_root), 0, nf)
    law = raw['cal|law']
    cal = dict(bg=raw['cal|bg'], a=float(law[0]), b0=float(law[3]), x_ref=float(law[2]))
    idx = np.unique(np.linspace(0, nf - 1, n).astype(int))
    la0 = os.getloadavg()
    rows = bench_iss(frames, cal, json.loads(str(raw['iss|config'])), idx)
    la1 = os.getloadavg()
    raw['iss|bench'] = rows
    raw['meta|bench'] = np.array(json.dumps(dict(n=len(idx), loadavg_before=[round(x, 2) for x in la0], loadavg_after=[round(x, 2) for x in la1], cpu_count=os.cpu_count())))
    C.save_npz_deterministic(raw_path, raw)
    log('bench: %d frames, %.2f s wall and %.2f s CPU per frame (mean), load average %s -> %s' % (len(idx), rows[:, 1].mean(), rows[:, 2].mean(), la0, la1))
    return raw


def calib_check(raw_path=RAW, seg=100, data_root=None, log=print):
    """Add the stability of the noise-law scale along the sequence (``scale_by_segment``) to an existing raw file (key cal|scale_segments)."""
    with np.load(raw_path, allow_pickle=False) as z:
        raw = {k: z[k] for k in z.files}
    nf = int(raw['meta|shape'][0])
    frames = load_real_frames(data_path(str(raw['meta|set']), data_root), 0, nf)
    raw['cal|scale_segments'] = scale_by_segment(frames, seg)
    C.save_npz_deterministic(raw_path, raw)
    log('calibration check: scale by segment of %d frames: %s' % (seg, ', '.join('%.2f' % r[2] for r in raw['cal|scale_segments'])))
    return raw


# ----------------------------------------------------------------------------------------------------------------------------------
# run
# ----------------------------------------------------------------------------------------------------------------------------------
def run(n_frames=None, workers=2, cache_dir=None, picasso_python=None, out_raw=RAW, set_name=SET, data_root=None, log=print):
    """Localize the frames with Picasso (R0 and the low-threshold variation) and the ISS, calibrate, and write the raw file. Returns the raw dict."""
    path = data_path(set_name, data_root)
    frames = load_real_frames(path, 0, n_frames)
    nf, ny, nx = frames.shape
    log('%s: %d frames of %d x %d px, ADU %.0f .. %.0f' % (set_name, nf, ny, nx, frames.min(), frames.max()))
    py = picasso_python or os.environ.get('SMLM_PICASSO_PYTHON')
    pic, info = {}, {}
    for key, params in (('pic', None), ('pic_ng1000', dict(min_net_gradient=PIC_LOW_NG))):
        info[key] = {}
        t0 = time.time()
        pic[key] = C.run_picasso(frames, params=params, python=py, info=info[key])
        log('Picasso %s: %d localizations, %.4f s per frame (identification + fit, one process), %.0f s wall' % (
            key, sum(len(d) for d in pic[key]), info[key]['run_s'] / nf, time.time() - t0))
    cal = calibrate(frames, pic['pic'])
    law, sig = cal['law'], cal['sig']
    log('calibration: sigma0 = %.3f px (%d bright isolated of %d R0 localizations%s), noise law a = %.3f, variance at the dark level %.1f ADU^2 (x_ref %.1f ADU; the single-frame law times %.3f)' % (
        sig['sigma0'], sig['n_sel'], sig['n_all'], ', FALLBACK' if sig['fallback'] else '', law['a'], law['b0'], law['x_ref'], law['scale']))
    cfg = iss_config(sig['sigma0'])
    t0 = time.time()
    res = run_iss(frames, dict(bg=cal['bg'], a=law['a'], b0=law['b0'], x_ref=law['x_ref']), cfg, workers=workers, cache_dir=cache_dir, log=log)
    log('ISS: %d localizations, %.2f s per frame, %.0f s wall' % (sum(len(r['det']) for r in res), float(np.mean([r['secs'] for r in res])), time.time() - t0))
    edges = [230, 245, 255, 265, 275, 290, 310, 340, 380, 430, 500, 600, 800, 1300]
    tv = temporal_variance_by_level(frames, edges)
    raw = {'README': np.array(RAW_README)}
    raw['meta|set'] = np.array(set_name)
    raw['meta|shape'] = np.array([nf, ny, nx], dtype=np.int64)
    raw['meta|pixel_nm'] = np.array(PIXEL_NM)
    raw['meta|zip_sha256'] = np.array(file_sha256(path))
    raw['meta|versions'] = np.array(json.dumps(dict(python=platform.python_version(), numpy=np.__version__, picasso=str(info['pic'].get('picasso_version', '')),
                                                    cpu_count=os.cpu_count(), loadavg_at_end=[round(x, 2) for x in os.getloadavg()], workers=workers)))
    raw['cal|bg'] = cal['bg']
    raw['cal|mean_image'] = frames.mean(axis=0).astype(np.float32)
    raw['cal|frame_delta'] = np.array([r['delta'] for r in res])
    raw['cal|law'] = np.array([law['a'], law['c'], law['x_ref'], law['b0'], law['n_frames'], law['n_pixels'], law['scale'], law['a_hf'], law['b0_hf'], law['v_dark'], law['x_dark'], law['n_dark']])
    raw['cal|law_x'], raw['cal|law_v'], raw['cal|law_n'] = law['x_bins'], law['v_bins'], law['n_bins']
    raw['cal|temporal_levels'] = np.array(tv).reshape(-1, 3)
    raw['cal|sigma0'] = np.array(json.dumps(sig))
    raw['cal|noise_params'] = np.array(json.dumps(dict(step=NL_STEP, f_lo=NL_F_LO, f_hi=NL_F_HI, sigma_mu=NL_SIGMA_MU, border=NL_BORDER, quantiles=list(NL_QUANTILES))))
    off, flat = iss_table(res)
    raw['iss|config'] = np.array(json.dumps(cfg))
    raw['iss|off'], raw['iss|flat'] = off, flat
    raw['iss|stats'] = _pack_results(res)[2]
    for key in ('pic', 'pic_ng1000'):
        raw['%s|off' % key], raw['%s|flat' % key] = pic_table(pic[key])
        raw['%s|run_s' % key] = np.array([info[key]['run_s'], info[key]['n_frames']])
        raw['%s|config' % key] = np.array(json.dumps(C._merge_params(None if key == 'pic' else dict(min_net_gradient=PIC_LOW_NG)), sort_keys=True, default=C._json_default))
    if ny != nx:
        raise ValueError('the FRC field is square; got %d x %d px frames' % (ny, nx))
    win = densest_window([raw['iss|flat'][:, :2], raw['pic|flat'][:, :2]], ZOOM_NM, max(ny, nx) * PIXEL_NM)
    raw['zoom|window_nm'] = np.array(win)
    for key in ('iss', 'pic'):
        flat, o = raw['%s|flat' % key], raw['%s|off' % key]
        m = (flat[:, 0] >= win[0]) & (flat[:, 0] < win[2]) & (flat[:, 1] >= win[1]) & (flat[:, 1] < win[3])
        raw['zoom|%s' % key] = np.c_[flat[m, 0], flat[m, 1], flat[m, 2], frames_of(o)[m]]
    os.makedirs(os.path.dirname(out_raw), exist_ok=True)
    C.save_npz_deterministic(out_raw, raw)
    log('wrote %s (%.1f MB)' % (out_raw, os.path.getsize(out_raw) / 1e6))
    return raw


RAW_README = ('Raw output of experiments/exp_smlm_real.py (real High Density sequence of the 2013 challenge). meta|*: set, [frames, rows, columns], assumed pixel size (nm), '
              'sha256 of the zip, versions and machine load. cal|*: bg (the background surface, ADU), mean_image (ADU, float32), frame_delta (per-frame offset against bg, ADU), '
              'law = [a, c, x_ref, b0, frames used, pixels used, scale, a_hf, b0_hf, v_dark, x_dark, n_dark]: the noise law var = a x + c (x in ADU) = a (x - x_ref) + b0 actually used = scale x the '
              'single-frame high-frequency law (a_hf, b0_hf), the scale being the temporal variance of the darkest pixels (v_dark ADU^2 at the mean level x_dark, n_dark pixels) over that law there; '
              'law_x / law_v / law_n (the bins the single-frame law was fitted to), temporal_levels = rows [temporal median level, temporal variance, pixels], '
              'sigma0 (JSON: the label-free PSF width and its sample), noise_params (JSON). Localization sets "iss", "pic" (Picasso R0) and "pic_ng1000" (Picasso with minimum net '
              'gradient 1000): <set>|off (frames + 1 int64 offsets: the localizations of frame i are rows off[i]:off[i + 1] of <set>|flat), <set>|flat columns '
              'iss: ' + ', '.join(ISS_COLS) + '; pic*: ' + ', '.join(PIC_COLS) + ' (positions in nm, corner convention (0, 0) = outer corner of the first pixel, x along the columns; '
              'intensity_adu = fitted integrated counts above the background; Picasso\'s = its photons x 6, the harness camera defaults, which are NOT this camera\'s, so only proportional); '
              'iss|stats rows [entries of the path, hit the cap (1/0), seconds, frame offset, weighted residual chi2 of the chosen block] per frame; iss|config, pic*|config (JSON); pic*|run_s = [seconds of identification + fit, frames]. '
              'cal|scale_segments (optional, added by --calib-check) rows [first frame, stop frame, scale, temporal dark variance, dark level] of the noise-law scale on consecutive segments; '
              'iss|bench (optional, added by --bench) rows [frame, wall s, CPU s, localizations] of the single-process timing, meta|bench its JSON (load averages, CPUs); '
              'zoom|window_nm = [x0, y0, x1, y1] of the densest 3 x 3 um window of the two sets together; zoom|iss, zoom|pic = the localizations inside it (x_nm, y_nm, intensity_adu, frame).')


# ----------------------------------------------------------------------------------------------------------------------------------
# analysis of one localization set and the report (from the raw alone)
# ----------------------------------------------------------------------------------------------------------------------------------
def md_table(header, rows):
    out = ['| ' + ' | '.join(header) + ' |', '|' + '|'.join(['---'] * len(header)) + '|']
    return '\n'.join(out + ['| ' + ' | '.join(str(c) for c in r) + ' |' for r in rows])


def nm(x, nd=1):
    return 'no crossing' if x is None or not np.isfinite(x) else '%.*f' % (nd, x)


def rng_str(v, nd=1):
    v = np.asarray(v, dtype=np.float64)
    v = v[np.isfinite(v)]
    return 'n/a' if len(v) == 0 else '%.*f [%.*f, %.*f]' % (nd, np.median(v), nd, v.min(), nd, v.max())


def nanmedian(v):
    """Median of the finite values (NaN if there are none; no warning)."""
    v = np.asarray(v, dtype=np.float64)
    v = v[np.isfinite(v)]
    return float(np.median(v)) if len(v) else float('nan')


def nanpercentile(v, q):
    v = np.asarray(v, dtype=np.float64)
    v = v[np.isfinite(v)]
    return [float(x) for x in np.percentile(v, q)] if len(v) else [float('nan')] * len(np.atleast_1d(q))


def max_abs_diff(a, b):
    """Largest |a - b| over the entries where both are finite (0.0 if there are none)."""
    d = np.abs(np.asarray(a, dtype=np.float64) - np.asarray(b, dtype=np.float64))
    d = d[np.isfinite(d)]
    return float(d.max()) if len(d) else 0.0


def pct_str(v, nd=1):
    """median [2.5 %, 97.5 %] of a bootstrap sample (NaN ignored; the number of NaN is reported separately)."""
    v = np.asarray(v, dtype=np.float64)
    v = v[np.isfinite(v)]
    return 'n/a' if len(v) == 0 else '%.*f [%.*f, %.*f]' % (nd, np.median(v), nd, np.percentile(v, 2.5), nd, np.percentile(v, 97.5))


def get_set(raw, key):
    """(xy_nm (n, 2), frame (n,), flat) of a stored localization set (``key``: iss, pic, pic_ng1000, or iss_interior: the ISS units whose fitted width is not at a bound)."""
    base = 'iss' if key.startswith('iss') else key
    off, flat = raw['%s|off' % base], raw['%s|flat' % base]
    fr = frames_of(off)
    if key == 'iss_interior':
        lo, hi = json.loads(str(raw['iss|config']))['widths']
        m = (flat[:, 3] > lo + WIDTH_TOL) & (flat[:, 3] < hi - WIDTH_TOL)
        flat, fr = flat[m], fr[m]
    ok = np.isfinite(flat[:, :2]).all(axis=1)
    return flat[ok, :2], fr[ok], flat[ok]


def analyse_set(xy, fr, n_frames, field_nm, n_splits=N_SPLITS, n_boot=N_BOOT):
    """Everything the report prints about one localization set (xy (n, 2) nm, fr (n,) frame indices): the FRC resolution by every split, the block bootstrap
    (``n_frames``: the same resamples for every set of the sequence) and the repeat precision."""
    inside = FRC.inside_mask(xy, field_nm)
    res = dict(n=len(xy), n_inside=int(inside.sum()), n_frames=n_frames, empty_frames=int(n_frames - len(np.unique(fr[inside]))))
    alt = FRC.frc_localizations(xy, fr, field_nm, FRC_PIXEL_NM, 'alternating', 0, BLOCK)
    blocks = [FRC.frc_localizations(xy, fr, field_nm, FRC_PIXEL_NM, 'blocks', s, BLOCK) for s in range(n_splits)]
    rand = [FRC.frc_localizations(xy, fr, field_nm, FRC_PIXEL_NM, 'random', s, BLOCK) for s in range(n_splits)]
    odd = FRC.frc_localizations(xy, fr, field_nm, FRC_PIXEL_NM, 'oddeven', 0, BLOCK)
    res['alt'] = alt['resolution_nm']
    res['alt_crossed'] = alt['crossed']
    res['alt_last'] = alt['last_nm']
    res['alt_0143'] = FRC.frc_resolution(alt, FRC_PIXEL_NM, FRC.SMOOTH, 0.143)['resolution_nm']
    res['blocks'] = np.array([r['resolution_nm'] for r in blocks])
    res['blocks_crossed'] = np.array([r['crossed'] for r in blocks])
    res['random'] = np.array([r['resolution_nm'] for r in rand])
    res['oddeven'] = odd['resolution_nm']
    res['oddeven_crossed'] = odd['crossed']
    res['max_0143_diff'] = max_abs_diff([res['alt_0143']] + [FRC.frc_resolution(r, FRC_PIXEL_NM, FRC.SMOOTH, 0.143)['resolution_nm'] for r in blocks],
                                        [alt['resolution_nm']] + [r['resolution_nm'] for r in blocks])
    res['curve'] = alt
    res['halves'] = (alt['n_a'], alt['n_b'])
    boot = FRC.bootstrap_blocks(xy, fr, n_boot=n_boot, seed=BOOT_SEED, block=BLOCK, field_nm=field_nm, pixel_nm=FRC_PIXEL_NM, n_frames=n_frames) if n_boot else None
    res['boot'] = boot['resolution_nm'] if boot else np.zeros(0)
    res['boot_crossed'] = boot['crossed'] if boot else np.zeros(0, bool)
    res['repeat'] = repeat_precision(xy, fr, n_frames, D_MAX_NM, n_boot=n_boot, seed=BOOT_SEED, block=BLOCK)
    return res


def ratio_str(num, den, nd=2):
    """median [2.5 %, 97.5 %] of the elementwise ratio of two paired bootstrap samples (NaN ignored)."""
    num, den = np.asarray(num, dtype=np.float64), np.asarray(den, dtype=np.float64)
    if len(num) == 0 or len(num) != len(den):
        return 'n/a'
    r = num / den
    r = r[np.isfinite(r)]
    if len(r) == 0:
        return 'n/a'
    return '%.*f [%.*f, %.*f]' % (nd, np.median(r), nd, np.percentile(r, 2.5), nd, np.percentile(r, 97.5))


def rel_phrase(x, ref, tol=0.10):
    """Words for the relation of a resolution x (nm, smaller is finer) to a reference: the same within ``tol``, or finer / coarser by a percentage."""
    if not (np.isfinite(x) and np.isfinite(ref) and ref > 0):
        return 'not comparable'
    r = x / ref
    if abs(r - 1) <= tol:
        return 'the same within %.0f %% (%.2fx)' % (100 * tol, r)
    return ('finer by %.0f %% (%.2fx)' if r < 1 else 'coarser by %.0f %% (%.2fx)') % (100 * abs(1 - r), r)


def thinned_resolution(xy, fr, n_keep, field_nm, seeds=THIN_SEEDS):
    """Alternating-block FRC resolution (nm) of random subsets of ``n_keep`` localizations of a set (without replacement; the localizations inside the field only; one subset
    per seed, ``numpy.random.default_rng(seed)``): the resolution at a given number of localizations, to separate the effect of how many a method reports from that of how
    good they are. NaN entries mean no crossing."""
    xy = np.asarray(xy, dtype=np.float64)
    fr = np.asarray(fr)
    ok = FRC.inside_mask(xy, field_nm)
    xy, fr = xy[ok], fr[ok]
    if n_keep > len(xy):
        raise ValueError('cannot keep %d of %d localizations' % (n_keep, len(xy)))
    out = []
    for sd in seeds:
        idx = np.sort(np.random.default_rng(sd).choice(len(xy), size=n_keep, replace=False))
        out.append(FRC.frc_localizations(xy[idx], fr[idx], field_nm, FRC_PIXEL_NM, 'alternating', 0, BLOCK)['resolution_nm'])
    return np.array(out)


def build_md(raw, n_boot=N_BOOT, n_splits=N_SPLITS, raw_sha=None, log=lambda s: None):
    """The report (a string) from the raw arrays alone: the same raw and arguments always give the same text."""
    nf, ny, nx = (int(v) for v in raw['meta|shape'])
    pixel = float(raw['meta|pixel_nm'])
    field = max(ny, nx) * pixel
    law = raw['cal|law']
    sig = json.loads(str(raw['cal|sigma0']))
    cfg = json.loads(str(raw['iss|config']))
    versions = json.loads(str(raw['meta|versions']))
    stats = raw['iss|stats']
    an = {}
    for key, label, kind, main in ROWS:
        xy, fr, flat = get_set(raw, key)
        log('analysing %s (%d localizations)' % (key, len(xy)))
        an[key] = dict(res=analyse_set(xy, fr, nf, field, n_splits, n_boot), xy=xy, fr=fr, flat=flat, label=label, kind=kind, main=main)
    n_ref = an['pic']['res']['n_inside']
    thin = {k: thinned_resolution(an[k]['xy'], an[k]['fr'], n_ref, field) for k in ('iss', 'pic_ng1000')}
    L = []
    P = L.append
    P('# Real experimental SMLM data: FRC resolution of the sliding inverse scale space against Picasso\n')
    P('Generated by `experiments/exp_smlm_real.py --from-raw` from `results/smlm/real_raw.npz`%s; every number below is computed from that file by `experiments/frc.py` and the functions of the script '
      '(the file regenerates byte for byte). The protocol and the a priori configuration are in the docstring of the script, the fetch script and the facts read on the pages in '
      '`data/smlm/fetch_real.py`.\n' % (' (sha256 `%s`)' % raw_sha if raw_sha else ''))
    P('**What this measures and what it does not.** The data are the real High Density sequence of the 2013 localization-microscopy challenge (Tubulin Alexa647 in a fixed cell, '
      '%d frames of %d x %d px, courtesy of N. Olivier, D. Keller and S. Manley, https://srm.epfl.ch/srm/dataset/challenge-2D-real/Real_High_Density/index.html). **There is no ground truth '
      '("Ground-truth localization unknown"). The experiment measures resolution (Fourier ring correlation) and a precision proxy (repeat-localization precision) only; it does not '
      'measure detection accuracy**: how many localizations of either method lie on a real filament, or on the same filament as the other method\'s, is not knowable here, and a method '
      'that reports more (or split) localizations is not thereby more accurate. FRC rewards localizations that two independent halves of the data agree on, and it can be flattered by '
      'consistent artefacts.\n' % (nf, ny, nx))

    # ------------------------------------------------------------------ summary
    ra, rb = an['pic']['res'], an['iss']['res']
    rl, rpl = an['pic_ng1000']['res'], an['iss_interior']['res']

    def rec(r):
        bt = r['boot'][np.isfinite(r['boot'])] if len(r['boot']) else np.zeros(0)
        if not len(bt) or not np.isfinite(r['alt']):
            return (float('nan'), float('nan'))
        med = np.median(bt)
        lo, hi = np.percentile(bt, [2.5, 97.5])
        return (r['alt'] * lo / med, r['alt'] * hi / med)
    ia, ib = rec(ra), rec(rb)
    cpu = ('%.2f s CPU per frame' % raw['iss|bench'][:, 2].mean()) if 'iss|bench' in raw else '%.2f s wall per frame' % stats[:, 2].mean()
    P('## Summary\n')
    P('* Localizations (%d frames of %d x %d px): Picasso R0 **%d** (%.1f per frame), sliding ISS **%d** (%.1f per frame, %.1fx).' % (
        nf, ny, nx, ra['n'], ra['n'] / nf, rb['n'], rb['n'] / nf, rb['n'] / max(ra['n'], 1)))
    P('* FRC resolution of the half sets (alternating blocks of %d frames, %g nm pixels, threshold 1/7; recentred bootstrap interval): Picasso R0 **%s nm** [%s, %s], ISS **%s nm** [%s, %s]; '
      'ratio ISS / Picasso R0 %s on the full sample, %s in the paired block bootstrap.' % (
          BLOCK, FRC_PIXEL_NM, nm(ra['alt']), nm(ia[0]), nm(ia[1]), nm(rb['alt']), nm(ib[0]), nm(ib[1]),
          '%.2f' % (rb['alt'] / ra['alt']) if np.isfinite(ra['alt']) and np.isfinite(rb['alt']) else 'n/a', ratio_str(rb['boot'], ra['boot'])))
    P('* **Count-controlled comparison**: the ISS thinned at random to Picasso R0\'s count (%d localizations; 5 thinnings) has an FRC resolution of %s nm against R0\'s %s nm, i.e. %s; Picasso with its minimum net gradient '
      'lowered to 1000 reaches %s nm with %d localizations (%s nm when thinned to R0\'s count). The resolution of each method improves with the number of localizations it reports (it cannot tell structure from '
      'consistent artefacts), so the difference between the two headline rows is %s.' % (
          n_ref, rng_str(thin['iss']), nm(ra['alt']), rel_phrase(nanmedian(thin['iss']), ra['alt']), nm(rl['alt']), rl['n_inside'], rng_str(thin['pic_ng1000']),
          'mostly the difference in the number of localizations' if abs(nanmedian(thin['iss']) / ra['alt'] - 1) <= 0.10 else 'not only the difference in the number of localizations'))
    P('* Repeat-localization precision (per axis): Picasso R0 %s nm, ISS %s nm. Time per frame: Picasso R0 %.4f s, ISS %s (a priori configuration, one core, shared machine).' % (
        nm(ra['repeat']['fit']['sigma_nm']), nm(rb['repeat']['fit']['sigma_nm']), raw['pic|run_s'][0] / raw['pic|run_s'][1], cpu))
    P('* Both methods localize the same filaments (%.1f %% of Picasso R0 localizations have an ISS partner within %g nm in the same frame; the reconstructions correlate with r = %.3f), but the ISS fit is '
      'far from the data model (chi2 %.1f, %.0f %% of its units at the lower width bound, PSF elongated by ~%.0f %% against its circular kernel). **No ground truth: resolution and precision proxies only, not detection accuracy.**\n' % (
          100 * agreement(an['pic']['xy'], an['pic']['fr'], an['iss']['xy'], an['iss']['fr'])[0], AGREE_NM,
          float(np.corrcoef(FRC.render(an['pic']['xy'], field, 100.0).ravel(), FRC.render(an['iss']['xy'], field, 100.0).ravel())[0, 1]),
          stats[:, 4].mean(), 100 * np.mean(an['iss']['flat'][:, 3] <= json.loads(str(raw['iss|config']))['widths'][0] + WIDTH_TOL),
          100 * (sig['median_sy'] / sig['median_sx'] - 1)))

    # ------------------------------------------------------------------ assumptions
    P('## What was read, what was assumed\n')
    P('* **Read on the hub\'s pages (2026-11-07):** the sample (fixed cell, mouse anti-alpha-tubulin primary, Alexa647 secondary), "High density settings, 500 frames", 16-bit TIFF, the '
      'courtesy, "Ground-truth localization unknown". **Read from the files:** %d frames of %d x %d px (one 16-bit TIFF per frame), sha256 of the zip `%s`.' % (nf, ny, nx, str(raw['meta|zip_sha256'])))
    P('* **Not given anywhere, therefore ASSUMED:** the **pixel size, taken as %g nm** (all nm below scale with it; a different pixel size multiplies every resolution and precision in nm by the '
      'same factor and leaves the comparison of the methods unchanged); consistent with the field of 12.8 um x 12.8 um, with a fitted PSF width of ~150 nm and with the density arithmetic of arXiv:1801.09631 (if their 500-frame experimental dataset is this sequence, which their text does not state), not proven by any of them. '
      'The camera is not given either: Picasso received the harness defaults (baseline 100, EM gain 300, 45 e-/ADU, QE 0.9, those of the simulated sets, not of this camera; they affect only '
      'its photon column, not its identification or positions) and the ISS works in raw ADU with the noise law below.' % pixel)
    P('* The frames are **not circular**: the bright isolated Picasso fits have median sx = %.2f px and sy = %.2f px (the PSF is elongated along y by ~%.0f %%); the ISS kernel is circular and the '
      'packaged model has no ellipticity parameter. This is a model mismatch of the a priori configuration, not something that was tuned away.\n' % (
          sig['median_sx'], sig['median_sy'], 100 * (sig['median_sy'] / sig['median_sx'] - 1)))

    # ------------------------------------------------------------------ calibration
    P('## Label-free calibration of the ISS (no ground truth, nothing tuned on the FRC)\n')
    P('* PSF width at entry: **sigma0 = %.3f px** (= %.0f nm at the assumed pixel), the median of sqrt(sx sy) over %d bright (integrated counts >= %.0f ADU, the median of all R0 localizations) isolated (no other localization of the frame within %g px) '
      'Picasso-R0 localizations out of %d (%d isolated, %d bright); for comparison the median over all R0 localizations is %.3f px and the simulated sets used 1.4 px. Entry width, sliding-width bounds '
      '(%.2f, %.2f) px and resolution limit r_ex = %.2f px follow from it (the simulated sets\' (1.0, 3.2) px and 2.0 px at 1.4 px).' % (
          sig['sigma0'], sig['sigma0'] * pixel, sig['n_sel'], sig['photon_threshold'] * C.DN_PER_PHOTON, ISO_PX, sig['n_all'], sig['n_iso'], sig['n_bright'],
          sig['median_all'], cfg['widths'][0], cfg['widths'][1], cfg['r_ex']))
    P('* Noise law, shape from the high-frequency content of every %dth frame (%d frames, %d pixels): var = a_hf x + c_hf with x the raw ADU level, **a_hf = %.3f, c_hf = %.1f ADU^2**, i.e. with the reference '
      'dark level x_ref = %.1f ADU (the mean of the darkest 1 %% of the smoothed pixel values; a reference of the regression, NOT the camera baseline) var = a_hf (x - x_ref) + b0_hf with b0_hf = %.1f ADU^2. '
      'The bins the line was fitted to (local mean ADU: robust variance ADU^2, pixels): %s.' % (
          NL_STEP, int(law[4]), int(law[5]), law[7], law[8] - law[7] * law[2], law[2], law[8],
          '; '.join('%.0f: %.0f (%d)' % (x, v, n) for x, v, n in zip(raw['cal|law_x'], raw['cal|law_v'], raw['cal|law_n']))))
    P('* Noise law, absolute scale from the temporal variance of the darkest pixels (lowest %g %% of the temporal median, %d pixels, level %.0f ADU, each frame\'s median removed, robust): %.0f ADU^2, against %.0f ADU^2 '
      'from the single-frame law at that level, i.e. the single-frame law (which sees the white part of the noise only) is multiplied by **s = %.3f**. **Law used: var = a (x - x_ref) + b0 with a = %.3f and b0 = %.1f ADU^2** '
      '(in raw ADU var = a x %+.1f). For purely white noise s would be 1. The temporal variance of the pixels by level of their temporal median, level: variance ADU^2 (pixels): %s - it explodes above the dark level '
      'because the pixels on filaments blink, which is why the temporal estimate of the simulated sets cannot give the dependence on the level here and only the darkest pixels anchor the scale.' % (
          100 * TD_QUANTILE, int(law[11]), law[10], law[9], law[7] * (law[10] - law[2]) + law[8], law[6], law[0], law[3], law[1],
          '; '.join('%.0f: %.0f (%d)' % (a, b, int(c)) for a, b, c in raw['cal|temporal_levels'])))
    if 'cal|scale_segments' in raw:
        ss = raw['cal|scale_segments']
        P('* Stability of that scale along the sequence (`--calib-check`: the same rule applied to consecutive segments of %d frames; first frame - scale s, temporal dark variance ADU^2): %s; the scale used (all %d frames) is %.3f. '
          'The temporal anchor is sensitive to the drift of the background in the first frames: a calibration on a short early stretch would have given a larger scale (a noisier law, fewer ISS units).' % (
              int(ss[0, 1] - ss[0, 0]), '; '.join('%d - %.2f, %.0f' % (r[0], r[2], r[3]) for r in ss), nf, law[6]))
    delta = raw['cal|frame_delta']
    bg = raw['cal|bg']
    P('* Background: the robust quadratic surface through the per-pixel temporal median (%.0f to %.0f ADU over the field) plus a per-frame offset (median of the frame minus the surface: %.0f ADU in frame 0, '
      '%.0f in the last frame, range %.0f to %.0f).\n' % (bg.min(), bg.max(), delta[0], delta[-1], delta.min(), delta.max()))
    P('ISS configuration actually run: ' + '`' + json.dumps(cfg, sort_keys=True) + '`. In words: configuration "%s" of `exp_smlm_iss.CONFIGS` (high-pass %g px, Gaussian core and a halo of 3x the width), start lattice %g px, '
      'stop rule significance kappa = %g, no width post-filter, no offset, crowded-field local refits (%g px radius, %g px window), no final global polish, cap of %d entries per frame. '
      'Picasso R0: `%s`.\n' % (cfg['cfg'], cfg['hp'], cfg['lat_h'], cfg['kappa'], cfg['fit_radius'], cfg['window'], cfg['cap'], str(raw['pic|config'])))

    # ------------------------------------------------------------------ counts and time
    P('## Localizations and time\n')
    rows = []
    for key, label, kind, main in ROWS:
        r = an[key]['res']
        fl = an[key]['flat']
        if kind == 'iss':
            sg = fl[:, 3]
            extra = 'fitted width (median) %.2f px' % np.median(sg)
        else:
            sg = np.sqrt(fl[:, 3] * fl[:, 4])
            extra = 'sqrt(sx sy) (median) %.2f px' % np.median(sg)
        rows.append([label, r['n'], r['n_inside'], '%.1f' % (r['n'] / nf), r['empty_frames'], '%.0f' % np.median(fl[:, 2]), extra])
    P(md_table(['row', 'localizations', 'inside the field', 'per frame', 'frames without any', 'median intensity (ADU, integrated)', 'width'], rows))
    secs = stats[:, 2]
    sg_all = an['iss']['flat'][:, 3]
    lo_w, hi_w = cfg['widths']

    def near_border(xy, d=400.0):
        xy = xy[FRC.inside_mask(xy, field)]
        return float(np.mean((xy.min(axis=1) < d) | (xy.max(axis=1) > field - d))) if len(xy) else float('nan')
    P('')
    P('Border: %.1f %% of the Picasso-R0 localizations and %.1f %% of the ISS localizations lie within 400 nm of the field border (Picasso cannot identify within box // 2 = 3 px at the top and left and 4 px at the bottom and right '
      'of an unpadded frame, the ISS covers the field; %d ISS and %d Picasso-R0 localizations have a fitted centre outside the field and are not used by the FRC). The FRC window tapers the outer 1.6 um of every side.\n' % (
          100 * near_border(an['pic']['xy']), 100 * near_border(an['iss']['xy']), an['iss']['res']['n'] - an['iss']['res']['n_inside'], an['pic']['res']['n'] - an['pic']['res']['n_inside']))
    P('Seconds per frame (one core): **Picasso R0 %.4f** (identification + fit of the %d frames in one process, JIT compilation excluded; the low-threshold variation %.4f); **sliding ISS %.2f mean, '
      '%.2f median** (min %.2f, max %.2f; wall clock inside a worker process on a shared machine: %d CPUs, other jobs were running, load average %s at the end of the run; the load rose during the run, '
      'see the benchmark below for CPU time). The ISS path ended at the cap of %d entries with every entrant still significant in %d of %d frames; the path had on average %.0f entries per frame and the chosen block %.0f '
      'localizations.\n' % (raw['pic|run_s'][0] / raw['pic|run_s'][1], int(raw['pic|run_s'][1]), raw['pic_ng1000|run_s'][0] / raw['pic_ng1000|run_s'][1], secs.mean(), np.median(secs), secs.min(),
                            secs.max(), versions['cpu_count'], versions['loadavg_at_end'], cfg['cap'], int(stats[:, 1].sum()), nf, stats[:, 0].mean(), an['iss']['res']['n'] / nf))
    if 'iss|bench' in raw:
        bn = raw['iss|bench']
        bi = json.loads(str(raw['meta|bench']))
        P('Timing benchmark (`--bench`, after the main run): the ISS on %d evenly spaced frames in ONE process, **%.2f s CPU per frame (mean; median %.2f)** and %.2f s wall clock (mean), %.0f localizations per frame; '
          'machine load average %s before and %s after (%d CPUs): the CPU time is much less affected by the other jobs than the wall clock of the main run but not immune to them (shared caches, hyper-threads).\n' % (
              len(bn), bn[:, 2].mean(), np.median(bn[:, 2]), bn[:, 1].mean(), bn[:, 3].mean(), bi['loadavg_before'], bi['loadavg_after'], bi['cpu_count']))
    P('Goodness of the ISS fit, label free: the weighted mean squared residual (chi2) of the high-passed frame after the chosen block is %.2f on average (median %.2f, 5 %%-95 %% %.2f to %.2f over the frames); '
      'it is about 1 for a complete model with the right noise law, so the excess is structure the model does not capture (the circular kernel on an elongated PSF, unresolved overlapping '
      'emitters, out-of-focus light) and/or noise above the law. The fitted widths of the ISS units: median %.2f px, **%.0f %% at the lower bound (%.2f px) and %.1f %% at the upper bound (%.2f px) '
      'of the sliding widths** (the entry width is %.2f px).\n' % (stats[:, 4].mean(), np.median(stats[:, 4]), np.percentile(stats[:, 4], 5), np.percentile(stats[:, 4], 95), np.median(sg_all),
                                                              100 * np.mean(sg_all <= lo_w + WIDTH_TOL), lo_w, 100 * np.mean(sg_all >= hi_w - WIDTH_TOL), hi_w, cfg['sigma0']))

    # ------------------------------------------------------------------ FRC
    nblk = -(-nf // BLOCK)
    P('## Fourier ring correlation resolution\n')
    P('FRC of two half-sets rendered as histograms on a %g nm grid over the %.1f um field (%d x %d pixels), Tukey window (0.25), window-weighted mean subtracted, 5-ring moving average, '
      'resolution = 1/f at the first crossing below **1/7** (`experiments/frc.py`; the "0.143" variant differs by at most %.3f nm from it over all the block splits below). '
      'The halves hold half of the localizations each, so these are the resolutions of HALF sets (the full set is somewhat better). **Splits.** The headline split is by BLOCKS: blocks of %d consecutive frames '
      'dealt to the two halves ("alternating": odd block numbers to one half, no randomness; "random blocks": %d random assignments, seeds 0-%d). A fluorophore that stays on for several frames is localized '
      'several times with independent errors, and a random localization split or an odd / even frame split puts those repeats into BOTH halves, which correlates the halves and makes the resolution '
      'look finer than it is (both are reported: the random localization split reads %.0f to %.0f %% finer than the alternating-block split in the rows of the table). **Bootstrap** (as specified): the %d blocks are resampled with replacement and the unique drawn blocks are dealt at random '
      'into the two halves, so that no localization is in both (%d replicates, seed %d). **Its replicates are coarser than the estimate on the full sample and its percentile interval does not bracket it**: '
      'a replicate holds only ~63 %% of the distinct blocks, hence fewer localizations per half, and the resolution depends on that number (the effect is the one shown by the rows of different numbers of localizations); '
      'the validity of the interval as an interval for the full-sample estimate is therefore NOT claimed. What is used from the bootstrap is (i) the RELATIVE spread of the replicates, applied to the full-sample '
      'estimate ("recentred": estimate x [2.5 %% / median, 97.5 %% / median]; expected to be conservative because a smaller sample has a larger relative spread: in the synthetic check of `tests/test_frc.py` the bootstrap spread exceeds the spread across independent datasets) and (ii) the PAIRED ratio of two methods resampled together '
      '(identical resamples and halves), in which the common loss of data cancels to first order.\n' % (
          FRC_PIXEL_NM, field / 1000.0, int(round(field / FRC_PIXEL_NM)), int(round(field / FRC_PIXEL_NM)), max(an[k]['res']['max_0143_diff'] for k, *_ in ROWS), BLOCK, n_splits, n_splits - 1,
          min(100 * (1 - nanmedian(an[k]['res']['random']) / an[k]['res']['alt']) for k, *_ in ROWS), max(100 * (1 - nanmedian(an[k]['res']['random']) / an[k]['res']['alt']) for k, *_ in ROWS),
          nblk, n_boot, BOOT_SEED))
    P('Resolution by split (nm; "no crossing" = better than the pixel limit):\n')
    rows = []
    for key, label, kind, main in ROWS:
        r = an[key]['res']
        rows.append([label, r['n_inside'], '%d + %d' % r['halves'], nm(r['alt']), rng_str(r['blocks']), '%d of %d' % (int(np.isfinite(r['blocks']).sum()), len(r['blocks'])),
                     rng_str(r['random']), nm(r['oddeven'])])
    P(md_table(['row', 'localizations inside the field', 'per half (alternating blocks)', '**alternating blocks (headline)**', 'random blocks, seeds 0-%d: median [min, max]' % (n_splits - 1),
                'random block splits with a crossing', 'random localization split, seeds 0-%d: median [min, max]' % (n_splits - 1), 'odd / even frames'], rows))
    P('')
    P('Block bootstrap (%d replicates; resolution in nm):\n' % n_boot)
    rows = []
    for key, label, kind, main in ROWS:
        r = an[key]['res']
        bt = r['boot'][np.isfinite(r['boot'])] if len(r['boot']) else np.zeros(0)
        if len(bt):
            med = np.median(bt)
            lo, hi = np.percentile(bt, [2.5, 97.5])
            rel = '[%.2f, %.2f]' % (lo / med, hi / med)
            rec = '[%s, %s]' % (nm(r['alt'] * lo / med), nm(r['alt'] * hi / med)) if np.isfinite(r['alt']) else 'n/a'
            rows.append([label, nm(r['alt']), '%.1f [%.1f, %.1f]' % (med, lo, hi), '%.2f' % (med / r['alt']) if np.isfinite(r['alt']) else 'n/a', rel, rec, '%d' % int((~r['boot_crossed']).sum())])
        else:
            rows.append([label, nm(r['alt']), 'n/a', 'n/a', 'n/a', 'n/a', 'n/a'])
    P(md_table(['row', 'alternating-blocks estimate (full sample)', 'bootstrap replicates: median [2.5 %, 97.5 %]', 'replicate median / estimate', 'relative spread [2.5 % / median, 97.5 % / median]',
                '**recentred interval**: estimate x relative spread', 'replicates without a crossing'], rows))
    P('')
    P('Paired comparison, ratio of resolutions (first row of the pair over the second; **below 1: the first has the finer resolution**), from the same %d resamples and halves (seed %d), and on the full sample (alternating blocks):\n' % (n_boot, BOOT_SEED))
    pairs = [('iss', 'pic'), ('iss', 'pic_ng1000'), ('iss_interior', 'pic_ng1000'), ('pic_ng1000', 'pic')]
    labels = dict((k, l) for k, l, _, _ in ROWS)
    rows = []
    for ka, kb in pairs:
        ra_, rb_ = an[ka]['res'], an[kb]['res']
        full = ra_['alt'] / rb_['alt'] if np.isfinite(ra_['alt']) and np.isfinite(rb_['alt']) else float('nan')
        rows.append(['%s / %s' % (labels[ka], labels[kb]), ra_['n_inside'], rb_['n_inside'], '%.2f' % full, ratio_str(ra_['boot'], rb_['boot'])])
    P(md_table(['pair', 'localizations (first)', 'localizations (second)', 'ratio on the full sample', 'paired bootstrap ratio: median [2.5 %, 97.5 %]'], rows))
    P('')
    P('Count-controlled comparison (alternating blocks; random subsets WITHOUT replacement of the localizations inside the field, seeds %s; nm):\n' % ', '.join(str(x) for x in THIN_SEEDS))
    rows = [[an['pic']['label'], an['pic']['res']['n_inside'], 'all', nm(an['pic']['res']['alt'])],
            [an['iss']['label'] + ', thinned', n_ref, 'random subset', rng_str(thin['iss'])],
            [an['pic_ng1000']['label'] + ', thinned', n_ref, 'random subset', rng_str(thin['pic_ng1000'])],
            [an['pic_ng1000']['label'], an['pic_ng1000']['res']['n_inside'], 'all', nm(an['pic_ng1000']['res']['alt'])],
            [an['iss']['label'], an['iss']['res']['n_inside'], 'all', nm(an['iss']['res']['alt'])]]
    P(md_table(['row', 'localizations used', 'which', 'resolution: median [min, max] over the thinnings'], rows))
    P('')
    ra_, rb_ = an['pic']['res'], an['iss']['res']
    P('The first and the last crossing of the smoothed curve of the alternating split are %s nm and %s nm for %s, %s nm and %s nm for %s%s.\n' % (
        nm(ra_['alt']), nm(ra_['alt_last']), an['pic']['label'], nm(rb_['alt']), nm(rb_['alt_last']), an['iss']['label'],
        ' (one crossing each)' if (abs(ra_['alt_last'] - ra_['alt']) < 0.5 and abs(rb_['alt_last'] - rb_['alt']) < 0.5)
        else ': the curve wanders around the threshold in at least one of them, so the crossing is poorly defined to a few nm'))
    P('FRC at fixed spatial frequencies (alternating blocks, smoothed curve; the threshold is 1/7 = 0.1429):\n')
    rows = []
    for key, label, kind, main in ROWS:
        c = an[key]['res']['curve']
        rows.append([label] + ['%.3f' % np.interp(f / 1000.0, c['freq'], c['smoothed']) for f in FREQS_UM])
    P(md_table(['row'] + ['%g /um (%.0f nm)' % (f, 1000.0 / f) for f in FREQS_UM], rows))
    P('')

    # ------------------------------------------------------------------ precision proxy
    P('## Repeat-localization precision (the second, ground-truth-free indicator)\n')
    P('For every localization, the distance to the nearest localization of the NEXT frame; those within %g nm are fitted by maximum likelihood to the mixture w Rayleigh(s) + (1 - w) "nearest of an unrelated '
      'random field" (`repeat_precision`; both truncated at %g nm). The same molecule seen in two frames with independent per-axis errors sigma has a Rayleigh distance with s = sqrt(2) sigma, so '
      '**sigma_rep = s / sqrt(2)** is the per-axis localization precision as seen by the repeats; w is the share of the distances in the same-molecule component. It describes the molecules that are '
      'localized in two consecutive frames (brighter than average), includes any position change from one frame to the next, and is biased by whatever else is nearest. Interval: the block bootstrap '
      '(blocks of %d frames, %d replicates).\n' % (D_MAX_NM, D_MAX_NM, BLOCK, n_boot))
    rows = []
    for key, label, kind, main in ROWS:
        rp = an[key]['res']['repeat']
        f = rp['fit']
        rows.append([label, rp['n_pairs'], '%.1f' % rp['median_nm'], nm(f['sigma_nm']), pct_str(rp['boot']) if len(rp['boot']) else 'n/a', '%.2f' % f['w'] if np.isfinite(f['w']) else 'n/a'])
    P(md_table(['row', 'consecutive-frame pairs within %g nm' % D_MAX_NM, 'median distance nm', 'sigma_rep nm (per axis)', 'bootstrap median [2.5 %, 97.5 %] nm', 'w (same-molecule share)'], rows))
    P('')

    # ------------------------------------------------------------------ agreement
    P('## Do the two methods localize the same structure? (consistency, not accuracy)\n')
    xa, fa = an['pic']['xy'], an['pic']['fr']
    xb, fb = an['iss']['xy'], an['iss']['fr']
    pa, na = agreement(xa, fa, xb, fb)
    pb, nb_ = agreement(xb, fb, xa, fa)
    img_a = FRC.render(xa, field, 100.0)
    img_b = FRC.render(xb, field, 100.0)
    corr = float(np.corrcoef(img_a.ravel(), img_b.ravel())[0, 1])
    P('Within the same frame, **%.1f %% of the %d Picasso-R0 localizations have an ISS localization within %g nm**, and **%.1f %% of the %d ISS localizations have a Picasso-R0 one**; the two pooled '
      'reconstructions rendered on 100 nm pixels correlate with Pearson r = %.3f. This says the two methods put localizations on the same bright structures in the same frames; it does not say that '
      'any of them is right, and the ISS localizations without a Picasso partner (%.1f %% of ISS) are faint spots, extra units in overlapping blobs and extra units of bright elongated blobs, '
      'which this experiment cannot tell apart.\n' % (100 * pa, na, AGREE_NM, 100 * pb, nb_, corr, 100 * (1 - pb)))

    # ------------------------------------------------------------------ zoom export
    w = raw['zoom|window_nm']
    P('## Figure-ready export\n')
    P('`results/smlm/real_raw.npz` holds both full localization sets (`iss|off`, `iss|flat`: x_nm, y_nm, intensity_adu, sigma_px, own_significance; `pic|off`, `pic|flat`: x_nm, y_nm, intensity_adu, sx_px, '
      'sy_px, net_gradient; offsets by frame), the mean image (`cal|mean_image`) and the localizations inside the densest 3 x 3 um window of the two sets together, window [%.0f, %.0f, %.0f, %.0f] nm '
      '(`zoom|window_nm`, `zoom|iss`: %d localizations, `zoom|pic`: %d). No figure is made here.\n' % (w[0], w[1], w[2], w[3], len(raw['zoom|iss']), len(raw['zoom|pic'])))

    # ------------------------------------------------------------------ verdict
    P('## Reading of the numbers\n')
    splits_ordered = [('alternating blocks', rb['alt'] < ra['alt']), ('random blocks (median)', nanmedian(rb['blocks']) < nanmedian(ra['blocks'])),
                      ('random localization split (median)', nanmedian(rb['random']) < nanmedian(ra['random'])), ('odd / even frames', rb['oddeven'] < ra['oddeven'])]
    n_fin = sum(1 for _, v in splits_ordered if v)
    P('* **Headline (alternating blocks, half sets):** Picasso R0 %s nm (recentred interval [%s, %s]), sliding ISS %s nm ([%s, %s]); paired bootstrap ratio ISS / Picasso R0 %s. The ISS has the finer resolution '
      'in %d of the 4 splits tabulated (%s).' % (
          nm(ra['alt']), nm(ia[0]), nm(ia[1]), nm(rb['alt']), nm(ib[0]), nm(ib[1]), ratio_str(rb['boot'], ra['boot']), n_fin,
          ', '.join('%s: %s' % (n, 'ISS finer' if v else 'Picasso finer') for n, v in splits_ordered)))
    P('* **Counts.** The ISS reports %.1fx as many localizations as Picasso R0 (%d against %d inside the field); lowering Picasso\'s minimum net gradient to 1000 raises its count '
      'to %d and moves its resolution from %s to %s nm (paired ratio %s). At comparable counts, the ISS without its units at a width bound (%d localizations, %s nm) is %s than Picasso at net gradient 1000 (%d, %s nm; '
      'paired ratio %s). The ISS thinned at random to the count of Picasso R0 (%d) has %s nm against Picasso R0\'s %s nm (%s), and Picasso at net gradient 1000 thinned to the same count %s nm. '
      'So the comparison "ISS against Picasso R0" is %s, and an FRC does not tell structure from consistent artefacts.' % (
          rb['n_inside'] / max(ra['n_inside'], 1), rb['n_inside'], ra['n_inside'], rl['n_inside'], nm(ra['alt']), nm(rl['alt']), ratio_str(rl['boot'], ra['boot']), rpl['n_inside'], nm(rpl['alt']),
          'finer' if rpl['alt'] < rl['alt'] else 'coarser', rl['n_inside'], nm(rl['alt']), ratio_str(rpl['boot'], rl['boot']),
          n_ref, rng_str(thin['iss']), nm(ra['alt']), rel_phrase(nanmedian(thin['iss']), ra['alt']), rng_str(thin['pic_ng1000']),
          'largely a comparison of how many localizations each reports' if abs(nanmedian(thin['iss']) / ra['alt'] - 1) <= 0.10 else 'not explained by the number of localizations alone'))
    P('* **Repeat precision** (per-axis, from consecutive frames): Picasso R0 %s nm, ISS %s nm, Picasso at net gradient 1000 %s nm, ISS without width-bound units %s nm (bootstrap intervals in the table above): %s.' % (
        nm(ra['repeat']['fit']['sigma_nm']), nm(rb['repeat']['fit']['sigma_nm']), nm(rl['repeat']['fit']['sigma_nm']), nm(rpl['repeat']['fit']['sigma_nm']),
        ('the ISS is less repeatable than Picasso R0' if rb['repeat']['fit']['sigma_nm'] > ra['repeat']['fit']['sigma_nm'] else 'the ISS is at least as repeatable as Picasso R0') +
        ' (the set of repeated molecules differs between the rows: Picasso R0 keeps only the localizations above its threshold, the brighter ones; the lower-threshold Picasso row is the like-for-like reference)'))
    P('* **Caveats.** Half-set resolutions of THIS sequence under the stated assumptions (pixel size 100 nm; sigma0 and the noise law from the data; the a priori ISS configuration with its circular kernel on an '
      'elongated PSF, a chi2 of %.1f and %.0f %% of its units at the lower width bound); no ground truth; nothing here is a detection accuracy; the seconds per frame were measured on a shared, heavily loaded machine (see above).\n' % (
          stats[:, 4].mean(), 100 * np.mean(sg_all <= lo_w + WIDTH_TOL)))
    P('## Reproduce\n')
    P('```\npython data/smlm/fetch_real.py                       # the zip, checksummed, into $SMLM_DATA_REAL\n'
      'SMLM_PICASSO_PYTHON=<interpreter with picasso> python experiments/exp_smlm_real.py --run --workers 2   # localizes (about 1 h on 2 busy cores), writes the raw and this file\n'
      'python experiments/exp_smlm_real.py --bench 50        # optional: single-process CPU timing of the ISS, added to the raw (the keys iss|bench, meta|bench)\n'
      'python experiments/exp_smlm_real.py --from-raw        # this file, from the raw alone (byte for byte)\n'
      'python -m pytest tests/test_frc.py tests/test_smlm_real.py -q\n```\n')
    return '\n'.join(L) + '\n'


def write_report(raw_path=RAW, md_path=MD, n_boot=N_BOOT, n_splits=N_SPLITS, log=print):
    with np.load(raw_path, allow_pickle=False) as z:
        raw = {k: z[k] for k in z.files}
    text = build_md(raw, n_boot, n_splits, file_sha256(raw_path), log=log)
    os.makedirs(os.path.dirname(md_path), exist_ok=True)
    with open(md_path, 'w') as f:
        f.write(text)
    log('wrote %s' % md_path)
    return text


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('--run', action='store_true', help='localize (Picasso, ISS), calibrate, write the raw and the md')
    ap.add_argument('--from-raw', action='store_true', help='regenerate the md from the raw alone')
    ap.add_argument('--calib-check', action='store_true', help='add the stability of the noise-law scale along the sequence to the raw (then the md is regenerated)')
    ap.add_argument('--bench', type=int, default=0, metavar='N', help='add the single-process CPU timing of the ISS on N evenly spaced frames to the raw (then the md is regenerated)')
    ap.add_argument('--frames', type=int, default=None, help='only the first N frames (smoke runs)')
    ap.add_argument('--workers', type=int, default=2, help='ISS worker processes (at most 2 on the shared machine)')
    ap.add_argument('--cache-dir', default=None, help='keep finished chunks of ISS frames here (resumable)')
    ap.add_argument('--boot', type=int, default=N_BOOT, help='bootstrap replicates')
    ap.add_argument('--splits', type=int, default=N_SPLITS, help='split seeds')
    ap.add_argument('--picasso-python', default=os.environ.get('SMLM_PICASSO_PYTHON'))
    ap.add_argument('--raw', default=RAW)
    ap.add_argument('--md', default=MD)
    a = ap.parse_args(argv)
    if not (a.run or a.from_raw or a.bench or a.calib_check):
        ap.error('give --run, --bench N, --calib-check or --from-raw')
    if a.run:
        run(a.frames, a.workers, a.cache_dir, a.picasso_python, a.raw)
    if a.calib_check:
        calib_check(a.raw)
    if a.bench:
        bench(a.raw, a.bench)
    write_report(a.raw, a.md, a.boot, a.splits)


if __name__ == '__main__':
    main()
