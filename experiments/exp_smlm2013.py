"""
The sliding inverse scale space and the Picasso ladder on the 2013 challenge's simulated 2D training sets (a second structure, a theoretical Gaussian PSF)
==================================================================================================================================

The four 2D simulated sets of the 2016 challenge that the earlier experiments use (``exp_smlm_picasso.py``, ``exp_smlm_iss.py``) are ONE structure (microtubules) at two noise
levels and two densities, with one experimentally derived PSF (a Gaussian core with a long halo). This experiment evaluates the same two methods, with the same harness, on the four
TRAINING sets of the 2013 challenge (``experiments/smlm2013.py``, data from ``data/smlm/fetch_2013.py``): a different structure (the 8 tubulins; bundles of 8 tubes), a different
pixel size (150 nm, 256 x 256 px; 100 nm, 64 x 64 px), a different camera (QE 1, 1 DN per electron, baseline 100) and a theoretical 2D Gaussian PSF of about 100 nm (no halo), on a
structured, time-varying autofluorescence background. ``smlm2013.py`` documents what the pages say, what the files say and what is assumed; the points that matter here:
the ground truth has ONE-FRAME activations (no emitter is on in two frames; the 2016 sets: 2.3 frames), the sets are 2401 / 2401 / 361 / 12,000 frames, and the densities are
41.6 emitters per 256 x 256 frame (Tubulins; about 2.6 per 64 x 64 tile, the density of the 2016 sparse sets), 6.75 (Bundled_Tubes_Long_Sequence) and 224.5 (Bundled_Tubes_High_Density,
one emitter per 18 pixels) per 64 x 64 frame.

Protocol (DECLARED BEFORE any test frame was scored; nothing is tuned or chosen on test frames)
------------------------------------------------------------------------------------------------
Development and test frame ranges: ``smlm2013.PROTOCOL`` (Tubulins: dev 0-39, test 200-299; Bundled long sequence: dev 0-299, test 400-899; Bundled high density: dev 0-29, test
100-159; 0-based stack indices). The long sequences are used in part only (the ISS costs seconds per frame; reasons in ``smlm2013.py``). Before the protocol was declared, the
following was looked at: the pages and the files (frame counts, the ground truth, 100 frames of each set for pixel statistics), Picasso on frames 0-199 of the long bundle set,
0-39 / 0-29 of the others (convention checks, candidate net gradients), the noise law and the PSF width on development frames, the ISS on ONE 64 x 64 tile of
frames 100-105 of Tubulins_I (timing and functioning), on frames 10-15 of the long bundle set and on frame 40 of the dense set (a path without early stop, to see what the entry
significances look like in a dense frame; frames 30-99 of the dense set and 40-199 of the Tubulins are outside both ranges), and the PSF file of each set (`psf.zip`). No frame of a
test range (Tubulins 200-299, long bundles 400-899, dense bundles 100-159) was looked at, and none was loaded for any other purpose than the test runs.

Both methods are tuned on the development frames and scored once on the test frames. Everything is scored with the harness (``smlm_common.score`` with the field of the set: Hungarian
matching at 250 nm, pooled Jaccard / RMSE / efficiency; the challenge setting of the 2016 harness (emitters within 450 nm of the border and the dimmest quarter of the whole set's
ground truth dropped after pairing) is reported as a second setting, its rules being from memory, see ``smlm_common.py``), with 95 % intervals of a bootstrap over blocks of 10
consecutive frames and paired differences (the same resampled blocks for both methods).

Picasso: the ladder of ``exp_smlm_picasso.py`` (its functions are imported; its grids are adapted below)
    R0  Picasso's defaults (minimum net gradient 5000, box 7, no padding, gausslq), with the camera of these sets passed to it (baseline 100, gain 1, sensitivity 1, QE 1, pixel size
        of the set: they change its photon column, not its identification or positions);
    R1E / R1J  the minimum net gradient alone tuned on the development frames (20 values) for efficiency E / Jaccard J;
    R2E / R2J  every option tuned: the structural grid of ``exp_smlm_picasso.STRUCT_GRID`` (3 x 3 x 3 x 2 x 2 = 108 settings: identification box, fit box, Gaussian filter, median padding,
        fit method) and, for each, the coordinate ascent of the post-filters (net gradient, minimum photons, minimum width, maximum width) for E and for J;
    +c  every rung also with the constant position offset (mean signed error of its matched development true positives) subtracted.
    GRIDS (declared here, before any run): the net-gradient grids are those of ``exp_smlm_picasso`` for its N2 sets (R1: ``GRID_R1['N2']``, 100-8000) extended at the top (28 values,
    100-15000) for the R2 line search; they were chosen from a low-threshold Picasso run on development frames of each set (true positives have net gradients from a few hundred to
    10,000 with a median of 2500-2900 in the sparse sets, false positives stay below 600-1600), the photon grid is the N1 grid of that module plus 2000 and 3000 (the dense set's
    photon counts are inflated by overlap), the width windows are the N1 windows scaled by 100 / 140 (the PSF of these sets is about 100 nm, that of the 2016 sets 140 nm).

The sliding inverse scale space: ``src/core/imaging.py`` through the pipeline of ``exp_smlm_iss.py`` (CONFIGS, KAPPAS, DEFAULT, the operating-point selection are imported), with
these changes, all forced by the sets and declared here:
  * the frame size is not 64 x 64: the Tubulin frames (256 x 256) are processed in tiles of 64 x 64 px cores with 8 px margins (72-80 px tiles; every detection whose centre is in the
    core of a tile is kept, and two detections of neighbouring tiles closer than the resolution limit r_ex are merged), which keeps the regime the pipeline was built for (a few emitters per
    tile, the same false-alarm count per unit area at the same kappa); the 64 x 64 sets are single tiles (the frame);
  * lengths tied to the PSF width scale with it: the entry width sigma0 (below), the bounds of the sliding widths (1.0, 3.2) px x sigma0 / 1.4, the resolution limit r_ex = 2.0 px x sigma0 / 1.4
    and the width bounds of the post-filter (2.0, 2.5, 3.2) px x sigma0 / 1.4 (all equal exp_smlm_iss's values at sigma0 = 1.4 px); the high-pass widths (3, 5 px) of CONFIGS and the lattice
    h = 0.5 px are in pixels as there;
  * sigma0 is the median width of a pixel-integrated Gaussian (plus a plane) fitted to isolated (no other emitter of the frame within 8 px) bright (at or above the median photons of the
    development frames) ground-truth emitters of the development frames, an 11 x 11 px window (the procedure of ``calibration.md``, labels used); the dense set borrows that of the long
    bundle set (it has no isolated emitters, and shares its PSF and camera), as the dense 2016 sets borrow the camera of their sparse set;
  * the noise law var = a x + b0 (x in photons = DN - 100): the slope a = 1.41 is the page's EM factor (it is checked, not tuned: the free slope fitted to the tubulin sets is 1.5-1.6),
    the constant b0 is calibrated on the development frames WITHOUT labels from the single-frame high-frequency noise (``noise_law``); the page's read-out sds (5, 10, 30) give 25, 100,
    900, the measured b0 are about 127 / 210 / 1016 (the dense set borrows that of the long bundle set); the background surface is the temporal median of the development frames
    smoothed by a 2 px Gaussian (the 2016 sets used a quadratic surface, which cannot follow the structured autofluorescence);
  * the path stops at the first entrant below min(KAPPAS) = 4 (the 2016 pipeline ran three entries further; the significance rule uses only the leading entries above kappa >= 4, so the
    result is identical for every kappa of the grid and the run is cheaper); at most 40 entries per tile.
  * CROWDED recipe, used for Bundled_Tubes_High_Density only (224 emitters in a 64 x 64 frame: the transitive cluster of a refit is the whole frame, and the first entrant below the
    significance stops the path after about 50 of them, see the diagnostic table in the results): refits only of the kernels within 4 px of an entrant, in a window of 10 px beyond them
    (``fit_radius``, ``rtrunc = 10 / reach``: the recipe of ``exp_astro_m2.py`` and ``exp_smlm_real.py``), no final polish of all kernels, at most 400 entries and NO early stop: the final
    block of the path is cleaned by the units' own significance (``imaging.unit_significance`` >= kappa, rule 'clean'). The protocol's early-stop rule ('snr') is stored from the same
    path and reported for the default only, as a diagnostic.
  * The ISS time per frame is measured in the test runs on a shared machine (one core per process): wall-clock and process CPU seconds are both reported.
  * Two descriptions that are NOT rungs of the protocol and were added after the test runs, not used to change anything: the recall / RMSE of the matched pairs by the distance of an emitter to
    its nearest neighbour (every set, test frames), and, for the dense set only, a diagnostic run of the E-tuned configuration with a final polish of the cleaned block on four development frames
    (``part_diag``; the crowded recipe has none), to see whether that omission costs accuracy.

DECLARED a priori default of the ISS, run once on the test frames with no choice made on any frames except the calibrations above: configuration 'hp5h' (Gaussian core 0.75 + halo 0.25 at 3
times the width, high-pass 5 px), kappa = 5, no width bound, no position offset correction (the recipe of the set: tiles / the 'snr' rule for the sparse sets, the crowded recipe with the
'clean' rule for the dense set). The OPERATING POINTS E and J are chosen on the development frames from the 4 configurations of exp_smlm_iss.CONFIGS x the 7 kappas x the 3 width bounds
(84 points, each with the position offset estimated on the same frames and removed), by the efficiency and by the Jaccard index (``exp_smlm_iss.choose_op``); the second round of the 2016
development search (entry width and resolution limit variants) is NOT run (cost). The ISS search is far smaller than Picasso's (84 points against ~10,000 scored configurations per set).

Usage:  python experiments/exp_smlm2013.py --run [--sets S ...] [--only meta cal psf pic isdev diag istest] [--workers 2] [--parts-dir DIR]
                                                                                                     (computes the parts, assembles smlm2013_raw.npz, writes smlm2013.md)
        python experiments/exp_smlm2013.py --assemble                                              (parts -> raw only)
        python experiments/exp_smlm2013.py --from-raw [--boot 1000]                               (regenerates smlm2013.md from the raw alone, byte for byte)
Environment: SMLM_DATA_2013, SMLM_PICASSO_PYTHON (an interpreter with picasso), SMLM2013_OUT (directory of the raw and the md, default results/smlm).
"""
import argparse
import collections
import hashlib
import io
import json
import multiprocessing
import os
import platform
import sys
import time
import zipfile

os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import numpy as np                                                  # noqa: E402
from scipy import ndimage as ndi                                    # noqa: E402
from scipy.optimize import least_squares                            # noqa: E402
from scipy.special import erf                                       # noqa: E402
from scipy.spatial import cKDTree                                   # noqa: E402

HERE = os.path.abspath(os.path.dirname(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
sys.path.insert(0, REPO)
import smlm_common as C                                             # noqa: E402
import smlm2013 as S                                                # noqa: E402
import exp_smlm_picasso as PIC                                      # noqa: E402
import exp_smlm_iss as ISS                                          # noqa: E402
from src.core import imaging as IM                                  # noqa: E402

OUT = os.environ.get('SMLM2013_OUT', os.path.join(REPO, 'results', 'smlm'))
RAW = os.path.join(OUT, 'smlm2013_raw.npz')
MD = os.path.join(OUT, 'smlm2013.md')
PARTS = os.path.join(OUT, 'smlm2013_parts')
SETS = S.SETS
BOOT, BLOCK = 1000, C.DEFAULT_BLOCK
CHALLENGE = dict(border_nm=450.0, keep=0.75)

# ---- the ISS declarations ------------------------------------------------------------------------------------------------------------
KAPPAS = ISS.KAPPAS                                                 # (4, 4.5, 5, 6, 7, 8, 10)
CONFIGS = ISS.CONFIGS
DEFAULT = dict(ISS.DEFAULT)                                         # hp5h, kappa 5
OPS = ISS.OPS                                                       # ('E', 'J')
K_MIN = min(KAPPAS)
SIGMA_REF = 1.4                                                     # px: sigma0 of the 2016 sets, at which the lengths below equal those of exp_smlm_iss
WIDTH_REL = (1.0 / SIGMA_REF, 3.2 / SIGMA_REF)                      # bounds of the sliding widths in units of sigma0
REX_REL = ISS.R_EX / SIGMA_REF                                      # resolution limit in units of sigma0
SMAX_REL = tuple(v / SIGMA_REF for v in ISS.SMAX)                   # width post-filter bounds in units of sigma0 (the last: the bound of the widths = no filter)
LAT_H = ISS.LAT_H                                                   # px
EM_FACTOR = S.EM_FACTOR                                             # a of the noise law (page)
BG_SMOOTH = 2.0                                                     # px: Gaussian smoothing of the temporal-median background
NL_SMOOTH, NL_BINS, NL_STEP = 3.0, 16, 1                            # the noise law: Gaussian of the high-pass residual, quantile bins, every NL_STEP-th development frame
S0_HALF, S0_ISO_PX, S0_BRIGHT = 5, 8.0, 0.5                         # sigma0: 11 x 11 window, isolation radius, quantile of the photons above which an emitter is bright
RECIPES = {
    'sparse': dict(tile=(64, 8), cap=40, fit_radius=None, rtrunc=4.0, refine_box=1.0, refine_rounds=2, refine_nfev=100, snr_min=K_MIN, kinds=('snr',)),
    'crowded': dict(tile=None, cap=400, fit_radius=4.0, rtrunc=None, refine_box=0.0, refine_rounds=2, refine_nfev=100, snr_min=None, kinds=('clean', 'snr')),
    # diagnostic only (``part_diag``): the crowded recipe with a final polish of the cleaned block (one round, at most 30 evaluations): NOT a rung of the protocol
    'crowded_polish': dict(tile=None, cap=400, fit_radius=4.0, rtrunc=None, refine_box=1.0, refine_rounds=1, refine_nfev=30, snr_min=None, kinds=('clean',)),
}
DIAG_FRAMES = (3, 11, 19, 27)                                       # development frames of the dense set on which the polish is looked at
RECIPE_OF = {'Tubulins_I': 'sparse', 'Tubulins_II': 'sparse', 'Bundled_Tubes_Long_Sequence': 'sparse', 'Bundled_Tubes_High_Density': 'crowded'}
RULE_KIND = {'sparse': 'snr', 'crowded': 'clean'}                   # the rule family of the operating points of a recipe

# ---- the Picasso grids (declared above) ---------------------------------------------------------------------------------------------
GRID_R1 = list(PIC.GRID_R1['N2'])
FILTER_GRID = collections.OrderedDict([
    ('min_net_gradient', [100, 125, 150, 175, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 1000, 1250, 1500, 2000, 2500, 3000, 3500, 4000, 5000, 6000, 8000, 10000, 12000, 15000]),
    ('min_photons', [None, 100, 200, 300, 400, 600, 800, 1000, 1500, 2000, 3000]),
    ('min_sigma_nm', [None, 50, 60, 70, 80, 90]),
    ('max_sigma_nm', [None, 140, 175, 210, 250, 300, 350, 500])])
PIC_NG_LOW = float(min(FILTER_GRID['min_net_gradient']))
CONV_NG = 1000.0                                                    # the fixed net gradient of the convention check (a few frames of the development range)
CONV_FRAMES = 20


def camera(name):
    """The Picasso camera options of a set (what the pages give: baseline 100 DN, 1 DN per electron, QE 1; pixel size)."""
    return dict(baseline=S.BASELINE_DN, gain=1.0, sensitivity=1.0, qe=1.0, pixel_nm=S.spec(name)['pixel_nm'])


# ==================================================================================================================================
# calibration (labels not used, except sigma0)
# ==================================================================================================================================
def background_surface(dev_photons, smooth=BG_SMOOTH):
    """The background of a set: the per-pixel temporal median of the development frames (n, ny, nx) in photons, smoothed by a Gaussian of ``smooth`` px."""
    return ndi.gaussian_filter(np.median(dev_photons, axis=0), smooth, mode='nearest')


_MED_CHI2 = 0.4549364231195724                                      # median of the chi-square distribution with one degree of freedom


def noise_law(dev_photons, a=EM_FACTOR, step=NL_STEP, nbins=NL_BINS, smooth=NL_SMOOTH):
    """The constant b0 of the noise law var = a x + b0 (x in photons), from the single-frame high-frequency noise of the development frames, without labels.

    The residual r = x - G x of every ``step``-th frame (G: a Gaussian of ``smooth`` px, mirrored borders) has, for noise that is white and independent of the neighbours, the variance
    (1 - 2 g0 + sum g^2) var with g the kernel (g0 its centre weight); its ROBUST variance (median(r^2) / 0.4549, resistant to the few pixels of emitters) is taken in ``nbins`` quantile bins
    of the smooth level G x. Returns (rows, b0): rows (nbins, 3) = (mean level, robust variance, pixels) and b0 = the median over the bins at or below the 90th percentile of the level of
    (variance - a max(level, 0)) (the bins of the brightest pixels of a sparse frame hold the emitters). The slope is the page's EM factor ``a``; ``free_slope`` fits it."""
    lev, res = [], []
    for f in np.asarray(dev_photons)[::step]:
        s = ndi.gaussian_filter(f, smooth, mode='mirror')
        lev.append(s.ravel())
        res.append((f - s).ravel())
    lev, res = np.concatenate(lev), np.concatenate(res)
    g0, sg2 = 1.0 / (2 * np.pi * smooth ** 2), 1.0 / (4 * np.pi * smooth ** 2)
    cf = 1.0 - 2 * g0 + sg2
    q = np.quantile(lev, np.linspace(0.0, 1.0, nbins + 1))
    rows = []
    for lo, hi in zip(q[:-1], q[1:]):
        m = (lev >= lo) & (lev <= hi)
        rows.append((float(lev[m].mean()), float(np.median(res[m] ** 2) / (_MED_CHI2 * cf)), int(m.sum())))
    rows = np.array(rows)
    use = rows[:, 0] <= np.quantile(lev, 0.9)
    b0 = float(np.median(rows[use, 1] - a * np.maximum(rows[use, 0], 0.0)))
    return rows, b0


def free_slope(rows):
    """Weighted least-squares line var = a x + b0 through the bins of ``noise_law`` (all bins, a check of the page's slope; meaningless where the levels of the bins span little, e.g. the
    long bundle set, or where emitters dominate the top bins): returns (a, b0, span) with span = the range of the bin levels."""
    x, v, w = np.maximum(rows[:, 0], 0.0), rows[:, 1], rows[:, 2]
    A = np.c_[x, np.ones(len(x))]
    sol = np.linalg.lstsq(A * np.sqrt(w)[:, None], v * np.sqrt(w), rcond=None)[0]
    return float(sol[0]), float(sol[1]), float(rows[:, 0].max() - rows[:, 0].min())


def _profile(ax, c, s):
    k = 1.0 / (np.sqrt(2.0) * s)
    return 0.5 * (erf((ax + 0.5 - c) * k) - erf((ax - 0.5 - c) * k))


def fit_width(patch, half, x0=0.0, y0=0.0):
    """Least-squares fit of a pixel-integrated Gaussian plus a plane to a (2 half + 1)^2 patch (photons; centre coordinates relative to the central pixel): returns (parameters
    (N, x, y, sigma, c0, cx, cy), rms residual). The plane absorbs the slope of the autofluorescence background."""
    ax = np.arange(-half, half + 1, dtype=float)
    U, V = np.meshgrid(ax, ax)

    def model(t):
        N, xc, yc, s, c0, cx, cy = t
        return c0 + cx * U + cy * V + N * np.outer(_profile(ax, yc, s), _profile(ax, xc, s))

    med = float(np.median(patch))
    t0 = [max(float(patch.sum()) - med * patch.size, 100.0), x0, y0, 1.0, med, 0.0, 0.0]
    lo = [0.0, -3.0, -3.0, 0.3, -np.inf, -np.inf, -np.inf]
    hi = [np.inf, 3.0, 3.0, 4.0, np.inf, np.inf, np.inf]
    r = least_squares(lambda t: (model(t) - patch).ravel(), t0, bounds=(lo, hi), x_scale='jac')
    return r.x, float(np.sqrt(np.mean(r.fun ** 2)))


def sigma0_estimate(frames_photons, gt_dev, first_frame, pixel_nm, half=S0_HALF, iso_px=S0_ISO_PX, bright=S0_BRIGHT):
    """The PSF width of a set from its development frames, labels used (the procedure of ``calibration.md``): the pixel-integrated Gaussian + plane (``fit_width``) fitted in an
    (2 half + 1)^2 window to every ground-truth emitter of the frames that is bright (photons at or above the ``bright`` quantile of the frames' ground truth), isolated (no other
    ground-truth emitter of its frame within ``iso_px``) and has its whole window inside the frame. Returns the array (m, 7) of (sigma px, dx nm, dy nm (fit - truth), fitted photons N,
    ground-truth photons, |z| nm, rms) - sigma0 is the median of its first column."""
    n, ny, nx = frames_photons.shape
    thr = np.quantile(gt_dev['photons'], bright)
    out = []
    for i in range(n):
        gf = gt_dev[gt_dev['frame'] == first_frame + i]
        if len(gf) == 0:
            continue
        xs, ys = gf['x_nm'] / pixel_nm - 0.5, gf['y_nm'] / pixel_nm - 0.5
        for j in range(len(gf)):
            if gf['photons'][j] < thr:
                continue
            d = np.hypot(xs - xs[j], ys - ys[j])
            d[j] = np.inf
            if d.min() < iso_px:
                continue
            ix, iy = int(round(xs[j])), int(round(ys[j]))
            if ix < half or ix > nx - 1 - half or iy < half or iy > ny - 1 - half:
                continue
            t, rms = fit_width(frames_photons[i, iy - half:iy + half + 1, ix - half:ix + half + 1], half, xs[j] - ix, ys[j] - iy)
            out.append((t[3], (ix + t[1] + 0.5) * pixel_nm - gf['x_nm'][j], (iy + t[2] + 0.5) * pixel_nm - gf['y_nm'][j], t[0], gf['photons'][j], abs(gf['z_nm'][j]), rms))
    return np.array(out).reshape(-1, 7)


def calibrate_set(name, log=print):
    """Calibration of a set from its development frames (a dense set from those of its sparse set, ``smlm2013.CAL_FROM``): the background surface, the noise law (b0 and the bins), the PSF
    width sigma0 with its sample. Returns a dict of arrays / scalars (the part 'cal')."""
    src = S.CAL_FROM.get(name, name)
    a, b = S.PROTOCOL[src]['dev']
    ph = S.to_photons(S.load_frames(S.data_path(src, 'stack'), a, b))
    gt = S.gt_in_frames(S.load_gt(S.data_path(src, 'gt')), a, b)
    bg = background_surface(ph)
    rows, b0 = noise_law(ph)
    af, bf, span = free_slope(rows)
    sig = sigma0_estimate(ph, gt, a, S.spec(src)['pixel_nm'])
    log('%s (calibration from %s frames %d-%d): background mean %.2f photons (range %.1f .. %.1f), b0 = %.1f (a = %.2f), free slope %.2f (span %.0f), sigma0 = %.3f px from %d emitters'
        % (name, src, a, b - 1, bg.mean(), bg.min(), bg.max(), b0, EM_FACTOR, af, span, np.median(sig[:, 0]), len(sig)))
    return dict(bg=bg, a=np.array(EM_FACTOR), b0=np.array(b0), nl_rows=rows, nl_free=np.array([af, bf, span]), sigma0=np.array(float(np.median(sig[:, 0]))), sigma_sample=sig,
                src=np.array(src), n_dev=np.array(b - a))


# ==================================================================================================================================
# the sliding inverse scale space on one frame (tiles), for the sets of 2013
# ==================================================================================================================================
_W = {}


def tile_grid(n, core, margin):
    """The tiles along one axis of n pixels: a list of (first, stop, core_index) with cores of ``core`` px (the last may be shorter) widened by ``margin`` px inside the frame. A detection with
    centre (pixel-index coordinates, pixel centres at integers) at u = x + 1/2 belongs to the tile with index clip(floor(u / core))."""
    nt = max(1, -(-n // core))
    return [(max(0, k * core - margin), min(n, (k + 1) * core + margin), k) for k in range(nt)]


def owner_index(x, core, nt):
    """Tile index along an axis of the pixel-index coordinate x (see ``tile_grid``)."""
    return np.clip(np.floor((np.asarray(x) + 0.5) / core).astype(int), 0, nt - 1)


def _init(name, cfg_name, cal, opts):
    cfg = ISS.get_cfg(cfg_name)
    sig0 = float(cal['sigma0'])
    rec = RECIPES[opts.get('recipe') or RECIPE_OF[name]]
    _W.clear()
    _W.update(name=name, cfg=cfg, sigma0=sig0, rec=rec, bg=np.asarray(cal['bg'], float), a=float(cal['a']), b0=float(cal['b0']), px=S.spec(name)['pixel_nm'], psf={},
              widths=(WIDTH_REL[0] * sig0, WIDTH_REL[1] * sig0), r_ex=REX_REL * sig0, kinds=opts['kinds'], size=S.spec(name)['size'])


def _psf(ny, nx):
    key = (ny, nx)
    if key not in _W['psf']:
        cfg, s0 = _W['cfg'], _W['sigma0']
        P = IM.PixelPSF(ny, nx, s0, parts=cfg['parts'], hp=cfg['hp'] or None)
        _W['psf'][key] = (P, IM.PixelPSF(ny, nx, s0, parts=cfg['parts']), IM.Lattice(P, h=LAT_H, margin=1.0))
    return _W['psf'][key]


def _rules(kinds):
    return [(k, kap) for k in kinds for kap in KAPPAS]


def run_tile(frame, bg):
    """The ISS on one tile (photons, ny x nx; bg the background of the tile): returns (dict rule -> (m, 4) detections (x px, y px in the tile, photons, sigma px), n_entries, capped)."""
    ny, nx = frame.shape
    P, Pp, L = _psf(ny, nx)
    rec, a, b0, widths = _W['rec'], _W['a'], _W['b0'], _W['widths']

    def sw_fn(c, w):
        m = Pp.synth(c, w) if len(w) else 0.0
        return 1.0 / np.sqrt(a * (np.maximum(bg, 0.0) + np.maximum(m, 0.0)) + b0)

    y0 = frame - bg
    y = P.highpass(y0 - np.median(y0))
    sw0 = sw_fn(np.zeros((0, 3)), np.zeros(0))
    rt = rec['rtrunc'] if rec['rtrunc'] is not None else 10.0 / P.reach(widths[1])
    R = IM.sliding_path(P, y, sw0, L, max_entries=rec['cap'], sw_fn=sw_fn, sw_on='params', sig_bounds=widths, tol=1e-8, passes=1, rtrunc=rt, r_ex=_W['r_ex'], snr_min=rec['snr_min'],
                        snr_extra=0, jump=None, extra=6, order='dual', fit_radius=rec['fit_radius'])
    out, cache = {}, {}
    for kind in _W['kinds']:
        for kap in KAPPAS:
            if kind == 'snr':
                i = IM.stop_snr(R['hist'], R['taus'], R['snr'], kap)[0]
            else:
                i = len(R['hist']) - 1
            h = R['hist'][i]
            cc, ww = h['c'], h['w']
            if kind == 'clean' and len(ww):
                keep = IM.unit_significance(P, sw0, cc, ww) >= kap
                cc, ww = cc[keep], ww[keep]
            key = (kind, i, kap if kind == 'clean' else None)
            if key not in cache:
                if rec['refine_box'] and len(ww):
                    cc, ww, _ = IM.refine(P, y, sw0, cc, ww, sig_bounds=widths, box=rec['refine_box'], sw_fn=sw_fn, sw_on='params', rounds=rec['refine_rounds'], max_nfev=rec['refine_nfev'])
                cache[key] = np.c_[cc[:, 0], cc[:, 1], ww, cc[:, 2]] if len(ww) else np.zeros((0, 4))
            out[(kind, kap)] = cache[key]
    capped = bool(len(R['taus']) >= rec['cap'] and (rec['snr_min'] is None or np.all(R['snr'] >= rec['snr_min'])))
    return out, len(R['taus']), capped


def merge_tiles(parts, r_ex):
    """Detections of one rule from all tiles of a frame ((m, 4) arrays in frame pixel-index coordinates, in tile order) with two detections of different tiles closer than ``r_ex`` px
    merged (the one of the larger photons kept). Returns the array in the order of the tiles."""
    if not parts:
        return np.zeros((0, 4))
    tid = np.concatenate([np.full(len(p), k) for k, p in enumerate(parts)])
    d = np.concatenate(parts)
    if len(d) < 2 or len(parts) < 2:
        return d
    keep = np.ones(len(d), bool)
    for i, j in sorted(cKDTree(d[:, :2]).query_pairs(r_ex)):
        if tid[i] != tid[j] and keep[i] and keep[j]:
            keep[i if d[i, 2] < d[j, 2] else j] = False
    return d[keep]


def run_frame(args):
    """One frame (photons, the full field): the path on every tile, the detections of every rule (kind, kappa) as (m, 4) arrays x_nm, y_nm, photons, sigma_px (corner convention, the
    centres inside the cores of their tiles, neighbouring duplicates merged). Returns (dict rule -> array, entries, capped tiles, wall seconds, process CPU seconds)."""
    k, frame = args
    t0, c0_ = time.perf_counter(), time.process_time()
    ny, nx = frame.shape
    rec, px, bg = _W['rec'], _W['px'], _W['bg']
    tl = [((0, ny, 0), (0, nx, 0))] if rec['tile'] is None else [(r, c) for r in tile_grid(ny, *rec['tile']) for c in tile_grid(nx, *rec['tile'])]
    core = rec['tile'][0] if rec['tile'] else max(ny, nx)
    nty, ntx = (1, 1) if rec['tile'] is None else (len(tile_grid(ny, *rec['tile'])), len(tile_grid(nx, *rec['tile'])))
    per = {key: [] for key in _rules(_W['kinds'])}
    n_ent = n_cap = 0
    for (r0, r1, ty), (c0, c1, tx) in tl:
        out, ne, cap = run_tile(frame[r0:r1, c0:c1], bg[r0:r1, c0:c1])
        n_ent += ne
        n_cap += int(cap)
        for key, d in out.items():
            if len(d):
                d = d.copy()
                d[:, 0] += c0
                d[:, 1] += r0
                if rec['tile'] is not None:
                    own = (owner_index(d[:, 0], core, ntx) == tx) & (owner_index(d[:, 1], core, nty) == ty)
                    d = d[own]
            per[key].append(d)
    res = {}
    for key, parts in per.items():
        d = merge_tiles(parts, _W['r_ex']) if len(parts) > 1 else (parts[0] if parts else np.zeros((0, 4)))
        res[key] = np.c_[(d[:, 0] + 0.5) * px, (d[:, 1] + 0.5) * px, d[:, 2], d[:, 3]] if len(d) else np.zeros((0, 4))
    return res, n_ent, n_cap, time.perf_counter() - t0, time.process_time() - c0_


def run_frames(name, frames_photons, cfg_name, cal, workers=1, recipe=None):
    """The ISS (recipe of the set, configuration ``cfg_name``) on frames (n, ny, nx) in photons. Returns (dets[rule][frame] arrays, entries per frame, capped tiles per frame, wall seconds per
    frame, process CPU seconds per frame); rules (kind, kappa) as ``_rules``."""
    kinds = RECIPES[recipe or RECIPE_OF[name]]['kinds']
    jobs = list(enumerate(frames_photons))
    opts = dict(kinds=kinds, recipe=recipe)
    if workers <= 1:
        _init(name, cfg_name, cal, opts)
        res = [run_frame(j) for j in jobs]
    else:
        with multiprocessing.get_context('fork').Pool(workers, initializer=_init, initargs=(name, cfg_name, cal, opts)) as pool:
            res = pool.map(run_frame, jobs, chunksize=1)
    rules = _rules(kinds)
    dets = {r: [x[0][r] for x in res] for r in rules}
    return dets, np.array([x[1] for x in res]), np.array([x[2] for x in res]), np.array([x[3] for x in res]), np.array([x[4] for x in res])


# ==================================================================================================================================
# raw storage: per-set parts (npz files, resumable) assembled into one deterministic file
# ==================================================================================================================================
SMAX_TOL = 1e-4                                                     # px: tolerance of the width post-filter (the widths are stored as float32)


def pack(dets):
    """List of per-frame (k, 4) arrays -> (offsets (F + 1,) int64, flat (M, 4) float32)."""
    off = np.r_[0, np.cumsum([len(d) for d in dets])].astype(np.int64)
    flat = np.concatenate([np.asarray(d).reshape(-1, 4) for d in dets]).astype(np.float32) if off[-1] else np.zeros((0, 4), np.float32)
    return off, flat


def unpack(off, flat):
    """Inverse of ``pack``: the per-frame arrays as float64."""
    return [np.asarray(flat[off[i]:off[i + 1]], dtype=np.float64) for i in range(len(off) - 1)]


def post(dets, smax):
    """Detections whose fitted width (column 3, px) is at most smax (``ISS.post`` with a tolerance for the float32 storage of the widths)."""
    return [d[d[:, 3] <= smax + SMAX_TOL] if len(d) else d for d in dets]


def save_npz(path, arrays):
    """A compressed npz that is the same bytes for the same arrays: sorted keys, a fixed zip date, deflate (``smlm_common.save_npz_deterministic`` without the compression)."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    tmp = path + '.tmp'
    with zipfile.ZipFile(tmp, 'w', compression=zipfile.ZIP_DEFLATED) as z:
        for k in sorted(arrays):
            zi = zipfile.ZipInfo(k + '.npy', date_time=(1980, 1, 1, 0, 0, 0))
            zi.compress_type = zipfile.ZIP_DEFLATED
            zi.external_attr = 0o644 << 16
            buf = io.BytesIO()
            np.lib.format.write_array(buf, np.asanyarray(arrays[k]), allow_pickle=False)
            z.writestr(zi, buf.getvalue(), compress_type=zipfile.ZIP_DEFLATED, compresslevel=6)
    os.replace(tmp, path)


def load_npz(path):
    with np.load(path, allow_pickle=False) as z:
        return {k: z[k] for k in z.files}


def part_file(name, part, parts_dir=None):
    return os.path.join(parts_dir or PARTS, '%s.%s.npz' % (name, part))


def have_part(name, part, parts_dir=None):
    return os.path.exists(part_file(name, part, parts_dir))


def save_part(name, part, arrays, parts_dir=None):
    save_npz(part_file(name, part, parts_dir), {k: np.asarray(v) for k, v in arrays.items()})


def load_part(name, part, parts_dir=None):
    return load_npz(part_file(name, part, parts_dir))


def part_names(name, parts_dir=None):
    d = parts_dir or PARTS
    pre = name + '.'
    return sorted(f[len(pre):-4] for f in os.listdir(d) if f.startswith(pre) and f.endswith('.npz')) if os.path.isdir(d) else []


def load_set_parts(name, parts_dir=None):
    """All parts of a set as one dict with keys 'part|key' (the keys of the assembled raw without the set name)."""
    out = {}
    for part in part_names(name, parts_dir):
        for k, v in load_part(name, part, parts_dir).items():
            out['%s|%s' % (part, k)] = v
    return out


def assemble(sets=SETS, parts_dir=None, raw_path=RAW, log=print):
    """Merge the part files of the sets into the raw file (keys 'set|part|key'); a set without parts is skipped. Deterministic (same parts, same bytes)."""
    have = [s for s in sets if part_names(s, parts_dir)]
    if not have:
        raise SystemExit('no part files of %s in %s: nothing to assemble (the raw file is not touched)' % (', '.join(sets), parts_dir or PARTS))
    raw = {'README': np.array(RAW_README), 'sets': np.array(json.dumps(have))}
    for name in sets:
        for k, v in load_set_parts(name, parts_dir).items():
            raw['%s|%s' % (name, k)] = v
    save_npz(raw_path, raw)
    log('assembled %s (%d arrays, %.1f MB)' % (raw_path, len(raw), os.path.getsize(raw_path) / 1e6))


def load_raw(path=RAW):
    """The raw file as {set: {'part|key': array}} plus the entries 'README' and 'sets'."""
    flat = load_npz(path)
    out = {k: flat[k] for k in ('README', 'sets')}
    for k, v in flat.items():
        if '|' in k:
            name, rest = k.split('|', 1)
            out.setdefault(name, {})[rest] = v
    return out


RAW_README = ('Raw results of experiments/exp_smlm2013.py. Keys "<set>|<part>|<key>": part meta (ground truth of the development and test frames, set statistics, data digests), cal (calibration: '
              'background, noise law, PSF width), pic (the Picasso ladder: tuning curves, search results, chosen options, detections of rows R0 R1E R1J R2E R2J on the development and test frames as '
              '"<split>|<row>|off" (F + 1 offsets) and "<split>|<row>|flat" (M x 4 float32: x_nm, y_nm, photons, sx_nm)), conv (convention checks), isdev.<config> and istest.<config> (the sliding ISS: '
              'detections of every rule "<kind><kappa>|off", "<kind><kappa>|flat" (M x 4 float32: x_nm, y_nm, photons, sigma_px), entries, capped tiles, seconds per frame). Positions in nm, corner '
              'convention, frames 0-based stack indices. results/smlm/smlm2013.md is generated from this file alone (python experiments/exp_smlm2013.py --from-raw).')


def _jsonable(o):
    return json.dumps(o, sort_keys=True, default=C._json_default)


# ==================================================================================================================================
# parts
# ==================================================================================================================================
def _sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def part_meta(name, log=print):
    """Facts of a set that the report needs from the data: ground truth of the development and test frames, whole-set statistics, the frame-base evidence, data digests."""
    (a, b), (c, d) = S.PROTOCOL[name]['dev'], S.PROTOCOL[name]['test']
    gt = S.load_gt(S.data_path(name, 'gt'))
    rl = C.run_lengths(gt)
    n_frames = S.n_frames_in(S.data_path(name, 'stack'))
    n_ex = min(12, b - a)                                            # development frames only
    ex = S.excess_at_ground_truth(S.load_frames(S.data_path(name, 'stack'), a, a + n_ex), gt, a, pixel_nm=S.spec(name)['pixel_nm'], min_photons=float(np.median(gt['photons'])))
    book = {}
    try:
        book = json.load(open(os.path.join(REPO, 'data', 'smlm', 'checksums_2013.json')))
    except Exception:
        pass
    sha = {k: dict(sha256=_sha256(S.data_path(name, k)), recorded=book.get('%s/%s' % (name, f), {}).get('sha256'))
           for k, f in (('stack', 'sequence.zip'), ('gt', 'fluorophores.zip'))}
    arrays = dict(gt_dev=S.gt_in_frames(gt, a, b), gt_test=S.gt_in_frames(gt, c, d), n_gt=len(gt), n_frames=n_frames, pmin=C.photon_threshold(gt, CHALLENGE['keep']),
                  runs=np.array([rl['runs'], rl['mean']]), photon_q=np.quantile(gt['photons'], [0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99]), z_range=np.array([gt['z_nm'].min(), gt['z_nm'].max()]),
                  z_median=np.median(np.abs(gt['z_nm'])), per_frame=np.bincount(gt['frame'], minlength=n_frames), frame_ranges=np.array([a, b, c, d]),
                  excess_offsets=np.array(sorted(ex)), excess_mean=np.array([ex[o][0] for o in sorted(ex)]), excess_n=np.array([ex[o][1] for o in sorted(ex)]),
                  sha=np.array(_jsonable(sha)), versions=np.array(_jsonable(dict(python=platform.python_version(), numpy=np.__version__))))
    log('%s meta: %d frames, %d ground-truth rows, dev %d / test %d emitters' % (name, n_frames, len(gt), len(arrays['gt_dev']), len(arrays['gt_test'])))
    return arrays


def part_cal(name, log=print):
    cal = calibrate_set(name, log)
    return {k: np.asarray(v) for k, v in cal.items()}


def get_cal(R):
    """The calibration dict of a set from its raw part 'cal'."""
    return dict(bg=R['cal|bg'], a=float(R['cal|a']), b0=float(R['cal|b0']), sigma0=float(R['cal|sigma0']))


def part_is(name, split, cfg_name, R, workers=1, log=print):
    """The ISS with configuration ``cfg_name`` on the development or test frames: detections of every rule of the recipe, with entries / capped tiles / seconds per frame."""
    a, b = S.PROTOCOL[name][split]
    frames = S.to_photons(S.load_frames(S.data_path(name, 'stack'), a, b))
    t0 = time.time()
    dets, nent, ncap, secs, cpu = run_frames(name, frames, cfg_name, get_cal(R), workers)
    arrays = {'frames': np.array([a, b]), 'entries': nent, 'capped': ncap, 'secs': secs, 'cpu': cpu}
    if split == 'test':
        arrays['bg_test'] = background_surface(frames).astype(np.float32)                 # label-free: the background of the test frames, to compare with the calibration's (stability)
    for (kind, kap), ds in dets.items():
        off, flat = pack(ds)
        arrays['%s%g|off' % (kind, kap)], arrays['%s%g|flat' % (kind, kap)] = off, flat
    gt = S.gt_in_frames(R['meta|gt_%s' % split], a, b) if 'meta|gt_%s' % split in R else None
    kind = RULE_KIND[RECIPE_OF[name]]
    msg = ''
    if gt is not None and (kind, 5.0) in dets:
        r = S.score(dets[(kind, 5.0)], gt, np.arange(a, b), name)
        msg = '; kappa 5: J %.3f RMSE %.1f E %.1f' % (r['jaccard'], r['rmse_nm'], r['efficiency'])
    log('ISS %s %s %s: %d frames in %.0f s (%.2f s wall, %.2f s CPU per frame, %.1f entries, %d capped tiles)%s' % (name, split, cfg_name, b - a, time.time() - t0, secs.mean(), cpu.mean(),
                                                                                                                   nent.mean(), ncap.sum(), msg))
    return arrays


# ---------------------------------------------------------------- the operating points of the ISS (development frames)
def cfg_list(R):
    have = {k[len('isdev.'):].split('|')[0] for k in R if k.startswith('isdev.')}
    return [c for c in CONFIGS if c in have]


def is_dets(R, split, cfg, kind, kap):
    key = 'is%s.%s|%s%g|' % (split, cfg, kind, kap)
    return unpack(R[key + 'off'], R[key + 'flat'])


def op_table(R, name):
    """Development-frame metrics of every operating point (config, kappa, width bound) of the rule family of the recipe, with the position offset estimated on the same frames and removed:
    a list of dict(cfg, kappa, smax, off, m = [J, P, R, RMSE, E, J_challenge, E_challenge]) (``exp_smlm_iss.op_table`` for these sets)."""
    kind = RULE_KIND[RECIPE_OF[name]]
    a, b = (int(v) for v in R['isdev.%s|frames' % cfg_list(R)[0]])
    idx = np.arange(a, b)
    gt, pmin, sig0 = R['meta|gt_dev'], float(R['meta|pmin']), float(R['cal|sigma0'])
    rows = []
    for cfg in cfg_list(R):
        for kap in KAPPAS:
            base = is_dets(R, 'dev', cfg, kind, kap)
            for rel in SMAX_REL:
                smax = rel * sig0
                dets = post(base, smax)
                off = ISS.offset_of(dets, gt, idx)
                dc = ISS.shift(dets, off)
                rows.append(dict(cfg=cfg, kappa=kap, smax=smax, off=off,
                                 m=ISS.metrics_row(S.score(dc, gt, idx, name), S.score(dc, gt, idx, name, border_nm=CHALLENGE['border_nm'], photon_min=pmin))))
    return rows


def select_ops(R, name):
    rows = op_table(R, name)
    return rows, {w: ISS.choose_op(rows, w) for w in OPS}


def op_test_dets(R, name, o):
    """Test-frame detections of an operating point: its configuration's run at its kappa, width-filtered, shifted by its development offset."""
    kind = RULE_KIND[RECIPE_OF[name]]
    return ISS.shift(post(is_dets(R, 'test', o['cfg'], kind, o['kappa']), o['smax']), o['off'])


def part_diag(name, R, workers=1, log=print):
    """DIAGNOSTIC of the dense set (not a rung of the protocol): the configuration of the E-tuned operating point run with the crowded recipe PLUS a final polish of the cleaned block
    (``RECIPES['crowded_polish']``) on the development frames ``DIAG_FRAMES``, to see whether the missing final polish of the crowded recipe costs position accuracy."""
    rows, ops = select_ops(R, name)
    o = ops['E']
    a, b = S.PROTOCOL[name]['dev']
    frames = S.to_photons(S.load_frames(S.data_path(name, 'stack'), a, b))[list(DIAG_FRAMES)]
    t0 = time.time()
    dets, nent, ncap, secs, cpu = run_frames(name, frames, o['cfg'], get_cal(R), workers, recipe='crowded_polish')
    arrays = {'frames': np.array([a + k for k in DIAG_FRAMES]), 'secs': secs, 'cpu': cpu, 'cfg': np.array(o['cfg']), 'kappa': np.array(o['kappa']), 'smax': np.array(o['smax']), 'off': np.asarray(o['off'])}
    for (kind, kap), ds in dets.items():
        arrays['%s%g|off' % (kind, kap)], arrays['%s%g|flat' % (kind, kap)] = pack(ds)
    log('%s diagnostic (polish, %s kappa %g) on dev frames %s: %.0f s (%.1f s CPU per frame)' % (name, o['cfg'], o['kappa'], list(arrays['frames']), time.time() - t0, cpu.mean()))
    return arrays


# ---------------------------------------------------------------- the PSF file shipped with the set (psf.zip), against the images
PSF_Z_NM = (0.0, 50.0, 100.0, 200.0, 290.0)                         # planes (10 nm voxels) looked at
PSF_RENDER_SHIFTS, PSF_RENDER_PHOTONS = 40, 2000.0


def psf_planes(path):
    """The planes of the PSF stack in a ``psf.zip`` (float32 TIFF 'psf-10.0nm.tif': 10 nm voxels, one plane per 10 nm of z, every plane summing to 1) as (n, ny, nx) float64."""
    from PIL import Image
    with zipfile.ZipFile(str(path)) as z:
        names = [i.filename for i in z.infolist() if i.filename.lower().endswith(('.tif', '.tiff'))]
        if len(names) != 1:
            raise ValueError('expected exactly one TIFF in %s' % path)
        im = Image.open(io.BytesIO(z.read(names[0])))
        out = []
        for k in range(im.n_frames):
            im.seek(k)
            out.append(np.asarray(im).astype(np.float64))
    return np.stack(out)


def render_psf(plane, pixel_nm, voxel_nm=10.0, shifts=PSF_RENDER_SHIFTS, photons=PSF_RENDER_PHOTONS, seed=0):
    """Noise-free camera patches (11 x 11 pixels of ``pixel_nm``) of an emitter with the PSF plane (voxels of ``voxel_nm``) at random sub-pixel positions, the voxels summed into pixels; returns
    the list of patches."""
    rng = np.random.default_rng(seed)
    k = int(round(pixel_nm / voxel_nm))
    n = 11 * k
    ny, nx = plane.shape
    pad = max(ny, nx)
    out = []
    for _ in range(shifts):
        can = np.zeros((n + 2 * pad, n + 2 * pad))
        cx, cy = int(rng.integers(0, k)), int(rng.integers(0, k))
        oy, ox = pad + n // 2 - k // 2 + cy - ny // 2, pad + n // 2 - k // 2 + cx - nx // 2
        can[oy:oy + ny, ox:ox + nx] += plane
        out.append(can[pad:pad + n, pad:pad + n].reshape(11, k, 11, k).sum(axis=(1, 3)) * photons)
    return out


def part_psf(name, log=print):
    """What the PSF file shipped with a set says (second-moment width and the width a pixel-integrated Gaussian fit finds in a noise-free camera rendering at the pixel size of the set, per
    plane), to be compared with the width measured on the images. None when ``psf.zip`` is not downloaded."""
    path = S.data_path(name, 'psf')
    if not os.path.exists(path):
        return None
    st = psf_planes(path)
    ny, nx = st.shape[1:]
    yy, xx = np.mgrid[0:ny, 0:nx]
    px = S.spec(name)['pixel_nm']
    rows = []
    for z in PSF_Z_NM:
        k = int(round(z / 10.0))
        if k >= len(st):
            continue
        pl = st[k]
        tot = pl.sum()
        cx, cy = (pl * xx).sum() / tot, (pl * yy).sum() / tot
        sm = np.sqrt((pl * ((xx - cx) ** 2 + (yy - cy) ** 2)).sum() / tot / 2.0)
        fits = [fit_width(p, 5, 0.0, 0.0)[0][3] for p in render_psf(pl, px)]
        rows.append((z, 10.0 * sm, float(np.median(fits)) * px))
    log('%s psf.zip: %d planes of %dx%d; second-moment sigma / rendered pixel-integrated fit (nm): %s' % (name, len(st), ny, nx, '; '.join('z %g: %.0f / %.0f' % r for r in rows)))
    return dict(rows=np.array(rows), n_planes=np.array(len(st)))


# ---------------------------------------------------------------- the Picasso ladder
def picasso_search(session, gt_dev, dev_idx, cam, log=print):
    """R2: every structural setting of ``exp_smlm_picasso.STRUCT_GRID`` run once at the lowest net gradient, then the coordinate ascent of the post-filters (``FILTER_GRID``) twice (for the
    efficiency and for the Jaccard index; one shared cache). ``exp_smlm_picasso.search_structures`` for the grids and the camera of these sets. Development frames only."""
    out = collections.OrderedDict()
    structs = PIC.struct_list()
    t_all = time.time()
    for i, st in enumerate(structs):
        t0 = time.time()
        cand = session.run(dict(st, min_net_gradient=PIC_NG_LOW, **cam))
        cache = {}
        rE = C.tune_postfilters(cand, gt_dev, dev_idx, FILTER_GRID, objective='efficiency', max_rounds=PIC.MAX_ROUNDS, cache=cache)
        rJ = C.tune_postfilters(cand, gt_dev, dev_idx, FILTER_GRID, objective='jaccard', max_rounds=PIC.MAX_ROUNDS, cache=cache)
        n_bad = int(sum(int((~np.isfinite(d[:, :2])).any(axis=1).sum()) for d in cand if len(d)))
        out[PIC.skey(st)] = dict(struct=dict(st), n_cand=int(sum(len(d) for d in cand)), n_nonfinite=n_bad, n_eval=len(cache), E=PIC._compact(rE), J=PIC._compact(rJ))
        log('  [%d/%d] %s: %d candidates, %d configurations scored, dev E %.1f (rounds %d), dev J %.1f (rounds %d)  [%.0f s, total %.0f s]' % (
            i + 1, len(structs), PIC.skey(st), out[PIC.skey(st)]['n_cand'], len(cache), rE['metrics']['efficiency'], rE['rounds'], 100 * rJ['metrics']['jaccard'], rJ['rounds'],
            time.time() - t0, time.time() - t_all))
    return out


def part_pic(name, py, cache_dir=None, log=print):
    """Rungs R0, R1E, R1J, R2E, R2J of the Picasso ladder on a set: tuned on the development frames only, detections of the development and test frames saved (packed), the convention check
    (Picasso at a fixed net gradient on the first frames of the development range), Picasso's seconds per frame (R0, test frames)."""
    (a, b), (c, d) = S.PROTOCOL[name]['dev'], S.PROTOCOL[name]['test']
    stack = S.data_path(name, 'stack')
    gt_all = S.load_gt(S.data_path(name, 'gt'))
    dev, test = S.load_frames(stack, a, b), S.load_frames(stack, c, d)
    dev_idx, test_idx = np.arange(a, b), np.arange(c, d)
    gt_dev = S.gt_in_frames(gt_all, a, b)
    cam = camera(name)
    arrays = {}

    # ---- the convention check: Picasso at a fixed threshold on the first development frames
    t_ = time.time()
    nconv = min(CONV_FRAMES, len(dev))
    dconv = C.run_picasso(dev[:nconv], dict(cam, min_net_gradient=CONV_NG), python=py)
    conv = S.check_conventions(name, dconv, gt_all, np.arange(a, a + nconv))
    arrays['conv'] = np.array(_jsonable(dict(conv, n_frames=nconv, net_gradient=CONV_NG)))
    log('%s: convention check (%d frames, net gradient > %g): corner TP %d bias (%.1f, %.1f) nm; centre bias (%.1f, %.1f); swapped TP %d; frame shifts %s [%.0f s]' % (
        name, nconv, CONV_NG, conv['corner']['tp'], conv['corner']['bias_x_nm'], conv['corner']['bias_y_nm'], conv['centre']['bias_x_nm'], conv['centre']['bias_y_nm'], conv['swapped']['tp'],
        conv['frames'], time.time() - t_))

    # ---- R1: the net gradient alone
    t_ = time.time()
    best_e, curve = C.tune_picasso(dev, gt_all, GRID_R1, params=cam, frame_indices=dev_idx, objective='efficiency', python=py, verify_frames=min(60, len(dev)))
    best_j = max(curve, key=lambda q: (q['jaccard'], -q['min_net_gradient']))['min_net_gradient']
    log('%s: R1 tuning %.0f s, best efficiency at %g, best Jaccard at %g' % (name, time.time() - t_, best_e, best_j))
    arrays['r1_curve'], arrays['r1_best'] = np.array(_jsonable(curve)), np.array([best_e, best_j])

    # ---- R2: the search over all options
    cache_file = os.path.join(cache_dir, '%s.search.json' % name) if cache_dir else None
    sig = hashlib.sha256(_jsonable(dict(name=name, proto=S.PROTOCOL[name], struct=[(k, list(v)) for k, v in PIC.STRUCT_GRID.items()], filt=[(k, v) for k, v in FILTER_GRID.items()],
                                        rounds=PIC.MAX_ROUNDS, defaults=C.PICASSO_DEFAULTS, cam=cam)).encode()).hexdigest()[:16]
    results = None
    if cache_file and os.path.exists(cache_file):
        blob = json.load(open(cache_file))
        if blob.get('signature') == sig:
            results = collections.OrderedDict(blob['results'])
            log('%s: R2 search results read from %s' % (name, cache_file))
    if results is None:
        t_ = time.time()
        with C.PicassoSession(dev, py) as sess:
            results = picasso_search(sess, gt_dev, dev_idx, cam, log=lambda s: log('%s %s' % (name, s)))
        log('%s: R2 search %.0f s' % (name, time.time() - t_))
        if cache_file:
            os.makedirs(cache_dir, exist_ok=True)
            json.dump(dict(signature=sig, results=results), open(cache_file, 'w'))
    keyE, keyJ = PIC.choose(results)
    chosen = collections.OrderedDict([('R2E', (keyE, 'E')), ('R2J', (keyJ, 'J'))])
    base = {'R0': dict(cam), 'R1E': dict(cam, min_net_gradient=float(best_e)), 'R1J': dict(cam, min_net_gradient=float(best_j))}
    for r, (k, o) in chosen.items():
        base[r] = dict(cam, **PIC.params_of(results[k]['struct'], results[k][o]['filters']))
        base[r]['min_net_gradient'] = float(base[r]['min_net_gradient'])
    arrays['r2'] = np.array(_jsonable(results))
    arrays['r2_keys'] = np.array(_jsonable([keyE, keyJ]))
    arrays['params'] = np.array(_jsonable({r: C._merge_params(p) for r, p in base.items()}))

    # ---- direct runs of every rung on the development and the test frames (and the self-check of the filtered low-threshold run)
    info0 = {}
    with C.PicassoSession(dev, py) as sd, C.PicassoSession(test, py) as st_:
        for r in PIC.RUNGS:
            dd = sd.run(base[r])
            dt = st_.run(base[r], info0 if r == 'R0' else None)
            for split, ds in (('dev', dd), ('test', dt)):
                arrays['%s|%s|off' % (split, r)], arrays['%s|%s|flat' % (split, r)] = pack([x[:, :4] if len(x) else np.zeros((0, 4)) for x in ds])
        for r, (k, o) in chosen.items():
            cand = sd.run(dict(results[k]['struct'], min_net_gradient=PIC_NG_LOW, **cam))
            direct = unpack(arrays['dev|%s|off' % r], arrays['dev|%s|flat' % r])
            for x, y in zip(direct, C.filter_detections(cand, base[r])):
                y = y[:, :4] if len(y) else np.zeros((0, 4))
                if x.shape != y.shape or not np.allclose(x, y.astype(np.float32).astype(np.float64), rtol=1e-6, atol=1e-3, equal_nan=True):
                    raise RuntimeError('%s %s: the filtered low-threshold run differs from the direct run' % (name, r))
    arrays['time_r0_test'] = np.array([info0['run_s'], info0['n_frames']])
    arrays['picasso_version'] = np.array(str(info0.get('picasso_version', 'unknown')))
    arrays['fit_schedule'] = np.array([info0.get('fit_tolerance', np.nan), info0.get('fit_max_iterations', -1)], dtype=float)
    arrays['ng_low'] = np.array(PIC_NG_LOW)
    log('%s: Picasso rungs done (R0 %.4f s per frame on %d test frames)' % (name, info0['run_s'] / info0['n_frames'], info0['n_frames']))
    return arrays


# ==================================================================================================================================
# orchestration of one set
# ==================================================================================================================================
def run_set(name, args, log=print):
    """All parts of a set, each saved when done and skipped when its file exists (unless --force): meta, cal, pic, the development runs of the four configurations, then the test runs of the
    configurations the declared default and the two operating points chosen on the development frames need."""
    pd, force, only = args.parts_dir, args.force, set(args.only) if args.only else None

    def want(tag):
        return only is None or tag in only

    def todo(part, tag):
        return want(tag) and (force or not have_part(name, part, pd))

    if todo('meta', 'meta'):
        save_part(name, 'meta', part_meta(name, log), pd)
    if todo('cal', 'cal'):
        save_part(name, 'cal', part_cal(name, log), pd)
    if todo('psf', 'psf'):
        pp = part_psf(name, log)
        if pp is not None:
            save_part(name, 'psf', pp, pd)
    if todo('pic', 'pic'):
        save_part(name, 'pic', part_pic(name, args.picasso_python, args.cache_dir, log), pd)
    R = load_set_parts(name, pd)
    for cfg in CONFIGS:
        if todo('isdev.' + cfg, 'isdev'):
            save_part(name, 'isdev.' + cfg, part_is(name, 'dev', cfg, R, args.workers, log), pd)
    if RECIPE_OF[name] == 'crowded' and todo('diag', 'diag') and all(have_part(name, 'isdev.' + c, pd) for c in CONFIGS):
        save_part(name, 'diag', part_diag(name, load_set_parts(name, pd), args.workers, log), pd)
    if want('istest'):
        R = load_set_parts(name, pd)
        if all(have_part(name, 'isdev.' + c, pd) for c in CONFIGS):
            rows, ops = select_ops(R, name)
            cfgs = sorted({DEFAULT['cfg']} | {o['cfg'] for o in ops.values()})
            log('%s: operating points chosen on the development frames: %s; test configurations %s' % (
                name, '; '.join('%s: %s kappa %g width <= %.3f px offset (%.1f, %.1f) dev E %.1f J %.1f' % (w, o['cfg'], o['kappa'], o['smax'], o['off'][0], o['off'][1], o['m'][4], 100 * o['m'][0])
                                for w, o in ops.items()), cfgs))
            for cfg in (cfgs[::-1] if getattr(args, 'reverse', False) else cfgs):
                if force or not have_part(name, 'istest.' + cfg, pd):
                    save_part(name, 'istest.' + cfg, part_is(name, 'test', cfg, R, args.workers, log), pd)
        else:
            log('%s: the development runs of all configurations are needed before the test runs' % name)


# ==================================================================================================================================
# the report (from the raw alone)
# ==================================================================================================================================
SHORT = {'Tubulins_I': 'Tubulins_I', 'Tubulins_II': 'Tubulins_II', 'Bundled_Tubes_High_Density': 'Bundled HD', 'Bundled_Tubes_Long_Sequence': 'Bundled LS'}
VARIANTS = PIC.VARIANTS
LADDER_VS = ISS.LADDER_VS
fmt_ci, table = ISS.fmt_ci, ISS.table


def _json(R, key):
    return json.loads(str(R[key]))


def _off(dets_dev, gt_dev, idx, name):
    r = S.score(dets_dev, gt_dev, idx, name)
    return np.nan_to_num(np.array([r['bias_x_nm'], r['bias_y_nm']]))


def evaluate_set(R, name, boot):
    """Everything scored on the test frames of a set, from the raw: dict(idx, res (primary), res2 (challenge), ci, ci2, ops, rows, off, dev (dev results of the Picasso rungs), meths (names))."""
    c, d = (int(v) for v in R['meta|frame_ranges'][2:4])
    a, b = (int(v) for v in R['meta|frame_ranges'][0:2])
    idx, idx_dev = np.arange(c, d), np.arange(a, b)
    gt, gt_dev, pmin = R['meta|gt_test'], R['meta|gt_dev'], float(R['meta|pmin'])
    kind = RULE_KIND[RECIPE_OF[name]]
    dets, off, dev_res = collections.OrderedDict(), {}, {}
    for r in PIC.RUNGS:
        dd = unpack(R['pic|dev|%s|off' % r], R['pic|dev|%s|flat' % r])
        dt = unpack(R['pic|test|%s|off' % r], R['pic|test|%s|flat' % r])
        off[r] = _off(dd, gt_dev, idx_dev, name)
        dev_res[r] = S.score(dd, gt_dev, idx_dev, name)
        dets[r], dets[r + '+c'] = dt, C.shift_detections(dt, *off[r])
    rows, ops = select_ops(R, name)
    dets['ISS-0'] = is_dets(R, 'test', DEFAULT['cfg'], kind, DEFAULT['kappa'])
    for w in OPS:
        dets['ISS-' + w] = op_test_dets(R, name, ops[w])
    res, res2, ci, ci2 = {}, {}, {}, {}
    for k, ds in dets.items():
        res[k] = S.score(ds, gt, idx, name)
        res2[k] = S.score(ds, gt, idx, name, border_nm=CHALLENGE['border_nm'], photon_min=pmin)
        ci[k] = C.bootstrap_ci(res[k], B=boot, seed=0, block=BLOCK, keys=C.BOOT_KEYS_EXT)
        ci2[k] = C.bootstrap_ci(res2[k], B=boot, seed=0, block=BLOCK)
    return dict(idx=idx, res=res, res2=res2, ci=ci, ci2=ci2, ops=ops, rows=rows, off=off, dev=dev_res, dets=dets)


NN_BINS = ((0.0, 150.0), (150.0, 250.0), (250.0, 400.0), (400.0, 1e9))


def nn_distance(gt, frame_indices):
    """For every row of ``gt`` the distance (nm) to the nearest OTHER row of the same frame (inf if it is alone), for the frames in ``frame_indices``; rows of other frames get inf."""
    out = np.full(len(gt), np.inf)
    for f in np.unique(np.asarray(frame_indices)):
        m = np.flatnonzero(gt['frame'] == f)
        if len(m) > 1:
            xy = np.c_[gt['x_nm'][m], gt['y_nm'][m]]
            d = np.hypot(xy[:, None, 0] - xy[None, :, 0], xy[:, None, 1] - xy[None, :, 1])
            np.fill_diagonal(d, np.inf)
            out[m] = d.min(axis=1)
    return out


def diff_ci(ev, A, B, boot):
    d = C.bootstrap_diff_ci(ev['res'][A]['per_frame'], ev['res'][B]['per_frame'], B=boot, seed=0, block=BLOCK)
    d2 = C.bootstrap_diff_ci(ev['res2'][A]['per_frame'], ev['res2'][B]['per_frame'], B=boot, seed=0, block=BLOCK)
    return d, d2


def method_label(key, ev, name):
    kind = RULE_KIND[RECIPE_OF[name]]
    if key == 'ISS-0':
        return 'ISS, a priori default (%s, kappa %g, no width bound, no offset correction)' % (DEFAULT['cfg'], DEFAULT['kappa'])
    if key in ('ISS-E', 'ISS-J'):
        o = ev['ops'][key[-1]]
        return 'ISS, %s-tuned (%s, kappa %g, width <= %.2f px, offset (%.1f, %.1f) removed)' % (key[-1], o['cfg'], o['kappa'], o['smax'], o['off'][0], o['off'][1])
    return 'Picasso %s' % key


def set_section(R, name, ev, boot):
    """The markdown lines of one set. Returns (lines, summary tuples, tuning tuple)."""
    L = []
    P = L.append
    sp = S.spec(name)
    px = sp['pixel_nm']
    a, b, c, d = (int(v) for v in R['meta|frame_ranges'])
    n_dev, n_test = b - a, d - c
    cal = get_cal(R)
    src = str(R['cal|src'])
    rec = RECIPE_OF[name]
    P('## %s\n' % name)
    P('%d x %d px of %g nm, %d frames, %d ground-truth rows in the whole set (%.2f per frame on average, %d to %d). Development frames %d-%d (%d frames, %d ground-truth emitters), test frames %d-%d '
      '(%d frames, %d emitters). Photon threshold of the challenge setting (25 %% quantile of the whole set\'s ground truth): %.0f photons. Persistence: %d runs of one fluorophore in consecutive '
      'frames, mean run length %.2f frames (every activation lasts one frame).\n' % (
          sp['size'], sp['size'], px, int(R['meta|n_frames']), int(R['meta|n_gt']), int(R['meta|n_gt']) / int(R['meta|n_frames']), int(R['meta|per_frame'].min()), int(R['meta|per_frame'].max()),
          a, b - 1, n_dev, len(R['meta|gt_dev']), c, d - 1, n_test, len(R['meta|gt_test']), float(R['meta|pmin']), int(R['meta|runs'][0]), float(R['meta|runs'][1])))
    # ---- calibration
    sig = R['cal|sigma_sample']
    af, bf, span = (float(v) for v in R['cal|nl_free'])
    rows_nl = R['cal|nl_rows']
    pick = [0, len(rows_nl) // 4, len(rows_nl) // 2, 3 * len(rows_nl) // 4, len(rows_nl) - 1]
    P('### Calibration and the ISS recipe\n')
    P('Calibration (development frames %d-%d of %s%s; no labels except for sigma0): background surface mean %.2f photons (range %.1f to %.1f); noise law var = a x + b0 with a = %.2f (the page\'s EM '
      'factor) and **b0 = %.1f photon^2** (the page\'s read-out sd %g gives %g); robust high-frequency variance against the smooth level (level: variance, model; the bins above the 90th percentile '
      'of the level hold the emitters of a sparse frame and are not used for b0): %s; %s. **sigma0 = %.3f px (%.0f nm)**: median over %d isolated bright emitters of the '
      'development frames of %s of the width of a pixel-integrated Gaussian plus a plane (IQR %.3f to %.3f; 5 to 95 %% %.3f to %.3f); the page\'s description (FWHM %.2f nm at the focal plane, '
      'growing exponentially to twice that at %g nm) gives %.3f px at the median |z| of the sample (%.0f nm) and %.3f px in focus; the fitted photons are %.2f (median; IQR %.2f to %.2f) of the '
      'ground-truth photons; the fit minus truth offset is (%.1f, %.1f) nm (median). The widths slide in [%.3f, %.3f] px, r_ex = %.3f px, width bounds of the post-filter %s px.\n' % (
          *S.PROTOCOL[src]['dev'][:1], S.PROTOCOL[src]['dev'][1] - 1, src, '' if src == name else '; the dense set uses the calibration of its sparse set', R['cal|bg'].mean(), R['cal|bg'].min(),
          R['cal|bg'].max(), cal['a'], cal['b0'], sp['read_sd'], sp['read_sd'] ** 2, '; '.join('%.0f: %.0f, %.0f' % (rows_nl[i, 0], rows_nl[i, 1], cal['a'] * max(rows_nl[i, 0], 0) + cal['b0']) for i in pick),
          ('the free slope through all bins is %.2f (span of the levels %.0f photons)' % (af, span)) if span >= 100.0 else ('the slope is not determinable from the bins here (the levels span only %.0f photons)' % span),
          cal['sigma0'], cal['sigma0'] * px, len(sig), src, *np.percentile(sig[:, 0], [25, 75]), *np.percentile(sig[:, 0], [5, 95]), S.PSF_FWHM_NM, S.spec(src)['z_double_nm'],
          float(np.median(S.sigma_page_nm(sig[:, 5], src)) / S.spec(src)['pixel_nm']), float(np.median(sig[:, 5])), float(S.sigma_page_nm(0.0, src) / S.spec(src)['pixel_nm']),
          float(np.median(sig[:, 3] / sig[:, 4])), *np.percentile(sig[:, 3] / sig[:, 4], [25, 75]), float(np.median(sig[:, 1])), float(np.median(sig[:, 2])), WIDTH_REL[0] * cal['sigma0'],
          WIDTH_REL[1] * cal['sigma0'], REX_REL * cal['sigma0'], ', '.join('%.3f' % (r_ * cal['sigma0']) for r_ in SMAX_REL)))
    if 'psf|rows' in R:
        P('The PSF file shipped with the set (`psf.zip`, a stack of planes of 10 nm voxels, each summing to 1; the page\'s PSF description is its z = 0 plane, FWHM %.2f nm): second-moment sigma and the width '
          'that the same pixel-integrated Gaussian fit finds in a noise-free camera rendering of it at %g nm pixels (nm), by plane: %s. The widths measured on the images above (sigma0 = %.0f nm, '
          'median over the emitters at |z| = %.0f nm) are smaller than any rendering of this PSF at the emitters\' depths (%.0f nm at z = 50 nm), the fitted photons are below the ground-truth photons '
          'while the PSF planes are normalized: the images do not follow the shipped PSF exactly (the cause was not investigated); the experiment uses the measured width.\n' % (
              S.PSF_FWHM_NM, px, '; '.join('z %g nm: %.0f / %.0f' % tuple(r) for r in R['psf|rows']), cal['sigma0'] * px, float(np.median(sig[:, 5])),
              float(np.interp(50.0, R['psf|rows'][:, 0], R['psf|rows'][:, 2]))))
    bt = R['istest.%s|bg_test' % DEFAULT['cfg']].astype(float) if (src == name and 'istest.%s|bg_test' % DEFAULT['cfg'] in R) else None      # (a dense set has no own calibration to compare)
    if bt is not None:
        dbg = bt - R['cal|bg']
        P('Stability of the background (label-free, measured after the test run, not used anywhere): the smoothed temporal median of the test frames minus that of the calibration frames has rms %.2f '
          'photons (|max| %.1f; mean %.2f), against the pixel noise sd of %.1f photons at the mean background.\n' % (float(np.sqrt(np.mean(dbg ** 2))), float(np.abs(dbg).max()), float(dbg.mean()),
                                                                                                                  float(np.sqrt(cal['a'] * max(R['cal|bg'].mean(), 0) + cal['b0']))))
    if rec == 'sparse':
        P('Recipe: %s; entries per frame %.1f on average (test, default configuration), tiles at the cap of %d entries in the test run: %d.\n' % (
            'tiles of 64 x 64 px cores with 8 px margins, the path stopped at the first entrant below kappa = %g' % K_MIN if sp['size'] > 64 else 'the whole 64 x 64 frame is one tile, the path stopped at the first entrant below kappa = %g' % K_MIN,
            float(R['istest.%s|entries' % DEFAULT['cfg']].mean()), RECIPES[rec]['cap'], int(R['istest.%s|capped' % DEFAULT['cfg']].sum())))
    else:
        P('Recipe (crowded): refits of the kernels within %g px of an entrant only, window 10 px beyond, no final polish, at most %d entries and no early stop, the final block cleaned by the units\' own '
          'significance (rule "clean"); entries per frame %.0f on average (the cap binds in every frame by design: the path runs on into the noise).\n' % (
              RECIPES[rec]['fit_radius'], RECIPES[rec]['cap'], float(R['istest.%s|entries' % DEFAULT['cfg']].mean())))

    # ---- Picasso on the development frames
    params = _json(R, 'pic|params')
    results = _json(R, 'pic|r2')
    keyE, keyJ = _json(R, 'pic|r2_keys')
    n_struct = len(results)
    n_eval = sum(r['n_eval'] for r in results.values())
    n_cand = sum(r['n_cand'] for r in results.values())
    n_conv = sum(int(r[o]['converged']) for r in results.values() for o in ('E', 'J'))
    P('### Picasso ladder: tuning on the development frames\n')
    P('Picasso %s, camera passed: baseline 100, gain 1, sensitivity 1, QE 1, pixel %g nm. **Configurations tried on the development frames: %d structural settings run, %d distinct (structure, '
      'post-filter) configurations scored** (R1: %d net-gradient values, all contained in them); %d candidate localizations fitted in total; the coordinate ascent converged within %d rounds in %d of %d '
      'runs. Net-gradient grid of R1: %s; R2 post-filter grids: net gradient %d values from %g to %g, minimum photons %s, minimum width %s nm, maximum width %s nm.\n' % (
          str(R['pic|picasso_version']), px, n_struct, n_eval, len(GRID_R1), n_cand, PIC.MAX_ROUNDS, n_conv, 2 * n_struct, ', '.join('%g' % g for g in GRID_R1), len(FILTER_GRID['min_net_gradient']),
          min(FILTER_GRID['min_net_gradient']), max(FILTER_GRID['min_net_gradient']), '/'.join('off' if v is None else '%g' % v for v in FILTER_GRID['min_photons']),
          '/'.join('off' if v is None else '%g' % v for v in FILTER_GRID['min_sigma_nm']), '/'.join('off' if v is None else '%g' % v for v in FILTER_GRID['max_sigma_nm'])))
    curve = _json(R, 'pic|r1_curve')
    be, bj = (float(v) for v in R['pic|r1_best'])
    P('R1 tuning curve of the minimum net gradient (development frames; bold: best E and best J):\n')
    P(table(['min net gradient', 'J %', 'precision %', 'recall %', 'RMSE nm', 'efficiency', 'detections / frame'],
            [[('**%g**' if q['min_net_gradient'] in (be, bj) else '%g') % q['min_net_gradient'], '%.1f' % (100 * q['jaccard']), '%.1f' % (100 * q['precision']), '%.1f' % (100 * q['recall']),
              '%.1f' % q['rmse_nm'], '%.1f' % q['efficiency'], '%.2f' % (q['n_det'] / n_dev)] for q in curve]))
    P('\nChosen configurations and their scores ON THE DEVELOPMENT FRAMES (all ground-truth emitters); offset = mean signed error (detection - truth) of the matched pairs, the correction of the "+c" rows:\n')
    P(table(['rung', 'configuration', 'dev J %', 'dev precision %', 'dev recall %', 'dev RMSE nm', 'dev E', 'offset (dx, dy) nm', 'dev TP'],
            [[r, PIC.describe(params[r]), '%.1f' % (100 * ev['dev'][r]['jaccard']), '%.1f' % (100 * ev['dev'][r]['precision']), '%.1f' % (100 * ev['dev'][r]['recall']),
              '%.1f' % ev['dev'][r]['rmse_nm'], '%.1f' % ev['dev'][r]['efficiency'], '(%.1f, %.1f)' % tuple(ev['off'][r]), ev['dev'][r]['tp']] for r in PIC.RUNGS]))
    P('')

    # ---- ISS on the development frames
    rows, ops = ev['rows'], ev['ops']
    kind = RULE_KIND[rec]
    P('### ISS: operating points on the development frames\n')
    P('%d operating points (%d configurations x %d significance thresholds x %d width bounds; rule "%s"), each with the position offset estimated on the same frames and removed. Bold: the point of the '
      'largest efficiency (E) and the one of the largest Jaccard index (J), the two frozen for the test frames.\n' % (len(rows), len(cfg_list(R)), len(KAPPAS), len(SMAX_REL), kind))
    smax_nf = SMAX_REL[-1] * cal['sigma0']                                   # the largest width bound = the bound of the sliding widths: no filter
    shown = [r for r in rows if abs(r['smax'] - smax_nf) < 1e-9 or r is ops['E'] or r is ops['J']]
    P(table(['config', 'kappa', 'width <= px', 'offset (dx, dy) nm', 'J %', 'precision %', 'recall %', 'RMSE nm', 'efficiency', 'J % (challenge)', 'efficiency (challenge)'],
            [[('**%s**' if (r is ops['E'] or r is ops['J']) else '%s') % r['cfg'], '%g' % r['kappa'], '%.3f' % r['smax'], '(%.1f, %.1f)' % tuple(r['off']), '%.1f' % (100 * r['m'][0]), '%.1f' % (100 * r['m'][1]),
              '%.1f' % (100 * r['m'][2]), '%.1f' % r['m'][3], '%.1f' % r['m'][4], '%.1f' % (100 * r['m'][5]), '%.1f' % r['m'][6]] for r in shown]))
    P('')

    # ---- test frames
    keys = list(VARIANTS) + ['ISS-0', 'ISS-E', 'ISS-J']
    P('### Test frames %d-%d\n' % (c, d - 1))
    P('All ground-truth emitters (primary setting); %d frames, %.3f true emitters per frame. Rows R0, R1E / R1J (net gradient tuned), R2E / R2J (every option tuned) of the Picasso ladder with and without the '
      'offset correction "+c"; ISS rows as defined above. 95 %% intervals of a block bootstrap (blocks of %d frames, %d resamples, seed 0; the test range holds %d blocks%s).\n' % (
          n_test, ev['res']['R0']['n_gt'] / n_test, BLOCK, boot, n_test // BLOCK, ': few for a percentile interval, treat the intervals as rough' if n_test // BLOCK < 10 else ''))
    P(table(['method', 'J %', 'precision %', 'recall %', 'RMSE nm', 'efficiency', 'detected / frame'],
            [[method_label(k, ev, name) if k.startswith('ISS') else k, fmt_ci(ev['res'][k]['jaccard'], ev['ci'][k]['jaccard']), fmt_ci(ev['res'][k]['precision'], ev['ci'][k]['precision']),
              fmt_ci(ev['res'][k]['recall'], ev['ci'][k]['recall']), fmt_ci(ev['res'][k]['rmse_nm'], ev['ci'][k]['rmse_nm'], 1.0), fmt_ci(ev['res'][k]['efficiency'], ev['ci'][k]['efficiency'], 1.0),
              '%.3f' % (ev['res'][k]['n_det'] / n_test)] for k in keys]))
    P('\nTP / FP / FN: %s.\n' % '; '.join('%s %d / %d / %d' % (k, ev['res'][k]['tp'], ev['res'][k]['fp'], ev['res'][k]['fn']) for k in ['R0', 'R2E', 'R2J', 'ISS-0', 'ISS-E', 'ISS-J']))
    P('Challenge assessment setting (emitters within %g nm of the border of the field and the dimmest quarter of the whole set\'s activations dropped after the pairing; the rules are those of '
      '`smlm_common.score`, reconstructed from the 2016 challenge description, with the field of the set):\n' % CHALLENGE['border_nm'])
    P(table(['method', 'J %', 'RMSE nm', 'efficiency'],
            [[method_label(k, ev, name) if k.startswith('ISS') else k, fmt_ci(ev['res2'][k]['jaccard'], ev['ci2'][k]['jaccard']), fmt_ci(ev['res2'][k]['rmse_nm'], ev['ci2'][k]['rmse_nm'], 1.0),
              fmt_ci(ev['res2'][k]['efficiency'], ev['ci2'][k]['efficiency'], 1.0)] for k in keys]))
    P('\nPaired differences (ISS minus Picasso, same frames, same resampled blocks), 95 % intervals; the offset-corrected rung is the like-for-like comparison:\n')
    rows_d, summary = [], []
    pairs = [('ISS-0', 'R0', 'a priori default ISS vs Picasso defaults (R0)')]
    for w in OPS:
        pairs += [('ISS-' + w, pr, '%s-tuned ISS vs Picasso %s' % (w, pr)) for pr in LADDER_VS[w]]
    for A, B, label in pairs:
        dd, dd2 = diff_ci(ev, A, B, boot)
        ra, rb, ra2, rb2 = ev['res'][A], ev['res'][B], ev['res2'][A], ev['res2'][B]
        rows_d.append([label, fmt_ci(ra['jaccard'] - rb['jaccard'], dd['jaccard']), fmt_ci(ra['rmse_nm'] - rb['rmse_nm'], dd['rmse_nm'], 1.0), fmt_ci(ra['efficiency'] - rb['efficiency'], dd['efficiency'], 1.0),
                       fmt_ci(ra2['efficiency'] - rb2['efficiency'], dd2['efficiency'], 1.0)])
        if B.endswith('+c'):
            summary.append((name, A[-1], ra['jaccard'] - rb['jaccard'], dd['jaccard'], ra['rmse_nm'] - rb['rmse_nm'], dd['rmse_nm'], ra['efficiency'] - rb['efficiency'], dd['efficiency'],
                            ra2['efficiency'] - rb2['efficiency'], dd2['efficiency']))
    P(table(['comparison', 'dJ % (all)', 'dRMSE nm (all)', 'dE (all)', 'dE (challenge)'], rows_d))
    # ---- the RMSE on the emitters both methods of a pair match (free of the survivor bias of the stricter method)
    P('\nRMSE on the ground-truth emitters that BOTH methods of a pair match (a method that loses its hard emitters is not flattered; all emitters, primary setting; paired bootstrap difference):\n')
    crow = []
    me_cache = {}
    for A, B in [('ISS-0', 'R0'), ('ISS-E', 'R2E+c'), ('ISS-J', 'R2J+c')]:
        for k in (A, B):
            if k not in me_cache:
                me_cache[k] = C.matched_errors(ev['dets'][k], R['meta|gt_test'], ev['idx'])
        common = me_cache[A]['matched'] & me_cache[B]['matched']
        ca, cb = PIC.common_columns(me_cache[A], common, n_test), PIC.common_columns(me_cache[B], common, n_test)
        if common.sum() == 0:
            crow.append([A + ' vs ' + B, 0, 'nan', 'nan', 'nan'])
            continue
        ra, rb = float(PIC.rmse_stat({k: v.sum() for k, v in ca.items()})), float(PIC.rmse_stat({k: v.sum() for k, v in cb.items()}))
        dci = C.bootstrap_pooled_diff_ci(ca, cb, PIC.rmse_stat, B=boot, seed=0, block=BLOCK)
        crow.append([A + ' vs ' + B, int(common.sum()), '%.1f' % ra, '%.1f' % rb, '%+.1f [%+.1f, %+.1f]' % (ra - rb, *dci)])
    P(table(['pair', 'common true positives', 'RMSE nm, first', 'RMSE nm, second', 'difference nm [95 % interval]'], crow))
    # ---- recall and RMSE by the distance to the nearest other emitter of the same frame (where the methods differ)
    gtt = R['meta|gt_test']
    nn = nn_distance(gtt, ev['idx'])
    P('\nRecall and RMSE of the matched pairs by the distance of the ground-truth emitter to the nearest OTHER ground-truth emitter of its frame (test frames, all emitters; a post-hoc description of where '
      'the methods differ, not a result of the protocol; RMSE in nm over the matched pairs of the bin):\n')
    brow = []
    for lo, hi in NN_BINS:
        sel = (nn >= lo) & (nn < hi)
        cells = ['%d-%s nm' % (lo, 'inf' if hi > 1e8 else '%d' % hi), int(sel.sum())]
        for k in ('R2E+c', 'ISS-0', 'ISS-E'):
            if k not in me_cache:
                me_cache[k] = C.matched_errors(ev['dets'][k], gtt, ev['idx'])
            me = me_cache[k]
            nnr = nn[me['gt_row']]
            ss = (nnr >= lo) & (nnr < hi)
            mt = ss & me['matched']
            cells += ['%.2f' % (mt.sum() / ss.sum()) if ss.sum() else 'nan', '%.1f' % np.sqrt(np.mean(me['dx'][mt] ** 2 + me['dy'][mt] ** 2)) if mt.sum() else 'nan']
        brow.append(cells)
    P(table(['nearest other emitter', 'emitters', 'recall R2E+c', 'RMSE R2E+c', 'recall ISS default', 'RMSE ISS default', 'recall ISS E-tuned', 'RMSE ISS E-tuned'], brow))
    tuning = (name, ev['res']['R0']['efficiency'], ev['res']['R2E+c']['efficiency'], ev['res']['ISS-0']['efficiency'], ev['res']['ISS-E']['efficiency'], ev['res']['R0']['jaccard'],
              ev['res']['R2J+c']['jaccard'], ev['res']['ISS-0']['jaccard'], ev['res']['ISS-J']['jaccard'])
    # ---- seconds per frame
    oE = ops['E']
    sE = R['istest.%s|secs' % oE['cfg']]
    cE = R['istest.%s|cpu' % oE['cfg']]
    t0r = R['pic|time_r0_test']
    P('\nSeconds per frame (one core of a shared machine; wall clock, and process CPU time for the ISS): E-tuned ISS (%s) %.2f wall (median %.2f), %.2f CPU; Picasso with its defaults %.4f wall (%d frames).\n' % (
        oE['cfg'], sE.mean(), np.median(sE), cE.mean(), t0r[0] / t0r[1], int(t0r[1])))
    # ---- the early stop in the dense set
    if rec == 'crowded':
        P('#### The protocol\'s early-stop rule against the crowded recipe (dense set)\n')
        P('Same path, same default configuration (%s); "snr" = the first-entrant-below-kappa rule of the 2016 pipeline (the block before it, no final polish), "clean" = the crowded recipe (final block, units '
          'of own significance >= kappa). Development and test frames; the detections per frame show how early the "snr" rule stops (true emitters per frame: %.1f dev, %.1f test).\n' % (
              DEFAULT['cfg'], len(R['meta|gt_dev']) / n_dev, len(R['meta|gt_test']) / n_test))
        drows = []
        for split, gts, fr in (('dev', R['meta|gt_dev'], np.arange(a, b)), ('test', R['meta|gt_test'], ev['idx'])):
            for kind2 in ('snr', 'clean'):
                for kap in (DEFAULT['kappa'], 8.0):
                    ds = is_dets(R, split, DEFAULT['cfg'], kind2, kap)
                    r = S.score(ds, gts, fr, name)
                    drows.append([split, kind2, '%g' % kap, '%.1f' % (100 * r['jaccard']), '%.1f' % (100 * r['precision']), '%.1f' % (100 * r['recall']), '%.1f' % r['rmse_nm'], '%.1f' % r['efficiency'],
                                  '%.1f' % (r['n_det'] / len(fr))])
        P(table(['frames', 'rule', 'kappa', 'J %', 'precision %', 'recall %', 'RMSE nm', 'efficiency', 'detected / frame'], drows))
        P('')
        if 'diag|frames' in R:
            dfr = R['diag|frames']
            gd = R['meta|gt_dev']
            gsel = gd[np.isin(gd['frame'], dfr)]
            nn_d = nn_distance(gsel, dfr)
            cfg_d, kap_d, smax_d, off_d = str(R['diag|cfg']), float(R['diag|kappa']), float(R['diag|smax']), R['diag|off']
            plain = [x for i_, x in enumerate(is_dets(R, 'dev', cfg_d, 'clean', kap_d)) if (a + i_) in set(dfr.tolist())]
            polished = unpack(R['diag|clean%g|off' % kap_d], R['diag|clean%g|flat' % kap_d])
            P('#### Diagnostic: does the missing final polish of the crowded recipe cost accuracy? (development frames %s only, E-tuned operating point: %s, kappa %g, width <= %.2f px, offset removed; not a rung of the '
              'protocol, decided after the test runs and not used to change them)\n' % (', '.join(str(int(v)) for v in dfr), cfg_d, kap_d, smax_d))
            drows = []
            for lab, ds in (('no polish (the crowded recipe as run)', plain), ('polish of the cleaned block (1 round, at most 30 evaluations)', polished)):
                ds = ISS.shift(post(ds, smax_d), off_d)
                r = S.score(ds, gsel, dfr, name)
                me = C.matched_errors(ds, gsel, dfr)
                iso = me['matched'] & (nn_d[me['gt_row']] > 300.0)
                drows.append([lab, '%.1f' % (100 * r['jaccard']), '%.1f' % (100 * r['precision']), '%.1f' % (100 * r['recall']), '%.1f' % r['rmse_nm'], '%.1f' % r['efficiency'],
                              '%d: %.1f' % (iso.sum(), np.sqrt(np.mean(me['dx'][iso] ** 2 + me['dy'][iso] ** 2)) if iso.sum() else float('nan'))])
            P(table(['variant', 'J %', 'precision %', 'recall %', 'RMSE nm', 'efficiency', 'matched emitters with no other emitter within 300 nm: n, RMSE nm'], drows))
            P('\nThe polish refits all kept units together against the whole frame; the diagnostic run polishes the block of each of the 7 thresholds (%.0f s CPU per frame, against %.0f s for the run without '
              'polish). The sample is small (%d frames, %d emitters in the last column), so this is an indication, not a measurement.\n' % (float(R['diag|cpu'].mean()), float(R['istest.%s|cpu' % cfg_d].mean()) if 'istest.%s|cpu' % cfg_d in R else float('nan'), len(dfr), int(drows[0][-1].split(':')[0])))
    return L, summary, tuning


def report(R, boot=BOOT, sets=SETS):
    """The markdown report from the raw alone (``load_raw``)."""
    sets = [s for s in sets if s in R and 'istest.%s|frames' % DEFAULT['cfg'] in R[s]]
    L = ['# The sliding inverse scale space and the Picasso ladder on the 2013 challenge\'s simulated 2D training sets', '']
    ver = _json(R[sets[0]], 'meta|versions')
    shas = ['%s: stack %s, ground truth %s' % (s, _json(R[s], 'meta|sha')['stack']['sha256'][:12], _json(R[s], 'meta|sha')['gt']['sha256'][:12]) for s in sets]
    bad = [s for s in sets for k in ('stack', 'gt') if _json(R[s], 'meta|sha')[k]['sha256'] != _json(R[s], 'meta|sha')[k]['recorded']]
    L += ['Generated by `experiments/exp_smlm2013.py` from `smlm2013_raw.npz`; every number below is printed by the script. Protocol, calibration, the Picasso grids and the declared ISS recipes: the '
          'script\'s docstring; data, conventions and frame ranges: `experiments/smlm2013.py`; data (not in the repository): `data/smlm/fetch_2013.py`, sha256 of the files used (first 12 digits; %s the digests recorded in `data/smlm/checksums_2013.json`): %s. '
          'Frame indices are 0-based stack indices, positions in nm (corner convention), metrics pooled over the frames with 95 %% block-bootstrap intervals (blocks of %d consecutive frames, %d resamples, seed 0); '
          'differences are paired (the same resampled blocks for both methods). Python %s, numpy %s, Picasso %s.' % ('all equal to' if not bad else 'NOT EQUAL for ' + ', '.join(bad) + ' to', '; '.join(shas), BLOCK, boot, ver['python'], ver['numpy'],
                                                                                                                 str(R[sets[0]]['pic|picasso_version'])), '']
    L += ['**What these sets are.** The four training sets of the 2013 challenge with public ground truth (the hub\'s page: "generated with a 2D theoretical PSF"): a realistic tubulin structure at low (Tubulins_I) '
          'and high (Tubulins_II) read-out noise, 256 x 256 px of 150 nm, and bundles of 8 tubes at high (Bundled HD) and low (Bundled LS) emitter density, 64 x 64 px of 100 nm. They differ from the '
          'four 2016 sets of `iss.md` in the structure, the pixel size, the camera (QE 1, 1 DN per electron), the PSF (a Gaussian of about 100 nm without halo, against 140 nm with a long halo) and the background '
          '(structured, time-varying autofluorescence), and in the ground truth (one-frame activations). The pages state "168 frames" for the high-density bundles and "2400 frames" for the tubulins; the files '
          'hold 361 and 2401 (`smlm2013.py`). The ground truth is `fluorophores/activation.csv` of each `fluorophores.zip` (the `oracle.zip` files hold only projection images); '
          'it was checked equal to the per-frame files. The ISS recipe (tiles, crowded fit) and the Picasso grids are adapted to the sets as declared in the script\'s docstring; **the long sequences are used in part** '
          '(frame ranges below), the dense set in 30 development and 60 test frames.', '']
    # ---- data table
    rows = []
    for s in sets:
        r_, sp = R[s], S.spec(s)
        pf = r_['meta|per_frame']
        pq = r_['meta|photon_q']
        rows.append([s, '%d x %d px of %g nm' % (sp['size'], sp['size'], sp['pixel_nm']), int(r_['meta|n_frames']), int(r_['meta|n_gt']), '%.2f (%d to %d)' % (pf.mean(), pf.min(), pf.max()), '%.0f / %.0f' % (pq[3], pq[2]),
                     '%.0f to %.0f' % tuple(r_['meta|z_range']), '%d-%d / %d-%d' % (r_['meta|frame_ranges'][0], r_['meta|frame_ranges'][1] - 1, r_['meta|frame_ranges'][2], r_['meta|frame_ranges'][3] - 1),
                     '%d / %d' % (len(r_['meta|gt_dev']), len(r_['meta|gt_test'])), RECIPE_OF[s]])
    L += ['## The sets', '', table(['set', 'frame', 'frames', 'ground-truth rows', 'emitters per frame', 'photons median / 25 % quantile', 'z (nm)', 'dev / test frames', 'dev / test emitters', 'ISS recipe'], rows), '']
    # ---- conventions
    L += ['## Conventions checked against the images', '',
          'Picasso (R0 settings with the camera of the set, minimum net gradient %g, first %d development frames; positions in the corner convention (p + 0.5) x pixel) matched to the ground truth, '
          'and the same detections under the alternatives that would be wrong if the files followed another convention. TP = matched pairs (250 nm). The corner convention without any shift is the files\' '
          'convention: the mean signed error is within a few nm of zero, the "centre" convention (pixel centre at j x pixel) would give -p/2, an exchange of the axes or a flipped axis destroys the matching, and '
          'a shift of the ground-truth frame by one or two destroys it in the sparse sets. In the dense set (one emitter per 18 pixels, a static structure) a detection finds some ground-truth emitter within '
          '250 nm after an exchange of the axes, a flip or a frame shift as well, so the matched counts do not discriminate there: the evidence is the mean error (corner: within a nm or two of zero; centre: tens of nm '
          'negative, less than -p/2 because the pairing then finds a nearer neighbour) and the excess table below.' % (CONV_NG, CONV_FRAMES), '']
    crow = []
    for s in sets:
        cv = _json(R[s], 'pic|conv')
        crow.append([s, '%d' % cv['corner']['tp'], '(%.1f, %.1f)' % (cv['corner']['bias_x_nm'], cv['corner']['bias_y_nm']), '(%.1f, %.1f)' % (cv['centre']['bias_x_nm'], cv['centre']['bias_y_nm']),
                     '%d / %d / %d' % (cv['swapped']['tp'], cv['flip_x']['tp'], cv['flip_y']['tp']), ' / '.join('%d' % cv['frames'][k] for k in ('-2', '-1', '0', '1', '2'))])
    L += [table(['set', 'TP (corner)', 'mean error (dx, dy) nm, corner', 'mean error, centre', 'TP with x<->y / x flipped / y flipped', 'TP with the ground-truth frame shifted by -2 / -1 / 0 / +1 / +2'], crow), '']
    erow = []
    for s in sets:
        r_ = R[s]
        erow.append([s] + ['%.1f' % v for v in r_['meta|excess_mean']])
    L += ['Mean excess (DN) of the 3 x 3 pixels around the bright ground-truth positions over the frame median, 12 development frames, the ground-truth frame looked up in the image with the given offset '
          '(0 = frame k of the stack with the ground-truth rows of file frame k + 1):', '', table(['set'] + ['offset %+d' % o for o in R[sets[0]]['meta|excess_offsets']], erow), '']
    # ---- per-set evaluation
    evs = {s: evaluate_set(R[s], s, boot) for s in sets}
    summary, tuning, secs = [], [], []
    sections = []
    for s in sets:
        lines, sm, tn = set_section(R[s], s, evs[s], boot)
        sections += lines
        summary += sm
        tuning.append(tn)
        oE = evs[s]['ops']['E']
        secs.append((s, R[s]['istest.%s|secs' % oE['cfg']], R[s]['istest.%s|cpu' % oE['cfg']], R[s]['pic|time_r0_test'], oE['cfg']))
    # ---- headline tables
    L += ['## Summary: test frames, the headline rows', '',
          'Point estimate [95 % block-bootstrap interval] on the held-out test frames. R0: Picasso\'s defaults; R2E+c / R2J+c: every Picasso option tuned on the development frames for E / J, with the position '
          'offset removed; ISS default: the a priori configuration; ISS E / J: the operating points chosen on the development frames (offset removed). All ground-truth emitters unless noted.', '']
    heads = [('R0', 'R0'), ('R2E+c', 'Picasso R2E+c'), ('R2J+c', 'Picasso R2J+c'), ('ISS-0', 'ISS default'), ('ISS-E', 'ISS E-tuned'), ('ISS-J', 'ISS J-tuned')]
    for label, mk, scale, nd, which in (('efficiency E', 'efficiency', 1.0, 1, 'all'), ('Jaccard %', 'jaccard', 100.0, 1, 'all'), ('RMSE nm', 'rmse_nm', 1.0, 1, 'all'),
                                        ('efficiency E, challenge setting', 'efficiency', 1.0, 1, 'ch')):
        rr = []
        for k, kl in heads:
            if which == 'all':
                rr.append([kl] + [fmt_ci(evs[s]['res'][k][mk], evs[s]['ci'][k][mk], scale, nd) for s in sets])
            else:
                rr.append([kl] + [fmt_ci(evs[s]['res2'][k][mk], evs[s]['ci2'][k][mk], scale, nd) for s in sets])
        L += ['%s (test frames):' % label, '', table(['row'] + [SHORT[s] for s in sets], rr), '']
    L += ISS.summary_lines(summary)
    L += ISS.tuning_lines(tuning)
    L += ['## Seconds per frame', '', 'One core of a shared machine (the load average was 10 to 12 on 4 cores during the runs, so wall-clock seconds are inflated; the ISS process CPU time is given too). E-tuned '
          'ISS configuration of the set; Picasso with its defaults (R0, identification and fit, JIT compilation excluded, wall clock).', '']
    L.append(table(['set', 'ISS config', 'ISS wall s / frame (mean, median)', 'ISS CPU s / frame', 'Picasso R0 s / frame', 'test frames (ISS)'],
                   [[SHORT[s], cfg, '%.2f, %.2f' % (sc.mean(), np.median(sc)), '%.2f' % cp.mean(), '%.4f' % (t[0] / t[1]), len(sc)] for s, sc, cp, t, cfg in secs]))
    L.append('')
    neval = [sum(r['n_eval'] for r in _json(R[s_], 'pic|r2').values()) for s_ in sets]
    L += ['## What this experiment does not establish', '',
          '* **The PSF family.** The simulated PSF is a Gaussian, which the ISS kernels (pixel-integrated Gaussians of sliding width) contain exactly, while Picasso fits a point-sampled Gaussian: on the '
          '2016 sets, whose PSF has a halo that Picasso\'s Gaussian does not contain, the ISS\'s kernels (a core plus a halo part, worth a point or two of efficiency on the development frames) match the PSF better than Picasso\'s, so the match is not specific to the 2013 sets; '
          'whether it decides the 2013 wins was not tested in this experiment. A real PSF is neither. The measured width and photon scale of the images also differ from the shipped '
          'PSF file (see each set\'s calibration paragraph); the cause was not investigated.',
          '* **Four sets, two structures, two noise levels.** The two tubulin sets share the ground truth and differ in the noise only; the two bundle sets share the 81,049 fluorophore positions and differ in the '
          'density. They are not four independent samples, and the dense set (224 emitters per 64 x 64 frame) needed its own recipe (declared in the docstring), so its row is a statement about that recipe.',
          '* **Frame ranges.** The long sequences are used in part (declared ranges); the dense set has 60 test frames (6 blocks of the block bootstrap). The intervals cover the sampling of frames only.',
          '* **The efficiency of the dense set** is dominated by its RMSE (pairs within 250 nm of an overlapped neighbour), not by detection; a method that matches more emitters is penalised on RMSE.',
          '* **The challenge setting** (450 nm border, dimmest quarter dropped after pairing) is the harness\'s reconstruction of the 2016 rules; whether the 2013 assessment used the same is not known.',
          '* **The noise law\'s slope** (1.41, the page\'s EM factor) was not calibrated on the data (the tubulin sets give 1.5-1.6 by a free fit); the sensitivity of the result to it was not tested. Only the constant b0 '
          'is measured. The background is a static surface from the development frames (its stability is reported per set).',
          '* **Tuning effort.** The ISS searches 84 operating points per set (no second round: entry width and resolution-limit variants were not run), Picasso %d to %d scored configurations; both on the same '
          'development frames.' % (min(neval), max(neval)) +
          ' The ISS recipe changes (tiles, scaled lengths, the dense recipe) were forced by the data and declared before the test runs, but they are choices.',
          '* **Time.** The ISS is two to three orders of magnitude slower per frame than Picasso (table above); the timings come from a shared machine and mix wall-clock and CPU time as stated.', '']
    return '\n'.join(L + sections) + '\n'


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('--run', action='store_true', help='compute the parts of the sets, assemble the raw file, write the report')
    ap.add_argument('--assemble', action='store_true', help='assemble the raw file from the parts only')
    ap.add_argument('--from-raw', action='store_true', help='regenerate the report from the raw file alone')
    ap.add_argument('--sets', nargs='+', default=list(SETS), choices=list(SETS))
    ap.add_argument('--only', nargs='+', choices=['meta', 'cal', 'psf', 'pic', 'isdev', 'diag', 'istest'], help='run only these parts')
    ap.add_argument('--force', action='store_true', help='recompute parts that exist')
    ap.add_argument('--reverse', action='store_true', help='run the test configurations in reverse order (two processes working on one set from both ends)')
    ap.add_argument('--workers', type=int, default=1, help='processes of the ISS runs (the machine allows 2 busy processes in all)')
    ap.add_argument('--boot', type=int, default=BOOT)
    ap.add_argument('--picasso-python', default=os.environ.get('SMLM_PICASSO_PYTHON'))
    ap.add_argument('--cache-dir', default=None, help='development only: keep / reuse the result of the R2 search per set here')
    ap.add_argument('--parts-dir', default=PARTS)
    ap.add_argument('--raw', default=RAW)
    ap.add_argument('--md', default=MD)
    ap.add_argument('--smoke', action='store_true', help='development of the script only: tiny frame ranges and a 2-setting structural grid; needs explicit --parts-dir, --raw and --md')
    A = ap.parse_args(argv)
    if A.smoke:
        if A.parts_dir == PARTS or A.raw == RAW or A.md == MD:
            ap.error('--smoke needs explicit --parts-dir, --raw and --md (it must never write the real results)')
        for nm in SETS:
            S.PROTOCOL[nm] = dict(dev=(0, 4), test=(10, 14))
        PIC.STRUCT_GRID = collections.OrderedDict([('box', (7,)), ('fit_box', (7,)), ('gauss_sigma', (None, 1.0)), ('pad', (0,)), ('fit', ('gausslq',))])
    if A.run:
        for name in A.sets:
            run_set(name, A, log=lambda s: print(s, flush=True))
        assemble(SETS, A.parts_dir, A.raw)
        R = load_raw(A.raw)
        done = [s for s in SETS if s in R and 'istest.%s|frames' % DEFAULT['cfg'] in R[s]]
        A.from_raw = bool(done) and all(s in done for s in A.sets) and not A.only
        if not A.from_raw:
            print('parts of %s are complete; the report is written when every requested set has all its parts (run without --only)' % (', '.join(done) or 'no set'))
    elif A.assemble:
        assemble(SETS, A.parts_dir, A.raw)
    if A.from_raw:
        R = load_raw(A.raw)
        with open(A.md, 'w') as f:
            f.write(report(R, A.boot, SETS))
        print('wrote', A.md)
    if not (A.run or A.assemble or A.from_raw):
        ap.error('give --run, --assemble or --from-raw')


if __name__ == '__main__':
    main()
