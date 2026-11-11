"""
Controls for the astronomy results: the PSF stamp against the stars, exact-stamp synthetic worlds, and the splitting and flux bias of bright stars (POST HOC DIAGNOSTICS)
=======================================================================================================================================================================

Why. Referee round IS1b-A (`review/REVIEW_roundIS1b_A.md`, findings 1, 2 and 8) built three controls that lived only in the referee's private working directory. This script ports them into the repository, in
its own structure (the referee's scripts were read as data; in scratch runs the renderer was checked bit for bit against the referee's frames with the referee's seed, and the selection and the statistics against
the figures he quotes; the tests pin the definitions), so that the numbers a draft quotes come from a committed, tested script. They are DIAGNOSTICS, run AFTER the transfer test and the kernel study, and they
use the truth list (HUGS) and the PSF stamps. Nothing in the project's methods, kernel rule, baselines or frozen configuration (`results/astro/freeze.json`, `freeze_repeat.json`) changes on their strength, no file
that the freeze records hash is touched, and no ISS run in this script sees a real frame: the ISS and the baselines run only on frames RENDERED here from a stamp and the truth list.

DESIGN (declared before the synthetic runs; every configuration is run once; nothing below depends on a result)

Patches.  Item 1 (stamps): all 18 patches with a stamp and a truth list: M2 (the four test windows (310,630), (310,730), (360,1040), (100,1100) and the two development patches (110,630), (210,630)), M53 (four),
M92 (two) and the six repeat patches of `astro_repeat.REPEAT` (frames m2r1, m2r2, m2r3). Item 2 (worlds): M2 (310,630), M53 (660,420), M92 (1100,220), M92 (1430,550) (the referee's four) and the six repeat
patches: ten worlds. Item 3 on the real catalogues: the ten test patches of M2, M53, M92 and the six repeat patches; on the synthetic ones: the ten worlds.

Item 2, the exact-stamp world of a patch.  Every truth source with F606W < 30 inside the patch box grown by 20 px is rendered with the patch's own 51 x 51 stamp (`psf51`: negative pixels clipped, unit sum, the
stamp the methods get), moved to its sub-pixel offset from the nearest pixel by a cubic spline (`scipy.ndimage.shift`, order 3, zero outside the stamp); flux = `astro_common.mag_to_counts(m)` (the law the scorer
inverts: zero point 22.5, 0.00546689 nanomaggies per count; the frame's own calibration differs from that constant by 0.00 to 0.08 mag, which the real data carry and the worlds do not); the pipeline sky of the
patch (`sky_dn`, per pixel) is added; Poisson noise at the patch gain (electrons = DN x gain) and Gaussian noise of variance `dark_variance` (DN^2). So the PSF is exactly the stamp, the noise model the methods
assume is exactly right, and the injected truth (the HUGS sources of the box, scored by `astro_common.truth_for_patch` as for the real data) is exact: no astrometric offset, no colour term, no missing truth, no
source outside the list, no PSF variation over the patch, no wing beyond 25 px, no saturation. The noise realisation of each world is `default_rng([1108, crc32(set name), x0, y0])`: one realisation per world.
Run on each world: (a) the frozen ISS, `exp_astro_m2.run_path` (isotropic kernel) at the configuration of `freeze.json` (sigma_b 3, r_ex 1.2, kappa 2, with the BIC key); (b) the ISS with the kernel that
`psf_general.fit_adequate` chooses from the stamp, `exp_astro_psf.run_path_psf` (where the rule keeps the isotropic fit the two are the same program on the same input: not run again, shown as "same");
(c) crowdsource at its default (`exp_astro_baselines.run_crowdsource`, then `exp_astro_m2.cs_dets` with the default filter, and with the frozen post-filter as a second row); then, after those, (d) photutils at the
variant that `freeze.json` names (`photutils_ext_all_local_i6`) and (e) the ISS with both kernels at the a priori default (sigma_b 5, r_ex 1.2, kappa 5). Scored against the injected truth with
`astro_common.score`: F of completeness and purity at mag < 22, 20 x 20 px blocks, intervals as in `exp_astro_m2` (95 % percentile intervals, 2000 resamples, seed 0, the same resamples for every row of a table,
stratified by patch when pooled), conditional on the patch and on the one noise realisation. Beside them the real-data F of the same patches, copied from the frozen results (`m2_raw.npz`: ISS as frozen and the
baselines; `psf_raw.npz`: the adequate kernel on M53 and M92; `repeat_raw.npz`: the repeat patches), and the residual scale (`exp_astro_psf.residual_sd`, robust sd of the standardized residual) in both worlds.

Item 1, the stamp against the stars (label-using).  Stars of a patch: truth sources with F606W < 20 whose position is more than 12 px from the patch edge (so that a 25 x 25 cut-out around the nearest pixel lies in
the patch) and that have no other truth source (any magnitude < 30, the box grown by 20 px) within 4.5 px that is brighter than m + 1.5. Each star: the 25 x 25 cut-out around the nearest pixel, minus the median
of its pixels farther than 8.5 px from the middle, then moved by a cubic spline (`order=3`, edge values repeated) so that its flux centroid (iterated from the middle pixel, radius 3 px, 6 iterations, negative
pixels counted as zero) is the middle pixel, cropped to 19 x 19 and divided by its sum within 6 px. The stack is the pixelwise MEDIAN over the stars; the stamp goes through the same centring and normalisation
(25 x 25 core, no sky). Reported: the central pixel of the stack against the stamp's, relative difference; the FWHM (`astro_second.psf_fwhm`: twice the radius at which the azimuthal mean falls to half the
middle pixel). The same procedure on the synthetic frame of the patch (the same stars) shows what the centring does to a star whose PSF IS the stamp, so the statistic is real minus synthetic (the
difference of the two relative differences, and the ratio of the central pixels minus 1, which is the figure the referee quoted). Uncertainty: 95 % percentile intervals of a bootstrap over the stars (2000
resamples, the same stars for the real and the synthetic stack, stratified by patch when stacks are pooled; a pooled stack is compared with the stamps of its patches averaged with the stars' numbers as weights).
A paired statistic was added after the first run of this part (which is cheap and was run before the worlds): the median over the stars of the ratio of the star's two central pixels (real over synthetic), because
the bootstrap interval of a pixelwise median of a few tens of stars is wide and skewed by a few stars with a companion just beyond 4.5 px or a crowd-inflated sky ring; the first statistics did not change.

Item 3, splitting and flux bias of bright stars (label-using).  A bright isolated star: a truth source in the scored area (at least 2.5 px from the patch edge) with F606W < 20.5 and no other truth source (the
box grown by 20 px) within 3 px that is brighter than m_star + 2.5. Its units: the detections of positive flux within 2 px of its truth position; a significant unit has at least 10 % of the star's flux
(`mag_to_counts(m_star)`); a split is at least 2 significant units, a miss none. The flux ratio is the units' summed flux over the star's. The bright pairs: a one-to-one position-only matching (`astro_common.match`,
0.7 px, no magnitude condition) of the detections in the scored area to the scored truth, the pairs with truth F606W < 19.5; their median detection-minus-truth magnitude and position offset. For the ISS (tuned
and a priori default, isotropic and chosen kernel), crowdsource (default and frozen post-filter), photutils (tuned): the real catalogues beside the synthetic ones. Intervals: bootstrap over the stars (pairs) as above.

Stalled paths (added to the design after the synthetic runs had begun, on a report of the project's review; the runs themselves are unchanged). An ISS path whose entry count reached the cap
(`exp_astro_m2.MAX_ENTRIES`; its `bic` record, block of the BIC minimum / blocks of the path / flag, ends in cap + 1 blocks) has stalled: `imaging.sliding_path` can enter and exit one lattice point again and again
at the same entry time until the cap, so that the run stops after a few hundred distinct entries whatever its stopping rule says. It is a defective run, not a censored valid one. `src/core/imaging.py` is hash-frozen
and is not touched. Here: the `bic` records stay in the raw (and are copied from the frozen results for the real runs; item 3 stores a `stalled` flag beside every catalogue); the report marks every row computed from
such a path STALLED, leaves it out of every pool and difference and names it, and states how many paths of the worlds and of the frozen results reached the cap, also when none did.

These controls identify nothing that was not put into the world: where the synthetic world reproduces a real number, the ingredients of the world are enough for it; where it does not, something that is not in
the world (and not in the stamp) is at work. The limitations are printed in the report.

    python experiments/exp_astro_diag.py --stamps                  item 1, all patches (minutes) -> results/astro/diag_raw.npz
    <domain python> -I experiments/exp_astro_diag.py --synth [--workers 2] [--only iso,adq,cs,phot,iso_def,adq_def] [--sets m2,m53,...]
                                                                   item 2 (hours; the domain environment holds crowdsource and photutils; a run that has a result is not run again; refuses unless the freeze matches)
    python experiments/exp_astro_diag.py --real                    the frozen results of the patches of the worlds copied into the raw (--synth does it too; run it again if the frozen raws change)
    python experiments/exp_astro_diag.py --bright                  item 3 from the stored catalogues (real: the frozen raws and baseline catalogues; synthetic: the worlds' catalogues)
    python experiments/exp_astro_diag.py --from-raw                regenerates results/astro/diag.md from diag_raw.npz alone (byte-identically)

Record of the committed raw.  `diag_raw.npz` was produced on a shared machine by single-worker calls of `--synth` on disjoint `--sets` (the first call ran the first stage of the referee's four worlds; later
calls, two processes at once, ran `--sets m2,m53,m92` and `--sets m2r1,m2r2,m2r3`), the keys of the repeat frames being copied into the raw of the first before `--real`, `--stamps` and `--bright` were run on
it. The runs are independent of each other (own seed, own input, no shared state), so one call `--synth --workers 2` is meant to give the same results; only the timings (`secs`) differ between calls. The
production calls did not use `--workers`; the pool is unit-tested with a fake pool, and was run once in a scratch directory (`--only cs --sets m2,m53 --workers 2`, domain environment), where the 40
result arrays of the two crowdsource runs, the 80 world headers and the 190 arrays of the real copy came out identical to those of the committed raw (NaN equal to NaN).
"""
import argparse
import hashlib
import json
import os
import sys
import zlib

os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
from scipy.ndimage import shift as nd_shift
from scipy.spatial import cKDTree

HERE = os.path.abspath(os.path.dirname(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)
import astro_common as A                  # noqa: E402
import astro_repeat as R                  # noqa: E402
import astro_second as S2                 # noqa: E402
import exp_astro_m2 as X                  # noqa: E402
import exp_astro_psf as PS                # noqa: E402
from src.core import imaging as IM        # noqa: E402
from src.core import psf_general as G     # noqa: E402

OUT = os.environ.get('DIAG_OUT', os.path.join(REPO, 'results', 'astro'))
RAW, MD = os.path.join(OUT, 'diag_raw.npz'), os.path.join(OUT, 'diag.md')
RES = os.path.join(REPO, 'results', 'astro')
M2RAW = os.environ.get('DIAG_M2_RAW', X.RAW)                                    # the frozen results, read only
PSFRAW = os.environ.get('DIAG_PSF_RAW', os.path.join(RES, 'psf_raw.npz'))
REPEATRAW = os.environ.get('DIAG_REPEAT_RAW', os.path.join(RES, 'repeat_raw.npz'))
BASE_REPEAT = os.path.join(RES, 'baseline_catalogs_repeat.npz')

DIAG_SEED = 1108                         # the noise of a world is default_rng([DIAG_SEED, crc32(set name), x0, y0])
N = A.CELL * 5                           # the patch side, 100 px
MARGIN, MAG_MAX = 20, 30.0               # px of truth rendered outside the box; every truth source brighter than this is rendered
STACK_MAG, STACK_RAD, STACK_DM = 20.0, 4.5, 1.5          # item 1: the stars (F606W limit; no source within STACK_RAD px brighter than m + STACK_DM)
HALF, PAD = 9, 3                         # the stacked cut-out is (2 HALF + 1)^2; PAD more pixels for the sky ring and the shift
SKY_RING, CEN_RAD, CEN_ITER, NORM_RAD = 8.5, 3.0, 6, 6.0
BRIGHT_MAG, BRIGHT_RAD, BRIGHT_DM = 20.5, 3.0, 2.5       # item 3: bright isolated stars
UNIT_RAD, UNIT_FRAC = 2.0, 0.1           # the units of a star: detections within UNIT_RAD px; significant: at least UNIT_FRAC of its flux
PAIR_RAD, PAIR_MAG = 0.7, 19.5           # bright pairs: position-only matching radius, truth magnitude limit

STAMP_SETS = ('m2', 'm53', 'm92', 'm2r1', 'm2r2', 'm2r3')
REPEAT_SETS = ('m2r1', 'm2r2', 'm2r3')
SYNTH = [('m2', (310, 630)), ('m53', (660, 420)), ('m92', (1100, 220)), ('m92', (1430, 550))] + [(n, tuple(p)) for n in REPEAT_SETS for p in R.REPEAT[n]['patches']]
POOLS = [('M2 (310, 630)', [p for p in SYNTH if p[0] == 'm2']), ('M53 (660, 420)', [p for p in SYNTH if p[0] == 'm53']), ('M92 (two windows of one frame)', [p for p in SYNTH if p[0] == 'm92']),
         ('M2 repeat frame m2r1', [p for p in SYNTH if p[0] == 'm2r1']), ('M2 repeat frame m2r2', [p for p in SYNTH if p[0] == 'm2r2']), ('M2 repeat frame m2r3', [p for p in SYNTH if p[0] == 'm2r3']),
         ('M2 repeat frames, all six patches', [p for p in SYNTH if p[0] in REPEAT_SETS])]
# the methods of a table: short name -> key suffix (the same in the worlds 'W|' and in the copies of the real results 'RF|'); a suffix names the results of one run and its readout
METHODS = {'iso': 'iso|ours|k2', 'adq': 'adq|ours|k2', 'iso_def': 'iso|default|k5', 'adq_def': 'adq|default|k5', 'cs': 'cs|default', 'cs_tuned': 'cs|tuned', 'phot': 'phot|tuned'}
METHOD_LABEL = {'iso': 'ISS, isotropic kernel, tuned (kappa 2)', 'adq': 'ISS, chosen kernel, tuned (kappa 2)', 'iso_def': 'ISS, isotropic kernel, a priori default (kappa 5)',
                'adq_def': 'ISS, chosen kernel, a priori default (kappa 5)', 'cs': 'crowdsource (default)', 'cs_tuned': 'crowdsource (frozen post-filter)', 'phot': 'photutils (tuned)'}
JOB_TAG = {'iso': 'iso|ours', 'adq': 'adq|ours', 'iso_def': 'iso|default', 'adq_def': 'adq|default', 'cs': 'cs', 'phot': 'phot'}      # a run's results in a world: 'W|<set>|x<x0>_y<y0>|<tag>|...'
TAG_METHOD = {'iso|ours': 'iso', 'adq|ours': 'adq', 'iso|default': 'iso_def', 'adq|default': 'adq_def'}                           # the ISS runs: tag -> method name
STAGES = [(('iso', 'adq', 'cs'), [p for p in SYNTH if p[0] not in REPEAT_SETS]), (('iso', 'adq', 'cs'), [p for p in SYNTH if p[0] in REPEAT_SETS]), (('phot',), SYNTH), (('iso_def', 'adq_def'), SYNTH)]
BRIGHT_ORDER = ('iso', 'iso_def', 'adq', 'adq_def', 'cs', 'cs_tuned', 'phot')


def set_patches(name):
    """The patches (x0, y0) of a set: 'm2' (the four test windows and the two development patches), 'm53', 'm92', and the repeat frames 'm2r1', 'm2r2', 'm2r3'."""
    if name == 'm2':
        return [tuple(p) for p in X.TESTS] + [tuple(p) for p in X.DEV]
    if name in R.REPEAT:
        return [tuple(p) for p in R.REPEAT[name]['patches']]
    return [tuple(p) for p in X.SECOND[name]['patches']]


def set_title(name):
    if name == 'm2':
        return 'M2 (NGC 7089), SDSS run 2583 camcol 2 field 136'
    return R.REPEAT[name]['title'] if name in R.REPEAT else X.SECOND[name]['title']


def truth_of(name):
    """dict(px, py, m606) of the truth list of the frame of a set."""
    return R.truth_arrays(name) if name in R.REPEAT else X.truth_arrays(name)


def patch_file(name, pxy):
    return R.patch_path(pxy[0], pxy[1], name) if name in R.REPEAT else X.patch_path(pxy[0], pxy[1], name)


def load_inputs(name, pxy):
    """(patch dict(dn, psf, gain, dark) as `exp_astro_m2.load_patch`, sky_dn (the pipeline sky of the patch), truth dict(px, py, m606) of the frame, psField sky sigma) of a patch."""
    path = patch_file(name, pxy)
    z = np.load(path)
    patch = dict(dn=z['dn'].astype(float), psf=z['psf51'].astype(float), gain=float(z['gain']), dark=float(z['dark_variance']))
    with open(os.path.join(os.path.dirname(path), 'meta.json')) as f:
        skysig = float(json.load(f)['psfield_skysig'])
    return patch, z['sky_dn'].astype(float), truth_of(name), skysig


# ----------------------------------------------------------------------------------------------------------------------------------
# the exact-stamp world
# ----------------------------------------------------------------------------------------------------------------------------------
def sources_near(truth, x0, y0, margin=0, n=N, mag_max=MAG_MAX):
    """(x, y, mag) in the patch's pixel coordinates (column, row; pixel centres at the integers) of the truth sources brighter than mag_max inside the box [x0 - 0.5 - margin, x0 + n - 0.5 + margin)^2
    (margin 0: the patch's own box, the sources `astro_common.truth_for_patch` scores)."""
    px, py, m = truth['px'] - x0, truth['py'] - y0, truth['m606']
    sel = (px >= -0.5 - margin) & (px < n - 0.5 + margin) & (py >= -0.5 - margin) & (py < n - 0.5 + margin) & (m < mag_max)
    return px[sel], py[sel], m[sel]


def render_sources(stamp, x, y, flux, shape, margin=MARGIN):
    """The noise-free image (`shape`) of sources at (x, y) (patch pixel coordinates) of flux `flux` (DN): each is the unit-sum `stamp` (its middle pixel is the source's centre) moved to the source's
    sub-pixel offset from the nearest pixel by a cubic spline (`scipy.ndimage.shift`, order 3, zero outside) and added; sources up to `margin` px outside the image contribute their wings."""
    ny, nx = shape
    canvas = np.zeros((ny + 2 * margin, nx + 2 * margin))
    h = stamp.shape[0] // 2
    for xi, yi, fi in zip(x, y, flux):
        ix, iy = int(np.round(xi)), int(np.round(yi))
        fx, fy = xi - ix, yi - iy
        st = nd_shift(stamp, (fy, fx), order=3, mode='constant', cval=0.0) if (abs(fx) > 1e-9 or abs(fy) > 1e-9) else stamp
        ya, xa = iy + margin - h, ix + margin - h
        y0_, x0_ = max(ya, 0), max(xa, 0)
        y1_, x1_ = min(ya + stamp.shape[0], canvas.shape[0]), min(xa + stamp.shape[1], canvas.shape[1])
        if y1_ <= y0_ or x1_ <= x0_:
            continue
        canvas[y0_:y1_, x0_:x1_] += fi * st[y0_ - ya: y1_ - ya, x0_ - xa: x1_ - xa]
    return canvas[margin: margin + ny, margin: margin + nx]


def add_noise(clean, gain, dark, rng):
    """DN with the noise model of the methods: Poisson electrons (DN x gain) and Gaussian noise of variance `dark` (DN^2)."""
    return rng.poisson(np.maximum(clean, 0.0) * gain) / gain + rng.normal(0.0, np.sqrt(dark), clean.shape)


def world_seed(name, pxy):
    return [DIAG_SEED, zlib.crc32(name.encode()), int(pxy[0]), int(pxy[1])]


def render_world(truth, x0, y0, stamp, gain, dark, sky, seed, margin=MARGIN, n=N, noise=True):
    """The exact-stamp world of a patch: dict(dn (the frame), clean (the frame without noise), x, y, mag (every source rendered: the box grown by `margin`))."""
    x, y, mag = sources_near(truth, x0, y0, margin, n)
    clean = render_sources(stamp, x, y, A.mag_to_counts(mag), (n, n), margin) + sky
    dn = add_noise(clean, gain, dark, np.random.default_rng(seed)) if noise else clean.copy()
    return dict(dn=dn, clean=clean, x=x, y=y, mag=mag)


def frame_stats(dn):
    """median, 5 and 95 percentile and maximum of a frame."""
    return np.array([np.median(dn), np.percentile(dn, 5), np.percentile(dn, 95), dn.max()])


# ----------------------------------------------------------------------------------------------------------------------------------
# item 1: the stack of isolated bright stars against the stamp
# ----------------------------------------------------------------------------------------------------------------------------------
DISK = np.hypot(*np.mgrid[-HALF:HALF + 1, -HALF:HALF + 1]) <= NORM_RAD


def isolated(x, y, m, cand, rad, dm):
    """Indices (into x, y, m) of the candidates (boolean mask `cand`) that have no OTHER source within `rad` px (distance at most rad) with a magnitude below m + dm (brighter than m + dm)."""
    tree = cKDTree(np.c_[x, y])
    out = []
    for i in np.flatnonzero(cand):
        if not any(j != i and m[j] < m[i] + dm for j in tree.query_ball_point([x[i], y[i]], rad)):
            out.append(i)
    return np.array(out, int)


def stack_stars(x, y, m, n=N, mag=STACK_MAG, rad=STACK_RAD, dm=STACK_DM, edge=HALF + PAD):
    """The stars of item 1 among the sources (x, y, m) of a patch and its margin: m < mag, more than `edge` px from every edge of the n x n patch, isolated (`isolated`)."""
    cand = (m < mag) & (x > edge) & (x < n - 1 - edge) & (y > edge) & (y < n - 1 - edge)
    return isolated(x, y, m, cand, rad, dm)


def centroid(img, rad=CEN_RAD, iters=CEN_ITER):
    """The flux centroid (column, row) of `img` inside `rad` pixels of the current estimate, iterated from the middle pixel; negative pixels count as zero."""
    n = img.shape[0]
    yy, xx = np.mgrid[0:n, 0:n]
    cx = cy = (n - 1) / 2.0
    for _ in range(iters):
        w = np.clip(img * (np.hypot(yy - cy, xx - cx) <= rad), 0.0, None)
        s = w.sum()
        if s <= 0:
            break
        cx, cy = (w * xx).sum() / s, (w * yy).sum() / s
    return cx, cy


def aligned(sub):
    """The (2 (HALF + PAD) + 1)^2 cut-out `sub` with its flux centroid moved to the middle pixel (cubic spline), cropped to (2 HALF + 1)^2 and divided by its sum within NORM_RAD."""
    n = sub.shape[0]
    cx, cy = centroid(sub)
    sh = nd_shift(sub, ((n - 1) / 2.0 - cy, (n - 1) / 2.0 - cx), order=3, mode='nearest')[PAD:-PAD, PAD:-PAD]
    return sh / sh[DISK].sum()


def cutouts(dn, x, y, idx):
    """The aligned cut-outs (len(idx), 2 HALF + 1, 2 HALF + 1) of the stars idx of the frame `dn`, the sky (median of the pixels beyond SKY_RING px of the middle) subtracted."""
    r = HALF + PAD
    yy, xx = np.mgrid[-r:r + 1, -r:r + 1]
    ring = np.hypot(yy, xx) > SKY_RING
    out = []
    for i in idx:
        ix, iy = int(np.round(x[i])), int(np.round(y[i]))
        sub = dn[iy - r: iy + r + 1, ix - r: ix + r + 1].copy()
        sub -= np.median(sub[ring])
        out.append(aligned(sub))
    return np.array(out).reshape(len(idx), 2 * HALF + 1, 2 * HALF + 1)


def stamp_core(stamp):
    """The stamp (middle pixel = the star's centre) through the same centring and normalisation as a star, without sky."""
    c, r = stamp.shape[0] // 2, HALF + PAD
    return aligned(stamp[c - r: c + r + 1, c - r: c + r + 1])


def fwhm(img):
    """FWHM (px) of a centred cut-out: `astro_second.psf_fwhm` (twice the radius at which the azimuthal mean falls to half the middle pixel); NaN where that is not a number in (0, 2 HALF)."""
    v = S2.psf_fwhm(img)
    return float(v) if np.isfinite(v) and 0.0 < v < 2 * HALF else float('nan')


STAT_NAMES = ('centre, real / stamp - 1', 'centre, synthetic / stamp - 1', 'centre, real minus synthetic', 'centre, real / synthetic - 1', 'FWHM real', 'FWHM synthetic', 'FWHM real - synthetic',
              'FWHM real / synthetic - 1', 'centre, median over the stars of (real / synthetic) - 1')


def stack_estimates(real, syn, ref):
    """The statistics of item 1 (STAT_NAMES; relative differences as fractions) from the aligned cut-outs (n, 19, 19) of the same n stars in the real and in the synthetic frame and the reference (the
    aligned stamp): from their median stacks, and the median over the stars of the ratio of the star's two central pixels (a paired estimator that a few stars with extreme values do not move)."""
    med_real, med_syn = np.median(real, 0), np.median(syn, 0)
    c = med_real.shape[0] // 2
    cr, cs, c0 = med_real[c, c], med_syn[c, c], ref[c, c]
    fr, fs = fwhm(med_real), fwhm(med_syn)
    rr, rs = cr / c0 - 1.0, cs / c0 - 1.0
    return np.array([rr, rs, rr - rs, cr / cs - 1.0, fr, fs, fr - fs, fr / fs - 1.0, np.median(real[:, c, c] / syn[:, c, c]) - 1.0])


def stack_bootstrap(real, syn, ref, sizes, B=X.BOOT, seed=0):
    """(estimate (9,), 95 % percentile interval (9, 2), sd (9,), number of resamples with a NaN) of `stack_estimates` for the cut-outs real, syn (n, 19, 19; the same n stars) and the reference `ref`;
    the stars are resampled with replacement within each patch (`sizes`: the stars per patch, in the order of the cut-outs), the same resamples for the real and the synthetic stack."""
    est = stack_estimates(real, syn, ref)
    idx = X.boot_idx_strat([int(s) for s in sizes if s > 0], B, seed)
    draws = np.array([stack_estimates(real[i], syn[i], ref) for i in idx])
    return est, np.nanpercentile(draws, [2.5, 97.5], axis=0).T, np.nanstd(draws, axis=0), int(np.isnan(draws).any(axis=1).sum())


def stack_arrays(real, syn, ref, sizes, seed, B=X.BOOT):
    """The arrays stored for one patch or one pooled set: stars, the median stacks and the aligned stamp, the estimates, their intervals and sd."""
    est, ci, sd, nbad = stack_bootstrap(real, syn, ref, sizes, B, seed)
    return dict(n=np.array([len(real)]), sizes=np.array(sizes), est=est, ci=ci, sd=sd, nan=np.array([nbad]), real=np.median(real, 0), syn=np.median(syn, 0), ref=ref)


def part_stamps(out_path=RAW, B=X.BOOT):
    """Item 1 for every patch of STAMP_SETS and for each set pooled (the stars of its patches); stored under 'S|<set>|x<x0>_y<y0>' and 'S|<set>|all'."""
    for si, name in enumerate(STAMP_SETS):
        reals, syns, refs, sizes = [], [], [], []
        for pi, pxy in enumerate(set_patches(name)):
            patch, sky, truth, _ = load_inputs(name, pxy)
            w = render_world(truth, pxy[0], pxy[1], patch['psf'], patch['gain'], patch['dark'], sky, world_seed(name, pxy))
            idx = stack_stars(w['x'], w['y'], w['mag'])
            real, syn, ref = cutouts(patch['dn'], w['x'], w['y'], idx), cutouts(w['dn'], w['x'], w['y'], idx), stamp_core(patch['psf'])
            if len(idx) < 3:
                print('S %s %s: %d stars, not stacked' % (name, pxy, len(idx)), flush=True)
                continue
            arrays = stack_arrays(real, syn, ref, [len(idx)], 100 * si + pi, B)
            X.save_part(out_path, 'S|%s|x%d_y%d' % (name, pxy[0], pxy[1]), arrays)
            reals.append(real), syns.append(syn), refs.append(ref), sizes.append(len(idx))
            print('S %s %s: %d stars; centre real / stamp %+.3f, synthetic / stamp %+.3f; FWHM real %.2f synthetic %.2f stamp %.2f' % (
                name, pxy, len(idx), arrays['est'][0], arrays['est'][1], arrays['est'][4], arrays['est'][5], fwhm(ref)), flush=True)
        if sizes:
            ref = sum(n * r for n, r in zip(sizes, refs)) / sum(sizes)
            arrays = stack_arrays(np.concatenate(reals), np.concatenate(syns), ref, sizes, 100 * si + 99, B)
            X.save_part(out_path, 'S|%s|all' % name, arrays)
            print('S %s pooled: %d stars; centre real minus synthetic %+.3f' % (name, arrays['n'][0], arrays['est'][2]), flush=True)


# ----------------------------------------------------------------------------------------------------------------------------------
# item 3: splitting and flux bias of bright stars
# ----------------------------------------------------------------------------------------------------------------------------------
def inside_scored(x, y, n=N):
    """True where a position is in the scored area (at least `astro_common.EDGE` px from the patch edge)."""
    return (x >= A.EDGE - 0.5) & (x <= n - 0.5 - A.EDGE) & (y >= A.EDGE - 0.5) & (y <= n - 0.5 - A.EDGE)


def bright_isolated(x, y, m, n=N):
    """The bright isolated stars among the sources (x, y, m) of a patch and its margin: in the scored area, m < BRIGHT_MAG, no other source within BRIGHT_RAD px brighter than m + BRIGHT_DM."""
    return isolated(x, y, m, inside_scored(x, y, n) & (m < BRIGHT_MAG), BRIGHT_RAD, BRIGHT_DM)


def star_units(det, x, y, m, idx):
    """For each star idx (indices into x, y, m) the detections `det` = (x, y, flux) of positive flux within UNIT_RAD px of it: (number of units, number of significant units (flux at least UNIT_FRAC
    of the star's, `mag_to_counts(m)`), summed flux of the units over the star's flux)."""
    dx, dy, f = [np.asarray(a, float) for a in det]
    k = f > 0
    dx, dy, f = dx[k], dy[k], f[k]
    nu, ns, fr = np.zeros(len(idx), int), np.zeros(len(idx), int), np.zeros(len(idx))
    if len(dx):
        tree = cKDTree(np.c_[dx, dy])
        for a, i in enumerate(idx):
            near = np.array(tree.query_ball_point([x[i], y[i]], UNIT_RAD), int)
            fl = A.mag_to_counts(m[i])
            nu[a], ns[a], fr[a] = len(near), int((f[near] >= UNIT_FRAC * fl).sum()), f[near].sum() / fl
    return nu, ns, fr


def bright_pairs(det, tx, ty, tm, n=N):
    """(n, 3) detection minus truth (dx, dy, dmag) of the pairs of truth magnitude below PAIR_MAG in a one-to-one, position-only (PAIR_RAD px) matching of the detections with finite magnitude in the
    scored area to the scored truth (tx, ty, tm: the sources of the patch's box)."""
    dx, dy, f = [np.asarray(a, float) for a in det]
    dm = A.counts_to_mag(f)
    k = np.isfinite(dm) & inside_scored(dx, dy, n)
    dx, dy, dm = dx[k], dy[k], dm[k]
    ti = np.flatnonzero(inside_scored(tx, ty, n))
    rd, rt = A.match(np.c_[dx, dy], dm, np.c_[tx[ti], ty[ti]], tm[ti], dpos=PAIR_RAD, dmag=np.inf)
    b = tm[ti][rt] < PAIR_MAG
    return np.c_[dx[rd][b] - tx[ti][rt][b], dy[rd][b] - ty[ti][rt][b], dm[rd][b] - tm[ti][rt][b]].reshape(-1, 3)


def bright_arrays(det, universe, box):
    """The per-star and per-pair arrays of one catalogue of one patch: `universe` = (x, y, m) of the sources of the box grown by MARGIN (the neighbours), `box` = (x, y, m) of the sources of the
    box. Returns dict(nunits, nsig, fratio, pairs)."""
    ux, uy, um = universe
    idx = bright_isolated(ux, uy, um)
    nu, ns, fr = star_units(det, ux, uy, um, idx)
    return dict(nunits=nu, nsig=ns, fratio=fr, pairs=bright_pairs(det, *box))


def boot_median(parts, B=X.BOOT, seed=0):
    """(median, 95 % percentile interval) of the median of the concatenated arrays `parts` (one per patch) under a bootstrap that resamples within each patch (a patch with no value is skipped)."""
    parts = [np.asarray(p, float) for p in parts if len(p)]
    if not parts:
        return float('nan'), np.full(2, np.nan)
    idx = X.boot_idx_strat([len(p) for p in parts], B, seed)
    v = np.concatenate(parts)
    return float(np.median(v)), np.percentile(np.median(v[idx], axis=1), [2.5, 97.5])


def boot_fraction(parts, pred, B=X.BOOT, seed=0):
    """(fraction, 95 % percentile interval) of the values of the concatenated arrays `parts` satisfying `pred`, bootstrap as `boot_median`."""
    parts = [np.asarray(p) for p in parts if len(p)]
    if not parts:
        return float('nan'), np.full(2, np.nan)
    idx = X.boot_idx_strat([len(p) for p in parts], B, seed)
    v = pred(np.concatenate(parts)).astype(float)
    return float(v.mean()), np.percentile(v[idx].mean(axis=1), [2.5, 97.5])


# ----------------------------------------------------------------------------------------------------------------------------------
# item 2: the runs on a world, the copy of the real results
# ----------------------------------------------------------------------------------------------------------------------------------
def iss_arrays(dets, secs, n, info, t, dk, rsd):
    """The arrays of one ISS run on a world (names relative to the run): the score of every readout ('k<kappa>|...', the BIC key is k0), the catalogue at the reported kappa dk ('k<dk>|det'), the run
    time and entrants ('secs'), the BIC block ('bic': index, blocks, 1 if the entry cap was hit) and the residual scale at dk ('rsd')."""
    arrays = {}
    for k, d in dets.items():
        arrays.update(X._pack('k%g' % k, X.score_dets(d, t)))
    arrays['secs'] = np.array([secs, n])
    arrays['bic'] = np.array([info['bic_index'], info['n_blocks'], int(n >= X.MAX_ENTRIES)])
    arrays['k%g|det' % dk] = np.c_[dets[dk][0], dets[dk][1], dets[dk][2]]
    arrays['rsd'] = np.array([rsd])
    return arrays


def cs_dict(cat, pxy):
    """A crowdsource catalogue in the key format of baseline_catalogs.npz, so that `exp_astro_m2.cs_dets` filters it."""
    return {'x%d_y%d|crowdsource|%s' % (pxy[0], pxy[1], k): v for k, v in cat.items()}


def cs_filter(name):
    """The post-filter (min flux / dflux, min fracflux, min qf, max spread) of `exp_astro_m2.cs_name` equal to `name`."""
    return [f for f in X.CS_FILTERS if X.cs_name(f) == name][0]


def cs_arrays(cat, pxy, t, tuned_name, secs):
    """The arrays of crowdsource on a world: the default filter and the frozen post-filter `tuned_name` (as freeze.json names it), each scored, with the catalogues."""
    b = cs_dict(cat, pxy)
    arrays = {'secs': np.array([secs, len(cat['x'])])}
    for label, f in (('default', (0.0, 0.0, 0.0, np.inf)), ('tuned', cs_filter(tuned_name))):
        det = X.cs_dets(b, pxy, f)
        arrays.update(X._pack(label, X.score_dets(det, t)))
        arrays[label + '|det'] = np.c_[det[0], det[1], det[2]]
    return arrays


def phot_arrays(det, t, secs):
    """The arrays of photutils (the tuned variant) on a world."""
    arrays = {'secs': np.array([secs, len(det[0])]), 'tuned|det': np.c_[det[0], det[1], det[2]]}
    arrays.update(X._pack('tuned', X.score_dets(det, t)))
    return arrays


def world_job(args):
    """One run on the world of a patch (a function of its arguments only, so that jobs can run in separate processes): (kind, set, patch, arrays)."""
    kind, name, pxy, cfg = args
    patch, sky, truth, skysig = load_inputs(name, pxy)
    w = render_world(truth, pxy[0], pxy[1], patch['psf'], patch['gain'], patch['dark'], sky, world_seed(name, pxy))
    world = dict(dn=w['dn'], psf=patch['psf'], gain=patch['gain'], dark=patch['dark'])
    t = A.truth_for_patch(truth, pxy[0], pxy[1])
    if kind in ('iso', 'iso_def', 'adq', 'adq_def'):
        (h, r), dk = ((cfg['hp'], cfg['rex']), cfg['kap']) if kind in ('iso', 'adq') else ((X.DEFAULT[0], X.DEFAULT[1]), X.DEFAULT[2])
        if kind.startswith('iso'):
            dets, secs, n, info = X.run_path(world, h, r)
            sig, parts, _ = IM.fit_psf_parts(world['psf'], n_parts=3, sigma0=1.0)
            P_iso = IM.PixelPSF(*world['dn'].shape, sig, parts=parts, hp=h or None)
            rsd = PS.residual_sd(world, P_iso, np.c_[dets[dk][0], dets[dk][1], dets[dk][2]])
        else:
            dets, secs, n, info = PS.run_path_psf(world, h, r, G.fit_adequate(world['psf']))
            rsd = info['rsd'][dk]
        arrays = iss_arrays(dets, secs, n, info, t, dk, rsd)
    elif kind == 'cs':
        cat, secs = X.AB.run_crowdsource(world['dn'], float(sky.mean()), world['psf'], world['gain'], world['dark'])
        arrays = cs_arrays(cat, pxy, t, cfg['tuned']['crowdsource'], secs)
    elif kind == 'phot':
        spec = [d for d in X.AB.PHOT_EXT if d['tag'] == cfg['tuned']['photutils']][0]
        cat, secs, _ = X.AB.run_photutils_ext(world['dn'], world['psf'], world['gain'], world['dark'], skysig, spec)
        arrays = phot_arrays((cat['x'], cat['y'], cat['flux']), t, secs)
    else:
        raise ValueError(kind)
    return kind, name, pxy, arrays


def world_header(name, pxy):
    """The arrays that describe the world of a patch (no run): the stamp fit (kind of kernel, rms, parts), the frame statistics of the real and the synthetic frame, the numbers of sources."""
    patch, sky, truth, _ = load_inputs(name, pxy)
    w = render_world(truth, pxy[0], pxy[1], patch['psf'], patch['gain'], patch['dark'], sky, world_seed(name, pxy))
    t = A.truth_for_patch(truth, pxy[0], pxy[1])
    arrays, fit = PS.stamp_arrays(patch)
    arrays['stats'] = np.array([frame_stats(patch['dn']), frame_stats(w['dn'])])
    arrays['counts'] = np.array([len(w['x']), len(t['x']), int((t['mag'][t['inner']] < 22.0).sum())])
    arrays['fwhm'] = np.array([S2.psf_fwhm(patch['psf'])])
    return arrays


def run_cfg(m2raw=M2RAW):
    """The frozen configuration and tuned baselines, read from the frozen results after `exp_astro_m2.check_freeze`: dict(hp, rex, kap, tuned)."""
    froz = X.load_raw(m2raw)
    X.check_freeze(froz)
    hp, rex, kap, f = X.choose(froz)
    return dict(hp=float(hp), rex=float(rex), kap=float(kap), tuned=X.choose_baselines(froz))


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


CODE_FILES = ('src/core/imaging.py', 'src/core/psf_general.py', 'experiments/astro_common.py', 'experiments/exp_astro_m2.py', 'experiments/exp_astro_psf.py', 'experiments/exp_astro_baselines.py')


def meta_arrays(cfg):
    """What the raw records of its provenance: the configuration read from the freeze and the hashes of the code the runs used."""
    return dict(cfg=np.array([cfg['hp'], cfg['rex'], cfg['kap']]), tuned=np.array([cfg['tuned']['crowdsource'], cfg['tuned']['photutils']]),
                files=np.array(['%s:%s' % (f, sha256_file(os.path.join(REPO, f))) for f in CODE_FILES]))


def jobs_todo(cur, kinds=None, sets=None, stages=None):
    """[(kind, set, patch)] of the runs of `stages` that have no result in `cur`, in the declared order (kinds, sets: optional filters); an 'adq' run is not needed where the rule keeps the isotropic kernel."""
    out = []
    for ks, patches in (STAGES if stages is None else stages):
        for name, pxy in patches:
            if sets is not None and name not in sets:
                continue
            pre = 'W|%s|x%d_y%d|' % (name, pxy[0], pxy[1])
            kind = cur.get('WH|%s|x%d_y%d|kind' % (name, pxy[0], pxy[1]))
            for kd in ks:
                if kinds is not None and kd not in kinds:
                    continue
                if kd.startswith('adq') and kind is not None and str(kind[0]) == 'isotropic':
                    continue
                if pre + JOB_TAG[kd] + '|secs' not in cur:
                    out.append((kd, name, pxy))
    return out


def part_synth(out_path=RAW, workers=1, kinds=None, sets=None, m2raw=M2RAW, psfraw=PSFRAW, repraw=REPEATRAW):
    """Item 2: the provenance, the copy of the real results, the worlds' headers, and every run of STAGES that has no result yet (the freeze is checked first; results are stored by this process only)."""
    cfg = run_cfg(m2raw)
    X.save_part(out_path, 'META', meta_arrays(cfg))
    part_real(out_path, m2raw, psfraw, repraw)
    cur = X.load_raw(out_path)
    for name, pxy in SYNTH:
        if 'WH|%s|x%d_y%d|kind' % (name, pxy[0], pxy[1]) not in cur:
            X.save_part(out_path, 'WH|%s|x%d_y%d' % (name, pxy[0], pxy[1]), world_header(name, pxy))
            print('W %s %s: header stored' % (name, pxy), flush=True)
    todo = jobs_todo(X.load_raw(out_path), kinds, sets)
    jobs = [(kd, name, pxy, cfg) for kd, name, pxy in todo]
    print('W: %d runs to do, in this order: %s' % (len(jobs), ', '.join('%s %s %s' % j[:3] for j in jobs)), flush=True)

    def store(res):
        kd, name, pxy, arrays = res
        X.save_part(out_path, 'W|%s|x%d_y%d|%s' % (name, pxy[0], pxy[1], JOB_TAG[kd]), arrays)
        print('W %s %s %s: stored (%.0f s)' % (kd, name, pxy, arrays['secs'][0]), flush=True)
    if workers > 1 and len(jobs) > 1:
        from multiprocessing import Pool
        with Pool(min(workers, len(jobs))) as pool:
            for res in pool.imap_unordered(world_job, jobs, chunksize=1):
                store(res)
    else:
        for j in jobs:
            store(world_job(j))


REAL_FRAME = {'iso|ours|k2': 'ours|k2', 'iso|default|k5': 'default|k5', 'cs|default': 'crowdsource (default)', 'cs|tuned': 'crowdsource (tuned)', 'phot|tuned': 'photutils (tuned)'}
REAL_REPEAT = dict({s: s for s in ('iso|ours|k2', 'iso|default|k5', 'adq|ours|k2', 'adq|default|k5')}, **{'cs|default': 'crowdsource (default)', 'cs|tuned': 'crowdsource (tuned)', 'phot|tuned': 'photutils (tuned)'})
RSD_RUNS = (('iso', 'ours'), ('iso', 'default'), ('adq', 'ours'), ('adq', 'default'))


def real_sources(m2raw=M2RAW, psfraw=PSFRAW, repraw=REPEATRAW):
    """The frozen results that the real columns are copied from: dict(m2raw, psf, rep) (rep: None if repeat_raw.npz does not exist)."""
    return dict(m2raw=X.load_raw(m2raw), psf=X.load_raw(psfraw), rep=X.load_raw(repraw) if os.path.exists(repraw) else None)


def real_kind(src, name, pxy):
    """The kernel the rule chose for a patch ('isotropic' or 'general') from the frozen results, or None."""
    if name in R.REPEAT:
        k = 'R|%s|x%d_y%d|kind' % (name, pxy[0], pxy[1])
        return None if src['rep'] is None or k not in src['rep'] else src['rep'][k]
    return src['psf'].get('S|%s|x%d_y%d|kind' % (name, pxy[0], pxy[1]))


def real_arrays(src, name, pxy):
    """The real-data results of a patch under the keys of the worlds (names relative to 'RF|<set>|x<x0>_y<y0>'): 'blk22' and 'tot' of every method of METHODS the frozen results hold, the residual
    scale ('<kernel>|<configuration>|rsd', one value per ISS run), the 'bic' record of each ISS run ('<kernel>|<configuration>|bic': block of the BIC minimum, blocks of the path, 1 if the entry cap was
    reached: a stalled path) and the kernel the rule chose ('kind'). Where the rule keeps the isotropic kernel the 'adq' rows are not copied (read the 'iso' rows)."""
    kind = real_kind(src, name, pxy)
    if kind is None:
        return {}
    arrays = {'kind': kind}
    same = str(kind[0]) == 'isotropic'
    if name in R.REPEAT:
        raw, pre = src['rep'], 'R|%s|x%d_y%d|' % (name, pxy[0], pxy[1])
        for s, k in REAL_REPEAT.items():
            if (same and s.startswith('adq')) or pre + k + '|blk22' not in raw:
                continue
            arrays[s + '|blk22'], arrays[s + '|tot'] = raw[pre + k + '|blk22'], raw[pre + k + '|tot']
        for kern, lab in RSD_RUNS:
            for what in ('rsd', 'bic'):
                if (same and kern == 'adq') or pre + '%s|%s|%s' % (kern, lab, what) not in raw:
                    continue
                arrays['%s|%s|%s' % (kern, lab, what)] = raw[pre + '%s|%s|%s' % (kern, lab, what)]
        return arrays
    raw, psf = src['m2raw'], src['psf']
    pre = ('T|x%d_y%d|' % pxy) if name == 'm2' else ('T2|%s|x%d_y%d|' % ((name,) + pxy))
    for s, k in REAL_FRAME.items():
        arrays[s + '|blk22'], arrays[s + '|tot'] = raw[pre + k + '|blk22'], raw[pre + k + '|tot']
    for lab in ('ours', 'default'):
        arrays['iso|%s|bic' % lab] = raw[pre + lab + '|bic']
    if name == 'm2':                                                    # the rule keeps the isotropic kernel on every M2 stamp: the residual scale is that of `exp_astro_psf --resid`
        for lab in ('ours', 'default'):
            arrays['iso|%s|rsd' % lab] = psf['R|m2|x%d_y%d|%s|rsd' % (pxy[0], pxy[1], lab)][:1]
        return arrays
    pp = 'P|%s|x%d_y%d|' % (name, pxy[0], pxy[1])
    for s, k in (('adq|ours|k2', 'ours|k2'), ('adq|default|k5', 'default|k5')):
        arrays[s + '|blk22'], arrays[s + '|tot'] = psf[pp + k + '|blk22'], psf[pp + k + '|tot']
    for lab in ('ours', 'default'):
        arrays['iso|%s|rsd' % lab], arrays['adq|%s|rsd' % lab] = psf[pp + lab + '|rsd'][:1], psf[pp + lab + '|rsd'][1:]
        arrays['adq|%s|bic' % lab] = psf[pp + lab + '|bic']
    return arrays


def part_real(out_path=RAW, m2raw=M2RAW, psfraw=PSFRAW, repraw=REPEATRAW):
    """The real-data results of the patches of SYNTH, copied from the frozen raws into the raw of this script (so that the report needs this raw alone), with the hashes of the sources."""
    src = real_sources(m2raw, psfraw, repraw)
    for name, pxy in SYNTH:
        arrays = real_arrays(src, name, pxy)
        if arrays:
            X.save_part(out_path, 'RF|%s|x%d_y%d' % (name, pxy[0], pxy[1]), arrays)
    X.save_part(out_path, 'RFsrc', dict(files=np.array(['%s:%s' % (os.path.basename(p), sha256_file(p)) for p in (m2raw, psfraw, repraw) if os.path.exists(p)])))


# ----------------------------------------------------------------------------------------------------------------------------------
# item 3: the catalogues
# ----------------------------------------------------------------------------------------------------------------------------------
BRIGHT_SETS = dict({'m2': [tuple(p) for p in X.TESTS], 'm53': [tuple(p) for p in X.SECOND['m53']['patches']], 'm92': [tuple(p) for p in X.SECOND['m92']['patches']]},
                   **{n: [tuple(p) for p in R.REPEAT[n]['patches']] for n in REPEAT_SETS})


def is_stalled(bic):
    """True if an ISS path reached the entry cap (`exp_astro_m2.MAX_ENTRIES`): the number of blocks of its 'bic' record (index of the BIC block, blocks of the path, flag) is the cap + 1 and the flag is
    set. Such a path has stalled: `imaging.sliding_path` can enter and exit one lattice point again and again at the same entry time until the cap, so that the run stops after a few
    hundred distinct entries whatever its stopping rule says. It is a defective run, not a censored valid one: its rows are marked STALLED and are never pooled."""
    return bic is not None and (int(bic[2]) == 1 or int(bic[1]) >= X.MAX_ENTRIES + 1)


def path_tag(s):
    """The ISS run a method suffix belongs to ('iso|ours' of 'iso|ours|k2'); None for the baselines, which are no paths."""
    return '|'.join(s.split('|')[:2]) if s.startswith(('iso', 'adq')) else None


def real_bic(src, name, pxy, s):
    """The 'bic' record of the frozen ISS run behind the method suffix `s` of a patch, or None (a baseline, or not held)."""
    tag = path_tag(s)
    if tag is None:
        return None
    kern, lab = tag.split('|')
    if name in R.REPEAT:
        return None if src['rep'] is None else src['rep'].get('R|%s|x%d_y%d|%s|bic' % (name, pxy[0], pxy[1], tag))
    if kern == 'iso':
        pre = ('T|x%d_y%d|' % pxy) if name == 'm2' else ('T2|%s|x%d_y%d|' % ((name,) + pxy))
        return src['m2raw'].get(pre + lab + '|bic')
    return src['psf'].get('P|%s|x%d_y%d|%s|bic' % (name, pxy[0], pxy[1], lab))


def alias_suffix(kind, s):
    """The key suffix to read for method suffix `s` of a patch whose kernel rule chose `kind`: the 'adq' runs are the 'iso' runs where the rule keeps the isotropic kernel."""
    return s.replace('adq', 'iso', 1) if kind == 'isotropic' and s.startswith('adq') else s


def real_dets(src, bases, name, pxy):
    """{method suffix: (x, y, flux)} of the real catalogues of a patch (the ISS at its reported readout, crowdsource, photutils tuned); only those the frozen results hold."""
    out = {}
    if name in R.REPEAT:
        raw, pre = src['rep'], 'R|%s|x%d_y%d|' % (name, pxy[0], pxy[1])
        if raw is None or bases['repeat'] is None or pre + 'kind' not in raw:
            return out
        for s in ('iso|ours', 'iso|default', 'adq|ours', 'adq|default'):
            if pre + s + '|det' in raw:
                out[s.replace('ours', 'ours|k2').replace('default', 'default|k5')] = tuple(raw[pre + s + '|det'].T)
        base = bases['repeat']
    else:
        raw = src['m2raw']
        pre = ('T|x%d_y%d|' % pxy) if name == 'm2' else ('T2|%s|x%d_y%d|' % ((name,) + pxy))
        out['iso|ours|k2'], out['iso|default|k5'] = tuple(raw[pre + 'ours|det'].T), tuple(raw[pre + 'default|det'].T)
        if name != 'm2':
            for s, k in (('adq|ours|k2', 'ours|det'), ('adq|default|k5', 'default|det')):
                out[s] = tuple(src['psf']['P|%s|x%d_y%d|%s' % (name, pxy[0], pxy[1], k)].T)
        base = bases['m2'] if name == 'm2' else bases['second']
    out['cs|default'], out['cs|tuned'] = X.cs_dets(base, pxy, (0.0, 0.0, 0.0, np.inf)), X.cs_dets(base, pxy, cs_filter(bases['tuned']['crowdsource']))
    out['phot|tuned'] = tuple(base['x%d_y%d|%s|%s' % (pxy[0], pxy[1], bases['tuned']['photutils'], c)] for c in ('x', 'y', 'flux'))
    return out


def load_bases(src):
    """The baseline catalogue files of the frozen results (M2, the second clusters, the repeat patches: None if absent) and the baselines' tuned choices."""
    return dict(m2=np.load(X.BASE), second=np.load(X.BASE2), repeat=np.load(BASE_REPEAT) if os.path.exists(BASE_REPEAT) else None, tuned=X.choose_baselines(src['m2raw']))


def part_bright(out_path=RAW, m2raw=M2RAW, psfraw=PSFRAW, repraw=REPEATRAW, bases=None):
    """Item 3 for every patch of BRIGHT_SETS: the real catalogues and, where the world has been run, the synthetic ones; stored under 'B|real|<set>|x<x0>_y<y0>' and 'B|synth|...'."""
    src = real_sources(m2raw, psfraw, repraw)
    cur = X.load_raw(out_path)
    bases = load_bases(src) if bases is None else bases
    for name, plist in BRIGHT_SETS.items():
        truth = truth_of(name)
        for pxy in plist:
            kind = real_kind(src, name, pxy)
            if kind is None:
                continue
            universe, box = sources_near(truth, pxy[0], pxy[1], MARGIN), sources_near(truth, pxy[0], pxy[1], 0)
            tag = '%s|x%d_y%d' % (name, pxy[0], pxy[1])
            arrays = {'kind': kind}
            for s, det in real_dets(src, bases, name, pxy).items():
                arrays.update({s + '|' + k: v for k, v in bright_arrays(det, universe, box).items()})
                arrays[s + '|stalled'] = np.array([int(is_stalled(real_bic(src, name, pxy, s)))])
            X.save_part(out_path, 'B|real|' + tag, arrays)
            sarr = {}
            for s in METHODS.values():
                key = 'W|%s|%s|det' % (tag, alias_suffix(str(kind[0]), s))
                if key in cur:
                    sarr.update({s + '|' + k: v for k, v in bright_arrays(tuple(cur[key].T), universe, box).items()})
                    run = path_tag(alias_suffix(str(kind[0]), s))
                    sarr[s + '|stalled'] = np.array([int(run is not None and is_stalled(cur.get('W|%s|%s|bic' % (tag, run))))])
            if sarr:
                X.save_part(out_path, 'B|synth|' + tag, sarr)
            print('B %s: %d real and %d synthetic catalogues' % (tag, sum(k.endswith('|nsig') for k in arrays), sum(k.endswith('|nsig') for k in sarr)), flush=True)


# ----------------------------------------------------------------------------------------------------------------------------------
# report
# ----------------------------------------------------------------------------------------------------------------------------------
def pc(v, nd=1):
    return 'nan' if not np.isfinite(v) else '%+.*f' % (nd, 100.0 * v)


def pci(v, ci, nd=1):
    return '%s [%s, %s]' % (pc(v, nd), pc(ci[0], nd), pc(ci[1], nd))


def pcu(v, ci, nd=0):
    """A fraction and its interval as unsigned percentages."""
    return '%.*f [%.*f, %.*f]' % (nd, 100.0 * v, nd, 100.0 * ci[0], nd, 100.0 * ci[1])


def num(v, ci, nd=2, sign=True):
    f = '%+.*f' if sign else '%.*f'
    return (f + ' [' + f + ', ' + f + ']') % (nd, v, nd, ci[0], nd, ci[1])


def kind_of(raw, name, pxy):
    """The kernel the rule chose for a patch ('isotropic' / 'general'), from the world's header, the copy of the real results or the bright-star arrays; '' if unknown."""
    for k in ('WH|%s|x%d_y%d|kind' % (name, pxy[0], pxy[1]), 'RF|%s|x%d_y%d|kind' % (name, pxy[0], pxy[1]), 'B|real|%s|x%d_y%d|kind' % (name, pxy[0], pxy[1])):
        if k in raw:
            return str(raw[k][0])
    return ''


def stamp_section(raw):
    L = ['## 1. The stamp against the stars (label-using, post hoc)', '',
         'Stars: truth sources with F606W < %g, more than %d px from the patch edge, with no other truth source within %g px that is brighter than m + %g. Each star: 25 x 25 cut-out around its nearest pixel, '
         'minus the median beyond %g px, moved by a cubic spline so that its flux centroid (radius %g px, %d iterations) is the middle pixel, cropped to 19 x 19, divided by its sum within %g px; the stack is '
         'the pixelwise MEDIAN. The stamp goes through the same centring and normalisation. The synthetic frame is the exact-stamp world of the patch (section 2: the same stars, a PSF that IS the stamp), so its '
         'stack shows what centring does to such a star. **Real minus synthetic** is the difference of the two relative differences of the central pixel (percentage points); "real / synthetic - 1" is the ratio of '
         'the central pixels of the two stacks minus 1 (the figure the referee quoted); the last column is the median over the stars of the ratio of the star\'s two central pixels (real over synthetic) minus 1, '
         'a paired estimator. FWHM: `astro_second.psf_fwhm` (px; the referee measured it in rings of 0.5 px, so his figures differ). Intervals: 95 %% percentile intervals of a bootstrap over the stars (%d '
         'resamples, the same stars for the real and the synthetic stack, stratified by patch when pooled; a pooled stack is compared with the stamps of its patches averaged with the numbers of stars as '
         'weights). A few stars (a companion just beyond %g px, or a sky level over-estimated by the crowd in the ring) have extreme central pixels in both frames; the median stack is robust to them, but '
         'the bootstrap interval of a median over a few tens of stars is wide and skewed by them, and the paired column is the sharper statement of the difference between the real and the synthetic stars.'
         % (STACK_MAG, HALF + PAD, STACK_RAD, STACK_DM, SKY_RING, CEN_RAD, CEN_ITER, NORM_RAD, X.BOOT, STACK_RAD), '',
         'Sets: ' + '; '.join('%s = %s' % (n, set_title(n)) for n in STAMP_SETS) + '.', '']
    hdr = ['set', 'patch', 'stars', 'FWHM stamp', 'FWHM real', 'FWHM synthetic', 'FWHM real - synthetic', 'centre real / stamp - 1 (%)', 'centre synthetic / stamp - 1 (%)',
           'centre real minus synthetic (points)', 'centre real / synthetic - 1 (%)', 'centre, median over the stars of real / synthetic - 1 (%)']
    rows, prow = [], []
    for name in STAMP_SETS:
        for pxy in set_patches(name) + ['all']:
            k = 'S|%s|%s|' % (name, 'all' if pxy == 'all' else 'x%d_y%d' % pxy)
            if k + 'est' not in raw:
                continue
            e, ci, n = raw[k + 'est'], raw[k + 'ci'], int(raw[k + 'n'][0])
            row = [name, 'all' if pxy == 'all' else '(%d, %d)' % pxy, '%d' % n, '%.2f' % fwhm(raw[k + 'ref']), num(e[4], ci[4], sign=False), num(e[5], ci[5], sign=False), num(e[6], ci[6]),
                   pci(e[0], ci[0]), pci(e[1], ci[1]), pci(e[2], ci[2]), pci(e[3], ci[3]), pci(e[8], ci[8])]
            (prow if pxy == 'all' else rows).append(row)
    L += ['### By patch\n', X.table(hdr, rows), '', '### Each set pooled\n', X.table(hdr, prow), '']
    return L


def blocks(raw, world, name, pxy, s):
    """blk22 (25, 4) of the method with key suffix `s` on a patch in the real ('RF') or synthetic ('W') results (the 'adq' rows of an isotropic-kernel patch are its 'iso' rows), or None."""
    k = '%s|%s|x%d_y%d|%s|blk22' % ('RF' if world == 'real' else 'W', name, pxy[0], pxy[1], alias_suffix(kind_of(raw, name, pxy), s))
    return raw.get(k)


def totals(raw, world, name, pxy, s):
    """The totals [true22, matched22, det22, det22_matched, true_all, matched_all, det_all] of the method with key suffix `s` on a patch in the real or synthetic results (alias as `blocks`), or None."""
    return raw.get('%s|%s|x%d_y%d|%s|tot' % ('RF' if world == 'real' else 'W', name, pxy[0], pxy[1], alias_suffix(kind_of(raw, name, pxy), s)))


STALLED = 'STALLED'


def stalled(raw, world, name, pxy, s):
    """True if the ISS path behind the method `s` of a patch in the real or synthetic results reached the entry cap (False for a baseline and where the record is absent)."""
    tag = path_tag(alias_suffix(kind_of(raw, name, pxy), s))
    return tag is not None and is_stalled(raw.get('%s|%s|x%d_y%d|%s|bic' % ('RF' if world == 'real' else 'W', name, pxy[0], pxy[1], tag)))


def split_pool(raw, patches, cols, get=None):
    """(valid, excluded) of the patches (set, origin) of a pool for the columns `cols` ([(world, suffix, header)]): `valid` have every column (`get`: `blocks` or `totals`) and no stalled path among them;
    `excluded` = [(patch, [headers of the columns whose path stalled])] of the patches that have every column but a stalled path (a stalled run is defective: it is left out of the pools, and the tables
    say so)."""
    get = blocks if get is None else get
    have = [q for q in patches if all(get(raw, w, q[0], q[1], s) is not None for w, s, _ in cols)]
    valid = [q for q in have if not any(stalled(raw, w, q[0], q[1], s) for w, s, _ in cols)]
    return valid, [(q, [h for w, s, h in cols if stalled(raw, w, q[0], q[1], s)]) for q in have if q not in valid]


def excluded_text(excluded):
    return '; '.join('%s (%d, %d): %s' % (q[0], q[1][0], q[1][1], ', '.join(hs)) for q, hs in excluded) if excluded else '-'


def pooled_blocks(raw, patches, world, s):
    parts = [blocks(raw, world, n, p, s) for n, p in patches]
    return None if not parts or any(b is None for b in parts) else np.concatenate(parts)


def cols_of(*specs):
    """[(world, method suffix, header)] from specs (world, method short name, header)."""
    return [(w, METHODS[m], h) for w, m, h in specs]


COLS_MAIN = cols_of(('real', 'iso', 'real iso'), ('real', 'adq', 'real adequate'), ('real', 'cs', 'real crowdsource'), ('synth', 'iso', 'synthetic iso'), ('synth', 'adq', 'synthetic adequate'),
                    ('synth', 'cs', 'synthetic crowdsource'))
COLS_DEFAULT = cols_of(('real', 'iso_def', 'real iso'), ('real', 'adq_def', 'real adequate'), ('real', 'cs', 'real crowdsource'), ('synth', 'iso_def', 'synthetic iso'),
                       ('synth', 'adq_def', 'synthetic adequate'), ('synth', 'cs', 'synthetic crowdsource'))
COLS_PHOT = cols_of(('real', 'iso', 'real iso'), ('real', 'phot', 'real photutils'), ('real', 'cs', 'real crowdsource'), ('synth', 'iso', 'synthetic iso'), ('synth', 'phot', 'synthetic photutils'),
                    ('synth', 'cs', 'synthetic crowdsource'))


def f_tables(raw, cols, title, note, need):
    """The F table by patch and the pooled F table of the columns `cols` ([(world, suffix, header)]): cells 'F [95 % interval]' ("(same)": the adequate kernel is the isotropic one; STALLED: the path reached
    the entry cap, a defective run); a pool counts the patches that have every column and no stalled path, the last column names the patches left out. Nothing is printed unless some world has a result for
    one of the method suffixes `need`."""
    present = [(n, p) for n, p in SYNTH if any(blocks(raw, 'synth', n, p, s) is not None for s in need)]
    if not present:
        return []
    rows = []
    idx = X.boot_idx(25)
    for n, p in present:
        kind = kind_of(raw, n, p)
        cells = []
        for w, s, _ in cols:
            b = blocks(raw, w, n, p, s)
            cells.append('-' if b is None else X.fmt_iv(*X.boot_cpf(b, idx)[2]) + (' (same)' if s.startswith('adq') and kind == 'isotropic' else '') + (' ' + STALLED if stalled(raw, w, n, p, s) else ''))
        rows.append([n, '(%d, %d)' % p, kind or '-'] + cells)
    L = ['### %s: by patch\n' % title, note, X.table(['set', 'patch', 'kernel chosen'] + [h for _, _, h in cols], rows), '']
    prow = []
    for label, patches in POOLS:
        full, excluded = split_pool(raw, patches, cols)
        if not full and not excluded:
            continue
        same = all(kind_of(raw, n, p) == 'isotropic' for n, p in full)
        idx = X.boot_idx_strat([25] * len(full)) if full else None
        cells = [X.fmt_iv(*X.boot_cpf(pooled_blocks(raw, full, w, s), idx)[2]) + (' (same)' if s.startswith('adq') and same else '') for w, s, _ in cols] if full else ['-'] * len(cols)
        prow.append([label, '%d of %d' % (len(full), len(patches))] + cells + [excluded_text(excluded)])
    if prow:
        L += ['### %s: pooled\n' % title, 'The blocks of the patches pooled (a pool counts the patches that have every column and no stalled path; the last column names the patches left out because a path stalled); the bootstrap '
              'resamples the blocks within each patch; the same resamples for every column.\n', X.table(['pool', 'patches'] + [h for _, _, h in cols] + ['left out (stalled path)'], prow), '']
    return L


def cp_table(raw, cols):
    """Completeness, purity (both at mag < 22) and the number of unmatched detections of any magnitude of the columns `cols`, pooled over the patches of each pool that have every column and no stalled path:
    'C / P / unmatched'; the last column names the patches left out."""
    rows = []
    for label, patches in POOLS:
        full, excluded = split_pool(raw, patches, cols, totals)
        if not full and not excluded:
            continue
        cells = []
        for w, s_, _ in cols:
            if not full:
                cells.append('-')
                continue
            tot = sum(totals(raw, w, n, p, s_) for n, p in full)
            c, pu, f = X.cpf(tot[:4])
            cells.append('%.3f / %.3f / %d' % (c, pu, tot[6] - tot[5]))
        rows.append([label, '%d of %d' % (len(full), len(patches))] + cells + [excluded_text(excluded)])
    return X.table(['pool', 'patches'] + [h for _, _, h in cols] + ['left out (stalled path)'], rows) if rows else ''


def lead_point(B, world, m):
    """F of the method m minus F of crowdsource in a world, on the pooled blocks B[(world, method)]."""
    return X.cpf(B[(world, m)].sum(0))[2] - X.cpf(B[(world, 'cs')].sum(0))[2]


def lead_draws(B, world, m, idx):
    """The same for each resample of the blocks (idx)."""
    return X.cpf(B[(world, m)][idx].sum(1))[2] - X.cpf(B[(world, 'cs')][idx].sum(1))[2]


COMPARISONS = (('synthetic iso - real iso', ('synth', 'iso'), ('real', 'iso')), ('synthetic adequate - real adequate', ('synth', 'adq'), ('real', 'adq')),
               ('synthetic crowdsource - real crowdsource', ('synth', 'cs'), ('real', 'cs')), ('real iso - real crowdsource', ('real', 'iso'), ('real', 'cs')),
               ('synthetic iso - synthetic crowdsource', ('synth', 'iso'), ('synth', 'cs')), ('real adequate - real crowdsource', ('real', 'adq'), ('real', 'cs')),
               ('synthetic adequate - synthetic crowdsource', ('synth', 'adq'), ('synth', 'cs')))
LEADS = (('lead of iso over crowdsource, real minus synthetic', 'iso'), ('lead of adequate over crowdsource, real minus synthetic', 'adq'))


def words(v, ci, up, down, flat):
    """A verdict from the sign of an interval: `up` if it excludes 0 from above, `down` from below, else `flat`."""
    return up if ci[0] > 0 else (down if ci[1] < 0 else flat)


def difference_table(raw):
    """Paired differences of F on the pooled blocks (same resamples for every row of a pool): between the worlds, between methods within a world, and the change of the lead of the ISS over crowdsource;
    a pool counts the patches that have all six results and no stalled path (the last column says how many, and names those left out)."""
    rows = []
    for label, patches in POOLS:
        full, excluded = split_pool(raw, patches, COLS_MAIN)
        if not full and not excluded:
            continue
        pooled = '%d of %d' % (len(full), len(patches)) + ('; left out (stalled path): ' + excluded_text(excluded) if excluded else '')
        if not full:
            rows.append([label, '-', '-', '-', pooled])
            continue
        idx = X.boot_idx_strat([25] * len(full))
        allgen = all(kind_of(raw, n, p) == 'isotropic' for n, p in full)
        B = {(w, m): pooled_blocks(raw, full, w, METHODS[m]) for w in ('real', 'synth') for m in ('iso', 'adq', 'cs')}
        for name, a, b in COMPARISONS:
            if allgen and 'adequate' in name:
                continue
            dc, dp, df = X.boot_diff(B[a], B[b], idx)
            between = a[0] != b[0]
            rows.append([label, name, X.fmt_iv(df[0], df[1]), words(*df, 'synthetic higher', 'synthetic lower', 'no difference') if between else X.verdict(*df), pooled])
        for name, m in LEADS:
            if allgen and m == 'adq':
                continue
            v = lead_point(B, 'real', m) - lead_point(B, 'synth', m)
            ci = np.percentile(lead_draws(B, 'real', m, idx) - lead_draws(B, 'synth', m, idx), [2.5, 97.5])
            rows.append([label, name, X.fmt_iv(v, ci), words(v, ci, 'real lead larger', 'real lead smaller', 'no change'), pooled])
    return X.table(['pool', 'difference of F (first minus second, paired blocks)', 'dF [95 % interval]', 'verdict', 'patches pooled'], rows) if rows else ''


def rsd_table(raw):
    rows = []
    for n, p in SYNTH:
        kind = kind_of(raw, n, p)
        row = [n, '(%d, %d)' % p, kind or '-']
        have = False
        for lab in ('ours', 'default'):
            for kern in ('iso', 'adq'):
                cells = []
                for pre in ('RF', 'W'):
                    k = '%s|%s|x%d_y%d|%s|%s|rsd' % (pre, n, p[0], p[1], alias_suffix(kind, kern), lab)
                    cells.append(('%.2f' % raw[k][0] + (' ' + STALLED if stalled(raw, 'real' if pre == 'RF' else 'synth', n, p, '%s|%s' % (kern, lab)) else '')) if k in raw else '-')
                    have = have or k in raw
                row += cells
        if have:
            rows.append(row)
    hdr = ['set', 'patch', 'kernel chosen'] + ['%s %s, %s' % (w, k, lab) for lab in ('tuned (kappa 2)', 'a priori default (kappa 5)') for k in ('isotropic kernel', 'chosen kernel') for w in ('real', 'synthetic')]
    return X.table(hdr, rows) if rows else ''


def frames_table(raw):
    rows = []
    for n, p in SYNTH:
        k = 'WH|%s|x%d_y%d|' % (n, p[0], p[1])
        if k + 'stats' not in raw:
            continue
        st, c = raw[k + 'stats'], raw[k + 'counts']
        rows.append([n, '(%d, %d)' % p, '%d' % c[0], '%d' % c[1], '%d' % c[2], '%.1f / %.1f' % (st[0, 0], st[1, 0]), '%.0f - %.0f / %.0f - %.0f' % (st[0, 1], st[0, 2], st[1, 1], st[1, 2]),
                     '%.0f / %.0f' % (st[0, 3], st[1, 3]), '%.2f' % raw[k + 'fwhm'][0]])
    return X.table(['set', 'patch', 'sources rendered', 'truth sources in the box', 'scored truth < 22', 'median DN, real / synthetic', '5 to 95 percentile DN, real / synthetic', 'maximum DN, real / synthetic',
                    'FWHM of the stamp (px)'], rows) if rows else ''


def runs_table(raw):
    rows = []
    for n, p in SYNTH:
        for kd in ('iso', 'adq', 'iso_def', 'adq_def', 'cs', 'phot'):
            k = 'W|%s|x%d_y%d|%s|' % (n, p[0], p[1], JOB_TAG[kd])
            if k + 'secs' not in raw:
                continue
            s = raw[k + 'secs']
            bic = raw.get(k + 'bic')
            rows.append([n, '(%d, %d)' % p, METHOD_LABEL[kd], '%d' % s[1], '%.0f' % s[0],
                         '-' if bic is None else '%d of %d' % (bic[0], bic[1]), (STALLED + ': the entry cap (%d) was reached, a defective run' % X.MAX_ENTRIES) if is_stalled(bic) else '-'])
    return X.table(['set', 'patch', 'run', 'entrants (detections for the baselines)', 'seconds', 'BIC block', 'stalled path'], rows) if rows else ''


def stalled_section(raw):
    """The paths that reached the entry cap, synthetic and real: listed, or the statement that there are none."""
    L = ['### Defective runs: paths that reached the entry cap\n',
         'A path that reaches the entry cap (%d) has stalled: `imaging.sliding_path` can enter and exit one lattice point again and again at the same entry time until the cap (a limit cycle reported after '
         'these runs were declared), so that such a run stops after a few hundred distinct entries whatever its stopping rule says. It is a **defective run, not a censored valid one**: every row computed from it is '
         'marked %s, no table counts it as valid, and the pooled tables and differences leave it out and name it. `src/core/imaging.py` is hash-frozen and was not touched; the rows stay in the raw with their '
         'record (`bic`).\n' % (X.MAX_ENTRIES, STALLED)]
    ran = [(n, p, kd) for n, p in SYNTH for kd in ('iso', 'adq', 'iso_def', 'adq_def') if 'W|%s|x%d_y%d|%s|bic' % (n, p[0], p[1], JOB_TAG[kd]) in raw]
    bad = [(n, p, kd) for n, p, kd in ran if is_stalled(raw['W|%s|x%d_y%d|%s|bic' % (n, p[0], p[1], JOB_TAG[kd])])]
    L.append('Synthetic worlds: ' + ('none of the %d ISS path runs reached the cap.' % len(ran) if not bad else '%d of the %d ISS path runs reached the cap and are defective: %s.' % (
        len(bad), len(ran), '; '.join('%s (%d, %d) %s' % (n, p[0], p[1], METHOD_LABEL[kd]) for n, p, kd in bad))))
    rran = [(n, p, tag) for n, p in SYNTH for tag in ('iso|ours', 'iso|default', 'adq|ours', 'adq|default') if 'RF|%s|x%d_y%d|%s|bic' % (n, p[0], p[1], tag) in raw]
    rbad = [(n, p, tag) for n, p, tag in rran if is_stalled(raw['RF|%s|x%d_y%d|%s|bic' % (n, p[0], p[1], tag)])]
    L.append('Real data (the frozen results of the same patches): ' + ('none of the %d ISS path runs copied reached the cap.' % len(rran) if not rbad else '%d of the %d ISS path runs copied reached the cap and are '
             'defective: %s.' % (len(rbad), len(rran), '; '.join('%s (%d, %d) %s' % (n, p[0], p[1], METHOD_LABEL[TAG_METHOD[tag]]) for n, p, tag in rbad))))
    return L + ['']


def world_section(raw):
    L = ['## 2. Exact-stamp synthetic worlds', '',
         'The world of a patch: every truth source of F606W < %g in the box grown by %d px, rendered with the patch\'s own 51 x 51 stamp at its sub-pixel position (cubic spline), flux from the HUGS magnitude by '
         '`astro_common.mag_to_counts`, the pipeline sky of the patch, Poisson noise at the patch gain, Gaussian noise of the dark variance; one noise realisation per patch (`default_rng([%d, crc32(set), x0, y0])`). '
         'Scored against the injected truth (the HUGS sources of the box) by `astro_common.score`: F of completeness and purity at mag < 22. The ISS at the frozen configuration (sigma_b 3, r_ex 1.2, kappa 2) with '
         'the isotropic kernel (**iso**) and with the kernel that the stamp rule chooses (**adequate**; where the rule keeps the isotropic kernel the two are one program on one input and the adequate column shows '
         'the isotropic result), crowdsource at its default (**crowdsource**). The real columns are the frozen results of the same patches (`m2_raw.npz`, `psf_raw.npz`, `repeat_raw.npz`). Intervals: 95 %% percentile '
         'intervals of the block bootstrap (20 x 20 px blocks, %d resamples), conditional on the patch and on the one noise realisation of the world. Nothing was tuned on a world.' % (MAG_MAX, MARGIN, DIAG_SEED, X.BOOT), '']
    if not any(k.startswith('W|') for k in raw):
        return L + ['No world has been run.', '']
    L += stalled_section(raw)
    L += f_tables(raw, COLS_MAIN, 'F at mag < 22, tuned configuration (kappa 2)', 'Real: the frozen results. Synthetic: the exact-stamp world of the same patch. %s: the path reached the entry cap (a defective run, see above).\n' % STALLED, [METHODS['iso']])
    ct = cp_table(raw, COLS_MAIN)
    if ct:
        L += ['### Completeness / purity at mag < 22 / unmatched detections of any magnitude, tuned configuration (kappa 2): pooled\n', 'The totals of the patches pooled, no interval (F over the blocks is in the table above). Unmatched '
              'detections: detections of any magnitude without a truth source within 0.5 px and 0.5 mag; the purity at mag < 22 does not see the faint ones (referee finding 7).\n', ct, '']
    dt = difference_table(raw)
    if dt:
        L += ['### Paired differences of F on the pooled blocks (tuned configuration)\n', 'Differences between the worlds (synthetic minus real, same method), between methods within a world, and the change of the lead of the ISS over '
              'crowdsource between the real data and the world (real lead minus synthetic lead; negative: the real data cost the ISS more than they cost crowdsource). Verdict: win / loss if the interval excludes 0.\n', dt, '']
    L += f_tables(raw, COLS_DEFAULT, 'F at mag < 22, a priori default (sigma_b 5, r_ex 1.2, kappa 5)', '%s: the path reached the entry cap (a defective run, see above); its F is not a valid row.\n' % STALLED, [METHODS['iso_def']])
    L += f_tables(raw, COLS_PHOT, 'F at mag < 22, photutils (the variant tuned on the development patches) beside the ISS and crowdsource', '', [METHODS['phot']])
    rt = rsd_table(raw)
    if rt:
        L += ['### Residual scale (robust sd of the standardized residual; 1 where the noise model and the kernel are right)\n', rt, '']
    ft = frames_table(raw)
    if ft:
        L += ['### The frames\n', 'The synthetic frame against the real one: median, 5 and 95 percentile and maximum of the DN.\n', ft, '']
    rn = runs_table(raw)
    if rn:
        L += ['### Run records\n', 'Entrants, seconds of the path to its last entrant (one core, shared machine), the BIC block and whether the path stalled (reached the entry cap: a defective run).\n', rn, '']
    return L


def bright_stats(raw, world, name, plist, s):
    """The statistics of item 1 of the bright stars of a method pooled over the patches `plist` of a set in the real or synthetic results, or None if a patch lacks the catalogue: dict of (value, interval) and
    counts, and `excluded`, the patches whose ISS path stalled (a defective run: left out of the pool). If every patch stalled, `n` is None."""
    kind = [kind_of(raw, name, p) for p in plist]
    keys = ['B|%s|%s|x%d_y%d|%s|' % (world, name, p[0], p[1], alias_suffix(kd, s)) for p, kd in zip(plist, kind)]
    if not plist or any(k + 'nsig' not in raw for k in keys):
        return None
    bad = [bool(raw[k + 'stalled'][0]) if k + 'stalled' in raw else False for k in keys]
    excluded = [p for p, b in zip(plist, bad) if b]
    keep = [k for k, b in zip(keys, bad) if not b]
    if not keep:
        return dict(n=None, excluded=excluded)
    ns, fr, pr = [raw[k + 'nsig'] for k in keep], [raw[k + 'fratio'] for k in keep], [raw[k + 'pairs'] for k in keep]
    return dict(n=sum(len(v) for v in ns), pairs=sum(len(v) for v in pr), split=boot_fraction(ns, lambda v: v >= 2), none=boot_fraction(ns, lambda v: v == 0), fratio=boot_median(fr),
                dmag=boot_median([v[:, 2] for v in pr]), dx=boot_median([v[:, 0] for v in pr]), dy=boot_median([v[:, 1] for v in pr]), excluded=excluded)


def stats_vector(st):
    """The numbers of a result of `bright_stats` as one vector."""
    out = [st['n'], st['pairs']]
    for k in ('split', 'none', 'fratio', 'dmag', 'dx', 'dy'):
        out += [st[k][0], st[k][1][0], st[k][1][1]]
    return np.array(out, float)


def same_stats(a, b):
    """True if two results of `bright_stats` are equal in every number and in the patches left out (None: not equal)."""
    return a is not None and b is not None and a['n'] is not None and b['n'] is not None and a['excluded'] == b['excluded'] and bool(np.array_equal(stats_vector(a), stats_vector(b), equal_nan=True))


def bright_section(raw):
    L = ['## 3. Splitting and flux bias of bright isolated stars (label-using, post hoc)', '',
         'A bright isolated star: a truth source in the scored area with F606W < %g and no other truth source within %g px that is brighter than m_star + %g. Its units: the detections of positive flux within %g px '
         'of its truth position; significant: at least %g %% of the star\'s flux. **Split**: at least 2 significant units. **No unit**: none. The flux ratio is the units\' summed flux over the star\'s. The bright pairs: '
         'one-to-one, position-only matching (%g px) of the detections in the scored area to the scored truth, truth F606W < %g; the median detection-minus-truth magnitude (dmag) and position (dx, dy, px). '
         'Intervals: 95 %% percentile intervals of a bootstrap over the stars (pairs) resampled within each patch. "Real, all patches": every test patch of the set; "real, the patches of the worlds": the patches that have '
         'a synthetic world (the same stars); "synthetic": the worlds. A catalogue of the adequate kernel on a patch where the rule keeps the isotropic kernel is the isotropic catalogue (a set whose patches all keep it has no adequate rows); '
         'the rows of the frozen crowdsource post-filter are printed only where they differ from those of the default. The catalogue of an ISS path that stalled (reached the entry cap: a defective run, section 2) is '
         'left out of the pool and the last column names the patch; a row counts the patches whose path did not stall, so the real and the synthetic row of a catalogue rest on the same patches only where '
         'nothing was left out (the numbers of stars show it).'
         % (BRIGHT_MAG, BRIGHT_RAD, BRIGHT_DM, UNIT_RAD, 100 * UNIT_FRAC, PAIR_RAD, PAIR_MAG), '']
    for name, plist in BRIGHT_SETS.items():
        wl = [p for p in plist if (name, p) in SYNTH and any(k.startswith('B|synth|%s|x%d_y%d|' % (name, p[0], p[1])) for k in raw)]
        iso_only = all(kind_of(raw, name, p) == 'isotropic' for p in plist)
        rows = []
        for m in BRIGHT_ORDER:
            s = METHODS[m]
            if iso_only and m.startswith('adq'):
                continue                                                  # the rule keeps the isotropic kernel on every patch of the set: the rows would repeat the isotropic ones
            worlds = [('real', plist, 'real, all %d patches' % len(plist))] if len(wl) != len(plist) else []
            worlds += [('real', wl, 'real, the patches of the worlds'), ('synth', wl, 'synthetic')]
            for world, pl, lab in worlds:
                st = bright_stats(raw, world, name, pl, s) if pl else None
                if st is None or (m == 'cs_tuned' and same_stats(st, bright_stats(raw, world, name, pl, METHODS['cs']))):
                    continue                                              # the frozen post-filter is printed only where it differs from the default
                ex = '; '.join('(%d, %d)' % p for p in st['excluded']) if st['excluded'] else '-'
                if st['n'] is None:
                    rows.append([METHOD_LABEL[m], lab, '%s: every path stalled' % STALLED] + ['-'] * 7 + [ex])
                    continue
                rows.append([METHOD_LABEL[m], lab, '%d' % st['n'], pcu(st['split'][0], st['split'][1]), pcu(st['none'][0], st['none'][1]), num(st['fratio'][0], st['fratio'][1], 2, False), '%d' % st['pairs'],
                             num(st['dmag'][0], st['dmag'][1]), num(st['dx'][0], st['dx'][1]), num(st['dy'][0], st['dy'][1]), ex])
        if rows:
            L += ['### %s: %s\n' % (name, set_title(name)), X.table(['catalogue', 'sample', 'bright isolated stars', 'split (%)', 'no unit (%)', 'flux ratio, median', 'bright pairs', 'dmag, median', 'dx, median', 'dy, median', 'patches left out (stalled path)'], rows), '']
    return L


LIMITATIONS = [
    'One noise realisation per world; the intervals are bootstrap intervals over 20 x 20 px blocks (stars) of that realisation and of the patch, conditional on both. No world was repeated, so a difference that is small against the '
    'intervals is not resolved, and the realisation-to-realisation scatter of a world is not measured here (the referee\'s frames of four of these worlds, `review/REVIEW_roundIS1b_A.md`, are another realisation of them; the '
    'F of the same method on the same patch differs between the two realisations).',
    'The worlds contain only what is put into them. The stamp is exact (no PSF variation over the patch, none beyond 25 px; the stamp is the psField reconstruction at the patch centre, clipped), the noise is exactly Poisson '
    'plus Gaussian at the patch gain and dark variance (the real noise has correlations, cosmic rays, flat-field and sky errors), the truth is the HUGS list (sources not in it, and the stars saturated in F606W, are absent from '
    'both the world and the scoring), positions and magnitudes are exact (no proper motion or astrometric offset: the image-based truth-to-frame shifts of `second_cluster_spec.md` are not applied; no colour term; '
    'the frame\'s calibration constant differs from the benchmark\'s, `repeat_spec.md`, which the real data carry and the worlds do not), and nothing saturates: a bright truth star has an unclipped profile in the world '
    'and a clipped one in the real frame (maximum DN in the table of the frames), and the world lacks the unresolved light and the sources the HUGS list does not hold (its median DN is lower than the frame\'s).',
    'The ISS in a world is given the stamp and fits its own kernel to it as in the real runs (`imaging.fit_psf_parts` or the adequacy rule), so a world tests the kernel misfit; it cannot test a stamp that differs from the stars. '
    'The baselines get the exact stamp as in the real runs.',
    'The stamp stack uses stars brighter than F606W = 20 (a few tens per patch, the table gives the number): it measures the PSF core at high signal to noise, not its wings or its dependence on magnitude or colour, and the interval of a patch is wide. '
    'Centring on a noisy centroid and the median stack are estimators with their own bias; the synthetic stack cancels it only to the extent that a rendered star and a real one have the same signal to noise and the same neighbours '
    '(the same truth stars are used for both).',
    'The definitions of section 3 (20.5 mag, 3 px, 2.5 mag, 10 % of the flux, 2 px, 19.5 mag, 0.7 px) are the referee\'s and were not tuned; they use the truth. A split says that more than one significant unit lies within 2 px of '
    'a bright star; it does not say that the units are wrong (a star with a lumpy core may be legitimately fitted by two units): the flux ratio and the magnitude offsets measure the consequence. Pairs are found by position '
    'alone, so a split star contributes the unit nearest to it.',
    'A path that reaches the entry cap has stalled (section 2): its rows are defective, marked STALLED, left out of every pool and named; the paths of the frozen results that stalled are among the real columns, '
    'the cause is in `imaging.sliding_path` (hash-frozen, untouched here), and a guarded path is a separate, labelled rerun.',
    'These controls do not identify a cause of the real M92 or M53 results beyond what the worlds reproduce; a quantity that a world does not reproduce is not thereby explained.']


def provenance(raw):
    L = []
    if 'META|cfg' in raw:
        c = raw['META|cfg']
        L.append('Configuration read from `freeze.json` (checked by `exp_astro_m2.check_freeze` before the runs): sigma_b %g, r_ex %g, kappa %g; the a priori default %s; tuned baselines: %s. Code the runs used (SHA-256, first 12 digits): %s.' % (
            c[0], c[1], c[2], X.DEFAULT, ', '.join(str(v) for v in raw['META|tuned']), '; '.join('%s %s' % (f.split(':')[0], f.split(':')[1][:12]) for f in raw['META|files'])))
    if 'RFsrc|files' in raw:
        L.append('Frozen results the real columns were copied from (SHA-256, first 12 digits): %s.' % '; '.join('%s %s' % (f.split(':')[0], f.split(':')[1][:12]) for f in raw['RFsrc|files']))
    return L


def report(raw):
    """The markdown report from the raw alone."""
    L = ['# Controls for the astronomy results: the stamp against the stars, exact-stamp synthetic worlds, bright-star splitting and flux bias', '',
         'Generated by `experiments/exp_astro_diag.py` from `diag_raw.npz`; every number below is printed by the script. **POST HOC DIAGNOSTICS, NOT TESTS**: they follow the frozen transfer test and the kernel study, '
         'they use the truth list (HUGS) and the PSF stamps, and nothing in the methods or their configuration changed on their strength. They port the controls of referee round IS1b-A '
         '(`review/REVIEW_roundIS1b_A.md`, findings 1, 2 and 8); the synthetic frames are rendered by the script and not stored.', ''] + [p + '\n' for p in provenance(raw)]
    L += stamp_section(raw)
    L += world_section(raw)
    L += bright_section(raw)
    L += ['## Limitations\n'] + ['* ' + t for t in LIMITATIONS] + ['']
    return '\n'.join(L) + '\n'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--stamps', action='store_true')
    ap.add_argument('--synth', action='store_true')
    ap.add_argument('--real', action='store_true', help='copy the real-data results of the patches of the worlds into the raw (done by --synth too)')
    ap.add_argument('--bright', action='store_true')
    ap.add_argument('--from-raw', action='store_true')
    ap.add_argument('--only', help='comma-separated kinds of runs of --synth: iso, adq, cs, phot, iso_def, adq_def')
    ap.add_argument('--sets', help='comma-separated sets of --synth: m2, m53, m92, m2r1, m2r2, m2r3')
    ap.add_argument('--workers', type=int, default=int(os.environ.get('ASTRO_WORKERS', '1')))
    A_ = ap.parse_args()
    if A_.stamps:
        part_stamps()
    elif A_.synth:
        part_synth(workers=A_.workers, kinds=A_.only.split(',') if A_.only else None, sets=A_.sets.split(',') if A_.sets else None)
    elif A_.real:
        part_real()
    elif A_.bright:
        part_bright()
    elif A_.from_raw:
        with open(MD, 'w') as f:
            f.write(report(X.load_raw(RAW)))
        print('wrote', MD)
    else:
        ap.error('give --stamps, --synth, --bright or --from-raw')


if __name__ == '__main__':
    main()
