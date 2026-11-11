"""
Post hoc sensitivity analyses of the held-out repeat-frame test: the controls of referee round REP, ported into the repository
===============================================================================================================================

STATUS: POST HOC. The repeat-frame test (`exp_astro_repeat.py`, `repeat_plan.md`, `repeat_raw.npz`, `repeat.md`, `freeze_repeat.json`) is finished and is not touched. Everything below is
done AFTER its result was seen, by an adversarial referee first (`review/REVIEW_roundREP.md`, findings 1, 3, 4, 5, 6, 9, 10) and here as committed, tested code. None of it changes `repeat.md`, the
plan or a freeze, and nothing in the project's methods or configuration (the ISS, its kernel rule, its configuration, the baselines and their tuning, the scoring of `astro_common.py`) may
change on its strength. The analyses move the TRUTH, the matching tolerances, the depth cut, the magnitudes of the detections and the block size and re-count; they never re-run a path.
This docstring was written to disk before any analysis below was run; a later substantive change of it is listed in the report's section 11.

    python experiments/exp_astro_repeat_sens.py --run [--cyclefix PATH]   re-scores the stored detections, measures the star offsets, writes results/astro/repeat_sens_raw.npz (needs the data of
                                                               $ASTRO_DATA, $ASTRO_DATA_2, $ASTRO_DATA_3; under a minute on one core; run it with `python -I`); the option is the hook of D8
    python experiments/exp_astro_repeat_sens.py --from-raw     regenerates results/astro/repeat_sens.md from repeat_sens_raw.npz alone (output directory: $REPEAT_SENS_OUT; about 15 s)
    Environment (all optional): REPEAT_SENS_OUT (directory of the raw and the report), REPEAT_SENS_REPEAT_RAW, REPEAT_SENS_M2_RAW, REPEAT_SENS_PSF_RAW, REPEAT_SENS_BASE (the read-only inputs).

D0. Inputs and the re-scoring rule (used by everything below)
-------------------------------------------------------------
Stored detections (x, y in patch pixels, flux in DN): `R|<frame>|x<x0>_y<y0>|<kern>|<cfg>|det` of `repeat_raw.npz`, kern in {iso, adq, scl} (isotropic kernel; the kernel the rule chose, which is
the isotropic one on four of the six patches; the chosen kernel with the lengths scaled by the seeing), cfg in {ours, default} (the tuned configuration, reported at kappa 2; the a priori default,
reported at kappa 5). The baselines: `baseline_catalogs_repeat.npz`, crowdsource through `exp_astro_m2.cs_dets` with the filter whose name is stored in `crowdsource (tuned)|name` (the default
is the filter (0, 0, 0, inf)), photutils with the tag stored in `photutils (tuned)|name` (default: tag `photutils`). Truth of a patch: `astro_common.truth_for_patch` of `astro_repeat.truth_arrays`.
The BIC readouts (key `k0`) have no stored detections: they are read from their stored counts and blocks (offset 0, tolerances 0.5 px / 0.5 mag, depth 22, blocks of 20 px) and are never re-scored.
`rescore(det, truth_patch, dpos, dmag, shift, moff, cuts, cells)` is our own scorer. It follows `astro_common.score` and differs from it in three respects only: the tolerances are arguments,
the matching is solved per connected component of the admissible pairs (a sparse formulation of the same problem), and it can count at several depth cuts and block sizes. Definitions:
  * scored detections: finite magnitude 22.5 - 2.5 log10(0.00546689 flux) minus `moff`, at least EDGE = 2.5 px from the patch edge (x, y in [2.0, 97.0]); scored truth: the `inner` truth of the patch;
  * admissible pair: distance <= dpos AND |detection magnitude - truth magnitude| <= dmag (both inclusive); the matching is one to one, of maximum size and, among those, of minimum total distance
    (maximum weight with weight BIG - distance, BIG = 1e4), exactly the objective of `astro_common.match`; the depth cut does not enter the matching;
  * counts at a cut c (c in 20, 21, 22): true_c = scored truth with magnitude < c; matched_c = matched truth with truth magnitude < c; det_c = scored detections with (corrected) magnitude < c;
    detm_c = matched detections with magnitude < c; and without a cut true_all, matched_all, det_all: the vector [true_c, matched_c, det_c, detm_c, true_all, matched_all, det_all] of
    `exp_astro_m2._pack`; F, completeness and purity from it (`exp_astro_m2.cpf`);
  * blocks of `cell` px (cell = 20: the report's 25 blocks per patch): the counts per block, a truth source and its match in the block of the truth source, a detection and its match in the block of
    the detection (`astro_common.score`); the block index is clip(floor((position + 0.5) / cell), 0, 100 / cell - 1), row-major;
  * a SHIFT of the truth by (sx, sy) moves the positions the matching sees; the scored set, the denominators and the block of every truth source are those of the UNSHIFTED truth, and the
    detections are never moved (so the denominators are the same for every shift and the comparison across shifts is paired). The offset of a method is the median over its matched pairs with truth
    magnitude < 21 of (detection - moved truth), as the `offset` of `astro_common.score`;
  * a magnitude offset `moff` (per patch and method) is estimated in a first pass at the same shift and tolerances without correction: the median of (detection magnitude - truth magnitude) of the
    matched pairs with 17 < truth magnitude < 20, if there are at least 5 such pairs, else 0; the second pass subtracts it from the detections' magnitudes (matching, depth cut and counts).
At the default tolerance, no shift, no offset and cut 22 the scorer must reproduce the stored `tot` and `blk22` of every stored row (60 rows: ISS iso, adq, scl at both configurations and the
four baseline rows, on six patches); `--run` refuses to write a raw unless it does.

D1. Truth-to-image offsets (the input of D2 and D4(iii))
---------------------------------------------------------
Two sources, each an offset (ox, oy) in px that is ADDED to the truth positions (the sign of the spec's "shift to add to the truth positions"):
  * scan: the rigid-scan offsets of `repeat_spec.md` section 6.2 read from the machine-readable copy `repeat_validation.json` ($ASTRO_DATA_3): key `<frame>["<x0>,<y0>"]["masked"]["asTrans"]["dx"|"dy"]`
    (the "masked scan, asTrans" column of the spec table; the truth lists of the repeat frames are asTrans-mapped); per patch. For the M2 windows the key `["masked"]["header"]` (the M2 truth file is
    header-mapped); the plain scan (`["image_shift_px"]`) is listed for reference.
  * own: template fits of bright isolated stars with the patch's own PSF stamp (`psf51`, unit sum). Star selection, declared before the fits: truth rows of the patch with 15.5 <= F606W < 19.5, not
    HST-saturated (nf606 != 99), at least 15 px from every patch edge (15 <= x, y <= 84 in patch pixels), with no OTHER truth row (any magnitude, saturated ones included) within 8 px that is
    brighter than m + 2. Fit: the 11 x 11 pixel cut around the rounded truth position, model A S(x - cx, y - cy) + B (S the stamp rendered by cubic-spline interpolation, origin its middle pixel;
    B a constant), weights 1 / sqrt(max(DN, 1) / gain + dark variance), `scipy.optimize.least_squares`, start (truth position, flux of the cut above its median, median); the star's offset is
    (cx, cy) - truth position; a fit that did not converge or with |offset| > 1 px is discarded (counted). Per patch: n, the median of the offsets, the robust scatter 1.4826 MAD. Per frame: the median
    over the stars of all the frame's patches with a 95 % percentile interval of a bootstrap over the stars (2000 resamples, seed 0). The offset used by every `own` scenario is the FRAME
    median, applied to each patch of the frame (a patch has 3 to 14 stars; the offset is a property of the frame's truth-to-image mapping); the per-patch medians are reported and used in the one
    scenario `ownp` (a patch with fewer than 3 stars uses the frame median there). A frame with fewer than 3 stars has no `own` rows. The same fits are made on the M2 test windows and the
    M53 and M92 patches and reported (no re-scoring of those).
    AMENDMENT (written after the NUMBER of stars of the selection above, here selection A, was seen: 3, 2 and 3 stars on m2r1, m2r2 and m2r3, 1 to 2 per patch, too few for a median; no offset
    of any other selection had been seen): selection B keeps the fit and changes the geometry to that of the referee's `template_fit.py`: 15.5 <= F606W < 19.5, not HST-saturated, at least 8 px from
    every patch edge (8 <= x, y <= 91), no other truth row within 5 px that is brighter than m + 2.5. The scenarios `own` and `ownp` use selection B; both selections are reported side by side. B is
    not a bias-free selection: it allows neighbours fainter than m + 2.5 at any distance and any neighbour at 5 px or more, whose pull on the fitted position a noise-free simulation with the stamp of m2r2
    (133, 281) gives by distance and flux ratio (printed in the report); the pull is toward the neighbour, in a direction that differs from star to star, which the median over the stars averages down.

D2. (a) P1 against the truth offset
-----------------------------------
The patches of P1 are those where the rule picks a general kernel (the two patches of m2r2). Readouts: tuned (kappa 2: `adq|ours|k2` against `iso|ours|k2`), a priori default (kappa 5) and the BIC stop
(`adq|default|k0` against `iso|default|k0`; alpha = 0 only, see D0). Scenarios: the truth shifted by alpha x offset, alpha in {0, 0.25, 0.5, 0.75, 1}, for the sources scan and own (alpha = 0 is
the stored test) and, at alpha = 1, `ownp` (per-patch own medians). For every scenario: the pooled dF (the two patches, stratified block bootstrap, 20 px blocks, 2000 resamples, seed 0: the report's
interval, `exp_astro_m2.boot_idx_strat`), dF per patch (2000, seed 0), F of both kernels, the median position offset of the bright matches of both kernels. For the three readouts at alpha = 0 and for the
tuned and default readouts at every alpha: the pooled interval from 40000 resamples (seed 0) and, over the seeds 1 to 200 (2000 resamples each), the share of seeds for which the lower end of the pooled
interval is above 0 (and its mean, minimum and maximum). Verdict words as the plan: win (lower end > 0), loss (upper end < 0), tie.

D3. (b) P1 robustness: block size
---------------------------------
Blocks of 10, 20, 25 and 50 px (100, 25, 16 and 4 blocks per patch), stratified bootstrap within patch (2000, seed 0): tuned and default readouts at alpha = 0 (the stored detections), pooled over the two
m2r2 patches, per patch, and over all six patches (the row is marked DILUTED: four of the six patches have dF = 0 by construction).

D4. (c) P2 robustness
---------------------
The pairs of the declared analysis: ISS with the chosen kernel (`adq|ours|k2`) and with the isotropic kernel (`iso|ours|k2`), each against tuned crowdsource and tuned photutils (four pairs), per patch
(verdict by the paired interval of one patch's blocks) and pooled over the six patches (stratified), under: (i) the match tolerances 0.5 px / 0.5 mag (the stored scoring), 0.35 / 0.5, 0.5 / 0.3, 0.75 / 0.75;
(ii) the depth cuts m < 20, 21, 22 (also F per method per patch and pooled; counts at the cut, matching unchanged); (iii) the truth shifted by alpha x offset, both sources, alpha in the grid of D2, every patch
by its own offset (scan: per patch; own: the frame median); (iv) the per-method median magnitude offset removed before matching (D0), without a shift and with the truth shifted by the full scan and the full own offsets.

D5. (d) Scaled lengths by frame
-------------------------------
The stored rows `scl|ours|k2` and `scl|default|k5`, `scl|default|k0` against the pixel-unit rows `adq|...` of the same configuration: F of both, dF (scaled - pixel) with the paired interval
(2000, seed 0), per frame (the frame's two patches pooled) and over the six patches, with the defective flag of D8 (the pixel-unit default rows of m2r3 are the capped runs).

D6. (e) F ratios and the seeing ladder
--------------------------------------
F(ISS, chosen kernel, tuned configuration) / F(tuned crowdsource), and the same against tuned photutils, per frame from `repeat_raw.npz` (the frame's two patches pooled) and for M2 (the four test windows of
`m2_raw.npz`, `T|x<x0>_y<y0>|ours|k2`), with the 95 % interval of the ratio over the same stratified resamples (2000, seed 0); the ratio of M2 on the windows that each repeat frame repeats ("the M2 windows'
own ratio"); and the seeing ladder as `repeat.md` prints it (F per window and frame, ISS / crowdsource / photutils), with the stamp FWHM of each frame.

D7. (f) Plan items missing from `repeat.md`
-------------------------------------------
  * the paired differences of the ISS default (kappa 5) and BIC readouts, with the chosen and with the isotropic kernel, against the baselines' DEFAULT photutils (and, for comparison with repeat.md, default
    crowdsource): per patch (verdict, dF), per frame, pooled over the six patches with and without the defective patches (D8);
  * the P3 verdict per readout: on the patches where the rule picks a general kernel, the robust standard deviation of the standardized residual (stored `rsd`, 1.4826 MAD) of the isotropic kernel is above
    that of the chosen kernel on every such patch ("supported") or not ("not supported", the plan's wording: the diagnosis K is not supported on these frames); readouts tuned and default (the BIC readout has no
    stored residual scale); beside the patchwise P1 dF;
  * the descriptive item "stars brighter than 19.5 mag matched": matched / in the truth for the first nine magnitude bins (15.0 to 19.5) of the stored rows, per patch and method.

D8. Defective (capped) paths
----------------------------
The a priori default runs of the isotropic kernel on m2r2 (133, 281), m2r3 (370, 441) and m2r3 (110, 501) hit the entry cap of 2500 entrants; this is a defective path of `imaging.sliding_path` (an enter/exit cycle of one
lattice point, `review/REVIEW_roundREP.md` finding 2), not a censored one. The flag of a run is the third element of the stored `bic = [block, hist length, cap flag]` of its kern and cfg (the chosen-kernel run
of a patch where the rule keeps the isotropic kernel is that same run, so it carries the same flag). Every analysis row that contains an ISS run whose flag is 1 is marked "defective path: yes"; a pooled row is
never given without the same row computed from the patches whose compared runs are not defective. No re-run is made here. Optional hook: `--run --cyclefix PATH` reads the guarded re-runs of `experiments/exp_astro_cyclefix.py` (`results/astro/cyclefix_raw.npz`)
from a raw that has them; the key format is the one that script's docstring and code give: `G|` + the stored key of the run + the relative key, for a repeat run `G|R|<frame>|x<x0>_y<y0>|<kern>|<cfg>|k<kappa>|tot`
(7 counts) and `|blk22` (25 x 4), `|bic` and `|secs`, kern iso, adq or scl (a run stored under both kernel names is there under iso only). The rows are copied into this raw (`gr|...`) and the report gains
section 8.1: the capped runs beside their guarded re-runs, and P1 and the default pairs at alpha 0 with the guarded rows in place of the defective ones. Without the option (the committed report) there is no such section.

D9. The camera-column confound
------------------------------
One table of the frames M2, m2r1, m2r2, m2r3, M53, M92: SDSS run-camcol-field, the camera column, gain and dark variance of the patches (`patch_arrays.npz`), the stamp FWHM, the rms of the isotropic stamp fit
(`psf_raw.npz` for M2, M53, M92; `repeat_raw.npz` for the repeat patches), the ellipticity of the stamp and the kernel the rule chooses.

D10. Positions of the two kernels on m2r2
-----------------------------------------
For the two m2r2 patches, the stored tuned catalogues of both kernels: the median (detection - truth) of their bright matches (stored `offset`), and the median (chosen - isotropic) position difference of
the flux-consistent unit pairs (mutual nearest neighbours within 0.4 px, flux ratio in [0.8, 1.25]) with detected magnitude < 19, 20, 21 and 22. It describes the 0.1 px difference; it does not identify its cause.

Intervals: percentile intervals (2.5, 97.5) of a block bootstrap, stratified within patch, `exp_astro_m2.boot_idx_strat`, 2000 resamples, seed 0, the same resamples for every row of a comparison; a single
patch is resampled alone. Both kernels, all methods and all shifts of a comparison see the same resamples. Nothing here is a test of a hypothesis declared before the data: no multiplicity correction is
made and the verdicts are descriptive.
"""
import argparse
import hashlib
import json
import os
import sys

os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
from scipy.ndimage import map_coordinates, spline_filter
from scipy.optimize import least_squares, linear_sum_assignment
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

HERE = os.path.abspath(os.path.dirname(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)
import astro_common as A                  # noqa: E402
import astro_repeat as R                  # noqa: E402
import astro_second as S2                 # noqa: E402
import exp_astro_m2 as X                  # noqa: E402
import exp_astro_repeat as E              # noqa: E402

OUT = os.environ.get('REPEAT_SENS_OUT', os.path.join(REPO, 'results', 'astro'))
RAW, MD = os.path.join(OUT, 'repeat_sens_raw.npz'), os.path.join(OUT, 'repeat_sens.md')
RES = os.path.join(REPO, 'results', 'astro')
REPEAT_RAW = os.environ.get('REPEAT_SENS_REPEAT_RAW', os.path.join(RES, 'repeat_raw.npz'))          # the finished test, read only
M2_RAW = os.environ.get('REPEAT_SENS_M2_RAW', os.path.join(RES, 'm2_raw.npz'))                      # the frozen M2 test, read only
PSF_RAW = os.environ.get('REPEAT_SENS_PSF_RAW', os.path.join(RES, 'psf_raw.npz'))                   # the stamp fits of every patch, read only
BASE = os.environ.get('REPEAT_SENS_BASE', os.path.join(RES, 'baseline_catalogs_repeat.npz'))        # the baselines' catalogues on the repeat patches, read only
VALID = os.path.join(R.ROOT3, 'repeat_validation.json')                                             # the machine-readable copy of repeat_spec.md section 6
PATCHES = list(E.PATCHES)                                                                           # [(frame, (x0, y0))] of the six repeat patches
FRAMES = sorted(R.REPEAT)
N = 100                                                                                             # side of a patch (px)
ALPHAS = (0.0, 0.25, 0.5, 0.75, 1.0)
CUTS = (20.0, 21.0, 22.0)                                                                           # depth cuts; the report's is the last
ICUT = CUTS.index(22.0)
CELLS = (10, 20, 25, 50)                                                                            # block sizes (px); the report's is A.CELL
TOLS = ((0.5, 0.5), (0.35, 0.5), (0.5, 0.3), (0.75, 0.75))                                          # (dpos px, dmag mag)
BIG = 1.0e4                                                                                         # weight of an admissible pair is BIG - distance
BOOT = X.BOOT                                                                                       # 2000 resamples, seed 0
BIG_B = 40000
SEEDS = tuple(range(1, 201))
MAG_LO, MAG_HI, MAG_MIN_PAIRS = 17.0, 20.0, 5                                                       # the magnitude-offset estimate
SELECTIONS = {'A': dict(mlo=15.5, mhi=19.5, margin=15.0, riso=8.0, dmag_iso=2.0),                  # the declared star selection (D1)
              'B': dict(mlo=15.5, mhi=19.5, margin=8.0, riso=5.0, dmag_iso=2.5)}                    # the amendment of D1, the one the own scenarios use
OWN_SEL = 'B'
STAR = dict(half=5, max_offset=1.0, min_stars=3)                                                    # the cut of the star fit (11 x 11 px), the largest offset kept, the fewest stars of a frame
BIAS_D, BIAS_Q = (3.0, 4.0, 5.0, 6.0, 7.0, 8.0), (0.063, 0.1, 0.16, 0.25, 0.4, 1.0)                # the distances (px) and flux ratios of the neighbour-bias simulation
ISS_KEYS = ('iso|ours|k2', 'adq|ours|k2', 'iso|default|k5', 'adq|default|k5')                       # the ISS rows re-scored in every scenario
TUNED_BASE = ('crowdsource (tuned)', 'photutils (tuned)')
SC_KEYS = ISS_KEYS + TUNED_BASE
BASE_EXTRA = ('scl|ours|k2', 'scl|default|k5', 'crowdsource (default)', 'photutils (default)')      # re-scored at the stored scoring only
STORED_ISS = tuple('%s|%s' % (kern, c) for kern in ('iso', 'adq', 'scl') for c in ('ours|k2', 'default|k5', 'default|k0'))
STORED_BASE = ('crowdsource (default)', 'crowdsource (tuned)', 'photutils (default)', 'photutils (tuned)')
STORED_ROWS = STORED_ISS + STORED_BASE
M2_ROWS = ('ours|k2', 'crowdsource (tuned)', 'photutils (tuned)')
CYCLEFIX_RAW = os.path.join(RES, 'cyclefix_raw.npz')          # the guarded re-runs of another script (optional, `--cyclefix`); read only


def _scenarios():
    out = {'base': dict(dpos=0.5, dmag=0.5, src=None, alpha=0.0, moff=False),
           'tol35': dict(dpos=0.35, dmag=0.5, src=None, alpha=0.0, moff=False), 'mag30': dict(dpos=0.5, dmag=0.3, src=None, alpha=0.0, moff=False),
           'tol75': dict(dpos=0.75, dmag=0.75, src=None, alpha=0.0, moff=False)}
    for src in ('scan', 'own'):
        for a in ALPHAS[1:]:
            out['%s%d' % (src, round(100 * a))] = dict(dpos=0.5, dmag=0.5, src=src, alpha=a, moff=False)
    out['ownp100'] = dict(dpos=0.5, dmag=0.5, src='ownp', alpha=1.0, moff=False)
    out['mo'] = dict(dpos=0.5, dmag=0.5, src=None, alpha=0.0, moff=True)
    out['mo+scan100'] = dict(dpos=0.5, dmag=0.5, src='scan', alpha=1.0, moff=True)
    out['mo+own100'] = dict(dpos=0.5, dmag=0.5, src='own', alpha=1.0, moff=True)
    return out


SCENARIOS = _scenarios()


def scn_name(src, alpha):
    """The scenario of a shift source and alpha: 'base' for alpha 0 (any source), else e.g. 'scan50', 'own100'."""
    return 'base' if alpha == 0 else '%s%d' % (src, round(100 * alpha))


# ----------------------------------------------------------------------------------------------------------------------------------
# the scorer: the matching with tolerances as arguments, the counts at several depth cuts, the blocks of several sizes
# ----------------------------------------------------------------------------------------------------------------------------------
def match_pairs(det_xy, det_mag, tr_xy, tr_mag, dpos=A.DPOS, dmag=A.DMAG):
    """One-to-one matching of detections to truth sources, of maximum size and, among those, of minimum total distance, over the pairs within dpos pixels AND dmag magnitudes (both inclusive): the
    objective of `astro_common.match`, solved by maximum weight (weight BIG - distance of an admissible pair) per connected component of the admissible pairs. Returns (detection index, truth index)
    arrays, sorted by detection."""
    det_xy, tr_xy = np.asarray(det_xy, float).reshape(-1, 2), np.asarray(tr_xy, float).reshape(-1, 2)
    det_mag, tr_mag = np.asarray(det_mag, float), np.asarray(tr_mag, float)
    nd, nt = len(det_xy), len(tr_xy)
    empty = (np.zeros(0, int), np.zeros(0, int))
    if nd == 0 or nt == 0:
        return empty
    cand = cKDTree(tr_xy).query_ball_point(det_xy, dpos + 1e-6)             # candidates only; the tests below are the exact ones
    di = np.repeat(np.arange(nd), [len(c) for c in cand])
    ti = np.fromiter((j for c in cand for j in c), dtype=int, count=len(di))
    dx, dy = det_xy[di, 0] - tr_xy[ti, 0], det_xy[di, 1] - tr_xy[ti, 1]
    d = np.sqrt(dx * dx + dy * dy)
    with np.errstate(invalid='ignore'):
        ok = (d <= dpos) & (np.abs(det_mag[di] - tr_mag[ti]) <= dmag)
    di, ti, d = di[ok], ti[ok], d[ok]
    if len(di) == 0:
        return empty
    _, lab = connected_components(coo_matrix((np.ones(len(di)), (di, nd + ti)), shape=(nd + nt, nd + nt)), directed=False)
    order = np.argsort(lab[di], kind='stable')
    di, ti, d = di[order], ti[order], d[order]
    cl = lab[di]
    bounds = np.flatnonzero(np.r_[True, cl[1:] != cl[:-1], True])
    md, mt = [], []
    for a, b in zip(bounds[:-1], bounds[1:]):
        cd, ct, cdist = di[a:b], ti[a:b], d[a:b]
        ud, ui = np.unique(cd, return_inverse=True)
        ut, uj = np.unique(ct, return_inverse=True)
        W = np.zeros((len(ud), len(ut)))
        W[ui, uj] = BIG - cdist
        r, c = linear_sum_assignment(W, maximize=True)
        keep = W[r, c] > 0
        md += list(ud[r[keep]])
        mt += list(ut[c[keep]])
    md, mt = np.array(md, int), np.array(mt, int)
    o = np.argsort(md, kind='stable')
    return md[o], mt[o]


def block_index(xy, cell, n=N):
    """Row-major index of the block of `cell` px of each position (x, y) in patch pixels: clip(floor((position + 0.5) / cell), 0, n / cell - 1) (as `astro_common.block_of` for cell = 20)."""
    xy = np.atleast_2d(np.asarray(xy, float))
    nb = n // cell
    ix = np.clip(np.floor((xy[:, 0] + 0.5) / cell).astype(int), 0, nb - 1)
    iy = np.clip(np.floor((xy[:, 1] + 0.5) / cell).astype(int), 0, nb - 1)
    return iy * nb + ix


def rescore(det, tp, dpos=A.DPOS, dmag=A.DMAG, shift=(0.0, 0.0), moff=0.0, cuts=CUTS, cells=(A.CELL,), n=N):
    """Our scorer (D0): the catalogue det = (x, y, flux) against the truth patch tp (`astro_common.truth_for_patch`), tolerances dpos, dmag, the truth moved by `shift` for the matching only, the
    magnitudes of the detections minus `moff` (a number, or 'auto': the median offset of the matched pairs with MAG_LO < truth magnitude < MAG_HI, at least MAG_MIN_PAIRS of them, else 0, estimated
    in a first pass without correction). Returns dict(tot (len(cuts), 7) int, blk {cell: (len(cuts), (n / cell)^2, 4)}, off (2,) the median (detection - moved truth) of the matches brighter than 21,
    moff the offset used, md and mt the matched detection and truth indices)."""
    x, y, f = (np.asarray(v, float) for v in det)
    mag = A.counts_to_mag(f)
    inside = (x >= A.EDGE - 0.5) & (x <= n - 0.5 - A.EDGE) & (y >= A.EDGE - 0.5) & (y <= n - 0.5 - A.EDGE)
    di = np.flatnonzero(np.isfinite(mag) & inside)
    ti = np.flatnonzero(tp['inner'])
    sx, sy = float(shift[0]), float(shift[1])
    tx, ty, tm = tp['x'], tp['y'], tp['mag']
    if isinstance(moff, str):
        assert moff == 'auto'
        r0, t0 = match_pairs(np.c_[x[di], y[di]], mag[di], np.c_[tx[ti] + sx, ty[ti] + sy], tm[ti], dpos, dmag)
        tmm = tm[ti][t0]
        sel = (tmm > MAG_LO) & (tmm < MAG_HI)
        moff = float(np.median(mag[di][r0][sel] - tmm[sel])) if sel.sum() >= MAG_MIN_PAIRS else 0.0
    mc = mag - moff
    rd, rt = match_pairs(np.c_[x[di], y[di]], mc[di], np.c_[tx[ti] + sx, ty[ti] + sy], tm[ti], dpos, dmag)
    md, mt = di[rd], ti[rt]
    tot = np.zeros((len(cuts), 7), int)
    blk = {c_: np.zeros((len(cuts), (n // c_) ** 2, 4)) for c_ in cells}
    for k, c in enumerate(cuts):
        s_t, s_d, s_mt, s_md = tm[ti] < c, mc[di] < c, tm[mt] < c, mc[md] < c
        tot[k] = [s_t.sum(), s_mt.sum(), s_d.sum(), s_md.sum(), len(ti), len(mt), len(di)]
        for cell in cells:
            np.add.at(blk[cell][k][:, 0], block_index(np.c_[tx[ti[s_t]], ty[ti[s_t]]], cell, n), 1)
            np.add.at(blk[cell][k][:, 1], block_index(np.c_[tx[mt[s_mt]], ty[mt[s_mt]]], cell, n), 1)
            np.add.at(blk[cell][k][:, 2], block_index(np.c_[x[di[s_d]], y[di[s_d]]], cell, n), 1)
            np.add.at(blk[cell][k][:, 3], block_index(np.c_[x[md[s_md]], y[md[s_md]]], cell, n), 1)
    bright = tm[mt] < 21.0
    off = np.median(np.c_[x[md] - (tx[mt] + sx), y[md] - (ty[mt] + sy)][bright], axis=0) if bright.any() else np.full(2, np.nan)
    return dict(tot=tot, blk=blk, off=off, moff=float(moff), md=md, mt=mt)


# ----------------------------------------------------------------------------------------------------------------------------------
# intervals: the report's block bootstrap, in chunks, with other numbers of resamples and seeds
# ----------------------------------------------------------------------------------------------------------------------------------
def _chunks(n, size):
    return [(i, min(i + size, n)) for i in range(0, n, size)]


def resample_sums(blk, idx):
    """The sums over the resampled blocks, (B, 4) for the block counts blk (nb, 4) and the resamples idx (B, m) of block indices: identical to blk[idx].sum(1) (the counts are integers, so the sums are
    exact in any order), computed as the matrix product of the (B, nb) matrix of how often each block was drawn with blk."""
    B, nb = idx.shape[0], blk.shape[0]
    cnt = np.bincount((np.arange(B)[:, None] * nb + idx).ravel(), minlength=B * nb).reshape(B, nb)
    return cnt @ blk


def boot_stat(fn, blks, sizes, B=BOOT, seed=0, chunk=5000):
    """(value, (lo, hi)): fn(*[blk.sum(0) for blk in blks]) on the blocks and its 2.5 and 97.5 percentiles over B stratified resamples (`exp_astro_m2.boot_idx_strat`, seed): the same resamples for
    every array of blks, so that differences and ratios are paired. fn takes the summed counts (..., 4) of each array and returns an array over the resamples."""
    idx = X.boot_idx_strat(list(sizes), B, seed)
    v0 = fn(*[b.sum(0) for b in blks])
    vs = np.concatenate([fn(*[resample_sums(b, idx[i:j]) for b in blks]) for i, j in _chunks(B, chunk)])
    return float(v0), np.percentile(vs, [2.5, 97.5])


def _dF(a, b):
    return X.cpf(a)[2] - X.cpf(b)[2]


def _ratio(a, b):
    return X.cpf(a)[2] / X.cpf(b)[2]


def boot_dF(blk_a, blk_b, sizes, B=BOOT, seed=0):
    """The paired dF = F(a) - F(b) with its 95 % interval; at B = BOOT, seed 0 it is `exp_astro_m2.boot_diff(a, b, boot_idx_strat(sizes))[2]`."""
    return boot_stat(_dF, [blk_a, blk_b], sizes, B, seed)


def boot_ratio(blk_a, blk_b, sizes, B=BOOT, seed=0):
    """The paired ratio F(a) / F(b) with its 95 % interval."""
    return boot_stat(_ratio, [blk_a, blk_b], sizes, B, seed)


def seed_lowers(blk_a, blk_b, sizes, seeds=None, B=BOOT):
    """The lower end of the 95 % interval of dF for each seed (B stratified resamples each); the seeds are SEEDS (1 to 200) unless given."""
    return np.array([boot_dF(blk_a, blk_b, sizes, B, s)[1][0] for s in (SEEDS if seeds is None else seeds)])


# ----------------------------------------------------------------------------------------------------------------------------------
# the offset between the truth and the image: template fits of bright isolated stars with the patch's own stamp (D1)
# ----------------------------------------------------------------------------------------------------------------------------------
def box_rows(truth, x0, y0, n=N):
    """The truth rows inside the box [x0 - 0.5, x0 + n - 0.5)^2 of a patch, in patch pixels: dict(x, y, mag, nf) (`truth`: arrays px, py, m606, nf606 of a whole truth list)."""
    sel = (truth['px'] >= x0 - 0.5) & (truth['px'] < x0 + n - 0.5) & (truth['py'] >= y0 - 0.5) & (truth['py'] < y0 + n - 0.5)
    return dict(x=truth['px'][sel] - x0, y=truth['py'][sel] - y0, mag=truth['m606'][sel], nf=truth['nf606'][sel])


def select_stars(rows, n=N, mlo=15.5, mhi=19.5, margin=15.0, riso=8.0, dmag_iso=2.0):
    """The indices (into `rows`, see `box_rows`) of the stars of the offset fit: mlo <= F606W < mhi, not HST-saturated (nf606 != 99), at least `margin` px from every edge (margin <= x, y <= n - 1 - margin)
    and no OTHER row of the box (any magnitude, saturated ones included) closer than `riso` px that is brighter than m + dmag_iso."""
    x, y, m = rows['x'], rows['y'], rows['mag']
    ok = (m >= mlo) & (m < mhi) & (rows['nf'] != 99) & (x >= margin) & (x <= n - 1 - margin) & (y >= margin) & (y <= n - 1 - margin)
    keep = []
    for k in np.flatnonzero(ok):
        d = np.hypot(x - x[k], y - y[k])
        d[k] = np.inf
        if not ((d < riso) & (m < m[k] + dmag_iso)).any():
            keep.append(k)
    return np.array(keep, int)


def pixel_weights(cut, gain, dark):
    """The weights of the star fit: 1 / sqrt(max(DN, 1) / gain + dark variance), from the counts of the cut itself."""
    return 1.0 / np.sqrt(np.maximum(cut, 1.0) / gain + dark)


def fit_star(dn, coef, gain, dark, tx, ty, half=STAR['half']):
    """Fit A S(x - cx, y - cy) + B to the (2 half + 1)^2 cut of the image `dn` (DN, rows are y) around the rounded position (tx, ty): S is the stamp whose spline coefficients are `coef`
    (`scipy.ndimage.spline_filter` of the unit-sum stamp; origin: its middle pixel), weights 1 / sqrt(max(DN, 1) / gain + dark). Returns (cx, cy, A, B, converged)."""
    ix, iy = int(np.floor(tx + 0.5)), int(np.floor(ty + 0.5))
    cut = dn[iy - half: iy + half + 1, ix - half: ix + half + 1]
    ys, xs = np.mgrid[iy - half: iy + half + 1, ix - half: ix + half + 1].astype(float)
    w = pixel_weights(cut, gain, dark)
    c = (coef.shape[0] - 1) / 2.0

    def resid(q):
        s = map_coordinates(coef, [ys - q[1] + c, xs - q[0] + c], order=3, mode='constant', cval=0.0, prefilter=False)
        return ((q[2] * s + q[3] - cut) * w).ravel()
    a0 = float(cut.sum() - cut.size * np.median(cut))
    sol = least_squares(resid, [tx, ty, a0, float(np.median(cut))], x_scale=[0.2, 0.2, max(abs(a0), 1.0) * 0.1, 5.0])
    return float(sol.x[0]), float(sol.x[1]), float(sol.x[2]), float(sol.x[3]), bool(sol.success)


def measure_patch(dn, psf, gain, dark, rows, selection='B', n=N):
    """The offsets (fitted - truth position) of the stars of a patch selected by SELECTIONS[selection]: (array (k, 6) = x, y, mag, dx, dy, flux of each star that was kept, number selected). A fit that
    did not converge or with an offset above STAR['max_offset'] px is dropped."""
    st = np.asarray(psf, float)
    coef = spline_filter(st / st.sum(), order=3)
    sel = select_stars(rows, n, **SELECTIONS[selection])
    out = []
    for k in sel:
        cx, cy, a, b, ok = fit_star(np.asarray(dn, float), coef, gain, dark, rows['x'][k], rows['y'][k])
        dx, dy = cx - rows['x'][k], cy - rows['y'][k]
        if ok and np.hypot(dx, dy) <= STAR['max_offset']:
            out.append([rows['x'][k], rows['y'][k], rows['mag'][k], dx, dy, a])
    return np.array(out, float).reshape(-1, 6), len(sel)


def rob_sd(v):
    """The robust standard deviation 1.4826 MAD about the median (NaN for an empty sample)."""
    v = np.asarray(v, float)
    return float(1.4826 * np.median(np.abs(v - np.median(v)))) if len(v) else float('nan')


def median_interval(v, B=BOOT, seed=0):
    """(median, (lo, hi)): the median of v and the 2.5 and 97.5 percentiles of the medians of B resamples of v (seed)."""
    v = np.asarray(v, float)
    idx = np.random.default_rng(seed).integers(0, len(v), size=(B, len(v)))
    return float(np.median(v)), np.percentile(np.median(v[idx], axis=1), [2.5, 97.5])


def frame_offsets(stars_by_patch, min_stars=STAR['min_stars']):
    """The `own` offsets of a frame from {patch: stars (k, 6)}: (frame median (dx, dy) or None if fewer than min_stars stars in all, {patch: median (dx, dy) of the patch, the frame's where the patch
    has fewer than min_stars})."""
    allst = np.concatenate([s for s in stars_by_patch.values()]) if stars_by_patch else np.zeros((0, 6))
    if len(allst) < min_stars:
        return None, {}
    fm = np.median(allst[:, 3:5], axis=0)
    return fm, {p: (np.median(s[:, 3:5], axis=0) if len(s) >= min_stars else fm) for p, s in stars_by_patch.items()}


def unit_pair_difference(a, b, mmax, rad=0.4, lo=0.8, hi=1.25):
    """The flux-consistent unit pairs of two catalogues (x, y, flux) a and b with detected magnitude < mmax: mutual nearest neighbours closer than rad px whose flux ratio a / b is in [lo, hi];
    returns (number of pairs, median of (a - b) x, median y), NaN medians for no pair."""
    a, b = np.asarray(a, float).reshape(-1, 3), np.asarray(b, float).reshape(-1, 3)
    a, b = a[A.counts_to_mag(a[:, 2]) < mmax], b[A.counts_to_mag(b[:, 2]) < mmax]
    if len(a) == 0 or len(b) == 0:
        return 0, float('nan'), float('nan')
    da, ja = cKDTree(b[:, :2]).query(a[:, :2])
    _, ib = cKDTree(a[:, :2]).query(b[:, :2])
    mutual = (ib[ja] == np.arange(len(a))) & (da < rad)
    ratio = a[:, 2] / b[ja, 2]
    good = mutual & (ratio >= lo) & (ratio <= hi)
    if not good.any():
        return 0, float('nan'), float('nan')
    d = a[good, :2] - b[ja[good], :2]
    return int(good.sum()), float(np.median(d[:, 0])), float(np.median(d[:, 1]))


# ----------------------------------------------------------------------------------------------------------------------------------
# the raw: everything the report prints, from the stored results of the finished test and the re-scored variants
# ----------------------------------------------------------------------------------------------------------------------------------
def pk(name, pxy):
    """The key fragment of a patch: 'm2r2|x133_y281'."""
    return '%s|x%d_y%d' % (name, pxy[0], pxy[1])


def catalogues(rep, base, name, pxy):
    """{row key: (x, y, flux)} of the rows that can be re-scored on a repeat patch: the stored detections of the ISS (SC_KEYS and BASE_EXTRA) and the baselines' catalogues through the frozen
    post-filter (`exp_astro_m2.cs_dets`) and tags."""
    pre = E.pre_of(name, pxy)
    cats = {}
    for key in ISS_KEYS + BASE_EXTRA[:2]:
        kern, cfg, _ = key.split('|')
        d = rep[pre + '%s|%s|det' % (kern, cfg)]
        cats[key] = (d[:, 0], d[:, 1], d[:, 2])
    fmap = {X.cs_name(f): f for f in X.CS_FILTERS}
    cats['crowdsource (default)'] = X.cs_dets(base, pxy, (0.0, 0.0, 0.0, np.inf))
    cats['crowdsource (tuned)'] = X.cs_dets(base, pxy, fmap[str(rep[pre + 'crowdsource (tuned)|name'][0])])
    for lab, tag in (('photutils (default)', 'photutils'), ('photutils (tuned)', str(rep[pre + 'photutils (tuned)|name'][0]))):
        k = 'x%d_y%d|%s|' % (pxy[0], pxy[1], tag)
        cats[lab] = (base[k + 'x'], base[k + 'y'], base[k + 'flux'])
    return cats


def shift_of(spec, name, pxy, scan, own, ownp):
    """The shift (px) of the truth of one patch in a scenario: alpha x the offset of the scenario's source (scan: of the patch; own: the frame's median; ownp: the patch's own median)."""
    a, s = spec['alpha'], spec['src']
    if s is None or a == 0:
        return (0.0, 0.0)
    o = {'scan': lambda: scan[(name, pxy)], 'own': lambda: own[name], 'ownp': lambda: ownp[(name, pxy)]}[s]()
    return (a * float(o[0]), a * float(o[1]))


def put(out, prefix, r):
    """Stores a `rescore` result under prefix: tot, blk (cell A.CELL), blkNN for the other cells, off, moff."""
    out[prefix + '|tot'] = r['tot']
    out[prefix + '|blk'] = r['blk'][A.CELL]
    for c_, b in r['blk'].items():
        if c_ != A.CELL:
            out[prefix + '|blk%d' % c_] = b
    out[prefix + '|off'] = r['off']
    out[prefix + '|moff'] = np.array([r['moff']])


def build_scenarios(rep, base, truths, scan, own, ownp, out):
    """Re-scores every scenario (SCENARIOS) on every repeat patch and stores it under 'sc|<scenario>|<frame>|x.._y..|<row>|...'. The base scenario must reproduce the stored rows exactly."""
    for name, pxy in PATCHES:
        tp = A.truth_for_patch(truths[name], *pxy)
        cats = catalogues(rep, base, name, pxy)
        pre = E.pre_of(name, pxy)
        for scn, spec in SCENARIOS.items():
            shift = shift_of(spec, name, pxy, scan, own, ownp)
            out['shift|%s|%s' % (scn, pk(name, pxy))] = np.array(shift, float)
            for key in SC_KEYS + (BASE_EXTRA if scn == 'base' else ()):
                cells = CELLS if (scn == 'base' and key in ISS_KEYS) else (A.CELL,)
                r = rescore(cats[key], tp, spec['dpos'], spec['dmag'], shift, 'auto' if spec['moff'] else 0.0, CUTS, cells)
                if scn == 'base':
                    if not (np.array_equal(r['tot'][ICUT], rep[pre + key + '|tot']) and np.array_equal(r['blk'][A.CELL][ICUT], rep[pre + key + '|blk22'])
                            and np.array_equal(r['off'], rep[pre + key + '|offset'], equal_nan=True)):
                        raise RuntimeError('the re-scored base row %s %s differs from the stored row' % (pk(name, pxy), key))
                put(out, 'sc|%s|%s|%s' % (scn, pk(name, pxy), key), r)


def build_stored(rep, m2raw, out):
    """Copies of the stored rows the report uses (tot, blk22, offset, the bright-star counts below 19.5 mag, the run records, the stamp fits), so that the report needs this raw alone."""
    for name, pxy in PATCHES:
        pre, p = E.pre_of(name, pxy), 'st|' + pk(name, pxy)
        for row in STORED_ROWS:
            out['%s|%s|tot' % (p, row)] = rep[pre + row + '|tot']
            out['%s|%s|blk22' % (p, row)] = rep[pre + row + '|blk22']
            out['%s|%s|offset' % (p, row)] = rep[pre + row + '|offset']
            out['%s|%s|bright' % (p, row)] = np.array([rep[pre + row + '|n_matched'][:9].sum(), rep[pre + row + '|n_true'][:9].sum()])
        for kern in ('iso', 'adq', 'scl'):
            for cfg in ('ours', 'default'):
                for k in ('bic', 'rsd', 'secs'):
                    out['%s|%s|%s|%s' % (p, kern, cfg, k)] = rep[pre + '%s|%s|%s' % (kern, cfg, k)]
        for k in ('kind', 'rms', 'shape', 'fwhm', 'scale', 'kernel|same', 'cfg'):
            out['%s|%s' % (p, k)] = rep[pre + k]
        out[p + '|parts'] = np.array([len(rep[pre + 'parts5']) if str(rep[pre + 'kind'][0]) == 'general' else 3])
        out[p + '|photutils (tuned)|name'] = rep[pre + 'photutils (tuned)|name']
        out[p + '|crowdsource (tuned)|name'] = rep[pre + 'crowdsource (tuned)|name']
    for w in X.TESTS:
        for row in M2_ROWS:
            out['m2|x%d_y%d|%s|tot' % (w[0], w[1], row)] = m2raw['T|x%d_y%d|%s|tot' % (w[0], w[1], row)]
            out['m2|x%d_y%d|%s|blk22' % (w[0], w[1], row)] = m2raw['T|x%d_y%d|%s|blk22' % (w[0], w[1], row)]
        out['m2|x%d_y%d|cfg' % w] = m2raw['T|x%d_y%d|cfg' % w]


def datasets():
    """The data sets whose offsets are measured: [(label, patches, loader (x0, y0) -> patch dict, path of the truth list)]: the three repeat frames, the M2 test windows, M53 and M92."""
    out = [(n, list(R.REPEAT[n]['patches']), (lambda x0, y0, n=n: R.load_patch(x0, y0, n)), os.path.join(R.REPEAT[n]['root'], R.REPEAT[n]['truth'])) for n in FRAMES]
    out.append(('M2', list(X.TESTS), lambda x0, y0: X.load_patch(x0, y0, 'm2'), os.path.join(X.ROOT, 'truth_hugs_meth1.npz')))
    for ds, lab in (('m53', 'M53'), ('m92', 'M92')):
        out.append((lab, list(X.SECOND[ds]['patches']), (lambda x0, y0, ds=ds: X.load_patch(x0, y0, ds)), os.path.join(X.SECOND[ds]['root'], X.SECOND[ds]['truth'])))
    return out


def measure_all(dsets):
    """{label: {patch: dict(stars={selection: (k, 6)}, nsel={selection: number selected}, gain, dark, fwhm)}} for the data sets of `datasets()` and both selections."""
    res = {}
    for label, plist, loader, tpath in dsets:
        z = np.load(tpath)
        truth = {k: z[k] for k in ('px', 'py', 'm606', 'nf606')}
        res[label] = {}
        for pxy in plist:
            patch = loader(*pxy)
            rows = box_rows(truth, pxy[0], pxy[1])
            m = {s_: measure_patch(patch['dn'], patch['psf'], patch['gain'], patch['dark'], rows, s_) for s_ in SELECTIONS}
            res[label][tuple(pxy)] = dict(stars={s_: v[0] for s_, v in m.items()}, nsel={s_: v[1] for s_, v in m.items()}, gain=patch['gain'], dark=patch['dark'], fwhm=S2.psf_fwhm(patch['psf']))
    return res


def neighbour_bias(psf, gain, dark, ds=BIAS_D, qs=BIAS_Q, angles=8):
    """(len(ds), len(qs)) array: the position error (px, noise-free) of the single-star fit `fit_star` on an image made of the stamp itself at a known position and ONE neighbour of flux ratio q at
    distance d, the largest over `angles` directions of the neighbour; the justification of the isolation radius of selection B."""
    st = np.asarray(psf, float)
    coef = spline_filter(st / st.sum(), order=3)
    yy, xx = np.mgrid[0:60, 0:60].astype(float)
    c = (st.shape[0] - 1) / 2.0

    def render(cx, cy):
        return map_coordinates(coef, [yy - cy + c, xx - cx + c], order=3, mode='constant', cval=0.0, prefilter=False)
    cx, cy, amp = 30.3, 29.8, 5000.0
    out = np.zeros((len(ds), len(qs)))
    for i, d in enumerate(ds):
        for j, q in enumerate(qs):
            for th in np.arange(angles) * 2 * np.pi / angles:
                img = 150.0 + amp * render(cx, cy) + q * amp * render(cx + d * np.cos(th), cy + d * np.sin(th))
                fx, fy, _, _, _ = fit_star(img, coef, gain, dark, cx, cy)
                out[i, j] = max(out[i, j], np.hypot(fx - cx, fy - cy))
    return out


def build_offsets(meas, scan, out):
    """Stores the star fits per patch ('own|<label>|x.._y..' (k, 6), '|nsel'), and the patch information of the confound table ('pt|...|info' = gain, dark variance, FWHM)."""
    for label, d in meas.items():
        for pxy, v in d.items():
            p = '%s|x%d_y%d' % (label, pxy[0], pxy[1])
            for s_ in SELECTIONS:
                out['own|%s|%s' % (s_, p)] = v['stars'][s_]
                out['own|%s|%s|nsel' % (s_, p)] = np.array([v['nsel'][s_]])
            out['pt|%s|info' % p] = np.array([v['gain'], v['dark'], v['fwhm']])
    for (label, pxy), v in scan.items():
        out['scan|%s|x%d_y%d' % (label, pxy[0], pxy[1])] = np.array(v)


def build_positions(rep, out):
    """The position difference of the two kernels' unit pairs on the patches where the rule picks a general kernel (D10), at the detected-magnitude limits 19, 20, 21, 22."""
    for name, pxy in PATCHES:
        pre = E.pre_of(name, pxy)
        if rep[pre + 'kernel|same'][0] == 1:
            continue
        for m in (19.0, 20.0, 21.0, 22.0):
            n_, dx, dy = unit_pair_difference(rep[pre + 'adq|ours|det'], rep[pre + 'iso|ours|det'], m)
            out['pos|%s|m%d' % (pk(name, pxy), m)] = np.array([n_, dx, dy])


def guarded_rows(path):
    """The guarded re-runs of the repeat paths in a cyclefix raw (key format: see the docstring, D8): {(frame, patch, kern, cfg): {'k2'|'k5'|'k0'...: (tot (7,), blk22 (25, 4)), 'bic': (3,), 'secs': (2,)}}; {} if the
    file is absent or holds no repeat run."""
    out = {}
    if not path or not os.path.exists(path):
        return out
    with np.load(path, allow_pickle=False) as z:
        for k in z.files:
            q = k.split('|')
            if len(q) < 7 or q[0] != 'G' or q[1] != 'R' or q[2] not in R.REPEAT or q[4] not in ('iso', 'adq', 'scl') or q[5] not in ('ours', 'default'):
                continue
            x0, y0 = q[3][1:].split('_y')
            run = out.setdefault((q[2], (int(x0), int(y0)), q[4], q[5]), {})
            if len(q) == 8 and q[7] in ('tot', 'blk22') and q[6].startswith('k'):
                run.setdefault(q[6], {})['tot' if q[7] == 'tot' else 'blk'] = z[k]
            elif len(q) == 7 and q[6] in ('bic', 'secs'):
                run[q[6]] = z[k]
    return {k: v for k, v in out.items() if any(isinstance(x, dict) and 'tot' in x and 'blk' in x for x in v.values())}


def build_guarded(g, out):
    """Copies the guarded rows into the raw: 'gr|<frame>|x.._y..|<kern>|<cfg>|<k>|tot', '|blk22', '|bic', '|secs'."""
    for (name, pxy, kern, cfg), run in sorted(g.items()):
        p = 'gr|%s|%s|%s|' % (pk(name, pxy), kern, cfg)
        for k, v in sorted(run.items()):
            if isinstance(v, dict):
                if 'tot' in v and 'blk' in v:
                    out[p + k + '|tot'], out[p + k + '|blk22'] = np.asarray(v['tot']), np.asarray(v['blk'])
            else:
                out[p + k] = np.asarray(v)


def read_scan(path=VALID):
    """(the masked asTrans scan offsets {(frame, patch): (dx, dy)} of the repeat patches, the same for the M2 windows {('M2', window): ...} (masked, header mapping), the plain asTrans scan of the repeat
    patches), all read from the machine-readable copy of repeat_spec.md section 6.2."""
    with open(path) as f:
        d = json.load(f)
    scan, scan_m2, plain = {}, {}, {}
    for name, pxy in PATCHES:
        e = d[name]['%d,%d' % pxy]
        scan[(name, pxy)] = (float(e['masked']['asTrans']['dx']), float(e['masked']['asTrans']['dy']))
        plain[(name, pxy)] = (float(e['image_shift_px']['asTrans']['dx']), float(e['image_shift_px']['asTrans']['dy']))
    for w in X.TESTS:
        e = d['m2']['%d,%d' % w]
        scan_m2[('M2', tuple(w))] = (float(e['masked']['header']['dx']), float(e['masked']['header']['dy']))
    return scan, scan_m2, plain


def own_offsets(meas):
    """({frame: (dx, dy) frame median or absent}, {(frame, patch): (dx, dy)}) from the star fits of the repeat frames."""
    own, ownp = {}, {}
    for name in FRAMES:
        fm, per = frame_offsets({p: v['stars'][OWN_SEL] for p, v in meas[name].items()})
        if fm is not None:
            own[name] = (float(fm[0]), float(fm[1]))
            ownp.update({(name, p): (float(o[0]), float(o[1])) for p, o in per.items()})
    return own, ownp


def build_stamps(psfraw, out):
    """The stamp fits of the patches of M2 (test windows), M53 and M92 from `psf_raw.npz` (the repeat patches' are in the repeat raw): 'pt|<label>|x.._y..|rms' (isotropic rms, chosen rms),
    '|kind', '|shape', '|parts'."""
    for label, ds, plist in (('M2', 'm2', X.TESTS), ('M53', 'm53', X.SECOND['m53']['patches']), ('M92', 'm92', X.SECOND['m92']['patches'])):
        for pxy in plist:
            s = 'S|%s|x%d_y%d|' % (ds, pxy[0], pxy[1])
            p = 'pt|%s|x%d_y%d|' % (label, pxy[0], pxy[1])
            out[p + 'rms'], out[p + 'kind'], out[p + 'shape'] = psfraw[s + 'rms'], psfraw[s + 'kind'], psfraw[s + 'shape']
            out[p + 'parts'] = np.array([len(psfraw[s + 'parts5']) if str(psfraw[s + 'kind'][0]) == 'general' else 3])


def build_raw(inp):
    """The whole raw from the loaded inputs `inp` (see `load_inputs`): dict key -> array."""
    out = {}
    own, ownp = own_offsets(inp['meas'])
    lacking = [n for n in FRAMES if n not in own]
    if lacking:
        raise RuntimeError('fewer than %d fitted stars in the frames %s: no own offset' % (STAR['min_stars'], lacking))
    build_scenarios(inp['rep'], inp['base'], inp['truth'], inp['scan'], own, ownp, out)
    build_stored(inp['rep'], inp['m2raw'], out)
    build_stamps(inp['psfraw'], out)
    build_offsets(inp['meas'], {**inp['scan'], **inp['scan_m2']}, out)
    for (name, pxy), v in inp['plain'].items():
        out['plain|%s' % pk(name, pxy)] = np.array(v)
    build_positions(inp['rep'], out)
    build_guarded(inp.get('guarded', {}), out)
    out['meta|cyclefix'] = np.array(['applied' if inp.get('guarded') else 'not applied'])
    out['bias|table'] = neighbour_bias(inp['bias_psf'], inp['bias_gain'], inp['bias_dark'])
    for k, v in sorted(inp['sha256'].items()):
        out['meta|sha256|' + k] = np.array([v])
    out['meta|scan_key'] = np.array(['%s[<frame>]["<x0>,<y0>"]["masked"]["asTrans"]["dx"|"dy"] (M2 windows: [m2]["<x0>,<y0>"]["masked"]["header"]); the plain scan: ["image_shift_px"]["asTrans"]' % os.path.basename(VALID)])
    return out


def _sha256_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def load_inputs(cyclefix=None):
    """Reads the finished test's results and the data (all read only): the stored raws and catalogues, the truth lists, the scan offsets, the star fits, and the sha256 of the input files; with `cyclefix`
    (the path of a cyclefix raw) also its guarded re-runs of the repeat paths."""
    scan, scan_m2, plain = read_scan()
    bp = R.load_patch(133, 281, 'm2r2')                                                      # the stamp of the neighbour-bias simulation
    files = {'repeat_raw.npz': REPEAT_RAW, 'baseline_catalogs_repeat.npz': BASE, 'm2_raw.npz': M2_RAW, 'psf_raw.npz': PSF_RAW, 'repeat_validation.json': VALID}
    files.update({os.path.basename(R.REPEAT[n]['truth']): os.path.join(R.REPEAT[n]['root'], R.REPEAT[n]['truth']) for n in FRAMES})
    guarded = guarded_rows(cyclefix)
    if guarded:
        files['cyclefix_raw.npz'] = cyclefix
    return dict(guarded=guarded, rep=X.load_raw(REPEAT_RAW), base=np.load(BASE, allow_pickle=False), m2raw=X.load_raw(M2_RAW), psfraw=X.load_raw(PSF_RAW), truth={n: R.truth_arrays(n) for n in FRAMES},
                scan=scan, scan_m2=scan_m2, plain=plain, meas=measure_all(datasets()), sha256={k: _sha256_file(p) for k, p in files.items()},
                bias_psf=bp['psf'], bias_gain=bp['gain'], bias_dark=bp['dark'])


def run(out_path=RAW, cyclefix=None):
    """`--run`: builds the raw from the data and writes it (it is a deterministic function of its inputs, so an existing file is replaced); `cyclefix`: the path of a cyclefix raw whose guarded re-runs of the
    repeat paths are copied in (the optional hook)."""
    raw = build_raw(load_inputs(cyclefix))
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    np.savez_compressed(out_path, **raw)
    print('wrote %s (%d arrays)' % (out_path, len(raw)))
    return raw


def load_raw(path=RAW):
    return X.load_raw(path)


# ----------------------------------------------------------------------------------------------------------------------------------
# report
# ----------------------------------------------------------------------------------------------------------------------------------
GEN_PATCHES = [('m2r2', (133, 281)), ('m2r2', (195, 689))]               # the patches where the rule picks a general kernel (checked against the raw by `general_patches`)
P2_PAIRS = [('chosen kernel vs crowdsource (tuned)', 'adq|ours|k2', 'crowdsource (tuned)'), ('chosen kernel vs photutils (tuned)', 'adq|ours|k2', 'photutils (tuned)'),
            ('isotropic kernel vs crowdsource (tuned)', 'iso|ours|k2', 'crowdsource (tuned)'), ('isotropic kernel vs photutils (tuned)', 'iso|ours|k2', 'photutils (tuned)')]
READOUTS = [('tuned configuration (kappa 2)', 'adq|ours|k2', 'iso|ours|k2'), ('a priori default (kappa 5)', 'adq|default|k5', 'iso|default|k5'),
            ('a priori default with the BIC stop', 'adq|default|k0', 'iso|default|k0')]
KEY_LABEL = {'iso|ours|k2': 'ISS isotropic, tuned', 'adq|ours|k2': 'ISS chosen kernel, tuned', 'scl|ours|k2': 'ISS chosen kernel, scaled lengths, tuned', 'iso|default|k5': 'ISS isotropic, default',
             'adq|default|k5': 'ISS chosen kernel, default', 'scl|default|k5': 'ISS chosen kernel, scaled lengths, default', 'iso|default|k0': 'ISS isotropic, BIC stop',
             'adq|default|k0': 'ISS chosen kernel, BIC stop', 'scl|default|k0': 'ISS chosen kernel, scaled lengths, BIC stop', 'crowdsource (tuned)': 'crowdsource (tuned)',
             'photutils (tuned)': 'photutils (tuned)', 'crowdsource (default)': 'crowdsource (default)', 'photutils (default)': 'photutils (default)'}
DEFECT_NOTE = 'D: a compared ISS run hit the entry cap (a defective path, section 8)'


def ptag(name, pxy):
    return '%s (%d, %d)' % (name, pxy[0], pxy[1])


def blocks(raw, plist, key, scn='base', cut=22.0, cell=A.CELL):
    """(concatenated blocks (nb, 4), [blocks per patch]) of a row over the patches: the scenario's re-scored row, or the stored row for the BIC readout (key ending in k0: base scenario, cut 22, 20 px)."""
    if key.endswith('|k0'):
        assert scn == 'base' and cut == 22.0 and cell == A.CELL, 'the BIC readout has no stored detections: it cannot be re-scored'
        parts = [raw['st|%s|%s|blk22' % (pk(n, p), key)] for n, p in plist]
    else:
        f = 'blk' if cell == A.CELL else 'blk%d' % cell
        parts = [raw['sc|%s|%s|%s|%s' % (scn, pk(n, p), key, f)][CUTS.index(cut)] for n, p in plist]
    return np.concatenate(parts), [len(q) for q in parts]


def total(raw, plist, key, scn='base', cut=22.0):
    """The summed counts [true_c, matched_c, det_c, detm_c] over the patches."""
    return blocks(raw, plist, key, scn, cut)[0].sum(0)


def capped(raw, name, pxy, key):
    """True if the ISS run behind the row `key` of the patch hit the entry cap (the stored bic = [block, hist length, cap flag]); a baseline row is never capped."""
    kern = key.split('|')[0]
    if kern not in ('iso', 'adq', 'scl'):
        return False
    return bool(raw['st|%s|%s|%s|bic' % (pk(name, pxy), kern, key.split('|')[1])][2])


def defective(raw, plist, keys):
    """The patches of plist where a run behind one of the rows `keys` hit the entry cap."""
    return [(n, p) for n, p in plist if any(capped(raw, n, p, k) for k in keys)]


def clean(raw, plist, keys):
    return [np_ for np_ in plist if np_ not in defective(raw, plist, keys)]


def yn(flag):
    return 'yes' if flag else 'no'


def general_patches(raw):
    """The patches where the rule picks a general kernel, from the stored flags."""
    return [(n, p) for n, p in PATCHES if raw['st|%s|kernel|same' % pk(n, p)][0] == 0]


def dF_stats(raw, plist, ka, kb, scn='base', cut=22.0, cell=A.CELL, B=BOOT, seed=0):
    ba, sizes = blocks(raw, plist, ka, scn, cut, cell)
    bb, _ = blocks(raw, plist, kb, scn, cut, cell)
    return boot_dF(ba, bb, sizes, B, seed)


def F_of(raw, plist, key, scn='base', cut=22.0):
    return float(X.cpf(total(raw, plist, key, scn, cut))[2])


def vcell(d, ci):
    return '%s (%+.3f)' % (X.verdict(d, ci), d)


def sgn(v):
    return '%+.3f' % v


def own_stats(raw, sel, label, plist):
    """(number of stars, median (dx, dy), stars (k, 6)) of the patches `plist` of a data set (selection sel)."""
    st = [raw['own|%s|%s|x%d_y%d' % (sel, label, p[0], p[1])] for p in plist]
    allst = np.concatenate(st) if st else np.zeros((0, 6))
    return len(allst), (np.median(allst[:, 3:5], axis=0) if len(allst) else np.full(2, np.nan)), allst


def shift_text(raw, scn, plist):
    """The shifts (px) a scenario applied to the truth of the patches, as text."""
    return ' \\| '.join('(%+.3f, %+.3f)' % tuple(raw['shift|%s|%s' % (scn, pk(n, p))]) for n, p in plist)


def sec_header(raw):
    sha = [(k[len('meta|sha256|'):], str(raw[k][0])) for k in sorted(raw) if k.startswith('meta|sha256|')]
    L = ['# Post hoc sensitivity analyses of the held-out repeat-frame test (finished test: nothing here changes it)', '',
         'Generated by `experiments/exp_astro_repeat_sens.py` from `repeat_sens_raw.npz`; every number below is printed by the script. **POST HOC.** The held-out repeat-frame test (`repeat.md`, `repeat_plan.md`, '
         '`repeat_raw.npz`, `freeze_repeat.json`) was finished and its result seen before any analysis below was made; the controls are those of the adversarial referee of round REP '
         '(`review/REVIEW_roundREP.md`, findings 1, 3, 4, 5, 6, 9, 10), ported here as committed, tested code. None of it changes `repeat.md`, the plan or a freeze, and nothing in the project\'s methods or '
         'configuration may change on its strength. The analyses re-count the STORED detections with the truth moved, other tolerances, other depth cuts, other magnitudes or other block sizes; no path is re-run. '
         'The exact definitions were written into the script\'s docstring before the first analysis was run (the one change after the first star counts is listed in section 11). Scoring: '
         '`astro_common.score` at its default (0.5 px, 0.5 mag, one to one, truth within 2.5 px of the edge not scored, F of completeness and purity at mag < 22), re-implemented with the tolerances as '
         'arguments (`rescore`); at the default it reproduces all %d stored rows of `repeat_raw.npz` (counts, 25-block arrays, median offsets) exactly. Intervals: 95 %% percentile intervals of a block '
         'bootstrap over the 20 x 20 px blocks (%d resamples, seed 0, stratified within each patch, the same resamples for every row of a comparison, so differences are paired), conditional on the '
         'patches; other settings are labelled. Verdicts (win / tie / loss: the paired interval is above, around, below 0) are descriptive: no hypothesis was declared for them and no multiplicity '
         'correction is made. A row marked "defective path: yes" (or D) contains an a priori default run of the isotropic kernel that hit the entry cap, an enter/exit limit cycle of `imaging.sliding_path` '
         '(section 8): it is never pooled silently.' % (len(PATCHES) * (len(SC_KEYS) + len(BASE_EXTRA)), X.BOOT), '']
    L += ['Input files (sha256 at the time of the run):', '', X.table(['file', 'sha256'], [['`%s`' % n, '`%s`' % h] for n, h in sha]), '']
    L += ['The scan offsets are read from `%s`.' % str(raw['meta|scan_key'][0]), '']
    return L


def all_labels():
    """[(label, patches)] of the data sets whose offsets are measured: the repeat frames, M2 (test windows), M53, M92."""
    return [(n, list(R.REPEAT[n]['patches'])) for n in FRAMES] + [('M2', list(X.TESTS))] + [(lab, list(X.SECOND[ds]['patches'])) for ds, lab in (('m53', 'M53'), ('m92', 'M92'))]


def xy(v):
    return '(%+.3f, %+.3f)' % (v[0], v[1])


def sec_offsets(raw):
    L = ['## 1. Truth-to-image offsets', '',
         'The offset (px) is ADDED to the truth positions to align them with the image (the sign of the spec\'s "shift to add to the truth positions"). Two sources: the rigid-scan offsets of `repeat_spec.md` '
         'section 6.2 (the masked scan with the asTrans mapping for the repeat patches, the header mapping for the M2 windows; the plain scan, which fails where HST-saturated stars set the flux scale, is '
         'listed for reference) and the template fits of this script (D1): bright isolated stars fitted with the patch\'s own PSF stamp. Selection A is the one declared before the fits (%g <= F606W < %g, '
         '%g px from every edge, no other truth source within %g px brighter than m + %g); it leaves %d to %d stars per patch of the repeat frames. Selection B is the amendment of section 11 (%g px from every edge, no '
         'other source within %g px brighter than m + %g) and the one the `own` scenarios use. Stars: kept after the fit (selected before the fit, in brackets).' % (
             SELECTIONS['A']['mlo'], SELECTIONS['A']['mhi'], SELECTIONS['A']['margin'], SELECTIONS['A']['riso'], SELECTIONS['A']['dmag_iso'],
             min(len(raw['own|A|%s|x%d_y%d' % (n, p[0], p[1])]) for n, p in PATCHES), max(len(raw['own|A|%s|x%d_y%d' % (n, p[0], p[1])]) for n, p in PATCHES),
             SELECTIONS['B']['margin'], SELECTIONS['B']['riso'], SELECTIONS['B']['dmag_iso']), '']
    rows = []
    for label, plist in all_labels():
        for pxy in plist:
            p = '%s|x%d_y%d' % (label, pxy[0], pxy[1])
            cells = [label, '(%d, %d)' % pxy, xy(raw['scan|' + p]) if 'scan|' + p in raw else '-', xy(raw['plain|' + p]) if 'plain|' + p in raw else '-']
            for sel in ('A', 'B'):
                st = raw['own|%s|%s' % (sel, p)]
                cells += ['%d (%d)' % (len(st), raw['own|%s|%s|nsel' % (sel, p)][0]), xy(np.median(st[:, 3:5], axis=0)) if len(st) else '-']
            st = raw['own|B|' + p]
            cells.append('(%.3f, %.3f)' % (rob_sd(st[:, 3]), rob_sd(st[:, 4])) if len(st) > 1 else '-')
            rows.append(cells)
    L += [X.table(['data set', 'patch', 'scan (masked)', 'scan (plain)', 'A: stars', 'A: median (dx, dy)', 'B: stars', 'B: median (dx, dy)', 'B: robust sd (dx, dy)'], rows), '']
    L += ['Per data set: the stars of all its patches pooled (the median of the offsets with the 95 %% percentile interval of a bootstrap over the stars, %d resamples, seed 0, dx and dy separately), against the '
          'mean of the scan offsets of its patches. The `own` scenarios apply the selection-B frame median to every patch of the repeat frame.' % X.BOOT, '']
    rows = []
    for label, plist in all_labels():
        cells = [label]
        for sel in ('A', 'B'):
            n_, med, st = own_stats(raw, sel, label, plist)
            if n_ >= 2:
                (mx, cx), (my, cy) = median_interval(st[:, 3]), median_interval(st[:, 4])
                cells += ['%d' % n_, '%+.3f [%+.3f, %+.3f]' % (mx, cx[0], cx[1]), '%+.3f [%+.3f, %+.3f]' % (my, cy[0], cy[1])]
            else:
                cells += ['%d' % n_, '-', '-']
        sc_ = [raw['scan|%s|x%d_y%d' % (label, p[0], p[1])] for p in plist if 'scan|%s|x%d_y%d' % (label, p[0], p[1]) in raw]
        if sc_:
            msc = np.mean(sc_, axis=0)
            cells.append(xy(msc))
            n_, med, _ = own_stats(raw, 'B', label, plist)
            cells.append(xy(med - msc) if n_ else '-')
        else:
            cells += ['-', '-']
        rows.append(cells)
    L += [X.table(['data set', 'A: stars', 'A: dx', 'A: dy', 'B: stars', 'B: dx', 'B: dy', 'scan: mean of the patches (dx, dy)', 'B median minus scan'], rows), '']
    bt = raw['bias|table']
    L += ['Contamination of a star fit by one neighbour (noise-free simulation with the PSF stamp of m2r2 (133, 281), FWHM %.2f px, single-star fit on the image of the stamp itself at a known position plus one neighbour; '
          'largest position error over 8 directions of the neighbour, px). Rows: distance of the neighbour (px); columns: its flux relative to the star (a neighbour 2.5 mag fainter has 0.100). Selection B rejects '
          'a neighbour closer than %g px that is brighter than m + %g; it allows fainter neighbours at any distance and any neighbour at %g px or more, so one neighbour pulls the fitted position by up to %.3f px (the '
          'flux 0.100 column at 3 px) or by up to %.3f px (an equally bright neighbour at 5 px; %.3f px at 6 px), toward the neighbour, in a direction that differs from star to star.' % (
              raw['pt|m2r2|x133_y281|info'][2], SELECTIONS['B']['riso'], SELECTIONS['B']['dmag_iso'], SELECTIONS['B']['riso'], raw['bias|table'][BIAS_D.index(3.0), BIAS_Q.index(0.1)],
              raw['bias|table'][BIAS_D.index(5.0), BIAS_Q.index(1.0)], raw['bias|table'][BIAS_D.index(6.0), BIAS_Q.index(1.0)]), '',
          X.table(['distance'] + ['flux %.3f' % q for q in BIAS_Q], [['%g px' % d] + ['%.3f' % v for v in bt[i]] for i, d in enumerate(BIAS_D)]), '']
    L += ['Alignment check: the median (detection - truth) position offset (px) of the matches brighter than 21 mag of the tuned rows, with the truth as it is (alpha 0), moved by the full scan offset and by the full '
          '`own` offset (alpha 1). A row that is centred after the shift says that the shift is the offset of that row.', '']
    rows = []
    for n, p in PATCHES:
        for key in ('iso|ours|k2', 'adq|ours|k2', 'crowdsource (tuned)'):
            rows.append([ptag(n, p), KEY_LABEL[key] + (' (the same run)' if key == 'adq|ours|k2' and raw['st|%s|kernel|same' % pk(n, p)][0] == 1 else '')] +
                        [xy(raw['sc|%s|%s|%s|off' % (scn, pk(n, p), key)]) for scn in ('base', 'scan100', 'own100')])
    L += [X.table(['patch', 'catalogue', 'alpha 0', 'scan, alpha 1', 'own, alpha 1'], rows), '']
    gen = general_patches(raw)
    L += ['On %s the median (detection - truth) offsets of the bright matches of tuned crowdsource, a catalogue that has nothing to do with the ISS kernels, are %s against the unshifted truth: the displacement is '
          'shared by catalogues of unrelated methods, as a property of the truth-to-image mapping of the frame would be and a property of a kernel would not. Whether it is the frame\'s astrometry, the proper motion '
          'over the epoch difference or both is not separated here.' % (
              ' and '.join(ptag(*g) for g in gen), ' and '.join(xy(raw['sc|base|%s|crowdsource (tuned)|off' % pk(*g)]) for g in gen)), '']
    return L


def scenario_list(include_ownp=True):
    """[(source, alpha, scenario)] of the shift scenarios in the order of the report: alpha 0 (the stored test), then scan and own by alpha, then ownp at alpha 1."""
    out = [('stored', 0.0, 'base')]
    for src in ('scan', 'own'):
        out += [(src, a, scn_name(src, a)) for a in ALPHAS[1:]]
    return out + ([('ownp', 1.0, 'ownp100')] if include_ownp else [])


SRC_LABEL = {'stored': 'none (the stored test)', 'scan': 'scan offsets of each patch', 'own': 'own: frame median of the star fits', 'ownp': 'own: median of each patch'}


def alphas_where(rows, cond):
    return ', '.join('%g' % a for a, d, lo, hi in rows if cond(d, lo, hi)) or 'none'


def sec_p1(raw):
    gen = general_patches(raw)
    assert gen == GEN_PATCHES, gen
    L = ['## 2. (a) P1 against the truth offset', '',
         'P1 is the pooled paired difference dF = F(chosen kernel) - F(isotropic kernel) on the patches where the rule picks a general kernel (%s), 20 px blocks, stratified bootstrap (%d resamples, seed 0). '
         'Here the HUGS truth is moved by alpha times the offset of section 1 before the matching (the scored set, the denominators and the blocks are those of the unshifted truth; the detections are never moved), '
         'for alpha in %s; alpha 0 is the stored test. Both kernels, and every row of the table, are scored against the same moved truth.' % (', '.join(ptag(*g) for g in gen), X.BOOT, ', '.join('%g' % a for a in ALPHAS)), '']
    summary = {}
    for src in ('scan', 'own'):
        rows, cells = [], []
        for s_, a, scn in [x for x in scenario_list() if x[0] in ('stored', src) or (src == 'own' and x[0] == 'ownp')]:
            if s_ == 'ownp' and src == 'scan':
                continue
            ba, sizes = blocks(raw, gen, 'adq|ours|k2', scn)
            bb, _ = blocks(raw, gen, 'iso|ours|k2', scn)
            d, ci = boot_dF(ba, bb, sizes)
            per = [dF_stats(raw, [g], 'adq|ours|k2', 'iso|ours|k2', scn) for g in gen]
            rows.append([('%g' % a) if s_ != 'ownp' else '1 (per-patch medians)', shift_text(raw, scn, gen), '%.3f' % X.cpf(bb.sum(0))[2], '%.3f' % X.cpf(ba.sum(0))[2], X.fmt_iv(d, ci), X.verdict(d, ci)] +
                        [X.fmt_iv(*q) for q in per])
            if s_ != 'ownp':
                cells.append((a, d, ci[0], ci[1]))
        summary[src] = cells
        L += ['### 2.%d Tuned configuration (kappa 2), source: %s\n' % (1 if src == 'scan' else 2, SRC_LABEL[src]),
              X.table(['alpha', 'shift of the truth (px): %s' % ' \\| '.join(ptag(*g) for g in gen), 'F isotropic', 'F chosen kernel', 'pooled dF', 'verdict'] + ['dF %s' % ptag(*g) for g in gen], rows), '']
    L += ['### 2.3 Intervals of the pooled dF: more resamples, other seeds\n',
          'The pooled interval of the report (%d resamples, seed 0) against %d resamples (seed 0) and, over the seeds 1 to %d (%d resamples each), the share of seeds for which the lower end is above 0 '
          'and the mean (minimum, maximum) of the lower end. The BIC readout is stored as counts and blocks only: alpha 0.' % (X.BOOT, BIG_B, len(SEEDS), X.BOOT), '']
    rows = []
    for ro_label, ka, kb in READOUTS:
        for s_, a, scn in scenario_list():
            if ka.endswith('|k0') and scn != 'base':
                continue
            cl = clean(raw, gen, [ka, kb])
            for pl, extra in [(gen, '')] + ([(cl, ', without the defective patch (%s)' % ', '.join(ptag(*g) for g in cl))] if (scn == 'base' and cl and cl != gen) else []):
                ba, sizes = blocks(raw, pl, ka, scn)
                bb, _ = blocks(raw, pl, kb, scn)
                d, ci = boot_dF(ba, bb, sizes)
                ci4 = boot_dF(ba, bb, sizes, BIG_B, 0)[1]
                lows = seed_lowers(ba, bb, sizes)
                rows.append([ro_label, SRC_LABEL[s_] + ('' if s_ == 'stored' else ', alpha %g' % a) + extra, yn(defective(raw, pl, [ka, kb])), '%.4f' % d, '[%.4f, %.4f]' % tuple(ci), '[%.4f, %.4f]' % tuple(ci4),
                             '%d of %d' % ((lows > 0).sum(), len(lows)), '%.4f (%.4f, %.4f)' % (lows.mean(), lows.min(), lows.max())])
    L += [X.table(['readout', 'truth shift', 'defective path', 'dF', 'interval, %d resamples' % X.BOOT, 'interval, %d resamples' % BIG_B, 'seeds with lower end > 0', 'lower end over the seeds: mean (min, max)'], rows), '']
    L += ['### 2.4 The a priori default and the BIC stop, with and without the defective path\n',
          'The isotropic default run of m2r2 (133, 281) hit the entry cap (a defective path): every pooled row of the two m2r2 patches contains it. The rows are given for both patches (defective: yes), for '
          'the patch whose runs are both complete (defective: no) and for the defective patch alone. The BIC readout cannot be moved (alpha 0 only). dF with its interval (%d resamples, seed 0), 20 px blocks.' % X.BOOT, '']
    rows = []
    for ro_label, ka, kb in READOUTS[1:]:
        for src in ('scan', 'own'):
            if ka.endswith('|k0') and src == 'own':
                continue
            for lab, pl in (('both patches', gen), ('without the defective patch', clean(raw, gen, [ka, kb])), ('the defective patch alone', defective(raw, gen, [ka, kb]))):
                if not pl:
                    continue
                cells = []
                for a in ALPHAS:
                    scn = scn_name(src, a)
                    cells.append(X.fmt_iv(*dF_stats(raw, pl, ka, kb, scn)) if not (ka.endswith('|k0') and a > 0) else '-')
                rows.append([ro_label, (src + ' (alpha 0: the stored test)') if not ka.endswith('|k0') else 'stored', lab + ' (%s)' % ', '.join(ptag(*g) for g in pl), yn(defective(raw, pl, [ka, kb]))] + cells)
    L += [X.table(['readout', 'source', 'patches', 'defective path'] + ['alpha %g' % a for a in ALPHAS], rows), '']
    L += ['### 2.5 With the per-method magnitude offsets removed\n',
          'The median (detected - true) magnitude of the matched pairs with true 17 < m < 20 is removed per patch and method before the matching (D0), without a shift and with the truth shifted by the full scan '
          'and the full own offsets (alpha 1). Pooled over the two m2r2 patches, 20 px blocks.', '']
    rows = []
    for ro_label, ka, kb in READOUTS[:2]:
        for scn, lab in (('base', 'none'), ('mo', 'offsets removed'), ('mo+scan100', 'offsets removed, scan shift (alpha 1)'), ('mo+own100', 'offsets removed, own shift (alpha 1)'), ('scan100', 'scan shift (alpha 1), no offsets removed')):
            d, ci = dF_stats(raw, gen, ka, kb, scn)
            rows.append([ro_label, lab, yn(defective(raw, gen, [ka, kb])), X.fmt_iv(d, ci), X.verdict(d, ci)] + [X.fmt_iv(*dF_stats(raw, [g], ka, kb, scn)) for g in gen])
    L += [X.table(['readout', 'truth and magnitudes', 'defective path', 'pooled dF', 'verdict'] + ['dF %s' % ptag(*g) for g in gen], rows), '']
    sc_, ow = summary['scan'], summary['own']
    v0, v1 = X.verdict(sc_[0][1], (sc_[0][2], sc_[0][3])), X.verdict(sc_[-1][1], (sc_[-1][2], sc_[-1][3]))
    hyp = ' and '.join('%.2f' % np.hypot(*raw['scan|%s|x%d_y%d' % (n, p[0], p[1])]) for n, p in gen)
    L += ['**What section 2 shows.** For the tuned readout the pooled dF is %+.3f at alpha 0 (the stored test, verdict %s) and %+.3f at alpha 1 with the scan offsets (verdict %s; %+.3f with the own offsets). The lower end of '
          'the pooled interval is above 0 at alpha = %s with the scan offsets and at alpha = %s with the own offsets; the point estimate is negative at alpha = %s (scan) and %s (own). The patchwise dF at alpha 1 with the '
          'scan offsets are %s. %s The distance between the truth and the image on the two patches is %s px (the scan offsets). ' % (
              sc_[0][1], v0, sc_[-1][1], v1, ow[-1][1], alphas_where(sc_, lambda d, lo, hi: lo > 0), alphas_where(ow, lambda d, lo, hi: lo > 0), alphas_where(sc_, lambda d, lo, hi: d < 0),
              alphas_where(ow, lambda d, lo, hi: d < 0), ' and '.join(sgn(dF_stats(raw, [g], 'adq|ours|k2', 'iso|ours|k2', 'scan100')[0]) for g in gen),
              ('The verdict of the declared P1 rule is therefore conditional on the truth-to-image offset that the plan left uncorrected.' if v0 != v1 else 'The verdict of the declared P1 rule does not change.'), hyp) +
          '**What it does not show.** Which alignment is right: both sources estimate the NET offset, which contains the SDSS astrometry and the proper motion over the epoch difference, and alpha 1 is a sensitivity, not a '
          'corrected test. That the chosen kernel is worse than the isotropic one: a lost interval is a loss of evidence for the gain, not evidence for a loss. Why the two kernels place the same stars at positions '
          'that differ systematically on m2r2 (section 10): the cause is not identified. Anything about other frames: the rule picks a general kernel on one frame only (two patches), so these intervals contain '
          'no between-patch and no between-frame variance. The stored test and its verdict (`repeat.md`) are not revised by this analysis.', '']
    return L


def sec_blocksize(raw):
    gen = general_patches(raw)
    L = ['## 3. (b) P1 robustness: block size', '',
         'dF of the chosen against the isotropic kernel with the stored detections (alpha 0), the blocks of the bootstrap %s px wide (stratified within patch, %d resamples, seed 0). The row over all six patches '
         'is DILUTED: %d of the %d patches keep the isotropic kernel, dF = 0 by construction there, so it narrows the interval and says nothing more than the %d patches of the rule\'s general kernel do.' % (
             ', '.join('%d' % c for c in CELLS), X.BOOT, len(PATCHES) - len(gen), len(PATCHES), len(gen)), '']
    rows, lowers = [], {}
    for ro_label, ka, kb in READOUTS[:2]:
        for lab, pl in (('both m2r2 patches pooled', gen), (ptag(*gen[0]), [gen[0]]), (ptag(*gen[1]), [gen[1]]), ('all six patches (DILUTED)', PATCHES)):
            cells = []
            for cell in CELLS:
                d, ci = dF_stats(raw, pl, ka, kb, 'base', 22.0, cell)
                cells.append(X.fmt_iv(d, ci))
                if lab.startswith('both'):
                    lowers.setdefault(ro_label, []).append((cell, ci[0]))
            rows.append([ro_label, lab, yn(defective(raw, pl, [ka, kb]))] + cells)
    L += [X.table(['readout', 'patches', 'defective path'] + ['%d px (%d blocks per patch)' % (c, (N // c) ** 2) for c in CELLS], rows), '']
    L += ['**What section 3 shows.** The lower end of the pooled interval of the two m2r2 patches is above 0 at the block sizes %s px for the tuned readout and %s px for the a priori default (a row containing a '
          'defective path). The interval is stable to the block size within the two patches. **What it does not show.** Variation between patches or between frames: the two patches are windows of one frame and '
          'the intervals are conditional on them; a different block size changes the resampling unit inside the patches only.' % (
              ', '.join('%d' % c for c, lo in lowers[READOUTS[0][0]] if lo > 0) or 'none', ', '.join('%d' % c for c, lo in lowers[READOUTS[1][0]] if lo > 0) or 'none'), '']
    return L


def pair_cells(raw, scn, ka, kb, cut=22.0):
    """(per patch [(dF, interval)], pooled (dF, interval) over the six patches, [verdicts]) of the pair ka - kb in a scenario at a depth cut."""
    per = [dF_stats(raw, [np_], ka, kb, scn, cut) for np_ in PATCHES]
    return per, dF_stats(raw, PATCHES, ka, kb, scn, cut), [X.verdict(d, ci) for d, ci in per]


def pair_rows(raw, settings):
    """The rows of the P2 pair table for the settings [(label, scenario, cut)]: [setting, comparison, per-patch cells, pooled, verdict, wins / ties / losses]."""
    rows = []
    for label, scn, cut in settings:
        for name, ka, kb in P2_PAIRS:
            per, pooled, v = pair_cells(raw, scn, ka, kb, cut)
            rows.append([label, name] + [vcell(d, ci) for d, ci in per] + [X.fmt_iv(*pooled), X.verdict(*pooled), '%d / %d / %d' % (v.count('win'), v.count('tie'), v.count('loss'))])
    return rows


PAIR_HEADER = ['setting', 'comparison'] + [ptag(n, p) for n, p in PATCHES] + ['pooled dF (six patches)', 'verdict', 'wins / ties / losses']
F_KEYS = ['iso|ours|k2', 'adq|ours|k2', 'crowdsource (tuned)', 'photutils (tuned)']


def f_rows(raw, settings):
    """F of the ISS rows and the tuned baselines per patch and pooled, for the settings [(label, scenario, cut)]."""
    rows = []
    for label, scn, cut in settings:
        for key in F_KEYS:
            rows.append([label, KEY_LABEL[key]] + ['%.3f' % F_of(raw, [np_], key, scn, cut) for np_ in PATCHES] + ['%.3f' % F_of(raw, PATCHES, key, scn, cut)])
    return rows


F_HEADER = ['setting', 'catalogue'] + [ptag(n, p) for n, p in PATCHES] + ['pooled F (six patches)']
TOL_LABEL = {'base': '0.5 px / 0.5 mag (stored)', 'tol35': '0.35 px / 0.5 mag', 'mag30': '0.5 px / 0.3 mag', 'tol75': '0.75 px / 0.75 mag'}


def sec_p2(raw):
    L = ['## 4. (c) P2 robustness', '',
         'The pairs of the declared P2 analysis at the tuned configuration: the ISS with the chosen kernel and with the isotropic kernel, each against tuned crowdsource and tuned photutils, per patch (verdict by '
         'the paired interval of the patch\'s %d blocks; in brackets dF) and pooled over the six patches (stratified bootstrap, %d resamples, seed 0). The detections of the baselines are the stored catalogues '
         'through the frozen post-filter. Wins / ties / losses count the six patches.' % ((N // A.CELL) ** 2, X.BOOT), '']
    tol = [(TOL_LABEL[s], s, 22.0) for s in ('base', 'tol35', 'mag30', 'tol75')]
    L += ['### 4.1 (i) Match tolerances\n', X.table(PAIR_HEADER, pair_rows(raw, tol)), '', 'F of the catalogues under the tolerances:\n', X.table(F_HEADER, f_rows(raw, tol)), '']
    rows, cp = [], {}
    for key in F_KEYS:
        t7 = sum(raw['st|%s|%s|tot' % (pk(n, p), key)] for n, p in PATCHES)
        c_, pu_, f_ = X.cpf(t7[:4])
        cp[key] = (c_, pu_, int(t7[6] - t7[5]))
        rows.append([KEY_LABEL[key], '%d / %d' % (t7[1], t7[0]), '%.3f' % c_, '%d / %d' % (t7[3], t7[2]), '%.3f' % pu_, '%.3f' % f_, '%.3f (%d)' % (t7[5] / max(t7[6], 1), t7[6] - t7[5])])
    L += ['Completeness and purity of the catalogues with the stored scoring (0.5 px / 0.5 mag, m < 22), pooled over the six patches:\n',
          X.table(['catalogue', 'matched / truth < 22', 'completeness', 'matched / detected < 22', 'purity', 'F', 'purity over all magnitudes (unmatched detections)'], rows), '']
    dep = [('m < %g' % c, 'base', c) for c in CUTS]
    L += ['### 4.2 (ii) Depth cuts\n', 'The cut applies to the counts (truth, matched truth, detections, matched detections); the matching does not depend on it.\n', X.table(PAIR_HEADER, pair_rows(raw, dep)), '',
          'F of the catalogues at the depth cuts:\n', X.table(F_HEADER, f_rows(raw, dep)), '']
    L += ['### 4.3 (iii) The truth moved by alpha times the offsets of section 1\n',
          'Every patch by its own offset (scan: of the patch; own: the frame median), alpha in %s; the same shift for all catalogues of a patch.\n' % ', '.join('%g' % a for a in ALPHAS)]
    for src in ('scan', 'own'):
        rows_d, rows_w, rows_f = [], [], []
        for s_, a, scn in [x for x in scenario_list(False) if x[0] in ('stored', src)]:
            cells_d, cells_w = [], []
            for name, ka, kb in P2_PAIRS:
                per, pooled, v = pair_cells(raw, scn, ka, kb)
                cells_d.append(X.fmt_iv(*pooled))
                cells_w.append('%d / %d / %d, smallest %s' % (v.count('win'), v.count('tie'), v.count('loss'), sgn(min(d for d, ci in per))))
            rows_d.append(['%g' % a] + cells_d)
            rows_w.append(['%g' % a] + cells_w)
            rows_f.append(['%g' % a] + ['%.3f' % F_of(raw, PATCHES, key, scn) for key in F_KEYS])
        L += ['Source: %s. Pooled dF [95 %%] of the four pairs by alpha:\n' % SRC_LABEL[src], X.table(['alpha'] + [p[0] for p in P2_PAIRS], rows_d), '',
              'Wins / ties / losses over the six patches and the smallest patchwise dF:\n', X.table(['alpha'] + [p[0] for p in P2_PAIRS], rows_w), '',
              'Pooled F of the catalogues:\n', X.table(['alpha'] + [KEY_LABEL[k] for k in F_KEYS], rows_f), '']
    mo = [('none', 'base', 22.0), ('offsets removed', 'mo', 22.0), ('offsets removed, scan shift (alpha 1)', 'mo+scan100', 22.0), ('offsets removed, own shift (alpha 1)', 'mo+own100', 22.0)]
    L += ['### 4.4 (iv) The per-method median magnitude offset removed before the matching\n',
          'Per patch and method the median (detected - true) magnitude of the matched pairs with true 17 < m < 20 (at least %d pairs, else 0) is subtracted from the detected magnitudes (matching, depth cut and counts), '
          'without a shift and with the full scan and own shifts. The frames\' own calibrations differ from the M2 zero point used for every frame by the amount given in `repeat_spec.md` section 5.3.\n' % MAG_MIN_PAIRS,
          X.table(PAIR_HEADER, pair_rows(raw, mo)), '', 'F of the catalogues:\n', X.table(F_HEADER, f_rows(raw, mo)), '', 'The magnitude offsets removed (mag, detected - true, no shift):\n',
          X.table(['patch'] + [KEY_LABEL[k] for k in F_KEYS], [[ptag(n, p)] + ['%+.3f' % raw['sc|mo|%s|%s|moff' % (pk(n, p), k)][0] for k in F_KEYS] for n, p in PATCHES]), '']
    allset = ([(TOL_LABEL[s], s, 22.0) for s in ('tol35', 'mag30', 'tol75')] + [('m < %g' % c, 'base', c) for c in (20.0, 21.0)] + [('%s, alpha %g' % (src, a), scn, 22.0) for src in ('scan', 'own') for s_, a, scn in scenario_list(False)
              if s_ == src] + [(lab, scn, 22.0) for lab, scn, _ in mo[1:]])
    exc, ncell = [], 0
    for label, scn, cut in [('stored', 'base', 22.0)] + allset:
        for name, ka, kb in P2_PAIRS:
            per, pooled, v = pair_cells(raw, scn, ka, kb, cut)
            for (n, p), (d, ci), vv in zip(PATCHES, per, v):
                ncell += 1
                if vv != 'win':
                    exc.append([label, name, ptag(n, p), X.fmt_iv(d, ci), vv])
    L += ['### 4.5 Every cell of sections 4.1 to 4.4 and the stored scoring where the ISS does not win\n',
          'Of %d (setting, comparison, patch) cells %d are not a win.' % (ncell, len(exc)), '']
    if exc:
        L += [X.table(['setting', 'comparison', 'patch', 'dF', 'verdict'], exc), '']
    losses = ['%s, %s, %s' % (r_[0], r_[1], r_[2]) for r_ in exc if r_[4] == 'loss']
    L += ['**What section 4 shows.** Under the match tolerances, the depth cuts, the shifted truth and the magnitude offsets removed, the pooled dF of the four pairs and the number of patches won are those of the tables '
          'above; the cells that are not a win are listed in 4.5, and those where the paired interval is below 0 (a baseline ahead) are: %s. **What it does not show.** That the ISS '
          'wins on other fields or with a stronger baseline: the baselines are the tuned family of the freeze (post-filters of the output and a photutils threshold, no high-pass front end); the six patches are four '
          'sky windows of one cluster and the truth is shared with M2; F is the harmonic mean of completeness and purity: the pooled completeness of the chosen-kernel ISS is %.3f against %.3f (crowdsource) and %.3f (photutils), its purity %.3f against %.3f and %.3f, with %d unmatched '
          'detections of any magnitude against %d and %d (section 4.1).' % ('; '.join(losses) or 'none', cp['adq|ours|k2'][0], cp['crowdsource (tuned)'][0], cp['photutils (tuned)'][0], cp['adq|ours|k2'][1],
                                                                          cp['crowdsource (tuned)'][1], cp['photutils (tuned)'][1], cp['adq|ours|k2'][2], cp['crowdsource (tuned)'][2], cp['photutils (tuned)'][2]), '']
    return L


SCALED = [('tuned configuration (kappa 2)', 'adq|ours|k2', 'scl|ours|k2'), ('a priori default (kappa 5)', 'adq|default|k5', 'scl|default|k5'), ('a priori default with the BIC stop', 'adq|default|k0', 'scl|default|k0')]


def entrants(raw, name, pxy, key):
    """Number of entrants of the run behind a row key (the stored secs = [seconds, entrants])."""
    return int(raw['st|%s|%s|%s|secs' % (pk(name, pxy), key.split('|')[0], key.split('|')[1])][1])


def sec_scaled(raw):
    L = ['## 5. (d) Scaled lengths by frame', '',
         'EXPLORATORY rows of `repeat.md`, by frame: the chosen kernel with sigma_b and r_ex multiplied by the seeing (FWHM of the stamp over the mean FWHM of the stamps of the two M2 development patches), kappa unchanged, '
         'against the same kernel with the pixel-unit constants. dF = F(scaled) - F(pixel units), paired intervals (%d resamples, seed 0). The pixel-unit a priori default of the chosen kernel hit the entry cap on %s '
         '(where the rule keeps the isotropic kernel the chosen-kernel run is the isotropic run); the scaled default runs hit it on %s. A row whose pixel-unit run is defective compares a defective run with the '
         'scaled run: its dF measures the cap as much as the scaling. Kappa was not re-tuned for the scaled lengths.' % (
             X.BOOT, ', '.join(ptag(*np_) for np_ in defective(raw, PATCHES, ['adq|default|k5'])) or 'no patch', ', '.join(ptag(*np_) for np_ in defective(raw, PATCHES, ['scl|default|k5'])) or 'no patch'), '']
    frames = [(n, [np_ for np_ in PATCHES if np_[0] == n]) for n in FRAMES]
    rows = []
    for lab, ka, kb in SCALED:
        groups = frames + [('all six patches', PATCHES)]
        wo = clean(raw, PATCHES, [ka, kb])
        if len(wo) != len(PATCHES):
            groups.append(('all patches without a defective run', wo))
        for gl, pl in groups:
            if not pl:
                continue
            d, ci = dF_stats(raw, pl, kb, ka)
            rows.append([lab, gl, yn(defective(raw, pl, [ka, kb])), '%.3f' % F_of(raw, pl, ka), '%.3f' % F_of(raw, pl, kb), X.fmt_iv(d, ci), X.verdict(d, ci)])
    L += [X.table(['readout', 'patches', 'defective path', 'F pixel units (chosen kernel)', 'F scaled lengths', 'dF (scaled - pixel)', 'verdict'], rows), '']
    rows = []
    for lab, ka, kb in SCALED:
        for n, p in PATCHES:
            d, ci = dF_stats(raw, [(n, p)], kb, ka)
            rows.append([lab, ptag(n, p), '%.2f' % raw['st|%s|scale' % pk(n, p)][0], yn(defective(raw, [(n, p)], [ka, kb])), '%.3f' % F_of(raw, [(n, p)], ka), '%.3f' % F_of(raw, [(n, p)], kb), X.fmt_iv(d, ci), X.verdict(d, ci),
                         '%d / %d' % ((entrants(raw, n, p, ka), entrants(raw, n, p, kb)) if not ka.endswith('k0') else (entrants(raw, n, p, ka), entrants(raw, n, p, kb)))])
    L += ['By patch (scale: the factor s applied to sigma_b and r_ex; entrants of the pixel-unit and the scaled path):\n',
          X.table(['readout', 'patch', 'scale s', 'defective path', 'F pixel units', 'F scaled', 'dF (scaled - pixel)', 'verdict', 'entrants (pixel / scaled)'], rows), '']
    d1 = {lab: dF_stats(raw, PATCHES, kb, ka) for lab, ka, kb in SCALED}
    L += ['**What section 5 shows.** At the tuned configuration the pooled dF (scaled - pixel units) is %s, by frame %s; at the a priori default it is %s, by frame %s (the default rows of m2r3 contain the capped runs). '
          '**What it does not show.** That scaling is a method: kappa was not re-tuned for it and nothing else was; the comparison at the default is not a comparison of valid runs on m2r3; the frame is the unit of '
          'replication for the seeing and there are three of them.' % (
              X.fmt_iv(*d1[SCALED[0][0]]), ', '.join('%s %s' % (n, sgn(dF_stats(raw, pl, SCALED[0][2], SCALED[0][1])[0])) for n, pl in frames), X.fmt_iv(*d1[SCALED[1][0]]),
              ', '.join('%s %s' % (n, sgn(dF_stats(raw, pl, SCALED[1][2], SCALED[1][1])[0])) for n, pl in frames)), '']
    return L


def sec_ratios(raw):
    L = ['## 6. (e) F ratios and the seeing ladder', '',
         'F(ISS, chosen kernel, tuned configuration) over F(tuned crowdsource) and over F(tuned photutils) by frame, from the stored rows of `repeat_raw.npz` and, for M2, of `m2_raw.npz` (the four test windows, the frozen '
         'isotropic ISS at kappa 2, which is also the chosen kernel there), with 95 %% intervals of the ratio over the same stratified resamples (%d resamples, seed 0). "M2 windows" is the ratio in M2 on the windows '
         'that the frame repeats (the patches\' windows, pooled): the ratio of the same sky at the seeing of M2.' % X.BOOT, '']
    m2w = [tuple(w) for w in X.TESTS]

    def m2_blocks(windows, row):
        parts = [raw['m2|x%d_y%d|%s|blk22' % (w[0], w[1], row)] for w in windows]
        return np.concatenate(parts), [len(q) for q in parts]

    def fwhm_of(label, plist):
        return float(np.mean([raw['pt|%s|x%d_y%d|info' % (label, p[0], p[1])][2] for p in plist]))
    rows = []
    groups = [('M2 (4 test windows)', 'M2', m2w, None)] + [(n, n, [p for nn, p in PATCHES if nn == n], None) for n in FRAMES] + [('the six repeat patches', None, None, None)]
    for gl, label, plist, _ in groups:
        if label is None:
            pl = PATCHES
            ba, sizes = blocks(raw, pl, 'adq|ours|k2')
            bc, _ = blocks(raw, pl, 'crowdsource (tuned)')
            bp, _ = blocks(raw, pl, 'photutils (tuned)')
            wins = [tuple(R.REPEAT[n]['windows'][R.REPEAT[n]['patches'].index(tuple(p))]) for n, p in pl]
            fw = float(np.mean([raw['pt|%s|x%d_y%d|info' % (n, p[0], p[1])][2] for n, p in pl]))
        elif label == 'M2':
            ba, sizes = m2_blocks(plist, 'ours|k2')
            bc, _ = m2_blocks(plist, 'crowdsource (tuned)')
            bp, _ = m2_blocks(plist, 'photutils (tuned)')
            wins = None
            fw = float(np.mean([raw['pt|M2|x%d_y%d|info' % (w[0], w[1])][2] for w in plist]))
        else:
            pl = [(n, p) for n, p in PATCHES if n == label]
            ba, sizes = blocks(raw, pl, 'adq|ours|k2')
            bc, _ = blocks(raw, pl, 'crowdsource (tuned)')
            bp, _ = blocks(raw, pl, 'photutils (tuned)')
            wins = [tuple(R.REPEAT[label]['windows'][R.REPEAT[label]['patches'].index(tuple(p))]) for n, p in pl]
            fw = fwhm_of(label, plist)
        fa, fc, fp = (float(X.cpf(b.sum(0))[2]) for b in (ba, bc, bp))
        rc, rp = boot_ratio(ba, bc, sizes), boot_ratio(ba, bp, sizes)
        cells = [gl, '%.2f' % fw, '%.3f' % fa, '%.3f' % fc, X.fmt_iv(rc[0], rc[1], 2), '%.3f' % fp, X.fmt_iv(rp[0], rp[1], 2)]
        if wins:
            ma, msz = m2_blocks(wins, 'ours|k2')
            mc, _ = m2_blocks(wins, 'crowdsource (tuned)')
            mp, _ = m2_blocks(wins, 'photutils (tuned)')
            mrc, mrp = boot_ratio(ma, mc, msz), boot_ratio(ma, mp, msz)
            cells += [X.fmt_iv(mrc[0], mrc[1], 2), X.fmt_iv(mrp[0], mrp[1], 2)]
        else:
            cells += ['-', '-']
        rows.append(cells)
    L += [X.table(['frame', 'stamp FWHM (px)', 'F ISS', 'F crowdsource', 'ratio ISS / crowdsource', 'F photutils', 'ratio ISS / photutils', 'M2 windows: ratio ISS / crowdsource', 'M2 windows: ratio ISS / photutils'], rows), '']
    L += ['The seeing ladder as `repeat.md` prints it: F at mag < 22 of the same M2 window in the frames that contain it (tuned configuration; ISS with the kernel chosen from the stamp, the isotropic kernel in brackets where '
          'it differs; tuned baselines), with the stamp FWHM of the patch and the ratio ISS / crowdsource.', '']
    rows = []
    for w in m2w:
        cells = []
        fa, fc, fp = (float(X.cpf(raw['m2|x%d_y%d|%s|blk22' % (w[0], w[1], r_)].sum(0))[2]) for r_ in M2_ROWS)
        cells.append('M2 (%.2f px): ISS %.3f, crowdsource %.3f, photutils %.3f, ratio %.2f' % (raw['pt|M2|x%d_y%d|info' % w][2], fa, fc, fp, fa / fc))
        for n, p in PATCHES:
            if tuple(R.REPEAT[n]['windows'][R.REPEAT[n]['patches'].index(tuple(p))]) != w:
                continue
            fa, fi, fc, fp = (F_of(raw, [(n, p)], k) for k in ('adq|ours|k2', 'iso|ours|k2', 'crowdsource (tuned)', 'photutils (tuned)'))
            same = raw['st|%s|kernel|same' % pk(n, p)][0] == 1
            cells.append('%s (%.2f px): ISS %.3f%s, crowdsource %.3f, photutils %.3f, ratio %.2f' % (n, raw['pt|%s|x%d_y%d|info' % (n, p[0], p[1])][2], fa, '' if same else ' (%.3f)' % fi, fc, fp, fa / fc))
        rows.append(['(%d, %d)' % w] + cells + ['-'] * (4 - len(cells)))
    L += [X.table(['M2 window', 'frame 1', 'frame 2', 'frame 3', 'frame 4'], rows), '']
    L += ['**What section 6 shows.** The ratio of the ISS to tuned crowdsource and to tuned photutils by frame, against the ratio in M2 on the same windows; a ratio that fell with the seeing would show as a decrease '
          'along the frames in order of FWHM. **What it does not show.** That the F level of a frame is controlled for the seeing: the cross-frame F depends on the window, the truth offset (section 1), the zero point '
          '(each frame is scored with the M2 zero point) and the epoch; the same-seeing control m2r1 (FWHM close to M2\'s) differs from M2 on the same windows (the F columns).', '']
    return L


DEFAULT_PAIRS = [('%s, %s vs %s' % (kl, rl, bl), '%s|default|%s' % (kern, kk), bk)
                 for kern, kl in (('adq', 'chosen kernel'), ('iso', 'isotropic kernel')) for kk, rl in (('k5', 'a priori default'), ('k0', 'BIC stop'))
                 for bk, bl in (('photutils (default)', 'photutils (default)'), ('crowdsource (default)', 'crowdsource (default)'))]


def sec_plan_items(raw):
    gen = general_patches(raw)
    L = ['## 7. (f) Plan items missing from `repeat.md`', '',
         '### 7.1 The default and BIC readouts against the baselines\' DEFAULT catalogues\n',
         'The plan asks for the a priori default and the BIC stop of the ISS against the baselines\' defaults; `repeat.md` has the crowdsource-default pairs only. Per patch: verdict of the paired interval and dF; %s. '
         'Pooled over the six patches (%d resamples, seed 0, stratified): over all six (the defective flag says whether a compared ISS run hit the entry cap) and over the patches without a defective run. The BIC readout '
         'is read from its stored blocks.' % (DEFECT_NOTE, X.BOOT), '']
    rows = []
    for name, ka, kb in DEFAULT_PAIRS:
        cells, vs = [], []
        for np_ in PATCHES:
            d, ci = dF_stats(raw, [np_], ka, kb)
            cells.append(vcell(d, ci) + (' D' if defective(raw, [np_], [ka, kb]) else ''))
        d, ci = dF_stats(raw, PATCHES, ka, kb)
        cl = clean(raw, PATCHES, [ka, kb])
        dc, cic = dF_stats(raw, cl, ka, kb) if cl else (float('nan'), (float('nan'), float('nan')))
        rows.append([name] + cells + [X.fmt_iv(d, ci), X.verdict(d, ci), yn(defective(raw, PATCHES, [ka, kb])), '%d' % len(cl), X.fmt_iv(dc, cic) if cl else '-', X.verdict(dc, cic) if cl else '-'])
    L += [X.table(['comparison'] + [ptag(n, p) for n, p in PATCHES] + ['pooled, all six', 'verdict', 'defective path', 'patches without a defective run', 'pooled without the defective runs', 'verdict'], rows), '']
    rows = []
    for name, ka, kb in DEFAULT_PAIRS:
        if not name.startswith('chosen kernel'):
            continue
        for fn in FRAMES:
            pl = [np_ for np_ in PATCHES if np_[0] == fn]
            d, ci = dF_stats(raw, pl, ka, kb)
            rows.append([name, fn, yn(defective(raw, pl, [ka, kb])), X.fmt_iv(d, ci), X.verdict(d, ci)])
    L += ['By frame (the chosen kernel; the frame\'s two patches pooled):\n', X.table(['comparison', 'frame', 'defective path', 'dF', 'verdict'], rows), '']
    L += ['### 7.2 P3: the residual scale, a verdict per readout\n',
          'The plan\'s rule: if the robust standard deviation of the standardized residual (1.4826 MAD; 1 where the noise model and the kernel are right) of the isotropic kernel is not above that of the chosen kernel on '
          'the patches where the rule picks a general kernel, the diagnosis (K) is not supported on these frames. The residual scale is stored for the tuned configuration (kappa 2) and the a priori default '
          '(kappa 5); the BIC readout has none. Beside it the patchwise P1 dF of the same readout. %s' % ('A readout whose isotropic run hit the entry cap has a residual scale of a defective catalogue.'), '']
    rows, verd = [], []
    for ro_label, cfg, ka, kb in (('tuned configuration (kappa 2)', 'ours', 'adq|ours|k2', 'iso|ours|k2'), ('a priori default (kappa 5)', 'default', 'adq|default|k5', 'iso|default|k5')):
        above = []
        for n, p in gen:
            ri, ra = float(raw['st|%s|iso|%s|rsd' % (pk(n, p), cfg)][0]), float(raw['st|%s|adq|%s|rsd' % (pk(n, p), cfg)][0])
            d, ci = dF_stats(raw, [(n, p)], ka, kb)
            above.append((ri > ra, defective(raw, [(n, p)], [ka, kb])))
            rows.append([ro_label, ptag(n, p), '%.2f' % ri, '%.2f' % ra, 'yes' if ri > ra else 'no', yn(defective(raw, [(n, p)], [ka, kb])), vcell(d, ci)])
        k = sum(1 for a_, dfc in above if a_)
        kc = sum(1 for a_, dfc in above if a_ and not dfc)
        nc = sum(1 for a_, dfc in above if not dfc)
        verd.append([ro_label, '%d of %d' % (k, len(above)), 'supported' if k == len(above) else 'not supported', '%d of %d' % (kc, nc) if nc else '-',
                     ('supported' if kc == nc else 'not supported') if nc else '-'])
    L += [X.table(['readout', 'patch', 'isotropic kernel', 'chosen kernel', 'isotropic above', 'defective path', 'P1 dF of the patch'], rows), '',
          X.table(['readout', 'patches with the isotropic value above', 'P3 verdict (the plan\'s rule)', 'among the patches without a defective run', 'verdict without the defective runs'], verd), '']
    L += ['### 7.3 Stars brighter than 19.5 mag matched\n',
          'The descriptive item of the plan: matched / in the truth for the first nine magnitude bins (15.0 to 19.5) of the stored rows (match 0.5 px / 0.5 mag), per patch and catalogue; pooled by frame and over the six patches. '
          '%s' % DEFECT_NOTE, '']
    keys = ['iso|ours|k2', 'adq|ours|k2', 'iso|default|k5', 'adq|default|k5', 'crowdsource (tuned)', 'crowdsource (default)', 'photutils (tuned)', 'photutils (default)']

    def cell(pl, k):
        v = sum(raw['st|%s|%s|bright' % (pk(n, p), k)] for n, p in pl)
        return '%d / %d%s' % (v[0], v[1], ' D' if defective(raw, pl, [k]) else '')
    rows = [[ptag(n, p)] + [cell([(n, p)], k) for k in keys] for n, p in PATCHES]
    rows += [['frame %s' % fn] + [cell([np_ for np_ in PATCHES if np_[0] == fn], k) for k in keys] for fn in FRAMES]
    rows.append(['the six patches'] + [cell(PATCHES, k) for k in keys])
    rows.append(['the patches without a defective run (their number in brackets)'] + [(lambda cl_: '%d / %d (%d)' % (*sum(raw['st|%s|%s|bright' % (pk(n, p), k)] for n, p in cl_), len(cl_)) if cl_ else '-')(clean(raw, PATCHES, [k])) for k in keys])
    L += [X.table(['patch'] + [KEY_LABEL[k] for k in keys], rows), '']
    L += ['**What section 7 shows.** The two missing pairs of the plan and the P3 verdict per readout, with the defective runs marked and the pooled rows given with and without them. **What it does not show.** A comparison of '
          'valid runs for the isotropic default on %d patches: those runs are defective, and a guarded re-run (another script\'s) is the way to replace them, not this script; P3 for the BIC readout (no stored '
          'residual); that the diagnosis (K) holds or fails beyond these frames.' % len(defective(raw, PATCHES, ['iso|default|k5'])), '']
    return L


def sec_capped(raw):
    L = ['## 8. Defective (capped) paths', '',
         'The flag is the third element of the stored `bic = [block, hist length, cap flag]` of each run. A flagged run stopped at the entry cap of %d entrants; the referee of round REP showed for m2r3 (110, 501) that '
         'this is an enter/exit limit cycle of one lattice point in `imaging.sliding_path` (`review/REVIEW_roundREP.md`, finding 2), not a censored path; the readouts of such a run are not the stated procedure. The runs '
         'below are the stored ones; none is re-run here. Where the rule keeps the isotropic kernel the chosen-kernel run is the isotropic run (the same flag).' % X.MAX_ENTRIES, '']
    rows = []
    for n, p in PATCHES:
        for kern in ('iso', 'adq', 'scl'):
            if kern == 'adq' and raw['st|%s|kernel|same' % pk(n, p)][0] == 1:
                continue
            for cfg in ('ours', 'default'):
                b, s_ = raw['st|%s|%s|%s|bic' % (pk(n, p), kern, cfg)], raw['st|%s|%s|%s|secs' % (pk(n, p), kern, cfg)]
                rows.append([ptag(n, p), {'iso': 'isotropic', 'adq': 'chosen kernel', 'scl': 'chosen kernel, scaled lengths'}[kern], cfg, '%d' % s_[1], '%d of %d' % (b[0], b[1]), 'capped' if b[2] else '-', yn(b[2])])
    L += [X.table(['patch', 'kernel', 'configuration', 'entrants', 'BIC block', 'cap', 'defective path'], rows), '']
    flagged = [(n, p, kern, cfg) for n, p in PATCHES for kern in ('iso', 'adq', 'scl') for cfg in ('ours', 'default') if raw['st|%s|%s|%s|bic' % (pk(n, p), kern, cfg)][2]]
    L += ['Flagged runs: %s. The analysis rows that contain them are marked "defective path: yes" (or D) in sections 2.3, 2.4, 3, 5 and 7, and every pooled row of those sections is also given without them.' % (
        '; '.join('%s %s %s' % (ptag(n, p), kern, cfg) for n, p, kern, cfg in flagged) or 'none'), '']
    return L


def gkey(raw, name, pxy, key, field='blk22'):
    """The key of the guarded row (field blk22 or tot) of the run behind a row key 'kern|cfg|kappa' of a patch, or None if the raw has none; the chosen-kernel run of a patch where the rule keeps the isotropic
    kernel is the isotropic run."""
    if key.count('|') != 2 or key.split('|')[0] not in ('iso', 'adq', 'scl'):
        return None
    kern, cfg, k = key.split('|')
    if kern == 'adq' and raw['st|%s|kernel|same' % pk(name, pxy)][0] == 1:
        kern = 'iso'
    g = 'gr|%s|%s|%s|%s|%s' % (pk(name, pxy), kern, cfg, k, field)
    return g if g in raw else None


def guarded_blocks(raw, plist, key):
    """(blocks, sizes) of a row over the patches with the guarded row in place of the stored one where the raw has it."""
    parts = []
    for n, p in plist:
        g = gkey(raw, n, p, key)
        parts.append(raw[g] if g else raw['st|%s|%s|blk22' % (pk(n, p), key)])
    return np.concatenate(parts), [len(q) for q in parts]


def still_defective(raw, plist, keys):
    """The patches where a run behind one of the rows is capped and not replaced by a guarded run that is itself not capped."""
    out = []
    for n, p in plist:
        for k in keys:
            if not capped(raw, n, p, k):
                continue
            kern, cfg = k.split('|')[:2]
            kern = 'iso' if kern == 'adq' and raw['st|%s|kernel|same' % pk(n, p)][0] == 1 else kern
            b = raw.get('gr|%s|%s|%s|bic' % (pk(n, p), kern, cfg))
            if gkey(raw, n, p, k) is None or b is None or b[2]:
                out.append((n, p))
                break
    return out


def sec_guarded(raw):
    """Section 8.1, only if the raw holds guarded re-runs (`--cyclefix`): the capped runs beside their guarded re-runs, and P1 and the default pairs at alpha 0 with the guarded rows in place of the defective ones."""
    runs = sorted({tuple(k[3:].split('|')[:1] + k[3:].split('|')[1:4]) for k in raw if k.startswith('gr|') and k.endswith('|tot')})
    if not runs:
        return []
    L = ['### 8.1 With the guarded re-runs', '',
         'The raw of `exp_astro_cyclefix.py` (read through the optional hook `--cyclefix`, sha256 above) holds guarded re-runs of the capped paths. Beside them: the stored run, the guarded run, and the readouts at '
         'kappa 5 and with the BIC stop (F of the patch, 20 px blocks, m < 22). The rows below replace a stored row by the guarded row of the same run where there is one, at alpha 0 (the guarded detections are not '
         'stored here, so no shift is applied); a pooled row is "defective: yes" while a capped run remains in it.', '']
    rows = []
    for (name, pxy, kern, cfg) in sorted({(k.split('|')[1], tuple(int(v) for v in k.split('|')[2][1:].split('_y')), k.split('|')[3], k.split('|')[4]) for k in raw if k.startswith('gr|') and k.endswith('|tot')}):
        key_r = 'gr|%s|%s|%s|' % (pk(name, pxy), kern, cfg)
        skey = '%s|%s|' % (kern, cfg)
        cells = [ptag(name, pxy), kern, cfg, '%d / %d' % (raw['st|%s|%s|%s|secs' % (pk(name, pxy), kern, cfg)][1], raw[key_r + 'secs'][1]) if key_r + 'secs' in raw else '-',
                 '%d / %d' % (raw['st|%s|%s|%s|bic' % (pk(name, pxy), kern, cfg)][2], raw[key_r + 'bic'][2]) if key_r + 'bic' in raw else '-']
        for k in ('k5', 'k0'):
            if cfg == 'default':
                cells += ['%.3f' % X.cpf(raw['st|%s|%s|blk22' % (pk(name, pxy), skey + k)].sum(0))[2], ('%.3f' % X.cpf(raw[key_r + k + '|blk22'].sum(0))[2]) if key_r + k + '|blk22' in raw else '-']
            else:
                cells += ['-', '-']
        rows.append(cells)
    L += [X.table(['patch', 'kernel', 'configuration', 'entrants (stored / guarded)', 'cap flag (stored / guarded)', 'F kappa 5, stored', 'F kappa 5, guarded', 'F BIC, stored', 'F BIC, guarded'], rows), '']
    gen = general_patches(raw)
    rows = []
    for ro_label, ka, kb in READOUTS[1:]:
        for lab, fn in (('stored rows', lambda pl, k: blocks(raw, pl, k)), ('guarded rows in place of the capped ones', lambda pl, k: guarded_blocks(raw, pl, k))):
            ba, sizes = fn(gen, ka)
            bb, _ = fn(gen, kb)
            d, ci = boot_dF(ba, bb, sizes)
            rows.append([ro_label, lab, yn(defective(raw, gen, [ka, kb]) if lab.startswith('stored') else still_defective(raw, gen, [ka, kb])), X.fmt_iv(d, ci), X.verdict(d, ci)] +
                        [X.fmt_iv(*boot_dF(*(fn([g], k)[0] for k in (ka, kb)), [25])) for g in gen])
    L += ['P1 at alpha 0 with the guarded rows (dF of the chosen against the isotropic kernel):\n', X.table(['readout', 'rows', 'defective path'] + ['pooled dF', 'verdict'] + ['dF %s' % ptag(*g) for g in gen], rows), '']
    rows = []
    for name, ka, kb in DEFAULT_PAIRS:
        cells = [name]
        for lab, fn, dfl in (('stored', lambda pl, k: blocks(raw, pl, k), lambda: defective(raw, PATCHES, [ka, kb])),
                             ('guarded', lambda pl, k: guarded_blocks(raw, pl, k), lambda: still_defective(raw, PATCHES, [ka, kb]))):
            ba, sizes = fn(PATCHES, ka)
            bb, _ = fn(PATCHES, kb)
            d, ci = boot_dF(ba, bb, sizes)
            cells += [X.fmt_iv(d, ci), X.verdict(d, ci), yn(dfl())]
        rows.append(cells)
    L += ['The default pairs against the baselines\' defaults (section 7.1), pooled over the six patches:\n',
          X.table(['comparison', 'stored: dF', 'verdict', 'defective path', 'guarded: dF', 'verdict', 'defective path'], rows), '']
    return L


def frame_id(label):
    """(run, camcol, field) of the SDSS frame of a data set."""
    if label in R.REPEAT:
        return tuple(R.REPEAT[label]['frame'])
    if label == 'M2':
        return tuple(R.M2_FRAME)
    ds = {'M53': 'm53', 'M92': 'm92'}[label]
    run, camcol, field = X.SECOND[ds]['prefix'].replace('patch_r_', '').split('_')
    return int(run), int(camcol), int(field)


def stamp_info(raw, label, pxy):
    """dict(gain, dark, fwhm, rms_iso, rms, ellipticity, kind, parts) of a patch from the raw."""
    p = '%s|x%d_y%d' % (label, pxy[0], pxy[1])
    info = raw['pt|%s|info' % p]
    if label in R.REPEAT:
        s = 'st|%s|' % pk(label, pxy)
        rms, shape, kind, parts = raw[s + 'rms'], raw[s + 'shape'], str(raw[s + 'kind'][0]), int(raw[s + 'parts'][0])
    else:
        s = 'pt|%s|' % p
        rms, shape, kind, parts = raw[s + 'rms'], raw[s + 'shape'], str(raw[s + 'kind'][0]), int(raw[s + 'parts'][0])
    return dict(gain=float(info[0]), dark=float(info[1]), fwhm=float(info[2]), rms_iso=float(rms[0]), rms=float(rms[1]), ell=float(shape[0]), kind=kind, parts=parts)


def sec_confound(raw):
    L = ['## 9. Frames, camera columns and stamp fits (the confound of the kernel rule)', '',
         'The kernel rule keeps the isotropic fit where its rms against the stamp\'s peak is at most %g. Gain and dark variance are those of the patches (`patch_arrays.npz`), the stamp quantities those of the stored stamp fits '
         '(`psf_raw.npz` for M2, M53, M92; `repeat_raw.npz` for the repeat patches).' % E.G.RMS_GOAL, '']
    rows, byframe = [], {}
    for label, plist in all_labels():
        run, cam, fld = frame_id(label)
        for pxy in plist:
            s = stamp_info(raw, label, pxy)
            rows.append([label, '%d-%d-%d' % (run, cam, fld), '%d' % cam, '(%d, %d)' % pxy, '%.3f' % s['gain'], '%.4f' % s['dark'], '%.2f' % s['fwhm'], '%.5f' % s['rms_iso'], '%.3f' % s['ell'], s['kind'], '%d' % s['parts'],
                         '%.5f' % s['rms']])
            byframe.setdefault(label, []).append(s['kind'])
    L += [X.table(['data set', 'run-camcol-field', 'camera column', 'patch', 'gain (e-/DN)', 'dark variance', 'stamp FWHM (px)', 'isotropic fit rms (of the peak)', 'stamp ellipticity', 'kernel the rule chooses', 'parts',
                   'rms of that kernel'], rows), '']
    gen_f = [lab for lab, ks in byframe.items() if 'general' in ks]
    iso_f = [lab for lab, ks in byframe.items() if 'general' not in ks]
    cam_g, cam_i = {frame_id(f)[1] for f in gen_f}, {frame_id(f)[1] for f in iso_f}
    L += ['The data sets where the rule picks a general kernel on at least one patch: %s (camera columns %s); where it keeps the isotropic kernel on every patch: %s (camera columns %s).' % (
        ', '.join(gen_f), ', '.join('%d' % frame_id(f)[1] for f in gen_f), ', '.join(iso_f), ', '.join('%d' % frame_id(f)[1] for f in iso_f)), '']
    if len(cam_g) == 1 and not cam_g & cam_i:
        L += ['**What section 9 shows.** The kernel rule\'s action coincides with the camera column in these frames: it acts in every frame of camera column %d and in no frame of the other camera column%s, so the stamp '
              'shape cannot be separated from the camera column (gain, dark variance, saturation level, the chip\'s PSF reconstruction) with these data. **What it does not show.** That the camera column causes the '
              'elliptical stamps, or that the stamp shape causes the gain of the chosen kernel: no frame of the other camera column%s with an elliptical stamp, and no frame of this camera column with a round one, is in '
              'the sample.' % (list(cam_g)[0], 's' if len(cam_i) > 1 else '', 's' if len(cam_i) > 1 else ''), '']
    else:
        L += ['**What section 9 shows.** The kernel rule\'s action does not coincide with a single camera column in these frames (the camera columns of the frames where it acts: %s; of the others: %s). **What it does '
              'not show.** That the stamp shape, rather than another property of the frames, causes the gain of the chosen kernel.' % (
                  ', '.join('%d' % c for c in sorted(cam_g)) or 'none', ', '.join('%d' % c for c in sorted(cam_i)) or 'none'), '']
    return L


def sec_positions(raw):
    L = ['## 10. Positions of the two kernels on m2r2', '',
         'The stored tuned catalogues (kappa 2) of the two kernels on the patches where the rule picks a general kernel. (i) The median (detection - truth) position offset of the matches brighter than 21 mag, truth '
         'as it is (stored `offset`). (ii) The median (chosen - isotropic) position difference of the flux-consistent unit pairs: mutual nearest neighbours closer than 0.4 px whose flux ratio is in [0.8, 1.25], '
         'with detected magnitude below the limit; pairs in brackets.', '']
    rows = []
    for n, p in general_patches(raw):
        cells = [ptag(n, p), xy(raw['st|%s|iso|ours|k2|offset' % pk(n, p)]), xy(raw['st|%s|adq|ours|k2|offset' % pk(n, p)])]
        for m in (19, 20, 21, 22):
            v = raw['pos|%s|m%d' % (pk(n, p), m)]
            cells.append('%s (%d)' % (xy(v[1:]), v[0]))
        rows.append(cells)
    L += [X.table(['patch', 'offset of the isotropic run', 'offset of the chosen-kernel run', 'chosen - isotropic, m < 19', 'm < 20', 'm < 21', 'm < 22'], rows), '']
    gen = general_patches(raw)
    fi0, fi1 = F_of(raw, gen, 'iso|ours|k2', 'base'), F_of(raw, gen, 'iso|ours|k2', 'scan100')
    fa0, fa1 = F_of(raw, gen, 'adq|ours|k2', 'base'), F_of(raw, gen, 'adq|ours|k2', 'scan100')
    further = all(np.hypot(*raw['st|%s|iso|ours|k2|offset' % pk(n, p)]) > np.hypot(*raw['st|%s|adq|ours|k2|offset' % pk(n, p)]) for n, p in gen)
    L += ['**What section 10 shows.** The two kernels place the same units at positions that differ by the amounts above, and the isotropic run\'s matches sit %s from the unshifted truth than the chosen kernel\'s '
          '(the offsets in (i)); with the truth moved by the full scan offsets the pooled F of the two m2r2 patches goes from %.3f to %.3f for the isotropic run and from %.3f to %.3f for the chosen kernel (section 2.1). '
          '**What it does not show.** Why the positions differ: the cause is NOT identified. The referee\'s noise-free test (a star built from the stamp is recovered by both kernels, so the difference is not a '
          'kernel-origin convention; `review/REVIEW_roundREP.md`, finding 1(b)) was not repeated here. Nor which kernel\'s position is nearer to the truth: the offset between the truth and the image is the one of '
          'section 1.' % ('further' if further else 'not further on every patch', fi0, fi1, fa0, fa1), '']
    return L


def sec_limits(raw):
    L = ['## 11. Limits, deviations, what was not done', '',
         '* **Post hoc.** The same patches, the same detections and the same truth as the finished test; the analyses were chosen after the result was seen, by a referee. They cannot confirm or revise the declared '
         'outcomes P1 to P3; they say how far those verdicts depend on a convention (the alignment of the truth, the tolerances, the depth, the block size, the zero point, the defective runs).',
         '* **One frame, two patches for P1; four sky windows, three exposures for P2.** The intervals are conditional on the patches. The truth (HUGS) and its defects (HST-saturated stars, stars without an '
         'F606W) are shared by all frames and with M2.',
         '* **The BIC readout is not re-scored.** The detections of the BIC readouts are not stored; their rows are the stored counts (alpha 0, 0.5 px / 0.5 mag, m < 22, 20 px blocks). The tuned and the a priori '
         'default readouts are re-scored from the stored detections at kappa 2 and 5.',
         '* **The offsets are estimates.** The scan offsets are the rigid shifts of the PSF-convolved truth model that minimise chi^2 (per patch); the own offsets are medians of %d to %d star fits per frame '
         '(selection B), with the neighbour contamination of section 1; neither separates the SDSS astrometry from the proper motion over the epoch difference, and neither is a measurement of the kernels\' own '
         'positional bias.' % (min(len(np.concatenate([raw['own|B|%s|x%d_y%d' % (n, p[0], p[1])] for p in R.REPEAT[n]['patches']])) for n in FRAMES),
                              max(len(np.concatenate([raw['own|B|%s|x%d_y%d' % (n, p[0], p[1])] for p in R.REPEAT[n]['patches']])) for n in FRAMES)),
         '* **Defective runs.** %d a priori default runs of the isotropic kernel are defective paths; this script marks them and gives the pooled rows with and without them, and does not replace them. The guarded '
         're-runs are another script\'s work (`exp_astro_cyclefix.py`). The optional hook `--run --cyclefix PATH` copies their rows from its raw (the key format is in the script\'s docstring); %s' % (
             len(defective(raw, PATCHES, ['iso|default|k5'])), 'it was used for this raw: section 8.1.' if str(raw['meta|cyclefix'][0]) == 'applied' else 'it was not used for this raw (no section 8.1).'),
         '* **Deviation from the declared design.** Selection A of the star fits, as declared, left %s stars on m2r1, m2r2, m2r3 (kept after the fit; %s selected); the frames below the declared minimum of %d stars '
         'are: %s. Selection B (the geometry of the referee\'s `template_fit.py`) was written into the docstring as an amendment after those NUMBERS were seen (and before any selection-B offset was computed), and the '
         '`own` scenarios use it; selection A is reported beside it. The amendment also states that a patch with fewer than %d selected stars uses the frame median in the scenario `ownp`.' % (
             ', '.join('%d' % len(np.concatenate([raw['own|A|%s|x%d_y%d' % (n, p[0], p[1])] for p in R.REPEAT[n]['patches']])) for n in FRAMES),
             ', '.join('%d' % sum(int(raw['own|A|%s|x%d_y%d|nsel' % (n, p[0], p[1])][0]) for p in R.REPEAT[n]['patches']) for n in FRAMES), STAR['min_stars'],
             ', '.join(n for n in FRAMES if len(np.concatenate([raw['own|A|%s|x%d_y%d' % (n, p[0], p[1])] for p in R.REPEAT[n]['patches']])) < STAR['min_stars']) or 'none', STAR['min_stars']),
         '* **Not done.** No path was re-run (the BIC catalogues at a shifted truth and the guarded default runs need re-runs); the post hoc M53 and M92 kernel gains were not tested for the same offset sensitivity (only '
         'their offsets are measured, section 1); proper motion, `fpM` masks and cosmic rays were not examined; the cause of the position difference of the two kernels (section 10) is open.',
         '']
    return L


def report(raw):
    """The markdown report from the raw alone."""
    L = []
    for sec in (sec_header, sec_offsets, sec_p1, sec_blocksize, sec_p2, sec_scaled, sec_ratios, sec_plan_items, sec_capped, sec_guarded, sec_confound, sec_positions, sec_limits):
        L += sec(raw)
    return '\n'.join(L).rstrip('\n') + '\n'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', action='store_true', help='re-score the stored detections, measure the star offsets, write repeat_sens_raw.npz')
    ap.add_argument('--from-raw', action='store_true', help='regenerate repeat_sens.md from repeat_sens_raw.npz')
    ap.add_argument('--cyclefix', metavar='PATH', help='with --run: also copy the guarded re-runs of the repeat paths from this cyclefix raw (optional hook; default: none)')
    A_ = ap.parse_args()
    if A_.run:
        run(cyclefix=A_.cyclefix)
    elif A_.from_raw:
        with open(MD, 'w') as f:
            f.write(report(load_raw()))
        print('wrote', MD)
    else:
        ap.error('give --run or --from-raw')


if __name__ == '__main__':
    main()
