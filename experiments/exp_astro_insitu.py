"""
An in-situ PSF for every method: does the stamp-versus-star mismatch explain the ISS's loss to crowdsource on M53 and M92? (POST HOC, label-using validation, label-free stamp)
==========================================================================================================================================================================

Why. `results/astro/baselines_hp2.md`: given the ISS's high-pass front end, crowdsource beats the frozen ISS on M53 and M92. `results/astro/diag.md` (`exp_astro_diag.py`): the SDSS psField stamp does not match the
stars there (the stack of isolated truth-selected stars against the stamp, real minus exact-stamp synthetic: central pixel +14.6 points on M53, -15.2 on M92, +0.8 on M2), and in synthetic worlds built from each
stamp exactly the ISS is level with crowdsource on M92 (-0.007 [-0.033, 0.019]) against -0.104 on the image. The hypothesis tested here, as stated by the coordinator: the loss comes from the stamp-versus-star
PSF mismatch, which all methods share in some way, and an IN-SITU PSF, built from the image itself without the HST truth, fixes it for every method. This script builds that PSF and gives it to every method
IN PLACE OF the psField stamp. It is a POST HOC study (the transfer test, the kernel study, the diagnostics and the baselines' rounds 1 and 2 are all known), labelled as such; nothing frozen changes.

THE DECLARATION. This docstring was written before any ISS, crowdsource or photutils run on an in-situ stamp and before any score of such a run. What was seen before it (development, all label-free except where
stated): (i) the stack statistics of the in-situ procedure on the two M2 DEVELOPMENT patches (110, 630), (210, 630), real and in synthetic worlds built (`exp_astro_diag.render_world`) from three planted stamps
(the M2, M53 and M92 psField stamps) on those patches' truth lists, to choose the estimator (below) among a ring-median stack, a clipped-mean stack and the flux-weighted stack with outlier rejection that was kept;
(ii) the number of selected stars and the central pixel / FWHM of the first, plain version of the stack on the ten target patches, compared with their psField stamps (label-free; it fixed the isolation radii
below, not any F); (iii) the truth-based stack numbers of `diag.md` for the target patches, which were known. No F of any method on an in-situ stamp had been computed. The numbers of the validation part
(`--stamps`) are label-using (HST truth) but are statistics of the stamp, not scores.

THE PATCHES. Targets: M53 (760, 590), (640, 680), (400, 660), (660, 420); M92 (1100, 220), (1430, 550); the four M2 TEST patches (310, 630), (310, 730), (360, 1040), (100, 1100) as the control (the in-situ
stamp must not hurt there). Development patches of M2 are not scored here (Task B uses them).

THE IN-SITU STAMP (label-free: the HST truth is never an input; inputs are `dn`, the gain, the dark variance and, for the wings, the patch's psField stamp)
  1. First-pass detection.  y = H (dn - median), H the ISS's high-pass of width 3 px (`exp_astro_baselines_hp.highpass_frame`); sm = gaussian_filter(y, 1 px); the local maxima of sm in a 5 x 5 footprint; their
     peak values pk. The noise unit is the median pixel sigma of the ISS's variance model (`exp_astro_baselines_hp.iss_weights`), as for the DAOStarFinder threshold of the baselines.
  2. Stars.  A maximum is a star if pk >= 20 sigma; the star is isolated if no other maximum of pk >= 5 sigma has pk >= 0.5 pk_star within 12 px and none has pk >= 0.1 pk_star within 5 px (a faint neighbour
     inside the core is excluded, a neighbour at the edge of the 25 x 25 cut-out only if it is comparable); it is more than 12 px from every edge of the patch; and no pixel of the 7 x 7 box around it exceeds
     30000 DN (SDSS saturates near 55000 DN; the 30000 DN cut is conservative and removes the stars whose cores are in the non-linear range). The stars are ranked by pk and the 40 brightest are kept.
  3. Cut-outs.  For each star: the 25 x 25 cut-out around the maximum pixel, minus the sigma-clipped mean (3 sigma, 5 passes) of its pixels farther than 8.5 px from the middle (a MEAN, not the median of
     `exp_astro_diag`: in the dense fields the median of the ring is below the mean of the unresolved light, which inflates the 6 px flux; the development worlds showed a bias of -15 to -20 % in the
     central pixel for the median and the per-star normalisation); moved by a cubic spline (`scipy.ndimage.shift`, order 3, edge values repeated) so that its flux centroid (`exp_astro_diag.centroid`: radius 3 px,
     6 iterations) is the middle pixel; cropped to 19 x 19. The measured centroid offset of the star from its maximum pixel is kept (the sub-pixel phase).
  4. Stack.  Initial shape: the median of the stars' cut-outs each divided by its sum within 6 px, normalised to sum 1 within 6 px. Then four iterations of: f_i = <c_i, s> / <s, s> (the fitted amplitude
     of star i); e_i = rms of (c_i - f_i s) within 6 px divided by the fitted PEAK f_i max(s) (the misfit as a fraction of the star's peak; typically 0.01 to 0.07 for the stars of the target patches); a star is rejected when f_i <= 0 or
     e_i > max(3 median(e over the kept stars), 0.10) (a neighbour or defect in the cut-out). [Amended before any F was computed: the first version divided by f_i, the flux, which made the floor 0.15 inert (a
     unit test showed it); the stamps were recomputed with the amended rule.]; the stack is
     the flux-weighted mean s = sum f_i c_i / sum f_i^2 over the kept stars, normalised to sum 1 within 6 px. At least 3 kept stars are needed; otherwise the patch has no in-situ stamp (reported, the psField
     stamp is NOT substituted silently; the run of that patch is skipped).
  5. Alignment bias.  The cubic-spline centring of an undersampled star is not the identity (noise-free, M2 stamp, random phases: the stack's central pixel -4.7 %, M53 +3.0 %, M92 -1.7 %). Three fixed-point
     steps correct it by analysis-by-synthesis, with no truth: core <- core + (s_obs - T(core)), where T(core) is the same procedure (shift by the star's measured sub-pixel offset, ring mean, centroid centring,
     crop, normalise) applied noise-free to a star of profile = the stamp built from `core`, weighted by f_i^2 / sum f^2 like the data. Crowding and noise are NOT in T (they are not corrected).
  6. The 51 x 51 stamp.  The stack is the core only (radius 9 px), so the stamp is a HYBRID: the in-situ core, scaled to carry the same flux within 6 px as the psField stamp (the in-situ stamp changes the SHAPE,
     not the aperture correction), inside radius 6.5 px; the psField stamp beyond 9 px; a linear blend between; negative pixels clipped to 0; unit sum (the flux convention of every method). The in-situ stamp
     is centred on the flux centroid of the stack (the psField stamp's own centroid offset is not removed from it: the stamp is the ISS's `fit_psf_parts` input, which assumes the middle pixel).
  SCOPES.  'patch' (PRIMARY): the stars of the patch itself. 'frame' (secondary, same analyses): the stars of all the patches of the same frame in the target list pooled (M53 four, M92 two, M2 four), the
  stamp used on each of them; it shows how much of a patch stamp's difference is the few stars (typically 4 to 10 per patch).

VALIDATION OF THE STAMP (`--stamps`; label-using where marked; seeds fixed; no F is computed)
  V1  In-situ against psField, per patch and scope: the central pixel of the centred, 6 px-normalised 19 x 19 core (`exp_astro_diag.stamp_core`) relative difference; FWHM (`exp_astro_diag.fwhm`); rms of the
      difference of the cores over the peak of the psField core; ellipticity (`exp_astro_psf.stamp_shape`). 95 % percentile intervals of a bootstrap over the SELECTED stars (1000 resamples, seed 0, stratified
      by patch for the frame scope, each resample re-estimated end to end).
  V2  (label-using) Against stars selected with the truth: the truth stars of `diag` item 1 (HUGS F606W < 20, more than 12 px from the edge, no source within 4.5 px brighter than m + 1.5, from
      `exp_astro_diag.stack_stars`, here also without a pixel above 30000 DN) go through the SAME estimator (steps 3 to 6): the 'truth-star stamp'. Reported: in-situ against the truth-star stamp and psField
      against the truth-star stamp (centre, FWHM, rms as V1), and the diag-style PLAIN median stack of the truth stars and of the in-situ stars against the psField core (the quantity of `diag.md`, with its
      alignment/crowding bias, to link this table to that one). Also how many of the in-situ stars are truth-isolated stars in diag's sense (nearest truth source within 1 px; no other truth source within 4.5 px
      brighter than m + 1.5) and their median F606W.
  V3  (label-using, for the world only) Recovery in the exact-stamp world of each patch (`exp_astro_diag.render_world`, the diag seed): the in-situ builder on the rendered frame against the planted stamp
      (centre, FWHM, rms), and the plain stack. It measures what the estimator does to a stamp that is known, in the patch's own crowding.
  V4  The mismatch in the front-end domain (label-using for the held-out stars only): the truth-isolated stars of V2 that are farther than 5 px from every in-situ star are 'held out'; each is fitted in the
      high-passed image (`highpass_frame`, 3 px) by amplitude and sub-pixel shift (least squares, weights `iss_weights`) with H applied to the shifted stamp (`highpass_stamp`), for the psField stamp and the in-situ
      stamp; reported: the median over the held-out stars of the rms of the standardized residual in the 11 x 11 box, per stamp, and the number of stars for which the in-situ stamp has the smaller rms. V4 is
      reported for the patch scope only.

THE RUNS (`--run`; every ISS path through `src/core/imaging_guard.sliding_path_guarded` as `exp_astro_cyclefix.py` runs it, the frozen run functions otherwise untouched; one process, nice, sequential)
  For every target patch, with the in-situ stamp of the 'patch' scope (then of the 'frame' scope) as the patch's `psf` (the psField stamp reads from the frozen results, below):
   iso     the ISS at the FROZEN configuration (`freeze.json`: sigma_b 3, r_ex 1.2, kappa 2, plus the BIC readout) with the isotropic kernel fitted to the stamp (`exp_astro_m2.run_path`);
   adq     the ISS with the kernel the rule `psf_general.fit_adequate` chooses from the in-situ stamp (`exp_astro_psf.run_path_psf`); where the rule keeps the isotropic kernel it is the same program on the
           same input and is not run again (shown as 'same'); the rule is NOT re-tuned (rms goal 2e-3 of the peak);
   cs_hp   crowdsource + the front-end high-pass, the variant the round-2 selection holds (`exp_astro_hp2.verify_selection2`: `cs2_hp2_t3_pm`: 2 px high-pass, threshold 3, miniter 8, maxiter 20, fewstars 20;
           the stamp goes through `highpass_stamp` exactly as there);
   phu_hp  photutils + the front-end high-pass, the round-2 selection (`phu2_hp3_all_i6_b1`);
   cs_def  crowdsource at its default, no high-pass (`exp_astro_baselines.run_crowdsource`, the default filter), a reference row (the stamp enters through `SimplePSF`).
  The psField-stamp rows come from the frozen results, READ ONLY and copied into the raw by `--ref`: iso and adq from `cyclefix_raw.npz` (the guarded runs, C| keys: `C|T|` M2, `C|T2|` M53/M92 iso, `C|P|`
  the kernel study's adq on M53/M92; on M2 the rule keeps the isotropic kernel so adq = iso), cs_hp and phu_hp from `baselines_hp2_raw.npz`, cs_def from `m2_raw.npz`. A path whose entry count reached the cap
  (`exp_astro_m2.MAX_ENTRIES`) is STALLED, marked, left out of pools and named (the guard removes the known cycle; the flag stays as a safeguard).
  Scoring: `astro_common.score`, F of completeness and purity at mag < 22, 20 x 20 px blocks.

THE ANALYSES (all declared; nothing else is read into the conclusion)
  F1  Pooled F (completeness, purity) at mag < 22 per cluster group (M2 test x4, M53 x4, M92 x2) and per patch, for every method with the psField stamp and with the in-situ stamp (patch and frame scopes).
  F2  Paired differences on the pooled blocks, 95 % percentile intervals of the block bootstrap (20 x 20 px blocks, 2000 resamples, seed 0, `exp_astro_m2.boot_idx_strat`: blocks resampled within each patch, the
      same resamples for every row): ISS iso minus cs_hp, ISS adq minus cs_hp, ISS iso / adq minus phu_hp, and minus cs_def, each with the psField stamp and with the in-situ stamp (the SAME stamp for
      both methods of a difference).
  F3  The stamp effect per method, in-situ minus psField (paired, same blocks), and the change of the ISS-minus-crowdsource gap, (gap with the in-situ stamp) minus (gap with psField) (paired).
  F4  The question, answered by the script's verdict cells and stated plainly in the report: on M53 and on M92, with the in-situ stamp for all methods, does the ISS-minus-crowdsource(+HP) gap close
      (its interval no longer excludes 0 from below, or turns positive) or reverse (excludes 0 from above, in the ISS's favour), and does the M2 control not get worse?
  The verdicts: win / loss if the 95 % interval excludes 0. The groups are windows of one frame each (not independent trials); intervals are conditional on the patches.

AMENDMENT 1 (post hoc, added after the first crowdsource runs on the patch-scope stamps and before the ISS runs; the declared analyses above are unchanged and stay the primary ones). The first runs
showed that crowdsource+HP loses a lot of F with a patch-scope in-situ stamp on M2 (pooled 0.393 with the psField stamp against 0.230) although the stamp differs from psField by a few per cent in the core, and a
quick scratch check on one patch (not stored) showed a large loss for the psField stamp plus a pixel perturbation of the same size in a random order. To separate the sensitivity of a method to ANY perturbation of
the stamp from the effect of the in-situ stamp's particular shape, a control scope 'scr' is added (F5): the psField stamp plus the pixels of (patch-scope in-situ stamp - psField stamp) permuted at random (seeded
by the patch; `scrambled_stamp`), clipped, unit sum: the amplitude distribution of the in-situ perturbation, none of its structure. It is run for the same methods (the baselines on all patches; the ISS where time
allowed: the report lists what ran) and enters the tables as a stamp column next to the psField and the in-situ columns. F5: the effect of the control stamp on each method, scr minus psField (paired); a method whose
F changes by about as much under the scrambled perturbation as under the in-situ stamp is not measuring the stamp's shape.

TASK B (secondary, DEVELOPMENT patches only; `--cs-dev`; the crowdsource options that the round-2 plan lists as not searched, `baselines_hp2_plan.md` section 7, plus others `fit_im` accepts)
  Declared set, run on the two M2 development patches (110, 630), (210, 630) with the round-2 selected variant as the base B0 (high-pass 2 px, threshold 3, miniter 8, maxiter 20, fewstars 20, psfderiv True,
  blendthreshu 2, titer_thresh 2, psfvalsharpcutfac 0.7, psfsharpsat 0.7, the ISS weights, no bright-star list, no mask; the base must reproduce the stored `cs2_hp2_t3_pm` development catalogue counts):
  stage 1, one option at a time: blendthreshu in {1, 3, 4}; titer_thresh in {1, 3, 5}; psfderiv False; maxiter in {30, 40}; miniter in {4, 12}; fewstars 5; psfvalsharpcutfac in {0.5, 0.9}; psfsharpsat in
  {0.5, 1.0}; threshold in {2.5, 4}; blist (the bright-star list of crowdsource: the 5 or the 15 brightest maxima of the step-1 detection, magnitude from the peak) ; weights: 'flat' (the constant median
  inverse noise), 'smooth3' (the ISS weights with the 3 px smoothing of the expected counts instead of 1.5 px), 'satmask' (weight 0 where dn > 30000 DN). 23 variants + the base.
  stage 2, declared now: one combination variant 'combo' = every stage-1 option whose pooled development F exceeds the base's by 0.003 or more (in the order of the list above, each family at its best level;
  none if no option reaches it).
  Selection: the largest pooled F at mag < 22 over the two development patches (counts added, rounded to 6 decimals; a tie to the base, among other tied variants the first in the declared order). Report the best against the round-2 selection. A variant that beats the
  selection by MORE THAN 0.010 pooled development F is a finding; only then (`--cs-test`, clearly post hoc, once) it is scored on the four M2 TEST patches beside the round-2 selection. A run that fails is
  recorded as failed, not retried with other settings.

    python experiments/exp_astro_insitu.py --stamps                  validation and the stamps of both scopes (minutes) -> results/astro/insitu_raw.npz
    python experiments/exp_astro_insitu.py --scr                     the control stamps of Amendment 1 (scope 'scr'), from the patch-scope stamps
    python experiments/exp_astro_insitu.py --ref                     the psField-stamp rows copied from the frozen results (read only) into the raw
    <domain python> -I experiments/exp_astro_insitu.py --run [--scope patch|frame] [--methods iso,adq,...] [--sets m2,m53,m92]   the runs (hours; one at a time)
    <domain python> -I experiments/exp_astro_insitu.py --cs-dev [--stage 1|2]     Task B on the development patches
    <domain python> -I experiments/exp_astro_insitu.py --cs-test                  only if a variant is a finding
    python experiments/exp_astro_insitu.py --from-raw                regenerates results/astro/insitu.md from insitu_raw.npz alone (byte-identically)

Record of the committed raw. `insitu_raw.npz` was produced on a shared machine (up to five other jobs) in single-process calls, one heavy process at a time, in this order: `--stamps` (twice: the first
version's rejection scale was amended, see step 4, and the stamps recomputed before any run), `--ref`, the baselines on the patch scope, `--scr` and the baselines on the control scope, the ISS (iso, adq) on the
patch scope for M53, M92 and then M2, `--run --scope patch --methods phu_hp`, Task B (`--cs-dev` stages 1 and 2), the baselines on the frame scope, and the ISS on the frame scope for M53 and M92 only (the M2
frame-scope ISS runs, the ISS and photutils on the control scope were NOT run: the report shows dashes). A container restart interrupted the last M2 patch-scope run, which was resumed (a run whose keys exist is
not run again). The seconds in the raw depend on the load of the machine (the same run took from 1.5 to 5 times longer under load). `insitu_raw.npz` is a snapshot of that state; the report regenerates from it
byte for byte.

Environment: INSITU_OUT (directory of the raw and the md; default results/astro), INSITU_CYCLEFIX_RAW / INSITU_HP2_RAW / INSITU_M2_RAW (the frozen raws read by --ref).
Raw keys: 'ST|<scope>|<set>|x<x0>_y<y0>|<field>' (stamps and validation), 'RUN|<scope>|<set>|x..|<method>|<field>' (scores in the layout of `exp_astro_m2._pack` under '.../k2' and '.../k0', '|secs', '|det', '|trace'),
'REF|<set>|x..|<method>|<field>' (psField rows), 'CSD|<tag>|x<x0>_y<y0>|<field>' and 'CSD|SEL|...' (Task B), 'meta|...' (hashes, versions). Write-once, atomic, under a lock.
"""
import argparse
import contextlib
import hashlib
import io
import json
import zlib
import os
import sys
import time

os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
from scipy.ndimage import gaussian_filter, maximum_filter
from scipy.ndimage import shift as nd_shift
from scipy.optimize import least_squares
from scipy.spatial import cKDTree

HERE = os.path.abspath(os.path.dirname(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)
import astro_common as A                  # noqa: E402
import exp_astro_baselines_hp as HP       # noqa: E402
import exp_astro_baselines_hp2 as H2      # noqa: E402
import exp_astro_cyclefix as CF           # noqa: E402
import exp_astro_diag as D                # noqa: E402
import exp_astro_hp2 as HS                # noqa: E402
import exp_astro_m2 as X                  # noqa: E402
import exp_astro_psf as PS                # noqa: E402
from src.core import imaging as IM        # noqa: E402
from src.core import psf_general as PG    # noqa: E402

ASTRO = os.path.join(REPO, 'results', 'astro')
OUT = os.environ.get('INSITU_OUT', ASTRO)
RAW, MD = os.path.join(OUT, 'insitu_raw.npz'), os.path.join(OUT, 'insitu.md')
CYCLEFIX_RAW = os.environ.get('INSITU_CYCLEFIX_RAW', os.path.join(ASTRO, 'cyclefix_raw.npz'))
HP2_RAW = os.environ.get('INSITU_HP2_RAW', os.path.join(ASTRO, 'baselines_hp2_raw.npz'))
M2_RAW = os.environ.get('INSITU_M2_RAW', X.RAW)

# ---- constants of the in-situ stamp (declared in the docstring)
SB_DET, SMOOTH, FOOT = 3.0, 1.0, 5          # first-pass detection: high-pass width, smoothing sigma, maximum footprint
T_DET, T_STAR = 5.0, 20.0                   # sigma: a neighbour counts from T_DET, a star from T_STAR
ISO_R1, ISO_F1, ISO_R2, ISO_F2 = 5.0, 0.1, 12.0, 0.5
SAT_DN, SAT_BOX = 30000.0, 3                # no pixel of the (2 SAT_BOX + 1)^2 box above SAT_DN
N_MAX, N_MIN = 40, 3
HALF, PAD, SKY_RING, NORM_RAD = D.HALF, D.PAD, D.SKY_RING, D.NORM_RAD
EDGE = HALF + PAD
CUT = 2 * EDGE + 1                          # the cut-out side (25)
CORE = 2 * HALF + 1                         # the core side (19)
KEEP_K, KEEP_FLOOR, N_STACK_IT, N_FIX_IT = 3.0, 0.10, 4, 3
R_FULL, R_OUT = 6.5, 9.0                    # the blend of the in-situ core into the psField wings
SCR_SEED = 2611
B_STAMP = 1000                              # bootstrap resamples over stars (V1)
HELD_RAD, FIT_BOX = 5.0, 5                  # V4: held-out stars farther than this from every in-situ star; the fit box is (2 FIT_BOX + 1)^2
RING = np.hypot(*np.mgrid[-EDGE:EDGE + 1, -EDGE:EDGE + 1]) > SKY_RING
DISK = D.DISK

TARGETS = [('m2', tuple(p)) for p in X.TESTS] + [('m53', tuple(p)) for p in X.SECOND['m53']['patches']] + [('m92', tuple(p)) for p in X.SECOND['m92']['patches']]
GROUPS = (('M2 test patches (control)', 'm2'), ('M53 (NGC 5024)', 'm53'), ('M92 (NGC 6341)', 'm92'))
SCOPES = ('patch', 'frame', 'scr')
METHODS = ('iso', 'adq', 'cs_hp', 'phu_hp', 'cs_def')
LABEL = {'iso': 'ISS, isotropic kernel (frozen config)', 'adq': 'ISS, kernel from the stamp (rule)', 'cs_hp': 'crowdsource + HP (round 2)', 'phu_hp': 'photutils + HP (round 2)',
         'cs_def': 'crowdsource default (no HP)'}
STAMPS = ('psf', 'patch', 'frame', 'scr')
STAMP_LABEL = {'psf': 'psField stamp', 'patch': 'in-situ stamp (patch stars)', 'frame': 'in-situ stamp (frame stars)', 'scr': 'control: psField + scrambled in-situ difference'}


def pkey(name, pxy):
    return '%s|x%d_y%d' % (name, pxy[0], pxy[1])


# ----------------------------------------------------------------------------------------------------------------------------------
# step 1 and 2: detection and the selection of stars
# ----------------------------------------------------------------------------------------------------------------------------------
def noise_unit(dn, gain, dark):
    """The median pixel sigma of the ISS's variance model (DN)."""
    return float(np.median(1.0 / HP.iss_weights(dn, gain, dark)))


def detect_maxima(dn, gain, dark):
    """(x, y, pk, noise): the local maxima (5 x 5 footprint, column x and row y of the maximum pixel) of the high-passed (3 px), 1 px-smoothed frame, their values and the noise unit."""
    y = HP.highpass_frame(dn, SB_DET)
    sm = gaussian_filter(y, SMOOTH)
    rows, cols = np.nonzero(sm == maximum_filter(sm, size=FOOT))
    return cols.astype(float), rows.astype(float), sm[rows, cols], noise_unit(dn, gain, dark)


def select_stars(x, y, pk, noise, dn, n=100, t_star=T_STAR, r1=ISO_R1, f1=ISO_F1, r2=ISO_R2, f2=ISO_F2, sat=SAT_DN, edge=EDGE, n_max=N_MAX):
    """(xs, ys, pks, idx): the maxima of at least T_DET sigma (every neighbour that counts) and the indices into them of the stars, brightest first: pk >= t_star sigma, more than `edge` px from every edge
    of the n x n patch, no other maximum with pk >= f2 pk_star within r2 px nor with pk >= f1 pk_star within r1 px (distances at most the radius), no pixel above `sat` in the 7 x 7 box; at most n_max."""
    keep = pk >= T_DET * noise
    xs, ys, pks = x[keep], y[keep], pk[keep]
    tree = cKDTree(np.c_[xs, ys])
    out = []
    for i in np.flatnonzero(pks >= t_star * noise):
        if not (xs[i] > edge and xs[i] < n - 1 - edge and ys[i] > edge and ys[i] < n - 1 - edge):
            continue
        near = [j for j in tree.query_ball_point([xs[i], ys[i]], r2) if j != i]
        if any(pks[j] >= f2 * pks[i] for j in near):
            continue
        if any(pks[j] >= f1 * pks[i] and np.hypot(xs[j] - xs[i], ys[j] - ys[i]) <= r1 for j in near):
            continue
        ix, iy = int(xs[i]), int(ys[i])
        if dn[max(iy - SAT_BOX, 0): iy + SAT_BOX + 1, max(ix - SAT_BOX, 0): ix + SAT_BOX + 1].max() > sat:
            continue
        out.append(i)
    out = np.array(out, int)
    out = out[np.argsort(-pks[out], kind='stable')][:n_max] if len(out) else out
    return xs, ys, pks, out


# ----------------------------------------------------------------------------------------------------------------------------------
# step 3: cut-outs
# ----------------------------------------------------------------------------------------------------------------------------------
def sky_clip(v, k=3.0, passes=5):
    """The sigma-clipped mean of v (k sigma about the current mean, `passes` passes)."""
    v = np.asarray(v, float)
    m = np.ones(len(v), bool)
    for _ in range(passes):
        mu, sd = v[m].mean(), v[m].std()
        m2 = np.abs(v - mu) <= k * sd
        if m2.sum() == m.sum() or m2.sum() < 3:
            break
        m = m2
    return float(v[m].mean())


def align_unnorm(sub):
    """The cut-out `sub` (CUT x CUT) with its flux centroid moved to the middle pixel (cubic spline, edge values repeated), cropped to CORE x CORE, NOT normalised; and the centroid offset (ox, oy) from the
    middle pixel before the move."""
    n = sub.shape[0]
    cx, cy = D.centroid(sub)
    sh = nd_shift(sub, ((n - 1) / 2.0 - cy, (n - 1) / 2.0 - cx), order=3, mode='nearest')[PAD:-PAD, PAD:-PAD]
    return sh, cx - (n - 1) / 2.0, cy - (n - 1) / 2.0


def star_cores(dn, x, y, idx):
    """(c (n, CORE, CORE), off (n, 2)): the aligned, ring-mean-subtracted, unnormalised cut-outs of the stars idx at (x, y) in the frame `dn`, and their centroid offsets (ox, oy)."""
    cs, off = [], []
    for i in idx:
        ix, iy = int(np.round(x[i])), int(np.round(y[i]))
        sub = dn[iy - EDGE: iy + EDGE + 1, ix - EDGE: ix + EDGE + 1].copy()
        sub -= sky_clip(sub[RING])
        c, ox, oy = align_unnorm(sub)
        cs.append(c)
        off.append((ox, oy))
    return np.array(cs).reshape(len(idx), CORE, CORE), np.array(off).reshape(len(idx), 2)


# ----------------------------------------------------------------------------------------------------------------------------------
# steps 4 to 6: stack, alignment bias, the hybrid stamp
# ----------------------------------------------------------------------------------------------------------------------------------
def robust_stack(c):
    """(s, f, e, use): the flux-weighted stack of the aligned cut-outs c (n, CORE, CORE), sum 1 within NORM_RAD px, with the fitted amplitudes f, the rms residuals e within NORM_RAD px as a fraction of the
    fitted peak f max(s), and the mask of the kept stars (see the docstring, step 4). s is NaN if fewer than N_MIN stars are kept."""
    n = len(c)
    ssum = c[:, DISK].sum(1)
    ok0 = ssum > 0
    if ok0.sum() < N_MIN:
        return np.full((CORE, CORE), np.nan), np.zeros(n), np.full(n, np.inf), np.zeros(n, bool)
    s = np.median(c[ok0] / ssum[ok0][:, None, None], 0)
    s = s / s[DISK].sum()
    use = ok0.copy()
    f, e = np.zeros(n), np.full(n, np.inf)
    for _ in range(N_STACK_IT):
        f = (c * s).sum((1, 2)) / (s * s).sum()
        r = c - f[:, None, None] * s
        e = np.sqrt((r[:, DISK] ** 2).mean(1)) / np.maximum(np.abs(f) * s.max(), 1e-12)
        thr = max(KEEP_K * float(np.median(e[use])), KEEP_FLOOR)
        use = (f > 0) & (e <= thr)
        if use.sum() < N_MIN:
            return np.full((CORE, CORE), np.nan), f, e, np.zeros(n, bool)
        s = (f[use, None, None] * c[use]).sum(0) / (f[use] ** 2).sum()
        s = s / s[DISK].sum()
    return s, f, e, use


def hybrid_stamp(core, psf):
    """The 51 x 51 stamp of a 19 x 19 core (sum 1 within NORM_RAD px, centred): the core scaled to the psField flux within NORM_RAD px inside R_FULL px, the psField stamp beyond R_OUT px, a linear blend
    between, negative pixels clipped, unit sum."""
    n = psf.shape[0]
    c = n // 2
    r = np.hypot(*np.mgrid[-c:c + 1, -c:c + 1])
    a = float(psf[r <= NORM_RAD].sum())
    big = np.zeros_like(psf)
    big[c - HALF: c + HALF + 1, c - HALF: c + HALF + 1] = a * core
    w = np.clip((R_OUT - r) / (R_OUT - R_FULL), 0.0, 1.0)
    s = np.clip(w * big + (1.0 - w) * psf, 0.0, None)
    return s / s.sum()


def forward_core(stamp, off, w):
    """T: the aligned, normalised core that stars of profile `stamp` at the centroid offsets off (n, 2) = (ox, oy) would give through the procedure of steps 3 and 4 (noise-free, no neighbours), weighted
    by w (n,) and normalised to sum 1 within NORM_RAD px."""
    c0 = stamp.shape[0] // 2
    acc = np.zeros((CORE, CORE))
    for (ox, oy), wi in zip(off, w):
        sub = nd_shift(stamp, (oy, ox), order=3, mode='constant', cval=0.0)[c0 - EDGE: c0 + EDGE + 1, c0 - EDGE: c0 + EDGE + 1]
        sub = sub - sky_clip(sub[RING])
        core, _, _ = align_unnorm(sub)
        acc += wi * core / core[DISK].sum()
    return acc / acc[DISK].sum()


def estimate(c, off, psf, n_fix=N_FIX_IT):
    """The in-situ stamp from the aligned cut-outs c (n, CORE, CORE) with centroid offsets off (n, 2) and the psField stamp for the wings: dict(stamp (or NaN), core, core_obs, f, e, use, resid (rms of the
    last fixed-point residual over the peak of the stack), ok)."""
    s_obs, f, e, use = robust_stack(c)
    if use.sum() < N_MIN:
        return dict(stamp=np.full(psf.shape, np.nan), core=np.full((CORE, CORE), np.nan), core_obs=s_obs, f=f, e=e, use=use, resid=np.nan, ok=False)
    w = f[use] ** 2 / (f[use] ** 2).sum()
    core = s_obs.copy()
    for _ in range(n_fix):
        t = forward_core(hybrid_stamp(core, psf), off[use], w)
        core = core + (s_obs - t)
        core = core / core[DISK].sum()
    stamp = hybrid_stamp(core, psf)
    t = forward_core(stamp, off[use], w)
    resid = float(np.sqrt(np.mean((s_obs - t) ** 2)) / s_obs.max())
    return dict(stamp=stamp, core=core, core_obs=s_obs, f=f, e=e, use=use, resid=resid, ok=True)


def collect(frames, psf):
    """The star table and aligned cut-outs of a list of frames (dict(dn, gain, dark)): dict(c, off, frame (index), x, y, snr (pk / noise), n_by_frame)."""
    cs, offs, fr, xs_, ys_, snr, nby = [], [], [], [], [], [], []
    for k, fm in enumerate(frames):
        x, y, pk, noise = detect_maxima(fm['dn'], fm['gain'], fm['dark'])
        xs, ys, pks, idx = select_stars(x, y, pk, noise, fm['dn'], n=fm['dn'].shape[0])
        c, off = star_cores(fm['dn'], xs, ys, idx)
        cs.append(c)
        offs.append(off)
        fr += [k] * len(idx)
        xs_.append(xs[idx])
        ys_.append(ys[idx])
        snr.append(pks[idx] / noise)
        nby.append(len(idx))
    return dict(c=np.concatenate(cs), off=np.concatenate(offs), frame=np.array(fr, int), x=np.concatenate(xs_), y=np.concatenate(ys_), snr=np.concatenate(snr), n_by_frame=np.array(nby, int))


def insitu_stamp(frames, psf):
    """The in-situ stamp of the stars of `frames` (a list of dict(dn, gain, dark)): the `estimate` dict plus the star table under 'stars'."""
    st = collect(frames, psf)
    est = estimate(st['c'], st['off'], psf) if len(st['c']) >= N_MIN else dict(stamp=np.full(psf.shape, np.nan), core=np.full((CORE, CORE), np.nan), core_obs=np.full((CORE, CORE), np.nan),
                                                                              f=np.zeros(len(st['c'])), e=np.full(len(st['c']), np.inf), use=np.zeros(len(st['c']), bool), resid=np.nan, ok=False)
    est['stars'] = st
    return est


# ----------------------------------------------------------------------------------------------------------------------------------
# comparisons of stamps
# ----------------------------------------------------------------------------------------------------------------------------------
def stamp_core(stamp):
    """The centred, NORM_RAD-normalised CORE x CORE core of a 51 x 51 stamp (`exp_astro_diag.stamp_core`)."""
    return D.stamp_core(stamp)


def compare_cores(a, b):
    """[central pixel of a over b minus 1, FWHM a, FWHM b, rms (a - b) over the peak of b] of two cores (CORE x CORE)."""
    c = CORE // 2
    return np.array([a[c, c] / b[c, c] - 1.0, D.fwhm(a), D.fwhm(b), np.sqrt(np.mean((a - b) ** 2)) / b[c, c]])


def shape_of(stamp):
    """[ellipticity, orientation (degrees)] of a stamp inside 5 px (`exp_astro_psf.stamp_shape`)."""
    return PS.stamp_shape(stamp)[:2]


def truth_stars(truth, x0, y0, dn):
    """(x, y, m, idx) of the truth stars of a patch in diag's sense (`exp_astro_diag.stack_stars` on the sources of the box grown by MARGIN) without a pixel above SAT_DN in the 7 x 7 box."""
    x, y, m = D.sources_near(truth, x0, y0, D.MARGIN)
    idx = D.stack_stars(x, y, m)
    keep = [i for i in idx if dn[max(int(np.round(y[i])) - SAT_BOX, 0): int(np.round(y[i])) + SAT_BOX + 1, max(int(np.round(x[i])) - SAT_BOX, 0): int(np.round(x[i])) + SAT_BOX + 1].max() <= SAT_DN]
    return x, y, m, np.array(keep, int)


def nearest_truth(x, y, tx, ty, tm, rad=1.0):
    """For each position (x, y) the magnitude of the nearest truth source within `rad` px (NaN if none) and the index of that source (-1 if none): by brute force."""
    mag, ind = np.full(len(x), np.nan), np.full(len(x), -1, int)
    for i in range(len(x)):
        d = np.hypot(tx - x[i], ty - y[i])
        if len(d) and d.min() <= rad:
            j = int(np.argmin(d))
            mag[i], ind[i] = tm[j], j
    return mag, ind


def truth_isolated_mask(ind, tx, ty, tm):
    """True where the truth source `ind` is isolated in diag's sense: no other truth source within 4.5 px (at most) brighter than m + 1.5 (`exp_astro_diag.STACK_RAD`, `STACK_DM`)."""
    out = np.zeros(len(ind), bool)
    for k, j in enumerate(ind):
        if j < 0:
            continue
        d = np.hypot(tx - tx[j], ty - ty[j])
        out[k] = not np.any((d <= D.STACK_RAD) & (np.arange(len(tx)) != j) & (tm < tm[j] + D.STACK_DM))
    return out


def stat_vector(est_stamp, ref_core):
    """[central pixel over the reference's minus 1, FWHM, rms difference over the reference's peak] of an in-situ stamp against a reference core; NaN if the stamp is NaN."""
    if not np.all(np.isfinite(est_stamp)):
        return np.full(3, np.nan)
    r = compare_cores(stamp_core(est_stamp), ref_core)
    return np.array([r[0], r[1], r[3]])


def star_bootstrap(st, psf_refs, psfs, sizes, B=B_STAMP, seed=0):
    """The 95 % percentile intervals (2 x 3 per reference) of `stat_vector` of the in-situ stamp over a bootstrap of the SELECTED stars (resampled within each stratum of `sizes`, each resample re-estimated end
    to end with the wing stamp psfs[k]) against each reference core of psf_refs: array (len(refs), 3, 2); and the number of resamples without an estimate."""
    idx = X.boot_idx_strat([s for s in sizes if s > 0], B, seed)
    c, off = st['c'], st['off']
    draws = np.full((B, len(psf_refs), 3), np.nan)
    for b in range(B):
        e = estimate(c[idx[b]], off[idx[b]], psfs[0])
        if e['ok']:
            for k, ref in enumerate(psf_refs):
                draws[b, k] = stat_vector(e['stamp'], ref)
    ci = np.nanpercentile(draws, [2.5, 97.5], axis=0) if np.isfinite(draws).any() else np.full((2, len(psf_refs), 3), np.nan)
    return np.transpose(ci, (1, 2, 0)), int((~np.isfinite(draws[:, 0, 0])).sum())


# ----------------------------------------------------------------------------------------------------------------------------------
# V4: held-out stars in the high-passed image
# ----------------------------------------------------------------------------------------------------------------------------------
def hp_fit_rms(dn, gain, dark, stamp, positions, sb=SB_DET):
    """For each position (x, y) of a star, [rms of the standardized residual, rms of the residual over the fitted peak of the model] in the (2 FIT_BOX + 1)^2 box of the high-passed frame after a least-squares
    fit of an amplitude and a sub-pixel shift of `stamp` (the model: H applied to the shifted stamp, `highpass_stamp`); NaN where the box is not inside the frame. Array (n, 2)."""
    yh = HP.highpass_frame(dn, sb)
    sw = HP.iss_weights(dn, gain, dark)
    n = stamp.shape[0]
    c0 = n // 2
    out = []
    for (px, py) in positions:
        ix, iy = int(np.round(px)), int(np.round(py))
        ya, xa = iy - FIT_BOX, ix - FIT_BOX
        if ya < 0 or xa < 0 or iy + FIT_BOX + 1 > dn.shape[0] or ix + FIT_BOX + 1 > dn.shape[1]:
            out.append((np.nan, np.nan))
            continue
        d = yh[ya: iy + FIT_BOX + 1, xa: ix + FIT_BOX + 1]
        w = sw[ya: iy + FIT_BOX + 1, xa: ix + FIT_BOX + 1]

        def model(p):
            sh = nd_shift(stamp, (p[2], p[1]), order=3, mode='constant', cval=0.0)
            return HP.highpass_stamp(sh, sb)[c0 - FIT_BOX: c0 + FIT_BOX + 1, c0 - FIT_BOX: c0 + FIT_BOX + 1]

        def resid(p):
            return ((d - p[0] * model(p)) * w).ravel()
        peak = float(d.max())
        p0 = np.array([max(peak / max(float(HP.highpass_stamp(stamp, sb)[c0, c0]), 1e-9), 1.0), px - ix, py - iy])
        sol = least_squares(resid, p0, x_scale=np.array([max(p0[0], 1.0), 0.3, 0.3]))
        m = sol.x[0] * model(sol.x)
        out.append((float(np.sqrt(np.mean(sol.fun ** 2))), float(np.sqrt(np.mean((d - m) ** 2)) / m.max())))
    return np.array(out).reshape(-1, 2)


# ----------------------------------------------------------------------------------------------------------------------------------
# data of a patch
# ----------------------------------------------------------------------------------------------------------------------------------
_CACHE = {}


def inputs(name, pxy):
    """(patch dict(dn, psf, gain, dark), sky_dn array, truth dict of the frame) of a target patch, cached."""
    k = (name, tuple(pxy))
    if k not in _CACHE:
        patch, sky, truth, _ = D.load_inputs(name, pxy)
        _CACHE[k] = (patch, sky, truth)
    return _CACHE[k]


def frame_patches(name):
    return [p for n, p in TARGETS if n == name]


def stamp_arrays(prefix, est, ref_core, ci=None, nfail=0):
    """The arrays stored for one in-situ stamp: the stamp, the cores, the star table and the comparison with the psField core."""
    st = est['stars']
    out = {'stamp': est['stamp'], 'core': est['core'], 'core_obs': est['core_obs'], 'n_sel': np.array([len(st['c'])]), 'n_used': np.array([int(est['use'].sum())]), 'n_by_frame': st['n_by_frame'],
           'star_x': st['x'], 'star_y': st['y'], 'star_frame': st['frame'], 'star_snr': st['snr'], 'star_f': est['f'], 'star_e': est['e'], 'star_use': est['use'].astype(int), 'resid': np.array([est['resid']])}
    if est['ok']:
        out['vs_psf'] = compare_cores(stamp_core(est['stamp']), ref_core)
        out['shape'] = shape_of(est['stamp'])
    else:
        out['vs_psf'] = np.full(4, np.nan)
        out['shape'] = np.full(2, np.nan)
    if ci is not None:
        out['vs_psf_ci'] = ci
        out['boot_fail'] = np.array([nfail])
    return {prefix + '|' + k: np.asarray(v) for k, v in out.items()}


def plain_stack_stat(dn, x, y, idx, ref_core):
    """The diag-style plain stack (`exp_astro_diag.cutouts`, pixelwise median) of the stars idx at (x, y): [central pixel over the reference's minus 1, FWHM, number of stars]."""
    if len(idx) < 1:
        return np.full(3, np.nan)
    core = np.median(D.cutouts(dn, x, y, idx), 0)
    return np.array([core[CORE // 2, CORE // 2] / ref_core[CORE // 2, CORE // 2] - 1.0, D.fwhm(core), float(len(idx))])


def validate_patch(name, pxy, scopes_est, B=B_STAMP):
    """V2, V3, V4 for one patch (label-using where marked): arrays under the patch's own keys (prefix 'ST|patch|<name>|<xy>|') for the patch scope estimate, and 'ST|frame|...' for V2 of the frame scope."""
    patch, sky, truth = inputs(name, pxy)
    ref = stamp_core(patch['psf'])
    out = {}
    tx, ty, tm, tidx = truth_stars(truth, pxy[0], pxy[1], patch['dn'])
    c_t, off_t = star_cores(patch['dn'], tx, ty, tidx)
    est_t = estimate(c_t, off_t, patch['psf']) if len(tidx) >= N_MIN else None
    out['truth|n'] = np.array([len(tidx)])
    out['truth|plain'] = plain_stack_stat(patch['dn'], tx, ty, tidx, ref)
    if est_t is not None and est_t['ok']:
        out['truth|stamp'] = est_t['stamp']
        out['truth|vs_psf'] = compare_cores(stamp_core(est_t['stamp']), ref)
        out['truth|n_used'] = np.array([int(est_t['use'].sum())])
        core_t = stamp_core(est_t['stamp'])
    else:
        core_t = None
        out['truth|vs_psf'] = np.full(4, np.nan)
        out['truth|n_used'] = np.array([0])
    res = {}
    for scope, est in scopes_est.items():
        o = {}
        o['vs_truth'] = compare_cores(stamp_core(est['stamp']), core_t) if (core_t is not None and est['ok']) else np.full(4, np.nan)
        st = est['stars']
        if scope == 'patch':
            sx, sy = st['x'], st['y']
            o['plain'] = plain_stack_stat(patch['dn'], sx, sy, np.arange(len(sx)), ref)
            mag, ind = nearest_truth(sx, sy, *D.sources_near(truth, pxy[0], pxy[1], D.MARGIN)) if len(sx) else (np.zeros(0), np.zeros(0, int))
            tx_, ty_, tm_ = D.sources_near(truth, pxy[0], pxy[1], D.MARGIN)
            iso = truth_isolated_mask(ind, tx_, ty_, tm_)
            o['truthmatch'] = np.array([int(np.isfinite(mag).sum()), int(iso.sum()), float(np.nanmedian(mag)) if np.isfinite(mag).any() else np.nan])
            # V4: held-out truth stars
            if est['ok'] and len(tidx):
                far = np.array([np.min(np.hypot(sx - tx[i], sy - ty[i])) if len(sx) else np.inf for i in tidx])
                held = tidx[far > HELD_RAD]
                pos = [(tx[i], ty[i]) for i in held]
                r_psf = hp_fit_rms(patch['dn'], patch['gain'], patch['dark'], patch['psf'], pos)
                r_ins = hp_fit_rms(patch['dn'], patch['gain'], patch['dark'], est['stamp'], pos)
                o['held'] = np.hstack([r_psf, r_ins])                # columns: standardized rms (psField), fractional rms (psField), standardized rms (in-situ), fractional rms (in-situ)
            else:
                o['held'] = np.zeros((0, 4))
        res[scope] = o
    # V3: recovery in the exact-stamp world
    w = D.render_world(truth, pxy[0], pxy[1], patch['psf'], patch['gain'], patch['dark'], sky, D.world_seed(name, pxy))
    frames = [dict(dn=w['dn'], gain=patch['gain'], dark=patch['dark'])]
    ew = insitu_stamp(frames, patch['psf'])
    out['world|n_sel'] = np.array([len(ew['stars']['c'])])
    out['world|vs_planted'] = compare_cores(stamp_core(ew['stamp']), ref) if ew['ok'] else np.full(4, np.nan)
    sx, sy = ew['stars']['x'], ew['stars']['y']
    out['world|plain'] = plain_stack_stat(w['dn'], sx, sy, np.arange(len(sx)), ref)
    return out, res


def frame_truth_pooled(name, est_f):
    """V2 pooled by frame (label-using): the truth stars of all the target patches of the frame through the same estimator (the wings of the first patch's psField stamp, as `est_f`, the pooled in-situ
    estimate), against the in-situ frame stamp and against the mean of the patches' psField cores. Arrays under '<name>|all|'."""
    plist = frame_patches(name)
    cs, offs, plain = [], [], []
    for p in plist:
        patch, sky, truth = inputs(name, p)
        tx, ty, tm, tidx = truth_stars(truth, p[0], p[1], patch['dn'])
        c, off = star_cores(patch['dn'], tx, ty, tidx)
        cs.append(c)
        offs.append(off)
    c, off = np.concatenate(cs), np.concatenate(offs)
    ref = np.mean([stamp_core(inputs(name, p)[0]['psf']) for p in plist], 0)
    out = {'truth_n': np.array([len(c)]), 'psf_mean_core': ref}
    et = estimate(c, off, inputs(name, plist[0])[0]['psf']) if len(c) >= N_MIN else None
    nan4 = np.full(4, np.nan)
    if et is not None and et['ok'] and est_f['ok']:
        ct = stamp_core(et['stamp'])
        out['truth_n_used'] = np.array([int(et['use'].sum())])
        out['truth_vs_psfmean'] = compare_cores(ct, ref)
        out['insitu_vs_psfmean'] = compare_cores(stamp_core(est_f['stamp']), ref)
        out['insitu_vs_truth'] = compare_cores(stamp_core(est_f['stamp']), ct)
    else:
        out['truth_n_used'] = np.array([0])
        out['truth_vs_psfmean'] = out['insitu_vs_psfmean'] = out['insitu_vs_truth'] = nan4
    return {'ST|frame|%s|all|%s' % (name, k): v for k, v in out.items()}


# ----------------------------------------------------------------------------------------------------------------------------------
# raw
# ----------------------------------------------------------------------------------------------------------------------------------
def load_raw(path=None):
    path = path or RAW
    return CF.load_raw(path) if os.path.exists(path) else {}


def save_part(arrays, path=None):
    """Merge arrays into the raw (write-once, atomic, locked: `exp_astro_cyclefix.save_part`)."""
    CF.save_part(path or RAW, arrays)


def sha256_file(path):
    return CF.sha256_file(path)


CODE_FILES = ('src/core/imaging.py', 'src/core/imaging_guard.py', 'src/core/psf_general.py', 'experiments/astro_common.py', 'experiments/exp_astro_m2.py', 'experiments/exp_astro_psf.py',
              'experiments/exp_astro_diag.py', 'experiments/exp_astro_baselines_hp.py', 'experiments/exp_astro_baselines_hp2.py', 'experiments/exp_astro_cyclefix.py')


def meta_arrays():
    """The hashes of the code this script imports and of the freeze record, and the versions, recorded when the first result is saved."""
    import scipy
    files = [f for f in CODE_FILES + ('results/astro/freeze.json',) if os.path.exists(os.path.join(REPO, f))]
    return {'meta|files': np.array(files), 'meta|sha256': np.array([sha256_file(os.path.join(REPO, f)) for f in files]),
            'meta|versions': np.array(['python %d.%d.%d' % sys.version_info[:3], 'numpy ' + np.__version__, 'scipy ' + scipy.__version__])}


def save_with_meta(arrays, path=None):
    cur = load_raw(path)
    if 'meta|files' not in cur:
        arrays = dict(arrays, **meta_arrays())
    save_part(arrays, path)


def check_hashes(path=None):
    """Raises if a file hashed in the raw by the first run differs now."""
    cur = load_raw(path)
    if 'meta|files' not in cur:
        return
    now = meta_arrays()
    rec = dict(zip((str(f) for f in cur['meta|files']), (str(h) for h in cur['meta|sha256'])))
    cur_h = dict(zip((str(f) for f in now['meta|files']), (str(h) for h in now['meta|sha256'])))
    diff = sorted(f for f, h in rec.items() if cur_h.get(f) != h)
    if diff:
        raise RuntimeError('files differ from what the first run recorded: %s' % ', '.join(diff))


# ----------------------------------------------------------------------------------------------------------------------------------
# --stamps
# ----------------------------------------------------------------------------------------------------------------------------------
def part_stamps(path=None, B=B_STAMP):
    """V1 to V4 and the stamps of both scopes for the ten target patches, saved patch by patch (a patch whose keys exist is not recomputed)."""
    cur = load_raw(path)
    done = lambda scope, name, pxy: 'ST|%s|%s|n_sel' % (scope, pkey(name, pxy)) in cur           # noqa: E731
    est_frame = {}
    names = list(dict.fromkeys(n for n, _ in TARGETS))
    for name in names:
        plist = frame_patches(name)
        frames = [dict(dn=inputs(name, p)[0]['dn'], gain=inputs(name, p)[0]['gain'], dark=inputs(name, p)[0]['dark']) for p in plist]
        est_frame[name] = insitu_stamp(frames, inputs(name, plist[0])[0]['psf'])         # the wings of the frame stamp: the first patch's psField stamp is replaced per patch below
    for name in names:
        if 'ST|frame|%s|all|truth_n' % name not in cur:
            save_with_meta(frame_truth_pooled(name, est_frame[name]), path)
    for name, pxy in TARGETS:
        if done('patch', name, pxy) and done('frame', name, pxy):
            print('ST %s %s: exists' % (name, pxy), flush=True)
            continue
        t0 = time.time()
        patch, sky, truth = inputs(name, pxy)
        ref = stamp_core(patch['psf'])
        fr = dict(dn=patch['dn'], gain=patch['gain'], dark=patch['dark'])
        est_p = insitu_stamp([fr], patch['psf'])
        # the frame scope: the pooled cut-outs, the wings of THIS patch's psField stamp
        plist = frame_patches(name)
        st = est_frame[name]['stars']
        est_f = estimate(st['c'], st['off'], patch['psf']) if len(st['c']) >= N_MIN else est_frame[name]
        est_f['stars'] = st
        arrays = {}
        ci_p, nf_p = star_bootstrap(est_p['stars'], [ref], [patch['psf']], est_p['stars']['n_by_frame'], B, 0) if est_p['ok'] else (np.full((1, 3, 2), np.nan), B)
        arrays.update(stamp_arrays('ST|patch|' + pkey(name, pxy), est_p, ref, ci_p, nf_p))
        sizes_f = st['n_by_frame']
        # the frame bootstrap against every patch of the frame at once: the stamp's wings differ per patch, so a resample is re-estimated with this patch's psField stamp
        ci_f, nf_f = star_bootstrap(st, [ref], [patch['psf']], sizes_f, B, 0) if est_f['ok'] else (np.full((1, 3, 2), np.nan), B)
        arrays.update(stamp_arrays('ST|frame|' + pkey(name, pxy), est_f, ref, ci_f, nf_f))
        out, res = validate_patch(name, pxy, {'patch': est_p, 'frame': est_f})
        for k, v in out.items():
            arrays['ST|patch|%s|%s' % (pkey(name, pxy), k)] = np.asarray(v)
        for scope, o in res.items():
            for k, v in o.items():
                arrays['ST|%s|%s|%s' % (scope, pkey(name, pxy), k)] = np.asarray(v)
        arrays['ST|psf|%s|stamp' % pkey(name, pxy)] = patch['psf']
        arrays['ST|psf|%s|shape' % pkey(name, pxy)] = shape_of(patch['psf'])
        arrays['ST|psf|%s|core' % pkey(name, pxy)] = ref
        save_with_meta(arrays, path)
        print('ST %s %s: patch scope %d stars (%d used), centre %+.3f fwhm %.2f (psField %.2f); frame scope %d stars (%d used), centre %+.3f; %.0f s' % (
            name, pxy, arrays['ST|patch|%s|n_sel' % pkey(name, pxy)][0], arrays['ST|patch|%s|n_used' % pkey(name, pxy)][0], arrays['ST|patch|%s|vs_psf' % pkey(name, pxy)][0],
            arrays['ST|patch|%s|vs_psf' % pkey(name, pxy)][1], arrays['ST|patch|%s|vs_psf' % pkey(name, pxy)][2], arrays['ST|frame|%s|n_sel' % pkey(name, pxy)][0],
            arrays['ST|frame|%s|n_used' % pkey(name, pxy)][0], arrays['ST|frame|%s|vs_psf' % pkey(name, pxy)][0], time.time() - t0), flush=True)


def scrambled_stamp(psf, insitu, seed):
    """The control stamp of the sensitivity analysis: the psField stamp plus the pixels of (in-situ - psField) in a random order (a permutation of the 51 x 51 pixels, `default_rng(seed)`), clipped at 0, unit sum. It has
    the amplitude distribution of the in-situ perturbation and none of its structure (no change of the shape, the width, the centring)."""
    d = (insitu - psf).ravel()
    s = np.clip(psf + np.random.default_rng(seed).permutation(d).reshape(psf.shape), 0.0, None)
    return s / s.sum()


def scr_seed(name, pxy):
    return [SCR_SEED, zlib.crc32(name.encode()), int(pxy[0]), int(pxy[1])]


def part_scr(path=None):
    """The control stamps (scope 'scr') of the target patches that have a patch-scope in-situ stamp, stored under 'ST|scr|<patch>|' (stamp, the comparison with psField, the rms distance from the in-situ stamp)."""
    cur = load_raw(path)
    arrays = {}
    for name, pxy in TARGETS:
        k = 'ST|patch|%s|stamp' % pkey(name, pxy)
        if k not in cur or not np.all(np.isfinite(cur[k])):
            continue
        patch = inputs(name, pxy)[0]
        st = scrambled_stamp(patch['psf'], cur[k], scr_seed(name, pxy))
        p = 'ST|scr|%s|' % pkey(name, pxy)
        arrays[p + 'stamp'] = st
        arrays[p + 'vs_psf'] = compare_cores(stamp_core(st), stamp_core(patch['psf']))
        arrays[p + 'rms_to_psf'] = np.array([np.sqrt(np.mean((st - patch['psf']) ** 2)) / patch['psf'].max(), np.sqrt(np.mean((cur[k] - patch['psf']) ** 2)) / patch['psf'].max()])
    save_with_meta(arrays, path)
    print('SCR: %d control stamps' % (len(arrays) // 3))


# ----------------------------------------------------------------------------------------------------------------------------------
# --ref: the psField-stamp rows
# ----------------------------------------------------------------------------------------------------------------------------------
def ref_sources(name, pxy):
    """{method: (raw name, key prefix of the scores, key prefix of the run)} of the psField-stamp rows of a target patch in the frozen raws."""
    x0, y0 = pxy
    if name == 'm2':
        iso = 'C|T|x%d_y%d|ours' % pxy
        adq = iso
        cd = 'T|x%d_y%d|crowdsource (default)' % pxy
    else:
        iso = 'C|T2|%s|x%d_y%d|ours' % (name, x0, y0)
        adq = 'C|P|%s|x%d_y%d|ours' % (name, x0, y0)
        cd = 'T2|%s|x%d_y%d|crowdsource (default)' % (name, x0, y0)
    return dict(iso=('cyc', iso), adq=('cyc', adq), cs_hp=('hp2', 'T|%s|x%d_y%d|cs2' % (name, x0, y0)), phu_hp=('hp2', 'T|%s|x%d_y%d|phu2' % (name, x0, y0)), cs_def=('m2', cd))


def part_ref(path=None):
    """Copies the psField-stamp rows (tot, blk22, detections, seconds, entries) of the frozen results into the raw under 'REF|<set>|x..|<method>|' (read only; write-once)."""
    raws = dict(cyc=CF.load_raw(CYCLEFIX_RAW), hp2=CF.load_raw(HP2_RAW), m2=CF.load_raw(M2_RAW))
    arrays = {}
    for name, pxy in TARGETS:
        for meth, (rn, pre) in ref_sources(name, pxy).items():
            z = raws[rn]
            base = 'REF|%s|%s' % (pkey(name, pxy), meth)
            if meth in ('iso', 'adq'):
                for k in ('k2', 'k0'):
                    arrays['%s|%s|tot' % (base, k)] = z['%s|%s|tot' % (pre, k)]
                    arrays['%s|%s|blk22' % (base, k)] = z['%s|%s|blk22' % (pre, k)]
                arrays[base + '|secs'] = z[pre + '|secs']
                arrays[base + '|det'] = z[pre + '|det']
                arrays[base + '|barred'] = np.array([len(z[pre + '|trace|barred'])])
                arrays[base + '|bic'] = z[pre + '|bic']
            else:
                arrays[base + '|k2|tot'] = z[pre + '|tot']
                arrays[base + '|k2|blk22'] = z[pre + '|blk22']
                if pre + '|secs' in z:
                    arrays[base + '|secs'] = z[pre + '|secs']
    arrays['REF|sources|sha256'] = np.array([sha256_file(p) for p in (CYCLEFIX_RAW, HP2_RAW, M2_RAW)])
    arrays['REF|sources|names'] = np.array(['cyclefix_raw.npz', 'baselines_hp2_raw.npz', 'm2_raw.npz'])
    save_with_meta(arrays, path)
    print('REF: %d arrays' % len(arrays))


# ----------------------------------------------------------------------------------------------------------------------------------
# --run
# ----------------------------------------------------------------------------------------------------------------------------------
def frozen_config():
    """(hp, rex, kap) of the frozen configuration, after the freeze check (`exp_astro_m2.check_freeze`)."""
    raw = X.load_raw(M2_RAW)
    X.check_freeze(raw)
    hp, rex, kap, _ = X.choose(raw)
    return hp, rex, kap


def with_guard(fn, *args, **kw):
    """Runs fn(*args, **kw) with `imaging.sliding_path` replaced by the guarded function and its side channel (`exp_astro_cyclefix.install`); returns (result, trace of the single path run)."""
    CF.install()
    CF.TRACE.clear()
    try:
        res = fn(*args, **kw)
        if len(CF.TRACE) != 1:
            raise RuntimeError('%d paths were run, expected one' % len(CF.TRACE))
        return res, CF.TRACE[0]
    finally:
        CF.uninstall()


def pack_iss(dets, secs, n, info, t, kap, trace):
    """The arrays of an ISS run, keys relative to the run: scores of the readouts kappa and BIC ('k2|...', 'k0|...'), the catalogue at kappa, seconds and entries, the BIC record, the guard's trace."""
    arrays = {}
    for k in (kap, 0.0):
        arrays.update(X._pack('k%g' % k, X.score_dets(dets[k], t)))
    d = dets[kap]
    arrays['det'] = np.c_[d[0], d[1], d[2]]
    arrays['secs'] = np.array([secs, n])
    arrays['bic'] = np.array([info['bic_index'], info['n_blocks'], int(n >= X.MAX_ENTRIES)])
    arrays['trace|counts'] = trace['counts']
    arrays['trace|barred'] = trace['barred']
    return arrays


def run_iss(kind, patch, truth_patch, cfg, fit=None):
    """The ISS on `patch` (its 'psf' is the stamp) at the frozen configuration, guarded: kind 'iso' (isotropic kernel, `exp_astro_m2.run_path`) or 'adq' (the kernel `fit` of `psf_general.fit_adequate`,
    `exp_astro_psf.run_path_psf`). Returns the arrays relative to the run key."""
    hp, rex, kap = cfg
    if kind == 'iso':
        (dets, secs, n, info), tr = with_guard(X.run_path, patch, hp, rex)
    else:
        (dets, secs, n, info), tr = with_guard(PS.run_path_psf, patch, hp, rex, fit)
    return pack_iss(dets, secs, n, info, truth_patch, kap, tr)


def round2_specs():
    """(crowdsource spec, photutils spec) of the round-2 selection (verified)."""
    sel = HS.verify_selection2(CF.load_raw(HP2_RAW))
    return next(d for d in H2.CS2 if d['tag'] == sel['crowdsource']), next(d for d in H2.PHOT2 if d['tag'] == sel['photutils'])


def run_baseline(meth, patch, sky_mean, truth_patch, specs):
    """A baseline on `patch` ('psf' = the stamp): arrays relative to the run key ('k2|' scores, 'det', 'secs')."""
    if meth == 'cs_hp':
        cat, secs = H2.run_crowdsource_hp2(patch['dn'], patch['psf'], patch['gain'], patch['dark'], specs[0])
        det = (cat['x'], cat['y'], cat['flux'])
    elif meth == 'phu_hp':
        cat, secs = H2.run_photutils_hp2(patch['dn'], patch['psf'], patch['gain'], patch['dark'], specs[1])
        det = (cat['x'], cat['y'], cat['flux'])
    elif meth == 'cs_def':
        with contextlib.redirect_stdout(io.StringIO()):
            cat, secs = HP.AB.run_crowdsource(patch['dn'], sky_mean, patch['psf'], patch['gain'], patch['dark'])
        key = 'x0_y0|crowdsource|'
        b = {key + k: cat[k] for k in ('x', 'y', 'flux', 'dflux', 'fracflux', 'qf', 'spread')}
        det = X.cs_dets(b, (0, 0), (0.0, 0.0, 0.0, np.inf))
    else:
        raise ValueError(meth)
    arrays = X._pack('k2', X.score_dets(det, truth_patch))
    arrays['det'] = np.c_[det[0], det[1], det[2]]
    arrays['secs'] = np.array([secs])
    return {k: np.asarray(v) for k, v in arrays.items()}


def part_run(scope, methods, sets, path=None):
    """The runs of `methods` on the target patches of `sets` with the in-situ stamp of `scope` (from the raw; `--stamps` first); a run whose keys exist is not run again; results are saved as they arrive."""
    check_hashes(path)
    cfg = frozen_config()
    specs = round2_specs()
    for name, pxy in TARGETS:
        if name not in sets:
            continue
        cur = load_raw(path)
        k0 = 'ST|%s|%s|' % (scope, pkey(name, pxy))
        if k0 + 'stamp' not in cur or not np.all(np.isfinite(cur[k0 + 'stamp'])):
            print('RUN %s %s %s: no in-situ stamp, skipped' % (scope, name, pxy), flush=True)
            continue
        patch0, sky, truth = inputs(name, pxy)
        patch = dict(patch0, psf=cur[k0 + 'stamp'])
        t = A.truth_for_patch(truth, pxy[0], pxy[1])
        fit = PG.fit_adequate(patch['psf'])
        for meth in methods:
            pre = 'RUN|%s|%s|%s' % (scope, pkey(name, pxy), meth)
            if pre + '|wall' in cur or pre + '|same' in cur:
                print('RUN %s %s %s %s: exists' % (scope, name, pxy, meth), flush=True)
                continue
            t0 = time.time()
            if meth == 'adq' and fit['kind'] == 'isotropic':
                save_with_meta({pre + '|same': np.array([1]), pre + '|kernel': np.array([fit['kind']]), pre + '|rms': np.array([fit['rms_isotropic'], fit['rms']])}, path)
                print('RUN %s %s %s adq: the rule keeps the isotropic kernel (same as iso)' % (scope, name, pxy), flush=True)
                continue
            if meth in ('iso', 'adq'):
                arr = run_iss(meth, patch, t, cfg, fit if meth == 'adq' else None)
                arr['kernel'] = np.array([fit['kind'] if meth == 'adq' else 'isotropic'])
                arr['rms'] = np.array([fit['rms_isotropic'], fit['rms']])
                arr['same'] = np.array([0])
            else:
                arr = run_baseline(meth, patch, float(sky.mean()), t, specs)
            out = {pre + '|' + k: v for k, v in arr.items()}
            out[pre + '|wall'] = np.array([time.time() - t0])
            save_with_meta(out, path)
            print('RUN %s %s %s %s: F %.3f, %.0f s' % (scope, name, pxy, meth, X.f_measure(out[pre + '|k2|tot'])[0], time.time() - t0), flush=True)


# ----------------------------------------------------------------------------------------------------------------------------------
# Task B: crowdsource options on the development patches
# ----------------------------------------------------------------------------------------------------------------------------------
CS_BASE = dict(sb=2.0, thr=3.0, miniter=8, maxiter=20, fewstars=20, psfderiv=True, blendthreshu=2, titer_thresh=2, psfvalsharpcutfac=0.7, psfsharpsat=0.7, weights='iss', blist=0)
CS_BASE_TAG = 'cs2_hp2_t3_pm'                # the round-2 selection that the base must reproduce
# stage 1: option family -> levels (declared); the tag is 'cs3_<family><level>'
CS_LEVELS = (('blendthreshu', (1, 3, 4)), ('titer_thresh', (1, 3, 5)), ('psfderiv', (False,)), ('maxiter', (30, 40)), ('miniter', (4, 12)), ('fewstars', (5,)), ('psfvalsharpcutfac', (0.5, 0.9)),
             ('psfsharpsat', (0.5, 1.0)), ('thr', (2.5, 4.0)), ('blist', (5, 15)), ('weights', ('flat', 'smooth3', 'satmask')))
COMBO_GAIN = 0.003                           # a stage-1 option enters the combination if its pooled development F exceeds the base's by at least this
FINDING_GAIN = 0.010                         # a variant is a finding if it beats the round-2 selection by more than this pooled development F


def cs3_tag(fam, level):
    return 'cs3_%s%s' % (fam, str(level).replace('.', 'p'))


def cs3_variants():
    """[(tag, family, level, spec)] of stage 1 (the base first): spec = CS_BASE with one option changed."""
    out = [('cs3_base', None, None, dict(CS_BASE))]
    for fam, levels in CS_LEVELS:
        for lv in levels:
            out.append((cs3_tag(fam, lv), fam, lv, dict(CS_BASE, **{fam: lv})))
    return out


def cs3_weights(dn, gain, dark, kind):
    """The inverse-noise weights of a Task B variant: 'iss' (`exp_astro_baselines_hp.iss_weights`), 'flat' (the constant median of them), 'smooth3' (the ISS formula with a 3 px smoothing of the expected counts),
    'satmask' (the ISS weights, 0 where dn > SAT_DN)."""
    sw = HP.iss_weights(dn, gain, dark)
    if kind == 'iss':
        return sw
    if kind == 'flat':
        return np.full(dn.shape, float(np.median(sw)))
    if kind == 'smooth3':
        expected = np.maximum(gaussian_filter(dn, 3.0, mode='reflect'), np.percentile(dn, 1))
        return 1.0 / np.sqrt(expected / gain + dark)
    if kind == 'satmask':
        return np.where(dn > SAT_DN, 0.0, sw)
    raise ValueError(kind)


def cs3_blist(dn, gain, dark, n):
    """crowdsource's bright-star list (x along the first array axis, y along the second, magnitude) of the n brightest maxima of the step-1 detection (the first-axis coordinate is our row)."""
    x, y, pk, noise = detect_maxima(dn, gain, dark)
    o = np.argsort(-pk, kind='stable')[:n]
    return (y[o].astype('f4'), x[o].astype('f4'), np.zeros(len(o), 'f4'))


def run_cs3(dn, psf51, gain, dark, spec):
    """One Task B crowdsource variant (needs the domain environment): the round-1/2 high-passed call with the options of `spec`; returns the catalogue dict (x, y swapped as there) and the CPU seconds."""
    import crowdsource.psf as psfmod
    from crowdsource import crowdsource_base
    t0 = time.process_time()
    sw = cs3_weights(dn, gain, dark, spec['weights'])
    base = HP.cs_stamp(psf51)
    s19 = float(psfmod.central_stamp(base, 19).sum())
    psf = psfmod.SimplePSF(HP.highpass_stamp(base / s19, spec['sb']), normalize=0)
    kw = dict(maxstars=10000000, fewstars=spec['fewstars'], miniter=spec['miniter'], maxiter=spec['maxiter'], ntilex=1, ntiley=1, refit_psf=False, verbose=False, threshold=spec['thr'],
              psfderiv=spec['psfderiv'], blendthreshu=spec['blendthreshu'], titer_thresh=spec['titer_thresh'], psfvalsharpcutfac=spec['psfvalsharpcutfac'], psfsharpsat=spec['psfsharpsat'])
    if spec['blist']:
        kw['blist'] = cs3_blist(dn, gain, dark, spec['blist'])
    flag = np.zeros(dn.shape, dtype=np.int64)
    with contextlib.redirect_stdout(io.StringIO()):
        res = crowdsource_base.fit_im(HP.highpass_frame(dn, spec['sb']), psf, weights=np.asarray(sw, dtype=np.float64), dq=flag, **kw)
    cat = res['stars']
    return dict(x=np.asarray(cat['y'], float), y=np.asarray(cat['x'], float), flux=np.asarray(cat['flux'], float)), time.process_time() - t0


def cs_dev_patches():
    return [tuple(p) for p in X.DEV]


def part_cs_variant(tag, spec, patches, truth, prefix, path=None):
    """Runs one variant on the patches (a patch whose result exists is not run again) and saves the scores under '<prefix>|<tag>|x<x0>_y<y0>|'; a failure is recorded as 'failed'."""
    cur = load_raw(path)
    for p in patches:
        k = '%s|%s|x%d_y%d|' % (prefix, tag, p[0], p[1])
        if k + 'tot' in cur or k + 'failed' in cur:
            continue
        a = HP.load_arrays('m2', p[0], p[1])
        t = A.truth_for_patch(truth, p[0], p[1])
        try:
            cat, secs = run_cs3(a['dn'], a['psf51'], a['gain'], a['dark'], spec)
        except Exception as exc:                                                     # noqa: BLE001  (recorded, not retried)
            save_with_meta({k + 'failed': np.array([('%s: %s' % (type(exc).__name__, exc))[:300]])}, path)
            print('CSD %s %s: FAILED %s' % (tag, p, exc), flush=True)
            continue
        arrays = {k + f.lstrip('|'): v for f, v in X._pack('', X.score_dets((cat['x'], cat['y'], cat['flux']), t)).items()}
        arrays[k + 'secs'] = np.array([secs])
        arrays[k + 'ndet'] = np.array([len(cat['x'])])
        save_with_meta(arrays, path)
        print('CSD %s %s: F %.3f (%d detections, %.1f s)' % (tag, p, X.f_measure(arrays[k + 'tot'])[0], len(cat['x']), secs), flush=True)


def dev_f(raw, tag, prefix='CSD'):
    """Pooled F at mag < 22 of a variant over the two development patches (NaN if a patch is missing or failed)."""
    ks = ['%s|%s|x%d_y%d|tot' % (prefix, tag, p[0], p[1]) for p in cs_dev_patches()]
    if not all(k in raw for k in ks):
        return float('nan')
    return float(X.f_measure(sum(raw[k] for k in ks))[0])


def combo_spec(raw):
    """(spec, options): the stage-2 combination declared in the docstring: for every family the level with the largest pooled development F among those that beat the base by COMBO_GAIN or more."""
    base = dev_f(raw, 'cs3_base')
    spec, picked = dict(CS_BASE), []
    for fam, levels in CS_LEVELS:
        cand = [(dev_f(raw, cs3_tag(fam, lv)), lv) for lv in levels]
        cand = [(f, lv) for f, lv in cand if np.isfinite(f) and f >= base + COMBO_GAIN - 1e-12]
        if cand:
            f, lv = max(cand, key=lambda fl: (round(fl[0], 6), -levels.index(fl[1])))
            spec[fam] = lv
            picked.append('%s=%s' % (fam, lv))
    return spec, picked


def part_cs_check(path=None, hp2_raw=None):
    """Stores whether the base reproduces the stored development totals of the round-2 selection (`baselines_hp2_raw.npz`, read only): 'CSD|check|<patch>' = 1 if the totals are equal."""
    hp2_raw = hp2_raw or HP2_RAW
    if not os.path.exists(hp2_raw):
        return
    raw, z = load_raw(path), CF.load_raw(hp2_raw)
    out = {}
    for p in cs_dev_patches():
        a, b = 'CSD|cs3_base|x%d_y%d|tot' % p, 'DEV|x%d_y%d|%s|tot' % (p[0], p[1], CS_BASE_TAG)
        if a in raw and b in z:
            out['CSD|check|x%d_y%d' % p] = np.array([int(np.array_equal(raw[a], z[b]))])
    if out:
        save_with_meta(out, path)


def part_cs_dev(stage, path=None):
    """Task B on the development patches: stage 1 (the base and every one-option variant), stage 2 (the combination; then the selection, written once)."""
    check_hashes(path)
    truth = X.truth_arrays()
    if stage == 1:
        for tag, fam, lv, spec in cs3_variants():
            part_cs_variant(tag, spec, cs_dev_patches(), truth, 'CSD', path)
        part_cs_check(path)
        return
    raw = load_raw(path)
    if 'CSD|SEL|best' in raw:
        raise RuntimeError('a selection exists: it is never made twice')
    spec, picked = combo_spec(raw)
    if picked:
        part_cs_variant('cs3_combo', spec, cs_dev_patches(), truth, 'CSD', path)
        raw = load_raw(path)
    tags = [t for t, _, _, _ in cs3_variants()] + (['cs3_combo'] if picked else [])
    scores = {t: dev_f(raw, t) for t in tags}
    ok = [t for t in tags if np.isfinite(scores[t])]
    best = max(ok, key=lambda t: (round(scores[t], 6), t == 'cs3_base'))
    base_f = scores['cs3_base']
    finding = round((round(scores[best], 6) - round(base_f, 6)) * 1e6) > round(FINDING_GAIN * 1e6)            # integer micro-units: 0.450 - 0.440 is not more than 0.010
    save_with_meta({'CSD|SEL|best': np.array([best]), 'CSD|SEL|best_F': np.array([scores[best]]), 'CSD|SEL|base_F': np.array([base_f]), 'CSD|SEL|combo': np.array(picked or ['none']),
                    'CSD|SEL|finding': np.array([int(finding)])}, path)
    print('CSD selection: best %s F %.4f against the base %.4f; finding %s' % (best, scores[best], base_f, finding))


def part_cs_test(path=None):
    """Only if the selection is a finding: the best variant and the base (the round-2 selection) on the four M2 test patches, once."""
    check_hashes(path)
    raw = load_raw(path)
    if 'CSD|SEL|finding' not in raw or int(raw['CSD|SEL|finding'][0]) != 1:
        raise RuntimeError('the best development variant is not a finding (more than %.3f above the round-2 selection): no test scoring' % FINDING_GAIN)
    best = str(raw['CSD|SEL|best'][0])
    spec_by_tag = {t: s for t, _, _, s in cs3_variants()}
    spec_by_tag['cs3_combo'] = combo_spec(raw)[0]
    truth = X.truth_arrays()
    for tag in ('cs3_base', best):
        part_cs_variant(tag, spec_by_tag[tag], [tuple(p) for p in X.TESTS], truth, 'CSTEST', path)


# ----------------------------------------------------------------------------------------------------------------------------------
# report: accessors and statistics (from the raw alone)
# ----------------------------------------------------------------------------------------------------------------------------------
def nan_ok(v):
    return v is not None and np.all(np.isfinite(v))


def pct(v, nd=1):
    return '-' if not np.isfinite(v) else '%+.*f' % (nd, 100.0 * v)


def pct_ci(v, ci, nd=1):
    if not np.isfinite(v):
        return '-'
    return '%+.*f [%+.*f, %+.*f]' % (nd, 100.0 * v, nd, 100.0 * ci[0], nd, 100.0 * ci[1])


def num(v, nd=2):
    return '-' if not np.isfinite(v) else '%.*f' % (nd, v)


def run_prefix(raw, stamp, name, pxy, meth):
    """The key prefix of the results of a method on a patch with a stamp ('psf': the frozen rows copied by --ref; else the run). Where the rule kept the isotropic kernel, adq is the iso run."""
    if stamp == 'psf':
        return 'REF|%s|%s' % (pkey(name, pxy), meth)
    pre = 'RUN|%s|%s|%s' % (stamp, pkey(name, pxy), meth)
    if meth == 'adq' and pre + '|same' in raw and int(raw[pre + '|same'][0]) == 1:
        pre = 'RUN|%s|%s|iso' % (stamp, pkey(name, pxy))
    return pre


def blk_of(raw, stamp, name, pxy, meth, ro='k2'):
    k = '%s|%s|blk22' % (run_prefix(raw, stamp, name, pxy, meth), ro)
    return raw[k] if k in raw else None


def is_stalled(raw, stamp, name, pxy, meth):
    """True where the ISS path of the row reached the entry cap (a defective run)."""
    if meth not in ('iso', 'adq'):
        return False
    k = run_prefix(raw, stamp, name, pxy, meth) + '|bic'
    return k in raw and len(raw[k]) >= 3 and int(raw[k][2]) == 1


def pool(raw, name, combos):
    """The patches of the set `name` that have a valid (present, not stalled) row for every (stamp, method) of `combos`, their blocks pooled: dict(patches, blk {combo: (n, 4)}, sizes, left [(patch, [combo])])."""
    used, left = [], []
    for n, p in TARGETS:
        if n != name:
            continue
        bad = [c for c in combos if blk_of(raw, c[0], n, p, c[1]) is None or is_stalled(raw, c[0], n, p, c[1])]
        (left if bad else used).append((p, bad) if bad else p)
    blk = {c: np.concatenate([blk_of(raw, c[0], name, p, c[1]) for p in used]) if used else np.zeros((0, 4)) for c in combos}
    return dict(patches=used, blk=blk, sizes=[25] * len(used), left=left)


def boot_dd(a1, b1, a0, b0, idx):
    """((a1 - b1) - (a0 - b0)) in F with its 95 % interval over the resamples idx: the change of the gap a - b from the rows 0 to the rows 1."""
    f = lambda blk: X.cpf(blk.sum(0))[2]                      # noqa: E731
    fs = lambda blk: X.cpf(blk[idx].sum(1))[2]                # noqa: E731
    d0 = (f(a1) - f(b1)) - (f(a0) - f(b0))
    dd = (fs(a1) - fs(b1)) - (fs(a0) - fs(b0))
    return d0, np.percentile(dd, [2.5, 97.5])


def diff_cell(d, ci):
    return '%s %s' % (X.fmt_iv(d, ci), X.verdict(d, ci))


PAIRS = (('iso', 'cs_hp'), ('adq', 'cs_hp'), ('iso', 'phu_hp'), ('adq', 'phu_hp'), ('iso', 'cs_def'), ('adq', 'cs_def'))
PAIR_NAME = {('iso', 'cs_hp'): 'ISS iso minus crowdsource+HP', ('adq', 'cs_hp'): 'ISS adq minus crowdsource+HP', ('iso', 'phu_hp'): 'ISS iso minus photutils+HP',
             ('adq', 'phu_hp'): 'ISS adq minus photutils+HP', ('iso', 'cs_def'): 'ISS iso minus crowdsource default', ('adq', 'cs_def'): 'ISS adq minus crowdsource default'}


def gap_sentence(v0, v1, vd):
    """The plain statement of what happened to the gap between the psField row (verdict v0) and the in-situ row (v1); vd the verdict of the change of the gap."""
    base = {('loss', 'loss'): 'the ISS still loses', ('loss', 'tie'): 'the loss closes (the interval now contains 0)', ('loss', 'win'): 'the loss REVERSES (the ISS wins)',
            ('tie', 'loss'): 'a tie becomes a loss', ('tie', 'tie'): 'still a tie', ('tie', 'win'): 'a tie becomes a win', ('win', 'loss'): 'a win becomes a LOSS', ('win', 'tie'): 'the win closes to a tie',
            ('win', 'win'): 'the ISS still wins'}[(v0, v1)]
    ch = {'win': 'the gap changes in the ISS\'s favour (interval of the change excludes 0)', 'loss': 'the gap changes AGAINST the ISS (interval of the change excludes 0)', 'tie': 'the change of the gap is not distinguishable from 0'}[vd]
    return '%s; %s' % (base, ch)


def patch_label(p):
    return '(%d, %d)' % p


# ----------------------------------------------------------------------------------------------------------------------------------
# report sections
# ----------------------------------------------------------------------------------------------------------------------------------
def stamps_section(raw):
    L = ['## 1. The in-situ stamp: against the psField stamp, against stars selected with the truth, in the exact-stamp world, and on held-out stars', '']
    if not any(k.startswith('ST|patch|') for k in raw):
        return L + ['(no stamps in the raw yet)', '']
    for scope in ('patch', 'frame'):
        L += ['### 1.%s Scope: %s' % ('a' if scope == 'patch' else 'b', 'PATCH (primary: the stars of the patch itself)' if scope == 'patch' else 'FRAME (secondary: the stars of all the target patches of the frame pooled)'), '']
        rows = []
        for name, pxy in TARGETS:
            k = 'ST|%s|%s|' % (scope, pkey(name, pxy))
            if k + 'n_sel' not in raw:
                continue
            vs, ci = raw[k + 'vs_psf'], raw.get(k + 'vs_psf_ci', np.full((1, 3, 2), np.nan))
            sh = raw[k + 'shape']
            sh0 = raw['ST|psf|%s|shape' % pkey(name, pxy)]
            rows.append([name, patch_label(pxy), '%d / %d' % (raw[k + 'n_sel'][0], raw[k + 'n_used'][0]), pct_ci(vs[0], ci[0, 0]), '%s [%s, %s]' % (num(vs[1]), num(ci[0, 1, 0]), num(ci[0, 1, 1])), num(vs[2]),
                         '%s [%s, %s]' % (num(100 * vs[3], 1), num(100 * ci[0, 2, 0], 1), num(100 * ci[0, 2, 1], 1)), '%s / %s' % (num(sh[0], 3), num(sh0[0], 3)), num(raw[k + 'resid'][0] * 100, 2)])
        L += [X.table(['set', 'patch', 'stars selected / kept', 'central pixel, in-situ against psField (%)', 'FWHM in-situ [95 %]', 'FWHM psField', 'rms difference of the cores / psField peak (%)',
                       'ellipticity in-situ / psField', 'last fixed-point residual / peak (%)'], rows), '']
        for _, name in GROUPS:
            vals = [raw['ST|%s|%s|vs_psf' % (scope, pkey(name, p))] for n, p in TARGETS if n == name and 'ST|%s|%s|vs_psf' % (scope, pkey(name, p)) in raw]
            vals = [v for v in vals if np.isfinite(v[0])]
            if vals:
                m = np.mean(vals, 0)
                L += ['Mean over the %d %s patches with a stamp: central pixel %s %%, FWHM %s against %s.' % (len(vals), name, pct(m[0]), num(m[1]), num(m[2])), '']
    L += ['### 1.c Against the stars selected with the truth, patch scope (label-using)', '',
          'The truth stars (diag item 1: HUGS F606W < 20, more than 12 px from the edge, no source within 4.5 px brighter than m + 1.5, no pixel above %d DN) go through the SAME estimator (steps 3 to 6): the truth-star stamp. '
          '"Plain" is the pixelwise median stack of `exp_astro_diag` (centred, 6 px-normalised), whose alignment and crowding bias is the one `diag.md` reports (against the psField core).' % int(SAT_DN), '']
    rows = []
    for name, pxy in TARGETS:
        k = 'ST|patch|%s|' % pkey(name, pxy)
        if k + 'truth|n' not in raw:
            continue
        vt, vi, pl, pli = raw[k + 'truth|vs_psf'], raw[k + 'vs_truth'], raw[k + 'truth|plain'], raw[k + 'plain']
        tm = raw[k + 'truthmatch']
        rows.append([name, patch_label(pxy), '%d / %d' % (raw[k + 'truth|n'][0], raw[k + 'truth|n_used'][0]), pct(vt[0]), num(vt[1]), pct(vi[0]), num(vi[3] * 100), pct(pl[0]), pct(pli[0]),
                     '%d of %d (median F606W %s)' % (tm[1], tm[0], num(tm[2]))])
    L += [X.table(['set', 'patch', 'truth stars / kept', 'truth-star stamp against psField, central pixel (%)', 'FWHM of the truth-star stamp', 'in-situ against truth-star stamp, central pixel (%)',
                   'rms difference of the cores (% of peak)', 'plain stack, truth stars against psField (%)', 'plain stack, in-situ stars against psField (%)', 'in-situ stars that are truth-isolated, of those matched to a truth source (median F606W)'], rows), '']
    rows = []
    for _, name in GROUPS:
        k = 'ST|frame|%s|all|' % name
        if k + 'truth_n' not in raw:
            continue
        vt, vi, vv = raw[k + 'truth_vs_psfmean'], raw[k + 'insitu_vs_psfmean'], raw[k + 'insitu_vs_truth']
        rows.append([name, '%d / %d' % (raw[k + 'truth_n'][0], raw[k + 'truth_n_used'][0]), pct(vt[0]), num(vt[1]), pct(vi[0]), num(vi[1]), pct(vv[0]), num(vv[3] * 100), num(vt[3] * 100), num(vi[3] * 100)])
    L += ['Pooled by frame (the truth stars of all the target patches of the frame through the same estimator, the wings of the first patch\'s psField stamp; the reference is the mean of the patches\' psField cores):', '',
          X.table(['set', 'truth stars / kept', 'truth-star stamp against psField, central pixel (%)', 'FWHM truth-star stamp', 'in-situ (frame scope) against psField, central pixel (%)', 'FWHM in-situ (frame scope)',
                   'in-situ against truth-star stamp, central pixel (%)', 'rms difference in-situ - truth-star (% of peak)', 'rms difference truth-star - psField (% of peak)', 'rms difference in-situ - psField (% of peak)'], rows), '']
    L += ['### 1.d Recovery in the exact-stamp world and the held-out stars, patch scope', '',
          'World: the in-situ builder run on the exact-stamp synthetic frame of the patch (`exp_astro_diag.render_world`, the diag seed; the wings are the psField stamp that the world was built with), against the planted stamp; '
          'plain: the diag-style stack of the same stars. Held-out: the truth stars farther than %g px from every in-situ star, fitted in the high-passed image (amplitude and shift; the model H applied to the shifted stamp): '
          'median over the stars of the rms of the residual in the %d x %d box over the fitted peak of the model (the misfit, %%), and of the standardized rms, with the psField stamp and with the in-situ stamp, and the number of stars for which the in-situ misfit is lower.' % (HELD_RAD, 2 * FIT_BOX + 1, 2 * FIT_BOX + 1), '']
    rows = []
    for name, pxy in TARGETS:
        k = 'ST|patch|%s|' % pkey(name, pxy)
        if k + 'world|n_sel' not in raw:
            continue
        vw, pw = raw[k + 'world|vs_planted'], raw[k + 'world|plain']
        h = raw[k + 'held']
        hh = h[np.isfinite(h).all(1)] if len(h) else h.reshape(0, 4)
        rows.append([name, patch_label(pxy), '%d' % raw[k + 'world|n_sel'][0], pct(vw[0]), num(vw[1]), num(vw[2]), num(vw[3] * 100), pct(pw[0]), '%d' % len(hh),
                     num(float(np.median(hh[:, 1])) * 100, 1) if len(hh) else '-', num(float(np.median(hh[:, 3])) * 100, 1) if len(hh) else '-', '%d of %d' % (int((hh[:, 3] < hh[:, 1]).sum()), len(hh)) if len(hh) else '-',
                     num(float(np.median(hh[:, 0])), 1) if len(hh) else '-', num(float(np.median(hh[:, 2])), 1) if len(hh) else '-'])
    L += [X.table(['set', 'patch', 'world: stars', 'recovered against planted, central pixel (%)', 'FWHM recovered', 'FWHM planted', 'rms difference (% of peak)', 'plain stack against planted (%)',
                   'held-out stars', 'median misfit / peak (%), psField', 'median misfit / peak (%), in-situ', 'in-situ lower', 'median standardized rms, psField', 'median standardized rms, in-situ'], rows), '']
    L += reading_lines(raw)
    return L


def reading_lines(raw):
    """### 1.e: the validation summarized per group by the script (patch scope, with the pooled-by-frame truth comparison)."""
    L = ['### 1.e Reading of the validation (computed from the tables above)', '']
    for title, name in GROUPS:
        ps = [p for n, p in TARGETS if n == name and 'ST|patch|%s|vs_psf' % pkey(n, p) in raw]
        if not ps:
            continue
        vs = np.array([raw['ST|patch|%s|vs_psf' % pkey(name, p)] for p in ps])
        wd = np.array([raw['ST|patch|%s|world|vs_planted' % pkey(name, p)] for p in ps])
        held = [raw['ST|patch|%s|held' % pkey(name, p)] for p in ps]
        held = np.concatenate([h[np.isfinite(h).all(1)] for h in held]) if held else np.zeros((0, 4))
        k = 'ST|frame|%s|all|' % name
        txt = '* %s: in-situ stamp against psField, central pixel %s %% on average over the patches (FWHM %s against %s px)' % (title, pct(float(np.nanmean(vs[:, 0]))), num(float(np.nanmean(vs[:, 1]))), num(float(np.nanmean(vs[:, 2]))))
        if k + 'truth_n' in raw and np.isfinite(raw[k + 'truth_vs_psfmean'][0]):
            txt += '; the truth-star stamp pooled over the frame (%d stars) %s %%, the in-situ frame stamp %s %%, in-situ against truth-star %s %%' % (
                raw[k + 'truth_n'][0], pct(raw[k + 'truth_vs_psfmean'][0]), pct(raw[k + 'insitu_vs_psfmean'][0]), pct(raw[k + 'insitu_vs_truth'][0]))
        txt += '; recovery of a planted stamp in the same crowding: central pixel %s %% on average, at worst %s %%' % (pct(float(np.nanmean(wd[:, 0]))), pct(float(wd[np.nanargmax(np.abs(wd[:, 0])), 0])))
        if len(held):
            txt += '; held-out stars: the in-situ stamp has the lower misfit on %d of %d' % (int((held[:, 3] < held[:, 1]).sum()), len(held))
        L.append(txt + '.')
    return L + ['']


def f_rows_section(raw):
    L = ['## 2. F at mag < 22 with each stamp', '']
    have = {s: any(k.startswith('RUN|%s|' % s) for k in raw) for s in ('patch', 'frame', 'scr')}
    if not have['patch'] and not have['frame']:
        return L + ['(no runs in the raw yet)', '']
    stamps = ['psf'] + [s for s in ('patch', 'frame', 'scr') if have[s]]
    L += ['The psField rows are the frozen results (guarded runs of `cyclefix_raw.npz` for the ISS, `baselines_hp2_raw.npz` for crowdsource+HP and photutils+HP, `m2_raw.npz` for the crowdsource default), copied by `--ref`. '
          'A cell is F [95 % interval of the 20 x 20 px block bootstrap, 2000 resamples, seed 0]; a row computed from a path that reached the entry cap is left out of every pool and named under the pool.', '']
    L += ['### 2.1 Pooled by group', '']
    rows, stalled, missing = [], [], {}
    for title, name in GROUPS:
        for meth in METHODS:
            cells = []
            for s in stamps:
                g = pool(raw, name, [(s, meth)])
                for p, bad in g['left']:
                    if blk_of(raw, s, name, p, meth) is not None:
                        stalled.append('%s %s, %s, %s' % (name, patch_label(p), STAMP_LABEL[s], LABEL[meth]))
                    else:
                        missing[(STAMP_LABEL[s], LABEL[meth])] = missing.get((STAMP_LABEL[s], LABEL[meth]), 0) + 1
                if not g['patches']:
                    cells.append('-')
                    continue
                idx = X.boot_idx_strat(g['sizes'])
                f, ci = X.boot_cpf(g['blk'][(s, meth)], idx)[2]
                cells.append('%s (%d)' % (X.fmt_iv(f, ci), len(g['patches'])))
            rows.append([title if meth == METHODS[0] else '', LABEL[meth]] + cells)
    L += [X.table(['group', 'method'] + [STAMP_LABEL[s] + ' (patches)' for s in stamps], rows), '']
    L += ['Left out of a pool because the ISS path reached the entry cap (a defective run): %s.' % ('; '.join(stalled) if stalled else 'none'), '',
          'Rows with no run (patches missing from the pool): %s.' % ('; '.join('%s, %s: %d' % (st_, me_, n_) for (st_, me_), n_ in missing.items()) if missing else 'none'), '']
    L += ['### 2.2 Per patch', '']
    for meth in METHODS:
        rows = []
        for name, pxy in TARGETS:
            cells = []
            for s in stamps:
                b = blk_of(raw, s, name, pxy, meth)
                if b is None:
                    cells.append('-')
                    continue
                f, ci = X.boot_cpf(b, X.boot_idx_strat([25]))[2]
                cells.append(X.fmt_iv(f, ci) + (' STALLED' if is_stalled(raw, s, name, pxy, meth) else ''))
            rows.append([name, patch_label(pxy)] + cells)
        L += ['**%s**' % LABEL[meth], '', X.table(['set', 'patch'] + [STAMP_LABEL[s] for s in stamps], rows), '']
    L += ['### 2.3 Paired differences ISS minus baseline, the same stamp for both methods, pooled blocks', '',
          'First minus second on the blocks of the patches that have both rows valid; win / loss when the 95 % interval excludes 0.', '']
    for title, name in GROUPS:
        rows = []
        for a, b in PAIRS:
            cells = []
            for s in stamps:
                g = pool(raw, name, [(s, a), (s, b)])
                if not g['patches']:
                    cells.append('-')
                    continue
                d, dp, df = X.boot_diff(g['blk'][(s, a)], g['blk'][(s, b)], X.boot_idx_strat(g['sizes']))
                cells.append(diff_cell(df[0], df[1]) + ' (%d)' % len(g['patches']))
            rows.append([PAIR_NAME[(a, b)]] + cells)
        L += ['**%s**' % title, '', X.table(['comparison'] + [STAMP_LABEL[s] + ' (patches)' for s in stamps], rows), '']
    L += ['### 2.4 The effect of the stamp, in-situ minus psField, per method (paired, same blocks)', '',
          'The control column (if present) is the psField stamp plus a random permutation of the in-situ perturbation (Amendment 1): the same pixel noise, none of the shape.', '']
    for title, name in GROUPS:
        rows = []
        for meth in METHODS:
            cells = []
            for s in stamps[1:]:
                g = pool(raw, name, [(s, meth), ('psf', meth)])
                if not g['patches']:
                    cells.append('-')
                    continue
                d, dp, df = X.boot_diff(g['blk'][(s, meth)], g['blk'][('psf', meth)], X.boot_idx_strat(g['sizes']))
                cells.append(diff_cell(df[0], df[1]) + ' (%d)' % len(g['patches']))
            rows.append([LABEL[meth]] + cells)
        L += ['**%s**' % title, '', X.table(['method'] + [STAMP_LABEL[s] + ' minus psField (patches)' for s in stamps[1:]], rows), '']
    L += ['### 2.5 The ISS-minus-crowdsource gap: psField against in-situ, and the change of the gap', '',
          'Gap = F(ISS) minus F(crowdsource+HP) with the same stamp for both. "Change" = (gap with the in-situ stamp) minus (gap with the psField stamp), paired on the same blocks and resamples.', '']
    verdicts = []
    for title, name in GROUPS:
        rows = []
        for a, b in (('iso', 'cs_hp'), ('adq', 'cs_hp'), ('iso', 'phu_hp'), ('adq', 'phu_hp'), ('iso', 'cs_def'), ('adq', 'cs_def')):
            for s in stamps[1:]:
                g = pool(raw, name, [(s, a), (s, b), ('psf', a), ('psf', b)])
                if not g['patches']:
                    continue
                idx = X.boot_idx_strat(g['sizes'])
                d0 = X.boot_diff(g['blk'][('psf', a)], g['blk'][('psf', b)], idx)[2]
                d1 = X.boot_diff(g['blk'][(s, a)], g['blk'][(s, b)], idx)[2]
                dc, dci = boot_dd(g['blk'][(s, a)], g['blk'][(s, b)], g['blk'][('psf', a)], g['blk'][('psf', b)], idx)
                v0, v1, vd = X.verdict(*d0), X.verdict(*d1), X.verdict(dc, dci)
                rows.append([PAIR_NAME[(a, b)], STAMP_LABEL[s], len(g['patches']), diff_cell(*d0), diff_cell(*d1), diff_cell(dc, dci), gap_sentence(v0, v1, vd)])
                verdicts.append((name, a, b, s, v0, v1, vd))
        L += ['**%s**' % title, '', X.table(['comparison', 'in-situ stamp', 'patches', 'gap, psField', 'gap, in-situ', 'change of the gap', 'reading'], rows), '']
    return L, verdicts


def answer_section(verdicts):
    L = ['## 3. The answer, as the tables give it', '']
    if not verdicts:
        return L + ['(no runs yet)', '']
    for title, name in GROUPS:
        for s in ('patch', 'frame'):
            for a, b in (('iso', 'cs_hp'), ('adq', 'cs_hp')):
                v = [x for x in verdicts if x[0] == name and x[1] == a and x[2] == b and x[3] == s]
                if v:
                    _, _, _, _, v0, v1, vd = v[0]
                    L.append('* %s, %s stamp, %s: psField verdict %s, in-situ verdict %s; %s.' % (title, STAMP_LABEL[s].split('(')[1].rstrip(')'), PAIR_NAME[(a, b)], v0, v1, gap_sentence(v0, v1, vd)))
    return L + ['']


def runs_section(raw):
    rows = []
    for scope in ('patch', 'frame'):
        for name, pxy in TARGETS:
            for meth in METHODS:
                pre = 'RUN|%s|%s|%s' % (scope, pkey(name, pxy), meth)
                if pre + '|wall' in raw:
                    secs = raw[pre + '|secs']
                    if meth in ('iso', 'adq'):
                        tr = raw[pre + '|trace|counts']
                        rows.append([scope, name, patch_label(pxy), LABEL[meth], str(raw[pre + '|kernel'][0]), '%d' % tr[0], '%d' % tr[6], 'STALLED' if is_stalled(raw, scope, name, pxy, meth) else '-',
                                     '%d of %d' % tuple(raw[pre + '|bic'][:2]), '%.0f' % secs[0]])
                    else:
                        rows.append([scope, name, patch_label(pxy), LABEL[meth], '-', '%d' % raw[pre + '|det'].shape[0], '-', '-', '-', '%.0f' % secs[0]])
                elif pre + '|same' in raw:
                    rows.append([scope, name, patch_label(pxy), LABEL[meth], str(raw[pre + '|kernel'][0]) + ' (same as iso)', '-', '-', '-', '-', '-'])
    L = ['## 4. Run records', '', 'Entries (detections for the baselines), the guard\'s activations (lattice points barred), the BIC block, seconds (one core, shared machine).', '']
    return L + ([X.table(['scope', 'set', 'patch', 'run', 'kernel', 'entries / detections', 'guard activations', 'stalled', 'BIC block', 'seconds'], rows)] if rows else ['(none yet)']) + ['']


def csb_section(raw):
    L = ['## 5. Task B: crowdsource options on the development patches', '']
    tags = [t for t, _, _, _ in cs3_variants()] + ['cs3_combo']
    if not any(k.startswith('CSD|') for k in raw):
        return L + ['(not run yet)', '']
    L += ['Base: the round-2 selection `%s` (high-pass 2 px, threshold 3, miniter 8, maxiter 20, fewstars 20, psfderiv, blendthreshu 2, titer_thresh 2, psfvalsharpcutfac 0.7, psfsharpsat 0.7, the ISS weights). '
          'One option at a time (stage 1), then the declared combination (stage 2: every option that beats the base by %.3f or more). Pooled F at mag < 22 on the development patches %s.' % (
              CS_BASE_TAG, COMBO_GAIN, ', '.join('(%d, %d)' % p for p in cs_dev_patches())), '']
    base = dev_f(raw, 'cs3_base')
    chk = [int(raw['CSD|check|x%d_y%d' % p][0]) for p in cs_dev_patches() if 'CSD|check|x%d_y%d' % p in raw]
    if chk:
        L += ['Check: the base reproduces the stored development totals of `%s` (`baselines_hp2_raw.npz`) on %d of %d patches.' % (CS_BASE_TAG, sum(chk), len(chk)), '']
    rows = []
    desc = {t: (f, lv) for t, f, lv, _ in cs3_variants()}
    for t in tags:
        f = dev_f(raw, t)
        failed = [p for p in cs_dev_patches() if 'CSD|%s|x%d_y%d|failed' % (t, p[0], p[1]) in raw]
        if not np.isfinite(f) and not failed:
            continue
        ks = [sum(raw['CSD|%s|x%d_y%d|tot' % (t, p[0], p[1])] for p in cs_dev_patches())] if np.isfinite(f) else None
        c, pu = (X.f_measure(ks[0])[1], X.f_measure(ks[0])[2]) if ks else (np.nan, np.nan)
        fam, lv = desc.get(t, (None, None))
        what = 'base' if t == 'cs3_base' else ('combination: ' + ', '.join(str(x) for x in raw['CSD|SEL|combo']) if t == 'cs3_combo' and 'CSD|SEL|combo' in raw else '%s = %s' % (fam, lv))
        rows.append([t, what, num(c, 3), num(pu, 3), num(f, 3), ('%+.3f' % (f - base)) if np.isfinite(f) and np.isfinite(base) else '-', ('FAILED: ' + str(raw['CSD|%s|x%d_y%d|failed' % (t, failed[0][0], failed[0][1])][0])) if failed else ''])
    L += [X.table(['variant', 'option', 'completeness (mag < 22)', 'purity (mag < 22)', 'F', 'F minus base', 'note'], rows), '']
    if 'CSD|SEL|best' in raw:
        best, bf, base_f = str(raw['CSD|SEL|best'][0]), float(raw['CSD|SEL|best_F'][0]), float(raw['CSD|SEL|base_F'][0])
        L += ['Selected (largest pooled development F; a tie to the base): **%s**, F %.3f against the round-2 selection %.3f (difference %+.3f). Combination used in stage 2: %s.' % (
            best, bf, base_f, bf - base_f, ', '.join(str(x) for x in raw['CSD|SEL|combo'])), '']
        if int(raw['CSD|SEL|finding'][0]):
            L += ['**A finding** (more than %.3f above the round-2 selection on the development patches). Scored once on the M2 test patches, post hoc:' % FINDING_GAIN, '']
            rows = []
            for p in X.TESTS:
                cells = []
                for t in ('cs3_base', best):
                    k = 'CSTEST|%s|x%d_y%d|tot' % (t, p[0], p[1])
                    cells.append(num(X.f_measure(raw[k])[0], 3) if k in raw else '-')
                rows.append(['(%d, %d)' % tuple(p)] + cells)
            tot = {}
            for t in ('cs3_base', best):
                ks = ['CSTEST|%s|x%d_y%d|tot' % (t, p[0], p[1]) for p in X.TESTS]
                tot[t] = num(X.f_measure(sum(raw[k] for k in ks))[0], 3) if all(k in raw for k in ks) else '-'
            rows.append(['pooled counts', tot['cs3_base'], tot[best]])
            L += [X.table(['M2 test patch', 'round-2 selection (base) F', '%s F' % best], rows), '']
        else:
            L += ['Not a finding (the best variant does not exceed the round-2 selection by more than %.3f pooled development F): the test patches were not scored.' % FINDING_GAIN, '']
    return L


def limits_section():
    return ['## 6. Limitations', '',
            '* POST HOC. The study follows the transfer test, the kernel study, the diagnostics and the baselines\' two rounds; the hypothesis came from the diagnostics. The estimator was chosen on the M2 development patches and their synthetic worlds, '
            'and its star-selection radii after a label-free look at the star counts and the stack on the target patches (docstring); no F on an in-situ stamp had been computed. The truth-based numbers of section 1 are validation, not inputs.',
            '* FEW STARS. A 100 x 100 patch has 4 to 10 stars that pass the isolation rules in these crowded fields; the patch-scope stamp is noisy (intervals in 1.a) and its noise enters every method. The frame scope pools the patches of one frame, '
            'at the price of assuming one PSF across the frame.',
            '* The estimator corrects the alignment bias of the centring but NOT crowding or noise; section 1.d shows what it does to a stamp that is known, in the patch\'s own crowding. The in-situ stamp changes the core (radius 9 px), the flux within 6 px and the wings stay psField\'s.',
            '* The stamp for the ISS goes through `fit_psf_parts` (a three-component isotropic mixture) or, where the rule `fit_adequate` (2e-3 of the peak, not re-tuned) rejects it, the general parts; a noisy in-situ stamp can push the rule to a general kernel that fits the noise.',
            '* SENSITIVITY. A method whose F moves by about as much under the scrambled-perturbation control as under the in-situ stamp (section 2.4) is not measuring the stamp\'s shape: its F on these patches is a fragile function of the stamp\'s pixels.',
            '* The ISS configuration, the baselines\' selected variants and the high-pass are those frozen/selected with the psField stamp (nothing is re-tuned for the in-situ stamp); crowdsource and photutils were given the stamp through the same high-pass as the ISS.',
            '* Windows of one frame each (M53 four, M92 two, M2 four): intervals are conditional on the patches, groups are not independent trials. The M92 group has two patches.',
            '* The psField rows are the frozen results: the ISS rows are the guarded re-runs of `cyclefix_raw.npz` (identical to the stored ones where the guard never acted), the in-situ ISS rows are guarded runs of the same functions.', '']


def report(raw):
    meta = ''
    if 'meta|files' in raw:
        meta = ' Code the runs used (SHA-256, first 12 digits): %s. Versions: %s.' % ('; '.join('%s %s' % (f, h[:12]) for f, h in zip((str(x) for x in raw['meta|files']), (str(x) for x in raw['meta|sha256']))), ', '.join(str(v) for v in raw['meta|versions']))
    L = ['# An in-situ PSF for every method: does the stamp-versus-star mismatch explain the ISS\'s loss to crowdsource on M53 and M92?', '',
         'Generated by `experiments/exp_astro_insitu.py` from `insitu_raw.npz` alone; every number below is printed by the script. **POST HOC**: the design and every analysis were declared in the docstring of the script before any '
         'method was run on an in-situ stamp; the stamp is label-free (the HST truth is never an input), its validation (section 1) uses the truth. The frozen results, the freeze records and the selections of the baselines are untouched.%s' % meta, '']
    L += stamps_section(raw)
    fr = f_rows_section(raw)
    if isinstance(fr, tuple):
        L += fr[0]
        L += answer_section(fr[1])
    else:
        L += fr
    L += runs_section(raw)
    L += csb_section(raw)
    L += limits_section()
    return '\n'.join(L).replace('\n\n\n', '\n\n').rstrip('\n') + '\n'


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--stamps', action='store_true')
    ap.add_argument('--scr', action='store_true')
    ap.add_argument('--ref', action='store_true')
    ap.add_argument('--run', action='store_true')
    ap.add_argument('--scope', default='patch', choices=SCOPES)
    ap.add_argument('--methods', default=','.join(METHODS))
    ap.add_argument('--sets', default='m2,m53,m92')
    ap.add_argument('--cs-dev', action='store_true')
    ap.add_argument('--stage', type=int, default=1, choices=(1, 2))
    ap.add_argument('--cs-test', action='store_true')
    ap.add_argument('--from-raw', action='store_true')
    a = ap.parse_args()
    if a.stamps:
        part_stamps()
    elif a.scr:
        part_scr()
    elif a.ref:
        part_ref()
    elif a.run:
        part_run(a.scope, [m for m in a.methods.split(',') if m], a.sets.split(','))
    elif a.cs_dev:
        part_cs_dev(a.stage)
    elif a.cs_test:
        part_cs_test()
    elif a.from_raw:
        os.makedirs(OUT, exist_ok=True)
        with open(MD, 'w') as f:
            f.write(report(load_raw()))
        print('wrote', MD)
    else:
        ap.error('give --stamps, --scr, --ref, --run, --cs-dev, --cs-test or --from-raw')


if __name__ == '__main__':
    main()
