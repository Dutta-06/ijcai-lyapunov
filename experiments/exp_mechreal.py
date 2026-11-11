"""
A pre-registered REAL-DATA test of the mechanism: the factorised network trained on SDSS star-field windows, scored against HST, with the two conditions switched off in turn
===========================================================================================================================================================================

Plan (committed before any test window was processed by any arm): results/astro/mechreal_plan.md. Freeze: results/astro/mechreal_freeze.json. Report: results/astro/mechreal.md, printed
from results/astro/mechreal_raw.npz alone (`--from-raw`).

What is trained. `experiments/exp_network_imaging.py` trains the paper's network on SIMULATED frames. Here the same network (same kernel family, same loss, same solver; `mechreal_core.py`,
pinned to the harness by tests) is trained on 20 x 20 px windows cut from the SDSS r-band patches of the astronomy experiments, whose truth is the HST (HUGS) catalogue:
    kernel   the pixel-integrated three-Gaussian fit of the frame's psField stamp (`imaging.fit_psf_parts(stamp, 3, 1.0)`, as `exp_astro_m2.run_path`), high-passed with the
             border-consistent filter of width 3 px; start lattice of spacing 1 px with a margin of 1 px (22 x 22 = 484 units);
    loss     F = 1/2 sum_p W_p (H y - sum_i w_i H phi_i)_p^2 with the pipeline's noise law W = 1 / (smoothed data / gain + dark variance), computed on the 100 x 100 patch and FIXED;
    units    the data are divided by the flux limit f_lim = Z_LIM / scale, scale = `imaging.noise_scale` of the window (the flux of matched-filter significance Z_LIM): a counted
             star has amplitude >= 1, the noise scale of the dual is Z_LIM and the time stop of significance kappa is tau = 1 / (kappa Z_LIM).
Arms (same dictionary, same loss): A the paper's unit by gradient flow from eps (the mechanism); B depth 2 (w = u^2 - v^2); C depth 1 (w = u); D the paper's unit by Adam (default
settings); E the paper's unit by Adam with one shared second moment; F the computed sliding path (`imaging_guard.sliding_path_guarded`, global refits, no exclusion zone, fixed widths
and weights: Definition 1 restricted to non-negative amplitudes), the reference.

Truth. The HUGS stars of the patch (F606W as r, `astro_common.mag_to_counts`), shifted by the patch's registered rigid offset where the registry has one; stars closer than R_BLEND
are one source (summed flux, flux-weighted centre); a source counts when its flux is >= f_lim and its centre is in the window interior (EDGE px from every side).

Outcomes of a snapshot (all arms alike): the formed units |w| > THETA (in f_lim), same-sign formed units within R_MERGE merged (single linkage, summed amplitude, |w|-weighted centre), those
in the interior counted: K_hat; count error |K_hat - K_true|; detection F at R_MATCH (one-to-one, maximum cardinality then minimum distance).

    python experiments/exp_mechreal.py --select dev|test                    the windows (truth geometry only), printed and stored in the raw
    python experiments/exp_mechreal.py --path dev|test                      arm F and the snapshot times of every window
    python experiments/exp_mechreal.py --pack dev|test --out DIR [--arms A,B] [--eps 1e-8] [--chunks 4]     the job files of the network arms (small derived arrays only)
    python experiments/exp_mechreal.py --work JOBFILE --out RESULTFILE [--procs 1]      run a job file (local; the remote worker runs the same function)
    python experiments/exp_mechreal.py --collect dev|test RESULTFILE[,...]  merge result files into the raw (write-once for the test stage)
    python experiments/exp_mechreal.py --freeze                             write the freeze record from the development raw
    python experiments/exp_mechreal.py --from-raw                           regenerate results/astro/mechreal.md
Environment: MECHREAL_OUT (directory of the raw, the md and the freeze; default results/astro).
"""
import argparse
import base64
import hashlib
import io
import json
import os
import sys
import time

os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.sparse.csgraph import connected_components
from scipy.spatial.distance import cdist

HERE = os.path.abspath(os.path.dirname(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
sys.path.insert(0, REPO)
import mechreal_core as MC                      # noqa: E402
from src.core import imaging as IM              # noqa: E402
from src.core import imaging_guard as IG        # noqa: E402

ASTRO = os.path.join(REPO, 'results', 'astro')
OUT = os.environ.get('MECHREAL_OUT', ASTRO)
RAW, MD, FREEZE = os.path.join(OUT, 'mechreal_raw.npz'), os.path.join(OUT, 'mechreal.md'), os.path.join(OUT, 'mechreal_freeze.json')
PLAN = os.path.join(ASTRO, 'mechreal_plan.md')
FILES = ('experiments/exp_mechreal.py', 'experiments/mechreal_core.py', 'src/core/imaging.py', 'src/core/imaging_guard.py', 'experiments/astro_common.py', 'experiments/astro_fresh.py',
         'experiments/exp_astro_m2.py', 'results/astro/mechreal_plan.md', 'results/astro/clusterdep_patches.json', 'results/astro/fresh_patches.json')

# ---- the registry of everything the plan fixes ------------------------------------------------------------------------------------------------------------------
NW, HP, EDGE = 20, 3.0, 2.5                       # window side (px), high-pass width (px), border excluded from the scoring (px; `astro_common.EDGE`)
Z_LIM = 40.0                                      # the flux limit in matched-filter significance (f_lim = Z_LIM / noise scale of the window)
Z_BRIGHT = 50.0                                   # no truth star brighter than Z_BRIGHT f_lim inside the window padded by PAD px (unsaturated, moderate dynamic range)
PAD = 5.0
K_MIN, K_MAX = 2, 12                              # counted truth sources in the interior of an eligible window
R_BLEND, Z_FAINT = 0.5, 0.1                       # truth stars closer than R_BLEND px are one source; stars fainter than Z_FAINT f_lim are ignored in the blend
CAND_ORIGINS = (10, 30, 50, 70)                   # the candidate windows of a 100 x 100 patch: the 16 non-overlapping windows with these origins in x and y
SEED_SELECT = 20261110
N_DEV, N_TEST_PER_PATCH = 8, 5
DEV_PATCHES = (('M2', 'm2', (110, 630)), ('M2', 'm2', (210, 630)))
FAILED_REG = (('ngc5904', (640, 240)), ('ngc6254', (1900, 80)), ('ngc6934', (370, 1230)))      # clusterdep.md section 3.3: the registration fit failed
R_MERGE, R_MATCH = 0.5, 0.5                       # px
R_MERGE_SENS = (0.25, 1.0, 1.5)
KAPPAS = (40.0, 30.0, 20.0, 15.0, 10.0)           # the grid of the time stop (development)
THETAS = (0.25, 0.5, 0.7, 0.85, 1.0)              # the grid of the formed threshold (development)
EPS_LIST = (1e-6, 1e-8)
STEPS = (1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000)      # the snapshots of the discrete optimizers
LR_ADAM = 1e-3
S_MAX = 100000                                    # the training budget of the discrete optimizers D and E (steps): their primary snapshot
W_AWAKE = 0.05                                    # a unit is 'awake' when |w| > W_AWAKE (the smallest amplitude every arm's stored snapshot keeps)
C_EARLY = (1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1)   # arm C: extra early snapshots at these fractions of the stop time (it has no wall, hence no clock of its own)
PATH_TOL, PATH_NFEV, PATH_MAX = 1e-12, 2000, 150
N_BOOT, BOOT_SEED = 20000, 20261110


def utc():
    return time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


# ---------------------------------------------------------------------------------------------------------------------------------------------------------------
# pure parts: truth sources, the formed set, the outcomes
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------
def in_interior(c, nw=NW, edge=EDGE):
    """Mask of the positions (k, 2) = (x, y) at least `edge` px from every side of an nw x nw window (pixel centres at the integers; `astro_common.truth_for_patch`'s rule)."""
    c = np.atleast_2d(np.asarray(c, float)).reshape(-1, 2)
    lo, hi = edge - 0.5, nw - 0.5 - edge
    return (c[:, 0] >= lo) & (c[:, 0] <= hi) & (c[:, 1] >= lo) & (c[:, 1] <= hi)


def blend_sources(x, y, z, r_blend=R_BLEND, z_faint=Z_FAINT):
    """The truth sources of a star list (positions x, y; fluxes z in units of the flux limit): stars with z >= z_faint closer than r_blend (single linkage) are one source with
    the summed flux at the flux-weighted centre. Returns (c (k, 2), z (k,)) sorted by decreasing flux."""
    x, y, z = np.asarray(x, float), np.asarray(y, float), np.asarray(z, float)
    keep = z >= z_faint
    c, z = np.c_[x[keep], y[keep]], z[keep]
    if len(z) == 0:
        return np.zeros((0, 2)), np.zeros(0)
    n_comp, lab = connected_components(cdist(c, c) < r_blend, directed=False)
    zs = np.bincount(lab, weights=z, minlength=n_comp)
    cs = np.c_[np.bincount(lab, weights=z * c[:, 0], minlength=n_comp), np.bincount(lab, weights=z * c[:, 1], minlength=n_comp)] / zs[:, None]
    o = np.argsort(-zs, kind='stable')
    return cs[o], zs[o]


def counted_truth(src_c, src_z, nw=NW, edge=EDGE):
    """The counted truth: the sources with flux >= 1 (the flux limit) whose centre is in the interior. Returns their centres (k, 2)."""
    src_c = np.asarray(src_c, float).reshape(-1, 2)
    return src_c[(np.asarray(src_z, float) >= 1.0) & in_interior(src_c, nw, edge)]


def formed_units(c, w, theta, r_merge=R_MERGE):
    """The formed set of a snapshot: the units with |w| > theta; same-sign formed units within r_merge of each other (single linkage) are one unit with the summed amplitude at
    the |w|-weighted centre. Returns (c (k, 2), w (k,))."""
    c, w = np.asarray(c, float).reshape(-1, 2), np.asarray(w, float)
    f = np.flatnonzero(np.isfinite(w) & (np.abs(w) > theta))
    if len(f) == 0:
        return np.zeros((0, 2)), np.zeros(0)
    cf, wf = c[f], w[f]
    near = (cdist(cf, cf) < r_merge) & (np.sign(wf)[:, None] == np.sign(wf)[None, :])
    n_comp, lab = connected_components(near, directed=False)
    a = np.abs(wf)
    asum = np.bincount(lab, weights=a, minlength=n_comp)
    cm = np.c_[np.bincount(lab, weights=a * cf[:, 0], minlength=n_comp), np.bincount(lab, weights=a * cf[:, 1], minlength=n_comp)] / asum[:, None]
    return cm, np.bincount(lab, weights=wf, minlength=n_comp)


def match_pairs(det, tru, r_match=R_MATCH):
    """The number of one-to-one pairs (detection, truth) within r_match: maximum cardinality, then minimum total distance (`astro_common.match` without the magnitude gate)."""
    det, tru = np.asarray(det, float).reshape(-1, 2), np.asarray(tru, float).reshape(-1, 2)
    if len(det) == 0 or len(tru) == 0:
        return 0
    D = cdist(det, tru)
    ok = D <= r_match
    if not ok.any():
        return 0
    big = 4.0 * (min(D.shape) + 1) * r_match + 1.0
    r, q = linear_sum_assignment(np.where(ok, D, big))
    return int(ok[r, q].sum())


OUT_FIELDS = ('k_hat', 'k_true', 'err', 'signed', 'tp', 'f1', 'n_formed_all', 'n_neg', 'k_pos', 'n_awake')


def outcomes(c, w, truth_c, theta, r_merge=R_MERGE, r_match=R_MATCH, nw=NW, edge=EDGE):
    """The OUT_FIELDS vector of a snapshot (unit centres c (d, 2), amplitudes w (d) in units of the flux limit) against the counted truth centres: K_hat (merged formed units of
    either sign in the interior), K_true, |K_hat - K_true|, K_hat - K_true, matched pairs, F = 2 TP / (K_hat + K_true) (1 if both are 0), the merged formed units anywhere on
    the lattice, the negative ones among K_hat, the positive ones among K_hat, and the units (unmerged, anywhere) with |w| > W_AWAKE."""
    cm, wm = formed_units(c, w, theta, r_merge)
    inn = in_interior(cm, nw, edge) if len(wm) else np.zeros(0, bool)
    det, wd = cm[inn], wm[inn]
    truth_c = np.asarray(truth_c, float).reshape(-1, 2)
    kh, kt = len(det), len(truth_c)
    tp = match_pairs(det, truth_c, r_match)
    f1 = 1.0 if kh + kt == 0 else 2.0 * tp / (kh + kt)
    w = np.asarray(w, float)
    return np.array([kh, kt, abs(kh - kt), kh - kt, tp, f1, len(wm), int((wd < 0).sum()), int((wd > 0).sum()), int((np.abs(w[np.isfinite(w)]) > W_AWAKE).sum())], float)


def boot_mean(x, n_boot=N_BOOT, seed=BOOT_SEED):
    """(mean, lo, hi): the mean of x over the windows and its 95 % percentile bootstrap interval (resampling windows with replacement)."""
    x = np.asarray(x, float)
    if len(x) == 0:
        return float('nan'), float('nan'), float('nan')
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(x), (n_boot, len(x)))
    m = x[idx].mean(1)
    return float(x.mean()), float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def boot_mean_cluster(x, groups, n_boot=N_BOOT, seed=BOOT_SEED):
    """(mean, lo, hi) with whole groups (patches) resampled: the cluster bootstrap of the mean over windows."""
    x, groups = np.asarray(x, float), np.asarray(groups)
    labs = sorted(set(groups.tolist()))
    sums = np.array([x[groups == g].sum() for g in labs])
    cnts = np.array([(groups == g).sum() for g in labs], float)
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(labs), (n_boot, len(labs)))
    m = sums[idx].sum(1) / cnts[idx].sum(1)
    return float(x.mean()), float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def clopper_pearson(k, n, level=0.95):
    """The exact binomial interval of k successes in n trials."""
    from scipy.stats import beta
    if n == 0:
        return float('nan'), float('nan')
    a = (1.0 - level) / 2.0
    lo = 0.0 if k == 0 else float(beta.ppf(a, k, n - k + 1))
    hi = 1.0 if k == n else float(beta.ppf(1.0 - a, k + 1, n - k))
    return lo, hi


# ---------------------------------------------------------------------------------------------------------------------------------------------------------------
# data: patches, windows, problems
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------
def stage_patches(stage):
    """The patches of a stage: list of dict(label, cluster, pxy, shift, load, truth) with `load()` -> dict(dn, psf, gain, dark) and `truth()` -> dict(px, py, m606) (frame
    coordinates). dev: the two M2 development patches; test: the held-out patches of `clusterdep_patches.json` whose registration did not fail, and M3 / M15."""
    import astro_fresh as F
    import exp_astro_m2 as X
    out = []
    if stage == 'dev':
        for lab, ds, p in DEV_PATCHES:
            out.append(dict(label=lab, cluster='ngc7089', pxy=tuple(p), shift=(0.0, 0.0), load=(lambda p=p, ds=ds: X.load_patch(p[0], p[1], ds)), truth=(lambda ds=ds: X.truth_arrays(ds))))
        return out
    assert stage == 'test'
    for fn in ('clusterdep_patches.json', 'fresh_patches.json'):
        reg = F.load_registry(os.path.join(ASTRO, fn))
        for c in sorted(reg):
            for i, p in enumerate(reg[c]['patches']):
                if (c, tuple(p)) in FAILED_REG:
                    continue
                e = reg[c]
                out.append(dict(label=c, cluster=c, pxy=tuple(p), shift=tuple(float(s) for s in e['shifts'][i]), load=(lambda p=p, e=e: F.load_patch(p[0], p[1], e)),
                                truth=(lambda e=e: F.truth_arrays(e))))
    return out


def patch_model(pa):
    """The pipeline's kernel and noise law of a patch (`exp_astro_m2.run_path`): (sigma, parts (3, 2), sw (100, 100))."""
    from scipy.ndimage import gaussian_filter
    sigma, parts, _ = IM.fit_psf_parts(pa['psf'], n_parts=3, sigma0=1.0)
    dn = pa['dn']
    expected = np.maximum(gaussian_filter(dn, 1.5, mode='reflect'), np.percentile(dn, 1))
    return float(sigma), np.asarray(parts, float), 1.0 / np.sqrt(expected / pa['gain'] + pa['dark'])


def window_arrays(dn, sw, sigma, parts, wx, wy, nw=NW, hp=HP, z_lim=Z_LIM):
    """The normalised problem of the window with origin (wx, wy) of a patch: dict(y (nw, nw) = H(data) / f_lim, sw = sqrt weights x f_lim, f_lim (DN), scale (1 / DN))."""
    d = np.asarray(dn, float)[wy:wy + nw, wx:wx + nw]
    s = np.asarray(sw, float)[wy:wy + nw, wx:wx + nw]
    P = IM.PixelPSF(nw, nw, sigma, parts=tuple(map(tuple, np.asarray(parts, float))), hp=hp)
    lat = IM.Lattice(P, h=MC.LAT_H, margin=MC.LAT_MARGIN)
    scale = IM.noise_scale(lat, s ** 2)
    f_lim = z_lim / scale
    return dict(y=P.highpass(d - np.median(d)) / f_lim, sw=s * f_lim, f_lim=float(f_lim), scale=float(scale))


def window_truth(tr, pxy, shift, wx, wy, f_lim, nw=NW):
    """The truth of a window: the stars of the frame list `tr` (px, py, m606) shifted by `shift`, in window coordinates, fluxes in units of f_lim, within the window padded by
    PAD px. Returns dict(x, y, z)."""
    import astro_common as A
    x = np.asarray(tr['px'], float) + shift[0] - pxy[0] - wx
    y = np.asarray(tr['py'], float) + shift[1] - pxy[1] - wy
    m = np.asarray(tr['m606'], float)
    sel = (x >= -0.5 - PAD) & (x < nw - 0.5 + PAD) & (y >= -0.5 - PAD) & (y < nw - 0.5 + PAD) & (m < 30.0)
    return dict(x=x[sel], y=y[sel], z=A.mag_to_counts(m[sel]) / f_lim)


def eligible(truth, nw=NW):
    """(eligible, K_true, brightest z): the selection rule on the truth of a window (truth geometry only)."""
    c, z = blend_sources(truth['x'], truth['y'], truth['z'])
    k = len(counted_truth(c, z, nw))
    zmax = float(truth['z'].max()) if len(truth['z']) else 0.0
    return bool(K_MIN <= k <= K_MAX and zmax <= Z_BRIGHT), k, zmax


def select_windows(stage):
    """The windows of a stage, by the seeded geometric rule: per patch the eligible ones among the 16 candidates, in the order of a seeded permutation; dev: the first N_DEV / 2
    of each patch, test: the first N_TEST_PER_PATCH. Returns (windows, table lines); a window is dict(name, stage, patch (index), label, cluster, pxy, wxy, k_true, zmax, f_lim,
    scale, sigma, parts, y, sw, src_c, src_z, truth_c)."""
    patches = stage_patches('dev' if stage == 'devx' else stage)
    per = N_DEV // len(patches) if stage in ('dev', 'devx') else N_TEST_PER_PATCH
    wins, lines = [], []
    for ip, pt in enumerate(patches):
        pa, tr = pt['load'](), pt['truth']()
        sigma, parts, sw = patch_model(pa)
        rng = np.random.default_rng([SEED_SELECT, 1 if stage == 'test' else 0, ip])
        cands = [(wx, wy) for wy in CAND_ORIGINS for wx in CAND_ORIGINS]
        order = rng.permutation(len(cands))
        n_el, taken = 0, 0
        for j in order:
            wx, wy = cands[j]
            arr = window_arrays(pa['dn'], sw, sigma, parts, wx, wy)
            tw = window_truth(tr, pt['pxy'], pt['shift'], wx, wy, arr['f_lim'])
            ok, k, zmax = eligible(tw)
            n_el += ok
            if ok and stage == 'devx' and n_el <= per:
                continue                                                # devx: the eligible development candidates that 'dev' did not take (arm F only: tuning of the stop)
            if ok and (taken < per or stage == 'devx'):
                taken += 1
                sc, sz = blend_sources(tw['x'], tw['y'], tw['z'])
                name = '%s%s_p%d_x%d_y%d' % ({'dev': 'D', 'devx': 'X', 'test': 'T'}[stage], pt['label'], ip, wx, wy)
                wins.append(dict(name=name, stage=stage, patch=ip, label=pt['label'], cluster=pt['cluster'], pxy=pt['pxy'], wxy=(wx, wy), k_true=k, zmax=zmax, f_lim=arr['f_lim'],
                                 scale=arr['scale'], sigma=sigma, parts=parts, y=arr['y'], sw=arr['sw'], src_c=sc, src_z=sz, truth_c=counted_truth(sc, sz)))
        lines.append('%s patch %s: sigma %.3f px, %d of 16 candidates eligible, %d taken' % (pt['label'], pt['pxy'], sigma, n_el, taken))
    return wins, lines


def problem_of(win):
    return MC.make_problem(win['y'], win['sw'], win['sigma'], win['parts'], HP)


# ---------------------------------------------------------------------------------------------------------------------------------------------------------------
# arm F: the computed path
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------
def tau_of(kappa, z_lim=Z_LIM):
    """The time stop of significance kappa in the normalised problem (noise scale Z_LIM, wall 1): tau = 1 / (kappa Z_LIM)."""
    return 1.0 / (float(kappa) * z_lim)


def run_path(prob, tau_max, max_entries=PATH_MAX):
    """Arm F: `sliding_path_guarded` (dual order, global refits, no exclusion zone, fixed widths and weights) up to and including the first entry after tau_max. Returns
    dict(hist, taus, n_bar, capped)."""
    R = _path_until(prob, tau_max, max_entries)
    return dict(hist=R['hist'], taus=np.asarray(R['taus'], float), n_bar=len(R['barred']), capped=bool(len(R['taus']) >= max_entries))


def _path_until(prob, tau_max, max_entries):
    """The guarded path run with growing max_entries until its last entry is later than tau_max (so that the block at tau_max and the next entry time are both known), or it ends."""
    M = 30
    while True:
        R = IG.sliding_path_guarded(prob['P'], prob['y'], prob['sw'], prob['lat'], rho=1.0, local=False, r_ex=0.0, order='dual', max_entries=min(M, max_entries), tol=PATH_TOL,
                                    max_nfev=PATH_NFEV, passes=1, cycle_guard=True)
        taus = np.asarray(R['taus'], float)
        if len(taus) < min(M, max_entries) or (len(taus) and taus[-1] > tau_max) or M >= max_entries:
            return R
        M *= 2


def block_loss(prob, c, w):
    """F = 1/2 sum_p W (y - model)^2 of a block (the training loss of the arms)."""
    return 0.5 * float(np.sum(prob['W'] * (prob['y'] - prob['P'].synth(np.asarray(c, float).reshape(-1, 2), np.asarray(w, float))) ** 2))


def path_block(path, tau):
    """The block of the path at time tau (the block after the last entry with tau_entry <= tau): (c (k, 2), w (k), tau_K, tau_next); tau_K = 0 for the empty block, tau_next = inf
    if the path has no later entry."""
    taus = path['taus']
    i, K = IM.stop_time(path['hist'], taus, tau)
    h = path['hist'][i]
    later = taus[taus > tau]
    return h['c'][:, :2].copy(), h['w'].copy(), (float(taus[K - 1]) if K else 0.0), (float(later.min()) if len(later) else float('inf'))


TAU_CAP, Q_MID = 4.0, 0.5                         # the harness's snapshot between two entries: tau_mid = tau_K + Q_MID (min(tau_next, TAU_CAP tau_K) - tau_K)
REF_FIELDS = ('k', 'tau_K', 'tau_next', 'tau_mid', 'thr', 'polished', 'k_signed', 'n_neg_signed', 'thr_signed', 'polished_signed')


def path_refs(prob, path, kappa):
    """The references of the path agreement at the stop of significance kappa, as the harness takes them: the path's block at tau = `tau_of(kappa)` polished to its stationary
    point (`exp_network_imaging.polish_block`), the harness's formed threshold of that block, the snapshot tau_mid between the block's last entry and the next one, and the
    SIGNED path's block (`exp_network_imaging.signed_path`, Definition 1) at tau_mid, polished. Returns dict(c, w, sc, sw_, ref = REF_FIELDS vector)."""
    import exp_network_imaging as NI
    c, w, tK, tn = path_block(path, tau_of(kappa))
    nan = float('nan')
    if len(w) == 0:
        return dict(c=c, w=w, sc=np.zeros((0, 2)), sw_=np.zeros(0), ref=np.array([0, 0.0, tn, nan, nan, 0, 0, 0, nan, 0], float))
    tau_end = min(tn, TAU_CAP * tK)
    tau_mid = tK + Q_MID * (tau_end - tK)
    cp, wp, acc, _, _ = NI.polish_block(prob['P'], prob['y'], prob['sw'], c, w)
    R = NI.signed_path(prob['P'], prob['y'], prob['sw'], prob['lat'], max_entries=PATH_MAX, tau_stop=tau_mid * (1.0 + 1e-9))
    h = R['hist'][max(j for j, q in enumerate(R['hist']) if q['tau'] <= tau_mid * (1.0 + 1e-12))]
    sc, sw_, sacc = h['c'], h['w'], False
    if len(sw_):
        sc, sw_, sacc, _, _ = NI.polish_block(prob['P'], prob['y'], prob['sw'], sc, sw_)
    ref = [len(wp), tK, tn, tau_mid, NI.formed_threshold(wp), float(acc), len(sw_), int((sw_ < 0).sum()), NI.formed_threshold(sw_) if len(sw_) else nan, float(sacc)]
    return dict(c=cp, w=wp, sc=sc, sw_=sw_, ref=np.array(ref, float))


AGREE_FIELDS = ('n_formed', 'n_unmatched', 'max_d', 'rel_model', 'amp_max', 'dormant', 'n_neg', 'c1', 'c2', 'c3', 'c4', 'set_ok', 'all_ok')


def agreement(prob, c_net, w_net, c_ref, w_ref, thr):
    """The AGREE_FIELDS vector of a network state against a reference block, by the harness's `compare`: C1 (same count, every formed unit within 0.5 px of a block unit), C2
    (one-to-one within 0.25 px), C3 (model and amplitudes within 1e-3), C4 (dormant units below 1e-3); set_ok = C1 and C2 (the formed SET is the block), all_ok = C1..C4."""
    import exp_network_imaging as NI
    if len(w_ref) == 0 or not np.isfinite(thr) or not np.all(np.isfinite(w_net)):
        return np.full(len(AGREE_FIELDS), np.nan)
    g = NI.compare(prob['P'], prob['W'], c_net, w_net, c_ref, w_ref, thr)
    g['set_ok'] = float(g['c1'] and g['c2'])
    g['all_ok'] = g['ok']
    return np.array([g[f] for f in AGREE_FIELDS], float)


# ---------------------------------------------------------------------------------------------------------------------------------------------------------------
# the raw file
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------
WIN_META = ('patch', 'wx', 'wy', 'k_true', 'zmax', 'f_lim', 'scale', 'sigma', 'px', 'py')


def save_part(path, arrays, write_once=False):
    """Merge `arrays` into the raw at `path` through a temporary file. write_once: refuse to replace an existing key (the test stage is run once)."""
    old = {}
    if os.path.exists(path):
        with np.load(path, allow_pickle=False) as z:
            old = {k: z[k] for k in z.files}
    if write_once:
        clash = sorted(set(old) & set(arrays))
        if clash:
            raise RuntimeError('write-once: %d keys exist already (first: %s)' % (len(clash), clash[0]))
    old.update({k: np.asarray(v) for k, v in arrays.items()})
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    tmp = path[:-4] + '.tmp.npz'
    np.savez_compressed(tmp, **old)
    os.replace(tmp, path)


def load_raw(path=None):
    path = path or RAW
    if not os.path.exists(path):
        return {}
    with np.load(path, allow_pickle=False) as z:
        return {k: z[k] for k in z.files}


def window_arrays_raw(win):
    pre = 'W|%s|' % win['name']
    meta = dict(patch=win['patch'], wx=win['wxy'][0], wy=win['wxy'][1], k_true=win['k_true'], zmax=win['zmax'], f_lim=win['f_lim'], scale=win['scale'], sigma=win['sigma'],
                px=win['pxy'][0], py=win['pxy'][1])
    return {pre + 'meta': np.array([meta[f] for f in WIN_META], float), pre + 'label': np.array([win['label']]), pre + 'y': win['y'], pre + 'sw': win['sw'],
            pre + 'parts': np.asarray(win['parts'], float), pre + 'src_c': win['src_c'], pre + 'src_z': win['src_z'], pre + 'truth_c': win['truth_c']}


def window_from_raw(raw, name):
    pre = 'W|%s|' % name
    m = dict(zip(WIN_META, raw[pre + 'meta']))
    return dict(name=name, stage={'D': 'dev', 'X': 'devx', 'T': 'test'}[name[0]], patch=int(m['patch']), label=str(raw[pre + 'label'][0]), pxy=(int(m['px']), int(m['py'])),
                wxy=(int(m['wx']), int(m['wy'])), k_true=int(m['k_true']), zmax=float(m['zmax']), f_lim=float(m['f_lim']), scale=float(m['scale']), sigma=float(m['sigma']),
                parts=raw[pre + 'parts'], y=raw[pre + 'y'], sw=raw[pre + 'sw'], src_c=raw[pre + 'src_c'], src_z=raw[pre + 'src_z'], truth_c=raw[pre + 'truth_c'])


def stage_names(raw, stage):
    k = 'names|' + stage
    return [str(n) for n in raw[k]] if k in raw else []


def eps_tag(eps):
    return '%.0e' % float(eps)


def run_key(name, arm, eps):
    return '%s|%s|%s' % (name, arm, eps_tag(eps))


# ---------------------------------------------------------------------------------------------------------------------------------------------------------------
# the stages
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------
def cmd_select(stage, raw_path):
    raw = load_raw(raw_path)
    if stage_names(raw, stage):
        raise RuntimeError('the %s windows are selected already' % stage)
    wins, lines = select_windows(stage)
    arrays = {'names|' + stage: np.array([w['name'] for w in wins]), 'select|%s|log' % stage: np.array(lines), 'select|%s|utc' % stage: np.array([utc()])}
    for w in wins:
        arrays.update(window_arrays_raw(w))
    save_part(raw_path, arrays, write_once=True)
    print('\n'.join(lines))
    for w in wins:
        print('%-28s K_true %2d  brightest %5.1f f_lim  f_lim %.0f DN' % (w['name'], w['k_true'], w['zmax'], w['f_lim']))


def path_arrays(name, prob, path, kappas, ref_kappas):
    """The raw keys of arm F for a window: the entry times, the block at every kappa of `kappas` ('P|name|k<kappa>|c', 'w', 'tt' = (tau_K, tau_next)), the references of the
    path agreement at `ref_kappas` ('P|name|ref<kappa>|c', 'w', 'sc', 'sw', 'ref')."""
    pre = 'P|%s|' % name
    out = {pre + 'taus': path['taus'], pre + 'info': np.array([path['n_bar'], float(path['capped'])], float)}
    for k in kappas:
        c, w, tK, tn = path_block(path, tau_of(k))
        out.update({pre + 'k%g|c' % k: c, pre + 'k%g|w' % k: w, pre + 'k%g|tt' % k: np.array([tK, tn]), pre + 'k%g|loss' % k: np.array([block_loss(prob, c, w)])})
    for k in ref_kappas:
        r = path_refs(prob, path, k)
        out.update({pre + 'ref%g|c' % k: r['c'], pre + 'ref%g|w' % k: r['w'], pre + 'ref%g|sc' % k: r['sc'], pre + 'ref%g|sw' % k: r['sw_'], pre + 'ref%g|ref' % k: r['ref']})
    return out


def cmd_path(stage, raw_path, kappas, ref_kappas):
    raw = load_raw(raw_path)
    for name in stage_names(raw, stage):
        if 'P|%s|taus' % name in raw:
            miss = [k for k in kappas if 'P|%s|k%g|c' % (name, k) in raw and 'P|%s|k%g|loss' % (name, k) not in raw]
            if miss:                                                    # blocks stored before the loss was: the loss of the stored block
                prob = problem_of(window_from_raw(raw, name))
                save_part(raw_path, {'P|%s|k%g|loss' % (name, k): np.array([block_loss(prob, raw['P|%s|k%g|c' % (name, k)], raw['P|%s|k%g|w' % (name, k)])]) for k in miss})
            continue
        t0 = time.time()
        prob = problem_of(window_from_raw(raw, name))
        path = run_path(prob, tau_of(min(kappas)))
        save_part(raw_path, path_arrays(name, prob, path, kappas, ref_kappas), write_once=(stage == 'test'))
        print('F %-28s %3d entries, capped %d, barred %d, %.0f s' % (name, len(path['taus']), path['capped'], path['n_bar'], time.time() - t0), flush=True)


def flow_grid(raw, name, arm, kappas, ref_kappas):
    """The snapshot times (tau, ascending, unique) of a flow arm on a window: the time stops of `kappas`; for arm A also tau_mid of each reference kappa; for arm C also the early
    fractions C_EARLY of the latest time stop of the development grid (absolute times, the same on every window and stage)."""
    g = [tau_of(k) for k in kappas]
    if arm == 'A':
        for k in ref_kappas:
            tm = float(raw['P|%s|ref%g|ref' % (name, k)][REF_FIELDS.index('tau_mid')])
            if np.isfinite(tm):
                g.append(tm)
    if arm == 'C':
        g += [f * tau_of(min(DEV_FLOW_KAPPAS)) for f in C_EARLY]
    g = np.sort(np.array(g, float))
    return g[np.r_[True, np.diff(g) > 1e-12 * g[1:]]]


def cmd_pack(stage, raw_path, out_dir, arms, eps_list, kappas, ref_kappas, budget, chunks, steps=None):
    """Write the job files of the network arms: `chunks` files with the jobs dealt out round-robin (longest arms first), and store every run's snapshot grid in the raw."""
    raw = load_raw(raw_path)
    specs, grids = [], {}
    for name in stage_names(raw, stage):
        win = window_from_raw(raw, name)
        for arm in arms:
            for eps in eps_list:
                key = run_key(name, arm, eps)
                if 'R|%s|run' % key in raw:
                    continue
                grid = flow_grid(raw, name, arm, kappas, ref_kappas) if arm in MC.ARM_READOUT else np.array((steps or {}).get(arm, STEPS), float)
                grids['G|%s' % key] = grid
                specs.append(dict(key=key, arm=arm, eps=eps, grid=grid, budget=budget, lr=LR_ADAM, y=win['y'], sw=win['sw'], sigma=win['sigma'], parts=win['parts'], hp=HP))
    save_part(raw_path, grids)
    os.makedirs(out_dir, exist_ok=True)
    order = sorted(range(len(specs)), key=lambda i: ('ABCDE'.index(specs[i]['arm']), i))
    for c in range(chunks):
        MC.save_jobs(os.path.join(out_dir, 'jobs_%s_%d.npz' % (stage, c)), [specs[i] for i in order[c::chunks]])
    print('%d jobs in %d files under %s' % (len(specs), chunks, out_dir))


def cmd_collect(stage, raw_path, files):
    raw = load_raw(raw_path)
    names = set(stage_names(raw, stage))
    arrays = {}
    for f in files:
        with np.load(f, allow_pickle=False) as z:
            for k in z.files:
                if k.split('|')[1] in names:
                    arrays[k] = z[k]
    new = {k: v for k, v in arrays.items() if k not in raw}
    save_part(raw_path, new, write_once=(stage == 'test'))
    raw = load_raw(raw_path)
    fix = {}
    for k in [k for k in new if k.endswith('|n')]:                       # a run holds the snapshots of ITS job file: its grid is the first len(n) values of the packed grid
        g = 'G|' + k[2:-2]
        if g in raw and len(raw[g]) > len(new[k]):
            if stage == 'test':
                raise RuntimeError('test run %s has fewer snapshots than its grid' % k)
            fix[g] = raw[g][:len(new[k])]
    save_part(raw_path, fix)
    raw.update(fix)
    q = {}
    for name in sorted(names):
        for k in [k for k in new if k.startswith('R|%s|A|' % name) and k.endswith('|run')]:
            q.update(agreement_arrays(raw, name, float(k.split('|')[3])))
    save_part(raw_path, {k: v for k, v in q.items() if k not in raw}, write_once=(stage == 'test'))
    print('%d keys read, %d new, %d agreement vectors' % (len(arrays), len(new), len(q)))


# ---------------------------------------------------------------------------------------------------------------------------------------------------------------
# the agreement with the path (computed when results are collected; stored), the outcome tables, the development choices, the freeze
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------
M_E = 1.0                                         # P-O (E): the mean count error of E is within M_E of A's (the 95 % interval of the mean difference inside [-M_E, M_E])
SHARE_PATH = 0.60                                 # P-path: arm A forms the path's block on at least this share of the test windows
RATIO_C = 2.0                                     # P-R: "many more units": mean K_hat of C at least RATIO_C times that of A (and the interval of the mean difference above 0)
GATE = 0.5                                        # the test is informative iff the mean count error of the reference path F is at most GATE x mean K_true
EPS_MAIN = 1e-8
BUDGET_DEV, BUDGET_TEST = 3600.0, 3600.0             # wall-clock budget of one run (s): development, test (the same rule)
BIC_PAR = 3                                       # parameters per unit of the BIC stop (position and amplitude)


def grid_index(grid, value):
    """The index of `value` in a snapshot grid (-1 if it is not there to 1e-9 relative)."""
    grid = np.asarray(grid, float)
    if len(grid) == 0:
        return -1
    i = int(np.argmin(np.abs(grid - value)))
    return i if abs(grid[i] - value) <= 1e-9 * max(abs(value), 1e-300) else -1


def full_state(prob, snap):
    """The (c (d, 2), w (d)) of all units from a stored sparse snapshot (idx, w, c): the units not kept are at their lattice start with w = 0."""
    idx, w, c = snap
    cf, wf = prob['c0'].copy(), np.zeros(len(prob['c0']))
    cf[idx], wf[idx] = c, w
    return cf, wf


def agreement_arrays(raw, name, eps):
    """The stored path-agreement vectors of arm A on a window: for every reference kappa in the raw, 'Q|key|nn<kappa>' (against arm F's polished block) and 'Q|key|sg<kappa>'
    (against the signed path's block), each AGREE_FIELDS + (done,) at the snapshot tau_mid."""
    key = run_key(name, 'A', eps)
    r = MC.unpack_result(raw, key)
    out = {}
    if r is None:
        return out
    prob = None
    for k in sorted(q for q in raw if q.startswith('P|%s|ref' % name) and q.endswith('|ref')):
        kap = float(k.split('|')[2][3:])
        ref = dict(zip(REF_FIELDS, raw[k]))
        i = grid_index(raw['G|' + key], ref['tau_mid']) if np.isfinite(ref['tau_mid']) else -1
        if i < 0:
            continue
        prob = prob or problem_of(window_from_raw(raw, name))
        c, w = full_state(prob, r['snaps'][i])
        pre = 'P|%s|ref%g|' % (name, kap)
        out['Q|%s|nn%g' % (key, kap)] = np.r_[agreement(prob, c, w, raw[pre + 'c'], raw[pre + 'w'], ref['thr']), r['done'][i]]
        out['Q|%s|sg%g' % (key, kap)] = np.r_[agreement(prob, c, w, raw[pre + 'sc'], raw[pre + 'sw'], ref['thr_signed']), r['done'][i]]
    return out


def arm_value(arm, kappa, steps):
    """The primary snapshot of an arm: the time stop of significance kappa for the flows, the step budget for the discrete optimizers."""
    return tau_of(kappa) if arm in MC.ARM_READOUT else float(steps[arm])


def equal_fit_values(raw, names, arm, eps, kappa, values):
    """The equal-fit stop of an arm (label-free): per window the first reached snapshot among `values` whose training loss is at most the loss of arm F's block at the time stop
    of significance kappa; if none is, the last reached one. Returns dict name -> value (absent if the run has no reached snapshot)."""
    out = {}
    for name in names:
        key = run_key(name, arm, eps)
        r = MC.unpack_result(raw, key)
        if r is None or 'P|%s|k%g|loss' % (name, kappa) not in raw:
            continue
        target = float(raw['P|%s|k%g|loss' % (name, kappa)][0])
        last = None
        for v in sorted(values):
            i = grid_index(raw['G|' + key], v)
            if i < 0 or r['done'][i] != 1.0:
                continue
            last = float(v)
            if r['loss'][i] <= target:
                break
        if last is not None:
            out[name] = last
    return out


def run_table(raw, names, arm, eps, value, theta, r_merge=R_MERGE):
    """(O (n, len(OUT_FIELDS)), done (n,), present (n,)) of an arm at the snapshot `value` (a number or a dict name -> number). An absent run or snapshot gives a nan row; a
    snapshot the run did not reach holds the last accepted state (done = 0)."""
    O, done, present = np.full((len(names), len(OUT_FIELDS)), np.nan), np.zeros(len(names)), np.zeros(len(names), bool)
    for j, name in enumerate(names):
        tc = raw['W|%s|truth_c' % name]
        if arm == 'F':
            k = 'P|%s|k%g|' % (name, value)
            if k + 'c' in raw:
                O[j], done[j], present[j] = outcomes(raw[k + 'c'], raw[k + 'w'], tc, theta, r_merge), 1.0, True
            continue
        key = run_key(name, arm, eps)
        r = MC.unpack_result(raw, key)
        if r is None:
            continue
        i = grid_index(raw['G|' + key], value[name] if isinstance(value, dict) else value)
        if i < 0:
            continue
        _, w, c = r['snaps'][i]
        O[j], done[j], present[j] = outcomes(c, w, tc, theta, r_merge), r['done'][i], True
    return O, done, present


def col(O, field):
    return O[:, OUT_FIELDS.index(field)]


def tune_stop(raw, names):
    """The stop and the threshold from arm F on the development windows: e(kappa, theta) = mean |K_hat - K_true|. Among the configurations within one standard error (of the
    best configuration's errors) of the smallest e: the largest kappa (the earliest stop), then the smallest e, then the largest mean F, then the largest theta. Returns
    dict(kappa, theta, err (n_kappa, n_theta), f (n_kappa, n_theta), se, best (kappa, theta))."""
    E, Fm, SD = (np.zeros((len(KAPPAS), len(THETAS))) for _ in range(3))
    for a, k in enumerate(KAPPAS):
        for b, th in enumerate(THETAS):
            O, _, _ = run_table(raw, names, 'F', None, k, th)
            E[a, b], Fm[a, b], SD[a, b] = col(O, 'err').mean(), col(O, 'f1').mean(), col(O, 'err').std(ddof=1)
    a0, b0 = np.unravel_index(int(np.argmin(E)), E.shape)
    se = SD[a0, b0] / np.sqrt(len(names))
    cands = [(-KAPPAS[a], round(E[a, b], 9), -round(Fm[a, b], 9), -THETAS[b], a, b) for a in range(len(KAPPAS)) for b in range(len(THETAS)) if E[a, b] <= E[a0, b0] + se]
    a1, b1 = sorted(cands)[0][4:]
    return dict(kappa=float(KAPPAS[a1]), theta=float(THETAS[b1]), err=E, f=Fm, se=float(se), best=(float(KAPPAS[a0]), float(THETAS[b0])))


def tune_values(raw, names, arm, eps, values, thetas):
    """The (value, theta) of the smallest mean count error of an arm over the development windows, among snapshot values present on every window: ties to the earlier value,
    then the larger theta. Returns (value, theta, mean error, table (n_values, n_thetas), the values used)."""
    vals, rows = [], []
    for v in values:
        row, ok = [], True
        for th in thetas:
            O, _, present = run_table(raw, names, arm, eps, v, th)
            ok = ok and bool(present.all())
            row.append(col(O, 'err').mean() if ok else np.nan)
        if ok:
            vals.append(float(v))
            rows.append(row)
    T = np.array(rows, float).reshape(len(vals), len(thetas))
    if len(vals) == 0:
        return float('nan'), float('nan'), float('nan'), T, vals
    cands = sorted((round(T[a, b], 9), vals[a], -thetas[b], a, b) for a in range(len(vals)) for b in range(len(thetas)))
    a, b = cands[0][3:]
    return vals[a], float(thetas[b]), float(T[a, b]), T, vals


def step_values():
    """The step snapshots of the discrete optimizers within the training budget."""
    return [float(s) for s in STEPS if s <= S_MAX]


def flow_values(arm, kappa):
    """The snapshot values of a flow arm that the test runs will hold (the time stops >= the frozen kappa; for arm C also its early snapshots)."""
    v = [tau_of(k) for k in DEV_FLOW_KAPPAS if k >= kappa]
    if arm == 'C':
        v = [f * tau_of(min(DEV_FLOW_KAPPAS)) for f in C_EARLY] + v
    return sorted(v)


def bic_values(raw, names, arm, eps, theta, values):
    """The BIC stop of an arm (label-free): per window the snapshot, among `values` present and reached, of the smallest 2 F + BIC_PAR K ln N, K the merged formed units
    (anywhere) at theta and N the pixels. Returns dict name -> value (absent if the run has none)."""
    out = {}
    for name in names:
        key = run_key(name, arm, eps)
        r = MC.unpack_result(raw, key)
        if r is None:
            continue
        best = None
        for v in values:
            i = grid_index(raw['G|' + key], v)
            if i < 0 or r['done'][i] != 1.0 or not np.isfinite(r['loss'][i]):
                continue
            _, w, c = r['snaps'][i]
            b = 2.0 * r['loss'][i] + BIC_PAR * len(formed_units(c, w, theta)[1]) * np.log(NW * NW)
            if best is None or b < best[0]:
                best = (b, float(v))
        if best is not None:
            out[name] = best[1]
    return out


def freeze_record(raw):
    """Everything the test needs, from the development raw: kappa and theta (arm F on the dev and devx windows), the step budget of D and E, the best case of every arm
    (its own best snapshot and threshold on the dev windows), the constants, the hashes."""
    dev, devx = stage_names(raw, 'dev'), stage_names(raw, 'devx')
    ts = tune_stop(raw, dev + devx)
    kappa, theta = ts['kappa'], ts['theta']
    steps, best = {'D': float(S_MAX), 'E': float(S_MAX)}, {}
    for arm in MC.ARMS:
        v, th, e, _, _ = tune_values(raw, dev, arm, EPS_MAIN, flow_values(arm, kappa) if arm in MC.ARM_READOUT else step_values(), THETAS)
        best[arm] = [v, th, e]
    return dict(utc=utc(), kappa=kappa, theta=theta, eps=EPS_MAIN, steps=steps, best=best, r_merge=R_MERGE, r_match=R_MATCH, z_lim=Z_LIM, margin_e=M_E, share_path=SHARE_PATH,
                ratio_c=RATIO_C, gate=GATE, budget_s=BUDGET_TEST, n_dev=len(dev), n_devx=len(devx), sha256={f: sha256_file(os.path.join(REPO, f)) for f in FILES})


def write_freeze(raw_path=None, path=None):
    raw_path, path = raw_path or RAW, path or FREEZE
    raw = load_raw(raw_path)
    if stage_names(raw, 'test'):
        raise RuntimeError('test windows exist already: the freeze is closed')
    rec = freeze_record(raw)
    with open(path, 'w') as f:
        json.dump(rec, f, indent=1, sort_keys=True)
        f.write('\n')
    save_part(raw_path, {'freeze|json': np.array([json.dumps(rec, sort_keys=True)])})
    print(json.dumps({k: v for k, v in rec.items() if k != 'sha256'}, indent=1, sort_keys=True))


# ---------------------------------------------------------------------------------------------------------------------------------------------------------------
# the predictions and the report (from the raw alone)
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------
ARM_LABEL = {'A': 'A: paper unit, gradient flow', 'B': 'B: depth 2, gradient flow', 'C': 'C: depth 1, gradient flow', 'D': 'D: paper unit, Adam',
             'E': 'E: paper unit, Adam with one shared second moment', 'F': 'F: computed sliding path (reference)'}


def fnum(x, nd=2):
    return 'n/a' if not np.isfinite(x) else '%.*f' % (nd, x)


def fiv(m, lo, hi, nd=2):
    return '%s [%s, %s]' % (fnum(m, nd), fnum(lo, nd), fnum(hi, nd))


def table(header, rows):
    L = ['| ' + ' | '.join(header) + ' |', '|' + '---|' * len(header)]
    L += ['| ' + ' | '.join(str(x) for x in r) + ' |' for r in rows]
    return L


def primary_tables(raw, names, fz, theta=None, r_merge=R_MERGE, values=None):
    """{arm: (O, done, present)} of the six arms at the frozen primary snapshots (or at `values`: arm -> value or dict)."""
    theta = fz['theta'] if theta is None else theta
    out = {}
    for arm in MC.ARMS + ('F',):
        v = (values or {}).get(arm)
        if v is None:
            v = fz['kappa'] if arm == 'F' else arm_value(arm, fz['kappa'], fz['steps'])
        out[arm] = run_table(raw, names, arm, fz['eps'], v, theta[arm] if isinstance(theta, dict) else theta, r_merge)
    return out


def paired(T, a, b, field, groups=None):
    """(mean, lo, hi, n) of field(a) - field(b) over the windows where both arms are present; the cluster bootstrap over `groups` if given."""
    ok = T[a][2] & T[b][2]
    d = col(T[a][0], field)[ok] - col(T[b][0], field)[ok]
    m, lo, hi = boot_mean(d) if groups is None else boot_mean_cluster(d, np.asarray(groups)[ok])
    return m, lo, hi, int(ok.sum())


def verdicts(T, agree_ok, n_win, fz, agree_signed=None):
    """The registered predictions from the primary tables T and the path-agreement flags of arm A (one per test window; a missing or unconverged run is a failure).
    `agree_signed`: the same flags against the signed path (the second registered reading of P-path). Returns dict of name -> dict(held, numbers)."""
    V = {}
    for name, x in (('R1', 'B'), ('R2', 'C')):
        m, lo, hi, n = paired(T, x, 'A', 'err')
        V[name] = dict(held=bool(lo > 0.0), m=m, lo=lo, hi=hi, n=n)
    m, lo, hi, n = paired(T, 'C', 'A', 'k_hat')
    ok = T['C'][2] & T['A'][2]
    ka, kc = col(T['A'][0], 'k_hat')[ok].mean() if ok.any() else np.nan, col(T['C'][0], 'k_hat')[ok].mean() if ok.any() else np.nan
    V['R3'] = dict(held=bool(lo > 0.0 and kc >= fz['ratio_c'] * ka), m=m, lo=lo, hi=hi, n=n, ka=float(ka), kc=float(kc))
    m, lo, hi, n = paired(T, 'C', 'A', 'n_awake')
    wa, wc = col(T['A'][0], 'n_awake')[ok].mean() if ok.any() else np.nan, col(T['C'][0], 'n_awake')[ok].mean() if ok.any() else np.nan
    V['R3b'] = dict(held=bool(lo > 0.0 and wc >= fz['ratio_c'] * wa), m=m, lo=lo, hi=hi, n=n, ka=float(wa), kc=float(wc))
    V['R'] = dict(held=bool(V['R1']['held'] and V['R2']['held'] and V['R3']['held']))
    mk, lok, hik, n = paired(T, 'D', 'A', 'k_hat')
    me, loe, hie, _ = paired(T, 'D', 'A', 'err')
    V['O1'] = dict(held=bool(lok > 0.0 and loe > 0.0), mk=mk, lok=lok, hik=hik, me=me, loe=loe, hie=hie, n=n)
    m, lo, hi, n = paired(T, 'E', 'A', 'err')
    V['O2'] = dict(held=bool(lo >= -fz['margin_e'] and hi <= fz['margin_e']), m=m, lo=lo, hi=hi, n=n)
    m, lo, hi, n = paired(T, 'E', 'A', 'k_hat')
    V['O2b'] = dict(held=bool(hi <= fz['margin_e']), m=m, lo=lo, hi=hi, n=n)
    V['O'] = dict(held=bool(V['O1']['held'] and V['O2']['held']))
    for name, flags in (('path', agree_ok), ('pathS', agree_signed if agree_signed is not None else np.zeros(0))):
        k = int(np.sum(flags))
        lo, hi = clopper_pearson(k, n_win)
        V[name] = dict(held=bool(n_win > 0 and k / n_win >= fz['share_path']), k=k, n=n_win, lo=lo, hi=hi)
    okf = T['F'][2]
    ef, kt = col(T['F'][0], 'err')[okf].mean() if okf.any() else np.nan, col(T['F'][0], 'k_true')[okf].mean() if okf.any() else np.nan
    V['gate'] = dict(held=bool(ef <= fz['gate'] * kt), ef=float(ef), kt=float(kt))
    return V


def agree_flags(raw, names, eps, kappa, which='nn', field='set_ok'):
    """Per window: (flag, present, done): the path-agreement field of arm A at tau_mid of the reference `kappa` against arm F's block ('nn') or the signed path's ('sg')."""
    flag, present, done = np.zeros(len(names)), np.zeros(len(names), bool), np.zeros(len(names))
    for j, name in enumerate(names):
        k = 'Q|%s|%s%g' % (run_key(name, 'A', eps), which, kappa)
        if k in raw:
            q = raw[k]
            present[j], done[j] = True, q[-1]
            v = q[AGREE_FIELDS.index(field)]
            flag[j] = float(np.isfinite(v) and v == 1.0 and q[-1] == 1.0)
    return flag, present, done


def held(b):
    return 'HELD' if b else 'FAILED'


def stage_summary(raw, names, fz, title, groups=None):
    """The outcome table of a stage at the primary snapshots: one row per arm."""
    T = primary_tables(raw, names, fz)
    rows = []
    for arm in ('A', 'B', 'C', 'D', 'E', 'F'):
        O, done, present = T[arm]
        ok = present
        if not ok.any():
            rows.append((ARM_LABEL[arm], '0 of %d' % len(names), '-', '-', '-', '-', '-', '-', '-', '-'))
            continue
        rows.append((ARM_LABEL[arm], '%d of %d' % (int(present.sum()), len(names)), '%d' % int((done[ok] == 1.0).sum()), fnum(col(O, 'k_true')[ok].mean()), fiv(*boot_mean(col(O, 'k_hat')[ok])),
                     fiv(*boot_mean(col(O, 'err')[ok])), fiv(*boot_mean(col(O, 'signed')[ok])), fiv(*boot_mean(col(O, 'f1')[ok]), nd=3), fnum(col(O, 'n_neg')[ok].mean()), fnum(col(O, 'n_awake')[ok].mean(), 1)))
    return T, table(('arm', 'runs present', 'reached the snapshot', 'mean K_true', 'mean K_hat', 'mean count error', 'mean signed error', 'mean F', 'mean negative units counted',
                     'mean units with |w| > %g (of %d)' % (W_AWAKE, (NW + 2) ** 2)), rows)


def report(raw):
    L = ['# The mechanism on real star-field windows: the factorised network trained on SDSS pixels, scored against HST, with the two conditions switched off', '',
         'Generated by `experiments/exp_mechreal.py --from-raw` from `mechreal_raw.npz`; every number below is printed by the script. Plan: `mechreal_plan.md` (committed before any test '
         'window was processed by any arm); freeze: `mechreal_freeze.json`. Windows of %d x %d px, lattice spacing 1 px (%d units), high-pass %g px, data in units of the flux limit '
         '(matched-filter significance %g); truth = HST (HUGS) sources above the limit in the window interior (%.1f px border excluded), stars closer than %.1f px are one source. '
         'K_hat: formed units (|w| > theta, either sign) after merging same-sign units within %.1f px, in the interior. Count error = |K_hat - K_true|. F at a match radius of %.1f px. '
         'Intervals: 95 %% percentile bootstrap over windows (%d resamples, seed %d).' % (NW, NW, (NW + 2) ** 2, HP, Z_LIM, EDGE, R_BLEND, R_MERGE, R_MATCH, N_BOOT, BOOT_SEED), '']
    if 'freeze|json' not in raw:
        return '\n'.join(L + ['No freeze record in the raw: development is not finished.', '']) + '\n'
    fz = json.loads(str(raw['freeze|json'][0]))
    dev, devx, test = stage_names(raw, 'dev'), stage_names(raw, 'devx'), stage_names(raw, 'test')
    L += ['Frozen (UTC %s) on the development windows only: time stop kappa = %g (tau = 1 / (kappa x %g) in the clock of each flow), formed threshold theta = %g, eps = %s; D and E are scored after their '
          'training budget of %d steps (step size %g).' % (fz['utc'], fz['kappa'], Z_LIM, fz['theta'], eps_tag(fz['eps']), fz['steps']['D'], LR_ADAM), '']

    # ---- 1. test: verdicts
    if test:
        groups = [raw['W|%s|label' % n][0] + '_%d' % int(raw['W|%s|meta' % n][0]) for n in test]
        T, tab = stage_summary(raw, test, fz, 'test')
        flag, present, done = agree_flags(raw, test, fz['eps'], fz['kappa'])
        V = verdicts(T, flag, len(test), fz, agree_flags(raw, test, fz['eps'], fz['kappa'], 'sg')[0])
        L += ['## 1. Test windows (n = %d; selected %s UTC): the registered predictions' % (len(test), str(raw['select|test|utc'][0])), '']
        L += ['* Informative? The reference path F has a mean count error of %s against a mean K_true of %s; gate (at most %g x mean K_true): %s.' % (
            fnum(V['gate']['ef']), fnum(V['gate']['kt']), fz['gate'], 'passed' if V['gate']['held'] else 'NOT passed: every arm fails on these windows and the predictions below carry no information'),
              '* **P-R (readout depth): %s.**' % held(V['R']['held']),
              '    * R1, count error of B (depth 2) minus A: %s (n = %d); interval above 0: %s.' % (fiv(V['R1']['m'], V['R1']['lo'], V['R1']['hi']), V['R1']['n'], held(V['R1']['held'])),
              '    * R2, count error of C (depth 1) minus A: %s (n = %d); interval above 0: %s.' % (fiv(V['R2']['m'], V['R2']['lo'], V['R2']['hi']), V['R2']['n'], held(V['R2']['held'])),
              '    * R3, K_hat of C minus A: %s; mean K_hat %s (C) against %s (A), required ratio %g: %s.' % (fiv(V['R3']['m'], V['R3']['lo'], V['R3']['hi']), fnum(V['R3']['kc']), fnum(V['R3']['ka']),
                                                                                                           fz['ratio_c'], held(V['R3']['held'])),
              '    * R3b (registered companion, not part of the verdict: no dormancy), units with |w| > %g of C minus A: %s; mean %s (C) against %s (A), required ratio %g: %s.' % (
                  W_AWAKE, fiv(V['R3b']['m'], V['R3b']['lo'], V['R3b']['hi'], 1), fnum(V['R3b']['kc'], 1), fnum(V['R3b']['ka'], 1), fz['ratio_c'], held(V['R3b']['held'])),
              '* **P-O (optimizer): %s.**' % held(V['O']['held']),
              '    * O1, Adam (D) against A: K_hat difference %s, count-error difference %s (n = %d); both intervals above 0: %s.' % (
                  fiv(V['O1']['mk'], V['O1']['lok'], V['O1']['hik']), fiv(V['O1']['me'], V['O1']['loe'], V['O1']['hie']), V['O1']['n'], held(V['O1']['held'])),
              '    * O2, shared second moment (E) against A: count-error difference %s (n = %d); interval inside [-%g, %g]: %s.' % (
                  fiv(V['O2']['m'], V['O2']['lo'], V['O2']['hi']), V['O2']['n'], fz['margin_e'], fz['margin_e'], held(V['O2']['held'])),
              '    * O2b (registered companion, one-sided: E forms no excess units), K_hat of E minus A: %s; interval below +%g: %s.' % (
                  fiv(V['O2b']['m'], V['O2b']['lo'], V['O2b']['hi']), fz['margin_e'], held(V['O2b']['held'])),
              '* **P-path (arm A forms the path\'s block): %s.** %d of %d windows (share %s, exact 95 %% interval [%s, %s]; required share %g). A window counts when the run reached the snapshot '
              'tau_mid and the formed set is the block of arm F (same count, every formed unit within 0.5 px of a block unit, one-to-one within 0.25 px; the harness\'s C1 and C2).' % (
                  held(V['path']['held']), V['path']['k'], V['path']['n'], fnum(V['path']['k'] / max(V['path']['n'], 1), 3), fnum(V['path']['lo'], 3), fnum(V['path']['hi'], 3), fz['share_path']),
              '* **P-path-S (second registered reading: against the SIGNED path, Definition 1, of which arm F is the non-negative restriction): %s.** %d of %d windows (share %s, exact 95 %% '
              'interval [%s, %s]).' % (held(V['pathS']['held']), V['pathS']['k'], V['pathS']['n'], fnum(V['pathS']['k'] / max(V['pathS']['n'], 1), 3), fnum(V['pathS']['lo'], 3),
                                       fnum(V['pathS']['hi'], 3)), '']
        L += ['### 1.1 Outcomes by arm at the frozen stop', ''] + tab + ['']
        # the path agreement in detail
        rows = []
        for which, lab in (('nn', 'arm F (non-negative path)'), ('sg', 'signed path (Definition 1)')):
            for field, flab in (('set_ok', 'formed set (C1 and C2)'), ('all_ok', 'C1 to C4 (the harness\'s criterion)'), ('c1', 'C1 alone')):
                f, pr, dn = agree_flags(raw, test, fz['eps'], fz['kappa'], which, field)
                lo, hi = clopper_pearson(int(f.sum()), len(test))
                rows.append((lab, flab, '%d of %d' % (int(f.sum()), len(test)), '[%s, %s]' % (fnum(lo, 3), fnum(hi, 3)), '%d' % int(pr.sum()), '%d' % int((dn == 1.0).sum())))
        L += ['### 1.2 Arm A against the path at tau_mid (eps = %s)' % eps_tag(fz['eps']), ''] + table(('reference', 'criterion', 'windows that agree', 'exact 95 % interval', 'runs present', 'reached tau_mid'), rows) + ['']
        cnt = {}
        for j, name in enumerate(test):
            k = 'Q|%s|nn%g' % (run_key(name, 'A', fz['eps']), fz['kappa'])
            if k not in raw:
                why = 'no run'
            else:
                q = dict(zip(AGREE_FIELDS + ('done',), raw[k]))
                kref = raw['P|%s|ref%g|ref' % (name, fz['kappa'])][0]
                why = ('not converged' if q['done'] != 1.0 else 'agrees' if q['set_ok'] == 1.0 else 'empty block' if not np.isfinite(q['n_formed']) else
                       'more formed units than the block, a negative one among them' if q['n_formed'] > kref and q['n_neg'] > 0 else 'more formed units than the block' if q['n_formed'] > kref else
                       'fewer formed units than the block' if q['n_formed'] < kref else 'same count, positions differ')
            cnt[why] = cnt.get(why, 0) + 1
        L += ['Windows by outcome against arm F: ' + '; '.join('%s %d' % (k, cnt[k]) for k in sorted(cnt)) + '.', '']
        # convergence
        rows = []
        for arm in MC.ARMS:
            st, secs, steps = [], [], []
            for name in test:
                k = 'R|%s|run' % run_key(name, arm, fz['eps'])
                if k in raw:
                    st.append(raw[k][0])
                    secs.append(raw[k][1])
                    steps.append(raw[k][5])
            st = np.array(st)
            rows.append((ARM_LABEL[arm], '%d' % len(st), '%d' % int((st == 0).sum()), '%d' % int((st == 1).sum()), '%d' % int((st >= 2).sum()), fnum(np.median(secs) if secs else np.nan, 0),
                         fnum(np.max(secs) if secs else np.nan, 0), fnum(np.median(steps) if steps else np.nan, 0)))
        L += ['### 1.3 Runs', '', 'A run that did not reach a snapshot within its wall-clock budget (%g s) is scored at the last accepted state and is a failure for P-path.' % fz['budget_s'], ''] + table(
            ('arm', 'runs', 'finished', 'timeout', 'failed', 'median seconds', 'max seconds', 'median steps'), rows) + ['']
        # ---- secondary analyses
        L += ['## 2. Secondary analyses (registered in the plan; none decides a verdict)', '']
        rows = []

        def sec_row(label, Tx):
            def d(a, b, f):
                m, lo, hi, n = paired(Tx, a, b, f)
                return fiv(m, lo, hi)
            rows.append((label, fnum(col(Tx['A'][0], 'err')[Tx['A'][2]].mean() if Tx['A'][2].any() else np.nan), d('B', 'A', 'err'), d('C', 'A', 'err'), d('D', 'A', 'err'), d('E', 'A', 'err'), d('C', 'A', 'k_hat'),
                         d('D', 'A', 'k_hat'), d('E', 'A', 'k_hat'), d('A', 'F', 'err')))
        sec_row('primary', T)
        okall = np.ones(len(test), bool)
        for arm in MC.ARMS:
            okall &= T[arm][2] & (T[arm][1] == 1.0)
        if okall.any():
            sec_row('windows where every arm reached its snapshot (n = %d)' % int(okall.sum()), {a: (T[a][0][okall], T[a][1][okall], T[a][2][okall]) for a in T})
        best_theta = {a: fz['best'][a][1] for a in MC.ARMS}
        best_theta['F'] = fz['theta']
        Tb = primary_tables(raw, test, fz, theta=best_theta, values={a: fz['best'][a][0] for a in ('B', 'C', 'D', 'E')})
        Tb['A'] = T['A']
        sec_row('best case of B, C, D, E (each at its own best development snapshot and threshold; A as frozen)', Tb)
        Te = primary_tables(raw, test, fz, values={a: equal_fit_values(raw, test, a, fz['eps'], fz['kappa'], step_values()) for a in ('D', 'E')})
        sec_row('equal fit for D and E (first snapshot whose training loss is at most that of arm F\'s block; label-free)', Te)
        Tc = primary_tables(raw, test, fz, values={a: bic_values(raw, test, a, fz['eps'], fz['theta'], flow_values(a, fz['kappa']) if a in MC.ARM_READOUT else step_values()) for a in MC.ARMS})
        sec_row('BIC stop over each arm\'s snapshots (label-free, the same rule for A to E)', Tc)
        for r in R_MERGE_SENS:
            sec_row('merge radius %g px' % r, primary_tables(raw, test, fz, r_merge=r))
        L += table(('analysis', 'mean count error of A', 'B - A, count error', 'C - A, count error', 'D - A, count error', 'E - A, count error', 'C - A, K_hat', 'D - A, K_hat', 'E - A, K_hat',
                    'A - F, count error'), rows) + ['']
        rows = []
        for a, b, f in (('B', 'A', 'err'), ('C', 'A', 'err'), ('D', 'A', 'err'), ('E', 'A', 'err'), ('C', 'A', 'k_hat'), ('D', 'A', 'k_hat'), ('E', 'A', 'k_hat')):
            m, lo, hi, n = paired(T, a, b, f, groups)
            rows.append(('%s - %s, %s' % (a, b, 'count error' if f == 'err' else 'K_hat'), fiv(m, lo, hi), '%d' % len(set(groups))))
        L += ['Cluster bootstrap (whole patches resampled) of the primary differences:', ''] + table(('difference', 'mean [95 % interval]', 'patches'), rows) + ['']
        rows = []
        for arm in ('A', 'B', 'C', 'D', 'E', 'F'):
            O, _, pr = T[arm]
            rows.append((ARM_LABEL[arm], fiv(*boot_mean(col(O, 'k_pos')[pr])), fiv(*boot_mean(np.abs(col(O, 'k_pos')[pr] - col(O, 'k_true')[pr]))), fiv(*boot_mean(col(O, 'n_formed_all')[pr]))))
        L += ['Positive units only (a negative unit is not a star), and the merged formed units anywhere on the lattice:', ''] + table(
            ('arm', 'mean positive K_hat', 'mean count error of the positive units', 'mean formed units on the whole lattice'), rows) + ['']
        labs = sorted(set(str(raw['W|%s|label' % n][0]) for n in test))
        rows = []
        for lab in labs:
            sel = np.array([str(raw['W|%s|label' % n][0]) == lab for n in test])
            rows.append((lab, '%d' % int(sel.sum())) + tuple(fnum(col(T[a][0], 'err')[sel & T[a][2]].mean() if (sel & T[a][2]).any() else np.nan) for a in ('A', 'B', 'C', 'D', 'E', 'F')) + (
                '%d of %d' % (int(flag[sel].sum()), int(sel.sum())),))
        L += ['By cluster (mean count error; windows where arm A forms the path\'s block):', ''] + table(('cluster', 'windows', 'A', 'B', 'C', 'D', 'E', 'F', 'A = block of F'), rows) + ['']

    # ---- development
    L += ['## 3. Development (the two M2 development patches; nothing here is a test)', '']
    if dev:
        ts = tune_stop(raw, dev + devx)
        L += ['### 3.1 The stop and the threshold, from arm F on the %d eligible development candidates (mean count error; mean F in brackets)' % (len(dev) + len(devx)), '']
        L += table(('kappa',) + tuple('theta %g' % t for t in THETAS), [('%g' % k,) + tuple('%s (%s)' % (fnum(ts['err'][a, b]), fnum(ts['f'][a, b], 3)) for b in range(len(THETAS))) for a, k in enumerate(KAPPAS)])
        L += ['', 'Smallest error at kappa %g, theta %g; one standard error %s; rule (plan): among the configurations within one standard error, the largest kappa, then the smallest error, then the '
              'largest F, then the largest theta: kappa = %g, theta = %g.' % (ts['best'][0], ts['best'][1], fnum(ts['se'], 3), ts['kappa'], ts['theta']), '']
        Td, tab = stage_summary(raw, dev, fz, 'dev')
        L += ['### 3.2 The arms on the %d development windows at the frozen stop' % len(dev), ''] + tab + ['']
        rows = []
        for arm in ('D', 'E'):
            for s_ in STEPS:
                O, d_, p_ = run_table(raw, dev, arm, EPS_MAIN, float(s_), fz['theta'])
                if p_.all():
                    lo_ = [float(raw['R|%s|loss' % run_key(n, arm, EPS_MAIN)][grid_index(raw['G|' + run_key(n, arm, EPS_MAIN)], float(s_))]) for n in dev]
                    lf_ = [float(raw['P|%s|k%g|loss' % (n, fz['kappa'])][0]) for n in dev]
                    rows.append((ARM_LABEL[arm], '%d' % s_, fnum(col(O, 'k_hat').mean()), fnum(col(O, 'err').mean()), fnum(col(O, 'n_awake').mean(), 1), fnum(float(np.median(np.array(lo_) / np.array(lf_))), 3)))
        L += ['The discrete optimizers along their training on the development windows (theta = %g; K_true mean %s). They have no clock; the registered primary snapshot is the end of the '
              'training budget, %d steps:' % (fz['theta'], fnum(np.mean([len(raw['W|%s|truth_c' % n]) for n in dev])), S_MAX), ''] + table(
            ('arm', 'steps', 'mean K_hat', 'mean count error', 'mean units with |w| > %g' % W_AWAKE, 'median training loss / loss of the block of F'), rows) + ['']
        rows = [(ARM_LABEL[a], '%.6g' % fz['best'][a][0], '%g' % fz['best'][a][1], fnum(fz['best'][a][2])) for a in MC.ARMS]
        L += ['Best case of every arm on the development windows (its own snapshot, tau or steps, and threshold):', ''] + table(('arm', 'snapshot', 'theta', 'mean count error'), rows) + ['']
        rows = []
        for eps in EPS_LIST:
            for kap in DEV_REF_KAPPAS:
                f, pr, dn = agree_flags(raw, dev, eps, kap)
                f2 = agree_flags(raw, dev, eps, kap, 'sg')[0]
                f3 = agree_flags(raw, dev, eps, kap, 'nn', 'all_ok')[0]
                O, d_, p_ = run_table(raw, dev, 'A', eps, tau_of(kap), fz['theta'])
                if pr.any():
                    rows.append((eps_tag(eps), '%g' % kap, '%d' % int(pr.sum()), '%d' % int((dn == 1.0).sum()), '%d' % int(f.sum()), '%d' % int(f3.sum()), '%d' % int(f2.sum()),
                                 fnum(col(O, 'err')[p_].mean() if p_.any() else np.nan)))
        L += ['Arm A on the development windows by eps and stop (path agreement at tau_mid; count error at the stop):', ''] + table(
            ('eps', 'kappa', 'runs', 'reached tau_mid', 'formed set = block of F', 'C1 to C4 against F', 'formed set = signed block', 'mean count error'), rows) + ['']

    # ---- windows
    for stage, names in (('test', test), ('dev', dev)):
        if not names:
            continue
        T = primary_tables(raw, names, fz)
        flag = agree_flags(raw, names, fz['eps'], fz['kappa'])[0]
        rows = []
        for j, n in enumerate(names):
            m = dict(zip(WIN_META, raw['W|%s|meta' % n]))
            rows.append((n, str(raw['W|%s|label' % n][0]), '(%d, %d)' % (m['px'], m['py']), '(%d, %d)' % (m['wx'], m['wy']), fnum(m['sigma'], 3), fnum(m['f_lim'], 0), fnum(m['zmax'], 1), '%d' % m['k_true']) + tuple(
                ('%d%s' % (col(T[a][0], 'k_hat')[j], '' if T[a][1][j] == 1.0 else '*') if T[a][2][j] else '-') for a in ('A', 'B', 'C', 'D', 'E', 'F')) + ('yes' if flag[j] else 'no',))
        L += ['## %s. The %s windows (K_hat by arm at the frozen stop; * = the run did not reach the snapshot)' % ('4' if stage == 'test' else '5', stage), ''] + table(
            ('window', 'cluster', 'patch', 'origin', 'PSF sigma (px)', 'f_lim (DN)', 'brightest star (f_lim)', 'K_true', 'A', 'B', 'C', 'D', 'E', 'F', 'A = block of F'), rows) + ['']
        if 'select|%s|log' % stage in raw:
            L += ['Selection: ' + '; '.join(str(x) for x in raw['select|%s|log' % stage]) + '.', '']
    return '\n'.join(L) + '\n'


def bundle_script(jobfile, procs=4):
    """A self-contained worker script (text) for a remote CPU machine: it carries `src/core/imaging.py`, `mechreal_core.py` and the job file, runs `mechreal_core.work` and leaves
    result.npz in its working directory. Only small derived arrays travel (the normalised windows, their weights, the kernel parameters)."""
    def b64(path):
        with open(path, 'rb') as f:
            return base64.b64encode(f.read()).decode()
    return ('import base64, os, sys\n'
            'FILES = {"imaging.py": "%s", "mechreal_core.py": "%s", "jobs.npz": "%s"}\n'
            'os.makedirs("wk", exist_ok=True)\n'
            'for n, b in FILES.items():\n'
            '    open(os.path.join("wk", n), "wb").write(base64.b64decode(b))\n'
            'os.environ["OMP_NUM_THREADS"] = "1"\nos.environ["OPENBLAS_NUM_THREADS"] = "1"\nos.environ["MKL_NUM_THREADS"] = "1"\n'
            'sys.path.insert(0, "wk")\n'
            'import mechreal_core as MC\n'
            'if __name__ == "__main__":\n'
            '    MC.work("wk/jobs.npz", "result.npz", %d)\n' % (b64(os.path.join(REPO, 'src', 'core', 'imaging.py')), b64(os.path.join(HERE, 'mechreal_core.py')), b64(jobfile), procs))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--select')
    ap.add_argument('--path')
    ap.add_argument('--pack')
    ap.add_argument('--work')
    ap.add_argument('--bundle')
    ap.add_argument('--collect', nargs=2)
    ap.add_argument('--freeze', action='store_true')
    ap.add_argument('--from-raw', action='store_true')
    ap.add_argument('--out')
    ap.add_argument('--arms', default='A,B,C,D,E')
    ap.add_argument('--eps', default='1e-8')
    ap.add_argument('--kappas')
    ap.add_argument('--ref-kappas')
    ap.add_argument('--budget', type=float)
    ap.add_argument('--chunks', type=int, default=1)
    ap.add_argument('--procs', type=int, default=1)
    a = ap.parse_args()
    fz = load_freeze()
    kappas = [float(k) for k in a.kappas.split(',')] if a.kappas else None
    refk = [float(k) for k in a.ref_kappas.split(',')] if a.ref_kappas else None
    for stage in (a.select, a.path, a.pack, a.collect[0] if a.collect else None):
        if stage == 'test':
            check_freeze(fz)
    if a.select:
        cmd_select(a.select, RAW)
    if a.path:
        dev = a.path != 'test'
        cmd_path(a.path, RAW, kappas or (KAPPAS if dev else (fz['kappa'],)), refk or ((DEV_REF_KAPPAS if a.path == 'dev' else ()) if dev else (fz['kappa'],)))
    if a.pack:
        dev = a.pack == 'dev'
        cmd_pack(a.pack, RAW, a.out, a.arms.split(','), [float(e) for e in a.eps.split(',')], kappas or (DEV_FLOW_KAPPAS if dev else test_flow_kappas(fz)),
                 refk or (DEV_REF_KAPPAS if dev else (fz['kappa'],)), a.budget or (BUDGET_DEV if dev else fz['budget_s']), a.chunks,
                 None if dev else {arm: test_steps(fz, arm) for arm in ('D', 'E')})
    if a.work:
        MC.work(a.work, a.out, a.procs)
    if a.bundle:
        with open(a.out, 'w') as f:
            f.write(bundle_script(a.bundle, a.procs))
    if a.collect:
        cmd_collect(a.collect[0], RAW, a.collect[1].split(','))
    if a.freeze:
        write_freeze()
    if a.from_raw:
        text = report(load_raw(RAW))
        with open(MD, 'w') as f:
            f.write(text)
        print('wrote', MD)


DEV_FLOW_KAPPAS, DEV_REF_KAPPAS = (40.0, 30.0, 20.0), (30.0, 20.0)      # development: the time stops the network arms are integrated to, and those with path-agreement references


def test_steps(fz, arm):
    """The step snapshots of a discrete optimizer in the test runs: the grid up to its step budget."""
    return np.array([s for s in STEPS if s <= float(fz['steps'][arm])], float)


def test_flow_kappas(fz):
    """The time stops of the test runs: the frozen kappa and the earlier ones of the development grid (free snapshots on the way)."""
    return tuple(k for k in DEV_FLOW_KAPPAS if k >= fz['kappa'])


def load_freeze(path=None):
    path = path or FREEZE
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)


def check_freeze(fz):
    """The test stage runs only with a freeze record whose hashed files are unchanged."""
    if fz is None:
        raise RuntimeError('no freeze record: the test stage is closed')
    bad = [f for f, h in fz['sha256'].items() if sha256_file(os.path.join(REPO, f)) != h]
    if bad:
        raise RuntimeError('frozen files changed: %s' % ', '.join(bad))


if __name__ == '__main__':
    main()
