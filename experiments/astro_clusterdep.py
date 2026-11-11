"""
Cluster dependence of the ISS lead: the merge of units the kernel cannot tell apart (the repair tested on held-out clusters) and a noise model with a model-error floor (the first attempt, dropped)
===========================================================================================================================================================================================

A NEW module (nothing frozen is edited; `imaging.py`, `psf_general.py`, `exp_astro_m2.py`, `exp_astro_psf.py` are imported unchanged). It holds

* the repair: `coherence_table` (the normalised autocorrelation of the kernel of a frame, tabulated on offsets of +-4 px), `coherence_at` (bilinear lookup) and `merge_coherent`: after the path has stopped, two units whose
  kernels have coherence at least rho are one source (the pair with the largest coherence first; the merged unit sits at the flux-weighted centroid and carries the summed flux; coherences with the neighbours are recomputed).
  A post-processing of the finished catalogue; label-free (pixels, stamp, catalogue; the truth is never read). The threshold rho* = 0.8 is declared in `results/astro/clusterdep_plan.md`.
* the first attempt, reported as a development result and not carried into the test (`clusterdep.md`, section 2b): `calibrate` / `calibrated_weights` fit  var = s0^2 (E/gain + dark) + eta^2 qd^2  to the residual of a first-pass
  catalogue (qd: the smoothed absolute value of its high-passed model; a robust moment fit in bins of qd; s0 >= 1, eta >= 0) and `run_path_weighted` runs the frozen program with those weights. Where the first-pass residual
  is white at the nominal noise (s0 = 1, eta = 0) the weights are the frozen ones and the path is the frozen path.
"""
import os
import sys
import time

os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
from scipy.ndimage import gaussian_filter

HERE = os.path.abspath(os.path.dirname(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
sys.path.insert(0, REPO)
import exp_astro_m2 as X                  # noqa: E402
from src.core import imaging as IM        # noqa: E402
from src.core import psf_general as PG    # noqa: E402

QD_SIGMA = 1.0                            # px: the smoothing of |model| that sets where the model error lives
N_BINS = 8                                # bins of qd (equal counts) of the moment fit
BORDER = 5                                # px: the interior of the residual image used by the fit (as `exp_astro_psf.residual_sd`)
S0_MIN = 1.0                              # the calibration never lowers the nominal noise


def nominal_weights(patch):
    """(sw, var0, expected): the frozen weights of `exp_astro_m2.run_path` and the variance image they are the inverse square root of."""
    dn = patch['dn']
    expected = np.maximum(gaussian_filter(dn, 1.5, mode='reflect'), np.percentile(dn, 1))
    var0 = expected / patch['gain'] + patch['dark']
    return 1.0 / np.sqrt(var0), var0, expected


def make_kernel(patch, hp, fit):
    """The kernel of the frozen programs: `psf_general.make_psf(fit, ...)` of `fit` = `psf_general.fit_adequate(stamp)`, whose isotropic case is the mixture of `imaging.fit_psf_parts` that
    `exp_astro_m2.run_path` builds and whose general case is the kernel of `exp_astro_psf.run_path_psf`; `fit` None means the isotropic mixture."""
    ny, nx = patch['dn'].shape
    if fit is None:
        fit = PG.fit_adequate(patch['psf'], rms_goal=np.inf)
    return PG.make_psf(fit, ny, nx, hp)


def robust_var(v):
    """(1.4826 MAD)^2 of the values v about their median."""
    v = np.asarray(v, float)
    return float((1.4826 * np.median(np.abs(v - np.median(v)))) ** 2)


def model_error_image(P, det, sigma=QD_SIGMA):
    """qd: the smoothed absolute value (DN) of the high-passed model image of the catalogue det = (n, 3) of x, y, flux."""
    det = np.asarray(det, float).reshape(-1, 3)
    M = P.synth(det[:, :2], det[:, 2])
    return gaussian_filter(np.abs(M), sigma, mode='reflect'), M


def calibrate(patch, P, det, n_bins=N_BINS, border=BORDER):
    """Fit s0^2 and eta^2 of  var = s0^2 var0 + eta^2 qd^2  to the residual of the first-pass catalogue `det` on the high-passed patch: in n_bins bins of qd (equal counts, the interior
    of the image) the robust variance (1.4826 MAD)^2 of the residual r = H y - model (DN) is regressed on (mean var0, mean qd^2) by least squares weighted by the bin counts; s0^2 is
    clipped to >= S0_MIN^2 and eta^2 to >= 0 (a constrained refit if either clips). Returns dict(s0, eta, qd, M, var0, n_bins used). Label-free: pixels, kernel and catalogue only."""
    sw, var0, _ = nominal_weights(patch)
    qd, M = model_error_image(P, det)
    y = P.highpass(patch['dn'] - np.median(patch['dn']))
    r = y - M
    sl = (slice(border, -border), slice(border, -border))
    q, rr, v0 = qd[sl].ravel(), r[sl].ravel(), var0[sl].ravel()
    order = np.argsort(q, kind='stable')
    bins = np.array_split(order, n_bins)
    V = np.array([robust_var(rr[b]) for b in bins])
    A0 = np.array([v0[b].mean() for b in bins])
    A1 = np.array([(q[b] ** 2).mean() for b in bins])
    w = np.array([len(b) for b in bins], float)
    # weighted non-negative least squares of V on (A0, A1), with s0^2 >= S0_MIN^2: V - S0_MIN^2 A0 = d0 A0 + eta2 A1, d0, eta2 >= 0
    T = V - S0_MIN ** 2 * A0
    best = None
    cand = []
    X2 = np.c_[A0, A1] * np.sqrt(w)[:, None]
    yv = T * np.sqrt(w)
    sol = np.linalg.lstsq(X2, yv, rcond=None)[0]
    cand.append(sol)
    cand.append(np.array([0.0, max(float(X2[:, 1] @ yv) / max(float(X2[:, 1] @ X2[:, 1]), 1e-300), 0.0)]))
    cand.append(np.array([max(float(X2[:, 0] @ yv) / max(float(X2[:, 0] @ X2[:, 0]), 1e-300), 0.0), 0.0]))
    cand.append(np.zeros(2))
    for s in cand:
        if s[0] < 0 or s[1] < 0:
            continue
        c = float(np.sum((X2 @ s - yv) ** 2))
        if best is None or c < best[0]:
            best = (c, s)
    d0, eta2 = best[1]
    return dict(s0=float(np.sqrt(S0_MIN ** 2 + d0)), eta=float(np.sqrt(eta2)), qd=qd, M=M, var0=var0, bins=(V, A0, A1, w))


def calibrated_weights(patch, P, det):
    """(sw, params): the weights 1 / sqrt(s0^2 var0 + eta^2 qd^2) of the calibration fitted to det, and its parameters (dict s0, eta)."""
    c = calibrate(patch, P, det)
    sw = 1.0 / np.sqrt(c['s0'] ** 2 * c['var0'] + c['eta'] ** 2 * c['qd'] ** 2)
    return sw, dict(s0=c['s0'], eta=c['eta'])


def run_path_weighted(patch, hp, rex, sw, fit=None, kappas=X.KAPPAS):
    """`exp_astro_m2.run_path` / `exp_astro_psf.run_path_psf` (time rule, dual order, clean-up by the units' own significance, key 0.0 the BIC block) with the weight image `sw` in place of
    the frozen 1 / sqrt(E / gain + dark); the kernel is that of `make_kernel(patch, hp, fit)`. With `sw` the frozen weights this is the frozen program. Returns (kappa -> (x, y, flux), seconds, entrants, info)."""
    ny, nx = patch['dn'].shape
    P = make_kernel(patch, hp, fit)
    L = IM.Lattice(P, h=X.LAT_H, margin=X.MARGIN)
    y = P.highpass(patch['dn'] - np.median(patch['dn']))
    reach = P.reach()
    scale = IM.noise_scale(L, sw ** 2)
    t0 = time.time()
    R = IM.sliding_path(P, y, sw, L, max_entries=X.MAX_ENTRIES, tol=1e-6, passes=1, rtrunc=10.0 / reach, r_ex=rex, fit_radius=X.FIT_RADIUS, order='dual', tau_max=1.0 / (min(kappas) * scale))
    secs = time.time() - t0
    out = {}
    for k in kappas:
        i, K = IM.stop_time(R['hist'], R['taus'], 1.0 / (k * scale))
        h = R['hist'][i]
        keep = IM.unit_significance(P, sw, h['c'], h['w']) >= k
        out[k] = (h['c'][keep, 0].copy(), h['c'][keep, 1].copy(), h['w'][keep].copy())
    i, K = IM.stop_bic(R['hist'], ny * nx, n_par=3)
    h = R['hist'][i]
    kb = float(np.sqrt(3.0 * np.log(ny * nx)))
    keep = IM.unit_significance(P, sw, h['c'], h['w']) >= kb
    out[0.0] = (h['c'][keep, 0].copy(), h['c'][keep, 1].copy(), h['w'][keep].copy())
    return out, secs, len(R['taus']), dict(bic_index=int(i), n_blocks=len(R['hist']))


# ----------------------------------------------------------------------------------------------------------------------------------
# merging of units that the kernel cannot tell apart (coherence of the dictionary)
# ----------------------------------------------------------------------------------------------------------------------------------
COH_GRID, COH_STEP, COH_FRAME = 4.0, 0.25, 41        # the table of the kernel's normalised autocorrelation: offsets in [-4, 4]^2 px, step 0.25 px, on a 41 x 41 px frame


def coherence_table(patch, fit):
    """rho(dx, dy) = <phi(0), phi(d)> / <phi(0), phi(0)> of the unit-sum kernel of `fit` (no high-pass; its pixel-integrated form on a COH_FRAME-px frame), tabulated on the grid of offsets
    [-COH_GRID, COH_GRID]^2 with step COH_STEP. Returns (offsets, table (n, n)) with table[iy, ix] = rho(offsets[ix], offsets[iy])."""
    P = make_kernel(dict(dn=np.zeros((COH_FRAME, COH_FRAME)), psf=patch['psf']), None, fit)
    c0 = (COH_FRAME - 1) / 2.0
    off = np.arange(-COH_GRID, COH_GRID + 1e-9, COH_STEP)
    ph0 = P.phi(np.array([[c0, c0]]))[:, 0]
    tab = np.empty((len(off), len(off)))
    for iy, dy in enumerate(off):
        ph = P.phi(np.c_[c0 + off, np.full(len(off), c0 + dy)])
        tab[iy] = (ph0 @ ph) / (ph0 @ ph0)
    return off, tab


def coherence_at(off, tab, dx, dy):
    """Bilinear interpolation of the table at the offsets (dx, dy) (arrays); zero outside the grid."""
    from scipy.ndimage import map_coordinates
    shape = np.shape(dx)
    fx = np.atleast_1d((np.asarray(dx, float) - off[0]) / (off[1] - off[0]))
    fy = np.atleast_1d((np.asarray(dy, float) - off[0]) / (off[1] - off[0]))
    out = map_coordinates(tab, [fy, fx], order=1, mode='constant', cval=0.0)
    return out.reshape(shape) if shape else float(out[0])


def merge_coherent(det, off, tab, rho):
    """Merges the units of the catalogue det = (n, 3) of x, y, flux whose kernels have coherence at least `rho` (the pair with the largest coherence first, until none is left): the merged unit
    sits at the flux-weighted centroid and carries the summed flux (a unit of non-positive flux keeps the position of the other). Deterministic (ties: smaller indices). Returns the new (m, 3) array,
    in the order of the surviving first members, and the number of merges."""
    from scipy.spatial import cKDTree
    det = np.asarray(det, float).reshape(-1, 3).copy()
    n = len(det)
    alive = np.ones(n, bool)
    reach = COH_GRID - 1e-6

    def pair_r(i, js):
        js = np.asarray(js, int)
        if len(js) == 0:
            return np.zeros(0)
        dx, dy = det[i, 0] - det[js, 0], det[i, 1] - det[js, 1]
        ok = (np.abs(dx) <= reach) & (np.abs(dy) <= reach)
        return np.where(ok, coherence_at(off, tab, np.where(ok, dx, 0.0), np.where(ok, dy, 0.0)), 0.0)

    pairs = {}
    if n > 1:
        tree = cKDTree(det[:, :2])
        for i, js in enumerate(tree.query_ball_point(det[:, :2], r=np.sqrt(2.0) * reach)):
            js = [j for j in js if j > i]
            for j, r in zip(js, pair_r(i, js)):
                if r >= rho:
                    pairs[(i, int(j))] = float(r)
    merges = 0
    while pairs:
        (i, j), r = max(pairs.items(), key=lambda kv: (kv[1], -kv[0][0], -kv[0][1]))
        fi, fj = det[i, 2], det[j, 2]
        if fi > 0 and fj > 0:
            det[i, 0] = (det[i, 0] * fi + det[j, 0] * fj) / (fi + fj)
            det[i, 1] = (det[i, 1] * fi + det[j, 1] * fj) / (fi + fj)
        elif fj > fi:
            det[i, 0], det[i, 1] = det[j, 0], det[j, 1]
        det[i, 2] = fi + fj
        alive[j] = False
        merges += 1
        pairs = {k: v for k, v in pairs.items() if i not in k and j not in k}
        others = np.flatnonzero(alive)
        others = others[others != i]
        near = others[(np.abs(det[others, 0] - det[i, 0]) <= reach) & (np.abs(det[others, 1] - det[i, 1]) <= reach)]
        for k, rr in zip(near, pair_r(i, near)):
            if rr >= rho:
                pairs[(min(i, int(k)), max(i, int(k)))] = float(rr)
    return det[alive], merges
