"""
Does the factorised RBF network, trained by gradient flow on an IMAGING loss, form the units that the sliding path (and the packaged algorithm) returns?
================================================================================================================================================

The question (a referee's: "No factorized network is ever trained on imaging data"). `src/core/imaging.py` computes the sliding inverse scale space
(`sliding_path`) of a camera frame: the small-initialisation limit of the paper's network, plus adaptations. Here the NETWORK ITSELF is trained, by its
gradient flow from a small start, on small simulated frames with the SAME loss the path uses (pixel-integrated Gaussian PSF kernels, optionally with the
border-consistent high-pass `hp`, pixel weights W = 1/variance of the noise law at the true model, fixed: no reweighting), and its units are compared with the
path's.

The network (the paper's Eqs. 1-2; `factorized_rbf.PaperLaw`; constants of `exp_parametrization.LAW`: m = 12, G = 1, gamma = 0.1, a_max = 1, b_i = 1)
    units i = 1..d on the start lattice of the path (`Lattice`, spacing 1 px, margin 1 px; d = 484 on 20 x 20 px, 676 on 24 x 24), kernel phi(c_i) (H phi if the model
    is high-passed), w_i = nu_i^2 a_i / m, a_i = a_max tanh(alpha_i / a_max), F = 1/2 sum_p W_p r_p^2, r = H y - sum_i w_i H phi_i,
    alpha' = -G m dF/dalpha, nu' = -(gamma m / 4) dF/dnu, c' = -G m dF/dc,   alpha_i(0) = 0, nu_i(0)^2 = eps  (equal walls rho_i = 1).
Equivalent form used for the right-hand side: alpha' = G beta a' g, beta' = gamma beta a g, c' = G m w <r grad phi>, g_i = sum_p W_p r_p phi_i(p) (`PixelPSF.comps`
once per evaluation, then the same sums as `PixelPSF.synth` and `PixelPSF.corr`; a test pins them to those). The conserved Q_i = beta_i - v_s cosh(2 alpha_i / a_max),
v_s = gamma a_max^2 / (4G), is monitored (not imposed). The path time is tau = t / Lambda(eps), Lambda = (pi a_max / (4G)) (2 / (v_s eps))^(1/2)
(`rich_transit.transit_scale`): 7025 at eps = 1e-6, 70248 at 1e-8.

Solver. The flow is integrated in the raw coordinates (alpha, nu, c) of ALL d units (4d = 1936 / 2704 unknowns) by Radau (BDF from the last accepted step if Radau
fails), rtol 1e-8, atol 1e-10 sqrt(eps) on (alpha, nu) and 1e-12 on the centres, from t = 0 to the snapshot times, a wall-clock budget of an hour per run. The Newton
matrix is a SPARSE APPROXIMATE Jacobian (it only affects the convergence of the Newton iteration, not the solution): the analytic own 4 x 4 block of every unit with
the drive frozen, and for the active units (beta_i > 1000 eps) the exact columns by forward differences with the kernel factors of the other units reused; a dense
Jacobian of 2704 unknowns costs seconds per factorisation, this one milliseconds. Part A reruns three frames at rtol 1e-11 and prints the difference of the final states.

The limit objects.
  Definition 1 as the brief fixes it: `imaging.sliding_path(P, y, sw, lat, rho=1, local=False, r_ex=0, order='dual')` with the same P, lattice and weights (global refits,
  widths fixed, no reweighting, no clean-up), non-negative amplitudes; its stop block (below) polished (unbounded Levenberg-Marquardt from the path's block, accepted if the
  cost does not rise, the signs hold and no centre moves more than 0.05 px) so that the trust-region fit's own error (1e-7 px) is not what is measured; the centre box
  (2 px) is checked inactive. The network's readout is SIGNED and `sliding_path`'s is not: a unit whose dual reaches the wall -1 forms a negative unit in the network that
  the non-negative path cannot have. So the network is also compared with
  the signed path (`signed_path`, the same algorithm with |u| = rho and a sign constraint per unit; Definition 1 of the paper, of which `sliding_path` is the restriction to
  non-negative amplitudes), whose block at the snapshot times is taken. All emitters are positive. The runner-up margin of the path (the smallest 1 - u_j of a kernel that
  is not in the block, over the window up to tau_mid) is recorded: the limit is not uniform in the margin.

Frames (fixed seeds, declared here and not changed after the first run; geometry seed = 1000 + index (S) or 2000 + index (A), noise seed = geometry seed + 50000)
  SMLM-like (S): 20 x 20 px, PSF sigma 1.2 px, background 100 photons/px, EM-CCD (`imaging.simulate_emccd`; weights from `imaging.emccd_variance` at the true model).
      K = 3..6 emitters, fluxes uniform in 500-1000 (L), 1000-2000 (M), 2000-3000 (H) photons, centres in [3, 17]^2 at least 4 px apart, except a close pair at 2.5 sigma
      (S06, level M) and at 2.0 sigma (S10, level H). S01-S04 L (K = 3, 4, 5, 6), S05-S08 M (3, 4 with the pair, 5, 6), S09-S12 H (3, 5 with the pair, 4, 6).
      Noise-free copies (same emitters and weights, data = the true model): S03n, S06n, S10n, S12n.
  Astronomy-like (A): 24 x 24 px, sigma 1.2 px, sky 300 DN, gain 4.6, dark variance 1 DN^2, DN = Poisson(mean gain)/gain + N(0, dark), weights 1/(mean/gain + dark) at the
      true model (the `exp_astro_m2` law sw = 1/sqrt(expected/gain + dark)). K = 6, 7, 8, 8, 9, 10 sources (A01 ... A06), fluxes log-uniform 400-25000 DN; in A03 and A06 a
      bright star (15000-25000 DN) with a faint satellite (3-6 % of its flux) 2.5 sigma away. Noise-free copies: A03n, A06n.
  Two variants of each: hp = 3 px (the model H phi, the frame H y; the main one) and none (the plain Gaussian dictionary; "the high-pass is a change of dictionary").

Stop rule of the path (fixed before looking at any network). Noisy frames: the block of `imaging.stop_bic` (3 parameters per unit) on the path run until the minimum is
followed by at least three more blocks (or the path ends); noise-free frames: the first block with exactly the emitters' number of units (reported: whether it recovers them
to 1e-3 px). K_stop = its number of units, tau_K = the time of the entry that made it, tau_next = the time of the next entry (infinite if none), tau_end = min(tau_next,
4 tau_K), tau_mid = tau_K + (tau_end - tau_K) / 2 (the midpoint of the brief; the cap 4 tau_K acts on noise-free frames, whose next entry is the numerical-zero one, and on a
gap larger than 3 tau_K), tau_q1 = tau_K + (tau_end - tau_K) / 4 (an earlier snapshot, a diagnostic). The network is integrated to t = Lambda(eps) tau at both snapshots and
compared at both; the criteria are declared for tau_mid.

ACCEPTANCE (declared before the first run). A network run "forms the path's units" at tau_mid iff all of
  C1 count      the formed set F = {i : |w_i| > THR}, THR = min(0.01 max|w_path|, 0.1 min|w_path|) (the 1 % of the brief; the second term as `factorized_rbf.analyse`, so that a
                faint source still counts), has exactly K_stop members and every one has a block unit within 0.5 px;
  C2 position   the one-to-one (Hungarian) match of F to the block has its largest centre distance <= 0.25 px;
  C3 model      the relative rms difference of the model images (all d units of the network against the block) <= 1e-3 and the largest relative difference of the matched
                amplitudes <= 1e-3;
  C4 dormant    the largest |w| among the units NOT in F, relative to max|w|, <= 1e-3 (theory: O(eps^{3/2}) times a power of the margin).
against the non-negative block of Definition 1 (the brief's criterion); the same four against the signed path's block at the same time are reported besides. Also reported:
the signed counts, whether the formed units are the path's lattice units (identity), the order and the times of formation against the path's entry times, the stationarity of
the formed units, the drift of Q, the dormant units against the theory (`FreeLaw.core_tail` at the path's runner-up margin), steps, evaluations and seconds, and for every failure its first
reason. Failures are reported and not repaired: nothing is tuned after the first run; if a bug is found it is fixed and everything is rerun.
ADDED AFTER THE FIRST RUN, none of it a criterion and each labelled as such in the report: (1) eps = 1e-10 for the frames that failed at 1e-8 (hp: S11, A02, A04, A05, A06, A06n;
plain: S02, S04, S06, S10n, A06, A06n); (2) the duplicates diagnostic (post hoc, from the stored states): same-sign formed units within 0.05 px of each other counted as one unit,
for the runs that have more formed units than the block; (3) the stalled-path check (part S, requested after the declared runs): a path of `imaging.sliding_path` can stall in an
enter/exit limit cycle (one lattice point enters and exits again and again at an identical tau until max_entries); the part recomputes the paths of the declared runs and counts
their events; a stalled path would be a defective path, not a path-versus-network result.
The packaged algorithm (`fit_frame`: r_ex = 2 px, local refits, widths sliding in [0.7, 2.3] sigma, weights iterated from the fit, BIC stop with 4 parameters per unit) is
compared with the non-negative block at its own BIC stop and at the same number of units: equal iff the same number of units, largest matched distance <= 0.25 px, largest
relative amplitude difference <= 5 %; then each adaptation is switched off (exclusion zone, local refits, sliding widths, reweighting), alone and cumulatively, and the
configuration where it becomes equal names the cause. The script's replica of `fit_frame` (the switches) is checked against `fit_frame` itself.

Parts and usage (a part writes its keys into the raw frame by frame and skips what is already there; two processes use different NI_OUT directories, merged by --merge):
    python experiments/exp_network_imaging.py --only P                 the paths (Definition 1, signed path, margins) of all frames and both variants (seconds)
    python experiments/exp_network_imaging.py --only H --eps 1e-6      the network, hp variant
    python experiments/exp_network_imaging.py --only Q --eps 1e-6      the network, plain variant
    python experiments/exp_network_imaging.py --only B                 the bridge to `fit_frame` and the ablations of its adaptations (hp frames)
    python experiments/exp_network_imaging.py --only A --eps 1e-8      the solver check (rtol 1e-11 on S02, S10, A03)
    python experiments/exp_network_imaging.py --only S                 the stalled-path check (recomputes the paths of the declared runs and reads their events)
    python experiments/exp_network_imaging.py --merge DIR/network_imaging_raw.npz[,...]    merge other raws into this one
    python experiments/exp_network_imaging.py --from-raw               regenerate results/synthetic/network_imaging.md from the raw alone
Environment: NI_OUT (directory of the raw and the md, default results/synthetic). Options: --frames S01,S02 (default all), --eps 1e-6,1e-8, --force.
"""
import argparse
import collections
import contextlib
import os
import sys
import time

os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import numpy as np
import scipy.sparse as sp
from scipy.integrate import BDF, Radau
from scipy.optimize import least_squares, linear_sum_assignment, nnls
from scipy.sparse.csgraph import connected_components
from scipy.spatial.distance import cdist

HERE = os.path.abspath(os.path.dirname(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
sys.path.insert(0, REPO)
import exp_parametrization as EP                 # noqa: E402  (the constants of the free law)
from src.core import imaging as IM               # noqa: E402
from src.core.rich_transit import transit_scale  # noqa: E402

OUT = os.environ.get('NI_OUT', os.path.join(REPO, 'results', 'synthetic'))
RAW, MD = os.path.join(OUT, 'network_imaging_raw.npz'), os.path.join(OUT, 'network_imaging.md')

# ---- the law and the numerics -------------------------------------------------------------------------------------------------------------------------------
LAW = EP.law()                                   # m = 12, G = 1, gamma = 0.1, a_max = 1 (FreeLaw; v_s = gamma a_max^2 / (4G) = LAW.c)
EPS_LIST = (1e-6, 1e-8)
RTOL, ATOL_P, ATOL_C = 1e-8, 1e-10, 1e-12        # Radau: rtol; atol = ATOL_P sqrt(eps) on (alpha, nu), ATOL_C on the centres
ACT_FACTOR = 1e3                                 # a unit is active (finite-difference Jacobian columns) when beta_i > ACT_FACTOR eps
FD_REL, FD_FLOOR = 1e-7, (1e-3, 1e-3, 1.0, 1.0)  # finite-difference step: FD_REL max(|z|, floor) for (alpha, nu, x, y)
BUDGET_S = 3600.0                                # wall-clock budget of one run (s)
PATH_TOL, PATH_NFEV = 1e-12, 2000                # the path's block fits
EXTRA_NOISY, EXTRA_FREE, MIN_AFTER = 6, 2, 3     # entries beyond K_true the noisy path runs to; the BIC minimum must be followed by MIN_AFTER blocks
TAU_CAP, Q_EARLY, Q_MID = 4.0, 0.25, 0.5         # tau_next is capped at TAU_CAP tau_K; the snapshots at these fractions of (tau_end - tau_K)
POLISH_MAX = 0.05                                # px: a polished block may move no centre farther than this
THETA, THETA_MIN = 0.01, 0.1                     # formed: |w| > min(THETA max|w_path|, THETA_MIN min|w_path|)
D_NEAR, D_MATCH = 0.5, 0.25                      # px: C1 (a path unit within D_NEAR of every formed unit), C2 (matched within D_MATCH)
REL_TOL, DORM_TOL = 1e-3, 1e-3                   # C3, C4
BR_AMP = 0.05                                    # the bridge: equal iff same count, largest distance <= D_MATCH, largest relative amplitude difference <= BR_AMP
R_EX, WIDTHS_REL = 2.0, (0.7, 2.3)               # the packaged algorithm's exclusion zone (px) and the bounds of the sliding widths (x sigma)

# ---- the frames --------------------------------------------------------------------------------------------------------------------------------------------------
SIGMA, HP, LAT_H, LAT_MARGIN = 1.2, 3.0, 1.0, 1.0
SMLM = dict(ny=20, nx=20, bg=100.0)
ASTRO = dict(ny=24, nx=24, sky=300.0, gain=4.6, dark=1.0)
FLUX = dict(L=(500.0, 1000.0), M=(1000.0, 2000.0), H=(2000.0, 3000.0))      # photons per emitter at the SNR levels low, mid, high
MIN_SEP, SAT_SIGMAS = 4.0, 2.5                                              # px between unrelated emitters; the satellite's distance in sigma
S_TABLE = [('L', 3, None), ('L', 4, None), ('L', 5, None), ('L', 6, None), ('M', 3, None), ('M', 4, 2.5), ('M', 5, None), ('M', 6, None),
           ('H', 3, None), ('H', 5, 2.0), ('H', 4, None), ('H', 6, None)]    # S01 ... S12: (level, K, close pair separation in sigma)
A_TABLE = [(6, False), (7, False), (8, True), (9, False), (10, False), (8, True)]      # A01 ... A06: (K, bright star with a faint satellite)
FREE_COPIES = ('S03', 'S06', 'S10', 'S12', 'A03', 'A06')                    # the frames that also run noise-free, as <name>n


def _specs():
    out = {}
    for k, (level, K, pair) in enumerate(S_TABLE, 1):
        out['S%02d' % k] = dict(name='S%02d' % k, kind='smlm', noisy=True, source='S%02d' % k, seed=1000 + k, level=level, K=K, pair=pair, sat=False)
    for k, (K, sat) in enumerate(A_TABLE, 1):
        out['A%02d' % k] = dict(name='A%02d' % k, kind='astro', noisy=True, source='A%02d' % k, seed=2000 + k, level='-', K=K, pair=None, sat=sat)
    for nm in FREE_COPIES:
        out[nm + 'n'] = dict(out[nm], name=nm + 'n', noisy=False)
    order = ['S%02d' % k for k in range(1, 13)] + [n + 'n' for n in FREE_COPIES if n[0] == 'S'] + ['A%02d' % k for k in range(1, 7)] + [n + 'n' for n in FREE_COPIES if n[0] == 'A']
    return {n: out[n] for n in order}


SPECS = _specs()
FRAMES = tuple(SPECS)


def draw_geometry(spec):
    """The emitters of a frame, from its seed: centres (K, 2) = (x, y) in pixel units and fluxes (K,) (photons for SMLM, DN for astronomy). The pair / satellite are
    the first two emitters; the others are drawn at least MIN_SEP from all earlier ones."""
    rng = np.random.default_rng(spec['seed'])
    smlm = spec['kind'] == 'smlm'
    lo, hi = 3.0, (SMLM['nx'] if smlm else ASTRO['nx']) - 3.0
    K = spec['K']
    sep = spec['pair'] * SIGMA if spec['pair'] else (SAT_SIGMAS * SIGMA if spec['sat'] else None)
    pts = []
    if sep is not None:
        for _ in range(100000):
            p = rng.uniform(lo, hi, 2)
            th = rng.uniform(0.0, 2.0 * np.pi)
            q = p + sep * np.array([np.cos(th), np.sin(th)])
            if (q >= lo).all() and (q <= hi).all():
                break
        pts = [p, q]
    for _ in range(100000):
        if len(pts) >= K:
            break
        p = rng.uniform(lo, hi, 2)
        if all(np.linalg.norm(p - q) >= MIN_SEP for q in pts):
            pts.append(p)
    if len(pts) < K:
        raise RuntimeError('draw_geometry: could not place %d emitters' % K)
    c = np.array(pts)
    if smlm:
        f = rng.uniform(FLUX[spec['level']][0], FLUX[spec['level']][1], K)
    else:
        f = 10.0 ** rng.uniform(np.log10(400.0), np.log10(25000.0), K)
        if spec['sat']:
            f[0] = rng.uniform(15000.0, 25000.0)
            f[1] = f[0] * rng.uniform(0.03, 0.06)
    return c, f


def simulate(spec):
    """The frame of a spec: dict(c, f, frame (with its background), y0 (background removed; the true model if noise-free), bg, var (the oracle variance map at the true model),
    var_a, var_b0 (the noise law var = var_a (bg + model) + var_b0 for `fit_frame`), ny, nx)."""
    c, f = draw_geometry(spec)
    noisy = spec['noisy']
    rng = np.random.default_rng(spec['seed'] + 50000)
    if spec['kind'] == 'smlm':
        ny, nx, bg = SMLM['ny'], SMLM['nx'], SMLM['bg']
        P0 = IM.PixelPSF(ny, nx, SIGMA)
        sig = P0.synth(c, f)
        var = IM.emccd_variance(bg + sig)
        if noisy:
            frame = IM.simulate_emccd(P0, c, f, bg, rng)
            y0 = frame - bg
        else:
            frame, y0 = bg + sig, sig
        var_a, var_b0 = 2.0 / IM.EMCCD['qe'], float(IM.emccd_variance(0.0))
    else:
        ny, nx, bg = ASTRO['ny'], ASTRO['nx'], ASTRO['sky']
        P0 = IM.PixelPSF(ny, nx, SIGMA)
        sig = P0.synth(c, f)
        mean = bg + sig
        var = mean / ASTRO['gain'] + ASTRO['dark']
        if noisy:
            frame = rng.poisson(mean * ASTRO['gain']).astype(float) / ASTRO['gain'] + rng.normal(0.0, np.sqrt(ASTRO['dark']), mean.shape)
            y0 = frame - bg
        else:
            frame, y0 = mean, sig
        var_a, var_b0 = 1.0 / ASTRO['gain'], ASTRO['dark']
    return dict(c=c, f=f, frame=frame, y0=y0, bg=bg, var=var, var_a=var_a, var_b0=var_b0, ny=ny, nx=nx)


def build(spec, hp):
    """The imaging problem of a frame: the PSF (high-passed with width hp, or plain if hp is None), the start lattice and its centres c0 (d, 2) (unit u = b nlx + a at
    (xs[a], ys[b])), the data y = H y0, the weights W = sw^2."""
    sim = simulate(spec)
    P = IM.PixelPSF(sim['ny'], sim['nx'], SIGMA, hp=hp)
    lat = IM.Lattice(P, h=LAT_H, margin=LAT_MARGIN)
    X, Y = np.meshgrid(lat.xs, lat.ys)
    sw = 1.0 / np.sqrt(sim['var'])
    return dict(spec=spec, sim=sim, P=P, lat=lat, c0=np.c_[X.ravel(), Y.ravel()], y=P.highpass(sim['y0']), sw=sw, W=sw ** 2, hp=hp)


# ---------------------------------------------------------------------------------------------------------------------------------------------------------------
# the path (Definition 1) and its stop
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------
def formed_threshold(w_path):
    """THR = min(THETA max|w_path|, THETA_MIN min|w_path|): a unit has formed when |w| exceeds it (the 1 % of the brief, and a tenth of the faintest unit of the block)."""
    a = np.abs(np.asarray(w_path, float))
    return float(min(THETA * a.max(), THETA_MIN * a.min()))


def block_stationarity(P, W, y, c, w):
    """How far the block (c (k, 2), w (k)) is from a stationary point of F: (largest |g_i| / sqrt(sum_p W phi_i^2), the largest Newton displacement of a centre in px,
    |g_x| / (|w| sum_p W (d_x phi)^2)-style estimate with the Gaussian curvature w^2 nrm^2 / (2 sigma^2))."""
    c = np.atleast_2d(np.asarray(c, float))
    w = np.asarray(w, float)
    if len(w) == 0:
        return 0.0, 0.0
    Rw = W * (y - P.synth(c, w))
    v, vx, vy = P.corr(c, Rw)
    nrm = np.sqrt(W.ravel() @ P.phi(c) ** 2)
    curv = np.abs(w) * 0.5 * nrm ** 2 / P.sigma ** 2
    return float(np.max(np.abs(v) / nrm)), float(max(np.max(np.abs(vx) / curv), np.max(np.abs(vy) / curv)))


def polish_block(P, y, sw, c, w, max_shift=POLISH_MAX):
    """The block (c (k, 2), w (k)) polished to the stationary point of the weighted fit: unbounded Levenberg-Marquardt with the analytic Jacobian, three rounds. Accepted
    iff the cost does not rise, every amplitude keeps its sign and no centre moves more than max_shift px. Returns (c, w, accepted, largest centre shift, largest relative
    amplitude shift)."""
    k = len(w)
    swv = np.asarray(sw, float) * np.ones((P.ny, P.nx))

    def res(z):
        return (swv * (y - P.synth(z[k:].reshape(2, k).T, z[:k]))).ravel()

    def jac(z):
        return P.jac(z[k:].reshape(2, k).T, z[:k], swv)
    z0 = np.r_[w, np.asarray(c, float).T.ravel()]
    z = z0.copy()
    for _ in range(3):
        z = least_squares(res, z, jac=jac, method='lm', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=2000).x
    c2, w2 = z[k:].reshape(2, k).T.copy(), z[:k].copy()
    dc = float(np.abs(c2 - c).max())
    dw = float(np.abs(w2 - w).max() / np.abs(w).max())
    ok = bool(np.sum(res(z) ** 2) <= np.sum(res(z0) ** 2) * (1.0 + 1e-12) and (w2 * w > 0).all() and dc <= max_shift)
    return (c2, w2, True, dc, dw) if ok else (np.array(c, float), np.array(w, float), False, dc, dw)


def signed_block_fit(psf, y, swv, c, w, s, c_lo, c_hi, max_nfev=PATH_NFEV, tol=PATH_TOL):
    """`imaging.block_fit` with a sign for every unit (s_i = +1 or -1): min 1/2 sum (sw (y - sum_i w_i phi_i))^2 over s_i w_i >= 0 and the centres in [c_lo, c_hi]: the coefficients
    first by non-negative least squares at the centres in place, then all together by the trust-region reflective solver with the analytic Jacobian. Returns (c, w)."""
    k = len(w)
    s = np.asarray(s, float)
    u1, _ = nnls(swv.reshape(-1, 1) * psf.phi(c) * s[None, :], (swv * y).ravel(), maxiter=200 * k)
    lo = np.r_[np.zeros(k), c_lo.T.ravel()]
    hi = np.r_[np.full(k, np.inf), c_hi.T.ravel()]
    z0 = np.r_[u1, np.clip(c, c_lo, c_hi).T.ravel()]

    def res(z):
        return (swv * (y - psf.synth(z[k:].reshape(2, k).T, s * z[:k]))).ravel()

    def jac(z):
        J = psf.jac(z[k:].reshape(2, k).T, s * z[:k], swv)
        J[:, :k] *= s[None, :]
        return J
    sol = least_squares(res, z0, jac=jac, bounds=(lo, hi), method='trf', x_scale='jac', xtol=tol, ftol=tol, gtol=tol, max_nfev=max_nfev)
    return sol.x[k:].reshape(2, k).T.copy(), s * sol.x[:k]


def signed_path(psf, y, sw, lat, rho=1.0, max_entries=40, box=2.0, tie_tol=1e-9, vtol=1e-12, tau_stop=None):
    """The SIGNED sliding inverse scale space (Definition 1 of the paper; `imaging.sliding_path` is its restriction to non-negative amplitudes): every lattice kernel carries a dual
    u_j with du_j/dtau = <W r phi_j>, it enters when |u_j| reaches the wall rho (with the sign of u_j), the block relaxes to a stationary point of the fit under s_i w_i >= 0 with the
    centres sliding (`signed_block_fit`), a unit at 0 whose correlation points back leaves; global refits, no exclusion zone, widths fixed, weights fixed. Stops when nothing
    grows, after max_entries entries, or before the first entry later than tau_stop. Returns dict(events = [(tau, (b, a), 'enter' / 'exit', sign)], hist = [dict(tau, c, w, s, ids,
    chi2, n_ent)] (hist[0] the empty block))."""
    ny, nx = psf.ny, psf.nx
    swv = np.asarray(sw, float) * np.ones((ny, nx))
    U = np.zeros(lat.shape)
    taken = np.zeros(lat.shape, bool)
    ids, c, w, c0, sg = [], np.zeros((0, 2)), np.zeros(0), np.zeros((0, 2)), np.zeros(0)
    tau, n_ent, events, vscale = 0.0, 0, [], None
    model = np.zeros((ny, nx))

    def chi2():
        return float(np.mean((swv * (y - model)) ** 2))
    hist = [dict(tau=0.0, c=c.copy(), w=w.copy(), s=sg.copy(), ids=[], chi2=chi2(), n_ent=0)]
    for _ in range(11 * max_entries):
        V = lat.corr(swv ** 2 * (y - model))
        if vscale is None:
            vscale = float(np.abs(V).max())
            if vscale <= 0.0:
                break
        v = np.where(taken, 0.0, V)
        up, dn = v > vtol * vscale, v < -vtol * vscale
        dt = np.full(lat.shape, np.inf)
        dt[up] = (rho - U[up]) / v[up]
        dt[dn] = (-rho - U[dn]) / v[dn]
        if not np.isfinite(dt).any():
            break
        dmin = float(dt.min())
        if tau_stop is not None and tau + dmin > tau_stop:
            break
        ent = np.argwhere(dt <= dmin * (1.0 + tie_tol) + tie_tol * 1e-300)
        U[~taken] += dmin * v[~taken]
        tau += dmin
        for b, a in ent:
            s_ = 1.0 if v[b, a] > 0 else -1.0
            U[b, a] = s_ * rho
            taken[b, a] = True
            ids.append((int(b), int(a)))
            p0 = lat.pos(b, a)
            c0, c, w, sg = np.vstack([c0, p0[None, :]]), np.vstack([c, p0[None, :]]), np.r_[w, 0.0], np.r_[sg, s_]
            events.append((tau, (int(b), int(a)), 'enter', s_))
            n_ent += 1
        c, w = signed_block_fit(psf, y, swv, c, w, sg, c0 - box, c0 + box)
        model = psf.synth(c, w)
        g = psf.corr(c, swv ** 2 * (y - model))[0]
        gone = (np.abs(w) <= 1e-12 * max(np.abs(w).max(), 1e-300)) & (sg * g <= vtol * vscale)
        for i in np.flatnonzero(gone)[::-1]:
            events.append((tau, ids[i], 'exit', sg[i]))
            taken[ids[i]] = False
            U[ids[i]] = sg[i] * rho
            del ids[i]
            c, w, c0, sg = np.delete(c, i, 0), np.delete(w, i), np.delete(c0, i, 0), np.delete(sg, i)
        if gone.any():
            model = psf.synth(c, w)
        hist.append(dict(tau=tau, c=c.copy(), w=w.copy(), s=sg.copy(), ids=list(ids), chi2=chi2(), n_ent=n_ent))
        if n_ent >= max_entries:
            break
    return dict(events=events, hist=hist)


SPV_FIELDS = ('thr', 'K', 'same_nn', 'tau_neg', 'id_neg', 'n_neg', 'polished', 'polish_dc', 'polish_dw', 'kkt_amp', 'kkt_pos', 'tau_snap')
SREF_KEYS = ('c', 'w', 'ids', 'tau', 'pv')


def signed_refs(prob, path):
    """The signed path's block at the snapshot times of the non-negative path (the last block whose entry time is <= tau_q1, tau_mid), polished: {'sq1', 'smid'} -> dict(c, w, ids
    (flat lattice index), tau (entry times), pv = the SPV_FIELDS vector: thr, K, whether its units are the non-negative block's, the time of the first negative entry (inf if none),
    the number of negative units, ...)."""
    P, y, sw, lat = prob['P'], prob['y'], prob['sw'], prob['lat']
    pv = dict(zip(PV_FIELDS, path['pv']))
    nlx = lat.shape[1]
    R = signed_path(P, y, sw, lat, max_entries=int(max(pv['K_true'], pv['K_stop'])) + 12, tau_stop=pv['tau_mid'] * (1.0 + 1e-9))
    ent = {tuple(ij): t for t, ij, kind, s_ in R['events'] if kind == 'enter'}
    negs = [(t, ij) for t, ij, kind, s_ in R['events'] if kind == 'enter' and s_ < 0]
    tau_neg, id_neg = (negs[0][0], float(negs[0][1][0] * nlx + negs[0][1][1])) if negs else (float('inf'), -1.0)
    out = {}
    for name in ('q1', 'mid'):
        ts = pv['tau_' + name]
        h = R['hist'][max(j for j, q in enumerate(R['hist']) if q['tau'] <= ts * (1.0 + 1e-12))]
        K = len(h['w'])
        c, w, acc, dc, dw = polish_block(P, y, sw, h['c'], h['w']) if K else (h['c'], h['w'], False, 0.0, 0.0)
        ids = np.array([b * nlx + a for b, a in h['ids']], float)
        ka, kp = block_stationarity(P, prob['W'], y, c, w)
        spv = dict(thr=formed_threshold(w) if K else float('nan'), K=K, same_nn=float(set(ids.astype(int).tolist()) == set(np.asarray(path['blk_ids']).astype(int).tolist())),
                   tau_neg=tau_neg, id_neg=id_neg, n_neg=int((w < 0).sum()), polished=float(acc), polish_dc=dc, polish_dw=dw, kkt_amp=ka, kkt_pos=kp, tau_snap=h['tau'])
        out['s' + name] = dict(c=c, w=w, ids=ids, tau=np.array([ent[tuple(ij)] for ij in h['ids']]), pv=np.array([spv[f] for f in SPV_FIELDS], float))
    return out


PV_FIELDS = ('K_true', 'K_stop', 'i_stop', 'n_hist', 'n_entries', 'interior', 'tau_K', 'tau_next', 'tau_end', 'tau_q1', 'tau_mid', 'cap_active', 'min_gap', 'gap_next',
             'box_active', 'polished', 'polish_dc', 'polish_dw', 'kkt_amp', 'kkt_pos', 'chi2', 'exact_ok', 'n_exits', 'thr', 'secs', 'margin_min', 'margin_mid', 'margin_id')


def snapshot_times(tau_K, tau_next):
    """(tau_end, tau_q1, tau_mid, capped): tau_end = min(tau_next, TAU_CAP tau_K), the snapshots at Q_EARLY and Q_MID of the way from tau_K to tau_end."""
    tau_end = min(float(tau_next), TAU_CAP * float(tau_K))
    return tau_end, tau_K + Q_EARLY * (tau_end - tau_K), tau_K + Q_MID * (tau_end - tau_K), bool(tau_next > TAU_CAP * tau_K)


def dual_margin(prob, hist, tau_stop, final_ids):
    """The runner-up margins of the non-negative path on the window [0, tau_stop]: the duals u_j (du_j/dtau = <W r phi_j> of the current block, pinned at the wall rho = 1 for the
    block's units) of the lattice kernels that are NOT in the final block `final_ids` (flat lattice indices), at the entry times and at tau_stop. Returns (smallest 1 - u_j over the
    window, 1 - u_j at tau_stop of the smallest, the flat index of the kernel with the smallest margin): the limit is not uniform in the margin (a kernel within the entry shift of
    the wall when a neighbour enters is captured by the network)."""
    P, y, W, lat = prob['P'], prob['y'], prob['W'], prob['lat']
    nlx = lat.shape[1]
    U = np.zeros(lat.shape)
    final = np.zeros(lat.shape, bool)
    for u in np.asarray(final_ids, int):
        final[u // nlx, u % nlx] = True
    best, best_id, mid_best = float('inf'), -1, float('inf')
    for j, h in enumerate(hist):
        taken = np.zeros(lat.shape, bool)
        for b, a in h['ids']:
            U[b, a] = 1.0
            taken[b, a] = True
        t0 = h['tau']
        if j > 0 and t0 >= tau_stop:
            break
        t1 = min(hist[j + 1]['tau'], tau_stop) if j + 1 < len(hist) else tau_stop
        V = lat.corr(W * (y - (P.synth(h['c'], h['w']) if len(h['w']) else 0.0)))
        U = U + (t1 - t0) * np.where(taken, 0.0, V)
        m = np.where(final | taken, np.inf, 1.0 - U)
        k = int(np.argmin(m))
        if m.flat[k] < best:
            best, best_id = float(m.flat[k]), k
        mid_best = float(m.flat[k]) if t1 >= tau_stop else mid_best
    return best, mid_best, float(best_id)


def run_path(prob):
    """Definition 1 on a problem: the path (`sliding_path`, global refits, no exclusion zone, widths fixed, fixed weights), its stop block (BIC for a noisy frame, the first
    block with K_true units for a noise-free one), the block polished, the snapshot times. Returns dict(pv = the PV_FIELDS vector, taus (the entry times), blk_c (K, 2),
    blk_w (K), blk_ids (K,) the lattice units (flat index), blk_tau (K,) their entry times, true_c, true_f, hist_k (units per block), hist_tau)."""
    t0 = time.time()
    P, y, sw, lat, spec = prob['P'], prob['y'], prob['sw'], prob['lat'], prob['spec']
    noisy, K_true = spec['noisy'], spec['K']
    M = K_true + (EXTRA_NOISY if noisy else EXTRA_FREE)
    while True:
        R = IM.sliding_path(P, y, sw, lat, rho=1.0, local=False, r_ex=0.0, order='dual', max_entries=M, tol=PATH_TOL, max_nfev=PATH_NFEV, passes=1)
        hist = R['hist']
        ended = len(R['taus']) < M
        if noisy:
            i, _ = IM.stop_bic(hist, P.N, n_par=3)
            done = ended or len(hist) - 1 - i >= MIN_AFTER
        else:
            idx = [j for j, h in enumerate(hist) if len(h['w']) == K_true]
            i = idx[0] if idx else IM.stop_bic(hist, P.N, n_par=3)[0]
            done = True
        if done or M >= K_true + 30:
            break
        M += EXTRA_NOISY
    h = hist[i]
    K = len(h['w'])
    interior = bool(ended or len(hist) - 1 - i >= MIN_AFTER) if noisy else True
    tau_K = float(h['tau'])
    tau_next = float(hist[i + 1]['tau']) if i + 1 < len(hist) else float('inf')
    tau_end, tau_q1, tau_mid, capped = snapshot_times(tau_K, tau_next)
    taus = np.asarray(R['taus'], float)
    tk = taus[:h['n_ent']]
    min_gap = float(np.min(np.diff(tk) / tk[:-1])) if len(tk) > 1 else float('nan')
    gap_next = (tau_next - tau_K) / tau_K if np.isfinite(tau_next) else float('inf')
    c0_blk = np.array([lat.pos(b, a) for b, a in h['ids']]).reshape(-1, 2)
    box_active = float(np.abs(h['c'] - c0_blk).max() >= 2.0 - 1e-6) if K else 0.0
    ent = {}
    for tau, ij, kind in R['events']:
        if kind == 'enter':
            ent[tuple(ij)] = tau
    blk_tau = np.array([ent[tuple(ij)] for ij in h['ids']])
    c_ref, w_ref, acc, dc, dw = (polish_block(P, y, sw, h['c'], h['w']) if K else (h['c'], h['w'], False, 0.0, 0.0))
    kkt_amp, kkt_pos = block_stationarity(P, prob['W'], y, c_ref, w_ref)
    exact_ok = float('nan')
    if not noisy:
        D = cdist(c_ref, prob['sim']['c']) if K else np.zeros((0, 0))
        exact_ok = float(K == K_true and K > 0 and D.min(0).max() < 1e-3 and D.min(1).max() < 1e-3)
    nlx = lat.shape[1]
    ids_flat = np.array([b * nlx + a for b, a in h['ids']], float)
    margin_min, margin_mid, margin_id = dual_margin(prob, hist, tau_mid, ids_flat) if K else (float('nan'), float('nan'), -1.0)
    pv = dict(K_true=K_true, K_stop=K, i_stop=i, n_hist=len(hist), n_entries=len(taus), interior=float(interior), tau_K=tau_K, tau_next=tau_next, tau_end=tau_end, tau_q1=tau_q1,
              tau_mid=tau_mid, cap_active=float(capped), min_gap=min_gap, gap_next=gap_next, box_active=box_active, polished=float(acc), polish_dc=dc, polish_dw=dw,
              kkt_amp=kkt_amp, kkt_pos=kkt_pos, chi2=float(h['chi2']), exact_ok=exact_ok, n_exits=float(sum(1 for e in R['events'] if e[2] == 'exit')),
              thr=formed_threshold(w_ref) if K else float('nan'), secs=time.time() - t0, margin_min=margin_min, margin_mid=margin_mid, margin_id=margin_id)
    path = dict(pv=np.array([pv[f] for f in PV_FIELDS], float), taus=taus, blk_c=c_ref, blk_w=w_ref, blk_ids=ids_flat,
                blk_tau=blk_tau, true_c=prob['sim']['c'], true_f=prob['sim']['f'], hist_k=np.array([len(q['w']) for q in hist], float),
                hist_tau=np.array([q['tau'] for q in hist], float))
    if K:
        path.update(signed_refs(prob, path))
    return path


STALL_KINDS = ('d1', 'signed', 'bridge')                  # Definition 1 (`imaging.sliding_path`), this module's signed path, the packaged pipeline of the bridge (`sliding_path` again)
STALL_FIELDS = ('calls', 'exits', 'max_repeat', 'stalled', 'fill')


def stall_stats(events, max_entries, n_entries):
    """The enter/exit limit-cycle check of one path run (events = [(tau, (b, a), kind), ...]): (number of exit events, the largest number of events (entries and exits) of one lattice
    point at one identical tau, whether the path is stalled, the entries over max_entries). A STALLED path (a defect of `imaging.sliding_path` reported after the declared runs: one
    lattice point enters and exits again and again at the same entry time until max_entries is reached, because the exit test reads the correlation at the slid centre and the re-entry
    test the one at the lattice start, so that dt = 0 forever) has more than two events of one lattice point at one tau."""
    cnt = collections.Counter((tuple(ij), float('%.12g' % t)) for t, ij, kind in events)
    rep = max(cnt.values()) if cnt else 0
    return float(sum(1 for e in events if e[2] == 'exit')), float(rep), float(rep > 2), n_entries / max_entries


@contextlib.contextmanager
def record_path_calls():
    """Record every call of `imaging.sliding_path` (also inside `fit_frame`) and of this module's `signed_path` while the block runs: yields the list of (kind, events, max_entries,
    number of entries) with kind 'sliding_path' or 'signed_path' and the events as (tau, (b, a), 'enter' / 'exit'). The functions are put back on exit."""
    calls = []
    real_sp, real_sg = IM.sliding_path, signed_path

    def sp(psf, y, sw, lat, *a, **k):
        out = real_sp(psf, y, sw, lat, *a, **k)
        calls.append(('sliding_path', [(t, tuple(ij), kind) for t, ij, kind in out['events']], k.get('max_entries', 40), len(out['taus'])))
        return out

    def sg(psf, y, sw, lat, *a, **k):
        out = real_sg(psf, y, sw, lat, *a, **k)
        calls.append(('signed_path', [(t, tuple(ij), kind) for t, ij, kind, s in out['events']], k.get('max_entries', 40), sum(1 for e in out['events'] if e[2] == 'enter')))
        return out
    IM.sliding_path, globals()['signed_path'] = sp, sg
    try:
        yield calls
    finally:
        IM.sliding_path, globals()['signed_path'] = real_sp, real_sg


def stall_rows(calls_by_kind):
    """The STALL_KINDS x STALL_FIELDS array of the recorded calls {kind: [(call kind, events, max_entries, n_entries), ...]}: calls, exit events, largest repeat, stalled calls, largest fill."""
    out = np.zeros((len(STALL_KINDS), len(STALL_FIELDS)))
    for i, kd in enumerate(STALL_KINDS):
        st = [stall_stats(ev, me, ne) for _, ev, me, ne in calls_by_kind.get(kd, [])]
        if st:
            out[i] = [len(st), sum(s[0] for s in st), max(s[1] for s in st), sum(s[2] for s in st), max(s[3] for s in st)]
    return out


# ---------------------------------------------------------------------------------------------------------------------------------------------------------------
# the network and its flow
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------
def sech2(x):
    """1/cosh^2 x without the cancellation of 1 - tanh^2 (the saturation's derivative a' = sech^2(alpha / a_max))."""
    x = np.minimum(np.abs(np.asarray(x, float)), 300.0)
    return 1.0 / np.cosh(x) ** 2


def lam_of(eps, law=LAW):
    """Lambda(eps) = (pi a_max / (4G)) (2 / (v_s eps))^(1/2): the path time unit, t = Lambda tau (`rich_transit.transit_scale`)."""
    return float(transit_scale(law, eps))


class _Budget(Exception):
    pass


class Network:
    """The paper's network on an imaging problem: units i = 1..d at the start centres c0 (d, 2), state z = (alpha (d), nu (d), x (d), y (d)), w_i = nu_i^2 a_i / m, a_i = a_max
    tanh(alpha_i / a_max), loss F = 1/2 sum_p W_p (y - sum_i w_i phi_i)_p^2 (P: the PixelPSF, high-passed or not). `rhs(t, z)` is -(group rate) dF/dz:
        alpha' = G beta a' g,   nu' = (gamma / 2) nu a g,   c' = G m w (d_x g, d_y g)   (g_i = sum_p W r phi_i, the centre derivative of the same sum),
    `jac(t, z)` a sparse APPROXIMATE Jacobian for the Newton iteration of an implicit solver (see the module docstring); `nfd` counts the evaluations it spends."""

    def __init__(self, P, y, W, c0, eps, law=LAW, act_factor=None):
        self.P, self.y, self.W = P, y, np.asarray(W, float) * np.ones((P.ny, P.nx))
        self.c0 = np.asarray(c0, float).reshape(-1, 2)
        self.d, self.eps = len(self.c0), float(eps)
        self.m, self.G, self.gam, self.A = law.m, law.G, law.gamma, law.A
        self.act_factor = float(ACT_FACTOR if act_factor is None else act_factor)
        self.nrhs = self.nfd = 0
        self.deadline = None

    def init(self):
        """The start: alpha = 0, nu^2 = eps (b_i = 1), the centres on the lattice."""
        return np.concatenate([np.zeros(self.d), np.full(self.d, np.sqrt(self.eps)), self.c0[:, 0], self.c0[:, 1]])

    def split(self, z):
        d = self.d
        return z[:d], z[d:2 * d], z[2 * d:3 * d], z[3 * d:]

    def w_of(self, z):
        al, nu, _, _ = self.split(z)
        return nu * nu * self.A * np.tanh(al / self.A) / self.m

    def centres(self, z):
        return np.c_[z[2 * self.d:3 * self.d], z[3 * self.d:]]

    def pieces(self, z, comps=None):
        """Everything the flow needs at z: (comps, a, a', beta, w, model, v, vx, vy); comps = P.comps(c) (reused if given), model = sum_i w_i phi_i (= P.synth), v_i = sum_p W r
        phi_i and its centre derivatives (= P.corr(c, W r))."""
        al, nu, cx, cy = self.split(z)
        a = self.A * np.tanh(al / self.A)
        ap = sech2(al / self.A)
        beta = nu * nu
        w = beta * a / self.m
        if comps is None:
            comps = self.P.comps(np.c_[cx, cy])
        model = 0.0
        for co, gy, gx, dgy, dgx in comps:
            model = model + co * ((gy * w) @ gx.T)
        Rw = self.W * (self.y - model)
        v = vx = vy = 0.0
        for co, gy, gx, dgy, dgx in comps:
            A_, B_ = gy.T @ Rw, dgy.T @ Rw
            v, vx, vy = v + co * (A_ * gx.T).sum(1), vx + co * (A_ * dgx.T).sum(1), vy + co * (B_ * gx.T).sum(1)
        return comps, a, ap, beta, w, model, v, vx, vy

    def deriv(self, z, pc):
        _, a, ap, beta, w, _, v, vx, vy = pc
        nu = z[self.d:2 * self.d]
        return np.concatenate([self.G * beta * ap * v, 0.5 * self.gam * nu * a * v, self.G * self.m * w * vx, self.G * self.m * w * vy])

    def rhs(self, t, z):
        self.nrhs += 1
        if self.deadline is not None and (self.nrhs & 63) == 0 and time.time() > self.deadline:
            raise _Budget()
        return self.deriv(z, self.pieces(z))

    def _swap_unit(self, comps, u, c_new):
        """comps with the factors of unit u recomputed at the centre c_new (2,): the other units' columns are reused."""
        new = self.P.comps(np.asarray(c_new, float).reshape(1, 2))
        out = []
        for (co, gy, gx, dgy, dgx), (_, ngy, ngx, ndgy, ndgx) in zip(comps, new):
            gy, gx, dgy, dgx = gy.copy(), gx.copy(), dgy.copy(), dgx.copy()
            gy[:, u], gx[:, u], dgy[:, u], dgx[:, u] = ngy[:, 0], ngx[:, 0], ndgy[:, 0], ndgx[:, 0]
            out.append((co, gy, gx, dgy, dgx))
        return out

    def active(self, z):
        """The active units: beta_i > act_factor eps (the formed ones and those on their way)."""
        return np.flatnonzero(self.split(z)[1] ** 2 > self.act_factor * self.eps)

    def jac(self, t, z):
        """The sparse approximate Jacobian d rhs / d z: the own 4 x 4 block of every unit with the drive (g, d_x g, d_y g) frozen (analytic), and for the active units the
        full columns by forward differences (the kernel factors of the other units reused: exact)."""
        d = self.d
        pc = self.pieces(z)
        comps, a, ap, beta, w, _, v, vx, vy = pc
        f0 = self.deriv(z, pc)
        al, nu, _, _ = self.split(z)
        app = -2.0 * a * ap / self.A ** 2
        idx = np.arange(d)
        G, gam, m = self.G, self.gam, self.m
        blocks = [(idx, idx, G * beta * app * v), (idx, d + idx, 2.0 * G * nu * ap * v), (d + idx, idx, 0.5 * gam * nu * ap * v), (d + idx, d + idx, 0.5 * gam * a * v),
                  (2 * d + idx, idx, G * beta * ap * vx), (2 * d + idx, d + idx, 2.0 * G * nu * a * vx), (3 * d + idx, idx, G * beta * ap * vy),
                  (3 * d + idx, d + idx, 2.0 * G * nu * a * vy)]
        R = np.concatenate([b[0] for b in blocks])
        C = np.concatenate([b[1] for b in blocks])
        V = np.concatenate([b[2] for b in blocks])
        act = self.active(z)
        if len(act):
            colmask = np.zeros(4 * d, bool)
            for k in range(4):
                colmask[k * d + act] = True
            keep = ~colmask[C]
            R, C, V = R[keep], C[keep], V[keep]
            rr, cc, vv = [R], [C], [V]
            for u in act:
                for k in range(4):
                    j = k * d + u
                    h = FD_REL * max(abs(z[j]), FD_FLOOR[k])
                    zp = z.copy()
                    zp[j] += h
                    cp = comps if k < 2 else self._swap_unit(comps, u, (zp[2 * d + u], zp[3 * d + u]))
                    col = (self.deriv(zp, self.pieces(zp, comps=cp)) - f0) / h
                    self.nfd += 1
                    nz = np.flatnonzero(col)
                    rr.append(nz)
                    cc.append(np.full(len(nz), j))
                    vv.append(col[nz])
            R, C, V = np.concatenate(rr), np.concatenate(cc), np.concatenate(vv)
        return sp.coo_matrix((V, (R, C)), shape=(4 * d, 4 * d)).tocsc()

    def q_drift(self, z):
        """The drift of the conserved Q_i = beta_i - v_s cosh(2 alpha_i / a_max) (v_s = LAW.c) from its start eps - v_s, relative to the unit's mass: |Q_i - Q_i(0)| / max(beta_i, eps)."""
        al, nu, _, _ = self.split(z)
        q = nu * nu - LAW.c * np.cosh(np.minimum(2.0 * np.abs(al) / self.A, 700.0))
        return np.abs(q - (self.eps - LAW.c)) / np.maximum(nu * nu, self.eps)


class Tracker:
    """The first time each unit's |w| reaches thr, by log-linear interpolation between the accepted steps (nan if never), and the sign of w then."""

    def __init__(self, d, thr):
        self.thr = float(thr)
        self.tform = np.full(d, np.nan)
        self.sform = np.zeros(d)

    def update(self, t0, t1, w0, w1):
        a0, a1 = np.abs(w0), np.abs(w1)
        for i in np.flatnonzero(np.isnan(self.tform) & (a1 >= self.thr)):
            lo = max(a0[i], 1e-300)
            fr = (np.log(self.thr) - np.log(lo)) / (np.log(a1[i]) - np.log(lo)) if (lo < self.thr <= a1[i]) else 1.0
            self.tform[i] = t0 + min(max(fr, 0.0), 1.0) * (t1 - t0)
            self.sform[i] = np.sign(w1[i])


def advance(net, z0, t0, t1, tracker, rtol=RTOL, atol_p=ATOL_P, atol_c=ATOL_C):
    """Integrate the network from (t0, z0) to t1 by Radau with the sparse approximate Jacobian, restarting from the last accepted state with BDF if Radau fails; the tracker
    sees every accepted step. Returns (z, t, status, counts) with status 'ok', 'timeout' or 'failed: <message>' and counts = dict(nfev, njev, nlu, steps, fallback)."""
    d = net.d
    atol = np.concatenate([np.full(2 * d, atol_p * np.sqrt(net.eps)), np.full(2 * d, atol_c)])
    cnt = dict(nfev=0, njev=0, nlu=0, steps=0, fallback=0)
    z, t = np.array(z0, float), float(t0)
    status = 'ok'
    for meth in (Radau, BDF):
        if meth is BDF:
            cnt['fallback'] = 1
        try:
            solver = meth(net.rhs, t, z, t1, rtol=rtol, atol=atol, jac=net.jac)
            wp = net.w_of(z)
            while solver.status == 'running':
                tp = solver.t
                solver.step()
                if solver.status == 'failed':
                    break
                wn = net.w_of(solver.y)
                tracker.update(tp, solver.t, wp, wn)
                wp = wn
                z, t = solver.y.copy(), float(solver.t)
                cnt['steps'] += 1
        except _Budget:
            return z, t, 'timeout', cnt
        cnt['nfev'] += solver.nfev
        cnt['njev'] += solver.njev
        cnt['nlu'] += solver.nlu
        if solver.status == 'finished':
            return z, t1, 'ok', cnt
        status = 'failed: ' + str(solver.message)[:60]
    return z, t, status, cnt


# ---------------------------------------------------------------------------------------------------------------------------------------------------------------
# the comparison of a network state with the path's block
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------
SN_FIELDS = ('n_formed', 'n_unmatched', 'max_d', 'rms_d', 'rel_model', 'nsd_model', 'amp_rms', 'amp_max', 'dormant', 'n_neg', 'neg_max', 'c1', 'c2', 'c3', 'c4', 'ok',
             'id_equal', 'order_ok', 'tdev', 'kkt_amp', 'kkt_pos', 'qd_formed', 'qd_all', 'F_net', 'F_path')
RUN_FIELDS = ('status', 'secs', 'nfev', 'njev', 'nlu', 'steps', 'nfd', 'fallback', 'n_active', 'lam', 'rtol', 'tau_neg_net', 'id_neg_net')
STATUS = {'ok': 0.0, 'timeout': 1.0, 'failed': 2.0, 'error': 3.0}


def compare(P, W, c_net, w_net, c_path, w_path, thr):
    """The network's units (c_net (d, 2), w_net (d)) against the path's block (c_path (K, 2), w_path (K)) at the formed threshold thr: dict with the SN_FIELDS' comparison
    entries: n_formed; n_unmatched (formed units with no block unit within D_NEAR); max_d, rms_d (centre distances of the one-to-one Hungarian match of the formed units to the
    block); rel_model (rms of the model-image difference over the rms of the block's image, all d units of the network), nsd_model (that difference in noise standard
    deviations: rms of sqrt(W) times it); amp_rms, amp_max (relative amplitude differences of the matched pairs); dormant (largest |w| outside the formed set over max|w|);
    n_neg (formed units with w < 0), neg_max (largest negative |w| over max|w|); the criteria c1..c4 and ok."""
    c_net, c_path = np.atleast_2d(np.asarray(c_net, float)), np.atleast_2d(np.asarray(c_path, float))
    w_net, w_path = np.asarray(w_net, float), np.asarray(w_path, float)
    K, d = len(w_path), len(w_net)
    formed = np.flatnonzero(np.abs(w_net) > thr)
    nf = len(formed)
    nan = float('nan')
    out = dict(n_formed=nf, n_unmatched=0, max_d=nan, rms_d=nan, amp_rms=nan, amp_max=nan)
    if nf and K:
        D = cdist(c_net[formed], c_path)
        out['n_unmatched'] = int((D.min(1) >= D_NEAR).sum())
        r, q = linear_sum_assignment(D)
        dd = D[r, q]
        wn, wp = w_net[formed[r]], w_path[q]
        out.update(max_d=float(dd.max()), rms_d=float(np.sqrt(np.mean(dd ** 2))), amp_rms=float(np.linalg.norm(wn - wp) / np.linalg.norm(wp)),
                   amp_max=float(np.max(np.abs(wn - wp) / np.abs(wp))))
    M_net, M_path = P.synth(c_net, w_net), P.synth(c_path, w_path)
    out['rel_model'] = float(np.linalg.norm(M_net - M_path) / np.linalg.norm(M_path))
    out['nsd_model'] = float(np.sqrt(np.mean(np.asarray(W) * (M_net - M_path) ** 2)))
    rest = np.ones(d, bool)
    rest[formed] = False
    wmax = float(np.abs(w_net).max()) if d else nan
    out['dormant'] = float(np.abs(w_net[rest]).max() / wmax) if rest.any() and wmax > 0 else nan
    neg = w_net < 0
    out['n_neg'] = int((neg & (np.abs(w_net) > thr)).sum())
    out['neg_max'] = float(np.abs(w_net[neg]).max() / wmax) if neg.any() and wmax > 0 else 0.0
    out['c1'] = float(nf == K and out['n_unmatched'] == 0)
    out['c2'] = float(nf == K and out['max_d'] <= D_MATCH)
    out['c3'] = float(out['rel_model'] <= REL_TOL and out['amp_max'] <= REL_TOL)
    out['c4'] = float(out['dormant'] <= DORM_TOL)
    out['ok'] = float(out['c1'] and out['c2'] and out['c3'] and out['c4'])
    return out


D_DUP = 0.05                                     # px: POST HOC (not a declared criterion): same-sign formed units closer than this are one unit in the duplicate diagnostic


def merge_formed(c, w, thr, d_dup=D_DUP):
    """POST HOC diagnostic: the network's units (c (d, 2), w (d)) with the same-sign formed units (|w| > thr) that lie within d_dup of each other (single linkage) merged into one unit
    at their |w|-weighted centre with the summed amplitude (the merged-away units get w = 0). Two units at one position are one function: only the sum of their amplitudes is
    determined. Returns (c, w, number of units merged away)."""
    c, w = np.array(c, float), np.array(w, float)
    f = np.flatnonzero(np.abs(w) > thr)
    if len(f) < 2:
        return c, w, 0
    near = (cdist(c[f], c[f]) < d_dup) & (np.sign(w[f])[:, None] == np.sign(w[f])[None, :])
    n_comp, lab = connected_components(near, directed=False)
    gone = 0
    for k in range(n_comp):
        m = f[lab == k]
        if len(m) < 2:
            continue
        a = np.abs(w[m])
        cm, wm = (a[:, None] * c[m]).sum(0) / a.sum(), w[m].sum()
        w[m] = 0.0
        w[m[0]], c[m[0]] = wm, cm
        gone += len(m) - 1
    return c, w, gone


def formation_vs_path(tform, lam, blk_ids, blk_tau):
    """(order_ok, tdev) of the units of the path's block: order_ok = the block's lattice units formed in the order of the path's entry times; tdev = the largest
    |t_formed / (Lambda tau_entry) - 1| (nan if one of them has not formed)."""
    tf = np.asarray(tform, float)[np.asarray(blk_ids, int)]
    tau = np.asarray(blk_tau, float)
    if np.isnan(tf).any():
        return 0.0, float('nan')
    order = float(np.all(np.argsort(tf, kind='stable') == np.argsort(tau, kind='stable')))
    return order, float(np.max(np.abs(tf / (lam * tau) - 1.0)))


def pv_get(path, name):
    """One field of a path dict's PV_FIELDS vector."""
    return float(path['pv'][PV_FIELDS.index(name)])


def ref_nn(path):
    """The non-negative Definition-1 reference of a path dict: the polished BIC (or exact-count) block."""
    return dict(c=path['blk_c'], w=path['blk_w'], ids=path['blk_ids'], tau=path['blk_tau'], thr=pv_get(path, 'thr'))


def ref_signed(path, name):
    """The signed path's block at the snapshot `name` ('q1' or 'mid') as a reference."""
    r = path['s' + name]
    return dict(c=r['c'], w=r['w'], ids=r['ids'], tau=r['tau'], thr=float(r['pv'][SPV_FIELDS.index('thr')]))


def snapshot_stats(prob, ref, net, z, tform, lam):
    """The SN_FIELDS vector of the network state z against a reference block `ref` (c, w, ids, tau, thr): `compare` plus identity, formation, stationarity of the formed units,
    drift of Q (all nan if the reference block is empty)."""
    if len(ref['w']) == 0 or not np.isfinite(ref['thr']):
        return np.full(len(SN_FIELDS), np.nan)
    if not np.all(np.isfinite(z)):
        return np.full(len(SN_FIELDS), np.nan)
    w = net.w_of(z)
    c = net.centres(z)
    cmp_ = compare(prob['P'], prob['W'], c, w, ref['c'], ref['w'], ref['thr'])
    formed = np.flatnonzero(np.abs(w) > ref['thr'])
    ids = np.asarray(ref['ids'], int)
    id_equal = float(set(formed.tolist()) == set(ids.tolist()))
    order_ok, tdev = formation_vs_path(tform, lam, ids, ref['tau'])
    kkt_amp, kkt_pos = block_stationarity(prob['P'], prob['W'], prob['y'], c[formed], w[formed]) if len(formed) else (float('nan'), float('nan'))
    qd = net.q_drift(z)

    def F(cc, ww):
        return 0.5 * float(np.sum(prob['W'] * (prob['y'] - prob['P'].synth(cc, ww)) ** 2))
    vals = dict(cmp_, id_equal=id_equal, order_ok=order_ok, tdev=tdev, kkt_amp=kkt_amp, kkt_pos=kkt_pos, qd_formed=float(qd[formed].max()) if len(formed) else float('nan'),
                qd_all=float(qd.max()), F_net=F(c, w), F_path=F(ref['c'], ref['w']))
    return np.array([vals[f] for f in SN_FIELDS], float)


def run_network(prob, path, eps, rtol=RTOL, budget=BUDGET_S):
    """Train the network on the problem from eps to the path's snapshots tau_q1 and tau_mid (t = Lambda(eps) tau). Returns dict(run = RUN_FIELDS vector; s_q1, s_mid = SN_FIELDS
    vectors against the non-negative block, t_q1, t_mid against the signed path's block at the same times; w_q1, c_q1, w_mid, c_mid; tform (d; the first time |w| reached the formed
    threshold, in units of Lambda, i.e. in tau), sform (d; the sign of w then)). A solver failure or timeout is a recorded result (the snapshots then hold the last state reached)."""
    t0 = time.time()
    net = Network(prob['P'], prob['y'], prob['W'], prob['c0'], eps)
    net.deadline = t0 + budget
    lam = lam_of(eps)
    refs = dict(s_q1=ref_nn(path), s_mid=ref_nn(path), t_q1=ref_signed(path, 'q1'), t_mid=ref_signed(path, 'mid'))
    tr = Tracker(net.d, min(r['thr'] for r in refs.values() if np.isfinite(r['thr'])))
    z, t = net.init(), 0.0
    cnt = dict(nfev=0, njev=0, nlu=0, steps=0, fallback=0)
    status, snaps, n_active = 'ok', {}, 0
    try:
        for name in ('q1', 'mid'):
            if status == 'ok':
                z, t, status, c_ = advance(net, z, t, lam * pv_get(path, 'tau_' + name), tr, rtol=rtol)
                for k in ('nfev', 'njev', 'nlu', 'steps'):
                    cnt[k] += c_[k]
                cnt['fallback'] = max(cnt['fallback'], c_['fallback'])
            snaps[name] = z.copy()
            n_active = max(n_active, len(net.active(z)))
        kind = status.split(':')[0]
    except Exception as e:                                                                  # noqa: BLE001 - recorded as a failed run, not a crash of the part
        kind = 'error'
        for name in ('q1', 'mid'):
            snaps.setdefault(name, z.copy())
        print('network run failed: %s: %s' % (type(e).__name__, e), flush=True)
    neg = np.flatnonzero((tr.sform < 0) & ~np.isnan(tr.tform))
    first = neg[np.argmin(tr.tform[neg])] if len(neg) else -1
    run = dict(status=STATUS.get(kind, 2.0), secs=time.time() - t0, nfev=cnt['nfev'], njev=cnt['njev'], nlu=cnt['nlu'], steps=cnt['steps'], nfd=net.nfd, fallback=cnt['fallback'],
               n_active=n_active, lam=lam, rtol=rtol, tau_neg_net=tr.tform[first] / lam if first >= 0 else float('nan'), id_neg_net=float(first))
    out = dict(run=np.array([run[f] for f in RUN_FIELDS], float), w_q1=net.w_of(snaps['q1']), c_q1=net.centres(snaps['q1']), w_mid=net.w_of(snaps['mid']),
               c_mid=net.centres(snaps['mid']), tform=tr.tform / lam, sform=tr.sform)
    for key, ref in refs.items():
        out[key] = snapshot_stats(prob, ref, net, snaps[key[2:]], tr.tform, lam)
    return out


def classify(s, status, K):
    """The first reason a snapshot (SN_FIELDS vector s, the run's status code, K = the block's number of units) fails the criteria: 'solver', 'extra negative unit', 'extra unit',
    'missing unit', 'unmatched unit', 'position', 'model/amplitude', 'dormant'; 'ok' if it passes."""
    g = dict(zip(SN_FIELDS, s))
    if status != 0.0:
        return 'solver'
    if g['ok'] == 1.0:
        return 'ok'
    if g['c1'] != 1.0:
        if g['n_formed'] > K:
            return 'extra negative unit' if g['n_neg'] > 0 else 'extra unit'
        return 'missing unit' if g['n_formed'] < K else 'unmatched unit'
    if g['c2'] != 1.0:
        return 'position'
    if g['c3'] != 1.0:
        return 'model/amplitude'
    return 'dormant'


# ---------------------------------------------------------------------------------------------------------------------------------------------------------------
# the bridge: the packaged algorithm (`imaging.fit_frame`) against Definition 1, and the switches of its adaptations
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------
BR_CONFIGS = (('full', dict(r_ex=R_EX, local=True, widths=True, fixed=False)), ('no_rex', dict(r_ex=0.0, local=True, widths=True, fixed=False)),
              ('no_local', dict(r_ex=R_EX, local=False, widths=True, fixed=False)), ('no_widths', dict(r_ex=R_EX, local=True, widths=False, fixed=False)),
              ('no_reweight', dict(r_ex=R_EX, local=True, widths=True, fixed=True)), ('off2', dict(r_ex=0.0, local=True, widths=False, fixed=False)),
              ('off3', dict(r_ex=0.0, local=False, widths=False, fixed=False)), ('off4', dict(r_ex=0.0, local=False, widths=False, fixed=True)))
BR_LABEL = {'full': 'full', 'no_rex': 'no exclusion zone', 'no_local': 'no local refits', 'no_widths': 'no width sliding', 'no_reweight': 'no reweighting',
            'off2': 'no widths, no exclusion', 'off3': 'no widths, no exclusion, global refits', 'off4': 'all four off (Definition 1 in the packaged pipeline)'}
BR_FIELDS = ('n', 'K_ref', 'same_count', 'max_d', 'amp_max', 'equal')


def packaged_fit(sim, psf, lat, cfg, k_target):
    """`imaging.fit_frame` with switches (the same steps in the same order: variance model or fixed map, median-subtracted high-passed frame, `sliding_path` stopped three entries
    after the first entrant below 3 sigma, BIC block, `refine`): r_ex, local refits (`fit_frame` has no switch for it), sliding widths (bounds WIDTHS_REL x sigma; BIC with 4
    parameters per unit instead of 3), reweighting (cfg['fixed']: the oracle variance map). Returns {'bic': (c, w), 'cnt': (c, w)}: the refined block at the BIC minimum and the first block
    with k_target units (the same count as Definition 1's block; (zeros, zeros) if the path has none)."""
    ny, nx = sim['ny'], sim['nx']
    fixed = cfg['fixed']
    widths = tuple(r * SIGMA for r in WIDTHS_REL) if cfg['widths'] else None
    nc = 3 if widths else 2
    Pp = IM.PixelPSF(ny, nx, psf.sigma, parts=psf.parts)
    bgv = np.asarray(sim['bg'], float) * np.ones((ny, nx))
    var_a = sim['var'] if fixed else sim['var_a']
    var_b0 = 0.0 if fixed else sim['var_b0']

    def sw_fn(c, w):
        if fixed:
            return 1.0 / np.sqrt(np.asarray(var_a, float) + var_b0)
        m = Pp.synth(c, w) if len(w) else 0.0
        return 1.0 / np.sqrt(var_a * (bgv + np.maximum(m, 0.0)) + var_b0)
    y0 = np.asarray(sim['frame'], float) - bgv
    y = psf.highpass(y0 - np.median(y0))
    sw0 = sw_fn(np.zeros((0, nc)), np.zeros(0))
    R = IM.sliding_path(psf, y, sw0, lat, max_entries=40, sw_fn=None if fixed else sw_fn, sw_on='params', sig_bounds=widths, tol=1e-8, passes=1, rtrunc=4.0, r_ex=cfg['r_ex'],
                        fit_radius=None, order='dual', snr_min=3.0, snr_extra=3, local=cfg['local'])
    out = {}
    for kind in ('bic', 'cnt'):
        if kind == 'bic':
            i = IM.stop_bic(R['hist'], ny * nx, n_par=4 if widths else 3)[0]
        else:
            idx = [j for j, h in enumerate(R['hist']) if len(h['w']) == k_target]
            i = idx[0] if idx else None
        if i is None:
            out[kind] = (np.zeros((0, nc)), np.zeros(0))
            continue
        c, w = R['hist'][i]['c'], R['hist'][i]['w']
        if len(w):
            c, w, _ = IM.refine(psf, y, sw0, c, w, sig_bounds=widths, box=1.0, sw_fn=None if fixed else sw_fn, sw_on='params')
        out[kind] = (c, w)
    return out


def packaged_reference(sim, psf, lat, cfg):
    """`imaging.fit_frame` itself for a configuration it supports (local refits on; reweighting or the fixed map; widths or not; BIC stop): (c, w)."""
    widths = tuple(r * SIGMA for r in WIDTHS_REL) if cfg['widths'] else None
    fixed = cfg['fixed']
    r = IM.fit_frame(sim['frame'], psf, lat, sim['bg'], sim['var'] if fixed else sim['var_a'], 0.0 if fixed else sim['var_b0'], widths=widths, stop=('bic', 4 if widths else 3),
                     max_entries=40, r_ex=cfg['r_ex'], fit_radius=None, rtrunc=4.0, refine_box=1.0, snr_extra=3, tol=1e-8, clean=None, order='dual')
    return r['c'], r['w']


def replica_error(a, b):
    """The largest absolute difference between two (c, w) results (1.0 if their shapes differ)."""
    (c1, w1), (c2, w2) = a, b
    if len(w1) != len(w2) or np.shape(c1) != np.shape(c2):
        return 1.0
    return float(max(np.abs(c1 - c2).max(), np.abs(w1 - w2).max())) if len(w1) else 0.0


def bridge_compare(c_alg, w_alg, c_ref, w_ref):
    """The BR_FIELDS vector of a returned block (c_alg (n, 2 or 3), w_alg (n)) against Definition 1's block: n, K_ref, same count, largest matched centre distance (px, Hungarian
    match on x, y), largest relative amplitude difference of the matched pairs, equal (same count, max_d <= D_MATCH, amp_max <= BR_AMP)."""
    n, K = len(w_alg), len(w_ref)
    nan = float('nan')
    if n == 0 or K == 0:
        return np.array([n, K, float(n == K), nan, nan, float(n == K and n == 0)], float)
    D = cdist(np.atleast_2d(c_alg)[:, :2], c_ref)
    r, q = linear_sum_assignment(D)
    max_d = float(D[r, q].max())
    amp = float(np.max(np.abs(np.asarray(w_alg)[r] - np.asarray(w_ref)[q]) / np.abs(np.asarray(w_ref)[q])))
    return np.array([n, K, float(n == K), max_d, amp, float(n == K and max_d <= D_MATCH and amp <= BR_AMP)], float)


def bridge_frame(prob, path):
    """The bridge on one frame (the hp problem and its Definition-1 path): dict(config -> (2, len(BR_FIELDS)) array [BIC stop; the same count], 'replica' -> the largest difference
    between `packaged_fit` and `fit_frame` over the configurations `fit_frame` supports, 'secs')."""
    t0 = time.time()
    sim, psf, lat = prob['sim'], prob['P'], prob['lat']
    K = int(path['pv'][PV_FIELDS.index('K_stop')])
    res, err = {}, 0.0
    for name, cfg in BR_CONFIGS:
        o = packaged_fit(sim, psf, lat, cfg, K)
        res[name] = np.array([bridge_compare(o[k][0], o[k][1], path['blk_c'], path['blk_w']) for k in ('bic', 'cnt')])
        if cfg['local']:
            err = max(err, replica_error(o['bic'], packaged_reference(sim, psf, lat, cfg)))
    res['replica'] = np.array([err])
    res['secs'] = np.array([time.time() - t0])
    return res


# ---------------------------------------------------------------------------------------------------------------------------------------------------------------
# the raw file and the parts
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------
VARIANTS = (('hp', HP), ('plain', None))
ACC_FRAMES, RTOL_CHECK = ('S02', 'S10', 'A03'), 1e-11


def eps_tag(eps):
    return '%.0e' % float(eps)


def key_path(var, frame, what):
    return 'P|%s|%s|%s' % (var, frame, what)


def key_net(var, eps, frame, what):
    return 'N|%s|%s|%s|%s' % (var, eps_tag(eps), frame, what)


def key_acc(var, eps, frame, what):
    return 'A|%s|%s|%s|%s' % (var, eps_tag(eps), frame, what)


def key_stall(var, frame):
    return 'S|%s|%s' % (var, frame)


def save_part(path, arrays):
    """Merge `arrays` (key -> array) into the raw at `path` (existing keys replaced); written through a temporary file so that an interrupted write loses nothing."""
    old = {}
    if os.path.exists(path):
        with np.load(path, allow_pickle=False) as z:
            old = {k: z[k] for k in z.files}
    old.update({k: np.asarray(v) for k, v in arrays.items()})
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    tmp = path[:-4] + '.tmp.npz'
    np.savez_compressed(tmp, **old)
    os.replace(tmp, path)


def load_raw(path=RAW):
    if not os.path.exists(path):
        return {}
    with np.load(path, allow_pickle=False) as z:
        return {k: z[k] for k in z.files}


def merge_raws(dest, sources):
    """Merge the raws `sources` into `dest` (keys of later sources replace earlier ones); returns the number of keys written."""
    arrays = {}
    for s in sources:
        arrays.update(load_raw(s))
    save_part(dest, arrays)
    return len(arrays)


PATH_KEYS = ('pv', 'taus', 'blk_c', 'blk_w', 'blk_ids', 'blk_tau', 'true_c', 'true_f', 'hist_k', 'hist_tau')


def path_arrays(var, frame, path):
    """The raw keys of a path dict (`run_path`), the signed references' arrays under 'sq1_<key>' and 'smid_<key>'."""
    out = {key_path(var, frame, k): path[k] for k in PATH_KEYS}
    for name in ('sq1', 'smid'):
        if name in path:
            out.update({key_path(var, frame, name + '_' + k): path[name][k] for k in SREF_KEYS})
    return out


def load_path(raw, var, frame):
    """The path dict of `run_path` from the raw (None if absent)."""
    if key_path(var, frame, 'pv') not in raw:
        return None
    path = {k: raw[key_path(var, frame, k)] for k in PATH_KEYS}
    for name in ('sq1', 'smid'):
        if key_path(var, frame, name + '_pv') in raw:
            path[name] = {k: raw[key_path(var, frame, name + '_' + k)] for k in SREF_KEYS}
    return path


NET_KEYS = ('run', 's_q1', 's_mid', 't_q1', 't_mid', 'w_q1', 'c_q1', 'w_mid', 'c_mid', 'tform', 'sform')


def net_arrays(var, eps, frame, res):
    return {key_net(var, eps, frame, k): res[k] for k in NET_KEYS}


def part_P(frames, raw_path, force=False):
    """The paths (Definition 1) of the frames, both variants."""
    raw = load_raw(raw_path)
    for var, hp in VARIANTS:
        for fr in frames:
            if not force and key_path(var, fr, 'pv') in raw:
                continue
            t0 = time.time()
            path = run_path(build(SPECS[fr], hp))
            save_part(raw_path, path_arrays(var, fr, path))
            print('P %s %s: K_stop %d of %d, tau_K %.4f, next %.4f, %.1f s' % (var, fr, pv_get(path, 'K_stop'), SPECS[fr]['K'], pv_get(path, 'tau_K'), pv_get(path, 'tau_next'),
                                                                               time.time() - t0), flush=True)


def part_network(var, eps_list, frames, raw_path, force=False, rtol=RTOL):
    """The network on the frames (variant var), one run per (eps, frame); each run is saved at once and skipped on a rerun."""
    hp = dict(VARIANTS)[var]
    for eps in eps_list:
        for fr in frames:
            raw = load_raw(raw_path)
            if not force and key_net(var, eps, fr, 'run') in raw:
                continue
            prob = build(SPECS[fr], hp)
            path = load_path(raw, var, fr) or run_path(prob)
            arrays = path_arrays(var, fr, path)
            if pv_get(path, 'K_stop') < 1:
                print('N %s %s eps %s: the path stops at the empty block, skipped' % (var, fr, eps_tag(eps)), flush=True)
                save_part(raw_path, arrays)
                continue
            res = run_network(prob, path, eps, rtol=rtol)
            arrays.update(net_arrays(var, eps, fr, res))
            save_part(raw_path, arrays)
            g = dict(zip(SN_FIELDS, res['s_mid']))
            print('N %s %s eps %s: K_stop %d formed %d ok %d rel_model %.2e max_d %.2e status %d steps %d nfev %d %.0f s' % (
                var, fr, eps_tag(eps), pv_get(path, 'K_stop'), g['n_formed'], g['ok'], g['rel_model'], g['max_d'], res['run'][0], res['run'][5], res['run'][2], res['run'][1]), flush=True)


def part_bridge(frames, raw_path, force=False):
    raw = load_raw(raw_path)
    for fr in frames:
        if not force and 'B|%s|full' % fr in raw:
            continue
        prob = build(SPECS[fr], HP)
        path = load_path(raw, 'hp', fr) or run_path(prob)
        res = bridge_frame(prob, path)
        arrays = path_arrays('hp', fr, path)
        arrays.update({'B|%s|%s' % (fr, k): v for k, v in res.items()})
        save_part(raw_path, arrays)
        print('B %s: full config equal at BIC %d, at the same count %d; replica error %.1e; %.1f s' % (fr, res['full'][0, -1], res['full'][1, -1], res['replica'][0], res['secs'][0]), flush=True)


def part_check(frames, eps_list, raw_path, force=False):
    """The solver check: the hp network again at rtol RTOL_CHECK on the frames, to compare with the main run."""
    for eps in eps_list:
        for fr in frames:
            raw = load_raw(raw_path)
            if not force and key_acc('hp', eps, fr, 'run') in raw:
                continue
            prob = build(SPECS[fr], HP)
            path = load_path(raw, 'hp', fr) or run_path(prob)
            res = run_network(prob, path, eps, rtol=RTOL_CHECK)
            save_part(raw_path, {key_acc('hp', eps, fr, k): res[k] for k in ('run', 's_mid', 't_mid', 'w_mid', 'c_mid')})
            print('A %s eps %s: rtol %.0e, %d steps, %.0f s' % (fr, eps_tag(eps), RTOL_CHECK, res['run'][5], res['run'][1]), flush=True)


def part_stall(frames, raw_path, force=False):
    """The stalled-path check (added at the coordinator's request after the declared runs; it recomputes the paths of the declared runs and reads their events): for each frame and
    variant the calls of `imaging.sliding_path` of `run_path` (Definition 1), of this module's `signed_path`, and (hp) of the bridge's packaged pipeline, as STALL_KINDS x STALL_FIELDS."""
    for var, hp in VARIANTS:
        for fr in frames:
            raw = load_raw(raw_path)
            if not force and key_stall(var, fr) in raw:
                continue
            prob = build(SPECS[fr], hp)
            by_kind = {}
            with record_path_calls() as calls:
                path = run_path(prob)
            by_kind['d1'] = [c for c in calls if c[0] == 'sliding_path']
            by_kind['signed'] = [c for c in calls if c[0] == 'signed_path']
            if var == 'hp':
                with record_path_calls() as calls:
                    bridge_frame(prob, path)
                by_kind['bridge'] = [c for c in calls if c[0] == 'sliding_path']
            rows = stall_rows(by_kind)
            save_part(raw_path, {key_stall(var, fr): rows})
            print('S %s %s: %d path calls, %d exit events, largest repeat of one point at one tau %d, stalled %d' % (var, fr, rows[:, 0].sum(), rows[:, 1].sum(), rows[:, 2].max(),
                                                                                                                     rows[:, 3].sum()), flush=True)


# ---------------------------------------------------------------------------------------------------------------------------------------------------------------
# the report (everything below is printed from the raw alone)
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------
GROUPS = (('SMLM, noisy', [n for n in FRAMES if n[0] == 'S' and not n.endswith('n')]), ('SMLM, noise-free', [n for n in FRAMES if n[0] == 'S' and n.endswith('n')]),
          ('astronomy, noisy', [n for n in FRAMES if n[0] == 'A' and not n.endswith('n')]), ('astronomy, noise-free', [n for n in FRAMES if n[0] == 'A' and n.endswith('n')]),
          ('all frames', list(FRAMES)))
VAR_LABEL = {'hp': 'high-passed (hp = %g px)' % HP, 'plain': 'plain Gaussian dictionary (no high-pass)'}


def _x(v, f='%.2g'):
    """A number for a table cell: 'n/a' for nan, 'inf' for infinity."""
    v = float(v)
    if np.isnan(v):
        return 'n/a'
    if np.isinf(v):
        return 'inf'
    return f % v


def _tab(head, rows):
    return ['| ' + ' | '.join(head) + ' |', '|' + '|'.join(['---'] * len(head)) + '|'] + ['| ' + ' | '.join(str(c) for c in r) + ' |' for r in rows]


def _sn(raw, var, eps, fr, which):
    """(SN_FIELDS dict, run dict) of a stored run, or (None, None)."""
    k = key_net(var, eps, fr, which)
    if k not in raw:
        return None, None
    return dict(zip(SN_FIELDS, raw[k])), dict(zip(RUN_FIELDS, raw[key_net(var, eps, fr, 'run')]))


def _flags(g):
    return ''.join('Y' if g[k] == 1.0 else 'n' for k in ('c1', 'c2', 'c3', 'c4'))


def _frame_desc(fr):
    s = SPECS[fr]
    if s['kind'] == 'smlm':
        return 'SMLM %s%s%s' % (s['level'], ', noise-free' if not s['noisy'] else '', ', pair %.1f sigma' % s['pair'] if s['pair'] else '')
    return 'astro%s%s' % (', noise-free' if not s['noisy'] else '', ', satellite %.1f sigma' % SAT_SIGMAS if s['sat'] else '')


def eps_present(raw):
    """The eps of the network runs in the raw, largest first."""
    return sorted({float(k.split('|')[2]) for k in raw if k.startswith('N|')}, reverse=True)


def _frames_have(raw, var, kind='pv'):
    return [fr for fr in FRAMES if key_path(var, fr, kind) in raw]


def sec_paths(raw, var):
    frs = _frames_have(raw, var)
    L = ['### Definition 1 (non-negative path), %s' % VAR_LABEL[var], '']
    if not frs:
        return L + ['(not run)', '']
    rows = []
    for fr in frs:
        p = load_path(raw, var, fr)
        g = dict(zip(PV_FIELDS, p['pv']))
        sq = dict(zip(SPV_FIELDS, p['smid']['pv'])) if 'smid' in p else None
        rows.append([fr, _frame_desc(fr), int(g['K_true']), int(g['K_stop']), _x(g['tau_K'], '%.4g'), _x(g['tau_next'] / g['tau_K'], '%.3f'), _x(g['min_gap'], '%.3g'),
                     'cap' if g['cap_active'] else '-', 'yes' if g['interior'] else 'no', _x(g['exact_ok'], '%d') if not SPECS[fr]['noisy'] else '-',
                     _x(g['margin_min'], '%.2e'), _x(g['margin_mid'], '%.2e'), _x(sq['tau_neg'] / g['tau_K'], '%.3f') if sq else 'n/a', 'yes' if g['box_active'] else 'no',
                     _x(g['polish_dc'], '%.1e'), _x(g['kkt_pos'], '%.1e')])
    L += _tab(['frame', 'kind', 'K true', 'K stop', 'tau_K', 'tau_next / tau_K', 'smallest relative gap between the first K entries', 'cap', 'BIC minimum interior',
               'recovers the emitters (noise-free)', 'runner-up margin: smallest 1 - u over the window', 'runner-up margin at tau_mid', 'first negative entry of the signed path / tau_K',
               'centre box active', 'polish: largest centre shift (px)',
               'block stationarity: Newton displacement (px)'], rows)
    L.append('')
    return L


def sec_results(raw, var, eps, which='mid', ref='s'):
    """The per-frame table of one (variant, eps, snapshot, reference): ref 's' the non-negative block, 't' the signed path's block at the snapshot."""
    frs = [fr for fr in FRAMES if key_net(var, eps, fr, ref + '_' + which) in raw]
    ttl = 'tau_mid' if which == 'mid' else 'tau_q1'
    L = ['### %s, eps = %s, at %s, against %s' % (VAR_LABEL[var], eps_tag(eps), ttl, 'Definition 1 (non-negative path)' if ref == 's' else 'the signed path at the same time'), '']
    if not frs:
        return L + ['(not run)', '']
    rows = []
    for fr in frs:
        g, run = _sn(raw, var, eps, fr, ref + '_' + which)
        p = load_path(raw, var, fr)
        K = int(p['pv'][PV_FIELDS.index('K_stop')]) if ref == 's' else int(p['s' + which]['pv'][SPV_FIELDS.index('K')])
        if np.isnan(g['n_formed']):
            rows.append([fr, K, 'n/a'] + [''] * 11)
            continue
        rows.append([fr, K, '%d (%d)' % (g['n_formed'], g['n_neg']), _x(g['n_unmatched'], '%d'), _x(g['max_d'], '%.2e'), _x(g['rel_model'], '%.2e'), _x(g['amp_max'], '%.2e'),
                     _x(g['dormant'], '%.1e'), _x(g['kkt_pos'], '%.1e'), _flags(g), 'pass' if g['ok'] == 1.0 else classify([g[f] for f in SN_FIELDS], run['status'], K), _x(g['tdev'], '%.1e'),
                     _x(run['steps'], '%d'), _x(run['secs'], '%.0f')])
    L += _tab(['frame', 'K block', 'formed (negative)', 'formed units with no block unit within 0.5 px', 'largest matched centre distance (px)', 'relative rms model difference',
               'largest relative amplitude difference', 'dormant max / max', 'stationarity of the formed units: Newton displacement (px)', 'C1 C2 C3 C4', 'result', 'largest relative formation-time offset', 'steps', 'seconds'], rows)
    L.append('')
    return L


def _count(raw, var, eps, frames, key):
    n = k = 0
    for fr in frames:
        g, _ = _sn(raw, var, eps, fr, key)
        if g is None or np.isnan(g['ok']):
            continue
        n += 1
        k += int(g['ok'] == 1.0)
    return k, n


def sec_summary(raw, epss):
    L = ['## Summary: the network forms the block of Definition 1 (C1-C4 all hold)', '',
         'Counts of runs in which all four criteria hold, "k of n", n the runs present. Columns: the non-negative path (the brief) at tau_mid and at tau_q1, and the signed path (the '
         'limit of the signed network) at the same times.', '']
    rows = []
    for var, _ in VARIANTS:
        for eps in epss:
            for gname, frames in GROUPS:
                cells = []
                for key in ('s_mid', 's_q1', 't_mid', 't_q1'):
                    k, n = _count(raw, var, eps, frames, key)
                    cells.append('%d of %d' % (k, n) if n else 'n/a')
                if any(c != 'n/a' for c in cells):
                    rows.append([var, eps_tag(eps), gname] + cells)
    L += _tab(['variant', 'eps', 'frames', 'non-negative path, tau_mid', 'non-negative path, tau_q1', 'signed path, tau_mid', 'signed path, tau_q1'], rows)
    L.append('')
    return L


def sec_eps(raw, var, epss):
    """The discrepancy against eps: for the frames run at the eps present, the relative model difference, the largest centre distance and the result."""
    L = ['### The discrepancy against eps, %s, non-negative path at tau_mid' % VAR_LABEL[var], '']
    if len(epss) < 2:
        return L + ['(needs two eps)', '']
    rows = []
    for fr in FRAMES:
        gs = [_sn(raw, var, e, fr, 's_mid')[0] for e in epss]
        if sum(g is not None for g in gs) < 2:
            continue
        have = [g for g in gs if g is not None]
        rows.append([fr] + [_x(g['rel_model'], '%.2e') if g is not None else 'n/a' for g in gs] + [_x(g['max_d'], '%.2e') if g is not None else 'n/a' for g in gs] +
                    [('pass' if g['ok'] == 1.0 else 'fail') if g is not None else 'n/a' for g in gs] + [_x(have[0]['rel_model'] / have[-1]['rel_model'], '%.3g')])
    if not rows:
        return L + ['(needs two eps on the same frame)', '']
    L += _tab(['frame'] + ['rel. model, eps = %s' % eps_tag(e) for e in epss] + ['max centre distance (px), eps = %s' % eps_tag(e) for e in epss] + ['result, eps = %s' % eps_tag(e) for e in epss] +
              ['ratio of the relative model difference, largest eps / smallest eps'], rows)
    L.append('')
    return L


PRED_THRESHOLDS = (1.0, 2.0, 4.0, 8.0, 16.0, 32.0)


def sec_predictor(raw, var, epss):
    """Does the runner-up margin predict the failures? For each eps and threshold t: the frames whose smallest margin is <= t eps^(1/3) (the scale of the entry shift) and how many of
    them fail at tau_mid, the same for the frames above."""
    L = ['### Failures against the runner-up margin, %s' % VAR_LABEL[var], '',
         'The limit is not uniform in the margin: a lattice kernel that is within the entry shift of the wall when a neighbour enters is captured by the network. The margin is the '
         'smallest 1 - u_j of the kernels outside the block over the window (Definition 1 alone); the entry shift scales as eps^(1/3). "Fail": not all of C1-C4 at tau_mid.', '']
    rows = []
    for eps in epss:
        data = []
        for fr in FRAMES:
            g, _ = _sn(raw, var, eps, fr, 's_mid')
            if g is None or np.isnan(g['ok']):
                continue
            data.append((pv_get(load_path(raw, var, fr), 'margin_min') / eps ** (1.0 / 3.0), g['ok'] == 1.0))
        for t in PRED_THRESHOLDS:
            lo, hi = [ok for r, ok in data if r <= t], [ok for r, ok in data if r > t]
            rows.append([eps_tag(eps), '%g' % t, len(lo), sum(not ok for ok in lo), len(hi), sum(not ok for ok in hi)])
    if not rows:
        return L + ['(not run)', '']
    L += _tab(['eps', 't', 'frames with margin <= t eps^(1/3)', 'of which fail', 'frames above', 'of which fail'], rows)
    L.append('')
    return L


def sec_failures(raw, var, eps):
    """The runs that fail at tau_mid against the non-negative block: the reason, and what predicts it (the entry gaps, the first negative entry of the signed path)."""
    L = ['### Failures at tau_mid, %s, eps = %s' % (VAR_LABEL[var], eps_tag(eps)), '']
    rows = []
    for fr in FRAMES:
        g, run = _sn(raw, var, eps, fr, 's_mid')
        if g is None or g['ok'] == 1.0 or np.isnan(g['ok']):
            continue
        p = load_path(raw, var, fr)
        pv = dict(zip(PV_FIELDS, p['pv']))
        gq, _ = _sn(raw, var, eps, fr, 's_q1')
        gt, _ = _sn(raw, var, eps, fr, 't_mid')
        sp = dict(zip(SPV_FIELDS, p['smid']['pv'])) if 'smid' in p else None
        rows.append([fr, classify([g[f] for f in SN_FIELDS], run['status'], int(pv['K_stop'])), _x(pv['min_gap'], '%.3g'), _x(pv['margin_min'] / eps ** (1.0 / 3.0), '%.3g'), _x(pv['gap_next'], '%.3g'),
                     _x(sp['tau_neg'] / pv['tau_K'], '%.3f') if sp else 'n/a', 'pass' if gq['ok'] == 1.0 else 'fail', 'pass' if gt['ok'] == 1.0 else 'fail'])
    if not rows:
        return L + ['none', '']
    L += _tab(['frame', 'reason (first that applies)', 'smallest relative entry gap among the first K', 'runner-up margin / eps^(1/3)', 'relative gap to the next entry',
               'signed path: first negative entry / tau_K',
               'non-negative path at tau_q1', 'signed path at tau_mid'], rows)
    L.append('')
    return L


def sec_duplicates(raw, var, eps):
    """POST HOC (not a declared criterion): the runs that fail at tau_mid against the non-negative block with MORE formed units than the block, after merging the same-sign formed units
    within D_DUP of each other (two units at one position are one function); is the declared criterion met then?"""
    L = ['### Duplicated units (post hoc), %s, eps = %s' % (VAR_LABEL[var], eps_tag(eps)), '']
    rows = []
    hp = dict(VARIANTS)[var]
    for fr in FRAMES:
        g, run = _sn(raw, var, eps, fr, 's_mid')
        if g is None or np.isnan(g['ok']) or g['ok'] == 1.0:
            continue
        p = load_path(raw, var, fr)
        K = len(p['blk_w'])
        if g['n_formed'] <= K:
            continue
        prob = build(SPECS[fr], hp)
        thr = pv_get(p, 'thr')
        c, w, gone = merge_formed(raw[key_net(var, eps, fr, 'c_mid')], raw[key_net(var, eps, fr, 'w_mid')], thr)
        s = compare(prob['P'], prob['W'], c, w, p['blk_c'], p['blk_w'], thr)
        f = np.flatnonzero(np.abs(raw[key_net(var, eps, fr, 'w_mid')]) > thr)
        cf = raw[key_net(var, eps, fr, 'c_mid')][f]
        dm = cdist(cf, cf) + np.eye(len(f)) * 1e9
        rows.append([fr, K, int(g['n_formed']), gone, int(s['n_formed']), _x(dm.min(), '%.2e'), _x(s['max_d'], '%.2e'), _x(s['rel_model'], '%.2e'), _x(s['amp_max'], '%.2e'),
                     'pass' if s['ok'] == 1.0 else 'fail'])
    if not rows:
        return L + ['none', '']
    L += _tab(['frame', 'K block', 'formed', 'merged away', 'formed after the merge', 'closest pair of formed units (px)', 'largest matched centre distance (px)',
               'relative rms model difference', 'largest relative amplitude difference', 'C1-C4 after the merge'], rows)
    L.append('')
    return L


def sec_dormant(raw, var, epss):
    """(iii) of the brief: the dormant units against the theory. The largest |w| outside the formed set (absolute), the coefficient the theory gives to an atom whose dual is at the
    runner-up margin mu of the path at tau_mid (`FreeLaw.core_tail`: a_max eps^(3/2) cos(pi mu / 2) / (m (2 v_s)^(1/2) sin^3(pi mu / 2)), leading order in eps), and their ratio."""
    L = ['### Dormant units against the theory, %s' % VAR_LABEL[var], '']
    head, rows = ['frame', 'runner-up margin at tau_mid'], []
    for fr in FRAMES:
        row = [fr]
        for eps in epss:
            g, _ = _sn(raw, var, eps, fr, 's_mid')
            if g is None or np.isnan(g['dormant']):
                row += ['n/a'] * 3
                continue
            mu = pv_get(load_path(raw, var, fr), 'margin_mid')
            obs = g['dormant'] * float(np.abs(raw[key_net(var, eps, fr, 'w_mid')]).max())
            pred = float(LAW.core_tail(eps, mu)) if 0.0 < mu < 1.0 else float('nan')
            row += [_x(obs, '%.2e'), _x(pred, '%.2e'), _x(obs / pred, '%.2g')]
        if any(c != 'n/a' for c in row[1:]):
            rows.append([row[0], _x(pv_get(load_path(raw, var, fr), 'margin_mid'), '%.2e')] + row[1:])
    if not rows:
        return L + ['(not run)', '']
    for eps in epss:
        head += ['largest dormant |w|, eps = %s' % eps_tag(eps), 'theory, eps = %s' % eps_tag(eps), 'ratio, eps = %s' % eps_tag(eps)]
    L += _tab(head, rows)
    L.append('')
    return L


def sec_scaling(raw, var, epss):
    """The size of the discrepancy of the runs that pass (C1-C4 at tau_mid against the non-negative block) against eps: medians and maxima over the noisy and the noise-free frames, the
    identity of the formed lattice units and the order of formation, and the formation-time offset in units of eps^(1/2)."""
    L = ['### Size of the discrepancy of the runs that pass, %s' % VAR_LABEL[var], '']
    groups = (('noisy', [n for n in FRAMES if SPECS[n]['noisy']]), ('noise-free', [n for n in FRAMES if not SPECS[n]['noisy']]))
    rows = []
    for eps in epss:
        for gname, frames in groups:
            gs = [g for g in (_sn(raw, var, eps, fr, 's_mid')[0] for fr in frames) if g is not None and g['ok'] == 1.0]
            if not gs:
                continue
            col = lambda k: np.array([g[k] for g in gs])
            td = col('tdev')
            rows.append([eps_tag(eps), gname, len(gs), _x(np.median(col('rel_model')), '%.1e'), _x(col('rel_model').max(), '%.1e'), _x(np.median(col('max_d')), '%.1e'), _x(col('max_d').max(), '%.1e'),
                         _x(np.median(col('amp_max')), '%.1e'), _x(col('amp_max').max(), '%.1e'), _x(np.median(col('dormant')), '%.1e'), _x(col('dormant').max(), '%.1e'),
                         '%d of %d' % (col('id_equal').sum(), len(gs)), '%d of %d' % (col('order_ok').sum(), len(gs)), _x(np.nanmedian(td), '%.1e') if np.isfinite(td).any() else 'n/a',
                         _x(np.nanmedian(td) / eps ** 0.5, '%.3g') if np.isfinite(td).any() else 'n/a'])
    if not rows:
        return L + ['(none passes)', '']
    L += _tab(['eps', 'frames', 'runs that pass', 'median rel. model difference', 'max', 'median largest centre distance (px)', 'max', 'median largest rel. amplitude difference', 'max',
               'median dormant max / max', 'max', 'same lattice units as the block', 'same formation order', 'median largest formation-time offset', 'that / eps^(1/2)'], rows)
    L.append('')
    return L


def sec_negative(raw, var, eps):
    """The negative entries: the signed path enters a negative unit before tau_mid in these frames; does the network form the same lattice unit, with its sign, at the same time?"""
    L = ['### Negative units, %s, eps = %s' % (VAR_LABEL[var], eps_tag(eps)), '']
    rows = []
    for fr in FRAMES:
        k = key_net(var, eps, fr, 'run')
        if k not in raw:
            continue
        p = load_path(raw, var, fr)
        sp = dict(zip(SPV_FIELDS, p['smid']['pv'])) if 'smid' in p else None
        run = dict(zip(RUN_FIELDS, raw[k]))
        if sp is None or not np.isfinite(sp['tau_neg']):
            continue
        tp, tn = sp['tau_neg'], run['tau_neg_net']
        rows.append([fr, _x(tp, '%.4f'), int(sp['id_neg']), _x(tn, '%.4f'), int(run['id_neg_net']), _x(abs(tn - tp) / tp, '%.1e') if np.isfinite(tn) else 'n/a',
                     'yes' if sp['id_neg'] == run['id_neg_net'] else 'no'])
    if not rows:
        return L + ['no frame has a negative entry of the signed path before tau_mid', '']
    L += _tab(['frame', 'signed path: first negative entry (tau)', 'lattice unit', 'network: first negative unit formed (tau)', 'lattice unit', 'relative time difference', 'same unit'], rows)
    L.append('')
    return L


def sec_stall(raw):
    """The stalled-path check: one sentence with the totals, and a table of the defective paths if there are any (a stalled path is not a path-versus-network result)."""
    L = ['## Stalled paths (the enter/exit limit cycle of `imaging.sliding_path`)', '']
    keys = [(var, fr) for var, _ in VARIANTS for fr in FRAMES if key_stall(var, fr) in raw]
    if not keys:
        return L + ['(not run)', '']
    tot = sum(raw[key_stall(v, f)] for v, f in keys)
    kinds = ', '.join('%d %s' % (tot[i, 0], STALL_LABEL[k]) for i, k in enumerate(STALL_KINDS))
    L.append('Check of the paths of the declared runs (%d frame-variant pairs; %d path calls: %s): %d exit events in all, the largest number of events of one lattice point at one '
             'identical tau is %d (a stalled path has more than 2), stalled paths: %d%s' % (len(keys), tot[:, 0].sum(), kinds, tot[:, 1].sum(), max(raw[key_stall(v, f)][:, 2].max() for v, f in keys),
                                                                                              tot[:, 3].sum(), '.' if tot[:, 3].sum() == 0 else ' (listed below; defective paths, not path-versus-network results).'))
    L.append('')
    bad = [(v, f, STALL_KINDS[i]) for v, f in keys for i in range(len(STALL_KINDS)) if raw[key_stall(v, f)][i, 3] > 0]
    if bad:
        L += _tab(['variant', 'frame', 'path', 'stalled calls', 'largest repeat'], [[v, f, STALL_LABEL[k], int(raw[key_stall(v, f)][STALL_KINDS.index(k), 3]),
                                                                                  int(raw[key_stall(v, f)][STALL_KINDS.index(k), 2])] for v, f, k in bad])
        L.append('')
    return L


STALL_LABEL = {'d1': 'of Definition 1', 'signed': 'of the signed path', 'bridge': 'of the bridge\'s packaged pipeline'}


def sec_solver(raw, epss):
    L = ['## The solver', '']
    rows = []
    for var, _ in VARIANTS:
        for eps in epss:
            runs = [dict(zip(RUN_FIELDS, raw[key_net(var, eps, fr, 'run')])) for fr in FRAMES if key_net(var, eps, fr, 'run') in raw]
            if not runs:
                continue
            col = lambda k: np.array([r[k] for r in runs])
            qd = [max(dict(zip(SN_FIELDS, raw[key_net(var, eps, fr, 's_mid')]))['qd_all'], dict(zip(SN_FIELDS, raw[key_net(var, eps, fr, 's_q1')]))['qd_all'])
                  for fr in FRAMES if key_net(var, eps, fr, 'run') in raw]
            rows.append([var, eps_tag(eps), len(runs), int((col('status') == 0).sum()), int(col('fallback').sum()), _x(np.median(col('secs')), '%.0f'), _x(col('secs').max(), '%.0f'),
                         _x(np.median(col('steps')), '%.0f'), _x(col('steps').max(), '%.0f'), _x(np.median(col('nfev')), '%.0f'), _x(np.median(col('njev')), '%.0f'),
                         _x(np.median(col('nlu')), '%.0f'), _x(np.median(col('nfd')), '%.0f'), _x(col('n_active').max(), '%d'), _x(np.nanmax(qd), '%.1e')])
    L += _tab(['variant', 'eps', 'runs', 'finished', 'BDF fallbacks', 'median seconds', 'max seconds', 'median steps', 'max steps', 'median right-hand sides', 'median Jacobians',
               'median factorisations', 'median finite-difference evaluations', 'most active units', 'largest drift of Q (relative)'], rows) if rows else ['(not run)']
    L.append('')
    acc = [(eps, fr) for eps in epss for fr in FRAMES if key_acc('hp', eps, fr, 'run') in raw and key_net('hp', eps, fr, 'run') in raw]
    L += ['### Tolerance check (rtol %.0e against the main runs at rtol %.0e)' % (RTOL_CHECK, RTOL), '']
    rows = []
    for eps, fr in acc:
        P = build(SPECS[fr], HP)['P']
        wa, ca = raw[key_acc('hp', eps, fr, 'w_mid')], raw[key_acc('hp', eps, fr, 'c_mid')]
        wb, cb = raw[key_net('hp', eps, fr, 'w_mid')], raw[key_net('hp', eps, fr, 'c_mid')]
        Ma, Mb = P.synth(ca, wa), P.synth(cb, wb)
        rows.append([fr, eps_tag(eps), _x(np.linalg.norm(Ma - Mb) / np.linalg.norm(Mb), '%.1e'), _x(np.abs(wa - wb).max() / np.abs(wb).max(), '%.1e'), _x(np.abs(ca - cb).max(), '%.1e'),
                     _x(raw[key_acc('hp', eps, fr, 'run')][RUN_FIELDS.index('steps')], '%d'), _x(raw[key_net('hp', eps, fr, 'run')][RUN_FIELDS.index('steps')], '%d')])
    L += _tab(['frame', 'eps', 'relative rms model difference of the final states', 'largest |dw| / max|w|', 'largest |dc| (px)', 'steps at the tight tolerance', 'steps at the main'],
              rows) if rows else ['(not run)']
    L.append('')
    return L


def sec_bridge(raw):
    L = ['## The packaged algorithm (`imaging.fit_frame`) against Definition 1', '']
    frs = [fr for fr in FRAMES if 'B|%s|full' % fr in raw]
    if not frs:
        return L + ['(not run)', '']
    names = [n for n, _ in BR_CONFIGS]
    L += ['The high-passed frames. `full`: r_ex = %g px, local refits, widths sliding in [%g, %g] sigma, weights iterated from the fit (the packaged settings, BIC stop with 4 parameters per unit); '
          'the other configurations switch adaptations off (`off2`, `off3`, `off4` cumulatively). Equal: same count as the block, largest matched distance <= %g px, largest relative amplitude '
          'difference <= %g. Largest difference between this script\'s replica of `fit_frame` and `fit_frame` itself: %.1e.' % (
              R_EX, WIDTHS_REL[0], WIDTHS_REL[1], D_MATCH, BR_AMP, max(float(raw['B|%s|replica' % fr][0]) for fr in frs)), '']
    for k, ttl in ((0, 'at the algorithm\'s own BIC stop'), (1, 'at the same number of units as Definition 1\'s block')):
        L += ['### Equal to Definition 1\'s block, %s' % ttl, '']
        rows = []
        for fr in frs:
            row = [fr, int(raw['B|%s|full' % fr][k][BR_FIELDS.index('K_ref')])]
            for n in names:
                g = dict(zip(BR_FIELDS, raw['B|%s|%s' % (fr, n)][k]))
                row.append('equal' if g['equal'] == 1.0 else ('n=%d' % g['n'] if g['same_count'] != 1.0 else ('d=%s' % _x(g['max_d'], '%.2g') if not g['max_d'] <= D_MATCH else 'a=%s' % _x(g['amp_max'], '%.2g'))))
            rows.append(row)
        tot = ['equal in', '']
        for j, n in enumerate(names):
            tot.append('%d of %d' % (sum(int(raw['B|%s|%s' % (fr, n)][k][BR_FIELDS.index('equal')]) for fr in frs), len(frs)))
        L += _tab(['frame', 'K block'] + [BR_LABEL[n] for n in names], rows + [tot])
        L.append('')
    L += ['### The full packaged algorithm against the block, numbers', '']
    rows = []
    for fr in frs:
        a, b = (dict(zip(BR_FIELDS, raw['B|%s|full' % fr][k])) for k in (0, 1))
        rows.append([fr, int(a['K_ref']), int(a['n']), _x(a['max_d'], '%.2e'), _x(a['amp_max'], '%.2e'), int(b['n']), _x(b['max_d'], '%.2e'), _x(b['amp_max'], '%.2e')])
    L += _tab(['frame', 'K block', 'units at its BIC stop', 'largest distance (px)', 'largest relative amplitude difference', 'units at the same count', 'largest distance (px)',
               'largest relative amplitude difference'], rows)
    L.append('')
    return L


def report(path=RAW):
    """The markdown report, from the raw alone."""
    raw = load_raw(path)
    L = ['# Does the factorised RBF network, trained by gradient flow on an imaging loss, form the units that the sliding path returns?', '',
         'Generated by `experiments/exp_network_imaging.py` from `network_imaging_raw.npz`; every number below is printed by the script. Design, frames, stop rule and acceptance '
         'criteria: the script\'s docstring. The network: m = %g, G = %g, gamma = %g, a_max = %g, equal walls; Radau (BDF fallback), rtol %.0e, atol %.0e sqrt(eps) on (alpha, nu) and %.0e on '
         'the centres, sparse approximate Jacobian; eps = %s. PSF sigma = %g px, lattice spacing %g px and margin %g px. "Formed": |w| > min(%g max|w_path|, %g min|w_path|). '
         'C1: exactly K formed units, each within %g px of a block unit; C2: one-to-one match within %g px; C3: relative rms model difference and largest relative amplitude '
         'difference <= %g; C4: dormant max / max <= %g. tau_q1 and tau_mid are the quarter and the midpoint of the way from tau_K to min(tau_next, %g tau_K).\n' % (
             LAW.m, LAW.G, LAW.gamma, LAW.A, RTOL, ATOL_P, ATOL_C, ', '.join(eps_tag(e) for e in eps_present(raw)), SIGMA, LAT_H, LAT_MARGIN, THETA, THETA_MIN, D_NEAR, D_MATCH, REL_TOL,
             DORM_TOL, TAU_CAP)]
    epss = eps_present(raw)
    L += sec_summary(raw, epss)
    L += ['## The paths', ''] + sec_paths(raw, 'hp') + sec_paths(raw, 'plain')
    for var, _ in VARIANTS:
        if not any(key_net(var, eps, fr, 'run') in raw for eps in epss for fr in FRAMES):
            continue
        L += ['## The network, %s' % VAR_LABEL[var], '']
        for eps in epss:
            L += sec_results(raw, var, eps, 'mid', 's') + sec_failures(raw, var, eps) + sec_duplicates(raw, var, eps) + sec_negative(raw, var, eps)
        L += sec_eps(raw, var, epss) + sec_scaling(raw, var, epss) + sec_predictor(raw, var, epss) + sec_dormant(raw, var, epss)
        for eps in epss:
            L += sec_results(raw, var, eps, 'q1', 's') + sec_results(raw, var, eps, 'mid', 't')
    L += sec_stall(raw) + sec_solver(raw, epss) + sec_bridge(raw)
    return '\n'.join(L) + '\n'


# ---------------------------------------------------------------------------------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--only', help='parts to run: P (paths), H (network, hp), Q (network, plain), B (bridge), A (solver check), S (stalled-path check); comma separated')
    ap.add_argument('--eps', default=','.join(eps_tag(e) for e in EPS_LIST), help='the eps of the network parts (default: %(default)s)')
    ap.add_argument('--frames', default=','.join(FRAMES), help='the frames to run (default: all)')
    ap.add_argument('--force', action='store_true', help='recompute what is already in the raw')
    ap.add_argument('--merge', help='comma-separated raws to merge into this one')
    ap.add_argument('--from-raw', action='store_true', help='regenerate the md from the raw')
    A = ap.parse_args()
    frames = [f for f in A.frames.split(',') if f]
    for f in frames:
        if f not in SPECS:
            ap.error('unknown frame %r' % f)
    eps_list = [float(e) for e in A.eps.split(',') if e]
    if A.merge:
        print('merged %d keys' % merge_raws(RAW, A.merge.split(',')))
    parts = set(A.only.split(',')) if A.only else set()
    if 'P' in parts:
        part_P(frames, RAW, A.force)
    if 'H' in parts:
        part_network('hp', eps_list, frames, RAW, A.force)
    if 'Q' in parts:
        part_network('plain', eps_list, frames, RAW, A.force)
    if 'B' in parts:
        part_bridge(frames, RAW, A.force)
    if 'A' in parts:
        part_check([f for f in ACC_FRAMES if f in frames], eps_list, RAW, A.force)
    if 'S' in parts:
        part_stall(frames, RAW, A.force)
    if A.from_raw:
        with open(MD, 'w') as f:
            f.write(report())
        print('wrote', MD)
    if not (parts or A.merge or A.from_raw):
        ap.error('give --only, --merge or --from-raw')


if __name__ == '__main__':
    main()
