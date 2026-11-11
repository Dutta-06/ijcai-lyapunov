"""
Baselines for choosing the number of kernels of a sliding-kernel fit (WP3), on the same inputs as the sliding inverse scale space
====================================================================================================================================

Everything works on the dictionary interface of `rich_transit_nd.BumpsND` (`phi`, `dphi`, `synth`, `corr`, `N`, `n`): data y on the
samples, a start grid of centres, Gaussian kernels of one width. The methods a practitioner would try instead of (or beside) the training
of Paper 17, each with the stopping or selection rule it comes with:

  omp_slide        matching pursuit with sliding: the grid kernel of the largest |⟨rφ⟩| enters, the whole block (coefficients and centres)
                   is refitted by unconstrained least squares — the greedy cousin of the sliding path (same candidates, same block
                   fits; the entry time is the instantaneous correlation instead of its integral) and of ADCG / sliding Frank–Wolfe
  stepwise_bic     forward selection on the fixed grid (coefficients only, centres stay on the grid): ordinary least squares with BIC
  lasso_grid       ℓ1-penalised least squares on the fixed grid (coordinate descent), the penalty by BIC along the path: the convex relaxation
  penalised_net    the network trained from a standard (not small) start with an ℓ1 penalty on the coefficients and trainable centres,
                   the penalty chosen by BIC: "penalised training"
  prune_net        the over-parameterised network (every grid kernel a unit) trained without penalty from a random start of O(1) and pruned
                   by a threshold on |w|: "random initialisation and prune"
  sfw              sliding Frank–Wolfe (Denoyelle, Duval, Peyré, Soubies 2020) for the BLASSO  min ½⟨(y − Σ w_iφ(c_i))²⟩ + λ‖w‖₁: the kernel of the
                   largest |⟨rφ⟩| enters while that exceeds λ, the coefficients are refitted by the lasso on the support and the whole block (coefficients
                   and centres) is slid on the penalised objective; a penalty λ is its one parameter
  cpgd             conic particle gradient descent (Chizat 2022) for the same problem: a particle for every grid kernel and sign, weights updated
                   multiplicatively, positions by gradient steps independent of the weights (the Wasserstein–Fisher–Rao flow), the penalty λ the one
                   parameter, the particles of the result clustered into spikes

and the rules that turn a sequence of nested fits into a number of kernels: `select_discrepancy` (needs the noise level), `select_bic`,
`select_significance` (the entrant's correlation against the noise scale), and for the sliding path `select_gap`. Nothing here is tuned
to make the sliding path look good: the experiment that uses it (`experiments/exp_noise_order.py`) gives each method its own rule.
"""
import numpy as np
from scipy.optimize import least_squares, linear_sum_assignment


# ---------------------------------------------------------------------------------------------------------------------------
# block fits
# ---------------------------------------------------------------------------------------------------------------------------
def refit_block(K, y, w, c, max_nfev=200):
    """Unconstrained least squares of the block (w (k), c (k, n)) from where it is (Levenberg–Marquardt, analytic Jacobian): the fit a matching
    pursuit with sliding performs after each addition. Returns (w, c)."""
    k, n = len(w), K.n
    if k == 0:
        return w, c

    def res(z):
        return (y - K.synth(z[k:].reshape(k, n), z[:k])) / np.sqrt(K.N)

    def jac(z):
        cc = z[k:].reshape(k, n)
        return -np.c_[K.phi(cc), (K.dphi(cc) * z[:k][None, :, None]).reshape(K.N, k * n)] / np.sqrt(K.N)
    sol = least_squares(res, np.r_[w, np.asarray(c, float).ravel()], jac=jac, method='lm', xtol=1e-12, ftol=1e-12, gtol=1e-12, max_nfev=max_nfev)
    return sol.x[:k], sol.x[k:].reshape(k, n)


def coef_lsq(K, y, c):
    """Least-squares coefficients of the kernels at the centres c (k, n)."""
    if len(c) == 0:
        return np.zeros(0)
    return np.linalg.lstsq(K.phi(c), y, rcond=None)[0]


def entry_significance(K, r, c, s):
    """The entrant's correlation against its noise scale: |⟨rφ_c⟩|·N/(s‖φ_c‖) — the number of noise standard deviations (white noise of sd s)."""
    g = K.corr(c, r)[0]
    return np.abs(g) * K.N / (s * np.linalg.norm(K.phi(c), axis=0))


# ---------------------------------------------------------------------------------------------------------------------------
# greedy pursuit with sliding
# ---------------------------------------------------------------------------------------------------------------------------
def omp_slide(K, y, c_grid, max_units=12, slide=True, s=None):
    """Matching pursuit with sliding on the start grid c_grid (d, n). Returns dict(blocks, order, sig): blocks[k] = dict(w, c, rss) after k units
    (blocks[0] the empty fit), order = the grid indices in the order of entry, sig = the entrant's significance against noise level s (if given)."""
    c_grid = np.atleast_2d(np.asarray(c_grid, float))
    if c_grid.shape[0] == 1 and K.n == 1 and c_grid.shape[1] != 1:
        c_grid = c_grid.T
    n = K.n
    c, w = np.zeros((0, n)), np.zeros(0)
    r = np.array(y, float)
    blocks = [dict(w=w, c=c, rss=float(r @ r))]
    order, sig, used = [], [], np.zeros(len(c_grid), bool)
    for _ in range(max_units):
        g = np.abs(K.corr(c_grid, r)[0])
        g[used] = -1.0
        j = int(np.argmax(g))
        if g[j] <= 0:
            break
        if s is not None:
            sig.append(float(entry_significance(K, r, c_grid[j:j + 1], s)[0]))
        used[j] = True
        order.append(j)
        c = np.vstack([c, c_grid[j][None, :]])
        w = coef_lsq(K, y, c)
        if slide:
            w, c = refit_block(K, y, w, c)
        r = y - K.synth(c, w)
        blocks.append(dict(w=w, c=c, rss=float(r @ r)))
    return dict(blocks=blocks, order=order, sig=np.array(sig))


def stepwise_bic(K, y, c_grid, max_units=12):
    """Forward selection of grid kernels by least squares (the centres stay on the grid): the sequence of nested fits."""
    c_grid = np.atleast_2d(np.asarray(c_grid, float))
    if K.n == 1 and c_grid.shape[0] == 1:
        c_grid = c_grid.T
    Phi = K.phi(c_grid)
    chosen, blocks = [], [dict(w=np.zeros(0), c=np.zeros((0, K.n)), rss=float(y @ y))]
    for _ in range(max_units):
        best = None
        for j in range(Phi.shape[1]):
            if j in chosen:
                continue
            idx = chosen + [j]
            wj, rs, *_ = np.linalg.lstsq(Phi[:, idx], y, rcond=None)
            rss = float(np.sum((y - Phi[:, idx] @ wj) ** 2))
            if best is None or rss < best[0]:
                best = (rss, j, wj)
        chosen.append(best[1])
        blocks.append(dict(w=best[2], c=c_grid[chosen], rss=best[0]))
    return dict(blocks=blocks, order=chosen)


# ---------------------------------------------------------------------------------------------------------------------------
# the sliding inverse scale space in its limit form, on the same interface
# ---------------------------------------------------------------------------------------------------------------------------
def iss_path(K, y, c_grid, max_entries=12, rho=1.0, tie_tol=1e-9, vtol=1e-12, s=None):
    """The sliding inverse scale space of Paper 17 in its limit form on the dictionary interface (the polish-only version of
    `rich_transit_nd.sliding_iss_nd`, which reached the same path on the benchmark): the inactive kernels' duals grow at ⟨rφ⟩ until |u_j| = ρ,
    the entrant (a tie together) joins with the sign of its dual, the block is refitted by sign-constrained least squares with sliding centres
    (`polish_block_nd`), a kernel whose coefficient falls to 0 against its dual's push leaves. Returns dict(blocks, taus, order, sig):
    blocks[k] = dict(w, c, rss) after k entries (blocks[0] the empty fit), taus the entry times, order the grid indices, sig the entrants'
    significance against the noise level s (if given)."""
    from src.core.rich_transit_nd import polish_block_nd
    c_grid = np.atleast_2d(np.asarray(c_grid, float))
    if K.n == 1 and c_grid.shape[0] == 1:
        c_grid = c_grid.T
    d, n = c_grid.shape
    c, w, u = c_grid.copy(), np.zeros(d), np.zeros(d)
    sg = np.zeros(d)
    act = np.zeros(d, bool)
    tau, taus, order, sig = 0.0, [], [], []
    blocks = [dict(w=np.zeros(0), c=np.zeros((0, n)), rss=float(y @ y))]
    vscale = None
    for _ in range(max_entries):
        r = y - K.synth(c[act], w[act]) if act.any() else np.array(y, float)
        v = K.corr(c, r)[0]
        if vscale is None:
            vscale = float(np.abs(v).max())
            if vscale <= 0:
                break
        dt = np.full(d, np.inf)
        ina = ~act & (np.abs(v) > vtol * vscale)
        dt[ina] = (rho - np.sign(v[ina]) * u[ina]) / np.abs(v[ina])
        dmin = float(dt.min())
        if not np.isfinite(dmin):
            break
        ent = np.flatnonzero(dt <= dmin * (1 + tie_tol) + 1e-300)
        u[~act] += dmin * v[~act]
        tau += dmin
        for j in ent:
            u[j] = np.sign(v[j]) * rho
            sg[j] = np.sign(v[j])
            act[j] = True
            taus.append(tau)
            order.append(int(j))
            if s is not None:
                sig.append(float(entry_significance(K, r, c[j:j + 1], s)[0]))
        idx = np.flatnonzero(act)
        w0 = coef_lsq(K, y, c[idx])
        w0 = sg[idx] * np.maximum(sg[idx] * w0, 1e-6 * max(np.abs(w0).max(), 1e-12))
        wp, cp = polish_block_nd(K, y, w0, c[idx], sg[idx])
        c[idx], w[idx] = cp, wp
        r = y - K.synth(c[act], w[act])
        g = K.corr(c[idx], r)[0]
        out = (np.abs(wp) <= 1e-9 * max(np.abs(wp).max(), 1e-300)) & (sg[idx] * g <= vtol * vscale)
        for i in idx[out]:
            act[i] = False
            w[i] = 0.0
            c[i] = c_grid[i]
            u[i] = sg[i] * rho
        idx = np.flatnonzero(act)
        r = y - K.synth(c[idx], w[idx]) if idx.size else np.array(y, float)
        blocks.append(dict(w=w[idx].copy(), c=c[idx].copy(), rss=float(r @ r)))
    return dict(blocks=blocks, taus=np.array(taus), order=order, sig=np.array(sig))


# ---------------------------------------------------------------------------------------------------------------------------
# the convex relaxation on the grid
# ---------------------------------------------------------------------------------------------------------------------------
def lasso_cd(Phi, y, lam, w0=None, iters=2000, tol=1e-10):
    """min ½‖y − Φw‖² + λ‖w‖₁ by cyclic coordinate descent (warm start w0)."""
    d = Phi.shape[1]
    w = np.zeros(d) if w0 is None else np.array(w0, float)
    col2 = (Phi ** 2).sum(0)
    r = y - Phi @ w
    for _ in range(iters):
        dmax = 0.0
        for j in range(d):
            if col2[j] == 0:
                continue
            rho = Phi[:, j] @ r + col2[j] * w[j]
            wj = np.sign(rho) * max(abs(rho) - lam, 0.0) / col2[j]
            if wj != w[j]:
                r -= Phi[:, j] * (wj - w[j])
                dmax = max(dmax, abs(wj - w[j]))
                w[j] = wj
        if dmax < tol:
            break
    return w


def lasso_grid(K, y, c_grid, n_lam=40, lam_min_ratio=1e-4):
    """The ℓ1 path on the grid (a geometric ladder of penalties from the smallest that gives w = 0 downwards, warm starts): a list of
    dict(lam, w (d), rss, nnz)."""
    c_grid = np.atleast_2d(np.asarray(c_grid, float))
    if K.n == 1 and c_grid.shape[0] == 1:
        c_grid = c_grid.T
    Phi = K.phi(c_grid)
    lam_max = float(np.abs(Phi.T @ y).max()) * (1.0 + 1e-9)      # the smallest penalty that gives w = 0 (strictly above, against rounding)
    out, w = [], None
    for lam in lam_max * np.geomspace(1.0, lam_min_ratio, n_lam):
        w = lasso_cd(Phi, y, lam, w)
        out.append(dict(lam=float(lam), w=w.copy(), rss=float(np.sum((y - Phi @ w) ** 2)), nnz=int((w != 0).sum())))
    return out


def merge_adjacent(w, c_grid, width):
    """Number of clusters of the nonzero coefficients when atoms closer than `width` (consecutive on the sorted grid) count as one spike."""
    idx = np.flatnonzero(w != 0)
    if idx.size == 0:
        return 0
    pos = np.sort(np.asarray(c_grid, float).ravel()[idx])
    return int(1 + np.sum(np.diff(pos) > width))


# ---------------------------------------------------------------------------------------------------------------------------
# trained networks: penalised, and random start + prune
# ---------------------------------------------------------------------------------------------------------------------------
def train_net(K, y, c0, w0, lam=0.0, steps=4000, lr_w=None, lr_c=None, train_c=True, prox=True):
    """Full-batch gradient descent on ½⟨(y − Σ w_iφ(c_i))²⟩ + λ‖w‖₁ (the ℓ1 term by a proximal step on w, the centres by plain gradient):
    a standard start (not small), the network of the baselines "penalised training" and "random initialisation and prune". The step for w is
    1/λ_max of the Gram matrix at the start unless given, the centres' step is lr_c·(that) (default 0.1)."""
    c = np.array(c0, float).reshape(len(w0), K.n)
    w = np.array(w0, float)
    Phi = K.phi(c)
    L = max(float(np.linalg.eigvalsh(Phi.T @ Phi / K.N).max()), 1e-12)
    lw = (1.0 / L) if lr_w is None else lr_w
    lc = (0.1 * lw) if lr_c is None else lr_c
    for _ in range(steps):
        r = y - K.synth(c, w)
        g = K.corr(c, r)
        gw = -g[0]
        if train_c:
            c = c + lc * (w[:, None] * g[1])
        w = w - lw * gw
        if lam > 0 and prox:
            w = np.sign(w) * np.maximum(np.abs(w) - lw * lam, 0.0)
    return w, c


def penalised_net(K, y, c_grid, lams, steps=3000, seed=0):
    """Penalised training: one network of d units (one per grid kernel, coefficients drawn N(0, 0.1²), centres trainable) per penalty in `lams`.
    Returns a list of dict(lam, w, c, rss, nnz)."""
    c_grid = np.atleast_2d(np.asarray(c_grid, float))
    if K.n == 1 and c_grid.shape[0] == 1:
        c_grid = c_grid.T
    rng = np.random.default_rng(seed)
    w0 = 0.1 * rng.standard_normal(len(c_grid))
    out = []
    for lam in lams:
        w, c = train_net(K, y, c_grid, w0, lam=lam, steps=steps)
        r = y - K.synth(c, w)
        out.append(dict(lam=float(lam), w=w, c=c, rss=float(r @ r), nnz=int((np.abs(w) > 1e-6 * max(np.abs(w).max(), 1e-12)).sum())))
    return out


def prune_net(K, y, c_grid, steps=6000, seed=0, thr=0.1, scale=0.3):
    """Random start of O(1) coefficients (N(0, scale²)) for every grid kernel, trained without penalty, then the units with |w| < thr·max|w|
    removed. Returns dict(w, c, kept (bool d), n_kept, rss)."""
    c_grid = np.atleast_2d(np.asarray(c_grid, float))
    if K.n == 1 and c_grid.shape[0] == 1:
        c_grid = c_grid.T
    rng = np.random.default_rng(seed)
    w, c = train_net(K, y, c_grid, scale * rng.standard_normal(len(c_grid)), lam=0.0, steps=steps)
    kept = np.abs(w) >= thr * np.abs(w).max()
    r = y - K.synth(c[kept], w[kept])
    return dict(w=w, c=c, kept=kept, n_kept=int(kept.sum()), rss=float(r @ r))


# ---------------------------------------------------------------------------------------------------------------------------
# the BLASSO by greedy and by particle methods, a penalty λ their parameter
# ---------------------------------------------------------------------------------------------------------------------------
def _slide_l1(K, y, w, c, lam, lo, hi, max_iter=300):
    """Local minimisation of J(w, c) = ½⟨(y − Σ w_iφ(c_i))²⟩ + λ‖w‖₁ over the coefficients and the centres (c within [lo, hi]), the signs of w kept
    (so the penalty is linear in |w|): L-BFGS-B on (|w|, c) with the analytic gradient. Returns (w, c)."""
    from scipy.optimize import minimize
    k, n = len(w), K.n
    s = np.sign(w)

    def fun(z):
        b, cc = z[:k], z[k:].reshape(k, n)
        r = y - K.synth(cc, s * b)
        g0, g1 = K.corr(cc, r)
        return 0.5 * float(r @ r) / K.N + lam * float(b.sum()), np.r_[lam - s * g0, (-(s * b)[:, None] * g1).ravel()]

    bounds = [(0.0, None)] * k + [(float(lo[j]), float(hi[j])) for _ in range(k) for j in range(n)]
    sol = minimize(fun, np.r_[np.abs(w), np.asarray(c, float).ravel()], jac=True, method='L-BFGS-B', bounds=bounds,
                   options=dict(maxiter=max_iter, ftol=1e-16, gtol=1e-11, maxcor=20))
    return s * sol.x[:k], sol.x[k:].reshape(k, n)


def sfw(K, y, c_grid, lam, max_units=12, max_iter=60, tol=1e-3):
    """Sliding Frank–Wolfe for min ½⟨(y − Σ w_iφ(c_i))²⟩ + λ‖w‖₁ with the candidate positions the start grid c_grid: while the certificate
    max_grid |⟨rφ⟩|/λ exceeds 1 + tol the grid kernel that attains it joins the support, the coefficients are refitted by the lasso on the support
    (centres fixed), the block (coefficients and centres) is slid on the penalised objective and the kernels whose coefficient fell to zero are dropped.
    Stops when the certificate is within 1 + tol (a solution of the BLASSO up to the grid of candidates), at max_units kernels or at max_iter additions.
    Returns dict(w, c, rss, obj, cert (the final certificate), n_iter, converged)."""
    c_grid = np.atleast_2d(np.asarray(c_grid, float))
    if K.n == 1 and c_grid.shape[0] == 1:
        c_grid = c_grid.T
    lo, hi = c_grid.min(0), c_grid.max(0)
    c, w = np.zeros((0, K.n)), np.zeros(0)
    converged, it = False, 0
    for it in range(max_iter):
        r = y - K.synth(c, w) if len(w) else np.array(y, float)
        g = np.abs(K.corr(c_grid, r)[0])
        j = int(np.argmax(g))
        if g[j] <= lam * (1.0 + tol):
            converged = True
            break
        if len(w) >= max_units:
            break
        c = np.vstack([c, c_grid[j][None, :]])
        w = lasso_cd(K.phi(c), y, lam * K.N, np.r_[w, 0.0])
        keep = w != 0
        if keep.any():
            w, c = _slide_l1(K, y, w[keep], c[keep], lam, lo, hi)
            keep = np.abs(w) > 1e-9 * max(np.abs(w).max(), 1e-300)
            w, c = w[keep], c[keep]
        else:
            w, c = w[keep], c[keep]
    r = y - K.synth(c, w) if len(w) else np.array(y, float)
    cert = float(np.abs(K.corr(c_grid, r)[0]).max() / lam)
    return dict(w=w, c=c, rss=float(r @ r), obj=0.5 * float(r @ r) / K.N + lam * float(np.abs(w).sum()), cert=cert, n_iter=it, converged=converged)


def threshold_units(K, y, res, lam, thr=0.1):
    """A BLASSO solution `res` (dict with w, c, as `sfw` returns) without the units of |w| below thr times the largest, rss and objective recomputed: the
    post-processing `cpgd` applies to its spikes, so that the solutions of the two methods are read alike (the weak spikes of a small penalty fit noise)."""
    w, c = res['w'], res['c']
    if len(w):
        keep = np.abs(w) >= thr * np.abs(w).max()
        w, c = w[keep], c[keep]
    r = y - K.synth(c, w) if len(w) else np.array(y, float)
    return dict(res, w=w, c=c, rss=float(r @ r), obj=0.5 * float(r @ r) / K.N + lam * float(np.abs(w).sum()))


def _cluster(x, s, radius):
    """Labels of the single-linkage clusters of the particles at x (m, n) with the signs s (m): two particles of one sign closer than radius share one."""
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import connected_components
    d = np.linalg.norm(x[:, None, :] - x[None, :, :], axis=2)
    adj = (d <= radius) & (s[:, None] == s[None, :])
    return connected_components(csr_matrix(adj), directed=False)[1]


def cpgd(K, y, c_grid, lam, steps=3000, w_init=1e-3, link=None, thr=0.1):
    """Conic particle gradient descent for min ½⟨(y − Σ w_iφ(c_i))²⟩ + λ‖w‖₁ (Chizat 2022) with signed weights by the doubling trick: one particle of
    weight w_i ≥ 0 and sign +1 and one of sign −1 on every grid kernel, w_init each. With r the residual of Σ s_iw_iφ(x_i), V_i = λ − s_i⟨rφ(x_i)⟩ the
    potential of particle i: w_i ← w_i exp(−η V_i) (the multiplicative, conic update) and x_i ← x_i + η_x s_i⟨r∇φ(x_i)⟩ (the position step does not scale
    with the weight: the Wasserstein–Fisher–Rao flow), η = 1/(2 max|⟨yφ⟩|) and η_x = 0.2 σ² η (σ the kernel width), both halved while a step would
    raise the objective. After `steps` steps the particles of weight above 10 w_init and of one sign closer than `link` (default 1.5 times the grid spacing) are one
    spike (mass: their total, position: their mass-weighted mean) and spikes of mass below thr times the largest are dropped.
    Returns dict(w, c, n_kept, rss, obj, w_part, x_part, s_part)."""
    c_grid = np.atleast_2d(np.asarray(c_grid, float))
    if K.n == 1 and c_grid.shape[0] == 1:
        c_grid = c_grid.T
    d, n = c_grid.shape
    lo, hi = c_grid.min(0), c_grid.max(0)
    x = np.vstack([c_grid, c_grid])
    s = np.r_[np.ones(d), -np.ones(d)]
    w = np.full(2 * d, float(w_init))
    G = max(float(np.abs(K.corr(c_grid, y)[0]).max()), 1e-300)
    eta = 0.5 / G
    ratio_x = 0.2 * getattr(K, 'sig', 1.0) ** 2                       # η_x = 0.2 σ² η: a position step of about a tenth of a width at the start
    r = y - K.synth(x, s * w)
    obj = 0.5 * float(r @ r) / K.N + lam * float(w.sum())
    for _ in range(steps):
        g0, g1 = K.corr(x, r)
        V = lam - s * g0
        for _try in range(30):
            w_new = np.maximum(w * np.exp(np.clip(-eta * V, -30.0, 30.0)), 1e-300)
            x_new = np.clip(x + eta * ratio_x * (s[:, None] * g1), lo, hi)
            r_new = y - K.synth(x_new, s * w_new)
            obj_new = 0.5 * float(r_new @ r_new) / K.N + lam * float(w_new.sum())
            if obj_new <= obj:
                break
            eta *= 0.5
        else:
            break
        w, x, r, obj = w_new, x_new, r_new, obj_new
        eta = min(eta * 1.05, 0.5 / G)
    if link is None:
        dd = np.linalg.norm(c_grid[:, None, :] - c_grid[None, :, :], axis=2)
        link = 1.5 * float(dd[dd > 0].min()) if (dd > 0).any() else 1.0
    live = w > 10.0 * w_init                                           # a particle that has not grown to ten times its start is not part of a spike
    xl, sl, wl = x[live], s[live], w[live]
    lab = _cluster(xl, sl, link) if live.any() else np.zeros(0, int)       # only the live particles are clustered: dormant ones must not bridge two spikes
    W, C = [], []
    for k in np.unique(lab):
        m = lab == k
        W.append(sl[m][0] * wl[m].sum())
        C.append((xl[m] * wl[m][:, None]).sum(0) / wl[m].sum())
    W, C = np.array(W), np.array(C).reshape(-1, n)
    if len(W):
        big = np.abs(W) >= thr * np.abs(W).max()
        W, C = W[big], C[big]
    rr = y - K.synth(C, W) if len(W) else np.array(y, float)
    return dict(w=W, c=C, n_kept=len(W), rss=float(rr @ rr), obj=0.5 * float(rr @ rr) / K.N + lam * float(np.abs(W).sum()), w_part=w, x_part=x, s_part=s)


# ---------------------------------------------------------------------------------------------------------------------------
# rules
# ---------------------------------------------------------------------------------------------------------------------------
def bic_value(rss, N, n_params):
    return N * np.log(max(rss, 1e-300) / N) + n_params * np.log(N)


def select_bic(blocks, N, per_unit):
    """The block of the smallest BIC, N ln(rss/N) + per_unit·k ln N (per_unit = 2 when centres are fitted too, 1 on the grid). Returns k."""
    vals = [bic_value(b['rss'], N, per_unit * len(b['w'])) for b in blocks]
    return int(np.argmin(vals))


def select_discrepancy(blocks, N, s, z=2.0):
    """The first block whose residual mean square is within the noise: rss/N ≤ s²(1 + z√(2/N)) (the discrepancy principle; needs s)."""
    thr = s * s * (1.0 + z * np.sqrt(2.0 / N))
    for k, b in enumerate(blocks):
        if b['rss'] / N <= thr:
            return k
    return len(blocks) - 1


def select_significance(sig, kappa=4.5):
    """The number of leading entrants whose significance is ≥ kappa."""
    bad = np.flatnonzero(np.asarray(sig) < kappa)
    return int(bad[0]) if bad.size else len(sig)


def select_gap(taus, k_min=1):
    """The largest-relative-gap rule on the entry times τ_1 ≤ τ_2 ≤ …: K̂ = argmax_{k ≥ k_min} τ_{k+1}/τ_k (0 if there are fewer than k_min + 1)."""
    taus = np.asarray(taus, float)
    if taus.size <= k_min:
        return 0
    return int(np.argmax(taus[k_min:] / np.maximum(taus[k_min - 1:-1], 1e-300))) + k_min


# ---------------------------------------------------------------------------------------------------------------------------
# scoring a recovered set against the truth
# ---------------------------------------------------------------------------------------------------------------------------
def match_units(c_est, w_est, c_true, w_true, tol):
    """Hungarian pairing of estimated and true kernels (centres in the kernel space, within tol): returns (n_matched, position rms over the
    pairs, relative amplitude rms over the pairs); nan where nothing matched."""
    c_est, c_true = np.atleast_2d(np.asarray(c_est, float)), np.atleast_2d(np.asarray(c_true, float))
    if len(w_est) == 0 or len(w_true) == 0:
        return 0, float('nan'), float('nan')
    if c_est.shape[0] != len(w_est):
        c_est = c_est.T
    if c_true.shape[0] != len(w_true):
        c_true = c_true.T
    d = np.linalg.norm(c_est[:, None, :] - c_true[None, :, :], axis=2)
    big = 1e6
    rr, cc = linear_sum_assignment(np.where(d <= tol, d, big))
    ok = d[rr, cc] <= tol
    if not ok.any():
        return 0, float('nan'), float('nan')
    rr, cc = rr[ok], cc[ok]
    pos = float(np.sqrt(np.mean(d[rr, cc] ** 2)))
    amp = float(np.sqrt(np.mean(((np.asarray(w_est)[rr] - np.asarray(w_true)[cc]) / np.asarray(w_true)[cc]) ** 2)))
    return int(ok.sum()), pos, amp
