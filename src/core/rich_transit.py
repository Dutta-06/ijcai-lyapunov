"""
Proposition D‴.48 (PARADIGM.md §5d): the rich transit — from small masses the kernels enter one at a time
====================================================================================================================================

The free mass–amplitude law of Proposition D‴.45 (dα_i = Gβ_ia′_ig_i, dβ_i = γβ_ia_ig_i, w_i = β_ia_i/m, Q_i conserved) started at α = 0
with masses β⁰_i = εb_i follows, in the time τ = t/Λ(ε), Λ = (πa_max/(4G))(2/(cε))^{1/2}, the inverse scale space flow of
Σ_i|w_i|/√b_i (`mirror_free.iss_path`): the kernels enter one at a time. With the kernels' centres moving by the paradigm's law
(Proposition D‴.34's positions, Euclidean: ċ_i = Gβ_ia_i⟨r∂_cφ_i⟩ = Gmw_i⟨r∂_cφ_i⟩) an inactive kernel — coefficient O(ε^{3/2}) — does
not move, and an entered one slides: the limit is the sliding inverse scale space (`sliding_iss`): the inactive kernels' dual coordinates
integrate the residual's correlation at their starting positions, an entrant joins the active block, the block relaxes on the fast time
(on the critical curve Q = −c, an ε-independent flow) to a stationary point of the fit over its coefficients and centres. With a
Wasserstein–Fisher–Rao metric on the positions instead (the mass dropped from the centre's speed, ċ_i = Ga_i⟨r∂_cφ_i⟩) an inactive kernel's
amplitude is O(ε^{1/2}) and its centre moves at O(1) in τ — dc/dτ = (πa_max²/(4cρ))tan(πu/(2ρ))⟨r∂_cφ⟩, du/dτ = ⟨rφ⟩ (`wfr_inactive_rhs`):
the kernels climb the correlation's peaks before entering — condensation.

This module: Gaussian kernels on a sample grid (`Bumps`; on a reference orbit with the closed loop's error filter, `OrbitBumps` —
Proposition D‴.49's leading order —, and the closed loop's exact averaged drive along the frozen periodic orbit, `FrozenOrbit` with
`frozen_orbit_rhs`), the free law with moving centres under either metric (`moving_rhs`) and with moving widths too
(`moving_rhs_widths`), the active block's fast flow on the critical curve (`critical_rhs`), the sliding inverse scale space
(`sliding_iss`) and the inactive kernels' limit under the Wasserstein–Fisher–Rao metric (`wfr_inactive_rhs`). Proposition D‴.51 (the
transit at a finite mass): the critical curve's inner profile near the wall (`eta_crit`, `coef_crit`, `inner_coef`, `inner_xi`), the
entrant's overshoot (`overshoot`), the ε-free inner problem of a near tie (`near_tie_inner`), the sliding path's duals (`path_duals`), the
creep of a kernel left near its wall along the path (`creep_reaches`) and the mass below which the loop forms exactly the path's kernels,
predicted from the path alone (`recovery_threshold`).
"""
import numpy as np
from scipy.integrate import solve_ivp


class Bumps:
    """Gaussian kernels exp(−(x − c)²/(2σ²)) of one width on the samples xs; ⟨f g⟩ is the mean over the samples."""

    def __init__(self, xs, sig):
        self.xs, self.sig = np.asarray(xs, float), float(sig)

    def phi(self, c):
        return np.exp(-(self.xs[:, None] - np.atleast_1d(c)[None, :]) ** 2 / (2 * self.sig ** 2))

    def dphi(self, c):
        """∂_cφ."""
        c = np.atleast_1d(c)
        return self.phi(c) * (self.xs[:, None] - c[None, :]) / self.sig ** 2

    def mean(self, v):
        return np.asarray(v).sum(axis=0) / self.xs.size

    def filt(self, r):
        """The residual as the drive sees it: itself on a dictionary (`OrbitBumps`: the closed loop's error filter)."""
        return r


class OrbitBumps(Bumps):
    """The closed loop's averaged drive (Proposition D‴.49): the samples x_m(t_s) of a T-periodic reference over one period (M points) and
    the error filter L = 1/(s + k_e) — the periodic steady state of ė = −k_ee + r, by FFT —, so that ⟨φ_i filt(r)⟩ = T⁻¹∮φ_i(x_m)L[r(x_m)]dt
    is the averaged drive ⟨eφ_i(x)⟩ of the slow adaptation."""

    def __init__(self, xm, T, M, sig, ke):
        ts = T * np.arange(M) / M
        super().__init__(xm(ts), sig)
        om = 2 * np.pi * np.fft.fftfreq(M, d=T / M)
        self.H, self.ke, self.T = 1.0 / (1j * om + ke), float(ke), float(T)

    def filt(self, r):
        return np.real(np.fft.ifft(np.fft.fft(r) * self.H))


class FrozenOrbit:
    """The closed loop's exact averaged drive (Proposition D‴.49): with the parameters frozen, the loop ė = −k_ee + r(x_m + e),
    r = f − f̂, has for k_e > Lip(r) a unique, exponentially stable T-periodic solution e_p = L[r(x_m + e_p)] (L = 1/(s + k_e); the map is
    a contraction of the periodic functions with factor Lip(r)/k_e) — computed on M samples of a period by the fixed point through the FFT,
    warm-started from the last call. The slow adaptation averages the regressor along the frozen orbit x_p = x_m + e_p: g_i = ⟨φ_i(x_p)e_p⟩,
    h_i = ⟨∂_cφ_i(x_p)e_p⟩ (the period's mean). `OrbitBumps` is its truncation at x_p = x_m, e_p = L[r(x_m)]."""

    def __init__(self, xm, T, M, sig, ke, f, tol=1e-12, it_max=2000):
        self.ts = T * np.arange(M) / M
        self.xm = xm(self.ts)
        om = 2 * np.pi * np.fft.fftfreq(M, d=T / M)
        self.H = 1.0 / (1j * om + ke)
        self.sig, self.ke, self.f, self.tol, self.it_max = float(sig), float(ke), f, tol, it_max
        self.e = np.zeros(M)
        self.iters = 0

    def orbit(self, c, w):
        """The frozen orbit's error e_p (raises if the fixed point does not converge — k_e too small against Lip(r))."""
        e = self.e
        for k in range(self.it_max):
            x = self.xm + e
            ph = np.exp(-(x[:, None] - c[None, :]) ** 2 / (2 * self.sig ** 2))
            en = np.real(np.fft.ifft(np.fft.fft(self.f(x) - ph @ w) * self.H))
            dif = np.abs(en - e).max()
            e = en
            if dif <= self.tol * max(1.0, np.abs(e).max()):
                break
        else:
            raise RuntimeError('FrozenOrbit: the periodic orbit did not converge (k_e ≤ Lip(f − f̂)?)')
        self.iters += k + 1
        self.e = e
        return e

    def drive(self, c, w):
        """(g, h, e_p): the averaged correlations ⟨φ_i(x_p)e_p⟩, ⟨∂_cφ_i(x_p)e_p⟩ and the orbit's error."""
        e = self.orbit(c, w)
        x = self.xm + e
        dx = x[:, None] - c[None, :]
        ph = np.exp(-dx ** 2 / (2 * self.sig ** 2))
        return ph.T @ e / e.size, (ph * dx / self.sig ** 2).T @ e / e.size, e


def frozen_orbit_rhs(L, O, b0, metric='paradigm'):
    """The averaged closed loop on the free law's curves from α = 0 with the masses b0 (β = b0 + 2c sinh²(α/a_max), Q kept), z = (α, c):
    α̇ = Gβa′g, ċ = Gβa h (metric 'paradigm'), Ga h ('wfr') or 0 ('fixed'); (g, h) the frozen orbit's drive (`FrozenOrbit`)."""
    if metric not in ('paradigm', 'wfr', 'fixed'):
        raise ValueError(metric)
    b0 = np.asarray(b0, float)

    def f(t, z):
        d = z.size // 2
        al, c = z[:d], z[d:]
        b = b0 + 2 * L.c * np.sinh(al / L.A) ** 2
        a = L.A * np.tanh(al / L.A)
        g, h, _ = O.drive(c, b * a / L.m)
        speed = L.G * b * a if metric == 'paradigm' else (L.G * a if metric == 'wfr' else 0.0 * a)
        return np.r_[curve_mobility(L, al, b0) * g, speed * h]
    return f


def transit_scale(L, eps):
    """Λ(ε) = (πa_max/(4G))(2/(cε))^{1/2}: λ₁ = Λ·b^{−1/2}(1 + O(ε)) for the mass εb (Proposition D‴.45(ii))."""
    return (np.pi * L.A / (4 * L.G)) * np.sqrt(2 / (L.c * eps))


def moving_rhs(L, K, fstar, metric='paradigm'):
    """The free law with moving centres, z = (α, log β, c): α̇ = Gβa′g, (log β)˙ = γag, ċ = Gβa⟨r∂_cφ⟩ (metric 'paradigm') or Ga⟨r∂_cφ⟩
    ('wfr'); g = ⟨rφ⟩, r = f* − Σ_i(β_ia_i/m)φ_i."""
    if metric not in ('paradigm', 'wfr'):
        raise ValueError(metric)

    def f(t, z):
        d = z.size // 3
        al, lb, c = z[:d], z[d:2 * d], z[2 * d:]
        b = np.exp(lb); a, da = L.amp(al)
        Ph = K.phi(c)
        r = K.filt(fstar - Ph @ (b * a / L.m))
        g = Ph.T @ r / K.xs.size
        h = K.dphi(c).T @ r / K.xs.size
        speed = L.G * b * a if metric == 'paradigm' else L.G * a
        return np.r_[L.G * b * da * g, L.gamma * a * g, speed * h]
    return f


def moving_rhs_widths(L, K, fstar):
    """The free law with moving centres and widths (the paradigm's metric), z = (α, log β, c, σ): ċ = Gβa⟨r∂_cφ⟩, σ̇ = Gβa⟨r∂_σφ⟩,
    φ_i = exp(−(x − c_i)²/(2σ_i²)) on K's samples (K's own width unused)."""
    def f(t, z):
        d = z.size // 4
        al, lb, c, sg = z[:d], z[d:2 * d], z[2 * d:3 * d], z[3 * d:]
        b = np.exp(lb); a, da = L.amp(al)
        dx = K.xs[:, None] - c[None, :]
        Ph = np.exp(-dx ** 2 / (2 * sg ** 2))
        r = K.filt(fstar - Ph @ (b * a / L.m))
        g = Ph.T @ r / K.xs.size
        h = (Ph * dx / sg ** 2).T @ r / K.xs.size
        hs = (Ph * dx ** 2 / sg ** 3).T @ r / K.xs.size
        sp = L.G * b * a
        return np.r_[L.G * b * da * g, L.gamma * a * g, sp * h, sp * hs]
    return f


def critical_w(L, al, b0=0.0):
    """The coefficient on the curve from α = 0 with the mass b0 (b0 = 0: the critical curve Q = −c): β = b0 + 2c sinh²(α/a_max) — Q kept
    exactly, without cancellation —, w = βa/m (≈ 2cα³/(ma_max²) at small α on the critical curve)."""
    return (b0 + 2 * L.c * np.sinh(al / L.A) ** 2) * L.A * np.tanh(al / L.A) / L.m


def curve_mobility(L, al, b0=0.0):
    """dα/dη = Gβa′ on the curve from α = 0 with the mass b0: G(b0 sech²(α/a_max) + 2c tanh²(α/a_max)), no cancellation at large α."""
    x = np.minimum(np.abs(al) / L.A, 350.0)
    return L.G * (b0 / np.cosh(x) ** 2 + 2 * L.c * np.tanh(x) ** 2)


def critical_rhs(L, K, fstar, metric='paradigm', b0=0.0):
    """The free law with moving centres on its curves, z = (α, c): α̇ = Gβa′g, ċ = Gβa⟨r∂_cφ⟩ (metric 'paradigm') or Ga⟨r∂_cφ⟩ ('wfr'),
    β = b0 + 2c sinh²(α/a_max) — each kernel kept on its curve exactly (round 37: integrating (α, log β) drifts off Q at small masses).
    b0 = 0: the active block's fast flow on the critical curve; b0 = εb: the whole loop from α = 0 with the masses εb."""
    if metric not in ('paradigm', 'wfr'):
        raise ValueError(metric)

    def f(t, z):
        d = z.size // 2
        al, c = z[:d], z[d:]
        b = b0 + 2 * L.c * np.sinh(al / L.A) ** 2
        a = L.A * np.tanh(al / L.A)
        Ph = K.phi(c)
        r = K.filt(fstar - Ph @ (b * a / L.m))
        g = Ph.T @ r / K.xs.size
        h = K.dphi(c).T @ r / K.xs.size
        speed = L.G * b * a if metric == 'paradigm' else L.G * a
        return np.r_[curve_mobility(L, al, b0) * g, speed * h]
    return f


def polish_block(K, fstar, w, c, s, cbox=None, max_nfev=5000, info=False, gtol=1e-15):
    """The active block's fit from (w, c) by a sign-constrained least-squares solve (trust-region reflective): min ½⟨(f* − Φ(c)w)²⟩ over
    (w, c) with s_iw_i ≥ 0. On a dictionary (K.filt the identity) the fast flow is a gradient flow of this fit, so its end is a stationary
    point of it; the solve reaches the point where the flow crawls along the block's slow modes (two kernels Δ apart: the fit's Hessian has
    eigenvalues ∝ 1, 1, Δ², Δ⁶ — Proposition D‴.52). cbox: (lo, hi), the centres kept in [lo, hi] (on noisy data the fit can have no
    minimiser: a kernel leaving the samples' window fits the noise at the edge with a coefficient growing without bound — referee round
    45). max_nfev: the solve's budget of evaluations (its end where the budget binds is not stationary: referee round 54). gtol: the
    solve's absolute gradient test (TRF's ‖g·v‖∞ < gtol, g·v the gradient scaled to the bounds); None disables it — at a fit error of
    10⁻⁹–10⁻¹⁰ the gradient along the block's slow modes is 10⁻¹⁶–10⁻¹⁷, and 10⁻¹⁵ stops the solve above the optimum (referee round 60).
    Returns (w, c), and with info=True also dict(status, nfev) of the solve (status 0: the budget spent; 1: the gradient test)."""
    from scipy.optimize import least_squares
    n, nx = w.size, K.xs.size
    clo, chi = (-np.inf, np.inf) if cbox is None else cbox
    lo = np.r_[np.where(s > 0, 0.0, -np.inf), np.full(n, clo)]
    hi = np.r_[np.where(s > 0, np.inf, 0.0), np.full(n, chi)]
    z0 = np.r_[np.where(s > 0, np.maximum(w, 1e-300), np.minimum(w, -1e-300)), c]
    sol = least_squares(lambda z: (fstar - K.phi(z[n:]) @ z[:n]) / np.sqrt(nx), z0,
                        jac=lambda z: -np.c_[K.phi(z[n:]), K.dphi(z[n:]) * z[:n][None, :]] / np.sqrt(nx),
                        bounds=(lo, hi), method='trf', xtol=1e-15, ftol=1e-15, gtol=gtol, max_nfev=max_nfev)
    if info:
        return sol.x[:n], sol.x[n:], dict(status=int(sol.status), nfev=int(sol.nfev))
    return sol.x[:n], sol.x[n:]


def ascent_target(K, r, c, s, step=None, max_steps=4000):
    """Where a kernel with a near-zero coefficient of sign s at c is carried in the residual r: while its coefficient is small r barely
    changes and its centre climbs s⟨rφ(c)⟩ (ċ ∝ w⟨r∂_cφ⟩, w of sign s), so the end of that ascent — the local maximum of s⟨rφ(c)⟩ it
    reaches. Plain gradient ascent with steps of σ/40 in the gradient's direction, stopped where the gradient changes direction."""
    c = np.array(c, float).reshape(-1)
    step = K.sig / 40 if step is None else step
    nd = c.size > 1 or np.ndim(K.xs) != 1
    nx = K.N if hasattr(K, 'N') else K.xs.size
    grad = lambda cc: s * (K.dphi(cc[None, :] if nd else cc).reshape(nx, -1).T @ r) / nx
    g = grad(c)
    for _ in range(max_steps):
        gn = float(np.linalg.norm(g))
        if gn == 0:
            break
        c_new = c + step * g / gn
        g2 = grad(c_new)
        if float(np.dot(g, g2)) <= 0:
            break
        c, g = c_new, g2
    return c


STIFF_BUDGET = 200000                       # right-hand-side evaluations of one LSODA chunk before `relax_chunk` redoes it by Radau
RELAX_SPLITS = 6                            # a chunk whose solve fails even by Radau is halved, to this depth


class _Stalled(Exception):
    pass


class _OutOfBudget(Exception):
    pass


def relax_chunk(rhs, z, chunk, jac=None, rtol=1e-10, atol=1e-13, budget=None, _depth=0):
    """One chunk of a block's fast flow from z over (0, chunk): LSODA, unless it spends more than `budget` (STIFF_BUDGET) right-hand-side
    evaluations — its stiffness detection can fail on a block whose coefficients span many decades (fast modes at −44 and slow ones at
    10⁻¹², integrated by the non-stiff method: 6.7·10⁶ evaluations against Radau's 100) —, then the chunk is redone from z by Radau (the
    exact Jacobian when given). A solve that fails (LSODA's "unexpected istate" from relaxed states on chunks beyond ~10¹¹, its last state
    the start — returned before as if the chunk had run: referee round 48) is also redone by Radau, and if Radau fails too, as two half
    chunks, to the depth RELAX_SPLITS. Below the budget a successful LSODA call is the one `sliding_iss` always made, so paths that never
    stall or fail are unchanged. Returns (the end state, how: 0 LSODA, 1 redone by Radau after the stall, 2 redone after a failed solve,
    3 failed even halved — the state the last successful piece reached)."""
    budget = STIFF_BUDGET if budget is None else budget
    kw = {} if jac is None else dict(jac=jac)
    n = [0]

    def counted(s, y):
        n[0] += 1
        if n[0] > budget:
            raise _Stalled
        return rhs(s, y)
    try:
        sol = solve_ivp(counted, (0, chunk), z, method='LSODA', rtol=rtol, atol=atol, **kw)
        if sol.success:
            return sol.y[:, -1], 0
        how = 2
    except _Stalled:
        how = 1
    sol = solve_ivp(rhs, (0, chunk), z, method='Radau', rtol=rtol, atol=atol, **kw)
    if sol.success:
        return sol.y[:, -1], how
    if _depth >= RELAX_SPLITS:
        return z, 3
    z1, h1 = relax_chunk(rhs, z, chunk / 2, jac, rtol, atol, budget, _depth + 1)
    if h1 == 3:
        return z1, 3
    z2, h2 = relax_chunk(rhs, z1, chunk / 2, jac, rtol, atol, budget, _depth + 1)
    return z2, 3 if h2 == 3 else 2


def leaves(si, g, wi, big, own=None):
    """The block's kernels that leave: a small coefficient (not `big`) and a correlation pointing into its box (s_i⟨rφ_i⟩ < 0) — and, with
    own (⟨φ_iφ_i⟩, the exit_kkt rule), still inward with the kernel's own contribution removed, s_i(⟨rφ_i⟩ + w_i⟨φ_iφ_i⟩) < 0."""
    out = ~big & (si * g < 0)
    return out if own is None else out & (si * (g + wi * own) < 0)


def sliding_iss(L, K, fstar, c0, rho, seed=1e-7, t_fast=1e5, tol=1e-9, exit_tol=1e-6, max_events=60, tie_tol=1e-9, gtol=1e-9,
                max_doublings=14, polish=False, max_entries=None, polish_step=1.0, polish_wfrac=1e-2, box=None, w_scale=None,
                max_rhs=None, exit_kkt=False, polish_nfev=5000, repolish=0, polish_gtol=1e-15):
    """The sliding inverse scale space (the paradigm's metric). The inactive kernels stay where they are (their start, or where they
    left the block), their u_j growing at ⟨rφ(c_j)⟩ (r the active block's residual); the entrants — every kernel whose |u_j| reaches ρ_j
    within the relative tie_tol of the first, so that exact ties enter together — join the active block with the coefficient `seed` (the
    sign of u_j); the block relaxes by `critical_rhs` (the fast flow on the critical curve) over t_fast, then over the elapsed time again
    (doubling it) until it is stationary — on every kernel of the block with a coefficient above exit_tol·max, |⟨rφ_i⟩| and
    |w_i⟨r∂_cφ_i⟩| below gtol·max_j|⟨f*φ_j⟩|; on the others, no outward correlation (s_i⟨rφ_i⟩ ≤ gtol·max) — or max_doublings is reached.
    A kernel of the block whose coefficient has fallen below exit_tol·max leaves, frozen where it is, only if its correlation points into
    its box (s_i⟨rφ_i⟩ < 0); one with an outward correlation stays on its wall. Stops when the rms residual is below tol. Returns
    dict(events = [(τ, j, 'enter'/'exit', c_j)], hist = [(τ, c, w, rms)], c, w, act, relaxed = [whether each relaxation became
    stationary]). polish=True (a dictionary only): after the first t_fast the block is moved to its fit's stationary point by
    `polish_block` — the flow's end when no other stationary point lies on its way — and checked there, the flow resuming if needed.
    The polish is accepted only if every kernel whose coefficient is still below polish_wfrac·max|w| (near its seed) lands within
    polish_step kernel widths of where the flow would carry it — its start, or the end of its correlation's ascent (`ascent_target`):
    such a kernel's centre is undetermined by the fit, and the solve can carry it into another feature's basin (referee round 43); a kernel
    with a coefficient moves where the fit takes it. Otherwise the flow goes on and the polish is tried again after each doubling. The
    count of rejected polishes is returned as 'polish_rejected', that of the flow's chunks redone by Radau (`relax_chunk`) as 'stiff_redone'.
    max_entries: stop as soon as that many kernels have entered, before the block relaxes (the path's events up to that entry).
    Past the first kernel that fits noise the fit can have no minimiser — a kernel leaves the samples' window and its coefficient grows
    without bound, fitting the noise at the edge — and thresholds relative to max|w| then take the features' kernels for small ones
    (referee round 45). box: (lo, hi) — the centres stay in [lo, hi] in the polish, and a relaxation that carries one outside stops the
    path ('runaway' True; that block's state is returned as c, w but is not a fit of hist). w_scale: the exits', the stationarity check's
    and the polish's thresholds as multiples of w_scale (the data's scale) instead of max|w|. max_rhs: a budget of right-hand-side
    evaluations for the whole path (deterministic, unlike a wall-clock limit); spent, the path stops before the block that spent it is
    relaxed ('out_of_budget' True). 'ev_iter' gives each event's iteration i: hist[i] is the state before the iteration's entries and
    hist[i + 1] the block relaxed after them (ties enter in one iteration; an exit and a re-entry at the same τ in two). exit_kkt: a
    kernel leaves only if its correlation stays inward with its own contribution w_i⟨φ_iφ_i⟩ removed — otherwise the removal turns it
    outward and the kernel re-enters at the same τ, a cycle of entries and exits (referee round 48); False keeps the paths of the results
    committed before. The chunks redone after a failed solve and those that failed even halved (`relax_chunk`) are counted as
    'solver_redone' and 'solver_failed' (a failed chunk ends its relaxation unrelaxed, and no kernel exits on that block: referee round
    54). polish_nfev: the polish's budget of evaluations, polish_gtol its gradient test (`polish_block`'s gtol; None: none, referee
    round 60); repolish: how many more polishes a block may take once it passes the
    stationarity test — each accepted if it lowers the fit's error by more than a relative 1e-9 and passes the seed rule, the flow then
    resuming from it (the flow alone crawls along a block's slow modes, and a single polish can spend its budget before the optimum:
    referee round 54). 'polish_log' lists per iteration its polishes (status, nfev, accepted; e0, e1 the rms before and after for the
    repeated ones) and 'entry_info' per entrant its dual and rate before the entry: u0, the rate ⟨rφ_j⟩, the largest inactive rate, its
    time to the wall and the least (the entrant is the kernel of the least time to the wall, not of the largest rate)."""
    d = len(c0)
    rho = np.broadcast_to(np.asarray(rho, float), (d,))
    u = np.zeros(d); c = np.array(c0, float); al = np.zeros(d)
    act = np.zeros(d, bool)
    tau, ev, hist, relaxed, n_rej, n_stiff, n_redo, n_fail = 0.0, [], [], [], 0, 0, 0, 0
    polish_log, entry_info = [], []
    nx = K.xs.size
    scale = float(np.abs(K.phi(c).T @ K.filt(fstar)).max()) / nx
    rhs = critical_rhs(L, K, fstar)
    nfev = [0]
    if max_rhs is not None:
        rhs_ = rhs

        def rhs(s_, y_):
            nfev[0] += 1
            if nfev[0] > max_rhs:
                raise _OutOfBudget
            return rhs_(s_, y_)
    runaway = out_of_budget = False
    ev_iter = []
    if polish and type(K).filt.__qualname__ != 'Bumps.filt':
        raise ValueError('polish needs a dictionary (the identity filter)')
    for it in range(max_events):
        w = np.where(act, critical_w(L, al), 0.0)
        r = fstar - K.phi(c) @ w
        rms = float(np.sqrt(np.mean(r ** 2)))
        hist.append((tau, c.copy(), w.copy(), rms))
        if rms < tol:
            break
        v = K.phi(c).T @ K.filt(r) / nx
        dt = np.full(d, np.inf)
        ina = ~act & (v != 0)
        dt[ina] = (rho[ina] - np.sign(v[ina]) * u[ina]) / np.abs(v[ina])
        dmin = dt.min()
        if not np.isfinite(dmin):
            break
        ent = np.flatnonzero(dt <= dmin * (1 + tie_tol) + tie_tol * np.finfo(float).tiny)
        rmax = float(np.abs(v[~act]).max())
        entry_info += [dict(it=it, j=int(j), u0=float(u[j]), rate=float(v[j]), rate_max=rmax, t_wall=float(dt[j]), t_min=float(dmin))
                       for j in ent]
        u[~act] += dmin * v[~act]
        tau += dmin
        for j in ent:
            u[j] = np.sign(v[j]) * rho[j]
            act[j] = True
            al[j] = np.sign(v[j]) * (seed * L.m * L.A ** 2 / (2 * L.c)) ** (1 / 3)
            ev.append((tau, int(j), 'enter', float(c[j]))); ev_iter.append(it)
        if max_entries is not None and sum(e[2] == 'enter' for e in ev) >= max_entries:
            break
        idx = np.flatnonzero(act)
        n = idx.size
        z = np.r_[al[idx], c[idx]]
        elapsed, chunk, ok, polished, failed = 0.0, float(t_fast), False, not polish, False
        plog, n_rep = [], 0
        polish_log.append(plog)
        try:
            for _k in range(max_doublings + 1 + repolish):           # a repeated polish does not spend a doubling
                z, how = relax_chunk(rhs, z, chunk)
                n_stiff += how == 1
                n_redo += how == 2
                elapsed += chunk
                if how == 3:
                    n_fail += 1
                    failed = True
                    break
                if box is not None and np.any((z[n:] < box[0]) | (z[n:] > box[1])):
                    runaway = True
                    break
                if not polished:
                    w0 = critical_w(L, z[:n])
                    wp, cp, pinfo = polish_block(K, fstar, w0, z[n:], np.sign(u[idx]), cbox=box, max_nfev=polish_nfev, info=True,
                                                  gtol=polish_gtol)
                    seedy = np.flatnonzero(np.abs(w0) <= polish_wfrac * (np.abs(w0).max() if w_scale is None else w_scale))
                    r0 = K.filt(fstar - K.phi(z[n:]) @ w0)
                    good = all(abs(cp[i] - z[n + i]) <= polish_step * K.sig or
                               abs(cp[i] - ascent_target(K, r0, z[n + i], np.sign(u[idx[i]]))[0]) <= polish_step * K.sig for i in seedy)
                    plog.append(dict(pinfo, accepted=bool(good)))
                    if good:
                        z = np.r_[np.sign(wp) * _alpha_of(L, np.abs(wp), coef_crit), cp]
                        polished = True
                    else:
                        n_rej += 1
                wi = critical_w(L, z[:n])
                ri = K.filt(fstar - K.phi(z[n:]) @ wi)
                g = K.phi(z[n:]).T @ ri / nx
                h = K.dphi(z[n:]).T @ ri / nx
                big = np.abs(wi) > exit_tol * (np.abs(wi).max() if w_scale is None else w_scale)
                si = np.sign(u[idx])
                if (np.all(np.abs(g[big]) <= gtol * scale) and np.all(np.abs(wi[big] * h[big]) <= gtol * scale)
                        and np.all(si[~big] * g[~big] <= gtol * scale)):
                    if polish and n_rep < repolish:                    # the test passed: polish again (referee round 54)
                        n_rep += 1
                        e0 = float(np.sqrt(np.mean(ri ** 2)))
                        wp, cp, pinfo = polish_block(K, fstar, wi, z[n:], si, cbox=box, max_nfev=polish_nfev, info=True,
                                                      gtol=polish_gtol)
                        e1 = float(np.sqrt(np.mean(K.filt(fstar - K.phi(cp) @ wp) ** 2)))
                        seedy = np.flatnonzero(np.abs(wi) <= polish_wfrac * (np.abs(wi).max() if w_scale is None else w_scale))
                        good = e1 < e0 * (1 - 1e-9) and all(
                            abs(cp[i] - z[n + i]) <= polish_step * K.sig or
                            abs(cp[i] - ascent_target(K, ri, z[n + i], si[i])[0]) <= polish_step * K.sig for i in seedy)
                        plog.append(dict(pinfo, accepted=bool(good), e0=e0, e1=e1))
                        if good:
                            z = np.r_[np.sign(wp) * _alpha_of(L, np.abs(wp), coef_crit), cp]
                            chunk = elapsed
                            continue
                    ok = True
                    break
                chunk = elapsed
        except _OutOfBudget:
            out_of_budget = True
            relaxed.append(False)
            break
        relaxed.append(ok)
        al[idx] = z[:n]; c[idx] = z[n:]
        if runaway:
            break
        if failed:                    # a chunk failed: the block is left unrelaxed and no kernel exits on it (referee round 54)
            continue
        own = None
        if exit_kkt:
            Pz = K.phi(z[n:])
            own = np.array([Pz[:, q] @ K.filt(Pz[:, q]) for q in range(n)]) / nx
        for k, i in enumerate(idx):
            if leaves(si, g, wi, big, own)[k]:
                act[i] = False; al[i] = 0.0
                ev.append((tau, int(i), 'exit', float(c[i]))); ev_iter.append(it)
    return dict(events=ev, hist=hist, c=c, w=np.where(act, critical_w(L, al), 0.0), act=act, relaxed=relaxed, polish_rejected=n_rej,
                stiff_redone=n_stiff, ev_iter=ev_iter, runaway=runaway, out_of_budget=out_of_budget, nfev=nfev[0], solver_redone=n_redo,
                polish_log=polish_log, entry_info=entry_info,
                solver_failed=n_fail)


def wfr_inactive_rhs(L, K, fstar, rho):
    """The inactive kernels' limit under the Wasserstein–Fisher–Rao metric before any entry (r = f*), z = (u, c), masses εb (ρ = b^{−1/2}):
    in the core a ≈ α = a_max(εb/(2c))^{1/2}tan(πη/(2λ₁)), so du/dτ = ⟨f*φ(c)⟩, dc/dτ = (πa_max²/(4cρ))tan(πu/(2ρ))⟨f*∂_cφ(c)⟩."""
    rho = np.asarray(rho, float)
    k = np.pi * L.A ** 2 / (4 * L.c)

    def f(t, z):
        d = z.size // 2
        u, c = z[:d], z[d:]
        g = K.phi(c).T @ fstar / K.xs.size
        h = K.dphi(c).T @ fstar / K.xs.size
        return np.r_[g, k * np.tan(np.pi * u / (2 * rho)) * h / rho]
    return f


# ------------------------------------------------------------------------------------------------------------------------------------
# Proposition D‴.51 (PARADIGM.md §5d): the transit at a finite mass — near ties, the entrant's overshoot, the creep, the predicted ε
# ------------------------------------------------------------------------------------------------------------------------------------

def eta_crit(L, al):
    """The critical curve's mirror coordinate relative to the wall, η_c(α) = lim_{ε→0}(η_ε(α) − λ₁) = (α − a coth(α/a))/(2Gc):
    −a²/(2Gcα) near α = 0 (the core's tangent), (α − a)/(2Gc) for α ≫ a (the exponential branch)."""
    x = np.asarray(al, float) / L.A
    return (L.A / (2 * L.G * L.c)) * (x - 1 / np.tanh(x))


def coef_crit(L, al):
    """The coefficient w = βa/m on the critical curve β = 2c sinh²(α/a), a = a_max tanh(α/a) (|α|, overflow-safe)."""
    x = np.minimum(np.abs(np.asarray(al, float)) / L.A, 350.0)
    return 2 * L.c * np.sinh(x) ** 2 * L.A * np.tanh(x) / L.m


def _alpha_of(L, target, fn, lo=1e-12, hi=1e3, it=200):
    """α > 0 with fn(L, α) = target (fn increasing), by bisection in log α, vectorised."""
    t = np.asarray(target, float)
    a, b = np.full(t.shape, np.log(lo)), np.full(t.shape, np.log(hi))
    for _ in range(it):
        m_ = 0.5 * (a + b)
        up = fn(L, np.exp(m_)) < t
        a, b = np.where(up, m_, a), np.where(up, b, m_)
    return np.exp(0.5 * (a + b))


def inner_coef(L, xi, rho=1.0):
    """The inner profile: a kernel whose dual coordinate is u = ρ(1 − ε^{1/2}ξ) carries, as ε → 0 at fixed ξ, the coefficient
    𝒲(ξ) = w_c(α) with η_c(α) = −Λ₀ρξ, Λ₀ = Λ(ε)ε^{1/2} = (πa/(4G))(2/c)^{1/2}: ≈ K'³ρ⁻³ξ⁻³ for ξ ≫ 1 (the core's profile,
    K'³ = (8/π³)a/(m(2c)^{1/2})), ∝ exp(−πρ(2c)^{1/2}ξ) for ξ ≪ −1 (the exponential branch)."""
    lam0 = transit_scale(L, 1.0)
    al = _alpha_of(L, -lam0 * rho * np.asarray(xi, float), eta_crit)
    return coef_crit(L, al)


def inner_xi(L, w, rho=1.0):
    """The inverse of `inner_coef`: the inner dual coordinate ξ of a kernel carrying |w| near its wall."""
    al = _alpha_of(L, np.abs(np.asarray(w, float)), coef_crit)
    return -eta_crit(L, al) / (transit_scale(L, 1.0) * rho)


def overshoot(L, A, rho=1.0):
    """The entrant's overshoot: an entered kernel carrying the coefficient A has its dual coordinate beyond its wall by
    η_c(α_A)/Λ = ε^{1/2}·(−inner_xi(A)) — what a runner-up's dual is pulled along by while the entrant grows (Proposition D‴.51(i))."""
    return -inner_xi(L, A, rho)


def near_tie_inner(L, K, A, cstar, c_e, c_r, x, rho=1.0, xi0=400.0, s_end=5000.0, both=True, tol=1e-9, max_doublings=8):
    """The ε-free inner problem of a near tie (Proposition D‴.51(i)): an entrant at c_e and a runner-up at c_r, a bump A at cstar, in
    s = (τ − τ_e)/ε^{1/2} and ξ_i = (ρ − s_Au_i)/(ρε^{1/2}) (s_A the bump's sign): dξ_i/ds = −s_A⟨rφ_i⟩/ρ, dc_i/ds = Λ₀Gm w_i⟨r∂_cφ_i⟩,
    w_i = s_A𝒲(ξ_i), r = Aφ(c*) − Σ_iw_iφ(c_i); the runner-up starts x behind the entrant (ξ_r − ξ_e = x, from ξ_e = xi0). both=False
    freezes the runner-up (its coefficient dropped): the entrant alone, the runner-up's ξ then recording the pull. The start xi0 behind the
    wall truncates the entrant's approach, where its core coefficient 𝒲 ∝ ξ⁻³ already acts: an error ∝ 1/xi0 in the pull — the entrant's
    centre creeps by O(xi0⁻²) in that far field, which shifts its arrival by O(1/xi0) (the coefficient's direct action on the runner-up
    is O(xi0⁻²)); Richardson extrapolation in xi0 removes it. Integrated over the entrant's outer arrival time s₀ = xi0ρ/|⟨fφ_e⟩| plus s_end, then over the elapsed time again (doubling) until
    stationary — every rate below tol·|A| (a small |A| relaxes slowly) — or max_doublings is reached. Returns the end state with the
    runner-up's share |w_r|/(|w_e| + |w_r|), whether it became stationary, and (both=False) the matched pull — the runner-up's outer margin
    at the entrant's arrival, x + s₀(d_e − d_r) (d the drifts before the entry), less ξ_r at the end, its drift after the entry (the
    end residual's) subtracted over s − s₀ (`entry_pull`'s definition; = x − ξ_r(∞) at an exact tie onto a fitted bump)."""
    lam0 = transit_scale(L, 1.0)
    sA = np.sign(A)
    f_ = A * K.phi(np.array([cstar]))[:, 0]
    nx = K.xs.size

    def W(xi):
        return float(inner_coef(L, xi, rho))

    def rhs(s, z):
        xe, xr, ce, cr = z
        we = sA * W(xe)
        wr = sA * W(xr) if both else 0.0
        pe, pr = K.phi(np.array([ce]))[:, 0], K.phi(np.array([cr]))[:, 0]
        r = f_ - we * pe - wr * pr
        de, dr = K.dphi(np.array([ce]))[:, 0], K.dphi(np.array([cr]))[:, 0]
        return [-sA * (pe @ r) / nx / rho, -sA * (pr @ r) / nx / rho, lam0 * L.G * L.m * we * (de @ r) / nx,
                lam0 * L.G * L.m * wr * (dr @ r) / nx]
    pe0, pr0 = K.phi(np.array([c_e]))[:, 0], K.phi(np.array([c_r]))[:, 0]
    d_e, d_r = sA * (pe0 @ f_) / nx / rho, sA * (pr0 @ f_) / nx / rho
    s_arr = xi0 / d_e
    z, s0, chunk, ok = np.array([xi0, xi0 + x, c_e, c_r], float), 0.0, s_arr + float(s_end), False
    for _ in range(max_doublings + 1):
        sol = solve_ivp(rhs, (s0, s0 + chunk), z, method='LSODA', rtol=1e-10, atol=1e-12)
        z, s0 = sol.y[:, -1], s0 + chunk
        if np.all(np.abs(rhs(s0, z)) <= tol * abs(A)):
            ok = True
            break
        chunk = s0
    xe, xr, ce, cr = z
    we, wr = W(xe), (W(xr) if both else 0.0)
    d_end = -rhs(s0, z)[1]
    pull = (x + s_arr * (d_e - d_r)) - xr - d_end * (s0 - s_arr) if not both else np.nan
    return dict(xi_e=float(xe), xi_r=float(xr), c_e=float(ce), c_r=float(cr), w_e=we, w_r=wr, share=wr / (we + wr),
                status=int(sol.status), stationary=ok, s=s0, pull=float(pull), s_arrival=float(s_arr))


def path_duals(K, fstar, R):
    """The sliding path's (`sliding_iss`) dual coordinates reconstructed at its event times: on each interval the block's residual r_k is
    fixed, the inactive kernels' u grows at ⟨r_kφ(c_j)⟩, an entrant is pinned at its wall (±1 in units of ρ — equal masses). Returns
    dict(taus, res (r_k after the k-th relaxation), centres, U (the duals after the events at taus[k]), ever (the kernels that enter))."""
    hist = R['hist']
    taus = [h[0] for h in hist]
    res = [fstar - K.phi(h[1]) @ h[2] for h in hist]
    cs = [h[1] for h in hist]
    d = len(cs[0])
    nx = K.xs.size
    ever = sorted({j for t, j, k, c in R['events'] if k == 'enter'})
    u = np.zeros(d); act = np.zeros(d, bool); U = [u.copy()]
    ev = R['events']; e = 0
    for k in range(len(taus) - 1):
        v = K.phi(cs[k]).T @ K.filt(res[k]) / nx
        u = u + np.where(act, 0.0, (taus[k + 1] - taus[k]) * v)
        while e < len(ev) and ev[e][0] <= taus[k + 1] + 1e-12:
            t, j, kind, c = ev[e]
            if kind == 'enter':
                act[j] = True; u[j] = np.sign(u[j])
            else:
                act[j] = False
            e += 1
        U.append(u.copy())
    return dict(taus=taus, res=res, centres=cs, U=np.array(U), ever=ever)


def fit_hessian(K, fstar, w, c):
    """The fit E = ½⟨(f* − Φ(c)w)²⟩ at (w, c): its Jacobian J = [Φ(c), Φ′(c)diag(w)] on the samples and its Hessian in (w, c),
    H = (JᵀJ − Σ_x r(x)∇²f(x))/N — the residual's curvature included (∂²f/∂w_i∂c_i = ∂_cφ_i, ∂²f/∂c_i² = w_i∂²_cφ_i). Returns (J, H)."""
    n, nx = w.size, K.xs.size
    P, D = K.phi(c), K.dphi(c)
    D2 = P * (((K.xs[:, None] - c[None, :]) / K.sig ** 2) ** 2 - 1 / K.sig ** 2)
    r = fstar - P @ w
    J = np.c_[P, D * w[None, :]]
    C = np.zeros((2 * n, 2 * n))
    i = np.arange(n)
    C[i, n + i] = C[n + i, i] = r @ D
    C[n + i, n + i] = w * (r @ D2)
    return J, (J.T @ J - C) / nx


def block_response(K, fstar, w, c):
    """The residual's first-order response to a perturbation ξ of the data at a stationary point (w, c) of the block's fit
    ½⟨(f* − Φ(c)w)²⟩: δr = Rξ, R = I − JH⁻¹Jᵀ/N (`fit_hessian`; a nondegenerate stationary point). At an exact fit R is the projection off
    the tangent space. A dictionary only (the identity filter)."""
    nx = K.xs.size
    if w.size == 0:
        return np.eye(nx)
    J, H = fit_hessian(K, fstar, w, c)
    return np.eye(nx) - J @ np.linalg.solve(H, J.T) / nx


def block_rate(L, K, fstar, w, c, rho=1.0):
    """The slowest relaxation rate, in the inner time s, of the block's fast flow at its stationary point (w, c) (Lemma D‴.53a: a gradient
    flow of the fit with the mobilities |𝒲′(ξ_i)|/ρ for the coefficients and Λ₀Gm for the centres): the smallest eigenvalue of
    M^{1/2}HM^{1/2}."""
    J, H = fit_hessian(K, fstar, w, c)
    xi = inner_xi(L, w, rho)
    h = 1e-6 * np.maximum(1.0, np.abs(xi))
    dW = np.abs(inner_coef(L, xi + h, rho) - inner_coef(L, xi - h, rho)) / (2 * h)
    m = np.sqrt(np.r_[dW / rho, np.full(w.size, transit_scale(L, 1.0) * L.G * L.m)])
    return float(np.linalg.eigvalsh(m[:, None] * H * m[None, :]).min())


def path_response(K, fstar, R, c0):
    """Proposition D‴.56: the sliding path's first-order response to a perturbation ξ of the data (measurement noise). Along the path R
    (`sliding_iss` on f* with polish: entries only, at distinct times, each block relaxed to a nondegenerate stationary point of its fit)
    block i's residual moves by R_iξ (`block_response`), the inactive kernels' velocities v_j = ⟨r_iφ_j⟩ by δv_j = ⟨R_iξ, φ_j⟩, and so
    the entry times and the duals: the entrant j* reaches its wall at δτ_{i+1} = δτ_i − (δu_{j*}(τ_i) + (τ_{i+1} − τ_i)δv_{j*})/v_{j*},
    the others' duals move by δu_j(τ_{i+1}) = δu_j(τ_i) + (δτ_{i+1} − δτ_i)v_j + (τ_{i+1} − τ_i)δv_j. Returns dict(taus (the entry
    times), entrants, u (the duals at the last entry), U (d, N): ξ ↦ δu there, as rows — zero on the entered kernels —, u_all, U_all
    (the same at every entry: the duals are piecewise linear in τ, so their extremes over the path lie at the entries), T (M, N): the
    entry times' rows, Rlast: the residual's response at the last block relaxed, None if R stopped before relaxing it)."""
    if type(K).filt.__qualname__ != 'Bumps.filt':
        raise ValueError('path_response needs a dictionary (the identity filter)')
    ev, hist, nx = R['events'], R['hist'], K.xs.size
    if any(e[2] != 'enter' for e in ev):
        raise ValueError('path_response: the path has exits')
    taus, ent = np.array([e[0] for e in ev]), [e[1] for e in ev]
    if np.any(np.diff(taus) <= 0):
        raise ValueError('path_response: simultaneous entries')
    Phi0 = K.phi(np.asarray(c0, float))                                 # the inactive kernels never moved
    d, M = Phi0.shape[1], len(ev)
    u, U, T = np.zeros(d), np.zeros((d, nx)), np.zeros((M, nx))
    u_all, U_all = [], []
    t_prev, tau_prev = np.zeros(nx), 0.0
    for i in range(M):
        c, w, A = hist[i][1], hist[i][2], ent[:i]
        v = Phi0.T @ (fstar - K.phi(c) @ w) / nx
        dv = Phi0.T @ block_response(K, fstar, w[A], c[A]) / nx
        dt, j = taus[i] - tau_prev, ent[i]
        T[i] = t_prev - (U[j] + dt * dv[j]) / v[j]
        on = np.ones(d, bool); on[A] = False
        U[on] += np.outer(v[on], T[i] - t_prev) + dt * dv[on]
        u[on] += dt * v[on]
        U[j] = 0.0; u[j] = np.sign(v[j])
        u_all.append(u.copy()); U_all.append(U.copy())
        t_prev, tau_prev = T[i], taus[i]
    Rlast = block_response(K, fstar, hist[M][2][ent], hist[M][1][ent]) if len(hist) > M else None
    return dict(taus=taus, entrants=ent, u=u, U=U, u_all=u_all, U_all=U_all, T=T, Rlast=Rlast)


def creep_reaches(L, K, P, j, eps, c0j, b=1.0):
    """Whether kernel j (never entering on the path), started at c0j and integrated along the path's residuals with its creep law
    dc/dτ = ε(πa²/(4c))b^{3/2}cos(πμ/2)/sin³(πμ/2)·s⟨r∂_cφ⟩ (Proposition D‴.48(iii)(a)) and du/dτ = ⟨rφ(c)⟩, comes within the inner
    layer 3ε^{1/2} of its wall before the path's last event (P from `path_duals`)."""
    pref = np.pi * L.A ** 2 / (4 * L.c) * b ** 1.5
    nx = K.xs.size
    mu_in = 3 * np.sqrt(eps)
    z = np.array([0.0, c0j])
    for k in range(len(P['taus']) - 1):
        r = K.filt(P['res'][k])

        def f(t, z):
            u, c = z
            mu = max(1 - abs(u), mu_in)
            ph, dph = K.phi(np.array([c]))[:, 0], K.dphi(np.array([c]))[:, 0]
            return [ph @ r / nx, eps * pref * np.cos(np.pi * mu / 2) / np.sin(np.pi * mu / 2) ** 3 * np.sign(u) * (dph @ r) / nx]

        def hit(t, z):
            return 1 - abs(z[0]) - mu_in
        hit.terminal = True
        sol = solve_ivp(f, (P['taus'][k], P['taus'][k + 1]), z, method='LSODA', rtol=1e-8, atol=1e-11, events=hit)
        if sol.status == 1:
            return True
        z = sol.y[:, -1]
    return False


def entry_pull(L, K, fstar, R, P, k, i, tracers, rho=1.0, xi0=400.0, s_end=200.0, tol=1e-10, max_doublings=12, cross=None,
               rtol=1e-9, method='LSODA', scale_alpha=False, deviation=False):
    """The pull an entry exerts on kernels that do not enter, measured on the entry's ε-free inner problem (Proposition D‴.51(iv)): in
    s = (τ − τ_e)/ε^{1/2} the active block of the path before the event k (index into P's times) and the entrant i follow the loop on the
    critical curve — dα/ds = Λ₀·2Gc tanh²(α/a)⟨rφ⟩, dc/ds = Λ₀Gm w⟨r∂_cφ⟩, w = 2c sinh²(α/a)·a tanh(α/a)/m —, the block from its
    coefficients, the entrant from the inner coordinate xi0 behind its wall; the tracers j stay where they are (inactive kernels do not
    move), their inner coordinates ξ_j = μ_j/ε^{1/2} moving at dξ_j/ds = −s_j⟨rφ_j⟩/ρ. The pull Δ_j = lim_s (ξ_j^out(s) − ξ_j(s)), ξ_j^out
    the outer solution — the drift of the residual before the entry until the entrant's outer arrival at its wall, s₀ = xi0ρ/(s_i⟨r⁻φ_i⟩),
    the inner problem's final residual's after: kernel j forms at the mass ε iff μ_j/ε^{1/2} − Δ_j < ξ_θ (𝒲(ξ_θ) = θ). Integrated over
    s₀ + s_end, then doubling until the block's rates fall below tol (or max_doublings); rtol the integrator's relative tolerance (its
    absolute one 10⁻¹²), method its scheme. scale_alpha: the coefficients' coordinates α get the absolute tolerance 10⁻¹² times the
    entrant's starting |α| — in the far field α ≈ (m𝒲/(2c))^{1/3}a^{2/3} is small (≈ 2·10⁻³ at xi0 = 1600) and ξ ∝ 1/α there, so a fixed
    10⁻¹² limits ξ's accuracy to ~xi0·10⁻⁹ per step. deviation=True integrates instead the entrant's and the tracers' departures from
    their outer solutions — ζ_e = ξ_e − (xi0 − v₀min(s, s₀)) (the entrant's α from ξ_e by Newton on x − coth x = −(2Gc/a)Λ₀ρξ) and
    ζ_j = ξ_j − ξ_j^out(s), the outer rates after the entry the path's (⟨r_kφ_j⟩, the limit of the inner problem's by (H3)) — so that
    the pulls, Δ_j = −lim ζ_j, are O(1) quantities rather than differences of O(xi0) ones: the relative tolerance then bounds their error
    by ~rtol·|Δ| instead of rtol·xi0 accumulated over the steps (s_θ moved by 10⁻³ between rtol 10⁻⁹ and 10⁻¹⁰ in the plain form:
    referee round 41's unconverged predictions); [0, s₀] and the rest are integrated separately (the outer rates jump at s₀). Returns dict(pull = {j: Δ_j}, stationary, s, w_end, c_end (the block's and the entrant's coefficients and centres
    at the end)); with cross = θ also s_cross, the time at which the entrant's |w| first reaches θ, measured from its outer arrival s₀
    (Theorem D‴.53(c): the loop's crossing at τ_e + ε^{1/2}s_cross)."""
    lam0 = transit_scale(L, 1.0)
    nx = K.xs.size
    hist = R['hist']
    cb, wb = hist[k - 1][1], hist[k - 1][2]
    blk = [int(j) for j in np.flatnonzero(wb != 0) if j != i]
    idx = blk + [int(i)]
    n = len(idx)
    si = float(np.sign(P['U'][k][i])) or 1.0
    tr = [int(j) for j in tracers]
    st = np.array([np.sign(P['U'][k][j]) or 1.0 for j in tr])
    ctr = np.array([cb[j] for j in tr])
    Ptr = K.phi(ctr) if tr else np.zeros((nx, 0))
    al0 = np.r_[[np.sign(wb[j]) * _alpha_of(L, abs(wb[j]), coef_crit) for j in blk],
                si * _alpha_of(L, -lam0 * rho * xi0, eta_crit)]
    z0 = np.r_[al0, [cb[j] for j in idx], np.zeros(len(tr))]
    r_pre = K.filt(fstar - K.phi(cb) @ wb)
    d_pre = st * (Ptr.T @ r_pre) / nx / rho
    s0 = xi0 * rho / (si * (K.phi(np.array([cb[i]]))[:, 0] @ r_pre) / nx)

    def resid(z):
        al, c = z[:n], z[n:2 * n]
        w = coef_crit(L, al) * np.sign(al)
        return K.filt(fstar - K.phi(c) @ w), al, c, w

    def rhs(s, z):
        r, al, c, w = resid(z)
        g = K.phi(c).T @ r / nx
        h = K.dphi(c).T @ r / nx
        x = np.minimum(np.abs(al) / L.A, 350.0)
        return np.r_[lam0 * 2 * L.G * L.c * np.tanh(x) ** 2 * g, lam0 * L.G * L.m * w * h, -st * (Ptr.T @ r) / nx / rho]
    def hit(s, z):
        return coef_crit(L, z[n - 1]) - cross
    hit.direction = 1
    if deviation:
        return _entry_pull_deviation(L, K, fstar, R, k, i, idx, blk, si, tr, st, Ptr, cb, al0, r_pre, d_pre, s0, rho, xi0, s_end, tol,
                                     max_doublings, cross, rtol, method)
    atol = np.r_[np.full(n, 1e-12 * (min(1.0, abs(al0[-1])) if scale_alpha else 1.0)), np.full(n + len(tr), 1e-12)]
    z, s_, chunk, ok, s_cross = z0, 0.0, s0 + float(s_end), False, None
    for _ in range(max_doublings + 1):
        sol = solve_ivp(rhs, (s_, s_ + chunk), z, method=method, rtol=rtol, atol=atol, events=hit if cross is not None else None)
        if cross is not None and s_cross is None and sol.t_events[0].size:
            s_cross = float(sol.t_events[0][0]) - s0
        z, s_ = sol.y[:, -1], s_ + chunk
        if np.all(np.abs(rhs(s_, z)[:2 * n]) <= tol):
            ok = True
            break
        chunk = s_
    r_end, al, c, w = resid(z)
    d_post = st * (Ptr.T @ r_end) / nx / rho
    xi_out = -(d_pre * min(s_, s0) + d_post * max(s_ - s0, 0.0))
    pull = {j: float(xi_out[q] - z[2 * n + q]) for q, j in enumerate(tr)}
    return dict(pull=pull, stationary=ok, s=s_, w_end=dict(zip(idx, w.tolist())), c_end=dict(zip(idx, c.tolist())), s_cross=s_cross)


def _alpha_from_xi(L, xi, rho, lam0):
    """The critical curve's α ≥ 0 at the inner dual coordinate ξ (η_c(α) = −Λ₀ρξ): x = α/a solves x − coth x = y = −(2Gc/a)Λ₀ρξ;
    Newton from the root of x − 1/x = y, below the root (x − coth x is increasing and concave), converging monotonically."""
    y = -(2 * L.G * L.c / L.A) * lam0 * rho * np.asarray(xi, float)
    x = 0.5 * (y + np.sqrt(y * y + 4.0))
    for _ in range(60):
        th = np.tanh(x)
        dx = (x - 1.0 / th - y) * th * th
        x = x - dx
        if np.all(np.abs(dx) <= 1e-15 * np.maximum(x, 1e-300)):
            break
    return L.A * x


def _entry_pull_deviation(L, K, fstar, R, k, i, idx, blk, si, tr, st, Ptr, cb, al0, r_pre, d_pre, s0, rho, xi0, s_end, tol,
                          max_doublings, cross, rtol, method):
    """`entry_pull` with deviation=True: the state (α of the block, ζ_e, the centres, ζ_tracers)."""
    lam0 = transit_scale(L, 1.0)
    nx, n, nb = K.xs.size, len(idx), len(idx) - 1
    hk = R['hist'][k]
    d_post = st * (Ptr.T @ K.filt(fstar - K.phi(hk[1]) @ hk[2])) / nx / rho
    v0 = xi0 / s0

    def unpack(s, z, pre):
        xi_e = z[nb] + (xi0 - v0 * s if pre else 0.0)
        al = np.r_[z[:nb], si * _alpha_from_xi(L, xi_e, rho, lam0)]
        c = z[n:2 * n]
        w = coef_crit(L, al) * np.sign(al)
        return al, c, w, K.filt(fstar - K.phi(c) @ w)

    def rates(s, z, pre):
        al, c, w, r = unpack(s, z, pre)
        g = K.phi(c).T @ r / nx
        h = K.dphi(c).T @ r / nx
        x = np.minimum(np.abs(al) / L.A, 350.0)
        return al, c, w, r, g, lam0 * 2 * L.G * L.c * np.tanh(x) ** 2 * g, lam0 * L.G * L.m * w * h

    def make(pre):                                                     # the outer rates jump at s₀: one right-hand side per side
        def rhs(s, z):
            al, c, w, r, g, dal, dc = rates(s, z, pre)
            return np.r_[dal[:nb], -si * g[nb] / rho + (v0 if pre else 0.0), dc,
                         -st * (Ptr.T @ r) / nx / rho + (d_pre if pre else d_post)]

        def hit(s, z):
            return float(coef_crit(L, _alpha_from_xi(L, z[nb] + (xi0 - v0 * s if pre else 0.0), rho, lam0))) - cross
        hit.direction = 1
        return rhs, hit
    z = np.r_[al0[:nb], 0.0, [cb[j] for j in idx], np.zeros(len(tr))]
    atol = np.r_[np.full(nb, 1e-12), 1e-12, np.full(n + len(tr), 1e-12)]
    s_, ok, s_cross = 0.0, False, None
    spans = [s0, s0 + float(s_end)]
    for q in range(max_doublings + 2):
        s1 = spans[q] if q < 2 else 2 * s_
        rhs, hit = make(q == 0)
        sol = solve_ivp(rhs, (s_, s1), z, method=method, rtol=rtol, atol=atol, events=hit if cross is not None else None)
        if sol.status < 0:
            raise RuntimeError(f'entry_pull (deviation): {sol.message}')
        if cross is not None and s_cross is None and sol.t_events[0].size:
            s_cross = float(sol.t_events[0][0]) - s0
        z, s_ = sol.y[:, -1], s1
        if q >= 1:
            al, c, w, r, g, dal, dc = rates(s_, z, False)
            if np.all(np.abs(np.r_[dal, dc]) <= tol):
                ok = True
                break
    al, c, w, r = unpack(s_, z, False)
    pull = {j: float(-z[2 * n + q]) for q, j in enumerate(tr)}
    return dict(pull=pull, stationary=ok, s=s_, w_end=dict(zip(idx, w.tolist())), c_end=dict(zip(idx, c.tolist())), s_cross=s_cross)


def recovery_threshold(L, K, fstar, c0, R, th_abs, n_creep=4, eps_lo=1e-14, eps_hi=1e-3, n_bisect=14, measured_pull=True):
    """An estimate of the mass below which the loop forms exactly the sliding path's kernels, from the path alone (Proposition D‴.51(iv);
    the leading order — referee round 39 measured the loop's own threshold at 0.77 of it, median over 48 random targets; equal masses only,
    the path's duals pinned at ±1):
    (a) at each entry i, every kernel j that never enters is pulled toward its wall by the entrant's overshoot — it is formed when its
    margin μ_j(τ_i) < ε^{1/2}x_θ(A_i)κ_ij, x_θ(A) = ln(|A|/θ − 1)/κ (κ = πρ(2c)^{1/2} the inner layer's Gibbs rate, θ the formation
    threshold th_abs), κ_ij = s_j⟨Δr_iφ_j⟩/(s_i⟨Δr_iφ_i⟩) (Δr_i the residual the entry removes, φ_i the entrant after its relaxation) —
    the near-tie law's estimate, right for the entrant's near-tie partner, which slides onto the entrant's feature and shares it. The pull
    on a kernel that stays where it is is measured on the entry's ε-free inner problem (`entry_pull`, Δ_ij): x_ij = Δ_ij + ξ_θ
    (𝒲(ξ_θ) = θ), raised to x_θκ_ij when j is pulled (Δ_ij > 0) and receives at most the entrant's own share of the residual change
    (0 < κ_ij ≤ 1: a near-tie partner); j can form only if x_ij > 0. Where an entry rearranges the block — an active kernel sent on to
    another feature, the residual changing far from the entrant, ⟨Δr_iφ_i⟩ small — κ_ij is large and the near-tie estimate would
    overstate the pull many times. ε_entry = min (μ_j/x_ij)²; (b) the n_creep never-entering kernels closest to their walls, integrated with their creep along the
    path (`creep_reaches`): ε_creep the smallest ε (bisection in log ε) at which one reaches the inner layer. Returns dict(eps_entry,
    entry=(i, j, μ, κ_ij, x_θ, Δ_ij), eps_creep, creep={j: ε_j}, eps_c = min). measured_pull=False: the near-tie estimate alone."""
    P = path_duals(K, fstar, R)
    taus, res, cs, U = P['taus'], P['res'], P['centres'], P['U']
    d = len(c0)
    nx = K.xs.size
    never = [j for j in range(d) if j not in P['ever']]
    kap = np.pi * np.sqrt(2 * L.c)
    xi_th = float(inner_xi(L, th_abs))
    eps_entry, arg = np.inf, None
    for t, i, kind, c in R['events']:
        if kind != 'enter':
            continue
        k = int(np.argmin(np.abs(np.array(taus) - t)))
        dr = K.filt(res[k - 1] - res[k])
        gi = K.phi(np.array([cs[k][i]]))[:, 0] @ dr / nx
        Ai = abs(float(R['w'][i]))                                   # the entrant's coefficient at the path's end (0 if it left)
        if gi == 0 or Ai <= th_abs:
            continue
        xth = np.log(Ai / th_abs - 1) / kap if Ai > 2 * th_abs else 0.0
        pull = entry_pull(L, K, fstar, R, P, k, i, never)['pull'] if (measured_pull and never) else {}
        for j in never:
            sj = np.sign(U[k][j]) if U[k][j] != 0 else 1.0
            kij = sj * (K.phi(np.array([c0[j]]))[:, 0] @ dr / nx) / abs(gi)
            if j in pull:
                x_ij = pull[j] + xi_th
                if pull[j] > 0 and 0 < kij <= 1:
                    x_ij = max(x_ij, xth * kij)
            else:
                x_ij = xth * kij if (kij > 0 and xth > 0) else -np.inf
            if x_ij > 0:
                e_ = ((1 - abs(U[k][j])) / x_ij) ** 2
                if e_ < eps_entry:
                    eps_entry, arg = float(e_), (int(i), int(j), float(1 - abs(U[k][j])), float(kij), float(xth),
                                                 float(pull.get(j, np.nan)))
    closest = sorted(never, key=lambda j: min(1 - abs(U[k][j]) for k in range(len(taus))))[:n_creep]
    creep = {}
    for j in closest:
        if not creep_reaches(L, K, P, j, eps_hi, c0[j]):
            creep[j] = np.inf
            continue
        lo, hi = np.log(eps_lo), np.log(eps_hi)
        if creep_reaches(L, K, P, j, eps_lo, c0[j]):
            creep[j] = eps_lo
            continue
        for _ in range(n_bisect):
            mid = 0.5 * (lo + hi)
            if creep_reaches(L, K, P, j, np.exp(mid), c0[j]):
                hi = mid
            else:
                lo = mid
        creep[j] = float(np.exp(hi))
    eps_creep = min(creep.values()) if creep else np.inf
    return dict(eps_entry=eps_entry, entry=arg, eps_creep=eps_creep, creep=creep, eps_c=min(eps_entry, eps_creep))


def loop_with_duals(L, K, fstar, c0, eps, taus, b=None, levels=None, P=None, rtol=1e-11, atol=1e-12):
    """The free law with moving centres from α = 0 with the masses εb (Theorem D‴.53's loop), the dual coordinates integrated alongside:
    the state (α, c, δη), each kernel on its curve β = εb_i + 2c sinh²(α_i/a_max). Without P, δη = η (η̇_i = ⟨rφ(c_i)⟩ exactly); with the
    sliding path's duals P (`path_duals`), δη_i = η_i − Λu_i^path(t/Λ), integrated at ⟨rφ(c_i)⟩ minus the path's rate on each of its
    intervals (⟨r_kφ(c_i^k)⟩ for its inactive kernels, 0 on its walls) — a well-scaled deviation (O(1) where η is O(Λ)), split at the path's
    events. The tolerances are per component: α's absolute tolerance atol times the core's scale a_max(εb/(2c))^{1/2}, the centres' atol·10⁻²,
    δη's atol. Records on the grid `taus` (in τ = t/Λ(ε)) the coefficients, centres, u = u^path + δη/Λ and du = δη/Λ; with
    levels = {i: θ_i} the first time each |w_i| reaches θ_i. Returns dict(taus, w, c, u, du, cross = {i: τ or None}, Lam). The path's duals
    (`path_duals`) are pinned at ±1, i.e. in units of ρ for equal masses: with P, b must be uniform."""
    d = len(c0)
    bb = np.ones(d) if b is None else np.asarray(b, float)
    if P is not None and np.ptp(bb) > 0:
        raise ValueError('loop_with_duals: the path\'s duals are for equal masses')
    b0 = eps * bb
    Lam = transit_scale(L, eps)
    taus = np.asarray(taus, float)
    base = critical_rhs(L, K, fstar, 'paradigm', b0=b0)
    nx = K.xs.size
    al_core = L.A * np.sqrt(b0 / (2 * L.c))
    tol = np.r_[atol * al_core, np.full(d, 1e-2 * atol), np.full(d, atol)]
    if P is not None:
        T = np.asarray(P['taus'], float)
        rates = []
        for k in range(len(T)):
            act = np.abs(np.abs(P['U'][k]) - 1) < 1e-12
            v = K.phi(np.asarray(P['centres'][k], float)).T @ K.filt(P['res'][k]) / nx
            rates.append(np.where(act, 0.0, v))
        cuts = [t for t in T[1:] if t < taus[-1]]
    else:
        T, rates, cuts = np.array([0.0]), [np.zeros(d)], []

    def make(rate):
        def rhs(t, z):
            al, c = z[:d], z[d:2 * d]
            r = K.filt(fstar - K.phi(c) @ critical_w(L, al, b0))
            return np.r_[base(t, z[:2 * d]), K.phi(c).T @ r / nx - rate]
        return rhs
    lev = dict(levels or {})
    evs = []
    for i, th in lev.items():
        def ev(t, z, i=i, th=th):
            return abs(float(critical_w(L, z[i], b0[i]))) - th
        ev.direction = 1
        evs.append(ev)
    bounds = [0.0] + [Lam * t for t in cuts] + [Lam * taus[-1]]
    z = np.r_[np.zeros(d), np.asarray(c0, float), np.zeros(d)]
    Y = np.full((taus.size, 3 * d), np.nan)
    cross = {i: None for i in lev}
    for q in range(len(bounds) - 1):
        t0, t1 = bounds[q], bounds[q + 1]
        k = int(np.searchsorted(T, t0 / Lam + 1e-12, side='right') - 1) if P is not None else 0
        sel = (Lam * taus >= t0) & (Lam * taus < t1)
        if q == len(bounds) - 2:
            sel = (Lam * taus >= t0)
        te = np.unique(np.r_[Lam * taus[sel], t1])                       # the segment's end always evaluated: the next one starts there
        sol = solve_ivp(make(rates[k]), (t0, t1), z, method='LSODA', rtol=rtol, atol=tol, t_eval=te, events=evs or None)
        if sol.status < 0:
            raise RuntimeError(f'loop_with_duals: the segment [{t0 / Lam:.6g}, {t1 / Lam:.6g}] failed: {sol.message}')
        if sel.any():
            Y[sel] = sol.y.T[np.searchsorted(te, Lam * taus[sel])]
        for qq, i in enumerate(lev):
            if cross[i] is None and sol.t_events[qq].size:
                cross[i] = float(sol.t_events[qq][0]) / Lam
        z = sol.y[:, -1]
    al, c, deta = Y[:, :d], Y[:, d:2 * d], Y[:, 2 * d:]
    if P is not None:
        up = np.array([np.asarray(P['U'][int(np.searchsorted(T, t + 1e-12, side='right') - 1)])
                       + (t - T[int(np.searchsorted(T, t + 1e-12, side='right') - 1)])
                       * rates[int(np.searchsorted(T, t + 1e-12, side='right') - 1)] for t in taus])
    else:
        up = np.zeros_like(deta)
    return dict(taus=taus, w=critical_w(L, al, b0), c=c, u=up + deta / Lam, du=deta / Lam, cross=cross, Lam=Lam)


def creep_constants(L, K, P, j, c0j, taus, b=1.0):
    """The leading-order displacement of a kernel that does not enter, and its dual's, along the sliding path (Theorem D‴.53(a), (d)):
    c_j = c0j + εC_j(τ), u_j = u_j^path(τ) + (pulls) + εD_j(τ), with dC_j/dτ = κ(μ_j)s_j⟨r_k∂_cφ(c0j)⟩ (the creep law of D‴.48(iii)(a),
    κ(μ) = (πa_max²/(4c))b^{3/2}cos(πμ/2)/sin³(πμ/2), μ_j = 1 − |u_j| the path's margin, ρ = 1) and dD_j/dτ = C_j⟨r_k∂_cφ(c0j)⟩ on each of
    the path's intervals (P from `path_duals`). Returns (C, D) at `taus` (which must precede j's entry, if it enters)."""
    pref = np.pi * L.A ** 2 / (4 * L.c) * b ** 1.5
    nx = K.xs.size
    ph, dph = K.phi(np.array([c0j]))[:, 0], K.dphi(np.array([c0j]))[:, 0]
    T = list(P['taus'])
    taus = np.asarray(taus, float)
    out = np.full((taus.size, 2), np.nan)
    z = np.zeros(2)
    for k in range(len(T)):
        t0 = T[k]
        t1 = T[k + 1] if k + 1 < len(T) else max(taus.max(), t0) + 1e-12
        r = K.filt(P['res'][k])
        v, gp = ph @ r / nx, dph @ r / nx
        u0 = P['U'][k][j]

        def f(t, z, u0=u0, v=v, gp=gp, t0=t0):
            u = u0 + (t - t0) * v
            mu = 1 - abs(u)
            return [pref * np.cos(np.pi * mu / 2) / np.sin(np.pi * mu / 2) ** 3 * np.sign(u) * gp, z[0] * gp]
        sel = (taus >= t0) & (taus <= t1)
        if t1 > t0:
            sol = solve_ivp(f, (t0, t1), z, method='LSODA', rtol=1e-10, atol=1e-14, dense_output=True)
            if sel.any():
                out[sel] = sol.sol(taus[sel]).T
            z = sol.y[:, -1]
        if t1 >= taus.max():
            break
    return out[:, 0], out[:, 1]
