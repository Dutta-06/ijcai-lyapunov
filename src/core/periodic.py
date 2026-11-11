"""
Periodic orbits of the finite-m closed loop and the averaged system (PARADIGM.md Theorem J)
==========================================================================================

The real closed loop with saturated features (langevin.SaturatedLiftedNetwork), no dither and the leak σ = σ₀/k_e:

    ẋ = f(x) + ẋ_m − k_e e − f̂_Θ(x),     e := x − x_m(t),
    ϑ̇_i = G ∇_ϑ[eᵀΦ̃(x; ϑ_i)] − (σ₀/k_e)(ϑ_i − ϑ⁰),      i = 1..m,                                   (CL)

with f ∈ C² autonomous and x_m T₀-periodic.  Its averaged system is the gradient flow

    Θ̇ = −(1/k_e) G_m ∇J_m(Θ),      G_m := I_m ⊗ G,
    J_m(Θ) := m L̄_m(Θ) + (σ₀/2) Σ_i ‖ϑ_i − ϑ⁰‖²_{G⁻¹},     L̄_m(Θ) := ½ ∫ ‖f − (1/m)Σ_j Φ̃(·; ϑ_j)‖² dρ,          (AV)

ρ the occupation measure of x_m.  Theorem J: for k_e large, every critical point Θ* of J_m with ∇²J_m(Θ*)
non-singular has a unique T₀-periodic orbit of (CL) within O(1/k_e), with the same stability type; the period mean
of the orbit is Θ* − J⁻¹ b + O(k_e⁻²), with J the Jacobian of the leading averaged drift and b the first-order
correction of Proposition H.7 and Remark H.7a (frozen-measure terms plus the intra-period oscillation term).

This module provides: the closed-loop vector field, J_m with its analytic gradient and a finite-difference Hessian, a
critical-point finder, the Poincaré map of (CL) over one period (high-order adaptive integration), Newton shooting
for the periodic orbit together with the monodromy matrix (Floquet multipliers), and the first-order offset
prediction.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
from . import averaging as A


class ClosedLoop:
    """(CL) for a saturated network `net`, plant f_vec (vectorised, (N,n)→(N,n)), reference (x_ref, x_ref_dot, period)."""

    def __init__(self, net, f_vec, x_ref, x_ref_dot, period, ke, G_diag, sigma0, theta0=None):
        self.net, self.f_vec, self.xr, self.xrd = net, f_vec, x_ref, x_ref_dot
        self.T0, self.ke = float(period), float(ke)
        self.G = np.asarray(G_diag, float)
        self.sigma0 = float(sigma0)
        self.p = 2 * net.n + 1
        self.theta0 = np.zeros(self.p) if theta0 is None else np.asarray(theta0, float)

    # -------------------------------------------------------------- the vector field
    def rhs(self, t, z, m):
        n = self.net.n
        x, Th = z[:n], z[n:].reshape(m, self.p)
        xr, xrd = self.xr(t), self.xrd(t)
        e = x - xr
        fx = self.f_vec(x[None, :])[0]
        xdot = fx + xrd - self.ke * e - self.net.f_hat(x, Th)
        Thdot = self.G[None, :] * self.net.grad_g(x, e, Th) - (self.sigma0 / self.ke) * (Th - self.theta0[None, :])
        return np.concatenate([xdot, Thdot.ravel()])

    def integrate(self, z0, t0, t1, m, rtol=1e-11, atol=1e-13):
        sol = solve_ivp(lambda t, z: self.rhs(t, z, m), (t0, t1), z0, method='DOP853', rtol=rtol, atol=atol)
        if not sol.success:
            raise RuntimeError(sol.message)
        return sol.y[:, -1]

    def poincare(self, z0, m, **kw):
        """One period of (CL) from t = 0."""
        return self.integrate(z0, 0.0, self.T0, m, **kw)

    def simulate_periods(self, z0, m, n_periods, **kw):
        """Iterate the Poincaré map; returns the list of states at t = 0, T₀, 2T₀, …"""
        out = [np.array(z0, float)]
        for _ in range(n_periods):
            out.append(self.poincare(out[-1], m, **kw))
        return out

    # -------------------------------------------------------------- the averaged objective
    def orbit(self, N=2000, sub=1):
        return A.Orbit(self.xr, self.xrd, self.T0, N=N, sub=sub)

    def J_m(self, Th, orb):
        """J_m(Θ) = m L̄_m(Θ) + (σ₀/2) Σ_i ‖ϑ_i − ϑ⁰‖²_{G⁻¹}  (ρ by the orbit's grid)."""
        m = Th.shape[0]
        r = self.f_vec(orb.xm) - self.net.f_hat_batch(orb.xm, Th)
        Lbar = 0.5 * float(np.mean(np.sum(r ** 2, axis=1)))
        leak = 0.5 * self.sigma0 * float(np.sum((Th - self.theta0[None, :]) ** 2 / self.G[None, :]))
        return m * Lbar + leak

    def grad_J_m(self, Th, orb):
        """∇_{ϑ_i} J_m = −⟨∇_ϑ(r_Θᵀ Φ̃(·; ϑ_i))⟩_ρ + σ₀ G⁻¹(ϑ_i − ϑ⁰)   (m, p)."""
        r = self.f_vec(orb.xm) - self.net.f_hat_batch(orb.xm, Th)                 # (N, n)
        g = np.zeros_like(Th)
        for x, rr in zip(orb.xm, r):
            g += self.net.grad_g(x, rr, Th)
        g /= orb.N
        return -g + self.sigma0 * (Th - self.theta0[None, :]) / self.G[None, :]

    def hess_J_m(self, Th, orb, h=1e-5):
        """Finite-difference Hessian of J_m (central differences of the analytic gradient), (mp, mp), symmetrised."""
        m, p = Th.shape
        d = m * p
        H = np.zeros((d, d))
        base = Th.ravel()
        for j in range(d):
            e = np.zeros(d); e[j] = h
            H[:, j] = (self.grad_J_m((base + e).reshape(m, p), orb) - self.grad_J_m((base - e).reshape(m, p), orb)).ravel() / (2 * h)
        return 0.5 * (H + H.T)

    def averaged_drift(self, Th, orb):
        """The leading averaged drift (1/k_e)·(−G_m ∇J_m)(Θ), (m, p)."""
        return -(self.G[None, :] * self.grad_J_m(Th, orb)) / self.ke

    def find_critical_point(self, Th0, orb, gtol=1e-12):
        """Minimise J_m from Th0 (L-BFGS with the analytic gradient), then Newton-polish; returns (Θ*, ∇²J_m(Θ*))."""
        m, p = Th0.shape
        fun = lambda v: self.J_m(v.reshape(m, p), orb)
        jac = lambda v: self.grad_J_m(v.reshape(m, p), orb).ravel()
        res = minimize(fun, Th0.ravel(), jac=jac, method='L-BFGS-B', options=dict(maxiter=5000, gtol=1e-10, ftol=1e-15))
        v = res.x
        for _ in range(20):                                              # Newton polish to machine precision
            g = jac(v)
            if np.linalg.norm(g) < gtol:
                break
            H = self.hess_J_m(v.reshape(m, p), orb)
            v = v - np.linalg.solve(H, g)
        Th_star = v.reshape(m, p)
        return Th_star, self.hess_J_m(Th_star, orb)

    # -------------------------------------------------------------- the periodic orbit of (CL)
    def shoot(self, z0, m, max_iter=12, tol=1e-11, fd=1e-7, verbose=False):
        """
        Newton on P(z) − z = 0 with a finite-difference Jacobian.  Returns (z_p, DP, residual) where DP is the Jacobian
        of the Poincaré map at the fixed point (the monodromy matrix of the periodic orbit).
        """
        z = np.array(z0, float)
        d = len(z)
        DP = None
        for it in range(max_iter):
            Pz = self.poincare(z, m)
            res = Pz - z
            if verbose:
                print(f"  shoot iter {it}: residual {np.linalg.norm(res):.3e}")
            DP = np.zeros((d, d))
            for j in range(d):
                ej = np.zeros(d); ej[j] = fd * max(1.0, abs(z[j]))
                DP[:, j] = (self.poincare(z + ej, m) - self.poincare(z - ej, m)) / (2 * ej[j])
            if np.linalg.norm(res) < tol:
                return z, DP, float(np.linalg.norm(res))
            z = z - np.linalg.solve(DP - np.eye(d), res)
        Pz = self.poincare(z, m)
        return z, DP, float(np.linalg.norm(Pz - z))

    def period_mean(self, z_p, m, N=400):
        """Period mean of Θ along the periodic orbit through z_p (dense output on N points)."""
        sol = solve_ivp(lambda t, z: self.rhs(t, z, m), (0.0, self.T0), z_p, method='DOP853', rtol=1e-11, atol=1e-13,
                        t_eval=np.linspace(0, self.T0, N, endpoint=False))
        n = self.net.n
        Th = sol.y[n:, :].T.reshape(N, m, self.p)
        E = sol.y[:n, :].T - np.array([self.xr(t) for t in sol.t])
        return Th.mean(axis=0), Th, E, sol.t

    # -------------------------------------------------------------- first-order prediction (H.7 + H.7a)
    def offset_prediction(self, Th_star, orb, jac_vec, h=1e-5):
        """
        Y* − Θ* ≈ −J⁻¹ (b_frozen + b_osc), with J the Jacobian of the leading averaged drift at Θ* (in time t, so it
        carries 1/k_e), b_frozen,i := (1/k_e²) G ∇_ϑ(U_even + U_odd)[μ*](ϑ_i) (Proposition H.7) and b_osc the oscillation
        term of Remark H.7a (also with its 1/k_e²).  Returns (full prediction, frozen-only prediction), each (m, p).
        """
        m, p = Th_star.shape
        fm = A.FrozenMeasure(self.net, Th_star)

        def grad_at_particles(fun):
            E = np.vstack([Th_star + h * np.eye(p)[k][None, :] for k in range(p)] + [Th_star - h * np.eye(p)[k][None, :] for k in range(p)])
            v = fun(E).reshape(2, p, m)
            return (v[0] - v[1]).T / (2 * h)

        g = grad_at_particles(lambda E: sum(A.expansion_terms(self.f_vec, jac_vec, fm, orb, E)))
        b_frozen = (self.G[None, :] * g / self.ke ** 2).ravel()
        b_osc = (A.oscillation_correction(self.net, self.f_vec, orb, Th_star, self.G) / self.ke ** 2).ravel()
        d = m * p
        J = np.zeros((d, d))
        base = Th_star.ravel()
        for j in range(d):
            e = np.zeros(d); e[j] = h
            J[:, j] = (self.averaged_drift((base + e).reshape(m, p), orb) - self.averaged_drift((base - e).reshape(m, p), orb)).ravel() / (2 * h)
        full = (-np.linalg.solve(J, b_frozen + b_osc)).reshape(m, p)
        frozen = (-np.linalg.solve(J, b_frozen)).reshape(m, p)
        return full, frozen


def floquet_exponents(DP, period, n_state):
    """
    Floquet exponents −log|λ|/T₀ of the monodromy matrix DP, split into the n_state fastest (error channel) and the
    rest (parameters); sorted ascending within each group.
    """
    lam = np.linalg.eigvals(DP)
    ex = np.sort(-np.log(np.abs(lam)) / period)
    return ex[-n_state:], ex[:-n_state]


def predicted_parameter_exponents(H, G_diag, m, ke):
    """Eigenvalues of (1/k_e) G_m^{1/2} ∇²J_m G_m^{1/2}, the exponents of the averaged linearisation (ascending)."""
    sq = np.sqrt(np.tile(np.asarray(G_diag, float), m))
    M = (sq[:, None] * H * sq[None, :]) / ke
    return np.sort(np.linalg.eigvalsh(M))


# ---------------------------------------------------------------- collapsed configurations and the orientation of the orbit (Proposition K)
def separation_projector(m, p):
    """Orthogonal projector of ℝ^{mp} onto the separation subspace {(u_1..u_m) : Σ_i u_i = 0} (dimension (m − 1)p)."""
    P = np.eye(m) - np.ones((m, m)) / m
    return np.kron(P, np.eye(p))


def restricted_exponents(DP, period, m, p, n_state):
    """
    Floquet exponents of the parameter block of the monodromy matrix DP, split by the permutation symmetry of a
    collapsed periodic orbit into the separation modes (Σ_i u_i = 0) and the collective modes (u_i equal).  The state
    block is removed by restricting to the parameter coordinates (valid up to the O(e^{−k_e T₀}) coupling).
    Returns (separation exponents, collective exponents), each sorted ascending, as −log|λ|/T₀.
    """
    M = DP[n_state:, n_state:]
    Psep = separation_projector(m, p)
    Pcol = np.eye(m * p) - Psep
    # eigenvalues of the compressions; the blocks are invariant up to the numerical symmetry breaking
    def ex(Pr, dim):
        # orthonormal basis of the range of Pr
        U, s, _ = np.linalg.svd(Pr)
        B = U[:, :dim]
        lam = np.linalg.eigvals(B.T @ M @ B)
        return np.sort(-np.log(np.abs(lam)) / period)
    return ex(Psep, (m - 1) * p), ex(Pcol, p)


def collapsed_critical_point(cl, orb, theta_guess):
    """The single-kernel critical point ϑ̄ of J_1 (with the leak) — a collapsed critical point of J_m for every m."""
    Th_star, H = cl.find_critical_point(np.atleast_2d(theta_guess), orb)
    return Th_star[0], H


def first_variation_hessian(cl, orb, theta_bar, h=1e-5):
    """
    R₂ := ∇²_ϑ⟨r_μ̄ᵀΦ̃(·;ϑ)⟩_ρ at ϑ = ϑ̄ with the residual r_μ̄ = f − Φ̃(·;ϑ̄) FROZEN (the Hessian of −V_L[μ̄] at the atom).
    The leading-order separation Hessian of J_m at the collapsed configuration is σ₀G⁻¹ − R₂, independent of m.
    """
    p = len(theta_bar)
    Th = theta_bar[None, :]
    r = cl.f_vec(orb.xm) - cl.net.f_hat_batch(orb.xm, Th)                       # frozen residual (N, n)

    def grad_V(th):                                                               # ∇_ϑ⟨rᵀΦ̃(·;ϑ)⟩ at frozen r
        g = np.zeros(p)
        for x, rr in zip(orb.xm, r):
            g += cl.net.grad_g(x, rr, th[None, :])[0]
        return g / orb.N
    R2 = np.zeros((p, p))
    for j in range(p):
        e = np.zeros(p); e[j] = h
        R2[:, j] = (grad_V(theta_bar + e) - grad_V(theta_bar - e)) / (2 * h)
    return 0.5 * (R2 + R2.T)


def sigma_crit(cl, orb, theta_guess, iters=30):
    """
    The leak gain at which the softest separation mode of the collapsed configuration is marginal at leading order:
    σ_crit solves σ = λ_max(G^{1/2} R₂(ϑ̄(σ)) G^{1/2}), with ϑ̄(σ) the single-kernel critical point for the leak σ.
    Secant iteration from the current cl.sigma0; returns (σ_crit, ϑ̄, R₂).
    """
    sq = np.sqrt(cl.G)

    def gap(s):
        cl.sigma0 = s
        tb, _ = collapsed_critical_point(cl, orb, theta_guess)
        R2 = first_variation_hessian(cl, orb, tb)
        return s - np.linalg.eigvalsh(sq[:, None] * R2 * sq[None, :])[-1], tb, R2
    s0 = cl.sigma0
    g0, tb, R2 = gap(s0)
    s1 = s0 - g0
    g1, tb, R2 = gap(s1)
    for _ in range(iters):
        if abs(g1) < 1e-12 or abs(g1 - g0) < 1e-15:
            break
        s2 = s1 - g1 * (s1 - s0) / (g1 - g0)
        s0, g0 = s1, g1
        s1 = s2
        g1, tb, R2 = gap(s1)
    cl.sigma0 = s1
    return s1, tb, R2


def second_order_drift(cl, Th, orb, jac_vec, h=1e-5):
    """The averaged drift through order k_e⁻²: (1/k_e)⟨H⟩(Θ) + b_frozen(Θ) + b_osc(Θ), (m, p)."""
    m, p = Th.shape
    fm = A.FrozenMeasure(cl.net, Th)
    E = np.vstack([Th + h * np.eye(p)[k][None, :] for k in range(p)] + [Th - h * np.eye(p)[k][None, :] for k in range(p)])
    v = sum(A.expansion_terms(cl.f_vec, jac_vec, fm, orb, E)).reshape(2, p, m)
    g = (v[0] - v[1]).T / (2 * h)
    b_frozen = cl.G[None, :] * g / cl.ke ** 2
    b_osc = A.oscillation_correction(cl.net, cl.f_vec, orb, Th, cl.G) / cl.ke ** 2
    return cl.averaged_drift(Th, orb) + b_frozen + b_osc


def jacobian_fd(fun, Th, h=1e-5):
    """Finite-difference Jacobian (mp × mp) of a drift Θ ↦ (m, p)."""
    m, p = Th.shape
    d = m * p
    J = np.zeros((d, d))
    base = Th.ravel()
    for j in range(d):
        e = np.zeros(d); e[j] = h
        J[:, j] = (fun((base + e).reshape(m, p)) - fun((base - e).reshape(m, p))).ravel() / (2 * h)
    return J
