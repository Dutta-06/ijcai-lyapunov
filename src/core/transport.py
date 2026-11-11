"""
Transport adaptive control: lifted particles, optimal assignment, OT Lyapunov function
=======================================================================================

Implements PARADIGM.md.  A kernel is a lifted particle ϑ = (a, c, s) ∈ ℝⁿ×ℝⁿ×ℝ,
the network is the empirical measure μ = (1/m)Σδ_{ϑ_i}, f̂(x) = (1/m)Σ a_i φ(x;c_i,s_i).

    law (L):   ϑ̇_i = Proj_K(ϑ_i, G ∇g_t(ϑ_i)),   g_t(ϑ) = eᵀ a φ(x; c, s)
    law (L_σ): ϑ̇_i = Proj_K(ϑ_i, G ∇g_t(ϑ_i) − σ_m (ϑ_i − ϑ⁰_{τ(i)})),  τ = optimal assignment to the prior

    V_OT = ½‖e‖² + ½ W(μ, μ*)²,   W² = min_π (1/m)Σ ‖ϑ_i − ϑ*_{π(i)}‖²_{G⁻¹}  (linear assignment)
"""
import numpy as np
from scipy.optimize import linear_sum_assignment
from .projection import ProjectionOperator, BandwidthProjection


# ------------------------------------------------------------------ Gaussian constants
def gaussian_constants(s_min: float, a_max: float):
    """
    Closed-form suprema over u = ‖x−c‖²/s² ≥ 0 (φ = e^{−u/2}) of the derivative magnitudes, scaled by s:

        h_cc = sup e^{−u/2} max(|u−1|, 1) = 1                    (at u = 0)          ‖∂²_cc φ‖ s²
        h_cs = sup √u e^{−u/2} |u−2|      at u = (5 − √17)/2                          ‖∂²_cs φ‖ s²
        h_ss = sup e^{−u/2} |u² − 3u|     = 2e^{−1/2}  (at u = 1)                     |∂²_ss φ| s²
        J_c  = e^{−1/2}/s_min ≥ ‖∇_c φ‖,   J_s = 2e^{−1}/s_min ≥ |∂_s φ|,   J := (J_c² + J_s²)^{1/2}

    H_φ := (max(h_cc, h_ss) + h_cs)/s_min² ≥ ‖∇²_{(c,s)} φ‖ (an upper bound, not the supremum),
    H_K := J + a_max H_φ ≥ sup_{x, ϑ∈Θ_K, ‖u‖=1} ‖∇²_ϑ (uᵀΦ(x;ϑ))‖,  K := H_K/2   (PARADIGM.md eq. (1)).
    """
    u_cs = (5.0 - np.sqrt(17.0)) / 2.0
    h_cc = 1.0
    h_cs = float(np.sqrt(u_cs) * np.exp(-u_cs / 2) * abs(u_cs - 2.0))
    h_ss = float(2.0 * np.exp(-0.5))
    J_c = np.exp(-0.5) / s_min
    J_s = 2.0 * np.exp(-1.0) / s_min
    J = float(np.sqrt(J_c ** 2 + J_s ** 2))
    H_phi = (max(h_cc, h_ss) + h_cs) / s_min ** 2
    H_K = J + a_max * H_phi
    return dict(J_c=J_c, J_s=J_s, J=J, H_phi=float(H_phi), H_K=float(H_K), K=float(0.5 * H_K),
                h_cc=h_cc, h_cs=h_cs, h_ss=h_ss)


# ------------------------------------------------------------------ features
def phi_and_grads(x, P):
    """P: (m, p) lifted particles [a (n), c (n), s]. Returns phi (m,), grad of φ wrt (c, s): gc (m,n), gs (m,)."""
    n = (P.shape[1] - 1) // 2
    c = P[:, n:2 * n]; s = P[:, 2 * n]
    diff = x[None, :] - c
    sq = np.sum(diff ** 2, axis=1)
    phi = np.exp(-sq / (2.0 * s ** 2))
    gc = (phi / s ** 2)[:, None] * diff
    gs = phi * sq / s ** 3
    return phi, gc, gs


def f_hat(x, P):
    n = (P.shape[1] - 1) // 2
    phi, _, _ = phi_and_grads(x, P)
    return (phi[:, None] * P[:, :n]).sum(axis=0) / P.shape[0]


def grad_g(x, e, P):
    """∇_ϑ g_t(ϑ_i) for all particles, g_t(ϑ) = eᵀ a φ(x;c,s).  Returns (m, p)."""
    n = (P.shape[1] - 1) // 2
    phi, gc, gs = phi_and_grads(x, P)
    ae = P[:, :n] @ e                                  # (m,)
    out = np.empty_like(P)
    out[:, :n] = phi[:, None] * e[None, :]
    out[:, n:2 * n] = ae[:, None] * gc
    out[:, 2 * n] = ae * gs
    return out


def grad_Phi(x, P):
    """∇_ϑ Φ(x;ϑ_i) ∈ ℝ^{n×p} for each particle: returns (m, n, p)."""
    m, p = P.shape; n = (p - 1) // 2
    phi, gc, gs = phi_and_grads(x, P)
    out = np.zeros((m, n, p))
    for i in range(m):
        out[i, :, :n] = phi[i] * np.eye(n)
        out[i, :, n:2 * n] = np.outer(P[i, :n], gc[i])
        out[i, :, 2 * n] = P[i, :n] * gs[i]
    return out


# ------------------------------------------------------------------ optimal assignment
def assignment_cost_matrix(P, Q, Ginv_diag):
    d = P[:, None, :] - Q[None, :, :]
    return np.einsum('ijk,k,ijk->ij', d, Ginv_diag, d)


def optimal_assignment(P, Q, Ginv_diag):
    """π minimising (1/m)Σ ‖P_i − Q_{π(i)}‖²_{G⁻¹}; returns (pi, W²)."""
    C = assignment_cost_matrix(P, Q, Ginv_diag)
    rows, cols = linear_sum_assignment(C)
    pi = np.empty(P.shape[0], dtype=int); pi[rows] = cols
    return pi, float(C[rows, cols].sum() / P.shape[0])


def w2(P, Q, Ginv_diag):
    return np.sqrt(optimal_assignment(P, Q, Ginv_diag)[1])


def assignment_cost(P, Q, pi, Ginv_diag):
    d = P - Q[pi]
    return float(np.einsum('ik,k,ik->', d, Ginv_diag, d) / P.shape[0])


def all_assignment_costs(P, Q, Ginv_diag):
    """Brute force over S_m (m ≤ 7) — for gap computation in tests."""
    from itertools import permutations
    return np.array([assignment_cost(P, Q, np.array(pi), Ginv_diag) for pi in permutations(range(P.shape[0]))])


def assignment_gap(P_star, P0, Ginv_diag):
    costs = np.sort(all_assignment_costs(P_star, P0, Ginv_diag))
    b2 = costs[0]
    above = costs[costs > b2 + 1e-12]
    return (float(above[0] - b2) if len(above) else np.inf), float(b2)


# ------------------------------------------------------------------ optimal plans and barycentric projection (Theorem C′)
def optimal_plan(P, Q, Ginv_diag):
    """
    Optimal transport plan γ (m × m₀) between μ = (1/m)Σδ_{P_i} and ν = (1/m₀)Σδ_{Q_j} for the cost ‖P_i − Q_j‖²_{G⁻¹}:
    a linear program on the transportation polytope (row sums 1/m, column sums 1/m₀), solved with the HiGHS dual simplex
    (so the plan is a vertex, as the law of Theorem C′ requires).  Returns
    (γ, W²).  For m = m₀ a vertex is a permutation matrix divided by m (Birkhoff), so this reduces to optimal_assignment.
    """
    from scipy.optimize import linprog
    m, m0 = P.shape[0], Q.shape[0]
    C = assignment_cost_matrix(P, Q, Ginv_diag)
    A_eq = np.zeros((m + m0, m * m0))
    for i in range(m):
        A_eq[i, i * m0:(i + 1) * m0] = 1.0
    for j in range(m0):
        A_eq[m + j, j::m0] = 1.0
    b_eq = np.concatenate([np.full(m, 1.0 / m), np.full(m0, 1.0 / m0)])
    res = linprog(C.ravel(), A_eq=A_eq[:-1], b_eq=b_eq[:-1], bounds=(0, None), method='highs-ds')   # one constraint is redundant;
    if not res.success:                                                                                 # dual simplex returns a vertex
        raise RuntimeError(res.message)
    gamma = np.clip(res.x, 0.0, None).reshape(m, m0)                                                   # round-off can give −1e-9
    return gamma, float(np.sum(gamma * C))


def barycentric_projection(P, Q, Ginv_diag):
    """
    T̄_i := m Σ_j γ_{ij} Q_j, the barycentre of the prior mass that the optimal plan assigns to particle i.
    Returns (T̄ (m, p), γ, W², Var_γ) with Var_γ := Σ_{ij} γ_{ij} ‖Q_j − T̄_i‖²_{G⁻¹}, the within-cell variance, so that
    W² = (1/m)Σ_i‖P_i − T̄_i‖²_{G⁻¹} + Var_γ   (the Pythagorean decomposition of the plan cost).
    """
    gamma, W2 = optimal_plan(P, Q, Ginv_diag)
    Tbar = P.shape[0] * (gamma @ Q)
    d = Q[None, :, :] - Tbar[:, None, :]
    Var = float(np.einsum('ij,ijk,k,ijk->', gamma, d, Ginv_diag, d))
    return Tbar, gamma, W2, Var


def centroid_defect(P, Q, Ginv_diag):
    """V_bar² := (1/m)Σ_i‖P_i − T̄_i‖²_{G⁻¹} = W² − Var_γ: zero exactly at the equal-mass centroidal configurations of ν."""
    Tbar, gamma, W2, Var = barycentric_projection(P, Q, Ginv_diag)
    d = P - Tbar
    return float(np.einsum('ik,k,ik->', d, Ginv_diag, d) / P.shape[0]), Tbar, gamma, W2, Var


def lloyd_equal_mass(P0, Q, Ginv_diag, n_iter=200, tol=1e-12):
    """Equal-mass Lloyd iteration P ← T̄(P) (the time-discrete form of the pure barycentric leak); returns (P, W²)."""
    P = np.array(P0, float)
    for _ in range(n_iter):
        Tbar, _, W2, _ = barycentric_projection(P, Q, Ginv_diag)
        if np.max(np.abs(Tbar - P)) < tol:
            P = Tbar
            break
        P = Tbar
    return P, barycentric_projection(P, Q, Ginv_diag)[2]


# ------------------------------------------------------------------ controller
class TransportController:
    """
    Lifted-particle adaptive controller (PARADIGM.md law (L) / (L_σ)).

    Parameters
    ----------
    m, n, ke : as usual.
    g_a, g_c, g_s : lifted gains (metric G = diag(g_a I, g_c I, g_s)).
    a_max, c_max, s_min, s_max, eps : per-particle box Θ_K and boundary layer.
    leak : None | 'ot' | 'labelled' | 'ot_bary' ; sigma_m : leak gain ; prior : (m, p) prior particles, or (m₀, p) with
           m₀ ≠ m for 'ot_bary' (the leak toward the barycentric projection of the optimal plan, Theorem C′).
    """

    def __init__(self, m, n, ke=5.0, g_a=1.0, g_c=1.0, g_s=1.0,
                 a_max=10.0, c_max=5.0, s_min=0.3, s_max=3.0, eps=0.1,
                 leak=None, sigma_m=0.0, prior=None):
        self.m, self.n, self.ke = int(m), int(n), float(ke)
        self.p = 2 * n + 1
        self.G_diag = np.concatenate([np.full(n, float(g_a)), np.full(n, float(g_c)), [float(g_s)]])
        self.Ginv_diag = 1.0 / self.G_diag
        self.lam_G = float(self.G_diag.max())
        self.a_max, self.c_max, self.s_min, self.s_max, self.eps = map(float, (a_max, c_max, s_min, s_max, eps))
        self.proj_a = ProjectionOperator(a_max, eps)
        self.proj_c = ProjectionOperator(c_max, eps)
        self.proj_s = BandwidthProjection(s_min, s_max, eps)
        self.leak, self.sigma_m = leak, float(sigma_m)
        self.prior = None if prior is None else np.array(prior, float)
        cst = gaussian_constants(s_min, a_max)
        self.K, self.J, self.H_phi = cst['K'], cst['J'], cst['H_phi']
        self.L0 = 2.0 * self.a_max                                   # Lemma B.2, linear branch
        self.L1 = float(np.sqrt(1.0 + (self.a_max * self.J) ** 2))

    def remainder_bound(self, e_norm, sq_matched_each):
        """Per-particle bound of Lemma B.2: ‖e‖·min(K‖ϑ̃‖², L0 + L1‖ϑ̃‖)."""
        d = np.sqrt(sq_matched_each)
        return e_norm * np.minimum(self.K * sq_matched_each, self.L0 + self.L1 * d)

    def excess_bound(self, e_norm, W):
        """Averaged bound of Theorem B: ‖e‖·min(K λ_G W², L0 + L1 λ_G^{1/2} W)."""
        return e_norm * min(self.K * self.lam_G * W ** 2, self.L0 + self.L1 * np.sqrt(self.lam_G) * W)

    # geometry helpers
    def diameter(self):
        """diam_{G⁻¹}(Θ_K)."""
        da = 2 * self.a_max; dc = 2 * self.c_max; ds = self.s_max - self.s_min
        n = self.n
        return float(np.sqrt(da ** 2 / self.G_diag[0] + dc ** 2 / self.G_diag[n] + ds ** 2 / self.G_diag[2 * n]))

    def inside(self, P, tol=1e-9):
        n = self.n
        return bool(np.all(np.linalg.norm(P[:, :n], axis=1) <= self.a_max + tol)
                    and np.all(np.linalg.norm(P[:, n:2 * n], axis=1) <= self.c_max + tol)
                    and np.all(P[:, 2 * n] >= self.s_min - tol) and np.all(P[:, 2 * n] <= self.s_max + tol))

    def admissible_ideal(self, P, strict=False):
        """Ideal particles in the inner set Θ_K^ε (closed); strict=True requires its interior (Theorem E)."""
        n, eps = self.n, self.eps
        le = (lambda a, b: a < b) if strict else (lambda a, b: a <= b)
        return bool(np.all(le(np.linalg.norm(P[:, :n], axis=1), self.a_max - eps))
                    and np.all(le(np.linalg.norm(P[:, n:2 * n], axis=1), self.c_max - eps))
                    and np.all(le(self.s_min + eps, P[:, 2 * n])) and np.all(le(P[:, 2 * n], self.s_max - eps)))

    # dynamics
    def control_law(self, x, x_ref, x_ref_dot, P):
        e = x - x_ref
        return x_ref_dot - self.ke * e - f_hat(x, P)

    def raw_velocity(self, x, e, P):
        """G ∇g_t(ϑ_i) (+ leak), before projection.  Returns (vel (m,p), tau or None)."""
        v = self.G_diag[None, :] * grad_g(x, e, P)
        tau = None
        if self.leak is not None:
            if self.prior is None:
                raise ValueError("leak requires a prior configuration")
            if self.leak == 'ot':
                tau, _ = optimal_assignment(P, self.prior, self.Ginv_diag)
                v = v - self.sigma_m * (P - self.prior[tau])
            elif self.leak == 'ot_bary':
                Tbar, gamma, _, _ = barycentric_projection(P, self.prior, self.Ginv_diag)
                tau = gamma
                v = v - self.sigma_m * (P - Tbar)
            else:
                tau = np.arange(self.m)
                v = v - self.sigma_m * (P - self.prior[tau])
        return v, tau

    def project(self, P, v):
        n = self.n
        out = np.empty_like(v)
        for i in range(self.m):
            out[i, :n] = self.proj_a.project_vector(P[i, :n], v[i, :n])
            out[i, n:2 * n] = self.proj_c.project_vector(P[i, n:2 * n], v[i, n:2 * n])
        out[:, 2 * n] = self.proj_s.project_vector(P[:, 2 * n], v[:, 2 * n])
        return out

    def velocity(self, x, e, P):
        v, tau = self.raw_velocity(x, e, P)
        return self.project(P, v), tau

    def hard_clip(self, P):
        n = self.n
        Q = P.copy()
        for i in range(self.m):
            Q[i, :n] = self.proj_a.hard_clip_vector(Q[i, :n])
            Q[i, n:2 * n] = self.proj_c.hard_clip_vector(Q[i, n:2 * n])
        Q[:, 2 * n] = self.proj_s.hard_clip(Q[:, 2 * n])
        return Q

    # Lyapunov quantities (oracle: needs the ideal configuration)
    def ot_distance(self, P, P_star):
        return optimal_assignment(P, P_star, self.Ginv_diag)

    def lyapunov_report(self, x, t, plant, P, P_star, P_star_dot=None):
        """
        Exact evaluation of Theorem B at one state: V_OT, W, active π, V̇_π (chain rule),
        the bound (B), per-particle remainders, cross term.  P_star constant unless P_star_dot given.
        """
        n = self.n
        xr, xrd = plant.x_ref(t), plant.x_ref_dot(t)
        e = x - xr
        u = self.control_law(x, xr, xrd, P)
        edot = plant.f(x, t) + u - xrd
        vel, _ = self.velocity(x, e, P)
        pi, W2 = optimal_assignment(P, P_star, self.Ginv_diag)
        Pm = P_star[pi]
        tilde = P - Pm                                            # matched errors
        Ginv = self.Ginv_diag
        V = 0.5 * float(e @ e) + 0.5 * W2
        Vdot = float(e @ edot) + float(np.einsum('ik,k,ik->', tilde, Ginv, vel)) / self.m
        if P_star_dot is not None:
            Vdot -= float(np.einsum('ik,k,ik->', tilde, Ginv, P_star_dot[pi])) / self.m
        # decomposition pieces
        g_now = (P[:, :n] @ e) * phi_and_grads(x, P)[0]          # g_t(ϑ_i)
        g_star = (Pm[:, :n] @ e) * phi_and_grads(x, Pm)[0]       # g_t(ϑ*_{π(i)})
        gg = grad_g(x, e, P)
        rem = g_star - g_now - np.einsum('ik,ik->i', gg, Pm - P)  # matched Taylor remainders
        cross = float(np.einsum('ik,k,ik->', tilde, Ginv, vel) - np.einsum('ik,ik->', tilde, gg)) / self.m
        eps_vec = plant.f(x, t) - f_hat(x, Pm)
        enorm = float(np.linalg.norm(e))
        sq_each = np.einsum('ik,ik->i', tilde, tilde)
        W = float(np.sqrt(W2))
        bound_B = -self.ke * enorm ** 2 + enorm * float(np.linalg.norm(eps_vec)) + self.excess_bound(enorm, W)
        V_star_t = 0.0
        if P_star_dot is not None:                    # drift term of (B): + V_*(t)·W, V_*(t) = RMS ideal speed at t
            V_star_t = float(np.sqrt(np.einsum('ik,k,ik->', P_star_dot, Ginv, P_star_dot) / self.m))
            bound_B += V_star_t * W
        return dict(t=t, V=V, W=W, Vdot=Vdot, bound_B=bound_B, cross=cross, V_star_t=V_star_t,
                    rem=rem, rem_bound=self.remainder_bound(enorm, sq_each),
                    e_norm=enorm, eps_norm=float(np.linalg.norm(eps_vec)), pi=pi,
                    identity_err=float(np.sqrt(np.einsum('ik,k,ik->', P - P_star, Ginv, P - P_star) / self.m)),
                    inside=self.inside(P),
                    Vdot_identity_check=Vdot - (-self.ke * enorm ** 2 + float(e @ eps_vec) + float(rem.sum()) / self.m + cross
                                                - (float(np.einsum('ik,k,ik->', tilde, Ginv, P_star_dot[pi])) / self.m if P_star_dot is not None else 0.0)))


# ------------------------------------------------------------------ integration
def rk4_closed_loop(plant, ctrl: TransportController, x0, P0, T, dt, record_every=5, callback=None):
    """Fixed-step RK4 on (x, P) with hard clip after every stage.  Returns dict of arrays."""
    def deriv(x, P, t):
        xr, xrd = plant.x_ref(t), plant.x_ref_dot(t)
        e = x - xr
        u = ctrl.control_law(x, xr, xrd, P)
        v, _ = ctrl.velocity(x, e, P)
        return plant.f(x, t) + u, v
    x, P, t = np.array(x0, float), ctrl.hard_clip(np.array(P0, float)), 0.0
    n_steps = int(round(T / dt))
    ts, xs, Ps, es = [], [], [], []
    for k in range(n_steps):
        if k % record_every == 0:
            ts.append(t); xs.append(x.copy()); Ps.append(P.copy()); es.append(np.linalg.norm(x - plant.x_ref(t)))
            if callback is not None:
                callback(x, P, t)
        k1x, k1P = deriv(x, P, t)
        k2x, k2P = deriv(x + 0.5 * dt * k1x, ctrl.hard_clip(P + 0.5 * dt * k1P), t + 0.5 * dt)
        k3x, k3P = deriv(x + 0.5 * dt * k2x, ctrl.hard_clip(P + 0.5 * dt * k2P), t + 0.5 * dt)
        k4x, k4P = deriv(x + dt * k3x, ctrl.hard_clip(P + dt * k3P), t + dt)
        x = x + dt / 6 * (k1x + 2 * k2x + 2 * k3x + k4x)
        P = ctrl.hard_clip(P + dt / 6 * (k1P + 2 * k2P + 2 * k3P + k4P))
        t = (k + 1) * dt
    return dict(t=np.array(ts), x=np.array(xs), P=np.array(Ps), e_norm=np.array(es))


# ------------------------------------------------------------------ representable plant in lifted form
class LiftedRepresentablePlant:
    """f(x) = (1/m)Σ a*_j φ(x; c*_j, s*_j) with a multi-sine reference (ε ≡ 0)."""

    def __init__(self, P_star, amplitude=1.0, freqs=None, phases=None):
        self.P_star = np.array(P_star, float)
        self.n = (self.P_star.shape[1] - 1) // 2
        self.A = amplitude
        self.freqs = np.array(freqs if freqs is not None else 1.0 + 0.3 * np.arange(self.n), float)
        self.phases = np.array(phases if phases is not None else 0.5 * np.arange(self.n), float)

    def f(self, x, t):
        return f_hat(x, self.P_star)

    def x_ref(self, t):
        return self.A * np.sin(self.freqs * t + self.phases)

    def x_ref_dot(self, t):
        return self.A * self.freqs * np.cos(self.freqs * t + self.phases)


class RichReferencePlant(LiftedRepresentablePlant):
    """Same plant, reference = sum of several incommensurate sines per channel (for excitation)."""

    def __init__(self, P_star, amplitude=1.0, n_tones=3):
        super().__init__(P_star, amplitude)
        rng = np.random.default_rng(7)
        self.tones = rng.uniform(0.4, 2.5, size=(self.n, n_tones))
        self.tph = rng.uniform(0, 2 * np.pi, size=(self.n, n_tones))
        self.amps = self.A / n_tones * np.ones((self.n, n_tones))

    def x_ref(self, t):
        return np.sum(self.amps * np.sin(self.tones * t + self.tph), axis=1)

    def x_ref_dot(self, t):
        return np.sum(self.amps * self.tones * np.cos(self.tones * t + self.tph), axis=1)


def sample_lifted_ideal(m, n, ctrl: TransportController, rng, a_scale=1.0, c_spread=1.0, min_sep=0.5):
    """Random admissible ideal configuration (strictly inside the inner box, centres separated)."""
    eps = ctrl.eps
    c = rng.uniform(-c_spread, c_spread, size=(m, n))
    for _ in range(500):
        moved = False
        for i in range(m):
            for j in range(i + 1, m):
                d = c[i] - c[j]; dist = np.linalg.norm(d)
                if dist < min_sep:
                    moved = True; step = 0.05 * (d / dist if dist > 1e-9 else rng.normal(size=n))
                    c[i] += step; c[j] -= step
        if not moved:
            break
    nc = np.linalg.norm(c, axis=1); c[nc > ctrl.c_max - 2 * eps] *= ((ctrl.c_max - 2 * eps) / nc[nc > ctrl.c_max - 2 * eps])[:, None]
    a = rng.normal(0, a_scale, size=(m, n))
    na = np.linalg.norm(a, axis=1); a[na > ctrl.a_max - 2 * eps] *= ((ctrl.a_max - 2 * eps) / na[na > ctrl.a_max - 2 * eps])[:, None]
    s = rng.uniform(ctrl.s_min + 2 * eps, min(ctrl.s_max - 2 * eps, ctrl.s_min + 2 * eps + 1.0), size=m)
    return np.hstack([a, c, s[:, None]])


def perturb_and_permute(P_star, rng, noise, permute=True):
    m = P_star.shape[0]
    pi = rng.permutation(m) if permute else np.arange(m)
    return P_star[pi] + rng.normal(0, noise, size=P_star.shape), pi
