"""
Frozen-measure averaging of the adaptive closed loop (PARADIGM.md §6c, Proposition H.7)
=======================================================================================

At a FROZEN parameter measure μ = Σ_j w_j δ_{ϑ_j} (the weights may be any reals, so that first and second
variations can be taken by perturbing a weight) the closed loop

    ẋ = f(x) + ẋ_m − k_e e − f̂_μ(x),    e = x − x_m(t),    f̂_μ = Σ_j w_j Φ̃(·; ϑ_j),

reduces to the error equation  ė = r_μ(x_m + e) − k_e e,  r_μ := f − f̂_μ.  For a T₀-periodic reference and
k_e large it has a unique periodic solution e_μ NEAR 0 (globally unique only under a one-sided Lipschitz condition
on f; PARADIGM.md H.7(a)), and the period average of the adaptive drift
G∇_ϑ[eᵀΦ̃(x;ϑ)] is G∇_ϑ P_μ(ϑ), with the AVERAGED POTENTIAL

    P_μ(ϑ) := ⟨ e_μ(t)ᵀ Φ̃(x_m(t) + e_μ(t); ϑ) ⟩_t .

Proposition H.7 states  P_μ = (1/k_e)⟨r_μᵀΦ̃_ϑ⟩_ρ + (1/k_e²)(U_even + U_odd) + O(k_e⁻³)  with

    U_even[μ](ϑ) = ⟨ (r_μ·∇_x)(r_μᵀΦ̃_ϑ) ⟩_ρ             (state shift x = x_m + e; a functional of ρ)
    U_odd[μ](ϑ)  = −⟨ (d/dt r_μ(x_m(t)))ᵀ Φ̃_ϑ(x_m(t)) ⟩_t   (filter lag; changes sign with the orbit)

and U_odd = −⟨ḟᵀΦ̃_ϑ⟩ + ∫K(ϑ,ϑ')dμ(ϑ'), where K(ϑ,ϑ') = ⟨Φ̃(x_m;ϑ)ᵀ (d/dt)Φ̃(x_m;ϑ')⟩_t is the
antisymmetric circulation kernel.

Autonomous plants only: the caller supplies f and its Jacobian as vectorised callables
f_vec(X) -> (N, n) and jac_vec(X) -> (N, n, n) with jac[t, k, l] = ∂f_k/∂x_l.
The features are those of langevin.SaturatedLiftedNetwork.
"""
import numpy as np


class Orbit:
    """A T₀-periodic reference sampled on N grid points (for averages) and on an RK4 grid of N·sub steps."""

    def __init__(self, x_ref, x_ref_dot, period, N=2000, sub=2):
        self.period, self.N, self.sub = float(period), int(N), int(sub)
        self.ts = np.arange(self.N) * self.period / self.N
        self.xm = np.array([x_ref(t) for t in self.ts], dtype=float)
        self.xmd = np.array([x_ref_dot(t) for t in self.ts], dtype=float)
        M = self.N * self.sub
        self.dt = self.period / M
        self.fine = np.array([x_ref(k * self.dt) for k in range(M + 1)], dtype=float)
        self.half = np.array([x_ref((k + 0.5) * self.dt) for k in range(M)], dtype=float)


def features(net, X, Th):
    """Φ̃(X; ϑ_j) as (N, k, n) and ∂_xΦ̃ as (N, k, n, n) with [t, j, k, l] = ∂Φ̃_k/∂x_l."""
    a, c, s = net.physical(np.atleast_2d(Th))
    D = X[:, None, :] - c[None, :, :]                                        # (N, k, n)
    ph = np.exp(-np.sum(D ** 2, axis=2) / (2 * s[None, :] ** 2))             # (N, k)
    Phi = ph[:, :, None] * a[None, :, :]
    dPhi = -(ph[:, :, None, None] * a[None, :, :, None] * D[:, :, None, :]) / (s[None, :, None, None] ** 2)
    return Phi, dPhi


class FrozenMeasure:
    """μ = Σ_j w_j δ_{ϑ_j}; f̂_μ(x) = Σ_j w_j Φ̃(x; ϑ_j)."""

    def __init__(self, net, Th, w=None):
        self.net = net
        self.Th = np.atleast_2d(np.asarray(Th, float))
        m = len(self.Th)
        self.w = np.full(m, 1.0 / m) if w is None else np.asarray(w, float)

    def with_atom(self, theta, weight):
        return FrozenMeasure(self.net, np.vstack([self.Th, np.atleast_2d(theta)]), np.append(self.w, weight))

    def f_hat(self, X):
        Phi, _ = features(self.net, X, self.Th)
        return np.einsum('j,tjn->tn', self.w, Phi)

    def jac(self, X):
        _, dPhi = features(self.net, X, self.Th)
        return np.einsum('j,tjkl->tkl', self.w, dPhi)


def periodic_error(f_vec, fm, orbit, ke, n_periods=2):
    """The periodic solution e_μ near 0 of ė = r_μ(x_m + e) − k_e e on the grid orbit.ts (RK4 from e = 0, which lies
    in its basin for k_e large; the transient decays like exp(−(k_e − Lip r_μ)T₀), so the last of n_periods periods is
    periodic to machine precision)."""
    n = orbit.xm.shape[1]
    e = np.zeros(n)

    def r(x):
        X = x[None, :]
        return f_vec(X)[0] - fm.f_hat(X)[0]

    M, dt = orbit.N * orbit.sub, orbit.dt
    out = np.zeros((orbit.N, n))
    for rep in range(n_periods):
        last = rep == n_periods - 1
        for k in range(M):
            if last and k % orbit.sub == 0:
                out[k // orbit.sub] = e
            x0, xh, x1 = orbit.fine[k], orbit.half[k], orbit.fine[k + 1]
            k1 = r(x0 + e) - ke * e
            k2 = r(xh + e + 0.5 * dt * k1) - ke * (e + 0.5 * dt * k1)
            k3 = r(xh + e + 0.5 * dt * k2) - ke * (e + 0.5 * dt * k2)
            k4 = r(x1 + e + dt * k3) - ke * (e + dt * k3)
            e = e + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return out


def averaged_potential(f_vec, fm, orbit, ke, Th_eval):
    """P_μ(ϑ) = ⟨e_μᵀ Φ̃(x_m + e_μ; ϑ)⟩ at each row of Th_eval (one periodic solve for all rows)."""
    e = periodic_error(f_vec, fm, orbit, ke)
    Phi, _ = features(fm.net, orbit.xm + e, Th_eval)
    return np.einsum('tn,tjn->j', e, Phi) / orbit.N


def leading_potential(f_vec, fm, orbit, ke, Th_eval):
    """The (𝔄) potential (1/k_e)⟨r_μᵀΦ̃_ϑ⟩_ρ = −(1/k_e) V_L[μ](ϑ)."""
    r = f_vec(orbit.xm) - fm.f_hat(orbit.xm)
    Phi, _ = features(fm.net, orbit.xm, Th_eval)
    return np.einsum('tn,tjn->j', r, Phi) / (orbit.N * ke)


def expansion_terms(f_vec, jac_vec, fm, orbit, Th_eval):
    """(U_even, U_odd) of Proposition H.7 at each row of Th_eval."""
    X, Xd = orbit.xm, orbit.xmd
    r = f_vec(X) - fm.f_hat(X)
    Jr = jac_vec(X) - fm.jac(X)                                              # ∂_l r_k
    rdot = np.einsum('tkl,tl->tk', Jr, Xd)
    Phi, dPhi = features(fm.net, X, Th_eval)
    U_even = (np.einsum('tkl,tl,tjk->j', Jr, r, Phi) + np.einsum('tk,tjkl,tl->j', r, dPhi, r)) / orbit.N
    U_odd = -np.einsum('tk,tjk->j', rdot, Phi) / orbit.N
    return U_even, U_odd


def circulation_kernel(net, orbit, ThA, ThB):
    """K[i, j] = ⟨Φ̃(x_m; ϑ_i)ᵀ (d/dt) Φ̃(x_m; ϑ'_j)⟩_t (loop integral (1/T₀)∮ Φ̃_ϑᵀ dΦ̃_ϑ')."""
    PhiA, _ = features(net, orbit.xm, ThA)
    _, dPhiB = features(net, orbit.xm, ThB)
    dB = np.einsum('tjkl,tl->tjk', dPhiB, orbit.xmd)
    return np.einsum('tik,tjk->ij', PhiA, dB) / orbit.N


def second_variation_formula(f_vec, jac_vec, fm, orbit, Th):
    """Predicted k_e²-coefficient of the second variation δP_μ(ϑ_i)/δμ(ϑ_j): S_even + S_odd, with S_odd = K."""
    X = orbit.xm
    r = f_vec(X) - fm.f_hat(X)
    Jr = jac_vec(X) - fm.jac(X)
    Phi, dPhi = features(fm.net, X, Th)
    grad_rPhi = np.einsum('tkl,tik->til', Jr, Phi) + np.einsum('tk,tikl->til', r, dPhi)    # ∇_x(rᵀΦ̃_i)
    term1 = np.einsum('tjl,til->ij', Phi, grad_rPhi)                                       # (Φ̃_j·∇)(rᵀΦ̃_i)
    term2 = np.einsum('tl,tjkl,tik->ij', r, dPhi, Phi) + np.einsum('tl,tjk,tikl->ij', r, Phi, dPhi)  # (r·∇)(Φ̃_jᵀΦ̃_i)
    S_even = -(term1 + term2) / orbit.N
    return S_even, circulation_kernel(fm.net, orbit, Th, Th)


def exact_second_variation(f_vec, fm, orbit, ke, Th, h=1e-3):
    """S[i, j] = d/dε P_{μ + ε δ_{ϑ_j}}(ϑ_i) by central differences on the exact averaged potential."""
    k = len(Th)
    S = np.zeros((k, k))
    for j in range(k):
        Pp = averaged_potential(f_vec, fm.with_atom(Th[j], +h), orbit, ke, Th)
        Pm = averaged_potential(f_vec, fm.with_atom(Th[j], -h), orbit, ke, Th)
        S[:, j] = (Pp - Pm) / (2 * h)
    return S


def cyclic_defect(S):
    """Cyclic sum of the antisymmetric part over the first three points: zero for every triple iff the
    antisymmetric part is a coboundary a(ϑ) − a(ϑ'), i.e. iff the interaction is reciprocal (a first variation)."""
    A = S - S.T
    return float(A[0, 1] + A[1, 2] + A[2, 0])


def oscillation_correction(net, f_vec, orbit, Th, G_diag, rel_step=1e-5):
    """
    Second-order averaging term of the m-particle closed loop (coefficient of 1/k_e²), with weights 1/m:

        ⟨ D_Θ H(Θ, t) · u(Θ, t) ⟩_t ,   H_i(Θ, t) := G ∇_ϑ[ r_Θ(x_m(t))ᵀ Φ̃(x_m(t); ϑ_i) ],

    where u is the zero-mean antiderivative in t of H − ⟨H⟩ (so that Θ = Y + u/k_e + O(k_e⁻²) and the
    period-mean of Θ is Y + O(k_e⁻²)).  This term comes from the intra-period oscillation of the parameters,
    which the frozen-measure expansion ignores; it is odd under reversal of the orbit.  Returns (m, p).
    The leak −σ(ϑ − ϑ⁰) is constant in t and does not contribute.
    """
    Th = np.atleast_2d(np.asarray(Th, float))
    m, p = Th.shape
    G = np.asarray(G_diag, float)
    X, N = orbit.xm, orbit.N

    def H_at(Thx, k):
        x = X[k]
        r = f_vec(x[None, :])[0] - net.f_hat(x, Thx)
        return G[None, :] * net.grad_g(x, r, Thx)

    Hs = np.array([H_at(Th, k) for k in range(N)])                    # (N, m, p)
    F = np.fft.rfft(Hs - Hs.mean(0), axis=0)
    kk = np.arange(F.shape[0])
    U = np.zeros_like(F)
    U[1:] = F[1:] / (1j * kk[1:, None, None] * (2 * np.pi / orbit.period))
    u = np.fft.irfft(U, n=N, axis=0)                                  # zero-mean antiderivative
    out = np.zeros((m, p))
    for k in range(N):
        uk = u[k]
        d = rel_step / (np.sqrt(np.sum(uk ** 2)) + 1e-30)
        out += (H_at(Th + d * uk, k) - H_at(Th - d * uk, k)) / (2 * d)
    return out / N
