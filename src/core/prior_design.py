"""
Proposition D‴.9 (PARADIGM.md §5d): the leak toward a measure-valued prior, and the design of the window
=====================================================================================================

The σ-modification leaks every kernel toward one prior ϑ⁰. With a prior *measure* ν = Σ_j q_jδ_{ϑ⁰_j} the natural leak is the
G-gradient of the distance-to-support potential

    V_ν(ϑ) := ½ min_j ‖ϑ − ϑ⁰_j‖²_{G⁻¹},        −G∇V_ν(ϑ) = −(ϑ − ϑ⁰_{j(ϑ)}),    j(ϑ) := argmin_j‖ϑ − ϑ⁰_j‖_{G⁻¹}

(the nearest-atom leak; the barycentric OT leak of Theorem C′ is its mass-constrained relative).  Inside each Voronoi cell of ν
the leak is the ordinary σ-modification toward that cell's atom, so Theorem D‴, Proposition D‴.5 and Corollary D‴.7 hold
verbatim for the averaged slow loop

    dϑ_i/ds = G⟨∇_ϑ(r_Θ Φ̃(·;ϑ_i))⟩_ρ − σ₀(ϑ_i − ϑ⁰_{j(ϑ_i)})                                                              (𝔄_m^ν)

whenever the ϱ′-neighbourhood of the support of the mean-field flow stays inside the cells (the assignment is locally
constant).  What changes is the stationary residual: a prior whose atoms sit near the plant's structure removes the leak's
bias, lowers κ_R[μ*] (the window's lower end) and raises the margin κ, and *raises* the weight sensitivity S (the floor of
Proposition D‴.5): the stationary atoms carry the plant's full amplitudes, to which the sensitivity is proportional — bias
against variance.  This module implements (𝔄_m^ν) as a subclass of `chaos_learning.AveragedSlowLoop` and the design table of
`experiments/exp_prior_design.py`.
"""
import numpy as np

from .chaos_learning import AveragedSlowLoop, _sig


class MixturePriorLoop(AveragedSlowLoop):
    """(𝔄_m^ν): the averaged slow loop with the nearest-atom leak toward the prior atoms `prior` (K₀ × 3, unconstrained
    coordinates).  With a single atom it is `AveragedSlowLoop` with theta0 = prior[0]."""

    def __init__(self, net, f_vec, x_ref, period, G, sigma0, prior, n_quad=256):
        prior = np.atleast_2d(np.asarray(prior, float))
        super().__init__(net, f_vec, x_ref, period, G, sigma0, theta0=prior[0], n_quad=n_quad)
        self.prior = prior

    # ------------------------------------------------------------------ the leak
    def assignment(self, Th):
        """j(ϑ_i): the nearest prior atom in the G⁻¹ metric, (m,)."""
        d2 = np.sum((Th[:, None, :] - self.prior[None, :, :]) ** 2 / self.G[None, None, :], axis=2)
        return np.argmin(d2, axis=1)

    def leak_target(self, Th):
        """ϑ⁰_{j(ϑ_i)}, (m, 3)."""
        return self.prior[self.assignment(Th)]

    def leak_potential(self, Th):
        """V_ν(ϑ_i) = ½ min_j‖ϑ_i − ϑ⁰_j‖²_{G⁻¹}, (m,)."""
        d2 = np.sum((Th[:, None, :] - self.prior[None, :, :]) ** 2 / self.G[None, None, :], axis=2)
        return 0.5 * d2.min(axis=1)

    def cell_margin(self, Th):
        """For every row, the gap between the distances to the nearest and the second-nearest prior atom (G⁻¹ metric); a
        gap above 2ϱ′ on the support is sufficient for the 'assignment locally constant' hypothesis of Proposition D‴.9 (the G⁻¹-distance
        to the nearest Voronoi boundary is at least half the gap, with equality on the segment between the two atoms)."""
        d = np.sqrt(np.sum((Th[:, None, :] - self.prior[None, :, :]) ** 2 / self.G[None, None, :], axis=2))
        d.sort(axis=1)
        return d[:, 1] - d[:, 0] if d.shape[1] > 1 else np.full(Th.shape[0], np.inf)

    # ------------------------------------------------------------------ the objective and its gradient
    def J(self, Th):
        """J_m^ν(Θ) = mL̄(μ_Θ) + σ₀Σ_iV_ν(ϑ_i)."""
        return Th.shape[0] * self.L_bar(Th) + self.sigma0 * float(np.sum(self.leak_potential(Th)))

    def grad_J(self, Th):
        """∇_{ϑ_i}J_m^ν = ∇V_L[μ_Θ](ϑ_i) + σ₀G⁻¹(ϑ_i − ϑ⁰_{j(ϑ_i)})  (away from the Voronoi boundaries)."""
        return self.grad_V_L(self.residual(Th), Th) + self.sigma0 * (Th - self.leak_target(Th)) / self.G[None, :]

    def drift(self, Th):
        a, c, s, diff, phi = self.features(Th)
        m = Th.shape[0]
        r = self.F - phi @ a / m
        net = self.net
        da = 1.0 - (a / net.a_max) ** 2
        dc = 1.0 - (c / net.c_max) ** 2
        sk = _sig(Th[:, 2]); ds = (net.s_max - net.s_min) * sk * (1 - sk)
        g = np.empty_like(Th)
        g[:, 0] = (phi.T @ r) * da
        g[:, 1] = a * dc / s ** 2 * ((phi * diff).T @ r)
        g[:, 2] = a * ds / s ** 3 * ((phi * diff ** 2).T @ r)
        return self.G[None, :] * g / self.N - self.sigma0 * (Th - self.leak_target(Th))

    # ------------------------------------------------------------------ the three terms (the leak of the reference's cell)
    def decomposition(self, Th, Thb, fhat_mf=None):
        """As in AveragedSlowLoop.decomposition, with the leak of the REFERENCE's assignment in T_A (the proof's convention when
        the assignment is locally constant: both ϑ_i and ϑ̄_i lie in the cell of ϑ̄_i)."""
        d = Th - Thb
        rT, rB = self.residual(Th), self.residual(Thb)
        r_mf = rB if fhat_mf is None else self.F - np.asarray(fhat_mf, float)
        D = float(np.mean(np.sum(d ** 2 / self.G[None, :], axis=1)))
        tgt = self.leak_target(Thb)
        gA = self.grad_V_L(rB, Th) - self.grad_V_L(rB, Thb) + self.sigma0 * d / self.G[None, :]
        T_A = -float(np.mean(np.sum(d * gA, axis=1)))
        gB = self.grad_V_L(rT, Th) - self.grad_V_L(rB, Th)
        T_B = -float(np.mean(np.sum(d * gB, axis=1)))
        gC = self.grad_V_L(rB, Thb) - self.grad_V_L(r_mf, Thb)
        T_C = -float(np.mean(np.sum(d * gC, axis=1)))
        drift_m = -self.G[None, :] * (self.grad_V_L(rT, Th) + self.sigma0 * (Th - tgt) / self.G[None, :])
        drift_mf = -self.G[None, :] * (self.grad_V_L(r_mf, Thb) + self.sigma0 * (Thb - tgt) / self.G[None, :])
        half_dD = float(np.mean(np.sum(d * (drift_m - drift_mf) / self.G[None, :], axis=1)))
        h = self.f_hat(Th) - self.f_hat(Thb)
        return dict(T_A=T_A, T_B=T_B, T_C=T_C, total=T_A + T_B + T_C, half_dD=half_dD, D=D,
                    h_norm=float(np.sqrt(np.mean(h ** 2))), bound_B_coeff=1.0 / 16.0, same_cell=bool(np.all(self.assignment(Th) == self.assignment(Thb))))


# ---------------------------------------------------------------------- Proposition D‴.5 with a measure-valued prior
def weighted_energy_grad_prior(loop, w, Th):
    """∇_ϑ𝒱^ν[μ_w](ϑ_k) for the K-atom measure with weights w and the nearest-atom leak: ∇V_L[μ_w](ϑ_k) + σ₀G⁻¹(ϑ_k − ϑ⁰_{j(ϑ_k)})."""
    a, c, s, diff, phi = loop.features(Th)
    r = loop.F - phi @ (np.asarray(w, float) * a)
    return loop.grad_V_L(r, Th) + loop.sigma0 * (Th - loop.leak_target(Th)) / loop.G[None, :]


def weighted_stationary_atoms_prior(loop, w, Th0, iters=100, h=1e-6, tol=1e-13, w0=None, n_cont=1):
    """Newton with backtracking (and optional continuation in the weights) on ∇𝒱^ν[μ_w](ϑ_k) = 0; returns (Θ*(w), residual)."""
    K = Th0.shape[0]; d = 3 * K
    v = Th0.ravel().copy()
    w = np.asarray(w, float)
    path = [w] if w0 is None else [np.asarray(w0, float) + (w - np.asarray(w0, float)) * (j + 1) / n_cont for j in range(n_cont)]
    for wj in path:
        F = lambda v, wj=wj: weighted_energy_grad_prior(loop, wj, v.reshape(K, 3)).ravel()
        for _ in range(iters):
            g = F(v); gn = np.abs(g).max()
            if gn < tol:
                break
            J = np.zeros((d, d))
            for j in range(d):
                e = np.zeros(d); e[j] = h
                J[:, j] = (F(v + e) - F(v - e)) / (2 * h)
            step = np.linalg.solve(J, g)
            lam = 1.0
            while lam > 1e-6:
                vn = v - lam * step
                if np.abs(F(vn)).max() < gn:
                    break
                lam *= 0.5
            v = vn
    return v.reshape(K, 3), float(np.abs(weighted_energy_grad_prior(loop, w, v.reshape(K, 3))).max())


def weight_sensitivity_prior(loop, w, atoms, h=1e-4):
    """S[k, :, j] = ∂ϑ*_k/∂w_j for the nearest-atom leak (central differences of the weighted stationary point), (K, 3, K)."""
    K = atoms.shape[0]
    S = np.zeros((K, 3, K))
    for j in range(K):
        wp = np.array(w, float); wp[j] += h
        wm = np.array(w, float); wm[j] -= h
        Tp, rp = weighted_stationary_atoms_prior(loop, wp, atoms)
        Tm, rm = weighted_stationary_atoms_prior(loop, wm, atoms)
        assert rp < 1e-9 and rm < 1e-9
        S[:, :, j] = (Tp - Tm) / (2 * h)
    return S


# ---------------------------------------------------------------------- the function-level law
def function_sensitivity(loop, w, atoms, S):
    """F[:, j] = ∂f̂_w/∂w_j on the quadrature grid, for f̂_w(x) = Σ_k w_kΦ̃(x;ϑ*_k(w)):
    ∂f̂/∂w_j = Φ̃(·;ϑ*_j) + Σ_k w_k ∇_ϑΦ̃(·;ϑ*_k)·S[k, :, j]   (N, K)."""
    K = atoms.shape[0]
    a, c, s, diff, phi = loop.features(atoms)                            # phi (N, K) without the amplitude
    Phi = phi * a[None, :]                                                # Φ̃(x; ϑ*_k)
    net = loop.net
    da = 1.0 - (a / net.a_max) ** 2
    dc = 1.0 - (c / net.c_max) ** 2
    sk = _sig(atoms[:, 2]); ds = (net.s_max - net.s_min) * sk * (1 - sk)
    # ∇_ϑΦ̃(x;ϑ_k) in the unconstrained coordinates, (N, K, 3)
    gradPhi = np.stack([phi * da[None, :], Phi * diff / s[None, :] ** 2 * dc[None, :], Phi * diff ** 2 / s[None, :] ** 3 * ds[None, :]], axis=2)
    F = np.empty((loop.N, K))
    for j in range(K):
        F[:, j] = Phi[:, j] + np.einsum('k,nkc,kc->n', np.asarray(w, float), gradPhi, S[:, :, j])
    return F


def predicted_stationary_function_D(F, wbar, m, M):
    """E‖f̂_m − f̂_M‖²_{L²(ρ)} ≈ (1/m − 1/M)·tr(F C̃ Fᵀ)/N with C̃ = diag(w̄) − w̄w̄ᵀ (the function-level form of Proposition D‴.5(ii))."""
    wbar = np.asarray(wbar, float)
    Ct = np.diag(wbar) - np.outer(wbar, wbar)
    return (1.0 / m - 1.0 / M) * float(np.trace(F @ Ct @ F.T)) / F.shape[0]


def splitting_margin_prior(loop, w, atoms, wbar=None, n_cont=4):
    """The splitting margin σ₀ − max_k κ_R of the replicated weighted stationary configuration under the nearest-atom leak (the
    analogue of chaos_learning.splitting_margin); Θ*(w) by continuation from wbar.  Returns (margin, Θ*(w), κ_R at the atoms)."""
    K = atoms.shape[0]
    wbar = np.ones(K) / K if wbar is None else np.asarray(wbar, float)
    Tw, res = weighted_stationary_atoms_prior(loop, w, atoms, w0=wbar, n_cont=n_cont)
    a, c, s, diff, phi = loop.features(Tw)
    r = loop.F - phi @ (np.asarray(w, float) * a)
    k = loop.kappa_R(r, Tw)
    return float(loop.sigma0 - k.max()), Tw, k
