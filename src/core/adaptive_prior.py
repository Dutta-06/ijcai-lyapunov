"""
Theorem D‴.23 (PARADIGM.md §5d): adaptive priors — the OT leak toward a learned prior
====================================================================================

The prior's atoms Z = (z_1, …, z_{K₀}) — its masses q stay fixed — move toward the barycentres of the particles that the plan assigns to them,
while the particles follow the OT leak toward ν_Z := Σ_jq_jδ_{z_j}:

    dϑ_i/ds = G⟨∇_ϑ(r_ΘΦ̃(·;ϑ_i))⟩_ρ − σ₀(ϑ_i − T̄_i),    dz_j/ds = −η(z_j − T̂_j),    T̂_j := (1/q_j)Σ_iγ_{ij}ϑ_i,          (𝔄_m^{OT,η})

γ an optimal plan from μ_Θ to ν_Z and T̄_i = mΣ_jγ_{ij}z_j its barycentric projection. By the envelope theorem
∇_{z_j}W²(μ_Θ, ν_Z) = 2q_jG⁻¹(z_j − T̂_j), so (𝔄_m^{OT,η}) is the gradient flow of the joint energy 𝒥(μ_Θ, Z) := L̄(μ_Θ) + (σ₀/2)W²_{G⁻¹}(μ_Θ, ν_Z)
in the metric (1/m)G⁻¹ per particle and (σ₀q_j/η)G⁻¹ on z_j (η > 0). At a replicated stationary configuration with z_j at the class atom the
leak force vanishes: the replicated stationary points are the critical points of the unregularised K₀-atom loss with weights q, while a
non-replicated one may hold a class apart by the leak toward its barycentre (the K₀-cluster bias of the fast limit). As η → ∞ the prior is
the class barycentres, z_j = T̂_j, and the flow is the gradient flow of L̄ + (σ₀/2)𝒬_q with 𝒬_q(μ) := min_Z W²(μ, ν_Z) the
capacity-constrained K₀-means distortion. The leak toward an adaptive prior is an internal force, so it no longer bounds the parameters.

This module implements the joint drift (a vertex plan frozen inside a step, the prior moving with the particles), the joint energy, the
integrator of `experiments/exp_adaptive_prior.py` (η finite or infinite), and the linear theory at an exact representation: the whitened
weighted Gauss–Newton block H_L, the 2 × 2 collective blocks [[−(h + σ₀), σ₀], [η, −η]] and their rates, the first-order bias of a fixed prior,
and the projections of a K₀-atom deviation on the eigenvectors of H_L.
"""
import numpy as np
from scipy.optimize import linear_sum_assignment

from .ot_leak import _cost_matrix, transport_plan


def plan_to(Th, Z, q, G):
    """A vertex optimal plan γ (m × K₀; row sums 1/m, column sums q) from μ_Θ to ν_Z and W²_{G⁻¹}; integral (the assignment route) when m is
    commensurate with q, the transportation LP otherwise."""
    m, K0 = Th.shape[0], Z.shape[0]
    C = _cost_matrix(Th, Z, G)
    counts = m * q
    if np.allclose(counts, np.round(counts), atol=1e-9):
        n = np.round(counts).astype(int)
        rep = np.repeat(np.arange(K0), n)
        rows, cols = linear_sum_assignment(C[:, rep])
        gamma = np.zeros((m, K0)); gamma[rows, rep[cols]] = 1.0 / m
    else:
        gamma = transport_plan(C, np.full(m, 1.0 / m), q)
    return gamma, float(np.sum(gamma * C))


def barycentres(Th, gamma, q):
    """T̂_j = (1/q_j)Σ_iγ_{ij}ϑ_i for a fixed plan γ, (K₀, 3)."""
    return (gamma.T @ Th) / q[:, None]


def nearest_plan(Th, Z, G):
    """The nearest-atom (cell) plan from μ_Θ to the atoms Z: γ_{ij} = 1/m for the nearest atom j of ϑ_i (G⁻¹ metric), and the distortion
    (1/m)Σ_i min_j‖ϑ_i − z_j‖²_{G⁻¹} — the minimum of W²(μ_Θ, Σ_jq_jδ_{z_j}) over the masses q, attained at the cell masses (semi-relaxed
    transport: with free target masses every particle goes to its nearest atom)."""
    m = Th.shape[0]
    C = _cost_matrix(Th, Z, G)
    lab = np.argmin(C, axis=1)
    gamma = np.zeros((m, Z.shape[0])); gamma[np.arange(m), lab] = 1.0 / m
    return gamma, float(np.mean(C[np.arange(m), lab]))


class AdaptivePriorLoop:
    """(𝔄_m^{OT,η}) on top of an `ot_leak.OTPriorLoop` (whose drift(Th, T) supplies the particles' learning term and the leak toward the targets
    T; its own prior is not used), with the prior mobility η (np.inf: the prior is the class barycentres). plan='nearest' (Proposition D‴.24):
    the nearest-atom leak toward the adaptive prior — the cell plan, the cell masses in place of q (an empty cell keeps its atom), no transit —
    the gradient flow of L̄ + (σ₀/2)(1/m)Σ_i min_j‖ϑ_i − z_j‖²_{G⁻¹} = min over the masses of the OT energy."""

    def __init__(self, ot_loop, eta, plan='ot', box=None):
        assert plan in ('ot', 'nearest')
        self.loop = ot_loop
        self.eta = float(eta)
        self.q = ot_loop.q
        self.G = ot_loop.G
        self.sigma0 = ot_loop.sigma0
        self.plan_kind = plan
        # Proposition D‴.27: the prior's atoms confined to the box 𝒵 = [lo, hi] (unconstrained coordinates, the same for every atom); the
        # prior's flow is then the projected gradient flow of 𝒥 on Θ × 𝒵^{K₀} (its velocity projected onto the tangent cone of the box), and at
        # η = ∞ the prior is the projection of the class barycentres onto the box (the minimiser of W² over 𝒵 for the plan, G being diagonal).
        self.box = None if box is None else (np.asarray(box[0], float), np.asarray(box[1], float))

    def project(self, Z):
        """Z clipped to the box (unchanged without a box)."""
        return Z if self.box is None else np.clip(Z, self.box[0][None, :], self.box[1][None, :])

    def tangent(self, Z, V):
        """The velocity V of the atoms projected onto the tangent cone of the box at Z (unchanged without a box)."""
        if self.box is None:
            return V
        lo, hi = self.box
        V = np.array(V, float)
        V[(Z <= lo[None, :] + 1e-14) & (V < 0)] = 0.0
        V[(Z >= hi[None, :] - 1e-14) & (V > 0)] = 0.0
        return V

    def plan(self, Th, Z):
        """(γ, the transport cost): the optimal plan to ν_Z with the masses q (OT) or the cell plan (nearest)."""
        return plan_to(Th, Z, self.q, self.G) if self.plan_kind == 'ot' else nearest_plan(Th, Z, self.G)

    def bary(self, Th, gamma, Z):
        """The class barycentres of the plan (OT: divided by q; nearest: by the cell masses, an empty cell keeping z_j)."""
        if self.plan_kind == 'ot':
            return barycentres(Th, gamma, self.q)
        w = gamma.sum(axis=0); out = np.array(Z, float).copy(); nz = w > 0
        out[nz] = (gamma.T @ Th)[nz] / w[nz, None]
        return out

    def drift(self, Th, Z, gamma, tangent=True):
        """(dΘ/ds, dZ/ds) for the plan γ (finite η: with a box, the prior's velocity projected onto the tangent cone of the box at Z — the
        projected flow — unless tangent=False); for η = ∞ the prior is replaced by the barycentres and dZ/ds is not used."""
        m = Th.shape[0]
        if np.isinf(self.eta):
            Zb = self.project(self.bary(Th, gamma, Z))
            return self.loop.drift(Th, m * (gamma @ Zb)), np.zeros_like(Z)
        V = -self.eta * (Z - self.bary(Th, gamma, Z))
        return self.loop.drift(Th, m * (gamma @ Z)), (self.tangent(Z, V) if tangent else V)

    def energy(self, Th, Z):
        """𝒥(μ_Θ, Z) = L̄(μ_Θ) + (σ₀/2)W²(μ_Θ, ν_Z) with an optimal plan (for η = ∞ evaluate it at Z = the barycentres of an optimal plan)."""
        return self.loop.L_bar(Th) + 0.5 * self.sigma0 * self.plan(Th, Z)[1]

    def in_transit(self, Th, Z, gamma):
        """True if some particle is assigned by the plan to an atom of ν_Z that is not its nearest one (the mass in transit; never for the
        cell plan)."""
        if self.plan_kind == 'nearest':
            return False
        lab = np.argmax(gamma, axis=1)
        return bool(np.any(lab != np.argmin(_cost_matrix(Th, Z, self.G), axis=1)))

    def fixed_prior_drift(self, Th, Z, gamma):
        """The largest whitened particle drift of the fixed-prior OT leak toward ν_Z with the plan γ."""
        d = self.loop.drift(Th, Th.shape[0] * (gamma @ Z))
        return float(np.sqrt(np.sum(d ** 2 / self.G[None, :], axis=1)).max())

    def run(self, Th, Z, S_end, dt, plan_every=10, rec_every=250, callback=None, gate=False, settle_tol=1e-5):
        """Midpoint rule in slow time; the plan to ν_Z is recomputed every `plan_every` steps and frozen in between (the vertex-selection
        integrator), the prior moving with the particles inside a step; for η = ∞ the prior is set to the class barycentres of the plan.
        gate = True (or 'transit') freezes the prior while mass is in transit (checked at every plan refresh): the class barycentres of the
        plan then contain particles of other clusters, and an adaptive prior would be dragged toward them. gate = 'settled' also waits, the
        first time, until the fixed-prior drift is below settle_tol, so that the adaptation starts from the fixed-prior stationary point.
        `callback(s, Th, Z, adapting)` at the records. With a box (Proposition D‴.27) the step projects positions, not velocities:
        Z_mid = Π(Z + (dt/2)v(Z)), Z_new = Π(Z + dt·v(Z_mid)), v the unprojected velocity of the prior; its fixed points on a face are the
        projected flow's stationary points. (Projecting the velocity onto the tangent cone at Z_mid instead stalls the prior just inside a
        face whenever the midpoint lies past it: referee round 19, D3n-2.) Without a box the step is the plain midpoint rule. Returns
        (Θ, Z, the first time the prior adapted)."""
        Th = np.array(Th, float); Z = np.array(Z, float)
        n = int(round(S_end / dt))
        gamma, adapting, t_open = None, True, float('nan')
        eta = self.eta
        settled = gate != 'settled'
        for k in range(n + 1):
            if k % plan_every == 0:
                gamma, _ = self.plan(Th, Z)
                free = not self.in_transit(Th, Z, gamma)
                if free and not settled and self.fixed_prior_drift(Th, Z, gamma) < settle_tol:
                    settled = True
                adapting = (not gate) or (free and settled)
                if adapting and not np.isfinite(t_open):
                    t_open = k * dt
                if adapting and np.isinf(eta):
                    Z = self.project(self.bary(Th, gamma, Z))
                    gamma, _ = self.plan(Th, Z)
            if callback is not None and k % rec_every == 0:
                callback(k * dt, Th, Z, adapting)
            if k == n:
                break
            self.eta = eta if adapting else 0.0
            if adapting or not np.isinf(eta):
                k1T, k1Z = self.drift(Th, Z, gamma, tangent=False)
                k2T, k2Z = self.drift(Th + 0.5 * dt * k1T, self.project(Z + 0.5 * dt * k1Z), gamma, tangent=False)
            else:                                                          # η = ∞ but gated: the fixed-prior OT leak
                m = Th.shape[0]
                k1T = self.loop.drift(Th, m * (gamma @ Z)); k1Z = np.zeros_like(Z)
                k2T = self.loop.drift(Th + 0.5 * dt * k1T, m * (gamma @ Z)); k2Z = np.zeros_like(Z)
            self.eta = eta
            Th, Z = Th + dt * k2T, self.project(Z + dt * k2Z)
            if adapting and np.isinf(eta):
                Z = self.project(self.bary(Th, gamma, Z))
        return Th, Z, t_open


# ---------------------------------------------------------------------- the linear theory at an exact representation
def gauss_newton_block(loop, atoms, q, h=1e-6):
    """H_L: the whitened, weighted Hessian of the loss part of the class potentials at the K₀-atom measure Σ_kq_kδ_{atoms_k} —
    diag(√q ⊗ √G)·∂(∇V_L[μ](ϑ_k))/∂ϑ_l·diag(√G/√q), symmetrised (central differences). At an exact representation (zero residual) it is the
    Gauss–Newton matrix [√(q_kq_l)G^{1/2}⟨∇Φ̃_kᵀ∇Φ̃_l⟩_ρG^{1/2}]. Returns (H_L, eigenvalues ascending, eigenvectors as columns)."""
    K = atoms.shape[0]
    d = 3 * K
    sqG = np.sqrt(loop.G)

    def grad_loss(V):
        Th = V.reshape(K, 3)
        a, c, s, diff, phi = loop.features(Th)
        r = loop.F - phi @ (q * a)
        return loop.grad_V_L(r, Th).ravel()

    v = np.asarray(atoms, float).ravel()
    Jm = np.zeros((d, d))
    for j in range(d):
        e = np.zeros(d); e[j] = h
        Jm[:, j] = (grad_loss(v + e) - grad_loss(v - e)) / (2 * h)
    sw = np.repeat(np.sqrt(q), 3); sg = np.tile(sqG, K)
    HL = (sw * sg)[:, None] * Jm * (sg / sw)[None, :]
    HL = 0.5 * (HL + HL.T)
    ev, V = np.linalg.eigh(HL)
    return HL, ev, V


def pair_rates(h, sigma0, eta):
    """The rates (λ⁻, λ⁺) of the collective pair [[−(h + σ₀), σ₀], [η, −η]]: the roots of λ² − (h + σ₀ + η)λ + ηh = 0; for η = ∞, (h, ∞)."""
    h = np.asarray(h, float)
    if np.isinf(eta):
        return h, np.full_like(h, np.inf)
    tr = h + sigma0 + eta
    disc = np.sqrt(np.maximum(tr ** 2 - 4 * eta * h, 0.0))
    lm = 2 * eta * h / (tr + disc)                                             # the smaller root, written stably
    return lm, tr - lm


def whitened(dev, q, G):
    """y = (√q_k G^{−1/2}dev_k)_k stacked, the whitened weighted coordinates of a K₀-atom deviation (K₀, 3) → (3K₀,)."""
    return (np.sqrt(q)[:, None] * dev / np.sqrt(G)[None, :]).ravel()


def pair_solution(h, sigma0, eta, y0, z0, s):
    """The exact solution of the pair ẏ = −(h + σ₀)y + σ₀z, ż = η(y − z) (η finite) or of ẏ = −hy (η = ∞, z = y) at times s."""
    s = np.asarray(s, float)
    if np.isinf(eta):
        return y0 * np.exp(-h * s), y0 * np.exp(-h * s)
    A = np.array([[-(h + sigma0), sigma0], [eta, -eta]])
    w, P = np.linalg.eig(A)
    c = np.linalg.solve(P, np.array([y0, z0]))
    sol = (P[:, None, :] * (c[None, None, :] * np.exp(w[None, None, :] * s[None, :, None]))).sum(axis=2)
    return sol[0].real, sol[1].real
