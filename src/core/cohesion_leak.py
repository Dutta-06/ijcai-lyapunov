"""
Proposition D‴.32 (PARADIGM.md §5d): the cohesion leak — the leak on the unexcited subspace, without a prior
==========================================================================================================

The particles attract each other through a radial interaction of range R (in the G⁻¹ metric) instead of leaking toward a prior:

    dϑ_i/ds = G⟨∇_ϑ(r_ΘΦ̃(·;ϑ_i))⟩_ρ − (σ₀/m)Σ_jψ(‖ϑ_i − ϑ_j‖²_{G⁻¹}/R²)(ϑ_i − ϑ_j),                          (𝔄_m^coh)

with ψ(z) = (1 − z)³₊ (kind 'compact': the force vanishes beyond the distance R) or ψ(z) = e^{−z/2} (kind 'gauss': the scale R), ψ(0) = 1.
It is the G-gradient flow of

    J^coh_m(Θ) := mL̄(μ_Θ) + (σ₀m/2)∫∫W dμ_Θdμ_Θ,    W(u) = R²ω(‖u‖²_{G⁻¹}/R²),    ω(z) = (1 − (1 − z)⁴₊)/8  resp.  1 − e^{−z/2},

since G∇W(u) = ψ(‖u‖²_{G⁻¹}/R²)u (ω′ = ψ/2) and ∂_{ϑ_i}[(σ₀/(2m))Σ_{j,l}W(ϑ_j − ϑ_l)] = (σ₀/m)Σ_j∇W(ϑ_i − ϑ_j); ∇²W(0) = G⁻¹. The force on a
particle is σ₀ times the local kernel mass (1/m)Σ_jψ_ij times the mean-shift vector Σ_jψ_ijϑ_j/Σ_jψ_ij − ϑ_i: each particle climbs the
kernel density of the cloud (the shadow kernel of ψ), so the leak pulls the network toward the modes of its own parameter distribution — no
prior, no atom count, no plan.

At a configuration replicated by classes whose atoms are pairwise at least R apart ('compact') the interaction force vanishes — inside a
class ψ(0)·0, across classes ψ = 0 — so such a configuration is stationary iff its atoms are a critical point of the unregularised weighted
K-atom loss L_w (w the class masses), for every K and every w: no leak bias, the property of the adaptive prior (Theorem D‴.23,
Proposition D‴.24) without the prior. The linearisation there acts on the separations of class k as G(R_k − σ₀w_kG⁻¹) (whitened:
−(σ₀w_kI − R^G_k)) and on the class means as the unregularised K-atom flow's (whitened: −Ã_w, Ã_w = H_L(w) − R̃(w)); the separations are
exactly the directions the plant does not see, so the leak acts exactly on the unexcited subspace. With the Gaussian kernel the classes
attract each other with the force σ₀w_kw_le^{−d²/(2R²)}d, an exponentially small bias.

This module: the kernel and its energy, the averaged loop with the cohesion leak on top of an `AveragedSlowLoop` (its learning term only),
the clusters of a configuration (single linkage), the linear theory at a replicated configuration (the predicted spectrum and the Jacobian
by central differences), and the real closed loop (CL) with the cohesion leak, σ = σ₀/k_e, with the lift of a K-atom configuration.
"""
import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage

from .periodic import ClosedLoop


def psi(z, kind):
    """ψ(z) = (1 − z)³₊ ('compact') or e^{−z/2} ('gauss')."""
    if kind == 'compact':
        return np.clip(1.0 - z, 0.0, None) ** 3
    if kind == 'gauss':
        return np.exp(-0.5 * z)
    raise ValueError(kind)


def omega(z, kind):
    """ω(z) with ω′ = ψ/2, ω(0) = 0: (1 − (1 − z)⁴₊)/8 or 1 − e^{−z/2}."""
    if kind == 'compact':
        return (1.0 - np.clip(1.0 - z, 0.0, None) ** 4) / 8.0
    if kind == 'gauss':
        return 1.0 - np.exp(-0.5 * z)
    raise ValueError(kind)


def pair_d2(Th, G):
    """(differences (m, m, p), squared G⁻¹ distances (m, m))."""
    D = Th[:, None, :] - Th[None, :, :]
    return D, np.sum(D ** 2 / np.asarray(G, float)[None, None, :], axis=2)


def cohesion_force(Th, G, R, sigma0, kind='compact'):
    """−(σ₀/m)Σ_jψ(‖ϑ_i − ϑ_j‖²_{G⁻¹}/R²)(ϑ_i − ϑ_j), (m, p): minus the G-gradient of (σ₀m/2)∫∫W dμ_Θdμ_Θ."""
    D, d2 = pair_d2(np.asarray(Th, float), G)
    return -(sigma0 / Th.shape[0]) * np.einsum('ij,ijk->ik', psi(d2 / R ** 2, kind), D)


def interaction_energy(Th, G, R, kind='compact'):
    """∫∫W dμ_Θdμ_Θ = (1/m²)Σ_{i,j}W(ϑ_i − ϑ_j)."""
    D, d2 = pair_d2(np.asarray(Th, float), G)
    return float(R ** 2 * omega(d2 / R ** 2, kind).sum() / Th.shape[0] ** 2)


def clusters(Th, G, thr):
    """Single-linkage clusters at the G⁻¹ distance thr: labels (m,) in 0..K−1."""
    Th = np.asarray(Th, float)
    if Th.shape[0] == 1:
        return np.zeros(1, int)
    return fcluster(linkage(Th / np.sqrt(np.asarray(G, float))[None, :], method='single'), thr, criterion='distance') - 1


class CohesionLoop:
    """(𝔄_m^coh) on top of an `AveragedSlowLoop` `base` (its residual and grad_weighted give the learning term; its own leak is not used),
    with the range R, the kernel kind and the strength σ₀ (default base.sigma0)."""

    def __init__(self, base, R, kind='compact', sigma0=None):
        self.base = base
        self.G = base.G
        self.R = float(R)
        self.kind = kind
        self.sigma0 = base.sigma0 if sigma0 is None else float(sigma0)

    def learning(self, Th):
        """G⟨∇_ϑ(r_ΘΦ̃(·;ϑ_i))⟩_ρ, (m, p)."""
        return self.G[None, :] * self.base.grad_weighted(self.base.residual(Th), Th)

    def drift(self, Th):
        return self.learning(Th) + cohesion_force(Th, self.G, self.R, self.sigma0, self.kind)

    def J(self, Th):
        """J^coh_m(Θ) = mL̄(μ_Θ) + (σ₀m/2)∫∫W dμ_Θdμ_Θ."""
        m = Th.shape[0]
        return m * self.base.L_bar(Th) + 0.5 * self.sigma0 * m * interaction_energy(Th, self.G, self.R, self.kind)

    def run(self, Th, S_end, dt, rec_every=200, callback=None):
        """Midpoint rule in slow time; callback(s, Θ) at the records. Returns Θ at S_end."""
        Th = np.array(Th, float)
        n = int(round(S_end / dt))
        for k in range(n + 1):
            if callback is not None and k % rec_every == 0:
                callback(k * dt, Th)
            if k == n:
                break
            k1 = self.drift(Th)
            Th = Th + dt * self.drift(Th + 0.5 * dt * k1)
        return Th

    def jacobian(self, Th, h=1e-6):
        """The Jacobian of the drift at Θ by central differences, (mp, mp), in the unconstrained coordinates."""
        x = np.asarray(Th, float).ravel()
        d = x.size
        out = np.zeros((d, d))
        for j in range(d):
            e = np.zeros(d); e[j] = h
            out[:, j] = (self.drift((x + e).reshape(Th.shape)) - self.drift((x - e).reshape(Th.shape))).ravel() / (2 * h)
        return out


def predicted_spectrum(base, atoms, counts, sigma0, h=1e-5):
    """The eigenvalues of the cohesion drift's Jacobian predicted at the replicated configuration (counts n_k copies of atoms[k], the
    classes pairwise at least R apart): on the separations of class k, −σ₀w_k + eig(G^{1/2}R_kG^{1/2}) with multiplicity n_k − 1 each
    (R_k the residual curvature at the atom, w = n/m), and on the class means the eigenvalues of the linearised unregularised weighted
    K-atom flow (the whitened weighted Hessian of L_w with the sign changed), computed by central differences of its drift.
    Returns (separation eigenvalues, collective eigenvalues), each sorted."""
    atoms = np.asarray(atoms, float)
    counts = np.asarray(counts, int)
    w = counts / counts.sum()
    G = base.G
    a, c, s, diff, phi = base.features(atoms)
    r = base.F - phi @ (w * a)
    Rk = base.residual_curvature(r, atoms)
    sq = np.sqrt(G)
    sep = []
    for k in range(len(atoms)):
        ev = np.linalg.eigvalsh(sq[:, None] * Rk[k] * sq[None, :])
        sep += list(np.repeat(-sigma0 * w[k] + ev, counts[k] - 1))

    def coll(A):
        a_, c_, s_, d_, ph = base.features(A)
        return G[None, :] * base.grad_weighted(base.F - ph @ (w * a_), A)
    d = atoms.size
    Jm = np.zeros((d, d))
    for j in range(d):
        e = np.zeros(d); e[j] = h
        Jm[:, j] = (coll((atoms.ravel() + e).reshape(atoms.shape)) - coll((atoms.ravel() - e).reshape(atoms.shape))).ravel() / (2 * h)
    return np.sort(np.array(sep)), np.sort(np.linalg.eigvals(Jm).real)


class CohesionRealLoop(ClosedLoop):
    """(CL) with the cohesion leak at σ = σ₀/k_e: ϑ̇_i = G∇_ϑ(eᵀΦ̃(x;ϑ_i)) − (σ₀/(k_em))Σ_jψ(‖ϑ_i − ϑ_j‖²_{G⁻¹}/R²)(ϑ_i − ϑ_j)."""

    def __init__(self, net, f_vec, x_ref, x_ref_dot, period, ke, G_diag, sigma0, R, kind='compact'):
        super().__init__(net, f_vec, x_ref, x_ref_dot, period, ke, G_diag, sigma0)
        self.R = float(R)
        self.kind = kind

    def rhs(self, t, z, m):
        n, p = self.net.n, self.p
        x, Th = z[:n], z[n:n + m * p].reshape(m, p)
        e = x - self.xr(t)
        xdot = self.f_vec(x[None, :])[0] + self.xrd(t) - self.ke * e - self.net.f_hat(x, Th)
        Thdot = self.G[None, :] * self.net.grad_g(x, e, Th) + cohesion_force(Th, self.G, self.R, self.sigma0 / self.ke, self.kind)
        return np.concatenate([xdot, Thdot.ravel()])

    def lift(self, atoms, counts):
        """The state at t = 0 of perfect tracking at the K-atom configuration: x = x_m(0) and n_k copies of atoms[k]."""
        return np.concatenate([np.atleast_1d(self.xr(0.0)), np.repeat(np.asarray(atoms, float), np.asarray(counts, int), axis=0).ravel()])
