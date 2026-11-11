"""
The averaged slow loop for an n-dimensional plant (Proposition D‴.24's planar check)
=====================================================================================

`chaos_learning.AveragedSlowLoop` and the classes built on it are written for n = 1. This module gives the averaged slow loop of Theorem J
for the saturated network on ℝⁿ (n ≥ 1), with a particle ϑ = (α, γ, κ) ∈ ℝⁿ × ℝⁿ × ℝ (p = 2n + 1) and the feature
Φ̃(x; ϑ) = a φ(x; c, s) ∈ ℝⁿ, a = a_max tanh(α/a_max), c = c_max tanh(γ/c_max), s = s_min + (s_max − s_min)sig(κ) (`langevin`):

    dϑ_i/ds = G⟨r_Θ(x)ᵀ∇_ϑΦ̃(x; ϑ_i)⟩_ρ − σ₀(ϑ_i − T_i),    r_Θ = f − f̂_Θ,    f̂_Θ = Σ_iw_iΦ̃(·; ϑ_i)  (w_i = 1/m by default),

ρ the uniform measure of one period of the reference orbit. It is the average over a period of the real loop's law
ϑ̇ = G∇_ϑ(eᵀΦ̃(x; ϑ)) − σ(ϑ − T) at the slaved error e = r(x_m)/k_e, in the slow time s = t/k_e with σ₀ = k_eσ (Theorem J), and the
gradient flow of L̄ = ½⟨|r|²⟩_ρ plus the leak in the metric (1/m)G⁻¹ — the same objects as in one dimension, for every n. With weights w
on K atoms it is the replicated (weighted K-atom) system of Proposition D‴.5.

On top of it: the nearest-atom leak toward atoms Z (the cell plan in the G⁻¹ metric), the adaptive prior of Proposition D‴.24(ii) (the atoms
move toward the means of their cells at the rate η; η = ∞: the atoms are the cell means), a midpoint integrator, the weighted stationary
point of the fixed prior by Newton, and the whitened weighted Jacobian Ĵ and Gauss–Newton block H_L = ĴᵀĴ of Theorem D‴.23(ii).
"""
import numpy as np
from scipy.optimize import root
from scipy.special import expit as _sig


class AveragedLoopN:
    """The averaged slow loop on ℝⁿ: plant f (f_vec: (N, n) → (N, n)), reference x_ref(t) → (n,), period, gains G (p,), leak σ₀; ρ the
    uniform grid of n_quad times of one period."""

    def __init__(self, net, f_vec, x_ref, period, G, sigma0, n_quad=400):
        self.net = net
        self.n = int(net.n)
        self.p = 2 * self.n + 1
        self.G = np.asarray(G, float)
        assert self.G.shape == (self.p,)
        self.sigma0 = float(sigma0)
        self.period = float(period)
        self.N = int(n_quad)
        ts = np.arange(self.N) * self.period / self.N
        self.X = np.array([np.atleast_1d(x_ref(t)) for t in ts]).reshape(self.N, self.n)
        self.F = np.asarray(f_vec(self.X), float).reshape(self.N, self.n)
        self.F2 = float(np.mean(np.sum(self.F ** 2, axis=1)))

    # ------------------------------------------------------------------ features
    def features(self, Th):
        """(a (m, n), c (m, n), s (m,), diff (N, m, n), sq (N, m), φ (N, m))."""
        a, c, s = self.net.physical(np.asarray(Th, float))
        diff = self.X[:, None, :] - c[None, :, :]
        sq = np.sum(diff ** 2, axis=2)
        phi = np.exp(-sq / (2 * s[None, :] ** 2))
        return a, c, s, diff, sq, phi

    def weights(self, Th, w):
        return np.full(Th.shape[0], 1.0 / Th.shape[0]) if w is None else np.asarray(w, float)

    def f_hat(self, Th, w=None):
        a, c, s, diff, sq, phi = self.features(Th)
        return phi @ (self.weights(Th, w)[:, None] * a)

    def residual(self, Th, w=None):
        return self.F - self.f_hat(Th, w)

    def L_bar(self, Th, w=None):
        """½⟨|r|²⟩_ρ."""
        return 0.5 * float(np.mean(np.sum(self.residual(Th, w) ** 2, axis=1)))

    def fit_error(self, Th, w=None):
        """(⟨|r|²⟩_ρ/⟨|f|²⟩_ρ)^{1/2}."""
        return float(np.sqrt(np.mean(np.sum(self.residual(Th, w) ** 2, axis=1)) / self.F2))

    def learning(self, Th, w=None, r=None):
        """⟨r(x)ᵀ∇_ϑΦ̃(x; ϑ_i)⟩_ρ in the unconstrained coordinates (the chain rule through the saturations), (m, p); r the residual of Θ with
        the weights w unless given."""
        Th = np.asarray(Th, float)
        n, net = self.n, self.net
        a, c, s, diff, sq, phi = self.features(Th)
        if r is None:
            r = self.F - phi @ (self.weights(Th, w)[:, None] * a)
        ra = r @ a.T                                                        # (N, m): r(x)·a_i
        da = 1.0 - (a / net.a_max) ** 2
        dc = 1.0 - (c / net.c_max) ** 2
        sk = _sig(Th[:, 2 * n]); ds = (net.s_max - net.s_min) * sk * (1 - sk)
        out = np.empty_like(Th)
        out[:, :n] = (phi.T @ r) * da
        out[:, n:2 * n] = np.einsum('xi,xic->ic', ra * phi, diff) / s[:, None] ** 2 * dc
        out[:, 2 * n] = np.sum(ra * phi * sq, axis=0) / s ** 3 * ds
        return out / self.N

    def drift(self, Th, T, w=None):
        """dΘ/ds = G⟨rᵀ∇_ϑΦ̃(ϑ_i)⟩_ρ − σ₀(ϑ_i − T_i) for the leak targets T (m, p)."""
        return self.G[None, :] * self.learning(Th, w) - self.sigma0 * (np.asarray(Th, float) - np.asarray(T, float))


def cost_matrix(Th, Z, G):
    """‖ϑ_i − z_j‖²_{G⁻¹}, (m, K)."""
    return np.sum((np.asarray(Th, float)[:, None, :] - np.asarray(Z, float)[None, :, :]) ** 2 / np.asarray(G, float)[None, None, :], axis=2)


def cells(Th, Z, G):
    """The nearest atom of every particle (G⁻¹ metric), (m,)."""
    return np.argmin(cost_matrix(Th, Z, G), axis=1)


def cell_means(Th, lab, Z):
    """The mean of every cell; an empty cell keeps its atom."""
    out = np.array(Z, float).copy()
    for j in range(out.shape[0]):
        if np.any(lab == j):
            out[j] = Th[lab == j].mean(axis=0)
    return out


def run(loop, Th, Z, S_end, dt, eta=0.0, cell_every=4, rec_every=10 ** 9, callback=None):
    """The nearest-atom leak toward the atoms Z, fixed (η = 0) or adaptive (dZ/ds = −η(Z − cell means); η = ∞: Z = the cell means, set at
    every refresh of the cells), midpoint rule in slow time; the cells are recomputed every `cell_every` steps and frozen in between.
    `callback(s, Θ, Z)` at the records. Returns (Θ, Z)."""
    Th = np.array(Th, float); Z = np.array(Z, float)
    n_steps = int(round(S_end / dt))
    lab = None
    for k in range(n_steps + 1):
        if k % cell_every == 0:
            lab = cells(Th, Z, loop.G)
            if np.isinf(eta):
                Z = cell_means(Th, lab, Z)
        if callback is not None and k % rec_every == 0:
            callback(k * dt, Th, Z)
        if k == n_steps:
            break
        T = Z[lab]
        k1 = loop.drift(Th, T)
        Thm = Th + 0.5 * dt * k1
        if np.isfinite(eta) and eta > 0:
            Zm = Z - 0.5 * dt * eta * (Z - cell_means(Th, lab, Z))
            Tm = Zm[lab]
            Z = Z - dt * eta * (Zm - cell_means(Thm, lab, Zm))
        else:
            Tm = T
        Th = Th + dt * loop.drift(Thm, Tm)
    return Th, Z


def weighted_stationary(loop, w, prior, init=None, tol=1e-13):
    """The weighted K-atom stationary point of the fixed prior with every atom leaking to its own prior atom: G⟨r_w∇Φ̃(ϑ_k)⟩ = σ₀(ϑ_k − ϑ⁰_k),
    r_w = f − Σ_lw_lΦ̃(ϑ_l) (scipy's hybrid Powell method from `init`, default the prior). Returns (atoms, residual, whether every atom is in
    its own prior cell)."""
    K, p = prior.shape
    x0 = (prior if init is None else np.asarray(init, float)).ravel()
    fun = lambda v: loop.drift(v.reshape(K, p), prior, w).ravel()
    sol = root(fun, x0, method='hybr', tol=tol)
    atoms = sol.x.reshape(K, p)
    return atoms, float(np.abs(fun(sol.x)).max()), bool(np.array_equal(cells(atoms, prior, loop.G), np.arange(K)))


def whitened_jacobian(loop, w, atoms, h=1e-6):
    """Ĵ: the derivative of f̂ = Σ_kw_kΦ̃(·; ϑ_k) on the grid in the whitened weighted coordinates y_k = √w_kG^{−1/2}(ϑ_k − ϑ*_k), as an
    ((N·n) × (K·p)) matrix (central differences); H_L = ĴᵀĴ/N."""
    K, p = atoms.shape
    w = np.asarray(w, float)
    J = np.zeros((loop.N * loop.n, K * p))
    for k in range(K):
        for c in range(p):
            e = np.zeros_like(atoms); e[k, c] = h
            J[:, k * p + c] = ((loop.f_hat(atoms + e, w) - loop.f_hat(atoms - e, w)) / (2 * h)).ravel()
    return J * (np.tile(np.sqrt(loop.G), K) / np.repeat(np.sqrt(w), p))[None, :]


def gauss_newton(loop, w, atoms):
    """(H_L, eigenvalues ascending, eigenvectors) of the whitened weighted Gauss–Newton block at the K-atom measure Σ_kw_kδ_{ϑ_k}."""
    J = whitened_jacobian(loop, w, atoms)
    H = J.T @ J / loop.N
    ev, V = np.linalg.eigh(H)
    return H, ev, V


def whitened(D, w, G):
    """√w_kG^{−1/2}D_k, flattened."""
    return (np.sqrt(np.asarray(w, float))[:, None] * np.asarray(D, float) / np.sqrt(np.asarray(G, float))[None, :]).ravel()
