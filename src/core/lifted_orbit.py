"""
Corollary D‴.14 and Proposition D‴.15 (PARADIGM.md §5d): the lifted periodic orbit of the real closed loop and its transverse spectrum
===================================================================================================================================

On a configuration collapsed by class — n_k = mw_k particles at ϑ_k — the real closed loop (CL) of Theorem J with a leak whose target is
the same for every particle of a class (the nearest-atom leak inside the cells; the OT leak with the diagonal plan at masses w = q) is the
weighted K₀-atom loop: the collapsed subspace is invariant, so the K₀-atom loop's T₀-periodic orbit lifts to a T₀-periodic orbit of the
m-loop for every commensurate m, with the same orbit measure, learned function and tracking error (Corollary D‴.14).  The linearisation
of (CL) along a lifted orbit leaves each class's separation subspace S_k = {δ : δ_i = 0 off the class, Σ_{i∈k}δ_i = 0} invariant and acts
on it as M_k ⊗ I_{n_k−1}, M_k the monodromy of one particle of the class driven by the orbit's plant state with weight zero in f̂ (a *ghost*
particle): the Floquet spectrum of the lift is the K₀-atom orbit's together with the ghosts' exponents, each of multiplicity n_k − 1, all
independent of m (Proposition D‴.15).  This module implements the per-class-target loop, the lift, the class-wise separation bases, the
compressions of the monodromy, and the ghost monodromy.
"""
import numpy as np
from scipy.integrate import solve_ivp

from .periodic import ClosedLoop


class ClassLeakLoop(ClosedLoop):
    """(CL) with a per-particle leak target (m, p) in place of the single prior ϑ⁰: the nearest-atom leak with its assignment frozen, or
    the OT leak with its plan frozen — exact where no particle changes cell and the plan does not switch."""

    def __init__(self, net, f_vec, x_ref, x_ref_dot, period, ke, G_diag, sigma0, targets, weights=None):
        super().__init__(net, f_vec, x_ref, x_ref_dot, period, ke, G_diag, sigma0)
        self.targets = np.atleast_2d(np.asarray(targets, float))
        self.weights = None if weights is None else np.asarray(weights, float)      # f̂ = Σ_i ω_iΦ̃(·;ϑ_i); None: uniform 1/m

    def f_hat(self, x, Th):
        if self.weights is None:
            return self.net.f_hat(x, Th)
        a, c, s = self.net.physical(Th)
        phi = np.exp(-np.sum((x[None, :] - c) ** 2, axis=1) / (2 * s ** 2))
        return (self.weights[:, None] * phi[:, None] * a).sum(axis=0)

    def rhs(self, t, z, m):
        n = self.net.n
        x, Th = z[:n], z[n:].reshape(m, self.p)
        xr, xrd = self.xr(t), self.xrd(t)
        e = x - xr
        fx = self.f_vec(x[None, :])[0]
        xdot = fx + xrd - self.ke * e - self.f_hat(x, Th)
        Thdot = self.G[None, :] * self.net.grad_g(x, e, Th) - (self.sigma0 / self.ke) * (Th - self.targets)
        return np.concatenate([xdot, Thdot.ravel()])


def lift(z_K, labels, n, p):
    """The lift of the K₀-atom state z_K = (x, ϑ_1, …, ϑ_{K₀}) to the collapsed m-state (x, ϑ_{labels[1]}, …, ϑ_{labels[m]})."""
    z_K = np.asarray(z_K, float)
    ThK = z_K[n:].reshape(-1, p)
    return np.concatenate([z_K[:n], ThK[np.asarray(labels)].ravel()])


def class_separation_bases(labels, n, p):
    """Orthonormal bases (columns) of the class separation subspaces S_k ⊂ ℝ^{n+mp} — zero plant component, zero mean within the class,
    zero off the class — as a dict k → B_k of shape (n + mp, (n_k − 1)p); classes with one particle have no separation subspace."""
    labels = np.asarray(labels)
    m = len(labels)
    out = {}
    for k in np.unique(labels):
        idx = np.flatnonzero(labels == k)
        nk = len(idx)
        if nk < 2:
            continue
        Q, _ = np.linalg.qr(np.column_stack([np.ones(nk), np.eye(nk)[:, :nk - 1]]))
        V = Q[:, 1:]                                                       # orthonormal basis of {v : Σv = 0} in ℝ^{n_k}
        B = np.zeros((n + m * p, (nk - 1) * p))
        for a in range(nk - 1):
            for j in range(p):
                col = np.zeros(n + m * p)
                col[n + idx * p + j] = V[:, a]
                B[:, a * p + j] = col
        out[int(k)] = B
    return out


def collective_basis(labels, n, p):
    """Orthonormal basis of the complement of ⊕_k S_k: the plant coordinates and the class means, (n + mp, n + K₀p)."""
    labels = np.asarray(labels)
    m = len(labels)
    cols = [np.eye(n + m * p)[:, i] for i in range(n)]
    for k in np.unique(labels):
        idx = np.flatnonzero(labels == k)
        for j in range(p):
            col = np.zeros(n + m * p)
            col[n + idx * p + j] = 1.0 / np.sqrt(len(idx))
            cols.append(col)
    return np.column_stack(cols)


def compression_spectrum(DP, B, period):
    """Eigenvalues of Bᵀ DP B (the monodromy compressed to the span of B), their Floquet exponents −log|λ|/T₀ sorted ascending, and the
    invariance defect ‖DP B − B(Bᵀ DP B)‖_F/‖DP B‖_F (zero iff span B is invariant)."""
    C = B.T @ DP @ B
    lam = np.linalg.eigvals(C)
    defect = float(np.linalg.norm(DP @ B - B @ C) / max(np.linalg.norm(DP @ B), 1e-300))
    return lam, np.sort(-np.log(np.abs(lam)) / period), defect


def ghost_monodromy(cl_K, z_p, K0, k, target_k, fd=1e-7, rtol=1e-11, atol=1e-13):
    """M_k: the one-period map's Jacobian of a ghost particle of class k — it feels the plant state (x, e) of the K₀-atom loop's periodic
    orbit through z_p and leaks toward target_k, but has weight zero in f̂ — by central differences from the class atom; returns
    (M_k (p × p), its Floquet exponents sorted ascending, the ghost's return residual)."""
    n, p = cl_K.net.n, cl_K.p
    tgt = np.asarray(target_k, float)

    def rhs(t, z):
        x = z[:n]; ThK = z[n:n + K0 * p].reshape(K0, p); g = z[n + K0 * p:].reshape(1, p)
        xr, xrd = cl_K.xr(t), cl_K.xrd(t)
        e = x - xr
        fx = cl_K.f_vec(x[None, :])[0]
        xdot = fx + xrd - cl_K.ke * e - cl_K.f_hat(x, ThK)
        ThKdot = cl_K.G[None, :] * cl_K.net.grad_g(x, e, ThK) - (cl_K.sigma0 / cl_K.ke) * (ThK - cl_K.targets)
        gdot = cl_K.G[None, :] * cl_K.net.grad_g(x, e, g) - (cl_K.sigma0 / cl_K.ke) * (g - tgt[None, :])
        return np.concatenate([xdot, ThKdot.ravel(), gdot.ravel()])

    def period_map(g0):
        sol = solve_ivp(rhs, (0.0, cl_K.T0), np.concatenate([z_p, g0]), method='DOP853', rtol=rtol, atol=atol)
        if not sol.success:
            raise RuntimeError(sol.message)
        return sol.y[-p:, -1]

    g_star = z_p[n + k * p:n + (k + 1) * p]
    residual = float(np.linalg.norm(period_map(g_star) - g_star))
    M = np.zeros((p, p))
    for j in range(p):
        ej = np.zeros(p); ej[j] = fd
        M[:, j] = (period_map(g_star + ej) - period_map(g_star - ej)) / (2 * fd)
    lam = np.linalg.eigvals(M)
    return M, np.sort(-np.log(np.abs(lam)) / cl_K.T0), residual


def monodromy_fd(cl, z_p, m, fd=1e-7):
    """The Jacobian of the Poincaré map at z_p by central differences (the monodromy matrix when z_p is a fixed point)."""
    d = len(z_p)
    DP = np.zeros((d, d))
    for j in range(d):
        ej = np.zeros(d); ej[j] = fd * max(1.0, abs(z_p[j]))
        DP[:, j] = (cl.poincare(z_p + ej, m) - cl.poincare(z_p - ej, m)) / (2 * ej[j])
    return DP
