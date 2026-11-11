"""
Proposition D‴.59: the transit with moving widths — a Gaussian feature of any width formed by one kernel
========================================================================================================

The paradigm's kernel law moves the widths beside the centres (σ̇ = Gβa⟨r∂_σφ⟩, ċ = Gβa⟨r∂_cφ⟩: the same mobility, Proposition D‴.48(v)).
On one-dimensional samples a kernel φ(x; c, σ) = exp(−(x − c)²/(2σ²)) is then a point θ = (c, σ) of the plane, and the sliding path,
the loop and the block's polish of `rich_transit_nd` run unchanged on `BumpsCW` — θ in place of a centre in ℝ², ⟨r∇_θφ⟩ the pair
(⟨r∂_cφ⟩, ⟨r∂_σφ⟩) —, the inactive kernels staying at their starting points (c_j, σ₀).

One feature, one kernel, in the continuum: with the coefficient at its least-squares value the fit is the reduced energy
J(c, σ) = ½‖f*‖² − ½⟨f*, φ_{c,σ}⟩²/‖φ_{c,σ}‖² (`reduced_energy`), and ⟨f*, φ_{c,σ}⟩²/‖φ_{c,σ}‖² ∝ (σσ*²/(σ² + σ*²))e^{−d²/(σ² + σ*²)}
(d = c − c*) has the single critical point (c*, σ*) — its maximum —, so the kernel's flow converges to the exact fit from every start
whose coefficient has the feature's sign (on the free law's curve the coefficient cannot change sign), whatever the width mismatch
(`one_kernel_grad` for checks).
"""
import numpy as np


class BumpsCW:
    """Gaussian kernels with their own widths on the one-dimensional samples xs; a kernel is θ = (c, σ) (rows of an array (d, 2)); ⟨f g⟩
    the mean over the samples. The interface of `rich_transit_nd.BumpsND` with n = 2 (sig: the starting width, the scale of a step)."""

    def __init__(self, xs, sig):
        self.x = np.asarray(xs, float).reshape(-1)
        self.xs = self.x[:, None]
        self.N, self.n, self.sig = self.x.size, 2, float(sig)

    def _parts(self, th):
        th = np.atleast_2d(np.asarray(th, float))
        c, s = th[:, 0], th[:, 1]
        u = self.x[:, None] - c[None, :]
        return u, s, np.exp(-u ** 2 / (2 * s ** 2))

    def phi(self, th):
        """φ_j at the samples: (N, d)."""
        return self._parts(th)[2]

    def dphi(self, th):
        """(∂_cφ_j, ∂_σφ_j) at the samples: (N, d, 2)."""
        u, s, P = self._parts(th)
        return np.stack([P * u / s ** 2, P * u ** 2 / s ** 3], -1)

    def synth(self, th, w):
        return self.phi(th) @ w

    def corr(self, th, r, second=False):
        """⟨rφ_j⟩ (d), ⟨r∇_θφ_j⟩ (d, 2) and, with second, ⟨r∇²_θφ_j⟩ (d, 2, 2)."""
        u, s, P = self._parts(th)
        rP = r[:, None] * P
        a, b = u / s ** 2, u ** 2 / s ** 3                             # ∂_c log φ, ∂_σ log φ
        g = rP.sum(0) / self.N
        h = np.stack([(rP * a).sum(0), (rP * b).sum(0)], -1) / self.N
        if not second:
            return g, h
        Q = np.empty((g.size, 2, 2))
        Q[:, 0, 0] = (rP * (a * a - 1 / s ** 2)).sum(0) / self.N
        Q[:, 0, 1] = Q[:, 1, 0] = (rP * (a * b - 2 * u / s ** 3)).sum(0) / self.N
        Q[:, 1, 1] = (rP * (b * b - 3 * u ** 2 / s ** 4)).sum(0) / self.N
        return g, h, Q

    def grams(self, th):
        """P = ⟨φ_iφ_j⟩ (d, d), D[i, j] = ⟨φ_i∇_θφ_j⟩ (d, d, 2), H[i, j] = ⟨∇_θφ_i∇_θφ_jᵀ⟩ (d, d, 2, 2)."""
        P, dP = self.phi(th), self.dphi(th)
        return (P.T @ P / self.N, np.einsum('si,sjk->ijk', P, dP) / self.N, np.einsum('sik,sjl->ijkl', dP, dP) / self.N)


def overlap(c, s, c2, s2):
    """∫φ(x; c, σ)φ(x; c₂, σ₂)dx = (2π)^{1/2}σσ₂/(σ² + σ₂²)^{1/2}·exp(−(c − c₂)²/(2(σ² + σ₂²)))."""
    v = s ** 2 + s2 ** 2
    return np.sqrt(2 * np.pi) * s * s2 / np.sqrt(v) * np.exp(-(c - c2) ** 2 / (2 * v))


def reduced_energy(c, s, A, cs, ss):
    """The one-kernel fit of Aφ(·; c*, σ*) in the continuum with the coefficient at its least-squares value, ½‖f*‖² − ½⟨f*, φ⟩²/‖φ‖²
    (‖φ_σ‖² = π^{1/2}σ); zero exactly at (c*, σ*)."""
    return 0.5 * A ** 2 * (np.sqrt(np.pi) * ss - overlap(c, s, cs, ss) ** 2 / (np.sqrt(np.pi) * s))


def one_kernel_grad(c, s, A, cs, ss):
    """∇_(c, σ) of `reduced_energy` in closed form: with v = σ² + σ*² and ⟨f*, φ⟩²/‖φ‖² = 2π^{1/2}A²σσ*²e^{−d²/v}/v,
    ∂_c = (2d/v)·E₁ and ∂_σ = −(1/σ − 2σ/v + 2σd²/v²)·E₁, E₁ = π^{1/2}A²σσ*²e^{−d²/v}/v."""
    v = s ** 2 + ss ** 2
    d = c - cs
    E1 = np.sqrt(np.pi) * A ** 2 * s * ss ** 2 * np.exp(-d ** 2 / v) / v
    return np.array([2 * d / v * E1, -(1 / s - 2 * s / v + 2 * s * d ** 2 / v ** 2) * E1])
