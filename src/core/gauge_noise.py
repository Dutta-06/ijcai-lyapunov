"""
Proposition D‴.40 (PARADIGM.md §5d): the gauge under general noise
================================================================

At an exact fit with the masses held fixed, white noise of covariance Σ (per unit slow time) on the positions makes them fluctuate about
the fit: to leading order a Gaussian of covariance C, the solution of the Lyapunov equation

    (−GH)C + C(−GH)ᵀ + Σ = 0,   H = ∇²_ΘJ at the fit, G the positions' metric (dϑ = −G∇J ds + noise),

displaced by the mean shift μ = −½H⁻¹∇_Θ tr(H(Θ)C)|_C (the cubic term of J). The masses' fitness u_k = −∂_{β_k}J (zero at the fit)
averaged over these fluctuations is

    ū_k = ∇u_k·μ + ½ tr(∇²u_k C),

and in the slow-replicator limit the masses follow the replicator flow of ū (dW_k/ds = γW_k(ū_k − Σ_lW_lū_l)) — for Σ = 2τ diag(G) the
descent of the entropic free energy (τ/2) log det H (Proposition D‴.37(ii)); for any Σ, at this order, ū_k = −(τ_k/m)(1 + (a_k/a_max)²)/
(W_ka′_k) with τ_k = Σ_{α_kα_k}/(2G_α) only (Proposition D‴.40: the gauge identity makes ū_k the generator's image of a function of α_k),
the Shahshahani gradient flow descending Σ_kτ_k log(W_k a′_k); beyond this order the averages depend on all of Σ.

This module: the averaged fitness by central differences (`averaged_fitness`), from callables for ∇_ΘJ and the fitness vector; for an
Ornstein–Uhlenbeck noise η of correlation time τ_c with the same integrated covariance Σ (covariance Σe^{−|s−s′|/τ_c}/(2τ_c)), the
linearised positions' cross-covariance with the noise, E[δϑ ηᵀ] = (I + τ_cGH)⁻¹Σ/2 (`coloured_cross`; the white noise's Σ/2 at τ_c = 0),
which sets the coloured noise's level −E[g_kη_{α_k}]/G_α = −g′_k[(I + τ_cGH)⁻¹Σ]_{α_kα_k}/(2G_α) at leading order (referee round 147), and
the same from the Lyapunov equation of the system augmented with η (`coloured_cross_augmented`, the check).
"""
import numpy as np
from scipy.linalg import solve_continuous_lyapunov


def averaged_fitness(gradJ, u, theta, Sigma, Gvec, h=1e-4, hH=1e-5):
    """ū (the fitness averaged over the positions' stationary fluctuations at fixed masses), C and μ at the exact fit `theta` (n,): `gradJ(t)`
    the gradient of J in the positions (n,), `u(t)` the fitness vector (K,), `Sigma` the positions' noise covariance (n × n), `Gvec` the
    metric's diagonal (n,). Central differences: steps `hH` for the Hessian of J, `h` for the third derivatives and those of u."""
    theta = np.asarray(theta, float)
    n = len(theta)
    I = np.eye(n)
    H = np.array([(gradJ(theta + hH * e) - gradJ(theta - hH * e)) / (2 * hH) for e in I]).T
    H = 0.5 * (H + H.T)
    A = -(np.asarray(Gvec, float)[:, None] * H)
    C = solve_continuous_lyapunov(A, -np.asarray(Sigma, float))
    C = 0.5 * (C + C.T)
    lam, V = np.linalg.eigh(C)

    def trHC(t):                                                              # tr(H(t)C) = Σ_m λ_m ∂²_{v_m}J(t)
        return sum(l * v @ (gradJ(t + h * v) - gradJ(t - h * v)) / (2 * h) for l, v in zip(lam, V.T) if l != 0.0)
    g = np.array([(trHC(theta + h * e) - trHC(theta - h * e)) / (2 * h) for e in I])
    mu = -0.5 * np.linalg.solve(H, g)
    gu = np.array([(u(theta + h * e) - u(theta - h * e)) / (2 * h) for e in I]).T          # K × n
    K = gu.shape[0]
    Hu = np.zeros((K, n, n))
    for i in range(n):
        for j in range(i, n):
            ei, ej = I[i] * h, I[j] * h
            d = (u(theta + ei + ej) - u(theta + ei - ej) - u(theta - ei + ej) + u(theta - ei - ej)) / (4 * h * h)
            Hu[:, i, j] = d
            Hu[:, j, i] = d
    return gu @ mu + 0.5 * np.einsum('kij,ij->k', Hu, C), C, mu


def replicator_flow(W, ubar):
    """dW_k/ds per unit γ for the fitness ū: W_k(ū_k − Σ_lW_lū_l)."""
    W = np.asarray(W, float)
    return W * (ubar - W @ ubar)


def coloured_cross(GH, Sigma, tc):
    """E[δϑ ηᵀ] in the stationary law of dδϑ = −GH δϑ ds + η ds, η an Ornstein–Uhlenbeck noise of correlation time `tc` whose integrated
    covariance is `Sigma` (η's covariance Σe^{−|s−s′|/tc}/(2tc)): ∫₀^∞ e^{−GHt}e^{−t/tc}dt Σ/(2tc) = (I + tc·GH)⁻¹Σ/2."""
    GH = np.asarray(GH, float)
    return 0.5 * np.linalg.solve(np.eye(len(GH)) + tc * GH, np.asarray(Sigma, float))


def coloured_cross_augmented(GH, Sigma, tc):
    """The same from the stationary covariance of the augmented linear system (δϑ, η): drift [[−GH, I], [0, −I/tc]], noise on η of
    covariance Σ/tc² per unit time; the (δϑ, η) block of the Lyapunov solution."""
    GH = np.asarray(GH, float); n = len(GH)
    M = np.block([[-GH, np.eye(n)], [np.zeros((n, n)), -np.eye(n) / tc]])
    Qn = np.zeros((2 * n, 2 * n)); Qn[n:, n:] = np.asarray(Sigma, float) / tc ** 2
    C = solve_continuous_lyapunov(M, -Qn)
    return C[:n, n:]
