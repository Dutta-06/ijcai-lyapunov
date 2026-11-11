"""
Proposition D‴.37 (PARADIGM.md §5d): the dither's drift along a manifold of equilibria
=====================================================================================

Let ẋ = F(x) have a manifold M of equilibria, parametrised by a chart w ↦ γ(w), that attracts its neighbourhood, with the linearisation
A = ∂F diagonalisable at M: the eigenvalue 0 on the tangent T = range ∂γ, the others with negative real parts. The limit map Φ (the point of
M the flow converges to) has ∂Φ = P, the spectral projection onto T along range A, and, from differentiating ∂Φ(x)F(x) ≡ 0 twice and
Φ∘γ = γ, in an eigenbasis (e_i, λ_i) of A completed on T by the chart's columns,

    ∂²Φ[e_i, e_j] = −P ∂²F[e_i, e_j] / (λ_i + λ_j)        unless both e_i, e_j ∈ T,
    ∂²Φ[t, t′]    = (I − P) ∂²γ[c, c′]                     for t = ∂γ c, t′ = ∂γ c′ ∈ T.

Under the dither dx = F ds + √Σ dB with Σ = O(τ), Katzenberger's limit on the time scale 1/τ is dY = P√Σ dB + ½∂²Φ[Σ] ds on M, and an
observable w(x) (the chart's coordinates, read off by `coords`) moves by dw = ∂w·dY + ½∂²w[PΣPᵀ] ds: its drift is
∂w ½∂²Φ[Σ] + ½∂²w[PΣPᵀ], its diffusion ∂w PΣPᵀ ∂wᵀ. A deterministic perturbation of the same order, dx = (F + F₁) ds + √Σ dB with
F₁ = O(τ) (Proposition D‴.38's correction), adds ∂w P F₁ to the drift: on the slow time scale ∂Φ F₁ = P F₁ is its part along M.

This module: that drift and diffusion from F, the chart and Σ, by finite differences (`limit_drift`), and the limit map itself by integrating
the flow (`limit_map`, for checks).
"""
import numpy as np


def _jac(f, x, h):
    return np.array([(f(x + h * e) - f(x - h * e)) / (2 * h) for e in np.eye(len(x))]).T


def _d2(f, x, u, v, h):
    """∂²f(x)[u, v] by polarisation."""
    return (f(x + h * (u + v)) - f(x + h * (u - v)) - f(x - h * (u - v)) + f(x - h * (u + v))) / (4 * h * h)


def spectral_projection(A, Jg):
    """P, the projection onto T = range Jg along range A, without A's eigenvectors: P = Jg(WᵀJg)⁻¹Wᵀ, W spanning ker Aᵀ = (range A)^⊥ (the
    left singular vectors of A for its dim T smallest singular values); and L = (WᵀJg)⁻¹Wᵀ, the chart coordinates of P (P = Jg L)."""
    k = Jg.shape[1]
    W = np.linalg.svd(A)[0][:, -k:]
    L = np.linalg.solve(W.T @ Jg, W.T)
    return Jg @ L, L


def limit_drift(F, chart, coords, w0, Sigma, h=1e-4, h2=2e-3, by_noise=False, F1=None, method='eig', cond_max=1e8):
    """The drift and the diffusion of the coordinates w = coords(x) of Katzenberger's limit at the point chart(w0) of M, for the noise
    covariance Sigma (per unit time) and an optional perturbation F1 of the drift of the same order (its part ∂w P F1(x0) is added).
    Returns (drift (dim w), diffusion (dim w × dim w), the nonzero eigenvalues of A); with `by_noise` the drift is also split by the noise
    coordinates (Sigma diagonal), as a (n, dim w) array (F1's part not included in the split).
    method 'eig': A's eigenbasis, as in the module's formulas — used unless the computed eigenvector matrix has a condition number above
    cond_max (a degenerate eigenvalue, where it is ill-determined); method 'schur' (then): the same formulas without eigenvectors —
    ½∂²Φ[S] = ½(I − P)∂²γ[LSLᵀ] − ½P∂²F[X], X the symmetric solution of AX + XAᵀ = S − PSPᵀ (X = X_N + PSQᵀA_Dᵀ + A_DQSPᵀ, Q = I − P,
    A_D = (A − σP)⁻¹Q, X_N from the Lyapunov equation of A − σP with QSQᵀ; Bartels–Stewart), ∂²F[X] = Σ_k μ_k∂²F[q_k, q_k] over X's
    eigenpairs; the eigenvalues then lose their imaginary parts of the size of the differences' error."""
    w0 = np.asarray(w0, float)
    x0 = chart(w0)
    n, k = len(x0), len(w0)
    A = _jac(F, x0, h)
    Jg = _jac(chart, w0, h)                                                   # n × k: the tangent
    lam, V = np.linalg.eig(A)
    o = np.argsort(np.abs(lam))
    lam, V = lam[o], V[:, o]
    if np.abs(lam[:k]).max() > 1e-6 * np.abs(lam).max() or np.abs(A @ Jg).max() > 1e-6 * np.abs(A).max() * np.abs(Jg).max():
        raise ValueError('the chart is not a manifold of equilibria of F here')
    if method == 'eig':
        Vc = np.c_[Jg, V[:, k:].real]
        if np.abs(lam[k:].imag).max() > 1e-12 * np.abs(lam).max() or np.linalg.cond(Vc) > cond_max:
            method = 'schur'
    if method == 'schur':
        if np.abs(lam[k:].imag).max() > 1e-6 * np.abs(lam).max():
            raise ValueError('complex eigenvalues off M (not handled)')
        return _limit_drift_schur(F, chart, coords, w0, x0, A, Jg, np.sort(lam[k:].real)[::-1], Sigma, h, h2, by_noise, F1)
    if np.abs(lam[k:].imag).max() > 1e-12 * np.abs(lam).max():
        raise ValueError('complex eigenvalues off M (not handled)')
    lam = np.r_[np.zeros(k), lam[k:].real]
    V = np.c_[Jg, V[:, k:].real]
    Vi = np.linalg.inv(V)
    P = V[:, :k] @ Vi[:k, :]
    Jw = _jac(coords, x0, h)

    def drift_for(S):
        C = Vi @ S @ Vi.T                                                    # the noise in the eigenbasis
        d2P = np.zeros(n)
        for i in range(n):
            for j in range(i, n):
                c = C[i, j] if i == j else C[i, j] + C[j, i]
                if c == 0.0:
                    continue
                if i < k and j < k:
                    ci, cj = np.eye(k)[i], np.eye(k)[j]                       # e_i = ∂γ c_i on T
                    g2 = _d2(chart, w0, ci, cj, h2)
                    B = g2 - P @ g2
                else:
                    B = -P @ _d2(F, x0, V[:, i], V[:, j], h2) / (lam[i] + lam[j])
                d2P += 0.5 * c * B
        Dy = P @ S @ P.T
        ito = np.zeros(k)
        for i in range(n):
            for j in range(n):
                if Dy[i, j] != 0.0:
                    ito += 0.5 * Dy[i, j] * _d2(coords, x0, np.eye(n)[i], np.eye(n)[j], h2)
        return Jw @ d2P + ito

    drift = drift_for(Sigma)
    if F1 is not None:
        drift = drift + Jw @ P @ F1(x0)
    diff = Jw @ P @ Sigma @ P.T @ Jw.T
    if by_noise:
        parts = np.array([drift_for(np.diag(np.eye(n)[i] * Sigma[i, i])) for i in range(n)])
        return drift, diff, lam[k:], parts
    return drift, diff, lam[k:]


def _limit_drift_schur(F, chart, coords, w0, x0, A, Jg, lam, Sigma, h, h2, by_noise, F1):
    """`limit_drift`'s method 'schur' (its docstring)."""
    from scipy.linalg import solve_continuous_lyapunov
    n, k = len(x0), len(w0)
    P, L = spectral_projection(A, Jg)
    Q = np.eye(n) - P
    As = A - np.abs(lam).max() * P                                            # nonsingular, the same as A on range Q
    AD = np.linalg.solve(As, Q)
    Jw = _jac(coords, x0, h)

    def drift_for(S):
        C = L @ S @ L.T
        d2P = np.zeros(n)
        for i in range(k):
            for j in range(i, k):
                c = C[i, j] if i == j else C[i, j] + C[j, i]
                if c != 0.0:
                    g2 = _d2(chart, w0, np.eye(k)[i], np.eye(k)[j], h2)
                    d2P += 0.5 * c * (g2 - P @ g2)
        X = solve_continuous_lyapunov(As, Q @ S @ Q.T) + P @ S @ Q.T @ AD.T + AD @ Q @ S @ P.T
        mu, qv = np.linalg.eigh(0.5 * (X + X.T))
        d2P -= 0.5 * P @ sum(mu[i] * _d2(F, x0, qv[:, i], qv[:, i], h2) for i in range(n))
        Dy = P @ S @ P.T
        ito = np.zeros(k)
        for i in range(n):
            for j in range(n):
                if Dy[i, j] != 0.0:
                    ito += 0.5 * Dy[i, j] * _d2(coords, x0, np.eye(n)[i], np.eye(n)[j], h2)
        return Jw @ d2P + ito

    drift = drift_for(Sigma)
    if F1 is not None:
        drift = drift + Jw @ P @ F1(x0)
    diff = Jw @ P @ Sigma @ P.T @ Jw.T
    if by_noise:
        parts = np.array([drift_for(np.diag(np.eye(n)[i] * Sigma[i, i])) for i in range(n)])
        return drift, diff, lam, parts
    return drift, diff, lam


def limit_map(F, x, S_end, rtol=1e-12, atol=1e-14):
    """Φ(x) by integrating ẋ = F(x) to S_end (LSODA); returns (the end point, the drift's size there)."""
    from scipy.integrate import solve_ivp
    sol = solve_ivp(lambda t, y: F(y), (0.0, S_end), np.asarray(x, float), method='LSODA', rtol=rtol, atol=atol)
    y = sol.y[:, -1]
    return y, float(np.abs(F(y)).max())
