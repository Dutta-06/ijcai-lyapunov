"""
Perfect tracking with an adaptive prior in the real closed loop: the Floquet spectrum (Corollary D‴.26, PARADIGM.md §5d)
=========================================================================================================================

The real closed loop (CL) of Theorem J with the nearest-atom leak toward an adaptive prior: the particles leak at the rate σ = σ₀/k_e toward
the atom z_{π(i)} of their class, and the atoms move toward the class means at the rate η/k_e,

    ẋ = f(x) + ẋ_m − k_e(x − x_m) − f̂_Θ(x),   ϑ̇_i = G∇_ϑ(eᵀΦ̃(x; ϑ_i)) − (σ₀/k_e)(ϑ_i − z_{π(i)}),   ż_j = −(η/k_e)(z_j − ϑ̄_j),

ϑ̄_j the mean of class j (η = ∞: z_j = ϑ̄_j, eliminated). The classes are frozen: exact near a replicated configuration whose particles lie
inside their cells. With Θ the lift of an exact representation ϑ^e(w) of the sampled weights and Z = ϑ^e(w), the plant on the reference is
a T₀-periodic solution with e ≡ 0 (Proposition D‴.24(vi)); this module gives the loop with the prior in the state, the lift, and the
predicted slow Floquet exponents (1/k_e)·{σ₀ on the separations, λ_j^±(w) on the collective pairs} of Theorem D‴.23(ii) at the weights w;
`monodromy` (central differences of the period map) and `periodic.floquet_exponents` compute the actual ones.

`CollapsedVariational` is the variational equation along the lift restricted to the collapsed subspace (the weighted K₀-atom loop's
linearisation) in the whitened weighted coordinates, where it is the adjoint pair δė = −k_eδe − Ĵ(t)y, ẏ = Ĵ(t)ᵀδe − (σ₀/k_e)(y − ζ),
ζ̇ = −(η/k_e)(ζ − y): its monodromy by integrating the variational equation (no differences), the error channel's exponent by Liouville's
formula, and the averaged generator of the slow block to second order in 1/k_e (Corollary D‴.26(c)): after the slaving of δe,
dY/ds = (B₀ + εB₁ + ε²B₂)Y in slow time (ε = 1/k_e), B₀ = C − J̄ᵀJ̄, B₁ = J̄ᵀJ̄̇, B₂ = −J̄ᵀJ̄̈ + J̄ᵀJ̄C − (J̄ᵀJ̄)² (J̄ = Ĵ on the y block, C the
leak and the prior), and the near-identity averaging gives Ā₀ = ⟨B₀⟩, Ā₁ = ⟨B₁⟩ + ⟨B̃₀U₁⟩ (skew in the metric diag(I, (σ₀/η)I)), Ā₂ =
⟨B₂⟩ + ⟨B₁U₁⟩ + ⟨B₀U₂⟩ (U₁′ = B̃₀, U₂′ = B₁ + B₀U₁ − U₁Ā₀ − Ā₁, zero means).

Lemma D‴.26a: `CollapsedVariational.chang` gives the slaving exactly — the periodic P(t; ε) whose graph δe = εPY is invariant, by the
contraction P = 𝒢F(P) on a Fourier grid, for either sign of ε — and `slow_generator` the slow block's generator log Φ(ε)/(εT₀), the real
parts of whose eigenvalues are exactly minus the slow Floquet exponents (in slow time) and whose spectrum at −ε is that at ε (the parity as
an identity).
"""
import numpy as np
from scipy.integrate import solve_ivp

from .periodic import ClosedLoop
from .adaptive_prior import pair_rates


class AdaptiveRealLoop(ClosedLoop):
    """(CL) with the nearest-atom leak toward an adaptive prior; the state is (x, Θ, Z) for finite η and (x, Θ) for η = ∞; `labels` (m,)
    the frozen classes, K0 the number of atoms; `weights` (m,) the weights of the particles in f̂ (default 1/m; the weighted K₀-atom loop has
    one particle per class with the weights w)."""

    def __init__(self, net, f_vec, x_ref, x_ref_dot, period, ke, G_diag, sigma0, labels, K0, eta, weights=None):
        super().__init__(net, f_vec, x_ref, x_ref_dot, period, ke, G_diag, sigma0)
        self.labels = np.asarray(labels, int)
        self.K0 = int(K0)
        self.eta = float(eta)
        self.counts = np.bincount(self.labels, minlength=self.K0)
        assert np.all(self.counts > 0)
        self.weights = None if weights is None else np.asarray(weights, float)

    def f_hat(self, x, Th):
        if self.weights is None:
            return self.net.f_hat(x, Th)
        a, c, s = self.net.physical(Th)
        phi = np.exp(-np.sum((x[None, :] - c) ** 2, axis=1) / (2 * s ** 2))
        return (self.weights[:, None] * phi[:, None] * a).sum(axis=0)

    def class_means(self, Th):
        out = np.zeros((self.K0, self.p))
        np.add.at(out, self.labels, Th)
        return out / self.counts[:, None]

    def state_dim(self, m):
        return self.net.n + m * self.p + (0 if np.isinf(self.eta) else self.K0 * self.p)

    def rhs(self, t, z, m):
        n, p = self.net.n, self.p
        x, Th = z[:n], z[n:n + m * p].reshape(m, p)
        Z = self.class_means(Th) if np.isinf(self.eta) else z[n + m * p:].reshape(self.K0, p)
        xr, xrd = self.xr(t), self.xrd(t)
        e = x - xr
        fx = self.f_vec(x[None, :])[0]
        xdot = fx + xrd - self.ke * e - self.f_hat(x, Th)
        Thdot = self.G[None, :] * self.net.grad_g(x, e, Th) - (self.sigma0 / self.ke) * (Th - Z[self.labels])
        if np.isinf(self.eta):
            return np.concatenate([xdot, Thdot.ravel()])
        Zdot = -(self.eta / self.ke) * (Z - self.class_means(Th))
        return np.concatenate([xdot, Thdot.ravel(), Zdot.ravel()])

    def lift(self, atoms):
        """The state of perfect tracking at t = 0: x = x_m(0), the lift of the atoms by the classes, Z = the atoms."""
        atoms = np.asarray(atoms, float)
        parts = [np.atleast_1d(self.xr(0.0)), atoms[self.labels].ravel()]
        if not np.isinf(self.eta):
            parts.append(atoms.ravel())
        return np.concatenate(parts)


def predicted_slow_exponents(h, sigma0, eta, ke, m, K0, p):
    """The predicted slow Floquet exponents (rates per unit time, ascending) at a lifted exact representation: σ₀/k_e with multiplicity
    (m − K₀)p (the separations), and λ_j^±/k_e for the K₀p collective pairs (η = ∞: h_j/k_e only, the prior eliminated)."""
    h = np.asarray(h, float)
    lm, lp = pair_rates(h, sigma0, eta)
    coll = list(lm) if np.isinf(eta) else list(lm) + list(lp)
    sep = [sigma0] * ((m - K0) * p)
    return np.sort(np.array(coll + sep)) / ke


def period_map(cl, z0, m, rtol=1e-13, atol=1e-15):
    """One period of the loop from t = 0 (DOP853 at tight tolerances)."""
    sol = solve_ivp(lambda t, z: cl.rhs(t, z, m), (0.0, cl.T0), np.asarray(z0, float), method='DOP853', rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.y[:, -1]


def monodromy(cl, z_p, m, fd=1e-5, rtol=1e-13, atol=1e-15):
    """The monodromy matrix at the periodic state z_p by central differences of the period map, with the relative step fd·max(1, |z_j|)
    (truncation O(fd²); the integrator's error enters divided by the step, so a larger step and tight tolerances resolve the multipliers
    near 1 far better than lifted_orbit.monodromy_fd's 10⁻⁷ step at rtol 10⁻¹¹)."""
    d = len(z_p)
    DP = np.zeros((d, d))
    for j in range(d):
        h = fd * max(1.0, abs(z_p[j]))
        e = np.zeros(d); e[j] = h
        DP[:, j] = (period_map(cl, z_p + e, m, rtol, atol) - period_map(cl, z_p - e, m, rtol, atol)) / (2 * h)
    return DP



# ---------------------------------------------------------------------- the variational equation on the collapsed subspace
def whitened_jacobian(net, x, atoms, w, G):
    """Ĵ(x), (n, K₀p): the whitened weighted Jacobian of f̂ on the collapsed subspace at the K₀-atom configuration `atoms` (K₀, p) with the
    weights w, Ĵ[r, (k, l)] = √w_k √G_l ∂_{ϑ_kl}Φ̃_r(x; ϑ_k). With y_k = √w_k G^{−1/2}δϑ_k the linearised error injection is Ĵy and the
    learning law's linearisation at e = 0 is ẏ = Ĵᵀδe: the two are exactly adjoint."""
    x = np.atleast_1d(np.asarray(x, float))
    n = net.n
    atoms = np.asarray(atoms, float)
    sc = np.sqrt(np.asarray(w, float))[:, None] * np.sqrt(np.asarray(G, float))[None, :]
    out = np.empty((n, atoms.size))
    for r in range(n):
        e = np.zeros(n); e[r] = 1.0
        out[r] = (sc * net.grad_g(x, e, atoms)).ravel()
    return out


class CollapsedVariational:
    """The variational equation of (CL) along the lift of an exact representation ϑ^e(w) (perfect tracking, e ≡ 0, the prior at its atoms),
    restricted to the collapsed subspace — the weighted K₀-atom loop's linearisation — in the whitened weighted coordinates (δe, y, ζ):

        δė = −k_eδe − Ĵ(t)y,   ẏ = Ĵ(t)ᵀδe − (σ₀/k_e)(y − ζ),   ζ̇ = −(η/k_e)(ζ − y)      (η = ∞: ζ = y and the leak vanishes),

    and at η = 0 the fixed prior at the atoms (no ζ: ẏ = Ĵ(t)ᵀδe − (σ₀/k_e)y). Ĵ(t) = whitened_jacobian at x_m(t). It assumes ∂_xf̂ = ∂_xf on the orbit (f̂ = f as functions at an exact representation); otherwise the
    error equation also carries ∂_x(f − f̂)(x_m(t))δe. `x_ref` returns x_m(t) as an array of shape (n,)."""

    def __init__(self, net, x_ref, period, G, sigma0, eta, atoms, w):
        self.net, self.xr, self.T0 = net, x_ref, float(period)
        self.G = np.asarray(G, float)
        self.sigma0, self.eta = float(sigma0), float(eta)
        self.atoms = np.asarray(atoms, float)
        self.w = np.asarray(w, float)
        self.n = net.n
        self.d = self.atoms.size
        self.fixed = self.eta == 0.0                                          # the fixed prior at the atoms: no ζ, the leak −(σ₀/k_e)y

    def J(self, t):
        return whitened_jacobian(self.net, self.xr(t), self.atoms, self.w, self.G)

    @property
    def dim(self):
        return self.n + (self.d if (np.isinf(self.eta) or self.fixed) else 2 * self.d)

    def generator(self, t, ke):
        """The matrix of the variational equation at time t (real time)."""
        n, d = self.n, self.d
        Jt = self.J(t)
        M = np.zeros((self.dim, self.dim))
        M[:n, :n] = -ke * np.eye(n)
        M[:n, n:n + d] = -Jt
        M[n:n + d, :n] = Jt.T
        if self.fixed:
            M[n:, n:] = -(self.sigma0 / ke) * np.eye(d)
        elif not np.isinf(self.eta):
            I = np.eye(d)
            M[n:n + d, n:n + d] = -(self.sigma0 / ke) * I
            M[n:n + d, n + d:] = (self.sigma0 / ke) * I
            M[n + d:, n:n + d] = (self.eta / ke) * I
            M[n + d:, n + d:] = -(self.eta / ke) * I
        return M

    def monodromy(self, ke, rtol=1e-13, atol=1e-16):
        """Ψ(T₀) of Ψ̇ = M(t)Ψ, Ψ(0) = I (DOP853)."""
        D = self.dim
        f = lambda t, y: (self.generator(t, ke) @ y.reshape(D, D)).ravel()
        sol = solve_ivp(f, (0.0, self.T0), np.eye(D).ravel(), method='DOP853', rtol=rtol, atol=atol)
        if not sol.success:
            raise RuntimeError(sol.message)
        return sol.y[:, -1].reshape(D, D)

    def exponents(self, ke, **kw):
        """(the slow Floquet exponents — rates per unit time, ascending — and the error channel's): the slow ones from the multipliers of the
        monodromy (all but the n smallest, ≈ e^{−k_eT₀}, below the resolution), the error channel's by Liouville's formula, the exponents
        summing to −⟨tr M⟩ = n·k_e + (σ₀ + η)d/k_e (η = ∞: n·k_e; η = 0: n·k_e + σ₀d/k_e)."""
        lam = np.linalg.eigvals(self.monodromy(ke, **kw))
        lam = lam[np.argsort(-np.abs(lam))][:self.dim - self.n]
        slow = np.sort(-np.log(np.abs(lam)) / self.T0)
        total = self.n * ke + (0.0 if np.isinf(self.eta) else (self.sigma0 + self.eta) * self.d / ke)
        if self.fixed:
            total = self.n * ke + self.sigma0 * self.d / ke
        return slow, (total - slow.sum()) / self.n

    def averaged(self, N=1024):
        """(Ā₀, Ā₁, Ā₂, H, the filter matrix ⟨Ĵ̇ᵀĴ̇⟩) of the slow block in slow time, from N samples of Ĵ over a period (spectral
        derivatives and antiderivatives; Ĵ is analytic and periodic). H = ⟨ĴᵀĴ⟩ is the whitened weighted Gauss–Newton matrix h_j(w)."""
        T0, d, n = self.T0, self.d, self.n
        t = np.arange(N) * T0 / N
        J = np.array([self.J(s) for s in t])                                  # (N, n, d)
        k = np.fft.fftfreq(N, d=1.0 / N) * (2 * np.pi / T0)                   # angular frequencies

        def deriv(X, order):
            return np.real(np.fft.ifft(((1j * k) ** order)[:, None, None] * np.fft.fft(X, axis=0), axis=0))

        def antider(X):                                                       # the zero-mean antiderivative of a zero-mean array
            kk = k.copy(); kk[0] = 1.0
            Xf = np.fft.fft(X, axis=0) / (1j * kk)[:, None, None]; Xf[0] = 0.0
            return np.real(np.fft.ifft(Xf, axis=0))
        Jd, Jdd = deriv(J, 1), deriv(J, 2)
        if np.isinf(self.eta):
            E = np.eye(d); C = np.zeros((d, d))
        elif self.fixed:
            E = np.eye(d); C = -self.sigma0 * np.eye(d)
        else:
            I = np.eye(d)
            E = np.hstack([I, np.zeros((d, d))]); C = np.block([[-self.sigma0 * I, self.sigma0 * I], [self.eta * I, -self.eta * I]])
        Jb, Jbd, Jbdd = J @ E, Jd @ E, Jdd @ E
        JtJ = np.einsum('nri,nrj->nij', Jb, Jb)
        B0 = C[None] - JtJ
        B1 = np.einsum('nri,nrj->nij', Jb, Jbd)
        B2 = -np.einsum('nri,nrj->nij', Jb, Jbdd) + JtJ @ C[None] - JtJ @ JtJ
        A0 = B0.mean(axis=0)
        U1 = antider(B0 - A0[None])
        A1 = B1.mean(axis=0) + (B0 @ U1).mean(axis=0)
        U2 = antider(B1 + B0 @ U1 - U1 @ A0[None] - A1[None])
        A2 = B2.mean(axis=0) + (B1 @ U1).mean(axis=0) + (B0 @ U2).mean(axis=0)
        H = np.einsum('nri,nrj->nij', J, J).mean(axis=0)
        filt = np.einsum('nri,nrj->nij', Jd, Jd).mean(axis=0)
        return A0, A1, A2, H, filt

    def predicted(self, ke, order=2, N=1024, gens=None):
        """The collective exponents (rates per unit time, ascending) predicted by the averaged generator to the given order in 1/k_e:
        order 1: eig(−Ā₀)/k_e (= λ_j^±(w)/k_e); order 2: eig(−(Ā₀ + Ā₁/k_e + Ā₂/k_e²))/k_e."""
        A0, A1, A2 = (gens or self.averaged(N))[:3]
        A = A0 if order == 1 else A0 + A1 / ke + A2 / ke ** 2
        return np.sort(np.linalg.eigvals(-A).real) / ke

    def second_order(self, N=1024, gens=None):
        """(λ_j, c_j, the filter part of c_j): the slow-time rates are λ_j + c_j/k_e² + O(k_e⁻⁴) (Corollary D‴.26(c)), with
        c_j = −ℓ_jĀ₂r_j − Σ_{i≠j}(ℓ_jĀ₁r_i)(ℓ_iĀ₁r_j)/(λ_i − λ_j) (r_j, ℓ_j the right and left eigenvectors of −Ā₀, ℓ_jr_j = 1), and the error
        filter's part −ℓ_jE^T⟨Ĵ̇ᵀĴ̇⟩Er_j ≤ 0 (the attenuation of the regressor's harmonics by the error channel, k_e/(k_e + iω))."""
        A0, A1, A2, H, filt = gens or self.averaged(N)
        mu, R = np.linalg.eig(-A0)
        o = np.argsort(mu.real); lam = mu.real[o]; R = R[:, o].real
        Lm = np.linalg.inv(R)
        K = Lm @ A1 @ R; Q = Lm @ A2 @ R
        D = len(lam)
        c = np.array([-Q[j, j] - sum(K[j, i] * K[i, j] / (lam[i] - lam[j]) for i in range(D) if i != j) for j in range(D)])
        d = self.d
        E = np.eye(d) if (np.isinf(self.eta) or self.fixed) else np.hstack([np.eye(d), np.zeros((d, d))])
        F = Lm @ (E.T @ filt @ E) @ R
        return lam, c, -np.diag(F).real

    # ------------------------------------------------------------------ Lemma D‴.26a: the exact slow subspace (Chang's decoupling)
    def chang_grid(self, N=512):
        """(t, K, Λ, ω): N samples of K = [Ĵ, 0] (0 < η < ∞) or Ĵ (η = ∞, and η = 0), shape (N, n, q), on t_k = kT₀/N; the leak and the
        prior's constant block Λ in slow time (Λ = [[σ₀I, −σ₀I], [−ηI, ηI]]; η = ∞: 0; η = 0, the fixed prior: σ₀I); the angular frequencies of
        the grid."""
        t = np.arange(N) * self.T0 / N
        J = np.array([self.J(s) for s in t])
        if np.isinf(self.eta):
            K, Lam = J, np.zeros((self.d, self.d))
        elif self.fixed:
            K, Lam = J, self.sigma0 * np.eye(self.d)
        else:
            I = np.eye(self.d)
            K = np.concatenate([J, np.zeros_like(J)], axis=2)
            Lam = np.block([[self.sigma0 * I, -self.sigma0 * I], [-self.eta * I, self.eta * I]])
        om = np.fft.fftfreq(N, d=1.0 / N) * (2 * np.pi / self.T0)
        return t, K, Lam, om

    def chang(self, eps, N=512, mismatch=None, tol=1e-16, maxit=500):
        """Lemma D‴.26a(i): the T₀-periodic P(t; ε) whose graph w = PY is invariant under εẇ = −(1 − εa_mis)w − KY, Ẏ = ε(Kᵀw − ΛY) (δe = εw,
        ε = 1/k_e): εṖ + P = −K + εa_misP − ε²P(KᵀP − Λ), by the contraction P ← 𝒢F(P), 𝒢 the periodic solution of εu̇ + u = g (a division by
        1 + iεω on the Fourier grid; anti-causal for ε < 0, where the same equation defines the slow block at negative ε). `mismatch`: the
        samples of a derivative mismatch a_mis(t) in the error channel (default none). Returns (P (N, n, q), A = KᵀP − Λ (N, q, q), the
        iterations, the invariance residual max|εṖ + P − F(P)|)."""
        t, K, Lam, om = self.chang_grid(N)
        F = lambda P: (-K - eps ** 2 * np.einsum('tij,tjk->tik', P, np.einsum('tji,tjk->tik', K, P) - Lam[None])
                       + (0.0 if mismatch is None else eps * np.asarray(mismatch)[:, None, None] * P))
        P = -K.copy()
        for it in range(maxit):
            Pn = np.real(np.fft.ifft(np.fft.fft(F(P), axis=0) / (1 + 1j * eps * om)[:, None, None], axis=0))
            dP = np.abs(Pn - P).max(); P = Pn
            if dP < tol:
                break
        Pd = np.real(np.fft.ifft((1j * om)[:, None, None] * np.fft.fft(P, axis=0), axis=0))
        res = float(np.abs(eps * Pd + P - F(P)).max())
        A = np.einsum('tji,tjk->tik', K, P) - Lam[None]
        return P, A, it, res

    def slow_generator(self, eps, N=512, mismatch=None, rtol=1e-13):
        """Lemma D‴.26a(iii): Ā(ε) = log Φ(ε)/(εT₀), Φ the monodromy of the slow block Ẏ = εA(t)Y (A from `chang`, interpolated
        trigonometrically), whose eigenvalues' real parts are minus the slow-time rates; and the Magnus parameter εT₀max_t‖A‖₂ (the full
        Magnus series converges below π). Returns (Ā, the parameter)."""
        from scipy.linalg import logm
        t, K, Lam, om = self.chang_grid(N)
        A = self.chang(eps, N, mismatch)[1]
        q = A.shape[1]
        Ah = np.fft.fft(A, axis=0) / N
        Aat = lambda s: np.real(np.einsum('k,kij->ij', np.exp(1j * om * s), Ah))
        sol = solve_ivp(lambda s, y: (eps * Aat(s) @ y.reshape(q, q)).ravel(), (0.0, self.T0), np.eye(q).ravel(), method='DOP853',
                        rtol=rtol, atol=1e-16)
        Phi = sol.y[:, -1].reshape(q, q)
        return np.real(logm(Phi)) / (eps * self.T0), float(abs(eps) * self.T0 * max(np.linalg.norm(a, 2) for a in A))

