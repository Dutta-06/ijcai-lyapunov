"""
The stationary law of the dithered finite-m closed loop (PARADIGM.md Theorem J′)
================================================================================

Closed loop (CL_τ) of Theorem H's design at finite m, in the scaling σ = σ₀/k_e, τ = τ₀/k_e, ε := 1/k_e:

    ẋ = f(x) + ẋ_m − k_e e − f̂_Θ(x),    dϑ_i = [G∇_ϑ(eᵀΦ̃(x;ϑ_i)) − (σ₀/k_e)(ϑ_i − ϑ⁰)]dt + √(2τ₀/k_e) G^{1/2} dB_i.

In slow time s = εt the parameters follow, to leading order, the Langevin dynamics dΘ = −G_m∇J_m ds + √(2τ₀)G_m^{1/2}dW
whose invariant measure is the Gibbs measure π ∝ exp(−J_m/τ₀).  Theorem J′: the stationary Θ-marginal m_ε of the real
loop satisfies m_ε(φ) = π(φ) + ε∫φ g dπ + O(ε²), where g solves the weighted Poisson equation τ₀∇·(πG∇g) = ∇·(πb) with
the first-order drift b = b_fr + b_osc of Theorem J (Proposition H.7 + Remark H.7a) as a vector field on Θ-space; the
first-order stationary current is j₁ = π(b − τ₀G∇g), nonzero exactly when b is not a G-gradient (Theorem I), and the
circulation observable Λ := E[Θ ∧ k_e·(drift of Θ)] of the real loop equals ε[∫Θ∧j₁ + π(⟨u∧u̇⟩_θ)] + O(ε²).

This module implements the single-kernel instance used for the numerics (n = 1, m = 1, ϑ = (α, γ): saturated amplitude
and centre, frozen width), the averaged data (J, ∇J, b_fr, b_osc, ⟨u∧u̇⟩) on arbitrary sets of parameter points, the
grid solver for π, g, the Poisson equation, the current and the entropy-production coefficient, and a vectorised
simulator of the real loop (Strang splitting: exact Gaussian half-kicks around an RK4 step of the deterministic part)
that returns time averages of the observables with batch-mean standard errors.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


class SingleKernelInstance:
    """
    n = 1, m = 1.  ϑ = (α, γ) ∈ ℝ², a = a_max tanh(α/a_max), c = c_max tanh(γ/c_max), width s fixed:
    Φ̃(x; ϑ) = a φ(x; c, s), φ = exp(−(x−c)²/(2s²)).  Reference x_m(θ) = X[sin(ωθ′) + r₂ sin(2ωθ′)], θ′ = orientation·θ
    (T₀ = 2π/ω): the occupation measure does not depend on the orientation (±1), which reverses the direction of traversal.
    Plant f with derivative df (vectorised callables).  Leak σ₀ toward ϑ⁰ in the G⁻¹ metric; dither temperature τ₀.
    """

    def __init__(self, f, df, a_max=3.0, c_max=2.5, s=0.6, G=(1.0, 1.0), sigma0=0.2, theta0=(0.0, 0.0), tau0=0.1,
                 X=1.5, omega=1.0, r2=0.0, orientation=+1, N_theta=512):
        self.f, self.df = f, df
        self.a_max, self.c_max, self.s = float(a_max), float(c_max), float(s)
        self.G = np.asarray(G, float)
        self.sigma0, self.tau0 = float(sigma0), float(tau0)
        self.theta0 = np.asarray(theta0, float)
        self.X, self.omega, self.r2, self.orientation = float(X), float(omega), float(r2), int(orientation)
        self.T0 = 2 * np.pi / self.omega
        self.N_theta = int(N_theta)
        self.ths = np.arange(self.N_theta) * self.T0 / self.N_theta
        self.xm = self.x_ref(self.ths)
        self.xmd = self.x_ref_dot(self.ths)

    # ------------------------------------------------------------------ reference and saturations
    def x_ref(self, t):
        """x_m(t) = X[sin(ωt′) + r₂ sin(2ωt′)], t′ := orientation·t.  With r₂ ≠ 0 the reversed reference is not a time shift
        of the original (for r₂ = 0 it is, and every orientation-odd quantity vanishes identically)."""
        w = self.omega * self.orientation * np.asarray(t, float)
        return self.X * (np.sin(w) + self.r2 * np.sin(2 * w))

    def x_ref_dot(self, t):
        w = self.omega * self.orientation * np.asarray(t, float)
        return self.X * self.omega * self.orientation * (np.cos(w) + 2 * self.r2 * np.cos(2 * w))

    def phys(self, al, ga):
        """a, c and the derivatives a′ = da/dα, a″, c′, c″ (broadcasting)."""
        ta = np.tanh(al / self.a_max); tc = np.tanh(ga / self.c_max)
        a = self.a_max * ta; c = self.c_max * tc
        a1 = 1 - ta ** 2; c1 = 1 - tc ** 2
        a2 = -2 * ta * a1 / self.a_max; c2 = -2 * tc * c1 / self.c_max
        return a, c, a1, a2, c1, c2

    # ------------------------------------------------------------------ the real loop's vector field (physical time)
    def rhs(self, t, e, al, ga, ke):
        """(ė, α̇, γ̇) for arrays e, α, γ at physical time t (deterministic part, leak σ₀/k_e)."""
        xm = self.x_ref(t)
        x = xm + e
        a, c, a1, _, c1, _ = self.phys(al, ga)
        d = x - c
        phi = np.exp(-d ** 2 / (2 * self.s ** 2))
        fh = a * phi
        edot = self.f(x) - fh - ke * e
        # ∇_ϑ(e Φ̃) = e (a′ φ, a φ d/s² c′)
        aldot = self.G[0] * e * a1 * phi - (self.sigma0 / ke) * (al - self.theta0[0])
        gadot = self.G[1] * e * a * phi * d / self.s ** 2 * c1 - (self.sigma0 / ke) * (ga - self.theta0[1])
        return edot, aldot, gadot

    # ------------------------------------------------------------------ averaged data at parameter points
    def averaged(self, al, ga, with_hessian=False):
        """
        For arrays α, γ (any common shape) return dict with J (= J_1 = L̄ + leak), Lbar, gradJ (2,...), b_fr (2,...),
        b_osc (2,...), uu (⟨u∧u̇⟩_θ = ⟨u_α u̇_γ − u_γ u̇_α⟩, …), and optionally hessJ (2,2,...).  All θ-averages use the
        N_theta-point grid of the reference period (trapezoidal = exact for trigonometric polynomials of degree < N/2).
        """
        shape = np.broadcast(al, ga).shape
        al = np.broadcast_to(al, shape)[..., None]; ga = np.broadcast_to(ga, shape)[..., None]
        xm, xmd = self.xm, self.xmd                                       # (N,)
        a, c, a1, a2, c1, c2 = self.phys(al, ga)                            # (..., 1)
        d = xm - c                                                           # (..., N)
        s2 = self.s ** 2
        phi = np.exp(-d ** 2 / (2 * s2))
        fx, dfx = self.f(xm), self.df(xm)
        Phi = a * phi
        dPhi_x = -a * phi * d / s2                                           # ∂_xΦ̃
        r = fx - Phi                                                         # residual on the orbit
        Dr = dfx - dPhi_x                                                    # D_x r
        rdot = Dr * xmd                                                      # d/dθ r_Θ(x_m(θ))
        # ∇_ϑΦ̃ = (a′φ, a φ d/s² c′);  ∂_x∇_ϑΦ̃ = (a′ ∂_xφ, a c′ ∂_c∂_xφ) with ∂_xφ = −dφ/s², ∂_c∂_xφ = φ(1 − d²/s²)/s²
        gA = a1 * phi
        gC = a * phi * d / s2 * c1
        gxA = -a1 * phi * d / s2
        gxC = a * c1 * phi * (1 - d ** 2 / s2) / s2
        mean = lambda v: v.mean(axis=-1)
        Lbar = 0.5 * mean(r ** 2)
        lk = self.sigma0 / self.G
        J = Lbar + 0.5 * (lk[0] * (al[..., 0] - self.theta0[0]) ** 2 + lk[1] * (ga[..., 0] - self.theta0[1]) ** 2)
        # ∇J = −⟨r ∇_ϑΦ̃⟩ + σ₀G⁻¹(ϑ − ϑ⁰)
        gradJ = np.stack([-mean(r * gA) + lk[0] * (al[..., 0] - self.theta0[0]), -mean(r * gC) + lk[1] * (ga[..., 0] - self.theta0[1])])
        # b_fr = G ∇_{ϑ'}(U_even + U_odd)|_{ϑ'=ϑ}, residual frozen:
        #   U_odd(ϑ') = −⟨ṙ Φ̃(x_m;ϑ')⟩,  U_even(ϑ') = ⟨r (Dr Φ̃(x_m;ϑ') + r ∂_xΦ̃(x_m;ϑ'))⟩
        bfrA = self.G[0] * (-mean(rdot * gA) + mean(r * (Dr * gA + r * gxA)))
        bfrC = self.G[1] * (-mean(rdot * gC) + mean(r * (Dr * gC + r * gxC)))
        # H(ϑ,θ) = G r ∇_ϑΦ̃ − σ₀(ϑ−ϑ⁰) = −G∇_ϑ𝒥, 𝒥 = ½r² + (σ₀/2)|ϑ−ϑ⁰|²_{G⁻¹};  D_ϑH = −G∇²𝒥, ∇²𝒥 = ∇Φ̃∇Φ̃ᵀ − r∇²Φ̃ + σ₀G⁻¹
        HA = self.G[0] * r * gA - self.sigma0 * (al - self.theta0[0])
        HC = self.G[1] * r * gC - self.sigma0 * (ga - self.theta0[1])
        # ∇²_ϑΦ̃: αα: a″φ; αγ: a′c′ φ d/s²; γγ: a[c″ φ d/s² + c′² φ(−1/s² + d²/s⁴)]
        hAA = a2 * phi
        hAC = a1 * c1 * phi * d / s2
        hCC = a * (c2 * phi * d / s2 + c1 ** 2 * phi * (-1 / s2 + d ** 2 / s2 ** 2))
        JAA = gA * gA - r * hAA + lk[0]
        JAC = gA * gC - r * hAC
        JCC = gC * gC - r * hCC + lk[1]
        # zero-mean antiderivative in θ (FFT) of H − ⟨H⟩
        def zma(v):
            F = np.fft.rfft(v - v.mean(axis=-1, keepdims=True), axis=-1)
            k = np.arange(F.shape[-1]) * (2 * np.pi / self.T0)
            F[..., 1:] = F[..., 1:] / (1j * k[1:])
            F[..., 0] = 0
            return np.fft.irfft(F, n=self.N_theta, axis=-1)
        uA, uC = zma(HA), zma(HC)
        # b_osc = ⟨(D_ϑH) u⟩ = −G⟨∇²𝒥 u⟩
        boscA = -self.G[0] * mean(JAA * uA + JAC * uC)
        boscC = -self.G[1] * mean(JAC * uA + JCC * uC)
        HAd, HCd = HA - HA.mean(axis=-1, keepdims=True), HC - HC.mean(axis=-1, keepdims=True)    # u̇ = H − ⟨H⟩
        uu = mean(uA * HCd - uC * HAd)                                       # ⟨u_α u̇_γ − u_γ u̇_α⟩
        out = dict(J=J, Lbar=Lbar, gradJ=gradJ, b_fr=np.stack([bfrA, bfrC]), b_osc=np.stack([boscA, boscC]), uu=uu)
        if with_hessian:
            out['hessJ'] = np.stack([np.stack([mean(JAA), mean(JAC)]), np.stack([mean(JAC), mean(JCC)])])
        return out


class GibbsGrid:
    """
    Uniform grid on [α_lo, α_hi] × [γ_lo, γ_hi] with nα × nγ nodes.  Finite-volume discretisation of the weighted
    Laplacian 𝔇g := τ₀∇·(πG∇g) with zero flux through the boundary, face weights π_face = exp(−(J_i + J_j)/(2τ₀)).
    Solves 𝔇g = ∇·(πb) (first-order correction g = π₁/π) and the Poisson equation 𝔇ψ = π(φ − π(φ)) (L̄₀ψ = φ − π(φ)).
    """

    def __init__(self, inst, al_lim, ga_lim, n_al=161, n_ga=161):
        self.inst = inst
        self.al = np.linspace(al_lim[0], al_lim[1], n_al)
        self.ga = np.linspace(ga_lim[0], ga_lim[1], n_ga)
        self.ha, self.hg = self.al[1] - self.al[0], self.ga[1] - self.ga[0]
        self.AL, self.GA = np.meshgrid(self.al, self.ga, indexing='ij')         # (nα, nγ)
        self.shape = self.AL.shape
        self.data = inst.averaged(self.AL, self.GA, with_hessian=True)
        J = self.data['J']
        self.Jmin = J.min()
        w = np.exp(-(J - self.Jmin) / inst.tau0)
        self.cell = self.ha * self.hg
        self.Z = w.sum() * self.cell
        self.pi = w / self.Z                                                   # density on the grid
        self.b = self.data['b_fr'] + self.data['b_osc']                        # (2, nα, nγ)
        self._build_operator()
        self.g = self.solve_correction()

    # ------------------------------------------------------------------ operator
    def _build_operator(self):
        inst, J = self.inst, self.data['J']
        na, ng = self.shape
        N = na * ng
        idx = np.arange(N).reshape(na, ng)
        tau0, G = inst.tau0, inst.G
        rows, cols, vals = [], [], []
        # faces in α
        pf = np.exp(-(J[1:, :] + J[:-1, :] - 2 * self.Jmin) / (2 * tau0)) / self.Z           # π at α-faces, (na−1, ng)
        wA = tau0 * G[0] * pf / self.ha ** 2
        i0, i1 = idx[:-1, :].ravel(), idx[1:, :].ravel()
        for (r_, c_, v_) in ((i0, i0, -wA.ravel()), (i1, i1, -wA.ravel()), (i0, i1, wA.ravel()), (i1, i0, wA.ravel())):
            rows.append(r_); cols.append(c_); vals.append(v_)
        pfg = np.exp(-(J[:, 1:] + J[:, :-1] - 2 * self.Jmin) / (2 * tau0)) / self.Z
        wG = tau0 * G[1] * pfg / self.hg ** 2
        j0, j1 = idx[:, :-1].ravel(), idx[:, 1:].ravel()
        for (r_, c_, v_) in ((j0, j0, -wG.ravel()), (j1, j1, -wG.ravel()), (j0, j1, wG.ravel()), (j1, j0, wG.ravel())):
            rows.append(r_); cols.append(c_); vals.append(v_)
        self.D = sp.csc_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(N, N))
        self.pf_al, self.pf_ga = pf, pfg

    def divergence_pi_b(self, b):
        """∇·(πb) by the same finite-volume stencil: face values of πb are averages of the nodal values."""
        pb = self.pi[None] * b
        out = np.zeros(self.shape)
        fa = 0.5 * (pb[0, 1:, :] + pb[0, :-1, :])                              # flux through the α-face between i and i+1
        out[:-1, :] += fa / self.ha; out[1:, :] -= fa / self.ha                 # (∇·F)_i = (F_{i+1/2} − F_{i−1/2})/h
        fg = 0.5 * (pb[1, :, 1:] + pb[1, :, :-1])
        out[:, :-1] += fg / self.hg; out[:, 1:] -= fg / self.hg
        return out

    def _solve(self, rhs):
        """Solve 𝔇u = rhs (rhs must have zero sum; the constant mode is pinned at the mode of π) and return u with π(u) = 0."""
        N = self.D.shape[0]
        k = int(np.argmax(self.pi))
        A = self.D.tolil()
        A[k, :] = 0; A[k, k] = 1.0
        r = rhs.ravel().copy(); r[k] = 0.0
        u = spla.spsolve(A.tocsc(), r).reshape(self.shape)
        return u - self.expect(u)

    def solve_correction(self):
        return self._solve(self.divergence_pi_b(self.b))

    def poisson(self, phi):
        """ψ with L̄₀ψ = φ − π(φ), i.e. 𝔇ψ = π(φ − π(φ))."""
        return self._solve(self.pi * (phi - self.expect(phi)))

    # ------------------------------------------------------------------ expectations
    def expect(self, phi):
        return float(np.sum(phi * self.pi) * self.cell)

    def first_order(self, phi):
        """ε-coefficient of m_ε(φ): ∫φ g dπ."""
        return float(np.sum(phi * self.g * self.pi) * self.cell)

    def grad(self, u):
        gA = np.gradient(u, self.ha, axis=0); gG = np.gradient(u, self.hg, axis=1)
        return np.stack([gA, gG])

    def current(self):
        """j₁ = π(b − τ₀G∇g) (2, nα, nγ) and b_⊥ = b − τ₀G∇g."""
        bperp = self.b - self.inst.tau0 * self.inst.G[:, None, None] * self.grad(self.g)
        return self.pi[None] * bperp, bperp

    def epr_coefficient(self):
        """τ₀⁻¹ π(|b_⊥|²_{G⁻¹}): the ε² coefficient of the entropy production rate of the first-order averaged dynamics."""
        _, bperp = self.current()
        return self.expect(np.sum(bperp ** 2 / self.inst.G[:, None, None], axis=0)) / self.inst.tau0

    def circulation_prediction(self):
        """ε-coefficient of Λ = E[Θ ∧ k_e·drift]: ∫Θ∧j₁ + π(⟨u∧u̇⟩), with Θ∧v := α v_γ − γ v_α."""
        j1, _ = self.current()
        wedge = float(np.sum((self.AL * j1[1] - self.GA * j1[0])) * self.cell)
        return wedge + self.expect(self.data['uu']), wedge, self.expect(self.data['uu'])

    def stationary_density(self, eps):
        """
        Stationary density m of the first-order averaged SDE dΘ = (−G∇J + εb)ds + √(2τ₀)G^{1/2}dW on the grid, by the
        finite-volume discretisation ∇·(m v) − τ₀∇·(G∇m) = 0 with central face values of m and zero flux through the
        boundary (cell Péclet numbers ≲ 1 on the grids used), normalised to ∫m = 1; and its current j = m v − τ₀G∇m
        evaluated at the cell faces and averaged to the nodes.  Returns (m, j) with j of shape (2, nα, nγ).
        """
        inst = self.inst
        na, ng = self.shape
        N = na * ng
        idx = np.arange(N).reshape(na, ng)
        v = -inst.G[:, None, None] * self.data['gradJ'] + eps * self.b                      # drift (2, nα, nγ)
        tau0, G = inst.tau0, inst.G
        rows, cols, vals = [], [], []
        # α-faces between (i, ·) and (i+1, ·): flux F = v_face (m_i + m_j)/2 − τ₀G_α (m_j − m_i)/h, outward from i
        vf = 0.5 * (v[0, 1:, :] + v[0, :-1, :])
        i0, i1 = idx[:-1, :].ravel(), idx[1:, :].ravel()
        d = tau0 * G[0] / self.ha
        # (∇·F)_i gets +F/h, (∇·F)_j gets −F/h
        for (r_, c_, val) in ((i0, i0, (0.5 * vf.ravel() + d) / self.ha), (i0, i1, (0.5 * vf.ravel() - d) / self.ha),
                              (i1, i0, -(0.5 * vf.ravel() + d) / self.ha), (i1, i1, -(0.5 * vf.ravel() - d) / self.ha)):
            rows.append(r_); cols.append(c_); vals.append(val)
        vg = 0.5 * (v[1, :, 1:] + v[1, :, :-1])
        j0, j1 = idx[:, :-1].ravel(), idx[:, 1:].ravel()
        d = tau0 * G[1] / self.hg
        for (r_, c_, val) in ((j0, j0, (0.5 * vg.ravel() + d) / self.hg), (j0, j1, (0.5 * vg.ravel() - d) / self.hg),
                              (j1, j0, -(0.5 * vg.ravel() + d) / self.hg), (j1, j1, -(0.5 * vg.ravel() - d) / self.hg)):
            rows.append(r_); cols.append(c_); vals.append(val)
        A = sp.csc_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(N, N)).tolil()
        k = int(np.argmax(self.pi))
        A[k, :] = self.cell                                                                   # normalisation row ∫m = 1
        rhs = np.zeros(N); rhs[k] = 1.0
        m = spla.spsolve(A.tocsc(), rhs).reshape(self.shape)
        # current at the faces, averaged to nodes
        Fa = vf * 0.5 * (m[1:, :] + m[:-1, :]) - tau0 * G[0] * (m[1:, :] - m[:-1, :]) / self.ha
        Fg = vg * 0.5 * (m[:, 1:] + m[:, :-1]) - tau0 * G[1] * (m[:, 1:] - m[:, :-1]) / self.hg
        jA = np.zeros(self.shape); jG = np.zeros(self.shape)
        jA[1:-1, :] = 0.5 * (Fa[1:, :] + Fa[:-1, :]); jA[0, :] = 0.5 * Fa[0, :]; jA[-1, :] = 0.5 * Fa[-1, :]
        jG[:, 1:-1] = 0.5 * (Fg[:, 1:] + Fg[:, :-1]); jG[:, 0] = 0.5 * Fg[:, 0]; jG[:, -1] = 0.5 * Fg[:, -1]
        return m, np.stack([jA, jG])

    def sample(self, K, rng, eps=0.0):
        """K samples from the grid density π (eps = 0) or from the first-order law π(1 + εg), clipped at zero and
        renormalised (cell chosen by weight, uniform jitter inside the cell)."""
        dens = self.pi * (1 + eps * self.g) if eps else self.pi
        p = (np.maximum(dens, 0.0) * self.cell).ravel(); p = p / p.sum()
        k = rng.choice(p.size, size=K, p=p)
        ia, ig = np.unravel_index(k, self.shape)
        return self.al[ia] + (rng.random(K) - 0.5) * self.ha, self.ga[ig] + (rng.random(K) - 0.5) * self.hg


# ---------------------------------------------------------------------- the real loop, vectorised
OBS_NAMES = ('alpha', 'gamma', 'alpha2', 'gamma2', 'alpha_gamma', 'eta2', 'circ')


def observables(inst, t, e, al, ga, ke):
    """Pointwise observables of the full state: α, γ, α², γ², αγ, η² = (k_e e)², and the circulation
    α V_γ − γ V_α with V := k_e × (deterministic drift of (α, γ) in physical time)."""
    _, ad, gd = inst.rhs(t, e, al, ga, ke)
    VA, VG = ke * ad, ke * gd
    return np.stack([al, ga, al ** 2, ga ** 2, al * ga, (ke * e) ** 2, al * VG - ga * VA])


def simulate_stationary(inst, ke, K, T_slow, burn_slow, seed, al0, ga0, e0=None, steps_per_unit_ke=2.0):
    """
    Strang splitting in physical time with h = T₀/n_per, n_per = ceil(steps_per_unit_ke·k_e·T₀) (so k_e h ≤ 1/steps_per_unit_ke):
    half Gaussian kick √(τ₀h/k_e)G^{1/2}ξ on (α, γ), RK4 step of the deterministic (e, α, γ) system, half kick.
    Time-averages the observables over all steps after the burn-in, for K independent replicas started at (e0, al0, ga0);
    the burn-in and the averaging window are whole numbers of reference periods (at least burn_slow and T_slow slow-time
    units), so that the phase is averaged uniformly.  Returns (means, standard errors, n_steps): the mean is the average of
    the per-replica time averages and the standard error their standard deviation over √K (the replicas are independent;
    batch means in time would be correlated through the slowest mode of L̄₀).
    """
    rng = np.random.default_rng(seed)
    T0 = inst.T0
    n_per = int(np.ceil(steps_per_unit_ke * ke * T0))
    h = T0 / n_per
    e = np.zeros(K) if e0 is None else np.array(e0, float)
    al, ga = np.array(al0, float), np.array(ga0, float)
    kick = np.sqrt(inst.tau0 * h / ke) * np.sqrt(inst.G)
    n_burn = n_per * int(np.ceil(burn_slow * ke / T0)); n_avg = n_per * int(np.ceil(T_slow * ke / T0))
    acc = np.zeros((len(OBS_NAMES), K))
    t = 0.0
    for k in range(n_burn + n_avg):
        xi = rng.standard_normal((2, K))
        al += kick[0] * xi[0]; ga += kick[1] * xi[1]
        k1 = inst.rhs(t, e, al, ga, ke)
        k2 = inst.rhs(t + 0.5 * h, e + 0.5 * h * k1[0], al + 0.5 * h * k1[1], ga + 0.5 * h * k1[2], ke)
        k3 = inst.rhs(t + 0.5 * h, e + 0.5 * h * k2[0], al + 0.5 * h * k2[1], ga + 0.5 * h * k2[2], ke)
        k4 = inst.rhs(t + h, e + h * k3[0], al + h * k3[1], ga + h * k3[2], ke)
        e = e + h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        al = al + h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        ga = ga + h / 6 * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2])
        xi = rng.standard_normal((2, K))
        al += kick[0] * xi[0]; ga += kick[1] * xi[1]
        t = (k + 1) * h
        if k >= n_burn:
            acc += observables(inst, t, e, al, ga, ke)
    per_rep = acc / n_avg                                                     # per-replica time averages, (n_obs, K)
    mean = per_rep.mean(axis=1)
    se = per_rep.std(axis=1, ddof=1) / np.sqrt(K)
    return mean, se, n_burn + n_avg
