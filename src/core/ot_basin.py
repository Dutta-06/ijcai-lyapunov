"""
Theorem D‴.21 (PARADIGM.md §5d): a measure-level basin for the barycentric OT leak — a Polyak–Łojasiewicz inequality uniform in m
==============================================================================================================================

For m commensurate with ν = Σ_j q_jδ_{ϑ⁰_j} every optimal plan of μ_Θ → ν is an assignment a with |a⁻¹(j)| = mq_j (Proposition D‴.12), and

    J_m^{OT} = min_a J_a,    J_a(Θ) := mL̄(μ_Θ) + (σ₀/2)Σ_i‖ϑ_i − ϑ⁰_{a(i)}‖²_{G⁻¹}.

Against a replicated stationary measure μ* = Σ_j q_jδ_{ϑ*_j} (every atom critical for its class potential, the diagonal plan optimal; residual
r*, energy 𝒥* := L̄(μ*) + (σ₀/2)Σ_jq_j‖ϑ*_j − ϑ⁰_j‖²_{G⁻¹}) write u_i := ϑ_i − ϑ*_{a(i)}, D_a := (1/m)Σ_i‖u_i‖²_{G⁻¹} and
Δf̂ := f̂_Θ − f̂[μ*] = (1/m)Σ_i[Φ̃(·;ϑ_i) − Φ̃(·;ϑ*_{a(i)})] — the last equality is where the counts mq_j enter (with other counts both
identities acquire the defects of `counts_mismatch`; (E) holds for every coupling with the column masses q, (P) needs an assignment, one
test direction per particle). Two identities are exact:

  (E)  E_a := J_a(Θ)/m − 𝒥* = ½‖Δf̂‖² + (σ₀/2)D_a − ⟨r*, (1/m)Σ_i R₂,i⟩,              R₂,i := Φ̃(ϑ_i) − Φ̃(ϑ*_{a(i)}) − ∇Φ̃(ϑ*_{a(i)})u_i,
  (P)  Q_a := (1/m)Σ_i⟨∇_{ϑ_i}J_a, u_i⟩
            = ‖Δf̂‖² + σ₀D_a − ⟨r*, (1/m)Σ_i(∇Φ̃(ϑ_i) − ∇Φ̃(ϑ*_{a(i)}))u_i⟩ + ⟨Δf̂, (1/m)Σ_i R̃₂,i⟩,    R̃₂,i := ∇Φ̃(ϑ_i)u_i − Φ̃(ϑ_i) + Φ̃(ϑ*_{a(i)}),

⟨·,·⟩ the L²(ρ) pairing; the remainders are the integral forms of the plan-segment curvatures,
⟨r*, R₂,i⟩ = ∫₀¹(1−t)u_iᵀR[μ*](ϑ*_{a(i)} + tu_i)u_i dt and the last two terms of (P) = −∫₀¹u_iᵀ[(1−t)R[μ*] + tR[μ_Θ]](ϑ*_{a(i)} + tu_i)u_i dt.
Cauchy–Schwarz gives |∇J_a|²_G ≥ mQ_a²/D_a, and with the plan-segment margin α := σ₀ − (κ* + κ^Θ)/2 > 0 (κ*, κ^Θ the largest eigenvalue of
the whitened curvature at the stationary and at the current residual along the segments [ϑ*_{a(i)}, ϑ_i]) the PL inequality
|∇J_a|²_G ≥ λ·mE_a with λ = α·min{2, 2α/(σ₀ + κ*_−)}.

`BasinDiagnostics` computes both sides of (E) and (P), the exact effective margins (κ_E, κ_Q: E_a = ½‖Δf̂‖² + ((σ₀ − κ_E)/2)D_a and
Q_a = ‖Δf̂‖² + (σ₀ − κ_Q)D_a), the observed PL ratio |∇J_a|²_G/(mE_a) and its lower bounds Q_a²/(D_aE_a) ≥ λ, and the segment curvatures on a
grid of t (a lower estimate of the suprema).
"""
import numpy as np

from .chaos_learning import _sig


def feature_tables(loop, Th):
    """Φ̃(x; ϑ_i) (N, m) and ∇_ϑΦ̃(x; ϑ_i) (N, m, 3) on the quadrature grid, in the unconstrained coordinates (chain rule through the
    saturations, as in AveragedSlowLoop.grad_weighted)."""
    a, c, s, diff, phi = loop.features(Th)
    net = loop.net
    da = 1.0 - (a / net.a_max) ** 2
    dc = 1.0 - (c / net.c_max) ** 2
    sk = _sig(Th[:, 2]); ds = (net.s_max - net.s_min) * sk * (1 - sk)
    Phi = phi * a[None, :]
    dPhi = np.stack([phi * da[None, :], Phi * diff / s[None, :] ** 2 * dc[None, :], Phi * diff ** 2 / s[None, :] ** 3 * ds[None, :]], axis=2)
    return Phi, dPhi


class BasinDiagnostics:
    """The identities (E), (P) and the PL constants of Theorem D‴.21 for an `ot_leak.OTPriorLoop` at commensurate sizes, against the
    replicated stationary measure with atoms `atoms` (K₀ × 3) and the loop's masses q."""

    def __init__(self, loop, atoms):
        self.loop = loop
        self.atoms = np.asarray(atoms, float)
        self.q = loop.q
        self.G = loop.G
        self.sigma0 = loop.sigma0
        Phi_s, dPhi_s = feature_tables(loop, self.atoms)
        self.Phi_s, self.dPhi_s = Phi_s, dPhi_s
        self.fhat_s = Phi_s @ self.q                                         # f̂[μ*] on the grid
        self.r_s = loop.F - self.fhat_s
        self.Jstar = 0.5 * float(np.mean(self.r_s ** 2)) + 0.5 * self.sigma0 * float(np.sum(self.q * np.sum((self.atoms - loop.prior) ** 2 / self.G[None, :], axis=1)))
        # stationarity of every atom for its class potential: ∇V_L[μ*](ϑ*_j) + σ₀G⁻¹(ϑ*_j − ϑ⁰_j)
        self.stationarity = float(np.abs(loop.grad_V_L(self.r_s, self.atoms) + self.sigma0 * (self.atoms - loop.prior) / self.G[None, :]).max())

    # ------------------------------------------------------------------ plans
    def labels(self, Th):
        """The assignment a of an optimal (integral) plan μ_Θ → ν; m must be commensurate."""
        assert self.loop.commensurate(Th.shape[0]), "Theorem D‴.21 is stated for commensurate m"
        gamma, _ = self.loop.plan(Th)
        return np.argmax(gamma, axis=1)

    # ------------------------------------------------------------------ the identities
    def terms(self, Th, lab=None, check_counts=True):
        """Both sides of (E) and (P) for the assignment `lab` (default: an optimal one), the gradient norm and the PL ratios. The identities
        need the counts mq (check_counts=False evaluates both sides anyway: with counts N ≠ mq, (E) is off by Σ_k(N_k/m − q_k)c_k with
        c_k = −⟨r*, Φ̃(ϑ*_k)⟩ + (σ₀/2)‖ϑ*_k − ϑ⁰_k‖²_{G⁻¹}, and (P) by −⟨Δf̂, f̂[μ_N] − f̂[μ*]⟩, μ_N := Σ_k(N_k/m)δ_{ϑ*_k})."""
        Th = np.asarray(Th, float)
        m = Th.shape[0]
        lab = self.labels(Th) if lab is None else np.asarray(lab)
        counts = np.bincount(lab, minlength=self.atoms.shape[0])
        if check_counts:
            assert np.all(counts == np.round(m * self.q)), "the assignment must have the counts mq"
        G, s0 = self.G, self.sigma0
        u = Th - self.atoms[lab]
        D = float(np.mean(np.sum(u ** 2 / G[None, :], axis=1)))
        Phi, dPhi = feature_tables(self.loop, Th)
        Phi_s, dPhi_s = self.Phi_s[:, lab], self.dPhi_s[:, lab, :]
        fhat = Phi.mean(axis=1)
        r = self.loop.F - fhat
        df = fhat - self.fhat_s
        ndf2 = float(np.mean(df ** 2))
        # direct evaluations
        leak = 0.5 * s0 * float(np.mean(np.sum((Th - self.loop.prior[lab]) ** 2 / G[None, :], axis=1)))
        E_direct = 0.5 * float(np.mean(r ** 2)) + leak - self.Jstar
        grad = self.loop.grad_V_L(r, Th) + s0 * (Th - self.loop.prior[lab]) / G[None, :]          # ∇_{ϑ_i}J_a, (m, 3)
        Q_direct = float(np.mean(np.sum(grad * u, axis=1)))
        gnorm2 = float(np.sum(G[None, :] * grad ** 2))                                             # |∇J_a|²_G
        # the identities
        du_s = np.einsum('nmc,mc->nm', dPhi_s, u)                                                   # ∇Φ̃(ϑ*_{a(i)})u_i
        du = np.einsum('nmc,mc->nm', dPhi, u)                                                       # ∇Φ̃(ϑ_i)u_i
        R2 = (Phi - Phi_s - du_s).mean(axis=1)
        R2t = (du - Phi + Phi_s).mean(axis=1)
        dgrad = (du - du_s).mean(axis=1)
        E_ident = 0.5 * ndf2 + 0.5 * s0 * D - float(np.mean(self.r_s * R2))
        Q_ident = ndf2 + s0 * D - float(np.mean(self.r_s * dgrad)) + float(np.mean(df * R2t))
        kE = 2.0 * float(np.mean(self.r_s * R2)) / D if D > 0 else 0.0                            # E = ½‖Δf̂‖² + ((σ₀ − κ_E)/2)D
        kQ = (float(np.mean(self.r_s * dgrad)) - float(np.mean(df * R2t))) / D if D > 0 else 0.0  # Q = ‖Δf̂‖² + (σ₀ − κ_Q)D
        E = E_direct
        ratio = gnorm2 / (m * E) if E > 0 else float('nan')
        lam_eff = Q_direct ** 2 / (D * E) if (D > 0 and E > 0) else float('nan')
        return dict(m=m, D=D, ndf2=ndf2, E=E, E_ident=E_ident, Q=Q_direct, Q_ident=Q_ident, gnorm2=gnorm2, kappa_E=kE, kappa_Q=kQ,
                    ratio=ratio, lam_eff=lam_eff, lab=lab, u=u, r=r, df=df, counts=counts)

    def counts_mismatch(self, counts, df):
        """The defects of (E) and (P) for an assignment with the counts N (zero iff N = mq): (Σ_k(N_k/m − q_k)c_k, −⟨Δf̂, f̂[μ_N] − f̂[μ*]⟩)."""
        counts = np.asarray(counts, float); m = counts.sum()
        dw = counts / m - self.q
        c = -np.mean(self.r_s[:, None] * self.Phi_s, axis=0) + 0.5 * self.sigma0 * np.sum((self.atoms - self.loop.prior) ** 2 / self.G[None, :], axis=1)
        return float(dw @ c), -float(np.mean(df * (self.Phi_s @ dw)))

    # ------------------------------------------------------------------ the plan-segment curvatures and the theorem's rate
    def segment_curvatures(self, Th, lab=None, n_t=9, r=None):
        """κ* = max λ_max(R^G[μ*]), κ*_− = max (−λ_min(R^G[μ*]))₊ and κ^Θ = max λ_max(R^G[μ_Θ]) over the points ϑ*_{a(i)} + t u_i,
        t on n_t equispaced points of [0, 1] (a lower estimate of the suprema over the segments)."""
        Th = np.asarray(Th, float)
        lab = self.labels(Th) if lab is None else np.asarray(lab)
        u = Th - self.atoms[lab]
        ts = np.linspace(0.0, 1.0, n_t)
        P = (self.atoms[lab][None, :, :] + ts[:, None, None] * u[None, :, :]).reshape(-1, 3)
        P = np.unique(np.round(P, 13), axis=0)
        r = self.loop.residual(Th) if r is None else r
        sq = np.sqrt(self.G)
        out = {}
        for key, res in (('star', self.r_s), ('theta', r)):
            R = self.loop.residual_curvature(res, P)
            Rw = sq[None, :, None] * R * sq[None, None, :]
            ev = np.linalg.eigvalsh(Rw)
            out[key] = (float(ev[:, -1].max()), float(max(0.0, -ev[:, 0].min())))
        return dict(kappa_star=out['star'][0], kappa_star_minus=out['star'][1], kappa_theta=out['theta'][0])

    def theorem_rate(self, kappa_star, kappa_star_minus, kappa_theta):
        """α = σ₀ − (κ* + κ^Θ)/2 and λ = α·min{2, 2α/(σ₀ + κ*_−)} (0 when α ≤ 0)."""
        alpha = self.sigma0 - 0.5 * (kappa_star + kappa_theta)
        lam = alpha * min(2.0, 2.0 * alpha / (self.sigma0 + kappa_star_minus)) if alpha > 0 else 0.0
        return alpha, lam
