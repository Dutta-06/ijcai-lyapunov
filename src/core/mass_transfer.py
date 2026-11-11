"""
Proposition D‴.34 (PARADIGM.md §5d): mass transfer — the replicator weight law with the cohesion leak
=====================================================================================================

The particles carry weights β_i > 0, the network f̂ = (1/m)Σ_iβ_iΦ̃(·; ϑ_i) and its measure μ_{Θ,β} = (1/m)Σ_iβ_iδ_{ϑ_i}, and the weights
adapt by the replicator law

    dβ_i/ds = γβ_i(u_i − λ),   u_i := ⟨r Φ̃(·; ϑ_i)⟩_ρ,   λ := Σ_jβ_ju_j/Σ_jβ_j,                                    (𝔉)

(real time: β̇_i = γβ_ieᵀ(Φ̃(x; ϑ_i) − f̂(x)) when Σβ = m), which keeps Σ_iβ_i = m — the network's output bounded by the saturation as for
equal weights —, beside the transport law with the cohesion leak (Proposition D‴.32),

    dϑ_i/ds = β_iG⟨∇_ϑ(r Φ̃(·; ϑ_i))⟩_ρ − (σ₀/m)Σ_jψ(‖ϑ_i − ϑ_j‖²_{G⁻¹}/R²)(ϑ_i − ϑ_j).

The pair is the gradient flow of J(Θ, β) = mL̄(μ_{Θ,β}) + (σ₀m/2)∫∫W dμ_Θdμ_Θ (the cohesion energy of the unweighted cloud) in the metric
Σ_i‖δϑ_i‖²_{G⁻¹} + Σ_iδβ_i²/(γβ_i) restricted to Σ_iδβ_i = 0 (Euclidean on the positions, Fisher–Rao on the weights): ∂_{β_i}(mL̄) = −u_i and
∂_{ϑ_i}(mL̄) = −β_i⟨r∇_ϑΦ̃(·; ϑ_i)⟩_ρ. With β ≡ 1 and γ = 0 it is the loop of Proposition D‴.32. In real time with the positions fixed at
the atoms of an exact representation (weights β*, Σβ* = m) the replicator law has the Lyapunov function ½|e|² + (1/(γm))Σ_iβ*_i log(β*_i/β_i),
whose derivative is −k_e|e|².

This module: the weighted averaged loop on top of an `AveragedSlowLoop` (its quadrature, features and gradients), the energy, the run (the
midpoint rule on (Θ, log β)), the class masses, and the real-loop law with the KL Lyapunov function.
"""
import numpy as np

from . import cohesion_leak as CO


class ReplicatorCohesionLoop:
    """The cohesion leak (range R, kind, strength σ₀) with the replicator weight law of rate γ on top of an `AveragedSlowLoop` `base`;
    `balanced` keeps Σβ = m (λ the β-mean of u), otherwise λ = 0."""

    def __init__(self, base, R, kind='gauss', sigma0=None, gamma=1.0, balanced=True, eps=0.0, fixman=0.0, fixman_tau=0.0, fixman_A=None,
                 fixman_TA=None, box=None, box_sigma=1.0, fixman_floor=0.01, fixman_classes=None):
        self.base = base
        self.G = base.G
        self.R, self.kind = float(R), kind
        self.sigma0 = base.sigma0 if sigma0 is None else float(sigma0)
        self.gamma, self.balanced = float(gamma), bool(balanced)
        self.eps = float(eps)                     # Proposition D‴.36: the entropic leak on the masses, energy (ε/γ)Σβ_i log β_i
        self.fixman = float(fixman)                                       # Proposition D‴.38: the gain κ_F,
        self.fixman_tau = np.asarray(fixman_tau, float)                   # the dither's τ, or the classes' amplitude-noise levels τ_k (D‴.40)
        self.fixman_A = None if fixman_A is None else np.array(fixman_A, float)   # the classes' contributions Â_k (the masses form)
        self.fixman_TA = fixman_TA                                        # their running average's time constant (None: fixed)
        self.fixman_floor = float(fixman_floor)                           # the denominator's floor, a fraction of W_i
        self.fixman_classes = None if fixman_classes is None else np.asarray(fixman_classes, int)   # Proposition D‴.42: the particles'
        # classes; the correction is then the class's (W_k = Σ_{i∈k}β_i/m, Â_k per class), given to each of its particles
        self.box = None if box is None else np.asarray(box, float)          # Proposition D‴.39: the dead zone's half-widths (L_α, L_γ, L_κ)
        self.box_sigma = float(box_sigma)                                   # and its stiffness σ_B

    def f_hat(self, Th, beta):
        a, c, s, diff, phi = self.base.features(Th)
        return phi @ (beta * a) / Th.shape[0]

    def residual(self, Th, beta):
        return self.base.F - self.f_hat(Th, beta)

    def fitness(self, r, Th):
        """u_i = ⟨r Φ̃(·; ϑ_i)⟩_ρ, (m,)."""
        a, c, s, diff, phi = self.base.features(Th)
        return a * (phi.T @ r) / self.base.N

    def drift(self, Th, beta, xi=None):
        """The drift of (𝔉) with the transport law; `xi` (N,) is added to the residual (Proposition D‴.41: the measured error r + ξ drives
        both laws)."""
        r = self.residual(Th, beta)
        if xi is not None:
            r = r + xi
        dTh = beta[:, None] * self.G[None, :] * self.base.grad_weighted(r, Th) + CO.cohesion_force(Th, self.G, self.R, self.sigma0, self.kind)
        if self.box is not None:
            dTh = dTh - self.box_sigma * self.box_excess(Th)
        u = self.fitness(r, Th)
        if self.eps:                              # −∂_β of (ε/γ)Σβ log β, up to the multiplier
            u = u - (self.eps / self.gamma) * np.log(beta)
        if self.fixman:
            u = u + self.fixman_correction(Th, beta)
        lam = float(beta @ u / beta.sum()) if self.balanced else 0.0
        return dTh, self.gamma * beta * (u - lam)

    def fixman_correction(self, Th, beta):
        """Proposition D‴.38: κ_F(τ/2)∂_{β_i}ℓ along the manifold of exact fits, ℓ = log det ∇²_ΘJ = 2Σ log(W_i − W_min,i²/W_i) + const for
        one particle per class (Proposition D‴.37), W_i = β_i/m: κ_F(τ/m)(1 + W_min,i²/W_i²)/(W_i − W_min,i²/W_i). The masses form takes
        W_min,i = |Â_i|/a_max from the classes' contributions Â_i (`fixman_A`; on M, Â_i = W_i a_i = A_i), so it depends on the masses only;
        the plug-in form (no `fixman_A`) takes W_min,i = W_i|a_i|/a_max at the current amplitudes, i.e. κ_F(τ/m)(1 + (a_i/a_max)²)/(W_i a′_i) —
        equal on M, but off M it reads an amplitude's excursion as a mass at its saturation. Both vanish without dither and act on the masses
        only. With Â fixed the masses form is the gradient of a separable function of W (of −κ_F(τ/2)ℓ where the floor is inactive); the
        plug-in form depends on the amplitudes. The denominator is floored at `fixman_floor`·W_i (0.01: the cap binds within
        δ_cap = W_min,i(1/√0.99 − 1) of the saturation, where the correction stops growing — Proposition D‴.38(iii) with the floor)."""
        amax = self.base.net.a_max
        m = Th.shape[0]
        if self.fixman_classes is not None:                                 # Proposition D‴.42: per class
            lab = self.fixman_classes; K = int(lab.max()) + 1
            Wc = np.bincount(lab, beta, K) / m
            Ac = np.bincount(lab, beta * self.base.features(Th)[0], K) / m if self.fixman_A is None else self.fixman_A
            q = (np.abs(Ac) / amax / Wc) ** 2
            ck = self.fixman * (self.fixman_tau / m) * (1 + q) / np.maximum(Wc * (1 - q), self.fixman_floor * Wc)
            return ck[lab]
        W = beta / m
        if self.fixman_A is None:
            a = self.base.features(Th)[0]
            q = (a / amax) ** 2                                               # (W_min/W)² at the current amplitudes
        else:
            q = (np.abs(self.fixman_A) / amax / W) ** 2
        return self.fixman * (self.fixman_tau / m) * (1 + q) / np.maximum(W * (1 - q), self.fixman_floor * W)

    def box_excess(self, Th):
        """Proposition D‴.39: ϑ_i − Π_B(ϑ_i), the excess of the unconstrained coordinates over the box B = Π_j[−L_j, L_j]; the dead zone
        V = (σ_B/2)Σ_i‖ϑ_i − Π_B(ϑ_i)‖²_{G⁻¹} has the G-gradient σ_B(ϑ_i − Π_B(ϑ_i)), zero inside B."""
        return Th - np.clip(Th, -self.box[None, :], self.box[None, :])

    def J(self, Th, beta):
        m = Th.shape[0]
        r = self.residual(Th, beta)
        out = 0.5 * m * float(np.mean(r ** 2)) + 0.5 * self.sigma0 * m * CO.interaction_energy(Th, self.G, self.R, self.kind)
        if self.box is not None:
            out += 0.5 * self.box_sigma * float(np.sum(self.box_excess(Th) ** 2 / self.G[None, :]))
        if self.eps:
            out += (self.eps / self.gamma) * float(np.sum(beta * np.log(beta)))
        return out

    def fit_error(self, Th, beta):
        r = self.residual(Th, beta)
        return float(np.sqrt(np.mean(r ** 2) / np.mean(self.base.F ** 2)))

    def run(self, Th, beta, S_end, dt, rec_every=200, callback=None):
        """The midpoint rule in slow time on (Θ, log β) (the weights stay positive; with `balanced` they are renormalised to Σβ = m after
        each step, which the exact flow conserves); callback(s, Θ, β) at the records. Returns (Θ, β) at S_end."""
        Th, lb = np.array(Th, float), np.log(np.array(beta, float))
        m = Th.shape[0]
        n = int(round(S_end / dt))
        for k in range(n + 1):
            b = np.exp(lb)
            if callback is not None and k % rec_every == 0:
                callback(k * dt, Th, b)
            if k == n:
                break
            d1, g1 = self.drift(Th, b)
            Tm, lm = Th + 0.5 * dt * d1, lb + 0.5 * dt * g1 / b
            d2, g2 = self.drift(Tm, np.exp(lm))
            Th, lb = Th + dt * d2, lb + dt * g2 / np.exp(lm)
            if self.balanced:
                lb = lb - np.log(np.exp(lb).sum() / m)
        return Th, np.exp(lb)


    def run_dithered(self, Th, beta, S_end, dt, tau, rng, rec_every=4000, callback=None, mask=None, noise_L=None, label_v=0.0,
                     label_scheme='ito', label_rho=0.0):
        """Proposition D‴.36: the positions dithered (Euler–Maruyama, the noise √(2τG) per particle), the weights by the midpoint-free Euler
        step of (𝔉) on log β, renormalised to Σβ = m; callback(s, Θ, β) at the records. `mask` (p,) scales the noise per coordinate
        (Proposition D‴.40: (0, 1, 1) dithers the centres and widths only); the random draws are the same with or without it. `noise_L`
        (mp × mp) replaces the noise by √dt·L ξ, covariance Σ = LLᵀ per unit time on the flattened positions (τ and mask then unused).
        Proposition D‴.41: `label_v` > 0 adds measurement noise ξ to the residual that drives both laws, independent over the quadrature
        nodes, of intensity v (⟨ξh⟩_ρ has the variance v⟨h²⟩_ρ per unit time): 'ito' draws it afresh at each step (variance vN/dt), 'ar1'
        correlates successive draws by ρ = `label_rho` at the same intensity (per-step variance vN(1 − ρ)/((1 + ρ)dt), started stationary),
        'heun' is Heun's predictor–corrector with the same draw in both stages (the Stratonovich interpretation; with `label_rho` ≠ 0 the
        draws are the AR(1) sequence, a coloured noise integrated as an ordinary differential equation), the predictor renormalised to
        Σβ = m like every step when `balanced`."""
        Th, lb = np.array(Th, float), np.log(np.array(beta, float))
        m = Th.shape[0]
        sq = np.sqrt(2.0 * tau * dt * self.G)[None, :]
        if mask is not None:
            sq = sq * np.asarray(mask, float)[None, :]
        if noise_L is not None:
            Ld = np.sqrt(dt) * np.asarray(noise_L, float)
        if label_v:
            if label_scheme not in ('ito', 'ar1', 'heun'):
                raise ValueError(label_scheme)
            if label_scheme == 'heun' and self.fixman_A is not None and self.fixman_TA:
                raise NotImplementedError("the tracked contributions with Heun's scheme")
            rho = float(label_rho) if label_scheme in ('ar1', 'heun') else 0.0
            s_xi = np.sqrt(label_v * self.base.N / dt * (1 - rho) / (1 + rho))
            xi = s_xi * rng.standard_normal(self.base.N)
        n = int(round(S_end / dt))
        for k in range(n + 1):
            b = np.exp(lb)
            if callback is not None and k % rec_every == 0:
                callback(k * dt, Th, b)
            if k == n:
                break
            if label_v and k > 0:
                xi = rho * xi + np.sqrt(1 - rho ** 2) * s_xi * rng.standard_normal(self.base.N)
            dT, db = self.drift(Th, b, xi if label_v else None)
            if self.fixman_A is not None and self.fixman_TA:                  # the contributions' running average (Proposition D‴.38)
                ba = b / m * self.base.features(Th)[0]
                if self.fixman_classes is not None:
                    ba = np.bincount(self.fixman_classes, ba, len(self.fixman_A))
                self.fixman_A = self.fixman_A + (dt / self.fixman_TA) * (ba - self.fixman_A)
            dW = sq * rng.standard_normal(Th.shape) if noise_L is None else (Ld @ rng.standard_normal(Th.size)).reshape(Th.shape)
            if label_v and label_scheme == 'heun':
                Tp, lp = Th + dt * dT + dW, lb + dt * db / b
                if self.balanced:
                    lp = lp - np.log(np.exp(lp).sum() / m)
                dT2, db2 = self.drift(Tp, np.exp(lp), xi)
                Th, lb = Th + 0.5 * dt * (dT + dT2) + dW, lb + 0.5 * dt * (db / b + db2 / np.exp(lp))
            else:
                Th = Th + dt * dT + dW
                lb = lb + dt * db / b
            if self.balanced:
                lb = lb - np.log(np.exp(lb).sum() / m)
        return Th, np.exp(lb)


def class_masses(Th, beta, G, thr):
    """The clusters of Θ (single linkage at thr) with their masses Σ_{i∈k}β_i/m and their β-weighted atoms: (labels, masses (K,), atoms (K, p))."""
    lab = CO.clusters(Th, G, thr)
    K = int(lab.max() + 1)
    m = Th.shape[0]
    W = np.array([beta[lab == k].sum() / m for k in range(K)])
    A = np.array([(beta[lab == k, None] * Th[lab == k]).sum(axis=0) / beta[lab == k].sum() for k in range(K)])
    return lab, W, A


def replicator_real_rhs(t, z, net, f_vec, x_ref, x_ref_dot, ke, Th, gamma, m):
    """The real closed loop with the positions Θ frozen and the replicator law on the weights: z = (x, β); ẋ = f(x) + ẋ_m − k_ee − f̂(x),
    β̇_i = γβ_ie(Φ̃(x; ϑ_i) − f̂(x)) (n = 1, Σβ = m conserved)."""
    x, beta = z[:1], z[1:]
    e = x - x_ref(t)
    a, c, s = net.physical(Th)
    phi = np.exp(-np.sum((x[None, :] - c) ** 2, axis=1) / (2 * s ** 2))
    Phi = phi * a[:, 0]
    fh = float(beta @ Phi) / m
    xdot = f_vec(x[None, :])[0] + x_ref_dot(t) - ke * e - fh
    return np.concatenate([xdot, gamma * beta * float(e[0]) * (Phi - fh)])


def kl_lyapunov(z, x_ref, t, beta_star, gamma, m):
    """½e² + (1/(γm))Σ_iβ*_i log(β*_i/β_i) (Σβ = Σβ* = m)."""
    e = z[0] - x_ref(t)[0]
    beta = z[1:]
    return 0.5 * float(e ** 2) + float(np.sum(beta_star * np.log(beta_star / beta))) / (gamma * m)
