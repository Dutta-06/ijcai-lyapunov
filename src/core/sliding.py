"""
Second-order plants: the sliding-variable reduction and the averaged slow loop on the phase plane
==================================================================================================

The plant ẍ = f(x, ẋ) + u (scalar x, unknown f) tracks a T₀-periodic reference x_m with the filtered error s := ė + λe (e := x − x_m, λ > 0)
and the controller u = ẍ_m − λė − k_s s − f̂_Θ(x, ẋ). Then

    ṡ = −k_s s + r_Θ(z),    ė = −λe + s,    r_Θ := f − f̂_Θ,    z := (x, ẋ) = z_m + (e, ė),    z_m := (x_m, ẋ_m),

the error channel of Theorem J with s in place of e and the phase-plane point z in place of x; the adaptive law
ϑ̇_i = G∇_ϑ(sΦ̃(z; ϑ_i)) − (σ₀/k_s)(ϑ_i − T_i) is the law of Theorem J. With k_s → ∞ at fixed λ both s and e are fast (s = r_Θ(z_m)/k_s + O(k_s⁻²),
e = O(1/k_s)), and in the slow time τ = t/k_s the parameters follow the averaged slow loop of Theorem J along the phase-plane orbit,

    dϑ_i/dτ = G⟨r_Θ∇_ϑΦ̃(·; ϑ_i)⟩_ρ − σ₀(ϑ_i − T_i),    ρ the occupation measure of z_m over one period,

the same objects as in one dimension with a curve in the plane in place of an interval. This module gives the network on the phase plane
(saturated Gaussian kernels with a scalar amplitude, ϑ = (α, γ₁, γ₂, κ), p = 4), the averaged loop with the interface the adaptive-prior and
OT machinery uses (drift(Θ, T), L_bar, G, σ₀), the residual curvature, the whitened weighted Gauss–Newton block, and the real closed loop.
"""
import numpy as np
from scipy.special import expit as _sig


class PhasePlaneNetwork:
    """Φ̃(z; ϑ) = a φ(z; c, s), z ∈ ℝ^d, a = a_max tanh(α/a_max), c_j = c_max,j tanh(γ_j/c_max,j), s = s_min + (s_max − s_min)sig(κ);
    ϑ = (α, γ₁, …, γ_d, κ), p = d + 2. c_max may differ per coordinate (position and velocity have different scales)."""

    def __init__(self, d=2, a_max=4.0, c_max=(3.0, 3.0), s_min=0.2, s_max=2.0):
        self.d = int(d)
        self.p = self.d + 2
        self.a_max = float(a_max)
        self.c_max = np.broadcast_to(np.asarray(c_max, float), (self.d,)).copy()
        self.s_min, self.s_max = float(s_min), float(s_max)

    def physical(self, Th):
        Th = np.asarray(Th, float)
        a = self.a_max * np.tanh(Th[:, 0] / self.a_max)
        c = self.c_max[None, :] * np.tanh(Th[:, 1:1 + self.d] / self.c_max[None, :])
        s = self.s_min + (self.s_max - self.s_min) * _sig(Th[:, 1 + self.d])
        return a, c, s

    def unconstrained(self, a, c, s):
        a = np.asarray(a, float); c = np.atleast_2d(np.asarray(c, float)); s = np.asarray(s, float)
        al = self.a_max * np.arctanh(np.clip(a / self.a_max, -1 + 1e-12, 1 - 1e-12))
        ga = self.c_max[None, :] * np.arctanh(np.clip(c / self.c_max[None, :], -1 + 1e-12, 1 - 1e-12))
        u = np.clip((s - self.s_min) / (self.s_max - self.s_min), 1e-12, 1 - 1e-12)
        return np.column_stack([al, ga, np.log(u / (1 - u))])

    def chain(self, Th):
        """da/dα, dc/dγ (m, d), ds/dκ."""
        a, c, s = self.physical(Th)
        sk = _sig(np.asarray(Th, float)[:, 1 + self.d])
        return 1.0 - (a / self.a_max) ** 2, 1.0 - (c / self.c_max[None, :]) ** 2, (self.s_max - self.s_min) * sk * (1 - sk)


class AveragedSlidingLoop:
    """The averaged slow loop along the phase-plane orbit z_m(t) = (x_m(t), ẋ_m(t)): plant f (callable on (N, 2) arrays → (N,)), the reference
    and its derivative, period, gains G (p,), leak σ₀; ρ the uniform grid of n_quad times of one period. f̂_Θ = Σ_iw_iΦ̃(·; ϑ_i), w_i = 1/m by
    default."""

    def __init__(self, net, f, x_ref, x_ref_dot, period, G, sigma0, n_quad=400):
        self.net = net
        self.p = net.p
        self.G = np.asarray(G, float)
        assert self.G.shape == (self.p,)
        self.sigma0 = float(sigma0)
        self.period = float(period)
        self.N = int(n_quad)
        ts = np.arange(self.N) * self.period / self.N
        self.Zq = np.column_stack([[float(x_ref(t)) for t in ts], [float(x_ref_dot(t)) for t in ts]])     # (N, 2)
        self.F = np.asarray(f(self.Zq), float).reshape(self.N)
        self.F2 = float(np.mean(self.F ** 2))

    def features(self, Th):
        """(a (m,), c (m, d), s (m,), diff (N, m, d), sq (N, m), φ (N, m))."""
        a, c, s = self.net.physical(Th)
        diff = self.Zq[:, None, :] - c[None, :, :]
        sq = np.sum(diff ** 2, axis=2)
        return a, c, s, diff, sq, np.exp(-sq / (2 * s[None, :] ** 2))

    def weights(self, Th, w):
        return np.full(np.asarray(Th).shape[0], 1.0 / np.asarray(Th).shape[0]) if w is None else np.asarray(w, float)

    def f_hat(self, Th, w=None):
        a, c, s, diff, sq, phi = self.features(Th)
        return phi @ (self.weights(Th, w) * a)

    def residual(self, Th, w=None):
        return self.F - self.f_hat(Th, w)

    def L_bar(self, Th, w=None):
        """½⟨r²⟩_ρ."""
        return 0.5 * float(np.mean(self.residual(Th, w) ** 2))

    def fit_error(self, Th, w=None):
        return float(np.sqrt(np.mean(self.residual(Th, w) ** 2) / self.F2))

    def grad_features(self, Th):
        """∇_ϑΦ̃(z_n; ϑ_i) in the unconstrained coordinates, (N, m, p)."""
        a, c, s, diff, sq, phi = self.features(Th)
        da, dc, ds = self.net.chain(Th)
        d = self.net.d
        out = np.empty((self.N, Th.shape[0], self.p))
        out[:, :, 0] = phi * da[None, :]
        out[:, :, 1:1 + d] = (a * phi)[:, :, None] * diff / (s ** 2)[None, :, None] * dc[None, :, :]
        out[:, :, 1 + d] = a[None, :] * phi * sq / s[None, :] ** 3 * ds[None, :]
        return out

    def learning(self, Th, w=None, r=None):
        """⟨r∇_ϑΦ̃(·; ϑ_i)⟩_ρ, (m, p); r the residual of Θ with the weights w unless given."""
        Th = np.asarray(Th, float)
        if r is None:
            r = self.residual(Th, w)
        return np.einsum('n,nip->ip', r, self.grad_features(Th)) / self.N

    def grad_V_L(self, r, Th):
        """∇_ϑV_L[μ](ϑ_i) = −⟨r∇_ϑΦ̃(·;ϑ_i)⟩_ρ for a frozen residual r."""
        return -self.learning(Th, r=r)

    def drift(self, Th, T, w=None):
        """dΘ/dτ = G⟨r∇_ϑΦ̃(ϑ_i)⟩_ρ − σ₀(ϑ_i − T_i)."""
        return self.G[None, :] * self.learning(Th, w) - self.sigma0 * (np.asarray(Th, float) - np.asarray(T, float))

    def residual_curvature(self, r, Th, h=1e-6):
        """R(ϑ_i) = ∇²_ϑ⟨rΦ̃(·;ϑ)⟩_ρ at ϑ_i with r frozen, (m, p, p), by central differences of the learning force."""
        Th = np.asarray(Th, float)
        R = np.zeros((Th.shape[0], self.p, self.p))
        for j in range(self.p):
            E = np.zeros_like(Th); E[:, j] = h
            R[:, :, j] = (self.learning(Th + E, r=r) - self.learning(Th - E, r=r)) / (2 * h)
        return 0.5 * (R + np.transpose(R, (0, 2, 1)))

    def kappa_R(self, r, Th):
        """λ_max(G^{1/2}R(ϑ_i)G^{1/2}) per row: Proposition D‴.5's splitting margin is σ₀ minus it."""
        sq = np.sqrt(self.G)
        R = self.residual_curvature(r, Th)
        return np.linalg.eigvalsh(sq[None, :, None] * R * sq[None, None, :])[:, -1]

    def whitened_jacobian(self, w, atoms):
        """Ĵ (N, K·p): the columns √w_k∇_ϑΦ̃(·;ϑ_k)G^{1/2}, k-major; H_L = ĴᵀĴ/N."""
        g = self.grad_features(atoms)
        return (np.sqrt(np.asarray(w, float))[None, :, None] * g * np.sqrt(self.G)[None, None, :]).reshape(self.N, -1)

    def gauss_newton(self, w, atoms):
        J = self.whitened_jacobian(w, atoms)
        H = J.T @ J / self.N
        ev, V = np.linalg.eigh(0.5 * (H + H.T))
        return H, ev, V


def weighted_critical_point(loop, w, Th0, prior=None, iters=80, tol=1e-13, h=1e-7):
    """Newton (finite-difference Jacobian, backtracking) on the weighted K-atom stationarity G⟨r_w∇Φ̃(ϑ_k)⟩ − σ₀(ϑ_k − ϑ⁰_k) = 0 (prior=None:
    the unregularised critical point). Returns (atoms, sup-norm residual)."""
    w = np.asarray(w, float)
    K, p = np.asarray(Th0).shape; d = K * p
    T0 = None if prior is None else np.asarray(prior, float)
    F = (lambda v: (loop.G[None, :] * loop.learning(v.reshape(K, p), w)).ravel()) if T0 is None else \
        (lambda v: loop.drift(v.reshape(K, p), T0, w).ravel())
    v = np.asarray(Th0, float).ravel().copy()
    for _ in range(iters):
        g = F(v); gn = np.abs(g).max()
        if gn < tol:
            break
        J = np.zeros((d, d))
        for i in range(d):
            e = np.zeros(d); e[i] = h
            J[:, i] = (F(v + e) - F(v - e)) / (2 * h)
        step = np.linalg.lstsq(J, g, rcond=None)[0]; lam = 1.0
        while lam > 1e-10:
            vn = v - lam * step
            if np.abs(F(vn)).max() < gn:
                break
            lam *= 0.5
        v = vn
    return v.reshape(K, p), float(np.abs(F(v)).max())


class SlidingRealLoop:
    """The real closed loop under the sliding controller: ẍ = f(x, ẋ) + u, u = ẍ_m − λė − k_s s − f̂_Θ(x, ẋ), s = ė + λe, and the adaptive law
    ϑ̇_i = G∇_ϑ(sΦ̃(z; ϑ_i)) − (σ₀/k_s)(ϑ_i − T_i) with leak targets T (frozen between refreshes)."""

    def __init__(self, net, f, x_ref, x_ref_dot, x_ref_ddot, ks, lam, G, sigma0):
        self.net = net
        self.f = f
        self.xr, self.xrd, self.xrdd = x_ref, x_ref_dot, x_ref_ddot
        self.ks, self.lam = float(ks), float(lam)
        self.G = np.asarray(G, float)
        self.sigma0 = float(sigma0)

    def f_hat_point(self, z, Th):
        a, c, s = self.net.physical(Th)
        phi = np.exp(-np.sum((np.asarray(z, float)[None, :] - c) ** 2, axis=1) / (2 * s ** 2))
        return float(np.mean(a * phi)), a, c, s, phi

    def rhs(self, t, x, v, Th, T):
        """(ẋ, v̇, Θ̇) and the error s."""
        e = x - self.xr(t); ed = v - self.xrd(t)
        sv = ed + self.lam * e
        z = np.array([x, v])
        fh, a, c, s, phi = self.f_hat_point(z, Th)
        u = self.xrdd(t) - self.lam * ed - self.ks * sv - fh
        vdot = float(self.f(z[None, :])[0]) + u
        da, dc, ds = self.net.chain(Th)
        d = self.net.d
        diff = z[None, :] - c
        sq = np.sum(diff ** 2, axis=1)
        g = np.empty_like(Th)
        g[:, 0] = sv * phi * da
        g[:, 1:1 + d] = (sv * a * phi / s ** 2)[:, None] * diff * dc
        g[:, 1 + d] = sv * a * phi * sq / s ** 3 * ds
        Thdot = self.G[None, :] * g - (self.sigma0 / self.ks) * (Th - T)
        return v, vdot, Thdot, sv


class SlidingServo:
    """The real closed loop with a regressor of one coordinate: ẍ = f(x) + u (f a function of the position, e.g. a cogging torque),
    u = ẍ_m − λė − k_s s − f̂_Θ(x), s = ė + λe, and the adaptive law ϑ̇_i = G∇_ϑ(sΦ̃(x; ϑ_i)) − (σ₀/k_s)(ϑ_i − T_i) with the one-dimensional
    saturated network of `langevin` (kernels on the position). RK4 with the parameters in the state; `targets(Θ)` gives the leak targets,
    refreshed every `refresh` steps. The averaged slow loop is Theorem J's along the position orbit x_m(t) (the module docstring)."""

    def __init__(self, net, f, x_ref, x_ref_dot, x_ref_ddot, ks, lam, G, sigma0):
        self.net = net
        self.f = f
        self.xr, self.xrd, self.xrdd = x_ref, x_ref_dot, x_ref_ddot
        self.ks, self.lam = float(ks), float(lam)
        self.G = np.asarray(G, float)
        self.sigma0 = float(sigma0)

    def field(self, t, x, v, Th, T):
        e = x - self.xr(t); ed = v - self.xrd(t)
        s = ed + self.lam * e
        fh = float(self.net.f_hat(np.array([x]), Th)[0])
        u = self.xrdd(t) - self.lam * ed - self.ks * s - fh
        vd = float(np.asarray(self.f(np.array([x])), float).ravel()[0]) + u
        return v, vd, self.G[None, :] * self.net.grad_g(np.array([x]), np.array([s]), Th) - (self.sigma0 / self.ks) * (Th - T)

    def run(self, Th0, targets, T_end, dt=2e-3, refresh=100, x0=None, v0=None, rec_every=100):
        """Integrate from the reference's state at t = 0 (or x0, v0); returns (Θ, times, e² at the records)."""
        x = float(self.xr(0.0)) if x0 is None else float(x0)
        v = float(self.xrd(0.0)) if v0 is None else float(v0)
        Th = np.array(Th0, float); t = 0.0
        T = targets(Th)
        ts, e2 = [], []
        for i in range(int(round(T_end / dt))):
            if i % refresh == 0:
                T = targets(Th)
            if i % rec_every == 0:
                ts.append(t); e2.append((x - float(self.xr(t))) ** 2)
            k1 = self.field(t, x, v, Th, T)
            k2 = self.field(t + dt / 2, x + dt / 2 * k1[0], v + dt / 2 * k1[1], Th + dt / 2 * k1[2], T)
            k3 = self.field(t + dt / 2, x + dt / 2 * k2[0], v + dt / 2 * k2[1], Th + dt / 2 * k2[2], T)
            k4 = self.field(t + dt, x + dt * k3[0], v + dt * k3[1], Th + dt * k3[2], T)
            x += dt / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]); v += dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
            Th = Th + dt / 6 * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]); t += dt
        return Th, np.array(ts), np.array(e2)


def filtered_error(r_orbit, period, ks, lam, exact=False):
    """The averaged prediction of the position error along one period: s = r(x_m)/k_s (the slaved value; exact=True: the periodic solution
    of ṡ = −k_s s + r, i.e. (d/dt + k_s)⁻¹r, which keeps the error filter's attenuation of r's harmonics) and e the periodic solution of
    ė = −λe + s (FFT on the uniform grid of the orbit). Returns (s, e) on the grid."""
    r_orbit = np.asarray(r_orbit, float)
    N = r_orbit.size
    k = np.fft.fftfreq(N, d=period / N) * 2 * np.pi
    s = np.real(np.fft.ifft(np.fft.fft(r_orbit) / (1j * k + ks))) if exact else r_orbit / ks
    e = np.real(np.fft.ifft(np.fft.fft(s) / (1j * k + lam)))
    return s, e
