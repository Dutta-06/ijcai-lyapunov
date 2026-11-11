"""
Proposition D‴.44: the replicator's own masses on a kernel dictionary under measurement noise — the flat and the paradigm's update, the
free masses' gauge equilibrium and logarithmic selection, the conserved mass's inflation, the null kernel under white noise
====================================================================================================================================

The loop of Proposition D‴.34 at fixed positions on a dictionary (Φ, N × d): masses β (Σβ = m, or free), amplitudes a = a_max tanh(α/a_max),
w = β⊙a/m, label noise on the measured error; the flat update's selection along the exact fits is the gradient flow of tr Ψ(K) in the flat
coordinates (α, √β), K = (m/N)Φ diag(D)Φᵀ (− (γ/N)yyᵀ with the conserved mass), D_i = (β_i/m²)(Gβ_ia′_i² + γa_i²); the paradigm's update
on log β adds the replicator drift of the filtered second moment of the log-masses' increment, and is not a gradient (`mass_dictionary`).

A. Identities: the closed forms against Katzenberger's drift of the scheme (`gauge_drift.limit_drift` on a chart of M) for the flat and the
   paradigm's updates, fresh and AR(1) draws at a finite step, conserved and free masses; the cancellation of the masses' own selection;
   the free law's gauge invariant Q: its drift against Katzenberger's, its quadratic variation; loop integrals of the drifts on M (the
   flat drift a gradient, the paradigm's not); K's spectrum against the flow's normal rates; the excess-mass bound with the conserved K;
   D* against its two laws and its split against the paradigm's gauge equilibrium; the null kernel's potential against the free one's
   (white: a constant; slow: not); along the white flows of Part B's instance, the Itô term against the flat drift, the conserved law's
   cross term against its net drift and their pushes on ‖w‖₁, the free masses' gauge ratio.
B. The paradigm's noisy scheme against the reduced flow, projected: from the flow's state at s = S0_B, every path's state at the window's
   ends mapped onto M by the noiseless flow (Katzenberger's limit map), the projected advance over the window against the flow's
   (ratio ± s.e.), the raw mean's offset from the projected one; fresh draws, and AR(1) draws (τ_c = 10³) at σ²τ_c from 0.05 to 0.9, the
   steps 1 and 0.5 — the reduced flow is the limit σ²τ_c/h → 0.
C. The selection's ends on random instances (N = 5, d = 12 Gaussian designs, 2- and 3-sparse truths of ‖w*‖₁ < a_max): the conserved, free
   and null-kernel laws under white and slow (AR(1), τ_c = 10³, 10⁴) noise — ‖w‖₁/a_max, the excess, the atoms used and alive, the truth
   recovered and the agreement with basis pursuit, split by whether basis pursuit is the truth; the null kernel against free masses per
   instance; the free masses' gauge ratio at the ends; the flat update's end and the minimiser of Σ_i|φ_i|²D*(|w_i|).
D. The same on Gaussian kernel dictionaries (16 centres, 8 samples on [0, 1]).
E. The conserved law's white end for a_max = 0.5, 1, 2, 4: ‖w‖₁/a_max, the excess against (iii)'s bound, the amplitudes' saturation.
F. Katzenberger's range (referee round 146): along Part B's instance's white flows (the paradigm's update) to s = 10⁸, K's spectrum (the
   normal rates), the masses, Q/c on the free law's dying and alive atoms, and λ_min(K) over the selection's fastest relative mass rate at
   Part B's noise (σ = SIG_B, h = 1) — the reduced flow is the scheme's limit only while that ratio is large; the paradigm's scheme
   against the reduced flow, projected as in B, from the free flow's state at s = 10⁵ and 2·10⁵ at two σ; Q's drift's sign on the free
   law's curves for both updates (κ = 1, 2) and its zeros against their closed forms; Q/c and the Itô/flat ratio at D*'s split and at the
   paradigm's gauge equilibrium; the null kernel's potential at a finite step (fresh draws, Ψ = −(2/h)log(1 − hκ/2)) at Part A's points.

Output: results/chaos/mass_dictionary.md and mass_dictionary_raw.npz (MD_WORKERS; `--from-raw` regenerates the md; `--only A,B,C,D,E,F`).
"""
import os
import sys
import time
import warnings
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
from src.core import mass_dictionary as MD
from src.core import colour_selection as CS
from src.core import gauge_drift as GDR
from src.core.runcache import cached
import exp_cohesion_leak as E

OUT = os.environ.get('MD_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
N_WORKERS = int(os.environ.get('MD_WORKERS', 2))
N_I, D_I, M_I, G_I, GAM_I, AMAX = 5, 12, 12.0, 1.0, 0.1, 1.0
MODES = ('conserved', 'free', 'null kernel')
TAUS = (0.0, 1e3, 1e4)
S_END = 1e8
FLOW_SECS = 900.0                                                      # a flow's wall-clock limit (one rbf instance's piece ran for hours)
KS_C, N_INST_C = (2, 3), 16
KS_D, N_INST_D, N_D, D_D, WIDTH_D = (2, 3), 8, 8, 16, 0.08
AMAXES_E, N_INST_E = (0.5, 1.0, 2.0, 4.0), 12
TOL_BP = 1e-3                                                          # an end "at basis pursuit": within this (largest coefficient's distance)
# Part B: the window [S0_B + burn-in, … + DS_B] of the reduced flow (s units) — the late stretch where the conserved law's net drift is a
# small difference of large terms (referee round 34) —, the burn-in BURN_B in time units (10τ_c, about eight of the slowest normal
# relaxation times), P_B paths; the runs (τ_c, σ²τ_c, h), fresh draws (τ_c = 0) at σ = SIG_B
SEED_B, S0_B, DS_B, BURN_B, P_B, SIG_B = 4, 1500.0, 450.0, 1e4, 64, 0.1
B_RUNS = ((0.0, None, 1.0), (1e3, 0.9, 1.0), (1e3, 0.225, 1.0), (1e3, 0.05, 1.0), (1e3, 0.45, 0.5), (1e3, 0.1, 0.5))
LOOP_EPS = (1e-2, 5e-3, 2.5e-3)
S_ALONG = (0.0, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7)
TAU_NULL = 1e4
# Part F (referee round 146): the flows' sample points S_F (the flows to S_END, as C's), Part B's step H_F; the masses' rates over the atoms
# with β > LIVE_F; the scheme's runs from the free flow's state at S0_F, at the noise levels SIGS_F (fresh draws), windows of DS_F in s after
# a burn-in of BURN_RELAX_F of the slowest normal relaxation times (1/λ_min(K)), the limit map run for SETTLE_RELAX_F of them (at least
# settle's 2·10⁵); Q's drift on the curves Q/c ∈ QC_F at the amplitudes TS_F; the null kernel's finite steps H_NULL_F; RATIO_F the target
# ratio for the σ column
S_F = (0.0, 1e3, 1500.0, 3e3, 1e4, 3e4, 1e5, 2e5, 3e5, 5e5, 1e6, 3e6, 1e7, 3e7, 1e8)
H_F, LIVE_F, RATIO_F = 1.0, 1e-12, 10.0
S0_F, SIGS_F, DS_F, BURN_RELAX_F, SETTLE_RELAX_F = (1e5, 2e5), (0.1, 0.05), 2000.0, 10.0, 100.0
QC_F = (1.0, 0.0, -0.5, -0.9, -1.0, -1.1, -2.0, -10.0)
TS_F = np.linspace(0.005, 0.995, 199)
H_NULL_F = (1.0, 10.0)


def instance(k, i, kind='gauss'):
    """A design with unit columns and a k-sparse truth (signs ±1, magnitudes in [0.15, 0.3], so ‖w*‖₁ < a_max = 1) — Gaussian (N_I × D_I)
    or Gaussian kernels on [0, 1] (N_D samples, D_D centres) — and the start's random amplitudes."""
    rng = np.random.default_rng((1000 if kind == 'gauss' else 5000) + 100 * k + i)
    if kind == 'gauss':
        Phi = rng.standard_normal((N_I, D_I)); Phi /= np.linalg.norm(Phi, axis=0)
    else:
        Phi = CS.rbf_dictionary(np.linspace(0.0, 1.0, N_D), np.linspace(0.0, 1.0, D_D), WIDTH_D)
    d = Phi.shape[1]
    ws = np.zeros(d); idx = rng.choice(d, k, replace=False); ws[idx] = rng.choice([-1, 1], k) * rng.uniform(0.15, 0.3, k)
    return Phi, ws, 0.05 * rng.standard_normal(d + 1)


def law(Phi, ws, mode, a_max=AMAX):
    """The MassDictionary for a mode (the null kernel: a zero column appended) and its design."""
    P2 = np.c_[Phi, np.zeros(Phi.shape[0])] if mode == 'null kernel' else Phi
    return MD.MassDictionary(P2, Phi @ ws, M_I, a_max=a_max, G=G_I, gamma=GAM_I, balanced=(mode != 'free')), P2


def start(X, al0):
    """The noiseless flow's end from the amplitudes al0 and equal masses m/(atoms) (the loop's own law, in (α, log β)), on M."""
    return X.settle(al0[:X.d], np.full(X.d, X.m / X.d))


def filt_of(tau):
    return (lambda k: np.ones_like(k)) if tau == 0 else (lambda k, t=tau: 1.0 / (1.0 + t * k))


def summary(X, th, ws, wbp, d0):
    """The end's ‖w‖₁/a_max, excess/m over the real atoms (1 − ‖w‖₁/a_max: the null atom's coefficient is not part of the fit), atoms used
    (|w| > 10⁻³), alive (β > 10⁻³), the truth recovered (10⁻⁴), distances, the null mass, the largest amplitude, the coefficients, and
    for free masses the gauge ratio Gβ(1 − t⁴)/(γa²) on the alive atoms."""
    w = X.w(th)[:d0]; b = X.parts(th)[2]
    l1 = float(np.abs(w).sum() / X.a_max)
    alive = b[:d0] > 1e-3
    return dict(l1=l1, exc=1.0 - l1 if X.balanced else np.nan, used=int((np.abs(w) > 1e-3).sum()), alive=int(alive.sum()),
                ok=bool(np.abs(w - ws).max() < 1e-4 * np.abs(ws).max()), d_truth=float(np.abs(w - ws).max()),
                d_bp=float(np.abs(w - wbp).max()), b0=float(b[d0]) if len(b) > d0 else np.nan,
                amax=float(np.abs(X.parts(th)[3][:d0]).max() / X.a_max), w=w.copy(),
                gr=(X.gauge_ratio(th)[:d0][alive].copy() if not X.balanced else None))


def dstar_minimiser(X, d0):
    """The minimiser of Σ_i|φ_i|²D*(|w_i|) on {Φw = y} (free masses' white selection after the gauge, the flat update), by SLSQP with D*'s
    envelope gradient; D* from the inner minimisation on a log grid, interpolated."""
    from scipy.optimize import minimize, minimize_scalar
    G, g, am, m = X.G, X.gamma, X.a_max, X.m
    xs = np.geomspace(1e-9, 2.0, 400)
    ts = np.array([minimize_scalar(lambda t: G * x * x * (1 - t * t) ** 2 / (am * am * t * t) + g * x * am * t / m, bounds=(1e-9, 1.0),
                                   method='bounded', options=dict(xatol=1e-13)).x for x in xs])
    tstar = lambda x: np.interp(np.log(np.maximum(x, 1e-9)), np.log(xs), ts)
    Phi, y = X.Phi[:, :d0], X.y
    c = (Phi ** 2).sum(axis=0)

    def f(w):
        x = np.abs(w); t = tstar(x)
        val = G * x * x * (1 - t * t) ** 2 / (am * am * t * t) + g * x * am * t / m
        dv = (2 * G * x * (1 - t * t) ** 2 / (am * am * t * t) + g * am * t / m) * np.sign(w)
        return float(c @ val), c * dv
    w0 = np.linalg.lstsq(Phi, y, rcond=None)[0]
    res = minimize(f, w0, jac=True, method='SLSQP', constraints=[dict(type='eq', fun=lambda w: Phi @ w - y, jac=lambda w: Phi)],
                   options=dict(ftol=1e-14, maxiter=2000))
    return res.x


def flat_field(X, x):
    """The noiseless flow in the flat coordinates (α, ϱ) (the replicator's λ for the conserved mass)."""
    d = X.d
    _, rho, b, a, da, _ = X.parts(x)
    g = X.Phi.T @ (X.y - X.Phi @ (b * a / X.m)) / X.N
    u = a * g
    lam = b @ u / b.sum() if X.balanced else 0.0
    return np.r_[X.G * b * da * g, 0.5 * rho * X.gamma * (u - lam)]


def own_removed(X, th, filt):
    """The flat drift without the masses' own selection (the γ-term of D_i differentiated in √β): −P_g(P∇R) with that term dropped."""
    al, rho, b, a, da, dda = X.parts(th)
    om = (X.m / X.N) * X.weights(th, filt); m2 = X.m ** 2
    gr = np.r_[om * (b / m2) * (2 * X.G * b * da * dda + 2 * X.gamma * a * da), om * (2 * X.G * b * da ** 2 / m2) * 2 * rho]
    return -X.project(th, X.Pm * gr)


def mnorm(X, q):
    """A velocity's length in the flow's metric (|dα|²/G + (4/γ)|dϱ|²)^{1/2}."""
    return float(np.sqrt(q @ (q / X.Pm)))


def l1_rate(X, th, q):
    """d‖w‖₁/ds along the velocity q."""
    _, rho, b, a, da, _ = X.parts(th)
    sg = np.sign(X.w(th))
    return float(sg @ (b * da / X.m * q[:X.d] + 2 * rho * a / X.m * q[X.d:]))


def loop_integral(X, th0, u, v, eps, ito, n=64):
    """∮ ⟨P⁻¹q, dθ⟩ around θ(φ) = retract(θ0 + ε(cos φ u + sin φ v)) (u, v orthonormal tangent vectors in the metric), q the white
    selection's drift — zero for a gradient flow in that metric, ε²·(a 2-form) otherwise — by the trapezoid rule with spectral tangents."""
    ph = 2 * np.pi * np.arange(n) / n
    P_ = np.array([X.retract(th0 + eps * (np.cos(p) * u + np.sin(p) * v)) for p in ph])
    k = np.fft.fftfreq(n, 1.0 / n)
    dP = np.real(np.fft.ifft(1j * k[:, None] * np.fft.fft(P_, axis=0), axis=0))
    one = lambda kk: np.ones_like(kk)
    return float(sum((X.drift(p, one, ito=ito) / X.Pm) @ dp for p, dp in zip(P_, dP)) * 2 * np.pi / n)


def part_A():
    """The identities at random exact fits (a 5 × 12 Gaussian design) and along Part B's white flows."""
    warnings.filterwarnings('ignore')
    rng = np.random.default_rng(3)
    Phi = rng.standard_normal((N_I, D_I)); Phi /= np.linalg.norm(Phi, axis=0)
    ws = np.zeros(D_I); ws[[2, 7]] = [0.2, -0.17]
    rows, canc, spec, bound, qd, loops = [], [], [], [], [], []
    for mode in ('conserved', 'free'):
        X, _ = law(Phi, ws, mode)
        for rep in range(3):
            beta = rng.uniform(0.5, 1.5, D_I); beta *= M_I / beta.sum()
            th = X.retract(X.theta(0.3 * rng.standard_normal(D_I), beta))
            B = np.argsort(-np.abs(X.parts(th)[3]))[:N_I]
            for h, tc in ((0.0, 0.0), (0.5, 200.0)):
                corr = ('fresh',) if tc == 0 else ('ar1', float(np.exp(-h / tc)))
                filt = (lambda k: np.ones_like(k)) if tc == 0 else (lambda k, h_=h, c_=corr: MD.LN.lr_factor(k, h_, c_))
                for coords in ('flat', 'log'):
                    g, Sf = X.scheme_drift(th, B, h=h, corr=corr, coords=coords)
                    q = X.drift(th, filt, ito=(coords == 'log')) * X.m / 4
                    rows.append(dict(mode=mode, rep=rep, tc=tc, coords=coords, rel=float(np.abs(g - q).max() / np.abs(g).max())))
                    if mode == 'free':                                  # Q's drift (per unit v) and its quadratic variation
                        al, rho, b, a, da, _ = X.parts(th)
                        c = X.gamma * X.a_max ** 2 / (4 * X.G)
                        dQ = -c * (2 / X.a_max) * np.sinh(2 * al / X.a_max) * g[:X.d] + 2 * rho * g[X.d:]
                        cf = X.q_drift(th, filt, ito=(coords == 'log')) * X.m / 4
                        rec = dict(mode=mode, rep=rep, tc=tc, coords=coords, rel=float(np.abs(dQ - cf).max() / np.abs(cf).max()))
                        if coords == 'log':
                            gz = np.r_[-c * (2 / X.a_max) * np.sinh(2 * al / X.a_max), b]
                            blk = lambda i: [i, X.d + i]
                            qv = [abs(gz[blk(i)] @ Sf[np.ix_(blk(i), blk(i))] @ gz[blk(i)]) / (gz[blk(i)] @ gz[blk(i)] * np.abs(Sf[np.ix_(blk(i), blk(i))]).max())
                                  for i in range(X.d)]
                            rec['qv'] = float(max(qv))
                        qd.append(rec)
                if coords == 'log':
                    full = X.drift(th, filt, ito=True)
                    canc.append(dict(mode=mode, rep=rep, tc=tc, rel=float(np.abs(full - own_removed(X, th, filt)).max() / np.abs(full).max())))
            A = GDR._jac(lambda x: flat_field(X, x), th, 1e-6)
            ev = np.sort(-np.linalg.eigvals(A).real)[::-1][:N_I]
            kap = np.sort(np.linalg.eigvalsh(X.gram(th)))[::-1]
            spec.append(dict(mode=mode, rep=rep, rel=float(np.abs(ev - kap).max() / kap.max())))
            if X.balanced:
                E_ = X.excess(th)
                bound.append(dict(rep=rep, E=float(E_ / X.m), trK=float(X.potential(th, lambda k: k)),
                                  trKf=float(X.potential(th, lambda k: k) + X.gamma / X.N * X.y @ X.y),
                                  lb=float(X.white_bound(E_)), sba=float(X.parts(th)[2] @ X.parts(th)[4] / X.m)))
            if rep == 0:                                                # loop integrals of the white drifts on M
                q1, q2 = X.project(th, rng.standard_normal(2 * D_I)), X.project(th, rng.standard_normal(2 * D_I))
                ip = lambda x, y: x @ (y / X.Pm)
                u = q1 / np.sqrt(ip(q1, q1)); q2 = q2 - ip(q2, u) * u; v = q2 / np.sqrt(ip(q2, q2))
                for ito in (False, True):
                    loops.append(dict(mode=mode, ito=ito, I=[loop_integral(X, th, u, v, e, ito) for e in LOOP_EPS]))
    X, _ = law(Phi, ws, 'free')
    thr = X.gamma * X.a_max ** 3 / (X.G * X.m)
    xs = np.array([1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3]) * 1.0
    ds = X.relevance_min(xs)
    t_flat, t_par = X.gauge_t(xs, paradigm=False), X.gauge_t(xs, paradigm=True)
    am, m = X.a_max, X.m
    d_at = X.G * xs ** 2 * (1 - t_flat ** 2) ** 2 / (am * am * t_flat ** 2) + X.gamma * xs * am * t_flat / m       # D on the orbit at t_flat
    # the null kernel's potential against the free masses' at six points of M_null (white and slow Ψ)
    Y, _ = law(Phi, ws, 'null kernel')
    Xf, _ = law(Phi, ws, 'free')
    nul = []
    for rep in range(6):
        beta = rng.uniform(0.5, 1.5, D_I + 1); beta *= M_I / beta.sum()
        th = Y.retract(Y.theta(0.3 * rng.standard_normal(D_I + 1), beta))
        thr_ = np.r_[th[:D_I], th[D_I + 1:2 * D_I + 1]]
        nul.append(dict(white=Y.potential(th, lambda k: k) - Xf.potential(thr_, lambda k: k),
                        slow=Y.potential(th, lambda k: np.log1p(TAU_NULL * k) / TAU_NULL) - Xf.potential(thr_, lambda k: np.log1p(TAU_NULL * k) / TAU_NULL),
                        on_M=float(np.abs(Xf.constraints(thr_)).max())))
    return dict(rows=rows, canc=canc, spec=spec, bound=bound, qd=qd, loops=loops, thr=thr, xs=xs, ds=ds, t_flat=t_flat, t_par=t_par,
                d_at=d_at, lin=X.gamma * X.a_max * xs / X.m, four3=1.5 * (2 * X.G) ** (1 / 3) * (X.gamma / X.m) ** (2 / 3) * xs ** (4 / 3),
                null=nul, null_const=-X.gamma / X.N * float(X.y @ X.y), along=along_flows())


def along_flows():
    """Along Part B's instance's white flows (the paradigm's update): the Itô term against the flat drift (lengths in the metric), the
    conserved law's cross term against its net drift and the two pushes on ‖w‖₁ (the cross term's and the rest's), the free masses'
    gauge ratio on the alive atoms (β > 10⁻³)."""
    Phi, ws, al0 = instance(2, SEED_B)
    one = filt_of(0.0)
    out = []
    for mode in ('conserved', 'free'):
        X, _ = law(Phi, ws, mode)
        _, Th = X.flow(start(X, al0), S_ALONG[-1], one, ito=True, t_eval=np.array(S_ALONG), chunks=60)
        for s, th in zip(S_ALONG, Th):
            flat = X.drift(th, one, ito=False)
            ito = X.project(th, X.ito_term(th, filt=one) * (4.0 / X.m))
            net = flat + ito
            r = dict(mode=mode, s=s, l1=float(np.abs(X.w(th)).sum() / X.a_max), ito=mnorm(X, ito) / mnorm(X, flat))
            if X.balanced:
                cr = X.project(th, X.cross_term(th, one) * (4.0 / X.m))
                r.update(cross=mnorm(X, cr) / mnorm(X, net), push_cross=l1_rate(X, th, cr), push_rest=l1_rate(X, th, net - cr))
            else:
                b = X.parts(th)[2]
                gr = X.gauge_ratio(th)[b > 1e-3]
                r.update(gr_lo=float(gr.min()), gr_hi=float(gr.max()), alive=int((b > 1e-3).sum()))
            out.append(r)
    return out


def scheme_window(mode, tc, s2tc, h, s0, sig_fresh, ds_win, burn=None, settle_T=2e5):
    """The paradigm's scheme from the reduced flow's state at s = s0 (Part B's instance), a burn-in of `burn` time units (None: BURN_RELAX_F
    slowest normal relaxation times at s0), then a window of ds_win (s units); every path's state at the window's ends mapped onto M by the
    noiseless flow run for settle_T (`settle`; None: SETTLE_RELAX_F slowest relaxation times, at least 2·10⁵); the projected advance against
    the flow's over the same window (the flow at the step's filter). Fresh draws (τ_c = 0) at σ = sig_fresh, AR(1) at σ²τ_c = s2tc."""
    Phi, ws, al0 = instance(2, SEED_B)
    X, _ = law(Phi, ws, mode)
    th0 = start(X, al0)
    rho = 0.0 if tc == 0 else float(np.exp(-h / tc))
    corr = ('fresh',) if tc == 0 else ('ar1', rho)
    filt = lambda k: MD.LN.lr_factor(k, h, corr)
    sig = sig_fresh if tc == 0 else float(np.sqrt(s2tc / tc))
    ds = X.m * (h * sig ** 2 / X.N) / 4 * h                                    # s per step
    adaptive = burn is None or settle_T is None
    if adaptive:                                                               # the flow to s0, then on from its state there
        th_s0 = X.flow(th0, s0, filt, ito=True, t_eval=np.array([s0]), chunks=8)[1][-1]
        kap0 = np.sort(np.linalg.eigvalsh(X.gram(th_s0)))[::-1]
        lam0 = float(kap0[-1])
        burn = BURN_RELAX_F / lam0 if burn is None else burn
        settle_T = max(2e5, SETTLE_RELAX_F / lam0) if settle_T is None else settle_T
    n_burn, n_win = int(round(burn / h)), int(round(ds_win / ds))
    s_a = s0 + n_burn * ds; s_b = s_a + n_win * ds
    if adaptive:
        _, Th = X.flow(th_s0, s_b - s0, filt, ito=True, t_eval=np.array([0.0, s_a - s0, s_b - s0]), chunks=8)
    else:
        _, Th = X.flow(th0, s_b, filt, ito=True, t_eval=np.array([s0, s_a, s_b]), chunks=8)
    wa, wb = X.w(Th[1]), X.w(Th[2])
    b0 = X.parts(Th[0])[2]
    rec, _, _, kept = X.noisy(Th[0][:X.d], b0, h, sig, rho, n_burn + n_win, np.random.default_rng(91), P_B, n_burn + n_win,
                              keep=(n_burn, n_burn + n_win))
    proj = [np.array([X.w(X.settle(al[j], b[j], T=settle_T)) for j in range(P_B)]) for al, b in kept]
    raw = [b * (X.a_max * np.tanh(al / X.a_max)) / X.m for al, b in kept]
    df = wb - wa
    ratio = (proj[1] - proj[0]) @ df / (df @ df)
    off = raw[1] - proj[1]
    out = dict(mode=mode, tc=tc, s2tc=s2tc, h=h, sig=sig, n_burn=n_burn, n_win=n_win, s_a=s_a, s_b=s_b, df=float(np.abs(df).max()),
               ratio=float(ratio.mean()), ratio_se=float(ratio.std(ddof=1) / np.sqrt(P_B)),
               start_off=float(np.abs(proj[0].mean(axis=0) - wa).max()),
               off=float(np.abs(off.mean(axis=0)).max()), off_z=float(np.abs(off.mean(axis=0) / (off.std(axis=0, ddof=1) / np.sqrt(P_B) + 1e-300)).max()),
               dev_raw=float(np.abs(raw[1].mean(axis=0) - wb).max()), dev_proj=float(np.abs(proj[1].mean(axis=0) - wb).max()))
    if adaptive:
        out.update(s0=s0, lam0=lam0, settle_T=settle_T, kap0=kap0, rate0=X.drift(th_s0, filt, ito=True)[X.d:] * 2 / np.sqrt(X.parts(th_s0)[2]),
                   beta0=X.parts(th_s0)[2])
    return out


@cached
def b_job(args):
    """Part B, one (mode, τ_c, σ²τ_c, h): `scheme_window` from s = S0_B, a burn-in of BURN_B time units, a window of DS_B."""
    warnings.filterwarnings('ignore')
    mode, tc, s2tc, h = args
    t0 = time.time()
    out = scheme_window(mode, tc, s2tc, h, S0_B, SIG_B, DS_B, burn=BURN_B)
    print(f"[{time.time() - t0:.0f}s] B {mode} tau_c {tc:g} s2tc {s2tc} h {h:g}", flush=True)
    return out


def run_flow(X, th0, filt, ito, ws, wbp, d0):
    """One selection flow to S_END → its end's summary, or, if it stopped (over the wall-clock limit or failed), the stop and the last
    state's summary."""
    try:
        th = X.flow(th0, S_END, filt, ito=ito, chunks=60, max_secs=FLOW_SECS)[1][-1]
        return summary(X, th, ws, wbp, d0), th
    except MD.FlowStopped as e:
        return dict(stopped='over the limit' if isinstance(e, MD.FlowTimeout) else 'failed', s=e.s,
                    last=summary(X, e.th, ws, wbp, d0) if e.th is not None else None), None


@cached
def c_job(args):
    """Parts C and D, one (kind, k, instance): the ends of the reduced flows (the paradigm's update) for every mode and colour; for white
    noise also the flat update's end; basis pursuit and the minimiser of Σ|φ_i|²D*(|w_i|)."""
    warnings.filterwarnings('ignore')
    kind, k, i = args
    t0 = time.time()
    Phi, ws, al0 = instance(k, i, kind)
    d0 = Phi.shape[1]
    wbp = CS.basis_pursuit(Phi, Phi @ ws, (Phi ** 2).sum(axis=0))
    out = dict(kind=kind, k=k, i=i, bp_ok=bool(np.abs(wbp - ws).max() < 1e-4 * np.abs(ws).max()), l1_truth=float(np.abs(ws).sum()),
               l1_bp=float(np.abs(wbp).sum()))
    for mode in MODES:
        X, P2 = law(Phi, ws, mode)
        th0 = start(X, al0)
        out[(mode, 'start')] = summary(X, th0, ws, wbp, d0)
        for tau in TAUS:
            out[(mode, tau)] = run_flow(X, th0, filt_of(tau), True, ws, wbp, d0)[0]
        if mode in ('conserved', 'free'):
            out[(mode, 'flat')], th = run_flow(X, th0, filt_of(0.0), False, ws, wbp, d0)
            if mode == 'free':
                wd = dstar_minimiser(X, d0)
                out['dstar_vs_flat'] = float(np.abs(X.w(th)[:d0] - wd).max()) if th is not None else np.nan
                out['dstar_vs_bp'] = float(np.abs(wd - wbp).max())
    print(f"[{time.time() - t0:.0f}s] C {kind} k {k} instance {i}", flush=True)
    return out


@cached
def e_job(args):
    """Part E, one (instance, a_max): the conserved law's white end against the representability boundary, for several a_max (the
    instance kept when its basis-pursuit solution is representable, ‖w_bp‖₁ < a_max)."""
    warnings.filterwarnings('ignore')
    i, am = args
    Phi, ws, al0 = instance(2, i)
    wbp = CS.basis_pursuit(Phi, Phi @ ws, (Phi ** 2).sum(axis=0))
    if np.abs(wbp).sum() >= am:
        return dict(i=i, am=am, skip=True)
    X, _ = law(Phi, ws, 'conserved', a_max=am)
    th0 = start(X, al0 * am)
    end, th = run_flow(X, th0, filt_of(0.0), True, ws, wbp, Phi.shape[1])
    if th is None:
        return dict(i=i, am=am, skip=True, stopped=end)
    c2 = (Phi ** 2).sum(axis=0)
    return dict(i=i, am=am, skip=False, l1_bp=float(np.abs(wbp).sum()), start=summary(X, th0, ws, wbp, Phi.shape[1]), end=end,
                bound=float(am * np.sqrt(X.gamma * X.m * np.sum(1.0 / c2) * c2.max() / X.G) / X.m))      # (iii)'s bound on E/m


def f_flows():
    """Part F: along Part B's instance's white flows (the paradigm's update; the flows to S_END with C's pieces, so the free end is C's), at
    s ∈ S_F: K's spectrum (the normal rates, decreasing), the masses, the coefficients, Q/c (free masses) and the masses' relative rates
    d log β_i/ds."""
    Phi, ws, al0 = instance(2, SEED_B)
    one = filt_of(0.0)
    out = []
    for mode in ('free', 'conserved'):
        X, _ = law(Phi, ws, mode)
        _, Th = X.flow(start(X, al0), S_F[-1], one, ito=True, t_eval=np.array(S_F), chunks=60)
        c = X.gamma * X.a_max ** 2 / (4 * X.G)
        for s, th in zip(S_F, Th):
            b = X.parts(th)[2]
            out.append(dict(mode=mode, s=float(s), kap=np.sort(np.linalg.eigvalsh(X.gram(th)))[::-1], beta=b.copy(), w=X.w(th),
                            qc=None if X.balanced else X.gauge_invariant(th) / c, rate=X.drift(th, one, ito=True)[X.d:] * 2 / np.sqrt(b)))
    return out


def ds_dt(sig, h=H_F):
    """s per unit time, s = (m v_h/4)t, v_h = hσ²/N (Part A's instance's m and N)."""
    return M_I * (h * sig ** 2 / N_I) / 4


def range_ratio(kap, beta, rate, sig=SIG_B, h=H_F):
    """λ_min(K) (the slowest normal rate, per unit time) over the selection's fastest relative mass rate per unit time at the noise σ and
    step h — max |d log β_i/ds| over the atoms with β_i > LIVE_F, times ds/dt. Katzenberger's limit needs it large."""
    live = np.asarray(beta) > LIVE_F
    return float(np.min(kap) / (np.abs(np.asarray(rate)[live]).max() * ds_dt(sig, h)))


def crossing(s, r, level=1.0):
    """The first s at which r falls below `level` (log–log interpolation between the samples), or None."""
    s, r = np.asarray(s, float), np.asarray(r, float)
    for j in range(1, len(s)):
        if r[j - 1] >= level > r[j] and s[j - 1] > 0:
            x0, x1, y0, y1 = np.log(s[j - 1]), np.log(s[j]), np.log(r[j - 1]), np.log(r[j])
            return float(np.exp(x0 + (np.log(level) - y0) * (x1 - x0) / (y1 - y0)))
    return None


def curve_state(X, qc, t):
    """θ on the free law's curve Q = qc·c at the amplitudes t (t = a/a_max; β = Q + c cosh(2α/a_max), cosh(2α/a_max) = (1 + t²)/(1 − t²))."""
    c = X.gamma * X.a_max ** 2 / (4 * X.G)
    t = np.asarray(t, float)
    return X.theta(X.a_max * np.arctanh(t), c * (qc + (1 + t * t) / (1 - t * t)))


def q_zero_closed(qc, kappa):
    """The closed-form zero t > 0 of Q's drift on the curve Q/c = qc: κ = 2 at t = 1/κ′ = ((1 + qc)/(1 − qc))^{1/2} (−1 < qc < 0), κ = 1 at
    t = ((−1 − qc)/(1 − qc))^{1/4} (qc < −1, on a branch); None elsewhere."""
    if kappa == 2 and -1 < qc < 0:
        return float(np.sqrt((1 + qc) / (1 - qc)))
    if kappa == 1 and qc < -1:
        return float(((-1 - qc) / (1 - qc)) ** 0.25)
    return None


def f_signs():
    """Part F: Q's drift (`q_drift`, white, κ = 1 flat, 2 the paradigm's) on the free law's curves Q/c ∈ QC_F at the amplitudes TS_F within
    the curve's range (β > 0), on states built on each curve (a one-row design of ones: ω = 1); its zeros by bisection on the library's
    drift, against `q_zero_closed`."""
    from scipy.optimize import brentq
    X1 = MD.MassDictionary(np.ones((1, 1)), np.zeros(1), M_I, a_max=AMAX, G=G_I, gamma=GAM_I, balanced=False)
    one = filt_of(0.0)
    out = []
    for qc in QC_F:
        t = TS_F[qc + (1 + TS_F ** 2) / (1 - TS_F ** 2) > 0]                       # the curve's range: β > 0
        Xt = MD.MassDictionary(np.ones((1, t.size)), np.zeros(1), M_I, a_max=AMAX, G=G_I, gamma=GAM_I, balanced=False)
        th = curve_state(Xt, qc, t)
        rec = dict(qc=qc, t=t)
        for kappa in (1, 2):
            v = Xt.q_drift(th, one, ito=(kappa == 2))
            g = lambda x: float(X1.q_drift(curve_state(X1, qc, [x]), one, ito=(kappa == 2))[0])
            ch = np.where(np.sign(v[:-1]) != np.sign(v[1:]))[0]
            zeros = [brentq(g, t[j], t[j + 1], xtol=1e-14) for j in ch]
            rec[kappa] = dict(sign=np.sign(v), zeros=zeros, closed=q_zero_closed(qc, kappa))
        out.append(rec)
    return out


def f_splits(xs):
    """Part F: at |w| ∈ xs, D*'s split (the flat update's rest) and the paradigm's gauge equilibrium (`gauge_t`): Q/c on states built
    there (`gauge_invariant`) against −(1 + t⁴)/(1 − t⁴) and −(1 − t²)/(1 + t²); at the paradigm's equilibrium the Itô term's mass
    component over the flat drift's before the projection (`ito_term`, `grad`) against (1 + t²)/(3 − t²)."""
    xs = np.asarray(xs, float)
    X = MD.MassDictionary(np.ones((1, xs.size)), np.zeros(1), M_I, a_max=AMAX, G=G_I, gamma=GAM_I, balanced=False)
    c = X.gamma * X.a_max ** 2 / (4 * X.G)
    one = filt_of(0.0)
    out = dict(xs=xs)
    for key, par in (('flat', False), ('par', True)):
        t = X.gauge_t(xs, paradigm=par)
        th = X.theta(X.a_max * np.arctanh(t), X.m * xs / (X.a_max * t))
        out[key] = dict(t=t, qc=X.gauge_invariant(th) / c, w=X.w(th))
        if par:
            ito = X.ito_term(th, filt=one)[X.d:] * (4.0 / X.m)
            flat = -(X.Pm * X.grad(th, one))[X.d:]
            out[key]['ito_flat'] = ito / -flat
    return out


def a_null_points():
    """Part A's six points of the null kernel's manifold, replayed: Part A's generator's draws up to them (the design; per law three
    points, the first with its loop's two tangent vectors), then the six points (the null law and the free law on the real atoms)."""
    rng = np.random.default_rng(3)
    Phi = rng.standard_normal((N_I, D_I)); Phi /= np.linalg.norm(Phi, axis=0)
    ws = np.zeros(D_I); ws[[2, 7]] = [0.2, -0.17]
    for _ in ('conserved', 'free'):
        for rep in range(3):
            rng.uniform(0.5, 1.5, D_I); rng.standard_normal(D_I)
            if rep == 0:
                rng.standard_normal(2 * D_I); rng.standard_normal(2 * D_I)
    Y, _ = law(Phi, ws, 'null kernel')
    Xf, _ = law(Phi, ws, 'free')
    pts = []
    for rep in range(6):
        beta = rng.uniform(0.5, 1.5, D_I + 1); beta *= M_I / beta.sum()
        th = Y.retract(Y.theta(0.3 * rng.standard_normal(D_I + 1), beta))
        pts.append((th, np.r_[th[:D_I], th[D_I + 1:2 * D_I + 1]]))
    return Y, Xf, pts


def psi_step(h):
    """Ψ for fresh draws at the step h: Ψ′ = 2/(2 − hκ) (`label_noise.lr_factor`), Ψ = −(2/h)log(1 − hκ/2)."""
    return lambda k: -(2.0 / h) * np.log1p(-h * k / 2.0)


def f_null():
    """Part F: R_null − R_free at Part A's six points under white noise, slow noise (τ = TAU_NULL) and fresh draws at the steps H_NULL_F."""
    Y, Xf, pts = a_null_points()
    slow = lambda k: np.log1p(TAU_NULL * k) / TAU_NULL
    out = []
    for th, thr in pts:
        r = dict(white=Y.potential(th, lambda k: k) - Xf.potential(thr, lambda k: k), slow=Y.potential(th, slow) - Xf.potential(thr, slow))
        for h in H_NULL_F:
            r[h] = Y.potential(th, psi_step(h)) - Xf.potential(thr, psi_step(h))
        out.append(r)
    return out


@cached
def f_job(args):
    """Part F, one (σ, s0): `scheme_window` for the free law, fresh draws at σ and h = H_F, from the reduced flow's state at s0."""
    warnings.filterwarnings('ignore')
    sig, s0 = args
    t0 = time.time()
    out = scheme_window('free', 0.0, None, H_F, s0, sig, DS_F, burn=None, settle_T=None)
    print(f"[{time.time() - t0:.0f}s] F sigma {sig:g} s0 {s0:g}", flush=True)
    return out


def part_F():
    """Part F's deterministic pieces (the scheme's runs, `f_job`, are added by `main`)."""
    warnings.filterwarnings('ignore')
    return dict(along=f_flows(), signs=f_signs(), splits=f_splits([1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3]), null=f_null())


def _fmt_range(x):
    x = np.asarray(x, float)
    return f"{np.median(x):.4f} [{x.min():.4f}, {x.max():.4f}]" if x.size else "—"


def sign_runs(t, sg):
    """A sign pattern along t as runs: '− on [t0, t1], + on [t2, t3]' ('0' for exact zeros)."""
    sym = {-1.0: '−', 0.0: '0', 1.0: '+'}
    runs, j = [], 0
    for k in range(1, len(sg) + 1):
        if k == len(sg) or sg[k] != sg[j]:
            runs.append(f"{sym[float(sg[j])]} on [{t[j]:.3f}, {t[k - 1]:.3f}]")
            j = k
    return ', '.join(runs)


def curve_name(qc):
    return ('a free curve, Q ≥ 0' if qc >= 0 else 'a free curve, −c < Q < 0' if qc > -1 else 'the critical curve' if qc == -1
            else 'a branch, Q < −c')


def report_F(F, A=None, C=None):
    """Part F's section of the md."""
    L_ = ["## F. Katzenberger's range along the flows, Q's sign maps, the null kernel at a finite step (referee round 146)", "",
          f"Along Part B's instance's white flows (the paradigm's update; to s = {S_F[-1]:.0e} with Part C's pieces): K's spectrum (the normal "
          f"rates, decreasing), the atoms alive (β > 10⁻³), the dying atoms' masses (β ≤ 10⁻³: smallest, largest), Q/c on the dying and on "
          "the alive atoms (free masses; c = γa_max²/(4G)), the fastest relative mass rate max|d log β_i/ds| over the atoms with β_i > "
          f"10⁻¹², and λ_min(K) over that rate per unit time at Part B's noise (σ = {SIG_B:g}, h = {H_F:g}: ds/dt = {ds_dt(SIG_B):.1e}) — "
          f"Katzenberger's limit needs it large —, with the σ at which it would be {RATIO_F:g} (the ratio ∝ 1/σ² at fixed h):", "",
          "| masses | s | ‖w‖₁/a_max | K's spectrum | alive | dying masses: min, max | Q/c dying: min, max | Q/c alive: min, max | "
          f"max\\|d log β/ds\\| | λ_min(K)/rate at σ = {SIG_B:g} | σ for a ratio of {RATIO_F:g} |", "|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in F['along']:
        b, kap = r['beta'], r['kap']
        al = b > 1e-3
        rr = range_ratio(kap, b, r['rate'])
        live = b > LIVE_F
        dm = f"{b[~al].min():.2e}, {b[~al].max():.2e}" if (~al).any() else "—"
        if r['qc'] is None:
            qd = qa = "—"
        else:
            qd = f"{r['qc'][~al].min():.5f}, {r['qc'][~al].max():.5f}" if (~al).any() else "—"
            qa = f"{r['qc'][al].min():.4f}, {r['qc'][al].max():.4f}"
        L_.append(f"| {r['mode']} | {r['s']:g} | {np.abs(r['w']).sum() / AMAX:.4f} | " + ', '.join(f"{k:.1e}" for k in kap) + f" | {al.sum()} | "
                  f"{dm} | {qd} | {qa} | {np.abs(r['rate'][live]).max():.2e} | {rr:.3g} | {SIG_B * np.sqrt(rr / RATIO_F):.2g} |")
    L_.append("")
    for mode in ('free', 'conserved'):
        rs = [r for r in F['along'] if r['mode'] == mode]
        ss = [r['s'] for r in rs]; rr = [range_ratio(r['kap'], r['beta'], r['rate']) for r in rs]
        x1, x10 = crossing(ss, rr, 1.0), crossing(ss, rr, RATIO_F)
        e = rs[-1]
        txt = (f"{mode.capitalize()} masses: the ratio is {rr[0]:.3g} at s = 0 and {rr[ss.index(S0_B)]:.3g} at Part B's s = {S0_B:g}; "
               + (f"it falls below {RATIO_F:g} at s ≈ {x10:.2g} and below 1 at s ≈ {x1:.2g} (log–log interpolation), and is {rr[-1]:.2g} at "
                  f"s = {ss[-1]:.0e}" if x1 is not None else f"it stays above 1 (at least {min(rr):.3g})")
               + f"; at s = {ss[-1]:.0e} {int((e['beta'] > 1e-3).sum())} of {e['beta'].size} atoms are alive, the smallest mass "
               f"{e['beta'].min():.2e}, λ_min(K) = {e['kap'].min():.1e}")
        if mode == 'free' and C is not None:
            rc = [r for r in C if r['k'] == 2 and r['i'] == SEED_B and 'stopped' not in r[('free', 0.0)]]
            if rc:
                txt += (f"; the end is Part C's (k = 2, instance {SEED_B}, free, white) to {np.abs(e['w'] - rc[0][('free', 0.0)]['w']).max():.1e}")
        L_.append(txt + ".")
        L_.append("")
    if F.get('scheme'):
        L_ += ["The paradigm's scheme against the free law's reduced flow, projected as in B (64 paths, fresh draws, h = 1, the flow at the "
               f"step's filter), from the flow's state at s0: a burn-in of {BURN_RELAX_F:g}/λ_min(K) time units, a window of {DS_F:g} in s, the "
               f"limit map (`settle`) run for max(2·10⁵, {SETTLE_RELAX_F:g}/λ_min(K)); λ_min(K) at s0 and the ratio there at the run's σ:", "",
               "| σ | s0 | λ_min(K) at s0 | λ_min(K)/rate at s0 | steps (burn-in + window) | the flow's advance (largest component) | projected "
               "advance / the flow's | the raw mean's offset (z) | projected mean − flow | raw mean − flow | start: projected − flow |",
               "|---|---|---|---|---|---|---|---|---|---|---|"]
        for r in F['scheme']:
            L_.append(f"| {r['sig']:g} | {r['s0']:g} | {r['lam0']:.1e} | {range_ratio(r['kap0'], r['beta0'], r['rate0'], sig=r['sig']):.3g} | "
                      f"{r['n_burn']} + {r['n_win']} | {r['df']:.2e} | {r['ratio']:.3f} ± {r['ratio_se']:.3f} | {r['off']:.1e} ({r['off_z']:.1f}) | "
                      f"{r['dev_proj']:.1e} | {r['dev_raw']:.1e} | {r['start_off']:.1e} |")
        L_.append("")
    L_ += ["Free masses, Q's drift dQ_i/ds = (ω_iγβ_i/(Nm))(κγa_i² − 2Gβ_i(1 − t_i⁴)) (`q_drift`, white) on the free law's curves Q/c, at "
           f"{TS_F.size} amplitudes t = a/a_max from {TS_F[0]:g} to {TS_F[-1]:g} within each curve's range (β > 0): its sign along the curve "
           "for the flat update (κ = 1) and the paradigm's (κ = 2), and its zero (bisection on `q_drift`) against the closed form — κ = 2: "
           "t = 1/κ′ = ((c + Q)/(c − Q))^{1/2} on −c < Q < 0; κ = 1: t = ((−Q − c)/(c − Q))^{1/4} on a branch:", "",
           "| Q/c | the curve | the flat update (κ = 1) | its zero (closed form) | the paradigm's (κ = 2) | its zero (closed form) |",
           "|---|---|---|---|---|---|"]
    for r in F['signs']:
        z = lambda k: (', '.join(f"{x:.6f}" for x in r[k]['zeros']) or '—') + (f" ({r[k]['closed']:.6f})" if r[k]['closed'] is not None else "")
        L_.append(f"| {r['qc']:g} | {curve_name(r['qc'])} | {sign_runs(r['t'], r[1]['sign'])} | {z(1)} | {sign_runs(r['t'], r[2]['sign'])} | {z(2)} |")
    sp = F['splits']
    L_ += ["", "Where each update's Q-drift vanishes on the orbit β|a| = m|w|: D*'s split (the flat update, 2Gm|w|(1 − t⁴) = γa_max³t³) and the "
           "paradigm's gauge equilibrium (Gm|w|(1 − t⁴) = γa_max³t³) — Q/c on states built there (`gauge_invariant`) against the closed forms "
           "−(1 + t⁴)/(1 − t⁴) (a branch for t > 0) and −(1 − t²)/(1 + t²); at the paradigm's equilibrium the Itô term's mass component over "
           "the flat drift's before the projection (`ito_term`, `grad`) against (1 + t²)/(3 − t²):", "",
           "| \\|w\\| | D*'s split: t | Q/c | −(1 + t⁴)/(1 − t⁴) | the paradigm's: t | Q/c | −(1 − t²)/(1 + t²) | Itô/flat | (1 + t²)/(3 − t²) |",
           "|---|---|---|---|---|---|---|---|---|"]
    for j, x in enumerate(sp['xs']):
        tf, tp = sp['flat']['t'][j], sp['par']['t'][j]
        L_.append(f"| {x:g} | {tf:.4f} | {sp['flat']['qc'][j]:.4f} | {-(1 + tf ** 4) / (1 - tf ** 4):.4f} | {tp:.4f} | {sp['par']['qc'][j]:.4f} | "
                  f"{-(1 - tp ** 2) / (1 + tp ** 2):.4f} | {sp['par']['ito_flat'][j]:.4f} | {(1 + tp ** 2) / (3 - tp ** 2):.4f} |")
    nl = F['null']
    rng_ = lambda k: (min(r[k] for r in nl), max(r[k] for r in nl))
    chk = ""
    if A is not None:
        chk = (f" (Part A's points replayed: its white and slow values reproduced to "
               f"{max(max(abs(a['white'] - r['white']), abs(a['slow'] - r['slow'])) for a, r in zip(A['null'], nl)):.1e})")
    L_ += ["", f"The null kernel's potential against the free masses' at Part A's six points{chk}: R_null − R_free from {rng_('white')[0]:.6e} to "
           f"{rng_('white')[1]:.6e} under white noise (spread {np.ptp([r['white'] for r in nl]):.1e}); under fresh draws at a finite step, "
           "Ψ = −(2/h)log(1 − hκ/2) (Ψ′ = 2/(2 − hκ)), " + '; '.join(
               f"h = {h:g} from {rng_(h)[0]:.4e} to {rng_(h)[1]:.4e} (spread {np.ptp([r[h] for r in nl]):.2e})" for h in H_NULL_F)
           + f"; slow noise (τ = {TAU_NULL:g}) from {rng_('slow')[0]:.4e} to {rng_('slow')[1]:.4e} (spread {np.ptp([r['slow'] for r in nl]):.2e}).", ""]
    return L_


def report(A=None, B=None, C=None, D=None, Ee=None, F=None):
    L_ = ["# Proposition D‴.44: the replicator's own masses on a kernel dictionary under measurement noise (`exp_mass_dictionary.py`)", "",
          f"The loop at fixed positions: masses β (Σβ = m = {M_I:g}, or free), amplitudes a = a_max tanh(α/a_max) (a_max = {AMAX:g}), w = β⊙a/m, "
          f"G = {G_I:g}, γ = {GAM_I:g}; label noise on the measured error; the selection along the exact fits (`mass_dictionary`), in s = (m v/4)t. "
          "The paradigm's update: Euler on (α, log β), renormalised; the flat one: Euler on (α, √β), renormalised radially.", ""]
    if A is not None:
        L_ += ["## A. Identities", "",
               "At random exact fits (5 × 12, three points per law): the closed-form drift against Katzenberger's drift of the scheme "
               "(`gauge_drift.limit_drift` on a chart of M), largest relative difference:", "",
               "| masses | the draws | the flat update | the paradigm's update (with its Itô term) |", "|---|---|---|---|"]
        for mode in ('conserved', 'free'):
            for tc in (0.0, 200.0):
                r = lambda c: max(x['rel'] for x in A['rows'] if x['mode'] == mode and x['tc'] == tc and x['coords'] == c)
                L_.append(f"| {mode} | {'fresh, h → 0' if tc == 0 else 'AR(1), τ_c = 200, h = 0.5'} | {r('flat'):.1e} | {r('log'):.1e} |")
        L_ += ["", "The paradigm's drift against the flat drift without the masses' own selection (the γ-term of D_i differentiated in √β), "
               "largest relative difference (the conserved law's leftover is the centring's cross term, `cross_term`; its size along the flows "
               "below): " + '; '.join(
                   f"{mode}, {'fresh' if tc == 0 else 'AR(1)'} {max(x['rel'] for x in A['canc'] if x['mode'] == mode and x['tc'] == tc):.1e}"
                   for mode in ('free', 'conserved') for tc in (0.0, 200.0)) + ".", "",
               "Free masses, the gauge invariant Q_i = β_i − c cosh(2α_i/a_max) (Proposition D‴.45): ∇Q·(Katzenberger's drift) against the "
               "closed form (ω_iγβ_i/(Nm))(κγa_i² − 2Gβ_i(1 − t_i⁴)) (κ = 1 flat, 2 paradigm), largest relative difference: " + '; '.join(
                   f"{'fresh' if tc == 0 else 'AR(1)'}, {'flat' if c == 'flat' else 'paradigm'} "
                   f"{max(x['rel'] for x in A['qd'] if x['tc'] == tc and x['coords'] == c):.1e}" for tc in (0.0, 200.0) for c in ('flat', 'log'))
               + "; Q's quadratic variation under the paradigm's update relative to its scale: "
               + ', '.join(f"{'fresh' if tc == 0 else 'AR(1)'} {max(x['qv'] for x in A['qd'] if x['tc'] == tc and 'qv' in x):.1e}"
                           for tc in (0.0, 200.0)) + ".", "",
               "Loop integrals ∮⟨P⁻¹q, dθ⟩ of the white drifts on M around loops of radius ε (in the metric) at one point per law — zero "
               "for a gradient flow in the flow's metric, ε² times a 2-form otherwise:", "",
               "| masses | the update | " + ' | '.join(f'ε = {e:g}' for e in LOOP_EPS) + " | ratios |", "|---|---|" + "---|" * len(LOOP_EPS) + "---|"]
        for r in A['loops']:
            I = r['I']
            L_.append(f"| {r['mode']} | {'the paradigm' if r['ito'] else 'flat'} | " + ' | '.join(f'{x:.2e}' for x in I) + " | "
                      + ', '.join(f'{I[j] / I[j + 1]:.2f}' for j in range(len(I) - 1)) + " |")
        L_ += ["", "K's spectrum against the flow's normal rates (−eigenvalues of its Jacobian in (α, √β)): largest relative difference "
               + ', '.join(f"{mode} {max(x['rel'] for x in A['spec'] if x['mode'] == mode):.1e}" for mode in ('conserved', 'free')) + ".", "",
               "The excess-mass bound at the conserved law's points (tr K the conserved law's, K_free + its −(γ/N)yyᵀ):", "",
               "| E/m | Σβa′/m | tr K | tr K_free | (G/(Nm))E²/Σ|φ_i|⁻² |", "|---|---|---|---|---|"]
        for r in A['bound']:
            L_.append(f"| {r['E']:.4f} | {r['sba']:.4f} | {r['trK']:.4e} | {r['trKf']:.4e} | {r['lb']:.4e} |")
        L_ += ["", f"The gauge-minimal relevance D*(|w|) against (γa_max/m)|w| and the 4/3 law (the threshold γa_max³/(Gm) = {A['thr']:.4g}); "
               "its minimiser's t = |a|/a_max (the root of 2Gm|w|(1 − t⁴) = γa_max³t³; D at that t against D*) and the paradigm's gauge "
               "equilibrium (Gm|w|(1 − t⁴) = γa_max³t³):", "",
               "| |w| | D* | over (γa_max/m)|w| | over (3/2)(2G)^{1/3}(γ/m)^{2/3}|w|^{4/3} | D*'s t | D(t)/D* − 1 | the paradigm's t |",
               "|---|---|---|---|---|---|---|"]
        for x, dsv, li, f3, tf, da_, tp in zip(A['xs'], A['ds'], A['lin'], A['four3'], A['t_flat'], A['d_at'], A['t_par']):
            L_.append(f"| {x:g} | {dsv:.4e} | {dsv / li:.4f} | {dsv / f3:.4f} | {tf:.4f} | {da_ / dsv - 1:.1e} | {tp:.4f} |")
        wv = [r['white'] for r in A['null']]; sv = [r['slow'] for r in A['null']]
        L_ += ["", f"The null kernel's potential against the free masses' at six points of its manifold (the real atoms' state an exact fit "
               f"of the free law to {max(r['on_M'] for r in A['null']):.1e}): white noise, R_null − R_free from {min(wv):.6e} to {max(wv):.6e} "
               f"(−(γ/N)|y|² = {A['null_const']:.6e}); slow noise (τ = {TAU_NULL:g}), from {min(sv):.4e} to {max(sv):.4e}.", "",
               f"Along Part B's instance's white flows (the paradigm's update), at s = {', '.join(f'{s:g}' for s in S_ALONG)}: the Itô term's "
               "length against the flat drift's; for the conserved law the cross term's against the net drift's and the pushes on ‖w‖₁ "
               "(d‖w‖₁/ds from the cross term, from the rest); for free masses the gauge ratio Gβ(1 − t⁴)/(γa²) on the alive atoms:", "",
               "| masses | s | ‖w‖₁/a_max | Itô/flat | cross/net | push: cross, rest | gauge ratio (alive) |", "|---|---|---|---|---|---|---|"]
        for r in A['along']:
            if r['mode'] == 'conserved':
                L_.append(f"| conserved | {r['s']:g} | {r['l1']:.4f} | {r['ito']:.3f} | {r['cross']:.3g} | {r['push_cross']:.2e}, {r['push_rest']:.2e} | — |")
            else:
                L_.append(f"| free | {r['s']:g} | {r['l1']:.4f} | {r['ito']:.3f} | — | — | {r['gr_lo']:.4f}–{r['gr_hi']:.4f} ({r['alive']}) |")
        L_.append("")
    if B is not None:
        L_ += [f"## B. The paradigm's noisy scheme against the reduced flow, projected (5 × 12, the 2-sparse instance {SEED_B}, {P_B} paths)", "",
               f"From the reduced flow's state at s = {S0_B:g}, a burn-in of {BURN_B:g} time units, then a window of {DS_B:g} in s; every path's "
               "state at the window's ends mapped onto M by the noiseless flow (`settle`, Katzenberger's limit map); the projected advance's "
               "component along the flow's advance over the same window, relative to it (± s.e. over the paths); the raw paths' mean "
               "offset from the projected one at the window's end (largest component; in s.e.); the projected and raw means against the "
               f"flow's state there. Fresh draws at σ = {SIG_B:g}; AR(1) draws at τ_c = 10³ (in time units) with σ set by σ²τ_c.", "",
               "| masses | the draws | h | σ²τ_c/h | steps (burn-in + window) | the flow's advance (largest component) | projected advance / the flow's | "
               "the raw mean's offset (z) | projected mean − flow | raw mean − flow | start: projected − flow |", "|---|---|---|---|---|---|---|---|---|---|---|"]
        for r in B:
            noise = 'fresh' if r['tc'] == 0 else f"AR(1), σ²τ_c = {r['s2tc']:g}"
            q = '—' if r['tc'] == 0 else f"{r['s2tc'] / r['h']:g}"
            L_.append(f"| {r['mode']} | {noise} | {r['h']:g} | {q} | {r['n_burn']} + {r['n_win']} | {r['df']:.2e} | {r['ratio']:.3f} ± {r['ratio_se']:.3f} | "
                      f"{r['off']:.1e} ({r['off_z']:.1f}) | {r['dev_proj']:.1e} | {r['dev_raw']:.1e} | {r['start_off']:.1e} |")
        L_.append("")
    for name, R, ks in (('C', C, KS_C), ('D', D, KS_D)):
        if R is None:
            continue
        stopped = []
        hdr = (f"## C. The selection's ends on Gaussian designs ({N_I} × {D_I}, {N_INST_C} instances per k)" if name == 'C' else
               f"## D. The selection's ends on Gaussian kernel dictionaries ({D_D} centres, {N_D} samples on [0, 1], width {WIDTH_D:g}, "
               f"{N_INST_D} instances per k)")
        L_ += [hdr, "", f"The paradigm's update unless marked; reduced flows to s = {S_END:.0e}; the start the noiseless flow's end from equal "
               "masses and amplitudes N(0, 0.05²). Each row over the instances whose flow finished (n; a stopped flow is listed below with "
               "its last state): means of ‖w‖₁/a_max (the truth's over the same instances), the excess E/m over the real atoms (conserved "
               "law and null kernel), the atoms used (|w| > 10⁻³; the largest count in brackets), alive (β > 10⁻³); the fractions recovering "
               f"the truth (within 10⁻⁴·max|w*|) and within {TOL_BP:g} of basis pursuit's solution, on the instances where basis pursuit is "
               "the truth / where it is not; the largest distance to basis pursuit; the null kernel's mass.", "",
               "| k | masses | the noise | n | ‖w‖₁/a_max (truth's) | E/m | used (max) | alive | the truth: bp ok / bp fails | at basis pursuit: "
               "bp ok / bp fails | largest distance to basis pursuit | the null kernel's mass |", "|---|---|---|---|---|---|---|---|---|---|---|---|"]
        for k in ks:
            rs = [r for r in R if r['k'] == k]
            if not rs:
                continue
            nok = sum(r['bp_ok'] for r in rs)
            L_.append(f"| {k} | — | basis pursuit | {len(rs)} | {np.mean([r['l1_bp'] for r in rs]):.3f} ({np.mean([r['l1_truth'] for r in rs]):.3f}) | | | | "
                      f"{nok}/{len(rs)} | | | |")
            for mode in MODES:
                for key in ['start'] + list(TAUS) + (['flat'] if mode != 'null kernel' else []):
                    fin = [r for r in rs if 'stopped' not in r[(mode, key)]]
                    for r in rs:
                        if 'stopped' in r[(mode, key)]:
                            x = r[(mode, key)]; la = x['last']
                            stopped.append(f"k = {k}, instance {r['i']}, {mode}, {key}: {x['stopped']} at s = {x['s']:.3g}"
                                           + (f" (its last state: ‖w‖₁/a_max {la['l1']:.3f}, the truth {'recovered' if la['ok'] else 'not recovered'} "
                                              f"(distance {la['d_truth']:.1e}), distance to basis pursuit {la['d_bp']:.1e}, {la['used']} atoms used, "
                                              f"{la['alive']} alive)" if la is not None else ""))
                    q = [r[(mode, key)] for r in fin]
                    lab = {'start': 'the start', 'flat': 'white, the flat update', 0.0: 'white'}.get(key, None) or f"slow, τ_c = {key:.0e}"
                    if not q:
                        L_.append(f"| {k} | {mode} | {lab} | 0 | | | | | | | | |"); continue
                    ok = [x['ok'] for r, x in zip(fin, q) if r['bp_ok']]; fl = [x['ok'] for r, x in zip(fin, q) if not r['bp_ok']]
                    at_ok = [x['d_bp'] < TOL_BP for r, x in zip(fin, q) if r['bp_ok']]; at_fl = [x['d_bp'] < TOL_BP for r, x in zip(fin, q) if not r['bp_ok']]
                    L_.append(f"| {k} | {mode} | {lab} | {len(q)} | {np.mean([x['l1'] for x in q]):.3f} ({np.mean([r['l1_truth'] for r in fin]):.3f}) | "
                              + (f"{np.mean([x['exc'] for x in q]):.3f}" if mode != 'free' else "—") + f" | {np.mean([x['used'] for x in q]):.1f} "
                              f"({max(x['used'] for x in q)}) | {np.mean([x['alive'] for x in q]):.1f} | {sum(ok)}/{len(ok)}, {sum(fl)}/{len(fl)} | "
                              f"{sum(at_ok)}/{len(at_ok)}, {sum(at_fl)}/{len(at_fl)} | {max(x['d_bp'] for x in q):.3f} | "
                              + (f"{np.mean([x['b0'] for x in q]):.2f} |" if mode == 'null kernel' else "— |"))
        L_ += ["", "The null kernel against free masses, instance by instance (the paradigm's update, the same colour): the instances whose ends "
               f"differ (largest coefficient's distance above {TOL_BP:g}); the instances where the null kernel misses the truth while basis "
               "pursuit and free masses recover it:", "", "| k | the noise | ends differ | the null kernel misses the truth where both recover it |",
               "|---|---|---|---|"]
        for k in ks:
            rs = [r for r in R if r['k'] == k]
            if not rs:
                continue
            for tau in TAUS:
                both = [r for r in rs if 'stopped' not in r[('null kernel', tau)] and 'stopped' not in r[('free', tau)]]
                diff = [r['i'] for r in both if np.abs(r[('null kernel', tau)]['w'] - r[('free', tau)]['w']).max() > TOL_BP]
                miss = [r['i'] for r in both if r['bp_ok'] and r[('free', tau)]['ok'] and not r[('null kernel', tau)]['ok']]
                lab = 'white' if tau == 0 else f"slow, τ_c = {tau:.0e}"
                L_.append(f"| {k} | {lab} | {len(diff)} of {len(both)} ({', '.join(map(str, diff)) or '—'}) | "
                          f"{len(miss)} of {sum(r['bp_ok'] and r[('free', tau)]['ok'] for r in both)} ({', '.join(map(str, miss)) or '—'}) |")
        L_ += ["", "Free masses: the gauge ratio Gβ(1 − t⁴)/(γa²) on the alive atoms at the ends (1 at the paradigm's gauge equilibrium, 1/2 at "
               "D*'s split), median [range] over the atoms of all finished flows; the flat update's end against the minimiser of "
               "Σ|φ_i|²D*(|w_i|), that minimiser against basis pursuit's solution (median, largest), and the flat ends on the truth:", "",
               "| k | white | slow, τ_c = 1e+03 | slow, τ_c = 1e+04 | white, the flat update | flat end − D* minimiser | D* minimiser − "
               "basis pursuit | flat ends on the truth |", "|---|---|---|---|---|---|---|---|"]
        for k in ks:
            rs = [r for r in R if r['k'] == k]
            if not rs:
                continue
            g = lambda key: np.concatenate([r[('free', key)]['gr'] for r in rs if 'stopped' not in r[('free', key)]] or [np.zeros(0)])
            dvf = [r['dstar_vs_flat'] for r in rs if np.isfinite(r['dstar_vs_flat'])] or [np.nan]
            fo = [r for r in rs if 'stopped' not in r[('free', 'flat')]]
            L_.append(f"| {k} | " + ' | '.join(_fmt_range(g(key)) for key in list(TAUS) + ['flat']) + f" | {np.median(dvf):.1e}, {max(dvf):.1e} | "
                      f"{np.median([r['dstar_vs_bp'] for r in rs]):.1e}, {max(r['dstar_vs_bp'] for r in rs):.1e} | "
                      f"{sum(r[('free', 'flat')]['ok'] for r in fo)} of {len(fo)} |")
        if stopped:
            L_ += ["", f"Flows stopped (over the {FLOW_SECS:g}-s limit, or failed), left out of the rows: " + '; '.join(stopped) + "."]
        L_.append("")
    if Ee is not None:
        L_ += [f"## E. The conserved law's white end against the representability boundary ({N_I} × {D_I}, 2-sparse truths, {N_INST_E} instances; "
               "those whose basis-pursuit solution is not representable, ‖w_bp‖₁ ≥ a_max, left out)", "",
               "(iii)'s bound E/m ≤ (a_max/m)(γmΣ|φ_i|⁻²max|φ_i|²/G)^{1/2} (vacuous at 1 and above); the largest amplitude |a|/a_max at the end "
               "(median over the instances); the saturation threshold γa_max³/(Gm) against the median used coefficient |w_i| at the end:", "",
               "| a_max | instances | ‖w‖₁/a_max at the start | at the end (min) | E/m at the end | the bound | max|a|/a_max | γa_max³/(Gm), median used |w_i| "
               "| the atoms used, alive | ‖w_bp‖₁/a_max |", "|---|---|---|---|---|---|---|---|---|---|"]
        for am in AMAXES_E:
            rs = [r for r in Ee if r['am'] == am and not r['skip']]
            n_st = sum('stopped' in r for r in Ee if r['am'] == am)
            if n_st:
                L_.append(f"| {am:g} | ({n_st} stopped, left out) | | | | | | | | |")
            if not rs:
                L_.append(f"| {am:g} | 0 | | | | | | | | |"); continue
            uw = np.concatenate([np.abs(r['end']['w'])[np.abs(r['end']['w']) > 1e-3] for r in rs])
            L_.append(f"| {am:g} | {len(rs)} | {np.mean([r['start']['l1'] for r in rs]):.3f} | {np.mean([r['end']['l1'] for r in rs]):.3f} "
                      f"({np.min([r['end']['l1'] for r in rs]):.3f}) | {np.mean([r['end']['exc'] for r in rs]):.4f} | "
                      f"{np.median([r['bound'] for r in rs]):.3f} | {np.median([r['end']['amax'] for r in rs]):.3f} | "
                      f"{GAM_I * am ** 3 / (G_I * M_I):.4f}, {np.median(uw):.3f} | "
                      f"{np.mean([r['end']['used'] for r in rs]):.1f}, {np.mean([r['end']['alive'] for r in rs]):.1f} | "
                      f"{np.mean([r['l1_bp'] / am for r in rs]):.3f} |")
        L_.append("")
    if F is not None:
        L_ += report_F(F, A, C)
    L_.append("Numbers are produced by `experiments/exp_mass_dictionary.py`; nothing is edited by hand.")
    return '\n'.join(L_) + '\n'


def main():
    raw = os.path.join(OUT, 'mass_dictionary_raw.npz')
    md = os.path.join(OUT, 'mass_dictionary.md')
    old = E.load(raw)
    get = lambda k: old[k] if k in old else None
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(*[get(k) for k in 'ABCDEF']))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else set('ABCDEF')
    from multiprocessing import Pool
    res = dict(old)
    t0 = time.time()
    log = lambda k: print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
    with Pool(N_WORKERS) as pool:
        h = {}
        if 'B' in parts:                                                # the longest jobs first
            jobs = [(mo,) + r for r in B_RUNS for mo in MODES]
            cost = lambda j: 0 if j[1] == 0 else 1.0 / (j[2] / j[1] * j[3] ** 2)
            h['B'] = pool.map_async(b_job, sorted(jobs, key=cost, reverse=True), chunksize=1)
        if 'C' in parts:
            h['C'] = pool.map_async(c_job, [('gauss', k, i) for k in KS_C for i in range(N_INST_C)], chunksize=1)
        if 'D' in parts:
            h['D'] = pool.map_async(c_job, [('rbf', k, i) for k in KS_D for i in range(N_INST_D)], chunksize=1)
        if 'E' in parts:
            h['E'] = pool.map_async(e_job, [(i, am) for i in range(N_INST_E) for am in AMAXES_E], chunksize=1)
        if 'F' in parts:
            h['F'] = pool.map_async(f_job, [(sg, s0) for s0 in S0_F for sg in SIGS_F], chunksize=1)
        if 'A' in parts:
            res['A'] = part_A(); E.save(raw, **res); log('A')
        if 'F' in parts:
            fd = part_F(); log('F (flows, signs, splits, null kernel)')
        for k in 'BCDEF':
            if k in h:
                res[k] = h[k].get()
                if k == 'B':
                    res[k] = sorted(res[k], key=lambda r: (MODES.index(r['mode']), [x[:3] for x in B_RUNS].index((r['tc'], r['s2tc'], r['h']))))
                if k == 'F':
                    res[k] = dict(fd, scheme=res[k])
                E.save(raw, **res); log(k)
    open(md, 'w').write(report(*[res.get(k) for k in 'ABCDEF']))


if __name__ == '__main__':
    main()
