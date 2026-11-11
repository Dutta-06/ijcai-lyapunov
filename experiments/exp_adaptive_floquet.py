"""
Corollary D‴.26: perfect tracking with an adaptive prior in the real closed loop — the Floquet spectrum, uniformly in m
========================================================================================================================

The three-bump instance (σ₀ = 0.15σ_crit, the reference x_m = 2 sin t, T₀ = 2π), the nearest-atom leak toward an adaptive prior in the real
closed loop of Theorem J (src/core/adaptive_floquet.py), and the sampled weights w = (6, 5, 4)/15, representable (and a saddle of the fixed
prior, Proposition D‴.24's numerics): the lift of their exact
representation ϑ^e(w) (the amplitudes A_k/w_k) with the prior at its atoms and the plant on the reference is a T₀-periodic solution with
e ≡ 0 (Proposition D‴.24(vi)). For the weighted K₀-atom loop (one particle per class, the weights w) and the lifts m = 15 and 30 (the
counts 6, 5, 4 times 1 and 2), η = σ₀ and ∞, and k_e = 25 and 50: the return residual of one period, the monodromy matrix by central
differences of the period map (step 10⁻⁵ relative, DOP853 at rtol 10⁻¹³) and its Floquet exponents (rates −log|λ|/T₀). Checked: the
separations' exponents are σ₀/k_e exactly ((m − K₀)p of them), the lifts' other exponents are the K₀-atom loop's, and the K₀-atom loop's are
the averaged prediction
(1/k_e)λ_j^±(w) of Theorem D‴.23(ii) at the weights w (η = ∞: h_j(w)/k_e) up to the averaging correction, h_j(w) the whitened weighted
Gauss–Newton spectrum at ϑ^e(w); the correction's scaling from k_e = 25 to 50.

Part C (`--only-C`): the variational equation along the lift on the collapsed subspace (src/core/adaptive_floquet.CollapsedVariational: the
whitened weighted coordinates, the monodromy by integrating it, no differences) — checked against the real loop's linearisation and against
Part A's differences — at k_e = 12.5 … 400, for x_m = 2 sin t (where the first-order averaging correction Ā₁ vanishes identically: in one
dimension and by the symmetry t → π − t) and for the asymmetric x_m = 2 sin t + 0.8 sin 2t (Ā₁ ≠ 0, skew): the deviations of the collective
exponents from the first-order prediction λ_j^±(w)/k_e (O(k_e⁻²) relative) and from the second-order one (O(k_e⁻⁴): the odd orders vanish),
the second-order coefficients c_j and their error-filter part, and the error channel's exponent by Liouville's formula.

Part D (`--add-D`, Lemma D‴.26a): on the asymmetric reference at η = 0 (the fixed prior at the atoms: K = Ĵ, Λ = σ₀I), σ₀ and ∞, k_e = 12.5
… 1600 — Chang's slow subspace P(t; ε) by its contraction (the iterations and the invariance residual, the truncations |P − P₂| against ε³ and
|P − P₍₃₎| against ε⁴), the slow block's exponents against the full variational equation's, the slow generator at ±ε (the exponents equal:
the parity identity; their first-order deviation for comparison, relative to the largest rate and per mode), the same with a derivative
mismatch in the error channel (the parity broken at first order), the Magnus parameter; and under white measurement noise the second
transformation's decoupling residual, the slow noise b against b₀ = −(K̇ + aK)ᵀ and its next term εb₁, the exact periodic covariance at
t = 0 in the coordinates (v, Y′) against the lemma's leading terms, the weakest mode's variance against Σ̄ split into the slow noise's and the
rate's corrections, Cov(Y) − Cov(Y′) against the jitter and its next term, and the floor's period average and range over the period against
Proposition D‴.30(ii)'s formula (referee round 22, D3s-1 to D3s-3, D3s-6, D3s-8).

Output: results/chaos/adaptive_floquet.md and adaptive_floquet_raw.npz (AFL_WORKERS; `--from-raw` regenerates the md).
"""
import os
import sys
import time
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.core import langevin as L
from src.core import ot_leak as OT
from src.core import adaptive_prior as AP
from src.core import adaptive_floquet as AF
from src.core import measurement_noise as MN
from src.core.periodic import floquet_exponents

OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos')


def save_raw(path, **arrays):
    """np.savez to a temporary file, then an atomic rename: an interrupted save never truncates the raw file."""
    tmp = path + '.tmp.npz'
    np.savez(tmp, **arrays)
    os.replace(tmp, path)
NET_KW = dict(a_max=4.0, c_max=3.0, s_min=0.3, s_max=2.0)
G = np.array([1.0, 0.25, 0.25])
BUMPS = [(1.0, -1.5, 0.5), (-0.7, 0.0, 0.4), (0.8, 1.5, 0.5)]
PER = 2 * np.pi
SIGMA_CRIT = 0.047162567584665846
SIGMA0 = 0.15 * SIGMA_CRIT
COUNTS = np.array([6, 5, 4])
MULTS = (0, 1, 2)                                                          # 0: the weighted K₀-atom loop
ETAS = (SIGMA0, np.inf)
KES = (25.0, 50.0)
KES_C = (12.5, 25.0, 50.0, 100.0, 200.0, 400.0)                           # Part C
REFS_C = ('sym', 'asym')
N_WORKERS = int(os.environ.get('AFL_WORKERS', 3))


def x_ref_C(key):
    if key == 'sym':
        return lambda t: np.array([2.0 * np.sin(t)])
    return lambda t: np.array([2.0 * np.sin(t) + 0.8 * np.sin(2.0 * t)])


def x_ref_dot_C(key):
    if key == 'sym':
        return lambda t: np.array([2.0 * np.cos(t)])
    return lambda t: np.array([2.0 * np.cos(t) + 1.6 * np.cos(2.0 * t)])


REF_NAME = {'sym': '2 sin t', 'asym': '2 sin t + 0.8 sin 2t'}


def f_vec(X):
    X = np.asarray(X, float)
    return sum(A * np.exp(-(X - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS)


def setup():
    net = L.SaturatedLiftedNetwork(1, **NET_KW)
    w = COUNTS / COUNTS.sum()
    A = np.array([b[0] for b in BUMPS]); C = np.array([b[1] for b in BUMPS]); S = np.array([b[2] for b in BUMPS])
    exact_w = net.unconstrained((A / w)[:, None], C[:, None], S)
    q_exact = net.unconstrained((3 * A)[:, None], C[:, None], S)
    loop = OT.OTPriorLoop(net, lambda X: f_vec(X), lambda t: np.array([2.0 * np.sin(t)]), PER, G, SIGMA0, q_exact, n_quad=256)
    return net, w, exact_w, loop


def job(args):
    mult, eta, ke = args
    t0 = time.time()
    net, w, exact_w, loop = setup()
    if mult == 0:
        labels, weights = np.arange(3), w
    else:
        labels, weights = np.repeat(np.arange(3), COUNTS * mult), None
    m = len(labels)
    cl = AF.AdaptiveRealLoop(net, lambda X: f_vec(X), lambda t: np.array([2.0 * np.sin(t)]), lambda t: np.array([2.0 * np.cos(t)]), PER, ke, G,
                             SIGMA0, labels, 3, eta, weights=weights)
    z_p = cl.lift(exact_w)
    residual = float(np.abs(AF.period_map(cl, z_p, m) - z_p).max())
    DP = AF.monodromy(cl, z_p, m)
    fast, slow = floquet_exponents(DP, PER, 1)
    h = AP.gauss_newton_block(loop, exact_w, w)[1]
    pred = AF.predicted_slow_exponents(h, SIGMA0, eta, ke, m, 3, 3)
    print(f"[{time.time() - t0:.0f}s] m {m} eta {eta:g} ke {ke:g}", flush=True)
    return dict(mult=mult, m=m, eta=eta, ke=ke, residual=residual, fast=fast, slow=slow, pred=pred, h=h)


def job_C(args):
    """Part C: one reference and one η — the variational exponents at KES_C against the first- and second-order predictions, the coefficients
    c_j, the norms of Ā₁ (and of its self-adjoint part in the metric diag(I, (σ₀/η)I)), and the variational generator against the central
    differences of the real loop's field (AdaptiveRealLoop.rhs, the weighted K₀-atom loop with the prior in the state) at six times."""
    key, eta = args
    t0 = time.time()
    net, w, exact_w, loop = setup()
    cv = AF.CollapsedVariational(net, x_ref_C(key), PER, G, SIGMA0, eta, exact_w, w)
    gens = cv.averaged()
    A0, A1, A2, H, filt = gens
    d = cv.d
    Wm = np.eye(A0.shape[0]) if np.isinf(eta) else np.diag(np.r_[np.ones(d), np.full(d, SIGMA0 / eta)])
    lam, c, cf = cv.second_order(gens=gens)
    rows = []
    for ke in KES_C:
        slow, err = cv.exponents(ke)
        rows.append(dict(ke=ke, slow=slow, err=err, p1=cv.predicted(ke, 1, gens=gens), p2=cv.predicted(ke, 2, gens=gens)))
    # the generator against the real loop's linearisation (the K₀-atom loop with the weights w; whitened coordinates)
    cl = AF.AdaptiveRealLoop(net, lambda X: f_vec(X), x_ref_C(key), x_ref_dot_C(key), PER, KES[0], G, SIGMA0, np.arange(3), 3, eta, weights=w)
    Dw = (np.sqrt(w)[:, None] / np.sqrt(G)[None, :]).ravel()
    T = np.concatenate([[1.0], Dw] + ([] if np.isinf(eta) else [Dw]))
    lin = 0.0
    for t in np.linspace(0.0, PER, 6, endpoint=False):
        z = cl.lift(exact_w); z[0] = x_ref_C(key)(t)[0]
        Dz = len(z); Jr = np.zeros((Dz, Dz))
        for j in range(Dz):
            hh = 1e-6 * max(1.0, abs(z[j])); e = np.zeros(Dz); e[j] = hh
            Jr[:, j] = (cl.rhs(t, z + e, 3) - cl.rhs(t, z - e, 3)) / (2 * hh)
        Mw = T[:, None] * Jr / T[None, :]
        M = cv.generator(t, KES[0])
        lin = max(lin, float(np.abs(Mw - M).max() / np.abs(M).max()))
    print(f"[{time.time() - t0:.0f}s] C {key} eta {eta:g}", flush=True)
    return dict(key=key, eta=eta, h=np.linalg.eigvalsh(H), lam=lam, c=c, cf=cf, rows=rows, A1=float(np.abs(A1).max()),
                A1_selfadj=float(np.abs(Wm @ A1 + (Wm @ A1).T).max()), lin=lin, trH=float(np.trace(H)))


def derived_C(key, eta, G_=None):
    """The η-independent mean squared frequencies of the modes in the regressor, ω_j² = v_jᵀ⟨Ĵ̇ᵀĴ̇⟩v_j/h_j (v_j, h_j the eigenpairs of H), and
    Ā₀'s self-adjointness: |WĀ₀ − (WĀ₀)ᵀ| in the metric W = diag(I, (σ₀/η)I) and |Ā₀ − Ā₀ᵀ| in I."""
    net, w, exact_w, loop = setup()
    cv = AF.CollapsedVariational(net, x_ref_C(key), PER, G if G_ is None else G_, SIGMA0, eta, exact_w, w)
    A0, A1, A2, H, filt = cv.averaged()
    h, V = np.linalg.eigh(H)
    om2 = np.einsum('ij,ik,kj->j', V, filt, V) / h
    d = cv.d
    Wm = np.eye(A0.shape[0]) if np.isinf(eta) else np.diag(np.r_[np.ones(d), np.full(d, SIGMA0 / eta)])
    return dict(omega2=om2, A0_metric=float(np.abs(Wm @ A0 - (Wm @ A0).T).max()), A0_plain=float(np.abs(A0 - A0.T).max()), A0_norm=float(np.abs(A0).max()),
                A1_skew_metric=float(np.abs(Wm @ A1 + (Wm @ A1).T).max()), A1_norm=float(np.abs(A1).max()))


GAIN_C = 30.0
KES_GAIN = (100.0, 200.0, 400.0)


def job_C_gain(args):
    """The adaptation gain scaled, G → GAIN_C·G (the symmetric reference, η = ∞): the second-order coefficients c_j (both signs) and the
    variational exponents at KES_GAIN against the first-order prediction — the sign of c_j is not fixed by Corollary D‴.26(c)."""
    key, eta = args
    t0 = time.time()
    net, w, exact_w, loop = setup()
    cv = AF.CollapsedVariational(net, x_ref_C(key), PER, GAIN_C * G, SIGMA0, eta, exact_w, w)
    gens = cv.averaged()
    lam, c, cf = cv.second_order(gens=gens)
    rows = []
    for ke in KES_GAIN:
        slow, err = cv.exponents(ke)
        rows.append(dict(ke=ke, slow=slow, p1=cv.predicted(ke, 1, gens=gens), p2=cv.predicted(ke, 2, gens=gens)))
    print(f"[{time.time() - t0:.0f}s] C gain {GAIN_C:g} {key} eta {eta:g}", flush=True)
    return dict(key=key, eta=eta, gain=GAIN_C, lam=lam, c=c, cf=cf, rows=rows)


def split(r):
    """The (m − K₀)p exponents nearest σ₀/k_e (the separations) and the others (the collective ones), both sorted."""
    s = np.sort(r['slow']); n_sep = (r['m'] - 3) * 3
    order = np.argsort(np.abs(s - SIGMA0 / r['ke']))
    sep = np.sort(s[order[:n_sep]]); coll = np.sort(s[order[n_sep:]])
    return sep, coll


# ---------------------------------------------------------------------------- D: Lemma D‴.26a
KES_D = (12.5, 25.0, 50.0, 100.0, 200.0, 400.0, 800.0, 1600.0)


def fprime(x):
    return sum(-A * (x - C) / S ** 2 * np.exp(-(x - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS)


def job_D(args):
    """Part D at one η on the asymmetric reference: per k_e the Chang quantities, the slow block against the full equation, the parity
    identity (relative to the largest rate and per mode) and its breaking by a mismatch, the Magnus parameter; and under white noise the
    decoupling, the slow noise to first order, the covariances at t = 0, the weakest mode's variance split into its corrections, the jitter
    with its next term, and the floor's period average and range against Proposition D‴.30(ii)'s formula."""
    eta = args
    t0 = time.time()
    net, w, exact_w, loop = setup()
    cv = AF.CollapsedVariational(net, x_ref_C('asym'), PER, G, SIGMA0, eta, exact_w, w)
    tg, K, Lam, om = cv.chang_grid(512)
    der = lambda X, k=1: np.real(np.fft.ifft(((1j * om) ** k).reshape(-1, 1, 1) * np.fft.fft(X, axis=0), axis=0))
    a_mis = np.array([fprime(float(cv.xr(t)[0])) for t in tg])
    A0 = cv.averaged()[0]
    nu0 = np.sort(np.linalg.eigvals(-A0).real)
    Kd, Kdd, Kddd = der(K), der(K, 2), der(K, 3)
    KKtK = lambda X, Y, Z: np.einsum('tij,tkj,tkl->til', X, Y, Z)
    P3 = Kddd + 2 * KKtK(Kd, K, K) + KKtK(K, Kd, K) + 2 * KKtK(K, K, Kd) + 2 * Kd @ Lam
    H, Q = MN.noise_matrices(cv, fprime)
    d = cv.d
    hv, Vh = np.linalg.eigh(H)
    ew = np.zeros(A0.shape[0]); ew[:d] = Vh[:, 0]                           # the weakest mode on the y block
    trH2 = float(np.trace(H @ H))
    rows = []
    for ke in KES_D:
        eps = 1.0 / ke
        P, A, it, res = cv.chang(eps)
        P2 = -K + eps * Kd - eps ** 2 * (Kdd + KKtK(K, K, K) + K @ Lam)
        rate = lambda Ab: np.sort(np.linalg.eigvals(-Ab).real)
        Ap, mag = cv.slow_generator(eps); Am = cv.slow_generator(-eps)[0]
        sp, sm = rate(Ap), rate(Am)
        full = np.sort(cv.exponents(ke)[0]) * ke
        mp_, mm_ = rate(cv.slow_generator(eps, mismatch=a_mis)[0]), rate(cv.slow_generator(-eps, mismatch=a_mis)[0])
        row = dict(ke=ke, it=it, res=res, trunc=float(np.abs(P - P2).max() / eps ** 3), trunc3=float(np.abs(P - P2 - eps ** 3 * P3).max() / eps ** 4),
                   vs_full=float(np.abs(sp / full - 1).max()),
                   even=float(np.abs(sp - sm).max() / np.abs(sp).max()), first=float(np.abs(sp - nu0).max() / np.abs(nu0).max()),
                   even_mode=float(np.max(np.abs(sp - sm) / np.abs(sp))), first_mode=float(np.max(np.abs(sp - nu0) / np.abs(nu0))),
                   first_slow=float((sp[0] - nu0[0]) / nu0[0]),
                   even_mis=float(np.abs(mp_ - mm_).max() / np.abs(mp_).max()), even_mis_mode=float(np.max(np.abs(mp_ - mm_) / np.abs(mp_))),
                   magnus=mag)
        P_, S_, b, b0, coup, b1 = MN.noise_decoupling(cv, fprime, eps, with_b1=True)
        o = MN.decoupled_covariance(cv, fprime, ke, period_average=True)
        Sig = o['Sigma']
        Qbar = np.zeros_like(A0); Qbar[:d, :d] = Q
        bb = np.einsum('tin,tjn->ij', b, b) / b.shape[0]
        Sig_noise = MN.solve_continuous_lyapunov(A0, -eps * bb)
        Sig_rate = MN.solve_continuous_lyapunov(Ap, -eps * Qbar)
        wv = lambda X: float(ew @ X @ ew)
        formula = float(MN.floor_closed_form(H, Q, SIGMA0, eta, 1.0, ke)[0]) + trH2 / (2 * ke)
        jit = o['fast']; dYY = o['CY'] - o['CYp']
        row.update(coupling=coup, b_rel=float(np.abs(b - b0).max() / np.abs(b0).max()), b1=float(np.abs(b - b0 - eps * b1).max() / eps ** 2),
                   Vv=float(o['Vv'][0, 0] * 2 * eps ** 3 - 1),
                   C_rel=float(np.abs(o['C'] + o['b'].T).max() / np.abs(o['b']).max()),
                   CYp_rel=float(np.abs(o['CYp'] - Sig).max() / np.abs(Sig).max()),
                   weak_exact=wv(o['CYp']) / wv(Sig) - 1, weak_noise=wv(Sig_noise) / wv(Sig) - 1, weak_rate=wv(Sig_rate) / wv(Sig) - 1,
                   jit_rel=float(np.abs(dYY - jit).max() / np.abs(jit).max()), jit_next=float(np.abs(dYY - jit - o['next']).max() / np.abs(jit).max()),
                   floor0=float(np.trace(H @ o['CY'][:d, :d])), floor0_pred=float(np.trace(H @ (Sig + jit)[:d, :d])),
                   fast_share=float(np.trace(H @ jit[:d, :d]) / np.trace(H @ (Sig + jit)[:d, :d])),
                   floor_avg=o['floor_avg'] / formula, floor_min=o['floor_min'] / formula, floor_max=o['floor_max'] / formula,
                   pred_min=o['pred_min'] / formula, pred_max=o['pred_max'] / formula)
        rows.append(row)
        print(f"[{time.time() - t0:.0f}s] D eta {eta:g} ke {ke:g}", flush=True)
    return dict(eta=eta, rows=rows)


def report(R, t0=None, RC=None, RG=None, RD=None):
    eta_s = lambda e: '∞' if np.isinf(e) else f'{e:.5f}'
    net, w, exact_w, loop = setup()
    lines = ["# Corollary D‴.26: perfect tracking with an adaptive prior in the real loop — the Floquet spectrum (`exp_adaptive_floquet.py`)", "",
             f"The three-bump instance at σ₀ = 0.15σ_crit = {SIGMA0:.5f}, x_m = 2 sin t, the nearest-atom leak toward an adaptive prior in the real closed loop; the sampled "
             f"weights w = {tuple(int(c) for c in COUNTS)}/{int(COUNTS.sum())} (the amplitudes A_k/w_k = " + ", ".join(f"{b[0] / wk:.3f}" for b, wk in zip(BUMPS, w)) + "), the "
             "lift of their exact representation with the prior at its atoms and the plant on the reference. Floquet exponents: rates −log|λ|/T₀ of the monodromy matrix (central "
             "differences of the period map, step 10⁻⁵ relative, DOP853 at rtol 10⁻¹³); the error channel's multiplier (≈ e^{−k_eT₀}) is below the resolution of the differences and is left out.", "",
             "The whitened weighted Gauss–Newton spectrum at ϑ^e(w): h_j(w) = " + ", ".join(f"{x:.4e}" for x in R[0]['h']) + ".", "",
             "## The weighted K₀-atom loop against the averaged prediction", "",
             "| η | k_e | return residual | collective exponents × k_e (actual / predicted (1/k_e)λ_j^±(w) × k_e) | largest relative deviation |", "|---|---|---|---|---|"]
    K = {(r['eta'], r['ke']): r for r in R if r['mult'] == 0}
    for (eta, ke), r in sorted(K.items(), key=lambda t: (np.isinf(t[0][0]), t[0][1])):
        a = np.sort(r['slow']); p = np.sort(r['pred'])
        lines.append(f"| {eta_s(eta)} | {ke:g} | {r['residual']:.1e} | " + ", ".join(f"{x * ke:.4e}/{y * ke:.4e}" for x, y in zip(a, p)) + f" | {np.max(np.abs(a - p) / p):.2e} |")
    lines += [""]
    for eta in ETAS:
        if (eta, KES[0]) in K and (eta, KES[1]) in K:
            d = [np.max(np.abs(np.sort(K[(eta, ke)]['slow']) - np.sort(K[(eta, ke)]['pred'])) * ke) for ke in KES]
            rl = [np.max(np.abs(np.sort(K[(eta, ke)]['slow']) - np.sort(K[(eta, ke)]['pred'])) / np.sort(K[(eta, ke)]['pred'])) for ke in KES]
            lines.append(f"η = {eta_s(eta)}: the largest deviation of the slow-time rates (exponent × k_e) from λ_j^±(w) is {d[0]:.2e} at k_e = {KES[0]:g} and {d[1]:.2e} at "
                         f"k_e = {KES[1]:g} (ratio {d[0] / d[1]:.2f}, i.e. k_e^{{−{np.log(d[0] / d[1]) / np.log(KES[1] / KES[0]):.2f}}}; relative {rl[0]:.2e} and {rl[1]:.2e}) — the "
                         "averaging correction is of second order in 1/k_e in slow time, so the exponents are (1/k_e)λ_j^±(w) + O(k_e⁻³) (Corollary D‴.26(c)).")
    lines += ["", "## The lifts: separations at σ₀/k_e, the collective exponents the K₀-atom loop's", "",
              "| m | η | k_e | return residual | separations: count, largest relative distance from σ₀/k_e | collective exponents: count, largest relative difference from the K₀-atom loop's |",
              "|---|---|---|---|---|---|"]
    for r in sorted([r for r in R if r['mult'] > 0], key=lambda r: (np.isinf(r['eta']), r['ke'], r['m'])):
        sep, coll = split(r)
        ref = np.sort(K[(r['eta'], r['ke'])]['slow'])
        lines.append(f"| {r['m']} | {eta_s(r['eta'])} | {r['ke']:g} | {r['residual']:.1e} | {len(sep)}, {np.max(np.abs(sep * r['ke'] / SIGMA0 - 1)):.1e} | "
                     f"{len(coll)}, {np.max(np.abs(coll - ref) / ref):.1e} |")
    if RC:
        lines += ["", "## C. The variational equation: the remainder of (c), to k_e = 400", "",
                  "The variational equation along the lift on the collapsed subspace, in the whitened weighted coordinates (δe, y, ζ): δė = −k_eδe − Ĵ(t)y, ẏ = Ĵ(t)ᵀδe − "
                  "(σ₀/k_e)(y − ζ), ζ̇ = −(η/k_e)(ζ − y) (η = ∞: ζ = y); its monodromy by integrating it (DOP853, rtol 10⁻¹³; no differences), the error channel's "
                  "exponent by Liouville's formula (the exponents sum to k_e + (σ₀ + η)·9/k_e; η = ∞: k_e). Two references: x_m = 2 sin t, where the first-order "
                  "averaging correction Ā₁ vanishes identically (in one dimension ⟨ĴᵀĴ̇⟩ is the integral of an exact differential, and the commutator term vanishes "
                  "by the symmetry t → π − t), and the asymmetric x_m = 2 sin t + 0.8 sin 2t, where it does not.", ""]
        for r in sorted(RC, key=lambda r: (r['key'] != 'sym', np.isinf(r['eta']))):
            lines.append(f"- {REF_NAME[r['key']]}, η = {eta_s(r['eta'])}: ‖Ā₁‖_max = {r['A1']:.2e}, its self-adjoint part {r['A1_selfadj']:.1e} (Ā₁ is skew); "
                         f"the generator against the central differences of the real loop's field at six times: {r['lin']:.1e} (relative); h_j = "
                         + ", ".join(f"{x:.4e}" for x in r['h']) + ".")
        fdd = []
        for r in RC:
            if r['key'] != 'sym':
                continue
            for row in r['rows']:
                if (r['eta'], row['ke']) in K:
                    fdd.append(np.max(np.abs(np.sort(K[(r['eta'], row['ke'])]['slow']) - row['slow']) / row['slow']))
        if fdd:
            lines += ["", f"Against the finite differences of the first table (x_m = 2 sin t, k_e = 25 and 50, both η) the variational exponents differ by at most "
                      f"{max(fdd):.1e} (relative): the resolution of the differences.", ""]
        lines += ["| reference | η | k_e | slowest mode: (actual − first order)/first order | largest over the modes | largest \|actual − second order\|/first order | "
                  "error channel (Liouville) | k_e − tr H/k_e |", "|---|---|---|---|---|---|---|---|"]
        for r in sorted(RC, key=lambda r: (r['key'] != 'sym', np.isinf(r['eta']))):
            for row in r['rows']:
                d1 = (row['slow'] - row['p1']) / row['p1']; d2 = np.abs(row['slow'] - row['p2']) / row['p1']
                lines.append(f"| {REF_NAME[r['key']]} | {eta_s(r['eta'])} | {row['ke']:g} | {d1[0]:+.4e} | {np.max(np.abs(d1)):.3e} | {np.max(d2):.2e} | "
                             f"{row['err']:.6f} | {row['ke'] - r['trH'] / row['ke']:.6f} |")
        lines += [""]
        for r in sorted(RC, key=lambda r: (r['key'] != 'sym', np.isinf(r['eta']))):
            d1 = [np.max(np.abs(row['slow'] - row['p1']) / row['p1']) for row in r['rows']]
            d2 = [np.max(np.abs(row['slow'] - row['p2']) / row['p1']) for row in r['rows']]
            neg = all(np.all(row['slow'] < row['p1']) for row in r['rows'])
            lines.append(f"- {REF_NAME[r['key']]}, η = {eta_s(r['eta'])}: per doubling of k_e the largest first-order deviation falls by "
                         + ", ".join(f"{a / b:.2f}" for a, b in zip(d1[:-1], d1[1:])) + " (→ 4: O(k_e⁻²) relative) and the second-order residual by "
                         + ", ".join(f"{a / b:.2f}" for a, b in zip(d2[:-1], d2[1:])) + " (→ 16: O(k_e⁻⁴), the odd orders vanishing); every exponent below its first-order "
                         f"prediction at every k_e: {neg} (on this instance: every c_j < 0 here; see the gain-scaled case below)."
                         + (f" Ā₀ is self-adjoint in the metric diag(I, (σ₀/η)I) to {r['A0_metric']:.1e} (in I: {r['A0_plain']:.1e}; |Ā₀| = {r['A0_norm']:.1e})."
                            if 'A0_metric' in r else ""))
        lines += ["", "The second-order coefficients (the slow-time rates are λ_j + c_j/k_e² + O(k_e⁻⁴)), c_j/λ_j in the order of λ_j; the error filter's part "
                  "of the slowest mode's, −κ₁ω₁² with ω₁² = v₁ᵀ⟨Ĵ̇ᵀĴ̇⟩v₁/h₁ = Σ_n n²|Ĵ_nv₁|²/Σ_n|Ĵ_nv₁|² — the mode's mean squared frequency in the regressor, "
                  "independent of η — and κ₁ = h₁(dλ₁/dh)/λ₁ (1 at η = ∞); and the measured (actual·k_e − λ₁)k_e²/λ₁ of the slowest mode at k_e = 100, 200, 400:", ""]
        for r in sorted(RC, key=lambda r: (r['key'] != 'sym', np.isinf(r['eta']))):
            meas = [(row['slow'][0] * row['ke'] - r['lam'][0]) * row['ke'] ** 2 / r['lam'][0] for row in r['rows'] if row['ke'] >= 100]
            if 'omega2' in r:
                om2 = float(r['omega2'][0]); kap = -(r['cf'][0] / r['lam'][0]) / om2
                ftxt = f"of which the filter's {r['cf'][0] / r['lam'][0]:.2f} = −{kap:.3f}ω₁², ω₁ = {np.sqrt(om2):.4f}"
            else:
                ftxt = f"of which the filter's {r['cf'][0] / r['lam'][0]:.2f}"
            lines.append(f"- {REF_NAME[r['key']]}, η = {eta_s(r['eta'])}: c_j/λ_j = " + ", ".join(f"{x:.2f}" for x in r['c'] / r['lam'])
                         + f"; the slowest c₁/λ₁ = {r['c'][0] / r['lam'][0]:.2f}, " + ftxt + "; measured " + ", ".join(f"{x:.2f}" for x in meas) + ".")
        if RG:
            for g in RG:
                if g.get('metric_check'):
                    lines += ["", f"At η = 10σ₀ ({REF_NAME[g['key']]}), where the metric W = diag(I, (σ₀/η)I) is not the identity: |WĀ₀ − (WĀ₀)ᵀ| = {g['A0_metric']:.1e} "
                              f"against |Ā₀ − Ā₀ᵀ| = {g['A0_plain']:.1e} (|Ā₀| = {g['A0_norm']:.1e}): Ā₀ is self-adjoint in W only; |WĀ₁ + (WĀ₁)ᵀ| = "
                              f"{g['A1_skew_metric']:.1e} (|Ā₁| = {g['A1_norm']:.1e})."]
                    continue
                pos = g['c'] > 0
                lines += ["", f"The adaptation gain scaled by {g['gain']:g} ({REF_NAME[g['key']]}, η = {eta_s(g['eta'])}): c_j/λ_j = " + ", ".join(f"{x:+.2f}" for x in g['c'] / g['lam'])
                          + f" — {int(np.sum(pos))} positive (the slaving's −(J̄ᵀJ̄)² and the level repulsion outweigh the filter where |Ĵ|² is large against "
                          "the frequencies). The relative deviations (actual − first order)/first order of those modes at k_e = "
                          + ", ".join(f"{row['ke']:g}" for row in g['rows']) + ": " + "; ".join(", ".join(f"{x:+.3e}" for x in ((row['slow'] - row['p1']) / row['p1'])[pos]) for row in g['rows'])
                          + ", against c_j/(λ_jk_e²) = " + "; ".join(", ".join(f"{x:+.3e}" for x in (g['c'] / g['lam'])[pos] / row['ke'] ** 2) for row in g['rows'])
                          + ": those exponents lie above their first-order prediction, and the slowest below it ("
                          + ", ".join(f"{((row['slow'] - row['p1']) / row['p1'])[0]:+.3e}" for row in g['rows']) + ")."]
    if RD:
        lines += ["", "## D. Lemma D‴.26a: the exact slow subspace (x_m = 2 sin t + 0.8 sin 2t)", "",
                  "Chang's P(t; ε) by its contraction on a 512-point Fourier grid: the iterations to 10⁻¹⁶ and the invariance residual; the truncations "
                  "|P − P₂|_∞/ε³ (P₂ = −K + εK̇ − ε²(K̈ + KKᵀK + KΛ)) and |P − P₍₃₎|_∞/ε⁴ (P₃ = K⃛ + 2K̇KᵀK + KK̇ᵀK + 2KKᵀK̇ + 2K̇Λ); the slow block's slow-time "
                  "rates (from log Φ(ε)/(εT₀)) against the full variational equation's; the rates at ±ε (the slow block at −ε defined by the same equation) "
                  "— max|ν(ε) − ν(−ε)| relative to the largest rate and per mode — against their first-order deviation (relative to the largest rate, per "
                  "mode, and the slowest mode's signed); the same at ±ε with the plant's slope as a mismatch in the error channel; the Magnus parameter "
                  "εT₀max‖A‖₂ (the series converges below π).", ""]
        for rd in RD:
            lines += [f"η = {eta_s(rd['eta'])}{' (the fixed prior at the atoms: K = Ĵ, Λ = σ₀I)' if rd['eta'] == 0 else ''}:", "",
                      "| k_e | iterations | invariance residual | |P − P₂|/ε³ | |P − P₍₃₎|/ε⁴ | slow block vs full | ±ε (largest; per mode) | first-order deviation (largest; per mode; slowest) | ±ε with a mismatch (largest; per mode) | Magnus |",
                      "|---|---|---|---|---|---|---|---|---|---|"]
            for r in rd['rows']:
                lines.append(f"| {r['ke']:g} | {r['it']} | {r['res']:.1e} | {r['trunc']:.4f} | {r['trunc3']:.0f} | {r['vs_full']:.1e} | {r['even']:.1e}; {r['even_mode']:.1e} | "
                             f"{r['first']:.2e}; {r['first_mode']:.2e}; {r['first_slow']:+.3e} | {r['even_mis']:.2e}; {r['even_mis_mode']:.2e} | {r['magnus']:.3f} |")
            lines += ["", "Under white measurement noise (N₀ = 1): the second transformation's decoupling residual; the slow noise b against "
                      "b₀ = −(K̇ + aK)ᵀ (max over the period) and |b − b₀ − εb₁|/ε² (b₁ = −(K̈ + KKᵀK + KΛᵀ + aK̇)ᵀ); the exact periodic covariance at "
                      "t = 0 in the coordinates (v, Y′) — Cov(v)·2ε³ − 1, ‖E[vY′ᵀ] + bᵀ‖/‖b‖, ‖Cov(Y′) − Σ̄‖/‖Σ̄‖ —; the weakest mode's variance (the y block "
                      "in H's eigenbasis) relative to Σ̄: exact, with only the slow noise's attenuation (Σ with ⟨bbᵀ⟩ in place of Q̄), with only the rate's "
                      "correction (Σ with the exact slow generator in place of Ā₀); Cov(Y) − Cov(Y′) against the jitter (ε/2)K(0)ᵀK(0), without and with "
                      "its next term −(ε²/2)(d(KᵀK)/dt + 2aKᵀK)(0); the floor at t = 0, tr(H Cov(y)), against tr(H(Σ̄ + (ε/2)KᵀK)) (the jitter's share); "
                      "and the floor tr(H Cov(y(t))) over the period — its average and range — and the predicted tr(H(Σ̄ + (ε/2)Ĵ(t)ᵀĴ(t)))'s range, both "
                      "relative to Proposition D‴.30(ii)'s formula.", "",
                      "| k_e | decoupling | ‖b − b₀‖/‖b₀‖ | |b − b₀ − εb₁|/ε² | Cov(v)·2ε³ − 1 | ‖E[vY′ᵀ] + bᵀ‖/‖b‖ | ‖Cov(Y′) − Σ̄‖/‖Σ̄‖ | weakest mode: exact; noise only; rate only | jitter: rel. error; with the next term | floor at t = 0: ratio; fast share | floor over the period / formula: average; range | predicted range |",
                      "|---|---|---|---|---|---|---|---|---|---|---|---|"]
            for r in rd['rows']:
                lines.append(f"| {r['ke']:g} | {r['coupling']:.1e} | {r['b_rel']:.3e} | {r['b1']:.0f} | {r['Vv']:+.3e} | {r['C_rel']:.3e} | {r['CYp_rel']:.3e} | "
                             f"{r['weak_exact']:+.3e}; {r['weak_noise']:+.3e}; {r['weak_rate']:+.3e} | {r['jit_rel']:.3f}; {r['jit_next']:.3f} | "
                             f"{r['floor0'] / r['floor0_pred']:.4f}; {r['fast_share']:.4f} | {r['floor_avg']:.5f}; [{r['floor_min']:.4f}, {r['floor_max']:.4f}] | "
                             f"[{r['pred_min']:.4f}, {r['pred_max']:.4f}] |")
            lines += [""]
    lines += ["", "Numbers are produced by `experiments/exp_adaptive_floquet.py`; nothing is edited by hand."]
    with open(os.path.join(OUT, 'adaptive_floquet.md'), 'w', encoding='utf-8') as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    if t0 is not None:
        print(f"[{time.time() - t0:.0f}s]")


def main():
    import multiprocessing as mp
    t0 = time.time()
    raw = os.path.join(OUT, 'adaptive_floquet_raw.npz')
    if '--from-raw' in sys.argv:
        with np.load(raw, allow_pickle=True) as z:
            R = list(z['R'])
            RC = list(z['RC']) if 'RC' in z.files else None
            RG = list(z['RG']) if 'RG' in z.files else None
            RD = list(z['RD']) if 'RD' in z.files else None
        report(R, t0, RC, RG, RD)
        return
    if '--add-D' in sys.argv:                                               # Lemma D‴.26a; the other parts kept
        with np.load(raw, allow_pickle=True) as z:
            save = {k: z[k] for k in z.files}
        with mp.Pool(min(N_WORKERS, 3)) as pool:
            save['RD'] = np.array(pool.map(job_D, [0.0, SIGMA0, np.inf], chunksize=1), dtype=object)
        save_raw(raw, **save)
        report(list(save['R']), t0, list(save['RC']), list(save['RG']), list(save['RD']))
        return
    if '--add-C' in sys.argv:                                               # the derived quantities and the gain-scaled case; R and RC kept
        with np.load(raw, allow_pickle=True) as z:
            R = list(z['R']); RC = list(z['RC'])
        for r in RC:
            r.update(derived_C(r['key'], r['eta']))
        with mp.Pool(min(N_WORKERS, 2)) as pool:
            RG = pool.map(job_C_gain, [('sym', np.inf)], chunksize=1)
        for key in REFS_C:                                                  # the metric at η = 10σ₀, where it is not the identity
            g = derived_C(key, 10 * SIGMA0); g.update(key=key, eta=10 * SIGMA0, metric_check=True); RG.append(g)
        save_raw(raw, R=np.array(R, dtype=object), RC=np.array(RC, dtype=object), RG=np.array(RG, dtype=object))
        report(R, t0, RC, RG)
        return
    if '--only-C' in sys.argv:
        with np.load(raw, allow_pickle=True) as z:
            R = list(z['R'])
        with mp.Pool(min(N_WORKERS, 4)) as pool:
            RC = pool.map(job_C, [(key, eta) for key in REFS_C for eta in ETAS], chunksize=1)
        save_raw(raw, R=np.array(R, dtype=object), RC=np.array(RC, dtype=object))
        report(R, t0, RC)
        return
    jobs = [(mult, eta, ke) for mult in MULTS for eta in ETAS for ke in KES]
    jobs.sort(key=lambda j: -j[0])                                          # the large systems first
    with mp.Pool(N_WORKERS) as pool:
        R = pool.map(job, jobs, chunksize=1)
        RC = pool.map(job_C, [(key, eta) for key in REFS_C for eta in ETAS], chunksize=1)
    save_raw(raw, R=np.array(R, dtype=object), RC=np.array(RC, dtype=object))
    report(R, t0, RC)


if __name__ == '__main__':
    main()
