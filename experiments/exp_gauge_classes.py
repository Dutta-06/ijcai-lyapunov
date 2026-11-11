"""
Proposition D‴.42: several particles per class — the class masses carry the gauge's laws, the split within a class is frozen
=========================================================================================================================

Proposition D‴.36's instance with n particles per class (m = 3n), the particles of a class at the class's exact-fit position (the class's
amplitude A_k/W_k, W_k = B_k/m, B_k = Σ_{i∈k}β_i), the class's mass split among them in the fractions s_i = β_i/B_k. Coordinates on the
manifold of exact fits M_m: the class masses (W_1, W_2) and the fractions (the first n − 1 of each class).

A. The determinant: log det ∇²_ΘJ at m = 6 and 9, at three class-mass points and three splits, minus the one-particle-per-class value at
   the same class masses — a constant (the cohesion's within-class Laplacian): independent of the class masses and of the split.
B. Katzenberger's limit at m = 6 and 9 (γ = 10⁻⁵, 10⁻²): the class masses' drift and diffusion against m = 3's at the same class masses
   (the ratio 3/m), the fractions' drift and diffusion.
C. Proposition D‴.38's correction per class (the masses form at the class masses and contributions, `fixman_classes`, given to each of
   the class's particles) at κ_F = 1 and γ = 10⁻⁵: the corrected class drift against m = 3's; the fractions' drift.
D. Measurement noise (Proposition D‴.41) at m = 6: the Itô drift, the Wong–Zakai term and their sum at equal and unequal splits, per
   component — the class masses' Stratonovich residual over their own Itô drift and the fractions' over theirs; the class masses' Itô
   drift over m = 3's at the same class masses (1 at the equal split: there co-located particles of equal mass get identical drifts and
   identical noise, an invariant subspace on which the loop is the m = 3 loop); the fractions' Itô drift against the replicator flow of the
   full trace's part (v/4)Σ_k⟨‖∇_ϑΦ̃_k‖²_G⟩_ρΣ_{i∈k}β_i² and of its amplitudes' part alone.
E. The full loop at m = 6 (γ = 10⁻², τ = 10⁻⁴, SEEDS_E, S_E slow units) from equal fractions and from (0.2, 0.7, 0.5): the class masses'
   spread against Proposition D‴.37's Part H (m = 3, the same γ, τ) at half the time (B's ratio 3/m = 1/2), the fractions' change
   against B's drift at the start.
F. Anisotropic dither (Proposition D‴.40's masks) at m = 6: Katzenberger's drift at γ = 10⁻³ … 10⁻⁶ with the noise on the amplitudes only
   (1, 0, 0), on the centres and widths only (0, 1, 1), and isotropic (1, 1, 1), at the equal split and at (0.2, 0.7, 0.5): the class drift
   per γτ against the per-class law (m = 3's with the same mask, times 3/m) and the fractions' drift per γτ — the frozen split and the
   per-class law need the isotropic dither (the positions' law at fixed masses Gibbs); an anisotropic one correlates a class's barycentre
   with its internal differences at unequal splits, and the fractions drift at O(γτ).

Output: results/chaos/gauge_classes.md and gauge_classes_raw.npz (GC_WORKERS; `--from-raw` regenerates the md; `--only A,B,…,F`).
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
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
from src.core import mass_transfer as MT
from src.core import gauge_drift as GDR
from src.core import label_noise as LN
from src.core.runcache import cached
import exp_cohesion_leak as E
import exp_gauge_dither as GD
import exp_gauge_drift as X

OUT = os.environ.get('GC_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
N_WORKERS = int(os.environ.get('GC_WORKERS', 4))
TAU, A_MAX, DT, REC = GD.TAU, 4.0, GD.DT, GD.REC
WC_A = ((1 / 3, 1 / 3, 1 / 3), (0.30, 0.35, 0.35), (0.40, 0.25, 0.35))
SPLITS_A = {2: ((0.5, 0.5, 0.5), (0.2, 0.7, 0.5), (0.05, 0.95, 0.3)), 3: ((1 / 3,) * 6, (0.2, 0.3, 0.6, 0.1, 0.4, 0.4), (0.05, 0.05, 0.9, 0.05, 0.2, 0.2))}
POINTS_B = (((1 / 3, 1 / 3), (0.5, 0.5, 0.5)), ((1 / 3, 1 / 3), (0.2, 0.7, 0.5)), ((0.30, 0.35), (0.3, 0.6, 0.8)))
POINTS_B3 = (((1 / 3, 1 / 3), (0.30, 0.30, 0.35, 0.35, 0.32, 0.36)), ((0.30, 0.35), (0.2, 0.3, 0.6, 0.1, 0.4, 0.4)),
             ((1 / 3, 1 / 3), (1 / 3,) * 6), ((0.30, 0.35), (1 / 3,) * 6))   # at n = 3's equal split the cohesion's within-class modes are
# exactly degenerate: `limit_drift` then leaves A's eigenbasis for its Schur method
MASKS_F, GAMMAS_F, POINTS_F = ((1, 0, 0), (0, 1, 1), (1, 1, 1)), (1e-3, 1e-4, 1e-5, 1e-6), (((1 / 3, 1 / 3), (0.5, 0.5, 0.5)), ((1 / 3, 1 / 3), (0.2, 0.7, 0.5)))
GAMMAS_B = (1e-5, 1e-2)
GAMMA_E, S_E, SEEDS_E = 1e-2, 4e5, tuple(range(1000, 1008))
STARTS_E = ((0.5, 0.5, 0.5), (0.2, 0.7, 0.5))


def split_beta(Wc, s, n):
    """The masses β (m = 3n) for the class masses Wc and the fractions s (n − 1 per class, the last one the complement)."""
    m = 3 * n
    s = np.asarray(s, float).reshape(3, n - 1)
    fr = np.c_[s, 1 - s.sum(axis=1)]
    return (m * np.asarray(Wc, float)[:, None] * fr).ravel()


def exact_m(net, Wc, s, n):
    return np.repeat(E.exact_rep(net, np.asarray(Wc, float)), n, axis=0), split_beta(Wc, s, n)


def log_det(loop, Th, beta, h=1e-4):
    x0 = Th.ravel(); N = x0.size
    J = lambda t: loop.J(t.reshape(Th.shape), beta)
    H = np.zeros((N, N)); I = np.eye(N) * h
    for i in range(N):
        for j in range(i, N):
            H[i, j] = H[j, i] = (J(x0 + I[i] + I[j]) - J(x0 + I[i] - I[j]) - J(x0 - I[i] + I[j]) + J(x0 - I[i] - I[j])) / (4 * h * h)
    sg, ld = np.linalg.slogdet(H)
    assert sg > 0
    return ld


def a_job(args):
    n, Wc = args
    net, base = E.setup()
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2)
    l1 = log_det(loop, E.exact_rep(net, np.asarray(Wc)), 3 * np.asarray(Wc))
    out = []
    for s in SPLITS_A[n]:
        Th, beta = exact_m(net, Wc, s, n)
        out.append(dict(split=s, diff=log_det(loop, Th, beta) - l1, resid=float(np.abs(loop.residual(Th, beta)).max())))
    return dict(n=n, Wc=Wc, l1=l1, rows=out)


def instance(n, gamma, **kw):
    """The flow in (Θ, log(β_j/β_m)), the chart of M_m by (W_1, W_2, the fractions), the coordinates, the dither's covariance."""
    m = 3 * n
    net, base = E.setup()
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=gamma, **kw)

    def beta_of(z):
        e = np.exp(np.r_[z, 0.0]); return m * e / e.sum()

    def F(x):
        T = x[:3 * m].reshape(m, 3); b = beta_of(x[3 * m:]); dT, db = loop.drift(T, b); dl = db / b
        return np.concatenate([dT.ravel(), dl[:m - 1] - dl[m - 1]])

    def chart(w):
        Wc = np.r_[w[:2], 1 - w[:2].sum()]
        Th, beta = exact_m(net, Wc, w[2:], n)
        return np.concatenate([Th.ravel(), np.log(beta[:m - 1] / beta[m - 1])])

    def coords(x):
        b = beta_of(x[3 * m:]).reshape(3, n); Bc = b.sum(axis=1)
        return np.r_[Bc[:2] / m, (b[:, :n - 1] / Bc[:, None]).ravel()]
    S = np.zeros((4 * m - 1, 4 * m - 1)); S[np.arange(3 * m), np.arange(3 * m)] = 2 * TAU * np.tile(E.G, m)
    return F, chart, coords, S, loop, beta_of


def b_job(args):
    n, gamma, Wc, s = args
    F, chart, coords, S, loop, beta_of = instance(n, gamma)
    w = np.r_[Wc, s]
    d, D, lam = GDR.limit_drift(F, chart, coords, w, S, h=X.H_ID, h2=X.H2_ID)
    F3, chart3, coords3, S3 = X.instance(gamma)
    d3, D3, _ = GDR.limit_drift(F3, chart3, coords3, np.asarray(Wc, float), S3, h=X.H_ID, h2=X.H2_ID)
    return dict(n=n, gamma=gamma, Wc=Wc, s=s, cls=d[:2] / (gamma * TAU), cls3=d3 / (gamma * TAU), frac=d[2:] / (gamma * TAU),
                Dcls=np.diag(D)[:2] / (gamma ** 2 * TAU), Dcls3=np.diag(D3) / (gamma ** 2 * TAU), Dfrac=float(np.abs(D[2:, 2:]).max()) / (gamma ** 2 * TAU),
                slow=float(np.abs(lam).min()))


def c_job(args):
    n, gamma, Wc, s = args
    m = 3 * n
    lab = np.repeat(np.arange(3), n)
    kw = dict(fixman=1.0, fixman_tau=TAU, fixman_A=E.A_BUMP.copy(), fixman_classes=lab)
    F, chart, coords, S, loop, beta_of = instance(n, gamma)
    Fc = instance(n, gamma, **kw)[0]
    w = np.r_[Wc, s]
    d0 = GDR.limit_drift(F, chart, coords, w, S, h=X.H_ID, h2=X.H2_ID)[0]
    dc = GDR.limit_drift(F, chart, coords, w, S, h=X.H_ID, h2=X.H2_ID, F1=lambda x: Fc(x) - F(x))[0]
    F3, chart3, coords3, S3 = X.instance(gamma)
    loop3c = MT.ReplicatorCohesionLoop(E.setup()[1], 1.0, 'compact', gamma=gamma, fixman=1.0, fixman_tau=TAU, fixman_A=E.A_BUMP.copy())

    def F3c(x):
        T = x[:9].reshape(3, 3); e = np.exp(np.r_[x[9:], 0.0]); b = 3 * e / e.sum(); dT, db = loop3c.drift(T, b); dl = db / b
        return np.concatenate([dT.ravel(), dl[:2] - dl[2]])
    d3 = GDR.limit_drift(F3, chart3, coords3, np.asarray(Wc, float), S3, h=X.H_ID, h2=X.H2_ID)[0]
    d3c = GDR.limit_drift(F3, chart3, coords3, np.asarray(Wc, float), S3, h=X.H_ID, h2=X.H2_ID, F1=lambda x: F3c(x) - F3(x))[0]
    return dict(n=n, gamma=gamma, Wc=Wc, s=s, cls0=d0[:2] / (gamma * TAU), clsc=dc[:2] / (gamma * TAU), frac0=d0[2:] / (gamma * TAU),
                fracc=dc[2:] / (gamma * TAU), cls3=d3 / (gamma * TAU), cls3c=d3c / (gamma * TAU))


def d_job(args):
    gamma, Wc, s = args
    n = 2
    F, chart, coords, S, loop, beta_of = instance(n, gamma)
    w = np.r_[Wc, s]; x0 = chart(w)
    U0 = LN.flat_fields(loop, x0, beta_of, m=3 * n)
    Sig = U0.T @ U0 / U0.shape[0]
    ito, D, lam = GDR.limit_drift(F, chart, coords, w, Sig, h=X.H_ID, h2=X.H2_ID)
    A = GDR._jac(F, x0, X.H_ID); Jg = GDR._jac(chart, w, X.H_ID); k = len(w)
    ev, V = np.linalg.eig(A); o = np.argsort(np.abs(ev)); V = np.c_[Jg, V[:, o][:, k:].real]
    P = V[:, :k] @ np.linalg.inv(V)[:k, :]
    Jw = GDR._jac(coords, x0, X.H_ID)
    wz = Jw @ P @ (0.5 * LN.wong_zakai(loop, x0, beta_of, m=3 * n))
    F3, chart3, coords3, S3 = X.instance(gamma)                               # m = 3 at the same class masses
    W3 = np.asarray(Wc, float); x3 = chart3(W3)
    net, base = E.setup()
    loop3 = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=gamma)
    U3 = LN.flat_fields(loop3, x3, lambda z: (lambda e: 3 * e / e.sum())(np.exp(np.r_[z, 0.0])))
    ito3 = GDR.limit_drift(F3, chart3, coords3, W3, U3.T @ U3 / U3.shape[0], h=X.H_ID, h2=X.H2_ID)[0]
    # the fractions' Itô drift from the replicator flow of (v/4)Σ_kC_kΣ_{i∈k}β_i²: ds_i/ds = −γ(v/2)C_kB_ks_i(s_i − Σ_{j∈k}s_j²)
    Wf = np.r_[Wc, 1 - sum(Wc)]
    UT = LN.fields(loop3, E.exact_rep(net, Wf), np.ones(3))[0]                # β = 1: G∇_ϑΦ̃ per class
    C_full = np.mean(np.sum(UT ** 2 / E.G[None, None, :], axis=2), axis=0)    # ⟨‖∇_ϑΦ̃_k‖²_G⟩_ρ
    C_amp = np.mean(UT[:, :, 0] ** 2 / E.G[0], axis=0)                        # its amplitudes' part
    fr = np.asarray(s, float)
    pred = lambda Ck: np.array([-gamma * 0.5 * Ck[k] * 6 * Wf[k] * fr[k] * (fr[k] - fr[k] ** 2 - (1 - fr[k]) ** 2) for k in range(3)])
    return dict(gamma=gamma, Wc=Wc, s=s, ito=ito, wz=wz, D=float(np.abs(D).max()), ito3=ito3, frac_full=pred(C_full), frac_amp=pred(C_amp))


H2_F = (X.H2_ID, 2 * X.H2_ID, 4 * X.H2_ID)                                    # Part F: the second derivatives' steps, extrapolated


def f_job(args):
    """Part F: Katzenberger's drift at m = 6 under the dither masked per coordinate (α, centre, width), and m = 3's with the same mask —
    Richardson-extrapolated in the second derivatives' step h2 (the central differences' error is O(h2²), and the drifts here are small
    enough for it to matter: referee round 30), from (h2, 2h2); the same from (2h2, 4h2) gives the resolution."""
    mask, gamma, Wc, s = args
    n, m = 2, 6
    F, chart, coords, S, loop, beta_of = instance(n, gamma)
    Sm = S.copy(); Sm[np.arange(3 * m), np.arange(3 * m)] *= np.tile(np.asarray(mask, float), m)
    F3, chart3, coords3, S3 = X.instance(gamma)
    S3m = S3.copy(); S3m[np.arange(9), np.arange(9)] *= np.tile(np.asarray(mask, float), 3)
    d = np.array([GDR.limit_drift(F, chart, coords, np.r_[Wc, s], Sm, h=X.H_ID, h2=h2)[0] for h2 in H2_F]) / (gamma * TAU)
    d3 = np.array([GDR.limit_drift(F3, chart3, coords3, np.asarray(Wc, float), S3m, h=X.H_ID, h2=h2)[0] for h2 in H2_F]) * (3 / m) / (gamma * TAU)
    rich = lambda x, i: (4 * x[i] - x[i + 1]) / 3
    dR, dR2, d3R, d3R2 = rich(d, 0), rich(d, 1), rich(d3, 0), rich(d3, 1)
    return dict(mask=mask, gamma=gamma, Wc=Wc, s=s, cls=dR[:2], law=d3R, frac=dR[2:], cls_h2=d[0, :2], law_h2=d3[0],
                res_cls=np.abs(dR - dR2)[:2], res_law=np.abs(d3R - d3R2), res_frac=np.abs(dR - dR2)[2:])


@cached
def e_series(args):
    """The full loop at m = 6 from the fractions s0: records (s, the class masses W_1..3, the fractions s_1..3, the largest |α|, the fit error)."""
    s0, seed, gamma, S = args
    t0 = time.time()
    net, base = E.setup()
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=gamma)
    Th, beta = exact_m(net, np.ones(3) / 3, s0, 2)
    rec = []

    def cb(s_, T, b):
        bb = b.reshape(3, 2); Bc = bb.sum(axis=1)
        rec.append(np.r_[s_, Bc / 6, bb[:, 0] / Bc, np.abs(T[:, 0]).max(), loop.fit_error(T, b)])
    loop.run_dithered(Th, beta, S, DT, TAU, np.random.default_rng(seed), rec_every=REC, callback=cb)
    print(f"[{time.time() - t0:.0f}s] E s0 {s0} seed {seed}", flush=True)
    return dict(s0=s0, seed=seed, rec=np.array(rec))


def summarise_E(runs, H_rec):
    out = []
    for s0 in STARTS_E:
        R = np.array([r['rec'] for r in runs if tuple(r['s0']) == tuple(s0)])
        s = R[0, :, 0]
        i_half = int(np.argmin(np.abs(H_rec[0, :, 0] - S_E / 2)))
        W6 = R[:, -1, 1:4]; W3 = np.asarray(H_rec, float)[:, i_half, 4:7]
        fr = R[:, :, 4:7]
        out.append(dict(s0=s0, n=len(R), S=s[-1], std6=W6.std(axis=0, ddof=1), mean6=W6.mean(axis=0), std3=W3.std(axis=0, ddof=1),
                        mean3=W3.mean(axis=0), s3=float(H_rec[0, i_half, 0]), dfrac=(fr[:, -1] - fr[:, 0]).mean(axis=0),
                        dfrac_se=(fr[:, -1] - fr[:, 0]).std(axis=0, ddof=1) / np.sqrt(len(R)), frac_std=fr.std(axis=1).mean(axis=0),
                        almax=float(R[:, :, 7].max()), fit=float(R[:, R.shape[1] // 2:, 8].mean())))
    return out


def report(A=None, B=None, C=None, D=None, Ee=None, Ff=None):
    v = lambda x, p=4: ', '.join(f'{y:+.{p}f}' for y in x)
    L_ = ["# Proposition D‴.42: several particles per class (`exp_gauge_classes.py`)", "",
          f"Proposition D‴.36's instance with n particles per class (m = 3n) at the class's exact-fit position, τ = {TAU:g}; the class masses W_k "
          "and the fractions s_i = β_i/B_k within each class.", ""]
    if A is not None:
        L_ += ["## A. The determinant: log det ∇²_ΘJ(m) − log det ∇²_ΘJ(3) at the same class masses", "",
               "| n | the class masses | the splits (fractions per class) → the difference | the exact fit's residual |", "|---|---|---|---|"]
        for r in A:
            L_.append(f"| {r['n']} | {', '.join(f'{x:.3f}' for x in r['Wc'])} | "
                      + '; '.join(f"({', '.join(f'{x:.2f}' for x in q['split'])}) → {q['diff']:.6f}" for q in r['rows'])
                      + f" | {max(q['resid'] for q in r['rows']):.1e} |")
        for n in sorted({r['n'] for r in A}):
            d = [q['diff'] for r in A if r['n'] == n for q in r['rows']]
            L_.append(f"\nn = {n}: the difference {np.mean(d):.6f}, its range over the {len(d)} points {np.ptp(d):.1e}.")
        L_.append("")
    if B is not None:
        L_ += ["## B. Katzenberger's limit: the class masses against one particle per class, and the fractions", "",
               "| n | γ | W_1, W_2 | the fractions | class drift / γτ | m = 3's | ratio | class diffusion / γ²τ (diagonal) | m = 3's | ratio | "
               "the fractions' drift / γτ | their diffusion / γ²τ (largest) | the slowest normal rate |",
               "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
        for r in B:
            L_.append(f"| {r['n']} | {r['gamma']:g} | {', '.join(f'{x:.3f}' for x in r['Wc'])} | {', '.join(f'{x:.2f}' for x in r['s'])} | "
                      f"{v(r['cls'])} | {v(r['cls3'])} | {', '.join(f'{a / b:.4f}' for a, b in zip(r['cls'], r['cls3']))} | "
                      f"{', '.join(f'{x:.4f}' for x in r['Dcls'])} | {', '.join(f'{x:.4f}' for x in r['Dcls3'])} | "
                      f"{', '.join(f'{a / b:.4f}' for a, b in zip(r['Dcls'], r['Dcls3']))} | {', '.join(f'{x:+.2e}' for x in r['frac'])} | "
                      f"{r['Dfrac']:.1e} | {r['slow']:.3e} |")
        L_.append("")
    if C is not None:
        L_ += ["## C. The correction per class (κ_F = 1, the masses form at the class masses)", "",
               "| n | γ | W_1, W_2 | the fractions | class drift / γτ, uncorrected; corrected | m = 3's, uncorrected; corrected | "
               "the fractions' drift / γτ, uncorrected; corrected |", "|---|---|---|---|---|---|---|"]
        for r in C:
            L_.append(f"| {r['n']} | {r['gamma']:g} | {', '.join(f'{x:.3f}' for x in r['Wc'])} | {', '.join(f'{x:.2f}' for x in r['s'])} | "
                      f"{v(r['cls0'])}; {v(r['clsc'], 6)} | {v(r['cls3'])}; {v(r['cls3c'], 6)} | {', '.join(f'{x:+.2e}' for x in r['frac0'])}; "
                      f"{', '.join(f'{x:+.2e}' for x in r['fracc'])} |")
        L_.append("")
    if D is not None:
        per = 'ito3' in D[0]
        L_ += ["## D. Measurement noise at m = 6 (per unit v): the Itô drift, the Wong–Zakai term, their sum (W_1, W_2, s_1, s_2, s_3)", "",
               "| γ | W_1, W_2 | the fractions | the Itô drift | the Wong–Zakai term | their sum | its size over the Itô drift's | " +
               ("the class masses' sum over their Itô drift (W_1, W_2) | the fractions' sum over their Itô drift (largest) | the class Itô drift over "
                "m = 3's (W_1, W_2) | " if per else "") + "the diffusion along M |",
               "|---|---|---|---|---|---|---|" + ("---|---|---|" if per else "") + "---|"]
        for r in D:
            tot = r['ito'] + r['wz']
            pc = ""
            if per:
                fi = np.abs(r['ito'][2:]).max()
                pc = (f"{', '.join(f'{a / b:+.2e}' for a, b in zip(tot[:2], r['ito'][:2]))} | "
                      + (f"{np.abs(tot[2:]).max() / fi:.2e} | " if fi > 1e-6 * np.abs(r['ito'][:2]).max() else "— (no Itô drift) | ")
                      + f"{', '.join(f'{a / b:.4f}' for a, b in zip(r['ito'][:2], r['ito3']))} | ")
            L_.append(f"| {r['gamma']:g} | {', '.join(f'{x:.3f}' for x in r['Wc'])} | {', '.join(f'{x:.2f}' for x in r['s'])} | "
                      f"{', '.join(f'{x:+.3e}' for x in r['ito'])} | {', '.join(f'{x:+.3e}' for x in r['wz'])} | {', '.join(f'{x:+.2e}' for x in tot)} | "
                      f"{np.abs(tot).max() / np.abs(r['ito']).max():.2e} | " + pc + f"{r['D']:.1e} |")
        L_.append("")
        if per:
            L_ += ["The fractions' Itô drift against the replicator flow of (v/4)Σ_kC_kΣ_{i∈k}β_i² — C_k = ⟨‖∇_ϑΦ̃_k‖²_G⟩_ρ (the full trace's part: the "
                   "fractions' tangent moves no position, so no amplitude correction enters; this flow is minus the Wong–Zakai term's fractions, and "
                   "the Itô drift differs from it by the cohesion's Stratonovich residual) and C_k = G_α⟨(∂_αΦ̃_k)²⟩_ρ (its amplitudes' part alone):", "",
                   "| γ | W_1, W_2 | the fractions | the fractions' Itô drift | the full trace's flow | the amplitudes' part's flow |", "|---|---|---|---|---|---|"]
            for r in D:
                if np.abs(r['ito'][2:]).max() > 1e-6 * np.abs(r['ito'][:2]).max():
                    L_.append(f"| {r['gamma']:g} | {', '.join(f'{x:.3f}' for x in r['Wc'])} | {', '.join(f'{x:.2f}' for x in r['s'])} | "
                              f"{', '.join(f'{x:+.3e}' for x in r['ito'][2:])} | {', '.join(f'{x:+.3e}' for x in r['frac_full'])} | "
                              f"{', '.join(f'{x:+.3e}' for x in r['frac_amp'])} |")
            L_.append("")
    if Ee is not None:
        L_ += [f"## E. The full loop at m = 6 (γ = {GAMMA_E:g}, τ = {TAU:g}, {Ee[0]['n']} seeds, {S_E:g} slow units)", "",
               "| the start's fractions | the class masses at the end: mean; std over the seeds | m = 3 at half the time (Proposition D‴.37's Part H): "
               "mean; std | the fractions' change (± s.e.) | their within-run std | the largest |α| | the fit error |", "|---|---|---|---|---|---|---|"]
        for r in Ee:
            L_.append(f"| {', '.join(f'{x:.2f}' for x in r['s0'])} | {', '.join(f'{x:.4f}' for x in r['mean6'])}; {', '.join(f'{x:.4f}' for x in r['std6'])} | "
                      f"{', '.join(f'{x:.4f}' for x in r['mean3'])}; {', '.join(f'{x:.4f}' for x in r['std3'])} (at {r['s3']:g}) | "
                      f"{', '.join(f'{a:+.4f} ± {b:.4f}' for a, b in zip(r['dfrac'], r['dfrac_se']))} | {', '.join(f'{x:.4f}' for x in r['frac_std'])} | "
                      f"{r['almax']:.2f} | {E.fmt(r['fit'])} |")
        L_.append("")
    if Ff is not None:
        rich = 'res_cls' in Ff[0]
        L_ += ["## F. Anisotropic dither at m = 6: the class drift against the per-class law, the fractions' drift (per γτ)", "",
               ("Richardson-extrapolated in the second derivatives' step h2 from (h2, 2h2), h2 = " + f"{H2_F[0]:g}" + "; the resolution: the same "
                "from (2h2, 4h2), its difference; a ratio is printed only where both drifts exceed ten times their resolution. At h2 alone the "
                "central differences' O(h2²) error is a few percent of these drifts (referee round 30)." if rich else ""), "",
               "| the mask (α, centre, width) | the fractions | γ | the class drift (W_1, W_2) | the per-class law (m = 3's × 3/m) | their ratio | "
               "the resolution (class drift; law) | the ratio at h2 alone | the fractions' drift (± resolution) |", "|---|---|---|---|---|---|---|---|---|"]
        for r in Ff:
            if rich:
                ok = (np.abs(r['cls']) > 10 * r['res_cls']) & (np.abs(r['law']) > 10 * r['res_law'])
                rat = ', '.join(f'{a / b:.4f}' if o else 'unresolved' for a, b, o in zip(r['cls'], r['law'], ok))
                res_ = f"{', '.join(f'{x:.1e}' for x in r['res_cls'])}; {', '.join(f'{x:.1e}' for x in r['res_law'])}"
                rat_h2 = ', '.join(f'{a / b:.4f}' for a, b in zip(r['cls_h2'], r['law_h2']))
                fr = ', '.join(f'{x:+.3e} ± {e:.0e}' for x, e in zip(r['frac'], r['res_frac']))
            else:
                rat = ', '.join(f'{a / b:.4f}' for a, b in zip(r['cls'], r['law'])); res_ = rat_h2 = '—'
                fr = ', '.join(f'{x:+.3e}' for x in r['frac'])
            ve = lambda x: ', '.join(f'{y:+.4e}' for y in x)
            L_.append(f"| {tuple(r['mask'])} | {', '.join(f'{x:.2f}' for x in r['s'])} | {r['gamma']:g} | {ve(r['cls'])} | {ve(r['law'])} | "
                      f"{rat} | {res_} | {rat_h2} | {fr} |")
        L_.append("")
    L_.append("Numbers are produced by `experiments/exp_gauge_classes.py`; nothing is edited by hand.")
    return '\n'.join(L_) + '\n'


def main():
    raw = os.path.join(OUT, 'gauge_classes_raw.npz')
    md = os.path.join(OUT, 'gauge_classes.md')
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
        if 'E' in parts:
            h['E'] = pool.map_async(e_series, [(s0, sd, GAMMA_E, S_E) for s0 in STARTS_E for sd in SEEDS_E], chunksize=1)
        if 'A' in parts:
            h['A'] = pool.map_async(a_job, [(n, Wc) for n in (2, 3) for Wc in WC_A], chunksize=1)
        if 'B' in parts:
            h['B'] = pool.map_async(b_job, [(2, g, Wc, s) for g in GAMMAS_B for Wc, s in POINTS_B]
                                    + [(3, g, Wc, s) for g in GAMMAS_B for Wc, s in POINTS_B3], chunksize=1)
        if 'C' in parts:
            h['C'] = pool.map_async(c_job, [(2, 1e-5, Wc, s) for Wc, s in POINTS_B], chunksize=1)
        if 'D' in parts:
            h['D'] = pool.map_async(d_job, [(g, Wc, s) for g in GAMMAS_B for Wc, s in POINTS_B], chunksize=1)
        if 'F' in parts:
            h['F'] = pool.map_async(f_job, [(mk, g, Wc, s) for mk in MASKS_F for Wc, s in POINTS_F for g in GAMMAS_F], chunksize=1)
        for k in 'ABCDF':
            if k in h:
                res[k] = h[k].get(); E.save(raw, **res); log(k)
        if 'E' in h:
            runs = h['E'].get()
            res['E_runs'] = runs
            z = E.load(os.path.join(OUT, 'gauge_drift_raw.npz'))
            res['E'] = summarise_E(runs, np.asarray(z['H_rec'], float))
            E.save(raw, **res); log('E')
    open(md, 'w').write(report(*[res.get(k) for k in 'ABCDEF']))


if __name__ == '__main__':
    main()
