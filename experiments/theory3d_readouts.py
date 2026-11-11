"""
Round theory3, direction D: which readouts select - theorems for the readouts of depth p (exact scalar unit), their numerical checks and the mismatch margin of the 200 benchmark targets.

The note (outside the repository) proves, for a unit of depth p trained by gradient flow from the start s^2 = eps b: its coefficient is an explicit function of the dual eta = int g
(depth 1: w = ell eta; depth 2: 2 eps b sinh(4 ell eta); depth p >= 3: s^p[(1-eta/lam)^-q - (1+eta/lam)^-q], q = p/(p-2)), hence (a) depth p >= 3: |w| <= (eps b)^(p/2) mu^(-p/(p-2)) at
the margin mu of the wall and the unit is formed iff its dual is within d_theta ~ theta^(-(p-2)/p)/(ell p (p-2)) of the wall, i.e. at a margin below ~ theta^(-(p-2)/p)(eps b)^((p-2)/2);
(b) depth 2: exact formation margin ln(2/(theta + sqrt(theta^2 + 4 (eps b)^2)))/ln(1/(eps b)); (c) depth 1: w = ell eta, no wall (fixed centres: the ridge-like gradient flow). This script
verifies them against the exact ODE (A), in ball arithmetic (B), against the repository's network flows with moving centres (C), the depth-1 flow (D), the sliding unit's first integral and
its competitor-overtaking margin (E), and applies them to the 200 targets of `exp_parametrization.py` Part C (F) and to the per-target scan of `theory3d_scan.py` (G); H extrapolates the
overtaking model to the 1034 one-dimensional targets of the benchmark and explains the observed exponent of the mismatch margin; L the fixed-centre limit (the duals of D1, D2, D3 and the paper's unit against the weighted-l1 inverse scale space as eps -> 0); K is a one-parameter family (one feature between two kernels) where the
critical eps of each member is found by bisection: the margin of the tie against the critical eps, at fixed geometry.

Usage: python experiments/theory3d_readouts.py [--run [--only A,B,C,D,E,L,K,P] ] [--analysis] [--from-raw]
   --run       parts A-E and P (the sliding paths of the 1034 targets, ~10 min) and then the analyses F, G, H;  --analysis: F, G, H only (from the stored P, the stored runs, the scan);
   --from-raw  rewrites results/synthetic/theory3d.md from the stored raws (theory3d_raw.npz, theory3d_scan_raw.npz).
Environment: T3D_OUT (output directory; default results/synthetic), T3D_PARAM_RAW (the stored runs of exp_parametrization), T3D_CERT_RAW (the stored targets of exp_transit_certificate).
"""
import os
import sys
import time
import warnings

os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
from scipy.integrate import solve_ivp

HERE = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..'))
import theory3d_unit as U            # noqa: E402
import theory3d_scan as SC           # noqa: E402
import exp_parametrization as P      # noqa: E402
from src.core import factorized_rbf as FR      # noqa: E402
from src.core import rich_transit_nd as RN     # noqa: E402

OUT = os.environ.get('T3D_OUT', os.path.join(HERE, '..', 'results', 'synthetic'))
RAW, MD = os.path.join(OUT, 'theory3d_raw.npz'), os.path.join(OUT, 'theory3d.md')
SCAN_RAW = os.path.join(OUT, 'theory3d_scan_raw.npz')
PARAM_RAW = os.environ.get('T3D_PARAM_RAW', os.path.join(HERE, '..', 'results', 'chaos', 'parametrization_raw.npz'))
CERT_RAW = os.environ.get('T3D_CERT_RAW', os.path.join(HERE, '..', 'results', 'chaos', 'transit_certificate_raw.npz'))
PARTS = ('A', 'B', 'C', 'D', 'E', 'L', 'K', 'P', 'F', 'G', 'H')
MU_C = P.LR_C                                  # the centres' rate G m = 12 (centres in the units of the samples)
A_RO = {'D3': U.lead_const('D3'), 'P0': U.lead_const('P0')}
EPS_H = tuple(10.0 ** (-k / 2.0) for k in range(6, 25))      # 1e-3 ... 1e-12
BOOT_N, BOOT_SIZE = 2000, 200                 # H: bootstrap samples of the population, of the size of Part C's sample
SCORE_FLIP = 1.0                               # overtaking margin / the competitor's margin at the entry >= 1: the model predicts a flip


def save(path, res):
    np.savez_compressed(path, **{k: np.array([v], dtype=object) for k, v in res.items() if v is not None})


def load(path):
    if not os.path.exists(path):
        return {}
    z = np.load(path, allow_pickle=True)
    return {k: z[k][0] for k in z.files}


def formed_theta(As):
    """the formed threshold of a target: min(0.01 max|A|, 0.1 min|A|)"""
    return FR.formed_threshold(As)


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# A. the exact scalar unit against the coordinate ODE, for a prescribed drive
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
def part_A():
    rows = []
    for p in (2, 3, 4, 5, 6):
        eps, b, ell = 1e-3, 1.0, 0.7
        lam = U.lam_2(eps, b, ell) if p == 2 else U.lam_p(p, eps, b, ell)
        om, T = 0.7, 40.0
        tt = np.linspace(0, T, 4001)
        e_unit = (1 - np.cos(om * tt)) / om + 0.3 * tt          # int_0^t (sin(om s) + 0.3) ds
        Amp = 0.9 * lam / np.abs(e_unit).max()
        g = lambda t: Amp * (np.sin(om * t) + 0.3)
        s = np.sqrt(eps * b)

        def rhs(t, y):
            u, v, e = y
            gg = g(t)
            return [ell * p * u ** (p - 1) * gg, -ell * p * v ** (p - 1) * gg, gg]
        sol = solve_ivp(rhs, (0, T), [s, s, 0.0], t_eval=np.linspace(0, T, 401), method='DOP853', rtol=1e-12, atol=1e-18)
        u, v, e = sol.y
        w_ode = u ** p - v ** p
        w_cf = U.W_2(e, eps, b, ell) if p == 2 else U.W_p(e, p, eps, b, ell)
        e_exact = Amp * ((1 - np.cos(om * sol.t)) / om + 0.3 * sol.t)
        rel = np.abs(w_ode - w_cf) / np.maximum(np.abs(w_cf), 1e-300)
        sel = np.abs(w_cf) > 1e-9 * np.abs(w_cf).max()
        rows.append(dict(p=p, eps=eps, ell=ell, lam=float(lam), eta_max_over_lam=float(np.abs(e).max() / lam), err_eta=float(np.abs(e - e_exact).max() / lam),
                         err_w=float(rel[sel].max()), wmax=float(np.abs(w_cf).max())))
    # the paper's unit: the raw (alpha, nu) ODE against Lemma 2.1's dual
    eps = 1e-4
    pu = U.PaperUnit(eps)
    L = U.LAW
    om, T = 0.7, 40.0
    e_unit = np.abs((1 - np.cos(om * np.linspace(0, T, 4001))) / om + 0.3 * np.linspace(0, T, 4001)).max()
    Amp = 0.5 * pu.lam / e_unit
    g = lambda t: Amp * (np.sin(om * t) + 0.3)
    a, m, G, gam = L['a_max'], L['m'], L['G'], L['gamma']

    def rhs(t, y):
        al, nu, e = y
        th = np.tanh(al / a)
        dw_da = nu ** 2 * (1 - th ** 2) / m
        dw_dn = 2 * nu * a * th / m
        return [G * m * g(t) * dw_da, gam * m / 4 * g(t) * dw_dn, g(t)]
    sol = solve_ivp(rhs, (0, T), [0.0, np.sqrt(eps), 0.0], t_eval=np.linspace(0, T, 401), method='DOP853', rtol=1e-12, atol=1e-18)
    al, nu, e = sol.y
    x = al / a
    eta_cf = pu.eta_of_x(x)
    Qc = nu ** 2 - pu.vs * np.cosh(2 * x)
    rows.append(dict(p=0, eps=eps, ell=1.0, lam=float(pu.lam), eta_max_over_lam=float(np.abs(e).max() / pu.lam), err_eta=float(np.abs(e - eta_cf).max() / pu.lam),
                     err_w=float(np.abs(Qc - Qc[0]).max() / pu.vs), wmax=float(np.abs(nu ** 2 * a * np.tanh(x) / m).max())))
    # Lemma 2.2 of the T1 note: the exact dormant sup of the paper's unit, its cos/sin^3 form and its mu^-3 form
    chain = []
    for eps in (1e-4, 1e-6, 1e-8, 1e-10):
        pu = U.PaperUnit(eps)
        for mu in (0.5, 0.2, 0.05, 0.025, 0.01, 1e-3):
            ex = pu.dormant_bound(mu)
            f1, f2 = pu.dormant_bound_T1(mu)
            chain.append(dict(eps=eps, mu=mu, applicable=bool(1.0 / (np.tan(np.pi * mu / 2) * pu.kap) < 1), exact=float(ex), cos_sin=float(f1), mu3=float(f2), lead=float(8 / np.pi ** 3 * L['a_max'] / (m * np.sqrt(2 * pu.vs)) * eps ** 1.5 * mu ** -3)))
    return dict(rows=rows, chain=chain)


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# B. ball arithmetic
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
def part_B():
    return dict(counts=U.ball_checks())


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# C. the repository's network flows (moving centres): the dual is the integral of the drive, and the coefficient is W of it
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
READOUT_CELLS = (('D2', 1e-6), ('D3', 1e-6), ('D4', 1e-4), ('D5', 1e-2), ('P0', 1e-6))


def eta_inverse(name, w, eps):
    """the dual carrying the coefficients w (array) and the wall"""
    if name == 'D2':
        return U.eta_of_w_2(w, eps), U.lam_2(eps)
    if name == 'P0':
        pu = U.PaperUnit(eps)
        return np.array([pu.eta_of_w(x) for x in np.atleast_1d(w)]), pu.lam
    p = int(name[1:])
    return U.eta_of_w_p(w, p, eps), U.lam_p(p, eps)


def run_with_dual(name, eps, tau_end=P.TAU_B, n_out=121):
    """the design target trained by the plain flow of `factorized_rbf.make_rhs` with the duals eta_j = int g_j as extra states; returns the record and the dual from the closed-form inverse"""
    L, K, c0 = P.setup()
    f = K.synth(np.array(P.CSTAR)[:, None], np.array(P.ASTAR))
    ro = P.make_readout(name)
    d, n, q = len(c0), 1, ro.q
    c0 = FR.pts(c0)
    scale = ro.scale(eps)
    base = FR.make_rhs(ro, K, f, P.LR_C, d, n)

    def rhs(t, z):
        zz = z[:d * (q + n)]
        out = base(t, zz)
        c = zz[d * q:].reshape(d, n)
        p_ = zz[:d * q].reshape(d, q)
        w = ro.w(p_)
        r = f - K.synth(c, w)
        g, _ = K.corr(c, r)
        return np.r_[out, g]
    p0 = ro.init(d, eps)
    z0 = np.r_[p0.ravel(), c0.ravel(), np.zeros(d)]
    atol = np.r_[1e-10 * np.maximum(np.abs(p0).ravel(), np.sqrt(eps)), np.full(d * n, 1e-12), np.full(d, 1e-12 * scale)]
    taus = np.linspace(0.0, tau_end, n_out)
    sol = solve_ivp(rhs, (0.0, tau_end * scale), z0, t_eval=taus * scale, method='LSODA', rtol=1e-10, atol=atol)
    if sol.status < 0 or sol.y.shape[1] != n_out:
        sol = solve_ivp(rhs, (0.0, tau_end * scale), z0, t_eval=taus * scale, method='Radau', rtol=1e-10, atol=atol)
    Z = sol.y.T
    Pm = Z[:, :d * q].reshape(-1, d, q)
    W = np.stack([ro.w(Pm[k]) for k in range(Pm.shape[0])])
    eta_int = Z[:, d * (q + n):]
    eta_inv = np.stack([eta_inverse(name, W[k], eps)[0] for k in range(W.shape[0])])
    lam = eta_inverse(name, np.array([0.0]), eps)[1]
    return dict(taus=taus, W=W, eta_int=eta_int, eta_inv=eta_inv, lam=float(lam), status=int(sol.status))


def part_C():
    rows = []
    for name, eps in READOUT_CELLS:
        r = run_with_dual(name, eps)
        err = np.abs(r['eta_int'] - r['eta_inv']) / r['lam']
        W, eta, lam = r['W'], r['eta_int'], r['lam']
        thr = FR.formed_threshold(P.ASTAR)
        formed = np.abs(W[-1]) > thr
        row = dict(name=name, eps=eps, lam=lam, err_max=float(err.max()), err_formed=float(err[:, formed].max()) if formed.any() else float('nan'),
                   err_dormant=float(err[:, ~formed].max()) if (~formed).any() else float('nan'), n_formed=int(formed.sum()), status=r['status'])
        if name != 'P0' and name != 'D2':
            p = int(name[1:])
            mu_obs = 1.0 - np.abs(eta[:, ~formed]).max(axis=0) / lam           # the margin of each dormant unit over the run, from the integrated dual
            wmax = np.abs(W[:, ~formed]).max(axis=0)
            bound = U.dormant_bound_p(mu_obs, p, eps)
            simple = U.dormant_bound_simple_p(mu_obs, p, eps)
            row.update(ratio_min=float((wmax / bound).min()), ratio_max=float((wmax / bound).max()), simple_over_exact_max=float((simple / bound).max()),
                       mu_min=float(mu_obs.min()), wmax_over_simple=float((wmax / simple).max()))
        elif name == 'D2':
            mu_obs = 1.0 - np.abs(eta[:, ~formed]).max(axis=0) / lam
            wmax = np.abs(W[:, ~formed]).max(axis=0)
            bound = U.dormant_bound_2(mu_obs, eps)
            row.update(ratio_min=float((wmax / bound).min()), ratio_max=float((wmax / bound).max()), mu_min=float(mu_obs.min()), simple_over_exact_max=1.0,
                       wmax_over_simple=float((wmax / (eps ** mu_obs)).max()))
        else:
            pu = U.PaperUnit(eps)
            mu_obs = 1.0 - np.abs(eta[:, ~formed]).max(axis=0) / lam
            wmax = np.abs(W[:, ~formed]).max(axis=0)
            bound = np.array([pu.dormant_bound(m) for m in mu_obs])
            t1 = np.array([pu.dormant_bound_T1(m)[1] for m in mu_obs])
            row.update(ratio_min=float((wmax / bound).min()), ratio_max=float((wmax / bound).max()), simple_over_exact_max=float((t1 / bound).max()), mu_min=float(mu_obs.min()),
                       wmax_over_simple=float((wmax / t1).max()))
        rows.append(row)
    return dict(rows=rows)


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# D. depth 1
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
def part_D():
    L, K, c0 = P.setup()
    f = K.synth(np.array(P.CSTAR)[:, None], np.array(P.ASTAR))
    c0p = FR.pts(c0)
    Phi = K.phi(c0p)
    N = Phi.shape[0]
    H = Phi.T @ Phi / N
    b0 = Phi.T @ f / N
    ell = P.LR_C
    thr = FR.formed_threshold(P.ASTAR)
    sv = np.linalg.svd(Phi, compute_uv=False)
    wls = np.linalg.lstsq(Phi, f, rcond=None)[0]
    normH = float(np.linalg.eigvalsh(H).max())
    # fixed centres: the flow of the repository (lr_c = 0) against the matrix exponential
    ro = FR.Linear(ell)
    T = 3000.0
    res = FR.train_flow(ro, K, f, c0, 1e-6, 1.0, 0.0, n_t=601, t_scale=T)
    wt, gt = U.d1_fixed_centres(H, b0, ell, res['taus'] * T)
    err = float(np.abs(res['W'] - wt).max())
    # the formation times: first t with |w_j| > theta on a fine log grid of the closed form (the flow agrees with it to err_fixed), and the leading prediction theta/(ell |g_j(0)|)
    tt = np.geomspace(1e-3, T, 3000)
    Wf, _ = U.d1_fixed_centres(H, b0, ell, tt)
    tj_sim, tj_pred, g0 = [], [], b0
    for j in range(len(c0)):
        k = np.flatnonzero(np.abs(Wf[:, j]) > thr)
        tj_sim.append(float(tt[k[0]]) if k.size else float('nan'))
        tj_pred.append(float(thr / (ell * abs(g0[j]))))
    cert, tj_c, Tj_c = U.d1_formation_certificate(g0, normH, ell, thr)
    ok = [bool(np.isfinite(tj_sim[j]) and tj_sim[j] <= tj_c[j] * (1 + 1e-9)) for j in np.flatnonzero(cert)]
    # the pure-ridge reading: w(t) in the long time and the LS limit
    return dict(sigma_min=float(sv.min()), sigma_max=float(sv.max()), n_ls_formed=int((np.abs(wls) > thr).sum()), n=len(c0), wls_min=float(np.abs(wls).min()),
                wls_absmax=float(np.abs(wls).max()), err_fixed=err, normH=normH, normg0=float(np.linalg.norm(g0)), thr=float(thr), n_cert=int(cert.sum()), cert_ok=ok,
                tj_sim=tj_sim, tj_pred=tj_pred, g0=[float(x) for x in g0], n_formed_T=int((np.abs(res['W'][-1]) > thr).sum()))


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# E. the sliding unit: first integral and overtaking margin
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
def part_E():
    eps, mu_c = 1e-6, 12.0
    W = lambda e: U.W_p(e, 3, eps)
    lam = U.lam_p(3, eps)
    from numpy.polynomial.legendre import leggauss
    xg, wg = leggauss(60)
    # (1) a Gaussian correlation landscape; the unit starts at the margin x_a on the flank and is followed until its centre is 0.5 sigma from the maximum (or its dual at 0.9999 wall)
    sig, cs_, A = 0.12, 0.3, 0.06
    G = lambda c: A * np.exp(-(c - cs_) ** 2 / (2 * sig ** 2))
    dG = lambda c: -A * (c - cs_) / sig ** 2 * np.exp(-(c - cs_) ** 2 / (2 * sig ** 2))
    c0 = cs_ - 1.3 * sig
    rows = []
    for xa in (0.2, 0.05, 0.02):
        eta_a = (1 - xa) * lam
        t, eta, c = U.integrate_sliding_unit(G, dG, eta_a, c0, mu_c, W, 3.0 * lam * xa / G(c0), lam, c_stop=(cs_, 0.5 * sig))

        def Theta(cend):
            cc = 0.5 * (cend - c0) * xg + 0.5 * (cend + c0)
            return float(0.5 * (cend - c0) * (wg * G(cc) / (mu_c * dG(cc))).sum())
        Ps = lambda e: float(U.Psi_p(e, 3, eps))
        k = len(t)
        sel = list(range(0, k, max(1, k // 60))) + [k - 1]
        Ecur = [Ps(eta[i]) - Ps(eta_a) - Theta(c[i]) for i in sel]
        scale = max(abs(Ps(eta[k - 1]) - Ps(eta_a)), 1e-300)
        rows.append(dict(x_a=xa, steps=k, drift=float(np.abs(Ecur).max() / scale), dc_over_sig=float(abs(c[k - 1] - c0) / sig), x_end=float(1 - eta[k - 1] / lam),
                         w_end=float(W(eta[k - 1])), G_ratio=float(G(c[k - 1]) / G(c0))))
    # (2) the linear landscape from the start of the run (eta = 0, c = 0): the exact first integral g^2 = g0^2 + 2 kappa Psi*(eta) and the first-order lead against the ODE at the frozen margin 0.02
    lin = []
    for h_rel in (0.3, 1.0, 3.0, 10.0):                    # h/g0 in 1/length (c in actual units, sigma = 0.12)
        g0 = 0.02
        h = h_rel * g0
        Gl = lambda c, g0=g0, h=h: g0 + h * c
        dGl = lambda c, h=h: h
        x_fr = 0.02
        T = lam * (1 - x_fr) / g0                        # the frozen dual reaches the margin x_fr at this time
        t, eta, c = U.integrate_sliding_unit(Gl, dGl, 0.0, 0.0, mu_c, W, T, lam, wall_frac=0.9999)
        kappa = mu_c * h ** 2
        sel = list(range(0, len(t), max(1, len(t) // 60))) + [len(t) - 1]
        gap = max(abs(Gl(c[i]) ** 2 - (g0 ** 2 + 2 * kappa * float(U.Psi_p(eta[i], 3, eps)))) / Gl(c[i]) ** 2 for i in sel)
        i_end = len(t) - 1
        reached = bool(abs(t[i_end] - T) < 1e-9 * T)
        sg_fr = g0 * t[i_end] / lam
        lead = (eta[i_end] - g0 * t[i_end]) / lam
        pred1 = U.lead_first_order_D3(kappa, g0, eps, sg_fr)
        lin.append(dict(h_rel=h_rel, kappa_over_g2=float(kappa / g0 ** 2), first_integral_err=float(gap), x_end=float(1 - sg_fr), lead=float(lead), pred=float(pred1),
                        pred_lead_form=float(U.drift_lead(kappa, g0, eps, 1 - sg_fr, 1.0, A_RO['D3'])), c_over_sig=float(c[i_end] / 0.12), reached=reached))
    return dict(landscape=rows, linear=lin)


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# L. fixed centres: the duals of the network against the inverse scale space of the weighted l1 norm (Corollary 5)
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
L_CELLS = (('D3', (1e-4, 1e-6, 1e-8, 1e-10)), ('P0', (1e-4, 1e-6, 1e-8, 1e-10)), ('D2', (1e-6, 1e-12, 1e-24, 1e-48, 1e-96)), ('D1', (1.0,)))
L_TAU = 25.0                                    # before the fourth entry of the fixed-centre path (27.4)


def part_L():
    from src.core import mirror_free as MF
    L, K, c0 = P.setup()
    f = K.synth(np.array(P.CSTAR)[:, None], np.array(P.ASTAR))
    c0p = FR.pts(c0)
    Phi = K.phi(c0p)
    d = len(c0)
    iss = MF.iss_path(Phi, f, np.ones(d))
    tk, Uk = np.array(iss['taus']), np.array(iss['U'])
    Wk = np.array(iss['W'])
    tg = np.linspace(0.0, L_TAU, 501)
    u_iss = np.stack([np.interp(tg, tk, Uk[:, j]) for j in range(d)], 1)
    seg = np.clip(np.searchsorted(tk, tg, side='right') - 1, 0, len(Wk) - 1)
    w_iss = Wk[seg]
    fit_iss = w_iss @ Phi.T
    ents = [(float(t), int(j)) for t, j, kind in iss['events'] if kind == 'enter' and t <= L_TAU]
    rows = []
    for name, epss in L_CELLS:
        for eps in epss:
            ro = P.make_readout(name)
            if name == 'D1':
                res = FR.train_flow(ro, K, f, c0, 1e-6, L_TAU, 0.0, n_t=501, t_scale=float(FR.PaperLaw(P.law()).scale(1e-6)))
                u_net = None
                sup = float('nan')
            else:
                res = FR.train_flow(ro, K, f, c0, eps, L_TAU, 0.0, n_t=501)
                lam = eta_inverse(name, np.array([0.0]), eps)[1]
                u_net = np.stack([eta_inverse(name, res['W'][k], eps)[0] / lam for k in range(501)])
                sup = float(np.abs(u_net - u_iss).max())
            fit_net = res['W'] @ Phi.T
            l2 = float(np.sqrt(np.mean((fit_net - fit_iss) ** 2)) / np.sqrt(np.mean(f ** 2)))
            thr = FR.formed_threshold(P.ASTAR)
            rows.append(dict(name=name, eps=float(eps), sup_dual=sup, l2=l2, n_formed=int((np.abs(res['W'][-1]) > thr).sum()), status=str(res['status'])))
    return dict(rows=rows, entries=ents, n_events=len(tk), tau=L_TAU)


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# K. a one-parameter family of ties: one feature between two kernels; the critical eps by bisection
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
K_FAMILY = dict(D3=dict(deltas=(0.0003, 0.0006, 0.0012, 0.0025, 0.005, 0.01, 0.02, 0.04), lo=-10.0, hi=-2.0, nit=8),
                P0=dict(deltas=(0.0006, 0.0012, 0.0025, 0.005, 0.01, 0.02), lo=-9.5, hi=-2.5, nit=7))


def part_K():
    """feature of amplitude 0.8 at 0.05 + delta (the midpoint of the kernels at 0 and 0.1 of the dictionary of 21), delta small: the path enters the kernel at 0.1 and the kernel at 0 is the
    runner-up, margin ~ 0.0035 delta/0.001. The set differs from the path's for eps >= eps* (bisection on log10 eps between lo and hi, nit steps; flip tested by the plain flow)."""
    L, K, c0 = P.setup()
    A = 0.8
    out = {}
    for name, sp in K_FAMILY.items():
        rows = []
        for dlt in sp['deltas']:
            cs, As = np.array([0.05 + dlt]), np.array([A])
            f = K.synth(cs[:, None], As)
            rec = path_record(L, K, c0, cs, As)
            ref = P.reference(L, K, c0, f)
            tau_end = 1.05 * max(ref['tau']) + 1.0
            ro = P.make_readout(name)

            def flips(le):
                c = P.cell(ro, K, f, c0, cs, As, 10.0 ** le, ref, tau_end)
                return not (c['set_equal'] and c['order_equal'])
            lo, hi = sp['lo'], sp['hi']
            f_lo, f_hi = flips(lo), flips(hi)
            n_runs = 2
            if f_hi and not f_lo:
                for _ in range(sp['nit']):
                    mid = 0.5 * (lo + hi)
                    n_runs += 1
                    if flips(mid):
                        hi = mid
                    else:
                        lo = mid
            pairs = stage_pairs(rec)
            rows.append(dict(delta=float(dlt), margin=float(rec['margin']), flip_lo=bool(f_lo), flip_hi=bool(f_hi), log_eps_lo=float(lo), log_eps_hi=float(hi), pred=critical_eps(pairs, A_RO[name]),
                             rho=rho_stat(pairs), runs=n_runs))
        out[name] = rows
    return out


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# P. the sliding paths of the 1034 one-dimensional benchmark targets: the stagewise drives, the stages' duals, the margins
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
def population_targets(cert_raw=CERT_RAW):
    z = np.load(cert_raw, allow_pickle=True)
    out = []
    for fam in ('A', 'B', 'D'):
        for idx, r in enumerate(z[fam][0]):
            if not r.get('skip'):
                out.append((fam, idx, np.asarray(r['cs'], float), np.asarray(r['As'], float)))
    return out


def path_record(L, K, c0, cs, As):
    """the stagewise data of the sliding path of a target: for every stage k (between event k and k+1): the drives g, the centre-slopes h = <r d_c phi> at the starting centres under the
    stage's residual, the duals at both ends; the events; the runner-up margins"""
    f = K.synth(cs[:, None], As)
    c0a = np.asarray(c0, float)[:, None]
    d = len(c0)
    R = RN.sliding_iss_nd(L, K, f, c0a, np.ones(d), polish=True)
    m_all, m_nev = P.margins(K, f, R, d)
    Pd = RN.path_duals_nd(K, f, R)
    taus = np.array(Pd['taus'])
    Uabs = np.abs(Pd['U'])
    ns = len(taus) - 1
    G = np.zeros((ns, d))
    H = np.zeros((ns, d))
    for k in range(ns):
        g, h = K.corr(c0a, Pd['res'][k])
        G[k], H[k] = g, h[:, 0]
    ev = [(float(t), int(j), 1 if kind == 'enter' else 0) for t, j, kind, c in R['events']]
    return dict(cs=cs, As=As, taus=taus, G=G, H=H, U0=Uabs[:-1], U1=Uabs[1:], ev=ev, margin=m_all, margin_never=m_nev, exits=int(sum(1 for e in ev if e[2] == 0)),
                relaxed=bool(all(R['relaxed'])), failed=int(R.get('solver_failed', 0)) + int(bool(R.get('out_of_budget', False))))


def part_P(limit=None):
    warnings.filterwarnings('ignore')
    L, K, c0 = P.setup()
    pop = population_targets()
    recs = []
    for i, (fam, idx, cs, As) in enumerate(pop[:limit]):
        r = path_record(L, K, c0, cs, As)
        r.update(fam=fam, idx=idx)
        recs.append(r)
    return dict(recs=recs)


def stage_pairs(rec, mu_c=MU_C):
    """the (stage, competitor) pairs of a path record: arrays (x0, x1, q) with x0, x1 the competitor's margins at the start and the end of the stage, q = mu_c h^2/g^2: the competitor is
    an inactive kernel other than the stage's entrant, with a positive margin at the end of the stage"""
    ev = rec['ev']
    ns = rec['G'].shape[0]
    act = np.zeros(rec['G'].shape[1], bool)
    e = 0
    x0s, x1s, qs, ids = [], [], [], []
    for k in range(ns):
        while e < len(ev) and ev[e][0] <= rec['taus'][k] + 1e-12:
            act[ev[e][1]] = bool(ev[e][2])
            e += 1
        ent = [ev[i][1] for i in range(len(ev)) if abs(ev[i][0] - rec['taus'][k + 1]) <= 1e-12 and ev[i][2] == 1]
        g, h = rec['G'][k], rec['H'][k]
        x0, x1 = 1.0 - rec['U0'][k], 1.0 - rec['U1'][k]
        ok = ~act & (x1 > 1e-12) & (np.abs(g) > 1e-300)
        for j in ent:
            ok[j] = False
        j_ = np.flatnonzero(ok)
        x0s.append(np.maximum(x0[j_], x1[j_] * 1.0001))
        x1s.append(x1[j_])
        qs.append(mu_c * h[j_] ** 2 / g[j_] ** 2)
        ids.append(np.stack([np.full(j_.size, k), j_], 1))
    cat = lambda a, dt=float: np.concatenate(a) if a else np.zeros(0, dt)
    return cat(x0s), cat(x1s), cat(qs), (np.concatenate(ids) if ids else np.zeros((0, 2), int))


def overtake_score(pairs, eps, a_ro):
    """max over the pairs of (the overtaking margin of the competitor)/(its margin at the entry): >= 1 predicts that the competitor overtakes (flip)"""
    x0, x1, q, _ = pairs
    if x1.size == 0:
        return 0.0
    return float((U.overtake_margin(q * a_ro * eps, x0) / x1).max())


def critical_eps(pairs, a_ro, lo=1e-20, hi=1.0):
    """the eps at which the score reaches 1 (the score increases with eps): bisection in log eps; lo if it is below 1 at hi... returns (eps*, flag): flag 0 inside, -1 never within [lo, hi]"""
    if pairs[1].size == 0:
        return float('nan')
    f = lambda e: overtake_score(pairs, e, a_ro) - 1.0
    if f(hi) < 0:
        return float('nan')
    if f(lo) > 0:
        return lo
    a, b = np.log(lo), np.log(hi)
    for _ in range(60):
        m_ = 0.5 * (a + b)
        if f(np.exp(m_)) > 0:
            b = m_
        else:
            a = m_
    return float(np.exp(0.5 * (a + b)))


def rho_stat(pairs):
    """the leading-order critical square root: min over pairs of x1/sqrt(q a_ro) (without a_ro): the competitor flips when sqrt(eps a_ro) > x1/sqrt(q)"""
    x0, x1, q, _ = pairs
    if x1.size == 0:
        return float('nan')
    return float((x1 / np.sqrt(np.maximum(q, 1e-300))).min())


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# F. the 200 targets: what the scalar theorem predicts, what the overtaking model adds, against the stored runs
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
def match_population(recs, T):
    """index of each of the 200 targets of Part C in the population records (by the centres and amplitudes)"""
    idx = {}
    for i, r in enumerate(recs):
        idx[(tuple(np.round(r['cs'], 12)), tuple(np.round(r['As'], 12)))] = i
    return [idx.get((tuple(np.round(np.asarray(t['cs']), 12)), tuple(np.round(np.asarray(t['As']), 12)))) for t in T]


def cells_of(t, name, eps):
    return [c for c in t['cells'] if c['name'] == name and c['eps'] == eps][0]


def confusion(pred, obs):
    pred, obs = np.asarray(pred, bool), np.asarray(obs, bool)
    return dict(tp=int((pred & obs).sum()), fp=int((pred & ~obs).sum()), fn=int((~pred & obs).sum()), tn=int((~pred & ~obs).sum()))


def part_F(Pres):
    T = SC.load_targets(PARAM_RAW)
    recs = Pres['recs']
    mi = match_population(recs, T)
    out = dict(n_matched=int(sum(m is not None for m in mi)), cells=[])
    CELLS = (('D3', 1e-6), ('D3', 1e-8), ('P0', 1e-6), ('P0', 1e-8))
    for name, eps in CELLS:
        theta_t = np.array([formed_theta(t['As']) for t in T])
        mnev = np.array([t['ref']['margin_never'] for t in T])
        mall = np.array([t['ref']['margin'] for t in T])
        if name == 'D3':
            mu_th = np.array([U.mu_theta_p(th, 3, eps) for th in theta_t])
        else:
            pu = U.PaperUnit(eps)
            mu_th = np.array([pu.mu_theta(th) for th in theta_t])
        setm = np.array([not cells_of(t, name, eps)['set_equal'] for t in T])
        anym = np.array([(not cells_of(t, name, eps)['set_equal']) or (not cells_of(t, name, eps)['order_equal']) for t in T])
        pred_scalar = np.where(np.isnan(mnev), False, mnev < mu_th)
        sc = np.array([overtake_score(stage_pairs(recs[m]), eps, A_RO[name]) if m is not None else np.nan for m in mi])
        pred_drift = np.where(np.isnan(sc), False, sc >= SCORE_FLIP)
        mm_obs = mall[setm]
        mm_pred = mall[pred_drift]
        # dormant coefficients against the bound at the path's margin (units that never enter): the ratio dormant_wmax/bound(mu_never), for targets with mu_never >= 100 sqrt(eps)
        ratios = []
        for t, mn in zip(T, mnev):
            if np.isnan(mn) or mn < 100 * np.sqrt(eps):
                continue
            c = cells_of(t, name, eps)
            if not (c['n_formed'] == len(t['ref']['set'])):
                continue
            b = U.dormant_bound_p(mn, 3, eps) if name == 'D3' else U.PaperUnit(eps).dormant_bound(mn)
            ratios.append(c['dormant_wmax'] / b)
        out['cells'].append(dict(name=name, eps=eps, n_set_mismatch=int(setm.sum()), n_any_flip=int(anym.sum()), conf_scalar=confusion(pred_scalar, setm),
                                 conf_drift_any=confusion(pred_drift, anym), conf_drift_set=confusion(pred_drift, setm), max_margin_obs=float(mm_obs.max()) if mm_obs.size else float('nan'),
                                 max_margin_pred=float(mm_pred.max()) if mm_pred.size else float('nan'), n_pred_scalar=int(pred_scalar.sum()), n_pred_drift=int(pred_drift.sum()),
                                 mu_theta_min=float(mu_th.min()), mu_theta_max=float(mu_th.max()), ratios=np.array(ratios), pred_scalar_margins_max=float(mnev[pred_scalar].max()) if pred_scalar.any() else float('nan'),
                                 score=sc, pred_scalar_margins=np.sort(mnev[pred_scalar]), obs_set_margins=np.sort(mnev[setm]), theta_t=theta_t))
    # observed exponents (the paper's statistics): max margin of a mismatch and the count, between 1e-6 and 1e-8
    # depth 2: the exact formation margin against the rule, the stored cells (200 targets x 2 eps, and Part D's 21 x 6)
    d2 = []
    for eps in (1e-6, 1e-12):
        pe, pr, ob = [], [], []
        for t in T:
            th = formed_theta(t['As'])
            mn = t['ref']['margin_never']
            mn_ = np.inf if np.isnan(mn) else mn
            pe.append(mn_ > U.mu_theta_2(th, eps))
            pr.append(mn_ > U.rule_depth2(0.01, eps))
            ob.append(cells_of(t, 'D2', eps)['set_equal'])
        d2.append(dict(eps=eps, exact=confusion(pe, ob), rule=confusion(pr, ob), n_obs_pass=int(np.sum(ob))))
    zp = np.load(PARAM_RAW, allow_pickle=True)
    Dp = list(zp['D'][0])                         # Part D: the design, then the first 20 targets of Part C
    agree_e = agree_r = tot = 0
    for k, r in enumerate(Dp):
        th = formed_theta(P.ASTAR) if k == 0 else formed_theta(T[k - 1]['As'])
        mn = r['margin_never']
        mn_ = np.inf if (mn is None or np.isnan(mn)) else mn
        for row in r['rows']:
            tot += 1
            ob = bool(row['set_equal'])
            agree_e += int((mn_ > U.mu_theta_2(th, row['eps'])) == ob)
            agree_r += int((mn_ > U.rule_depth2(0.01, row['eps'])) == ob)
    out['d2'] = d2
    out['d2_partD'] = dict(total=tot, exact_agree=agree_e, rule_agree=agree_r)
    # the design (Part B): the scalar prediction 'path reached iff mu > mu_theta' for every readout and eps
    Bd = zp['B'][0]
    mu_d = Bd['ref']['margin_never']
    th_d = formed_theta(P.ASTAR)
    des = []
    for row in Bd['rows']:
        ro, eps = row['readout'], row['eps']
        if ro in ('D3', 'D4', 'D5'):
            mt = U.mu_theta_p(th_d, int(ro[1:]), eps)
        elif ro == 'D2':
            mt = U.mu_theta_2(th_d, eps)
        elif ro == 'P0':
            mt = U.PaperUnit(eps).mu_theta(th_d)
        else:
            continue
        des.append(dict(readout=ro, eps=eps, mu_theta=float(mt), pred_pass=bool(mu_d > mt), obs_pass=bool(row['set_equal']), n_formed=int(row['n_formed']), ok=(row['status'] == 'ok')))
    disp = []
    for row in Bd['rows']:
        ro, eps = row['readout'], row['eps']
        if ro in ('D3', 'D4', 'D5') and eps <= 1e-6 and row['status'] == 'ok':
            pp = int(ro[1])
            ps = float(U.Psi_p((1.0 - mu_d) * U.lam_p(pp, eps), pp, eps))
            disp.append(dict(readout=ro, eps=eps, disp=float(row['dormant_disp']), psi=ps, ratio=float(row['dormant_disp']) / ps))
    out['design'] = dict(mu=float(mu_d), theta=float(th_d), rows=des, disp=disp)
    return out


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# G. the per-target critical eps of the scan against the model's
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
def observed_critical(cells_by_eps):
    """from a target's cells {eps: cell}: (eps_hi, eps_lo) with a flip (set or order differs) at eps_hi and none at the next lower grid value eps_lo; or (None, eps_min) if no flip at all, or
    (eps_min, None) if it flips at the lowest value (censored); the smallest flipped eps counts"""
    es = sorted(cells_by_eps)
    fl = {e: (not cells_by_eps[e]['set_equal']) or (not cells_by_eps[e]['order_equal']) for e in es}
    flipped = [e for e in es if fl[e]]
    if not flipped:
        return None, es[0]
    e_lo_f = flipped[0]
    below = [e for e in es if e < e_lo_f]
    if not below:
        return e_lo_f, None
    return e_lo_f, below[-1]


def part_G(Pres, scan):
    T = SC.load_targets(PARAM_RAW)
    recs = Pres['recs']
    mi = match_population(recs, T)
    out = dict(rows=[])
    for name in ('D3', 'P0'):
        tab = scan.get(name, {})
        targets = sorted({i for (i, _) in tab})
        rows = []
        for i in targets:
            cells = {e: c for (ii, e), c in tab.items() if ii == i}
            hi, lo = observed_critical(cells)
            m = mi[i]
            pairs = stage_pairs(recs[m])
            ec = critical_eps(pairs, A_RO[name])
            rows.append(dict(i=int(i), margin=float(T[i]['ref']['margin']), eps_hi=hi if hi is not None else float('nan'), eps_lo=lo if lo is not None else float('nan'),
                             pred=ec, censored_low=bool(lo is None), no_flip=bool(hi is None), n_eps=len(cells)))
        out[name] = rows
    return out


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# H. the population: the model's counts and maxima as functions of eps, the exponent
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
def local_exponent(e1, e2, v1, v2):
    return float(np.log(v1 / v2) / np.log(e1 / e2))


def part_H(Pres):
    recs = Pres['recs']
    good = [r for r in recs if not r['failed'] and r['exits'] == 0 and r['G'].shape[0] > 0]
    out = dict(n=len(recs), n_good=len(good), eps=np.array(EPS_H), by_name={})
    pairs = [stage_pairs(r) for r in good]
    marg = np.array([r['margin'] for r in good])
    rng = np.random.default_rng(20261110)
    boots = rng.integers(0, len(good), size=(BOOT_N, BOOT_SIZE))
    for name in ('D3', 'P0'):
        cnt, mx, q95 = [], [], []
        flips = {}
        for eps in EPS_H:
            sc = np.array([overtake_score(p_, eps, A_RO[name]) for p_ in pairs])
            fl = sc >= SCORE_FLIP
            flips[eps] = fl
            cnt.append(int(fl.sum()))
            mx.append(float(marg[fl].max()) if fl.any() else float('nan'))
            q95.append(float(np.percentile(marg[fl], 95)) if fl.sum() >= 5 else float('nan'))
        rho = np.array([rho_stat(p_) for p_ in pairs])
        cc = np.concatenate([np.sqrt(A_RO[name] * p_[2])[p_[1] < 0.1] for p_ in pairs])         # c = sqrt(a_ro q): the competitor flips when x < c sqrt(eps)
        # bootstrap: samples of BOOT_SIZE targets (like the 200 of Part C): the exponents between 1e-6 and 1e-8 of the count and of the largest margin among the flips
        f6, f8 = flips[10.0 ** -6], flips[10.0 ** -8]
        ex_n, ex_m = [], []
        for idx in boots:
            n6, n8 = int(f6[idx].sum()), int(f8[idx].sum())
            if n6 > 0 and n8 > 0:
                ex_n.append(local_exponent(1e-6, 1e-8, n6, n8))
                ex_m.append(local_exponent(1e-6, 1e-8, marg[idx][f6[idx]].max(), marg[idx][f8[idx]].max()))
        out['by_name'][name] = dict(count=np.array(cnt), maxmargin=np.array(mx), q95=np.array(q95), rho=rho, c_pct=np.percentile(cc, [10, 50, 90, 99]), c_n=len(cc), c_frac30=float((cc > 30).mean()), boot_n=len(ex_n), boot_count_exp=np.percentile(ex_n, [5, 25, 50, 75, 95]),
                                    boot_max_exp=np.percentile(ex_m, [5, 25, 50, 75, 95]))
    return out


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# the report
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
def g(x, f='.3g'):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return 'n/a'
    return format(x, f)


def table(head, rows):
    out = ['| ' + ' | '.join(head) + ' |', '|' + '|'.join(['---'] * len(head)) + '|']
    out += ['| ' + ' | '.join(str(c) for c in r) + ' |' for r in rows]
    return '\n'.join(out) + '\n'


def conf_str(c):
    return f"TP {c['tp']}, FP {c['fp']}, FN {c['fn']}, TN {c['tn']}"


def report(A=None, B=None, C=None, D=None, E=None, Lf=None, Kf=None, Pr=None, F=None, G=None, H=None, scan=None):
    L = ['# Which readouts select: theorems for the readouts of depth p, their checks, and the mismatch margin (round theory3, direction D)\n',
         'The statements (with the proofs, in the note outside the repository) concern the SCALAR unit: whatever the other units and its own centre do, a unit of depth p enters the network only through '
         'its dual eta = int g, g = <r phi(.;c)> the drive, and its coefficient is an explicit function of eta. Everything below tests this against the exact ODE (A), ball arithmetic (B), the '
         'repository\'s network flows with moving centres (C), the depth-1 flow (D), the sliding unit (E), the fixed-centre limit (L), a family of ties at fixed geometry (K), the paths of the benchmark population (P), the stored runs of 200 benchmark targets (F), the per-target scan (G) and the population model (H). '
         'Status of each statement is in the last section.\n']
    if A:
        L.append('## A. The closed forms against the coordinate ODE (prescribed drive)\n')
        L.append('Drive g(t) = Amp (sin 0.7t + 0.3), Amp so that max|eta| = 0.9 wall (0.5 for the paper\'s unit); DOP853 with rtol 1e-12; w from (u, v) against W(eta) with eta from the same integration; '
                 'row p = 0 is the paper\'s unit ((alpha, nu) ODE against Lemma 2.1\'s dual; the last column is the drift of the conserved Q = nu^2 - v_s cosh(2 alpha/a_max) over v_s).\n')
        L.append(table(['p', 'eps', 'ell', 'wall', 'max eta/wall', 'err eta/wall', 'max rel err w (p=0: Q drift)', 'max w'],
                       [[r['p'], g(r['eps']), g(r['ell']), g(r['lam']), g(r['eta_max_over_lam']), g(r['err_eta']), g(r['err_w']), g(r['wmax'])] for r in A['rows']]))
        L.append('The paper\'s unit: exact dormant supremum, the cos/sin^3 form and the mu^-3 form of Lemma 2.2 (T1 note), and the leading constant 8/pi^3 a/(m sqrt(2 v_s)) eps^(3/2) mu^-3:\n')
        L.append(table(['eps', 'mu', 'exact sup', 'cos/sin^3 form', 'mu^-3 form', 'leading', 'exact <= cos/sin^3 <= mu^-3'],
                       [[g(r['eps']), g(r['mu']), g(r['exact']), g(r['cos_sin']), g(r['mu3']), g(r['lead']), ('yes' if r['exact'] <= r['cos_sin'] * (1 + 1e-12) <= r['mu3'] * (1 + 1e-12) else 'NO') if r['applicable'] else 'n/a (s_eps >= 1: the lemma needs cot(pi mu/2) < kappa)'] for r in A['chain']]))
    if B:
        L.append('## B. Ball arithmetic (python-flint, 400 bits)\n')
        names = {'B1': 'dormant supremum s^p[mu^-q - (2-mu)^-q] is exact and below s^p mu^-q', 'B2': 'distance form (ell p (p-2) d)^-q[1 - (x/(2-x))^q]',
                 'B3': 'u\' = ell p u^(p-1) g, v\' = -ell p v^(p-1) g (arb_series derivative of the closed forms)', 'B4': 'mass identity (rigorous quadrature) and Psi_p <= eps b mu^(-2/(p-2))/(2 ell p)',
                 'B5': 'formation distance bracket d^0 (1 - r^q)^(1/q) < d_theta < d^0, d^0 = theta^(-1/q)/(ell p (p-2))', 'B6': 'depth 2: exact dormant value, distance form, the bracket of mu_theta and the exact d_theta'}
        L.append(table(['check', 'statement', 'cases', 'failed'], [[k, names[k], v[0], v[1]] for k, v in B['counts'].items()]))
        L.append('Grid: p in 3..6, eps in {1e-4, 1e-6, 1e-8, 1e-12}, b in {1, 3/2}, ell in {1, 1/4}, mu in {1/2, 1/10, 1/100, 1/1000}, theta in {1/100, 1/1000}; inputs exact rationals; "failed" = a comparison not certified.\n')
    if C:
        L.append('## C. The repository\'s network flows (moving centres): the dual is the integral of the drive\n')
        L.append('Design target (3 features, 21 kernels), plain flow of `factorized_rbf.make_rhs` with the 21 duals eta_j = int g_j dt integrated as extra states (LSODA, rtol 1e-10), horizon tau = 20; w_j(t) is inverted '
                 'through W^-1 and compared with the integrated eta_j. "ratio" = max over the unformed units of (max_t |w_j|)/(the exact dormant supremum at the margin mu_j = 1 - max_t|eta_j|/wall actually reached): '
                 'the theorem says = 1 up to the sampling of the time grid (121 points) and never above 1.\n')
        L.append(table(['readout', 'eps', 'wall', 'max |eta_int - eta_W^-1|/wall', 'formed units', 'smallest dormant margin', 'ratio min', 'ratio max', 'bound-form/exact (max)', 'max w / bound-form'],
                       [[r['name'], g(r['eps']), g(r['lam']), g(r['err_max']), r['n_formed'], g(r['mu_min']), g(r['ratio_min']), g(r['ratio_max']), g(r['simple_over_exact_max']), g(r['wmax_over_simple'])] for r in C['rows']]))
        L.append('(bound-form: (eps)^(p/2) mu^(-p/(p-2)) for depth p; eps^mu for depth 2; Lemma 2.2\'s mu^-3 form for the paper\'s unit.)\n')
    if D:
        L.append('## D. Depth 1 (w = ell eta, no wall)\n')
        L.append(f"Design target, 21 kernels at fixed centres: sigma_min(Phi) = {g(D['sigma_min'])}, sigma_max = {g(D['sigma_max'])} (full column rank, so the flow converges to the least-squares solution); "
                 f"the LS solution has {D['n_ls_formed']} of {D['n']} coefficients above the formed threshold {g(D['thr'])} (smallest |w_LS| = {g(D['wls_min'])}, largest {g(D['wls_absmax'])}). "
                 f"The repository's depth-1 flow with lr_c = 0 against H^+(I - exp(-ell H t))b over t <= 3000: max abs difference {g(D['err_fixed'])}. At t = 3000 {D['n_formed_T']} units are above the threshold. "
                 f"Formation-time certificate (Theorem 4): ||H|| = {g(D['normH'])}, ||g(0)|| = {g(D['normg0'])}; units with g_j(0)^2 > 2 theta ||H|| ||g(0)||: {D['n_cert']}, "
                 f"all formed by 2 theta/(ell |g_j(0)|): {all(D['cert_ok']) if D['cert_ok'] else 'vacuous (none certified)'}.\n")
        sim, pr, g0 = np.array(D['tj_sim']), np.array(D['tj_pred']), np.array(D['g0'])
        order = np.argsort(-np.abs(g0))[:8]
        L.append('Formation time of the eight units with the largest initial drive (first t with |w_j| > theta, against theta/(ell |g_j(0)|)):\n')
        L.append(table(['unit', 'g_j(0)', 't formed (flow)', 'theta/(ell |g_j(0)|)', 'ratio'], [[int(j), g(g0[j]), g(sim[j]), g(pr[j]), g(sim[j] / pr[j])] for j in order]))
    if E:
        L.append('## E. The sliding unit: first integral and overtaking margin\n')
        L.append('Model: eta\' = G(c), c\' = mu_c W(eta) G\'(c) (frozen residual; depth 3, eps = 1e-6, mu_c = 12). First integral Psi*(eta) - Theta(c), Theta\' = G/(mu_c G\'). '
                 'Gaussian landscape G(c) = 0.06 exp(-(c - 0.3)^2/(2 sigma^2)), sigma = 0.12, start 1.3 sigma from the maximum, margin x_a of the start:\n')
        L.append(table(['x_a', 'steps', 'drift of the first integral (relative)', '|c - c0|/sigma', 'margin at the end', 'w at the end', 'G(c)/G(c0)'],
                       [[g(r['x_a']), r['steps'], g(r['drift']), g(r['dc_over_sig']), g(r['x_end']), g(r['w_end']), g(r['G_ratio'])] for r in E['landscape']]))
        L.append('Linear landscape G = g0 + h c (g0 = 0.02) from the start of the run (eta = 0, c = 0), followed until the frozen dual is at margin 0.02: the exact relation g^2 = g0^2 + 2 kappa Psi*(eta), kappa = mu_c h^2, '
                 'and the lead (eta - g0 t)/wall against the first order in kappa, kappa eps/(6 g0^2)[1/(1-s) - 1/(1+s) - 2s], and against its leading form kappa a_ro eps/g0^2 (1/x - 1) (a_ro = 1/(6 ell)); '
                 'the displacement is in widths of the Gaussian kernel (a linear landscape beyond ~1 width is not a landscape of the benchmark; the large-kappa rows show where the first order stops):\n')
        L.append(table(['h/g0 (1/length)', 'kappa/g0^2', 'first-integral error (relative, max)', 'frozen margin', 'lead (ODE)', 'first order', 'leading form', 'displacement/sigma'],
                       [[g(r['h_rel']), g(r['kappa_over_g2']), g(r['first_integral_err']), g(r['x_end']), g(r['lead']), g(r['pred']), g(r['pred_lead_form']), g(r['c_over_sig'])] for r in E['linear']]))
    if Lf:
        L.append('## L. Fixed centres: the network against the inverse scale space of the weighted l1 norm (Corollary 5)\n')
        L.append(f"Design target, 21 kernels, centres frozen (lr_c = 0), time tau = t/wall, window [0, {g(Lf['tau'])}] ({len(Lf['entries'])} entries of the fixed-centre path: "
                 f"{', '.join(f'kernel {j} at {g(t)}' for t, j in Lf['entries'])}). Reference: `mirror_free.iss_path` (equal walls). sup dual = max over units and 501 times of |u_net - u_ISS|, u_net = (dual from the coefficient through the closed form)/wall; "
                 f"L2 = rms over time and samples of the fitted function's difference over the rms of the target. Depth 1 has no dual in this sense (its flow is the l2 flow); it is run to the same window of tau at the paper's wall.\n")
        rows = []
        for r in Lf['rows']:
            ex = ''
            if r['name'] in ('D3', 'P0') and np.isfinite(r['sup_dual']):
                ex = g(r['sup_dual'] / np.sqrt(r['eps']), '.2f')
            elif r['name'] == 'D2':
                ex = g(r['sup_dual'] * np.log(1 / r['eps']), '.2f')
            rows.append([r['name'], g(r['eps']), g(r['sup_dual']), ex, g(r['l2']), r['n_formed']])
        L.append(table(['readout', 'eps', 'sup dual error', 'sup/sqrt(eps) (depth 3, paper) or sup x ln(1/eps) (depth 2)', 'L2 fit error', 'units formed at the end'], rows))
        for name in ('D3', 'P0', 'D2'):
            rr = [r for r in Lf['rows'] if r['name'] == name]
            sl = float(np.polyfit(np.log([r['eps'] for r in rr]), np.log([r['sup_dual'] for r in rr]), 1)[0])
            L.append(f"{name}: slope of log sup dual error on log eps: {g(sl, '.3f')} (1/2 for a sqrt(eps) rate; for depth 2 the rate is 1/ln(1/eps), slope ~ 0).\n")
    if Kf:
        L.append('## K. A one-parameter family of ties (one feature between two kernels): the margin against the critical eps at fixed geometry\n')
        L.append('Feature of amplitude 0.8 at 0.05 + delta (between the kernels at 0 and 0.1 of the 21), delta small; the path enters one kernel and the other is the runner-up. For each delta the critical eps is '
                 'found by bisection on log10 eps (flip: the formed set or the entry order differs from the path\'s); "model eps*" is the overtaking model\'s. A pure sqrt(eps) law at fixed geometry means '
                 'a slope 1/2 of log margin on log eps*.\n')
        for name, rows in Kf.items():
            ok = [r for r in rows if r['flip_hi'] and not r['flip_lo']]
            L.append(f'{name}:\n')
            L.append(table(['delta', 'margin', 'eps* bracket (from bisection)', 'geometric mean', 'model eps*', 'model/obs'],
                           [[g(r['delta']), g(r['margin']), f"{g(10 ** r['log_eps_lo'])} - {g(10 ** r['log_eps_hi'])}" if (r['flip_hi'] and not r['flip_lo']) else ('flips at both ends' if r['flip_lo'] else 'never flips'),
                             g(10 ** (0.5 * (r['log_eps_lo'] + r['log_eps_hi']))) if (r['flip_hi'] and not r['flip_lo']) else 'n/a', g(r['pred']),
                             g(r['pred'] / 10 ** (0.5 * (r['log_eps_lo'] + r['log_eps_hi']))) if (r['flip_hi'] and not r['flip_lo']) else 'n/a'] for r in rows]))
            if len(ok) >= 3:
                le = np.array([0.5 * (r['log_eps_lo'] + r['log_eps_hi']) for r in ok])
                lm = np.log10([r['margin'] for r in ok])
                sl = float(np.polyfit(le, lm, 1)[0])
                lp = np.log10([r['pred'] for r in ok])
                slp = float(np.polyfit(lp, lm, 1)[0])
                L.append(f"Slope of log10 margin on log10 eps* (observed): {g(sl, '.3f')} over {len(ok)} members; the same slope with the model's eps*: {g(slp, '.3f')}. (1/2: pure sqrt(eps).)\n")
    if Pr:
        recs = Pr['recs']
        L.append('## P. Sliding paths of the benchmark population\n')
        L.append(f"{len(recs)} one-dimensional targets of `exp_transit_certificate` (families A, B, D), equal masses; per target the stagewise drives and centre-slopes of the dictionary at the starting centres under "
                 f"each block's residual, the duals at the stages' ends, the runner-up margins. Failed or with exits: {sum(1 for r in recs if r['failed'] or r['exits'])}; "
                 f"median stage count {int(np.median([r['G'].shape[0] for r in recs]))}.\n")
    if F:
        L.append('## F. The 200 targets of Part C: what the theorem predicts, what the overtaking model adds\n')
        L.append(f"Matched to the population records: {F['n_matched']} of 200. For each cell: *scalar*: the target is predicted to mismatch iff its never-entering margin is below the exact formation margin mu_theta(eps) "
                 f"(theta = the target's formed threshold; depth 3: d_theta/wall, the paper's unit: its Lemma 2.1 distance); *overtaking*: the model of Section E applied to every (stage, competitor) of the path "
                 f"with mu_c = 12 and the readout's constant a_ro (no fitted parameter): a flip is predicted iff some competitor's overtaking margin reaches its margin at the stage's end. "
                 f"'Mismatch' in the first column is the paper's (the formed set differs from the path's); 'flip' also counts an entry order different from the path's. For the paper's unit mu_theta is negative when theta exceeds the coefficient that the unit carries at its wall (0.0079): such a unit is formed only after its dual has passed the wall, and the scalar prediction is empty for those targets.\n")
        rows = []
        for c in F['cells']:
            rows.append([c['name'], g(c['eps']), c['n_set_mismatch'], c['n_any_flip'], f"{g(c['mu_theta_min'])} - {g(c['mu_theta_max'])}", c['n_pred_scalar'], conf_str(c['conf_scalar']),
                         c['n_pred_drift'], conf_str(c['conf_drift_any']), conf_str(c['conf_drift_set'])])
        L.append(table(['readout', 'eps', 'set mismatches', 'flips', 'mu_theta (range over targets)', 'scalar: predicted', 'scalar vs set mismatch', 'overtaking: predicted', 'overtaking vs flip', 'overtaking vs set mismatch'], rows))
        for c in F['cells']:
            ps_, os_ = c['pred_scalar_margins'], c['obs_set_margins']
            L.append(f"{c['name']}, eps = {g(c['eps'])}: the scalar theorem predicts a mismatch for {len(ps_)} targets, whose never-entering margins are {', '.join(g(float(v), '.5f') for v in ps_[:40])}; "
                     f"the {len(os_)} observed set mismatches have margins {', '.join(g(float(v), '.5f') for v in os_[:40])}.\n")
        L.append('Largest runner-up margin among the mismatches (the statistic of the paper), observed and predicted by the overtaking model, and the exponent between 1e-6 and 1e-8:\n')
        rows = []
        by = {(c['name'], c['eps']): c for c in F['cells']}
        for name in ('D3', 'P0'):
            c6, c8 = by[(name, 1e-6)], by[(name, 1e-8)]
            rows.append([name, g(c6['max_margin_obs'], '.4f'), g(c8['max_margin_obs'], '.4f'), g(local_exponent(1e-6, 1e-8, c6['max_margin_obs'], c8['max_margin_obs']), '.3f'),
                         g(c6['max_margin_pred'], '.4f'), g(c8['max_margin_pred'], '.4f'), g(local_exponent(1e-6, 1e-8, c6['max_margin_pred'], c8['max_margin_pred']), '.3f'),
                         c6['n_set_mismatch'], c8['n_set_mismatch'], g(local_exponent(1e-6, 1e-8, c6['n_set_mismatch'], c8['n_set_mismatch']), '.3f')])
        L.append(table(['readout', 'max margin obs 1e-6', 'obs 1e-8', 'exponent obs', 'max margin pred 1e-6', 'pred 1e-8', 'exponent pred', 'mismatches 1e-6', '1e-8', 'count exponent'], rows))
        L.append('The theorem-level bound (a) on the targets: (max dormant coefficient of the run)/(the exact dormant supremum at the path\'s never-entering margin), over the targets with mu_never >= 100 sqrt(eps) '
                 'whose run formed the path\'s number of units (the actual dual differs from the path\'s by O(sqrt(eps)), so the ratio is 1 up to ~3 dev sqrt(eps)/mu):\n')
        rows = []
        for c in F['cells']:
            r = c['ratios']
            rows.append([c['name'], g(c['eps']), len(r), g(float(np.median(r))) if len(r) else 'n/a', g(float(np.percentile(r, 10))) if len(r) else 'n/a', g(float(np.percentile(r, 90))) if len(r) else 'n/a', g(float(r.max())) if len(r) else 'n/a'])
        L.append(table(['readout', 'eps', 'targets', 'median ratio', '10th pct', '90th pct', 'max'], rows))
        L.append('Depth 2 against the exact formation margin mu_theta(eps) = ln(2/(theta + sqrt(theta^2 + 4 eps^2)))/ln(1/eps) (theta per target) and the paper\'s rule ln(100)/ln(1/eps): predicted outcome "formed set = path\'s" iff mu_never > the margin, against the stored D2 runs:\n')
        L.append(table(['eps', 'D2 passes (observed)', 'exact margin', 'rule'], [[g(d['eps']), d['n_obs_pass'], conf_str(d['exact']), conf_str(d['rule'])] for d in F['d2']]))
        pd_ = F['d2_partD']
        L.append(f"Part D of the stored runs (D2 on the design and 20 targets at eps = 1e-6 ... 1e-144): outcome predicted correctly by the exact margin in {pd_['exact_agree']} of {pd_['total']} runs, by the rule in {pd_['rule_agree']}.\n")
        ds = F['design']
        L.append(f"The design (margin {g(ds['mu'])}, theta = {g(ds['theta'])}): the scalar theorem's \"path reached iff mu > mu_theta(readout, eps)\" against the stored runs (set equal to the path's):\n")
        _design_table = table(['readout', 'eps', 'mu_theta', 'predicted pass', 'observed pass', 'formed units', 'agree'], [[r['readout'], g(r['eps']), g(r['mu_theta']), r['pred_pass'], r['obs_pass'], r['n_formed'], 'yes' if (r['pred_pass'] == r['obs_pass'] or not r['ok']) else 'NO'] for r in ds['rows']])
        L.append(_design_table)
        dr = [r['ratio'] for r in ds['disp']]
        L.append('The displacement of the dormant centres on the design (stored runs; largest over the unformed units, in kernel widths) against the mass identity: by Theorem 2(iii) it is (centre rate x drive-slope/drive) x Psi_p((1 - mu) wall) with the same target-dependent factor for every depth and eps, Psi_p = eps b mu^(-2/(p-2))/(2 ell p) at the leading order; the ratio disp/Psi_p over the ' + str(len(dr)) + ' runs (D3, D4, D5 at eps = 1e-6, 1e-8, 1e-10; the failed run is excluded) has median ' + g(float(np.median(dr))) + ' and spread (max/min) ' + g(max(dr) / min(dr), '.2f') + ':\n')
        L.append(table(['readout', 'eps', 'dormant displacement', 'Psi_p at the design margin', 'ratio'], [[r['readout'], g(r['eps']), g(r['disp']), g(r['psi']), g(r['ratio'])] for r in ds['disp']]))
    if G:
        L.append('## G. Per-target critical eps (scan of `theory3d_scan.py`) against the overtaking model\n')
        for name in ('D3', 'P0'):
            rows_ = G.get(name, [])
            if not rows_:
                continue
            ok = [r for r in rows_ if not r['no_flip'] and not r['censored_low'] and np.isfinite(r['pred'])]
            if len(ok) >= 3:
                eo = np.array([np.sqrt(r['eps_hi'] * r['eps_lo']) for r in ok])
                ep = np.array([r['pred'] for r in ok])
                lr_ = np.log10(ep / eo)
                rk = lambda a: np.argsort(np.argsort(a)).astype(float)
                sp = float(np.corrcoef(rk(eo), rk(ep))[0, 1])
                slope = float(np.polyfit(np.log10(eo), np.log10(ep), 1)[0])
                mg = np.array([r['margin'] for r in ok])
                ex_obs = float(np.polyfit(np.log10(eo), np.log10(mg), 1)[0])
                s = (f"{name}: {len(rows_)} targets scanned (grid {g(min(SC.SPEC[name]['eps']))} ... {g(max(SC.SPEC[name]['eps']))}, targets of margin below {g(SC.SPEC[name]['margin_max'])}{(' and ' + str(sum(r['margin'] >= SC.SPEC[name]['margin_max'] for r in rows_)) + ' of larger margin from an earlier pass') if any(r['margin'] >= SC.SPEC[name]['margin_max'] for r in rows_) else ''}); {sum(r['no_flip'] for r in rows_)} do not flip even at the largest eps of the grid (critical eps above it), {sum(r['censored_low'] for r in rows_)} already flip at the smallest (critical eps below it; censored); "
                     f"{len(ok)} with a bracketed critical eps. log10(predicted/observed critical eps): median {g(float(np.median(lr_)), '.2f')}, quartiles {g(float(np.percentile(lr_, 25)), '.2f')} / {g(float(np.percentile(lr_, 75)), '.2f')}; "
                     f"Spearman correlation of the two critical eps {g(sp, '.3f')}; slope of log10(pred) on log10(obs) {g(slope, '.2f')}; slope of log10(margin) on log10(observed critical eps) {g(ex_obs, '.2f')} (1/2 would be a pure sqrt(eps) law, 1/4 the paper's pooled exponent). ")
                L.append(s + '\n')
                a_ = np.log10(mg)
                b_ = a_ - 0.5 * np.log10(eo)                                   # log10 c, c = margin/sqrt(observed critical eps)
                va, vb = float(np.var(a_)), float(np.var(b_))
                cab = float(np.mean((a_ - a_.mean()) * (b_ - b_.mean())))
                s_id = (va - cab) / (2.0 * (va + vb - 2.0 * cab))
                L.append(f"Decomposition: log10 margin = log10 c + (1/2) log10 eps*, c = margin/sqrt(eps*) the observed constant of the target. Var(log10 margin) = {g(va, '.3f')}, Var(log10 c) = {g(vb, '.3f')}, correlation "
                         f"{g(cab / np.sqrt(va * vb), '.2f')}: the regression slope of log margin on log eps* is (Var a - Cov(a, b))/(2 (Var a + Var b - 2 Cov(a, b))) = {g(s_id, '.3f')} (identity; if log c were independent of "
                         f"the margin it would be (1/2) Var a/(Var a + Var b) = {g(0.5 * va / (va + vb), '.3f')}): the spread of c (a factor {g(10 ** np.sqrt(vb), '.1f')} per standard deviation) attenuates the slope from 1/2 to the observed value.\n")
            else:
                L.append(f"{name}: {len(rows_)} targets scanned; too few bracketed critical eps for statistics.\n")
        L.append('The columns of the table below are, per scanned target: the runner-up margin, the observed bracket (flip at eps_hi, none at eps_lo), the model\'s critical eps and the ratio of the margin to sqrt(observed eps).\n')
        for name in ('D3', 'P0'):
            rows_ = [r for r in G.get(name, []) if not r['no_flip'] and not r['censored_low']]
            if rows_:
                rows_ = sorted(rows_, key=lambda r: -r['margin'])[:15]
                L.append(f'{name}, the 15 flipping targets of largest margin:\n')
                L.append(table(['target', 'margin', 'eps_hi (flip)', 'eps_lo (none)', 'model eps*', 'margin/sqrt(eps_hi)'],
                               [[r['i'], g(r['margin'], '.4f'), g(r['eps_hi']), g(r['eps_lo']), g(r['pred']), g(r['margin'] / np.sqrt(r['eps_hi']), '.1f')] for r in rows_]))
    if H:
        L.append('## H. The population: counts and the largest flipped margin of the model as functions of eps, and the exponent\n')
        L.append(f"{H['n_good']} of {H['n']} population targets with a usable path (no exit, no solver failure). The model's count of flips N(eps) and the largest runner-up margin among them, for the two readouts; "
                 f"in the model N(eps) = #(targets with rho_t < sqrt(a_ro eps)) to leading order, rho_t = min over (stage, competitor) of x/sqrt(q), so the local exponent of N in eps is half the local exponent of the distribution function of rho_t.\n")
        for name in ('D3', 'P0'):
            b = H['by_name'][name]
            eps = H['eps']
            L.append(f'{name}:\n')
            L.append(table(['eps', 'flips N(eps)', 'largest margin among flips', 'its value / sqrt(eps)', '95th pct margin'], [[g(e), int(n), g(m, '.4f'), g(m / np.sqrt(e), '.1f') if np.isfinite(m) else 'n/a', g(q, '.4f')] for e, n, m, q in zip(eps, b['count'], b['maxmargin'], b['q95'])]))
            L.append(f"Bootstrap ({BOOT_N} samples of {BOOT_SIZE} targets from the population, the model's flips at 1e-6 and 1e-8, {b['boot_n']} samples with flips at both): exponent between 1e-6 and 1e-8 of the count N(eps): 5/25/50/75/95th percentiles {', '.join(g(float(v), '.2f') for v in b['boot_count_exp'])}; of the largest margin among the flips: {', '.join(g(float(v), '.2f') for v in b['boot_max_exp'])}.\n")
            L.append(f"The constant of the flip margin c = sqrt(a_ro mu_c) |h/g| (a competitor flips when its margin is below c sqrt(eps)), over the {b['c_n']} (stage, competitor) pairs with margin below 0.1: 10/50/90/99th percentiles {', '.join(g(float(v), '.1f') for v in b['c_pct'])}; fraction above 30: {g(b['c_frac30'], '.3f')}.\n")
            if name == 'D3':
                rho = b['rho']
                rho = rho[np.isfinite(rho)]
                rows = []
                for lo, hi in ((1e-2, 1e-3), (1e-3, 1e-4), (1e-4, 1e-5)):
                    a, b_ = int((rho < lo).sum()), int((rho < hi).sum())
                    rows.append([f'{g(lo)} -> {g(hi)}', a, b_, g(np.log(a / b_) / np.log(lo / hi), '.2f') if b_ > 0 and a > 0 else 'n/a'])
                L.append('Distribution function of rho_t = min over (stage, competitor) of x |g|/(sqrt(mu_c)|h|) (the same for both readouts; the flip condition is rho_t < sqrt(a_ro eps)): count of targets below, local index between consecutive decades '
                         '(1 would mean a pure sqrt(eps) law for the counts; the first row is saturated: every target has rho_t below 0.01):\n')
                L.append(table(['rho range', 'count below the larger', 'count below the smaller', 'local index'], rows))
    L.append('## Status of the statements (the note has the proofs)\n')
    L.append('* PROVED (exact, scalar): the closed forms of the unit of depth p (Theorem 1); the dormant coefficient bound |w| <= (eps b)^(p/2) mu^(-p/(p-2)), sharp, and the same in terms of a drive bound (Theorem 2); the distance form and the formation '
             'criterion "formed iff the dual is within d_theta of the wall", d_theta ~ theta^(-(p-2)/p)/(ell p (p-2)) (Theorem 2); depth 2: w = 2 eps b sinh(4 ell eta), formation margin exactly '
             'ln(2/(theta + sqrt(theta^2 + 4 (eps b)^2)))/ln(1/(eps b)) and the bracket by ln(1/theta)/ln(1/(eps b)) (Theorem 3); depth 1: w = ell eta, the fixed-centre flow and its formation times (Theorem 4). Verified here: A (ODE, 1e-10), B (ball arithmetic), C (the repository\'s flows with moving centres), D.\n'
             '* PROVED up to the cited abstract steps (those of the paper\'s fixed-centre theorem): the fixed-centre limits (Corollary 5): weighted l1 inverse scale space for depth p >= 3, unweighted l1 for depth 2 (at a logarithmic rate), the l2 flow for depth 1; the epi-convergence of the potentials is proved. Verified in L.\n'
             '* PROVED (exact, model): the first integral of the frozen-environment sliding unit (E) and the proportionality of the dormant centre\'s displacement to the mass identity Psi_p (checked on the stored runs of D3, D4, D5: ratio constant to 6%).\n'
             '* COMPUTER-ASSISTED (a first-order model without a fitted parameter): the competitor\'s overtaking margin c sqrt(eps), c = sqrt(mu_c a_ro)|h/g|, its agreement with the 200 targets (F), with the per-target critical eps (G) and with the tie family (K, slope 0.53 at fixed geometry), and the explanation of the pooled exponent '
             '(the attenuation of the slope by the spread of c; the bootstrap of H). The explanation is not a theorem.\n'
             '* NOT proved: that the actual duals of the network follow the path\'s to the accuracy that the thresholds need (the hypothesis of the dormant bound is on the actual dual); the sliding-limit theorem for depth p != the paper\'s unit; the layer analysis for D_p; the rate of depth 2\'s dual error (observed to fall faster than 1/ln(1/eps) on the design). '
             'The scalar theorem alone does not predict the paper\'s mismatch margins (it accounts for 12 of the 24 depth-3 set mismatches at 1e-6 and 2 of the 34 of the paper\'s unit); the dual deviation does.\n')
    return '\n'.join(L)


def main():
    os.makedirs(OUT, exist_ok=True)
    args = sys.argv[1:]
    if '--from-raw' in args:
        res = load(RAW)
        scan = SC.load_scan(SCAN_RAW)
        open(MD, 'w').write(report(*[res.get(k) for k in PARTS], scan=scan))
        return
    res = load(RAW)
    only = set(args[args.index('--only') + 1].split(',')) if '--only' in args else set(PARTS)
    analysis = '--analysis' in args
    t0 = time.time()
    warnings.filterwarnings('ignore')
    fns = dict(A=part_A, B=part_B, C=part_C, D=part_D, E=part_E, L=part_L, K=part_K, P=part_P)
    if '--run' in args:
        for k in ('A', 'B', 'C', 'D', 'E', 'L', 'K', 'P'):
            if k in only:
                res[k] = fns[k]()
                save(RAW, res)
                print(f'[{time.time() - t0:.0f}s] {k}', flush=True)
    if '--run' in args or analysis:
        scan = SC.load_scan(SCAN_RAW)
        if res.get('P'):
            res['F'] = part_F(res['P'])
            res['G'] = part_G(res['P'], scan)
            res['H'] = part_H(res['P'])
            save(RAW, res)
            print(f'[{time.time() - t0:.0f}s] F G H', flush=True)
    open(MD, 'w').write(report(*[res.get(k) for k in PARTS], scan=SC.load_scan(SCAN_RAW)))


if __name__ == '__main__':
    main()
