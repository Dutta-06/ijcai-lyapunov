"""
Proposition D‴.47: the free law's branches — an exponentiated gradient with frozen signs; the limit for every start is the Bregman
projection onto the sign-constrained least-squares fits, the wrong-signed particles going extinct; positive kernel dictionaries
====================================================================================================================================

The free mass–amplitude law of Proposition D‴.45 with particles on branches (Q_i < −c: a start with β_i(0) < 2c sinh²(α_i(0)/a_max)),
represented without cancellation near their branch points (`mirror_free.ParticleSet`).

A. Identities: the slope dη/dlog|w| = a/(γa² + Gβa′²) at random states (finite differences); the temperature 1/(γa₀) at branch points;
   B(w̄, w) → β/(γm) as w̄ → 0 on a branch; the fraction of small random starts on branches; a ± pair's coefficient against
   2u sinh(γa₀Δ) and its potential against the hyperbolic entropy; the critical curve's power law.
B. The limit on random Gaussian designs (6 × 12) from small random starts (all on branches, random or positive signs; and a mixed start):
   the flow to t = 3·10⁷ against the Bregman projection onto S (computed independently, `sign_ls_set` and `projection`), J − min_K J,
   the extinct masses and their decay rates against γa₀|g_i(∞)|, B(w̄, w(t)) nonincreasing.
C. The rich limits on one design (8 × 20): branch starts at masses ε → the sign-constrained basis pursuit weighted by 1/|a(0)|; starts at
   α = 0 → basis pursuit; a mixed start (half at α = 0).
D. A pair per atom: the rich limit of ± pairs → basis pursuit weighted by 1/a₀.
E. Positive kernel dictionaries (Gaussian atoms on a grid, few samples, a positive 3-sparse target with two adjacent atoms): the
   nonnegative certificate's margin; basis pursuit with equal and with column-norm weights; the loop from positive branches at four scales
   and from α = 0 at two, to t = 10⁹; the other atoms' mass against the certificate's bound.
F. The real loop on positive branches: V̇ = −k_ee² with w* vanishing on unused atoms (random states); a run (twelve atoms, a positive
   target on three): V, the error, the unused atoms' masses by quarters.
G. (round 144) The mixed start of C resolved on both scales. The branch atoms' own vertex has |S| = 3 < N atoms, so the dual directions
   off range Φ_S couple only to the free and off-support coefficients, whose contribution to Φw lies below the residual at which
   `projection` stops once ε ≤ 10⁻⁸. Here the dual equations are solved on range Φ_S and, scaled separately, on its complement (free
   particles exact for ε ≥ 10⁻¹², by their linear core below), down to ε = 10⁻¹²⁸, against the leading order z = log(1/ε)ζ*, ζ* the
   minimiser of Σ_free(Φᵀζ)_j² subject to Φ_Sᵀζ = ρ_S·s_S and s_j(Φᵀζ)_j ≤ −ρ_j on the off-support branch atoms (ρ = 1/(γ|a(0)|)).

Output: results/chaos/mirror_branch.md and mirror_branch_raw.npz (`--from-raw` regenerates the md; `--only A,B,C,D,E,F,G`).
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
from scipy.integrate import solve_ivp, quad

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
from src.core import mirror_free as MF
import exp_cohesion_leak as E

OUT = os.environ.get('MB_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
LAW_B = dict(m=12.0, G=1.0, gamma=0.1, a_max=1.0)                 # c = 0.025
LAW_E = dict(m=12.0, G=0.25, gamma=1.0, a_max=1.0)                # c = 1
SEEDS_B, NB, DB, T_B = range(8), 6, 12, 3e7
ND_B, DD_B = 5, 4                                                 # B, duplicated: 5 samples, 4 atoms twice (round 36)
NC, DC, EPS_BR, EPS_FREE = 8, 20, (1e-4, 1e-8, 1e-16, 1e-32, 1e-64, 1e-128), (1e-2, 1e-4, 1e-6, 1e-8, 1e-10)
U_PAIR = (1e-4, 1e-8, 1e-16, 1e-32, 1e-64, 1e-128, 1e-256)
CASES_E = ((21, 0.15, 7), (41, 0.15, 7), (41, 0.1, 11))
STARTS_E = (('branch', 1.0, 1e-6), ('branch', 1.0, 1e-2), ('branch', 1.0, 1.0), ('branch', 2.0, 10.0), ('alpha=0', 0.0, 1e-6),
            ('alpha=0', 0.0, 1e2))
T_E = 1e9
Q_E, GEN_START_E = (0.25, 0.5, 0.75, 0.9), 1                         # E: margin quantiles; the start (branch, 1, 10⁻²)
LAW_F = dict(m=12.0, G=1.0, gamma=1.0, a_max=1.0)                # c = 0.25
KE_F, T_F = 5.0, 4000.0
EPS_G = (1e-4, 1e-8, 1e-12, 1e-16, 1e-24, 1e-32, 1e-48, 1e-64, 1e-96, 1e-128)
EPS_G_EXACT = 1e-12                                              # G: free particles exact down to here, their linear core below
XR = (lambda t: np.sin(t) + 0.4 * np.sin(2.7 * t))
XRD = (lambda t: np.cos(t) + 1.08 * np.cos(2.7 * t))


def part_A():
    L = MF.FreeLaw(**LAW_B)
    rng = np.random.default_rng(0)
    # A1: the slope at random states, free and on branches (finite differences along each particle's state)
    al_i = rng.uniform(-1.5, 1.5, 40); b_i = 10 ** rng.uniform(-6, 1, 40)
    P = MF.ParticleSet(L, al_i, b_i)
    x = P.x0.copy(); h = 1e-6
    err = 0.0
    for j in range(P.d):
        xp, xm = x.copy(), x.copy(); xp[j] += h; xm[j] -= h
        num = (P.eta(xp)[j] - P.eta(xm)[j]) / (np.log(abs(P.w(xp)[j])) - np.log(abs(P.w(xm)[j])))
        err = max(err, abs(num / P.temperature(x)[j] - 1))
    # A2: the temperature at branch points
    tmp = []
    for a0 in (0.05, 0.3, 1.0, 2.0):
        Pb = MF.ParticleSet(L, np.array([a0 + 1e-30]), np.array([2 * L.c * np.sinh(a0) * np.sinh(1e-30) * 2]))
        xb = np.array([np.log(1e-40)])
        Pb.al0[:] = a0
        tmp.append(dict(a0=a0, T=float(Pb.temperature(xb)[0]), ref=1 / (L.gamma * L.A * np.tanh(a0 / L.A))))
    # A3: B(w̄, w) → β/(γm) as w̄ → 0 on a branch
    Pc = MF.ParticleSet(L, np.array([0.4]), np.array([1e-3]))
    xw = np.array([np.log(0.05)])
    b0 = []
    for dbar in (1e-6, 1e-12, 1e-30, 1e-100):
        xbar = np.array([np.log(dbar)])
        Bv = Pc.bdiv(Pc.w(xbar), xbar, xw)
        b0.append(dict(dbar=dbar, wbar=float(Pc.w(xbar)[0]), B=Bv, ref=float(Pc.beta(xw)[0] / (L.gamma * L.m))))
    # A4: the fraction of small random starts on branches
    frac = []
    for b in (1e-3, 1e-5, 1e-7, 1e-9):
        al = rng.uniform(-0.3, 0.3, 20000)
        frac.append(dict(b=b, frac=float(np.mean(b < 2 * L.c * np.sinh(al / L.A) ** 2))))
    # A5: a ± pair: W(Δ) against 2u sinh(γa₀Δ), the potential against the hyperbolic entropy (relative errors by |W|/u)
    al0, bb = 0.3, 1e-7
    Pp = MF.ParticleSet(L, np.array([al0, -al0]), np.array([bb, bb]))
    e0 = Pp.eta(Pp.x0); u = float(Pp.w(Pp.x0)[0]); a0 = L.A * np.tanh(Pp.al0[0] / L.A)
    pair = []
    for D in (1.0, 10.0, 30.0, 60.0, 90.0, 120.0):
        xs = Pp.x_of_eta(e0 + D)
        W = float(Pp.w(xs).sum())
        ref = 2 * u * np.sinh(L.gamma * a0 * D)
        pair.append(dict(D=D, ratio=W / u, rel=abs(W / ref - 1)))
    # A6: the critical curve Q = −c: η·w^{1/3} → −(a_max²/(2Gc))(2c/(ma_max²))^{1/3} as α → 0⁺ (η from α to 1 by quadrature, shifted)
    A, G, c, m = L.A, L.G, L.c, L.m
    beta_c = lambda a_: 2 * c * np.sinh(a_ / A) ** 2
    integrand = lambda a_: 1 / (G * beta_c(a_) * (1 - np.tanh(a_ / A) ** 2))
    crit = []
    for a_ in (1e-2, 1e-3, 1e-4):
        eta_ = -quad(integrand, a_, 1.0, limit=200, epsabs=0, epsrel=1e-12)[0]
        w_ = beta_c(a_) * A * np.tanh(a_ / A) / m
        crit.append(dict(alpha=a_, val=eta_ * w_ ** (1 / 3), ref=-(A ** 2 / (2 * G * c)) * (2 * c / (m * A ** 2)) ** (1 / 3)))
    return dict(slope=err, temp=tmp, b0=b0, frac=frac, pair=pair, u=u, a0=a0, crit=crit)


def gauss_design(seed, N, d):
    return np.random.default_rng(seed).standard_normal((N, d))


def flow_B(P, Phi, y, T, n=60):
    ts = np.concatenate([[0], np.geomspace(1e-1, T, n)])
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        sol = solve_ivp(P.rhs(Phi, y), (0, T), P.x0, method='LSODA', rtol=1e-10, atol=1e-12, t_eval=ts)
    return sol


def b_case(seed, kind):
    L = MF.FreeLaw(**LAW_B)
    rng = np.random.default_rng(100 + seed)
    if kind == 'duplicated':
        # round 36: each atom twice (the least-squares set S not a point), one atom's frozen sign against its unconstrained least-squares
        # coefficient (Z ≠ ∅: the fit inexact), negative branches alive — the Bregman selection with extinction
        r2 = np.random.default_rng(950 + seed)
        Phi0 = r2.standard_normal((ND_B, DD_B)); y = r2.standard_normal(ND_B)
        s0 = np.sign(np.linalg.lstsq(Phi0, y, rcond=None)[0]); s0[-1] = -s0[-1]
        Phi = np.hstack([Phi0, Phi0])
        al_i, b_i = np.r_[s0, s0] * r2.uniform(0.4, 0.8, 2 * DD_B), np.full(2 * DD_B, 1e-2)
    else:
        Phi = gauss_design(seed, NB, DB)
        ws = np.zeros(DB); ws[[0, 3, 5, 8]] = [0.9, -0.7, 0.6, -0.5]; y = Phi @ ws
        amp = rng.uniform(0.2, 0.4, DB)
        if kind == 'positive':
            al_i, b_i = amp, np.full(DB, 1e-3)
        elif kind == 'random signs':
            al_i, b_i = amp * rng.choice([-1, 1], DB), np.full(DB, 1e-3)
        else:                                                                # mixed: half the particles at α = 0 (free signs)
            al_i = amp * rng.choice([-1, 1], DB); al_i[::2] = 0.0; b_i = np.full(DB, 1e-3)
    N, d = Phi.shape
    P = MF.ParticleSet(L, al_i, b_i)
    eta0 = P.eta(P.x0)
    yh, Z, _ = MF.sign_ls_set(Phi, y, P.s)
    keep = np.setdiff1d(np.arange(d), Z)
    wbar, it, res = MF.projection(P, Phi, yh, eta0, keep)
    sol = flow_B(P, Phi, y, T_B)
    W = np.array([P.w(sol.y[:, k]) for k in range(sol.y.shape[1])])
    Bt = np.array([P.beta(sol.y[:, k]) for k in range(sol.y.shape[1])])
    # B(w̄, w(t)) along the flow (w̄ = the projection; its zeros on Z read as β/(γm)); w̄'s own state off Z by root finding
    xb = P.state_of_w(wbar, keep)
    Bv = np.array([P.bdiv(wbar, xb, sol.y[:, k]) for k in range(sol.y.shape[1])])
    incr = float(max(0.0, np.max(np.diff(Bv)) / max(Bv[0], 1e-300)))
    J = 0.5 * N ** -1 * L.m * np.sum((W @ Phi.T - y) ** 2, axis=1)
    Jstar = 0.5 * L.m / N * np.sum((yh - y) ** 2)
    g_inf = -(Phi.T @ (Phi @ wbar - y)) / N
    a0 = L.A * np.tanh(P.al0 / L.A)
    rates = []
    for i in Z:
        sel = (Bt[:, i] > 1e-250) & (Bt[:, i] < 1e-12)
        if sel.sum() >= 3:
            sl = np.polyfit(sol.t[sel], np.log(Bt[sel, i]), 1)[0]
            rates.append(dict(i=int(i), fit=float(-sl), ref=float(L.gamma * a0[i] * abs(g_inf[i])), sg=float(P.s[i] * g_inf[i])))
    return dict(seed=seed, kind=kind, nbr=int(P.br.sum()), yres=float(np.linalg.norm(y - yh)), Z=[int(i) for i in Z], it=it, res=res,
                status=int(sol.status), dist=[float(np.abs(W[k] - wbar).max()) for k in range(W.shape[0])], t=sol.t.tolist(),
                dJ=[float(J[k] - Jstar) for k in range(W.shape[0])], bZ=[float(Bt[k, Z].max()) if len(Z) else 0.0 for k in range(W.shape[0])],
                Bincr=incr, rates=rates, nneg=int(np.sum(wbar[keep] < 0)),
                split=float(np.abs(wbar - np.r_[(wbar[:d // 2] + wbar[d // 2:]) / 2, (wbar[:d // 2] + wbar[d // 2:]) / 2]).max())
                if kind == 'duplicated' else float('nan'))


def part_B():
    return [b_case(s, kind) for kind in ('positive', 'random signs', 'mixed', 'duplicated') for s in SEEDS_B]


def part_C():
    L = MF.FreeLaw(**LAW_B)
    rng = np.random.default_rng(7)
    Phi = rng.standard_normal((NC, DC))
    ws = np.zeros(DC); ws[[1, 5, 9]] = [0.8, -0.6, 0.5]; y = Phi @ ws
    sg = rng.choice([-1.0, 1.0], DC); sg[[1, 5, 9]] = np.sign(ws[[1, 5, 9]])
    amp0 = rng.uniform(0.2, 0.6, DC)
    a_init = L.A * np.tanh(amp0 / L.A)
    w_scbp = MF.weighted_bp(Phi, y, 1 / a_init, signs=sg)
    w_bp = MF.weighted_bp(Phi, y, np.ones(DC))
    S_sc = np.flatnonzero(np.abs(w_scbp) > 1e-9); S_bp = np.flatnonzero(np.abs(w_bp) > 1e-9)
    off_sc = np.setdiff1d(np.arange(DC), S_sc); off_bp = np.setdiff1d(np.arange(DC), S_bp)
    mu_sc = MF.certificate_slacks(Phi, w_scbp, S_sc, 1 / (L.gamma * a_init), frozen=sg)[0] if S_sc.size == NC else np.full(DC, np.nan)
    mu_bp, xi_bp = MF.certificate_slacks(Phi, w_bp, S_bp, np.ones(DC)) if S_bp.size == NC else (np.full(DC, np.nan), None)
    br = []
    for eps in EPS_BR:
        P = MF.ParticleSet(L, sg * amp0, np.full(DC, eps))
        wb, it, res = MF.projection(P, Phi, y, P.eta(P.x0), exact=True)
        w1 = MF.rich_leading_order(P, Phi, y, w_scbp, S_sc) if S_sc.size == NC else np.full(DC, np.nan)
        br.append(dict(eps=eps, all_branch=bool(P.br.all()), it=it, res=res, d=float(np.abs(wb - w_scbp).max()),
                       d1=float(np.abs(w1 - w_scbp).max()), agree=float(np.abs(wb - w1).max())))
    fr = []
    for eps in EPS_FREE:
        P = MF.ParticleSet(L, np.zeros(DC), np.full(DC, eps))
        wb, it, res = MF.projection(P, Phi, y, P.eta(P.x0), exact=True)
        w1 = MF.rich_leading_order(P, Phi, y, w_bp, S_bp) if S_bp.size == NC else np.full(DC, np.nan)
        cf = np.sign(Phi.T @ xi_bp) * L.core_tail(eps, mu_bp)                  # the closed form off the support, the certificate's sign
        fr.append(dict(eps=eps, it=it, res=res, d=float(np.abs(wb - w_bp).max()), d1=float(np.abs(w1 - w_bp).max()),
                       agree=float(np.abs(wb - w1).max()), cf=float(np.max(np.abs(wb[off_bp] - cf[off_bp]) / np.abs(cf[off_bp])))))
    # a mixed start: the even atoms at α = 0 (free), the odd ones on branches with the truth's signs where it has them
    mx = []
    odd = np.arange(DC) % 2 == 1
    w_scbp_odd = np.zeros(DC)
    sub = MF.weighted_bp(Phi[:, odd], y, 1 / a_init[odd], signs=sg[odd])
    w_scbp_odd[odd] = sub
    for eps in (1e-4, 1e-8, 1e-12):
        al_i = np.where(odd, sg * amp0, 0.0)
        P = MF.ParticleSet(L, al_i, np.full(DC, eps))
        wb, it, res = MF.projection(P, Phi, y, P.eta(P.x0), exact=True)
        mx.append(dict(eps=eps, it=it, res=res, free=float(np.abs(wb[~odd]).sum()), d=float(np.abs(wb - w_scbp_odd).max())))
    # round 36: the odd atoms' frozen signs against the truth's on its atoms — they cannot fit y; the free atoms then enter on the power
    # scale: the limit is basis pursuit over the free atoms with the branch atoms sign-constrained and (to leading order) free of charge,
    # their own weights 1/(γa₀) breaking ties at the next order
    sg2 = sg.copy(); sg2[[1, 5, 9]] = -np.sign(ws[[1, 5, 9]])
    sub2 = MF.weighted_bp(Phi[:, odd], y, 1 / a_init[odd], signs=sg2[odd])
    w_mix2 = MF.weighted_bp(Phi, y, np.where(odd, 1e-6 / (L.gamma * a_init), 1.0), signs=np.where(odd, sg2, 0.0))
    mx2 = []
    for eps in (1e-4, 1e-8, 1e-12):
        P = MF.ParticleSet(L, np.where(odd, sg2 * amp0, 0.0), np.full(DC, eps))
        wb, it, res = MF.projection(P, Phi, y, P.eta(P.x0), exact=True)
        mx2.append(dict(eps=eps, it=it, res=res, free=float(np.abs(wb[~odd]).sum()), d=float(np.abs(wb - w_mix2).max())))
    tol_r = 1e-8 * (1 + np.linalg.norm(y))
    nfail = sum(1 for r in br + fr + mx + mx2 if not r['res'] < tol_r)
    return dict(scbp_supp=[int(i) for i in S_sc], scbp_is_truth=bool(np.abs(w_scbp - ws).max() < 1e-8),
                bp_supp=[int(i) for i in S_bp], br=br, fr=fr, mx=mx, odd_feasible=bool(np.all(np.isfinite(sub))),
                mx2=mx2, odd2_feasible=bool(np.all(np.isfinite(sub2))), mix2_free=float(np.abs(w_mix2[~odd]).sum()), nfail=nfail, tol_r=tol_r,
                mu_sc=np.sort(mu_sc[off_sc]).tolist(), mu_bp=np.sort(mu_bp[off_bp]).tolist())


def part_D():
    L = MF.FreeLaw(**LAW_B)
    rng = np.random.default_rng(7)
    Phi = rng.standard_normal((NC, DC))
    ws = np.zeros(DC); ws[[1, 5, 9]] = [0.8, -0.6, 0.5]; y = Phi @ ws
    al = 0.4
    Phi2 = np.hstack([Phi, Phi])
    w_bp = MF.weighted_bp(Phi, y, np.ones(DC))                                 # equal amplitudes: the weights 1/a₀ are equal
    S_bp = np.flatnonzero(np.abs(w_bp) > 1e-9)
    off_bp = np.setdiff1d(np.arange(DC), S_bp)
    mu_bp = MF.certificate_slacks(Phi, w_bp, S_bp, np.ones(DC))[0] if S_bp.size == NC else np.full(DC, np.nan)
    S2 = np.where(w_bp[S_bp] > 0, S_bp, DC + S_bp)                              # the member of each support atom with its sign
    rows = []
    for u in U_PAIR:
        b = u * L.m / (L.A * np.tanh(al / L.A))                               # w(0) = u for each member
        P = MF.ParticleSet(L, np.r_[np.full(DC, al), np.full(DC, -al)], np.full(2 * DC, b))
        assert P.br.all(), u
        wb, it, res = MF.projection(P, Phi2, y, P.eta(P.x0), exact=True)
        W = wb[:DC] + wb[DC:]
        a0 = L.A * np.tanh(P.al0[0] / L.A)
        wref = MF.weighted_bp(Phi, y, np.full(DC, 1 / a0))
        w2 = np.zeros(2 * DC); w2[S2] = wref[S_bp]
        w1 = MF.rich_leading_order(P, Phi2, y, w2, S2) if S_bp.size == NC else np.full(2 * DC, np.nan)
        W1 = w1[:DC] + w1[DC:]
        rows.append(dict(u=u, it=it, res=res, d=float(np.abs(W - wref).max()), d1=float(np.abs(W1 - wref).max()),
                         agree=float(np.abs(W - W1).max()),
                         both=int(np.sum((np.abs(wb[:DC]) > 1e-3 * np.abs(W).max()) & (np.abs(wb[DC:]) > 1e-3 * np.abs(W).max())))))
    tol_r = 1e-8 * (1 + np.linalg.norm(y))
    return dict(rows=rows, bp_supp=[int(i) for i in S_bp], mu_bp=np.sort(mu_bp[off_bp]).tolist(), tol_r=tol_r,
                nfail=sum(1 for r in rows if not r['res'] < tol_r))


def e_case(d, wid, N):
    L = MF.FreeLaw(**LAW_E)
    cen = np.linspace(-1, 1, d); xs = np.linspace(-1.1, 1.1, N)
    Phi = np.exp(-(xs[:, None] - cen[None, :]) ** 2 / (2 * wid ** 2))
    i0 = d // 4; supp = np.array([i0, i0 + 1, 3 * d // 4])
    ws = np.zeros(d); ws[supp] = [1.0, 0.7, 0.5]; y = Phi @ ws
    mu, v = MF.nonneg_certificate(Phi, supp)
    w_bp = MF.weighted_bp(Phi, y, np.ones(d)); w_bpn = MF.weighted_bp(Phi, y, np.linalg.norm(Phi, axis=0))
    rows = []
    for lab, al_v, b_v in STARTS_E:
        P = MF.ParticleSet(L, np.full(d, al_v), np.full(d, b_v))
        ts = np.concatenate([[0], np.geomspace(1.0, T_E, 36)])
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            sol = solve_ivp(P.rhs(Phi, y), (0, T_E), P.x0, method='LSODA', rtol=1e-10, atol=1e-12, t_eval=ts)
        Wt = np.array([P.w(sol.y[:, k]) for k in range(sol.y.shape[1])])
        off = np.setdiff1d(np.arange(d), supp)
        row = dict(lab=lab, al=al_v, b=b_v, branch=bool(P.br.all()), t=sol.t.tolist(), err=np.abs(Wt - ws).max(axis=1).tolist(),
                   offmass=np.abs(Wt[:, off]).sum(axis=1).tolist(), res=float(np.linalg.norm(Phi @ Wt[-1] - y)),
                   supp=[int(i) for i in np.flatnonzero(np.abs(Wt[-1]) > 1e-4 * np.abs(Wt[-1]).max())])
        if P.br.all() and v is not None:
            # the certificate's bound (|v|/μ)(2N·B(w*, w(0))/t)^{1/2}: B(w*, w(0)) with w*'s zeros read as β/(γm)
            xs_ = P.state_of_w(ws, supp)
            B0 = P.bdiv(ws, xs_, P.x0)
            row['bound'] = [float(np.linalg.norm(v) / mu * np.sqrt(2 * N * B0 / t)) if t > 0 else float('inf') for t in sol.t]
        rows.append(row)
    # round 36: the margin over every 3-support (the target's rank among them), and the loop from one start on supports at the margins'
    # quartiles (amplitudes 1.0, 0.7, 0.5) — the convergence's dependence on the margin
    from itertools import combinations
    trip = list(combinations(range(d), 3))
    mus = np.array([MF.nonneg_certificate(Phi, list(t_))[0] for t_ in trip]) + 0.0          # + 0.0: no signed zeros
    order = np.argsort(np.nan_to_num(mus, nan=-1.0))
    gen = []
    for q_ in Q_E:
        t_ = np.array(trip[order[int(q_ * (len(trip) - 1))]])
        wg = np.zeros(d); wg[t_] = [1.0, 0.7, 0.5]; yg = Phi @ wg
        lab, al_v, b_v = STARTS_E[GEN_START_E]
        P = MF.ParticleSet(L, np.full(d, al_v), np.full(d, b_v))
        ts = np.concatenate([[0], np.geomspace(1.0, T_E, 19)])
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            sol = solve_ivp(P.rhs(Phi, yg), (0, T_E), P.x0, method='LSODA', rtol=1e-10, atol=1e-12, t_eval=ts)
        Wt = np.array([P.w(sol.y[:, k]) for k in range(sol.y.shape[1])])
        offg = np.setdiff1d(np.arange(d), t_)
        gen.append(dict(q=q_, supp=t_.tolist(), mu=float(mus[order[int(q_ * (len(trip) - 1))]]), t=sol.t.tolist(),
                        err=np.abs(Wt - wg).max(axis=1).tolist(), offmass=np.abs(Wt[:, offg]).sum(axis=1).tolist(),
                        res=float(np.linalg.norm(Phi @ Wt[-1] - yg))))
    fin = mus[np.isfinite(mus)]
    return dict(d=d, wid=wid, N=N, supp=supp.tolist(), mu=float(mu), bp=float(np.abs(w_bp - ws).max()), bpn=float(np.abs(w_bpn - ws).max()),
                rows=rows, n_supp=len(trip), mu_q=np.quantile(fin, [0.1, 0.25, 0.5, 0.75, 0.9]).tolist(), mu_zero=float(np.mean(fin <= 1e-9)),
                mu_rank=float(np.mean(fin < mu)), gen=gen)


def part_E():
    return [e_case(*c) for c in CASES_E]


def part_F():
    L = MF.FreeLaw(**LAW_F)
    rng = np.random.default_rng(5)
    d = 12
    cen, wid = np.linspace(-1.4, 1.4, d), np.full(d, 0.25)
    wstar = np.zeros(d); wstar[[2, 6, 9]] = [0.12, 0.08, 0.1]
    P = MF.ParticleSet(L, np.full(d, 0.5), np.full(d, 1e-3))
    assert P.br.all()
    f = P.real_rhs(cen, wid, wstar, XR, XRD, KE_F)
    # the identity V̇ = −k_ee² at random states (w* on the curves of the used atoms, β/(γm) on the others)
    used = np.flatnonzero(wstar)
    xs_ = P.state_of_w(wstar, used)
    V = lambda t, z: 0.5 * (z[0] - XR(t)) ** 2 + P.bdiv(wstar, xs_, z[1:])
    worst = 0.0
    for k in range(12):
        t = rng.uniform(0, 50); z = np.r_[XR(t) + rng.normal(0, 0.3), P.x0 + rng.normal(0, 0.5, d)]
        dz = f(t, z); h = 1e-6
        num = (V(t + h, z + h * dz) - V(t - h, z - h * dz)) / (2 * h)
        e = z[0] - XR(t)
        worst = max(worst, abs(num + KE_F * e ** 2) / max(abs(KE_F * e ** 2), 1e-12))
    z0 = np.r_[XR(0.0), P.x0]
    ts = np.linspace(0, T_F, 8001)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        sol = solve_ivp(f, (0, T_F), z0, method='LSODA', rtol=1e-9, atol=1e-12, t_eval=ts)
    Vt = np.array([V(t, sol.y[:, k]) for k, t in enumerate(sol.t)])
    e = sol.y[0] - XR(sol.t)
    unused = np.setdiff1d(np.arange(d), used)
    q = np.array_split(np.arange(sol.t.size), 4)
    Bu = np.array([P.beta(sol.y[1:, k])[unused].max() for k in range(sol.t.size)])
    Wend = P.w(sol.y[1:, -1])
    return dict(dV=worst, V=[float(Vt[i].max()) for i in q], Vinc=float(max(0.0, np.max(np.diff(Vt)))), erms=[float(np.sqrt(np.mean(e[i] ** 2))) for i in q],
                bunused=[float(Bu[i].max()) for i in q], err_end=float(np.abs(Wend - wstar).max()), V0=float(Vt[0]))


def mixed_C():
    """Part C's design and its mixed start (the even atoms at α = 0, the odd ones on branches with the truth's signs on its atoms): the
    same random draws as `part_C`. Returns a dict with the law, Φ, y, the frozen signs, the amplitudes and the branch atoms' own
    sign-constrained weighted basis pursuit (support S)."""
    L = MF.FreeLaw(**LAW_B)
    rng = np.random.default_rng(7)
    Phi = rng.standard_normal((NC, DC))
    ws = np.zeros(DC); ws[[1, 5, 9]] = [0.8, -0.6, 0.5]; y = Phi @ ws
    sg = rng.choice([-1.0, 1.0], DC); sg[[1, 5, 9]] = np.sign(ws[[1, 5, 9]])
    amp0 = rng.uniform(0.2, 0.6, DC)
    a_init = L.A * np.tanh(amp0 / L.A)
    odd = np.arange(DC) % 2 == 1
    wsc = np.zeros(DC); wsc[odd] = MF.weighted_bp(Phi[:, odd], y, 1 / a_init[odd], signs=sg[odd])
    S = np.flatnonzero(np.abs(wsc) > 1e-9)
    return dict(L=L, Phi=Phi, y=y, sg=sg, amp0=amp0, rho=1 / (L.gamma * a_init), odd=odd, wsc=wsc, S=S)


def mixed_qp(Phi, S, sg, rho, ob, fr):
    """The leading order of the dual off range Φ_S: ζ* = argmin ½Σ_{j∈fr}(Φᵀζ)_j² subject to Φ_Sᵀζ = ρ_S·s_S and s_j(Φᵀζ)_j ≤ −ρ_j for
    j ∈ ob (the off-support branch atoms), by enumerating the active sets (KKT: feasible, multipliers ≥ 0). Returns (ζ*, active atoms,
    the largest KKT violation)."""
    from itertools import combinations
    N = Phi.shape[0]
    M = Phi[:, fr] @ Phi[:, fr].T
    for k in range(ob.size + 1):
        for act in combinations(ob.tolist(), k):
            act = list(act)
            Aeq = np.vstack([Phi[:, S].T] + [sg[j] * Phi[:, j][None, :] for j in act])
            beq = np.r_[rho[S] * sg[S], [-rho[j] for j in act]]
            K = np.block([[M, Aeq.T], [Aeq, np.zeros((Aeq.shape[0], Aeq.shape[0]))]])
            try:
                sol = np.linalg.solve(K, np.r_[np.zeros(N), beq])
            except np.linalg.LinAlgError:
                continue
            z, lam = sol[:N], sol[N + S.size:]
            viol = max([0.0] + [sg[j] * Phi[:, j] @ z + rho[j] for j in ob] + [-v for v in lam])
            if viol <= 1e-9 * (1 + np.abs(rho).max()):
                return z, act, float(viol)
    return None, [], float('inf')


def mixed_resolved(m, eps, z0, exact_free):
    """The Bregman projection of the mixed start at masses ε, its dual z solved on both scales: Newton on Uᵀ(Φw − y) = 0 (U spanning
    range Φ_S) and VᵀΦ_{S^c}w_{S^c}/s = 0 (V its complement, s = Σ_{j∉S}|w_j||φ_j| at the pass's start, so both blocks are O(1)), with
    backtracking on the residual's norm, from z0; a second pass rescales s at the first pass's end. Branch particles by
    `ParticleSet`; free particles by `ParticleSet` (exact_free) or by their linear core w = G(β⁰)²η/m, dw/dη = G(β⁰)²/m (below
    ε ≈ 10⁻¹⁶ their Q = ε − c is not representable; at ε = 10⁻¹² the core's relative error is O(ςGεη²) ≈ 10⁻⁶). Returns a dict."""
    L, Phi, y, sg, amp0, odd, S = m['L'], m['Phi'], m['y'], m['sg'], m['amp0'], m['odd'], m['S']
    oi, fr = np.flatnonzero(odd), np.flatnonzero(~odd)
    off = np.setdiff1d(np.arange(DC), S)
    Uq = np.linalg.qr(Phi[:, S], mode='complete')[0]
    U, V = Uq[:, :S.size], Uq[:, S.size:]
    Pb = MF.ParticleSet(L, sg[oi] * amp0[oi], np.full(oi.size, eps))
    Pf = MF.ParticleSet(L, np.zeros(fr.size), np.full(fr.size, eps)) if exact_free else None
    e0 = np.zeros(DC); e0[oi] = Pb.eta(Pb.x0)

    nphi = np.linalg.norm(Phi[:, off], axis=0)

    def at(z, sc):
        e = e0 + Phi.T @ z
        xb = Pb.x_of_eta(e[oi])
        w = np.zeros(DC); mD = np.zeros(DC)
        w[oi], mD[oi] = Pb.w(xb), Pb.mD(xb)
        if exact_free:
            xf = Pf.x_of_eta(e[fr]); w[fr], mD[fr] = Pf.w(xf), Pf.mD(xf)
        else:
            w[fr], mD[fr] = L.G * eps ** 2 * e[fr] / L.m, L.G * eps ** 2 / L.m
        return np.r_[U.T @ (Phi @ w - y), V.T @ (Phi[:, off] @ w[off]) / sc], w, mD, e
    z = np.array(z0, float)
    it = 0
    for pas in range(2):
        sc = float(np.sum(np.abs(at(z, 1.0)[1][off]) * nphi))
        f, w, mD, e = at(z, sc)
        for k in range(400):
            it += 1
            J = np.vstack([U.T @ Phi @ (mD[:, None] * Phi.T), V.T @ Phi[:, off] @ (mD[off, None] * Phi[:, off].T) / sc])
            dz = -np.linalg.lstsq(J, f, rcond=None)[0]
            n0, lam, acc = np.linalg.norm(f), 1.0, None
            while lam > 1e-14:
                try:
                    with np.errstate(over='ignore', invalid='ignore'):
                        f2, w2, mD2, e2 = at(z + lam * dz, sc)
                        ok = np.all(np.isfinite(f2)) and np.linalg.norm(f2) < (1 - 1e-4 * lam) * n0
                    if ok:
                        acc = (z + lam * dz, f2, w2, mD2, e2)
                        break
                except ValueError:
                    pass
                lam /= 2
            if acc is None:                                                   # no decrease: the rounding floor
                break
            z, f, w, mD, e = acc
    lg = np.log(1 / eps)
    ob = np.array([j for j in off if odd[j]])
    scale = float(np.sum(np.abs(w[off]) * np.linalg.norm(Phi[:, off], axis=0)))
    with np.errstate(divide='ignore'):
        expo = np.log(np.abs(w[ob])) / np.log(eps)
    return dict(eps=eps, exact=bool(exact_free), it=it, res=float(np.linalg.norm(Phi @ w - y)),
                perp=float(np.linalg.norm(V.T @ (Phi[:, off] @ (w[off] / scale)))) if scale > 0 else float('nan'),
                free=float(np.abs(w[fr]).sum()), free_e2=float(np.abs(w[fr]).sum() / eps ** 2),
                free_e2L=float(np.abs(w[fr]).sum() / eps ** 2 / lg), eta_free=float(np.abs(e[fr]).max() / lg),
                offout=float(scale * np.linalg.norm(Phi[:, off] @ (w[off] / scale))) if scale > 0 else 0.0,
                dist=float(np.abs(w - m['wsc']).max()),
                expo={int(j): float(x) for j, x in zip(ob, expo)}, z=z.tolist())


def part_G():
    m = mixed_C()
    L, Phi, y, sg, odd, S, rho = m['L'], m['Phi'], m['y'], m['sg'], m['odd'], m['S'], m['rho']
    fr = np.flatnonzero(~odd); off = np.setdiff1d(np.arange(DC), S)
    ob = np.array([j for j in off if odd[j]])
    from scipy.optimize import linprog
    V = np.linalg.qr(Phi[:, S], mode='complete')[0][:, S.size:]
    A_ = np.array([sg[j] * V.T @ Phi[:, j] for j in ob]).T
    lp = linprog(np.zeros(ob.size), A_eq=np.vstack([A_, np.ones(ob.size)]), b_eq=np.r_[np.zeros(V.shape[1]), 1.0],
                 bounds=[(0, None)] * ob.size, method='highs')
    zs, act, viol = mixed_qp(Phi, S, sg, rho, ob, fr)
    const = float(L.G / L.m * np.abs(Phi[:, fr].T @ zs).sum())
    # C's projection (`projection`, exact=True) at the same masses, from which the continuation starts
    cproj = []
    z = None
    for eps in (1e-4, 1e-8, 1e-12):
        P = MF.ParticleSet(L, np.where(odd, sg * m['amp0'], 0.0), np.full(DC, eps))
        e0 = P.eta(P.x0)
        wb, it, res = MF.projection(P, Phi, y, e0, exact=True)
        cproj.append(dict(eps=eps, res=res, free=float(np.abs(wb[fr]).sum()), dist=float(np.abs(wb - m['wsc']).max())))
        if z is None:
            st = P.state_of_w(wb, np.flatnonzero(wb != 0))
            z = np.linalg.lstsq(Phi.T, P.eta(st) - e0, rcond=None)[0]
    rows, lg_prev = [], None
    for eps in EPS_G:
        lg = np.log(1 / eps)
        z = z if lg_prev is None else z * lg / lg_prev                     # continuation: z ≈ log(1/ε)ζ*
        r = mixed_resolved(m, eps, z, eps >= EPS_G_EXACT)
        r['zdev'] = float(np.linalg.norm(np.array(r['z']) / lg - zs))
        rows.append(r)
        z, lg_prev = np.array(r['z']), lg
    chk = mixed_resolved(m, EPS_G_EXACT, rows[[r['eps'] for r in rows].index(EPS_G_EXACT)]['z'], False)
    core_err = abs(chk['free'] / rows[[r['eps'] for r in rows].index(EPS_G_EXACT)]['free'] - 1)
    return dict(S=S.tolist(), N=NC, ob=ob.tolist(), gordan=bool(lp.status != 0), zstar=zs.tolist(), act=[int(j) for j in act],
                qp_viol=viol, const=const, cproj=cproj, tol_proj=float(1e-12 * (1 + np.linalg.norm(y))), rows=rows, core_err=core_err)


def report(A=None, B=None, C=None, D=None, Ee=None, F=None, G=None):
    fmt = E.fmt
    L_ = ["# Proposition D‴.47: the free law's branches (`exp_mirror_branch.py`)", "",
          "The free mass–amplitude law (Proposition D‴.45) with particles on branches (Q < −c), represented without cancellation near their "
          "branch points (`mirror_free.ParticleSet`).", ""]
    if A is not None:
        L_ += ["## A. Identities", "",
               f"dη/dlog|w| against a/(γa² + Gβa′²) at 40 random states (free and on branches; finite differences): largest relative "
               f"difference {fmt(A['slope'])}.", "", "The temperature at the branch point against 1/(γa₀):", "", "| α₀ | a/(γa² + Gβa′²) | 1/(γa₀) |",
               "|---|---|---|"]
        L_ += [f"| {r['a0']:g} | {r['T']:.10g} | {r['ref']:.10g} |" for r in A['temp']]
        L_ += ["", "B(w̄, w) on a branch as w̄ → 0 (δ̄ the distance to the branch point) against β/(γm):", "",
               "| δ̄ | w̄ | B(w̄, w) | β/(γm) |", "|---|---|---|---|"]
        L_ += [f"| {r['dbar']:g} | {fmt(r['wbar'])} | {r['B']:.12g} | {r['ref']:.12g} |" for r in A['b0']]
        L_ += ["", "The fraction of starts on branches, α(0) uniform on [−0.3, 0.3]: "
               + ", ".join(f"{r['frac']:.4f} at β(0) = {r['b']:g}" for r in A['frac']) + ".", "",
               f"A ± pair (α(0) = ±0.3, u = w₊(0) = {fmt(A['u'])}, a₀ = {A['a0']:.6g}): the atom's coefficient W(Δ) against 2u sinh(γa₀Δ):", "",
               "| Δ | W/u | relative difference |", "|---|---|---|"]
        L_ += [f"| {r['D']:g} | {fmt(r['ratio'])} | {fmt(r['rel'])} |" for r in A['pair']]
        L_ += ["", "The critical curve Q = −c: η·w^{1/3} (η by quadrature, up to its constant) against −(a_max²/(2Gc))(2c/(ma_max²))^{1/3} = "
               f"{A['crit'][0]['ref']:.6g}: " + ", ".join(f"{r['val']:.6g} at α = {r['alpha']:g}" for r in A['crit']) + ".", ""]
    if B is not None:
        L_ += [f"## B. The limit on random Gaussian designs ({NB} × {DB}, the target on four atoms with signs +, −, +, −; masses 10⁻³, amplitudes "
               "0.2–0.4)", "",
               "The flow (LSODA) to t = 3·10⁷ against the Bregman projection of its start onto S (computed independently); |y − ŷ| the "
               "sign-constrained residual, |Z| the coordinates vanishing on S; J − min_K J at the end; the largest mass on Z at the end; "
               "B(w̄, w(t))'s largest relative increase; the extinction rates (fitted where 10⁻²⁵⁰ < β < 10⁻¹²) against γa₀|g(∞)|:", "",
               "| start | seed | branches | |y − ŷ| | |Z| | projection (it, residual) | max |w(T) − w̄| | J − J* | max β on Z | B's increase | rates: largest relative error (n) |",
               "|---|---|---|---|---|---|---|---|---|---|---|"]
        for r in B:
            if r['kind'] == 'duplicated':
                continue
            rr = max((abs(x['fit'] / x['ref'] - 1) for x in r['rates']), default=float('nan'))
            L_.append(f"| {r['kind']} | {r['seed']} | {r['nbr']} | {fmt(r['yres'])} | {len(r['Z'])} | {r['it']}, {fmt(r['res'])} | {fmt(r['dist'][-1])} | "
                      f"{fmt(r['dJ'][-1])} | {fmt(r['bZ'][-1])} | {fmt(r['Bincr'])} | {fmt(rr)} ({len(r['rates'])}) |")
        dup = [r for r in B if r['kind'] == 'duplicated']
        if dup:
            L_ += ["", f"Duplicated atoms (round 36; {ND_B} samples, {DD_B} Gaussian atoms each carried by two particles — the least-squares set S "
                   "is not a point —, y random (the fit inexact) and one atom's frozen sign against its unconstrained least-squares coefficient (Z ≠ ∅ where it "
                   "binds), "
                   "amplitudes ±U(0.4, 0.8), masses 10⁻²): the flow against the Bregman projection, and the projection's distance from the "
                   "equal split of the duplicates (another point of S):", "",
                   "| seed | branches | |y − ŷ| | Z | negative coefficients alive | projection vs the equal split | max |w(T) − w̄| | max β on Z |",
                   "|---|---|---|---|---|---|---|---|"]
            L_ += [f"| {r['seed']} | {r['nbr']} | {fmt(r['yres'])} | {r['Z']} | {r['nneg']} | {fmt(r['split'])} | {fmt(r['dist'][-1])} | "
                   f"{fmt(r['bZ'][-1])} |" for r in dup]
        L_.append("")
    if C is not None:
        L_ += [f"## C. The rich limits on one design ({NC} × {DC}; the target on atoms [1, 5, 9])", "",
               f"Sign-constrained basis pursuit weighted by 1/|a(0)| (the frozen signs: the target's on its atoms, random elsewhere): support "
               f"{C['scbp_supp']} (the target: {'yes' if C['scbp_is_truth'] else 'no'}); basis pursuit: support {C['bp_supp']}.", "",
               "The vertices' certificates (|S| = N, ξ from Φ_Sᵀξ = ρ_S·sign): the smallest relative slacks off the support μ_j = 1 − s_j(Φᵀξ)_j/ρ_j — "
               f"sign-constrained, weights 1/(γ|a(0)|): {', '.join(f'{v:.4f}' for v in C['mu_sc'][:3])}; basis pursuit (equal weights): "
               f"{', '.join(f'{v:.4f}' for v in C['mu_bp'][:3])}.", "",
               "The projection's distance to its limit against the leading order computed from the certificate (the support's mirror coordinates "
               "at the vertex, the off-support coefficients through each particle's mirror map: ≈ e^{O(1)}ε^{μ_j} on branches, "
               "a_maxε^{3/2}cos(πμ_j/2)/(m(2c)^{1/2}sin³(πμ_j/2)) from α(0) = 0 — the closed form's largest relative error off the support in the "
               "last column):", "",
               "| start | ε | all on branches | projection (it, residual) | distance to its limit | the leading order's distance | "
               "|w − w₁|∞ | the closed form |", "|---|---|---|---|---|---|---|---|"]
        L_ += [f"| branches (amplitudes 0.2–0.6) | {r['eps']:g} | {'yes' if r['all_branch'] else 'no'} | {r['it']}, {fmt(r['res'])} | {fmt(r['d'])} | "
               f"{fmt(r['d1'])} | {fmt(r['agree'])} | — |" for r in C['br']]
        L_ += [f"| α(0) = 0 | {r['eps']:g} | no | {r['it']}, {fmt(r['res'])} | {fmt(r['d'])} | {fmt(r['d1'])} | {fmt(r['agree'])} | {fmt(r['cf'])} |"
               for r in C['fr']]
        L_ += ["", "A mixed start (the even atoms at α = 0, the odd ones on branches): the free particles' total |w| and the distance to the "
               "branch atoms' own sign-constrained weighted basis pursuit:", "", "| ε | projection (it, residual) | Σ|w| on free particles | distance |",
               "|---|---|---|---|"]
        L_ += [f"| {r['eps']:g} | {r['it']}, {fmt(r['res'])} | {fmt(r['free'])} | {fmt(r['d'])} |" for r in C['mx']]
        failed = lambda r: '' if r['res'] < C.get('tol_r', np.inf) else ' (the projection failed)'
        if 'mx2' in C:
            L_ += ["", f"The same with the odd atoms' frozen signs against the truth's on its atoms (round 36; their own sign-constrained fit "
                   f"{'exists' if C['odd2_feasible'] else 'does not exist'}): the free particles enter on the power scale; the distance to basis "
                   "pursuit over the free atoms with the branch atoms sign-constrained and free of charge (Σ|w| on the free atoms there: "
                   f"{fmt(C['mix2_free'])}):", "", "| ε | projection (it, residual) | Σ|w| on free particles | distance |", "|---|---|---|---|"]
            L_ += [f"| {r['eps']:g} | {r['it']}, {fmt(r['res'])}{failed(r)} | {fmt(r['free'])} | {fmt(r['d'])} |" for r in C['mx2']]
        if 'nfail' in C:
            L_ += ["", f"Projections whose residual |Φw − y| exceeds {fmt(C['tol_r'])}: {C['nfail']}."]
        L_.append("")
    if D is not None:
        L_ += ["## D. A pair per atom (± 0.4, equal masses; the design of C)", "",
               "The atom's coefficient W = w₊ + w₋ of the projection against basis pursuit weighted by 1/a₀; the atoms where both members "
               f"carry more than 10⁻³ of the largest |W| (basis pursuit's support {D['bp_supp']}; its certificate's smallest relative slacks "
               f"{', '.join(f'{v:.4f}' for v in D['mu_bp'][:3])}); the leading order from the certificate as in C:", "",
               "| u = w(0) | projection (it, residual) | distance | the leading order's distance | |W − W₁|∞ | atoms with both members active |",
               "|---|---|---|---|---|---|"]
        L_ += [f"| {r['u']:g} | {r['it']}, {fmt(r['res'])} | {fmt(r['d'])} | {fmt(r['d1'])} | {fmt(r['agree'])} | {r['both']} |" for r in D['rows']]
        if 'nfail' in D:
            L_ += ["", f"Projections whose residual exceeds {fmt(D['tol_r'])}: {D['nfail']}."]
        L_.append("")
    if Ee is not None:
        L_ += ["## E. Positive kernel dictionaries (Gaussian atoms of one width on a grid of [−1, 1], samples evenly on [−1.1, 1.1], the target "
               "1.0, 0.7, 0.5 on two adjacent atoms and one more)", ""]
        for c in Ee:
            L_ += [f"### {c['d']} atoms, width {c['wid']:g}, {c['N']} samples (the target on {c['supp']})", "",
                   f"The nonnegative certificate's margin μ = {fmt(c['mu'])}; basis pursuit's max |w − w*|: {fmt(c['bp'])} with equal weights, "
                   f"{fmt(c['bpn'])} with the columns' norms.", "",
                   "| start | α(0) | β(0) | on branches | max |w − w*| at the times nearest 10³, 10⁵, 10⁷, 10⁹ | the other atoms' mass at 10⁹ (the "
                   "certificate's bound) | |Φw − y| | support at 10⁹ |", "|---|---|---|---|---|---|---|---|"]
            for r in c['rows']:
                ts = np.array(r['t'])
                pick = [int(np.argmin(np.abs(np.log10(np.maximum(ts, 1e-300)) - k))) for k in (3, 5, 7, 9)]
                if r is c['rows'][0]:
                    L_[-2] = L_[-2].replace('the times nearest 10³, 10⁵, 10⁷, 10⁹', 't = ' + ', '.join(f'{ts[k]:.3g}' for k in pick))
                bd = f" ({fmt(r['bound'][-1])})" if 'bound' in r else ""
                L_.append(f"| {r['lab']} | {r['al']:g} | {r['b']:g} | {'yes' if r['branch'] else 'no'} | "
                          + ", ".join(fmt(r['err'][k]) for k in pick) + f" | {fmt(r['offmass'][-1])}{bd} | {fmt(r['res'])} | {r['supp']} |")
            if 'gen' in c:
                L_ += ["", f"The margin over all {c['n_supp']} supports of three atoms (round 36): the 10, 25, 50, 75, 90% quantiles "
                       f"{', '.join(fmt(v) for v in c['mu_q'])}; {c['mu_zero']:.0%} at most 10⁻⁹ (numerically no certificate); the target above "
                       f"{c['mu_rank']:.0%} of them. The loop from ({STARTS_E[GEN_START_E][1]:g}, {STARTS_E[GEN_START_E][2]:g}) on positive branches "
                       "to t = 10⁹ on the supports at the margin's quantiles (amplitudes 1.0, 0.7, 0.5):", "",
                       "| quantile | support | margin μ | max |w − w*| at t = 10⁵, 10⁷, 10⁹ | the other atoms' mass at 10⁹ | |Φw − y| |",
                       "|---|---|---|---|---|---|"]
                for g_ in c['gen']:
                    tg = np.array(g_['t'])
                    pk = [int(np.argmin(np.abs(np.log10(np.maximum(tg, 1e-300)) - k))) for k in (5, 7, 9)]
                    L_.append(f"| {g_['q']:g} | {g_['supp']} | {fmt(g_['mu'])} | " + ", ".join(fmt(g_['err'][k]) for k in pk)
                              + f" | {fmt(g_['offmass'][-1])} | {fmt(g_['res'])} |")
            L_.append("")
    if F is not None:
        L_ += [f"## F. The real loop on positive branches (k_e = {KE_F:g}, G = γ = 1, twelve atoms of width 0.25, the target on atoms [2, 6, 9], "
               f"α(0) = 0.5, β(0) = 10⁻³, {T_F:g} time units)", "",
               f"V̇ against −k_ee² at 12 random states (w* vanishing on nine atoms, their terms β/(γm)): largest relative difference {fmt(F['dV'])}.", "",
               "By quarters: max V, the rms error, the largest mass on the unused atoms; V's largest increase between samples "
               f"{fmt(F['Vinc'])} (V(0) = {fmt(F['V0'])}); max |w − w*| at the end {fmt(F['err_end'])}:", "",
               "| | quarter 1 | quarter 2 | quarter 3 | quarter 4 |", "|---|---|---|---|---|",
               "| max V | " + " | ".join(fmt(x) for x in F['V']) + " |",
               "| rms e | " + " | ".join(fmt(x) for x in F['erms']) + " |",
               "| max β, unused atoms | " + " | ".join(fmt(x) for x in F['bunused']) + " |", ""]
    if G is not None:
        act = G['act']
        L_ += ["## G. The mixed start of C resolved on both scales (round 144)", "",
               f"The branch atoms' own sign-constrained weighted basis pursuit (C's mixed start) has the support {G['S']}: {len(G['S'])} atoms "
               f"against N = {G['N']}, so the dual directions off range Φ_S couple only to the free particles and to the off-support branch "
               f"atoms {G['ob']}. C's `projection` stops once |Φw − y| < {fmt(G['tol_proj'])}; the off-support coefficients' contribution "
               "to Φw lies below that from ε = 10⁻⁸ on, so C's free totals there are not resolved:", "",
               "| ε | C's projection: |Φw − y| | C's Σ|w| on free particles | resolved Σ|w| on free particles | resolved |Φ_{S^c}w_{S^c}| |",
               "|---|---|---|---|---|"]
        byeps = {r['eps']: r for r in G['rows']}
        L_ += [f"| {c['eps']:g} | {fmt(c['res'])} | {fmt(c['free'])} | {fmt(byeps[c['eps']]['free'])} | {fmt(byeps[c['eps']]['offout'])} |"
               for c in G['cproj']]
        L_ += ["", "Resolved: Newton on Uᵀ(Φw − y) = 0 (U spanning range Φ_S) and VᵀΦ_{S^c}w_{S^c} = 0 (V its complement), scaled by "
               "Σ_{j∉S}|w_j||φ_j|, continued in ε "
               f"from z ∝ log(1/ε); free particles exact for ε ≥ {EPS_G_EXACT:g}, by their linear core below (at {EPS_G_EXACT:g} the core changes "
               f"their total by {fmt(G['core_err'])}, relative). The off-support branch atoms' signed columns projected off range Φ_S have "
               f"{'no' if G['gordan'] else 'a'} vanishing nonnegative combination"
               + (" (a direction shrinks them all: they balance the free particles)" if G['gordan'] else "") + ". The leading order "
               "z = log(1/ε)ζ*, ζ* = argmin Σ_free(Φᵀζ)_j² subject to Φ_Sᵀζ = ρ_S·s_S and s_j(Φᵀζ)_j ≤ −ρ_j off S (ρ = 1/(γ|a(0)|)): "
               f"active on {act} (KKT violation {fmt(G['qp_viol'])}), the free particles' total (G/m)Σ|Φᵀζ*| = {G['const']:.4g} times "
               "ε²log(1/ε), the active atoms' exponents 2.", "",
               "| ε | |Φw − y| | complement residual (relative) | Σ|w| on free particles | Σ|w|/ε² | Σ|w|/(ε²log(1/ε)) | "
               "max |η| on free particles / log(1/ε) | exponents log|w_j|/log ε on the active atoms | |z/log(1/ε) − ζ*| | distance to the branch "
               "atoms' basis pursuit |", "|---|---|---|---|---|---|---|---|---|---|"]
        for r in G['rows']:
            ex = [r['expo'][j] for j in act]
            L_.append(f"| {r['eps']:g} | {fmt(r['res'])} | {fmt(r['perp'])} | {fmt(r['free'])} | {r['free_e2']:.4g} | {r['free_e2L']:.4g} | "
                      f"{r['eta_free']:.4g} | {min(ex):.3f}–{max(ex):.3f} | {r['zdev']:.3g} | {fmt(r['dist'])} |")
        L_.append("")
    L_.append("Numbers are produced by `experiments/exp_mirror_branch.py`; nothing is edited by hand.")
    return '\n'.join(L_) + '\n'


def main():
    raw = os.path.join(OUT, 'mirror_branch_raw.npz')
    md = os.path.join(OUT, 'mirror_branch.md')
    old = E.load(raw)
    get = lambda k: old[k] if k in old else None
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(*[get(k) for k in 'ABCDEFG']))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else set('ABCDEFG')
    res = dict(old)
    t0 = time.time()
    for k, fn in (('A', part_A), ('C', part_C), ('D', part_D), ('F', part_F), ('B', part_B), ('E', part_E), ('G', part_G)):
        if k in parts:
            res[k] = fn(); E.save(raw, **res); print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
    open(md, 'w').write(report(*[res.get(k) for k in 'ABCDEFG']))


if __name__ == '__main__':
    main()
