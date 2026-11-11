"""
The planar (n = 2) certificate of Theorem T1: limit recursion, `cert` with Euclidean geometry, the integrated planar network, the ingredient checks.

`cert2d` is the source of `exp_theory_t1.cert` with its geometry-dependent expressions replaced (distances become Euclidean norms, the 1-D interval
[c - dc, c + dc] becomes the disc of radius dc, G2_sup includes the tangential eigenvalue Gamma/2 of the Hessian of Gamma); every other line, in particular all
constants of the box (n-independent: radial offsets, operator norms of the cross blocks of the Hessian), is the same.  See THEORY_T1.md section 12.2.
"""
import inspect
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root

HERE = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..'))
import exp_theory_t1 as T      # noqa: E402

Gam = T.Gam


def _nrm(x):
    return np.linalg.norm(np.asarray(x, float), axis=-1)


def gam2(d):
    return np.exp(-np.asarray(d, float) ** 2 / 4.0)


def G2_sup2(x0):
    """sup over y >= x0 of the operator norm of the Hessian of Gamma(|x|): max(|Gamma''(y)|, Gamma(y)/2) (radial and tangential eigenvalues)."""
    x0 = np.maximum(np.asarray(x0, float), 0.0)
    return np.maximum(T.G2_sup(x0), np.exp(-x0 ** 2 / 4) / 2)


# ----------------------------------------------------------------------------------------------------------------------------------
# limit object in the plane
# ----------------------------------------------------------------------------------------------------------------------------------
def crit_block2(cs, As, feats):
    k = len(feats)

    def Fg(z):
        w, c = z[:k], z[k:].reshape(k, 2)
        D = c[:, None, :] - cs[None, :, :]
        E = c[:, None, :] - c[None, :, :]
        g = gam2(_nrm(D)) @ As - gam2(_nrm(E)) @ w
        gradT = -(D / 2) * gam2(_nrm(D))[..., None]
        gradU = -(E / 2) * gam2(_nrm(E))[..., None]
        for i in range(k):
            gradU[i, i] = 0.0
        h = np.einsum('ikd,k->id', gradT, As) - np.einsum('ijd,j->id', gradU, w)
        return np.r_[g, (w[:, None] * h).ravel()]
    s = root(Fg, np.r_[As[feats], cs[feats].ravel()], method='hybr', tol=1e-14)
    return s.x[:k], s.x[k:].reshape(k, 2)


def limit_path2(cs, As, c0):
    K, J = len(cs), len(c0)
    u = np.zeros(J)
    tau = 0.0
    active, assigned = [], []
    w, c = np.zeros(0), np.zeros((0, 2))
    ev, unfit = [], list(range(K))
    tau_br, U_br, P = [0.0], [u.copy()], []
    for k in range(K):
        p = gam2(_nrm(c0[:, None, :] - cs[None, :, :])) @ As - (gam2(_nrm(c0[:, None, :] - c[None, :, :])) @ w if len(w) else 0.0)
        p = np.where(np.isin(np.arange(J), active), 0.0, p)
        dt = np.full(J, np.inf)
        for j in range(J):
            if j in active or abs(p[j]) < 1e-14:
                continue
            dt[j] = (1 - np.sign(p[j]) * u[j]) / abs(p[j])
        order = np.argsort(dt)
        j1 = int(order[0])
        gap = (dt[order[1]] - dt[order[0]]) / dt[order[0]] if np.isfinite(dt[order[1]]) else np.inf
        P.append(p.copy())
        u = u + dt[j1] * p
        tau += dt[j1]
        f = min(unfit, key=lambda m: _nrm(c0[j1] - cs[m]))
        unfit.remove(f)
        active.append(j1)
        assigned.append(f)
        w, c = crit_block2(cs, As, assigned)
        u[j1] = np.sign(p[j1])
        tau_br.append(tau)
        U_br.append(u.copy())
        ev.append(dict(k=k + 1, kernel=j1, feature=f, tau=tau, gap=float(gap), d0=float(_nrm(c0[j1] - cs[f])), p_e=float(abs(p[j1])),
                       sign_ok=bool(np.sign(p[j1]) == np.sign(As[f])), w=w.copy(), c=c.copy()))
    pend = gam2(_nrm(c0[:, None, :] - cs[None, :, :])) @ As - gam2(_nrm(c0[:, None, :] - c[None, :, :])) @ w
    P.append(np.where(np.isin(np.arange(J), active), 0.0, pend))
    return dict(events=ev, tau_br=np.array(tau_br), U_br=np.array(U_br), P=np.array(P), active=active, K=K, J=J)


def prepare2(cs, As, c0):
    cs, As, c0 = np.asarray(cs, float), np.asarray(As, float), np.asarray(c0, float)
    K = len(cs)
    a = np.abs(As)
    amin = float(a.min())
    D = _nrm(cs[:, None, :] - cs[None, :, :]) + np.eye(K) * 1e9
    Delta = float(D.min()) if K > 1 else np.inf
    L = limit_path2(cs, As, c0)
    ent = set(L['active'])
    non = [j for j in range(L['J']) if j not in ent]
    mu = float(min(1 - np.abs(L['U_br'][:, non]).max(axis=1))) if non else np.inf
    return dict(cs=cs, As=As, c0=c0, K=K, a=a, amin=amin, R=float(a.max() / amin), Delta=Delta, L=L, mu_min=mu,
                gap_min=float(min(e['gap'] for e in L['events'])), d0=[e['d0'] for e in L['events']], sign_ok=all(e['sign_ok'] for e in L['events']))


def make_target2(K, sep, R, seed, rad=3.5):
    """K features in the plane with pairwise distances >= sep (first at the origin, the others by rejection in an annulus), amplitudes (1, ..., R)*0.7 in random
    order, random signs; the dictionary is a square lattice of spacing 0.1/sigma at a random phase, restricted to the kernels within rad of a feature."""
    rng = np.random.default_rng(seed)
    pts = [np.zeros(2)]
    while len(pts) < K:
        cand = pts[rng.integers(len(pts))] + rng.uniform(sep, 1.3 * sep) * np.array([np.cos(a := rng.uniform(0, 2 * np.pi)), np.sin(a)])
        if all(np.linalg.norm(cand - p) >= sep for p in pts):
            pts.append(cand)
    cs = np.array(pts)
    cs = cs - cs.mean(0)
    a = np.array([1.0]) if K == 1 else (np.array([1.0, R]) if K == 2 else R ** (np.arange(K) / (K - 1.0)))
    a = a[rng.permutation(K)] * 0.7
    As = a * rng.choice([-1.0, 1.0], K)
    h = T.H_LAT
    lo, hi = cs.min(0) - rad - h, cs.max(0) + rad + h
    ph = rng.uniform(0, h, 2)
    xs, ys = np.arange(lo[0], hi[0], h) + ph[0], np.arange(lo[1], hi[1], h) + ph[1]
    G = np.array([[x, y] for x in xs for y in ys])
    keep = (_nrm(G[:, None, :] - cs[None, :, :]) <= rad).any(1)
    return cs, As, G[keep]


# ----------------------------------------------------------------------------------------------------------------------------------
# Euclidean geometry for the certificate
# ----------------------------------------------------------------------------------------------------------------------------------
def rgs2(prep, stage, center, radius, n=21):
    """sup |grad ghat| and oscillation of ghat^(stage)(c) over the disc |c - center| <= radius (grid + Lipschitz margin)."""
    cs, As, L = prep['cs'], prep['As'], prep['L']
    if stage == 0:
        w, cc = np.zeros(0), np.zeros((0, 2))
    else:
        ev = L['events'][stage - 1]
        w, cc = ev['w'], ev['c']
    g1 = np.linspace(-radius, radius, n)
    G = np.array([[x, y] for x in g1 for y in g1])
    G = G[_nrm(G) <= radius * (1 + 1e-12)] + np.asarray(center, float)
    D = G[:, None, :] - cs[None, :, :]
    gT = gam2(_nrm(D))
    g = gT @ As
    grad = np.einsum('nkd,k->nd', -(D / 2) * gT[..., None], As)
    if len(w):
        E = G[:, None, :] - cc[None, :, :]
        gE = gam2(_nrm(E))
        g = g - gE @ w
        grad = grad - np.einsum('nkd,k->nd', -(E / 2) * gE[..., None], w)
    step = 2 * radius / (n - 1) * np.sqrt(2) / 2 if radius > 0 else 0.0
    sec = 0.5 * (np.abs(As).sum() + np.abs(w).sum())
    gn = float(_nrm(grad).max())
    return gn + step * sec, float(g.max() - g.min() + step * gn)


_RG2 = {}


def rgs2_cached(prep, stage, center, radius):
    k = (id(prep), stage, tuple(np.round(np.asarray(center, float), 12)), round(float(radius), 12))
    if k not in _RG2:
        _RG2[k] = rgs2(prep, stage, center, radius)
    return _RG2[k]


def _build_cert2():
    src = inspect.getsource(T.cert)
    rep = [("np.abs(c0[:, None] - cs[featarr][None, :])", "_nrm(c0[:, None, :] - cs[featarr][None, :, :])"),
           ("abs(cs[feat[k]] - cs[feat[i]])", "_nrm(cs[feat[k]] - cs[feat[i]])"),
           ("resid_grad_sup(prep, k, cs[f_] - u0, cs[f_] + u0)", "resid_grad_sup(prep, k, cs[f_], u0)"),
           ("resid_grad_sup(prep, k, c0[j] - dc, c0[j] + dc)", "resid_grad_sup(prep, k, c0[j], dc)"),
           ("np.abs(c0[j] - cs[featarr])", "_nrm(c0[j][None, :] - cs[featarr])")]
    for a, b in rep:
        assert a in src, a
        src = src.replace(a, b)
    ns = dict(T.__dict__)
    ns.update(_nrm=_nrm, resid_grad_sup=rgs2_cached, G2_sup=G2_sup2)
    exec(compile(src, 'cert2d', 'exec'), ns)
    return ns['cert']


cert2d = _build_cert2()


def best_cert2(prep, eps, pars=None):
    pars = T.PARS if pars is None else pars
    o = None
    for i, par in enumerate(pars):
        o = cert2d(prep, eps, par)
        if o['ok']:
            return o, i
    return o, -1


# ----------------------------------------------------------------------------------------------------------------------------------
# the integrated planar network
# ----------------------------------------------------------------------------------------------------------------------------------
def simulate2(prep, eps, Tfin_tau=None, rtol=1e-9, nlocal=1500, nglob=3000):
    cs, As, c0 = prep['cs'], prep['As'], prep['c0']
    J = len(c0)
    L = prep['L']
    ev = L['events']
    U = T.Unit(eps)
    lam, muc = U.lam, U.mu_c
    tau_end = (1.05 * ev[-1]['tau'] + 1.0) if Tfin_tau is None else Tfin_tau
    Tfin = lam * tau_end
    tk = lam * L['tau_br'][1:]

    def f(t, z):
        x, c = z[:J], z[J:].reshape(J, 2)
        w = T._wfun(x, eps, U)
        D = c[:, None, :] - cs[None, :, :]
        E = c[:, None, :] - c[None, :, :]
        gD, gE = gam2(_nrm(D)), gam2(_nrm(E))
        g = gD @ As - gE @ w
        gradE = -(E / 2) * gE[..., None]
        idx = np.arange(J)
        gradE[idx, idx] = 0.0
        h = np.einsum('jkd,k->jd', -(D / 2) * gD[..., None], As) - np.einsum('jid,i->jd', gradE, w)
        return np.r_[(U.G / U.amax) * (eps + (2 * U.vs - eps) * np.tanh(x) ** 2) * g, (muc * w[:, None] * h).ravel()]
    loc = [np.linspace(max(0, t - 600), min(Tfin, t + 900), nlocal) for t in tk]
    ts = np.unique(np.concatenate([np.linspace(0, Tfin, nglob)] + loc))
    sol = solve_ivp(f, (0, Tfin), np.r_[np.zeros(J), c0.ravel()], method='LSODA', rtol=rtol, atol=np.r_[np.full(J, 1e-14), np.full(2 * J, 1e-12)], t_eval=ts)
    return sol, U


def follows2(prep, eps, sol, U):
    As, cs, c0 = prep['As'], prep['cs'], prep['c0']
    J = len(c0)
    ev = prep['L']['events']
    ent = sorted(e['kernel'] for e in ev)
    x, c = sol.y[:J, -1], sol.y[J:, -1].reshape(J, 2)
    w = T._wfun(x, eps, U)
    thr = min(T.THETA_REL * np.abs(As).max(), 0.1 * np.abs(As).min())
    formed = sorted(np.flatnonzero(np.abs(w) > thr).tolist())
    ok = (formed == ent) and sol.status == 0
    ex = bool(ok and all(abs(abs(w[e['kernel']]) - abs(As[e['feature']])) < 1e-3 * abs(As[e['feature']]) and _nrm(c[e['kernel']] - cs[e['feature']]) < 0.05 for e in ev))
    return bool(ok), ex


def diagnostics2(prep, eps, par):
    """The main ingredient inequalities on the planar trajectory: returns [cert_ok, D, E/D max, shift_ok, T1_ok, T1_ratio, Ee_ok, Ee_ratio, mono, d_ok, w_ok, omega ratio, formed_ok, final coefficient error]."""
    o = cert2d(prep, eps, par)
    if not o['ok']:
        return None
    cs, As, c0 = prep['cs'], prep['As'], prep['c0']
    J, K, L = len(c0), prep['K'], prep['L']
    ev = L['events']
    amin = prep['amin']
    sol, U = simulate2(prep, eps)
    lam = U.lam
    t = sol.t
    x, c = sol.y[:J], sol.y[J:].reshape(J, 2, -1)
    w = T._wfun(x, eps, U)
    eta = np.where(x >= 0, lam - U.s_of_x(np.abs(x)), -(lam - U.s_of_x(np.abs(x))))
    ent = [e['kernel'] for e in ev]
    non = [j for j in range(J) if j not in ent]
    tk = o['tk']
    tb = np.r_[0.0, tk, t[-1]]

    def ubar(j):
        u = np.zeros_like(t)
        for i in range(K + 1):
            m = (t >= tb[i]) & (t <= tb[i + 1])
            u[m] = L['U_br'][i][j] * lam + L['P'][i][j] * (t[m] - tb[i])
        return u
    ratio = max(np.abs(eta[j] - ubar(j)).max() / o['Dv'][j] for j in non)
    r, vbar, Bx = par['r'], o['vbar'], o['Bx']
    sh_ok, T1_ok, T1r, Ee_ok, Eer, mono, d_ok, w_ok = True, True, 0.0, True, 0.0, True, True, True
    for k, e in enumerate(ev):
        j, f = e['kernel'], e['feature']
        wj = np.abs(w[j])
        ta = t[np.argmax(wj >= o['w1'])]
        t1 = t[np.argmax(wj >= o['vsharp'])]
        d = np.linalg.norm(c[j] - cs[f][:, None], axis=0)
        inD = (d <= Bx['dhat'] + 1e-12) & (wj >= Bx['vD_lo'] * amin) & (wj <= Bx['vD_hi'] * amin)
        bad = np.flatnonzero(~inD & (t >= t1))
        trel = t1 if not len(bad) else (t[bad[-1] + 1] if bad[-1] + 1 < len(t) else np.inf)
        shift = ta - tk[k]
        sh_ok &= bool(-o['Sm'][k] - 1e-6 <= shift <= o['Sp'][k] + 1e-6)
        T1_ok &= bool(t1 - ta <= o['T1'] + 1e-6)
        T1r = max(T1r, (t1 - ta) / o['T1'])
        Ee = abs(-par['delta1'] - e['p_e'] * (ta - tk[k]))
        Ee_ok &= bool(Ee <= o['De'][k] + 1e-9)
        Eer = max(Eer, Ee / o['De'][k])
        mono &= bool(np.all(np.diff(wj[t <= trel]) >= -1e-13 * vbar))
        d_ok &= bool(d.max() <= r + 1e-12)
        w_ok &= bool(wj.max() <= vbar + 1e-12)
    om = float(np.abs(w[non]).sum(0).max() / amin)
    thr = min(T.THETA_REL * np.abs(As).max(), 0.1 * np.abs(As).min())
    formed = sorted(np.flatnonzero(np.abs(w[:, -1]) > thr).tolist())
    ferr = max(abs(abs(w[ent[k], -1]) - abs(As[ev[k]['feature']])) for k in range(K))
    return [1.0, o['D'], ratio, float(sh_ok), float(T1_ok), T1r, float(Ee_ok), Eer, float(mono), float(d_ok), float(w_ok), om / max(o['omega_peak'], 1e-300), float(formed == sorted(ent)), ferr]
