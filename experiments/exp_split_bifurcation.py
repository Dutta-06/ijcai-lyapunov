"""
Proposition D‴.35: the split of a class under the cohesion leak is an S_n-symmetric transcritical bifurcation
==========================================================================================================

Proposition D‴.32's resolution instance (the fourth bump, ε = 0.3, the three-atom approximant ϑ^b(q) lifted with the counts (10, 10, 10),
m = 30, the compact kernel at R = 1): the middle class splits first, at σ = λ⁺/w_k. Along the splitting direction v (the top whitened
eigenvector of its residual curvature) the separations ξ_i (Σξ = 0) obey to second order

    dξ_i/ds = μξ_i + (c₃/2)(ξ_i² − S₂/n),   μ = λ⁺ − σw_k,   c₃ = ∂³_t⟨r Φ̃(·; ϑ_k + tG^{1/2}v)⟩,

the standard representation of S_n with its quadratic equivariant; the cohesion, even, enters at third order.

A. The reduced coefficients at the lift: per class λ⁺, v, c₃ (central differences at two steps) and the thresholds λ⁺/w_k.
B. The branches: at μ = 2·10⁻⁶ and 5·10⁻⁷ the S_q × S_{n−q} split states (q = 1, …, 4) by Newton on the full particle drift from the
   predicted t* = −2μ/(c₃(n − 2q)), their positions against the prediction and the Jacobian's eigenvalues on the class's separations
   along v against −μn/(n − 2q) (q − 1), +μn/(n − 2q) (n − q − 1) and −μ (radial); the number of unstable eigenvalues of the full Jacobian.
C. Who leaves, when and where: from the lift with a jitter 10⁻³ (seeds 0–5) at μ = 1.39·10⁻⁴ (σ₀/16, Part G's first split) and at
   μ = 5·10⁻⁵: the ranking of sign(c₃)ξ_i(0); the order in which particles leave the class (a particle leaves when it is more than 0.1
   from the median of the class's ξ); the first departure's time against the reduced quadratic equation integrated from the measured
   ξ(0), without the quadratic term (linear) and with c₃'s sign flipped; its side against sign(c₃); the class's pieces beyond the kernel's
   range at the end (single linkage at R, round 25); the margin of the class's main piece at the end.
D. The cascade and two changes of the law (Part G's schedule: σ halved every 10⁵ slow units from σ₀ to σ₀/64, the jitter renewed at each
   stage, seeds 0–2): the cohesion leak with the replicator weight law (Proposition D‴.34, γ = 0.1), and the mass-free leak σ(m_i − ϑ_i)
   (the normalised mean shift: separation stiffness σ instead of σw_k), against Part G's counts.
E. Part D's schedule with the replicator at γ = 1 (seeds 0–2).
F. The mass-free leak continued to σ₀/256 (two more stages, seeds 0–2; round 25), replacing Part D's mass-free rows.

Output: results/chaos/split_bifurcation.md and split_bifurcation_raw.npz (SB_WORKERS; `--from-raw` regenerates the md; `--only A,B,…`).
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
from src.core import cohesion_leak as CO
from src.core import mass_transfer as MT
import exp_cohesion_leak as E

OUT = os.environ.get('SB_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
N_WORKERS = int(os.environ.get('SB_WORKERS', 3))
COUNTS = np.array([10, 10, 10])
MID = 1                                                                      # the middle class splits first
MUS_B = (2e-6, 5e-7)
SEEDS_C = tuple(range(6))
MUS_C = (5e-5, None)                                                         # None: σ₀/16
S_C, LEAVE = 100000.0, 0.1
GAMMA_D, SEEDS_D = 0.1, (0, 1, 2)
STAGES = (1, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625)
STAGES_MF_EXTRA = (0.0078125, 0.00390625)                                    # round 25 (D3x-1): the mass-free leak's thresholds are w_k lower


def setup():
    net, base, bq = E.approximant()
    w = COUNTS / COUNTS.sum()
    at = E.gauge(net, bq, w)
    Th0 = np.repeat(at, COUNTS, axis=0)
    return net, base, at, w, Th0


def coefficients(base, at, Th0, k, hs=(2e-3, 1e-3)):
    """λ⁺, the whitened unit eigenvector v (sign fixed: first nonzero centre component positive) and c₃ at two steps for class k."""
    sqG = np.sqrt(E.G)
    r0 = base.residual(Th0)
    Rk = base.residual_curvature(r0, at[k:k + 1])[0]
    ev, V = np.linalg.eigh(sqG[:, None] * Rk * sqG[None, :])
    v = V[:, -1] * np.sign(V[1, -1])
    g = lambda t: float(base.grad_weighted(r0, (at[k] + t * sqG * v)[None, :])[0] @ (sqG * v))
    c3 = [(g(h) - 2 * g(0.0) + g(-h)) / h ** 2 for h in hs]
    return float(ev[-1]), v, c3, ev


def part_A():
    net, base, at, w, Th0 = setup()
    rows = []
    for k in range(3):
        lam, v, c3, ev = coefficients(base, at, Th0, k)
        rows.append(dict(k=k, w=float(w[k]), lam=lam, v=v, c3=c3, ev=ev, sigma_c=lam / w[k]))
    return dict(rows=rows, sigma0=E.SIGMA0)


def drift_fn(loop):
    return lambda x: loop.drift(x.reshape(-1, 3)).ravel()


def jac(F, x, h=1e-6):
    d = x.size
    J = np.zeros((d, d))
    for j in range(d):
        e = np.zeros(d); e[j] = h
        J[:, j] = (F(x + e) - F(x - e)) / (2 * h)
    return J


def b_job(mu_target):
    t0 = time.time()
    net, base, at, w, Th0 = setup()
    n = int(COUNTS[MID])
    lam, v, c3s, _ = coefficients(base, at, Th0, MID)
    c3 = c3s[0]
    sqG = np.sqrt(E.G)
    sigma = (lam - mu_target) / w[MID]
    loop = CO.CohesionLoop(base, E.R_WARM, 'compact', sigma0=sigma)
    mu = lam - sigma * w[MID]
    F = drift_fn(loop)
    idx = np.arange(COUNTS[:MID].sum(), COUNTS[:MID + 1].sum())
    Q = np.linalg.qr(np.vstack([np.ones(n), np.eye(n)[:-1]]).T)[0][:, 1:]    # an orthonormal basis of {Σζ = 0}
    B = np.zeros((Th0.size, n - 1)); Wl = np.zeros((n - 1, Th0.size))
    for j in range(n - 1):
        U = np.zeros_like(Th0); U[idx] = Q[:, j][:, None] * (sqG * v)[None, :]; B[:, j] = U.ravel()
        U = np.zeros_like(Th0); U[idx] = Q[:, j][:, None] * (v / sqG)[None, :]; Wl[j] = U.ravel()
    out = []
    for q in (1, 2, 3, 4):
        t = -2 * mu / (c3 * (n - 2 * q))
        xi = np.r_[np.full(q, (n - q) * t), np.full(n - q, -q * t)]
        Th = Th0.copy(); Th[idx] += xi[:, None] * (sqG * v)[None, :]
        x = Th.ravel()
        for _ in range(40):
            f = F(x)
            if np.abs(f).max() < 1e-17:
                break
            x = x - np.linalg.lstsq(jac(F, x), f, rcond=None)[0]
        J = jac(F, x)
        T = x.reshape(-1, 3)
        xn = ((T[idx] - T[idx].mean(axis=0)) / sqG[None, :]) @ v
        evr = np.sort(np.linalg.eigvals(Wl @ J @ B).real)
        pred = np.sort(np.r_[np.full(q - 1, -mu * n / (n - 2 * q)), np.full(n - q - 1, mu * n / (n - 2 * q)), [-mu]])
        out.append(dict(q=q, t=t, drift=float(np.abs(F(x)).max()), groups=(float(xn[:q].mean()), float(xn[q:].mean())),
                        pred_groups=((n - q) * t, -q * t), ev=evr, pred=pred, n_unstable=int(np.sum(np.linalg.eigvals(J).real > 1e-12))))
    print(f"[{time.time() - t0:.0f}s] B mu {mu_target:g}", flush=True)
    return dict(mu=mu, sigma=sigma, c3=c3, rows=out)


def reduced_crossing(xi0, mu, c3, level, s_max, dt=5.0):
    """Integrate dξ_i/ds = μξ_i + (c₃/2)(ξ_i² − S₂/n) (RK4) from ξ(0); the first time a particle is `level` from the median of the ξ, and which."""
    xi = np.array(xi0, float); n = xi.size
    f = lambda x: mu * x + 0.5 * c3 * (x ** 2 - np.sum(x ** 2) / n)
    s = 0.0
    while s < s_max:
        dev = np.abs(xi - np.median(xi))
        if dev.max() > level:
            return s, int(np.argmax(dev))
        k1 = f(xi); k2 = f(xi + 0.5 * dt * k1); k3 = f(xi + 0.5 * dt * k2); k4 = f(xi + dt * k3)
        xi = xi + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
        if not np.all(np.isfinite(xi)):
            break
        s += dt
    return np.nan, -1


def c_job(args):
    mu_target, seed = args
    t0 = time.time()
    net, base, at, w, Th0 = setup()
    lam, v, c3s, _ = coefficients(base, at, Th0, MID)
    c3 = c3s[0]
    sqG = np.sqrt(E.G)
    sigma = E.SIGMA0 / 16 if mu_target is None else (lam - mu_target) / w[MID]
    mu = lam - sigma * w[MID]
    rng = np.random.default_rng(seed)
    Th = Th0 + 1e-3 * rng.standard_normal(Th0.shape) * sqG[None, :]
    idx = np.arange(COUNTS[:MID].sum(), COUNTS[:MID + 1].sum())
    xi_of = lambda T: ((T[idx] - T[idx].mean(axis=0)) / sqG[None, :]) @ v
    xi0 = xi_of(Th)
    left = {}

    def cb(s_, T):
        x = xi_of(T)
        dev = np.abs(x - np.median(x))
        for i in np.where(dev > LEAVE)[0]:
            if int(i) not in left:
                left[int(i)] = (s_, float(np.sign(x[i] - np.median(x))))
    T = CO.CohesionLoop(base, E.R_WARM, 'compact', sigma0=sigma).run(Th, S_C, E.DT, rec_every=40, callback=cb)
    s_pred, i_pred = reduced_crossing(xi0, mu, c3, LEAVE, 4 * S_C)
    # the cascade: the margin σw − λ_max(R^G) of the largest remaining piece of the middle class at the end, with the current residual
    lab = CO.clusters(T, E.G, E.CLUSTER_THR)
    main = np.bincount(lab[idx]).argmax()
    keep = idx[lab[idx] == main]
    rem_margin = float(sigma * len(keep) / T.shape[0] - base.kappa_R(base.residual(T), T[keep].mean(axis=0, keepdims=True))[0])
    # round 25 (D3x-3): the pieces of the class at the kernel's range (single linkage at R), with their distances from the main piece
    labR = CO.clusters(T[idx], E.G, E.R_WARM)
    mainR = np.bincount(labR).argmax()
    cm = T[idx][labR == mainR].mean(axis=0)
    pieces = [(np.where(labR == k)[0].tolist(), float(np.sqrt(np.sum((T[idx][labR == k].mean(axis=0) - cm) ** 2 / E.G))))
              for k in range(labR.max() + 1) if k != mainR]
    order = sorted(left.items(), key=lambda kv: kv[1][0])
    print(f"[{time.time() - t0:.0f}s] C mu {mu_target} seed {seed}", flush=True)
    return dict(mu_target=mu_target, mu=mu, seed=seed, c3=c3, xi0=xi0, rank=np.argsort(-np.sign(c3) * xi0), left=order,
                s_pred=s_pred, i_pred=i_pred, K_end=int(lab.max() + 1), fit_end=base.fit_error(T), n_main=int(len(keep)), rem_margin=rem_margin, pieces=pieces)


class MassFree(CO.CohesionLoop):
    """The normalised mean shift σ(m_i − ϑ_i), m_i the kernel-weighted mean: separation stiffness σ at a clustered configuration."""

    def drift(self, Th):
        D, d2 = CO.pair_d2(Th, self.G)
        P = CO.psi(d2 / self.R ** 2, self.kind)
        return self.learning(Th) - self.sigma0 * np.einsum('ij,ijk->ik', P, D) / P.sum(axis=1)[:, None]


def d_job(args):
    law, seed = args
    t0 = time.time()
    net, base, at, w, Th = setup()
    rng = np.random.default_rng(seed)
    beta = np.ones(Th.shape[0])
    stages = []
    for f in (STAGES + STAGES_MF_EXTRA if law == 'mass-free' else STAGES):
        Th = Th + 1e-3 * rng.standard_normal(Th.shape) * np.sqrt(E.G)[None, :]
        if law.startswith('replicator'):
            gam = GAMMA_D if law == 'replicator' else float(law.split(' ')[1])
            loop = MT.ReplicatorCohesionLoop(base, E.R_WARM, 'compact', sigma0=E.SIGMA0 * f, gamma=gam)
            Th, beta = loop.run(Th, beta, E.S_SPLIT, E.DT, rec_every=10 ** 9)
            fit = loop.fit_error(Th, beta)
        else:
            Th = MassFree(base, E.R_WARM, 'compact', sigma0=E.SIGMA0 * f).run(Th, E.S_SPLIT, E.DT)
            fit = base.fit_error(Th)
        lab, W, A = MT.class_masses(Th, beta, E.G, E.CLUSTER_THR)
        a, c, s = net.physical(A)
        o = np.argsort(c[:, 0])
        stages.append(dict(f=f, K=int(len(W)), fit=fit, masses=W[o], n=np.array([np.sum(lab == k) for k in o]), c=c[o, 0], a=a[o, 0]))
    print(f"[{time.time() - t0:.0f}s] D {law} seed {seed}", flush=True)
    return dict(law=law, seed=seed, stages=stages)


# ---------------------------------------------------------------------------- report
def report(A=None, B=None, C=None, D=None, G_ref=None):
    fmt = E.fmt
    L_ = ["# Proposition D‴.35: the split of a class is an S_n-symmetric transcritical bifurcation (`exp_split_bifurcation.py`)", "",
          "Proposition D‴.32's resolution instance: the fourth bump (ε = 0.3), the lift of ϑ^b(q) with the counts (10, 10, 10), m = 30, the compact "
          "kernel at R = 1; ξ_i the whitened separations of the middle class along the top eigenvector v of its residual curvature.", ""]
    if A is not None:
        L_ += ["## A. The reduced coefficients at the lift", "",
               "| class | w_k | λ⁺ = λ_max(R^G_k) | the other eigenvalues | v (whitened: amplitude, centre, width) | c₃ (steps 2·10⁻³, 10⁻³) | threshold λ⁺/w_k (/σ₀) |",
               "|---|---|---|---|---|---|---|"]
        for r in A['rows']:
            L_.append(f"| {r['k'] + 1} | {r['w']:.3f} | {fmt(r['lam'])} | {fmt(r['ev'][0])}, {fmt(r['ev'][1])} | ({', '.join(f'{x:+.4f}' for x in r['v'])}) | "
                      f"{r['c3'][0]:.6e}, {r['c3'][1]:.6e} | {r['sigma_c'] / A['sigma0']:.4f} |")
        L_.append("")
    if B is not None:
        L_ += ["## B. The split states S_q × S_{n−q} near the threshold (the middle class, n = 10)", "",
               "| μ | q | t* | Newton's drift | group positions (ξ) | predicted | eigenvalues on the separations along v (/μ) | predicted (/μ) | unstable eigenvalues of the full Jacobian |",
               "|---|---|---|---|---|---|---|---|---|"]
        for b in sorted(B, key=lambda b: -b['mu']):
            for r in b['rows']:
                L_.append(f"| {fmt(b['mu'])} | {r['q']} | {fmt(r['t'])} | {fmt(r['drift'])} | {r['groups'][0]:+.4e}, {r['groups'][1]:+.4e} | "
                          f"{r['pred_groups'][0]:+.4e}, {r['pred_groups'][1]:+.4e} | {', '.join(f'{x / b['mu']:+.4f}' for x in r['ev'])} | "
                          f"{', '.join(f'{x / b['mu']:+.3f}' for x in r['pred'])} | {r['n_unstable']} |")
        L_.append("")
    if C is not None:
        L_ += ["## C. Who leaves, when and where (jitter 10⁻³, 10⁵ slow units)", "",
               f"An excursion: a particle more than {LEAVE:g} from the median of the class's ξ (far inside the kernel's range R = 1). The predictions integrate "
               "the reduced equation from the measured ξ(0) with the quadratic term c₃, without it (linear) and with c₃'s sign flipped; the pieces at the "
               "end are the class's single-linkage pieces at the range R, with their distances from the main piece.", "",
               "| μ | seed | ranking by sign(c₃)ξ(0) (first three) and their values | excursions in order (side) | first at s | quadratic: first (particle) | linear | flipped c₃ | pieces beyond the range at the end (distance) | the main piece: particles, margin |",
               "|---|---|---|---|---|---|---|---|---|---|"]
        tally = dict(q=0, lin=0, flip=0, n=0, err=[], err_lin=[], side=0)
        for r in sorted(C, key=lambda r: (-r['mu'], r['seed'])):
            sc = np.sign(r['c3'])
            rk = r['rank'][:3]
            left = ', '.join(f"{i}({'+' if sd > 0 else '−'})" for i, (s_, sd) in r['left']) or '—'
            s_lin, i_lin = reduced_crossing(r['xi0'], r['mu'], 0.0, LEAVE, 4 * S_C)
            s_fl, i_fl = reduced_crossing(r['xi0'], r['mu'], -r['c3'], LEAVE, 4 * S_C)
            pcs = '; '.join(f"{{{', '.join(str(i) for i in p_)}}} ({fmt(d)})" for p_, d in r.get('pieces', [])) or '—'
            if r['left']:
                i1, (s1, sd1) = r['left'][0]
                tally['n'] += 1; tally['q'] += int(i1 == r['i_pred']); tally['lin'] += int(i1 == i_lin); tally['flip'] += int(i1 == i_fl)
                tally['err'].append(abs(r['s_pred'] / s1 - 1)); tally['err_lin'].append(abs(s_lin / s1 - 1)); tally['side'] += int(sd1 == sc)
            L_.append(f"| {fmt(r['mu'])} | {r['seed']} | {', '.join(str(int(i)) for i in rk)} ({', '.join(f'{sc * r['xi0'][i]:.2e}' for i in rk)}) | {left} | "
                      f"{fmt(r['left'][0][1][0]) if r['left'] else '—'} | {fmt(r['s_pred'])} ({r['i_pred']}) | {fmt(s_lin)} ({i_lin}) | {fmt(s_fl)} ({i_fl}) | {pcs} | "
                      f"{r['n_main']}, {fmt(r['rem_margin'])} |")
        L_ += ["", f"The first excursion's particle: the quadratic equation's in {tally['q']} of {tally['n']} runs, the linear one's in {tally['lin']}, the flipped one's in "
               f"{tally['flip']}; its time to {100 * max(tally['err']):.1f}% (quadratic) and {100 * max(tally['err_lin']):.1f}% (linear) at most; on the side sign(c₃)v in "
               f"{tally['side']} of {tally['n']}.", ""]
    if D is not None:
        L_ += ["## D. The cascade under Part G's schedule: the replicator weight law and the mass-free leak", "",
               "| law | seed | σ/σ₀: clusters (fit error) | at the schedule's end (σ₀/64; the mass-free leak σ₀/256, Part F): particles per class | centres | masses |", "|---|---|---|---|---|---|"]
        rows = list(D)
        if G_ref is not None:
            rows = [dict(law='cohesion leak (Part G)', seed=g['seed'], stages=g['stages']) for g in G_ref] + rows
        for r in rows:
            st = r['stages']
            L_.append(f"| {r['law']} | {r['seed']} | " + ', '.join(f"{s_['f']:g}: {s_['K']} ({fmt(s_['fit'])})" for s_ in st) + " | "
                      + f"{', '.join(str(int(x)) for x in st[-1]['n'])} | {', '.join(f'{x:+.2f}' for x in st[-1]['c'])} | "
                      + (', '.join(f'{x:.3f}' for x in st[-1]['masses']) if 'masses' in st[-1] else '—') + " |")
        L_.append("")
    L_.append("Numbers are produced by `experiments/exp_split_bifurcation.py`; nothing is edited by hand.")
    return '\n'.join(L_) + '\n'


def d_rows(d):
    """Part D's rows with the mass-free ones replaced by F's (to σ₀/256) when present, and E's."""
    rows = [r for r in (list(d['D']) if 'D' in d else []) if not ('F' in d and r['law'] == 'mass-free')]
    return rows + (list(d['F']) if 'F' in d else []) + (list(d['E']) if 'E' in d else []) or None


def g_reference():
    """Part G of exp_cohesion_leak.py (the cohesion leak itself under the same schedule), from its raw file."""
    z = E.load(os.path.join(E.OUT, 'cohesion_leak_raw.npz'))
    if 'G' not in z:
        return None
    out = []
    for g in z['G']:
        st = []
        for s_ in g['stages']:
            cl = np.asarray(s_['clusters'])
            o = np.argsort(cl[:, 2])
            st.append(dict(f=float(s_['frac']), K=int(s_['K_end']), fit=float(s_['fit']), n=np.rint(cl[o, 0] * COUNTS.sum()).astype(int), c=cl[o, 2],
                           masses=cl[o, 0]))
        out.append(dict(seed=int(g['seed']), stages=st))
    return out


def main():
    raw = os.path.join(OUT, 'split_bifurcation_raw.npz')
    md = os.path.join(OUT, 'split_bifurcation.md')
    args = sys.argv[1:]
    old = E.load(raw)
    if '--from-raw' in args:
        open(md, 'w').write(report(old.get('A'), old.get('B'), old.get('C'), d_rows(old), g_reference()))
        return
    parts = set('ABCDEF')
    if '--only' in args:
        parts = set(args[args.index('--only') + 1].split(','))
    from multiprocessing import Pool
    res = dict(old)
    with Pool(N_WORKERS) as pool:
        asyncs = {}
        if 'B' in parts:
            asyncs['B'] = pool.map_async(b_job, MUS_B, chunksize=1)
        if 'C' in parts:
            asyncs['C'] = pool.map_async(c_job, [(mu, s) for mu in MUS_C for s in SEEDS_C], chunksize=1)
        if 'D' in parts:
            asyncs['D'] = pool.map_async(d_job, [(law, s) for law in ('replicator', 'mass-free') for s in SEEDS_D], chunksize=1)
        if 'E' in parts:
            asyncs['E'] = pool.map_async(d_job, [('replicator 1', s) for s in SEEDS_D], chunksize=1)
        if 'F' in parts:                                                     # round 25: the mass-free leak to σ₀/256 (replaces D's mass-free rows)
            asyncs['F'] = pool.map_async(d_job, [('mass-free', s) for s in SEEDS_D], chunksize=1)
        if 'A' in parts:
            res['A'] = part_A(); E.save(raw, **res)
        for k in 'BCDEF':
            if k in parts:
                res[k] = asyncs[k].get(); E.save(raw, **res)
    open(md, 'w').write(report(res.get('A'), res.get('B'), res.get('C'), d_rows(res), g_reference()))


if __name__ == '__main__':
    main()
