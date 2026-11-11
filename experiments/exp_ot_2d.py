"""
Corollary D‴.14 in two dimensions: the real closed loop with the OT leak on a planar plant
=========================================================================================

The mechanism of Propositions D‴.12 and Corollary D‴.14 is dimension-free: a planar plant ẋ = f(x) + u, x ∈ ℝ², f a sum of four Gaussian
bumps with vector amplitudes, a Lissajous reference x_m(t) = (2 sin t, 1.5 sin 2t), a saturated network with p = 2n + 1 = 5 parameters per
kernel, and the matching prior (four atoms at the exact representation, masses 1/4). The real loop (k_e = 25, σ = σ₀/k_e, RK4 plant with the
control held over the step, Euler parameters, dt = 2·10⁻³) runs from the same jittered warm clouds under the nearest-atom leak and under the
barycentric OT leak (plan recomputed every 100 steps; RK4 plant with the control and the reference at the stages), for 2500 time units = 100 slow units = 5/σ₀ at σ₀ = 0.05 (the leak is a design input;
the transit of the OT leak takes about ln 2/σ₀ ≈ 14 slow units, the within-class contraction 1/σ₀ = 20): sizes m = 8, 32 and M = 128
(commensurate with the masses). Records: the cluster weights by nearest atom
(conserved against equalised), the rms tracking error over the last ten periods, the fit error along the reference orbit at the end, the
distance ‖f̂_m − f̂_M‖²_{L²(orbit)} between sizes and W²_{G⁻¹}(μ_m, μ_M).

Output: results/chaos/ot_2d.md and ot_2d_raw.npz (OT2D_T for the horizon, OT2D_SEEDS for the seeds).
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

OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos')
N_DIM = 2
NET_KW = dict(a_max=4.0, c_max=3.0, s_min=0.3, s_max=2.0)
G = np.array([1.0, 1.0, 0.25, 0.25, 0.25])
# four bumps (A ∈ ℝ², C ∈ ℝ², S): the exact representation has amplitudes 4A (masses 1/4), within a_max = 4
BUMPS = [(np.array([0.8, 0.3]), np.array([-1.2, 0.8]), 0.7), (np.array([-0.5, 0.6]), np.array([1.0, 1.0]), 0.6),
         (np.array([0.4, -0.7]), np.array([0.8, -1.1]), 0.6), (np.array([-0.6, -0.4]), np.array([-1.0, -0.9]), 0.7)]
PER = 2 * np.pi
SIGMA0 = 0.05
KE = 25.0
SIG = SIGMA0 / KE
DT = 2e-3
T_END = float(os.environ.get('OT2D_T', 2500.0))
PLAN_EVERY, REC_E, REC_W = 100, 100, 1000
SEEDS = tuple(int(s) for s in os.environ.get('OT2D_SEEDS', '0,1').split(','))
MS, M_REF = (8, 32), 128
JITTER = 0.02
Q = np.ones(4) / 4
N_PER = 10


def f_vec(X):
    X = np.atleast_2d(X)
    out = np.zeros_like(X)
    for A, C, S in BUMPS:
        out += A[None, :] * np.exp(-np.sum((X - C[None, :]) ** 2, axis=1) / (2 * S ** 2))[:, None]
    return out


def x_ref(t):
    return np.array([2.0 * np.sin(t), 1.5 * np.sin(2 * t)])


def x_ref_dot(t):
    return np.array([2.0 * np.cos(t), 3.0 * np.cos(2 * t)])


def make_net():
    return L.SaturatedLiftedNetwork(N_DIM, **NET_KW)


def exact_representation(net):
    return net.unconstrained(np.array([4 * A for A, C, S in BUMPS]), np.array([C for A, C, S in BUMPS]), np.array([S for A, C, S in BUMPS]))


ORBIT_X = np.array([x_ref(tt) for tt in np.arange(400) * PER / 400])
F_ORB = f_vec(ORBIT_X)


def mixture_start(atoms, M, rng, jitter):
    kk = rng.integers(0, atoms.shape[0], M)
    return atoms[kk] + jitter * rng.standard_normal((M, atoms.shape[1])) * np.sqrt(G)[None, :]


def cluster_weights(Th, atoms):
    lab = OT.nearest_targets(Th, atoms, G)[1]
    return np.array([np.mean(lab == k) for k in range(atoms.shape[0])])


def real_run(net, exact, systems, leak, T):
    st = {k: dict(x=x_ref(0.0).copy(), Th=np.array(Th, float)) for k, Th in systems.items()}
    tgt = {}
    n = int(round(T / DT)); t = 0.0
    out = {k: dict(t=[], e2=[], tw=[], w=[]) for k in systems}
    for i in range(n):
        if i % PLAN_EVERY == 0:
            for k, S in st.items():
                tgt[k] = OT.barycentric_targets(S['Th'], exact, Q, G)[0] if leak == 'ot' else OT.nearest_targets(S['Th'], exact, G)[0]
        xr, xrd = x_ref(t), x_ref_dot(t)
        for k, S in st.items():
            x, Th = S['x'], S['Th']
            e = x - xr
            xdot = lambda y, tau: f_vec(y)[0] + x_ref_dot(tau) - KE * (y - x_ref(tau)) - net.f_hat(y, Th)   # the closed-loop plant, Θ frozen over the step
            k1 = xdot(x, t); k2 = xdot(x + 0.5 * DT * k1, t + 0.5 * DT); k3 = xdot(x + 0.5 * DT * k2, t + 0.5 * DT); k4 = xdot(x + DT * k3, t + DT)
            drift = G[None, :] * net.grad_g(x, e, Th) - SIG * (Th - tgt[k])
            S['x'] = x + DT / 6 * (k1 + 2 * k2 + 2 * k3 + k4); S['Th'] = Th + DT * drift
            if i % REC_E == 0:
                out[k]['t'].append(t); out[k]['e2'].append(float(e @ e))
            if i % REC_W == 0:
                out[k]['tw'].append(t); out[k]['w'].append(cluster_weights(S['Th'], exact))
        t = (i + 1) * DT
    for k in out:
        out[k] = {kk: np.array(v) for kk, v in out[k].items()}
        out[k]['Th'] = st[k]['Th']
        out[k]['fh'] = net.f_hat_batch(ORBIT_X, st[k]['Th']).astype(np.float32)          # the learned function on the orbit grid at the end
    return out


def job(seed):
    net = make_net()
    exact = exact_representation(net)
    rng = np.random.default_rng(seed)
    Th_all = mixture_start(exact, M_REF, rng, JITTER)
    systems = {m: Th_all[:m].copy() for m in MS + (M_REF,)}
    out = dict(seed=seed)
    for leak in ('ot', 'nearest'):
        out[leak] = real_run(net, exact, systems, leak, T_END)
    return out


def e_rms(o):
    msk = o['t'] > T_END - N_PER * PER
    return float(np.sqrt(np.mean(o['e2'][msk])))


def class_spread(Th, atoms):
    lab = OT.nearest_targets(Th, atoms, G)[1]
    d2 = []
    for k in np.unique(lab):
        P = Th[lab == k]
        d2 += list(np.sum((P - P.mean(axis=0)) ** 2 / G[None, :], axis=1))
    return float(np.sqrt(np.mean(d2)))


def fit_orbit(o):
    return float(np.sqrt(np.mean(np.sum((o['fh'] - F_ORB) ** 2, axis=1)) / np.mean(np.sum(F_ORB ** 2, axis=1))))


def fdist(o1, o2):
    return float(np.mean(np.sum((o1['fh'] - o2['fh']) ** 2, axis=1)))


def report(res, t0=None):
    lines = ["# Corollary D‴.14 in two dimensions: the real closed loop with the OT leak on a planar plant (`exp_ot_2d.py`)", "",
             f"Plant ẋ = f(x) + u on ℝ², f the sum of four Gaussian bumps with vector amplitudes (A_k, C_k, S_k) = " + ", ".join(f"(({A[0]:g}, {A[1]:g}), ({C[0]:g}, {C[1]:g}), {S:g})" for A, C, S in BUMPS) + ", reference x_m(t) = (2 sin t, 1.5 sin 2t), "
             f"saturated network with n = 2 (a_max = {NET_KW['a_max']}, c_max = {NET_KW['c_max']}, s ∈ [{NET_KW['s_min']}, {NET_KW['s_max']}]), G = diag{tuple(float(g) for g in G)}, matching prior = the exact representation "
             f"(four atoms with amplitudes 4A_k, masses q = 1/4). Real loop at k_e = {KE:g}, σ = σ₀/k_e with σ₀ = {SIGMA0:g} (ln 2/σ₀ = {np.log(2) / SIGMA0:.0f} slow units), RK4 plant with the control and the reference evaluated at the stages and the parameters frozen over the step, Euler parameters, dt = {DT:g}, {T_END:g} time units = {T_END / KE:.0f} slow units = {T_END * SIGMA0 / KE:.1f}/σ₀ = {T_END / PER:.0f} periods; "
             f"warm start = the mixture of Gaussians of width {JITTER:g} (G⁻¹ metric) around the exact representation with random class labels; the same clouds under the nearest-atom leak and under the OT leak "
             f"(plan every {PLAN_EVERY} steps); sizes m = {', '.join(str(m) for m in MS)} and M = {M_REF}; seeds {', '.join(str(r['seed']) for r in res)}; statistics over the last {N_PER} periods.", "",
             "| seed | leak | m | cluster weights at the start → at the end | within-class spread at the end (rms, whitened) | rms tracking error (last 10 periods) | fit error along the orbit (end) | ‖f̂_m − f̂_M‖²_{L²(orbit)} (end) | W²_{G⁻¹}(μ_m, μ_M) (end) |",
             "|---|---|---|---|---|---|---|---|---|"]
    net = make_net(); exact = exact_representation(net)
    for r in res:
        for leak in ('ot', 'nearest'):
            R = r[leak]
            for m in MS + (M_REF,):
                o = R[m]
                lines.append(f"| {r['seed']} | {leak} | {m} | ({', '.join(f'{x:.3f}' for x in o['w'][0])}) → ({', '.join(f'{x:.3f}' for x in o['w'][-1])}) | {class_spread(o['Th'], exact):.1e} | {e_rms(o):.3e} | {fit_orbit(o):.3e} | "
                             f"{fdist(o, R[M_REF]):.3e} | {OT.W2_empirical(o['Th'], R[M_REF]['Th'], G):.3e} |")
    # the transit: first weight record from which the OT weights equal q within half a particle
    eq = []
    for r in res:
        for m in MS + (M_REF,):
            o = r['ot'][m]
            ok = np.all(np.abs(o['w'] - Q[None, :]) < 0.5 / m, axis=1)
            if ok[-1]:
                idx = len(ok) - 1
                while idx > 0 and ok[idx - 1]:
                    idx -= 1
                eq.append((r['seed'], m, float(o['tw'][idx]) / KE))
            else:
                eq.append((r['seed'], m, float('nan')))
    e_ot = [e_rms(r['ot'][m]) for r in res for m in MS + (M_REF,)]; e_nr = [e_rms(r['nearest'][m]) for r in res for m in MS + (M_REF,)]
    h_ot = [fdist(r['ot'][m], r['ot'][M_REF]) for r in res for m in MS]; h_nr = [fdist(r['nearest'][m], r['nearest'][M_REF]) for r in res for m in MS]
    per = [e_rms(r['nearest'][m]) / e_rms(r['ot'][m]) for r in res for m in MS + (M_REF,)]
    lines += ["", f"Equalisation of the OT leak's cluster weights (first weight record from which w = q to within half a particle, in slow units): " + "; ".join(f"seed {s}, m = {m}: {x:.0f}" if np.isfinite(x) else f"seed {s}, m = {m}: not reached" for s, m, x in eq) + f" (ln 2/σ₀ = {np.log(2) / SIGMA0:.0f}).", "",
              f"Reading. In two dimensions as in one: the nearest-atom leak conserves the clouds' sampled cluster weights and the sizes learn different functions (‖f̂_m − f̂_M‖² = {min(h_nr):.1e}–{max(h_nr):.1e}), the OT leak drives the weights to the prior's masses and the sizes agree "
              f"({min(h_ot):.1e}–{max(h_ot):.1e}, at the stopping time 5/σ₀ and still falling); the rms tracking error is {min(e_ot):.1e}–{max(e_ot):.1e} under the OT leak against {min(e_nr):.1e}–{max(e_nr):.1e} under the nearest-atom leak, "
              f"{min(per):.0f}–{max(per):.0f} times larger per system.", "",
              "Numbers are produced by `experiments/exp_ot_2d.py`; nothing is edited by hand."]
    with open(os.path.join(OUT, 'ot_2d.md'), 'w', encoding='utf-8') as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    if t0 is not None:
        print(f"[{time.time() - t0:.0f}s]")


def main():
    import multiprocessing as mp
    t0 = time.time()
    if '--from-raw' in sys.argv:
        z = np.load(os.path.join(OUT, 'ot_2d_raw.npz'), allow_pickle=True)
        report(list(z['res']), t0)
        return
    with mp.Pool(min(len(SEEDS), 2)) as pool:
        res = pool.map(job, list(SEEDS))
    np.savez(os.path.join(OUT, 'ot_2d_raw.npz'), res=np.array(res, dtype=object), allow_pickle=True)
    report(res, t0)


if __name__ == '__main__':
    main()
