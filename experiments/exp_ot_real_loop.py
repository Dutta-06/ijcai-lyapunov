"""
Corollary D‴.14: the real closed loop with the OT leak — the periodic orbit is the lift of the weighted K₀-atom loop's, the same for every
commensurate size
====================================================================================================================================

The real closed loop of Theorem J (k_e = 25, σ = σ₀/k_e, RK4 plant with the control and the reference evaluated at the stages and the
parameters frozen over the step, Euler parameters, dt = 10⁻³) on the
three-bump instance, with the nearest-atom leak and with the barycentric OT leak toward the matching prior (plan recomputed every 100 steps),
started from the clouds at the end of the averaged runs of `exp_ot_leak.py` (commensurate set: m = 15, 60 and M = 510; the OT
runs' clouds have the prior's masses, the nearest-atom runs' clouds their conserved sampled weights) and from the collapsed lifts — the
three atoms of the exact representation for the OT leak (the weighted 3-atom loop with q = (1/3, 1/3, 1/3)), and the replicated weighted
stationary configurations Θ*(w^m) for the nearest-atom leak. Over the last five periods: the learned function along the reference orbit
against the reference's and against the lift's, the fit error along the orbit, the rms tracking error, W²_{G⁻¹} between configurations.

Output: results/chaos/ot_real_loop.md and ot_real_loop_raw.npz (ORL_T for the horizon, ORL_SEEDS for the seeds).
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
from src.core import chaos_learning as CL
from src.core import prior_design as PD
from src.core import ot_leak as OT

OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos')
NET_KW = dict(a_max=4.0, c_max=3.0, s_min=0.3, s_max=2.0)
G = np.array([1.0, 0.25, 0.25])
BUMPS = [(1.0, -1.5, 0.5), (-0.7, 0.0, 0.4), (0.8, 1.5, 0.5)]
PER = 2 * np.pi
N_Q = 256
SIGMA_CRIT = 0.047162567584665846
SIGMA0 = 0.15 * SIGMA_CRIT
KE = 25.0
SIG = SIGMA0 / KE
DT = 1e-3
T_END = float(os.environ.get('ORL_T', 100.0))
PLAN_EVERY, REC_EVERY = 100, 100
SEEDS = tuple(int(s) for s in os.environ.get('ORL_SEEDS', '0,1').split(','))
MS, M_REF = (15, 60), 510
N_PER = 5                                                           # the last N_PER periods are the statistics window
Q = np.ones(3) / 3


def f_vec(X):
    return sum(A * np.exp(-(X - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS)


def x_ref(t):
    return np.array([2.0 * np.sin(t)])


def x_ref_dot(t):
    return np.array([2.0 * np.cos(t)])


def make_net():
    return L.SaturatedLiftedNetwork(1, **NET_KW)


def exact_representation(net):
    return net.unconstrained(np.array([[3 * A] for A, C, S in BUMPS]), np.array([[C] for A, C, S in BUMPS]), np.array([S for A, C, S in BUMPS]))


ORBIT_X = np.array([x_ref(tt)[0] for tt in np.arange(400) * PER / 400])[:, None]
F_ORB = f_vec(ORBIT_X[:, 0])


def real_run(net, loop, systems, leak, T):
    """The real closed loop for a dict name → Θ₀ (all sharing the plant's initial state x_ref(0)); returns name → records (t, f̂ on the orbit
    grid, e², the final configuration)."""
    st = {k: dict(x=x_ref(0.0).copy(), Th=np.array(Th, float)) for k, Th in systems.items()}
    tgt = {}
    n = int(round(T / DT)); t = 0.0
    out = {k: dict(t=[], fh=[], e2=[]) for k in systems}
    for i in range(n):
        if i % PLAN_EVERY == 0:
            for k, S in st.items():
                tgt[k] = loop.leak_target(S['Th']) if leak == 'ot' else loop.nearest.leak_target(S['Th'])
        xr, xrd = x_ref(t), x_ref_dot(t)
        for k, S in st.items():
            x, Th = S['x'], S['Th']
            e = x - xr
            xdot = lambda y, tau: f_vec(y) + x_ref_dot(tau) - KE * (y - x_ref(tau)) - net.f_hat(y, Th)   # the closed-loop plant, Θ frozen over the step
            k1 = xdot(x, t); k2 = xdot(x + 0.5 * DT * k1, t + 0.5 * DT); k3 = xdot(x + 0.5 * DT * k2, t + 0.5 * DT); k4 = xdot(x + DT * k3, t + DT)
            drift = G[None, :] * net.grad_g(x, e, Th) - SIG * (Th - tgt[k])
            S['x'] = x + DT / 6 * (k1 + 2 * k2 + 2 * k3 + k4); S['Th'] = Th + DT * drift
            if i % REC_EVERY == 0:
                out[k]['t'].append(t); out[k]['e2'].append(float(e[0] ** 2))
                out[k]['fh'].append(net.f_hat_batch(ORBIT_X, S['Th'])[:, 0])
        t = (i + 1) * DT
    for k in out:
        out[k] = {kk: np.array(v) for kk, v in out[k].items()}
        out[k]['Th'] = st[k]['Th']
        keep = out[k]['t'] > T - N_PER * PER                               # f̂ on the orbit grid is kept on the statistics window only (raw-file size)
        out[k]['fh'] = out[k]['fh'][keep].astype(np.float32); out[k]['t_fh'] = out[k]['t'][keep]      # float32: 10⁻⁷ relative, below the printed digits
    return out


def lift_nearest(loop, exact, w, m):
    """The collapsed lift of the nearest-atom loop with the realised weights w: round(m w_k) particles at ϑ*_k(w)."""
    atoms, res = PD.weighted_stationary_atoms_prior(loop.nearest, w, exact, w0=Q, n_cont=4)
    n = np.round(m * np.asarray(w)).astype(int)
    assert n.sum() == m
    return np.repeat(atoms, n, axis=0), atoms, res


def job(seed):
    net = make_net()
    exact = exact_representation(net)
    loop = OT.OTPriorLoop(net, f_vec, x_ref, PER, G, SIGMA0, exact, n_quad=N_Q)
    z = np.load(os.path.join(OUT, 'ot_leak_raw.npz'), allow_pickle=True)
    runs = list(z['runs'])
    pick = lambda leak: next(r for r in runs if r['set'] == 'commensurate' and r['leak'] == leak and r['seed'] == seed)
    R = {leak: pick(leak) for leak in ('ot', 'nearest')}
    out = dict(seed=seed)
    # (A) the OT leak from the equalised, contracted clouds, and the lift (the exact representation, q = 1/3)
    sysA = {f'm{m}': R['ot'][f'final{m}'] for m in MS}; sysA[f'm{M_REF}'] = R['ot'][f'final{M_REF}']; sysA['lift'] = exact.copy()
    out['ot'] = real_run(net, loop, sysA, 'ot', T_END)
    out['ot_w0'] = {m: CL.cluster_weights(R['ot'][f'final{m}'], exact, G)[1] for m in MS + (M_REF,)}
    # (B) the nearest-atom leak from the contracted clouds with their conserved weights, and the collapsed lifts at Θ*(w^m)
    sysB, lifts = {}, {}
    for m in MS + (M_REF,):
        Th = R['nearest'][f'final{m}']; sysB[f'm{m}'] = Th
        w = CL.cluster_weights(Th, exact, G)[1]
        sysB[f'lift{m}'], atoms, res = lift_nearest(loop, exact, w, m)
        lifts[m] = dict(w=w, atoms=atoms, res=res)
    out['nearest'] = real_run(net, loop, sysB, 'nearest', T_END)
    out['nearest_lifts'] = lifts
    return out


def last_mask(t):
    return t > T_END - N_PER * PER


def fdist(o1, o2):
    msk = last_mask(o1['t_fh'])
    return float(np.mean(np.mean((o1['fh'][msk] - o2['fh'][msk]) ** 2, axis=1)))


def fit_orbit(o):
    msk = last_mask(o['t_fh'])
    return float(np.sqrt(np.mean(np.mean((o['fh'][msk] - F_ORB[None, :]) ** 2, axis=1)) / np.mean(F_ORB ** 2)))


def e_rms(o):
    msk = last_mask(o['t'])
    return float(np.sqrt(np.mean(o['e2'][msk])))


def class_spread(Th, exact):
    """rms whitened distance of the particles to their class mean (classes by nearest prior atom): the within-class spread."""
    lab = CL.cluster_weights(Th, exact, G)[0]
    d2 = []
    for k in np.unique(lab):
        P = Th[lab == k]
        d2 += list(np.sum((P - P.mean(axis=0)) ** 2 / G[None, :], axis=1))
    return float(np.sqrt(np.mean(d2)))


def report(res, t0=None):
    net = make_net()
    exact = exact_representation(net)
    loop = OT.OTPriorLoop(net, f_vec, x_ref, PER, G, SIGMA0, exact, n_quad=N_Q)
    lines = ["# Corollary D‴.14: the real closed loop with the OT leak — the orbit is the lift of the weighted K₀-atom loop's, the same for every commensurate size (`exp_ot_real_loop.py`)", "",
             f"The real closed loop of Theorem J on the three-bump instance (k_e = {KE:g}, σ = σ₀/k_e with σ₀ = 0.15σ_crit = {SIGMA0:.5f}, RK4 plant with the control and the reference evaluated at the stages and the parameters frozen over the step, Euler parameters, dt = {DT:g}), "
             f"{T_END:g} time units = {T_END / PER:.1f} periods = {T_END / KE:.1f} slow units; statistics over the last {N_PER} periods on the 400-point grid of the reference orbit. "
             f"Starts: the clouds at the end of `exp_ot_leak.py`'s averaged runs (contracted where the realised weights are admissible) (commensurate set, m = {', '.join(str(m) for m in MS)}, M = {M_REF}; the OT runs' clouds carry the prior's masses q = (1/3, 1/3, 1/3), "
             f"the nearest-atom runs' clouds their conserved sampled weights), and the collapsed lifts: the three atoms of the exact representation for the OT leak (the weighted 3-atom loop with weights q, "
             f"itself the m = 3 loop), and the replicated weighted stationary configurations Θ*(w^m) for the nearest-atom leak. Seeds {', '.join(str(r['seed']) for r in res)}; the leak target (the nearest atom, or the barycentric "
             f"projection of the vertex plan) recomputed every {PLAN_EVERY} steps.", ""]
    lines += ["## 1. The OT leak: every size on the same orbit", "",
              "| seed | system | cluster weights at the start | within-class spread at the end (rms, whitened) | ‖f̂_m − f̂_M‖²_{L²(orbit)} (M = 510) | ‖f̂_m − f̂_lift‖² | fit error along the orbit | rms tracking error | W²_{G⁻¹}(μ_m, μ_M) at the end | W² to the lift |",
              "|---|---|---|---|---|---|---|---|---|---|"]
    for r in res:
        A = r['ot']
        for k in [f'm{m}' for m in MS + (M_REF,)] + ['lift']:
            m = int(k[1:]) if k != 'lift' else 3
            w = r['ot_w0'][m] if k != 'lift' else Q
            lines.append(f"| {r['seed']} | {k if k != 'lift' else 'lift (m = 3)'} | ({', '.join(f'{x:.4f}' for x in w)}) | {class_spread(A[k]['Th'], exact):.1e} | {fdist(A[k], A[f'm{M_REF}']):.3e} | {fdist(A[k], A['lift']):.3e} | {fit_orbit(A[k]):.3e} | {e_rms(A[k]):.3e} | "
                         f"{OT.W2_empirical(A[k]['Th'], A[f'm{M_REF}']['Th'], G):.3e} | {OT.W2_empirical(A[k]['Th'], A['lift']['Th'], G):.3e} |")
    lines += ["", "## 2. The nearest-atom leak: every size with admissible weights on the orbit of its own weights", "",
              "| seed | m | conserved weights w^m | splitting margin σ₀ − max κ_R at w^m (Proposition D‴.5(i); ≤ 0: the lift is a saddle orbit) | within-class spread at the end (rms, whitened) | ‖f̂_m − f̂_M‖² | ‖f̂_m − f̂_{lift(w^m)}‖² | ‖f̂_{lift(w^m)} − f̂_{lift(w^M)}‖² | fit error along the orbit (m-system / its lift) | rms tracking error (m-system / its lift) | W²(μ_m, μ_M) | W² to its lift |",
              "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    margins = {}
    for r in res:
        B = r['nearest']
        for m in MS + (M_REF,):
            k, lk = f'm{m}', f'lift{m}'
            mg = PD.splitting_margin_prior(loop.nearest, r['nearest_lifts'][m]['w'], exact)[0]
            margins[(r['seed'], m)] = mg
            lines.append(f"| {r['seed']} | {m} | ({', '.join(f'{x:.4f}' for x in r['nearest_lifts'][m]['w'])}) | {mg:+.5f}{' (saddle)' if mg <= 0 else ''} | {class_spread(B[k]['Th'], exact):.1e} | {fdist(B[k], B[f'm{M_REF}']):.3e} | {fdist(B[k], B[lk]):.3e} | {fdist(B[lk], B[f'lift{M_REF}']):.3e} | "
                         f"{fit_orbit(B[k]):.3e} / {fit_orbit(B[lk]):.3e} | {e_rms(B[k]):.3e} / {e_rms(B[lk]):.3e} | {OT.W2_empirical(B[k]['Th'], B[f'm{M_REF}']['Th'], G):.3e} | {OT.W2_empirical(B[k]['Th'], B[lk]['Th'], G):.3e} |")
    # summary numbers
    ot_mM = [fdist(r['ot'][f'm{m}'], r['ot'][f'm{M_REF}']) for r in res for m in MS]
    nr_mM = [fdist(r['nearest'][f'm{m}'], r['nearest'][f'm{M_REF}']) for r in res for m in MS]
    nr_ll = [fdist(r['nearest'][f'lift{m}'], r['nearest'][f'lift{M_REF}']) for r in res for m in MS]
    nr_ml = [fdist(r['nearest'][f'm{m}'], r['nearest'][f'lift{m}']) for r in res for m in MS + (M_REF,)]
    ot_ml = [fdist(r['ot'][f'm{m}'], r['ot']['lift']) for r in res for m in MS + (M_REF,)]
    e_ot = [e_rms(r['ot'][f'm{m}']) for r in res for m in MS + (M_REF,)]; e_nr = [e_rms(r['nearest'][f'm{m}']) for r in res for m in MS + (M_REF,)]
    ratio_ll = [a / b for a, b in zip(nr_mM, nr_ll)]
    adm = [(r['seed'], m) for r in res for m in MS if margins[(r['seed'], m)] > 0 and margins[(r['seed'], M_REF)] > 0]
    sad = [(r['seed'], m) for r in res for m in MS + (M_REF,) if margins[(r['seed'], m)] <= 0]
    ratio_adm = [fdist(r['nearest'][f'm{m}'], r['nearest'][f'm{M_REF}']) / fdist(r['nearest'][f'lift{m}'], r['nearest'][f'lift{M_REF}']) for r in res for m in MS if (r['seed'], m) in adm]
    ml_adm = [fdist(r['nearest'][f'm{m}'], r['nearest'][f'lift{m}']) for r in res for m in MS + (M_REF,) if (r['seed'], m) not in sad]
    ml_sad = [fdist(r['nearest'][f'm{m}'], r['nearest'][f'lift{m}']) for r in res for m in MS + (M_REF,) if (r['seed'], m) in sad]
    sp_sad = [class_spread(r['nearest'][f'm{m}']['Th'], exact) for r in res for m in MS + (M_REF,) if (r['seed'], m) in sad]
    sp_adm = [class_spread(r['nearest'][f'm{m}']['Th'], exact) for r in res for m in MS + (M_REF,) if (r['seed'], m) not in sad]
    lines += ["", f"Reading. Under the OT leak the learned functions of m = 15, 60 and M = 510 agree along the orbit to ‖f̂_m − f̂_M‖² = {min(ot_mM):.1e}–{max(ot_mM):.1e} and with the 3-atom lift to {min(ot_ml):.1e}–{max(ot_ml):.1e} "
              f"(the within-class spread of the clouds, second order); under the nearest-atom leak they disagree by {min(nr_mM):.1e}–{max(nr_mM):.1e}. Where both realised weight vectors are admissible this disagreement is the distance "
              f"between the lifted orbits of the conserved weights to a ratio of {', '.join(f'{x:.4f}' for x in ratio_adm)} ({len(ratio_adm)} pair{'s' if len(ratio_adm) != 1 else ''}), and each admissible m-system is within {max(ml_adm):.1e} of its own lift (within-class spread {min(sp_adm):.1e}–{max(sp_adm):.1e}); "
              + (f"the {len(sad)} systems with inadmissible realised weights — " + ", ".join(f"seed {s}, m = {m} (margin {margins[(s, m)]:+.5f})" for s, m in sad) + f" — are not on their lifts, which are saddle orbits: their over-weighted class is split (within-class spread {min(sp_sad):.1e}–{max(sp_sad):.1e}), "
                 f"the distance to the lift {min(ml_sad):.1e}–{max(ml_sad):.1e} being the split's second-order effect on f̂, and the lift-to-lift ratio at those pairs is {', '.join(f'{a / b:.3f}' for a, b, (s, m) in [(fdist(r['nearest'][f'm{m}'], r['nearest'][f'm{M_REF}']), fdist(r['nearest'][f'lift{m}'], r['nearest'][f'lift{M_REF}']), (r['seed'], m)) for r in res for m in MS] if (s, m) not in adm)}. " if sad else "all realised weights are admissible. ")
              + f"The rms tracking error is {min(e_ot):.1e}–{max(e_ot):.1e} under the OT leak (the m-systems' remaining fit error over k_e; the collapsed lift tracks to {e_rms(res[0]['ot']['lift']):.1e}, the scheme's floor) "
              f"against {min(e_nr):.1e}–{max(e_nr):.1e} under the nearest-atom leak — the weight bias of the fit is a tracking error, which the mass constraint removes for every size at once.", "",
              "Numbers are produced by `experiments/exp_ot_real_loop.py`; nothing is edited by hand."]
    with open(os.path.join(OUT, 'ot_real_loop.md'), 'w', encoding='utf-8') as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    if t0 is not None:
        print(f"[{time.time() - t0:.0f}s]")


def main():
    t0 = time.time()
    if '--from-raw' in sys.argv:
        z = np.load(os.path.join(OUT, 'ot_real_loop_raw.npz'), allow_pickle=True)
        report(list(z['res']), t0)
        return
    res = [job(s) for s in SEEDS]
    np.savez(os.path.join(OUT, 'ot_real_loop_raw.npz'), res=np.array(res, dtype=object), allow_pickle=True)
    report(res, t0)


if __name__ == '__main__':
    main()
