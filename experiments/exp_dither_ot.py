"""
Remark D‴.16: dither and the mass constraint — the dithered averaged slow loop under the nearest-atom and the OT leak
====================================================================================================================

The averaged slow loop of Theorem D‴ with dither in slow time (Theorem J′'s leading order: dϑ_i = −G∇_{ϑ_i}J_m ds + √(2τ₀) G^{1/2} dW_i, whose
stationary law is the Gibbs measure ∝ exp(−J_m/τ₀)), on the three-bump instance with the matching prior (three atoms at the exact
representation, masses q = 1/3), m = 15 particles started with the cluster masses (6, 5, 4)/15, under the nearest-atom leak (𝔄_m^ν) and the
barycentric OT leak (𝔄_m^{OT}), at τ₀ = 0.005, 0.01, 0.02 and σ₀ = 0.15σ_crit.

Records the cluster masses by nearest atom along long runs: under the nearest-atom leak the masses change only by a particle hopping across a
Voronoi ridge (Arrhenius), and their stationary distribution is the Laplace approximation of the Gibbs measure over the count vectors n,
P(n) ∝ multinomial(n)·exp(−m𝒢*_{n/m}/τ₀)·Π_k det(σ₀I − R^G_k(n))^{−(n_k−1)/2}·det(H_n)^{−1/2} (the replicated weighted stationary configuration's
energy and its Hessian's block determinants, Proposition D‴.5(i)); under the OT leak the plan's masses are q by construction and a particle
that crosses a ridge stays assigned to its own atom, so the counts by nearest atom leave q only during reversible excursions, which the
report classifies (the crossing particle returns, or two particles exchange their assignments).

Output: results/chaos/dither_ot.md and dither_ot_raw.npz (DOT_SEEDS, DOT_SCALE to shorten the horizons).
"""
import os
import sys
import time
from math import lgamma, log
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
N_Q = 128
SIGMA_CRIT = 0.047162567584665846
SIGMA0 = 0.15 * SIGMA_CRIT
M = 15
COUNTS0 = (6, 5, 4)
TAUS = (0.005, 0.01, 0.02)
SCALE = float(os.environ.get('DOT_SCALE', 1.0))
S_END = {0.005: 50000.0 * SCALE, 0.01: 200000.0 * SCALE, 0.02: 50000.0 * SCALE}
DT, PLAN_EVERY, REC_EVERY = 0.05, 10, 200
SEEDS = tuple(int(s) for s in os.environ.get('DOT_SEEDS', '0,1').split(','))
Q = np.ones(3) / 3


def f_vec(X):
    return sum(A * np.exp(-(X - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS)


def x_ref(t):
    return np.array([2.0 * np.sin(t)])


def make_net():
    return L.SaturatedLiftedNetwork(1, **NET_KW)


def exact_representation(net):
    return net.unconstrained(np.array([[3 * A] for A, C, S in BUMPS]), np.array([[C] for A, C, S in BUMPS]), np.array([S for A, C, S in BUMPS]))


def make_loop(net):
    return OT.OTPriorLoop(net, f_vec, x_ref, PER, G, SIGMA0, exact_representation(net), n_quad=N_Q)


def run(leak, tau0, seed):
    net = make_net(); exact = exact_representation(net); loop = make_loop(net)
    rng = np.random.default_rng(seed)
    sqG = np.sqrt(G)
    Th = np.repeat(exact, COUNTS0, axis=0) + 0.02 * rng.standard_normal((M, 3)) * sqG[None, :]
    n = int(round(S_END[tau0] / DT)); labs = []; T = None
    for k in range(n):
        if k % PLAN_EVERY == 0:
            T = loop.leak_target(Th) if leak == 'ot' else loop.nearest.leak_target(Th)
        Th = Th + DT * loop.drift(Th, T) + np.sqrt(2 * tau0 * DT) * sqG[None, :] * rng.standard_normal((M, 3))
        if k % REC_EVERY == 0:
            labs.append(CL.cluster_weights(Th, exact, G)[0].astype(np.int8))
    return dict(leak=leak, tau0=tau0, seed=seed, labs=np.array(labs), s=np.arange(len(labs)) * REC_EVERY * DT)


def laplace_table(net, loop, exact):
    """For every count vector n (n_k ≥ 1): the weighted stationary point's energy 𝒢*_n, the per-class separation blocks' log-determinants and the
    collective block's; returns dict n → (log multinomial, 𝒢*, log det factor)."""
    out = {}
    sqG = np.sqrt(G)
    for n1 in range(1, M - 1):
        for n2 in range(1, M - n1):
            n3 = M - n1 - n2
            if n3 < 1:
                continue
            n = (n1, n2, n3); w = np.array(n) / M
            Tw, res = PD.weighted_stationary_atoms_prior(loop.nearest, w, exact, w0=Q, n_cont=4)
            a, c, s, diff, phi = loop.features(Tw); r = loop.F - phi @ (w * a)
            Gw = 0.5 * float(np.mean(r ** 2)) + 0.5 * SIGMA0 * float(np.sum(w * np.sum((Tw - exact) ** 2 / G[None, :], axis=1)))
            R = loop.residual_curvature(r, Tw); Rw = sqG[None, :, None] * R * sqG[None, None, :]
            logdet_sep = 0.0
            for k in range(3):
                ev = np.linalg.eigvalsh(SIGMA0 * np.eye(3) - Rw[k])
                if np.any(ev <= 0):
                    logdet_sep = float('inf'); break                     # a saddle: not a Gibbs mode (the Laplace term is absent)
                logdet_sep += (n[k] - 1) * float(np.sum(np.log(ev)))
            H = hess_weighted(loop, w, Tw)
            evH = np.linalg.eigvalsh(H)
            logdet_col = float(np.sum(np.log(evH))) if np.all(evH > 0) else float('inf')
            logmult = lgamma(M + 1) - sum(lgamma(x + 1) for x in n)
            out[n] = (logmult, Gw, logdet_sep, logdet_col)
    return out


def hess_weighted(loop, w, atoms, h=1e-6):
    K = atoms.shape[0]; d = 3 * K; v = atoms.ravel(); Jm = np.zeros((d, d))
    grad_fn = lambda V: PD.weighted_energy_grad_prior(loop.nearest, w, V.reshape(K, 3)).ravel()
    for j in range(d):
        e = np.zeros(d); e[j] = h
        Jm[:, j] = (grad_fn(v + e) - grad_fn(v - e)) / (2 * h)
    sw = np.repeat(np.sqrt(np.asarray(w, float)), 3); sg = np.tile(np.sqrt(G), K)
    Hs = (sw * sg)[:, None] * Jm * (sg / sw)[None, :]
    return 0.5 * (Hs + Hs.T)


def gibbs_prediction(table, tau0, with_hessian=True):
    logp = {}
    for n, (logmult, Gw, lds, ldc) in table.items():
        if with_hessian and (lds == float('inf') or ldc == float('inf')):
            continue
        logp[n] = logmult - M * Gw / tau0 - (0.5 * (lds + ldc) if with_hessian else 0.0)
    mx = max(logp.values()); Z = sum(np.exp(v - mx) for v in logp.values())
    return {n: float(np.exp(v - mx) / Z) for n, v in logp.items()}


def excursions(labs, balanced=(5, 5, 5)):
    """After the first balanced record: the maximal runs of records whose counts are not `balanced`, classified by the per-particle labels of
    the balanced records before and after the run — 'return' (the same labels) or 'exchange' (balanced again with some labels changed).
    Returns (n_return, n_exchange, off-balance fraction after the first balanced record, index of the first balanced record)."""
    counts = np.array([[int(np.sum(l == k)) for k in range(3)] for l in labs])
    bal = np.all(counts == np.array(balanced)[None, :], axis=1)
    if not bal.any():
        return 0, 0, float('nan'), -1
    i0 = int(np.argmax(bal))
    n_ret = n_exc = 0
    i = i0
    while i < len(bal):
        if bal[i]:
            i += 1
            continue
        j = i
        while j < len(bal) and not bal[j]:
            j += 1
        if j < len(bal):
            if np.array_equal(labs[i - 1], labs[j]):
                n_ret += 1
            else:
                n_exc += 1
        i = j
    return n_ret, n_exc, float(1.0 - bal[i0:].mean()), i0


def undone_changes(labs):
    """Record-level label changes, and how many of them are undone by the same particle at the next change (recrossings)."""
    L = labs.astype(int)
    ch = [i for i in range(1, len(L)) if np.any(L[i] != L[i - 1])]
    undone = 0
    for a, b in zip(ch[:-1], ch[1:]):
        moved = np.flatnonzero(L[a] != L[a - 1])
        if np.array_equal(L[b][moved], L[a - 1][moved]) and np.array_equal(np.flatnonzero(L[b] != L[b - 1]), moved):
            undone += 1
    return len(ch), undone


def count_stats(labs):
    counts = np.array([[int(np.sum(l == k)) for k in range(3)] for l in labs])
    hops = int(np.sum(np.any(np.diff(labs.astype(int), axis=0) != 0, axis=1)))
    frac = {}
    for c in counts:
        c = tuple(int(x) for x in c)
        frac[c] = frac.get(c, 0) + 1
    frac = {k: v / len(counts) for k, v in frac.items()}
    return counts, hops, frac


def report(res, table, t0=None):
    lines = ["# Remark D‴.16: dither and the mass constraint — the dithered averaged slow loop under the nearest-atom and the OT leak (`exp_dither_ot.py`)", "",
             f"The three-bump instance with the matching prior (three atoms at the exact representation, masses 1/3), σ₀ = 0.15σ_crit = {SIGMA0:.5f}, m = {M} particles started at the collapsed "
             f"configuration with counts {COUNTS0} (jitter 0.02), the Langevin dynamics dϑ_i = −G∇J_m ds + √(2τ₀)G^{{1/2}}dW_i in slow time (Euler–Maruyama, ds = {DT}; the leak target recomputed every "
             f"{PLAN_EVERY} steps), cluster masses by nearest atom recorded every {REC_EVERY * DT:g} slow units; horizons {', '.join(f'{S_END[t]:g}' for t in TAUS)} slow units at τ₀ = {', '.join(str(t) for t in TAUS)}; "
             f"seeds {', '.join(str(s) for s in SEEDS)}. The Gibbs prediction for the nearest-atom leak is the Laplace approximation over the count vectors, "
             "P(n) ∝ multinomial(n)·exp(−m𝒢*_{n/m}/τ₀)·Π_k det(σ₀I − R^G_k)^{−(n_k−1)/2}·det(H_n)^{−1/2} (count vectors whose replicated configuration is a saddle carry no Laplace mode and are omitted); "
             "'hop records' counts the records at which some particle's nearest atom changed.", ""]
    arr = {}
    for tau0 in TAUS:
        pred = gibbs_prediction(table, tau0); pred0 = gibbs_prediction(table, tau0, with_hessian=False)
        top = sorted(pred0.items(), key=lambda kv: -kv[1])[:5]; topH = sorted(pred.items(), key=lambda kv: -kv[1])[:3]
        lines += [f"## τ₀ = {tau0} (horizon {S_END[tau0]:g} slow units; the leak's part of the lowest ridge σ₀Δ_ν²/8 = {SIGMA0 * 6.263 ** 2 / 8:.3f} = {SIGMA0 * 6.263 ** 2 / 8 / tau0:.1f}τ₀)", "",
                  "The energy–multinomial law over the count vectors, multinomial(n)·exp(−m𝒢*_n/τ₀) (count vectors with n/m ∈ W_adm and the saddles alike; not an order of an expansion: the Laplace determinant factors are of the same order and are omitted), P(n): " + ", ".join(f"{n}: {p:.3f}" for n, p in top)
                  + "; with the Laplace determinant factors (saddle count vectors omitted): " + ", ".join(f"{n}: {p:.3f}" for n, p in topH) + ".", "",
                  "| leak | seed | records with a label change / records (undone at the next change) | time fraction by counts (top four) | mean |w − q|₁ | final counts | excursions off the balanced counts after the first balanced record: returns / exchanges; off-balance fraction |", "|---|---|---|---|---|---|---|"]
        for leak in ('nearest', 'ot'):
            for r in [r for r in res if r['tau0'] == tau0 and r['leak'] == leak]:
                counts, hops, frac = count_stats(r['labs'])
                nch, und = undone_changes(r['labs'])
                nret, nexc, offb, i0 = excursions(r['labs'])
                arr[(leak, tau0, r['seed'])] = (hops, len(r['labs']) - 1, float(np.mean(np.abs(counts / M - Q[None, :]).sum(axis=1))), frac.get((5, 5, 5), 0.0), nret, nexc, offb, und)
                top4 = sorted(frac.items(), key=lambda kv: -kv[1])[:4]
                lines.append(f"| {leak} | {r['seed']} | {hops}/{len(r['labs']) - 1} ({und}) | " + ", ".join(f"{k}: {v:.3f}" for k, v in top4) + f" | {arr[(leak, tau0, r['seed'])][2]:.3f} | {tuple(int(x) for x in counts[-1])} | "
                             + (f"{nret} / {nexc}; {offb:.2e}" if leak == 'ot' and np.isfinite(offb) else "—") + " |")
        lines.append("")
    # Arrhenius estimate from the nearest-atom hop rates at the two higher temperatures
    rates = {tau0: np.mean([arr[('nearest', tau0, s)][0] / S_END[tau0] for s in SEEDS]) for tau0 in TAUS}
    est = []
    for ta, tb in ((0.01, 0.02), (0.005, 0.01)):
        if rates[ta] > 0 and rates[tb] > 0:
            est.append((ta, tb, float(np.log(rates[tb] / rates[ta]) / (1 / ta - 1 / tb))))
    bal = {tau0: (np.mean([arr[('nearest', tau0, s)][3] for s in SEEDS]), np.mean([arr[('ot', tau0, s)][3] for s in SEEDS]), gibbs_prediction(table, tau0, with_hessian=False)[(5, 5, 5)]) for tau0 in TAUS}
    ret = {t: (sum(arr[('ot', t, s)][4] for s in SEEDS), sum(arr[('ot', t, s)][5] for s in SEEDS)) for t in TAUS}
    offb = {t: np.mean([arr[('ot', t, s)][6] for s in SEEDS]) for t in TAUS}
    eff = float(np.log(offb[0.02] / offb[0.01]) / (1 / 0.01 - 1 / 0.02)) if offb[0.01] > 0 and offb[0.02] > 0 else float('nan')
    ratios = [arr[('ot', t, s)][0] / arr[('nearest', t, s)][0] for t in (0.01, 0.02) for s in SEEDS if arr[('nearest', t, s)][0] > 0]
    lines += ["## Reading", "",
              "Nearest-atom leak: record-level label changes (every 10 slow units, about half undone at the next change — recrossings) at " + ", ".join(f"{rates[t]:.2e} per slow unit at τ₀ = {t}" for t in TAUS)
              + (" — Arrhenius barriers ΔV = " + ", ".join(f"{e:.3f} from τ₀ = {a}, {b}" for a, b, e in est) + " (the second pair rests on one event at τ₀ = 0.005)" if est else "") + f"; the leak's part of the lowest ridge is σ₀Δ_ν²/8 = {SIGMA0 * 6.263 ** 2 / 8:.3f}. "
              + "Time fraction at the balanced counts (5, 5, 5), nearest-atom leak against the OT leak against the energy–multinomial law: " + "; ".join(f"τ₀ = {t}: {bal[t][0]:.3f} / {bal[t][1]:.3f} / {bal[t][2]:.3f}" for t in TAUS)
              + " (two seeds; the balanced indicator's autocorrelation is long at τ₀ = 0.01, so these fractions carry standard errors of order 0.1 there). "
              + "OT leak: the plan's masses are q by construction; the counts by nearest atom leave the balanced counts in excursions that end with the crossing particle back in its cluster (returns) or, rarely, with two particles' assignments exchanged — "
              + "; ".join(f"τ₀ = {t}: {ret[t][0]} returns, {ret[t][1]} exchanges, off-balance fraction {offb[t]:.2e}" for t in TAUS)
              + f" — an effective barrier of {eff:.3f} between τ₀ = 0.01 and 0.02, and label changes as frequent as the nearest-atom leak's (ratios {min(ratios):.2f}–{max(ratios):.2f}): the crossing barrier is the nearest-atom ridge barrier, and what the mass constraint changes is that a crossing is reversible — "
              "the particle stays assigned to its own atom and is pulled back. Dither does on the Arrhenius time scale, and only in distribution, what the mass constraint does deterministically in ln 2/σ₀.", "",
              "Numbers are produced by `experiments/exp_dither_ot.py`; nothing is edited by hand."]
    with open(os.path.join(OUT, 'dither_ot.md'), 'w', encoding='utf-8') as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    if t0 is not None:
        print(f"[{time.time() - t0:.0f}s]")


def job(args):
    leak, tau0, seed = args
    return run(leak, tau0, seed)


def main():
    import multiprocessing as mp
    t0 = time.time()
    if '--from-raw' in sys.argv:
        z = np.load(os.path.join(OUT, 'dither_ot_raw.npz'), allow_pickle=True)
        report(list(z['res']), z['table'].item(), t0)
        return
    net = make_net(); loop = make_loop(net); exact = exact_representation(net)
    table = laplace_table(net, loop, exact)
    jobs = [(leak, tau0, seed) for tau0 in TAUS for leak in ('nearest', 'ot') for seed in SEEDS]
    with mp.Pool(2) as pool:
        res = pool.map(job, jobs, chunksize=1)
    np.savez(os.path.join(OUT, 'dither_ot_raw.npz'), res=np.array(res, dtype=object), table=np.array(table, dtype=object), allow_pickle=True)
    report(res, table, t0)


if __name__ == '__main__':
    main()
