"""
Is there a regime in which the dual order of entry matters? The sliding inverse scale space against a sliding matching pursuit (workstream C)
=============================================================================================================================================

The sliding inverse scale space (ISS, `order='dual'`: a candidate enters when its ACCUMULATED correlation reaches a common wall) and sliding matching pursuit
(`order='greedy'`: the candidate of the largest INSTANTANEOUS correlation enters) differ only in the order of entry; on SMLM frames and on astronomy patches the
controls (results/smlm/controls.md, results/astro/controls.md) found no gain from the dual order.  This script is a one-dimensional model of the whole pipeline
with every other ingredient identical for the two orders (dictionary, resolution limit, sliding, non-negative block fit, stopping rule), so that the order of entry is the only
difference, and scans the regimes in which the theory says the two could differ (a strong feature next to a weak one, close pairs, unequal widths, large noise, many features).

The one-dimensional model: samples x_i = i dx (dx = 0.5, NS = 96), a signal sum_k a_k exp(-(x - c_k)^2/(2 s_k^2)) with a_k > 0, white noise of sd `noise`
(known), the significance of a feature the matched-filter value a_k |phi|/noise.  Dictionary: the lattice of spacing H = 0.25 with the nominal width 1; resolution limit
`r_ex` (lattice points within r_ex of an active centre do not grow, as `imaging.sliding_path`), block fit by non-negative least squares then the trust-region fit with the
centres within `box` of their lattice starts (and, with `widths`, the widths in the bounds), the dormant kernel with w = 0 and a non-positive correlation leaves.
`path(order=...)` is the one function; the stopping rules are applied to the stored path (`bic`: the block of least chi^2 + n_par k ln N, known noise; `sig`: the leading entrants
of significance >= kappa; `disc`: the discrepancy rule).  Metrics: exact = K-hat = K and a one-to-one match within TOL of the truth; position rms of the exact ones.

Comparators (same data, same rules): the greedy order; each order followed by a backward step (`prune_block`: remove the unit whose removal raises chi^2 least while
that is below the penalty); sliding Frank-Wolfe for the BLASSO (`baselines.sfw`, a fixed ladder of eight penalties, the penalty chosen by the same BIC / penalty-16 rule).
Regimes (see results/synthetic/dualregime_plan.md): P close pairs, Q crowded windows, W four features of unequal widths with the widths free, E separation at or beyond the
theory's entry threshold (equivalence), part T (noise-free entry sequences), the two-dimensional arm (fresh scenes of `exp_order_pairs`), and a POST HOC reading of the stored
SMLM test detections by the distance to the nearest neighbour.  The development scan (`--scan`) covers grids A (pairs), B (random clusters), C (unequal widths), D (crowds).

    python experiments/exp_dualregime.py --scan  [--grids C,D,A,B] [--workers 1]   development scan (seeds from DEV_SEED0), writes dualregime_dev_raw.npz
    python experiments/exp_dualregime.py --run   [--parts T,P,E,Q,W,2d]            the pre-registered held-out run (seeds from TEST_SEED0), writes dualregime_raw.npz
    python experiments/exp_dualregime.py --smlm                                     the POST HOC SMLM analysis (reads results/smlm/controls_raw.npz), into dualregime_raw.npz
    python experiments/exp_dualregime.py --from-raw                                 regenerates results/synthetic/dualregime.md from the raws
    (DR_OUT redirects the raws and the md, DR_WORKERS the processes of a grid or of the two-dimensional arm)
"""
import argparse
import os
import sys

os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
from scipy.optimize import least_squares, linear_sum_assignment, nnls

HERE = os.path.abspath(os.path.dirname(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, REPO)

OUT = os.environ.get('DR_OUT', os.path.join(REPO, 'results', 'synthetic'))
RAW, MD = os.path.join(OUT, 'dualregime_raw.npz'), os.path.join(OUT, 'dualregime.md')

DX, NS, H, MARGIN = 0.5, 96, 0.25, 2.0
XS = DX * np.arange(NS)
LEN = DX * (NS - 1)
LAT = np.arange(-MARGIN, LEN + MARGIN + 1e-9, H)
TOL = 0.7                      # match radius (nominal widths), as exp_order_pairs
SAME_TOL = 0.05                # two blocks are the same set when the sorted centres agree within this
NORM2 = float(np.sum(np.exp(-(XS - XS[NS // 2]) ** 2)))     # |phi|^2 of the nominal kernel (interior)


def gauss(c, s=1.0):
    """Kernel columns (NS, k): exp(-(x - c_k)^2 / (2 s_k^2))."""
    c = np.atleast_1d(np.asarray(c, float))
    s = np.broadcast_to(np.asarray(s, float), c.shape)
    return np.exp(-0.5 * ((XS[:, None] - c[None, :]) / s[None, :]) ** 2)


_PHI_LAT = gauss(LAT)


def block_fit(y, c, s, w, c0, widths, box=1.5, max_nfev=100, tol=1e-10):
    """The block's fit from (c, s, w): non-negative least squares at the centres in place, then the trust-region fit of all (w, c[, s]) with the centres in c0 +- box
    and, if widths = (lo, hi), the widths inside.  Returns (c, s, w)."""
    k = len(w)
    if k == 0:
        return c, s, w
    w1, _ = nnls(gauss(c, s), y, maxiter=200 * k)
    lo = np.r_[np.zeros(k), c0 - box] if widths is None else np.r_[np.zeros(k), c0 - box, np.full(k, widths[0])]
    hi = np.r_[np.full(k, np.inf), c0 + box] if widths is None else np.r_[np.full(k, np.inf), c0 + box, np.full(k, widths[1])]
    z0 = np.clip(np.r_[w1, c] if widths is None else np.r_[w1, c, s], lo, hi)

    def parts(z):
        return z[:k], z[k:2 * k], (np.ones(k) if widths is None else z[2 * k:])

    def res(z):
        ww, cc, ss = parts(z)
        return y - gauss(cc, ss) @ ww

    def jac(z):
        ww, cc, ss = parts(z)
        G = gauss(cc, ss)
        d = XS[:, None] - cc[None, :]
        cols = [G, G * d / ss ** 2 * ww]
        if widths is not None:
            cols.append(G * d ** 2 / ss ** 3 * ww)
        return -np.concatenate(cols, axis=1)
    sol = least_squares(res, z0, jac=jac, bounds=(lo, hi), method='trf', x_scale='jac', xtol=tol, ftol=tol, gtol=tol, max_nfev=max_nfev)
    ww, cc, ss = parts(sol.x)
    return cc.copy(), ss.copy(), ww.copy()


def path(y, noise, order='dual', r_ex=1.5, widths=None, box=1.5, max_entries=10, rho=1.0, vtol=1e-12, tie_tol=1e-9):
    """The sliding path of the data y (white noise of sd `noise`) with the entry order `order` ('dual': the accumulated correlation reaches the wall rho; 'greedy': the
    largest instantaneous correlation enters, tau = the step count).  Returns dict(blocks = [dict(c, s, w, rss)] (blocks[0] the empty fit), taus, snr (each entrant's
    matched-filter significance v / (noise |phi|)), ids (lattice indices of the active units), entered (the lattice indices of every entry, in the order of entry))."""
    D = len(LAT)
    U = np.zeros(D)
    taken = np.zeros(D, bool)
    c, s, w, c0 = np.zeros(0), np.zeros(0), np.zeros(0), np.zeros(0)
    ids, taus, snrs, tau, entered = [], [], [], 0.0, []
    blocks = [dict(c=c.copy(), s=s.copy(), w=w.copy(), rss=float(y @ y))]
    model = np.zeros(NS)
    vscale = None
    for _ in range(11 * max_entries):
        V = _PHI_LAT.T @ (y - model)
        if vscale is None:
            vscale = float(np.abs(V).max())
            if vscale <= 0.0:
                break
        blocked = taken.copy()
        if r_ex > 0.0 and len(c):
            blocked |= (np.abs(LAT[:, None] - c[None, :]) < r_ex).any(1)
        v = np.where(blocked, 0.0, V)
        grow = v > vtol * vscale
        if not grow.any():
            break
        if order == 'greedy':
            j = int(np.argmax(np.where(grow, v, -np.inf)))
            tau += 1.0
        else:
            dt = np.full(D, np.inf)
            dt[grow] = (rho - U[grow]) / v[grow]
            dmin = float(dt.min())
            j = int(np.argmax(dt <= dmin * (1.0 + tie_tol) + tie_tol * 1e-300))
            U[~blocked] += dmin * v[~blocked]
            tau += dmin
        U[j] = rho
        taken[j] = True
        ids.append(j)
        entered.append(j)
        snrs.append(float(v[j] / (noise * np.sqrt(NORM2))))
        taus.append(tau)
        c0 = np.r_[c0, LAT[j]]
        c, s, w = np.r_[c, LAT[j]], np.r_[s, 1.0], np.r_[w, 0.0]
        c, s, w = block_fit(y, c, s, w, c0, widths, box)
        model = gauss(c, s) @ w
        g = gauss(c, s).T @ (y - model)
        gone = (w <= 1e-12 * max(w.max(), 1e-300)) & (g <= vtol * vscale)
        for i in np.flatnonzero(gone)[::-1]:
            taken[ids[i]] = False
            U[ids[i]] = rho
            del ids[i]
            c, s, w, c0 = np.delete(c, i), np.delete(s, i), np.delete(w, i), np.delete(c0, i)
        if gone.any():
            model = gauss(c, s) @ w if len(w) else np.zeros(NS)
        r = y - model
        blocks.append(dict(c=c.copy(), s=s.copy(), w=w.copy(), rss=float(r @ r)))
        if len(taus) >= max_entries:
            break
    return dict(blocks=blocks, taus=np.array(taus), snr=np.array(snrs), ids=ids, entered=entered)


# ---- the stopping rules, applied to the stored path --------------------------------------------------------------------------------------------------------------
def select_pen(blocks, noise, pen):
    """Index of the block of least rss/noise^2 + pen k (the first minimum): a unit is kept if it lowers chi^2 by more than pen.  BIC is pen = n_par ln N."""
    cost = [b['rss'] / noise ** 2 + pen * len(b['w']) for b in blocks]
    return int(np.argmin(cost))


def select_time(blocks, taus, kappa, noise):
    """The stopping time of the inverse scale space: the last block with at most K-hat entrants, K-hat = the entrants with tau <= 1/(kappa noise |phi|) (a feature of
    significance kappa correlates at kappa noise |phi| and enters at tau = 1 / that)."""
    K = int(np.sum(np.asarray(taus) <= 1.0 / (kappa * noise * np.sqrt(NORM2))))
    return max(i for i, b in enumerate(blocks) if len(b['w']) <= K) if K else 0


def select_sig(blocks, snr, kappa):
    """Index of the last block with at most K-hat entrants, K-hat = the number of leading entrants of significance >= kappa."""
    bad = np.flatnonzero(np.asarray(snr) < kappa)
    K = int(bad[0]) if bad.size else len(snr)
    return max(i for i, b in enumerate(blocks) if len(b['w']) <= K) if K else 0


def select_disc(blocks, noise, z=2.0):
    """The first block whose mean square residual is within the noise: rss/N <= noise^2 (1 + z sqrt(2/N))."""
    thr = noise ** 2 * (1.0 + z * np.sqrt(2.0 / NS))
    for i, b in enumerate(blocks):
        if b['rss'] / NS <= thr:
            return i
    return len(blocks) - 1


def match(c_est, c_true, tol=TOL):
    """One-to-one match within tol: (number matched, position rms over the matched pairs)."""
    if len(c_est) == 0 or len(c_true) == 0:
        return 0, float('nan')
    d = np.abs(np.asarray(c_est)[:, None] - np.asarray(c_true)[None, :])
    rr, cc = linear_sum_assignment(np.where(d <= tol, d, 1e6))
    ok = d[rr, cc] <= tol
    if not ok.any():
        return 0, float('nan')
    return int(ok.sum()), float(np.sqrt(np.mean(d[rr, cc][ok] ** 2)))


def score(block, c_true):
    """(exact, position rms) of a block: exact = as many units as truth and all matched."""
    n, rms = match(block['c'], c_true)
    ex = int(len(block['w']) == len(c_true) and n == len(c_true))
    return ex, (rms if ex else float('nan'))


# ---- the scenes ----------------------------------------------------------------------------------------------------------------------------------------------
def feat_norm(s):
    """|phi| of a kernel of width s (interior)."""
    return np.sqrt(np.sum(np.exp(-(XS - XS[NS // 2]) ** 2 / s ** 2)))


def scene(rng, fam):
    """A scene from the family dict `fam`: returns (c, a, s, noise).  kind 'pair': two features at separation fam['d'] (nominal widths), amplitude ratio fam['ratio'] (the
    stronger first), the weaker of significance fam['snr'], widths 1 and fam['wr'] (the weaker's, default 1), a random position and orientation.  kind 'cluster' (or 'crowd': the centres in a window of width K fam['gap']): fam['K'] features,
    separations at least fam['dmin'] (rejection), amplitudes log-uniform over a factor fam['ratio'], widths log-uniform in [1/fam['wr'], fam['wr']] (default 1), the weakest of
    significance fam['snr']; the centres in [6, LEN - 6]."""
    lo, hi = 6.0, LEN - 6.0
    if fam['kind'] == 'pair':
        m = rng.uniform(lo + 2.0, hi - 2.0)
        sgn = rng.choice([-1.0, 1.0])
        c = np.array([m - 0.5 * sgn * fam['d'], m + 0.5 * sgn * fam['d']])
        a = np.array([fam['ratio'], 1.0])
        s = np.array([1.0, fam.get('wr', 1.0)])
    else:
        K = fam['K']
        if fam['kind'] == 'crowd':                         # the K centres inside a window of width K * fam['gap'] at a random place
            lo = rng.uniform(lo, hi - K * fam['gap'])
            hi = lo + K * fam['gap']
        for _ in range(100000):
            c = np.sort(rng.uniform(lo, hi, K))
            if K == 1 or np.diff(c).min() >= fam['dmin']:
                break
        else:
            raise RuntimeError('no placement')
        rng.shuffle(c)
        a = np.exp(rng.uniform(0.0, np.log(fam['ratio']), K)) if fam['ratio'] > 1 else np.ones(K)
        a = a / a.min()
        wr = fam.get('wr', 1.0)
        s = np.exp(rng.uniform(-np.log(wr), np.log(wr), K)) if wr > 1 else np.ones(K)
    sig = a * np.array([feat_norm(x) for x in s])           # significance units
    noise = sig.min() / fam['snr']
    return c, a, s, noise


def prune_block(y, b, noise, pen, widths=None):
    """Backward elimination of one block: while removing some unit (the others refitted from where they are) raises chi^2 by less than pen, remove the unit whose removal
    raises it least.  Returns the pruned block (the same dict form)."""
    c, w, sg = b['c'].copy(), b['w'].copy(), b['s'].copy()
    rss = b['rss']
    while len(w) > 0:
        best = None
        for i in range(len(w)):
            keep = np.arange(len(w)) != i
            cc, ss, ww = block_fit(y, c[keep], sg[keep], w[keep], c[keep], widths)
            r = y - gauss(cc, ss) @ ww if len(ww) else y
            q = float(r @ r)
            if best is None or q < best[0]:
                best = (q, cc, ss, ww)
        if (best[0] - rss) / noise ** 2 < pen:
            rss, c, sg, w = best
        else:
            break
    return dict(c=c, s=sg, w=w, rss=float(rss))


RULES = ('bic', 'pen', 'disc', 'nat', 'bicbw', 'penbw', 'oracleK', 'path')
PEN_BIC, PEN_HI = 2.0 * np.log(NS), 16.0
KAPPA_NAT = 4.5
PRUNE_EXTRA, EXTRA_ENTRIES = 1, 3     # backward elimination of the blocks of at most K + 1 units; a path runs to K + 3 entries
METHODS = ('dual', 'greedy', 'sfw')
LAM_REL = tuple(np.geomspace(0.5, 0.003, 8))


def sfw_blocks(y, noise, min_sig=3.0):
    """Sliding Frank-Wolfe for the BLASSO (baselines.sfw) on the ladder of penalties lam = f max|<y phi>|, f = LAM_REL (fixed a priori); the spikes whose own significance
    w |phi| / noise is below min_sig are dropped (the weak spikes of a small penalty fit noise), and the ladder's solutions are returned as a list of blocks (from the empty fit,
    chi^2 recomputed after the drop)."""
    from src.core import baselines as Bs
    from src.core import rich_transit_nd as RN
    Kd = RN.BumpsND(XS[:, None], 1.0)
    grid = LAT[:, None]
    G = float(np.abs(Kd.corr(grid, y)[0]).max())
    out = [dict(c=np.zeros(0), s=np.zeros(0), w=np.zeros(0), rss=float(y @ y))]
    for f in LAM_REL:
        r = Bs.sfw(Kd, y, grid, f * G, max_units=10)
        w, c = np.asarray(r['w']), np.asarray(r['c']).ravel()
        keep = (w * np.sqrt(NORM2) / noise) >= min_sig
        w, c = w[keep], c[keep]
        res = y - gauss(c) @ w if len(w) else y
        out.append(dict(c=c.copy(), s=np.ones(len(w)), w=w.copy(), rss=float(res @ res)))
    return out


def apply_rules(y, B, c_true, noise, widths, nat=None, bw=True):
    """The rules' (exact, position rms) for one path or ladder B (list of blocks): returns array (len(RULES), 2)."""
    out = np.full((len(RULES), 2), np.nan)
    npar_pen = PEN_BIC * (2 if widths is None else 3) / 2.0
    pick = {'bic': select_pen(B, noise, npar_pen), 'pen': select_pen(B, noise, PEN_HI), 'disc': select_disc(B, noise)}
    for k, i in pick.items():
        out[RULES.index(k)] = score(B[i], c_true)
    if nat is not None:
        out[RULES.index('nat')] = score(B[nat], c_true)
    for k, pen in (('bicbw', npar_pen), ('penbw', PEN_HI)) if bw else ():
        P = [prune_block(y, b, noise, pen, widths) if len(b['w']) <= len(c_true) + PRUNE_EXTRA else b for b in B]
        out[RULES.index(k)] = score(P[select_pen(P, noise, pen)], c_true)
    kk = [i for i, b in enumerate(B) if len(b['w']) == len(c_true)]
    out[RULES.index('oracleK')] = score(B[kk[0]], c_true) if kk else (0.0, np.nan)
    ex = [q for q in (score(b, c_true) for b in B) if q[0] == 1]
    out[RULES.index('path')] = (1.0, min(q[1] for q in ex)) if ex else (0.0, np.nan)
    return out


def same_block(b1, b2, tol=SAME_TOL):
    """The two blocks have the same number of units and the sorted centres agree within tol."""
    return int(len(b1['w']) == len(b2['w']) and (len(b1['w']) == 0 or float(np.abs(np.sort(b1['c']) - np.sort(b2['c'])).max()) <= tol))


def trial_full(fam, seed, r_ex=1.5, widths=None, methods=METHODS):
    """One scene, every method on the same data.  Returns (out, agree): out = array (len(methods), len(RULES), 2) = (exact, position rms); agree = (the dual's and
    the greedy path's blocks selected by 'bic' are the same set, ... by 'pen').  'nat' is each method's own natural rule: the stopping time (kappa 4.5) for the dual, the
    entrant's significance (4.5) for matching pursuit (sfw: none).  'oracleK': the first block with exactly K units; 'path': some block is exact (the position rms of the best).
    'bicbw', 'penbw': the rule after backward elimination of every block."""
    rng = np.random.default_rng(seed)
    c, a, s, noise = scene(rng, fam)
    y = gauss(c, s) @ a + noise * rng.standard_normal(NS)
    out = np.full((len(methods), len(RULES), 2), np.nan)
    npar_pen = PEN_BIC * (2 if widths is None else 3) / 2.0
    sel = {}
    for m, name in enumerate(methods):
        if name == 'sfw':
            out[m] = apply_rules(y, sfw_blocks(y, noise), c, noise, widths, bw=False)
            continue
        R = path(y, noise, name, r_ex=r_ex, widths=widths, max_entries=len(c) + EXTRA_ENTRIES)
        nat = select_time(R['blocks'], R['taus'], KAPPA_NAT, noise) if name == 'dual' else select_sig(R['blocks'], R['snr'], KAPPA_NAT)
        out[m] = apply_rules(y, R['blocks'], c, noise, widths, nat=nat)
        sel[name] = (R['blocks'][select_pen(R['blocks'], noise, npar_pen)], R['blocks'][select_pen(R['blocks'], noise, PEN_HI)])
    agree = np.array([same_block(sel['dual'][0], sel['greedy'][0]), same_block(sel['dual'][1], sel['greedy'][1])]) if 'dual' in sel and 'greedy' in sel else np.zeros(2, int)
    return out, agree


def trial(fam, seed, r_ex=1.5, widths=None, methods=METHODS):
    """`trial_full` without the agreement flags."""
    return trial_full(fam, seed, r_ex, widths, methods)[0]


# ---- the cells ----------------------------------------------------------------------------------------------------------------------------------------------
def run_cell(args):
    """(fam, seed0, n_trials, r_ex, widths) -> array (n_trials, len(METHODS), len(RULES), 2); trial i uses seed seed0 + i."""
    fam, seed0, n, r_ex, widths = args
    return np.array([trial(fam, seed0 + i, r_ex=r_ex, widths=widths) for i in range(n)])


def diff_stats(E, a, b, rule):
    """Paired difference in exact-recovery rate, methods a minus b under `rule` (points): (mean, standard error) over the trials of E (n, methods, rules, 2)."""
    d = (E[:, METHODS.index(a), RULES.index(rule), 0] - E[:, METHODS.index(b), RULES.index(rule), 0]) * 100.0
    return float(d.mean()), float(d.std(ddof=1) / np.sqrt(len(d)))


# ---- part T: noise-free, the entry sequence ------------------------------------------------------------------------------------------------------------------
T_REX, T_DS, T_RATIOS, T_N, T_SEED = (0.0, 1.5), (1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0), (1.0, 2.0, 5.0, 20.0), 12, 31


def part_t():
    """Noise-free pairs (y = a1 phi(c1) + a2 phi(c1 + d), amplitude ratio a1/a2, the first centre at 20 + u, u uniform in [0, 1), so that the pair sits at a random phase of the
    lattice), both orders run for two entries.  Returns array (len(T_REX), len(T_DS), len(T_RATIOS), T_N, 4) of 0/1: the entry sequences are the same; the blocks after two
    entries are the same set; the dual's block is exact; the greedy block is exact."""
    out = np.zeros((len(T_REX), len(T_DS), len(T_RATIOS), T_N, 4))
    for i, r_ex in enumerate(T_REX):
        for j, d in enumerate(T_DS):
            for k, ratio in enumerate(T_RATIOS):
                rng = np.random.default_rng(T_SEED + 1000 * i + 10 * j + k)
                for n in range(T_N):
                    c = 20.0 + rng.uniform(0.0, 1.0) + np.array([0.0, d])
                    y = gauss(c) @ np.array([ratio, 1.0])
                    R = {o: path(y, 1e-6, o, r_ex=r_ex, max_entries=2) for o in ('dual', 'greedy')}
                    out[i, j, k, n] = (R['dual']['entered'] == R['greedy']['entered'], same_block(R['dual']['blocks'][-1], R['greedy']['blocks'][-1]),
                                       score(R['dual']['blocks'][-1], c)[0], score(R['greedy']['blocks'][-1], c)[0])
    return out


# ---- the grids: development scan and held-out run ------------------------------------------------------------------------------------------------------------
DEV_SEED0, TEST_SEED0 = 1000000, 7000000
WIDTHS = (0.4, 3.0)
DEV_GRIDS = {}
DEV_GRIDS['A'] = (60, [(dict(kind='pair', d=d, ratio=r, snr=q), None) for d in (1.2, 1.6, 2.0, 2.5, 3.5, 5.0) for r in (1, 3, 10, 30) for q in (5, 8, 12, 20)])
DEV_GRIDS['B'] = (40, [(dict(kind='cluster', K=K, dmin=dm, ratio=r, snr=q), None) for K in (4, 7) for dm in (1.5, 2.2, 3.2) for r in (1, 5, 25) for q in (5, 8, 14)])
DEV_GRIDS['C'] = (40, [(dict(kind='cluster', K=K, dmin=dm, ratio=r, snr=q, wr=w), WIDTHS) for w in (1.5, 2.0) for K in (2, 4) for dm in (1.5, 2.5) for r in (1, 5) for q in (6, 12)])
DEV_GRIDS['D'] = (50, [(dict(kind='crowd', K=K, gap=g, dmin=1.0, ratio=r, snr=q), None) for K in (3, 5) for g in (1.6, 2.2, 3.0) for r in (1, 6) for q in (8, 14)])
TEST_GRIDS = {}
TEST_GRIDS['P'] = (250, [(dict(kind='pair', d=d, ratio=r, snr=q), None) for d in (1.2, 1.6) for r in (1, 3) for q in (8, 12)])
TEST_GRIDS['E'] = (100, [(dict(kind='pair', d=d, ratio=r, snr=q), None) for d in (3.5, 5.0) for r in (1, 10) for q in (8, 12)])
TEST_GRIDS['Q'] = (100, list(DEV_GRIDS['D'][1]))          # the whole crowd grid D, fresh seeds
TEST_GRIDS['W'] = (120, [c for c in DEV_GRIDS['C'][1] if c[0]['K'] == 4])     # unequal widths, free in WIDTHS: the K = 4 cells of grid C, fresh seeds


def grid_seeds(grids, seed0, step):
    """The first seed of every cell of every grid, in the order of the dict: seed0 + step * (running cell index)."""
    out, k = {}, 0
    for g, (n, cells) in grids.items():
        out[g] = []
        for _ in cells:
            out[g].append(seed0 + step * k)
            k += 1
    return out


def run_grid(grids, seed0, step, which, workers=1):
    """Run the grids named in `which`; returns the dict of raw arrays: '<g>|<j>|fam' (json), '<g>|<j>|widths' (json), '<g>|<j>|E' (n, methods, rules, 2), '<g>|<j>|A' (n, 2)."""
    import json
    seeds = grid_seeds(grids, seed0, step)
    raw = {}
    for g in which:
        n, cells = grids[g]
        jobs = [(fam, seeds[g][j] + i, 1.5, wd) for j, (fam, wd) in enumerate(cells) for i in range(n)]
        if workers > 1:
            from multiprocessing import Pool
            with Pool(workers) as pool:
                res = pool.map(_job, jobs, chunksize=4)
        else:
            res = [_job(a) for a in jobs]
        for j, (fam, wd) in enumerate(cells):
            r = res[j * n:(j + 1) * n]
            raw['%s|%d|fam' % (g, j)] = np.array(json.dumps(fam, sort_keys=True))
            raw['%s|%d|widths' % (g, j)] = np.array(json.dumps(wd))
            raw['%s|%d|E' % (g, j)] = np.array([q[0] for q in r])
            raw['%s|%d|A' % (g, j)] = np.array([q[1] for q in r])
        print('grid %s done' % g, flush=True)
    return raw


def _job(a):
    fam, seed, r_ex, wd = a
    return trial_full(fam, seed, r_ex=r_ex, widths=wd)


# ---- statistics ----------------------------------------------------------------------------------------------------------------------------------------------
BOOT, BOOT_SEED = 2000, 0


def paired_boot(d, boot=BOOT):
    """(mean, lo, hi) of the per-scene differences d (a 1-D array) with the 95 % percentile interval of a bootstrap over scenes (BOOT resamples, seed BOOT_SEED)."""
    d = np.asarray(d, float)
    d = d[np.isfinite(d)]
    if len(d) == 0:
        return float('nan'), float('nan'), float('nan')
    idx = np.random.default_rng(BOOT_SEED).integers(0, len(d), size=(boot, len(d)))
    lo, hi = np.percentile(d[idx].mean(1), [2.5, 97.5])
    return float(d.mean()), float(lo), float(hi)


def verdict(lo, hi):
    """win: the interval excludes 0 above; loss: below; tie: it contains 0."""
    return 'win' if lo > 0 else 'loss' if hi < 0 else 'tie'


def clopper_pearson(k, n, alpha=0.05):
    """Two-sided Clopper-Pearson interval of a proportion k / n."""
    from scipy.stats import beta
    lo = 0.0 if k == 0 else float(beta.ppf(alpha / 2, k, n - k + 1))
    hi = 1.0 if k == n else float(beta.ppf(1 - alpha / 2, k + 1, n - k))
    return lo, hi


def cell_keys(raw, g):
    """The cell indices of grid g present in the raw, in order."""
    j = 0
    while '%s|%d|E' % (g, j) in raw:
        j += 1
    return list(range(j))


def exact_of(E, method, rule):
    """Per-scene exact flags (n,) of a method under a rule from a cell array E (n, methods, rules, 2)."""
    return E[:, METHODS.index(method), RULES.index(rule), 0]


def pooled(raw, g, method, rule, js=None):
    """The exact flags of a method and rule concatenated over the cells of grid g (or the cells js), in scene order."""
    js = cell_keys(raw, g) if js is None else js
    return np.concatenate([exact_of(raw['%s|%d|E' % (g, j)], method, rule) for j in js])


def pooled_rms(raw, g, method, rule):
    js = cell_keys(raw, g)
    return np.concatenate([raw['%s|%d|E' % (g, j)][:, METHODS.index(method), RULES.index(rule), 1] for j in js])


def fam_of(raw, g, j):
    import json
    return json.loads(str(raw['%s|%d|fam' % (g, j)]))


def fam_label(fam):
    return ' '.join('%s=%g' % (k, v) for k, v in fam.items() if k != 'kind') if fam['kind'] != 'pair' else 'd=%g ratio=%g snr=%g' % (fam['d'], fam['ratio'], fam['snr'])


# ---- the two-dimensional arm: fresh scenes of exp_order_pairs ------------------------------------------------------------------------------------------------
TWOD_SEED0, TWOD_SCENES, TWOD_STOP = 9000000, 100, ('bic', 4.0)
TWOD_SEPS, TWOD_PHOTONS = (1.2, 1.6, 2.5, 3.0), (0, 1)           # separations (widths); photon levels = indices into exp_order_pairs.PHOTONS (2000, 6000)
TWOD_CLOSE, TWOD_MID = (1.2, 1.6), (2.5, 3.0)                    # the two groups of the plan: expected dual ahead / dual behind


def part_2d(workers=1):
    """Scenes of exp_order_pairs (the same generator, EM-CCD frames, equal fluxes) with fresh seeds TWOD_SEED0 + 1000 * cell + s; the classes under the dual and the greedy
    order (0 mislocated, 1 resolved, 2 merged, 3 split) at the stop TWOD_STOP.  Returns array (len(TWOD_SEPS), len(TWOD_PHOTONS), TWOD_SCENES, 2)."""
    sys.path.insert(0, HERE)
    import exp_order_pairs as OP
    jobs = [(d, 1.0, ip, TWOD_SEED0 + 1000 * (len(TWOD_PHOTONS) * i + k) + s, TWOD_STOP, None) for i, d in enumerate(TWOD_SEPS) for k, ip in enumerate(TWOD_PHOTONS) for s in range(TWOD_SCENES)]
    if workers > 1:
        from multiprocessing import Pool
        with Pool(workers, initializer=OP._init) as pool:
            res = pool.map(OP.trial, jobs, chunksize=4)
    else:
        OP._init()
        res = [OP.trial(a) for a in jobs]
    return np.array([[o[0] for o in r] for r in res]).reshape(len(TWOD_SEPS), len(TWOD_PHOTONS), TWOD_SCENES, 2)


# ---- the SMLM frames, POST HOC: pair members in the stored test detections -----------------------------------------------------------------------------------
SMLM_BANDS_NM = ((0.0, 200.0), (200.0, 300.0), (300.0, 450.0), (450.0, 1e9))      # nearest-neighbour distance of a ground-truth emitter in its frame (the PSF sd is about 150 nm)
SMLM_TOL_NM, SMLM_BLOCK = 150.0, 10


def part_smlm():
    """POST HOC.  From results/smlm/controls_raw.npz (read only): for each set the test frames' ground-truth emitters, each with its nearest-neighbour distance in its frame, and whether
    it is matched (one to one, Hungarian, SMLM_TOL_NM) by a detection of the dual ISS (`ISS-E`) and of the greedy order at the dual's operating point (`GR-E`).  Returns the dict of
    raw arrays '<set>|nn', '<set>|frame', '<set>|hit_iss', '<set>|hit_gr'."""
    sys.path.insert(0, HERE)
    import exp_smlm_controls as C
    import smlm_common as S
    raw = C.load_raw()
    out = {}
    for name in C.sets_in(raw):
        a, b = (int(v) for v in raw['ref|%s|frames' % name][2:4])
        gt = raw['ref|%s|gt_test' % name]
        dets = {m: C.dets_of(raw, name, m) for m in ('ISS-E', 'GR-E')}
        if any(v is None for v in dets.values()):
            continue
        nn, fr, hi, hg = [], [], [], []
        for k, f in enumerate(range(a, b)):
            g = gt[gt['frame'] == f]
            xy = np.c_[g['x_nm'], g['y_nm']]
            if len(xy) == 0:
                continue
            D = np.linalg.norm(xy[:, None] - xy[None], axis=2) + 1e9 * np.eye(len(xy))
            hits = []
            for m in ('ISS-E', 'GR-E'):
                det = np.asarray(dets[m][k], float)
                det = det[:, :2] if len(det) else np.zeros((0, 2))
                M = S.match(det, xy, tol_nm=SMLM_TOL_NM)
                h = np.zeros(len(xy), bool)
                h[M.pairs[:, 1]] = True
                hits.append(h)
            nn.append(D.min(1))
            fr.append(np.full(len(xy), f))
            hi.append(hits[0])
            hg.append(hits[1])
        out[name + '|nn'], out[name + '|frame'] = np.concatenate(nn), np.concatenate(fr)
        out[name + '|hit_iss'], out[name + '|hit_gr'] = np.concatenate(hi), np.concatenate(hg)
    return out


# ---- the report ----------------------------------------------------------------------------------------------------------------------------------------------
DEV_RAW = os.path.join(OUT, 'dualregime_dev_raw.npz')
COMPARISONS = (              # (id, label, (method, rule) of the first arm, (method, rule) of the second arm); the difference is first minus second, in points
    ('1', 'dual + backward step vs greedy + backward step (order only)', ('dual', 'bicbw'), ('greedy', 'bicbw')),
    ('2', 'dual vs greedy (as the pipelines stand)', ('dual', 'bic'), ('greedy', 'bic')),
    ('3', 'dual vs greedy + backward step', ('dual', 'bic'), ('greedy', 'bicbw')),
    ('4', 'dual vs sliding Frank-Wolfe (BLASSO)', ('dual', 'bic'), ('sfw', 'bic')),
    ('5', 'dual + backward step vs sliding Frank-Wolfe', ('dual', 'bicbw'), ('sfw', 'bic')),
    ('6', 'dual vs greedy, penalty 16 (sensitivity)', ('dual', 'pen'), ('greedy', 'pen')),
    ('7', 'dual + backward step vs greedy + backward step, penalty 16 (sensitivity)', ('dual', 'penbw'), ('greedy', 'penbw')),
    ('8', 'dual vs sliding Frank-Wolfe, penalty 16 (sensitivity)', ('dual', 'pen'), ('sfw', 'pen')),
    ('9', 'dual vs greedy, first block with K units (order only, no stopping rule)', ('dual', 'oracleK'), ('greedy', 'oracleK')),
    ('10', 'dual vs greedy, own natural rule (stopping time vs entrant significance)', ('dual', 'nat'), ('greedy', 'nat')),
)


BW_IDS = ('1', '3', '5', '7')            # the comparisons with a backward-step arm (the backward step is capped at K_true + PRUNE_EXTRA units: K leaks into these arms)


def comparison_status(g, cid):
    """How a pooled comparison may be quoted: the order-only ones, the oracle-K ones, the ones against the sliding Frank-Wolfe solver (a different algorithm), and those that are not an order comparison."""
    st = []
    if cid == '1':
        st.append('order effect, everything else equal')
    if cid == '9':
        st.append('oracle K (order only, uses the true number of units)')
    if cid == '10':
        st.append('not an order comparison (different stopping rules), not to be quoted')
    if cid in ('4', '5', '8'):
        st.append('different algorithm (sliding Frank-Wolfe), not an order comparison' + ('; misspecified here (fixed width), not to be quoted' if g == 'W' else ''))
    if g == 'W' and cid in ('6', '7'):
        st.append('penalty-16 sensitivity in the regime W rows W4-W8, not to be quoted')
    if cid in BW_IDS:
        st.append('backward step capped at K_true + %d units' % PRUNE_EXTRA)
    return '; '.join(st)


def md_table(head, rows):
    return ['| ' + ' | '.join(head) + ' |', '|' + '---|' * len(head)] + ['| ' + ' | '.join(str(c) for c in r) + ' |' for r in rows]


def pct(x):
    return '%.0f' % (100.0 * x)


def dev_lines(dev):
    """The development-scan tables: one row per cell, the exact-recovery rate under BIC of the dual / greedy / sFW and the paired differences (points, standard error)."""
    L = []
    flags = []
    for g in DEV_GRIDS:
        if '%s|0|E' % g not in dev:
            continue
        n, cells = DEV_GRIDS[g]
        rows = []
        for j in cell_keys(dev, g):
            E = dev['%s|%d|E' % (g, j)]
            fam = fam_of(dev, g, j)
            ds = [diff_stats(E, a, b, r) for a, b, r in (('dual', 'greedy', 'bic'), ('dual', 'greedy', 'pen'), ('dual', 'greedy', 'oracleK'), ('dual', 'greedy', 'bicbw'), ('dual', 'sfw', 'bic'))]
            flags.append((g, j, fam, ds))
            rows.append([fam_label(fam)] + [pct(np.mean(exact_of(E, m, 'bic'))) for m in METHODS] + ['%+.0f (%.0f)' % d for d in ds])
        L += ['### Grid %s (%d scenes per cell, seeds from %d)' % (g, n, grid_seeds(DEV_GRIDS, DEV_SEED0, 1000)[g][0]), '']
        L += md_table(['cell', 'dual', 'greedy', 'sFW', 'dual - greedy, BIC', 'pen 16', 'first block with K units', 'both + backward step', 'dual - sFW, BIC'], rows) + ['']
    return L, flags


def smlm_lines(sm):
    """The POST HOC table of the SMLM sets: recall of ground-truth emitters by the distance to their nearest neighbour in the frame, the dual ISS and the greedy order at the dual's E-operating point, paired difference with
    the 95 % interval of a block bootstrap over the test frames (blocks of SMLM_BLOCK consecutive frames, BOOT resamples, seed BOOT_SEED)."""
    sets = sorted(set(k.split('|')[0] for k in sm))
    order = ['MT0.N1.LD', 'MT0.N1.HD', 'MT0.N2.LD', 'MT0.N2.HD', 'Tubulins_I', 'Tubulins_II', 'Bundled_Tubes_High_Density', 'Bundled_Tubes_Long_Sequence']
    sets = [s for s in order if s in sets]
    rows, tally = [], {}
    for s in sets:
        nn, fr, hi, hg = sm[s + '|nn'], sm[s + '|frame'], sm[s + '|hit_iss'], sm[s + '|hit_gr']
        frames = np.unique(fr)
        starts = np.arange(0, len(frames), SMLM_BLOCK)
        blocks = [np.flatnonzero(np.isin(fr, frames[a:a + SMLM_BLOCK])) for a in starts]
        rng = np.random.default_rng(BOOT_SEED)
        pick = rng.integers(0, len(blocks), size=(BOOT, len(blocks)))
        for lo, hi_ in SMLM_BANDS_NM:
            band = (nn >= lo) & (nn < hi_)
            if band.sum() < 30:
                continue
            ni = np.array([band[b].sum() for b in blocks], float)
            ci = np.array([hi[b][band[b]].sum() for b in blocks], float)
            cg = np.array([hg[b][band[b]].sum() for b in blocks], float)
            bs = np.array([(ci[p].sum() - cg[p].sum()) / max(ni[p].sum(), 1.0) for p in pick])
            d = (hi[band].mean() - hg[band].mean())
            l_, h_ = np.percentile(bs, [2.5, 97.5])
            v = verdict(l_, h_)
            tally[(lo, v)] = tally.get((lo, v), 0) + 1
            rows.append([s, '[%g, %s)' % (lo, 'inf' if hi_ > 1e8 else '%g' % hi_), int(band.sum()), '%.3f' % hi[band].mean(), '%.3f' % hg[band].mean(), '%+.3f [%+.3f, %+.3f]' % (d, l_, h_), v])
    L = md_table(['set', 'nearest neighbour (nm)', 'emitters', 'recall dual ISS', 'recall greedy', 'dual - greedy', 'verdict'], rows) + ['']
    L += ['Tally of the verdicts (dual ahead = win, greedy ahead = loss) by band: ' + '; '.join(
        '[%g, %s): %s' % (lo, 'inf' if hi_ > 1e8 else '%g' % hi_, ', '.join('%d %s' % (tally.get((lo, v), 0), v) for v in ('win', 'tie', 'loss'))) for lo, hi_ in SMLM_BANDS_NM), '']
    return L


def t_lines(T):
    """Part T: noise-free pairs, per resolution limit and separation, over the ratios: the number of the T_N phases (of each ratio) in which the entry sequences are the same / the blocks after
    two entries are the same set / the dual's block is exact / the greedy block is exact."""
    L = []
    for i, r_ex in enumerate(T_REX):
        rows = []
        for j, d in enumerate(T_DS):
            cells = []
            for k, ratio in enumerate(T_RATIOS):
                q = T[i, j, k]
                cells.append('%d / %d / %d / %d' % tuple(int(q[:, m].sum()) for m in range(4)))
            rows.append(['%g' % d] + cells)
        L += ['Resolution limit r_ex = %g (each cell: same entry sequence / same block / dual exact / greedy exact, out of %d phases):' % (r_ex, T_N), '']
        L += md_table(['separation'] + ['ratio %g' % r for r in T_RATIOS], rows) + ['']
    return L


def test_lines(raw):
    """The held-out run: for each of the regimes P (close pairs) and Q (crowded windows) the pooled comparisons, with the plan's verdict rule, and the equivalence in E."""
    L, res = [], {}
    for g, title in (('P', 'close pairs'), ('Q', 'crowded windows'), ('W', 'four features of unequal widths, widths free')):
        if '%s|0|E' % g not in raw:
            continue
        js = cell_keys(raw, g)
        rows = []
        for j in js:
            E, fam = raw['%s|%d|E' % (g, j)], fam_of(raw, g, j)
            rows.append([fam_label(fam), len(E)] + [pct(np.mean(exact_of(E, m, r))) for m, r in (('dual', 'bic'), ('greedy', 'bic'), ('greedy', 'bicbw'), ('sfw', 'bic'), ('dual', 'oracleK'), ('greedy', 'oracleK'))])
        L += ['### Regime %s (%s): exact recovery (%%) under the BIC rule, by cell' % (g, title), '']
        L += md_table(['cell', 'scenes', 'dual', 'greedy', 'greedy + backward step', 'sFW', 'dual, first block with K units', 'greedy, first block with K units'], rows) + ['']
        rows, rg, num = [], {}, {}
        for cid, label, (m1, r1), (m2, r2) in COMPARISONS:
            a, b = pooled(raw, g, m1, r1), pooled(raw, g, m2, r2)
            mean, lo, hi = paired_boot((a - b) * 100.0)
            v = verdict(lo, hi)
            rg[cid] = v
            num[cid] = (mean, lo, hi)
            rows.append([g + cid, label, pct(np.mean(a)), pct(np.mean(b)), '%+.1f [%+.1f, %+.1f]' % (mean, lo, hi), v, comparison_status(g, cid)])
        L += ['### Regime %s: pooled paired comparisons (%d scenes; difference in points of exact recovery, first arm minus second; 95 %% bootstrap interval over scenes, %d resamples, seed %d)' % (g, len(a), BOOT, BOOT_SEED), '']
        L += md_table(['id', 'comparison', 'first arm (%)', 'second arm (%)', 'difference [95 % interval]', 'verdict', 'status'], rows) + ['']
        rd, rgr = pooled_rms(raw, g, 'dual', 'bic'), pooled_rms(raw, g, 'greedy', 'bic')
        both = np.isfinite(rd) & np.isfinite(rgr)
        mean, lo, hi = paired_boot(rd[both] - rgr[both])
        L += ['Localisation (%s11): over the %d scenes in which both orders are exact under BIC, the position rms (nominal widths) is %.3f for the dual and %.3f for the greedy order; paired difference %+.3f [%+.3f, %+.3f] (%s; lower is better, so a negative difference favours the dual).' % (
            g, both.sum(), np.mean(rd[both]), np.mean(rgr[both]), mean, lo, hi, 'win' if hi < 0 else 'loss' if lo > 0 else 'tie'), '']
        need = ('1', '3') if g == 'W' else ('1', '3', '4')
        rg['regime'] = 'yes' if all(rg[k] == 'win' for k in need) else 'no'
        L += ['**Plan rule: the dual order matters in %s iff %s are all wins%s.  Result: %s** (%s).' % (
            g, ', '.join(g + k for k in need), '' if g != 'W' else ' (the sliding Frank-Wolfe comparator has a fixed width and is shown, not counted)', rg['regime'], ', '.join('%s%s %s' % (g, k, rg[k]) for k in need)), '']
        rg['num'] = num
        res[g] = rg
    if 'E|0|E' in raw:
        rows, k_all, n_all = [], 0, 0
        for j in cell_keys(raw, 'E'):
            E, A = raw['E|%d|E' % j], raw['E|%d|A' % j]
            fam = fam_of(raw, 'E', j)
            k, n = int(A[:, 0].sum()), len(A)
            k_all, n_all = k_all + k, n_all + n
            rows.append([fam_label(fam), n, pct(np.mean(exact_of(E, 'dual', 'bic'))), pct(np.mean(exact_of(E, 'greedy', 'bic'))), '%d / %d' % (k, n), '%d / %d' % (int(A[:, 1].sum()), n)])
        lo, hi = clopper_pearson(k_all, n_all)
        res['equiv'] = 'yes' if lo >= 0.98 else 'no'
        res['equiv_n'] = (k_all, n_all)
        L += ['### Regime E (separation at or beyond the theory\'s entry threshold): do the two orders return the same block?', '']
        L += md_table(['cell', 'scenes', 'dual exact (%)', 'greedy exact (%)', 'same block, BIC', 'same block, penalty 16'], rows) + ['']
        L += ['Pooled: the BIC blocks of the two orders are the same set in %d of %d scenes (%.1f %%, Clopper-Pearson 95 %% interval [%.1f, %.1f]).  **Plan rule: equivalent iff the lower end is at least 98 %%.  Result: %s.**' % (
            k_all, n_all, 100.0 * k_all / n_all, 100 * lo, 100 * hi, res['equiv']), '']
    return L, res


def twod_sep(twod, d):
    """(mean, lo, hi) in points of the paired difference dual minus greedy in resolved fraction at the single separation d (widths), the photon levels and the scenes pooled; None if d is not a separation of the arm."""
    if d not in TWOD_SEPS:
        return None
    r = (twod[TWOD_SEPS.index(d)] == 1)
    return paired_boot(((r[..., 0].astype(float) - r[..., 1]) * 100.0).ravel())


def twod_lines(twod):
    """The two-dimensional arm: resolved fractions (class 1) per cell and the pooled paired differences of the close and the mid groups."""
    sys.path.insert(0, HERE)
    import exp_order_pairs as OP
    rows = []
    for i, d in enumerate(TWOD_SEPS):
        for k, ip in enumerate(TWOD_PHOTONS):
            r = twod[i, k] == 1
            diff = (r[:, 0].astype(float) - r[:, 1]) * 100.0
            mean, lo, hi = paired_boot(diff)
            rows.append(['%g' % d, '%g' % OP.PHOTONS[ip], len(r), pct(r[:, 0].mean()), pct(r[:, 1].mean()), '%+.1f [%+.1f, %+.1f]' % (mean, lo, hi), verdict(lo, hi)])
    L = md_table(['separation (widths)', 'photons (brighter)', 'scenes', 'resolved, dual (%)', 'resolved, greedy (%)', 'dual - greedy [95 % interval]', 'verdict'], rows) + ['']
    out = {}
    for name, grp in (('close', TWOD_CLOSE), ('mid', TWOD_MID)):
        ii = [TWOD_SEPS.index(d) for d in grp]
        r = (twod[ii] == 1)
        diff = ((r[..., 0].astype(float) - r[..., 1]) * 100.0).ravel()
        mean, lo, hi = paired_boot(diff)
        out[name] = verdict(lo, hi)
        L += ['Group %s (separations %s widths, %d scenes pooled): resolved %.1f %% (dual) against %.1f %% (greedy), difference %+.1f [%+.1f, %+.1f] points: **%s**.' % (
            name, ', '.join('%g' % d for d in grp), diff.size, 100 * r[..., 0].mean(), 100 * r[..., 1].mean(), mean, lo, hi, out[name]), '']
    t3 = twod_sep(twod, 3.0)
    if t3 is not None:
        L += ['Separation 3 widths alone (both photon levels pooled): dual minus greedy %+.1f [%+.1f, %+.1f] points: **%s**.  Two dimensions at 3.5 widths and more were not run.' % (t3 + (verdict(t3[1], t3[2]),)), '']
    return L, out


def dev_tally(dev):
    """Counts over the development cells: for each comparison (dual minus greedy under BIC; both with the backward step; dual minus sFW) the cells with the difference at least two standard errors above 0, below 0, and neither,
    over the cells whose standard error is positive (a zero standard error is a cell in which the two never differ)."""
    out = {}
    for name, (a, b, rule) in (('dual - greedy, BIC', ('dual', 'greedy', 'bic')), ('dual + backward step - greedy + backward step', ('dual', 'greedy', 'bicbw')), ('dual - sFW, BIC', ('dual', 'sfw', 'bic'))):
        up = dn = nl = tot = 0
        for g in DEV_GRIDS:
            for j in cell_keys(dev, g):
                E = dev['%s|%d|E' % (g, j)]
                m, se = diff_stats(E, a, b, rule)
                tot += 1
                if se == 0.0:
                    nl += 1
                elif m >= 2 * se:
                    up += 1
                elif m <= -2 * se:
                    dn += 1
        out[name] = (up, dn, nl, tot)
    return out


def summary_lines(raw, res, dev):
    """The generated summary at the top of the report."""
    L = []
    for g, title in (('P', 'close pairs'), ('Q', 'crowded windows'), ('W', 'unequal widths')):
        if g in res:
            r = res[g]
            need = ('1', '3') if g == 'W' else ('1', '3', '4')
            order_ok, sfw_ok = all(r[k] == 'win' for k in need if k != '4'), (g == 'W' or r['4'] == 'win')
            cause = ('met' if r['regime'] == 'yes' else 'the failing part is (ii) alone' if order_ok else 'the failing part is (i) alone' if sfw_ok else 'neither (i) nor (ii) is met')
            part2 = ('(ii) the comparison with sliding Frank-Wolfe, a different algorithm and not an order comparison (%s4): a %s.' % (g, r['4']) if g != 'W' else
                     '(ii) the comparison with sliding Frank-Wolfe (W4-W8) is not an order comparison and, with its fixed width, misspecified here: not to be quoted.')
            L.append('* Regime %s (%s): (i) the order effect with everything else equal (same dictionary, resolution limit, block fit and stopping rule; the same backward step on both, %s1): a %s; as the pipelines stand (%s2) a %s; '
                     'against greedy + backward step (%s3) a %s; first block with K units (%s9, oracle K) a %s; %s  Registered plan rule (%s all wins): **%s**, %s.  %s10 (each order\'s own stopping rule) is not an order comparison and is not to be quoted.' % (
                         g, title, g, r['1'], g, r['2'], g, r['3'], g, r['9'], part2, ', '.join(g + k for k in need), r['regime'], cause, g))
    if 'equiv' in res:
        k_, n_ = res.get('equiv_n', (None, None))
        t3 = twod_sep(raw['twod'], 3.0) if 'twod' in raw else None
        L.append('* Regime E (separation 3.5 and 5 widths): the two orders return the same block in the BIC rule: equivalent (plan rule) **%s**%s.  Scope: one-dimensional pairs with two features, the final BIC blocks (sorted centres within %g widths), not the '
                 'entry sequences and not images%s.' % (res['equiv'], '' if k_ is None else ' (%d of %d scenes)' % (k_, n_), SAME_TOL,
                                                        '' if t3 is None else ': in two dimensions at 3 widths the dual minus greedy difference is %+.1f [%+.1f, %+.1f] points (a %s) and 3.5 widths and more were not run' % (t3 + (verdict(t3[1], t3[2]),))))
    if 'twod' in res:
        L.append('* Two-dimensional arm (fresh scenes of exp_order_pairs, BIC stop): close separations (1.2, 1.6 widths) dual vs greedy **%s**, mid separations (2.5, 3.0 widths) **%s**.' % (res['twod']['close'], res['twod']['mid']))
    if 'T' in raw:
        T = raw['T']
        for i, r_ex in enumerate(T_REX):
            far = [T_DS.index(d) for d in T_DS if d >= 5.0]
            near = [T_DS.index(d) for d in T_DS if d <= 3.5]
            L.append('* Part T (noise-free, r_ex %g): the entry sequences of the two orders are the same in %d of %d phases at separations of 5 widths and more, and in %d of %d phases at separations of 3.5 widths and less.' % (
                r_ex, int(T[i, far, :, :, 0].sum()), T[i, far, :, :, 0].size, int(T[i, near, :, :, 0].sum()), T[i, near, :, :, 0].size))
    if dev:
        for name, (up, dn, nl, tot) in dev_tally(dev).items():
            L.append('* Development scan, %s: of %d cells, %d with the difference at least two standard errors above 0, %d at least two below, %d never differ (zero standard error), %d neither.' % (name, tot, up, dn, nl, tot - up - dn - nl))
    return L + ['']


def conclusion_lines(res):
    """The conclusions, as sentences whose verdict words are taken from the results; the only numbers stated are the pooled differences of P1, P2, Q1, Q2 (res[g]['num']) and the referee\'s rerun of the backward-step cap, which is quoted and labelled as the referee\'s."""
    L = []
    if res.get('equiv') == 'yes':
        L.append('1. **Equivalence.**  For separations of 3.5 widths and more (regime E) the two orders select the same block in at least 98 %% of the scenes (plan rule met), and in part T the noise-free entry sequences coincide from about 5 widths, the order of the continuum '
                 'entry threshold of the Gram (4.7 to 5.8 widths).  This is what the structure of the two orders predicts: when the Gram matrix of the features is nearly diagonal, the instantaneous correlation of every feature not yet explained is constant in time (fitting '
                 'the others does not change it), so the accumulated dual of a feature is its constant correlation times the elapsed time, and both orders sort the same numbers, the amplitudes.  They can differ only through the coherence between features.  Scope: this is for one-dimensional pairs with two features (amplitude ratios 1 and 10, separations 3.5 and 5 widths)%s and for the FINAL BIC blocks, '
                 'not for the entry sequences (in the noise-free part T they differ at 3.5 widths and less, see the summary) and not for images: in the two-dimensional arm the dual order against the greedy order is a %s at the mid separations (2.5, 3 widths), and two dimensions at 3.5 widths and more were not run.' %
                 (('' if 'equiv_n' not in res else ' (%d of %d scenes returned the same block)' % res['equiv_n']), res.get('twod', {}).get('mid', 'not run')))
    elif 'equiv' in res:
        L.append('1. **Equivalence.**  The plan rule for equivalence in regime E was NOT met (see its table).')
    for g, name in (('P', 'close pairs'), ('Q', 'crowded windows')):
        if g not in res:
            continue
        r = res[g]
        L.append('2%s. **Regime %s (%s).**  As the pipelines stand the dual order against the bare greedy order is a %s (%s2); with the same backward step on both, the order alone is a %s (%s1); against greedy + backward step a %s (%s3); with no stopping rule at all (the first block with K units) a %s (%s9).  '
                 'This is the order effect, everything else equal, (%s1) and (%s2); the comparison with the sliding Frank-Wolfe solver of the BLASSO is a comparison with a different algorithm, not an order comparison, and the dual is a %s there (%s4).  '
                 'The oracle-K comparison (%s9) uses the true number of units.  The plan rule for the claim that the dual order matters in %s is %s met.' % (
                     'ab'[g == 'Q'], g, name, r['2'], g, r['1'], g, r['3'], g, r['9'], g, g, g, r['4'], g, g, g, '' if r['regime'] == 'yes' else 'NOT'))
    if 'W' in res:
        r = res['W']
        L.append('3. **Regime W (four features of unequal widths, widths free).**  Dual against greedy as the pipelines stand: %s (W2); with the same backward step: %s (W1); dual against greedy + backward step: %s (W3); first block with K units: %s (W9).  The plan rule (W1 and W3) is %s met.' % (
            r['2'], r['1'], r['3'], r['9'], '' if r['regime'] == 'yes' else 'NOT'))
    if 'twod' in res:
        L.append('4. **Two dimensions (fresh scenes of the SMLM close-pair experiment).**  Close separations: dual against greedy a %s; mid separations: a %s.' % (res['twod']['close'], res['twod']['mid']))
    L.append('5. **Mechanism (a hypothesis, not a finding).**  The data are CONSISTENT WITH early correlation memory being what lets the dual order lead where it does; this was not tested: there is no ablation (a leaky or reset accumulator was not run), and the evidence is one planted '
             'scene and an inspection of the first merged two-dimensional scenes (not scripted).  The picture: when two features are close the first atom is fitted to their union and absorbs part of the second; the residual correlation at the second feature is then no larger than the noise maxima '
             'elsewhere in the window, and matching pursuit would enter a noise peak first (the planted scene of the tests: seed 100002 of the pair at 1.6 widths).  The accumulated dual would still carry the early correlation, which both features contributed to, and so enter the second feature first.  '
             'A backward step (removing a unit whose removal costs less than the penalty) removes most of the noise units that matching pursuit then carries, which is consistent with the lead shrinking from the pipelines-as-they-stand comparison to the same-backward-step comparison.  The same memory '
             'could be a liability when a unit is fitted exactly: the duals of its neighbours have accumulated to nearly the wall and, without a resolution limit, the dual enters ghosts next to it (the planted wide-Gaussian test shows the ghosts; the resolution limit of the pipelines suppresses '
             'them).  A second possible way the memory hurts, suggested by the two-dimensional arm at 2.5 to 3 widths (an inspection of the first merged scenes at 3 widths, hence a hypothesis), is that '
             'the first atom, entered at the midpoint of a resolved pair, slides towards one emitter while the accumulated (possibly stale) dual of the other emitter\'s lattice points keeps its head start: the dual would enter the already-covered emitter next and leave the first atom as a weak third unit '
             '(a split), where matching pursuit, which sees the residual after the slide, enters the uncovered emitter.  The post hoc SMLM table below is the only check on real frames available here and is not a confirmation.')
    nP, nQ = res.get('P', {}).get('num', {}), res.get('Q', {}).get('num', {})
    red = ''
    if all(k in d for d in (nP, nQ) for k in ('1', '2')):
        v = [nP['1'][0], nQ['1'][0]]
        red = ('  The backward step reduces the lead over the bare greedy order but does not remove it: from %+.1f (P2) and %+.1f (Q2) points to %+.1f (P1) and %+.1f (Q1) points, that is to about %d-%d points (P1 %s, Q1 %s).' % (
            nP['2'][0], nQ['2'][0], nP['1'][0], nQ['1'][0], round(min(v)), round(max(v)), res['P']['1'], res['Q']['1']))
    L.append('6. **What the paper may claim.**  The dual order of entry is not what makes the method better than the classical solvers: in the regimes where it leads the bare greedy order, the lead is modest once a backward step is allowed and it is exceeded by the sliding Frank-Wolfe '
             'solver of the BLASSO (a different algorithm: that is not an order comparison).%s  The accumulated-dual order is the object that the small-initialisation training dynamics computes and that the theory analyses (entry times as a statistic, the stopping time); it is an analysis tool and a modest advantage in crowded, near-resolution scenes, '
             'not a source of accuracy over the best classical comparator.' % red)
    L.append('7. **Where the plan\'s stated expectations were not met.**  The comparisons P1 and P3 were predicted ties and are small wins; in part T a few cells fall short of the counts the plan quoted for 5 widths and above (and one cell at 3.5 widths and below exceeds them), the pooled counts being in the summary; the plan\'s other stated expectations were met.')
    L.append('8. **Design caveats (the referee\'s review; not corrected in the stored runs).**  The backward step is applied only to blocks of at most K_true + %d units (`apply_rules` uses `len(c_true)`), so the true number of features leaks into the backward-step arms (comparisons %s).  '
             'The referee\'s rerun without the cap on held-out scenes (180 P, 120 W, 80 Q; not recomputed here) left P identical and moved greedy + backward step from 78 to 82 of 120 in W (dual + backward step from 81 to 82) and from 8 to 10 of 80 in Q: P1, Q1 and W1 are inflated by about 0 to 3 points.  '
             'The comparisons marked oracle K use the true number of units by construction, and the comparisons against the sliding Frank-Wolfe solver and the comparisons with each method\'s own stopping rule are not order comparisons (see the status column).' % (PRUNE_EXTRA, ', '.join(BW_IDS)))
    return L + ['']


def report(path=None, dev_path=None):
    path, dev_path = path or RAW, dev_path or DEV_RAW
    raw, dev = load_raw(path), load_raw(dev_path)
    import platform
    import scipy
    L = ['# Is there a regime in which the dual order of entry matters?  The sliding inverse scale space against sliding matching pursuit and the classical comparators', '',
         'Generated by `experiments/exp_dualregime.py` from `dualregime_raw.npz` (the pre-registered held-out run, part T and the two-dimensional arm) and `dualregime_dev_raw.npz` (the development scan); every number below is printed by the script.  Setting, '
         'engine and metrics: the script\'s docstring.  The plan is `dualregime_plan.md`, committed before the held-out run was scored.  Python %s, numpy %s, scipy %s.  Paired differences are first arm minus second, in points of '
         'exact-recovery rate; a **win** (for the first arm) is an interval that excludes 0 above, a **loss** below, a **tie** otherwise.' % (platform.python_version(), np.__version__, scipy.__version__), '']
    res, tl, tw = {}, [], []
    if any('%s|0|E' % g in raw for g in TEST_GRIDS):
        tl, res = test_lines(raw)
    if 'twod' in raw:
        tw, res['twod'] = twod_lines(raw['twod'])
    L += ['## Summary (generated)', ''] + summary_lines(raw, res, dev)
    if 'T' in raw:
        L += ['## Part T. Noise-free pairs: when do the two orders enter the same atoms?', '',
              'Two Gaussians of the dictionary\'s width at separation d (widths), amplitude ratio a1/a2, at %d random phases of the lattice, no noise; both orders run for two entries (the resolution limit r_ex as stated).' % T_N, '']
        L += t_lines(raw['T'])
    if any('%s|0|E' % g in raw for g in TEST_GRIDS):
        L += ['## Held-out run (pre-registered; seeds from %d)' % TEST_SEED0, '']
        L += tl
    if 'twod' in raw:
        L += ['## Two-dimensional arm (pre-registered; fresh scenes of `exp_order_pairs`, EM-CCD frames, equal fluxes, stop %s %g, seeds from %d)' % (TWOD_STOP[0], TWOD_STOP[1], TWOD_SEED0), '']
        L += tw
    sm = {k[5:]: v for k, v in raw.items() if k.startswith('smlm|')}
    if res:
        L += ['## Conclusions (generated; verdict words from the tables)', ''] + conclusion_lines(res)
    if sm:
        L += ['## POST HOC: the SMLM test frames, by the distance to the nearest neighbour', '',
              'Nothing here was pre-registered: the stored test detections of `results/smlm/controls_raw.npz` (the dual ISS at its E-tuned point and the greedy order at the same point, read only) are scored on the ground-truth emitters of the '
              'test frames, each with the distance to its nearest neighbour in the frame; an emitter is recalled if a detection is matched to it one to one within %g nm.  Intervals: block bootstrap over the test frames (blocks of %d frames, %d resamples, seed %d).' % (
                  SMLM_TOL_NM, SMLM_BLOCK, BOOT, BOOT_SEED), '']
        L += smlm_lines(sm)
    if dev:
        dl, flags = dev_lines(dev)
        L += ['## Development scan (exploratory: every cell is a look; nothing here is a test)', '',
              'Cells of the grids A (pairs), B (random clusters), C (random clusters with unequal widths, widths free in [%g, %g]), D (crowded windows), %s.  Each cell: exact-recovery rate under the BIC rule of the dual, the greedy order and the sliding Frank-Wolfe '
              'comparator, then the paired difference dual minus the other with its standard error over scenes (points; in brackets).  "both + backward step" is the dual and the greedy order each followed by backward elimination of every block.' % (
                  WIDTHS[0], WIDTHS[1], 'DEV seeds from %d' % DEV_SEED0), '']
        L += dl
    return '\n'.join(L) + '\n'


# ---- raw files ----------------------------------------------------------------------------------------------------------------------------------------------
def load_raw(path):
    """The raw file as a dict (empty if absent)."""
    if not os.path.exists(path):
        return {}
    with np.load(path, allow_pickle=False) as z:
        return {k: z[k] for k in z.files}


def save_raw(path, new):
    """Merge `new` into the raw file at `path` (atomic: temporary file, then rename)."""
    d = load_raw(path)
    d.update(new)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + '.tmp.npz'
    np.savez_compressed(tmp, **d)
    os.replace(tmp, path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--scan', action='store_true', help='the development scan (DEV_GRIDS) into the dev raw')
    ap.add_argument('--run', action='store_true', help='the held-out run: parts T, P, E, Q, W, 2d (smlm is the POST HOC part, --smlm)')
    ap.add_argument('--smlm', action='store_true', help='the POST HOC analysis of the SMLM test detections (reads results/smlm/controls_raw.npz)')
    ap.add_argument('--parts', default='T,P,E,Q,W,2d', help='with --run: which parts, comma separated')
    ap.add_argument('--grids', default='C,D,A,B', help='with --scan: which grids, comma separated')
    ap.add_argument('--workers', type=int, default=int(os.environ.get('DR_WORKERS', '1')))
    ap.add_argument('--from-raw', action='store_true', help='regenerate the md from the raw files')
    A = ap.parse_args()
    if A.scan:
        for g in A.grids.split(','):
            save_raw(DEV_RAW, run_grid(DEV_GRIDS, DEV_SEED0, 1000, [g], A.workers))
    if A.run:
        parts = A.parts.split(',')
        if 'T' in parts:
            save_raw(RAW, {'T': part_t()})
        for g in ('P', 'E', 'Q', 'W'):
            if g in parts:
                save_raw(RAW, run_grid(TEST_GRIDS, TEST_SEED0, 10000, [g], A.workers))
        if '2d' in parts:
            save_raw(RAW, {'twod': part_2d(A.workers)})
    if A.smlm:
        save_raw(RAW, {'smlm|' + k: v for k, v in part_smlm().items()})
    if A.from_raw:
        with open(MD, 'w') as f:
            f.write(report())
        print('wrote', MD)


if __name__ == '__main__':
    main()
