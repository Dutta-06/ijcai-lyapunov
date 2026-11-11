"""
The separated regime for two features of equal amplitude (round theory3, direction Y; statement S6' of the round note): a closed-form sufficient condition on (Delta, h, a_1, a_2) for the
hypotheses (a), (b) of Theorem K, checked against the exact recursion of the limit object.

Setting: K = 2, n = 1, continuum Gram Z Gamma(c - c'), equal walls, |A_1| = |A_2| (any signs), the starting kernels on a lattice of spacing h (a sub-dictionary that contains the nearest kernel
of each feature is allowed).  a_m in [0, h/2] is the distance of feature m to the lattice, ordered a_1 <= a_2 (a_1 < a_2 < h/2: the generic case; a_1 = a_2 is the tie of the two entry times,
a_2 = h/2 the tie of two kernels of one feature).  Write g = Gamma, phi(D) = g(D) + g(Delta - D).  Constants of Prop. A for K = 2, R = 1, Omega = 0 at the radius r = h/2 + 0.02
(which exceeds every entrant offset, so the same box serves both layers): d-tilde, E_v with Q = 1 - g(d-tilde) + E_v, and
    rho(x) = d-tilde psi_sup(x - d-tilde) + Q g(x - d-tilde)     (x >= d-tilde; nonincreasing)
bounding |g(x) - v g(x - d)| for the first block (|d| <= d-tilde, |1 - v| <= Q).  The sufficient condition (all conditions are explicit numbers):
  (S0)  (L), (C), (N) of Prop. A at r = h/2 + 0.02 (ball arithmetic), and phi nonincreasing on [h, Delta/2] (Lemma phi);
  (S1)  P_- := g(a_1) - g(Delta - a_1)  >  L1 := max{phi(a_2), phi(h - a_1), phi(h - a_2), phi(h)}                      [first entry: the nearest kernel of feature 1, uniquely];
  (S2)  q_2 := g(a_2) - rho(Delta - a_2) > 0;
  (S3)  T_* := (1 - (g(a_2) - g(Delta - a_2))/P_+) / q_2   <   min over the competitors C1..C4 of (1 - u_bar)/r_bar          [second entry: the nearest kernel of feature 2, uniquely],
        with P_+ = g(a_1) + g(Delta - a_1) and  C1: (u_bar, r_bar) = (phi(h - a_2)/P_-, g(h - a_2) + rho(Delta - h + a_2)),  C2: (phi(h - a_1)/P_-, g(Delta - h + a_1) + rho(h - a_1)),
        C3: (phi(h)/P_-, g(h) + rho(Delta/2)),  C4: (phi(h)/P_-, g(Delta/2) + rho(h)).
Then Theorem K (a), (b) hold: the limit object has exactly the two entries (the nearest kernels, the one of smaller offset first), no exit, final block the target.
Checked here against `exp_theory_t1.limit_path` (the exact recursion) on random lattice phases, separations and signs: the condition never holds when the recursion fails.
"""
import os
import sys

os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np

HERE = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..'))
import exp_theory_t1 as T          # noqa: E402

DELTAS = (4.7, 4.9, 5.2, 5.6, 6.0, 7.0, 8.0)
HS = (0.8333, 0.5)
N_TRIALS = 300
N_GRID = 120
SEED = 424242
AMP = 0.7


def g(x):
    return T.Gam(x)


def phi(D, Delta):
    return g(D) + g(Delta - D)


def box_for(Delta, h):
    """Prop. A's box for K = 2, R = 1, Omega = 0 at r = h/2 + 0.02: (L), (C), (N) certain (ball arithmetic), or None."""
    r = h / 2 + 0.02
    b = T.box_t1(2, 1.0, Delta, r, 0.0, 0.05, 0.01)
    if 'nu' in b and b['nu'] > 0 and b['vlo'] > 0 and b['Cmargin'] > 0 and b['vD_lo'] > 0:
        b = dict(b)
        b['r'] = r
        return b
    return None


def lemma_phi(Delta, h):
    """phi is nonincreasing on [h, Delta/2]: psi(D) >= psi(Delta - D) there.  For D >= sqrt 2 both arguments are beyond the maximum of psi and D <= Delta - D; for D < sqrt 2 it suffices that
    psi(h) >= psi(Delta - sqrt 2) (psi increases on [0, sqrt 2] and decreases beyond; Delta >= 2 sqrt 2).  Returns that inequality (certain in floats with a 1e-12 margin)."""
    if Delta < 2 * np.sqrt(2.0) or h >= np.sqrt(2.0):
        return False
    return bool(T.psi(h) > T.psi(Delta - np.sqrt(2.0)) + 1e-12)


def rho_bar(x, dt, Q):
    """|g(x) - v g(x - d)| <= rho(x) for |d| <= dt, |1 - v| <= Q, x >= dt."""
    xx = np.maximum(np.asarray(x, float) - dt, 0.0)
    return dt * T.psi_sup(xx) + Q * g(xx)


def conditions(Delta, h, a1, a2, box=None):
    """The sufficient condition of the note for the target with offsets a1 <= a2 (both in [0, h/2]) at the separation Delta (the lattice spacing h).  Returns a dict of booleans and margins."""
    if box is None:
        box = box_for(Delta, h)
    out = dict(s0=False, s1=False, s2=False, s3=False, ok=False, P_minus=np.nan, L1=np.nan, q2=np.nan, T_star=np.nan, T_comp=np.nan)
    if box is None or not lemma_phi(Delta, h) or not (0 <= a1 < a2 < h / 2):
        return out
    out['s0'] = True
    dt, Ev = box['dtilde'], box['Ev']
    Q = 1 - float(g(dt)) + Ev
    Pm, Pp = float(g(a1) - g(Delta - a1)), float(g(a1) + g(Delta - a1))
    L1 = float(max(phi(a2, Delta), phi(h - a1, Delta), phi(h - a2, Delta), phi(h, Delta)))
    out.update(P_minus=Pm, L1=L1)
    out['s1'] = bool(Pm > L1)
    rho = lambda x: float(rho_bar(x, dt, Q))
    q2 = float(g(a2) - rho(Delta - a2))
    out['q2'] = q2
    out['s2'] = bool(q2 > 0)
    if not (out['s1'] and out['s2']):
        return out
    Tstar = float((1 - (g(a2) - g(Delta - a2)) / Pp) / q2)
    comps = [((phi(h - a2, Delta)) / Pm, float(g(h - a2)) + rho(Delta - h + a2)),
             (float(phi(h - a1, Delta)) / Pm, float(g(Delta - h + a1)) + rho(h - a1)),
             (float(phi(h, Delta)) / Pm, float(g(h)) + rho(Delta / 2)),
             (float(phi(h, Delta)) / Pm, float(g(Delta / 2)) + rho(h))]
    Tc = min((1 - u) / r for (u, r) in comps)
    out.update(T_star=Tstar, T_comp=float(Tc))
    out['s3'] = bool(Tstar < Tc)
    out['ok'] = bool(out['s1'] and out['s2'] and out['s3'])
    return out


def build_target(Delta, h, a1, a2, sg1, sg2, sgn_sep, flip):
    """Two equal-amplitude features at lattice offsets +-a1, +-a2 (signs sg1, sg2) and the lattice h Z; the features are N = round(Delta/h) cells apart; flip swaps which feature is the left one.
    Returns (cs, As, c0, actual separation, (a_left, a_right))."""
    N = int(round(Delta / h))
    left_off, right_off = (sg1 * a1, sg2 * a2) if not flip else (sg2 * a2, sg1 * a1)
    cs = np.array([left_off, N * h + right_off])
    sep = float(cs[1] - cs[0])
    As = AMP * np.array([1.0, sgn_sep])
    c0 = np.arange(np.floor((cs.min() - 8) / h) * h, cs.max() + 8, h)
    return cs, As, c0, sep


def oracle(cs, As, c0, a1, a2):
    """The exact recursion of Theorem K on the target: True if the limit path has exactly the two entries (the nearest kernels, the one of the smaller offset first), both signs right, both entrants
    unique (relative gap > 0), the final block the target."""
    L = T.limit_path(cs, As, c0)
    ev = L['events']
    if len(ev) != 2:
        return False
    near = [int(np.argmin(np.abs(c0 - c))) for c in cs]
    off = [abs(c0[near[m]] - cs[m]) for m in range(2)]
    first = int(np.argmin(off))
    want = [near[first], near[1 - first]]
    if [e['kernel'] for e in ev] != want:
        return False
    if not all(e['sign_ok'] for e in ev) or not all(e['gap'] > 1e-9 for e in ev):
        return False
    fin = ev[1]
    return bool(np.allclose(np.sort(fin['c']), np.sort(cs), atol=1e-8) and np.allclose(np.sort(np.abs(fin['w'])), np.sort(np.abs(As)), atol=1e-8))


def run_part(quick=False):
    """Part S: for each lattice spacing h and separation Delta, (i) the region of the condition in the offset square (grid), (ii) random targets: the condition against the exact recursion
    (violations = condition true and recursion not as the theorem says: must be 0), the share of recursion-successes covered by the condition, (iii) the benchmark's K = 2, R = 1 targets."""
    rng = np.random.default_rng(SEED)
    deltas = DELTAS[::3] if quick else DELTAS
    hs = HS[:1] if quick else HS
    nt = 30 if quick else N_TRIALS
    ng = 30 if quick else N_GRID
    rows, reg, bench = [], [], []
    for h in hs:
        for Delta in deltas:
            box = box_for(Delta, h)
            # region of the condition on a grid of offsets (a1 < a2 in (0, h/2))
            ax = (np.arange(ng) + 0.5) / ng * (h / 2)
            cnt, tot = 0, 0
            if box is not None and lemma_phi(Delta, h):
                for i in range(ng):
                    for j in range(i + 1, ng):
                        tot += 1
                        cnt += int(conditions(Delta, h, ax[i], ax[j], box)['ok'])
            else:
                tot = ng * (ng - 1) // 2
            reg.append([h, Delta, float(box is not None), cnt / tot if tot else np.nan, tot])
            # random targets against the oracle
            n_cond = n_ok = n_viol = n_both = n_trial = 0
            for _ in range(nt):
                a1, a2 = np.sort(rng.uniform(0.0, h / 2, 2))
                sg1, sg2 = rng.choice([-1.0, 1.0], 2)
                sgn_sep = float(rng.choice([-1.0, 1.0]))
                flip = bool(rng.integers(0, 2))
                cs, As, c0, sep = build_target(Delta, h, a1, a2, sg1, sg2, sgn_sep, flip)
                c = conditions(sep, h, a1, a2)
                ok = oracle(cs, As, c0, a1, a2)
                n_trial += 1
                n_cond += int(c['ok'])
                n_ok += int(ok)
                n_both += int(c['ok'] and ok)
                n_viol += int(c['ok'] and not ok)
            rows.append([h, Delta, n_trial, n_cond, n_ok, n_both, n_viol])
    # the benchmark's K = 2, R = 1 targets (exp_theory_t1.C_TARGETS with R = 1)
    for tgt in [t for t in T.C_TARGETS if t[0] == 2 and t[2] == 1.0][:2 if quick else 10]:
        cs, As, c0 = T.make_target(*tgt)
        h = float(c0[1] - c0[0])
        a = [float(np.min(np.abs(c0 - c))) for c in cs]
        a1, a2 = sorted(a)
        sep = float(abs(cs[1] - cs[0]))
        c = conditions(sep, h, a1, a2)
        bench.append([tgt[1], tgt[3], h, sep, a1, a2, float(c['ok']), float(oracle(cs, As, c0, a1, a2)), c['P_minus'], c['L1'], c['T_star'], c['T_comp']])
    return dict(S_rows=np.array(rows), S_reg=np.array(reg), S_bench=np.array(bench))


def report_part(raw):
    """Markdown lines of part S from the raw arrays."""
    R, G, B = raw['S_rows'], raw['S_reg'], raw['S_bench']
    L = ['## S. The separated regime, K = 2, n = 1, equal amplitudes: the closed-form condition of Theorem S6\'', '',
         'The sufficient condition (S0)-(S3) of the round note (module docstring of `theory3y_sep.py`) depends on the target only through the two lattice offsets a_1 < a_2 of its features and on the '
         'separation Delta and lattice spacing h.  (S0): Prop. A\'s (L), (C), (N) at r = h/2 + 0.02 (ball arithmetic), and phi = g(D) + g(Delta - D) nonincreasing on [h, Delta/2].', '',
         '**Region of the condition.**  Share of the offset pairs 0 < a_1 < a_2 < h/2 (uniform midpoint grid, the number of pairs is in the last column) for which (S0)-(S3) all hold; `S0` = the box conditions hold at all:', '',
         '| h | Delta | S0 | share of offset pairs covered | pairs |', '|---|---|---|---|---|']
    for r in G:
        L.append('| %g | %g | %d | %s | %d |' % (r[0], r[1], int(r[2]), ('%.4f' % r[3]) if np.isfinite(r[3]) else 'nan', int(r[4])))
    L += ['', '**Against the exact recursion.**  Random targets (offsets uniform, random signs of the features, equal or opposite amplitudes\' signs, random side), the exact recursion of Theorem K '
          '(`exp_theory_t1.limit_path`, hull +- 8 widths): the theorem\'s conclusion = the entries are the nearest kernels, the smaller offset first, right signs, unique entrants, final block = target.  '
          '`violations` = condition true and conclusion false (the theorem would be wrong): ', '',
          '| h | Delta | targets | condition holds | conclusion holds | both | violations |', '|---|---|---|---|---|---|---|']
    for r in R:
        L.append('| %g | %g | %d | %d | %d | %d | **%d** |' % (r[0], r[1], r[2], r[3], r[4], r[5], r[6]))
    L += ['', 'Total violations: **%d** of %d targets (condition true on %d).' % (int(R[:, 6].sum()), int(R[:, 2].sum()), int(R[:, 3].sum())), '']
    if len(B):
        L += ['**The benchmark\'s K = 2, R = 1 targets** (`exp_theory_t1.C_TARGETS`, h = 0.833 widths): offsets, the condition and the recursion:', '',
              '| sep (nominal) | seed | h | Delta | a_1 | a_2 | condition | recursion | P_- | L1 | T_* | min competitor |', '|---|---|---|---|---|---|---|---|---|---|---|---|']
        for r in B:
            L.append('| %g | %d | %.4f | %.3f | %.4f | %.4f | %d | %d | %.4g | %.4g | %.4g | %.4g |' % (r[0], int(r[1]), r[2], r[3], r[4], r[5], int(r[6]), int(r[7]), r[8], r[9], r[10], r[11]))
        L += ['']
    return L
