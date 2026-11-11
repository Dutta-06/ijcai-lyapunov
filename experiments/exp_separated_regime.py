"""
Exact recovery of separated features: the path against its decoupled limit (the perturbative theorem of the separated regime)
=============================================================================================================================

The claim. A target of K Gaussian features of the kernels' width, separated by Delta, has a *decoupled limit* (Gamma(Delta) -> 0, Gamma(d) =
exp(-d^2/(4 sigma^2)) the kernels' overlap): K one-feature paths running in parallel. In that limit feature k is entered by its best starting
kernel j_k (the one of largest correlation a_j = <A_k phi(c*_k) phi_j>) at the hitting time T_k = rho_{j_k} / |a_{j_k}|, the entered unit slides onto the
feature exactly, the residual near the feature is zero and the other kernels near it stop at |u_j| / rho_j = T_k / T_j < 1. If the decoupled limit is
generic (each feature's best kernel unique, the T_k distinct) the path at a large enough Delta has the same K entries in the same order, no exit,
and ends at the exact fit: the blocks' critical points are nondegenerate (implicit function theorem), the duals are piecewise linear in tau, the
entries are transversal, and a coupled entry converges (supplement, Lemma "inner problem" (iii)); the error of every entry time is O(Gamma(Delta)).

What this script measures (1-D, sigma = 0.12, 21 starting kernels of spacing 0.1 on [-1, 1], 61 samples, equal masses; the path is
`rich_transit_nd.sliding_iss_nd`, the same as the 1154-target record of `exp_transit_certificate.py`):
  A. K = 2 targets, separation 2..8 sigma, 60 per separation, amplitudes +-U(0.4, 1), positions uniform: how often the path is exact with exactly K
     units and no exit; how often its entries are the decoupled prediction's (set and order); where they are, the relative error of the entry times against
     the first-order cross-talk Gamma(Delta - h/2) of a kernel half a lattice cell from its feature (the best kernel of a feature is at most h/2 away, so the
     error of a first entry is (A_l/A_k) Gamma(Delta - d_1)/Gamma(d_1) with d_1 <= h/2: it is NOT O(Gamma(Delta)) for a fixed offset d_1 > 0); the agreement as a
     function of the decoupled limit's own margins (mu0 = the smallest over features of 1 - T_(1)/T_(2), gap = the smallest relative gap between two features'
     hitting times). The theorem needs both to exceed a constant times its overlap parameter, whose constants put Delta_0 above every separation tested here
     (referee round IS2); the experiment tests the STRUCTURE (the decoupled prediction and the size of its error), and the gate min(mu0, gap) > Gamma(Delta - h/2) is
     an empirical, lattice-specific law, not a theorem.
  B. K = 3 targets with minimal separation 3, 4 and 5 sigma, 60 per separation (rejection sampling as in `exp_transit_certificate.targets`), the same columns;
     the separation of a target is its actual minimal separation.

  C. (added after the full-draft review) K = 2 at Delta = 10, 12 and 14 sigma, 60 per separation, in a wider window (centres in [-0.95, 0.95], the lattice spans [-1, 1]): the separations at which the
     perturbative theorem's own constants put Delta_0 (about 12 sigma), where the overlap Gamma(Delta) is below 1e-5; run by `--run-extra`, which adds the keys 'C|<Delta>' to the raw of A and B
     (their keys are not touched).

Usage: python experiments/exp_separated_regime.py --run [--workers N] | --run-extra [--workers N] | --from-raw          -> results/chaos/separated_regime.md (+ _raw.npz)
"""
import argparse
import os
import sys
import warnings

os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np

HERE = os.path.abspath(os.path.dirname(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)

OUT = os.environ.get('SEP_OUT', os.path.join(REPO, 'results', 'chaos'))
RAW, MD = os.path.join(OUT, 'separated_regime_raw.npz'), os.path.join(OUT, 'separated_regime.md')
SEPS_A = (2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0)
SEPS_B = (3.0, 4.0, 5.0)
SEPS_C = (10.0, 12.0, 14.0)
WIN_C = 0.95                      # part C: the centres lie in [-WIN_C, WIN_C]
N_PER = 60
SEED_A, SEED_B, SEED_C = 1000, 2000, 3000
RMS_TOL = 1e-6                    # an exact fit: the final rms residual
H1_MARGIN = 0.05                  # the margin at which Paper 17 counts a path's hypotheses (H1) as holding on random targets
_SETUP = {}


def setup():
    """The 1-D samples, the starting kernels and the kernels at the samples (cached per process)."""
    if not _SETUP:
        import exp_transit_mass as TM
        L, K, c0 = TM.setup()
        _SETUP.update(L=L, K=K, c0=np.asarray(c0, float), Phi0=K.phi(np.asarray(c0, float)))
    return _SETUP['L'], _SETUP['K'], _SETUP['c0'], _SETUP['Phi0']


def gamma_of(sep, sig):
    """The kernels' overlap Gamma(d) = exp(-d^2 / (4 sigma^2))."""
    return float(np.exp(-sep ** 2 / (4 * sig ** 2)))


def decoupled(Phi0, bump, cs, As):
    """Per feature the decoupled hitting times of all starting kernels, T_j = 1 / |a_j| with a_j the mean over the samples of A phi(c*) phi_j (equal walls rho = 1).
    Returns [(best kernel j_k, T_k, margin_k)] with margin_k = 1 - T_k / T_second (the second-best kernel's hitting time)."""
    out = []
    N = Phi0.shape[0]
    for c, A in zip(cs, As):
        a = np.abs(Phi0.T @ (A * bump(float(c)))) / N
        order = np.argsort(-a)
        j, j2 = int(order[0]), int(order[1])
        out.append((j, 1.0 / a[j], 1.0 - a[j2] / a[j]))
    return out


def decoupled_summary(dec):
    """(j_k in the order of the hitting times, the sorted T_k, mu0 = the smallest margin, gap = the smallest relative gap between two features' T_k)."""
    srt = sorted(dec, key=lambda x: x[1])
    Ts = np.array([t for _, t, _ in srt])
    gap = float(np.min(np.diff(Ts)) / Ts[-1]) if len(Ts) > 1 else 1.0
    return [j for j, _, _ in srt], Ts, float(min(m for _, _, m in dec)), gap


def one(cs, As, sig_scale):
    """The path on one target and its comparison with the decoupled prediction."""
    import exp_mirror_transit as X
    import exp_transit_certificate as XC
    L, K, c0, Phi0 = setup()
    warnings.filterwarnings('ignore')
    f = X.bumps_target(K, cs, As)
    dec = decoupled(Phi0, lambda c: X.bumps_target(K, [c], [1.0]), cs, As)
    order_dec, Ts, mu0, gap = decoupled_summary(dec)
    R = XC.path_1d(L, K, f, c0)
    ev = [(float(t), int(j)) for t, j, k, c in R['events'] if k == 'enter']
    n_exit = sum(1 for t, j, k, c in R['events'] if k == 'exit')
    act = sorted(int(j) for j in np.flatnonzero(R['act']))
    kf = len(cs)
    rms = float(R['hist'][-1][3])
    cf = np.asarray(R['c'], float).reshape(len(c0), -1)[act, 0] if act else np.zeros(0)
    wf = np.asarray(R['w'], float).reshape(-1)[act] if act else np.zeros(0)
    srt = np.argsort(cs)
    if len(act) == kf:
        o = np.argsort(cf)
        cerr, aerr = float(np.abs(cf[o] - np.asarray(cs)[srt]).max()), float(np.abs(wf[o] - np.asarray(As)[srt]).max())
    else:
        cerr = aerr = np.nan
    same_seq = len(ev) == kf and [j for t, j in ev] == order_dec
    terr = float(np.max([abs(t - T) / T for (t, j), T in zip(ev, Ts)])) if same_seq else np.nan      # only where the entries ARE the decoupled prediction's: entry k against its own T_k
    sepmin = float(np.min(np.abs(np.subtract.outer(cs, cs))[np.triu_indices(kf, 1)]) / sig_scale)
    return np.array([float(rms < RMS_TOL and len(act) == kf and n_exit == 0), float(len(act)), float(n_exit),
                     float(sorted(order_dec) == act), float(same_seq), terr, cerr, aerr, mu0, gap, sepmin])


COLS = ('exact', 'units', 'exits', 'same_set', 'same_order', 'terr', 'cerr', 'aerr', 'mu0', 'gap', 'sepmin')


def target_a(seed, sep, sig):
    rng = np.random.default_rng(seed)
    c1 = rng.uniform(-0.7, 0.7 - sep * sig)
    return np.array([c1, c1 + sep * sig]), rng.choice([-1.0, 1.0], 2) * rng.uniform(0.4, 1.0, 2)


def target_b(seed, dmin, sig):
    rng = np.random.default_rng(seed)
    while True:
        cs = rng.uniform(-0.7, 0.7, 3)
        if min(abs(cs[i] - cs[j]) for i in range(3) for j in range(i)) >= dmin * sig:
            return cs, rng.choice([-1.0, 1.0], 3) * rng.uniform(0.4, 1.0, 3)


def target_c(seed, sep, sig):
    """Two features at exactly sep sigma in the wider window [-WIN_C, WIN_C]."""
    rng = np.random.default_rng(seed)
    c1 = rng.uniform(-WIN_C, WIN_C - sep * sig)
    return np.array([c1, c1 + sep * sig]), rng.choice([-1.0, 1.0], 2) * rng.uniform(0.4, 1.0, 2)


def job(args):
    part, sep, i = args
    _, K, _, _ = setup()
    cs, As = {'A': lambda: target_a(SEED_A + i, sep, K.sig), 'B': lambda: target_b(SEED_B + i, sep, K.sig), 'C': lambda: target_c(SEED_C + i, sep, K.sig)}[part]()
    return one(cs, As, K.sig)


def run(workers):
    import multiprocessing as mp
    _, K, _, _ = setup()
    jobs = [('A', s, i) for s in SEPS_A for i in range(N_PER)] + [('B', s, i) for s in SEPS_B for i in range(N_PER)]
    if workers > 1:
        with mp.Pool(workers) as pool:
            rows = pool.map(job, jobs, chunksize=4)
    else:
        rows = [job(j) for j in jobs]
    raw = {}
    for (part, s) in [('A', s) for s in SEPS_A] + [('B', s) for s in SEPS_B]:
        raw['%s|%g' % (part, s)] = np.array([r for (p, s_, i), r in zip(jobs, rows) if p == part and s_ == s])
    raw['sig'] = np.array([K.sig])
    raw['h'] = np.array([float(setup()[2][1] - setup()[2][0])])           # the lattice spacing of the starting kernels
    os.makedirs(OUT, exist_ok=True)
    np.savez_compressed(RAW, **raw)
    return raw


def run_extra(workers):
    """Part C only: the rows 'C|<Delta>' are added to the raw of A and B, whose keys stay as they are."""
    import multiprocessing as mp
    raw = load_raw()
    jobs = [('C', s, i) for s in SEPS_C for i in range(N_PER)]
    if workers > 1:
        with mp.Pool(workers) as pool:
            rows = pool.map(job, jobs, chunksize=4)
    else:
        rows = [job(j) for j in jobs]
    for s in SEPS_C:
        raw['C|%g' % s] = np.array([r for (p, s_, i), r in zip(jobs, rows) if s_ == s])
    np.savez_compressed(RAW, **raw)
    return raw


def load_raw():
    with np.load(RAW) as z:
        return {k: z[k] for k in z.files}


def slope(x, y):
    """Least-squares slope of log y on log x over the positive finite pairs."""
    m = (x > 0) & (y > 0) & np.isfinite(x) & np.isfinite(y)
    return float(np.polyfit(np.log(x[m]), np.log(y[m]), 1)[0]) if m.sum() >= 2 else np.nan


def table(header, rows):
    out = ['| ' + ' | '.join(header) + ' |', '|' + '|'.join(['---'] * len(header)) + '|']
    return '\n'.join(out + ['| ' + ' | '.join(str(c) for c in r) + ' |' for r in rows])


def report(raw):
    sig, h = float(raw['sig'][0]), float(raw['h'][0])
    ci = {c: i for i, c in enumerate(COLS)}
    gam_of = lambda sep_sig: np.exp(-np.maximum(sep_sig * sig - h / 2, 0.0) ** 2 / (4 * sig ** 2))          # Gamma(Delta - h/2): the cross-talk of a kernel half a lattice cell from its feature
    L = ['# Exact recovery of separated features: the path against its decoupled limit', '',
         'Generated by `experiments/exp_separated_regime.py` from `separated_regime_raw.npz`; every number below is printed by the script (protocol: its docstring). '
         '1-D, sigma = %g, 21 starting kernels of spacing h = %g, equal masses, `sliding_iss_nd`; %d targets per row. Gamma(d) = exp(-d^2 / (4 sigma^2)) is the overlap of two '
         'kernels at the distance d; Delta is the separation of the features (for K = 3 the target\'s actual minimal separation; the row\'s Delta is the lower bound it was drawn with). '
         '"Decoupled" is the prediction by K one-feature paths running in parallel (the best starting kernel of each feature, entered at T_k = 1 / |a_j|). '
         'Exact: final rms residual below %g, exactly K units, no exit. "Gate": min(mu0, gap) > Gamma(Delta - h/2), where mu0 is the smallest margin of a feature\'s best kernel over its '
         'runner-up and gap the smallest relative gap between two features\' hitting times of the decoupled limit: an empirical, lattice-specific law, not the theorem\'s condition '
         '(the theorem\'s constants put its Delta_0 above these separations, referee round IS2). The entry-time error is the largest over the entries of |tau_k - T_k| / T_k, '
         'taken over the targets whose entries ARE the decoupled prediction\'s (same kernels, same order); it is divided by Gamma(Delta - h/2) of the same target.' % (sig, h, N_PER, RMS_TOL), '']
    tot_gate = tot_gate_ok = tot_non = tot_non_ok = 0
    ratios_far = []
    slopes = {}
    parts = [('A', 'K = 2, separation exactly Delta (positions uniform, amplitudes +-U(0.4, 1))', SEPS_A, 'Delta/sigma'),
             ('B', 'K = 3, minimal separation at least Delta', SEPS_B, 'min Delta/sigma')]
    if all('C|%g' % s_ in raw for s_ in SEPS_C):
        parts.append(('C', 'K = 2, separation exactly Delta, wider window (centres in [-%g, %g]; the separations of the perturbative theorem\'s own constants)' % (WIN_C, WIN_C), SEPS_C, 'Delta/sigma'))
    for part, title, seps, unit in parts:
        rows, gam_row, err_max = [], [], []
        for s in seps:
            d = raw['%s|%g' % (part, s)]
            n = len(d)
            gt = gam_of(d[:, ci['sepmin']])
            te = d[:, ci['terr']]
            cnt = lambda col: int(round(d[:, ci[col]].sum()))
            gate = np.minimum(d[:, ci['mu0']], d[:, ci['gap']]) > gt
            ss = d[:, ci['same_set']] > 0
            ok = np.isfinite(te)
            ratio = te[ok] / gt[ok]
            rows.append(['%g' % s, '%.2e' % gam_of(np.array(s)), '%d / %d' % (cnt('exact'), n), '%d' % cnt('exits'), '%d / %d' % (cnt('same_set'), n), '%d / %d' % (cnt('same_order'), n),
                         '%d / %d' % (int((ss & gate).sum()), int(gate.sum())), '%d / %d' % (int((ss & ~gate).sum()), int((~gate).sum())), '%d' % int(ok.sum()),
                         ('%.2f' % np.median(ratio)) if (ok.any() and part != 'C') else '-', ('%.2f' % ratio.max()) if (ok.any() and part != 'C') else '-', '%.1e' % np.nanmax(d[:, ci['cerr']]),
                         '%.1e' % np.nanmax(d[:, ci['aerr']]), '%d / %d' % (int((np.minimum(d[:, ci['mu0']], d[:, ci['gap']]) >= H1_MARGIN).sum()), n)])
            gam_row.append(float(gam_of(np.array(s))))
            err_max.append(float(te[ok].max()) if ok.any() else np.nan)
            if s >= 5.0 and ok.any() and part != 'C':
                ratios_far.append(ratio.max())
            tot_gate += int(gate.sum())
            tot_gate_ok += int((ss & gate).sum())
            tot_non += int((~gate).sum())
            tot_non_ok += int((ss & ~gate).sum())
        L += ['## %s' % title, '',
              table([unit, 'Gamma(Delta - h/2)', 'exact, K units', 'exits (all targets)', 'same entry set as decoupled', 'same order', 'same set where gate', 'same set elsewhere',
                     'targets with the decoupled entries', 'entry-time error / Gamma(Delta - h/2), median', 'max', 'max centre error', 'max coefficient error',
                     'min(mu0, gap) >= %g' % H1_MARGIN], rows), '']
        if part == 'C':           # at these separations the overlap is below the solver's precision: the entry-time error is not first order in Gamma any more
            terr_all = np.concatenate([raw['C|%g' % s_][:, ci['terr']] for s_ in seps])
            L += ['In C the overlap Gamma(Delta - h/2) is below %.0e: the largest relative entry-time error over the targets with the decoupled entries is %.1e (the solver\'s precision, not first order '
                  'in Gamma), so the ratio columns are not given.' % (max(gam_of(np.array(s_)) for s_ in seps), np.nanmax(terr_all)), '']
        if part == 'A':           # the separation is exact only in A (in B it is a lower bound), so the scaling is read there
            big = np.array([s >= 4.0 for s in seps])
            gam_row, err_max = np.array(gam_row), np.array(err_max)
            L += ['The log-log slope of the largest entry-time error against Gamma(Delta - h/2) over the separations of at least 4 sigma: %.2f (first-order cross-talk: 1).' % slope(
                gam_row[big], err_max[big]), '']
    allexact = sum(int(round(raw[k][:, ci['exact']].sum())) for k in raw if k[:2] in ('A|', 'B|', 'C|'))
    alln = sum(len(raw[k]) for k in raw if k[:2] in ('A|', 'B|', 'C|'))
    L += ['Over all %d targets of %s the path is exact with exactly K units and no exit on %d.' % (alln, 'A, B and C' if len(parts) == 3 else 'A and B', allexact),
          'Where the decoupled limit\'s margins exceed the gate (min(mu0, gap) > Gamma(Delta - h/2): %d targets) the path\'s entry set is the decoupled prediction\'s on %d; elsewhere (%d targets) on %d.' % (
              tot_gate, tot_gate_ok, tot_non, tot_non_ok),
          'At the separations of at least 5 sigma the largest entry-time error / Gamma(Delta - h/2) over the rows is %.2f.' % (max(ratios_far) if ratios_far else float('nan')), '']
    return '\n'.join(L) + '\n'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', action='store_true')
    ap.add_argument('--run-extra', action='store_true')
    ap.add_argument('--from-raw', action='store_true')
    ap.add_argument('--workers', type=int, default=int(os.environ.get('SEP_WORKERS', '2')))
    a = ap.parse_args()
    if a.run:
        run(a.workers)
    if a.run_extra:
        run_extra(a.workers)
    if a.run or a.run_extra or a.from_raw:
        with open(MD, 'w') as f:
            f.write(report(load_raw()))
        print('wrote', MD)
    else:
        ap.error('give --run, --run-extra or --from-raw')


if __name__ == '__main__':
    main()
