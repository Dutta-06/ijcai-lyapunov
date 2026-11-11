"""Proposition D‴.56: the transit under measurement noise — the sliding path's linear response (`rich_transit.path_response`), the window
and the stopped fit (`experiments/exp_noisy_certificate.py`)."""
import os
import sys
import warnings

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'experiments')))
warnings.filterwarnings('ignore')

from src.core import rich_transit as RT
import exp_noisy_certificate as NC
import exp_transit_mass as TM
import exp_mirror_transit as X

CS, AS = np.array([-0.52, 0.02, 0.61]), np.array([0.8, -0.35, 0.9])


def test_block_response_at_an_exact_fit_is_the_projection_off_the_tangent_space():
    L, K, c0 = TM.setup()
    f = X.bumps_target(K, CS, AS)
    Rm = RT.block_response(K, f, AS, CS)
    J = np.c_[K.phi(CS), K.dphi(CS) * AS[None, :]]
    assert np.abs(Rm @ Rm - Rm).max() < 1e-10 and np.abs(Rm - Rm.T).max() < 1e-10
    assert np.abs(Rm @ J).max() < 1e-10
    assert abs(np.trace(Rm) - (K.xs.size - 2 * CS.size)) < 1e-8                 # rank N − 2K


def test_sliding_iss_stops_at_max_entries_before_relaxing():
    L, K, c0 = TM.setup()
    f = X.bumps_target(K, CS, AS)
    R = RT.sliding_iss(L, K, f, c0, np.ones(c0.size), polish=True, max_entries=2)
    assert [e[2] for e in R['events']] == ['enter', 'enter'] and len(R['relaxed']) == 1 and len(R['hist']) == 2


def test_path_response_is_the_derivative_of_the_noisy_path():
    """Central differences of the noisy path's entry times and of its duals at the K-th entry against T·ξ and U·ξ."""
    L, K, c0 = TM.setup()
    f = X.bumps_target(K, CS, AS)
    R0 = RT.sliding_iss(L, K, f, c0, np.ones(c0.size), polish=True)
    P = RT.path_response(K, f, R0, c0)
    xi = np.random.default_rng(3).standard_normal(K.xs.size)
    h, nK = 1e-5, CS.size
    kw = dict(polish=True, max_events=nK + 3, tol=0.0, max_entries=nK + 1)
    Rp = RT.sliding_iss(L, K, f + h * xi, c0, np.ones(c0.size), **kw)
    Rm = RT.sliding_iss(L, K, f - h * xi, c0, np.ones(c0.size), **kw)
    tp, tm = (np.array([e[0] for e in R_['events'][:nK]]) for R_ in (Rp, Rm))
    assert np.abs((tp - tm) / (2 * h) - P['T'] @ xi).max() < 1e-5 * np.abs(P['T'] @ xi).max()
    off = [j for j in range(c0.size) if j not in P['entrants']]
    up, um = (RT.path_duals(K, f + s_ * h * xi, R_)['U'][nK] for s_, R_ in ((1, Rp), (-1, Rm)))
    fd = ((up - um) / (2 * h))[off]
    assert np.abs(fd - (P['U'] @ xi)[off]).max() < 1e-5 * np.abs((P['U'] @ xi)[off]).max()
    assert np.abs(P['u_all'][-1] - P['u']).max() == 0 and np.abs(P['U_all'][-1] - P['U']).max() == 0


def test_window_next_kernel_and_oracle_at_small_noise():
    """s = 0.001 on a separated target: the path enters exactly E, the first noise kernel enters after the first-order window W (to 1%:
    the relative error is O(s), the residual R_Kξ + O(|ξ|²)) and is the predicted one, and the fit stopped at the K-th entry is the oracle
    — its error is the noise's projection on the tangent space (to 1%)."""
    L, K, c0, f, cert, R0, P = NC.noiseless(CS, AS)
    nK, nx = CS.size, K.xs.size
    off = [j for j in range(c0.size) if j not in cert['e']]
    for seed in (1, 2):
        xi = 0.001 * np.random.default_rng(seed).standard_normal(nx)
        R = NC.noisy(L, K, c0, f, xi, nK)
        ent = [e for e in R['events'] if e[2] == 'enter']
        assert sorted(e[1] for e in ent[:nK]) == sorted(cert['e']) and len(ent) == nK + 1
        v = K.phi(c0).T @ (P['Rlast'] @ xi) / nx
        u = P['u'] + P['U'] @ xi
        W = {j: (1 - np.sign(v[j]) * u[j]) / abs(v[j]) for j in off}
        jn = min(W, key=W.get)
        assert ent[nK][1] == jn and abs((ent[nK][0] - ent[nK - 1][0]) / W[jn] - 1) < 1e-2
        hK = R['hist'][nK]
        err = np.sqrt(np.mean((K.phi(hK[1]) @ hK[2] - f) ** 2))
        err_T = np.sqrt(np.mean(((np.eye(nx) - P['Rlast']) @ xi) ** 2))
        assert abs(err / err_T - 1) < 1e-2


def _runaway_draw():
    """Referee round 45's runaway path: B's target 16 at s = 0.1, the old Part E stream's draw 3 — past the first noise kernel the edge
    kernel 20 leaves the samples' window and the fit has no minimiser."""
    cs, As = NC.targets()[16]
    L, K, c0, f, cert, R0, P = NC.noiseless(cs, As)
    rng = np.random.default_rng(3016)
    xs = [s * rng.standard_normal(K.xs.size) for s in (0.01, 0.03, 0.1, 0.3) for _ in range(5)]
    return L, K, c0, f + xs[13]


def test_sliding_iss_box_flags_a_runaway_kernel_and_the_budget_stops_the_path():
    L, K, c0, y = _runaway_draw()
    kw = dict(polish=True, max_events=14, tol=0.0, max_entries=7)
    R = RT.sliding_iss(L, K, y, c0, np.ones(c0.size), **kw)
    assert not R['runaway'] and np.abs(R['hist'][-1][2]).max() > 1e3            # without the box: a "fit" with |w| ~ 10⁶
    assert len(R['ev_iter']) == len(R['events'])
    ex = [i for i, e in zip(R['ev_iter'], R['events']) if e[2] == 'exit']
    assert ex and all(i == ex[0] for i in ex)                                   # the exits burst in one iteration
    Rb = RT.sliding_iss(L, K, y, c0, np.ones(c0.size), box=NC.box_of(K), w_scale=float(np.abs(y).max()), **kw)
    lo, hi = NC.box_of(K)
    a = Rb['w'] != 0
    assert Rb['runaway'] and (Rb['c'][a].max() > hi or Rb['c'][a].min() < lo)
    assert len(Rb['hist']) == len(Rb['relaxed']) and not Rb['relaxed'][-1]       # the runaway block is no fit of hist
    assert all(np.abs(h[2]).max() < 1.0 for h in Rb['hist'])
    assert [e[2] for e in Rb['events']] == ['enter'] * 5
    Rw = RT.sliding_iss(L, K, y, c0, np.ones(c0.size), w_scale=float(np.abs(y).max()), **kw)
    assert all(e[2] == 'enter' for e in Rw['events'])                          # on the data's scale the features' kernels stay big …
    assert max(np.abs(h[2]).max() for h in Rw['hist']) > 1e3 and all(Rw['relaxed'])   # … but only the box catches the runaway
    Rc = RT.sliding_iss(L, K, y, c0, np.ones(c0.size), max_rhs=1000, **kw)
    assert Rc['out_of_budget'] and Rc['nfev'] == 1001 and not Rc['relaxed'][-1] and len(Rc['events']) < len(R['events'])


def test_the_duals_at_every_entry_are_the_derivative_of_the_noisy_path():
    """Central differences of the noisy path's duals at each of its first K entries against u_all/U_all for every kernel not yet entered
    there — off the path and the path's later entrants (the order criterion's rows)."""
    L, K, c0, f, cert, R0, P = NC.noiseless(CS, AS)
    nK, nx = CS.size, K.xs.size
    xi = np.random.default_rng(5).standard_normal(nx)
    h = 1e-5
    Ds = []
    for sg in (1, -1):
        R = NC.noisy(L, K, c0, f, sg * h * xi, nK)
        ent, exits, es, eo = NC.classify(R, P['entrants'])
        assert eo
        D = RT.path_duals(K, f + sg * h * xi, R)['U']
        Ds.append(np.array([D[NC.fit_after(ent, k)] for k in range(1, nK + 1)]))
    fd = (Ds[0] - Ds[1]) / (2 * h)
    rs, ro = NC.crit_rows(P, c0.size)
    lin = np.array([P['U_all'][k][j] @ xi for k, j in ro])
    num = np.array([fd[k, j] for k, j in ro])
    assert np.abs(num - lin).max() < 1e-5 * np.abs(lin).max()
    later = [(k, j) for k, j in ro if j in P['entrants']]
    assert later and set(rs) | set(later) == set(ro)
    assert all(abs(P['u_all'][k][j]) < 1 for k, j in ro)


def test_criteria_rows_radius_and_union_bound():
    L, K, c0, f, cert, R0, P = NC.noiseless(CS, AS)
    ent = list(P['entrants'])
    rs, ro = NC.crit_rows(P, c0.size)
    assert all(j not in ent for k, j in rs) and all(j not in ent[:k + 1] for k, j in ro)
    assert {k for k, j in rs} == {k for k, j in ro} == set(range(len(ent)))     # every entry, the K-th included
    assert {(0, ent[1]), (0, ent[2]), (1, ent[2])} <= set(ro) and (1, ent[1]) not in ro
    u, U = NC.crit_arrays(P, ro)
    r1 = NC.radius(u, U)
    i = int(np.argmin((1 - np.abs(u)) / np.linalg.norm(U, axis=1)))
    xi = r1 * U[i] / np.linalg.norm(U[i]) * np.sign(u[i] if u[i] else 1.0)     # the worst direction: the row reaches its wall
    assert abs(abs(u[i] + U[i] @ xi) - 1) < 1e-12
    assert np.abs(u + U @ (0.999 * xi)).max() < 1
    from scipy.stats import norm
    s = 0.05
    assert abs(NC.union_bound(u, U, s) - np.sum(2 * norm.sf((1 - np.abs(u)) / (s * np.linalg.norm(U, axis=1))))) < 1e-15


def test_classify_and_the_stopping_rules():
    def R_(evs):
        return dict(events=[(t, j, k, 0.0) for t, j, k, i in evs], ev_iter=[i for t, j, k, i in evs])
    order = [5, 3, 9]
    ent, ex, es, eo = NC.classify(R_([(1, 5, 'enter', 0), (2, 3, 'enter', 1), (3, 9, 'enter', 2), (9, 7, 'enter', 3)]), order)
    assert es and eo and ex == 0 and NC.fit_after(ent, 3) == 3
    ent, ex, es, eo = NC.classify(R_([(1, 5, 'enter', 0), (2, 9, 'enter', 1), (3, 3, 'enter', 2)]), order)
    assert es and not eo                                                        # a swap: the set, not the order
    ent, ex, es, eo = NC.classify(R_([(1, 5, 'enter', 0), (2, 3, 'enter', 1), (2, 3, 'exit', 1), (2, 3, 'enter', 2),
                                      (3, 9, 'enter', 3)]), order)
    assert not es and ex == 1                                                   # a re-entry among the first K entries
    ent, ex, es, eo = NC.classify(R_([(1, 5, 'enter', 0), (2, 3, 'enter', 1), (3, 9, 'enter', 2), (3, 3, 'exit', 2),
                                      (9, 7, 'enter', 3)]), order)
    assert not es and not eo and ex == 1                                        # the first K entries E, an exit before the (K + 1)-th
    ent, ex, es, eo = NC.classify(R_([(1, 5, 'enter', 0), (2, 3, 'enter', 1), (3, 9, 'enter', 2), (9, 7, 'enter', 3),
                                      (9, 7, 'exit', 3)]), order)
    assert es and eo and ex == 0                                                # a later exit does not count
    assert NC.disc_stop([0.5, 0.2, 0.012, 0.0105, 0.009], 0.01) == 4
    assert NC.disc_stop([0.5, 0.2], 0.01) is None
    taus = [1.0, 2.0, 2.5, 30.0, 33.0]
    assert NC.gap_stop(taus, 5) == 3 and NC.gap_stop(taus, 3) == 1 and NC.gap_stop(taus, 1) is None


def test_b_draw_on_a_seeded_draw():
    """One draw of Part B against an independent computation of its classification, criteria, window, stop and fit."""
    L, K, c0, f, cert, R0, P = NC.noiseless(CS, AS)
    nK, nx = CS.size, K.xs.size
    order = list(P['entrants'])
    off = [j for j in range(c0.size) if j not in order]
    rs, ro = NC.crit_rows(P, c0.size)
    rows = (NC.crit_arrays(P, rs), NC.crit_arrays(P, ro))
    s = 0.03
    xi = NC.noise(2000, 0, 2, 0, s, nx)
    d = NC.b_draw(L, K, c0, f, P, order, off, rows, s, xi)
    R = NC.noisy(L, K, c0, f, xi, nK)
    ent = [e for e in R['events'] if e[2] == 'enter']
    assert d['exact'] == (sorted(e[1] for e in ent[:nK]) == sorted(order)) and d['exact_order'] == ([e[1] for e in ent[:nK]] == order)
    assert d['exact_order']
    crit = max(abs(P['u_all'][k][j] + P['U_all'][k][j] @ xi) for k in range(nK) for j in off)
    crit_ord = max(abs(P['u_all'][k][j] + P['U_all'][k][j] @ xi) for k in range(nK) for j in range(c0.size) if j not in order[:k + 1])
    assert abs(d['crit'] - crit) < 1e-14 and abs(d['crit_ord'] - crit_ord) < 1e-14 and crit_ord >= crit
    v = K.phi(c0).T @ (P['Rlast'] @ xi) / nx
    u = P['u'] + P['U'] @ xi
    W = min((1 - np.sign(v[j]) * u[j]) / abs(v[j]) for j in off)
    assert abs(d['W_pred'] - W) < 1e-12 * W and abs(d['W_obs'] - (ent[nK][0] - ent[nK - 1][0])) < 1e-12
    taus = [e[0] for e in ent]
    assert np.allclose(d['gaps'], [(taus[i + 1] - taus[i]) / taus[i] for i in range(len(taus) - 1)], rtol=0, atol=1e-15)
    it = [i for e, i in zip(R['events'], R['ev_iter']) if e[2] == 'enter']
    rms = [R['hist'][it[k] + 1][3] for k in range(nK)]
    assert d['k_disc'] == next((k + 1 for k in range(nK) if rms[k] <= 1.1 * s), -1) and d['rms_K'] == rms[-1]
    hK = R['hist'][it[nK - 1] + 1]
    assert abs(d['err'] - np.sqrt(np.mean((K.phi(hK[1]) @ hK[2] - f) ** 2))) < 1e-15
    assert d['lin_err'] < 10 * s ** 2 and abs(d['crit_path'] - crit) < d['lin_err'] + 1e-12


def _fake_path(s, taus, fit_rms, fit_err, oracle, exact=True):
    n = len(fit_rms)
    return dict(s=s, exact=exact, exact_order=exact, exits_before=0, exits=0, runaway=False, budget=False, nfev=0, stiff=0,
                n_ent=len(taus), taus=taus, ent_j=list(range(len(taus))), ent_iter=list(range(len(taus))), relaxed=[True] * n,
                fit_rms=fit_rms, fit_err=fit_err, fit_out=[0.0] * n, fit_maxw=[1.0] * n, oracle=oracle, x_end=0.5, maxw_end=1.0)


def test_report_E_on_synthetic_paths():
    """Part E's table from two hand-built paths: one where both rules stop at the K-th fit, one where a later noise gap beats W/τ_K over
    K + 4 entries and the discrepancy rule stops one fit late."""
    pa = _fake_path(0.01, [1, 1.5, 2, 40, 45, 50, 52], [0.3, 0.1, 0.0105, 0.0102, 0.0100, 0.0098], [0.3, 0.1, 0.004, 0.0045, 0.005, 0.0055], 0.004)
    pb = _fake_path(0.01, [1, 1.5, 2, 10, 100, 101, 102], [0.3, 0.1, 0.0115, 0.0105, 0.0100, 0.0098], [0.3, 0.1, 0.004, 0.006, 0.007, 0.008], 0.004)
    E_ = [dict(fam='base', q=0, K=3, more=4, order=[0, 1, 2], amps=[1, 1, 1], levels=[dict(s=0.01, draws=[pa, pb])])]
    md = NC.report(E_=E_)
    row = next(l for l in md.splitlines() if l.startswith('| base | 0.01 |'))
    assert row == '| base | 0.01 | 2 (0) | 0, 0 | 2 | 2/2, 1/2 | 2, 1 of 2 | 1, 1, 0, 0 | 1.250 [1.5] | 1.250 [1.5] | 0, 0 |'
    row2 = [l for l in md.splitlines() if l.startswith('| base | 0.01 |')][1]
    assert row2 == '| base | 0.01 | 0.5 [0.5] | 4.56 [9] | 11.5 [4] |'                  # later gaps: max 0.125 and 9
