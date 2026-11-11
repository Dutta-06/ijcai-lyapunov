"""Tests of src/core/baselines.py on the one-dimensional design of Paper 17 / the parametrisation study: 61 samples of [-1, 1], Gaussian kernels of
width 0.12, a 21-kernel start grid, three features (-0.52, 0.07, 0.61) with amplitudes (0.9, -0.6, 0.75)."""
import numpy as np
import pytest

from src.core import baselines as B
from src.core import rich_transit_nd as RN

XS = np.linspace(-1, 1, 61)
K = RN.BumpsND(XS[:, None], 0.12)
GRID = np.linspace(-1, 1, 21)[:, None]
CS = np.array([[-0.52], [0.07], [0.61]])
AS = np.array([0.9, -0.6, 0.75])


def data(sd=0.0, seed=0):
    y = K.synth(CS, AS)
    return y + sd * np.random.default_rng(seed).standard_normal(y.size)


def test_refit_block_recovers_noise_free_features_from_displaced_start():
    y = data()
    c0 = CS + np.array([[0.03], [-0.04], [0.02]])
    w, c = B.refit_block(K, y, np.array([0.7, -0.4, 0.6]), c0)
    assert np.abs(c - CS).max() < 1e-7 and np.abs(w - AS).max() < 1e-7
    assert B.refit_block(K, y, np.zeros(0), np.zeros((0, 1)))[0].size == 0


def test_omp_slide_noise_free_finds_the_three_features_and_the_rules_agree():
    y = data()
    R = B.omp_slide(K, y, GRID, max_units=6)
    blocks = R['blocks']
    assert len(blocks) == 7 and blocks[0]['w'].size == 0 and abs(blocks[0]['rss'] - y @ y) < 1e-12
    b3 = blocks[3]
    srt = np.argsort(b3['c'][:, 0])
    assert np.abs(b3['c'][srt] - CS).max() < 1e-6 and np.abs(b3['w'][srt] - AS).max() < 1e-6 and b3['rss'] < 1e-16
    rss = [b['rss'] for b in blocks]
    assert all(rss[i + 1] <= rss[i] + 1e-12 for i in range(6))
    assert B.select_discrepancy(blocks, K.N, 1e-3) == 3
    assert len(set(R['order'])) == len(R['order'])                          # a grid kernel enters once


def test_omp_slide_in_noise_selects_three_by_bic_discrepancy_and_significance():
    sd = 0.02
    hits = 0
    for seed in range(6):
        y = data(sd, seed)
        R = B.omp_slide(K, y, GRID, max_units=8, s=sd)
        k_bic = B.select_bic(R['blocks'], K.N, 2)
        k_dis = B.select_discrepancy(R['blocks'], K.N, sd)
        k_sig = B.select_significance(R['sig'], 4.5)
        hits += (k_bic == 3) + (k_dis == 3) + (k_sig == 3)
    assert hits >= 15                                                      # of 18: the three rules each right in most of six draws


def test_entry_significance_scales_with_the_signal_and_the_noise_level():
    y = data()
    s1 = B.entry_significance(K, y, CS[:1], 0.1)[0]
    s2 = B.entry_significance(K, y, CS[:1], 0.2)[0]
    assert abs(s1 / s2 - 2) < 1e-12
    # white noise of sd s: the significance of a fixed kernel has unit variance
    rng = np.random.default_rng(1)
    z = np.array([B.entry_significance(K, 0.3 * rng.standard_normal(K.N), CS[:1], 0.3)[0] for _ in range(4000)])        # |signed significance|
    assert abs(np.sqrt(np.mean(z ** 2)) - 1) < 0.05                                                                      # unit second moment


def test_stepwise_bic_rss_decreases_and_the_grid_cannot_fit_off_grid_features_exactly():
    y = data()
    R = B.stepwise_bic(K, y, GRID, max_units=8)
    rss = [b['rss'] for b in R['blocks']]
    assert all(rss[i + 1] <= rss[i] + 1e-12 for i in range(len(rss) - 1))
    assert R['blocks'][3]['rss'] > 1e-6                                    # three grid kernels are not the three off-grid features
    assert len(set(R['order'])) == len(R['order'])


def test_lasso_cd_is_the_soft_threshold_for_an_orthonormal_design_and_satisfies_the_kkt_conditions():
    rng = np.random.default_rng(2)
    Q = np.linalg.qr(rng.standard_normal((30, 5)))[0]
    y = rng.standard_normal(30)
    lam = 0.4
    w = B.lasso_cd(Q, y, lam)
    z = Q.T @ y
    assert np.abs(w - np.sign(z) * np.maximum(np.abs(z) - lam, 0)).max() < 1e-9
    Phi = K.phi(GRID)
    yy = data(0.05, 3)
    lam = 0.5
    w = B.lasso_cd(Phi, yy, lam, iters=20000, tol=1e-13)
    g = Phi.T @ (yy - Phi @ w)
    act = w != 0
    assert np.abs(g[act] - lam * np.sign(w[act])).max() < 1e-6 and np.abs(g[~act]).max() <= lam + 1e-6


def test_lasso_grid_path_starts_at_zero_and_the_support_grows():
    y = data(0.02, 4)
    P = B.lasso_grid(K, y, GRID, n_lam=25)
    assert P[0]['nnz'] == 0 and P[-1]['nnz'] >= 4
    rss = [p['rss'] for p in P]
    assert all(rss[i + 1] <= rss[i] + 1e-9 for i in range(len(rss) - 1))
    assert B.merge_adjacent(np.array([0, 1.0, 2.0, 0, 0, 3.0]), np.arange(6) * 0.1, 0.15) == 2
    assert B.merge_adjacent(np.zeros(4), np.arange(4) * 0.1, 0.15) == 0


def test_train_net_decreases_the_loss_and_a_large_penalty_empties_the_network():
    y = data()
    rng = np.random.default_rng(5)
    w0 = 0.2 * rng.standard_normal(21)
    loss = lambda w, c: 0.5 * np.mean((y - K.synth(c, w)) ** 2)
    l0 = loss(w0, GRID)
    w1, c1 = B.train_net(K, y, GRID, w0, steps=300)
    w2, c2 = B.train_net(K, y, GRID, w0, steps=2000)
    assert loss(w2, c2) < loss(w1, c1) < l0
    w3, c3 = B.train_net(K, y, GRID, w0, lam=1.0, steps=500)
    assert np.abs(w3).max() == 0
    # the proximal step solves the same ℓ1 problem as coordinate descent when the centres stay fixed
    lam = 0.02
    wp, _ = B.train_net(K, y, GRID, w0, lam=lam, steps=30000, train_c=False)
    wc = B.lasso_cd(K.phi(GRID) / np.sqrt(K.N), y / np.sqrt(K.N), lam, iters=50000, tol=1e-14)
    assert np.abs(wp - wc).max() < 5e-3


def test_penalised_and_pruned_networks_fit_and_report_their_size():
    y = data(0.02, 6)
    P = B.penalised_net(K, y, GRID, [1e-4, 1e-2, 1e-1], steps=1500)
    assert [p['lam'] for p in P] == [1e-4, 1e-2, 1e-1] and P[0]['nnz'] >= P[1]['nnz'] >= P[2]['nnz']
    R = B.prune_net(K, y, GRID, steps=3000)
    assert 1 <= R['n_kept'] <= 21 and R['rss'] < 0.05 * (y @ y) and R['kept'].dtype == bool


def test_selection_rules_on_hand_made_sequences():
    N = 100
    blocks = [dict(w=np.zeros(k), rss=r) for k, r in ((0, 50.0), (1, 20.0), (2, 1.05), (3, 1.0), (4, 0.99))]
    assert B.select_bic(blocks, N, 2) == 2
    assert B.select_discrepancy(blocks, N, 0.1) == 2                       # rss/N = 0.0105 ≤ 0.01(1 + 2·0.1414) = 0.01283
    assert B.select_discrepancy(blocks, N, 0.001) == 4                     # nothing within the noise: the last block
    assert B.select_significance([9, 7, 5, 3, 8], 4.5) == 3 and B.select_significance([9, 9], 4.5) == 2 and B.select_significance([1, 9], 4.5) == 0
    assert B.select_gap([1.0, 1.2, 1.3, 9.0, 9.5]) == 3 and B.select_gap([1.0]) == 0 and B.select_gap([1.0, 4.0, 4.1], k_min=2) == 2
    assert abs(B.bic_value(1.0, np.e, 2) - (np.e * np.log(1 / np.e) + 2)) < 1e-12


def test_match_units_pairs_within_the_tolerance_only():
    c_t, w_t = np.array([[0.0], [1.0], [2.0]]), np.array([1.0, 2.0, 3.0])
    n, pos, amp = B.match_units(np.array([[0.01], [1.2], [5.0]]), np.array([1.1, 2.0, 9.0]), c_t, w_t, 0.05)
    assert n == 1 and abs(pos - 0.01) < 1e-12 and abs(amp - 0.1) < 1e-12
    assert B.match_units(np.zeros((0, 1)), np.zeros(0), c_t, w_t, 0.1)[0] == 0
    n, pos, amp = B.match_units(np.array([0.0, 1.0, 2.0]), np.array([1.0, 2.0, 3.0]), c_t, w_t, 0.1)    # 1-D arrays are accepted
    assert n == 3 and pos == 0.0 and amp == 0.0


def test_iss_path_matches_the_repositorys_sliding_iss_on_the_design_and_selects_three_in_noise():
    from experiments import exp_parametrization as EP
    L, Kb, c0 = EP.setup()
    f = Kb.synth(CS, AS)
    ref = RN.sliding_iss_nd(L, Kb, f, c0[:, None], np.ones(len(c0)), polish=True)
    ent = [(t, j) for t, j, kd, c in ref['events'] if kd == 'enter']
    R = B.iss_path(Kb, f, c0, max_entries=3)
    assert R['order'] == [j for t, j in ent] and np.allclose(R['taus'], [t for t, j in ent], rtol=1e-6)
    b3 = R['blocks'][3]
    srt = np.argsort(b3['c'][:, 0])
    assert np.abs(b3['c'][srt] - CS).max() < 1e-6 and np.abs(b3['w'][srt] - AS).max() < 1e-6 and b3['rss'] < 1e-14
    sd = 0.02
    hits = 0
    for seed in range(6):
        y = data(sd, seed)
        R = B.iss_path(K, y, GRID[:, 0], max_entries=8, s=sd)
        hits += (B.select_significance(R['sig'], 4.5) == 3) + (B.select_discrepancy(R['blocks'], K.N, sd) == 3) + (B.select_gap(R['taus']) == 3)
        assert len(R['taus']) == len(R['sig']) == len(R['order']) and np.all(np.diff(R['taus']) >= 0)
    assert hits >= 15


# ---------------------------------------------------------------------------------------------------------------------------
# the BLASSO by sliding Frank–Wolfe and by conic particle gradient descent
# ---------------------------------------------------------------------------------------------------------------------------
FINE = np.linspace(-1, 1, 801)[:, None]


def lam_max(y):
    return float(np.abs(K.corr(GRID, y)[0]).max())


def test_sfw_noise_free_satisfies_the_optimality_conditions_of_the_blasso():
    y = data()
    lam = 0.03 * lam_max(y)
    R = B.sfw(K, y, GRID, lam)
    assert R['converged'] and len(R['w']) == 3
    srt = np.argsort(R['c'][:, 0])
    assert np.abs(R['c'][srt] - CS).max() < 0.01 and np.all(np.sign(R['w'][srt]) == np.sign(AS))
    assert np.all(np.abs(R['w'][srt]) < np.abs(AS)) and np.abs(np.abs(R['w'][srt]) - np.abs(AS)).max() < 0.2        # shrunk by the penalty, not by much
    r = y - K.synth(R['c'], R['w'])
    eta = K.corr(FINE, r)[0] / lam                                                           # the dual certificate on a fine grid
    assert np.abs(eta).max() <= 1.0 + 2e-3                                                    # |eta| <= 1 everywhere
    at = K.corr(R['c'], r)[0] / lam
    assert np.abs(at * np.sign(R['w']) - 1).max() < 2e-3                                      # eta = sign(w) at the support
    assert abs(R['obj'] - (0.5 * (r @ r) / K.N + lam * np.abs(R['w']).sum())) < 1e-15 and abs(R['rss'] - r @ r) < 1e-15


def test_sfw_has_empty_support_above_lambda_max_and_a_growing_one_below():
    y = data()
    G = lam_max(y)
    R = B.sfw(K, y, GRID, 1.01 * G)
    assert len(R['w']) == 0 and R['converged'] and R['cert'] <= 1.0
    sizes = [len(B.sfw(K, y, GRID, f * G)['w']) for f in (1.5, 0.99, 0.5, 0.3, 0.01)]
    assert sizes[0] == 0 and sizes[1] >= 1 and sizes[-1] == 3 and all(sizes[i] <= sizes[i + 1] for i in range(len(sizes) - 1))
    R = B.sfw(K, y, GRID, 0.01 * G, max_units=2)
    assert len(R['w']) == 2 and not R['converged']                                            # stopped by the cap, not by the certificate


def test_cpgd_reaches_the_same_blasso_solution_as_sliding_frank_wolfe():
    y = data()
    lam = 0.1 * lam_max(y)
    R = B.sfw(K, y, GRID, lam)
    Q = B.cpgd(K, y, GRID, lam, steps=3000)
    assert Q['n_kept'] == len(R['w']) == 3
    a, b = np.argsort(R['c'][:, 0]), np.argsort(Q['c'][:, 0])
    assert np.abs(R['c'][a] - Q['c'][b]).max() < 5e-3 and np.abs(R['w'][a] - Q['w'][b]).max() < 2e-2
    assert abs(R['obj'] - Q['obj']) < 2e-3 * R['obj']                                         # the clustering drops the particles' tails below ten start weights


def test_cpgd_particles_stay_non_negative_the_objective_falls_and_a_negative_spike_gets_a_negative_mass():
    c1, a1 = np.array([[0.2]]), np.array([-0.8])
    y = K.synth(c1, a1)
    lam = 0.05 * lam_max(y)
    Q = B.cpgd(K, y, GRID, lam, steps=1500)
    assert np.all(Q['w_part'] >= 0) and Q['w_part'].shape == Q['s_part'].shape == (42,) and Q['x_part'].shape == (42, 1)
    assert Q['obj'] < 0.5 * (y @ y) / K.N and Q['n_kept'] == 1
    assert Q['w'][0] < 0 and abs(Q['c'][0, 0] - 0.2) < 0.02 and 0.7 < -Q['w'][0] < 0.8
    sgn = Q['s_part'][Q['w_part'] > 0.5 * np.abs(Q['w_part']).max()]
    assert np.all(sgn == -1)                                                                   # the dominant particles are the negative ones
    zero = B.cpgd(K, y, GRID, 1.5 * lam_max(y), steps=300)                                       # above lambda_max every particle decays: nothing is kept
    assert zero['n_kept'] == 0 and zero['w'].size == 0


def test_threshold_units_drops_weak_spikes_and_recomputes_the_fit():
    y = data(0.02, 3)
    lam = 0.003 * lam_max(y)
    R = B.sfw(K, y, GRID, lam)
    assert len(R['w']) > 3                                                                    # a small penalty fits noise with weak spikes
    T = B.threshold_units(K, y, R, lam)
    assert len(T['w']) < len(R['w']) and np.abs(T['w']).min() >= 0.1 * np.abs(T['w']).max() and np.abs(R['w']).max() == np.abs(T['w']).max()
    gone = np.abs(R['w']) < 0.1 * np.abs(R['w']).max()
    assert gone.any() and len(T['w']) == int((~gone).sum()) and np.allclose(np.sort(T['w']), np.sort(R['w'][~gone]))        # exactly the weak spikes are dropped
    r = y - K.synth(T['c'], T['w'])
    assert abs(T['rss'] - r @ r) < 1e-15 and abs(T['obj'] - (0.5 * (r @ r) / K.N + lam * np.abs(T['w']).sum())) < 1e-15
    assert T['rss'] >= R['rss'] - 1e-15 and T['converged'] == R['converged']
    E = B.threshold_units(K, y, dict(w=np.zeros(0), c=np.zeros((0, 1)), rss=0.0, obj=0.0), lam)                # the empty solution
    assert E['w'].size == 0 and abs(E['rss'] - y @ y) < 1e-15


# ---------------------------------------------------------------------------------------------------------------------------
# added after the first mutation run of baselines.py
# ---------------------------------------------------------------------------------------------------------------------------
def test_stepwise_bic_picks_the_best_grid_kernel_first_never_twice_and_finds_on_grid_features():
    on = np.array([3, 10, 16])
    y = K.synth(GRID[on], np.array([1.0, -0.5, 0.8]))
    R = B.stepwise_bic(K, y, GRID, max_units=8)                            # past the three features all candidates tie at rounding level: a kernel must still not be chosen twice
    blocks = R['blocks']
    Phi = K.phi(GRID)
    one_fit = [float(np.sum((y - Phi[:, j] * (Phi[:, j] @ y) / (Phi[:, j] @ Phi[:, j])) ** 2)) for j in range(len(GRID))]
    assert R['order'][0] == int(np.argmin(one_fit)) and abs(blocks[1]['rss'] - min(one_fit)) < 1e-12      # by the smallest residual, not the largest
    assert len(set(R['order'])) == len(R['order'])                                                       # a kernel is chosen once
    assert sorted(R['order'][:3]) == sorted(on.tolist()) and blocks[3]['rss'] < 1e-20
    rss = [b['rss'] for b in blocks]
    assert all(rss[i + 1] <= rss[i] + 1e-12 for i in range(len(rss) - 1)) and blocks[2]['rss'] > 1e-3          # noise-free: BIC is degenerate at rss ~ 1e-30, the residual is what is tested


def test_iss_path_ties_enter_together_with_equal_times():
    y = K.synth(np.array([[-0.45], [0.45]]), np.array([0.8, 0.8]))         # mirror-symmetric target on a symmetric grid: two duals reach the wall at the same time
    R = B.iss_path(K, y, GRID[:, 0], max_entries=2)
    assert len(R['taus']) == 2 and abs(R['taus'][0] - R['taus'][1]) < 1e-9 * R['taus'][0]
    assert len(R['blocks'][1]['w']) == 2                                   # both are in the first block: a tie is one event


def test_sfw_tolerance_stops_on_the_right_side_and_reports_the_certificate_in_units_of_lambda():
    y = data()
    G = lam_max(y)
    R = B.sfw(K, y, GRID, 0.7 * G, tol=0.5)                                # the first candidate's certificate is 1/0.7 = 1.43 <= 1.5: the tolerance stops it
    assert len(R['w']) == 0 and R['converged']
    assert len(B.sfw(K, y, GRID, 0.7 * G, tol=0.1)['w']) >= 1
    R = B.sfw(K, y, GRID, 1.25 * G)
    assert len(R['w']) == 0 and abs(R['cert'] - 0.8) < 1e-9                # max|<r phi>| / lambda with r = y


def test_sfw_slide_keeps_signs_and_bounds_the_centres():
    y = K.synth(np.array([[0.07]]), np.array([-0.8]))
    lam = 0.01 * lam_max(y)
    w, c = B._slide_l1(K, y, np.array([0.3]), np.array([[0.07]]), lam, np.array([-1.0]), np.array([1.0]))
    assert abs(w[0]) < 1e-9                                                # a positive coefficient cannot fit a negative feature: it goes to 0, it does not change sign
    y2 = K.synth(np.array([[0.5]]), np.array([0.8]))
    w, c = B._slide_l1(K, y2, np.array([0.6]), np.array([[0.0]]), 1e-4, np.array([-1.0]), np.array([0.2]))
    assert abs(c[0, 0] - 0.2) < 1e-6 and w[0] > 0                          # the centre stops at its bound, short of the feature at 0.5
    w, c = B._slide_l1(K, y2, np.array([0.6]), np.array([[0.0]]), 1e-4, np.array([-1.0]), np.array([1.0]))
    assert abs(c[0, 0] - 0.5) < 1e-4                                       # unbounded in practice: it reaches the feature


def test_cluster_joins_close_particles_of_one_sign_only():
    x = np.array([[0.0], [0.05], [0.08], [0.5]])
    lab = B._cluster(x, np.array([1.0, -1.0, 1.0, 1.0]), 0.1)
    assert lab[0] == lab[2] and lab[1] != lab[0] and lab[3] != lab[0] and len(set(lab)) == 3
    lab = B._cluster(x, np.ones(4), 0.1)
    assert len(set(lab)) == 2 and lab[0] == lab[1] == lab[2]


def test_cpgd_spikes_are_the_total_mass_at_the_mass_weighted_position_and_weak_spikes_are_dropped():
    y = K.synth(np.array([[-0.5], [0.5]]), np.array([1.0, 0.04]))
    lam = 0.01 * lam_max(y)
    Q = B.cpgd(K, y, GRID, lam, steps=1500)
    assert Q['n_kept'] == 1 and abs(Q['c'][0, 0] + 0.5) < 0.02            # the spike of 0.04 is below a tenth of the largest and is dropped ...
    Q2 = B.cpgd(K, y, GRID, lam, steps=1500, thr=0.01)
    assert Q2['n_kept'] == 2                                               # ... and kept when the threshold allows it
    yy = data()
    w0 = 0.01
    P = B.cpgd(K, yy, GRID, 0.1 * lam_max(yy), steps=14, w_init=w0)       # unconverged: the particles of a spike have different weights and positions
    live = P['w_part'] > 10 * w0
    xl, sl, wl = P['x_part'][live], P['s_part'][live], P['w_part'][live]
    lab = B._cluster(xl, sl, 1.5 * 0.1)
    assert live.sum() == 8 and len(set(lab)) == 3 and P['n_kept'] == 3     # eight live particles, three spikes (not bridged by the dormant particles between them)
    for w, c in zip(P['w'], P['c']):
        hit = [k for k in np.unique(lab) if abs(np.sum((sl * wl)[lab == k]) - w) < 1e-12]
        assert len(hit) == 1
        m = lab == hit[0]
        assert abs(c[0] - np.sum(xl[m, 0] * wl[m]) / wl[m].sum()) < 1e-12                                   # the mass-weighted mean, not the plain one
    assert len({round(float(x), 6) for x in P['c'][:, 0]}) == 3
    far = B.cpgd(K, yy, GRID, 0.1 * lam_max(yy), steps=14, w_init=w0, link=0.02)
    assert far['n_kept'] == 8                                              # a smaller link radius cuts the neighbours of a spike apart


def test_cpgd_objective_never_rises_even_from_a_heavy_start_and_recovers_its_step():
    y = data()
    lam = 0.05 * lam_max(y)
    w0 = 1.0                                                               # a heavy start: the first multiplicative steps overshoot unless they are rejected and halved
    d = len(GRID)

    def fit(q):
        return 0.5 * float(np.sum((y - K.synth(q['x_part'], q['s_part'] * q['w_part'])) ** 2)) / K.N + lam * float(q['w_part'].sum())
    f_start = 0.5 * float(np.sum((y - K.synth(np.vstack([GRID, GRID]), np.r_[np.ones(d), -np.ones(d)] * w0)) ** 2)) / K.N + lam * 2 * d * w0
    fs = [fit(B.cpgd(K, y, GRID, lam, steps=n, w_init=w0)) for n in (1, 2, 3, 5, 10, 30, 100)]
    assert fs[0] <= f_start + 1e-12 and all(fs[i + 1] <= fs[i] + 1e-12 for i in range(len(fs) - 1))      # never rises (without the acceptance test it explodes to 1e25)
    assert fs[-1] < 0.05                                                   # the step recovers after the rejections (a step that never grows back is still at 0.10 after 100 steps)
    R = B.sfw(K, y, GRID, lam)
    assert fit(B.cpgd(K, y, GRID, lam, steps=3000, w_init=w0)) < R['obj'] * 1.05 + 1e-6     # 3000 steps reach the BLASSO optimum


def test_sfw_drops_a_unit_whose_slide_ends_at_zero_weight():
    """A noisy two-feature target found by a random search (the 369th draw of default_rng(1), noise sd 0.1) at a penalty of 1 % of the largest correlation: eight units
    enter and the slide of the block drives one coefficient to exactly zero; the returned support must not keep it."""
    rng = np.random.default_rng(1)
    for _ in range(369):
        k = rng.integers(2, 6)
        c = np.sort(rng.uniform(-0.8, 0.8, k))
        a = rng.uniform(0.2, 1.0, k) * rng.choice([-1, 1], k)
        sd = rng.choice([0.0, 0.02, 0.1])
        y = K.synth(c[:, None], a) + sd * rng.standard_normal(K.N)
    assert k == 2 and sd == 0.1
    R = B.sfw(K, y, GRID, 0.01 * lam_max(y))
    assert len(R['w']) == 7 and np.abs(R['w']).min() > 1e-3
