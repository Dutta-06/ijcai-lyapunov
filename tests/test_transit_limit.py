"""Theorem D‴.53: the sliding limit — the inner problem is a gradient flow of the fit (an isolated entrant converges to its feature), and
the loop reaches the sliding path with the predicted first-order correction (the entries' pulls and shifts, the creep)."""
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'experiments'))
import exp_transit_limit as TL
from src.core import rich_transit as RT


def test_isolated_entrant_inner_problem_descends_the_fit_to_its_feature():
    """One bump A·φ(c*) and one entrant from ξ = 200: dξ/ds = −⟨rφ⟩, dc/ds = Λ₀Gm w⟨r∂_cφ⟩, w = 𝒲(ξ) — the fit E = ½⟨r²⟩ never increases
    (Lemma D‴.53a: a gradient flow in the mirror coordinate's metric) and the entrant ends on the feature, w = A, c = c*."""
    L, K, c0 = TL.setup()
    A, cs, ce = 0.8, 0.051, 0.1
    f = A * K.phi(np.array([cs]))[:, 0]
    lam0 = RT.transit_scale(L, 1.0)
    nx = K.xs.size

    def wc(z):
        xi, c = z
        w = float(RT.inner_coef(L, xi))
        return w, f - w * K.phi(np.array([c]))[:, 0]

    def rhs(s, z):
        w, r = wc(z)
        c = z[1]
        return [-(K.phi(np.array([c]))[:, 0] @ r) / nx, lam0 * L.G * L.m * w * (K.dphi(np.array([c]))[:, 0] @ r) / nx]
    ss = np.linspace(0, 6000, 601)
    sol = solve_ivp(rhs, (0, ss[-1]), [200.0, ce], method='LSODA', rtol=1e-10, atol=1e-12, t_eval=ss)
    E = np.array([0.5 * np.mean(wc(z)[1] ** 2) for z in sol.y.T])
    assert np.all(np.diff(E) <= 1e-12 * E[0])
    assert E[-1] < 1e-14 * E[0]
    w_end, _ = wc(sol.y[:, -1])
    assert abs(w_end - A) < 1e-6 and abs(sol.y[1, -1] - cs) < 1e-6


def run(n, eps):
    """D‴.51's target n: its constants (`b_job`: converged pulls, s_θ, (H3)) and the loop against the path at the mass ε (`a_job`)."""
    cs, As = TL.TM.random_targets(TL.SEED_A, TL.N_POOL)[n]
    C = TL.b_job((n, cs, As))
    saved = TL.EPS_A
    try:
        TL.EPS_A = (eps,)
        r = TL.a_job(C)
    finally:
        TL.EPS_A = saved
    return cs, As, C, r['rows'][0]


def test_loop_reaches_the_path_with_its_first_order_correction():
    """A three-bump target whose path leaves margins ≥ 0.15 (D‴.51's target 39), the loop at ε = 10⁻¹⁰: exactly the path's kernels; the
    never-entering kernels on the core's profile; their centres displaced by the creep's εC_j; their duals off the path by ε^{1/2}·(the
    entries' pulls and shifts) + εD_j to 2·10⁻³ε^{1/2}; the entries' crossings at τ_k + ε^{1/2}(s_θ + σ_k) to 0.02ε^{1/2}."""
    cs, As, C, row = run(39, 1e-10)
    assert C['ok'] and C['mu0'] > 0.15
    eps = row['eps']
    assert row['n_form'] == 3
    assert row['block'] < 1e-10
    disp, creep = np.array(row['disp']), np.array(row['creep'])
    big = np.abs(creep) > 1e-3 * np.abs(creep).max()
    assert np.abs(disp[big] / creep[big] - 1).max() < 5e-3
    res = np.array([d[2] - d[3] - d[4] for d in row['dev']])
    assert np.abs(res).max() < 2e-3
    alone = np.array([d[2] - d[5] for d in row['dev']])
    assert np.abs(alone).max() > 0.03                                    # the shifts' term is needed
    assert max(abs(a - b) for _, a, b in row['cross']) < 0.02
    # the never-entering kernels' coefficients: the core's profile at their path margins, a(εb)^{3/2}/(m(2c)^{1/2})·cos(πμ/2)/sin³(πμ/2)
    L, K, c0, f, R, P = TL.path_of(cs, As)
    never = [j for j in range(c0.size) if j not in P['ever']]
    t_end = 1.3 * C['entries'][-1][0]
    k = len(P['taus']) - 1
    U = P['U'][k] + (t_end - P['taus'][k]) * (K.phi(P['centres'][k]).T @ P['res'][k] / K.xs.size)
    m = 1 - np.abs(U[never])
    prof = L.A * eps ** 1.5 / (L.m * np.sqrt(2 * L.c)) * np.cos(np.pi * m / 2) / np.sin(np.pi * m / 2) ** 3
    assert abs(row['w_never'] / prof.max() - 1) < 1e-2


def test_sign_per_entry_and_the_shifts_on_a_sliding_target():
    """D‴.51's target 41 at ε = 10⁻⁹ (referee round 41: its second entrant slides, its third arrives 4.8ε^{1/2} late): the duals to
    1.7·10⁻³ε^{1/2} and the crossings to 0.045ε^{1/2} (about 3× the observed 5.7·10⁻⁴ and 0.015) with converged predictions — the
    extrapolation in xi0 spread below 10⁻⁴ —; (H3): every entry's inner problem ends at the path's S_k; the block within 2·10⁻¹¹."""
    cs, As, C, row = run(41, 1e-9)
    assert max(C['s_spread']) < 1e-4 and max(C['p_spread']) < 1e-5
    assert max(C['h3w']) < 1e-9 and max(C['h3c']) < 1e-9 and all(C['stationary'])
    assert max(abs(s) for s in C['shifts']) > 1.0                       # an entry's shift matters here
    assert row['n_form'] == 3 and row['block'] < 2e-11
    res = np.array([d[2] - d[3] - d[4] for d in row['dev']])
    assert np.abs(res).max() < 1.7e-3
    assert max(abs(a - b) for _, a, b in row['cross']) < 0.045


def test_entry_shift_carries_the_entrants_mass():
    """σ_k = −s_k·dev/(ρ_{j_k}v_k) (Theorem D‴.53(c), referee round 140): a dual s(u + ε^{1/2}dev) growing at the absolute rate ρv reaches
    its wall ρ at τ_k + ε^{1/2}σ_k, for either sign; at ρ = 1 it is B's −s·dev/v."""
    for s, rho, v, dev in ((1.0, 0.5, 0.07, 0.3), (-1.0, 1.7, 0.02, -0.11), (1.0, 1.0, 0.05, 0.2)):
        tau_k, eps = 10.0, 1e-8
        u = lambda t: s * rho - s * rho * v * (tau_k - t)                # the path's dual: s·u reaches ρ at τ_k at the rate ρv
        t_wall = tau_k - eps ** 0.5 * s * dev / (rho * v)                 # s(u(t) + ε^{1/2}dev) = ρ
        assert np.isclose(s * (u(t_wall) + eps ** 0.5 * dev), rho, rtol=1e-14, atol=0)
        assert np.isclose(TL.entry_shift(dev, s, v, rho), (t_wall - tau_k) / eps ** 0.5, rtol=1e-6, atol=0)
    assert TL.entry_shift(0.2, 1.0, 0.05) == -0.2 / 0.05


def test_unequal_masses_the_shift_needs_one_over_rho():
    """Part U (referee round 140): D‴.51's target 41 with masses log-uniform in [0.5, 2]; the path exact with simple entries; the loop at
    ε = 10⁻¹⁰ crosses θ at s_θ + σ_k to 0.03 (observed 0.011) with σ_k = −s_k·dev/(ρ_{j_k}v_k), while without the factor 1/ρ_{j_k} it misses
    by more than 0.1 (observed 0.19); the never-entering duals at ε^{1/2}dev to 0.02 (observed 0.0056)."""
    saved = TL.EPS_U
    try:
        TL.EPS_U = (1e-10,)
        U = TL.u_job()
    finally:
        TL.EPS_U = saved
    assert U['simple'] and U['mu0'] > 0.05 and len(U['entries']) == 3
    rho = np.array([e['rho'] for e in U['entries']])
    assert np.ptp(rho) > 0.5                                             # the entrants' masses differ
    e3 = U['entries'][2]
    assert abs(e3['sigma'] - e3['sigma_wrong']) > 0.15
    assert np.isclose(e3['sigma'], e3['sigma_wrong'] / e3['rho'], rtol=1e-12, atol=0)
    row = U['rows'][0]
    assert row['n_form'] == 3
    cr = [c for _, c in row['cross']]
    assert max(abs(c - e['scross'] - e['sigma']) for c, e in zip(cr, U['entries'])) < 0.03
    assert max(abs(c - e['scross'] - e['sigma_wrong']) for c, e in zip(cr, U['entries'])) > 0.1
    assert max(abs(a - b) for _, _, a, b in row['du']) < 0.02
    # the stored Part U reproduced (mutation testing after round 140): the masses, the path's margin, every entry's constants and the
    # loop's row at ε = 10⁻¹⁰, so that a change in `u_masses`, `u_job`'s ρ, v_k, pulls or σ_k shows here and not only in the thresholds
    S = _raw()['U']
    assert tuple(r['eps'] for r in S['rows']) == TL.EPS_U and S['n'] == TL.U_TARGET
    assert np.array_equal(TL.u_masses(len(S['b'])), np.array(S['b'])) and U['b'] == list(S['b'])
    assert U['simple'] == S['simple'] and np.isclose(U['mu0'], S['mu0'], rtol=1e-10, atol=0)
    for e, f in zip(U['entries'], S['entries']):
        assert (e['k'], e['j'], e['stationary']) == (f['k'], f['j'], f['stationary'])
        for key in ('tau', 'rho', 'v', 'theta', 'scross', 'dev', 'sigma', 'sigma_wrong'):
            assert np.isclose(e[key], f[key], rtol=1e-9, atol=1e-12), key
        assert np.isclose(e['s_spread'], f['s_spread'], rtol=1e-6, atol=1e-12)
        assert np.allclose([e['devs_after'][j] for j in sorted(e['devs_after'])], [f['devs_after'][j] for j in sorted(f['devs_after'])],
                           rtol=1e-9, atol=1e-12)
    srow = S['rows'][-1]
    assert srow['eps'] == 1e-10 and row['n_form'] == srow['n_form']
    assert [j for j, _ in row['cross']] == [j for j, _ in srow['cross']]
    assert np.allclose([c for _, c in row['cross']], [c for _, c in srow['cross']], rtol=1e-9, atol=0)
    assert [d[:2] for d in row['du']] == [tuple(d[:2]) for d in srow['du']]
    assert np.allclose([d[2:] for d in row['du']], [d[2:] for d in srow['du']], rtol=1e-8, atol=1e-10)


RAW = os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos', 'transit_limit_raw.npz')


def _raw():
    return TL.E.load(RAW)


def test_path_duals_rho_pins_the_entrants_at_their_walls():
    """`path_duals_rho` (Part U): from its entry on, each entrant's dual sits at its own wall sign·ρ_j (not ±1), from the first entry
    on; the never-entering kernels' duals are `RT.path_duals`' (mutation testing after round 140)."""
    L, K, c0 = TL.setup()
    rho = TL.u_masses(c0.size) ** -0.5
    cs, As = TL.TM.random_targets(TL.SEED_A, TL.N_POOL)[TL.U_TARGET]
    f = TL.X.bumps_target(K, cs, As)
    R = RT.sliding_iss(L, K, f, c0, rho, polish=True)
    P0, P = RT.path_duals(K, f, R), TL.path_duals_rho(K, f, R, rho)
    T = P['taus']
    first = {}
    for t, j, kind, c in R['events']:
        first.setdefault(j, min(k for k in range(len(T)) if abs(T[k] - t) < 1e-12))
    assert len(first) == 3 and min(first.values()) == 1
    for k in range(len(T)):
        for j in range(c0.size):
            if j in first and k >= first[j]:
                assert abs(abs(P['U'][k][j]) - rho[j]) < 1e-15 and np.sign(P['U'][k][j]) == np.sign(P0['U'][k][j]), (k, j)
            else:
                assert P['U'][k][j] == P0['U'][k][j], (k, j)
    assert np.ptp([rho[j] for j in first]) > 0.5


def test_part_W_reproduces_the_stored_windows():
    """Part W (referee round 140) against its raw: the targets are A's (`chosen`, in order, `part_W`'s jobs B's rows), λ_k from B; target 4
    (its second interval's window skipped at ε ≥ 10⁻⁷) rerun at every ε gives the stored windows and blocks (mutation testing after round
    140: the start τ_k + (C_W/λ_k)ε^{1/2}log(1/ε), the skip at τ_{k+1} − ε^{1/4}, the block over the coefficients and the centres)."""
    import multiprocessing
    raw = _raw()
    B, W = list(raw['B']), list(raw['W'])
    byB = {b['n']: b for b in B}
    assert [t['n'] for t in W] == TL.chosen(B)
    assert all(list(t['lam']) == list(byB[t['n']]['lam']) for t in W)
    assert all(tuple(r['eps'] for r in t['rows']) == TL.EPS_A for t in W)

    class FakePool:
        def __init__(self, *a, **k):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

        def map(self, fn, jobs, chunksize=1):
            return [(fn.__name__, j['n']) for j in jobs]
    saved = multiprocessing.Pool
    try:
        multiprocessing.Pool = FakePool
        assert TL.part_W(B) == [('w_job', t['n']) for t in W]
    finally:
        multiprocessing.Pool = saved
    t = {x['n']: x for x in W}[4]
    w = TL.w_job(byB[4])
    assert w['n'] == 4 and w['lam'] == list(t['lam'])
    assert [q for q, _ in t['rows'][0]['wins']] == [0, 2] and [q for q, _ in t['rows'][2]['wins']] == [0, 1, 2]
    for r, s in zip(w['rows'], t['rows']):
        assert r['eps'] == s['eps'] and [q for q, _ in r['wins']] == [q for q, _ in s['wins']]
        assert np.allclose([d for _, d in r['wins']], [d for _, d in s['wins']], rtol=1e-6, atol=0), r['eps']
    # planted: λ_2 such that at ε = 10⁻⁶ the second window starts 0.01 before τ_3 — inside τ_3 − ε^{1/4} = τ_3 − 0.0316, so skipped
    tk = [e[0] for e in byB[4]['entries']]
    eps = 1e-6
    C = dict(byB[4], lam=[byB[4]['lam'][0], TL.C_W * np.sqrt(eps) * np.log(1 / eps) / (tk[2] - tk[1] - 0.01), byB[4]['lam'][2]])
    saved = TL.EPS_A
    try:
        TL.EPS_A = (eps,)
        w = TL.w_job(C)
    finally:
        TL.EPS_A = saved
    assert [q for q, _ in w['rows'][0]['wins']] == [0, 2]



def test_report_W_takes_the_transient_over_the_intervals_present():
    """W's column 'A: C_B·λ_k' is C_B times the least λ_k over the intervals whose window (A's start τ_k + C_Bε^{1/2}log(1/ε), its end
    τ_{k+1} − ε^{1/4}) is non-empty at every 10⁻⁹ ≤ ε ≤ 10⁻⁶. On the stored targets λ_k decreases, so the last interval decides; planted
    (target 39): the first interval the slowest, 0.424 after its entry the second entry, so that its window is empty at ε = 10⁻⁶
    (0.414 + 10⁻¹·⁵ > 0.424) though not if it ended at τ_{k+1} − ε^{1/2} — the column stays C_B·λ_3 (mutation testing after round 140)."""
    import copy
    raw = copy.deepcopy(_raw())
    B, A, W, U = (list(raw[k]) if k != 'U' else raw[k] for k in 'BAWU')
    t = next(x for x in W if x['n'] == 39)
    b = next(x for x in B if x['n'] == 39)
    row = lambda md: next(ln for ln in md.splitlines() if ln.startswith('| 39 | '))
    assert row(TL.report(B, A, W, U)).endswith(f" | {TL.C_B * t['lam'][2]:.3f} |")
    t['lam'] = [0.04] + list(t['lam'][1:])
    e = b['entries']
    b['entries'] = [e[0], (e[0][0] + 0.424,) + tuple(e[1][1:])] + list(e[2:])
    assert row(TL.report(B, A, W, U)).endswith(f" | {TL.C_B * t['lam'][2]:.3f} |")


if __name__ == '__main__':
    for name, fn in list(globals().items()):
        if name.startswith('test_'):
            fn(); print(name, 'ok')
