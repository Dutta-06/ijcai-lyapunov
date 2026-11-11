"""Proposition D‴.51: the transit at a finite mass — the inner profile, the entrant's overshoot, the near-tie inner problem, the path's
duals and the predicted mass."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'experiments'))
from src.core import mirror_free as MF
from src.core import rich_transit as RT

L = MF.FreeLaw(m=12.0, G=1.0, gamma=0.1, a_max=1.0)


def grid():
    return RT.Bumps(np.linspace(-1, 1, 61), 0.12), np.linspace(-1, 1, 21)


def test_inner_profile_is_the_small_mass_limit_of_the_mirror_map():
    """A kernel of mass ε at α has its dual (λ₁ − η(α))/(ε^{1/2}Λ) below the wall; as ε → 0 this is inner_xi(w(α)) — the critical
    curve's η_c, not the core's power law (which 𝒲 falls short of by 6.7% at ξ = 15 and 65% at ξ = 3)."""
    al = np.array([0.05, 0.2, 0.5, 1.0, 3.3])
    for eps, tol in ((1e-7, 5e-3), (1e-9, 5e-4)):
        q = L.Q(np.zeros(1), np.array([eps]))
        Lam = RT.transit_scale(L, eps)
        lam1 = Lam * (1 + 0.0)
        xi_lib = np.array([(lam1 - float(L.eta(np.array([a]), q)[0])) / (np.sqrt(eps) * Lam) for a in al])
        w = np.array([float(L.w(np.array([a]), q)[0]) for a in al])
        assert np.allclose(RT.inner_xi(L, w), xi_lib, rtol=0, atol=tol)
    xi = np.array([-9.0, -1.0, 0.0, 2.0, 40.0, 500.0])
    assert np.allclose(RT.inner_xi(L, RT.inner_coef(L, xi)), xi, rtol=0, atol=1e-9)
    Kp3 = (8 / np.pi ** 3) * L.A / (L.m * np.sqrt(2 * L.c))
    assert abs(RT.inner_coef(L, 500.0) * 500.0 ** 3 / Kp3 - 1) < 2e-3                      # the core's power law far below the wall
    assert abs(RT.inner_coef(L, 15.0) * 15.0 ** 3 / Kp3 - 1) > 0.05


def test_overshoot_on_the_exponential_branch():
    """The overshoot −inner_xi(A) = η_c(α_A)/(Λ₀ρ); for large A the branch is exponential, w ≈ (ca/(2m))e^{2α/a}, so the overshoot grows like
    (a/(4GcΛ₀))·log A, and two kernels whose duals differ by Δξ carry coefficients in the ratio e^{−κΔξ}, κ = πρ(2c)^{1/2}."""
    A = np.array([0.3, 0.8, 2.0, 20.0, 200.0])
    ov = RT.overshoot(L, A)
    assert np.allclose(ov, -RT.inner_xi(L, A))
    lam0 = RT.transit_scale(L, 1.0)
    slope = (ov[-1] - ov[-2]) / np.log(10.0)
    assert abs(slope / (L.A / (4 * L.G * L.c * lam0)) - 1) < 1e-3
    kap = np.pi * np.sqrt(2 * L.c)
    w1, w2 = RT.inner_coef(L, -6.0), RT.inner_coef(L, -6.0 + 1.5)
    assert abs(np.log(w1 / w2) / (kap * 1.5) - 1) < 2e-3


def test_near_tie_inner_problem():
    """At an exact tie (x = 0) the inner problem splits the bump equally; at small x the split is the Gibbs one 1/(1 + e^{κx}); far behind,
    the runner-up keeps the frozen coefficient K'³ρ⁻³(x − Δ)⁻³, Δ the pull a frozen runner-up feels (the overshoot less I∞ > 0)."""
    K, c0 = grid()
    kap = np.pi * np.sqrt(2 * L.c)
    Kp3 = (8 / np.pi ** 3) * L.A / (L.m * np.sqrt(2 * L.c))
    r0 = RT.near_tie_inner(L, K, 0.8, 0.05, 0.1, 0.0, 0.0)
    assert abs(r0['share'] - 0.5) < 1e-6 and abs(r0['c_e'] - 0.05) < 1e-6 and abs(r0['c_r'] - 0.05) < 1e-6
    r1 = RT.near_tie_inner(L, K, 0.8, 0.05, 0.1, 0.0, 1.0)
    assert abs(r1['share'] * (1 + np.exp(kap)) - 1) < 0.03
    fr = RT.near_tie_inner(L, K, 0.8, 0.05, 0.1, 0.0, 200.0, both=False)
    delta = 200.0 - fr['xi_r']
    ov = float(RT.overshoot(L, 0.8))
    assert abs(fr['xi_e'] + ov) < 1e-4 and 0.5 < ov - delta < 2.0                         # the entrant ends at its overshoot; I∞ ∈ (0.5, 2)
    r3 = RT.near_tie_inner(L, K, 0.8, 0.05, 0.1, 0.0, 300.0)
    assert abs(r3['w_r'] * (300.0 - delta) ** 3 / Kp3 - 1) < 0.01
    rn = RT.near_tie_inner(L, K, -0.8, 0.05, 0.1, 0.0, 3.0)
    rp = RT.near_tie_inner(L, K, 0.8, 0.05, 0.1, 0.0, 3.0)
    assert abs(rn['share'] / rp['share'] - 1) < 1e-6


def test_near_tie_loop_against_the_inner_problem():
    """The loop at ε = 10⁻⁹ on a bump 10⁻³ off the midpoint (x ≈ 110): the runner-up's share is the inner problem's at the same x to 2%;
    at ε = 10⁻⁶ (x ≈ 3.5) both kernels form, the runner-up with the Gibbs share to 10%."""
    import exp_mirror_transit as X
    K, c0 = grid()
    f = X.bumps_target(K, (0.051,), (0.8,))
    g = K.phi(c0).T @ f / K.xs.size
    mu = 1 - abs(g[10]) / abs(g[11])
    kap = np.pi * np.sqrt(2 * L.c)
    for eps, tol in ((1e-9, 0.02), (1e-6, 0.10)):
        taus, W, cc, Lam = X.run_moving(L, K, f, c0, eps, 1.32 / abs(g[11]), 'paradigm', n_t=111)
        share = abs(W[-1][10]) / (abs(W[-1][10]) + abs(W[-1][11]))
        x = mu / np.sqrt(eps)
        inner = RT.near_tie_inner(L, K, 0.8, 0.05, 0.1, 0.0, x)
        assert abs(share / inner['share'] - 1) < tol
        if eps == 1e-6:
            assert abs(share * (1 + np.exp(kap * x)) - 1) < 0.1


def three_bumps(n, seed=21):
    rng = np.random.default_rng(seed)
    for k in range(n + 1):
        while True:
            cs = np.sort(rng.uniform(-0.8, 0.8, 3))
            if np.diff(cs).min() >= 0.3:
                break
        As = rng.choice([-1, 1], 3) * rng.uniform(0.4, 1.0, 3)
    return cs, As


def test_path_duals_reach_the_walls_at_the_entries():
    """The reconstructed duals: every entrant's dual reaches its wall exactly at its entry (it is pinned there by the path's own entry
    rule), and the kernels that never enter stay inside their boxes all along."""
    import exp_mirror_transit as X
    K, c0 = grid()
    cs, As = three_bumps(10)
    f = X.bumps_target(K, cs, As)
    R = RT.sliding_iss(L, K, f, c0, np.ones(21))
    P = RT.path_duals(K, f, R)
    taus = np.array(P['taus'])
    nx = K.xs.size
    for t, j, kind, c in R['events']:
        if kind != 'enter':
            continue
        k = int(np.argmin(np.abs(taus - t)))
        u_before = P['U'][k - 1][j] + (taus[k] - taus[k - 1]) * (K.phi(P['centres'][k - 1]).T @ P['res'][k - 1] / nx)[j]
        assert abs(abs(u_before) - 1) < 1e-8
    never = [j for j in range(21) if j not in P['ever']]
    assert np.abs(P['U'][:, never]).max() < 1


def test_recovery_threshold_on_the_lost_near_tie():
    """D‴.48's target 10: the path leaves kernel 15 0.0065 inside its wall beside an unfitted bump — the predicted mass is the creep's,
    ≈ 5·10⁻⁹ (the loop forms the path's kernels at 2·10⁻⁹ and not at 10⁻⁸, `exp_transit_mass` C); the creep is monotone in ε."""
    import exp_mirror_transit as X
    K, c0 = grid()
    cs, As = three_bumps(10)
    f = X.bumps_target(K, cs, As)
    R = RT.sliding_iss(L, K, f, c0, np.ones(21))
    p = RT.recovery_threshold(L, K, f, c0, R, 1e-2 * np.abs(As).max())
    assert 3e-9 < p['eps_creep'] < 8e-9 and min(p['creep'], key=p['creep'].get) == 15
    assert p['eps_c'] == p['eps_creep'] and p['eps_entry'] > 1e-7
    P = RT.path_duals(K, f, R)
    assert RT.creep_reaches(L, K, P, 15, 2e-8, c0[15]) and not RT.creep_reaches(L, K, P, 15, 1e-9, c0[15])


def test_entry_pull_is_the_matched_frozen_pull():
    """`entry_pull` on a path's first entry (one bump 10⁻³ off the midpoint) gives the frozen runner-up's matched pull of the two-kernel
    inner problem — its outer margin at the entrant's arrival less its final ξ — whatever the start, which it converges in like 1/xi0."""
    import exp_mirror_transit as X
    K, c0 = grid()
    f = X.bumps_target(K, (0.051,), (0.8,))
    R = RT.sliding_iss(L, K, f, c0, np.ones(21), polish=True)
    P = RT.path_duals(K, f, R)
    k = int(np.argmin(np.abs(np.array(P['taus']) - R['events'][0][0])))
    pulls = []
    for xi0 in (100.0, 200.0):
        ep = RT.entry_pull(L, K, f, R, P, k, 11, [10], xi0=xi0)['pull'][10]
        fr = RT.near_tie_inner(L, K, 0.8, 0.051, 0.1, 0.0, 200.0, both=False, xi0=xi0)
        assert abs(ep - fr['pull']) < 1e-3
        pulls.append(ep)
    assert 0 < pulls[0] - pulls[1] < 0.2                                 # the far-field creep's 1/xi0 tail


def test_rearranging_entry_pull_is_measured():
    """Part E's target 327: the third entry sends an active kernel on to another feature, so κ_ij = 444 and the near-tie estimate alone
    predicts a lost near tie up to ε ≈ 7·10⁻¹⁰; the measured pull puts the entry threshold above 10⁻⁷, and the loop at 2·10⁻⁹ forms exactly
    the path's kernels."""
    import exp_mirror_transit as X
    import exp_transit_mass as T
    Lx, K, c0 = T.setup()
    cs, As = T.random_targets(T.SEED_E, T.N_E)[327]
    f = X.bumps_target(K, cs, As)
    R = RT.sliding_iss(L, K, f, c0, np.ones(21), polish=True)
    th = X.TH_ACT * np.abs(As).max()
    p0 = RT.recovery_threshold(L, K, f, c0, R, th, n_creep=0, measured_pull=False)
    p1 = RT.recovery_threshold(L, K, f, c0, R, th, n_creep=0)
    assert p0['entry'][3] > 100 and p0['eps_entry'] < 1e-9 < 1e-7 < p1['eps_entry']
    r = T.d_loop((327, cs, As, 2e-9, 1.3 * R['events'][-1][0], X.net_entries(R)))
    assert r['ok']



def test_inner_profile_next_orders_and_its_rho_scaling():
    """(i): below the wall 𝒲(ξ) = K'³ρ⁻³ξ⁻³(1 − 8(κξ)⁻² + 64(κξ)⁻⁴ + O(ξ⁻⁶)) — z = α/a = t(1 − 2t²/3 + 13t⁴/15), t = 2/(κξ), from
    coth's series —; and 𝒲 depends on ρ and ξ only through ρξ (η_c(α) = −Λ₀ρξ)."""
    K3 = (8 / np.pi ** 3) * L.A / (L.m * np.sqrt(2 * L.c))
    kap = np.pi * np.sqrt(2 * L.c)
    xi = np.array([80.0, 160.0])
    r = RT.inner_coef(L, xi) * xi ** 3 / K3 - (1 - 8 / (kap * xi) ** 2)
    assert np.abs(r * (kap * xi) ** 4 / 64 - 1).max() < 0.005
    assert abs(RT.inner_coef(L, 15.0) * 15.0 ** 3 / K3 - 1 + 0.067) < 0.002          # 6.7% short of the power law at ξ = 15
    xs = np.array([-3.0, 0.5, 4.0, 30.0])
    for rho in (0.5, 2.0):
        assert np.allclose(RT.inner_coef(L, xs, rho), RT.inner_coef(L, rho * xs, 1.0), rtol=1e-12, atol=0)


def test_inner_threshold_lies_below_the_gibbs_value():
    """(ii): at the exact tie (A = 0.8) the inner problem's runner-up ends above θ = 0.01|A| at x = 6.2 and below it at x = 6.4 — its
    threshold x_θ ≈ 6.30, below the Gibbs value ln(|A|/θ − 1)/κ = 6.54: kernels apart creep toward the bump before the entry, the runner-up
    less, and it falls further behind."""
    K, c0 = grid()
    th = 0.01 * 0.8
    assert RT.near_tie_inner(L, K, 0.8, 0.05, 0.1, 0.0, 6.2, max_doublings=10)['w_r'] > th
    assert RT.near_tie_inner(L, K, 0.8, 0.05, 0.1, 0.0, 6.4, max_doublings=10)['w_r'] < th
    assert np.log(1 / 0.01 - 1) / (np.pi * np.sqrt(2 * L.c)) > 6.5


def test_the_partners_approach_creep_raises_the_threshold():
    """(ii)(a): one bump (A = 0.8) 0.03 off the kernels' midpoint, μ* = 0.0989: at x = μ*/ε^{1/2} = 7.4 the free loop forms the runner-up
    (its threshold 8.19: the runner-up, farther from the bump, gains on the entrant by its approach creep), while with both centres held
    fixed until the entrant's arrival it does not (6.52, the Gibbs value)."""
    import exp_transit_mass as T
    mu = T.near_tie_case(0.8, 0.03)[-1]
    assert abs(mu - 0.0989) < 5e-4
    eps = (mu / 7.4) ** 2
    assert T.h_formed(0.03, 'free', eps) and not T.h_formed(0.03, 'both', eps)


def test_unequal_masses_break_the_tie():
    """(ii)(b): an exact tie made by unequal masses (the bump at 0.08 between kernels 11 and 12, kernel 12's mass b = 1.626 so that both
    reach their walls together) does not split ½: the loop at ε = 10⁻⁸ forms kernel 12 alone, with 99.7% of the bump, where the path forms
    both."""
    import exp_transit_mass as T
    r = T.h_masses((0.08, 11, 12, 1e-8))
    assert abs(r['b'] - 1.626) < 2e-3 and r['share'] > 0.99 and r['formed'] == [12] and sorted(r['path']) == [11, 12]


def test_a_near_side_runner_up_is_passed_over():
    """(ii)(d): E's target 117 at ε = 10⁻⁴ (x = 0.52): kernel 14, between the third entrant (13) and its bump, never forms — its largest
    coefficient below θ/10 — as the entrant slides through it; the loop forms the path's kernels."""
    import exp_transit_mass as T
    r = T.h_near_side(1e-4)
    assert r['w14'] < 0.1 * r['theta'] and r['gap'] < 1e-3 and sorted(r['order']) == sorted(r['path'])


def test_the_estimate_overstates_the_loops_threshold_at_a_large_margin():
    """(iv): D's target 7, bound by an entry at the margin of a near tie (ε_c = 1.7·10⁻⁴): the loop does not form exactly the path's kernels
    at 0.5ε_c, where the estimate says it would — its own threshold is 0.27ε_c (G) — and does at 0.15ε_c."""
    import exp_transit_mass as T
    import exp_mirror_transit as X
    Lx, K, c0 = T.setup()
    cs, As = T.random_targets(T.SEED_D, T.N_D)[7]
    f = X.bumps_target(K, cs, As)
    R = RT.sliding_iss(Lx, K, f, c0, np.ones(21), polish=True)
    p = RT.recovery_threshold(Lx, K, f, c0, R, X.TH_ACT * np.abs(As).max())
    ec = min(p['eps_entry'], p['eps_creep'])
    assert p['eps_entry'] < p['eps_creep'] and 1.5e-4 < ec < 2e-4
    te, path = 1.3 * R['events'][-1][0], X.net_entries(R)
    assert not T.d_loop((7, cs, As, 0.5 * ec, te, path))['ok'] and T.d_loop((7, cs, As, 0.15 * ec, te, path))['ok']
