"""Theorem D‴ (PARADIGM.md §5d): the averaged slow loop, the three terms of the coupling identity, the residual-weighted
curvature against Proposition K, the deterministic neighbourhood search, the constants of the theorem (the bootstrap's
constraints), and a coupled run from a near-stationary warm start where hypothesis (H) is verified."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.core import langevin as L
from src.core import periodic as PO
from src.core import chaos_learning as CL

PER = 2 * np.pi
BUMPS = [(1.0, -1.5, 0.5), (-0.7, 0.0, 0.4), (0.8, 1.5, 0.5)]
G = np.array([1.0, 0.25, 0.25])
F_VEC = lambda X: sum(A * np.exp(-(X - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS)
X_REF = lambda t: np.array([2.0 * np.sin(t)])
X_REF_DOT = lambda t: np.array([2.0 * np.cos(t)])
SIGMA_CRIT = 0.047162567584665846          # Proposition K for this instance (asserted in test 3)
LC = dict(L_th_G=4.0435377, L_2_G=11.1111111)   # searched whitened constants (asserted in test 4)


def _instance(sigma0=0.014, n_quad=400):
    net = L.SaturatedLiftedNetwork(1, a_max=4.0, c_max=3.0, s_min=0.3, s_max=2.0)
    loop = CL.AveragedSlowLoop(net, F_VEC, X_REF, PER, G, sigma0, n_quad=n_quad)
    cl = PO.ClosedLoop(net, lambda Xa: F_VEC(Xa[:, 0])[:, None], X_REF, X_REF_DOT, PER, 25.0, G, sigma0)
    return net, loop, cl


def test_averaged_slow_loop_is_the_G_gradient_flow_of_J_m():
    """(𝔄_m): the drift equals −G∇J_m with J_m of periodic.ClosedLoop (same quadrature grid), to machine precision."""
    net, loop, cl = _instance()
    orb = cl.orbit(N=400)
    rng = np.random.default_rng(3)
    Th = CL.cold_start(net, 12, rng)
    assert np.abs(loop.grad_J(Th) - cl.grad_J_m(Th, orb)).max() < 1e-13
    assert np.abs(loop.drift(Th) + G[None, :] * cl.grad_J_m(Th, orb)).max() < 1e-13
    J = 12 * loop.L_bar(Th) + 0.5 * loop.sigma0 * float(np.sum(Th ** 2 / G[None, :]))
    assert abs(J - cl.J_m(Th, orb)) < 1e-12


def test_coupling_identity_three_terms_and_the_sign_of_the_interaction():
    """½dD/ds = T_A + T_B + T_C exactly; the interaction term T_B equals −‖f̂_Θ − f̂_Θ̄‖² up to a relative O(D^{1/2}) correction
    (flat convexity: non-positive at linear order) and obeys (3.1) T_B ≤ (L₂^G)²D²/16 at every scale of δ."""
    net, loop, _ = _instance()
    rng = np.random.default_rng(3)
    Th = CL.cold_start(net, 12, rng)
    Thb = Th + 0.05 * rng.standard_normal(Th.shape)
    fh_mf = loop.f_hat(CL.cold_start(net, 500, rng))                     # a stand-in for f̂[μ_s]
    d = loop.decomposition(Th, Thb, fh_mf)
    assert abs(d['total'] - d['half_dD']) < 1e-15 * max(1.0, abs(d['half_dD']))
    ratios = []
    for eps in (1.0, 0.3, 0.1, 0.03, 0.01):
        dd = loop.decomposition(Thb + eps * (Th - Thb), Thb)
        assert dd['T_C'] == 0.0
        assert dd['T_B'] <= LC['L_2_G'] ** 2 * dd['D'] ** 2 / 16                     # (3.1)
        ratios.append(dd['T_B'] / (-dd['h_norm'] ** 2))
    assert abs(ratios[-1] - 1) < 1e-3 and abs(ratios[0] - 1) < 0.05          # T_B/(−‖h‖²) → 1 as δ → 0 (hence T_B < 0 there)
    assert all(abs(r - 1) <= abs(ratios[0] - 1) + 1e-12 for r in ratios)


def test_residual_curvature_at_a_collapsed_configuration_is_proposition_K_and_the_sphere_search_is_first_order_exact():
    """κ_R at the atoms of a collapsed configuration with the collapsed residual equals λ_max(G^{1/2}R₂G^{1/2}) of Proposition K,
    which at σ₀ = σ_crit equals σ₀; the collapsed critical point has zero drift.  The deterministic neighbourhood search over
    the sphere of radius ρ agrees with the first-order prediction λ_max(atom) + ρ‖G^{1/2}∇_ϑλ_max‖ to O(ρ²) (relative error
    decreasing with ρ) — it finds the supremum, not a random point."""
    net, loop, cl = _instance()
    orb = cl.orbit(N=400)
    g = net.unconstrained(np.array([[1.0]]), np.array([[0.0]]), np.array([0.8]))[0]
    s_crit, tb, R2 = PO.sigma_crit(cl, orb, g)
    assert abs(s_crit - SIGMA_CRIT) < 1e-9
    loop_c = CL.AveragedSlowLoop(net, F_VEC, X_REF, PER, G, s_crit, n_quad=400)
    Thc = np.tile(tb[None, :], (5, 1))
    r = loop_c.residual(Thc)
    kr = loop_c.kappa_R(r, Thc)
    lam = np.linalg.eigvalsh(np.sqrt(G)[:, None] * R2 * np.sqrt(G)[None, :])[-1]
    assert np.abs(kr - lam).max() < 1e-8 and abs(lam - s_crit) < 1e-9
    assert np.abs(loop_c.drift(Thc)).max() < 1e-12
    # the sphere search against the first-order expansion of ϑ ↦ λ_max(R^G(ϑ)) at the atom (frozen residual)
    h = 1e-6
    grad = np.array([(loop_c.kappa_R(r, tb[None, :] + h * np.sqrt(G[j]) * np.eye(3)[j][None, :])[0]
                      - loop_c.kappa_R(r, tb[None, :] - h * np.sqrt(G[j]) * np.eye(3)[j][None, :])[0]) / (2 * h) for j in range(3)])
    errs = []
    for rho in (0.04, 0.02, 0.01):
        found = loop_c.kappa_R_neighbourhood(r, tb[None, :], rho, n_dirs=400, radii=(1.0,))
        pred = lam + rho * np.linalg.norm(grad)
        errs.append(abs(found - pred) / (rho * np.linalg.norm(grad)))
        assert found >= lam
    assert errs[0] < 0.05 and errs[-1] < errs[0]                              # O(ρ) relative error, decreasing with ρ


def test_whitened_constants_and_the_bootstrap_constraints_of_theorem_constants():
    """The searched L_ϑ^G, L₂^G reproduce the values of the document and lie below the closed-form bounds; the constants
    returned by theorem_constants satisfy the four constraints the bootstrap of the proof needs: (L₂^G)²ϱ²/16 ≤ κ/8,
    L₂^Gθ ≤ κ/8, 4e^{2Λ₁}(L_ϑ^G)²θ²/κ² ≤ ϱ²/4 and 2e^{Λ₁}κ⁻¹[(L_ϑ^G)²ϱ/2 + L_ϑ^Gθ] ≤ 3ϱ′/4; Λ₁ of a profile integrates exactly
    to s₁ (referee round 8, D3-6)."""
    net, loop, _ = _instance()
    Lc = CL.whitened_lipschitz_constants(net, G, 2.5, n_samples=20000)
    assert abs(Lc['L_th_G'] - LC['L_th_G']) < 1e-3 and abs(Lc['L_2_G'] - LC['L_2_G']) < 1e-3
    assert Lc['L_th_G'] <= Lc['L_th_G_closed'] * (1 + 1e-9) and Lc['L_2_G'] <= Lc['L_2_G_closed'] * (1 + 1e-9)
    for kappa, Lam, rp in ((0.005, 0.0, 0.3), (0.002, 1.5, 0.1), (0.02, 0.0, 1.0)):
        c = CL.theorem_constants(loop, Lc['L_th_G'], Lc['L_2_G'], R0=3.0, kappa=kappa, Lambda1=Lam, rho_prime=rp)
        Lt, L2, rho, th = Lc['L_th_G'], Lc['L_2_G'], c['rho'], c['theta']
        assert L2 ** 2 * rho ** 2 / 16 <= kappa / 8 * (1 + 1e-12)
        assert L2 * th <= kappa / 8 * (1 + 1e-12)
        assert 4 * np.exp(2 * Lam) * Lt ** 2 * th ** 2 / kappa ** 2 <= rho ** 2 / 4 * (1 + 1e-12)
        assert 2 * np.exp(Lam) / kappa * (Lt ** 2 * rho / 2 + Lt * th) <= 0.75 * rp * (1 + 1e-12)
        assert c['Delta'] > 0 and c['m0'] == 16 * net.K_Phi ** 2 / th ** 2 and c['R_Theta'] >= 3.0
    s = np.linspace(0, 10, 11)
    assert CL.Lambda1_from_profile(s, 0.01 + 0 * s, 0.005, 5.0) == 0.0
    assert abs(CL.Lambda1_from_profile(s, -0.005 + 0 * s, 0.005, 5.5) - 0.055) < 1e-12


def test_hypothesis_H_at_the_three_atom_stationary_measure_and_the_weight_sensitivity_law():
    """At σ₀ = 0.15σ_crit the three-atom critical point of J₃ is a stationary measure of the mean-field flow (zero drift) at which
    (H) holds with ϱ′ = 0.1: κ_R over the sphere of radius ϱ′ around the atoms is below σ₀ (margin κ > 0), while σ₀ < σ_crit (the
    collapsed configuration is unstable: the window of Remark D‴.1).  Along a coupled run (m = 16, 64 against M = 512) from the
    jittered mixture μ₀, κ_R[μ*] over the neighbourhood of the reference's support stays below σ₀ (hypothesis (H) with Λ₁ = 0),
    the two curvature statistics differ by at most the proof's bound (0.2) L₂^G‖r_ref − r*‖_{L¹(ρ)}, D₆₄ < D₁₆, the function-level
    distance obeys ‖f̂_Θ − f̂_Θ̄‖ ≤ L_ϑ^G·mean‖δ_i‖, and the (vacuous) bound C₁/m holds.  Proposition D‴.5: the weighted K-atom
    stationary point ϑ*(w) replicated with multiplicities n_k = m·w_k is a critical point of J_m; the linearisation ϑ*(w̄ + dw) −
    ϑ*(w̄) = S·dw has a remainder of order |dw|² (checked at two scales); the splitting margin separates admissible weights
    (a strict local minimiser) from a saddle; and on the run the per-cluster mean deviations of the
    m-system from its reference point along the predicted vectors S(w^m − w^M) (cosine ≥ 0.8) with D within a factor 2 of
    Σ_k w_k‖S_k(w^m − w^M)‖²_{G⁻¹}, although at s = 60 ≪ 1/κ the systems are still relaxing toward their stationary points."""
    sigma0 = 0.15 * SIGMA_CRIT
    net, loop, _ = _instance(sigma0=sigma0, n_quad=256)
    Th0 = net.unconstrained(np.array([[3 * A] for A, C, S in BUMPS]), np.array([[C] for A, C, S in BUMPS]), np.array([S for A, C, S in BUMPS]))
    atoms, gnorm, eig = CL.stationary_atoms(loop, Th0)
    assert gnorm < 1e-12 and np.abs(loop.drift(atoms)).max() < 1e-12 and eig.min() > 0
    r_star = loop.residual(atoms)
    rho_p = 0.1
    k_nb = loop.kappa_R_neighbourhood(r_star, atoms, rho_p, n_dirs=200, radii=(1.0, 0.5))
    kappa = sigma0 - k_nb
    assert kappa > 0                                                         # (H) at μ* (σ₀ = 0.15σ_crit is inside the learning regime by construction)
    # Proposition D‴.5, structural part: the weighted stationary point and its replication
    wbar = np.ones(3) / 3
    A3, res3 = CL.weighted_stationary_atoms(loop, wbar, atoms)
    assert res3 < 1e-12 and np.abs(A3 - atoms).max() < 1e-9
    w16 = np.array([6, 7, 3]) / 16.0
    A16, res16 = CL.weighted_stationary_atoms(loop, w16, atoms)
    Th16 = np.repeat(A16, [6, 7, 3], axis=0)                                 # 16 particles with multiplicities n_k
    assert res16 < 1e-12 and np.abs(loop.grad_J(Th16)).max() < 1e-10         # a critical point of J₁₆
    S = CL.weight_sensitivity(loop, wbar, atoms)
    mg_ok, _, _ = CL.splitting_margin(loop, np.array([7, 5, 4]) / 16.0, atoms)
    mg_bad, _, _ = CL.splitting_margin(loop, np.array([9, 4, 3]) / 16.0, atoms)
    assert mg_ok > 0 > mg_bad                                                # (7,5,4)/16 admissible, (9,4,3)/16 a saddle (referee round 9, D3b-1)
    for eps in (0.05, 0.025):
        dw = eps * np.array([1.0, -0.5, -0.5])
        Aw, _ = CL.weighted_stationary_atoms(loop, wbar + dw, atoms)
        lin = np.einsum('kij,j->ki', S, dw)
        err = np.abs(Aw - atoms - lin).max()
        assert err < 0.2 * eps ** 2 * 100                                      # second-order remainder (|dw|² scale)
        if eps == 0.05:
            err05 = err
    assert err < err05 / 3                                                    # … decreasing like |dw|²
    # the coupled run
    rng = np.random.default_rng(11)
    Th_all = CL.mixture_start(atoms, 512, rng, jitter=0.02, G=G)
    rec = loop.coupled_run(Th_all, [16, 64], S_end=60.0, dt=0.02, rec_every=50, rho_nbhd=rho_p, n_dirs=16, rng=rng, nb_rows_max=128,
                           r_fixed=r_star)
    assert rec['kR_ref_nbhd_fixed'].max() < sigma0                           # (H) along the run with the stationary residual (lower estimate of the sup)
    assert np.all(np.abs(rec['kR_ref_nbhd'] - rec['kR_ref_nbhd_fixed']) <= LC['L_2_G'] * rec['resid_dev_L1'] * (1 + 1e-9))   # (0.2) (a loose bound: slack ≥ 10², referee round 9)
    kappa_run = min(kappa, sigma0 - rec['kR_ref_nbhd_fixed'].max())
    c = CL.theorem_constants(loop, LC['L_th_G'], LC['L_2_G'], R0=float(np.sqrt(np.sum(Th_all ** 2 / G[None, :], axis=1)).max()),
                             kappa=kappa_run, Lambda1=0.0, rho_prime=rho_p)
    labM, wM = CL.cluster_weights(rec['final'][512], atoms, G)
    for m in (16, 64):
        D = rec['D'][m]
        assert D.max() < c['C1'] / m                                         # the theorem's bound (vacuous here: C₁ is astronomically large)
        assert np.all(rec['hnorm'][m] <= LC['L_th_G'] * rec['dev1'][m] * (1 + 1e-9))
        lab, wm = CL.cluster_weights(rec['final'][m], atoms, G)
        dw = wm - wM
        pred = np.einsum('kij,j->ki', S, dw) / np.sqrt(G)[None, :]           # predicted per-cluster displacement, whitened
        dev = (rec['final'][m] - rec['final'][512][:m]) / np.sqrt(G)[None, :]
        for k in range(3):
            if np.any(lab == k):
                meas = dev[lab == k].mean(0)
                assert np.dot(meas, pred[k]) / (np.linalg.norm(meas) * np.linalg.norm(pred[k])) > 0.8
        D_pred = float(np.sum(wm * np.sum(pred ** 2, axis=1)))
        assert 0.5 < D[-1] / D_pred < 2.0
    assert rec['D'][64].max() < rec['D'][16].max()
