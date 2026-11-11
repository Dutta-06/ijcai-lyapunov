"""Proposition D‴.31 (PARADIGM.md §5d): the adaptive prior is a washout leak — its leak target lags the class barycentre by the barycentre's
exponentially forgotten increments, so the leak acts on the separations in full and on the barycentre through the high-pass iω/(iω + η) —
across the plan's switches too — and along every affine time-invariant prior the stationary response R to a persistent force and the bias
B toward a misplaced anchor satisfy RH + B = I (HR + B = I only where the zero-frequency leak commutes with H); the function-level response
to a force is at most the norm of the smallest target perturbation exerting it, under every leak S ≥ 0; and the desired regressor (the
adaptive law's regressor at the reference) leaves the lift, the linearisation and the first-order noise response unchanged while turning
the noise's own second-order force into the least-squares force of g = k_eE[e] — the parameters' fluctuation adding a force of the same
order that is not of that form."""
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.core import langevin as L
from src.core import adaptive_floquet as AF
from src.core import adaptive_prior as AP
from src.core import ot_leak as OT
from src.core import measurement_noise as MN

PER = 2 * np.pi
BUMPS = [(1.0, -1.5, 0.5), (-0.7, 0.0, 0.4), (0.8, 1.5, 0.5)]
G = np.array([1.0, 0.25, 0.25])
SIGMA0 = 0.15 * 0.047162567584665846
Q3 = np.ones(3) / 3
F_VEC = lambda X: sum(A * np.exp(-(np.asarray(X, float) - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS)
FP = lambda x: sum(-A * (x - C) / S ** 2 * np.exp(-(x - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS)
FPP = lambda x: sum(A * ((x - C) ** 2 / S ** 4 - 1 / S ** 2) * np.exp(-(x - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS)
XR = lambda t: np.array([2.0 * np.sin(t)])
XRD = lambda t: np.array([2.0 * np.cos(t)])
KE = 25.0


def _setup():
    net = L.SaturatedLiftedNetwork(1, a_max=4.0, c_max=3.0, s_min=0.3, s_max=2.0)
    atoms = net.unconstrained(np.array([[3 * A] for A, C, S in BUMPS]), np.array([[C] for A, C, S in BUMPS]), np.array([S for A, C, S in BUMPS]))
    cv = AF.CollapsedVariational(net, XR, PER, G, SIGMA0, np.inf, atoms, Q3)
    return net, atoms, cv


def test_the_prior_lags_the_barycentre_by_its_forgotten_increments():
    """(i) on the real loop (finite η, a perturbed start, classes (0, 0, 1, 2)): T̂ − z = ∫₀ᵗe^{−(η/k_e)(t−τ)}dT̂(τ) + e^{−(η/k_e)t}(T̂(0) − z(0)),
    the convolution integrated along the solution from the barycentre's velocity (the class mean of Θ̇)."""
    net, atoms, cv = _setup()
    eta = 5 * SIGMA0
    labels = np.array([0, 0, 1, 2]); m, p = 4, 3
    cl = AF.AdaptiveRealLoop(net, F_VEC, XR, XRD, PER, KE, G, SIGMA0, labels, 3, eta, weights=np.array([1, 1, 2, 2]) / 6)
    z0 = cl.lift(atoms)
    rng = np.random.default_rng(3)
    z0[1:1 + m * p] += 0.05 * rng.standard_normal(m * p)
    z0[1 + m * p:] += 0.05 * rng.standard_normal(3 * p)
    c = eta / KE
    T = 3 * PER

    def aug(t, y):
        z = y[:cl.state_dim(m)]
        dz = cl.rhs(t, z, m)
        dT = cl.class_means(dz[1:1 + m * p].reshape(m, p))
        return np.concatenate([dz, (dT - c * y[cl.state_dim(m):].reshape(3, p)).ravel()])
    y0 = np.concatenate([z0, np.zeros(3 * p)])
    sol = solve_ivp(aug, (0, T), y0, method='DOP853', rtol=1e-11, atol=1e-13)
    yT = sol.y[:, -1]
    Th0 = z0[1:1 + m * p].reshape(m, p); Z0 = z0[1 + m * p:].reshape(3, p)
    ThT = yT[1:1 + m * p].reshape(m, p); ZT = yT[1 + m * p:1 + m * p + 3 * p].reshape(3, p); W = yT[cl.state_dim(m):].reshape(3, p)
    lhs = cl.class_means(ThT) - ZT
    rhs = W + np.exp(-c * T) * (cl.class_means(Th0) - Z0)
    assert np.abs(lhs - rhs).max() < 1e-9 and np.abs(W).max() > 1e-4


def test_the_prior_lags_the_barycentre_across_a_plan_switch():
    """(i) on the averaged loop with the nearest-atom plan recomputed every step (η = 1, slow time): a particle placed past the bisector of
    two atoms switches its cell once the prior moves (T̂ jumps by about 0.5), and T̂ − z equals the forgotten increments — the smooth ones
    and the jump, each at its own time — plus e^{−ηs}(T̂(0) − z(0)), to the midpoint rule's accuracy."""
    net, atoms, cv = _setup()
    loop = OT.OTPriorLoop(net, F_VEC, XR, PER, G, SIGMA0, atoms, n_quad=128)
    eta, dt = 1.0, 0.005
    a0, a1, a2 = atoms
    Th0 = np.array([a0, a0, a1, a1, a2 + 0.4 * (a0 - a2), a2 + 0.55 * (a0 - a2)])
    apl = AP.AdaptivePriorLoop(loop, eta, plan='nearest')
    rec = []
    apl.run(Th0, atoms.copy(), 2.0, dt, plan_every=1, rec_every=1, callback=lambda s_, Th, Z, ad: rec.append((s_, Th.copy(), Z.copy())))
    bary = lambda Th, Z, gam: apl.bary(Th, gam, Z)
    gams = [AP.nearest_plan(Th, Z, G)[0] for _, Th, Z in rec]
    smooth = [bary(rec[k + 1][1], rec[k + 1][2], gams[k]) - bary(rec[k][1], rec[k][2], gams[k]) for k in range(len(rec) - 1)]
    jumps = [bary(rec[k + 1][1], rec[k + 1][2], gams[k + 1]) - bary(rec[k + 1][1], rec[k + 1][2], gams[k]) for k in range(len(rec) - 1)]
    s_ = np.array([r[0] for r in rec])
    assert max(np.abs(j).max() for j in jumps) > 0.1                                   # the plan switched
    T0 = bary(rec[0][1], rec[0][2], gams[0])
    for n in (200, 300, 400):                                                         # s = 1, 1.5, 2: after the switch at s ≈ 0.6
        lhs = bary(rec[n][1], rec[n][2], gams[n]) - rec[n][2]
        rhs = sum(np.exp(-eta * (s_[n] - 0.5 * (s_[k] + s_[k + 1]))) * smooth[k] + np.exp(-eta * (s_[n] - s_[k + 1])) * jumps[k] for k in range(n))
        rhs = rhs + np.exp(-eta * s_[n]) * (T0 - rec[0][2])
        jc = sum(np.exp(-eta * (s_[n] - s_[k + 1])) * jumps[k] for k in range(n))
        assert np.abs(lhs - rhs).max() < 2e-5 and np.abs(jc).max() > 0.05             # the identity (2.4·10⁻⁶ here, the midpoint rule's), the jump's share large


def test_the_pair_block_is_the_washed_out_leak_and_the_stationary_identity_holds_along_the_prior_family():
    """(i) linearised: the transfer function of the averaged collective block from a force on y to y is, mode by mode,
    1/(iω + h_j + σ₀(iω + κ)/(iω + η + κ)); (ii) at zero frequency RH + B = I with R = (H + S₀)⁻¹, S₀ = σ₀κ/(η + κ), for the adaptive prior
    (κ = 0: S₀ = 0, B = 0, R = H⁻¹), the anchored priors and the fixed prior (a scalar S₀: HR + B = I as well); for priors whose
    zero-frequency leak does not commute with H — fixed on the amplitudes (η = σ₀ elsewhere), one active face of the projected adaptive
    prior (η = ∞), a matrix anchor gain ηP_amp, a random anchor gain with K + Kᵀ ⪰ 0 — R and B from the dynamic block equal (H + S₀)⁻¹ and
    (H + S₀)⁻¹S₀ with S₀ by elimination, so RH + B = I and the per-mode h_jr_jj + b_jj = 1 hold, and HR + B = I fails (round 22, ID-5:
    derived from the block, not from the definitions)."""
    net, atoms, cv = _setup()
    H, Q = MN.noise_matrices(cv, FP)
    h, V = np.linalg.eigh(H)
    d = len(h); eta = SIGMA0
    for kap in (0.0, 0.1 * eta, eta, 10 * eta):
        A = MN.collective_block(H, SIGMA0, eta, kap)
        for om in (0.0, 0.3 * eta, 3 * eta, 30 * eta):
            Rw = np.linalg.solve(1j * om * np.eye(A.shape[0]) - A, np.vstack([np.eye(d), np.zeros((d, d))]))[:d]
            pred = V @ np.diag(MN.washout_transfer(h, SIGMA0, eta, om, kap)) @ V.T
            assert np.abs(Rw - pred).max() < 1e-9 * np.abs(pred).max()
        R, Bm = MN.stationary_operators(H, SIGMA0, eta, kap)
        S0 = MN.zero_frequency_leak(SIGMA0, eta, kap)
        assert np.abs(H @ R + Bm - np.eye(d)).max() < 1e-12
        assert np.abs(R - np.linalg.inv(H + S0 * np.eye(d))).max() < 1e-10 * np.abs(R).max()
    R, Bm = MN.stationary_operators(H, SIGMA0, 0.0)
    assert np.abs(H @ R + Bm - np.eye(d)).max() < 1e-12 and np.abs(Bm - SIGMA0 * np.linalg.inv(H + SIGMA0 * np.eye(d))).max() < 1e-12
    R, Bm = MN.stationary_operators(H, SIGMA0, eta, 0.0)
    assert np.abs(Bm).max() < 1e-12 and np.abs(R - np.linalg.inv(H)).max() < 1e-8 * np.abs(R).max()
    amp = np.isin(np.arange(d), [0, 3, 6])                                            # whitened (k, l) ↦ 3k + l, l = 0 the amplitude
    Mr = np.random.default_rng(3).standard_normal((d, d))
    Kr = eta * (Mr @ Mr.T / np.linalg.norm(Mr @ Mr.T, 2) + 0.5 * (Mr - Mr.T) / np.linalg.norm(Mr - Mr.T, 2))
    for e_, K, fixed in ((eta, None, amp), (np.inf, None, np.arange(d) == 0), (eta, eta * np.diag(amp.astype(float)), None), (eta, Kr, None)):
        R, Bm = MN.anchored_operators(H, SIGMA0, e_, K, fixed)
        S0m = MN.anchored_leak(H, SIGMA0, e_, K, fixed)
        R0, B0 = MN.zero_frequency_operators(H, S0m)
        assert np.abs(R - R0).max() < 1e-10 * np.abs(R0).max() and np.abs(Bm - B0).max() < 1e-10
        assert np.abs(R @ H + Bm - np.eye(d)).max() < 1e-12
        assert np.abs(H @ S0m - S0m @ H).max() > 1e-5 and np.abs(H @ R + Bm - np.eye(d)).max() > 0.1
        assert np.abs(h * np.einsum('ij,ik,kj->j', V, R, V) + np.einsum('ij,ik,kj->j', V, Bm, V) - 1).max() < 1e-10
    A, Bin = MN.anchored_block(H, SIGMA0, eta, eta * np.eye(d))                      # a scalar anchor gain: the prior family's block
    assert np.abs(A - MN.collective_block(H, SIGMA0, eta, eta)).max() == 0.0


def test_the_response_to_a_force_is_at_most_its_smallest_equivalent_target_perturbation_under_every_leak():
    """(iii): Fᵀ(H + S)⁻¹H(H + S)⁻¹F ≤ FᵀH⁻¹F for S ≥ 0, with equality at S = 0; FᵀH⁻¹F is the squared norm of g = ĴH⁻¹F, the smallest g with
    ⟨Ĵᵀg⟩ = F; a least-squares force ⟨Ĵᵀg⟩ has a response at most ‖g‖²_ρ, and at S = 0 the response Ĵy is the projection of g."""
    net, atoms, cv = _setup()
    N = 512
    t = np.arange(N) * PER / N
    J = np.array([cv.J(s)[0] for s in t])
    H = J.T @ J / N
    rng = np.random.default_rng(0)
    d = H.shape[0]
    for trial in range(20):
        F = rng.standard_normal(d)
        M = rng.standard_normal((d, d))
        S = [0.0, SIGMA0, 10 * SIGMA0, 0.01 * M @ M.T][trial % 4]
        assert MN.response_size(H, F, S) <= MN.equivalent_perturbation(H, F) * (1 + 1e-10)
        g_min = J @ np.linalg.solve(H, F)
        assert abs(np.mean(g_min ** 2) / MN.equivalent_perturbation(H, F) - 1) < 1e-9 and np.abs(J.T @ g_min / N - F).max() < 1e-9 * np.abs(F).max()
        g = rng.standard_normal(N)
        Fg = J.T @ g / N
        assert MN.response_size(H, Fg, S) <= np.mean(g ** 2)
        y = np.linalg.solve(H, Fg)
        Pg = J @ np.linalg.solve(J.T @ J, J.T @ g)
        assert np.abs(J @ y - Pg).max() < 1e-8 * np.abs(Pg).max()
    assert abs(MN.response_size(H, F, 0.0) / MN.equivalent_perturbation(H, F) - 1) < 1e-10


def test_the_desired_regressor_keeps_the_lift_and_the_first_order_response():
    """(iv)(a)–(b): without noise the lift stays exact under DR and DR*; on a common noise path DR's run differs from MR's by the second-order
    terms only — the difference falls like the noise's variance while each deviation falls like its standard deviation."""
    net, atoms, cv = _setup()
    grid = 2 * np.sin(np.arange(64) * PER / 64)
    runs = {}
    for reg, comp in (('measured', 'measured'), ('reference', 'measured'), ('reference', 'reference')):
        cl = MN.NoisyRealLoop(net, F_VEC, XR, XRD, PER, KE, G, SIGMA0, np.arange(3), 3, np.inf, weights=Q3, regressor=reg, compensation=comp)
        z0 = cl.lift(atoms)
        rec = lambda t, z: np.array([cl.f_hat(np.array([x]), z[1:10].reshape(3, 3))[0] for x in grid]) - F_VEC(grid)
        clean = MN.simulate_noisy(cl, z0, 3, 2.0, 1e-3, 100.0, 0.0, np.random.default_rng(0), 2000, rec)
        assert np.abs(clean[-1]).max() < 1e-12                                       # the lift is a solution
        for std in (1e-3, 1e-2):
            runs[(reg, comp, std)] = MN.simulate_noisy(cl, z0, 3, 5.0, 1e-3, 100.0, std ** 2, np.random.default_rng(7), 5000, rec)[-1]
    dev = {s: np.sqrt(np.mean(runs[('measured', 'measured', s)] ** 2)) for s in (1e-3, 1e-2)}
    dif = {s: np.sqrt(np.mean((runs[('measured', 'measured', s)] - runs[('reference', 'measured', s)]) ** 2)) for s in (1e-3, 1e-2)}
    assert 7 < dev[1e-2] / dev[1e-3] < 14                                              # first order: ∝ std
    assert 70 < dif[1e-2] / dif[1e-3] < 140 and dif[1e-3] < 0.05 * dev[1e-3]            # MR − DR: ∝ std², small at small noise
    dif_star = np.sqrt(np.mean((runs[('measured', 'measured', 1e-3)] - runs[('reference', 'reference', 1e-3)]) ** 2))
    assert dif_star < 0.1 * dev[1e-3]                                                 # DR*: the same response up to O(a/k_e)
    # the regressor's point (ν = 0, x ≠ x_m): DR's law is G∇_ϑ(eΦ̃(x_m; ϑ)), MR's G∇_ϑ(eΦ̃(x; ϑ)); with ν and e = 0, DR's is G∇_ϑ(νΦ̃(x_m; ϑ))
    for reg in ('measured', 'reference'):
        cl = MN.NoisyRealLoop(net, F_VEC, XR, XRD, PER, KE, G, SIGMA0, np.arange(3), 3, np.inf, weights=Q3, regressor=reg)
        z = cl.lift(atoms); t = 0.7; xm = XR(t); z[0] = xm[0] + 0.3
        pt = xm if reg == 'reference' else z[:1]
        law = cl.rhs_noisy(t, z, 3, np.zeros(1))[1:10].reshape(3, 3)
        assert np.abs(law - G[None, :] * net.grad_g(pt, np.array([0.3]), atoms)).max() < 1e-12
    z = cl.lift(atoms); z[0] = xm[0]
    law = cl.rhs_noisy(t, z, 3, np.array([0.02]))[1:10].reshape(3, 3)
    assert np.abs(law - G[None, :] * net.grad_g(xm, np.array([0.02]), atoms)).max() < 1e-12


def test_the_desired_regressor_turns_the_errors_in_variables_force_into_a_least_squares_force():
    """(iv)(c): DR's and DR*'s forces are ⟨Ĵᵀg⟩/k_e (real time) with g = k_eE[e], |g| ≤ ½|f″|Eν²(1 + 2|a|/k_e) along the orbit (DR) and
    ≤ ½|f″|Eν²(k_e/(k_e − a))³ (DR*), their adaptive-prior responses at most ‖g‖²_ρ; MR's force exceeds its least-squares part by
    F₁ = k_e⟨Ĵ′ᵀE[ε_m²]⟩, whose projected mean tracking error F₁ᵀH⁻¹F₁ is three orders of magnitude above DR's bias here."""
    net, atoms, cv = _setup()
    H, Q = MN.noise_matrices(cv, FP)
    N = 1024
    t = np.arange(N) * PER / N
    J = np.array([cv.J(s)[0] for s in t])
    x = 2 * np.sin(t); a = FP(x); b2 = FPP(x)
    for ke, B in ((25.0, 100.0), (100.0, 1000.0)):
        Fr, g = MN.eiv_force(cv, FP, FPP, ke, B, 1.0, variant='reference', return_g=True)
        assert np.abs(Fr - J.T @ (g / ke) / N).max() < 1e-12
        assert np.all(np.abs(g) <= 0.5 * np.abs(b2) * (1 + 2 * np.abs(a) / ke) + 1e-15)
        Fb, gb = MN.eiv_force(cv, FP, FPP, ke, B, 1.0, variant='reference_both', return_g=True)
        assert np.all(np.abs(gb) <= 0.5 * np.abs(b2) * (ke / (ke - a)) ** 3 + 1e-15)
        for FF, gg in ((Fr, g), (Fb, gb)):
            assert MN.equivalent_perturbation(H, ke * FF) <= np.mean(gg ** 2) * (1 + 1e-9)
        Fm = MN.eiv_force(cv, FP, FPP, ke, B, 1.0, variant='measured')
        F1 = ke * (Fm - Fr)
        assert MN.equivalent_perturbation(H, F1) > 1e3 * MN.equivalent_perturbation(H, ke * Fr)


def test_the_noise_part_of_the_force_is_the_hessian_contraction_and_the_parameters_fluctuation_is_not_least_squares():
    """The quasi-static forces of `eiv_force` against the full second-order mean of `second_order_mean` (the periodic covariance contracted
    with the true field's second derivatives, k_e = 25, B = 100, η = ∞): the noise parts agree for MR, DR and DR* (the closed forms of E[ε_m²]
    and E[e] right, their signs and the regressor's part); the parameters' fluctuation changes MR's bias by less than 10⁻⁵ but DR's and DR*'s
    by more than 1% (round 21, D3r-1) — and DR's full bias stays three orders of magnitude below MR's."""
    net, atoms, cv = _setup()
    H, Q = MN.noise_matrices(cv, FP)
    full = {}
    for reg, comp in (('measured', 'measured'), ('reference', 'measured'), ('reference', 'reference')):
        v = 'measured' if reg == 'measured' else ('reference' if comp == 'measured' else 'reference_both')
        cl = MN.NoisyRealLoop(net, F_VEC, XR, XRD, PER, KE, G, SIGMA0, np.arange(3), 3, np.inf, weights=Q3, regressor=reg, compensation=comp)
        out = MN.second_order_mean(cl, cv, FP, 100.0, 1.0, N=128)
        Fq = KE * MN.eiv_force(cv, FP, FPP, KE, 100.0, 1.0, variant=v)
        assert np.linalg.norm(out['F_noise'] - Fq) < 3e-3 * np.linalg.norm(Fq)
        b = MN.mean_biases(H, out)
        full[v] = b['full']
        if v == 'measured':
            assert abs(b['full'] / b['noise'] - 1) < 1e-5
        else:
            assert abs(b['full'] / b['noise'] - 1) > 0.01
    assert full['reference'] < 1e-3 * full['measured'] and full['reference_both'] < 1e-3 * full['measured']


def test_with_the_regressor_at_the_reference_the_separations_drift_is_the_mean_errors_term_alone():
    """Corollary D‴.33: under DR the separations' second-order drift is ⟨E[e]G^{1/2}∇²_ϑΦ̃G^{1/2}⟩ alone — MR's drift minus its regressor term,
    since DR's mean error is MR's at the lift — negative definite in every class on the instance (no noise threshold), while MR's is 0.30
    in class 3; under DR* it is positive but O(1/k_e): k_e times its largest eigenvalue is the same at k_e = 25 and 50 (B = 4k_e), so its
    threshold does not fall with k_e, while MR's largest eigenvalue does not decrease; and the averaged noise-only residual (I − P)g/k_e at
    the adaptive prior's stationary point is below 1% of g/k_e — so it is the parameters' part, not this residual, that sets the stationary
    mean error (the next tests; referee round 23, D3v-1)."""
    net, atoms, cv = _setup()
    fpp = lambda x: sum(A * ((x - C) ** 2 / S ** 4 - 1 / S ** 2) * np.exp(-(x - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS)
    S_mr = MN.separation_drift(cv, FP, fpp, 25.0, 100.0, 1.0, N=128)
    S_dr = MN.separation_drift(cv, FP, fpp, 25.0, 100.0, 1.0, N=128, variant='reference')
    # MR's E[e] term alone, computed through the series with E[ε_m²] = 0 and MR's quasi-static mean error
    t = np.arange(128) * PER / 128
    Ee = np.array([MN.mean_error(FP(2 * np.sin(s)), fpp(2 * np.sin(s)), 25.0, 100.0, 1.0) for s in t])
    S_e = MN.separation_drift_series(cv, t, Ee, np.zeros(128)).mean(axis=0)
    assert np.abs(S_dr - S_e).max() < 1e-12 * np.abs(S_mr).max()
    assert all(np.linalg.eigvalsh(s)[-1] < 0 for s in S_dr) and np.linalg.eigvalsh(S_mr[2])[-1] > 0.25
    lam = {ke: (max(np.linalg.eigvalsh(s)[-1] for s in MN.separation_drift(cv, FP, fpp, ke, 4 * ke, 1.0, N=128, variant='reference_both')),
                max(np.linalg.eigvalsh(s)[-1] for s in MN.separation_drift(cv, FP, fpp, ke, 4 * ke, 1.0, N=128))) for ke in (25.0, 50.0)}
    assert lam[25.0][0] > 0 and abs(50.0 * lam[50.0][0] / (25.0 * lam[25.0][0]) - 1) < 0.01
    assert lam[50.0][1] > lam[25.0][1] > 0.25
    H = MN.noise_matrices(cv, FP)[0]
    F, g = MN.eiv_force(cv, FP, fpp, 25.0, 100.0, 1.0, N=256, variant='reference', return_g=True)
    t = np.arange(256) * PER / 256
    J = np.array([AF.whitened_jacobian(net, np.array([2 * np.sin(s)]), cv.atoms, cv.w, cv.G)[0] for s in t])
    res = g - J @ np.linalg.solve(H, 25.0 * F)
    assert np.sqrt(np.mean(res ** 2)) < 0.01 * np.sqrt(np.mean(g ** 2))


def test_the_desired_compensation_laws_drift_is_the_quadrature_of_its_mean_error():
    """Corollary D‴.33 (round 23, D3v-7): DR*'s separation drift at the lift against an independent quadrature — the Hessian of the
    single-particle feature by second differences of f̂ itself, DR*'s mean error ½∂²_xf·Var(e)/(k_e − a) with Var(e) =
    k_e²Eν²/((k_e − a)(k_e − a + B)) typed here — to 10⁻⁴, so a rescaled DR* drift (the doubled one passed every test of round 22) is caught."""
    net, atoms, cv = _setup()
    ke, B, N = 25.0, 100.0, 64
    S = MN.separation_drift(cv, FP, FPP, ke, B, 1.0, N=N, variant='reference_both')
    phi = lambda x, th: net.f_hat(np.array([x]), th[None, :])[0]
    E = 1e-3 * np.eye(3)
    hess = lambda x, th: np.array([[(phi(x, th + E[i] + E[j]) - phi(x, th + E[i] - E[j]) - phi(x, th - E[i] + E[j]) + phi(x, th - E[i] - E[j]))
                                    / 4e-6 for j in range(3)] for i in range(3)])
    g = np.sqrt(G)
    Sq = np.zeros_like(S)
    for s_ in np.arange(N) * PER / N:
        x = 2 * np.sin(s_); c = ke - FP(x)
        Ee = 0.5 * FPP(x) * (ke ** 2 / (c * (c + B))) / c
        for k, th in enumerate(atoms):
            Sq[k] += Ee * g[:, None] * hess(x, th) * g[None, :] / N
    assert np.abs(S - Sq).max() < 1e-4 * np.abs(S).max()
    assert all(np.linalg.eigvalsh(sk)[0] > 0 for sk in S)


def test_at_the_stationary_point_the_mean_error_carries_the_parameters_part():
    """Corollary D‴.33(iii) (round 23, D3v-1): the exact periodic second-order mean of DR's K₀-atom loop under the adaptive prior is the sum
    of its noise part and its parameters' part (the forcing is linear in the covariance); the parameters' part dominates the stationary mean
    error (its rms more than twice the noise part's) and lies largely in the tangent space, which the balance does not remove; and the
    drift there, the cross term included, is indefinite."""
    net, atoms, cv = _setup()
    cl = MN.NoisyRealLoop(net, F_VEC, XR, XRD, PER, KE, G, SIGMA0, np.arange(3), 3, np.inf, weights=Q3, regressor='reference')
    o = MN.second_order_mean(cl, cv, FP, 100.0, 1.0, N=128, periodic=True)
    e, en, ep = o['m'][:, 0], o['m_noise'][:, 0], o['m_param'][:, 0]
    rms = lambda v: np.sqrt(np.mean(v ** 2))
    assert np.abs(e - en - ep).max() < 1e-8 * np.abs(e).max()
    assert rms(ep) > 2 * rms(en)
    H = MN.noise_matrices(cv, FP)[0]
    Fe = o['J'].T @ ep / len(ep)
    assert np.sqrt(Fe @ np.linalg.solve(H, Fe)) > 0.5 * rms(ep)
    S = MN.separation_drift_series(cv, o['t'], e, np.zeros_like(e), Eeth=o['Eeth']).mean(axis=0)
    lam = np.linalg.eigvalsh(S[0])
    assert lam[0] < 0 < lam[-1]


def test_the_period_maps_continue_the_direct_transient_and_drs_drift_turns_indefinite():
    """`long_transient` (round 23, D3v-2): over the first periods its period-averaged drift (the mean error's two parts and the cross term)
    is `second_order_transient`'s on the same grid to 10⁻⁷; and along DR's transient the drift, negative definite near the lift, turns
    indefinite — in class 3 after 300 time units and within 800 periods."""
    net, atoms, cv = _setup()
    cl = MN.NoisyRealLoop(net, F_VEC, XR, XRD, PER, KE, G, SIGMA0, np.arange(3), 3, np.inf, weights=Q3, regressor='reference')
    lt = MN.long_transient(cl, cv, FP, 100.0, 1.0, 800, N=126, keep=np.arange(4))
    o = MN.second_order_transient(cl, cv, FP, 100.0, 1.0, 504 * (PER / 126), dt=PER / 126)
    Sd = MN.separation_drift_series(cv, o['t'], o['Ee'], np.zeros_like(o['Ee']), Eeth=o['Eeth'])
    for n in range(4):
        r = lt['rec'][n]
        assert np.abs(Sd[n * 126:(n + 1) * 126].mean(axis=0) - (r['S_noise'] + r['S_param'] + r['C'])).max() < 1e-7 * np.abs(Sd).max()
    lam, t = lt['lam_tot'][:, 2], lt['t']
    assert np.all(lam[t < 300] < 0) and np.any(lam > 0)
