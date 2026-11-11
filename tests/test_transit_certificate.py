"""Theorem D‴.55: the transit's certificate (`experiments/exp_transit_certificate.py`)."""
import os
import sys
import warnings

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'experiments')))
warnings.filterwarnings('ignore')

import exp_transit_certificate as TC
import exp_transit_mass as TM


def test_one_feature_certificate_is_the_normalised_profile():
    L, K, c0 = TM.setup()
    cs, As = np.array([0.137]), np.array([-0.7])
    cert = TC.certificate(K, c0, cs, As)
    e = cert['e'][0]
    assert e == int(np.argmin(np.abs(c0 - 0.137)))
    # η(c_e) = 1, and η elsewhere is the sample-Gram profile normalised by the nearest kernel's, signed
    g = K.phi(c0).T @ K.phi(cs)[:, 0] / K.xs.size
    assert abs(cert['eta'][e] + 1.0) < 1e-12                           # sign(A) = −1
    assert np.allclose(cert['eta'], -g / g[e])
    assert 0 < cert['m'] < 1 and cert['q'] == 0.0


def test_certificate_depends_on_signs_and_positions_only():
    L, K, c0 = TM.setup()
    cs = np.array([-0.5, 0.05, 0.62])
    a = TC.certificate(K, c0, cs, np.array([0.9, -0.6, 0.75]))
    b = TC.certificate(K, c0, cs, np.array([1e-3, -50.0, 2.0]))
    assert np.allclose(a['eta'], b['eta']) and a['e'] == b['e'] and a['m'] == b['m']


def test_separated_features_are_recovered_by_their_nearest_kernels():
    r = TC.job_1d((np.array([-0.52, 0.02, 0.61]), np.array([0.8, -0.05, 0.9])))
    assert r['m'] > 0.8 * r['q'] and r['out'] == 'nearest' and r['rms'] < TC.RMS_TOL


def test_a_point_between_two_same_sign_features_intrudes():
    # two same-sign features 0.27 (2.25σ) apart, a starting point midway: the certificate there exceeds one and that kernel enters
    r = TC.job_1d((np.array([-0.135, 0.135, 0.7]), np.array([0.8, 0.8, 0.5])))
    assert r['m'] < -0.8 * r['q'] and r['out'] != 'nearest'
    L, K, c0 = TM.setup()
    assert r['argm'] == int(np.argmin(np.abs(c0)))                    # the kernel at 0, midway
    assert r['argm'] in r['order']


def test_skeleton_for_one_feature_is_the_single_feature_margin():
    """K = 1: the skeleton enters the feature's kernel and its margin is 1 − max_{j≠e}⟨φφ_j⟩/⟨φφ_e⟩, the certificate's and the
    single-feature margin alike."""
    L, K, c0 = TM.setup()
    cs, As = np.array([0.137]), np.array([-0.7])
    cert = TC.certificate(K, c0, cs, As)
    sk = TC.skeleton(K, c0, cs, As)
    assert sk['ok'] and sk['entrants'] == cert['e']
    assert abs(sk['margin'] - cert['m']) < 1e-12 and abs(cert['m_vor'] - cert['m']) < 1e-12 and abs(cert['m_rel'] - cert['m']) < 1e-12


def test_skeleton_follows_the_amplitude_ratio_where_the_certificate_cannot():
    """Round 43's counterexample: two opposite-sign features on starting points (10⁻⁶ off), 0.3 apart. The amplitude-free certificate is
    the same at every ratio (m = 1.057q); the path is exact with E at equal amplitudes and not when the second is halved (the kernel at 0.3
    slides onto the feature at 0.2) — and the skeleton, which follows the order of the entries and the screening, says the same."""
    L, K, c0 = TM.setup()
    cs = np.array([-0.1, 0.2]) + 1e-6
    out = {}
    for a2 in (-1.0, -0.5):
        As = np.array([1.0, a2])
        r = TC.job_1d((cs, As))
        sk = TC.skeleton(K, c0, cs, As)
        assert abs(r['m'] / r['q'] - 1.057) < 2e-3
        assert (sk['margin'] > 0) == (r['out'] == 'nearest') == sk['ok']
        out[a2] = r['out']
    assert out == {-1.0: 'nearest', -0.5: 'other'}


def test_fitted_skeleton_resolves_the_linear_skeletons_near_tie():
    """Part D(a)'s lattice pair 0.3 apart, same signs, the amplitude ratio 0.02 (and 50): the path is exact with E, the linear skeleton
    lets the intruder in by a near tie (m_sk = −0.0030), the fitted skeleton — each block at its own least-squares fit — the rightful
    kernel by +10⁻⁴: the linear skeleton's misses are its linearisation."""
    L, K, c0 = TM.setup()
    cs = np.array([-0.1, 0.2]) + TC.LAT_OFF
    for As in (np.array([1.0, 0.02]), np.array([1.0, 50.0])):
        lin, fit = TC.skeleton(K, c0, cs, As), TC.skeleton(K, c0, cs, As, fit=True)
        assert not lin['ok'] and -4e-3 < lin['margin'] < -2e-3
        assert fit['ok'] and 0 < fit['margin'] < 1e-3 and sorted(fit['entrants']) == [9, 12]
    assert TC.job_1d((cs, np.array([1.0, 0.02])))['out'] == 'nearest'


def test_relative_margin_compares_with_the_entrants_own_values():
    """m_rel = min_k s_kη(c_{e_k}) − max_{j∉E}|η(c_j)|: for an opposite-sign pair η(c_{e_k}) < 1, so m_rel < m."""
    L, K, c0 = TM.setup()
    cs, As = np.array([-0.1, 0.2]) + 1e-6, np.array([1.0, -1.0])
    c = TC.certificate(K, c0, cs, As)
    off = [j for j in range(c0.size) if j not in c['e']]
    own = min(np.sign(As[k]) * c['eta'][c['e'][k]] for k in range(2))
    assert own < 1 and abs(c['m_rel'] - (own - np.abs(c['eta'][off]).max())) < 1e-14 and c['m_rel'] < c['m']


def test_outcome_needs_the_nearest_kernels_and_no_exits():
    ev = lambda *xs: [(float(i), j, k, 0.0) for i, (j, k) in enumerate(xs)]
    act = np.zeros(21, bool); act[[3, 9]] = True
    hist = [(0.0, None, None, 1.0), (1.0, None, None, 1e-12)]
    assert TC.outcome(dict(events=ev((3, 'enter'), (9, 'enter')), act=act, hist=hist), [3, 9], 2) == 'nearest'
    assert TC.outcome(dict(events=ev((3, 'enter'), (5, 'enter'), (5, 'exit'), (9, 'enter')), act=act, hist=hist), [3, 9], 2) == 'other'
    assert TC.outcome(dict(events=ev((3, 'enter'), (9, 'enter')), act=act, hist=hist[:1]), [3, 9], 2) == 'not'


def test_bounded_polish_keeps_a_seed_level_entrant_from_jumping():
    """Round 43's target 267 of Part A with a short fast flow (t_fast = 10³): unbounded, the polish carries the second entrant, still at
    its seed coefficient, 6.7σ onto another feature and the path ends with other kernels; bounded (a near-seed kernel may move by one
    kernel width at most), it ends exactly with E, as with the default fast flow."""
    rng = np.random.default_rng(TC.SEED_A)
    tg = [(cs, s * a) for dmin in TC.DMIN_A for cs, s, a, u in TC.targets(rng, TC.N_A, dmin)]
    cs, As = tg[267]
    L, K, c0 = TM.setup()
    cert = TC.certificate(K, c0, cs, As)
    f = TC.X.bumps_target(K, cs, As)
    R_free = TC.RT.sliding_iss(L, K, f, c0, np.ones(c0.size), polish=True, t_fast=1e3, polish_step=np.inf)
    R_bnd = TC.RT.sliding_iss(L, K, f, c0, np.ones(c0.size), polish=True, t_fast=1e3)
    assert TC.outcome(R_free, cert['e'], 3) == 'other' and TC.outcome(R_bnd, cert['e'], 3) == 'nearest'
    assert R_bnd['polish_rejected'] > 0


def test_band_counts_the_targets_outside_it():
    rows = [dict(m=0.5, q=0.1, out='nearest'), dict(m=-0.5, q=0.1, out='other'), dict(m=0.01, q=0.1, out='other'),
            dict(m=0.2, q=0.1, out='other')]
    txt = TC.band(rows)
    assert 'C = 0.5: 2 of 3' in txt and 'C = 1: 2 of 3' in txt




def test_fitted_skeleton_fits_negative_features():
    """D(a)'s opposite-sign pair 0.3 apart at the ratio 1 (exact with E): the fitted skeleton's boxes follow the features' signs — with every
    box positive (referee round 46's mutant K7) the negative feature could not be fitted."""
    L, K, c0 = TM.setup()
    cs, As = np.array([-0.1, 0.2]) + TC.LAT_OFF, np.array([1.0, -1.0])
    fit = TC.skeleton(K, c0, cs, As, fit=True)
    assert fit['ok'] and fit['margin'] > 0 and sorted(fit['entrants']) == sorted(TC.entrants(K, c0, cs)[0])
    assert TC.job_1d((cs, As))['out'] == 'nearest'


def test_atoms_baseline_drops_the_tangent_directions():
    """The atoms-only baseline subtracts the fitted features' atoms and nothing else: for one feature its residual after the entry is zero,
    and on D(a)'s lattice pair at the ratio 0.02 its margin differs from the linear skeleton's (which also removes the tangent directions)."""
    L, K, c0 = TM.setup()
    a = TC.skeleton(K, c0, np.array([0.137]), np.array([-0.7]), atoms=True)
    assert a['ok'] and a['margin'] > 0
    cs, As = np.array([-0.1, 0.2]) + TC.LAT_OFF, np.array([1.0, 0.02])
    assert abs(TC.skeleton(K, c0, cs, As, atoms=True)['margin'] - TC.skeleton(K, c0, cs, As)['margin']) > 1e-4


def test_entrants_are_the_largest_correlations_near_the_window_edge():
    """Near the window's edge the samples' largest correlation is not the nearest kernel (referee round 46): a feature at 0.97 is entered by
    the kernel at 0.9 (index 19), which E must name for the path to count as exact with E."""
    L, K, c0 = TM.setup()
    cs, As = np.array([0.97]), np.array([0.8])
    e = TC.entrants(K, c0, cs)[0]
    assert e == [19] and int(np.argmin(np.abs(c0 - 0.97))) == 20
    assert TC.job_1d((cs, As))['out'] == 'nearest'


def test_ascent_target_ends_at_the_feature():
    """`ascent_target`: a seed of the feature's sign started 0.08 away climbs its correlation with the feature to the feature's centre
    (within the step σ/40); of the opposite sign it climbs away."""
    L, K, c0 = TM.setup()
    r = TC.X.bumps_target(K, np.array([0.137]), np.array([0.7]))
    c = TC.RT.ascent_target(K, r, np.array([0.137 + 0.08]), 1.0)
    assert abs(float(c[0]) - 0.137) < K.sig / 40 + 1e-12
    assert abs(float(TC.RT.ascent_target(K, r, np.array([0.137 + 0.08]), -1.0)[0]) - 0.137) > 0.08


def test_outcome_tolerances():
    """`outcome`: exact with E needs the residual below RMS_TOL, the active set E, no exits and exactly K events."""
    act = np.zeros(21, bool); act[[3, 9]] = True
    ev = [(0.0, 3, 'enter', 0.0), (1.0, 9, 'enter', 0.0)]
    R = lambda rms, events=ev: dict(events=events, act=act, hist=[(0.0, None, None, 1.0), (1.0, None, None, rms)])
    assert TC.outcome(R(0.5 * TC.RMS_TOL), [3, 9], 2) == 'nearest'
    assert TC.outcome(R(2 * TC.RMS_TOL), [3, 9], 2) == 'not'
    assert TC.outcome(R(0.5 * TC.RMS_TOL, ev + [(2.0, 5, 'enter', 0.0)]), [3, 9], 2) == 'other'
    assert TC.outcome(R(0.5 * TC.RMS_TOL), [3, 8], 2) == 'other'


def test_skeletons_on_the_targets_that_break_their_mutants():
    """Round 46's surviving mutants, each on a stored target it gets wrong: the linear skeleton's tangent directions at the features (not
    the starting points; B's configuration 31 at r = 3, exact with E); the fitted skeleton's boxes with the features' signs (D(a)'s
    opposite-sign pair 0.5 apart at the ratio 50) and its fit started at the features (B's configuration 11 at r = 3)."""
    L, K, c0 = TM.setup()
    cs = np.array([0.628354315580911, -0.04456980956766876, -0.6300330902222451])
    As = np.array([-0.5978787114809792, -0.31619518893970555, -0.0021571116965809])
    sk = TC.skeleton(K, c0, cs, As)
    assert sk['ok'] and sk['margin'] > 0.01                            # tangents at the starting points: −0.97
    assert TC.job_1d((cs, As))['out'] == 'nearest'
    fit = TC.skeleton(K, c0, np.array([-0.1, 0.4]) + TC.LAT_OFF, np.array([1.0, -50.0]), fit=True)
    assert fit['ok'] and fit['margin'] > 0.1                           # every box positive: −0.99
    cs = np.array([-0.26146688489923686, 0.6985887704301363, -0.6172569386174308])
    As = np.array([0.9608083993684658, 0.001283451133584242, 0.0018132196297868347])
    fit = TC.skeleton(K, c0, cs, As, fit=True)
    assert fit['ok'] and fit['margin'] > 0.01                          # started at the starting points: −0.34
    assert TC.job_1d((cs, As))['out'] == 'nearest'


def test_outcome_counts_kernels_and_its_tolerance_is_absolute():
    """`outcome`: the residual's tolerance is 10⁻⁸ (absolute), and an exact fit with more kernels than features is neither exact with E
    nor with other kernels (Part F's four-kernel fit)."""
    act = np.zeros(21, bool); act[[3, 9]] = True
    ev = [(0.0, 3, 'enter', 0.0), (1.0, 9, 'enter', 0.0)]
    R = lambda rms, a=act, events=ev: dict(events=events, act=a, hist=[(0.0, None, None, 1.0), (1.0, None, None, rms)])
    assert TC.outcome(R(5e-9), [3, 9], 2) == 'nearest' and TC.outcome(R(2e-8), [3, 9], 2) == 'not'
    a3 = act.copy(); a3[5] = True
    assert TC.outcome(R(1e-12, a3, ev + [(2.0, 5, 'enter', 0.0)]), [3, 9], 2) == 'not'


def test_the_polish_safeguard_on_stored_paths():
    """Round 46's surviving safeguard mutants, on stored paths (both `sliding_iss` and `sliding_iss_nd`, the experiment's): A's target 114
    accepts a near-seed kernel's polish 3.8σ from its start because it lands at its correlation's ascent (bounded by the start only, by
    the ascent of the wrong sign or by a quarter width it is rejected three times); A's target 12 and B's configuration 28 at r = 3 reject
    2 and 4 polishes, each retried (B's configuration 14 at r = 2 ends exact with E only with the near-seed scale relative to the block)."""
    L, K, c0 = TM.setup()
    T = {'A114': ([0.6541004699596475, 0.14717182832824738, -0.1741542573559407],
                  [-0.9993392739038046, -0.5218652653923312, -0.8701696976682551], 'other', 0),
         'A12': ([-0.3482759573852491, -0.0002953631766601239, -0.6856589204161185],
                 [0.8834332947114262, 0.9508709095889194, 0.40545651538789224], 'other', 2),
         'B115': ([0.23946960475117063, -0.19182582990096742, -0.6226746645823192],
                  [-0.017550238780768175, 0.0041097263697728054, 0.0032936667939692126], 'nearest', 4),
         'B58': ([-0.05502623411688401, -0.4293485291166097, 0.5138963933303842],
                 [0.044839474500752985, 0.08864282863781996, 0.03386534213286307], 'nearest', 1)}
    for name, (cs, As, out, rej) in T.items():
        cs, As = np.array(cs), np.array(As)
        e = TC.entrants(K, c0, cs)[0]
        R1 = TC.RT.sliding_iss(L, K, TC.X.bumps_target(K, cs, As), c0, np.ones(c0.size), polish=True)
        Rn = TC.job_1d((cs, As))
        assert (TC.outcome(R1, e, 3), int(R1['polish_rejected'])) == (out, rej), name
        assert (Rn['out'], Rn['rej']) == (out, rej), name


def test_eta_v_interpolates_signs_with_vanishing_derivative():
    """Part H (referee round 140): η_V = Σ_sφ_s(·)q_s with q = Γ(ΓᵀΓ)⁻¹(sign A, 0) takes the value sign A_k at each feature with zero
    derivative there (on the samples); for two separated features it is nondegenerate, its margin on the starting points off E positive."""
    L, K0, c0 = TM.setup()
    K = TC._k1(K0)
    cs, As = np.array([[-0.3], [0.35]]), np.array([0.8, -0.5])
    Gm = np.c_[K.phi(cs), K.dphi(cs).reshape(K.N, -1)]
    q = Gm @ np.linalg.solve(Gm.T @ Gm, np.r_[np.sign(As), 0.0, 0.0])
    assert np.allclose(K.phi(cs).T @ q, np.sign(As), rtol=0, atol=1e-10)
    assert np.allclose(K.dphi(cs).reshape(K.N, -1).T @ q, 0.0, rtol=0, atol=1e-8)
    e = TC.entrants(K0, c0, cs[:, 0])[0]
    ev = TC.eta_v(K, cs, As, c0, e)
    assert ev['nondeg'] and ev['hess'] < 0 and 0 < ev['m_v'] < 1 and ev['sup1'] < 1


def test_basis_pursuit_certificate_and_omp_on_one_feature():
    """One feature: basis pursuit on the grid recovers it (its least certificate t < 1), and continuous OMP finds it in one step, exact."""
    L, K0, c0 = TM.setup()
    K = TC._k1(K0)
    cs, As = np.array([[0.1234]]), np.array([0.7])
    t = TC.bp_certificate(K, cs, As)
    assert np.isfinite(t) and t < 1
    f = K.phi(cs) @ As
    rms, c, w = TC.omp(K, f, 1, 1)
    assert rms < TC.RMS_TOL and abs(c[0, 0] - 0.1234) < 1e-6 and abs(w[0] - 0.7) < 1e-6


def test_omp_needs_the_slide_on_overlapping_gaussians():
    """Two features 2.5σ apart: continuous OMP selects the correlation's peaks, shifted off the features by the overlap, and is not exact
    after two steps (the Gaussian kernel is not admissible in Elvira et al.'s sense); with the joint slide after each step it is exact."""
    L, K0, c0 = TM.setup()
    K = TC._k1(K0)
    cs, As = np.array([[-0.15], [0.15]]), np.array([0.9, 0.6])
    f = K.phi(cs) @ As
    rms, c, w = TC.omp(K, f, 2, 1)
    assert rms > 1e-4
    rms_s, c_s, w_s = TC.omp(K, f, 2, 1, slide=True)
    assert rms_s < TC.RMS_TOL
    assert np.allclose(np.sort(c_s[:, 0]), [-0.15, 0.15], rtol=0, atol=1e-6)


RAW = os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos', 'transit_certificate_raw.npz')


def _raw():
    import exp_cohesion_leak as E
    return E.load(RAW)


class _FakePool:
    """multiprocessing.Pool's stand-in: `map` returns each job's function name and arguments (the job lists of Parts G and H)."""
    def __init__(self, *a, **k):
        pass

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False

    def map(self, fn, jobs, chunksize=1):
        return [dict(fn=fn.__name__, args=j) for j in jobs]


def test_parts_G_and_H_job_lists_are_the_stored_rows():
    """Parts G and H's jobs (mutation testing after round 140): G runs `job_d` with the flow on B and D(a), D(c) — not D(b), whose rows
    have it — and `job_2d_flow` on C; H runs `h_job` on every target of A, B, D in one dimension and C in the plane; both in the stored
    rows' order."""
    import multiprocessing
    raw = _raw()
    saved = multiprocessing.Pool
    try:
        multiprocessing.Pool = _FakePool
        g = TC.part_G(raw['B'], raw['C'], raw['D'])
        h = TC.part_H(raw['A'], raw['B'], raw['C'], raw['D'])
    finally:
        multiprocessing.Pool = saved
    assert all(x['fn'] == 'job_d' and x['args'][3] == TC.T_MAX and x['args'][4] is True for x in g['rows'])
    assert [x['args'][0] for x in g['rows']] == [tuple(r['tag']) for r in raw['G']['rows']]
    assert all(raw['D'][i]['tag'][0] in 'ac' for p, i in (x['args'][0] for x in g['rows']) if p == 'D')
    assert all(x['fn'] == 'job_2d_flow' for x in g['plane'])
    assert [tuple(x['src']) for x in g['plane']] == [tuple(r['src']) for r in raw['G']['plane']]
    assert all(x['fn'] == 'h_job' for x in h)
    assert [x['args'][0] for x in h] == [tuple(r['src']) for r in raw['H']]
    assert all(x['args'][3] == (2 if x['args'][0][0] == 'C' else 1) for x in h)


def test_h_job_reproduces_stored_baselines():
    """Part H's baselines rerun on stored targets of every family (two of D(b)'s twenty on which η_V is degenerate, a plane target)
    against the raw: η_V's margin, sups and Hessian, basis pursuit's t, continuous OMP's residual and distance, with the slide."""
    raw = _raw()
    hb = {tuple(h['src']): h for h in raw['H']}
    dc = next(i for i, r in enumerate(raw['D']) if r['tag'][0] == 'c')
    keys = [('A', 0), ('A', 1), ('B', 0), ('B', 115), ('D', 0), ('D', 74), ('D', 75), ('D', dc), ('C', 0)]
    assert not hb[('D', 74)]['nondeg'] and not hb[('D', 75)]['nondeg']
    for key in keys:
        r = raw[key[0]][key[1]]
        h = TC.h_job((key, r['cs'], r['As'], 2 if key[0] == 'C' else 1))
        s = hb[key]
        assert tuple(h['src']) == key and h['nondeg'] == s['nondeg'], key
        for f in ('m_v', 'sup', 'sup1', 'hess', 'omp_rms', 'omp_dmax'):
            assert np.isclose(h[f], s[f], rtol=1e-6, atol=0), (key, f)
        assert h['somp_rms'] < 1e-12 and s['somp_rms'] < 1e-12
        assert (np.isnan(h['t_bp']) and np.isnan(s['t_bp'])) if key[0] == 'C' else np.isclose(h['t_bp'], s['t_bp'], rtol=1e-6, atol=0), key


def test_eta_v_nondegenerate_needs_the_hessian():
    """`eta_v`'s nondegeneracy needs both |η_V| < 1 off the features and s_k∇²η_V(c*_k) negative definite: with the second derivatives'
    sign planted reversed (the samples' `corr` wrapped) the same target is degenerate."""
    L, K0, c0 = TM.setup()
    K = TC._k1(K0)
    cs, As = np.array([[-0.3], [0.35]]), np.array([0.8, -0.5])
    e = TC.entrants(K0, c0, cs[:, 0])[0]
    ev = TC.eta_v(K, cs, As, c0, e)

    class Flip:
        def __init__(self, K):
            self.K = K

        def __getattr__(self, a):
            return getattr(self.K, a)

        def corr(self, *a, **k):
            out = list(self.K.corr(*a, **k))
            out[2] = -out[2]
            return tuple(out)
    fv = TC.eta_v(Flip(K), cs, As, c0, e)
    assert ev['nondeg'] and fv['sup'] == ev['sup'] and fv['hess'] > 0 and not fv['nondeg']


def test_job_2d_flow_reports_the_unpolished_flow():
    """Part G's plane job against its raw (C's first target), and with the unpolished flow planted off the path (its residual large, one
    block not relaxed): `out_flow` and `flow_relaxed` are the unpolished flow's, `out` the polished path's."""
    raw = _raw()
    r, s = raw['C'][0], raw['G']['plane'][0]
    x = TC.job_2d_flow((np.asarray(r['cs']), np.asarray(r['As'])))
    assert (x['out'], x['out_flow'], x['flow_relaxed'], list(x['e'])) == (s['out'], s['out_flow'], s['flow_relaxed'], list(s['e']))
    orig = TC.RN.sliding_iss_nd

    def planted(*a, **k):
        R = orig(*a, **k)
        if k.get('polish', True) is False:
            R = dict(R, relaxed=list(R['relaxed']) + [False], hist=list(R['hist'][:-1]) + [tuple(R['hist'][-1][:3]) + (1.0,)])
        return R
    try:
        TC.RN.sliding_iss_nd = planted
        y = TC.job_2d_flow((np.asarray(r['cs']), np.asarray(r['As'])))
    finally:
        TC.RN.sliding_iss_nd = orig
    assert y['out'] == s['out'] != 'not' and y['out_flow'] == 'not' and y['flow_relaxed'] is False


def test_job_d_with_the_flow_reproduces_stored_G_rows():
    """Part G's one-dimensional rows rerun (B's configuration 115, which rejects four polishes, and D(a)'s first) against the raw."""
    raw = _raw()
    rows = {tuple(g['tag']): g for g in raw['G']['rows']}
    for key in (('B', 115), ('D', 0)):
        r, s = raw[key[0]][key[1]], rows[key]
        x = TC.job_d((key, np.asarray(r['cs']), np.asarray(r['As']), TC.T_MAX, True))
        for f in ('out', 'out_flow', 'flow_relaxed', 'rej', 'order', 'order_flow'):
            assert x[f] == s[f], (key, f)
        assert x['out'] == r['out']


def test_report_G_and_H_count_planted_rows():
    """The G and H sections on the raw with planted rows: two of G's B rows whose unpolished flow ends elsewhere and a third whose blocks
    are not all relaxed; a timeout's rej = −1 on a B row (not counted as a rejection); a D(a) row that timed out (H drops it)."""
    import copy
    raw = copy.deepcopy(_raw())
    G = raw['G']
    ix = {tuple(g['tag']): i for i, g in enumerate(G['rows'])}
    for i in (0, 1):
        G['rows'][ix[('B', i)]]['out_flow'] = 'not' if raw['B'][i]['out'] != 'not' else 'nearest'
    G['rows'][ix[('B', 2)]]['flow_relaxed'] = False
    raw['B'][3]['rej'] = -1
    raw['D'][0]['out'] = 'timeout'
    assert raw['D'][0]['tag'][0] == 'a' and raw['B'][3].get('rej', 0) == -1
    md = TC.report(*[raw.get(k) for k in 'ABCDEFGH'])
    g = md[md.index('## G.'):md.index('## H.')]
    h = md[md.index('## H.'):]
    assert '| B | 160 of 160 | 158 | 159 | 160 | 8 (4) |' in g
    assert 'on 1151 of 1154 targets, its blocks relaxed on 1153.' in g             # the timed-out D(a) row's flow differs too
    assert '\n| D(a) | 53 |' in h and '\n| all | 1153 |' in h


def test_omp_selects_by_the_normalised_correlation():
    """Continuous OMP's selection is |⟨rφ(c)⟩|/‖φ(c)‖ (its maximum at a feature by Cauchy–Schwarz), not |⟨rφ(c)⟩|: two separated
    features, one near the samples' edge where ‖φ(c)‖ is smaller, with amplitudes A₂/A₁ = (‖φ(c₁)‖/‖φ(c₂)‖)^{3/2} between the two rules'
    break-even ratios — the first step selects the edge feature (unnormalised it would take the inner one; mutation testing after round
    140)."""
    L, K0, c0 = TM.setup()
    K = TC._k1(K0)
    cs = np.array([[0.0], [0.97]])
    n0, n1 = np.sqrt((K.phi(cs) ** 2).sum(0))
    assert n1 < 0.95 * n0
    As = np.array([1.0, (n0 / n1) ** 1.5])
    f = K.phi(cs) @ As
    rms, c, w = TC.omp(K, f, 1, 1)
    assert abs(c[0, 0] - 0.97) < 1e-4


# ---- Parts I (BLASSO at a finite regularisation) and J (basis pursuit by a jet certificate), added after round 140 -----------------------

MD = os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos', 'transit_certificate.md')


def _lasso_brute(G, b, lam):
    """The LASSO by enumeration of supports and signs: the point satisfying the KKT conditions (unique for G positive definite)."""
    from itertools import product
    n = len(b)
    for sg in product((-1, 0, 1), repeat=n):
        sg = np.array(sg, float)
        S = np.flatnonzero(sg)
        x = np.zeros(n)
        if len(S):
            x[S] = np.linalg.solve(G[np.ix_(S, S)], b[S] - lam * sg[S])
            if np.any(np.sign(x[S]) != sg[S]):
                continue
        c = b - G @ x
        off = sg == 0
        if np.all(np.abs(c[off]) <= lam * (1 + 1e-12)):
            return x
    raise AssertionError('no KKT point')


def test_lasso_fs_is_the_lasso():
    """Feature-sign search against the soft threshold (G = I) and against enumeration on random positive definite G (several λ, warm and
    cold starts); its reported KKT violation is the true one."""
    b = np.array([0.3, -1.2, 0.05, 0.8])
    x, kkt, it = TC.lasso_fs(np.eye(4), b, 0.25)
    assert np.allclose(x, np.sign(b) * np.maximum(np.abs(b) - 0.25, 0), atol=1e-15) and kkt < 1e-12
    rng = np.random.default_rng(7)
    for trial in range(12):
        A = rng.normal(size=(9, 7)) + 2.5 * rng.normal(size=(9, 1))        # correlated columns
        G, b = A.T @ A / 9, A.T @ rng.normal(size=9) / 9
        x0 = None
        for l in (0.3, 0.1, 0.02):
            lam = l * np.abs(b).max()
            x, kkt, it = TC.lasso_fs(G, b, lam, x0)
            xb = _lasso_brute(G, b, lam)
            assert np.allclose(x, xb, atol=1e-10), (trial, l)
            c = b - G @ x
            S = np.flatnonzero(x)
            off = np.setdiff1d(np.arange(len(b)), S)
            true = max(np.abs(c[off]).max() / lam - 1 if len(off) else -1, np.abs(c[S] - lam * np.sign(x[S])).max() / lam if len(S) else 0)
            assert abs(kkt - true) < 1e-14 and kkt < 1e-10
            x0 = x
            xc, _, _ = TC.lasso_fs(G, b, lam)
            assert np.allclose(xc, xb, atol=1e-10)


def test_blasso_verdict_plants():
    """Recovered: every nonzero within the tolerance of a feature and each feature's cluster of its amplitude's sign; exact: the support
    the features; lost: an empty cluster; the amplitude error the clusters' sums against the amplitudes."""
    cs, As = np.array([[-0.4], [0.1], [0.5]]), np.array([0.8, -0.5, 0.3])
    v = TC.blasso_verdict(dict(pos=cs.copy(), x=As.copy()), cs, As, 0.01)
    assert v['rec'] and v['exact'] and v['lost'] == 0 and v['far'] == 0 and v['err'] == 0.0 and v['n'] == 3
    pos = np.r_[cs, [[0.108]]]                                          # a neighbour within the tolerance shares the second feature
    v = TC.blasso_verdict(dict(pos=pos, x=np.array([0.8, -0.3, 0.3, -0.15])), cs, As, 0.01)
    assert v['rec'] and not v['exact'] and abs(v['err'] - 0.1) < 1e-12 and v['n'] == 4
    v = TC.blasso_verdict(dict(pos=np.r_[cs, [[0.111]]], x=np.array([0.8, -0.3, 0.3, -0.15])), cs, As, 0.01)
    assert not v['rec'] and v['far'] == 1 and v['lost'] == 0
    v = TC.blasso_verdict(dict(pos=cs[:2], x=As[:2]), cs, As, 0.01)                 # the third feature lost
    assert not v['rec'] and v['lost'] == 1 and v['far'] == 0 and abs(v['err'] - 1.0) < 1e-12
    v = TC.blasso_verdict(dict(pos=cs, x=np.array([0.8, 0.5, 0.3])), cs, As, 0.01)  # the second cluster of the wrong sign
    assert not v['rec'] and v['lost'] == 0 and abs(v['err'] - 2.0) < 1e-12
    v = TC.blasso_verdict(dict(pos=np.r_[cs, [[0.105]]], x=np.array([0.8, -0.5, 0.3, 1e-9])), cs, As, 0.01)
    assert v['rec'] and not v['exact']                                  # a tiny same-cluster neighbour of the other sign: the sum decides
    v = TC.blasso_verdict(dict(pos=np.zeros((0, 1)), x=np.zeros(0)), cs, As, 0.01)
    assert not v['rec'] and v['lost'] == 3 and v['n'] == 0


def test_bp_grid_replaces_the_features_neighbours():
    cs = np.array([[0.1012], [-0.3]])
    g = TC.bp_grid(TC._k1(TM.setup()[1]), cs, 0.005)
    assert np.allclose(g[-2:], cs) and len(g) == len(TC._grid(1, 0.005)) - 2 + 2
    assert np.min(np.abs(g[:-2] - cs.T)) > 0.0025


def test_parts_I_and_J_job_lists_are_the_stored_rows():
    raw = _raw()
    ij = TC.i_jobs(raw['A'], raw['B'], raw['D'])
    jj = TC.j_jobs(raw['A'], raw['B'], raw['C'], raw['D'])
    assert [j[0] for j in ij] == [tuple(r['src']) for r in raw['I']] and len(ij) == 1034
    assert [j[0] for j in jj] == [tuple(r['src']) for r in raw['J']] and len(jj) == 1154
    assert all(j[3] == (2 if j[0][0] == 'C' else 1) for j in jj)
    assert [tuple(r['src']) for r in raw['H'] if r['src'][0] != 'C'] == [j[0] for j in ij]


def test_run_chunks_resumes_and_keeps_the_job_order():
    import multiprocessing

    class Pool1(_FakePool):
        def map(self, fn, jobs, chunksize=1):
            return [fn(j) for j in jobs]
    calls, saves = [], []
    fn = lambda j: calls.append(j[0]) or dict(src=j[0], v=j[1])
    jobs = [(('A', i), 10 * i) for i in range(7)]
    done = [dict(src=('A', 4), v=40), dict(src=('A', 1), v=10)]
    saved, chunk = multiprocessing.Pool, TC.I_CHUNK
    try:
        multiprocessing.Pool = Pool1
        out = TC.run_chunks(fn, jobs, done, lambda: saves.append(len(done)), chunk=2)
    finally:
        multiprocessing.Pool = saved
    assert calls == [('A', i) for i in (0, 2, 3, 5, 6)] and saves == [4, 6, 7]
    assert out is done and [r['src'] for r in out] == [j[0] for j in jobs] and [r['v'] for r in out] == [10 * i for i in range(7)]


def test_i_job_reproduces_stored_blasso():
    """Part I rerun on stored targets: one of A, three of B, two of D(a) (one at the ratio 10⁻³) and a D(b) pair on which η_V is
    degenerate — the λ_max, the supports' positions and coefficients, the iterations, the KKT violations below 10⁻¹¹."""
    raw = _raw()
    ib = {tuple(r['src']): r for r in raw['I']}
    hb = {tuple(h['src']): h for h in raw['H']}
    db = next(k for k in ib if k[0] == 'D' and not hb[k]['nondeg'])
    da = next(i for i, r in enumerate(raw['D']) if r['tag'][0] == 'a' and not r['skip'] and abs(r['tag'][3] - 1e-3) < 1e-12)
    for k in [('A', 0), ('B', 13), ('B', 83), ('B', 115), ('D', 21), ('D', da), db]:   # B 13, B 83, D 21: a full step flips a sign
        r = raw[k[0]][k[1]]
        x = TC.i_job((k, r['cs'], r['As']))
        s = ib[k]
        assert abs(x['lmax'] - s['lmax']) < 1e-14
        for a, b in zip(x['sols'], s['sols']):
            assert a['l'] == b['l'] and np.allclose(a['pos'], b['pos'], atol=0) and np.allclose(a['x'], b['x'], rtol=1e-9, atol=1e-13)
            assert a['it'] == b['it']                                   # the warm start's iteration counts (printed in the md)
            assert a['kkt'] < 1e-11
    assert len(raw['I'][0]['sols']) == len(TC.I_LAMS)


def test_part_I_printed_rows_recount_from_the_raw():
    """The md's row "all | 0.0001" and the KKT sentence recounted from the raw with an independent recovery rule (every coefficient within
    two spacings of a feature, each feature's cluster sum of its amplitude's sign)."""
    raw = _raw()
    md = open(MD, encoding='utf-8').read()
    tol = 2 * 0.005
    for n_l, l in enumerate(TC.I_LAMS):
        rec = ex = lost = 0
        errs = []
        for r in raw['I']:
            t = raw[r['src'][0]][r['src'][1]]
            cs, As = np.asarray(t['cs'], float).ravel(), np.asarray(t['As'], float)
            s = r['sols'][n_l]
            p, x = np.asarray(s['pos'], float).ravel(), np.asarray(s['x'], float)
            d = np.abs(p[:, None] - cs[None, :])
            near, ok = d.argmin(1), d.min(1) <= tol
            sums = np.array([x[(near == k) & ok].sum() for k in range(len(As))])
            cnt = np.array([np.sum((near == k) & ok) for k in range(len(As))])
            good = ok.all() and (cnt > 0).all() and (np.sign(sums) == np.sign(As)).all()
            rec += good
            ex += len(x) == len(As) and bool(np.all(d.min(1) < 1e-12))
            lost += (cnt == 0).any()
            if good:
                errs.append(np.max(np.abs(sums - As) / np.abs(As)))
        import exp_cohesion_leak as E
        assert f"| all | {l:g} | 1034 | {rec} | {ex} | {lost} | {E.fmt(float(np.median(errs)))} [{E.fmt(max(errs))}] |" in md
    kkt = max(s['kkt'] for r in raw['I'] for s in r['sols'])
    assert kkt < 1e-11 and f"over all 4136 solves {E.fmt(kkt)} relative to λ" in md


def test_jets_and_the_jet_certificate():
    """`jets` are φ, ∇φ and ∇²φ at the features (the second derivative by finite differences); the jet certificate bounds Part H's t from
    above and is bounded by η_V's sup on the grid; for one feature in one dimension its p of the span of φ, φ′, φ″ beats η_V."""
    L, K0, c0 = TM.setup()
    K = TC._k1(K0)
    cs = np.array([[0.13], [-0.4]])
    J = TC.jets(K, cs)
    assert J.shape == (K.N, 6)
    h = 1e-4
    for k in range(2):
        e = np.zeros((2, 1))
        e[k] = h
        assert np.allclose(J[:, 2 + k], (K.phi(cs + e)[:, k] - K.phi(cs - e)[:, k]) / (2 * h), atol=1e-6)
        assert np.allclose(J[:, 4 + k], (K.phi(cs + e)[:, k] - 2 * K.phi(cs)[:, k] + K.phi(cs - e)[:, k]) / h ** 2, atol=1e-4)
    K2 = TC.TN.full()[0]
    J2 = TC.jets(K2, np.array([[0.1, -0.2]]))
    assert J2.shape == (K2.N, 6)
    raw = _raw()
    hb = {tuple(x['src']): x for x in raw['H']}
    jb = {tuple(x['src']): x for x in raw['J']}
    for k in [('A', 3), ('D', 74), ('C', 0)]:
        r = raw[k[0]][k[1]]
        dim = 2 if k[0] == 'C' else 1
        x = TC.j_job((k, r['cs'], r['As'], dim))
        assert x['dim'] == len(r['As']) * (6 if dim == 2 else 3) == jb[k]['dim']
        assert abs(x['t_jet'] - jb[k]['t_jet']) < 1e-9 and abs(x['t_etav'] - jb[k]['t_etav']) < 1e-12
        assert x['t_jet'] <= x['t_etav'] + 1e-9
        if dim == 1:
            assert x['t_jet'] >= hb[k]['t_bp'] - 1e-9
    one = TC.bp_jet_certificate(K, np.array([[0.0123]]), np.array([-1.0]))
    assert one['t_jet'] < one['t_etav'] - 1e-6 and one['t_jet'] < 1


def test_part_J_printed_rows_recount_from_the_raw():
    raw = _raw()
    md = open(MD, encoding='utf-8').read()
    pl = [r for r in raw['J'] if r['src'][0] == 'C']
    tj = [r['t_jet'] for r in pl]
    assert f"| C | 120 | {sum(t < 1 for t in tj)} | {np.median(tj):.5f} [{max(tj):.5f}] |" in md
    hb = {tuple(h['src']): h for h in raw['H']}
    one = [r for r in raw['J'] if r['src'][0] != 'C']
    agree = sum((r['t_jet'] < 1) == (hb[tuple(r['src'])]['t_bp'] < 1) for r in one)
    assert f"In one dimension ({len(one)} targets) t_J < 1 exactly where Part H's t < 1 on {agree}" in md
    assert f"certified exact on {sum(t < 1 for t in tj)} of {len(pl)} targets" in md


def test_report_I_and_J_on_planted_rows():
    """The I and J sections on planted rows: an I row whose smallest-λ solution gains a far atom (one fewer recovered, as path 'other'),
    a J plane row with t_J ≥ 1 (one fewer certified), a one-dimensional J row above 1 where Part H's t is below (one fewer agreement) and
    one where both are above 1 (still an agreement)."""
    import copy
    raw = copy.deepcopy(_raw())
    k = next(i for i, r in enumerate(raw['I']) if raw[r['src'][0]][r['src'][1]]['out'] == 'other'
             and TC.blasso_verdict(r['sols'][-1], raw[r['src'][0]][r['src'][1]]['cs'], raw[r['src'][0]][r['src'][1]]['As'], 0.01)['rec'])
    s = raw['I'][k]['sols'][-1]
    s['pos'] = np.r_[np.asarray(s['pos'], float).reshape(-1, 1), [[1.19]]]
    s['x'] = np.r_[s['x'], 1e-6]
    jc = next(i for i, r in enumerate(raw['J']) if r['src'][0] == 'C')
    raw['J'][jc]['t_jet'] = 1.0
    raw['J'][0]['t_jet'] = 1.5
    k1 = tuple(raw['J'][1]['src'])                                      # both above 1: an agreement that is not a certification
    next(h for h in raw['H'] if tuple(h['src']) == k1)['t_bp'] = 1.2
    raw['J'][1]['t_jet'] = 1.3
    md = TC.report(*[raw.get(x) for x in 'ABCDEFGHIJ'])
    i_ = md[md.index('## I.'):md.index('## J.')]
    j_ = md[md.index('## J.'):]
    assert '| all | 0.0001 | 1034 | 1011 |' in i_ and '680 of 680, 331 of 354 |' in i_
    assert '| C | 120 | 119 |' in j_ and 'certified exact on 119 of 120' in j_
    assert "t_J < 1 exactly where Part H's t < 1 on 1033" in j_


def test_k_classify_plants():
    """Part K's types on planted solutions: a nonzero beyond the tolerance but within K_RD·σ is displaced (and keeps its feature, so the
    feature is not lost although its cluster is empty), beyond K_RD·σ spurious; a feature without a nonzero within K_RD·σ is lost; the
    class is the first of lost, spurious, displaced, wrong sign; the distances and the lost features' relative amplitudes."""
    cs, As = np.array([[-0.4], [0.1], [0.5]]), np.array([0.8, -0.5, 0.04])
    tol, rd = 0.01, 0.06
    v = TC.k_classify(dict(pos=cs.copy(), x=As.copy()), cs, As, tol, rd)
    assert v['cls'] == 'rec' and v['lost'] == 0 and v['lost_tol'] == 0 and v['disp'] == 0 and v['spur'] == 0 and not v['sign']
    v = TC.k_classify(dict(pos=np.array([[-0.4], [0.13], [0.5]]), x=As.copy()), cs, As, tol, rd)   # the second feature's atom 0.03 off
    assert v['cls'] == 'displaced' and v['lost'] == 0 and v['lost_tol'] == 1 and v['disp'] == 1 and abs(v['ddisp'] - 0.03) < 1e-12
    v = TC.k_classify(dict(pos=np.array([[-0.4], [0.1], [0.58]]), x=As.copy()), cs, As, tol, rd)   # 0.08 off: lost and spurious
    assert v['cls'] == 'lost' and v['lost'] == 1 and v['spur'] == 1 and abs(v['dspur'] - 0.08) < 1e-12
    assert np.allclose(v['lost_rel'], [0.05], rtol=1e-12)
    v = TC.k_classify(dict(pos=np.r_[cs, [[-0.15]]], x=np.r_[As, 1e-4]), cs, As, tol, rd)          # an atom 0.25 from every feature
    assert v['cls'] == 'spurious' and v['lost'] == 0 and v['spur'] == 1 and np.isnan(v['ddisp']) and abs(v['dspur'] - 0.25) < 1e-12
    v = TC.k_classify(dict(pos=np.r_[cs, [[-0.15]], [[0.115]]], x=np.r_[As, 1e-4, -1e-4]), cs, As, tol, rd)
    assert v['cls'] == 'spurious' and v['disp'] == 1                                         # spurious before displaced
    v = TC.k_classify(dict(pos=cs.copy(), x=np.array([0.8, 0.5, 0.04])), cs, As, tol, rd)         # a cluster of the wrong sign
    assert v['cls'] == 'sign' and v['sign']
    v = TC.k_classify(dict(pos=cs[:2], x=As[:2]), cs, As, tol, rd)                              # the small feature dropped
    assert v['cls'] == 'lost' and v['lost'] == 1 and v['lost_tol'] == 1 and np.allclose(v['lost_rel'], [0.05], rtol=1e-12)
    v = TC.k_classify(dict(pos=np.array([[0.3]]), x=np.array([0.1])), cs, As, tol, rd)
    assert v['cls'] == 'lost' and v['lost'] == 3 and v['spur'] == 1


def test_part_K_reproduces_the_stored_rows():
    """Part K recomputed from the stored Parts A, B, D, G and I equals the raw's K, and its classes agree with Part I's verdict."""
    raw = _raw()
    k = TC.part_K(raw['A'], raw['B'], raw['D'], raw['G'], raw['I'])
    assert k['tol'] == raw['K']['tol'] == 0.01 and abs(k['rd'] - 0.06) < 1e-12
    assert len(k['rows']) == len(raw['K']['rows']) == 1034 and len(k['seq']) == 1034 and len(k['plane']) == 120
    for a, b in zip(k['rows'], raw['K']['rows']):
        assert a['src'] == tuple(b['src'])
        for x, y in zip(a['types'], b['types']):
            assert x['cls'] == y['cls'] and x['lost'] == y['lost'] and x['lost_tol'] == y['lost_tol']
            assert x['disp'] == y['disp'] and x['spur'] == y['spur'] and x['lost_rel'] == y['lost_rel']
    assert [x['same'] for x in k['seq']] == [x['same'] for x in raw['K']['seq']]
    for r, t in zip(raw['I'], k['rows']):
        p = raw[r['src'][0]][r['src'][1]]
        for s, x in zip(r['sols'], t['types']):
            assert (x['cls'] == 'rec') == TC.blasso_verdict(s, p['cs'], p['As'], 0.01)['rec']
    import copy                                                         # planted: a permuted sequence, a plane class that differs
    A2, G2 = copy.deepcopy(raw['A']), copy.deepcopy(raw['G'])
    A2[0]['order_flow'] = list(A2[0]['order'])[::-1]
    G2['plane'][0]['out_flow'] = 'not'
    k2 = TC.part_K(A2, raw['B'], raw['D'], G2, raw['I'][:1])
    assert k2['seq'][0]['same'] is False and sum(x['same'] for x in k2['seq']) == 1033 and k2['plane'][0]['same'] is False


def test_part_K_printed_rows_recount_from_the_raw():
    """The md's K rows "all" recounted from the stored solutions with an independent rule (a feature lost: no nonzero within 0.06 whose
    nearest feature it is; spurious: a nonzero beyond 0.06 of every feature; displaced: one between 0.01 and 0.06), and the entry
    sequences recounted from the Parts' rows."""
    raw = _raw()
    md = open(MD, encoding='utf-8').read()
    for n_l, l in enumerate(TC.I_LAMS):
        c = dict(bad=0, lost=0, spur=0, disp=0, empty=0)
        for r in raw['I']:
            t = raw[r['src'][0]][r['src'][1]]
            cs, As = np.asarray(t['cs'], float).ravel(), np.asarray(t['As'], float)
            s = r['sols'][n_l]
            p, x = np.asarray(s['pos'], float).ravel(), np.asarray(s['x'], float)
            d = np.abs(p[:, None] - cs[None, :])
            near, dm = d.argmin(1), d.min(1)
            ok = dm <= 0.01
            cnt = np.array([np.sum((near == k) & ok) for k in range(len(As))])
            sums = np.array([x[(near == k) & ok].sum() for k in range(len(As))])
            if ok.all() and (cnt > 0).all() and (np.sign(sums) == np.sign(As)).all():
                continue
            c['bad'] += 1
            c['empty'] += (cnt == 0).any()
            if any(not np.any((near == k) & (dm <= 0.06)) for k in range(len(As))):
                c['lost'] += 1
            elif (dm > 0.06).any():
                c['spur'] += 1
            elif (dm > 0.01).any():
                c['disp'] += 1
        assert f"| all | {l:g} | 1034 | {c['bad']} | {c['lost']} | {c['spur']} | {c['disp']} | 0 | {c['empty']} |" in md
    same = sum(list(r['order']) == list(r['order_flow']) for p in 'AD' for r in raw[p] if not r['skip'] and 'order_flow' in r)
    same += sum(list(raw[g['tag'][0]][g['tag'][1]]['order']) == list(g['order_flow']) for g in raw['G']['rows'])
    cls = sum(g['out'] == g['out_flow'] for g in raw['G']['plane'])
    own = sum(not r['skip'] and 'order_flow' in r for p in 'AD' for r in raw[p])
    assert own == 640 and f"in the same order as the polished path on {same} of the 1034 one-dimensional targets ({own} of A and D(b)" in md
    assert f"its outcome class is the polished path's on {cls} of 120." in md


def test_report_K_on_planted_rows():
    """The K section on planted rows: one A failure moved from displaced to spurious, a sequence that differs, a plane class that
    differs; the section starts after J's and the report's last line stays last."""
    import copy
    raw = copy.deepcopy(_raw())
    i = next(n for n, r in enumerate(raw['K']['rows']) if r['src'][0] == 'A' and r['types'][0]['cls'] == 'displaced')
    raw['K']['rows'][i]['types'][0].update(cls='spurious', dspur=0.3)
    raw['K']['seq'][0]['same'] = False
    raw['K']['plane'][0]['same'] = False
    md = TC.report(*[raw.get(x) for x in 'ABCDEFGHIJK'])
    k_ = md[md.index('## K.'):]
    assert md.index('## J.') < md.index('## K.') and md.rstrip().endswith('nothing is edited by hand.')
    assert '| A | 0.1 | 400 | 63 | 0 | 1 | 62 | 0 |' in k_ and '| all | 0.1 | 1034 | 366 | 159 | 1 | 206 | 0 | 306 |' in k_
    assert 'on 1033 of the 1034 one-dimensional targets (640 of A and D(b)' in k_ and "polished path's on 119 of 120" in k_
    assert TC.report(*[raw.get(x) for x in 'ABCDEFGHIJ']) == md[:md.index('## K.')] + md[md.rindex('Numbers are produced'):]


# ---- L and M (after round 149) --------------------------------------------------------------------------------------------------------

def _plane_lp(K, cs, As, h, cut):
    """Part L's linear program built column by column from `K.phi` (the grid's columns on the samples, p in the truncated span of the
    tensor basis written out as a full 1681 × dim matrix), for HiGHS: an independent assembly of the same program."""
    from scipy.optimize import linprog
    U1, s1, V1, ti, tj = TC._kron_basis(K, h, cut)
    n = len(K.x1)
    Ub = np.stack([np.outer(U1[:, i], U1[:, j]).ravel() / (s1[i] * s1[j]) for i, j in zip(ti, tj)], 1)   # p = Ub w
    g = TC._grid(2, h)
    g = g[np.min(np.linalg.norm(g[:, None, :] - cs[None, :, :], axis=2), axis=1) > 0.5 * h]
    Qg, Qs = K.phi(g).T @ Ub, K.phi(cs).T @ Ub
    m = Ub.shape[1]
    sol = linprog(np.r_[np.zeros(m), 1.0], A_ub=np.r_[np.c_[Qg, -np.ones(len(g))], np.c_[-Qg, -np.ones(len(g))]], b_ub=np.zeros(2 * len(g)),
                  A_eq=np.c_[Qs, np.zeros(len(As))], b_eq=np.sign(As), bounds=[(None, None)] * m + [(0, None)], method='highs')
    assert n * n == K.N and sol.status == 0
    return float(sol.x[-1]), m


def test_kron_basis_is_the_tensor_grids_span():
    """The pairs kept at Part H's truncation are 1094, the dimension Part J's md reports for the span of the 6561 grid columns; in the
    basis, the grid values of p = (u_i ⊗ u_j)/(s_is_j) are v_i ⊗ v_j (orthonormal); the coarser truncations keep fewer pairs."""
    K = TC.TN.full()[0]
    U1, s1, V1, ti, tj = TC._kron_basis(K, TC.H_BP[1])
    assert len(ti) == 1094 and V1.shape == (81, 41) and np.all(s1[ti] * s1[tj] > 1e-13 * s1[0] ** 2)
    assert np.allclose(V1.T @ V1, np.eye(41), atol=1e-12)
    g = TC._grid(2, TC.H_BP[1])
    for i, j in [(0, 0), (3, 7), (10, 5)]:
        p = np.outer(U1[:, i], U1[:, j]).ravel() / (s1[i] * s1[j])
        eta = K.phi(g).T @ p
        assert np.allclose(eta, np.outer(V1[:, i], V1[:, j]).ravel(), rtol=0, atol=1e-9)
    dims = [len(TC._kron_basis(K, TC.H_BP[1], c)[3]) for c in TC.L_CUTS]
    assert dims[0] == 1094 and dims[0] > dims[1] > dims[2]


def test_bp_certificate_kron_is_the_linear_program():
    """With 11 × 11 samples and a grid of spacing 0.1 (finer than the samples, as in Part L) the interior-point solution equals HiGHS on
    the same program assembled independently from `K.phi`,
    the bracket holds it, the certificate interpolates the signs at the features and its sup on the grid is t; a smaller span (a coarser
    truncation) has a larger least t."""
    xs = np.arange(-1, 1 + 1e-9, 0.2)
    K = TC.RN.BumpsGrid2(xs, xs, 0.25)
    cs, As = np.array([[0.31, 0.18], [-0.27, 0.42], [-0.47, -0.06]]), np.array([-0.45, 0.88, -0.73])
    out = {}
    for c in (1e-7, 1e-4):
        x = TC.bp_certificate_kron(K, cs, As, h=0.1, cut=c)
        t_h, m = _plane_lp(K, cs, As, 0.1, c)
        g = TC._grid(2, 0.1)
        assert x['dim'] == m and x['n_con'] == int(np.sum(np.min(np.linalg.norm(g[:, None, :] - cs[None], axis=2), axis=1) > 0.05))
        assert x['t_low'] - 1e-9 <= t_h <= x['t'] + 1e-9 and x['t'] - x['t_low'] < 1e-6 and abs(x['t'] - t_h) < 1e-6
        assert x['re'] < 1e-12 and x['feq'] < 1e-4 and abs(x['t_dir'] - x['t']) < 1e-3
        out[c] = x['t']
    assert out[1e-7] < out[1e-4]


def test_dual_lower_makes_the_dual_feasible():
    """`_dual_lower` on a small explicit program (B 40 × 6, two equalities, four grid rows dropped) with planted multipliers whose
    stationarity residual is large: the corrected multipliers are nonnegative and satisfy B_keepᵀ(a − b) + Cᵀν = 0, and the bound lies
    below the program's least t (HiGHS)."""
    from scipy.optimize import linprog
    rng = np.random.default_rng(7)
    B, C, s = rng.normal(size=(40, 6)), rng.normal(size=(2, 6)), np.array([1.0, -1.0])
    keep = np.ones(40, bool)
    keep[[3, 11, 20, 33]] = False
    Bk = B[keep]
    l1, l2 = rng.uniform(0, 0.05, 40) * keep, rng.uniform(0, 0.05, 40) * keep
    nu = rng.normal(size=2)
    low, a, b = TC._dual_lower(lambda q: B.T @ q, lambda w: B @ w, lambda v: np.linalg.solve(Bk.T @ Bk, v), keep, l1, l2, nu, C, s)
    assert np.all(a >= 0) and np.all(b >= 0) and np.all(a[~keep] == 0) and np.all(b[~keep] == 0)
    r0 = Bk.T @ (l1 - l2)[keep] + C.T @ nu
    assert np.abs(r0).max() > 0.1 and np.abs(Bk.T @ (a - b)[keep] + C.T @ nu).max() < 1e-12
    assert np.isclose(low, -(s @ nu) / (a.sum() + b.sum()), rtol=1e-15, atol=0)
    n = int(keep.sum())
    sol = linprog(np.r_[np.zeros(6), 1.0], A_ub=np.r_[np.c_[Bk, -np.ones(n)], np.c_[-Bk, -np.ones(n)]], b_ub=np.zeros(2 * n),
                  A_eq=np.c_[C, np.zeros(2)], b_eq=s, bounds=[(None, None)] * 6 + [(0, None)], method='highs')
    assert sol.status == 0 and low <= sol.fun + 1e-12

def test_bp_certificate_kron_below_the_jet_certificate():
    """One feature and a pair, on the plane's grid at Part H's truncation: the least certificate's bracket lies below Part J's jet
    certificate and η_V's sup (both certificates of the same program, up to the truncation) and below one."""
    K = TC.TN.full()[0]
    for cs, As in ((np.array([[0.0123, -0.031]]), np.array([1.0])), (np.array([[0.0123, -0.031], [0.3311, 0.2]]), np.array([1.0, -0.4]))):
        x = TC.bp_certificate_kron(K, cs, As)
        j = TC.bp_jet_certificate(K, cs, As)
        assert x['t_low'] <= x['t'] < j['t_jet'] <= j['t_etav'] + 1e-9 and x['t'] < 1 and x['t'] - x['t_low'] < 1e-5


def test_part_L_job_list_and_stored_rows():
    """L's jobs are the 120 plane targets of C; l_job on two stored targets reproduces the raw's rows (the iteration is deterministic)."""
    raw = _raw()
    jobs = TC.l_jobs(raw['C'])
    assert [j[0] for j in jobs] == [('C', i) for i, r in enumerate(raw['C']) if not r['skip']] and len(jobs) == 120
    lb = {tuple(r['src']): r for r in raw['L']}
    assert set(lb) == {j[0] for j in jobs}
    for j in (jobs[0], jobs[57]):
        x = TC.l_job(j)
        for a, b in zip(x['cuts'], lb[j[0]]['cuts']):
            assert a['cut'] == b['cut'] and a['it'] == b['it'] and a['dim'] == b['dim'] and a['nact'] == b['nact']
            assert abs(a['t'] - b['t']) < 1e-12 and abs(a['t_low'] - b['t_low']) < 1e-12


def test_part_L_printed_rows_recount_from_the_raw():
    """The md's L rows recounted from the raw: per truncation the targets with t < 1, t's median, those with t ≤ t_J, the widest bracket."""
    raw = _raw()
    md = open(MD, encoding='utf-8').read()
    jb = {tuple(r['src']): r['t_jet'] for r in raw['J']}
    for n_l, c in enumerate(TC.L_CUTS):
        t = [r['cuts'][n_l]['t'] for r in raw['L']]
        le = sum(r['cuts'][n_l]['t'] <= jb[tuple(r['src'])] for r in raw['L'])
        dim = raw['L'][0]['cuts'][n_l]['dim']
        assert f"| {c:g} | {dim} | 120 | {sum(v < 1 for v in t)} | {np.median(t):.5f} [{min(t):.5f}, {max(t):.5f}] | {le} |" in md
        w = max(r['cuts'][n_l]['t'] - r['cuts'][n_l]['t_low'] for r in raw['L'])
        assert w > 0 and w < 1e-4
    lo = [r['cuts'][0]['t'] < 1 for r in raw['L']]
    assert f"is exact (t < 1) on {sum(lo)} of 120 plane targets" in md
    rep = []                                                            # every feature within half a spacing of its nearest grid point
    for r in raw['L']:
        cs = np.asarray(raw['C'][r['src'][1]]['cs'], float).reshape(-1, 2)
        near = np.clip(np.round((cs + 1) / 0.025), 0, 80) * 0.025 - 1
        rep.append(bool(np.all(np.hypot(*(cs - near).T) <= 0.0125)))
    t_r = [r['cuts'][0]['t'] for r, x in zip(raw['L'], rep) if x]
    t_n = [r['cuts'][0]['t'] for r, x in zip(raw['L'], rep) if not x]
    assert (f"t has the median [range] {np.median(t_r):.5f} [{min(t_r):.5f}, {max(t_r):.5f}] ({len(t_r)} targets), where some feature "
            f"replaced none {np.median(t_n):.5f} [{min(t_n):.5f}, {max(t_n):.5f}] ({len(t_n)})") in md


def test_part_L_three_digit_values_are_the_raws_rounded():
    """The three-digit values that the paper and supplement quote (D3eu-2) are printed by the report and equal the raw's values rounded,
    recomputed here independently: t's range and medians, t_J − t's range, the grid-geometry and path-outcome medians."""
    raw = _raw()
    md = open(MD, encoding='utf-8').read()
    l_ = md[md.index('## L.'):md.index('## M.')]
    line = next(x for x in l_.splitlines() if x.startswith('At three digits'))
    t0 = np.array([r['cuts'][0]['t'] for r in raw['L']])
    jb = {tuple(r['src']): r['t_jet'] for r in raw['J']}
    g0 = np.array([jb[tuple(r['src'])] - r['cuts'][0]['t'] for r in raw['L']])
    med = [np.median([r['cuts'][n]['t'] for r in raw['L']]) for n in range(3)]
    rep = []
    for r in raw['L']:
        cs = np.asarray(raw['C'][r['src'][1]]['cs'], float).reshape(-1, 2)
        near = np.clip(np.round((cs + 1) / 0.025), 0, 80) * 0.025 - 1
        rep.append(bool(np.all(np.hypot(*(cs - near).T) <= 0.0125)))
    rep = np.array(rep)
    oth = np.array([raw['C'][r['src'][1]]['out'] == 'other' for r in raw['L']])
    want = [round(float(x), 3) for x in (t0.min(), np.median(t0), t0.max(), *med, g0.min(), g0.max(),
                                         np.median(t0[rep]), np.median(t0[~rep]), np.median(t0[~oth]), np.median(t0[oth]))]
    import re
    got = [float(x) for x in re.findall(r'(?<![\d.e-])(\d\.\d{3})(?!\d)', line)]
    assert got == want
    assert want == [0.567, 0.620, 0.785, 0.620, 0.696, 0.780, 0.205, 0.409, 0.602, 0.724, 0.712, 0.613]


def test_kron_gram_is_the_gram():
    """`KronGram` against the dense ⟨φ_iφ_j⟩ on the plane's samples: products with sparse and dense vectors, from either side, blocks."""
    K = TC.TN.full()[0]
    rng = np.random.default_rng(3)
    pos = rng.uniform(-1, 1, (40, 2))
    P = K.phi(pos)
    G0 = P.T @ P / K.N
    G = TC.KronGram(np.exp(-(K.x1[:, None] - pos[None, :, 0]) ** 2 / (2 * K.sig ** 2)),
                    np.exp(-(K.x2[:, None] - pos[None, :, 1]) ** 2 / (2 * K.sig ** 2)), K.N)
    x = np.zeros(40)
    x[[2, 9, 31]] = [0.5, -1.2, 2.0]
    y = rng.normal(size=40)
    assert np.allclose(G @ x, G0 @ x, rtol=0, atol=1e-14) and np.allclose(G @ y, G0 @ y, rtol=0, atol=1e-13)
    assert np.isclose(y @ G @ y, y @ G0 @ y, rtol=1e-13, atol=0) and np.isclose(0.5 * y @ G @ y, 0.5 * y @ G0 @ y, rtol=1e-13, atol=0)
    S = np.array([1, 4, 30])
    assert np.allclose(G[np.ix_(S, S)], G0[np.ix_(S, S)], rtol=0, atol=1e-15) and np.allclose(G[np.ix_(S, S[:2])], G0[np.ix_(S, S[:2])])


def test_blasso_kron_is_blasso_in_the_plane():
    """On a coarse plane grid (spacing 0.1) `blasso_kron` gives `blasso`'s solutions (the dense Gram matrix) at every λ."""
    K = TC.TN.full()[0]
    cs, As = np.array([[0.31, 0.18], [-0.27, 0.42], [-0.47, -0.06]]), np.array([-0.45, 0.88, -0.73])
    a, b = TC.blasso_kron(K, cs, As, h=0.1), TC.blasso(K, cs, As, h=0.1)
    assert np.isclose(a['lmax'], b['lmax'], rtol=1e-13, atol=0)
    for x, y in zip(a['sols'], b['sols']):
        assert x['l'] == y['l'] and np.array_equal(x['pos'], y['pos']) and np.allclose(x['x'], y['x'], rtol=1e-9, atol=0)
        assert x['kkt'] < 1e-11 and y['kkt'] < 1e-11


def test_part_M_reproduces_stored_rows():
    """m_job on three stored plane targets reproduces the raw's M rows; the jobs are L's (the 120 plane targets)."""
    raw = _raw()
    mb = {tuple(r['src']): r for r in raw['M']}
    jobs = TC.l_jobs(raw['C'])
    assert set(mb) == {j[0] for j in jobs}
    for j in (jobs[0], jobs[40], jobs[119]):
        x = TC.m_job(j)
        assert np.isclose(x['lmax'], mb[j[0]]['lmax'], rtol=1e-14, atol=0)
        for a, b in zip(x['sols'], mb[j[0]]['sols']):
            assert np.array_equal(a['pos'], np.asarray(b['pos'])) and np.allclose(a['x'], b['x'], rtol=1e-10, atol=0) and a['it'] == b['it']


def test_part_M_printed_rows_recount_from_the_raw():
    """The md's M rows recounted from the stored solutions with an independent rule (tolerance 0.05, radius 0.075): recovered, the
    support exactly the features, a cluster empty, and the failures by type."""
    raw = _raw()
    md = open(MD, encoding='utf-8').read()
    for n_l, l in enumerate(TC.I_LAMS):
        c = dict(rec=0, exact=0, empty=0, lost=0, spur=0, disp=0, sign=0)
        for r in raw['M']:
            t = raw['C'][r['src'][1]]
            cs, As = np.asarray(t['cs'], float).reshape(-1, 2), np.asarray(t['As'], float)
            s = r['sols'][n_l]
            p, x = np.asarray(s['pos'], float).reshape(-1, 2), np.asarray(s['x'], float)
            d = np.linalg.norm(p[:, None, :] - cs[None], axis=2)
            near, dm = d.argmin(1), d.min(1)
            ok = dm <= 0.05
            cnt = np.array([np.sum((near == k) & ok) for k in range(3)])
            sums = np.array([x[(near == k) & ok].sum() for k in range(3)])
            c['exact'] += len(x) == 3 and bool(np.all(dm < 1e-12))
            c['empty'] += bool((cnt == 0).any())
            if ok.all() and (cnt > 0).all() and (np.sign(sums) == np.sign(As)).all():
                c['rec'] += 1
            elif any(not np.any((near == k) & (dm <= 0.075)) for k in range(3)):
                c['lost'] += 1
            elif (dm > 0.075).any():
                c['spur'] += 1
            elif (dm > 0.05).any():
                c['disp'] += 1
            else:
                c['sign'] += 1
        assert (f"| {l:g} | 120 | {c['rec']} | {c['exact']} | {c['empty']} | {c['lost']}, {c['spur']}, {c['disp']}, {c['sign']} |") in md
    kk = max(s['kkt'] for r in raw['M'] for s in r['sols'])
    assert kk < 1e-10 and f"over all 480 solves {TC.E.fmt(kk)} relative to λ" in md


def test_report_L_and_M_on_planted_rows():
    """The L and M sections on planted rows: an L row above one at Part H's truncation and above t_J, an M solution with a far atom
    (recovered one fewer, a spurious failure); the sections follow K's and the report without them is the old one."""
    import copy
    raw = copy.deepcopy(_raw())
    full = TC.report(*[raw.get(x) for x in 'ABCDEFGHIJKLM'])
    assert full.index('## K.') < full.index('## L.') < full.index('## M.') and full.rstrip().endswith('nothing is edited by hand.')
    assert TC.report(*[raw.get(x) for x in 'ABCDEFGHIJK']) == full[:full.index('## L.')] + full[full.rindex('Numbers are produced'):]
    l0 = full[full.index('## L.'):full.index('## M.')]
    n0 = sum(r['cuts'][0]['t'] < 1 for r in raw['L'])
    le0 = sum(r['cuts'][0]['t'] <= next(j['t_jet'] for j in raw['J'] if tuple(j['src']) == tuple(r['src'])) for r in raw['L'])
    assert f"| 1e-13 | 1094 | 120 | {n0} |" in l0
    mono = sum(r['cuts'][0]['t'] < r['cuts'][1]['t'] < r['cuts'][2]['t'] for r in raw['L'])                 # recounted, strictly
    assert f"(a smaller span, a larger least t) on {mono} of 120" in l0
    assert "radius 0.5σ = 0.075 (σ = 0.15)" in full and "every nonzero within 2 spacings (0.05) of a feature" in full
    raw['L'][0]['cuts'][0]['t'] = 1.5
    m = next(i for i, r in enumerate(raw['M']) if TC.blasso_verdict(r['sols'][-1], raw['C'][r['src'][1]]['cs'], raw['C'][r['src'][1]]['As'], 0.05)['rec'])
    s = raw['M'][m]['sols'][-1]
    s['pos'] = np.r_[np.asarray(s['pos'], float).reshape(-1, 2), [[0.95, 0.95]]]
    s['x'] = np.r_[s['x'], 1e-6]
    md = TC.report(*[raw.get(x) for x in 'ABCDEFGHIJKLM'])
    l_, m_ = md[md.index('## L.'):md.index('## M.')], md[md.index('## M.'):]
    assert f"| 1e-13 | 1094 | 120 | {n0 - 1} |" in l_ and f"[{min(r['cuts'][0]['t'] for r in raw['L']):.5f}, 1.50000] | {le0 - 1} |" in l_
    rec0 = full[full.index('| 0.0001 | 120 | '):].split(' | ')[2]
    assert f"| 0.0001 | 120 | {int(rec0) - 1} |" in m_
