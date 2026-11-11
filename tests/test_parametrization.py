"""Tests of experiments/exp_parametrization.py: the target sample, the reference path, the reformulation's checks (part A), the cell
records, and that the committed report is the report of the committed raw."""
import os
import sys
import warnings

import numpy as np
import pytest

ROOT = os.path.join(os.path.dirname(__file__), '..')
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'experiments'))
import exp_cohesion_leak as E
import exp_parametrization as P

warnings.filterwarnings('ignore')
RAW = os.path.join(ROOT, 'results', 'chaos', 'parametrization_raw.npz')
MD = os.path.join(ROOT, 'results', 'chaos', 'parametrization.md')


def test_target_sample_is_deterministic_and_stratified():
    tg = P.sample_targets()
    assert len(tg) == 200
    fams = {k: sum(t['fam'] == k for t in tg) for k in ('A', 'B', 'Da', 'Db', 'Dc')}
    assert fams == dict(A=70, B=30, Da=20, Db=50, Dc=30)
    assert sum(t['out'] == 'nearest' for t in tg) == 118 and sum(t['out'] == 'other' for t in tg) == 82
    assert np.allclose(tg[0]['cs'], [0.59799952, -0.03874132, -0.54359771], rtol=0, atol=1e-8)
    assert np.allclose(tg[0]['As'], [0.99698898, -0.63429886, -0.65616795], rtol=0, atol=1e-8)
    assert tg[0]['out'] == 'nearest' and tg[0]['order'] == [16, 5, 10]
    again = P.sample_targets()
    assert all(np.array_equal(a['cs'], b['cs']) and np.array_equal(a['As'], b['As']) for a, b in zip(tg, again))


def test_the_sample_has_two_or_three_features_inside_the_window():
    """The designed families D(a), D(b) are pairs (K = 2); A, B and D(c) are triples."""
    for t in P.sample_targets():
        K = 2 if t['fam'] in ('Da', 'Db') else 3
        assert t['cs'].size == K and t['As'].size == K and np.all(np.abs(t['cs']) <= 0.8) and np.all(np.abs(t['As']) > 0), t['fam']
    wide = [t for t in P.sample_targets() if np.abs(t['As']).max() / np.abs(t['As']).min() > 100]
    assert len(wide) >= 10                                  # the targets for which the formed threshold must follow the smallest feature


def test_reference_path_of_the_design():
    L, K, c0 = P.setup()
    f = K.synth(np.array(P.CSTAR)[:, None], np.array(P.ASTAR))
    ref = P.reference(L, K, c0, f)
    assert ref['order'] == [5, 16, 11] and ref['set'] == [5, 11, 16] and ref['exits'] == 0
    assert np.allclose(ref['tau'], [10.7217, 12.8471, 16.3858], rtol=0, atol=2e-4)       # Paper 17: 10.72, 12.85, 16.39
    assert ref['rms'] < 1e-12


def test_reference_path_reproduces_the_stored_order_on_a_few_targets():
    L, K, c0 = P.setup()
    for t in P.sample_targets()[:6]:
        ref = P.reference(L, K, c0, K.synth(t['cs'][:, None], t['As']))
        assert ref['order'] == t['order']                           # the stored entry sequence of `exp_transit_certificate`


def test_part_A_identities():
    A = P.part_A()
    assert max(A['identity']) < 1e-14                                # the law is the Euclidean flow, to rounding
    assert [r['eps'] for r in A['flows']] == [1e-3, 1e-4, 1e-6]
    assert all(r['dw'] < 5e-6 and r['dc'] < 5e-6 for r in A['flows'])   # raw (α, ν, c) flow against the reduced law on its curves
    assert {r['readout'] for r in A['conserved']} == {'P0', 'D2', 'D3'}
    assert all(r['drift'] < 1e-6 for r in A['conserved'])


def test_cell_record_is_flat_and_complete():
    L, K, c0 = P.setup()
    cs, As = np.array(P.CSTAR), np.array(P.ASTAR)
    f = K.synth(cs[:, None], As)
    ref = P.reference(L, K, c0, f)
    a = P.cell(P.make_readout('D3'), K, f, c0, cs, As, 1e-6, ref, 20.0)
    for k in ('n_formed', 'n_formed_sweep', 'gap', 'dormant_disp', 'dormant_wmax', 'fit_rms', 'exact', 'set_equal', 'order_equal', 'tau_rel',
              'secs', 'status', 'nfev', 'eps', 'readout'):
        assert k in a
    assert a['n_formed'] == 3 and a['readout'] == 'D3' and a['status'] == 'ok' and len(a['n_formed_sweep']) == 4
    # the centres' rate is a parameter: with fixed centres the network is a dictionary and one unit per off-grid feature cannot fit it
    b = P.cell(P.make_readout('D3'), K, f, c0, cs, As, 1e-6, ref, 20.0, lr_c=0.0)
    assert not b['exact'] and b['fit_rms'] > 1e-3


def test_readout_factory():
    assert [P.make_readout(n).name for n in ('P0', 'P0u', 'P4', 'D3', 'D2', 'D1')] == ['P0', 'P0u', 'P0', 'D3', 'D2', 'D1']
    assert np.allclose(P.make_readout('P4').lr, [12.0, 12.0], rtol=1e-14, atol=0)              # equal learning rates
    assert np.allclose(P.make_readout('P0').lr, [12.0, 0.3], rtol=1e-14, atol=0)
    with pytest.raises(ValueError):
        P.make_readout('D9')


def test_span_formats_the_range_of_the_finite_values():
    assert P.span([0.5, float('nan'), 0.25, 0.4]) == '0.25 – 0.5' and P.span([float('nan')]) == 'nan' and P.span([]) == 'nan'
    assert P.span([3.14159, 2.71828], '.2f') == '2.72 – 3.14' and P.span(np.array([1e-5, 2e-5])) == '1e-05 – 2e-05'


def test_report_formatting_helpers():
    assert P.g(float('nan')) == 'nan' and P.g(None) == 'nan' and P.g(0.000123456) == '0.000123' and P.g(12.34, '.1f') == '12.3'
    t = P.table(['a', 'b'], [[1, 'x'], [2, 'y']])
    assert t[0] == '| a | b |' and t[1] == '|---|---|' and t[2] == '| 1 | x |' and t[-1] == ''
    assert P.pct(3, 5) == '3 of 5'


@pytest.mark.skipif(not os.path.exists(RAW), reason='the raw is committed with the report')
def test_the_committed_report_is_the_report_of_the_committed_raw():
    raw = E.load(RAW)
    assert set('ABDEFGH') <= set(raw) <= set(P.PARTS)                   # C (the 200 targets), I (centres' rate) and J (the plane): long runs, merged when finished
    text = P.report(*[raw.get(k) for k in P.PARTS])
    assert text == open(MD, encoding='utf-8').read()
    assert P.report(*[raw.get(k) for k in P.PARTS]) == text         # deterministic


@pytest.mark.skipif(not os.path.exists(RAW), reason='the raw is committed with the report')
def test_the_raw_pins_the_claims():
    """The conclusions of WP1, pinned against the committed raw (a change of the engine or the analysis that moves them fails here)."""
    raw = E.load(RAW)
    rows = {(r['name'], r['eps']): r for r in raw['B']['rows']}
    for n in ('P0', 'P0u', 'P4', 'D3'):
        for e in (1e-6, 1e-8, 1e-10):
            r = rows[(n, e)]
            assert r['n_formed'] == 3 and r['set_equal'] and r['exact'] and r['order_equal'], (n, e)
    for e in P.EPS_B:
        r = rows[('D2', e)]
        assert r['n_formed'] >= 6 and not r['exact'], e           # depth 2: no clean selection at any practical scale
        assert rows[('D1', e)]['n_formed'] >= 8 and not rows[('D1', e)]['set_equal']
    for r in raw['F']:                                              # fixed-step descent: no unstable iteration is reported as evidence
        assert r['h'] * r['lam_max'] < 2.0, (r['i'], r['name'], r['h'])
    f6 = [r for r in raw['F'] if r['eps'] == 1e-6 and r['name'] in ('P0', 'D3')]
    assert len(f6) >= 5 and all(r['exact'] and r['set_equal'] and r['n_formed'] == 3 for r in f6)      # stable descent follows the flow
    if 'C' in raw:                                                  # the 200 benchmark targets (a long run: merged when finished)
        C = raw['C']
        assert len(C) == 200
        for k, (nm, ep) in enumerate(P.CELLS_C):
            done = [r['cells'][k] for r in C if r['cells'][k]['status'] == 'ok']
            if nm == 'D2' and ep >= 1e-12:
                assert sum(c['set_equal'] for c in done) <= 0.2 * len(done), (nm, ep)    # depth 2 at a practical scale: rarely the path's set
                assert not any(c['exact'] for c in done), (nm, ep)
            if nm in ('P0', 'D3'):                                  # depth 3: K units on the features on nearly every target, the path's set on most
                assert len(done) == 200 and sum(c['exact'] for c in done) >= (0.9 if ep == 1e-6 else 0.95) * len(done), (nm, ep)
                assert sum(c['set_equal'] for c in done) >= (0.8 if ep == 1e-6 else 0.9) * len(done), (nm, ep)
    for part in ('I', 'J'):                                         # the centres' rate: fixed centres never fit exactly; depth 2 never does
        if part in raw and part == 'I':
            rows = raw['I'][0]['rows']
            for j, r0 in enumerate(rows):
                ok = [t['rows'][j] for t in raw['I'] if t['rows'][j]['status'] == 'ok']
                if r0['mult'] == 0.0:
                    assert not any(r['exact'] for r in ok), r0['name']
                if r0['name'] == 'D3' and r0['mult'] > 0:
                    assert all(r['exact'] for r in ok), r0['mult']
                if r0['name'] == 'D2':
                    assert not any(r['exact'] for r in ok), r0['mult']
    rows = {(r['name'], r['eps']): r for r in raw['B']['rows']}                      # depth 4 and 5 (the referee's check of "depth >= 3")
    for nm in ('D4', 'D5'):
        for e in (1e-4, 1e-6, 1e-8):
            assert rows[(nm, e)]['n_formed'] == 3 and rows[(nm, e)]['set_equal'] and rows[(nm, e)]['exact'], (nm, e)
    assert rows[('D5', 1e-10)]['status'].startswith('failed')                       # a recorded solver failure at t ~ 1e13, not a success
    D = raw['D']                                                                     # the depth-2 rule, with no fitted parameter
    for j, e in enumerate(P.EPS_D):
        thr = np.log(1 / P.FR.THETA) / np.log(1 / e)
        assert sum((t['margin'] > thr) == t['rows'][j]['set_equal'] for t in D) >= 19, e
    assert [sum(t['rows'][j]['exact_topk'] for t in D) for j in range(len(P.EPS_D))] == [0, 1, 10, 13, 20, 21]     # the referee's counts
    if 'C' in raw:                                                                   # the margin law: no mismatch where the runner-up margin is >= 0.02
        C = raw['C']
        hi = [r for r in C if r['ref']['margin'] >= 0.02]
        assert len(hi) == 98
        for k, (nm, ep) in enumerate(P.CELLS_C):
            if nm in ('P0', 'D3'):
                assert all(r['cells'][k]['set_equal'] and r['cells'][k]['exact'] and r['cells'][k]['order_equal'] for r in hi), (nm, ep)
        assert sum(c.get('radau', 0) > 0 for r in C for c in r['cells']) == 104      # cells that used the Radau fallback are recorded
        assert all(r['ref']['relaxed'] and r['ref']['solver_failed'] == 0 and not r['ref']['out_of_budget'] for r in C)
    d2 = raw['D']
    prod = [np.nanmedian([t['rows'][j]['tau_rel'] for t in d2]) * np.log(1 / P.EPS_D[j]) for j in range(len(P.EPS_D))]
    assert max(prod) / min(prod) < 1.05                           # τ error × ln(1/ε) is constant: the logarithmic law


def test_lsoda_failure_on_a_wide_ratio_pair_is_redone_by_radau():
    """Target 108 of the sample (a pair with amplitudes 1 and −1000): LSODA reports "Unexpected istate" on a relaxed state, which crashed the
    first long run; the flow is now redone by Radau and finishes with the path's two units."""
    L, K, c0 = P.setup()
    t = P.sample_targets()[108]
    assert t['fam'] == 'Da' and np.abs(t['As']).max() / np.abs(t['As']).min() > 500
    f = K.synth(t['cs'][:, None], t['As'])
    ref = P.reference(L, K, c0, f)
    a = P.cell(P.make_readout('P0'), K, f, c0, t['cs'], t['As'], 1e-6, ref, 1.05 * max(ref['tau']) + 1.0)
    assert a['status'] == 'ok' and a['n_formed'] == 2 and a['set_equal'] and a['exact']


def test_a_failing_run_is_a_recorded_result_not_a_crash(monkeypatch):
    L, K, c0 = P.setup()
    cs, As = np.array(P.CSTAR), np.array(P.ASTAR)
    f = K.synth(cs[:, None], As)
    ref = P.reference(L, K, c0, f)

    def boom(*a, **k):
        raise RuntimeError('solver exploded')
    monkeypatch.setattr(P.FR, 'train_flow', boom)
    a = P.cell(P.make_readout('D3'), K, f, c0, cs, As, 1e-6, ref, 20.0)
    assert a['status'].startswith('error') and not a['exact'] and not a['set_equal'] and a['n_formed'] == 0
    assert np.isnan(a['fit_rms']) and a['readout'] == 'D3'


def test_runner_up_margins_of_the_design_and_of_a_target_the_referee_recomputed():
    """`margins`: the smallest 1 − |u_j| over the kernels not yet active at an entry. The values were recomputed independently by the referee
    (target 59: 0.044); the margin is the distance from the limit that every finite-ε statement has to be read against."""
    L, K, c0 = P.setup()
    r = P.reference(L, K, c0, K.synth(np.array(P.CSTAR)[:, None], np.array(P.ASTAR)))
    assert np.isclose(r['margin'], 0.0664, atol=1e-3, rtol=0) and r['margin_never'] >= r['margin']
    assert r['relaxed'] and r['solver_failed'] == 0 and not r['out_of_budget']
    t = P.sample_targets()[59]
    r = P.reference(L, K, c0, K.synth(t['cs'][:, None], t['As']))
    assert np.isclose(r['margin'], 0.0437, atol=1e-3, rtol=0) and 0 < r['margin'] < 1
    # unequal masses: the duals are reconstructed for equal masses only, so no margin is reported
    rho = np.ones(len(c0))
    rho[3] = 1.05
    r2 = P.reference(L, K, c0, K.synth(t['cs'][:, None], t['As']), rho=rho)
    assert np.isnan(r2['margin']) and np.isnan(r2['margin_never'])


@pytest.mark.skipif(not os.path.exists(RAW), reason='the raw is committed with the report')
def test_the_report_states_the_status_of_every_cell():
    import copy
    raw = E.load(RAW)
    B = copy.deepcopy(raw['B'])
    n = len(B['rows'])
    bad0 = sum(r['status'] != 'ok' for r in B['rows'])           # B holds one recorded solver failure (D5 at 1e-10)
    text = P.report(B=B)
    assert f'B: {n} / {bad0}' in text and 'never counted as a success' in text
    first_ok = next(r for r in B['rows'] if r['status'] == 'ok')
    first_ok['status'] = 'timeout'                               # a further non-ok cell is reported, not hidden
    assert f'B: {n} / {bad0 + 1}' in P.report(B=B)



# -- parts K (branch pairs) and L (the dictionary's spacing) -------------------------------------------------------------------------
def test_cell_atoms_controls_go_through_the_same_analysis_as_cell():
    """`cell_atoms` with a one-unit-per-atom readout is `cell` (the same record for the keys they share) plus the atoms' counts."""
    L, K, c0 = P.setup()
    cs, As = np.array(P.CSTAR), np.array(P.ASTAR)
    f = K.synth(cs[:, None], As)
    ref = P.reference(L, K, c0, f)
    a = P.cell(P.make_readout('D3'), K, f, c0, cs, As, 1e-6, ref, 20.0)
    b = P.cell_atoms(P.make_readout('D3'), K, f, c0, cs, As, 1e-6, ref, 20.0)
    for k in ('n_formed', 'n_formed_sweep', 'gap', 'gap_k', 'dormant_disp', 'dormant_wmax', 'fit_rms', 'exact', 'exact_topk', 'set_equal', 'order_equal',
              'tau_rel', 'pos_err', 'amp_err', 'form_order', 'status', 'nfev', 'eps', 'readout', 'radau', 'aborted'):
        assert a[k] == b[k], k
    assert b['n_atoms'] == b['n_formed'] == 3 and b['atoms'] == [5, 11, 16] and b['one_each'] and b['tau_ratio'] < 1e-3


def test_cell_atoms_branch_pairs_on_the_design():
    L, K, c0 = P.setup()
    cs, As = np.array(P.CSTAR), np.array(P.ASTAR)
    f = K.synth(cs[:, None], As)
    ref = P.reference(L, K, c0, f)
    ro = P.FR.PaperPairs(P.law(), 0.1)
    a = P.cell_atoms(ro, K, f, c0, cs, As, 1e-6, ref, 20.0)
    assert a['status'] == 'ok' and a['n_formed'] == 3 and a['n_atoms'] == 3 and a['atoms'] == [5, 11, 16] and a['one_each']
    assert a['set_equal'] and a['order_equal'] and a['exact'] and a['readout'] == 'P0p' and len(a['n_formed_sweep']) == 4
    assert 1.8 < a['gap_k'] < 2.0 and a['gap_k'] == a['gap']
    assert abs(a['tau_rel'] - (1 - ro.form_ratio(1e-6, P.FR.formed_threshold(As)))) < 5e-3 and a['tau_ratio'] < 5e-3      # the derived factor: here every atom's (separated features)
    assert a['form_order'] == [5, 16, 11]                                  # the atoms in the path's order
    # a₀ = 0.3 at 10⁻⁶ forms six atoms: not the path's set, not exact (the pair start is not depth 3)
    b = P.cell_atoms(P.FR.PaperPairs(P.law(), 0.3), K, f, c0, cs, As, 1e-6, ref, 20.0)
    assert b['n_formed'] == 6 and b['n_atoms'] == 6 and not b['set_equal'] and not b['exact'] and b['gap_k'] < 1.0


def test_cell_atoms_exactness_needs_one_unit_per_atom(monkeypatch):
    """A fourth unit on a formed atom's partner (both units of the pair above the threshold) is not 'one unit each': three atoms, four units."""
    L, K, c0 = P.setup()
    cs, As = np.array(P.CSTAR), np.array(P.ASTAR)
    f = K.synth(cs[:, None], As)
    ref = P.reference(L, K, c0, f)
    real = P.FR.train_flow

    def with_partner(ro, K_, f_, c0u, eps, tau_end, lr_c, **kw):
        res = real(ro, K_, f_, c0u, eps, tau_end, lr_c, **kw)
        res['W'] = res['W'].copy()
        res['W'][-1, 26] = -0.3                                              # the − unit of atom 5, started with the + unit (5) that has formed
        res['C'] = res['C'].copy()
        res['C'][-1, 26, 0] = 0.07                                           # (on another feature: it would not be the nearest unit of the −0.52 feature)
        return res
    monkeypatch.setattr(P.FR, 'train_flow', with_partner)
    a = P.cell_atoms(P.FR.PaperPairs(P.law(), 0.1), K, f, c0, cs, As, 1e-6, ref, 20.0)
    assert a['n_formed'] == 4 and a['n_atoms'] == 3 and not a['one_each'] and not a['exact'] and a['set_equal']


def test_cell_atoms_failing_run_is_a_recorded_result(monkeypatch):
    L, K, c0 = P.setup()
    cs, As = np.array(P.CSTAR), np.array(P.ASTAR)
    f = K.synth(cs[:, None], As)
    ref = P.reference(L, K, c0, f)

    def boom(*a, **k):
        raise RuntimeError('solver exploded')
    monkeypatch.setattr(P.FR, 'train_flow', boom)
    a = P.cell_atoms(P.FR.PaperPairs(P.law(), 0.3), K, f, c0, cs, As, 1e-6, ref, 20.0)
    assert a['status'].startswith('error') and not a['exact'] and not a['set_equal'] and a['n_formed'] == 0 and a['n_atoms'] == 0
    assert np.isnan(a['fit_rms']) and np.isnan(a['tau_ratio']) and a['readout'] == 'P0p' and a['eps'] == 1e-6


def test_cell_atoms_exactness_needs_k_atoms_not_only_k_units(monkeypatch):
    """The units' test (`analyse`: K units formed, one per feature, within the tolerances) can hold with only K − 1 ATOMS: the two units of one atom matched
    to two features (here atom 5's pair: unit 5 on feature 0, its partner unit 26 on feature 1, and atom 16's unit 37 on feature 2). `cell_atoms` must
    then not call the run exact (the clause n_atoms == K; a mutant without it passed every test, since the other test had four units)."""
    L, K, c0 = P.setup()
    cs, As = np.array(P.CSTAR), np.array(P.ASTAR)
    f = K.synth(cs[:, None], As)
    ref = P.reference(L, K, c0, f)

    def synthetic(ro, K_, f_, c0u, eps, tau_end, lr_c, **kw):
        T = 5
        W = np.zeros((T, 42))
        C = np.tile(np.asarray(c0u, float)[None], (T, 1, 1)).copy()
        for u, k in ((5, 0), (26, 1), (37, 2)):
            W[2:, u] = As[k]
            C[-1, u, 0] = cs[k]
        return dict(taus=np.linspace(0.0, tau_end, T), W=W, C=C, status='ok', nfev=0, secs=0.0, radau_pieces=0, aborted_nfev=0)
    monkeypatch.setattr(P.FR, 'train_flow', synthetic)
    a = P.cell_atoms(P.FR.PaperPairs(P.law(), 0.1), K, f, c0, cs, As, 1e-6, ref, 20.0)
    res = synthetic(None, K, f, np.r_[c0, c0][:, None], 1e-6, 20.0, 12.0)
    units = P.FR.analyse(res, K, f, np.r_[c0, c0][:, None], cs, As)
    assert units['exact'] and units['n_formed'] == 3                                    # as units: three formed, one on each feature
    assert a['n_formed'] == 3 and a['n_atoms'] == 2 and a['atoms'] == [5, 16] and not a['one_each'] and not a['exact']      # as atoms: two


def test_the_k_job_predicts_with_the_targets_own_formed_threshold(monkeypatch):
    """`job_k`: the pairs' predicted earliness is 1 − form_ratio at the target's own formed threshold min(0.01 max|A*|, 0.1 min|A*|), not at the design's
    0.009: a target with A = (0.5, −0.25, 0.3) has 0.005 (a mutant with 0.009 fixed passed every test: invisible on the design and within the 0.02
    tolerance of the raw test elsewhere). The cells are stubbed; the grid of rows is the job's."""
    tgt = dict(cs=np.array([-0.4, 0.1, 0.6]), As=np.array([0.5, -0.25, 0.3]))
    thr = P.FR.formed_threshold(tgt['As'])
    assert np.isclose(thr, 0.005, rtol=1e-14, atol=0)
    seen = []
    monkeypatch.setattr(P, 'reference', lambda L, K, c0, f: dict(order=[1], tau=[1.0], set=[1], margin=0.1, margin_never=0.2))

    def fake(ro, K, f, c0, cs, As, eps, ref, tau_end, lr_c=None):
        seen.append((ro.name, float(eps), float(tau_end)))
        return dict(status='ok', eps=float(eps), readout=ro.name)
    monkeypatch.setattr(P, 'cell_atoms', fake)
    out = P.job_k((3, tgt))
    assert out['i'] == 3 and out['kind'] == 'std' and out['margin'] == 0.1 and np.isclose(out['sep'], 0.5 / 0.12, rtol=1e-12, atol=0)
    pairs = [r for r in out['rows'] if r['start'].startswith('pairs')]
    assert [(r['a0'], r['eps']) for r in pairs] == [(a0, e) for a0 in P.A0_K for e in P.EPS_K + P.EPS_K_FAR]
    for r in pairs:
        ro = P.FR.PaperPairs(P.law(), r['a0'])
        assert r['tau_pred'] == 1.0 - ro.form_ratio(r['eps'], thr)
        assert abs(r['tau_pred'] - (1.0 - ro.form_ratio(r['eps'], 0.009))) > 1e-4             # the design's threshold would give another prediction (by 6·10⁻⁴ at 10⁻⁴⁸)
    ctl = [(r['start'], r['eps']) for r in out['rows'] if not r['start'].startswith('pairs')]
    assert ctl == [('P0, α(0) = 0', e) for e in P.EPS_K] + [('D3', e) for e in P.EPS_K] + [('D2', e) for e in P.EPS_K + P.EPS_K_FAR]
    assert all(np.isnan(r['tau_pred']) for r in out['rows'] if not r['start'].startswith('pairs'))
    assert len(seen) == len(out['rows']) and all(np.isclose(t, 2.05, rtol=1e-14, atol=0) for _, _, t in seen)       # τ_end = 1.05 max τ + 1 of the stubbed path


@pytest.mark.skipif(not os.path.exists(RAW), reason='the raw is committed with the report')
def test_part_K_job_of_a_sample_target_recomputes_the_raw(monkeypatch):
    """A target whose formed threshold is not the design's (sample target 0: 0.00997), on a reduced grid (a₀ = 0.1 at 10⁻⁶ and D2), equals the raw's cells and
    carries the prediction at its own threshold."""
    raw = E.load(RAW)
    if 'K' not in raw:
        pytest.skip('part K is merged when finished')
    tgt = P.sample_targets()[0]
    assert not np.isclose(P.FR.formed_threshold(tgt['As']), 0.009, rtol=1e-3, atol=0)
    for name, val in (('A0_K', (0.1,)), ('EPS_K', (1e-6,)), ('EPS_K_FAR', ()), ('NAMES_K', ('D2',))):
        monkeypatch.setattr(P, name, val)
    got = P.job_k((1, tgt))
    want = next(t for t in raw['K'] if t.get('kind', 'std') == 'std' and t['i'] == 1)
    assert [(r['start'], r['eps']) for r in got['rows']] == [('pairs a₀ = 0.1', 1e-6), ('D2', 1e-6)]
    for r in got['rows']:
        _close(r, next(x for x in want['rows'] if (x['start'], x['eps']) == (r['start'], r['eps'])))
    assert np.isclose(got['rows'][0]['tau_pred'], 1 - P.FR.PaperPairs(P.law(), 0.1).form_ratio(1e-6, P.FR.formed_threshold(tgt['As'])), rtol=1e-14, atol=0)
    assert np.isclose(got['margin'], want['margin'], rtol=1e-12, atol=0) and np.isclose(got['sep'], want['sep'], rtol=1e-12, atol=0)


def test_the_d_job_stores_both_margins_of_the_path(monkeypatch):
    """`job_d` (PZ2-6): the rule of part D reads μ never as well as μ; both come from the target's own reference path (cells stubbed)."""
    tgt = dict(cs=np.array(P.CSTAR), As=np.array(P.ASTAR), fam='design', out='nearest', order=[], e=[])
    monkeypatch.setattr(P, 'reference', lambda L, K, c0, f: dict(order=[1], tau=[2.0], set=[1], margin=0.11, margin_never=0.22))
    seen = []

    def fake(ro, K, f, c0, cs, As, eps, ref, tau_end, **kw):
        seen.append((ro.name, float(eps), float(tau_end)))
        return dict(set_equal=False, eps=float(eps))
    monkeypatch.setattr(P, 'cell', fake)
    out = P.job_d((7, tgt))
    assert out['i'] == 7 and out['margin'] == 0.11 and out['margin_never'] == 0.22
    assert [r['eps'] for r in out['rows']] == list(P.EPS_D) and [r['eps'] for r in out['rows_long']] == list(P.EPS_D[:2])
    assert np.allclose([t for _, _, t in seen[:6]], 3.1, rtol=1e-14, atol=0) and np.allclose([t for _, _, t in seen[6:]], 9.3, rtol=1e-14, atol=0)      # τ_end, 3τ_end
    assert {n for n, _, _ in seen} == {'D2'} and len(seen) == 8


def test_the_coupled_targets_are_the_sixteen_close_pairs_of_family_Db():
    """K = 2 targets of C's sample whose features are at most 2σ apart: family D(b) 120 … 135, amplitudes (1, 1) or (1, −1), 1.5σ and 2σ apart."""
    cp = P.coupled_targets()
    assert [i for i, t in cp] == list(range(120, 136))
    same = [bool(t['As'][0] * t['As'][1] > 0) for i, t in cp]
    assert sum(same) == 8 and len(same) - sum(same) == 8
    assert {i for (i, t), sm in zip(cp, same) if sm} == {120, 121, 122, 123, 124, 125, 130, 131}
    seps = np.array([abs(t['cs'][1] - t['cs'][0]) / P.SIG for i, t in cp])
    assert np.allclose(seps[:10], 1.5, rtol=0, atol=1e-9) and np.allclose(seps[10:], 2.0, rtol=0, atol=1e-9)
    assert all(t['fam'] == 'Db' and np.allclose(np.abs(t['As']), 1.0) for i, t in cp)


def test_the_first_entrant_is_early_by_the_derived_factor_the_second_atom_of_a_same_sign_pair_by_more():
    """The PZ2 referee's finding, recomputed (a₀ = 0.1, ε = 10⁻²⁴): target 121 (features 1.5σ apart, A = (1, 1)): the first entrant is early by 1 − form_ratio
    (0.1122), the second atom by 0.48, four times that, the entry-time ratios off by 0.42; target 126 (A = (1, −1)): both atoms by the prediction. The
    derived factor holds for the first entrant and for later atoms only while their duals are straight (formed is not entered)."""
    L, K, c0 = P.setup()
    tg = P.sample_targets()
    for i, same in ((121, True), (126, False)):
        t = tg[i]
        assert (t['As'][0] * t['As'][1] > 0) == same
        f = K.synth(t['cs'][:, None], t['As'])
        ref = P.reference(L, K, c0, f)
        ro = P.FR.PaperPairs(P.law(), 0.1)
        a = P.cell_atoms(ro, K, f, c0, t['cs'], t['As'], 1e-24, ref, 1.05 * max(ref['tau']) + 1.0)
        pred = 1 - ro.form_ratio(1e-24, P.FR.formed_threshold(t['As']))
        e = a['earliness']
        assert a['status'] == 'ok' and a['set_equal'] and len(e) == 2 and abs(e[0] - pred) < 5e-3, i             # the first entrant: the derived factor
        if same:
            assert e[1] > 3.5 * pred and a['tau_ratio'] > 0.3, i
        else:
            assert abs(e[1] - pred) < 0.01 and a['tau_ratio'] < 0.01, i
        assert np.isclose(max(e), a['tau_rel'], rtol=1e-12, atol=0) and len(a['form_times']) == 2


def test_the_rule_reads_nan_as_infinity_and_a_three_kernel_dictionary_has_no_kernel_that_never_enters():
    """μ never is NaN when every kernel enters (nothing is left to violate the rule): the rule reads +∞ there. With the raw comparison `NaN > threshold`
    is False and the rule would predict a failure that depth 2 does not have (3 kernels, 3 atoms)."""
    nan = float('nan')
    assert P.rule_margin(nan) == float('inf') and P.rule_margin(None) == float('inf') and P.rule_margin(0.25) == 0.25
    assert np.isclose(P.rule_threshold(1e-6), np.log(100) / np.log(1e6), rtol=1e-14, atol=0) and P.rule_threshold(1e-12) < P.rule_threshold(1e-6)
    assert np.isclose(P.rule_threshold(1e-6), 1 / 3, rtol=1e-12, atol=0) and np.isclose(P.rule_threshold(1e-12), 1 / 6, rtol=1e-12, atol=0)
    assert P.rule_passes(nan, 1e-6) and P.rule_passes(float('inf'), 1e-96) and not P.rule_passes(0.2, 1e-6) and P.rule_passes(0.2, 1e-12)
    assert not P.rule_passes(P.rule_threshold(1e-6), 1e-6)                                  # strict: the margin must exceed the threshold
    assert P.rule_passes(np.nextafter(P.rule_threshold(1e-6), 1.0), 1e-6)
    L, K, c0 = P.setup_nd(3)
    cs, As = np.array(P.CSTAR), np.array(P.ASTAR)
    f = K.synth(cs[:, None], As)
    ref = P.reference(L, K, c0, f)
    assert np.isnan(ref['margin_never']) and 0 < ref['margin'] < 1                          # every one of the 3 kernels enters
    a = P.cell(P.make_readout('D2'), K, f, c0, cs, As, 1e-6, ref, 1.05 * max(ref['tau']) + 1.0)
    assert a['set_equal'] and not (ref['margin_never'] > P.rule_threshold(1e-6))            # the unconverted comparison would predict a failure
    assert P.rule_passes(ref['margin_never'], 1e-6) == a['set_equal']


@pytest.mark.skipif(not os.path.exists(RAW), reason='the raw is committed with the report')
def test_the_report_reads_a_nan_margin_never_as_infinity_in_parts_D_and_L():
    import copy
    raw = E.load(RAW)
    D = copy.deepcopy(raw['D'])
    for t in D:
        t['margin_never'] = float('nan')
    text = P.report(D=D)
    for j, e in enumerate(P.EPS_D):                                                         # +∞ > threshold: the rule predicts the path's set on every target
        n_set = sum(t['rows'][j]['set_equal'] for t in D)
        assert f"| {e:.0e} | {len(D)} of {len(D)} | {P.g(P.rule_threshold(e))} | {n_set} of {len(D)} |" in text
    Ls = copy.deepcopy(raw['L'])
    t21 = next(t for t in Ls if t['variant'] == 'off' and t['nd'] == 21)
    t21['ref']['margin_never'] = float('nan')
    text = P.report(Ls=Ls)
    rows = [l for l in text.splitlines() if l.startswith('| 21 | 0.1 | 0.0664 | nan | D2 |')]
    assert len(rows) == 2 and '| 0.333: yes |' in rows[0] and '| 0.167: yes |' in rows[1]


def test_the_dictionaries_of_part_L():
    """21 … 7 kernels at spacing 2/(n − 1); the on-grid arm — a control that removes the sliding — is run at EVERY spacing (the first version ran it
    only where neighbouring kernels overlap by at most 0.5, which are exactly the spacings where depth 2 passes); its features are the kernel centres
    nearest to the design's; the scan runs D2 over ε = 10^(−6) … 10^(−90) on six dictionaries."""
    assert P.ND_L == (21, 16, 13, 11, 9, 7) and P.ND_ON_L == P.ND_L and P.EPS_L == (1e-6, 1e-12) and P.NAMES_L == ('P0', 'D3', 'D2')
    assert P.ND_SCAN_L == (11, 13, 16, 21, 25, 31) and P.K_SCAN_L == tuple(range(6, 91, 2)) and len(P.K_SCAN_L) == 43
    L, K, c0 = P.setup()
    L21, K21, c21 = P.setup_nd(21)
    assert np.array_equal(c21, c0) and K21.sig == K.sig == 0.12 and K21.xs.shape == (61, 1)
    low_overlap = []
    for n in P.ND_L:
        _, Kn, cn = P.setup_nd(n)
        sp = 2.0 / (n - 1)
        assert cn.size == n and np.isclose(cn[1] - cn[0], sp, rtol=1e-12, atol=0) and cn[0] == -1.0 and cn[-1] == 1.0
        overlap = np.exp(-sp ** 2 / (4 * P.SIG ** 2))                                   # ⟨φ_iφ_{i+1}⟩/⟨φ²⟩
        if overlap <= 0.5 + 1e-12:
            low_overlap.append(n)
    assert low_overlap == [11, 9, 7]                                                    # the first version's on-grid arm; 21, 16, 13 overlap by more
    for n, want in ((11, [-0.6, 0.0, 0.6]), (9, [-0.5, 0.0, 0.5]), (7, [-2 / 3, 0.0, 2 / 3]), (21, [-0.5, 0.1, 0.6]), (16, [-1 + 8 / 15, -1 + 16 / 15, -1 + 24 / 15]),
                    (13, [-0.5, 0.0, 2 / 3])):
        assert np.allclose(P.on_grid(np.linspace(-1, 1, n)), want, rtol=0, atol=1e-14), n
    for n in P.ND_ON_L:
        c = np.linspace(-1, 1, n)
        g_ = P.on_grid(c)
        assert all(np.isclose(c, x, rtol=0, atol=1e-14).any() for x in g_) and len(set(np.round(g_, 9))) == 3


def test_part_K_has_the_design_and_ten_targets_and_the_controls():
    assert P.A0_K == (0.1, 0.3) and P.EPS_K == (1e-6, 1e-12) and P.N_K == 10 and P.NAMES_K == ('P0', 'D3', 'D2')
    # the PZ2 round: the pairs and D2 further down; the coupled targets; the scan of a₀
    assert P.EPS_K_FAR == (1e-24, 1e-48, 1e-96) and P.A0_KC == 0.1 and P.EPS_KC == (1e-6, 1e-12, 1e-24, 1e-48) and np.isclose(P.SEP_KC, 0.24, rtol=1e-14, atol=0)
    assert P.A0_SCAN_K == (0.003, 0.01, 0.03, 0.1, 0.2, 0.3, 0.5, 1.0, 3.0) and P.EPS_SCAN_K == (1e-6, 1e-12)


@pytest.mark.skipif(not os.path.exists(RAW), reason='the raw is committed with the report')
def test_the_raw_pins_part_K_the_branch_pairs():
    raw = E.load(RAW)
    if 'K' not in raw:
        pytest.skip('part K is merged when finished')
    Ks = [t for t in raw['K'] if t.get('kind', 'std') == 'std']                      # the design and the first 10 targets (the coupled targets and the a₀ scan follow)
    assert len(Ks) == 1 + P.N_K and [t['i'] for t in Ks] == list(range(11))
    assert [t.get('kind', 'std') for t in raw['K']] == ['std'] * 11 + ['coupled'] * 16 + ['a0scan']
    keys = [(r['start'], r['eps']) for r in Ks[0]['rows']]
    assert keys == ([(f'pairs a₀ = {a0:g}', e) for a0 in P.A0_K for e in P.EPS_K + P.EPS_K_FAR]
                    + [(n + (', α(0) = 0' if n == 'P0' else ''), e) for n in P.NAMES_K for e in P.EPS_K + (P.EPS_K_FAR if n == 'D2' else ())])
    assert all(r['status'] == 'ok' for t in raw['K'] for r in t['rows'])
    design = {(r['start'], r['eps']): r for r in Ks[0]['rows']}
    # the referee's validation cells, on the design
    r = design[('pairs a₀ = 0.3', 1e-6)]
    assert r['n_formed'] == 6 and r['n_atoms'] == 6 and not r['set_equal'] and not r['exact'] and np.isclose(r['gap_k'], 0.71, atol=0.01, rtol=0)
    r = design[('pairs a₀ = 0.3', 1e-12)]
    assert r['n_formed'] == 3 and r['set_equal'] and r['exact'] and np.isclose(r['gap_k'], 2.16, atol=0.01, rtol=0)
    r = design[('pairs a₀ = 0.1', 1e-6)]
    assert r['n_formed'] == 3 and r['set_equal'] and r['exact'] and np.isclose(r['gap_k'], 1.89, atol=0.01, rtol=0)
    assert np.isclose(design[('D3', 1e-6)]['gap_k'], 5.27, atol=0.01, rtol=0) and design[('D2', 1e-6)]['n_formed'] == 9
    # the derived factor 1 − c/L_ε, at the first two ε: every pair run in which the atoms formed is early by it to 0.02 (31 of 44 runs), and the entry times
    # are the path's up to that factor to 0.05. (The sample: the design and 10 family-A targets whose features are at least 2.2σ apart.)
    n_fin = 0
    for t in Ks:
        for r in t['rows']:
            if r['start'].startswith('pairs') and r['eps'] in P.EPS_K and not np.isnan(r['tau_rel']):
                n_fin += 1
                assert abs(r['tau_rel'] - r['tau_pred']) < 0.02 and r['tau_ratio'] < 0.05, (r['start'], r['eps'])
    assert n_fin == 31 and min(t['sep'] for t in Ks) > 2.2
    # over the design and the ten targets: depth 2 never selects at 10⁻⁶ and 10⁻¹², depth 3 always gives K units; the pairs lie between (the gap K ordered by depth)
    med = {}
    for j, (st, e) in enumerate(keys):
        rs = [t['rows'][j] for t in Ks]
        med[(st, e)] = np.median([x['gap_k'] for x in rs])
        if st == 'D2' and e in P.EPS_K:
            assert not any(x['set_equal'] or x['exact'] for x in rs)
        if st.startswith('pairs a₀ = 0.3') and e == 1e-6:
            assert not any(x['set_equal'] or x['exact'] for x in rs) and np.median([x['n_formed'] for x in rs]) >= 5
        if st in ('P0, α(0) = 0', 'D3') and e == 1e-12:
            assert all(x['set_equal'] and x['exact'] and x['n_formed'] == 3 for x in rs)
    for e in P.EPS_K:
        assert med[('D2', e)] < med[('pairs a₀ = 0.3', e)] < med[('pairs a₀ = 0.1', e)] < med[('D3', e)], e
    # P0 and D3 dormant units: the maxima (not only the medians) over the 11 targets
    mx = lambda st, e: max(x['dormant_disp'] for x in (t['rows'][keys.index((st, e))] for t in Ks))              # noqa: E731
    assert 0.5 < mx('pairs a₀ = 0.1', 1e-6) < 3.0 and mx('P0, α(0) = 0', 1e-6) < 0.6 and mx('D3', 1e-6) < 1.0


@pytest.mark.skipif(not os.path.exists(RAW), reason='the raw is committed with the report')
def test_the_raw_pins_part_K_further_down_in_eps():
    """The pairs at ε = 10⁻²⁴, 10⁻⁴⁸, 10⁻⁹⁶ on the 11 targets (the referee's own cells: set 8, 8, 9 and exact 11, 11, 11 of 11 at a₀ = 0.1; set 5, 6, 8 and exact 9,
    10, 11 at a₀ = 0.3): the selection improves as 1/L_ε, the derived factor is exact to 0.004 from 10⁻²⁴ on (0.012 at 10⁻⁶), and the gap K grows."""
    raw = E.load(RAW)
    if 'K' not in raw:
        pytest.skip('part K is merged when finished')
    Ks = [t for t in raw['K'] if t.get('kind', 'std') == 'std']
    rows = lambda st, e: [next(r for r in t['rows'] if (r['start'], r['eps']) == (st, e)) for t in Ks]       # noqa: E731
    want = {0.1: {1e-24: (8, 11), 1e-48: (8, 11), 1e-96: (9, 11)}, 0.3: {1e-24: (5, 9), 1e-48: (6, 10), 1e-96: (8, 11)}}
    for a0, w in want.items():
        for e, (n_set, n_exact) in w.items():
            rs = rows(f'pairs a₀ = {a0:g}', e)
            assert (sum(r['set_equal'] for r in rs), sum(r['exact'] for r in rs)) == (n_set, n_exact), (a0, e)
            fin = [r for r in rs if not np.isnan(r['tau_rel'])]
            assert len(fin) >= 6 and all(abs(r['tau_rel'] - r['tau_pred']) < 0.004 and r['tau_ratio'] < 0.007 for r in fin), (a0, e)
    for a0 in P.A0_K:
        gaps = [np.median([r['gap_k'] for r in rows(f'pairs a₀ = {a0:g}', e)]) for e in P.EPS_K + P.EPS_K_FAR]
        assert np.all(np.diff(gaps) > 0), (a0, gaps)                              # the gap K grows with L_ε
    assert 8.0 < np.median([r['gap_k'] for r in rows('pairs a₀ = 0.1', 1e-96)]) < 9.5
    d2 = [sum(r['set_equal'] for r in rows('D2', e)) for e in P.EPS_K + P.EPS_K_FAR]                       # depth 2 also improves, from none at 10⁻⁶ and 10⁻¹²
    assert d2[:2] == [0, 0] and d2[-1] > 0 and all(r['status'] == 'ok' for e in P.EPS_K_FAR for r in rows('D2', e))


@pytest.mark.skipif(not os.path.exists(RAW), reason='the raw is committed with the report')
def test_the_raw_pins_part_K_the_a0_scan_and_the_coupled_targets():
    raw = E.load(RAW)
    if 'K' not in raw:
        pytest.skip('part K is merged when finished')
    Ka = [t for t in raw['K'] if t.get('kind') == 'a0scan']
    Kc = [t for t in raw['K'] if t.get('kind') == 'coupled']
    assert len(Ka) == 1 and [(r['a0'], r['eps']) for r in Ka[0]['rows']] == [(a0, e) for a0 in P.A0_SCAN_K for e in P.EPS_SCAN_K]
    # a₀ interpolates between P0 and D2: the gap K falls with a₀ at every step up to a₀ = 1, from near P0's (6.15 at 10⁻⁶, 15.3 at 10⁻¹²) to near D2's (0.14, 0.68);
    # beyond, at the plateau, it moves by less than 0.01 up (0.121 → 0.125 at 10⁻⁶) or keeps falling (0.85 → 0.65 at 10⁻¹²)
    design = {(r['start'], r['eps']): r for t in raw['K'] if t.get('kind', 'std') == 'std' and t['i'] == 0 for r in t['rows']}
    for e in P.EPS_SCAN_K:
        g_ = [next(r['gap_k'] for r in Ka[0]['rows'] if (r['a0'], r['eps']) == (a0, e)) for a0 in P.A0_SCAN_K]
        k1 = P.A0_SCAN_K.index(1.0)
        assert np.all(np.diff(g_[:k1 + 1]) < 0) and np.all(np.diff(g_[k1:]) < 0.01), (e, g_)
        p0, d2 = design[('P0, α(0) = 0', e)]['gap_k'], design[('D2', e)]['gap_k']
        assert 0.5 * p0 < g_[0] < p0 and abs(g_[-1] - d2) < 0.3 and g_[0] > 5 * g_[k1]
    # the a₀ = 0.1, 0.3 cells of the scan are the design's cells
    for a0 in (0.1, 0.3):
        for e in P.EPS_SCAN_K:
            _close(next(r for r in Ka[0]['rows'] if (r['a0'], r['eps']) == (a0, e)), design[(f'pairs a₀ = {a0:g}', e)])
    # the coupled targets: 16 of C's family D(b), eight of each sign
    assert [t['i'] for t in Kc] == list(range(120, 136)) and sum(t['same_sign'] for t in Kc) == 8
    assert all(r['status'] == 'ok' and r['a0'] == 0.1 and len(r['earliness']) == 2 for t in Kc for r in t['rows'])
    for t in Kc:
        assert [r['eps'] for r in t['rows']] == list(P.EPS_KC) and len(t['path_tau']) == 2 and np.isclose(t['slope'][0] * t['path_tau'][0], 1.0, rtol=1e-3, atol=0)
        for r in t['rows']:
            pr = 1 - P.FR.PaperPairs(P.law(), 0.1).form_ratio(r['eps'], P.FR.formed_threshold([1.0, 1.0]))
            assert np.isclose(r['tau_pred'], pr, rtol=1e-12, atol=0)
            assert abs(r['earliness'][0] - pr) < 0.005, (t['i'], r['eps'])                       # the first entrant: the derived factor, on every target at every ε
            if t['same_sign']:
                assert r['earliness'][1] > 1.3 * pr and r['tau_ratio'] > 0.1, (t['i'], r['eps'])     # the second atom of a same-sign pair: early by far more
            else:
                assert abs(r['earliness'][1] - pr) < 0.02 and r['tau_ratio'] < 0.02, (t['i'], r['eps'])    # of an opposite-sign pair: by the prediction
    # target 121 (the referee's example): formed at τ = 3.05 and 3.35, entered by the path at 5.53 and 19.86 (ε = 10⁻⁶)
    t121 = next(t for t in Kc if t['i'] == 121)
    assert np.allclose(t121['path_tau'], [5.526, 19.864], atol=2e-3, rtol=0)
    assert np.allclose(t121['rows'][0]['form_times'], [3.047, 3.354], atol=2e-3, rtol=0)
    assert np.isclose(t121['slope'][1] * t121['path_tau'][1], 3.32, atol=0.02, rtol=0)           # the path enters it 3.3 times later than a straight dual would
    t126 = next(t for t in Kc if t['i'] == 126)
    assert np.isclose(t126['slope'][1] * t126['path_tau'][1], 1.0, atol=0.01, rtol=0)             # an opposite-sign pair's second dual is straight to its wall


@pytest.mark.skipif(not os.path.exists(RAW), reason='the raw is committed with the report')
def test_the_raw_pins_part_L_the_dictionary_spacing():
    raw = E.load(RAW)
    if 'L' not in raw:
        pytest.skip('part L is merged when finished')
    assert [(t['nd'], t['variant']) for t in raw['L']] == ([(n, 'off') for n in P.ND_L] + [(n, 'on') for n in P.ND_ON_L] + [(n, 'scan') for n in P.ND_SCAN_L])
    assert all(r['status'] == 'ok' for t in raw['L'] for r in t['rows'] + t.get('rows_long', []))
    Ls = [t for t in raw['L'] if t['variant'] in ('off', 'on')]
    assert [(r['readout'], r['eps']) for r in Ls[0]['rows']] == [(n, e) for n in P.NAMES_L for e in P.EPS_L]
    off = {t['nd']: t for t in Ls if t['variant'] == 'off'}
    on = {t['nd']: t for t in Ls if t['variant'] == 'on'}
    row = lambda t, nm, e: next(r for r in t['rows'] if r['readout'] == nm and r['eps'] == e)             # noqa: E731
    # the referee's validation: D2 forms 9, 8, 6, 5, 4, 4 units at 1e-6 and passes at 1e-12 only at spacing 0.33
    assert [row(off[n], 'D2', 1e-6)['n_formed'] for n in P.ND_L] == [9, 8, 6, 5, 4, 4]
    assert not any(row(off[n], 'D2', 1e-6)['set_equal'] or row(off[n], 'D2', 1e-6)['exact'] for n in P.ND_L)
    assert [row(off[n], 'D2', 1e-12)['set_equal'] for n in P.ND_L] == [False] * 5 + [True]
    assert [row(off[n], 'D2', 1e-12)['exact'] for n in P.ND_L] == [False] * 5 + [True]
    # P0 and D3 pass at every spacing, with a gap of 5 to 8 decades at 1e-6
    for n in P.ND_L:
        for nm in ('P0', 'D3'):
            for e in P.EPS_L:
                r = row(off[n], nm, e)
                assert r['set_equal'] and r['order_equal'] and r['exact'] and r['n_formed'] == 3, (n, nm, e)
            assert 5.0 < row(off[n], nm, 1e-6)['gap_k'] < 8.5, (n, nm)
    # with the features on kernel centres where neighbouring kernels overlap by at most 0.5 (11, 9, 7 kernels) depth 2 forms the path's set at both ε;
    # exact once the amplitudes have converged
    for n in (11, 9, 7):
        for nm in P.NAMES_L:
            for e in P.EPS_L:
                r = row(on[n], nm, e)
                assert r['set_equal'] and r['order_equal'] and r['n_formed'] == 3, (n, nm, e)
                assert r['exact'] == (nm != 'D2' or e == 1e-12), (n, nm, e)
        assert 2.5 < row(on[n], 'D2', 1e-6)['gap_k'] < 5.0
    # the on-grid arm at ALL six spacings (PZ2-5), a control that removes the sliding: P0 and D3 pass everywhere; D2 passes iff the margin rule holds, which
    # is not the case with 21 kernels (9 and 6 units at 1e-6 and 1e-12) nor with 16 at 1e-6 (9 units) — the referee's cells
    for n in P.ND_ON_L:
        for nm in ('P0', 'D3'):
            for e in P.EPS_L:
                r = row(on[n], nm, e)
                assert r['set_equal'] and r['order_equal'] and r['exact'] and r['n_formed'] == 3, (n, nm, e)
        for e in P.EPS_L:
            assert row(on[n], 'D2', e)['set_equal'] == P.rule_passes(on[n]['ref']['margin_never'], e), (n, e)
    assert [row(on[n], 'D2', 1e-6)['set_equal'] for n in P.ND_ON_L] == [False, False, True, True, True, True]
    assert [row(on[n], 'D2', 1e-12)['set_equal'] for n in P.ND_ON_L] == [False, True, True, True, True, True]
    assert [row(on[21], 'D2', 1e-6)['n_formed'], row(on[21], 'D2', 1e-12)['n_formed'], row(on[16], 'D2', 1e-6)['n_formed']] == [9, 6, 9]
    assert np.isclose(on[21]['ref']['margin_never'], 0.152, atol=1e-3, rtol=0) and np.isclose(on[16]['ref']['margin_never'], 0.262, atol=1e-3, rtol=0)
    # the margin that decides depth 2: the rule with the margin of the kernels that never enter predicts every one of the 24 D2 cells, with the later entrants' it does not
    d2 = [(t, r) for t in Ls for r in t['rows'] if r['readout'] == 'D2']
    assert len(d2) == 24
    assert all(P.rule_passes(t['ref']['margin_never'], r['eps']) == r['set_equal'] for t, r in d2)
    assert sum((t['ref']['margin'] > P.rule_threshold(r['eps'])) == r['set_equal'] for t, r in d2) == 18                         # the later entrants' margin: 6 of the 24 disagree
    assert all(t['ref']['margin_never'] >= t['ref']['margin'] for t in Ls)
    assert np.isclose(off[21]['ref']['margin'], 0.0664, atol=1e-3, rtol=0)                                  # the design's margin (the test above)
    # the on-grid depth-2 cells at 1e-6 that form the path's set are not exact at the study's horizon (the amplitudes 11% off, 5% allowed); at three times the
    # horizon every cell that forms the set is exact, and every D2 cell is as at the study's horizon
    for n in (13, 11, 9, 7):
        r = row(on[n], 'D2', 1e-6)
        assert not r['exact'] and 0.10 < r['amp_err'] < 0.12 and r['pos_err'] < 0.01
    for n in (16, 13, 11, 9, 7):
        assert 0.02 < row(on[n], 'D2', 1e-12)['amp_err'] < 0.03
    for t in Ls:
        assert [r['eps'] for r in t['rows_long']] == list(P.EPS_L) and all(r['readout'] == 'D2' and r['status'] == 'ok' for r in t['rows_long'])
        for r in t['rows_long']:
            r1 = row(t, 'D2', r['eps'])
            assert (r['n_formed'] == r1['n_formed']) or t['variant'] == 'on' and t['nd'] in (21, 16)               # the dictionaries on which D2 fails: 9 → 9, 6 → 7, 9 → 9
            assert r['set_equal'] == r1['set_equal']
            assert r['exact'] == (r1['set_equal'] and (t['variant'] == 'on' or (t['nd'] == 7 and r['eps'] == 1e-12)))
            if r['exact'] and t['variant'] == 'on':
                assert r['amp_err'] < 0.015 and r['pos_err'] < 0.01


@pytest.mark.skipif(not os.path.exists(RAW), reason='the raw is committed with the report')
def test_the_raw_pins_part_L_the_rule_out_of_sample_the_scan_of_epsilon():
    """D2 over ε = 10⁻⁶ … 10⁻⁹⁰ on 11, 13, 16, 21, 25, 31 kernels (the referee's scan): the ε from which it passes is where the rule puts it, to 0.90 … 1.02
    of L = ln(1/θ)/μ never, with no fitted parameter."""
    raw = E.load(RAW)
    if 'L' not in raw:
        pytest.skip('part L is merged when finished')
    scan = [t for t in raw['L'] if t['variant'] == 'scan']
    assert [t['nd'] for t in scan] == list(P.ND_SCAN_L)
    want = {11: (14, 1.02), 13: (18, 0.90), 16: (36, 1.00), 21: (28, 0.93), 25: (54, 0.97), 31: (42, 1.02)}            # k* and L_obs/L_rule
    for t in scan:
        assert [r['eps'] for r in t['rows']] == [float(f'1e-{k}') for k in P.K_SCAN_L] and all(r['readout'] == 'D2' and r['status'] == 'ok' for r in t['rows'])
        ks = [int(round(-np.log10(r['eps']))) for r in t['rows']]
        ok = [bool(r['set_equal']) for r in t['rows']]
        kstar = next(k for j, k in enumerate(ks) if all(ok[j:]))                            # D2 passes from 10^(−k*) on, and at every smaller ε of the scan
        l_rule = np.log(100) / P.rule_margin(t['ref']['margin_never'])
        assert kstar == want[t['nd']][0] and np.isclose(kstar * np.log(10) / l_rule, want[t['nd']][1], atol=0.006, rtol=0), t['nd']
        assert 0.85 < kstar * np.log(10) / l_rule < 1.1 and not ok[0]                       # not at 10⁻⁶
        for r in t['rows']:                                                                 # far from the transition the rule is exactly right
            if abs(np.log(1 / r['eps']) - l_rule) > 0.2 * l_rule:
                assert r['set_equal'] == P.rule_passes(t['ref']['margin_never'], r['eps']), (t['nd'], r['eps'])


@pytest.mark.skipif(not os.path.exists(RAW), reason='the raw is committed with the report')
def test_parts_K_and_L_are_counted_in_the_status_header():
    import copy
    raw = E.load(RAW)
    for key, kw in (('K', 'Kb'), ('L', 'Ls')):
        if key not in raw:
            continue
        part = copy.deepcopy(raw[key])
        n = sum(len(t['rows']) + len(t.get('rows_long', [])) for t in part)                   # (L also holds D2 at three times the horizon)
        assert f'{key}: {n} / 0' in P.report(**{kw: part})
        part[0]['rows'][0]['status'] = 'timeout'
        assert f'{key}: {n} / 1' in P.report(**{kw: part})


def _close(a, b, path=''):
    """Recursive equality of two cell records: flags, counts, strings and lists equal, floats to 1e-6 (secs and nfev are not compared)."""
    if isinstance(a, dict):
        assert set(a) == set(b), path
        for k in a:
            if k not in ('secs', 'nfev'):
                _close(a[k], b[k], f'{path}.{k}')
    elif isinstance(a, (list, tuple, np.ndarray)):
        assert len(a) == len(b), path
        for i, (x, y) in enumerate(zip(a, b)):
            _close(x, y, f'{path}[{i}]')
    elif isinstance(a, (float, np.floating)):
        assert (np.isnan(a) and np.isnan(b)) or np.isclose(a, b, rtol=1e-6, atol=1e-12), (path, a, b)
    else:
        assert a == b, (path, a, b)


@pytest.mark.skipif(not os.path.exists(RAW), reason='the raw is committed with the report')
def test_part_K_design_job_recomputes_the_raw():
    raw = E.load(RAW)
    if 'K' not in raw:
        pytest.skip('part K is merged when finished')
    got = P.job_k((0, dict(cs=np.array(P.CSTAR), As=np.array(P.ASTAR))))
    _close(got, raw['K'][0])


@pytest.mark.skipif(not os.path.exists(RAW), reason='the raw is committed with the report')
def test_part_L_seven_kernels_job_recomputes_the_raw():
    raw = E.load(RAW)
    if 'L' not in raw:
        pytest.skip('part L is merged when finished')
    got = P.job_l((7, 'off'))
    _close(got, next(t for t in raw['L'] if (t['nd'], t['variant']) == (7, 'off')))


@pytest.mark.skipif(not os.path.exists(RAW), reason='the raw is committed with the report')
def test_part_K_coupled_and_a0_scan_jobs_recompute_the_raw(monkeypatch):
    """The coupled job (target 121 at ε = 10⁻²⁴ only) and the a₀-scan job (a₀ = 0.3 at 10⁻⁶ only) on reduced grids equal the raw's cells, with the path's
    entries, the features' separation and sign, and the duals' initial slopes the coupled job stores."""
    raw = E.load(RAW)
    if 'K' not in raw or not any(t.get('kind') == 'coupled' for t in raw['K']):
        pytest.skip('part K is merged when finished')
    cp = dict(P.coupled_targets())
    monkeypatch.setattr(P, 'EPS_KC', (1e-24,))
    got = P.job_kc((121, cp[121]))
    want = next(t for t in raw['K'] if t.get('kind') == 'coupled' and t['i'] == 121)
    assert got['kind'] == 'coupled' and got['same_sign'] is True and np.isclose(got['sep'], 1.5, rtol=1e-12, atol=0)
    assert [r['eps'] for r in got['rows']] == [1e-24] and got['path_order'] == want['path_order'] and len(got['path_tau']) == 2
    _close(got['rows'][0], next(r for r in want['rows'] if r['eps'] == 1e-24))
    for k in ('margin', 'path_tau', 'slope', 'sep'):
        _close(got[k], want[k])
    assert got['slope'][0] > got['slope'][1] > 0 and np.isclose(got['slope'][0] * got['path_tau'][0], 1.0, rtol=1e-3, atol=0)     # the first entrant: a straight dual to its wall
    monkeypatch.setattr(P, 'A0_SCAN_K', (0.3,))
    monkeypatch.setattr(P, 'EPS_SCAN_K', (1e-6,))
    ga = P.job_ka((0, dict(cs=np.array(P.CSTAR), As=np.array(P.ASTAR))))
    wa = next(t for t in raw['K'] if t.get('kind') == 'a0scan')
    assert ga['kind'] == 'a0scan' and ga['i'] == 0 and [(r['start'], r['eps']) for r in ga['rows']] == [('pairs a₀ = 0.3', 1e-6)]
    _close(ga['rows'][0], next(r for r in wa['rows'] if (r['start'], r['eps']) == ('pairs a₀ = 0.3', 1e-6)))


@pytest.mark.skipif(not os.path.exists(RAW), reason='the raw is committed with the report')
def test_part_L_scan_job_recomputes_the_raw(monkeypatch):
    """The scan job of an 11-kernel dictionary on a reduced grid (ε = 10⁻⁶, 10⁻¹², 10⁻²⁴) equals the raw's cells; the scan's features are the design's."""
    raw = E.load(RAW)
    if 'L' not in raw or not any(t['variant'] == 'scan' for t in raw['L']):
        pytest.skip('part L is merged when finished')
    monkeypatch.setattr(P, 'K_SCAN_L', (6, 12, 24))
    got = P.job_l((11, 'scan'))
    want = next(t for t in raw['L'] if (t['nd'], t['variant']) == (11, 'scan'))
    assert got['variant'] == 'scan' and [r['eps'] for r in got['rows']] == [1e-6, 1e-12, 1e-24] and np.array_equal(got['cs'], np.array(P.CSTAR))
    for r in got['rows']:
        _close(r, next(x for x in want['rows'] if x['eps'] == r['eps']))
    _close(got['ref'], want['ref'])
    assert 'rows_long' not in got and np.isclose(got['spacing'], 0.2, rtol=1e-12, atol=0)


def test_the_job_lists_of_parts_K_and_L_and_the_dispatcher(monkeypatch):
    """Part K: the 11 main jobs, then the 16 coupled targets, then the a₀ scan; part L: 6 off-grid, 6 on-grid and 6 scan dictionaries."""
    seen = {}

    class FakePool:
        def __init__(self, n):
            seen['n'] = n

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

        def map(self, fn, jobs, chunksize=None):
            seen['fn'], seen['jobs'], seen['chunksize'] = fn, list(jobs), chunksize
            return ['result']
    import multiprocessing
    monkeypatch.setattr(multiprocessing, 'Pool', FakePool)
    assert P.part_K() == ['result'] and seen['fn'] is P.run_k and seen['chunksize'] == 1 and seen['n'] == P.N_WORKERS
    kinds = [k for k, _ in seen['jobs']]
    assert kinds == ['std'] * 11 + ['coupled'] * 16 + ['a0scan']
    assert [a[0] for k, a in seen['jobs'] if k == 'std'] == list(range(11)) and [a[0] for k, a in seen['jobs'] if k == 'coupled'] == list(range(120, 136))
    assert np.allclose(seen['jobs'][0][1][1]['cs'], P.CSTAR) and seen['jobs'][-1][1][0] == 0 and np.allclose(seen['jobs'][-1][1][1]['As'], P.ASTAR)
    assert P.part_L() == ['result'] and seen['fn'] is P.job_l and seen['chunksize'] == 1
    assert seen['jobs'] == [(n, 'off') for n in P.ND_L] + [(n, 'on') for n in P.ND_L] + [(n, 'scan') for n in P.ND_SCAN_L]
    monkeypatch.setattr(P, 'job_k', lambda a: ('k', a))
    monkeypatch.setattr(P, 'job_kc', lambda a: ('kc', a))
    monkeypatch.setattr(P, 'job_ka', lambda a: ('ka', a))
    assert P.run_k(('std', 1)) == ('k', 1) and P.run_k(('coupled', 2)) == ('kc', 2) and P.run_k(('a0scan', 3)) == ('ka', 3)


@pytest.mark.skipif(not os.path.exists(RAW), reason='the raw is committed with the report')
def test_the_md_header_gives_the_radau_counts_of_the_parts_that_record_them():
    """PZ2-10: the header's Radau count was Part C's alone (104); K's and L's (and B's, D's, E's, J's) are in it, counted from the raw."""
    raw = E.load(RAW)
    head = next(l for l in open(MD, encoding='utf-8').read().splitlines() if l.startswith('Engine:'))
    assert 'Cells with a piece redone by Radau' in head and 'C: 104 of 1200' in head
    cells = dict(C=[c for r in raw['C'] for c in r['cells']], D=[c for t in raw['D'] for c in t['rows'] + t['rows_long']],
                 K=[r for t in raw['K'] for r in t['rows']], L=[r for t in raw['L'] for r in t['rows'] + t.get('rows_long', [])])
    for part, cs_ in cells.items():
        assert f"{part}: {sum(c['radau'] > 0 for c in cs_)} of {len(cs_)}" in head, part
    assert sum(c['radau'] > 0 for c in cells['K']) > 0 and sum(c['radau'] > 0 for c in cells['L']) > 0
    assert P.report(B=raw['B']).count('Radau (the raw') == 1 and 'B: 2 of 40' in P.report(B=raw['B'])          # a part alone: only its own count


@pytest.mark.skipif(not os.path.exists(RAW), reason='the raw is committed with the report')
def test_the_report_and_the_engine_do_not_claim_the_common_factor_for_every_entry():
    """PZ2-1: the factor 1 − c/L_ε is the first entrant's (exact), a later atom's only while its dual is straight; part K's intro says that the sliding
    extension is tested, not proved, and that ε = 10⁻¹² is pre-asymptotic; the K sample's figures are stated with their sample."""
    text = open(MD, encoding='utf-8').read()
    k = text[text.index('## K.'):text.index('## L.')]
    assert 'tested, not proved' in k and 'is pre-asymptotic' in k and 'claimed for the first entrant' in k
    assert 'is early by the same factor only while its dual is straight until its wall' in k and 'both figures belong to this sample' in k
    assert 'Coupled targets: where the common factor holds and where it fails' in k and 'same sign' in k
    d = text[text.index('## D.'):text.index('## E.')]
    assert 'first entrant' in d and 'only while its dual is straight until its wall' in d and 'μ never' in d
    assert 'by which every entry is early' not in text and 'every entry is early' not in text
    doc = P.FR.PaperPairs.form_ratio.__doc__
    assert 'FIRST entrant' in doc and 'NOT the earliness of every entry' in doc and 'every entry of the network is early by this factor' not in doc
    assert 'FIRST entrant' in P.FR.PaperPairs.__doc__


def test_the_parts_and_the_report_arguments_correspond():
    import inspect
    assert P.PARTS == 'ABCDEFGHIJKL' and len(inspect.signature(P.report).parameters) == len(P.PARTS)
    assert "zip(PARTS, (part_A, part_B, part_C, part_D, part_E, part_F, part_G, part_H, part_I, part_J, part_K, part_L))" in inspect.getsource(P.main)
    assert [P.part_K.__name__, P.part_L.__name__] == ['part_K', 'part_L']
