"""Tests of experiments/exp_kan_merge.py (the POST HOC merge analysis of the KAN spike): the headline numbers of kan.md reproduced from the stored raw, the merge
sweep invariants, the strata, the duplicate statistics, a planted duplicate pair, the bootstrap, and the byte-identical regeneration of kan_merge.md."""
import json
import os
import subprocess
import sys

import numpy as np
import pytest

HERE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'experiments')
sys.path.insert(0, HERE)
import exp_kan_merge as M  # noqa: E402
import exp_kan_spike as K  # noqa: E402

SCRIPT = M.__file__
RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results', 'ml')
RAW = os.path.join(RES, 'kan_raw.npz')


@pytest.fixture(scope='module')
def raw():
    meta, recs = K.load_raw(RAW)
    return meta, recs


@pytest.fixture(scope='module')
def comp(raw):
    return M.compute(raw[1])


# ---------------------------------------------------------------------------------------------------------------------------------
# registries
# ---------------------------------------------------------------------------------------------------------------------------------
def test_registries():
    assert M.DELTAS == (0.0, 0.1, 0.2, 0.3, 0.4, 0.5)
    assert M.GAP_THRESH == (0.1, 0.3) and M.MATCH_TOLS == (0.2, 0.5)
    assert (M.N_BOOT, M.BOOT_SEED) == (20000, 0)
    assert M.CKEY == 'centres_0.02' and M.NKEY == 'counts_0.02' and K.TAU == 0.02


# ---------------------------------------------------------------------------------------------------------------------------------
# merging, planted
# ---------------------------------------------------------------------------------------------------------------------------------
def test_merged_count_planted_pair():
    pair = [1.0, 1.25]                                            # a duplicate pair, gap 0.25
    assert M.merged_count(pair, 0.0) == 2 and M.merged_count(pair, 0.1) == 2 and M.merged_count(pair, 0.2) == 2
    assert M.merged_count(pair, 0.3) == 1 and M.merged_count(pair, 0.5) == 1
    assert M.merged_count([1.25, 1.0], 0.3) == 1                  # order of the stored centres is irrelevant
    assert M.merged_count([1.25, 1.0], 0.2) == 2
    assert M.merged_count([1.0, 1.0], 0.0) == 2                   # delta = 0 merges nothing, even identical centres
    assert M.merged_count([1.0, 1.0], 0.1) == 1
    assert M.merged_count([], 0.5) == 0 and M.merged_count([0.3], 0.5) == 1


def test_merged_count_single_linkage_chain():
    chain = [0.0, 0.2, 0.4, 2.0]                                  # 0.0-0.2-0.4 are one cluster by single linkage although 0.0 and 0.4 are 0.4 apart
    assert M.merged_count(chain, 0.25) == 2
    assert M.merged_count(chain, 0.15) == 4
    assert M.merged_count(chain, 0.45) == 2
    assert M.merged_count(chain, 2.0) == 1
    assert M.gaps([0.4, 0.0, 0.2]) == [0.2, 0.2]


# ---------------------------------------------------------------------------------------------------------------------------------
# planted runs
# ---------------------------------------------------------------------------------------------------------------------------------
def rec(method, centres, tid, noise, seed, truth, truth_c, conf='c', val=0.1, rmse=0.05, P=None, steps=100):
    r = dict(kind='synth', tid=tid, noise=noise, seed=seed, method=method, conf=conf, val_mse=val, test_rmse=rmse, n_val=500, secs=1.0, steps=steps, arch='A',
             truth=list(truth), truth_c=truth_c)
    r['centres_0.02'] = [list(c) for c in centres]
    r['counts_0.02'] = [len(c) for c in centres]
    n = sum(len(c) for c in centres)
    for tau in K.TAUS:
        r['P_%g' % tau] = 3 * n + 5 if P is None else P
    return r


TRUTH = [0, 1, 2]
TRUTH_C = [[], [0.0], [-1.5, 1.5]]


def planted_run(tid, noise, seed, ours_c, ctrl_cands, rmse_o=0.07, rmse_c=0.04):
    out = [rec('ours', ours_c, tid, noise, seed, TRUTH, TRUTH_C, conf='o', rmse=rmse_o),
           rec('fixed', [[], [0.0], [-1.5, 1.5]], tid, noise, seed, TRUTH, TRUTH_C, conf='f', val=0.5)]
    for i, (c, val) in enumerate(ctrl_cands):
        out.append(rec('ctrl', c, tid, noise, seed, TRUTH, TRUTH_C, conf='lam%d' % i, val=val, rmse=rmse_c))
    return out


EXACT = [[], [0.0], [-1.5, 1.5]]
DUP = [[], [0.0, 0.05], [-1.5, -1.45, 1.5, 1.6]]              # every non-empty edge doubled by a pair 0.05 / 0.05 / 0.1 apart


def planted(n_targets=3):
    recs = []
    for t in range(1, n_targets + 1):
        for seed in (0, 1):
            recs += planted_run(t, 0.1, seed, EXACT, [(DUP, 0.10), (DUP, 0.101)])
    return recs


def test_planted_sweep_and_strata():
    recs = planted()
    runs = M.select_runs(recs)
    assert len(runs) == 6
    T = M.truth_matrix(runs)
    assert T.shape == (6, 3) and T.sum() == 18
    assert (M.counts_matrix(runs, 'ours', 0.0) == T).all()
    c0 = M.counts_matrix(runs, 'ctrl', 0.0)
    assert c0.tolist()[0] == [0, 2, 4]
    assert (M.counts_matrix(runs, 'ctrl', 0.06)[:, 1] == 1).all() and (M.counts_matrix(runs, 'ctrl', 0.06)[:, 2] == 3).all()   # gaps 0.05, 0.05, 0.1
    assert (M.counts_matrix(runs, 'ctrl', 0.11)[:, 2] == 2).all()                                                                # the 0.1 gap merged
    assert M.counts_matrix(runs, 'ctrl', 0.1).tolist()[0] == [0, 1, 3]                                                          # gap 0.1 is not below 0.1
    st = M.strata(T, c0)
    assert st['empty'] == dict(n=6, exact_n=6, share=1.0, over_run=0.0, under_run=0.0, over_edge=0.0)
    ne = st['nonempty']
    assert ne['n'] == 12 and ne['exact_n'] == 0 and ne['share'] == 0.0 and ne['over_run'] == 3.0 and ne['over_edge'] == 1.5 and ne['under_run'] == 0.0
    # the strata add up to the all-edge share
    assert st['empty']['exact_n'] + ne['exact_n'] == int((c0 == T).sum())
    assert st['empty']['n'] + ne['n'] == T.size
    s_ours = M.strata(T, M.counts_matrix(runs, 'ours', 0.0))
    assert s_ours['nonempty']['share'] == 1.0 and s_ours['nonempty']['over_run'] == 0.0
    # an under-counted edge
    st2 = M.strata(np.array([[0, 2, 1]]), np.array([[0, 1, 3]]))
    assert st2['nonempty']['over_run'] == 2.0 and st2['nonempty']['under_run'] == 1.0 and st2['empty']['over_run'] == 0.0
    st3 = M.strata(np.array([[0, 1]]), np.array([[2, 1]]))
    assert st3['empty']['over_run'] == 2.0 and st3['empty']['share'] == 0.0 and st3['nonempty']['share'] == 1.0


def test_planted_compute_gain_and_oracle():
    res = M.compute(planted())
    sw = res['sweep']
    # unmerged: ours exact everywhere, the control exact only on the empty edge
    assert np.allclose(sw[0.0]['share']['ours'], 1.0) and np.allclose(sw[0.0]['share']['ctrl'], 1 / 3)
    g = sw[0.0]['gain']['ctrl']
    assert g[0] == pytest.approx(2 / 3) and g[1] == pytest.approx(2 / 3) and g[2] == pytest.approx(2 / 3)      # identical gain in every target: degenerate bootstrap
    # at delta 0.2 every pair is merged: no gain left
    assert np.allclose(sw[0.2]['share']['ctrl'], 1.0) and sw[0.2]['gain']['ctrl'] == (0.0, 0.0, 0.0)
    # at delta 0.1 the 0.1 pair (gap exactly 0.1) is not yet merged, the two 0.05 pairs are
    assert np.allclose(sw[0.1]['share']['ctrl'], 2 / 3)
    assert sw[0.1]['all5']['ctrl'] == 0.0 and sw[0.2]['all5']['ctrl'] == 1.0
    assert sw[0.0]['units']['ours'] == 3.0 and sw[0.0]['units']['ctrl'] == 6.0
    assert sw[0.1]['nonempty']['ctrl'] == 0.5 and sw[0.1]['nonempty']['ours'] == 1.0


def test_gain_ci_cluster_bootstrap():
    runs = [dict(key=(t, 0.1, s)) for t in (1, 2, 3, 4) for s in (0, 1)]
    a = np.array([1.0, 1.0, 0.5, 0.5, 0.0, 0.0, 1.0, 1.0])
    b = np.zeros(8)
    m, lo, hi = M.gain_ci(runs, a, b)
    assert m == pytest.approx(0.625) and lo <= m <= hi
    by = [1.0, 0.5, 0.0, 1.0]
    assert (m, lo, hi) == K.boot_ci(by, n_boot=20000, seed=0)
    assert M.gain_ci(runs, a, b) == (m, lo, hi)                                                # deterministic
    # clustered over targets, not over runs: the pair of runs of a target moves together
    rng = np.random.default_rng(0).integers(0, 4, size=(20000, 4))
    assert lo == pytest.approx(float(np.quantile(np.array(by)[rng].mean(1), 0.025)))


def test_oracle_merged_selection():
    good = [[], [0.0], [-1.5, 1.5]]
    recs = planted_run(1, 0.1, 0, EXACT, [(DUP, 0.10), ([[], [0.0, 0.5], [-1.5, 1.5]], 0.2), (good, 0.3)])
    runs = M.select_runs(recs)
    # the validation control (one-SE rule: fewest active parameters among val <= best (1 + sqrt(2/500)) = 0.1 * 1.063) is the first one (DUP; the second has val 0.2)
    assert runs[0]['ctrl']['conf'] == 'lam0'
    # at delta 0 the oracle picks the exactly right penalty, in the merged world the first two both reach it, the tie is broken by fewer active parameters
    assert runs[0]['ctrl_or']['conf'] == 'lam2'
    C = M.oracle_merged(runs, 0.0)
    assert C.tolist() == [[0, 1, 2]]
    C = M.oracle_merged(runs, 0.2)
    assert C.tolist() == [[0, 1, 2]]
    # a run whose best penalty is the one with the extra unit merged at 0.6 only
    recs = planted_run(1, 0.1, 0, EXACT, [(DUP, 0.10), ([[], [0.0, 0.5], [-1.5, 1.5]], 0.2)])
    runs = M.select_runs(recs)
    assert M.oracle_merged(runs, 0.1).tolist() == [[0, 2, 2]]                                      # delta 0.1: DUP -> [0,1,3], other -> [0,2,2]: both 2/3 exact, tie broken by fewer active parameters
    assert M.oracle_merged(runs, 0.11).tolist() == [[0, 1, 2]]                                     # DUP is fully merged and exact
    assert M.oracle_merged(runs, 0.6).tolist() == [[0, 1, 2]]


def test_select_runs_mirrors_spike_module(raw):
    recs = raw[1]
    runs = M.select_runs(recs)
    rows = K.synth_summary(recs)
    assert [r['key'] for r in runs] == [r['key'] for r in rows] and len(runs) == 60
    for r, row in zip(runs, rows):
        assert r['ours']['conf'] == row['ours']['conf'] and r['ctrl']['conf'] == row['ctrl']['conf'] and r['ctrl_or']['conf'] == row['ctrl_or']['conf']
        assert len(r['ctrl_all']) == 6
        assert sum(1 for x in raw[1] if x['kind'] == 'synth' and x['method'] == 'ours' and K.run_key(x) == r['key']) == 1
    # an incomplete run (no fixed-grid record) is not scored, as in the module
    inc = planted()
    inc = [r for r in inc if not (r['method'] == 'fixed' and r['tid'] == 1 and r['seed'] == 0)]
    assert len(M.select_runs(inc)) == 5 and len(K.synth_summary(inc)) == 5


# ---------------------------------------------------------------------------------------------------------------------------------
# the stored raw: kan.md headline and the sweep invariants
# ---------------------------------------------------------------------------------------------------------------------------------
def test_headline_numbers_of_kan_md(comp, raw):
    sw = comp['sweep']
    assert round(float(np.mean(sw[0.0]['share']['ours'])), 3) == 0.723
    assert round(float(np.mean(sw[0.0]['share']['ctrl'])), 3) == 0.283
    assert round(float(np.mean(sw[0.0]['share']['ctrl_or'])), 3) == 0.287
    kan_md = open(os.path.join(RES, 'kan.md')).read()
    for s in ('| ours (factorised, small start) | 0.723 |', '| control: L1 + prune, penalty by validation | 0.283 |', '| control: penalty by the truth (oracle) | 0.287 |'):
        assert s in kan_md
    # exactly the published numbers, not only rounded
    rows = K.synth_summary(raw[1])
    for nm in ('ours', 'ctrl', 'ctrl_or'):
        assert float(np.mean(sw[0.0]['share'][nm])) == pytest.approx(float(np.mean([r[nm]['share'] for r in rows])), abs=1e-12)
    # counts from the stored centres are the stored counts
    for r in comp['runs']:
        assert M.counts_matrix([r], 'ours', 0.0).tolist()[0] == r['ours']['counts_0.02']
        assert M.counts_matrix([r], 'ctrl', 0.0).tolist()[0] == r['ctrl']['counts_0.02']
        assert M.counts_matrix([r], 'ctrl_or', 0.0).tolist()[0] == r['ctrl_or']['counts_0.02']
    assert comp['n_checked'] == 60 * 7


def test_referee_numbers(comp):
    T = comp['T']
    assert T.shape == (60, 5) and int((T == 0).sum()) == 84 and int((T > 0).sum()) == 216
    sw = comp['sweep']
    st = {nm: M.strata(T, sw[0.0]['C'][nm]) for nm in ('ours', 'ctrl')}
    assert st['ours']['empty']['exact_n'] == 84 and st['ctrl']['empty']['exact_n'] == 84
    assert st['ctrl']['nonempty']['exact_n'] == 1 and st['ours']['nonempty']['exact_n'] == 133
    assert round(st['ctrl']['nonempty']['over_run'], 2) == 8.67 and st['ctrl']['empty']['over_run'] == 0.0
    assert round(st['ours']['nonempty']['over_run'], 2) == 1.78
    assert round(sw[0.0]['all5']['ours'] * 60) == 14 and sw[0.0]['all5']['ctrl'] == 0.0
    n, sh = M.dup_stats(comp['runs'], 'ctrl')
    assert n == 784 and round(sh[0], 3) == 0.406 and round(sh[1], 3) == 0.541
    n, sh = M.dup_stats(comp['runs'], 'ours')
    assert round(sh[0], 3) == 0.049 and round(sh[1], 3) == 0.126
    rc, pc = M.match_stats(comp['runs'], 'ctrl', 0.2)
    assert round(rc, 3) == 0.986 and round(pc, 3) == 0.882
    ro, po = M.match_stats(comp['runs'], 'ours', 0.2)
    assert round(ro, 3) == 0.930 and round(po, 3) == 0.832
    # the sweep as published in the referee report
    exp = {0.1: (0.750, 0.563), 0.2: (0.783, 0.687), 0.3: (0.810, 0.767), 0.4: (0.830, 0.840), 0.5: (0.863, 0.910)}
    for d, (o, c) in exp.items():
        assert round(float(np.mean(sw[d]['share']['ours'])), 3) == o and round(float(np.mean(sw[d]['share']['ctrl'])), 3) == c
    assert round(float(np.mean(sw[0.5]['share']['ctrl_best'])), 3) == 0.957


def test_sweep_invariants(comp):
    sw, T = comp['sweep'], comp['T']
    for nm in ('ours', 'ctrl', 'ctrl_or'):
        assert (sw[0.0]['C'][nm] == np.array([r[nm]['counts_0.02'] for r in comp['runs']])).all()         # delta = 0 is the unmerged count
        prev = sw[0.0]['C'][nm]
        for d in M.DELTAS[1:]:
            cur = sw[d]['C'][nm]
            assert (cur <= prev).all()                                                                    # merging never increases the count
            assert ((cur == 0) == (prev == 0)).all()                                                      # and never removes the last unit of an edge
            prev = cur
        assert (sw[0.5]['C'][nm] <= sw[0.0]['C'][nm]).all() and (sw[0.5]['C'][nm] >= (sw[0.0]['C'][nm] > 0)).all()
    # the best-penalty control is at least as exact as the others at every delta
    for d in M.DELTAS:
        assert np.all(sw[d]['share']['ctrl_best'] >= sw[d]['share']['ctrl'] - 1e-12) and np.all(sw[d]['share']['ctrl_best'] >= sw[d]['share']['ctrl_or'] - 1e-12)
        for nm in ('ours', 'ctrl'):
            assert np.allclose(sw[d]['share'][nm], (sw[d]['C'][nm] == T).mean(1))
            assert sw[d]['units'][nm] == pytest.approx(float(sw[d]['C'][nm].sum(1).mean()))
    # strata add up at every delta
    for d in M.DELTAS:
        for nm in ('ours', 'ctrl'):
            C = sw[d]['C'][nm]
            st = M.strata(T, C)
            assert st['empty']['n'] + st['nonempty']['n'] == T.size == 300
            assert st['empty']['exact_n'] + st['nonempty']['exact_n'] == int((C == T).sum())
            assert (st['empty']['exact_n'] + st['nonempty']['exact_n']) / 300 == pytest.approx(float(np.mean(sw[d]['share'][nm])))
            assert st['nonempty']['share'] == pytest.approx(sw[d]['nonempty'][nm])
    # per-K strata add up as well
    ex = sw[0.0]['C']['ctrl'] == T
    assert sum(int(ex[T == k].sum()) for k in range(5)) == int(ex.sum())
    assert sum(int((T == k).sum()) for k in range(5)) == 300


# ---------------------------------------------------------------------------------------------------------------------------------
# duplicate statistics, planted
# ---------------------------------------------------------------------------------------------------------------------------------
def test_dup_and_match_planted():
    runs = M.select_runs(planted(1))
    n, sh = M.dup_stats(runs, 'ours')
    assert n == 2 * 1 and sh == [0.0, 0.0]                                                            # gaps 3.0 only (K=2 edge): one per run, two runs
    # ctrl = the DUP rows: gaps 0.05, then 0.05, 3.05 (-1.45 -> 1.5), 0.1 within the K=2 edge  [-1.5,-1.45,1.5,1.6]: 0.05, 2.95, 0.1
    runs_c = [dict(r, ctrl=dict(r['ours'], **{'centres_0.02': DUP})) for r in runs]
    n, sh = M.dup_stats(runs_c, 'ctrl')
    assert n == 2 * 4 and sh[0] == pytest.approx(4 / 8) and sh[1] == pytest.approx(6 / 8)              # gaps per run: 0.05, 0.05, 2.95, 0.1 -> below 0.1: 2; below 0.3: 3 (0.1 < 0.1 false)
    rc, pc = M.match_stats(runs, 'ours', 0.2)
    assert rc == 1.0 and pc == 1.0
    rc, pc = M.match_stats(runs_c, 'ctrl', 0.2)
    assert rc == 1.0 and pc == 1.0                                                                    # DUP: units 0.0, 0.05, -1.5, -1.45, 1.5, 1.6 are all within 0.2 of a bump
    assert M.match_stats(runs_c, 'ctrl', 0.04)[1] == pytest.approx(3 / 6)                             # within 0.04: 0.0, -1.5, 1.5 only
    assert M.match_stats(runs_c, 'ctrl', 0.04)[0] == 1.0


# ---------------------------------------------------------------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------------------------------------------------------------
def test_first_delta():
    sw = {d: dict(v=d) for d in M.DELTAS}
    assert M.first_delta(sw, lambda s: s['v'] >= 0.3) == 0.3
    assert M.first_delta(sw, lambda s: False) is None


def test_report_content_and_determinism(raw):
    meta, recs = raw
    txt = M.report(meta, recs)
    assert txt == M.report(meta, recs)
    for s in ('POST HOC', 'Reproduction of the published headline', '(match: yes)', '0.440 [0.363, 0.517]', '| 0.3 | 0.810 | 0.767 | 0.043 [-0.030, 0.120] |',
              'Clean-test RMSE: ours 0.068, control (validation) 0.039', 'ours reached the 12000-step cap in 32 of 60 runs',
              'point estimate first is not positive at delta = 0.4', 'first includes 0 at delta = 0.3', 'shrinks to 0.187 at delta = 0.1 and to 0.097 at delta = 0.2',
              'verdict KILL, non-inferior on 1 of 8 sets', '84 of the 300 edges have K = 0', '0.957'):
        assert s in txt, s


def test_report_byte_identity(tmp_path, raw):
    meta, recs = raw
    txt = M.report(meta, recs)
    repo_md = open(os.path.join(RES, 'kan_merge.md')).read()
    assert repo_md == txt + '\n'
    out1, out2 = tmp_path / 'a.md', tmp_path / 'b.md'
    for o in (out1, out2):
        subprocess.run([sys.executable, '-I', SCRIPT, '--from-raw', RAW, '--md', str(o)], check=True)
    assert out1.read_bytes() == out2.read_bytes() == (txt + '\n').encode()
    # the harness form: no path, the output directory from the environment
    d = tmp_path / 'o'
    d.mkdir()
    import shutil
    shutil.copy(RAW, d / 'kan_raw.npz')
    subprocess.run([sys.executable, '-I', SCRIPT, '--from-raw'], check=True, env=dict(os.environ, KAN_OUT=str(d)))
    assert (d / 'kan_merge.md').read_bytes() == (txt + '\n').encode()


def test_consistency_check_catches_mismatch():
    recs = planted(1)
    runs = M.select_runs(recs)
    assert M.check_consistency(runs) == 2 * 3
    bad = [dict(r) for r in recs]
    for r in bad:
        if r['method'] == 'ours':
            r['counts_0.02'] = [0, 1, 3]
    with pytest.raises(AssertionError):
        M.check_consistency(M.select_runs(bad))
