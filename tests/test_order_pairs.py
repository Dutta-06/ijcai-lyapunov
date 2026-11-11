"""Tests of experiments/exp_order_pairs.py: the scene, the classification of an outcome, the Wilson interval, a tiny end-to-end run and its report."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'experiments'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import exp_order_pairs as X      # noqa: E402


def test_scene_has_the_requested_separation_ratio_and_is_deterministic():
    p1, w1 = X.scene(np.random.default_rng(3), 2.0, 3.0, 6000.0)
    p2, w2 = X.scene(np.random.default_rng(3), 2.0, 3.0, 6000.0)
    assert np.array_equal(p1, p2) and np.array_equal(w1, w2)
    assert abs(np.linalg.norm(p1[0] - p1[1]) - 2.0 * X.SIG) < 1e-12 and np.allclose(w1, [6000.0, 2000.0])
    assert np.all((p1 > 8) & (p1 < 24))


def test_classify_the_outcomes():
    pos = np.array([[10.0, 10.0], [13.0, 10.0]])
    assert X.classify(np.array([[10.2, 10.1], [12.9, 9.9]]), np.ones(2), pos) == 1                 # resolved
    assert X.classify(np.array([[12.9, 9.9], [10.2, 10.1]]), np.ones(2), pos) == 1                 # in the other order
    assert X.classify(np.array([[11.5, 10.0]]), np.ones(1), pos) == 2                              # merged
    assert X.classify(np.array([[10.0, 10.0], [13.0, 10.0], [11.0, 10.0]]), np.ones(3), pos) == 3  # split
    assert X.classify(np.array([[10.0, 10.0], [18.0, 10.0]]), np.ones(2), pos) == 0                # one unit far off
    assert X.classify(np.zeros((0, 2)), np.zeros(0), pos) == 0
    assert X.classify(np.array([[10.0, 10.0], [10.5, 10.0]]), np.ones(2), pos) == 0                # both on the first emitter: not one-to-one


def test_wilson_interval_contains_the_proportion_and_handles_edges():
    lo, hi = X.wilson(30, 100)
    assert lo < 0.3 < hi and 0.2 < lo and hi < 0.4
    lo0, hi0 = X.wilson(0, 50)
    assert lo0 == 0.0 or abs(lo0) < 1e-12
    assert hi0 > 0 and np.isnan(X.wilson(0, 0)[0])


def test_tiny_run_and_report_end_to_end(tmp_path):
    X.SEPS, X.RATIOS, X.PHOTONS = (3.0,), (1.0,), (6000.0,)
    path = str(tmp_path / 'raw.npz')
    X.run(1, 2, ('bic', 4.0), out_path=path)
    z = np.load(path)
    assert z['cls'].shape == (1, 1, 1, 2, 2)
    md = X.report(path, boot=50)
    assert md == X.report(path, boot=50) and '| separation (widths)' in md and 'resolved' in md


def test_clean_up_option_is_stored_and_reported_and_changes_nothing_when_absent(tmp_path):
    path = str(tmp_path / 'raw_clean.npz')
    X.run(1, 1, ('bic', 4.0), out_path=path, clean=5.27)
    z = dict(np.load(path))
    assert abs(float(z['clean'][0]) - 5.27) < 1e-12
    assert 'then the units of own significance below 5.27 discarded' in X.report(path)
    path2 = str(tmp_path / 'raw_noclean.npz')
    X.run(1, 1, ('bic', 4.0), out_path=path2)
    z2 = dict(np.load(path2))
    assert np.isnan(float(z2['clean'][0])) and 'discarded' not in X.report(path2)
