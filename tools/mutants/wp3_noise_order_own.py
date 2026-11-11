"""Single-point mutants of experiments/exp_noise_order.py (the target generator, the scoring, the noise draw, the Wilson interval and the report's rates, biases and strata), for tools/mutants.py:

    venv/bin/python tools/mutants.py --exp exp_noise_order --tests tests/test_noise_order.py --raw results/chaos/noise_order_raw.npz --md results/chaos/noise_order.md \\
        --out-env NO_OUT --spec tools/mutants/wp3_noise_order_own.py --work <scratch dir> [--slot a --i0 0 --step 2]

A tuple is (id, what, old, new[, k]): `old` replaced by `new` at its k-th occurrence in experiments/exp_noise_order.py. The kill criterion is the tests file alone.
"""
M = [
    ('N01', 'make_target: no minimal separation', "if K == 1 or np.diff(c).min() >= MIN_SEP:", "if True:"),
    ('N02', 'make_target: amplitudes down to zero', "a = rng.uniform(0.3, 1.0, K) * rng.choice([-1.0, 1.0], K)", "a = rng.uniform(0.0, 1.0, K) * rng.choice([-1.0, 1.0], K)"),
    ('N03', 'make_target: the largest amplitude not normalised', "a[int(np.argmax(np.abs(a)))] = np.sign(a[int(np.argmax(np.abs(a)))])", "pass"),
    ('N04', 'make_target: all amplitudes positive', "* rng.choice([-1.0, 1.0], K)\n    a[int", "\n    a[int"),
    ('N05', '_score: an empty result counts as matched', "return [Khat, 0, np.nan, np.nan]", "return [Khat, Khat, np.nan, np.nan]"),
    ('N06', '_trial: the noise is not scaled by the largest amplitude', "s = NOISE[ni] * np.abs(at).max()\n    y = f + s", "s = NOISE[ni]\n    y = f + s"),
    ('N07', '_trial: the noise draw does not depend on the draw index', "np.random.default_rng(seed * 1000 + draw + 1)", "np.random.default_rng(seed * 1000 + 1)"),
    ('J01', 'run: target seeds', "jobs = [(K, ni, 1000 * K + t, d) for K in KS", "jobs = [(K, ni, 100 * K + t, d) for K in KS"),
    ('W01', 'wilson: z = 1', "def wilson(k, n, z=1.96):", "def wilson(k, n, z=1.0):"),
    ('W02', 'wilson: centred on the proportion', "mid = (p + z * z / (2 * n)) / den", "mid = p"),
    ('W03', 'wilson: half width without the second term', "half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den", "half = z * np.sqrt(p * (1 - p) / n) / den"),
    ('R01', 'report: exact against K + 1', "exact = res[:, :, 0] == Kt[:, None]", "exact = res[:, :, 0] == Kt[:, None] + 1"),
    ('R02', 'report: bias with the wrong sign', "d = res[:, j, 0] - Kt", "d = Kt - res[:, j, 0]"),
    ('R03', 'report: over- and under-estimates swapped', "'%.0f' % (100 * (d > 0).mean()), '%.0f' % (100 * (d < 0).mean())", "'%.0f' % (100 * (d < 0).mean()), '%.0f' % (100 * (d > 0).mean())"),
    ('R04', 'report: failed trials kept', "res, meta, secs = res[~bad], meta[~bad], secs[~bad]", "pass"),
    ('R05', 'report: the table by K groups by the noise index', "exact[Kt == k, j].mean()", "exact[ni == k, j].mean()"),
    ('R06', 'report: median cost replaced by the mean', "float(np.median(secs[:, i]))", "float(np.mean(secs[:, i]))"),
    ('R07', 'report: the all-noise column uses the first noise level', "k = int(exact[:, j].sum())\n        lo, hi = wilson(k, len(res))", "k = int(exact[ni == 0, j].sum())\n        lo, hi = wilson(k, len(res))"),
    ('S01', 'order_strata: identical and different swapped', "out[(Kn, seed)] = 0 if a == b else (1 if sorted(a) == sorted(b) else 2)", "out[(Kn, seed)] = 2 if a == b else (1 if sorted(a) == sorted(b) else 0)"),
    ('S02', 'order_strata: only the first entry compared', "a = [int(j) for j in B.iss_path(K, f, grid, max_entries=Kn)['order'][:Kn]]", "a = [int(j) for j in B.iss_path(K, f, grid, max_entries=Kn)['order'][:1]]"),
    ('S03', 'order_strata: a different set counts as the same kernels in another order', "(1 if sorted(a) == sorted(b) else 2)", "(1 if len(a) == len(b) else 2)"),
    ('S04', 'report: the strata table pairs the sliding ISS with itself', "('discrepancy', 'iss|disc', 'omp|disc')", "('discrepancy', 'iss|disc', 'iss|disc')"),
    ('S05', 'report: the strata count trials as targets', "row = [STRATA[g], str(len({(int(k), int(sd)) for k, sd in zip(Kt[m], meta[m, 2])})), str(int(m.sum()))]", "row = [STRATA[g], str(int(m.sum())), str(int(m.sum()))]"),
    ('N60', 'cluster diffs: the trials are resampled, not the targets', "ids = np.unique(targets)\n    groups = [np.flatnonzero(targets == t) for t in ids]", "ids = np.arange(len(targets))\n    groups = [np.flatnonzero(ids == t) for t in ids]"),
    ('N61', 'cluster diffs: the difference has the wrong sign', "d = np.array([exact[g, base].mean() - exact[g, jo].mean() for g in groups])", "d = np.array([exact[g, jo].mean() - exact[g, base].mean() for g in groups])"),
    ('N62', 'cluster diffs: a 90 % interval', "float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))", "float(np.percentile(boot, 5.0)), float(np.percentile(boot, 95.0))"),
    ('N63', 'cluster diffs: another seed', "np.random.default_rng(seed).integers(0, len(ids), size=(B, len(ids)))", "np.random.default_rng(seed + 1).integers(0, len(ids), size=(B, len(ids)))"),
    ('N64', 'report: the base column is compared with itself', "oth = [j for j in range(len(cols)) if j != jb]", "oth = [j for j in range(len(cols))]"),
    ('N65', 'report: the base is the discrepancy rule', "jb = cols.index('iss|bic')", "jb = cols.index('iss|disc')"),
    ('N66', 'target ids: the seed is ignored', "return np.asarray(Kt).astype(int) * 10 ** 6 + np.asarray(seeds).astype(int)", "return np.asarray(Kt).astype(int) * 10 ** 6"),
]
