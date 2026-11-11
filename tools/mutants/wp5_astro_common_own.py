"""Single-point mutants of experiments/astro_common.py (the scorer of the astronomy experiment), for tools/mutants.py; the referee's round IS1a mutants A1-A10 (A10 is an equivalent control):

    venv/bin/python tools/mutants.py --exp astro_common --tests tests/test_astro_selection.py --raw results/smlm/calibration.md --md results/smlm/calibration.md \\
        --out-env MUT_OUT --spec tools/mutants/wp5_astro_common_own.py --work <scratch dir> [--slot a --i0 0 --step 2]

(only `tests` means anything for these mutants). A tuple is (id, what, old, new[, k]): `old` replaced by `new` at its k-th occurrence in experiments/astro_common.py.
"""
M = [
    ('A1', 'match: position AND magnitude -> OR', "ok = (d <= dpos) & (dm <= dmag)", "ok = (d <= dpos) | (dm <= dmag)"),
    ('A2', 'score: detections in the edge zone are scored', "d_keep = ok & inside", "d_keep = ok"),
    ('A3', 'score: matched truth counted in the detection block', "np.add.at(blk[:, 1], block_of(t_xy[mt[sel_mt]], n), 1)", "np.add.at(blk[:, 1], block_of(det_xy[md[sel_mt]], n), 1)"),
    ('A4', 'score: matched detections selected by the truth magnitude below 22', "dmag_ < 22.0, det_mag[md] < 22.0", "dmag_ < 22.0, t['mag'][mt] < 22.0"),
    ('A5', 'match: any maximum matching (no minimum distance)', "r, c = linear_sum_assignment(np.where(ok, d, big))", "r, c = linear_sum_assignment(np.where(ok, 0.0, big))"),
    ('A6', 'truth_for_patch: upper y margin 0.5 px too small', "(y <= n - 0.5 - EDGE)", "(y <= n - EDGE)"),
    ('A7', 'score: non-positive-flux detections not dropped', "    ok = np.isfinite(det_mag)\n    inside", "    ok = np.ones(len(det_mag), bool)\n    inside"),
    ('A8', 'block_of: last block column merged into the one before', "ix = np.clip(np.floor((xy[:, 0] + 0.5) / cell).astype(int), 0, nb - 1)", "ix = np.clip(np.floor((xy[:, 0] + 0.5) / cell).astype(int), 0, nb - 2)"),
    ('A9', 'score: truth depth 22.5 instead of 22 for the totals', "sel22t, sel22d, sel22m = tmag < 22.0, dmag_ < 22.0,", "sel22t, sel22d, sel22m = tmag < 22.5, dmag_ < 22.0,"),
    ('A10', 'control (equivalent): must survive', "return np.zeros(0, int), np.zeros(0, int)\n    d = np.linalg.norm", "return np.zeros(0, int), np.zeros(0, int)\n    tr_xy = tr_xy + 0.0; d = np.linalg.norm"),
]
