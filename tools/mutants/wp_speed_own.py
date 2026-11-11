"""Single-point mutants of experiments/exp_smlm_speed.py, for tools/mutants.py:

    venv/bin/python tools/mutants.py --exp exp_smlm_speed --tests tests/test_smlm_speed.py --raw results/smlm/speed_raw.npz \\
        --md results/smlm/speed.md --out-env SPEED_OUT --spec tools/mutants/wp_speed_own.py --work <scratch dir>

A tuple is (id, what, old, new[, k]): `old` replaced by `new` at its k-th occurrence in experiments/exp_smlm_speed.py. The kill criterion is the tests file alone.
"""
M = [
    ('E01', 'frames: more frames than the range holds are asked for', "n = min(n, b - a)", "n = n"),
    ('E02', 'frames: the last index is one past the range', "np.linspace(0, b - a - 1, n)", "np.linspace(0, b - a, n)"),
    ('E03', 'frames: the range start is dropped', "return a + np.round(", "return np.round("),
    ('C01', 'compare: only frames with fewer full detections are skipped', "        if len(a) != len(b):\n            continue", "        if len(a) < len(b):\n            continue"),
    ('C02', 'compare: the full detections are sorted by y then x', "oa, ob = np.lexsort((a[:, 1], a[:, 0])), np.lexsort((b[:, 1], b[:, 0]))", "oa, ob = np.lexsort((a[:, 0], a[:, 1])), np.lexsort((b[:, 1], b[:, 0]))"),
    ('C03', 'compare: the smallest position difference of a frame', "np.abs(a[oa][:, :2] - b[ob][:, :2]).max()", "np.abs(a[oa][:, :2] - b[ob][:, :2]).min()"),
    ('C04', 'compare: the smallest position difference over frames', "diff_pos = max(diff_pos, ", "diff_pos = min(diff_pos, "),
    ('C05', 'compare: the smallest photon difference over frames', "diff_ph = max(diff_ph, ", "diff_ph = min(diff_ph, "),
    ('C06', 'compare: the photon difference relative to the other setting', "np.maximum(np.abs(a[oa][:, 2]), 1e-9)", "np.maximum(np.abs(b[ob][:, 2]), 1e-9)"),
    ('C07', 'compare: a frame counts twice', "same += 1", "same += 2"),
    ('C08', 'compare: a single detection is not compared', "        if len(a):\n", "        if len(a) > 1:\n"),
    ('C09', 'compare: the photon difference is signed', "(np.abs(a[oa][:, 2] - b[ob][:, 2]) / ", "((a[oa][:, 2] - b[ob][:, 2]) / "),
    ('T01', 'timed: wall-clock instead of process CPU time', "t0 = time.process_time()", "t0 = time.time()"),
    ('T02', 'timed: two worker processes', "opts, workers=1)", "opts, workers=2)"),
    ('T03', 'timed: another significance level', "rules=[('snr', ISS.DEFAULT['kappa'])]", "rules=[('snr', 3.0)]"),
    ('T04', 'timed: all rules returned, not the first', "return dets[0], time.process_time() - t0", "return dets, time.process_time() - t0"),
    ('T05', 'timed: another configuration', "ISS.run_frames(frames, ISS.DEFAULT['cfg'], ", "ISS.run_frames(frames, 'hp3', "),
    ('S01', 'settings: full tolerance 1e-6', "dict(tol=1e-8, max_nfev=200)", "dict(tol=1e-6, max_nfev=200)"),
    ('S02', 'settings: full evaluations 100', "dict(tol=1e-8, max_nfev=200)", "dict(tol=1e-8, max_nfev=100)"),
    ('S03', 'settings: fast tolerance 1e-3', "dict(tol=1e-4, max_nfev=25)", "dict(tol=1e-3, max_nfev=25)"),
    ('S04', 'settings: fast evaluations 30', "dict(tol=1e-4, max_nfev=25)", "dict(tol=1e-4, max_nfev=30)"),
    ('U01', 'run: no alternation of the settings', "for label, setting in (SETTINGS if rep == 0 else SETTINGS[::-1]):", "for label, setting in SETTINGS:"),
    ('U02', 'run: the larger of the two repeats is kept', "best[label] = min(best.get(label, np.inf), cpu)", "best[label] = max(best.get(label, np.inf), cpu)"),
    ('U03', 'run: the times stored in the other order', "np.array([best['full'], best['fast']])", "np.array([best['fast'], best['full']])"),
    ('U04', 'run: the fast detections compared with themselves', "compare(dets['full'], dets['fast'])", "compare(dets['fast'], dets['fast'])"),
    ('U05', 'run: one more frame than asked', "idx = evenly_spaced(a, b, n_frames)", "idx = evenly_spaced(a, b, n_frames + 1)"),
    ('U06', 'run: development frames instead of test frames', "ISS.dev_test_frames(name, 'test')", "ISS.dev_test_frames(name, 'dev')"),
    ('U07', 'run: detection counts stored in the other order', "np.array([sum(len(d) for d in dets['full']), sum(len(d) for d in dets['fast'])])", "np.array([sum(len(d) for d in dets['fast']), sum(len(d) for d in dets['full'])])"),
    ('U08', 'run: a single repeat', "for rep in range(2):", "for rep in range(1):"),
    ('R01', 'report: full time not per frame', "'%.3f' % (cpu[0] / n)", "'%.3f' % cpu[0]"),
    ('R02', 'report: fast time not per frame', "'%.3f' % (cpu[1] / n)", "'%.3f' % cpu[1]"),
    ('R03', 'report: speed-up inverted', "'%.2f' % (cpu[0] / cpu[1])", "'%.2f' % (cpu[1] / cpu[0])"),
    ('R04', 'report: equal-count frames swapped with the frames', "'%d / %d' % (cmp_[0], cmp_[1])", "'%d / %d' % (cmp_[1], cmp_[0])"),
    ('R05', 'report: detection counts swapped', "'%d / %d' % (nd[0], nd[1])", "'%d / %d' % (nd[1], nd[0])"),
    ('R06', 'report: position column prints the photon difference', "'%.4f' % cmp_[2]", "'%.4f' % cmp_[3]"),
    ('R07', 'report: photon column prints the position difference', "'%.2e' % cmp_[3]", "'%.2e' % cmp_[2]"),
]
