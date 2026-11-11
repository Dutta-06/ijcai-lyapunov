"""Single-point mutants of experiments/exp_astro_m2.py (the selection of the frozen configuration and of the tuned baselines, the statistics, the verdicts, the freeze record and the
write-once test part), for tools/mutants.py; the referee's round IS1a mutants E1-E12, E14, E15 and ours for the freeze, parts D, DB and T:

    venv/bin/python tools/mutants.py --exp exp_astro_m2 --tests tests/test_astro_selection.py --raw results/smlm/calibration.md --md results/smlm/calibration.md \\
        --out-env MUT_OUT --spec tools/mutants/wp5_astro_m2_own.py --work <scratch dir> [--slot a --i0 0 --step 2]

(only `tests` means anything for these mutants: the fast tests of test_astro_selection.py). E13 (the clean-up threshold of `run_path`) needs the slow synthetic-path test of
tests/test_astro_common.py and is checked there. A tuple is (id, what, old, new[, k]): `old` replaced by `new` at its k-th occurrence in experiments/exp_astro_m2.py.
"""
M = [
    ('E1', 'choose: selects the smallest F', "if best is None or (round(f, 6), k) > (round(best[0], 6), best[3]):", "if best is None or (round(f, 6), k) < (round(best[0], 6), best[3]):"),
    ('E2', 'choose: ties go to the smaller kappa', "if best is None or (round(f, 6), k) > (round(best[0], 6), best[3]):", "if best is None or (round(f, 6), -k) > (round(best[0], 6), -best[3]):"),
    ('E3', 'f_measure: purity denominator = true22', "pur = tot[3] / max(tot[2], 1)\n    return 2 * comp", "pur = tot[3] / max(tot[0], 1)\n    return 2 * comp"),
    ('E4', 'choose: only the first development patch enters F', "tot = sum(raw['D|x%d_y%d|hp%g|rex%g|k%g|tot' % (p[0], p[1], hp, rex, k)] for p in DEV)\n                f = f_measure(tot)[0]\n                if best",
     "tot = sum(raw['D|x%d_y%d|hp%g|rex%g|k%g|tot' % (p[0], p[1], hp, rex, k)] for p in DEV[:1])\n                f = f_measure(tot)[0]\n                if best"),
    ('E5', 'choose_baselines: photutils picks the worst variant', "out['photutils'] = max(tags, key=lambda t: (round(_dev_F(raw, t), 6), t))", "out['photutils'] = min(tags, key=lambda t: (round(_dev_F(raw, t), 6), t))"),
    ('E5b', 'choose_baselines: crowdsource picks the worst filter', "out['crowdsource'] = max(names, key=lambda n: (round(_dev_F(raw, n), 6), n))", "out['crowdsource'] = min(names, key=lambda n: (round(_dev_F(raw, n), 6), n))"),
    ('E6', 'boot_diff: interval of b - a beside the point estimate a - b', "        d = aa[i] - bb[i]\n        out.append((a0[i] - b0[i], np.percentile(d, [2.5, 97.5])))", "        d = bb[i] - aa[i]\n        out.append((a0[i] - b0[i], np.percentile(d, [2.5, 97.5])))"),
    ('E7', 'verdict: win if the upper bound is positive', "return 'win' if lohi[0] > 0 else ('loss' if lohi[1] < 0 else 'tie')", "return 'win' if lohi[1] > 0 else ('loss' if lohi[0] < 0 else 'tie')"),
    ('E8', 'cs_dets: the crowdsource S/N filter inverted', "m = (b[key + 'flux'] > 0) & (b[key + 'flux'] >= snr * b[key + 'dflux'])", "m = (b[key + 'flux'] > 0) & (b[key + 'flux'] <= snr * b[key + 'dflux'])"),
    ('E9', 'boot_idx: the last block is never resampled', "return np.random.default_rng(seed).integers(0, nb, size=(B, nb))", "return np.random.default_rng(seed).integers(0, nb - 1, size=(B, nb))"),
    ('E10', 'boot_cpf: 90 % interval labelled 95 %', "return [(v0, np.percentile(v, [2.5, 97.5])) for v0, v in ((c0, cs), (p0, ps), (f0, fs))]", "return [(v0, np.percentile(v, [5.0, 95.0])) for v0, v in ((c0, cs), (p0, ps), (f0, fs))]"),
    ('E11', 'cpf: purity denominator = true22', "p = tot[..., 3] / np.maximum(tot[..., 2], 1)", "p = tot[..., 3] / np.maximum(tot[..., 0], 1)"),
    ('E12', 'part_T: no tuned baselines run', "    tuned = choose_baselines(raw)\n    todo = []", "    tuned = {}\n    todo = []"),
    ('E14', 'part_T: stored default catalogue at the tuned kappa', "dk = kap if label == 'ours' else DEFAULT[2]", "dk = kap if label == 'ours' else kap"),
    ('E15', 'REXS grid enlarged without the D part', "REXS = (0.8, 1.2)", "REXS = (0.8, 1.2, 0.0)"),
    ('F1', 'freeze: the hashes of the input files are constants', "files={f_: _sha256_file(os.path.join(REPO, f_)) for f_ in FREEZE_FILES}", "files={f_: 'x' for f_ in FREEZE_FILES}"),
    ('F2', 'freeze: check_freeze ignores differences', "    if diff:\n        raise RuntimeError('the freeze record differs", "    if False:\n        raise RuntimeError('the freeze record differs"),
    ('F3', 'freeze: the development results are not hashed', "frozen_code=frozen_code_hash(), dev_keys=dev_keys_hash(raw))", "frozen_code=frozen_code_hash(), dev_keys='x')"),
    ('F4', 'freeze: the constants are not hashed', "consts = repr((TESTS, DEV, DEFAULT, KAPPAS, HPS, REXS, LAT_H, MARGIN, FIT_RADIUS, MAX_ENTRIES, PHOT_TAGS, CS_FILTERS))", "consts = ''"),
    ('F5', 'part_T: a patch with a result is run again', "        if key_fmt % pxy + '|cfg' in raw:\n            print('T patch %s %s: a result exists, not run again' % (ds, pxy), flush=True)\n        else:\n            todo.append(pxy)", "        todo.append(pxy)"),
    ('F6', 'part_T: no check of the freeze', "    raw = load_raw(out_path)\n    check_freeze(raw)\n    hp, rex, kap, f = choose(raw)\n    tuned = choose_baselines(raw)", "    raw = load_raw(out_path)\n    hp, rex, kap, f = choose(raw)\n    tuned = choose_baselines(raw)"),
    ('F7', 'part_T: the as-shipped DAOPHOT row stores the shifted catalogue', "arrays.update(_pack('daophot (as shipped)', bs['daophot_asis']))", "arrays.update(_pack('daophot (as shipped)', bs['daophot']))"),
    ('F8', 'part_DB: keys under D', "    save_part(out_path, 'DB', arrays)", "    save_part(out_path, 'D', arrays)"),
    ('F9', 'part_D: the cap flag never set', "int(n >= MAX_ENTRIES)])\n        print('D patch", "0])\n        print('D patch"),
    ('F10', 'boot_idx_strat: offsets dropped', "out.append(off + rng.integers(0, n, size=(B, n)))", "out.append(rng.integers(0, n, size=(B, n)))"),
    ('F11', 'baseline_scores: the unshifted catalogue is the shifted one', "    for name in ('daophot', 'daophot_asis'):", "    for name in ('daophot', 'daophot'):"),
    ('F12', 'dev_keys_hash: only the keys, not the values', "        h.update(np.ascontiguousarray(raw[k]).tobytes())", "        h.update(b'')"),
]
