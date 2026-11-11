"""Single-point mutants of experiments/exp_astro_repeat.py (the held-out repeat-frame test: the per-patch bookkeeping, the freeze, write-once running, the report), for tools/mutants.py:

    venv/bin/python tools/mutants.py --exp exp_astro_repeat --tests tests/test_exp_astro_repeat.py --raw results/astro/repeat_raw.npz --md results/astro/repeat.md --out-env REPEAT_OUT \\
        --extra results/astro/m2_raw.npz --extra results/astro/freeze.json --extra results/astro/freeze_repeat.json --extra results/astro/repeat_plan.md \\
        --extra results/astro/baseline_catalogs_repeat.npz --extra data/astro/checksums_repeat.json --spec tools/mutants/wp_astrorepeat_own.py --work <scratch dir>

A tuple is (id, what, old, new[, k]): `old` replaced by `new` at its k-th occurrence (0-based) in experiments/exp_astro_repeat.py. The kill criterion is the tests file alone.
"""
M = [
    ('A01', 'patch_arrays: the kernel flag is inverted', "arrays['kernel|same'] = np.array([int(fit['kind'] == 'isotropic')])", "arrays['kernel|same'] = np.array([int(fit['kind'] != 'isotropic')])"),
    ('A02', 'patch_arrays: the seeing scale is inverted', "s = fwhm / fw_ref", "s = fw_ref / fwhm"),
    ('A03', 'patch_arrays: the scaled run does not scale sigma_b', "dets3, secs3, n3, info3 = PS.run_path_psf(patch, h * s, r * s, fit)", "dets3, secs3, n3, info3 = PS.run_path_psf(patch, h, r * s, fit)"),
    ('A04', 'patch_arrays: the scaled run does not scale r_ex', "dets3, secs3, n3, info3 = PS.run_path_psf(patch, h * s, r * s, fit)", "dets3, secs3, n3, info3 = PS.run_path_psf(patch, h * s, r, fit)"),
    ('A05', 'patch_arrays: the isotropic run is not stored under the chosen-kernel name', "            _alias(arrays, 'iso|' + label, 'adq|' + label)", "            pass"),
    ('A06', 'patch_arrays: the a priori default runs at the tuned configuration', "('default', (X.DEFAULT[0], X.DEFAULT[1]), X.DEFAULT[2])):\n        dets, secs, n, info = X.run_path(patch, h, r)", "('default', (hp, rex), X.DEFAULT[2])):\n        dets, secs, n, info = X.run_path(patch, h, r)"),
    ('A07', 'patch_arrays: the default is reported at the tuned significance', "('default', (X.DEFAULT[0], X.DEFAULT[1]), X.DEFAULT[2])):\n        dets, secs, n, info = X.run_path(patch, h, r)", "('default', (X.DEFAULT[0], X.DEFAULT[1]), kap)):\n        dets, secs, n, info = X.run_path(patch, h, r)"),
    ('A08', 'patch_arrays: the tuned baseline is stored as the default', "arrays.update(X._pack('%s (tuned)' % fam, bs[nm]))", "arrays.update(X._pack('%s (default)' % fam, bs[nm]))"),
    ('A09', 'patch_arrays: the chosen kernel\'s residual scale is read at significance 0', "arrays['adq|%s|rsd' % label] = np.array([info2['rsd'][dk]])", "arrays['adq|%s|rsd' % label] = np.array([info2['rsd'][0.0]])"),
    ('A10', 'patch_arrays: the seeing pair is stored in the wrong order', "arrays['fwhm'] = np.array([fwhm, fw_ref])", "arrays['fwhm'] = np.array([fw_ref, fwhm])"),
    ('A11', 'patch_arrays: the isotropic residual scale has no high-pass', "P_iso = IM.PixelPSF(*patch['dn'].shape, sig, parts=parts, hp=h or None)", "P_iso = IM.PixelPSF(*patch['dn'].shape, sig, parts=parts, hp=None)"),
    ('A12', 'patch_arrays: the adequate-kernel run is made for every patch', "        if fit['kind'] == 'isotropic':\n            _alias(arrays, 'iso|' + label, 'adq|' + label)\n        else:", "        if False:\n            _alias(arrays, 'iso|' + label, 'adq|' + label)\n        else:"),
    ('F01', 'check_freeze: any difference is accepted', "diff = [k for k in now if now[k] != old.get(k)]", "diff = [k for k in now if now[k] == old.get(k) and False]"),
    ('F02', 'check_freeze: a missing record is accepted', "    if not os.path.exists(path):\n        raise RuntimeError('no freeze record", "    if False:\n        raise RuntimeError('no freeze record"),
    ('F03', 'frozen_code_hash: the patch list is not hashed', "consts = repr((PATCHES, X.DEFAULT,", "consts = repr((X.DEFAULT,"),
    ('F04', 'freeze_record: the truth files are not hashed', "        data[os.path.join(R.REPEAT[name]['root'], R.REPEAT[name]['truth'])] = None", "        pass"),
    ('F05', 'freeze_record: the reference seeing is not recorded', "fwhm_ref=round(ref_fwhm(), 6)", "fwhm_ref=0.0"),
    ('F06', 'load_frozen: another default is accepted', "    if tuple(rec['default']) != tuple(X.DEFAULT):", "    if False:"),
    ('F07', 'load_frozen: the tuned baselines are not read', "tuned=dict(rec['tuned_baselines']))", "tuned={})"),
    ('F08', 'frozen_code_hash: the kernel rule constants are not hashed', "G.RMS_GOAL, G.K_MAX, X.DEV, LABELS, FILES))", "X.DEV, LABELS, FILES))"),
    ('P01', 'part_run: finished patches are run again', "if 'R|%s|x%d_y%d|cfg' % (n, p[0], p[1]) not in cur]", "if 'R|%s|x%d_y%d|cfg' % (n, p[0], p[1]) in cur]"),
    ('P02', 'part_run: the freeze is not checked', "    rec = check_freeze(freeze_path)\n    fr = load_frozen()", "    rec = dict(fwhm_ref=2.15)\n    fr = load_frozen()"),
    ('P03', 'part_run: the stored key swaps the origin', "            for name, pxy, arrays in pool.imap_unordered(patch_arrays, jobs, chunksize=1):\n                X.save_part(out_path, 'R|%s|x%d_y%d' % (name, pxy[0], pxy[1]), arrays)", "            for name, pxy, arrays in pool.imap_unordered(patch_arrays, jobs, chunksize=1):\n                X.save_part(out_path, 'R|%s|x%d_y%d' % (name, pxy[1], pxy[0]), arrays)"),
    ('P04', 'part_run: the reference seeing is not passed to the jobs', "jobs = [(n, p, fr['hp'], fr['rex'], fr['kap'], fr['tuned'], rec['fwhm_ref'], base) for n, p in todo]", "jobs = [(n, p, fr['hp'], fr['rex'], fr['kap'], fr['tuned'], 2.0, base) for n, p in todo]"),
    ('T01', 'purity_all: the purity at mag < 22', "return tot[5] / max(tot[6], 1), int(tot[6] - tot[5])", "return tot[3] / max(tot[2], 1), int(tot[6] - tot[5])"),
    ('T02', 'purity_all: the unmatched count is the detections', "return tot[5] / max(tot[6], 1), int(tot[6] - tot[5])", "return tot[5] / max(tot[6], 1), int(tot[6])"),
    ('T03', 'p1_text: the patches where the kernels are the same', "gen = [np_ for np_ in plist if raw[pre_of(*np_) + 'kernel|same'][0] == 0]", "gen = [np_ for np_ in plist if raw[pre_of(*np_) + 'kernel|same'][0] == 1]"),
    ('T04', 'p1_text: the verdict wording is swapped', "{'win': 'confirmed on held-out frames', 'tie': 'not confirmed', 'loss': 'contradicted'}", "{'win': 'contradicted', 'tie': 'not confirmed', 'loss': 'confirmed on held-out frames'}"),
    ('T05', 'pair_table: the difference is taken the other way', "r, df = diff_row(name, blk[a], blk[b], idx)", "r, df = diff_row(name, blk[b], blk[a], idx)"),
    ('T06', 'report: a frame pools every patch', "            fl = [np_ for np_ in plist if np_[0] == fname]", "            fl = plist"),
    ('T07', 'report: the ladder lists the other windows', "if tuple(R.REPEAT[n]['windows'][R.REPEAT[n]['patches'].index(tuple(p))]) == tuple(w))", "if tuple(R.REPEAT[n]['windows'][R.REPEAT[n]['patches'].index(tuple(p))]) != tuple(w))"),
    ('T08', 'report: the exploratory rows are never shown', "if all(pre_of(*np_) + 'scl|ours|k2|blk22' in raw for np_ in plist):", "if False:"),
    ('T09', 'report: the run records list the chosen kernel where it is the same', "                if kern == 'adq' and raw[pre + 'kernel|same'][0] == 1:", "                if kern == 'adq' and raw[pre + 'kernel|same'][0] == 0:"),
    ('T10', 'report: the cap flag is inverted', "'capped' if b_[2] else '-'", "'-' if b_[2] else 'capped'"),
    ('T11', 'report: the residual table says same for the wrong patches', "            same = raw[pre + 'kernel|same'][0] == 1\n            rows.append([name, '(%d, %d)' % p, '%.2f' % raw[pre + 'iso|ours|rsd'][0]", "            same = raw[pre + 'kernel|same'][0] == 0\n            rows.append([name, '(%d, %d)' % p, '%.2f' % raw[pre + 'iso|ours|rsd'][0]"),
    ('T12', 'method_table: the bootstrap is not stratified', "    idx = X.boot_idx_strat(sizes)\n    rows = []\n    for label, key in defs:", "    idx = X.boot_idx(sum(sizes))\n    rows = []\n    for label, key in defs:"),
    ('T13', 'report: the patch table reads the chosen kernel\'s rms twice', "'%.5f' % rms[0], kind,", "'%.5f' % rms[1], kind,"),
    ('T14', 'report: the truth count is read from the default run', "tr = raw[pre + 'iso|ours|k2|tot']", "tr = raw[pre + 'iso|default|k5|tot']"),
    ('T15', 'report: the record by patch uses the pooled resamples', "                    idx1 = X.boot_idx(raw[pre + a + '|blk22'].shape[0])", "                    idx1 = X.boot_idx(raw[pre + a + '|blk22'].shape[0], seed=1)"),
]
