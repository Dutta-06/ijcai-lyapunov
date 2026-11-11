"""Single-point mutants of experiments/exp_astro_hp2.py (round 2 of the baselines given the ISS's front end), for tools/mutants.py:

    venv/bin/python tools/mutants.py --exp exp_astro_hp2 --tests tests/test_exp_astro_hp2.py --raw results/astro/baselines_hp2_raw.npz --md results/astro/baselines_hp2.md --out-env HP2_OUT \\
        --extra results/astro/m2_raw.npz --extra results/astro/psf_raw.npz --extra results/astro/repeat_raw.npz --extra results/astro/baselines_hp_raw.npz --extra results/astro/baselines_hp_plan.md \\
        --extra results/astro/baselines_hp2_plan.md --extra tests/test_astro_report.py --extra tests/test_exp_astro_hp.py --extra tests/test_exp_astro_baselines_hp.py --extra tests/test_exp_astro_baselines_hp2.py \\
        --spec tools/mutants/wp_astrohp2_own.py --work <scratch dir> [--slot a --i0 0 --step 2]

A tuple is (id, what, old, new[, k]): `old` replaced by `new` at its k-th occurrence in experiments/exp_astro_hp2.py. The kill criterion is the tests file alone. The slot's report is regenerated with the
round-1 raw copied by --extra to results/astro/ (the default of HP2_RAW1).

Survivor on the final run: D10 (the guard `set(sel) != set(G2.FAMILIES2)` in part_dev2 removed). Equivalent on every reachable input: the round-1 selected variant of each family is always a candidate and is scored on both
development patches before the selection (its totals must reproduce the round-1 ones, else part_dev2 has already raised), so `choose_hp2` returns both families; with the guard removed a missing family would still stop
part_dev2 with a KeyError at `sel['photutils']` before anything is stored, i.e. only the exception type would differ.
"""
M = [
    # ---- the selection rule
    ('S01', 'choose_hp2: F compared after rounding to 3 decimals', "max(elig, key=lambda t: (round(E1.dev_F(raw, t), 6), t))", "max(elig, key=lambda t: (round(E1.dev_F(raw, t), 3), t))"),
    ('S02', 'choose_hp2: F not rounded', "max(elig, key=lambda t: (round(E1.dev_F(raw, t), 6), t))", "max(elig, key=lambda t: (E1.dev_F(raw, t), t))"),
    ('S03', 'choose_hp2: a tie goes to the smaller tag', "max(elig, key=lambda t: (round(E1.dev_F(raw, t), 6), t))", "max(elig, key=lambda t: (round(E1.dev_F(raw, t), 6), tuple(-ord(c) for c in t)))"),
    ('S04', 'choose_hp2: the smallest F', "out[fam] = max(elig,", "out[fam] = min(elig,"),
    ('S05', 'choose_hp2: a candidate with a score on one patch is eligible', "elig = [t for t in tags if all('DEV|x%d_y%d|%s|tot' % (p[0], p[1], t) in raw for p in X.DEV)]", "elig = [t for t in tags if any('DEV|x%d_y%d|%s|tot' % (p[0], p[1], t) in raw for p in X.DEV)]"),
    ('S06', 'candidates: the round-1 selection is not a candidate', "return {fam: G2.FAMILIES2[fam] + [sel1[fam]] for fam in G2.FAMILIES2}", "return {fam: G2.FAMILIES2[fam] for fam in G2.FAMILIES2}"),
    ('S07', 'candidates: the round-1 selection of the other family', "return {fam: G2.FAMILIES2[fam] + [sel1[fam]] for fam in G2.FAMILIES2}", "return {fam: G2.FAMILIES2[fam] + [sel1['photutils']] for fam in G2.FAMILIES2}"),
    ('S08', 'candidates: the round-1 selection first', "return {fam: G2.FAMILIES2[fam] + [sel1[fam]] for fam in G2.FAMILIES2}", "return {fam: [sel1[fam]] + G2.FAMILIES2[fam] for fam in G2.FAMILIES2}"),
    # ---- the freeze
    ('V01', 'verify_selection2: one missing key is enough to pass', "if any(k not in raw for k in keys):", "if all(k not in raw for k in keys):"),
    ('V02', 'verify_selection2: the round-1 selection is not verified', "    sel1 = round1_selection(raw1)\n    if {f:", "    sel1 = {f: str(raw['SEL|r1_' + f][0]) for f in G2.FAMILIES2}\n    if {f:"),
    ('V03', 'verify_selection2: the round-1 selection used is not compared', "if {f: str(raw['SEL|r1_' + f][0]) for f in G2.FAMILIES2} != sel1:", "if False:"),
    ('V04', 'verify_selection2: the stored selection is not compared with the rule', "if stored != choose_hp2(raw, sel1):", "if False:"),
    ('V05', 'verify_selection2: the development hash is not checked', "if str(raw['SEL|dev_hash'][0]) != E1.dev_hash(raw):", "if False:"),
    ('V06', 'verify_selection2: the plan hash is not checked', "if str(raw['SEL|plan_sha256'][0]) != E1.sha256_file(plan):", "if False:"),
    ('V07', 'verify_selection2: the default plan ignored', "plan = plan or PLAN2", "plan = PLAN2"),
    ('V08', 'verify_selection2: the raw1 argument ignored', "    sel1 = round1_selection(raw1)\n    if {f:", "    sel1 = round1_selection()\n    if {f:"),
    ('V09', 'round1_selection: the raw1 argument ignored', "return E1.verify_selection(raw1 if raw1 is not None else X.load_raw(RAW1))", "return E1.verify_selection(X.load_raw(RAW1))"),
    ('V10', 'round1_selection: not verified', "return E1.verify_selection(raw1 if raw1 is not None else X.load_raw(RAW1))", "r = raw1 if raw1 is not None else X.load_raw(RAW1)\n    return {f: str(r['SEL|' + f][0]) for f in G2.FAMILIES2}"),
    # ---- part_dev2
    ('D01', 'part_dev2: the first development patch only', "    for p in X.DEV:\n        t = A.truth_for_patch(truth, p[0], p[1])\n        for b, tags in", "    for p in X.DEV[:1]:\n        t = A.truth_for_patch(truth, p[0], p[1])\n        for b, tags in"),
    ('D02', 'part_dev2: the seconds not stored', "                    if k + 'secs' in b.files:\n                        arrays[k + 'secs'] = np.asarray(b[k + 'secs'])", "                    pass"),
    ('D03', 'part_dev2: a second selection allowed', "if 'SEL|photutils' in cur:", "if False:"),
    ('D04', 'part_dev2: the round-1 variants need not reproduce their totals', "if not np.array_equal(arrays[k], raw1['DEV|' + k]):", "if False:"),
    ('D05', 'part_dev2: the plan hash of nothing', "plan_sha256=np.array([E1.sha256_file(PLAN2)])", "plan_sha256=np.array(['0' * 64])"),
    ('D06', 'part_dev2: the round-1 selection stored exchanged', "r1_photutils=np.array([sel1['photutils']]),\n                                       r1_crowdsource=np.array([sel1['crowdsource']])", "r1_photutils=np.array([sel1['crowdsource']]),\n                                       r1_crowdsource=np.array([sel1['photutils']])"),
    ('D07', 'part_dev2: the truth of the wrong patch', "t = A.truth_for_patch(truth, p[0], p[1])\n        for b, tags in", "t = A.truth_for_patch(truth, p[1], p[0])\n        for b, tags in"),
    ('D08', 'part_dev2: the round-1 crowdsource selection is not scored', "(b1, [sel1['photutils'], sel1['crowdsource']])", "(b1, [sel1['photutils']])"),
    ('D09', 'part_dev2: the round-2 crowdsource variants are not scored', "(b2, G2.FAMILIES2['photutils'] + G2.FAMILIES2['crowdsource'])", "(b2, G2.FAMILIES2['photutils'])"),
    ('D10', 'part_dev2: a family without a candidate tolerated', "if set(sel) != set(G2.FAMILIES2):", "if False:"),
    ('D11', 'part_dev2: the round-1 variants read from the round-2 file', "(b1, [sel1['photutils'], sel1['crowdsource']])", "(b2, [sel1['photutils'], sel1['crowdsource']])"),
    ('D12', 'part_dev2: the development hash of nothing', "dev_hash=np.array([E1.dev_hash(raw)])", "dev_hash=np.array(['0' * 64])"),
    # ---- part_test2
    ('T01', 'part_test2: the selection not verified', "sel = verify_selection2(raw, raw1)", "sel = {f: str(raw['SEL|' + f][0]) for f in G2.FAMILIES2}"),
    ('T02', 'part_test2: the repeat frames without the flag', "if name in G.REPEAT_SETS and not (flag and os.path.exists(flag)):", "if False:"),
    ('T03', 'part_test2: the flag need not exist', "not (flag and os.path.exists(flag))", "not flag"),
    ('T04', 'part_test2: a patch with a result scored again', "if pre + '|phu2|tot' in raw:", "if False:"),
    ('T05', 'part_test2: the origin of the key swapped', "pre = 'T|%s|x%d_y%d' % (name, p[0], p[1])\n            if pre", "pre = 'T|%s|x%d_y%d' % (name, p[1], p[0])\n            if pre"),
    ('T06', 'part_test2: a repeat frame scored against the M2 truth', "truth = R.truth_arrays(name) if name in G.REPEAT_SETS else X.truth_arrays(ds)", "truth = X.truth_arrays(ds)"),
    ('T07', 'part_test2: a second cluster scored against the M2 truth', "else X.truth_arrays(ds)", "else X.truth_arrays('m2')"),
    ('T08', 'part_test2: the tag not stored', "                arrays[row + '|name'] = np.array([tag])", "                pass"),
    ('T09', 'part_test2: the seconds of the other variant', "arrays[row + '|secs'] = np.asarray(b['x%d_y%d|%s|secs' % (p[0], p[1], tag)])", "arrays[row + '|secs'] = np.asarray(b['x%d_y%d|%s|secs' % (p[0], p[1], tags['phu2'])])"),
    ('T10', 'part_test2: the two rows exchanged', "tags = {'phu2': sel['photutils'], 'cs2': sel['crowdsource']}", "tags = {'phu2': sel['crowdsource'], 'cs2': sel['photutils']}"),
    ('T11', 'part_test2: an unknown set accepted', "if name not in G.TEST_SETS:", "if False:"),
    ('T12', 'part_test2: the truth of the wrong patch', "t = A.truth_for_patch(truth, p[0], p[1])\n            arrays = {}", "t = A.truth_for_patch(truth, p[1], p[0])\n            arrays = {}"),
    # ---- the rows of the report
    ('P01', 'patch_rows2: no stamp-kernel ISS row', "if r in ('iss', 'iss_adq', 'phu_old_t', 'cs_old_t')]", "if r in ('iss', 'phu_old_t', 'cs_old_t')]"),
    ('P02', 'patch_rows2: the round-1 rows from each other\'s keys', "for row, src in (('phu1', 'phu_sel'), ('cs1', 'cs_sel')):", "for row, src in (('phu1', 'cs_sel'), ('cs1', 'phu_sel')):"),
    ('P03', 'patch_rows2: the round-1 rows missing', "if raw1 is not None and pre1 + src + '|tot' in raw1:", "if False:"),
    ('P04', 'patch_rows2: the round-2 rows from the round-1 raw', "        if pre1 + row + '|tot' in raw:\n            rows.append((row, '%s [%s]' % (LABELS2[row], E1.shown(raw[pre1 + row + '|name'][0])), raw[pre1 + row + '|tot'], raw[pre1 + row + '|blk22'], raw[pre1 + row + '|offset']))",
     "        if pre1 + row + '|tot' in raw:\n            rows.append((row, '%s [%s]' % (LABELS2[row], E1.shown(raw[pre1 + row + '|name'][0])), raw[pre1 + row + '|tot'], raw[pre1 + row + '|blk22'][::-1], raw[pre1 + row + '|offset']))"),
    ('P05', 'patch_rows2: the tag printed with its pipes', "(row, '%s [%s]' % (LABELS2[row], E1.shown(raw[pre1 + row + '|name'][0])),", "(row, '%s [%s]' % (LABELS2[row], str(raw[pre1 + row + '|name'][0])),"),
    ('P06', 'patch_rows2: the round-1 tag printed with its pipes', "'%s [%s]' % (LABELS2[row], E1.shown(raw1[pre1 + src + '|name'][0]))", "'%s [%s]' % (LABELS2[row], str(raw1[pre1 + src + '|name'][0]))"),
    ('P07', 'group_data2: rows present on one patch are pooled', "rows = [r for r in ORDER2 if r in per and len(per[r]) == len(done)]", "rows = [r for r in ORDER2 if r in per]"),
    ('P08', 'group_data2: the strata are the pooled blocks', "sizes=[len(b) for l, t, b in per[rows[0]]] if rows else [])", "sizes=[sum(len(b) for l, t, b in per[rows[0]])] if rows else [])"),
    ('P09', 'group_data2: the totals of the first patch only', "tot={r: sum(t for l, t, b in per[r]) for r in rows}", "tot={r: per[r][0][1] for r in rows}"),
    ('P10', 'group_data2: a patch without round-2 rows is a patch', "if 'T|%s|x%d_y%d|phu2|tot' % (name, p[0], p[1]) not in raw:", "if False:"),
    ('P11', 'ORDER2: the round-2 rows first', "ORDER2 = ('iss', 'iss_adq', 'phu_old_t', 'cs_old_t', 'phu1', 'cs1') + NEW2_ROWS", "ORDER2 = NEW2_ROWS + ('iss', 'iss_adq', 'phu_old_t', 'cs_old_t', 'phu1', 'cs1')"),
    ('P12', 'diff_table2: second minus first', "dc, dp, df = X.boot_diff(g['blk'][a], g['blk'][b], idx)", "dc, dp, df = X.boot_diff(g['blk'][b], g['blk'][a], idx)"),
    ('P13', 'diff_table2: completeness and purity columns exchanged', "X.fmt_iv(dc[0], dc[1]), X.fmt_iv(dp[0], dp[1]), X.fmt_iv(df[0], df[1]), X.verdict(*df)]", "X.fmt_iv(dp[0], dp[1]), X.fmt_iv(dc[0], dc[1]), X.fmt_iv(df[0], df[1]), X.verdict(*df)]"),
    ('P14', 'diff_table2: the verdict of the completeness', "X.fmt_iv(df[0], df[1]), X.verdict(*df)]", "X.fmt_iv(df[0], df[1]), X.verdict(*dc)]"),
    ('P15', 'COMPARISONS2_PATCH: round 2 against the old default instead of the round-1 selection', "('phu2', 'phu1'), ('cs2', 'cs1'),", "('phu2', 'phu1'), ('cs2', 'phu1'),"),
    ('P16', 'COMPARISONS2_POOLED: the ISS against the round-1 photutils twice', "(('iss', 'phu1'), ('iss', 'cs1'), ('iss_adq', 'phu1'), ('iss_adq', 'cs1'))", "(('iss', 'phu1'), ('iss', 'phu1'), ('iss_adq', 'phu1'), ('iss_adq', 'cs1'))"),
    # ---- the development table
    ('E01', 'dev_table2: the selected candidate not marked', "mark = '**' if sel.get(fam) == tag else ''", "mark = ''"),
    ('E02', 'dev_table2: the round-1 selection not marked', "' (round 1 selection)' if sel1.get(fam) == tag else ''", "''"),
    ('E03', 'dev_table2: the round-1 candidate not listed', "for tag in tags + ([sel1[fam]] if fam in sel1 else []):", "for tag in tags:"),
    ('E04', 'dev_table2: the seconds of the first patch only', "for p in X.DEV if 'DEV|x%d_y%d|%s|secs' % (p[0], p[1], tag) in raw)", "for p in X.DEV[:1] if 'DEV|x%d_y%d|%s|secs' % (p[0], p[1], tag) in raw)"),
    ('E05', 'dev_table2: the detections of all magnitudes', "'%d' % tot[2],\n                         '%d' % (tot[6] - tot[5]), '%.0f' % secs])", "'%d' % tot[6],\n                         '%d' % (tot[6] - tot[5]), '%.0f' % secs])"),
    ('E06', 'dev_table2: the old rows pooled over one patch', "if all('OLD|x%d_y%d|%s|tot' % (p[0], p[1], label) in raw1 for p in X.DEV):", "if any('OLD|x%d_y%d|%s|tot' % (p[0], p[1], label) in raw1 for p in X.DEV):"),
    ('E07', 'dev_table2: the old rows without their names', "E1.shown(raw1['OLD|%s|name' % label][0])", "label"),
    ('E08', 'describe2: the bound of a photutils variant', "return '%s, xy_bounds %g px' % (E1.describe(d['parent']), d['bounds'])", "return '%s, xy_bounds %g px' % (E1.describe(d['parent']), d['bounds'] + 1)"),
    ('E09', 'describe2: the passes of a crowdsource variant', "p['miniter'], p['maxiter'], p['fewstars'])", "p['miniter'], p['maxiter'], p['fewstars'] + 1)"),
    ('E10', 'describe2: the width of a crowdsource variant', "(d['sb'], d['thr'], d['passes'],", "(d['sb'] + 1, d['thr'], d['passes'],"),
    # ---- the verdict counts and the summary
    ('C01', 'patch_verdicts: the stamp-kernel row is always the isotropic one', "('stamp', 'iss_adq' if 'iss_adq' in g['blk'] else 'iss')", "('stamp', 'iss')"),
    ('C02', 'patch_verdicts: the stamp-kernel row is required', "('stamp', 'iss_adq' if 'iss_adq' in g['blk'] else 'iss')", "('stamp', 'iss_adq')"),
    ('C03', 'patch_verdicts: the bootstrap of another seed', "idx = X.boot_idx_strat(g['sizes'])\n            for kind, row in", "idx = X.boot_idx_strat(g['sizes'], seed=1)\n            for kind, row in"),
    ('C04', 'patch_verdicts: second minus first', "X.verdict(*X.boot_diff(g['blk'][row], g['blk'][base], idx)[2])", "X.verdict(*X.boot_diff(g['blk'][base], g['blk'][row], idx)[2])"),
    ('C05', 'patch_verdicts: the baselines of round 1 left out', "for base in ('phu2', 'cs2', 'phu1', 'cs1'):\n                    if row in", "for base in ('phu2', 'cs2'):\n                    if row in"),
    ('C06', 'verdict_table: wins and losses exchanged in the counts', "for w in ('win', 'tie', 'loss')) if n in per else '-' for n in names]", "for w in ('loss', 'tie', 'win')) if n in per else '-' for n in names]"),
    ('C07', 'verdict_table: the total counts ties as wins', "'%d / %d / %d (%d patches)' % (tuple(sum(1 for v in allv if v == w) for w in ('win', 'tie', 'loss')) + (len(allv),))]", "'%d / %d / %d (%d patches)' % (tuple(sum(1 for v in allv if v == w) for w in ('tie', 'win', 'loss')) + (len(allv),))]"),
    ('C08', 'verdict_table: the stamp-kernel row labelled as the isotropic one', "'tuned ISS, isotropic kernel as frozen' if kind == 'iso' else", "'tuned ISS, isotropic kernel as frozen' if kind != 'iso' else"),
    ('C09', 'summary_tables2: the difference of the crowdsource against the photutils round 1', "('phu2', 'phu1'), ('cs2', 'cs1'), ('phu2', 'phu_old_t'), ('cs2', 'cs_old_t'))])", "('phu2', 'phu1'), ('cs2', 'phu1'), ('phu2', 'phu_old_t'), ('cs2', 'cs_old_t'))])"),
    ('C10', 'summary_tables2: unmatched detections of mag < 22', "return '%.3f [%.3f, %.3f] (%d)' % (ff, ci[0], ci[1], g['tot'][r][6] - g['tot'][r][5])", "return '%.3f [%.3f, %.3f] (%d)' % (ff, ci[0], ci[1], g['tot'][r][2] - g['tot'][r][3])"),
    ('C11', 'summary_tables2: the F columns in another order', "t1.append([title + ' (%d patches)' % len(g['done'])] + [cell(r) for r in ORDER2])", "t1.append([title + ' (%d patches)' % len(g['done'])] + [cell(r) for r in ORDER2[::-1]])"),
    # ---- the report
    ('R01', 'pooled_block2: the patches pooled without stratification', "def pooled_block2(heading, g, raw_summary, title):\n    idx = X.boot_idx_strat(g['sizes'])", "def pooled_block2(heading, g, raw_summary, title):\n    idx = X.boot_idx(sum(g['sizes']))"),
    ('R02', 'report2: the per-patch bootstrap with another seed', "idx = X.boot_idx_strat(gp['sizes'])", "idx = X.boot_idx_strat(gp['sizes'], seed=1)"),
    ('R03', 'report2: the pooled block shown for one patch too', "        if len(g['done']) > 1:\n            body += pooled_block2('### All %d test patches together'", "        if len(g['done']) > 0:\n            body += pooled_block2('### All %d test patches together'"),
    ('R04', 'report2: the repeat patches pooled only from four', "if len(gr['done']) > 1:", "if len(gr['done']) > 3:"),
    ('R05', 'report2: the repeat patches pooled over two frames', "gr = group_data2(raw, raw1, [(n, G.SETS[n][1]) for n in sorted(R.REPEAT)], frozen)", "gr = group_data2(raw, raw1, [(n, G.SETS[n][1]) for n in sorted(R.REPEAT)[:2]], frozen)"),
    ('R06', 'report2: the truth count of the patch from the wrong column', "gp['tot']['phu2'][0]", "gp['tot']['phu2'][1]"),
    ('R07', 'report2: the plan hash not shown', "% (str(raw['SEL|plan_sha256'][0]) if 'SEL|plan_sha256' in raw else 'none yet', X.BOOT), '']", "% ('none yet', X.BOOT), '']"),
    ('R08', 'report2: the shas of the frozen results exchanged', "(('m2_raw.npz', shas.get('m2')), ('psf_raw.npz', shas.get('psf')), ('repeat_raw.npz', shas.get('repeat')))", "(('m2_raw.npz', shas.get('psf')), ('psf_raw.npz', shas.get('m2')), ('repeat_raw.npz', shas.get('repeat')))"),
    ('R09', 'report2: the selection shown exchanged', "% (raw['SEL|photutils'][0], raw['SEL|r1_photutils'][0], raw['SEL|crowdsource'][0], raw['SEL|r1_crowdsource'][0])", "% (raw['SEL|crowdsource'][0], raw['SEL|r1_photutils'][0], raw['SEL|photutils'][0], raw['SEL|r1_crowdsource'][0])"),
    ('R10', 'report2: the check of the round-1 variant always says identical', "'identical' if same else 'DIFFERENT'", "'identical'"),
    ('R11', 'report2: the check compares one patch only', "same = all(np.array_equal(raw[k % p], raw['DEV|x%d_y%d|%s|tot' % (p[0], p[1], raw['SEL|r1_crowdsource'][0])]) for p in X.DEV)", "same = all(np.array_equal(raw[k % p], raw['DEV|x%d_y%d|%s|tot' % (p[0], p[1], raw['SEL|r1_crowdsource'][0])]) for p in X.DEV[:1])"),
    ('R12', 'report2: the 95 % text', "Paired differences (first minus second, same blocks), 95 % intervals", "Paired differences (first minus second, same blocks), 95 %% intervals"),
    ('R13', 'report2: no verdict counts', "        vd = patch_verdicts(raw, raw1, frozen)\n        if vd:", "        vd = {}\n        if vd:"),
    ('R14', 'report2: the single-patch group missing from the summary', "        else:\n            idx = X.boot_idx_strat(g['sizes'])\n            summary.append(", "        else:\n            idx = X.boot_idx_strat(g['sizes'])\n            (lambda *a: None)("),
    ('B01', 'main: the test sets default to the M2 patches only', "ap.add_argument('--sets', nargs='*', default=['m2', 'm53', 'm92'])", "ap.add_argument('--sets', nargs='*', default=['m2'])"),
]
