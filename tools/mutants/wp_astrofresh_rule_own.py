"""Single-point mutants of experiments/astro_fresh.py (the rule that chooses the fresh clusters, the registry loaders), for tools/mutants.py:

    venv/bin/python tools/mutants.py --exp astro_fresh --tests tests/test_exp_astro_fresh.py --raw results/astro/fresh_raw.npz --md results/astro/fresh.md --out-env FRESH_OUT \\
        --spec tools/mutants/wp_astrofresh_rule_own.py --work <scratch dir> [--slot a --i0 0 --step 2]

A tuple is (id, what, old, new[, k]): `old` replaced by `new` at its k-th occurrence in experiments/astro_fresh.py. The kill criterion is the tests file alone (the same file as for `wp_astrofresh_own.py`:
`tests/test_exp_astro_fresh.py`). `astro_fresh.py` has no `--from-raw`, so the md column of the tool is meaningless for this spec (always unchanged, the regeneration exits with an argument error).
The code paths that need astropy and the downloaded data (`survey_fields`, `cmd_fields`, `cmd_select`, `cmd_select_frames`, `cmd_build`, `cmd_registry`) are not exercised by the tests and are not mutated: they were run once on the real
data and their outputs are the stored files `fresh_survey.txt`, `fresh_fields.json`, `fresh_selection.txt`, `fresh_selection_A1.txt`, `fresh_patches.json` (a limit, stated in the report). The constant `RULE['min_patches']` is 2 (version 1 of the rule); amendment A2 (a cluster with
one patch stays) was applied by hand when the clusters were built, the code was not changed after the freeze record was written, so the test pins the constant as it is.

Final run (complete spec, final tests file): 78 mutants, 78 killed, no survivor. A first pass against the first version of the tests left 7 survivors (S09, A02, A12, A13, H02, H04, L06), each killed by a test added for it.
"""
M = [
    # ---- the constants
    ('C01', 'RULE: the PSF status up to 3', "status_max=2,", "status_max=3,"),
    ('C02', 'RULE: no PSF star needed', "nstar_min=1,", "nstar_min=0,"),
    ('C03', 'RULE: the lower seeing limit 2.0', "fwhm_lo=2.3,", "fwhm_lo=2.0,"),
    ('C04', 'RULE: the upper seeing limit 4.0', "fwhm_hi=3.6,", "fwhm_hi=4.0,"),
    ('C05', 'RULE: 50 windows are enough', "windows_min=100,", "windows_min=50,"),
    ('C06', 'RULE: three fresh clusters', "n_fresh=2)", "n_fresh=3)"),
    ('C07', 'RULE: a single patch is enough (amendment A2 as a constant)', "min_patches=2,", "min_patches=1,"),
    ('C08', 'RULE: the frame centre elsewhere', "centre=(1023.5, 744.0)", "centre=(1024.0, 744.0)"),
    ('C09', 'UNIVERSE: M53 included', "UNIVERSE = ('ngc5053',", "UNIVERSE = ('ngc5024', 'ngc5053',"),
    ('C10', 'UNIVERSE: M92 included', "'ngc7078')\nRULE", "'ngc7078', 'ngc6341')\nRULE"),
    ('C11', 'UNIVERSE: M15 left out', "'ngc6981', 'ngc7078')\nRULE", "'ngc6981')\nRULE"),
    # ---- the survey parser
    ('S01', 'SURVEY_LINE: the camcol of two digits', r"\|\s+(\d+)\s+(\d)\s+(\d+)\s+\|", r"\|\s+(\d+)\s+(\d\d)\s+(\d+)\s+\|"),
    ('S02', 'SURVEY_LINE: the density column must be a number', r"(\d+)\s+(\S+)\s+\|\s+psp_status", r"(\d+)\s+(\d+)\s+\|\s+psp_status"),
    ('S03', 'parse_survey: the field and the camcol exchanged', "c, run, cam, fld, see, frac, nw, dens, ps, nst, sat = m.groups()", "c, run, fld, cam, see, frac, nw, dens, ps, nst, sat = m.groups()"),
    ('S04', 'parse_survey: the status and the PSF stars exchanged', "status=int(ps), nstar=int(nst)", "status=int(nst), nstar=int(ps)"),
    ('S05', 'parse_survey: the seeing and the footprint fraction exchanged', "seeing=float(see), frac=float(frac)", "seeing=float(frac), frac=float(see)"),
    ('S06', 'parse_survey: the windows from the status', "windows=int(nw)", "windows=int(ps)"),
    ('S07', 'parse_survey: the saturation level as an integer string', "sat=float(sat)", "sat=sat"),
    ('S08', 'parse_survey: the density string dropped', "density=dens,", "density='',"),
    ('S09', 'SURVEY_LINE: matches lines without a field box', r"\|\s+psp_status\s+(\d+)\s+psf_nstar\(r\)\s+(\d+)\s+sat\s+(\d+)')", r"\|?\s*psp_status?\s*(\d*)\s*psf_nstar\(r\)?\s*(\d*)\s*sat?\s*(\d*)')"),
    # ---- the rule
    ('A01', 'apply_rule: PSF status below the limit only', "r['status'] <= rule['status_max']", "r['status'] < rule['status_max']"),
    ('A02', 'apply_rule: PSF stars above the minimum only', "r['nstar'] >= rule['nstar_min']", "r['nstar'] > rule['nstar_min']"),
    ('A03', 'apply_rule: status or stars', "r['status'] <= rule['status_max'] and r['nstar'] >= rule['nstar_min']", "r['status'] <= rule['status_max'] or r['nstar'] >= rule['nstar_min']"),
    ('A04', 'apply_rule: the lower seeing limit exclusive', "rule['fwhm_lo'] <= f <= rule['fwhm_hi']", "rule['fwhm_lo'] < f <= rule['fwhm_hi']"),
    ('A05', 'apply_rule: the upper seeing limit exclusive', "rule['fwhm_lo'] <= f <= rule['fwhm_hi']", "rule['fwhm_lo'] <= f < rule['fwhm_hi']"),
    ('A06', 'apply_rule: an unmeasured stamp passes', "r['ok_seeing'] = f is not None and rule['fwhm_lo']", "r['ok_seeing'] = f is None or rule['fwhm_lo']"),
    ('A07', 'apply_rule: the windows limit exclusive', "r['ok_windows'] = r['windows'] >= rule['windows_min']", "r['ok_windows'] = r['windows'] > rule['windows_min']"),
    ('A08', 'apply_rule: a cluster passes with any two of the three tests', "r['passes'] = bool(r['ok_psf'] and r['ok_seeing'] and r['ok_windows'])", "r['passes'] = bool(r['ok_psf'] + r['ok_seeing'] + r['ok_windows'] >= 2)"),
    ('A09', 'apply_rule: the seeing test not needed', "r['passes'] = bool(r['ok_psf'] and r['ok_seeing'] and r['ok_windows'])", "r['passes'] = bool(r['ok_psf'] and r['ok_windows'])"),
    ('A10', 'apply_rule: the fewest windows first', "(not r['passes'], -r['windows'] if r['passes'] else 0, ngc(r))", "(not r['passes'], r['windows'] if r['passes'] else 0, ngc(r))"),
    ('A11', 'apply_rule: a tie goes to the larger NGC number', "(not r['passes'], -r['windows'] if r['passes'] else 0, ngc(r))", "(not r['passes'], -r['windows'] if r['passes'] else 0, -ngc(r))"),
    ('A12', 'apply_rule: the failing clusters first', "(not r['passes'], -r['windows']", "(r['passes'], -r['windows']"),
    ('A13', 'apply_rule: the NGC number read as text', "ngc = lambda r: int(r['cluster'][3:])\n    out.sort", "ngc = lambda r: r['cluster']\n    out.sort"),
    ('A14', 'apply_rule: the clusters without a row are not reported', "        if c not in byc:\n            missing.append(c)\n            continue", "        if c not in byc:\n            continue"),
    ('A15', 'apply_rule: the stamp of another cluster', "f = fwhm.get(c)", "f = fwhm.get('ngc5272')"),
    ('A16', 'apply_rule: clusters outside the universe ranked', "    for c in universe:\n        if c not in byc:", "    for c in byc:\n        if c not in byc:"),
    # ---- the choice among the passing clusters
    ('H01', 'choose_fresh: a third cluster taken', "if not r['passes'] or len(chosen) >= rule['n_fresh']:", "if not r['passes'] or len(chosen) > rule['n_fresh']:"),
    ('H02', 'choose_fresh: the failing clusters considered', "if not r['passes'] or len(chosen) >= rule['n_fresh']:", "if len(chosen) >= rule['n_fresh']:"),
    ('H03', 'choose_fresh: exactly the minimum of patches is too few', "(chosen if n >= rule['min_patches'] else dropped)", "(chosen if n > rule['min_patches'] else dropped)"),
    ('H04', 'choose_fresh: a cluster not yet built is skipped instead of stopping the walk', "        if n is None:\n            break", "        if n is None:\n            continue"),
    ('H05', 'choose_fresh: dropped clusters counted as chosen', "(chosen if n >= rule['min_patches'] else dropped).append(r['cluster'])", "chosen.append(r['cluster'])"),
    ('H06', 'choose_fresh: the default rule ignored', "    rule = rule or RULE\n    chosen, dropped = [], []", "    rule = RULE\n    chosen, dropped = [], []"),
    # ---- the frames (amendment A1)
    ('B01', 'apply_rule_frames: the frame with the fewest windows', "key = (-r['windows'], r['seeing'], r['run'], r['camcol'], r['field'])", "key = (r['windows'], r['seeing'], r['run'], r['camcol'], r['field'])"),
    ('B02', 'apply_rule_frames: a tie goes to the larger PSF width', "key = (-r['windows'], r['seeing'], r['run'], r['camcol'], r['field'])", "key = (-r['windows'], -r['seeing'], r['run'], r['camcol'], r['field'])"),
    ('B03', 'apply_rule_frames: a failing frame can be the cluster\'s frame', "    for r in allf:\n        if r['passes']:", "    for r in allf:\n        if True:"),
    ('B04', 'apply_rule_frames: the clusters ordered by their smallest frame', "ranked = sorted((v[1] for v in best.values()), key=lambda r: (-r['windows'], int(r['cluster'][3:])))", "ranked = sorted((v[1] for v in best.values()), key=lambda r: (r['windows'], int(r['cluster'][3:])))"),
    ('B05', 'apply_rule_frames: frames keyed by the cluster alone are not told apart by the stamp', "f = fwhm.get(frame_id(r))", "f = fwhm.get(r['cluster'])"),
    ('B06', 'frame_id: the camcol and the field exchanged', "return '%s:%d:%d:%d' % (r['cluster'], r['run'], r['camcol'], r['field'])", "return '%s:%d:%d:%d' % (r['cluster'], r['run'], r['field'], r['camcol'])"),
    ('B07', 'need_psfield: the windows test dropped', "r['nstar'] >= rule['nstar_min'] and r['windows'] >= rule['windows_min']]", "r['nstar'] >= rule['nstar_min']]"),
    ('B08', 'need_psfield: the PSF test dropped', "r['status'] <= rule['status_max'] and r['nstar'] >= rule['nstar_min'] and r['windows']", "r['windows']".join(["True and ", ""]) if False else "r['status'] <= rule['status_max'] and r['windows']"),
    ('B09', 'need_psfield: the universe ignored', "if r['cluster'] in universe and r['status']", "if r['status']"),
    ('B10', 'apply_rule_frames: the universe ignored', "        if r['cluster'] not in universe:\n            continue\n        r = dict(r)", "        r = dict(r)"),
    # ---- the patches' stamp filter and the loaders
    ('L01', 'fwhm_filter: the lower limit exclusive', "if rule['fwhm_lo'] <= f <= rule['fwhm_hi']:", "if rule['fwhm_lo'] < f <= rule['fwhm_hi']:"),
    ('L02', 'fwhm_filter: the upper limit exclusive', "if rule['fwhm_lo'] <= f <= rule['fwhm_hi']:", "if rule['fwhm_lo'] <= f < rule['fwhm_hi']:"),
    ('L03', 'fwhm_filter: the stamp read at the origin instead of the centre of the patch', "f = stamp_fwhm(kl, q['x0'] + (S2.N_PATCH - 1) / 2.0, q['y0'] + (S2.N_PATCH - 1) / 2.0)", "f = stamp_fwhm(kl, q['x0'], q['y0'])"),
    ('L04', 'fwhm_filter: x and y exchanged', "f = stamp_fwhm(kl, q['x0'] + (S2.N_PATCH - 1) / 2.0, q['y0'] + (S2.N_PATCH - 1) / 2.0)", "f = stamp_fwhm(kl, q['y0'] + (S2.N_PATCH - 1) / 2.0, q['x0'] + (S2.N_PATCH - 1) / 2.0)"),
    ('L05', 'fwhm_filter: the FWHM not attached', "kept.append(dict(q, fwhm=f))", "kept.append(dict(q))"),
    ('L06', 'fwhm_filter: the removed count wrong', "return kept, len(cands) - len(kept)", "return kept, len(kept)"),
    ('L07', 'patch_path: the origin swapped', "'%s_x%d_y%d_100x100' % (entry['prefix'], x0, y0)", "'%s_x%d_y%d_100x100' % (entry['prefix'], y0, x0)"),
    ('L08', 'patch_path: the explicit root ignored', "return os.path.join(root or entry['root'],", "return os.path.join(entry['root'],"),
    ('L09', 'load_patch: the PSF read from the 25 px stamp', "psf=z['psf51'].astype(float)", "psf=z['psf25'].astype(float)"),
    ('L10', 'load_patch: gain and dark variance exchanged', "gain=float(z['gain']), dark=float(z['dark_variance']))", "gain=float(z['dark_variance']), dark=float(z['gain']))"),
    ('L11', 'truth_arrays: x and y exchanged', "return dict(px=z['px'], py=z['py'], m606=z['m606'])", "return dict(px=z['py'], py=z['px'], m606=z['m606'])"),
    ('L12', 'load_registry: the registry path ignored', "with open(path or REGISTRY) as f:", "with open(REGISTRY) as f:"),
    ('L13', 'psfield_path: the field written with three digits', "'psField-%06d-%d-%04d.fit' % (run, camcol, field))", "'psField-%06d-%d-%03d.fit' % (run, camcol, field))"),
    ('L14', 'frame_prefix: the camcol and the field exchanged', "return 'patch_r_%06d_%d_%04d' % (run, camcol, field)", "return 'patch_r_%06d_%d_%04d' % (run, field, camcol)"),
    # the same tests inside apply_rule_frames (k = 1) and need_psfield (k = 2): the lines are repeated in the source
    ('X01', 'apply_rule_frames: PSF status below the limit only', "r['status'] <= rule['status_max']", "r['status'] < rule['status_max']", 1),
    ('X02', 'apply_rule_frames: PSF stars above the minimum only', "r['nstar'] >= rule['nstar_min']", "r['nstar'] > rule['nstar_min']", 1),
    ('X03', 'apply_rule_frames: status or stars', "r['status'] <= rule['status_max'] and r['nstar'] >= rule['nstar_min']", "r['status'] <= rule['status_max'] or r['nstar'] >= rule['nstar_min']", 1),
    ('X04', 'apply_rule_frames: the lower seeing limit exclusive', "rule['fwhm_lo'] <= f <= rule['fwhm_hi']", "rule['fwhm_lo'] < f <= rule['fwhm_hi']", 1),
    ('X05', 'apply_rule_frames: the upper seeing limit exclusive', "rule['fwhm_lo'] <= f <= rule['fwhm_hi']", "rule['fwhm_lo'] <= f < rule['fwhm_hi']", 1),
    ('X06', 'apply_rule_frames: an unmeasured stamp passes', "r['ok_seeing'] = f is not None and rule['fwhm_lo']", "r['ok_seeing'] = f is None or rule['fwhm_lo']", 1),
    ('X07', 'apply_rule_frames: the windows limit exclusive', "r['ok_windows'] = r['windows'] >= rule['windows_min']", "r['ok_windows'] = r['windows'] > rule['windows_min']", 1),
    ('X08', 'apply_rule_frames: a frame passes with any two of the three tests', "r['passes'] = bool(r['ok_psf'] and r['ok_seeing'] and r['ok_windows'])", "r['passes'] = bool(r['ok_psf'] + r['ok_seeing'] + r['ok_windows'] >= 2)", 1),
    ('X09', 'apply_rule_frames: the seeing test not needed', "r['passes'] = bool(r['ok_psf'] and r['ok_seeing'] and r['ok_windows'])", "r['passes'] = bool(r['ok_psf'] and r['ok_windows'])", 1),
    ('X10', 'need_psfield: PSF status below the limit only', "r['status'] <= rule['status_max']", "r['status'] < rule['status_max']", 2),
    ('X11', 'need_psfield: PSF stars above the minimum only', "r['nstar'] >= rule['nstar_min']", "r['nstar'] > rule['nstar_min']", 2),
    ('X12', 'need_psfield: the windows limit exclusive', "r['windows'] >= rule['windows_min']]", "r['windows'] > rule['windows_min']]"),
]
