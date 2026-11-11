"""Single-point mutants of experiments/exp_smlm_decode.py, for tools/mutants.py:

    venv/bin/python tools/mutants.py --exp exp_smlm_decode --tests tests/test_exp_smlm_decode.py --raw results/smlm/decode_raw.npz --md results/smlm/decode.md --out-env DECODE_OUT \\
        --extra results/smlm/learned_raw.npz --extra results/smlm/decode.md --spec tools/mutants/wp_smlmdecode_own.py --work <scratch dir> [--slot a --i0 0 --step 2]

A tuple is (id, what, old, new[, k]): `old` replaced by `new` at its k-th occurrence in experiments/exp_smlm_decode.py. The kill criterion is the tests file alone (`pytest -x`; the stored-md test reads the copies of the
raw and the md that `--extra` puts in the slot, everything else is planted). Ids: K the constants, C the DECODE-to-ground-truth convention, P the packing of the candidates, R the raw, X the consistency check of the raw,
V the verdict and the formatting, M the method table, S the summary (and the tally), T the per-set report, B the binned descriptive tables, H the threshold table, N the ingest of a run (the text composed from the
returned files), L the comparison rows read from learned_raw.npz, I the command line.

Result on the 16 tests of tests/test_exp_smlm_decode.py (the unmutated source passes the control with them and regenerates the md unchanged): 97 mutants, 95 killed, 2 survivors that change no output (equivalent
mutants, left in on purpose):
  K03  `OPS = ('J', 'E')` for `('E', 'J')`: the operating points are looped over to check the stored detections and to build the raw; the order changes nothing that is stored or printed;
  P11  the ground-truth `id` field set to zeros in `gt_array`: neither the matching nor the scoring of `smlm_common` reads it.
Order of work: a first complete pass of 88 mutants, on a planted raw (the stored run was not there yet), left 26 survivors; the killable ones got tests (an E point that differs from the J point, so that the stored
arrays, the offsets, the method table, the summary and the threshold table can be told apart; two sets for the mean row and the tally; the boundaries of the verdict and of the frame range; the photon threshold
and the ground truth that the raw keeps; a planted learned raw for the CNN rows; the text composed from the returned files), 9 more mutants were added for the tally and that text, and the 97 were run again on
the stored run (94 killed; N11, the second GPU line, was then killed by a test with two GPUs and re-run).

Not mutated, on purpose: `ingest`'s reading of the SMLM data (`C.load_gt(C.data_path(...))`, the default `gt_loader`) and of `results/smlm/learned_raw.npz` (the tests plant both), `main`'s argument handling beyond
the two switches above, and the wording of the protocol bullets of `compose_meta` (quoted from `decode_run_plan.md`; the tests check the computed parts of the text, not the fixed sentences).
"""
M = [
    # ---- K: constants
    ('K01', 'the bootstrap size', "BOOT, BLOCK = 1000, C.DEFAULT_BLOCK", "BOOT, BLOCK = 500, C.DEFAULT_BLOCK"),
    ('K02', 'the bootstrap block', "BOOT, BLOCK = 1000, C.DEFAULT_BLOCK", "BOOT, BLOCK = 1000, 5"),
    ('K03', 'the operating points in another order', "OPS = ('E', 'J')", "OPS = ('J', 'E')"),
    ('K04', 'the raw key tag', "TAG = 'DECODE'", "TAG = 'DECODE2'"),
    ('K05', 'the candidate columns exchanged', "CAND_COLS = ('frame', 'x_nm', 'y_nm', 'photons', 'prob')", "CAND_COLS = ('frame', 'y_nm', 'x_nm', 'photons', 'prob')"),
    ('K06', 'the probability and photons exchanged', "CAND_COLS = ('frame', 'x_nm', 'y_nm', 'photons', 'prob')", "CAND_COLS = ('frame', 'x_nm', 'y_nm', 'prob', 'photons')"),
    # ---- C: convention
    ('C01', 'no axis swap in the conversion (x from x_px)', "return (np.asarray(y_px, dtype=np.float64) + 0.5) * C.PIXEL_NM, (np.asarray(x_px, dtype=np.float64) + 0.5) * C.PIXEL_NM",
     "return (np.asarray(x_px, dtype=np.float64) + 0.5) * C.PIXEL_NM, (np.asarray(y_px, dtype=np.float64) + 0.5) * C.PIXEL_NM"),
    ('C02', 'no half-pixel shift', "return (np.asarray(y_px, dtype=np.float64) + 0.5) * C.PIXEL_NM", "return (np.asarray(y_px, dtype=np.float64) + 0.0) * C.PIXEL_NM"),
    ('C03', 'no half-pixel shift in y', "(np.asarray(x_px, dtype=np.float64) + 0.5) * C.PIXEL_NM\n", "(np.asarray(x_px, dtype=np.float64) + 0.0) * C.PIXEL_NM\n"),
    ('C04', 'a pixel of 50 nm', "(np.asarray(y_px, dtype=np.float64) + 0.5) * C.PIXEL_NM,", "(np.asarray(y_px, dtype=np.float64) + 0.5) * 50.0,"),
    # ---- P: packing
    ('P01', 'the reshape of the candidates', "rows = np.asarray(rows, dtype=np.float64).reshape(-1, 5)\n    out = np.zeros", "rows = np.asarray(rows, dtype=np.float64).reshape(-1, 5)[:-1]\n    out = np.zeros"),
    ('P02', 'to_struct drops a column', "        out[c] = rows[:, i]", "        out[c] = rows[:, 0]"),
    ('P03', 'a negative frame accepted', "if len(fr) and (fr.min() < 0 or fr.max() >= n_frames):", "if len(fr) and (fr.min() < -1 or fr.max() >= n_frames):"),
    ('P04', 'the last frame index outside', "if len(fr) and (fr.min() < 0 or fr.max() >= n_frames):", "if len(fr) and (fr.min() < 0 or fr.max() > n_frames):"),
    ('P05', 'the frame is not sorted', "order = np.argsort(fr, kind='stable')\n    rows, fr = rows[order], fr[order]", "order = np.arange(len(fr))\n    rows, fr = rows[order], fr[order]"),
    ('P06', 'an unstable sort', "order = np.argsort(fr, kind='stable')\n    rows, fr = rows[order], fr[order]", "order = np.argsort(-fr, kind='stable')\n    rows, fr = rows[order], fr[order]"),
    ('P07', 'the offsets count the frame itself', "np.add.at(off, fr + 1, 1)", "np.add.at(off, fr, 1)"),
    ('P08', 'the photons dropped from the flat array', "return np.cumsum(off), rows[:, 1:5].copy()", "return np.cumsum(off), rows[:, 1:4].copy()"),
    ('P09', 'the x column dropped', "return np.cumsum(off), rows[:, 1:5].copy()", "return np.cumsum(off), rows[:, 2:6].copy()"),
    ('P10', 'the frame index of unpack shifted', "fr = np.repeat(np.arange(len(off) - 1), np.diff(off))", "fr = np.repeat(np.arange(1, len(off)), np.diff(off))"),
    ('P11', 'the ground-truth id', "a['id'] = np.arange(1, len(a) + 1)", "a['id'] = np.zeros(len(a), dtype=np.int64)"),
    ('P12', 'the ground truth without z', "    for k in ('frame', 'x_nm', 'y_nm', 'z_nm', 'photons'):\n        a[k] = g[k]", "    for k in ('frame', 'x_nm', 'y_nm', 'photons'):\n        a[k] = g[k]"),
    # ---- R: raw
    ('R01', 'the development and test ranges exchanged', "(da, db), (ta, tb) = protocol(name)['dev'], protocol(name)['test']\n        arrays['%s|dev_frames'", "(ta, tb), (da, db) = protocol(name)['dev'], protocol(name)['test']\n        arrays['%s|dev_frames'"),
    ('R02', 'the development candidates packed with the test range', "for split, (a, b) in (('dev', (da, db)), ('test', (ta, tb))):", "for split, (a, b) in (('dev', (ta, tb)), ('test', (ta, tb))):"),
    ('R03', 'the stored E detections are the J ones', "off, flat = IN.pack(o['test_dets'])\n            arrays['%s|%s-%s|test|off'", "off, flat = IN.pack(res['J']['test_dets'])\n            arrays['%s|%s-%s|test|off'"),
    ('R04', 'the photon threshold of the challenge setting lost', "arrays['%s|pmin' % name] = np.array(s['pmin'])", "arrays['%s|pmin' % name] = np.array(0.0)"),
    ('R05', 'the comparison rows packed with fewer columns', "off, flat = C.pack_detections(dets, columns=(0, 1, 2, 3))", "off, flat = C.pack_detections(dets, columns=(0, 1, 2))"),
    ('R06', 'the ground truth stored without the dev frames of another key', "arrays['%s|gt|%s' % (name, k)] = np.asarray(s['gt'][k])", "arrays['%s|gt|%s' % (name, k)] = np.asarray(s['gt'][k])[::-1]"),
    ('R07', 'the candidates unpacked from the test key for the development frames', "def cand_of(raw, name, split):\n    return unpack_cand(raw['%s|%s|%s|cand|off' % (name, TAG, split)], raw['%s|%s|%s|cand|flat' % (name, TAG, split)])",
     "def cand_of(raw, name, split):\n    return unpack_cand(raw['%s|%s|test|cand|off' % (name, TAG)], raw['%s|%s|test|cand|flat' % (name, TAG)])"),
    # ---- X: consistency check
    ('X01', 'the check passes when either array matches', "if not (np.array_equal(off, raw['%s|%s-%s|test|off' % (name, TAG, op)]) and np.array_equal(flat, raw['%s|%s-%s|test|flat' % (name, TAG, op)])):",
     "if not (np.array_equal(off, raw['%s|%s-%s|test|off' % (name, TAG, op)]) or np.array_equal(flat, raw['%s|%s-%s|test|flat' % (name, TAG, op)])):"),
    ('X02', 'the offsets are not checked', "if not (np.array_equal(off, raw['%s|%s-%s|test|off' % (name, TAG, op)]) and np.array_equal(flat,", "if not (True and np.array_equal(flat,"),
    ('X03', 'the J point is not checked', "        for op in OPS:\n            off, flat = IN.pack(res[op]['test_dets'])", "        for op in OPS[:1]:\n            off, flat = IN.pack(res[op]['test_dets'])"),
    ('X04', 'the check never raises', "raise ValueError('%s %s: the stored test detections differ from those recomputed from the candidates' % (name, op))", "pass"),
    # ---- V: verdict and formatting
    ('V01', 'a win at the lower bound 0', "return 'win' if ci[0] > 0 else ('loss' if ci[1] < 0 else 'tie')", "return 'win' if ci[0] >= 0 else ('loss' if ci[1] < 0 else 'tie')"),
    ('V02', 'a loss at the upper bound 0', "return 'win' if ci[0] > 0 else ('loss' if ci[1] < 0 else 'tie')", "return 'win' if ci[0] > 0 else ('loss' if ci[1] <= 0 else 'tie')"),
    ('V03', 'win and loss exchanged', "return 'win' if ci[0] > 0 else ('loss' if ci[1] < 0 else 'tie')", "return 'loss' if ci[0] > 0 else ('win' if ci[1] < 0 else 'tie')"),
    ('V04', 'the interval formatted with the wrong scale', "return '%.*f [%.*f, %.*f]' % (nd, scale * v, nd, scale * ci[0], nd, scale * ci[1])", "return '%.*f [%.*f, %.*f]' % (nd, scale * v, nd, ci[0], nd, scale * ci[1])"),
    ('V05', 'percentages scaled by 10', "def pct(v, ci, nd=1):\n    return fmt_ci(v, ci, 100.0, nd)", "def pct(v, ci, nd=1):\n    return fmt_ci(v, ci, 10.0, nd)"),
    ('V06', 'the table separator', "'|' + '|'.join(['---'] * len(header)) + '|'", "'|' + '|'.join(['--'] * len(header)) + '|'"),
    # ---- M: methods
    ('M01', 'the J point taken from the E result', "('DECODE-J', 'DECODE, probability threshold and offset chosen for J on the development frames', res['J']['test_dets'])]", "('DECODE-J', 'DECODE, probability threshold and offset chosen for J on the development frames', res['E']['test_dets'])]"),
    ('M02', 'the E point taken from the J result', "ms = [('DECODE-E', 'DECODE, probability threshold and offset chosen for E on the development frames', res['E']['test_dets']),", "ms = [('DECODE-E', 'DECODE, probability threshold and offset chosen for E on the development frames', res['J']['test_dets']),"),
    ('M03', 'the comparison rows missing', "        d = cmp_dets(raw, name, tag)\n        if d is not None:\n            ms.append((tag, lab, d))", "        d = cmp_dets(raw, name, tag)\n        if d is None:\n            ms.append((tag, lab, d))"),
    # ---- S: summary
    ('S01', 'the challenge border of the report', "r2 = C.score(dets, gt, idx, border_nm=450.0, photon_min=pmin)", "r2 = C.score(dets, gt, idx, border_nm=400.0, photon_min=pmin)"),
    ('S02', 'no photon threshold in the challenge setting', "r2 = C.score(dets, gt, idx, border_nm=450.0, photon_min=pmin)", "r2 = C.score(dets, gt, idx, border_nm=450.0, photon_min=None)"),
    ('S03', 'the bootstrap seed', "R[name][key] = (lab, r, r2, C.bootstrap_ci(r, B=boot, seed=0, block=BLOCK), C.bootstrap_ci(r2, B=boot, seed=0, block=BLOCK))", "R[name][key] = (lab, r, r2, C.bootstrap_ci(r, B=boot, seed=1, block=BLOCK), C.bootstrap_ci(r2, B=boot, seed=0, block=BLOCK))"),
    ('S04', 'the challenge-setting seed', "R[name][key] = (lab, r, r2, C.bootstrap_ci(r, B=boot, seed=0, block=BLOCK), C.bootstrap_ci(r2, B=boot, seed=0, block=BLOCK))", "R[name][key] = (lab, r, r2, C.bootstrap_ci(r, B=boot, seed=0, block=BLOCK), C.bootstrap_ci(r2, B=boot, seed=1, block=BLOCK))"),
    ('S05', 'the efficiency shown as the Jaccard index', "row.append('%.1f (%.1f)' % (R[s][k][1]['efficiency'], 100 * R[s][k][1]['jaccard']) if k in R[s] else '-')", "row.append('%.1f (%.1f)' % (R[s][k][1]['jaccard'], 100 * R[s][k][1]['jaccard']) if k in R[s] else '-')"),
    ('S06', 'the summary difference of the challenge setting', "ra, rb = R[s]['DECODE-E'][1], R[s][kb][1]", "ra, rb = R[s]['DECODE-E'][2], R[s][kb][2]"),
    ('S07', 'the summary compares the J point', "ra, rb = R[s]['DECODE-E'][1], R[s][kb][1]", "ra, rb = R[s]['DECODE-J'][1], R[s][kb][1]"),
    ('S08', 'the summary interval of another metric', "d = C.bootstrap_diff_ci(ra['per_frame'], rb['per_frame'], B=boot, seed=0, block=BLOCK)['efficiency']", "d = C.bootstrap_diff_ci(ra['per_frame'], rb['per_frame'], B=boot, seed=0, block=BLOCK)['jaccard']"),
    ('S09', 'the summary pairs lose the CNN row', "('CNN-a', 'like-for-like (both use the development labels)'),\n                     ('ISS-0'", "('ISS-E', 'like-for-like (both use the development labels)'),\n                     ('ISS-0'"),
    ('S10', 'the mean over the sets drops DECODE', "for k in ('DECODE-E', 'ISS-E', 'R2E+c'))]", "for k in ('ISS-E', 'ISS-E', 'R2E+c'))]"),
    ('S11', 'the mean over sets needs only two sets', "if len(full) > 1:", "if len(full) > 2:"),
    ('S12', 'the tally compares ISS-E twice', "    for kb in ('ISS-E', 'R2E+c'):\n        v = [verdict(0,", "    for kb in ('ISS-E', 'ISS-E'):\n        v = [verdict(0,"),
    ('S13', 'the tally uses the J point', "v = [verdict(0, C.bootstrap_diff_ci(R[s]['DECODE-E'][1]['per_frame']", "v = [verdict(0, C.bootstrap_diff_ci(R[s]['DECODE-J'][1]['per_frame']"),
    ('S14', 'the tally counts losses as wins', "(kb, v.count('win'), v.count('tie'), v.count('loss')))", "(kb, v.count('loss'), v.count('tie'), v.count('win')))"),
    ('S15', 'the tally counts ties as wins', "(kb, v.count('win'), v.count('tie'), v.count('loss')))", "(kb, v.count('tie'), v.count('win'), v.count('loss')))"),
    # ---- T: per-set
    ('T01', 'the test ground truth range', "tg = gt[(gt['frame'] >= a) & (gt['frame'] < b)]", "tg = gt[(gt['frame'] >= a) & (gt['frame'] <= b)]"),
    ('T02', 'the emitters per frame', "len(tg) / float(b - a)", "len(tg) / float(b - a + 1)"),
    ('T03', 'the E point described with the J result', "op_line(name, 'E', res['E']), op_line(name, 'J', res['J']))]", "op_line(name, 'E', res['J']), op_line(name, 'J', res['J']))]"),
    ('T04', 'the dev J shown as a fraction', "100 * o['dev']['jaccard'], o['dev']['rmse_nm'], o['dev']['efficiency'])", "o['dev']['jaccard'], o['dev']['rmse_nm'], o['dev']['efficiency'])"),
    ('T05', 'the snapshot table of the wrong convention', "for r in sel['all'] if r['swap'] == sel['swap_chosen']]", "for r in sel['all'] if r['swap'] != sel['swap_chosen']]"),
    ('T06', 'the convention shown wrongly', "'x_nm from DECODE y (rows are DECODE x)' if sel.get('swap_chosen') else 'plain'", "'plain' if sel.get('swap_chosen') else 'x_nm from DECODE y (rows are DECODE x)'"),
    ('T07', 'the detected per frame', "'%.3f' % (r['n_det'] / len(idx))", "'%.3f' % (r['n_det'] / (len(idx) + 1))"),
    ('T08', 'the bias columns exchanged', "'(%.1f, %.1f)' % (r['bias_x_nm'], r['bias_y_nm'])", "'(%.1f, %.1f)' % (r['bias_y_nm'], r['bias_x_nm'])"),
    ('T09', 'the J difference compared with the DECODE-J row only', "for ka in ('DECODE-E', 'DECODE-J'):", "for ka in ('DECODE-E',):"),
    ('T10', 'the challenge difference uses the all-emitter rows', "f1(ra[2]['efficiency'] - rb[2]['efficiency'], d2['efficiency'])])", "f1(ra[1]['efficiency'] - rb[1]['efficiency'], d2['efficiency'])])"),
    ('T11', 'the RMSE difference with the wrong sign', "f1(ra[1]['rmse_nm'] - rb[1]['rmse_nm'], d['rmse_nm'])", "f1(rb[1]['rmse_nm'] - ra[1]['rmse_nm'], d['rmse_nm'])"),
    ('T12', 'the J difference of the challenge setting', "d2 = C.bootstrap_diff_ci(ra[2]['per_frame'], rb[2]['per_frame'], B=boot, seed=0, block=BLOCK)", "d2 = C.bootstrap_diff_ci(ra[1]['per_frame'], rb[1]['per_frame'], B=boot, seed=0, block=BLOCK)"),
    # ---- B: binned tables
    ('B01', 'the photon bins', "edges = np.quantile(ph, [0.0, 0.125, 0.25, 0.5, 0.75, 1.0])", "edges = np.quantile(ph, [0.0, 0.25, 0.5, 0.75, 0.9, 1.0])"),
    ('B02', 'the last photon bin closed', "    edges[-1] = np.inf\n    cols = ['%.0f-%s photons", "    cols = ['%.0f-%s photons"),
    ('B03', 'recall of the unmatched', "rows.append([m] + ['%.1f' % (100.0 * mm['matched'][(p_ >= lo) & (p_ < hi)].mean()) for lo, hi in zip(edges[:-1], edges[1:])])", "rows.append([m] + ['%.1f' % (100.0 * (1 - mm['matched'][(p_ >= lo) & (p_ < hi)].mean())) for lo, hi in zip(edges[:-1], edges[1:])])"),
    ('B04', 'the z bins', "edges = [0.0, 100.0, 200.0, 350.0, np.inf]", "edges = [0.0, 100.0, 200.0, 400.0, np.inf]"),
    ('B05', 'z signed', "z = np.abs(np.asarray(gt['z_nm']))", "z = np.asarray(gt['z_nm'])"),
    ('B06', 'the RMSE of the bins without y', "math.sqrt(float(np.mean(mm['dx'][ok] ** 2 + mm['dy'][ok] ** 2)))", "math.sqrt(float(np.mean(mm['dx'][ok] ** 2)))"),
    ('B07', 'the z table needs z', "if not np.any(gt['z_nm'] != 0):\n        return []", "if np.any(gt['z_nm'] != 0):\n        return []"),
    # ---- H: threshold table
    ('H01', 'the offset not removed in the threshold table', "d = IN.shift(IN.per_frame(to_struct(cand), len(idx), t), off)", "d = IN.per_frame(to_struct(cand), len(idx), t)"),
    ('H02', 'the J offset used', "off = res['E']['off']\n    rows, es = [], []", "off = res['J']['off']\n    rows, es = [], []"),
    ('H03', 'the best threshold by the worst E', "best = max(range(len(es)), key=lambda i: -np.inf if math.isnan(es[i]) else es[i])", "best = min(range(len(es)), key=lambda i: np.inf if math.isnan(es[i]) else es[i])"),
    ('H04', 'the bold threshold', "('**%.1f**' if t == res['E']['thr'] else '%.1f')", "('**%.1f**' if t == res['J']['thr'] else '%.1f')"),
    # ---- N: ingest
    ('N01', 'the development ground truth range of the ingest', "g = gt[(gt['frame'] >= da) & (gt['frame'] < tb)]", "g = gt[(gt['frame'] >= ta) & (gt['frame'] < tb)]"),
    ('N02', 'the test CSV read for the development candidates', "dev = read_candidates(os.path.join(run_dir, '%s_dev.csv' % name))", "dev = read_candidates(os.path.join(run_dir, '%s_test.csv' % name))"),
    ('N03', 'the photon threshold of another set', "pmin=float(z['%s|pmin' % name])", "pmin=0.0"),
    ('N04', 'the GPU name taken from the wrong line', "gpu_name = gpu[1].split(',')[0] if len(gpu) > 1 else 'unknown GPU'", "gpu_name = gpu[0].split(',')[0] if len(gpu) > 1 else 'unknown GPU'"),
    ('N05', 'the DECODE commit shortened wrongly', "decode_commit=DECODE_COMMIT[:8] + ' (' + DECODE_COMMIT + ')'", "decode_commit=DECODE_COMMIT[:7] + ' (' + DECODE_COMMIT + ')'"),
    ('N06', 'the budget text without the epochs', "'/'.join(str(eps[s_]) for s_ in sets if s_ in eps)", "'/'.join(str(0) for s_ in sets if s_ in eps)"),
    ('N07', 'the minimum of the training seconds', "(min(secs), max(secs),", "(max(secs), max(secs),"),
    ('N08', 'the conda file suffix kept in the versions', ".replace('.json', '')", ".replace('.jsn', '')"),
    ('N09', 'the packages joined with commas', "'; '.join(' '.join(l.split())", "', '.join(' '.join(l.split())"),
    ('N10', 'the lines added by the person who ran it dropped', "meta['versions_text'] += list(load_json(os.path.join(run_dir, 'run_extra.json')).get('versions_text', []))", "meta['versions_text'] += []"),
    ('N11', 'the GPU line of the versions without the second GPU', "' / '.join(gpu[1:3]) if len(gpu) > 2", "' / '.join(gpu[1:2]) if len(gpu) > 2"),
    ('N12', 'the budget text without the unit', "'%.0f-%.0f s (%s epochs of 10,000 simulated frames)'", "'%.0f-%.0f (%s epochs of 10,000 simulated frames)'"),
    # ---- L: comparison rows
    ('L01', 'the CNN (a) threshold not applied', "d = [c[c[:, 3] >= sel['thr']] if len(c) else c for c in C.unpack_detections(z[key + '|off'], z[key + '|flat'])]", "d = [c for c in C.unpack_detections(z[key + '|off'], z[key + '|flat'])]"),
    ('L02', 'the offset applied to the simulation-trained CNN too', "if sel['tag'].startswith('sup'):", "if True:"),
    ('L03', 'the comparison rows exchanged', "for tag in ('ISS-0', 'ISS-E', 'R0', 'R2E+c'):\n        out[tag] = C.unpack_detections(z['%s|cmp|%s|off' % (name, tag)], z['%s|cmp|%s|flat' % (name, tag)])",
     "for tag in ('ISS-0', 'ISS-E', 'R0', 'R2E+c'):\n        out[tag] = C.unpack_detections(z['%s|cmp|%s|off' % (name, 'ISS-0')], z['%s|cmp|%s|flat' % (name, 'ISS-0')])"),
    # ---- I: command line
    ('I01', 'from-raw does not write the md', "    if A.from_raw or A.ingest:", "    if A.ingest:"),
    ('I02', 'the md is written without the final newline', "return '\\n'.join(L) + '\\n'", "return '\\n'.join(L)"),
]
