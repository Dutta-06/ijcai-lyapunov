"""
Round 2 of the baselines given the ISS's front end: scoring, the selection on the development patches (the round-1 selection a candidate), the test patches and the report
======================================================================================================================================================================

Declared before any round-2 catalogue was scored in `results/astro/baselines_hp2_plan.md` (read it first; its SHA-256 is stored with the selection). `exp_astro_baselines_hp2.py` writes the catalogues of the two
round-2 families (photutils with `xy_bounds`, 6 variants; crowdsource with a 2 or 3 px high-pass, threshold 3 or 5, two pass settings, 8 variants). This script scores them with `astro_common.score` (as round 1), selects the
variant of each family with the largest pooled F at mag < 22 on the two M2 development patches among the round-2 variants AND the round-1 selected variant of the family (so that round 2 cannot return something worse on
the development patches), scores only the two selected variants on the 16 test patches and reports them beside the round-1 selections, the old tuned baselines and the tuned ISS rows of the frozen results (READ, never written).
Round 1 (`exp_astro_hp.py`, `baselines_hp_raw.npz`, `baseline_catalogs_hp.npz`) is imported / read only.

    python experiments/exp_astro_hp2.py --dev                         score the development catalogues, select (refuses if a selection exists)
    python experiments/exp_astro_hp2.py --test [--sets m2 m53 m92 m2r1 ...] [--flag FILE]   the two selected variants on the test patches (the repeat frames only if FILE exists)
    python experiments/exp_astro_hp2.py --from-raw                    regenerates results/astro/baselines_hp2.md from baselines_hp2_raw.npz (and the round-1 raw and the frozen raws, read only)

Raw keys: 'DEV|x<x0>_y<y0>|<tag>|<field>' (round-2 tags and the round-1 selected tags), 'SEL|photutils', 'SEL|crowdsource', 'SEL|r1_photutils', 'SEL|r1_crowdsource', 'SEL|plan_sha256', 'SEL|dev_hash',
'T|<set>|x<x0>_y<y0>|<row>|<field>' with row in phu2, cs2 ('|name' = the tag).
"""
import argparse
import os
import sys

os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np

HERE = os.path.abspath(os.path.dirname(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)
import astro_common as A                  # noqa: E402
import astro_repeat as R                  # noqa: E402
import exp_astro_baselines_hp as G        # noqa: E402  (round 1, read only)
import exp_astro_baselines_hp2 as G2      # noqa: E402
import exp_astro_hp as E1                 # noqa: E402  (round 1, read only)
import exp_astro_m2 as X                  # noqa: E402

ASTRO = os.path.join(REPO, 'results', 'astro')
OUT = os.environ.get('HP2_OUT', ASTRO)
RAW2, MD2 = os.path.join(OUT, 'baselines_hp2_raw.npz'), os.path.join(OUT, 'baselines_hp2.md')
CATS2 = os.environ.get('HP2_CATALOGS', os.path.join(ASTRO, 'baseline_catalogs_hp2.npz'))
CATS1 = os.environ.get('HP2_CATALOGS1', os.path.join(ASTRO, 'baseline_catalogs_hp.npz'))
RAW1 = os.environ.get('HP2_RAW1', os.path.join(ASTRO, 'baselines_hp_raw.npz'))
PLAN2 = G2.PLAN2
NEW2_ROWS = ('phu2', 'cs2')
ORDER2 = ('iss', 'iss_adq', 'phu_old_t', 'cs_old_t', 'phu1', 'cs1') + NEW2_ROWS


# ----------------------------------------------------------------------------------------------------------------------------------
# the selection
# ----------------------------------------------------------------------------------------------------------------------------------
def round1_selection(raw1=None):
    """The round-1 selection {family: tag}, verified (it is what the round-1 rule gives, the round-1 plan is unchanged)."""
    return E1.verify_selection(raw1 if raw1 is not None else X.load_raw(RAW1))


def candidates(sel1):
    """{family: [tags]}: the round-2 variants and the round-1 selected variant of the family."""
    return {fam: G2.FAMILIES2[fam] + [sel1[fam]] for fam in G2.FAMILIES2}


def choose_hp2(raw, sel1):
    """{family: tag}: among the candidates of the family the one with the largest pooled F at mag < 22 on the two development patches, compared after rounding to 6 decimals, a tie to the larger tag (the round-1
    rule, `exp_astro_hp.choose_hp`); a candidate needs a score on both patches. The round-1 tags sort above the round-2 tags, so the round-1 variant keeps a tie."""
    out = {}
    for fam, tags in candidates(sel1).items():
        elig = [t for t in tags if all('DEV|x%d_y%d|%s|tot' % (p[0], p[1], t) in raw for p in X.DEV)]
        if elig:
            out[fam] = max(elig, key=lambda t: (round(E1.dev_F(raw, t), 6), t))
    return out


def verify_selection2(raw, raw1=None, plan=None):
    """The stored round-2 selection {family: tag} after checking that the round-1 selection it used is the verified round-1 selection, that the stored selection is what `choose_hp2` gives on the stored development
    scores, that these scores are unchanged and that the round-2 plan is unchanged. Raises RuntimeError otherwise."""
    plan = plan or PLAN2
    keys = ('SEL|photutils', 'SEL|crowdsource', 'SEL|r1_photutils', 'SEL|r1_crowdsource', 'SEL|plan_sha256', 'SEL|dev_hash')
    if any(k not in raw for k in keys):
        raise RuntimeError('no selection in the raw: run `exp_astro_hp2.py --dev` first')
    sel1 = round1_selection(raw1)
    if {f: str(raw['SEL|r1_' + f][0]) for f in G2.FAMILIES2} != sel1:
        raise RuntimeError('the round-1 selection is not the one this selection was made with')
    stored = {f: str(raw['SEL|' + f][0]) for f in G2.FAMILIES2}
    if stored != choose_hp2(raw, sel1):
        raise RuntimeError('the stored selection %s is not what the rule gives on the stored development scores: %s' % (stored, choose_hp2(raw, sel1)))
    if str(raw['SEL|dev_hash'][0]) != E1.dev_hash(raw):
        raise RuntimeError('the development scores changed after the selection')
    if str(raw['SEL|plan_sha256'][0]) != E1.sha256_file(plan):
        raise RuntimeError('the plan file changed after the selection')
    return stored


# ----------------------------------------------------------------------------------------------------------------------------------
# scoring
# ----------------------------------------------------------------------------------------------------------------------------------
def part_dev2(cats2=CATS2, cats1=CATS1, raw1_path=RAW1, out_path=RAW2):
    """The round-2 variants and the round-1 selected variants on the two development patches, scored against their truth, and the selection (refuses if a selection exists in the raw)."""
    cur = X.load_raw(out_path) if os.path.exists(out_path) else {}
    if 'SEL|photutils' in cur:
        raise RuntimeError('a selection exists in %s: it is never made twice' % out_path)
    raw1 = X.load_raw(raw1_path)
    sel1 = round1_selection(raw1)
    truth = X.truth_arrays()
    b2, b1 = np.load(cats2), np.load(cats1)
    arrays = {}
    for p in X.DEV:
        t = A.truth_for_patch(truth, p[0], p[1])
        for b, tags in ((b2, G2.FAMILIES2['photutils'] + G2.FAMILIES2['crowdsource']), (b1, [sel1['photutils'], sel1['crowdsource']])):
            for tag in tags:
                k = 'x%d_y%d|%s|' % (p[0], p[1], tag)
                if k + 'x' in b.files:
                    arrays.update(X._pack(k[:-1], E1._score(b, p[0], p[1], tag, t)))
                    if k + 'secs' in b.files:
                        arrays[k + 'secs'] = np.asarray(b[k + 'secs'])
    for fam, tag in sel1.items():                                        # the round-1 selected variants are scored again: they must reproduce the round-1 development totals
        for p in X.DEV:
            k = 'x%d_y%d|%s|tot' % (p[0], p[1], tag)
            if not np.array_equal(arrays[k], raw1['DEV|' + k]):
                raise RuntimeError('the round-1 selected variant %s does not reproduce its round-1 development totals on %s' % (tag, p))
    X.save_part(out_path, 'DEV', arrays)
    raw = X.load_raw(out_path)
    sel = choose_hp2(raw, sel1)
    if set(sel) != set(G2.FAMILIES2):
        raise RuntimeError('a family has no candidate with a score on both development patches')
    X.save_part(out_path, 'SEL', dict(photutils=np.array([sel['photutils']]), crowdsource=np.array([sel['crowdsource']]), r1_photutils=np.array([sel1['photutils']]),
                                       r1_crowdsource=np.array([sel1['crowdsource']]), plan_sha256=np.array([E1.sha256_file(PLAN2)]), dev_hash=np.array([E1.dev_hash(raw)])))
    return sel


def part_test2(sets, cats2=CATS2, out_path=RAW2, flag=None, raw1=None):
    """The two selected variants on the patches of `sets`, scored against the truth of the frame; a patch whose result exists is not scored again. The selection is verified first; the repeat frames need the flag file."""
    raw = X.load_raw(out_path) if os.path.exists(out_path) else {}
    sel = verify_selection2(raw, raw1)
    tags = {'phu2': sel['photutils'], 'cs2': sel['crowdsource']}
    b = np.load(cats2)
    for name in sets:
        if name not in G.TEST_SETS:
            raise RuntimeError('unknown set %s' % name)
        if name in G.REPEAT_SETS and not (flag and os.path.exists(flag)):
            raise RuntimeError('the repeat frames are scored only after the flag file of the ISS repeat run exists')
        ds, plist = G.SETS[name]
        truth = R.truth_arrays(name) if name in G.REPEAT_SETS else X.truth_arrays(ds)
        for p in plist:
            pre = 'T|%s|x%d_y%d' % (name, p[0], p[1])
            if pre + '|phu2|tot' in raw:
                print('T %s %s: a result exists, not scored again' % (name, p), flush=True)
                continue
            t = A.truth_for_patch(truth, p[0], p[1])
            arrays = {}
            for row, tag in tags.items():
                arrays.update(X._pack(row, E1._score(b, p[0], p[1], tag, t)))
                arrays[row + '|name'] = np.array([tag])
                arrays[row + '|secs'] = np.asarray(b['x%d_y%d|%s|secs' % (p[0], p[1], tag)])
            X.save_part(out_path, pre, arrays)
            print('T %s %s scored: %s' % (name, p, ', '.join('%s=%s' % kv for kv in tags.items())), flush=True)


# ----------------------------------------------------------------------------------------------------------------------------------
# the report
# ----------------------------------------------------------------------------------------------------------------------------------
def describe2(tag):
    """The settings of a variant (round 2 or the round-1 selected ones) as printed in the development table."""
    for d in G2.PHOT2:
        if d['tag'] == tag:
            return '%s, xy_bounds %g px' % (E1.describe(d['parent']), d['bounds'])
    for d in G2.CS2:
        if d['tag'] == tag:
            p = G2.PASSES[d['passes']]
            return 'crowdsource, image and stamp high-passed %g px, the ISS weights, threshold %g, passes %s (miniter %d, maxiter %d, fewstars %d)' % (d['sb'], d['thr'], d['passes'], p['miniter'], p['maxiter'], p['fewstars'])
    return E1.describe(tag)


LABELS2 = {'phu1': 'photutils, round 1 selection', 'cs1': 'crowdsource, round 1 selection', 'phu2': 'photutils, round 2 selection', 'cs2': 'crowdsource, round 2 selection'}
SHORT2 = {'iss': 'tuned ISS', 'iss_adq': 'tuned ISS (kernel from the stamp)', 'phu_old_t': 'photutils old tuned', 'cs_old_t': 'crowdsource old tuned', 'phu1': 'photutils round 1', 'cs1': 'crowdsource round 1',
          'phu2': 'photutils round 2', 'cs2': 'crowdsource round 2'}
COMPARISONS2_PATCH = (('iss', 'phu2'), ('iss', 'cs2'), ('iss_adq', 'phu2'), ('iss_adq', 'cs2'), ('phu2', 'phu1'), ('cs2', 'cs1'), ('phu2', 'phu_old_t'), ('cs2', 'cs_old_t'))
COMPARISONS2_POOLED = COMPARISONS2_PATCH + (('iss', 'phu1'), ('iss', 'cs1'), ('iss_adq', 'phu1'), ('iss_adq', 'cs1'))


def patch_rows2(raw, raw1, name, pxy, frozen):
    """[(row, label, tot, blk22, offset)] of a test patch: the tuned ISS rows and the old tuned baselines (frozen results), the round-1 selections (round-1 raw), the round-2 selections (this raw), in the order of ORDER2."""
    rows = [(r, lab, s[p + '|tot'], s[p + '|blk22'], s[p + '|offset']) for r, lab, s, p in E1.frozen_rows(name, pxy, frozen) if r in ('iss', 'iss_adq', 'phu_old_t', 'cs_old_t')]
    pre1 = 'T|%s|x%d_y%d|' % (name, pxy[0], pxy[1])
    for row, src in (('phu1', 'phu_sel'), ('cs1', 'cs_sel')):
        if raw1 is not None and pre1 + src + '|tot' in raw1:
            rows.append((row, '%s [%s]' % (LABELS2[row], E1.shown(raw1[pre1 + src + '|name'][0])), raw1[pre1 + src + '|tot'], raw1[pre1 + src + '|blk22'], raw1[pre1 + src + '|offset']))
    for row in NEW2_ROWS:
        if pre1 + row + '|tot' in raw:
            rows.append((row, '%s [%s]' % (LABELS2[row], E1.shown(raw[pre1 + row + '|name'][0])), raw[pre1 + row + '|tot'], raw[pre1 + row + '|blk22'], raw[pre1 + row + '|offset']))
    return rows


def group_data2(raw, raw1, name_list, frozen):
    """The patches of `name_list` = [(set, patches)] with round-2 rows in the raw, pooled: dict(done, rows (present on every one of them, in the order of ORDER2), label, blk, tot, sizes)."""
    done, per = [], {}
    for name, plist in name_list:
        for p in plist:
            if 'T|%s|x%d_y%d|phu2|tot' % (name, p[0], p[1]) not in raw:
                continue
            done.append((name, tuple(p)))
            for row, label, tot, blk, off in patch_rows2(raw, raw1, name, p, frozen):
                per.setdefault(row, []).append((label, tot, blk))
    rows = [r for r in ORDER2 if r in per and len(per[r]) == len(done)]
    return dict(done=done, rows=rows, label={r: per[r][0][0] for r in rows}, blk={r: np.concatenate([b for l, t, b in per[r]]) for r in rows},
                tot={r: sum(t for l, t, b in per[r]) for r in rows}, sizes=[len(b) for l, t, b in per[rows[0]]] if rows else [])


def diff_table2(pairs, g, idx):
    """The paired differences first minus second for the pairs of rows present in the group: (table or '', {pair: (d, interval) of F})."""
    rows, dfs = [], {}
    for a, b in pairs:
        if a in g['blk'] and b in g['blk']:
            dc, dp, df = X.boot_diff(g['blk'][a], g['blk'][b], idx)
            rows.append(['%s minus %s' % (SHORT2[a], SHORT2[b]), X.fmt_iv(dc[0], dc[1]), X.fmt_iv(dp[0], dp[1]), X.fmt_iv(df[0], df[1]), X.verdict(*df)])
            dfs[(a, b)] = df
    return (E1.table(['comparison', 'd completeness', 'd purity', 'd F', 'verdict on F'], rows) if rows else ''), dfs


def dev_table2(raw, raw1):
    """The round-2 development table: every candidate (round-2 variants, then the round-1 selected variants) with its pooled counts on the two development patches, the selected ones marked, the old baselines of round 1's OLD keys."""
    sel = {f: str(raw['SEL|' + f][0]) for f in G2.FAMILIES2} if 'SEL|photutils' in raw else {}
    sel1 = {f: str(raw['SEL|r1_' + f][0]) for f in G2.FAMILIES2} if 'SEL|r1_photutils' in raw else {}
    rows = []
    for fam, tags in G2.FAMILIES2.items():
        for tag in tags + ([sel1[fam]] if fam in sel1 else []):
            if not all('DEV|x%d_y%d|%s|tot' % (p[0], p[1], tag) in raw for p in X.DEV):
                continue
            tot = sum(raw['DEV|x%d_y%d|%s|tot' % (p[0], p[1], tag)] for p in X.DEV)
            ff, c, pu = X.f_measure(tot)
            secs = sum(float(raw['DEV|x%d_y%d|%s|secs' % (p[0], p[1], tag)][0]) for p in X.DEV if 'DEV|x%d_y%d|%s|secs' % (p[0], p[1], tag) in raw)
            mark = '**' if sel.get(fam) == tag else ''
            rows.append([fam, '%s%s%s%s' % (mark, tag, mark, ' (round 1 selection)' if sel1.get(fam) == tag else ''), describe2(tag), '%.3f' % c, '%.3f' % pu, '%s%.3f%s' % (mark, ff, mark), '%d' % tot[2],
                         '%d' % (tot[6] - tot[5]), '%.0f' % secs])
    if raw1 is not None:
        for label in ('photutils (default)', 'photutils (tuned)', 'crowdsource (default)', 'crowdsource (tuned)'):
            if all('OLD|x%d_y%d|%s|tot' % (p[0], p[1], label) in raw1 for p in X.DEV):
                tot = sum(raw1['OLD|x%d_y%d|%s|tot' % (p[0], p[1], label)] for p in X.DEV)
                ff, c, pu = X.f_measure(tot)
                rows.append(['old family', label, E1.shown(raw1['OLD|%s|name' % label][0]), '%.3f' % c, '%.3f' % pu, '%.3f' % ff, '%d' % tot[2], '%d' % (tot[6] - tot[5]), '-'])
    return E1.table(['family', 'variant', 'settings', 'completeness (mag < 22)', 'purity (mag < 22)', 'F', 'detections < 22 mag', 'unmatched detections (all magnitudes)', 'CPU seconds (two patches)'], rows)


def patch_verdicts(raw, raw1, frozen):
    """The per-patch verdicts of the paired F difference (95 % interval of the patch's own 25 blocks, the round-1 bootstrap) for the test patches in the raw: {(isskind, baseline): {set: [verdicts]}} with isskind 'iso'
    (the isotropic-kernel row as frozen) or 'stamp' (the row with the kernel from the stamp where the patch has one: M53 and M92 the post hoc row, the repeat patches the `adq` row; M2 the isotropic row, which is that kernel there),
    baseline in phu2, cs2, phu1, cs1."""
    out = {}
    for title, name in E1.GROUPS:
        for pxy in G.SETS[name][1]:
            g = group_data2(raw, raw1, [(name, [pxy])], frozen)
            if not g['done']:
                continue
            idx = X.boot_idx_strat(g['sizes'])
            for kind, row in (('iso', 'iss'), ('stamp', 'iss_adq' if 'iss_adq' in g['blk'] else 'iss')):
                for base in ('phu2', 'cs2', 'phu1', 'cs1'):
                    if row in g['blk'] and base in g['blk']:
                        out.setdefault((kind, base), {}).setdefault(name, []).append(X.verdict(*X.boot_diff(g['blk'][row], g['blk'][base], idx)[2]))
    return out


def verdict_table(vd):
    """The table of the per-patch verdict counts: wins / ties / losses of the ISS over the 16 patches and per group."""
    names = [n for _, n in E1.GROUPS]
    rows = []
    for (kind, base), per in sorted(vd.items(), key=lambda kv: (kv[0][0], ('phu2', 'cs2', 'phu1', 'cs1').index(kv[0][1]))):
        allv = [v for n in names for v in per.get(n, [])]
        cells = ['%d / %d / %d' % tuple(sum(1 for v in per.get(n, []) if v == w) for w in ('win', 'tie', 'loss')) if n in per else '-' for n in names]
        rows.append(['%s minus %s' % ('tuned ISS, isotropic kernel as frozen' if kind == 'iso' else 'tuned ISS, kernel from the stamp (M2: the isotropic one)', SHORT2[base]),
                     '%d / %d / %d (%d patches)' % (tuple(sum(1 for v in allv if v == w) for w in ('win', 'tie', 'loss')) + (len(allv),))] + cells)
    return E1.table(['comparison (ISS better = win)', 'wins / ties / losses'] + ['%s' % n for n in names], rows)


def summary_tables2(summary):
    """Two tables, one row per pooled group: the F with its 95 % interval of every row, and the pooled paired differences with their verdicts."""
    t1, t2 = [], []
    for title, g, fs, dfs in summary:
        def cell(r):
            if r not in fs:
                return '-'
            ff, ci = fs[r]
            return '%.3f [%.3f, %.3f] (%d)' % (ff, ci[0], ci[1], g['tot'][r][6] - g['tot'][r][5])

        def vd(a, b):
            return '-' if (a, b) not in dfs else '%+.3f %s' % (dfs[(a, b)][0], X.verdict(*dfs[(a, b)]))
        t1.append([title + ' (%d patches)' % len(g['done'])] + [cell(r) for r in ORDER2])
        t2.append([title] + [vd(a, b) for a, b in (('iss', 'phu2'), ('iss', 'cs2'), ('iss_adq', 'phu2'), ('iss_adq', 'cs2'), ('phu2', 'phu1'), ('cs2', 'cs1'), ('phu2', 'phu_old_t'), ('cs2', 'cs_old_t'))])
    return (E1.table(['group'] + [SHORT2[r] for r in ORDER2], t1),
            E1.table(['group', 'ISS minus photutils r2', 'ISS minus crowdsource r2', 'ISS (stamp kernel) minus photutils r2', 'ISS (stamp kernel) minus crowdsource r2', 'photutils r2 minus r1',
                      'crowdsource r2 minus r1', 'photutils r2 minus old tuned', 'crowdsource r2 minus old tuned'], t2))


def pooled_block2(heading, g, raw_summary, title):
    idx = X.boot_idx_strat(g['sizes'])
    fs = {r: X.boot_cpf(g['blk'][r], idx)[2] for r in g['rows']}
    dtab, dfs = diff_table2(COMPARISONS2_POOLED, g, idx)
    raw_summary.append((title, g, fs, dfs))
    return [heading + '\n', 'The blocks of the %d patches pooled (%d blocks); counts are sums; the bootstrap resamples the blocks WITHIN each patch (stratified), the same resamples for every row. The patches '
            'are windows of frames of one cluster, not independent trials.\n' % (len(g['done']), sum(g['sizes'])), E1.method_table(g, idx), '', 'Paired differences on the pooled blocks (first minus second):\n', dtab, '']


def report2(raw, frozen=None, raw1=None):
    """The markdown report of round 2 from its raw, the round-1 raw and the frozen results (all read only)."""
    frozen = frozen if frozen is not None else E1.load_frozen()
    raw1 = raw1 if raw1 is not None else X.load_raw(RAW1)
    L = ['# The baselines given the ISS\'s front end, round 2: photutils with `xy_bounds`, crowdsource with more passes and a 2 px high-pass', '',
         'Generated by `experiments/exp_astro_hp2.py` from `baselines_hp2_raw.npz` (and the round-1 `baselines_hp_raw.npz`, the frozen `m2_raw.npz`, `psf_raw.npz`, `repeat_raw.npz`, all read only); every number below is printed by the script. '
         'Round 2 was declared in `baselines_hp2_plan.md` before any round-2 catalogue was scored (SHA-256 of the plan stored with the selection: %s). The family was written after the referee\'s development probes (plan, section 0): '
         'the development F of the winners is optimistic, and the round-1 test results were known when round 2 was run, so the test scores are results of a declared procedure, not a blind confirmation. '
         'Scoring and intervals as round 1 (match at 0.5 px and 0.5 mag, F at mag < 22, 95 %% percentile intervals of a bootstrap over the 20 x 20 px blocks, %d resamples, seed 0, the same resampled blocks for every row of a table; '
         'the pooled tables resample the blocks within each patch). The selection reads the two M2 development patches only; the tuned ISS rows are READ from the frozen results and influenced nothing. '
         'On the repeat frames these baselines are added later, as in round 1.' % (str(raw['SEL|plan_sha256'][0]) if 'SEL|plan_sha256' in raw else 'none yet', X.BOOT), '']
    shas = frozen.get('sha', {})
    L += ['Frozen results read: %s.' % '; '.join('`%s` %s' % (n, (s[:16] + '...') if s else 'absent') for n, s in (('m2_raw.npz', shas.get('m2')), ('psf_raw.npz', shas.get('psf')), ('repeat_raw.npz', shas.get('repeat')))), '']
    if any(k.startswith('DEV|') for k in raw):
        L += ['## Development patches %s: every candidate\n' % ', '.join('(%d, %d)' % p for p in X.DEV),
              'Pooled counts of the two M2 development patches (the only data the selection sees); bold: the candidate of each family with the largest F (a tie to the larger tag; the round-1 tags sort above the round-2 tags). '
              'The round-1 selected variants are candidates (scored again from the round-1 catalogues); the old baselines are for reference.\n', dev_table2(raw, raw1), '']
        if 'SEL|photutils' in raw:
            L += ['Selected: photutils `%s` (round 1: `%s`), crowdsource `%s` (round 1: `%s`).' % (raw['SEL|photutils'][0], raw['SEL|r1_photutils'][0], raw['SEL|crowdsource'][0], raw['SEL|r1_crowdsource'][0]), '']
            k = 'DEV|x%d_y%d|cs2_hp3_t5_pd|tot'
            if all(k % p in raw and 'DEV|x%d_y%d|%s|tot' % (p[0], p[1], raw['SEL|r1_crowdsource'][0]) in raw for p in X.DEV):
                same = all(np.array_equal(raw[k % p], raw['DEV|x%d_y%d|%s|tot' % (p[0], p[1], raw['SEL|r1_crowdsource'][0])]) for p in X.DEV)
                L += ['Check: the round-2 runner reproduces the round-1 selected crowdsource variant (`cs2_hp3_t5_pd` against `%s`, development totals of both patches): %s.' % (raw['SEL|r1_crowdsource'][0], 'identical' if same else 'DIFFERENT'), '']
    summary, body = [], []
    for title, name in E1.GROUPS:
        g = group_data2(raw, raw1, [(name, G.SETS[name][1])], frozen)
        if not g['done']:
            continue
        body += ['## %s\n' % title]
        for _, p in g['done']:
            gp = group_data2(raw, raw1, [(name, [p])], frozen)
            idx = X.boot_idx_strat(gp['sizes'])
            body += ['### Patch (%d, %d)\n' % p, 'Truth sources brighter than 22 mag in the scored area: %d.\n' % gp['tot']['phu2'][0], E1.method_table(gp, idx), '',
                     'Paired differences (first minus second, same blocks), 95 % intervals:\n', diff_table2(COMPARISONS2_PATCH, gp, idx)[0], '']
        if len(g['done']) > 1:
            body += pooled_block2('### All %d test patches together' % len(g['done']), g, summary, title.split(',')[0])
        else:
            idx = X.boot_idx_strat(g['sizes'])
            summary.append((title.split(',')[0], g, {r: X.boot_cpf(g['blk'][r], idx)[2] for r in g['rows']}, diff_table2(COMPARISONS2_POOLED, g, idx)[1]))
    gr = group_data2(raw, raw1, [(n, G.SETS[n][1]) for n in sorted(R.REPEAT)], frozen)
    if len(gr['done']) > 1:
        body += pooled_block2('## All the repeat patches together', gr, summary, 'All the repeat patches')
    if summary:
        t1, t2 = summary_tables2(summary)
        L += ['## Summary: pooled F at mag < 22 by cluster\n', 'F with its 95 % interval and, in brackets, the unmatched detections of all magnitudes. "Tuned ISS" is the isotropic-kernel row as frozen; "tuned ISS (kernel from the stamp)" '
              'is the post hoc row of `psf_raw.npz` on M53 and M92 and the repeat plan\'s `adq` row on the repeat patches (equal to the isotropic row where the rule keeps the isotropic kernel); it does not exist for M2.\n', t1, '',
              'Pooled paired differences dF (first minus second) with the verdict of the 95 % interval.\n', t2, '']
        vd = patch_verdicts(raw, raw1, frozen)
        if vd:
            L += ['## Verdict counts over the test patches\n', 'For each test patch in the raw the verdict of the paired difference of F on the patch\'s own 25 blocks (95 % interval; win: the ISS is better), counted by the script. Two ISS rows: '
                  'the isotropic-kernel row as frozen, and the row with the kernel from the stamp (M2 patches: the isotropic row; M53 and M92: the post hoc stamp-kernel row; repeat patches: the `adq` row of `repeat_raw.npz`). '
                  'The columns are the groups (wins / ties / losses).\n', verdict_table(vd), '']
    return '\n'.join(L + body) + '\n'


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--dev', action='store_true')
    ap.add_argument('--test', action='store_true')
    ap.add_argument('--sets', nargs='*', default=['m2', 'm53', 'm92'])
    ap.add_argument('--flag')
    ap.add_argument('--from-raw', action='store_true')
    A_ = ap.parse_args()
    if A_.dev:
        print(part_dev2())
    elif A_.test:
        part_test2(A_.sets, flag=A_.flag)
    elif A_.from_raw:
        with open(MD2, 'w') as f:
            f.write(report2(X.load_raw(RAW2)))
        print('wrote', MD2)
    else:
        ap.error('give --dev, --test or --from-raw')


if __name__ == '__main__':
    main()
