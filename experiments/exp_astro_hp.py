"""
The baselines given the ISS's front end and a wider pipeline search: scoring, the selection on the development patches, the test patches and the report
==================================================================================================================================================

Declared before any catalogue was scored in `results/astro/baselines_hp_plan.md` (read it first; its SHA-256 is stored with the selection). `exp_astro_baselines_hp.py` writes the catalogues of two families
(photutils on the ISS's high-passed patch, 8 variants; crowdsource with its detection threshold and the same front end, 8 variants; the current defaults are the references); this script scores them with
`astro_common.score` (match 0.5 px and 0.5 mag, one to one; F-measure of completeness and purity at mag < 22, as `exp_astro_m2`), selects the variant of each family with the largest pooled F on the two M2
development patches (`exp_astro_m2.choose_baselines`' rule: counts of the two patches added, F rounded to 6 decimals, a tie to the larger tag), scores only the selected variants and the references on the test
patches (the four M2 test patches, the four M53 and two M92 patches, the six repeat patches), and reports them beside the old tuned baselines and the tuned ISS rows of the frozen results, which are READ (never
written) from `m2_raw.npz`, `psf_raw.npz` and `repeat_raw.npz`: they influence nothing that is declared or selected here.

    python experiments/exp_astro_hp.py --dev                       score the development catalogues, select, write the DEV| and SEL| keys (refuses if a selection exists)
    python experiments/exp_astro_hp.py --old                       the old baselines on the development patches, for the development table (after the selection; no part in it)
    python experiments/exp_astro_hp.py --test [--sets m2 m53 m92 m2r1 ...] [--flag FILE]   the selected variants and the references on the test patches (the repeat frames only if FILE exists)
    python experiments/exp_astro_hp.py --from-raw                  regenerates results/astro/baselines_hp.md from baselines_hp_raw.npz (and the frozen raws, read only)

Raw keys: 'DEV|x<x0>_y<y0>|<tag>|<field>' (the fields of `exp_astro_m2._pack`, '|secs'), 'OLD|x..|<label>|<field>', 'SEL|photutils', 'SEL|crowdsource', 'SEL|plan_sha256', 'SEL|dev_hash',
'T|<set>|x<x0>_y<y0>|<row>|<field>' with row in phu_sel, phu_ref, cs_sel, cs_hp, cs_default ('|name' = the tag).
"""
import argparse
import hashlib
import json
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
import exp_astro_baselines_hp as G        # noqa: E402  (numpy only at import: the registries)
import exp_astro_m2 as X                  # noqa: E402

ASTRO = os.path.join(REPO, 'results', 'astro')
OUT = os.environ.get('HP_OUT', ASTRO)
RAW, MD = os.path.join(OUT, 'baselines_hp_raw.npz'), os.path.join(OUT, 'baselines_hp.md')
CATS = os.environ.get('HP_CATALOGS', os.path.join(ASTRO, 'baseline_catalogs_hp.npz'))
PLAN = G.PLAN
M2RAW, PSFRAW, REPEATRAW = (os.environ.get(e, os.path.join(ASTRO, f)) for e, f in (('HP_M2_RAW', 'm2_raw.npz'), ('HP_PSF_RAW', 'psf_raw.npz'), ('HP_REPEAT_RAW', 'repeat_raw.npz')))
FREEZE = os.environ.get('HP_FREEZE', os.path.join(ASTRO, 'freeze.json'))
NEW_ROWS = ('phu_sel', 'phu_ref', 'cs_sel', 'cs_hp', 'cs_default')            # the catalogues scored on every test patch
GROUPS = (('M2 (NGC 7089), SDSS run 2583', 'm2'),) + tuple((X.SECOND[k]['title'], k) for k in ('m53', 'm92')) + tuple((R.REPEAT[n]['title'], n) for n in sorted(R.REPEAT))


# ----------------------------------------------------------------------------------------------------------------------------------
# hashes, the selection
# ----------------------------------------------------------------------------------------------------------------------------------
def sha256_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def dev_hash(raw):
    """The hash of the development scores of a raw (every key 'DEV|...', in sorted order)."""
    h = hashlib.sha256()
    for k in sorted(k for k in raw if k.startswith('DEV|')):
        h.update(k.encode())
        h.update(np.ascontiguousarray(raw[k]).tobytes())
    return h.hexdigest()


def dev_F(raw, tag):
    """Pooled F at mag < 22 of a variant on the two development patches (the counts added)."""
    return X.f_measure(sum(raw['DEV|x%d_y%d|%s|tot' % (p[0], p[1], tag)] for p in X.DEV))[0]


def choose_hp(raw):
    """{family: tag}: the variant of each family with the largest pooled F at mag < 22 on the two development patches, compared after rounding to 6 decimals, a tie to the larger tag
    (`exp_astro_m2.choose_baselines`); a variant needs a score on both patches to be eligible."""
    out = {}
    for fam, tags in G.FAMILIES.items():
        elig = [t for t in tags if all('DEV|x%d_y%d|%s|tot' % (p[0], p[1], t) in raw for p in X.DEV)]
        if elig:
            out[fam] = max(elig, key=lambda t: (round(dev_F(raw, t), 6), t))
    return out


def verify_selection(raw, plan=None):
    """The stored selection {family: tag}, after checking that it is what `choose_hp` gives on the stored development scores, that these scores are the ones it was made from and that the plan file is
    unchanged. Raises RuntimeError otherwise (a test number is never produced from another selection)."""
    plan = plan or PLAN
    if 'SEL|photutils' not in raw or 'SEL|crowdsource' not in raw:
        raise RuntimeError('no selection in the raw: run `exp_astro_hp.py --dev` first')
    stored = {f: str(raw['SEL|' + f][0]) for f in G.FAMILIES}
    if stored != choose_hp(raw):
        raise RuntimeError('the stored selection %s is not what the rule gives on the stored development scores: %s' % (stored, choose_hp(raw)))
    if str(raw['SEL|dev_hash'][0]) != dev_hash(raw):
        raise RuntimeError('the development scores changed after the selection')
    if str(raw['SEL|plan_sha256'][0]) != sha256_file(plan):
        raise RuntimeError('the plan file changed after the selection')
    return stored


def test_tags(sel):
    """{row: tag} of the catalogues scored on every test patch: the selected variants and the references."""
    return {'phu_sel': sel['photutils'], 'phu_ref': G.PHOT_REF, 'cs_sel': sel['crowdsource'], 'cs_hp': 'cs_hp5_t5', 'cs_default': 'cs_default'}


# ----------------------------------------------------------------------------------------------------------------------------------
# scoring
# ----------------------------------------------------------------------------------------------------------------------------------
def _score(b, x0, y0, tag, t):
    k = 'x%d_y%d|%s|' % (x0, y0, tag)
    return X.score_dets((b[k + 'x'], b[k + 'y'], b[k + 'flux']), t)


def part_dev(cats=CATS, out_path=RAW):
    """Every variant on the two development patches, scored against their truth, and the selection (refuses to run if a selection exists in the raw)."""
    cur = X.load_raw(out_path) if os.path.exists(out_path) else {}
    if 'SEL|photutils' in cur:
        raise RuntimeError('a selection exists in %s: it is never made twice' % out_path)
    truth = X.truth_arrays()
    b = np.load(cats)
    arrays = {}
    for p in X.DEV:
        t = A.truth_for_patch(truth, p[0], p[1])
        for tags in G.FAMILIES.values():
            for tag in tags:
                k = 'x%d_y%d|%s|' % (p[0], p[1], tag)
                if k + 'x' in b.files:
                    arrays.update(X._pack(k[:-1], _score(b, p[0], p[1], tag, t)))
                    if k + 'secs' in b.files:
                        arrays[k + 'secs'] = np.asarray(b[k + 'secs'])
    X.save_part(out_path, 'DEV', arrays)
    raw = X.load_raw(out_path)
    sel = choose_hp(raw)
    if set(sel) != set(G.FAMILIES):
        raise RuntimeError('a family has no variant with a score on both development patches')
    X.save_part(out_path, 'SEL', dict(photutils=np.array([sel['photutils']]), crowdsource=np.array([sel['crowdsource']]), plan_sha256=np.array([sha256_file(PLAN)]),
                                       dev_hash=np.array([dev_hash(raw)])))
    return sel


def part_old(out_path=RAW, freeze=FREEZE, base=None):
    """The old baselines on the development patches (the default and the tuned variant of each family, the tuned names from freeze.json), for the development table. Run after the selection: they take no part in it."""
    with open(freeze) as f:
        tuned = json.load(f)['tuned_baselines']
    truth = X.truth_arrays()
    arrays = {}
    for p in X.DEV:
        bs = X.baseline_scores(p, A.truth_for_patch(truth, p[0], p[1]), base)
        for label, name in (('photutils (default)', 'photutils'), ('crowdsource (default)', X.cs_name((0.0, 0.0, 0.0, np.inf))), ('photutils (tuned)', tuned['photutils']),
                            ('crowdsource (tuned)', tuned['crowdsource'])):
            arrays.update(X._pack('x%d_y%d|%s' % (p[0], p[1], label), bs[name]))
            arrays['%s|name' % label] = np.array([name])
    X.save_part(out_path, 'OLD', arrays)


def part_test(sets, cats=CATS, out_path=RAW, flag=None):
    """The selected variants and the references on the patches of `sets` (names of `exp_astro_baselines_hp.SETS`), scored against the truth of the frame; a patch whose result exists is not scored again. The
    selection is verified first; the repeat frames need the flag file."""
    raw = X.load_raw(out_path) if os.path.exists(out_path) else {}
    tags = test_tags(verify_selection(raw))
    b = np.load(cats)
    for name in sets:
        if name not in G.TEST_SETS:
            raise RuntimeError('unknown set %s' % name)
        if name in G.REPEAT_SETS and not (flag and os.path.exists(flag)):
            raise RuntimeError('the repeat frames are scored only after the flag file of the ISS repeat run exists')
        ds, plist = G.SETS[name]
        truth = R.truth_arrays(name) if name in G.REPEAT_SETS else X.truth_arrays(ds)
        for p in plist:
            pre = 'T|%s|x%d_y%d' % (name, p[0], p[1])
            if pre + '|phu_sel|tot' in raw:
                print('T %s %s: a result exists, not scored again' % (name, p), flush=True)
                continue
            t = A.truth_for_patch(truth, p[0], p[1])
            arrays = {}
            for row, tag in tags.items():
                arrays.update(X._pack(row, _score(b, p[0], p[1], tag, t)))
                arrays[row + '|name'] = np.array([tag])
                arrays[row + '|secs'] = np.asarray(b['x%d_y%d|%s|secs' % (p[0], p[1], tag)])
            X.save_part(out_path, pre, arrays)
            print('T %s %s scored: %s' % (name, p, ', '.join('%s=%s' % kv for kv in tags.items())), flush=True)


# ----------------------------------------------------------------------------------------------------------------------------------
# the frozen results (read only)
# ----------------------------------------------------------------------------------------------------------------------------------
def load_frozen(paths=None):
    """{'m2': raw, 'psf': raw, 'repeat': raw (None where the file is absent), 'sha': {name: sha256 or None}} of the frozen results, loaded read only."""
    paths = paths or dict(m2=M2RAW, psf=PSFRAW, repeat=REPEATRAW)
    out = {k: (X.load_raw(p) if os.path.exists(p) else None) for k, p in paths.items()}
    out['sha'] = {k: (sha256_file(p) if os.path.exists(p) else None) for k, p in paths.items()}
    return out


def kappa_of(src, pre):
    """The reported significance of a frozen patch (the third element of its 'cfg')."""
    return float(src[pre + 'cfg'][2])


def frozen_rows(name, pxy, frozen):
    """[(row, label, source raw, key prefix)] of the frozen rows of a test patch that exist: the tuned ISS (isotropic kernel as frozen; the kernel chosen from the stamp where a post hoc or repeat-plan row
    exists) and the old baselines (tuned and default)."""
    out = []
    if name in G.REPEAT_SETS:
        src, pre = frozen.get('repeat'), 'R|%s|x%d_y%d|' % (name, pxy[0], pxy[1])
        if src is None or pre + 'cfg' not in src:
            return out
        k = 'k%g' % kappa_of(src, pre)
        iss = [('iss', 'tuned ISS, isotropic kernel (as frozen)', src, pre + 'iso|ours|' + k), ('iss_adq', 'tuned ISS, kernel chosen from the stamp', src, pre + 'adq|ours|' + k)]
        base = pre
    else:
        src = frozen.get('m2')
        pre = 'T|x%d_y%d|' % pxy if name == 'm2' else 'T2|%s|x%d_y%d|' % (name, pxy[0], pxy[1])
        if src is None or pre + 'cfg' not in src:
            return out
        k = 'k%g' % kappa_of(src, pre)
        iss = [('iss', 'tuned ISS, isotropic kernel (as frozen)', src, pre + 'ours|' + k)]
        if name != 'm2' and frozen.get('psf') is not None:                                                   # (a patch the post hoc run did not cover has no such row: dropped below)
            iss.append(('iss_adq', 'tuned ISS, kernel chosen from the stamp (post hoc)', frozen['psf'], 'P|%s|x%d_y%d|ours|%s' % (name, pxy[0], pxy[1], k)))
        base = pre
    out += iss
    for row, label, key in (('phu_old_t', 'photutils, old tuned', 'photutils (tuned)'), ('cs_old_t', 'crowdsource, old tuned', 'crowdsource (tuned)'),
                            ('phu_old_d', 'photutils, old default', 'photutils (default)'), ('cs_old_d', 'crowdsource, old default', 'crowdsource (default)')):
        if src is not None and base + key + '|tot' in src:
            nm = src.get(base + key + '|name')
            out.append((row, label + (' [%s]' % shown(nm[0]) if nm is not None else ''), src, base + key))
    return [r for r in out if r[3] + '|tot' in r[2]]


NEW_LABELS = {'phu_sel': 'photutils on the high-passed patch, selected', 'phu_ref': 'photutils on the 5 px high-pass, default settings (reference)',
              'cs_sel': 'crowdsource, selected', 'cs_hp': 'crowdsource, default pipeline on the 5 px high-pass with the ISS weights (reference)', 'cs_default': 'crowdsource, default (re-run)'}


def patch_rows(raw, name, pxy, frozen):
    """[(row, label, tot, blk22, offset)] of a test patch in the table order: the frozen rows, then the new ones (present in the raw)."""
    rows = [(r, lab, s[p + '|tot'], s[p + '|blk22'], s[p + '|offset']) for r, lab, s, p in frozen_rows(name, pxy, frozen)]
    pre = 'T|%s|x%d_y%d|' % (name, pxy[0], pxy[1])
    for r in NEW_ROWS:
        if pre + r + '|tot' in raw:
            rows.append((r, NEW_LABELS[r] + (' [%s]' % shown(raw[pre + r + '|name'][0]) if r.endswith('sel') else ''), raw[pre + r + '|tot'], raw[pre + r + '|blk22'], raw[pre + r + '|offset']))
    return rows


# ----------------------------------------------------------------------------------------------------------------------------------
# the report
# ----------------------------------------------------------------------------------------------------------------------------------
def shown(tag):
    """A tag as it is printed in a table cell: the tags of the old baselines contain pipes ('crowdsource|snr5|ff0|qf0|spinf') that would split the cell, so they print with slashes. Only the printed
    text changes: the stored keys, the names in the raw and the freeze records are untouched."""
    return str(tag).replace('|', '/')


def table(header, rows):
    """`exp_astro_m2.table` with the guarantee that a cell stays one cell: a pipe in a cell is refused (print tags with `shown`)."""
    for r in [header] + list(rows):
        for c in r:
            if '|' in str(c):
                raise ValueError('a pipe in the table cell %r would split it' % str(c))
    return X.table(header, rows)


def describe(tag):
    for d in G.PHOT_HP:
        if d['tag'] == tag:
            return 'photutils, high-pass %g px, mode %s, %d iterations%s' % (d['sb'], d['mode'], d['iters'], ', grouper %g px, relaxed finder cuts' % d['group'] if d['group'] else '')
    for d in G.CS_HP:
        if d['tag'] == tag:
            fe = {'default': 'raw image, the default weights', 'raw': 'raw image, the ISS weights', 'hp': 'image and stamp high-passed %g px, the ISS weights' % d['sb']}[d['kind']]
            return 'crowdsource, %s, threshold %g' % (fe, d['thr'])
    raise KeyError(tag)


def purity_all(tot):
    """(matched detections over all detections of any magnitude, unmatched detections) from the totals [true22, matched22, det22, det22_matched, true_all, matched_all, det_all]."""
    return tot[5] / max(tot[6], 1), int(tot[6] - tot[5])


ORDER = ('iss', 'iss_adq', 'phu_old_t', 'cs_old_t', 'phu_old_d', 'cs_old_d') + NEW_ROWS
COMPARISONS_PATCH = (('iss', 'phu_sel'), ('iss', 'cs_sel'), ('iss_adq', 'phu_sel'), ('iss_adq', 'cs_sel'), ('phu_sel', 'phu_old_t'), ('cs_sel', 'cs_old_t'))
COMPARISONS_POOLED = (('iss', 'phu_sel'), ('iss', 'cs_sel'), ('iss', 'phu_ref'), ('iss', 'cs_hp'), ('iss_adq', 'phu_sel'), ('iss_adq', 'cs_sel'), ('iss_adq', 'phu_ref'), ('iss_adq', 'cs_hp'),
                      ('phu_sel', 'phu_old_t'), ('cs_sel', 'cs_old_t'), ('phu_ref', 'phu_old_d'), ('cs_hp', 'cs_old_d'), ('phu_sel', 'phu_ref'), ('cs_sel', 'cs_hp'))
SHORT = {'iss': 'tuned ISS', 'iss_adq': 'tuned ISS (kernel from the stamp)', 'phu_sel': 'photutils HP selected', 'cs_sel': 'crowdsource selected', 'phu_ref': 'photutils HP reference',
         'cs_hp': 'crowdsource HP reference', 'phu_old_t': 'photutils old tuned', 'cs_old_t': 'crowdsource old tuned', 'phu_old_d': 'photutils old default', 'cs_old_d': 'crowdsource old default',
         'cs_default': 'crowdsource default re-run'}


def group_data(raw, name_list, frozen):
    """The patches of `name_list` = [(set, patches)] that have new rows in the raw, pooled: dict(done = [(set, patch)], rows = the rows present on every one of them in the table order, label, blk (blocks of
    the patches stacked), tot (sums of the totals), sizes (blocks per patch, the strata of the bootstrap))."""
    done, per = [], {}
    for name, plist in name_list:
        for p in plist:
            if 'T|%s|x%d_y%d|phu_sel|tot' % (name, p[0], p[1]) not in raw:
                continue
            done.append((name, tuple(p)))
            for row, label, tot, blk, off in patch_rows(raw, name, p, frozen):
                per.setdefault(row, []).append((label, tot, blk))
    rows = [r for r in ORDER if r in per and len(per[r]) == len(done)]
    return dict(done=done, rows=rows, label={r: per[r][0][0] for r in rows}, blk={r: np.concatenate([b for l, t, b in per[r]]) for r in rows},
                tot={r: sum(t for l, t, b in per[r]) for r in rows}, sizes=[len(b) for l, t, b in per[rows[0]]] if rows else [])


def method_table(g, idx):
    """The table of a group: completeness, purity, F with 95 % intervals over the resamples idx, purity over all magnitudes and the unmatched detections."""
    out = []
    for row in g['rows']:
        tot = g['tot'][row]
        (c, cci), (pu, pci), (ff, fci) = X.boot_cpf(g['blk'][row], idx)
        pa, un = purity_all(tot)
        out.append([g['label'][row], '%d / %d' % (tot[1], tot[0]), X.fmt_iv(c, cci), '%d / %d' % (tot[3], tot[2]), X.fmt_iv(pu, pci), X.fmt_iv(ff, fci), '%.3f (%d)' % (pa, un)])
    return table(['catalogue', 'matched / truth < 22', 'completeness', 'matched / detected < 22', 'purity', 'F', 'purity over all magnitudes (unmatched detections)'], out)


def diff_table(pairs, g, idx):
    """The paired differences first minus second (same blocks, same resamples) for the pairs of rows present in the group: (table or '', {pair: (d, interval) of F})."""
    rows, dfs = [], {}
    for a, b in pairs:
        if a in g['blk'] and b in g['blk']:
            dc, dp, df = X.boot_diff(g['blk'][a], g['blk'][b], idx)
            rows.append(['%s minus %s' % (SHORT[a], SHORT[b]), X.fmt_iv(dc[0], dc[1]), X.fmt_iv(dp[0], dp[1]), X.fmt_iv(df[0], df[1]), X.verdict(*df)])
            dfs[(a, b)] = df
    return (table(['comparison', 'd completeness', 'd purity', 'd F', 'verdict on F'], rows) if rows else ''), dfs


def dev_table(raw):
    """The development table: every variant of both families with its pooled counts on the two development patches (selection marked), and the old baselines where scored."""
    sel = {f: str(raw['SEL|' + f][0]) for f in G.FAMILIES} if 'SEL|photutils' in raw else {}
    rows = []
    for fam, tags in G.FAMILIES.items():
        for tag in tags:
            if not all('DEV|x%d_y%d|%s|tot' % (p[0], p[1], tag) in raw for p in X.DEV):
                continue
            tot = sum(raw['DEV|x%d_y%d|%s|tot' % (p[0], p[1], tag)] for p in X.DEV)
            ff, c, pu = X.f_measure(tot)
            secs = sum(float(raw['DEV|x%d_y%d|%s|secs' % (p[0], p[1], tag)][0]) for p in X.DEV if 'DEV|x%d_y%d|%s|secs' % (p[0], p[1], tag) in raw)
            mark = '**' if sel.get(fam) == tag else ''
            ref = ' (reference)' if tag in (G.PHOT_REF,) + G.CS_REFS else ''
            rows.append([fam, '%s%s%s%s' % (mark, tag, mark, ref), describe(tag), '%.3f' % c, '%.3f' % pu, '%s%.3f%s' % (mark, ff, mark), '%d' % tot[2], '%d' % (tot[6] - tot[5]), '%.0f' % secs])
    for label in ('photutils (default)', 'photutils (tuned)', 'crowdsource (default)', 'crowdsource (tuned)'):
        if all('OLD|x%d_y%d|%s|tot' % (p[0], p[1], label) in raw for p in X.DEV):
            tot = sum(raw['OLD|x%d_y%d|%s|tot' % (p[0], p[1], label)] for p in X.DEV)
            ff, c, pu = X.f_measure(tot)
            rows.append(['old family', label, shown(raw['OLD|%s|name' % label][0]), '%.3f' % c, '%.3f' % pu, '%.3f' % ff, '%d' % tot[2], '%d' % (tot[6] - tot[5]), '-'])
    return table(['family', 'variant', 'settings', 'completeness (mag < 22)', 'purity (mag < 22)', 'F', 'detections < 22 mag', 'unmatched detections (all magnitudes)', 'CPU seconds (two patches)'], rows)


def summary_table(summary):
    """One row per pooled group: F (95 % interval, unmatched detections of all magnitudes in brackets) of the tuned ISS, the old tuned baselines and the selected new baselines, and the pooled paired
    differences (tuned ISS minus each new baseline; each new baseline minus its old tuned counterpart) with their verdicts."""
    rows = []
    for title, g, fs, dfs in summary:
        def cell(r):
            if r not in fs:
                return '-'
            (ff, fci) = fs[r]
            tot = g['tot'][r]
            return '%.3f [%.3f, %.3f] (%d)' % (ff, fci[0], fci[1], tot[6] - tot[5])

        def vd(a, b):
            return '-' if (a, b) not in dfs else '%+.3f %s' % (dfs[(a, b)][0], X.verdict(*dfs[(a, b)]))
        rows.append([title + ' (%d patches)' % len(g['done']), cell('iss'), cell('phu_old_t'), cell('cs_old_t'), cell('phu_sel'), cell('cs_sel'), vd('iss', 'phu_sel'), vd('iss', 'cs_sel'),
                     vd('phu_sel', 'phu_old_t'), vd('cs_sel', 'cs_old_t')])
    return ('F at mag < 22 with its 95 % interval and, in brackets, the unmatched detections of all magnitudes; the last four columns are the pooled paired differences dF and the verdict (win: the first '
            'catalogue is better by the 95 % interval, tie: the interval contains 0, loss).\n\n' + table(['group', 'tuned ISS', 'photutils old tuned', 'crowdsource old tuned', 'photutils HP selected',
                                                                                                          'crowdsource selected', 'ISS minus photutils HP', 'ISS minus crowdsource',
                                                                                                          'photutils HP minus old tuned', 'crowdsource minus old tuned'], rows))


def report(raw, frozen=None):
    """The markdown report from the raw of this script and the frozen results (read only)."""
    frozen = frozen if frozen is not None else load_frozen()
    L = ['# The baselines given the ISS\'s front end and a wider pipeline search', '',
         'Generated by `experiments/exp_astro_hp.py` from `baselines_hp_raw.npz` (and the frozen `m2_raw.npz`, `psf_raw.npz`, `repeat_raw.npz`, read only); every number below is printed by the script. The family, the '
         'selection rule and what is reported were declared in `baselines_hp_plan.md` before any catalogue was scored (SHA-256 of the plan stored with the selection: %s). Question (referee): does a baseline '
         'given the ISS\'s front end (the border-consistent high-pass of the frame and of the PSF stamp, the ISS\'s weights) and a search of its own pipeline close the gap to the sliding ISS? Scoring and intervals as '
         '`exp_astro_m2`: match at 0.5 px and 0.5 mag, F-measure of completeness and purity at mag < 22, 95 %% percentile intervals of a bootstrap over the 20 x 20 px blocks (%d resamples, seed 0, the same resampled '
         'blocks for every row of a table, so differences are paired; the pooled tables resample the blocks within each patch). The tuned ISS rows and the old baselines are READ from the frozen results and '
         'influenced nothing declared or selected here; nothing is re-tuned on a test patch. Deviations of the plan from the brief (plan, section 10): the threshold of photutils is in the variance model\'s pixel '
         'noise, not the MAD of the crowded high-passed image; the iteration counts are {3, 10, 30} (mode new) and 6 (mode all); the crowdsource family searches the detection threshold and the front end, with '
         'no post-filters. **On the repeat frames these baselines are added later** (`repeat_plan.md`, section 2: chosen from a declared family on the two M2 development patches only, their catalogues written to a '
         'file of their own, `baseline_catalogs_hp.npz`, before they were scored); they belong to none of the declared outcomes P1 to P3 of that plan, which concern the ISS.'
         % (str(raw['SEL|plan_sha256'][0]) if 'SEL|plan_sha256' in raw else 'none yet', X.BOOT), '']
    shas = frozen.get('sha', {})
    L += ['Frozen results read: %s.' % '; '.join('`%s` %s' % (n, (s[:16] + '...') if s else 'absent') for n, s in (('m2_raw.npz', shas.get('m2')), ('psf_raw.npz', shas.get('psf')), ('repeat_raw.npz', shas.get('repeat')))), '']
    if any(k.startswith('DEV|') for k in raw):
        L += ['## Development patches %s: every variant\n' % ', '.join('(%d, %d)' % p for p in X.DEV),
              'Pooled counts of the two M2 development patches (the only data the selection sees); bold: the variant of each family with the largest F (a tie to the larger tag). The old family\'s rows are for '
              'reference: they take no part in the selection.\n', dev_table(raw), '']
        if 'SEL|photutils' in raw:
            L += ['Selected: photutils `%s`, crowdsource `%s`.' % (raw['SEL|photutils'][0], raw['SEL|crowdsource'][0]), '']
    summary, body = [], []
    for title, name in GROUPS:
        g = group_data(raw, [(name, G.SETS[name][1])], frozen)
        if not g['done']:
            continue
        body += ['## %s\n' % title]
        for _, p in g['done']:
            gp = group_data(raw, [(name, [p])], frozen)
            idx = X.boot_idx_strat(gp['sizes'])
            body += ['### Patch (%d, %d)\n' % p, 'Truth sources brighter than 22 mag in the scored area: %d.\n' % gp['tot']['phu_sel'][0], method_table(gp, idx), '',
                     'Paired differences (first minus second, same blocks), 95 % intervals:\n', diff_table(COMPARISONS_PATCH, gp, idx)[0], '']
        if len(g['done']) > 1:
            body += pooled_block('### All %d test patches together' % len(g['done']), g, summary, title.split(',')[0])
        else:
            summary_entry(g, summary, title.split(',')[0])
    gr = group_data(raw, [(n, G.SETS[n][1]) for n in sorted(R.REPEAT)], frozen)
    if len(gr['done']) > 1:
        body += pooled_block('## All the repeat patches together', gr, summary, 'All the repeat patches')
    if summary:
        L += ['## Summary: pooled F at mag < 22 by cluster\n', summary_table(summary), '']
    return '\n'.join(L + body) + '\n'


def summary_entry(g, summary, title):
    """Appends (title, group, {row: (F, interval)}, {pair: dF}) for the summary table."""
    idx = X.boot_idx_strat(g['sizes'])
    summary.append((title, g, {r: X.boot_cpf(g['blk'][r], idx)[2] for r in g['rows']}, diff_table(COMPARISONS_POOLED, g, idx)[1]))


def pooled_block(heading, g, summary, title):
    """The lines of the pooled tables of a group of patches (stratified block bootstrap); the group also goes into the summary."""
    idx = X.boot_idx_strat(g['sizes'])
    summary_entry(g, summary, title)
    return [heading + '\n', 'The blocks of the %d patches pooled (%d blocks); counts are sums; the bootstrap resamples the blocks WITHIN each patch (stratified), the same resamples for every row. The patches '
            'are windows of frames of one cluster, not independent trials.\n' % (len(g['done']), sum(g['sizes'])), method_table(g, idx), '', 'Paired differences on the pooled blocks (first minus second):\n',
            diff_table(COMPARISONS_POOLED, g, idx)[0], '']


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--dev', action='store_true')
    ap.add_argument('--old', action='store_true')
    ap.add_argument('--test', action='store_true')
    ap.add_argument('--sets', nargs='*', default=['m2', 'm53', 'm92'])
    ap.add_argument('--flag', help='the file that says the ISS repeat run is complete (needed to score the repeat frames)')
    ap.add_argument('--from-raw', action='store_true')
    A_ = ap.parse_args()
    if A_.dev:
        print(part_dev())
    elif A_.old:
        part_old()
    elif A_.test:
        part_test(A_.sets, flag=A_.flag)
    elif A_.from_raw:
        with open(MD, 'w') as f:
            f.write(report(X.load_raw(RAW)))
        print('wrote', MD)
    else:
        ap.error('give --dev, --old, --test or --from-raw')


if __name__ == '__main__':
    main()
