"""
Held-out test on repeat imaging of M2: the frozen configuration and the kernel rule on the same four windows in other SDSS runs, at other seeing
=================================================================================================================================================

Declared before the run in `results/astro/repeat_plan.md` (read it first; its hash is in the freeze record). The data are the six repeat patches of `experiments/astro_repeat.py` (frames m2r1, m2r2,
m2r3: SDSS runs 4187, 7778, 3388, psField FWHM 2.28, 3.05, 3.52 px; M2's frame has 2.1 px), the baselines' catalogues at their frozen settings (`baseline_catalogs_repeat.npz`, written before this
script existed) and the configuration frozen on the two M2 development patches (`freeze.json`): sigma_b 3 px, r_ex 1.2 px, kappa 2; the a priori default (5, 1.2, 5); the BIC stop.

Per patch, once: the stamp fit (label-free); the sliding ISS with the isotropic mixture of `imaging.fit_psf_parts` (as in the frozen M2 test, `exp_astro_m2.run_path`) and with the kernel that
`psf_general.fit_adequate` chooses from the stamp (`exp_astro_psf.run_path_psf`; where the rule keeps the isotropic fit the two are the same program on the same input, the run is done once and stored
under both names, flag `same`), each at the tuned configuration (kappa 2), the a priori default (kappa 5) and the default with the BIC stop (key 0); an EXPLORATORY run with the lengths sigma_b and r_ex
multiplied by the seeing (FWHM of the stamp over the mean FWHM of the stamps of the two M2 development patches); the baselines as tuned and at their defaults. Scored as everywhere
(`astro_common.score`: 0.5 px and 0.5 mag, one to one, F of completeness and purity at mag < 22). Nothing is re-tuned: the configuration, the kernel rule and the baselines' choices are read from the
freeze records and checked.

    python experiments/exp_astro_repeat.py --freeze            write results/astro/freeze_repeat.json (hashes of the code, the inputs, the plan) -- before the run
    python experiments/exp_astro_repeat.py --run [--workers 2]   the ISS on the six patches, once each (refuses unless the freeze record matches; a patch whose result exists is not run again)
    python experiments/exp_astro_repeat.py --from-raw          regenerate results/astro/repeat.md from repeat_raw.npz (reads the frozen m2_raw.npz for the seeing ladder)
"""
import argparse
import hashlib
import inspect
import json
import os
import sys
import time

os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np

HERE = os.path.abspath(os.path.dirname(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)
import astro_common as A                  # noqa: E402
import astro_repeat as R                  # noqa: E402
import astro_second as S2                 # noqa: E402
import exp_astro_m2 as X                  # noqa: E402
import exp_astro_psf as PS                # noqa: E402
from src.core import imaging as IM        # noqa: E402
from src.core import psf_general as G     # noqa: E402

OUT = os.environ.get('REPEAT_OUT', os.path.join(REPO, 'results', 'astro'))
RAW, MD = os.path.join(OUT, 'repeat_raw.npz'), os.path.join(OUT, 'repeat.md')
FREEZE = os.path.join(OUT, 'freeze_repeat.json')
PLAN = os.path.join(REPO, 'results', 'astro', 'repeat_plan.md')
M2FREEZE = os.path.join(REPO, 'results', 'astro', 'freeze.json')
BASE = os.path.join(REPO, 'results', 'astro', 'baseline_catalogs_repeat.npz')
M2RAW = os.environ.get('REPEAT_M2_RAW', X.RAW)           # the frozen M2 results, read only
PATCHES = [(name, pxy) for name in sorted(R.REPEAT) for pxy in R.REPEAT[name]['patches']]
FILES = ('src/core/imaging.py', 'src/core/psf_general.py', 'experiments/astro_common.py', 'experiments/astro_repeat.py', 'experiments/astro_second.py', 'experiments/exp_astro_baselines.py',
         'results/astro/baseline_catalogs_repeat.npz', 'results/astro/repeat_plan.md', 'results/astro/freeze.json', 'data/astro/checksums_repeat.json')
LABELS = ('ours', 'default')              # the tuned configuration (kappa from the freeze) and the a priori default (kappa DEFAULT[2]); each also with the BIC stop (key 0)


def load_frozen(path=M2FREEZE):
    """The configuration frozen on the M2 development patches: dict(hp, rex, kap, default, tuned) from `freeze.json`; raises unless the a priori default is the code's `exp_astro_m2.DEFAULT`."""
    with open(path) as f:
        rec = json.load(f)
    c = rec['config']
    if tuple(rec['default']) != tuple(X.DEFAULT):
        raise RuntimeError('the default of freeze.json is not exp_astro_m2.DEFAULT')
    return dict(hp=float(c['sigma_b']), rex=float(c['r_ex']), kap=float(c['kappa']), default=tuple(rec['default']), tuned=dict(rec['tuned_baselines']))


def ref_fwhm():
    """The FWHM (px) of the stamps of the two M2 development patches, their mean: the seeing against which the exploratory runs scale the lengths (label-free)."""
    return float(np.mean([S2.psf_fwhm(X.load_patch(*p)['psf']) for p in X.DEV]))


def _store(arrays, prefix, dets, secs, n, info, t, dk):
    """Scores of the catalogues of one run (every kappa and the BIC key 0), the run time and entrants, the BIC block, and the catalogue at the reported significance dk."""
    for k, d in dets.items():
        arrays.update(X._pack('%s|k%g' % (prefix, k), X.score_dets(d, t)))
    arrays['%s|secs' % prefix] = np.array([secs, n])
    arrays['%s|bic' % prefix] = np.array([info['bic_index'], info['n_blocks'], int(n >= X.MAX_ENTRIES)])
    arrays['%s|det' % prefix] = np.c_[dets[dk][0], dets[dk][1], dets[dk][2]]


def _alias(arrays, src, dst):
    """The arrays stored under `src|...` also under `dst|...` (the same run, used for two rows)."""
    for k in [k for k in arrays if k.startswith(src + '|')]:
        arrays[dst + k[len(src):]] = arrays[k]


def patch_arrays(args):
    """The arrays of one repeat patch: stamp fit, the ISS with both kernels at both configurations (and the BIC stop), the exploratory scaled runs, the baselines. A function of its arguments only."""
    name, pxy, hp, rex, kap, tuned, fw_ref, base = args
    patch = R.load_patch(pxy[0], pxy[1], name)
    t = A.truth_for_patch(R.truth_arrays(name), pxy[0], pxy[1])
    arrays, fit = PS.stamp_arrays(patch)
    fwhm = S2.psf_fwhm(patch['psf'])
    arrays['fwhm'] = np.array([fwhm, fw_ref])
    arrays['cfg'] = np.array([hp, rex, kap])
    arrays['kernel|same'] = np.array([int(fit['kind'] == 'isotropic')])
    sig, parts, _ = IM.fit_psf_parts(patch['psf'], n_parts=3, sigma0=1.0)
    for label, (h, r), dk in (('ours', (hp, rex), kap), ('default', (X.DEFAULT[0], X.DEFAULT[1]), X.DEFAULT[2])):
        dets, secs, n, info = X.run_path(patch, h, r)
        _store(arrays, 'iso|' + label, dets, secs, n, info, t, dk)
        P_iso = IM.PixelPSF(*patch['dn'].shape, sig, parts=parts, hp=h or None)
        arrays['iso|%s|rsd' % label] = np.array([PS.residual_sd(patch, P_iso, np.c_[dets[dk][0], dets[dk][1], dets[dk][2]])])
        print('R patch %s %s %s iso: sigma_b %g r_ex %g, %d entries in %.0f s, BIC block %d of %d' % (name, pxy, label, h, r, n, secs, info['bic_index'], info['n_blocks']), flush=True)
        if fit['kind'] == 'isotropic':
            _alias(arrays, 'iso|' + label, 'adq|' + label)
        else:
            dets2, secs2, n2, info2 = PS.run_path_psf(patch, h, r, fit)
            _store(arrays, 'adq|' + label, dets2, secs2, n2, info2, t, dk)
            arrays['adq|%s|rsd' % label] = np.array([info2['rsd'][dk]])
            print('R patch %s %s %s adequate (%d parts): %d entries in %.0f s, BIC block %d of %d' % (name, pxy, label, len(fit['parts5']), n2, secs2, info2['bic_index'], info2['n_blocks']), flush=True)
    s = fwhm / fw_ref                                              # exploratory: lengths in units of the seeing, the kernel the rule chose
    arrays['scale'] = np.array([s])
    for label, (h, r), dk in (('ours', (hp, rex), kap), ('default', (X.DEFAULT[0], X.DEFAULT[1]), X.DEFAULT[2])):
        dets3, secs3, n3, info3 = PS.run_path_psf(patch, h * s, r * s, fit)
        _store(arrays, 'scl|' + label, dets3, secs3, n3, info3, t, dk)
        arrays['scl|%s|rsd' % label] = np.array([info3['rsd'][dk]])
        print('R patch %s %s %s scaled x%.2f: %d entries in %.0f s' % (name, pxy, label, s, n3, secs3), flush=True)
    bs = X.baseline_scores(pxy, t, base)
    for label, nm in (('photutils (default)', 'photutils'), ('crowdsource (default)', X.cs_name((0.0, 0.0, 0.0, np.inf)))):
        if nm in bs:
            arrays.update(X._pack(label, bs[nm]))
    for fam, nm in tuned.items():
        if nm in bs:
            arrays.update(X._pack('%s (tuned)' % fam, bs[nm]))
            arrays['%s (tuned)|name' % fam] = np.array([nm])
    return name, pxy, arrays


# ----------------------------------------------------------------------------------------------------------------------------------
# freeze
# ----------------------------------------------------------------------------------------------------------------------------------
def _sha256_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def frozen_code_hash():
    """The hash of what decides a number of this test: the source of this script's run functions, of the functions of `exp_astro_m2` and `exp_astro_psf` that run and score the path, and the constants
    (patches, grids, default, baseline families, kernel rule). Edits of the report do not change it."""
    src = ''.join(inspect.getsource(fn) for fn in (load_frozen, ref_fwhm, _store, _alias, patch_arrays, X.run_path, X.load_patch, X.score_dets, X._pack, X.f_measure, X.baseline_scores,
                                                  X.cs_dets, X.cs_name, PS.run_path_psf, PS.residual_sd, PS.stamp_arrays, PS.stamp_shape, R.load_patch, R.truth_arrays, R.patch_path))
    consts = repr((PATCHES, X.DEFAULT, X.KAPPAS, X.LAT_H, X.MARGIN, X.FIT_RADIUS, X.MAX_ENTRIES, X.PHOT_TAGS, X.CS_FILTERS, G.RMS_GOAL, G.K_MAX, X.DEV, LABELS, FILES))
    return hashlib.sha256((src + consts).encode()).hexdigest()


def freeze_record():
    """What is frozen: the configuration and baseline choices read from freeze.json, the hashes of the code and of every input (code files, catalogues, plan, patch and truth files)."""
    fr = load_frozen()
    data = {}
    for name, pxy in PATCHES:
        data[R.patch_path(pxy[0], pxy[1], name)] = None
    for name in R.REPEAT:
        data[os.path.join(R.REPEAT[name]['root'], R.REPEAT[name]['truth'])] = None
    data_files = {os.path.relpath(p, R.ROOT3): _sha256_file(p) for p in sorted(data)}
    return dict(config=dict(sigma_b=fr['hp'], r_ex=fr['rex'], kappa=fr['kap']), default=list(fr['default']), tuned_baselines=fr['tuned'], fwhm_ref=round(ref_fwhm(), 6),
                files={f: _sha256_file(os.path.join(REPO, f)) for f in FILES}, data_files=data_files, frozen_code=frozen_code_hash(), patches=[[n, list(p)] for n, p in PATCHES])


def write_freeze(path=None):
    path = path or FREEZE
    rec = freeze_record()
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as f:
        json.dump(rec, f, indent=1, sort_keys=True)
        f.write('\n')
    return rec


def check_freeze(path=None):
    """Raises unless the freeze record exists and equals what the code, the inputs and the plan are now."""
    path = path or FREEZE
    if not os.path.exists(path):
        raise RuntimeError('no freeze record %s: run `exp_astro_repeat.py --freeze` first' % path)
    with open(path) as f:
        old = json.load(f)
    now = freeze_record()
    diff = [k for k in now if now[k] != old.get(k)]
    if diff:
        raise RuntimeError('the freeze record differs from the present code / inputs in: %s' % ', '.join(diff))
    return old


def part_run(out_path=RAW, workers=1, freeze_path=None, base=BASE):
    """Each repeat patch, once, at the frozen configuration; a patch whose result exists is not run again."""
    rec = check_freeze(freeze_path)
    fr = load_frozen()
    cur = X.load_raw(out_path) if os.path.exists(out_path) else {}
    todo = [(n, p) for n, p in PATCHES if 'R|%s|x%d_y%d|cfg' % (n, p[0], p[1]) not in cur]
    jobs = [(n, p, fr['hp'], fr['rex'], fr['kap'], fr['tuned'], rec['fwhm_ref'], base) for n, p in todo]
    print('R: %d patches to run, %d exist' % (len(jobs), len(PATCHES) - len(jobs)), flush=True)
    if workers > 1 and len(jobs) > 1:
        from multiprocessing import Pool
        with Pool(min(workers, len(jobs))) as pool:
            for name, pxy, arrays in pool.imap_unordered(patch_arrays, jobs, chunksize=1):
                X.save_part(out_path, 'R|%s|x%d_y%d' % (name, pxy[0], pxy[1]), arrays)
    else:
        for j in jobs:
            name, pxy, arrays = patch_arrays(j)
            X.save_part(out_path, 'R|%s|x%d_y%d' % (name, pxy[0], pxy[1]), arrays)


# ----------------------------------------------------------------------------------------------------------------------------------
# report
# ----------------------------------------------------------------------------------------------------------------------------------
KERN = {'iso': 'isotropic kernel (as frozen)', 'adq': 'kernel chosen from the stamp', 'scl': 'kernel chosen, lengths scaled by the seeing (exploratory)'}
CFG = {'ours': ('tuned configuration (sigma_b 3, r_ex 1.2, kappa 2)', 'k2'), 'default': ('a priori default (5, 1.2, kappa 5)', 'k5')}
BASE_ROWS = [('crowdsource (default)', 'crowdsource (default)'), ('crowdsource (tuned)', 'crowdsource (tuned)'), ('photutils (default)', 'photutils (default)'), ('photutils (tuned)', 'photutils (tuned)')]


def rows_def():
    """[(label, key)] of the rows of every table: the ISS with the isotropic and the chosen kernel at the three readouts, the baselines; the exploratory rows are separate."""
    out = []
    for kern in ('iso', 'adq'):
        out += [('ISS, %s, %s' % (KERN[kern], CFG['ours'][0]), '%s|ours|k2' % kern), ('ISS, %s, %s' % (KERN[kern], CFG['default'][0]), '%s|default|k5' % kern),
                ('ISS, %s, a priori default with the BIC stop' % KERN[kern], '%s|default|k0' % kern)]
    return out + BASE_ROWS


EXPLORATORY = [('ISS, %s, %s' % (KERN['scl'], CFG[c][0]), 'scl|%s|%s' % (c, CFG[c][1])) for c in ('ours', 'default')]


def pre_of(name, pxy):
    return 'R|%s|x%d_y%d|' % (name, pxy[0], pxy[1])


def done_patches(raw):
    return [(n, p) for n, p in PATCHES if pre_of(n, p) + 'cfg' in raw]


def purity_all(tot):
    """Matched detections over all detections of any magnitude, and the number of unmatched detections, from the totals [true22, matched22, det22, det22_matched, true_all, matched_all, det_all]."""
    return tot[5] / max(tot[6], 1), int(tot[6] - tot[5])


def pooled_blocks(raw, plist, key):
    parts = [raw[pre_of(*np_) + key + '|blk22'] for np_ in plist]
    return np.concatenate(parts), [len(q) for q in parts]


def diff_row(name, a, b, idx):
    dc, dp, df = X.boot_diff(a, b, idx)
    return [name, X.fmt_iv(dc[0], dc[1]), X.fmt_iv(dp[0], dp[1]), X.fmt_iv(df[0], df[1]), X.verdict(*df)], df


def method_table(raw, plist, defs, title_note=''):
    """The rows `defs` pooled over the patches `plist` (stratified block bootstrap): completeness, purity, F, purity over all magnitudes."""
    blk, sizes = {}, None
    for label, key in defs:
        if all(pre_of(*np_) + key + '|blk22' in raw for np_ in plist):
            blk[key], sizes = pooled_blocks(raw, plist, key)
    if not blk:
        return '', {}, None
    idx = X.boot_idx_strat(sizes)
    rows = []
    for label, key in defs:
        if key not in blk:
            continue
        tot = sum(raw[pre_of(*np_) + key + '|tot'] for np_ in plist)
        (c, cci), (pu, pci), (ff, fci) = X.boot_cpf(blk[key], idx)
        pa, un = purity_all(tot)
        rows.append([label, '%d / %d' % (tot[1], tot[0]), X.fmt_iv(c, cci), '%d / %d' % (tot[3], tot[2]), X.fmt_iv(pu, pci), X.fmt_iv(ff, fci), '%.3f (%d)' % (pa, un)])
    return X.table(['catalogue', 'matched / truth < 22', 'completeness', 'matched / detected < 22', 'purity', 'F', 'purity over all magnitudes (unmatched detections)'], rows), blk, idx


PAIRS_MAIN = [('adequate kernel vs isotropic, tuned configuration', 'adq|ours|k2', 'iso|ours|k2'), ('adequate kernel vs isotropic, a priori default', 'adq|default|k5', 'iso|default|k5'),
              ('adequate kernel vs isotropic, BIC stop', 'adq|default|k0', 'iso|default|k0'),
              ('adequate kernel (tuned configuration) vs crowdsource (tuned)', 'adq|ours|k2', 'crowdsource (tuned)'), ('adequate kernel (tuned configuration) vs photutils (tuned)', 'adq|ours|k2', 'photutils (tuned)'),
              ('adequate kernel (a priori default) vs crowdsource (default)', 'adq|default|k5', 'crowdsource (default)'), ('adequate kernel (BIC stop) vs crowdsource (default)', 'adq|default|k0', 'crowdsource (default)'),
              ('isotropic kernel (tuned configuration) vs crowdsource (tuned)', 'iso|ours|k2', 'crowdsource (tuned)'), ('isotropic kernel (tuned configuration) vs photutils (tuned)', 'iso|ours|k2', 'photutils (tuned)')]


def pair_table(blk, idx):
    rows, dfs = [], {}
    for name, a, b in PAIRS_MAIN:
        if a in blk and b in blk:
            r, df = diff_row(name, blk[a], blk[b], idx)
            rows.append(r)
            dfs[(a, b)] = df
    return X.table(['comparison', 'd completeness', 'd purity', 'd F', 'verdict on F'], rows), dfs


def p1_text(raw, plist):
    """The declared primary outcome P1: on the patches where the rule picks a general kernel, the pooled paired dF of the adequate against the isotropic kernel."""
    gen = [np_ for np_ in plist if raw[pre_of(*np_) + 'kernel|same'][0] == 0]
    if not gen:
        return 'P1: the rule keeps the isotropic kernel on every patch run; the two are the same program, dF = 0 by construction, and P1 has nothing to test.'
    lines = ['P1 (kernel rule on held-out frames), the %d patches where the rule picks a general kernel (%s):' % (len(gen), ', '.join('%s (%d, %d)' % (n, p[0], p[1]) for n, p in gen))]
    for label, a, b in (('tuned configuration (kappa 2)', 'adq|ours|k2', 'iso|ours|k2'), ('a priori default (kappa 5)', 'adq|default|k5', 'iso|default|k5'), ('BIC stop', 'adq|default|k0', 'iso|default|k0')):
        ba, sizes = pooled_blocks(raw, gen, a)
        bb, _ = pooled_blocks(raw, gen, b)
        idx = X.boot_idx_strat(sizes)
        df = X.boot_diff(ba, bb, idx)[2]
        v = {'win': 'confirmed on held-out frames', 'tie': 'not confirmed', 'loss': 'contradicted'}[X.verdict(*df)]
        lines.append('* %s: dF = %s, %s' % (label, X.fmt_iv(df[0], df[1]), v))
    return '\n'.join(lines)


def report(raw, m2raw=None):
    plist = done_patches(raw)
    L = ['# Held-out test on repeat imaging of M2: the frozen configuration and the kernel rule at other seeing', '',
         'Generated by `experiments/exp_astro_repeat.py` from `repeat_raw.npz` (and the frozen `m2_raw.npz`, read only); every number below is printed by the script. The analysis was declared in `repeat_plan.md` '
         'before any run; the freeze record `freeze_repeat.json` hashes the plan, the code and the inputs, and the run refused to start unless it matched. **Nothing was re-tuned on these patches**: the configuration '
         '(sigma_b 3, r_ex 1.2, kappa 2), the a priori default (5, 1.2, 5), the BIC stop, the kernel rule (`psf_general.fit_adequate`, rms goal %g of the peak) and the baselines\' choices come from the M2 development '
         'patches and the M53/M92 episode. The windows are the M2 test windows (same sky, same HUGS stars): images, PSFs and noise are new, the truth and its defects are shared. Six patches from three frames: the '
         'unit of replication for the seeing is the frame. Intervals: 95 %% percentile intervals of a bootstrap over the 20 x 20 px blocks (%d resamples, seed 0, stratified within each patch, the same resamples '
         'for every row of a table, so differences are paired); they are conditional on the patches. Magnitudes from counts use M2\'s fixed zero point for every frame (0.03 to 0.08 mag offset, `repeat_spec.md`).'
         % (G.RMS_GOAL, X.BOOT), '']
    rows = []
    for name, p in plist:
        pre = pre_of(name, p)
        rms, kind, p5 = raw[pre + 'rms'], str(raw[pre + 'kind'][0]), raw[pre + 'parts5']
        sh = raw[pre + 'shape']
        tr = raw[pre + 'iso|ours|k2|tot']
        rows.append([name, '(%d, %d)' % p, '(%d, %d)' % R.REPEAT[name]['windows'][R.REPEAT[name]['patches'].index(tuple(p))], '%.2f' % raw[pre + 'fwhm'][0], '%d' % tr[0], '%.5f' % rms[0], kind,
                     '%d' % (len(p5) if kind == 'general' else 3), '%.5f' % rms[1], '%.3f' % sh[0]])
    L += ['## The patches and their stamp fits (label-free)\n', 'Frame, patch origin, the M2 window it repeats, FWHM of the psField stamp (px; the M2 development stamps have %.2f), truth sources brighter than 22 mag in the scored area, rms of the isotropic fit '
          'against the stamp\'s peak, the kernel the rule chose, its number of parts and rms, the ellipticity of the stamp.\n' % (raw[pre_of(*plist[0]) + 'fwhm'][1] if plist else float('nan')),
          X.table(['frame', 'patch', 'M2 window', 'FWHM', 'truth < 22', 'isotropic rms', 'kernel chosen', 'parts', 'its rms', 'ellipticity'], rows), '']
    if plist:
        L += ['## All %d patches together\n' % len(plist), 'The blocks of the patches pooled (%d blocks); counts are sums.\n' % (25 * len(plist))]
        tbl, blk, idx = method_table(raw, plist, rows_def())
        L += [tbl, '', 'Paired differences on the pooled blocks:\n']
        ptbl, dfs = pair_table(blk, idx)
        L += [ptbl, '', '### Primary outcome P1\n', p1_text(raw, plist), '']
        for fname in sorted(R.REPEAT):
            fl = [np_ for np_ in plist if np_[0] == fname]
            if not fl:
                continue
            L += ['## Frame %s: %s\n' % (fname, R.REPEAT[fname]['title']), 'The %d patches of this frame pooled.\n' % len(fl)]
            tbl, blk, idx = method_table(raw, fl, rows_def())
            ptbl, _ = pair_table(blk, idx)
            L += [tbl, '', ptbl, '']
        L += ['## Record by patch (verdict on F, paired intervals)\n']
        cmp_rows = []
        for name_, a, b in PAIRS_MAIN:
            cells = []
            for np_ in plist:
                pre = pre_of(*np_)
                if pre + a + '|blk22' in raw and pre + b + '|blk22' in raw:
                    idx1 = X.boot_idx(raw[pre + a + '|blk22'].shape[0])
                    df = X.boot_diff(raw[pre + a + '|blk22'], raw[pre + b + '|blk22'], idx1)[2]
                    cells.append('%s (%+.3f)' % (X.verdict(*df), df[0]))
                else:
                    cells.append('-')
            cmp_rows.append([name_] + cells)
        L += [X.table(['comparison'] + ['%s (%d, %d)' % (n, p[0], p[1]) for n, p in plist], cmp_rows), '']
        L += ['## P3: the residual scale of the catalogues (label-free)\n', 'Robust standard deviation of the standardized residual (1: the noise model and the kernel are right), at the reported significance; the isotropic and the chosen kernel; '
              '"same" where the rule keeps the isotropic kernel.\n']
        rows = []
        for name, p in plist:
            pre = pre_of(name, p)
            same = raw[pre + 'kernel|same'][0] == 1
            rows.append([name, '(%d, %d)' % p, '%.2f' % raw[pre + 'iso|ours|rsd'][0], 'same' if same else '%.2f' % raw[pre + 'adq|ours|rsd'][0], '%.2f' % raw[pre + 'iso|default|rsd'][0],
                         'same' if same else '%.2f' % raw[pre + 'adq|default|rsd'][0]])
        L += [X.table(['frame', 'patch', 'isotropic, tuned', 'chosen kernel, tuned', 'isotropic, default', 'chosen kernel, default'], rows), '']
        L += ['## Run records\n', 'Entrants, run time of the path to its last entrant (one core, shared machine), the BIC block of the path and whether the entry cap (%d) was hit (a capped run is censored).\n' % X.MAX_ENTRIES]
        rows = []
        for name, p in plist:
            pre = pre_of(name, p)
            for kern in ('iso', 'adq'):
                if kern == 'adq' and raw[pre + 'kernel|same'][0] == 1:
                    continue
                for c in LABELS:
                    s_, b_ = raw[pre + '%s|%s|secs' % (kern, c)], raw[pre + '%s|%s|bic' % (kern, c)]
                    rows.append([name, '(%d, %d)' % p, KERN[kern], c, '%d' % s_[1], '%.0f' % s_[0], '%d of %d' % (b_[0], b_[1]), 'capped' if b_[2] else '-'])
        L += [X.table(['frame', 'patch', 'kernel', 'configuration', 'entrants', 'seconds', 'BIC block', 'cap'], rows), '']
        L += ['## Seeing ladder on the M2 windows\n', 'F at mag < 22 of the same window in the frames that contain it: M2 (the frozen test, isotropic kernel, which is also the chosen kernel there) and the repeat frames (the chosen kernel; the isotropic '
              'kernel in brackets where it differs). Tuned configuration; tuned baselines. Single patches: conditional on the patch.\n']
        froz = m2raw if m2raw is not None else X.load_raw(M2RAW)
        rows = []
        for w in X.TESTS:
            cells = []
            for fname, key in (('M2', 'T|x%d_y%d|' % w),) + tuple((n, pre_of(n, p)) for n, p in plist if tuple(R.REPEAT[n]['windows'][R.REPEAT[n]['patches'].index(tuple(p))]) == tuple(w)):
                src = froz if fname == 'M2' else raw
                row = []
                for mk in (('ours|k2',) if fname == 'M2' else ('adq|ours|k2',)) + ('crowdsource (tuned)', 'photutils (tuned)'):
                    b = src.get(key + mk + '|blk22')
                    row.append('-' if b is None else '%.3f' % X.cpf(b.sum(0))[2])
                if fname != 'M2':
                    bi = raw.get(key + 'iso|ours|k2|blk22')
                    row[0] += ' (%.3f)' % X.cpf(bi.sum(0))[2] if bi is not None and raw[key + 'kernel|same'][0] == 0 else ''
                cells.append('%s: ISS %s, crowdsource %s, photutils %s' % ((fname,) + tuple(row)))
            rows.append(['(%d, %d)' % w] + cells + ['-'] * (4 - len(cells)))
        L += [X.table(['M2 window', 'frame 1', 'frame 2', 'frame 3', 'frame 4'], rows), '']
        if all(pre_of(*np_) + 'scl|ours|k2|blk22' in raw for np_ in plist):
            L += ['## Exploratory: lengths scaled by the seeing\n', 'DESCRIPTIVE, not a result of the declared plan (a question: do the pixel-unit constants need to follow the seeing?). sigma_b and r_ex multiplied by s = FWHM of the stamp over '
                  'the mean FWHM of the stamps of the two M2 development patches (%s), kappa unchanged, the kernel the rule chose. Pooled over the six patches.\n'
                  % ', '.join('%s %.2f' % (n, raw[pre_of(n, p) + 'scale'][0]) for n, p in plist)]
            tbl, blk, idx = method_table(raw, plist, [('ISS, chosen kernel, tuned configuration (pixels)', 'adq|ours|k2')] + EXPLORATORY + [('ISS, chosen kernel, a priori default (pixels)', 'adq|default|k5')] + BASE_ROWS)
            L += [tbl, '']
    return '\n'.join(L) + '\n'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--freeze', action='store_true')
    ap.add_argument('--run', action='store_true')
    ap.add_argument('--from-raw', action='store_true')
    ap.add_argument('--workers', type=int, default=int(os.environ.get('ASTRO_WORKERS', '2')))
    A_ = ap.parse_args()
    if A_.freeze:
        print(json.dumps(write_freeze(), indent=1, sort_keys=True))
    elif A_.run:
        part_run(workers=A_.workers)
    elif A_.from_raw:
        with open(MD, 'w') as f:
            f.write(report(X.load_raw(RAW)))
        print('wrote', MD)
    else:
        ap.error('give --freeze, --run or --from-raw')


if __name__ == '__main__':
    main()
