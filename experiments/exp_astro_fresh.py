"""
Fresh-cluster test: the frozen configuration and the kernel rule on two globular clusters that no choice of the project has seen (M3 and M15), against the baselines as frozen and as given the ISS's front end
===========================================================================================================================================================================================================

Declared before the run in `results/astro/fresh_plan.md` (read it first; its hash is in the freeze record; the choice of the clusters is `fresh_rule.md`). The data are the patches of
`results/astro/fresh_patches.json` (`experiments/astro_fresh.py`: M3 one patch, M15 two), the baselines' catalogues (`fresh_baseline_catalogs.npz`, written by `exp_astro_fresh_baselines.py` before this script
was run, no truth read) and the configuration frozen on the two M2 development patches (`freeze.json`): sigma_b 3 px, r_ex 1.2 px, kappa 2; the a priori default (5, 1.2, 5); the BIC stop.

Per patch, once: the stamp fit (label-free); the sliding ISS with the isotropic mixture of `imaging.fit_psf_parts` (`exp_astro_m2.run_path`, the program frozen on M2) and with the kernel that
`psf_general.fit_adequate` chooses from the stamp (`exp_astro_psf.run_path_psf`; where the rule keeps the isotropic fit the two are the same program on the same input, the run is done once and stored under both
names, flag `same`), each at the tuned configuration (kappa 2), the a priori default (kappa 5) and the default with the BIC stop (key 0); the baselines (photutils and crowdsource: default, tuned as frozen in
`freeze.json`, and the two variants selected in round 2 of the baselines given the ISS's front end), all scored against the truth as it is ('' keys) and against the truth shifted by the rigid offset measured
beforehand ('reg|' keys, secondary). EVERY ISS run uses the guarded path (`imaging_guard.sliding_path_guarded`, through the wrapper of `exp_astro_cyclefix.install`); the run records the cap flag, the entry count and the
guard activations. Nothing is re-tuned: the configuration, the kernel rule and the baselines' choices are read from the freeze records and checked. The unit of replication is the CLUSTER (n = 2).

    python experiments/exp_astro_fresh.py --freeze             write results/astro/freeze_fresh.json (hashes of the code, the inputs, the plan) -- before the run
    python experiments/exp_astro_fresh.py --run [--workers 1]  the ISS on the patches, once each (refuses unless the freeze record matches; a patch whose result exists is not run again)
    python experiments/exp_astro_fresh.py --from-raw           regenerate results/astro/fresh.md from fresh_raw.npz (reads nothing else)

Environment (tests): FRESH_OUT (directory of the raw, the md and the freeze record), FRESH_PLAN, FRESH_BASE (the baseline catalogue file), FRESH_REGISTRY (the registry of `astro_fresh.py`).
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
import astro_fresh as F                   # noqa: E402
import astro_second as S2                 # noqa: E402
import exp_astro_baselines_repeat as BR   # noqa: E402  (parse_cs_name)
import exp_astro_cyclefix as CF           # noqa: E402  (the guarded path wrapper and its record)
import exp_astro_m2 as X                  # noqa: E402
import exp_astro_psf as PS                # noqa: E402
from src.core import imaging as IM        # noqa: E402
from src.core import psf_general as G     # noqa: E402

ASTRO = os.path.join(REPO, 'results', 'astro')
OUT = os.environ.get('FRESH_OUT', ASTRO)
RAW, MD = os.path.join(OUT, 'fresh_raw.npz'), os.path.join(OUT, 'fresh.md')
FREEZE = os.path.join(OUT, 'freeze_fresh.json')
PLAN = os.environ.get('FRESH_PLAN', os.path.join(ASTRO, 'fresh_plan.md'))
M2FREEZE = os.path.join(ASTRO, 'freeze.json')
BASE = os.environ.get('FRESH_BASE', os.path.join(ASTRO, 'fresh_baseline_catalogs.npz'))
FILES = ('src/core/imaging.py', 'src/core/psf_general.py', 'src/core/imaging_guard.py', 'experiments/astro_common.py', 'experiments/astro_second.py', 'experiments/astro_fresh.py',
         'experiments/exp_astro_m2.py', 'experiments/exp_astro_psf.py', 'experiments/exp_astro_baselines.py', 'experiments/exp_astro_baselines_repeat.py', 'experiments/exp_astro_cyclefix.py',
         'experiments/exp_astro_fresh_baselines.py', 'results/astro/fresh_baseline_catalogs.npz', 'results/astro/fresh_plan.md', 'results/astro/fresh_rule.md', 'results/astro/fresh_survey.txt',
         'results/astro/fresh_fields.json', 'results/astro/fresh_selection.txt', 'results/astro/fresh_selection_A1.txt', 'results/astro/fresh_patches.json', 'results/astro/freeze.json',
         'results/astro/baselines_hp2_raw.npz', 'data/astro/checksums_fresh.json')
LABELS = ('ours', 'default')              # the tuned configuration (kappa from freeze.json) and the a priori default (kappa DEFAULT[2]); each also with the BIC stop (key 0)
TRUTHS = ('', 'reg|')                     # key prefixes of the scoring: the truth as it is (primary), the truth shifted by the measured rigid offset (secondary)
# baseline rows: label, tag source ('photutils' = the default catalogue; 'tuned_phot'; 'cs_default'; 'cs_tuned'; 'fe_phot'; 'fe_cs')
BASE_ROWS = [('photutils (default)', 'photutils'), ('photutils (tuned)', 'tuned_phot'), ('photutils (front end)', 'fe_phot'),
             ('crowdsource (default)', 'cs_default'), ('crowdsource (tuned)', 'cs_tuned'), ('crowdsource (front end)', 'fe_cs')]
BASE_LABELS = tuple(l for l, _ in BASE_ROWS)


def registry():
    return F.load_registry()


def patches(reg=None):
    """[(cluster, (x0, y0))] of the registry, clusters in name order."""
    reg = reg if reg is not None else registry()
    return [(c, tuple(p)) for c in sorted(reg) for p in reg[c]['patches']]


try:
    PATCHES = patches()
except (OSError, KeyError, ValueError):
    PATCHES = []


def pre_of(cluster, pxy):
    return 'F|%s|x%d_y%d|' % (cluster, pxy[0], pxy[1])


def load_frozen(path=M2FREEZE):
    """The configuration frozen on the M2 development patches: dict(hp, rex, kap, default, tuned) from `freeze.json`; raises unless the a priori default is the code's `exp_astro_m2.DEFAULT`."""
    with open(path) as f:
        rec = json.load(f)
    c = rec['config']
    if tuple(rec['default']) != tuple(X.DEFAULT):
        raise RuntimeError('the default of freeze.json is not exp_astro_m2.DEFAULT')
    return dict(hp=float(c['sigma_b']), rex=float(c['r_ex']), kap=float(c['kappa']), default=tuple(rec['default']), tuned=dict(rec['tuned_baselines']))


# ----------------------------------------------------------------------------------------------------------------------------------
# the baselines' catalogues
# ----------------------------------------------------------------------------------------------------------------------------------
def pkey(cluster, pxy):
    return '%s|x%d_y%d' % (cluster, pxy[0], pxy[1])


def cs_filtered(b, key, f):
    """The crowdsource catalogue under `key` ('<cluster>|x<x0>_y<y0>|crowdsource|' prefix) with the post-filter f = (min flux / dflux, min fracflux, min qf, max |spread|): the rule of `exp_astro_m2.cs_dets`
    (flux > 0 and flux >= snr x dflux, fracflux >= ff, qf >= qf, |spread| <= sp); the default (no filter) is (0, 0, 0, inf)."""
    snr, ff, qf, sp = f
    m = (b[key + 'flux'] > 0) & (b[key + 'flux'] >= snr * b[key + 'dflux']) & (b[key + 'fracflux'] >= ff) & (b[key + 'qf'] >= qf) & (np.abs(b[key + 'spread']) <= sp)
    return b[key + 'x'][m], b[key + 'y'][m], b[key + 'flux'][m]


def baseline_tags(b, tuned):
    """The tags of the catalogue file `b` (np.load), checked against the frozen tuned choices: dict(tuned_phot, cs_tuned, fe_phot, fe_cs). Raises if the file's tuned tags are not the frozen ones."""
    meta = {k: str(b['meta|' + k][0]) for k in ('photutils_tuned', 'crowdsource_tuned', 'photutils_fe', 'crowdsource_fe')}
    if meta['photutils_tuned'] != tuned['photutils'] or meta['crowdsource_tuned'] != tuned['crowdsource']:
        raise RuntimeError('the tuned baselines of the catalogue file (%s, %s) are not the frozen ones %s' % (meta['photutils_tuned'], meta['crowdsource_tuned'], tuned))
    return dict(tuned_phot=meta['photutils_tuned'], cs_tuned=meta['crowdsource_tuned'], fe_phot=meta['photutils_fe'], fe_cs=meta['crowdsource_fe'])


def baseline_dets(b, cluster, pxy, tags):
    """{row label: (x, y, flux)} of the six baseline rows of a patch from the catalogue file `b`."""
    k = pkey(cluster, pxy) + '|'
    tri = lambda tag: (b[k + tag + '|x'], b[k + tag + '|y'], b[k + tag + '|flux'])
    return {'photutils (default)': tri('photutils'), 'photutils (tuned)': tri(tags['tuned_phot']), 'photutils (front end)': tri(tags['fe_phot']),
            'crowdsource (default)': cs_filtered(b, k + 'crowdsource|', (0.0, 0.0, 0.0, np.inf)), 'crowdsource (tuned)': cs_filtered(b, k + 'crowdsource|', BR.parse_cs_name(tags['cs_tuned'])),
            'crowdsource (front end)': tri(tags['fe_cs'])}


# ----------------------------------------------------------------------------------------------------------------------------------
# one patch
# ----------------------------------------------------------------------------------------------------------------------------------
def _store(arrays, prefix, dets, secs, n, info, truths, dk, trace):
    """Scores of the catalogues of one run (every kappa and the BIC key 0) against each truth of `truths` ({key prefix: truth patch}), the run time and entrants, the BIC block with the cap flag, the catalogue at
    the reported significance dk, and the record of the guarded path (counts, barred points, entry times)."""
    for tname, t in truths.items():
        for k, d in dets.items():
            arrays.update(X._pack('%s%s|k%g' % (tname, prefix, k), X.score_dets(d, t)))
    arrays['%s|secs' % prefix] = np.array([secs, n])
    arrays['%s|bic' % prefix] = np.array([info['bic_index'], info['n_blocks'], int(n >= X.MAX_ENTRIES)])
    arrays['%s|det' % prefix] = np.c_[dets[dk][0], dets[dk][1], dets[dk][2]]
    arrays['%s|trace|counts' % prefix] = trace['counts']
    arrays['%s|trace|barred' % prefix] = trace['barred']
    arrays['%s|trace|taus' % prefix] = trace['taus']


def _alias(arrays, src, dst):
    """The arrays stored under `src|...` and `reg|src|...` also under `dst|...` and `reg|dst|...` (the same run, used for two rows)."""
    for k in [k for k in arrays if k.startswith(src + '|') or k.startswith('reg|' + src + '|')]:
        if k.startswith('reg|'):
            arrays['reg|' + dst + k[len('reg|' + src):]] = arrays[k]
        else:
            arrays[dst + k[len(src):]] = arrays[k]


def _guarded_run(fn, *args):
    """One path through the guarded `imaging.sliding_path` (installed by the caller): the result of `fn` and the record of the single call of the guarded function."""
    if IM.sliding_path is not CF.guarded_traced:
        raise RuntimeError('the guarded path is not installed: no ISS run may use the unguarded path')
    del CF.TRACE[:]
    res = fn(*args)
    if len(CF.TRACE) != 1:
        raise RuntimeError('a run must call the guarded path exactly once, it called it %d times' % len(CF.TRACE))
    return res, CF.TRACE[0]


def patch_arrays(args):
    """The arrays of one fresh patch: stamp fit, the ISS with both kernels at both configurations (and the BIC stop) through the guarded path, the baselines, each scored against both truths. A function of its arguments only."""
    cluster, pxy, entry, shift, hp, rex, kap, tuned, base = args
    patch = F.load_patch(pxy[0], pxy[1], entry)
    truth = F.truth_arrays(entry)
    t = A.truth_for_patch(truth, pxy[0], pxy[1])
    t_reg = A.truth_for_patch(dict(px=truth['px'] + shift[0], py=truth['py'] + shift[1], m606=truth['m606']), pxy[0], pxy[1])
    truths = {TRUTHS[0]: t, TRUTHS[1]: t_reg}
    arrays, fit = PS.stamp_arrays(patch)
    arrays['fwhm'] = np.array([S2.psf_fwhm(patch['psf'])])
    arrays['cfg'] = np.array([hp, rex, kap])
    arrays['shift'] = np.array([shift[0], shift[1]])
    arrays['kernel|same'] = np.array([int(fit['kind'] == 'isotropic')])
    sig, parts, _ = IM.fit_psf_parts(patch['psf'], n_parts=3, sigma0=1.0)
    CF.install()
    try:
        for label, (h, r), dk in (('ours', (hp, rex), kap), ('default', (X.DEFAULT[0], X.DEFAULT[1]), X.DEFAULT[2])):
            (dets, secs, n, info), tr = _guarded_run(X.run_path, patch, h, r)
            _store(arrays, 'iso|' + label, dets, secs, n, info, truths, dk, tr)
            P_iso = IM.PixelPSF(*patch['dn'].shape, sig, parts=parts, hp=h or None)
            arrays['iso|%s|rsd' % label] = np.array([PS.residual_sd(patch, P_iso, np.c_[dets[dk][0], dets[dk][1], dets[dk][2]])])
            print('F patch %s %s %s iso: sigma_b %g r_ex %g, %d entries in %.0f s, BIC block %d of %d, guard activations %d' % (cluster, pxy, label, h, r, n, secs, info['bic_index'], info['n_blocks'], tr['counts'][6]), flush=True)
            if fit['kind'] == 'isotropic':
                _alias(arrays, 'iso|' + label, 'adq|' + label)
            else:
                (dets2, secs2, n2, info2), tr2 = _guarded_run(PS.run_path_psf, patch, h, r, fit)
                _store(arrays, 'adq|' + label, dets2, secs2, n2, info2, truths, dk, tr2)
                arrays['adq|%s|rsd' % label] = np.array([info2['rsd'][dk]])
                print('F patch %s %s %s adequate (%d parts): %d entries in %.0f s, BIC block %d of %d, guard activations %d' % (cluster, pxy, label, len(fit['parts5']), n2, secs2, info2['bic_index'], info2['n_blocks'], tr2['counts'][6]), flush=True)
    finally:
        CF.uninstall()
    b = np.load(base)
    tags = baseline_tags(b, tuned)
    for label, d in baseline_dets(b, cluster, pxy, tags).items():
        for tname, tt in truths.items():
            arrays.update(X._pack(tname + label, X.score_dets(d, tt)))
        arrays[label + '|n'] = np.array([len(d[0])])
    for k_, v_ in tags.items():
        arrays['tag|' + k_] = np.array([v_])
    return cluster, pxy, arrays


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
    """The hash of what decides a number of this test: the source of this script's run functions, of the functions of `exp_astro_m2`, `exp_astro_psf`, `astro_fresh` and `exp_astro_cyclefix` that run, guard and score
    the path, and the constants (patches, grids, default, kernel rule, rows). Edits of the report do not change it."""
    src = ''.join(inspect.getsource(fn) for fn in (load_frozen, patches, pkey, cs_filtered, baseline_tags, baseline_dets, _store, _alias, _guarded_run, patch_arrays, X.run_path, X.load_patch, X.score_dets, X._pack, X.f_measure,
                                                  PS.run_path_psf, PS.residual_sd, PS.stamp_arrays, PS.stamp_shape, F.load_patch, F.truth_arrays, F.patch_path, CF.guarded_traced, CF.trace_of, CF.install, CF.uninstall))
    consts = repr((PATCHES, X.DEFAULT, X.KAPPAS, X.LAT_H, X.MARGIN, X.FIT_RADIUS, X.MAX_ENTRIES, G.RMS_GOAL, G.K_MAX, X.DEV, LABELS, TRUTHS, BASE_ROWS, FILES))
    return hashlib.sha256((src + consts).encode()).hexdigest()


def freeze_record():
    """What is frozen: the configuration and baseline choices read from freeze.json and the catalogue file, the hashes of the code, the plan and every input (patch, truth and catalogue files)."""
    fr = load_frozen()
    reg = registry()
    data = {}
    for c, pxy in patches(reg):
        data[F.patch_path(pxy[0], pxy[1], reg[c])] = None
    for c in reg:
        data[os.path.join(reg[c]['root'], reg[c]['truth'])] = None
    data_files = {os.path.relpath(p, reg[next(iter(reg))]['root']): _sha256_file(p) for p in sorted(data)}
    with np.load(BASE) as b:
        tags = baseline_tags(b, fr['tuned'])
    return dict(config=dict(sigma_b=fr['hp'], r_ex=fr['rex'], kappa=fr['kap']), default=list(fr['default']), tuned_baselines=fr['tuned'], front_end_baselines=dict(photutils=tags['fe_phot'], crowdsource=tags['fe_cs']),
                files={f: _sha256_file(os.path.join(REPO, f)) for f in FILES}, data_files=data_files, frozen_code=frozen_code_hash(), patches=[[n, list(p)] for n, p in patches(reg)],
                shifts={'%s|%d,%d' % (c, p[0], p[1]): reg[c]['shifts'][reg[c]['patches'].index(list(p))] for c, p in patches(reg)})


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
        raise RuntimeError('no freeze record %s: run `exp_astro_fresh.py --freeze` first' % path)
    with open(path) as f:
        old = json.load(f)
    now = freeze_record()
    diff = [k for k in now if now[k] != old.get(k)]
    if diff:
        raise RuntimeError('the freeze record differs from the present code / inputs in: %s' % ', '.join(diff))
    return old


def part_run(out_path=RAW, workers=1, freeze_path=None, base=BASE):
    """Each fresh patch, once, at the frozen configuration; a patch whose result exists is not run again. Refuses unless the freeze record matches."""
    rec = check_freeze(freeze_path)
    fr = load_frozen()
    reg = registry()
    cur = X.load_raw(out_path) if os.path.exists(out_path) else {}
    todo = [(c, p) for c, p in patches(reg) if pre_of(c, p) + 'cfg' not in cur]
    shifts = {(c, tuple(p)): reg[c]['shifts'][reg[c]['patches'].index(list(p))] for c, p in patches(reg)}
    jobs = [(c, p, reg[c], shifts[(c, p)], fr['hp'], fr['rex'], fr['kap'], fr['tuned'], base) for c, p in todo]
    print('F: %d patches to run, %d exist' % (len(jobs), len(patches(reg)) - len(jobs)), flush=True)
    meta = {'meta|plan_sha256': np.array([_sha256_file(PLAN)]), 'meta|frozen_code': np.array([rec['frozen_code']]), 'meta|freeze_sha256': np.array([_sha256_file(freeze_path or FREEZE)]),
            'meta|n_clusters': np.array([len(reg)])}
    for c in reg:
        meta['title|' + c] = np.array([reg[c]['title']])
    if 'meta|plan_sha256' not in cur:
        X.save_part(out_path, 'meta', {k[len('meta|'):]: v for k, v in meta.items() if k.startswith('meta|')})
        X.save_part(out_path, 'title', {c: meta['title|' + c] for c in reg})
    if workers > 1 and len(jobs) > 1:
        from multiprocessing import Pool
        with Pool(min(workers, len(jobs))) as pool:
            for c, pxy, arrays in pool.imap_unordered(patch_arrays, jobs, chunksize=1):
                X.save_part(out_path, pre_of(c, pxy)[:-1], arrays)
    else:
        for j in jobs:
            c, pxy, arrays = patch_arrays(j)
            X.save_part(out_path, pre_of(c, pxy)[:-1], arrays)


# ----------------------------------------------------------------------------------------------------------------------------------
# report
# ----------------------------------------------------------------------------------------------------------------------------------
KERN = {'iso': 'isotropic kernel (as frozen)', 'adq': 'rule kernel (from the stamp)'}
CFG = {'ours': ('tuned (sigma_b 3, r_ex 1.2, kappa 2)', 'k2'), 'default': ('a priori default (5, 1.2, kappa 5)', 'k5')}


def rows_def():
    """[(label, key)] of the rows of every table: the ISS with the isotropic and the rule kernel at the three readouts, then the baselines."""
    out = []
    for kern in ('iso', 'adq'):
        out += [('ISS, %s, %s' % (KERN[kern], CFG['ours'][0]), '%s|ours|k2' % kern), ('ISS, %s, %s' % (KERN[kern], CFG['default'][0]), '%s|default|k5' % kern),
                ('ISS, %s, a priori default with the BIC stop' % KERN[kern], '%s|default|k0' % kern)]
    return out + [(l, l) for l in BASE_LABELS]


# primary and secondary comparisons: (id, name, first, second)
PRIMARY = [('O1', 'rule-kernel ISS (tuned) minus crowdsource (tuned)', 'adq|ours|k2', 'crowdsource (tuned)'),
           ('O2', 'rule-kernel ISS (tuned) minus crowdsource (front end)', 'adq|ours|k2', 'crowdsource (front end)'),
           ('O1i', 'isotropic ISS (tuned) minus crowdsource (tuned)', 'iso|ours|k2', 'crowdsource (tuned)'),
           ('O2i', 'isotropic ISS (tuned) minus crowdsource (front end)', 'iso|ours|k2', 'crowdsource (front end)')]
SECONDARY = [('rule-kernel ISS (tuned) minus photutils (tuned)', 'adq|ours|k2', 'photutils (tuned)'), ('rule-kernel ISS (tuned) minus photutils (front end)', 'adq|ours|k2', 'photutils (front end)'),
             ('isotropic ISS (tuned) minus photutils (tuned)', 'iso|ours|k2', 'photutils (tuned)'), ('isotropic ISS (tuned) minus photutils (front end)', 'iso|ours|k2', 'photutils (front end)'),
             ('rule-kernel ISS (default) minus crowdsource (default)', 'adq|default|k5', 'crowdsource (default)'), ('rule-kernel ISS (default) minus crowdsource (front end)', 'adq|default|k5', 'crowdsource (front end)'),
             ('rule-kernel ISS (default) minus photutils (default)', 'adq|default|k5', 'photutils (default)'),
             ('rule-kernel ISS (BIC stop) minus crowdsource (default)', 'adq|default|k0', 'crowdsource (default)'), ('rule-kernel ISS (BIC stop) minus crowdsource (front end)', 'adq|default|k0', 'crowdsource (front end)'),
             ('rule-kernel ISS (BIC stop) minus photutils (default)', 'adq|default|k0', 'photutils (default)'),
             ('rule kernel minus isotropic, tuned', 'adq|ours|k2', 'iso|ours|k2'), ('rule kernel minus isotropic, a priori default', 'adq|default|k5', 'iso|default|k5'),
             ('rule kernel minus isotropic, BIC stop', 'adq|default|k0', 'iso|default|k0')]
CONTEXT = [('M2 (4 test patches), isotropic kernel', '+0.118', 'win', '-'), ('M53 (4 patches)', '-0.058', 'loss', '+0.049 win'), ('M92 (2 patches)', '-0.215', 'loss', '-0.060 loss'),
           ('M2 repeat frames (6 patches)', '+0.054', 'win', '+0.063 win')]


def clusters_in(raw):
    """[(cluster, [patch, ...])] with a result in the raw, in name order."""
    out = {}
    for k in raw:
        if k.startswith('F|') and k.endswith('|cfg'):
            _, c, xy = k.split('|')[:3]
            x0, y0 = xy[1:].split('_y')
            out.setdefault(c, []).append((int(x0), int(y0)))
    return [(c, sorted(p)) for c, p in sorted(out.items())]


def blocks(raw, plist_c, key, tname=''):
    """The blocks of the patches (cluster, pxy) for the row `key` (truth `tname`) concatenated, and the sizes of the strata."""
    parts = [raw[pre_of(c, p) + tname + key + '|blk22'] for c, p in plist_c]
    return np.concatenate(parts), [len(q) for q in parts]


def has(raw, plist_c, key, tname=''):
    return all(pre_of(c, p) + tname + key + '|blk22' in raw for c, p in plist_c)


def tot_of(raw, plist_c, key, tname=''):
    return sum(raw[pre_of(c, p) + tname + key + '|tot'] for c, p in plist_c)


def purity_all(tot):
    """Matched detections over all detections of any magnitude, and the number of unmatched detections, from the totals [true22, matched22, det22, det22_matched, true_all, matched_all, det_all]."""
    return tot[5] / max(tot[6], 1), int(tot[6] - tot[5])


def method_table(raw, plist_c, defs, tname=''):
    """The rows `defs` pooled over the patches (stratified block bootstrap): matched / truth, completeness, purity, F, purity over all magnitudes. Returns (markdown table, blocks by key, resamples)."""
    blk, sizes = {}, None
    for label, key in defs:
        if has(raw, plist_c, key, tname):
            blk[key], sizes = blocks(raw, plist_c, key, tname)
    if not blk:
        return '', {}, None
    idx = X.boot_idx_strat(sizes)
    rows = []
    for label, key in defs:
        if key not in blk:
            continue
        tot = tot_of(raw, plist_c, key, tname)
        (c, cci), (pu, pci), (ff, fci) = X.boot_cpf(blk[key], idx)
        pa, un = purity_all(tot)
        rows.append([label, '%d / %d' % (tot[1], tot[0]), X.fmt_iv(c, cci), '%d / %d' % (tot[3], tot[2]), X.fmt_iv(pu, pci), X.fmt_iv(ff, fci), '%.3f (%d)' % (pa, un)])
    return X.table(['catalogue', 'matched / truth < 22', 'completeness', 'matched / detected < 22', 'purity', 'F', 'purity over all magnitudes (unmatched detections)'], rows), blk, idx


def diff_row(name, a, b, idx):
    dc, dp, df = X.boot_diff(a, b, idx)
    return [name, X.fmt_iv(dc[0], dc[1]), X.fmt_iv(dp[0], dp[1]), X.fmt_iv(df[0], df[1]), X.verdict(*df)], df


def pair_table(blk, idx, pairs):
    rows, dfs = [], {}
    for name, a, b in pairs:
        if a in blk and b in blk:
            r, df = diff_row(name, blk[a], blk[b], idx)
            rows.append(r)
            dfs[(a, b)] = df
    return X.table(['comparison', 'd completeness', 'd purity', 'd F', 'verdict on F'], rows), dfs


def cluster_verdicts(raw, clusters, tname=''):
    """{(id, cluster): (dF, interval, verdict)} of the primary comparisons per cluster (patches pooled, stratified bootstrap, seed 0)."""
    out = {}
    for c, pl in clusters:
        plist_c = [(c, p) for p in pl]
        for cid, name, a, b in PRIMARY:
            if has(raw, plist_c, a, tname) and has(raw, plist_c, b, tname):
                ba, sizes = blocks(raw, plist_c, a, tname)
                bb, _ = blocks(raw, plist_c, b, tname)
                df = X.boot_diff(ba, bb, X.boot_idx_strat(sizes))[2]
                out[(cid, c)] = (df[0], df[1], X.verdict(*df))
    return out


def reading(verdicts, n_win_text, none_text, mixed_text):
    """The declared reading of a list of per-cluster verdicts (fresh_plan.md section 4)."""
    n, w, l = len(verdicts), verdicts.count('win'), verdicts.count('loss')
    s = n_win_text if w == n else (none_text if w == 0 else mixed_text)
    if l:
        s += ' A loss at %d cluster%s is reported as a loss, whatever the other cluster shows.' % (l, '' if l == 1 else 's')
    return s


def primary_text(raw, clusters, tname=''):
    cv = cluster_verdicts(raw, clusters, tname)
    L = ['| outcome | ' + ' | '.join('%s (%d patch%s)' % (c, len(pl), '' if len(pl) == 1 else 'es') for c, pl in clusters) + ' | clusters won / tied / lost |', '|---|' + '---|' * (len(clusters) + 1)]
    for cid, name, a, b in PRIMARY:
        cells, vs = [], []
        for c, pl in clusters:
            if (cid, c) in cv:
                d, ci, v = cv[(cid, c)]
                cells.append('%s %s' % (X.fmt_iv(d, ci), v))
                vs.append(v)
            else:
                cells.append('-')
        L.append('| **%s** %s | %s | %d / %d / %d |' % (cid, name, ' | '.join(cells), vs.count('win'), vs.count('tie'), vs.count('loss')))
    out = '\n'.join(L)
    n = len(clusters)
    o2 = [cv[('O2', c)][2] for c, _ in clusters if ('O2', c) in cv]
    o1 = [cv[('O1', c)][2] for c, _ in clusters if ('O1', c) in cv]
    out += '\n\nReading of O2 (fresh_plan.md section 4): ' + reading(o2, 'the lead over the strongest baseline holds at all %d clusters never seen by any choice (n = %d clusters; one of them a single window).' % (n, n),
                                                                     'tie or loss at every cluster: the claim is restricted to the M2 field and seeing; against baselines that have the ISS\'s front end there is no demonstrated advantage outside it.',
                                                                     'the lead is not general: it is cluster-dependent; the per-cluster verdicts are the statement, with no general sentence.')
    out += '\n\nReading of O1: ' + reading(o1, 'the lead over the baseline as originally tuned holds at all %d fresh clusters.' % n, 'the lead over the baseline as originally tuned does not hold at any fresh cluster.',
                                              'the lead over the baseline as originally tuned holds at some fresh clusters only.')
    return out


def p1_text(raw, plist_all):
    """S1: on the patches where the rule picks a general kernel, the pooled paired dF of the rule kernel against the isotropic one."""
    gen = [(c, p) for c, p in plist_all if raw[pre_of(c, p) + 'kernel|same'][0] == 0]
    if not gen:
        return 'S1: the rule keeps the isotropic kernel on every patch; the two are the same program, dF = 0 by construction.'
    lines = ['S1 (kernel rule on a held-out cluster), the %d patches where the rule picks a general kernel (%s); M3\'s patch keeps the isotropic kernel (dF = 0 by construction) and does not enter:' % (len(gen), ', '.join('%s (%d, %d)' % (c, p[0], p[1]) for c, p in gen))]
    for label, a, b in (('tuned configuration (kappa 2)', 'adq|ours|k2', 'iso|ours|k2'), ('a priori default (kappa 5)', 'adq|default|k5', 'iso|default|k5'), ('BIC stop', 'adq|default|k0', 'iso|default|k0')):
        ba, sizes = blocks(raw, gen, a)
        bb, _ = blocks(raw, gen, b)
        df = X.boot_diff(ba, bb, X.boot_idx_strat(sizes))[2]
        v = {'win': 'confirmed on a held-out cluster', 'tie': 'not confirmed', 'loss': 'contradicted'}[X.verdict(*df)]
        lines.append('* %s: dF = %s, %s' % (label, X.fmt_iv(df[0], df[1]), v))
    return '\n'.join(lines)


MAG_EDGES = A.BIN_EDGES


def completeness_table(raw, plist_c, keys):
    """Completeness by magnitude bin (matched truth / truth) for the rows `keys` = [(label, key)], the patches pooled; bins with no truth left out."""
    nt = sum(raw[pre_of(c, p) + keys[0][1] + '|n_true'] for c, p in plist_c)
    rows = []
    for b in range(len(MAG_EDGES) - 1):
        if nt[b] == 0:
            continue
        rows.append(['%.1f-%.1f' % (MAG_EDGES[b], MAG_EDGES[b + 1]), '%d' % nt[b]] + ['%.2f' % (sum(raw[pre_of(c, p) + k + '|n_matched'][b] for c, p in plist_c) / nt[b]) for l, k in keys])
    return X.table(['mag bin', 'truth'] + [l for l, k in keys], rows)


def report(raw):
    clusters = clusters_in(raw)
    plist_all = [(c, p) for c, pl in clusters for p in pl]
    fr_cfg = raw[pre_of(*plist_all[0]) + 'cfg'] if plist_all else np.zeros(3)
    L = ['# Fresh-cluster test: the frozen configuration and the kernel rule at two globular clusters that no choice of the project has seen (SDSS r band against HST, M3 and M15)', '',
         'Generated by `experiments/exp_astro_fresh.py` from `fresh_raw.npz` alone; every number below is printed by the script. The analysis was declared in `fresh_plan.md` (SHA-256 stored in the raw: %s) before the '
         'baselines\' catalogues were built and before any ISS run; the freeze record `freeze_fresh.json` hashes the plan, the code and the inputs (SHA-256 of the record in the raw: %s), and the run refused to start unless it matched. '
         '**Nothing was re-tuned on these patches**: the configuration (sigma_b %g, r_ex %g, kappa %g), the a priori default (5, 1.2, 5), the BIC stop, the kernel rule (`psf_general.fit_adequate`, rms goal %g of the peak) and the baselines\' '
         'choices come from the M2 development patches and the earlier episodes. **The unit of replication is the CLUSTER: n = %d clusters (%s), %d patches in all; the patches of a cluster are windows of one frame, not independent '
         'trials, and the intervals are conditional on them. One cluster has a single patch.** Intervals: 95 %% percentile intervals of a bootstrap over the 20 x 20 px blocks (%d resamples, seed 0, stratified within each patch, the '
         'same resamples for every row of a table, so differences are paired). **Every ISS run in this report used the guarded path** (`imaging_guard.sliding_path_guarded` through `exp_astro_cyclefix.install`; the unguarded path has the '
         'limit-cycle defect of `cyclefix.md`); cap flags and guard activations are in the run records. The offset between truth and image is NOT corrected in the main tables (section "Registration-corrected truth" shows the corrected ones). '
         'Magnitudes from counts use M2\'s fixed zero point.' % (str(raw['meta|plan_sha256'][0])[:16] + '...' if 'meta|plan_sha256' in raw else 'n/a', str(raw['meta|freeze_sha256'][0])[:16] + '...' if 'meta|freeze_sha256' in raw else 'n/a',
                                                                  fr_cfg[0], fr_cfg[1], fr_cfg[2], G.RMS_GOAL, len(clusters), ', '.join('%s: %d' % (raw['title|' + c][0].split(',')[0] if 'title|' + c in raw else c, len(pl)) for c, pl in clusters),
                                                                  len(plist_all), X.BOOT), '']
    # the patches
    rows = []
    for c, p in plist_all:
        pre = pre_of(c, p)
        rms, kind, p5, sh = raw[pre + 'rms'], str(raw[pre + 'kind'][0]), raw[pre + 'parts5'], raw[pre + 'shape']
        rows.append([c, '(%d, %d)' % p, '%.2f' % raw[pre + 'fwhm'][0], '%d' % raw[pre + 'iso|ours|k2|tot'][0], '%.5f' % rms[0], kind, '%d' % (len(p5) if kind == 'general' else 3), '%.5f' % rms[1], '%.3f' % sh[0],
                     '(%+.2f, %+.2f)' % tuple(raw[pre + 'shift'])])
    L += ['## The clusters and patches, and their stamp fits (label-free)\n', 'Cluster, patch origin, FWHM of the psField stamp (px; the M2 development stamps have 2.15, M53 3.3 to 3.4, M92 2.9 to 3.0), truth sources brighter than 22 mag in the scored area, rms of the '
          'isotropic fit against the stamp\'s peak, the kernel the rule chose, its number of parts and rms, the ellipticity of the stamp, the rigid offset of the truth against the image (px added to the truth in the corrected scoring).\n',
          X.table(['cluster', 'patch', 'FWHM', 'truth < 22', 'isotropic rms', 'kernel chosen', 'parts', 'its rms', 'ellipticity', 'offset'], rows), '']
    for c, pl in clusters:
        if 'title|' + c in raw:
            L.append('* %s: patches %s' % (raw['title|' + c][0], ', '.join('(%d, %d)' % p for p in pl)))
    L.append('')
    # the primary outcomes
    L += ['## Primary outcomes O1 and O2 (declared in fresh_plan.md section 4), per cluster\n', 'dF = F of the first minus F of the second at mag < 22, the cluster\'s patches pooled; win = the 95 % interval above 0, loss = below, tie = contains 0. '
          'O1 and O2 use the rule-kernel row (the isotropic row where the rule keeps the isotropic kernel: M3); O1i and O2i the isotropic row (the program frozen on M2).\n', primary_text(raw, clusters), '',
          'For the record, the same difference to the front-end crowdsource at the earlier clusters (read from `baselines_hp2.md`, not recomputed): ' + '; '.join(
              '%s %s %s%s' % (n, d, v, (' (stamp kernel ' + s + ')') if s != '-' else '') for n, d, v, s in CONTEXT) + '.', '']
    # clusters
    for c, pl in clusters:
        plist_c = [(c, p) for p in pl]
        L += ['## %s\n' % (raw['title|' + c][0] if 'title|' + c in raw else c), 'The %d patch%s of the cluster pooled (%d blocks).\n' % (len(pl), '' if len(pl) == 1 else 'es', 25 * len(pl))]
        tbl, blk, idx = method_table(raw, plist_c, rows_def())
        L += [tbl, '', 'Paired differences (first minus second, same blocks), 95 % intervals:\n']
        pairs = [(cid + ' ' + name, a, b) for cid, name, a, b in PRIMARY] + [(n, a, b) for n, a, b in SECONDARY]
        ptbl, _ = pair_table(blk, idx, pairs)
        L += [ptbl, '']
        if len(pl) > 1:
            L += ['Record by patch (verdict on F of the primary comparisons, dF with its paired interval over the patch\'s own 25 blocks):\n']
            rows = []
            for cid, name, a, b in PRIMARY:
                cells = []
                for p in pl:
                    pre = pre_of(c, p)
                    idx1 = X.boot_idx(raw[pre + a + '|blk22'].shape[0])
                    df = X.boot_diff(raw[pre + a + '|blk22'], raw[pre + b + '|blk22'], idx1)[2]
                    cells.append('%s (%+.3f [%.3f, %.3f])' % (X.verdict(*df), df[0], df[1][0], df[1][1]))
                rows.append(['%s %s' % (cid, name)] + cells)
            L += [X.table(['comparison'] + ['(%d, %d)' % p for p in pl], rows), '']
        L += ['Completeness by magnitude (matched truth / truth) of the main rows:\n', completeness_table(raw, plist_c, [('ISS rule kernel, tuned', 'adq|ours|k2'), ('ISS isotropic, tuned', 'iso|ours|k2'), ('crowdsource (tuned)', 'crowdsource (tuned)'),
                                                                                                   ('crowdsource (front end)', 'crowdsource (front end)'), ('photutils (front end)', 'photutils (front end)')]), '']
    # S1
    L += ['## S1: the kernel rule on a held-out cluster\n', p1_text(raw, plist_all), '']
    # S3
    L += ['## S3: both clusters pooled (descriptive)\n', 'The %d patches pooled (%d blocks), weighted by patches (M3 one, M15 two); the stratified bootstrap keeps every patch\'s weight. n = %d clusters: no verdict rests on this table.\n' % (len(plist_all), 25 * len(plist_all), len(clusters))]
    tbl, blk, idx = method_table(raw, plist_all, rows_def())
    ptbl, _ = pair_table(blk, idx, [(cid + ' ' + name, a, b) for cid, name, a, b in PRIMARY])
    L += [tbl, '', 'Paired differences on the pooled blocks (first minus second):\n', ptbl, '']
    # S2
    if all(has(raw, [(c, p) for p in pl], 'adq|ours|k2', 'reg|') for c, pl in clusters):
        L += ['## S2: registration-corrected truth (secondary)\n', 'The same catalogues scored against the truth shifted by the rigid offset of each patch (fitted to the image and the truth before any catalogue existed; the offsets are in the first table). '
              'It never decides a verdict; if a verdict differs from the main table the difference is stated here.\n', primary_text(raw, clusters, 'reg|'), '']
        cvm, cvr = cluster_verdicts(raw, clusters), cluster_verdicts(raw, clusters, 'reg|')
        diffs = ['%s at %s: %s -> %s' % (cid, c, cvm[(cid, c)][2], cvr[(cid, c)][2]) for (cid, c) in sorted(cvm) if (cid, c) in cvr and cvm[(cid, c)][2] != cvr[(cid, c)][2]]
        L += ['Verdicts that change under the corrected truth: %s.' % ('; '.join(diffs) if diffs else 'none'), '']
        for c, pl in clusters:
            plist_c = [(c, p) for p in pl]
            L += ['### %s, corrected truth\n' % (raw['title|' + c][0] if 'title|' + c in raw else c)]
            tbl, blk, idx = method_table(raw, plist_c, rows_def(), 'reg|')
            L += [tbl, '']
    # P3
    L += ['## Residual scale of the catalogues (label-free)\n', 'Robust standard deviation of the standardised residual (1: the noise model and the kernel are right) at the reported significance; the isotropic and the rule kernel; "same" where the rule keeps the isotropic kernel.\n']
    rows = []
    for c, p in plist_all:
        pre = pre_of(c, p)
        same = raw[pre + 'kernel|same'][0] == 1
        rows.append([c, '(%d, %d)' % p, '%.2f' % raw[pre + 'iso|ours|rsd'][0], 'same' if same else '%.2f' % raw[pre + 'adq|ours|rsd'][0], '%.2f' % raw[pre + 'iso|default|rsd'][0], 'same' if same else '%.2f' % raw[pre + 'adq|default|rsd'][0]])
    L += [X.table(['cluster', 'patch', 'isotropic, tuned', 'rule kernel, tuned', 'isotropic, default', 'rule kernel, default'], rows), '']
    # run records
    L += ['## Run records (guarded path)\n', 'Entrants, run time of the path to its last entrant (one core, shared machine), the BIC block of the path, the entry cap (%d; a capped run is flagged and still reported), the activations of the guard (lattice points barred after two exits at one entry time), '
          'the exit events, the distinct lattice points that entered, the most entries of one point (5 or more: possible residual cycling that the two-strike guard does not remove) and the re-entries at the exit time.\n' % X.MAX_ENTRIES]
    rows = []
    for c, p in plist_all:
        pre = pre_of(c, p)
        for kern in ('iso', 'adq'):
            if kern == 'adq' and raw[pre + 'kernel|same'][0] == 1:
                continue
            for cf in LABELS:
                s_, b_, tc = raw[pre + '%s|%s|secs' % (kern, cf)], raw[pre + '%s|%s|bic' % (kern, cf)], raw[pre + '%s|%s|trace|counts' % (kern, cf)]
                rows.append([c, '(%d, %d)' % p, KERN[kern], cf, '%d' % s_[1], '%.0f' % s_[0], '%d of %d' % (b_[0], b_[1]), 'CAPPED' if b_[2] else '-', '%d' % tc[6], '%d' % tc[1], '%d' % tc[2], '%d' % tc[3], '%d' % tc[7]])
    L += [X.table(['cluster', 'patch', 'kernel', 'configuration', 'entrants', 'seconds', 'BIC block', 'cap', 'guard activations', 'exit events', 'distinct points', 'most entries of one point', 're-entries at the exit time'], rows), '']
    capped = [(c, p) for c, p in plist_all if any(raw[pre_of(c, p) + '%s|%s|bic' % (k, cf)][2] for k in ('iso', 'adq') for cf in LABELS)]
    acts = sum(int(raw[pre_of(c, p) + '%s|%s|trace|counts' % (k, cf)][6]) for c, p in plist_all for k in ('iso', 'adq') for cf in LABELS)
    L += ['Runs that reached the entry cap: %s. Guard activations in all runs: %d.' % ('none' if not capped else ', '.join('%s (%d, %d)' % (c, p[0], p[1]) for c, p in capped), acts), '']
    runs = [(c, p, k, cf) for c, p in plist_all for k in ('iso', 'adq') for cf in LABELS if pre_of(c, p) + '%s|%s|trace|counts' % (k, cf) in raw and not (k == 'adq' and raw[pre_of(c, p) + 'kernel|same'][0] == 1)]
    cyc = [(c, p, k, cf) for c, p, k, cf in runs if raw[pre_of(c, p) + '%s|%s|trace|counts' % (k, cf)][3] >= CF.CYCLING_FLAG]
    L += ['Runs in which one lattice point entered %d or more times (possible residual cycling that the two-strike guard does not remove): %d of %d%s.' % (
        CF.CYCLING_FLAG, len(cyc), len(runs), (' (' + '; '.join('%s (%d, %d), %s kernel, %s' % (c, p[0], p[1], 'isotropic' if k == 'iso' else 'rule', cf) for c, p, k, cf in cyc) + ')') if cyc else ''), '']
    return '\n'.join(L) + '\n'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--freeze', action='store_true')
    ap.add_argument('--run', action='store_true')
    ap.add_argument('--from-raw', action='store_true')
    ap.add_argument('--workers', type=int, default=int(os.environ.get('ASTRO_WORKERS', '1')))
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
