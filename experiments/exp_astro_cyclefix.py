"""
The astronomy paths that an enter/exit limit cycle of `imaging.sliding_path` made invalid: re-run with the guard, with a check of every other published path
===================================================================================================================================================================

Declared before any of it was run in `results/astro/cyclefix_plan.md` (read it first; the sha256 of the plan is in the raw). The defect (referee round REP, finding 2): `imaging.sliding_path` lets a lattice kernel
enter on the correlation at its lattice start and lets a unit leave on the correlation at its slid centre; when the block fit of a crowded cluster returns the entrant's coefficient 0 the unit leaves, its dual goes
back to the wall, and it re-enters at dt = 0, at the same entry time, until `max_entries` (2500). The time rule never acts (tau does not advance). The cap flag of the raw files (`bic = [block, hist length, 1]`) is
the signature of that cycle, not a censoring. The guard (`src/core/imaging_guard.sliding_path_guarded`): a lattice point that has exited twice at the same entry time is barred for the rest of the path.

This script imports the frozen modules read-only (`exp_astro_m2`, `exp_astro_psf`, `exp_astro_repeat`, `astro_repeat`, `astro_common`), replaces `imaging.sliding_path` in its own worker processes by a thin wrapper of the
guarded function (it keeps `barred` and the event counts in a side channel, because the frozen run functions drop the path's result), and calls the frozen run and scoring functions with the original configuration:

    python experiments/exp_astro_cyclefix.py --run [--workers N]         stage 1: the paths whose stored run hit the entry cap (11 by the `bic` flag, plus the paths whose stored entry count is the cap)
    python experiments/exp_astro_cyclefix.py --check [--only D,T,..] [--reverse]    check mode: every other published path, re-run with the guard and compared with its stored result (identical or not); families in
                                                                         the order given (default T, D, P, R, V); --reverse: the shortest stored runs first, for a second process on the same families
    python experiments/exp_astro_cyclefix.py --from-raw                  regenerate results/astro/cyclefix.md from cyclefix_raw.npz (reads the frozen raws, read only)
    python experiments/exp_astro_cyclefix.py --list                      the inventory of the stored path runs (what is affected, what is checked)

Raw `results/astro/cyclefix_raw.npz` (write-once: a run whose result exists is not run again, also when another process of this script saved it first; saving is atomic and under a lock; new file, the frozen raws
are never written). Key = `G|` (stage 1) or `C|` (check mode) + the original key of the run
(for example `G|D|x110_y630|hp5|rex0.8|k2|tot`), plus per run `|trace|counts`, `|trace|barred` (tau, b, a of every barred point), `|trace|taus` (the entry times; stage 1, and any run with an activation or at the cap)
and `|wall` (seconds). Environment: CYCLEFIX_OUT (directory of the raw and the md, default results/astro; for tests), CYCLEFIX_M2_RAW / CYCLEFIX_PSF_RAW / CYCLEFIX_REPEAT_RAW (the frozen raws, read only).
"""
import argparse
import collections
import fcntl
import hashlib
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
import exp_astro_m2 as X                  # noqa: E402
import exp_astro_psf as PS                # noqa: E402
import exp_astro_repeat as RE             # noqa: E402
from src.core import imaging as IM        # noqa: E402
from src.core import imaging_guard as G   # noqa: E402
from src.core import psf_general as PG    # noqa: E402

ASTRO = os.path.join(REPO, 'results', 'astro')
OUT = os.environ.get('CYCLEFIX_OUT', ASTRO)
RAW, MD = os.path.join(OUT, 'cyclefix_raw.npz'), os.path.join(OUT, 'cyclefix.md')
PLAN = os.environ.get('CYCLEFIX_PLAN', os.path.join(ASTRO, 'cyclefix_plan.md'))
M2RAW = os.environ.get('CYCLEFIX_M2_RAW', os.path.join(ASTRO, 'm2_raw.npz'))
PSFRAW = os.environ.get('CYCLEFIX_PSF_RAW', os.path.join(ASTRO, 'psf_raw.npz'))
REPRAW = os.environ.get('CYCLEFIX_REPEAT_RAW', os.path.join(ASTRO, 'repeat_raw.npz'))
GNS, CNS = 'G', 'C'                       # key namespaces: the guarded re-runs of the affected paths, the check-mode runs
LABEL = 'guarded re-run (post hoc correction of a defect found after the runs)'
ASRUN = 'as run (entry cap hit: an enter/exit cycle, invalid)'
CYCLING_FLAG = 5                          # a lattice point that entered this many times in one path is flagged as possible residual cycling
FROZEN = ('src/core/imaging.py', 'src/core/psf_general.py', 'experiments/astro_common.py', 'experiments/astro_second.py', 'experiments/astro_repeat.py', 'experiments/exp_astro_baselines.py',
          'experiments/exp_astro_m2.py', 'experiments/exp_astro_psf.py', 'experiments/exp_astro_repeat.py', 'results/astro/m2_raw.npz', 'results/astro/psf_raw.npz', 'results/astro/repeat_raw.npz',
          'results/astro/baseline_catalogs.npz', 'results/astro/baseline_catalogs_second.npz', 'results/astro/baseline_catalogs_repeat.npz', 'results/astro/repeat_plan.md',
          'results/astro/freeze.json', 'results/astro/freeze_repeat.json', 'data/astro/checksums_repeat.json')
HASHED = FROZEN + ('src/core/imaging_guard.py',)
COUNT_FIELDS = ('entries', 'exit events', 'distinct points entered', 'most entries of one point', 'points entered 3 or more times', 'blocks', 'guard activations', 're-entries at the exit time')
FAMILIES = ('D', 'T', 'V', 'P', 'R')
CHECK_ORDER = ('T', 'D', 'P', 'R', 'V')   # the priority of check mode (plan section 6): test patches, uncapped development paths, kernel study, repeat patches, design variants


def utc():
    return time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())


# ----------------------------------------------------------------------------------------------------------------------------------
# the guarded function in place of imaging.sliding_path, with a side channel
# ----------------------------------------------------------------------------------------------------------------------------------
ORIG_PATH = IM.sliding_path               # the unguarded function of the frozen module
TRACE = []                                # one record per call of the replacement (the frozen run functions call it exactly once per path)


def trace_of(res):
    """The side-channel record of a path result: counts (see COUNT_FIELDS), the barred points (tau, b, a) and the entry times."""
    ev = res['events']
    enters = [e for e in ev if e[2] == 'enter']
    per_point = collections.Counter(e[1] for e in enters)
    last_exit, again = {}, 0
    for tau, pid, kind in ev:
        if kind == 'exit':
            last_exit[pid] = tau
        elif pid in last_exit and abs(tau - last_exit[pid]) <= 1e-12 * max(1.0, abs(tau)):
            again += 1
    bars = np.array([[t, b, a] for t, (b, a) in res['barred']], float).reshape(-1, 3)
    taus = np.asarray(res['taus'], float)
    counts = np.array([len(taus), len(ev) - len(enters), len(per_point), max(per_point.values()) if per_point else 0, sum(1 for v in per_point.values() if v >= 3), len(res['hist']), len(bars), again])
    return dict(counts=counts, barred=bars, taus=taus)


def guarded_traced(*args, **kw):
    """`imaging_guard.sliding_path_guarded` (cycle_guard on) with the record of its result kept in TRACE."""
    res = G.sliding_path_guarded(*args, **kw)
    TRACE.append(trace_of(res))
    return res


def install():
    """The assignment `IM.sliding_path = <guarded function>` that every call site of the frozen modules sees (they look the function up at call time)."""
    IM.sliding_path = guarded_traced


def uninstall():
    IM.sliding_path = ORIG_PATH


# ----------------------------------------------------------------------------------------------------------------------------------
# raw files
# ----------------------------------------------------------------------------------------------------------------------------------
def load_raw(path):
    with np.load(path, allow_pickle=False) as z:
        return {k: z[k] for k in z.files}


def same_array(a, b):
    """Equal in shape, dtype kind and every element (NaN equal to NaN)."""
    return a.shape == b.shape and np.array_equal(a, b, equal_nan=a.dtype.kind == 'f' and b.dtype.kind == 'f')


def save_part(path, arrays):
    """Merge arrays into the raw file, atomically (written beside it and renamed) and under a lock (two processes of this script may save into the same raw); nothing already in it is changed."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path + '.lock', 'a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        old = load_raw(path) if os.path.exists(path) else {}
        clash = [k for k in arrays if k in old and not same_array(old[k], np.asarray(arrays[k]))]
        if clash:
            raise RuntimeError('write-once raw: %d keys exist with other values, for example %s' % (len(clash), clash[0]))
        old.update({k: np.asarray(v) for k, v in arrays.items()})
        tmp = '%s.tmp%d.npz' % (path, os.getpid())
        np.savez_compressed(tmp, **old)
        os.replace(tmp, path)


def frozen_raws():
    """The three frozen raws, read only: dict(m2, psf, rep)."""
    return dict(m2=load_raw(M2RAW), psf=load_raw(PSFRAW), rep=load_raw(REPRAW))


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def hashes(files=HASHED):
    return {f: sha256_file(os.path.join(REPO, f)) for f in files if os.path.exists(os.path.join(REPO, f))}


def meta_record():
    """What is recorded in the raw when the first result is saved: the hashes of the frozen files and of the guard, the hash of the plan, the versions."""
    import scipy
    hs = hashes()
    return {'meta|files': np.array(sorted(hs)), 'meta|sha256': np.array([hs[f] for f in sorted(hs)]), 'meta|plan_sha256': np.array([sha256_file(PLAN) if os.path.exists(PLAN) else 'absent']),
            'meta|versions': np.array(['python %d.%d.%d' % sys.version_info[:3], 'numpy ' + np.__version__, 'scipy ' + scipy.__version__])}


def preflight(out_path=None, freeze_checks=True):
    """Refuses to run unless the frozen modules still match their freeze records and the files hashed in the raw by earlier runs (and the plan) are unchanged. Returns (the current raw, what is to be
    recorded in it by this run: the whole record if the raw is new, else the hashes of the files that earlier runs did not record, if any)."""
    out_path = out_path or RAW
    if freeze_checks:
        X.check_freeze(load_raw(M2RAW))
        RE.check_freeze()
    cur = load_raw(out_path) if os.path.exists(out_path) else {}
    now = meta_record()
    if 'meta|sha256' not in cur:
        return cur, now
    now_h = dict(zip((str(f) for f in now['meta|files']), (str(h) for h in now['meta|sha256'])))
    rec = dict(zip((str(f) for f in cur['meta|files']), (str(h) for h in cur['meta|sha256'])))
    if 'meta|files_late' in cur:
        rec.update(zip((str(f) for f in cur['meta|files_late']), (str(h) for h in cur['meta|sha256_late'])))
    diff = sorted(f for f, h in rec.items() if now_h.get(f) != h)
    if diff or not np.array_equal(cur['meta|plan_sha256'], now['meta|plan_sha256']):
        raise RuntimeError('the frozen files, the guard or the plan differ from what the first run recorded: %s' % (', '.join(diff) or 'the plan'))
    late = [f for f in now_h if f not in rec]
    return cur, ({'meta|files_late': np.array(late), 'meta|sha256_late': np.array([now_h[f] for f in late])} if late else {})


# ----------------------------------------------------------------------------------------------------------------------------------
# the inventory of the stored path runs
# ----------------------------------------------------------------------------------------------------------------------------------
def frozen_config():
    """dict(hp, rex, kap, default) of the configuration frozen on the M2 development patches (freeze.json)."""
    fr = RE.load_frozen()
    return dict(hp=fr['hp'], rex=fr['rex'], kap=fr['kap'], default=tuple(fr['default']))


def test_groups():
    """[(data set, key format of the part, patches)] of the test patches stored in m2_raw."""
    return [('m2', 'T|x%d_y%d|', list(X.TESTS))] + [(ds, 'T2|%s|x%%d_y%%d|' % ds, list(v['patches'])) for ds, v in X.SECOND.items()]


def inventory(raws, cfg):
    """Every path run stored in the three raws, as dicts (rid = the key prefix of the run; fam D, T, V, P, R; the configuration h, r and the reported kappa dk; cost = the stored run time), in the
    order of the raws. A repeat run stored under both kernel names (the rule keeps the isotropic kernel) is listed once, under 'iso'."""
    m2, psf, rep = raws['m2'], raws['psf'], raws['rep']
    out = []

    def conf(label):
        return (cfg['hp'], cfg['rex'], cfg['kap']) if label == 'ours' else tuple(cfg['default'])

    for p in X.DEV:
        for hp in X.HPS:
            for rex in X.REXS:
                rid = 'D|x%d_y%d|hp%g|rex%g' % (p[0], p[1], hp, rex)
                if rid + '|bic' in m2:
                    out.append(dict(rid=rid, fam='D', raw='m2', pxy=tuple(p), h=float(hp), r=float(rex), dk=None))
    for ds, kfmt, plist in test_groups():
        for p in plist:
            for label in ('ours', 'default'):
                rid = (kfmt % tuple(p)) + label
                if rid + '|bic' in m2:
                    h, r, dk = conf(label)
                    out.append(dict(rid=rid, fam='T', raw='m2', ds=ds, pxy=tuple(p), label=label, h=float(h), r=float(r), dk=float(dk)))
    p0 = X.TESTS[0]
    for name in X.VARIANTS:
        rid = 'V|x%d_y%d|%s' % (p0[0], p0[1], name)
        if rid + '|secs' in m2:
            out.append(dict(rid=rid, fam='V', raw='m2', pxy=tuple(p0), variant=name, h=float(cfg['hp']), r=float(cfg['rex']), dk=None))
    for ds in PS.CLUSTERS:
        for p in X.SECOND[ds]['patches']:
            for label in ('ours', 'default'):
                rid = 'P|%s|x%d_y%d|%s' % (ds, p[0], p[1], label)
                if rid + '|bic' in psf:
                    h, r, dk = conf(label)
                    out.append(dict(rid=rid, fam='P', raw='psf', ds=ds, pxy=tuple(p), label=label, h=float(h), r=float(r), dk=float(dk),
                                    froz='T2|%s|x%d_y%d|%s|det' % (ds, p[0], p[1], label)))
    for name, p in RE.PATCHES:
        pre = RE.pre_of(name, p)
        if pre + 'kernel|same' not in rep:
            continue
        same = int(rep[pre + 'kernel|same'][0]) == 1
        for kern in ('iso', 'adq', 'scl'):
            if kern == 'adq' and same:
                continue
            for label in ('ours', 'default'):
                rid = pre + '%s|%s' % (kern, label)
                if rid + '|bic' in rep:
                    h, r, dk = conf(label)
                    s = float(rep[pre + 'scale'][0]) if kern == 'scl' else 1.0
                    out.append(dict(rid=rid, fam='R', raw='rep', name=name, pxy=tuple(p), kern=kern, label=label, h=float(h * s), r=float(r * s), dk=float(dk)))
    for s in out:
        z = raws[s['raw']]
        s['cost'] = float(z[s['rid'] + '|secs'][0]) if s['rid'] + '|secs' in z else 0.0
    return out


def capped(raws):
    """The stored path runs that hit the entry cap: {rid: how}, 'bic' where the stored `bic` flag (third element) is 1, 'entries' where there is no `bic` key and the stored entry count is the cap."""
    out = {}
    for z in raws.values():
        for k, v in z.items():
            if k.endswith('|bic') and len(v) >= 3 and int(v[2]) == 1:
                out[k[:-4]] = 'bic'
        for k, v in z.items():
            if k.endswith('|secs') and len(v) == 2 and int(v[1]) >= X.MAX_ENTRIES and k[:-5] + '|bic' not in z and k[:-5] not in out:
                out[k[:-5]] = 'entries'
    return out


def split_inventory(raws, cfg):
    """(affected specs, check specs, aliases): the runs of the inventory that hit the cap, the others, and the cap-flagged runs stored under the second kernel name (aliases of an isotropic run, not listed)."""
    inv = inventory(raws, cfg)
    cap = capped(raws)
    listed = {s['rid'] for s in inv}
    aff = [dict(s, how=cap[s['rid']]) for s in inv if s['rid'] in cap]
    chk = [s for s in inv if s['rid'] not in cap]
    return aff, chk, sorted(k for k in cap if k not in listed)


# ----------------------------------------------------------------------------------------------------------------------------------
# one path run (a function of its argument only, so that it can run in a worker process)
# ----------------------------------------------------------------------------------------------------------------------------------
_FROZEN_CACHE = {}


def _m2raw():
    if 'm2' not in _FROZEN_CACHE:
        _FROZEN_CACHE['m2'] = load_raw(M2RAW)
    return _FROZEN_CACHE['m2']


def store_run(arrays, dets, secs, n, info, t, dk, det=True, bic=True):
    """The arrays of a run in the layout of the frozen scripts, relative to the run's key: the scores of every readout (`exp_astro_m2._pack`), [seconds, entries], [BIC block, blocks, cap flag] and the
    catalogue at the reported kappa."""
    for k, d in dets.items():
        arrays.update(X._pack('k%g' % k, X.score_dets(d, t)))
    arrays['secs'] = np.array([secs, n])
    if bic:
        arrays['bic'] = np.array([info['bic_index'], info['n_blocks'], int(n >= X.MAX_ENTRIES)])
    if det:
        arrays['det'] = np.c_[dets[dk][0], dets[dk][1], dets[dk][2]]
    return arrays


def iso_rsd(patch, h, det):
    """The residual scale of a catalogue for the isotropic kernel of the patch's stamp at high-pass h (what `exp_astro_psf` and `exp_astro_repeat` store as `rsd` for the isotropic kernel)."""
    sig, parts, _ = IM.fit_psf_parts(patch['psf'], n_parts=3, sigma0=1.0)
    P_iso = IM.PixelPSF(*patch['dn'].shape, sig, parts=parts, hp=h or None)
    return PS.residual_sd(patch, P_iso, det)


def run_path_spec(spec):
    """The arrays (relative keys) of one path run, computed by the frozen run and scoring functions with the original configuration; the guard must be installed by the caller."""
    fam, pxy = spec['fam'], spec['pxy']
    arrays = {}
    if fam == 'D':
        patch, t = X.load_patch(*pxy), A.truth_for_patch(X.truth_arrays(), pxy[0], pxy[1])
        dets, secs, n, info = X.run_path(patch, spec['h'], spec['r'])
        store_run(arrays, dets, secs, n, info, t, None, det=False)
    elif fam == 'T':
        patch, t = X.load_patch(*pxy, ds=spec['ds']), A.truth_for_patch(X.truth_arrays(spec['ds']), pxy[0], pxy[1])
        dets, secs, n, info = X.run_path(patch, spec['h'], spec['r'])
        store_run(arrays, dets, secs, n, info, t, spec['dk'])
        arrays['rsd'] = np.array([iso_rsd(patch, spec['h'], arrays['det'])])               # extra (the frozen script stores none): the isotropic kernel's residual scale of the catalogue
    elif fam == 'V':
        patch, t = X.load_patch(*pxy), A.truth_for_patch(X.truth_arrays(), pxy[0], pxy[1])
        ov = dict(X.VARIANTS[spec['variant']])
        dets, secs, n, info = X.run_path(patch, ov.pop('hp', spec['h']), ov.pop('rex', spec['r']), **ov)
        store_run(arrays, dets, secs, n, info, t, None, det=False, bic=False)
        arrays['bic'] = np.array([info['bic_index'], info['n_blocks'], int(n >= X.MAX_ENTRIES)])      # extra (the frozen script stores none for the variants)
    elif fam == 'P':
        patch, t = X.load_patch(*pxy, ds=spec['ds']), A.truth_for_patch(X.truth_arrays(spec['ds']), pxy[0], pxy[1])
        fit = PG.fit_adequate(patch['psf'])
        dets, secs, n, info = PS.run_path_psf(patch, spec['h'], spec['r'], fit)
        store_run(arrays, dets, secs, n, info, t, spec['dk'])
        arrays['rsd'] = np.array([iso_rsd(patch, spec['h'], _m2raw()[spec['froz']]), info['rsd'][spec['dk']]])
    elif fam == 'R':
        name = spec['name']
        patch, t = R.load_patch(pxy[0], pxy[1], name), A.truth_for_patch(R.truth_arrays(name), pxy[0], pxy[1])
        if spec['kern'] == 'iso':
            dets, secs, n, info = X.run_path(patch, spec['h'], spec['r'])
            store_run(arrays, dets, secs, n, info, t, spec['dk'])
            arrays['rsd'] = np.array([iso_rsd(patch, spec['h'], arrays['det'])])
        else:                                                                               # the kernel the rule chose ('adq'), or the same with the lengths scaled by the seeing ('scl')
            fit = PG.fit_adequate(patch['psf'])
            dets, secs, n, info = PS.run_path_psf(patch, spec['h'], spec['r'], fit)
            store_run(arrays, dets, secs, n, info, t, spec['dk'])
            arrays['rsd'] = np.array([info['rsd'][spec['dk']]])
    else:
        raise ValueError('unknown family %r' % (fam,))
    return arrays


def run_job(job):
    """job = (spec, namespace). Installs the guard, runs the path, returns (key prefix of the arrays, arrays incl. the trace, wall seconds). A function of its argument only."""
    spec, ns = job
    install()
    TRACE.clear()
    t0 = time.time()
    print('%s start  %s' % (utc(), spec['rid']), flush=True)
    arrays = run_path_spec(spec)
    wall = time.time() - t0
    if len(TRACE) != 1:
        raise RuntimeError('%s: %d paths were run, expected exactly one' % (spec['rid'], len(TRACE)))
    tr = TRACE[0]
    arrays['trace|counts'] = tr['counts']
    arrays['trace|barred'] = tr['barred']
    if ns == GNS or len(tr['barred']) or tr['counts'][0] >= X.MAX_ENTRIES:
        arrays['trace|taus'] = tr['taus']
    arrays['wall'] = np.array([wall])
    c = tr['counts']
    print('%s done   %s: %d entries, %d guard activations, cap %s, %.0f s' % (utc(), spec['rid'], c[0], c[6], 'HIT' if c[0] >= X.MAX_ENTRIES else 'no', wall), flush=True)
    return ns + '|' + spec['rid'], arrays, wall


def run_jobs(specs, ns, out_path=None, workers=1, reverse=False):
    """Each spec once (a run whose result exists is not run again), results saved as they arrive; the longest stored runs first (the shortest first if `reverse`: a second process of this script working
    on the same paths meets the first one in the middle)."""
    out_path = out_path or RAW
    cur, now = preflight(out_path)
    todo = [s for s in specs if '%s|%s|wall' % (ns, s['rid']) not in cur]
    todo.sort(key=lambda s: -s.get('cost', 0.0), reverse=reverse)
    print('%s %s: %d paths to run, %d exist' % (utc(), ns, len(todo), len(specs) - len(todo)), flush=True)
    pending = dict(now)

    def save(prefix, arrays, wall):
        if os.path.exists(out_path) and prefix + '|wall' in load_raw(out_path):          # another process of this script finished the same path first: its result stands
            print('%s %s: already saved by another process' % (utc(), prefix), flush=True)
            return
        new = {prefix + '|' + k: v for k, v in arrays.items()}
        new.update(pending)
        pending.clear()
        save_part(out_path, new)

    if workers > 1 and len(todo) > 1:
        from multiprocessing import Pool
        with Pool(min(workers, len(todo))) as pool:
            for prefix, arrays, wall in pool.imap_unordered(run_job, [(s, ns) for s in todo], chunksize=1):
                save(prefix, arrays, wall)
    else:
        for s in todo:
            if os.path.exists(out_path) and '%s|%s|wall' % (ns, s['rid']) in load_raw(out_path):   # saved meanwhile by another process: not run again
                print('%s %s|%s: saved meanwhile by another process, not run' % (utc(), ns, s['rid']), flush=True)
                continue
            save(*run_job((s, ns)))
    return len(todo)


def stage1(workers=1, out_path=None):
    """The affected paths (the stored runs that hit the cap), re-run with the guard."""
    raws = frozen_raws()
    aff, chk, alias = split_inventory(raws, frozen_config())
    return run_jobs(aff, GNS, out_path, workers)


def check_mode(workers=1, only=None, out_path=None, reverse=False):
    """Every other stored path, re-run with the guard, in the order of CHECK_ORDER (or of `only`)."""
    raws = frozen_raws()
    aff, chk, alias = split_inventory(raws, frozen_config())
    order = list(only) if only else list(CHECK_ORDER)
    specs = []
    for fam in order:
        specs += [s for s in chk if s['fam'] == fam]
    cur, _ = preflight(out_path)
    tiers = []
    for fam in order:                                       # one tier per family, so that the priority order holds (a tier runs to its end before the next starts)
        tiers.append([s for s in specs if s['fam'] == fam])
    n = 0
    for tier in tiers:
        if tier:
            n += run_jobs(tier, CNS, out_path, workers, reverse)
    return n


# ----------------------------------------------------------------------------------------------------------------------------------
# comparison of a re-run with the stored run
# ----------------------------------------------------------------------------------------------------------------------------------
def stored_fields(old, rid):
    """The fields of the stored run `rid` in a raw: {relative key: array}."""
    pre = rid + '|'
    return {k[len(pre):]: v for k, v in old.items() if k.startswith(pre)}


def compare_run(new_fields, old_fields):
    """(identical, differing fields, missing fields, extra fields): the entry count (`secs[1]`), `bic` and every score array equal to the stored ones exactly (NaN equal to NaN); fields the frozen script did not
    store are extras (not compared); the run time is not compared."""
    diff, missing = [], []
    for f, ov in old_fields.items():
        if f == 'secs':
            nv = new_fields.get(f)
            if nv is None:
                missing.append(f)
            elif int(nv[1]) != int(ov[1]):
                diff.append(f)
            continue
        nv = new_fields.get(f)
        if nv is None:
            missing.append(f)
        elif not same_array(nv, ov):
            diff.append(f)
    extra = sorted(f for f in new_fields if f not in old_fields and not f.startswith('trace|') and f != 'wall')
    return (not diff and not missing), sorted(diff), sorted(missing), extra


def overlay(base, cf, ns=GNS):
    """A copy of the frozen raw `base` in which the keys of the re-runs of namespace ns replace the stored ones (trace keys and the wall time are not copied)."""
    out = dict(base)
    pre = ns + '|'
    for k, v in cf.items():
        if k.startswith(pre) and '|trace|' not in k and not k.endswith('|wall'):
            out[k[len(pre):]] = v
    return out


def alias_repeat(rep):
    """In a repeat raw, the isotropic runs of the patches where the rule keeps the isotropic kernel are also the chosen-kernel runs (`exp_astro_repeat._alias`): copy them over after a replacement."""
    out = dict(rep)
    for name, p in RE.PATCHES:
        pre = RE.pre_of(name, p)
        if pre + 'kernel|same' in out and int(out[pre + 'kernel|same'][0]) == 1:
            for lab in RE.LABELS:
                src = pre + 'iso|%s|' % lab
                for k in [k for k in out if k.startswith(src)]:
                    out[pre + 'adq|%s|' % lab + k[len(src):]] = out[k]
    return out


# ----------------------------------------------------------------------------------------------------------------------------------
# report
# ----------------------------------------------------------------------------------------------------------------------------------
RAWOF = {'D': 'm2', 'T': 'm2', 'T2': 'm2', 'V': 'm2', 'P': 'psf', 'R': 'rep'}
NAMES = {'ours': 'tuned configuration', 'default': 'a priori default'}
KERN = {'iso': 'isotropic kernel', 'adq': 'chosen kernel', 'scl': 'chosen kernel, lengths scaled by the seeing'}


def cfg_from_raw(m2):
    """The configuration stored with the test runs of the frozen m2 raw: dict(hp, rex, kap, default)."""
    for k in sorted(m2):
        if k.endswith('|cfg') and (k.startswith('T|') or k.startswith('T2|')):
            hp, rex, kap = (float(v) for v in m2[k])
            return dict(hp=hp, rex=rex, kap=kap, default=tuple(X.DEFAULT))
    raise ValueError('no test patch with a stored configuration in the m2 raw')


def _xy(s):
    a, b = s[1:].split('_y')
    return '%s, %s' % (a, b)


def describe(rid):
    """A short description of a stored path run from its key prefix."""
    p = rid.split('|')
    if p[0] == 'D':
        return 'development patch (%s), sigma_b %s, r_ex %s' % (_xy(p[1]), p[2][2:], p[3][3:])
    if p[0] == 'T':
        return 'M2 test patch (%s), %s' % (_xy(p[1]), NAMES[p[2]])
    if p[0] == 'T2':
        return '%s test patch (%s), %s' % (p[1].upper(), _xy(p[2]), NAMES[p[3]])
    if p[0] == 'V':
        return 'design variant "%s", M2 test patch (%s), tuned configuration' % (p[2], _xy(p[1]))
    if p[0] == 'P':
        return 'kernel study, %s patch (%s), %s' % (p[1].upper(), _xy(p[2]), NAMES[p[3]])
    return 'repeat frame %s patch (%s), %s, %s' % (p[1], _xy(p[2]), KERN[p[3]], NAMES[p[4]])


def kq(rid):
    """A run key as code in a table cell (the pipes of the key are escaped)."""
    return '`%s`' % rid.replace('|', '\\|')


def run_info(z, pre):
    """(entries, seconds, BIC block, blocks, capped) of the run stored under the key prefix `pre` (BIC block and blocks None where the frozen script stored no `bic`)."""
    s, b = z[pre + 'secs'], z.get(pre + 'bic')
    n = int(s[1])
    if b is None:
        return n, float(s[0]), None, None, n >= X.MAX_ENTRIES
    return n, float(s[0]), int(b[0]), int(b[1]), bool(b[2])


def f_txt(z, pre, k):
    """F (completeness, purity) at mag < 22 of the readout k of the run stored under `pre`, as text; '-' if the run has no such readout."""
    key = pre + 'k%g|tot' % k
    if key not in z:
        return '-'
    f, c, p = X.f_measure(z[key])
    return '%.3f (%.3f, %.3f)' % (f, c, p)


def readouts(rid, z, pre):
    """The readouts shown for a run: kappa 2, 5 and the BIC block for the development grid; the reported kappa and the BIC block otherwise (those the run has)."""
    p = rid.split('|')
    ks = (2.0, 5.0, 0.0) if p[0] == 'D' else ((5.0, 0.0) if 'default' in p else (2.0, 0.0))
    return [k for k in ks if pre + 'k%g|tot' % k in z]


def status_of(counts, bars):
    """How a guarded run ended, in words (plan section 3): cap not reached / cap reached without activation (censoring) / cap reached with activations; and a flag for a point entered many times."""
    n, act, most = int(counts[0]), int(counts[6]), int(counts[3])
    if n < X.MAX_ENTRIES:
        s = 'ended before the cap'
    elif act == 0:
        s = 'CAP REACHED, no activation: censoring of a long path'
    else:
        s = 'CAP REACHED with %d activations' % act
    if most >= CYCLING_FLAG:
        s += '; possible residual cycling (a point entered %d times)' % most
    return s


def bars_txt(bars, limit=12):
    if not len(bars):
        return '-'
    s = '; '.join('%.5g (%d, %d)' % (t, b, a) for t, b, a in bars[:limit])
    return s + ('; ... (%d in all)' % len(bars) if len(bars) > limit else '')


def have(cf, ns, rid):
    return '%s|%s|wall' % (ns, rid) in cf


def affected_section(fro, cf, aff):
    L = ['## 3. The affected paths: as run and guarded, side by side\n',
         'The paths whose stored run hit the entry cap (%d; the 11 flagged by the `bic` record and %d added by the scan of the stored entry counts), re-run with the guard and the original configuration. Entries count every entry, '
         're-entries included. "Activations" are the lattice points the guard barred (the entry time of their second exit at one entry time, then the lattice point (b, a)). **%s** for the guarded rows; **%s** '
         'for the stored rows.\n' % (len(aff), sum(1 for s in aff if s['how'] == 'entries'), LABEL, ASRUN)]
    rows, frows, notes = [], [], []
    for s in aff:
        rid, z = s['rid'], fro[RAWOF[s['rid'].split('|')[0]]]
        n0, t0, b0, nb0, c0 = run_info(z, rid + '|')
        head = '%s %s%s' % (describe(rid), kq(rid), '' if s['how'] == 'bic' else ' (added by the scan of the entry counts)')
        old_blocks = '-' if b0 is None else '%d of %d' % (b0, nb0)
        if not have(cf, GNS, rid):
            rows.append([head, '%d' % n0, old_blocks, '%.0f' % t0] + ['not re-run'] + [''] * 6)
            continue
        pre = '%s|%s|' % (GNS, rid)
        n1, t1, b1, nb1, c1 = run_info(cf, pre)
        cnt, bars = cf[pre + 'trace|counts'], cf[pre + 'trace|barred']
        rows.append([head, '%d' % n0, old_blocks, '%.0f' % t0, '%d' % n1, '%d of %d' % (b1, nb1), '%.0f' % t1, '%d' % cnt[6], '%d' % cnt[1], '%d' % cnt[2], '%d' % cnt[3], status_of(cnt, bars)])
        for k in sorted(set(readouts(rid, z, rid + '|')) | set(readouts(rid, cf, pre)), key=lambda k: (k == 0.0, k)):
            frows.append([describe(rid), 'kappa %g' % k if k else 'BIC block', f_txt(z, rid + '|', k), f_txt(cf, pre, k)])
        notes.append('* %s: guard activations (entry time, lattice point (b, a)): %s.' % (kq(rid), bars_txt(bars) if len(bars) else 'none'))
        if n1 >= X.MAX_ENTRIES:
            tau = cf[pre + 'trace|taus']
            q = np.percentile(tau, [0, 25, 50, 75, 100])
            notes.append('  Cap reached again. Entry times (of %d entries): min %.5g, quartiles %.5g, %.5g, %.5g, max %.5g.' % (len(tau), q[0], q[1], q[2], q[3], q[4]))
    L += [X.table(['path', 'as run: entries', 'as run: BIC block', 'as run: seconds', 'guarded: entries', 'guarded: BIC block', 'guarded: seconds', 'guard activations', 'exit events', 'distinct points entered',
                   'most entries of one point', 'how the guarded run ended'], rows), '']
    L += ['F at mag < 22 (completeness, purity) at the readouts, the stored run and the guarded re-run of the same path (kappa 2 is the frozen readout, kappa 5 the a priori default, BIC the BIC block):\n',
          X.table(['path', 'readout', 'as run', 'guarded re-run'], frows), '']
    L += notes + ['']
    chk_counts = [cf[k] for k in sorted(cf) if k.startswith(CNS + '|') and k.endswith('|trace|counts')]
    g_counts = [cf['%s|%s|trace|counts' % (GNS, s['rid'])] for s in aff if have(cf, GNS, s['rid'])]
    if chk_counts and g_counts:
        def stat(cs):
            most, frac = [int(c[3]) for c in cs], [c[1] / max(c[0], 1) for c in cs]
            return 'most entries of one point: median %g, maximum %d; exit events per entry: median %.3f, maximum %.3f' % (np.median(most), max(most), np.median(frac), max(frac))
        L += ['Added after the plan. For comparison with the guarded re-runs of the capped paths (%s), the %d check-mode paths run so far (%s). The guard removes only the repeats of a point at one entry time; a point that exits and enters '
              'again at different entry times is not touched, so these counts show how much of that behaviour is left in the guarded paths.\n' % (stat(g_counts), len(chk_counts), stat(chk_counts))]
    return L


def dev_section(fro, cf, aff):
    m2 = fro['m2']
    done = {s['rid'] for s in aff if s['fam'] == 'D' and have(cf, GNS, s['rid'])}
    ndev = sum(1 for s in aff if s['fam'] == 'D')
    corr = overlay(m2, cf, GNS)
    L = ['## 4. The development table recomputed, and the choice of the frozen configuration\n',
         'The twelve paths of the grid (two development patches x sigma_b 3, 5, 8 x r_ex 0.8, 1.2), the %d capped ones (%d of them re-run so far) replaced by their guarded re-runs, scored exactly as `exp_astro_m2.choose` '
         'does (largest F at mag < 22 of the two patches\' counts added, over (sigma_b, r_ex, kappa); ties on `round(F, 6)` go to the larger kappa; capped rows count as valid rows). The left block is the stored table '
         '(%s where a patch\'s path hit the cap), the right block is the table with the guarded re-runs in place (**%s**).\n' % (ndev, len(done), ASRUN.split(' (')[0], LABEL)]
    h0, r0, k0, f0 = X.choose(m2)
    h1, r1, k1, f1 = X.choose(corr)
    rows = []
    for h in X.HPS:
        for r in X.REXS:
            cap0 = sum(int(m2['D|x%d_y%d|hp%g|rex%g|bic' % (p[0], p[1], h, r)][2]) for p in X.DEV)
            cap1 = sum(int(corr['D|x%d_y%d|hp%g|rex%g|bic' % (p[0], p[1], h, r)][2]) for p in X.DEV)
            re_run = sum(1 for p in X.DEV if 'D|x%d_y%d|hp%g|rex%g' % (p[0], p[1], h, r) in done)
            for k in X.KAPPAS:
                t0 = sum(m2['D|x%d_y%d|hp%g|rex%g|k%g|tot' % (p[0], p[1], h, r, k)] for p in X.DEV)
                t1 = sum(corr['D|x%d_y%d|hp%g|rex%g|k%g|tot' % (p[0], p[1], h, r, k)] for p in X.DEV)
                (ff0, c0, p0), (ff1, c1, p1) = X.f_measure(t0), X.f_measure(t1)
                m0 = '**' if (h, r, k) == (h0, r0, k0) else ''
                m1 = '**' if (h, r, k) == (h1, r1, k1) else ''
                rows.append(['%g' % h, '%g' % r, '%g' % k, '%.3f' % c0, '%.3f' % p0, '%s%.3f%s' % (m0, ff0, m0), ('%d of %d' % (cap0, len(X.DEV))) if cap0 else '-',
                             '%.3f' % c1 if re_run else 'same', '%.3f' % p1 if re_run else 'same', ('%s%.3f%s' % (m1, ff1, m1)) if (re_run or m1) else 'same',
                             ('%d of %d' % (cap1, len(X.DEV))) if cap1 else '-'])
    L += [X.table(['sigma_b', 'r_ex', 'kappa', 'as run: completeness', 'as run: purity', 'as run: F', 'as run: capped patches', 'guarded: completeness', 'guarded: purity', 'guarded: F', 'guarded: capped patches'], rows), '',
          'Bold: the argmax of the table (`exp_astro_m2.choose`). Where "same" stands, the path was not re-run (not capped).\n']
    ranked = []
    for h in X.HPS:
        for r in X.REXS:
            cap = sum(int(corr['D|x%d_y%d|hp%g|rex%g|bic' % (p[0], p[1], h, r)][2]) for p in X.DEV)
            for k in X.KAPPAS:
                ff, c, pu = X.f_measure(sum(corr['D|x%d_y%d|hp%g|rex%g|k%g|tot' % (p[0], p[1], h, r, k)] for p in X.DEV))
                ranked.append((round(float(ff), 6), k, h, r, float(ff), cap))
    ranked.sort(reverse=True)                                                    # the order in which `choose` ranks the rows: F, then the larger kappa
    L += ['The ten best rows of the corrected table, in the order `choose` ranks them (F rounded to 6 decimals, then the larger kappa) (added after the plan):\n',
          X.table(['rank', 'sigma_b', 'r_ex', 'kappa', 'pooled F', 'capped patches'], [['%d' % (i + 1), '%g' % h, '%g' % r, '%g' % k, '%.4f' % ff, ('%d of %d' % (cap, len(X.DEV))) if cap else '-'] for i, (_, k, h, r, ff, cap) in enumerate(ranked[:10])]), '']
    cfg = cfg_from_raw(m2)
    frozen = (cfg['hp'], cfg['rex'], cfg['kap'])
    fz = X.f_measure(sum(corr['D|x%d_y%d|hp%g|rex%g|k%g|tot' % (p[0], p[1], *frozen)] for p in X.DEV))[0]
    L.append('* Argmax of the stored table: sigma_b %g, r_ex %g, kappa %g, pooled F %.4f (the frozen configuration).' % (h0, r0, k0, f0))
    L.append('* Argmax of the corrected table: sigma_b %g, r_ex %g, kappa %g, pooled F %.4f.' % (h1, r1, k1, f1))
    left = len(done) < ndev
    if left:
        out = 'INCOMPLETE: %d of the %d capped development paths are not re-run yet; no outcome is declared.' % (ndev - len(done), ndev)
    elif (h1, r1, k1) == frozen:
        out = '**Outcome (A)**: the argmax is the frozen configuration (sigma_b %g, r_ex %g, kappa %g); nothing else changes. Corrected pooled F of the frozen point: %.4f (stored: %.4f).' % (*frozen, fz, f0)
    else:
        out = ('**Outcome (B)**: another configuration wins: sigma_b %g, r_ex %g, kappa %g with pooled F %.4f, against %.4f for the frozen point (stored: %.4f). The plan declares that only now, as a labelled post hoc addition '
               'without any claim that it is held out, that configuration is run once with the guard on the 16 test patches and scored beside the as-run results; this report does not contain that run unless a section below says so.'
               % (h1, r1, k1, f1, fz, f0))
    L.append('* ' + out)
    # the robustness line: the argmax over the rows of the paths that did not reach the cap (a capped row is zeroed, so that `choose` itself decides)
    masked = dict(corr)
    for h in X.HPS:
        for r in X.REXS:
            for p in X.DEV:
                if int(corr['D|x%d_y%d|hp%g|rex%g|bic' % (p[0], p[1], h, r)][2]):
                    for k in X.KAPPAS:
                        key = 'D|x%d_y%d|hp%g|rex%g|k%g|tot' % (p[0], p[1], h, r, k)
                        masked[key] = np.zeros_like(corr[key])
    h2, r2, k2, f2 = X.choose(masked)
    L += ['* Robustness line (descriptive, not the outcome rule): the argmax over the rows of the paths that did not reach the cap, capped rows excluded for both patches of a configuration: sigma_b %g, r_ex %g, kappa %g, '
          'pooled F %.4f%s.' % (h2, r2, k2, f2, '' if (h2, r2, k2) == (h1, r1, k1) else ' (differs from the argmax of the full corrected table)'), '']
    return L


def dev_all_section(fro, cf):
    """Added after the plan: the development table with every one of the twelve paths a guarded run (the capped ones from stage 1, the others from check mode), because check mode shows that the guard acts
    also in paths that did not hit the cap."""
    m2 = fro['m2']
    dev_rids = ['D|x%d_y%d|hp%g|rex%g' % (p[0], p[1], h, r) for p in X.DEV for h in X.HPS for r in X.REXS]
    n_g = sum(1 for rid in dev_rids if have(cf, GNS, rid) or have(cf, CNS, rid))
    L = ['### 4b. The development table with all twelve paths guarded (added after the plan)\n',
         'Check mode shows that the guard acts also in paths whose run did not hit the cap (a double exit at one entry time that the unguarded path got over with a further entry), and that such a path differs from its stored run. '
         'The development table above replaces only the capped paths, as declared. This table takes every one of the twelve paths from a guarded run (the capped ones from stage 1, the others from check mode), '
         'so that all rows are rows of the same, guarded algorithm; it is scored by `exp_astro_m2.choose` as before. It is not the outcome rule of the plan.\n']
    if n_g < len(dev_rids):
        return L + ['INCOMPLETE: %d of the %d development paths have a guarded run so far.\n' % (n_g, len(dev_rids))]
    allg = overlay(overlay(m2, cf, GNS), {k: v for k, v in cf.items() if k.startswith(CNS + '|D|')}, CNS)
    h1, r1, k1, f1 = X.choose(allg)
    cfg = cfg_from_raw(m2)
    frozen = (cfg['hp'], cfg['rex'], cfg['kap'])
    ranked = []
    for h in X.HPS:
        for r in X.REXS:
            cap = sum(int(allg['D|x%d_y%d|hp%g|rex%g|bic' % (p[0], p[1], h, r)][2]) for p in X.DEV)
            for k in X.KAPPAS:
                ff = X.f_measure(sum(allg['D|x%d_y%d|hp%g|rex%g|k%g|tot' % (p[0], p[1], h, r, k)] for p in X.DEV))[0]
                ranked.append((round(float(ff), 6), k, h, r, float(ff), cap))
    ranked.sort(reverse=True)
    rank_frozen = 1 + [(h, r, k) for _, k, h, r, _f, _c in ranked].index(frozen)
    f_frozen = ranked[rank_frozen - 1][4]
    L += ['The ten best rows, in the order `choose` ranks them:\n',
          X.table(['rank', 'sigma_b', 'r_ex', 'kappa', 'pooled F', 'capped patches'], [['%d' % (i + 1), '%g' % h, '%g' % r, '%g' % k, '%.4f' % ff, ('%d of %d' % (cap, len(X.DEV))) if cap else '-'] for i, (_, k, h, r, ff, cap) in enumerate(ranked[:10])]), '']
    act = [(rid, int(cf['%s|%s|trace|counts' % (ns, rid)][6])) for rid in dev_rids for ns in (GNS, CNS) if have(cf, ns, rid)]
    L.append('* Guard activations in the twelve development paths: %s.' % ', '.join('%s %d' % (kq(rid), a) for rid, a in act))
    if (h1, r1, k1) == frozen:
        L.append('* Argmax of this table: sigma_b %g, r_ex %g, kappa %g, pooled F %.4f: **the frozen configuration is still the argmax** (rank 1; the second row has F %.4f).' % (h1, r1, k1, f1, ranked[1][4]))
    else:
        L.append('* Argmax of this table: sigma_b %g, r_ex %g, kappa %g, pooled F %.4f, **not the frozen configuration**, which has rank %d with F %.4f. The plan\'s outcome (A) / (B) is decided by the table of section 4, '
                 'not by this one; whether the configuration of this table is run on the test patches is for the coordinator to decide.' % (h1, r1, k1, f1, rank_frozen, f_frozen))
    L.append('')
    return L


def _blk(z, pre, key):
    return z[pre + key + '|blk22']


def _ci(v):
    return X.fmt_iv(v[0], v[1])


def m92_section(fro, cf):
    m2 = fro['m2']
    corr = overlay(m2, cf, GNS)
    plist = list(X.SECOND['m92']['patches'])
    L = ['### 5.1 M92 (NGC 6341): the a priori default and the BIC rows\n',
         'The two a priori default runs of M92 hit the entry cap (`T2|m92|*|default`); their guarded re-runs replace the stored arrays, everything else of the M92 report (the tuned run, the baselines) is the stored one. Intervals: '
         '95 %% percentile intervals of the block bootstrap of `exp_astro_m2.report` (%d resamples, seed 0, the same resamples for every row of a table, so differences are paired; pooled: stratified within each patch).\n' % X.BOOT]
    ok = all(have(cf, GNS, 'T2|m92|x%d_y%d|default' % tuple(p)) for p in plist)
    if not ok:
        return L + ['Not all M92 default paths are re-run yet.\n']
    keys = [('ISS a priori default (kappa 5)', 'default|k5'), ('ISS a priori default, BIC stop', 'default|k0')]
    rows, prow = [], []
    for scope, pl in [('(%d, %d)' % tuple(p), [p]) for p in plist] + [('pooled', plist)]:
        pre = ['T2|m92|x%d_y%d|' % tuple(p) for p in pl]
        for label, key in keys:
            cells = []
            for z in (m2, corr):
                blk = np.concatenate([_blk(z, q, key) for q in pre])
                idx = X.boot_idx(25) if len(pl) == 1 else X.boot_idx_strat([25] * len(pl))
                (c, cci), (pu, pci), (ff, fci) = X.boot_cpf(blk, idx)
                tot = sum(z[q + key + '|tot'] for q in pre)
                cells += [_ci((c, cci)), _ci((pu, pci)), _ci((ff, fci)), '%.3f (%d)' % (tot[5] / max(tot[6], 1), tot[6] - tot[5])]
            rows.append([scope, label] + cells)
        for ours_label, ours_key, base in (('default ISS', 'default|k5', ('photutils (default)', 'crowdsource (default)')), ('BIC-stop default ISS', 'default|k0', ('photutils (default)', 'crowdsource (default)', 'photutils (tuned)', 'crowdsource (tuned)'))):
            for b in base:
                cells = []
                for z in (m2, corr):
                    ba = np.concatenate([_blk(z, q, ours_key) for q in pre])
                    bb = np.concatenate([_blk(z, q, b) for q in pre])
                    idx = X.boot_idx(25) if len(pl) == 1 else X.boot_idx_strat([25] * len(pl))
                    df = X.boot_diff(ba, bb, idx)[2]
                    cells += [X.fmt_iv(df[0], df[1]), X.verdict(*df)]
                prow.append([scope, '%s vs %s' % (ours_label, b)] + cells)
    L += [X.table(['patch', 'row', 'as run: completeness', 'as run: purity', 'as run: F', 'as run: purity over all magnitudes (unmatched)', 'guarded: completeness', 'guarded: purity', 'guarded: F',
                   'guarded: purity over all magnitudes (unmatched)'], rows), '',
          'Paired differences (sliding ISS minus the baseline, same blocks), d F with its 95 % interval and the verdict on F:\n',
          X.table(['patch', 'comparison', 'as run: d F', 'as run: verdict', 'guarded: d F', 'guarded: verdict'], prow), '']
    for title, z in (('M92 report as `m2.md` prints it, stored arrays (%s)' % ASRUN.split(' (')[0], m2), ('M92 report as `m2.md` prints it, **%s**' % LABEL, corr)):
        sub = {k: v for k, v in z.items() if k.startswith('T2|m92|')}
        md = X.report(sub)
        body = md.split('\n## ', 1)[1]
        body = '## ' + body
        body = '\n'.join(('#### ' + l[3:]) if l.startswith('## ') else l for l in body.split('\n'))
        L += ['#### %s\n' % title, body.rstrip('\n'), '']
    return L


REPEAT_DEFAULT_ROWS = [('ISS, isotropic kernel, a priori default (5, 1.2, kappa 5)', 'iso|default|k5'), ('ISS, isotropic kernel, a priori default with the BIC stop', 'iso|default|k0'),
                       ('ISS, kernel chosen from the stamp, a priori default (5, 1.2, kappa 5)', 'adq|default|k5'), ('ISS, kernel chosen from the stamp, a priori default with the BIC stop', 'adq|default|k0'),
                       ('crowdsource (default)', 'crowdsource (default)'), ('photutils (default)', 'photutils (default)')]
REPEAT_PAIRS = [('adequate kernel vs isotropic, a priori default', 'adq|default|k5', 'iso|default|k5'), ('adequate kernel vs isotropic, BIC stop', 'adq|default|k0', 'iso|default|k0'),
                ('adequate kernel (a priori default) vs crowdsource (default)', 'adq|default|k5', 'crowdsource (default)'), ('adequate kernel (BIC stop) vs crowdsource (default)', 'adq|default|k0', 'crowdsource (default)')]
REPEAT_PAIRS_ADDED = [('isotropic kernel (a priori default) vs crowdsource (default)', 'iso|default|k5', 'crowdsource (default)'), ('isotropic kernel (BIC stop) vs crowdsource (default)', 'iso|default|k0', 'crowdsource (default)'),
                      ('adequate kernel (a priori default) vs photutils (default)', 'adq|default|k5', 'photutils (default)'), ('adequate kernel (BIC stop) vs photutils (default)', 'adq|default|k0', 'photutils (default)'),
                      ('isotropic kernel (a priori default) vs photutils (default)', 'iso|default|k5', 'photutils (default)'), ('isotropic kernel (BIC stop) vs photutils (default)', 'iso|default|k0', 'photutils (default)')]


def repeat_corrected(fro, cf):
    """The repeat raw with the guarded re-runs in place of the capped ones (the isotropic runs also under the chosen-kernel name where the rule keeps the isotropic kernel)."""
    return alias_repeat(overlay(fro['rep'], cf, GNS))


def _pair_tables(raw, plist, pairs):
    """Rows of the paired differences of `pairs` on the patches `plist` pooled (stratified bootstrap) for a repeat raw, with the blocks of the rows REPEAT_DEFAULT_ROWS and the tuned rows."""
    defs = REPEAT_DEFAULT_ROWS
    tbl, blk, idx = RE.method_table(raw, plist, defs)
    rows = []
    for name, a, b in pairs:
        if a in blk and b in blk:
            r, df = RE.diff_row(name, blk[a], blk[b], idx)
            rows.append(r)
    return tbl, rows


def p1_lines(raw, plist, B=X.BOOT):
    """P1's lines as `exp_astro_repeat.p1_text` computes them (dF of the chosen against the isotropic kernel, pooled over the patches where the rule picks a general kernel), with B resamples:
    [(readout, dF, interval, wording)]."""
    gen = [np_ for np_ in plist if raw[RE.pre_of(*np_) + 'kernel|same'][0] == 0]
    out = []
    for label, a, b in (('tuned configuration (kappa 2)', 'adq|ours|k2', 'iso|ours|k2'), ('a priori default (kappa 5)', 'adq|default|k5', 'iso|default|k5'), ('BIC stop', 'adq|default|k0', 'iso|default|k0')):
        ba, sizes = RE.pooled_blocks(raw, gen, a)
        bb, _ = RE.pooled_blocks(raw, gen, b)
        idx = X.boot_idx_strat(sizes, B=B)
        df = X.boot_diff(ba, bb, idx)[2]
        out.append((label, df[0], df[1], {'win': 'confirmed on held-out frames', 'tie': 'not confirmed', 'loss': 'contradicted'}[X.verdict(*df)]))
    return gen, out


def repeat_section(fro, cf):
    rep0 = fro['rep']
    rep1 = repeat_corrected(fro, cf)
    plist = [np_ for np_ in RE.done_patches(rep0)]
    L = ['### 5.2 The repeat frames: the a priori default and BIC rows\n',
         'The three isotropic a priori default runs that hit the cap (m2r2 (133, 281), m2r3 (370, 441), m2r3 (110, 501); the two m2r3 runs are also the chosen-kernel runs there, the rule keeps the isotropic kernel) are replaced by their guarded '
         're-runs. Pooled over the six patches and by frame (stratified block bootstrap, %d resamples, seed 0, the same resamples for every row of a table); the pairs are those `repeat.md` prints for these readouts '
         '(`PAIRS_MAIN`), and the pairs added after the plan are marked.\n' % X.BOOT]
    need = ['R|m2r2|x133_y281|iso|default', 'R|m2r3|x370_y441|iso|default', 'R|m2r3|x110_y501|iso|default']
    if not all(have(cf, GNS, rid) for rid in need):
        return L + ['Not all three repeat default paths are re-run yet.\n']
    scopes = [('All %d patches' % len(plist), plist)] + [('Frame %s' % fname, [np_ for np_ in plist if np_[0] == fname]) for fname in sorted(R.REPEAT)]
    for scope, pl in scopes:
        if not pl:
            continue
        for title, raw in (('%s' % ASRUN.split(' (')[0], rep0), ('**%s**' % LABEL, rep1)):
            tbl, rows = _pair_tables(raw, pl, REPEAT_PAIRS)
            _, rows_added = _pair_tables(raw, pl, REPEAT_PAIRS_ADDED)
            L += ['#### %s, %s\n' % (scope, title), tbl, '', X.table(['comparison', 'd completeness', 'd purity', 'd F', 'verdict on F'], rows + [[r[0] + ' (added after the plan)'] + r[1:] for r in rows_added]), '']
    L += ['#### Record by patch (verdict on F, paired intervals), the pairs of the default and BIC readouts\n']
    cells_rows = []
    for name_, a, b in REPEAT_PAIRS:
        for title, raw in (('as run', rep0), ('guarded', rep1)):
            cells = []
            for np_ in plist:
                pre = RE.pre_of(*np_)
                idx1 = X.boot_idx(raw[pre + a + '|blk22'].shape[0])
                df = X.boot_diff(raw[pre + a + '|blk22'], raw[pre + b + '|blk22'], idx1)[2]
                cells.append('%s (%+.3f)' % (X.verdict(*df), df[0]))
            cells_rows.append(['%s (%s)' % (name_, title)] + cells)
    L += [X.table(['comparison'] + ['%s (%d, %d)' % (n, p[0], p[1]) for n, p in plist], cells_rows), '']
    return L


def p1_section(fro, cf):
    rep0 = fro['rep']
    rep1 = repeat_corrected(fro, cf)
    plist = RE.done_patches(rep0)
    L = ['### 5.3 P1, the secondary lines, with the corrected isotropic arm on m2r2 (133, 281)\n']
    if not have(cf, GNS, 'R|m2r2|x133_y281|iso|default'):
        return L + ['The isotropic default path of m2r2 (133, 281) is not re-run yet.\n']
    gen, a2000 = p1_lines(rep0, plist)
    _, g2000 = p1_lines(rep1, plist)
    _, a40k = p1_lines(rep0, plist, B=40000)
    _, g40k = p1_lines(rep1, plist, B=40000)
    check = RE.p1_text(rep0, plist)
    same = all('dF = %s' % X.fmt_iv(d, ci) in check for _, d, ci, _w in a2000)
    L += ['P1 is the pooled dF = F(chosen kernel) - F(isotropic kernel) over the %d patches where the rule picks a general kernel (%s), 95 %% intervals of the block bootstrap with 2000 resamples (seed 0, as `repeat.md`) and, beside them, '
          'with 40000 resamples (seed 0). The wording is the plan of the repeat test\'s: a positive interval "confirmed on held-out frames", an interval containing 0 "not confirmed", a negative one "contradicted". The three arms of the '
          'stored table are reproduced by the function below from the stored raw: %s. Only the isotropic a priori default arm of m2r2 (133, 281) is capped; the tuned line does not change.\n' % (
              len(gen), ', '.join('%s (%d, %d)' % (n, p[0], p[1]) for n, p in gen), 'identical to the lines of `repeat.md`' if same else 'NOT identical to `repeat.md`: a defect of this report')]
    rows = []
    for (label, d0, ci0, w0), (_, d1, ci1, w1), (_, d2, ci2, w2), (_, d3, ci3, w3) in zip(a2000, g2000, a40k, g40k):
        rows.append([label, '%s, %s' % (X.fmt_iv(d0, ci0), w0), '%s, %s' % (X.fmt_iv(d1, ci1), w1), '%s, %s' % (X.fmt_iv(d2, ci2), w2), '%s, %s' % (X.fmt_iv(d3, ci3), w3)])
    L += [X.table(['readout', 'as run, 2000 resamples', 'guarded, 2000 resamples', 'as run, 40000 resamples', 'guarded, 40000 resamples'], rows), '',
          '(%s for the guarded columns.)\n' % LABEL]
    return L


def extras_section(fro, cf):
    L = ['### 5.4 Added after the plan: the residual scale of the corrected catalogues (label-free)\n']
    rows = []
    for rid, z0key in (('R|m2r2|x133_y281|iso|default', None), ('R|m2r3|x370_y441|iso|default', None), ('R|m2r3|x110_y501|iso|default', None)):
        if have(cf, GNS, rid):
            rows.append([describe(rid), '%.2f' % fro['rep'][rid + '|rsd'][0], '%.2f' % cf['%s|%s|rsd' % (GNS, rid)][0]])
    L += [X.table(['repeat default run (P3 of `repeat.md`, isotropic kernel, a priori default)', 'as run', 'guarded'], rows), ''] if rows else ['No repeat default run is re-run yet.\n']
    rows = []
    for p in X.SECOND['m92']['patches']:
        rid = 'T2|m92|x%d_y%d|default' % tuple(p)
        if have(cf, GNS, rid) and 'P|m92|x%d_y%d|default|rsd' % tuple(p) in fro['psf']:
            rows.append([describe(rid), '%.2f' % fro['psf']['P|m92|x%d_y%d|default|rsd' % tuple(p)][0], '%.2f' % cf['%s|%s|rsd' % (GNS, rid)][0]])
    if rows:
        L += ['`psf.md` prints the residual scale of the frozen isotropic catalogues next to the chosen kernel\'s; its frozen column for the M92 default is read from the capped catalogue:\n',
              X.table(['M92 default run (frozen isotropic kernel)', 'as run (psf.md)', 'guarded'], rows), '']
    return L


def psf_section(fro, cf):
    """`psf.md` compares the chosen kernel with the FROZEN rows of the same patches (the capped M92 default rows among them): its M92 section as `exp_astro_psf.report` prints it, with the stored and with the
    corrected frozen rows."""
    L = ['### 5.5 The kernel study (`psf.md`) against the corrected M92 default rows\n',
         'Added after the plan. `psf.md` compares the ISS with the chosen kernel (`psf_raw.npz`, not capped) with the frozen isotropic rows of the same patches; the frozen M92 a priori default and BIC rows are the capped ones. Its M92 sections, '
         'as `exp_astro_psf.report` prints them, with the stored frozen rows and with the guarded re-runs in place of the capped ones (the chosen-kernel rows are the same in both).\n']
    plist = list(X.SECOND['m92']['patches'])
    if not all(have(cf, GNS, 'T2|m92|x%d_y%d|default' % tuple(p)) for p in plist):
        return L + ['Not all M92 default paths are re-run yet.\n']
    corr = overlay(fro['m2'], cf, GNS)
    for title, froz in (('stored frozen rows (%s)' % ASRUN.split(' (')[0], fro['m2']), ('**%s**' % LABEL, corr)):
        md = PS.report(fro['psf'], froz)
        i = md.find('\n## M92')
        if i < 0:
            L += ['#### %s\n' % title, 'No M92 section in the kernel-study raw.\n']
            continue
        j = md.find('\n## ', i + 1)
        body = md[i + 1:j if j >= 0 else len(md)]
        body = '\n'.join(('#### ' + l[3:]) if l.startswith('## ') else ('##### ' + l[4:] if l.startswith('### ') else l) for l in body.split('\n'))
        L += ['#### %s\n' % title, body.rstrip('\n'), '']
    return L


def check_section(fro, cf, chk):
    L = ['## 6. Check mode: every other published path, re-run with the guard\n',
         'Each path below was re-run with the guard, the original configuration and the frozen functions; it is **identical** when its entry count, its BIC record and every score array of every readout, the stored catalogue '
         'and the stored residual scales are equal to the stored ones exactly (NaN equal to NaN); the run time is not compared. A path without any guard activation must be identical (the guard cannot act); a path with an '
         'activation had a double exit at one entry time and is expected to differ. Order of the families: test patches (T), uncapped development paths (D), kernel study (P), repeat patches (R), design variants (V).\n']
    rows, diffs = [], []
    counts = collections.OrderedDict((f, [0, 0, 0]) for f in CHECK_ORDER)
    chk = sorted(chk, key=lambda s: CHECK_ORDER.index(s['fam']))
    for s in chk:
        rid = s['rid']
        z0 = fro[s['raw']]
        if not have(cf, CNS, rid):
            counts[s['fam']][2] += 1
            rows.append([s['fam'], describe(rid) + ' ' + kq(rid), '-', '-', '-', 'not run'])
            continue
        old, new = stored_fields(z0, rid), stored_fields(cf, CNS + '|' + rid)
        ok, diff, missing, extra = compare_run(new, old)
        cnt = new['trace|counts']
        counts[s['fam']][0 if ok else 1] += 1
        rows.append([s['fam'], describe(rid) + ' ' + kq(rid), '%d' % int(old['secs'][1]), '%d' % int(new['secs'][1]), '%d' % cnt[6], 'identical' if ok else 'DIFFERENT'])
        if not ok:
            diffs.append((s, old, new, diff, missing, cnt))
    n_id = sum(c[0] for c in counts.values())
    n_df = sum(c[1] for c in counts.values())
    n_nr = sum(c[2] for c in counts.values())
    L += ['Summary: %d paths identical, %d different, %d not run (of %d). Per family (identical / different / not run): %s.\n' % (
        n_id, n_df, n_nr, len(chk), '; '.join('%s %d / %d / %d' % (f, *c) for f, c in counts.items()))]
    act = [(s['rid'], int(cf['%s|%s|trace|counts' % (CNS, s['rid'])][6])) for s in chk if have(cf, CNS, s['rid']) and int(cf['%s|%s|trace|counts' % (CNS, s['rid'])][6])]
    L.append('Paths with a guard activation: %s.\n' % (', '.join('%s (%d)' % (kq(a[0]), a[1]) for a in act) if act else 'none'))
    L += [X.table(['family', 'path', 'stored entries', 'guarded entries', 'guard activations', 'result'], rows), '']
    if diffs:
        L += ['### Paths that differ: F before and after\n']
        frows = []
        for s, old, new, diff, missing, cnt in diffs:
            rid = s['rid']
            for k in readouts(rid, new, ''):
                frows.append([describe(rid), 'kappa %g' % k if k else 'BIC block', f_txt(old, '', k), f_txt(new, '', k), '%d / %d' % (int(old['bic'][0]) if 'bic' in old else -1, int(old['bic'][1]) if 'bic' in old else -1),
                              '%d / %d' % (int(new['bic'][0]), int(new['bic'][1])), ', '.join(diff[:4]) + (' ...' if len(diff) > 4 else '') + (' missing: ' + ', '.join(missing) if missing else '')])
        L += [X.table(['path', 'readout', 'F before (completeness, purity)', 'F after', 'BIC block / blocks before', 'after', 'fields that differ'], frows), '']
    return L


def static_defect():
    return ['## 1. The defect and the guard\n',
            '`imaging.sliding_path` lets a lattice kernel enter when its dual reaches the wall; the entry test uses the correlation at the kernel\'s *lattice start* (`v`). A unit leaves when its coefficient is 0 and the correlation at its '
            '*slid centre* (`g`) is not positive; its dual then goes back to the wall (`U = rho`), so a kernel with `v > 0` re-enters at `dt = 0`, at the same entry time. If the block fit of the crowded cluster returns its coefficient '
            '0 each time, the path loops at one tau until `max_entries` (%d); the time rule (`tau_max`) never acts because tau does not advance. Such a run stops after about 150 distinct entries whatever its stopping rule says, and '
            'the cap flag of the raw files (`bic = [block, hist length, 1]`) is the signature of the cycle, not of a long path (found by referee round REP, finding 2; repeat patch m2r3 (110, 501), a priori default: the lattice point '
            '(110, 54) enters and exits 2345 times at tau = 2.18969).\n' % X.MAX_ENTRIES,
            'The guard (`src/core/imaging_guard.sliding_path_guarded`): a lattice point that has exited twice at the same entry time (the two times within 1e-12 x max(1, tau): the repeats of a cycle are at exactly the same tau, since the re-entry is at dt = 0) is barred from entering for the rest of the path; the result lists the barred '
            'points (`barred`). A first exit followed by a successful re-entry is untouched, so a path without a double exit at one entry time is the unguarded path bit for bit, and `cycle_guard=False` is `imaging.sliding_path` bit for bit '
            '(tests: `tests/test_imaging_guard.py`). The guard removes the repeated entry, not its cause (the entry test at the lattice start and the exit test at the slid centre disagree). The re-runs below call the frozen run functions '
            '(`exp_astro_m2.run_path`, `exp_astro_psf.run_path_psf`) with `imaging.sliding_path` replaced, in the worker process, by a wrapper of the guarded function that keeps the barred points and the event counts in a side channel.\n']


def frozen_section(cf):
    L = ['## 2. Frozen state\n']
    if 'meta|sha256' not in cf:
        return L + ['No run yet.\n']
    L += ['Recorded by the first run (sha256; the frozen modules were imported read only, and the run refuses to start if any differs from what the first run recorded; the freeze records `freeze.json` and `freeze_repeat.json` were '
          'checked by `exp_astro_m2.check_freeze` and `exp_astro_repeat.check_freeze` at every start). Plan `cyclefix_plan.md`: sha256 %s. Versions: %s.\n' % (cf['meta|plan_sha256'][0], ', '.join(str(v) for v in cf['meta|versions'])),
          X.table(['file', 'sha256'], [[str(f), str(h)] for f, h in zip(cf['meta|files'], cf['meta|sha256'])]), '']
    if 'meta|files_late' in cf:
        L += ['Files that the first run did not record (a file that `freeze_repeat.json` hashes was added to the list afterwards), recorded by a later run:\n',
              X.table(['file', 'sha256'], [[str(f), str(h)] for f, h in zip(cf['meta|files_late'], cf['meta|sha256_late'])]), '']
    return L


def limitations(fro, cf, aff, chk):
    n_run = sum(1 for s in chk if have(cf, CNS, s['rid']))
    nr = [s['rid'] for s in chk if not have(cf, CNS, s['rid'])]
    L = ['## 7. Limitations\n',
         '* **The corrected rows are post hoc.** The defect was found by a referee after the results had been seen; the guard was written after that; every corrected row says so and none is held out.',
         '* **The guard changes the path, it does not restore "what the unguarded algorithm would have done".** After a lattice point is barred, its neighbours enter in its place; the corrected rows are the rows of the guarded algorithm. '
         'The guard removes repeated entries at one entry time; a cycle whose turns advance tau, or a point that exits and re-enters at different times, is not removed. The largest number of entries of one lattice point and the number of points '
         'entered three times or more are printed for every guarded re-run (table in section 3); a point entered %d times or more is flagged.' % CYCLING_FLAG,
         '* **The guard keys on the entry time, so it never acts in `order=\'greedy\'`** (tau is the step count there); the design variant `greedy_snr` is the one greedy astronomy path.',
         '* **A cap reached again is not "the" path.** Where a guarded re-run reaches %d entries again its rows are cut by the cap; the report says whether it did so with or without activations (section 3), and the development table '
         'flags such rows in its capped column.' % X.MAX_ENTRIES,
         '* **Check mode: %d of %d paths were run%s.** %s' % (n_run, len(chk), '' if n_run == len(chk) else ', %d were not' % (len(chk) - n_run),
                                                              ('Not run: ' + ', '.join('`%s`' % r for r in nr)) if nr else 'All of them were run.'),
         '* **Identity means identity in this environment** (same Python, numpy, scipy and machine type as the stored runs; the stored runs were reproduced bit for bit wherever the guard did not act). The run times are those of a shared, loaded machine and '
         'are not comparable with the stored ones.',
         '* **Added by the scan, not by the brief:** the design variant `norex` stored 2500 entries but no `bic` record; it is re-run in stage 1 and labelled as an addition.',
         '* **Derived tables.** `psf.md` compares the adequate kernel with the *frozen* M92 default rows and prints their residual scale; those rows are the capped ones (section 5.4 gives the corrected residual scale and section 5.5 reprints the M92 part of `psf.md` '
         'with the corrected rows; the other tables of `psf.md` are not regenerated here). `m2.md` and `repeat.md` are not regenerated; the coordinator decides how the first outputs are kept.',
         '* **Intervals** are conditional on the patches (blocks of 20 x 20 px resampled inside fixed patches), as in the frozen reports; windows and truth are shared between the frames.', '']
    return L


def report(cf, fro):
    """The markdown report from the cyclefix raw `cf` and the frozen raws `fro` = dict(m2, psf, rep) (read only)."""
    cfg = cfg_from_raw(fro['m2'])
    aff, chk, alias = split_inventory(fro, cfg)
    L = ['# The astronomy paths that an enter/exit limit cycle of `imaging.sliding_path` made invalid: guarded re-runs and a check of every other path', '',
         'Generated by `experiments/exp_astro_cyclefix.py` from `cyclefix_raw.npz` and the frozen `m2_raw.npz`, `psf_raw.npz` and `repeat_raw.npz` (read only); every number below is printed by the script. Declared before any re-run in '
         '`cyclefix_plan.md`. **%s** is the label of every corrected row; the stored rows are not changed or replaced anywhere (**%s**). Frozen configuration: sigma_b %g, r_ex %g, kappa %g; a priori default %g, %g, %g.' % (
             LABEL, ASRUN, cfg['hp'], cfg['rex'], cfg['kap'], *cfg['default']), '']
    L += static_defect() + frozen_section(cf)
    L += affected_section(fro, cf, aff)
    L += dev_section(fro, cf, aff)
    L += dev_all_section(fro, cf)
    L += ['## 5. The corrected rows that the paper uses (**%s**)\n' % LABEL]
    L += m92_section(fro, cf) + repeat_section(fro, cf) + p1_section(fro, cf) + extras_section(fro, cf) + psf_section(fro, cf)
    L += check_section(fro, cf, chk)
    L += limitations(fro, cf, aff, chk)
    return '\n'.join(L) + '\n'


def list_inventory():
    raws = frozen_raws()
    aff, chk, alias = split_inventory(raws, frozen_config())
    print('stage 1 (%d paths: the stored run hit the entry cap):' % len(aff))
    for s in aff:
        print('  %-48s flagged by %-7s stored run time %6.0f s' % (s['rid'], s['how'], s['cost']))
    print('aliases of an isotropic run (not run): %s' % ', '.join(alias))
    print('check mode (%d paths), in the order %s:' % (len(chk), ','.join(CHECK_ORDER)))
    for fam in CHECK_ORDER:
        sub = [s for s in chk if s['fam'] == fam]
        print('  %s: %d paths, stored run time %.0f s' % (fam, len(sub), sum(s['cost'] for s in sub)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', action='store_true', help='stage 1: re-run the paths whose stored run hit the entry cap, with the guard')
    ap.add_argument('--check', action='store_true', help='check mode: re-run every other stored path with the guard and compare')
    ap.add_argument('--from-raw', action='store_true')
    ap.add_argument('--list', action='store_true')
    ap.add_argument('--only', default=None, help='check mode: families in the order to run them (letters of D,T,V,P,R)')
    ap.add_argument('--reverse', action='store_true', help='check mode: the shortest stored runs first (a second process meets the first one in the middle)')
    ap.add_argument('--workers', type=int, default=int(os.environ.get('ASTRO_WORKERS', '2')))
    A_ = ap.parse_args()
    if A_.list:
        list_inventory()
    elif A_.run:
        stage1(A_.workers)
    elif A_.check:
        check_mode(A_.workers, only=A_.only.split(',') if A_.only else None, reverse=A_.reverse)
    elif A_.from_raw:
        with open(MD, 'w') as f:
            f.write(report(load_raw(RAW), frozen_raws()))
        print('wrote', MD)
    else:
        ap.error('give --run, --check, --from-raw or --list')


if __name__ == '__main__':
    main()
