"""
Two POST HOC controls of the astronomy results: the grid edge of the development optimum, and the question "is the sliding ISS just a greedy fitter?"
===================================================================================================================================================

Everything below is POST HOC: it answers reviewer objections after the test numbers of the paper were seen. Nothing here changes the frozen configuration (`results/astro/freeze.json`,
`freeze_repeat.json`), a frozen file, a stored raw or an existing report; the frozen modules are imported read only and the scripts' own freeze checks run before every run. The
sliding path is always the guarded one (`src/core/imaging_guard.sliding_path_guarded`, put in place of `imaging.sliding_path` in the worker process exactly as `experiments/exp_astro_cyclefix.py`
does it: the defect and the guard are described there); the cap flag (2500 entries) and the guard activations of every path are recorded and printed. For `order='greedy'` the guard
never acts (tau is the step count, no two exits share an entry time: `imaging_guard`'s docstring), so the greedy paths are the unguarded paths, and the report says so with the (zero) activations.

CONTROL 1 -- the grid edge (DEVELOPMENT PATCHES ONLY)
-----------------------------------------------------
The development optimum (sigma_b 3, r_ex 1.2, kappa 2) lies on the edge of every axis of the development grid of `exp_astro_m2` (sigma_b {3, 5, 8}, r_ex {0.8, 1.2}, kappa {2, 2.5, ..., 8}):
the lowest sigma_b, the highest r_ex, the lowest kappa. The control extends the grid beyond each edge on the TWO DEVELOPMENT PATCHES (110, 630) and (210, 630) only (a test patch is never loaded
or scored by this control; `run_spec` refuses it): sigma_b in {2, 2.5} (new), r_ex in {1.6, 2.0} (new), kappa in {1.25, 1.5, 1.75} (new). The kappa axis costs no extra path: one path run to the
kappa = 1.25 time carries the readouts of every larger kappa (the shorter stopping time is a prefix of the same path; a test and the report compare the readouts of kappa >= 2 of the
re-run old paths with the guarded development rows of `cyclefix_raw.npz`, bit for bit). The set of paths, per development patch, one run each (13 paths, 26 in all; 9 are new
(sigma_b, r_ex) combinations and 4 are old combinations re-run to the kappa = 1.25 time so that the kappa axis can be extended there), every path read at the 11 kappas
{1.25, 1.5, 1.75, 2, 2.5, 3, 3.5, 4, 5, 6, 8}:

    (sigma_b, r_ex) in  {2, 2.5, 3} x {1.2, 1.6, 2.0}      the three lowest sigma_b against the three highest r_ex (9 paths; (3, 1.2) is the frozen combination, re-run)
                        (3, 0.8)                           the old lower r_ex at the frozen sigma_b (kappa extension only)
                        (5, 1.2), (5, 1.6)                 the middle sigma_b at the frozen r_ex and at the next r_ex
                        (8, 1.2)                           the widest sigma_b at the frozen r_ex (kappa extension only)

Not crossed (cost): r_ex 0.8 at sigma_b 2 and 2.5, r_ex 2.0 at sigma_b 5 and 8, r_ex 1.6 at sigma_b 8. The extended surface is the union of these paths and the old guarded development rows of
`cyclefix_raw.npz` (the six capped old paths from their guarded re-runs `G|`, the other six from check mode `C|`: all twelve are guarded paths); for the old combinations that were re-run the old rows
and the new rows of kappa >= 2 are the same numbers. Scoring is exactly `exp_astro_m2.choose`: the pooled F-measure of completeness and purity at mag < 22 (the two development patches' counts added),
the argmax over (sigma_b, r_ex, kappa), ties on round(F, 6) to the larger kappa (a test checks the function against `exp_astro_m2.choose` on the stored development rows). A second surface with the old
rows as stored (unguarded, `m2_raw.npz`; kappa >= 2) in place of the guarded ones is argmax-ed in the same way. Reported: whether the argmax moves, by how much F, the F surface around the frozen point
(axis profiles through it, the full surface, the ranking), the entry count, cap flag and guard activations of every path. Reading rule (fixed here): an extended point that beats the frozen point by more than
0.01 pooled development F, against the guarded and the stored value of the frozen point alike, is a FINDING and is stated as such; it is NOT run on any test patch. A row from a path that hit the
cap is flagged (the cap censors the path at its end; the flag is per path, so it is conservative for the large kappas). The maximum over a larger grid on the same two patches is optimistically
biased; the report says so and does not claim a held-out value.

CONTROL 2 -- is it just a greedy fitter? (all astronomy test patches the paper uses for the tuned ISS)
--------------------------------------------------------------------------------------------------------
`order='greedy'` is the same code as the dual order: sliding matching pursuit with the same PSF mixture, the same high-pass, the same weights, the same block fits, exits, exclusion zone, clean-up
by the units' own significance; the candidate of the largest matched-filter significance enters at each step instead of the first kernel whose dual reaches the wall. With the time stop it never
reaches the stop (the dual's entry times are the rule's input; the stored run went to the cap and was dropped, `exp_astro_m2.VARIANTS`), so the order is compared with the stop that both orders
share, the entrant's instantaneous significance (`stop='snr'`), at the frozen tuned configuration (sigma_b 3, r_ex 1.2, kappa 2) and the isotropic kernel (`exp_astro_m2.run_path`; the repeat patches
are run with the isotropic kernel as `iso|ours` of `repeat_raw.npz`). The greedy arm is `exp_astro_m2.VARIANTS['greedy_snr']`; it is run on the 16 test patches: the 4 M2 patches (`exp_astro_m2.TESTS`), the
4 M53 and 2 M92 patches (`exp_astro_m2.SECOND`) and the 6 repeat patches (`exp_astro_repeat.PATCHES`). HARNESS TEST: on the benchmark patch (310, 630) the run must reproduce the stored row
`V|x310_y630|greedy_snr` of `m2_raw.npz` (every readout, every stored array, and the entry count) bit for bit; the report prints the comparison (the same for the stored `snr` row and the
dual-order arm below). Added by this script, labelled so: the dual order with the SAME instantaneous-significance stop (`VARIANTS['snr']`, the stop-matched control of the order, about 20 s per patch) on
the 16 patches. Per group (M2 4, M53 4, M92 2, repeat 6, and all 16) the report pools the 20 x 20 px blocks of the patches and uses the stratified paired block bootstrap of the frozen scripts
(`exp_astro_m2.boot_idx_strat`, 2000 resamples, seed 0, the same resamples for every method of a table) to give completeness, purity and F with 95 % intervals and the paired differences
(ISS minus greedy, the question "does the dual order matter on frames"): the dual-order ISS tuned (the stored rows `T|`, `T2|`, `R|...|iso|ours|k2`) against the greedy arm; the dual order with the same
stop against the greedy arm; the greedy arm against crowdsource with the ISS front end (round 2, `baselines_hp2_raw.npz`, stored) and the ISS against it for reference; the guarded re-run of the tuned ISS
(`cyclefix_raw.npz`) in place of the stored one (added). Also descriptive: the greedy arm at the other kappa of the instantaneous stop, and the best kappa of the greedy arm CHOSEN ON THE TEST PATCHES (this favours
the greedy arm and is a bound, not a tuned result). Verdict on F: win if the 95 % interval of the difference is above 0, loss if below, tie otherwise. Nothing is tuned and no number of the paper changes.

    python experiments/exp_astro_controls.py --list                          the planned runs and their state
    python experiments/exp_astro_controls.py --run [--part c1|c2|all] [--first N]   run (one process at a time: a lock refuses a second one; write-once raw; a run whose result exists is not run again;
                                                                             c2 first, the harness test patch first)
    python experiments/exp_astro_controls.py --from-raw                      regenerate results/astro/controls.md from the raw and the stored raws (read only), byte for byte

Raw `controls_raw.npz` (`CONTROLS_OUT`, default results/astro): `E1|x<x0>_y<y0>|hp<s>|rex<r>|...` (control 1) and `E2|<set>|x<x0>_y<y0>|<arm>|...` (control 2, arm greedy or dual_snr), each with the scores of
every readout in the frozen scripts' layout, `secs` [path seconds, entries], `bic`, `trace|counts` (`exp_astro_cyclefix.COUNT_FIELDS`), `trace|barred`, `wall`; `meta|*`: the hashes of the frozen files, the sha256 of
this docstring (the declaration, checked at every run), versions. Inputs read only: `m2_raw.npz`, `repeat_raw.npz`, `cyclefix_raw.npz`, `baselines_hp2_raw.npz`
(`CONTROLS_M2_RAW`, `CONTROLS_REPEAT_RAW`, `CONTROLS_CYCLEFIX_RAW`, `CONTROLS_HP2_RAW`).
"""
import argparse
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
import exp_astro_cyclefix as CF           # noqa: E402
import exp_astro_m2 as X                  # noqa: E402
import exp_astro_repeat as RE             # noqa: E402

ASTRO = os.path.join(REPO, 'results', 'astro')
OUT = os.environ.get('CONTROLS_OUT', ASTRO)
RAW, MD = os.path.join(OUT, 'controls_raw.npz'), os.path.join(OUT, 'controls.md')
M2RAW = os.environ.get('CONTROLS_M2_RAW', os.path.join(ASTRO, 'm2_raw.npz'))
REPRAW = os.environ.get('CONTROLS_REPEAT_RAW', os.path.join(ASTRO, 'repeat_raw.npz'))
CFRAW = os.environ.get('CONTROLS_CYCLEFIX_RAW', os.path.join(ASTRO, 'cyclefix_raw.npz'))
HP2RAW = os.environ.get('CONTROLS_HP2_RAW', os.path.join(ASTRO, 'baselines_hp2_raw.npz'))
POSTHOC = 'POST HOC'

# control 1 (the declared set; the docstring above is the declaration)
EXT_KAPPAS = (1.25, 1.5, 1.75)
C1_KAPPAS = EXT_KAPPAS + tuple(X.KAPPAS)
C1_PATHS = ((3.0, 1.2), (3.0, 1.6), (3.0, 2.0), (2.5, 1.2), (2.5, 1.6), (2.5, 2.0), (2.0, 1.2), (2.0, 1.6), (2.0, 2.0), (3.0, 0.8), (5.0, 1.2), (5.0, 1.6), (8.0, 1.2))
FINDING = 0.01                                                    # an extended point beats the frozen one by more than this pooled development F: a finding
PACK_FIELDS = ('n_true', 'n_matched', 'n_det', 'n_det_matched', 'completeness', 'fp_fraction', 'offset', 'blk22', 'tot')

# control 2
ARMS = {'greedy': 'greedy_snr', 'dual_snr': 'snr'}              # arm -> exp_astro_m2.VARIANTS name
ARM_ORDER = ('greedy', 'dual_snr')
HARNESS = {'greedy': 'V|x310_y630|greedy_snr|', 'dual_snr': 'V|x310_y630|snr|'}      # the stored rows of the benchmark patch the harness must reproduce
HARNESS_PATCH = ('m2', (310, 630))
GROUP_TITLES = {'m2': 'M2 (NGC 7089), SDSS run 2583', 'm53': 'M53 (NGC 5024), SDSS run 5390', 'm92': 'M92 (NGC 6341), SDSS run 4682', 'rep': 'M2 repeat frames m2r1, m2r2, m2r3 (SDSS runs 4187, 7778, 3388)'}


def c2_patches():
    """The 16 test patches of the tuned ISS, as dicts(set, group, pxy), in the order M2, M53, M92, repeat: set = the data set name of the loaders and of `baselines_hp2_raw.npz`."""
    out = [dict(set='m2', group='m2', pxy=tuple(p)) for p in X.TESTS]
    for ds in ('m53', 'm92'):
        out += [dict(set=ds, group=ds, pxy=tuple(p)) for p in X.SECOND[ds]['patches']]
    out += [dict(set=name, group='rep', pxy=tuple(p)) for name, p in RE.PATCHES]
    return out


def iss_prefix(p):
    """The key prefix of the stored dual-order tuned ISS row (kappa 2) of a test patch, in `m2_raw.npz` / `repeat_raw.npz` and (with a leading 'C|') in `cyclefix_raw.npz`."""
    x, y = p['pxy']
    if p['group'] == 'm2':
        return 'T|x%d_y%d|ours|k2|' % (x, y)
    if p['group'] in ('m53', 'm92'):
        return 'T2|%s|x%d_y%d|ours|k2|' % (p['set'], x, y)
    return 'R|%s|x%d_y%d|iso|ours|k2|' % (p['set'], x, y)


def c2_pre(p, arm):
    return 'E2|%s|x%d_y%d|%s|' % (p['set'], p['pxy'][0], p['pxy'][1], arm)


def c1_pre(p, h, r):
    return 'E1|x%d_y%d|hp%g|rex%g|' % (p[0], p[1], h, r)


def utc():
    return time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())


def declaration_sha():
    return hashlib.sha256((__doc__ or '').encode()).hexdigest()


# ----------------------------------------------------------------------------------------------------------------------------------
# the runs
# ----------------------------------------------------------------------------------------------------------------------------------
def jobs(part='all'):
    """The planned runs, in the order they are run: control 2 first, the harness test patch first (the greedy arm, then the others), then control 1 (the frozen combination first)."""
    out = []
    if part in ('c2', 'all'):
        pats = c2_patches()
        pats.sort(key=lambda p: (p['set'], p['pxy']) != (HARNESS_PATCH[0], HARNESS_PATCH[1]))
        for arm in ARM_ORDER:
            out += [dict(kind='c2', arm=arm, **p) for p in pats]
    if part in ('c1', 'all'):
        for h, r in C1_PATHS:
            out += [dict(kind='c1', pxy=tuple(p), hp=h, rex=r) for p in X.DEV]
    for s in out:
        s['prefix'] = c1_pre(s['pxy'], s['hp'], s['rex']) if s['kind'] == 'c1' else c2_pre(s, s['arm'])
    return out


def load_for(spec):
    """(patch, truth of the patch) of a spec; control 1 refuses everything but the two development patches."""
    if spec['kind'] == 'c1':
        if tuple(spec['pxy']) not in [tuple(p) for p in X.DEV]:
            raise ValueError('control 1 runs the development patches only, not %s' % (spec['pxy'],))
        return X.load_patch(*spec['pxy']), A.truth_for_patch(X.truth_arrays(), spec['pxy'][0], spec['pxy'][1])
    x, y = spec['pxy']
    if spec['group'] == 'rep':
        return R.load_patch(x, y, spec['set']), A.truth_for_patch(R.truth_arrays(spec['set']), x, y)
    return X.load_patch(x, y, ds=spec['set']), A.truth_for_patch(X.truth_arrays(spec['set']), x, y)


def frozen_cfg():
    fr = RE.load_frozen()
    return fr['hp'], fr['rex'], fr['kap']


def run_spec(spec):
    """The arrays (keys relative to the prefix) of one run, computed by `exp_astro_m2.run_path` with the guard in place of `imaging.sliding_path`; the guard is installed here and removed afterwards."""
    hp, rex, kap = frozen_cfg()
    patch, t = load_for(spec)
    CF.install()
    CF.TRACE.clear()
    try:
        t0 = time.time()
        if spec['kind'] == 'c1':
            dets, secs, n, info = X.run_path(patch, spec['hp'], spec['rex'], kappas=C1_KAPPAS)
        else:
            ov = dict(X.VARIANTS[ARMS[spec['arm']]])
            dets, secs, n, info = X.run_path(patch, ov.pop('hp', hp), ov.pop('rex', rex), **ov)
        wall = time.time() - t0
        if len(CF.TRACE) != 1:
            raise RuntimeError('%d paths were run, expected one' % len(CF.TRACE))
        tr = CF.TRACE[0]
    finally:
        CF.uninstall()
        CF.TRACE.clear()
    arrays = CF.store_run({}, dets, secs, n, info, t, kap, det=spec['kind'] == 'c2', bic=True)
    arrays['trace|counts'], arrays['trace|barred'], arrays['trace|taus'], arrays['wall'] = tr['counts'], tr['barred'], tr['taus'], np.array([wall])
    return arrays


def meta_record():
    import scipy
    hs = CF.hashes()
    return {'meta|files': np.array(sorted(hs)), 'meta|sha256': np.array([hs[f] for f in sorted(hs)]), 'meta|declaration_sha256': np.array([declaration_sha()]),
            'meta|inputs': np.array(['m2_raw.npz', 'repeat_raw.npz', 'cyclefix_raw.npz', 'baselines_hp2_raw.npz']),
            'meta|inputs_sha256': np.array([CF.sha256_file(p) for p in (M2RAW, REPRAW, CFRAW, HP2RAW)]),
            'meta|first_run_utc': np.array([utc()]),
            'meta|versions': np.array(['python %d.%d.%d' % sys.version_info[:3], 'numpy ' + np.__version__, 'scipy ' + scipy.__version__])}


def preflight(out_path=None):
    """Refuses to run unless the frozen modules match their freeze records, the hashed files are those of the first run, and this docstring (the declaration) is the one of the first run. Returns
    (the current raw, what is to be recorded with the next save)."""
    out_path = out_path or RAW
    X.check_freeze(X.load_raw(M2RAW))
    RE.check_freeze()
    cur = CF.load_raw(out_path) if os.path.exists(out_path) else {}
    now = meta_record()
    if 'meta|sha256' not in cur:
        return cur, now
    rec = dict(zip((str(f) for f in cur['meta|files']), (str(h) for h in cur['meta|sha256'])))
    now_h = dict(zip((str(f) for f in now['meta|files']), (str(h) for h in now['meta|sha256'])))
    diff = sorted(f for f, h in rec.items() if now_h.get(f) != h)
    if diff:
        raise RuntimeError('frozen files differ from what the first run recorded: %s' % ', '.join(diff))
    if str(cur['meta|declaration_sha256'][0]) != declaration_sha():
        raise RuntimeError('the docstring (the declaration) differs from the one recorded by the first run')
    return cur, {}


def run(part='all', first=None, out_path=None):
    """Each planned run once (a run whose result exists is not run again), saved as it arrives; at most one process of this script at a time (a lock)."""
    out_path = out_path or RAW
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    with open(out_path + '.runlock', 'a') as lk:
        try:
            fcntl.flock(lk, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            raise RuntimeError('another process of this script is running (lock %s.runlock)' % out_path)
        cur, pending = preflight(out_path)
        planned = jobs(part)
        todo = [s for s in planned if s['prefix'] + 'wall' not in cur]
        n_exist = len(planned) - len(todo)
        if first is not None:
            todo = todo[:first]
        print('%s %d runs to do (%s), %d exist' % (utc(), len(todo), part, n_exist), flush=True)
        for s in todo:
            print('%s start %s' % (utc(), s['prefix']), flush=True)
            arrays = run_spec(s)
            c = arrays['trace|counts']
            print('%s done  %s: %d entries, %d guard activations, cap %s, %.0f s' % (utc(), s['prefix'], c[0], c[6], 'HIT' if c[0] >= X.MAX_ENTRIES else 'no', arrays['wall'][0]), flush=True)
            new = {s['prefix'] + k: v for k, v in arrays.items()}
            new.update(pending)
            pending = {}
            CF.save_part(out_path, new)
    return len(todo)


def list_jobs(out_path=None):
    out_path = out_path or RAW
    cur = CF.load_raw(out_path) if os.path.exists(out_path) else {}
    for s in jobs('all'):
        print('%-5s %-60s %s' % (s['kind'], s['prefix'], 'done' if s['prefix'] + 'wall' in cur else 'todo'))


# ----------------------------------------------------------------------------------------------------------------------------------
# control 1: the surface
# ----------------------------------------------------------------------------------------------------------------------------------
def _row(z, pre, k, cap, act):
    return dict(tot=np.asarray(z[pre + 'k%g|tot' % k], float), blk=np.asarray(z[pre + 'k%g|blk22' % k], float), cap=cap, act=act)


def old_surface(ctx, mode):
    """{(sigma_b, r_ex, kappa): {patch: row}} of the old development grid (kappa >= 2), row = dict(tot, blk, cap, act): mode 'guarded': the guarded paths of `cyclefix_raw` (`G|` for the six capped
    paths, `C|` for the others); mode 'stored': `m2_raw`, as run (unguarded). A combination needs both patches."""
    out = {}
    for h in X.HPS:
        for r in X.REXS:
            per = {}
            for p in X.DEV:
                rid = 'D|x%d_y%d|hp%g|rex%g|' % (p[0], p[1], h, r)
                if mode == 'stored':
                    if rid + 'bic' not in ctx['m2']:
                        break
                    z, pre, cap, act = ctx['m2'], rid, int(ctx['m2'][rid + 'bic'][2]), 0
                else:
                    ns = 'G' if 'G|' + rid + 'wall' in ctx['cf'] else 'C'
                    if ns + '|' + rid + 'wall' not in ctx['cf']:
                        break
                    z, pre = ctx['cf'], ns + '|' + rid
                    c = ctx['cf'][pre + 'trace|counts']
                    cap, act = int(c[0] >= X.MAX_ENTRIES), int(c[6])
                per[p] = (z, pre, cap, act)
            if len(per) == len(X.DEV):
                for k in X.KAPPAS:
                    out[(h, r, k)] = {p: _row(z, pre, k, cap, act) for p, (z, pre, cap, act) in per.items()}
    return out


def ext_surface(raw):
    """{(sigma_b, r_ex, kappa): {patch: row}} of the control-1 paths of the raw (both patches present)."""
    out = {}
    for h, r in C1_PATHS:
        if not all(c1_pre(p, h, r) + 'wall' in raw for p in X.DEV):
            continue
        for k in C1_KAPPAS:
            row = {}
            for p in X.DEV:
                pre = c1_pre(p, h, r)
                c = raw[pre + 'trace|counts']
                row[p] = _row(raw, pre, k, int(c[0] >= X.MAX_ENTRIES), int(c[6]))
            out[(h, r, k)] = row
    return out


def merge_surface(old, ext):
    """The old rows, plus the control-1 rows at the (sigma_b, r_ex, kappa) the old grid does not hold; each entry is a dict with tot (the two patches' counts added), per (by patch), blk, cap (patches capped), act."""
    out = {}
    for src in (old, ext):
        for key, per in src.items():
            if key not in out:
                out[key] = dict(tot=sum(per[p]['tot'] for p in X.DEV), per=[per[p]['tot'] for p in X.DEV], blk=[per[p]['blk'] for p in X.DEV], cap=sum(per[p]['cap'] for p in X.DEV),
                                act=sum(per[p]['act'] for p in X.DEV), src='control 1' if src is ext else 'old grid')
    return out


def choose_surface(surf, uncapped_only=False):
    """(sigma_b, r_ex, kappa, F) of the largest pooled F at mag < 22, by the rule of `exp_astro_m2.choose`: the keys in ascending order (sigma_b, r_ex, kappa), a row replaces the best if (round(F, 6), kappa)
    is larger; `uncapped_only` leaves out the rows of paths that hit the entry cap."""
    best = None
    for key in sorted(surf):
        if uncapped_only and surf[key]['cap']:
            continue
        f = X.f_measure(surf[key]['tot'])[0]
        if best is None or (round(f, 6), key[2]) > (round(best[0], 6), best[3]):
            best = (f,) + key
    return None if best is None else (best[1], best[2], best[3], best[0])


def F_of(row):
    return X.f_measure(row['tot'])[0]


def Fp(row, i):
    return X.f_measure(row['per'][i])[0]


def pivot_table(surf, combos, frozen, best, ks=C1_KAPPAS):
    rows = []
    for h, r in combos:
        cells = ['%g' % h, '%g' % r]
        for k in ks:
            v = surf.get((h, r, k))
            if v is None:
                cells.append('-')
                continue
            s = '%.4f' % F_of(v) + (' (cap)' if v['cap'] else '')
            if (h, r, k) == frozen:
                s = '**' + s + '**'
            elif (h, r, k) == best:
                s += ' (max)'
            cells.append(s)
        rows.append(cells)
    return X.table(['sigma_b px', 'r_ex px'] + ['kappa %g' % k for k in ks], rows)


def combos_of(surf):
    return sorted({(h, r) for h, r, k in surf})


def c1_checks(raw, ctx):
    """[(sigma_b, r_ex, readouts compared, identical)] for the old combinations re-run by control 1: the readouts of kappa >= 2 against the guarded rows of cyclefix_raw."""
    out = []
    for h, r in C1_PATHS:
        if (h, r) not in [(a, b) for a in X.HPS for b in X.REXS] or not all(c1_pre(p, h, r) + 'wall' in raw for p in X.DEV):
            continue
        n = same = 0
        for p in X.DEV:
            rid = 'D|x%d_y%d|hp%g|rex%g|' % (p[0], p[1], h, r)
            ns = 'G' if 'G|' + rid + 'wall' in ctx['cf'] else 'C'
            for k in X.KAPPAS:
                for f in PACK_FIELDS:
                    a, b = raw.get(c1_pre(p, h, r) + 'k%g|%s' % (k, f)), ctx['cf'].get(ns + '|' + rid + 'k%g|%s' % (k, f))
                    if a is None or b is None:
                        continue
                    n += 1
                    same += int(CF.same_array(a, b))
        out.append((h, r, n, same))
    return out


def c1_section(raw, ctx, summ=None):
    L = ['## Control 1: the development optimum against the edges of the grid (%s; development patches %s only)\n' % (POSTHOC, ', '.join('(%d, %d)' % p for p in X.DEV)),
         'The grid of `exp_astro_m2` (sigma_b {3, 5, 8}, r_ex {0.8, 1.2}, kappa {2, ..., 8}) is extended beyond every edge the frozen point (sigma_b 3, r_ex 1.2, kappa 2) touches: sigma_b {2, 2.5}, r_ex '
         '{1.6, 2.0}, kappa {1.25, 1.5, 1.75}, on the two development patches only (no test patch is loaded or scored by this control); guarded paths; scored by `exp_astro_m2.choose`\'s rule: pooled F at mag < 22, '
         'the two patches\' counts added, ties on round(F, 6) to the larger kappa. The set is declared in the script\'s docstring (13 paths per patch, each read at 11 kappas).\n']
    present = [(h, r) for h, r in C1_PATHS if all(c1_pre(p, h, r) + 'wall' in raw for p in X.DEV)]
    if not present:
        return L + ['No control-1 path has been run.\n']
    old_g, old_s, ext = old_surface(ctx, 'guarded'), old_surface(ctx, 'stored'), ext_surface(raw)
    sg, ss = merge_surface(old_g, ext), merge_surface(old_s, ext)
    frozen = (3.0, 1.2, 2.0)
    done = len(present)
    L.append('Paths run: %d of %d per patch (%s).\n' % (done, len(C1_PATHS), ', '.join('(%g, %g)' % c for c in present)))
    # run record
    rows = []
    for h, r in present:
        for p in X.DEV:
            pre = c1_pre(p, h, r)
            c, secs = raw[pre + 'trace|counts'], raw[pre + 'secs']
            rows.append(['(%d, %d)' % p, '%g' % h, '%g' % r, '%d' % c[0], 'HIT' if c[0] >= X.MAX_ENTRIES else 'no', '%d' % c[6], '%d' % c[1], '%d' % c[3], '%.0f' % secs[0], '%.0f' % raw[pre + 'wall'][0]])
    L += ['### Runs: entries, cap flag, guard activations\n', X.table(['patch', 'sigma_b', 'r_ex', 'entries (path to the kappa = 1.25 time)', 'entry cap (%d)' % X.MAX_ENTRIES, 'guard activations', 'exit events',
                                                                      'most entries of one lattice point', 'path seconds', 'wall seconds'], rows), '']
    chk = c1_checks(raw, ctx)
    if chk:
        L.append('Old combinations re-run to the kappa = 1.25 time, readouts of kappa >= 2 (every stored array of `exp_astro_m2._pack`) against the guarded development rows of `cyclefix_raw.npz` (the path to a larger time is a prefix-extension, '
                 'so these must be identical): %s.\n' % '; '.join('(%g, %g): %d of %d identical' % c for c in chk))
    # argmax
    L += ['### The argmax: does it move?\n']
    for name, surf in (('guarded old rows (all twelve old paths guarded) and control 1', sg), ('stored old rows (as run, unguarded) and control 1', ss)):
        if not surf:
            continue
        h, r, k, f = choose_surface(surf)
        ref = surf.get(frozen)
        fr = F_of(ref) if ref else float('nan')
        L.append('* Surface of %s: argmax sigma_b %g, r_ex %g, kappa %g, pooled F %.4f; frozen point (3, 1.2, 2) F %.4f; difference %+.4f; **%s**.' % (
            name, h, r, k, f, fr, f - fr, 'the argmax is the frozen point' if (h, r, k) == frozen else 'the argmax moves'))
        hu, ru, ku, fu = choose_surface(surf, uncapped_only=True)
        L.append('  Over the rows of paths that did not hit the cap: sigma_b %g, r_ex %g, kappa %g, pooled F %.4f.' % (hu, ru, ku, fu))
    L.append('')
    h, r, k, f = choose_surface(sg)
    best = (h, r, k)
    fg = F_of(sg[frozen]) if frozen in sg else float('nan')
    fs = F_of(ss[frozen]) if frozen in ss else float('nan')
    over = sorted(((F_of(v), key) for key, v in sg.items() if v['src'] == 'control 1' and F_of(v) - max(fg, fs) > FINDING), reverse=True)
    if summ is not None:
        summ['c1'] = dict(best=best, f=f, fg=fg, fs=fs, nover=len(over), n=len(present))
    L += ['### Findings (an extended point that beats the frozen point by more than %.2f pooled development F, against the guarded value %.4f and the stored value %.4f alike)\n' % (FINDING, fg, fs)]
    if over:
        L.append('**FINDING (%s, development patches only, not run on any test patch):** %d extended point(s) beat the frozen point by more than %.2f: %s.\n' % (
            POSTHOC, len(over), FINDING, '; '.join('(sigma_b %g, r_ex %g, kappa %g) F %.4f (%+.4f)' % (key + (v, v - fg)) for v, key in over[:10])))
    else:
        L.append('No extended point beats the frozen point by more than %.2f pooled development F. The largest excess over the guarded frozen value is %+.4f (at sigma_b %g, r_ex %g, kappa %g); the number of '
                 'extended points within 0.005 of the frozen value is %d.\n' % (FINDING, f - fg, h, r, k, sum(1 for key, v in sg.items() if key != frozen and abs(F_of(v) - fg) <= 0.005 and v['src'] == 'control 1')))
    # edges
    ext_ax = dict(sigma_b=sorted({kk[0] for kk in sg}), r_ex=sorted({kk[1] for kk in sg}), kappa=sorted({kk[2] for kk in sg}))
    edges = []
    for nm, v, g in (('sigma_b', best[0], ext_ax['sigma_b']), ('r_ex', best[1], ext_ax['r_ex']), ('kappa', best[2], ext_ax['kappa'])):
        if v in (g[0], g[-1]):
            edges.append('%s %g is at the %s edge of the extended grid %s' % (nm, v, 'lower' if v == g[0] else 'upper', tuple(g)))
    L.append('Position of the extended argmax on the extended grid: %s.\n' % ('; '.join(edges) if edges else 'interior on every axis'))
    # bootstrap of the best against the frozen
    if best != frozen:
        blk_b, blk_f = np.concatenate(sg[best]['blk']), np.concatenate(sg[frozen]['blk'])
        idx = X.boot_idx_strat([len(b) for b in sg[frozen]['blk']])
        dc, dp, df = X.boot_diff(blk_b, blk_f, idx)
        L.append('Paired block bootstrap on the development blocks (stratified within the two patches, %d resamples, seed 0) of the extended argmax minus the frozen point: dF %s (selected as the maximum of a larger '
                 'grid on these same blocks, so the interval is optimistic and is not a held-out statement).\n' % (X.BOOT, X.fmt_iv(df[0], df[1], 4)))
    # profiles
    L += ['### The F surface around the frozen point (guarded rows; per-patch F in brackets: (110, 630) / (210, 630))\n']
    prof, axis_top = [], {}
    for ax, name, keys in (('sigma_b', 'sigma_b', [(h_, frozen[1], frozen[2]) for h_ in ext_ax['sigma_b']]), ('r_ex', 'r_ex', [(frozen[0], r_, frozen[2]) for r_ in ext_ax['r_ex']]),
                           ('kappa', 'kappa', [(frozen[0], frozen[1], k_) for k_ in ext_ax['kappa']])):
        axis_top[name] = max([F_of(sg[key]) for key in keys if key in sg] or [float('nan')]) <= fg + 1e-12
        for key in keys:
            if key not in sg:
                continue
            v = sg[key]
            prof.append([name, '%g' % key[{'sigma_b': 0, 'r_ex': 1, 'kappa': 2}[ax]], '%.4f' % F_of(v), '%+.4f' % (F_of(v) - fg), '%.4f / %.4f' % (Fp(v, 0), Fp(v, 1)),
                         '%.3f / %.3f' % X.cpf(v['tot'])[:2], '%d' % v['tot'][2], 'cap' if v['cap'] else '-'])
    L += [X.table(['axis through the frozen point', 'value', 'pooled F', 'dF to frozen', 'F by patch', 'completeness / purity', 'detections < 22', 'cap'], prof), '']
    above = sorted((F_of(v) - fg, key) for key, v in sg.items() if F_of(v) > fg)
    L.append('Along each axis through the frozen point (the other two values fixed) the frozen value is the largest: sigma_b %s, r_ex %s, kappa %s. Points of the extended surface with a larger pooled F than the frozen point (guarded value %.4f): '
             '%d, of which %d beyond the old edge; the largest excess is %s.\n' % ('yes' if axis_top['sigma_b'] else 'no', 'yes' if axis_top['r_ex'] else 'no', 'yes' if axis_top['kappa'] else 'no', fg, len(above),
                                                                                    sum(1 for _, key in above if sg[key]['src'] == 'control 1'), ('%+.4f at sigma_b %g, r_ex %g, kappa %g' % ((above[-1][0],) + above[-1][1])) if above else 'none'))
    L += ['### The full extended surface: pooled F at mag < 22\n', 'Rows: (sigma_b, r_ex); columns: kappa. Guarded old rows (kappa >= 2) where the control did not run the combination. Bold: the frozen point; "(max)": the argmax of the table; '
          '"(cap)": the path hit the entry cap (censored). "-": not run.\n', pivot_table(sg, combos_of(sg), frozen, best), '']
    L += ['### The ten best rows (the order of `choose`: round(F, 6), then the larger kappa)\n']
    rk = sorted(sg, key=lambda key: (round(F_of(sg[key]), 6), key[2]), reverse=True)[:10]
    L += [X.table(['rank', 'sigma_b', 'r_ex', 'kappa', 'pooled F', 'F by patch', 'source', 'cap'], [[str(i + 1), '%g' % kk[0], '%g' % kk[1], '%g' % kk[2], '%.4f' % F_of(sg[kk]),
                                                                                                      '%.4f / %.4f' % (Fp(sg[kk], 0), Fp(sg[kk], 1)), sg[kk]['src'], 'cap' if sg[kk]['cap'] else '-'] for i, kk in enumerate(rk)]), '']
    L += ['The table of the stored (as run) old rows in place of the guarded ones:\n', pivot_table(ss, combos_of(ss), frozen, choose_surface(ss)[:3] if ss else None), '']
    L.append('Reading: the maximum over a larger grid on the same two patches is optimistically biased; the control shows whether the frozen point sits on a slope that continues beyond the old edge (a better point beyond '
             'it) or on a plateau/ridge, not a held-out value; no test patch was scored for it.\n')
    return L


# ----------------------------------------------------------------------------------------------------------------------------------
# control 2
# ----------------------------------------------------------------------------------------------------------------------------------
def harness_checks(raw, ctx):
    """[(arm, arrays compared, identical, entries equal, entries stored, entries now)] of the arms on the benchmark patch against the stored `V|...` rows of m2_raw."""
    out = []
    p = dict(set='m2', group='m2', pxy=HARNESS_PATCH[1])
    for arm in ARM_ORDER:
        pre = c2_pre(p, arm)
        st = HARNESS[arm]
        if pre + 'wall' not in raw or st + 'secs' not in ctx['m2']:
            continue
        n = same = 0
        for k in X.KAPPAS:
            for f in PACK_FIELDS:
                a, b = raw.get(pre + 'k%g|%s' % (k, f)), ctx['m2'].get(st + 'k%g|%s' % (k, f))
                if a is None or b is None:
                    continue
                n += 1
                same += int(CF.same_array(a, b))
        out.append((arm, n, same, int(raw[pre + 'secs'][1]) == int(ctx['m2'][st + 'secs'][1]), int(ctx['m2'][st + 'secs'][1]), int(raw[pre + 'secs'][1])))
    return out


def method_sources(ctx, raw):
    """[(label, getter(patch) -> prefix-source pair)]: each method is a function patch -> (blk22 array, tot array) or None if absent."""
    def stored(z, ns=''):
        def f(p):
            pre = ns + iss_prefix(p)
            return (z[pre + 'blk22'], z[pre + 'tot']) if pre + 'blk22' in z else None
        return f

    def mine(arm, k='k2'):
        def f(p):
            pre = c2_pre(p, arm)
            return (raw[pre + k + '|blk22'], raw[pre + k + '|tot']) if pre + k + '|blk22' in raw else None
        return f

    def hp2(row):
        def f(p):
            pre = 'T|%s|x%d_y%d|%s|' % (p['set'], p['pxy'][0], p['pxy'][1], row)
            return (ctx['hp2'][pre + 'blk22'], ctx['hp2'][pre + 'tot']) if pre + 'blk22' in ctx['hp2'] else None
        return f
    zz = dict(m2=ctx['m2'], rep=ctx['rep'])

    def iss(p):
        return stored(zz['rep'] if p['group'] == 'rep' else zz['m2'])(p)

    def iss_cf(p):
        return stored(ctx['cf'], 'C|')(p)
    return dict(iss=('ISS, dual order, tuned (stored rows; kappa 2)', iss), iss_cf=('ISS, dual order, tuned, guarded re-run (cyclefix; kappa 2)', iss_cf),
                greedy=('greedy order, instantaneous-significance stop, kappa 2 (this control)', mine('greedy')), dual_snr=('dual order, the same instantaneous-significance stop, kappa 2 (this control, added)', mine('dual_snr')),
                cs2=('crowdsource with the ISS front end (round 2, stored)', hp2('cs2')), phu2=('photutils with the ISS front end (round 2, stored)', hp2('phu2')))


def pool(meths, key, plist):
    """(blocks of the patches concatenated, sizes, totals summed) of a method over the patches, or None if any patch lacks it."""
    got = [meths[key][1](p) for p in plist]
    if not got or any(g is None for g in got):
        return None
    return np.concatenate([g[0] for g in got]), [len(g[0]) for g in got], sum(np.asarray(g[1], float) for g in got)


def c2_section(raw, ctx, summ=None):
    L = ['## Control 2: is the sliding ISS just a greedy fitter? (%s; the dual order against the greedy order on the 16 test patches)\n' % POSTHOC,
         'The same code, the same PSF, high-pass, weights, exclusion zone, clean-up and kappa; only the order of entry differs: the accumulated dual (the sliding inverse scale space) or the instantaneous '
         'matched-filter significance (sliding matching pursuit, `order=\'greedy\'`). The stop is the one both orders share, the entrant\'s instantaneous significance (the time stop of the ISS needs the dual\'s entry '
         'times, the greedy path\'s never reach it); the dual order with that same stop is added as the stop-matched control. Frozen tuned configuration (sigma_b 3, r_ex 1.2, kappa 2), isotropic kernel, guarded '
         'function (for `order=\'greedy\'` the guard never acts: tau is the step count, no two exits share an entry time; the activations below are those the guard reports, and for the dual arm they are real).\n']
    pats = c2_patches()
    have = [p for p in pats if c2_pre(p, 'greedy') + 'wall' in raw]
    if not have:
        return L + ['No control-2 path has been run.\n']
    meths = method_sources(ctx, raw)
    # harness
    L += ['### Harness test: the benchmark patch (310, 630) reproduces the stored design-variant rows\n']
    hc = harness_checks(raw, ctx)
    if hc:
        L.append(X.table(['arm', 'stored row', 'arrays compared (kappa 2 to 8, 9 fields each)', 'identical', 'entry count stored / now', 'bit for bit'],
                         [[a, HARNESS[a][:-1].replace('|', '/'), str(n), str(s), '%d / %d' % (e0, e1), 'YES' if (n > 0 and s == n and eq) else 'NO'] for a, n, s, eq, e0, e1 in hc]))
        L.append('')
    else:
        L.append('The benchmark patch has not been run yet.\n')
    # runs
    rows = []
    for p in have:
        for arm in ARM_ORDER:
            pre = c2_pre(p, arm)
            if pre + 'wall' not in raw:
                continue
            c = raw[pre + 'trace|counts']
            rows.append(['%s (%d, %d)' % (p['set'], p['pxy'][0], p['pxy'][1]), arm, '%d' % c[0], 'HIT' if c[0] >= X.MAX_ENTRIES else 'no', '%d' % c[6], '%d' % c[1], '%.0f' % raw[pre + 'secs'][0], '%.0f' % raw[pre + 'wall'][0]])
    L += ['### Runs: entries, cap flag, guard activations\n', X.table(['patch', 'arm', 'entries', 'entry cap (%d)' % X.MAX_ENTRIES, 'guard activations', 'exit events', 'path seconds', 'wall seconds'], rows), '']
    allp = []
    for g in ('m2', 'm53', 'm92', 'rep'):
        pl = [p for p in pats if p['group'] == g]
        allp.append((g, GROUP_TITLES[g], pl))
    allp.append(('all', 'all 16 test patches', pats))
    for _g, title, pl in allp:
        if any(c2_pre(p, 'greedy') + 'wall' not in raw for p in pl):
            L += ['### %s: not complete (%d of %d patches run)\n' % (title, sum(1 for p in pl if c2_pre(p, 'greedy') + 'wall' in raw), len(pl))]
            continue
        order = [m for m in ('iss', 'iss_cf', 'greedy', 'dual_snr', 'cs2', 'phu2') if pool(meths, m, pl) is not None]
        pooled = {m: pool(meths, m, pl) for m in order}
        sizes = pooled['iss'][1]
        idx = X.boot_idx_strat(sizes)
        L += ['### %s (%d patches, %d blocks): pooled, stratified paired block bootstrap (%d resamples, seed 0)\n' % (title, len(pl), sum(sizes), X.BOOT)]
        rows = []
        for m in order:
            blk, _, tot = pooled[m]
            (c, cci), (pu, pci), (ff, fci) = X.boot_cpf(blk, idx)
            rows.append([meths[m][0], '%d / %d' % (tot[1], tot[0]), X.fmt_iv(c, cci), '%d / %d' % (tot[3], tot[2]), X.fmt_iv(pu, pci), X.fmt_iv(ff, fci)])
        L += [X.table(['catalogue', 'matched / truth < 22', 'completeness', 'matched / detected < 22', 'purity', 'F'], rows), '']
        if summ is not None:
            summ.setdefault('f2', {})[_g] = {m: X.cpf(pooled[m][0].sum(0))[2] for m in pooled}
        rows = []
        pairs = [('iss', 'greedy', 'THE QUESTION: dual-order ISS tuned minus greedy order'), ('dual_snr', 'greedy', 'dual order minus greedy order, the same stop (order only)'), ('iss', 'dual_snr', 'tuned ISS (time stop) minus dual order with the instantaneous stop (the stop only)'),
                 ('iss_cf', 'greedy', 'guarded tuned ISS minus greedy order'), ('greedy', 'cs2', 'greedy order minus crowdsource with the front end'), ('iss', 'cs2', 'ISS tuned minus crowdsource with the front end (reference)')]
        for a, b, name in pairs:
            if a in pooled and b in pooled:
                dc, dp, df = X.boot_diff(pooled[a][0], pooled[b][0], idx)
                if summ is not None:
                    summ.setdefault('c2', {})[(_g, a, b)] = (df[0], df[1][0], df[1][1], X.verdict(*df))
                rows.append([name, X.fmt_iv(dc[0], dc[1]), X.fmt_iv(dp[0], dp[1]), X.fmt_iv(df[0], df[1]), X.verdict(*df)])
        L += ['Paired differences (a minus b, the same resampled blocks), 95 % intervals:\n', X.table(['comparison', 'd completeness', 'd purity', 'd F', 'verdict on F'], rows), '']
        # greedy at the other kappa
        kr = []
        gk = {}
        for k in X.KAPPAS:
            got = [(raw.get(c2_pre(p, 'greedy') + 'k%g|blk22' % k), raw.get(c2_pre(p, 'greedy') + 'k%g|tot' % k)) for p in pl]
            dk = [(raw.get(c2_pre(p, 'dual_snr') + 'k%g|blk22' % k), raw.get(c2_pre(p, 'dual_snr') + 'k%g|tot' % k)) for p in pl]
            if any(a is None for a, _ in got):
                continue
            bg = np.concatenate([a for a, _ in got])
            gk[k] = bg
            (c, cci), (pu, pci), (ff, fci) = X.boot_cpf(bg, idx)
            row = ['%g' % k, X.fmt_iv(c, cci), X.fmt_iv(pu, pci), X.fmt_iv(ff, fci)]
            if all(a is not None for a, _ in dk):
                bd = np.concatenate([a for a, _ in dk])
                row.append(X.fmt_iv(*X.boot_cpf(bd, idx)[2]))
            else:
                row.append('-')
            kr.append(row)
        if kr:
            L += ['Descriptive: the readouts at the other kappa of the instantaneous stop (no choice rests on these rows):\n', X.table(['kappa', 'greedy: completeness', 'greedy: purity', 'greedy: F', 'dual order, same stop: F'], kr), '']
            bk = max(gk, key=lambda k: (round(X.cpf(gk[k].sum(0))[2], 6), k))
            dc, dp, df = X.boot_diff(pooled['iss'][0], gk[bk], idx)
            if summ is not None:
                summ.setdefault('bound', {})[_g] = (bk, df[0], df[1][0], df[1][1], X.verdict(*df))
            L += ['Bound in favour of the greedy order: the greedy kappa CHOSEN ON THESE TEST PATCHES is %g (largest pooled F of the greedy arm); the tuned ISS (kappa 2, chosen on the development patches) minus the greedy order at that kappa: dF %s (%s).\n' % (
                bk, X.fmt_iv(df[0], df[1]), X.verdict(*df))]
    # per patch
    rows = []
    for p in have:
        res = {m: meths[m][1](p) for m in ('iss', 'greedy', 'dual_snr', 'cs2')}
        if any(v is None for v in res.values()):
            continue
        idx = X.boot_idx(len(res['iss'][0]))
        dc, dp, df = X.boot_diff(res['iss'][0], res['greedy'][0], idx)
        fs = {m: X.f_measure(np.asarray(v[1], float))[0] for m, v in res.items()}
        rows.append(['%s (%d, %d)' % (p['set'], p['pxy'][0], p['pxy'][1]), '%d' % res['iss'][1][0], '%.3f' % fs['iss'], '%.3f' % fs['greedy'], '%.3f' % fs['dual_snr'], '%.3f' % fs['cs2'], X.fmt_iv(df[0], df[1]), X.verdict(*df)])
    if rows:
        L += ['### By patch: F at mag < 22 (point values) and the paired difference of the dual-order ISS and the greedy order with its block-bootstrap interval (25 blocks, seed 0)\n',
              X.table(['patch', 'truth < 22', 'F ISS dual, tuned', 'F greedy', 'F dual, same stop', 'F crowdsource front end', 'dF ISS minus greedy', 'verdict'], rows), '']
    return L


# ----------------------------------------------------------------------------------------------------------------------------------
# report
# ----------------------------------------------------------------------------------------------------------------------------------
def load_ctx():
    return dict(m2=CF.load_raw(M2RAW), rep=CF.load_raw(REPRAW), cf=CF.load_raw(CFRAW), hp2=CF.load_raw(HP2RAW))


def summary_lines(summ):
    """The headline of the report, from the numbers the two sections computed: the argmax of control 1 and the ISS minus greedy rows of control 2."""
    L = ['## Summary (%s)\n' % POSTHOC] if summ else []
    if 'c2' in summ:
        rows = []
        for g in ('m2', 'm53', 'm92', 'rep', 'all'):
            if (g, 'iss', 'greedy') not in summ['c2']:
                continue
            f, d = summ['f2'][g], summ['c2'][(g, 'iss', 'greedy')]
            ds, bd = summ['c2'].get((g, 'dual_snr', 'greedy')), summ.get('bound', {}).get(g)
            rows.append([GROUP_TITLES.get(g, 'all 16 test patches').split(', SDSS')[0].split(' (SDSS')[0], '%.3f' % f['iss'], '%.3f' % f['greedy'], '%+.3f [%+.3f, %+.3f]' % d[:3], d[3],
                         ('%+.3f [%+.3f, %+.3f] %s' % (ds[0], ds[1], ds[2], ds[3])) if ds else '-', ('%+.3f [%+.3f, %+.3f] %s (kappa %g)' % (bd[1], bd[2], bd[3], bd[4], bd[0])) if bd else '-'])
        L += ['Control 2, pooled F at mag < 22 of the dual-order ISS tuned (stored rows) and of the greedy order (this control), and their paired difference (ISS minus greedy) with the 95 % interval '
              'of the stratified block bootstrap, the same difference with the stop matched (dual order with the instantaneous stop minus greedy), and the bound with the greedy kappa chosen on the test patches:\n',
              X.table(['group', 'F ISS tuned', 'F greedy', 'ISS minus greedy', 'verdict', 'dual order, same stop, minus greedy', 'ISS minus greedy at the greedy arm\'s best kappa'], rows), '']
    if 'c1' in summ:
        c = summ['c1']
        L += ['Control 1 (development patches only): over the extended grid (%d paths per patch run) the argmax is sigma_b %g, r_ex %g, kappa %g with pooled F %.4f; the frozen point (3, 1.2, 2) has %.4f (guarded) and %.4f (stored): '
              '**%s**; extended points that beat the frozen one by more than %.2f: %d.\n' % (
                  c['n'], c['best'][0], c['best'][1], c['best'][2], c['f'], c['fg'], c['fs'], 'the argmax is the frozen point' if c['best'] == (3.0, 1.2, 2.0) else 'the argmax moves', FINDING, c['nover'])]
    return L


def report(raw, ctx):
    L = ['# Two post hoc controls of the astronomy results: the grid edge, and the greedy order', '',
         'Generated by `experiments/exp_astro_controls.py` from `controls_raw.npz` and the stored raws `m2_raw.npz`, `repeat_raw.npz`, `cyclefix_raw.npz`, `baselines_hp2_raw.npz` (read only); every number below is printed '
         'by the script. **Everything here is %s**: both controls answer reviewer objections after the test numbers were seen; nothing changes the frozen configuration or any existing report. The set of runs and the rules '
         'were declared in the script\'s docstring before any run (its sha256 is stored in the raw as `meta|declaration_sha256`: %s; the present docstring %s).' % (
             POSTHOC, str(raw['meta|declaration_sha256'][0])[:16] + '...' if 'meta|declaration_sha256' in raw else 'no run yet',
             'matches it' if 'meta|declaration_sha256' in raw and str(raw['meta|declaration_sha256'][0]) == declaration_sha() else 'does NOT match it' if 'meta|declaration_sha256' in raw else 'n/a'),
         '', 'Method of every number: the frozen scripts\' own functions (`exp_astro_m2.run_path` with the guarded sliding path in place of `imaging.sliding_path`, `astro_common.score`, the 20 x 20 px block bootstrap '
         'with %d resamples, seed 0, stratified within the patches for pooled rows, the same resamples for every method of a table). A path that hit the entry cap (%d) is flagged; the guard activations of every path are printed.' % (X.BOOT, X.MAX_ENTRIES), '']
    if 'meta|first_run_utc' in raw:
        L += ['First run recorded: %s; %s.' % (str(raw['meta|first_run_utc'][0]), '; '.join(str(v) for v in raw['meta|versions'])), '']
    summ = {}
    s2, s1 = c2_section(raw, ctx, summ), c1_section(raw, ctx, summ)
    L += summary_lines(summ) + s2 + s1
    wall = [float(raw[k][0]) for k in sorted(raw) if k.endswith('|wall') and k.startswith('E')]
    if wall:
        L += ['## Run time\n', 'Wall seconds of the %d paths (one core, `nice 10`, shared machine): total %.0f s (%.1f h); control 1: %.0f s; control 2: %.0f s.\n' % (
            len(wall), sum(wall), sum(wall) / 3600.0, sum(float(raw[k][0]) for k in raw if k.endswith('|wall') and k.startswith('E1|')), sum(float(raw[k][0]) for k in raw if k.endswith('|wall') and k.startswith('E2|')))]
    return '\n'.join(L) + '\n'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', action='store_true')
    ap.add_argument('--part', choices=['c1', 'c2', 'all'], default='all')
    ap.add_argument('--first', type=int, default=None, help='run at most this many of the pending runs')
    ap.add_argument('--list', action='store_true')
    ap.add_argument('--from-raw', action='store_true')
    A_ = ap.parse_args()
    if A_.list:
        list_jobs()
    elif A_.run:
        run(A_.part, A_.first)
    elif A_.from_raw:
        with open(MD, 'w') as f:
            f.write(report(CF.load_raw(RAW), load_ctx()))
        print('wrote', MD)
    else:
        ap.error('give --run, --from-raw or --list')


if __name__ == '__main__':
    main()
