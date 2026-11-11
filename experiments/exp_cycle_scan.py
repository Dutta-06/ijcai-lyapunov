"""
SMLM spot check of the enter/exit limit cycle of `imaging.sliding_path`: does the defect act on the localisation frames?
=======================================================================================================================

The astronomy runs hit a limit cycle of the sliding path (`results/astro/cyclefix_plan.md`, `src/core/imaging_guard.py`): a lattice kernel enters on the correlation at its lattice start, the
block fit returns its coefficient 0, the unit leaves on the correlation at its slid centre, its dual goes back to the wall, and it re-enters at dt = 0, at the same entry time, until
`max_entries`. The localisation experiments reach `imaging.sliding_path` through `imaging.fit_frame`-like code (`exp_smlm_iss.run_frame`, `exp_smlm2013.run_tile`) with stopping rules (significance,
`snr_extra`) that may end a cycle after a few entries; a cycle at a significant entrant would waste the frame's remaining entries. This script counts, on the first N test frames of each of the eight
sets (the four 2016 sets, the four 2013 sets), at the tuned (E-optimal) operating point and at the a priori default, how many frames show a double exit at one entry time (the frame is fitted by the
experiments' own functions with the guarded path function in place of `imaging.sliding_path`, and the guard's `barred` is read), and compares the detections of the guarded run with those of the unguarded
run, frame by frame: for every frame with an activation and, as a control of the claim that a path without a double exit is the unguarded path bit for bit, for the first CONTROL frames without one.
It reports activations and any change of detections (frames, emitters); it does not recompute the E / J tables.

    python experiments/exp_cycle_scan.py --run [--sets NAME ...] [--frames 40] [--ops E default]    (adds to results/smlm/cycle_scan_raw.npz; a (set, operating point) whose result exists is not run again)
    python experiments/exp_cycle_scan.py --from-raw                                                  (regenerates results/smlm/cycle_scan.md from the raw alone, byte for byte)

Environment: CYCLE_SCAN_OUT (directory of the raw and the md, default results/smlm), SMLM_DATA, SMLM_DATA_2013 (the data, as the experiments read them).
"""
import argparse
import fcntl
import os
import sys
import time

os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import numpy as np

HERE = os.path.abspath(os.path.dirname(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
sys.path.insert(0, REPO)
import exp_astro_cyclefix as CF               # noqa: E402  (the guard's wrapper with its side channel)
import exp_smlm_iss as ISS                    # noqa: E402
import exp_smlm2013 as S13                    # noqa: E402
import smlm2013 as S                          # noqa: E402
import smlm_common as C                       # noqa: E402

OUT = os.environ.get('CYCLE_SCAN_OUT', os.path.join(REPO, 'results', 'smlm'))
RAW, MD = os.path.join(OUT, 'cycle_scan_raw.npz'), os.path.join(OUT, 'cycle_scan.md')
SETS = [('2016', n) for n in ISS.SETS] + [('2013', n) for n in S13.SETS]
OPS = ('E', 'default')
N_FRAMES, CONTROL = 40, 3
TOL_NM = 1e-6                                 # two detections closer than this are the same detection (bitwise, for practical purposes)
FINE_NM = 10.0                                # and the physical tolerance: two detections closer than this (nm) are the same emitter for the localisation scores (matching radius 250 nm)


def key(year, name, op, field):
    return '%s|%s|%s|%s' % (year, name, op, field)


def compare_dets(a, b, tol=TOL_NM):
    """(identical, n in a, n in b, only in a, only in b) for two detection arrays (k, 4): identical = equal bit for bit; the emitters only in one of them are those with no detection of the other
    within `tol` in x and y (nm)."""
    a, b = np.asarray(a, float).reshape(-1, 4), np.asarray(b, float).reshape(-1, 4)
    if a.shape == b.shape and np.array_equal(a, b):
        return True, len(a), len(b), 0, 0

    def unmatched(p, q):
        if not len(p):
            return 0
        if not len(q):
            return len(p)
        d = np.abs(p[:, None, :2] - q[None, :, :2]).max(axis=2)
        return int(np.sum(d.min(axis=1) > tol))
    return False, len(a), len(b), unmatched(a, b), unmatched(b, a)


def per_frame(records, nframes, per):
    """The side-channel records of consecutive path calls grouped by frame (`per` calls per frame: the tiles): [(entries, exit events, activations, barred (n, 4) = tile, tau, b, a)]."""
    if len(records) != nframes * per:
        raise RuntimeError('%d path calls for %d frames x %d tiles' % (len(records), nframes, per))
    out = []
    for f in range(nframes):
        rs = records[f * per:(f + 1) * per]
        bars = [np.c_[np.full(len(r['barred']), t), r['barred']] for t, r in enumerate(rs) if len(r['barred'])]
        out.append((sum(int(r['counts'][0]) for r in rs), sum(int(r['counts'][1]) for r in rs), sum(int(r['counts'][6]) for r in rs), np.concatenate(bars) if bars else np.zeros((0, 4))))
    return out


def tiles_per_frame(name, shape):
    rec = S13.RECIPES[S13.RECIPE_OF[name]]
    if rec['tile'] is None:
        return 1
    return len(S13.tile_grid(shape[0], *rec['tile'])) * len(S13.tile_grid(shape[1], *rec['tile']))


# ---- the fits: the experiments' own functions ---------------------------------------------------------------------------------------

def setup_2016(name, nframes, raw=None):
    """(frames, calibration, operating points {op: (config name, kappa)}) of a 2016 set: the first `nframes` test frames, the calibration and the E-optimal point of the experiment (`exp_smlm_iss`)."""
    raw = raw if raw is not None else ISS.load_raw()
    stack, gt = ISS.load_set(name)
    a, b = ISS.dev_test_frames(name, 'test', nframes)
    frames = C.to_photons(C.load_frames(stack, a, b))
    o = ISS.choose_op(ISS.op_table(raw, name, gt), 'E')
    return frames, ISS.get_cal(raw, name), {'E': (o['cfg'], float(o['kappa'])), 'default': (ISS.DEFAULT['cfg'], float(ISS.DEFAULT['kappa']))}, (a, b)


def setup_2013(name, nframes, raw=None):
    R = (raw if raw is not None else S13.load_raw())[name]
    a, b = S.PROTOCOL[name]['test']
    b = min(b, a + nframes)
    frames = S.to_photons(S.load_frames(S.data_path(name, 'stack'), a, b))
    rows, ops = S13.select_ops(R, name)
    return frames, S13.get_cal(R), {'E': (ops['E']['cfg'], float(ops['E']['kappa'])), 'default': (S13.DEFAULT['cfg'], float(S13.DEFAULT['kappa']))}, (a, b)


def fit(year, name, frames, cfg, cal, kappa, guard):
    """The experiments' own fit of the frames at (config, kappa), with the guarded path function in place (guard) or the frozen one: (detections per frame (k, 4), [per-frame records of the guard], seconds per frame)."""
    (CF.install if guard else CF.uninstall)()
    CF.TRACE.clear()
    try:
        if year == '2016':
            dets, taus, snr, secs = ISS.run_frames(frames, cfg, cal, dict(rules=[('snr', kappa)]), workers=1)
            d, per = dets[0], 1
        else:
            dets, nent, ncap, secs, cpu = S13.run_frames(name, frames, cfg, cal, workers=1)
            kind = S13.RULE_KIND[S13.RECIPE_OF[name]]
            d, per = dets[(kind, kappa)], tiles_per_frame(name, frames.shape[1:])
        recs = per_frame(list(CF.TRACE), len(frames), per) if guard else None
    finally:
        CF.uninstall()
        CF.TRACE.clear()
    return [np.asarray(x, float).reshape(-1, 4) for x in d], recs, np.asarray(secs, float)


def scan(year, name, op, frames, cal, cfg, kappa, rng_span, log=print):
    """The arrays of one (set, operating point): the guarded fit of all frames, then the unguarded fit of the frames with an activation and of the first CONTROL frames without one."""
    t0 = time.time()
    dg, recs, secs = fit(year, name, frames, cfg, cal, kappa, True)
    act = np.array([r[2] for r in recs])
    comp = list(np.flatnonzero(act > 0)) + list(np.flatnonzero(act == 0)[:CONTROL])
    why = ['activation'] * int(np.sum(act > 0)) + ['control'] * (len(comp) - int(np.sum(act > 0)))
    du, _, secs_u = fit(year, name, frames[comp], cfg, cal, kappa, False) if comp else ([], None, np.zeros(0))
    cmp_n = np.array([compare_dets(dg[f], du[j]) for j, f in enumerate(comp)], int).reshape(-1, 5)
    cmp_fine = np.array([compare_dets(dg[f], du[j], FINE_NM)[3:] for j, f in enumerate(comp)], int).reshape(-1, 2)
    bars = np.concatenate([np.c_[np.full(len(r[3]), f), r[3]] for f, r in enumerate(recs) if len(r[3])]) if act.sum() else np.zeros((0, 5))
    log('%s %s %s (%s, kappa %g): %d frames, %d with an activation (%d in all), %d compared, %d changed, %.0f s' % (
        year, name, op, cfg, kappa, len(frames), int(np.sum(act > 0)), int(act.sum()), len(comp), int(np.sum(cmp_n[:, 0] == 0)) if len(comp) else 0, time.time() - t0))
    return {key(year, name, op, 'cfg'): np.array([cfg, str(kappa)]), key(year, name, op, 'frames'): np.array(rng_span), key(year, name, op, 'act'): act,
            key(year, name, op, 'entries'): np.array([r[0] for r in recs]), key(year, name, op, 'exits'): np.array([r[1] for r in recs]), key(year, name, op, 'barred'): bars,
            key(year, name, op, 'secs'): np.array([secs.mean(), secs_u.mean() if len(secs_u) else 0.0]), key(year, name, op, 'cmp_frames'): np.array(comp, int), key(year, name, op, 'cmp_why'): np.array(why or ['']),
            key(year, name, op, 'cmp'): cmp_n, key(year, name, op, 'cmp10'): cmp_fine}


def load_raw(path=RAW):
    with np.load(path, allow_pickle=False) as z:
        return {k: z[k] for k in z.files}


def save_raw(path, arrays):
    """Merge arrays into the raw (atomically, under a lock: two scans of different sets may save into the same raw)."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path + '.lock', 'a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        old = load_raw(path) if os.path.exists(path) else {}
        old.update({k: np.asarray(v) for k, v in arrays.items()})
        tmp = '%s.tmp%d.npz' % (path, os.getpid())
        np.savez_compressed(tmp, **old)
        os.replace(tmp, path)


def run(sets=None, nframes=N_FRAMES, ops=OPS, out_path=RAW, log=print):
    """Each (set, operating point) once; a result that exists is not run again."""
    for year, name in SETS:
        if sets and name not in sets:
            continue
        cur = load_raw(out_path) if os.path.exists(out_path) else {}
        todo = [op for op in ops if key(year, name, op, 'act') not in cur]
        if not todo:
            log('%s %s: results exist' % (year, name))
            continue
        frames, cal, points, span = (setup_2016 if year == '2016' else setup_2013)(name, nframes)
        for op in todo:
            cfg, kappa = points[op]
            save_raw(out_path, scan(year, name, op, frames, cal, cfg, kappa, span, log))


# ---- the report ----------------------------------------------------------------------------------------------------------------------

def table(header, rows):
    return '\n'.join(['| ' + ' | '.join(header) + ' |', '|' + '|'.join(['---'] * len(header)) + '|'] + ['| ' + ' | '.join(str(c) for c in r) + ' |' for r in rows])


def report(raw):
    L = ['# SMLM spot check of the enter/exit limit cycle of `imaging.sliding_path`', '',
         'Generated by `experiments/exp_cycle_scan.py` from `cycle_scan_raw.npz`; every number below is printed by the script. The defect (`results/astro/cyclefix_plan.md`): a lattice kernel enters on the correlation at its lattice start, '
         'the block fit returns its coefficient 0, the unit leaves on the correlation at its slid centre, and it re-enters at dt = 0, at the same entry time, until the entry cap. The guard (`src/core/imaging_guard.py`) bars a lattice point '
         'that has exited twice at the same entry time. Here the frames are fitted by the experiments\' own functions (`exp_smlm_iss.run_frames`, `exp_smlm2013.run_frames`) at the tuned (E-optimal) operating point of the experiment and at '
         'its a priori default, with the guarded path function in place of `imaging.sliding_path`; **an activation** is a lattice point barred in a frame (a tile for the tiled 2013 sets). For every frame with an activation and for '
         'the first %d frames without one (a control: a path without a double exit at one entry time is the unguarded path bit for bit) the frame is fitted again with the frozen, unguarded function and the detections are compared. '
         'Detections: x, y (nm), photons, width (px) at the operating point\'s rule; "changed" means not equal bit for bit. The E and J tables are not recomputed.' % CONTROL, '']
    rows, det_rows, bar_rows = [], [], []
    for year, name in SETS:
        for op in OPS:
            if key(year, name, op, 'act') not in raw:
                continue
            act, ent, ex = raw[key(year, name, op, 'act')], raw[key(year, name, op, 'entries')], raw[key(year, name, op, 'exits')]
            cfg, kap = (str(v) for v in raw[key(year, name, op, 'cfg')])
            comp, why, cmp = raw[key(year, name, op, 'cmp_frames')], raw[key(year, name, op, 'cmp_why')], raw[key(year, name, op, 'cmp')]
            ncomp = len(comp)
            changed = [j for j in range(ncomp) if cmp[j, 0] == 0]
            ctrl_changed = [j for j in changed if str(why[j]) == 'control']
            sec = raw[key(year, name, op, 'secs')]
            rows.append([year, name, 'E-optimal' if op == 'E' else 'a priori default', '%s, kappa %s' % (cfg, kap), '%d' % len(act), '%d' % int(np.sum(act > 0)), '%d' % act.sum(), '%.1f' % ent.mean(), '%d' % ex.sum(),
                         '%d' % ncomp, '%d' % len(changed), '%d' % len(ctrl_changed), '%.2f' % sec[0]])
            fine = raw.get(key(year, name, op, 'cmp10'))
            for j in changed:
                f = int(comp[j])
                det_rows.append([year, name, op, '%d' % (int(raw[key(year, name, op, 'frames')][0]) + f), str(why[j]), '%d' % cmp[j, 1], '%d' % cmp[j, 2], '%d' % cmp[j, 3], '%d' % cmp[j, 4],
                                 '%d' % fine[j, 0] if fine is not None else '-', '%d' % fine[j, 1] if fine is not None else '-'])
            bars = raw[key(year, name, op, 'barred')]
            for f in np.flatnonzero(act > 0)[:10]:
                bb = bars[bars[:, 0] == f]
                bar_rows.append([year, name, op, '%d' % (int(raw[key(year, name, op, 'frames')][0]) + f), '%d' % act[f], '; '.join('%.5g (%d, %d)%s' % (r[2], r[3], r[4], ' tile %d' % r[1] if r[1] else '') for r in bb[:6])
                                 + ('; ...' if len(bb) > 6 else '')])
    L += ['## Activations and changes of detections, per set and operating point\n', table(['year', 'set', 'operating point', 'configuration', 'frames', 'frames with an activation', 'activations',
                                                                                          'entries per frame', 'exit events', 'frames compared', 'changed', 'changed among the controls', 'seconds per frame'], rows), '']
    n_act = sum(int(np.sum(raw[k] > 0)) for k in raw if k.endswith('|act'))
    n_frames = sum(len(raw[k]) for k in raw if k.endswith('|act'))
    n_chg = sum(int(np.sum(raw[k[:-len('|cmp_why')] + '|cmp'][:, 0] == 0)) for k in raw if k.endswith('|cmp_why') and len(raw[k[:-len('|cmp_why')] + '|cmp']))
    n_ctrl = sum(int(np.sum(raw[k] == 'control')) for k in raw if k.endswith('|cmp_why'))
    L += ['**%d frames scanned (all sets and both operating points), %d of them with at least one activation, %d frames whose detections changed (%d control frames compared).**\n' % (n_frames, n_act, n_chg, n_ctrl)]
    if det_rows:
        L += ['## Frames whose detections changed\n', 'Emitters detected only by the guarded or only by the unguarded run: with no detection of the other within %g nm (bitwise) and within %g nm (the position at which two detections are the '
              'same emitter for the scores of the localisation experiments, whose matching radius is 250 nm).\n' % (TOL_NM, FINE_NM),
              table(['year', 'set', 'operating point', 'frame (stack index)', 'why compared', 'detections guarded', 'detections unguarded', 'only guarded (%g nm)' % TOL_NM, 'only unguarded (%g nm)' % TOL_NM,
                     'only guarded (%g nm)' % FINE_NM, 'only unguarded (%g nm)' % FINE_NM], det_rows), '']
    else:
        L += ['## Frames whose detections changed\n', 'None among the frames compared.\n']
    if bar_rows:
        L += ['## The frames with an activation (the first 10 of each set and operating point)\n', table(['year', 'set', 'operating point', 'frame (stack index)', 'activations', 'barred points: entry time (b, a)'], bar_rows), '']
    L += ['## Limits\n', '* Only the first %d test frames of each set were scanned (the experiments\' test ranges are 100 to 900 frames), and only the first %d frames without an activation were compared with an unguarded run; for the others the '
          'guard\'s property (no double exit at one entry time: the path is the unguarded one) is relied on, and the controls test it.' % (N_FRAMES, CONTROL),
          '* The guard acts only on repeated entries at one entry time; a frame\'s cycle that advances tau is not seen here.',
          '* For the tiled 2013 sets an activation is counted per tile and the entries per frame are the sums over the tiles.', '']
    return '\n'.join(L) + '\n'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', action='store_true')
    ap.add_argument('--from-raw', action='store_true')
    ap.add_argument('--sets', nargs='*')
    ap.add_argument('--frames', type=int, default=N_FRAMES)
    ap.add_argument('--ops', nargs='*', default=list(OPS), choices=OPS)
    A_ = ap.parse_args()
    if A_.run:
        run(A_.sets, A_.frames, tuple(A_.ops))
    elif A_.from_raw:
        with open(MD, 'w') as f:
            f.write(report(load_raw()))
        print('wrote', MD)
    else:
        ap.error('give --run or --from-raw')


if __name__ == '__main__':
    main()
