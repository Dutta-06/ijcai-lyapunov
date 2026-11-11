"""
Controls for the SMLM comparison of the sliding inverse scale space (ISS) with Picasso: the order of entry, the symmetry of the Picasso tuning, and the development-versus-test gap
==============================================================================================================================================================================

POST HOC. Everything below was decided AFTER the results of `iss.md`, `smlm2013.md` and `picasso_baseline.md` were known, in answer to three objections of simulated reviewers; it
changes none of the existing results, reports, raws or scripts (they are read, never written). The analyses are declared here, before the first run of this script.

Objections
    (a) the imaging algorithm is a greedy find-fit-refit fitter, so the ISS's order of entry (the accumulated dual) may be irrelevant and the SMLM wins may come from the ingredients
        that any greedy method with them would share (the pixel-integrated PSF with a halo, the high-pass, the noise weights, the resolution limit, the width bound, the stopping rule);
    (b) the tuning is asymmetric: the ISS's operating point maximises E AFTER the position offset is estimated on the development frames and removed (`exp_smlm_iss.op_table`), Picasso's
        R2E / R2J maximise E / J BEFORE it and the offset is applied afterwards (the "+c" rows; `picasso_baseline.md`: "we did not examine the effect");
    (c) no development-versus-test gap is reported for either side.

Analyses (A1-A3; the sets are the four 2016 sets of `iss.md` and the four 2013 sets of `smlm2013.md`; 8 sets in all)
    A1  GREEDY CONTROL. The same code with only the order of entry changed: `imaging.sliding_path(order='greedy')` (matching pursuit with sliding: at each step the lattice candidate of the
        largest matched-filter significance v / sqrt(sum W phi^2) of the current residual enters; tau is the step count) instead of `order='dual'` (entry by the integral of the
        correlation). Everything else is the frozen pipeline of `exp_smlm_iss.py` (2016) / `exp_smlm2013.py` (2013): the PSF (core + halo parts), the high-pass, the weights 1 / var from
        the background and the current fit, the block fits, exits, the resolution limit r_ex, the sliding-width bounds, the stopping rule `stop_snr`, the refinement, the width post-filter
        and the position offset. The order is forced by replacing `imaging.sliding_path` for the duration of a run by a wrapper that sets `order` (nothing else) and records the events
        of each path (`OrderForcing`); the pipelines themselves are called unchanged.
        A1a  Greedy at the dual ISS's own development-chosen operating points: E-tuned and J-tuned (configuration, significance kappa, width bound and offset exactly as the dual's
             `choose_op` selected them on the development frames, recomputed from the frozen raws) and at the a priori default (hp5h, kappa 5, no width bound, no offset). NOTHING is
             re-tuned on any frame. Scored on the same test frames as the dual ISS and Picasso; E, J, RMSE with paired block-bootstrap intervals (blocks of 10 frames, 1000 resamples,
             seed 0, the same as the existing reports) for greedy minus dual and greedy minus Picasso (R2E+c for the E point, R2J+c for the J point, R0 for the default).
        A1b  (secondary, declared with A1) Greedy tuned on its own, on the development frames only, with the same search as the dual's (the four configurations x 7 kappas x 3 width bounds,
             and for the 2016 sets the second round around its E-optimal base configuration; offset estimated and removed on the same frames; selection by `choose_op`), scored once on the
             test frames. This answers "any greedy method with these ingredients and the same tuning effort", where A1a answers "the same operating point".
        Validation first: the dual path is re-run through the harness (original settings) and must reproduce a stored ISS row bit for bit (2016: MT0.N1.LD, E-tuned configuration, all
        seven significance levels, 500 test frames; 2013: Bundled_Tubes_Long_Sequence, both configurations, all keys of the stored part); a dual run with the truncated stop below must equal
        it as well.
        Speed-up that does not change a result (checked by the validation above): the 2016 pipeline stops its path 3 entries after the first entrant below significance 3 (`SNR_MIN`,
        `SNR_EXTRA`) although its rule `stop_snr` reads only the leading entrants above kappa >= 4; the runs here stop at the first entrant below 4 (`snr_min = 4`, `snr_extra = 0`; the
        2013 pipeline already does). The detections of every kappa >= 4 are those of the long path.
        Guard: `imaging_guard.sliding_path_guarded` supports `order='greedy'` but its guard keys on the entry time, which advances by 1 at every greedy entry, so it can never act there; a
        check on a sample of frames confirms that the guarded and the unguarded greedy path are identical with no barred point. The runs use the unguarded path, as all stored results do;
        for every path the events are scanned for the condition the guard would act on (a point that exits twice at one entry time) and the numbers are reported with the exit /
        re-entry counts and the paths that reached the cap on entries.
    A2  PICASSO OBJECTIVE SYMMETRY. The stored ladder does not hold what a re-selection needs (only the chosen configurations' detections, not the 108 x ~120 scored combinations), so the
        search is re-run on the DEVELOPMENT frames only, with the SAME grids and the SAME coordinate ascent as `exp_smlm_picasso.py` (2016) / `exp_smlm2013.py` (2013) (108 structural settings
        run once at the lowest net gradient; per setting the coordinate ascent over net gradient, minimum photons, minimum width, maximum width, at most 3 rounds, strict improvements, shared
        cache for the two objectives) but with the objective of the ISS: E (R2Eo) / J (R2Jo) AFTER the constant offset (mean signed error of the matched true positives of that very
        configuration on the development frames) has been removed and the pairing redone, evaluated for every candidate configuration. The structural setting is chosen by the best
        offset-corrected development value (earlier setting on ties). The chosen configuration is run once on the test frames, the development offset subtracted, and scored (the stored
        rows R2E+c / R2J+c are the original selection with the offset appended). Before the search the stored R2E / R2J configurations are re-run on the development frames and must reproduce
        the stored development detections. ISS-minus-Picasso paired differences are reported next to the stored ones, and the verdict (wins / ties / losses by the paired interval) before and
        after. The search does not contain the stored configurations as such (the ascent may land elsewhere) but the same space; the offset-corrected optimum over it is, by construction
        of the ascent, not guaranteed to be at least the stored configuration's corrected value (local search): the corrected development value of the stored configuration is reported next
        to the new one.
    A3  DEV-VS-TEST GAP. For the ISS (dual E-tuned, J-tuned; greedy A1b if run) and Picasso (R2E+c, R2J+c as stored; R2Eo, R2Jo of A2): E and J on the development frames and on the test frames
        of each set, their difference, with a 95 % interval of the difference from independent block bootstraps of the two frame ranges (seeds 1 and 2, 1000 resamples, blocks of 10 frames;
        the two ranges are different frames, so the difference of two noisy values has the noise of both). Development detections of the stored rows come from the frozen raws (ISS: the
        development runs with the chosen kappa, width bound and offset; Picasso: the ladder's development detections); a positive test-minus-development difference is no over-fitting.

Rules for reading the results (fixed now): a paired difference is a win / loss when its 95 % interval excludes 0 (by more than 1e-4 in the metric's units: two paths that fit the same kernels in another order
differ by about 1e-14, and float32 storage of the 2013 detections by about 1e-5 nm, which is rounding, not a result), a tie otherwise; headline metrics are the primary setting (all
ground-truth emitters); the challenge setting is reported for E. Intervals cover the sampling of frames only.

Usage:  python experiments/exp_smlm_controls.py --run [--parts REF G GB PS] [--sets ...] [--frames N] [--force] [--max-structs N] [--cache-dir D] [--picasso-python PY]
            (computes the parts, merges them into controls_raw.npz; the report is written when --run finishes; resumable: a part that is in the raw is skipped unless --force)
        python experiments/exp_smlm_controls.py --from-raw [--boot 1000]       (regenerates controls.md from controls_raw.npz alone, byte for byte)
Environment: CONTROLS_OUT (directory of the raw and the md, default results/smlm), SMLM_PICASSO_PYTHON (an interpreter with picasso for the part PS), SMLM_DATA, SMLM_DATA_2013.
Parts: REF the frozen references (ISS dual rows, Picasso rows, ground truth, operating points) copied into the raw so that the report needs nothing else; G greedy at the dual's points
(A1a) with the validation runs; GB greedy tuned on its own (A1b); PS Picasso with the offset-corrected objective (A2).
"""
import argparse
import hashlib
import json
import os
import sys
import time

os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import numpy as np                                                  # noqa: E402

HERE = os.path.abspath(os.path.dirname(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
sys.path.insert(0, REPO)
import smlm_common as C                                             # noqa: E402
import smlm2013 as S13                                              # noqa: E402
import exp_smlm_iss as ISS                                          # noqa: E402
import exp_smlm_picasso as PIC                                      # noqa: E402
import exp_smlm2013 as E13                                          # noqa: E402
from src.core import imaging as IM                                  # noqa: E402
from src.core import imaging_guard as IG                            # noqa: E402

OUT = os.environ.get('CONTROLS_OUT', os.path.join(REPO, 'results', 'smlm'))
RAW = os.path.join(OUT, 'controls_raw.npz')
MD = os.path.join(OUT, 'controls.md')
FROZEN = os.path.join(REPO, 'results', 'smlm')
ISS_RAW, LADDER, RAW13 = (os.path.join(FROZEN, f) for f in ('iss_raw.npz', 'picasso_ladder_raw.npz', 'smlm2013_raw.npz'))

SETS16, SETS13 = tuple(ISS.SETS), tuple(E13.SETS)
SETS = SETS16 + SETS13
SHORT = {'MT0.N1.LD': 'N1 LD', 'MT0.N1.HD': 'N1 HD', 'MT0.N2.LD': 'N2 LD', 'MT0.N2.HD': 'N2 HD', 'Tubulins_I': 'Tubulins I', 'Tubulins_II': 'Tubulins II',
         'Bundled_Tubes_High_Density': 'Bundled HD', 'Bundled_Tubes_Long_Sequence': 'Bundled LS'}
BOOT, BLOCK = 1000, C.DEFAULT_BLOCK
CHALLENGE = dict(border_nm=450.0)
OPS = ('E', 'J', '0')                                               # the dual's E-tuned and J-tuned points, and the a priori default
TRUNC = dict(snr_min=4.0, snr_extra=0)                              # the shortened stop of the 2016 runs (see the docstring); the 2013 pipeline already stops this way
GAP_SEEDS = (1, 2)                                                  # independent bootstraps of the development and of the test frames (A3)
PARTS = ('REF', 'G', 'GB', 'PS')

LABEL = {'ISS-E': 'dual ISS, E-tuned', 'ISS-J': 'dual ISS, J-tuned', 'ISS-0': 'dual ISS, a priori default',
         'GR-E': 'greedy at the E-tuned point', 'GR-J': 'greedy at the J-tuned point', 'GR-0': 'greedy at the a priori default',
         'GT-E': 'greedy tuned on its own (E)', 'GT-J': 'greedy tuned on its own (J)',
         'R0': 'Picasso R0', 'R2E+c': 'Picasso R2E+c (stored)', 'R2J+c': 'Picasso R2J+c (stored)', 'R2Eo': 'Picasso R2Eo (offset-corrected objective)',
         'R2Jo': 'Picasso R2Jo (offset-corrected objective)'}
PIC_OF = {'E': 'R2E+c', 'J': 'R2J+c', '0': 'R0'}                    # the like-for-like Picasso row of each operating point
PICO_OF = {'E': 'R2Eo', 'J': 'R2Jo'}


def fam(name):
    return '16' if name in SETS16 else '13'


# ==================================================================================================================================
# raw storage: one deterministic npz, parts merged by key prefix
# ==================================================================================================================================
RAW_README = ('Raw results of experiments/exp_smlm_controls.py (post hoc controls; the analyses are declared in the script\'s docstring). Keys "<part>|<set>|...": ref (the frozen references copied '
              'from iss_raw.npz, smlm2013_raw.npz and picasso_ladder_raw.npz: ground truth of the development and test frames, the dual ISS\'s operating points, the scored detections of the '
              'dual ISS, Picasso R0 / R2E+c / R2J+c as "<method>|<split>|off" (F + 1 offsets) and "<method>|<split>|flat" (M x 4: x_nm, y_nm, photons, width)), g (the greedy path at the '
              'dual\'s points: per configuration the detections of every significance rule before the width filter and the offset, "<cfg>|<rule>|off|flat", and the path statistics), gb (the '
              'greedy path tuned on its own: development and test detections of the configurations, its operating points), ps (Picasso with the offset-corrected objective: the search, the chosen '
              'configurations, the offset-corrected detections), val (the validation runs), meta (versions, digests of the frozen files). results/smlm/controls.md is generated from this file alone '
              '(python experiments/exp_smlm_controls.py --from-raw).')


def load_raw(path=None):
    path = path or RAW
    if not os.path.exists(path):
        return {}
    return E13.load_npz(path)


def save_part(path, prefix, arrays, replace=True):
    """Merge ``arrays`` into the raw file under the key prefix ``prefix + '|'`` (existing keys of that prefix are replaced); deterministic bytes (sorted keys, fixed zip date)."""
    old = load_raw(path)
    if replace:
        old = {k: v for k, v in old.items() if not k.startswith(prefix + '|')}
    old['README'] = np.array(RAW_README)
    old.update({prefix + '|' + k: np.asarray(v) for k, v in arrays.items()})
    E13.save_npz(path, old)


def jdump(o):
    return np.array(json.dumps(o, sort_keys=True, default=C._json_default))


def jload(a):
    return json.loads(str(a))


def pack(dets):
    return ISS.pack(dets)


def unpack(raw, key):
    return ISS.unpack(raw[key + '|off'], raw[key + '|flat'])


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


# ==================================================================================================================================
# the order of entry, forced from outside
# ==================================================================================================================================
def path_stats(events, n_entries, max_entries):
    """Facts of one path from its event list [(tau, (b, a), 'enter' | 'exit')]: entries, exits, re-entries (an entry of a point that had exited before), guard events (an exit of a point that
    had exited at the very same entry time before: the condition `imaging_guard` bars on) and whether the path reached the cap on entries. Returns the vector (calls, entries, exits,
    reentries, guard events, paths with a guard event, paths at the cap)."""
    exited, last, guard, exits, reent = set(), {}, 0, 0, 0
    for tau, pid, kind in events:
        if kind == 'enter':
            reent += int(pid in exited)
        else:
            exits += 1
            if pid in last and abs(tau - last[pid]) <= 1e-12 * max(1.0, abs(tau)):
                guard += 1
            last[pid] = tau
            exited.add(pid)
    return np.array([1, n_entries, exits, reent, guard, int(guard > 0), int(n_entries >= max_entries)], dtype=np.int64)


STAT_NAMES = ('paths', 'entries', 'exits', 're-entries', 'guard events', 'paths with a guard event', 'paths at the cap')


class OrderForcing:
    """Context manager: for its duration `imaging.sliding_path` (as called by the frozen pipelines through the module) forces ``order`` and the keyword overrides ``override``, and accumulates
    ``path_stats`` of every path in ``self.stats``. ``guarded=True`` calls `imaging_guard.sliding_path_guarded` (cycle_guard on) instead and counts its barred points in ``self.barred``.
    The original function is restored on exit, also after an exception."""

    def __init__(self, order, override=None, guarded=False):
        self.order, self.override, self.guarded = order, dict(override or {}), guarded
        self.stats, self.barred = np.zeros(len(STAT_NAMES), dtype=np.int64), 0
        self._orig = None

    def __enter__(self):
        self._orig = IM.sliding_path
        IM.sliding_path = self._call
        return self

    def __exit__(self, *exc):
        IM.sliding_path = self._orig
        return False

    def _call(self, *args, **kw):
        kw['order'] = self.order
        kw.update(self.override)
        R = (IG.sliding_path_guarded(*args, cycle_guard=True, **kw) if self.guarded else self._orig(*args, **kw))
        if self.guarded:
            self.barred += len(R['barred'])
        self.stats += path_stats(R['events'], len(R['taus']), kw.get('max_entries', 40))
        return R


# ==================================================================================================================================
# scoring helpers (both families)
# ==================================================================================================================================
def score_set(name, dets, gt, idx, pmin=None):
    """Score with the field of the set; pmin switches on the challenge setting (border 450 nm and the photon threshold)."""
    kw = dict(border_nm=CHALLENGE['border_nm'], photon_min=pmin) if pmin is not None else {}
    return C.score(dets, gt, idx, **kw) if fam(name) == '16' else S13.score(dets, gt, idx, name, **kw)


def equal_dets(a, b):
    return len(a) == len(b) and all(np.array_equal(x, y) for x, y in zip(a, b))


def ops_json(ops):
    return {w: dict(cfg=o['cfg'], kappa=float(o['kappa']), smax=float(o['smax']), off=[float(v) for v in o['off']], m=[float(v) for v in o['m']]) for w, o in ops.items()}


def cut(dets, n):
    return dets if n is None else dets[:n]


# ==================================================================================================================================
# part REF: the frozen references
# ==================================================================================================================================
def ref_keys(name):
    return 'ref|%s|' % name


def iss_ops16(iss, name, gt):
    rows = ISS.op_table(iss, name, gt)
    return rows, {w: ISS.choose_op(rows, w) for w in ISS.OPS}


def part_ref(name, n_frames=None, log=print):
    """The frozen references of a set, copied into the raw: ground truth of the development and test frames, the dual ISS's operating points (recomputed from its development runs with the
    frozen selection code) and the scored detections of the dual ISS (E, J, default) and Picasso (R0, R2E+c, R2J+c) on the test frames, of E, J and R2E+c / R2J+c on the development frames."""
    t0 = time.time()
    A = {}
    if fam(name) == '16':
        iss = ISS.load_raw(ISS_RAW)
        gt = C.load_gt(C.data_path(name, 'gt'))
        pmin = C.photon_threshold(gt, 0.75)
        (d0, d1), (a, b) = ISS.PROTOCOL[name]['dev'], ISS.PROTOCOL[name]['test']
        b = min(b, a + n_frames) if n_frames else b
        _, ops = iss_ops16(iss, name, gt)
        A.update(kind=np.array('snr'), pmin=np.array(pmin), frames=np.array([d0, d1, a, b]), gt_dev=gt[(gt['frame'] >= d0) & (gt['frame'] < d1)], gt_test=gt[(gt['frame'] >= a) & (gt['frame'] < b)])
        k0 = 'T0|%s|run|snr%g|' % (name, ISS.DEFAULT['kappa'])
        test = {'ISS-E': ISS.op_test_dets(iss, name, ops['E']), 'ISS-J': ISS.op_test_dets(iss, name, ops['J']), 'ISS-0': ISS.unpack(iss[k0 + 'off'], iss[k0 + 'flat'])}
        dev = {}
        for w in ISS.OPS:
            o = ops[w]
            key = 'D|%s|%s|' % (name, o['cfg'])
            dev['ISS-' + w] = ISS.shift(ISS.post(ISS.unpack(iss[key + '%g|off' % o['kappa']], iss[key + '%g|flat' % o['kappa']]), o['smax']), o['off'])
        for r in ('R0', 'R2E+c', 'R2J+c'):
            test[r] = C.load_ladder_detections(LADDER, name, r, 'test')[0]
        for r in ('R2E+c', 'R2J+c'):
            dev[r] = C.load_ladder_detections(LADDER, name, r, 'dev')[0]
    else:
        R = E13.load_raw(RAW13)[name]
        rows, ops = E13.select_ops(R, name)
        kind = E13.RULE_KIND[E13.RECIPE_OF[name]]
        (d0, d1), (a, b) = (int(v) for v in R['meta|frame_ranges'][0:2]), (int(v) for v in R['meta|frame_ranges'][2:4])
        b = min(b, a + n_frames) if n_frames else b
        gt_dev, gt_test = R['meta|gt_dev'], R['meta|gt_test']
        A.update(kind=np.array(kind), pmin=np.array(float(R['meta|pmin'])), frames=np.array([d0, d1, a, b]), gt_dev=gt_dev, gt_test=S13.gt_in_frames(gt_test, a, b))
        test = {'ISS-E': E13.op_test_dets(R, name, ops['E']), 'ISS-J': E13.op_test_dets(R, name, ops['J']), 'ISS-0': E13.is_dets(R, 'test', E13.DEFAULT['cfg'], kind, E13.DEFAULT['kappa'])}
        dev = {}
        for w in E13.OPS:
            o = ops[w]
            dev['ISS-' + w] = ISS.shift(E13.post(E13.is_dets(R, 'dev', o['cfg'], kind, o['kappa']), o['smax']), o['off'])
        idx_dev = np.arange(d0, d1)
        for r in ('R0', 'R2E', 'R2J'):
            dd, dt = E13.unpack(R['pic|dev|%s|off' % r], R['pic|dev|%s|flat' % r]), E13.unpack(R['pic|test|%s|off' % r], R['pic|test|%s|flat' % r])
            off = E13._off(dd, gt_dev, idx_dev, name)
            if r == 'R0':
                test['R0'] = dt
            else:
                test[r + '+c'], dev[r + '+c'] = C.shift_detections(dt, *off), C.shift_detections(dd, *off)
    A['ops'] = jdump(ops_json(ops))
    for m, ds in test.items():
        A['%s|test|off' % m], A['%s|test|flat' % m] = pack(cut(ds, b - a))
    for m, ds in dev.items():
        A['%s|dev|off' % m], A['%s|dev|flat' % m] = pack(ds)
    save_part(RAW, 'ref|' + name, A)
    log('REF %s: operating points E %s kappa %g, J %s kappa %g; %d test frames  [%.0f s]' % (name, ops['E']['cfg'], ops['E']['kappa'], ops['J']['cfg'], ops['J']['kappa'], b - a, time.time() - t0))


def part_meta(log=print):
    """Digests of the frozen files the references were copied from, and the versions."""
    import platform
    d = {os.path.basename(p): sha256(p) for p in (ISS_RAW, LADDER, RAW13) if os.path.exists(p)}
    save_part(RAW, 'meta', dict(sha=jdump(d), versions=jdump(dict(python=platform.python_version(), numpy=np.__version__))))
    log('meta: digests of %s' % ', '.join(sorted(d)))


# ==================================================================================================================================
# part G: the greedy path at the dual's operating points (A1a), with the validation runs
# ==================================================================================================================================
def get_ops(name, raw=None):
    """The dual ISS's operating points of a set as stored by REF (cfg, kappa, smax, off, dev metrics)."""
    raw = raw if raw is not None else load_raw()
    return jload(raw['ref|%s|ops' % name])


def cfgs_needed(ops, default_cfg):
    return sorted({ops['E']['cfg'], ops['J']['cfg'], default_cfg})


def run_cfg16(name, frames, cfg, cal, order, override=None, guarded=False):
    """The 2016 pipeline on ``frames`` for one configuration with the order forced; returns (arrays of detections of every rule, stats, seconds per frame, barred)."""
    with OrderForcing(order, override, guarded) as of:
        dets, _, _, secs = ISS.run_frames(frames, cfg, cal, dict(rules=ISS.rules_dev()), workers=1)
    arrays = {}
    ISS._store(arrays, cfg, ISS.rules_dev(), dets)
    return arrays, of.stats, secs, of.barred


def validate16(name, frames, ops, cal, iss, log):
    """Bit-for-bit: the dual path through this harness (original settings, and the truncated stop) against the stored T part of the E-tuned configuration, all seven rules."""
    cfg = ops['E']['cfg']
    out = {}
    for tag, order, ov in (('dual_original_stop', None, None), ('dual_truncated_stop', 'dual', TRUNC)):
        t0 = time.time()
        if order is None:
            dets, _, _, _ = ISS.run_frames(frames, cfg, cal, dict(rules=ISS.rules_dev()), workers=1)
            arr = {}
            ISS._store(arr, cfg, ISS.rules_dev(), dets)
        else:
            arr, _, _, _ = run_cfg16(name, frames, cfg, cal, order, ov)
        ok, nfr = [], 0
        for k in ISS.KAPPAS:
            sk = 'snr%g' % k
            ref_off, ref_flat = iss['T|%s|%s|%s|off' % (name, cfg, sk)], iss['T|%s|%s|%s|flat' % (name, cfg, sk)]
            n = len(frames)
            mine = ISS.unpack(arr['%s|%s|off' % (cfg, sk)], arr['%s|%s|flat' % (cfg, sk)])
            ref = ISS.unpack(ref_off, ref_flat)[:n]
            ok.append(equal_dets(mine, ref))
            nfr = len(ref)
        out[tag] = dict(set=name, cfg=cfg, frames=nfr, rules=len(ok), identical=bool(all(ok)), rules_identical=int(sum(ok)))
        log('VAL %s %s %s: %d frames, %d of %d rules bit-identical to the stored T part  [%.0f s]' % (name, cfg, tag, nfr, sum(ok), len(ok), time.time() - t0))
    return out


def validate_guard(name, frames, ops, cal, n, log):
    """Greedy path, guarded against unguarded, on the first n test frames: identical detections and no barred point."""
    cfg = ops['E']['cfg']
    a, st1, _, _ = run_cfg16(name, frames[:n], cfg, cal, 'greedy', TRUNC)
    b, st2, _, barred = run_cfg16(name, frames[:n], cfg, cal, 'greedy', TRUNC, guarded=True)
    same = all(np.array_equal(a[k], b[k]) for k in a)
    log('VAL guard %s %s: %d frames, guarded greedy path identical to the unguarded: %s, barred points %d' % (name, cfg, n, same, barred))
    return dict(set=name, cfg=cfg, frames=int(n), identical=bool(same), barred=int(barred))


def part_g16(name, n_frames=None, log=print, validate=True):
    raw = load_raw()
    ops = get_ops(name, raw)
    iss = ISS.load_raw(ISS_RAW)
    cal = ISS.get_cal(iss, name)
    stack = C.data_path(name, 'stack')
    a, b = ISS.dev_test_frames(name, 'test', n_frames)
    frames = C.to_photons(C.load_frames(stack, a, b))
    A = {'frames': np.array([a, b])}
    cfgs = cfgs_needed(ops, ISS.DEFAULT['cfg'])
    for cfg in cfgs:
        t0 = time.time()
        arr, stats, secs, _ = run_cfg16(name, frames, cfg, cal, 'greedy', TRUNC)
        A.update(arr)
        A['%s|stats' % cfg], A['%s|secs' % cfg] = stats, secs
        log('G %s greedy %s: %d frames, %.0f s (%.2f s/frame); paths %d, exits %d, re-entries %d, guard events %d, at the cap %d' % (name, cfg, b - a, time.time() - t0, secs.mean(), stats[0], stats[2], stats[3], stats[4], stats[6]))
    # the dual's own path statistics (for the same configurations, truncated stop; the detections equal the stored ones, see the validation)
    for cfg in cfgs:
        _, stats, secs, _ = run_cfg16(name, frames, cfg, cal, 'dual', TRUNC)
        A['%s|dual_stats' % cfg], A['%s|dual_secs' % cfg] = stats, secs
        log('G %s dual %s: paths %d, exits %d, re-entries %d, guard events %d, at the cap %d' % (name, cfg, stats[0], stats[2], stats[3], stats[4], stats[6]))
    save_part(RAW, 'g|' + name, A)
    if validate:
        V = {}
        if name == 'MT0.N1.LD':
            V.update(validate16(name, frames, ops, cal, iss, log))
        if name == 'MT0.N1.HD':
            V['guard'] = validate_guard(name, frames, ops, cal, min(20, b - a), log)
        if V:
            save_part(RAW, 'val|' + name, {k: jdump(v) for k, v in V.items()})


def part_g13(name, n_frames=None, log=print, validate=True):
    raw = load_raw()
    ops = get_ops(name, raw)
    R = E13.load_raw(RAW13)[name]
    a, b = (int(v) for v in R['meta|frame_ranges'][2:4])
    cfgs = cfgs_needed(ops, E13.DEFAULT['cfg'])
    A = {k[len('g|%s|' % name):]: v for k, v in raw.items() if k.startswith('g|%s|' % name) and k != 'g|%s|frames' % name}       # configurations finished before an interruption are kept
    for cfg in cfgs:
        if '%s|stats' % cfg in A:
            continue
        t0 = time.time()
        with OrderForcing('greedy') as of:
            arr = E13.part_is(name, 'test', cfg, R, 1, log)
        arr.pop('bg_test', None)
        for k, v in arr.items():
            A['%s|%s' % (cfg, k)] = v
        A['%s|stats' % cfg] = of.stats
        save_part(RAW, 'g|' + name, A)                                      # no 'frames' key yet: the set counts as done only at the end
        log('G %s greedy %s: %.0f s; tile paths %d, exits %d, re-entries %d, guard events %d, at the cap %d' % (name, cfg, time.time() - t0, of.stats[0], of.stats[2], of.stats[3], of.stats[4], of.stats[6]))
    if validate and name == 'Bundled_Tubes_Long_Sequence':                  # the cheap 2013 set: the dual path through this harness against the whole stored part
        for cfg in cfgs:
            with OrderForcing('dual') as of:
                arr = E13.part_is(name, 'test', cfg, R, 1, log)
            same = all(np.array_equal(arr[k], R['istest.%s|%s' % (cfg, k)]) for k in arr if k not in ('bg_test', 'secs', 'cpu'))
            A['%s|dual_identical' % cfg], A['%s|dual_stats' % cfg] = np.array(int(same)), of.stats
            log('VAL %s %s: dual path through this harness identical to the stored part: %s' % (name, cfg, same))
    A['frames'] = np.array([a, b])
    save_part(RAW, 'g|' + name, A)


# ==================================================================================================================================
# part GB: the greedy path tuned on its own (A1b)
# ==================================================================================================================================
def part_gb16(name, n_frames=None, log=print):
    """Greedy, the development search of the dual (four configurations, second round around the E-optimal base), the operating points by ``choose_op``, the test frames once."""
    iss = ISS.load_raw(ISS_RAW)
    cal = ISS.get_cal(iss, name)
    stack = C.data_path(name, 'stack')
    gt = C.load_gt(C.data_path(name, 'gt'))
    (d0, d1), (a, b) = ISS.dev_test_frames(name, 'dev'), ISS.dev_test_frames(name, 'test', n_frames)
    dev = C.to_photons(C.load_frames(stack, d0, d1))
    A = {'frames': np.array([d0, d1, a, b])}
    view = {}

    def run_dev(cfg):
        t0 = time.time()
        arr, stats, secs, _ = run_cfg16(name, dev, cfg, cal, 'greedy', TRUNC)
        for r in ISS.rules_dev():
            for s in ('off', 'flat'):
                A['D|%s|snr%g|%s' % (cfg, r[1], s)] = arr['%s|snr%g|%s' % (cfg, r[1], s)]
                view['D|%s|%s|%g|%s' % (name, cfg, r[1], s)] = arr['%s|snr%g|%s' % (cfg, r[1], s)]
        view['D|%s|%s|frames' % (name, cfg)] = np.array([d0, d1])
        log('GB %s dev greedy %s: %d frames, %.0f s' % (name, cfg, d1 - d0, time.time() - t0))
    for cfg in ISS.CONFIGS:
        run_dev(cfg)
    rows = [r for r in ISS.op_table(view, name, gt) if r['cfg'] in ISS.CONFIGS]
    base = ISS.choose_op(rows, 'E')['cfg']
    for suf in ISS.SUFFIX:
        run_dev(base + suf)
    rows = ISS.op_table(view, name, gt)
    ops = {w: ISS.choose_op(rows, w) for w in ISS.OPS}
    A['ops'] = jdump(ops_json(ops))
    test = C.to_photons(C.load_frames(stack, a, b))
    for cfg in sorted({o['cfg'] for o in ops.values()}):
        t0 = time.time()
        arr, stats, secs, _ = run_cfg16(name, test, cfg, cal, 'greedy', TRUNC)
        for k, v in arr.items():
            A['T|' + k] = v
        log('GB %s test greedy %s: %.0f s' % (name, cfg, time.time() - t0))
    save_part(RAW, 'gb|' + name, A)
    log('GB %s own points: E %s kappa %g width <= %g, J %s kappa %g width <= %g' % (name, ops['E']['cfg'], ops['E']['kappa'], ops['E']['smax'], ops['J']['cfg'], ops['J']['kappa'], ops['J']['smax']))


def part_gb13(name, n_frames=None, log=print):
    R = E13.load_raw(RAW13)[name]
    view = {k: v for k, v in R.items() if not k.startswith('isdev.')}
    A = {}
    (d0, d1), (a, b) = (int(v) for v in R['meta|frame_ranges'][0:2]), (int(v) for v in R['meta|frame_ranges'][2:4])
    A['frames'] = np.array([d0, d1, a, b])
    for cfg in E13.CONFIGS:
        t0 = time.time()
        with OrderForcing('greedy'):
            arr = E13.part_is(name, 'dev', cfg, R, 1, log)
        arr.pop('bg_test', None)
        for k, v in arr.items():
            A['D|%s|%s' % (cfg, k)] = v
            view['isdev.%s|%s' % (cfg, k)] = v
        log('GB %s dev greedy %s: %.0f s' % (name, cfg, time.time() - t0))
    rows = E13.op_table(view, name)
    ops = {w: ISS.choose_op(rows, w) for w in E13.OPS}
    A['ops'] = jdump(ops_json(ops))
    for cfg in sorted({o['cfg'] for o in ops.values()}):
        t0 = time.time()
        with OrderForcing('greedy'):
            arr = E13.part_is(name, 'test', cfg, R, 1, log)
        arr.pop('bg_test', None)
        for k, v in arr.items():
            A['T|%s|%s' % (cfg, k)] = v
        log('GB %s test greedy %s: %.0f s' % (name, cfg, time.time() - t0))
    save_part(RAW, 'gb|' + name, A)
    log('GB %s own points: E %s kappa %g width <= %.3f, J %s kappa %g width <= %.3f' % (name, ops['E']['cfg'], ops['E']['kappa'], ops['E']['smax'], ops['J']['cfg'], ops['J']['kappa'], ops['J']['smax']))


# ==================================================================================================================================
# part PS: Picasso with the offset-corrected objective (A2)
# ==================================================================================================================================
def digest_dets(dets):
    h = hashlib.blake2b(digest_size=16)
    h.update(np.array([len(d) for d in dets], dtype=np.int64).tobytes())
    for d in dets:
        if len(d):
            h.update(np.ascontiguousarray(d, dtype=np.float64).tobytes())
    return h.digest()


def offset_corrected_metrics(dets, gt, idx):
    """The ISS's objective for a set of detections: the constant offset (mean signed error of the matched pairs) is estimated on these very frames, subtracted, the pairing redone, and the
    metrics of the shifted detections returned with the offset (``exp_smlm_iss.op_table``: ``offset_of``, ``shift``, ``score``)."""
    r0 = C.score(dets, gt, idx)
    off = np.nan_to_num(np.array([r0['bias_x_nm'], r0['bias_y_nm']]))
    r = C.score(C.shift_detections(dets, off[0], off[1]), gt, idx)
    m = {k: r[k] for k in C._TUNE_KEYS}
    m['off_x'], m['off_y'] = float(off[0]), float(off[1])
    return m


def ascent_offset(candidates, gt, idx, grids, objective, cache, dcache, max_rounds=3, evaluator=None):
    """`smlm_common.tune_postfilters` (the same coordinate ascent: grids in order, strict improvements, at most ``max_rounds`` rounds, NaN as -inf) with the objective evaluated AFTER the
    offset correction. ``cache`` (filter combination -> metrics) is shared by the two objectives; ``dcache`` (digest of the filtered detections -> metrics) skips the scoring of combinations that
    keep the same detections (a speed-up that changes no value). ``evaluator(dets, gt, idx)`` is the metric function (default ``offset_corrected_metrics``; with ``plain_metrics`` it is
    `tune_postfilters` itself)."""
    ev_fn = evaluator or offset_corrected_metrics
    keys = list(grids)
    cur = {k: None for k in keys}

    def ev(f):
        key = tuple(f.get(k) for k in keys)
        if key not in cache:
            dets = C.filter_detections(candidates, f)
            dg = digest_dets(dets)
            if dg not in dcache:
                dcache[dg] = ev_fn(dets, gt, idx)
            cache[key] = dcache[dg]
        return cache[key]

    def val(m):
        v = m[objective]
        return -np.inf if v is None or (isinstance(v, float) and np.isnan(v)) else float(v)
    best, first, rounds, changed = val(ev(cur)), None, 0, True
    for rnd in range(int(max_rounds)):
        rounds = rnd + 1
        changed = False
        for ki, k in enumerate(keys):
            for v in grids[k]:
                if v == cur[k]:
                    continue
                q = dict(cur)
                q[k] = v
                x = val(ev(q))
                if x > best + 1e-12:
                    best, cur, changed = x, q, True
            if rnd == 0 and ki == 0:
                first = dict(filters=dict(cur), value=best, metrics=dict(ev(cur)))
        if not changed:
            break
    return dict(filters=dict(cur), value=best, metrics=dict(ev(cur)), rounds=rounds, converged=not changed, first=first, n_eval=len(cache))


def pic_setup(name):
    """What the Picasso search of a set needs: frames (DN), ground truth, indices, grids, the lowest net gradient, the camera options, the stored configurations."""
    if fam(name) == '16':
        noise = PIC.noise_of(name)
        (d0, d1), (t0, t1) = PIC.PROTOCOL[name]['dev'], PIC.PROTOCOL[name]['test']
        stack = C.data_path(name, 'stack')
        gt = C.load_gt(C.data_path(name, 'gt'))
        z = np.load(LADDER, allow_pickle=False)
        stored = {r: json.loads(str(z['%s|%s|config' % (name, r)])) for r in ('R2E', 'R2J')}
        return dict(dev=C.load_frames(stack, d0, d1), test=C.load_frames(stack, t0, t1), d=(d0, d1), t=(t0, t1), gt=gt, grids=PIC.FILTER_GRIDS[noise], ng_low=PIC.ng_low(noise), cam={}, stored=stored,
                    stored_dev={r: C.load_ladder_detections(LADDER, name, r, 'dev')[0] for r in ('R2E', 'R2J')})
    R = E13.load_raw(RAW13)[name]
    (d0, d1), (t0, t1) = S13.PROTOCOL[name]['dev'], S13.PROTOCOL[name]['test']
    stack = S13.data_path(name, 'stack')
    params = jload(R['pic|params'])
    return dict(dev=S13.load_frames(stack, d0, d1), test=S13.load_frames(stack, t0, t1), d=(d0, d1), t=(t0, t1), gt=S13.load_gt(S13.data_path(name, 'gt')), grids=E13.FILTER_GRID,
                ng_low=E13.PIC_NG_LOW, cam=E13.camera(name), stored={r: params[r] for r in ('R2E', 'R2J')},
                stored_dev={r: E13.unpack(R['pic|dev|%s|off' % r], R['pic|dev|%s|flat' % r]) for r in ('R2E', 'R2J')})


def search_offset(session, gt_dev, dev_idx, cam, ng_low, grids, structs, cache_file=None, sig='', max_structs=None, log=print):
    """A2: every structural setting once at the lowest net gradient, the coordinate ascent with the offset-corrected objective for E and for J (one shared cache). Development frames only.
    Resumable through ``cache_file`` (a JSON written after every setting; ``sig`` guards it)."""
    out = {}
    if cache_file and os.path.exists(cache_file):
        blob = json.load(open(cache_file))
        if blob.get('signature') == sig:
            out = blob['results']
            log('  search results of %d settings read from %s' % (len(out), cache_file))
    t_all = time.time()
    for i, st in enumerate(structs[:max_structs]):
        k = PIC.skey(st)
        if k in out:
            continue
        t0 = time.time()
        cand = session.run(dict(st, min_net_gradient=ng_low, **cam))
        cache, dcache = {}, {}
        rE = ascent_offset(cand, gt_dev, dev_idx, grids, 'efficiency', cache, dcache, PIC.MAX_ROUNDS)
        rJ = ascent_offset(cand, gt_dev, dev_idx, grids, 'jaccard', cache, dcache, PIC.MAX_ROUNDS)
        out[k] = dict(struct=dict(st), n_cand=int(sum(len(d) for d in cand)), n_eval=len(cache), E=rE, J=rJ)
        out[k]['E'].pop('first'), out[k]['J'].pop('first')
        log('  [%d/%d] %s: %d candidates, %d configurations scored, dev E %.1f (rounds %d), dev J %.1f (rounds %d)  [%.0f s, total %.0f s]' % (
            i + 1, len(structs), k, out[k]['n_cand'], len(cache), rE['metrics']['efficiency'], rE['rounds'], 100 * rJ['metrics']['jaccard'], rJ['rounds'], time.time() - t0, time.time() - t_all))
        if cache_file:
            os.makedirs(os.path.dirname(os.path.abspath(cache_file)), exist_ok=True)
            json.dump(dict(signature=sig, results=out), open(cache_file, 'w'), default=C._json_default)
    return out


def choose_offset(results):
    """The structural settings of the best offset-corrected development E (objective E) and J (objective J); the earlier setting on ties (strictly greater replaces)."""
    bestE = bestJ = None
    for k, r in results.items():
        if bestE is None or r['E']['metrics']['efficiency'] > results[bestE]['E']['metrics']['efficiency']:
            bestE = k
        if bestJ is None or r['J']['metrics']['jaccard'] > results[bestJ]['J']['metrics']['jaccard']:
            bestJ = k
    return bestE, bestJ


def plain_metrics(dets, gt, idx):
    """The metrics of ``smlm_common.tune_postfilters`` (no offset correction): with this evaluator ``ascent_offset`` is that function."""
    r = C.score(dets, gt, idx)
    return {k: r[k] for k in C._TUNE_KEYS}


def same_detections(dd, ref, float32=False):
    """(equal to rounding, bit-identical) of two lists of per-frame detection arrays on their first four columns (``float32``: the stored ones were cast to float32)."""
    if len(dd) != len(ref):
        return False, False
    close = bits = True
    dt = np.float32 if float32 else np.float64
    for x, y in zip(dd, ref):
        x, y = np.asarray(x), np.asarray(y)
        if len(x) != len(y):
            return False, False
        if len(x):
            close &= bool(np.allclose(x[:, :4], y[:, :4], rtol=1e-6, atol=1e-3))
            bits &= bool(np.array_equal(x[:, :4].astype(dt), y[:, :4].astype(dt)))
    return close, bits


def part_ps(name, py, cache_dir=None, max_structs=None, log=print):
    t_all = time.time()
    P = pic_setup(name)
    (d0, d1), (t0, t1) = P['d'], P['t']
    dev_idx = np.arange(d0, d1)
    gt = P['gt']
    gt_dev = gt[(gt['frame'] >= d0) & (gt['frame'] < d1)]
    structs = PIC.struct_list()
    sig = hashlib.sha256(json.dumps(dict(name=name, dev=P['d'], struct=[(k, list(v)) for k, v in PIC.STRUCT_GRID.items()], filt=[(k, v) for k, v in P['grids'].items()], ng=P['ng_low'], cam=P['cam'],
                                         rounds=PIC.MAX_ROUNDS, defaults=C.PICASSO_DEFAULTS, objective='offset-corrected'), sort_keys=True, default=str).encode()).hexdigest()[:16]
    cache_file = os.path.join(cache_dir, '%s.ps.json' % name) if cache_dir else None
    A = {}
    with C.PicassoSession(P['dev'], py) as sd:
        # ---- validation: the stored R2E / R2J configurations re-run on the development frames reproduce the stored development detections
        val, stored_corr = {}, {}
        for r in ('R2E', 'R2J'):
            dd = sd.run(dict(P['cam'], **P['stored'][r]))
            close, bits = same_detections(dd, P['stored_dev'][r], float32=fam(name) == '13')
            val[r] = dict(identical_to_rounding=close, bit_identical=bits, detections=int(sum(len(x) for x in dd)))
            stored_corr[r] = offset_corrected_metrics(dd, gt_dev, dev_idx)       # the corrected development value of the stored configuration
            log('PS %s validation %s: the stored configuration re-run on the development frames reproduces the stored detections: %s (bit-identical %s)' % (name, r, close, bits))
        # ---- the search
        results = search_offset(sd, gt_dev, dev_idx, P['cam'], P['ng_low'], P['grids'], structs, cache_file, sig, max_structs, log)
        keyE, keyJ = choose_offset(results)
        chosen = {'R2Eo': (keyE, 'E'), 'R2Jo': (keyJ, 'J')}
        params_of_row, dev_dets, info = {}, {}, {}
        for r, (k, o) in chosen.items():
            params = dict(P['cam'], **PIC.params_of(results[k]['struct'], results[k][o]['filters']))
            if params.get('min_net_gradient') is not None:
                params['min_net_gradient'] = float(params['min_net_gradient'])
            dd = sd.run(params)
            cand = sd.run(dict(results[k]['struct'], min_net_gradient=P['ng_low'], **P['cam']))
            for x, y in zip(dd, C.filter_detections(cand, params)):                      # the self-check of the original search: the filtered low-threshold run equals the direct run
                y = y[:, :4] if len(y) else np.zeros((0, 4))
                x = x[:, :4] if len(x) else np.zeros((0, 4))
                if x.shape != y.shape or not np.allclose(x, y, rtol=1e-6, atol=1e-3, equal_nan=True):
                    raise RuntimeError('%s %s: the filtered low-threshold run differs from the direct run' % (name, r))
            r0 = C.score(dd, gt_dev, dev_idx)
            off = np.nan_to_num(np.array([r0['bias_x_nm'], r0['bias_y_nm']]))
            params_of_row[r], dev_dets[r] = params, dd
            A['%s|config' % r], A['%s|offset' % r], A['%s|key' % r] = jdump(C._merge_params(params)), off, np.array(k)
            info[r] = dict(key=k, objective=o, offset=[float(v) for v in off])
    with C.PicassoSession(P['test'], py) as st_:                                          # the test frames, once per chosen configuration
        for r in chosen:
            dt = st_.run(params_of_row[r])
            off = A['%s|offset' % r]
            for split, ds in (('dev', C.shift_detections(dev_dets[r], *off)), ('test', C.shift_detections(dt, *off))):
                A['%s|%s|off' % (r, split)], A['%s|%s|flat' % (r, split)] = C.pack_detections(ds, columns=(0, 1, 2, 3))
    A['search'] = jdump(results)
    A['stored_config'] = jdump(P['stored'])
    A['stored_corrected'] = jdump(stored_corr)
    A['validation'] = jdump(val)
    A['summary'] = jdump(dict(info=info, structs=len(results), n_eval=int(sum(r['n_eval'] for r in results.values())), n_cand=int(sum(r['n_cand'] for r in results.values())),
                              secs=float(time.time() - t_all), max_structs=max_structs))
    save_part(RAW, 'ps|' + name, A)
    log('PS %s: R2Eo %s dev E %.1f, R2Jo %s dev J %.1f  [%.0f s]' % (name, keyE, results[keyE]['E']['metrics']['efficiency'], keyJ, 100 * results[keyJ]['J']['metrics']['jaccard'], time.time() - t_all))


# ==================================================================================================================================
# the report (from the raw alone)
# ==================================================================================================================================
fmt_ci, table = ISS.fmt_ci, ISS.table


def base_dets(raw, name, key):
    """The detections stored under ``key`` (``key|off``, ``key|flat``) of a set as float64 arrays (the 2013 pipeline stores float32; the dual rows were widened the same way before the shift)."""
    return [np.asarray(d, dtype=np.float64) for d in unpack(raw, key)]


def post_of(name):
    return ISS.post if fam(name) == '16' else E13.post


def dets_of(raw, name, method):
    """The scored test detections of a method of a set, or None when the raw does not hold them. Keys: ISS-*, R0, R2E+c, R2J+c (ref), GR-* (g), GT-* (gb), R2Eo / R2Jo (ps)."""
    ref = 'ref|%s|' % name
    if method in ('ISS-E', 'ISS-J', 'ISS-0', 'R0', 'R2E+c', 'R2J+c'):
        return unpack(raw, ref + method + '|test') if ref + method + '|test|off' in raw else None
    if method in ('R2Eo', 'R2Jo'):
        return unpack(raw, 'ps|%s|%s|test' % (name, method)) if 'ps|%s|%s|test|off' % (name, method) in raw else None
    w = method[-1]
    if method.startswith('GR-'):
        pre, ops_key = 'g|%s|' % name, ref + 'ops'
        if pre + 'frames' not in raw or ops_key not in raw:
            return None
        o = jload(raw[ops_key])[w] if w != '0' else None
        cfg = o['cfg'] if o else (ISS.DEFAULT['cfg'])
        kappa = o['kappa'] if o else ISS.DEFAULT['kappa']
        kind = str(raw[ref + 'kind'])
        key = '%s%s|%s%g|' % (pre, cfg, kind, kappa)
        if key + 'off' not in raw:
            return None
        base = base_dets(raw, name, key[:-1])
        if o is None:
            return base
        return ISS.shift(post_of(name)(base, o['smax']), o['off'])
    if method.startswith('GT-'):
        pre = 'gb|%s|' % name
        if pre + 'ops' not in raw:
            return None
        o = jload(raw[pre + 'ops'])[w]
        kind = str(raw[ref + 'kind'])
        key = '%sT|%s|%s%g' % (pre, o['cfg'], kind, o['kappa'])
        if key + '|off' not in raw:
            return None
        base = base_dets(raw, name, key)
        return ISS.shift(post_of(name)(base, o['smax']), o['off'])
    return None


def dev_dets_of(raw, name, method):
    ref = 'ref|%s|' % name
    if method in ('ISS-E', 'ISS-J', 'R2E+c', 'R2J+c'):
        return unpack(raw, ref + method + '|dev') if ref + method + '|dev|off' in raw else None
    if method in ('R2Eo', 'R2Jo'):
        return unpack(raw, 'ps|%s|%s|dev' % (name, method)) if 'ps|%s|%s|dev|off' % (name, method) in raw else None
    if method in ('GT-E', 'GT-J'):
        pre = 'gb|%s|' % name
        if pre + 'ops' not in raw:
            return None
        o = jload(raw[pre + 'ops'])[method[-1]]
        kind = str(raw[ref + 'kind'])
        key = '%sD|%s|%s%g' % (pre, o['cfg'], kind, o['kappa'])
        if key + '|off' not in raw:
            return None
        base = base_dets(raw, name, key)
        return ISS.shift(post_of(name)(base, o['smax']), o['off'])
    return None


class SetEval:
    """Everything scored on the test frames of one set, lazily and cached: the metrics of a method (primary and challenge setting), its intervals, and the paired differences."""

    def __init__(self, raw, name, boot):
        self.raw, self.name, self.boot = raw, name, boot
        a, b = (int(v) for v in raw['ref|%s|frames' % name][2:4])
        self.idx = np.arange(a, b)
        self.gt, self.pmin = raw['ref|%s|gt_test' % name], float(raw['ref|%s|pmin' % name])
        self._res, self._diff, self._ci, self._has = {}, {}, {}, {}

    def has(self, m):
        if m not in self._has:
            self._has[m] = dets_of(self.raw, self.name, m) is not None
        return self._has[m]

    def res(self, m):
        if m not in self._res:
            ds = dets_of(self.raw, self.name, m)
            if len(ds) != len(self.idx):
                raise ValueError('%s %s: %d detection lists for %d test frames' % (self.name, m, len(ds), len(self.idx)))
            self._res[m] = (score_set(self.name, ds, self.gt, self.idx), score_set(self.name, ds, self.gt, self.idx, self.pmin))
        return self._res[m]

    def ci(self, m):
        if m not in self._ci:
            self._ci[m] = C.bootstrap_ci(self.res(m)[0], B=self.boot, seed=0, block=BLOCK)
        return self._ci[m]

    def diff(self, A, B):
        """Paired differences A - B: values and intervals of J, RMSE, E (primary) and E (challenge)."""
        if (A, B) not in self._diff:
            (ra, ra2), (rb, rb2) = self.res(A), self.res(B)
            d = C.bootstrap_diff_ci(ra['per_frame'], rb['per_frame'], B=self.boot, seed=0, block=BLOCK)
            d2 = C.bootstrap_diff_ci(ra2['per_frame'], rb2['per_frame'], B=self.boot, seed=0, block=BLOCK)
            self._diff[(A, B)] = dict(J=(ra['jaccard'] - rb['jaccard'], d['jaccard']), R=(ra['rmse_nm'] - rb['rmse_nm'], d['rmse_nm']), E=(ra['efficiency'] - rb['efficiency'], d['efficiency']),
                                      Ec=(ra2['efficiency'] - rb2['efficiency'], d2['efficiency']))
        return self._diff[(A, B)]


def dcell(d, scale=1.0):
    return fmt_ci(d[0], d[1], scale)


VERDICT_EPS = 1e-4                                                  # below this (points, nm, or a fraction) a bound is rounding noise: two paths that fit the same kernels in another order differ by ~1e-14, float32 storage by ~1e-5


def verdict(d, higher_better=True):
    """'win' / 'loss' / 'tie' of a paired difference ``d`` = (value, (lo, hi)): the interval excludes 0 (by more than ``VERDICT_EPS``) in favour of / against the first method, else a tie."""
    lo, hi = (d[1][0], d[1][1]) if higher_better else (-d[1][1], -d[1][0])
    return 'win' if lo > VERDICT_EPS else ('loss' if hi < -VERDICT_EPS else 'tie')


def metric_cell(ev, m, key, scale=1.0):
    r = ev.res(m)[0]
    return fmt_ci(r[key], ev.ci(m)[key], scale)


def tally(triples, key='E', higher_better=True):
    """(wins, ties, losses) of the paired difference ``key`` over the (SetEval, A, B) triples whose methods exist."""
    w = t = l = 0
    for ev, A, B in triples:
        if ev.has(A) and ev.has(B):
            v = verdict(ev.diff(A, B)[key], higher_better)
            w, t, l = w + (v == 'win'), t + (v == 'tie'), l + (v == 'loss')
    return w, t, l


def sets_in(raw, kinds=('ref',)):
    return [s for s in SETS if all(('%s|%s|ops' % (k, s)) in raw for k in kinds)]


def report(raw, boot=BOOT):
    """The markdown report from the raw alone."""
    L = ['# Controls for the SMLM comparison: the order of entry, the symmetry of the Picasso tuning, the development-versus-test gap', '']
    sets = sets_in(raw)
    if not sets:
        return '\n'.join(L + ['No reference part in the raw: nothing to report.', '']) + '\n'
    meta = jload(raw['meta|sha']) if 'meta|sha' in raw else {}
    ver = jload(raw['meta|versions']) if 'meta|versions' in raw else {}
    L += ['Generated by `experiments/exp_smlm_controls.py` from `controls_raw.npz`; every number below is printed by the script. **Post hoc**: the analyses A1-A3 were declared in the script\'s docstring after the '
          'results of `iss.md`, `smlm2013.md` and `picasso_baseline.md` were known, in answer to three objections (the order of entry may be irrelevant; the ISS and Picasso are tuned on different '
          'objectives; no development-versus-test gap is reported). Nothing in the existing reports, raws or scripts is changed; they were read (digests of the files the references were copied from: %s). '
          'Frame indices are 0-based stack indices, positions in nm; metrics pooled over the test frames with 95 %% intervals of a bootstrap over blocks of %d consecutive frames (%d resamples, seed 0); '
          'differences are paired (the same resampled blocks for both methods). A paired difference is a **win** (for the first method) when its interval excludes 0 in its favour, a **loss** when it '
          'excludes 0 against it, a **tie** otherwise. Python %s, numpy %s.' % (', '.join('%s %s' % (k, v[:12]) for k, v in sorted(meta.items())) or 'not recorded', BLOCK, boot,
                                                                              ver.get('python', '?'), ver.get('numpy', '?')), '']
    evs = {s: SetEval(raw, s, boot) for s in sets}
    L += status_lines(raw, sets)
    L += validation_lines(raw)
    L += greedy_lines(raw, sets, evs)
    L += own_lines(raw, sets, evs)
    L += picasso_lines(raw, sets, evs)
    L += gap_lines(raw, sets, evs, boot)
    L += limits_lines(raw)
    return '\n'.join(L) + '\n'


def status_lines(raw, sets):
    rows = []
    for s in sets:
        a, b = (int(v) for v in raw['ref|%s|frames' % s][2:4])
        rows.append([SHORT[s], '%d-%d (%d frames)' % (a, b - 1, b - a), 'yes' if 'g|%s|frames' % s in raw else 'no', 'yes' if 'gb|%s|ops' % s in raw else 'no', 'yes' if 'ps|%s|summary' % s in raw else 'no'])
    return ['## What is in this report', '',
            table(['set', 'test frames', 'A1a greedy at the dual\'s points', 'A1b greedy tuned on its own', 'A2 Picasso, offset-corrected objective'], rows), '',
            'A part that is absent for a set was not run (time-box: the 2016 sets first, then the 2013 sets); its rows and tallies below cover the sets that have it.', '']


def validation_lines(raw):
    L = ['## Harness validation', '']
    rows = []
    for s in SETS:
        for k in ('dual_original_stop', 'dual_truncated_stop'):
            if 'val|%s|%s' % (s, k) in raw:
                v = jload(raw['val|%s|%s' % (s, k)])
                rows.append([SHORT[s], 'dual path through this harness, ' + k.replace('dual_', '').replace('_', ' '), v['cfg'], v['frames'], '%d of %d' % (v['rules_identical'], v['rules']), 'yes' if v['identical'] else 'NO'])
        if 'val|%s|guard' % s in raw:
            v = jload(raw['val|%s|guard' % s])
            rows.append([SHORT[s], 'greedy path, guarded against unguarded (`imaging_guard`), %d frames' % v['frames'], v['cfg'], v['frames'], 'barred points: %d' % v['barred'], 'yes' if v['identical'] else 'NO'])
        for k in sorted(k_ for k_ in raw if k_.startswith('g|%s|' % s) and k_.endswith('|dual_identical')):
            rows.append([SHORT[s], '2013 dual path through this harness, whole stored part (detections, entries, capped tiles)', k.split('|')[2], int(raw['g|%s|frames' % s][1] - raw['g|%s|frames' % s][0]), 'all keys', 'yes' if int(raw[k]) else 'NO'])
    if not rows:
        return L + ['No validation run in the raw.', '']
    L += ['Before any greedy number is read: the dual path is re-run through the harness (the frozen pipeline, `imaging.sliding_path` called unchanged) and compared array by array with the stored rows. '
          '"truncated stop": the 2016 runs of this script stop the path at the first entrant below significance 4 instead of 3 entries after the first below 3 (the rule `stop_snr` reads only the entrants above kappa >= 4), '
          'which must change no detection.', '',
          table(['set', 'check', 'configuration', 'frames', 'rules / result', 'identical'], rows), '']
    return L


def op_title(w):
    return {'E': 'E-tuned operating point', 'J': 'J-tuned operating point', '0': 'a priori default'}[w]


def greedy_lines(raw, sets, evs):
    L = ['## A1a. Greedy control: the same code, only the order of entry changed', '',
         'Greedy = `order=\'greedy\'` (matching pursuit with sliding: at every step the candidate of the largest matched-filter significance enters). Same PSF, high-pass, weights, halo, resolution limit, width '
         'bound, stopping rule, refinement, width post-filter and offset as the dual ISS, at the dual ISS\'s own development-chosen operating point (not re-tuned). "d" = greedy minus the row named; '
         'a win means greedy is better, a loss worse.', '']
    gsets = [s for s in sets if evs[s].has('GR-E') or evs[s].has('GR-0')]
    if not gsets:
        return L + ['The greedy part is not in the raw.', '']
    for w in OPS:
        g, d, p = 'GR-' + w, 'ISS-' + w, PIC_OF[w]
        rows = []
        for s in gsets:
            ev = evs[s]
            if not (ev.has(g) and ev.has(d)):
                continue
            o = jload(raw['ref|%s|ops' % s])[w] if w != '0' else None
            cfgs = (('%s kappa %g width <= %.3g' % (o['cfg'], o['kappa'], o['smax'])) if o else '%s kappa %g' % (ISS.DEFAULT['cfg'], ISS.DEFAULT['kappa']))
            dg, dp, dd_ = ev.diff(g, d), ev.diff(g, p), ev.diff(d, p)
            rows.append([SHORT[s], cfgs, metric_cell(ev, d, 'efficiency'), metric_cell(ev, g, 'efficiency'), dcell(dg['E']) + ' ' + verdict(dg['E']), dcell(dg['J'], 100.0) + ' ' + verdict(dg['J']),
                         dcell(dg['R']) + ' ' + verdict(dg['R'], False), dcell(dg['Ec']) + ' ' + verdict(dg['Ec']), dcell(dp['E']) + ' ' + verdict(dp['E']), dcell(dd_['E']) + ' ' + verdict(dd_['E'])])
        L += ['### %s' % op_title(w).capitalize(), '',
              'Efficiency E of the dual and of the greedy path (primary setting), then the paired differences greedy minus dual in E, J (points), RMSE (nm, lower is better) and E in the challenge setting, then greedy minus %s '
              'and, for reference, dual minus %s in E.' % (LABEL[p], LABEL[p]), '',
              table(['set', 'operating point', 'E dual', 'E greedy', 'dE greedy - dual', 'dJ % greedy - dual', 'dRMSE nm greedy - dual', 'dE (challenge) greedy - dual', 'dE greedy - ' + p, 'dE dual - ' + p], rows), '']
    # tallies
    L += ['### Record over the sets (wins / ties / losses by the paired interval, primary setting)', '']
    rows = []
    for w in OPS:
        g, d, p = 'GR-' + w, 'ISS-' + w, PIC_OF[w]
        for lab, A, B in (('greedy vs dual', g, d), ('greedy vs ' + p, g, p), ('dual vs ' + p, d, p)):
            cells = [op_title(w), lab]
            for key, hb in (('E', True), ('J', True), ('R', False)):
                wi, ti, lo = tally([(evs[s], A, B) for s in gsets], key, hb)
                cells.append('%d / %d / %d' % (wi, ti, lo))
            rows.append(cells)
    L += [table(['operating point', 'comparison', 'E', 'J', 'RMSE'], rows), '']
    # mean differences
    rows = []
    for w in OPS:
        g, d = 'GR-' + w, 'ISS-' + w
        ss = [s for s in gsets if evs[s].has(g)]
        if not ss:
            continue
        dE = [evs[s].diff(g, d)['E'][0] for s in ss]
        dJ = [100 * evs[s].diff(g, d)['J'][0] for s in ss]
        dR = [evs[s].diff(g, d)['R'][0] for s in ss]
        rows.append([op_title(w), len(ss), '%+.2f' % np.mean(dE), '%+.2f / %+.2f' % (min(dE), max(dE)), '%+.2f' % np.mean(dJ), '%+.2f' % np.mean(dR)])
    L += ['Mean over the sets of greedy minus dual (point estimates):', '', table(['operating point', 'sets', 'mean dE', 'dE min / max', 'mean dJ (points)', 'mean dRMSE (nm)'], rows), '']
    # how often the paths differ, and path statistics
    rows = []
    for s in gsets:
        for w in OPS:
            g, d = dets_of(raw, s, 'GR-' + w), dets_of(raw, s, 'ISS-' + w)
            if g is None or d is None:
                continue
            ndiff = sum(1 for x, y in zip(g, d) if x.shape != y.shape or (len(x) and not np.allclose(np.asarray(x)[:, :2], np.asarray(y)[:, :2], rtol=0, atol=1e-6)))
            nc = sum(1 for x, y in zip(g, d) if len(x) != len(y))
            rows.append([SHORT[s], op_title(w), len(g), ndiff, nc, sum(len(x) for x in g), sum(len(y) for y in d)])
    L += ['### How often the two paths give different detections', '',
          'Frames whose scored detections differ (a different count or a position differing by more than 1e-6 nm) between the greedy and the dual path at the same operating point, and those whose count differs.', '',
          table(['set', 'operating point', 'test frames', 'frames that differ', 'frames with a different count', 'detections greedy', 'detections dual'], rows), '']
    L += path_stat_lines(raw, gsets)
    L += time_lines(raw, gsets)
    return L


def time_lines(raw, gsets):
    rows = []
    for s in gsets:
        for cfg in sorted({k.split('|')[2] for k in raw if k.startswith('g|%s|' % s) and k.endswith('|secs')}):
            g, d, c = raw.get('g|%s|%s|secs' % (s, cfg)), raw.get('g|%s|%s|dual_secs' % (s, cfg)), raw.get('g|%s|%s|cpu' % (s, cfg))
            rows.append([SHORT[s], cfg, '%.3f' % g.mean(), '%.3f' % c.mean() if c is not None else '', '%.3f' % d.mean() if d is not None else '', len(g)])
    if not rows:
        return []
    return ['### Seconds per frame', '',
            'Wall-clock seconds per frame of the greedy runs and, for the 2016 sets, of the dual runs of the same configuration (both with the shortened stop, one core of a shared machine, so the figures are '
            'inflated and noisy); process CPU seconds for the 2013 runs. The 2016 pipeline stops its path at the first entrant below 4, which is why these are shorter than the times in `iss.md`.', '',
            table(['set', 'configuration', 'greedy wall s / frame', 'greedy CPU s / frame', 'dual wall s / frame', 'test frames'], rows), '']


def path_stat_lines(raw, gsets):
    rows = []
    for s in gsets:
        cfgs = sorted({k.split('|')[2] for k in raw if k.startswith('g|%s|' % s) and k.endswith('|stats')})
        for cfg in cfgs:
            for lab, key in (('greedy', 'g|%s|%s|stats' % (s, cfg)), ('dual', 'g|%s|%s|dual_stats' % (s, cfg))):
                if key in raw:
                    v = raw[key]
                    rows.append([SHORT[s], cfg, lab] + [int(x) for x in v])
    if not rows:
        return []
    return ['### Guard activations and cycles', '',
            '`imaging_guard.sliding_path_guarded` bars a point that exits twice at the same entry time (the enter / exit limit cycle). It supports `order=\'greedy\'` but cannot act there (the entry time is the step count, '
            'so no two exits share it; confirmed by the validation above), so the stored dual results and these greedy results are both from the unguarded path. Counts over all paths of the test frames '
            '(tiles for the 2013 sets) of each configuration run: guard events = exits of a point that had exited at the same entry time before (what the guard would bar); re-entries = entries of a point that '
            'had exited; paths at the cap = paths that ended because the cap on entries was reached.', '',
            table(['set', 'configuration', 'order'] + list(STAT_NAMES), rows), '']


def own_lines(raw, sets, evs):
    gsets = [s for s in sets if evs[s].has('GT-E')]
    L = ['## A1b. Greedy tuned on its own development frames', '']
    if not gsets:
        return L + ['Not run for any set (time-box).', '']
    L += ['The greedy path given the same development search as the dual (four configurations x 7 significance levels x 3 width bounds; for the 2016 sets also the second round around its E-optimal base; offset estimated '
          'and removed on the same frames; `choose_op`), scored once on the test frames. This is "any greedy method with these ingredients and the same tuning effort". "d" = greedy (own tuning) minus the row named.', '']
    rows = []
    for s in gsets:
        ev = evs[s]
        for w in ('E', 'J'):
            g, d, p = 'GT-' + w, 'ISS-' + w, PIC_OF[w]
            og, od = jload(raw['gb|%s|ops' % s])[w], jload(raw['ref|%s|ops' % s])[w]
            dg, dp = ev.diff(g, d), ev.diff(g, p)
            rows.append([SHORT[s], w, '%s kappa %g width <= %.3g' % (od['cfg'], od['kappa'], od['smax']), '%s kappa %g width <= %.3g' % (og['cfg'], og['kappa'], og['smax']), metric_cell(ev, d, 'efficiency'),
                         metric_cell(ev, g, 'efficiency'), dcell(dg['E']) + ' ' + verdict(dg['E']), dcell(dg['J'], 100.0) + ' ' + verdict(dg['J']), dcell(dg['R']) + ' ' + verdict(dg['R'], False),
                         dcell(dp['E']) + ' ' + verdict(dp['E'])])
    L += [table(['set', 'point', 'dual: chosen', 'greedy: chosen', 'E dual', 'E greedy', 'dE greedy - dual', 'dJ % greedy - dual', 'dRMSE nm greedy - dual', 'dE greedy - Picasso R2*+c'], rows), '']
    rows = []
    for w in ('E', 'J'):
        for lab, A, B in (('greedy (own tuning) vs dual', 'GT-', 'ISS-'), ('greedy (own tuning) vs Picasso R2*+c', 'GT-', None)):
            cells = [w + '-tuned', lab]
            for key, hb in (('E', True), ('J', True), ('R', False)):
                pairs = [(evs[s], A + w, (B + w) if B else PIC_OF[w]) for s in gsets]
                cells.append('%d / %d / %d' % tally(pairs, key, hb))
            rows.append(cells)
    L += ['Record over the %d sets (wins / ties / losses):' % len(gsets), '', table(['point', 'comparison', 'E', 'J', 'RMSE'], rows), '']
    return L


def picasso_lines(raw, sets, evs):
    L = ['## A2. Picasso objective symmetry: R2E / R2J re-selected with the ISS\'s objective (E / J after the offset is removed)', '']
    psets = [s for s in sets if evs[s].has('R2Eo')]
    if not psets:
        return L + ['Not run for any set (time-box).', '']
    L += ['Declared search (docstring): the stored ladder does not hold the 108 structural settings x their scored post-filter combinations, so the search was re-run on the development frames only, with the grids and '
          'the coordinate ascent of the original, the objective being E (R2Eo) / J (R2Jo) after the constant offset of that configuration is estimated on the development frames, removed, and the pairing redone, for every candidate '
          'configuration. The chosen configuration was run once on the test frames with the development offset subtracted. R2E+c / R2J+c are the stored rows (selected before the offset, offset appended).', '']
    rows = []
    for s in psets:
        sm = jload(raw['ps|%s|summary' % s])
        res = jload(raw['ps|%s|search' % s])
        sc = jload(raw['ps|%s|stored_corrected' % s])
        stc = jload(raw['ps|%s|stored_config' % s])
        val = jload(raw['ps|%s|validation' % s])
        for r, o, st, ob in (('R2Eo', 'E', 'R2E', 'efficiency'), ('R2Jo', 'J', 'R2J', 'jaccard')):
            k = sm['info'][r]['key']
            m = res[k][o]['metrics']
            same = PIC.describe(jload(raw['ps|%s|%s|config' % (s, r)])) == PIC.describe(stc[st])
            f = (lambda v: v * 100.0) if ob == 'jaccard' else (lambda v: v)
            rows.append([SHORT[s], r, 'yes' if val[st]['bit_identical'] else ('to rounding' if val[st]['identical_to_rounding'] else 'NO'), PIC.describe(stc[st]),
                         '%.1f' % f(sc[st][ob]), PIC.describe(jload(raw['ps|%s|%s|config' % (s, r)])), '%.1f' % f(m[ob]), '(%.1f, %.1f)' % tuple(sm['info'][r]['offset']), 'same' if same else 'different'])
    L += [table(['set', 'row', 'stored config re-run reproduces stored dev detections', 'stored configuration (R2E / R2J)', 'its dev value, offset-corrected', 'offset-corrected selection (R2Eo / R2Jo)',
                 'its dev value (offset-corrected)', 'offset (dx, dy) nm', 'configuration'], rows), '',
          'The development value is E for the E rows and J (in percent) for the J rows, after the offset is removed. The ascent is a local search over the same space, so the new value need not exceed the stored configuration\'s '
          'corrected value, but it normally does. Search sizes: ' + '; '.join('%s %d settings, %d configurations scored' % (SHORT[s], jload(raw['ps|%s|summary' % s])['structs'], jload(raw['ps|%s|summary' % s])['n_eval']) for s in psets) + '.', '']
    # test scores, the stored and the symmetric rows
    rows = []
    for s in psets:
        ev = evs[s]
        for lab in ('R2E+c', 'R2Eo', 'R2J+c', 'R2Jo'):
            r = ev.res(lab)[0]
            rows.append([SHORT[s], LABEL[lab], metric_cell(ev, lab, 'jaccard', 100.0), metric_cell(ev, lab, 'rmse_nm'), metric_cell(ev, lab, 'efficiency'), '%.3f' % (r['n_det'] / len(ev.idx))])
    L += ['### Test frames: the stored and the symmetric Picasso rows', '', table(['set', 'row', 'J %', 'RMSE nm', 'efficiency', 'detected / frame'], rows), '']
    # ISS minus Picasso, stored and symmetric
    L += ['### Dual ISS minus Picasso, with the stored selection and with the symmetric one', '',
          'Paired differences (dual ISS minus Picasso) on the test frames. "stored" = R2E+c / R2J+c, "symmetric" = R2Eo / R2Jo. The verdict is that of the interval of the difference in the metric of the column.', '']
    rows = []
    for s in psets:
        ev = evs[s]
        for w in ('E', 'J'):
            d = 'ISS-' + w
            a, b = ev.diff(d, PIC_OF[w]), ev.diff(d, PICO_OF[w])
            rows.append([SHORT[s], w + '-tuned', dcell(a['E']) + ' ' + verdict(a['E']), dcell(b['E']) + ' ' + verdict(b['E']), dcell(a['J'], 100.0) + ' ' + verdict(a['J']), dcell(b['J'], 100.0) + ' ' + verdict(b['J']),
                         dcell(a['R']) + ' ' + verdict(a['R'], False), dcell(b['R']) + ' ' + verdict(b['R'], False), dcell(a['Ec']) + ' ' + verdict(a['Ec']), dcell(b['Ec']) + ' ' + verdict(b['Ec'])])
    L += [table(['set', 'point', 'dE stored', 'dE symmetric', 'dJ % stored', 'dJ % symmetric', 'dRMSE nm stored', 'dRMSE nm symmetric', 'dE (challenge) stored', 'dE (challenge) symmetric'], rows), '']
    rows = []
    for w in ('E', 'J'):
        for key, lab, hb in (('E', 'E', True), ('J', 'J', True), ('R', 'RMSE', False), ('Ec', 'E (challenge)', True)):
            a = tally([(evs[s], 'ISS-' + w, PIC_OF[w]) for s in psets], key, hb)
            b = tally([(evs[s], 'ISS-' + w, PICO_OF[w]) for s in psets], key, hb)
            rows.append([w + '-tuned ISS', lab, '%d / %d / %d' % a, '%d / %d / %d' % b])
    L += ['Record over the %d sets (wins / ties / losses of the dual ISS):' % len(psets), '', table(['point', 'metric', 'vs the stored Picasso row', 'vs the symmetric Picasso row'], rows), '']
    mean_rows = []
    for w in ('E', 'J'):
        a = [evs[s].diff('ISS-' + w, PIC_OF[w])['E'][0] for s in psets]
        b = [evs[s].diff('ISS-' + w, PICO_OF[w])['E'][0] for s in psets]
        pe = [evs[s].diff(PICO_OF[w], PIC_OF[w])['E'][0] for s in psets]
        mean_rows.append([w + '-tuned', '%+.2f' % np.mean(a), '%+.2f' % np.mean(b), '%+.2f' % np.mean(pe), '%+.2f / %+.2f' % (min(pe), max(pe))])
    L += ['Mean over the sets (point estimates): ', '', table(['point', 'mean dE dual - stored Picasso', 'mean dE dual - symmetric Picasso', 'mean dE of the Picasso row (symmetric - stored)', 'its min / max'], mean_rows), '']
    # Picasso symmetric minus stored, paired
    rows = []
    for s in psets:
        ev = evs[s]
        for w in ('E', 'J'):
            dd = ev.diff(PICO_OF[w], PIC_OF[w])
            rows.append([SHORT[s], w, dcell(dd['E']) + ' ' + verdict(dd['E']), dcell(dd['J'], 100.0) + ' ' + verdict(dd['J']), dcell(dd['R']) + ' ' + verdict(dd['R'], False)])
    L += ['Symmetric minus stored Picasso row (paired, test frames):', '', table(['set', 'point', 'dE', 'dJ %', 'dRMSE nm'], rows), '']
    L += verdict_lines(psets, evs)
    return L


def verdict_lines(psets, evs):
    """The computed answer to objection (b): did the verdicts change?"""
    L = ['### Did the conclusions change?', '']
    changed = []
    for s in psets:
        for w in ('E', 'J'):
            for key, hb in (('E', True), ('J', True), ('R', False)):
                a = verdict(evs[s].diff('ISS-' + w, PIC_OF[w])[key], hb)
                b = verdict(evs[s].diff('ISS-' + w, PICO_OF[w])[key], hb)
                if a != b:
                    changed.append('%s %s-tuned %s: %s -> %s' % (SHORT[s], w, {'E': 'E', 'J': 'J', 'R': 'RMSE'}[key], a, b))
    n = len(psets) * 2 * 3
    L.append('Of the %d verdicts (dual ISS against Picasso: %d sets x E-tuned / J-tuned x E, J, RMSE), %d change when the stored Picasso row is replaced by the symmetric one%s' % (n, len(psets), len(changed), (': ' + '; '.join(changed) + '.') if changed else '.'))
    L.append('')
    return L


def boot_reps(per_frame, B, seed, keys=('jaccard', 'efficiency', 'rmse_nm')):
    cols = C._sum_columns(per_frame)
    reps = {k: np.empty(B) for k in keys}
    for b0, idx in C._boot_indices(len(per_frame), B, seed, BLOCK):
        for j, m in enumerate(C._replicate_metrics(cols, idx, C.ALPHA_PER_NM)):
            for k in keys:
                reps[k][b0 + j] = m[k]
    return reps


def gap_lines(raw, sets, evs, boot):
    L = ['## A3. Development-versus-test gap', '',
         'E and J of each tuned row on the frames it was tuned on (development) and on the held-out test frames (primary setting, all emitters), with the difference test minus development and a 95 %% interval of that difference '
         'from two independent block bootstraps (blocks of %d frames, %d resamples; seeds %d and %d). A negative difference is a loss on new frames, a positive one is no over-fitting (the two ranges are different frames of '
         'different difficulty: the development ranges of the HD sets are 30-100 frames). ISS rows: the development runs with the selected kappa, width bound and offset; Picasso rows: the detections of the selected configuration '
         'on the development frames with the development offset subtracted.' % (BLOCK, boot, *GAP_SEEDS), '']
    methods = ['ISS-E', 'ISS-J', 'R2E+c', 'R2J+c', 'R2Eo', 'R2Jo', 'GT-E', 'GT-J']
    allrows, summ = [], {}
    for s in sets:
        ev = evs[s]
        a, b, c0, c1 = (int(v) for v in raw['ref|%s|frames' % s])
        gt_dev = raw['ref|%s|gt_dev' % s]
        for m in methods:
            dd = dev_dets_of(raw, s, m)
            if dd is None or not ev.has(m):
                continue
            rd = score_set(s, dd, gt_dev, np.arange(a, b))
            rt = ev.res(m)[0]
            xd, xt = boot_reps(rd['per_frame'], boot, GAP_SEEDS[0]), boot_reps(rt['per_frame'], boot, GAP_SEEDS[1])
            cells = [SHORT[s], LABEL[m]]
            for key, scale, nd in (('efficiency', 1.0, 1), ('jaccard', 100.0, 1)):
                diff = xt[key] - xd[key]
                diff = diff[~np.isnan(diff)]
                lo, hi = (np.percentile(diff, 2.5), np.percentile(diff, 97.5)) if len(diff) else (np.nan, np.nan)
                cells += ['%.*f' % (nd, scale * rd[key]), '%.*f' % (nd, scale * rt[key]), '%+.*f [%+.*f, %+.*f]' % (nd, scale * (rt[key] - rd[key]), nd, scale * lo, nd, scale * hi)]
                summ.setdefault((m, key), []).append(scale * (rt[key] - rd[key]))
            cells.append('%d / %d' % (len(range(a, b)), len(ev.idx)))
            allrows.append(cells)
    if not allrows:
        return L + ['No development detections in the raw.', '']
    L += [table(['set', 'row', 'dev E', 'test E', 'test - dev E', 'dev J %', 'test J %', 'test - dev J (points)', 'dev / test frames'], allrows), '']
    rows = []
    for m in methods:
        if (m, 'efficiency') in summ:
            e, j = summ[(m, 'efficiency')], summ[(m, 'jaccard')]
            rows.append([LABEL[m], len(e), '%+.2f' % np.mean(e), '%+.2f / %+.2f' % (min(e), max(e)), '%+.2f' % np.mean(j), '%+.2f / %+.2f' % (min(j), max(j))])
    L += ['Mean of the test-minus-development difference over the sets (point estimates):', '', table(['row', 'sets', 'mean E', 'E min / max', 'mean J (points)', 'J min / max'], rows), '']
    return L


def limits_lines(raw):
    return ['## What this does not establish', '',
            '* **Same operating point is not the best greedy.** A1a keeps the dual\'s development-chosen kappa, width bound and offset, which were tuned for the dual path; A1b gives the greedy path its own tuning where it was run.',
            '* **Frames only.** Intervals cover the sampling of frames (blocks of 10), not other datasets or the simulation; the development and test frames share the structures. The 2016 sets share one structure, the 2013 '
            'sets two, and within each pair the sets are not independent samples.',
            '* **The Picasso re-selection is a local search** (a coordinate ascent over the original grids from the original start), not the global optimum of the offset-corrected objective, as the original was not for the raw one. '
            'It reproduces the stored search space, not its path: the stored R2E / R2J configurations are listed next to the new ones and their offset-corrected development values are given.',
            '* **The offset is estimated per configuration on the development frames** (mean signed error of the matched pairs), as in both existing pipelines; it is a constant, so it cannot change which emitters are '
            'detected except through the re-pairing at the 250 nm tolerance.',
            '* **Post hoc.** The analyses were chosen after the existing results were known; they were declared in the docstring before they were run, and the existing reports are unchanged.', '']


# ==================================================================================================================================
# command line
# ==================================================================================================================================
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('--run', action='store_true', help='compute the parts, merge them into the raw file, write the report')
    ap.add_argument('--from-raw', action='store_true', help='regenerate the report from the raw file alone')
    ap.add_argument('--parts', nargs='+', default=list(PARTS), choices=list(PARTS))
    ap.add_argument('--sets', nargs='+', default=list(SETS), choices=list(SETS))
    ap.add_argument('--frames', type=int, default=None, help='limit the test frames of REF / G / GB to the first N (smoke runs; use a scratch CONTROLS_OUT)')
    ap.add_argument('--force', action='store_true', help='recompute parts that are in the raw')
    ap.add_argument('--max-structs', type=int, default=None, help='PS: limit the structural settings (smoke runs)')
    ap.add_argument('--cache-dir', default=None, help='PS: keep / reuse the search per set here (resumable)')
    ap.add_argument('--no-validate', action='store_true')
    ap.add_argument('--boot', type=int, default=BOOT)
    ap.add_argument('--picasso-python', default=os.environ.get('SMLM_PICASSO_PYTHON'))
    A = ap.parse_args(argv)
    if A.run:
        for part in PARTS:
            if part not in A.parts:
                continue
            for name in A.sets:
                have = load_raw()
                if part == 'REF' and (A.force or 'ref|%s|ops' % name not in have):
                    part_ref(name, A.frames)
                elif part == 'G' and (A.force or 'g|%s|frames' % name not in have):
                    (part_g16 if fam(name) == '16' else part_g13)(name, A.frames, validate=not A.no_validate)
                elif part == 'GB' and (A.force or 'gb|%s|ops' % name not in have):
                    (part_gb16 if fam(name) == '16' else part_gb13)(name, A.frames)
                elif part == 'PS' and (A.force or 'ps|%s|summary' % name not in have):
                    part_ps(name, A.picasso_python, A.cache_dir, A.max_structs)
        if 'meta|sha' not in load_raw():
            part_meta()
        A.from_raw = True
    if A.from_raw:
        raw = load_raw()
        with open(MD, 'w') as f:
            f.write(report(raw, A.boot))
        print('wrote', MD)
    if not (A.run or A.from_raw):
        ap.error('give --run or --from-raw')


if __name__ == '__main__':
    main()
