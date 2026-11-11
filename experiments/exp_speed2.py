"""
The cost of the sliding inverse scale space fit: the original evaluation against `src/core/imaging_fast.py` (the same algorithm, a cheaper arithmetic, an exact early exit), against Picasso
=====================================================================================================================================================================================

Declared before it was run in `results/smlm/speed2_plan.md` (read it first). Development frames and development patches only. Parts (each saved in `results/smlm/speed2_raw.npz` under its own prefix):

    E   equivalence, SMLM: the a priori configuration, the original `exp_smlm_iss.run_frame` and `smlm_fast.run_frame_fast` on 150 (LD) / 40 (HD) development frames per set, detections packed,
        counts, largest differences, and the scored metrics of both (per-frame statistics, for the paired interval)
    V   the variants battery (other configuration, greedy order, no width sliding, two rules): 12 frames per set and variant
    T   cost per frame, SMLM: 24 development frames per set, original / new / new without the early exit, 3 interleaved repeats, CPU and wall seconds
    TT  cost per frame of the E-tuned configuration of each set (the operating point the main results use, chosen on the development frames by `exp_smlm_iss.choose_op`): 12 frames per set, the same engines
    ENV the machine: CPU model, cores, and the load average sampled every 30 s (a log written by a sampler running beside the timing parts, SPEED2_LOADLOG)
    P   Picasso (rung R2E of the ladder) on the same 24 frames, one frame at a time, 3 repeats (needs SMLM_PICASSO_PYTHON, an interpreter that can import picasso)
    S   scaling on simulated frames: S1 the number of emitters at 64 x 64, S2 the frame side at constant density
    A   astronomy development patches: paths truncated at 25 ... 400 entries, original (guarded) against new, equivalence and CPU seconds; AF the complete path of one patch

    python experiments/exp_speed2.py --run E,V,T,P,S,A [--sets ...] [--quick]          (merged into results/smlm/speed2_raw.npz)
    python experiments/exp_speed2.py --run AF                                              (the complete astronomy path, long)
    python experiments/exp_speed2.py --from-raw                                            (regenerates results/smlm/speed2.md byte for byte)

Environment: SPEED2_OUT (directory of the raw and the md, default results/smlm), SMLM_PICASSO_PYTHON, SMLM_DATA, ASTRO_DATA. One process, one thread (OMP_NUM_THREADS = 1), run under `nice`.
"""
import argparse
import json
import os
import subprocess
import sys
import tempfile
import time

os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import numpy as np
from scipy.optimize import linear_sum_assignment

HERE = os.path.abspath(os.path.dirname(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)
import exp_smlm_iss as ISS                  # noqa: E402
import smlm_common as C                     # noqa: E402
import smlm_fast as SF                      # noqa: E402
from src.core import imaging as IM          # noqa: E402
from src.core import imaging_fast as F      # noqa: E402
from src.core import imaging_guard as G     # noqa: E402

OUT = os.environ.get('SPEED2_OUT', os.path.join(REPO, 'results', 'smlm'))
RAW, MD = os.path.join(OUT, 'speed2_raw.npz'), os.path.join(OUT, 'speed2.md')
SETS = ISS.SETS
N_EQ = {'LD': 150, 'HD': 40}                 # equivalence frames per set
N_VAR = 12                                   # variants battery frames per set and variant
N_TIME = 24                                  # timed frames per set
N_TUNED = 12                                 # timed frames per set, E-tuned configuration
REPEATS = 3
BOOT = 1000
# tolerances of the plan (section 3)
TOL_POS_NM, TOL_PH_REL, TOL_SIG_PX, TOL_E = 0.1, 1e-3, 1e-3, 1e-3
TOL_ASTRO_POS, TOL_ASTRO_FLUX, TOL_ASTRO_TAU = 1e-3, 1e-3, 1e-6
VARIANTS = {                                 # name -> (configuration, options overriding the defaults)
    'hp3_k7': ('hp3', dict(rules=[('snr', 7.0)])),
    'greedy': (ISS.DEFAULT['cfg'], dict(rules=[('snr', 5.0)], order='greedy')),
    'nowidth': (ISS.DEFAULT['cfg'], dict(rules=[('snr', 5.0)], widths=None)),
    'two_rules': (ISS.DEFAULT['cfg'], dict(rules=[('snr', 5.0), ('bic', 4)])),
}
ASTRO_PATCHES = ((110, 630), (210, 630))
ASTRO_CAPS = (25, 50, 100, 200, 400)
ASTRO_HP, ASTRO_REX = 3.0, 1.2
S1_K = (1, 2, 4, 8, 12, 16, 24)
S2_SIDES = (32, 48, 64, 96, 128)
S_PHOTONS, S_DENSITY = 1500.0, 6.0           # photons per emitter; emitters per 64 x 64 px at constant density


def evenly_spaced(a, b, n):
    n = min(n, b - a)
    return a + np.round(np.linspace(0, b - a - 1, n)).astype(int)


def compare_dets(da, db):
    """(same count, largest |dx|,|dy| in nm, largest relative photon difference, largest |d sigma| in px) of two detection arrays (k, 4) = x_nm, y_nm, photons, sigma_px, the detections of
    db paired with those of da by the minimum-cost assignment on the distance; (False, nan, nan, nan) for different counts."""
    if len(da) != len(db):
        return False, np.nan, np.nan, np.nan
    if len(da) == 0:
        return True, 0.0, 0.0, 0.0
    cost = np.linalg.norm(da[:, None, :2] - db[None, :, :2], axis=2)
    r, c = linear_sum_assignment(cost)
    a, b = da[r], db[c]
    return True, float(np.abs(a[:, :2] - b[:, :2]).max()), float((np.abs(a[:, 2] - b[:, 2]) / np.maximum(np.abs(a[:, 2]), 1e-9)).max()), float(np.abs(a[:, 3] - b[:, 3]).max())


def summarize_compare(rows):
    """rows: list of compare_dets results. -> array [frames, frames with equal counts, max pos nm, max rel photon, max sigma px] (the maxima over the frames with equal counts)."""
    ok = [r for r in rows if r[0]]
    mx = lambda j: float(max([r[j] for r in ok])) if ok else float('nan')       # noqa: E731
    return np.array([len(rows), len(ok), mx(1), mx(2), mx(3)])


def compare_astro(da, db):
    """(same count, largest position difference px, largest relative flux difference) of two (x, y, flux) triples of an astronomy path."""
    a, b = np.c_[da[0], da[1], da[2], np.zeros(len(da[0]))], np.c_[db[0], db[1], db[2], np.zeros(len(db[0]))]
    if len(a) != len(b):
        return False, np.nan, np.nan
    if len(a) == 0:
        return True, 0.0, 0.0
    cost = np.linalg.norm(a[:, None, :2] - b[None, :, :2], axis=2)
    r, c = linear_sum_assignment(cost)
    return True, float(np.abs(a[r, :2] - b[c, :2]).max()), float((np.abs(a[r, 2] - b[c, 2]) / np.maximum(np.abs(a[r, 2]), 1e-9)).max())


def compare_traces(ta, tb):
    """(events identical, largest relative entry-time difference, barred sets identical) of two paths' traces (events (tau, (b, a), kind), taus, barred)."""
    ea = [(e[1], e[2]) for e in ta['events']]
    eb = [(e[1], e[2]) for e in tb['events']]
    same = ea == eb
    dt = float(np.abs(np.asarray(ta['taus']) - np.asarray(tb['taus'])).max() / max(np.abs(ta['taus']).max(), 1e-300)) if (same and len(ta['taus'])) else float('nan')
    return same, dt, [b for _, b in ta['barred']] == [b for _, b in tb['barred']]


def timed(fn, arg):
    c0, w0 = time.process_time(), time.perf_counter()
    out = fn(arg)
    return out, time.process_time() - c0, time.perf_counter() - w0


def pf_matrix(sc):
    p = sc['per_frame']
    return np.stack([p[k].astype(float) for k in ('frame', 'tp', 'fp', 'fn', 'n_det', 'n_gt', 'sse', 'sx', 'sy')], axis=1)


def pf_struct(m):
    p = np.zeros(len(m), dtype=C.PER_FRAME_DTYPE)
    for j, k in enumerate(('frame', 'tp', 'fp', 'fn', 'n_det', 'n_gt', 'sse', 'sx', 'sy')):
        p[k] = m[:, j]
    return p


def _setup(name, raw_iss, cfg=None, opts=None):
    cal = ISS.get_cal(raw_iss, name)
    o = dict(rules=[('snr', ISS.DEFAULT['kappa'])], tol=1e-8, max_nfev=200)
    o.update(opts or {})
    ISS._init(cfg or ISS.DEFAULT['cfg'], cal['bg'], cal['a'], o)
    return o


def load_frames_of(name, idx):
    stack, gt = ISS.load_set(name)
    return np.stack([C.to_photons(C.load_frames(stack, int(i), int(i) + 1))[0] for i in idx]), gt


# ----------------------------------------------------------------------------------------------------------------------------------
# E, V: equivalence on SMLM frames
# ----------------------------------------------------------------------------------------------------------------------------------
def _run_pair(frames, early=True):
    old, new, t_old, t_new = [], [], 0.0, 0.0
    for f in frames:
        a, ca, _ = timed(ISS.run_frame, (0, f))
        b, cb, _ = timed(lambda x: SF.run_frame_fast(x, early=early), (0, f))
        old.append(a[0])
        new.append(b[0])
        t_old += ca
        t_new += cb
    return old, new, t_old, t_new


def part_E(sets, quick=False, out_path=RAW):
    raw_iss = ISS.load_raw(ISS.RAW)
    for name in sets:
        a, b = ISS.dev_test_frames(name, 'dev')
        n = N_EQ[name[-2:]] // (10 if quick else 1)
        idx = evenly_spaced(a, b, n)
        frames, gt = load_frames_of(name, idx)
        _setup(name, raw_iss)
        old, new, t_old, t_new = _run_pair(frames)
        o0, f0 = ISS.pack([d[0] for d in old])
        o1, f1 = ISS.pack([d[0] for d in new])
        so, sn = C.score([d[0][:, :2] if len(d[0]) else d[0] for d in old], gt, idx), C.score([d[0][:, :2] if len(d[0]) else d[0] for d in new], gt, idx)
        cmp_ = summarize_compare([compare_dets(x[0], y[0]) for x, y in zip(old, new)])
        ISS.save_part(out_path, 'E|' + name, {'frames': idx, 'off_old': o0, 'flat_old': f0, 'off_new': o1, 'flat_new': f1, 'pf_old': pf_matrix(so), 'pf_new': pf_matrix(sn), 'cmp': cmp_,
                                              'cpu': np.array([t_old, t_new])})
        print('E %s: %d frames, counts equal %d, max pos %.3g nm' % (name, len(idx), cmp_[1], cmp_[2]), flush=True)


def part_V(sets, quick=False, out_path=RAW):
    raw_iss = ISS.load_raw(ISS.RAW)
    for name in sets:
        a, b = ISS.dev_test_frames(name, 'dev')
        idx = evenly_spaced(a, b, 3 if quick else N_VAR)
        frames, _ = load_frames_of(name, idx)
        for vname, (cfg, opts) in VARIANTS.items():
            _setup(name, raw_iss, cfg, opts)
            rows = []
            for f in frames:
                x, _, _ = timed(ISS.run_frame, (0, f))
                y, _, _ = timed(SF.run_frame_fast, (0, f))
                for r in range(len(x[0])):
                    rows.append(compare_dets(x[0][r], y[0][r]))
            cmp_ = summarize_compare(rows)
            ISS.save_part(out_path, 'V|%s|%s' % (vname, name), {'frames': idx, 'cmp': cmp_})
            print('V %s %s: %d, counts equal %d, max pos %.3g nm' % (vname, name, cmp_[0], cmp_[1], cmp_[2]), flush=True)


# ----------------------------------------------------------------------------------------------------------------------------------
# T, P: the cost per frame
# ----------------------------------------------------------------------------------------------------------------------------------
def time_frames(frames):
    """Per frame and per engine (original, new, new without the early exit) the minimum over REPEATS interleaved repeats of the CPU and of the wall seconds; returns (cpu, wall) (3, n) and the
    detection counts of the original and of the new (n, 2)."""
    n = len(frames)
    cpu, wall = np.full((3, n), np.inf), np.full((3, n), np.inf)
    nd = np.zeros((n, 2), int)
    fns = [ISS.run_frame, SF.run_frame_fast, lambda x: SF.run_frame_fast(x, early=False)]
    for rep in range(REPEATS):
        for i, f in enumerate(frames):
            for e, fn in enumerate(fns):
                out, c, w = timed(fn, (0, f))
                cpu[e, i], wall[e, i] = min(cpu[e, i], c), min(wall[e, i], w)
                if rep == 0 and e < 2:
                    nd[i, e] = len(out[0][0])
    return cpu, wall, nd


def part_T(sets, quick=False, out_path=RAW):
    raw_iss = ISS.load_raw(ISS.RAW)
    for name in sets:
        a, b = ISS.dev_test_frames(name, 'dev')
        idx = evenly_spaced(a, b, 4 if quick else N_TIME)
        frames, _ = load_frames_of(name, idx)
        _setup(name, raw_iss)
        cpu, wall, nd = time_frames(frames)
        ISS.save_part(out_path, 'T|' + name, {'frames': idx, 'cpu': cpu, 'wall': wall, 'ndet': nd})
        print('T %s: CPU s / frame old %.3f new %.3f new-no-early %.3f' % (name, cpu[0].mean(), cpu[1].mean(), cpu[2].mean()), flush=True)


def part_TT(sets, quick=False, out_path=RAW):
    raw_iss = ISS.load_raw(ISS.RAW)
    for name in sets:
        _, gt = ISS.load_set(name)
        op = ISS.choose_op(ISS.op_table(raw_iss, name, gt), 'E')
        a, b = ISS.dev_test_frames(name, 'dev')
        idx = evenly_spaced(a, b, 3 if quick else N_TUNED)
        frames, _ = load_frames_of(name, idx)
        _setup(name, raw_iss, op['cfg'], dict(rules=[('snr', op['kappa'])]))
        cpu, wall, nd = time_frames(frames)
        ISS.save_part(out_path, 'TT|' + name, {'frames': idx, 'cpu': cpu, 'wall': wall, 'ndet': nd, 'cfg': np.array(op['cfg']), 'kappa': np.array(float(op['kappa']))})
        print('TT %s (%s, kappa %g): CPU s / frame old %.3f new %.3f new-no-early %.3f' % (name, op['cfg'], op['kappa'], cpu[0].mean(), cpu[1].mean(), cpu[2].mean()), flush=True)


def part_ENV(out_path=RAW):
    model = [l.split(':', 1)[1].strip() for l in open('/proc/cpuinfo') if l.startswith('model name')][:1]
    arrays = {'cpu_model': np.array(model[0] if model else 'unknown'), 'cores': np.array(os.cpu_count())}
    log = os.environ.get('SPEED2_LOADLOG')
    if log and os.path.exists(log):
        arrays['load'] = np.array([[float(v) for v in l.split()] for l in open(log) if len(l.split()) == 4])       # unix time, 1, 5, 15 minute load averages
    ISS.save_part(out_path, 'ENV', arrays)


def part_P(sets, quick=False, out_path=RAW):
    py = os.environ.get('SMLM_PICASSO_PYTHON')
    if not py:
        raise SystemExit('set SMLM_PICASSO_PYTHON to an interpreter that can import picasso')
    lad = np.load(os.path.join(REPO, 'results', 'smlm', 'picasso_ladder_raw.npz'), allow_pickle=True)
    for name in sets:
        a, b = ISS.dev_test_frames(name, 'dev')
        idx = evenly_spaced(a, b, 4 if quick else N_TIME)
        stack, _ = ISS.load_set(name)
        frames = np.stack([C.load_frames(stack, int(i), int(i) + 1)[0] for i in idx])
        cfg = str(lad[name + '|R2E|config'])
        with tempfile.TemporaryDirectory(prefix='speed2_') as tmp:
            fin, fp, fout = os.path.join(tmp, 'f.npy'), os.path.join(tmp, 'p.json'), os.path.join(tmp, 'o.npz')
            np.save(fin, frames)
            open(fp, 'w').write(cfg)
            env = dict(os.environ, OMP_NUM_THREADS='1', NUMBA_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
            r = subprocess.run([py, '-I', os.path.join(HERE, 'speed2_picasso_worker.py'), fin, fp, fout, str(REPEATS)], env=env, capture_output=True, text=True)
            if r.returncode != 0:
                raise RuntimeError('Picasso worker failed:\n' + (r.stderr or r.stdout)[-3000:])
            z = np.load(fout)
            ISS.save_part(out_path, 'P|' + name, {'frames': idx, 'cpu': z['cpu'].min(0), 'wall': z['wall'].min(0), 'nloc': z['nloc'], 'config': np.array(cfg), 'version': np.array(str(z['picasso_version']))})
        print('P %s: Picasso CPU s / frame %.4f' % (name, z['cpu'].min(0).mean()), flush=True)


# ----------------------------------------------------------------------------------------------------------------------------------
# S: scaling on simulated frames
# ----------------------------------------------------------------------------------------------------------------------------------
def sim_frame(n, K, seed, bg, a_slope):
    """An n x n EM-CCD frame of K emitters (photons S_PHOTONS, PSF width SIGMA0, uniform positions >= 3 px from the border), the background bg (scalar photons per pixel)."""
    rng = np.random.default_rng(seed)
    cfg = ISS.get_cfg(ISS.DEFAULT['cfg'])
    Pp = IM.PixelPSF(n, n, ISS.SIGMA0, parts=cfg['parts'])
    c = rng.uniform(3.0, n - 4.0, size=(K, 2))
    return IM.simulate_emccd(Pp, c, np.full(K, S_PHOTONS), bg, rng)


def fit_generic(frame, bg, a_slope, engine, early=True, kappa=5.0):
    """The a priori fit (the configuration of `ISS.DEFAULT`, the settings of `ISS.run_frame`) of one frame of any size by the original or the new engine; returns the detections (k, 4)."""
    ny, nx = frame.shape
    cfg = ISS.get_cfg(ISS.DEFAULT['cfg'])
    P = IM.PixelPSF(ny, nx, ISS.SIGMA0, parts=cfg['parts'], hp=cfg['hp'])
    Pp = IM.PixelPSF(ny, nx, ISS.SIGMA0, parts=cfg['parts'])
    new = engine == 'new'
    L = (F.FastLattice if new else IM.Lattice)(P, h=ISS.LAT_H, margin=1.0)
    bgv = np.full((ny, nx), float(bg))
    synth = (F.Kern(Pp).synth if new else Pp.synth)

    def sw(c, w):
        m = synth(c, w) if len(w) else 0.0
        return 1.0 / np.sqrt(a_slope * (bgv + np.maximum(m, 0.0)) + ISS.B0)
    y0 = frame - bgv
    y = P.highpass(y0 - np.median(y0))
    sw0 = sw(np.zeros((0, 3)), np.zeros(0))
    kw = dict(max_entries=ISS.MAX_ENTRIES, sw_fn=sw, sw_on='params', sig_bounds=ISS.WIDTHS, tol=1e-8, max_nfev=200, passes=1, rtrunc=4.0, r_ex=cfg.get('r_ex', ISS.R_EX), snr_min=ISS.SNR_MIN,
              snr_extra=ISS.SNR_EXTRA, extra=6)
    if new:
        R = F.sliding_path_fast(P, y, sw0, L, cycle_guard=False, stop_kappa=kappa if early else None, **kw)
    else:
        R = IM.sliding_path(P, y, sw0, L, **kw)
    i = IM.stop_snr(R['hist'], R['taus'], R['snr'], kappa)[0]
    h = R['hist'][i]
    c, w = h['c'], h['w']
    if len(w):
        c, w, _ = (F.refine_fast if new else IM.refine)(P, y, sw0, c, w, sig_bounds=ISS.WIDTHS, box=1.0, sw_fn=sw, sw_on='params')
    return np.c_[(c[:, 0] + 0.5) * C.PIXEL_NM, (c[:, 1] + 0.5) * C.PIXEL_NM, w, c[:, 2]] if len(w) else np.zeros((0, 4))


def part_S(quick=False, out_path=RAW):
    raw_iss = ISS.load_raw(ISS.RAW)
    cal = ISS.get_cal(raw_iss, 'MT0.N1.LD')
    bg, a_slope = float(np.mean(cal['bg'])), float(cal['a'])
    rows = []                                          # (sweep, size, K, frame, engine, cpu)
    plan = [(1, 64, K, 3) for K in (S1_K[:3] if quick else S1_K)] + [(2, n, max(1, int(round(S_DENSITY * n * n / 4096.0))), 2 if n >= 128 else 3) for n in (S2_SIDES[:3] if quick else S2_SIDES)]
    for sweep, n, K, nfr in plan:
        for fr in range(nfr):
            frame = sim_frame(n, K, 1000 * sweep + 17 * n + 101 * K + fr, bg, a_slope)
            for e, eng in enumerate(('old', 'new')):
                best = np.inf
                for _ in range(2):
                    _, c, _ = timed(lambda f: fit_generic(f, bg, a_slope, eng), frame)
                    best = min(best, c)
                rows.append((sweep, n, K, fr, e, best))
            print('S sweep %d side %d K %d frame %d: old %.3f new %.3f' % (sweep, n, K, fr, rows[-2][5], rows[-1][5]), flush=True)
    ISS.save_part(out_path, 'S', {'rows': np.array(rows, float), 'bg': np.array(bg), 'a': np.array(a_slope)})


# ----------------------------------------------------------------------------------------------------------------------------------
# A: astronomy
# ----------------------------------------------------------------------------------------------------------------------------------
def run_astro(eng, xy, cap_n):
    import astro_fast as AF
    import exp_astro_m2 as X
    patch = X.load_patch(*xy)
    with AF.engine(eng), AF.cap(cap_n):
        c0, w0 = time.process_time(), time.perf_counter()
        out, _, n, info = X.run_path(patch, ASTRO_HP, ASTRO_REX)
        cpu, wall = time.process_time() - c0, time.perf_counter() - w0
        trace = dict(AF.TRACE)
    return out, cpu, wall, n, trace


def _astro_pair(xy, cap_n, repeats, engines=('old', 'new')):
    res = {}
    for eng, name in zip(engines, ('old', 'new')):
        best = np.inf
        for _ in range(repeats):
            out, cpu, wall, n, trace = run_astro(eng, xy, cap_n)
            best = min(best, cpu)
        res[name] = (out, best, wall, n, trace)
    return res


def _astro_arrays(res):
    ev_same, dtau, bar_same = compare_traces(res['old'][4], res['new'][4])
    keys = sorted(res['old'][0])
    cmp_ = [compare_astro(res['old'][0][k], res['new'][0][k]) for k in keys]
    arr = dict(cpu=np.array([res['old'][1], res['new'][1]]), wall=np.array([res['old'][2], res['new'][2]]), n=np.array([res['old'][3], res['new'][3]]),
               trace=np.array([ev_same, dtau, bar_same], float), kappas=np.array(keys),
               cmp=np.array([[float(c[0]), c[1], c[2]] for c in cmp_]), ndet=np.array([[len(res['old'][0][k][0]), len(res['new'][0][k][0])] for k in keys]))
    return arr


def part_A(quick=False, out_path=RAW):
    caps = ASTRO_CAPS[:2] if quick else ASTRO_CAPS
    for xy in ASTRO_PATCHES:
        for cap_n in caps:
            res = _astro_pair(xy, cap_n, 1 if quick else 2)
            ISS.save_part(out_path, 'A|x%d_y%d|%d' % (xy + (cap_n,)), _astro_arrays(res))
            print('A patch %s cap %d: entries %d, CPU old %.1f new %.1f, events identical %s' % (xy, cap_n, res['old'][3], res['old'][1], res['new'][1], res['old'][4]['events'] == res['new'][4]['events'] or
                                                                                              [(e[1], e[2]) for e in res['old'][4]['events']] == [(e[1], e[2]) for e in res['new'][4]['events']]), flush=True)


def part_AC(out_path=RAW):
    """POST HOC control (not in the plan): the original guarded path against itself with the block-fit tolerance 10 times smaller, on the lengths and patches where the new path failed the tolerances."""
    for xy in ASTRO_PATCHES:
        for cap_n in (200, 400):
            res = _astro_pair(xy, cap_n, 1, engines=('old', 'old_tight'))
            ISS.save_part(out_path, 'AC|x%d_y%d|%d' % (xy + (cap_n,)), _astro_arrays(res))
            print('AC patch %s cap %d: entries %d, CPU %.1f / %.1f' % (xy, cap_n, res['old'][3], res['old'][1], res['new'][1]), flush=True)


def part_AX(out_path=RAW):
    """POST HOC equivalence study of the COMPLETE paths (not in the plan; run after the complete path of patch (110, 630) turned out to differ in its events): on both development patches the
    complete path (about 1300 entries) by the original, the new engine and the original at a ten times smaller block-fit tolerance; the F-measure counts per 20 x 20 px block (mag < 22) at every
    threshold and the BIC block, for the paired block bootstrap."""
    import astro_common as A
    import exp_astro_m2 as X
    truth = X.truth_arrays()
    for xy in ASTRO_PATCHES:
        t = A.truth_for_patch(truth, xy[0], xy[1])
        arrays = {}
        for eng in ('old', 'new', 'old_tight'):
            out, cpu, wall, n, trace = run_astro(eng, xy, 2500)
            keys = sorted(out)
            arrays[eng + '|kappas'] = np.array(keys)
            arrays[eng + '|blk'] = np.stack([X.score_dets(out[k], t)['blk22'] for k in keys])
            arrays[eng + '|cpu'] = np.array(cpu)
            arrays[eng + '|n'] = np.array(n)
            print('AX patch %s %s: %d entries, %.0f CPU s' % (xy, eng, n, cpu), flush=True)
        ISS.save_part(out_path, 'AX|x%d_y%d' % xy, arrays)


def part_AF(out_path=RAW):
    xy = ASTRO_PATCHES[0]
    res = _astro_pair(xy, 2500, 1)
    arr = _astro_arrays(res)
    ISS.save_part(out_path, 'AF|x%d_y%d' % xy, arr)
    print('AF patch %s: entries %d, CPU old %.1f new %.1f' % (xy, res['old'][3], res['old'][1], res['new'][1]), flush=True)


# ----------------------------------------------------------------------------------------------------------------------------------
# report (from the raw alone)
# ----------------------------------------------------------------------------------------------------------------------------------
def boot_ratio(num, den, B=BOOT, seed=0):
    """(ratio of means, lower, upper): sum(num)/sum(den) over frames, the 95 % percentile interval of the frame bootstrap (resampling frames with replacement, B resamples, seed)."""
    num, den = np.asarray(num, float), np.asarray(den, float)
    n = len(num)
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, n, size=(B, n))
    r = num[idx].sum(1) / den[idx].sum(1)
    return float(num.sum() / den.sum()), float(np.percentile(r, 2.5)), float(np.percentile(r, 97.5))


def fmt_ratio(t):
    return '%.1f [%.1f, %.1f]' % t


def pooled_boot(cpus, pics, B=BOOT, seed=0):
    """The four pooled ratios of the cost table with a 95 % percentile interval: cpus = list over sets of (3, n) CPU seconds (rows: old, new, new without the early exit), pics = list over sets of (n,)
    Picasso seconds on the same frames. Pooled = ratio of the means over sets of the sets' mean seconds (as the pooled line of the table); the frames of each set are resampled with replacement (the
    same frame draw for every engine), B resamples, seed. Returns {'old/new': (point, lo, hi), 'new/pic': ..., 'nonearly/pic': ..., 'old/pic': ...}."""
    rng = np.random.default_rng(seed)
    draws = [rng.integers(0, c.shape[1], size=(B, c.shape[1])) for c in cpus]

    def pooled(idx):
        m = lambda get: sum(get(c, p, None if idx is None else idx[j]) for j, (c, p) in enumerate(zip(cpus, pics))) / len(cpus)      # noqa: E731
        o = m(lambda c, p, i: c[0].mean() if i is None else c[0][i].mean(1))
        n = m(lambda c, p, i: c[1].mean() if i is None else c[1][i].mean(1))
        ne = m(lambda c, p, i: c[2].mean() if i is None else c[2][i].mean(1))
        pc = m(lambda c, p, i: p.mean() if i is None else p[i].mean(1))
        return {'old/new': o / n, 'new/pic': n / pc, 'nonearly/pic': ne / pc, 'old/pic': o / pc}
    pt, bs = pooled(None), pooled(draws)
    return {k: (float(pt[k]), float(np.percentile(bs[k], 2.5)), float(np.percentile(bs[k], 97.5))) for k in pt}


def span(v, fmt='%.0f'):
    """'a-b' of the smallest and the largest value, 'a' when they print alike."""
    a, b = fmt % min(v), fmt % max(v)
    return a if a == b else a + '-' + b


def fit_slope(x, y, B=BOOT, seed=0):
    """Least-squares slope of log y on log x with a 95 % bootstrap interval over the (x, y) points."""
    lx, ly = np.log(np.asarray(x, float)), np.log(np.asarray(y, float))
    s = np.polyfit(lx, ly, 1)[0]
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(B):
        i = rng.integers(0, len(lx), len(lx))
        if np.ptp(lx[i]) == 0:
            continue
        out.append(np.polyfit(lx[i], ly[i], 1)[0])
    return float(s), float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))


def md_table(header, rows):
    return '| ' + ' | '.join(header) + ' |\n|' + '|'.join(['---'] * len(header)) + '|\n' + '\n'.join('| ' + ' | '.join(r) + ' |' for r in rows)


def report(raw):
    L = ['# The cost of the ISS frame fit: the original evaluation, the fast evaluation of the same algorithm, and Picasso', '',
         'Generated by `experiments/exp_speed2.py` from `speed2_raw.npz`; every number below is printed by the script. Declared in `speed2_plan.md` (committed before the runs). Development frames and '
         'development patches only. One process, one thread, a machine shared with other jobs: the seconds are the minimum over interleaved repeats of the PROCESS CPU time (the wall time is in the raw), '
         'so they are the cost with the least interference, not a benchmark of an idle machine. "old" = `exp_smlm_iss.run_frame` (`imaging.sliding_path`, `imaging.refine`, scipy `least_squares`), '
         '"new" = `smlm_fast.run_frame_fast` (`imaging_fast`), "new, no early exit" = the same without the exit at the first entrant below the significance threshold.', '']
    sets = [s for s in SETS if 'T|%s|cpu' % s in raw]
    if 'ENV|cpu_model' in raw:
        env = 'Platform: %s, %d cores' % (str(raw['ENV|cpu_model']), int(raw['ENV|cores']))
        if 'ENV|load' in raw and len(raw['ENV|load']):
            ld = raw['ENV|load'][:, 1]
            env += '; 1-minute load average sampled every 30 s during the timing parts (%d samples): minimum %.1f, median %.1f, maximum %.1f. The machine was shared with other jobs throughout, so absolute seconds are inflated and noisy; every ratio below compares engines measured interleaved in the same session on this one machine (never across machines).' % (len(ld), ld.min(), np.median(ld), ld.max())
        L += [env, '']
    # ---- equivalence
    L += ['## 1. Is the output the same? (SMLM, development frames, the a priori configuration)', '',
          'Tolerances fixed in the plan: equal detection counts on every frame; positions within %.1f nm, photons within a relative %.0e, widths within %.0e px; pooled efficiency E (in the challenge units) within %.0e '
          '(paired block-bootstrap interval of E_new - E_old must contain 0).' % (TOL_POS_NM, TOL_PH_REL, TOL_SIG_PX, TOL_E), '']
    rows = []
    for s in [s for s in SETS if 'E|%s|cmp' % s in raw]:
        cm, po, pn = raw['E|%s|cmp' % s], pf_struct(raw['E|%s|pf_old' % s]), pf_struct(raw['E|%s|pf_new' % s])
        so, sn = C.pooled_metrics(*[float(po[k].sum()) for k in ('tp', 'fp', 'fn', 'sse', 'sx', 'sy')]), C.pooled_metrics(*[float(pn[k].sum()) for k in ('tp', 'fp', 'fn', 'sse', 'sx', 'sy')])
        ci = C.bootstrap_diff_ci(pn, po, B=BOOT, keys=('efficiency', 'jaccard'))
        dE = sn['efficiency'] - so['efficiency']
        ok = cm[1] == cm[0] and cm[2] <= TOL_POS_NM and cm[3] <= TOL_PH_REL and cm[4] <= TOL_SIG_PX and abs(dE) <= TOL_E
        rows.append([s, '%d' % cm[0], '%d / %d' % (cm[1], cm[0]), '%.2e' % cm[2], '%.2e' % cm[3], '%.2e' % cm[4], '%.3f' % so['efficiency'], '%.3f' % sn['efficiency'], '%.2e' % dE,
                     '[%.2e, %.2e]' % (ci['efficiency'][0], ci['efficiency'][1]), 'pass' if ok else 'FAIL'])
    L.append(md_table(['set', 'frames', 'equal counts', 'max position diff (nm)', 'max relative photon diff', 'max width diff (px)', 'E old', 'E new', 'E new - old', 'interval', 'verdict'], rows))
    L.append('')
    vrows = []
    for v in VARIANTS:
        for s in [s for s in SETS if 'V|%s|%s|cmp' % (v, s) in raw]:
            cm = raw['V|%s|%s|cmp' % (v, s)]
            ok = cm[1] == cm[0] and cm[2] <= TOL_POS_NM and cm[3] <= TOL_PH_REL and cm[4] <= TOL_SIG_PX
            vrows.append([v, s, '%d' % cm[0], '%d / %d' % (cm[1], cm[0]), '%.2e' % cm[2], '%.2e' % cm[3], '%.2e' % cm[4], 'pass' if ok else 'FAIL'])
    if vrows:
        L += ['Variants battery (12 development frames per set; "two_rules" keeps the whole tail, the early exit is off; frames x rules counted):', '',
              md_table(['variant', 'set', 'detection lists', 'equal counts', 'max position diff (nm)', 'max relative photon diff', 'max width diff (px)', 'verdict'], vrows), '']
    # ---- cost per frame
    L += ['## 2. Cost per frame (SMLM, a priori configuration, 24 development frames per set)', '']
    rows, tot = [], {'old': 0.0, 'new': 0.0, 'ne': 0.0, 'pic': 0.0}
    per_set = {}                                                      # set -> (old/new, new/pic, no early exit/pic, old/pic, early exit alone, algebra alone) boot_ratio triples
    for s in sets:
        cpu = raw['T|%s|cpu' % s]
        pic = raw['P|%s|cpu' % s] if 'P|%s|cpu' % s in raw else None
        n = cpu.shape[1]
        per_set[s] = (boot_ratio(cpu[0], cpu[1]), boot_ratio(cpu[1], pic) if pic is not None else None, boot_ratio(cpu[2], pic) if pic is not None else None,
                      boot_ratio(cpu[0], pic) if pic is not None else None, boot_ratio(cpu[2], cpu[1]), boot_ratio(cpu[0], cpu[2]))
        rows.append([s, str(n), '%.4f' % cpu[0].mean(), '%.4f' % cpu[2].mean(), '%.4f' % cpu[1].mean(), '%.4f' % pic.mean() if pic is not None else '-', fmt_ratio(boot_ratio(cpu[0], cpu[1])),
                     fmt_ratio(boot_ratio(cpu[2], cpu[1])), fmt_ratio(boot_ratio(cpu[0], cpu[2])), '%.1f' % float(np.median(cpu[0] / cpu[1])),
                     '%.0f [%.0f, %.0f]' % boot_ratio(cpu[1], pic) if pic is not None else '-', '%.0f [%.0f, %.0f]' % boot_ratio(cpu[0], pic) if pic is not None else '-'])
        tot['old'] += cpu[0].mean() / len(sets)
        tot['new'] += cpu[1].mean() / len(sets)
        tot['ne'] += cpu[2].mean() / len(sets)
        tot['pic'] += (pic.mean() / len(sets)) if pic is not None else 0.0
    L.append(md_table(['set', 'frames', 'old CPU s / frame', 'new, no early exit', 'new', 'Picasso R2E', 'old / new', 'early exit alone (no-early / new)', 'algebra alone (old / no-early)', 'median per-frame old / new',
                       'new / Picasso', 'old / Picasso'], rows))
    L += ['', 'Intervals: 95 %% percentile bootstrap over the frames (%d resamples, seed 0) of the ratio of the sums over frames. Picasso: rung R2E of the ladder (`picasso_ladder_raw.npz`, options per set), identification and fit of '
          'one frame at a time after a JIT warm-up, Picasso %s, the same CPU-time minimum over %d repeats.' % (BOOT, str(raw['P|%s|version' % sets[0]]) if 'P|%s|version' % sets[0] in raw else '-', REPEATS), '']
    if tot['pic'] > 0:
        L += ['**Pooled over the four sets (mean of the sets\' means):** old %.4f s / frame, new %.4f s / frame, Picasso %.5f s / frame: old / new = %.1f, new / Picasso = %.0f, old / Picasso = %.0f.' %
              (tot['old'], tot['new'], tot['pic'], tot['old'] / tot['new'], tot['new'] / tot['pic'], tot['old'] / tot['pic']), '']
        pb = pooled_boot([raw['T|%s|cpu' % s] for s in sets], [raw['P|%s|cpu' % s] for s in sets])
        pt = lambda k: '%.0f [%.0f, %.0f]' % pb[k]                    # noqa: E731
        L += ['**Three cost bands against Picasso** (per set: ratio of the sums over the 24 frames; pooled: ratio of the means of the sets\' mean seconds, with a 95 %% percentile interval of a frame bootstrap '
              'within each set, %d resamples, seed 0, the same frame draw for every engine): the new evaluation with its early exit at one fixed threshold is %sx slower than Picasso per set (pooled %s); the new '
              'evaluation without the early exit (the whole path kept) is %sx per set (pooled %s); the original is %sx per set (pooled %s). Only the first band is the cost of a single fixed threshold; the other '
              'two are the cost of the path.' % (BOOT, span([v[1][0] for v in per_set.values()]), pt('new/pic'), span([v[2][0] for v in per_set.values()]), pt('nonearly/pic'),
                                                  span([v[3][0] for v in per_set.values()]), pt('old/pic')), '']
        lo10 = [s for s in sets if per_set[s][0][1] < 10.0]
        L += ['**What the speed-up consists of.** The pooled old / new = %.1f [%.1f, %.1f] mixes two effects: the exact arithmetic alone (old / new without the early exit) is %sx per set, and the rest is the '
              'early exit (%sx per set), which exists only for a single fixed `snr` threshold and is off for any other rule or for a kept path. Per set old / new is %sx. The plan states "at least 10x" for a '
              'set only if the lower end of its interval is >= 10: %s. The pooled lower end is %.1f (%s).' %
              (pb['old/new'] + (span([v[5][0] for v in per_set.values()], '%.1f'), span([v[4][0] for v in per_set.values()], '%.1f'), span([v[0][0] for v in per_set.values()]),
                                ('this does not hold for ' + ', '.join('%s (lower end %.1f)' % (s, per_set[s][0][1]) for s in lo10) + ', which does not reach the plan\'s own 10x rule') if lo10
                                else 'this holds for every set', pb['old/new'][1], 'at least 10' if pb['old/new'][1] >= 10.0 else 'below 10')), '']
        L += ['**The limit-cycle guard never triggered on real data.** The SMLM engines run with `cycle_guard` off (as the originals do), and on the astronomy patches the barred set was empty for both '
              'engines on the complete paths in the referee\'s rerun (the raw stores only whether the two barred lists agree, not whether they are empty, so "barred set identical" in section 4 is vacuous on real data). '
              'The equivalence of the guard is tested on planted cycles only (`tests/test_imaging_fast.py`).', '']
    trows = []
    for s in [s for s in SETS if 'TT|%s|cpu' % s in raw]:
        cpu = raw['TT|%s|cpu' % s]
        pic = raw['P|%s|cpu' % s] if 'P|%s|cpu' % s in raw else None
        trows.append([s, str(raw['TT|%s|cfg' % s]), '%g' % float(raw['TT|%s|kappa' % s]), str(cpu.shape[1]), '%.4f' % cpu[0].mean(), '%.4f' % cpu[2].mean(), '%.4f' % cpu[1].mean(), fmt_ratio(boot_ratio(cpu[0], cpu[1])), '%d / %d' % (int((raw['TT|%s|ndet' % s][:, 0] == raw['TT|%s|ndet' % s][:, 1]).sum()), cpu.shape[1]),
                      '%.0f [%.0f, %.0f]' % boot_ratio(cpu[1], pic.mean() * np.ones(cpu.shape[1])) if pic is not None else '-', '%.0f [%.0f, %.0f]' % boot_ratio(cpu[0], pic.mean() * np.ones(cpu.shape[1])) if pic is not None else '-'])
    if trows:
        L += ['The E-tuned operating point of each set (configuration and significance threshold chosen on the development frames by the main study; the width bound and the offset are post-processing and cost nothing), a single threshold, 12 development frames per set; Picasso as above (mean of its 24 frames):', '',
              md_table(['set', 'configuration', 'kappa', 'frames', 'old CPU s / frame', 'new, no early exit', 'new', 'old / new', 'frames with equal detection counts', 'new / Picasso', 'old / Picasso'], trows), '']
    # ---- scaling
    if 'S|rows' in raw:
        S = raw['S|rows']
        L += ['## 3. Scaling on simulated frames', '', 'Per setting the mean over frames of the CPU seconds (minimum of 2 repeats per frame); slopes are log-log least-squares exponents with a 95 % bootstrap interval over the individual frame timings.', '']
        for sweep, title, col in ((1, 'S1: 64 x 64 frames, the number of emitters K', 2), (2, 'S2: frame side n at constant density (6 emitters per 64 x 64)', 1)):
            sub = S[S[:, 0] == sweep]
            xs = sorted(set(sub[:, col]))
            rows = []
            for x in xs:
                a = sub[(sub[:, col] == x) & (sub[:, 4] == 0)][:, 5]
                b = sub[(sub[:, col] == x) & (sub[:, 4] == 1)][:, 5]
                k = int(sub[sub[:, col] == x][0, 2])
                rows.append([('%d' % x), ('%d' % k) if col == 1 else '', '%.3f' % a.mean(), '%.3f' % b.mean(), '%.1f' % (a.sum() / b.sum())])
            L += [title, '', md_table(['K' if col == 2 else 'side n', 'emitters' if col == 1 else '', 'old s', 'new s', 'old / new'], rows), '']
            xx = sub[sub[:, 4] == 0][:, col] if col == 2 else sub[sub[:, 4] == 0][:, col] ** 2
            so = fit_slope(xx, sub[sub[:, 4] == 0][:, 5])
            sn = fit_slope(xx, sub[sub[:, 4] == 1][:, 5])
            L += ['Exponent of time in %s: old %.2f [%.2f, %.2f], new %.2f [%.2f, %.2f].' % (('K' if col == 2 else 'the number of pixels',) + so + sn), '']
    # ---- astronomy
    arows = []
    for xy in ASTRO_PATCHES:
        for cap_n in ASTRO_CAPS:
            k = 'A|x%d_y%d|%d|' % (xy + (cap_n,))
            if k + 'cpu' not in raw:
                continue
            cpu, n, tr, cm, nd = raw[k + 'cpu'], raw[k + 'n'], raw[k + 'trace'], raw[k + 'cmp'], raw[k + 'ndet']
            ok = bool(tr[0]) and tr[1] <= TOL_ASTRO_TAU and bool(tr[2]) and bool(np.all(cm[:, 0] == 1)) and float(np.nanmax(cm[:, 1])) <= TOL_ASTRO_POS and float(np.nanmax(cm[:, 2])) <= TOL_ASTRO_FLUX
            arows.append(['(%d, %d)' % xy, str(int(n[0])), '%.2f' % cpu[0], '%.2f' % cpu[1], '%.2f' % (cpu[0] / cpu[1]), '%.1f' % (1e3 * cpu[0] / n[0]), '%.1f' % (1e3 * cpu[1] / n[1]),
                          'yes' if tr[0] else 'NO', '%.1e' % tr[1] if tr[0] else '-', '%d / %d' % (int((cm[:, 0] == 1).sum()), len(cm)), '%.1e' % np.nanmax(cm[:, 1]), '%.1e' % np.nanmax(cm[:, 2]), 'pass' if ok else 'FAIL'])
    if arows:
        L += ['## 4. Astronomy development patches (M2, hp 3, r_ex 1.2, the frozen `run_path`; guarded original against new)', '',
              md_table(['patch', 'entries', 'old CPU s', 'new CPU s', 'old / new', 'old ms / entry', 'new ms / entry', 'events identical', 'max rel. entry-time diff', 'thresholds with equal counts', 'max position diff (px)',
                        'max rel. flux diff', 'verdict'], arows), '']
    crows = []
    for xy in ASTRO_PATCHES:
        for cap_n in (200, 400):
            k = 'AC|x%d_y%d|%d|' % (xy + (cap_n,))
            ka = 'A|x%d_y%d|%d|' % (xy + (cap_n,))
            if k + 'cpu' in raw and ka + 'cpu' in raw:
                tr, cm, ta, cma = raw[k + 'trace'], raw[k + 'cmp'], raw[ka + 'trace'], raw[ka + 'cmp']
                crows.append(['(%d, %d)' % xy, str(cap_n), 'yes' if tr[0] else 'NO', '%.1e' % tr[1], '%.1e' % np.nanmax(cm[:, 1]), '%.1e' % np.nanmax(cm[:, 2]),
                              'yes' if ta[0] else 'NO', '%.1e' % ta[1], '%.1e' % np.nanmax(cma[:, 1]), '%.1e' % np.nanmax(cma[:, 2])])
    if crows:
        L += ['POST HOC control (not in the plan; run after the table above showed the failures): the differences of the ORIGINAL guarded path from itself when its block-fit tolerance (1e-6) is made ten times smaller, '
              'beside the differences of the new path from the original at the same tolerance. If the first are as large as the second, the discrepancy is the iteration tolerance of the fits, not an error of the new evaluation.', '',
              md_table(['patch', 'entries', 'events identical (old vs old at tol/10)', 'max rel. entry-time diff', 'max position diff (px)', 'max rel. flux diff', 'events identical (new vs old)', 'max rel. entry-time diff', 'max position diff (px)',
                        'max rel. flux diff'], crows), '']
    for xy in ASTRO_PATCHES:
        k = 'AF|x%d_y%d|' % xy
        if k + 'cpu' in raw:
            cpu, n, tr, cm = raw[k + 'cpu'], raw[k + 'n'], raw[k + 'trace'], raw[k + 'cmp']
            L += ['Complete path of patch (%d, %d): %d entries (old) / %d (new); CPU old %.0f s, new %.0f s, old / new = %.1f; events identical: %s; thresholds with equal counts %d / %d, max position diff %.1e px, '
                  'max relative flux diff %.1e.' % (xy + (int(n[0]), int(n[1]), cpu[0], cpu[1], cpu[0] / cpu[1], 'yes' if tr[0] else 'NO', int((cm[:, 0] == 1).sum()), len(cm), np.nanmax(cm[:, 1]), np.nanmax(cm[:, 2]))), '']
    if all('AX|x%d_y%d|new|blk' % xy in raw for xy in ASTRO_PATCHES):
        L += ['POST HOC (not in the plan): equivalence of the COMPLETE paths by the F-measure at mag < 22 (completeness and purity from the 20 x 20 px block counts), the two development patches pooled; '
              'old = the guarded original, new = `imaging_fast`, tight = the original with its block-fit tolerance ten times smaller (the reference for how much a result of this algorithm moves under a tiny perturbation). '
              'Differences with a 95 %% interval of the paired bootstrap over the %d blocks of the two patches (%d resamples, seed 0).' % (2 * 25, BOOT), '']
        rows = []
        kappas = raw['AX|x%d_y%d|old|kappas' % ASTRO_PATCHES[0]]
        rng = np.random.default_rng(0)
        nb = 25 * len(ASTRO_PATCHES)
        idx = rng.integers(0, nb, size=(BOOT, nb))
        dnew, dtight, exc_new, exc_tight = [], [], 0, 0

        for j, kap in enumerate(kappas):
            blk = {e: np.concatenate([raw['AX|x%d_y%d|%s|blk' % (xy + (e,))][j] for xy in ASTRO_PATCHES]) for e in ('old', 'new', 'old_tight')}

            def F_(b):
                tot = b.sum(0) if b.ndim == 2 else b
                comp, pur = tot[1] / max(tot[0], 1), tot[3] / max(tot[2], 1)
                return 2 * comp * pur / max(comp + pur, 1e-12)
            Fv = {e: F_(b) for e, b in blk.items()}
            rep = {e: np.array([F_(blk[e][i]) for i in idx]) for e in blk}
            d1, d2 = rep['new'] - rep['old'], rep['old_tight'] - rep['old']
            dnew.append(Fv['new'] - Fv['old'])
            dtight.append(Fv['old_tight'] - Fv['old'])
            lo1, hi1, lo2, hi2 = np.percentile(d1, 2.5), np.percentile(d1, 97.5), np.percentile(d2, 2.5), np.percentile(d2, 97.5)
            exc_new += int(lo1 > 0 or hi1 < 0)
            exc_tight += int(lo2 > 0 or hi2 < 0)
            rows.append(['BIC' if kap == 0 else '%g' % kap, '%.4f' % Fv['old'], '%.4f' % Fv['new'], '%.4f' % Fv['old_tight'], '%+.4f [%+.4f, %+.4f]' % (Fv['new'] - Fv['old'], np.percentile(d1, 2.5), np.percentile(d1, 97.5)),
                         '%+.4f [%+.4f, %+.4f]' % (Fv['old_tight'] - Fv['old'], np.percentile(d2, 2.5), np.percentile(d2, 97.5))])
        L += [md_table(['threshold', 'F old', 'F new', 'F tight', 'new - old', 'tight - old'], rows), '']
        sg = lambda v: '%+.4f' % v                                    # noqa: E731
        L += ['The %d new-minus-old F differences are %s (%s to %s); %d of %d intervals exclude 0. The tight-tolerance original moves F in %s (%s to %s; %d of %d of its intervals exclude 0). So "within the '
              'perturbation scale" holds in magnitude, not in the sign pattern.' %
              (len(dnew), 'all negative' if all(d < 0 for d in dnew) else 'not all negative', sg(min(dnew)), sg(max(dnew)), exc_new, len(dnew),
               'both directions' if min(dtight) < 0 < max(dtight) else 'one direction', sg(min(dtight)), sg(max(dtight)), exc_tight, len(dtight)), '']
        ents = ['(%d, %d): %d (old) / %d (new) / %d (tight)' % (xy + tuple(int(raw['AX|x%d_y%d|%s|n' % (xy + (e,))]) for e in ('old', 'new', 'old_tight'))) for xy in ASTRO_PATCHES]
        L += ['Entry counts of the complete paths can differ between the engines: ' + '; '.join(ents) + '.', '']
    return '\n'.join(L) + '\n'


def load_raw(path=RAW):
    with np.load(path, allow_pickle=False) as z:
        return {k: z[k] for k in z.files}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', default='')
    ap.add_argument('--from-raw', action='store_true')
    ap.add_argument('--sets', nargs='*', default=list(SETS))
    ap.add_argument('--quick', action='store_true')
    A = ap.parse_args()
    parts = [p for p in A.run.split(',') if p]
    for p in parts:
        t0 = time.time()
        if p == 'E':
            part_E(A.sets, A.quick)
        elif p == 'V':
            part_V(A.sets, A.quick)
        elif p == 'T':
            part_T(A.sets, A.quick)
        elif p == 'ENV':
            part_ENV()
        elif p == 'TT':
            part_TT(A.sets, A.quick)
        elif p == 'P':
            part_P(A.sets, A.quick)
        elif p == 'S':
            part_S(A.quick)
        elif p == 'A':
            part_A(A.quick)
        elif p == 'AC':
            part_AC()
        elif p == 'AX':
            part_AX()
        elif p == 'AF':
            part_AF()
        else:
            ap.error('unknown part %r' % p)
        print('part %s done in %.0f s' % (p, time.time() - t0), flush=True)
    if parts or A.from_raw:
        with open(MD, 'w') as f:
            f.write(report(load_raw()))
        print('wrote', MD)
    else:
        ap.error('give --run or --from-raw')


if __name__ == '__main__':
    main()
