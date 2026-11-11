"""
The sliding inverse scale space on the simulated 2D single-molecule benchmarks (WP4), against Picasso
======================================================================================================

`src/core/imaging.py` on the 64 x 64 EM-CCD frames of the 2016 localisation-microscopy challenge's 2D microtubule sets (public ground
truth), scored with the harness of `experiments/smlm_common.py` (Hungarian matching at 250 nm; Jaccard, RMSE, efficiency; frame bootstrap,
paired for the differences) and compared with the Picasso ladder of `experiments/exp_smlm_picasso.py` (Picasso 0.11.3: its defaults, every option tuned on the
development frames for E and for J, each also with a position-offset correction) on the SAME test frames, whose detections are read from
`results/smlm/picasso_ladder_raw.npz`.

Protocol (nothing is tuned on test frames)
    development frames: the first 1000 (LD) or 100 (HD) frames of a set; test frames: LD 1000-1499, HD 100-499 (0-based stack indices).
    Calibrated on the development frames WITHOUT labels: the background (a quadratic surface fitted robustly to the per-pixel temporal median)
    and the noise variance var = a x + b0 (x in photons; a from the clipped per-pixel temporal variance, b0 = 0.078 photon^2 the nominal read-out
    term); the PSF width sigma0 = 1.4 px at which every kernel enters is the median of `results/smlm/calibration.md` (isolated bright emitters,
    labels used there). Tuned on the development frames WITH labels: the model (part D: high-pass width, halo) and the significance threshold
    kappa of the stop rule `stop_snr`, a bound on the fitted width (the same kind of post-filter as Picasso's width window) and a constant position offset
    (estimated on the same frames and removed), by the efficiency E of the challenge and, separately, by the Jaccard index J (two operating points per set),
    as the ladder rungs R2E and R2J of Picasso are tuned. The search is far smaller than Picasso's (about 170 operating points against 12,000 to 15,000
    scored configurations per set, `picasso_baseline.md`).

Model and algorithm (all of it in `src/core/imaging.py`)
    kernels: the Gaussian PSF integrated over each pixel, a width of its own (1.0 .. 3.2 px, sliding with the centre: the moving widths of
    Proposition D'''.59), optionally with a halo part; the model and the frame are high-passed with the same border-consistent filter (removes the
    pedestal under a bright emitter, which a Gaussian cannot represent); sub-pixel start lattice h = 0.5 px, resolution limit r_ex = 2 px; weights
    1/var from the background and the current fit (iterated). The path is stopped at the first entrant whose matched-filter significance is below
    kappa (the number of noise standard deviations by which the residual correlates with its kernel). Parts: C calibration, D (D2) development grid,
    T test frames, T0 test frames with the a priori default, B the default with the BIC stop (no threshold) on the development and test frames, G the label-free PSF width, L the label-free default (that width and the BIC stop) on the test frames, V design variants, S stop rules on the whole path, P Picasso's seconds per frame.

Usage:  python experiments/exp_smlm_iss.py --only C|D|D2|T|T0|B|G|L|V|S|P [--sets ...] [--frames N]   (parts are merged into results/smlm/iss_raw.npz)
        python experiments/exp_smlm_iss.py --from-raw                                 (regenerates results/smlm/iss.md from the raw)
Environment: ISS_OUT (directory of the raw and the md, default results/smlm), ISS_WORKERS (processes, default 3), SMLM_PICASSO_PYTHON, SMLM_DATA.
"""
import argparse
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
import smlm_common as C                      # noqa: E402
from src.core import imaging as IM           # noqa: E402

OUT = os.environ.get('ISS_OUT', os.path.join(REPO, 'results', 'smlm'))
RAW = os.path.join(OUT, 'iss_raw.npz')
MD = os.path.join(OUT, 'iss.md')
WORKERS = int(os.environ.get('ISS_WORKERS', '3'))

SETS = ('MT0.N1.LD', 'MT0.N1.HD', 'MT0.N2.LD', 'MT0.N2.HD')
PROTOCOL = {s: dict(dev=(0, 1000), test=(1000, 1500)) if s.endswith('LD') else dict(dev=(0, 100), test=(100, 500)) for s in SETS}
SIGMA0 = 1.4                     # px: median width of isolated bright emitters (results/smlm/calibration.md)
WIDTHS = (1.0, 3.2)              # px: bounds of the sliding widths
B0 = 0.078                       # photon^2: the read-out term of the variance (nominal; below the shot noise by 3 orders)
LAT_H, R_EX, MAX_ENTRIES, SNR_MIN, SNR_EXTRA = 0.5, 2.0, 40, 3.0, 3
CONFIGS = {                      # name -> high-pass width (px), kernel parts (coefficient, width multiplier)
    'hp3': dict(hp=3.0, parts=((1.0, 1.0),)),
    'hp3h': dict(hp=3.0, parts=((0.75, 1.0), (0.25, 3.0))),
    'hp5': dict(hp=5.0, parts=((1.0, 1.0),)),
    'hp5h': dict(hp=5.0, parts=((0.75, 1.0), (0.25, 3.0))),
}
SUFFIX = {'_s12': dict(sigma0=1.2), '_s17': dict(sigma0=1.7), '_r15': dict(r_ex=1.5), '_r25': dict(r_ex=2.5)}     # round 2 of the development search


def get_cfg(name):
    """A configuration by name: one of CONFIGS, optionally with a suffix of SUFFIX (the entry width sigma0 or the resolution limit r_ex changed)."""
    for suf, ov in SUFFIX.items():
        if name.endswith(suf):
            return dict(CONFIGS[name[:-len(suf)]], **ov)
    return dict(CONFIGS[name])


KAPPAS = (4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.0)
# variants of the chosen method (parts V): the overrides of its options. 'tail' keeps the whole noise-fitting tail so that the other stop rules can be read.
VARIANTS = {
    'full': {},
    'greedy': dict(order='greedy'),
    'nowidth': dict(widths=None),
    'nohp': dict(hp=0.0),
    'norex': dict(r_ex=0.0),
    'lat1': dict(lat_h=1.0),
    'fast': dict(tol=1e-4, max_nfev=25),         # looser least-squares tolerance (1e-8 -> 1e-4, 200 -> 25 evaluations): the same detections to 0.01 nm at a third of the time
    'ols': dict(flat_var=True),                  # unweighted least squares: the variance of the empty frame for every pixel and every iterate (the model's photon noise dropped)
    'sig12': dict(sigma0=1.2),
    'sig17': dict(sigma0=1.7),
}
TAIL_RULES = [('disc', 2.0), ('disc', 3.0), ('gap', None), ('gapg', 5.0), ('gapg', 6.0)]
DEFAULT = dict(cfg='hp5h', kappa=5.0)          # the a priori configuration, declared before its test run: the halo model with the 5 px high-pass, the entry width of the
                                               # calibration (1.4 px), the 2 px resolution limit, the conventional 5-sigma significance, no width bound, no position offset
LADDER = os.path.join(REPO, 'results', 'smlm', 'picasso_ladder_raw.npz')          # the Picasso ladder's detections (exp_smlm_picasso.py)
LADDER_ROWS = [('R0', 'Picasso, defaults (R0)'), ('R1E+c', 'Picasso, only the net gradient tuned for E, + offset correction (R1E+c)'), ('R2E', 'Picasso, every option tuned for E (R2E)'), ('R2E+c', 'Picasso, R2E + offset correction'),
               ('R2J', 'Picasso, every option tuned for J (R2J)'), ('R2J+c', 'Picasso, R2J + offset correction')]
LADDER_VS = {'E': ('R0', 'R2E', 'R2E+c'), 'J': ('R0', 'R2J', 'R2J+c')}             # the rows each operating point is compared with
BOOT = 1000
BLOCK = 10                       # frames per resampled block: an emitter stays on for ~2.3 frames, so neighbouring frames are not independent


# ----------------------------------------------------------------------------------------------------------------------------------
# calibration (labels not used)
# ----------------------------------------------------------------------------------------------------------------------------------
def fit_surface(med, order=2, iters=6, reject=2.5):
    """A robust polynomial surface through the per-pixel temporal median `med` (ny, nx): least squares of the pixels below median + reject * MAD
    of the residual (the contamination by emitters is one-sided), iterated. Returns the surface (ny, nx)."""
    ny, nx = med.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    u, v = (xx - (nx - 1) / 2) / nx, (yy - (ny - 1) / 2) / ny
    cols = [u ** i * v ** j for i in range(order + 1) for j in range(order + 1 - i)]
    A = np.stack([c.ravel() for c in cols], axis=1)
    keep = np.ones(med.size, bool)
    for _ in range(iters):
        coef = np.linalg.lstsq(A[keep], med.ravel()[keep], rcond=None)[0]
        res = med.ravel() - A @ coef
        mad = 1.4826 * np.median(np.abs(res[keep] - np.median(res[keep])))
        keep = res < np.median(res[keep]) + reject * mad
    return (A @ coef).reshape(ny, nx)


def calibrate(dev_photons):
    """Label-free calibration from development frames in photons (n, ny, nx): background surface and the variance slope a of var = a x + B0.
    The temporal variance of each pixel is taken after clipping at 4 robust standard deviations of its own series (twice), and a is the median over
    the pixels of (var - B0) / background."""
    med = np.median(dev_photons, axis=0)
    bg = fit_surface(med)
    mad = 1.4826 * np.median(np.abs(dev_photons - med), axis=0)
    keep = np.ones(dev_photons.shape, bool)
    for _ in range(2):
        keep = np.abs(dev_photons - med) < 4.0 * np.maximum(mad, 1e-9)
    x = np.where(keep, dev_photons, np.nan)
    var = np.nanvar(x, axis=0, ddof=1)
    a = float(np.median((var - B0) / np.maximum(bg, 1e-9)))
    return dict(bg=bg, a=a, bg_mean=float(bg.mean()), sd=float(np.sqrt(np.median(var))))


def _pix_gauss(x0, y0, s, half):
    """The pixel-integrated unit-flux Gaussian on the (2 half + 1)^2 patch, centre offsets (x0, y0) from the middle pixel."""
    from scipy.special import erf
    p = np.arange(-half, half + 1, dtype=float)
    sq = np.sqrt(2.0) * s
    gx = 0.5 * (erf((p + 0.5 - x0) / sq) - erf((p - 0.5 - x0) / sq))
    gy = 0.5 * (erf((p + 0.5 - y0) / sq) - erf((p - 0.5 - y0) / sq))
    return np.outer(gy, gx)


def estimate_sigma0(dev_photons, bg, a, b0=B0, n_frames=300, snr_peak=8.0, iso_px=7.0, top=0.25, half=5, chi2_max=2.0):
    """The PSF core width (px) from development frames WITHOUT labels: in every frame (background removed) the local maxima of the frame smoothed by 1 px whose raw peak is at
    least snr_peak noise sd above the background and that have no other maximum within iso_px are isolated bright emitters; each is fitted by a pixel-integrated Gaussian plus an
    offset on a (2 half + 1)^2 window (weights from the noise law a x + b0, refitted once with the variance of the fitted model; fits of reduced chi^2 above chi2_max, e.g. unresolved blends, are dropped); the estimate is the median fitted width of the brightest `top` fraction. Returns (sigma0, number of spots
    fitted, 25 %, 75 % of the widths of the brightest); nan if no spot."""
    from scipy.ndimage import gaussian_filter, maximum_filter
    from scipy.optimize import least_squares
    fr = np.asarray(dev_photons, float)[:n_frames] - bg
    noise = np.sqrt(a * bg + b0)
    thr_sm = 0.28 * float(np.median(noise))                       # the sd of the 1 px-smoothed noise is 1 / (2 sqrt(pi)) of the pixel sd
    ny, nx = fr.shape[1:]
    spots = []
    for f in fr:
        sm = gaussian_filter(f, 1.0)
        mx = (sm == maximum_filter(sm, size=5)) & (sm > 3.0 * thr_sm)
        ys, xs = np.nonzero(mx)
        for y, x in zip(ys, xs):
            if y < half or x < half or y >= ny - half or x >= nx - half or f[y, x] < snr_peak * noise[y, x]:
                continue
            if np.sum((ys - y) ** 2 + (xs - x) ** 2 <= iso_px ** 2) > 1:      # another local maximum within iso_px: not isolated
                continue
            patch = f[y - half:y + half + 1, x - half:x + half + 1]
            bgp = bg[y - half:y + half + 1, x - half:x + half + 1]
            w = 1.0 / noise[y - half:y + half + 1, x - half:x + half + 1]
            z0, ok = np.array([float(patch.sum()), 0.0, 0.0, 1.5, 0.0]), False
            for _pass in range(2):                                      # the second pass weights by the variance of the fitted model (a bright spot is noisier than the background)
                sol = least_squares(lambda z: (((z[0] * _pix_gauss(z[1], z[2], z[3], half) + z[4]) - patch) * w).ravel(), z0,
                                    bounds=([0, -2, -2, 0.6, -50], [np.inf, 2, 2, 4.0, 50]), method='trf')
                z0 = sol.x
                w = 1.0 / np.sqrt(a * (bgp + np.maximum(sol.x[0] * _pix_gauss(sol.x[1], sol.x[2], sol.x[3], half) + sol.x[4], 0.0)) + b0)
            chi2 = float(np.sum(((sol.x[0] * _pix_gauss(sol.x[1], sol.x[2], sol.x[3], half) + sol.x[4] - patch) * w) ** 2)) / (patch.size - 5)
            if sol.success and 0.7 < sol.x[3] < 3.9 and chi2 < chi2_max:        # a blend of two emitters fits a Gaussian badly: rejected
                spots.append((sol.x[0], sol.x[3]))
    if not spots:
        return float('nan'), 0, float('nan'), float('nan')
    spots = np.array(spots)
    k = max(int(np.ceil(top * len(spots))), 1)
    wid = spots[np.argsort(-spots[:, 0])[:k], 1]
    return float(np.median(wid)), len(spots), float(np.percentile(wid, 25)), float(np.percentile(wid, 75))


# ----------------------------------------------------------------------------------------------------------------------------------
# the method on one frame
# ----------------------------------------------------------------------------------------------------------------------------------
_W = {}


def _init(cfg_name, bg, a, opts):
    cfg = get_cfg(cfg_name)
    cfg.update({k: opts[k] for k in ('hp', 'sigma0', 'r_ex') if k in opts})
    sig0 = cfg.get('sigma0', SIGMA0)
    P = IM.PixelPSF(64, 64, sig0, parts=cfg['parts'], hp=cfg['hp'] or None)
    _W.update(P=P, Pp=IM.PixelPSF(64, 64, sig0, parts=cfg['parts']), L=IM.Lattice(P, h=opts.get('lat_h', LAT_H), margin=1.0), bg=bg, a=a, opts=opts, r_ex=cfg.get('r_ex', R_EX))


def _sw(c, w):
    """The weights 1/sqrt(variance) of the pixels: var = a (bg + model) + b0 (the variant 'ols' drops the model: a fixed variance, the unweighted fit)."""
    m = 0.0 if _W['opts'].get('flat_var') else (_W['Pp'].synth(c, w) if len(w) else 0.0)
    return 1.0 / np.sqrt(_W['a'] * (_W['bg'] + np.maximum(m, 0.0)) + B0)


def _rule_block(rule, p):
    kind, par = rule
    if kind == 'snr':
        return IM.stop_snr(p['hist'], p['taus'], p['snr'], par)[0]
    if kind == 'disc':
        return IM.stop_discrepancy(p['hist'], 64 * 64, z=par)
    if kind == 'gap':
        return IM.stop_gap(p['hist'], p['taus'])[0]
    if kind == 'gapg':
        return IM.stop_gap(p['hist'], p['taus'], snr=p['snr'], gate=par)[0]
    if kind == 'bic':                                       # no threshold: par = the parameters of a unit (position, width, amplitude: 4)
        return IM.stop_bic(p['hist'], 64 * 64, n_par=int(par))[0]
    raise ValueError(rule)


def run_frame(args):
    """One frame (photons, background not yet removed): the path, then every rule's block refined against the whole frame.
    Returns (per rule: array (k, 4) = x_nm, y_nm, photons, sigma_px), the entries' times and significances, the seconds."""
    k, frame = args
    t0 = time.perf_counter()
    P, L, bg, o = _W['P'], _W['L'], _W['bg'], _W['opts']
    widths = o.get('widths', WIDTHS)
    nc = 3 if widths else 2
    y0 = frame - bg
    y = P.highpass(y0 - np.median(y0))
    sw0 = _sw(np.zeros((0, nc)), np.zeros(0))
    R = IM.sliding_path(P, y, sw0, L, max_entries=o.get('max_entries', MAX_ENTRIES), sw_fn=_sw, sw_on='params', sig_bounds=widths, tol=o.get('tol', 1e-8), max_nfev=o.get('max_nfev', 200),
                        passes=1, rtrunc=4.0, r_ex=_W['r_ex'], snr_min=o.get('snr_min', SNR_MIN), snr_extra=SNR_EXTRA,
                        jump=o.get('jump'), extra=6, order=o.get('order', 'dual'))
    p = dict(hist=R['hist'], taus=R['taus'], snr=R['snr'])
    cache, out = {}, []
    for rule in o['rules']:
        i = _rule_block(rule, p)
        if i not in cache:
            h = R['hist'][i]
            cc, ww = h['c'], h['w']
            if len(ww):
                cc, ww, _ = IM.refine(P, y, sw0, cc, ww, sig_bounds=widths, box=1.0, sw_fn=_sw, sw_on='params')
            sg = cc[:, 2] if nc == 3 else np.full(len(ww), P.sigma)
            cache[i] = np.c_[(cc[:, 0] + 0.5) * C.PIXEL_NM, (cc[:, 1] + 0.5) * C.PIXEL_NM, ww, sg] if len(ww) else np.zeros((0, 4))
        out.append(cache[i])
    return out, R['taus'], R['snr'], time.perf_counter() - t0


def run_frames(frames, cfg_name, cal, opts, workers=WORKERS):
    """The method on frames (n, 64, 64) in photons. Returns (dets[rule][frame] = (k, 4) arrays, taus list, snr list, seconds per frame)."""
    from multiprocessing import Pool
    jobs = list(enumerate(frames))
    if workers <= 1:
        _init(cfg_name, cal['bg'], cal['a'], opts)
        res = [run_frame(j) for j in jobs]
    else:
        with Pool(workers, initializer=_init, initargs=(cfg_name, cal['bg'], cal['a'], opts)) as pool:
            res = pool.map(run_frame, jobs, chunksize=4)
    nr = len(opts['rules'])
    dets = [[r[0][i] for r in res] for i in range(nr)]
    return dets, [r[1] for r in res], [r[2] for r in res], np.array([r[3] for r in res])


# ----------------------------------------------------------------------------------------------------------------------------------
# raw storage: ragged per-frame detections as offsets + one flat array
# ----------------------------------------------------------------------------------------------------------------------------------
def pack(dets):
    """list over frames of (k, 4) arrays -> (offsets (F+1,), flat (M, 4))."""
    off = np.r_[0, np.cumsum([len(d) for d in dets])].astype(np.int64)
    flat = np.concatenate([d.reshape(-1, 4) for d in dets]) if off[-1] else np.zeros((0, 4))
    return off, flat


def unpack(off, flat):
    return [flat[off[i]:off[i + 1]] for i in range(len(off) - 1)]


def save_part(path, part, arrays):
    """Merge the arrays of one part into the raw file (keys are prefixed with the part's name; existing keys of that part are replaced)."""
    old = {}
    if os.path.exists(path):
        with np.load(path, allow_pickle=False) as z:
            old = {k: z[k] for k in z.files if not k.startswith(part + '|')}
    old.update({part + '|' + k: np.asarray(v) for k, v in arrays.items()})
    os.makedirs(os.path.dirname(path), exist_ok=True)
    np.savez_compressed(path, **old)


def load_raw(path=RAW):
    with np.load(path, allow_pickle=False) as z:
        return {k: z[k] for k in z.files}


# ----------------------------------------------------------------------------------------------------------------------------------
# parts
# ----------------------------------------------------------------------------------------------------------------------------------
def dev_test_frames(name, which, n_max=None):
    d0, d1 = PROTOCOL[name][which]
    if n_max:
        d1 = min(d1, d0 + n_max)
    return d0, d1


def load_set(name):
    stack = C.data_path(name, 'stack')
    return stack, C.load_gt(C.data_path(name, 'gt'))


def part_C(sets, out_path=RAW):
    arrays = {}
    for name in sets:
        stack, _ = load_set(name)
        d0, d1 = PROTOCOL[name]['dev']
        cal = calibrate(C.to_photons(C.load_frames(stack, d0, d1)))
        print('%s: background %.2f photons (surface), a = %.3f, sd %.2f' % (name, cal['bg_mean'], cal['a'], cal['sd']), flush=True)
        arrays.update({name + '|bg': cal['bg'], name + '|a': [cal['a']], name + '|sd': [cal['sd']], name + '|n': [d1 - d0]})
    save_part(out_path, 'C', arrays)


CAL_FROM = {'MT0.N1.HD': 'MT0.N1.LD', 'MT0.N2.HD': 'MT0.N2.LD'}      # the dense sets borrow the camera of the sparse set of their noise class


def get_cal(raw, name):
    """The calibration used for a set. A dense set (HD) does not calibrate itself: the temporal median and variance of a pixel are contaminated by the
    emitters that cover a third of the field in every frame (the HD medians range from 94 to 242 photons against 97.5 to 104.7 in LD), so it uses the
    camera of the low-density set of its noise class (the same simulator parameters, `calibration.md`); `part_C` records both."""
    src = CAL_FROM.get(name, name)
    return dict(bg=raw['C|%s|bg' % src], a=float(raw['C|%s|a' % src][0]))


def rules_dev():
    return [('snr', k) for k in KAPPAS]


def part_D(sets, cfgs, n_frames=None, out_path=RAW):
    raw = load_raw(out_path)
    for name in sets:
        stack, gt = load_set(name)
        d0, d1 = dev_test_frames(name, 'dev', n_frames)
        frames = C.to_photons(C.load_frames(stack, d0, d1))
        cal = get_cal(raw, name)
        for cfg in cfgs:
            t0 = time.time()
            dets, taus, snr, secs = run_frames(frames, cfg, cal, dict(rules=rules_dev()))
            arrays = {}
            for r, rule in enumerate(rules_dev()):
                off, flat = pack(dets[r])
                arrays['%g|off' % rule[1]] = off
                arrays['%g|flat' % rule[1]] = flat
            arrays['secs'] = secs
            arrays['frames'] = [d0, d1]
            save_part(out_path, 'D|' + name + '|' + cfg, arrays)
            sc = C.score(dets[KAPPAS.index(5.0)], gt, np.arange(d0, d1))
            print('D %s %s: %d frames in %.0f s; kappa 5: J %.3f RMSE %.1f E %.1f' % (name, cfg, d1 - d0, time.time() - t0, sc['jaccard'], sc['rmse_nm'], sc['efficiency']),
                  flush=True)


def part_D2(sets, n_frames=None, out_path=RAW):
    """Round 2 of the development search: around the E-optimal base configuration of each set, the entry width sigma0 = 1.2, 1.7 and the resolution limit
    r_ex = 1.5, 2.5 (coordinate ascent; the same significance thresholds)."""
    raw = load_raw(out_path)
    for name in sets:
        stack, gt = load_set(name)
        rows = [r for r in op_table(raw, name, gt) if r['cfg'] in CONFIGS]
        base = choose_op(rows, 'E')['cfg']
        d0, d1 = dev_test_frames(name, 'dev', n_frames)
        frames = C.to_photons(C.load_frames(stack, d0, d1))
        for suf in SUFFIX:
            cfg = base + suf
            t0 = time.time()
            dets, taus, snr, secs = run_frames(frames, cfg, get_cal(raw, name), dict(rules=rules_dev()))
            arrays = {}
            for r, rule in enumerate(rules_dev()):
                off, flat = pack(dets[r])
                arrays['%g|off' % rule[1]] = off
                arrays['%g|flat' % rule[1]] = flat
            arrays['secs'] = secs
            arrays['frames'] = [d0, d1]
            save_part(out_path, 'D|' + name + '|' + cfg, arrays)
            print('D2 %s %s: %d frames in %.0f s' % (name, cfg, d1 - d0, time.time() - t0), flush=True)


def metrics_row(sc, sc2):
    return [sc['jaccard'], sc['precision'], sc['recall'], sc['rmse_nm'], sc['efficiency'], sc2['jaccard'], sc2['efficiency']]


SMAX = (2.0, 2.5, 3.2)           # post-filter on the fitted width (px), the same kind of filter as Picasso's; 3.2 is the bound of the widths: none
OPS = ('E', 'J')                 # two operating points per set, chosen on the development frames by the efficiency and by the Jaccard index


def post(dets, smax):
    """Detections whose fitted width is at most smax (px)."""
    return [d[d[:, 3] <= smax + 1e-9] if len(d) else d for d in dets]


def shift(dets, off):
    """Detections moved by -off (nm): the position offset of a symmetric fit against the simulated emitter (estimator- and z-dependent, see
    `calibration.md`) estimated on the development frames."""
    return [np.c_[d[:, 0] - off[0], d[:, 1] - off[1], d[:, 2:]] if len(d) else d for d in dets]


def offset_of(dets, gt, idx):
    """The mean signed error (dx, dy) in nm of the matched pairs on the given frames (detection minus truth); zero where nothing matched."""
    r = C.score(dets, gt, idx)
    return np.nan_to_num(np.array([r['bias_x_nm'], r['bias_y_nm']]))


def cfg_names(raw, name):
    """The configurations with development results for a set, in a fixed order (the four of CONFIGS, then their round-2 variants)."""
    pre, suf = 'D|%s|' % name, '|frames'
    have = {k[len(pre):-len(suf)] for k in raw if k.startswith(pre) and k.endswith(suf)}
    order = list(CONFIGS) + [b + x for b in CONFIGS for x in SUFFIX]
    return [c for c in order if c in have]


def op_table(raw, name, gt):
    """Development-frame metrics of every operating point (config, kappa, s_max) with the offset estimated on the same frames and removed: a list of
    dict(cfg, kappa, smax, off, m = [J, P, R, RMSE, E, J_challenge, E_challenge])."""
    pmin = C.photon_threshold(gt, 0.75)
    rows = []
    for cfg in cfg_names(raw, name):
        key = 'D|%s|%s|' % (name, cfg)
        d0, d1 = (int(v) for v in raw[key + 'frames'])
        idx = np.arange(d0, d1)
        for kap in KAPPAS:
            base = unpack(raw[key + '%g|off' % kap], raw[key + '%g|flat' % kap])
            for smax in SMAX:
                dets = post(base, smax)
                off = offset_of(dets, gt, idx)
                dc = shift(dets, off)
                rows.append(dict(cfg=cfg, kappa=kap, smax=smax, off=off, m=metrics_row(C.score(dc, gt, idx), C.score(dc, gt, idx, border_nm=450.0, photon_min=pmin))))
    return rows


def choose_op(rows, which):
    """The operating point of the largest development efficiency ('E'; ties: larger kappa, smaller width bound) or Jaccard index ('J'; ties: smaller RMSE)."""
    if which == 'E':
        return max(rows, key=lambda r: (round(r['m'][4], 6), r['kappa'], -r['smax']))
    return max(rows, key=lambda r: (round(r['m'][0], 6), -r['m'][3], r['kappa']))


def part_T(sets, n_frames=None, out_path=RAW):
    """The test frames, run once for every configuration that an operating point (E, J) needs, with all significance thresholds stored."""
    raw = load_raw(out_path)
    for name in sets:
        stack, gt = load_set(name)
        ops = {w: choose_op(op_table(raw, name, gt), w) for w in OPS}
        cfgs = sorted({o['cfg'] for o in ops.values()})
        a, b = dev_test_frames(name, 'test', n_frames)
        frames = C.to_photons(C.load_frames(stack, a, b))
        arrays = {'frames': [a, b]}
        for cfg in cfgs:
            t0 = time.time()
            dets, taus, snr, secs = run_frames(frames, cfg, get_cal(raw, name), dict(rules=rules_dev()))
            _store(arrays, cfg, rules_dev(), dets)
            arrays['%s|secs' % cfg] = secs
            print('T %s %s: %d frames, %.0f s, %.2f s/frame' % (name, cfg, b - a, time.time() - t0, secs.mean()), flush=True)
        save_part(out_path, 'T|' + name, arrays)


def part_G(sets, out_path=RAW):
    """The label-free PSF width of each camera: `estimate_sigma0` on the first 300 development frames of the sparse (LD) set of its noise class (the dense sets use their sparse set's
    camera, as in `get_cal`); stored as G|<set>|sigma0 = [estimate, spots, 25 %, 75 %]."""
    raw = load_raw(out_path)
    arrays = {}
    for name in sets:
        src = CAL_FROM.get(name, name)
        stack, _ = load_set(src)
        d0, d1 = PROTOCOL[src]['dev']
        cal = get_cal(raw, src)
        t0 = time.time()
        r = estimate_sigma0(C.to_photons(C.load_frames(stack, d0, min(d1, d0 + 300))), cal['bg'], cal['a'])
        arrays['%s|sigma0' % name] = np.array(r)
        print('G %s (camera of %s): sigma0 %.3f from %d isolated bright emitters (brightest quartile %.3f-%.3f), %.0f s' % (name, src, r[0], r[1], r[2], r[3], time.time() - t0), flush=True)
    save_part(out_path, 'G', arrays)


def part_L(sets, n_frames=None, out_path=RAW):
    """The label-free default on the test frames: the default configuration with the entry width of `part_G` (no labels anywhere) and the BIC stop (nothing to tune)."""
    raw = load_raw(out_path)
    for name in sets:
        stack, gt = load_set(name)
        sg = float(raw['G|%s|sigma0' % name][0])
        a, b = dev_test_frames(name, 'test', n_frames)
        frames = C.to_photons(C.load_frames(stack, a, b))
        rules = [('bic', 4)]
        t0 = time.time()
        dets, taus, snr, secs = run_frames(frames, DEFAULT['cfg'], get_cal(raw, name), dict(rules=rules, sigma0=sg))
        arrays = {'frames': [a, b], 'secs': secs, 'sigma0': [sg]}
        _store(arrays, 'run', rules, dets)
        save_part(out_path, 'L|' + name, arrays)
        print('L %s sigma0 %.3f: %d frames, %.0f s, %.2f s/frame' % (name, sg, b - a, time.time() - t0, secs.mean()), flush=True)


def part_B(sets, n_frames=None, out_path=RAW):
    """The a priori default (`DEFAULT`'s configuration) with the BIC stop instead of a significance level: a unit is kept if it lowers the weighted sum of squares of the frame by
    more than 4 ln(4096) = 33.3 (four parameters: position, width, amplitude), a rule with nothing to tune; on the development frames (to see it against the best level) and on the test frames."""
    raw = load_raw(out_path)
    for name in sets:
        stack, gt = load_set(name)
        arrays = {}
        for which in ('dev', 'test'):
            a, b = dev_test_frames(name, which, n_frames)
            frames = C.to_photons(C.load_frames(stack, a, b))
            rules = [('bic', 4)]
            t0 = time.time()
            dets, taus, snr, secs = run_frames(frames, DEFAULT['cfg'], get_cal(raw, name), dict(rules=rules))
            _store(arrays, which, rules, dets)
            arrays[which + '|frames'] = [a, b]
            arrays[which + '|secs'] = secs
            print('B %s %s: %d frames, %.0f s, %.2f s/frame' % (name, which, b - a, time.time() - t0, secs.mean()), flush=True)
        save_part(out_path, 'B|' + name, arrays)


def part_T0(sets, n_frames=None, out_path=RAW):
    """The test frames with the a priori default (`DEFAULT`): no operating point chosen on any frames, the counterpart of Picasso's defaults (R0)."""
    raw = load_raw(out_path)
    for name in sets:
        stack, gt = load_set(name)
        a, b = dev_test_frames(name, 'test', n_frames)
        frames = C.to_photons(C.load_frames(stack, a, b))
        t0 = time.time()
        rules = [('snr', DEFAULT['kappa'])]
        dets, taus, snr, secs = run_frames(frames, DEFAULT['cfg'], get_cal(raw, name), dict(rules=rules))
        arrays = {'frames': [a, b], 'secs': secs}
        _store(arrays, 'run', rules, dets)
        save_part(out_path, 'T0|' + name, arrays)
        print('T0 %s %s kappa %g: %d frames, %.0f s, %.2f s/frame' % (name, DEFAULT['cfg'], DEFAULT['kappa'], b - a, time.time() - t0, secs.mean()), flush=True)


def _store(arrays, tag, rules, dets):
    for r, rule in enumerate(rules):
        off, flat = pack(dets[r])
        nm = '%s|%s%s' % (tag, rule[0], '' if rule[1] is None else '%g' % rule[1])
        arrays[nm + '|off'] = off
        arrays[nm + '|flat'] = flat


def op_test_dets(raw, name, o):
    """The test-frame detections of an operating point: its configuration's run at its threshold, width-filtered, shifted by its dev offset."""
    key = 'T|%s|%s|snr%g|' % (name, o['cfg'], o['kappa'])
    return shift(post(unpack(raw[key + 'off'], raw[key + 'flat']), o['smax']), o['off'])


def part_V(sets, variants, n_frames=None, out_path=RAW):
    """The variants of the E-tuned method on the development frames (every threshold) and on the test frames."""
    raw = load_raw(out_path)
    for name in sets:
        stack, gt = load_set(name)
        cfg = choose_op(op_table(raw, name, gt), 'E')['cfg']
        cal = get_cal(raw, name)
        for vname in variants:
            t0 = time.time()
            ov = VARIANTS[vname]
            arrays = {}
            for which in ('dev', 'test'):
                a, b = dev_test_frames(name, which, n_frames)
                frames = C.to_photons(C.load_frames(stack, a, b))
                dets, taus, snr, secs = run_frames(frames, cfg, cal, dict(ov, rules=rules_dev()))
                _store(arrays, which, rules_dev(), dets)
                arrays[which + '|frames'] = [a, b]
                arrays[which + '|secs'] = secs
            save_part(out_path, 'V|%s|%s' % (name, vname), arrays)
            print('V %s %s (%s): %.0f s' % (name, vname, cfg, time.time() - t0), flush=True)


def part_S(sets, n_frames=None, out_path=RAW):
    """The stop rules on the test frames with the whole noise-fitting tail of the E-tuned method (no early stop): the significance rule at its tuned
    threshold against the discrepancy rule (z = 2, 3), the largest-gap rule and the gap rule gated by the first entrant's significance (5, 6)."""
    raw = load_raw(out_path)
    for name in sets:
        stack, gt = load_set(name)
        o = choose_op(op_table(raw, name, gt), 'E')
        a, b = dev_test_frames(name, 'test', n_frames if n_frames else (300 if name.endswith('LD') else 150))
        frames = C.to_photons(C.load_frames(stack, a, b))
        rules = [('snr', o['kappa'])] + TAIL_RULES
        t0 = time.time()
        dets, taus, snr, secs = run_frames(frames, o['cfg'], get_cal(raw, name), dict(rules=rules, snr_min=None, max_entries=25 if name.endswith('LD') else 40))
        arrays = {'frames': [a, b], 'secs': secs, 'kappa': [o['kappa']], 'smax': [o['smax']], 'off': o['off']}
        _store(arrays, 'rules', rules, dets)
        save_part(out_path, 'S|' + name, arrays)
        print('S %s (%s, kappa %g): %d frames, %.0f s' % (name, o['cfg'], o['kappa'], b - a, time.time() - t0), flush=True)


def part_P(sets, n_frames=None, out_path=RAW, python=None):
    """Picasso's seconds per frame with its defaults (one core; identification and fit, JIT compilation excluded) on the first 100 test frames of a set. The
    Picasso detections themselves come from the ladder (`LADDER`)."""
    for name in sets:
        stack, gt = load_set(name)
        a, b = dev_test_frames(name, 'test', n_frames)
        b = min(b, a + 100)
        info = {}
        C.run_picasso(C.load_frames(stack, a, b), params=None, python=python, info=info)
        save_part(out_path, 'P|' + name, {'secs': [info['run_s'], info['n_frames']]})
        print('P %s: %.4f s / frame (%d frames)' % (name, info['run_s'] / info['n_frames'], info['n_frames']), flush=True)


# ----------------------------------------------------------------------------------------------------------------------------------
# report
# ----------------------------------------------------------------------------------------------------------------------------------
def fmt_ci(v, ci, scale=100.0, nd=1):
    return '%.*f [%.*f, %.*f]' % (nd, scale * v, nd, scale * ci[0], nd, scale * ci[1])


def table(header, rows):
    out = ['| ' + ' | '.join(header) + ' |', '|' + '|'.join(['---'] * len(header)) + '|']
    return '\n'.join(out + ['| ' + ' | '.join(str(c) for c in r) + ' |' for r in rows])


def report(raw, boot=BOOT, sets=SETS):
    L = ['# The sliding inverse scale space on the 2D single-molecule benchmarks, against Picasso', '',
         'Generated by `experiments/exp_smlm_iss.py` from `iss_raw.npz`; every number below is printed by the script. Protocol, model and algorithm: the '
         'script\'s docstring. Frame indices are 0-based stack indices, positions in nm (corner convention), metrics pooled over the frames with 95 %% '
         'frame-bootstrap intervals in blocks of %d consecutive frames (%d resamples, seed 0; emitters stay on for about 2.3 frames, so frames are not independent); '
         'differences are paired (the same resampled blocks for both methods).' % (BLOCK, boot), '']
    summary, tuning = [], []
    for name in sets:
        if 'C|%s|bg' % name not in raw:
            continue
        gt = C.load_gt(C.data_path(name, 'gt'))
        pmin = C.photon_threshold(gt, 0.75)
        L += ['## %s' % name, '']
        cal = get_cal(raw, name)
        src = CAL_FROM.get(name, name)
        L.append('Calibration (development frames %d-%d of %s, no labels%s): background surface mean %.2f photons (range %.2f to %.2f), variance slope a = %.3f '
                 '(var = a x + %.3f), PSF width at entry sigma0 = %.2f px, widths slide in [%.1f, %.1f] px. Own estimate of this set: background %.2f photons, '
                 'a = %.3f.\n' % (PROTOCOL[src]['dev'][0], PROTOCOL[src]['dev'][1] - 1, src, '' if src == name else '; the dense set uses the camera of its sparse set',
                                   cal['bg'].mean(), cal['bg'].min(), cal['bg'].max(), cal['a'], B0, SIGMA0, *WIDTHS, raw['C|%s|bg' % name].mean(),
                                   float(raw['C|%s|a' % name][0])))
        rows = op_table(raw, name, gt)
        if not rows:
            continue
        ops = {w: choose_op(rows, w) for w in OPS}
        L.append('### Development frames: operating points\n')
        L.append('%d operating points (%d model configurations x %d significance thresholds x %d width bounds), each with the position offset estimated on the same '
                 'frames and removed. Bold: the point of the largest efficiency (E) and the one of the largest Jaccard index (J), the two frozen for the test frames.\n' % (
                     len(rows), len({r['cfg'] for r in rows}), len(KAPPAS), len(SMAX)))
        shown = [r for r in rows if r['smax'] == 3.2 or r is ops['E'] or r is ops['J']]
        L.append(table(['config', 'kappa', 'width <=', 'offset (dx, dy) nm', 'J %', 'precision %', 'recall %', 'RMSE nm', 'efficiency', 'J % (challenge)', 'efficiency (challenge)'],
                       [[('**%s**' if (r is ops['E'] or r is ops['J']) else '%s') % r['cfg'], '%g' % r['kappa'], '%g' % r['smax'], '(%.1f, %.1f)' % tuple(r['off']),
                         '%.1f' % (100 * r['m'][0]), '%.1f' % (100 * r['m'][1]), '%.1f' % (100 * r['m'][2]), '%.1f' % r['m'][3], '%.1f' % r['m'][4], '%.1f' % (100 * r['m'][5]),
                         '%.1f' % r['m'][6]] for r in shown]))
        if 'B|%s|dev|frames' % name in raw and 'D|%s|%s|frames' % (name, DEFAULT['cfg']) in raw:
            L += bic_dev_lines(raw, name, gt)
        tkey = 'T|%s|' % name
        if tkey + 'frames' in raw and os.path.exists(LADDER):
            a, b = (int(v) for v in raw[tkey + 'frames'])
            idx = np.arange(a, b)
            meths = []
            for row, label in LADDER_ROWS:
                pdets, pfr = C.load_ladder_detections(LADDER, name, row, 'test')
                if pfr[0] != a or pfr[-1] != b - 1:
                    raise ValueError('the ladder\'s test frames %d-%d differ from this run\'s %d-%d' % (pfr[0], pfr[-1], a, b - 1))
                meths.append((row, label, pdets))
            if 'T0|%s|frames' % name in raw:
                k0 = 'T0|%s|run|snr%g|' % (name, DEFAULT['kappa'])
                meths.append(('ISS-0', 'sliding ISS, a priori default (%s, kappa %g, no width bound, no offset correction)' % (DEFAULT['cfg'], DEFAULT['kappa']),
                              unpack(raw[k0 + 'off'], raw[k0 + 'flat'])))
            if 'L|%s|frames' % name in raw:
                kl = 'L|%s|run|bic4|' % name
                meths.append(('ISS-L', 'sliding ISS, label-free default (%s, BIC stop, width %.2f px estimated from the frames, no width bound, no offset correction)' % (DEFAULT['cfg'], raw['L|%s|sigma0' % name][0]),
                              unpack(raw[kl + 'off'], raw[kl + 'flat'])))
            if 'B|%s|test|frames' % name in raw:
                kb = 'B|%s|test|bic4|' % name
                meths.append(('ISS-B', 'sliding ISS, BIC stop (%s, nothing to tune, no width bound, no offset correction)' % DEFAULT['cfg'], unpack(raw[kb + 'off'], raw[kb + 'flat'])))
            for w in OPS:
                o = ops[w]
                meths.append(('ISS-' + w, 'sliding ISS, %s-tuned (%s, kappa %g, width <= %g px, offset removed)' % (w, o['cfg'], o['kappa'], o['smax']), op_test_dets(raw, name, o)))
            L.append('\n### Test frames %d-%d\n' % (a, b - 1))
            res = {}
            for key, mname, dets in meths:
                r = C.score(dets, gt, idx)
                r2 = C.score(dets, gt, idx, border_nm=450.0, photon_min=pmin)
                res[key] = (mname, r, r2, C.bootstrap_ci(r, B=boot, seed=0, block=BLOCK), C.bootstrap_ci(r2, B=boot, seed=0, block=BLOCK))
            L.append('All ground-truth emitters (primary setting):\n')
            L.append(table(['method', 'J %', 'precision %', 'recall %', 'RMSE nm', 'efficiency', 'detected / frame', 'true / frame'],
                           [[m, fmt_ci(r['jaccard'], ci['jaccard']), fmt_ci(r['precision'], ci['precision']), fmt_ci(r['recall'], ci['recall']),
                             fmt_ci(r['rmse_nm'], ci['rmse_nm'], 1.0), fmt_ci(r['efficiency'], ci['efficiency'], 1.0), '%.3f' % (r['n_det'] / len(idx)),
                             '%.3f' % (r['n_gt'] / len(idx))] for m, r, r2, ci, ci2 in res.values()]))
            L.append('\nChallenge assessment setting (emitters within 450 nm of the border and the dimmest quarter of the activations dropped after the pairing; the rules '
                     'are those of `smlm_common.score`, reconstructed from the challenge description):\n')
            L.append(table(['method', 'J %', 'RMSE nm', 'efficiency'],
                           [[m, fmt_ci(r2['jaccard'], ci2['jaccard']), fmt_ci(r2['rmse_nm'], ci2['rmse_nm'], 1.0), fmt_ci(r2['efficiency'], ci2['efficiency'], 1.0)]
                            for m, r, r2, ci, ci2 in res.values()]))
            L.append('\nPaired differences (sliding ISS minus Picasso, same frames, same resampled blocks), 95 % intervals; the offset-corrected rung is the like-for-like comparison:\n')
            rows_d = []
            if 'ISS-0' in res:
                ra, rb = res['ISS-0'], res['R0']
                d = C.bootstrap_diff_ci(ra[1]['per_frame'], rb[1]['per_frame'], B=boot, seed=0, block=BLOCK)
                d2 = C.bootstrap_diff_ci(ra[2]['per_frame'], rb[2]['per_frame'], B=boot, seed=0, block=BLOCK)
                rows_d.append(['a priori default ISS vs Picasso defaults (R0)', fmt_ci(ra[1]['jaccard'] - rb[1]['jaccard'], d['jaccard']),
                               fmt_ci(ra[1]['rmse_nm'] - rb[1]['rmse_nm'], d['rmse_nm'], 1.0), fmt_ci(ra[1]['efficiency'] - rb[1]['efficiency'], d['efficiency'], 1.0),
                               fmt_ci(ra[2]['efficiency'] - rb[2]['efficiency'], d2['efficiency'], 1.0)])
                # Picasso's default threshold is in camera units (net gradient 5000 DN) and is far above the signal of the low-photon sets: the defaults' gap is partly a scale effect, so the untuned ISS
                # is also compared with Picasso whose threshold alone is tuned (R1E+c) and with the fully tuned, offset-corrected rung (R2E+c)
                for other, olabel in (('R1E+c', 'Picasso with only its net gradient tuned for E, offset corrected (R1E+c)'), ('R2E+c', 'Picasso with every option tuned for E, offset corrected (R2E+c)')):
                    if other not in res:
                        continue
                    rb = res[other]
                    d = C.bootstrap_diff_ci(ra[1]['per_frame'], rb[1]['per_frame'], B=boot, seed=0, block=BLOCK)
                    d2 = C.bootstrap_diff_ci(ra[2]['per_frame'], rb[2]['per_frame'], B=boot, seed=0, block=BLOCK)
                    rows_d.append(['a priori default ISS vs %s' % olabel, fmt_ci(ra[1]['jaccard'] - rb[1]['jaccard'], d['jaccard']), fmt_ci(ra[1]['rmse_nm'] - rb[1]['rmse_nm'], d['rmse_nm'], 1.0),
                                   fmt_ci(ra[1]['efficiency'] - rb[1]['efficiency'], d['efficiency'], 1.0), fmt_ci(ra[2]['efficiency'] - rb[2]['efficiency'], d2['efficiency'], 1.0)])
            if 'ISS-L' in res:
                for other, olabel in (('R0', 'Picasso defaults (R0)'), ('ISS-0', 'the a priori default with kappa 5'), ('ISS-B', 'the BIC-stop default with the labelled width')):
                    if other not in res:
                        continue
                    ra, rb = res['ISS-L'], res[other]
                    d = C.bootstrap_diff_ci(ra[1]['per_frame'], rb[1]['per_frame'], B=boot, seed=0, block=BLOCK)
                    d2 = C.bootstrap_diff_ci(ra[2]['per_frame'], rb[2]['per_frame'], B=boot, seed=0, block=BLOCK)
                    rows_d.append(['label-free ISS vs %s' % olabel, fmt_ci(ra[1]['jaccard'] - rb[1]['jaccard'], d['jaccard']), fmt_ci(ra[1]['rmse_nm'] - rb[1]['rmse_nm'], d['rmse_nm'], 1.0),
                                   fmt_ci(ra[1]['efficiency'] - rb[1]['efficiency'], d['efficiency'], 1.0), fmt_ci(ra[2]['efficiency'] - rb[2]['efficiency'], d2['efficiency'], 1.0)])
            if 'ISS-B' in res:
                for other, olabel in (('R0', 'Picasso defaults (R0)'), ('ISS-0', 'the a priori default with kappa 5'), ('ISS-E', 'the E-tuned ISS')):
                    if other not in res:
                        continue
                    ra, rb = res['ISS-B'], res[other]
                    d = C.bootstrap_diff_ci(ra[1]['per_frame'], rb[1]['per_frame'], B=boot, seed=0, block=BLOCK)
                    d2 = C.bootstrap_diff_ci(ra[2]['per_frame'], rb[2]['per_frame'], B=boot, seed=0, block=BLOCK)
                    rows_d.append(['BIC-stop ISS vs %s' % olabel, fmt_ci(ra[1]['jaccard'] - rb[1]['jaccard'], d['jaccard']), fmt_ci(ra[1]['rmse_nm'] - rb[1]['rmse_nm'], d['rmse_nm'], 1.0),
                                   fmt_ci(ra[1]['efficiency'] - rb[1]['efficiency'], d['efficiency'], 1.0), fmt_ci(ra[2]['efficiency'] - rb[2]['efficiency'], d2['efficiency'], 1.0)])
            for w in OPS:
                for prow in LADDER_VS[w]:
                    ra, rb = res['ISS-' + w], res[prow]
                    d = C.bootstrap_diff_ci(ra[1]['per_frame'], rb[1]['per_frame'], B=boot, seed=0, block=BLOCK)
                    d2 = C.bootstrap_diff_ci(ra[2]['per_frame'], rb[2]['per_frame'], B=boot, seed=0, block=BLOCK)
                    rows_d.append(['%s-tuned ISS vs %s' % (w, rb[0].replace('Picasso, ', 'Picasso ')),
                                   fmt_ci(ra[1]['jaccard'] - rb[1]['jaccard'], d['jaccard']), fmt_ci(ra[1]['rmse_nm'] - rb[1]['rmse_nm'], d['rmse_nm'], 1.0),
                                   fmt_ci(ra[1]['efficiency'] - rb[1]['efficiency'], d['efficiency'], 1.0),
                                   fmt_ci(ra[2]['efficiency'] - rb[2]['efficiency'], d2['efficiency'], 1.0)])
                    if prow.endswith('+c'):
                        summary.append((name, w, ra[1]['jaccard'] - rb[1]['jaccard'], d['jaccard'], ra[1]['rmse_nm'] - rb[1]['rmse_nm'], d['rmse_nm'],
                                        ra[1]['efficiency'] - rb[1]['efficiency'], d['efficiency'], ra[2]['efficiency'] - rb[2]['efficiency'], d2['efficiency']))
            L.append(table(['comparison', 'dJ % (all)', 'dRMSE nm (all)', 'dE (all)', 'dE (challenge)'], rows_d))
            if 'ISS-0' in res:
                tuning.append((name, res['R0'][1]['efficiency'], res['R2E+c'][1]['efficiency'], res['ISS-0'][1]['efficiency'], res['ISS-E'][1]['efficiency'],
                               res['R0'][1]['jaccard'], res['R2J+c'][1]['jaccard'], res['ISS-0'][1]['jaccard'], res['ISS-J'][1]['jaccard']))
            sec = raw[tkey + '%s|secs' % ops['E']['cfg']]
            psec = ('; Picasso with its defaults: %.4f (%d frames)' % (raw['P|%s|secs' % name][0] / raw['P|%s|secs' % name][1], raw['P|%s|secs' % name][1])
                    if 'P|%s|secs' % name in raw else '')
            L.append('\nSeconds per frame of the E-tuned configuration (one core, shared machine): mean %.2f, median %.2f%s.\n' % (sec.mean(), np.median(sec), psec))
            L += report_variants(raw, name, gt, pmin, boot, ops, res['ISS-E'])
    if summary:
        L += summary_lines(summary)
    if tuning:
        L += tuning_lines(tuning)
    return '\n'.join(L) + '\n'


def tuning_lines(tuning):
    """What tuning on the development frames buys (point estimates on the test frames): the rise of E (and of J) from the defaults to the tuned, offset-corrected configuration."""
    L = ['## What tuning buys (test frames, point estimates)', '',
         'The rise from the untuned to the tuned configuration: Picasso from its defaults (R0) to the option-tuned, offset-corrected rung (R2E+c for E, R2J+c for J); the sliding ISS from its '
         'a priori default to the point tuned on the development frames (with the offset removed).', '']
    rows = []
    for name, e0, e2, ei0, ei1, j0, j2, ji0, ji1 in tuning:
        rows.append([name, '%.1f' % e0, '%.1f' % e2, '%+.1f' % (e2 - e0), '%.1f' % ei0, '%.1f' % ei1, '%+.1f' % (ei1 - ei0), '%.1f' % (100 * j0), '%.1f' % (100 * j2), '%+.1f' % (100 * (j2 - j0)),
                     '%.1f' % (100 * ji0), '%.1f' % (100 * ji1), '%+.1f' % (100 * (ji1 - ji0))])
    L.append(table(['set', 'E Picasso R0', 'E Picasso tuned', 'rise', 'E ISS default', 'E ISS tuned', 'rise', 'J Picasso R0', 'J Picasso tuned', 'rise', 'J ISS default', 'J ISS tuned', 'rise'], rows))
    L.append('')
    L.append('Mean rise of E over the sets: Picasso %+.1f, sliding ISS %+.1f; of J: Picasso %+.1f, sliding ISS %+.1f.' % (
        np.mean([t[2] - t[1] for t in tuning]), np.mean([t[4] - t[3] for t in tuning]), 100 * np.mean([t[6] - t[5] for t in tuning]), 100 * np.mean([t[8] - t[7] for t in tuning])))
    L.append('')
    return L


def bic_dev_lines(raw, name, gt):
    """The BIC stop against the significance levels on the development frames, for the same configuration, no width bound, no offset removed: E, J and the detections per frame."""
    key = 'D|%s|%s|' % (name, DEFAULT['cfg'])
    d0, d1 = (int(v) for v in raw[key + 'frames'])
    idx = np.arange(d0, d1)
    rows = []
    for kap in KAPPAS:
        sc = C.score(unpack(raw[key + '%g|off' % kap], raw[key + '%g|flat' % kap]), gt, idx)
        rows.append(['significance %g' % kap, '%.1f' % (100 * sc['jaccard']), '%.1f' % sc['rmse_nm'], '%.1f' % sc['efficiency'], '%.3f' % (sc['n_det'] / len(idx))])
    kb = 'B|%s|dev|bic4|' % name
    sc = C.score(unpack(raw[kb + 'off'], raw[kb + 'flat']), gt, idx)
    rows.append(['**BIC (no threshold)**', '%.1f' % (100 * sc['jaccard']), '%.1f' % sc['rmse_nm'], '%.1f' % sc['efficiency'], '%.3f' % (sc['n_det'] / len(idx))])
    return ['### The BIC stop against the significance levels on the development frames (%s, no width bound, offset not removed; frames %d-%d)\n' % (DEFAULT['cfg'], d0, d1 - 1),
            table(['stop', 'J %', 'RMSE nm', 'efficiency', 'detections / frame'], rows), '']


def wtl(d, ci, higher_better=True):
    """'win' if the paired interval of the difference excludes zero in the sliding ISS's favour, 'loss' if it excludes zero against it, else 'tie'."""
    lo, hi = (ci[0], ci[1]) if higher_better else (-ci[1], -ci[0])
    return 'win' if lo > 0 else ('loss' if hi < 0 else 'tie')


def summary_lines(summary):
    """The record of the sliding ISS against the offset-corrected, option-tuned Picasso (R2E+c for the E-tuned point, R2J+c for the J-tuned one) over the sets."""
    L = ['## Summary over the sets: sliding ISS against the tuned, offset-corrected Picasso', '',
         'Each cell: the paired difference (ISS minus Picasso) with its 95 % interval and the verdict (win: the interval excludes 0 in favour of the sliding ISS; '
         'loss: against it; tie: it contains 0). The E-tuned point is compared with R2E+c, the J-tuned point with R2J+c, both in the primary setting (all emitters) '
         'unless noted.', '']
    rows, tally = [], {}
    for name, w, dj, cij, dr, cir, de, cie, dec, ciec in summary:
        rows.append([name, w + '-tuned', '%s %s' % (fmt_ci(dj, cij), wtl(dj, cij)), '%s %s' % (fmt_ci(dr, cir, 1.0), wtl(dr, cir, False)),
                     '%s %s' % (fmt_ci(de, cie, 1.0), wtl(de, cie)), '%s %s' % (fmt_ci(dec, ciec, 1.0), wtl(dec, ciec))])
        for k, v in (('J', wtl(dj, cij)), ('RMSE', wtl(dr, cir, False)), ('E', wtl(de, cie)), ('E (challenge)', wtl(dec, ciec))):
            tally.setdefault((w, k), []).append(v)
    L.append(table(['set', 'operating point', 'dJ %', 'dRMSE nm (lower is better)', 'dE', 'dE (challenge)'], rows))
    L += ['', 'Record (wins / ties / losses over the %d sets):' % len({r[0] for r in summary}), '']
    L.append(table(['operating point', 'metric', 'wins', 'ties', 'losses'],
                   [[w + '-tuned', k, v.count('win'), v.count('tie'), v.count('loss')] for (w, k), v in sorted(tally.items())]))
    L.append('')
    return L


def _score_pair(dets, gt, idx, pmin):
    return C.score(dets, gt, idx), C.score(dets, gt, idx, border_nm=450.0, photon_min=pmin)


def _variant_op(raw, key, name, gt, frames_key='dev'):
    """For a variant stored as `key`: the (kappa, s_max) of the largest development efficiency with the dev offset removed; returns (kappa, smax, off, dev E)."""
    da, db = (int(x) for x in raw[key + 'dev|frames'])
    idx = np.arange(da, db)
    best = None
    for k in KAPPAS:
        base = unpack(raw[key + 'dev|snr%g|off' % k], raw[key + 'dev|snr%g|flat' % k])
        for smax in SMAX:
            dets = post(base, smax)
            off = offset_of(dets, gt, idx)
            e = C.score(shift(dets, off), gt, idx)['efficiency']
            cand = (round(e, 6), k, -smax)
            if best is None or cand > best[0]:
                best = (cand, k, smax, off, e)
    return best[1], best[2], best[3], best[4]


def report_variants(raw, name, gt, pmin, boot, ops, ref):
    """The design variants (part V) and the stop rules (part S) of one set, as table lines. `ref` = (label, score, score challenge, ci, ci challenge) of
    the main E-tuned result on the same test frames."""
    L = []
    vs = [v for v in VARIANTS if 'V|%s|%s|dev|frames' % (name, v) in raw]
    if vs:
        L.append('### Design variants of the E-tuned method (%s; threshold, width bound and offset re-chosen on the development frames for each variant)\n' % ops['E']['cfg'])
        va, vb = (int(x) for x in raw['V|%s|%s|test|frames' % (name, vs[0])])
        L.append('Test frames %d-%d; the differences are against the main result scored on the same frames.\n' % (va, vb - 1))
        rows = []
        for v in vs:
            key = 'V|%s|%s|' % (name, v)
            ta, tb = (int(x) for x in raw[key + 'test|frames'])
            kv, sv, offv, ev = _variant_op(raw, key, name, gt)
            dets = shift(post(unpack(raw[key + 'test|snr%g|off' % kv], raw[key + 'test|snr%g|flat' % kv]), sv), offv)
            sc, sc2 = _score_pair(dets, gt, np.arange(ta, tb), pmin)
            ra, rb = (int(x) for x in raw['T|%s|frames' % name])
            if (ta, tb) == (ra, rb):
                ref1, ref2 = ref[1], ref[2]
            else:                                                       # the variant ran on the first frames of the test range only: the main result is rescored on those
                ref1, ref2 = _score_pair(op_test_dets(raw, name, ops['E'])[ta - ra:tb - ra], gt, np.arange(ta, tb), pmin)
            dd = C.bootstrap_diff_ci(sc['per_frame'], ref1['per_frame'], B=boot, seed=0, block=BLOCK)
            dd2 = C.bootstrap_diff_ci(sc2['per_frame'], ref2['per_frame'], B=boot, seed=0, block=BLOCK)
            rows.append([v, '%g' % kv, '%g' % sv, '%.1f' % (100 * sc['jaccard']), '%.1f' % (100 * sc['precision']), '%.1f' % (100 * sc['recall']), '%.1f' % sc['rmse_nm'],
                         '%.1f' % sc['efficiency'], fmt_ci(sc['efficiency'] - ref1['efficiency'], dd['efficiency'], 1.0), fmt_ci(sc2['efficiency'] - ref2['efficiency'], dd2['efficiency'], 1.0),
                         '%.2f' % raw[key + 'test|secs'].mean()])
        L.append(table(['variant', 'kappa (dev)', 'width <= (dev)', 'J %', 'precision %', 'recall %', 'RMSE nm', 'efficiency', 'dE vs the main result (all)', 'dE (challenge)', 's / frame'], rows))
        L.append('\n`full` repeats the main E-tuned result through a second code path (same frames, same options, run separately): its row must equal the main result.\n')
    skey = 'S|%s|' % name
    if skey + 'frames' in raw:
        sa, sb = (int(v) for v in raw[skey + 'frames'])
        sidx = np.arange(sa, sb)
        L.append('### Stop rules on the whole path (test frames %d-%d, the path run to %d entries without an early stop)\n' % (sa, sb - 1, 25 if name.endswith('LD') else 40))
        kap = float(raw[skey + 'kappa'][0])
        off, smax = raw[skey + 'off'], float(raw[skey + 'smax'][0])
        names = [('significance, kappa %g' % kap, 'snr%g' % kap)] + [('discrepancy, z = %g' % z, 'disc%g' % z) for z in (2.0, 3.0)] + [
            ('largest relative gap', 'gap')] + [('gap gated by the first entrant (kappa %g)' % g, 'gapg%g' % g) for g in (5.0, 6.0)]
        res = []
        for label, tag in names:
            dets = shift(post(unpack(raw[skey + 'rules|%s|off' % tag], raw[skey + 'rules|%s|flat' % tag]), smax), off)
            sc, sc2 = _score_pair(dets, gt, sidx, pmin)
            res.append((label, sc, sc2))
        rows = []
        for label, sc, sc2 in res:
            dd = C.bootstrap_diff_ci(sc['per_frame'], res[0][1]['per_frame'], B=boot, seed=0, block=BLOCK)
            rows.append([label, '%.1f' % (100 * sc['jaccard']), '%.1f' % (100 * sc['precision']), '%.1f' % (100 * sc['recall']), '%.1f' % sc['rmse_nm'], '%.1f' % sc['efficiency'],
                         fmt_ci(sc['efficiency'] - res[0][1]['efficiency'], dd['efficiency'], 1.0), '%.2f' % (sc['n_det'] / len(sidx))])
        L.append(table(['rule', 'J %', 'precision %', 'recall %', 'RMSE nm', 'efficiency', 'dE vs significance (paired)', 'detections / frame'], rows))
        L.append('')
    return L


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--only', choices=['C', 'D', 'D2', 'T', 'T0', 'B', 'G', 'L', 'P', 'V', 'S'])
    ap.add_argument('--variants', nargs='*', default=list(VARIANTS))
    ap.add_argument('--sets', nargs='*', default=list(SETS))
    ap.add_argument('--configs', nargs='*', default=list(CONFIGS))
    ap.add_argument('--frames', type=int, default=None, help='limit the frames of each range (smoke runs)')
    ap.add_argument('--from-raw', action='store_true')
    ap.add_argument('--boot', type=int, default=BOOT)
    ap.add_argument('--picasso-python', default=os.environ.get('SMLM_PICASSO_PYTHON'))
    A = ap.parse_args()
    if A.from_raw:
        raw = load_raw()
        with open(MD, 'w') as f:
            f.write(report(raw, A.boot, A.sets))
        print('wrote', MD)
        return
    if A.only == 'C':
        part_C(A.sets)
    elif A.only == 'D':
        part_D(A.sets, A.configs, A.frames)
    elif A.only == 'D2':
        part_D2(A.sets, A.frames)
    elif A.only == 'T':
        part_T(A.sets, A.frames)
    elif A.only == 'T0':
        part_T0(A.sets, A.frames)
    elif A.only == 'B':
        part_B(A.sets, A.frames)
    elif A.only == 'G':
        part_G(A.sets)
    elif A.only == 'L':
        part_L(A.sets, A.frames)
    elif A.only == 'V':
        part_V(A.sets, A.variants, A.frames)
    elif A.only == 'S':
        part_S(A.sets, A.frames)
    elif A.only == 'P':
        part_P(A.sets, A.frames, python=A.picasso_python)
    else:
        ap.error('give --only PART or --from-raw')


if __name__ == '__main__':
    main()
