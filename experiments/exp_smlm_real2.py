"""
Second real-data study of the sliding inverse scale space against Picasso on nuclear-pore recordings (pre-registered, fresh frames)
==================================================================================================================================

Implements `results/smlm/real2_plan.md` (declared 2026-11-09 14:32 UTC) and writes `results/smlm/real2.md` from `results/smlm/real2_raw.npz`. Read the plan first: it is
the specification, this docstring only says where each piece of it lives. The functions of `exp_smlm_real_gt.py` (pore finder, pore model, radial-profile fit, J_s / E_s,
Picasso post-filters, the coordinate ascent, the bootstrap pieces) are imported unchanged; the first study's raw (`real_gt_raw.npz`, development frames, frozen Picasso
rows, calibrations) is READ, never written.

Test sets (fresh frames, plan section 1): T1 = D1 frames 15001-17000; T2a = D2 field 19 frames 1-1055; T2b = D2 field 17 frames 301-1054 in a window disjoint from the
development window. Parts, in order (each stores its results in the raw):
    ref   T_ref pores, windows, validity gates V1 / V4 / V3, expected yields (tables only; no raw frame, no localizer).
    reg   registration shift curves from the RAW frames against the published table of the same frames; gate R1.
    cal   label-free calibration of every test set (and of its merged frames, Part B).
    dev   the ISS4 development paths of 8 configurations per dataset, the four-dial coordinate ascent for E_s (Part A) and, for Part B, the Picasso and ISS post-filter
          ascents for J_B on merged development frames; freeze and hash.
    test  ONE run per test set (Picasso R0 / R1E / R2E and the ISS4 paths; Part B merged frames).
    diag  development-only diagnostics of the grid-edge optimum (plan section 0).
    (--from-raw regenerates results/smlm/real2.md from the raw alone, byte for byte.)
"""
import argparse
import hashlib
import json
import os
import sys
import time

os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import numpy as np                                                   # noqa: E402
from scipy import ndimage as ndi                                     # noqa: E402
from scipy.optimize import linear_sum_assignment                     # noqa: E402
from scipy.spatial import cKDTree                                    # noqa: E402

HERE = os.path.abspath(os.path.dirname(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
sys.path.insert(0, REPO)
import smlm_common as C                                              # noqa: E402
import exp_smlm_real_gt as G                                         # noqa: E402

OUT = os.environ.get('REAL2_OUT', os.path.join(REPO, 'results', 'smlm'))
RAW = os.path.join(OUT, 'real2_raw.npz')
MD = os.path.join(OUT, 'real2.md')
OLD_RAW = os.environ.get('REAL2_OLD_RAW', os.path.join(REPO, 'results', 'smlm', 'real_gt_raw.npz'))
DATA = G.DATA

# ==================================================================================================================================
# the specification of the plan (constants fixed in the plan before any test frame was touched)
# ==================================================================================================================================
BLOCK = G.BLOCK                      # frames per resampled block (10)
N_BOOT = G.N_BOOT                    # 1000
BOOT_SEED, BOOT_SEED_PORES = G.BOOT_SEED, G.BOOT_SEED_PORES
MARGIN_E = G.MARGIN_E                # 3 E points
MARGIN_J = 0.03                      # 3 J points (part B)
TESTS = {
    'T1': dict(ds='D1', field=None, frames=(15001, 17000), ref_ranges=[(2001, 10000), (12001, 15000), (17001, 20000)], group='D1'),
    'T2a': dict(ds='D2', field=19, frames=(1, 1055), ref_ranges=[(1057, 30000)], group='D2'),
    'T2b': dict(ds='D2', field=17, frames=(301, 1054), ref_ranges=[(1057, 30000)], group='D2'),
}
GROUPS = {'D1': ['T1'], 'D2': ['T2a', 'T2b']}
KAPPAS2 = (4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0, 20.0, 30.0, 50.0)
FW_GRID = (0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0)
SMIN_GRID = (None, 1.0, 1.02, 1.05, 1.1, 1.15, 1.2, 1.3)
SMAX_GRID = (1.2, 1.3, 1.4, 1.6, 1.8, 2.0, 2.5, 3.2)
ISS_START = dict(kappa=5.0, fw=0.0, smin=None, smax=3.2)
ISS_KEYS = ('kappa', 'fw', 'smin', 'smax')
ISS_ROUNDS = 3
CFG_SETS = {'D1': ('hp3', 'hp3h', 'hp5', 'hp5h', 'hp3h_s12', 'hp3h_s17', 'hp3h_r15', 'hp3h_r25'),
            'D2': ('hp3', 'hp3h', 'hp5', 'hp5h', 'hp5_s12', 'hp5_s17', 'hp5_r15', 'hp5_r25')}
SIGMA_RENDER_PX = {'D1': 1.6, 'D2': 1.33}
REG_BLOCK = 125                      # frames per block of the shift curve
REG_DEG = 3
REG_SEARCH_PX = 8
REG_HP_PX, REG_SMOOTH_PX = 6.0, 1.0
R1_MAX_DEV_NM, R1_MIN_Q = 20.0, 2.5
R2_BLOCK = 250
R2_MED_NM, R2_P90_NM = 15.0, 30.0
CAP = 400
# part B
K_MERGE, GAP_MERGE, B_BLOCK = 3, 7, 21
RM_REL = 0.5                         # matching radius in units of sigma0
CLOSE_REL = 2.0                      # close pairs: closer than 2 sigma0
BORDER_NM = 200.0
BOOT_BLOCK_B = 10
SPECIES = {'D1': 'zenodo_5068525', 'D2': 'zenodo_3540443'}


def fmt(x, nd=1):
    return G.fmt(x, nd)


jdump, jload = G.jdump, G.jload


def load_raw(path=None):
    return G.load_raw(path or RAW)


def save_raw(raw, path=None):
    G.save_raw(raw, path or RAW)


def old_raw():
    return G.load_raw(OLD_RAW)


# ==================================================================================================================================
# data access
# ==================================================================================================================================
def load_tab(ds, field=None):
    """The published table with the photon column: dict(frame (1-based), x, y (nm), ph)."""
    sp = G.DS[ds]
    path = os.path.join(DATA, sp['table'] % field if '%d' in sp['table'] else sp['table'])
    a = np.genfromtxt(path, delimiter=',', skip_header=sp['tab_skip'], usecols=[sp['tab_cols'][k] for k in ('frame', 'x', 'y')] + [5])
    a = a[np.isfinite(a[:, :3]).all(axis=1)]
    return dict(frame=a[:, 0].astype(np.int64), x=a[:, 1], y=a[:, 2], ph=a[:, 3])


def frames_path(ds, field=None):
    sp = G.DS[ds]
    if ds == 'D1':
        return os.path.join(DATA, 'zenodo_5068525', 'NPCData_glass_STORM.tif')
    return os.path.join(DATA, sp['prefix'] % field)


def read_stack(ds, field, a, b):
    """Frames a..b (1-based, inclusive) of a recording, uint16 (n, ny, nx), the whole field of view."""
    import tifffile
    with tifffile.TiffFile(frames_path(ds, field)) as t:
        return t.asarray(key=range(a - 1, b))


def crop(fr, win, n=None):
    if win is None:
        return fr
    n = n or G.DS['D2']['window_px']
    return fr[:, win[1]:win[1] + n, win[0]:win[0] + n]


def window_origin_nm(ds, win):
    return (0.0, 0.0) if win is None else (win[0] * G.DS[ds]['pixel_nm'], win[1] * G.DS[ds]['pixel_nm'])


def test_window(raw, name):
    return tuple(int(v) for v in raw['ref|%s|window' % name]) if TESTS[name]['ds'] == 'D2' else None


def n_frames_of(name):
    a, b = TESTS[name]['frames']
    return b - a + 1


def load_test_frames(raw, name):
    t = TESTS[name]
    a, b = t['frames']
    return crop(read_stack(t['ds'], t['field'], a, b), test_window(raw, name)).astype(np.float32)


# ==================================================================================================================================
# the pore model of a test set (the first study's Scorer, built from explicit arguments)
# ==================================================================================================================================
class Scorer2(G.Scorer):
    """`exp_smlm_real_gt.Scorer` for a set whose pores, ring radius, region and convention shift are given explicitly (same assignment, radial bins and areas)."""

    def __init__(self, centres, R, region, sp, conv, x_min=None, phases=None):
        self.set_name, self.ds = None, 'D1' if sp is G.DS['D1'] else 'D2'
        self.sp = sp
        self.conv = np.asarray(conv, float)
        cen = np.asarray(centres, float)
        self.allc = cen + self.conv
        self.R, self.region = float(R), tuple(region)
        sc = G.disc_inside(cen, self.R, sp['disc_margin'], self.region)
        if x_min is not None:
            sc &= cen[:, 0] > x_min
        sc &= G.disc_inside(self.allc, self.R, sp['disc_margin'], self.region)
        self.scored = sc
        self.phase_all = None if phases is None else np.asarray(phases, float)
        self.rmax = self.R + sp['disc_margin']
        self.nb = int(np.ceil(self.rmax / G.RADIAL_BIN_NM))
        self.edges = np.arange(self.nb + 1) * G.RADIAL_BIN_NM
        self.tree = cKDTree(self.allc)
        pitch = 10.0
        gx = np.arange(self.region[0] + pitch / 2, self.region[2], pitch)
        gy = np.arange(self.region[1] + pitch / 2, self.region[3], pitch)
        Gr = np.stack(np.meshgrid(gx, gy), axis=-1).reshape(-1, 2)
        r, j = self.tree.query(Gr, distance_upper_bound=self.rmax)
        ok = np.isfinite(r) & (j < len(self.allc))
        ok[ok] &= sc[j[ok]]
        b = np.minimum((r[ok] / G.RADIAL_BIN_NM).astype(int), self.nb - 1)
        self.area = np.bincount(b, minlength=self.nb) * pitch * pitch
        self.scored_idx = np.nonzero(sc)[0]
        pos = np.full(len(self.allc), -1, np.int64)
        pos[self.scored_idx] = np.arange(len(self.scored_idx))
        self.pos = pos
        ap = np.zeros((len(self.scored_idx), self.nb))
        np.add.at(ap, (pos[j[ok]], b), pitch * pitch)
        self.area_pore = ap
        self.n_pores = int(sc.sum())


def convention2(old, ds):
    """The first study's coordinate convention of the table against the localizers, measured on its development frames (D1: none, D2: half a pixel)."""
    return G.convention(old, ds)


def region_of(raw, name):
    t = TESTS[name]
    sp = G.DS[t['ds']]
    if t['ds'] == 'D1':
        return (0.0, 0.0, sp['shape'][1] * sp['pixel_nm'], sp['shape'][0] * sp['pixel_nm'])
    w = test_window(raw, name)
    n = sp['window_px']
    return (w[0] * sp['pixel_nm'], w[1] * sp['pixel_nm'], (w[0] + n) * sp['pixel_nm'], (w[1] + n) * sp['pixel_nm'])


def common_R(raw, name):
    """Ring radius of a test set: its own R_ref (D1), the mean of the two D2 fields' R_ref (so that the pooled radial bins agree)."""
    g = jload(raw['ref|gates'])
    if TESTS[name]['ds'] == 'D1':
        return g['T1']['R_ref']
    return float(np.mean([g['T2a']['R_ref'], g['T2b']['R_ref']]))


def scorer_of(raw, name, old=None):
    old = old if old is not None else old_raw()
    t = TESTS[name]
    sp = G.DS[t['ds']]
    return Scorer2(raw['ref|%s|centres' % name], common_R(raw, name), region_of(raw, name), sp, convention2(old, t['ds']),
                   x_min=sp['x_test_min'] if t['ds'] == 'D1' else None, phases=raw['ref|%s|phase' % name][:, 0] if sp['eight_fold'] else None)


def expected_yield2(raw, name, sc):
    """n_exp = the T_ref ring yield (the same pores) scaled to the number of scored frames."""
    t = TESTS[name]
    tab = load_tab(t['ds'], t['field'])
    ref = G.table_select(tab, t['ref_ranges'])
    r, j, keep = sc.assign(np.c_[ref['x'], ref['y']] + sc.conv)
    f = G.fit_profile(sc.hist(r, keep)[0], sc.area, sc.edges, sc.R)
    n_ref = sum(b - a + 1 for a, b in t['ref_ranges'])
    return f['N'] * n_frames_of(name) / n_ref, f


# ==================================================================================================================================
# step ref: pores, windows, gates (tables only)
# ==================================================================================================================================
def pick_window(cen, R, sp, exclude=None, step=10):
    """The 128 x 128 px window (origin in px, step 10) holding the most pores whose disc is inside it; windows overlapping `exclude` (an origin) are not allowed."""
    wpx = sp['window_px']
    best, bw = -1, None
    for x0 in range(0, 419 - wpx + 1, step):
        for y0 in range(0, 419 - wpx + 1, step):
            if exclude is not None and abs(x0 - exclude[0]) < wpx and abs(y0 - exclude[1]) < wpx:
                continue
            reg = (x0 * sp['pixel_nm'], y0 * sp['pixel_nm'], (x0 + wpx) * sp['pixel_nm'], (y0 + wpx) * sp['pixel_nm'])
            n = int(G.disc_inside(cen, R, sp['disc_margin'], reg).sum())
            if n > best:
                best, bw = n, (x0, y0)
    return bw, best


def step_ref(raw, old, log=print):
    gates = {}
    dev_win = tuple(int(v) for v in old['window|D2dev'])
    for name, t in TESTS.items():
        sp = G.DS[t['ds']]
        tab = load_tab(t['ds'], t['field'])
        ref = G.table_select(tab, t['ref_ranges'])
        region = (0.0, 0.0, sp['shape'][1] * sp['pixel_nm'], sp['shape'][0] * sp['pixel_nm']) if t['ds'] == 'D1' else (0.0, 0.0, sp['hist_size_nm'], sp['hist_size_nm'])
        t0 = time.time()
        cen, rad, cnt = G.find_pores(ref['x'], ref['y'], region, sp)
        R_ref = float(np.median(rad))
        g = dict(n_ref=int(len(ref['x'])), n_pores=int(len(cen)), R_ref=R_ref)
        for par, nm in ((1, 'odd'), (0, 'even')):
            h = G.block_split(ref, t['ref_ranges'], par)
            g['R_' + nm] = G.pore_radius(h['x'], h['y'], cen, sp)
        g['V1_abs_diff'] = abs(g['R_odd'] - g['R_even'])
        raw['ref|%s|centres' % name] = cen
        raw['ref|%s|radii' % name] = rad
        raw['ref|%s|counts' % name] = cnt
        if sp['eight_fold']:
            ph = G.pore_phases(ref['x'], ref['y'], cen, R_ref, sp)
            g['frac_8fold_significant'] = float(np.mean(ph[:, 2] < 0.05))
            raw['ref|%s|phase' % name] = ph
        gates[name] = g
        log('%s: T_ref %d localizations; %d accepted pores; R_ref %.2f nm (%.0f s)' % (name, g['n_ref'], len(cen), R_ref, time.time() - t0))
    for name, t in TESTS.items():
        if t['ds'] != 'D2':
            continue
        sp = G.DS['D2']
        excl = dev_win if (t['field'] == G.DS['D2']['fields']['dev']) else None
        bw, best = pick_window(raw['ref|%s|centres' % name], gates[name]['R_ref'], sp, exclude=excl)
        raw['ref|%s|window' % name] = np.array(bw, dtype=np.int64)
        gates[name]['window_px'] = list(bw)
        gates[name]['pores_in_window'] = best
        log('%s: window origin %s px holds %d accepted pores' % (name, bw, best))
    gates['V2'] = dict(diff=abs(gates['T2a']['R_ref'] - gates['T2b']['R_ref']))
    gates['V2']['pass'] = bool(gates['V2']['diff'] <= 3.0)
    raw['ref|gates'] = jdump(gates)
    # scored pores, expected yields
    for name in TESTS:
        sc = scorer_of(raw, name, old)
        n_exp, fref = expected_yield2(raw, name, sc)
        gates[name]['scored_pores'] = sc.n_pores
        gates[name]['n_exp'] = n_exp
        gates[name]['ref_fit'] = fref
        log('%s: %d scored pores, n_exp %.0f ring localizations (T_ref fit N %.0f mu %.1f w %.1f)' % (name, sc.n_pores, n_exp, fref['N'], fref['mu'], fref['w']))
    raw['ref|gates'] = jdump(gates)
    raw['ref|sha256'] = np.array(ref_hash(raw))
    return gates


def ref_hash(raw):
    h = hashlib.sha256()
    for k in sorted(raw):
        if k.startswith('ref|') and k != 'ref|sha256':
            h.update(k.encode())
            h.update(np.ascontiguousarray(raw[k]).tobytes())
    return h.hexdigest()


# ==================================================================================================================================
# step reg: registration shift curves from the raw frames (plan section 2)
# ==================================================================================================================================
def render_table(x, y, shape, pix, sigma_px):
    """The table localizations (nm, corner convention) on the pixel grid of the frames (pixel i covers [i pix, (i + 1) pix)), Gaussian-blurred."""
    ny, nx = shape
    H, _, _ = np.histogram2d(y, x, bins=[np.arange(ny + 1) * pix, np.arange(nx + 1) * pix])
    return ndi.gaussian_filter(H, sigma_px)


def sum_image(fr):
    """The positive part of the frames minus their temporal median, summed (the blinking emitters of the block)."""
    fr = np.asarray(fr, np.float32)
    return np.clip(fr - np.median(fr, axis=0)[None], 0.0, None).sum(0)


def xcorr_shift(A, B, pix, search=REG_SEARCH_PX, hp=REG_HP_PX, sm=REG_SMOOTH_PX):
    """The translation (nm) by which the image A is displaced against B (A(x) = B(x - d)): the peak of their circular cross-correlation after a Gaussian high-pass of hp px and
    a smoothing of sm px, within +-search px, refined by a parabola; returns (dx, dy, peak over the standard deviation of the correlation)."""
    A = ndi.gaussian_filter(A - ndi.gaussian_filter(A, hp), sm)
    B = ndi.gaussian_filter(B - ndi.gaussian_filter(B, hp), sm)
    F = np.fft.fftshift(np.fft.ifft2(np.fft.fft2(A) * np.conj(np.fft.fft2(B))).real)
    cy, cx = F.shape[0] // 2, F.shape[1] // 2
    sub = F[cy - search:cy + search + 1, cx - search:cx + search + 1]
    iy, ix = np.unravel_index(np.argmax(sub), sub.shape)

    def par(v, i):
        if 0 < i < len(v) - 1:
            d = v[i - 1] - 2 * v[i] + v[i + 1]
            return i + 0.5 * (v[i - 1] - v[i + 1]) / d if d != 0 else float(i)
        return float(i)
    py, px = par(sub[:, ix], iy), par(sub[iy, :], ix)
    return float((px - search) * pix), float((py - search) * pix), float(sub.max() / max(F.std(), 1e-12))


def reg_blocks(n):
    """Near-equal blocks of about REG_BLOCK frames: list of (first, stop) in 0-based positions of the set."""
    nb = max(1, int(round(n / REG_BLOCK)))
    cuts = np.linspace(0, n, nb + 1).astype(int)
    return [(int(cuts[i]), int(cuts[i + 1])) for i in range(nb)]


def fit_shift_curve(blocks, shifts, n):
    """Cubic (at most; fewer parameters than blocks minus one) polynomial in t = (frame + 0.5) / n per axis through the block shifts, one 2.5-sigma clip. Returns
    (coef (2, deg + 1), residual (nb, 2))."""
    t = np.array([(0.5 * (a + b)) / n for a, b in blocks])
    deg = int(min(REG_DEG, max(len(blocks) - 2, 0)))
    coef = np.zeros((2, REG_DEG + 1))
    res = np.zeros((len(blocks), 2))
    for ax in range(2):
        y = shifts[:, ax]
        keep = np.ones(len(y), bool)
        for _ in range(2):
            c = np.polyfit(t[keep], y[keep], deg)
            r = y - np.polyval(c, t)
            s = max(np.std(r[keep]), 1e-9)
            nk = np.abs(r) <= 2.5 * s
            if nk.sum() < deg + 2 or (nk == keep).all():
                break
            keep = nk
        c = np.polyfit(t[keep], y[keep], deg)
        coef[ax, REG_DEG - deg:] = c
        res[:, ax] = y - np.polyval(coef[ax], t)
    return coef, res


def shift_of(raw, name, n=None):
    """The per-frame translation (n, 2) nm of the raw frames against the table's coordinates (raw = table + s)."""
    n = n or n_frames_of(name)
    coef = raw['reg|%s|coef' % name]
    t = (np.arange(n) + 0.5) / n
    return np.c_[np.polyval(coef[0], t), np.polyval(coef[1], t)]


def step_reg(raw, old, log=print):
    for name, t in TESTS.items():
        ds, sp = t['ds'], G.DS[t['ds']]
        a, b = t['frames']
        n = b - a + 1
        fr = read_stack(ds, t['field'], a, b)
        tab = G.table_select(load_tab(ds, t['field']), [(a, b)])
        conv = convention2(old, ds)
        ny, nx = fr.shape[1:]
        rows = []
        blocks = reg_blocks(n)
        for lo, hi in blocks:
            A = sum_image(fr[lo:hi])
            m = (tab['frame'] >= a + lo) & (tab['frame'] < a + hi)
            B = render_table(tab['x'][m] + conv[0], tab['y'][m] + conv[1], (ny, nx), sp['pixel_nm'], SIGMA_RENDER_PX[ds])
            rows.append(xcorr_shift(A, B, sp['pixel_nm']))
        rows = np.array(rows)
        coef, res = fit_shift_curve(blocks, rows[:, :2], n)
        raw['reg|%s|blocks' % name] = np.c_[np.array(blocks, float), rows]
        raw['reg|%s|coef' % name] = coef
        raw['reg|%s|resid' % name] = res
        log('%s: %d blocks; shift (x, y) nm from (%.0f, %.0f) to (%.0f, %.0f); max residual %.1f nm; min quality %.2f' % (
            name, len(blocks), rows[0, 0], rows[0, 1], rows[-1, 0], rows[-1, 1], np.abs(res).max(), rows[:, 2].min()))
        del fr
    return raw


def gate_r1(raw, name):
    res = raw['reg|%s|resid' % name]
    q = raw['reg|%s|blocks' % name][:, 4]
    return dict(max_dev=float(np.abs(res).max()), min_q=float(q.min()), passed=bool(np.abs(res).max() <= R1_MAX_DEV_NM and q.min() >= R1_MIN_Q))


# ==================================================================================================================================
# step cal: label-free calibration of every test set (and of its merged frames)
# ==================================================================================================================================
def calibrate_store(raw, key, frames, ds, log=print):
    import exp_smlm_real as REAL
    r0 = C.run_picasso(frames, G.picasso_params(ds))
    cal = REAL.calibrate(frames, r0)
    law, sig = cal['law'], cal['sig']
    raw['cal|%s|bg' % key] = cal['bg']
    raw['cal|%s|law' % key] = np.array([law['a'], law['b0'], law['x_ref'], law['c'], law['scale'], law['a_hf'], law['b0_hf'], law['v_dark'], law['x_dark']])
    raw['cal|%s|sigma0' % key] = np.array(json.dumps(sig))
    raw['cal|%s|r0_per_frame' % key] = np.array(sum(len(x) for x in r0) / len(frames))
    log('%s: %d frames of %d x %d px; Picasso R0 %d localizations; sigma0 %.3f px (%d selected%s); noise law a %.3f b0 %.1f x_ref %.1f (scale %.3f)' % (
        key, len(frames), frames.shape[1], frames.shape[2], sum(len(x) for x in r0), sig['sigma0'], sig['n_sel'], ', FALLBACK' if sig['fallback'] else '', law['a'], law['b0'],
        law['x_ref'], law['scale']))
    return r0


def step_cal(raw, old, log=print):
    for name in TESTS:
        fr = load_test_frames(raw, name)
        calibrate_store(raw, name, fr, TESTS[name]['ds'], log)
    return raw


# ==================================================================================================================================
# the ISS paths (kappa read at ten values) and the four post-filters
# ==================================================================================================================================
def iss_frame2(frame, ctx, order, kappas):
    """`exp_smlm_real_gt.iss_frame` with a given list of kappas: returns (dict kappa -> (k, 4) x_nm, y_nm, w_adu, sigma_px), entries, seconds)."""
    from src.core import imaging as IM
    from src.core import imaging_guard as IG
    cfg, P = ctx['cfg'], ctx['P']
    ny, nx = frame.shape
    t0 = time.perf_counter()
    y = np.asarray(frame, dtype=np.float64) - ctx['x_ref']
    delta = float(np.median(y - ctx['bg']))
    bgv = ctx['bg'] + delta
    Pp = ctx.get('Pp')
    if Pp is None:
        Pp = ctx['Pp'] = IM.PixelPSF(ny, nx, cfg['sigma0'], parts=[tuple(p) for p in cfg['parts']])
    a, b0 = ctx['a'], ctx['b0']

    def sw_fn(c, w):
        m = Pp.synth(c, w) if len(w) else 0.0
        return 1.0 / np.sqrt(a * (bgv + np.maximum(m, 0.0)) + b0)
    y0 = y - bgv
    yh = P.highpass(y0 - np.median(y0))
    sw0 = sw_fn(np.zeros((0, 3)), np.zeros(0))
    R = IG.sliding_path_guarded(P, yh, sw0, ctx['L'], max_entries=cfg['cap'], sw_fn=sw_fn, sw_on='params', sig_bounds=tuple(cfg['widths']), tol=cfg['tol'], passes=1,
                                rtrunc=ctx['rtrunc'], r_ex=cfg['r_ex'], fit_radius=cfg['fit_radius'], order=order, snr_min=G.SNR_MIN_PATH, snr_extra=G.SNR_EXTRA_PATH, cycle_guard=True)
    out = {}
    for kap in kappas:
        i = IM.stop_snr(R['hist'], R['taus'], R['snr'], kap)[0]
        h = R['hist'][i]
        c, w = h['c'], h['w']
        pix = ctx['pixel_nm']
        out[kap] = np.c_[(c[:, 0] + 0.5) * pix, (c[:, 1] + 0.5) * pix, w, c[:, 2]] if len(w) else np.zeros((0, 4))
    return out, len(R['taus']), time.perf_counter() - t0


def iss_run2(frames, cal, cfg, pixel_nm, kappas=KAPPAS2, order='dual', log=print, tag=''):
    """The ISS path on every frame: (dets[kappa] = list of per-frame arrays, seconds per frame, entries per frame)."""
    import exp_smlm_real as REAL
    ctx = REAL.make_context(cal['bg'], cal['a'], cal['b0'], cal['x_ref'], cfg, frames.shape[1:])
    ctx['pixel_nm'] = pixel_nm
    dets = {k: [] for k in kappas}
    secs, ent = [], []
    t_start = time.time()
    for i, fr in enumerate(frames):
        d, n, sec = iss_frame2(fr, ctx, order, kappas)
        for k in kappas:
            dets[k].append(d[k])
        secs.append(sec)
        ent.append(n)
        if (i + 1) % 100 == 0:
            log('    %s %s: %d / %d frames, %.2f s per frame, %.0f s' % (tag, cfg['cfg'], i + 1, len(frames), float(np.mean(secs)), time.time() - t_start))
    return dets, np.array(secs), np.array(ent)


def pack_dets(dets):
    """dets[kappa] (per-frame arrays) -> flat arrays for the raw: {kappa: (frame (m,), data (m, 4))}."""
    out = {}
    for k, lst in dets.items():
        fr = np.concatenate([np.full(len(d), i, np.int32) for i, d in enumerate(lst)]) if len(lst) else np.zeros(0, np.int32)
        dat = np.concatenate([np.asarray(d, np.float64).reshape(len(d), 4) for d in lst if len(d)]) if any(len(d) for d in lst) else np.zeros((0, 4))
        out[k] = (fr, dat)
    return out


def put_paths(raw, key, dets, n_frames, secs, ent):
    for k, (fr, dat) in pack_dets(dets).items():
        raw['%s|k%g|frame' % (key, k)] = fr
        raw['%s|k%g|dat' % (key, k)] = dat
    raw['%s|n_frames' % key] = np.array(n_frames)
    raw['%s|secs' % key] = np.asarray(secs, float)
    raw['%s|entries' % key] = np.asarray(ent)


def get_paths(raw, key, kappas=KAPPAS2):
    n = int(raw['%s|n_frames' % key])
    dets = {}
    for k in kappas:
        fr, dat = raw['%s|k%g|frame' % (key, k)], raw['%s|k%g|dat' % (key, k)]
        order = np.searchsorted(fr, np.arange(n + 1))
        dets[k] = [dat[order[i]:order[i + 1]] for i in range(n)]
    return dets


def w_ref_of(dets):
    """The label-free amplitude scale: the median fitted amplitude of the kappa-10 detections of the frames."""
    w = np.concatenate([d[:, 2] for d in dets[10.0] if len(d)]) if any(len(d) for d in dets[10.0]) else np.zeros(1)
    return float(np.median(w))


def filter_dets(dets_list, f, w_ref, s_scale):
    """The ISS post-filters: amplitude >= fw x w_ref, sigma >= smin x s (if set), sigma <= smax x s; f = dict(fw, smin, smax)."""
    lo = -np.inf if f.get('smin') is None else f['smin'] * s_scale
    hi = f['smax'] * s_scale + 1e-9
    wmin = f['fw'] * w_ref
    out = []
    for d in dets_list:
        if len(d):
            m = (d[:, 2] >= wmin) & (d[:, 3] >= lo) & (d[:, 3] <= hi)
            out.append(d[m])
        else:
            out.append(d)
    return out


def ascent_iss(ev, grids, start=None, max_rounds=ISS_ROUNDS):
    """The coordinate ascent of the plan (section 3) for an arbitrary scorer ev(point) -> value: strict improvement, the dials in the order ISS_KEYS, at most max_rounds rounds."""
    cur = dict(start or ISS_START)
    best = ev(cur)
    n_ev = 1
    rounds = 0
    for rnd in range(max_rounds):
        rounds = rnd + 1
        changed = False
        for k in ISS_KEYS:
            for v in grids[k]:
                if v == cur[k]:
                    continue
                q = dict(cur)
                q[k] = v
                x = ev(q)
                n_ev += 1
                if x > best + 1e-12:
                    best, cur, changed = x, q, True
        if not changed:
            break
    return cur, best, rounds, n_ev


ISS_GRIDS = dict(kappa=KAPPAS2, fw=FW_GRID, smin=SMIN_GRID, smax=SMAX_GRID)


def choose_cfg(results):
    """The configuration of the largest development value (rounded to 6 places); ties: the earlier configuration. results: list of (name, value, point)."""
    best = None
    for r in results:
        if best is None or round(r[1], 6) > round(best[1], 6):
            best = r
    return best


# ==================================================================================================================================
# part B: merged frames, reference events, matching
# ==================================================================================================================================
def make_groups(n, K=K_MERGE, gap=GAP_MERGE):
    """Groups of K frames `gap` apart inside each block of K * gap consecutive positions: (g, K) integer array of 0-based positions."""
    blk = K * gap
    out = []
    for b0 in range(0, n - blk + 1, blk):
        for i in range(gap):
            out.append([b0 + i + gap * m for m in range(K)])
    return np.array(out, dtype=np.int64).reshape(-1, K)


def merge_frames(frames, groups, bg=None):
    """Merged frames: the sum of the K frames of a group minus (K - 1) times the per-pixel temporal median of the whole set."""
    frames = np.asarray(frames, np.float32)
    bg = np.median(frames, axis=0) if bg is None else bg
    K = groups.shape[1]
    return np.stack([frames[g].sum(0) - (K - 1) * bg for g in groups]).astype(np.float32)


def reference_events(tab, a, groups, origin_nm_, conv, shifts, region_nm, ph_bright, border=BORDER_NM):
    """Per group the table's events of its constituent frames, in the raw WINDOW coordinates (table + convention + the frame's shift - origin): arrays (m, 5) = x, y, photons,
    constituent index (0..K-1), bright (1.0 / 0.0). Events within `border` nm of the window edge are dropped. tab frames are 1-based, a = the set's first frame."""
    W = region_nm[2] - region_nm[0]
    H = region_nm[3] - region_nm[1]
    fr0 = tab['frame'] - a
    order = np.argsort(fr0, kind='stable')
    fs = fr0[order]
    out = []
    for g in groups:
        rows = []
        for m, f in enumerate(g):
            lo, hi = np.searchsorted(fs, f), np.searchsorted(fs, f + 1)
            idx = order[lo:hi]
            x = tab['x'][idx] + conv[0] + shifts[f, 0] - origin_nm_[0]
            y = tab['y'][idx] + conv[1] + shifts[f, 1] - origin_nm_[1]
            ok = (x >= border) & (x < W - border) & (y >= border) & (y < H - border)
            ph = tab['ph'][idx][ok]
            rows.append(np.c_[x[ok], y[ok], ph, np.full(ok.sum(), m, float), (ph >= ph_bright).astype(float)])
        out.append(np.concatenate(rows) if rows else np.zeros((0, 5)))
    return out


def match_group(D, R, rm):
    """One to one matching of detections D (n, >=2) to reference events R (m, >=2) closer than rm (the Hungarian assignment on the distances): returns (matched ref mask (m,),
    matched det mask (n,), index of the matched reference event of every detection (-1 if none))."""
    md = np.zeros(len(D), bool)
    mr = np.zeros(len(R), bool)
    which = np.full(len(D), -1, np.int64)
    if len(D) == 0 or len(R) == 0:
        return mr, md, which
    cost = np.hypot(D[:, None, 0] - R[None, :, 0], D[:, None, 1] - R[None, :, 1])
    big = np.where(cost <= rm, cost, 1e6)
    ri, ci = linear_sum_assignment(big)
    ok = big[ri, ci] < 1e5
    md[ri[ok]] = True
    mr[ci[ok]] = True
    which[ri[ok]] = ci[ok]
    return mr, md, which


def close_flags(R, close_nm):
    """For the bright reference events: whether another bright event of a different constituent frame lies within close_nm (and the number of such pairs)."""
    flag = np.zeros(len(R), bool)
    br = np.nonzero(R[:, 4] > 0.5)[0]
    pairs = []
    if len(br) > 1:
        tr = cKDTree(R[br, :2])
        for a_, b_ in tr.query_pairs(close_nm):
            i, j = br[a_], br[b_]
            if R[i, 3] != R[j, 3]:
                flag[i] = flag[j] = True
                pairs.append((i, j))
    return flag, pairs


def group_counts(D, R, rm, close_nm):
    """The count vector of one merged frame: [TP, FP, FN, close events, close matched, close pairs, pairs with both matched, detections]. A detection matched to a bright
    reference event is a TP, one matched to a dim one is neither TP nor FP, one matched to none is a FP, a bright event matched to none is a FN."""
    mr, md, which = match_group(D, R, rm)
    bright = R[:, 4] > 0.5 if len(R) else np.zeros(0, bool)
    tp = int(bright[which[which >= 0]].sum())
    fp = int((~md).sum())
    fn = int((bright & ~mr).sum())
    flag, pairs = close_flags(R, close_nm) if len(R) else (np.zeros(0, bool), [])
    both = sum(1 for i, j in pairs if mr[i] and mr[j])
    return np.array([tp, fp, fn, int(flag.sum()), int((flag & mr).sum()), len(pairs), both, len(D)], float)


def jb_of(counts):
    c = counts.sum(0)
    return float(c[0] / max(c[0] + c[1] + c[2], 1.0))


# ==================================================================================================================================
# step dev, part A: the ISS4 development paths and the four-dial ascent for E_s
# ==================================================================================================================================
def dev_frames_of(ds):
    a, b = G.DS[ds]['dev_frames']
    return b - a + 1


def dev_scorer(old, ds):
    """(Scorer, n_exp, origin) of a development set of the first study."""
    dev = ds + 'dev'
    sc = G.Scorer(old, dev)
    n_exp, _ = G.expected_yield(old, sc, dev_frames_of(ds))
    return sc, n_exp, G.origin_nm(old, dev)


def make_ev_E(dets, sc, n_exp, org, s_scale):
    """The E_s of an ISS operating point (kappa, fw, smin, smax) on a development path, cached; -inf for fewer than 30 detections."""
    w_ref = w_ref_of(dets)
    cache = {}

    def ev(p):
        key = tuple(p[k] for k in ISS_KEYS)
        if key not in cache:
            dl = filter_dets(dets[p['kappa']], p, w_ref, s_scale)
            if sum(len(d) for d in dl) < 30:
                cache[key] = -np.inf
            else:
                cache[key] = G.score_loc(sc, G.Loc.from_dets(dl, org), 0, n_exp)['E']
        return cache[key]
    return ev, w_ref


def step_dev_A(raw, old, ds, log=print, save=None):
    dev = ds + 'dev'
    sc, n_exp, org = dev_scorer(old, ds)
    cal, sigma0 = G.iss_cal(old, dev)
    pix = G.DS[ds]['pixel_nm']
    s_scale = sigma0 / 1.4
    frames = None
    res = []
    for name in CFG_SETS[ds]:
        key = 'devpath|%s|%s' % (ds, name)
        if key + '|n_frames' not in raw:
            if frames is None:
                frames = G.load_frames(dev, G.window_of(old, dev))
            t0 = time.time()
            dets, secs, ent = iss_run2(frames, cal, G.grid_cfg(name, sigma0), pix, log=log, tag=dev)
            put_paths(raw, key, dets, len(frames), secs, ent)
            log('%s dev path %s: %d frames, %.2f s per frame (%.0f s)' % (dev, name, len(frames), secs.mean(), time.time() - t0))
            if save:
                save(raw)
        dets = get_paths(raw, key)
        ev, w_ref = make_ev_E(dets, sc, n_exp, org, s_scale)
        cur, best, rounds, n_ev = ascent_iss(ev, ISS_GRIDS)
        log('%s dev ascent %s: E_s %.2f at %s (%d rounds, %d points scored)' % (dev, name, best, cur, rounds, n_ev))
        res.append((name, best, cur, rounds, n_ev, w_ref))
    name, best, cur, rounds, n_ev, w_ref = choose_cfg(res)
    dets = get_paths(raw, 'devpath|%s|%s' % (ds, name))
    ev, w_ref = make_ev_E(dets, sc, n_exp, org, s_scale)
    dl = filter_dets(dets[cur['kappa']], cur, w_ref, s_scale)
    r = G.score_loc(sc, G.Loc.from_dets(dl, org), 0, n_exp)
    raw['dev|%s|iss4' % ds] = jdump(dict(cfg=name, kappa=cur['kappa'], fw=cur['fw'], smin=cur['smin'], smax=cur['smax'], E=best, J=r['J'], N=r['fit']['N'], floor=r['fit']['floor'],
                                          w=r['fit']['w'], n_exp=n_exp, n_frames=dev_frames_of(ds), sigma0=sigma0, w_ref=w_ref, n_points_scored=int(sum(x[4] for x in res)),
                                          per_cfg={x[0]: dict(E=x[1], point=x[2], rounds=x[3], n_points=x[4]) for x in res}))
    log('%s dev ISS4 frozen: %s' % (ds, jload(raw['dev|%s|iss4' % ds])))
    return raw


# ==================================================================================================================================
# step dev, part B: merged development frames, the Picasso and ISS post-filter ascents for J_B
# ==================================================================================================================================
def pic_mask(arr, f):
    """Picasso's post-filters (strict inequalities) on a stacked (m, 6) detection array."""
    m = np.ones(len(arr), bool)
    if f.get('min_net_gradient') is not None:
        m &= arr[:, C.NG_COL] > f['min_net_gradient']
    if f.get('min_photons') is not None:
        m &= arr[:, C.PH_COL] > f['min_photons']
    sx, sy = np.abs(arr[:, C.SX_COL]), np.abs(arr[:, C.SY_COL])
    if f.get('min_sigma_nm') is not None:
        m &= (sx > f['min_sigma_nm']) & (sy > f['min_sigma_nm'])
    if f.get('max_sigma_nm') is not None:
        m &= (sx < f['max_sigma_nm']) & (sy < f['max_sigma_nm'])
    return m


def stack_dets(dets):
    """per-frame (n_i, k) arrays -> (group index (m,), array (m, k)) with finite positions."""
    fr = np.concatenate([np.full(len(d), i, np.int64) for i, d in enumerate(dets)]) if len(dets) else np.zeros(0, np.int64)
    k = max((d.shape[1] for d in dets if len(d)), default=2)
    arr = np.concatenate([np.asarray(d, float).reshape(len(d), -1) for d in dets if len(d)]) if any(len(d) for d in dets) else np.zeros((0, k))
    ok = np.isfinite(arr[:, :2]).all(axis=1) if len(arr) else np.zeros(0, bool)
    return fr[ok], arr[ok]


def split_by_group(fr, arr, n_groups):
    order = np.searchsorted(fr, np.arange(n_groups + 1))
    return [arr[order[i]:order[i + 1]] for i in range(n_groups)]


def counts_for(det_lists, refs, rm, close_nm):
    return np.array([group_counts(D, R, rm, close_nm) for D, R in zip(det_lists, refs)])


def b_constants(old, ds, sigma0_single):
    pix = G.DS[ds]['pixel_nm']
    return RM_REL * sigma0_single * pix, CLOSE_REL * sigma0_single * pix


def step_dev_B(raw, old, ds, log=print, save=None):
    import exp_smlm_picasso as PIC
    dev = ds + 'dev'
    sp = G.DS[ds]
    pix = sp['pixel_nm']
    frames = G.load_frames(dev, G.window_of(old, dev))
    n = len(frames)
    groups = make_groups(n)
    S = merge_frames(frames, groups)
    key = ds + 'devB'
    calibrate_store(raw, key, S, ds, log)
    sc = G.Scorer(old, dev)
    a, b = sp['dev_frames']
    tab = G.table_select(load_tab(ds, sp['fields']['dev'] if ds == 'D2' else None), [(a, b)])
    conv, org = convention2(old, ds), G.origin_nm(old, dev)
    all_ev = reference_events(tab, a, groups, org, conv, np.zeros((n, 2)), sc.region, np.inf)
    ph_bright = float(np.median(np.concatenate([e[:, 2] for e in all_ev])))
    refs = reference_events(tab, a, groups, org, conv, np.zeros((n, 2)), sc.region, ph_bright)
    sigma_single = G.iss_cal(old, dev)[1]
    rm, close_nm = b_constants(old, ds, sigma_single)
    log('%s dev B: %d merged frames, %d reference events per frame (%.1f bright), bright threshold %.0f photons, rm %.0f nm, close %.0f nm' % (
        dev, len(S), np.mean([len(e) for e in refs]), np.mean([e[:, 4].sum() for e in refs]), ph_bright, rm, close_nm))
    out = dict(ph_bright=ph_bright, rm=rm, close_nm=close_nm, n_groups=len(S), sigma_single=sigma_single)
    # Picasso: the post-filters re-ascended for J_B for two structures (the one of R2E and the defaults), the better kept
    pic = jload(old['dev|%s|picasso' % ds])
    grids, r1 = G.scaled_grids(pic['fac'])
    base = G.picasso_params(ds)
    best_pic = None
    for sname_, struct in (('R2E', pic['R2E']['struct']), ('default', {})):
        d0 = C.run_picasso(S, dict(base, **struct, min_net_gradient=grids['min_net_gradient'][0]))
        fr, arr = stack_dets(d0)
        cache = {}

        def evP(f, fr=fr, arr=arr, cache=cache):
            k = tuple(f.get(x) for x in grids)
            if k not in cache:
                m = pic_mask(arr, f)
                cache[k] = jb_of(counts_for(split_by_group(fr[m], arr[m], len(S)), refs, rm, close_nm))
            return cache[k]
        f_best, j_best, rounds = G.ascent(evP, grids, PIC.MAX_ROUNDS)
        log('%s dev B Picasso structure %s: J_B %.3f with %s (%d points)' % (dev, sname_, j_best, f_best, len(cache)))
        if best_pic is None or round(j_best, 6) > round(best_pic['J'], 6):
            best_pic = dict(structure=sname_, filters=f_best, J=j_best, rounds=rounds, n_points=len(cache), struct=struct, params=dict(base, **struct, **f_best))
        else:
            best_pic['n_points'] += len(cache)
    out['R2B'] = best_pic
    f_best, j_best = best_pic['filters'], best_pic['J']
    j0 = jb_of(counts_for(split_by_group(*stack_dets(C.run_picasso(S, base)), len(S)), refs, rm, close_nm))
    out['R0'] = dict(J=j0)
    log('%s dev B Picasso: R0 J_B %.3f; R2B J_B %.3f (%s structure) with %s (%d points)' % (dev, j0, j_best, best_pic['structure'], f_best, best_pic['n_points']))
    # ISS: the configuration of ISS4, the four dials re-ascended for J_B on the merged frames
    iss4 = jload(raw['dev|%s|iss4' % ds])
    cal, sigma0B = G.iss_cal(raw, key)
    kp = 'devpathB|%s' % ds
    if kp + '|n_frames' not in raw:
        dets, secs, ent = iss_run2(S, cal, G.grid_cfg(iss4['cfg'], sigma0B), pix, log=log, tag=dev + 'B')
        put_paths(raw, kp, dets, len(S), secs, ent)
        if save:
            save(raw)
    dets = get_paths(raw, kp)
    w_ref = w_ref_of(dets)
    s_scale = sigma0B / 1.4
    cacheI = {}

    def evI(p):
        k = tuple(p[x] for x in ISS_KEYS)
        if k not in cacheI:
            dl = filter_dets(dets[p['kappa']], p, w_ref, s_scale)
            cacheI[k] = jb_of(counts_for(dl, refs, rm, close_nm))
        return cacheI[k]
    cur, best, rounds, n_ev = ascent_iss(evI, ISS_GRIDS)
    out['ISSB'] = dict(cfg=iss4['cfg'], kappa=cur['kappa'], fw=cur['fw'], smin=cur['smin'], smax=cur['smax'], J=best, rounds=rounds, n_points=n_ev, w_ref=w_ref, sigma0=sigma0B)
    log('%s dev B ISS: J_B %.3f at %s (%d points)' % (dev, best, cur, n_ev))
    raw['dev|%s|B' % ds] = jdump(out)
    return raw


def freeze_hash2(raw):
    h = hashlib.sha256()
    merged = tuple('cal|%sB|' % n for n in TESTS)                    # the label-free calibrations of the merged TEST frames are made during the test run
    for k in sorted(raw):
        if k.startswith(('ref|', 'reg|', 'cal|', 'dev|')) and not k.startswith(merged):
            h.update(k.encode())
            h.update(np.ascontiguousarray(raw[k]).tobytes())
    return h.hexdigest()


def step_dev(raw, old, log=print, save=None):
    for ds in ('D1', 'D2'):
        if 'dev|%s|iss4' % ds not in raw:
            step_dev_A(raw, old, ds, log, save)
            if save:
                save(raw)
        if 'dev|%s|B' % ds not in raw:
            step_dev_B(raw, old, ds, log, save)
            if save:
                save(raw)
    raw['freeze|time_utc'] = np.array(time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()))
    raw['freeze|sha256'] = np.array(freeze_hash2(raw))
    if save:
        save(raw)
    log('frozen %s sha256 %s' % (str(raw['freeze|time_utc']), str(raw['freeze|sha256'])))
    return raw


# ==================================================================================================================================
# step test: ONE run per test set (Part A rows and Part B merged frames)
# ==================================================================================================================================
def table_loc2(raw, old, name):
    """The published table restricted to the frames of the set as a Loc (frames 0-based within the set, the convention applied)."""
    t = TESTS[name]
    a, b = t['frames']
    tab = G.table_select(load_tab(t['ds'], t['field']), [(a, b)])
    conv = convention2(old, t['ds'])
    return G.Loc(tab['frame'] - a, np.c_[tab['x'], tab['y']] + conv)


def step_test_set(raw, old, name, log=print):
    t = TESTS[name]
    ds = t['ds']
    pix = G.DS[ds]['pixel_nm']
    win = test_window(raw, name)
    org = window_origin_nm(ds, win)
    frames = load_test_frames(raw, name)
    n = len(frames)
    pic = jload(old['dev|%s|picasso' % ds])
    iss4 = jload(raw['dev|%s|iss4' % ds])
    devB = jload(raw['dev|%s|B' % ds])
    raw['test|%s|n_frames' % name] = np.array(n)
    for row in ('R0', 'R1E', 'R2E'):
        info = {}
        t0 = time.time()
        dets = C.run_picasso(frames, pic[row]['params'], info=info)
        G.put_loc(raw, 'test|%s|%s' % (name, row), G.Loc.from_dets(dets, org), secs=[info['run_s'] / n])
        log('%s test %s: %d localizations (%.0f s)' % (name, row, sum(len(d) for d in dets), time.time() - t0))
    cal, sigma0 = G.iss_cal(raw, name)
    t0 = time.time()
    dets, secs, ent = iss_run2(frames, cal, G.grid_cfg(iss4['cfg'], sigma0), pix, log=log, tag=name)
    put_paths(raw, 'testpath|%s|iss4' % name, dets, n, secs, ent)
    log('%s test ISS4 %s: %d frames, %.2f s per frame (%.0f s)' % (name, iss4['cfg'], n, secs.mean(), time.time() - t0))
    old_dual = jload(old['dev|%s|dual' % ds])
    if (old_dual['cfg'] != iss4['cfg'] or float(old_dual['kappa']) not in KAPPAS2) and ds == 'D2':
        t0 = time.time()
        d2, s2, e2 = iss_run2(frames, cal, G.grid_cfg(old_dual['cfg'], sigma0), pix, kappas=(float(old_dual['kappa']),), log=log, tag=name + 'old')
        put_paths(raw, 'testpath|%s|old' % name, d2, n, s2, e2)
    # part B
    groups = make_groups(n)
    S = merge_frames(frames, groups)
    calibrate_store(raw, name + 'B', S, ds, log)
    for row, params in (('R0', G.picasso_params(ds)), ('R2B', devB['R2B']['params'])):
        info = {}
        dd = C.run_picasso(S, params, info=info)
        G.put_loc(raw, 'testB|%s|%s' % (name, row), G.Loc.from_dets(dd, org), secs=[info['run_s'] / len(S)])
    calB, sigmaB = G.iss_cal(raw, name + 'B')
    t0 = time.time()
    dB, sB, eB = iss_run2(S, calB, G.grid_cfg(devB['ISSB']['cfg'], sigmaB), pix, log=log, tag=name + 'B')
    put_paths(raw, 'testpathB|%s' % name, dB, len(S), sB, eB)
    log('%s test B: %d merged frames, ISS %.2f s per frame (%.0f s)' % (name, len(S), sB.mean(), time.time() - t0))
    return raw


def step_test(raw, old, log=print, save=None):
    if 'freeze|sha256' not in raw or str(raw['freeze|sha256']) != freeze_hash2(raw):
        raise RuntimeError('the frozen development choices are missing or changed: the test run is refused')
    raw['test|start_utc'] = np.array(time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()))
    raw['test|meta'] = jdump(dict(python=sys.version.split()[0], numpy=np.__version__, cpu_count=os.cpu_count(), loadavg_at_start=[round(v, 2) for v in os.getloadavg()]))
    if save:
        save(raw)
    for name in TESTS:
        if 'test|%s|R0|xy' % name in raw and 'testpathB|%s|n_frames' % name in raw:
            continue
        step_test_set(raw, old, name, log)
        if save:
            save(raw)
    raw['test|end_utc'] = np.array(time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()))
    if save:
        save(raw)
    return raw


# ==================================================================================================================================
# the rows of a test set (from the raw alone), registration, gates
# ==================================================================================================================================
ROWS_A = ('R0', 'R1E', 'R2E', 'ISS4', 'ISS_OLD', 'PUB')
ROW_NAMES = {'R0': 'Picasso, defaults (R0)', 'R1E': 'Picasso, net gradient tuned for E_s (R1E)', 'R2E': 'Picasso, every option tuned for E_s (R2E)',
             'ISS4': 'sliding ISS, four-dial grid tuned for E_s (ISS4)', 'ISS_OLD': 'sliding ISS, the first study\'s two-dial choice (continuity)',
             'PUB': 'as published by the authors (context only)'}
ROWS_B = ('R0', 'R2B', 'ISSB')
ROW_NAMES_B = {'R0': 'Picasso, defaults (R0)', 'R2B': 'Picasso, post-filters re-tuned for J_B (R2B)', 'ISSB': 'sliding ISS, four dials re-tuned for J_B (ISS-B)'}


def reg_loc(loc, shifts):
    """A Loc translated by -s(frame): raw coordinates -> the table's coordinates."""
    return G.Loc(loc.frame, loc.xy - shifts[loc.frame], loc.extra)


def iss4_loc(raw, name, kappa=None):
    """The ISS4 row of a test set: the frozen filters (amplitude floor relative to the label-free w_ref of the test frames, widths relative to sigma0 / 1.4) at the frozen kappa
    (or at another kappa: the curve rows)."""
    t = TESTS[name]
    iss4 = jload(raw['dev|%s|iss4' % t['ds']])
    sigma0 = G.iss_cal(raw, name)[1]
    dets = get_paths(raw, 'testpath|%s|iss4' % name)
    kap = float(iss4['kappa'] if kappa is None else kappa)
    dl = filter_dets(dets[kap], dict(fw=iss4['fw'], smin=iss4['smin'], smax=iss4['smax']), w_ref_of(dets), sigma0 / 1.4)
    return G.Loc.from_dets(dl, window_origin_nm(t['ds'], test_window(raw, name)))


def iss_old_loc(raw, old, name):
    t = TESTS[name]
    ds = t['ds']
    iss4, od = jload(raw['dev|%s|iss4' % ds]), jload(old['dev|%s|dual' % ds])
    key = 'testpath|%s|iss4' % name if (od['cfg'] == iss4['cfg'] and float(od['kappa']) in KAPPAS2) else 'testpath|%s|old' % name
    if key + '|n_frames' not in raw:
        return None
    sigma0 = G.iss_cal(raw, name)[1]
    kap = float(od['kappa'])
    dets = get_paths(raw, key, kappas=(kap,))
    return G.Loc.from_dets(G.post_width(dets[kap], od['smax_rel'] * sigma0 / 1.4), window_origin_nm(ds, test_window(raw, name)))


def set_rows(raw, old, name, register=True, curve=True):
    """The localization sets of a test set: dict row -> Loc (registered unless register is False), plus the kappa-curve rows 'K<kappa>' of the frozen ISS4 point."""
    locs = {r: G.get_loc(raw, 'test|%s|%s' % (name, r)) for r in ('R0', 'R1E', 'R2E')}
    locs['ISS4'] = iss4_loc(raw, name)
    io = iss_old_loc(raw, old, name)
    if io is not None:
        locs['ISS_OLD'] = io
    if curve:
        for k in KAPPAS2:
            locs['K%g' % k] = iss4_loc(raw, name, k)
    if register:
        sh = shift_of(raw, name)
        locs = {k: reg_loc(v, sh) for k, v in locs.items()}
    locs['PUB'] = table_loc2(raw, old, name)
    return locs


def gate_r2(raw, old, name):
    """Gate R2: the translation that maps each registered row onto the published table, per block of R2_BLOCK frames: median and 90th percentile of its modulus."""
    n = n_frames_of(name)
    sc = scorer_of(raw, name, old)
    locs = set_rows(raw, old, name, register=True, curve=False)
    extent = float(max(sc.region[2], sc.region[3])) + 2000.0
    pub = locs['PUB']
    out = {}
    for row in ('R0', 'R2E', 'ISS4'):
        loc = locs[row]
        mods = []
        for k in range(int(np.ceil(n / R2_BLOCK))):
            m = (loc.frame >= k * R2_BLOCK) & (loc.frame < (k + 1) * R2_BLOCK)
            mp = (pub.frame >= k * R2_BLOCK) & (pub.frame < (k + 1) * R2_BLOCK)
            if m.sum() >= 50 and mp.sum() >= 50:
                dx, dy, _ = G.block_shift(loc.xy[m], pub.xy[mp], extent)
                mods.append(float(np.hypot(dx, dy)))
        med = float(np.median(mods)) if mods else float('nan')
        p90 = float(np.percentile(mods, 90)) if mods else float('nan')
        out[row] = dict(median=med, p90=p90, n_blocks=len(mods), passed=bool(len(mods) > 0 and med <= R2_MED_NM and p90 <= R2_P90_NM))
    return out


def set_valid(raw, old, name, r2=None):
    """The gates of a test set that decide whether it counts for the verdict: V1, V4, R1 and (given) R2."""
    g = jload(raw['ref|gates'])[name]
    r1 = gate_r1(raw, name)
    v = dict(V1=bool(g['V1_abs_diff'] <= 2.0), V4=bool(g['scored_pores'] >= 100), R1=r1['passed'])
    if r2 is not None:
        v['R2'] = all(x['passed'] for x in r2.values())
    v['valid'] = all(v.values())
    return v


# ==================================================================================================================================
# part A: the analysis (profile fits, paired block bootstrap, equal-off-ring yield)
# ==================================================================================================================================
def yield_at_floor(floors, Ns, F):
    """The ring yield of the curve (floor_k, N_k) at the off-ring count F, log-log interpolation, clamped at the ends of the curve."""
    f = np.asarray(floors, float)
    N = np.asarray(Ns, float)
    if len(f) < 2 or not np.isfinite(F) or F <= 0:
        return float('nan')
    o = np.argsort(f, kind='stable')
    return float(np.exp(np.interp(np.log(F), np.log(np.maximum(f[o], 1e-9)), np.log(np.maximum(N[o], 1e-9)))))


def pooled(members, fn):
    return sum(fn(m) for m in members)


def analyse_A(raw, old, gname, n_boot=N_BOOT, log=lambda s: None, register=True, with_curve=True):
    names = GROUPS[gname]
    gates = jload(raw['ref|gates'])
    mem = []
    for name in names:
        sc = scorer_of(raw, name, old)
        sc.sigma0_px = float(jload(raw['cal|%s|sigma0' % name])['sigma0'])
        n = n_frames_of(name)
        nb = int(np.ceil(n / BLOCK))
        locs = set_rows(raw, old, name, register=register, curve=with_curve)
        eight = G.DS[TESTS[name]['ds']]['eight_fold']
        prep = {r: G.prep_row(sc, l, nb, eight) for r, l in locs.items()}
        mem.append(dict(name=name, sc=sc, n=n, nb=nb, n_exp=gates[name]['n_exp'], prep=prep, eight=eight))
    sc0 = mem[0]['sc']
    rows = [r for r in ROWS_A if r in mem[0]['prep']]
    curve = ['K%g' % k for k in KAPPAS2] if with_curve else []
    n_exp = sum(m['n_exp'] for m in mem)
    area = sum(m['sc'].area for m in mem)
    n_pores = sum(m['sc'].n_pores for m in mem)
    n_tot = sum(m['n'] for m in mem)

    def hist(row, idxs=None):
        if idxs is None:
            return sum(m['prep'][row]['H'].sum(0) for m in mem)
        return sum(m['prep'][row]['H'][ix].sum(0) for m, ix in zip(mem, idxs))
    base = {}
    out = dict(group=gname, members=names, n_frames=n_tot, n_pores=n_pores, n_exp=n_exp, R_ref=sc0.R, rows={}, pairs={}, members_info={})
    for r in rows + curve:
        f, J, E = G.stats_of(sc0, hist(r), area, n_exp)
        base[r] = f
        if r in rows:
            cr = np.concatenate([m['prep'][r]['ring_per_pore'] for m in mem])
            d = dict(fit=f, J=J, E=E, n_loc=sum(m['prep'][r]['n_loc'] for m in mem), n_keep=sum(m['prep'][r]['n_keep'] for m in mem),
                     yield_per_pore_kframe=float(f['N'] / n_pores / (n_tot / 1000.0)), frac_pores_8=float(np.mean(cr >= 8)),
                     dup_frac=sum(m['prep'][r]['dup'] for m in mem) / max(sum(m['prep'][r]['n_loc'] for m in mem), 1), mu_bias=f['mu'] - sc0.R)
            if all(m['eight'] for m in mem):
                d['c8'] = float(sum(m['prep'][r]['c8_sum'].sum() for m in mem) / max(sum(m['prep'][r]['c8_n'].sum() for m in mem), 1.0))
                d['c8_null'] = float(np.mean([m['prep'][r]['c8_null'] for m in mem]))
            out['rows'][r] = d
        else:
            out.setdefault('curve', {})[r] = dict(N=f['N'], floor=f['floor'], w=f['w'], E=E, J=J)
    for m in mem:
        out['members_info'][m['name']] = dict(n_frames=m['n'], n_pores=m['sc'].n_pores, n_exp=m['n_exp'], R=m['sc'].R)

    def ratio_of(fits):
        fl = [fits[c][0]['floor'] for c in curve]
        Ns = [fits[c][0]['N'] for c in curve]
        return yield_at_floor(fl, Ns, fits['R2E'][0]['floor']) / max(fits['R2E'][0]['N'], 1e-9)
    fits0 = {r: (base[r], None, None) for r in base}
    out['off_ring_ratio'] = ratio_of(fits0) if with_curve else float('nan')
    pairs = [('ISS4', 'R2E'), ('ISS4', 'R0'), ('R2E', 'R0')] + ([('ISS_OLD', 'R2E'), ('ISS_OLD', 'R0')] if 'ISS_OLD' in rows else [])
    rng = np.random.default_rng(BOOT_SEED)
    rep = {r: dict(E=[], J=[], N=[], w=[]) for r in rows}
    dif = {pr: dict(dE=[], dJ=[], dw=[]) for pr in pairs}
    ratio_rep = []
    t0 = time.time()
    for it in range(n_boot):
        idxs = [rng.integers(0, m['nb'], m['nb']) for m in mem]
        fits = {}
        for r in rows + curve:
            Hs = hist(r, idxs)
            f, J, E = G.stats_of(sc0, Hs, area, n_exp, start=base[r]['x'])
            fits[r] = (f, J, E, Hs)
            if r in rep:
                rep[r]['E'].append(E), rep[r]['J'].append(J), rep[r]['N'].append(f['N']), rep[r]['w'].append(f['w'])
        for pr in pairs:
            a, b = pr
            dif[pr]['dE'].append(fits[a][2] - fits[b][2])
            dif[pr]['dJ'].append(fits[a][1] - fits[b][1])
        wa, wb = G.eq_yield_w(sc0, fits['ISS4'][0], fits['ISS4'][3], fits['R2E'][0], fits['R2E'][3], area, rng)
        dif[('ISS4', 'R2E')]['dw'].append(wa - wb)
        ratio_rep.append(ratio_of(fits) if with_curve else float('nan'))
    log('%s block bootstrap %d replicates: %.0f s' % (gname, n_boot, time.time() - t0))
    for r in rows:
        out['rows'][r]['ci'] = {k: G.percentile_ci(v) for k, v in rep[r].items()}
    for pr in pairs:
        a, b = pr
        d = dict(dE=out['rows'][a]['E'] - out['rows'][b]['E'], dE_ci=G.percentile_ci(dif[pr]['dE']), dJ=out['rows'][a]['J'] - out['rows'][b]['J'], dJ_ci=G.percentile_ci(dif[pr]['dJ']))
        if pr == ('ISS4', 'R2E'):
            wa, wb = G.eq_yield_w(sc0, base[a], hist(a), base[b], hist(b), area, np.random.default_rng(BOOT_SEED + 7))
            d.update(dw=wa - wb, dw_ci=G.percentile_ci(dif[pr]['dw']), w_a=wa, w_b=wb, off_ring_ratio=out['off_ring_ratio'], off_ring_ratio_ci=G.percentile_ci(ratio_rep) if with_curve else [float('nan')] * 2)
        out['pairs']['%s-%s' % pr] = d
    # sensitivity: resample the pores (all their frames); n_exp kept fixed
    rng = np.random.default_rng(BOOT_SEED_PORES)
    Hp = {r: np.concatenate([m['prep'][r]['Hp'] for m in mem]) for r in ('ISS4', 'R2E')}
    Ap = np.concatenate([m['sc'].area_pore for m in mem])
    dE_p = []
    for it in range(n_boot):
        ix = rng.integers(0, n_pores, n_pores)
        ar = Ap[ix].sum(0)
        E_ = {}
        for r in ('ISS4', 'R2E'):
            E_[r] = G.stats_of(sc0, Hp[r][ix].sum(0), ar, n_exp, start=base[r]['x'])[2]
        dE_p.append(E_['ISS4'] - E_['R2E'])
    out['pairs']['ISS4-R2E']['dE_ci_pores'] = G.percentile_ci(dE_p)
    return out


def read_A(res):
    """The reading of the primary outcome of one test dataset under the plan's rule A (section 4)."""
    pr = res['pairs']['ISS4-R2E']
    dE, dEp = pr['dE_ci'], pr['dE_ci_pores']
    better = dE[0] > 0 and dEp[0] > 0
    if better and pr['off_ring_ratio'] <= 1.0:
        return 'tradeoff'
    if better:
        return 'ISS better'
    if dE[0] > 0:
        return 'ISS better in the frame bootstrap only (not in the pore-resampling sensitivity: not claimed)'
    if dE[1] < 0:
        return 'Picasso better'
    if dE[0] > -MARGIN_E and dE[1] < MARGIN_E:
        return 'equivalent'
    return 'inconclusive'


def overall_A(reads, valid):
    """The rule A over the two test datasets (`reads` and `valid` keyed by group name 'D1', 'D2')."""
    vs = [g for g in ('D1', 'D2') if valid.get(g)]
    if len(vs) < 2:
        txt = 'NOT CONCLUDED in registered form: %s void under the plan\'s gates; a real-data advantage needs both datasets.' % ' and '.join(g for g in ('D1', 'D2') if not valid.get(g))
        if vs:
            txt += ' The valid dataset alone: %s.' % reads[vs[0]]
        return txt
    r1, r2 = reads['D1'], reads['D2']
    if r1 == 'ISS better' and r2 == 'ISS better':
        return 'ISS better on BOTH datasets: a real-data advantage may be written, with the caveats of the plan (section 1).'
    if (r1 == 'Picasso better' or r2 == 'Picasso better') and 'ISS better' not in (r1, r2):
        return 'ISS WORSE: Picasso better on at least one dataset and the ISS on neither.'
    if r1 == 'equivalent' and r2 == 'equivalent':
        return 'EQUIVALENT within +-%g E points on both datasets: no real-data advantage.' % MARGIN_E
    return 'MIXED: D1 %s; D2 %s.' % (r1, r2)


# ==================================================================================================================================
# part B: the analysis of the merged frames
# ==================================================================================================================================
def b_member(raw, old, name):
    """Counts per merged frame (n_groups, 8) of the three rows of Part B for a test set, and the bookkeeping."""
    t = TESTS[name]
    ds = t['ds']
    a, b = t['frames']
    n = n_frames_of(name)
    groups = make_groups(n)
    devB = jload(raw['dev|%s|B' % ds])
    win = test_window(raw, name)
    org = window_origin_nm(ds, win)
    tab = G.table_select(load_tab(ds, t['field']), [(a, b)])
    conv = convention2(old, ds)
    refs = reference_events(tab, a, groups, org, conv, shift_of(raw, name), region_of(raw, name), devB['ph_bright'])
    sigma_single = G.iss_cal(raw, name)[1]
    rm, close_nm = b_constants(old, ds, sigma_single)
    ng = len(groups)
    dets = {}
    for r in ('R0', 'R2B'):
        loc = G.get_loc(raw, 'testB|%s|%s' % (name, r))
        order = np.searchsorted(loc.frame, np.arange(ng + 1))
        dets[r] = [loc.xy[order[i]:order[i + 1]] - np.array(org) for i in range(ng)]
    pd = get_paths(raw, 'testpathB|%s' % name)
    sigmaB = G.iss_cal(raw, name + 'B')[1]
    iss = devB['ISSB']
    dl = filter_dets(pd[float(iss['kappa'])], dict(fw=iss['fw'], smin=iss['smin'], smax=iss['smax']), w_ref_of(pd), sigmaB / 1.4)
    dets['ISSB'] = [d[:, :2] for d in dl]
    counts = {r: counts_for(dets[r], refs, rm, close_nm) for r in ROWS_B}
    return dict(name=name, counts=counts, n_groups=ng, rm=rm, close_nm=close_nm, ref_per_frame=float(np.mean([len(e) for e in refs])),
                bright_per_frame=float(np.mean([e[:, 4].sum() for e in refs])))


def b_stats(c):
    s = c.sum(0)
    tp, fp, fn = s[0], s[1], s[2]
    return dict(J=float(tp / max(tp + fp + fn, 1.0)), recall=float(tp / max(tp + fn, 1.0)), fp_per_frame=float(fp / max(len(c), 1)), close_recall=float(s[4] / max(s[3], 1.0)),
                pair_resolution=float(s[6] / max(s[5], 1.0)), n_close=float(s[3]), n_pairs=float(s[5]), det_per_frame=float(s[7] / max(len(c), 1)))


def analyse_B(raw, old, gname, n_boot=N_BOOT, log=lambda s: None):
    mem = [b_member(raw, old, name) for name in GROUPS[gname]]
    out = dict(group=gname, members=GROUPS[gname], rows={}, pairs={})
    allc = {r: np.concatenate([m['counts'][r] for m in mem]) for r in ROWS_B}
    out['n_groups'] = int(sum(m['n_groups'] for m in mem))
    out['info'] = {m['name']: dict(n_groups=m['n_groups'], rm=m['rm'], close_nm=m['close_nm'], ref_per_frame=m['ref_per_frame'], bright_per_frame=m['bright_per_frame']) for m in mem}
    for r in ROWS_B:
        out['rows'][r] = b_stats(allc[r])
    nbs = [int(np.ceil(m['n_groups'] / BOOT_BLOCK_B)) for m in mem]
    pairs = [('ISSB', 'R2B'), ('ISSB', 'R0'), ('R2B', 'R0')]
    rng = np.random.default_rng(BOOT_SEED)
    rep = {r: {k: [] for k in ('J', 'recall', 'close_recall', 'pair_resolution')} for r in ROWS_B}
    dif = {pr: {k: [] for k in ('J', 'close_recall', 'recall', 'pair_resolution')} for pr in pairs}
    for it in range(n_boot):
        sel = {}
        idx = []
        for m, nb in zip(mem, nbs):
            bl = rng.integers(0, nb, nb)
            g = np.concatenate([np.arange(k * BOOT_BLOCK_B, min((k + 1) * BOOT_BLOCK_B, m['n_groups'])) for k in bl])
            idx.append(g)
        st = {}
        for r in ROWS_B:
            c = np.concatenate([m['counts'][r][g] for m, g in zip(mem, idx)])
            st[r] = b_stats(c)
            for k in rep[r]:
                rep[r][k].append(st[r][k])
        for pr in pairs:
            for k in dif[pr]:
                dif[pr][k].append(st[pr[0]][k] - st[pr[1]][k])
    for r in ROWS_B:
        out['rows'][r]['ci'] = {k: G.percentile_ci(v) for k, v in rep[r].items()}
    for pr in pairs:
        a, b = pr
        out['pairs']['%s-%s' % pr] = {k: dict(est=out['rows'][a][k] - out['rows'][b][k], ci=G.percentile_ci(v)) for k, v in dif[pr].items()}
    return out


def read_B(res):
    pr = res['pairs']['ISSB-R2B']
    dJ = pr['J']['ci']
    if dJ[0] > 0 and pr['close_recall']['est'] >= 0:
        return 'ISS better'
    if dJ[0] > 0:
        return 'ISS better in J_B but the close-pair recall difference is negative (not claimed)'
    if dJ[1] < 0:
        return 'Picasso better'
    if dJ[0] > -MARGIN_J and dJ[1] < MARGIN_J:
        return 'equivalent'
    return 'inconclusive'


def overall_B(reads, valid):
    vs = [g for g in ('D1', 'D2') if valid.get(g)]
    if len(vs) < 2:
        txt = 'NOT CONCLUDED in registered form: %s void under the plan\'s gates.' % ' and '.join(g for g in ('D1', 'D2') if not valid.get(g))
        if vs:
            txt += ' The valid dataset alone: %s.' % reads[vs[0]]
        return txt
    r1, r2 = reads['D1'], reads['D2']
    if r1 == 'ISS better' and r2 == 'ISS better':
        return 'ISS better on BOTH datasets on frames made by adding three real frames.'
    if (r1 == 'Picasso better' or r2 == 'Picasso better') and 'ISS better' not in (r1, r2):
        return 'ISS WORSE on merged frames.'
    if r1 == 'equivalent' and r2 == 'equivalent':
        return 'EQUIVALENT within +-%g J points on both datasets.' % MARGIN_J
    return 'MIXED: D1 %s; D2 %s.' % (r1, r2)


# ==================================================================================================================================
# step diag: development-only diagnostics of the grid-edge optimum (plan section 0)
# ==================================================================================================================================
DIAG = {'D1': dict(cfg='hp3h_s12', n=60), 'D2': dict(cfg='hp5_r25', n=100)}
DIAG_KAPPA = 4.0
DIAG_MASK_PX = 6.0


def robust_sd(a):
    a = np.asarray(a, float)
    return float(1.4826 * np.median(np.abs(a - np.median(a))))


def step_diag(raw, old, log=print):
    """On the first frames of the development sets: the robust standard deviation of the standardized high-passed frame and of the standardized residual of the kappa-4 block
    away from every fitted unit (is the label-free noise law too small?), the tails of that residual, and the fraction of the units inside the pore discs that lie on the
    ring, by amplitude fifth (what carries the ring?). Everything is computed on development frames and development pores."""
    import exp_smlm_real as REAL
    from src.core import imaging as IM
    from src.core import imaging_guard as IG
    from scipy.stats import norm
    for ds, d in DIAG.items():
        dev = ds + 'dev'
        sc, n_exp, org = dev_scorer(old, ds)
        frames = G.load_frames(dev, G.window_of(old, dev))[:d['n']]
        cal, sigma0 = G.iss_cal(old, dev)
        cfg = G.grid_cfg(d['cfg'], sigma0)
        ctx = REAL.make_context(cal['bg'], cal['a'], cal['b0'], cal['x_ref'], cfg, frames.shape[1:])
        ctx['pixel_nm'] = G.DS[ds]['pixel_nm']
        P, pix = ctx['P'], ctx['pixel_nm']
        Pp = IM.PixelPSF(*frames.shape[1:], cfg['sigma0'], parts=[tuple(p) for p in cfg['parts']])
        z0, zres, units = [], [], []
        yy, xx = np.mgrid[:frames.shape[1], :frames.shape[2]]
        for i, frame in enumerate(frames):
            y = np.asarray(frame, float) - ctx['x_ref']
            bgv = ctx['bg'] + float(np.median(y - ctx['bg']))
            a_, b0_ = ctx['a'], ctx['b0']

            def sw_fn(c, w):
                m = Pp.synth(c, w) if len(w) else 0.0
                return 1.0 / np.sqrt(a_ * (bgv + np.maximum(m, 0.0)) + b0_)
            y0 = y - bgv
            yh = P.highpass(y0 - np.median(y0))
            sw0 = sw_fn(np.zeros((0, 3)), np.zeros(0))
            R = IG.sliding_path_guarded(P, yh, sw0, ctx['L'], max_entries=cfg['cap'], sw_fn=sw_fn, sw_on='params', sig_bounds=tuple(cfg['widths']), tol=cfg['tol'], passes=1,
                                        rtrunc=ctx['rtrunc'], r_ex=cfg['r_ex'], fit_radius=cfg['fit_radius'], order='dual', snr_min=G.SNR_MIN_PATH, snr_extra=G.SNR_EXTRA_PATH,
                                        cycle_guard=True)
            h = R['hist'][IM.stop_snr(R['hist'], R['taus'], R['snr'], DIAG_KAPPA)[0]]
            c, w = h['c'], h['w']
            mh = P.synth(c, w) if len(w) else np.zeros_like(yh)
            res = (yh - mh) * sw_fn(c, w)
            mask = np.ones(yh.shape, bool)
            for cc in c:
                mask &= ((xx - cc[0]) ** 2 + (yy - cc[1]) ** 2) > DIAG_MASK_PX ** 2
            z0.append((yh * sw0).ravel())
            zres.append(res[mask])
            if len(w):
                units.append(np.c_[(c[:, 0] + 0.5) * pix + org[0], (c[:, 1] + 0.5) * pix + org[1], w])
        z0, zres = np.concatenate(z0), np.concatenate(zres)
        sd = robust_sd(zres)
        out = dict(cfg=d['cfg'], n_frames=len(frames), robust_sd_frame=robust_sd(z0), robust_sd_residual=sd, sd_residual=float(zres.std()))
        for t in (3, 4, 5, 6):
            out['tail_pos_%d' % t] = float(np.mean(zres > t * sd))
            out['tail_neg_%d' % t] = float(np.mean(zres < -t * sd))
            out['tail_gauss_%d' % t] = float(norm.sf(t))
        U = np.concatenate(units)
        r, j, keep = sc.assign(U[:, :2])
        ring = keep & (np.abs(r - sc.R) <= sc.sp['ring_half'])
        q = np.quantile(U[keep, 2], np.linspace(0, 1, 6))
        out['ring_fraction_by_amplitude_fifth'] = [float(ring[keep & (U[:, 2] >= lo) & (U[:, 2] <= hi)].mean()) for lo, hi in zip(q[:-1], q[1:])]
        out['ring_area_fraction'] = float((sc.R + sc.sp['ring_half']) ** 2 - max(sc.R - sc.sp['ring_half'], 0.0) ** 2) / sc.rmax ** 2
        out['units_per_frame'] = len(U) / len(frames)
        out['units_in_discs'] = int(keep.sum())
        # the label-free repair of the law: its variance multiplied by the residual's variance (robust sd squared), the four dials re-tuned on all development frames
        sd2 = sd ** 2
        cal2 = dict(cal, a=cal['a'] * sd2, b0=cal['b0'] * sd2)
        all_frames = G.load_frames(dev, G.window_of(old, dev))
        dets2, _, _ = iss_run2(all_frames, cal2, cfg, pix, log=log, tag=dev + 'rescaled')
        ev2, _ = make_ev_E(dets2, sc, n_exp, org, sigma0 / 1.4)
        cur2, best2, rounds2, n_ev2 = ascent_iss(ev2, ISS_GRIDS)
        stored = jload(raw['dev|%s|iss4' % ds])['per_cfg'].get(d['cfg']) if 'dev|%s|iss4' % ds in raw else None
        out['rescaled'] = dict(variance_factor=sd2, E=best2, point=cur2, n_points=n_ev2, E_stored=None if stored is None else stored['E'],
                               point_stored=None if stored is None else stored['point'])
        raw['diag|%s' % ds] = jdump(out)
        log('%s diag: %s' % (ds, out))
    return raw


# ==================================================================================================================================
# the report (from the raw alone)
# ==================================================================================================================================
md_table = G.md_table


def ci_str(v, ci, nd=1):
    return '%s [%s, %s]' % (fmt(v, nd), fmt(ci[0], nd), fmt(ci[1], nd))


def row_line_A(name, d):
    f = d['fit']
    return [name, d['n_loc'], fmt(f['N'], 0), '%s [%s, %s]' % (fmt(f['w']), fmt(d['ci']['w'][0]), fmt(d['ci']['w'][1])), fmt(f['floor'], 0),
            ci_str(d['J'], d['ci']['J'], 3), ci_str(d['E'], d['ci']['E']), fmt(d['mu_bias']), fmt(d['yield_per_pore_kframe'], 2), fmt(100 * d['dup_frac'], 1)]


def with_void(reading, valid):
    """The reading column of a table: the reading of the rule, marked '(void)' when the dataset's gates failed (the reading is then descriptive, not a verdict)."""
    return reading if valid else reading + ' (void)'


def edge_of(kappa):
    """Where a kappa lies on the extended grid KAPPAS2: the lower edge, the top edge, or the interior."""
    return 'at the lower edge' if kappa == KAPPAS2[0] else 'at the top edge' if kappa == KAPPAS2[-1] else 'in the interior'


def grid_edges_line(edge4, kB):
    """The sentence on the grid edges of the development optima: edge4 = {dataset: number of the 8 configurations of the ISS4 ascent whose optimum is at the top kappa}, kB = {dataset: the ISS-B optimum kappa}."""
    return ('Grid edges of the ISS tuning (development): Part A, the ISS4 optimum is at the top kappa (%g) for %d of 8 configurations of D1 and %d of 8 of D2; Part B, the ISS-B optimum is at kappa %g for D1 (%s of the kappa grid '
            '%g-%g) and at kappa %g for D2 (%s).' % (KAPPAS2[-1], edge4['D1'], edge4['D2'], kB['D1'], edge_of(kB['D1']), KAPPAS2[0], KAPPAS2[-1], kB['D2'], edge_of(kB['D2'])))


def rank_str(vals):
    """'a x > b y > ...' of {label: value}, largest first (ties in the order of the labels)."""
    return ' > '.join('%s %s' % (k, fmt(v, 1)) for k, v in sorted(vals.items(), key=lambda kv: -kv[1]))


def transfer_paragraph(raw, old, resA):
    """The sentences on the transfer of the development tuning to the test frames and on the continuity row, from the stored raw and the Part A results (void or not)."""
    out = []
    parts = []
    for g in ('D1', 'D2'):
        v = resA[g]['pairs']['R2E-R0']
        parts.append('%s %s' % (g, ci_str(v['dE'], v['dE_ci'])))
    out.append('**Tuning on the development frames did not transfer.** On the test frames untuned Picasso (R0) %s Picasso tuned on the development frames (R2E) on both datasets (paired Delta E_s, R2E - R0: %s). The '
               'registered comparison ISS4 - R2E (table above) therefore compares two points that are both %s than untuned Picasso (ISS4 - R0: %s).' %
               ('beats' if all(resA[g]['pairs']['R2E-R0']['dE_ci'][1] < 0 for g in resA) else 'does not beat', '; '.join(parts),
                'worse' if all(resA[g]['pairs']['ISS4-R0']['dE_ci'][1] < 0 for g in resA) else 'not always worse',
                '; '.join('%s %s' % (g, ci_str(resA[g]['pairs']['ISS4-R0']['dE'], resA[g]['pairs']['ISS4-R0']['dE_ci'])) for g in resA)))
    for g in ('D1', 'D2'):
        pic, od, i4 = jload(old['dev|%s|picasso' % g]), jload(old['dev|%s|dual' % g]), jload(raw['dev|%s|iss4' % g])
        dev = {'R0': pic['R0']['E'], 'R2E': pic['R2E']['E'], 'first-study ISS': od['E'], 'ISS4': i4['E']}
        te = {'R0': resA[g]['rows']['R0']['E'], 'R2E': resA[g]['rows']['R2E']['E'], 'ISS4': resA[g]['rows']['ISS4']['E']}
        if 'ISS_OLD' in resA[g]['rows']:
            te['first-study ISS'] = resA[g]['rows']['ISS_OLD']['E']
        keys = [k for k in ('R0', 'R2E', 'first-study ISS') if k in te]
        d_order, t_order = sorted(keys, key=lambda k: -dev[k]), sorted(keys, key=lambda k: -te[k])
        rev = 'exactly reversed' if t_order == d_order[::-1] else 'not reversed'
        out.append('%s: development ranking (dev E_s) %s; test ranking (E_s) %s; the order of %s is %s between the development and the test frames.' % (g, rank_str(dev), rank_str(te), ', '.join(keys), rev))
    if 'ISS_OLD-R2E' in resA['D2']['pairs']:
        a, b = resA['D2']['pairs']['ISS_OLD-R2E'], resA['D2']['pairs']['ISS_OLD-R0']
        out.append('Continuity row (D2; POST HOC, the first study\'s frozen two-dial point, outside the registered rule and on a void pooled set, so not evidence of an advantage): ISS_OLD - R2E %s and ISS_OLD - R0 %s; '
                   'against the registered ISS4 - R2E %s the sign of the D2 contrast depends on which frozen ISS point is used.' %
                   (ci_str(a['dE'], a['dE_ci']), ci_str(b['dE'], b['dE_ci']), ci_str(resA['D2']['pairs']['ISS4-R2E']['dE'], resA['D2']['pairs']['ISS4-R2E']['dE_ci'])))
    return out


def build_md(raw, n_boot=N_BOOT, log=lambda s: None):
    old = old_raw()
    gates = jload(raw['ref|gates'])
    lines = []
    A = lines.append
    A('# Second real-data pore study: fresh frames, a registration gate, a repaired ISS tuning grid, and an overlap test (pre-registered)')
    A('')
    A('Generated by `experiments/exp_smlm_real2.py --from-raw` from `results/smlm/real2_raw.npz`; every number below is computed from that raw (the file regenerates byte for byte). '
      'The specification is `results/smlm/real2_plan.md` (declared 2026-11-09 14:32 UTC; Amendment 1 declared before any localizer ran). **The frozen development choices were hashed '
      'at %s (sha256 `%s`) before the test run; the test run started at %s and ended at %s; one test run per test set.** The first study\'s raw (`real_gt_raw.npz`) is read, '
      'never written.' % (str(raw['freeze|time_utc']), str(raw['freeze|sha256']), str(raw['test|start_utc']), str(raw['test|end_utc'])))
    A('')
    # gates
    A('## Validity gates (plan sections 2 and 6)')
    A('')
    r2 = {n: gate_r2(raw, old, n) for n in TESTS}
    val = {n: set_valid(raw, old, n, r2[n]) for n in TESTS}
    rows = []
    for n in TESTS:
        g, r1 = gates[n], gate_r1(raw, n)
        rows.append([n, '%s (%s)' % (TESTS[n]['ds'], 'frames %d-%d' % TESTS[n]['frames']), fmt(g['V1_abs_diff'], 2), g['scored_pores'], fmt(g['R_ref'], 2), fmt(r1['max_dev']), fmt(r1['min_q']),
                     '; '.join('%s %s/%s nm' % (r, fmt(x['median']), fmt(x['p90'])) for r, x in r2[n].items()), 'yes' if val[n]['valid'] else 'NO (%s)' % ', '.join(k for k, v in val[n].items() if not v and k != 'valid')])
    A(md_table(['set', 'recording', 'V1: R_ref odd - even (nm, limit 2)', 'V4: scored pores (limit 100)', 'R_ref (nm)', 'R1: largest block deviation (nm, limit 20)', 'R1: smallest quality (limit 2.5)',
                'R2: median / 90th percentile of the residual shift per row', 'valid'], rows))
    A('')
    A('Gate V2 (D2): R_ref of fields 17 and 19 differ by %s nm (limit 3): %s. Gate V3 (D2, 8th harmonic significant in at least 50 %% of the pores; used for the C8 secondary only): '
      'field 19 %s, field 17 %s.' % (fmt(gates['V2']['diff'], 2), 'pass' if gates['V2']['pass'] else 'FAIL', fmt(100 * gates['T2a']['frac_8fold_significant'], 0) + ' %', fmt(100 * gates['T2b']['frac_8fold_significant'], 0) + ' %'))
    A('')
    A('Registration shift curves (raw minus table, nm; polynomial of the frame index through the blocks of about 125 frames):')
    A('')
    rows = []
    for n in TESTS:
        sh = shift_of(raw, n)
        rows.append([n, fmt(sh[0, 0], 0) + ', ' + fmt(sh[0, 1], 0), fmt(sh[len(sh) // 2, 0], 0) + ', ' + fmt(sh[len(sh) // 2, 1], 0), fmt(sh[-1, 0], 0) + ', ' + fmt(sh[-1, 1], 0), len(raw['reg|%s|blocks' % n])])
    A(md_table(['set', 'first frame (x, y)', 'middle frame', 'last frame', 'blocks'], rows))
    A('')
    # data
    A('## Data, splits and the reference')
    A('')
    rows = []
    for n, t in TESTS.items():
        g = gates[n]
        win = test_window(raw, n)
        rows.append([n, t['ds'], 'field %s' % t['field'] if t['field'] else 'one field', '%d-%d (%d)' % (t['frames'][0], t['frames'][1], n_frames_of(n)), 'whole frame' if win is None else 'window origin %s px' % (win,),
                     g['scored_pores'], fmt(g['R_ref'], 2), fmt(g['n_exp'], 0), ', '.join('%d-%d' % r for r in t['ref_ranges'])])
    A(md_table(['set', 'dataset', 'recording', 'frames (1-based)', 'scored image', 'scored pores', 'R_ref (nm)', 'n_exp (ring localizations)', 'T_ref frames'], rows))
    A('')
    A('D1 = Heil et al. 2018 (gp210, Zenodo 5068525); D2 = Thevathasan et al. 2019 (Nup96, Zenodo 3540443). Development sets: the first study\'s D1dev (frames 1-300) and D2dev (field 17 frames 1-300). '
      'The pooled D2 uses one ring radius R = %s nm (the mean of the two fields\' R_ref) so that the radial bins agree.' % fmt(common_R(raw, 'T2a'), 2))
    A('')
    # development
    A('## Development tuning (frames 1-300; E_s; Part B: J_B on merged development frames)')
    A('')
    rows = []
    for ds in ('D1', 'D2'):
        i4 = jload(raw['dev|%s|iss4' % ds])
        pic = jload(old['dev|%s|picasso' % ds])
        od = jload(old['dev|%s|dual' % ds])
        rows.append([ds, fmt(pic['R0']['E'], 2), fmt(pic['R2E']['E'], 2), pic['R2E']['n_configs'], fmt(od['E'], 2), '%s, kappa %g, width x%g' % (od['cfg'], od['kappa'], od['smax_rel']), fmt(i4['E'], 2),
                     '%s, kappa %g, fw %g, smin %s, smax %g' % (i4['cfg'], i4['kappa'], i4['fw'], 'none' if i4['smin'] is None else '%g' % i4['smin'], i4['smax']), i4['n_points_scored']])
    A(md_table(['set', 'Picasso R0 dev E_s', 'Picasso R2E dev E_s', 'R2E configurations scored', 'first-study ISS dev E_s', 'its point', 'ISS4 dev E_s', 'ISS4 frozen point', 'ISS4 operating points scored'], rows))
    A('')
    A('ISS4 per configuration (dev E_s after the four-dial ascent):')
    A('')
    rows = []
    for ds in ('D1', 'D2'):
        i4 = jload(raw['dev|%s|iss4' % ds])
        for c, v in i4['per_cfg'].items():
            p = v['point']
            rows.append([ds, c, fmt(v['E'], 2), '%g' % p['kappa'], '%g' % p['fw'], 'none' if p['smin'] is None else '%g' % p['smin'], '%g' % p['smax'], v['rounds'], v['n_points']])
    A(md_table(['set', 'configuration', 'dev E_s', 'kappa', 'fw', 'smin', 'smax', 'rounds', 'points'], rows))
    A('')
    A('Part B development (merged frames of the development frames; the Picasso structure kept is named in the log of the run):')
    A('')
    rows = []
    for ds in ('D1', 'D2'):
        b = jload(raw['dev|%s|B' % ds])
        rows.append([ds, b['n_groups'], fmt(b['ph_bright'], 0), fmt(b['rm'], 0), fmt(b['close_nm'], 0), fmt(b['R0']['J'], 3), '%s (%s)' % (fmt(b['R2B']['J'], 3), b['R2B']['structure']), b['R2B']['n_points'], fmt(b['ISSB']['J'], 3),
                     'kappa %g, fw %g, smin %s, smax %g' % (b['ISSB']['kappa'], b['ISSB']['fw'], 'none' if b['ISSB']['smin'] is None else '%g' % b['ISSB']['smin'], b['ISSB']['smax']), b['ISSB']['n_points']])
    A(md_table(['set', 'merged frames', 'bright threshold (table photons)', 'match radius (nm)', 'close-pair distance (nm)', 'R0 dev J_B', 'R2B dev J_B (structure)', 'R2B points', 'ISS-B dev J_B', 'ISS-B point', 'ISS-B points'], rows))
    A('')
    # part A
    resA = {g: analyse_A(raw, old, g, n_boot, log) for g in ('D1', 'D2')}
    A('## Part A (primary): registered localizations, tuned and frozen points')
    A('')
    reads = {g: read_A(resA[g]) for g in resA}
    validg = {'D1': val['T1']['valid'], 'D2': val['T2a']['valid'] and val['T2b']['valid']}
    rows = []
    for g in ('D1', 'D2'):
        pr = resA[g]['pairs']['ISS4-R2E']
        rows.append([g, '+'.join(resA[g]['members']), 'yes' if validg[g] else 'VOID', ci_str(pr['dE'], pr['dE_ci']), '[%s, %s]' % (fmt(pr['dE_ci_pores'][0]), fmt(pr['dE_ci_pores'][1])),
                     ci_str(pr['dJ'] * 100, [100 * pr['dJ_ci'][0], 100 * pr['dJ_ci'][1]]), ci_str(pr['dw'], pr['dw_ci']), '%s [%s, %s]' % (fmt(pr['off_ring_ratio'], 2), fmt(pr['off_ring_ratio_ci'][0], 2), fmt(pr['off_ring_ratio_ci'][1], 2)),
                     with_void(reads[g], validg[g])])
    A(md_table(['dataset', 'sets', 'valid', 'Delta E_s, ISS4 - Picasso R2E (E points) [95 % block-bootstrap interval]', 'pore-resampling interval', 'Delta J_s (points)', 'Delta w at equal yield (nm), secondary',
                'ISS ring yield at the off-ring count of R2E, over R2E\'s yield', 'reading'], rows))
    A('')
    A('**Registered reading of rule A: %s**' % overall_A(reads, validg))
    A('')
    A('Positive Delta E_s favours the ISS; the intervals are paired percentile intervals of %d resamples of blocks of %d consecutive frames (seed %d, the same blocks for every row; the members of the pooled D2 '
      'are resampled independently); the pore-resampling intervals resample the pores with all their frames (seed %d).' % (n_boot, BLOCK, BOOT_SEED, BOOT_SEED_PORES))
    A('')
    for t in transfer_paragraph(raw, old, resA):
        A(t)
        A('')
    edge4 = {ds: sum(1 for v in jload(raw['dev|%s|iss4' % ds])['per_cfg'].values() if v['point']['kappa'] == KAPPAS2[-1]) for ds in ('D1', 'D2')}
    kB = {ds: jload(raw['dev|%s|B' % ds])['ISSB']['kappa'] for ds in ('D1', 'D2')}
    A(grid_edges_line(edge4, kB))
    A('')
    for g in ('D1', 'D2'):
        r = resA[g]
        A('### %s (%s): all rows (registered)' % (g, ' + '.join(r['members'])))
        A('')
        A('%d scored frames, %d scored pores, expected ring yield n_exp %s, ring radius %s nm.' % (r['n_frames'], r['n_pores'], fmt(r['n_exp'], 0), fmt(r['R_ref'], 1)))
        A('')
        A(md_table(['row', 'localizations', 'ring N', 'w (nm) [interval]', 'floor', 'J_s [interval]', 'E_s [interval]', 'mu - R_ref (nm)', 'ring loc. per pore per 1000 frames', 'duplicates within sigma0 (%)'],
                   [row_line_A(ROW_NAMES[k], r['rows'][k]) for k in ROWS_A if k in r['rows']]))
        A('')
        pr_rows = []
        for k, v in r['pairs'].items():
            a_, b_ = k.split('-')
            pr_rows.append([ROW_NAMES[a_] + ' - ' + ROW_NAMES[b_], ci_str(v['dE'], v['dE_ci']), ci_str(100 * v['dJ'], [100 * v['dJ_ci'][0], 100 * v['dJ_ci'][1]])])
        A(md_table(['paired difference', 'Delta E_s (points)', 'Delta J_s (points)'], pr_rows))
        A('')
        cur = r['curve']
        A('The ISS4 curve in kappa (all other dials frozen): ' + '; '.join('kappa %g: N %s, floor %s, E_s %s' % (float(k[1:]), fmt(v['N'], 0), fmt(v['floor'], 0), fmt(v['E'])) for k, v in cur.items()))
        A('')
        if 'c8' in r['rows']['R0']:
            A('Corner contrast C8 (mean cos 8(theta - phi), pooled): ' + '; '.join('%s %s (null %s)' % (k, fmt(r['rows'][k]['c8'], 4), fmt(r['rows'][k]['c8_null'], 4)) for k in ROWS_A if k in r['rows']))
            A('')
    # unregistered side by side
    A('### The unregistered analysis (the shift curve NOT applied), side by side; Delta E_s only')
    A('')
    rows = []
    for g in ('D1', 'D2'):
        ru = analyse_A(raw, old, g, n_boot, log, register=False, with_curve=False)
        pr, pu = resA[g]['pairs']['ISS4-R2E'], ru['pairs']['ISS4-R2E']
        rows.append([g, ci_str(pr['dE'], pr['dE_ci']), ci_str(pu['dE'], pu['dE_ci']), fmt(resA[g]['rows']['ISS4']['E']), fmt(ru['rows']['ISS4']['E']), fmt(resA[g]['rows']['R2E']['E']), fmt(ru['rows']['R2E']['E'])])
    A(md_table(['dataset', 'Delta E_s registered', 'Delta E_s unregistered', 'ISS4 E_s registered', 'ISS4 E_s unregistered', 'R2E E_s registered', 'R2E E_s unregistered'], rows))
    A('')
    # per member
    A('### Per test set (descriptive): Delta E_s, ISS4 - R2E, and the gate status')
    A('')
    rows = []
    for n in TESTS:
        g = [k for k, v in GROUPS.items() if n in v][0]
        if len(GROUPS[g]) == 1:
            r = resA[g]
        else:
            r = analyse_A_single(raw, old, n, n_boot)
        pr = r['pairs']['ISS4-R2E']
        rows.append([n, 'yes' if val[n]['valid'] else 'VOID: ' + ', '.join(k for k, v in val[n].items() if not v and k != 'valid'), ci_str(pr['dE'], pr['dE_ci']), fmt(r['rows']['ISS4']['E']), fmt(r['rows']['R2E']['E']),
                     fmt(r['rows']['ISS4']['fit']['N'], 0), fmt(r['rows']['R2E']['fit']['N'], 0), fmt(r['rows']['ISS4']['fit']['floor'], 0), fmt(r['rows']['R2E']['fit']['floor'], 0),
                     fmt(r['rows']['PUB']['mu_bias']), fmt(r['rows']['PUB']['fit']['w'])])
    A(md_table(['set', 'valid', 'Delta E_s', 'ISS4 E_s', 'R2E E_s', 'ISS4 N', 'R2E N', 'ISS4 floor', 'R2E floor', 'published table: ring radius minus R_ref (nm)', 'published table: ring width (nm)'], rows))
    A('')
    A('(POST HOC observation, not a registered gate.) The published table of the same frames is the best available picture of the reference: where its own radial profile peaks far from R_ref the pore list '
      'of that window is not a clean set of single rings (touching pores, clusters), and a fit of the localizers\' profile there is degenerate (ring yield at its lower bound). The fit of a row at N = 0 '
      'means exactly that, not a measured absence of localizations.')
    A('')
    # part B
    resB = {g: analyse_B(raw, old, g, n_boot, log) for g in ('D1', 'D2')}
    readsB = {g: read_B(resB[g]) for g in resB}
    A('## Part B: merged frames (K = %d frames added; groups of frames %d apart in blocks of %d)' % (K_MERGE, GAP_MERGE, B_BLOCK))
    A('')
    rows = []
    for g in ('D1', 'D2'):
        pr = resB[g]['pairs']['ISSB-R2B']
        rows.append([g, '+'.join(resB[g]['members']), resB[g]['n_groups'], 'yes' if validg[g] else 'VOID', ci_str(100 * pr['J']['est'], [100 * pr['J']['ci'][0], 100 * pr['J']['ci'][1]]),
                     ci_str(100 * pr['recall']['est'], [100 * pr['recall']['ci'][0], 100 * pr['recall']['ci'][1]]), ci_str(100 * pr['close_recall']['est'], [100 * pr['close_recall']['ci'][0], 100 * pr['close_recall']['ci'][1]]),
                     ci_str(100 * pr['pair_resolution']['est'], [100 * pr['pair_resolution']['ci'][0], 100 * pr['pair_resolution']['ci'][1]]), with_void(readsB[g], validg[g])])
    A(md_table(['dataset', 'sets', 'merged frames', 'valid', 'Delta J_B, ISS-B - R2B (points) [95 % interval]', 'Delta recall of bright events (points)', 'Delta close-pair recall (points)',
                'Delta pair resolution (points)', 'reading'], rows))
    A('')
    A('**Registered reading of rule B: %s**' % overall_B(readsB, validg))
    A('')
    for g in ('D1', 'D2'):
        r = resB[g]
        A('### %s: rows of Part B' % g)
        A('')
        A(md_table(['row', 'J_B [interval]', 'recall [interval]', 'unconfirmed detections per merged frame', 'detections per merged frame', 'close-pair recall [interval]', 'pair resolution [interval]'],
                   [[ROW_NAMES_B[k], ci_str(100 * v['J'], [100 * x for x in v['ci']['J']]), ci_str(100 * v['recall'], [100 * x for x in v['ci']['recall']]), fmt(v['fp_per_frame'], 2), fmt(v['det_per_frame'], 1),
                      ci_str(100 * v['close_recall'], [100 * x for x in v['ci']['close_recall']]), ci_str(100 * v['pair_resolution'], [100 * x for x in v['ci']['pair_resolution']])] for k, v in r['rows'].items()]))
        A('')
        A('Bookkeeping: ' + '; '.join('%s: %d merged frames, %s reference events per frame (%s bright), match radius %s nm, close-pair distance %s nm' % (n, i['n_groups'], fmt(i['ref_per_frame'], 1), fmt(i['bright_per_frame'], 1), fmt(i['rm'], 0), fmt(i['close_nm'], 0))
                                for n, i in r['info'].items()) + '. Close events: %s, pairs: %s (ISS-B row).' % (fmt(r['rows']['ISSB']['n_close'], 0), fmt(r['rows']['ISSB']['n_pairs'], 0)))
        A('')
    # timing
    A('## Cost (CPU seconds per frame, one core, the machine shared with other jobs)')
    A('')
    rows = []
    for n in TESTS:
        ds = TESTS[n]['ds']
        rows.append([n, fmt(float(np.mean(raw['test|%s|R0|secs' % n])), 4), fmt(float(np.mean(raw['test|%s|R2E|secs' % n])), 4), fmt(float(np.mean(raw['testpath|%s|iss4|secs' % n])), 2),
                     fmt(float(np.mean(raw['testpath|%s|iss4|entries' % n])), 1), fmt(100 * float(np.mean(raw['testpath|%s|iss4|entries' % n] >= CAP)), 1)])
    A(md_table(['set', 'Picasso R0', 'Picasso R2E', 'ISS4 path', 'ISS4 mean entries per frame', 'frames at the cap of %d entries (%%)' % CAP], rows))
    A('')
    # diag
    A('## Development-only diagnostics of the grid-edge optimum of the first study')
    A('')
    rows = []
    for ds in ('D1', 'D2'):
        if 'diag|%s' % ds in raw:
            d = jload(raw['diag|%s' % ds])
            rows.append([ds, d['cfg'], d['n_frames'], fmt(d['robust_sd_frame'], 3), fmt(d['robust_sd_residual'], 3), '%.2e / %.2e / %.2e' % (d['tail_pos_4'], d['tail_neg_4'], d['tail_gauss_4']),
                         ' / '.join(fmt(x, 2) for x in d['ring_fraction_by_amplitude_fifth']), fmt(d['ring_area_fraction'], 2), fmt(d['units_per_frame'], 1)])
    A(md_table(['set', 'configuration', 'frames', 'robust sd of the standardized high-passed frame', 'robust sd of the standardized residual away from units', 'P(r > 4 sd) / P(r < -4 sd) / Gaussian',
                'ring fraction of the units in the discs, by amplitude fifth (low to high)', 'ring share of the disc area', 'units per frame (kappa 4)'], rows))
    A('')
    A('The label-free repair of the law (its variance multiplied by the residual\'s variance), the four dials re-tuned on the 300 development frames:')
    A('')
    rows2 = []
    for ds in ('D1', 'D2'):
        if 'diag|%s' % ds in raw and 'rescaled' in jload(raw['diag|%s' % ds]):
            r_ = jload(raw['diag|%s' % ds])['rescaled']
            rows2.append([ds, fmt(r_['variance_factor'], 3), fmt(r_['E_stored'], 2) if r_['E_stored'] is not None else 'n/a', fmt(r_['E'], 2), '%g' % r_['point']['kappa'], '%g' % r_['point']['fw'],
                          'none' if r_['point']['smin'] is None else '%g' % r_['point']['smin'], '%g' % r_['point']['smax']])
    A(md_table(['set', 'variance factor', 'dev E_s, noise law as calibrated (four dials)', 'dev E_s, variance rescaled by the residual', 'kappa', 'fw', 'smin', 'smax'], rows2))
    A('')
    A('How to read: a residual sd below 1 means the label-free noise law is too LARGE (conservative), above 1 too small; a tail on the positive side only is unfitted signal, not noise; a ring fraction '
      'that rises with amplitude is why the amplitude and width post-filters that Picasso\'s tuned pipeline has matter for the position of the optimum in kappa.')
    A('')
    A('## Not done, deviations, weaknesses')
    A('')
    A('* See plan sections 1, 4 and 5 for the weaknesses declared in advance: T1 is the same field (and pores) as the first D1 test, with new frames; D2 has three recordings of one preparation and no between-microscope variance; '
      'the reference is derived from third-party tables of other frames of the same recordings; Part B is a manufactured density with the table as reference.')
    edge = {ds: sum(1 for v in jload(raw['dev|%s|iss4' % ds])['per_cfg'].values() if v['point']['kappa'] == KAPPAS2[-1]) for ds in ('D1', 'D2')}
    A('* Grid edge: the development optimum of the ISS4 ascent is at the TOP of the extended kappa grid (%g) for %d of 8 configurations of D1 and %d of 8 of D2 (the plan did not extend the grid '
      'further; the frozen point is at kappa %g (D1) and %g (D2)).' % (KAPPAS2[-1], edge['D1'], edge['D2'], jload(raw['dev|D1|iss4'])['kappa'], jload(raw['dev|D2|iss4'])['kappa']))
    A('* The stored raw keeps every frozen choice, the per-configuration development results and the test paths; the development ISS paths (49 MB, `devpath|*`) were removed from it after the '
      'freeze to keep the file small (`--step dev` recomputes them; the freeze hash does not cover them).')
    A('* Computed locally on a shared 4-core machine (no Kaggle); the optional offer to move the long jobs to private CPU notebooks was not taken because the setup would have cost more than the jobs.')
    A('* No raw-frame DNA-origami or calibration-bead data with a known spacing were reachable (plan section 1).')
    A('* The greedy-order row was not run in this study (the controls of `results/smlm/controls.md` and `results/astro/controls.md` cover the order of entry).')
    A('* Everything labelled POST-AMENDMENT (the verdict for T1 under Amendment 1) is not part of the registered verdict.')
    A('')
    return '\n'.join(lines) + '\n'


def analyse_A_single(raw, old, name, n_boot):
    """Part A for one member of the pooled group (descriptive)."""
    saved = GROUPS.get('_one')
    GROUPS['_one'] = [name]
    try:
        return analyse_A(raw, old, '_one', n_boot, with_curve=True)
    finally:
        if saved is None:
            GROUPS.pop('_one', None)
        else:
            GROUPS['_one'] = saved


def write_report(raw_path=RAW, md_path=MD, n_boot=N_BOOT, log=print):
    raw = load_raw(raw_path)
    md = build_md(raw, n_boot, log)
    with open(md_path, 'w') as f:
        f.write(md)
    log('wrote %s' % md_path)
    return md


# ==================================================================================================================================
# command line
# ==================================================================================================================================
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('--step', choices=('ref', 'reg', 'cal', 'dev', 'test', 'diag'))
    ap.add_argument('--from-raw', action='store_true')
    a = ap.parse_args(argv)
    if a.from_raw:
        write_report()
        return
    raw, old = load_raw(), old_raw()
    save = save_raw
    if a.step == 'ref':
        step_ref(raw, old)
        save(raw)
    elif a.step == 'reg':
        step_reg(raw, old)
        save(raw)
    elif a.step == 'cal':
        step_cal(raw, old)
        save(raw)
    elif a.step == 'dev':
        step_dev(raw, old, save=save)
    elif a.step == 'test':
        step_test(raw, old, save=save)
    elif a.step == 'diag':
        step_diag(raw, old)
        save(raw)


if __name__ == '__main__':
    main()
