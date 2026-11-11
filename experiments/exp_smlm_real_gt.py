"""
Real SMLM recordings of nuclear pore complexes (known geometry): the sliding inverse scale space against Picasso (pre-registered)
=================================================================================================================================

Implements `results/smlm/real_gt_plan.md` (declared 2026-11-09 07:56 UTC; Amendments 1 and 2 declared before any pore finder or localizer ran on
development or test frames) and writes `results/smlm/real_gt.md` from `results/smlm/real_gt_raw.npz`. Read the plan first: it is the specification, this
docstring only says where each piece of it lives.

Data (fetched by `data/smlm/fetch_real_gt.py` into $SMLM_REAL_GT, default /home/user/smlm_real_data; the TIFF reading needs `tifffile`, the report from the
raw does not):
    D1  Heil et al. 2018, gp210 nuclear pores (Xenopus), STORM on glass: 20000 frames of 117 x 118 px, 102 nm, EM gain 100, 15.01 e-/ADU, offset 100 ADU.
        Reference table T_ref: the authors' ThunderSTORM table, frames 2001-10000 and 12001-20000. Development frames 1-300 (pores with x < 5800 nm), test
        frames 10001-12000 (pores with x > 6200 nm).
    D2  Thevathasan et al. 2019, Nup96-SNAP, 2D STORM, 106 nm, 419 x 419 px: development = field 17 frames 1-300, test = field 18 frames 1-1055, one
        128 x 128 px window per field; T_ref: the SMAP table of the same field, frames 1057-30000.
Parts (each stores its results in the raw; run in this order, one heavy process at a time):
    pores   step 2: the reference pores, the radius, the 8-fold phases and the validity gates from T_ref ONLY (no localizer is run); writes and prints the
            sha256 of the frozen pore lists.
    cal     step 3: label-free calibration of the ISS on the frames of each set (no labels exist), the Picasso scale factors from the development frames.
    dev     step 4: tuning of Picasso (R0, R1E, R2E) and of the ISS and greedy rows (168 operating points) on the development frames for E_s; freezes and hashes.
    test    step 5: ONE test run per dataset with the frozen configurations.
    (--from-raw regenerates results/smlm/real_gt.md from the raw alone, byte for byte.)
"""
import argparse
import hashlib
import json
import os
import platform
import sys
import time

os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import numpy as np                                                   # noqa: E402
from scipy import ndimage as ndi, signal                             # noqa: E402
from scipy.optimize import minimize                                  # noqa: E402
from scipy.spatial import cKDTree                                    # noqa: E402
from scipy.special import ndtr                                       # noqa: E402

HERE = os.path.abspath(os.path.dirname(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
sys.path.insert(0, REPO)
import smlm_common as C                                              # noqa: E402

DATA = os.environ.get('SMLM_REAL_GT', '/home/user/smlm_real_data')
OUT = os.environ.get('REALGT_OUT', os.path.join(REPO, 'results', 'smlm'))
RAW = os.path.join(OUT, 'real_gt_raw.npz')
MD = os.path.join(OUT, 'real_gt.md')

# ==================================================================================================================================
# the specification of the plan (constants fixed in the plan before any run)
# ==================================================================================================================================
HIST_NM = 10.0                       # histogram pixel of the pore finder
RADIAL_BIN_NM = 5.0                  # radial bin of the profile fit
BLOCK = 10                           # frames per resampled block
N_BOOT = 1000
BOOT_SEED = 0
BOOT_SEED_PORES = 1
N_THIN = 20
MARGIN_E = 3.0                       # E points: the equivalence margin
GATE_DEV_E = 10.0                    # E points: the development gate
ALPHA = 1.0                          # per nm, the efficiency of the challenge
MIN_ANNULUS = 15                     # localizations in the annulus for an accepted pore
RADIUS_TOL = 0.25                    # accepted radius within +-25 % of the median of the candidates
SETS = {
    'D1dev': dict(ds='D1', role='dev'), 'D1test': dict(ds='D1', role='test'),
    'D2dev': dict(ds='D2', role='dev'), 'D2test': dict(ds='D2', role='test'),
}
DS = {
    'D1': dict(name='Heil et al. 2018, gp210, Xenopus, STORM on glass (Zenodo 5068525)', pixel_nm=102.0, shape=(117, 118), n_frames=20000, r_start_px=8.0,
               annulus=(40.0, 125.0), search=130.0, pct=95.0, disc_margin=150.0, tau_fp=120.0, tau=40.0, ring_half=40.0,
               camera=dict(baseline=100.0, gain=100.0, sensitivity=15.01, qe=0.9), hist_size_nm=12100.0,
               ref_frames=[(2001, 10000), (12001, 20000)], dev_frames=(1, 300), test_frames=(10001, 12000), x_dev_max=5800.0, x_test_min=6200.0,
               table='zenodo_5068525/NPCData_glass_STORM_loc.csv', tab_cols=dict(frame=1, x=2, y=3), tab_skip=1, n_pores_min=100, eight_fold=False),
    'D2': dict(name='Thevathasan et al. 2019, Nup96-SNAP, 2D STORM (Zenodo 3540443)', pixel_nm=106.0, shape=(419, 419), n_frames=30000, r_start_px=5.3,
               annulus=(25.0, 80.0), search=90.0, pct=97.0, disc_margin=100.0, tau_fp=75.0, tau=25.0, ring_half=25.0,
               camera=dict(baseline=100.0, gain=1.0, sensitivity=1.0, qe=0.95), hist_size_nm=45000.0, window_px=128,
               ref_frames=[(1057, 30000)], dev_frames=(1, 300), test_frames=(1, 1055), n_pores_min=100, eight_fold=True,
               fields=dict(dev=17, test=18), table='zenodo_3540443/Cos7_Nup96_A647_COT_Au_%d_1_locs.csv', tab_cols=dict(frame=0, x=12, y=13), tab_skip=1,
               prefix='zenodo_3540443/Cos7_Nup96_A647_COT_Au_%d_1_MMStack_Pos0.ome.tif.part0-400MB', n_prefix_frames=1055),
}
R_LITERATURE_D2 = 53.5               # nm, from memory (nature.com is blocked); gate V2 only


def fmt(x, nd=1):
    return ('%.' + str(nd) + 'f') % x


# ==================================================================================================================================
# data access
# ==================================================================================================================================
def load_table(ds, field=None):
    """The published localization table of the dataset: dict(frame (1-based), x_nm, y_nm)."""
    sp = DS[ds]
    path = os.path.join(DATA, sp['table'] % field if '%d' in sp['table'] else sp['table'])
    a = np.genfromtxt(path, delimiter=',', skip_header=sp['tab_skip'], usecols=[sp['tab_cols'][k] for k in ('frame', 'x', 'y')])
    a = a[np.isfinite(a).all(axis=1)]
    return dict(frame=a[:, 0].astype(np.int64), x=a[:, 1], y=a[:, 2])


def table_select(tab, frame_ranges):
    m = np.zeros(len(tab['frame']), bool)
    for a, b in frame_ranges:
        m |= (tab['frame'] >= a) & (tab['frame'] <= b)
    return {k: v[m] for k, v in tab.items()}


def load_frames(set_name, win=None):
    """The frames of a set in ADU, float32 (n, ny, nx), cropped to the window win = (x0, y0) px (D2) and the frame range of the set."""
    import tifffile                                                    # only the raw data need it, not the report
    s = SETS[set_name]
    sp = DS[s['ds']]
    a, b = sp['dev_frames'] if s['role'] == 'dev' else sp['test_frames']
    if s['ds'] == 'D1':
        with tifffile.TiffFile(os.path.join(DATA, 'zenodo_5068525', 'NPCData_glass_STORM.tif')) as t:
            fr = t.asarray(key=range(a - 1, b))
        return fr.astype(np.float32)
    field = sp['fields'][s['role']]
    with tifffile.TiffFile(os.path.join(DATA, sp['prefix'] % field)) as t:
        fr = t.asarray(key=range(a - 1, b))
    n = sp['window_px']
    return fr[:, win[1]:win[1] + n, win[0]:win[0] + n].astype(np.float32)


def calibration_frames(set_name, win=None):
    """D2: the label-free calibration uses all 1055 frames of the field (the plan, section 2); D1: the frames of the set."""
    import tifffile
    s = SETS[set_name]
    sp = DS[s['ds']]
    if s['ds'] == 'D1':
        return load_frames(set_name)
    field = sp['fields'][s['role']]
    with tifffile.TiffFile(os.path.join(DATA, sp['prefix'] % field)) as t:
        fr = t.asarray(key=range(0, sp['n_prefix_frames']))
    n = sp['window_px']
    return fr[:, win[1]:win[1] + n, win[0]:win[0] + n].astype(np.float32)


# ==================================================================================================================================
# step 2: the reference pores (from T_ref only)
# ==================================================================================================================================
def ring_candidates(x, y, region, r_start_px, pix, pct):
    """Ring matched filter on a pix-nm histogram of the localizations inside region = (x0, y0, x1, y1): the local maxima above the pct-th percentile of the
    occupied pixels, as an (n, 2) array of centres (nm, absolute)."""
    x0, y0, x1, y1 = region
    nx, ny = int(np.ceil((x1 - x0) / pix)), int(np.ceil((y1 - y0) / pix))
    H, _, _ = np.histogram2d(y, x, bins=[y0 + pix * np.arange(ny + 1), x0 + pix * np.arange(nx + 1)])
    half = int(2.5 * r_start_px)
    yy, xx = np.mgrid[-half:half + 1, -half:half + 1]
    k = np.exp(-0.5 * ((np.hypot(xx, yy) - r_start_px) / (0.2 * r_start_px)) ** 2)
    k -= k.mean()
    k /= np.linalg.norm(k)
    c = signal.fftconvolve(H, k, mode='same')
    occupied = c[H > 0]
    if len(occupied) == 0:
        return np.zeros((0, 2))
    pk = (c == ndi.maximum_filter(c, size=int(2 * r_start_px) + 1)) & (c > np.percentile(occupied, pct))
    iy, ix = np.nonzero(pk)
    return np.c_[x0 + ix * pix + pix / 2, y0 + iy * pix + pix / 2]


def circle_fit(tree, x, y, a, b, r_in, r_out, search, iters=3, min_n=MIN_ANNULUS):
    """Three iterations of an algebraic circle fit of the localizations in the annulus (r_in, r_out) around (a, b); returns (a, b, R, n) or None."""
    n, R = 0, np.nan
    for _ in range(iters):
        idx = tree.query_ball_point([a, b], search)
        P = np.c_[x[idx], y[idx]] - [a, b]
        r = np.hypot(P[:, 0], P[:, 1])
        m = (r > r_in) & (r < r_out)
        n = int(m.sum())
        if n < min_n:
            return None
        Q = P[m]
        sol = np.linalg.lstsq(np.c_[2 * Q[:, 0], 2 * Q[:, 1], np.ones(n)], (Q ** 2).sum(1), rcond=None)[0]
        a, b = a + sol[0], b + sol[1]
        R = np.sqrt(max(sol[2] + sol[0] ** 2 + sol[1] ** 2, 0.0))
    return a, b, R, n


def find_pores(x, y, region, sp):
    """The pore finder of the plan (section 3.1) on the T_ref localizations (x, y) in nm: matched filter, circle fits, acceptance (>= 15 annulus localizations,
    radius within +-25 % of the median of the candidates), merging of centres closer than 1.5 R_start. Returns (centres (n, 2), radii (n), counts (n))."""
    pix = HIST_NM
    cand = ring_candidates(x, y, region, sp['r_start_px'], pix, sp['pct'])
    tree = cKDTree(np.c_[x, y])
    fits = []
    for a, b in cand:
        f = circle_fit(tree, x, y, a, b, sp['annulus'][0], sp['annulus'][1], sp['search'])
        if f is not None:
            fits.append(f)
    if not fits:
        return np.zeros((0, 2)), np.zeros(0), np.zeros(0, np.int64)
    fits = np.array(fits)
    med = np.median(fits[:, 2])
    ok = np.abs(fits[:, 2] - med) <= RADIUS_TOL * med
    fits = fits[ok]
    order = np.argsort(-fits[:, 3], kind='stable')
    keep = []
    dmin = 1.5 * sp['r_start_px'] * pix
    kt = []
    for i in order:
        c = fits[i, :2]
        if kt and np.min(np.hypot(*(np.array(kt) - c).T)) < dmin:
            continue
        keep.append(i)
        kt.append(c)
    keep = np.array(sorted(keep))
    f = fits[keep]
    return f[:, :2], f[:, 2], f[:, 3].astype(np.int64)


def pore_radius(x, y, centres, sp):
    """The median circle-fit radius of the pores with the given centres (kept fixed as the start) on localizations (x, y); NaN if none fits."""
    tree = cKDTree(np.c_[x, y])
    rs = []
    for a, b in centres:
        f = circle_fit(tree, x, y, a, b, sp['annulus'][0], sp['annulus'][1], sp['search'])
        if f is not None:
            rs.append(f[2])
    return float(np.median(rs)) if rs else float('nan')


def pore_phases(x, y, centres, R, sp):
    """The 8-fold phase of every pore (maximum of the 8th harmonic over phi in [0, 45 degrees)), the harmonic's modulus over n, the Rayleigh p-value, n."""
    tree = cKDTree(np.c_[x, y])
    out = np.zeros((len(centres), 4))
    for j, (a, b) in enumerate(centres):
        idx = tree.query_ball_point([a, b], R + sp['ring_half'])
        P = np.c_[x[idx], y[idx]] - [a, b]
        r = np.hypot(P[:, 0], P[:, 1])
        m = np.abs(r - R) <= sp['ring_half']
        th = np.arctan2(P[m, 1], P[m, 0])
        n = int(m.sum())
        Z = np.exp(8j * th).sum() if n else 0j
        out[j] = (np.mod(np.angle(Z) / 8.0, np.pi / 4), abs(Z) / max(n, 1), float(np.exp(-abs(Z) ** 2 / max(n, 1))), n)
    return out


def disc_inside(centres, R, margin, region):
    x0, y0, x1, y1 = region
    d = R + margin
    return (centres[:, 0] - d >= x0) & (centres[:, 0] + d <= x1) & (centres[:, 1] - d >= y0) & (centres[:, 1] + d <= y1)


def block_split(tab, ref_ranges, parity):
    """The T_ref localizations of the odd (parity 1) or even (0) 100-frame blocks of the reference frames, counted along the reference frames in order."""
    fr = tab['frame']
    cum = np.zeros(len(fr), np.int64)
    inside = np.zeros(len(fr), bool)
    off = 0
    for a, b in ref_ranges:
        m = (fr >= a) & (fr <= b)
        cum[m] = off + (fr[m] - a)
        inside |= m
        off += b - a + 1
    blk = (cum // 100) % 2
    return {k: v[inside & (blk == parity)] for k, v in tab.items()}


def step_pores(raw, log=print):
    """Step 2. Everything here uses T_ref only; no localizer output exists."""
    gates = {}
    for ds in ('D1', 'D2'):
        sp = DS[ds]
        fields = [None] if ds == 'D1' else [sp['fields']['dev'], sp['fields']['test']]
        for field in fields:
            tag = ds if field is None else '%s_f%d' % (ds, field)
            tab = load_table(ds, field)
            ref = table_select(tab, sp['ref_frames'])
            if ds == 'D1':
                region = (0.0, 0.0, sp['shape'][1] * sp['pixel_nm'], sp['shape'][0] * sp['pixel_nm'])
            else:
                region = (0.0, 0.0, sp['hist_size_nm'], sp['hist_size_nm'])
            t0 = time.time()
            cen, rad, cnt = find_pores(ref['x'], ref['y'], region, sp)
            R_ref = float(np.median(rad))
            log('%s: T_ref %d localizations; %d accepted pores; R_ref %.2f nm (%.0f s)' % (tag, len(ref['x']), len(cen), R_ref, time.time() - t0))
            g = dict(n_ref=int(len(ref['x'])), n_pores=int(len(cen)), R_ref=R_ref)
            for par, nm in ((1, 'odd'), (0, 'even')):
                h = block_split(ref, sp['ref_frames'], par)
                g['R_' + nm] = pore_radius(h['x'], h['y'], cen, sp)
            g['V1_abs_diff'] = abs(g['R_odd'] - g['R_even'])
            if sp['eight_fold']:
                ph = pore_phases(ref['x'], ref['y'], cen, R_ref, sp)
                g['frac_8fold_significant'] = float(np.mean(ph[:, 2] < 0.05))
                raw['pores|%s|phase' % tag] = ph
            raw['pores|%s|centres' % tag] = cen
            raw['pores|%s|radii' % tag] = rad
            raw['pores|%s|counts' % tag] = cnt
            gates[tag] = g
    # the windows (D2) and the scored pores
    sp = DS['D2']
    wpx = sp['window_px']
    for role in ('dev', 'test'):
        field = sp['fields'][role]
        tag = 'D2_f%d' % field
        cen = raw['pores|%s|centres' % tag]
        R = gates[tag]['R_ref']
        best, bw = -1, None
        step = 10
        for x0 in range(0, 419 - wpx + 1, step):
            for y0 in range(0, 419 - wpx + 1, step):
                reg = (x0 * sp['pixel_nm'], y0 * sp['pixel_nm'], (x0 + wpx) * sp['pixel_nm'], (y0 + wpx) * sp['pixel_nm'])
                n = int(disc_inside(cen, R, sp['disc_margin'], reg).sum())
                if n > best:
                    best, bw = n, (x0, y0)
        raw['window|D2%s' % role] = np.array(bw, dtype=np.int64)
        gates[tag]['window_px'] = list(bw)
        gates[tag]['pores_in_window'] = best
        log('D2 %s (field %d): window origin %s px holds %d accepted pores' % (role, field, bw, best))
    gates['V1'] = {k: gates[k]['V1_abs_diff'] <= 2.0 for k in gates if 'V1_abs_diff' in gates[k]}
    gates['V2'] = dict(diff=abs(gates['D2_f17']['R_ref'] - gates['D2_f18']['R_ref']), lit_dev=gates['D2_f17']['R_ref'] / R_LITERATURE_D2 - 1.0,
                       lit_test=gates['D2_f18']['R_ref'] / R_LITERATURE_D2 - 1.0)
    gates['V2']['pass'] = bool(gates['V2']['diff'] <= 3.0 and abs(gates['V2']['lit_dev']) <= 0.05 and abs(gates['V2']['lit_test']) <= 0.05)
    gates['V3'] = {k: gates[k]['frac_8fold_significant'] for k in gates if k.startswith('D2_f') and 'frac_8fold_significant' in gates[k]}
    raw['gates'] = np.array(json.dumps(gates, sort_keys=True))
    return gates


def scored_pores(raw, set_name):
    """The pores scored in a set: (centres (n, 2) nm in the table's coordinates, R_ref, region (x0, y0, x1, y1), phases (n,) or None, indices into the field's list)."""
    s = SETS[set_name]
    sp = DS[s['ds']]
    if s['ds'] == 'D1':
        tag = 'D1'
        cen = raw['pores|D1|centres']
        R = json.loads(str(raw['gates']))['D1']['R_ref']
        region = (0.0, 0.0, sp['shape'][1] * sp['pixel_nm'], sp['shape'][0] * sp['pixel_nm'])
        m = disc_inside(cen, R, sp['disc_margin'], region)
        m &= (cen[:, 0] < sp['x_dev_max']) if s['role'] == 'dev' else (cen[:, 0] > sp['x_test_min'])
        return cen[m], R, region, None, np.nonzero(m)[0]
    field = sp['fields'][s['role']]
    tag = 'D2_f%d' % field
    cen = raw['pores|%s|centres' % tag]
    R = json.loads(str(raw['gates']))[tag]['R_ref']
    x0, y0 = (int(v) for v in raw['window|D2' + s['role']])
    n = sp['window_px']
    region = (x0 * sp['pixel_nm'], y0 * sp['pixel_nm'], (x0 + n) * sp['pixel_nm'], (y0 + n) * sp['pixel_nm'])
    m = disc_inside(cen, R, sp['disc_margin'], region)
    return cen[m], R, region, raw['pores|%s|phase' % tag][m, 0], np.nonzero(m)[0]


def finish_gates(raw, gates):
    for sname in SETS:
        cen = scored_pores(raw, sname)[0]
        gates['V4_' + sname] = dict(n=int(len(cen)), ok=bool(len(cen) >= DS[SETS[sname]['ds']]['n_pores_min']))
    raw['gates'] = np.array(json.dumps(gates, sort_keys=True))
    return gates


def pore_list_hash(raw):
    h = hashlib.sha256()
    for k in sorted(raw):
        if k.startswith('pores|') or k.startswith('window|'):
            h.update(k.encode())
            h.update(np.ascontiguousarray(raw[k]).tobytes())
    return h.hexdigest()


# ==================================================================================================================================
# step 3: label-free calibration, the coordinate convention, the Picasso scale factors
# ==================================================================================================================================
N1_LOW_NG = 300.0                    # the lowest net gradient of the N1 grids of exp_smlm_picasso.FILTER_GRIDS
SIGMA_REF_NM = 140.0                 # the PSF width of the simulated sets: the width windows are scaled by sigma0 / 140 nm


def picasso_params(ds, **kw):
    sp = DS[ds]
    p = dict(pixel_nm=sp['pixel_nm'], **sp['camera'])
    p.update(kw)
    return p


def window_of(raw, set_name):
    s = SETS[set_name]
    return tuple(int(v) for v in raw['window|' + set_name[:2] + s['role']]) if s['ds'] == 'D2' else None


def origin_nm(raw, set_name):
    w = window_of(raw, set_name)
    return (0.0, 0.0) if w is None else (w[0] * DS['D2']['pixel_nm'], w[1] * DS['D2']['pixel_nm'])


def convention_offset(dets, tab, frame_range, region, origin, tol=250.0):
    """The median (Picasso - table) offset (nm, x and y) of the nearest-neighbour pairs (same frame, closer than tol) between Picasso R0 localizations (corner convention,
    window-relative; origin = the window's origin in nm) and the published table's localizations inside the window, frames frame_range (1-based, inclusive)."""
    d = []
    a, b = frame_range
    for i, dd in enumerate(dets):
        f = a + i
        m = (tab['frame'] == f) & (tab['x'] >= region[0]) & (tab['x'] < region[2]) & (tab['y'] >= region[1]) & (tab['y'] < region[3])
        if not m.any() or len(dd) == 0:
            continue
        t = np.c_[tab['x'][m], tab['y'][m]]
        p = dd[:, :2] + np.array(origin)
        dist, j = cKDTree(t).query(p)
        k = dist < tol
        d.append(p[k] - t[j[k]])
    d = np.concatenate(d)
    return np.median(d, axis=0), len(d)


def step_cal(raw, log=print):
    import exp_smlm_iss as ISS
    import exp_smlm_real as REAL
    # the N1 reference of the scale factors: Picasso at the lowest net gradient of the N1 grid on the development frames of MT0.N1.LD (defaults of the harness)
    fr = C.load_frames(C.data_path('MT0.N1.LD'), 0, 1000)
    d = C.run_picasso(fr, dict(min_net_gradient=N1_LOW_NG))
    allx = np.concatenate([x for x in d if len(x)])
    raw['cal|N1|ng_median'] = np.array(float(np.median(allx[:, C.NG_COL])))
    raw['cal|N1|ph_median'] = np.array(float(np.median(allx[:, C.PH_COL])))
    log('N1 reference (MT0.N1.LD frames 0-999, min net gradient %g): %d localizations, median net gradient %.0f, median photons %.0f' % (
        N1_LOW_NG, len(allx), np.median(allx[:, C.NG_COL]), np.median(allx[:, C.PH_COL])))
    for sname in SETS:
        sp = DS[SETS[sname]['ds']]
        win = window_of(raw, sname)
        cal_frames = calibration_frames(sname, win)
        t0 = time.time()
        r0 = C.run_picasso(cal_frames, picasso_params(SETS[sname]['ds']))
        lo = C.run_picasso(cal_frames, picasso_params(SETS[sname]['ds'], min_net_gradient=N1_LOW_NG))
        cal = REAL.calibrate(cal_frames, r0)
        law, sig = cal['law'], cal['sig']
        pix = sp['pixel_nm']
        log('%s: %d calibration frames of %d x %d px; Picasso R0 %d localizations; sigma0 = %.3f px (%.1f nm; %d selected%s); noise law a = %.3f b0 = %.1f x_ref %.1f (scale %.3f) (%.0f s)' % (
            sname, len(cal_frames), cal_frames.shape[1], cal_frames.shape[2], sum(len(x) for x in r0), sig['sigma0'], sig['sigma0'] * pix, sig['n_sel'],
            ', FALLBACK' if sig['fallback'] else '', law['a'], law['b0'], law['x_ref'], law['scale'], time.time() - t0))
        raw['cal|%s|bg' % sname] = cal['bg']
        raw['cal|%s|law' % sname] = np.array([law['a'], law['b0'], law['x_ref'], law['c'], law['scale'], law['a_hf'], law['b0_hf'], law['v_dark'], law['x_dark']])
        raw['cal|%s|sigma0' % sname] = np.array(json.dumps(sig))
        a_ng = np.concatenate([x for x in lo if len(x)])
        raw['cal|%s|low_run' % sname] = np.array([float(np.median(a_ng[:, C.NG_COL])), float(np.median(a_ng[:, C.PH_COL])), float(len(a_ng))])
        raw['cal|%s|r0_per_frame' % sname] = np.array(sum(len(x) for x in r0) / len(cal_frames))
        if SETS[sname]['role'] == 'dev':
            tab = load_table(SETS[sname]['ds'], sp['fields']['dev'] if SETS[sname]['ds'] == 'D2' else None)
            a, b = sp['dev_frames']
            region = scored_pores(raw, sname)[2] if SETS[sname]['ds'] == 'D2' else (0.0, 0.0, sp['shape'][1] * pix, sp['shape'][0] * pix)
            off, n = convention_offset(r0[:b - a + 1] if SETS[sname]['ds'] == 'D1' else r0[:b - a + 1], tab, (a, b), region, origin_nm(raw, sname))
            raw['cal|%s|convention_offset' % SETS[sname]['ds']] = np.array([off[0], off[1], n])
            log('%s: coordinate convention, median (Picasso R0 - published table) of %d pairs on frames %d-%d: dx %.1f dy %.1f nm' % (sname, n, a, b, off[0], off[1]))
    return raw


# ==================================================================================================================================
# the scoring of a localization set against the pore model (plan, section 4)
# ==================================================================================================================================
def convention(raw, ds):
    """The coordinate convention of the published table against the corner convention of the localizers (nm, added to the table's coordinates): the multiple of
    half a pixel nearest to the median (Picasso R0 - table) offset measured on the development frames (D1: 0, the measured -4 / -5 nm is a method difference; D2:
    -53 nm = half a pixel, SMAP counts from the pixel centre)."""
    off = raw['cal|%s|convention_offset' % ds][:2]
    half = DS[ds]['pixel_nm'] / 2.0
    return np.round(off / half) * half


class Scorer:
    """The pore model of one set: the scored pores (disc inside the scored region after the convention shift), the nearest-pore assignment, the radial bins and the
    area of the scored pore cells per bin (a grid of 10 nm), the 8-fold phases (D2)."""

    def __init__(self, raw, set_name):
        s = SETS[set_name]
        self.set_name, self.ds = set_name, s['ds']
        sp = DS[s['ds']]
        self.sp = sp
        cen_s, R, region, phases, idx = scored_pores(raw, set_name)
        tag = 'D1' if s['ds'] == 'D1' else 'D2_f%d' % sp['fields'][s['role']]
        self.conv = convention(raw, s['ds'])
        self.allc = raw['pores|%s|centres' % tag] + self.conv
        self.R, self.region = R, region
        sc = np.zeros(len(self.allc), bool)
        sc[idx] = True
        sc &= disc_inside(self.allc, R, sp['disc_margin'], region)
        self.scored = sc
        self.phase_all = raw['pores|%s|phase' % tag][:, 0] if sp['eight_fold'] else None
        self.rmax = R + sp['disc_margin']
        self.nb = int(np.ceil(self.rmax / RADIAL_BIN_NM))
        self.edges = np.arange(self.nb + 1) * RADIAL_BIN_NM
        self.tree = cKDTree(self.allc)
        pitch = 10.0
        gx = np.arange(region[0] + pitch / 2, region[2], pitch)
        gy = np.arange(region[1] + pitch / 2, region[3], pitch)
        G = np.stack(np.meshgrid(gx, gy), axis=-1).reshape(-1, 2)
        r, j = self.tree.query(G, distance_upper_bound=self.rmax)
        ok = np.isfinite(r) & (j < len(self.allc))
        ok[ok] &= sc[j[ok]]
        b = np.minimum((r[ok] / RADIAL_BIN_NM).astype(int), self.nb - 1)
        self.area = np.bincount(b, minlength=self.nb) * pitch * pitch
        self.scored_idx = np.nonzero(sc)[0]
        pos = np.full(len(self.allc), -1, np.int64)
        pos[self.scored_idx] = np.arange(len(self.scored_idx))
        self.pos = pos
        ap = np.zeros((len(self.scored_idx), self.nb))
        np.add.at(ap, (pos[j[ok]], b), pitch * pitch)
        self.area_pore = ap
        self.n_pores = int(sc.sum())

    def assign(self, xy):
        """(r, pore index j, keep): distance to the nearest pore, its index, and whether the localization is in the disc of a scored pore."""
        if len(xy) == 0:
            return np.zeros(0), np.zeros(0, np.int64), np.zeros(0, bool)
        r, j = self.tree.query(xy, distance_upper_bound=self.rmax)
        ok = np.isfinite(r) & (j < len(self.allc))
        keep = ok.copy()
        keep[ok] &= self.scored[j[ok]]
        return r, np.where(ok, j, 0), keep

    def hist(self, r, keep, block=None, nblocks=1):
        """Radial counts per block: (nblocks, nb)."""
        b = np.minimum((r[keep] / RADIAL_BIN_NM).astype(int), self.nb - 1)
        k = np.zeros(keep.sum(), np.int64) if block is None else np.asarray(block)[keep]
        H = np.zeros((nblocks, self.nb))
        np.add.at(H, (k, b), 1.0)
        return H


def fit_profile(counts, area, edges, R, w0=20.0, start=None):
    """Poisson maximum-likelihood fit of the pooled radial counts by N Normal(r; mu, w) + b x area: returns dict(N, mu, w, b, floor (= b x total area), n (counts), ok)."""
    counts = np.asarray(counts, float)
    lo, hi = edges[:-1], edges[1:]
    n = counts.sum()
    outer = slice(len(counts) * 2 // 3, None)
    b0 = max(counts[outer].sum() / max(area[outer].sum(), 1.0), 1e-9)
    N0 = max(n - b0 * area.sum(), 1.0)

    def nll(p):
        N, b, mu, w = np.exp(p[0]), np.exp(p[1]), p[2], np.exp(p[3])
        e = N * (ndtr((hi - mu) / w) - ndtr((lo - mu) / w)) + b * area
        e = np.maximum(e, 1e-12)
        return float(np.sum(e - counts * np.log(e)))
    x0 = np.array([np.log(N0), np.log(b0), R, np.log(w0)])
    bounds = [(np.log(1e-3), np.log(1e9)), (np.log(1e-12), np.log(1.0)), (R - 40.0, R + 40.0), (np.log(2.0), np.log(80.0))]
    best = None
    if start is not None:                                   # a warm start (the bootstrap and the thinnings): one run from the given parameters
        best = minimize(nll, np.asarray(start, float), method='L-BFGS-B', bounds=bounds)
    else:
        for w_start in (w0, 10.0, 35.0):
            x0[3] = np.log(w_start)
            r = minimize(nll, x0, method='L-BFGS-B', bounds=bounds)
            if best is None or r.fun < best.fun:
                best = r
    N, b, mu, w = np.exp(best.x[0]), np.exp(best.x[1]), best.x[2], np.exp(best.x[3])
    return dict(N=float(N), mu=float(mu), w=float(w), b=float(b), floor=float(b * area.sum()), n=float(n), ok=bool(best.success), x=[float(v) for v in best.x])


def e_score(fit, n_exp, alpha=ALPHA):
    """J_s and E_s of a profile fit (plan, section 4)."""
    miss = max(0.0, n_exp - fit['N'])
    J = fit['N'] / max(fit['N'] + fit['floor'] + miss, 1e-12)
    E = 100.0 - np.sqrt((100.0 - 100.0 * J) ** 2 + (alpha * fit['w']) ** 2)
    return float(J), float(E)


def expected_yield(raw, sc, n_frames_scored):
    """n_exp = A_ref x (scored frames / T_ref frames); A_ref = N of the profile fit of the T_ref localizations (the same pores, the convention shift applied)."""
    ds = sc.ds
    sp = sc.sp
    s = SETS[sc.set_name]
    field = sp['fields'][s['role']] if ds == 'D2' else None
    tab = table_select(load_table(ds, field), sp['ref_frames'])
    xy = np.c_[tab['x'], tab['y']] + sc.conv
    r, j, keep = sc.assign(xy)
    f = fit_profile(sc.hist(r, keep)[0], sc.area, sc.edges, sc.R)
    n_ref = sum(b - a + 1 for a, b in sp['ref_frames'])
    return f['N'] * n_frames_scored / n_ref, f


class Loc:
    """A localization set of the scored frames: frame (0-based within the set), x, y (nm, corner convention, absolute in the field), optional columns."""

    def __init__(self, frame, xy, extra=None):
        self.frame, self.xy = np.asarray(frame, np.int64), np.asarray(xy, float).reshape(-1, 2)
        self.extra = extra if extra is not None else np.zeros((len(self.frame), 0))

    @staticmethod
    def from_dets(dets, origin):
        """per-frame arrays whose first two columns are x_nm, y_nm (window-relative, corner convention) -> Loc (absolute)."""
        fr = np.concatenate([np.full(len(d), i, np.int64) for i, d in enumerate(dets)]) if len(dets) else np.zeros(0, np.int64)
        arr = np.concatenate([np.asarray(d, float).reshape(len(d), -1) for d in dets if len(d)]) if any(len(d) for d in dets) else np.zeros((0, 2))
        return Loc(fr, arr[:, :2] + np.array(origin), arr[:, 2:])

    def subset(self, m):
        return Loc(self.frame[m], self.xy[m], self.extra[m])


def score_loc(sc, loc, n_frames, n_exp, nblocks=1):
    """Profile fit, J_s, E_s of a Loc over all its frames."""
    r, j, keep = sc.assign(loc.xy)
    H = sc.hist(r, keep)
    f = fit_profile(H[0], sc.area, sc.edges, sc.R)
    J, E = e_score(f, n_exp)
    return dict(fit=f, J=J, E=E, n_loc=int(keep.sum()))


# ==================================================================================================================================
# step 4a: Picasso on the development frames (R0, R1E, R2E), scored for E_s
# ==================================================================================================================================
def ascent(ev, grids, max_rounds=3):
    """The coordinate ascent of `smlm_common.tune_postfilters` for an arbitrary scorer ev(filter dict) -> value (None = off; strict improvement; <= max_rounds rounds)."""
    keys = list(grids)
    cur = {k: None for k in keys}
    best = ev(cur)
    rounds = 0
    for rnd in range(max_rounds):
        rounds = rnd + 1
        changed = False
        for k in keys:
            for v in grids[k]:
                if v == cur[k]:
                    continue
                q = dict(cur)
                q[k] = v
                x = ev(q)
                if x > best + 1e-12:
                    best, cur, changed = x, q, True
        if not changed:
            break
    return cur, best, rounds


class PicCand:
    """One low-threshold Picasso run flattened (frame, x, y (nm, absolute), photons, sx, sy, net gradient) with the pore assignment done once."""

    def __init__(self, dets, origin, sc):
        fr = np.concatenate([np.full(len(d), i, np.int64) for i, d in enumerate(dets)]) if len(dets) else np.zeros(0, np.int64)
        arr = np.concatenate([np.asarray(d, float) for d in dets if len(d)]) if any(len(d) for d in dets) else np.zeros((0, 6))
        ok = np.isfinite(arr[:, :2]).all(axis=1) if len(arr) else np.zeros(0, bool)
        self.frame, self.arr = fr[ok], arr[ok]
        self.xy = self.arr[:, :2] + np.array(origin)
        self.r, self.j, self.keep = sc.assign(self.xy)
        self.sc = sc

    def mask(self, f):
        a = self.arr
        m = np.ones(len(a), bool)
        if f.get('min_net_gradient') is not None:
            m &= a[:, C.NG_COL] > f['min_net_gradient']
        if f.get('min_photons') is not None:
            m &= a[:, C.PH_COL] > f['min_photons']
        sx, sy = np.abs(a[:, C.SX_COL]), np.abs(a[:, C.SY_COL])
        if f.get('min_sigma_nm') is not None:
            m &= (sx > f['min_sigma_nm']) & (sy > f['min_sigma_nm'])
        if f.get('max_sigma_nm') is not None:
            m &= (sx < f['max_sigma_nm']) & (sy < f['max_sigma_nm'])
        return m

    def e_value(self, f, n_exp):
        m = self.mask(f) & self.keep
        if m.sum() < 30:
            return -np.inf
        b = np.minimum((self.r[m] / RADIAL_BIN_NM).astype(int), self.sc.nb - 1)
        h = np.bincount(b, minlength=self.sc.nb).astype(float)
        fit = fit_profile(h, self.sc.area, self.sc.edges, self.sc.R)
        return e_score(fit, n_exp)[1]


def scale_factors(raw, ds):
    """Picasso grid scale factors of a dataset from its DEVELOPMENT frames only (the plan, section 5)."""
    dev = ds + 'dev'
    low = raw['cal|%s|low_run' % dev]
    sig = json.loads(str(raw['cal|%s|sigma0' % dev]))
    return dict(f_ng=float(low[0] / raw['cal|N1|ng_median']), f_ph=float(low[1] / raw['cal|N1|ph_median']), f_sig=float(sig['sigma0'] * DS[ds]['pixel_nm'] / SIGMA_REF_NM))


def scaled_grids(fac):
    import exp_smlm_picasso as PIC
    g = PIC.FILTER_GRIDS['N1']
    sc = {'min_net_gradient': fac['f_ng'], 'min_photons': fac['f_ph'], 'min_sigma_nm': fac['f_sig'], 'max_sigma_nm': fac['f_sig']}
    grids = {k: [None if v is None else float(v) * sc[k] for v in vals] for k, vals in g.items()}
    r1 = [float(v) * fac['f_ng'] for v in PIC.GRID_R1['N1']]
    return grids, r1


def picasso_dev(raw, ds, log=print):
    """R0, R1E and R2E on the development frames of dataset ds for E_s; returns dict of the frozen configurations and the development scores."""
    import exp_smlm_picasso as PIC
    dev = ds + 'dev'
    sc = Scorer(raw, dev)
    frames = load_frames(dev, window_of(raw, dev))
    n = len(frames)
    n_exp, fref = expected_yield(raw, sc, n)
    org = origin_nm(raw, dev)
    fac = scale_factors(raw, ds)
    grids, r1 = scaled_grids(fac)
    log('%s dev: %d frames, %d scored pores, n_exp %.0f ring localizations (T_ref fit: N %.0f, mu %.1f, w %.1f); scale factors %s' % (ds, n, sc.n_pores, n_exp, fref['N'], fref['mu'], fref['w'], fac))
    base = picasso_params(ds)
    out = dict(fac=fac, n_exp=n_exp, n_frames=n)
    sess = C.PicassoSession(frames)
    try:
        d0 = sess.run(base)
        c0 = PicCand(d0, org, sc)
        out['R0'] = dict(params=base, E=c0.e_value({}, n_exp), n_loc=int(c0.keep.sum()))
        cand = PicCand(sess.run(dict(base, min_net_gradient=r1[0])), org, sc)
        vals = [(cand.e_value({'min_net_gradient': v}, n_exp), v) for v in r1]
        e1, v1 = max(vals, key=lambda t: (round(t[0], 9), -t[1]))
        out['R1E'] = dict(params=dict(base, min_net_gradient=v1), E=e1, grid=[list(t) for t in vals])
        log('%s dev: R0 E_s %.2f; R1E E_s %.2f at net gradient %.0f' % (ds, out['R0']['E'], e1, v1))
        best, n_cfg, t_all = None, 0, time.time()
        structs = PIC.struct_list()
        res = []
        for i, st in enumerate(structs):
            cand = PicCand(sess.run(dict(base, **dict(st), min_net_gradient=grids['min_net_gradient'][0])), org, sc)
            cache = {}

            def ev(f, cand=cand, cache=cache):
                key = tuple(f.get(k) for k in grids)
                if key not in cache:
                    cache[key] = cand.e_value(f, n_exp)
                return cache[key]
            f, val, rounds = ascent(ev, grids, PIC.MAX_ROUNDS)
            n_cfg += len(cache)
            res.append((val, i, dict(st), dict(f), rounds))
            if i % 12 == 0:
                log('  [%d/%d] %s: dev E_s %.2f (rounds %d), %d configurations so far, %.0f s' % (i + 1, len(structs), PIC.skey(st), val, rounds, n_cfg, time.time() - t_all))
        # the best structure by development E_s; ties: the earlier in the grid order
        val, i, st, f, rounds = max(res, key=lambda t: (round(t[0], 9), -t[1]))
        out['R2E'] = dict(params=dict(base, **st, **{k: v for k, v in f.items()}), E=val, struct=st, filters=f, rounds=rounds, n_structs=len(structs), n_configs=n_cfg)
        log('%s dev: R2E E_s %.2f with %s %s (%d configurations scored)' % (ds, val, PIC.skey(st), f, n_cfg))
    finally:
        sess.close()
    return out


# ==================================================================================================================================
# step 4b: the ISS (dual order, guarded) and the greedy order on the development frames: 168 operating points each
# ==================================================================================================================================
BASE_CFGS = ('hp3', 'hp3h', 'hp5', 'hp5h')
VARIANTS = ('_s12', '_s17', '_r15', '_r25')
KAPPAS = (4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.0)
SMAX_REL = (2.0, 2.5, 3.2)                       # px at sigma0 = 1.4 px; scaled by sigma0 / 1.4; the last is the bound of the widths (no filter)
DEFAULT_CFG = 'hp5h'                             # the label-free a priori configuration of real_frc.md (kappa 5, no width filter)
SNR_MIN_PATH, SNR_EXTRA_PATH = 4.0, 0            # the shortened stop of the controls (the detections of every kappa >= 4 are those of the long path)


def grid_cfg(name, sigma0):
    """The ISS configuration `name` (a base of ISS.CONFIGS, optionally with a suffix of the round-2 variants) for a label-free PSF width sigma0 (px)."""
    import exp_smlm_iss as ISS
    import exp_smlm_real as REAL
    base = name[:name.index('_')] if '_' in name else name
    suf = name[len(base):]
    s = float(sigma0) / ISS.SIGMA0
    cfg = REAL.iss_config(sigma0)
    b = ISS.CONFIGS[base]
    cfg.update(cfg=name, hp=float(b['hp']), parts=[list(p) for p in b['parts']])
    if suf == '_s12':
        cfg['sigma0'] = float(sigma0) * 1.2 / ISS.SIGMA0
    elif suf == '_s17':
        cfg['sigma0'] = float(sigma0) * 1.7 / ISS.SIGMA0
    elif suf == '_r15':
        cfg['r_ex'] = 1.5 * s
    elif suf == '_r25':
        cfg['r_ex'] = 2.5 * s
    elif suf:
        raise ValueError(name)
    return cfg


def iss_cal(raw, set_name):
    """The label-free calibration of a set (bg surface, noise law a, b0, x_ref; sigma0 px) as stored by step 3."""
    law = raw['cal|%s|law' % set_name]
    sig = json.loads(str(raw['cal|%s|sigma0' % set_name]))
    return dict(bg=raw['cal|%s|bg' % set_name], a=float(law[0]), b0=float(law[1]), x_ref=float(law[2])), float(sig['sigma0'])


def iss_frame(frame, ctx, order):
    """One frame (ADU): the path of the sliding ISS (order 'dual', through the cycle guard `imaging_guard.sliding_path_guarded`) or of the greedy order, computed once with
    the shortened stop, and for every kappa the block of `stop_snr` (no final polish, as in real_frc.md). Returns (dict kappa -> (k, 4) array x_nm, y_nm, w_adu, sigma_px),
    n entries, seconds)."""
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
                                rtrunc=ctx['rtrunc'], r_ex=cfg['r_ex'], fit_radius=cfg['fit_radius'], order=order, snr_min=SNR_MIN_PATH, snr_extra=SNR_EXTRA_PATH, cycle_guard=True)
    out = {}
    for kap in KAPPAS:
        i = IM.stop_snr(R['hist'], R['taus'], R['snr'], kap)[0]
        h = R['hist'][i]
        c, w = h['c'], h['w']
        pix = ctx['pixel_nm']
        out[kap] = np.c_[(c[:, 0] + 0.5) * pix, (c[:, 1] + 0.5) * pix, w, c[:, 2]] if len(w) else np.zeros((0, 4))
    return out, len(R['taus']), time.perf_counter() - t0


def iss_run(frames, cal, cfg, order, pixel_nm, log=print, tag=''):
    """The ISS / greedy path on every frame: returns (dets[kappa] = list of per-frame arrays, seconds per frame (n,), entries (n,))."""
    import exp_smlm_real as REAL
    ctx = REAL.make_context(cal['bg'], cal['a'], cal['b0'], cal['x_ref'], cfg, frames.shape[1:])
    ctx['pixel_nm'] = pixel_nm
    dets = {k: [] for k in KAPPAS}
    secs, ent = [], []
    t_start = time.time()
    for i, fr in enumerate(frames):
        d, n, sec = iss_frame(fr, ctx, order)
        for k in KAPPAS:
            dets[k].append(d[k])
        secs.append(sec)
        ent.append(n)
        if (i + 1) % 50 == 0:
            log('    %s %s %s: %d / %d frames, %.2f s per frame, %.0f s' % (tag, cfg['cfg'], order, i + 1, len(frames), float(np.mean(secs)), time.time() - t_start))
    return dets, np.array(secs), np.array(ent)


def post_width(dets, smax):
    return [d[d[:, 3] <= smax + 1e-9] if len(d) else d for d in dets]


def op_rows(dets_by_kappa, sc, origin, n_exp, s_scale, widths_hi):
    """E_s (and J_s, N) of the 7 x 3 operating points of one path: [(kappa, smax_rel, E, J, n_ring)]."""
    rows = []
    for kap in KAPPAS:
        for sm in SMAX_REL:
            smax = sm * s_scale
            dets = post_width(dets_by_kappa[kap], smax)
            loc = Loc.from_dets(dets, origin)
            r = score_loc(sc, loc, 0, n_exp)
            rows.append((kap, sm, r['E'], r['J'], r['fit']['N']))
    return rows


def choose_op(rows):
    """The operating point of the largest development E_s (ties: the larger kappa, then the smaller width bound)."""
    return max(rows, key=lambda r: (round(r[2], 6), r[0], -r[1]))


def iss_dev(raw, ds, order, log=print, cfgs=None, n_frames=None):
    """Development tuning of the ISS (order 'dual') or the greedy order: the four base configurations, then the four round-2 variants around the best base (the 8
    configurations x 7 kappas x 3 width bounds = 168 operating points of iss.md), for E_s. Returns the frozen choice and every point's score."""
    dev = ds + 'dev'
    sc = Scorer(raw, dev)
    frames = load_frames(dev, window_of(raw, dev))
    if n_frames:
        frames = frames[:n_frames]
    n = len(frames)
    n_exp, _ = expected_yield(raw, sc, n)
    org = origin_nm(raw, dev)
    cal, sigma0 = iss_cal(raw, dev)
    pix = DS[ds]['pixel_nm']
    s_scale = sigma0 / 1.4
    allrows, secs = {}, {}
    for stage, names in (('base', BASE_CFGS), ('variants', None)):
        if stage == 'variants':
            best_base = max(((choose_op(allrows[c]), c) for c in BASE_CFGS), key=lambda t: (round(t[0][2], 6), t[0][0], -t[0][1]))[1]
            names = [best_base + v for v in VARIANTS]
        for name in names:
            cfg = grid_cfg(name, sigma0)
            t0 = time.time()
            dets, sec, ent = iss_run(frames, cal, cfg, order, pix, log=log, tag=ds + 'dev')
            allrows[name] = op_rows(dets, sc, org, n_exp, s_scale, cfg['widths'][1])
            secs[name] = float(sec.mean())
            bo = choose_op(allrows[name])
            log('%s dev %s %s: %d frames, %.2f s per frame, best point kappa %g width %g: E_s %.2f J_s %.3f (%.0f s)' % (ds, order, name, n, sec.mean(), bo[0], bo[1], bo[2], bo[3], time.time() - t0))
    best = max(((choose_op(r), c) for c, r in allrows.items()), key=lambda t: (round(t[0][2], 6), t[0][0], -t[0][1]))
    (kap, sm, E, J, N), name = best
    dflt = [r for r in allrows[DEFAULT_CFG] if r[0] == 5.0 and r[1] == SMAX_REL[-1]][0]
    return dict(order=order, cfg=name, kappa=kap, smax_rel=sm, E=E, J=J, n_exp=n_exp, n_frames=n, default_E=dflt[2], rows={c: [list(map(float, r)) for r in v] for c, v in allrows.items()},
                secs=secs)


def jdump(o):
    return np.array(json.dumps(o, sort_keys=True, default=float))


def jload(a):
    return json.loads(str(a))


def freeze_hash(raw):
    """sha256 over the pore lists, the scale factors and the frozen development choices (everything the test run is allowed to depend on)."""
    h = hashlib.sha256()
    for k in sorted(raw):
        if k.startswith(('pores|', 'window|', 'cal|', 'dev|')):
            h.update(k.encode())
            h.update(np.ascontiguousarray(raw[k]).tobytes())
    return h.hexdigest()


def step_dev(raw, log=print):
    """Step 4: tune Picasso and the ISS / greedy rows on the development frames (checkpointed in the raw), then freeze and hash."""
    for ds in ('D1', 'D2'):
        if 'dev|%s|picasso' % ds not in raw:
            raw['dev|%s|picasso' % ds] = jdump(picasso_dev(raw, ds, log))
            save_raw(raw)
        for order in ('dual', 'greedy'):
            key = 'dev|%s|%s' % (ds, order)
            if key not in raw:
                raw[key] = jdump(iss_dev(raw, ds, order, log))
                save_raw(raw)
    raw['freeze|time_utc'] = np.array(time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()))
    raw['freeze|sha256'] = np.array(freeze_hash(raw))
    save_raw(raw)
    log('frozen %s sha256 %s' % (str(raw['freeze|time_utc']), str(raw['freeze|sha256'])))
    return raw


# ==================================================================================================================================
# step 5: the ONE test run (frozen configurations)
# ==================================================================================================================================
ROWS = ('R0', 'R1E', 'R2E', 'ISS', 'ISS0', 'GREEDY', 'PUB')
ROW_NAMES = {'R0': 'Picasso, defaults (R0)', 'R1E': 'Picasso, net gradient tuned for E_s (R1E)', 'R2E': 'Picasso, every option tuned for E_s (R2E)',
             'ISS': 'sliding ISS, tuned for E_s', 'ISS0': 'sliding ISS, label-free a priori configuration (hp5h, kappa 5, no width filter)',
             'GREEDY': 'greedy order (matching pursuit with sliding), tuned for E_s', 'PUB': 'as published by the authors (context only)'}


def put_loc(raw, key, loc, secs=None):
    raw[key + '|frame'] = loc.frame.astype(np.int32)
    raw[key + '|xy'] = loc.xy
    if secs is not None:
        raw[key + '|secs'] = np.asarray(secs, float)


def get_loc(raw, key):
    return Loc(raw[key + '|frame'].astype(np.int64), raw[key + '|xy'])


def table_loc(ds, set_name, conv):
    """The authors' table restricted to the frames of the set, as a Loc in the corner convention (frame 0-based within the set)."""
    s = SETS[set_name]
    sp = DS[ds]
    a, b = sp['dev_frames'] if s['role'] == 'dev' else sp['test_frames']
    tab = table_select(load_table(ds, sp['fields'][s['role']] if ds == 'D2' else None), [(a, b)])
    return Loc(tab['frame'] - a, np.c_[tab['x'], tab['y']] + conv)


def step_test(raw, ds, log=print):
    """Step 5 for dataset ds: every row once on the test frames with the frozen configurations (stored in the raw: `test|<ds>|<row>|...`)."""
    test = ds + 'test'
    sc = Scorer(raw, test)
    win = window_of(raw, test)
    frames = load_frames(test, win)
    n = len(frames)
    org = origin_nm(raw, test)
    pic, dual, greedy = jload(raw['dev|%s|picasso' % ds]), jload(raw['dev|%s|dual' % ds]), jload(raw['dev|%s|greedy' % ds])
    raw['test|%s|n_frames' % ds] = np.array(n)
    n_exp, fref = expected_yield(raw, sc, n)
    raw['test|%s|ref_fit' % ds] = jdump(dict(n_exp=n_exp, fit=fref))
    for row in ('R0', 'R1E', 'R2E'):
        info = {}
        t0 = time.time()
        dets = C.run_picasso(frames, pic[row]['params'], info=info)
        put_loc(raw, 'test|%s|%s' % (ds, row), Loc.from_dets(dets, org), secs=[info['run_s'] / n])
        log('%s test %s: %d localizations, %.4f s per frame (%.0f s)' % (ds, row, sum(len(d) for d in dets), info['run_s'] / n, time.time() - t0))
    cal, sigma0 = iss_cal(raw, test)
    pix = DS[ds]['pixel_nm']
    s_scale = sigma0 / 1.4
    cache = {}

    def path(order, name):
        if (order, name) not in cache:
            cfg = grid_cfg(name, sigma0)
            t0 = time.time()
            dets, secs, ent = iss_run(frames, cal, cfg, order, pix, log=log, tag=ds + 'test')
            cache[(order, name)] = (dets, secs, ent)
            log('%s test %s %s: %d frames, %.2f s per frame (%.0f s)' % (ds, order, name, n, secs.mean(), time.time() - t0))
        return cache[(order, name)]
    for row, order, ch in (('ISS', 'dual', dual), ('GREEDY', 'greedy', greedy), ('ISS0', 'dual', dict(cfg=DEFAULT_CFG, kappa=5.0, smax_rel=SMAX_REL[-1]))):
        dets, secs, ent = path(order, ch['cfg'])
        d = post_width(dets[float(ch['kappa'])], ch['smax_rel'] * s_scale)
        put_loc(raw, 'test|%s|%s' % (ds, row), Loc.from_dets(d, org), secs=secs)
        raw['test|%s|%s|entries' % (ds, row)] = ent
        raw['test|%s|%s|choice' % (ds, row)] = jdump(dict(cfg=ch['cfg'], kappa=float(ch['kappa']), smax_rel=float(ch['smax_rel']), sigma0=sigma0))
    put_loc(raw, 'test|%s|PUB' % ds, table_loc(ds, test, sc.conv))
    return raw


def step_test_all(raw, log=print):
    if 'freeze|sha256' not in raw or str(raw['freeze|sha256']) != freeze_hash(raw):
        raise RuntimeError('the frozen development choices are missing or changed: the test run is refused')
    raw['test|start_utc'] = np.array(time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()))
    raw['test|meta'] = jdump(dict(python=platform.python_version(), numpy=np.__version__, cpu_count=os.cpu_count(), loadavg_at_start=[round(v, 2) for v in os.getloadavg()]))
    save_raw(raw)
    for ds in ('D1', 'D2'):
        if 'test|%s|ISS|xy' % ds in raw:
            continue
        step_test(raw, ds, log)
        save_raw(raw)
    raw['test|end_utc'] = np.array(time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()))
    save_raw(raw)
    return raw


# ==================================================================================================================================
# the analysis of the test rows from the raw alone (plan sections 4, 6, 7)
# ==================================================================================================================================
def prep_row(sc, loc, nblocks, eight_fold):
    """Per-block radial histograms, per-pore histograms and the arrays of the secondary outcomes of one localization set."""
    r, j, keep = sc.assign(loc.xy)
    blk = np.minimum(loc.frame // BLOCK, nblocks - 1)
    b = np.minimum((np.where(keep, r, 0.0) / RADIAL_BIN_NM).astype(int), sc.nb - 1)
    H = np.zeros((nblocks, sc.nb))
    np.add.at(H, (blk[keep], b[keep]), 1.0)
    Hp = np.zeros((sc.n_pores, sc.nb))
    np.add.at(Hp, (sc.pos[j[keep]], b[keep]), 1.0)
    out = dict(H=H, Hp=Hp, n_loc=int(len(loc.frame)), n_keep=int(keep.sum()))
    ring = keep & (np.abs(r - sc.R) <= sc.sp['ring_half'])
    out['ring_per_pore'] = np.bincount(sc.pos[j[ring]], minlength=sc.n_pores)
    if eight_fold:
        c = sc.allc[j]
        th = np.arctan2(loc.xy[:, 1] - c[:, 1], loc.xy[:, 0] - c[:, 0])
        z = np.cos(8.0 * (th - sc.phase_all[j]))
        Zb = np.bincount(blk[ring], weights=z[ring], minlength=nblocks)
        nb_ = np.bincount(blk[ring], minlength=nblocks).astype(float)
        out['c8_sum'], out['c8_n'] = Zb, nb_
        rng = np.random.default_rng(0)
        nulls = []
        for _ in range(N_THIN):
            ph = rng.permutation(sc.phase_all)
            zz = np.cos(8.0 * (th - ph[j]))
            nulls.append(zz[ring].sum() / max(ring.sum(), 1))
        out['c8_null'] = float(np.mean(nulls))
    # duplicates: the share of localizations with another localization of the same frame closer than sigma0
    d = 0
    if len(loc.frame):
        t = cKDTree(np.c_[loc.xy, loc.frame * 1.0e7])
        pr = t.query_pairs(sc.sp['pixel_nm'] * sc.sigma0_px, output_type='ndarray') if hasattr(sc, 'sigma0_px') else np.zeros((0, 2), int)
        d = len(np.unique(pr)) if len(pr) else 0
    out['dup'] = d
    return out


def stats_of(sc, Hs, area, n_exp, start=None):
    f = fit_profile(Hs, area, sc.edges, sc.R, start=start)
    J, E = e_score(f, n_exp)
    return f, J, E


def eq_yield_w(sc, fa, Ha, fb, Hb, area, rng):
    """w of both sets after thinning the one with the larger ring yield (binomial thinning of its radial counts) to the other's N: mean over N_THIN thinnings."""
    if fa['N'] >= fb['N']:
        big, Hbig, small, swap = fa, Ha, fb, False
    else:
        big, Hbig, small, swap = fb, Hb, fa, True
    q = min(small['N'] / max(big['N'], 1e-9), 1.0)
    ws = []
    for _ in range(N_THIN):
        Ht = rng.binomial(Hbig.astype(np.int64), q).astype(float)
        st = list(big['x'])
        st[0] += np.log(max(q, 1e-9))
        ws.append(fit_profile(Ht, area, sc.edges, sc.R, start=st)['w'])
    wb = float(np.mean(ws))
    return (small['w'], wb) if swap else (wb, small['w'])


def percentile_ci(v, level=0.95):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    return [float(np.percentile(v, 100 * (1 - level) / 2)), float(np.percentile(v, 100 * (1 + level) / 2))]


REG_PIX_NM = 20.0
REG_BLOCK = 250
REG_SEARCH_NM = 800.0


def _reg_image(xy, extent):
    n = int(np.ceil(extent / REG_PIX_NM)) + 1
    H, _, _ = np.histogram2d(xy[:, 1], xy[:, 0], bins=[np.arange(n + 1) * REG_PIX_NM] * 2)
    return H


def block_shift(xy_a, xy_b, extent):
    """The translation (nm) that maps the localizations a onto b: the peak of the cross-correlation of their 20 nm histograms within +-800 nm, refined by a parabola
    through the peak and its neighbours; returns (dx, dy, peak over the standard deviation of the correlation)."""
    a, b = _reg_image(xy_a, extent), _reg_image(xy_b, extent)
    c = signal.fftconvolve(b - b.mean(), (a - a.mean())[::-1, ::-1], mode='same')
    cy, cx = np.array(c.shape) // 2
    r = int(REG_SEARCH_NM / REG_PIX_NM)
    sub = c[cy - r:cy + r + 1, cx - r:cx + r + 1]
    iy, ix = np.unravel_index(np.argmax(sub), sub.shape)
    off = []
    for k, i in ((0, iy), (1, ix)):
        if 0 < i < sub.shape[k] - 1:
            m1, p0, p1 = (sub[i - 1, ix], sub[i, ix], sub[i + 1, ix]) if k == 0 else (sub[iy, i - 1], sub[iy, i], sub[iy, i + 1])
            den = m1 - 2 * p0 + p1
            off.append(0.5 * (m1 - p1) / den if den != 0 else 0.0)
        else:
            off.append(0.0)
    dy, dx = (iy - r + off[0]) * REG_PIX_NM, (ix - r + off[1]) * REG_PIX_NM
    return float(dx), float(dy), float(sub.max() / max(c.std(), 1e-12))


def register_loc(loc, ref, n_frames, extent):
    """Post-hoc registration (Amendment 4): the localizations of every block of REG_BLOCK frames are translated by the shift that maps them onto the published table's
    localizations of the same frames; blocks with fewer than 50 localizations use the row's shift of all frames. Returns (Loc, shifts (nblocks, 3))."""
    nb = int(np.ceil(n_frames / REG_BLOCK))
    allsh = block_shift(loc.xy, ref.xy, extent)
    out = loc.xy.copy()
    sh = np.zeros((nb, 3))
    for k in range(nb):
        m = (loc.frame >= k * REG_BLOCK) & (loc.frame < (k + 1) * REG_BLOCK)
        mr = (ref.frame >= k * REG_BLOCK) & (ref.frame < (k + 1) * REG_BLOCK)
        sh[k] = block_shift(loc.xy[m], ref.xy[mr], extent) if m.sum() >= 50 and mr.sum() >= 50 else allsh
        out[m] += sh[k, :2]
    return Loc(loc.frame, out, loc.extra), sh


def analyse(raw, ds, n_boot=N_BOOT, log=lambda s: None, register=False):
    """Everything the report prints for one dataset: per-row point estimates and block-bootstrap intervals, the paired differences (primary: ISS - R2E; descriptive: ISS - greedy),
    and the pore-resampling sensitivity of the pairs."""
    test = ds + 'test'
    sc = Scorer(raw, test)
    sig = jload(raw['cal|%s|sigma0' % test])
    sc.sigma0_px = float(sig['sigma0'])
    n = int(raw['test|%s|n_frames' % ds])
    nblocks = int(np.ceil(n / BLOCK))
    ref = jload(raw['test|%s|ref_fit' % ds])
    n_exp = ref['n_exp']
    eight = DS[ds]['eight_fold']
    locs = {row: get_loc(raw, 'test|%s|%s' % (ds, row)) for row in ROWS}
    shifts = {}
    if register:
        extent = float(max(sc.region[2], sc.region[3])) + 2000.0
        for row in ROWS:
            if row != 'PUB':
                locs[row], shifts[row] = register_loc(locs[row], locs['PUB'], n, extent)
    P = {row: prep_row(sc, locs[row], nblocks, eight) for row in ROWS}
    out = dict(n_frames=n, nblocks=nblocks, n_pores=sc.n_pores, n_exp=n_exp, ref_fit=ref['fit'], R_ref=sc.R, rows={}, shifts={k: v.tolist() for k, v in shifts.items()})
    base = {}
    for row in ROWS:
        Hs = P[row]['H'].sum(0)
        f, J, E = stats_of(sc, Hs, sc.area, n_exp)
        base[row] = f
        cr = P[row]['ring_per_pore']
        d = dict(fit=f, J=J, E=E, n_loc=P[row]['n_loc'], n_keep=P[row]['n_keep'], yield_per_pore_kframe=float(f['N'] / sc.n_pores / (n / 1000.0)),
                 frac_pores_8=float(np.mean(cr >= 8)), dup_frac=P[row]['dup'] / max(P[row]['n_loc'], 1), mu_bias=f['mu'] - sc.R)
        if eight:
            d['c8'] = float(P[row]['c8_sum'].sum() / max(P[row]['c8_n'].sum(), 1.0))
            d['c8_null'] = P[row]['c8_null']
        sec = raw.get('test|%s|%s|secs' % (ds, row))
        d['secs_per_frame'] = float(np.mean(sec)) if sec is not None else float('nan')
        if 'test|%s|%s|entries' % (ds, row) in raw:
            e = raw['test|%s|%s|entries' % (ds, row)]
            d['mean_entries'], d['frac_at_cap'] = float(e.mean()), float(np.mean(e >= 400))
        out['rows'][row] = d
    pairs = (('ISS', 'R2E'), ('ISS', 'GREEDY'), ('ISS0', 'R2E'), ('R2E', 'R0'))
    # block bootstrap: the same resampled blocks for every row
    rng = np.random.default_rng(BOOT_SEED)
    rep = {row: dict(E=[], J=[], N=[], w=[], mu=[], c8=[]) for row in ROWS}
    dif = {pr: dict(dE=[], dJ=[], dw=[]) for pr in pairs}
    t0 = time.time()
    for it in range(n_boot):
        idx = rng.integers(0, nblocks, nblocks)
        fits = {}
        for row in ROWS:
            Hs = P[row]['H'][idx].sum(0)
            f, J, E = stats_of(sc, Hs, sc.area, n_exp, start=base[row]['x'])
            fits[row] = (f, J, E, Hs)
            r_ = rep[row]
            r_['E'].append(E), r_['J'].append(J), r_['N'].append(f['N']), r_['w'].append(f['w']), r_['mu'].append(f['mu'] - sc.R)
            if eight:
                r_['c8'].append(P[row]['c8_sum'][idx].sum() / max(P[row]['c8_n'][idx].sum(), 1.0))
        for pr in pairs:
            a, b = pr
            dif[pr]['dE'].append(fits[a][2] - fits[b][2])
            dif[pr]['dJ'].append(fits[a][1] - fits[b][1])
            if pr in pairs[:2]:                                     # the equal-yield width has an interval for the two pairs of Amendment 1 and the registered rule only
                wa, wb = eq_yield_w(sc, fits[a][0], fits[a][3], fits[b][0], fits[b][3], sc.area, rng)
                dif[pr]['dw'].append(wa - wb)
    log('%s block bootstrap %d replicates: %.0f s' % (ds, n_boot, time.time() - t0))
    for row in ROWS:
        out['rows'][row]['ci'] = {k: percentile_ci(v) for k, v in rep[row].items() if len(v)}
    out['pairs'] = {}
    for pr in pairs:
        a, b = pr
        wa, wb = eq_yield_w(sc, base[a], P[a]['H'].sum(0), base[b], P[b]['H'].sum(0), sc.area, np.random.default_rng(BOOT_SEED + 7))
        out['pairs']['%s-%s' % pr] = dict(dE=out['rows'][a]['E'] - out['rows'][b]['E'], dE_ci=percentile_ci(dif[pr]['dE']), dJ=out['rows'][a]['J'] - out['rows'][b]['J'],
                                           dJ_ci=percentile_ci(dif[pr]['dJ']), dw=wa - wb, dw_ci=percentile_ci(dif[pr]['dw']) if dif[pr]['dw'] else None, w_a=wa, w_b=wb)
    # sensitivity: resample the pores (clusters, all their frames); n_exp kept fixed
    rng = np.random.default_rng(BOOT_SEED_PORES)
    dsens = {pr: dict(dE=[], dw=[]) for pr in pairs[:2]}
    npor = sc.n_pores
    for it in range(n_boot):
        idx = rng.integers(0, npor, npor)
        area = sc.area_pore[idx].sum(0)
        fits = {}
        for row in ('ISS', 'R2E', 'GREEDY'):
            Hs = P[row]['Hp'][idx].sum(0)
            f, J, E = stats_of(sc, Hs, area, n_exp, start=base[row]['x'])
            fits[row] = (f, J, E, Hs)
        for pr in pairs[:2]:
            a, b = pr
            wa, wb = eq_yield_w(sc, fits[a][0], fits[a][3], fits[b][0], fits[b][3], area, rng)
            dsens[pr]['dE'].append(fits[a][2] - fits[b][2])
            dsens[pr]['dw'].append(wa - wb)
    for pr in pairs[:2]:
        out['pairs']['%s-%s' % pr]['dE_ci_pores'] = percentile_ci(dsens[pr]['dE'])
        out['pairs']['%s-%s' % pr]['dw_ci_pores'] = percentile_ci(dsens[pr]['dw'])
    return out


# ==================================================================================================================================
# the report (from the raw alone)
# ==================================================================================================================================
def md_table(header, rows):
    out = ['| ' + ' | '.join(header) + ' |', '|' + '|'.join('---' for _ in header) + '|']
    out += ['| ' + ' | '.join(str(c) for c in r) + ' |' for r in rows]
    return '\n'.join(out)


def ci_str(v, ci, nd=1):
    return '%s [%s, %s]' % (fmt(v, nd), fmt(ci[0], nd), fmt(ci[1], nd))


def verdict_ds(pr):
    """The reading of the primary outcomes of one dataset under the plan's section 7."""
    dE, dw = pr['dE_ci'], pr['dw_ci']
    better = dE[0] > 0 and dw[1] < 0
    better_pores = pr['dE_ci_pores'][0] > 0 and pr['dw_ci_pores'][1] < 0
    worse = dE[1] < 0
    equiv = (dE[0] > -MARGIN_E and dE[1] < MARGIN_E) and not (dE[0] > 0 or dE[1] < 0)
    if better and better_pores:
        return 'ISS better'
    if better and not better_pores:
        return 'ISS better in the frame bootstrap only (not in the pore-resampling sensitivity: not claimed)'
    if worse:
        return 'Picasso better'
    if dE[0] > 0 and not dw[1] < 0:
        return 'Delta E_s above 0 but the width at equal yield is not smaller: more localizations, not better ones (not an advantage)'
    if dw[1] < 0 and not dE[0] > 0:
        return 'narrower at equal yield, but Delta E_s does not exclude 0 (reported as such)'
    if equiv:
        return 'equivalent within +-%g E points' % MARGIN_E
    return 'no detectable difference (intervals include 0, wider than the margin)'


def overall_verdict(v1, v2, dev_gate_both):
    better = [v == 'ISS better' for v in (v1, v2)]
    worse = [v == 'Picasso better' for v in (v1, v2)]
    if all(better):
        txt = 'ISS better on BOTH datasets: a real-data advantage may be written, with the single-field and ring-model caveats of the plan (section 1).'
    elif any(worse) and not any(better):
        txt = 'DROP: Picasso better on at least one dataset and the ISS on neither; the real-data claim is dropped and the manuscript keeps only the FRC experiment, which makes no accuracy claim.'
    elif any(better) and any(worse):
        txt = 'MIXED: the ISS is better on one dataset and Picasso on the other; reported as such, no real-data advantage is claimed.'
    elif any(better):
        txt = 'MIXED: the ISS is better on one dataset only; the real-data sentence is limited to that dataset, the other result stated next to it.'
    else:
        txt = 'NO DETECTABLE DIFFERENCE on real data under the registered rule; the manuscript must not claim a real-data advantage.'
    if dev_gate_both:
        txt += ' The development gate of section 7 (iii) is also met (Picasso R2E exceeds the tuned ISS by more than %g E points on both development sets).' % GATE_DEV_E
    return txt


def build_md(raw, n_boot=N_BOOT, log=lambda s: None):
    gates = jload(raw['gates'])
    an = {ds: analyse(raw, ds, n_boot, log) for ds in ('D1', 'D2')}
    ar = {ds: analyse(raw, ds, n_boot, log, register=True) for ds in ('D1', 'D2')}
    dev = {ds: dict(pic=jload(raw['dev|%s|picasso' % ds]), dual=jload(raw['dev|%s|dual' % ds]), greedy=jload(raw['dev|%s|greedy' % ds])) for ds in ('D1', 'D2')}
    L = []
    w = L.append
    w('# Real nuclear-pore recordings with a known geometry: the sliding ISS against Picasso (pre-registered)')
    w('')
    w('Generated by `experiments/exp_smlm_real_gt.py --from-raw` from `results/smlm/real_gt_raw.npz`; every number below is computed from that raw (the file regenerates byte for byte). The specification is `results/smlm/real_gt_plan.md` (declared 2026-11-09 07:56 UTC; '
      'Amendments 1 and 2 declared before any pore finder or localizer ran on development or test frames). **The frozen development choices were hashed at %s (sha256 `%s`) before the '
      'test run; the test run started at %s and ended at %s; one test run per dataset.**' % (str(raw['freeze|time_utc']), str(raw['freeze|sha256']),
                                                                                           str(raw['test|start_utc']), str(raw['test|end_utc'])))
    w('')
    w('**What this measures and what it does not.** The ground truth is the geometry of a nuclear pore (a ring of radius R about 87 nm for gp210 in D1, about 53 nm for Nup96 in D2), '
      'not the emitters: which emitters are on in a frame is unknown, labelling is incomplete, and the pore centres, radii and 8-fold phases are taken from the authors\' localization '
      'tables of OTHER frames of the same field. A localizer is scored by how well its localizations of the held-out frames lie on the reference ring (a Poisson fit of the pooled '
      'radial profile: ring width w, ring yield N, uniform floor), by J_s = N / (N + floor + missed yield) and E_s = 100 - sqrt((100 - 100 J_s)^2 + (alpha w)^2), alpha = 1 per nm. '
      'It is **not** a measurement of detection accuracy per molecule, and it covers the frames and pores of ONE field of view per test set (no between-field variance). '
      'No raw DNA-origami movie with a designed spacing was publicly available, so no origami test was made; D3 (SIMFLUX) is not analysed in this report.')
    w('')
    # --- outcome
    v1, v2 = verdict_ds(an['D1']['pairs']['ISS-R2E']), verdict_ds(an['D2']['pairs']['ISS-R2E'])
    gate = all(dev[ds]['pic']['R2E']['E'] - dev[ds]['dual']['E'] > GATE_DEV_E for ds in ('D1', 'D2'))
    w('## Outcome (registered rule: the tuned ISS against Picasso R2E, tuned for E_s, on each test set)')
    w('')
    rows = []
    for ds, v in (('D1', v1), ('D2', v2)):
        pr = an[ds]['pairs']['ISS-R2E']
        rows.append([ds, ci_str(pr['dE'], pr['dE_ci']), ci_str(pr['dw'], pr['dw_ci']), ci_str(pr['dE'], pr['dE_ci_pores']).split(' ', 1)[1], ci_str(pr['dw'], pr['dw_ci_pores']).split(' ', 1)[1], v])
    w(md_table(['set', 'P1: Delta E_s, ISS - R2E (E points) [95 % block-bootstrap interval]', 'P2: Delta w at equal yield (nm) [interval]', 'P1, pore-resampling interval', 'P2, pore-resampling interval', 'reading'], rows))
    w('')
    w('**Literal reading of the registered rule: %s**' % overall_verdict(v1, v2, gate))
    w('')
    w('Negative Delta w favours the ISS; positive Delta E_s favours the ISS. The intervals are paired percentile intervals of %d resamples of blocks of %d consecutive frames (seed %d, the same blocks for every row); '
      'the pore-resampling intervals resample the pores with all their frames (seed %d, n_exp kept fixed).' % (n_boot, BLOCK, BOOT_SEED, BOOT_SEED_PORES))
    w('')
    # the registration finding
    sh = {ds: ar[ds]['shifts'] for ds in ('D1', 'D2')}
    rows = []
    for ds in ('D1', 'D2'):
        for row in ('R0', 'R2E', 'ISS', 'GREEDY'):
            a_ = np.array(sh[ds][row])
            rows.append([ds, ROW_NAMES[row], '%.0f / %.0f' % (np.median(a_[:, 0]), np.median(a_[:, 1])), '%.0f to %.0f / %.0f to %.0f' % (a_[:, 0].min(), a_[:, 0].max(), a_[:, 1].min(), a_[:, 1].max()), fmt(float(a_[:, 2].min()), 0)])
    w('### VALIDITY FINDING: the raw-frame localizations of D1 are not in the reference\'s coordinate system (found after the test result was seen)')
    w('')
    w('In D1 the localization rows of Picasso, ISS and greedy have a flat radial profile around the reference pores (ring yield N = 1 in the fit, J_s = 0), while the published table of the same frames has the ring. '
      'The cause, diagnosed afterwards from the stored localizations only (no localizer was run again): the published table of frames 10001-12000 is shifted against the raw frames by a nearly constant '
      'translation of the order of 200 to 300 nm (the authors registered / drift-corrected it; on the development frames 1-300 the shift is about zero, which is why the development tuning was not affected '
      'and the first-frame convention offset of step 3 did not see it). The pores were defined in the table\'s coordinates, so every raw-frame row is mis-registered by that shift. Translation (nm, x / y) that maps each '
      'row onto the table, by cross-correlation of 20 nm histograms per block of %d frames, with the parabolic refinement:' % REG_BLOCK)
    w('')
    w(md_table(['set', 'row', 'median shift (x / y, nm)', 'range over the blocks (x / y)', 'smallest correlation peak / s.d.'], rows))
    w('')
    w('**Reading.** The registered analysis above is, for D1, a test of whether the data are registered, not of the localizers: **D1 is void under the plan\'s section 7 ("no valid ground truth there"), by this '
      'gate that the plan did not list (registration of the frames to the reference).** D2 is registered (shifts of a few nm at most). The literal verdict above must therefore not be used. The coordinator '
      'decides whether the post-hoc analysis below, which is NOT registered, may be quoted.')
    w('')
    w('### Amendment 4 (post hoc, declared after the test result): the same rows after a per-block rigid registration to the table')
    w('')
    w('No localizer was re-run and nothing was re-tuned. Every row except the published one is translated, per block of %d frames, by the shift above (the same procedure for every row, the table of the same frames as the '
      'reference of the registration, the pores unchanged); the scoring, the intervals and the reading are the registered ones.' % REG_BLOCK)
    w('')
    rows = []
    vr = {}
    for ds in ('D1', 'D2'):
        pr = ar[ds]['pairs']['ISS-R2E']
        vr[ds] = verdict_ds(pr)
        rows.append([ds, ci_str(pr['dE'], pr['dE_ci']), ci_str(pr['dw'], pr['dw_ci']), ci_str(pr['dE'], pr['dE_ci_pores']).split(' ', 1)[1], ci_str(pr['dw'], pr['dw_ci_pores']).split(' ', 1)[1], vr[ds]])
    w(md_table(['set', 'P1: Delta E_s, ISS - R2E [interval]', 'P2: Delta w at equal yield (nm) [interval]', 'P1, pore-resampling interval', 'P2, pore-resampling interval', 'reading'], rows))
    w('')
    w('**Post-hoc reading (not registered): %s**' % overall_verdict(vr['D1'], vr['D2'], gate))
    w('')
    for ds in ('D1', 'D2'):
        a = ar[ds]
        rows = []
        for row in ROWS:
            d = a['rows'][row]
            c = d['ci']
            rows.append([ROW_NAMES[row], d['n_loc'], fmt(d['fit']['N'], 0), fmt(d['fit']['w'], 1) + ' [%s, %s]' % (fmt(c['w'][0], 1), fmt(c['w'][1], 1)), fmt(d['fit']['floor'], 0),
                         '%.3f [%.3f, %.3f]' % (d['J'], c['J'][0], c['J'][1]), ci_str(d['E'], c['E'], 1), fmt(d['mu_bias'], 1), fmt(d['yield_per_pore_kframe'], 2)])
        w('Post-hoc registered rows, %s (%d frames):' % (ds, a['n_frames']))
        w('')
        w(md_table(['row', 'localizations', 'ring N', 'w (nm) [interval]', 'floor', 'J_s [interval]', 'E_s [interval]', 'mu - R_ref (nm)', 'ring loc. per pore per 1000 frames'], rows))
        w('')
        rows = []
        for key, label in (('ISS-R2E', 'ISS tuned - Picasso R2E'), ('ISS-GREEDY', 'ISS tuned - greedy tuned (Amendment 1, descriptive)'), ('ISS0-R2E', 'ISS label-free default - R2E'), ('R2E-R0', 'R2E - R0')):
            p_ = a['pairs'][key]
            rows.append([label, ci_str(p_['dE'], p_['dE_ci']), ci_str(p_['dJ'] * 100, [100 * p_['dJ_ci'][0], 100 * p_['dJ_ci'][1]], 1), '%s (%s vs %s)' % (ci_str(p_['dw'], p_['dw_ci']) if p_['dw_ci'] else fmt(p_['dw'], 1) + ' (no interval)', fmt(p_['w_a'], 1), fmt(p_['w_b'], 1))])
        w(md_table(['paired difference', 'Delta E_s (points)', 'Delta J_s (points)', 'Delta w at equal yield (nm)'], rows))
        w('')
    # --- data and splits
    w('## Data, splits and the reference')
    w('')
    rows = []
    for sname in SETS:
        sc = Scorer(raw, sname)
        sp = DS[SETS[sname]['ds']]
        a, b = sp['dev_frames'] if SETS[sname]['role'] == 'dev' else sp['test_frames']
        rows.append([sname, SETS[sname]['ds'], '%d-%d' % (a, b), '%d x %d px at %g nm' % (sp['shape'][0] if SETS[sname]['ds'] == 'D1' else sp['window_px'], sp['shape'][1] if SETS[sname]['ds'] == 'D1' else sp['window_px'], sp['pixel_nm']),
                     sc.n_pores, fmt(sc.R, 1), '%g, %g' % (sc.conv[0], sc.conv[1])])
    w(md_table(['set', 'dataset', 'frames (1-based)', 'scored image', 'scored pores', 'R_ref (nm)', 'convention shift of the table (nm)'], rows))
    w('')
    w('D1 = %s; D2 = %s. The reference pores come from the authors\' tables restricted to the frames that are neither development nor test frames (D1: 2001-10000 and 12001-20000; D2: 1057-30000 of the '
      'same field), by the matched-filter finder of the plan (`find_pores`); the pore lists are in the raw (sha256 `%s`).' % (DS['D1']['name'], DS['D2']['name'], str(raw['pores|sha256'])))
    w('')
    w('### Validity gates (evaluated on the reference only, before any localizer ran)')
    w('')
    rows = [['V1: R_ref from the odd and even 100-frame blocks differ by at most 2 nm', 'D1 %.2f nm; D2 field 17 %.2f nm; field 18 %.2f nm' % (gates['D1']['V1_abs_diff'], gates['D2_f17']['V1_abs_diff'], gates['D2_f18']['V1_abs_diff']),
             'pass' if all(gates['V1'].values()) else 'FAIL'],
            ['V2: R_ref of the two D2 fields within 3 nm and within 5 % of 53.5 nm (literature value from memory)', 'fields 17 / 18: %.2f / %.2f nm; difference %.2f nm; vs 53.5 nm %+.1f %% / %+.1f %%' % (
                gates['D2_f17']['R_ref'], gates['D2_f18']['R_ref'], gates['V2']['diff'], 100 * gates['V2']['lit_dev'], 100 * gates['V2']['lit_test']), 'pass' if gates['V2']['pass'] else 'FAIL'],
            ['V3: at least 50 % of the D2 pores have a significant 8th harmonic (Rayleigh p < 0.05 in place of the permutation test)', 'fields 17 / 18: %.0f %% / %.0f %%' % (100 * gates['V3']['D2_f17'], 100 * gates['V3']['D2_f18']),
             'pass' if all(v >= 0.5 for v in gates['V3'].values()) else 'FAIL'],
            ['V4: at least 100 scored pores per set', ', '.join('%s %d' % (k[3:], v['n']) for k, v in sorted(gates.items()) if k.startswith('V4_')), 'pass' if all(v['ok'] for k, v in gates.items() if k.startswith('V4_')) else 'FAIL']]
    w(md_table(['gate', 'value', 'result'], rows))
    w('')
    w('R_ref (median circle-fit radius of the accepted pores): D1 %.2f nm over %d accepted pores; D2 field 17 %.2f nm (%d pores), field 18 %.2f nm (%d pores). The radial-profile fit of the reference '
      'localizations (nearest-pore assignment, uniform floor) gives a ring radius mu of %.1f nm (D1) and %.1f nm (D2 test), i.e. mu - R_ref = %+.1f / %+.1f nm: the difference between the two ways of '
      'measuring a radius is a property of the reference and of the model, common to all rows; the radius bias of a row (mu - R_ref in the tables) contains it and is a secondary outcome only.' % (
          gates['D1']['R_ref'], gates['D1']['n_pores'], gates['D2_f17']['R_ref'], gates['D2_f17']['n_pores'], gates['D2_f18']['R_ref'], gates['D2_f18']['n_pores'],
          an['D1']['ref_fit']['mu'], an['D2']['ref_fit']['mu'], an['D1']['ref_fit']['mu'] - an['D1']['R_ref'], an['D2']['ref_fit']['mu'] - an['D2']['R_ref']))
    w('')
    w('### Coordinate convention, label-free calibration and Picasso scale factors (development data, no labels)')
    w('')
    rows = []
    for sname in SETS:
        sig = jload(raw['cal|%s|sigma0' % sname])
        law = raw['cal|%s|law' % sname]
        rows.append([sname, '%.3f px (%.1f nm)' % (sig['sigma0'], sig['sigma0'] * DS[SETS[sname]['ds']]['pixel_nm']), '%d of %d R0 localizations%s' % (sig['n_sel'], sig['n_all'], ' (FALLBACK)' if sig['fallback'] else ''),
                     'a %.3f, b0 %.1f ADU^2, x_ref %.1f ADU, scale %.3f' % (law[0], law[1], law[2], law[4])])
    w(md_table(['set', 'sigma0 (label-free)', 'its sample', 'noise law var = a (x - x_ref) + b0'], rows))
    w('')
    for ds in ('D1', 'D2'):
        co = raw['cal|%s|convention_offset' % ds]
        fa = dev[ds]['pic']['fac']
        w('* %s: median (Picasso R0 - table) offset on %d pairs of the development frames: (%.1f, %.1f) nm, applied to the table as (%g, %g) nm (the nearest multiple of half a pixel); Picasso scale factors from '
          'the development frames: net gradient x %.3f, photons x %.3f, widths x %.3f (N1 reference: median net gradient %.0f and photons %.0f of Picasso at net gradient %g on MT0.N1.LD frames 0-999).' % (
              ds, co[2], co[0], co[1], convention(raw, ds)[0], convention(raw, ds)[1], fa['f_ng'], fa['f_ph'], fa['f_sig'], float(raw['cal|N1|ng_median']), float(raw['cal|N1|ph_median']), N1_LOW_NG))
    w('')
    # --- development
    w('## Development tuning (E_s on the development frames; same objective for every method)')
    w('')
    rows = []
    for ds in ('D1', 'D2'):
        d = dev[ds]
        gp = d['pic']
        rows.append([ds, d['dual']['n_frames'], fmt(gp['R0']['E'], 2), fmt(gp['R1E']['E'], 2), fmt(gp['R2E']['E'], 2), gp['R2E']['n_configs'], fmt(d['dual']['E'], 2),
                     '%s, kappa %g, width %g x sigma0/1.4' % (d['dual']['cfg'], d['dual']['kappa'], d['dual']['smax_rel']), fmt(d['dual']['default_E'], 2), fmt(d['greedy']['E'], 2),
                     '%s, kappa %g, width %g' % (d['greedy']['cfg'], d['greedy']['kappa'], d['greedy']['smax_rel'])])
    w(md_table(['set', 'dev frames', 'R0', 'R1E', 'R2E', 'Picasso configs scored (R2)', 'ISS tuned', 'ISS point', 'ISS label-free default', 'greedy tuned', 'greedy point'], rows))
    w('')
    w('Development E_s of the frozen configurations. The ISS search is 168 operating points (8 configurations x 7 kappas x 3 width bounds); Picasso R2 scores the numbers of configurations above '
      '(the simulated study scored 12,000 to 15,000 per set). R2E: %s.' % '; '.join('%s %s' % (ds, dev[ds]['pic']['R2E']['struct']) + ' with ' + ', '.join('%s=%s' % (k, 'off' if v is None else fmt(v, 1)) for k, v in sorted(dev[ds]['pic']['R2E']['filters'].items())) for ds in ('D1', 'D2')))
    w('')
    w('Development gate (plan section 7 iii: Picasso R2E exceeds the tuned ISS by more than %g E points on both sets): D1 %+.2f, D2 %+.2f points (R2E minus ISS): %s.' % (
        GATE_DEV_E, dev['D1']['pic']['R2E']['E'] - dev['D1']['dual']['E'], dev['D2']['pic']['R2E']['E'] - dev['D2']['dual']['E'], 'MET' if gate else 'not met'))
    w('')
    w('**Edge of the declared grid (stated after the freeze, nothing changed).** The coordinator noticed, from the development log, that the best significance threshold of many configurations is kappa = %g, the TOP of the declared '
      'grid %s, which suggests that a larger kappa might have scored better, and that the tuned ISS is below Picasso R2E and even below R0 on the development frames. The freeze (%s) was already written when this '
      'was noticed, so the grid was NOT extended and the registered configurations stand. The facts, from the development rows stored in the raw:' % (KAPPAS[-1], list(KAPPAS), str(raw['freeze|time_utc'])))
    w('')
    rows = []
    for ds in ('D1', 'D2'):
        for order in ('dual', 'greedy'):
            d = dev[ds][order]
            top = {c: max(r, key=lambda x: x[2]) for c, r in d['rows'].items()}
            at_edge = sum(1 for t in top.values() if t[0] == KAPPAS[-1])
            rows.append([ds, order, '%s, kappa %g' % (d['cfg'], d['kappa']), '%d of %d' % (at_edge, len(top)), ', '.join('%s %g' % (c, t[0]) for c, t in sorted(top.items())),
                         fmt(d['E'], 2), fmt(dev[ds]['pic']['R0']['E'], 2), fmt(dev[ds]['pic']['R2E']['E'], 2)])
    w(md_table(['set', 'order', 'frozen point', 'configurations with their best kappa at the top of the grid', 'best kappa of every configuration', 'dev E_s', 'Picasso R0 dev E_s', 'Picasso R2E dev E_s'], rows))
    w('')
    w('So the optimum is at the edge for D1 (both orders) but NOT for the frozen D2 points (kappa 4.5 and 5), although some D2 configurations also peak at kappa 10. A development optimum at the edge means that the ISS rows are '
      'probably tuned too little on the threshold relative to Picasso, whose coordinate ascent over the net gradient is not capped in the same way; the test result for the ISS may therefore be pessimistic. '
      'Descriptively (test frames, after the registered run; no label-free diagnostic such as the chi2 of the ISS fit was stored): the numbers of localizations per frame of the tuned ISS, of Picasso R2E and of the '
      'authors\' table are %s. A label-free noise law that is too small for these cameras (EM gain excess noise in D1; an unknown gain in D2) would make the ISS admit noise entrants and push the best kappa upward; '
      'the figures above are consistent with that but do not prove it.' % '; '.join('%s: %.2f / %.2f / %.2f' % (ds, an[ds]['rows']['ISS']['n_loc'] / an[ds]['n_frames'], an[ds]['rows']['R2E']['n_loc'] / an[ds]['n_frames'],
                                                                                          an[ds]['rows']['PUB']['n_loc'] / an[ds]['n_frames']) for ds in ('D1', 'D2')))
    w('')
    # --- test tables
    for ds in ('D1', 'D2'):
        a = an[ds]
        w('## Test set %s as run, NOT registered to the table (%d frames, %d scored pores, expected ring yield n_exp = %.0f)' % (ds, a['n_frames'], a['n_pores'], a['n_exp']))
        w('')
        rows = []
        for row in ROWS:
            d = a['rows'][row]
            c = d['ci']
            rows.append([ROW_NAMES[row], d['n_loc'], fmt(d['fit']['N'], 0), fmt(d['fit']['w'], 1) + ' [%s, %s]' % (fmt(c['w'][0], 1), fmt(c['w'][1], 1)), fmt(d['fit']['floor'], 0),
                         '%.3f [%.3f, %.3f]' % (d['J'], c['J'][0], c['J'][1]), ci_str(d['E'], c['E'], 1), fmt(d['mu_bias'], 1), fmt(d['yield_per_pore_kframe'], 2), fmt(100 * d['frac_pores_8'], 0),
                         fmt(100 * d['dup_frac'], 1), ('%.4f' % d['secs_per_frame']) if np.isfinite(d['secs_per_frame']) else 'n/a'])
        w(md_table(['row', 'localizations', 'ring N', 'w (nm) [interval]', 'floor', 'J_s [interval]', 'E_s [interval]', 'mu - R_ref (nm)', 'ring loc. per pore per 1000 frames', 'pores with >= 8 ring loc. (%)',
                    'duplicates within sigma0 (%)', 's per frame'], rows))
        w('')
        pr = a['pairs']
        rows = []
        for key, label in (('ISS-R2E', 'ISS tuned - Picasso R2E (registered P1 / P2)'), ('ISS-GREEDY', 'ISS tuned - greedy tuned (Amendment 1, descriptive)'), ('ISS0-R2E', 'ISS label-free default - R2E'), ('R2E-R0', 'R2E - R0')):
            p_ = pr[key]
            rows.append([label, ci_str(p_['dE'], p_['dE_ci']), ci_str(p_['dJ'] * 100, [100 * p_['dJ_ci'][0], 100 * p_['dJ_ci'][1]], 1), '%s (%s vs %s)' % (ci_str(p_['dw'], p_['dw_ci']) if p_['dw_ci'] else fmt(p_['dw'], 1) + ' (no interval)', fmt(p_['w_a'], 1), fmt(p_['w_b'], 1))])
        w(md_table(['paired difference', 'Delta E_s (points)', 'Delta J_s (points)', 'Delta w at equal yield (nm)'], rows))
        w('')
        if DS[ds]['eight_fold']:
            rows = [[ROW_NAMES[r], '%.4f [%.4f, %.4f]' % (a['rows'][r]['c8'], a['rows'][r]['ci']['c8'][0], a['rows'][r]['ci']['c8'][1]), '%.4f' % a['rows'][r]['c8_null']] for r in ROWS]
            w('Corner contrast C8 (mean cos(8 (theta - phi_j)) of the ring localizations, phases phi_j from the reference; null = phases permuted across pores, mean of %d permutations; not thinned to equal yield):' % N_THIN)
            w('')
            w(md_table(['row', 'C8 [interval]', 'null'], rows))
            w('')
        iss = a['rows']['ISS']
        if 'mean_entries' in iss:
            w('ISS paths: mean %.1f entries per frame, %.1f %% of the frames reached the cap of 400 entries; greedy: %.1f entries, %.1f %% at the cap.' % (
                iss['mean_entries'], 100 * iss['frac_at_cap'], a['rows']['GREEDY']['mean_entries'], 100 * a['rows']['GREEDY']['frac_at_cap']))
            w('')
    w('## Not done, deviations, weaknesses')
    w('')
    w('* ThunderSTORM was not run: the 1.3 jar (HohlbeinLab update site; the official update-site path returns 404) and ImageJ 1.54p were obtained, but in the Xvfb session the macro commands "Run analysis" and "Export results" '
      'returned without producing a results file after about one hour of attempts, and the cause was not found; the authors\' published tables stand in as context rows ("as published"), as the plan said.')
    w('* The permutation test of gate V3 was replaced by the Rayleigh p-value of the harmonic (the same null, analytic).')
    w('* The pore lists were frozen in the table\'s coordinates; the scored pores of a set are those whose disc lies inside the scored region AFTER the convention shift (D2: half a pixel), by a fixed rule, so that a few pores at the window border are dropped.')
    w('* One field of view per test set; intervals cover frames (and, in the sensitivity, pores) of that field only. The reference localizer (the authors\' software) and the pore model are shared by all rows; the "as published" rows are not comparators.')
    w('* D3 (SIMFLUX nanorulers) is outside the inference and not analysed here.')
    w('')
    return '\n'.join(L) + '\n'


def write_report(raw_path=RAW, md_path=MD, n_boot=N_BOOT, log=print):
    raw = load_raw(raw_path)
    md = build_md(raw, n_boot, log)
    with open(md_path, 'w') as f:
        f.write(md)
    log('wrote %s' % md_path)
    return md


def load_raw(path=None):
    path = path or RAW
    if not os.path.exists(path):
        return {}
    with np.load(path, allow_pickle=False) as z:
        return {k: z[k] for k in z.files}


def save_raw(raw, path=None):
    path = path or RAW
    os.makedirs(os.path.dirname(path), exist_ok=True)
    C.save_npz_deterministic(path, raw)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('--part', choices=('pores', 'cal', 'dev', 'test'))
    ap.add_argument('--from-raw', action='store_true')
    a = ap.parse_args(argv)
    raw = load_raw()
    if a.from_raw:
        write_report()
        return
    if a.part == 'dev':
        step_dev(raw)
    if a.part == 'test':
        step_test_all(raw)
    if a.part == 'cal':
        step_cal(raw)
        save_raw(raw)
    if a.part == 'pores':
        gates = step_pores(raw)
        gates = finish_gates(raw, gates)
        raw['pores|sha256'] = np.array(pore_list_hash(raw))
        save_raw(raw)
        print(json.dumps(gates, indent=1, sort_keys=True))
        print('pore list sha256', str(raw['pores|sha256']))


if __name__ == '__main__':
    main()
