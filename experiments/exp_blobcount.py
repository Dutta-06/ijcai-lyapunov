"""
Counting overlapping nuclei in fluorescence microscopy patches: the sliding inverse scale space against tuned classical detectors (a spike with a kill rule)
=====================================================================================================================================================

Question, data, split, methods, metrics and the GO / KILL rule are fixed in `results/ml/blobcount_plan.md` (committed before any patch was scored). In short: 128 x 128 patches of
the fluorescence images of the 2018 Data Science Bowl `stage1_train` (public Kaggle copy), split by image (60 development, 150 test, a seeded rule), ground-truth count = the masks
whose centroid lies in the patch; the target is the count itself. Methods: `ours` (the repository's sliding path, dual order, learned width per unit, BIC stop), `greedy` (the same with
the greedy order of entry: the control that matters), `log`, `dog`, `plm` (peak_local_max), `ws` (Otsu / adaptive threshold + distance-transform watershed). Every method has at most 24
configurations, chosen on the development patches by one objective, frozen, and run once on the test patches.

Stages (CLI):  bank   build the patch bank from the dataset archive (run inside a private notebook; no image leaves it except the derived 128 x 128 patches)
               dev    run every configuration of the chosen methods on the development patches (counts and matches only)
               freeze choose the configuration of every method from the development raw and write the freeze file
               test   run the frozen configurations once on the test patches
               report regenerate the markdown from a raw file (`--from-raw`), byte-identical
"""
import argparse
import hashlib
import io
import json
import os
import sys
import time

import numpy as np
from scipy import ndimage as ndi
from scipy.optimize import linear_sum_assignment

HERE = os.path.abspath(os.path.dirname(__file__))
REPO = os.environ.get('BLOBCOUNT_REPO') or os.path.abspath(os.path.join(HERE, '..'))
if REPO not in sys.path:
    sys.path.insert(0, REPO)

# ----------------------------------------------------------------------------------------------------------------------------------
# constants fixed by the plan
# ----------------------------------------------------------------------------------------------------------------------------------
SEED_SPLIT = 1847            # numpy default_rng seed of the image split
SEED_PATCH = 20181          # seed of the per-image tile choice
N_DEV, N_TEST = 60, 150
PATCH = 128
MIN_AREA = 10                # masks below this many pixels are debris
CROWD = 15                   # crowded stratum: at least this many ground-truth nuclei (Amendment 2 of the plan; 25 in the first draft)
CROWD_HI = 25                # the very crowded sub-stratum, reported separately
PER_IMAGE, BOOST = 2, 4      # base tiles per image, crowded-boost tiles per image
MED_MAX = 100                # fluorescence: grey and median grey level below this
R_MATCH = 6.0                # pixels: match radius of the detection F1
N_BOOT = 10000
MAX_CONFIGS = 24
M_REL = 0.05                 # C1: relative margin of the crowded gain
M_GREEDY = 0.5               # C2: nuclei per patch
M_NONINF = 0.10              # C3: relative non-inferiority margin
METHODS = ('ours', 'greedy', 'log', 'dog', 'plm', 'ws')
METHODS_NC = ('ours_nc', 'greedy_nc')           # POST HOC variants: the same grids without the early stop of the path (snr_min None)
SLIDING = ('ours', 'greedy') + METHODS_NC
CLASSICAL = ('log', 'dog', 'plm', 'ws')


# ----------------------------------------------------------------------------------------------------------------------------------
# data: eligibility, split, tiles, ground truth
# ----------------------------------------------------------------------------------------------------------------------------------
def is_fluorescence(rgb):
    """Mask-free eligibility: the colour channels are identical at every pixel and the median grey level is below MED_MAX."""
    rgb = np.asarray(rgb)
    if rgb.ndim == 2:
        return bool(np.median(rgb) < MED_MAX)
    c = rgb[..., :3]
    return bool(np.all(c[..., 0] == c[..., 1]) and np.all(c[..., 1] == c[..., 2]) and np.median(c[..., 0]) < MED_MAX)


def split_ids(ids):
    """(development ids, test ids): the sorted eligible ids permuted by default_rng(SEED_SPLIT); the first N_DEV, the next N_TEST."""
    ids = sorted(ids)
    perm = np.random.default_rng(SEED_SPLIT).permutation(len(ids))
    order = [ids[i] for i in perm]
    return order[:N_DEV], order[N_DEV:N_DEV + N_TEST]


def mask_stats(masks):
    """Centroids (x = mean column, y = mean row; pixel centres at integers) and areas of the masks with at least MIN_AREA pixels."""
    xy, area = [], []
    for m in masks:
        a = int(m.sum())
        if a < MIN_AREA:
            continue
        r, c = np.nonzero(m)
        xy.append((c.mean(), r.mean()))
        area.append(a)
    return np.array(xy, float).reshape(-1, 2), np.array(area, float)


def tile_grid(h, w):
    """Origins (row, column) of the non-overlapping PATCH x PATCH tiles of an h x w image (the remainder at the border is dropped)."""
    return [(r * PATCH, c * PATCH) for r in range(h // PATCH) for c in range(w // PATCH)]


def in_patch(xy, r0, c0):
    """Boolean mask of the points of xy (x, y) in image coordinates that lie in the tile at (r0, c0): patch coordinates in [-0.5, PATCH - 0.5)."""
    x, y = xy[:, 0] - c0, xy[:, 1] - r0
    return (x >= -0.5) & (x < PATCH - 0.5) & (y >= -0.5) & (y < PATCH - 0.5)


def choose_tiles(image_id, counts):
    """(base, boost): indices into the tile list. base = PER_IMAGE tiles by a per-image seeded permutation (all if fewer); boost = up to BOOST further tiles, in the same permutation
    order, with at least CROWD ground-truth nuclei."""
    seed = int(hashlib.sha256(image_id.encode()).hexdigest()[:8], 16)
    perm = np.random.default_rng([SEED_PATCH, seed]).permutation(len(counts))
    base = [int(i) for i in perm[:PER_IMAGE]]
    boost = [int(i) for i in perm[PER_IMAGE:] if counts[i] >= CROWD][:BOOST]
    return base, boost


def touching_ids(label):
    """Set of mask labels (> 0) that touch another mask in the 4-neighbourhood."""
    out = set()
    for a, b in ((label[:, :-1], label[:, 1:]), (label[:-1, :], label[1:, :])):
        m = (a > 0) & (b > 0) & (a != b)
        out.update(np.unique(a[m]).tolist())
        out.update(np.unique(b[m]).tolist())
    return out


def build_bank(zip_path, log=print):
    """Cut the patch bank from the dataset archive. Returns (arrays, info)."""
    import zipfile
    from PIL import Image
    z = zipfile.ZipFile(zip_path)
    names = z.namelist()
    ids = sorted({n.split('/')[2] for n in names if '/stage1_train/' in n and n.endswith('.png') and '/images/' in n})
    prefix = 'data-science-bowl-2018/stage1_train/'
    elig, imgs = [], {}
    for i in ids:
        im = np.array(Image.open(io.BytesIO(z.read('%s%s/images/%s.png' % (prefix, i, i)))))
        if is_fluorescence(im):
            elig.append(i)
            imgs[i] = im[..., 0] if im.ndim == 3 else im
    dev, test = split_ids(elig)
    log('images %d eligible %d dev %d test %d' % (len(ids), len(elig), len(dev), len(test)))
    rows, patches, gts, areas, pool = [], [], [], [], []
    for sp, group in ((0, dev), (1, test)):
        for i in group:
            g = imgs[i]
            h, w = g.shape
            mnames = sorted(n for n in names if n.startswith('%s%s/masks/' % (prefix, i)) and n.endswith('.png'))
            masks = [np.array(Image.open(io.BytesIO(z.read(n)))) > 0 for n in mnames]
            keep = [m for m in masks if m.sum() >= MIN_AREA]
            xy, area = mask_stats(keep)
            label = np.zeros(g.shape, np.int32)
            for k, m in enumerate(keep):
                label[m] = k + 1
            touch = touching_ids(label)
            tiles = tile_grid(h, w)
            counts = [int(in_patch(xy, r0, c0).sum()) for r0, c0 in tiles]
            base, boost = choose_tiles(i, counts)
            pool.append((sp, i, counts))
            for flag, sel in ((0, base), (1, boost)):
                for t in sel:
                    r0, c0 = tiles[t]
                    ins = in_patch(xy, r0, c0)
                    patches.append(g[r0:r0 + PATCH, c0:c0 + PATCH].astype(np.uint8))
                    gts.append((xy[ins] - np.array([c0, r0])).astype(np.float32))
                    areas.append(area[ins].astype(np.float32))
                    nt = int(sum(1 for k in np.flatnonzero(ins) if (k + 1) in touch))
                    rows.append((sp, i, r0, c0, flag, int(ins.sum()), nt))
    P = len(rows)
    mx = max(len(a) for a in gts)
    gt_xy = np.full((P, mx, 2), np.nan, np.float32)
    gt_area = np.full((P, mx), np.nan, np.float32)
    for k, (a, b) in enumerate(zip(gts, areas)):
        gt_xy[k, :len(a)] = a
        gt_area[k, :len(a)] = b
    uid = sorted({r[1] for r in rows})
    arrays = dict(patches=np.stack(patches), split=np.array([r[0] for r in rows], np.uint8), img=np.array([uid.index(r[1]) for r in rows], np.int32),
                  img_ids=np.array(uid), r0=np.array([r[2] for r in rows], np.int16), c0=np.array([r[3] for r in rows], np.int16),
                  boost=np.array([r[4] for r in rows], np.uint8), n_gt=np.array([r[5] for r in rows], np.int16), n_touch=np.array([r[6] for r in rows], np.int16),
                  gt_xy=gt_xy, gt_area=gt_area)
    info = dict(n_images=len(ids), n_eligible=len(elig), dev_ids=dev, test_ids=test, n_patches=P, pool_counts=[[sp, i, c] for sp, i, c in pool])
    return arrays, info


def load_bank(path):
    with np.load(path, allow_pickle=False) as z:
        return {k: z[k] for k in z.files}


def gt_of(bank, k):
    """Ground-truth centroids (n, 2) of patch k (patch coordinates)."""
    n = int(bank['n_gt'][k])
    return bank['gt_xy'][k, :n].astype(float)


def strata(bank, split):
    """Index arrays (all = base sample, crowded = any sampled patch with >= CROWD nuclei) of the patches of one split (0 development, 1 test)."""
    sel = bank['split'] == split
    base = np.flatnonzero(sel & (bank['boost'] == 0))
    crowd = np.flatnonzero(sel & (bank['n_gt'] >= CROWD))
    return base, crowd


# ----------------------------------------------------------------------------------------------------------------------------------
# the shared label-free front end
# ----------------------------------------------------------------------------------------------------------------------------------
def noise_model(g, nbin=8):
    """(a, b): pixel variance ~ a * mean + b, estimated from the patch itself. The squared horizontal first differences (halved) are binned by the local mean; the median of each bin over
    0.4549 (the median of chi-square with one degree of freedom) estimates the variance there, robust to the nuclei's edges; a weighted line is fitted (a >= 0, b >= 0.25)."""
    g = np.asarray(g, float)
    d2 = (g[:, 1:] - g[:, :-1]) ** 2 / 2.0
    m = 0.5 * (g[:, 1:] + g[:, :-1])
    edges = np.quantile(m, np.linspace(0, 1, nbin + 1))
    xs, vs, ws = [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        s = (m >= lo) & (m <= hi) if hi == edges[-1] else (m >= lo) & (m < hi)
        if s.sum() < 50:
            continue
        xs.append(float(np.median(m[s])))
        vs.append(float(np.median(d2[s])) / 0.4549)
        ws.append(float(s.sum()))
    xs, vs, ws = np.array(xs), np.array(vs), np.array(ws)
    if len(xs) < 2 or np.ptp(xs) < 1e-9:
        return 0.0, max(float(np.average(vs, weights=ws)) if len(vs) else 1.0, 0.25)
    W = np.sqrt(ws)
    A = np.c_[xs, np.ones_like(xs)] * W[:, None]
    a, b = np.linalg.lstsq(A, vs * W, rcond=None)[0]
    if a < 0:
        a, b = 0.0, float(np.average(vs, weights=ws))
    return float(a), float(max(b, 0.25))


def scale_estimate(p, sigmas=None, top=12, default=5.0):
    """The patch scale s_hat (pixels): the median width of the `top` strongest local maxima of the scale-normalised Laplacian of Gaussian -s^2 LoG_s p over space and scale."""
    sigmas = np.geomspace(2.0, 14.0, 9) if sigmas is None else sigmas
    R = np.stack([-(s * s) * ndi.gaussian_laplace(np.asarray(p, float), s, mode='reflect') for s in sigmas])
    pk = (R == ndi.maximum_filter(R, size=(3, 3, 3), mode='nearest')) & (R > 0)
    if not pk.any():
        return float(default)
    vals = R[pk]
    si = np.nonzero(pk)[0]
    o = np.argsort(-vals)[:top]
    return float(np.median(sigmas[si[o]]))


def prepare(patch):
    """The shared front end of one patch: background-subtracted image y, inverse-sd map sw, variance map var, patch scale s_hat."""
    g = np.asarray(patch, float)
    med = float(np.median(g))
    a, b = noise_model(g)
    var = np.maximum(a * ndi.gaussian_filter(g, 1.5) + b, 0.25)
    y = g - med
    return dict(g=g, y=y, var=var, sw=1.0 / np.sqrt(var), s_hat=scale_estimate(y), a=a, b=b, med=med)


# ----------------------------------------------------------------------------------------------------------------------------------
# ours and greedy: the repository's sliding path
# ----------------------------------------------------------------------------------------------------------------------------------
PATH_GRID = [
    dict(s0=1.3, lo=0.7, hi=3.0, hp=0.0, rex=0.0),
    dict(s0=1.3, lo=0.7, hi=3.0, hp=0.0, rex=0.8),
    dict(s0=1.8, lo=1.0, hi=4.0, hp=0.0, rex=0.0),
    dict(s0=1.8, lo=1.0, hi=4.0, hp=0.0, rex=0.8),
]
PATH_FIXED = dict(h=1.0, margin=1.0, box=1.0, fit=3.0, max_entries=120, snr_min=2.0, extra=6, rtrunc=4.0)
STOP_GRID = [dict(lam=l, clean=c, cal=1) for l, c in ((64.0, 0.0), (256.0, 0.0), (1024.0, 0.0), (64.0, 1.0), (256.0, 1.0), (1024.0, 1.0))]


def configs_sliding(nocap=False):
    fixed = dict(PATH_FIXED, snr_min=None) if nocap else PATH_FIXED
    return [dict(path=dict(p, **fixed), stop=s) for p in PATH_GRID for s in STOP_GRID]


def order_of(method):
    return 'dual' if method in ('ours', 'ours_nc') else 'greedy'


def sliding_run(prep, path, order):
    """The sliding path of one patch. Returns (psf, hist, n_entries)."""
    from src.core import imaging as IM
    from src.core import imaging_fast as F
    sh = prep['s_hat']
    P = IM.PixelPSF(PATCH, PATCH, path['s0'] * sh, hp=(path['hp'] * sh or None))
    L = F.FastLattice(P, h=path['h'], margin=path['margin'])
    y = P.highpass(prep['y'])
    R = F.sliding_path_fast(P, y, prep['sw'], L, max_entries=path['max_entries'], tol=1e-6, passes=1, rtrunc=path['rtrunc'], r_ex=path['rex'] * sh,
                            fit_radius=path['fit'] * sh, order=order, sig_bounds=(path['lo'] * sh, path['hi'] * sh), box=path['box'] * sh,
                            snr_min=path['snr_min'], snr_extra=path['extra'])
    return P, R


def sliding_stop(prep, P, R, stop):
    """The detections (x, y, width, amplitude) of the block chosen by the BIC stop with penalty `lam` per parameter (4 parameters per unit: centre, width, amplitude), optionally
    followed by the clean-up of units whose own significance is below `clean` * sqrt(lam ln N). The noise calibration (`cal` = 1): the penalty is multiplied by the weighted mean square
    residual of the path's last block (the residual variance of the over-complete fit, in units of the pixel-noise model; the nuclei's texture is not pixel noise, so the first-difference
    noise model alone underestimates it). Returns (array (n, 4), index of the block, hit_end)."""
    from src.core import imaging as IM
    N = PATCH * PATCH
    hist = R['hist']
    lam = stop['lam'] * (max(hist[-1]['chi2'], 1.0) if stop.get('cal', 1) else 1.0)
    i, _ = IM.stop_bic(hist, N, n_par=lam)
    h = hist[i]
    c, w = h['c'], h['w']
    if stop['clean'] and len(w):
        keep = IM.unit_significance(P, prep['sw'], c, w) >= stop['clean'] * np.sqrt(lam * np.log(N))
        c, w = c[keep], w[keep]
    return np.c_[c[:, :3], w].reshape(-1, 4), int(i), bool(i == len(hist) - 1)


def run_sliding_group(prep, cfgs, order):
    """All stops of the configurations that share one path. cfgs: list of (config index, config). Returns {index: (dets, hit_end)}."""
    P, R = sliding_run(prep, cfgs[0][1]['path'], order)
    out = {}
    for k, c in cfgs:
        d, _, hit = sliding_stop(prep, P, R, c['stop'])
        out[k] = (d, hit)
    return out


# ----------------------------------------------------------------------------------------------------------------------------------
# classical baselines
# ----------------------------------------------------------------------------------------------------------------------------------
def _norm01(img):
    """The image scaled so that its 99.9th percentile is 1 (negative values clipped at 0)."""
    img = np.maximum(np.asarray(img, float), 0.0)
    s = np.percentile(img, 99.9)
    return img / s if s > 0 else img


SCALES = {'rel_narrow': (1.0, 1.8, True), 'rel_wide': (0.7, 2.6, True), 'abs': (3.0, 12.0, False)}


BLOB_THR = {'log': (0.04, 0.08, 0.15, 0.3), 'dog': (0.15, 0.3, 0.6, 1.2)}


def configs_blob(kind='log'):
    return [dict(scale=sc, thr=t, overlap=o) for sc in ('abs', 'rel_wide') for t in BLOB_THR[kind] for o in (0.15, 0.3, 0.6)]


def run_blob(prep, cfg, kind):
    from skimage.feature import blob_dog, blob_log
    lo, hi, rel = SCALES[cfg['scale']]
    lo, hi = (lo * prep['s_hat'], hi * prep['s_hat']) if rel else (lo, hi)
    hi = max(hi, lo * 1.2)
    img = _norm01(prep['y'])
    if kind == 'log':
        b = blob_log(img, min_sigma=lo, max_sigma=hi, num_sigma=6, threshold=cfg['thr'], overlap=cfg['overlap'])
    else:
        b = blob_dog(img, min_sigma=lo, max_sigma=hi, sigma_ratio=1.4, threshold=cfg['thr'] / 4.0, overlap=cfg['overlap'])
    return np.c_[b[:, 1], b[:, 0]].reshape(-1, 2) if len(b) else np.zeros((0, 2))


def configs_plm():
    rel = [dict(sm=s, md=m, thr=0.12, rel=True) for s in (0.5, 0.8, 1.2) for m in (0.8, 1.3, 2.0)]
    ab = [dict(sm=s, md=m, thr=0.12, rel=False) for s in (2.0, 3.0, 4.0, 5.0) for m in (2, 3, 5)]
    th = [dict(sm=3.0, md=3, thr=0.05, rel=False), dict(sm=3.0, md=3, thr=0.25, rel=False), dict(sm=4.0, md=3, thr=0.25, rel=False)]
    return rel + ab + th


def run_plm(prep, cfg):
    from skimage.feature import peak_local_max
    sm, md = (cfg['sm'] * prep['s_hat'], cfg['md'] * prep['s_hat']) if cfg['rel'] else (cfg['sm'], cfg['md'])
    img = _norm01(ndi.gaussian_filter(prep['y'], sm))
    pk = peak_local_max(img, min_distance=max(int(round(md)), 1), threshold_rel=cfg['thr'], exclude_border=False)
    return np.c_[pk[:, 1], pk[:, 0]].reshape(-1, 2).astype(float) if len(pk) else np.zeros((0, 2))


def configs_ws():
    out = []
    for sm in (0.3, 0.6):
        for th in ('otsu', 'otsu0.7', 'local'):
            for q in (0.7, 1.2):
                out.append(dict(fam='peak', sm=sm, thr=th, q=q))
    for sm in (0.3, 0.6):
        for th in ('otsu0.7', 'local'):
            for q in (0.5, 1.0, 2.0):
                out.append(dict(fam='hmax', sm=sm, thr=th, q=q))
    return out


def run_ws(prep, cfg):
    from skimage.feature import peak_local_max
    from skimage.filters import threshold_local, threshold_otsu
    from skimage.morphology import h_maxima
    from skimage.segmentation import watershed
    s = ndi.gaussian_filter(prep['y'], cfg['sm'] * prep['s_hat'])
    if cfg['thr'] == 'local':
        mask = s > threshold_local(s, 51, offset=0) + 0.25 * max(float(s.std()), 1e-6)
        mask &= s > 0.3 * threshold_otsu(s)
    else:
        mask = s > threshold_otsu(s) * (0.7 if cfg['thr'] == 'otsu0.7' else 1.0)
    mask = ndi.binary_opening(mask, iterations=1)
    if not mask.any():
        return np.zeros((0, 2))
    edt = ndi.gaussian_filter(ndi.distance_transform_edt(mask), 1.0)
    if cfg['fam'] == 'peak':
        pk = peak_local_max(edt, min_distance=max(int(round(cfg['q'] * prep['s_hat'])), 1), labels=ndi.label(mask)[0], exclude_border=False)
        markers = np.zeros(mask.shape, np.int32)
        markers[pk[:, 0], pk[:, 1]] = np.arange(1, len(pk) + 1)
    else:
        markers = ndi.label(h_maxima(edt, cfg['q']))[0]
    lab = watershed(-edt, markers, mask=mask)
    ids = [k for k in range(1, int(lab.max()) + 1) if (lab == k).sum() >= 12]
    if not ids:
        return np.zeros((0, 2))
    cm = np.array(ndi.center_of_mass(mask, lab, ids))
    return np.c_[cm[:, 1], cm[:, 0]]


def configs_of(method):
    if method in ('ours', 'greedy'):
        return configs_sliding()
    if method in METHODS_NC:
        return configs_sliding(nocap=True)
    return {'log': lambda: configs_blob('log'), 'dog': lambda: configs_blob('dog'), 'plm': configs_plm, 'ws': configs_ws}[method]()


def run_group(method, prep, group):
    """Run one group of configurations of a method on a prepared patch. group: list of (config index, config). Returns {index: (dets (n, >=2), hit_end)}."""
    if method in SLIDING:
        return run_sliding_group(prep, group, order_of(method))
    fn = {'log': lambda c: run_blob(prep, c, 'log'), 'dog': lambda c: run_blob(prep, c, 'dog'), 'plm': lambda c: run_plm(prep, c), 'ws': lambda c: run_ws(prep, c)}[method]
    return {k: (fn(c), False) for k, c in group}


def groups_of(method, cfgs, only=None):
    """Partition the configurations into groups that share an expensive run (the sliding path), as lists of (index, config)."""
    sel = range(len(cfgs)) if only is None else only
    if method in SLIDING:
        g = {}
        for k in sel:
            g.setdefault(json.dumps(cfgs[k]['path'], sort_keys=True), []).append((k, cfgs[k]))
        return list(g.values())
    return [[(k, cfgs[k])] for k in sel]


# ----------------------------------------------------------------------------------------------------------------------------------
# scoring
# ----------------------------------------------------------------------------------------------------------------------------------
def inside(xy):
    """Rows of xy (patch coordinates) whose centre lies in the patch."""
    xy = np.asarray(xy, float).reshape(-1, xy.shape[-1] if np.ndim(xy) > 1 else 2)
    if len(xy) == 0:
        return xy
    return xy[(xy[:, 0] >= -0.5) & (xy[:, 0] < PATCH - 0.5) & (xy[:, 1] >= -0.5) & (xy[:, 1] < PATCH - 0.5)]


def match_count(det, gt, r=R_MATCH):
    """Number of true positives: the optimal one-to-one assignment of detections to ground-truth centres, pairs within r pixels only."""
    if len(det) == 0 or len(gt) == 0:
        return 0
    d = np.hypot(det[:, None, 0] - gt[None, :, 0], det[:, None, 1] - gt[None, :, 1])
    cost = np.where(d <= r, d, 1e6)
    i, j = linear_sum_assignment(cost)
    return int(np.sum(cost[i, j] <= r))


def score_patch(dets, gt):
    """(number of detections inside the patch, true positives) of a detection array (x, y, ...) against the ground-truth centres."""
    d = inside(np.asarray(dets, float)[:, :2]) if len(dets) else np.zeros((0, 2))
    return len(d), match_count(d, gt)


# ----------------------------------------------------------------------------------------------------------------------------------
# statistics
# ----------------------------------------------------------------------------------------------------------------------------------
def pooled(err, idx):
    """Mean absolute error over the patches idx."""
    return float(np.mean(np.abs(err[idx]))) if len(idx) else float('nan')


def boot_diff(ea, eb, img, idx, n_boot=N_BOOT, seed=0):
    """Paired bootstrap over images of mean|ea| - mean|eb| on the patches idx (ratio of sums within a resample; the unit of resampling is the image). Returns (point, lo, hi)."""
    idx = np.asarray(idx)
    d = np.abs(ea[idx]) - np.abs(eb[idx])
    im = img[idx]
    ids, inv = np.unique(im, return_inverse=True)
    s = np.bincount(inv, weights=d, minlength=len(ids))
    n = np.bincount(inv, minlength=len(ids)).astype(float)
    rng = np.random.default_rng(seed)
    k = rng.integers(0, len(ids), size=(n_boot, len(ids)))
    est = s[k].sum(1) / n[k].sum(1)
    return float(s.sum() / n.sum()), float(np.percentile(est, 2.5)), float(np.percentile(est, 97.5))


def verdict(table):
    """The plan's decision rule. table: dict with
    mae_c, mae_o, mae_g (crowded stratum; classical best, ours, greedy), diff_c (point, lo, hi of MAE(C*) - MAE(ours) on the crowded stratum),
    mae_all_c, mae_all_o, diff_all (point, lo, hi of MAE(ours) - MAE(C*) on all patches), diff_og (ours minus greedy on crowded).
    Returns dict(C1, C2, C3, GO, dual_contributes)."""
    c1 = (table['mae_c'] - table['mae_o'] >= M_REL * table['mae_c']) and table['diff_c'][1] > 0.0
    c2 = table['mae_o'] <= table['mae_g'] + M_GREEDY
    c3 = table['diff_all'][2] <= M_NONINF * table['mae_all_c']
    return dict(C1=bool(c1), C2=bool(c2), C3=bool(c3), GO=bool(c1 and c2 and c3), dual_contributes=bool(table['diff_og'][2] < 0.0))


# ----------------------------------------------------------------------------------------------------------------------------------
# running
# ----------------------------------------------------------------------------------------------------------------------------------
_BANK = {}


def _work(arg):
    method, k, group = arg
    t = time.time()
    prep = prepare(_BANK['patches'][k])
    try:
        res = run_group(method, prep, group)
        err = False
    except Exception as e:                      # a failed run counts as no detection and is recorded
        res = {i: (np.zeros((0, 4)), False) for i, _ in group}
        err = True
        sys.stderr.write('FAIL %s patch %d: %r\n' % (method, k, e))
    secs = time.time() - t
    out = []
    gt = gt_of(_BANK, k)
    for i, _ in group:
        d, hit = res[i]
        nd, tp = score_patch(d, gt)
        keep = np.asarray(d, float)[:, :2] if len(d) else np.zeros((0, 2))
        out.append((i, nd, tp, hit, keep))
    return k, out, secs, err


def run_method(method, bank, patch_idx, cfg_idx=None, workers=1, log=print):
    """Run the (selected) configurations of a method on the patches. Returns dict with arrays ndet, tp, hit (n_cfg, n_patch), secs (n_groups, n_patch), dets (per cfg, per patch)."""
    import multiprocessing as mp
    _BANK.clear()
    _BANK.update(bank)
    cfgs = configs_of(method)
    groups = groups_of(method, cfgs, cfg_idx)
    nc, npt = len(cfgs), len(patch_idx)
    ndet = np.zeros((nc, npt), np.int32)
    tp = np.zeros((nc, npt), np.int32)
    hit = np.zeros((nc, npt), np.uint8)
    secs = np.zeros((len(groups), npt))
    fails = np.zeros(npt, np.uint8)
    dets = {}
    pos = {int(k): j for j, k in enumerate(patch_idx)}
    for gi, group in enumerate(groups):
        tasks = [(method, int(k), group) for k in patch_idx]
        pool = mp.get_context('fork').Pool(workers) if workers > 1 else None
        it = pool.imap(_work, tasks, chunksize=1) if pool is not None else map(_work, tasks)
        t0 = time.time()
        for k, out, s, err in it:
            j = pos[k]
            secs[gi, j] = s
            fails[j] |= err
            for i, nd, t_, h, d in out:
                ndet[i, j], tp[i, j], hit[i, j] = nd, t_, h
                dets[(i, j)] = d
        if pool is not None:
            pool.close()
            pool.join()
        log('%s group %d/%d (%d cfg) %.0f s' % (method, gi + 1, len(groups), len(group), time.time() - t0))
    return dict(ndet=ndet, tp=tp, hit=hit, secs=secs, fails=fails, dets=dets, cfgs=cfgs, groups=[[i for i, _ in g] for g in groups])


def pack_dets(dets, cfg_ids, npt):
    """Flatten the detections of the given configurations into (xy (m, 2), offsets (len(cfg_ids) * npt + 1))."""
    chunks, off = [], [0]
    for i in cfg_ids:
        for j in range(npt):
            d = np.asarray(dets.get((i, j), np.zeros((0, 2))), float).reshape(-1, 2)
            chunks.append(d)
            off.append(off[-1] + len(d))
    return (np.concatenate(chunks) if chunks else np.zeros((0, 2))).astype(np.float32), np.array(off, np.int64)


def objective(ndet, n_gt, base, crowd):
    """Development objective of every configuration: the mean of the count MAE over the base sample and over the crowded stratum."""
    err = ndet.astype(float) - n_gt[None, :].astype(float)
    return 0.5 * (np.abs(err[:, base]).mean(1) + np.abs(err[:, crowd]).mean(1))


def stage_dev(a):
    bank = load_bank(a.bank)
    dev = np.flatnonzero(bank['split'] == 0)
    if a.limit:
        dev = dev[:a.limit]
    out = {'dev_patch': dev}
    for m in a.methods.split(','):
        r = run_method(m, bank, dev, workers=a.workers)
        for key in ('ndet', 'tp', 'hit', 'secs', 'fails'):
            out['%s|%s' % (m, key)] = r[key]
        sys.stdout.write('METHOD_DONE %s\n' % m)
        sys.stdout.flush()
        np.savez_compressed(os.path.join(a.out, 'blobcount_dev_raw_%s.npz' % m), **{k: v for k, v in out.items() if k == 'dev_patch' or k.startswith(m + '|')})


def stage_freeze(a):
    """Choose, for each method, the configuration of the lowest development objective from the assembled development raw file; write the freeze file."""
    with np.load(a.dev_raw, allow_pickle=False) as z:
        raw = {k: z[k] for k in z.files}
    ngt, boost = raw['n_gt'], raw['boost']
    base = np.flatnonzero(boost == 0)
    crowd = np.flatnonzero(ngt >= CROWD)
    res = {}
    for m in a.methods.split(','):
        obj = objective(raw['%s|ndet' % m], ngt, base, crowd)
        k = int(np.argmin(obj))
        res[m] = dict(index=k, config=configs_of(m)[k], objective=float(obj[k]), objective_all=[float(x) for x in obj], n_configs=len(obj))
    json.dump(dict(methods=res, n_dev_patches=int(len(ngt))), open(a.freeze, 'w'), indent=1, sort_keys=True)


def stage_test(a):
    bank = load_bank(a.bank)
    fr = json.load(open(a.freeze))['methods']
    test = np.flatnonzero(bank['split'] == 1)
    out = {'test_patch': test}
    for m in a.methods.split(','):
        k = fr[m]['index']
        assert configs_of(m)[k] == fr[m]['config'], 'freeze file does not match the configuration grid of the code'
        r = run_method(m, bank, test, cfg_idx=[k], workers=a.workers)
        xy, off = pack_dets(r['dets'], [k], len(test))
        out.update({'%s|ndet' % m: r['ndet'][k], '%s|tp' % m: r['tp'][k], '%s|hit' % m: r['hit'][k], '%s|secs' % m: r['secs'].sum(0), '%s|fails' % m: r['fails'],
                    '%s|dets' % m: xy, '%s|off' % m: off})
        sys.stdout.write('METHOD_DONE %s\n' % m)
        sys.stdout.flush()
        np.savez_compressed(os.path.join(a.out, 'blobcount_test_raw_%s.npz' % m), **{kk: v for kk, v in out.items() if kk == 'test_patch' or kk.startswith(m + '|')})


# ----------------------------------------------------------------------------------------------------------------------------------
# assembling the raw files and the reports (everything printed is computed from the raw files; `--from-raw` regenerates a report byte-identically)
# ----------------------------------------------------------------------------------------------------------------------------------
META_KEYS = ('n_gt', 'img', 'boost', 'n_touch')


def assemble_dev(meta_path, raw_dir, out, methods=METHODS):
    """Merge the per-method development pieces and the bank's metadata into one raw file."""
    meta = load_bank(meta_path)
    out_d = {}
    for m in methods:
        with np.load(os.path.join(raw_dir, 'dev_%s.npz' % m), allow_pickle=False) as z:
            dp = z['dev_patch']
            for k in ('ndet', 'tp', 'hit', 'secs', 'fails'):
                out_d['%s|%s' % (m, k)] = z['%s|%s' % (m, k)]
    out_d['dev_patch'] = dp
    for k in META_KEYS:
        out_d[k] = meta[k][dp]
    np.savez_compressed(out, **out_d)


def assemble_test(meta_path, raw_dir, freeze_path, info_path, out, methods=METHODS):
    meta = load_bank(meta_path)
    out_d = {}
    for m in methods:
        with np.load(os.path.join(raw_dir, 'test_%s.npz' % m), allow_pickle=False) as z:
            tp_ = z['test_patch']
            for k in ('ndet', 'tp', 'hit', 'secs', 'fails', 'dets', 'off'):
                out_d['%s|%s' % (m, k)] = z['%s|%s' % (m, k)]
    out_d['test_patch'] = tp_
    for k in META_KEYS:
        out_d[k] = meta[k][tp_]
    out_d['freeze'] = np.array(open(freeze_path).read())
    out_d['info'] = np.array(open(info_path).read())
    np.savez_compressed(out, **out_d)


def _f1(ndet, tp, ngt):
    p = tp.sum() / max(ndet.sum(), 1)
    r = tp.sum() / max(ngt.sum(), 1)
    return 2 * p * r / (p + r) if p + r > 0 else 0.0


def _row(name, nd, tp, ngt, idx):
    e = (nd[idx] - ngt[idx]).astype(float)
    pos = ngt[idx] >= 1
    rel = float(np.mean(np.abs(e[pos]) / ngt[idx][pos])) if pos.any() else float('nan')
    return '| %s | %.3f | %.3f | %+.3f | %.3f |' % (name, np.abs(e).mean(), rel, e.mean(), _f1(nd[idx], tp[idx], ngt[idx]))


def dev_report(raw):
    nm = [m for m in METHODS + METHODS_NC if '%s|ndet' % m in raw]
    ngt, boost = raw['n_gt'], raw['boost']
    base = np.flatnonzero(boost == 0)
    crowd = np.flatnonzero(ngt >= CROWD)
    L = ['# Nuclei counting, development run (60 development images)', '',
         'Generated by `experiments/exp_blobcount.py` from the stored development raw file (`--from-raw`). %d development patches: %d base, %d crowded (>= %d nuclei), from %d images.' % (
             len(ngt), len(base), len(crowd), CROWD, len(np.unique(raw['img']))), '',
         'Objective of every configuration: the mean of the count MAE over the base sample and over the crowded stratum (development patches only). The configuration of each method with the lowest objective is frozen.', '',
         '| method | configs | best index | objective | MAE base | MAE crowded | bias base | bias crowded | F1 base | share of paths ending at the cap (best) | seconds per patch (all configs) |', '|---|---|---|---|---|---|---|---|---|---|---|']
    for m in nm:
        nd, tp = raw['%s|ndet' % m], raw['%s|tp' % m]
        obj = objective(nd, ngt, base, crowd)
        k = int(np.argmin(obj))
        eb = nd[k, base] - ngt[base]
        ec = nd[k, crowd] - ngt[crowd]
        L.append('| %s | %d | %d | %.3f | %.3f | %.3f | %+.3f | %+.3f | %.3f | %.3f | %.2f |' % (m, nd.shape[0], k, obj[k], np.abs(eb).mean(), np.abs(ec).mean(), eb.mean(), ec.mean(),
                                                                                              _f1(nd[k, base], tp[k, base], ngt[base]), raw['%s|hit' % m][k].mean(), raw['%s|secs' % m].sum(0).mean()))
    L += ['', '## Objective of every configuration', '']
    for m in nm:
        obj = objective(raw['%s|ndet' % m], ngt, base, crowd)
        L.append('* %s: %s' % (m, ', '.join('%d:%.2f' % (i, o) for i, o in enumerate(obj))))
    return '\n'.join(L) + '\n'


BINS = ((1, 5), (6, 10), (11, 14), (15, 24), (25, 999))


def test_stats(raw):
    """All numbers of the test report, from the raw file: dict."""
    nm = [m for m in METHODS if '%s|ndet' % m in raw]
    ngt, img, boost = raw['n_gt'], raw['img'], raw['boost']
    idx = {'all': np.flatnonzero(boost == 0), 'crowded': np.flatnonzero(ngt >= CROWD), 'very': np.flatnonzero(ngt >= CROWD_HI)}
    err = {m: raw['%s|ndet' % m].astype(float) - ngt.astype(float) for m in nm}
    out = dict(methods=nm, idx=idx, err=err)
    best = {}
    for st in ('all', 'crowded'):
        cl = [m for m in CLASSICAL if m in nm]
        best[st] = min(cl, key=lambda m: pooled(err[m], idx[st]))
    out['best'] = best
    ea, eo = err[best['crowded']], err['ours']
    t = dict(mae_c=pooled(ea, idx['crowded']), mae_o=pooled(eo, idx['crowded']), mae_g=pooled(err['greedy'], idx['crowded']),
             diff_c=boot_diff(ea, eo, img, idx['crowded']), mae_all_c=pooled(err[best['all']], idx['all']), mae_all_o=pooled(eo, idx['all']),
             diff_all=boot_diff(eo, err[best['all']], img, idx['all']), diff_og=boot_diff(eo, err['greedy'], img, idx['crowded']))
    out['table'] = t
    out['verdict'] = verdict(t)
    return out


def test_report(raw):
    S = test_stats(raw)
    nm, idx, err, T, V = S['methods'], S['idx'], S['err'], S['table'], S['verdict']
    ngt, img = raw['n_gt'], raw['img']
    fr = json.loads(str(raw['freeze']))['methods']
    L = ['# Nuclei counting in fluorescence microscopy: the sliding inverse scale space against tuned classical detectors (test run)', '',
         'Generated by `experiments/exp_blobcount.py` from `results/ml/blobcount_raw.npz` (`--from-raw`). Rule and plan: `results/ml/blobcount_plan.md`. Frozen configurations: `results/ml/blobcount_freeze.json`.', '',
         'Test set: %d patches of 128 x 128 from %d images (150 test images; every image contributes up to %d base tiles and up to %d crowded-boost tiles). '
         '"All" = base sample (%d patches, %d images); "crowded" = every sampled patch with >= %d nuclei (%d patches, %d images); "very crowded" >= %d nuclei (%d patches, %d images).' % (
             len(ngt), len(np.unique(img)), PER_IMAGE, BOOST, len(idx['all']), len(np.unique(img[idx['all']])), CROWD, len(idx['crowded']), len(np.unique(img[idx['crowded']])), CROWD_HI,
             len(idx['very']), len(np.unique(img[idx['very']]))), '',
         'Mean true count: all %.2f, crowded %.2f; share of nuclei touching another mask in the crowded stratum %.2f.' % (
             ngt[idx['all']].mean(), ngt[idx['crowded']].mean(), raw['n_touch'][idx['crowded']].sum() / max(ngt[idx['crowded']].sum(), 1)), '']
    for st, title in (('all', 'All test patches (base sample)'), ('crowded', 'Crowded stratum (>= %d nuclei)' % CROWD), ('very', 'Very crowded (>= %d nuclei), descriptive only' % CROWD_HI)):
        L += ['## ' + title, '', '| method | count MAE | relative error | bias (detected - true) | detection F1 (6 px) |', '|---|---|---|---|---|']
        for m in nm:
            L.append(_row(m, raw['%s|ndet' % m], raw['%s|tp' % m], ngt, idx[st]))
        L.append('')
    L += ['## Count MAE by true count', '', '| true count | patches | ' + ' | '.join(nm) + ' |', '|---|---|' + '---|' * len(nm)]
    for lo, hi in BINS:
        s = np.flatnonzero((ngt >= lo) & (ngt <= hi))
        if len(s):
            L.append('| %d-%s | %d | ' % (lo, ('%d' % hi) if hi < 999 else 'up', len(s)) + ' | '.join('%.2f' % pooled(err[m], s) for m in nm) + ' |')
    L += ['', '## Paired bootstrap over images (10 000 resamples, 95 % percentile interval) of the difference in count MAE, ours minus the other method (negative: ours better)', '',
          '| other method | crowded: difference | interval | all: difference | interval |', '|---|---|---|---|---|']
    for m in nm:
        if m == 'ours':
            continue
        c = boot_diff(err['ours'], err[m], img, idx['crowded'])
        a = boot_diff(err['ours'], err[m], img, idx['all'])
        L.append('| %s | %+.3f | [%+.3f, %+.3f] | %+.3f | [%+.3f, %+.3f] |' % (m, c[0], c[1], c[2], a[0], a[1], a[2]))
    bc, ba = S['best']['crowded'], S['best']['all']
    L += ['', '## Verdict by the pre-registered rule', '',
          'Best classical method (lowest test MAE in the stratum, an oracle choice favouring the baselines): crowded %s, all %s.' % (bc, ba), '',
          '* C1 (crowded gain >= %d %% of the best classical MAE with the interval above 0): MAE(%s) = %.3f, MAE(ours) = %.3f, gain %.3f (%.1f %%), interval of MAE(%s) - MAE(ours) [%+.3f, %+.3f] -> **%s**' % (
              round(100 * M_REL), bc, T['mae_c'], T['mae_o'], T['mae_c'] - T['mae_o'], 100 * (T['mae_c'] - T['mae_o']) / T['mae_c'], bc, T['diff_c'][1], T['diff_c'][2], 'met' if V['C1'] else 'not met'),
          '* C2 (ours within %.1f nuclei of greedy on the crowded stratum): MAE(ours) = %.3f, MAE(greedy) = %.3f -> **%s**' % (M_GREEDY, T['mae_o'], T['mae_g'], 'met' if V['C2'] else 'not met'),
          '* C3 (all patches: upper end of the interval of MAE(ours) - MAE(%s) at most %d %% of MAE(%s)): MAE(ours) = %.3f, MAE(%s) = %.3f, difference %+.3f, interval [%+.3f, %+.3f], limit %.3f -> **%s**' % (
              ba, round(100 * M_NONINF), ba, T['mae_all_o'], ba, T['mae_all_c'], T['diff_all'][0], T['diff_all'][1], T['diff_all'][2], M_NONINF * T['mae_all_c'], 'met' if V['C3'] else 'not met'),
          '* Dual order contributes (ours minus greedy on the crowded stratum, interval below 0): %s (difference %+.3f, interval [%+.3f, %+.3f]).' % (
              'yes' if V['dual_contributes'] else 'no', T['diff_og'][0], T['diff_og'][1], T['diff_og'][2]), '',
          '**Verdict: %s**' % ('GO' if V['GO'] else 'KILL'), '', '## Diagnostics', '',
          '| method | seconds per patch (mean) | failed runs | paths whose BIC block is the last one |', '|---|---|---|---|']
    for m in nm:
        L.append('| %s | %.2f | %d | %s |' % (m, raw['%s|secs' % m].mean(), int(raw['%s|fails' % m].sum()), '%.3f' % raw['%s|hit' % m].mean() if m in SLIDING else '-'))
    L += ['', '## Frozen configurations (chosen on the development patches)', '']
    for m in nm:
        L.append('* %s (grid index %d of %d, development objective %.3f): `%s`' % (m, fr[m]['index'], fr[m]['n_configs'], fr[m]['objective'], json.dumps(fr[m]['config'], sort_keys=True)))
    return '\n'.join(L) + '\n'


# ----------------------------------------------------------------------------------------------------------------------------------
# POST HOC (after the test result): is the failure the path or the stop rule? An oracle stop along the frozen paths.
# ----------------------------------------------------------------------------------------------------------------------------------
def oracle_patch(prep, path, order, gt, stop):
    """Along the path of one patch: the count of every block (detections inside the patch), the true positives of every block, and the block of the frozen stop.
    Returns dict(counts, tps, chosen)."""
    P, R = sliding_run(prep, path, order)
    counts, tps = [], []
    for h in R['hist']:
        d = inside(h['c'][:, :2]) if len(h['w']) else np.zeros((0, 2))
        counts.append(len(d))
        tps.append(match_count(d, gt))
    _, i, _ = sliding_stop(prep, P, R, stop)
    return dict(counts=np.array(counts), tps=np.array(tps), chosen=i)


def _oracle_work(arg):
    method, k, cfg = arg
    prep = prepare(_BANK['patches'][k])
    o = oracle_patch(prep, cfg['path'], order_of(method), gt_of(_BANK, k), cfg['stop'])
    n = int(_BANK['n_gt'][k])
    err = np.abs(o['counts'] - n)
    j = int(np.argmin(err))
    return k, (int(err[j]), int(o['tps'][j]), int(o['tps'].max()), int(o['counts'][o['chosen']]), len(o['counts']), int(o['counts'][-1]), int(o['counts'][j]))


def stage_oracle(a):
    import multiprocessing as mp
    bank = load_bank(a.bank)
    fr = json.load(open(a.freeze))['methods']
    test = np.flatnonzero(bank['split'] == 1)
    _BANK.clear()
    _BANK.update(bank)
    out = {'test_patch': test}
    for m in ('ours', 'greedy'):
        pool = mp.get_context('fork').Pool(a.workers)
        res = dict(pool.imap(_oracle_work, [(m, int(k), fr[m]['config']) for k in test], chunksize=1))
        pool.close()
        pool.join()
        out['%s|oracle' % m] = np.array([res[int(k)] for k in test], np.int32)       # columns: oracle |error|, tp at that block, max tp, count at the frozen stop, blocks, count at the end, count at the oracle block
        sys.stdout.write('ORACLE_DONE %s\n' % m)
        sys.stdout.flush()
    np.savez_compressed(os.path.join(a.out, 'blobcount_posthoc_raw.npz'), **out)


def posthoc_report(raw, test_raw):
    """The POST HOC diagnostic, from the oracle raw (and the test raw for the strata and the actual counts)."""
    ngt, img, boost = test_raw['n_gt'], test_raw['img'], test_raw['boost']
    idx = {'all': np.flatnonzero(boost == 0), 'crowded': np.flatnonzero(ngt >= CROWD)}
    L = ['# POST HOC diagnostic (made after the test result): the frozen paths with an oracle stop', '',
         'The frozen path of `ours` (dual order) and of `greedy` is run on the test patches again; the block of the path whose count is closest to the true count (it reads the ground truth: a diagnostic, not a method) is taken. '
         'This separates "the path does not contain the right number of units" from "the stop rule picks the wrong block". Generated from `results/ml/blobcount_posthoc_raw.npz` and `results/ml/blobcount_raw.npz` (`--from-raw`).', '',
         '| method | stratum | count MAE at the frozen stop | count MAE at the oracle block | recall of true centres at the oracle block | best recall anywhere on the path | mean blocks on the path | MAE of the last block |', '|---|---|---|---|---|---|---|---|']
    for m in ('ours', 'greedy'):
        o = raw['%s|oracle' % m].astype(float)
        for st in ('all', 'crowded'):
            i = idx[st]
            act = np.abs(test_raw['%s|ndet' % m][i] - ngt[i]).mean()
            rec = o[i, 1].sum() / max(ngt[i].sum(), 1)
            recmax = o[i, 2].sum() / max(ngt[i].sum(), 1)
            L.append('| %s | %s | %.3f | %.3f | %.3f | %.3f | %.1f | %.3f |' % (m, st, act, o[i, 0].mean(), rec, recmax, o[i, 4].mean(), np.abs(o[i, 5] - ngt[i]).mean()))
    c = {m: raw['%s|oracle' % m][:, 0].astype(float) for m in ('ours', 'greedy')}
    d = boot_diff(c['ours'], c['greedy'], img, idx['crowded'])
    L += ['', 'Oracle-stop crowded MAE, ours minus greedy: %+.3f, interval [%+.3f, %+.3f] (paired bootstrap over images).' % d, '']
    return '\n'.join(L) + '\n'


# ----------------------------------------------------------------------------------------------------------------------------------
# REFERENCE (run after the test result, pretrained models with default parameters, nothing tuned, not a target of the verdict)
# ----------------------------------------------------------------------------------------------------------------------------------
REFS = ('cellpose', 'stardist')


def stage_ref(a):
    """Pretrained Cellpose ('nuclei', default diameter) and StarDist ('2D_versatile_fluo', default thresholds) on the test patches; the centroids of the instances are scored like any detector."""
    import importlib.metadata as md
    bank = load_bank(a.bank)
    test = np.flatnonzero(bank['split'] == 1)
    out = {'test_patch': test}
    versions = {}
    for name in REFS:
        det = []
        t0 = time.time()
        if name == 'cellpose':
            from cellpose import models
            versions[name] = md.version('cellpose')
            model = models.CellposeModel(gpu=False, model_type='nuclei')
            for k in test:
                im = bank['patches'][k]
                masks = model.eval(im, channels=[0, 0])[0]
                ids = [i for i in np.unique(masks) if i > 0]
                cm = np.array(ndi.center_of_mass(masks > 0, masks, ids)) if ids else np.zeros((0, 2))
                det.append(np.c_[cm[:, 1], cm[:, 0]] if len(cm) else np.zeros((0, 2)))
        else:
            from csbdeep.utils import normalize
            from stardist.models import StarDist2D
            versions[name] = md.version('stardist')
            model = StarDist2D.from_pretrained('2D_versatile_fluo')
            for k in test:
                lab, info = model.predict_instances(normalize(bank['patches'][k].astype(np.float32), 1, 99.8), verbose=False, show_tile_progress=False)
                pts = np.asarray(info['points'], float).reshape(-1, 2)
                det.append(np.c_[pts[:, 1], pts[:, 0]] if len(pts) else np.zeros((0, 2)))
        sc = [score_patch(d, gt_of(bank, k)) for d, k in zip(det, test)]
        out['%s|ndet' % name] = np.array([x[0] for x in sc], np.int32)
        out['%s|tp' % name] = np.array([x[1] for x in sc], np.int32)
        out['%s|secs' % name] = np.full(len(test), (time.time() - t0) / len(test))
        sys.stdout.write('REF_DONE %s %s %.0f s\n' % (name, versions[name], time.time() - t0))
        sys.stdout.flush()
    out['versions'] = np.array(json.dumps(versions))
    np.savez_compressed(os.path.join(a.out, 'blobcount_ref_raw.npz'), **out)


def ref_report(ref, test_raw):
    """REFERENCE table from the reference raw and the test raw."""
    ngt, img, boost = test_raw['n_gt'], test_raw['img'], test_raw['boost']
    idx = {'all': np.flatnonzero(boost == 0), 'crowded': np.flatnonzero(ngt >= CROWD), 'very crowded': np.flatnonzero(ngt >= CROWD_HI)}
    S = test_stats(test_raw)
    L = ['# Reference: pretrained instance-segmentation models on the same test patches (run after the test result)', '',
         'Pretrained models with their default parameters, nothing tuned or chosen on any patch, centroids of the instances scored exactly like the detectors of `blobcount.md`. Versions: %s. '
         'These are references, not targets: such models are reportedly trained on data that include nuclei images of this challenge (no check was made of any overlap with the test images), so their numbers are probably optimistic. '
         'Generated from `results/ml/blobcount_ref_raw.npz` and `results/ml/blobcount_raw.npz` (`--from-raw`).' % ', '.join('%s %s' % kv for kv in sorted(json.loads(str(ref['versions'])).items())), '']
    for st in ('all', 'crowded', 'very crowded'):
        L += ['## ' + st, '', '| method | count MAE | relative error | bias (detected - true) | detection F1 (6 px) |', '|---|---|---|---|---|']
        L.append(_row('ours', test_raw['ours|ndet'], test_raw['ours|tp'], ngt, idx[st]))
        L.append(_row(S['best']['crowded' if st != 'all' else 'all'] + ' (best classical)', test_raw['%s|ndet' % S['best']['crowded' if st != 'all' else 'all']],
                      test_raw['%s|tp' % S['best']['crowded' if st != 'all' else 'all']], ngt, idx[st]))
        for m in REFS:
            if '%s|ndet' % m in ref:
                L.append(_row(m, ref['%s|ndet' % m], ref['%s|tp' % m], ngt, idx[st]))
        L.append('')
    return '\n'.join(L)


def nocap_report(dev_raw, nc_raw, main_raw):
    """POST HOC: the sliding methods without the early stop of the path (grids otherwise identical, re-tuned on the development patches, frozen, run once on the test patches), put into the
    pre-registered rule in place of `ours` / `greedy`. dev_raw: the development raw of the variants; nc_raw: their test raw; main_raw: the test raw of the pre-registered run."""
    ngt, img, boost = main_raw['n_gt'], main_raw['img'], main_raw['boost']
    idx = {'all': np.flatnonzero(boost == 0), 'crowded': np.flatnonzero(ngt >= CROWD)}
    fr = json.loads(str(nc_raw['freeze']))['methods']
    L = ['# POST HOC: the sliding methods without the early stop of the path', '',
         'Made after the test result of the pre-registered run (`blobcount.md`, verdict KILL). The path of `ours` and `greedy` was cut short by a stop that was fixed in the plan and not tuned (a path stops a few entries after an entrant of matched-filter significance below 2); on crowded patches the cut path '
         'ended several units short of the true count (see `blobcount_posthoc.md`, column "MAE of the last block"). Here the same 24-configuration grids are run with that stop removed (the entry cap of 120 stays), re-tuned on the development patches with the same objective, frozen, and run once on the test patches. '
         'Everything here is POST HOC: it was not pre-registered and cannot support a claim without new held-out data.', '',
         '## Development (objective: mean of the base and crowded MAE)', '', '| method | best index | objective | MAE base | MAE crowded | bias base | bias crowded |', '|---|---|---|---|---|---|---|']
    dn, db, dc = dev_raw['n_gt'], np.flatnonzero(dev_raw['boost'] == 0), np.flatnonzero(dev_raw['n_gt'] >= CROWD)
    for m in METHODS_NC:
        nd = dev_raw['%s|ndet' % m]
        obj = objective(nd, dn, db, dc)
        k = int(np.argmin(obj))
        L.append('| %s | %d | %.3f | %.3f | %.3f | %+.3f | %+.3f |' % (m, k, obj[k], np.abs(nd[k, db] - dn[db]).mean(), np.abs(nd[k, dc] - dn[dc]).mean(), (nd[k, db] - dn[db]).mean(), (nd[k, dc] - dn[dc]).mean()))
    merged = dict(main_raw)
    for m, tag in (('ours_nc', 'ours'), ('greedy_nc', 'greedy')):
        for k in ('ndet', 'tp', 'hit', 'secs', 'fails'):
            merged['%s|%s' % (tag, k)] = nc_raw['%s|%s' % (m, k)]
    S = test_stats(merged)
    T, V = S['table'], S['verdict']
    L += ['', '## Test (the frozen configurations of the variants)', '', '| method | stratum | count MAE | relative error | bias | F1 |', '|---|---|---|---|---|---|']
    for st in ('all', 'crowded'):
        for m in ('ours', 'greedy'):
            L.append('| %s (pre-registered) | %s ' % (m, st) + _row('', main_raw['%s|ndet' % m], main_raw['%s|tp' % m], ngt, idx[st])[3:])
        for m in METHODS_NC:
            L.append('| %s | %s ' % (m, st) + _row('', nc_raw['%s|ndet' % m], nc_raw['%s|tp' % m], ngt, idx[st])[3:])
        b = S['best'][st]
        L.append('| %s (best classical) | %s ' % (b, st) + _row('', main_raw['%s|ndet' % b], main_raw['%s|tp' % b], ngt, idx[st])[3:])
    L += ['', 'Paired bootstrap over images of the difference in count MAE (variant minus other):', '', '| comparison | crowded | all |', '|---|---|---|']
    err = {m: nc_raw['%s|ndet' % m].astype(float) - ngt for m in METHODS_NC}
    for m in METHODS_NC:
        for o in S['methods']:
            if o in ('ours', 'greedy') or o in METHODS_NC:
                continue
            c = boot_diff(err[m], S['err'][o], img, idx['crowded'])
            a = boot_diff(err[m], S['err'][o], img, idx['all'])
            L.append('| %s - %s | %+.3f [%+.3f, %+.3f] | %+.3f [%+.3f, %+.3f] |' % (m, o, c[0], c[1], c[2], a[0], a[1], a[2]))
    c = boot_diff(err['ours_nc'], err['greedy_nc'], img, idx['crowded'])
    a = boot_diff(err['ours_nc'], err['greedy_nc'], img, idx['all'])
    L.append('| ours_nc - greedy_nc | %+.3f [%+.3f, %+.3f] | %+.3f [%+.3f, %+.3f] |' % (c + a))
    L += ['', 'The pre-registered rule with `ours_nc` and `greedy_nc` in place of `ours` and `greedy` (best classical: crowded %s, all %s): C1 %s, C2 %s, C3 %s, rule outcome **%s** (POST HOC, not the verdict).' % (
        S['best']['crowded'], S['best']['all'], 'met' if V['C1'] else 'not met', 'met' if V['C2'] else 'not met', 'met' if V['C3'] else 'not met', 'GO' if V['GO'] else 'KILL'), '',
        'Crowded MAE: ours_nc %.3f, best classical %.3f, greedy_nc %.3f; all-patch MAE: ours_nc %.3f, best classical %.3f. Share of paths whose BIC block is the last one: ours_nc %.3f, greedy_nc %.3f.' % (
            T['mae_o'], T['mae_c'], T['mae_g'], T['mae_all_o'], T['mae_all_c'], nc_raw['ours_nc|hit'].mean(), nc_raw['greedy_nc|hit'].mean()), '']
    for m in METHODS_NC:
        L.append('* %s (grid index %d, development objective %.3f): `%s`' % (m, fr[m]['index'], fr[m]['objective'], json.dumps(fr[m]['config'], sort_keys=True)))
    return '\n'.join(L) + '\n'


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--stage', choices=['bank', 'dev', 'freeze', 'test', 'assemble_dev', 'assemble_test', 'oracle', 'ref'])
    ap.add_argument('--from-raw', nargs='?', const='')
    ap.add_argument('--md')
    ap.add_argument('--meta')
    ap.add_argument('--info')
    ap.add_argument('--zip')
    ap.add_argument('--bank')
    ap.add_argument('--out', default='.')
    ap.add_argument('--methods', default=','.join(METHODS))
    ap.add_argument('--workers', type=int, default=1)
    ap.add_argument('--limit', type=int, default=0)
    ap.add_argument('--freeze')
    ap.add_argument('--dev_raw')
    ap.add_argument('--raw_dir', default='.')
    a = ap.parse_args(argv)
    if a.from_raw is not None:
        outd = os.environ.get('BLOBCOUNT_OUT', os.path.join(REPO, 'results', 'ml'))
        jobs = [(a.from_raw, a.md)] if a.from_raw else [(os.path.join(outd, 'blobcount_raw.npz'), os.path.join(outd, 'blobcount.md')),
                                                         (os.path.join(outd, 'blobcount_dev_raw.npz'), os.path.join(outd, 'blobcount_dev.md'))]
        for rawp, mdp in jobs:
            if not os.path.exists(rawp):
                continue
            with np.load(rawp, allow_pickle=False) as z:
                raw = {k: z[k] for k in z.files}
            open(mdp, 'w').write(test_report(raw) if 'test_patch' in raw else dev_report(raw))
        ph, tr = os.path.join(outd, 'blobcount_posthoc_raw.npz'), os.path.join(outd, 'blobcount_raw.npz')
        if not a.from_raw and os.path.exists(ph) and os.path.exists(tr):
            with np.load(ph, allow_pickle=False) as z, np.load(tr, allow_pickle=False) as y:
                open(os.path.join(outd, 'blobcount_posthoc.md'), 'w').write(posthoc_report({k: z[k] for k in z.files}, {k: y[k] for k in y.files}))
        nd_, nt_ = os.path.join(outd, 'blobcount_nocap_dev_raw.npz'), os.path.join(outd, 'blobcount_nocap_raw.npz')
        if not a.from_raw and os.path.exists(nd_) and os.path.exists(nt_) and os.path.exists(tr):
            with np.load(nd_, allow_pickle=False) as z, np.load(nt_, allow_pickle=False) as y, np.load(tr, allow_pickle=False) as w:
                open(os.path.join(outd, 'blobcount_nocap.md'), 'w').write(nocap_report({k: z[k] for k in z.files}, {k: y[k] for k in y.files}, {k: w[k] for k in w.files}))
        rf = os.path.join(outd, 'blobcount_ref_raw.npz')
        if not a.from_raw and os.path.exists(rf) and os.path.exists(tr):
            with np.load(rf, allow_pickle=False) as z, np.load(tr, allow_pickle=False) as y:
                open(os.path.join(outd, 'blobcount_ref.md'), 'w').write(ref_report({k: z[k] for k in z.files}, {k: y[k] for k in y.files}))
        return
    if a.stage == 'assemble_dev':
        assemble_dev(a.meta, a.raw_dir, a.out, methods=tuple(a.methods.split(',')))
    elif a.stage == 'assemble_test':
        assemble_test(a.meta, a.raw_dir, a.freeze, a.info, a.out, methods=tuple(a.methods.split(',')))
    elif a.stage == 'oracle':
        stage_oracle(a)
    elif a.stage == 'ref':
        stage_ref(a)
    elif a.stage == 'bank':
        arrays, info = build_bank(a.zip)
        np.savez_compressed(os.path.join(a.out, 'blobcount_bank.npz'), **arrays)
        json.dump(info, open(os.path.join(a.out, 'blobcount_bank_info.json'), 'w'), indent=1)
    elif a.stage == 'dev':
        stage_dev(a)
    elif a.stage == 'freeze':
        stage_freeze(a)
    else:
        stage_test(a)


if __name__ == '__main__':
    main()
