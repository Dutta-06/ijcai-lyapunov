"""
Calibration of the simulated 2D localization benchmark
======================================================

Measures what a localization method needs to know about the data of a set (default MT0.N1.LD, with MT0.N2.LD as the second noise point):
 (0) the conventions the ground truth follows (frame base, x/y axes, pixel-centre convention), checked against the images;
 (a) the empirical PSF width: a 2D Gaussian + constant fitted to isolated emitters (point-sampled and pixel-integrated; the fit window is printed with every
     width, and the dependence of the defocus widths on the window), and the offset between the centre of a symmetric fit and the ground-truth position,
     which depends on the estimator and on |z| (table of estimators);
 (b) the background level and noise on the QUIET frames (no ground-truth emitter in the frame or in either neighbour), the PSF halo that raises the
     background around emitters, and the radial excess profile around isolated single-emitter frames (cumulative fraction of the ground-truth photons
     within 3, 5, 10, 25 px; the abrupt cut-off);
 (c) the noise model var = a x + b in photon units x = (DN - 100) / 6: a and b from the two backgrounds (N1 and N2), and a check on emitter pixels;
 (d) the ground-truth emitters per frame, the photons per emitter and the persistence of an emitter over consecutive frames.

Frames: the PSF, offset and emitter statistics use the development frames of the protocol (default 0-999) only. The camera constants of (b) and (c) and the
halo profile use the quiet / single-emitter frames of the stack OUTSIDE the test frames (every frame except 1000-1499; the development-only numbers are printed
next to them): they are properties of the camera and the simulated PSF, the ground truth only selects the frames, and no test frame is used anywhere.

Output: results/smlm/calibration.md (regenerates byte for byte; no timings). Usage: python experiments/exp_smlm_calibration.py [--set MT0.N1.LD]
[--second-set MT0.N2.LD] [--out PATH]
"""
import argparse
import json
import math
import os
import sys

import numpy as np
from scipy.optimize import least_squares
from scipy.special import erf

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
import smlm_common as C

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
HALF = 5                       # patch half-size of the main fits (11 x 11 px)
WINDOW_HALVES = (3, 4, 5, 6, 9)                # fit windows of the window study: 7 x 7, 9 x 9, 11 x 11, 13 x 13, 19 x 19 px
ISO_PX = 8.0                   # isolation radius (px): no other ground-truth emitter of the frame closer than this
ISO_WINDOW_PX = 14.0           # isolation radius of the window study (the 19 x 19 window reaches 12.7 px at its corners)
MIN_PHOTONS = 2000.0           # emitters brighter than this enter the PSF fits (reported against photons below)
EXPECTED_SIGMA_PX = 0.21 * 660.0 / 1.49 / C.PIXEL_NM
Z_BINS = [(0, 100), (100, 200), (200, 350), (350, 800)]
HALO_EDGES = np.array(list(range(0, 13)) + [14, 16, 18, 20, 22, 23, 24, 25, 26, 28, 30, 32, 36, 40, 50, 100], dtype=float)      # radial annuli (px)
HALO_RADII = (3, 5, 10, 25)
CAL_EXCLUDE = C.TEST_FRAMES    # the test frames are never used
GROUPS = 10                    # contiguous groups of frames for the standard errors of the camera constants
SIM_PAGE = 'https://srm.epfl.ch/srm/dataset/challenge-3D-simulation/%s/index.html'


# ------------------------------------------------------------------ fits
_GRIDS = {}


def _grid(half):
    if half not in _GRIDS:
        _GRIDS[half] = np.mgrid[-half:half + 1, -half:half + 1].astype(float)
    return _GRIDS[half]


def _int_profile(ax, c, s):
    """Gaussian integrated over unit pixels centred at integer ``ax`` (1-D profile, sums to ~1)."""
    k = 1.0 / (np.sqrt(2.0) * s)
    return 0.5 * (erf((ax + 0.5 - c) * k) - erf((ax - 0.5 - c) * k))


def _model(theta, kind, half=HALF):
    yy, xx = _grid(half)
    if kind == 'point':                     # [A, x0, y0, s, bg], Gaussian sampled at the pixel centres
        A, x0, y0, s, bg = theta
        return bg + A * np.exp(-0.5 * ((xx - x0) ** 2 + (yy - y0) ** 2) / s ** 2)
    if kind == 'ellip':                     # [A, x0, y0, sx, sy, bg], point-sampled elliptical
        A, x0, y0, sx, sy, bg = theta
        return bg + A * np.exp(-0.5 * ((xx - x0) ** 2 / sx ** 2 + (yy - y0) ** 2 / sy ** 2))
    A, x0, y0, s, bg = theta                # 'int': [N, x0, y0, s, bg], integrated over each pixel, N = total counts
    return bg + A * np.outer(_int_profile(yy[:, 0], y0, s), _int_profile(xx[0], x0, s))


def fit_patch(patch, kind, half=HALF):
    """Least-squares fit (unweighted, DN) of a 2D Gaussian + constant to a (2 half + 1)^2 patch; returns (theta, model image, leverage image) or None."""
    yy, xx = _grid(half)
    bg0 = float(np.median(patch))
    w = np.clip(patch - bg0, 0, None)
    tot = w.sum()
    if tot <= 0:
        return None
    x0, y0 = (w * xx).sum() / tot, (w * yy).sum() / tot
    s0 = 1.2
    amp = float(patch.max() - bg0)
    if kind == 'int':
        th0, lo, hi = [amp * 2 * np.pi * s0 ** 2, x0, y0, s0, bg0], [0, -3, -3, 0.3, 0], [np.inf, 3, 3, 5, np.inf]
    elif kind == 'ellip':
        th0, lo, hi = [amp, x0, y0, s0, s0, bg0], [0, -3, -3, 0.3, 0.3, 0], [np.inf, 3, 3, 5, 5, np.inf]
    else:
        th0, lo, hi = [amp, x0, y0, s0, bg0], [0, -3, -3, 0.3, 0], [np.inf, 3, 3, 5, np.inf]
    try:
        r = least_squares(lambda t: (_model(t, kind, half) - patch).ravel(), th0, bounds=(lo, hi), method='trf', x_scale='jac')
        q, _ = np.linalg.qr(r.jac)
    except Exception:
        return None
    return r.x, _model(r.x, kind, half), (q ** 2).sum(axis=1).reshape(patch.shape)           # parameters, model image, leverage of every pixel


def _spread(v):
    v = np.asarray(v, float)
    q = np.percentile(v, [5, 25, 50, 75, 95])
    return dict(n=len(v), median=q[2], mean=float(v.mean()), std=float(v.std(ddof=1)), iqr=float(q[3] - q[1]), mad_sigma=float(1.4826 * np.median(np.abs(v - q[2]))),
                p5=q[0], p95=q[4])


def isolated_emitters(gt, n_frames, min_photons, half=HALF, iso_px=ISO_PX):
    """Indices into ``gt`` (frames 0-based from 0 .. n_frames - 1 of the array given) of emitters with a full (2 half + 1)^2 patch inside the image and no other
    emitter of the same frame within ``iso_px``."""
    keep = []
    xs, ys, fr = gt['x_nm'] / C.PIXEL_NM, gt['y_nm'] / C.PIXEL_NM, gt['frame']
    for i in range(len(gt)):
        if gt['photons'][i] < min_photons:
            continue
        ix, iy = int(np.floor(xs[i])), int(np.floor(ys[i]))          # pixel containing the emitter (corner convention)
        if ix < half or ix > 63 - half or iy < half or iy > 63 - half:
            continue
        same = np.flatnonzero((fr == fr[i]) & (np.arange(len(gt)) != i))
        if len(same) and np.min(np.hypot(xs[same] - xs[i], ys[same] - ys[i])) < iso_px:
            continue
        keep.append(i)
    return np.array(keep, dtype=int)


# ------------------------------------------------------------------ report helpers
def md_table(header, rows):
    out = ['| ' + ' | '.join(header) + ' |', '|' + '|'.join(['---'] * len(header)) + '|']
    out += ['| ' + ' | '.join(str(c) for c in r) + ' |' for r in rows]
    return '\n'.join(out)


def f(x, nd=2):
    return 'nan' if x is None or (isinstance(x, float) and np.isnan(x)) else ('%.*f' % (nd, x))


def ols(x, y, w=None):
    """Weighted least squares y = a x + b; returns a, b, SE(a), SE(b) from the residual variance."""
    X = np.stack([x, np.ones_like(x)], axis=1)
    w = np.ones_like(y) if w is None else w
    sw = np.sqrt(w)
    beta, *_ = np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)
    res = (y - X @ beta) * sw
    dof = max(len(y) - 2, 1)
    cov = np.linalg.inv((X * w[:, None]).T @ X) * (res @ res) / dof
    return beta[0], beta[1], float(np.sqrt(cov[0, 0])), float(np.sqrt(cov[1, 1]))


def load_stack(name):
    """(uint16 stack (n, 64, 64) of every frame, ground truth of the set)."""
    st = C.load_frames(C.data_path(name, 'stack'), 0, None).astype(np.uint16)
    return st, C.load_gt(C.data_path(name, 'gt'))


# ------------------------------------------------------------------ quiet frames, background, noise
def quiet_frames(gt, n):
    """Boolean mask of the frames with no ground-truth row in the frame or in its two neighbours (a missing neighbour counts as empty)."""
    cnt = np.bincount(gt['frame'], minlength=n)[:n]
    empty = cnt == 0
    q = empty.copy()
    q[1:] &= empty[:-1]
    q[:-1] &= empty[1:]
    return q


def calibration_frames(n):
    """Indices of every frame of the stack outside the test frames."""
    k = np.arange(n)
    return k[(k < CAL_EXCLUDE[0]) | (k >= CAL_EXCLUDE[1])]


def quiet_stats(stack, idx):
    """Background and noise from the frames ``idx`` (quiet frames), in photons x = (DN - 100) / 6. Standard errors from GROUPS contiguous groups of frames."""
    ph = (stack[idx].astype(np.float64) - C.BASELINE_DN) / C.DN_PER_PHOTON
    m = len(idx)
    mean = float(ph.mean())
    sd = float(ph.std())
    pv = ph.var(axis=0, ddof=1)
    var_t = float(pv.mean())                                                    # mean over pixels of the temporal variance
    fm = ph.mean(axis=(1, 2))
    yy, xx = np.mgrid[0:64, 0:64]
    rad = np.hypot(xx + 0.5 - 32, yy + 0.5 - 32)
    radial = []
    for lo, hi in ((0, 8), (8, 16), (16, 24), (24, 32), (32, 40), (40, 46)):
        sel = (rad >= lo) & (rad < hi)
        v = ph[:, sel]
        radial.append((lo, hi, float(v.mean()), float(v.std() / math.sqrt(v.size))))
    chunks = np.array_split(np.arange(m), GROUPS)
    gm = np.array([ph[c].mean() for c in chunks])
    gv = np.array([ph[c].var(axis=0, ddof=1).mean() for c in chunks])
    return dict(m=m, mean=mean, median=float(np.median(ph)), sd=sd, skew=float(((ph - mean) ** 3).mean() / sd ** 3), var=var_t, var_over_mean=var_t / mean,
                se_mean=float(gm.std(ddof=1) / math.sqrt(GROUPS)), se_var=float(gv.std(ddof=1) / math.sqrt(GROUPS)), frame_sd=float(fm.std()),
                expect_sd=math.sqrt(var_t / 4096.0), trend=float(np.polyfit(idx, fm, 1)[0] * 1000) if m > 2 else float('nan'),
                spatial_sd=float(pm_sd(ph)), spatial_expect=math.sqrt(var_t / m), radial=radial, gm=gm, gv=gv)


def pm_sd(ph):
    return ph.mean(axis=0).std()


def halo_leak(stack, gt, idx, base):
    """Frames with emitters: the mean of the 'emitter-free' pixels (farther than ISO_PX from every ground-truth emitter of the frame) against the quiet-frame
    background, and its regression on the total ground-truth photons of the frame (frames ``idx`` that have at least one emitter)."""
    yy, xx = np.mgrid[0:64, 0:64]
    xc, yc = (xx + 0.5) * C.PIXEL_NM, (yy + 0.5) * C.PIXEL_NM
    order = np.argsort(gt['frame'], kind='stable')
    g = gt[order]
    lo = np.searchsorted(g['frame'], idx, side='left')
    hi = np.searchsorted(g['frame'], idx, side='right')
    fm, tot, nfree = [], [], []
    for k, f_ in enumerate(idx):
        rows = g[lo[k]:hi[k]]
        if len(rows) == 0:
            continue
        free = np.ones((64, 64), bool)
        for r in rows:
            free &= np.hypot(xc - r['x_nm'], yc - r['y_nm']) >= ISO_PX * C.PIXEL_NM
        fm.append((stack[f_][free].astype(np.float64).mean() - C.BASELINE_DN) / C.DN_PER_PHOTON)
        tot.append(float(rows['photons'].sum()))
        nfree.append(int(free.sum()))
    fm, tot, nfree = np.array(fm), np.array(tot), np.array(nfree)
    slope, icpt = np.polyfit(tot, fm, 1)
    res = fm - (slope * tot + icpt)
    return dict(n=len(fm), mean=float(fm.mean()), excess=float(fm.mean() - base), sd=float(fm.std()), slope=float(slope), icpt=float(icpt), corr=float(np.corrcoef(tot, fm)[0, 1]),
                res_sd=float(res.std()), nfree=float(nfree.mean()), leaked=float(slope * nfree.mean()))


def halo_profile(stack, gt, idx, base, min_photons):
    """Radial excess profile around isolated single-emitter frames: frames of ``idx`` with exactly one ground-truth row of at least ``min_photons`` photons, more than
    16 px from every border. Returns, per class of emitters (all, in focus |z| < 100 nm, defocused 300 <= |z| < 700 nm), the mean excess (photons - base) per pixel
    in the annuli HALO_EDGES, the pixel counts, the mean ground-truth photons, the cumulative fraction of the ground-truth photons within HALO_RADII (annuli
    of the profile times their full area, i.e. assuming isotropy) with a standard error from GROUPS groups of emitters, and the directly summed excess of the
    whole image per ground-truth photon."""
    yy, xx = np.mgrid[0:64, 0:64]
    cnt = np.bincount(gt['frame'], minlength=len(stack))[:len(stack)]
    by_frame = {int(r['frame']): r for r in gt}
    classes = {'all': lambda z: np.ones_like(z, bool), 'in focus (|z| < 100 nm)': lambda z: np.abs(z) < 100, 'defocused (300 <= |z| < 700 nm)': lambda z: (np.abs(z) >= 300) & (np.abs(z) < 700)}
    nb = len(HALO_EDGES) - 1
    acc = {c: dict(S=np.zeros((GROUPS, nb)), N=np.zeros((GROUPS, nb)), P=np.zeros(GROUPS), n=np.zeros(GROUPS, int), tot=np.zeros(GROUPS)) for c in classes}
    k = 0
    for f_ in idx:
        if cnt[f_] != 1:
            continue
        r = by_frame[int(f_)]
        if r['photons'] < min_photons or not (1600 < r['x_nm'] < 4800 and 1600 < r['y_nm'] < 4800):
            continue
        ex = (stack[f_].astype(np.float64) - C.BASELINE_DN) / C.DN_PER_PHOTON - base
        d = np.hypot(xx + 0.5 - r['x_nm'] / C.PIXEL_NM, yy + 0.5 - r['y_nm'] / C.PIXEL_NM)
        b = np.clip(np.searchsorted(HALO_EDGES, d.ravel(), side='right') - 1, 0, nb - 1)
        S = np.bincount(b, weights=ex.ravel(), minlength=nb)
        N = np.bincount(b, minlength=nb)
        for c, sel in classes.items():
            if sel(np.array([r['z_nm']]))[0]:
                a = acc[c]
                g_ = k % GROUPS
                a['S'][g_] += S
                a['N'][g_] += N
                a['P'][g_] += r['photons']
                a['n'][g_] += 1
                a['tot'][g_] += ex.sum()
        k += 1
    area = np.pi * (HALO_EDGES[1:] ** 2 - HALO_EDGES[:-1] ** 2)
    out = {}
    for c, a in acc.items():
        n = int(a['n'].sum())
        if n < 5:
            continue
        S, N, P = a['S'].sum(axis=0), a['N'].sum(axis=0), a['P'].sum()
        prof = S / np.maximum(N, 1)                                              # mean excess photons per pixel in each annulus
        mean_p = P / n

        def cum(prof_, mean_p_):
            return np.array([np.sum((prof_ * area)[HALO_EDGES[1:] <= R]) / mean_p_ for R in HALO_RADII])
        full = cum(prof, mean_p)
        gcum = []
        for g_ in range(GROUPS):
            keep = np.arange(GROUPS) != g_                                       # jackknife: leave one group of emitters out
            pg = a['S'][keep].sum(axis=0) / np.maximum(a['N'][keep].sum(axis=0), 1)
            gcum.append(cum(pg, a['P'][keep].sum() / a['n'][keep].sum()))
        gcum = np.array(gcum)
        se = np.sqrt((GROUPS - 1) / GROUPS * ((gcum - gcum.mean(axis=0)) ** 2).sum(axis=0))
        out[c] = dict(n=n, mean_photons=mean_p, prof=prof, prof_per_1000=prof / mean_p * 1000.0, npix=N, cum=full, cum_se=se, image_total=float(a['tot'].sum() / P))
    return out


def illumination_profile(gt_all, summary_path):
    """The illumination of the simulation as the summary file lists it (centre, distance at half maximum) and the mean photons per activation by distance of the
    emitter to that centre: (centre (x, y) nm, half-maximum distance nm, rows (lo, hi, n, mean photons, median photons)) or None when the file cannot be read."""
    try:
        cen, half = [], None
        for line in open(summary_path):
            key, _, val = line.partition(',')
            key = key.strip()
            if key.startswith('Center of power illumination'):
                cen.append(float(val))
            elif key.startswith('Distance at half maximum illumination'):
                half = float(val)
        if len(cen) != 2 or half is None:
            return None
    except (OSError, ValueError):
        return None
    d = np.hypot(gt_all['x_nm'] - cen[0], gt_all['y_nm'] - cen[1])
    rows = []
    for lo, hi in ((0, 1000), (1000, 2000), (2000, 3000), (3000, 4000), (4000, 1e9)):
        m = (d >= lo) & (d < hi)
        if m.sum() >= 20:
            rows.append((lo, hi, int(m.sum()), float(gt_all['photons'][m].mean()), float(np.median(gt_all['photons'][m]))))
    return (cen[0], cen[1]), half, rows


# ------------------------------------------------------------------ offset estimators
def stacked_samples(frames, gt, rows, base, half=7):
    """Pixel samples around the emitters ``rows`` (indices into gt): for each emitter the (2 half + 1)^2 pixels around the pixel that contains it, as arrays
    (emitter number, dx, dy in px of the pixel centre minus the ground-truth position, excess photons per ground-truth photon, with the background ``base``)."""
    E, DX, DY, W = [], [], [], []
    for n, i in enumerate(rows):
        gx, gy = gt['x_nm'][i] / C.PIXEL_NM, gt['y_nm'][i] / C.PIXEL_NM
        ix, iy = int(gx), int(gy)
        p = (frames[int(gt['frame'][i]), iy - half:iy + half + 1, ix - half:ix + half + 1].astype(np.float64) - C.BASELINE_DN) / C.DN_PER_PHOTON - base
        jj, ii = np.meshgrid(np.arange(ix - half, ix + half + 1), np.arange(iy - half, iy + half + 1))
        E.append(np.full(p.size, n))
        DX.append((jj + 0.5 - gx).ravel())
        DY.append((ii + 0.5 - gy).ravel())
        W.append((p / gt['photons'][i]).ravel())
    return np.concatenate(E), np.concatenate(DX), np.concatenate(DY), np.concatenate(W)


def per_emitter_centroids(frames, gt, rows, base, radius=4.5, iters=4, half=7):
    """Centroid of the excess within ``radius`` px of the current estimate, started at the ground-truth position and iterated: (dx, dy) in nm, estimate minus truth."""
    out = []
    for i in rows:
        gx, gy = gt['x_nm'][i] / C.PIXEL_NM, gt['y_nm'][i] / C.PIXEL_NM
        cx, cy = gx, gy
        img = (frames[int(gt['frame'][i])].astype(np.float64) - C.BASELINE_DN) / C.DN_PER_PHOTON - base
        for _ in range(iters):
            ix, iy = int(np.clip(round(cx - 0.5), half, 63 - half)), int(np.clip(round(cy - 0.5), half, 63 - half))          # window kept inside the image
            p = img[iy - half:iy + half + 1, ix - half:ix + half + 1]
            jj, ii = np.meshgrid(np.arange(ix - half, ix + half + 1), np.arange(iy - half, iy + half + 1))
            w = np.where(np.hypot(jj + 0.5 - cx, ii + 0.5 - cy) < radius, p, 0.0)
            cx, cy = (w * (jj + 0.5)).sum() / w.sum(), (w * (ii + 0.5)).sum() / w.sum()
        out.append(((cx - gx) * C.PIXEL_NM, (cy - gy) * C.PIXEL_NM))
    return np.array(out)


def gaussian_profile_fit(dx, dy, w, radius=5.0):
    """Symmetric Gaussian + constant fitted to the pooled pixel samples within ``radius`` of the ground-truth position: (x0, y0 in nm, sigma in px)."""
    m = np.hypot(dx, dy) < radius
    dx, dy, w = dx[m], dy[m], w[m]
    r = least_squares(lambda t: t[4] + t[0] * np.exp(-0.5 * ((dx - t[1]) ** 2 + (dy - t[2]) ** 2) / t[3] ** 2) - w, [float(w.max()), 0.0, 0.0, 1.4, 0.0], x_scale='jac')
    return r.x[1] * C.PIXEL_NM, r.x[2] * C.PIXEL_NM, float(abs(r.x[3]))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--set', default='MT0.N1.LD')
    ap.add_argument('--second-set', default='MT0.N2.LD', help='the set of the second noise point (a different noise class of the same camera); "" to skip')
    ap.add_argument('--out', default=os.path.join(REPO, 'results', 'smlm', 'calibration.md'))
    ap.add_argument('--dev', type=int, nargs=2, default=list(C.DEV_FRAMES), help='development frames [start, stop)')
    a = ap.parse_args(argv)
    d0, d1 = a.dev
    stack, gt_all = load_stack(a.set)
    frames = stack[d0:d1]
    gt = gt_all[(gt_all['frame'] >= d0) & (gt_all['frame'] < d1)]
    nF, nAll = len(frames), len(stack)
    pmin_set = C.photon_threshold(gt_all)
    # the camera constants from the quiet frames (needed by the offset estimators of (a) and by (b), (c)): frames outside the test frames, and the development frames alone
    cal = calibration_frames(nAll)
    q_all = quiet_frames(gt_all, nAll)
    q_cal = cal[q_all[cal]]
    qs_cal = quiet_stats(stack, q_cal)
    q_dev = np.arange(d0, d1)[q_all[d0:d1]]
    qs_dev = quiet_stats(stack, q_dev)
    second = None
    if a.second_set:
        stack2, gt2 = load_stack(a.second_set)
        n2 = len(stack2)
        cal2 = calibration_frames(n2)
        q2 = cal2[quiet_frames(gt2, n2)[cal2]]
        second = dict(stack=stack2, gt=gt2, n=n2, cal=cal2, q=q2, qs=quiet_stats(stack2, q2))
    base = qs_cal['mean']
    L = []
    P = L.append
    P('# Calibration of %s from its development frames %d-%d\n' % (a.set, d0, d1 - 1))
    P('Produced by `experiments/exp_smlm_calibration.py`. The PSF, offset and emitter statistics use the first %d frames (the development frames of the protocol); the camera constants (b, c) '
      'and the halo profile use the quiet / single-emitter frames of the stack outside the test frames %d-%d, with the development-only numbers printed next to them; test frames are never used. '
      'Frame indices are 0-based stack indices, positions are in nm with the corner convention (see `experiments/smlm_common.py`). Camera, from the simulation parameters published on the '
      'dataset pages (%s, read 2026-11-07): baseline 100 DN, EM gain 300 (gamma-distributed), 45 e-/ADU, QE 0.90, read-out noise 74.4 e- (Gaussian), spurious (clock-induced) noise Poisson 0.002, '
      'wavelength 660 nm, NA 1.49, pixel 100 nm, "autofluorescence" Poisson at 100 nm pixel size; total gain 6.0 DN per photon (photons below are (DN - 100) / 6, i.e. incident photons, QE included).\n'
      % (nF, C.TEST_FRAMES[0], C.TEST_FRAMES[1] - 1, SIM_PAGE % a.set))
    try:
        ck = json.load(open(os.path.join(REPO, 'data', 'smlm', 'checksums.json')))
        sets = [a.set] + ([a.second_set] if a.second_set else [])
        P('Data: %s.\n' % ', '.join('`%s` (sha256 `%s...`), `activations.csv` (sha256 `%s...`)' % (os.path.basename(C.data_path(s, 'stack')),
                                                                                                    ck['%s/%s' % (s, os.path.basename(C.data_path(s, 'stack')))]['sha256'][:16],
                                                                                                    ck['%s/activations.csv' % s]['sha256'][:16]) for s in sets))
    except Exception:
        pass

    # ---------------- (0) conventions
    P('## 0. Conventions of the ground truth, checked against the images\n')
    bright = np.flatnonzero((gt['photons'] > 4000) & (gt['frame'] >= 3) & (gt['frame'] < nF - 3) & (gt['x_nm'] > 600) & (gt['x_nm'] < 5800) &
                            (gt['y_nm'] > 600) & (gt['y_nm'] < 5800))
    P('Mean excess of the 3 x 3 pixels around the ground-truth position over the frame median (DN), for %d bright (> 4000 photons) emitters, as a function of '
      'the offset between the stack index and (file frame - 1) and of the axis mapping. The correct mapping stands out sharply.\n' % len(bright))
    rows = []
    maps = [('x -> column, y -> row (used)', lambda x, y: (x, y)), ('x <-> y swapped', lambda x, y: (y, x)),
            ('x flipped', lambda x, y: (6400 - x, y)), ('y flipped', lambda x, y: (x, 6400 - y))]
    for name, fn in maps:
        vals = []
        for dlt in (-2, -1, 0, 1, 2):
            v = []
            for i in bright:
                k = int(gt['frame'][i]) + dlt
                if not 0 <= k < nF:
                    continue
                xx_, yy_ = fn(gt['x_nm'][i], gt['y_nm'][i])
                ix, iy = int(xx_ // C.PIXEL_NM), int(yy_ // C.PIXEL_NM)
                img = frames[k].astype(np.float64)
                v.append(img[max(iy - 1, 0):iy + 2, max(ix - 1, 0):ix + 2].mean() - np.median(img))
            vals.append(np.mean(v))
        rows.append([name] + [f(v, 1) for v in vals])
    P(md_table(['mapping'] + ['offset %+d' % d for d in (-2, -1, 0, 1, 2)], rows))
    P('\nOffset 0 is the stack index equal to (file frame - 1), i.e. the ground truth counts frames from 1. The pixel-centre convention is decided by the model-free centroid of the '
      'stacked PSF in section (a) (offset near zero under the corner convention, 50 nm in each axis under the other).\n')

    # ---------------- (a) PSF sigma
    iso = isolated_emitters(gt, nF, MIN_PHOTONS)
    iso_all = isolated_emitters(gt, nF, 0.0)
    P('## (a) Empirical PSF: 2D Gaussian + constant on isolated emitters\n')
    P('Selection: ground-truth emitters of frames %d-%d with no other ground-truth emitter of the same frame within %.0f px, the full %d x %d px fit window inside the '
      'image, and at least %.0f photons: **%d emitters** (of %d ground-truth activations in these frames, %d isolated and interior at any brightness). '
      'Least-squares fit in DN with amplitude, centre, width and constant free, **fit window %d x %d px unless another window is named**; "point" samples the Gaussian at the pixel centres (as Picasso does), "integrated" '
      'integrates it over each pixel. Expected from the optics (0.21 lambda / NA, 660 nm, NA 1.49): sigma = %.2f px = %.0f nm.\n'
      % (d0, d1 - 1, ISO_PX, 2 * HALF + 1, 2 * HALF + 1, MIN_PHOTONS, len(iso), len(gt), len(iso_all), 2 * HALF + 1, 2 * HALF + 1, EXPECTED_SIGMA_PX, EXPECTED_SIGMA_PX * C.PIXEL_NM))
    res = {k: [] for k in ('point', 'int', 'ellip')}
    keep_idx = []
    resid = []
    for i in iso:
        ix, iy = int(gt['x_nm'][i] // C.PIXEL_NM), int(gt['y_nm'][i] // C.PIXEL_NM)
        patch = frames[int(gt['frame'][i]), iy - HALF:iy + HALF + 1, ix - HALF:ix + HALF + 1].astype(np.float64)
        out = {k: fit_patch(patch, k) for k in res}
        if any(v is None for v in out.values()):
            continue
        keep_idx.append(i)
        for k, v in out.items():
            res[k].append(np.concatenate([v[0], [ix, iy]]))
        resid.append((patch, out['int'][1], out['int'][2]))
    keep_idx = np.array(keep_idx)
    R = {k: np.array(v) for k, v in res.items()}
    n_fit = len(keep_idx)
    zq = gt['z_nm'][keep_idx]
    ph = gt['photons'][keep_idx]
    pt, it, el = R['point'], R['int'], R['ellip']
    s_pt, s_int = pt[:, 3], it[:, 3]
    s_el = np.sqrt(el[:, 3] * el[:, 4])
    ell = np.abs(el[:, 3] - el[:, 4]) / (0.5 * (el[:, 3] + el[:, 4]))
    win = '%d x %d window' % (2 * HALF + 1, 2 * HALF + 1)
    rows = []
    for name, v in (('sigma, point-sampled circular (px, %s)' % win, s_pt), ('sigma, pixel-integrated circular (px, %s)' % win, s_int),
                    ('sigma, geometric mean of the elliptical fit (px, %s)' % win, s_el), ('ellipticity |sx - sy| / mean (%s)' % win, ell)):
        s = _spread(v)
        rows.append([name, f(s['median'], 3), f(s['mean'], 3), f(s['std'], 3), f(s['iqr'], 3), f(s['mad_sigma'], 3), '%s .. %s' % (f(s['p5'], 3), f(s['p95'], 3))])
    P(md_table(['quantity', 'median', 'mean', 'std', 'IQR', 'robust sd (1.48 MAD)', '5% .. 95%'], rows))
    s_pt_med, s_int_med = float(np.median(s_pt)), float(np.median(s_int))
    P('\n**Median sigma, all |z| (%s): %.3f px (point-sampled, %.0f nm), %.3f px (pixel-integrated, %.0f nm)**, against %.2f px expected from the optics (n = %d fits). '
      'The width is not a constant: the PSF broadens with defocus, and the emitters sit at z = %.0f ... %.0f nm, so the spread is partly physical (defocus), partly fit noise. '
      'The defocus widths depend on the fit window (table below).\n' % (win, s_pt_med, s_pt_med * C.PIXEL_NM, s_int_med, s_int_med * C.PIXEL_NM,
                                                                                                    EXPECTED_SIGMA_PX, n_fit, gt['z_nm'].min(), gt['z_nm'].max()))
    m0 = np.abs(zq) < 100
    P('**In focus (|z| < 100 nm, %d emitters, %s): median sigma %.3f px (pixel-integrated), %.3f px (point-sampled)**, i.e. %.0f nm, %.0f %% wider than the %.2f px of the ideal optics.\n'
      % (m0.sum(), win, np.median(s_int[m0]), np.median(s_pt[m0]), np.median(s_int[m0]) * C.PIXEL_NM, 100 * (np.median(s_int[m0]) / EXPECTED_SIGMA_PX - 1), EXPECTED_SIGMA_PX))
    P('Pixel-integrated minus point-sampled: median %.3f px (a point-sampled model must be wider by about sqrt(sigma^2 + 1/12) - sigma = %.3f px).\n'
      % (float(np.median(s_pt - s_int)), float(np.sqrt(s_int_med ** 2 + 1 / 12) - s_int_med)))
    rows = []
    for lo, hi in Z_BINS:
        m = (np.abs(zq) >= lo) & (np.abs(zq) < hi)
        if m.sum() >= 5:
            rows.append(['%d-%d nm' % (lo, hi), int(m.sum()), f(np.median(s_int[m]), 3), f(np.median(s_pt[m]), 3), f(np.percentile(s_int[m], 25), 3) + ' .. ' + f(np.percentile(s_int[m], 75), 3)])
    P('Width by axial position |z| of the emitter (%s):\n' % win)
    P(md_table(['|z|', 'n', 'median sigma_int (px)', 'median sigma_point (px)', 'IQR sigma_int'], rows))
    # window dependence
    iso_w = isolated_emitters(gt, nF, MIN_PHOTONS, half=max(WINDOW_HALVES), iso_px=ISO_WINDOW_PX)
    sig_w = {h: [] for h in WINDOW_HALVES}
    for i in iso_w:
        ix, iy = int(gt['x_nm'][i] // C.PIXEL_NM), int(gt['y_nm'][i] // C.PIXEL_NM)
        for h in WINDOW_HALVES:
            patch = frames[int(gt['frame'][i]), iy - h:iy + h + 1, ix - h:ix + h + 1].astype(np.float64)
            fit = fit_patch(patch, 'int', h)
            sig_w[h].append(fit[0][3] if fit is not None else np.nan)
    sig_w = {h: np.array(v) for h, v in sig_w.items()}
    zw = gt['z_nm'][iso_w]
    P('\nWindow dependence of the widths: the SAME %d emitters (isolated within %.0f px, at least %.0f photons, a 19 x 19 window inside the image) fitted with pixel-integrated Gaussians in '
      'windows of different size, median sigma in px (a window that is too small for a broad defocused PSF cuts its wings, a large one admits the halo and the neighbours):\n' % (len(iso_w), ISO_WINDOW_PX, MIN_PHOTONS))
    rows = []
    for lab, m in [('all |z|', np.ones(len(zw), bool))] + [('|z| %d-%d nm' % (lo, hi), (np.abs(zw) >= lo) & (np.abs(zw) < hi)) for lo, hi in Z_BINS]:
        if m.sum() >= 5:
            rows.append([lab, int(m.sum())] + [f(float(np.nanmedian(sig_w[h][m])), 3) for h in WINDOW_HALVES])
    P(md_table(['emitters', 'n'] + ['%d x %d' % (2 * h + 1, 2 * h + 1) for h in WINDOW_HALVES], rows))
    def window_spread(mask):
        meds = [float(np.nanmedian(sig_w[h][mask])) for h in WINDOW_HALVES]
        return 100.0 * (max(meds) - min(meds)) / min(meds)
    m_in, m_far = np.abs(zw) < 100, np.abs(zw) >= 350
    P('\nAcross the windows the median width of the in-focus emitters (|z| < 100 nm, n = %d) varies by %.0f %%, that of the strongly defocused ones (|z| >= 350 nm, n = %d) by %.0f %%: a defocus width is only meaningful '
      'together with the window it was fitted in, an in-focus width hardly depends on it.\n' % (m_in.sum(), window_spread(m_in) if m_in.sum() >= 5 else float('nan'), m_far.sum(), window_spread(m_far) if m_far.sum() >= 5 else float('nan')))
    rows = []
    for lo, hi in ((MIN_PHOTONS, 3500), (3500, 5000), (5000, 1e9)):
        m = (ph >= lo) & (ph < hi)
        if m.sum() >= 5:
            rows.append(['%.0f-%s photons' % (lo, '%.0f' % hi if hi < 1e8 else 'inf'), int(m.sum()), f(np.median(s_int[m]), 3), f(np.std(s_int[m], ddof=1), 3)])
    P('\nWidth by brightness (%s; to see that the selection does not drive the value):\n' % win)
    P(md_table(['photons', 'n', 'median sigma_int (px)', 'std'], rows))
    # centre offsets of the integrated fit under the two pixel conventions
    gx, gy = gt['x_nm'][keep_idx] / C.PIXEL_NM, gt['y_nm'][keep_idx] / C.PIXEL_NM
    fx, fy = it[:, 5] + it[:, 1], it[:, 6] + it[:, 2]              # fitted centre as a pixel INDEX (pixel centre = integer); columns: N, x0, y0, s, bg, ix, iy
    off = {}
    for conv, sh in (('corner (centre of pixel j at (j + 0.5) * 100 nm; used)', 0.5), ('centre (centre of pixel j at j * 100 nm)', 0.0)):
        off[conv] = ((fx - (gx - sh)) * C.PIXEL_NM, (fy - (gy - sh)) * C.PIXEL_NM)
    rows = []
    for conv, (ox_, oy_) in off.items():
        rows.append([conv, f(np.median(ox_), 1), f(np.median(oy_), 1), f(ox_.std(ddof=1), 1), f(oy_.std(ddof=1), 1)])
    P('\nFitted centre minus ground truth (nm, pixel-integrated fit, %s) under the two pixel conventions: the right one has a small median offset, the wrong one is displaced by a further 50 nm in each axis.\n' % win)
    P(md_table(['convention', 'median dx (nm)', 'median dy (nm)', 'std dx (nm)', 'std dy (nm)'], rows))
    ox, oy = off['corner (centre of pixel j at (j + 0.5) * 100 nm; used)']

    # ---------------- (a2) the offset is estimator- and z-dependent
    P('\n### Offset of a symmetric centre against the simulated emitter: it depends on the estimator and on |z|\n')
    P('The simulated PSF is not symmetric, so the centre that a symmetric estimator finds is displaced against the emitter position, and by how much depends on WHICH estimator and on |z|. '
      'It is not a constant of the data set: it must not be read as a convention error, and a correction estimated for one estimator (or one z range) does not carry over to another. Estimator by estimator '
      '(offset = estimate minus truth, nm, corner convention; isolated emitters of the development frames, at least %.0f photons, a 15 x 15 window inside the image; the background is the quiet-frame mean of section (b)):\n' % MIN_PHOTONS)
    iso7 = isolated_emitters(gt, nF, MIN_PHOTONS, half=7)
    E_, DX, DY, W = stacked_samples(frames, gt, iso7, base, half=7)
    cen = per_emitter_centroids(frames, gt, iso7, base)
    z7 = gt['z_nm'][iso7]
    # the per-emitter integrated fits (11 x 11), matched to the 15 x 15 sample
    pos_in_fit = {int(i): k for k, i in enumerate(keep_idx)}
    rows = []
    for lab, mz in [('all |z|', np.ones(len(iso7), bool))] + [('|z| %d-%d nm' % (lo, hi), (np.abs(z7) >= lo) & (np.abs(z7) < hi)) for lo, hi in Z_BINS]:
        if mz.sum() < 8:
            continue
        sel_em = np.flatnonzero(mz)
        mk = np.isin(E_, sel_em)
        r_ = np.hypot(DX[mk], DY[mk])
        cells = [lab, int(mz.sum())]
        cen_r = []
        for R_ in (3.0, 5.0, 6.0):
            w_ = W[mk][r_ < R_]
            cen_r.append((100 * (w_ * DX[mk][r_ < R_]).sum() / w_.sum(), 100 * (w_ * DY[mk][r_ < R_]).sum() / w_.sum()))
            cells.append('(%+.1f, %+.1f)' % cen_r[-1])
        gx0, gy0, gs = gaussian_profile_fit(DX[mk], DY[mk], W[mk])
        cells.append('(%+.1f, %+.1f), sigma %.3f' % (gx0, gy0, gs))
        cells.append('(%+.1f, %+.1f)' % (np.median(cen[mz, 0]), np.median(cen[mz, 1])))
        fi = [pos_in_fit[int(i)] for i in iso7[mz] if int(i) in pos_in_fit]
        cells.append('(%+.1f, %+.1f), n %d' % (np.median(ox[fi]), np.median(oy[fi]), len(fi)) if fi else 'n/a')
        rows.append(cells)
        if lab == 'all |z|':
            row_all = cells
            cen_all = cen_r
    P(md_table(['emitters', 'n', 'centroid of the stacked PSF, r < 3 px', 'r < 5 px', 'r < 6 px (whole PSF)', 'Gaussian fit to the stacked PSF, r < 5 px (nm; width in px)', 'per-emitter centroid (median, 4.5 px)',
                'per-emitter integrated Gaussian fit (11 x 11, median)'], rows))
    P('\n"Stacked PSF": the pixels of all emitters of the row at their position relative to the ground-truth position, each emitter normalised by its ground-truth photons; its centroid within a radius of the '
      'ground-truth position is the offset of the mean profile, model free. Pooled over |z| (row "all |z|") it is %s for r < 3 px, %s for r < 5 px and %s for the whole PSF (r < 6 px): the offset of a centre '
      'estimate depends on how much of the PSF it uses, and the whole-PSF centroid is within %.0f nm of the ground-truth position in each axis, which is the evidence for the corner convention (the other convention would be 50 nm '
      'off in each axis). The offsets of the Gaussian fits are those of the PSF core against the centre of mass; per-emitter estimators scatter around them; all of them change with |z| (rows). '
      'Picasso\'s own least-squares fits show their offsets in `picasso_baseline.md`.\n' % (row_all[2], row_all[3], row_all[4], max(abs(cen_all[2][0]), abs(cen_all[2][1]))))
    # per-emitter integrated fit by z (the earlier table)
    P('Residual offset of the pixel-integrated fit (%s, as above) under the corner convention by axial position:\n' % win)
    rows = []
    for lo, hi in ((-800, -400), (-400, -150), (-150, 150), (150, 400), (400, 800)):
        m = (zq >= lo) & (zq < hi)
        if m.sum() >= 5:
            rows.append(['%d .. %d nm' % (lo, hi), int(m.sum()), f(np.median(ox[m]), 1), f(np.median(oy[m]), 1)])
    P(md_table(['z', 'n', 'median dx (nm)', 'median dy (nm)'], rows))
    P('\nA method that fits a symmetric model therefore carries an offset of about (%.0f, %.0f) nm (this estimator, all |z|, median) into its RMSE, and the offset grows with |z|; `score` reports `bias_x_nm`, '
      '`bias_y_nm` and `rmse_centered_nm`, and the baseline ladder reports every row also with a constant offset estimated on the development frames subtracted (`picasso_baseline.md`).\n' % (np.median(ox), np.median(oy)))

    # ---------------- (b) background on the quiet frames
    P('## (b) Background and noise on the quiet frames, and the PSF halo\n')
    P('**Quiet frames**: frames with no ground-truth emitter in the frame or in either neighbouring frame (so that no emitter, and no signal of one that has just switched on or off, is in the frame). '
      'In %s there are %d such frames among the %d frames outside the test frames (%d among the %d development frames)%s. All values in photons, x = (DN - 100) / 6; standard errors from %d contiguous groups of frames.\n'
      % (a.set, len(q_cal), len(cal), len(q_dev), nF, '; in %s %d among the %d frames outside the test frames' % (a.second_set, len(second['q']), len(second['cal'])) if second else '', GROUPS))
    hdr = ['statistic', '%s, quiet frames of the development frames' % a.set, '%s, quiet frames outside the test frames' % a.set] + (['%s, quiet frames outside the test frames' % a.second_set] if second else [])
    cols = [qs_dev, qs_cal] + ([second['qs']] if second else [])
    rows = [['quiet frames'] + [c['m'] for c in cols],
            ['mean (photons)'] + ['%.3f +- %.3f' % (c['mean'], c['se_mean']) for c in cols],
            ['mean (DN)'] + ['%.2f' % (c['mean'] * C.DN_PER_PHOTON + C.BASELINE_DN) for c in cols],
            ['median (photons)'] + ['%.2f' % c['median'] for c in cols],
            ['std of the pixel values (photons)'] + ['%.3f' % c['sd'] for c in cols],
            ['skewness'] + ['%.3f' % c['skew'] for c in cols],
            ['temporal variance, mean over pixels (photon^2)'] + ['%.2f +- %.2f' % (c['var'], c['se_var']) for c in cols],
            ['variance / mean (photons)'] + ['%.3f' % c['var_over_mean'] for c in cols],
            ['sd of the frame means (photons)'] + ['%.3f' % c['frame_sd'] for c in cols],
            ['... expected from independent pixel noise alone'] + ['%.3f' % c['expect_sd'] for c in cols],
            ['linear trend of the frame mean (photons per 1000 frames)'] + ['%.3f' % c['trend'] for c in cols],
            ['sd of the per-pixel mean map (photons)'] + ['%.3f' % c['spatial_sd'] for c in cols],
            ['... expected from sampling alone'] + ['%.3f' % c['spatial_expect'] for c in cols]]
    P(md_table(hdr, rows))
    P('\nThe background of a quiet frame is **constant in space and in time**: radial means of the quiet-frame mean image (%s, frames outside the test frames), photons: %s; the sd of the frame means equals '
      'that of independent pixel noise (%.3f against %.3f), so there is no common-mode fluctuation from frame to frame, no spatial profile (the sd of the per-pixel mean map equals its sampling error, %.3f against %.3f) '
      'and no trend. The background of %s is %.2f photons (DN %.1f)%s.\n'
      % (a.set, '; '.join('r %d-%d px %.3f +- %.3f' % tuple(r_) for r_ in qs_cal['radial']), qs_cal['frame_sd'], qs_cal['expect_sd'], qs_cal['spatial_sd'], qs_cal['spatial_expect'],
         a.set, qs_cal['mean'], qs_cal['mean'] * C.DN_PER_PHOTON + C.BASELINE_DN, ', that of %s %.2f photons (DN %.1f)' % (a.second_set, second['qs']['mean'], second['qs']['mean'] * C.DN_PER_PHOTON + C.BASELINE_DN) if second else ''))
    hl = halo_leak(stack, gt_all, cal, qs_cal['mean'])
    P('**What the frames with emitters show instead is the PSF halo of the emitters.** In the %d frames outside the test frames that have at least one emitter, the pixels farther than %.0f px from every ground-truth emitter '
      'of the frame ("emitter-free" pixels) have a mean of %.2f photons, %.2f above the quiet-frame background, and the sd of the per-frame mean of these pixels is %.3f (quiet frames: %.3f). Regressing that frame mean on the '
      'total ground-truth photons of the frame gives a slope of %.2e photons per pixel per emitted photon (correlation %.3f); the residual sd after removing it is %.3f (pixel noise alone: %.3f). '
      'Over the %.0f emitter-free pixels per frame this is %.2f photons beyond %.0f px per emitted photon: signal far from the emitter, not illumination (the background of quiet frames is flat) and not a common-mode '
      'fluctuation (it is a deterministic function of the photons emitted in the frame).\n'
      % (hl['n'], ISO_PX, hl['mean'], hl['excess'], hl['sd'], qs_cal['frame_sd'], hl['slope'], hl['corr'], hl['res_sd'], math.sqrt(qs_cal['var'] / hl['nfree']), hl['nfree'], hl['leaked'], ISO_PX))
    # radial profile
    P('### Radial excess profile around isolated single-emitter frames\n')
    pm1 = pmin_set
    hp = {a.set: halo_profile(stack, gt_all, cal, qs_cal['mean'], pm1)}
    if second:
        hp[a.second_set] = halo_profile(second['stack'], second['gt'], second['cal'], second['qs']['mean'], C.photon_threshold(second['gt']))
    P('Frames outside the test frames with exactly one ground-truth emitter, brighter than the lower quartile of the whole set\'s activations (%s: %.0f photons%s), at least 16 px from every border; the image minus the quiet-frame '
      'background of its set is binned by the distance of the pixel centre to the ground-truth position. Cumulative fraction of the emitter\'s ground-truth photons within a radius = sum of the mean excess per pixel of the '
      'annuli times their full areas (isotropy assumed beyond the frame edge) divided by the mean ground-truth photons; +- from a leave-one-group-out jackknife over %d groups of emitters.\n'
      % (a.set, pm1, '; %s: %.0f' % (a.second_set, C.photon_threshold(second['gt'])) if second else '', GROUPS))
    rows = []
    for s_, hs in hp.items():
        for cname, h in hs.items():
            rows.append([s_, cname, h['n'], '%.0f' % h['mean_photons']] + ['%.3f +- %.3f' % (h['cum'][k], h['cum_se'][k]) for k in range(len(HALO_RADII))] + ['%.3f' % h['image_total']])
    P(md_table(['set', 'emitters', 'n', 'mean ground-truth photons'] + ['within %d px' % R_ for R_ in HALO_RADII] + ['whole image (64 x 64), summed directly'], rows))
    P('')
    rows = []
    sets_h = list(hp)
    for k in range(len(HALO_EDGES) - 1):
        if HALO_EDGES[k] >= 50:
            continue
        rows.append(['%g-%g' % (HALO_EDGES[k], HALO_EDGES[k + 1])] + sum([[('%.4f' % hp[s_]['all']['prof_per_1000'][k]), ('%.3f' % np.sum((hp[s_]['all']['prof'] * np.pi * (HALO_EDGES[1:] ** 2 - HALO_EDGES[:-1] ** 2))[:k + 1] / hp[s_]['all']['mean_photons']))]
                                                                         for s_ in sets_h if 'all' in hp[s_]], []))
    P('Mean excess per pixel per 1000 ground-truth photons in each annulus and the cumulative fraction (all emitters of the row "all" above):\n')
    P(md_table(['annulus (px)'] + sum([['%s excess / px / 1000 photons' % s_, '%s cumulative' % s_] for s_ in sets_h if 'all' in hp[s_]], []), rows))
    h1 = hp[a.set]['all']
    hin = hp[a.set].get('in focus (|z| < 100 nm)', h1)
    h2 = hp[a.second_set].get('all') if second else None
    if h2 is not None:
        same = max(abs(h1['cum'][k] - h2['cum'][k]) for k in range(len(HALO_RADII)))
        sentence = ('In %s the fractions are %s; they differ from those of %s by at most %.2f, so normalised by the photons the profile %s in both noise classes: it is a property of the simulated PSF, not of the noise'
                    % (a.second_set, ', '.join('%.2f at %d px' % (h2['cum'][k], HALO_RADII[k]) for k in range(len(HALO_RADII))), a.set, same, 'is the same' if same < 0.1 else 'is NOT the same'))
    else:
        sentence = 'No second set was given, so the noise dependence of the profile was not checked'
    edges = list(HALO_EDGES)

    def prof_at(lo):
        return h1['prof_per_1000'][edges.index(float(lo))]
    P('\nThe ground-truth "intensity" corresponds to the core of the PSF: the cumulative excess within 5 px is %.2f of it for the emitters in focus (%.2f over all emitters, whose PSFs are broader). But the PSF has a long halo: for the in-focus emitters the cumulative '
      'fraction keeps rising to %.2f at 25 px (%.0f %% of the ground-truth photons arrive beyond 5 px), and then the profile is cut off abruptly: the excess per pixel per 1000 ground-truth photons (all emitters) is %.3f at 20-22 px, %.3f at 22-23 px, %.3f at 24-25 px '
      'and %.3f (consistent with zero) at 26-28 px. %s. The cause of the halo was not investigated. Consequences: a "Gaussian + constant" patch model (Picasso\'s fit) cannot hold the halo, '
      'so its photon count is biased low (`picasso_baseline.md`), and any estimate of the background from the neighbourhood of an emitter (as the earlier version of this file did with an 8 px mask: 100.3 photons instead of %.2f) is raised by it.\n'
      % (hin['cum'][1], h1['cum'][1], hin['cum'][3], 100 * (hin['cum'][3] - hin['cum'][1]), prof_at(20), prof_at(22), prof_at(24), prof_at(26), sentence, qs_cal['mean']))

    # ---------------- (c) noise model
    P('## (c) Noise model: var = a x + b in photon units, x = (DN - 100) / 6\n')
    Q_E, EM, SIGR = C.QE, C.EM_GAIN, 74.4
    a_nom, b_nom = 2.0 / Q_E, (SIGR / (Q_E * EM)) ** 2 + (1.0 / 12.0) / C.DN_PER_PHOTON ** 2
    P('Nominal EMCCD model (gamma gain noise, excess noise factor^2 = 2, read noise %.1f e-, quantisation): a = 2 / QE = %.3f, b = (read noise / (QE EM))^2 + 1/(12 G^2) = %.4f photon^2.\n'
      % (SIGR, a_nom, b_nom))
    x1, v1 = qs_cal['mean'], qs_cal['var']
    P('**One operating point (quiet frames of %s).** Mean x = %.3f photons, temporal variance %.2f +- %.2f photon^2 (mean over the pixels of the per-pixel variance), variance / mean = %.3f, '
      'a_eff = (var - b) / x = %.3f with the nominal b, against a = %.3f expected.\n' % (a.set, x1, v1, qs_cal['se_var'], v1 / x1, (v1 - b_nom) / x1, a_nom))
    if second:
        qs2 = second['qs']
        x2, v2 = qs2['mean'], qs2['var']
        a2 = (v1 - v2) / (x1 - x2)
        b2 = v1 - a2 * x1
        ga = (qs_cal['gv'] - qs2['gv']) / (qs_cal['gm'] - qs2['gm'])                           # paired groups of frames: one estimate of a per group
        gb = qs_cal['gv'] - ga * qs_cal['gm']
        P('**Two operating points (the backgrounds of two noise classes of the same camera).** %s: x = %.3f photons, var = %.2f +- %.2f photon^2, var / mean = %.3f; %s: x = %.3f photons, var = %.3f +- %.3f photon^2, '
          'var / mean = %.3f. The line through the two points: **a = %.3f +- %.3f, b = %.2f +- %.2f photon^2** (standard errors from %d paired groups of frames), against the nominal a = %.3f, b = %.3f photon^2. '
          'The fitted b differs from the nominal b by %.2f photon^2 (%.1f standard errors of the group estimate; a systematic cause such as the spurious noise or an excess-noise factor that is not exactly 2 is not excluded). '
          'It is %.2f %% of the variance at the background of %s (%.0f photon^2), so it is irrelevant at that level and matters only for very dark pixels (var = a x + b at x = 0).\n'
          % (a.set, x1, v1, qs_cal['se_var'], v1 / x1, a.second_set, x2, v2, qs2['se_var'], v2 / x2, a2, ga.std(ddof=1) / math.sqrt(GROUPS), b2, gb.std(ddof=1) / math.sqrt(GROUPS), GROUPS, a_nom, b_nom,
             b2 - b_nom, abs(b2 - b_nom) / (gb.std(ddof=1) / math.sqrt(GROUPS)), 100 * abs(b2) / v1, a.set, v1))
    # emitter pixels: leverage-corrected residuals of the integrated fits, binned by the fitted mean
    mod = np.concatenate([m_.ravel() for _, m_, _ in resid])
    rsd = np.concatenate([(p_ - m_).ravel() for p_, m_, _ in resid])
    lev = np.concatenate([h_.ravel() for _, _, h_ in resid])
    xmod = (mod - C.BASELINE_DN) / C.DN_PER_PHOTON
    vres = rsd ** 2 / C.DN_PER_PHOTON ** 2 / (1.0 - lev)
    edges = np.array([50, 80, 100, 120, 150, 200, 300, 450, 700, 1100, 1800, 1e9])
    bx, bv, bn = [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (xmod >= lo) & (xmod < hi)
        if m.sum() > 200:
            bx.append(xmod[m].mean())
            bv.append(vres[m].mean())
            bn.append(int(m.sum()))
    bx, bv, bn = np.array(bx), np.array(bv), np.array(bn)
    a3, b3, sa3, sb3 = ols(bx, bv, w=bn / (a_nom * bx + b_nom) ** 2)
    a4 = float(np.sum(bn * bx * (bv - b_nom) / (a_nom * bx + b_nom) ** 2) / np.sum(bn * bx ** 2 / (a_nom * bx + b_nom) ** 2))
    P('**Emitter pixels (a rough check).** The %d pixels of the %d isolated-emitter patches (%s), residual of the pixel-integrated Gaussian fit divided by (1 - leverage), '
      'squared, averaged in bins of the fitted mean x. Weighted regression over the %d bins (x = %.0f .. %.0f): **a = %.2f +- %.2f, b = %.0f +- %.0f photon^2**; with b fixed at '
      'the nominal value a = %.2f. The bins are biased by the PSF misfit (upwards: the simulated PSF is not Gaussian, it has the halo of section (b)) and by the imperfect leverage correction of '
      'the unweighted fit (downwards at the brightest pixels), so treat them as a consistency check (a of the order of 2-3), not as a measurement.\n'
      % (len(mod), len(resid), win, len(bx), bx.min(), bx.max(), a3, sa3, b3, sb3, a4))
    P(md_table(['mean x (photons)', 'pixels', 'variance (photon^2)', 'nominal a x + b'], [[f(x, 1), n, f(v, 1), f(a_nom * x + b_nom, 1)] for x, v, n in zip(bx, bv, bn)]))
    P('\nThe regression of per-pixel variance on per-pixel mean across emitter-free pixels is not informative (the background is flat, so the means span only a few tenths of a photon and the sample mean and variance '
      'of a skewed distribution are correlated), which is why the two backgrounds, not a regression, fix a and b.\n')
    P('**Result: a = %.2f (quiet frames; nominal %.2f)%s, b = %s photon^2 (nominal %.3f).**\n\n'
      % ((v1 - b_nom) / x1, a_nom, '; two-point a = %.3f' % a2 if second else '', '%.2f +- %.2f' % (b2, gb.std(ddof=1) / math.sqrt(GROUPS)) if second else 'not measurable', b_nom))

    # ---------------- (d) emitters and photons
    P('## (d) Ground-truth emitters per frame and photons per emitter (frames %d-%d)\n' % (d0, d1 - 1))
    cnt = np.bincount(gt['frame'] - d0, minlength=nF)
    mu, var = cnt.mean(), cnt.var(ddof=1)
    kk = np.arange(0, cnt.max() + 1)
    P('%d ground-truth activations in %d frames: mean %.3f per frame, variance %.3f (variance / mean = %.2f; 1 for a Poisson count), %.1f %% of the frames are empty, maximum %d.\n'
      % (len(gt), nF, mu, var, var / mu, 100 * np.mean(cnt == 0), cnt.max()))
    from math import exp, factorial
    P(md_table(['emitters in frame'] + [str(k) for k in kk], [['frames'] + [str(int((cnt == k).sum())) for k in kk],
                                                                 ['Poisson, same mean'] + [f(nF * exp(-mu) * mu ** k / factorial(int(k)), 1) for k in kk]]))
    qq = np.percentile(gt['photons'], [1, 5, 25, 50, 75, 95, 99])
    P('\nPhotons per activation (the ground-truth `intensity`, photons emitted during the frame; partial-frame activations are included; development frames): mean %.0f, median %.0f, min %.0f, max %.0f. '
      'Whole set (%d activations): median %.0f, 25 %% quantile %.0f.\n' % (gt['photons'].mean(), np.median(gt['photons']), gt['photons'].min(), gt['photons'].max(), len(gt_all), np.median(gt_all['photons']), pmin_set))
    P(md_table(['1%', '5%', '25%', '50%', '75%', '95%', '99%'], [[f(v, 0) for v in qq]]))
    P('\nFraction of the development activations below 250 / 500 / 1000 / 2000 photons: %s. A quarter of the activations of the WHOLE set are dimmer than %.0f photons (the quartile used throughout the harness; '
      'the development frames alone give %.0f); the challenge discarded the dimmest quarter, after pairing, when it assessed software (rule from memory, `smlm_common.py`). A %.0f-photon emitter (the median of the whole set) of width %.2f px has a peak signal of about %.0f photons (%.0f DN) per '
      'pixel above the background of %.2f photons.\n'
      % (' / '.join('%.1f %%' % (100 * np.mean(gt['photons'] < t)) for t in (250, 500, 1000, 2000)), pmin_set, np.percentile(gt['photons'], 25), np.median(gt_all['photons']), s_int_med,
         np.median(gt_all['photons']) / (2 * np.pi * s_int_med ** 2), np.median(gt_all['photons']) / (2 * np.pi * s_int_med ** 2) * C.DN_PER_PHOTON, qs_cal['mean']))
    ill = illumination_profile(gt_all, C.data_path(a.set, 'summary'))
    if ill:
        P('**The illumination of the simulation acts on the emitters, not on the background.** The simulator\'s summary file (`stats_activation_sumary.csv`) lists an illumination profile (centre of power (%.1f, %.1f) nm, '
          'both coordinates are labelled "X" in the file; distance at half maximum %.0f nm). The brightness of the activations of the whole set falls with the distance of the emitter to that centre (mean photons per activation in the table), '
          'while the quiet frames of section (b) are flat:\n' % (ill[0][0], ill[0][1], ill[1]))
        P(md_table(['distance to the centre', 'activations', 'mean photons', 'median photons'], [['%d-%s nm' % (lo, '%d' % hi if hi < 1e8 else 'inf'), n, '%.0f' % mu_, '%.0f' % md_] for lo, hi, n, mu_, md_ in ill[2]]))
        P('')
    P('Axial positions: z from %.0f to %.0f nm (median |z| %.0f nm); %.0f %% of the activations lie within +-250 nm of the focal plane.\n'
      % (gt['z_nm'].min(), gt['z_nm'].max(), np.median(np.abs(gt['z_nm'])), 100 * np.mean(np.abs(gt['z_nm']) <= 250)))
    P('**Persistence of an emitter over consecutive frames.** A fluorophore stays on for several frames (runs of the same coordinates in consecutive frames of the whole ground truth), so neighbouring frames '
      'are not independent samples and frame statistics need a block bootstrap (`smlm_common.bootstrap_ci`, blocks of 10 frames):\n')
    rows = []
    for s_ in C.SETS:
        try:
            rl = C.run_lengths(gt_all if s_ == a.set else C.load_gt(C.data_path(s_, 'gt')))
        except Exception:
            continue
        rows.append([s_, rl['runs'], '%.3f' % rl['mean'], '%.1f %%' % (100 * rl['ge3'])] + [int(v) for v in rl['dist']])
    P(md_table(['set', 'runs', 'mean run length (frames)', 'runs of 3 or more frames'] + ['runs of %d' % k for k in range(1, 9)], rows))
    P('')
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(a.out, 'w') as fh:
        fh.write('\n'.join(L) + '\n')
    print('wrote', a.out)


if __name__ == '__main__':
    main()
