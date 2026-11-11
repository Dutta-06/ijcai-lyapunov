"""Usability check of the real known-geometry SMLM recordings fetched by ``data/smlm/fetch_real_gt.py`` (NO scoring of any localizer, nothing tuned).

What it prints (all of it is quoted in ``data/smlm/real_candidates.md`` and ``results/smlm/real_gt_plan.md``):
  * frame size, number of readable frames, dark level and range of the three recordings (from the files),
  * the metadata that is (and is not) stored in the files (Micro-Manager metadata block of the OME-TIFF, ThunderSTORM protocol of the Heil record),
  * Picasso 0.11.3 (``smlm_common.run_picasso``, its default pipeline R0: box 7, minimum net gradient 5000, ``gausslq``) on a few frames: the number of
    localizations per frame and, for the Nup96 frames, the fraction that lies within 150 nm of a localization of the same frame in the SMAP table of the
    authors (a sanity check that the frames are read with the right orientation and pixel size, NOT a score),
  * the NPC rings that a crude matched-filter finds in the published localization tables (centres, radii, localizations per ring), which show that the
    geometry the plan relies on is present in the data.
Needs numpy, scipy, tifffile and Picasso 0.11.3 (``pip install tifffile scipy picassosr==0.11.3``); run it with the interpreter of that environment:
    SMLM_PICASSO_PYTHON=$(which python) python experiments/real_gt_usability_check.py [--no-picasso]
Raw data are read from $SMLM_REAL_GT (default /home/user/smlm_real_data). Truncated prefixes make tifffile print 'invalid page offset' warnings: expected.
"""
import argparse
import json
import os
import sys
import warnings

import numpy as np
import tifffile
from scipy import ndimage as ndi, signal
from scipy.spatial import cKDTree

warnings.filterwarnings('ignore')
HERE = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, HERE)
ROOT = os.environ.get('SMLM_REAL_GT', '/home/user/smlm_real_data')
HEIL = os.path.join(ROOT, 'zenodo_5068525')
NUP = os.path.join(ROOT, 'zenodo_3540443')
SIM = os.path.join(ROOT, 'zenodo_7771932')
NUP_PREFIX = 'Cos7_Nup96_A647_COT_Au_%d_1_MMStack_Pos0.ome.tif.part0-400MB'


def read_pages(path, start, stop):
    with tifffile.TiffFile(path) as t:
        return t.asarray(key=range(start, stop)), len(t.pages)


def mm_meta(path):
    with tifffile.TiffFile(path) as t:
        m = t.pages[0].tags[51123].value
    keep = ('ROI', 'Exposure-ms', 'Camera', 'Camera-Name', 'Camera-ChipName', 'Camera-Gain', 'Camera-Offset', 'Camera-ReadoutRate', 'Camera-ActualInterval-ms',
            'PixelSizeUm', 'HamamatsuHam_DCAM-CameraName', 'Binning', 'Laser 642nm-Power Setpoint')
    out = {}
    for k in keep:
        v = m.get(k)
        out[k] = v.get('PropVal') if isinstance(v, dict) else v
    return out


def picasso(frames, params):
    import smlm_common as C
    info = {}
    dets = C.run_picasso(frames.astype(np.float32), params, info=info)
    n = sum(len(d) for d in dets)
    return dets, n, info


def rings(x, y, R_px, pix, size, r_in, r_out, search, pct):
    """Crude NPC finder on a localization table: ring matched filter on a histogram, then an algebraic circle fit (3 iterations) of the localizations
    in an annulus. Returns an (n, 4) array (cx, cy, radius, number of localizations in the annulus), nm."""
    H, _, _ = np.histogram2d(y, x, bins=[np.arange(0, size + pix, pix)] * 2)
    half = int(2.5 * R_px)
    yy, xx = np.mgrid[-half:half + 1, -half:half + 1]
    k = np.exp(-0.5 * ((np.hypot(xx, yy) - R_px) / (0.2 * R_px)) ** 2)
    k -= k.mean()
    k /= np.linalg.norm(k)
    c = signal.fftconvolve(H, k, mode='same')
    pk = (c == ndi.maximum_filter(c, size=int(2 * R_px) + 1)) & (c > np.percentile(c[H > 0], pct))
    iy, ix = np.nonzero(pk)
    tree, out = cKDTree(np.c_[x, y]), []
    for a, b in zip(ix * pix + pix / 2, iy * pix + pix / 2):
        ok = True
        for _ in range(3):
            idx = tree.query_ball_point([a, b], search)
            P = np.c_[x[idx], y[idx]] - [a, b]
            r = np.hypot(P[:, 0], P[:, 1])
            m = (r > r_in) & (r < r_out)
            if m.sum() < 15:
                ok = False
                break
            Q = P[m]
            sol = np.linalg.lstsq(np.c_[2 * Q[:, 0], 2 * Q[:, 1], np.ones(len(Q))], (Q ** 2).sum(1), rcond=None)[0]
            a, b, R = a + sol[0], b + sol[1], np.sqrt(sol[2] + sol[0] ** 2 + sol[1] ** 2)
        if ok:
            out.append((a, b, R, m.sum()))
    return np.array(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--no-picasso', action='store_true')
    a = ap.parse_args()
    res = {}

    # ---- Heil et al., gp210 NPC, glass, 20000 frames of 117 x 118 px ------------------------------------------------------------------------------
    p = os.path.join(HEIL, 'NPCData_glass_STORM.tif')
    fr, n = read_pages(p, 0, 20000)[0][::200], 20000
    res['heil'] = dict(frames=n, shape=list(fr.shape[1:]), median_adu=float(np.median(fr)), p99_9_adu=float(np.percentile(fr, 99.9)))
    tab = np.genfromtxt(os.path.join(HEIL, 'NPCData_glass_STORM_loc.csv'), delimiter=',', skip_header=1)     # id frame x y sigma intensity offset bkgstd chi2 unc det
    res['heil'].update(table_rows=len(tab), locs_per_frame=len(tab) / n, median_sigma_nm=float(np.median(tab[:, 4])), median_photons=float(np.median(tab[:, 5])),
                       median_uncertainty_nm=float(np.median(tab[:, 9])))
    for pct in (95.0, 99.0):                                    # the count of rings depends on the threshold of the matched filter (it is not a count of NPCs)
        r = rings(tab[:, 2], tab[:, 3], 8.0, 10.0, 12100, 40, 125, 130, pct)
        res['heil']['rings_matched_filter_above_percentile_%d' % pct] = dict(
            rings=len(r), radius_nm_q05_25_50_75_95=np.percentile(r[:, 2], [5, 25, 50, 75, 95]).round(1).tolist(), locs_in_annulus_median=float(np.median(r[:, 3])))
    if not a.no_picasso:
        fr2, _ = read_pages(p, 1000, 1200)
        dets, nloc, info = picasso(fr2, dict(pixel_nm=102.0, baseline=100.0, gain=100.0, sensitivity=15.01, qe=0.9))
        res['heil']['picasso_R0_frames_1000_1199'] = dict(locs=nloc, per_frame=nloc / 200.0)

    # ---- Thevathasan et al., Nup96-SNAP, FOV 17 and 18 ----------------------------------------------------------------------------------------------
    res['nup96'] = {}
    for k in (17, 18):
        p = os.path.join(NUP, NUP_PREFIX % k)
        with tifffile.TiffFile(p) as t:
            npg = len(t.pages) - 1                              # the last page of the prefix is cut
        fr, _ = read_pages(p, 0, npg)
        tab = np.genfromtxt(os.path.join(NUP, 'Cos7_Nup96_A647_COT_Au_%d_1_locs.csv' % k), delimiter=',', skip_header=1)
        ref = tab[tab[:, 0] > 1056]                             # reference = the frames NOT in the prefix
        rr = {}
        for pct in (90.0, 97.0):
            r = rings(ref[:, 12], ref[:, 13], 5.3, 10.0, 45000, 25, 80, 90, pct)
            rr['rings_from_frames_gt_1056_matched_filter_above_percentile_%d' % pct] = dict(
                rings=len(r), radius_nm_q05_25_50_75_95=np.percentile(r[:, 2], [5, 25, 50, 75, 95]).round(1).tolist(), locs_in_annulus_median=float(np.median(r[:, 3])))
        d = dict(frames_in_prefix=int(npg), shape=list(fr.shape[1:]), median_adu=float(np.median(fr)), meta=mm_meta(p), table_frames=int(tab[:, 0].max()),
                 table_rows=len(tab), table_locs_in_prefix_frames=int((tab[:, 0] <= 1056).sum()), table_x_nm_over_x_pix=float(np.median(tab[:, 12] / tab[:, 1])),
                 median_photons=float(np.median(tab[:, 5])), median_crlb_y_nm=float(np.nanmedian(tab[:, 15])), **rr)
        if not a.no_picasso and k == 17:
            dets, nloc, info = picasso(fr[400:500], dict(pixel_nm=106.0, baseline=100.0, gain=1.0, sensitivity=1.0, qe=0.95))
            hit = tot = 0
            for i, dd in enumerate(dets):
                s = tab[tab[:, 0] == 401 + i][:, 12:14]          # the table is 1-based
                if len(dd) and len(s):
                    dist, _ = cKDTree(s).query(dd[:, :2])
                    hit += int((dist < 150).sum())
                    tot += len(dd)
            d['picasso_R0_frames_400_499'] = dict(locs=nloc, per_frame=nloc / 100.0, fraction_within_150nm_of_table_loc_same_frame=hit / max(tot, 1),
                                                  smap_table_locs_per_frame_same_frames=float(((tab[:, 0] >= 401) & (tab[:, 0] <= 500)).sum() / 100.0))
        res['nup96']['fov%d' % k] = d

    # ---- SIMFLUX 80 nm nanorulers ----------------------------------------------------------------------------------------------------------------
    p = os.path.join(SIM, 'sim4_1_MMStack_Pos0.ome.tif.part0-300MB')
    with tifffile.TiffFile(p) as t:
        npg = len(t.pages) - 1
    fr, _ = read_pages(p, 0, npg)
    d = dict(frames_in_prefix=int(npg), shape=list(fr.shape[1:]), median_adu=float(np.median(fr)), meta=mm_meta(p), p99_99_adu=float(np.percentile(fr, 99.99)))
    if not a.no_picasso:
        for ng in (500, 1500, 5000):
            dets, nloc, info = picasso(fr[700:800], dict(pixel_nm=65.0, baseline=100.0, gain=1.0, sensitivity=0.46, qe=0.8, min_net_gradient=ng))
            d['picasso_frames_700_799_min_net_gradient_%d_locs_per_frame' % ng] = nloc / 100.0
    res['simflux'] = d
    print(json.dumps(res, indent=1, default=float))


if __name__ == '__main__':
    main()
