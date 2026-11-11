"""
Crowded-field benchmark helpers (the M2 / SDSS patch against HST): magnitudes, the truth list of a patch, one-to-one matching, completeness and purity by magnitude
=================================================================================================================================================================

The benchmark of Portillo et al. (2017) as defined by the code shipped with its follow-up (`results/astro/data_access.md`, section 3.3 and 5): a 100 x 100 pixel
patch of an SDSS r-band frame of the globular cluster M2, the truth an HST catalogue of the same stars (F606W used as r: zero point 22.5, 0.00546689 nanomaggies
per count); a detection matches a truth source when it is within 0.5 px and 0.5 mag; truth sources within 2.5 px of the patch edge are not scored; the figures are
by magnitude bins of 0.5 mag from 15.0 to 23.5. Positions are in the patch's pixel coordinates (0-based, pixel centre at the integer, x the column).

Only numpy and scipy; the truth list itself is built by `data/astro/make_truth.py` (it needs astropy).
"""
import numpy as np
from scipy.optimize import linear_sum_assignment

ZP = 22.5
NMGY_PER_COUNT = 0.00546689
BIN_EDGES = np.arange(15.0, 23.5 + 1e-9, 0.5)          # 17 bins
DPOS, DMAG, EDGE = 0.5, 0.5, 2.5


def counts_to_mag(flux):
    """SDSS-style magnitude of a flux in counts (DN): 22.5 - 2.5 log10(flux x 0.00546689); NaN for a non-positive flux."""
    flux = np.asarray(flux, float)
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(flux > 0, ZP - 2.5 * np.log10(np.where(flux > 0, flux, 1.0) * NMGY_PER_COUNT), np.nan)


def mag_to_counts(mag):
    return 10.0 ** (-(np.asarray(mag, float) - ZP) / 2.5) / NMGY_PER_COUNT


def truth_for_patch(truth, x0, y0, n=100, mag_max=30.0):
    """The truth sources inside the patch box [x0 - 0.5, x0 + n - 0.5) x [y0 - 0.5, y0 + n - 0.5): dict(x, y, mag) in patch pixel coordinates and
    `inner`, the boolean mask of those at least EDGE pixels from every edge (the scored ones)."""
    px, py, m = truth['px'], truth['py'], truth['m606']
    sel = (px >= x0 - 0.5) & (px < x0 + n - 0.5) & (py >= y0 - 0.5) & (py < y0 + n - 0.5) & (m < mag_max)
    x, y, mag = px[sel] - x0, py[sel] - y0, m[sel]
    inner = (x >= EDGE - 0.5) & (x <= n - 0.5 - EDGE) & (y >= EDGE - 0.5) & (y <= n - 0.5 - EDGE)
    return dict(x=x, y=y, mag=mag, inner=inner)


def match(det_xy, det_mag, tr_xy, tr_mag, dpos=DPOS, dmag=DMAG):
    """One-to-one matching of detections to truth sources: among all pairings of pairs within dpos pixels AND dmag magnitudes it takes one of maximum size and,
    among those, minimum total distance. Returns (det index, truth index) arrays of the matched pairs."""
    det_xy, tr_xy = np.atleast_2d(det_xy).reshape(-1, 2), np.atleast_2d(tr_xy).reshape(-1, 2)
    if len(det_xy) == 0 or len(tr_xy) == 0:
        return np.zeros(0, int), np.zeros(0, int)
    d = np.linalg.norm(det_xy[:, None, :] - tr_xy[None, :, :], axis=2)
    dm = np.abs(np.asarray(det_mag, float)[:, None] - np.asarray(tr_mag, float)[None, :])
    with np.errstate(invalid='ignore'):
        ok = (d <= dpos) & (dm <= dmag)
    if not ok.any():
        return np.zeros(0, int), np.zeros(0, int)
    big = 4.0 * (min(d.shape) + 1) * dpos + 1.0
    r, c = linear_sum_assignment(np.where(ok, d, big))
    keep = ok[r, c]
    return r[keep], c[keep]


CELL = 20                                               # px: the side of the square blocks of the block bootstrap (25 blocks per 100 x 100 patch)


def block_of(xy, n=100, cell=CELL):
    """The index (row-major over the n / cell x n / cell blocks of a patch) of the block of each position (x, y) in patch pixel coordinates (pixel centres at the integers)."""
    xy = np.atleast_2d(np.asarray(xy, float))
    nb = n // cell
    ix = np.clip(np.floor((xy[:, 0] + 0.5) / cell).astype(int), 0, nb - 1)
    iy = np.clip(np.floor((xy[:, 1] + 0.5) / cell).astype(int), 0, nb - 1)
    return iy * nb + ix


def score(det_xy, det_flux, truth_patch, n=100):
    """Completeness and purity of a catalogue against a patch's truth. Detections within EDGE pixels of the edge are not scored (as the truth). Returns a dict:
    per bin (BIN_EDGES): n_true, n_matched (truth sources matched), n_det, n_det_matched (detections matched) and the derived completeness and false-positive fraction,
    and totals for mag < 22 (the benchmark's depth) and for all magnitudes, plus the median position offset of the matched bright (mag < 21) pairs, and `blk22`
    (n_blocks, 4): the counts true22, matched22, det22, det22_matched per block of CELL x CELL pixels (a truth source and its match are counted in the truth's block, a
    detection and its match in the detection's), the units of the block bootstrap."""
    det_xy = np.atleast_2d(np.asarray(det_xy, float)).reshape(-1, 2)
    det_mag = counts_to_mag(det_flux)
    ok = np.isfinite(det_mag)
    inside = (det_xy[:, 0] >= EDGE - 0.5) & (det_xy[:, 0] <= n - 0.5 - EDGE) & (det_xy[:, 1] >= EDGE - 0.5) & (det_xy[:, 1] <= n - 0.5 - EDGE)
    d_keep = ok & inside
    t = truth_patch
    tk = t['inner']
    ti = np.flatnonzero(tk)
    di = np.flatnonzero(d_keep)
    rd, rt = match(det_xy[di], det_mag[di], np.c_[t['x'][ti], t['y'][ti]], t['mag'][ti])
    md, mt = di[rd], ti[rt]
    tmag, dmag_ = t['mag'][ti], det_mag[di]
    n_bins = len(BIN_EDGES) - 1
    nt = np.histogram(tmag, BIN_EDGES)[0]
    nm = np.histogram(t['mag'][mt], BIN_EDGES)[0]
    nd = np.histogram(dmag_, BIN_EDGES)[0]
    ndm = np.histogram(det_mag[md], BIN_EDGES)[0]
    with np.errstate(invalid='ignore', divide='ignore'):
        comp = nm / np.maximum(nt, 1)
        fpf = 1.0 - ndm / np.maximum(nd, 1)
    comp[nt == 0] = np.nan
    fpf[nd == 0] = np.nan
    bright = t['mag'][mt] < 21.0
    off = np.median(np.c_[det_xy[md, 0] - t['x'][mt], det_xy[md, 1] - t['y'][mt]][bright], axis=0) if bright.any() else np.full(2, np.nan)
    sel22t, sel22d, sel22m = tmag < 22.0, dmag_ < 22.0, det_mag[md] < 22.0
    nbl = (n // CELL) ** 2
    blk = np.zeros((nbl, 4))
    t_xy = np.c_[t['x'], t['y']]
    np.add.at(blk[:, 0], block_of(t_xy[ti[sel22t]], n), 1)
    sel_mt = t['mag'][mt] < 22.0
    np.add.at(blk[:, 1], block_of(t_xy[mt[sel_mt]], n), 1)
    np.add.at(blk[:, 2], block_of(det_xy[di[sel22d]], n), 1)
    np.add.at(blk[:, 3], block_of(det_xy[md[sel22m]], n), 1)
    return dict(n_true=nt, n_matched=nm, n_det=nd, n_det_matched=ndm, completeness=comp, fp_fraction=fpf,
                true22=int(sel22t.sum()), matched22=int((t['mag'][mt] < 22.0).sum()), det22=int(sel22d.sum()), det22_matched=int(sel22m.sum()),
                true_all=int(len(ti)), matched_all=int(len(mt)), det_all=int(len(di)), offset=off, blk22=blk)
