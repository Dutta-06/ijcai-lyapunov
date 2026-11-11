"""
A second globular cluster as an independent crowded-field test set: SDSS r-band patches against the HST (HUGS) catalogue, in the on-disk format of the M2 patches
==================================================================================================================================================================

The M2 benchmark (`exp_astro_m2.py`, `astro_common.py`, `data/astro/fetch.py`) lives on one cluster and one frame, and its configuration is frozen on two patches of that frame. This
module builds the test set that a statement about generalisation needs: 100 x 100 pixel patches of another cluster's SDSS r-band frame with the HUGS catalogue of the same stars as the
truth, in exactly the files `fetch.py extract-patch` writes (`patch_arrays.npz` with dn, electrons, nmgy, sky_dn, psf51, psf25, gain, dark_variance, x0, y0; `patch_dn.fits`; `meta.json`;
the FITS header carries one extra card, NMGY) and a truth list in the format of `truth_hugs_meth1.npz` (px, py, m606, e606, nf606 in frame pixels). Nothing here runs a cataloguing
method: only geometry, units, truth lists, counts and the choice of the patches. The findings are in `results/astro/second_cluster_spec.md`.

Conventions (those of M2). Pixel (x, y) = (column, row) of the frame, 0-based, the pixel centre at the integer. Truth magnitude = the HUGS F606W used as r (zero point 22.5, no colour
term). DN = frame / calibration vector + the bilinearly interpolated sky (bias-subtracted, flat-fielded counts); electrons = DN x gain; gain and dark variance from the psField (r); the PSF
is the psField reconstruction at the patch centre, 51 x 51, unit sum. A truth source belongs to a patch when it is inside [x0 - 0.5, x0 + n - 0.5) x [y0 - 0.5, y0 + n - 0.5) and is scored
when it is at least 2.5 px from the edges (`astro_common.truth_for_patch`).

Sky to pixel. The frame header carries a TAN projection (CRVAL, CRPIX, CD), a linearisation of the frame's full astrometric solution (the asTrans row in HDU 3) at the frame centre.
`tan_world2pix` is that header mapping (the mapping of `data/astro/make_truth.py` and of the M2 truth); `astrans_world2pix` is the full solution of the SDSS data model (great-circle
coordinates, cubic distortion in the column, optional colour term). The truth lists written here use the full solution (colour term off: it moves positions by < 0.03 px), because the
image prefers it away from the frame centre (spec, section 3); the header-mapped list is written next to it (`..._header.npz`) for parity with M2.

Choice of the patches (made by location, truth and image counts only; no detector output enters). Origins on a 10-pixel grid with 0.03 to 0.15 truth sources per pixel brighter than
F606W = 22 that also (1) keep 40 px from the frame edges, (2) lie inside the HUGS footprint (every 20 x 20 px cell of the box grown by 20 px holds a truth source), (3) have no pixel above
0.8 x the photoField saturation level within 10 px of the box, (4) have no truth source brighter than F606W = 14.5 within 10 px (SDSS saturates near 14; HUGS flags such stars saturated),
and (5) have no HUGS source with F814W < 17 but no F606W within 10 px (a bright star the truth lacks). From the four density classes [0.03, 0.05), [0.05, 0.08), [0.08, 0.11), [0.11, 0.15]
the candidate closest to the class target (0.04, 0.065, 0.095, 0.13) that does not overlap a patch already chosen, the scarcest class first, ties by (y0, x0).

Command line (needs astropy; run it as `python -I` because it reads downloaded, untrusted files; $ASTRO_DATA_2 holds the output of `data/astro/fetch_second.py fetch`):
    python -I experiments/astro_second.py survey                         the coverage table: clusters, SDSS runs, footprint fraction, candidate windows
    python -I experiments/astro_second.py validate-m2                    the alignment and units validation on M2 against the shipped M2 truth and patches
    python -I experiments/astro_second.py build --cluster ngc5024        truth lists, selected patches and their table (clusters in TARGETS; --dry-run writes nothing)
"""
import argparse
import hashlib
import json
import os
import warnings

import numpy as np
from scipy.ndimage import map_coordinates
from scipy.ndimage import shift as nd_shift
from scipy.spatial import cKDTree

ZP = 22.5
PIX = 0.396                                   # arcsec per SDSS pixel
FRAME_SHAPE = (1489, 2048)                    # rows, columns
EDGE = 2.5                                    # px, the scoring margin of astro_common
N_PATCH = 100
MAG_DEPTH = 22.0
CLASSES = ((0.03, 0.05), (0.05, 0.08), (0.08, 0.11), (0.11, 0.15))
TARGET_DENSITY = (0.04, 0.065, 0.095, 0.13)
STEP, FRAME_MARGIN, CELL, COLLAR = 10, 40, 20, 20
SAT_FRACTION, DN_MARGIN, BRIGHT_LIMIT, MISSING_LIMIT = 0.8, 10, 14.5, 17.0
M2_SNR22 = 6.33                               # peak S/N of a F606W = 22 point source in the M2 benchmark patch (validate-m2, section 6)
HUGS_COLS = (20, 21, 24, 33, 34, 26)          # 0-based: F606W magnitude, its rms, number of exposures it is found in, RA, Dec, F814W magnitude
TARGETS = {'ngc6341': (4682, 6, 208), 'ngc5024': (5390, 6, 162), 'ngc5904': (1458, 4, 700), 'ngc5272': (4646, 3, 81)}       # cluster -> (run, camcol, field) of the frame that holds the patches
ROOT2 = os.environ.get('ASTRO_DATA_2', '/home/user/astro_data_2')
ROOT = os.environ.get('ASTRO_DATA', '/home/user/astro_data')


# ----------------------------------------------------------------------------- sky geometry
def tan_params(hdr):
    """The TAN solution of a frame header (any mapping with CRVAL1/2, CRPIX1/2, CD1_1, CD1_2, CD2_1, CD2_2)."""
    g = lambda k: float(hdr[k])
    return dict(crval=(g('CRVAL1'), g('CRVAL2')), crpix=(g('CRPIX1'), g('CRPIX2')), cd=((g('CD1_1'), g('CD1_2')), (g('CD2_1'), g('CD2_2'))))


def tan_world2pix(ra, dec, p):
    """(RA, Dec) in degrees -> (x, y) = (column, row), 0-based with the pixel centre at the integer (FITS pixels are 1-based: x = u + CRPIX1 - 1 with (u, v) = CD^-1 (xi, eta))."""
    ra, dec = np.broadcast_arrays(np.asarray(ra, float), np.asarray(dec, float))
    a0, d0 = np.deg2rad(p['crval'][0]), np.deg2rad(p['crval'][1])
    a, d = np.deg2rad(ra.ravel()), np.deg2rad(dec.ravel())
    cosc = np.sin(d0) * np.sin(d) + np.cos(d0) * np.cos(d) * np.cos(a - a0)
    xi = np.rad2deg(np.cos(d) * np.sin(a - a0) / cosc)
    eta = np.rad2deg((np.cos(d0) * np.sin(d) - np.sin(d0) * np.cos(d) * np.cos(a - a0)) / cosc)
    uv = np.linalg.solve(np.array(p['cd'], float), np.vstack([xi, eta]))
    return (uv[0] + p['crpix'][0] - 1.0).reshape(ra.shape), (uv[1] + p['crpix'][1] - 1.0).reshape(ra.shape)


def tan_pix2world(x, y, p):
    """Inverse of `tan_world2pix`: (x, y) -> (RA, Dec) in degrees."""
    x, y = np.broadcast_arrays(np.asarray(x, float), np.asarray(y, float))
    a0, d0 = np.deg2rad(p['crval'][0]), np.deg2rad(p['crval'][1])
    xi_eta = np.deg2rad(np.array(p['cd'], float) @ np.vstack([x.ravel() + 1.0 - p['crpix'][0], y.ravel() + 1.0 - p['crpix'][1]]))
    xi, eta = xi_eta
    rho = np.hypot(xi, eta)
    c = np.arctan(rho)
    with np.errstate(invalid='ignore', divide='ignore'):
        sin_d = np.cos(c) * np.sin(d0) + np.where(rho > 0, eta * np.sin(c) * np.cos(d0) / np.where(rho > 0, rho, 1.0), 0.0)
    dec = np.arcsin(np.clip(sin_d, -1.0, 1.0))
    ra = a0 + np.arctan2(xi * np.sin(c), rho * np.cos(d0) * np.cos(c) - eta * np.sin(d0) * np.sin(c))
    return (np.rad2deg(ra) % 360.0).reshape(x.shape), np.rad2deg(dec).reshape(x.shape)


def eq_to_gc(ra, dec, node, incl):
    """Equatorial (RA, Dec) to the great-circle coordinates (mu, nu) of a scan with ascending node `node` and inclination `incl` (degrees)."""
    ra, dec, node, incl = (np.deg2rad(np.asarray(v, float)) for v in (ra, dec, node, incl))
    x = np.cos(ra - node) * np.cos(dec)
    y = np.sin(ra - node) * np.cos(dec) * np.cos(incl) + np.sin(dec) * np.sin(incl)
    z = -np.sin(ra - node) * np.cos(dec) * np.sin(incl) + np.sin(dec) * np.cos(incl)
    return (np.rad2deg(node) + np.rad2deg(np.arctan2(y, x))) % 360.0, np.rad2deg(np.arcsin(np.clip(z, -1.0, 1.0)))


def gc_to_eq(mu, nu, node, incl):
    """Inverse of `eq_to_gc`; mu may lie outside 0..360 (the asTrans constant `a` does)."""
    mu, nu, node, incl = (np.deg2rad(np.asarray(v, float)) for v in (mu, nu, node, incl))
    x, y, z = np.cos(mu - node) * np.cos(nu), np.sin(mu - node) * np.cos(nu), np.sin(nu)
    y2, z2 = y * np.cos(incl) - z * np.sin(incl), y * np.sin(incl) + z * np.cos(incl)
    return (np.rad2deg(np.arctan2(y2, x)) + np.rad2deg(node)) % 360.0, np.rad2deg(np.arcsin(np.clip(z2, -1.0, 1.0)))


ASTRANS_KEYS = ('NODE', 'INCL', 'A', 'B', 'C', 'D', 'E', 'F', 'DROW0', 'DROW1', 'DROW2', 'DROW3', 'DCOL0', 'DCOL1', 'DCOL2', 'DCOL3', 'CSROW', 'CSCOL', 'CCROW', 'CCCOL', 'RICUT')


def astrans_params(rec):
    """The constants of the asTrans row of a frame (HDU 3) as floats."""
    return {k: float(rec[k]) for k in ASTRANS_KEYS}


def _distortion(col, t, color):
    drow = t['DROW0'] + t['DROW1'] * col + t['DROW2'] * col ** 2 + t['DROW3'] * col ** 3
    dcol = t['DCOL0'] + t['DCOL1'] * col + t['DCOL2'] * col ** 2 + t['DCOL3'] * col ** 3
    if color is not None:                                    # the colour term of the data model: r - i < riCut: slope x colour; else a constant
        color = np.asarray(color, float)
        blue = color < t['RICUT']
        drow = drow + np.where(blue, t['CSROW'] * color, t['CCROW'])
        dcol = dcol + np.where(blue, t['CSCOL'] * color, t['CCCOL'])
    return drow, dcol


def astrans_pix2world(x, y, t, color=None):
    """The asTrans solution of the SDSS data model: (row, col) of the SDSS convention (pixel centres at half-integers: row = y + 0.5, col = x + 0.5) ->
    rowm = row + dRow(col), colm = col + dCol(col); mu = a + b rowm + c colm; nu = d + e rowm + f colm; then (mu, nu) -> (RA, Dec). Returns (RA, Dec) in degrees."""
    row, col = np.asarray(y, float) + 0.5, np.asarray(x, float) + 0.5
    drow, dcol = _distortion(col, t, color)
    rowm, colm = row + drow, col + dcol
    return gc_to_eq(t['A'] + t['B'] * rowm + t['C'] * colm, t['D'] + t['E'] * rowm + t['F'] * colm, t['NODE'], t['INCL'])


def astrans_world2pix(ra, dec, t, color=None, n_iter=8):
    """Inverse of `astrans_pix2world`: (mu, nu) -> (rowm, colm) by the 2 x 2 linear system, then the column distortion by fixed-point iteration (it is about 1e-3 per pixel)."""
    ra, dec = np.broadcast_arrays(np.asarray(ra, float), np.asarray(dec, float))
    mu, nu = eq_to_gc(ra.ravel(), dec.ravel(), t['NODE'], t['INCL'])
    dmu = (mu - t['A'] + 180.0) % 360.0 - 180.0
    rowm, colm = np.linalg.solve(np.array([[t['B'], t['C']], [t['E'], t['F']]]), np.vstack([dmu, nu - t['D']]))
    col = colm.copy()
    for _ in range(n_iter):
        col = colm - _distortion(col, t, color)[1]
    row = rowm - _distortion(col, t, color)[0]
    return (col - 0.5).reshape(ra.shape), (row - 0.5).reshape(ra.shape)


# ----------------------------------------------------------------------------- reading the downloaded files (astropy)
def read_frame(path):
    """The SDSS frame file (astropy reads .bz2 transparently): header, image (nanomaggies, rows x columns), calibration vector, sky table, asTrans row."""
    from astropy.io import fits
    with fits.open(path, memmap=False) as h:
        t = h[2].data
        return dict(header=h[0].header.copy(), nmgy=h[0].data.astype(np.float64), cal=h[1].data.astype(np.float64),
                    allsky=t['ALLSKY'][0].astype(np.float64), xinterp=t['XINTERP'][0].astype(np.float64), yinterp=t['YINTERP'][0].astype(np.float64),
                    astrans=astrans_params(h[3].data[0]))


def read_psfield(path):
    """The psField: the r-band KL eigen-images (HDU 3) and the r entries of gain, dark variance, sky and sky sigma (HDU 6)."""
    from astropy.io import fits
    with fits.open(path, memmap=False) as h:
        kl, row6 = h[3].data, h[6].data
        return dict(kl=kl, gain=float(row6['gain'][0][2]), dark_variance=float(row6['dark_variance'][0][2]), sky=float(row6['sky'][0][2]), skysig=float(row6['skysig'][0][2]))


def read_hugs(path):
    """The HUGS catalogue (method 1): arrays m606, e606, nf606, ra, dec, f814 of every row (m606 = -99.9999 where F606W is not measured, f814 likewise for F814W)."""
    a = np.loadtxt(path, comments='#', usecols=HUGS_COLS)
    return dict(m606=a[:, 0], e606=a[:, 1], nf606=a[:, 2], ra=a[:, 3], dec=a[:, 4], f814=a[:, 5])


# ----------------------------------------------------------------------------- units: the frame to DN, the psField to the PSF
def sky_image(allsky, xinterp, yinterp):
    """The sky (DN) on the full frame: bilinear interpolation of ALLSKY at (YINTERP, XINTERP), constant beyond the table (the data model)."""
    yy, xx = np.meshgrid(yinterp, xinterp, indexing='ij')
    return map_coordinates(allsky, [yy, xx], order=1, mode='nearest')


def frame_to_dn(nmgy, cal, sky):
    """dn = image / calibration + sky (SDSS data model): the bias-subtracted, flat-fielded counts of the original image."""
    return nmgy / cal[None, :] + sky


def dn_variance(dn, gain, dark_variance):
    """Variance of a DN pixel: electrons = DN x gain are Poisson, plus the dark variance (data model): DN / gain + dark variance."""
    return np.asarray(dn, float) / gain + dark_variance


def counts_to_mag(flux, nmgy_per_count):
    """22.5 - 2.5 log10(flux x nmgy per count); NaN for a non-positive flux. astro_common uses the M2 frame's 0.00546689; the frame's own NMGY differs by < 0.03 mag."""
    flux = np.asarray(flux, float)
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(flux > 0, ZP - 2.5 * np.log10(np.where(flux > 0, flux, 1.0) * nmgy_per_count), np.nan)


def sdss_psf_recon(kl, xpos, ypos):
    """KL-basis PSF image (51 x 51) at 0-based (column xpos, row ypos); the centre pixel is the star centre. The recipe of the SDSS psField data model (polynomials in
    (0.001 (row + 0.5), 0.001 (col + 0.5))), the same code as `data/astro/fetch.py`."""
    rc = 0.001
    rnrow, rncol = int(kl['RNROW'][0]), int(kl['RNCOL'][0])
    nr, nc = int(kl['nrow_b'].max()), int(kl['ncol_b'].max())
    coeffs = np.outer(((ypos + 0.5) * rc) ** np.arange(nr), ((xpos + 0.5) * rc) ** np.arange(nc))
    img = np.zeros(len(kl['RROWS'][0]), dtype=np.float64)
    for i in range(len(kl)):
        c = np.asarray(kl[i]['c'])[: int(kl[0]['nrow_b']), : int(kl[0]['ncol_b'])]
        img += np.asarray(kl[i]['RROWS'], dtype=np.float64) * (c * coeffs[: c.shape[0], : c.shape[1]].T).sum()
    return img[: rnrow * rncol].reshape(rncol, rnrow)


def unit_psf(psf):
    """Negative pixels to zero, unit sum (as `fetch.py extract-patch`)."""
    psf = np.clip(psf, 0, None)
    return psf / psf.sum()


def psf_fwhm(psf):
    """FWHM (pixels) of a PSF stamp whose centre pixel is the star centre: twice the radius at which the azimuthal mean of the PSF about the centre (32 directions, cubic interpolation of
    the stamp, steps of 0.02 px) falls to half of the centre pixel, linearly interpolated between the two steps."""
    c = psf.shape[0] // 2
    r = np.arange(0.0, c, 0.02)
    th = np.linspace(0.0, 2.0 * np.pi, 32, endpoint=False)
    prof = map_coordinates(psf, [c + r[:, None] * np.sin(th)[None, :], c + r[:, None] * np.cos(th)[None, :]], order=3).mean(axis=1) / psf[c, c]
    i = int(np.argmax(prof < 0.5))
    return float(2.0 * (r[i - 1] + (prof[i - 1] - 0.5) / (prof[i - 1] - prof[i]) * (r[i] - r[i - 1])))


def psf_vs_pydl(kl, xpos, ypos):
    """Largest absolute difference between `sdss_psf_recon` (normalised as `unit_psf`) and the PSF of the independent `pydl.photoop.image.sdss_psf_recon` at the same position;
    None when pydl is not installed."""
    try:
        from pydl.photoop.image import sdss_psf_recon as pydl_recon
    except ImportError:
        return None
    return float(np.abs(unit_psf(sdss_psf_recon(kl, xpos, ypos)) - unit_psf(np.asarray(pydl_recon(kl, xpos, ypos), float))).max())


def patch_name(run, camcol, field, x0, y0, n=N_PATCH):
    return 'patch_r_%06d_%d_%04d_x%d_y%d_%dx%d' % (run, camcol, field, x0, y0, n, n)


def extract_patch(frame, ps, x0, y0, n=N_PATCH, dn_full=None, sky=None):
    """The arrays of `patch_arrays.npz`, the content of `meta.json` and the extra header cards of `patch_dn.fits` for the patch with its lower-left pixel at (column x0, row y0).
    frame from `read_frame`, ps from `read_psfield`; dn_full / sky may be passed to avoid recomputing the full frame for several patches."""
    if sky is None:
        sky = sky_image(frame['allsky'], frame['xinterp'], frame['yinterp'])
    if dn_full is None:
        dn_full = frame_to_dn(frame['nmgy'], frame['cal'], sky)
    sl = (slice(y0, y0 + n), slice(x0, x0 + n))
    dn = dn_full[sl]
    psf = unit_psf(sdss_psf_recon(ps['kl'], x0 + (n - 1) / 2.0, y0 + (n - 1) / 2.0))
    c = psf.shape[0] // 2
    psf25 = psf[c - 12: c + 13, c - 12: c + 13]
    psf25 = psf25 / psf25.sum()
    gain, dark = ps['gain'], ps['dark_variance']
    arrays = dict(dn=dn, electrons=dn * gain, nmgy=frame['nmgy'][sl], sky_dn=sky[sl], psf51=psf, psf25=psf25, gain=gain, dark_variance=dark, x0=x0, y0=y0)
    meta = dict(x0=x0, y0=y0, size=n, gain=gain, dark_variance=dark, psfield_sky=ps['sky'], psfield_skysig=ps['skysig'],
                note='origin is the lower-left pixel, 0-based (column x0, row y0); array index [row, column]')
    return arrays, meta


def write_patch(outdir, arrays, meta, frame_header):
    """patch_arrays.npz, meta.json and patch_dn.fits exactly as `fetch.py extract-patch` writes them (the FITS header: the frame WCS moved to the patch, units, run/rerun/camcol/frame/filter,
    the origin PX0/PY0, GAINR; plus NMGY, the frame's calibration in nanomaggies per count, which the M2 header does not carry)."""
    from astropy.io import fits
    from astropy.wcs import WCS
    os.makedirs(outdir, exist_ok=True)
    x0, y0 = arrays['x0'], arrays['y0']
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')                      # the frame header's RADECSYS, MJD-OBS and the length of the BUNIT comment (as in fetch.py)
        w = WCS(frame_header)
        w.wcs.crpix = w.wcs.crpix - np.array([x0, y0])
        ph = fits.PrimaryHDU(arrays['dn'].astype(np.float32), header=w.to_header())
        ph.header['BUNIT'] = ('DN', 'bias-subtracted flat-fielded counts (frame/cal + sky)')
        for k in ('RUN', 'RERUN', 'CAMCOL', 'FRAME', 'FILTER'):
            ph.header[k] = frame_header[k]
        ph.header['PX0'], ph.header['PY0'] = x0, y0
        ph.header['GAINR'] = (arrays['gain'], 'e-/DN, psField HDU6')
        ph.header['NMGY'] = (float(frame_header['NMGY']), 'nMgy per count of the frame')
        fits.HDUList([ph]).writeto(os.path.join(outdir, 'patch_dn.fits'), overwrite=True)
    np.savez(os.path.join(outdir, 'patch_arrays.npz'), **arrays)
    with open(os.path.join(outdir, 'meta.json'), 'w') as f:
        json.dump(meta, f, indent=1)


def content_sha256(path):
    """Hash of the content of an .npz (names, dtypes, shapes, bytes of the arrays in name order): the zip container itself carries time stamps and is not reproducible."""
    h = hashlib.sha256()
    with np.load(path) as z:
        for k in sorted(z.files):
            a = np.ascontiguousarray(z[k])
            h.update(k.encode() + str(a.dtype).encode() + str(a.shape).encode() + a.tobytes())
    return h.hexdigest()


def file_sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for blk in iter(lambda: f.read(1 << 20), b''):
            h.update(blk)
    return h.hexdigest()


# ----------------------------------------------------------------------------- the truth list
def make_truth(cat, world2pix):
    """The HUGS sources with a valid F606W in frame pixels: dict(px, py, m606, e606, nf606), the format of `truth_hugs_meth1.npz` (`data/astro/make_truth.py`)."""
    ok = (cat['m606'] > -90) & (cat['m606'] < 90)
    px, py = world2pix(cat['ra'][ok], cat['dec'][ok])
    return dict(px=px, py=py, m606=cat['m606'][ok], e606=cat['e606'][ok], nf606=cat['nf606'][ok])


def unmeasured_f606w(cat, world2pix, f814_limit=22.0):
    """The HUGS sources whose F814W is measured (brighter than f814_limit) but whose F606W is not (saturated in F606W or lost): they are not in the truth list. dict(px, py, f814)."""
    ok = (cat['m606'] <= -90) & (cat['f814'] > -90) & (cat['f814'] < f814_limit)
    x, y = world2pix(cat['ra'][ok], cat['dec'][ok])
    return dict(px=x, py=y, f814=cat['f814'][ok])


def truth_for_patch(truth, x0, y0, n=N_PATCH, mag_max=30.0):
    """The truth sources inside [x0 - 0.5, x0 + n - 0.5) x [y0 - 0.5, y0 + n - 0.5) in patch pixel coordinates (x, y, mag) and `inner`, those at least EDGE px from every edge
    (the scored ones): the same rule as `astro_common.truth_for_patch`."""
    px, py, m = truth['px'], truth['py'], truth['m606']
    sel = (px >= x0 - 0.5) & (px < x0 + n - 0.5) & (py >= y0 - 0.5) & (py < y0 + n - 0.5) & (m < mag_max)
    x, y, mag = px[sel] - x0, py[sel] - y0, m[sel]
    inner = (x >= EDGE - 0.5) & (x <= n - 0.5 - EDGE) & (y >= EDGE - 0.5) & (y <= n - 0.5 - EDGE)
    return dict(x=x, y=y, mag=mag, inner=inner)


# ----------------------------------------------------------------------------- crowding
def patch_density(truth, x0, y0, n=N_PATCH, mag_max=MAG_DEPTH):
    """Truth sources brighter than mag_max per pixel in the patch box (count over n x n)."""
    t = truth_for_patch(truth, x0, y0, n, mag_max)
    return len(t['x']) / float(n * n)


def peak_snr(mag, psf_peak, nmgy_per_count, sky_dn, gain, dark_variance):
    """Peak signal to noise (per pixel, above the sky) of an isolated point source of magnitude `mag` (zero point 22.5): counts = 10^(-0.4 (mag - 22.5)) / nmgy per count, peak = counts x
    the peak of the unit-sum PSF, noise = sqrt((sky + peak) / gain + dark variance) in DN."""
    peak = 10.0 ** (-0.4 * (mag - ZP)) / nmgy_per_count * psf_peak
    return float(peak / np.sqrt((sky_dn + peak) / gain + dark_variance))


def equal_snr_magnitude(target_snr, psf_peak, nmgy_per_count, sky_dn, gain, dark_variance):
    """The magnitude at which `peak_snr` equals target_snr (e.g. the M2 value at F606W = 22, M2_SNR22): the depth of a patch expressed on the M2 scale."""
    from scipy.optimize import brentq
    return float(brentq(lambda m: peak_snr(m, psf_peak, nmgy_per_count, sky_dn, gain, dark_variance) - target_snr, 12.0, 28.0))


def blend_fraction(x, y, radius):
    """Fraction of the sources that have another source of the list closer than `radius` pixels (the crowding as seen by a point-spread function of that width)."""
    xy = np.c_[np.asarray(x, float), np.asarray(y, float)]
    if len(xy) < 2:
        return 0.0
    d, _ = cKDTree(xy).query(xy, k=2)
    return float(np.mean(d[:, 1] < radius))


def crowding(truth, x0, y0, fwhm, n=N_PATCH, mag_max=MAG_DEPTH):
    """The crowding descriptors of a patch: counts brighter than mag_max (box and scored), density per pixel, sources per PSF area (density x pi (FWHM / 2)^2) and the fraction of those
    sources with another one closer than the M2 PSF FWHM (2.22 px) and than the patch's own FWHM."""
    t = truth_for_patch(truth, x0, y0, n, mag_max)
    d = len(t['x']) / float(n * n)
    return dict(n=len(t['x']), n_scored=int(t['inner'].sum()), density=d, per_psf_area=d * np.pi * (fwhm / 2.0) ** 2, blend_m2_fwhm=blend_fraction(t['x'], t['y'], 2.22), blend_own_fwhm=blend_fraction(t['x'], t['y'], fwhm))


def block_counts(truth, x0, y0, n=N_PATCH, cell=CELL, mag_max=24.0):
    """Truth sources brighter than mag_max in each cell x cell block of the patch box, array (n // cell, n // cell). A block far below its neighbours is a hole of the truth (the box
    leaves the HUGS footprint): it makes every detection there a false positive by construction."""
    t = truth_for_patch(truth, x0, y0, n, mag_max)
    h = np.zeros((n // cell, n // cell), int)
    np.add.at(h, (np.clip(np.floor((t['y'] + 0.5) / cell).astype(int), 0, n // cell - 1), np.clip(np.floor((t['x'] + 0.5) / cell).astype(int), 0, n // cell - 1)), 1)
    return h


def window_counts(px, py, m, shape, mag_max=MAG_DEPTH, n=N_PATCH, step=STEP):
    """Truth sources brighter than mag_max inside every n x n box [x0 - 0.5, x0 + n - 0.5) x [y0 - 0.5, y0 + n - 0.5) whose origin is a multiple of `step` inside the frame:
    returns (x0 values, y0 values, counts[iy, ix]) by an integral image over bins of `step` pixels."""
    ny, nx = shape
    ex = np.arange(-0.5, nx + step, step)
    ey = np.arange(-0.5, ny + step, step)
    sel = m < mag_max
    h, _, _ = np.histogram2d(py[sel], px[sel], [ey, ex])
    cs = np.zeros((h.shape[0] + 1, h.shape[1] + 1))
    cs[1:, 1:] = h.cumsum(0).cumsum(1)
    k = n // step
    x0s = np.arange(0, nx - n + 1, step)
    y0s = np.arange(0, ny - n + 1, step)
    ix, iy = x0s // step, y0s // step
    c = cs[np.ix_(iy + k, ix + k)] - cs[np.ix_(iy, ix + k)] - cs[np.ix_(iy + k, ix)] + cs[np.ix_(iy, ix)]
    return x0s, y0s, c


def footprint_cells(px, py, shape, cell=CELL):
    """Boolean occupancy of the frame in cells of `cell` px (cell index floor((x + 0.5) / cell)): True where the truth list has at least one source."""
    ny, nx = shape
    occ = np.zeros((ny // cell + 2, nx // cell + 2), bool)
    ix = np.floor((px + 0.5) / cell).astype(int)
    iy = np.floor((py + 0.5) / cell).astype(int)
    ok = (ix >= 0) & (iy >= 0) & (ix < occ.shape[1]) & (iy < occ.shape[0])
    occ[iy[ok], ix[ok]] = True
    return occ


def inside_footprint(occ, x0, y0, n=N_PATCH, cell=CELL, collar=COLLAR):
    """True when every cell of the box grown by `collar` px on all sides holds at least one truth source (the box is inside the HUGS coverage with a margin)."""
    lo_x, hi_x = int(np.floor((x0 - collar + 0.5) / cell)), int(np.floor((x0 + n - 1 + collar + 0.5) / cell))
    lo_y, hi_y = int(np.floor((y0 - collar + 0.5) / cell)), int(np.floor((y0 + n - 1 + collar + 0.5) / cell))
    if lo_x < 0 or lo_y < 0 or hi_x >= occ.shape[1] or hi_y >= occ.shape[0]:
        return False
    return bool(occ[lo_y: hi_y + 1, lo_x: hi_x + 1].all())


def screen_patches(truth, dn_full, sat_level, shape=FRAME_SHAPE, n=N_PATCH, missing=None):
    """Every origin of the 10-pixel grid whose density is in the range of the classes, with the outcome of each rule: list of dict(x0, y0, count, density, dn_max, margin, footprint,
    saturation, bright, lost) (booleans: True = passes). `dn_full` is the frame in DN, truth the dict of `make_truth` (px, py, m606); missing: optional output of `unmeasured_f606w`."""
    x0s, y0s, c = window_counts(truth['px'], truth['py'], truth['m606'], shape, MAG_DEPTH, n)
    occ = footprint_cells(truth['px'], truth['py'], shape)
    bright = truth['m606'] < BRIGHT_LIMIT
    bx, by = truth['px'][bright], truth['py'][bright]
    if missing is not None:
        lost = missing['f814'] < MISSING_LIMIT
        lx, ly = missing['px'][lost], missing['py'][lost]
    else:
        lx, ly = np.zeros(0), np.zeros(0)

    def inside(px, py, x0, y0):                                       # a source within DN_MARGIN px of the box [x0 - 0.5, x0 + n - 0.5) x [y0 - 0.5, y0 + n - 0.5)
        return bool(np.any((px >= x0 - DN_MARGIN - 0.5) & (px < x0 + n + DN_MARGIN - 0.5) & (py >= y0 - DN_MARGIN - 0.5) & (py < y0 + n + DN_MARGIN - 0.5)))
    out = []
    for j, y0 in enumerate(y0s):
        for i, x0 in enumerate(x0s):
            d = c[j, i] / float(n * n)
            if not (CLASSES[0][0] <= d <= CLASSES[-1][1]):
                continue
            x0, y0 = int(x0), int(y0)
            dmax = float(dn_full[max(y0 - DN_MARGIN, 0): y0 + n + DN_MARGIN, max(x0 - DN_MARGIN, 0): x0 + n + DN_MARGIN].max())
            out.append(dict(x0=x0, y0=y0, count=int(c[j, i]), density=float(d), dn_max=dmax,
                            margin=not (x0 < FRAME_MARGIN or y0 < FRAME_MARGIN or x0 + n > shape[1] - FRAME_MARGIN or y0 + n > shape[0] - FRAME_MARGIN),
                            footprint=inside_footprint(occ, x0, y0, n), saturation=dmax < SAT_FRACTION * sat_level, bright=not inside(bx, by, x0, y0), lost=not inside(lx, ly, x0, y0)))
    return out


RULES = ('margin', 'footprint', 'saturation', 'bright', 'lost')


def candidate_patches(truth, dn_full, sat_level, shape=FRAME_SHAPE, n=N_PATCH, missing=None):
    """The origins of `screen_patches` that pass every rule of the module docstring: list of dict(x0, y0, count, density, dn_max)."""
    return [dict(x0=q['x0'], y0=q['y0'], count=q['count'], density=q['density'], dn_max=q['dn_max']) for q in screen_patches(truth, dn_full, sat_level, shape, n, missing)
            if all(q[r] for r in RULES)]


def funnel(screened):
    """How many origins of the density range survive the rules applied one after the other: [('density range', n), ('margin', n), ...]."""
    out, keep = [('density range 0.03-0.15', len(screened))], list(screened)
    for r in RULES:
        keep = [q for q in keep if q[r]]
        out.append((r, len(keep)))
    return out


def select_patches(cands, n=N_PATCH, classes=CLASSES, targets=TARGET_DENSITY):
    """One patch per density class: the candidate of the class whose density is closest to the class target and that does not overlap a patch chosen before; the classes are taken
    in order of the number of candidates they have (the scarcest first, ties by density) so that a class is not lost to the overlap of a more abundant one; ties between candidates by
    (y0, x0). A class without a usable candidate is skipped (reported by the caller). Returns the chosen candidates in order of increasing class density."""
    members = [[c for c in cands if lo <= c['density'] < hi or (hi == classes[-1][1] and c['density'] == hi)] for lo, hi in classes]
    chosen = {}
    for k in sorted(range(len(classes)), key=lambda k: (len(members[k]), k)):
        pool = [c for c in members[k] if all(abs(c['x0'] - o['x0']) >= n or abs(c['y0'] - o['y0']) >= n for o in chosen.values())]
        if pool:
            chosen[k] = min(pool, key=lambda c: (abs(c['density'] - targets[k]), c['y0'], c['x0']))
    return [chosen[k] for k in sorted(chosen)]


# ----------------------------------------------------------------------------- validation helpers
def robust_scatter(v):
    v = np.asarray(v, float)
    return float(1.4826 * np.median(np.abs(v - np.median(v)))) if v.size else float('nan')


def render_model(psf, x, y, flux, shape, pad=12):
    """Sum of flux_i x PSF(pixel - (x_i, y_i)) on a canvas of `shape` plus `pad` pixels on every side (canvas pixel (0, 0) is patch pixel (-pad, -pad)); `psf` is the 51 x 51 stamp
    whose centre pixel is the star centre, interpolated bicubically."""
    from scipy.interpolate import RectBivariateSpline
    r = psf.shape[0] // 2
    ax = np.arange(psf.shape[0]) - r
    spl = RectBivariateSpline(ax, ax, psf, kx=3, ky=3)
    ny, nx = shape[0] + 2 * pad, shape[1] + 2 * pad
    img = np.zeros((ny, nx))
    for xi, yi, fi in zip(x + pad, y + pad, flux):
        if xi < -r or yi < -r or xi > nx + r or yi > ny + r:
            continue
        k0, k1 = max(int(np.ceil(xi - r)), 0), min(int(np.floor(xi + r)), nx - 1)
        l0, l1 = max(int(np.ceil(yi - r)), 0), min(int(np.floor(yi + r)), ny - 1)
        if k1 < k0 or l1 < l0:
            continue
        img[l0: l1 + 1, k0: k1 + 1] += fi * spl(np.arange(l0, l1 + 1) - yi, np.arange(k0, k1 + 1) - xi)
    return img


def rigid_shift(dn, gain, dark_variance, psf, tx, ty, tmag, mag_faint=25.0, half_range=0.4, step=0.05, pad=12):
    """Image-based check of the sky-to-pixel mapping on one patch: the PSF-convolved model of the truth sources (relative flux 10^(-0.4 m), free scale, a free plane for the sky) is
    shifted rigidly over a grid and compared with the patch by weighted least squares (variance DN / gain + dark variance). Returns (sx, sy, chi2_min, chi2 at zero shift): the shift (px)
    that has to be ADDED to the truth positions to match the image, from a parabola through the chi^2 grid. A measurement of the net offset of the mapping (astrometry plus the proper
    motion between the two epochs), not applied to anything. tx, ty, tmag: the truth in patch pixel coordinates, including a margin of `pad` pixels around the patch."""
    sel = tmag < mag_faint
    canvas = render_model(psf, tx[sel], ty[sel], 10.0 ** (-0.4 * (tmag[sel] - ZP)), dn.shape, pad)
    sw = 1.0 / np.sqrt(dn_variance(dn, gain, dark_variance)).ravel()
    yy, xx = np.mgrid[0:dn.shape[0], 0:dn.shape[1]]
    plane = [np.ones(dn.size), (xx.ravel() - dn.shape[1] / 2.0) / dn.shape[1], (yy.ravel() - dn.shape[0] / 2.0) / dn.shape[0]]
    grid = np.arange(-half_range, half_range + 1e-9, step)
    chi = np.zeros((len(grid), len(grid)))
    for j, sy in enumerate(grid):
        for i, sx in enumerate(grid):
            mod = nd_shift(canvas, (sy, sx), order=3, mode='nearest')[pad: -pad, pad: -pad]
            A = np.vstack([mod.ravel()] + plane).T * sw[:, None]
            b = dn.ravel() * sw
            coef = np.linalg.lstsq(A, b, rcond=None)[0]
            chi[j, i] = float(((b - A @ coef) ** 2).sum())
    j0, i0 = np.unravel_index(np.argmin(chi), chi.shape)
    j0, i0 = min(max(j0, 1), len(grid) - 2), min(max(i0, 1), len(grid) - 2)

    def vertex(c_m, c_0, c_p):
        den = c_m - 2 * c_0 + c_p
        return 0.0 if den <= 0 else 0.5 * (c_m - c_p) / den
    sx = grid[i0] + step * vertex(chi[j0, i0 - 1], chi[j0, i0], chi[j0, i0 + 1])
    sy = grid[j0] + step * vertex(chi[j0 - 1, i0], chi[j0, i0], chi[j0 + 1, i0])
    return float(sx), float(sy), float(chi.min()), float(chi[len(grid) // 2, len(grid) // 2])


# ----------------------------------------------------------------------------- survey: which SDSS runs and fields cover a cluster
def run_candidates(runs, ra, dec, half_width=1.5, field_mu=1361 * PIX / 3600.0):
    """Runs of the photoRunAll table (dict of arrays RUN, NODE, INCL, MU_REF, FIELD_REF, STARTFIELD, ENDFIELD, STRIPE, STRIP) whose scan strip passes over (ra, dec): the great-circle
    latitude is within `half_width` degrees of the scan circle and the field coordinate along the run lies inside the run's fields. Returns a list of dict(run, strip, stripe, nu, field)."""
    mu, nu = eq_to_gc(ra, dec, runs['NODE'], runs['INCL'])
    dmu = (mu - runs['MU_REF'] + 180.0) % 360.0 - 180.0
    f = runs['FIELD_REF'] + dmu / field_mu
    ok = (np.abs(nu) < half_width) & (f > runs['STARTFIELD'] - 1) & (f < runs['ENDFIELD'] + 1)
    return [dict(run=int(runs['RUN'][i]), strip=str(runs['STRIP'][i]).strip(), stripe=int(runs['STRIPE'][i]), nu=float(nu[i]), field=float(f[i])) for i in np.flatnonzero(ok)]


def sources_in_field(ra, dec, node, incl, box):
    """Source positions in pixel units relative to a photoField box (dict MU_START, MU_END, NU_START, NU_END of one field, great-circle coordinates of the run; x from nu, y from
    mu, 0.396 arcsec pixels): returns (x, y, inside mask, width, height). The box is the 1361-row field proper; the frame carries 128 more rows beyond it."""
    mu, nu = eq_to_gc(ra, dec, node, incl)
    lo, hi = min(box['NU_START'], box['NU_END']), max(box['NU_START'], box['NU_END'])
    x = (nu - lo) * 3600.0 / PIX
    y = (((mu - box['MU_START'] + 180.0) % 360.0) - 180.0) * 3600.0 / PIX
    w, h = (hi - lo) * 3600.0 / PIX, ((box['MU_END'] - box['MU_START'] + 180.0) % 360.0 - 180.0) * 3600.0 / PIX
    return x, y, (x >= 0) & (x < w) & (y >= 0) & (y < h), w, h


# ----------------------------------------------------------------------------- command line
def _raw(cluster_dir, key, filename):
    return os.path.join(ROOT2, cluster_dir, key, filename)


def _m2_paths():
    return dict(frame=_raw('ngc7089', 'sdss_frame_r_002583_2_0136', 'frame-r-002583-2-0136.fits.bz2'), psfield=_raw('ngc7089', 'sdss_psfield_002583_2_0136', 'psField-002583-2-0136.fit'),
                catalog=_raw('ngc7089', 'hugs_ngc7089', 'hlsp_hugs_hst_wfc3-uvis-acs-wfc_ngc7089_multi_v1_catalog-meth1.txt'),
                photofield=_raw('ngc7089', 'sdss_photofield_002583_2', 'photoField-002583-2.fits'))


def _target_paths(cluster):
    run, cam, fld = TARGETS[cluster]
    return dict(frame=_raw(cluster, 'sdss_frame_r_%06d_%d_%04d' % (run, cam, fld), 'frame-r-%06d-%d-%04d.fits.bz2' % (run, cam, fld)),
                psfield=_raw(cluster, 'sdss_psfield_%06d_%d_%04d' % (run, cam, fld), 'psField-%06d-%d-%04d.fit' % (run, cam, fld)),
                catalog=_raw(cluster, 'hugs_' + cluster, 'hlsp_hugs_hst_wfc3-uvis-acs-wfc_%s_multi_v1_catalog-meth1.txt' % cluster),
                photofield=_raw('survey', 'sdss_photofield_%06d_%d' % (run, cam), 'photoField-%06d-%d.fits' % (run, cam)), run=run, camcol=cam, field=fld)


def _saturation_level(photofield_path, field):
    from astropy.io import fits
    t = fits.getdata(photofield_path, 1)
    return float(t['SATURATION_LEVEL'][t['FIELD'] == field][0][2])


def cmd_survey(args):
    from astropy.io import fits
    sums = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'astro', 'checksums_second.json')))
    centres = sums['cluster_centres']
    t = fits.getdata(_raw('survey', 'photorunall', 'photoRunAll-dr17.fits'), 1)
    runs = {k: np.asarray(t[k]) for k in ('RUN', 'NODE', 'INCL', 'MU_REF', 'FIELD_REF', 'STARTFIELD', 'ENDFIELD', 'STRIPE', 'STRIP')}
    print('photoRunAll-dr17: %d runs (rerun 301)' % len(t))
    print('%-8s %-9s %7s %8s %6s | %s' % ('cluster', 'name', 'RA', 'Dec', 'runs', 'best field: run camcol field | seeing"  footprint-fraction(F606W<22)  windows(0.03-0.15)  density-range'))
    rows = []
    for c, cen in centres.items():
        cand = run_candidates(runs, cen['ra'], cen['dec'])
        if c == 'ngc7089':
            print('%-8s %-9s %7.3f %+8.3f %6d | (the reference: the benchmark frame is run 2583 camcol 2 field 136; not part of the test-set selection)' % (c, cen['name'], cen['ra'], cen['dec'], len(cand)))
            continue
        if not cand:
            print('%-8s %-9s %7.3f %+8.3f %6d |' % (c, cen['name'], cen['ra'], cen['dec'], 0))
            continue
        cat_path = _raw(c, 'hugs_' + c, 'hlsp_hugs_hst_wfc3-uvis-acs-wfc_%s_multi_v1_catalog-meth1.txt' % c)
        cat = read_hugs(cat_path)
        ok = (cat['m606'] > -90) & (cat['m606'] < 90)
        ra, dec, m = cat['ra'][ok], cat['dec'][ok], cat['m606'][ok]
        best = None
        for rc in sorted({r['run'] for r in cand}):
            i = int(np.flatnonzero(runs['RUN'] == rc)[0])
            node, incl = float(runs['NODE'][i]), float(runs['INCL'][i])
            for cam in range(1, 7):
                pfp = _raw('survey', 'sdss_photofield_%06d_%d' % (rc, cam), 'photoField-%06d-%d.fits' % (rc, cam))
                if not os.path.exists(pfp):
                    continue
                pf = fits.getdata(pfp, 1)
                mu0, nu0 = eq_to_gc(cen['ra'], cen['dec'], node, incl)
                dm = (mu0 - pf['MU_START'] + 180) % 360 - 180
                near = np.flatnonzero((dm > -0.4) & (dm < 0.55) & (nu0 > np.minimum(pf['NU_START'], pf['NU_END']) - 0.5) & (nu0 < np.maximum(pf['NU_START'], pf['NU_END']) + 0.5))
                for k in near:
                    box = {q: float(pf[q][k]) for q in ('MU_START', 'MU_END', 'NU_START', 'NU_END')}
                    x, y, inside, w, h = sources_in_field(ra, dec, node, incl, box)
                    b = m < MAG_DEPTH
                    frac = float((inside & b).sum() / max(b.sum(), 1))
                    if frac < 0.02:
                        continue
                    sat = float(pf['SATURATION_LEVEL'][k][2])
                    tr = dict(px=x[inside], py=y[inside], m606=m[inside])
                    x0s, y0s, cnt = window_counts(tr['px'], tr['py'], tr['m606'], (int(h), int(w)))
                    occ = footprint_cells(tr['px'], tr['py'], (int(h), int(w)))
                    dens = []
                    for jj, yy in enumerate(y0s):
                        for ii, xx in enumerate(x0s):
                            d = cnt[jj, ii] / 1e4
                            if CLASSES[0][0] <= d <= CLASSES[-1][1] and xx >= 30 and yy >= 30 and xx + 100 <= w - 30 and yy + 100 <= h - 30 and inside_footprint(occ, int(xx), int(yy)):
                                dens.append(d)
                    rec = (frac, rc, cam, int(pf['FIELD'][k]), float(pf['PSF_WIDTH'][k][2]), len(dens), (min(dens), max(dens)) if dens else None, int(pf['PSP_STATUS'][k]), int(pf['PSF_NSTAR'][k][2]), sat)
                    if best is None or rec[0] > best[0] + 1e-9 or (abs(rec[0] - best[0]) < 1e-9 and rec[4] < best[4]):
                        best = rec
        if best is None:
            print('%-8s %-9s %7.3f %+8.3f %6d | (no field box of these runs holds >= 2%% of the footprint: the cluster falls in a camcol gap, or the photoField is not on disk)' % (c, cen['name'], cen['ra'], cen['dec'], len({r['run'] for r in cand})))
            continue
        frac, rc, cam, fld, see, nw, dr, ps, nst, sat = best
        print('%-8s %-9s %7.3f %+8.3f %6d | %5d %d %4d | %.3f  %.2f  %5d  %s | psp_status %d psf_nstar(r) %d sat %.0f' % (
            c, cen['name'], cen['ra'], cen['dec'], len({r['run'] for r in cand}), rc, cam, fld, see, frac, nw, ('%.2f-%.2f' % dr) if dr else '-', ps, nst, sat))
        rows.append((c, rc, cam, fld, see, frac, nw))
    return rows


def _fmt_stats(v, unit='px', fmt='%.4f'):
    v = np.asarray(v, float)
    return ('n %d  median ' + fmt + '  90th pct ' + fmt + '  max ' + fmt + ' %s') % (v.size, np.median(v), np.percentile(v, 90), v.max(), unit)


def world2pix_for(frame, mapping):
    """The sky-to-pixel mapping of a frame as a function (ra, dec) -> (x, y): 'header' (TAN of the frame header, the mapping of the M2 truth) or 'astrans' (the full solution, HDU 3)."""
    if mapping == 'header':
        p = tan_params(frame['header'])
        return lambda ra, dec: tan_world2pix(ra, dec, p)
    if mapping == 'astrans':
        return lambda ra, dec: astrans_world2pix(ra, dec, frame['astrans'])
    raise ValueError(mapping)


def compare_with_list(truth, lx, ly, lm, mag_lim, rmax=1.5, dmag=0.5):
    """Every source of a list (positions lx, ly in frame pixels, magnitudes lm) brighter than mag_lim against the nearest truth source within rmax px and dmag mag (the benchmark's
    magnitude tolerance): dict(n, n_matched, d (nearest distances), dx, dy (truth - list), dm)."""
    sel = lm < mag_lim
    rows = []
    for xx, yy, mm in zip(lx[sel], ly[sel], lm[sel]):
        near = np.flatnonzero((np.abs(truth['px'] - xx) < rmax) & (np.abs(truth['py'] - yy) < rmax) & (np.abs(truth['m606'] - mm) < dmag))
        if near.size:
            dd = np.hypot(truth['px'][near] - xx, truth['py'][near] - yy)
            k = near[np.argmin(dd)]
            if dd.min() <= rmax:
                rows.append((truth['px'][k] - xx, truth['py'][k] - yy, truth['m606'][k] - mm, dd.min()))
    a = np.array(rows).reshape(-1, 4)
    return dict(n=int(sel.sum()), n_matched=len(a), dx=a[:, 0], dy=a[:, 1], dm=a[:, 2], d=a[:, 3])


def cmd_validate_m2(args):
    """The validation on M2: the mappings against the M2 truth file and the independent shipped truth, the image-based check, the unit conversion against the shipped patches."""
    P = _m2_paths()
    fr, ps = read_frame(P['frame']), read_psfield(P['psfield'])
    cat = read_hugs(P['catalog'])
    ok = (cat['m606'] > -90) & (cat['m606'] < 90)
    th = make_truth(cat, world2pix_for(fr, 'header'))
    ta = make_truth(cat, world2pix_for(fr, 'astrans'))
    print('HUGS NGC 7089 rows %d, with a valid F606W %d' % (len(cat['m606']), len(th['m606'])))
    ref = np.load(os.path.join(ROOT, 'truth_hugs_meth1.npz'))
    print('[1] the TAN mapping of the header against truth_hugs_meth1.npz (data/astro/make_truth.py: the astropy WCS of the same header; 0-based pixel centres at integers)')
    d = np.hypot(th['px'] - ref['px'], th['py'] - ref['py'])
    b20 = th['m606'] < 20
    print('    all         : %s' % _fmt_stats(d, fmt='%.3g'))
    print('    F606W < 20  : %s' % _fmt_stats(d[b20], fmt='%.3g'))
    print('    magnitudes  : max |m606 - file m606| = %.3g ; e606 and nf606 identical to the file: %s' % (np.abs(th['m606'] - ref['m606']).max(), bool(np.array_equal(th['e606'], ref['e606']) and np.array_equal(th['nf606'], ref['nf606']))))
    print('[2] the asTrans solution of the frame (HDU 3, colour term off) against the TAN mapping of the header (the same catalogue rows)')
    da = np.hypot(ta['px'] - th['px'], ta['py'] - th['py'])
    box = (th['px'] >= 309.5) & (th['px'] < 409.5) & (th['py'] >= 629.5) & (th['py'] < 729.5)
    print('    benchmark patch box, F606W < 20 : %s' % _fmt_stats(da[box & b20]))
    print('    benchmark patch box, all        : %s' % _fmt_stats(da[box]))
    print('    benchmark patch box, signed median (asTrans - header): dx %+.4f dy %+.4f px' % (np.median(ta['px'][box] - th['px'][box]), np.median(ta['py'][box] - th['py'][box])))
    inframe = (th['px'] > 0) & (th['px'] < 2047) & (th['py'] > 0) & (th['py'] < 1488)
    print('    whole HUGS footprint in the frame: %s' % _fmt_stats(da[inframe]))
    for col in (-0.3, 0.3, 1.0):
        xc, yc = astrans_world2pix(cat['ra'][ok][box], cat['dec'][ok][box], fr['astrans'], color=col)
        print('    colour r-i = %+.1f against colour 0 (patch box): median shift (%+.4f, %+.4f) px' % (col, np.median(xc - ta['px'][box]), np.median(yc - ta['py'][box])))
    print('[3] the independent shipped truth (Data/sdss0921_true.txt of the public repository: ACS survey positions in the patch (310, 630), flux in SDSS counts) against both mappings')
    shipped = os.environ.get('ASTRO_SHIPPED', os.path.join(ROOT, 'shipped_inputs'))
    st = np.loadtxt(os.path.join(shipped, 'Data', 'sdss0921_true.txt'))
    sx, sy, sm = st[:, 0] + 310, st[:, 1] + 630, counts_to_mag(st[:, 2], 0.00546689)
    for name, t in (('header ', th), ('asTrans', ta)):
        for lim in (20, 22):
            c = compare_with_list(t, sx, sy, sm, lim)
            print('    %s, shipped sources brighter than %d: %d, matched (1.5 px, 0.5 mag) %d (%.1f%%); |HUGS - shipped|: median %.4f, 90th pct %.4f, max %.3f px ; median (dx, dy) = (%+.4f, %+.4f), scatter (%.4f, %.4f) ; dm median %+.4f, scatter %.4f, 90th pct |dm| %.4f' % (
                name, lim, c['n'], c['n_matched'], 100.0 * c['n_matched'] / c['n'], np.median(c['d']), np.percentile(c['d'], 90), c['d'].max(), np.median(c['dx']), np.median(c['dy']),
                robust_scatter(c['dx']), robust_scatter(c['dy']), np.median(c['dm']), robust_scatter(c['dm']), np.percentile(np.abs(c['dm']), 90)))
    print('[4] units: the benchmark patch (310, 630) rebuilt from the frame and the psField against the patch of fetch.py')
    arrays, meta = extract_patch(fr, ps, 310, 630)
    old = np.load(os.path.join(ROOT, 'patch_r_002583_2_0136_x310_y630_100x100', 'patch_arrays.npz'))
    print('    ' + '  '.join('%s %s' % (k, 'identical' if np.array_equal(np.asarray(arrays[k]), old[k]) else 'DIFFERS(max %.3g)' % np.abs(np.asarray(arrays[k], float) - old[k].astype(float)).max()) for k in old.files))
    print('    meta.json identical: %s' % (meta == json.load(open(os.path.join(ROOT, 'patch_r_002583_2_0136_x310_y630_100x100', 'meta.json')))))
    print('    PSF against pydl.photoop.image.sdss_psf_recon at the patch centre: max |difference| %s' % psf_vs_pydl(ps['kl'], 359.5, 679.5))
    import tempfile
    from astropy.io import fits
    with tempfile.TemporaryDirectory() as tmp:
        write_patch(tmp, arrays, meta, fr['header'])
        new_h, old_h = fits.getheader(os.path.join(tmp, 'patch_dn.fits')), fits.getheader(os.path.join(ROOT, 'patch_r_002583_2_0136_x310_y630_100x100', 'patch_dn.fits'))
        same_data = np.array_equal(fits.getdata(os.path.join(tmp, 'patch_dn.fits')), fits.getdata(os.path.join(ROOT, 'patch_r_002583_2_0136_x310_y630_100x100', 'patch_dn.fits')))
        zn, zo = np.load(os.path.join(tmp, 'patch_arrays.npz')), np.load(os.path.join(ROOT, 'patch_r_002583_2_0136_x310_y630_100x100', 'patch_arrays.npz'))
        print('    write_patch against the files of fetch.py: patch_dn.fits data identical %s, header cards: %d common with equal values %s, only in the new file %s, only in the old %s ; npz arrays and dtypes identical %s ; meta.json text identical %s' % (
            same_data, len([k for k in old_h if k in new_h]), all(old_h[k] == new_h[k] for k in old_h if k in new_h), [k for k in new_h if k not in old_h], [k for k in old_h if k not in new_h],
            zn.files == zo.files and all(np.array_equal(zn[k], zo[k]) and zn[k].dtype == zo[k].dtype for k in zo.files),
            open(os.path.join(tmp, 'meta.json')).read() == open(os.path.join(ROOT, 'patch_r_002583_2_0136_x310_y630_100x100', 'meta.json')).read()))
    cnt = np.loadtxt(os.path.join(shipped, 'Data', 'sdss0921_cntp.txt'))
    raw = cnt.reshape(100, 100) - 1044.0
    print('    shipped raw counts (bias 1044 removed) against dn: correlation %.5f, median difference (shipped - dn) %.2f DN' % (np.corrcoef(raw.ravel(), arrays['dn'].ravel())[0, 1], np.median(raw - arrays['dn'])))
    print('[5] image-based rigid shift of the HUGS positions on the benchmark patch (PSF-convolved model, free scale and sky plane); the shift to add to the positions')
    pad = 12
    for name, t in (('header ', th), ('asTrans', ta)):
        tp = truth_for_patch(t, 310 - pad, 630 - pad, N_PATCH + 2 * pad)
        sx_, sy_, chi_min, chi_0 = rigid_shift(arrays['dn'], ps['gain'], ps['dark_variance'], arrays['psf51'], tp['x'] - pad, tp['y'] - pad, tp['mag'])
        print('    %s: dx %+.3f dy %+.3f px ; chi2 %.0f at the best shift, %.0f at zero shift (%d pixels)' % (name, sx_, sy_, chi_min, chi_0, N_PATCH * N_PATCH))
    print('[6] the M2 patches by the same functions (truth: truth_hugs_meth1.npz); the last two were the first choice of test patches and were replaced after a referee round found them mostly outside the HUGS footprint')
    for (x0, y0, role) in ((110, 630, 'development'), (210, 630, 'development'), (310, 630, 'test (benchmark)'), (310, 730, 'test'), (360, 1040, 'test'), (100, 1100, 'test'), (410, 630, 'replaced'), (410, 730, 'replaced')):
        d_ = os.path.join(ROOT, 'patch_r_002583_2_0136_x%d_y%d_100x100' % (x0, y0), 'patch_arrays.npz')
        if not os.path.exists(d_):
            print('    (%4d, %4d) %-17s patch not on disk' % (x0, y0, role))
            continue
        z = np.load(d_)
        fw = psf_fwhm(z['psf51'])
        tr_ = dict(px=ref['px'], py=ref['py'], m606=ref['m606'])
        cr = crowding(tr_, x0, y0, fw)
        bc = block_counts(tr_, x0, y0)
        print('    (%4d, %4d) %-17s F606W<22 %4d (%.4f /px)  scored %4d  FWHM %.2f px  per PSF area %.2f  blend(2.22 px) %.2f  DNmax %6.0f  DNmedian %5.0f  sky %5.1f  peak S/N at 22 %.2f  blocks F606W<24: min %d median %.0f max %d, below 50: %d' % (
            x0, y0, role, cr['n'], cr['density'], cr['n_scored'], fw, cr['per_psf_area'], cr['blend_m2_fwhm'], z['dn'].max(), np.median(z['dn']), np.median(z['sky_dn']),
            peak_snr(22.0, z['psf51'].max(), float(fr['header']['NMGY']), float(np.median(z['sky_dn'])), float(z['gain']), float(z['dark_variance'])), bc.min(), np.median(bc), bc.max(), (bc < 50).sum()))


def cmd_build(args):
    c = args.cluster
    P = _target_paths(c)
    fr, ps = read_frame(P['frame']), read_psfield(P['psfield'])
    cat = read_hugs(P['catalog'])
    truth = make_truth(cat, world2pix_for(fr, 'astrans'))
    truth_h = make_truth(cat, world2pix_for(fr, 'header'))
    miss = unmeasured_f606w(cat, world2pix_for(fr, 'astrans'))
    print('frame run %d camcol %d field %d (rerun 301); NMGY %.8f; DATE-OBS %s; psField gain %.4f dark variance %.4f' % (P['run'], P['camcol'], P['field'], fr['header']['NMGY'], fr['header']['DATE-OBS'], ps['gain'], ps['dark_variance']))
    if not args.dry_run:
        for suffix, t in (('', truth), ('_header', truth_h)):
            np.savez(os.path.join(ROOT2, 'truth_hugs_%s_meth1%s.npz' % (c, suffix)), **t)
        print('wrote truth_hugs_%s_meth1.npz (asTrans mapping) and truth_hugs_%s_meth1_header.npz (header TAN mapping): %d sources with a valid F606W of %d rows' % (c, c, len(truth['m606']), len(cat['m606'])))
    sky = sky_image(fr['allsky'], fr['xinterp'], fr['yinterp'])
    dn_full = frame_to_dn(fr['nmgy'], fr['cal'], sky)
    sat = _saturation_level(P['photofield'], P['field'])
    screened = screen_patches(truth, dn_full, sat, missing=miss)
    cands = [dict(x0=q['x0'], y0=q['y0'], count=q['count'], density=q['density'], dn_max=q['dn_max']) for q in screened if all(q[r] for r in RULES)]
    print('saturation level (photoField, r) %.0f DN; origins on the 10 px grid after each rule: %s' % (sat, '; '.join('%s %d' % f for f in funnel(screened))))
    print('candidate origins satisfying all rules: %d' % len(cands))
    for (lo, hi) in CLASSES:
        print('   class [%.2f, %.2f]: %d candidates' % (lo, hi, sum(1 for q in cands if lo <= q['density'] <= hi)))
    chosen = select_patches(cands)
    if args.dry_run:                                                 # the selection only: nothing is written
        psf = unit_psf(sdss_psf_recon(ps['kl'], 1023.5, 744.0))
        print('dry run: PSF FWHM at the frame centre %.2f px; %d patches selected: %s' % (psf_fwhm(psf), len(chosen), [(q['x0'], q['y0'], round(q['density'], 4)) for q in chosen]))
        return []
    table = []
    pad = 12
    for q in chosen:
        x0, y0 = q['x0'], q['y0']
        arrays, meta = extract_patch(fr, ps, x0, y0, dn_full=dn_full, sky=sky)
        name = patch_name(P['run'], P['camcol'], P['field'], x0, y0)
        outdir = os.path.join(ROOT2, name)
        write_patch(outdir, arrays, meta, fr['header'])
        t, th_ = truth_for_patch(truth, x0, y0, N_PATCH), truth_for_patch(truth_h, x0, y0, N_PATCH)
        inbox = (miss['px'] >= x0 - 0.5) & (miss['px'] < x0 + N_PATCH - 0.5) & (miss['py'] >= y0 - 0.5) & (miss['py'] < y0 + N_PATCH - 0.5)
        files = ('patch_arrays.npz', 'patch_dn.fits', 'meta.json')
        row = dict(name=name, x0=x0, y0=y0, n22=int((t['mag'] < 22).sum()), density22=float((t['mag'] < 22).sum()) / 1e4, n21=int((t['mag'] < 21).sum()), n23=int((t['mag'] < 23).sum()),
                   n_all=int(len(t['mag'])), n22_scored=int(((t['mag'] < 22) & t['inner']).sum()), n22_header_mapping=int((th_['mag'] < 22).sum()), brightest=float(t['mag'].min()),
                   n_no_f606w_f814w_lt_20=int((inbox & (miss['f814'] < 20)).sum()), n_no_f606w_f814w_lt_22=int(inbox.sum()),
                   dn_max=float(arrays['dn'].max()), dn_median=float(np.median(arrays['dn'])), sky_dn_median=float(np.median(arrays['sky_dn'])), fwhm_px=psf_fwhm(arrays['psf51']), crowding=crowding(truth, x0, y0, psf_fwhm(arrays['psf51'])),
                   peak_snr_at_22=peak_snr(22.0, arrays['psf51'].max(), float(fr['header']['NMGY']), float(np.median(arrays['sky_dn'])), ps['gain'], ps['dark_variance']),
                   blocks_lt_24=dict(min=int(block_counts(truth, x0, y0).min()), median=float(np.median(block_counts(truth, x0, y0))), max=int(block_counts(truth, x0, y0).max()), below_50=int((block_counts(truth, x0, y0) < 50).sum())),
                   psf_vs_pydl=psf_vs_pydl(ps['kl'], x0 + (N_PATCH - 1) / 2.0, y0 + (N_PATCH - 1) / 2.0),
                   m2_equal_snr_mag=equal_snr_magnitude(M2_SNR22, arrays['psf51'].max(), float(fr['header']['NMGY']), float(np.median(arrays['sky_dn'])), ps['gain'], ps['dark_variance']),
                   sha256=dict(patch_arrays_content=content_sha256(os.path.join(outdir, files[0])), patch_dn_fits=file_sha256(os.path.join(outdir, files[1])), meta_json=file_sha256(os.path.join(outdir, files[2]))),
                   bytes={f: os.path.getsize(os.path.join(outdir, f)) for f in files})
        m = (truth['px'] >= x0 - 0.5) & (truth['px'] < x0 + N_PATCH - 0.5) & (truth['py'] >= y0 - 0.5) & (truth['py'] < y0 + N_PATCH - 0.5)
        mh = (truth_h['px'] >= x0 - 0.5) & (truth_h['px'] < x0 + N_PATCH - 0.5) & (truth_h['py'] >= y0 - 0.5) & (truth_h['py'] < y0 + N_PATCH - 0.5)
        both = np.flatnonzero(m & mh)
        dxy = np.c_[truth['px'][both] - truth_h['px'][both], truth['py'][both] - truth_h['py'][both]]
        dd = np.hypot(dxy[:, 0], dxy[:, 1])
        row['asTrans_minus_header_px'] = dict(median=float(np.median(dd)), p90=float(np.percentile(dd, 90)), max=float(dd.max()), signed_median_dx=float(np.median(dxy[:, 0])), signed_median_dy=float(np.median(dxy[:, 1])))
        shifts = {}
        for label, tt in (('asTrans', truth), ('header', truth_h)):
            tp = truth_for_patch(tt, x0 - pad, y0 - pad, N_PATCH + 2 * pad)
            sx, sy, chi_min, chi_0 = rigid_shift(arrays['dn'], ps['gain'], ps['dark_variance'], arrays['psf51'], tp['x'] - pad, tp['y'] - pad, tp['mag'])
            shifts[label] = dict(dx=sx, dy=sy, chi2_best=chi_min, chi2_zero=chi_0)
        row['image_shift_px'] = shifts
        # bright truth stars sit on bright pixels: the largest DN in the 3 x 3 pixels around every scored F606W < 18 truth source against the median DN of the patch
        br = (t['mag'] < 18) & t['inner']
        ix, iy = np.rint(t['x'][br]).astype(int), np.rint(t['y'][br]).astype(int)
        ratio = np.array([arrays['dn'][max(yy - 1, 0): yy + 2, max(xx - 1, 0): xx + 2].max() / np.median(arrays['dn']) for yy, xx in zip(iy, ix)])
        row['bright_on_bright'] = dict(n=int(br.sum()), median_ratio=float(np.median(ratio)) if br.any() else None, fraction_ratio_ge_2=float(np.mean(ratio >= 2.0)) if br.any() else None)
        table.append(row)
        bl = row['blocks_lt_24']
        print('patch %-44s x0 %4d y0 %4d  F606W<22: %4d (%.4f /px)  <21: %4d  <23: %4d  all: %5d  scored<22: %4d  brightest %.2f  DNmax %7.0f  FWHM %.2f px  per PSF area %.2f  blend(2.22 px) %.2f blend(FWHM) %.2f  peak S/N at 22: %.2f (same as M2 at 22.0: F606W %.2f)  20x20 blocks F606W<24: min %d median %.0f max %d' % (
            name, x0, y0, row['n22'], row['density22'], row['n21'], row['n23'], row['n_all'], row['n22_scored'], row['brightest'], row['dn_max'], row['fwhm_px'],
            row['crowding']['per_psf_area'], row['crowding']['blend_m2_fwhm'], row['crowding']['blend_own_fwhm'], row['peak_snr_at_22'], row['m2_equal_snr_mag'], bl['min'], bl['median'], bl['max']))
        a = row['asTrans_minus_header_px']
        print('      asTrans - header on the box: median %.4f p90 %.4f max %.4f px (signed median %+.3f, %+.3f) ; image shift asTrans (%+.3f, %+.3f) chi2 %.0f -> %.0f ; header (%+.3f, %+.3f) chi2 %.0f -> %.0f ; bright(<18) 3x3 max >= 2 x median: %.2f of %d' % (
            a['median'], a['p90'], a['max'], a['signed_median_dx'], a['signed_median_dy'], shifts['asTrans']['dx'], shifts['asTrans']['dy'], shifts['asTrans']['chi2_zero'], shifts['asTrans']['chi2_best'],
            shifts['header']['dx'], shifts['header']['dy'], shifts['header']['chi2_zero'], shifts['header']['chi2_best'], row['bright_on_bright']['fraction_ratio_ge_2'], row['bright_on_bright']['n']))
    with open(os.path.join(ROOT2, 'second_cluster_%s_patches.json' % c), 'w') as f:
        json.dump(dict(cluster=c, run=P['run'], camcol=P['camcol'], field=P['field'], nmgy=float(fr['header']['NMGY']), date_obs=fr['header']['DATE-OBS'], saturation_level_dn=sat, n_candidates=len(cands),
                       candidates_per_class=[sum(1 for q in cands if lo <= q['density'] <= hi) for lo, hi in CLASSES], funnel=funnel(screened), patches=table), f, indent=1)
    print('wrote', os.path.join(ROOT2, 'second_cluster_%s_patches.json' % c))
    return table


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    sub.add_parser('survey')
    sub.add_parser('validate-m2')
    p = sub.add_parser('build')
    p.add_argument('--cluster', choices=sorted(TARGETS), required=True)
    p.add_argument('--dry-run', action='store_true', help='print the candidate counts and the selection, write nothing')
    args = ap.parse_args()
    {'survey': cmd_survey, 'validate-m2': cmd_validate_m2, 'build': cmd_build}[args.cmd](args)


if __name__ == '__main__':
    main()
