"""Tests of experiments/astro_second.py: the sky-to-pixel mappings on made-up headers, the unit conversion of the frame to DN, the patch files, the truth selection and edge rules, the
density / crowding computation, the patch selection rules, and (skipped when the data are absent) the validation numbers on M2 and the patches of the second cluster that
results/astro/second_cluster_spec.md reports. Data tests read $ASTRO_DATA_2 (fetch_second.py), $ASTRO_DATA (the M2 files) and $ASTRO_SHIPPED (the shipped M2 inputs) and need astropy."""
import json
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'experiments'))
import astro_common as A      # noqa: E402
import astro_second as S      # noqa: E402

M2_PATCH = os.path.join(S.ROOT, 'patch_r_002583_2_0136_x310_y630_100x100')
SHIPPED = os.environ.get('ASTRO_SHIPPED', os.path.join(S.ROOT, 'shipped_inputs'))


def _header(ra0=150.0, dec0=30.0, theta=23.0, scale=0.396 / 3600.0):
    c, s = np.cos(np.deg2rad(theta)), np.sin(np.deg2rad(theta))
    return {'CRVAL1': ra0, 'CRVAL2': dec0, 'CRPIX1': 1025.0, 'CRPIX2': 745.0, 'CD1_1': scale * c, 'CD1_2': scale * s, 'CD2_1': scale * s, 'CD2_2': -scale * c}


def _astrans(**kw):
    t = dict(NODE=95.0, INCL=30.0, A=310.7, B=1.1e-4, C=2.8e-9, D=-0.2, E=1.1e-9, F=1.1e-4, DROW0=0.03, DROW1=-9e-5, DROW2=7e-8, DROW3=-1.3e-11, DCOL0=-0.16, DCOL1=7.8e-4, DCOL2=-8.5e-7,
             DCOL3=2.6e-10, CSROW=0.017, CSCOL=0.025, CCROW=1.7, CCCOL=2.5, RICUT=100.0)
    t.update(kw)
    return t


# ----------------------------------------------------------------------------- sky geometry
def test_tan_mapping_reference_pixel_axes_and_round_trip():
    p = S.tan_params(_header())
    x, y = S.tan_world2pix(150.0, 30.0, p)
    assert abs(x - 1024.0) < 1e-9 and abs(y - 744.0) < 1e-9                       # FITS pixel 1025 / 745 is 0-based 1024 / 744, the pixel centre at the integer
    ra, dec = S.tan_pix2world(1024.0, 744.0, p)
    assert abs(ra - 150.0) < 1e-9 and abs(dec - 30.0) < 1e-9
    # a source 10 arcsec north lands where CD^-1 (0, 10 arcsec) says; 10 arcsec east (xi = 10 arcsec cos(dec)) likewise
    cd = np.array(p['cd'])
    for dra, ddec in ((0.0, 10.0 / 3600.0), (10.0 / 3600.0 / np.cos(np.deg2rad(30.0)), 0.0)):
        x, y = S.tan_world2pix(150.0 + dra, 30.0 + ddec, p)
        exp = np.linalg.solve(cd, np.array([dra * np.cos(np.deg2rad(30.0)), ddec]))
        assert abs((x - 1024.0) - exp[0]) < 2e-3 and abs((y - 744.0) - exp[1]) < 2e-3
    # round trip over the frame, also across RA = 0
    for ra0 in (150.0, 359.99):
        p = S.tan_params(_header(ra0=ra0, dec0=-0.7))
        xs, ys = np.meshgrid(np.linspace(0, 2047, 7), np.linspace(0, 1488, 5))
        ra, dec = S.tan_pix2world(xs, ys, p)
        assert ra.min() >= 0 and ra.max() < 360
        x2, y2 = S.tan_world2pix(ra, dec, p)
        assert np.abs(x2 - xs).max() < 1e-7 and np.abs(y2 - ys).max() < 1e-7


def test_tan_mapping_matches_astropy_on_a_rotated_header():
    wcs = pytest.importorskip('astropy.wcs')
    fits = pytest.importorskip('astropy.io.fits')
    h = _header(ra0=259.25568, dec0=43.09338)
    hdr = fits.Header()
    for k, v in dict(CTYPE1='RA---TAN', CTYPE2='DEC--TAN', CUNIT1='deg', CUNIT2='deg', **h).items():
        hdr[k] = v
    w = wcs.WCS(hdr)
    rng = np.random.default_rng(1)
    x, y = rng.uniform(0, 2047, 200), rng.uniform(0, 1488, 200)
    ra, dec = w.all_pix2world(x, y, 0)
    xs, ys = S.tan_world2pix(ra, dec, S.tan_params(h))
    assert np.abs(xs - x).max() < 1e-6 and np.abs(ys - y).max() < 1e-6
    ra2, dec2 = S.tan_pix2world(x, y, S.tan_params(h))
    assert np.abs(((ra2 - ra) + 180) % 360 - 180).max() < 1e-9 and np.abs(dec2 - dec).max() < 1e-9


def test_astrans_round_trip_wrap_and_half_pixel_convention():
    rng = np.random.default_rng(2)
    x, y = rng.uniform(0, 2047, 300), rng.uniform(0, 1488, 300)
    for t in (_astrans(), _astrans(A=400.0), _astrans(INCL=0.0, NODE=48.4)):          # a > 360 is allowed by the data model
        ra, dec = S.astrans_pix2world(x, y, t)
        x2, y2 = S.astrans_world2pix(ra, dec, t)
        assert np.abs(x2 - x).max() < 1e-7 and np.abs(y2 - y).max() < 1e-7
        col = rng.uniform(-0.3, 1.2, 300)                                             # with the colour term
        ra, dec = S.astrans_pix2world(x, y, t, color=col)
        x2, y2 = S.astrans_world2pix(ra, dec, t, color=col)
        assert np.abs(x2 - x).max() < 1e-7 and np.abs(y2 - y).max() < 1e-7
    # the SDSS convention: pixel (row, col) has its centre at (row + 0.5, col + 0.5); the 0-based pixel (x, y) = (0, 0) is (row, col) = (0.5, 0.5)
    t = _astrans(DROW0=0, DROW1=0, DROW2=0, DROW3=0, DCOL0=0, DCOL1=0, DCOL2=0, DCOL3=0)
    ra, dec = S.astrans_pix2world(0.0, 0.0, t)
    mu, nu = t['A'] + t['B'] * 0.5 + t['C'] * 0.5, t['D'] + t['E'] * 0.5 + t['F'] * 0.5
    ra_e, dec_e = S.gc_to_eq(mu, nu, t['NODE'], t['INCL'])
    assert abs(ra - ra_e) < 1e-12 and abs(dec - dec_e) < 1e-12
    # (mu, nu) <-> (RA, Dec) are inverse and the scan circle's own latitude is nu = 0
    r, d = S.gc_to_eq(123.0, 0.0, 95.0, 30.0)
    mu, nu = S.eq_to_gc(r, d, 95.0, 30.0)
    assert abs(mu - 123.0) < 1e-9 and abs(nu) < 1e-9


def test_astrans_agrees_with_the_tan_mapping_for_an_equatorial_scan_without_distortion():
    """An equatorial scan (node 0, inclination 0) has mu = RA, nu = Dec: with the rows along RA and the columns along Dec, as in the SDSS frames, the full solution and the header TAN
    solution built for the same field centre must agree to a few thousandths of a pixel over +-1000 px (their difference is of third order in the field size)."""
    ra0, dec0, s = 100.0, 0.0, 0.396 / 3600.0
    t = _astrans(NODE=0.0, INCL=0.0, B=s, C=0.0, E=0.0, F=s, DROW0=0, DROW1=0, DROW2=0, DROW3=0, DCOL0=0, DCOL1=0, DCOL2=0, DCOL3=0)
    t['A'] = ra0 - s * 744.5                                  # the field centre is the SDSS (row, col) = (744.5, 1024.5), i.e. FITS pixel (1025, 745), 0-based (1024, 744)
    t['D'] = dec0 - s * 1024.5
    p = S.tan_params({'CRVAL1': ra0, 'CRVAL2': dec0, 'CRPIX1': 1025.0, 'CRPIX2': 745.0, 'CD1_1': 0.0, 'CD1_2': s, 'CD2_1': s, 'CD2_2': 0.0})
    rng = np.random.default_rng(3)
    x, y = rng.uniform(24, 2024, 200), rng.uniform(0, 1488, 200)
    ra, dec = S.tan_pix2world(x, y, p)
    xa, ya = S.astrans_world2pix(ra, dec, t)
    assert np.abs(xa - x).max() < 5e-3 and np.abs(ya - y).max() < 5e-3


def test_astrans_distortion_and_colour_terms_act_as_the_data_model_says():
    s = 0.396 / 3600.0
    z = dict(DROW0=0, DROW1=0, DROW2=0, DROW3=0, DCOL0=0, DCOL1=0, DCOL2=0, DCOL3=0)
    base = _astrans(NODE=0.0, INCL=0.0, B=s, C=0.0, E=0.0, F=s, A=100.0, D=0.0, CSROW=0.02, CSCOL=0.03, CCROW=1.5, CCCOL=2.5, RICUT=0.5, **z)
    r0, d0 = S.astrans_pix2world(500.0, 700.0, base)
    # blue object (r - i < riCut): rowm += csRow x colour, colm += csCol x colour; mu = a + b rowm, nu = d + f colm
    r1, d1 = S.astrans_pix2world(500.0, 700.0, base, color=0.2)
    assert abs((r1 - r0) - s * 0.02 * 0.2) < 1e-12 and abs((d1 - d0) - s * 0.03 * 0.2) < 1e-12
    # red object (r - i >= riCut): the constants ccRow, ccCol
    r2, d2 = S.astrans_pix2world(500.0, 700.0, base, color=0.8)
    assert abs((r2 - r0) - s * 1.5) < 1e-12 and abs((d2 - d0) - s * 2.5) < 1e-12
    # the distortion polynomials are in the column (the data model): rowm = row + dRow0 + dRow1 col + ..., colm = col + dCol0 + ...
    t = dict(base, DROW0=0.1, DROW1=1e-3, DROW2=1e-6, DROW3=1e-9, DCOL0=-0.2, DCOL1=2e-3, DCOL2=-1e-6, DCOL3=1e-9)
    col = 500.5
    r3, d3 = S.astrans_pix2world(500.0, 700.0, t)
    rowm = 700.5 + 0.1 + 1e-3 * col + 1e-6 * col ** 2 + 1e-9 * col ** 3
    colm = col - 0.2 + 2e-3 * col - 1e-6 * col ** 2 + 1e-9 * col ** 3
    assert abs(r3 - (100.0 + s * rowm)) < 1e-10 and abs(d3 - s * colm) < 1e-10


# ----------------------------------------------------------------------------- units
def test_sky_interpolation_dn_conversion_and_variance_invert_the_data_model():
    iy, ix = np.mgrid[0:10, 0:12]
    allsky = 100.0 + 2.0 * ix + 3.0 * iy                                            # a plane: bilinear interpolation reproduces it exactly
    xinterp, yinterp = np.linspace(0, 11, 60), np.linspace(0, 12, 40)                # y goes beyond the table: constant extrapolation (the data model)
    sky = S.sky_image(allsky, xinterp, yinterp)
    assert sky.shape == (40, 60)
    exp = 100.0 + 2.0 * xinterp[None, :] + 3.0 * np.minimum(yinterp, 9.0)[:, None]
    assert np.allclose(sky, exp, atol=1e-9)
    rng = np.random.default_rng(4)
    nmgy = rng.normal(0, 3, (40, 60))
    cal = rng.uniform(0.004, 0.006, 60)
    dn = S.frame_to_dn(nmgy, cal, sky)
    assert np.allclose((dn - sky) * cal[None, :], nmgy, rtol=1e-12, atol=1e-12)    # dn = image / calibration + sky
    assert np.allclose(S.dn_variance(dn, 4.6, 1.0), dn / 4.6 + 1.0)                 # electrons = dn x gain are Poisson, plus the dark variance
    f = np.array([100.0, 1e4, 3e5])
    assert np.allclose(S.counts_to_mag(f, 0.00546689), A.counts_to_mag(f), rtol=0, atol=1e-12)
    assert abs(S.counts_to_mag(1.0 / 0.00546689, 0.00546689) - 22.5) < 1e-12 and np.isnan(S.counts_to_mag(-1.0, 0.0055))


def _kl(rn=11, nb=3, seed=0):
    rng = np.random.default_rng(seed)
    dt = np.dtype([('RNROW', 'i4'), ('RNCOL', 'i4'), ('nrow_b', 'i4'), ('ncol_b', 'i4'), ('RROWS', 'f4', (rn * rn,)), ('c', 'f4', (nb, nb))])
    kl = np.zeros(4, dt)
    kl['RNROW'], kl['RNCOL'], kl['nrow_b'], kl['ncol_b'] = rn, rn, nb, nb
    kl['RROWS'] = rng.normal(0, 1, (4, rn * rn))
    kl['c'] = rng.normal(0, 1, (4, nb, nb))
    return kl


def test_psf_reconstruction_is_the_polynomial_sum_of_the_eigen_images():
    kl = _kl()
    x, y = 359.5, 679.5
    psf = S.sdss_psf_recon(kl, x, y)
    assert psf.shape == (11, 11)
    # sum_i RROWS_i x sum_jk c_i[j, k] (0.001 (x + 0.5))^j (0.001 (y + 0.5))^k
    exp = np.zeros(121)
    for i in range(4):
        w = sum(float(kl['c'][i][j, k]) * (0.001 * (x + 0.5)) ** j * (0.001 * (y + 0.5)) ** k for j in range(3) for k in range(3))
        exp += kl['RROWS'][i].astype(float) * w
    assert np.allclose(psf.ravel(), exp, rtol=1e-12, atol=1e-12)
    u = S.unit_psf(psf)
    assert u.min() >= 0 and abs(u.sum() - 1) < 1e-12


def test_peak_signal_to_noise_of_a_point_source():
    counts = 10.0 ** (0.4 * 0.5) / 0.00546689                                       # a 22nd magnitude source
    peak = 0.1422 * counts
    assert abs(S.peak_snr(22.0, 0.1422, 0.00546689, 149.2, 4.6, 1.0) - peak / np.sqrt((149.2 + peak) / 4.6 + 1.0)) < 1e-12
    assert S.peak_snr(21.0, 0.1, 0.0055, 100.0, 4.6, 1.0) > S.peak_snr(22.0, 0.1, 0.0055, 100.0, 4.6, 1.0) > S.peak_snr(22.0, 0.07, 0.0055, 100.0, 4.6, 1.0)       # brighter and sharper is better


def test_equal_snr_magnitude_inverts_peak_snr_and_is_fainter_for_a_deeper_image():
    m = S.equal_snr_magnitude(6.33, 0.0752, 0.00534669, 109.5, 4.895, 0.9025)
    assert abs(S.peak_snr(m, 0.0752, 0.00534669, 109.5, 4.895, 0.9025) - 6.33) < 1e-9 and 21.0 < m < 22.0       # a shallower image reaches M2's S/N at 22 only at a brighter magnitude
    assert S.equal_snr_magnitude(6.33, 0.1422, 0.00546689, 149.2, 4.6, 1.0) > 21.9


def test_psf_fwhm_of_a_gaussian_stamp():
    x = np.arange(51) - 25
    g = np.exp(-0.5 * (x[None, :] ** 2 + x[:, None] ** 2) / 1.5 ** 2)
    assert abs(S.psf_fwhm(g) - 2.3548 * 1.5) < 0.03 and abs(S.psf_fwhm(g[:, ::-1].T) - S.psf_fwhm(g)) < 1e-9


def _synthetic_frame(rows=160, cols=150):
    iy, ix = np.mgrid[0:rows, 0:cols]
    cal = 0.005 * (1.0 + 0.01 * np.arange(cols))
    return dict(nmgy=iy * 1000.0 + ix, cal=cal, allsky=100.0 + 0.5 * np.mgrid[0:12, 0:14][1] + np.mgrid[0:12, 0:14][0], xinterp=np.linspace(0, 13, cols), yinterp=np.linspace(0, 11, rows),
                header={'RUN': 4682, 'RERUN': '301', 'CAMCOL': 6, 'FRAME': 208, 'FILTER': 'r', 'NMGY': 0.00545525})


def test_extract_patch_has_the_m2_keys_shapes_and_row_column_slices():
    fr = _synthetic_frame()
    ps = dict(kl=_kl(rn=51), gain=4.895, dark_variance=0.9025, sky=100.0, skysig=6.0)
    arrays, meta = S.extract_patch(fr, ps, 20, 30)
    assert set(arrays) == {'dn', 'electrons', 'nmgy', 'sky_dn', 'psf51', 'psf25', 'gain', 'dark_variance', 'x0', 'y0'}          # the keys of patch_arrays.npz of the M2 patches
    assert set(meta) == {'x0', 'y0', 'size', 'gain', 'dark_variance', 'psfield_sky', 'psfield_skysig', 'note'}
    assert arrays['dn'].shape == (100, 100) and arrays['psf51'].shape == (51, 51) and arrays['psf25'].shape == (25, 25)
    sky = S.sky_image(fr['allsky'], fr['xinterp'], fr['yinterp'])
    for (r, c) in ((0, 0), (5, 7), (99, 99)):                                       # array index [row, column]; the frame pixel is (row y0 + r, column x0 + c)
        assert abs(arrays['dn'][r, c] - ((30 + r) * 1000.0 + 20 + c) / fr['cal'][20 + c] - sky[30 + r, 20 + c]) < 1e-9
        assert arrays['nmgy'][r, c] == (30 + r) * 1000.0 + 20 + c and arrays['sky_dn'][r, c] == sky[30 + r, 20 + c]
    assert np.allclose(arrays['electrons'], arrays['dn'] * 4.895)
    assert abs(arrays['psf51'].sum() - 1) < 1e-12 and abs(arrays['psf25'].sum() - 1) < 1e-12 and arrays['x0'] == 20 and arrays['y0'] == 30
    assert S.patch_name(4682, 6, 208, 1080, 220) == 'patch_r_004682_6_0208_x1080_y220_100x100' and S.patch_name(2583, 2, 136, 310, 630) == 'patch_r_002583_2_0136_x310_y630_100x100'


def test_write_patch_files_follow_the_m2_layout_and_the_content_hash_is_reproducible(tmp_path):
    fits = pytest.importorskip('astropy.io.fits')
    pytest.importorskip('astropy.wcs')
    fr = _synthetic_frame()
    hdr = fits.Header()
    for k, v in dict(CTYPE1='RA---TAN', CTYPE2='DEC--TAN', CUNIT1='deg', CUNIT2='deg', **_header()).items():
        hdr[k] = v
    for k, v in fr['header'].items():
        hdr[k] = v
    ps = dict(kl=_kl(rn=51), gain=4.895, dark_variance=0.9025, sky=100.0, skysig=6.0)
    arrays, meta = S.extract_patch(fr, ps, 20, 30)
    S.write_patch(str(tmp_path / 'a'), arrays, meta, hdr)
    S.write_patch(str(tmp_path / 'b'), arrays, meta, hdr)
    assert sorted(os.listdir(tmp_path / 'a')) == ['meta.json', 'patch_arrays.npz', 'patch_dn.fits']
    z = np.load(tmp_path / 'a' / 'patch_arrays.npz')
    assert set(z.files) == set(arrays) and np.array_equal(z['dn'], arrays['dn']) and z['x0'].dtype == np.int64
    h = fits.getheader(tmp_path / 'a' / 'patch_dn.fits')
    assert h['PX0'] == 20 and h['PY0'] == 30 and h['BUNIT'] == 'DN' and h['RUN'] == 4682 and h['FRAME'] == 208 and abs(h['GAINR'] - 4.895) < 1e-12 and abs(h['NMGY'] - 0.00545525) < 1e-12
    assert abs(h['CRPIX1'] - (1025.0 - 20)) < 1e-9 and abs(h['CRPIX2'] - (745.0 - 30)) < 1e-9          # the frame WCS moved to the patch origin
    assert fits.getdata(tmp_path / 'a' / 'patch_dn.fits').dtype == np.dtype('>f4')
    assert json.load(open(tmp_path / 'a' / 'meta.json')) == meta
    assert S.content_sha256(str(tmp_path / 'a' / 'patch_arrays.npz')) == S.content_sha256(str(tmp_path / 'b' / 'patch_arrays.npz'))
    assert S.file_sha256(str(tmp_path / 'a' / 'patch_dn.fits')) == S.file_sha256(str(tmp_path / 'b' / 'patch_dn.fits'))


# ----------------------------------------------------------------------------- the truth list
def test_hugs_reader_and_truth_list_have_the_format_of_the_m2_truth(tmp_path):
    rows = []
    for k, (m, e, nf, ra, dec) in enumerate([(18.5, 0.01, 5, 10.0, 20.0), (-99.9999, 99.9999, 0, 10.001, 20.001), (22.2, 0.1, 99, 10.002, 20.002), (99.0, 0.0, 0, 10.003, 20.003)]):
        f = ['0.0000'] * 37
        f[20], f[21], f[24], f[33], f[34], f[35], f[36], f[26] = '%.4f' % m, '%.4f' % e, '%d' % nf, '%.8f' % ra, '%.8f' % dec, 'R%07d' % (k + 1), '1', '%.4f' % (19.0 + k)
        rows.append('   ' + '  '.join(f))
    path = tmp_path / 'cat.txt'
    path.write_text('# header line\n# Cols. 1,2: X, Y\n' + '\n'.join(rows) + '\n')
    cat = S.read_hugs(str(path))
    assert cat['m606'].tolist() == [18.5, -99.9999, 22.2, 99.0] and cat['nf606'].tolist() == [5, 0, 99, 0]
    assert np.allclose(cat['ra'], [10.0, 10.001, 10.002, 10.003]) and np.allclose(cat['dec'], [20.0, 20.001, 20.002, 20.003]) and cat['f814'].tolist() == [19.0, 20.0, 21.0, 22.0]
    miss = S.unmeasured_f606w(cat, lambda ra, dec: ((ra - 10.0) * 1000.0, (dec - 20.0) * 1000.0), 22.0)       # F606W not measured (-99.9999), F814W measured and brighter than 22: row 1 only
    assert miss['f814'].tolist() == [20.0] and np.allclose(miss['px'], [1.0])
    tr = S.make_truth(cat, lambda ra, dec: ((ra - 10.0) * 1000.0, (dec - 20.0) * 1000.0))
    assert set(tr) == {'px', 'py', 'm606', 'e606', 'nf606'}                         # the arrays of truth_hugs_meth1.npz
    assert tr['m606'].tolist() == [18.5, 22.2] and np.allclose(tr['px'], [0.0, 2.0]) and np.allclose(tr['py'], [0.0, 2.0])      # -99.9999 (not measured) and 99 are dropped
    assert all(v.dtype == np.float64 for v in tr.values())


def test_truth_for_patch_follows_the_rules_of_astro_common():
    rng = np.random.default_rng(5)
    px, py = rng.uniform(300, 420, 3000), rng.uniform(620, 740, 3000)
    # sources exactly on the half-pixel box edges and on the scoring margin
    px[:6] = [309.5, 409.5, 309.4999, 409.4999, 312.0, 407.0]
    py[:6] = [680.0] * 6
    m = rng.uniform(14, 28, 3000)
    truth = dict(px=px, py=py, m606=m)
    for mag_max in (30.0, 22.5, 18.0):
        a, b = S.truth_for_patch(truth, 310, 630, 100, mag_max), A.truth_for_patch(truth, 310, 630, 100, mag_max)
        assert set(a) == set(b)
        for k in a:
            assert np.array_equal(a[k], b[k])
    t = S.truth_for_patch(truth, 310, 630, 100)
    assert 309.5 - 310 in t['x'] and (409.5 - 310) not in t['x'] and (409.4999 - 310) in t['x']       # half-open box [x0 - 0.5, x0 + n - 0.5)
    inner = {round(float(x), 4): bool(i) for x, i, y in zip(t['x'], t['inner'], t['y']) if abs(y - 50) < 1e-9}
    assert inner[2.0] is True and inner[97.0] is True and inner[-0.5] is False and inner[99.4999] is False


# ----------------------------------------------------------------------------- crowding
def test_density_windows_blend_fraction_and_footprint():
    rng = np.random.default_rng(6)
    shape = (300, 400)
    px, py = rng.uniform(-0.5, 399.5, 6000), rng.uniform(-0.5, 299.5, 6000)
    m = rng.uniform(15, 25, 6000)
    x0s, y0s, c = S.window_counts(px, py, m, shape, 22.0, 100, 10)
    assert c.shape == (len(y0s), len(x0s)) and x0s[-1] == 300 and y0s[-1] == 200
    for (x0, y0) in ((0, 0), (120, 40), (300, 200), (50, 170)):
        brute = int(((px >= x0 - 0.5) & (px < x0 + 99.5) & (py >= y0 - 0.5) & (py < y0 + 99.5) & (m < 22.0)).sum())
        assert c[y0 // 10, x0 // 10] == brute
        assert abs(S.patch_density(dict(px=px, py=py, m606=m), x0, y0) - brute / 1e4) < 1e-15
    # blends: two of three sources have a neighbour closer than 2 px
    assert abs(S.blend_fraction([0, 1, 10], [0, 0, 10], 2.0) - 2.0 / 3.0) < 1e-12 and S.blend_fraction([0.0], [0.0], 2.0) == 0.0
    # the footprint: sources fill x in [100, 300), y in [60, 260); a patch needs its box plus a 20 px collar inside
    fx, fy = rng.uniform(100, 300, 20000), rng.uniform(60, 260, 20000)
    occ = S.footprint_cells(fx, fy, shape)
    assert S.inside_footprint(occ, 140, 100) and not S.inside_footprint(occ, 100, 100) and not S.inside_footprint(occ, 140, 60) and not S.inside_footprint(occ, 220, 100)
    cr = S.crowding(dict(px=px, py=py, m606=m), 120, 40, 3.0)
    assert cr['n'] == c[4, 12] and abs(cr['density'] - c[4, 12] / 1e4) < 1e-15 and abs(cr['per_psf_area'] - cr['density'] * np.pi * 1.5 ** 2) < 1e-12 and 0.0 <= cr['blend_own_fwhm'] <= 1.0


def test_block_counts_find_a_hole_of_the_truth():
    rng = np.random.default_rng(9)
    px, py = rng.uniform(0, 200, 20000), rng.uniform(0, 200, 20000)
    keep = ~((px > 100) & (py < 100))                                                # the lower right quarter of the box (x in [100, 200), y in [0, 100)) is outside the footprint
    truth = dict(px=px[keep], py=py[keep], m606=rng.uniform(16, 26, keep.sum()))
    h = S.block_counts(truth, 50, 0)                                                 # the box x in [49.5, 149.5): its right half (x >= 100) is empty for y < 100
    assert h.shape == (5, 5) and h.sum() == len(S.truth_for_patch(truth, 50, 0, 100, 24.0)['x'])
    assert (h[:, 3:] == 0).all() and h[:, :2].min() > 0                               # the two right-hand columns of blocks are empty, the two left-hand ones are not
    full = S.block_counts(truth, 0, 100)                                             # a box inside the covered area: no empty block, similar counts
    assert full.min() > 0.5 * np.median(full)


def _synthetic_field(seed=7, shape=(500, 700)):
    """Truth whose density rises with x (0.02 to 0.20 per pixel brighter than 22 inside the footprint x in [60, 640), y in [60, 440))."""
    rng = np.random.default_rng(seed)
    n = 60000
    x, y = rng.uniform(60, 640, n), rng.uniform(60, 440, n)
    keep = rng.uniform(0, 1, n) < (x - 60) / 580.0
    x, y = x[keep], y[keep]
    m = np.where(rng.uniform(0, 1, len(x)) < 0.4, rng.uniform(17, 22, len(x)), rng.uniform(22, 26, len(x)))
    return dict(px=x, py=y, m606=m, nf606=np.full(len(x), 5.0)), shape


def test_candidates_obey_the_rules_and_selection_is_deterministic_one_per_class_and_non_overlapping():
    truth, shape = _synthetic_field()
    dn = np.full(shape, 100.0)
    dn[300, 300] = 60000.0                                                           # a saturated pixel
    truth['px'], truth['py'], truth['m606'], truth['nf606'] = (np.r_[v, w] for v, w in ((truth['px'], [500.0]), (truth['py'], [150.0]), (truth['m606'], [13.0]), (truth['nf606'], [99.0])))
    cands = S.candidate_patches(truth, dn, 52900.0, shape)
    assert len(cands) > 20
    hole = dict(px=np.array([350.0, 150.0]), py=np.array([250.0, 380.0]), f814=np.array([16.0, 18.5]))        # a star lost in F606W: the bright one excludes patches, the faint one does not
    c_hole = S.candidate_patches(truth, dn, 52900.0, shape, missing=hole)
    assert len(c_hole) < len(cands) and all(not (q['x0'] - 10 - 0.5 <= 350.0 < q['x0'] + 110 - 0.5 and q['y0'] - 10 - 0.5 <= 250.0 < q['y0'] + 110 - 0.5) for q in c_hole)
    assert len(S.candidate_patches(truth, dn, 52900.0, shape, missing=dict(hole, f814=np.array([18.0, 18.5])))) == len(cands)
    occ = S.footprint_cells(truth['px'], truth['py'], shape)
    for q in cands:
        x0, y0 = q['x0'], q['y0']
        assert S.CLASSES[0][0] <= q['density'] <= S.CLASSES[-1][1] and abs(q['density'] - S.patch_density(truth, x0, y0)) < 1e-12
        assert x0 >= S.FRAME_MARGIN and y0 >= S.FRAME_MARGIN and x0 + 100 <= shape[1] - S.FRAME_MARGIN and y0 + 100 <= shape[0] - S.FRAME_MARGIN
        assert S.inside_footprint(occ, x0, y0)
        assert not (x0 - 10 <= 300 < x0 + 110 and y0 - 10 <= 300 < y0 + 110)         # the saturated pixel (within 10 px of the box)
        assert not (x0 - 10 - 0.5 <= 500.0 < x0 + 110 - 0.5 and y0 - 10 - 0.5 <= 150.0 < y0 + 110 - 0.5)       # the 13th magnitude star
    # without the guards the same field offers strictly more candidates
    dn2 = np.full(shape, 100.0)
    t2 = dict(truth, m606=np.where(truth['m606'] < 14.5, 20.0, truth['m606']))
    assert len(S.candidate_patches(t2, dn2, 52900.0, shape)) > len(cands)
    chosen = S.select_patches(cands)
    assert chosen == S.select_patches(list(reversed(cands)))                         # deterministic, independent of the order of the list
    assert 3 <= len(chosen) <= 4
    members = [[c for c in cands if lo <= c['density'] < hi or (k == len(S.CLASSES) - 1 and c['density'] == hi)] for k, (lo, hi) in enumerate(S.CLASSES)]
    assert [c['density'] for c in chosen] == sorted(c['density'] for c in chosen)
    classes_of = [next(k for k, m in enumerate(members) if c in m) for c in chosen]
    assert len(set(classes_of)) == len(chosen)                                       # one per class
    for a in range(len(chosen)):
        for b in range(a):
            assert abs(chosen[a]['x0'] - chosen[b]['x0']) >= 100 or abs(chosen[a]['y0'] - chosen[b]['y0']) >= 100        # no overlap
    k_scarce = min((k for k in range(len(S.CLASSES)) if members[k]), key=lambda k: (len(members[k]), k))        # the scarcest class is taken first: its choice is the candidate closest to the target
    q = chosen[classes_of.index(k_scarce)]
    assert all(abs(c['density'] - S.TARGET_DENSITY[k_scarce]) >= abs(q['density'] - S.TARGET_DENSITY[k_scarce]) - 1e-15 for c in members[k_scarce])
    # a class without candidates is skipped, not filled from another class; overlapping candidates of two classes are not both taken, and the scarcer class wins the overlap
    mk = lambda x0, y0, d: dict(x0=x0, y0=y0, count=int(d * 1e4), density=d, dn_max=100.0)
    got = S.select_patches([mk(0, 0, 0.09), mk(300, 0, 0.13), mk(310, 50, 0.131)])
    assert [(g['x0'], g['y0']) for g in got] == [(0, 0), (300, 0)]
    got = S.select_patches([mk(0, 0, 0.045), mk(50, 0, 0.0451), mk(400, 0, 0.12)])
    assert [(g['x0'], g['y0']) for g in got] == [(0, 0), (400, 0)]                       # the closest to the 0.04 target; the second candidate overlaps it and is not a class member anyway
    got = S.select_patches([mk(0, 0, 0.095), mk(40, 40, 0.12), mk(200, 0, 0.094), mk(200, 200, 0.0965)])
    assert [(g['x0'], g['y0']) for g in got] == [(200, 0), (40, 40)]                     # the only dense candidate is kept; the 0.095 candidate it overlaps gives way to the next one
    assert S.select_patches([]) == []


def test_run_candidates_use_the_great_circle_latitude_and_the_field_range():
    runs = dict(RUN=np.array([100, 200, 300]), NODE=np.array([95.0, 95.0, 95.0]), INCL=np.array([0.0, 30.0, 0.0]), MU_REF=np.array([145.0, 145.0, 145.0]),
                FIELD_REF=np.array([11, 11, 11]), STARTFIELD=np.array([11, 11, 11]), ENDFIELD=np.array([100, 100, 30]), STRIPE=np.array([82, 82, 82]), STRIP=np.array(['S', 'N', 'S']))
    got = S.run_candidates(runs, 150.0, -0.8)
    assert [g['run'] for g in got] == [100]                                          # run 200 is inclined (its circle is far from the point), run 300 ends before the field
    assert abs(got[0]['nu'] + 0.8) < 1e-9 and abs(got[0]['field'] - (11 + 5.0 / (1361 * 0.396 / 3600.0))) < 1e-6
    assert S.run_candidates(runs, 150.0, 3.0) == []                                  # 3 degrees from the scan circle


# ----------------------------------------------------------------------------- the image-based check
def test_rigid_shift_recovers_a_known_shift_and_zero_for_aligned_truth():
    rng = np.random.default_rng(8)
    g = np.arange(51) - 25
    psf = np.exp(-0.5 * (g[None, :] ** 2 + g[:, None] ** 2) / 1.4 ** 2)
    psf /= psf.sum()
    n_src, pad = 60, 12
    tx, ty = rng.uniform(-pad, 100 + pad, n_src), rng.uniform(-pad, 100 + pad, n_src)
    tm = rng.uniform(17.0, 22.0, n_src)
    flux = 10.0 ** (-0.4 * (tm - S.ZP)) / 0.00546689
    for shift in ((0.0, 0.0), (0.20, -0.12)):
        img = S.render_model(psf, tx + shift[0], ty + shift[1], flux, (100, 100), pad) + 150.0
        dn = (img + rng.normal(0, 1, img.shape) * np.sqrt(img[pad:-pad, pad:-pad].mean() / 4.6 + 1.0))[pad:-pad, pad:-pad]
        sx, sy, chi_best, chi_zero = S.rigid_shift(dn, 4.6, 1.0, psf, tx, ty, tm)
        assert abs(sx - shift[0]) < 0.05 and abs(sy - shift[1]) < 0.05 and chi_best <= chi_zero + 1e-6


# ----------------------------------------------------------------------------- data tests (skipped when the data are absent)
def _have(*paths):
    return all(os.path.exists(p) for p in paths)


def test_m2_mappings_and_units_reproduce_the_shipped_m2_products():
    pytest.importorskip('astropy')
    P = S._m2_paths()
    ref_truth, ref_patch = os.path.join(S.ROOT, 'truth_hugs_meth1.npz'), os.path.join(M2_PATCH, 'patch_arrays.npz')
    if not _have(P['frame'], P['psfield'], P['catalog'], ref_truth, ref_patch):
        pytest.skip('the M2 frame, psField, HUGS catalogue, truth file or patch are not on disk ($ASTRO_DATA_2 / $ASTRO_DATA)')
    fr, ps, cat = S.read_frame(P['frame']), S.read_psfield(P['psfield']), S.read_hugs(P['catalog'])
    th, ta = S.make_truth(cat, S.world2pix_for(fr, 'header')), S.make_truth(cat, S.world2pix_for(fr, 'astrans'))
    ref = np.load(ref_truth)
    assert len(th['px']) == len(ref['px']) == 287653
    assert np.abs(th['px'] - ref['px']).max() < 1e-6 and np.abs(th['py'] - ref['py']).max() < 1e-6          # the header mapping is the one of data/astro/make_truth.py
    assert np.array_equal(th['m606'], ref['m606']) and np.array_equal(th['e606'], ref['e606']) and np.array_equal(th['nf606'], ref['nf606'])
    # the full solution against the header mapping in the benchmark patch (reported: median 0.0097 px, 90th percentile 0.0199 px)
    box = (th['px'] >= 309.5) & (th['px'] < 409.5) & (th['py'] >= 629.5) & (th['py'] < 729.5)
    d = np.hypot(ta['px'] - th['px'], ta['py'] - th['py'])[box]
    assert abs(np.median(d) - 0.0099) < 0.003 and abs(np.percentile(d, 90) - 0.0199) < 0.004 and d.max() < 0.04
    # the units: the benchmark patch rebuilt from the frame and the psField is the patch of fetch.py, bit for bit
    arrays, meta = S.extract_patch(fr, ps, 310, 630)
    old = np.load(ref_patch)
    assert set(old.files) == set(arrays) and all(np.array_equal(np.asarray(arrays[k]), old[k]) for k in old.files)
    assert meta == json.load(open(os.path.join(M2_PATCH, 'meta.json')))
    from astropy.io import fits
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:                                       # the files written by write_patch are those of fetch.py (plus the NMGY card)
        S.write_patch(tmp, arrays, meta, fr['header'])
        h_new, h_old = fits.getheader(os.path.join(tmp, 'patch_dn.fits')), fits.getheader(os.path.join(M2_PATCH, 'patch_dn.fits'))
        assert all(h_old[k] == h_new[k] for k in h_old) and [k for k in h_new if k not in h_old] == ['NMGY']
        assert np.array_equal(fits.getdata(os.path.join(tmp, 'patch_dn.fits')), fits.getdata(os.path.join(M2_PATCH, 'patch_dn.fits')))
        assert open(os.path.join(tmp, 'meta.json')).read() == open(os.path.join(M2_PATCH, 'meta.json')).read()
    # the independent shipped truth (reported for F606W < 20: median |HUGS - shipped| 0.0519 px (header) and 0.0435 px (asTrans), 90th percentile 0.060 / 0.057; magnitude offset +0.005)
    sp = os.path.join(SHIPPED, 'Data', 'sdss0921_true.txt')
    if _have(sp):
        st = np.loadtxt(sp)
        sx, sy, sm = st[:, 0] + 310, st[:, 1] + 630, S.counts_to_mag(st[:, 2], 0.00546689)
        for t, med in ((th, 0.0519), (ta, 0.0435)):
            c = S.compare_with_list(t, sx, sy, sm, 20)
            assert c['n'] == 234 and c['n_matched'] >= 232
            assert np.median(c['d']) < 0.15 and abs(np.median(c['d']) - med) < 0.006 and np.percentile(c['d'], 90) < 0.08         # the acceptance criterion: median below 0.15 px
            assert abs(np.median(c['dm'])) < 0.01 and S.robust_scatter(c['dm']) < 0.01


# the patches of the second cluster as reported in results/astro/second_cluster_spec.md: (x0, y0) -> (F606W < 22 in the box, scored, all valid F606W in the box)
REPORTED = {
    'ngc5024': dict(frame=(5390, 6, 162), patches={(760, 590): (402, 365, 4527), (640, 680): (671, 594, 6427), (400, 660): (950, 858, 8024), (660, 420): (1302, 1163, 10004)}),
    'ngc6341': dict(frame=(4682, 6, 208), patches={(1100, 220): (718, 638, 3141), (1430, 550): (1184, 1052, 4460)}),
}


@pytest.mark.parametrize('cluster', sorted(REPORTED))
def test_second_cluster_patches_on_disk_match_the_spec(cluster):
    pytest.importorskip('astropy')
    run, cam, fld = REPORTED[cluster]['frame']
    truth_path = os.path.join(S.ROOT2, 'truth_hugs_%s_meth1.npz' % cluster)
    dirs = {xy: os.path.join(S.ROOT2, S.patch_name(run, cam, fld, *xy)) for xy in REPORTED[cluster]['patches']}
    if not _have(truth_path, *[os.path.join(d, 'patch_arrays.npz') for d in dirs.values()]):
        pytest.skip('the patches of %s are not in $ASTRO_DATA_2 (python -I experiments/astro_second.py build --cluster %s)' % (cluster, cluster))
    tr = np.load(truth_path)
    assert set(tr.files) == {'px', 'py', 'm606', 'e606', 'nf606'}
    old_keys = {'dn', 'electrons', 'nmgy', 'sky_dn', 'psf51', 'psf25', 'gain', 'dark_variance', 'x0', 'y0'}
    for (x0, y0), (n22, scored, n_all) in REPORTED[cluster]['patches'].items():
        z = np.load(os.path.join(dirs[(x0, y0)], 'patch_arrays.npz'))
        assert set(z.files) == old_keys and z['dn'].shape == (100, 100) and z['psf51'].shape == (51, 51) and int(z['x0']) == x0 and int(z['y0']) == y0
        assert abs(z['psf51'].sum() - 1) < 1e-12 and np.allclose(z['electrons'], z['dn'] * float(z['gain']))
        t = S.truth_for_patch(tr, x0, y0)
        assert int((t['mag'] < 22).sum()) == n22 and int(((t['mag'] < 22) & t['inner']).sum()) == scored and len(t['mag']) == n_all
        assert 0.03 <= n22 / 1e4 <= 0.15 and t['mag'].min() >= S.BRIGHT_LIMIT and float(z['dn'].max()) < S.SAT_FRACTION * 52900.0
        assert sorted(os.listdir(dirs[(x0, y0)])) == ['meta.json', 'patch_arrays.npz', 'patch_dn.fits']
    # the units of one patch rebuilt from the downloaded frame (needs the frame and the psField)
    P = S._target_paths(cluster)
    if _have(P['frame'], P['psfield']):
        fr, ps = S.read_frame(P['frame']), S.read_psfield(P['psfield'])
        (x0, y0) = sorted(REPORTED[cluster]['patches'])[0]
        arrays, _ = S.extract_patch(fr, ps, x0, y0)
        z = np.load(os.path.join(dirs[(x0, y0)], 'patch_arrays.npz'))
        assert all(np.array_equal(np.asarray(arrays[k]), z[k]) for k in z.files)


# ----------------------------------------------------------------------------- data/astro/fetch_second.py (no network)
def _fetch_second():
    import importlib.util
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'astro', 'fetch_second.py')
    spec = importlib.util.spec_from_file_location('fetch_second', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_fetch_second_items_urls_directories_and_the_m2_hashes_of_fetch_py():
    import urllib.parse
    F = _fetch_second()
    items = F.build_items()
    assert items['sdss_frame_r_005390_6_0162'][2] == 'https://data.sdss.org/sas/dr17/eboss/photoObj/frames/301/5390/6/frame-r-005390-6-0162.fits.bz2'
    assert items['sdss_psfield_004682_6_0208'][2] == 'https://data.sdss.org/sas/dr17/eboss/photo/redux/301/4682/objcs/6/psField-004682-6-0208.fit'
    assert items['sdss_photofield_004682_6'][2] == 'https://data.sdss.org/sas/dr17/eboss/photoObj/301/4682/photoField-004682-6.fits'
    assert urllib.parse.parse_qs(urllib.parse.urlparse(items['hugs_ngc5024'][2]).query)['uri'] == ['mast:HLSP/hugs/ngc5024/hlsp_hugs_hst_wfc3-uvis-acs-wfc_ngc5024_multi_v1_catalog-meth1.txt']
    assert len({os.path.dirname(F.dest_path(items, k)) for k in items}) == len(items)        # every download goes to its own directory
    assert {v[0] for v in items.values()} == set(F.GROUPS)
    assert all(('ngc%d' % n) in [v[1] for v in items.values()] for n in F.CLUSTERS)
    for (_, _, run, cam, fld) in [(0, 0, *v) for v in F.TARGETS.values()]:
        assert 'sdss_frame_r_%06d_%d_%04d' % (run, cam, fld) in items and 'sdss_psfield_%06d_%d_%04d' % (run, cam, fld) in items and 'sdss_photofield_%06d_%d' % (run, cam) in items
    old_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'astro', 'checksums.json')
    if os.path.exists(old_path):                                                              # the M2 files are the bytes that fetch.py recorded
        old = json.load(open(old_path))['files']
        for k2, k1 in (('sdss_frame_r_002583_2_0136', 'sdss_frame_r'), ('sdss_psfield_002583_2_0136', 'sdss_psfield'), ('sdss_photofield_002583_2', 'sdss_photofield'), ('hugs_ngc7089', 'hugs_ngc7089_meth1')):
            assert F.M2_EXPECTED[k2] == (old[k1]['sha256'], old[k1]['bytes'])
    if os.path.exists(F.CHECKSUMS):
        assert set(json.load(open(F.CHECKSUMS))['files']) <= set(items)


def test_fetch_second_download_resumes_with_range_records_and_refuses_a_corrupted_file(tmp_path, monkeypatch):
    import hashlib
    import http.server
    import threading
    F = _fetch_second()
    payload = np.random.default_rng(0).bytes(300000)
    seen = []

    class H(http.server.BaseHTTPRequestHandler):
        def do_GET(self):
            rng = self.headers.get('Range')
            seen.append(rng)
            start = int(rng.split('=')[1].split('-')[0]) if rng else 0
            body = payload[start:]
            self.send_response(206 if rng else 200)
            if rng:
                self.send_header('Content-Range', 'bytes %d-%d/%d' % (start, len(payload) - 1, len(payload)))
            self.send_header('Content-Length', str(len(body)))
            self.send_header('Last-Modified', 'Wed, 01 Jan 2020 00:00:00 GMT')
            self.end_headers()
            self.wfile.write(body)

        def log_message(self, *a):
            pass
    srv = http.server.HTTPServer(('127.0.0.1', 0), H)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    for k in ('http_proxy', 'HTTP_PROXY'):
        monkeypatch.delenv(k, raising=False)
    monkeypatch.setenv('no_proxy', '127.0.0.1')
    monkeypatch.setenv('ASTRO_DATA_2', str(tmp_path))
    monkeypatch.setattr(F, 'CHECKSUMS', str(tmp_path / 'sums.json'))                          # never touch the real checksums_second.json
    try:
        items = {'k': ('target', 'd', 'http://127.0.0.1:%d/file.bin' % srv.server_address[1], 'file.bin')}
        dest = F.dest_path(items, 'k')
        os.makedirs(os.path.dirname(dest))
        with open(dest + '.part', 'wb') as f:                                                 # an interrupted download: the first 100,000 bytes
            f.write(payload[:100000])
        lastmod = F.download(items['k'][2], dest)
        assert seen == ['bytes=100000-'] and open(dest, 'rb').read() == payload and lastmod == 'Wed, 01 Jan 2020 00:00:00 GMT' and not os.path.exists(dest + '.part')
        sums = F.load_sums()
        F.check_or_record(items, 'k', dest, sums, lastmod)                                    # first download: recorded
        rec = json.load(open(F.CHECKSUMS))['files']['k']
        assert rec['bytes'] == len(payload) and rec['sha256'] == hashlib.sha256(payload).hexdigest() and rec['source_last_modified'] == lastmod
        F.check_or_record(items, 'k', dest, F.load_sums())                                    # verified against the record
        with open(dest, 'r+b') as f:                                                          # a flipped byte: refused and moved aside
            f.seek(1234)
            f.write(b'\x00\x01')
        with pytest.raises(SystemExit):
            F.check_or_record(items, 'k', dest, F.load_sums())
        assert os.path.exists(dest + '.bad') and not os.path.exists(dest)
    finally:
        srv.shutdown()
