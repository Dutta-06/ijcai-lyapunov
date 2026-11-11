"""
Catalogues of two standard crowded-field photometry tools on the M2 patches, for the comparison of `exp_astro_m2.py` (needs the DOMAIN environment: astropy, photutils, crowdsource)
=====================================================================================================================================================================================

Settings as set up and run by the data-access round (`results/astro/data_access.md`, section 4), no tuning, the SDSS psField PSF of the patch:

  photutils     DAOStarFinder(threshold = 5 x psField skysig, fwhm = that of the PSF) + IterativePSFPhotometry(ImagePSF(psf, origin = (25, 25)), fit_shape (7, 7),
                maxiters 3, no grouper, mode 'new', aperture_radius 3) on the patch minus its 1st percentile, errors sqrt(DN / gain + dark variance)
  crowdsource   crowdsource_base.fit_im(image - pipeline sky, SimplePSF(psf stamp), weights = 1 / sigma, maxstars 1e7, fewstars 60, ntilex = ntiley = 1, refit_psf False)

Writes results/astro/baseline_catalogs.npz (keys '<patch>|<method>|x', '|y', '|flux', in the patch's 0-based pixel coordinates, flux in DN) for every patch given.
For the benchmark's patch (310, 630) it also stores the 356-source DAOPHOT catalogue that ships with the benchmark's code (`daophot`; the file of the multiband PCAT repository,
Data/idR-002583-2-0136/truth/idR-002583-2-0136-tru.txt, SHA-256 below; columns x, y, r flux in counts, ...), with the file's y coordinate lowered by 1 pixel: the file's rows
are one pixel above the frame's 0-based rows (its own alignment with the HST truth: 303 of its 356 sources lie within 0.5 px of a truth source at that shift, 5 at none; a published
catalogue is given its best alignment, the convention of the file, not tuned per source).

Referee round IS1a (pre-test review) found that the 12 photutils variants above are 6 distinct catalogues (the threshold axis is inert: the patch minus its 1st percentile is far above
every threshold), that the family is under-tuned, and that the DAOPHOT row carries a shift fitted to the truth. So the family is extended by three documented-option variants
(`PHOT_EXT`: a 25 px Background2D with relaxed finder cuts and photutils' LocalBackground in the fit; a percentile-filter background with relaxed cuts and a 4.4 px SourceGrouper;
photutils' DAOPHOT-style mode 'all' with a grouper, LocalBackground and relaxed cuts), chosen from what the documentation offers and scored on the DEVELOPMENT patches only by
`exp_astro_m2.py --only DB`; and the catalogue as shipped (no shift) is stored beside the shifted one (`daophot_asis`). No truth enters any baseline run.

    python -I experiments/exp_astro_baselines.py [--patches 310,630 110,630 210,630] [--out FILE] [--only {base,ext}]     (patches are produced by `data/astro/fetch.py extract-patch --x0 X --y0 Y`)
    python -I experiments/exp_astro_baselines.py --merge FILE FILE ... --out FILE                                       (merge the catalogues of runs on different patches / families)
"""
import argparse
import json
import os
import time

import numpy as np

ROOT = os.environ.get('ASTRO_DATA', '/home/user/astro_data')
HERE = os.path.abspath(os.path.dirname(__file__))
OUT = os.path.join(HERE, '..', 'results', 'astro', 'baseline_catalogs.npz')
DAOPHOT = os.environ.get('DAOPHOT_FILE', os.path.join(ROOT, 'gh_RichardFeder-zz_multiband_pcat', 'Data', 'idR-002583-2-0136', 'truth', 'idR-002583-2-0136-tru.txt'))
DAOPHOT_SHA256 = '47a1cb145168630c14e197822ae37e4618ca429eff0fcacfc09c7e656aeb4602'
DAOPHOT_YSHIFT = -1.0
PHOT_VARIANTS = [(lb, ns, it) for lb in (False, True) for ns in (3.0, 5.0) for it in (3, 6, 10)]      # (background: global 1st percentile or a local 25 px median, DAOStarFinder threshold in sky sigmas, iterations)


# documented-option extensions (referee round IS1a); fields: tag, background ('b2d' = 25 px Background2D, 'lowmed' = 10th-percentile filter smoothed), relaxed finder cuts, SourceGrouper
# min_separation (px, 0 = none), LocalBackground in the fit, photutils mode, iterations, fit shape, threshold in sky sigmas
PHOT_EXT = [
    dict(tag='photutils_ext_b2d_local_i10', bkg='b2d', relax=True, group=0.0, local=True, mode='new', iters=10, fit=7, nsig=3.0),
    dict(tag='photutils_ext_lowmed_group_i6', bkg='lowmed', relax=True, group=4.4, local=False, mode='new', iters=6, fit=7, nsig=5.0),
    dict(tag='photutils_ext_all_local_i6', bkg='p1', relax=True, group=4.4, local=True, mode='all', iters=6, fit=7, nsig=5.0),
]
EXT_TAGS = tuple(d['tag'] for d in PHOT_EXT)


def phot_tag(lb, ns, it):
    return 'photutils' if (lb, ns, it) == (False, 5.0, 3) else 'photutils_%s_t%g_i%d' % ('local' if lb else 'global', ns, it)


def load_daophot(path=DAOPHOT, yshift=DAOPHOT_YSHIFT):
    """The shipped DAOPHOT catalogue of the benchmark patch as dict(x, y, flux) in the patch's 0-based pixel coordinates (y lowered by `yshift` = 1 px, see the docstring; yshift 0: as
    shipped); None if absent."""
    import hashlib
    if not os.path.exists(path):
        return None
    if hashlib.sha256(open(path, 'rb').read()).hexdigest() != DAOPHOT_SHA256:
        raise ValueError('%s is not the file with the recorded SHA-256' % path)
    d = np.loadtxt(path)
    return dict(x=d[:, 0].astype(float), y=d[:, 1].astype(float) + yshift, flux=d[:, 2].astype(float))


ROOT2 = os.environ.get('ASTRO_DATA_2', '/home/user/astro_data_2')
DATASETS = {'m2': (ROOT, 'patch_r_002583_2_0136'), 'm53': (ROOT2, 'patch_r_005390_6_0162'), 'm92': (ROOT2, 'patch_r_004682_6_0208')}       # (directory of the patches, name prefix): M2 and the second clusters


def patch_dir(x0, y0, n=100, ds='m2'):
    root, prefix = DATASETS[ds]
    return os.path.join(root, '%s_x%d_y%d_%dx%d' % (prefix, x0, y0, n, n))


def psf_fwhm(psf51):
    """The FWHM (pixels) of the 51 x 51 psField stamp from its radial profile in half-pixel rings."""
    yy, xx = np.indices(psf51.shape)
    rr = np.hypot(xx - 25, yy - 25)
    rb = np.arange(0, 12, 0.5)
    prof = np.array([psf51[(rr >= r - 0.25) & (rr < r + 0.25)].mean() if ((rr >= r - 0.25) & (rr < r + 0.25)).any() else np.nan for r in rb])
    half = psf51.max() / 2
    k = int(np.argmax(prof < half))
    return 2 * (rb[k - 1] + (prof[k - 1] - half) / (prof[k - 1] - prof[k]) * 0.5)


def run_photutils_ext(dn, psf51, gain, dark_var, skysig, spec):
    """One of the documented-option variants of `PHOT_EXT` (see the docstring): the same PSF stamp, errors and aperture radius as `run_photutils`."""
    from photutils.detection import DAOStarFinder
    from photutils.psf import IterativePSFPhotometry, ImagePSF, SourceGrouper
    from photutils.background import Background2D, MedianBackground, LocalBackground
    from scipy.ndimage import gaussian_filter, percentile_filter
    fwhm = psf_fwhm(psf51)
    if spec['bkg'] == 'p1':
        data = dn - float(np.percentile(dn, 1))
    elif spec['bkg'] == 'b2d':
        for excl in (10.0, 50.0, 100.0):
            try:
                data = dn - Background2D(dn, (25, 25), filter_size=(3, 3), bkg_estimator=MedianBackground(), exclude_percentile=excl).background
                break
            except ValueError:
                if excl == 100.0:
                    raise
    elif spec['bkg'] == 'lowmed':
        data = dn - gaussian_filter(percentile_filter(dn, 10, size=15, mode='reflect'), 5, mode='reflect')
    else:
        raise ValueError(spec['bkg'])
    err = np.sqrt(np.clip(dn, 1, None) / gain + dark_var)
    cuts = dict(sharplo=0.0, sharphi=2.0, roundlo=-2.0, roundhi=2.0) if spec['relax'] else {}
    finder = DAOStarFinder(threshold=spec['nsig'] * skysig, fwhm=fwhm, **cuts)
    model = ImagePSF(psf51, flux=1, x_0=0, y_0=0, origin=(25, 25), oversampling=1)
    grouper = SourceGrouper(min_separation=spec['group']) if spec['group'] else None
    lb = LocalBackground(6.0, 10.0, MedianBackground()) if spec['local'] else None
    fs = (spec['fit'], spec['fit'])
    phot = IterativePSFPhotometry(model, fit_shape=fs, finder=finder, grouper=grouper, maxiters=spec['iters'], mode=spec['mode'], aperture_radius=3.0, n_threads=1, localbkg_estimator=lb)
    t0 = time.time()
    res = phot(data, error=err)
    good = np.isfinite(res['flux_fit']) & (res['flux_fit'] > 0)
    return dict(x=np.asarray(res['x_fit'][good], float), y=np.asarray(res['y_fit'][good], float), flux=np.asarray(res['flux_fit'][good], float)), time.time() - t0, fwhm


def run_photutils(dn, psf51, gain, dark_var, skysig, nsigma=5.0, maxiters=3, local_bkg=False):
    from photutils.detection import DAOStarFinder
    from photutils.psf import IterativePSFPhotometry, ImagePSF
    fwhm = psf_fwhm(psf51)
    if local_bkg:
        from photutils.background import Background2D, MedianBackground
        for excl in (10.0, 50.0, 100.0):          # Background2D's default excludes boxes with more than 10 % of their pixels clipped; on a sparse patch every box can be excluded (patch (410, 630)): relax it
            try:
                data = dn - Background2D(dn, (25, 25), filter_size=(3, 3), bkg_estimator=MedianBackground(), exclude_percentile=excl).background
                break
            except ValueError:
                if excl == 100.0:
                    raise
    else:
        data = dn - float(np.percentile(dn, 1))
    err = np.sqrt(np.clip(dn, 1, None) / gain + dark_var)
    finder = DAOStarFinder(threshold=nsigma * skysig, fwhm=fwhm)
    model = ImagePSF(psf51, flux=1, x_0=0, y_0=0, origin=(25, 25), oversampling=1)
    phot = IterativePSFPhotometry(model, fit_shape=(7, 7), finder=finder, grouper=None, maxiters=maxiters, mode='new', aperture_radius=3.0, n_threads=1)
    t0 = time.time()
    res = phot(data, error=err)
    good = np.isfinite(res['flux_fit']) & (res['flux_fit'] > 0)
    return dict(x=np.asarray(res['x_fit'][good], float), y=np.asarray(res['y_fit'][good], float), flux=np.asarray(res['flux_fit'][good], float)), time.time() - t0, fwhm


def run_crowdsource(dn, sky_pipe, psf51, gain, dark_var):
    import crowdsource.psf as psfmod
    from crowdsource import crowdsource_base
    stamp = np.clip(psf51, 1e-10, np.inf)
    stamp = stamp / stamp.sum()
    psf = psfmod.SimplePSF(stamp)
    im = dn - sky_pipe
    sig = np.sqrt(np.clip(dn, 1, None) / gain + dark_var)
    flag = np.zeros(im.shape, dtype=np.int64)
    t0 = time.time()
    res = crowdsource_base.fit_im(im, psf, weights=(1.0 / sig).astype(np.float64), dq=flag, maxstars=10000000, fewstars=60, ntilex=1, ntiley=1, refit_psf=False, verbose=False)
    cat = res['stars']
    # crowdsource's 'x' runs along the array's first axis (rows) and 'y' along the second (columns): swapped against the patch's (column, row) convention [checked: 380 of 563
    # sources lie within 0.5 px of a truth source after the swap, 171 without it, which is the level of chance]
    return dict(x=np.asarray(cat['y'], float), y=np.asarray(cat['x'], float), flux=np.asarray(cat['flux'], float), dflux=np.asarray(cat['dflux'], float),
                qf=np.asarray(cat['qf'], float), fracflux=np.asarray(cat['fracflux'], float), spread=np.asarray(cat['spread_model'], float),
                rchi2=np.asarray(cat['rchi2'], float)), time.time() - t0


def merge(files, out):
    """Merge baseline catalogue files (keys are per patch and method, so the union is the merge; a key present twice must be equal)."""
    merged = {}
    for f in files:
        with np.load(f) as z:
            for k in z.files:
                if k in merged and not np.array_equal(merged[k], z[k]):
                    raise ValueError('key %s differs between files' % k)
                merged[k] = z[k]
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    np.savez_compressed(out, **merged)
    return merged


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--patches', nargs='*', default=['310,630', '110,630', '210,630'])
    ap.add_argument('--out', default=OUT)
    ap.add_argument('--dataset', choices=sorted(DATASETS), default='m2', help='m2 (NGC 7089), m53 (NGC 5024) or m92 (NGC 6341: the second clusters, results/astro/second_cluster_spec.md)')
    ap.add_argument('--only', choices=['base', 'ext'], help='base: the 12 original photutils variants, crowdsource and DAOPHOT; ext: the documented-option variants of PHOT_EXT')
    ap.add_argument('--merge', nargs='*')
    A = ap.parse_args()
    if A.merge:
        merge(A.merge, A.out)
        print('wrote', A.out)
        return
    out, meta = {}, {}
    for p in A.patches:
        x0, y0 = (int(v) for v in p.split(','))
        d = patch_dir(x0, y0, ds=A.dataset)
        z = np.load(os.path.join(d, 'patch_arrays.npz'))
        m = json.load(open(os.path.join(d, 'meta.json')))
        dn, psf51 = z['dn'].astype(float), z['psf51'].astype(float)
        gain, dark = float(z['gain']), float(z['dark_variance'])
        if A.only in (None, 'base'):
            for lb, ns, it in PHOT_VARIANTS:
                cat, secs, fwhm = run_photutils(dn, psf51, gain, dark, float(m['psfield_skysig']), nsigma=ns, maxiters=it, local_bkg=lb)
                tag = phot_tag(lb, ns, it)
                for k, v in cat.items():
                    out['x%d_y%d|%s|%s' % (x0, y0, tag, k)] = v
                meta['x%d_y%d|%s' % (x0, y0, tag)] = [len(cat['x']), secs]
            if A.dataset == 'm2' and (x0, y0) == (310, 630):
                for name, ysh in (('daophot', DAOPHOT_YSHIFT), ('daophot_asis', 0.0)):
                    dp = load_daophot(yshift=ysh)
                    if dp is not None:
                        for k, v in dp.items():
                            out['x%d_y%d|%s|%s' % (x0, y0, name, k)] = v
                        meta['x%d_y%d|%s' % (x0, y0, name)] = [len(dp['x']), 0.0]
            cat, secs = run_crowdsource(dn, float(z['sky_dn'].mean()), psf51, gain, dark)
            for k, v in cat.items():
                out['x%d_y%d|crowdsource|%s' % (x0, y0, k)] = v
            meta['x%d_y%d|crowdsource' % (x0, y0)] = [len(cat['x']), secs]
            print('patch (%d, %d): photutils (5 sigma) %d sources %.1f s; crowdsource %d sources %.1f s' % (
                x0, y0, meta['x%d_y%d|photutils' % (x0, y0)][0], meta['x%d_y%d|photutils' % (x0, y0)][1], meta['x%d_y%d|crowdsource' % (x0, y0)][0],
                meta['x%d_y%d|crowdsource' % (x0, y0)][1]), flush=True)
        if A.only in (None, 'ext'):
            for spec in PHOT_EXT:
                cat, secs, fwhm = run_photutils_ext(dn, psf51, gain, dark, float(m['psfield_skysig']), spec)
                for k, v in cat.items():
                    out['x%d_y%d|%s|%s' % (x0, y0, spec['tag'], k)] = v
                meta['x%d_y%d|%s' % (x0, y0, spec['tag'])] = [len(cat['x']), secs]
                print('patch (%d, %d): %s %d sources %.1f s' % (x0, y0, spec['tag'], len(cat['x']), secs), flush=True)
    os.makedirs(os.path.dirname(os.path.abspath(A.out)), exist_ok=True)
    np.savez_compressed(A.out, **out)
    print('wrote', A.out)


if __name__ == '__main__':
    main()
