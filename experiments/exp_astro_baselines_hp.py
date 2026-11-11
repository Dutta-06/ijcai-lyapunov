"""
Catalogues of the two standard crowded-field photometry tools given the ISS's front end (the border-consistent high-pass) and a wider pipeline search (needs the DOMAIN environment: astropy, photutils, crowdsource)
=====================================================================================================================================================================================================

The referee's remark behind it (round R3, finding F7): the ISS's largest effect, the 3 px border-consistent high-pass, was available to no baseline, and no pipeline parameter of the baselines was searched.
The family, its parameters, the selection rule and what is reported were declared BEFORE any catalogue of this study was scored, in `results/astro/baselines_hp_plan.md` (read it first; its SHA-256 is stored with
the selection). This script writes the catalogues; `exp_astro_hp.py` scores them, selects on the two M2 development patches and writes the report. It imports `exp_astro_baselines` (hash-frozen) and
`imaging` read-only; no truth enters any run.

The front end is exactly the ISS's (`exp_astro_m2.run_path`): data `H (dn - median)` with `H y = y - S y S^T`, `S = imaging.smoother(n, sigma_b)`, computed by `imaging.PixelPSF.highpass`; the weights
`sw = 1 / sqrt(max(gaussian_filter(dn, 1.5), 1st percentile) / gain + dark)`; the model of a star of unit flux `H` applied to the 51 x 51 psField stamp (embedded in a zero frame wide enough for the smoother's
4 sigma_b truncation, so the result is the interior operator, then cropped).

  photutils     8 variants (`PHOT_HP`): DAOStarFinder(threshold = 5 x median(1/sw), fwhm of the stamp) + IterativePSFPhotometry(ImagePSF(H stamp, origin (25, 25)), fit_shape 7, aperture radius 3, no local background),
                sigma_b in {3, 5}; mode 'new' with maxiters in {3, 10, 30}, and mode 'all' (SourceGrouper 4.4 px, relaxed sharpness / roundness cuts) with maxiters 6
  crowdsource   8 variants (`CS_HP`): `fit_im` as `exp_astro_baselines.run_crowdsource` calls it, the detection `threshold` in {3, 5, 8}, on the raw image with the ISS's weights, on the image high-passed with
                sigma_b = 3, and (threshold 5) with sigma_b = 5; the high-passed rows use the high-passed stamp in `SimplePSF(., normalize=0)` (the default's flux convention imposed beforehand); the current
                default (`cs_default`, `exp_astro_baselines.run_crowdsource` called unchanged) is the reference. No post-filter.

Writes (keys of `baseline_catalogs.npz`: '<patch>|<tag>|x', '|y', '|flux', patch = 'x<x0>_y<y0>' in the patch's 0-based pixel coordinates, flux in DN; crowdsource also '|dflux', '|qf', '|fracflux', '|spread',
'|rchi2'; every run also '|secs' = [CPU seconds]). The 18 patch origins of M2, M53, M92 and the repeat frames are distinct (checked).

    python -I experiments/exp_astro_baselines_hp.py --stage dev  --out FILE [--variants TAG ...] [--patches X,Y ...]        every variant on the development patches
    python -I experiments/exp_astro_baselines_hp.py --stage test --out FILE --raw baselines_hp_raw.npz [--sets m2 m53 m92 m2r1 ...]   only the SELECTED variants and the references, on the test patches
    python -I experiments/exp_astro_baselines_hp.py --merge FILE FILE ... --out FILE                                         (union of catalogue files; a key present twice must be equal)
"""
import argparse
import contextlib
import io
import os
import sys
import time
import warnings

os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np

HERE = os.path.abspath(os.path.dirname(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)
import exp_astro_baselines as AB      # noqa: E402  (numpy only at import; astropy, photutils and crowdsource are imported inside the runners)
import exp_astro_m2 as X              # noqa: E402
import astro_repeat as R              # noqa: E402
from src.core import imaging as IM    # noqa: E402

OUT = os.path.join(REPO, 'results', 'astro', 'baseline_catalogs_hp.npz')
PLAN = os.path.join(REPO, 'results', 'astro', 'baselines_hp_plan.md')

NSIG = 5.0                      # DAOStarFinder threshold in noise sigmas (the variance model's median pixel sigma); 3 gives the identical catalogue (plan, section 3 check 4)
FIT, APER, GROUP = 7, 3.0, 4.4  # fit shape, aperture radius of the initial flux, SourceGrouper min_separation of the 'all' variants (px)
RELAXED = dict(sharpness_range=(0.0, 2.0), roundness_range=(-2.0, 2.0))      # the relaxed cuts of `exp_astro_baselines.PHOT_EXT` (sharplo 0, sharphi 2, roundlo -2, roundhi 2)

# photutils on the high-passed patch: tag, sigma_b (px), mode, maxiters, grouper min_separation (0 = none), relaxed finder cuts
PHOT_HP = [dict(tag=t, sb=sb, mode=mode, iters=it, group=g, relax=rx) for t, sb, mode, it, g, rx in (
    ('phu_hp5_new_i3', 5.0, 'new', 3, 0.0, False),            # the reference: the settings of the default `photutils` baseline on the ISS's a priori high-pass
    ('phu_hp3_new_i3', 3.0, 'new', 3, 0.0, False),
    ('phu_hp3_new_i10', 3.0, 'new', 10, 0.0, False),
    ('phu_hp5_new_i10', 5.0, 'new', 10, 0.0, False),
    ('phu_hp3_new_i30', 3.0, 'new', 30, 0.0, False),
    ('phu_hp5_new_i30', 5.0, 'new', 30, 0.0, False),
    ('phu_hp3_all_i6', 3.0, 'all', 6, GROUP, True),         # the options of the best variant of the old family (mode 'all', grouper 4.4 px, relaxed cuts), no local background
    ('phu_hp5_all_i6', 5.0, 'all', 6, GROUP, True))]
PHOT_REF = 'phu_hp5_new_i3'

# crowdsource: tag, kind ('default': exp_astro_baselines.run_crowdsource unchanged; 'raw': raw image with the ISS's weights; 'hp': high-passed image and stamp with the ISS's weights), sigma_b, threshold
CS_HP = [dict(tag=t, kind=k, sb=sb, thr=thr) for t, k, sb, thr in (
    ('cs_default', 'default', 0.0, 5.0),                      # reference 1: the current default
    ('cs_w_t5', 'raw', 0.0, 5.0),
    ('cs_w_t3', 'raw', 0.0, 3.0),
    ('cs_w_t8', 'raw', 0.0, 8.0),
    ('cs_hp3_t5', 'hp', 3.0, 5.0),
    ('cs_hp3_t3', 'hp', 3.0, 3.0),
    ('cs_hp3_t8', 'hp', 3.0, 8.0),
    ('cs_hp5_t5', 'hp', 5.0, 5.0))]                           # reference 2: the default pipeline on the a priori high-pass
CS_REFS = ('cs_default', 'cs_hp5_t5')
FAMILIES = {'photutils': [d['tag'] for d in PHOT_HP], 'crowdsource': [d['tag'] for d in CS_HP]}

# the patch sets: name -> (data set for the loader, patches)
SETS = {'dev': ('m2', list(X.DEV)), 'm2': ('m2', list(X.TESTS)), 'm53': ('m53', list(X.SECOND['m53']['patches'])), 'm92': ('m92', list(X.SECOND['m92']['patches']))}
SETS.update({name: (name, [tuple(p) for p in R.REPEAT[name]['patches']]) for name in sorted(R.REPEAT)})
REPEAT_SETS = tuple(sorted(R.REPEAT))
TEST_SETS = ('m2', 'm53', 'm92') + REPEAT_SETS


def patch_key(x0, y0):
    return 'x%d_y%d' % (x0, y0)


def check_origins():
    """Raises if two patches of the sets share an origin (the key 'x<x0>_y<y0>' would collide); returns the number of origins."""
    origins = [tuple(p) for ds, plist in SETS.values() for p in plist]
    if len(set(origins)) != len(origins):
        raise SystemExit('two patches of the sets share an origin: the keys x<x0>_y<y0> would collide')
    return len(origins)


def load_arrays(ds, x0, y0):
    """dict(dn, psf51, gain, dark, sky) of a patch: `ds` is 'm2', 'm53', 'm92' or a repeat frame name (`astro_repeat.REPEAT`)."""
    path = os.path.join(AB.patch_dir(x0, y0, ds=ds), 'patch_arrays.npz') if ds in ('m2', 'm53', 'm92') else R.patch_path(x0, y0, ds)
    z = np.load(path)
    return dict(dn=z['dn'].astype(float), psf51=z['psf51'].astype(float), gain=float(z['gain']), dark=float(z['dark_variance']), sky=float(z['sky_dn'].mean()))


# ----------------------------------------------------------------------------------------------------------------------------------
# the ISS's front end
# ----------------------------------------------------------------------------------------------------------------------------------
def iss_weights(dn, gain, dark):
    """The ISS's inverse noise (`exp_astro_m2.run_path`): 1 / sqrt(max(gaussian_filter(dn, 1.5, 'reflect'), 1st percentile of dn) / gain + dark variance)."""
    from scipy.ndimage import gaussian_filter
    expected = np.maximum(gaussian_filter(dn, 1.5, mode='reflect'), np.percentile(dn, 1))
    return 1.0 / np.sqrt(expected / gain + dark)


def highpass_frame(dn, sb):
    """The ISS's data: H (dn - median(dn)), by `imaging.PixelPSF.highpass` itself (the border-consistent smoother of width sb)."""
    ny, nx = dn.shape
    return IM.PixelPSF(ny, nx, 1.0, hp=float(sb)).highpass(dn - np.median(dn))


def highpass_stamp(stamp, sb):
    """H applied to a square PSF stamp: the stamp is embedded in a zero frame of side n + 2 ceil(4 sb) (the smoother's truncation radius is 4 sb, so its border renormalisation never reaches a stamp pixel), the
    ISS's operator `x - S x S^T` is applied with `imaging.smoother`, and the middle n x n is returned: the interior operator, the model H phi of a unit-flux star."""
    n = stamp.shape[0]
    pad = int(np.ceil(4.0 * sb))
    m = n + 2 * pad
    big = np.zeros((m, m))
    big[pad:pad + n, pad:pad + n] = stamp
    S = IM.smoother(m, float(sb))
    return (big - S @ big @ S.T)[pad:pad + n, pad:pad + n]


def noise_sigma(dn, gain, dark):
    """The pixel noise of the ISS's variance model, the median of 1 / sw (DN): the unit of the DAOStarFinder threshold."""
    return float(np.median(1.0 / iss_weights(dn, gain, dark)))


# ----------------------------------------------------------------------------------------------------------------------------------
# the runners
# ----------------------------------------------------------------------------------------------------------------------------------
def run_photutils_hp(dn, psf51, gain, dark, spec):
    """One photutils variant of `PHOT_HP` on the high-passed patch; returns the catalogue dict(x, y, flux) of the sources with a finite positive fitted flux and the CPU seconds of the photometry."""
    from photutils.detection import DAOStarFinder
    from photutils.psf import IterativePSFPhotometry, ImagePSF, SourceGrouper
    from astropy.utils.exceptions import AstropyUserWarning
    sb = spec['sb']
    err = 1.0 / iss_weights(dn, gain, dark)
    data = highpass_frame(dn, sb)
    stamp = highpass_stamp(psf51, sb)
    finder = DAOStarFinder(threshold=NSIG * float(np.median(err)), fwhm=AB.psf_fwhm(psf51), **(RELAXED if spec['relax'] else {}))
    model = ImagePSF(stamp, flux=1, x_0=0, y_0=0, origin=(25, 25), oversampling=1)
    grouper = SourceGrouper(min_separation=spec['group']) if spec['group'] else None
    phot = IterativePSFPhotometry(model, fit_shape=(FIT, FIT), finder=finder, grouper=grouper, maxiters=spec['iters'], mode=spec['mode'], aperture_radius=APER, n_threads=1)
    t0 = time.process_time()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', AstropyUserWarning)          # the fits that did not converge are flagged in the table, not printed once per source
        res = phot(data, error=err)
    secs = time.process_time() - t0
    if res is None:
        return dict(x=np.zeros(0), y=np.zeros(0), flux=np.zeros(0)), secs
    good = np.isfinite(res['flux_fit']) & (res['flux_fit'] > 0)
    return dict(x=np.asarray(res['x_fit'][good], float), y=np.asarray(res['y_fit'][good], float), flux=np.asarray(res['flux_fit'][good], float)), secs


def cs_stamp(psf51):
    """The stamp of `exp_astro_baselines.run_crowdsource` before `SimplePSF`: clipped at 1e-10 and normalised to unit sum."""
    stamp = np.clip(psf51, 1e-10, np.inf)
    return stamp / stamp.sum()


def cs_fit(im, psf, weights, thr):
    """crowdsource's `fit_im` as `exp_astro_baselines.run_crowdsource` calls it, with the detection threshold `thr`; the printed centroid lines are discarded; the output with x and y swapped (crowdsource's 'x' runs
    along the array's first axis), as there."""
    from crowdsource import crowdsource_base
    flag = np.zeros(im.shape, dtype=np.int64)
    with contextlib.redirect_stdout(io.StringIO()):
        res = crowdsource_base.fit_im(im, psf, weights=np.asarray(weights, dtype=np.float64), dq=flag, maxstars=10000000, fewstars=60, ntilex=1, ntiley=1, refit_psf=False, verbose=False, threshold=thr)
    cat = res['stars']
    return dict(x=np.asarray(cat['y'], float), y=np.asarray(cat['x'], float), flux=np.asarray(cat['flux'], float), dflux=np.asarray(cat['dflux'], float), qf=np.asarray(cat['qf'], float),
                fracflux=np.asarray(cat['fracflux'], float), spread=np.asarray(cat['spread_model'], float), rchi2=np.asarray(cat['rchi2'], float))


def run_crowdsource_hp(dn, sky_pipe, psf51, gain, dark, spec):
    """One crowdsource variant of `CS_HP`; returns the catalogue dict and the CPU seconds."""
    t0 = time.process_time()
    if spec['kind'] == 'default':
        with contextlib.redirect_stdout(io.StringIO()):
            cat, _ = AB.run_crowdsource(dn, sky_pipe, psf51, gain, dark)
        return cat, time.process_time() - t0
    import crowdsource.psf as psfmod
    sw = iss_weights(dn, gain, dark)
    base = cs_stamp(psf51)
    if spec['kind'] == 'raw':
        cat = cs_fit(dn - sky_pipe, psfmod.SimplePSF(base.copy()), sw, spec['thr'])
    elif spec['kind'] == 'hp':
        s19 = float(psfmod.central_stamp(base, 19).sum())                   # the default's flux convention: SimplePSF(normalize=19) makes the central 19 x 19 pixels sum to 1
        psf = psfmod.SimplePSF(highpass_stamp(base / s19, spec['sb']), normalize=0)
        cat = cs_fit(highpass_frame(dn, spec['sb']), psf, sw, spec['thr'])
    else:
        raise ValueError(spec['kind'])
    return cat, time.process_time() - t0


def run_variant(tag, a):
    """The catalogue and CPU seconds of the variant `tag` on the arrays `a` (`load_arrays`)."""
    for spec in PHOT_HP:
        if spec['tag'] == tag:
            return run_photutils_hp(a['dn'], a['psf51'], a['gain'], a['dark'], spec)
    for spec in CS_HP:
        if spec['tag'] == tag:
            return run_crowdsource_hp(a['dn'], a['sky'], a['psf51'], a['gain'], a['dark'], spec)
    raise KeyError(tag)


def store(out, x0, y0, tag, cat, secs):
    k = '%s|%s|' % (patch_key(x0, y0), tag)
    for name, v in cat.items():
        out[k + name] = np.asarray(v)
    out[k + 'secs'] = np.array([secs])


def tags_of_test_stage(raw):
    """The tags run on every test patch: the selected variant of each family and the references, in a fixed order, without repeats."""
    import exp_astro_hp as H
    sel = H.verify_selection(raw)
    tags = [sel['photutils'], PHOT_REF, sel['crowdsource']] + list(CS_REFS)
    return [t for i, t in enumerate(tags) if t not in tags[:i]]


# ----------------------------------------------------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--stage', choices=['dev', 'test'])
    ap.add_argument('--out', required=True)
    ap.add_argument('--variants', nargs='*', help='dev stage: a subset of the tags (default: all 16)')
    ap.add_argument('--patches', nargs='*', help='dev stage: x,y origins (default: the two development patches); test stage: restrict a set to these origins')
    ap.add_argument('--sets', nargs='*', help='test stage: the sets (default m2 m53 m92); the repeat frames only with --flag')
    ap.add_argument('--raw', help='test stage: baselines_hp_raw.npz with the frozen selection')
    ap.add_argument('--flag', help='test stage, repeat frames: the file that says the ISS repeat run is complete')
    ap.add_argument('--merge', nargs='*')
    A_ = ap.parse_args()
    if A_.merge:
        AB.merge(A_.merge, A_.out)
        print('wrote', A_.out)
        return
    if not A_.stage:
        ap.error('give --stage or --merge')
    if os.path.exists(A_.out):
        raise SystemExit('%s exists: a catalogue file is never overwritten (merge partial files with --merge)' % A_.out)
    print('%d patch origins, all distinct' % check_origins(), flush=True)
    jobs = []                                                               # (set, data set, (x0, y0), tags)
    if A_.stage == 'dev':
        tags = A_.variants or FAMILIES['photutils'] + FAMILIES['crowdsource']
        unknown = [t for t in tags if t not in FAMILIES['photutils'] + FAMILIES['crowdsource']]
        if unknown:
            raise SystemExit('unknown variants: %s' % unknown)
        plist = [tuple(int(v) for v in p.split(',')) for p in A_.patches] if A_.patches else SETS['dev'][1]
        jobs = [('dev', 'm2', p, tags) for p in plist]
    else:
        if not A_.raw:
            raise SystemExit('--stage test needs --raw (the file with the frozen selection)')
        import exp_astro_m2 as _X
        raw = _X.load_raw(A_.raw)
        tags = tags_of_test_stage(raw)
        names = A_.sets or ['m2', 'm53', 'm92']
        for name in names:
            if name not in TEST_SETS:
                raise SystemExit('unknown set %s' % name)
            if name in REPEAT_SETS and not (A_.flag and os.path.exists(A_.flag)):
                raise SystemExit('the repeat frames are run only after the flag file of the ISS repeat run exists (--flag)')
        for name in names:
            ds, plist = SETS[name]
            if A_.patches:
                want = [tuple(int(v) for v in p.split(',')) for p in A_.patches]
                plist = [p for p in plist if tuple(p) in want]
            jobs += [(name, ds, tuple(p), tags) for p in plist]
        print('test stage: variants %s' % tags, flush=True)
    out = {}
    for name, ds, (x0, y0), tags in jobs:
        a = load_arrays(ds, x0, y0)
        for tag in tags:
            cat, secs = run_variant(tag, a)
            store(out, x0, y0, tag, cat, secs)
            print('%s patch (%d, %d) %s: %d sources, %.1f s CPU' % (name, x0, y0, tag, len(cat['x']), secs), flush=True)
    os.makedirs(os.path.dirname(os.path.abspath(A_.out)), exist_ok=True)
    np.savez_compressed(A_.out, **out)
    print('wrote', A_.out, '(%d arrays)' % len(out))


if __name__ == '__main__':
    main()
