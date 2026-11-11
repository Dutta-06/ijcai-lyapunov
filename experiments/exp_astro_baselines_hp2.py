"""
Round 2 of the baselines given the ISS's front end: photutils with `xy_bounds`, crowdsource with more passes and a 2 px high-pass (needs the DOMAIN environment: astropy, photutils, crowdsource)
===========================================================================================================================================================================================

Declared before any round-2 catalogue was made or scored in `results/astro/baselines_hp2_plan.md` (read it first; its SHA-256 is stored with the selection). Round 1 (`exp_astro_baselines_hp.py`, imported
read-only: the front end, the loaders, the patch sets, the registries) left out two kinds of options that matter on the development patches (referee round BHP, findings 1 and 2):

  photutils     the round-1 variants phu_hp3_all_i6, phu_hp3_new_i30, phu_hp5_new_i3, each with `IterativePSFPhotometry(xy_bounds=b)`, b in {1, 2} px (6 variants, `PHOT2`)
  crowdsource   the round-1 call with the high-pass width in {2, 3} px x threshold in {3, 5} x passes in {d: miniter 4, maxiter 10, fewstars 60; m: miniter 8, maxiter 20, fewstars 20} (8 variants, `CS2`)

Writes (keys of the round-1 catalogue file: '<patch>|<tag>|x', '|y', '|flux', crowdsource also '|dflux', '|qf', '|fracflux', '|spread', '|rchi2', every run '|secs'). In the test stage a selected variant that is a
round-1 variant is not run: its arrays are copied from the round-1 catalogue file (identical by construction), so that the round-2 file is complete.

    python -I experiments/exp_astro_baselines_hp2.py --stage dev  --out FILE [--variants TAG ...] [--patches X,Y ...]
    python -I experiments/exp_astro_baselines_hp2.py --stage test --out FILE --raw baselines_hp2_raw.npz [--sets m2 m53 m92 m2r1 ...] [--flag FILE] [--patches X,Y ...]
    python -I experiments/exp_astro_baselines_hp2.py --merge FILE FILE ... --out FILE
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
import exp_astro_baselines as AB          # noqa: E402
import exp_astro_baselines_hp as G        # noqa: E402  (round 1, read only)

PLAN2 = os.path.join(REPO, 'results', 'astro', 'baselines_hp2_plan.md')
CATS1 = os.path.join(REPO, 'results', 'astro', 'baseline_catalogs_hp.npz')

# photutils: a round-1 variant with xy_bounds (px)
PHOT2_PARENTS = ('phu_hp3_all_i6', 'phu_hp3_new_i30', 'phu_hp5_new_i3')
BOUNDS = (1, 2)
PHOT2 = [dict(next(d for d in G.PHOT_HP if d['tag'] == parent), tag='phu2_%s_b%d' % (parent[len('phu_'):], b), parent=parent, bounds=float(b)) for parent in PHOT2_PARENTS for b in BOUNDS]

# crowdsource: the round-1 high-passed call with the width, the threshold and the passes varied
PASSES = {'d': dict(miniter=4, maxiter=10, fewstars=60), 'm': dict(miniter=8, maxiter=20, fewstars=20)}
CS2 = [dict(tag='cs2_hp%d_t%d_p%s' % (sb, thr, ps), kind='hp', sb=float(sb), thr=float(thr), passes=ps) for sb in (2, 3) for thr in (3, 5) for ps in ('d', 'm')]
FAMILIES2 = {'photutils': [d['tag'] for d in PHOT2], 'crowdsource': [d['tag'] for d in CS2]}


def run_photutils_hp2(dn, psf51, gain, dark, spec):
    """One photutils variant of `PHOT2`: `exp_astro_baselines_hp.run_photutils_hp` with `xy_bounds=spec['bounds']`; the catalogue dict(x, y, flux) of the sources with a finite positive fitted flux and the CPU seconds."""
    from photutils.detection import DAOStarFinder
    from photutils.psf import IterativePSFPhotometry, ImagePSF, SourceGrouper
    from astropy.utils.exceptions import AstropyUserWarning
    sb = spec['sb']
    err = 1.0 / G.iss_weights(dn, gain, dark)
    data = G.highpass_frame(dn, sb)
    stamp = G.highpass_stamp(psf51, sb)
    finder = DAOStarFinder(threshold=G.NSIG * float(np.median(err)), fwhm=AB.psf_fwhm(psf51), **(G.RELAXED if spec['relax'] else {}))
    model = ImagePSF(stamp, flux=1, x_0=0, y_0=0, origin=(25, 25), oversampling=1)
    grouper = SourceGrouper(min_separation=spec['group']) if spec['group'] else None
    phot = IterativePSFPhotometry(model, fit_shape=(G.FIT, G.FIT), finder=finder, grouper=grouper, maxiters=spec['iters'], mode=spec['mode'], aperture_radius=G.APER, n_threads=1,
                                  xy_bounds=spec['bounds'])
    t0 = time.process_time()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', AstropyUserWarning)
        res = phot(data, error=err)
    secs = time.process_time() - t0
    if res is None:
        return dict(x=np.zeros(0), y=np.zeros(0), flux=np.zeros(0)), secs
    good = np.isfinite(res['flux_fit']) & (res['flux_fit'] > 0)
    return dict(x=np.asarray(res['x_fit'][good], float), y=np.asarray(res['y_fit'][good], float), flux=np.asarray(res['flux_fit'][good], float)), secs


def cs_fit2(im, psf, weights, thr, passes):
    """crowdsource's `fit_im` as the round-1 `cs_fit` calls it, with the passes of `PASSES[passes]` (`miniter`, `maxiter`, `fewstars`) and the detection threshold `thr`; x and y swapped as there."""
    from crowdsource import crowdsource_base
    p = PASSES[passes]
    flag = np.zeros(im.shape, dtype=np.int64)
    with contextlib.redirect_stdout(io.StringIO()):
        res = crowdsource_base.fit_im(im, psf, weights=np.asarray(weights, dtype=np.float64), dq=flag, maxstars=10000000, fewstars=p['fewstars'], miniter=p['miniter'], maxiter=p['maxiter'],
                                      ntilex=1, ntiley=1, refit_psf=False, verbose=False, threshold=thr)
    cat = res['stars']
    return dict(x=np.asarray(cat['y'], float), y=np.asarray(cat['x'], float), flux=np.asarray(cat['flux'], float), dflux=np.asarray(cat['dflux'], float), qf=np.asarray(cat['qf'], float),
                fracflux=np.asarray(cat['fracflux'], float), spread=np.asarray(cat['spread_model'], float), rchi2=np.asarray(cat['rchi2'], float))


def run_crowdsource_hp2(dn, psf51, gain, dark, spec):
    """One crowdsource variant of `CS2`: the high-passed image and stamp (round-1 flux convention) with the ISS weights; returns the catalogue dict and the CPU seconds."""
    import crowdsource.psf as psfmod
    t0 = time.process_time()
    sw = G.iss_weights(dn, gain, dark)
    base = G.cs_stamp(psf51)
    s19 = float(psfmod.central_stamp(base, 19).sum())
    psf = psfmod.SimplePSF(G.highpass_stamp(base / s19, spec['sb']), normalize=0)
    cat = cs_fit2(G.highpass_frame(dn, spec['sb']), psf, sw, spec['thr'], spec['passes'])
    return cat, time.process_time() - t0


def run_variant2(tag, a):
    """The catalogue and CPU seconds of the round-2 variant `tag` on the arrays `a` (`exp_astro_baselines_hp.load_arrays`)."""
    for spec in PHOT2:
        if spec['tag'] == tag:
            return run_photutils_hp2(a['dn'], a['psf51'], a['gain'], a['dark'], spec)
    for spec in CS2:
        if spec['tag'] == tag:
            return run_crowdsource_hp2(a['dn'], a['psf51'], a['gain'], a['dark'], spec)
    raise KeyError(tag)


def test_stage_tags(raw2):
    """The tags of the two selected variants of the round-2 raw (after `exp_astro_hp2.verify_selection2`), without repeats."""
    import exp_astro_hp2 as H2
    sel = H2.verify_selection2(raw2)
    tags = [sel['photutils'], sel['crowdsource']]
    return [t for i, t in enumerate(tags) if t not in tags[:i]]


def carry(out, x0, y0, tag, cats1):
    """Copies the arrays of a round-1 variant on a patch from the round-1 catalogue file into `out`."""
    k = '%s|%s|' % (G.patch_key(x0, y0), tag)
    found = [f for f in cats1.files if f.startswith(k)]
    if not found:
        raise SystemExit('the round-1 catalogue file has no arrays %s*' % k)
    for f in found:
        out[f] = np.asarray(cats1[f])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--stage', choices=['dev', 'test'])
    ap.add_argument('--out', required=True)
    ap.add_argument('--variants', nargs='*')
    ap.add_argument('--patches', nargs='*')
    ap.add_argument('--sets', nargs='*')
    ap.add_argument('--raw')
    ap.add_argument('--flag')
    ap.add_argument('--cats1', default=CATS1)
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
    print('%d patch origins, all distinct' % G.check_origins(), flush=True)
    jobs = []                                                                  # (set, data set, (x0, y0), tags)
    all_tags = FAMILIES2['photutils'] + FAMILIES2['crowdsource']
    if A_.stage == 'dev':
        tags = A_.variants or all_tags
        unknown = [t for t in tags if t not in all_tags]
        if unknown:
            raise SystemExit('unknown variants: %s' % unknown)
        plist = [tuple(int(v) for v in p.split(',')) for p in A_.patches] if A_.patches else G.SETS['dev'][1]
        jobs = [('dev', 'm2', p, tags) for p in plist]
    else:
        if not A_.raw:
            raise SystemExit('--stage test needs --raw (the file with the frozen selection)')
        import exp_astro_m2 as X
        tags = test_stage_tags(X.load_raw(A_.raw))
        names = A_.sets or ['m2', 'm53', 'm92']
        for name in names:
            if name not in G.TEST_SETS:
                raise SystemExit('unknown set %s' % name)
            if name in G.REPEAT_SETS and not (A_.flag and os.path.exists(A_.flag)):
                raise SystemExit('the repeat frames are run only after the flag file of the ISS repeat run exists (--flag)')
        for name in names:
            ds, plist = G.SETS[name]
            if A_.patches:
                want = [tuple(int(v) for v in p.split(',')) for p in A_.patches]
                plist = [p for p in plist if tuple(p) in want]
            jobs += [(name, ds, tuple(p), tags) for p in plist]
        print('test stage: variants %s' % tags, flush=True)
    cats1 = np.load(A_.cats1) if A_.stage == 'test' else None
    out = {}
    for name, ds, (x0, y0), tags in jobs:
        a = None
        for tag in tags:
            if tag not in all_tags:                                            # a selected round-1 variant: carried, not run again
                carry(out, x0, y0, tag, cats1)
                print('%s patch (%d, %d) %s: carried from the round-1 catalogue file' % (name, x0, y0, tag), flush=True)
                continue
            a = a or G.load_arrays(ds, x0, y0)
            cat, secs = run_variant2(tag, a)
            G.store(out, x0, y0, tag, cat, secs)
            print('%s patch (%d, %d) %s: %d sources, %.1f s CPU' % (name, x0, y0, tag, len(cat['x']), secs), flush=True)
    os.makedirs(os.path.dirname(os.path.abspath(A_.out)), exist_ok=True)
    np.savez_compressed(A_.out, **out)
    print('wrote', A_.out, '(%d arrays)' % len(out))


if __name__ == '__main__':
    main()
