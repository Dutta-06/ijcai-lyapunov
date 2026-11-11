"""
Catalogues of the two standard crowded-field photometry tools on the repeat patches of M2, at their frozen settings (needs the DOMAIN environment: astropy, photutils, crowdsource)
=============================================================================================================================================================================

The repeat patches (`astro_repeat.py`, `results/astro/repeat_spec.md`) are a held-out test set: the four M2 test windows in other SDSS runs at other seeing. This script writes the catalogues of
the baselines on them with the settings that `results/astro/freeze.json` records, tuned on nothing here:

  photutils     the default of `exp_astro_baselines.run_photutils` (DAOStarFinder at 5 sigma + IterativePSFPhotometry, 3 iterations, tag 'photutils') and the variant that
                `freeze.json` ('tuned_baselines', 'photutils') names, 'photutils_ext_all_local_i6' (`exp_astro_baselines.run_photutils_ext` with its entry of `PHOT_EXT`)
  crowdsource   `exp_astro_baselines.run_crowdsource` (fit_im, unfiltered, tag 'crowdsource', with dflux, qf, fracflux, spread and rchi2); the post-filter of `freeze.json`
                ('tuned_baselines', 'crowdsource', 'crowdsource|snr5|ff0|qf0|spinf': flux >= 5 dflux, no cut on fracflux and qf, no cut on the spread) is applied by the reader, as in
                `exp_astro_m2.cs_dets` (`frozen_catalogues` below does it)

The functions of `exp_astro_baselines.py` are imported, not copied (that file is hashed by the freeze record and is not modified); the PSF stamp, the errors, the sky and the units are those of
the patch files. Writes results/astro/baseline_catalogs_repeat.npz in the key format of baseline_catalogs.npz: '<patch>|<method>|x', '|y', '|flux' (crowdsource: also '|dflux', '|qf', '|fracflux',
'|spread', '|rchi2'), positions in the patch's 0-based pixel coordinates (column, row), flux in DN, <patch> = 'x<x0>_y<y0>' (the origins of the patches of all the frames of the registry are distinct,
which is checked). No truth enters any run. After the run it prints, as a check that the code ran and nothing else, the number of detections of each catalogue and how many lie within 0.5 px of a truth
source of F606W < 23 (positions only: no magnitudes, no completeness, no purity); no setting is looked at or changed on the strength of that count.

    python -I experiments/exp_astro_baselines_repeat.py [--frames m2r1 m2r2 ...] [--out FILE]
"""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, HERE)
import exp_astro_baselines as AB      # noqa: E402  (numpy only at import; astropy, photutils and crowdsource are imported inside its runners)
import astro_repeat as R              # noqa: E402

REPO = os.path.abspath(os.path.join(HERE, '..'))
OUT = os.path.join(REPO, 'results', 'astro', 'baseline_catalogs_repeat.npz')
FREEZE = os.path.join(REPO, 'results', 'astro', 'freeze.json')
DEFAULT_PHOTUTILS = 'photutils'


def frozen_choices(path=FREEZE):
    """(photutils variant, crowdsource post-filter name) of the freeze record, 'tuned_baselines'."""
    with open(path) as f:
        tb = json.load(f)['tuned_baselines']
    return tb['photutils'], tb['crowdsource']


def parse_cs_name(name):
    """'crowdsource|snr5|ff0|qf0|spinf' -> (min flux / dflux, min fracflux, min qf, max |spread_model|) = (5.0, 0.0, 0.0, inf), the inverse of `exp_astro_m2.cs_name`."""
    p = name.split('|')
    if p[0] != 'crowdsource' or len(p) != 5:
        raise ValueError(name)
    sp = p[4][2:]
    return float(p[1][3:]), float(p[2][2:]), float(p[3][2:]), np.inf if sp == 'inf' else float(sp)


def patch_key(x0, y0):
    return 'x%d_y%d' % (x0, y0)


def crowdsource_filtered(b, x0, y0, f):
    """The crowdsource catalogue of a patch with the post-filter f = (snr, ff, qf, spread): the rule of `exp_astro_m2.cs_dets` (flux > 0 and flux >= snr x dflux, fracflux >= ff, qf >= qf, |spread| <= sp)."""
    k = patch_key(x0, y0) + '|crowdsource|'
    snr, ff, qf, sp = f
    m = (b[k + 'flux'] > 0) & (b[k + 'flux'] >= snr * b[k + 'dflux']) & (b[k + 'fracflux'] >= ff) & (b[k + 'qf'] >= qf) & (np.abs(b[k + 'spread']) <= sp)
    return b[k + 'x'][m], b[k + 'y'][m], b[k + 'flux'][m]


def frozen_catalogues(b, x0, y0, freeze=FREEZE):
    """The catalogues of a patch at the frozen settings: {'photutils (default)': (x, y, flux), 'photutils_ext_all_local_i6 (frozen)': ..., 'crowdsource (frozen post-filter)': ...} from the
    loaded catalogue file b (np.load of baseline_catalogs_repeat.npz)."""
    tag, cs = frozen_choices(freeze)
    k = patch_key(x0, y0)
    out = {'photutils (default)': tuple(b['%s|%s|%s' % (k, DEFAULT_PHOTUTILS, c)] for c in ('x', 'y', 'flux')), 'photutils (frozen: %s)' % tag: tuple(b['%s|%s|%s' % (k, tag, c)] for c in ('x', 'y', 'flux'))}
    out['crowdsource (frozen: %s)' % cs] = crowdsource_filtered(b, x0, y0, parse_cs_name(cs))
    return out


def run_patch(name, x0, y0, tag_tuned):
    """The default and the frozen photutils variants and crowdsource on one patch: ({key: array}, {method: (n sources, seconds)})."""
    z = np.load(R.patch_path(x0, y0, name))
    with open(os.path.join(os.path.dirname(R.patch_path(x0, y0, name)), 'meta.json')) as f:
        meta = json.load(f)
    dn, psf51 = z['dn'].astype(float), z['psf51'].astype(float)
    gain, dark, skysig = float(z['gain']), float(z['dark_variance']), float(meta['psfield_skysig'])
    spec = [d for d in AB.PHOT_EXT if d['tag'] == tag_tuned]
    if len(spec) != 1:
        raise SystemExit('the frozen photutils variant %s is not one of exp_astro_baselines.PHOT_EXT' % tag_tuned)
    out, info = {}, {}
    k = patch_key(x0, y0)
    cat, secs, _ = AB.run_photutils(dn, psf51, gain, dark, skysig, nsigma=5.0, maxiters=3, local_bkg=False)
    for c, v in cat.items():
        out['%s|%s|%s' % (k, DEFAULT_PHOTUTILS, c)] = v
    info[DEFAULT_PHOTUTILS] = (len(cat['x']), secs)
    cat, secs, _ = AB.run_photutils_ext(dn, psf51, gain, dark, skysig, spec[0])
    for c, v in cat.items():
        out['%s|%s|%s' % (k, tag_tuned, c)] = v
    info[tag_tuned] = (len(cat['x']), secs)
    cat, secs = AB.run_crowdsource(dn, float(z['sky_dn'].mean()), psf51, gain, dark)
    for c, v in cat.items():
        out['%s|crowdsource|%s' % (k, c)] = v
    info['crowdsource'] = (len(cat['x']), secs)
    return out, info


def positional_check(det_x, det_y, truth_patch, radius=0.5, mag_max=23.0):
    """Number of detections with a truth source of F606W < mag_max within `radius` px (positions only): a check that the catalogue is in the patch's coordinates, nothing more."""
    from scipy.spatial import cKDTree
    sel = truth_patch['mag'] < mag_max
    if not sel.any() or len(det_x) == 0:
        return 0
    d, _ = cKDTree(np.c_[truth_patch['x'][sel], truth_patch['y'][sel]]).query(np.c_[det_x, det_y], k=1)
    return int((d <= radius).sum())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--frames', nargs='*', default=None, help='names of the registry (default: all)')
    ap.add_argument('--out', default=OUT)
    A_ = ap.parse_args()
    import astro_common as AC
    tag_tuned, cs_name = frozen_choices()
    names = A_.frames or sorted(R.REPEAT)
    origins = [p for n in sorted(R.REPEAT) for p in R.REPEAT[n]['patches']]
    if len(set(origins)) != len(origins):
        raise SystemExit('two patches of the registry share an origin: the key x<x0>_y<y0> would collide')
    print('frozen settings (results/astro/freeze.json): photutils %s ; crowdsource post-filter %s = %s' % (tag_tuned, cs_name, parse_cs_name(cs_name)))
    out = {}
    for name in names:
        truth = R.truth_arrays(name)
        for (x0, y0), w in zip(R.REPEAT[name]['patches'], R.REPEAT[name]['windows']):
            cat, info = run_patch(name, x0, y0, tag_tuned)
            out.update(cat)
            tp = AC.truth_for_patch(truth, x0, y0)
            k = patch_key(x0, y0)
            print('%s patch (%d, %d) = window %s: %s' % (name, x0, y0, tuple(w), '; '.join('%s %d sources %.1f s' % (m, n, s) for m, (n, s) in info.items())), flush=True)
            fz = frozen_catalogues(out, x0, y0)
            print('    code check (positions only: detections within 0.5 px of a truth source of F606W < 23, of n detections): %s' % '; '.join(
                '%s: %d of %d' % (lab, positional_check(c[0], c[1], tp), len(c[0])) for lab, c in fz.items()), flush=True)
    os.makedirs(os.path.dirname(os.path.abspath(A_.out)), exist_ok=True)
    if A_.frames and os.path.exists(A_.out):                    # a partial run is merged into the existing file
        with np.load(A_.out) as z:
            old = {key: z[key] for key in z.files}
        old.update(out)
        out = old
    np.savez_compressed(A_.out, **out)
    print('wrote', A_.out, '(%d arrays)' % len(out))


if __name__ == '__main__':
    main()
