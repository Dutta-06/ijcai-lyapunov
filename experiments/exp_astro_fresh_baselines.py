"""
Catalogues of the baselines on the fresh-cluster patches, at their frozen settings and with the ISS's front end as selected in round 2 (needs the DOMAIN environment: astropy, photutils, crowdsource)
=============================================================================================================================================================================================

The patches are those of `results/astro/fresh_patches.json` (`astro_fresh.py`). This script writes, once, the catalogues of

  photutils     the default of `exp_astro_baselines.run_photutils` (tag 'photutils') and the variant that `results/astro/freeze.json` names as tuned ('photutils_ext_all_local_i6',
                `exp_astro_baselines.run_photutils_ext` with its entry of `PHOT_EXT`)
  crowdsource   `exp_astro_baselines.run_crowdsource` (fit_im, unfiltered, tag 'crowdsource', with dflux, qf, fracflux, spread and rchi2); the tuned post-filter of `freeze.json`
                ('crowdsource|snr5|ff0|qf0|spinf') is applied by the reader (`exp_astro_fresh.cs_filtered`)
  front end     the two variants SELECTED in round 2 of the baselines given the ISS's front end (`baselines_hp2.md`: photutils `phu2_hp3_all_i6_b1`, crowdsource `cs2_hp2_t3_pm`), read from the
                verified selection of `baselines_hp2_raw.npz` (`exp_astro_hp2.verify_selection2`), run by `exp_astro_baselines_hp2.run_variant2` unchanged. NOTHING IS RE-SELECTED.

Every function of the earlier baseline scripts is imported, not copied (they are hash-frozen or part of the freeze record). No truth is read by this script (it imports none and compares nothing); it
prints the number of detections and the CPU seconds of each catalogue. Keys: '<cluster>|x<x0>_y<y0>|<tag>|<field>' (x, y, flux; crowdsource also dflux, qf, fracflux, spread, rchi2), positions in the
patch's 0-based pixel coordinates (column, row), flux in DN; '<tag>|secs' = CPU seconds; 'meta|photutils_tuned', 'meta|crowdsource_tuned', 'meta|photutils_fe', 'meta|crowdsource_fe' = the tags used.

    python -I experiments/exp_astro_fresh_baselines.py [--out FILE]      (refuses to overwrite an existing file)
"""
import argparse
import json
import os
import sys
import time

os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np

HERE = os.path.abspath(os.path.dirname(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)
import astro_fresh as F                       # noqa: E402
import exp_astro_baselines as AB              # noqa: E402
import exp_astro_baselines_hp2 as G2          # noqa: E402
import exp_astro_baselines_repeat as BR       # noqa: E402  (parse_cs_name, frozen_choices)
import exp_astro_hp2 as H2                    # noqa: E402
import exp_astro_m2 as X                      # noqa: E402

OUT = os.environ.get('FRESH_BASE', os.path.join(REPO, 'results', 'astro', 'fresh_baseline_catalogs.npz'))
FREEZE = os.path.join(REPO, 'results', 'astro', 'freeze.json')
DEFAULT_PHOTUTILS = 'photutils'


def pkey(cluster, x0, y0):
    return '%s|x%d_y%d' % (cluster, x0, y0)


def tags():
    """(tuned photutils tag, tuned crowdsource name, front-end photutils tag, front-end crowdsource tag): the first two from freeze.json, the last two from the verified round-2 selection."""
    tp, tc = BR.frozen_choices(FREEZE)
    sel = H2.verify_selection2(X.load_raw(H2.RAW2))
    return tp, tc, sel['photutils'], sel['crowdsource']


def run_patch(cluster, x0, y0, entry, tp, fe_p, fe_c, root=None):
    """The five catalogue families of one patch: ({key: array}, {tag: (n sources, seconds)})."""
    z = np.load(F.patch_path(x0, y0, entry, root))
    with open(os.path.join(os.path.dirname(F.patch_path(x0, y0, entry, root)), 'meta.json')) as f:
        meta = json.load(f)
    dn, psf51 = z['dn'].astype(float), z['psf51'].astype(float)
    gain, dark, skysig = float(z['gain']), float(z['dark_variance']), float(meta['psfield_skysig'])
    spec = [d for d in AB.PHOT_EXT if d['tag'] == tp]
    if len(spec) != 1:
        raise SystemExit('the frozen photutils variant %s is not one of exp_astro_baselines.PHOT_EXT' % tp)
    out, info = {}, {}
    k = pkey(cluster, x0, y0)

    def put(tag, cat, secs):
        for c, v in cat.items():
            out['%s|%s|%s' % (k, tag, c)] = np.asarray(v)
        out['%s|%s|secs' % (k, tag)] = np.array([secs])
        info[tag] = (len(cat['x']), secs)
    cat, secs, _ = AB.run_photutils(dn, psf51, gain, dark, skysig, nsigma=5.0, maxiters=3, local_bkg=False)
    put(DEFAULT_PHOTUTILS, cat, secs)
    cat, secs, _ = AB.run_photutils_ext(dn, psf51, gain, dark, skysig, spec[0])
    put(tp, cat, secs)
    cat, secs = AB.run_crowdsource(dn, float(z['sky_dn'].mean()), psf51, gain, dark)
    put('crowdsource', cat, secs)
    a = dict(dn=dn, psf51=psf51, gain=gain, dark=dark, sky=float(z['sky_dn'].mean()))
    for tag in (fe_p, fe_c):
        cat, secs = G2.run_variant2(tag, a)
        put(tag, cat, secs)
    return out, info


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', default=OUT)
    A_ = ap.parse_args()
    if os.path.exists(A_.out):
        raise SystemExit('%s exists: a catalogue file is never overwritten' % A_.out)
    tp, tc, fe_p, fe_c = tags()
    print('frozen settings: photutils tuned %s ; crowdsource tuned %s ; front end (round 2 selection): photutils %s, crowdsource %s' % (tp, tc, fe_p, fe_c), flush=True)
    reg = F.load_registry()
    out = {'meta|photutils_tuned': np.array([tp]), 'meta|crowdsource_tuned': np.array([tc]), 'meta|photutils_fe': np.array([fe_p]), 'meta|crowdsource_fe': np.array([fe_c])}
    for cluster in sorted(reg):
        for x0, y0 in reg[cluster]['patches']:
            t0 = time.time()
            cat, info = run_patch(cluster, x0, y0, reg[cluster], tp, fe_p, fe_c)
            out.update(cat)
            print('%s patch (%d, %d): %s  [wall %.0f s]' % (cluster, x0, y0, '; '.join('%s %d sources %.1f s' % (m, n, s) for m, (n, s) in info.items()), time.time() - t0), flush=True)
    os.makedirs(os.path.dirname(os.path.abspath(A_.out)), exist_ok=True)
    np.savez_compressed(A_.out, **out)
    print('wrote', A_.out, '(%d arrays)' % len(out))


if __name__ == '__main__':
    main()
