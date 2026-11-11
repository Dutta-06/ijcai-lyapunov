"""
Catalogues of the baselines on the held-out patches of the cluster-dependence test (needs the DOMAIN environment: astropy, photutils, crowdsource)
===================================================================================================================================================

The patches are those of `results/astro/clusterdep_patches.json` (`clusterdep_select.py build`). As `exp_astro_fresh_baselines.py` (whose `run_patch` and `tags` it imports, unchanged): photutils default and as
tuned on the M2 development patches (`freeze.json`), crowdsource default (the tuned post-filter is applied by the reader) and the two variants selected in round 2 of the baselines given the ISS's front end
(`baselines_hp2.md`: photutils `phu2_hp3_all_i6_b1`, crowdsource `cs2_hp2_t3_pm`). NOTHING IS RE-SELECTED. No truth is read (this script imports none and compares nothing); it prints the number of detections
and the CPU seconds of each catalogue. The file is written once, before the ISS run and before anything is scored.

    python -I experiments/exp_astro_clusterdep_baselines.py [--out FILE]      (refuses to overwrite an existing file)
"""
import argparse
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
import exp_astro_fresh_baselines as FB        # noqa: E402

OUT = os.environ.get('CLUSTERDEP_BASE', os.path.join(REPO, 'results', 'astro', 'clusterdep_baseline_catalogs.npz'))
REGISTRY = os.environ.get('CLUSTERDEP_REGISTRY', os.path.join(REPO, 'results', 'astro', 'clusterdep_patches.json'))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', default=OUT)
    a = ap.parse_args()
    if os.path.exists(a.out):
        raise SystemExit('%s exists: a catalogue file is never overwritten' % a.out)
    tp, tc, fe_p, fe_c = FB.tags()
    print('frozen settings: photutils tuned %s ; crowdsource tuned %s ; front end (round 2 selection): photutils %s, crowdsource %s' % (tp, tc, fe_p, fe_c), flush=True)
    reg = F.load_registry(REGISTRY)
    out = {'meta|photutils_tuned': np.array([tp]), 'meta|crowdsource_tuned': np.array([tc]), 'meta|photutils_fe': np.array([fe_p]), 'meta|crowdsource_fe': np.array([fe_c])}
    for cluster in sorted(reg):
        for x0, y0 in reg[cluster]['patches']:
            t0 = time.time()
            cat, info = FB.run_patch(cluster, x0, y0, reg[cluster], tp, fe_p, fe_c)
            out.update(cat)
            print('%s patch (%d, %d): %s  [wall %.0f s]' % (cluster, x0, y0, '; '.join('%s %d sources %.1f s' % (m, n, s) for m, (n, s) in info.items()), time.time() - t0), flush=True)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    np.savez_compressed(a.out, **out)
    print('wrote', a.out, '(%d arrays)' % len(out))


if __name__ == '__main__':
    main()
