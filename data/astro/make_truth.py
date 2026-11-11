#!/usr/bin/env python
"""The HST truth list for the M2 benchmark in the SDSS frame's pixel coordinates (needs numpy, pandas, astropy: the domain environment, not the repository's).

Reads the HUGS NGC 7089 catalogue (method 1; `fetch.py fetch --group hst`) and the SDSS r frame (`fetch.py fetch --group core`), converts every source with a valid
F606W magnitude from (RA, Dec) to the frame's 0-based pixel coordinates (pixel i centred at integer i, column = x, row = y) with the frame's own WCS — no further
shift: the image-based alignment of `results/astro/data_access.md` is 0.10 px in x and 0.02 px in y, far below the 0.5 px matching tolerance — and writes
$ASTRO_DATA/truth_hugs_meth1.npz with the arrays px, py, m606 (F606W, used as r: the 2017 benchmark's convention), e606 (its rms), nf606 (exposures it is found in).

    python -I data/astro/make_truth.py
"""
import os
import sys

import numpy as np

ROOT = os.environ.get('ASTRO_DATA', '/home/user/astro_data')
CAT = os.path.join(ROOT, 'hugs_ngc7089_meth1', 'hlsp_hugs_hst_wfc3-uvis-acs-wfc_ngc7089_multi_v1_catalog-meth1.txt')
FRAME = os.path.join(ROOT, 'sdss_frame_r', 'frame-r-002583-2-0136.fits.bz2')
OUT = os.path.join(ROOT, 'truth_hugs_meth1.npz')
COL = dict(m606=20, e606=21, nf606=24, ra=33, dec=34)           # 0-based columns of the numeric part (the string identifier, column 35, removed)


def main():
    try:
        import pandas as pd
        from astropy.io import fits
        from astropy.wcs import WCS
    except ImportError as e:
        sys.exit('needs pandas and astropy (%s): run it in the domain environment' % e)
    for p in (CAT, FRAME):
        if not os.path.exists(p):
            sys.exit('missing %s: run data/astro/fetch.py first' % p)
    df = pd.read_csv(CAT, sep=r'\s+', comment='#', header=None, dtype={35: str})
    a = df.drop(columns=[35]).to_numpy(dtype=float)
    m = a[:, COL['m606']]
    ok = (m > -90) & (m < 90)
    with fits.open(FRAME, memmap=False) as h:
        w = WCS(h[0].header)
    px, py = w.all_world2pix(a[ok, COL['ra']], a[ok, COL['dec']], 0)
    np.savez(OUT, px=px, py=py, m606=m[ok], e606=a[ok, COL['e606']], nf606=a[ok, COL['nf606']])
    print('wrote %s: %d sources with a valid F606W of %d rows' % (OUT, int(ok.sum()), len(a)))


if __name__ == '__main__':
    main()
