"""
Which 100 x 100 pixel patches of the M2 frame are fully inside the HST (HUGS) footprint? A survey from TRUTH GEOMETRY alone
================================================================================================================================

Referee round IS1a of the astronomy experiment found that two of the four test patches first chosen, (410, 630) and (410, 730), lay largely outside the footprint of the truth
catalogue (a block of 20 x 20 pixels with no truth source brighter than 22 mag while its pixels hold stars: every detection there is a false positive by construction). The
replacements were chosen by this survey before any test-patch detection existed. The rule (declared here, applied to truth positions only):

  * the patch origin (x0, y0) lies on a 20-pixel grid with 60 <= x0, x0 + 100 <= 1900, 60 <= y0, y0 + 100 <= 1380 (away from the frame's edges);
  * every one of the 25 blocks of 20 x 20 pixels holds at least MIN_BLOCK = 50 truth sources of magnitude < 24 ("fully covered");
  * no pixel overlap (and a margin of 5 pixels) with the development patches and the other test patches;
  * crowding levels spanning the benchmark's: the number of truth sources brighter than 22 mag in the scored area, 923 to 2057 for the patches chosen: the benchmark's 969, a denser one
    (2057) and two between and below (923, 1446).

    python experiments/astro_patch_survey.py [--min-block 50]          prints the candidates sorted by scored density (needs $ASTRO_DATA/truth_hugs_meth1.npz)
"""
import argparse
import os
import sys

import numpy as np

HERE = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, HERE)
import astro_common as A      # noqa: E402

ROOT = os.environ.get('ASTRO_DATA', '/home/user/astro_data')
MIN_BLOCK = 50
FRAME_X, FRAME_Y = (60, 1900), (60, 1380)
USED = [(110, 630), (210, 630), (310, 630), (310, 730)]


def coverage(truth, x0, y0):
    """(smallest number of truth sources of m < 24 over the 25 blocks, number of scored sources brighter than 22 mag) of the patch with origin (x0, y0)."""
    t = A.truth_for_patch(truth, x0, y0, mag_max=24.0)
    if len(t['x']) < 25:
        return 0, 0
    cnt = np.bincount(A.block_of(np.c_[t['x'], t['y']]), minlength=25)
    return int(cnt.min()), int(((t['mag'] < 22.0) & t['inner']).sum())


def overlaps(p, others, margin=5):
    return any(abs(p[0] - o[0]) < 100 + margin and abs(p[1] - o[1]) < 100 + margin for o in others)


def survey(truth, min_block=MIN_BLOCK, used=USED, step=20):
    """[(x0, y0, smallest block count, scored sources < 22 mag)] of the fully covered candidates, sorted by scored density."""
    rows = []
    for x0 in range(FRAME_X[0], FRAME_X[1] - 100 + 1, step):
        for y0 in range(FRAME_Y[0], FRAME_Y[1] - 100 + 1, step):
            if overlaps((x0, y0), used):
                continue
            mb, n22 = coverage(truth, x0, y0)
            if mb >= min_block:
                rows.append((x0, y0, mb, n22))
    return sorted(rows, key=lambda r: (r[3], r[0], r[1]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--min-block', type=int, default=MIN_BLOCK)
    a = ap.parse_args()
    z = np.load(os.path.join(ROOT, 'truth_hugs_meth1.npz'))
    truth = dict(px=z['px'], py=z['py'], m606=z['m606'])
    for p in USED + [(360, 1040), (100, 1100)]:
        mb, n22 = coverage(truth, *p)
        print('patch (%4d, %4d): smallest block count %3d, scored sources < 22 mag %4d' % (p[0], p[1], mb, n22))
    rows = survey(truth, a.min_block)
    print('%d fully covered candidates (every block >= %d truth sources of m < 24), not overlapping the development and the first two test patches' % (len(rows), a.min_block))
    for r in rows[::max(len(rows) // 25, 1)]:
        print('  x0 %4d y0 %4d  smallest block %3d  scored < 22 mag %5d' % r)


if __name__ == '__main__':
    main()
