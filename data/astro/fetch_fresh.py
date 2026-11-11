#!/usr/bin/env python3
"""Data access for the fresh-cluster test (two globular clusters that are not M2, M53 or M92): psField, frame and checksums.

Companion of `data/astro/fetch_second.py`, whose helpers it imports (`download`, `check_or_record`, `sha256_file`, the URL builders); only the checksum file
(`data/astro/checksums_fresh.json`) and the download root (`$ASTRO_DATA_4`, default /home/user/astro_data_4, never inside the repository) differ. The HUGS catalogues and the
photoField files of the survey runs are the ones that `fetch_second.py` already downloaded into `$ASTRO_DATA_2` and recorded in `checksums_second.json`; they are only
re-hashed here (`verify-shared`). Standard library only; the network hosts are data.sdss.org (SDSS SAS dr17) and nothing else; a refusal (401/403/404) is reported and never worked around.

Sub-commands
  list                         the files this script manages and whether they are on disk
  psfield   [--frames C:run:camcol:field ...]   the psField (0.4 MB) of every frame of the selection rule (`astro_fresh.py select` reads them); default: the frames of `FRAMES`
  frame     --frames ...       the r-band frame (3 MB) of the frames named (the selected clusters)
  verify                       re-hash the files of this script against checksums_fresh.json
  verify-shared                re-hash the shared HUGS and photoField files against checksums_second.json
"""
import argparse
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location('fetch_second', os.path.join(HERE, 'fetch_second.py'))
FS = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(FS)

CHECKSUMS = os.path.join(HERE, 'checksums_fresh.json')
DEFAULT_ROOT = '/home/user/astro_data_4'
FS.CHECKSUMS = CHECKSUMS                                              # the helpers of fetch_second read and write this module global
FS.root = lambda: os.environ.get('ASTRO_DATA_4', DEFAULT_ROOT)

# cluster -> (run, camcol, field) of the best frame of the coverage survey (`astro_second.py survey`, stored in results/astro/fresh_survey.txt); filled in by the rule's input, see `astro_fresh.py`
FRAMES = {
    'ngc5053': (5390, 5, 168), 'ngc5272': (4646, 3, 81), 'ngc5466': (4623, 1, 324), 'ngc5904': (1458, 4, 700), 'ngc6205': (3226, 5, 127), 'ngc6254': (5391, 5, 185),
    'ngc6535': (1463, 3, 61), 'ngc6838': (6895, 3, 52), 'ngc6934': (1738, 1, 29), 'ngc7078': (1749, 6, 141)}
HUGS_SHARED = ('ngc5053', 'ngc5272', 'ngc5466', 'ngc5904', 'ngc6205', 'ngc6254', 'ngc6535', 'ngc6838', 'ngc6934', 'ngc6981', 'ngc7078')


def items_for(frames, kinds=('psfield',)):
    """key -> (group, directory under the root, url, filename) for the psField and/or frame of each (cluster, run, camcol, field) of the list `frames`."""
    it = {}
    for c, run, cam, fld in frames:
        if 'psfield' in kinds:
            k, u, f = FS._psfield(run, cam, fld)
            it[k] = ('psfield', c, u, f)
        if 'frame' in kinds:
            k, u, f = FS._frame(run, cam, fld)
            it[k] = ('frame', c, u, f)
    return it


def parse_frames(specs):
    out = []
    for s in specs:
        c, run, cam, fld = s.split(':')
        out.append((c, int(run), int(cam), int(fld)))
    return out


def fetch(items):
    sums = FS.load_sums()
    for key, (g, d, url, fn) in items.items():
        dest = os.path.join(FS.root(), d, key, fn)
        print('[%s]' % key)
        lastmod = None
        if os.path.exists(dest):
            print('  already on disk; verifying')
        else:
            lastmod = FS.download(url, dest)
        FS.check_or_record(items, key, dest, sums, lastmod)


def verify(items, sums_path=CHECKSUMS, root=None):
    import json
    with open(sums_path) as f:
        sums = json.load(f)
    rc = 0
    for key, (g, d, url, fn) in items.items():
        p = os.path.join(root or FS.root(), d, key, fn)
        if not os.path.exists(p):
            print('[%s] missing: %s' % (key, p))
            rc = 1
            continue
        digest, size = FS.sha256_file(p)
        rec = sums.get('files', {}).get(key, {})
        ok = digest == rec.get('sha256') and size == rec.get('bytes')
        print('[%s] %s  %s  %d bytes' % (key, 'OK ' if ok else 'BAD', digest[:16], size))
        rc |= (not ok)
    return rc


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('cmd', choices=['list', 'psfield', 'frame', 'verify', 'verify-shared'])
    ap.add_argument('--frames', nargs='*', help='cluster:run:camcol:field, for example ngc5272:4646:3:81')
    a = ap.parse_args()
    frames = parse_frames(a.frames) if a.frames else [(c,) + f for c, f in FRAMES.items()]
    if a.cmd == 'list':
        for k, (g, d, u, f) in items_for(frames, ('psfield', 'frame')).items():
            print('%-34s %-8s %s %s' % (k, g, 'present' if os.path.exists(os.path.join(FS.root(), d, k, f)) else 'absent ', u))
    elif a.cmd in ('psfield', 'frame'):
        fetch(items_for(frames, (a.cmd,)))
    elif a.cmd == 'verify':
        sys.exit(verify(items_for(frames, ('psfield', 'frame'))))
    else:
        items = {'hugs_' + c: ('hugs', c, '', 'hlsp_hugs_hst_wfc3-uvis-acs-wfc_%s_multi_v1_catalog-meth1.txt' % c) for c in HUGS_SHARED}
        sys.exit(verify(items, FS.CHECKSUMS.replace('checksums_fresh', 'checksums_second'), os.environ.get('ASTRO_DATA_2', FS.DEFAULT_ROOT)))


if __name__ == '__main__':
    main()
