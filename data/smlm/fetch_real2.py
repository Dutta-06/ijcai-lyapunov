#!/usr/bin/env python
"""Resumable, checksummed download of the FRESH real recording of the second pore study (`results/smlm/real2_plan.md`).

One additional field of view of the Nup96 record already used by `fetch_real_gt.py` (Zenodo 3540443, Thevathasan et al., Nat Methods 16, 1045
(2019), CC BY 4.0): field 19, the first 400,000,000 bytes of piece 1 (the first 1055 whole frames of 419 x 419 px, an uncompressed multi-page TIFF
read with an HTTP Range request) and the SMAP localization table of the whole recording (55 MB, complete). Field 19 was never opened before the
second study (the first study used fields 17 and 18).

Raw third-party data never go into the repository: they are fetched into $SMLM_REAL_GT (default /home/user/smlm_real_data); the script refuses to
write inside the repository. Integrity as in `fetch_real_gt.py`: the table is checked against the md5 published by Zenodo; the prefix against the total
size announced in Content-Range, its byte count and the sha256 recorded in `checksums_real2.json`.

    python data/smlm/fetch_real2.py [--head-only | --verify-only | --record]
"""
import argparse
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fetch_real_gt as F                                                 # noqa: E402

CHECKSUMS = os.path.join(HERE, 'checksums_real2.json')
EXPECTED = [
    dict(record=3540443, name='Cos7_Nup96_A647_COT_Au_19_1_MMStack_Pos0.ome.tif', size=4303400016, md5='91c0877bd819217a49f1cd407580dddd', prefix=400000000),
    dict(record=3540443, name='Cos7_Nup96_A647_COT_Au_19_1_locs.csv', size=55176718, md5='556c8cc46abe603a6510cf04731cdb3e', prefix=None),
]
SUBDIR = 'zenodo_3540443'


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--head-only', action='store_true')
    ap.add_argument('--verify-only', action='store_true')
    ap.add_argument('--record', action='store_true')
    ap.add_argument('--data', default=os.environ.get('SMLM_REAL_GT', F.DEFAULT_DATA))
    a = ap.parse_args(argv)
    root = os.path.abspath(a.data)
    if os.path.commonpath([root, F.REPO]) == F.REPO:
        sys.exit('refusing to write raw data inside the repository (%s)' % root)
    rec = json.load(open(CHECKSUMS))['files'] if os.path.exists(CHECKSUMS) else {}
    new, bad = dict(rec), 0
    try:
        api = F.api_listing(3540443)
    except Exception as e:
        api = {}
        print('(Zenodo API unreachable: %r; the expected values are used)' % (e,))
    for it in EXPECTED:
        F.SUBDIR_OF[id(it)] = SUBDIR
        if it['name'] in api and api[it['name']] != (it['size'], it['md5']):
            sys.exit('Zenodo now lists %s as %s, expected %s: the record changed' % (it['name'], api[it['name']], (it['size'], it['md5'])))
        dest = os.path.join(root, SUBDIR, F.local_name(it))
        want = it['size'] if it['prefix'] is None else it['prefix']
        print(' - %s (%.1f MB)' % (it['name'], want / 1e6))
        if a.head_only:
            continue
        ok, msg, sha = F.check_one(it, root, rec) if os.path.exists(dest) else (False, 'missing', None)
        if not ok and not a.verify_only:
            part, total = F.fetch_range(F.url_of(it), dest, want)
            if total is not None and total != it['size']:
                sys.exit('server reports a total size of %d for %s, expected %d' % (total, it['name'], it['size']))
            os.replace(part, dest)
            ok, msg, sha = F.check_one(it, root, rec)
        print('   %s' % msg)
        if ok and a.record:
            new[it['name']] = dict(record='https://doi.org/10.5281/zenodo.%d' % it['record'], local_name=F.local_name(it), bytes=want, full_file_bytes=it['size'],
                                   full_file_md5_published_by_zenodo=it['md5'], prefix_only=it['prefix'] is not None, sha256=sha, url=F.url_of(it))
        bad += (not ok)
    if a.record and not bad:
        with open(CHECKSUMS, 'w') as f:
            json.dump({'recorded_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), 'files': new}, f, indent=1, sort_keys=True)
            f.write('\n')
    sys.exit(1 if bad else 0)


if __name__ == '__main__':
    main()
