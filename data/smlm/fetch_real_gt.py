#!/usr/bin/env python
"""Resumable, checksummed download of REAL single-molecule localization recordings with a known (or well-defined) geometry.

These are the candidates of the pre-registered real-data test (``data/smlm/real_candidates.md``, ``results/smlm/real_gt_plan.md``). Raw
third-party data NEVER go into the repository: everything is fetched into $SMLM_REAL_GT (default /home/user/smlm_real_data) and the script refuses
to write inside the repository; only this script, ``checksums_real_gt.json`` and derived results are committed.

Sets (all Zenodo, licence CC BY 4.0; the file list, sizes and md5 below are what the Zenodo API returned on 2026-11-09)
-----------------------------------------------------------------------------------------------------------------------
heil_npc    https://doi.org/10.5281/zenodo.5068525   Heil et al., Light Sci Appl 7, 99 (2018): Xenopus laevis nuclear pore complexes, gp210-Alexa Fluor
            647, STORM on glass; 20000 raw frames of 117 x 118 px (16 bit), 102 nm pixel, EM gain 100, ThunderSTORM camera settings (offset 100 ADU,
            15.01 photoelectrons per ADU) in the protocol file; the ThunderSTORM localization table of the same movie. Files are fetched COMPLETE and
            verified against the md5 that Zenodo publishes. (The nano-mirror movie of the same record, 557 MB, is NOT fetched.)
nup96_au    https://doi.org/10.5281/zenodo.3540443   Thevathasan et al., Nat Methods 16, 1045 (2019): Nup96-SNAP-Alexa Fluor 647 nuclear pores, 2D STORM,
            106 nm pixel, 419 x 419 px frames (16 bit), two fields of view (17 and 18). A recording is stored in several 4.3 GB OME-TIFF pieces (about
            30000 frames per field of view); ONLY THE FIRST 400,000,000 BYTES OF PIECE 1 (the first 1056 frames) of each field are fetched with an HTTP
            Range request (the file is an uncompressed multi-page TIFF, so a prefix is a readable TIFF of whole frames; tifffile warns about the cut
            page). Also fetched complete: the SMAP localization tables of the full recordings (locs.csv, 30000 frames: the reference positions of the
            plan) and the wide-field images.
simflux     https://doi.org/10.5281/zenodo.7771932   Cnossen et al., Nat Methods 17, 59 (2020) (SIMFLUX data set): DNA-PAINT on 80 nm nanorulers, 400 x
            400 px frames (16 bit, Hamamatsu ORCA-Flash4.0 V2), patterned illumination. ONLY THE FIRST 300,000,000 BYTES (923 whole frames) of
            sim4_1_MMStack_Pos0.ome.tif (4.30 GB in all) are fetched.

Integrity
---------
* complete files: size equal to the Content-Length of the server, md5 equal to the one published by Zenodo (EXPECTED below, cross-checked with the
  API when online), sha256 equal to the record in ``checksums_real_gt.json``;
* prefixes: Zenodo publishes the md5 of the whole file only, so a prefix cannot be verified against it. The script checks that the server's
  Content-Range reports the expected TOTAL size of the file, that the prefix has the expected number of bytes and the sha256 recorded in
  ``checksums_real_gt.json`` (computed on the first verified download, 2026-11-09), and that it starts with a TIFF header.
A local file that does not match its record is NOT taken to mean that the server changed: a fresh copy is downloaded to a '.part' file and
compared (as in data/smlm/fetch.py). Nothing is extracted or executed.

    python data/smlm/fetch_real_gt.py                         # all sets
    python data/smlm/fetch_real_gt.py --sets heil_npc         # one set
    python data/smlm/fetch_real_gt.py --head-only             # only compare the Zenodo API / HEAD answers with the expected sizes and md5
    python data/smlm/fetch_real_gt.py --verify-only           # no download: verify what is on disk
    python data/smlm/fetch_real_gt.py --record                # (first run only) write checksums_real_gt.json from verified local files

Standard library only (TLS verification always on; proxy and CA bundle from the environment, e.g. SSL_CERT_FILE=/root/.ccr/ca-bundle.crt).
"""
import argparse
import hashlib
import http.client
import json
import os
import ssl
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
CHECKSUMS = os.path.join(HERE, 'checksums_real_gt.json')
DEFAULT_DATA = '/home/user/smlm_real_data'

# set -> list of dict(record, name, size (full file), md5 (published for the full file), prefix (bytes fetched; None = the whole file))
EXPECTED = {
    'heil_npc': [
        dict(record=5068525, name='NPCData_glass_STORM.tif', size=562040008, md5='1573237de6ab5ec7c9a0508251faa996', prefix=None),
        dict(record=5068525, name='NPCData_glass_STORM_loc.csv', size=86046259, md5='e24bc3f79a5ab5b51f0825256a2f53c2', prefix=None),
        dict(record=5068525, name='NPCData_glass_STORM_loc-protocol.txt', size=2692, md5='d97187a4d9227f051e92d341d7c8aad5', prefix=None),
        dict(record=5068525, name='NPCData_glass_EPI.tif', size=2163862, md5='fcc31f485bc90bb4bf27dba9d0be36c8', prefix=None),
    ],
    'nup96_au': [
        dict(record=3540443, name='Cos7_Nup96_A647_COT_Au_17_1_MMStack_Pos0.ome.tif', size=4302893758, md5='5422b25b653fc8343e966d3696c780ac', prefix=400000000),
        dict(record=3540443, name='Cos7_Nup96_A647_COT_Au_18_1_MMStack_Pos0.ome.tif', size=4302909026, md5='b9e6a5fa7ee9581ca31205ce40eec03b', prefix=400000000),
        dict(record=3540443, name='Cos7_Nup96_A647_COT_Au_17_1_locs.csv', size=72404532, md5='38a816e0302ef102c38d645e3b856e0f', prefix=None),
        dict(record=3540443, name='Cos7_Nup96_A647_COT_Au_18_1_locs.csv', size=71802563, md5='daa88e73d3fc548509530ee5578c4444', prefix=None),
        dict(record=3540443, name='Cos7_Nup96_A647_COT_Au_WF17_MMStack_Pos0.ome.tif', size=625365, md5='c5766164271c597c3f1fe31e33802379', prefix=None),
        dict(record=3540443, name='Cos7_Nup96_A647_COT_Au_WF18_MMStack_Pos0.ome.tif', size=625363, md5='a1a715109bc7925128699eb3d40ebd93', prefix=None),
    ],
    'simflux': [
        dict(record=7771932, name='sim4_1_MMStack_Pos0.ome.tif', size=4301948112, md5='6d8e849e2f25d139a9785ae1809be30a', prefix=300000000),
    ],
}
SETS = tuple(EXPECTED)
SUBDIR = {'heil_npc': 'zenodo_5068525', 'nup96_au': 'zenodo_3540443', 'simflux': 'zenodo_7771932'}


def local_name(item):
    """File name on disk: a prefix gets a suffix that says what it is, so that nobody mistakes it for the whole file."""
    if item['prefix'] is None:
        return item['name']
    return '%s.part0-%dMB' % (item['name'], item['prefix'] // 1000000)


def url_of(item):
    return 'https://zenodo.org/api/records/%d/files/%s/content' % (item['record'], urllib.parse.quote(item['name']))


def _open(url, method='GET', headers=None, timeout=120):
    req = urllib.request.Request(url, method=method, headers=headers or {})
    return urllib.request.urlopen(req, timeout=timeout, context=ssl.create_default_context())


def api_listing(record):
    """{file name: (size, md5)} from the Zenodo API."""
    with _open('https://zenodo.org/api/records/%d' % record) as r:
        d = json.loads(r.read().decode('utf-8'))
    out = {}
    for f in d['files']:
        out[f['key']] = (int(f['size']), f['checksum'].split(':', 1)[-1])
    return out


def digests(path, chunk=1 << 20):
    h1, h2 = hashlib.sha256(), hashlib.md5()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(chunk), b''):
            h1.update(b)
            h2.update(b)
    return h1.hexdigest(), h2.hexdigest()


def fetch_range(url, dest, nbytes, retries=8, log=print):
    """Resumable GET of the first ``nbytes`` bytes of url into dest + '.part' (Range: bytes=have-(nbytes-1)); returns (part path, total size of the file
    as announced by Content-Range, or None when the whole file was requested and no Content-Range was given)."""
    part, total = dest + '.part', None
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    for attempt in range(retries):
        have = os.path.getsize(part) if os.path.exists(part) else 0
        if have > nbytes:
            os.remove(part)
            have = 0
        if have == nbytes:
            break
        headers = {'Range': 'bytes=%d-%d' % (have, nbytes - 1)}
        try:
            t0, got = time.time(), 0
            with _open(url, 'GET', headers) as r:
                cr = r.headers.get('Content-Range', '')
                if r.status != 206 or not cr.startswith('bytes %d-' % have):
                    raise ConnectionError('server did not answer the Range request as expected (status %s, Content-Range %r)' % (r.status, cr))
                total = int(cr.rsplit('/', 1)[-1])
                with open(part, 'ab' if have else 'wb') as f:
                    while True:
                        b = r.read1(1 << 20)
                        if not b:
                            break
                        f.write(b)
                        got += len(b)
                        if time.time() - t0 > 5:
                            log('    %.1f / %.1f MB' % ((have + got) / 1e6, nbytes / 1e6))
                            t0 = time.time()
        except urllib.error.HTTPError as e:
            log('    HTTP %s, retrying' % e.code)
            time.sleep(min(2 * (attempt + 1), 10))
        except (urllib.error.URLError, http.client.HTTPException, ConnectionError, TimeoutError, OSError) as e:
            log('    transfer interrupted (%r), resuming' % (e,))
            time.sleep(min(2 * (attempt + 1), 10))
    size = os.path.getsize(part) if os.path.exists(part) else -1
    if size != nbytes:
        raise RuntimeError('incomplete download of %s: %d of %d bytes (re-run to resume)' % (url, size, nbytes))
    return part, total


def load_checksums():
    if os.path.exists(CHECKSUMS):
        with open(CHECKSUMS) as f:
            return json.load(f)
    return {}


def check_one(item, root, rec, log=print):
    """Verify a local file against EXPECTED and the record; returns (ok, message, sha256)."""
    path = os.path.join(root, SUBDIR_OF[id(item)], local_name(item))
    if not os.path.exists(path):
        return False, 'missing', None
    want = item['size'] if item['prefix'] is None else item['prefix']
    size = os.path.getsize(path)
    if size != want:
        return False, 'size %d, expected %d' % (size, want), None
    sha, md5 = digests(path)
    if item['prefix'] is None and md5 != item['md5']:
        return False, 'md5 %s differs from the one published by Zenodo (%s)' % (md5, item['md5']), sha
    with open(path, 'rb') as f:
        head = f.read(4)
    if item['name'].endswith('.tif') and head[:2] not in (b'II', b'MM'):
        return False, 'does not start with a TIFF header', sha
    r = rec.get(item['name'])
    if r is not None and r['sha256'] != sha:
        return False, 'sha256 %s differs from the record %s' % (sha, r['sha256']), sha
    return True, 'ok (%s)' % ('complete, md5 equals the Zenodo md5' if item['prefix'] is None else 'prefix of %d bytes, sha256 equals the record' % want), sha


SUBDIR_OF = {}
for _s, _items in EXPECTED.items():
    for _it in _items:
        SUBDIR_OF[id(_it)] = SUBDIR[_s]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--sets', nargs='*', default=list(SETS), choices=SETS)
    ap.add_argument('--head-only', action='store_true')
    ap.add_argument('--verify-only', action='store_true')
    ap.add_argument('--record', action='store_true', help='write checksums_real_gt.json from the (verified) local files')
    ap.add_argument('--data', default=os.environ.get('SMLM_REAL_GT', DEFAULT_DATA))
    a = ap.parse_args(argv)
    root = os.path.abspath(a.data)
    if os.path.commonpath([root, REPO]) == REPO:
        sys.exit('refusing to write raw data inside the repository (%s); set --data or $SMLM_REAL_GT' % root)
    rec = load_checksums().get('files', {})
    new = dict(rec)
    bad = 0
    for s in a.sets:
        print('== %s -> %s' % (s, os.path.join(root, SUBDIR[s])))
        api = {}
        for r in sorted({it['record'] for it in EXPECTED[s]}):
            try:
                api.update(api_listing(r))
            except Exception as e:                                         # offline: the hard-coded values are used
                print('   (Zenodo API unreachable: %r; the expected values of this script are used)' % (e,))
        for it in EXPECTED[s]:
            if it['name'] in api and api[it['name']] != (it['size'], it['md5']):
                sys.exit('Zenodo now lists %s as %s, this script expects %s: the record changed' % (it['name'], api[it['name']], (it['size'], it['md5'])))
            dest = os.path.join(root, SUBDIR[s], local_name(it))
            want = it['size'] if it['prefix'] is None else it['prefix']
            print(' - %s (%.1f MB%s)' % (it['name'], want / 1e6, '' if it['prefix'] is None else ' = prefix of %.2f GB' % (it['size'] / 1e9)))
            if a.head_only:
                continue
            ok, msg, sha = check_one(it, root, rec) if os.path.exists(dest) else (False, 'missing', None)
            if not ok and not a.verify_only:
                if os.path.exists(dest) and msg != 'missing':
                    print('   local copy does not verify (%s); downloading again' % msg)
                    os.replace(dest, dest + '.suspect')
                part, total = fetch_range(url_of(it), dest, want)
                if total is not None and total != it['size']:
                    sys.exit('server reports a total size of %d for %s, expected %d' % (total, it['name'], it['size']))
                os.replace(part, dest)
                ok, msg, sha = check_one(it, root, rec)
                if ok and os.path.exists(dest + '.suspect'):
                    os.remove(dest + '.suspect')
            print('   %s' % msg)
            if ok and a.record:
                new[it['name']] = dict(set=s, record='https://doi.org/10.5281/zenodo.%d' % it['record'], local_name=local_name(it), bytes=want,
                                       full_file_bytes=it['size'], full_file_md5_published_by_zenodo=it['md5'], prefix_only=it['prefix'] is not None,
                                       sha256=sha, url=url_of(it))
            bad += (not ok)
    if a.record and not bad:
        with open(CHECKSUMS, 'w') as f:
            json.dump({'recorded_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), 'files': new}, f, indent=1, sort_keys=True)
            f.write('\n')
        print('recorded %d files in %s' % (len(new), CHECKSUMS))
    sys.exit(1 if bad else 0)


if __name__ == '__main__':
    main()
