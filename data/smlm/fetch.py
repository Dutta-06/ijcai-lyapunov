#!/usr/bin/env python
"""Resumable download of the 2D simulated single-molecule localization benchmark data.

The four 2D microtubule training sets of the 2016 localization-microscopy challenge have a public ground truth:

    MT0.N1.LD   19,996 frames, low density (0.2), high SNR
    MT0.N1.HD    2,500 frames, high density (2),  high SNR
    MT0.N2.LD   19,996 frames, low density,       low SNR
    MT0.N2.HD    2,500 frames, high density,      low SNR

For each requested set this script fetches (into $SMLM_DATA, default /home/user/smlm_data, NEVER into the repository)

    <set>/sequence-<set>-2D-Exp-as-stack.zip   one 16-bit TIFF stack of n 64 x 64 images (the 2D modality)
    <set>/activations.csv                      the ground truth: one row per emitter ON event (frame, x, y, z, photons)
    <set>/positions.csv                        the fluorophore sites (x, y, z in nm; no frame information)
    <set>/stats_activation_sumary.csv          the simulator's own summary (number of activations, ...)

Every file is checked against the server's Content-Length (HEAD request) and, where known, the size seen when the URLs were
verified; the sha256 of every file is recorded in checksums.json next to this script. A later run re-verifies the recorded
digest. A local file that no longer matches its record (a flipped byte keeps the size) is NOT taken to mean that the server
changed: the file is downloaded again to a '.part' file and compared. If the fresh copy matches the record the local file was
damaged and is replaced (repaired); only when the fresh copy does not match either does the script stop with "server content
changed", and --accept-new then records the digest of that FRESH download (never the digest of a file that was already on disk).
A downloaded zip that fails its CRC check is removed and nothing is recorded. Downloads resume from a '.part' file with a Range
request. Nothing is ever extracted or executed: downstream code reads the zip member in memory.

    python data/smlm/fetch.py                    # MT0.N1.LD stack + ground truth
    python data/smlm/fetch.py --sets MT0.N1.HD   # additionally / instead: the other sets
    python data/smlm/fetch.py --head-only        # only check URLs and sizes

Only the Python standard library is used (TLS verification is always on; the proxy and CA bundle come from the environment).
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
import urllib.request
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
CHECKSUMS = os.path.join(HERE, 'checksums.json')
DEFAULT_DATA = '/home/user/smlm_data'

ROOT = 'https://srm.epfl.ch/srm'
SETS = ('MT0.N1.LD', 'MT0.N1.HD', 'MT0.N2.LD', 'MT0.N2.HD')

# Sizes in bytes as returned by HEAD when the URLs were verified (2026-11-07); None = unknown, checked against HEAD only.
EXPECTED_BYTES = {
    'MT0.N1.LD': {'stack': 114729551, 'activations': 1317983, 'positions': 296172, 'summary': 665},
    'MT0.N1.HD': {'stack': 14921627, 'activations': 1670528, 'positions': 379001, 'summary': 667},
    'MT0.N2.LD': {'stack': 90625319, 'activations': 1294643, 'positions': 296172, 'summary': 662},
    'MT0.N2.HD': {'stack': 11851051, 'activations': 1640924, 'positions': 379001, 'summary': 664},
}


def files_for(name):
    """kind -> (relative local path, url) for one dataset."""
    stack = 'sequence-%s-2D-Exp-as-stack.zip' % name
    sample = '%s/dataset/challenge-3D-simulation/%s/sample/' % (ROOT, name)
    return {
        'stack': ('%s/%s' % (name, stack), '%s/Data/challenge-3D-simulation/%s/Data/%s' % (ROOT, name, stack)),
        'activations': ('%s/activations.csv' % name, sample + 'activations.csv'),
        'positions': ('%s/positions.csv' % name, sample + 'positions.csv'),
        'summary': ('%s/stats_activation_sumary.csv' % name, sample + 'stats_activation_sumary.csv'),
    }


def _open(url, method='GET', headers=None, timeout=60):
    req = urllib.request.Request(url, method=method, headers=headers or {})
    return urllib.request.urlopen(req, timeout=timeout, context=ssl.create_default_context())


def head(url, retries=4):
    """(content_length, last_modified, etag) from a HEAD request."""
    last = None
    for k in range(retries):
        try:
            with _open(url, 'HEAD') as r:
                n = r.headers.get('Content-Length')
                if n is None:
                    raise RuntimeError('no Content-Length for %s' % url)
                return int(n), r.headers.get('Last-Modified'), r.headers.get('ETag')
        except (urllib.error.URLError, http.client.HTTPException, ConnectionError, TimeoutError) as e:
            last = e
            time.sleep(2 * (k + 1))
    raise RuntimeError('HEAD failed for %s: %r' % (url, last))


def sha256_of(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(chunk), b''):
            h.update(b)
    return h.hexdigest()


def download(url, dest, total, retries=8, log=print, finalize=True):
    """Resumable GET of url into dest (via dest + '.part'); total = the size announced by HEAD. With ``finalize=False`` the complete
    '.part' file is left where it is (to be compared before it replaces anything) and its path is returned."""
    part = dest + '.part'
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    for attempt in range(retries):
        have = os.path.getsize(part) if os.path.exists(part) else 0
        if have > total:                                  # corrupt leftovers: start again
            os.remove(part)
            have = 0
        if have == total:
            break
        headers = {'Range': 'bytes=%d-' % have} if have else {}
        try:
            t0, got = time.time(), 0
            with _open(url, 'GET', headers, timeout=120) as r:
                if have and r.status != 206:              # the server ignored the Range header: restart from zero
                    have = 0
                elif have:                                # a partial answer must start where the part file ends, else the splice would be corrupt
                    cr = r.headers.get('Content-Range', '')
                    if not cr.startswith('bytes %d-' % have):
                        raise ConnectionError('unexpected Content-Range %r for a request starting at %d' % (cr, have))
                mode = 'ab' if have else 'wb'
                with open(part, mode) as f:
                    while True:
                        b = r.read1(1 << 20)                          # returns what has arrived; b'' at the end or when the connection drops
                        if not b:
                            break
                        f.write(b)
                        got += len(b)
                        if time.time() - t0 > 5:
                            log('    %.1f / %.1f MB' % ((have + got) / 1e6, total / 1e6))
                            t0 = time.time()
        except urllib.error.HTTPError as e:
            if e.code == 416:                             # nothing left to fetch (or a bad part file)
                if os.path.exists(part) and os.path.getsize(part) != total:
                    os.remove(part)
                continue
            log('    HTTP %s, retrying' % e.code)
            time.sleep(min(2 * (attempt + 1), 10))
        except (urllib.error.URLError, http.client.HTTPException, ConnectionError, TimeoutError, OSError) as e:
            log('    transfer interrupted (%r), resuming' % (e,))
            time.sleep(min(2 * (attempt + 1), 10))
    size = os.path.getsize(part) if os.path.exists(part) else -1
    if size != total:
        raise RuntimeError('incomplete download of %s: %d of %d bytes (re-run to resume)' % (url, size, total))
    if not finalize:
        return part
    os.replace(part, dest)
    return dest


def zip_members(path):
    """[(member name, size)] of a zip; RuntimeError when it is not a zip or a member fails its CRC check."""
    try:
        with zipfile.ZipFile(path) as z:
            bad = z.testzip()
            if bad is not None:
                raise RuntimeError('CRC error in %s member %s' % (path, bad))
            return [(i.filename, i.file_size) for i in z.infolist()]
    except zipfile.BadZipFile as e:
        raise RuntimeError('%s is not a valid zip file (%s)' % (path, e))


class ContentChanged(RuntimeError):
    """The file the server delivers now matches neither the local file nor the recorded digest."""


def ensure_file(url, dest, total, rec, accept_new=False, log=print):
    """Make ``dest`` the file the server delivers and return (sha256, action).

    rec: the recorded entry ({'sha256': ...}) or None. action: 'verified' (the local file has the right size and the recorded digest, nothing
    downloaded); 'downloaded' (fresh download, no record or the same digest as the record); 'recorded' (a local file without a record, taken as it is);
    'repaired' (the local file did not match its record, a fresh download did: the local file was damaged and has been replaced);
    'accepted-new' (the fresh download differs from the record and ``accept_new`` allowed recording it). A local file with a recorded digest that
    does not match is never trusted: it is compared with a fresh download made to '.part'. Raises ContentChanged when the fresh download differs
    from the record and ``accept_new`` is false (nothing is replaced then). --accept-new only ever accepts the digest of a FRESH download."""
    part = dest + '.part'
    have = os.path.exists(dest)
    if have and os.path.getsize(dest) != total:
        os.remove(dest)
        have = False
    local = sha256_of(dest) if have else None
    if have and rec is None:
        return local, 'recorded'                                    # a file without a record: recorded as it is (first use; delete it to force a download)
    if have and local == rec['sha256']:
        return local, 'verified'
    if not have:
        download(url, dest, total, log=log, finalize=False)         # resumes an interrupted '.part' of an earlier run
        fresh = sha256_of(part)
        action = 'downloaded'
    else:
        if os.path.exists(part):
            os.remove(part)                                         # a leftover could be damaged too: compare against a complete, fresh copy
        log('    local file does not match the recorded sha256; downloading a fresh copy to compare')
        download(url, dest, total, log=log, finalize=False)
        fresh = sha256_of(part)
        action = 'repaired'
    if rec is not None and fresh != rec['sha256']:
        if not accept_new:
            os.remove(part)
            raise ContentChanged('sha256 of a fresh download of %s (%s) differs from the recorded one (%s)%s: the server content changed. '
                                 'Re-run with --accept-new if that is intended.' % (os.path.basename(dest), fresh, rec['sha256'],
                                                                                    '; the local file (%s) is damaged as well' % local if have and local != fresh else ''))
        action = 'accepted-new'
    os.replace(part, dest)
    return fresh, action


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--sets', nargs='+', default=['MT0.N1.LD'], choices=SETS, help='datasets to fetch (default MT0.N1.LD)')
    ap.add_argument('--data-dir', default=os.environ.get('SMLM_DATA', DEFAULT_DATA), help='target directory (default $SMLM_DATA or %s)' % DEFAULT_DATA)
    ap.add_argument('--checksums', default=CHECKSUMS)
    ap.add_argument('--head-only', action='store_true', help='only check the URLs and their sizes')
    ap.add_argument('--gt-only', action='store_true', help='skip the image stacks')
    ap.add_argument('--accept-new', action='store_true', help='record the digest of a FRESH download that differs from the recorded one (never that of a file already on disk)')
    a = ap.parse_args(argv)

    data = os.path.abspath(a.data_dir)
    if os.path.commonpath([data, REPO]) == REPO:
        sys.exit('refusing to put raw data inside the repository: %s (set SMLM_DATA or --data-dir)' % data)
    os.makedirs(data, exist_ok=True)
    book = json.load(open(a.checksums)) if os.path.exists(a.checksums) else {}

    ok = True
    for name in a.sets:
        for kind, (rel, url) in files_for(name).items():
            if kind == 'stack' and a.gt_only:
                continue
            dest = os.path.join(data, rel)
            total, lastmod, etag = head(url)
            want = EXPECTED_BYTES[name][kind]
            flag = 'ok' if total == want else 'SIZE MISMATCH (expected %s)' % want
            print('%-52s %12d bytes  %s' % (rel, total, flag))
            if total != want:
                ok = False
                continue
            if a.head_only:
                continue
            rec = book.get(rel)
            try:
                digest, action = ensure_file(url, dest, total, rec, accept_new=a.accept_new)
            except ContentChanged as e:
                sys.exit(str(e))
            if action == 'verified':
                print('    present, sha256 matches the record')
                continue
            entry = {'url': url, 'bytes': total, 'sha256': digest, 'last_modified': lastmod, 'etag': etag}
            if kind == 'stack':
                try:
                    members = zip_members(dest)
                except RuntimeError as e:
                    os.remove(dest)                                  # nothing is recorded for a damaged archive; the next run downloads it again
                    sys.exit('%s: %s. The file was removed and no digest recorded; run again to download it anew.' % (rel, e))
                entry['zip_members'] = [m[0] for m in members]
                entry['zip_member_bytes'] = [m[1] for m in members]
            book[rel] = entry
            print('    %s, sha256 %s' % ({'downloaded': 'downloaded', 'recorded': 'present without a record, recorded as it is', 'repaired': 'repaired: the local file was damaged and has been replaced by a fresh download',
                                          'accepted-new': 'downloaded, new digest recorded (--accept-new)'}[action], digest))
            with open(a.checksums, 'w') as f:
                json.dump(book, f, indent=2, sort_keys=True)
                f.write('\n')
    if not ok:
        sys.exit('size mismatch: the remote file is not the one the sizes were recorded for')


if __name__ == '__main__':
    main()
