#!/usr/bin/env python
"""Resumable, checksummed download of the REAL 2D single-molecule localization sequences of the 2013 challenge (Tubulin Alexa647).

Raw data stays out of the repository: this script fetches into $SMLM_DATA_REAL (default /home/user/smlm_data_real) and refuses to write
inside the repository; only this script, the digests below and derived results (``results/smlm/real_*``) are committed.

What the pages say (https://srm.epfl.ch/srm/dataset/challenge-2D-real/index.html and its sub-pages Real_High_Density/index.html,
Real_Long_Sequence/index.html, Tubulin_AF647/index.html, all read on 2026-11-07; quoted or paraphrased closely)
----------------------------------------------------------------------------------------------------------------------------------
* Sample: "a fixed cell, stained with mouse anti-alpha-tubulin primary antibody and Alexa647 secondary antibody" (Tubulin AF647).
* "Datasets of the Challenge 2013", published 02 February 2013, courtesy of Nicolas Olivier, Debora Keller and Suliana Manley; reference
  Sage et al., Nature Methods 2015 (nmeth.3442). Each sub-page says "Real dataset of the challenge 2013. Ground-truth localization unknown."
* Real_High_Density: "High density settings, 500 frames"; sequence.zip, 11 Mb on the page (10,878,456 bytes by HEAD).
* Real_Long_Sequence: "Long sequence settings, 15000 frames" on its own page; sequence.zip, 71 Mb on the page (71,213,860 bytes by HEAD).
* Main page, entry Tubulin AF647 (the 218 Mb sequence.zip of the Tubulin_AF647 sub-page, NOT fetched here): "The dataset consists in 9990
  frames of 128x128 pixels stored as TIFF in 16 bits." The entries Real Long Sequence ("Long sequence settings, 15000 frames.") and Real High
  Density ("High density settings, 500 frames.") carry no frame size; a fourth entry, cell4 (not fetched), is "Basic exemple, 2000 frames, 256x256".
* Download links are relative to the sub-pages: ../../../Data/challenge-2D-real/<set>/Data/sequence.zip (-> URLS below). The results
  directories of other software packages that the pages list are not fetched (and the Data/ directory listing answers 403).

What NO page says (checked: the three sub-pages, the main page, the hub's home page; bigwww.epfl.ch redirects to the hub's home page)
----------------------------------------------------------------------------------------------------------------------------------
Pixel size, camera model, EM gain, electrons per count (ADU), baseline, exposure time, objective / NA, emission wavelength, laser power,
buffer. Anything built on one of them is an ASSUMPTION of the experiment that uses these data (``experiments/exp_smlm_real.py`` states
which); in particular the pixel size of 100 nm used to express positions in nm is not on the pages.

What the files themselves say (``--manifest`` reads it; not on the pages for the Long Sequence): the zips hold one 16-bit big-endian TIFF per
frame, sequence/00001.tif, sequence/00002.tif, ... Real_High_Density: 500 frames of 128 x 128 px (the frame size is not stated on its page).
Real_Long_Sequence: 15000 frames (agrees with its page) of 64 x 64 px, NOT the 128x128 of the main page's Tubulin AF647 entry.

Integrity (same discipline as data/smlm/fetch.py, whose helpers are reused)
----------------------------------------------------------------------------------------------------------------------------------
Every file is checked against the server's Content-Length (HEAD) and against the size recorded in EXPECTED_BYTES, and its sha256
against CHECKSUMS (recorded below after the first verified download, 2026-11-07). A local file that does not match its record is not taken
to mean that the server changed: a fresh copy is downloaded to a '.part' file and compared; only if that differs from the record as well
does the script stop with "server content changed" (``--accept-new`` then records the digest of that FRESH download). A zip that fails its
CRC check is removed. Nothing is ever extracted or executed: downstream code reads the zip member in memory.

    python data/smlm/fetch_real.py                       # Real_High_Density (500 frames)
    python data/smlm/fetch_real.py --sets Real_High_Density Real_Long_Sequence
    python data/smlm/fetch_real.py --head-only           # only check the URLs and sizes
    python data/smlm/fetch_real.py --manifest            # additionally print the facts above and what the local stacks contain

Only the Python standard library is needed to download (TLS verification always on; proxy and CA bundle from the environment);
``--manifest`` reads the TIFF header with Pillow when it is importable.
"""
import argparse
import io
import os
import sys
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, HERE)
import fetch as F                                    # noqa: E402  (the helpers of the simulated sets: head, ensure_file, zip_members, sha256_of)

DEFAULT_DATA = '/home/user/smlm_data_real'
ROOT = 'https://srm.epfl.ch/srm'
SETS = ('Real_High_Density', 'Real_Long_Sequence')
PAGES = {name: '%s/dataset/challenge-2D-real/%s/index.html' % (ROOT, name) for name in SETS}
URLS = {name: '%s/Data/challenge-2D-real/%s/Data/sequence.zip' % (ROOT, name) for name in SETS}

# Sizes in bytes as returned by HEAD when the URLs were verified (2026-11-07).
EXPECTED_BYTES = {'Real_High_Density': 10878456, 'Real_Long_Sequence': 71213860}

# sha256 of the downloaded files, recorded after the first verified download (2026-11-07; zip CRC checked, size equal to HEAD).
CHECKSUMS = {
    'Real_High_Density': '2138d15662900ad21be60d0c249851d50d94f8b90bddb1eabdca1e90d0780fcb',
    'Real_Long_Sequence': '7c0e1b597edee79b5f3aa68740bcef9c9822cef5928d06fe4558622e6a650fa5',
}

# The facts printed by --manifest (the page text, verbatim or closely paraphrased) and the facts the pages do NOT give.
PAGE_FACTS = [
    ('sample', 'fixed cell, mouse anti-alpha-tubulin primary antibody, Alexa647 secondary antibody (Tubulin AF647)'),
    ('provenance', 'challenge 2013, published 02 February 2013; courtesy of Nicolas Olivier, Debora Keller and Suliana Manley; ground truth unknown'),
    ('Real_High_Density', 'high density settings, 500 frames; sequence.zip 11 Mb'),
    ('Real_Long_Sequence', 'long sequence settings, 15000 frames; sequence.zip 71 Mb (no frame size on its page)'),
    ('Tubulin_AF647', '(not fetched) main page: "9990 frames of 128x128 pixels stored as TIFF in 16 bits"; sequence.zip 218 Mb'),
]
NOT_ON_THE_PAGES = ('pixel size', 'camera model', 'EM gain', 'electrons per ADU count', 'baseline / offset', 'exposure time', 'objective / NA',
                    'emission wavelength')


def local_path(name, data):
    return os.path.join(data, name, 'sequence.zip')


def tiff_members(members):
    """The (frame number, member name) pairs of the per-frame TIFFs of a zip listing [(name, size)], sorted by the integer in the file name
    (the sequences hold one file per frame, sequence/00001.tif, sequence/00002.tif, ...)."""
    out = []
    for name, size in members:
        base = os.path.basename(name)
        stem, ext = os.path.splitext(base)
        if size > 0 and ext.lower() in ('.tif', '.tiff') and stem.isdigit():
            out.append((int(stem), name))
    return sorted(out)


def summarise(members):
    """One line about the zip members: how many per-frame TIFFs, whether their numbers are 1..n without a gap, their sizes."""
    t = tiff_members(members)
    sizes = sorted({s for n, s in members if s > 0})
    if not t:
        return '%d members, no per-frame TIFF' % len(members)
    nums = [k for k, _ in t]
    gap = nums != list(range(nums[0], nums[0] + len(nums)))
    return '%d per-frame TIFFs %s .. %s (%s), member sizes in bytes %s' % (
        len(t), t[0][1], t[-1][1], 'numbers have a gap' if gap else 'consecutive numbers', ', '.join(map(str, sizes[:4])) + (' ...' if len(sizes) > 4 else ''))


def stack_facts(path):
    """(n_frames, height, width, mode, dtype) of the sequence in the zip at ``path``: the number of per-frame TIFF members and the header of the
    first one (Pillow, read in memory; nothing is extracted), or None when Pillow is not importable."""
    try:
        from PIL import Image
    except ImportError:
        return None
    import numpy as np
    with zipfile.ZipFile(path) as z:
        t = tiff_members([(i.filename, i.file_size) for i in z.infolist()])
        if not t:
            return (0, None, None, None, None)
        im = Image.open(io.BytesIO(z.read(t[0][1])))
        return (len(t), im.size[1], im.size[0], im.mode, str(np.asarray(im).dtype))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--sets', nargs='+', default=['Real_High_Density'], choices=SETS, help='datasets to fetch (default Real_High_Density)')
    ap.add_argument('--data-dir', default=os.environ.get('SMLM_DATA_REAL', DEFAULT_DATA),
                    help='target directory (default $SMLM_DATA_REAL or %s)' % DEFAULT_DATA)
    ap.add_argument('--head-only', action='store_true', help='only check the URLs and their sizes')
    ap.add_argument('--accept-new', action='store_true',
                    help='record the digest of a FRESH download that differs from the recorded one (never that of a file already on disk)')
    ap.add_argument('--manifest', action='store_true', help='print the facts read from the pages and the contents of the local stacks')
    a = ap.parse_args(argv)

    data = os.path.abspath(a.data_dir)
    if os.path.commonpath([data, REPO]) == REPO:
        sys.exit('refusing to put raw data inside the repository: %s (set SMLM_DATA_REAL or --data-dir)' % data)
    os.makedirs(data, exist_ok=True)

    if a.manifest:
        print('Facts read from the pages (2026-11-07; %s):' % ', '.join(PAGES[s] for s in SETS))
        for k, v in PAGE_FACTS:
            print('  %-18s %s' % (k, v))
        print('NOT given on any page: %s.' % ', '.join(NOT_ON_THE_PAGES))

    ok = True
    for name in a.sets:
        url, dest = URLS[name], local_path(name, data)
        total, lastmod, etag = F.head(url)
        want = EXPECTED_BYTES[name]
        flag = 'ok' if total == want else 'SIZE MISMATCH (expected %s)' % want
        print('%-34s %12d bytes  %s' % (os.path.relpath(dest, data), total, flag))
        if total != want:
            ok = False
            continue
        if a.head_only:
            continue
        sha = CHECKSUMS.get(name)
        rec = {'sha256': sha} if sha else None
        try:
            digest, action = F.ensure_file(url, dest, total, rec, accept_new=a.accept_new)
        except F.ContentChanged as e:
            sys.exit(str(e))
        try:
            members = F.zip_members(dest)
        except RuntimeError as e:
            os.remove(dest)
            sys.exit('%s: %s. The file was removed; run again to download it anew.' % (name, e))
        print('    %s; sha256 %s; zip: %s' % ({'verified': 'present, sha256 matches the record', 'downloaded': 'downloaded',
                                           'recorded': 'present (no record yet)', 'repaired': 'repaired (the local file was damaged and was replaced)',
                                           'accepted-new': 'downloaded, new digest accepted (--accept-new)'}[action], digest, summarise(members)))
        if sha is None:
            print('    NO RECORDED DIGEST for %s: after checking the line above, put  %r  into CHECKSUMS in this script.' % (name, digest))
        elif digest != sha:
            print('    NEW DIGEST (--accept-new): to keep it, put %r into CHECKSUMS in this script.' % digest)
        if a.manifest:
            sf = stack_facts(dest)
            if sf is None:
                print('    stack: Pillow is not importable, header not read')
            else:
                print('    stack (from the first TIFF itself): %s frames of %s x %s px (rows x columns), Pillow mode %s, numpy dtype %s' % sf)
    if not ok:
        sys.exit('size mismatch: the remote file is not the one the sizes were recorded for')


if __name__ == '__main__':
    main()
