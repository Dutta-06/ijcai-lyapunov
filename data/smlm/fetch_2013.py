#!/usr/bin/env python
"""Resumable, checksummed download of the four 2D TRAINING datasets (public ground truth) of the 2013 localization-microscopy challenge.

Raw data stays out of the repository: this script fetches into $SMLM_DATA_2013 (default /home/user/smlm_data_2013) and refuses to write inside the
repository; only this script, the digests (``checksums_2013.json``) and derived results (``results/smlm/smlm2013*``) are committed.

Source (read 2026-11-07): https://srm.epfl.ch/srm/dataset/challenge-2D-simulation/index.html, "Realistic Simulations 2D, Datasets of the Challenge
2013 ... generated with a 2D theoretical PSF", section "Training Datasets (with Ground-truth Localizations)", and the four sub-pages
Tubulins_I, Tubulins_II, Bundled_Tubes_High_Density, Bundled_Tubes_Long_Sequence (<set>/index.html). Each sub-page lists four files, linked
relative to it as ../../../Data/challenge-2D-simulation/<set>/Data/<file>.zip:

    sequence.zip      the image sequence: one 16-bit big-endian TIFF per frame, sequence/00001.tif ... (checked by ``--manifest``)
    fluorophores.zip  THE GROUND TRUTH: fluorophores/activation.csv (one row per activation: id, frame, xnano, ynano, znano, intensity; the same
                      header as the 2016 ``activations.csv``) and fluorophores/frames/NNNNN.csv (the same rows split by frame). NOT in the
                      page's text: the page calls the file "fluorophores" and shows "Views of the sample"; the content was read from the file.
    oracle.zip        NOT the localizations: four projections of the structure (oracle/projection-{avg,max}-{camera,high}-resolution.{tif,jpg}),
                      optional (``--extras``)
    psf.zip           the simulated PSF (a TIFF, psf/psf-10.0nm.tif and an illustration), optional (``--extras``)

Every file is checked against the server's Content-Length (HEAD) and against the size seen when the URLs were verified; the sha256 of every
file is recorded in ``checksums_2013.json`` next to this script after the first verified download and re-verified by every later run. The
integrity discipline is that of ``data/smlm/fetch.py`` (whose helpers are imported): a local file that no longer matches its record is
NOT taken to mean that the server changed; a fresh copy is downloaded to a '.part' file and compared, and only when that differs from the
record as well does the script stop with "server content changed" (``--accept-new`` then records the digest of that FRESH download, never that of
a file already on disk). A zip that fails its CRC check is removed and nothing is recorded. Downloads resume from a '.part' file. Nothing is ever
extracted or executed: downstream code (``experiments/smlm2013.py``) reads the zip members in memory.

    python data/smlm/fetch_2013.py                       # sequence + ground truth (fluorophores.zip) of all four sets (about 465 MB)
    python data/smlm/fetch_2013.py --sets Tubulins_I     # one set
    python data/smlm/fetch_2013.py --extras              # also oracle.zip and psf.zip
    python data/smlm/fetch_2013.py --head-only           # only check the URLs and sizes
    python data/smlm/fetch_2013.py --manifest            # additionally print the facts of the pages and what the local files contain

Only the Python standard library is needed to download (TLS verification is always on; proxy and CA bundle come from the environment);
``--manifest`` reads a TIFF header with Pillow and the member names / header of the ground truth with the standard library.
"""
import argparse
import io
import json
import os
import sys
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, HERE)
import fetch as F                                    # noqa: E402  (the helpers of the 2016 sets: head, ensure_file, zip_members, sha256_of)

DEFAULT_DATA = '/home/user/smlm_data_2013'
CHECKSUMS = os.path.join(HERE, 'checksums_2013.json')
ROOT = 'https://srm.epfl.ch/srm'
PAGE_ROOT = ROOT + '/dataset/challenge-2D-simulation'
SETS = ('Tubulins_I', 'Tubulins_II', 'Bundled_Tubes_High_Density', 'Bundled_Tubes_Long_Sequence')
KINDS = ('sequence', 'fluorophores', 'oracle', 'psf')
DEFAULT_KINDS = ('sequence', 'fluorophores')

# Sizes in bytes as returned by HEAD when the URLs were verified (2026-11-07); a different size is refused.
EXPECTED_BYTES = {
    'Tubulins_I': {'sequence': 190604305, 'oracle': 38077701, 'psf': 1190682, 'fluorophores': 5513600},
    'Tubulins_II': {'sequence': 202735580, 'oracle': 40777930, 'psf': 1190682, 'fluorophores': 5513600},
    'Bundled_Tubes_High_Density': {'sequence': 2077674, 'oracle': 4217912, 'psf': 4307112, 'fluorophores': 3817424},
    'Bundled_Tubes_Long_Sequence': {'sequence': 55752091, 'oracle': 4218346, 'psf': 4307112, 'fluorophores': 6923245},
}

# ----------------------------------------------------------------------------------------------------------------------------------
# What the pages say (read 2026-11-07; the numbers are copied from the pages' parameter tables, the quotes are verbatim)
# ----------------------------------------------------------------------------------------------------------------------------------
PAGE_FACTS = {
    'Tubulins_I': dict(
        title='Tubulins_I, 2D. "Low level of read-out noise autofluorescence background."',
        sample='7 tubulins of constant diameter 25 nm and 1 of 40 nm; depth of the sample 0 to 300 nm; 100000 fluorophores activated over 2400 frames '
               '(the parameter tables say "Frames From 1 to 2401")',
        frame='256 x 256 px, pixel 150 nm, field of view 38400 nm', optics='wavelength 723 nm, NA 1.40, diffraction limit 258.21 nm',
        psf='XY Gaussian, FWHM 258.21 nm (sigma 109.65 nm = 0.731 px); Z function exponential: 1 x FWHM at z = 0, 2 x FWHM at z = 500 nm',
        camera='QE 1.00 e-/photon; ADC gain 1.00 DN/e-, offset 0, baseline 100 DN, saturation 16383 DN, 14-bit, file TIFF 16-bit',
        noise='read-out Gaussian sd 10.0; dark noise Poisson 1.0; shot noise Poisson; "EM gain: multiplicative 1.41" (excess noise factor)',
        background='autofluorescence from 1000 sources (gain 50, size 200 nm, diffusion 1000 nm, displacement 500 nm, 5 scales, "Evolution in the sequence: '
                   'Dynamic", 20 % change), added at 75 nm pixels, constant background level mean 0 sd 0',
        emitters='41.64 activated fluorophores per frame on average (100000 over 2400 frames; activation rate 1.00 mol./frame; alive 40 ms of a 20 ms frame; '
                 'photons per molecule and frame: flux 1944.41 photons/mol. at the centre of a non-uniform "bumped" excitation, FWHM 30000 nm)',
        snr='peak signal 151..2141 (mean 502), background sd 10.35..281.02 (mean 26.67)'),
    'Tubulins_II': dict(
        title='Tubulins_II, 2D. "High level of read-out noise and autofluorescence background."',
        sample='same structure as Tubulins_I (same tube table, same 100000 fluorophores over 2400 frames)',
        frame='256 x 256 px, pixel 150 nm, field of view 38400 nm', optics='as Tubulins_I', psf='as Tubulins_I',
        camera='as Tubulins_I', noise='as Tubulins_I except read-out Gaussian sd 30.0',
        background='as Tubulins_I (1000 sources, gain 50)', emitters='as Tubulins_I', snr='not given on the page'),
    'Bundled_Tubes_High_Density': dict(
        title='Bundled_Tubes_High_Density, 2D. "Bundles of 8 tubes of 30 nm diameter. High-density, 81049 fluorophores on 168 frames"',
        sample='8 tubes of radius 15 nm (membrane 3 nm), total length 45.5 um, z between 30 and 68.5 nm; THE PAGE CONTRADICTS ITSELF on the length: '
               '"168 frames" in the headline, "Frames From 1 to 361" in the parameter tables (the stack decides: see --manifest)',
        frame='64 x 64 px, pixel 100 nm, field of view 6400 nm', optics='wavelength 723 nm, NA 1.40, diffraction limit 258.21 nm',
        psf='XY Gaussian, FWHM 258.21 nm (sigma 109.65 nm = 1.097 px); Z function exponential: 1 x FWHM at z = 0, 2 x FWHM at z = 300 nm',
        camera='QE 1.00 e-/photon; ADC gain 1.00 DN/e-, offset 0, baseline 100 DN, saturation 16383 DN, 14-bit, file TIFF 16-bit',
        noise='read-out Gaussian sd 5.0; dark noise Poisson 1.0; shot noise Poisson; "EM gain: multiplicative 1.41"',
        background='autofluorescence from 2000 sources (gain 30, size 200 nm, diffusion 1000 nm, displacement 500 nm, 5 scales, Dynamic, 20 % change), added at 100 nm pixels',
        emitters='81049 activated fluorophores', snr='peak signal 216..2051 (mean 868), background sd 20.50..450.39 (mean 184.87)'),
    'Bundled_Tubes_Long_Sequence': dict(
        title='Bundled_Tubes_Long_Sequence, 2D. "Bundles of 8 tubes of 30 nm diameter. Sparse density, 81049 molecules on 12000 frames"',
        sample='same bundles as the high-density set (page: "Frames From 1 to 12000")',
        frame='64 x 64 px, pixel 100 nm, field of view 6400 nm', optics='as Bundled_Tubes_High_Density', psf='as Bundled_Tubes_High_Density',
        camera='as Bundled_Tubes_High_Density', noise='as Bundled_Tubes_High_Density (read-out sd 5.0, EM multiplicative 1.41)',
        background='as Bundled_Tubes_High_Density (2000 sources, gain 30)', emitters='81049 activated molecules',
        snr='peak signal 160..1076 (mean 361), background sd 8.77..140.11 (mean 22.22)'),
}
PAGE_COMMON = ('"They are generated with a 2D theoretical PSF." Main page lists the reference as "Sage et al., Quantitative Evaluation of Software Packages for Single-Molecule '
               'Localization Microscopy ... Nature Methods 2019" with the link of the 2019 paper (s41592-019-0364-4), the sub-pages as "Sage et al., Nature Methods 2015" with the link '
               'nmeth.3442. Published 02 February 2013.')
NOT_ON_THE_PAGES = ('which of the four files holds the ground-truth localizations (fluorophores.zip: the pages call it "fluorophores" and show "Views of the sample"; oracle.zip holds projection images)',
                    'the format of the ground truth (columns, units, frame numbering, photons): read from the file',
                    'the frame numbering of the TIFF files (the pages say "From 1 to N", the files decide) and the true number of frames of two sets',
                    'what unit the Gaussian noise sds (5, 10, 30) and the "gain of sources" are in; whether the photon numbers are per frame or per second')


def files_for(name, kinds=DEFAULT_KINDS):
    """kind -> (relative local path, url) for one dataset."""
    return {k: ('%s/%s.zip' % (name, k), '%s/Data/challenge-2D-simulation/%s/Data/%s.zip' % (ROOT, name, k)) for k in kinds}


def tiff_members(members):
    """The (frame number, member name) pairs of the per-frame TIFFs of a zip listing [(name, size)], sorted by the integer in the file name."""
    out = []
    for name, size in members:
        stem, ext = os.path.splitext(os.path.basename(name))
        if size > 0 and ext.lower() in ('.tif', '.tiff') and stem.isdigit():
            out.append((int(stem), name))
    return sorted(out)


def sequence_facts(path):
    """(n_frames, first number, last number, consecutive?, rows, columns, Pillow mode) of the per-frame TIFFs in a sequence zip (the header of the first one)."""
    with zipfile.ZipFile(path) as z:
        t = tiff_members([(i.filename, i.file_size) for i in z.infolist()])
        if not t:
            return None
        nums = [k for k, _ in t]
        rows = cols = mode = None
        try:
            from PIL import Image
            im = Image.open(io.BytesIO(z.read(t[0][1])))
            cols, rows, mode = im.size[0], im.size[1], im.mode
        except ImportError:
            pass
        return (len(t), nums[0], nums[-1], nums == list(range(nums[0], nums[0] + len(nums))), rows, cols, mode)


def oracle_facts(path, nlines=3):
    """[(member, size, first lines)] of the members of an oracle zip (the first lines of text members, decoded as ASCII/UTF-8 with replacement)."""
    out = []
    with zipfile.ZipFile(path) as z:
        for i in z.infolist():
            if i.is_dir():
                continue
            with z.open(i) as f:
                head = f.read(2000).decode('utf-8', errors='replace').splitlines()[:nlines]
            out.append((i.filename, i.file_size, head))
    return out


def gt_facts(path):
    """Facts of a fluorophores.zip read from the file: (activation member, header, number of rows, first and last frame number, number of per-frame csv files, rows of the per-frame files,
    whether every frame number of the activation file occurs once per row, i.e. every activation lasts one frame, is not tested here: rows only)."""
    with zipfile.ZipFile(path) as z:
        acts = [i.filename for i in z.infolist() if not i.is_dir() and os.path.basename(i.filename) == 'activation.csv']
        if len(acts) != 1:
            return None
        lines = z.read(acts[0]).decode('utf-8', errors='replace').splitlines()
        frames = [int(l.split(',')[1]) for l in lines[1:] if l.strip()]
        per = [i.filename for i in z.infolist() if not i.is_dir() and '/frames/' in '/' + i.filename and i.filename.endswith('.csv')]
        return (acts[0], lines[0], len(lines) - 1, min(frames), max(frames), len(per))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--sets', nargs='+', default=list(SETS), choices=SETS, help='datasets to fetch (default: all four)')
    ap.add_argument('--extras', action='store_true', help='also fetch oracle.zip (projections of the structure) and psf.zip')
    ap.add_argument('--data-dir', default=os.environ.get('SMLM_DATA_2013', DEFAULT_DATA),
                    help='target directory (default $SMLM_DATA_2013 or %s)' % DEFAULT_DATA)
    ap.add_argument('--checksums', default=CHECKSUMS)
    ap.add_argument('--head-only', action='store_true', help='only check the URLs and their sizes')
    ap.add_argument('--accept-new', action='store_true',
                    help='record the digest of a FRESH download that differs from the recorded one (never that of a file already on disk)')
    ap.add_argument('--manifest', action='store_true', help='print the facts read from the pages and the contents of the local files')
    a = ap.parse_args(argv)

    data = os.path.abspath(a.data_dir)
    if os.path.commonpath([data, REPO]) == REPO:
        sys.exit('refusing to put raw data inside the repository: %s (set SMLM_DATA_2013 or --data-dir)' % data)
    os.makedirs(data, exist_ok=True)
    book = json.load(open(a.checksums)) if os.path.exists(a.checksums) else {}
    kinds = KINDS if a.extras else DEFAULT_KINDS

    if a.manifest:
        print('Facts read from the pages (2026-11-07; %s/index.html and <set>/index.html):' % PAGE_ROOT)
        print('  %s' % PAGE_COMMON)
        for s in a.sets:
            print('  == %s' % s)
            for k, v in PAGE_FACTS[s].items():
                print('     %-10s %s' % (k, v))
        print('NOT on the pages: %s.' % '; '.join(NOT_ON_THE_PAGES))

    ok = True
    for name in a.sets:
        for kind, (rel, url) in files_for(name, kinds).items():
            dest = os.path.join(data, rel)
            total, lastmod, etag = F.head(url)
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
                digest, action = F.ensure_file(url, dest, total, rec, accept_new=a.accept_new)
            except F.ContentChanged as e:
                sys.exit(str(e))
            if action == 'verified':
                print('    present, sha256 matches the record')
            else:
                try:
                    members = F.zip_members(dest)
                except RuntimeError as e:
                    os.remove(dest)                                  # nothing is recorded for a damaged archive; the next run downloads it again
                    sys.exit('%s: %s. The file was removed and no digest recorded; run again to download it anew.' % (rel, e))
                book[rel] = {'url': url, 'bytes': total, 'sha256': digest, 'last_modified': lastmod, 'etag': etag,
                             'zip_members': [m[0] for m in members[:3]] + (['...'] if len(members) > 3 else []), 'zip_n_members': len(members),
                             'zip_member_bytes_total': int(sum(m[1] for m in members))}
                print('    %s, sha256 %s' % ({'downloaded': 'downloaded', 'recorded': 'present without a record, recorded as it is',
                                              'repaired': 'repaired: the local file was damaged and has been replaced by a fresh download',
                                              'accepted-new': 'downloaded, new digest recorded (--accept-new)'}[action], digest))
                with open(a.checksums, 'w') as f:
                    json.dump(book, f, indent=2, sort_keys=True)
                    f.write('\n')
            if a.manifest:
                if kind == 'sequence':
                    sf = sequence_facts(dest)
                    print('    sequence (from the zip itself): ' + ('no per-frame TIFF' if sf is None else
                          '%d per-frame TIFFs numbered %d..%d (%s), first frame %s rows x %s columns, Pillow mode %s' % (
                              sf[0], sf[1], sf[2], 'consecutive' if sf[3] else 'WITH A GAP', sf[4], sf[5], sf[6])))
                elif kind == 'fluorophores':
                    gf = gt_facts(dest)
                    print('    ground truth (from the zip itself): ' + ('no activation.csv' if gf is None else '%s, header "%s", %d rows, frame numbers %d..%d, %d per-frame csv files' % gf))
                elif kind == 'oracle':
                    for m, size, head in oracle_facts(dest, 0):
                        print('    oracle member %s (%d bytes)' % (m, size))
    if not ok:
        sys.exit('size mismatch: the remote file is not the one the sizes were recorded for')


if __name__ == '__main__':
    main()
