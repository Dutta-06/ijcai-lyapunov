"""
The four 2D TRAINING datasets of the 2013 localization-microscopy challenge as benchmark sets for the harness of `smlm_common.py`
==================================================================================================================================

The EPFL Super-Resolution Microscopy Hub hosts the simulated 2D datasets of the 2013 challenge (https://srm.epfl.ch/srm/dataset/challenge-2D-simulation/index.html,
"generated with a 2D theoretical PSF"), with public ground truth for four "training" sets: a realistic structure of 8 tubulins (``Tubulins_I`` low,
``Tubulins_II`` high read-out noise; 256 x 256 px of 150 nm, 2401 frames, 100,000 fluorophores) and bundles of 8 tubes (``Bundled_Tubes_High_Density``, 361 frames, and
``Bundled_Tubes_Long_Sequence``, 12,000 frames, the same 81,049 fluorophores; 64 x 64 px of 100 nm). This module exposes each set the way `smlm_common.py` exposes the
2016 sets: frames in photons (``load_frames``, ``to_photons``), the ground truth as the harness's structured array (``load_gt``: 0-BASED frame index, positions in nm in the
harness's CORNER convention), scoring with the field of the set (``score``), and the declared development / test frame ranges (``PROTOCOL``). Files come from
``data/smlm/fetch_2013.py`` (into $SMLM_DATA_2013, default /home/user/smlm_data_2013; never into the repository).

What the pages say, what the files say, what is assumed (read 2026-11-07)
-------------------------------------------------------------------------
PAGES (``data/smlm/fetch_2013.py --manifest`` prints the full table): frame size and pixel (256 px / 150 nm for the tubulins, 64 px / 100 nm for the bundles); wavelength 723 nm,
NA 1.40, diffraction limit 258.21 nm; PSF "XY function: Gaussian, FWHM 258.21 nm; Z function: Exponential, focal plane 1 x FWHM at 0 nm, defocus plane 2 x FWHM at 500 nm
(tubulins) / 300 nm (bundles)"; camera "photon converter factor QE 1.00 e-/photon", "electron conversion gain 1.00 DN/e-, offset 0, baseline 100 DN, 14-bit, saturation
16383", noise "read-out Gaussian sd 10.0 (Tubulins_I) / 30.0 (Tubulins_II) / 5.0 (bundles), dark noise Poisson 1.0, shot noise Poisson, EM gain multiplicative 1.41"; a structured,
"dynamic" autofluorescence background; the activation counts (41.64 fluorophores per frame on average for the tubulins). The pages contradict themselves on the number of frames
of the tubulin sets ("2400 frames" in the text, "From 1 to 2401" in the tables) and of the high-density bundle set (headline "168 frames", table "From 1 to 361"): THE FILES DECIDE.
FILES (read from them, ``fetch_2013.py --manifest``): ``sequence.zip`` holds one 16-bit big-endian TIFF per frame (``sequence/00001.tif`` ... in the bundle sets, ``00001.tif`` ... in the
tubulin sets): 2401 frames of 256 x 256 (Tubulins_I, _II), 361 of 64 x 64 (Bundled_Tubes_High_Density), 12,000 of 64 x 64 (Bundled_Tubes_Long_Sequence). THE GROUND TRUTH IS IN
``fluorophores.zip``, NOT in ``oracle.zip`` (which holds four projection images of the structure): ``activation.csv`` (100,000 / 100,000 / 81,049 / 81,049 rows; header
``Ground-truth,frame,xnano,ynano,znano,intensity``, the header of the 2016 ``activations.csv``) and ``frames/NNNNN.csv`` (the same rows split by frame; checked to equal the
activation file). Every row is ONE activation in ONE frame (each fluorophore is on in exactly one frame: all run lengths are 1, unlike the 2016 sets, whose emitters stay on for 2.3
frames on average), so consecutive frames share no emitter; Tubulins_I and _II have the same ground truth (same activations, same positions) and differ in the noise only; the two
bundle sets have the same 81,049 fluorophore positions and differ in the frames and photons.
ASSUMED / DERIVED: photons = DN - 100 (QE 1.00 e-/photon and 1.00 DN/e- from the page; checked against the data only through the amplitude ratios in ``results/smlm/smlm2013.md``);
the Gaussian noise sds of the pages do NOT give the observed constant variance (``results/smlm/smlm2013.md`` measures it), so the experiment calibrates it; the PSF width is measured
on the development frames (the page's FWHM gives sigma = 109.65 nm at the focal plane; the measured widths are in the results).

Conventions (CHECKED against the images: ``check_conventions`` and the tests, with a Picasso run matched to the ground truth)
-----------------------------------------------------------------------------------------------------------------------------
* Positions are in nm, x along the image COLUMNS, y along the ROWS, (0, 0) the OUTER CORNER of the top-left pixel: the centre of the pixel with 0-based indices (row i, column j) is at
  ((j + 0.5) p, (i + 0.5) p) nm, p the pixel size of the set; the field is [0, size p] nm in both axes. THE FILES FOLLOW THIS (the harness's corner convention) WITHOUT ANY SHIFT: the
  mean signed error (detection - truth) of Picasso's symmetric Gaussian fit is within about 2 nm of zero in both axes on all four sets (the bias of 20 development frames: Tubulins_I (1.7, 0.6) nm, Tubulins_II (1.6, -1.0), dense bundles (0.5, -1.0), long bundles (-1.1, 1.9)), against -p/2 (50 or 75 nm) under the other
  convention; swapping x and y or flipping an axis destroys the matching.
* Frame indices are 0-BASED here: frame k of a set is the TIFF numbered k + 1 and the ground-truth rows with file frame k + 1 (the page says "From 1 to N", the files are numbered
  1..N). Checked: the matched count of a Picasso run collapses for any shift of the ground-truth frame by +-1 or +-2 (sparse sets) and the excess intensity at the ground-truth positions
  peaks at offset 0 (also for the dense set).
* Photons: ``to_photons`` = (DN - 100) / 1.0 (the pages' conversion factors).
* The field size enters the harness only through the border exclusion of the "challenge setting" (``smlm_common.FIELD_NM`` = 6400 nm for the 2016 sets): ``score`` here sets it to the
  field of the set (38,400 nm for the tubulins) for the duration of the call.

Protocol (DECLARED before any test frame was scored; nothing is tuned or chosen on test frames)
-----------------------------------------------------------------------------------------------
``PROTOCOL[set]`` = dev (first, last + 1) and test (first, last + 1) in 0-based frame indices. Development frames: tuning and calibration. The long sequences are not used in full: the
ISS costs seconds per frame (see ``experiments/exp_smlm2013.py``), so the development and test ranges are subsets, chosen so that each holds more ground-truth emitters than the
2016 sets' (development 1253, test about 600 for MT0.N1.LD) and, for the test frames, at least 10 blocks of 10 frames for the block bootstrap:
    Tubulins_I / _II            dev frames 0-39 (40 frames, about 1,670 emitters), test frames 200-299 (100 frames, about 4,170 emitters); frames 40-199 are not used
                                (frames 100-105 of Tubulins_I were run once through the ISS on one 64 x 64 tile as a timing benchmark before this protocol was written; they are
                                outside both ranges, as is everything below 200 except 0-39; see the docstring of ``exp_smlm2013.py``)
    Bundled_Tubes_Long_Sequence dev frames 0-299 (300 frames, about 2,000 emitters), test frames 400-899 (500 frames, about 3,400 emitters) of 12,000
    Bundled_Tubes_High_Density  dev frames 0-29 (30 frames, about 6,700 emitters), test frames 100-159 (60 frames, about 13,500 emitters) of 361

Dropping in another method (data from ``fetch_2013.py``, in $SMLM_DATA_2013)::

    import smlm2013 as S, smlm_common as C
    a, b = S.PROTOCOL['Tubulins_I']['test']
    frames = S.load_frames(S.data_path('Tubulins_I', 'stack'), a, b)        # (b - a, 256, 256) DN; S.to_photons(frames) for photons
    gt = S.load_gt(S.data_path('Tubulins_I', 'gt'))                          # harness structured array, 0-based frames, nm, corner convention
    dets = my_method(frames)                                                  # list of per-frame (n, 2+) arrays: x_nm, y_nm (corner convention)
    res = S.score(dets, gt, np.arange(a, b), 'Tubulins_I')                    # the metrics of smlm_common.score with the field of the set
    ci = C.bootstrap_ci(res)                                                  # 95 % block-bootstrap intervals (blocks of 10 frames)
"""
import contextlib
import io
import os
import re
import sys
import zipfile

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
import smlm_common as C                                  # noqa: E402

DEFAULT_DATA = '/home/user/smlm_data_2013'
SETS = ('Tubulins_I', 'Tubulins_II', 'Bundled_Tubes_High_Density', 'Bundled_Tubes_Long_Sequence')
BASELINE_DN = 100.0                  # page: "Baseline 100.00 DN"
DN_PER_PHOTON = 1.0                  # page: QE 1.00 e-/photon x 1.00 DN/e-
GT_FRAME_BASE = 1                    # files number frames from 1
GT_MEMBER = 'activation.csv'         # in fluorophores.zip (at the top level in the tubulin sets, under fluorophores/ in the bundle sets)
EM_FACTOR = 1.41                     # page: "EM Gain ... Distribution: Multiplicative 1.41" (the slope of the noise variance against the signal, checked on the data by the experiment)

# facts of the pages (read 2026-11-07) and of the files; n_frames from the files (the pages disagree with themselves for two sets)
SPEC = {
    'Tubulins_I': dict(size=256, pixel_nm=150.0, n_frames=2401, read_sd=10.0, z_double_nm=500.0, structure='tubulins', n_gt=100000),
    'Tubulins_II': dict(size=256, pixel_nm=150.0, n_frames=2401, read_sd=30.0, z_double_nm=500.0, structure='tubulins', n_gt=100000),
    'Bundled_Tubes_High_Density': dict(size=64, pixel_nm=100.0, n_frames=361, read_sd=5.0, z_double_nm=300.0, structure='bundles', n_gt=81049),
    'Bundled_Tubes_Long_Sequence': dict(size=64, pixel_nm=100.0, n_frames=12000, read_sd=5.0, z_double_nm=300.0, structure='bundles', n_gt=81049),
}
PSF_FWHM_NM = 258.21                 # page: XY function Gaussian, FWHM 258.21 nm (= wavelength / (2 NA) = 723 / 2.8)
FWHM_TO_SIGMA = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))

# the declared development / test frame ranges, [first, last + 1), 0-based (see the module docstring for the reasons)
PROTOCOL = {
    'Tubulins_I': dict(dev=(0, 40), test=(200, 300)),
    'Tubulins_II': dict(dev=(0, 40), test=(200, 300)),
    'Bundled_Tubes_High_Density': dict(dev=(0, 30), test=(100, 160)),
    'Bundled_Tubes_Long_Sequence': dict(dev=(0, 300), test=(400, 900)),
}
CAL_FROM = {'Bundled_Tubes_High_Density': 'Bundled_Tubes_Long_Sequence'}     # a dense set borrows the camera / PSF calibration of the sparse set of the same structure


def spec(name):
    if name not in SPEC:
        raise KeyError('unknown set %r (known: %s)' % (name, ', '.join(SETS)))
    return SPEC[name]


def field_nm(name):
    """Side of the (square) field of view in nm: 38,400 for the tubulin sets, 6,400 for the bundle sets."""
    s = spec(name)
    return s['size'] * s['pixel_nm']


def sigma_page_nm(z_nm, name):
    """The PSF width (nm) a reading of the page predicts at axial position z: sigma = FWHM / 2.3548 at the focal plane, growing exponentially to twice that at ``z_double_nm``
    (the page's 'Z function: Exponential, focal plane 1 x FWHM at 0 nm, defocus plane 2 x FWHM at 500 / 300 nm'). A READING of the page, to be compared with the measured widths."""
    return PSF_FWHM_NM * FWHM_TO_SIGMA * 2.0 ** (np.abs(np.asarray(z_nm, dtype=np.float64)) / spec(name)['z_double_nm'])


# ----------------------------------------------------------------------------------------------------------------------------------
# data access
# ----------------------------------------------------------------------------------------------------------------------------------
def data_path(name, kind='stack', root=None):
    """Path of a downloaded file: kind 'stack' (``sequence.zip``, one TIFF per frame), 'gt' (``fluorophores.zip``, holding ``activation.csv``), 'oracle', 'psf'."""
    spec(name)
    root = root or os.environ.get('SMLM_DATA_2013', DEFAULT_DATA)
    return os.path.join(root, name, {'stack': 'sequence.zip', 'gt': 'fluorophores.zip', 'oracle': 'oracle.zip', 'psf': 'psf.zip'}[kind])


def data_available(name, root=None):
    return os.path.exists(data_path(name, 'stack', root)) and os.path.exists(data_path(name, 'gt', root))


_FRAME_RE = re.compile(r'^(?:.*/)?(\d+)\.tiff?$', re.IGNORECASE)


def frame_members(zf):
    """{frame number: member name} of the per-frame TIFFs of an open zip (the integer of the file name; directories and other members ignored)."""
    out = {}
    for i in zf.infolist():
        m = _FRAME_RE.match(i.filename)
        if m and not i.is_dir():
            out[int(m.group(1))] = i.filename
    return out


def n_frames_in(path):
    """Number of frames of a zip of per-frame TIFFs; the numbers must be 1..n without a gap (else ValueError)."""
    with zipfile.ZipFile(path) as z:
        nums = sorted(frame_members(z))
    if not nums or nums != list(range(GT_FRAME_BASE, GT_FRAME_BASE + len(nums))):
        raise ValueError('the frame numbers of %s are not %d..n without a gap' % (path, GT_FRAME_BASE))
    return len(nums)


def load_frames(path, start=0, stop=None):
    """Frames ``start`` .. ``stop - 1`` (0-based, stop exclusive; None: to the end) of a zip of per-frame TIFFs (frame k = the TIFF numbered k + 1) as float64 (n, rows, cols) in DN.

    Each member is decoded in memory with Pillow (nothing is extracted); 16-bit big-endian TIFFs (Pillow mode I;16B) are handled. Raises IndexError when the archive has fewer than
    ``stop`` frames (no silent truncation), ValueError for a gap in the numbering or frames of different shapes."""
    from PIL import Image
    if start < 0 or (stop is not None and stop < start):
        raise ValueError('need 0 <= start <= stop, got start=%r stop=%r' % (start, stop))
    with zipfile.ZipFile(str(path)) as z:
        mem = frame_members(z)
        nums = sorted(mem)
        if not nums or nums != list(range(GT_FRAME_BASE, GT_FRAME_BASE + len(nums))):
            raise ValueError('the frame numbers of %s are not %d..n without a gap' % (path, GT_FRAME_BASE))
        n = len(nums)
        stop = n if stop is None else stop
        if stop > n:
            raise IndexError('stack has only %d frames, requested stop=%d' % (n, stop))
        out = []
        for k in range(start, stop):
            a = np.asarray(Image.open(io.BytesIO(z.read(mem[k + GT_FRAME_BASE]))))
            if a.ndim != 2:
                raise ValueError('expected single-channel frames, got array of shape %r' % (a.shape,))
            out.append(a.astype(np.float64))
    if not out:
        return np.empty((0, 0, 0))
    if any(o.shape != out[0].shape for o in out):
        raise ValueError('frames of different shapes in %s' % path)
    return np.stack(out)


def to_photons(frames):
    """Photons per pixel from DN: (DN - 100) / 1 (the pages: baseline 100 DN, QE 1.00 e-/photon, 1.00 DN/e-)."""
    return (np.asarray(frames, dtype=np.float64) - BASELINE_DN) / DN_PER_PHOTON


def parse_gt(text, frame_base=GT_FRAME_BASE):
    """The ground truth from the text of an ``activation.csv`` as the harness's structured array (``smlm_common.GT_DTYPE``), sorted by (frame, id): fields id, frame (0-BASED = file
    frame - ``frame_base``), x_nm, y_nm (corner convention, x column axis), z_nm, photons. The header must be ``Ground-truth,frame,xnano,ynano,znano,intensity`` (that of the 2016
    files). Same rules as ``smlm_common.load_gt``."""
    lines = text.splitlines()
    head = [s.strip().lower() for s in lines[0].split(',')] if lines else []
    if head != C._GT_HEADER:
        raise ValueError('unexpected ground-truth header %r (expected %r)' % (head, C._GT_HEADER))
    raw = np.loadtxt(io.StringIO('\n'.join(lines[1:])), delimiter=',', ndmin=2) if len(lines) > 1 else np.empty((0, 6))
    gt = np.empty(len(raw), dtype=C.GT_DTYPE)
    gt['id'] = raw[:, 0].astype(np.int64)
    fr = raw[:, 1]
    if not np.all(fr == np.round(fr)):
        raise ValueError('non-integer frame numbers')
    gt['frame'] = fr.astype(np.int64) - frame_base
    gt['x_nm'], gt['y_nm'], gt['z_nm'], gt['photons'] = raw[:, 2], raw[:, 3], raw[:, 4], raw[:, 5]
    if len(gt) and gt['frame'].min() < 0:
        raise ValueError('negative frame index after the base shift (frame_base=%r)' % frame_base)
    return gt[np.lexsort((gt['id'], gt['frame']))]


def load_gt(path, frame_base=GT_FRAME_BASE):
    """Ground truth of a set from ``fluorophores.zip`` (its ``activation.csv`` member, wherever it sits in the archive) or from a plain csv: see ``parse_gt``."""
    path = str(path)
    if path.lower().endswith('.zip'):
        with zipfile.ZipFile(path) as z:
            cands = [i.filename for i in z.infolist() if not i.is_dir() and os.path.basename(i.filename) == GT_MEMBER]
            if len(cands) != 1:
                raise ValueError('expected exactly one %s in %s, found %d' % (GT_MEMBER, path, len(cands)))
            text = z.read(cands[0]).decode('utf-8')
    else:
        with open(path) as f:
            text = f.read()
    return parse_gt(text, frame_base)


def gt_frame_files(path):
    """The per-frame ground-truth files of a ``fluorophores.zip`` (``frames/NNNNN.csv``) as one structured array (same fields, sorted by (frame, id))), for the cross-check against
    ``activation.csv``: the file ids restart at 1 in every frame, so only (frame, x, y, z, photons) are comparable."""
    parts = []
    with zipfile.ZipFile(str(path)) as z:
        for i in sorted(z.infolist(), key=lambda i: i.filename):
            if not i.is_dir() and re.search(r'(^|/)frames/\d+\.csv$', i.filename):
                g = parse_gt(z.read(i.filename).decode('utf-8'))
                if len(g):
                    parts.append(g)
    return np.concatenate(parts) if parts else np.zeros(0, dtype=C.GT_DTYPE)


def gt_in_frames(gt, a, b):
    """The ground-truth rows with frame in [a, b)."""
    return gt[(gt['frame'] >= a) & (gt['frame'] < b)]


# ----------------------------------------------------------------------------------------------------------------------------------
# scoring with the field of the set
# ----------------------------------------------------------------------------------------------------------------------------------
@contextlib.contextmanager
def field(name):
    """Within the block ``smlm_common.FIELD_NM`` is the field of the set (the harness uses it only for the border exclusion of the challenge setting) and is restored afterwards."""
    old = C.FIELD_NM
    C.FIELD_NM = field_nm(name)
    try:
        yield
    finally:
        C.FIELD_NM = old


def score(detections, gt, frame_indices, name, **kw):
    """``smlm_common.score`` with the field of the set (so that ``border_nm`` > 0 drops the emitters within that distance of the true border of the field); every argument and
    result as there."""
    with field(name):
        return C.score(detections, gt, frame_indices, **kw)


# ----------------------------------------------------------------------------------------------------------------------------------
# convention checks against the images
# ----------------------------------------------------------------------------------------------------------------------------------
def excess_at_ground_truth(frames_dn, gt, first_frame, offsets=(-2, -1, 0, 1, 2), pixel_nm=C.PIXEL_NM, half=1, min_photons=0.0):
    """Mean excess (DN) of the (2 half + 1)^2 pixels around the ground-truth positions over the median of their frame, as a function of the offset between the stack index and the
    ground-truth frame: for offset o the ground-truth rows of frame f are looked up in image f + o. A wrong frame base or a wrong axis mapping leaves no excess (the right
    mapping stands out sharply, see ``results/smlm/calibration.md`` section 0 for the 2016 sets). frames_dn: (n, rows, cols) DN of the frames ``first_frame`` .. ``first_frame + n - 1``;
    the pixel containing the emitter under the corner convention is (floor(y / p), floor(x / p)). Returns {offset: (mean excess, number of rows)}."""
    n, ny, nx = frames_dn.shape
    med = np.median(frames_dn.reshape(n, -1), axis=1)
    out = {}
    for o in offsets:
        vals = []
        for f in range(n):
            g = gt[(gt['frame'] == first_frame + f - o) & (gt['photons'] >= min_photons)]
            for x, y in zip(g['x_nm'], g['y_nm']):
                ix, iy = int(np.floor(x / pixel_nm)), int(np.floor(y / pixel_nm))
                if iy < half or iy >= ny - half or ix < half or ix >= nx - half:
                    continue
                vals.append(frames_dn[f, iy - half:iy + half + 1, ix - half:ix + half + 1].mean() - med[f])
        out[o] = (float(np.mean(vals)) if vals else float('nan'), len(vals))
    return out


def check_conventions(name, dets, gt, frame_indices):
    """Evidence for the position and frame conventions from a localization run (e.g. Picasso at a fixed low threshold on a few frames; ``dets`` per-frame (n, 2+) arrays in the corner
    convention as the harness produces them) against the ground truth: returns a dict of
      corner   (tp, bias_x_nm, bias_y_nm)  the matched count and the mean signed error under the harness's corner convention (the files' convention: bias near 0),
      centre   the same with the detections shifted by -p/2 (the convention in which the pixel centre is at j p: a bias of -p/2 in each axis),
      swapped  the matched count with x and y exchanged, flip_x / flip_y with an axis flipped (all collapse),
      frames   {offset: tp} for the ground-truth frame shifted by -2 .. +2 (the true base 0 stands out)."""
    p = spec(name)['pixel_nm']
    idx = np.asarray(frame_indices)

    def run(ds, g=gt):
        r = score(ds, g, idx, name)
        return dict(tp=r['tp'], bias_x_nm=r['bias_x_nm'], bias_y_nm=r['bias_y_nm'], rmse_nm=r['rmse_nm'])

    out = dict(corner=run(dets))
    out['centre'] = run([np.c_[d[:, 0] - 0.5 * p, d[:, 1] - 0.5 * p] if len(d) else d for d in dets])
    out['swapped'] = run([np.c_[d[:, 1], d[:, 0]] if len(d) else d for d in dets])
    fld = field_nm(name)
    out['flip_x'] = run([np.c_[fld - d[:, 0], d[:, 1]] if len(d) else d for d in dets])
    out['flip_y'] = run([np.c_[d[:, 0], fld - d[:, 1]] if len(d) else d for d in dets])
    fr = {}
    for o in (-2, -1, 0, 1, 2):
        g2 = gt.copy()
        g2['frame'] = g2['frame'] + o
        fr[o] = score(dets, g2, idx, name)['tp']
    out['frames'] = fr
    return out
