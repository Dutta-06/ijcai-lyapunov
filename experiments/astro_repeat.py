"""
Repeat imaging of M2: the four M2 test windows in other SDSS runs (Stripe 82), at other seeing, with the HUGS truth mapped into each frame
======================================================================================================================================

The sliding inverse scale space was tuned on two development patches and tested on four test windows of one frame of M2 (run 2583, camcol 2, field 136, psField FWHM about 2.1 px), then run
on M53 and M92 (FWHM 3.3 and 2.9 px) where it did not transfer. M2 lies on Stripe 82, which SDSS scanned many times, so the SAME stars (the same four windows, the same HUGS truth) are available
in other runs at other seeing. This module is the geometry, the selection of the frames, the patches and their validation; it runs no cataloguing method (the baselines are in
`exp_astro_baselines_repeat.py`). The findings, the declared selection rule and the numbers are in `results/astro/repeat_spec.md`.

Conventions are those of `astro_second.py` (pixel (x, y) = (column, row), 0-based, pixel centre at the integer; the box of a patch with origin (x0, y0) is [x0 - 0.5, x0 + 99.5)^2; truth = HUGS
F606W used as r; the full asTrans solution maps the sky to the frame; DN, electrons, PSF as `astro_second.extract_patch`). The window of a repeat frame is the 100 x 100 box centred (to 0.5 px) on
the sky position of the centre of an M2 test window.

Selection rule (declared before any repeat frame was looked at, refined twice before any patch was built; `results/astro/repeat_spec.md` sections 1, 2, 2b): a window is ELIGIBLE in a frame when
rules (1) and (3) of `astro_second` hold (40 px from the frame edges; no pixel above 0.8 x the photoField r saturation level within 10 px; the other three rules, footprint, bright truth star and
lost bright star, are computed and reported; rule (4) fails for every window in every frame, so no window passes all five); a frame is eligible with at least two eligible windows; its seeing is the mean FWHM of the psField stamp at its eligible windows; the frames chosen are,
by closest pair first (the pair of a target and a frame at the smallest distance quantised to 0.05 px, ties by the larger number of windows, then the smaller run, camcol, field, target),
for the targets 2.7, 3.1, 3.5 px. `select_frames` is that rule.

Command line (needs astropy; run it as `python -I` because it reads downloaded, untrusted files; $ASTRO_DATA_3 holds the output of `data/astro/fetch_repeat.py fetch`):
    python -I experiments/astro_repeat.py survey [--stage runs|geometry|frames|all]    the survey tables and the selection (stage `frames` needs the candidate frames on disk)
    python -I experiments/astro_repeat.py build                                         the truth lists and the patches of the frames of REPEAT
    python -I experiments/astro_repeat.py validate                                      the validation numbers of section 4 of the spec for every patch of REPEAT and for the M2 windows
"""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, HERE)
import astro_second as S      # noqa: E402

ROOT3 = os.environ.get('ASTRO_DATA_3', '/home/user/astro_data_3')
ROOT = os.environ.get('ASTRO_DATA', '/home/user/astro_data')                 # the M2 patches and truth (data/astro/fetch.py)
M2_FRAME = (2583, 2, 136)                                                    # run, camcol, field of the benchmark frame
M2_WINDOWS = ((310, 630), (310, 730), (360, 1040), (100, 1100))              # the M2 test windows (lower-left origins in the frame of M2_FRAME)
N_PATCH = S.N_PATCH
TARGETS_FWHM = (2.7, 3.1, 3.5)                                               # px, the FWHM of the psField stamp the three frames should have
QUANTUM = 0.05                                                               # px, the quantum of the distance to the target (ties below it)
MIN_WINDOWS = 2                                                              # eligible windows a frame needs
RUN_HALF_WIDTH = 1.5                                                         # deg, `astro_second.run_candidates`
GEOM_PAD = 0.06                                                              # deg, the prefilter around a photoField box (about 550 px)
GEOM_MARGIN = 30                                                             # px, frames with a window at least this far from the edge (photoField coefficients) are examined
ELIGIBLE_RULES = ('margin', 'saturation')                                   # the conditions of rule R3 as refined before building (spec section 2): the two that depend on the frame; all five of
                                                                             # `astro_second.RULES` are computed and reported ('eligible_literal' = all five, which no window satisfies in any frame)
SURVEY_JSON = os.path.join(ROOT3, 'repeat_survey.json')

# the frames chosen by rule R5' (filled in below from the output of `survey`); names in order of increasing seeing
REPEAT = {}


def _hugs_path():
    return os.path.join(ROOT3, 'ngc7089', 'hugs_ngc7089', 'hlsp_hugs_hst_wfc3-uvis-acs-wfc_ngc7089_multi_v1_catalog-meth1.txt')


def frame_paths(run, camcol, field):
    """Where `data/astro/fetch_repeat.py` puts the frame, the psField and the photoField of (run, camcol, field)."""
    if (run, camcol, field) == M2_FRAME:
        d = lambda key: os.path.join(ROOT3, 'ngc7089', key)
        return dict(frame=os.path.join(d('sdss_frame_r_002583_2_0136'), 'frame-r-002583-2-0136.fits.bz2'), psfield=os.path.join(d('sdss_psfield_002583_2_0136'), 'psField-002583-2-0136.fit'),
                    photofield=photofield_path(run, camcol))
    return dict(frame=os.path.join(ROOT3, 'frames', 'sdss_frame_r_%06d_%d_%04d' % (run, camcol, field), 'frame-r-%06d-%d-%04d.fits.bz2' % (run, camcol, field)),
                psfield=os.path.join(ROOT3, 'frames', 'sdss_psfield_%06d_%d_%04d' % (run, camcol, field), 'psField-%06d-%d-%04d.fit' % (run, camcol, field)), photofield=photofield_path(run, camcol))


def photofield_path(run, camcol):
    if (run, camcol) == M2_FRAME[:2]:
        return os.path.join(ROOT3, 'ngc7089', 'sdss_photofield_002583_2', 'photoField-002583-2.fits')
    return os.path.join(ROOT3, 'survey', 'sdss_photofield_%06d_%d' % (run, camcol), 'photoField-%06d-%d.fits' % (run, camcol))


def runall_path():
    return os.path.join(ROOT3, 'survey', 'photorunall', 'photoRunAll-dr17.fits')


# ----------------------------------------------------------------------------- the windows
def window_centre(x0, y0, n=N_PATCH):
    """The pixel position of the centre of the box with origin (x0, y0): (x0 + (n - 1) / 2, y0 + (n - 1) / 2)."""
    return x0 + (n - 1) / 2.0, y0 + (n - 1) / 2.0


def window_sky(astrans, windows=M2_WINDOWS, n=N_PATCH):
    """(RA, Dec) in degrees of the centres of the windows of a frame (rule R0): the asTrans solution of that frame, `astro_second.astrans_pix2world`."""
    c = np.array([window_centre(x0, y0, n) for x0, y0 in windows])
    ra, dec = S.astrans_pix2world(c[:, 0], c[:, 1], astrans)
    return np.atleast_1d(ra), np.atleast_1d(dec)


def origin_of_centre(xc, yc, n=N_PATCH):
    """The integer origin of the box centred on the pixel position (xc, yc) (rule R2): floor(xc - (n - 1) / 2 + 0.5), floor(yc - (n - 1) / 2 + 0.5)."""
    return int(np.floor(xc - (n - 1) / 2.0 + 0.5)), int(np.floor(yc - (n - 1) / 2.0 + 0.5))


def window_pixels(ra, dec, astrans, n=N_PATCH):
    """The windows of a frame: list of (xc, yc, x0, y0), (xc, yc) the pixel position of each sky position under the frame's asTrans solution and (x0, y0) the integer origin of the box."""
    xc, yc = S.astrans_world2pix(ra, dec, astrans)
    return [(float(a), float(b)) + origin_of_centre(a, b, n) for a, b in zip(np.atleast_1d(xc), np.atleast_1d(yc))]


def edge_margin(x0, y0, shape=S.FRAME_SHAPE, n=N_PATCH):
    """Pixels between the box with origin (x0, y0) and the nearest frame edge: min(x0, y0, columns - (x0 + n), rows - (y0 + n)); negative when the box leaves the frame."""
    return int(min(x0, y0, shape[1] - (x0 + n), shape[0] - (y0 + n)))


def frame_orientation(t, x, y, h=10.0):
    """Pixel scales (arcsec per px) and position angles (deg, east of north) of the +x (column) and +y (row) axes of a frame at pixel (x, y), from the asTrans solution `t`:
    [(scale_x, pa_x), (scale_y, pa_y)]. In the M2 frame +x points north (pa 0) and +y east (pa 90)."""
    ra0, dec0 = S.astrans_pix2world(x, y, t)
    out = []
    for dx, dy in ((h, 0.0), (0.0, h)):
        ra1, dec1 = S.astrans_pix2world(x + dx, y + dy, t)
        d_east, d_north = float(ra1 - ra0) * np.cos(np.deg2rad(float(dec0))) * 3600.0, float(dec1 - dec0) * 3600.0
        out.append((float(np.hypot(d_east, d_north) / h), float(np.rad2deg(np.arctan2(d_east, d_north)))))
    return out


# ----------------------------------------------------------------------------- the eligibility of a window in a frame (rules (1) to (5) of second_cluster_spec.md section 2.4)
def frame_context(truth, missing, dn_full, sat_level, shape=S.FRAME_SHAPE):
    """What the rules need of one frame: the HUGS footprint occupancy of its truth, the bright and the lost truth sources, the DN image and the photoField saturation level.
    truth: `astro_second.make_truth`; missing: `astro_second.unmeasured_f606w`."""
    bright = truth['m606'] < S.BRIGHT_LIMIT
    lost = missing['f814'] < S.MISSING_LIMIT
    return dict(occ=S.footprint_cells(truth['px'], truth['py'], shape), bx=truth['px'][bright], by=truth['py'][bright], lx=missing['px'][lost], ly=missing['py'][lost],
                dn=dn_full, sat=float(sat_level), shape=shape)


def _any_within(px, py, x0, y0, n, margin=S.DN_MARGIN):
    """True when a source lies within `margin` px of the box [x0 - 0.5, x0 + n - 0.5)^2 (the test of `astro_second.screen_patches`)."""
    return bool(np.any((px >= x0 - margin - 0.5) & (px < x0 + n + margin - 0.5) & (py >= y0 - margin - 0.5) & (py < y0 + n + margin - 0.5)))


def window_rules(ctx, x0, y0, n=N_PATCH):
    """The five rules for the box with origin (x0, y0): dict(margin_px, margin, footprint, saturation, bright, lost, dn_max, eligible_literal, eligible), the booleans True = passes;
    `dn_max` is the largest DN within 10 px of the box (NaN if that region is empty). The same tests as `astro_second.screen_patches`, for an arbitrary origin. `eligible_literal` needs all
    five; `eligible` needs `ELIGIBLE_RULES` (margin and saturation), the condition of rule R3 as refined before building."""
    m = edge_margin(x0, y0, ctx['shape'], n)
    sub = ctx['dn'][max(y0 - S.DN_MARGIN, 0): y0 + n + S.DN_MARGIN, max(x0 - S.DN_MARGIN, 0): x0 + n + S.DN_MARGIN]
    dmax = float(sub.max()) if sub.size else float('nan')
    r = dict(margin_px=m, margin=bool(m >= S.FRAME_MARGIN), footprint=bool(S.inside_footprint(ctx['occ'], x0, y0, n)), saturation=bool(sub.size and dmax < S.SAT_FRACTION * ctx['sat']),
             bright=not _any_within(ctx['bx'], ctx['by'], x0, y0, n), lost=not _any_within(ctx['lx'], ctx['ly'], x0, y0, n), dn_max=dmax)
    r['eligible_literal'] = all(r[k] for k in S.RULES)
    r['eligible'] = all(r[k] for k in ELIGIBLE_RULES)
    return r


def window_fwhm(kl, x0, y0, n=N_PATCH):
    """FWHM (px) of the psField r-band stamp (unit-sum 51 x 51 reconstruction) at the centre of the box with origin (x0, y0), as the patch's `psf51`."""
    xc, yc = window_centre(x0, y0, n)
    return S.psf_fwhm(S.unit_psf(S.sdss_psf_recon(kl, xc, yc)))


# ----------------------------------------------------------------------------- the selection (rule R5)
def quantised_distance(fwhm, target, quantum=QUANTUM):
    """floor(|fwhm - target| / quantum + 0.5): the distance to the target in units of `quantum` (two frames closer to each other than about one quantum are tied)."""
    return int(np.floor(abs(float(fwhm) - float(target)) / quantum + 0.5))


def select_key(frame, target, quantum=QUANTUM):
    """The sort key of a frame (dict with run, camcol, field, n_w, fwhm) for one target: quantised distance, then more eligible windows, then smaller run, camcol, field."""
    return (quantised_distance(frame['fwhm'], target, quantum), -int(frame['n_w']), int(frame['run']), int(frame['camcol']), int(frame['field']))


def eligible_frames(frames, min_windows=MIN_WINDOWS, exclude_runs=(M2_FRAME[0],)):
    """The frames (dicts with run, camcol, field, n_w, fwhm) with at least `min_windows` eligible windows and a finite FWHM, run 2583 excluded."""
    return [f for f in frames if f['n_w'] >= min_windows and f['run'] not in exclude_runs and np.isfinite(f['fwhm'])]


def rank_frames(frames, target, quantum=QUANTUM, min_windows=MIN_WINDOWS, exclude_runs=(M2_FRAME[0],)):
    """The eligible frames in the order of preference for one target."""
    return sorted(eligible_frames(frames, min_windows, exclude_runs), key=lambda f: select_key(f, target, quantum))


def _frame_id(f):
    return (int(f['run']), int(f['camcol']), int(f['field']))


def select_frames(frames, targets=TARGETS_FWHM, quantum=QUANTUM, min_windows=MIN_WINDOWS, exclude_runs=(M2_FRAME[0],)):
    """Rule R5' (closest pair first): among all pairs of a target not matched yet and an eligible frame not chosen yet, the pair with the smallest (`select_key`, target) is matched, and so on until
    every target has a frame or the frames run out. Returns a list of dict(target, q, frame) in the order of matching (the closest pair first)."""
    pool = eligible_frames(frames, min_windows, exclude_runs)
    out, used, free = [], set(), list(targets)
    while free:
        pairs = [(select_key(f, T, quantum), T, f) for T in free for f in pool if _frame_id(f) not in used]
        if not pairs:
            break
        key, T, f = min(pairs, key=lambda p: (p[0], p[1]))
        used.add(_frame_id(f))
        free.remove(T)
        out.append(dict(target=T, q=key[0], frame=f))
    return out


def select_frames_target_order(frames, targets=TARGETS_FWHM, quantum=QUANTUM, min_windows=MIN_WINDOWS, exclude_runs=(M2_FRAME[0],)):
    """Rule R5 as first declared and superseded before building (spec section 2b): the targets in the given order, each taking the best frame left. Kept for the comparison of the survey output."""
    out, used = [], set()
    for T in targets:
        left = [f for f in rank_frames(frames, T, quantum, min_windows, exclude_runs) if _frame_id(f) not in used]
        if not left:
            break
        used.add(_frame_id(left[0]))
        out.append(dict(target=T, q=quantised_distance(left[0]['fwhm'], T, quantum), frame=left[0]))
    return out


def frame_names(chosen):
    """{'m2r1': frame, 'm2r2': ..} in order of increasing FWHM of the chosen frames."""
    fr = sorted((c['frame'] for c in chosen), key=lambda f: (f['fwhm'], f['run'], f['camcol'], f['field']))
    return {'m2r%d' % (i + 1): f for i, f in enumerate(fr)}


# ----------------------------------------------------------------------------- the registry and the loaders (numpy only)
def frame_tag(run, camcol, field):
    return '%06d_%d_%04d' % (run, camcol, field)


def registry_entry(title, run, camcol, field, patches, windows, root=None):
    """An entry of REPEAT: the patches (x0, y0) of one frame, the M2 test window each belongs to, the directory prefix and the truth file."""
    return dict(title=title, root=root or ROOT3, prefix='patch_r_' + frame_tag(run, camcol, field), truth='truth_hugs_ngc7089_%s_meth1.npz' % frame_tag(run, camcol, field),
                patches=[tuple(p) for p in patches], windows=[tuple(w) for w in windows], frame=(run, camcol, field))


# The frames chosen by rule R5' (`survey`, spec section 3), named in order of increasing seeing; patches = the integer origins of the eligible windows in that frame, windows = the M2 test window of each.
REPEAT.update({
    'm2r1': registry_entry('M2 repeat 1 (NGC 7089), SDSS run 4187 camcol 2 field 105, 2003-10-18, FWHM 2.28 px', 4187, 2, 105, [(316, 124), (316, 224)], [(310, 630), (310, 730)]),
    'm2r2': registry_entry('M2 repeat 2 (NGC 7089), SDSS run 7778 camcol 6 field 72, 2008-10-30, FWHM 3.05 px', 7778, 6, 72, [(133, 281), (195, 689)], [(310, 630), (360, 1040)]),
    'm2r3': registry_entry('M2 repeat 3 (NGC 7089), SDSS run 3388 camcol 2 field 75, 2002-10-11, FWHM 3.52 px', 3388, 2, 75, [(370, 441), (110, 501)], [(360, 1040), (100, 1100)]),
})


def patch_path(x0, y0, name):
    """patch_arrays.npz of the patch with origin (x0, y0) of the repeat frame `name` (REPEAT)."""
    r = REPEAT[name]
    return os.path.join(r['root'], '%s_x%d_y%d_100x100' % (r['prefix'], x0, y0), 'patch_arrays.npz')


def load_patch(x0, y0, name):
    """dict(dn, psf, gain, dark) exactly like `exp_astro_m2.load_patch`."""
    z = np.load(patch_path(x0, y0, name))
    return dict(dn=z['dn'].astype(float), psf=z['psf51'].astype(float), gain=float(z['gain']), dark=float(z['dark_variance']))


def truth_arrays(name):
    """dict(px, py, m606) of the truth list of the repeat frame `name`, exactly like `exp_astro_m2.truth_arrays` (read it with `astro_common.truth_for_patch(truth, x0, y0)`)."""
    z = np.load(os.path.join(REPEAT[name]['root'], REPEAT[name]['truth']))
    return dict(px=z['px'], py=z['py'], m606=z['m606'])


# ----------------------------------------------------------------------------- the survey (astropy): runs, photoField geometry, frames
def read_runs():
    """The photoRunAll-dr17 table as dict of arrays (the columns `astro_second.run_candidates` needs, plus the date)."""
    from astropy.io import fits
    t = fits.getdata(runall_path(), 1)
    return {k: np.asarray(t[k]) for k in ('RUN', 'NODE', 'INCL', 'MU_REF', 'FIELD_REF', 'STARTFIELD', 'ENDFIELD', 'STRIPE', 'STRIP', 'DATESTRING')}


def read_photofield(run, camcol):
    from astropy.io import fits
    return fits.getdata(photofield_path(run, camcol), 1)


def _rband(pf, col, k):
    v = np.asarray(pf[col][k]).ravel()
    return float(v[2]) if v.size == 5 else float(v[0])


def photofield_astrans(pf, k, node, incl):
    """The asTrans constants of row k of a photoField table, r band, as the dict of `astro_second.astrans_params`. The table stores them in single precision (positions good to about
    0.1 px); the downloaded frame carries the double-precision row and decides in the final tables."""
    t = dict(NODE=float(node), INCL=float(incl))
    for key in S.ASTRANS_KEYS[2:]:
        t[key] = _rband(pf, 'FF' if key == 'F' else key, k)
    return t


def photofield_record(pf, k):
    """The r-band entries of row k of a photoField table that the survey reports."""
    return dict(psf_width=_rband(pf, 'PSF_WIDTH', k), psp_status=int(pf['PSP_STATUS'][k]), psf_nstar=int(_rband(pf, 'PSF_NSTAR', k)), sat_level=_rband(pf, 'SATURATION_LEVEL', k), quality=int(pf['QUALITY'][k]),
                nmgy=_rband(pf, 'NMGYPERCOUNT', k), gain=_rband(pf, 'GAIN', k), dark_variance=_rband(pf, 'DARK_VARIANCE', k), sky=_rband(pf, 'SKY', k))


def candidate_runs(runs, ra, dec, half_width=RUN_HALF_WIDTH):
    """Rule R1, first part: the runs of photoRunAll whose scan circle passes within `half_width` deg of any of the window positions and whose fields cover them
    (`astro_second.run_candidates`): list of dict(run, strip, stripe, node, incl, date, nu (per window), field (per window)), sorted by run."""
    found = {}
    for a, d in zip(ra, dec):
        for r in S.run_candidates(runs, a, d, half_width):
            found.setdefault(r['run'], []).append(r)
    out = []
    for run in sorted(found):
        i = int(np.flatnonzero(runs['RUN'] == run)[0])
        out.append(dict(run=run, strip=str(runs['STRIP'][i]).strip(), stripe=int(runs['STRIPE'][i]), node=float(runs['NODE'][i]), incl=float(runs['INCL'][i]), date=str(runs['DATESTRING'][i]).strip(),
                        nu=[r['nu'] for r in found[run]], field=[r['field'] for r in found[run]]))
    return out


def geometry_candidates(cands, ra, dec, pad=GEOM_PAD, min_margin=GEOM_MARGIN):
    """Rule R1, second part, from the photoField tables alone: for every candidate run and camcol the fields whose frame (the photoField r-band asTrans row) holds the centre of at least one
    window, with the window pixels and margins. Returns (frames, missing photoField paths); frames = list of dict(run, camcol, field, strip, date, windows=[dict(window, xc, yc, x0, y0,
    margin_px)]) for the frames in which some window has margin_px >= min_margin."""
    out, missing = [], []
    for c in cands:
        mu, nu = S.eq_to_gc(ra, dec, c['node'], c['incl'])
        for cam in range(1, 7):
            if not os.path.exists(photofield_path(c['run'], cam)):
                missing.append(photofield_path(c['run'], cam))
                continue
            pf = read_photofield(c['run'], cam)
            m0, m1 = np.asarray(pf['MU_START'], float), np.asarray(pf['MU_END'], float)
            n0, n1 = np.asarray(pf['NU_START'], float), np.asarray(pf['NU_END'], float)
            lo, hi, span = np.minimum(n0, n1), np.maximum(n0, n1), (m1 - m0 + 180.0) % 360.0 - 180.0
            near = np.zeros(len(m0), bool)
            for w in range(len(ra)):
                dm = (mu[w] - m0 + 180.0) % 360.0 - 180.0
                near |= (dm > -pad) & (dm < span + pad) & (nu[w] > lo - pad) & (nu[w] < hi + pad)
            for k in np.flatnonzero(near):
                t = photofield_astrans(pf, k, c['node'], c['incl'])
                wins = []
                for w, (xc, yc, x0, y0) in zip(M2_WINDOWS, window_pixels(ra, dec, t)):
                    if 0 <= xc < S.FRAME_SHAPE[1] and 0 <= yc < S.FRAME_SHAPE[0]:
                        wins.append(dict(window=list(w), xc=xc, yc=yc, x0=x0, y0=y0, margin_px=edge_margin(x0, y0)))
                if any(q['margin_px'] >= min_margin for q in wins):
                    out.append(dict(run=c['run'], camcol=cam, field=int(pf['FIELD'][k]), strip=c['strip'], date=c['date'], windows=wins))
    return out, missing


def _jsonable(o):
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, (np.floating, float)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return o


def evaluate_frame(run, camcol, field, ra, dec, hugs, strip='?'):
    """The survey record of one downloaded frame (rules R2 to R4): the photoField r-band entries, the header date, the window pixels from the frame's own asTrans, the five rules and the psField FWHM
    for every window whose centre is inside the frame, the number of windows inside with margin (`n_in`) and eligible (`n_w`), the frame's FWHM (mean over the eligible windows; NaN if none), and the
    largest disagreement (px) between the frame's asTrans and the photoField row at the window centres."""
    P = frame_paths(run, camcol, field)
    fr, ps = S.read_frame(P['frame']), S.read_psfield(P['psfield'])
    pf = read_photofield(run, camcol)
    k = int(np.flatnonzero(np.asarray(pf['FIELD']) == field)[0])
    rec = photofield_record(pf, k)
    w2p = lambda a, d: S.astrans_world2pix(a, d, fr['astrans'])
    truth, miss = S.make_truth(hugs, w2p), S.unmeasured_f606w(hugs, w2p)
    sky = S.sky_image(fr['allsky'], fr['xinterp'], fr['yinterp'])
    ctx = frame_context(truth, miss, S.frame_to_dn(fr['nmgy'], fr['cal'], sky), rec['sat_level'])
    pix = window_pixels(ra, dec, fr['astrans'])
    tpf = photofield_astrans(pf, k, fr['astrans']['NODE'], fr['astrans']['INCL'])
    xp, yp = S.astrans_world2pix(ra, dec, tpf)
    wins = []
    for w, (xc, yc, x0, y0) in zip(M2_WINDOWS, pix):
        if not (0 <= xc < S.FRAME_SHAPE[1] and 0 <= yc < S.FRAME_SHAPE[0]):
            continue
        r = window_rules(ctx, x0, y0)
        wins.append(dict(window=list(w), xc=xc, yc=yc, x0=x0, y0=y0, fwhm=window_fwhm(ps['kl'], x0, y0), **r))
    el = [q for q in wins if q['eligible']]
    n_literal = sum(1 for q in wins if q['eligible_literal'])
    inside = [q for q in wins if q['margin']]
    return dict(run=run, camcol=camcol, field=field, strip=strip, date_obs=str(fr['header']['DATE-OBS']), nmgy_header=float(fr['header']['NMGY']), psfield_gain=ps['gain'], psfield_dark_variance=ps['dark_variance'],
                **rec, windows=wins, n_in=sum(1 for q in wins if q['margin']), n_w=len(el), n_w_literal=n_literal, fwhm=float(np.mean([q['fwhm'] for q in el])) if el else float('nan'), fwhm_in=float(np.mean([q['fwhm'] for q in inside])) if inside else float('nan'),
                astrans_frame_vs_photofield_px=float(max(np.hypot(np.array([q[0] for q in pix]) - xp, np.array([q[1] for q in pix]) - yp))))


def run_coverage(c, ra, dec):
    """For one candidate run (dict of `candidate_runs`) and the first window: the camcol whose frame comes closest to holding the window centre, from the photoField r-band asTrans of the field of that
    camcol nearest to it along the scan. Returns dict(camcol, field, xc, yc, outside_px): the window centre in that frame's pixels and its distance (px) beyond the nearest column edge of the
    frame (negative: inside); used to say why a run does not cover the windows (the window falls in the gap between two camcols)."""
    mu, nu = S.eq_to_gc(ra, dec, c['node'], c['incl'])
    best = None
    for cam in range(1, 7):
        if not os.path.exists(photofield_path(c['run'], cam)):
            continue
        pf = read_photofield(c['run'], cam)
        m0, m1 = np.asarray(pf['MU_START'], float), np.asarray(pf['MU_END'], float)
        dm = (mu[0] - m0 + 180.0) % 360.0 - 180.0
        span = (m1 - m0 + 180.0) % 360.0 - 180.0
        k = int(np.nanargmin(np.abs(dm - span / 2.0)))
        t = photofield_astrans(pf, k, c['node'], c['incl'])
        xc, yc = S.astrans_world2pix(ra[0], dec[0], t)
        xc, yc = float(np.ravel(xc)[0]), float(np.ravel(yc)[0])
        out = max(-xc, xc - (S.FRAME_SHAPE[1] - 1.0))
        if best is None or out < best['outside_px']:
            best = dict(camcol=cam, field=int(pf['FIELD'][k]), xc=xc, yc=yc, outside_px=out)
    return best


def frame_lines(r):
    """The text lines that describe one evaluated frame (the record of `evaluate_frame`)."""
    out = ['run %d camcol %d field %d strip %s %s  seeing %.3f arcsec (%.2f px)  PSP_STATUS %d  PSF_NSTAR(r) %d  saturation %.0f DN  QUALITY %d  NMGY %.8f  gain %.3f dark %.4f' % (
        r['run'], r['camcol'], r['field'], r['strip'], r['date_obs'], r['psf_width'], r['psf_width'] / S.PIX, r['psp_status'], r['psf_nstar'], r['sat_level'], r['quality'], r['nmgy_header'], r['psfield_gain'], r['psfield_dark_variance'])]
    for q in r['windows']:
        out.append('    window (%4d, %4d): centre (%7.2f, %7.2f) origin (%4d, %4d) margin %5d | (2) footprint %-4s (3) saturation %-4s DNmax %6.0f (4) bright %-4s (5) lost %-4s | FWHM %.3f px | %s' % (
            q['window'][0], q['window'][1], q['xc'], q['yc'], q['x0'], q['y0'], q['margin_px'], 'ok' if q['footprint'] else 'FAIL', 'ok' if q['saturation'] else 'FAIL', q['dn_max'],
            'ok' if q['bright'] else 'FAIL', 'ok' if q['lost'] else 'FAIL', q['fwhm'], 'ELIGIBLE' + ('' if q['eligible_literal'] else ' (not under the literal rule)') if q['eligible'] else ('margin FAIL' if not q['margin'] else 'saturation FAIL')))
    out.append('    windows inside with margin: %d ; eligible windows n_w = %d (all five rules: %d) ; FWHM_f = %s px ; frame asTrans against photoField row at the window centres: %.3f px' % (
        r['n_in'], r['n_w'], r['n_w_literal'], ('%.3f' % r['fwhm']) if np.isfinite(r['fwhm']) else 'n/a', r['astrans_frame_vs_photofield_px']))
    return out


def cmd_survey(args):
    stage = args.stage
    fr0 = S.read_frame(frame_paths(*M2_FRAME)['frame'])
    ra, dec = window_sky(fr0['astrans'])
    xh, yh = S.tan_world2pix(ra, dec, S.tan_params(fr0['header']))
    print('R0 windows of the M2 frame (run %d camcol %d field %d): sky position of the box centre under the frame asTrans solution' % M2_FRAME)
    for w, a, d, x_, y_ in zip(M2_WINDOWS, ra, dec, xh, yh):
        xc, yc = window_centre(*w)
        print('    window (%4d, %4d): centre pixel (%.1f, %.1f) -> RA %.6f Dec %+.6f ; the header TAN solution puts that position at (%.3f, %.3f), %.3f px away' % (w[0], w[1], xc, yc, a, d, x_, y_, np.hypot(x_ - xc, y_ - yc)))
    runs = read_runs()
    cands = candidate_runs(runs, ra, dec)
    print('R1 runs of photoRunAll-dr17 (%d runs, rerun 301) whose scan circle passes within %.1f deg of a window and whose fields cover it: %d' % (len(runs['RUN']), RUN_HALF_WIDTH, len(cands)))
    cands_other = [c for c in cands if c['run'] != M2_FRAME[0]]
    geom, missing = geometry_candidates(cands_other, ra, dec)
    covered = {(g['run']) for g in geom}
    print('    run  stripe strip date        nu of window 1 (deg) | covering frames with a window at >= %d px from the edge, or the nearest camcol' % GEOM_MARGIN)
    for c in cands:
        if c['run'] == M2_FRAME[0]:
            print('    %4d %4d   %-3s %s %+8.3f             | the benchmark frame (camcol 2, field 136), not a candidate' % (c['run'], c['stripe'], c['strip'], c['date'][:10], c['nu'][0]))
            continue
        if c['run'] in covered:
            txt = 'camcol-field ' + ', '.join('%d-%d (%d windows)' % (g['camcol'], g['field'], sum(1 for q in g['windows'] if q['margin_px'] >= S.FRAME_MARGIN)) for g in geom if g['run'] == c['run'])
        else:
            b = run_coverage(c, ra, dec)
            txt = 'none: the window falls %.0f px %s the edge of camcol %d (centre x = %.0f in field %d)' % (abs(b['outside_px']), 'beyond' if b['outside_px'] > 0 else 'inside', b['camcol'], b['xc'], b['field']) if b else 'photoField missing'
        print('    %4d %4d   %-3s %s %+8.3f             | %s' % (c['run'], c['stripe'], c['strip'], c['date'][:10], c['nu'][0], txt))
    if missing:
        print('    photoField files missing on disk: %d (run `fetch_repeat.py fetch --group photofield`)' % len(missing))
    if stage == 'runs':
        return
    print('R1 geometric candidates from the photoField tables (frames with a window at >= %d px from the edge): %d' % (GEOM_MARGIN, len(geom)))
    for g in geom:
        print('    %d-%d-%d  strip %s %s : %s' % (g['run'], g['camcol'], g['field'], g['strip'], g['date'][:10], ' | '.join('(%d,%d)->(%d,%d) m%d' % (w['window'][0], w['window'][1], w['x0'], w['y0'], w['margin_px']) for w in g['windows'])))
    print('    FRAMES = %s' % [(g['run'], g['camcol'], g['field']) for g in geom])
    if stage == 'geometry':
        return
    absent = [(g['run'], g['camcol'], g['field']) for g in geom if not (os.path.exists(frame_paths(g['run'], g['camcol'], g['field'])['frame']) and os.path.exists(frame_paths(g['run'], g['camcol'], g['field'])['psfield']))]
    if absent:
        print('frames missing on disk (run `fetch_repeat.py fetch --group frames`): %s' % absent)
    hugs = S.read_hugs(_hugs_path())
    recs = []
    strip_of = {c['run']: c['strip'] for c in cands}
    for run, cam, fld in [M2_FRAME] + [(g['run'], g['camcol'], g['field']) for g in geom if (g['run'], g['camcol'], g['field']) not in absent]:
        rec = evaluate_frame(run, cam, fld, ra, dec, hugs, strip=strip_of.get(run, '?'))
        rec['reference'] = (run, cam, fld) == M2_FRAME
        recs.append(rec)
        for line in frame_lines(rec):
            print(line)
    pool = [r for r in recs if not r['reference']]
    print('\nR4 frames with at least %d eligible windows: %d of %d examined' % (MIN_WINDOWS, len(eligible_frames(pool)), len(pool)))
    print('    | run | camcol | field | strip | date | seeing (arcsec) | PSP_STATUS | PSF_NSTAR(r) | windows with margin | eligible windows | (all five rules) | FWHM_f (px) |')
    for r in pool:
        print('    | %d | %d | %d | %s | %s | %.3f | %d | %d | %d | %d | %d | %s |' % (r['run'], r['camcol'], r['field'], r['strip'], r['date_obs'][:10], r['psf_width'], r['psp_status'], r['psf_nstar'], r['n_in'], r['n_w'], r['n_w_literal'],
                                                                                 ('%.3f' % r['fwhm']) if np.isfinite(r['fwhm']) else 'n/a'))
    print('\nR5 ranking per target (quantised distance q in units of %.2f px, then more windows, then smaller run/camcol/field), first 5:' % QUANTUM)
    for T in TARGETS_FWHM:
        print('    target %.1f px: %s' % (T, '; '.join('%d-%d-%d (FWHM %.3f, n_w %d, q %d)' % (f['run'], f['camcol'], f['field'], f['fwhm'], f['n_w'], quantised_distance(f['fwhm'], T)) for f in rank_frames(pool, T)[:5])))
    chosen = select_frames(pool)
    names = frame_names(chosen)
    print('\nR5\' chosen (closest pair first): %s' % '; '.join('%.1f px -> %d-%d-%d (FWHM_f %.3f, q %d, n_w %d)' % (c['target'], c['frame']['run'], c['frame']['camcol'], c['frame']['field'], c['frame']['fwhm'], c['q'], c['frame']['n_w']) for c in chosen))
    for nm, f in names.items():
        print('    %s = run %d camcol %d field %d (%s), FWHM_f %.3f px, eligible windows %s ; windows contained but not eligible: %s' % (
            nm, f['run'], f['camcol'], f['field'], f['date_obs'], f['fwhm'], [tuple(q['window']) for q in f['windows'] if q['eligible']],
            [(tuple(q['window']), [k for k in S.RULES if not q[k]]) for q in f['windows'] if not q['eligible']]))
    if len(chosen) < len(TARGETS_FWHM):
        print('FEWER THAN %d ELIGIBLE FRAMES: the chosen set is the frames that exist' % len(TARGETS_FWHM))
    print('\nR5 as first declared (targets in order, each taking the closest frame left; superseded, spec section 2b): %s' % '; '.join(
        '%.1f px -> %d-%d-%d (FWHM_f %.3f)' % (c['target'], c['frame']['run'], c['frame']['camcol'], c['frame']['field'], c['frame']['fwhm']) for c in select_frames_target_order(pool)))
    print('\nSensitivity (not used): the same rule R5\' with other eligibility conditions')
    for label, key_n, key_f in (('literal rule (all five rules)', 'n_w_literal', None), ('margin rule alone (the saturation condition dropped)', 'n_in', 'fwhm_in')):
        alt = []
        for r in pool:
            if key_f is None:                                                   # the literal rule has no eligible window anywhere; the FWHM is that of the eligible windows
                alt.append(dict(r, n_w=r[key_n], fwhm=r['fwhm'] if r[key_n] else float('nan')))
            else:
                alt.append(dict(r, n_w=r[key_n], fwhm=r[key_f]))
        ch = select_frames(alt)
        print('    %-52s: %s' % (label, '; '.join('%.1f px -> %d-%d-%d (FWHM_f %.3f, n_w %d)' % (c['target'], c['frame']['run'], c['frame']['camcol'], c['frame']['field'], c['frame']['fwhm'], c['frame']['n_w']) for c in ch) or 'no eligible frame'))
    if not args.no_json:
        with open(SURVEY_JSON, 'w') as f:
            json.dump(_jsonable(dict(windows=[list(w) for w in M2_WINDOWS], window_ra=ra.tolist(), window_dec=dec.tolist(), runs=cands, frames=recs, chosen=[dict(target=c['target'], q=c['q'], run=c['frame']['run'], camcol=c['frame']['camcol'], field=c['frame']['field']) for c in chosen],
                                     names={k: [v['run'], v['camcol'], v['field']] for k, v in names.items()})), f, indent=1)
        print('wrote', SURVEY_JSON)
    return recs, chosen, names


# ----------------------------------------------------------------------------- build: truth lists and patches of the frames of REPEAT
def load_bundle(run, camcol, field, hugs):
    """Everything the build and the validation need of one frame: the frame and psField (astropy), the photoField r-band entries, the HUGS truth mapped with the asTrans solution
    (`truth`, the lists of this module) and with the header TAN solution (`truth_h`, the mapping of the M2 truth file), the HUGS sources without an F606W (`miss`), the sky and DN images and
    the context of the five rules."""
    P = frame_paths(run, camcol, field)
    fr, ps = S.read_frame(P['frame']), S.read_psfield(P['psfield'])
    pf = read_photofield(run, camcol)
    k = int(np.flatnonzero(np.asarray(pf['FIELD']) == field)[0])
    rec = photofield_record(pf, k)
    w2p = S.world2pix_for(fr, 'astrans')
    truth, truth_h, miss = S.make_truth(hugs, w2p), S.make_truth(hugs, S.world2pix_for(fr, 'header')), S.unmeasured_f606w(hugs, w2p)
    sky = S.sky_image(fr['allsky'], fr['xinterp'], fr['yinterp'])
    dn_full = S.frame_to_dn(fr['nmgy'], fr['cal'], sky)
    return dict(run=run, camcol=camcol, field=field, fr=fr, ps=ps, pf=rec, truth=truth, truth_h=truth_h, miss=miss, sky=sky, dn_full=dn_full, ctx=frame_context(truth, miss, dn_full, rec['sat_level']))


def patch_row(B, x0, y0, arrays):
    """The counts and descriptors of one patch (the quantities of sections 5.1 and 5.2 of the second-cluster spec, as `astro_second.cmd_build` computes them): truth counts in the box, scored
    counts, brightest truth, sources without F606W, DN statistics, FWHM, crowding, peak signal to noise, 20 x 20 block counts, PSF against pydl."""
    t, th = S.truth_for_patch(B['truth'], x0, y0, N_PATCH), S.truth_for_patch(B['truth_h'], x0, y0, N_PATCH)
    miss = B['miss']
    inbox = (miss['px'] >= x0 - 0.5) & (miss['px'] < x0 + N_PATCH - 0.5) & (miss['py'] >= y0 - 0.5) & (miss['py'] < y0 + N_PATCH - 0.5)
    fw = S.psf_fwhm(arrays['psf51'])
    nmgy = float(B['fr']['header']['NMGY'])
    bc = S.block_counts(B['truth'], x0, y0)
    return dict(x0=x0, y0=y0, n22=int((t['mag'] < 22).sum()), density22=float((t['mag'] < 22).sum()) / 1e4, n21=int((t['mag'] < 21).sum()), n23=int((t['mag'] < 23).sum()), n_all=int(len(t['mag'])),
                n22_scored=int(((t['mag'] < 22) & t['inner']).sum()), n22_header_mapping=int((th['mag'] < 22).sum()), brightest=float(t['mag'].min()),
                n_no_f606w_f814w_lt_20=int((inbox & (miss['f814'] < 20)).sum()), n_no_f606w_f814w_lt_22=int(inbox.sum()), dn_max=float(arrays['dn'].max()), dn_median=float(np.median(arrays['dn'])),
                sky_dn_median=float(np.median(arrays['sky_dn'])), fwhm_px=fw, crowding=S.crowding(B['truth'], x0, y0, fw),
                peak_snr_at_22=S.peak_snr(22.0, arrays['psf51'].max(), nmgy, float(np.median(arrays['sky_dn'])), B['ps']['gain'], B['ps']['dark_variance']),
                m2_equal_snr_mag=S.equal_snr_magnitude(S.M2_SNR22, arrays['psf51'].max(), nmgy, float(np.median(arrays['sky_dn'])), B['ps']['gain'], B['ps']['dark_variance']),
                blocks_lt_24=dict(min=int(bc.min()), median=float(np.median(bc)), max=int(bc.max()), below_50=int((bc < 50).sum())),
                psf_vs_pydl=S.psf_vs_pydl(B['ps']['kl'], x0 + (N_PATCH - 1) / 2.0, y0 + (N_PATCH - 1) / 2.0))


def print_patch_row(label, row):
    bl = row['blocks_lt_24']
    print('patch %-44s x0 %4d y0 %4d  F606W<22: %4d (%.4f /px)  <21: %4d  <23: %4d  all: %5d  scored<22: %4d  brightest %.2f  DNmax %7.0f  DNmedian %5.0f  sky %6.1f  FWHM %.2f px  per PSF area %.2f  '
          'blend(2.22 px) %.2f blend(FWHM) %.2f  peak S/N at 22: %.2f (as M2 at 22.0: F606W %.2f)  20x20 blocks F606W<24: min %d median %.0f max %d  no-F606W with F814W<20/<22: %d/%d' % (
              label, row['x0'], row['y0'], row['n22'], row['density22'], row['n21'], row['n23'], row['n_all'], row['n22_scored'], row['brightest'], row['dn_max'], row['dn_median'], row['sky_dn_median'],
              row['fwhm_px'], row['crowding']['per_psf_area'], row['crowding']['blend_m2_fwhm'], row['crowding']['blend_own_fwhm'], row['peak_snr_at_22'], row['m2_equal_snr_mag'], bl['min'], bl['median'],
              bl['max'], row['n_no_f606w_f814w_lt_20'], row['n_no_f606w_f814w_lt_22']))


def cmd_build(args):
    if not REPEAT:
        sys.exit('REPEAT is empty: run `survey`, then fill the registry of this module with the frames it selects')
    hugs = S.read_hugs(_hugs_path())
    table = {}
    for name, ent in REPEAT.items():
        run, camcol, field = ent['frame']
        B = load_bundle(run, camcol, field, hugs)
        print('%s: frame run %d camcol %d field %d (rerun 301); DATE-OBS %s; NMGY %.8f; psField gain %.4f dark variance %.4f; photoField r PSF_WIDTH %.3f arcsec, PSP_STATUS %d, PSF_NSTAR %d, saturation %.0f DN' % (
            name, run, camcol, field, B['fr']['header']['DATE-OBS'], B['fr']['header']['NMGY'], B['ps']['gain'], B['ps']['dark_variance'], B['pf']['psf_width'], B['pf']['psp_status'], B['pf']['psf_nstar'], B['pf']['sat_level']))
        (sx_, pax), (sy_, pay) = frame_orientation(B['fr']['astrans'], *window_centre(*ent['patches'][0]))
        print('  orientation at the first patch: +x (columns) scale %.4f arcsec/px, position angle %+.3f deg; +y (rows) scale %.4f arcsec/px, position angle %+.3f deg (M2 frame: +x north, +y east)' % (sx_, pax, sy_, pay))
        if not args.dry_run:
            os.makedirs(ent['root'], exist_ok=True)
            tp = os.path.join(ent['root'], ent['truth'])
            np.savez(tp, **B['truth'])
            np.savez(tp.replace('_meth1.npz', '_meth1_header.npz'), **B['truth_h'])
            print('  wrote %s (asTrans mapping) and its _header twin (header TAN mapping): %d sources with a valid F606W of %d rows' % (ent['truth'], len(B['truth']['m606']), len(hugs['m606'])))
        rows = []
        for (x0, y0), w in zip(ent['patches'], ent['windows']):
            rules = window_rules(B['ctx'], x0, y0)
            if not rules['eligible']:
                print('  WARNING: window %s with origin (%d, %d) is not eligible: %s' % (w, x0, y0, {k: rules[k] for k in S.RULES}))
            arrays, meta = S.extract_patch(B['fr'], B['ps'], x0, y0, dn_full=B['dn_full'], sky=B['sky'])
            row = patch_row(B, x0, y0, arrays)
            row.update(window=list(w), rules=rules, name=S.patch_name(run, camcol, field, x0, y0))
            if not args.dry_run:
                outdir = os.path.join(ent['root'], row['name'])
                S.write_patch(outdir, arrays, meta, B['fr']['header'])
                files = ('patch_arrays.npz', 'patch_dn.fits', 'meta.json')
                row['sha256'] = dict(patch_arrays_content=S.content_sha256(os.path.join(outdir, files[0])), patch_dn_fits=S.file_sha256(os.path.join(outdir, files[1])), meta_json=S.file_sha256(os.path.join(outdir, files[2])))
                row['bytes'] = {f: os.path.getsize(os.path.join(outdir, f)) for f in files}
            print_patch_row(row['name'], row)
            rows.append(row)
        table[name] = dict(frame=[run, camcol, field], date_obs=str(B['fr']['header']['DATE-OBS']), nmgy=float(B['fr']['header']['NMGY']), photofield=B['pf'], orientation=dict(x=[sx_, pax], y=[sy_, pay]), patches=rows)
        if not args.dry_run and not os.path.exists(tp):
            raise SystemExit('truth file not written')
    if not args.dry_run:
        path = os.path.join(ROOT3, 'repeat_patches.json')
        with open(path, 'w') as f:
            json.dump(_jsonable(table), f, indent=1)
        print('wrote', path)
    return table


# ----------------------------------------------------------------------------- validate: the checks of section 4 of the second-cluster spec for the repeat patches and for the M2 windows
def photometric_scale(dn, gain, dark_variance, psf, tx, ty, tmag, sx, sy, mag_faint=25.0, pad=12):
    """The flux scale (DN per nanomaggy of the model) and chi^2 of the weighted least-squares fit of the PSF-convolved truth model, shifted rigidly by (sx, sy) pixels, with a free plane for
    the sky (the model of `astro_second.rigid_shift` at one shift). For a perfect model, calibration and F606W = r the scale is 1 / NMGY of the frame, so `scale x NMGY` is 1."""
    from scipy.ndimage import shift as nd_shift
    sel = tmag < mag_faint
    canvas = S.render_model(psf, tx[sel], ty[sel], 10.0 ** (-0.4 * (tmag[sel] - S.ZP)), dn.shape, pad)
    mod = nd_shift(canvas, (sy, sx), order=3, mode='nearest')[pad: -pad, pad: -pad]
    sw = 1.0 / np.sqrt(S.dn_variance(dn, gain, dark_variance)).ravel()
    yy, xx = np.mgrid[0:dn.shape[0], 0:dn.shape[1]]
    plane = [np.ones(dn.size), (xx.ravel() - dn.shape[1] / 2.0) / dn.shape[1], (yy.ravel() - dn.shape[0] / 2.0) / dn.shape[0]]
    A = np.vstack([mod.ravel()] + plane).T * sw[:, None]
    b = dn.ravel() * sw
    coef = np.linalg.lstsq(A, b, rcond=None)[0]
    return float(coef[0]), float(((b - A @ coef) ** 2).sum())


MASK_RADIUS = 10.0                                                           # px, the radius around an HST-saturated or lost bright star that the masked scan leaves out


def flagged_sources(B, x0, y0, pad=12):
    """The sources whose HUGS magnitude is not a measurement or that the truth lacks, in the patch's pixel coordinates, within `pad` px of the box: the HST-saturated truth sources (`nf606` = 99)
    and the HUGS sources with F814W < 17 and no F606W (rule (5)). Returns (x, y) arrays."""
    lo_x, hi_x, lo_y, hi_y = x0 - pad - 0.5, x0 + N_PATCH + pad - 0.5, y0 - pad - 0.5, y0 + N_PATCH + pad - 0.5
    tr, ms = B['truth'], B['miss']
    a = (tr['nf606'] == 99) & (tr['px'] >= lo_x) & (tr['px'] < hi_x) & (tr['py'] >= lo_y) & (tr['py'] < hi_y)
    b = (ms['f814'] < S.MISSING_LIMIT) & (ms['px'] >= lo_x) & (ms['px'] < hi_x) & (ms['py'] >= lo_y) & (ms['py'] < hi_y)
    return np.r_[tr['px'][a], ms['px'][b]] - x0, np.r_[tr['py'][a], ms['py'][b]] - y0


def masked_scan(dn, gain, dark_variance, psf, truth, x0, y0, flagged, mask_radius=MASK_RADIUS, mag_faint=25.0, half_range=0.4, step=0.05, pad=12):
    """`astro_second.rigid_shift` with the sources whose HUGS magnitude is not a measurement left out of the model and the pixels within `mask_radius` of them (and of the lost bright sources,
    `flagged`: their patch coordinates) given zero weight. A star the HST saturated (F606W of 12 to 14.5 in the truth) can be wrong by magnitudes and sets the flux scale of the whole fit
    (the plain scan then returns a scale of 0.1 to 0.6 and a shift of 0.5 px). Returns dict(dx, dy, chi2_best, chi2_zero, scale, n_pixels, n_model): the shift to add to the truth positions,
    the chi^2 at the best and at zero shift, the flux scale (DN per nanomaggy of the model) at the best shift, the pixels used and the model sources."""
    from scipy.ndimage import shift as nd_shift
    x, y, m = truth['px'] - x0, truth['py'] - y0, truth['m606']
    sel = (x >= -pad - 0.5) & (x < N_PATCH + pad - 0.5) & (y >= -pad - 0.5) & (y < N_PATCH + pad - 0.5) & (m < mag_faint) & (truth['nf606'] != 99)
    canvas = S.render_model(psf, x[sel], y[sel], 10.0 ** (-0.4 * (m[sel] - S.ZP)), dn.shape, pad)
    yy, xx = np.mgrid[0:dn.shape[0], 0:dn.shape[1]]
    bad = np.zeros(dn.shape, bool)
    for xb, yb in zip(*flagged):
        bad |= (xx - xb) ** 2 + (yy - yb) ** 2 <= mask_radius ** 2
    sw = (1.0 / np.sqrt(S.dn_variance(dn, gain, dark_variance))) * (~bad)
    sw = sw.ravel()
    plane = [np.ones(dn.size), (xx.ravel() - dn.shape[1] / 2.0) / dn.shape[1], (yy.ravel() - dn.shape[0] / 2.0) / dn.shape[0]]
    grid = np.arange(-half_range, half_range + 1e-9, step)
    chi, scl = np.zeros((len(grid), len(grid))), np.zeros((len(grid), len(grid)))
    for j, sy in enumerate(grid):
        for i, sx in enumerate(grid):
            mod = nd_shift(canvas, (sy, sx), order=3, mode='nearest')[pad: -pad, pad: -pad]
            A = np.vstack([mod.ravel()] + plane).T * sw[:, None]
            b = dn.ravel() * sw
            coef = np.linalg.lstsq(A, b, rcond=None)[0]
            chi[j, i], scl[j, i] = float(((b - A @ coef) ** 2).sum()), float(coef[0])
    j0, i0 = np.unravel_index(np.argmin(chi), chi.shape)
    j0, i0 = min(max(j0, 1), len(grid) - 2), min(max(i0, 1), len(grid) - 2)

    def vertex(c_m, c_0, c_p):
        den = c_m - 2 * c_0 + c_p
        return 0.0 if den <= 0 else 0.5 * (c_m - c_p) / den
    return dict(dx=float(grid[i0] + step * vertex(chi[j0, i0 - 1], chi[j0, i0], chi[j0, i0 + 1])), dy=float(grid[j0] + step * vertex(chi[j0 - 1, i0], chi[j0, i0], chi[j0 + 1, i0])), chi2_best=float(chi.min()),
                chi2_zero=float(chi[len(grid) // 2, len(grid) // 2]), scale=float(scl[np.unravel_index(np.argmin(chi), chi.shape)]), n_pixels=int((sw > 0).sum()), n_model=int(sel.sum()))


def check_patch(B, x0, y0, arrays, pad=12):
    """The validation numbers of one patch: the asTrans against the header mapping on the box; the rigid-shift scan of the PSF-convolved truth against the patch for both mappings (shift to add to
    the truth, chi^2 at the best and at zero shift) and the flux scale at the best asTrans shift (and `phot_ratio` = scale x NMGY of the frame, 1 for a perfect calibration); bright truth stars on bright
    pixels; the PSF against pydl; the saturation figures."""
    out = {}
    t, th = S.truth_for_patch(B['truth'], x0, y0, N_PATCH), S.truth_for_patch(B['truth_h'], x0, y0, N_PATCH)
    m = (B['truth']['px'] >= x0 - 0.5) & (B['truth']['px'] < x0 + N_PATCH - 0.5) & (B['truth']['py'] >= y0 - 0.5) & (B['truth']['py'] < y0 + N_PATCH - 0.5)
    dxy = np.c_[B['truth']['px'][m] - B['truth_h']['px'][m], B['truth']['py'][m] - B['truth_h']['py'][m]]
    dd = np.hypot(dxy[:, 0], dxy[:, 1])
    out['asTrans_minus_header_px'] = dict(median=float(np.median(dd)), p90=float(np.percentile(dd, 90)), max=float(dd.max()), signed_median_dx=float(np.median(dxy[:, 0])), signed_median_dy=float(np.median(dxy[:, 1])))
    shifts = {}
    for label, tt in (('asTrans', B['truth']), ('header', B['truth_h'])):
        tp = S.truth_for_patch(tt, x0 - pad, y0 - pad, N_PATCH + 2 * pad)
        sx, sy, chi_min, chi_0 = S.rigid_shift(arrays['dn'], B['ps']['gain'], B['ps']['dark_variance'], arrays['psf51'], tp['x'] - pad, tp['y'] - pad, tp['mag'])
        shifts[label] = dict(dx=sx, dy=sy, chi2_best=chi_min, chi2_zero=chi_0)
        if label == 'asTrans':
            scale, chi_s = photometric_scale(arrays['dn'], B['ps']['gain'], B['ps']['dark_variance'], arrays['psf51'], tp['x'] - pad, tp['y'] - pad, tp['mag'], sx, sy)
            nm = float(B['fr']['header']['NMGY'])
            out['scale'] = dict(scale=scale, nmgy=nm, phot_ratio=scale * nm, chi2=chi_s, chi2_reduced=chi_s / (N_PATCH * N_PATCH - 4))
    out['image_shift_px'] = shifts
    fl = flagged_sources(B, x0, y0, pad)
    masked = {}
    for label, tt in (('asTrans', B['truth']), ('header', B['truth_h'])):
        r = masked_scan(arrays['dn'], B['ps']['gain'], B['ps']['dark_variance'], arrays['psf51'], tt, x0, y0, fl, pad=pad)
        r['phot_ratio'] = r['scale'] * float(B['fr']['header']['NMGY'])
        r['chi2_reduced'] = r['chi2_best'] / (r['n_pixels'] - 4)
        masked[label] = r
    out['masked'] = dict(n_flagged=len(fl[0]), n_masked_pixels=N_PATCH * N_PATCH - masked['asTrans']['n_pixels'], **masked)
    br = (t['mag'] < 18) & t['inner']
    ix, iy = np.rint(t['x'][br]).astype(int), np.rint(t['y'][br]).astype(int)
    ratio = np.array([arrays['dn'][max(yy - 1, 0): yy + 2, max(xx - 1, 0): xx + 2].max() / np.median(arrays['dn']) for yy, xx in zip(iy, ix)])
    out['bright_on_bright'] = dict(n=int(br.sum()), median_ratio=float(np.median(ratio)) if br.any() else None, fraction_ratio_ge_2=float(np.mean(ratio >= 2.0)) if br.any() else None)
    sat = B['pf']['sat_level']
    out['saturation'] = dict(level=sat, dn_max_patch=float(arrays['dn'].max()), dn_max_box_plus_10px=B['ctx']['dn'][max(y0 - S.DN_MARGIN, 0): y0 + N_PATCH + S.DN_MARGIN, max(x0 - S.DN_MARGIN, 0): x0 + N_PATCH + S.DN_MARGIN].max().item(),
                             fraction_of_level=float(arrays['dn'].max() / sat), pixels_above_0p8=int((arrays['dn'] > S.SAT_FRACTION * sat).sum()))
    out['psf_vs_pydl'] = S.psf_vs_pydl(B['ps']['kl'], x0 + (N_PATCH - 1) / 2.0, y0 + (N_PATCH - 1) / 2.0)
    return out


def print_check(label, c):
    a, s, b = c['asTrans_minus_header_px'], c['image_shift_px'], c['bright_on_bright']
    print('  %-44s asTrans - header on the box: median %.4f p90 %.4f max %.4f px (signed median %+.3f, %+.3f) ; image shift asTrans (%+.3f, %+.3f) chi2 %.0f -> %.0f ; header (%+.3f, %+.3f) chi2 %.0f -> %.0f' % (
        label, a['median'], a['p90'], a['max'], a['signed_median_dx'], a['signed_median_dy'], s['asTrans']['dx'], s['asTrans']['dy'], s['asTrans']['chi2_zero'], s['asTrans']['chi2_best'], s['header']['dx'], s['header']['dy'],
        s['header']['chi2_zero'], s['header']['chi2_best']))
    mk = c['masked']
    print('  %-44s masked scan (HST-saturated and lost bright stars left out, %d of them, %d pixels masked): asTrans shift (%+.3f, %+.3f) chi2 %.0f -> %.0f, scale x NMGY %.4f, reduced chi2 %.2f ; header (%+.3f, %+.3f) chi2 %.0f -> %.0f, scale x NMGY %.4f' % (
        '', mk['n_flagged'], mk['n_masked_pixels'], mk['asTrans']['dx'], mk['asTrans']['dy'], mk['asTrans']['chi2_zero'], mk['asTrans']['chi2_best'], mk['asTrans']['phot_ratio'], mk['asTrans']['chi2_reduced'],
        mk['header']['dx'], mk['header']['dy'], mk['header']['chi2_zero'], mk['header']['chi2_best'], mk['header']['phot_ratio']))
    q, z = c['scale'], c['saturation']
    print('  %-44s flux scale of the truth model %.3f DN/nMgy (NMGY %.8f, scale x NMGY = %.4f), reduced chi2 %.2f ; bright(<18) 3x3 max >= 2 x median: %s of %d (median ratio %s) ; PSF against pydl: max |diff| %s ; DNmax %.0f = %.3f of the saturation level %.0f, '
          'box+10 px max %.0f, pixels above 0.8 x level: %d' % ('', q['scale'], q['nmgy'], q['phot_ratio'], q['chi2_reduced'], ('%.2f' % b['fraction_ratio_ge_2']) if b['n'] else 'n/a', b['n'], ('%.1f' % b['median_ratio']) if b['n'] else 'n/a',
                                                                  ('%.2e' % c['psf_vs_pydl']) if c['psf_vs_pydl'] is not None else 'n/a (pydl absent)', z['dn_max_patch'], z['fraction_of_level'], z['level'], z['dn_max_box_plus_10px'], z['pixels_above_0p8']))


def cmd_validate(args):
    if not REPEAT:
        sys.exit('REPEAT is empty: run `survey`, then fill the registry of this module with the frames it selects')
    hugs = S.read_hugs(_hugs_path())
    result = {}
    # the M2 frame: the same code path, the four windows, against the patches of data/astro/fetch.py and the M2 truth file
    B0 = load_bundle(*M2_FRAME, hugs)
    ref = np.load(os.path.join(ROOT, 'truth_hugs_meth1.npz'))
    print('[1] M2 frame (run %d camcol %d field %d): the HUGS truth mapped with asTrans (the lists of this module) against truth_hugs_meth1.npz (header TAN mapping), per window' % M2_FRAME)
    for w in M2_WINDOWS:
        m = (B0['truth_h']['px'] >= w[0] - 0.5) & (B0['truth_h']['px'] < w[0] + N_PATCH - 0.5) & (B0['truth_h']['py'] >= w[1] - 0.5) & (B0['truth_h']['py'] < w[1] + N_PATCH - 0.5)
        d_hdr = np.hypot(B0['truth_h']['px'][m] - ref['px'][m], B0['truth_h']['py'][m] - ref['py'][m])
        d_as = np.hypot(B0['truth']['px'][m] - B0['truth_h']['px'][m], B0['truth']['py'][m] - B0['truth_h']['py'][m])
        print('    window (%4d, %4d): %d sources; header mapping against the file: max %.2e px; asTrans against header: median %.4f p90 %.4f max %.4f px' % (
            w[0], w[1], m.sum(), d_hdr.max(), np.median(d_as), np.percentile(d_as, 90), d_as.max()))
    (sx0, pax0), (sy0, pay0) = frame_orientation(B0['fr']['astrans'], *window_centre(*M2_WINDOWS[0]))
    print('    M2 frame orientation at the benchmark window: +x scale %.4f arcsec/px, position angle %+.3f deg; +y scale %.4f arcsec/px, position angle %+.3f deg' % (sx0, pax0, sy0, pay0))
    print('[2] units: `extract_patch` on the M2 frame and psField against the patches of data/astro/fetch.py, and the window rules in the M2 frame')
    m2 = {}
    for w in M2_WINDOWS:
        arrays, meta = S.extract_patch(B0['fr'], B0['ps'], w[0], w[1], dn_full=B0['dn_full'], sky=B0['sky'])
        old = np.load(os.path.join(ROOT, 'patch_r_002583_2_0136_x%d_y%d_100x100' % w, 'patch_arrays.npz'))
        same = all(np.array_equal(np.asarray(arrays[k]), old[k]) for k in old.files)
        meta_same = meta == json.load(open(os.path.join(ROOT, 'patch_r_002583_2_0136_x%d_y%d_100x100' % w, 'meta.json')))
        rules = window_rules(B0['ctx'], w[0], w[1])
        print('    window (%4d, %4d): every array of patch_arrays.npz identical to the fetch.py patch: %s ; meta.json identical: %s ; rules %s, DNmax within 10 px %.0f (0.8 x level %.0f)' % (
            w[0], w[1], same, meta_same, {k: rules[k] for k in S.RULES}, rules['dn_max'], S.SAT_FRACTION * B0['pf']['sat_level']))
        m2[w] = (arrays, rules)
    print('[3] image-based checks (rigid-shift scan of the PSF-convolved truth, flux scale, bright stars, PSF against pydl, saturation); the M2 windows first')
    scales, bundles = {}, {}
    result['m2'] = {}
    for w in M2_WINDOWS:
        arrays, rules = m2[w]
        c = check_patch(B0, w[0], w[1], arrays)
        c['rules'] = rules
        c['fwhm_px'] = S.psf_fwhm(arrays['psf51'])
        result['m2']['%d,%d' % w] = c
        print_check('M2 %s' % (w,), c)
        scales[('m2', w)] = c['masked']['asTrans']['phot_ratio']
    for name, ent in REPEAT.items():
        run, camcol, field = ent['frame']
        B = bundles[name] = load_bundle(run, camcol, field, hugs)
        print('%s: run %d camcol %d field %d, NMGY %.8f (%+.4f mag against the M2 frame), seeing PSF_WIDTH %.3f arcsec' % (
            name, run, camcol, field, B['fr']['header']['NMGY'], -2.5 * np.log10(B['fr']['header']['NMGY'] / B0['fr']['header']['NMGY']), B['pf']['psf_width']))
        result[name] = {}
        for (x0, y0), w in zip(ent['patches'], ent['windows']):
            z = np.load(patch_path(x0, y0, name))
            arrays = {k: z[k] for k in z.files}
            mem, _ = S.extract_patch(B['fr'], B['ps'], x0, y0, dn_full=B['dn_full'], sky=B['sky'])
            same = all(np.array_equal(np.asarray(mem[k]), arrays[k]) for k in arrays)
            c = check_patch(B, x0, y0, arrays)
            c['rules'] = window_rules(B['ctx'], x0, y0)
            c['fwhm_px'] = S.psf_fwhm(arrays['psf51'])
            c['rebuild_identical'] = bool(same)
            result[name]['%d,%d' % (x0, y0)] = c
            print_check('%s (%d, %d) = window %s' % (name, x0, y0, tuple(w)), c)
            print('  %-44s rebuilt in memory from the frame: every array identical to the file: %s ; FWHM %.3f px' % ('', same, c['fwhm_px']))
            scales[(name, tuple(w))] = c['masked']['asTrans']['phot_ratio']
    print('[4] the same stars: HUGS sources of F606W < 22 in the box of a repeat patch that are also in the box of the M2 window (the truth lists of all frames are row for row the same sources)')
    for name, ent in REPEAT.items():
        B = bundles[name]
        assert np.array_equal(B['truth']['m606'], ref['m606'])
        for (x0, y0), w in zip(ent['patches'], ent['windows']):
            inb = lambda tr_, a, b: (tr_['px'] >= a - 0.5) & (tr_['px'] < a + N_PATCH - 0.5) & (tr_['py'] >= b - 0.5) & (tr_['py'] < b + N_PATCH - 0.5) & (tr_['m606'] < 22)
            mine, theirs = inb(B['truth'], x0, y0), inb(ref, w[0], w[1])
            (_, pax), (_, pay) = frame_orientation(B['fr']['astrans'], *window_centre(x0, y0))
            print('    %s (%d, %d) = window %s: %d sources in the repeat box, %d in the M2 box, %d in both (%.1f %% of the repeat box, %.1f %% of the M2 box); frame axes rotated by %+.2f deg against the M2 frame' % (
                name, x0, y0, tuple(w), mine.sum(), theirs.sum(), (mine & theirs).sum(), 100.0 * (mine & theirs).sum() / mine.sum(), 100.0 * (mine & theirs).sum() / theirs.sum(), pax - pax0))
            result[name]['%d,%d' % (x0, y0)]['same_stars'] = dict(n_repeat=int(mine.sum()), n_m2=int(theirs.sum()), n_both=int((mine & theirs).sum()), rotation_deg=pax - pax0)
    print('[5] units: flux scale of the truth model in the masked scan (scale x NMGY, 1 for a perfect calibration) of the repeat patches against the same window in the M2 frame')
    for name in REPEAT:
        r = [scales[(name, tuple(w))] / scales[('m2', tuple(w))] for w in REPEAT[name]['windows']]
        print('    %s: ratios %s ; median %.4f' % (name, ' '.join('%.4f' % v for v in r), np.median(r)))
    path = os.path.join(ROOT3, 'repeat_validation.json')
    with open(path, 'w') as f:
        json.dump(_jsonable(result), f, indent=1)
    print('wrote', path)
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('survey')
    p.add_argument('--stage', choices=['runs', 'geometry', 'frames', 'all'], default='all')
    p.add_argument('--no-json', action='store_true', help='do not write the survey JSON')
    p = sub.add_parser('build')
    p.add_argument('--dry-run', action='store_true', help='print the counts, write nothing')
    sub.add_parser('validate')
    args = ap.parse_args()
    {'survey': cmd_survey, 'build': cmd_build, 'validate': cmd_validate}[args.cmd](args)


if __name__ == '__main__':
    main()
