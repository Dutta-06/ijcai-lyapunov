"""
The sliding inverse scale space on a crowded stellar field: the M2 patch of SDSS against HST (WP5)
==================================================================================================

The benchmark of Portillo et al. (2017) as defined by `experiments/astro_common.py`: a 100 x 100 pixel r-band patch of M2, about 0.1 sources per pixel brighter than 22 mag, truth the HST
(HUGS) catalogue of the same stars; detections are matched at 0.5 px and 0.5 mag; the figures are completeness (matched truth sources / truth sources, by magnitude) and the fraction of false
positives. `src/core/imaging.py` with the SDSS PSF of the patch (the psField reconstruction, 51 x 51, FWHM 2.2 px) represented as a sum of three Gaussians (`fit_psf_parts`), Poisson-plus-dark
noise (variance DN/gain + dark variance from the smoothed data), a border-consistent high-pass of width sigma_b for the sky, and the inverse-scale-space stopping time 1/(kappa x noise scale): the
entrants of time-averaged significance at least kappa that are also, once fitted, of own amplitude significance at least kappa (the instantaneous significance of the entrant, the rule of the SMLM frames, is fooled in a crowded field by candidates that built up their dual before a neighbour was fitted).

Protocol: tuned on two patches of the same frame inside the HST coverage ((110, 630) and (210, 630), non-overlapping with every test patch), by the F-measure of completeness and purity at
mag < 22 (the benchmark's depth); then run once on each test patch: the benchmark's patch (310, 630) and three more at other crowding levels, (310, 730), (360, 1040), (100, 1100) (truth
sources brighter than 22 mag in the scored area: 969, 2057, 923, 1446). Referee round IS1a (pre-test) showed that the first choice of test patches, (410, 630) and (410, 730), lay mostly
outside the HUGS footprint (28 of 100 blocks without a single truth source brighter than 22 mag while the pixels hold stars: every detection there is a false positive by construction); they
are replaced by patches whose 25 blocks all hold at least 50 truth sources of m < 24, chosen from truth geometry alone (`experiments/astro_patch_survey.py`),
before any test-patch detection was made. Also from that round: the kappa grid is extended to 2 (the optimum of the first grid lay at its edge), the path runs to the kappa = 2 time,
the baseline family is extended and scored on the development patches by part DB, and the freeze is mechanical (`--freeze` writes results/astro/freeze.json from the development result;
part T refuses to run if a hashed file differs or any test result exists). The baselines (photutils DAOStarFinder + IterativePSFPhotometry with 12 variants, crowdsource fit_im with 48 post-filters; `exp_astro_baselines.py`) run at their default settings and,
as a second row, with the variant tuned on the same two development patches by the same F-measure (the same tuning effort and objective as the path's). The path also runs with an a priori
default (`DEFAULT`): fixed in the code before the development grid ran and not chosen by it (a first smoke run of the path used it on development patch (110, 630)); its lengths are in pixels. Nothing is tuned on a test patch. Intervals: bootstrap over the 20 x 20 pixel blocks of the patch (25 per patch; seed 0), paired for the differences.

    python experiments/exp_astro_m2.py --only D|DB|T|V [--workers 3]    (parts merge into results/astro/m2_raw.npz; D the path on the development patches, DB the baselines on them; --freeze after D and DB; needs $ASTRO_DATA with the patches, the truth list and, for the baselines, results/astro/baseline_catalogs.npz)
    python experiments/exp_astro_m2.py --from-raw                       (regenerates results/astro/m2.md from the raw alone)
"""
import argparse
import json
import os
import sys
import time

os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np

HERE = os.path.abspath(os.path.dirname(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
sys.path.insert(0, REPO)
import astro_common as A                  # noqa: E402
from src.core import imaging as IM        # noqa: E402

ROOT = os.environ.get('ASTRO_DATA', '/home/user/astro_data')
OUT = os.environ.get('ASTRO_OUT', os.path.join(REPO, 'results', 'astro'))
RAW, MD = os.path.join(OUT, 'm2_raw.npz'), os.path.join(OUT, 'm2.md')
BASE = os.path.join(REPO, 'results', 'astro', 'baseline_catalogs.npz')
ROOT2 = os.environ.get('ASTRO_DATA_2', '/home/user/astro_data_2')
BASE2 = os.path.join(REPO, 'results', 'astro', 'baseline_catalogs_second.npz')
# the second clusters (results/astro/second_cluster_spec.md): test sets with the configuration and the baselines' tuning frozen on the M2 development patches, nothing re-tuned
SECOND = {'m53': dict(title='M53 (NGC 5024), SDSS run 5390 camcol 6 field 162', root=ROOT2, prefix='patch_r_005390_6_0162', truth='truth_hugs_ngc5024_meth1.npz',
                      patches=[(760, 590), (640, 680), (400, 660), (660, 420)]),
          'm92': dict(title='M92 (NGC 6341), SDSS run 4682 camcol 6 field 208', root=ROOT2, prefix='patch_r_004682_6_0208', truth='truth_hugs_ngc6341_meth1.npz',
                      patches=[(1100, 220), (1430, 550)])}
TESTS, DEV = [(310, 630), (310, 730), (360, 1040), (100, 1100)], [(110, 630), (210, 630)]      # test patches (the benchmark's first), development patches
DEFAULT = (5.0, 1.2, 5.0)                                                                       # a priori (sigma_b px, r_ex px, kappa): the middle high-pass, the 1.2 px resolution limit, the 5-sigma rule
BOOT = 2000
KAPPAS = (2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0)             # extended from 3 to 2 after referee round IS1a (the optimum of the first grid sat at its lower edge)
HPS = (3.0, 5.0, 8.0)
REXS = (0.8, 1.2)
LAT_H, MARGIN, FIT_RADIUS, MAX_ENTRIES = 0.5, 2.0, 4.0, 2500
import exp_astro_baselines as AB       # noqa: E402  (numpy only at import; astropy, photutils and crowdsource are imported inside its runners)
PHOT_TAGS = tuple('photutils' if (lb, ns, it) == (False, 5.0, 3) else 'photutils_%s_t%g_i%d' % ('local' if lb else 'global', ns, it)
                  for lb in (False, True) for ns in (3.0, 5.0) for it in (3, 6, 10)) + tuple(AB.EXT_TAGS)       # photutils variants (exp_astro_baselines.PHOT_VARIANTS and PHOT_EXT); 'photutils' is the default
FREEZE = os.path.join(OUT, 'freeze.json')
FREEZE_FILES = ('src/core/imaging.py', 'experiments/astro_common.py', 'experiments/exp_astro_baselines.py', 'results/astro/baseline_catalogs.npz', 'results/astro/baseline_catalogs_second.npz')      # hashed whole; this script by `frozen_code_hash`
CS_FILTERS = [(snr, ff, qf, sp) for snr in (0.0, 3.0, 5.0, 8.0) for ff in (0.0, 0.5, 0.8) for qf in (0.0, 0.5) for sp in (np.inf, 0.02)]      # crowdsource post-filters, tuned on the development patches


def cs_name(f):
    return 'crowdsource|snr%g|ff%g|qf%g|sp%g' % f


def cs_dets(b, patch, f):
    """The crowdsource catalogue of a patch with the post-filter f = (min flux/dflux, min fracflux, min qf, max spread_model); the default (no filter) is (0, 0, 0, inf)."""
    key = 'x%d_y%d|crowdsource|' % patch
    snr, ff, qf, sp = f
    m = (b[key + 'flux'] > 0) & (b[key + 'flux'] >= snr * b[key + 'dflux']) & (b[key + 'fracflux'] >= ff) & (b[key + 'qf'] >= qf) & (np.abs(b[key + 'spread']) <= sp)
    return b[key + 'x'][m], b[key + 'y'][m], b[key + 'flux'][m]


def baseline_scores(patch, t, base=None):
    """Scores of every baseline variant on a patch, {name: score dict} (none if the baseline catalogues are absent); `base`: the catalogue file (default: the M2 patches')."""
    out = {}
    base = base or BASE
    if not os.path.exists(base):
        return out
    b = np.load(base)
    for tag in PHOT_TAGS:
        key = 'x%d_y%d|%s|' % (patch[0], patch[1], tag)
        if key + 'x' in b.files:
            out[tag] = score_dets((b[key + 'x'], b[key + 'y'], b[key + 'flux']), t)
    if 'x%d_y%d|crowdsource|x' % patch in b.files:
        for f in CS_FILTERS:
            out[cs_name(f)] = score_dets(cs_dets(b, patch, f), t)
    for name in ('daophot', 'daophot_asis'):               # the shipped DAOPHOT catalogue with its y shift of -1 px (fitted to the HUGS truth of that patch) and as shipped
        if 'x%d_y%d|%s|x' % (patch + (name,)) in b.files:
            k = 'x%d_y%d|%s|' % (patch + (name,))
            out[name] = score_dets((b[k + 'x'], b[k + 'y'], b[k + 'flux']), t)
    return out


def patch_path(x0, y0, ds='m2'):
    if ds == 'm2':
        return os.path.join(ROOT, 'patch_r_002583_2_0136_x%d_y%d_100x100' % (x0, y0), 'patch_arrays.npz')
    return os.path.join(SECOND[ds]['root'], '%s_x%d_y%d_100x100' % (SECOND[ds]['prefix'], x0, y0), 'patch_arrays.npz')


def load_patch(x0, y0, ds='m2'):
    z = np.load(patch_path(x0, y0, ds))
    return dict(dn=z['dn'].astype(float), psf=z['psf51'].astype(float), gain=float(z['gain']), dark=float(z['dark_variance']))


def run_path(patch, hp, rex, kappas=KAPPAS, snr_min=None, order='dual', stop='time', cleanup=True):
    """The path on one patch; returns ({kappa: (x, y, flux)}, seconds, number of entries). order: 'dual' (the accumulated dual enters first: the sliding inverse scale space) or
    'greedy' (the instantaneous correlation: matching pursuit with sliding); stop: 'time' (the time rule, tau <= 1 / (kappa x noise scale)) or 'snr' (the entrant's instantaneous
    significance, the rule of the SMLM frames); cleanup: discard the units whose own amplitude significance is below kappa; hp 0: no high-pass; rex 0: no exclusion zone. With stop='time' the
    result also has key 0.0: the block of the smallest BIC along the same path (nothing to tune), cleaned up at sqrt(3 ln N). Returns (kappa -> (x, y, flux), seconds, number of entries,
    info), info = dict(bic_index, n_blocks): the BIC block's index along the path and the number of blocks (the BIC minimum is interior iff bic_index < n_blocks - 1)."""
    sigma, parts, rms = IM.fit_psf_parts(patch['psf'], n_parts=3, sigma0=1.0)
    ny, nx = patch['dn'].shape
    P = IM.PixelPSF(ny, nx, sigma, parts=parts, hp=hp or None)
    L = IM.Lattice(P, h=LAT_H, margin=MARGIN)
    from scipy.ndimage import gaussian_filter
    dn = patch['dn']
    expected = np.maximum(gaussian_filter(dn, 1.5, mode='reflect'), np.percentile(dn, 1))
    sw = 1.0 / np.sqrt(expected / patch['gain'] + patch['dark'])
    y = P.highpass(dn - np.median(dn))
    reach = P.reach()
    scale = IM.noise_scale(L, sw ** 2)
    t0 = time.time()
    R = IM.sliding_path(P, y, sw, L, max_entries=MAX_ENTRIES, tol=1e-6, passes=1, rtrunc=10.0 / reach, r_ex=rex, fit_radius=FIT_RADIUS, order=order,
                        tau_max=(1.0 / (min(kappas) * scale)) if stop == 'time' else None, snr_min=(min(kappas) if stop == 'snr' else snr_min))
    secs = time.time() - t0
    out = {}
    for k in kappas:
        i, K = IM.stop_time(R['hist'], R['taus'], 1.0 / (k * scale)) if stop == 'time' else IM.stop_snr(R['hist'], R['taus'], R['snr'], k)
        h = R['hist'][i]
        keep = (IM.unit_significance(P, sw, h['c'], h['w']) >= k) if cleanup else np.ones(len(h['w']), bool)      # the unit's own amplitude significance: removes ghosts and faint blends
        out[k] = (h['c'][keep, 0].copy(), h['c'][keep, 1].copy(), h['w'][keep].copy())
    if stop == 'time':                                  # key 0: the BIC block of the same path, the stop with nothing to tune (3 parameters per unit: position and flux)
        i, K = IM.stop_bic(R['hist'], ny * nx, n_par=3)
        h = R['hist'][i]
        kb = float(np.sqrt(3.0 * np.log(ny * nx)))      # the significance a single unit needs to pay for its parameters: the same penalty for the clean-up
        keep = (IM.unit_significance(P, sw, h['c'], h['w']) >= kb) if cleanup else np.ones(len(h['w']), bool)
        out[0.0] = (h['c'][keep, 0].copy(), h['c'][keep, 1].copy(), h['w'][keep].copy())
        return out, secs, len(R['taus']), dict(bic_index=int(i), n_blocks=len(R['hist']))
    return out, secs, len(R['taus']), dict(bic_index=-1, n_blocks=len(R['hist']))


def score_dets(det, truth_patch):
    x, y, f = det
    return A.score(np.c_[x, y], f, truth_patch)


def _pack(name, s):
    return {name + '|' + k: np.asarray(s[k]) for k in ('n_true', 'n_matched', 'n_det', 'n_det_matched', 'completeness', 'fp_fraction', 'offset', 'blk22')} | {
        name + '|tot': np.array([s['true22'], s['matched22'], s['det22'], s['det22_matched'], s['true_all'], s['matched_all'], s['det_all']])}


def save_part(path, part, arrays):
    old = {}
    if os.path.exists(path):
        with np.load(path, allow_pickle=False) as z:
            old = {k: z[k] for k in z.files if not k.startswith(part + '|')}
    old.update({part + '|' + k: v for k, v in arrays.items()})
    os.makedirs(os.path.dirname(path), exist_ok=True)
    np.savez_compressed(path, **old)


def load_raw(path=RAW):
    with np.load(path, allow_pickle=False) as z:
        return {k: z[k] for k in z.files}


def truth_arrays(ds='m2'):
    z = np.load(os.path.join(ROOT if ds == 'm2' else SECOND[ds]['root'], 'truth_hugs_meth1.npz' if ds == 'm2' else SECOND[ds]['truth']))
    return dict(px=z['px'], py=z['py'], m606=z['m606'])


def f_measure(tot):
    """F of completeness and purity at mag < 22 from the totals [true22, matched22, det22, det22_matched, ...]."""
    comp = tot[1] / max(tot[0], 1)
    pur = tot[3] / max(tot[2], 1)
    return 2 * comp * pur / max(comp + pur, 1e-12), comp, pur


def part_D(workers=3, out_path=RAW):
    """Development patches: every (sigma_b, r_ex) of the path with all thresholds (the baselines are scored by part DB)."""
    from multiprocessing import Pool
    truth = truth_arrays()
    jobs = [(p, hp, rex) for p in DEV for hp in HPS for rex in REXS]
    with Pool(workers) as pool:
        res = pool.map(_dev_job, [(j, truth) for j in jobs], chunksize=1)
    arrays = {}
    for (p, hp, rex), (scores, secs, n, info) in zip(jobs, res):
        for k, s_ in scores.items():
            arrays.update(_pack('x%d_y%d|hp%g|rex%g|k%g' % (p[0], p[1], hp, rex, k), s_))
        arrays['x%d_y%d|hp%g|rex%g|secs' % (p[0], p[1], hp, rex)] = np.array([secs, n])
        arrays['x%d_y%d|hp%g|rex%g|bic' % (p[0], p[1], hp, rex)] = np.array([info['bic_index'], info['n_blocks'], int(n >= MAX_ENTRIES)])
        print('D patch %s hp %g rex %g: %d entries, %.0f s, BIC block %d of %d' % (p, hp, rex, n, secs, info['bic_index'], info['n_blocks']), flush=True)
    save_part(out_path, 'D', arrays)


def part_DB(out_path=RAW):
    """The baselines on the development patches: every photutils variant (PHOT_TAGS), every crowdsource post-filter (CS_FILTERS), from results/astro/baseline_catalogs.npz. No truth enters the
    catalogues; the scores select the 'tuned' variant of each family (`choose_baselines`)."""
    truth = truth_arrays()
    arrays = {}
    for p in DEV:
        for name, s_ in baseline_scores(p, A.truth_for_patch(truth, p[0], p[1])).items():
            arrays.update(_pack('x%d_y%d|%s' % (p[0], p[1], name), s_))
    save_part(out_path, 'DB', arrays)


def _dev_job(args):
    (p, hp, rex), truth = args
    patch = load_patch(*p)
    t = A.truth_for_patch(truth, p[0], p[1])
    dets, secs, n, info = run_path(patch, hp, rex)
    return {k: score_dets(d, t) for k, d in dets.items()}, secs, n, info


def choose(raw):
    """(sigma_b, r_ex, kappa) of the largest F-measure at mag < 22, the counts of the two development patches added."""
    best = None
    for hp in HPS:
        for rex in REXS:
            for k in KAPPAS:
                tot = sum(raw['D|x%d_y%d|hp%g|rex%g|k%g|tot' % (p[0], p[1], hp, rex, k)] for p in DEV)
                f = f_measure(tot)[0]
                if best is None or (round(f, 6), k) > (round(best[0], 6), best[3]):
                    best = (f, hp, rex, k)
    return best[1], best[2], best[3], best[0]


def _dev_F(raw, name):
    return f_measure(sum(raw['DB|x%d_y%d|%s|tot' % (p[0], p[1], name)] for p in DEV))[0]


def choose_baselines(raw):
    """The development-tuned variant of each baseline: the photutils threshold and the crowdsource post-filter of the largest F at mag < 22 (counts of the two development patches added)."""
    out = {}
    tags = [t for t in PHOT_TAGS if all('DB|x%d_y%d|%s|tot' % (p[0], p[1], t) in raw for p in DEV)]
    if tags:
        out['photutils'] = max(tags, key=lambda t: (round(_dev_F(raw, t), 6), t))
    names = [cs_name(f) for f in CS_FILTERS if all('DB|x%d_y%d|%s|tot' % (p[0], p[1], cs_name(f)) in raw for p in DEV)]
    if names:
        out['crowdsource'] = max(names, key=lambda n: (round(_dev_F(raw, n), 6), n))
    return out


# 'greedy' with the time stop was dropped after its first run: the time rule needs the dual's entry times, the greedy path's never reach the stop (the path ran to the entry cap, 2500, with 21 detections kept),
# so that row measured the stop rule, not the order; the order is compared with the instantaneous-significance stop that both paths share ('greedy_snr' against 'snr').
VARIANTS = {'greedy_snr': dict(order='greedy', stop='snr'), 'nohp': dict(hp=0.0), 'norex': dict(rex=0.0), 'snr': dict(stop='snr'), 'nocleanup': dict(cleanup=False)}


def _variant_job(args):
    (name, hp, rex, pxy), truth = args
    patch = load_patch(*pxy)
    t = A.truth_for_patch(truth, pxy[0], pxy[1])
    ov = dict(VARIANTS[name])
    dets, secs, n, info = run_path(patch, ov.pop('hp', hp), ov.pop('rex', rex), **ov)
    return {k: score_dets(d, t) for k, d in dets.items()}, secs, n, info


def part_V(workers=3, out_path=RAW):
    """Design variants of the path on the benchmark's test patch (the frozen sigma_b, r_ex of the development patches): the accumulated dual replaced by the instantaneous
    correlation (greedy order, with the instantaneous-significance stop that it shares with the dual order's 'snr' row), no high-pass, no exclusion zone, the instantaneous significance as the
    stop, no clean-up by the units' own significance. DESCRIPTIVE: every variant is run
    at the frozen sigma_b, r_ex and kappa (kappa is not re-chosen per variant, unlike the SMLM ablations), on one patch: a drop measures "this component, at a kappa tuned without it", not
    a necessity (referee round IS1a)."""
    from multiprocessing import Pool
    raw = load_raw(out_path)
    check_freeze(raw)
    hp, rex, kap, f = choose(raw)
    truth = truth_arrays()
    jobs = [((name, hp, rex, TESTS[0]), truth) for name in VARIANTS]
    with Pool(workers) as pool:
        res = pool.map(_variant_job, jobs, chunksize=1)
    arrays = {}
    for ((name, _, _, _), _t), (scores, secs, n, info) in zip(jobs, res):
        for k, s_ in scores.items():
            arrays.update(_pack('%s|k%g' % (name, k), s_))
        arrays['%s|secs' % name] = np.array([secs, n])
        print('V %s: %d entries, %.0f s' % (name, n, secs), flush=True)
    save_part(out_path, 'V|x%d_y%d' % TESTS[0], arrays)


def _sha256_file(path):
    import hashlib
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def frozen_code_hash():
    """The hash of what decides a test number in this script: the source of the functions that run the path, score a catalogue and choose the configuration and the tuned baselines, and
    the constants (patches, grids, default, baseline families). Edits of the report or of the parts' bookkeeping do not change it."""
    import hashlib
    import inspect
    src = ''.join(inspect.getsource(fn) for fn in (run_path, load_patch, score_dets, _pack, choose, choose_baselines, f_measure, baseline_scores, cs_dets, cs_name))
    consts = repr((TESTS, DEV, DEFAULT, KAPPAS, HPS, REXS, LAT_H, MARGIN, FIT_RADIUS, MAX_ENTRIES, PHOT_TAGS, CS_FILTERS, SECOND))
    return hashlib.sha256((src + consts).encode()).hexdigest()


def dev_keys_hash(raw):
    """The hash of the development results in a raw (every key of the parts D and DB, in sorted order)."""
    import hashlib
    h = hashlib.sha256()
    for k in sorted(k for k in raw if k.startswith('D|') or k.startswith('DB|')):
        h.update(k.encode())
        h.update(np.ascontiguousarray(raw[k]).tobytes())
    return h.hexdigest()


def freeze_record(raw):
    """What is frozen when the development parts have run: the configuration, the tuned baselines, the hashes of the code and of the inputs that decide a test number."""
    hp, rex, kap, f = choose(raw)
    rec = dict(config=dict(sigma_b=hp, r_ex=rex, kappa=kap, dev_F=round(f, 6)), default=list(DEFAULT), tuned_baselines=choose_baselines(raw), tests=[list(p) for p in TESTS],
               files={f_: _sha256_file(os.path.join(REPO, f_)) for f_ in FREEZE_FILES}, truth=_sha256_file(os.path.join(ROOT, 'truth_hugs_meth1.npz')),
               truth_second={k: _sha256_file(os.path.join(v['root'], v['truth'])) for k, v in SECOND.items() if os.path.exists(os.path.join(v['root'], v['truth']))},
               frozen_code=frozen_code_hash(), dev_keys=dev_keys_hash(raw))
    return rec


def write_freeze(raw, path=None):
    import json
    path = path or FREEZE
    rec = freeze_record(raw)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as f:
        json.dump(rec, f, indent=1, sort_keys=True)
        f.write('\n')
    return rec


def check_freeze(raw, path=None):
    """Raises unless the freeze record exists and equals what the code, the inputs and the development results are now (a test number is never produced from a changed method)."""
    import json
    path = path or FREEZE
    if not os.path.exists(path):
        raise RuntimeError('no freeze record %s: run `exp_astro_m2.py --freeze` after the parts D and DB' % path)
    with open(path) as f:
        old = json.load(f)
    now = freeze_record(raw)
    diff = [k for k in now if now[k] != old.get(k)]
    if diff:
        raise RuntimeError('the freeze record differs from the present code / inputs / development results in: %s' % ', '.join(diff))
    return old


def _test_patch_arrays(args):
    """The arrays stored for one test patch (see `_run_test_patches`); a function of its arguments only, so that patches can run in separate processes."""
    pxy, ds, hp, rex, kap, tuned = args
    truth = truth_arrays(ds)
    patch = load_patch(*pxy, ds=ds)
    t = A.truth_for_patch(truth, pxy[0], pxy[1])
    arrays = {'cfg': np.array([hp, rex, kap])}
    runs = {}
    for label, (h, r) in (('ours', (hp, rex)), ('default', (DEFAULT[0], DEFAULT[1]))):
        if (h, r) not in runs:
            dets, secs, n, info = run_path(patch, h, r)
            runs[(h, r)] = (dets, secs, n, info)
            print('T patch %s %s sigma_b %g r_ex %g: %d entries in %.0f s, BIC block %d of %d' % (ds, pxy, h, r, n, secs, info['bic_index'], info['n_blocks']), flush=True)
        dets, secs, n, info = runs[(h, r)]
        arrays['%s|secs' % label] = np.array([secs, n])
        arrays['%s|bic' % label] = np.array([info['bic_index'], info['n_blocks'], int(n >= MAX_ENTRIES)])
        for k, d in dets.items():
            arrays.update(_pack('%s|k%g' % (label, k), score_dets(d, t)))
        dk = kap if label == 'ours' else DEFAULT[2]
        arrays['%s|det' % label] = np.c_[dets[dk][0], dets[dk][1], dets[dk][2]]            # the catalogue at the reported kappa: x, y (patch pixels), flux (DN)
    bs = baseline_scores(pxy, t, BASE if ds == 'm2' else BASE2)
    for label, name in (('photutils', 'photutils'), ('crowdsource', cs_name((0.0, 0.0, 0.0, np.inf)))):
        if name in bs:
            arrays.update(_pack(label + ' (default)', bs[name]))
    if 'daophot' in bs:
        arrays.update(_pack('daophot (published)', bs['daophot']))
    if 'daophot_asis' in bs:
        arrays.update(_pack('daophot (as shipped)', bs['daophot_asis']))
    for label, name in tuned.items():
        if name in bs:
            arrays.update(_pack(label + ' (tuned)', bs[name]))
            arrays['%s (tuned)|name' % label] = np.array([name])
    return arrays


def _run_test_patches(out_path, plist, key_fmt, ds, workers=1):
    """Each patch of `plist` of the data set `ds`, once, with the configuration frozen on the M2 development patches (and at DEFAULT) and the baselines at their defaults and as tuned there.
    Refuses to run unless the freeze record matches (`check_freeze`); a patch whose result exists is not run again. Results are stored under `key_fmt % patch`, by this process only
    (the patches may run in `workers` processes)."""
    raw = load_raw(out_path)
    check_freeze(raw)
    hp, rex, kap, f = choose(raw)
    tuned = choose_baselines(raw)
    todo = []
    for pxy in plist:
        if key_fmt % pxy + '|cfg' in raw:
            print('T patch %s %s: a result exists, not run again' % (ds, pxy), flush=True)
        else:
            todo.append(pxy)
    jobs = [(pxy, ds, hp, rex, kap, tuned) for pxy in todo]
    if workers > 1 and len(jobs) > 1:
        from multiprocessing import Pool
        with Pool(min(workers, len(jobs))) as pool:
            results = pool.map(_test_patch_arrays, jobs, chunksize=1)
    else:
        results = [_test_patch_arrays(j) for j in jobs]
    for pxy, arrays in zip(todo, results):
        save_part(out_path, key_fmt % pxy, arrays)
    print('T %s: sigma_b %g r_ex %g kappa %g (dev F %.3f), default %s' % (ds, hp, rex, kap, f, DEFAULT), flush=True)


def part_T(out_path=RAW, workers=1):
    """The M2 test patches (TESTS), once each, at the configuration frozen on the development patches."""
    _run_test_patches(out_path, TESTS, 'T|x%d_y%d', 'm2', workers)


def part_T2(ds, out_path=RAW, workers=1):
    """The patches of a second cluster (SECOND[ds]), once each, with the SAME frozen configuration and the baselines' tuning of the M2 development patches (nothing is re-tuned)."""
    _run_test_patches(out_path, SECOND[ds]['patches'], 'T2|' + ds + '|x%d_y%d', ds, workers)


# ----------------------------------------------------------------------------------------------------------------------------------
# report
# ----------------------------------------------------------------------------------------------------------------------------------
def table(header, rows):
    out = ['| ' + ' | '.join(header) + ' |', '|' + '|'.join(['---'] * len(header)) + '|']
    return '\n'.join(out + ['| ' + ' | '.join(str(c) for c in r) + ' |' for r in rows])


def cpf(tot):
    """Completeness, purity and their F-measure from totals [true22, matched22, det22, det22_matched] (any leading shape)."""
    tot = np.asarray(tot, float)
    c = tot[..., 1] / np.maximum(tot[..., 0], 1)
    p = tot[..., 3] / np.maximum(tot[..., 2], 1)
    return c, p, 2 * c * p / np.maximum(c + p, 1e-12)


def boot_idx(nb, B=BOOT, seed=0):
    """Block resamples (B, nb) shared by every method of a table: the differences are paired."""
    return np.random.default_rng(seed).integers(0, nb, size=(B, nb))


def boot_idx_strat(sizes, B=BOOT, seed=0):
    """Block resamples (B, sum(sizes)) that resample the blocks WITHIN each patch (a stratified bootstrap: every patch keeps its blocks' number and its weight)."""
    rng = np.random.default_rng(seed)
    out, off = [], 0
    for n in sizes:
        out.append(off + rng.integers(0, n, size=(B, n)))
        off += n
    return np.concatenate(out, axis=1)


def fmt_iv(v, lohi, nd=3):
    return '%.*f [%.*f, %.*f]' % (nd, v, nd, lohi[0], nd, lohi[1])


def boot_cpf(blk, idx):
    """((c, ci), (p, ci), (f, ci)) of a method's blocks (nb, 4) with 95 % percentile intervals over the resamples idx."""
    c0, p0, f0 = cpf(blk.sum(0))
    cs, ps, fs = cpf(blk[idx].sum(1))
    return [(v0, np.percentile(v, [2.5, 97.5])) for v0, v in ((c0, cs), (p0, ps), (f0, fs))]


def boot_diff(blk_a, blk_b, idx):
    """The paired differences (a - b) of completeness, purity and F with 95 % percentile intervals."""
    out = []
    a0, b0 = cpf(blk_a.sum(0)), cpf(blk_b.sum(0))
    aa, bb = cpf(blk_a[idx].sum(1)), cpf(blk_b[idx].sum(1))
    for i in range(3):
        d = aa[i] - bb[i]
        out.append((a0[i] - b0[i], np.percentile(d, [2.5, 97.5])))
    return out


def verdict(d, lohi):
    return 'win' if lohi[0] > 0 else ('loss' if lohi[1] < 0 else 'tie')


ROW_LABELS = {'daophot (published)': 'DAOPHOT catalogue of the benchmark, y shifted by -1 px (the shift was fitted to the HUGS truth of this patch)',
              'daophot (as shipped)': 'DAOPHOT catalogue of the benchmark as shipped (no shift)'}


def patch_methods(raw, pre, kap):
    """[(label, key)] of the rows of one test patch: the path (tuned, a priori default) and the baselines, those present in the raw."""
    dk = DEFAULT[2]
    rows = [('sliding ISS, tuned on the development patches (kappa %g)' % kap, 'ours|k%g' % kap),
            ('sliding ISS, a priori default (sigma_b %g, r_ex %g, kappa %g)' % DEFAULT, 'default|k%g' % dk),
            ('sliding ISS, a priori default with the BIC stop (nothing to tune)', 'default|k0'),
            ('sliding ISS, tuned sigma_b, r_ex with the BIC stop', 'ours|k0')]
    for m in ('photutils (default)', 'photutils (tuned)', 'crowdsource (default)', 'crowdsource (tuned)', 'daophot (published)', 'daophot (as shipped)'):
        if pre + m + '|tot' in raw:
            rows.append((ROW_LABELS.get(m, m) + (' [%s]' % str(raw[pre + m.split(' ')[0] + ' (tuned)|name'][0]) if '(tuned)' in m and pre + m.split(' ')[0] + ' (tuned)|name' in raw else ''), m))
    return rows


def report(raw):
    L = ['# The sliding inverse scale space on the M2 crowded-field patches (SDSS r band) against HST', '',
         'Generated by `experiments/exp_astro_m2.py` from `m2_raw.npz`; every number below is printed by the script. Protocol, model and the benchmark\'s definition: the script\'s docstring and '
         '`experiments/astro_common.py` (matching at 0.5 px and 0.5 mag; truth sources within 2.5 px of the edge are not scored; bins of 0.5 mag from 15.0 to 23.5). Truth: the HUGS NGC 7089 '
         'catalogue, F606W used as r. Completeness and purity at mag < 22 (the benchmark\'s depth) with 95 %% percentile intervals of a bootstrap over the 20 x 20 pixel blocks of a patch '
         '(%d resamples, seed 0, the same resampled blocks for every method of a table, so differences are paired). The first choice of test patches, (410, 630) and (410, 730), was replaced before any test patch was '
         'scored because they lay mostly outside the HUGS footprint (`astro_patch_survey.py`); no test result exists for them.' % BOOT, '']
    if 'D|%s' % ('x%d_y%d|hp%g|rex%g|k%g|tot' % (DEV[0][0], DEV[0][1], HPS[0], REXS[0], KAPPAS[0])) in raw:
        hp, rex, kap, f = choose(raw)
        L += ['## Development patches %s: sigma_b, r_ex, kappa\n' % ', '.join('(%d, %d)' % p for p in DEV),
              'F-measure of completeness and purity at mag < 22 (the two patches\' counts added); bold: the point frozen for the test patches.\n']
        rows = []
        for h in HPS:
            for r in REXS:
                capped = sum(int(raw['D|x%d_y%d|hp%g|rex%g|bic' % (p[0], p[1], h, r)][2]) for p in DEV if 'D|x%d_y%d|hp%g|rex%g|bic' % (p[0], p[1], h, r) in raw)
                for k in KAPPAS:
                    tot = sum(raw['D|x%d_y%d|hp%g|rex%g|k%g|tot' % (p[0], p[1], h, r, k)] for p in DEV)
                    ff, c, pu = f_measure(tot)
                    mark = '**' if (h, r, k) == (hp, rex, kap) else ''
                    rows.append(['%s%g%s' % (mark, h, mark), '%s%g%s' % (mark, r, mark), '%s%g%s' % (mark, k, mark), '%s%.3f%s' % (mark, c, mark), '%s%.3f%s' % (mark, pu, mark),
                                 '%s%.3f%s' % (mark, ff, mark), '%d' % tot[2], '%d' % (tot[6] - tot[5]), ('%d of %d' % (capped, len(DEV))) if capped else '-'])
        L.append(table(['sigma_b px', 'r_ex px', 'kappa', 'completeness (mag < 22)', 'purity (mag < 22)', 'F', 'detections < 22 mag', 'unmatched detections (all magnitudes)',
                        'patches whose path hit the entry cap (%d)' % MAX_ENTRIES], rows))
        L.append('\nWhere a path hit the entry cap its rows show the last block of a path cut at the cap, not the F of that significance level: they are censored and say little about the configuration.\n')
        edges = [('sigma_b', hp, HPS), ('r_ex', rex, REXS), ('kappa', kap, KAPPAS)]
        edges = ['%s %g is at the %s edge of its grid %s' % (nm, v, 'lower' if v == min(g) else 'upper', tuple(g)) for nm, v, g in edges if v in (min(g), max(g))]
        if edges:
            L.append('The frozen configuration is at the edge of the grid on every axis it was searched on where the next value was not tried: %s.\n' % '; '.join(edges))
        bic_keys = ['D|x%d_y%d|hp%g|rex%g|bic' % (p[0], p[1], hp, rex) for p in DEV]
        if all(k_ in raw for k_ in bic_keys):
            L.append('The BIC block of the path at the frozen (sigma_b, r_ex), per development patch: %s. The path was run to the kappa = %g time, so the BIC readout is the minimum of the path as '
                     'truncated there; "interior" means the minimum is before the last block.\n' % ('; '.join('block %d of %d (%s%s)' % (
                         raw[k_][0], raw[k_][1], 'interior' if raw[k_][0] < raw[k_][1] - 1 else 'at the end of the path', ', capped at %d entries' % MAX_ENTRIES if raw[k_][2] else '')
                         for k_ in bic_keys), min(KAPPAS)))
        tags = [t for t in PHOT_TAGS if all('DB|x%d_y%d|%s|tot' % (p[0], p[1], t) in raw for p in DEV)]
        if tags:
            L += ['### Baselines on the development patches\n', 'The same F-measure; the variant of each baseline family with the largest F (bold) is the one used as "tuned" on the test patches.\n']
            tuned = choose_baselines(raw)
            rows = []
            for name in tags + [cs_name(f_) for f_ in CS_FILTERS if all('DB|x%d_y%d|%s|tot' % (p[0], p[1], cs_name(f_)) in raw for p in DEV)]:
                tot = sum(raw['DB|x%d_y%d|%s|tot' % (p[0], p[1], name)] for p in DEV)
                ff, c, pu = f_measure(tot)
                mark = '**' if name in tuned.values() else ''
                rows.append([mark + name + mark, '%.3f' % c, '%.3f' % pu, '%.3f' % ff, '%d' % tot[2], '%d' % (tot[6] - tot[5])])
            L.append(table(['baseline variant', 'completeness (mag < 22)', 'purity (mag < 22)', 'F', 'detections < 22 mag', 'unmatched detections (all magnitudes)'], rows))
            L.append('')
    groups = [('M2 (NGC 7089), SDSS run 2583', 'T|x%d_y%d|', TESTS)] + [(v['title'], 'T2|%s|x%%d_y%%d|' % k, v['patches']) for k, v in SECOND.items()]
    for gi, (gtitle, kfmt, plist) in enumerate(groups):
        pooled, summary = {}, []
        for pxy in plist:
            pre = kfmt % pxy
            if pre + 'cfg' not in raw:
                continue
            hp, rex, kap = raw[pre + 'cfg']
            meths = patch_methods(raw, pre, kap)
            blks = {k: raw[pre + k + '|blk22'] for l, k in meths}
            idx = boot_idx(next(iter(blks.values())).shape[0])
            tr = raw[pre + 'ours|k%g|tot' % kap]
            empty = int((raw[pre + 'ours|k%g|blk22' % kap][:, 0] == 0).sum())
            bic = raw.get(pre + 'ours|bic')
            L += ['## %s: test patch (%d, %d)\n' % ((gtitle,) + pxy), 'Frozen from the development patches: sigma_b = %g px, r_ex = %g px, kappa = %g. Truth sources brighter than 22 mag in the scored area: %d (blocks of 20 x 20 px '
                  'without such a source: %d of %d). Run time of the path to its last entrant: %.0f s (%d entrants%s), one core, shared machine.%s\n' % (
                      hp, rex, kap, tr[0], empty, raw[pre + 'ours|k%g|blk22' % kap].shape[0], raw[pre + 'ours|secs'][0], raw[pre + 'ours|secs'][1],
                      ', capped at %d' % MAX_ENTRIES if (bic is not None and bic[2]) else '',
                      (' BIC readout: block %d of %d (%s).' % (bic[0], bic[1], 'interior' if bic[0] < bic[1] - 1 else 'the end of the path: the path was cut before the BIC minimum')) if bic is not None else '')]
            rows = []
            for label, key in meths:
                tot = raw[pre + key + '|tot']
                (c, cci), (pu, pci), (ff, fci) = boot_cpf(blks[key], idx)
                rows.append([label, '%d' % tot[1], fmt_iv(c, cci), '%d' % tot[2], '%d' % tot[3], fmt_iv(pu, pci), fmt_iv(ff, fci), '%d' % tot[6], '%d' % (tot[6] - tot[5]),
                             '(%.2f, %.2f)' % tuple(raw[pre + key + '|offset'])])
                pooled.setdefault(key, []).append(blks[key])
            L.append(table(['catalogue', 'matched < 22', 'completeness < 22', 'detections < 22', 'matched', 'purity < 22', 'F', 'detections (all mags)', 'unmatched (all mags)',
                            'median offset of bright matches (px)'], rows))
            L.append('\nPaired differences (sliding ISS minus the baseline, same blocks), 95 % intervals:\n')
            rows = []
            for ours_label, ours_key, base in (('tuned ISS', 'ours|k%g' % kap, ('photutils (tuned)', 'crowdsource (tuned)')), ('default ISS', 'default|k%g' % DEFAULT[2], ('photutils (default)', 'crowdsource (default)')),
                                               ('BIC-stop default ISS', 'default|k0', ('photutils (default)', 'crowdsource (default)', 'photutils (tuned)', 'crowdsource (tuned)'))):
                for b in base:
                    if b not in blks:
                        continue
                    dc, dp, df = boot_diff(blks[ours_key], blks[b], idx)
                    rows.append(['%s vs %s' % (ours_label, b), fmt_iv(dc[0], dc[1]), fmt_iv(dp[0], dp[1]), fmt_iv(df[0], df[1]), verdict(*df)])
                    summary.append((pxy, ours_label, b, df[0], df[1]))
            L.append(table(['comparison', 'd completeness', 'd purity', 'd F', 'verdict on F'], rows))
            L.append('\nThe path at the other significance levels (the same frozen sigma_b, r_ex; descriptive: no choice rests on these rows):\n')
            rows = []
            for k in KAPPAS:
                tot = raw[pre + 'ours|k%g|tot' % k]
                (c, cci), (pu, pci), (ff, fci) = boot_cpf(raw[pre + 'ours|k%g|blk22' % k], idx)
                rows.append(['%g' % k, '%d' % tot[1], fmt_iv(c, cci), '%d' % tot[2], fmt_iv(pu, pci), fmt_iv(ff, fci)])
            L.append(table(['kappa', 'matched < 22', 'completeness < 22', 'detections < 22', 'purity < 22', 'F'], rows))
            if pxy == plist[0]:
                edges = A.BIN_EDGES
                keys = [(l, k) for l, k in meths]
                L.append('\nCompleteness by magnitude (matched truth / truth; the truth counts are the first column):\n')
                hdr = ['mag bin', 'truth'] + [l.split(' [')[0].replace('sliding ISS, ', 'ISS ') for l, k in keys]
                rows = []
                for b in range(len(edges) - 1):
                    nt = int(raw[pre + '%s|n_true' % keys[0][1]][b])
                    if nt == 0:
                        continue
                    rows.append(['%.1f-%.1f' % (edges[b], edges[b + 1]), str(nt)] + ['%.2f' % raw[pre + '%s|completeness' % k][b] for l, k in keys])
                L.append(table(hdr, rows))
                L.append('\nFalse-positive fraction by detected magnitude (unmatched / detected, number detected in brackets):\n')
                hdr = ['mag bin'] + [l.split(' [')[0].replace('sliding ISS, ', 'ISS ') for l, k in keys]
                rows = []
                for b in range(len(edges) - 1):
                    nds = [int(raw[pre + '%s|n_det' % k][b]) for l, k in keys]
                    if sum(nds) == 0:
                        continue
                    rows.append(['%.1f-%.1f' % (edges[b], edges[b + 1])] + ['%.2f (%d)' % (raw[pre + '%s|fp_fraction' % k][b], nd) if nd else '-' for (l, k), nd in zip(keys, nds)])
                L.append(table(hdr, rows))
            L.append('')
            vpre = 'V|x%d_y%d|' % pxy
            if gi == 0 and pxy == plist[0] and any(k.startswith(vpre) for k in raw):
                L.append('### Design variants of the path on this patch (frozen sigma_b, r_ex and kappa; the same blocks)\n\nDescriptive only: kappa is NOT re-chosen for a variant (the SMLM ablations do re-choose it), the variants are run on this one patch, and removing a component changes what the frozen kappa means; a drop is "this component, at a kappa tuned without it", not a measure of necessity.\n')
                main_blk = raw[pre + 'ours|k%g|blk22' % kap]
                rows = []
                (c, cci), (pu, pci), (ff, fci) = boot_cpf(main_blk, idx)
                rows.append(['the path as frozen', fmt_iv(c, cci), fmt_iv(pu, pci), fmt_iv(ff, fci), '-', '%d' % raw[pre + 'ours|secs'][1], '%.0f' % raw[pre + 'ours|secs'][0]])
                for name in VARIANTS:
                    if vpre + '%s|k%g|blk22' % (name, kap) not in raw:
                        continue
                    vb = raw[vpre + '%s|k%g|blk22' % (name, kap)]
                    (c, cci), (pu, pci), (ff, fci) = boot_cpf(vb, idx)
                    dc, dp, df = boot_diff(vb, main_blk, idx)
                    rows.append([{'greedy_snr': 'instantaneous correlation instead of the accumulated dual, both with the instantaneous-significance stop (compare the row of that stop)', 'nohp': 'no high-pass', 'norex': 'no exclusion zone',
                                  'snr': 'instantaneous significance as the stop (the SMLM rule)', 'nocleanup': 'no clean-up by the units\' own significance'}[name], fmt_iv(c, cci), fmt_iv(pu, pci),
                                 fmt_iv(ff, fci), fmt_iv(df[0], df[1]), '%d' % raw[vpre + '%s|secs' % name][1], '%.0f' % raw[vpre + '%s|secs' % name][0]])
                L.append(table(['variant', 'completeness < 22', 'purity < 22', 'F', 'dF against the frozen path (paired)', 'entrants', 'path seconds'], rows))
                L.append('')
        reported = [p for p in plist if (kfmt % p) + 'cfg' in raw]
        # only the catalogues present for every reported patch are pooled (the shipped DAOPHOT catalogue covers the benchmark patch alone)
        pooled = {k: v for k, v in pooled.items() if len(v) == len(reported)}
        if len(reported) > 1 and pooled:
            npatch = len(reported)
            blk = {k: np.concatenate(v) for k, v in pooled.items()}
            nb = next(iter(blk.values())).shape[0]
            L += ['## %s: all %d test patches together\n' % (gtitle, npatch), 'The blocks of the patches pooled (%d blocks); the counts are the sums; the bootstrap resamples the blocks WITHIN each patch (stratified), '
                  'so that no patch changes its weight. The patches are windows of one frame of one cluster, not independent trials. Catalogues that exist for one patch only '
                  '(the DAOPHOT catalogue) are not pooled.\n' % nb]
            idx = boot_idx_strat([len(v_) for v_ in next(iter(pooled.values()))])
            kap = float(raw[(kfmt % reported[0]) + 'cfg'][2])
            rows = []
            for label, key in patch_methods(raw, kfmt % reported[0], kap):
                if key not in blk:
                    continue
                tot = blk[key].sum(0)
                (c, cci), (pu, pci), (ff, fci) = boot_cpf(blk[key], idx)
                tot_all = sum(raw[(kfmt % p_) + key + '|tot'] for p_ in reported)            # the totals of all magnitudes (the blocks count mag < 22 only)
                rows.append([label.split(' [')[0], '%d / %d' % (tot[1], tot[0]), fmt_iv(c, cci), '%d / %d' % (tot[3], tot[2]), fmt_iv(pu, pci), fmt_iv(ff, fci),
                             '%.3f (%d)' % (tot_all[5] / max(tot_all[6], 1), tot_all[6] - tot_all[5])])
            L.append(table(['catalogue', 'matched / truth < 22', 'completeness', 'matched / detected < 22', 'purity', 'F', 'purity over all magnitudes (unmatched detections)'], rows))
            L.append('\nPaired differences on the pooled blocks:\n')
            rows = []
            for ours_label, ours_key, base in (('tuned ISS', 'ours|k%g' % kap, ('photutils (tuned)', 'crowdsource (tuned)')), ('default ISS', 'default|k%g' % DEFAULT[2], ('photutils (default)', 'crowdsource (default)')),
                                               ('BIC-stop default ISS', 'default|k0', ('photutils (default)', 'crowdsource (default)', 'photutils (tuned)', 'crowdsource (tuned)'))):
                for b in base:
                    if b in blk and ours_key in blk:
                        dc, dp, df = boot_diff(blk[ours_key], blk[b], idx)
                        rows.append(['%s vs %s' % (ours_label, b), fmt_iv(dc[0], dc[1]), fmt_iv(dp[0], dp[1]), fmt_iv(df[0], df[1]), verdict(*df)])
            L.append(table(['comparison', 'd completeness', 'd purity', 'd F', 'verdict on F'], rows))
            L.append('\nRecord by patch (verdict on F, per comparison):\n')
            cmp_names = sorted({(a, b) for _, a, b, _, _ in summary})
            rows = [['%s vs %s' % (a, b)] + [verdict(d, ci) for (pxy, a_, b_, d, ci) in summary if (a_, b_) == (a, b)] +
                    ['%d / %d / %d' % tuple([sum(1 for (_, a_, b_, d, ci) in summary if (a_, b_) == (a, b) and verdict(d, ci) == v) for v in ('win', 'tie', 'loss')])] for a, b in cmp_names]
            L.append(table(['comparison'] + ['(%d, %d)' % p for p in plist if (kfmt % p) + 'cfg' in raw] + ['wins / ties / losses'], rows))
            L.append('')
    return '\n'.join(L) + '\n'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--only', choices=['D', 'DB', 'T', 'T2', 'V'])
    ap.add_argument('--freeze', action='store_true', help='write results/astro/freeze.json from the development parts D and DB')
    ap.add_argument('--workers', type=int, default=int(os.environ.get('ASTRO_WORKERS', '3')))
    ap.add_argument('--from-raw', action='store_true')
    A_ = ap.parse_args()
    if A_.from_raw:
        with open(MD, 'w') as f:
            f.write(report(load_raw()))
        print('wrote', MD)
    elif A_.freeze:
        print(json.dumps(write_freeze(load_raw()), indent=1, sort_keys=True))
    elif A_.only == 'D':
        part_D(A_.workers)
    elif A_.only == 'DB':
        part_DB()
    elif A_.only == 'T':
        part_T(workers=A_.workers)
    elif A_.only == 'T2':
        for ds_ in SECOND:
            part_T2(ds_, workers=A_.workers)
    elif A_.only == 'V':
        part_V(A_.workers)
    else:
        ap.error('give --only PART, --freeze or --from-raw')


if __name__ == '__main__':
    main()
