"""
Why does the ISS lead over the front-end crowdsource depend on the cluster? A post hoc diagnosis of the spent patches, one label-free repair, and its held-out test on six clusters no earlier test used
=========================================================================================================================================================================================================

Parts (all numbers of `results/astro/clusterdep.md` are printed from `clusterdep_raw.npz` alone: `--from-raw`):

  --diag     POST HOC. On the 19 spent test patches (M2 four windows and three repeat frames, M53, M92, M3, M15) the stored catalogues of the guarded ISS (`cyclefix_raw`, `repeat_raw`, `fresh_raw`: the isotropic-kernel row
             and the kernel-rule row at the tuned configuration, kappa 2) and the front-end crowdsource (`baseline_catalogs_hp2.npz`, `fresh_baseline_catalogs.npz`, the tag `cs2_hp2_t3_pm`) are read, never re-run. Per patch:
             covariates from the pixels, the PSF stamp and the truth list, F of every row, the completeness and purity by magnitude, the decomposition of the errors, and the F of the ISS catalogue after the merge of
             units that the kernel cannot tell apart (the repair, `astro_clusterdep.merge_coherent`) at a grid of coherence thresholds. The cross-patch statistics (rank correlations of the gap with each covariate, leave-one-
             cluster-out prediction) are computed in the report from the stored per-patch values. Nothing in this part is a test.
  --freeze   writes `results/astro/freeze_clusterdep.json` (hashes of the plan, the code, the registry, the selection and the baselines' catalogues) -- after the plan is committed, before `--run`.
  --run      the held-out test of `results/astro/clusterdep_plan.md`: the frozen ISS (isotropic kernel and the kernel the rule picks from the stamp, tuned configuration, guarded path) on the patches of
             `clusterdep_patches.json`, once each, write-once; the repair applied to the catalogues; every row scored against the truth shifted by the measured rigid offset (primary) and as it is (secondary), together
             with the baselines' catalogues of `clusterdep_baseline_catalogs.npz` (no truth was read when they were made).
  --r1       POST HOC (after the held-out run): the first attempt at a repair (a noise model with a model-error floor) run on seven spent cases; reported in section 2b, not part of the test.
  --from-raw regenerates `results/astro/clusterdep.md` from `clusterdep_raw.npz` (reads nothing else).

Environment (tests): CLUSTERDEP_OUT (directory of the raw, the md and the freeze record), CLUSTERDEP_PLAN, CLUSTERDEP_REGISTRY, CLUSTERDEP_BASE.
"""
import argparse
import hashlib
import json
import os
import sys
import time

os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np

HERE = os.path.abspath(os.path.dirname(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)
import astro_clusterdep as C              # noqa: E402
import astro_common as A                  # noqa: E402
import astro_fresh as F                   # noqa: E402
import astro_repeat as R                  # noqa: E402
import astro_second as S2                 # noqa: E402
import exp_astro_cyclefix as CF           # noqa: E402
import exp_astro_fresh as FR              # noqa: E402
import exp_astro_m2 as X                  # noqa: E402
import exp_astro_psf as PS                # noqa: E402
from src.core import imaging as IM        # noqa: E402
from src.core import psf_general as PG    # noqa: E402

ASTRO = os.path.join(REPO, 'results', 'astro')
OUT = os.environ.get('CLUSTERDEP_OUT', ASTRO)
RAW, MD = os.path.join(OUT, 'clusterdep_raw.npz'), os.path.join(OUT, 'clusterdep.md')
FREEZE = os.path.join(OUT, 'freeze_clusterdep.json')
PLAN = os.environ.get('CLUSTERDEP_PLAN', os.path.join(ASTRO, 'clusterdep_plan.md'))
REGISTRY = os.environ.get('CLUSTERDEP_REGISTRY', os.path.join(ASTRO, 'clusterdep_patches.json'))
BASE = os.environ.get('CLUSTERDEP_BASE', os.path.join(ASTRO, 'clusterdep_baseline_catalogs.npz'))
M2FREEZE = os.path.join(ASTRO, 'freeze.json')
CYCLEFIX_RAW, REPEAT_RAW, FRESH_RAW = (os.path.join(ASTRO, n) for n in ('cyclefix_raw.npz', 'repeat_raw.npz', 'fresh_raw.npz'))
HP2_CATS, FRESH_CATS = os.path.join(ASTRO, 'baseline_catalogs_hp2.npz'), os.path.join(ASTRO, 'fresh_baseline_catalogs.npz')
FILES = ('src/core/imaging.py', 'src/core/psf_general.py', 'src/core/imaging_guard.py', 'experiments/astro_common.py', 'experiments/astro_second.py', 'experiments/astro_fresh.py', 'experiments/exp_astro_m2.py',
         'experiments/exp_astro_psf.py', 'experiments/exp_astro_cyclefix.py', 'experiments/exp_astro_fresh.py', 'experiments/astro_clusterdep.py', 'experiments/clusterdep_select.py',
         'experiments/exp_astro_clusterdep_baselines.py', 'results/astro/clusterdep_plan.md', 'results/astro/clusterdep_selection.txt', 'results/astro/clusterdep_patches.json',
         'results/astro/clusterdep_baseline_catalogs.npz', 'results/astro/freeze.json', 'data/astro/checksums_clusterdep.json')

RHO_GRID = (0.6, 0.7, 0.75, 0.8, 0.85, 0.9)          # coherence thresholds of the merge scanned on the spent patches
RHO_STAR = 0.8                                       # the threshold of the repair (declared in the plan before the held-out run)
HP, REX, KAP = 3.0, 1.2, 2.0                         # the frozen tuned configuration (freeze.json)
ROWS = ('iso', 'rule')                               # the ISS rows: the isotropic kernel as frozen; the kernel of `psf_general.fit_adequate` (identical where the stamp is adequately isotropic)
CLUSTER_OF = {'M2': 'M2', 'M2r1': 'M2', 'M2r2': 'M2', 'M2r3': 'M2', 'M53': 'M53', 'M92': 'M92', 'M3': 'M3', 'M15': 'M15'}
COVARIATES = ('fwhm', 'stamp_rms', 'stamp_ell', 'stamp_centroid', 'sky', 'noise', 'snr22', 'bg_grad', 'dens22', 'n22', 'f15', 'f25', 'rsd_iso', 'rsd_rule', 'comp2122_cs', 'comp2122_iso')
COV_NOTE = {'fwhm': 'FWHM of the psField stamp (px)', 'stamp_rms': 'rms misfit of the round 3-Gaussian fit to the stamp / its peak', 'stamp_ell': 'ellipticity of the stamp (second moments, r <= 5 px)',
            'stamp_centroid': 'distance of the stamp centroid from its middle pixel (px)', 'sky': 'sky level (5th percentile of the patch, DN)', 'noise': 'model noise sd at the median pixel (DN)',
            'snr22': 'peak S/N of a 22 mag star', 'bg_grad': 'sd over the 25 blocks of the block 10th percentile / noise (background structure)', 'dens22': 'truth stars brighter than 22 mag per pixel (truth)',
            'n22': 'truth stars brighter than 22 mag in the scored area (truth)', 'f15': 'fraction of those with a brighter truth star within 1.5 FWHM (truth)',
            'f25': 'same within 2.5 FWHM (truth)', 'rsd_iso': 'residual scale (1.4826 MAD of the standardised residual) of the isotropic-kernel ISS catalogue (label-free)',
            'rsd_rule': 'same for the kernel-rule catalogue (label-free)', 'comp2122_cs': 'completeness of the front-end crowdsource at 21 to 22 mag (truth)', 'comp2122_iso': 'same for the ISS isotropic row (truth)'}


def utc():
    return time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


# ----------------------------------------------------------------------------------------------------------------------------------
# the spent patches
# ----------------------------------------------------------------------------------------------------------------------------------
def spent_patches():
    """The 19 spent test patches as dicts: label (M2, M2r1, ..., M15), pxy, patch (dn, psf, gain, dark), truth patch, the stored ISS catalogues 'iso' and 'rule' (n, 3 of x, y, flux, kappa 2, guarded path), 'cs'
    (the front-end crowdsource, `cs2_hp2_t3_pm`). Read only: the frozen raws and catalogue files."""
    cf, rp, fr = np.load(CYCLEFIX_RAW), np.load(REPEAT_RAW), np.load(FRESH_RAW)
    h2, fb = np.load(HP2_CATS), np.load(FRESH_CATS)
    reg = F.load_registry(os.path.join(ASTRO, 'fresh_patches.json'))
    L = []
    for p in X.TESTS:
        k = 'C|T|x%d_y%d|ours|det' % p
        L.append(dict(label='M2', pxy=p, ds='m2', iso=cf[k], rule=cf[k]))
    for ds, lab in (('m53', 'M53'), ('m92', 'M92')):
        for p in X.SECOND[ds]['patches']:
            L.append(dict(label=lab, pxy=p, ds=ds, iso=cf['C|T2|%s|x%d_y%d|ours|det' % ((ds,) + p)], rule=cf['C|P|%s|x%d_y%d|ours|det' % ((ds,) + p)]))
    for n in sorted(R.REPEAT):
        for p in R.REPEAT[n]['patches']:
            L.append(dict(label=n.replace('m2r', 'M2r'), pxy=tuple(p), ds=n, iso=rp['R|%s|x%d_y%d|iso|ours|det' % ((n,) + tuple(p))], rule=rp['R|%s|x%d_y%d|adq|ours|det' % ((n,) + tuple(p))]))
    for c in sorted(reg):
        for p in reg[c]['patches']:
            L.append(dict(label={'ngc5272': 'M3', 'ngc7078': 'M15'}[c], pxy=tuple(p), ds=c, iso=fr['F|%s|x%d_y%d|iso|ours|det' % ((c,) + tuple(p))], rule=fr['F|%s|x%d_y%d|adq|ours|det' % ((c,) + tuple(p))]))
    for d in L:
        p = d['pxy']
        if d['ds'] == 'm2':
            pa, tr = X.load_patch(*p), X.truth_arrays('m2')
        elif d['ds'] in X.SECOND:
            pa, tr = X.load_patch(*p, ds=d['ds']), X.truth_arrays(d['ds'])
        elif d['ds'] in R.REPEAT:
            pa, tr = R.load_patch(p[0], p[1], d['ds']), R.truth_arrays(d['ds'])
        else:
            pa, tr = F.load_patch(p[0], p[1], reg[d['ds']]), F.truth_arrays(reg[d['ds']])
        d['patch'], d['truth'] = pa, A.truth_for_patch(tr, p[0], p[1])
        if d['ds'] in reg:
            k = '%s|x%d_y%d|%s|' % (d['ds'], p[0], p[1], str(fb['meta|crowdsource_fe'][0]))
            d['cs'] = np.c_[fb[k + 'x'], fb[k + 'y'], fb[k + 'flux']]
        else:
            k = 'x%d_y%d|cs2_hp2_t3_pm|' % p
            d['cs'] = np.c_[h2[k + 'x'], h2[k + 'y'], h2[k + 'flux']]
    return L


def pkey(label, pxy):
    return 'D|%s|x%d_y%d|' % (label, pxy[0], pxy[1])


def score_det(det, truth):
    det = np.asarray(det, float).reshape(-1, 3)
    return X.score_dets((det[:, 0], det[:, 1], det[:, 2]), truth)


def f_of(s):
    """(F, completeness, purity) at mag < 22 of a score dict of `astro_common.score`."""
    c = s['matched22'] / max(s['true22'], 1)
    p = s['det22_matched'] / max(s['det22'], 1)
    return 2 * c * p / max(c + p, 1e-12), c, p


# ----------------------------------------------------------------------------------------------------------------------------------
# covariates and the decomposition of the errors
# ----------------------------------------------------------------------------------------------------------------------------------
def stamp_centroid(stamp, radius=6.0):
    st = np.asarray(stamp, float)
    n = st.shape[0]
    yy, xx = np.mgrid[0:n, 0:n]
    c = (n - 1) / 2.0
    w = st * (np.hypot(yy - c, xx - c) <= radius)
    w = w / w.sum()
    return float(np.hypot((w * (xx - c)).sum(), (w * (yy - c)).sum()))


def covariates(d):
    """The label-free covariates (pixels and stamp) and the covariates that read the truth list, of one spent patch, as a dict of floats (COVARIATES)."""
    from scipy.ndimage import gaussian_filter
    pa, t = d['patch'], d['truth']
    dn = pa['dn']
    fwhm = S2.psf_fwhm(pa['psf'])
    fit = PG.fit_adequate(pa['psf'])
    expected = np.maximum(gaussian_filter(dn, 1.5, mode='reflect'), np.percentile(dn, 1))
    noise = float(np.sqrt(np.median(expected) / pa['gain'] + pa['dark']))
    flux22 = A.mag_to_counts(22.0)
    sig, parts, rms_iso = IM.fit_psf_parts(pa['psf'], n_parts=3, sigma0=1.0)
    nb = 5
    q10 = np.array([np.percentile(dn[i * 20:(i + 1) * 20, j * 20:(j + 1) * 20], 10) for i in range(nb) for j in range(nb)])
    s22 = t['inner'] & (t['mag'] < 22)
    m, tx, ty = t['mag'], t['x'], t['y']
    D = np.hypot(tx[:, None] - tx[None, :], ty[:, None] - ty[None, :])
    br = m[None, :] < m[:, None]
    out = dict(fwhm=fwhm, stamp_rms=float(rms_iso), stamp_ell=float(PS.stamp_shape(pa['psf'])[0]), stamp_centroid=stamp_centroid(pa['psf']), sky=float(np.percentile(dn, 5)), noise=noise,
               snr22=float(flux22 * pa['psf'].max() / pa['psf'].sum() / noise), bg_grad=float(q10.std() / noise), dens22=float(s22.sum() / 9500.0), n22=float(s22.sum()),
               f15=float(((D < 1.5 * fwhm) & br & (D > 0)).any(1)[s22].mean()), f25=float(((D < 2.5 * fwhm) & br & (D > 0)).any(1)[s22].mean()))
    for row in ROWS:
        kfit = PG.fit_adequate(pa['psf'], rms_goal=np.inf) if row == 'iso' else fit
        P = C.make_kernel(pa, HP, kfit)
        out['rsd_' + row] = PS.residual_sd(pa, P, d[row])
    return out


def recall_bin(s, lo, hi):
    """Completeness of a score dict over the 0.5 mag bins between lo and hi mag."""
    edges = A.BIN_EDGES
    sel = (edges[:-1] >= lo - 1e-9) & (edges[1:] <= hi + 1e-9)
    return float(s['n_matched'][sel].sum() / max(s['n_true'][sel].sum(), 1))


def decompose(det, t, mag_max=22.0):
    """The errors of a catalogue against the truth at mag < mag_max, counts: [truth stars (scored), matched, unmatched with a detection within 0.5 px (so the magnitude differs by more than 0.5), unmatched with
    the nearest detection at 0.5 to 1.5 px, unmatched with none within 1.5 px, detections (scored, < mag_max), unmatched, of those: nearest truth star of any magnitude within 0.5 px, at 0.5 to 1.5 px, none within 1.5 px]."""
    det = np.asarray(det, float).reshape(-1, 3)
    mg = A.counts_to_mag(det[:, 2])
    ok = np.isfinite(mg) & (det[:, 0] >= A.EDGE - 0.5) & (det[:, 0] <= 99.5 - A.EDGE) & (det[:, 1] >= A.EDGE - 0.5) & (det[:, 1] <= 99.5 - A.EDGE)
    xy, mg = det[ok, :2], mg[ok]
    ti = np.flatnonzero(t['inner'] & (t['mag'] < mag_max))
    T = np.c_[t['x'], t['y']][ti]
    rd, rt = A.match(xy, mg, T, t['mag'][ti])
    out = np.zeros(10)
    out[0], out[1] = len(ti), len(rt)
    if len(xy):
        D = np.hypot(xy[:, None, 0] - T[None, :, 0], xy[:, None, 1] - T[None, :, 1])
        for j in np.setdiff1d(np.arange(len(ti)), rt):
            dm = D[:, j].min()
            out[2 + (0 if dm <= 0.5 else (1 if dm <= 1.5 else 2))] += 1
    else:
        out[4] = len(ti) - len(rt)
    sel = np.flatnonzero(mg < mag_max)
    out[5] = len(sel)
    um = np.setdiff1d(sel, rd)
    out[6] = len(um)
    allT = np.c_[t['x'], t['y']]
    for i in um:
        dm = np.hypot(allT[:, 0] - xy[i, 0], allT[:, 1] - xy[i, 1]).min() if len(allT) else np.inf
        out[7 + (0 if dm <= 0.5 else (1 if dm <= 1.5 else 2))] += 1
    return out


def bright_split(det, t, mag_max=19.0):
    """For the scored truth stars brighter than mag_max: [their number, mean number of detections within 1.5 px, fraction with a detection within 0.5 px, mean distance of the nearest detection,
    fraction matched (0.5 px, 0.5 mag)]."""
    det = np.asarray(det, float).reshape(-1, 3)
    mg = A.counts_to_mag(det[:, 2])
    ti = np.flatnonzero(t['inner'] & (t['mag'] < mag_max))
    if len(ti) == 0 or len(det) == 0:
        return np.array([len(ti), np.nan, np.nan, np.nan, np.nan])
    ok = np.isfinite(mg)
    D = np.hypot(det[ok, None, 0] - t['x'][None, ti], det[ok, None, 1] - t['y'][None, ti])
    rd, rt = A.match(det[ok, :2], mg[ok], np.c_[t['x'], t['y']][ti], t['mag'][ti])
    return np.array([len(ti), float((D < 1.5).sum(0).mean()), float((D.min(0) <= 0.5).mean()), float(D.min(0).mean()), len(rt) / len(ti)])


def diag_patch(d):
    """The arrays of one spent patch (prefix `pkey`): covariates, rows and their scores, decomposition, bright-star split, merged rows."""
    pre = pkey(d['label'], d['pxy'])
    t, pa = d['truth'], d['patch']
    arr = {}
    cov = covariates(d)
    rows = {'iso': d['iso'], 'rule': d['rule'], 'cs': d['cs']}
    sc = {r: score_det(v, t) for r, v in rows.items()}
    cov['comp2122_cs'], cov['comp2122_iso'] = recall_bin(sc['cs'], 21.0, 22.0), recall_bin(sc['iso'], 21.0, 22.0)
    for k in COVARIATES:
        arr[pre + 'cov|' + k] = np.array([cov[k]])
    arr[pre + 'kind'] = np.array([PG.fit_adequate(pa['psf'])['kind']])
    for r, s in sc.items():
        arr.update(X._pack(pre + r, s))
        arr[pre + r + '|dec'] = decompose(rows[r], t)
        arr[pre + r + '|split'] = bright_split(rows[r], t)
        arr[pre + r + '|bins'] = np.c_[s['n_true'], s['n_matched'], s['n_det'], s['n_det_matched']]
        arr[pre + r + '|n'] = np.array([len(rows[r])])
    for r in ROWS:
        kfit = PG.fit_adequate(pa['psf'], rms_goal=np.inf) if r == 'iso' else PG.fit_adequate(pa['psf'])
        off, tab = C.coherence_table(pa, kfit)
        for rho in RHO_GRID:
            m, nm = C.merge_coherent(rows[r], off, tab, rho)
            arr.update(X._pack('%s%s_m%g' % (pre, r, rho), score_det(m, t)))
            arr['%s%s_m%g|merges' % (pre, r, rho)] = np.array([nm])
    return arr


def part_diag(out_path=RAW):
    L = spent_patches()
    arrays = {'meta|labels': np.array(['%s|x%d_y%d' % (d['label'], d['pxy'][0], d['pxy'][1]) for d in L]),
              'meta|inputs_sha256': np.array([sha256_file(p) for p in (CYCLEFIX_RAW, REPEAT_RAW, FRESH_RAW, HP2_CATS, FRESH_CATS)]), 'meta|diag_utc': np.array([utc()])}
    for d in L:
        t0 = time.time()
        arrays.update(diag_patch(d))
        print('diag %s %s [%.0f s]' % (d['label'], d['pxy'], time.time() - t0), flush=True)
    _save_flat(out_path, arrays)


def _save_flat(path, arrays):
    """Adds `arrays` to the raw `path` (keys already complete); atomic; an existing key is never overwritten by a different value."""
    old = {}
    if os.path.exists(path):
        with np.load(path, allow_pickle=False) as z:
            old = {k: z[k] for k in z.files}
    for k, v in arrays.items():
        if k in old and not np.array_equal(old[k], v):
            raise SystemExit('key %s exists with another value: the raw is write-once' % k)
    old.update(arrays)
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    tmp = path + '.tmp.npz'
    np.savez_compressed(tmp, **old)
    os.replace(tmp, path)


# ----------------------------------------------------------------------------------------------------------------------------------
# held-out test: registry, baselines, ISS runs
# ----------------------------------------------------------------------------------------------------------------------------------
def registry():
    return F.load_registry(REGISTRY)


def heldout_patches(reg=None):
    reg = reg or registry()
    return [(c, tuple(p)) for c in sorted(reg) for p in reg[c]['patches']]


def hkey(cluster, pxy):
    return 'H|%s|x%d_y%d|' % (cluster, pxy[0], pxy[1])


def baseline_rows(b, cluster, pxy):
    """{label: (x, y, flux)}: the held-out comparators from the catalogue file `b`: crowdsource and photutils as tuned on the M2 development patches (freeze.json) and with the ISS's front end (round 2)."""
    meta = {k: str(b['meta|' + k][0]) for k in ('photutils_tuned', 'crowdsource_tuned', 'photutils_fe', 'crowdsource_fe')}
    return FR.baseline_dets(b, cluster, pxy, dict(tuned_phot=meta['photutils_tuned'], cs_tuned=meta['crowdsource_tuned'], fe_phot=meta['photutils_fe'], fe_cs=meta['crowdsource_fe']))


def heldout_arrays(args):
    """One held-out patch: the frozen ISS (isotropic kernel; kernel rule where the stamp is not adequately isotropic) through the guarded path at the tuned configuration, the repair, the baselines, each scored
    against the registration-corrected truth ('reg|', primary) and the truth as it is ('' , secondary). A function of its arguments only."""
    cluster, pxy, entry, shift, base = args
    patch = F.load_patch(pxy[0], pxy[1], entry)
    truth = F.truth_arrays(entry)
    t = A.truth_for_patch(truth, pxy[0], pxy[1])
    t_reg = A.truth_for_patch(dict(px=truth['px'] + shift[0], py=truth['py'] + shift[1], m606=truth['m606']), pxy[0], pxy[1])
    truths = {'': t, 'reg|': t_reg}
    fit = PG.fit_adequate(patch['psf'])
    arr = dict(kind=np.array([fit['kind']]), shift=np.array([shift[0], shift[1]]), fwhm=np.array([S2.psf_fwhm(patch['psf'])]), rms=np.array([fit['rms_isotropic'], fit['rms']]), shape=PS.stamp_shape(patch['psf']))
    CF.install()
    try:
        (dets, secs, n, info), tr = FR._guarded_run(X.run_path, patch, HP, REX)
        FR._store(arr, 'iso', dets, secs, n, info, truths, KAP, tr)
        if fit['kind'] == 'isotropic':
            FR._alias(arr, 'iso', 'rule')
        else:
            (dets2, secs2, n2, info2), tr2 = FR._guarded_run(PS.run_path_psf, patch, HP, REX, fit)
            FR._store(arr, 'rule', dets2, secs2, n2, info2, truths, KAP, tr2)
    finally:
        CF.uninstall()
    for r in ROWS:                                           # the tuned readout (kappa 2) also under the bare row name, for the tables
        for tname in truths:
            for f in ('n_true', 'n_matched', 'n_det', 'n_det_matched', 'completeness', 'fp_fraction', 'offset', 'blk22', 'tot'):
                arr['%s%s|%s' % (tname, r, f)] = arr['%s%s|k%g|%s' % (tname, r, KAP, f)]
    for r in ROWS:
        kfit = PG.fit_adequate(patch['psf'], rms_goal=np.inf) if r == 'iso' else fit
        off, tab = C.coherence_table(patch, kfit)
        det = arr[r + '|det']
        for rho in RHO_GRID:
            m, nm = C.merge_coherent(det, off, tab, rho)
            for tname, tt in truths.items():
                arr.update(X._pack('%s%s_m%g' % (tname, r, rho), score_det(m, tt)))
            arr['%s_m%g|merges' % (r, rho)] = np.array([nm])
            if rho == RHO_STAR:
                arr['%s_m|det' % r] = m
    b = np.load(base)
    for label, dd in baseline_rows(b, cluster, pxy).items():
        for tname, tt in truths.items():
            arr.update(X._pack(tname + label, score_det(np.c_[dd[0], dd[1], dd[2]], tt)))
        arr[label + '|n'] = np.array([len(dd[0])])
    return cluster, pxy, arr


def frozen_code_hash():
    """SHA-256 of the source of the functions and the constants that decide a held-out number (the run, the baseline rows, the scoring, the keys): a later edit of the report does not change it."""
    import inspect
    parts = [inspect.getsource(f) for f in (heldout_arrays, baseline_rows, part_run, score_det, hkey, heldout_patches, registry)]
    parts.append(repr((RHO_GRID, RHO_STAR, HP, REX, KAP, ROWS)))
    return hashlib.sha256('\n'.join(parts).encode()).hexdigest()


def freeze_record():
    rec = {}
    for f in FILES:
        p = os.path.join(REPO, f) if not os.path.isabs(f) else f
        if f.endswith('clusterdep_plan.md'):
            p = PLAN
        elif f.endswith('clusterdep_patches.json'):
            p = REGISTRY
        elif f.endswith('clusterdep_baseline_catalogs.npz'):
            p = BASE
        rec[f] = sha256_file(p) if os.path.exists(p) else None
    return dict(files=rec, code=frozen_code_hash(), rho_star=RHO_STAR, config=dict(sigma_b=HP, r_ex=REX, kappa=KAP), utc=utc())


def write_freeze(path=None):
    path = path or FREEZE
    if os.path.exists(path):
        raise SystemExit('%s exists: the freeze record is written once' % path)
    rec = freeze_record()
    if any(v is None for v in rec['files'].values()):
        raise SystemExit('missing files: %s' % [k for k, v in rec['files'].items() if v is None])
    with open(path, 'w') as f:
        json.dump(rec, f, indent=1, sort_keys=True)
        f.write('\n')
    print('wrote', path)


def check_freeze(path=None):
    path = path or FREEZE
    if not os.path.exists(path):
        raise SystemExit('no freeze record: run --freeze first')
    with open(path) as f:
        rec = json.load(f)
    now = freeze_record()
    bad = [k for k in now['files'] if now['files'][k] != rec['files'].get(k)]
    if bad or rec['rho_star'] != RHO_STAR or rec['config'] != now['config'] or rec.get('code') != now['code']:
        raise SystemExit('the freeze record does not match: %s' % (bad or 'rho_star, configuration or code'))
    return rec


def part_run(out_path=RAW, workers=1):
    check_freeze()
    reg = registry()
    shifts = {(c, tuple(p)): reg[c]['shifts'][i] for c in reg for i, p in enumerate(reg[c]['patches'])}
    todo = [(c, p) for c, p in heldout_patches(reg) if (hkey(c, p) + 'kind') not in (np.load(out_path).files if os.path.exists(out_path) else [])]
    for c, p in todo:
        t0 = time.time()
        _, _, arr = heldout_arrays((c, p, reg[c], shifts[(c, p)], BASE))
        _save_flat(out_path, {hkey(c, p) + k: v for k, v in arr.items()} | {'meta|run_utc|%s' % hkey(c, p): np.array([utc()])})
        print('H %s %s done [%.0f s]' % (c, p, time.time() - t0), flush=True)


R1_CASES = (('M92', (1100, 220), 'rule'), ('M92', (1100, 220), 'iso'), ('M92', (1430, 550), 'rule'), ('M53', (640, 680), 'rule'), ('M53', (660, 420), 'iso'), ('M2', (310, 630), 'iso'), ('M2r2', (133, 281), 'rule'))
R1_KAPPAS = (2.0, 3.0, 5.0)


def r1_key(label, pxy, row):
    return 'R1|%s|x%d_y%d|%s|' % (label, pxy[0], pxy[1], row)


def part_r1(out_path=RAW):
    """POST HOC development run of the first attempt at a repair (`astro_clusterdep.calibrate`): the weights 1 / sqrt(s0^2 var0 + eta^2 qd^2) fitted to the stored first-pass catalogue of the row, then the guarded path
    with them, read out at kappa 2, 3, 5. Seven spent cases; a case whose result exists is not run again."""
    have = set(np.load(out_path).files) if os.path.exists(out_path) else set()
    todo = [c for c in R1_CASES if r1_key(*c) + 'prm' not in have]
    if not todo:
        return
    L = {(d['label'], d['pxy']): d for d in spent_patches()}
    for label, pxy, row in todo:
        t0 = time.time()
        d = L[(label, pxy)]
        pa = d['patch']
        fit = PG.fit_adequate(pa['psf']) if row == 'rule' else None
        P = C.make_kernel(pa, HP, fit)
        sw, prm = C.calibrated_weights(pa, P, d[row])
        CF.install()
        try:
            (out, secs, n, info), tr = FR._guarded_run(C.run_path_weighted, pa, HP, REX, sw, fit)
        finally:
            CF.uninstall()
        k = r1_key(label, pxy, row)
        arr = {k + 'prm': np.array([prm['s0'], prm['eta']]), k + 'secs': np.array([secs, n]), k + 'counts': tr['counts']}
        for kap in R1_KAPPAS:
            arr.update(X._pack('%sk%g' % (k, kap), score_det(np.c_[out[kap][0], out[kap][1], out[kap][2]], d['truth'])))
        arr[k + 'split'] = bright_split(np.c_[out[KAP][0], out[KAP][1], out[KAP][2]], d['truth'])
        _save_flat(out_path, arr)
        print('R1 %s %s %s: s0 %.2f eta %.3f, %d entries [%.0f s]' % (label, pxy, row, prm['s0'], prm['eta'], n, time.time() - t0), flush=True)


# ----------------------------------------------------------------------------------------------------------------------------------
# statistics of the diagnosis
# ----------------------------------------------------------------------------------------------------------------------------------
def spearman(a, b):
    """Spearman rank correlation (average ranks for ties); nan for a constant input."""
    from scipy.stats import rankdata
    a, b = rankdata(a), rankdata(b)
    if a.std() == 0 or b.std() == 0:
        return float('nan')
    return float(np.corrcoef(a, b)[0, 1])


def loco(x, y, groups):
    """Leave-one-cluster-out prediction of y from the single covariate x: for each group, ordinary least squares on the patches of the other groups, prediction for the held-out patches. Returns (rmse of the
    prediction, rmse of the null predictor = the mean of the training y, Spearman correlation of prediction and truth, number of groups whose mean prediction has the sign of their mean y)."""
    x, y, groups = np.asarray(x, float), np.asarray(y, float), np.asarray(groups)
    pred, null = np.zeros_like(y), np.zeros_like(y)
    for g in np.unique(groups):
        tr, te = groups != g, groups == g
        xs = x[tr]
        if xs.std() > 0:
            b = np.cov(xs, y[tr], bias=True)[0, 1] / xs.var()
            a = y[tr].mean() - b * xs.mean()
        else:
            a, b = y[tr].mean(), 0.0
        pred[te], null[te] = a + b * x[te], y[tr].mean()
    sign = sum(1 for g in np.unique(groups) if np.sign(pred[groups == g].mean()) == np.sign(y[groups == g].mean()))
    return float(np.sqrt(np.mean((pred - y) ** 2))), float(np.sqrt(np.mean((null - y) ** 2))), spearman(pred, y), int(sign)


def spent_index(raw):
    """[(label, pxy)] of the spent patches in the raw, in stored order."""
    out = []
    for s in raw['meta|labels']:
        lab, xy = str(s).split('|')
        x, y = xy[1:].split('_y')
        out.append((lab, (int(x), int(y))))
    return out


def tot_f(raw, key):
    tot = raw[key + '|tot']
    return X.f_measure(tot)


def blocks_of(raw, keys):
    """The blocks (nb, 4) of the patches' rows concatenated and the strata sizes."""
    parts = [raw[k + '|blk22'] for k in keys]
    return np.concatenate(parts), [len(p) for p in parts]


# ----------------------------------------------------------------------------------------------------------------------------------
# the report
# ----------------------------------------------------------------------------------------------------------------------------------
def fmt(v, nd=3):
    return 'nan' if v is None or not np.isfinite(v) else '%.*f' % (nd, v)


def group_blocks(raw, idxs, row_fn):
    """Concatenated blocks of the rows row_fn(i) of the patches idxs, and the strata sizes."""
    return blocks_of(raw, [row_fn(i) for i in idxs])


def gap_row(raw, pk, a, b):
    """dF of the rows a and b of one spent patch (point value only)."""
    return tot_f(raw, pk + a)[0] - tot_f(raw, pk + b)[0]


def diag_section(raw):
    idx = spent_index(raw)
    pks = [pkey(l, p) for l, p in idx]
    labels = [l for l, p in idx]
    clusters = np.array([CLUSTER_OF[l] for l in labels])
    L = ['## 1. Diagnosis on the spent patches (POST HOC: nothing here is a test)', '',
         'The 19 test patches of the earlier rounds, every number from the stored catalogues (guarded ISS at the tuned configuration, kappa 2; front-end crowdsource `cs2_hp2_t3_pm`, the comparator of the fresh test), none re-run. '
         '"ISS iso" is the isotropic kernel as frozen, "ISS rule" the kernel that `psf_general.fit_adequate` chooses from the stamp (the same program as "iso" where the stamp is adequately round: M2 and its repeat frames 1 and 3, M3). '
         'The unit of replication of the dependence is the cluster: M2 (4 windows + 6 repeat patches), M53, M92, M3, M15 = 5 clusters.', '']
    rows = []
    for (l, p), pk in zip(idx, pks):
        fi, fr_, fc = tot_f(raw, pk + 'iso')[0], tot_f(raw, pk + 'rule')[0], tot_f(raw, pk + 'cs')[0]
        rows.append([l, '(%d, %d)' % p, str(raw[pk + 'kind'][0]), fmt(raw[pk + 'cov|n22'][0], 0), fmt(fc), fmt(fi), fmt(fr_), fmt(fi - fc), fmt(fr_ - fc), fmt(raw[pk + 'cov|rsd_iso'][0], 2), fmt(raw[pk + 'cov|rsd_rule'][0], 2)])
    L += [X.table(['frame', 'patch', 'kernel rule', 'truth < 22', 'F crowdsource (front end)', 'F ISS iso', 'F ISS rule', 'ISS iso - cs', 'ISS rule - cs', 'residual scale iso', 'residual scale rule'], rows), '']
    # per-cluster paired intervals
    L += ['### 1.1 The gap by cluster (patches pooled, 20 x 20 px block bootstrap, 2000 resamples, seed 0, stratified within patches)', '']
    rowsc = []
    for c in ['M2', 'M53', 'M92', 'M3', 'M15']:
        ii = [i for i in range(len(idx)) if clusters[i] == c]
        out = [c, str(len(ii))]
        for a in ('iso', 'rule'):
            ba, sizes = blocks_of(raw, [pks[i] + a for i in ii])
            bb, _ = blocks_of(raw, [pks[i] + 'cs' for i in ii])
            df = X.boot_diff(ba, bb, X.boot_idx_strat(sizes))[2]
            out.append('%s %s' % (X.fmt_iv(df[0], df[1]), X.verdict(*df)))
        rowsc.append(out)
    L += [X.table(['cluster', 'patches', 'ISS iso minus crowdsource (front end)', 'ISS rule minus crowdsource (front end)'], rowsc), '']
    # covariates
    L += ['### 1.2 Which covariate predicts the gap? (patch level, n = 19; rank correlation, and leave-one-cluster-out linear prediction over the 5 clusters)', '',
          'Gap = F(ISS) - F(crowdsource, front end) of the patch. Spearman rho on the 19 patches (the patches of one cluster are not independent: read it as description). LOCO: ordinary least squares on the patches of the 4 other clusters, '
          'prediction for the held-out cluster\'s patches; rmse of the prediction against the rmse of the null predictor (the training mean), over all 19 predictions; "signs" = clusters whose mean predicted gap has the sign of their mean gap (of 5). '
          'A covariate that really explains the dependence has LOCO rmse below the null rmse and 5 / 5 signs. NOTE (erratum to `clusterdep_plan.md` section 1, which was committed before this table was read closely): the plan says that sky, noise, crowding '
          'and background structure are at or above the null in LOCO rmse and wrong in the signs; the lines under the tables below list the covariates that are below the null in rmse, and some of those (sky, noise, density) are, with 3 / 5 signs or fewer, '
          'confounded with the cluster (the M2 frames are brighter and denser than the others); the plan\'s sentence is superseded by those lines.', '']
    for a, name in (('iso', 'ISS iso'), ('rule', 'ISS rule')):
        gap = np.array([gap_row(raw, pk, a, 'cs') for pk in pks])
        rowsx, better = [], []
        for k in COVARIATES:
            x = np.array([raw[pk + 'cov|' + k][0] for pk in pks])
            rho = spearman(x, gap)
            rm, rn, rs, sg = loco(x, gap, clusters)
            rowsx.append([k, COV_NOTE[k], fmt(rho, 2), fmt(rm, 3), fmt(rn, 3), fmt(rs, 2), '%d / 5' % sg])
            if rm < rn:
                better.append((k, rm, rn, sg))
        L += ['**%s minus crowdsource**' % name, '', X.table(['covariate', 'meaning', 'Spearman rho', 'LOCO rmse', 'null rmse', 'LOCO rho', 'signs'], rowsx), '']
        L += ['Covariates whose LOCO rmse is below the null rmse (printed by the script; best first): ' + ('; '.join('%s (%.3f against %.3f, signs %d / 5)' % (k, rm, rn, sg) for k, rm, rn, sg in sorted(better, key=lambda b: b[1])) if better else 'none') + '.', '']
    L += roundness_paragraph(raw, pks, clusters)
    # recall / purity by magnitude, pooled by cluster
    L += ['### 1.3 Completeness and purity by magnitude (pooled counts by cluster; ISS rule row and front-end crowdsource)', '']
    edges = A.BIN_EDGES
    bands = [(15.0, 18.0), (18.0, 19.0), (19.0, 20.0), (20.0, 21.0), (21.0, 22.0)]
    rowsm = []
    for c in ['M2', 'M53', 'M92', 'M3', 'M15']:
        ii = [i for i in range(len(idx)) if clusters[i] == c]
        for a, name in (('rule', 'ISS rule'), ('cs', 'crowdsource')):
            bins = sum(raw[pks[i] + a + '|bins'] for i in ii)
            cells = []
            for lo, hi in bands:
                sel = (edges[:-1] >= lo - 1e-9) & (edges[1:] <= hi + 1e-9)
                nt, nm, nd, ndm = (bins[sel, j].sum() for j in range(4))
                cells.append('%s / %s' % (fmt(nm / max(nt, 1), 2), fmt(ndm / max(nd, 1), 2)))
            rowsm.append([c, name] + cells)
    L += [X.table(['cluster', 'catalogue'] + ['%g-%g mag: recall / purity' % b for b in bands], rowsm), '']
    # error decomposition
    L += ['### 1.4 Where the ISS loses: decomposition of the errors at mag < 22 (pooled counts by cluster)', '',
          'Truth stars unmatched: a detection within 0.5 px whose magnitude is off by more than 0.5 (blend flux), the nearest detection at 0.5 to 1.5 px (position or split), none within 1.5 px (missed). '
          'Detections unmatched: a truth star of any magnitude within 0.5 px (magnitude off), at 0.5 to 1.5 px, none within 1.5 px (spurious). Bright-star split: the scored truth stars brighter than 19 mag: mean number of '
          'detections within 1.5 px, fraction with a detection within 0.5 px, mean distance of the nearest detection (px), fraction matched.', '']
    rowsd = []
    for c in ['M2', 'M53', 'M92', 'M3', 'M15']:
        ii = [i for i in range(len(idx)) if clusters[i] == c]
        for a, name in (('iso', 'ISS iso'), ('rule', 'ISS rule'), ('cs', 'crowdsource')):
            dec = sum(raw[pks[i] + a + '|dec'] for i in ii)
            spl = np.array([raw[pks[i] + a + '|split'] for i in ii])
            w = spl[:, 0]
            sp = (spl[:, 1:] * w[:, None]).sum(0) / max(w.sum(), 1)
            nt = max(dec[0], 1)
            rowsd.append([c, name, '%d' % dec[0], '%.2f' % (dec[1] / nt), '%.2f' % (dec[2] / nt), '%.2f' % (dec[3] / nt), '%.2f' % (dec[4] / nt), '%d' % dec[5], '%.2f' % (dec[6] / max(dec[5], 1)),
                          '%.2f / %.2f / %.2f' % tuple(dec[7:10] / max(dec[6], 1)), '%.2f, %.2f, %.2f, %.2f' % tuple(sp)])
    L += [X.table(['cluster', 'catalogue', 'truth < 22', 'matched', 'mag off (<= 0.5 px)', 'nearest 0.5-1.5 px', 'missed', 'detections < 22', 'unmatched', 'unmatched: <= 0.5 px / 0.5-1.5 px / none', 'bright split: dets within 1.5 px, dets within 0.5 px, nearest px, matched'], rowsd), '']
    return L, idx, pks, clusters


CLUSTERS5 = ['M2', 'M53', 'M92', 'M3', 'M15']


def loco_means(x, y):
    """Leave-one-cluster-out linear prediction at the level of the cluster means (n points, one per cluster): (rmse of the prediction, rmse of the null predictor = the mean of the other clusters, number of clusters whose
    prediction has the sign of their gap, Spearman correlation of x and y)."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    n = len(y)
    pred, null = np.zeros(n), np.zeros(n)
    for g in range(n):
        tr = np.arange(n) != g
        xs = x[tr]
        b = np.cov(xs, y[tr], bias=True)[0, 1] / xs.var() if xs.std() > 0 else 0.0
        a = y[tr].mean() - b * xs.mean()
        pred[g], null[g] = a + b * x[g], y[tr].mean()
    return float(np.sqrt(np.mean((pred - y) ** 2))), float(np.sqrt(np.mean((null - y) ** 2))), int(np.sum(np.sign(pred) == np.sign(y))), spearman(x, y)


def cluster_permutation(raw, pks, clusters, a):
    """The diagnosis at the level of the 5 clusters (the covariates are constant within a frame, so the 19 patches are pseudo-replicates): the cluster means of the gap F(row a) - F(front-end crowdsource) and of every covariate;
    the observed best covariate (smallest leave-one-cluster-out rmse of the cluster means) and the permutation p-value of that best: over all 5! = 120 assignments of the 5 gaps to the 5 clusters, the share in which the best of the
    16 covariates does at least as well (rmse not above the observed best). Also the value for the stamp ellipticity alone. Returns None unless the 5 clusters are present."""
    import itertools
    cl = [c for c in CLUSTERS5 if c in set(clusters)]
    if len(cl) != 5:
        return None
    gap = np.array([gap_row(raw, pk, a, 'cs') for pk in pks])
    y = np.array([gap[clusters == c].mean() for c in cl])
    X_ = {k: np.array([np.mean([raw[pk + 'cov|' + k][0] for pk, c_ in zip(pks, clusters) if c_ == c]) for c in cl]) for k in COVARIATES}
    obs = {k: loco_means(X_[k], y) for k in COVARIATES}
    best = min(obs, key=lambda k: (obs[k][0], k))
    n_best = n_ell = 0
    perms = list(itertools.permutations(range(5)))
    for pm in perms:
        yp = y[list(pm)]
        n_best += int(min(loco_means(X_[k], yp)[0] for k in COVARIATES) <= obs[best][0] + 1e-12)
        n_ell += int(loco_means(X_['stamp_ell'], yp)[0] <= obs['stamp_ell'][0] + 1e-12)
    return dict(best=best, obs=obs, n_best=n_best, n_ell=n_ell, n_perm=len(perms), gap=y)


def roundness_paragraph(raw, pks, clusters):
    """The statement of what the covariate diagnosis can bear, with the cluster-level numbers and the permutation p-values (iso and rule gaps)."""
    res = {a: cluster_permutation(raw, pks, clusters, a) for a in ('iso', 'rule')}
    if res['iso'] is None:
        return []
    L = ['### 1.2b How much the covariate diagnosis can bear (POST HOC)', '',
         'The roundness of the PSF stamp is a hypothesis at n = 5 clusters with 16 covariates, not a finding. The covariates are constant within a frame and the patches of a frame share them, so the 19 patches are pseudo-replicates of 5 clusters; the cluster-level '
         'numbers below use the 5 cluster means (gap and covariate), leave-one-cluster-out as before, and a permutation test of the cluster gaps: all %d assignments of the 5 gaps to the 5 clusters, and the share of them in which the BEST of the 16 covariates '
         'predicts at least as well as the observed best (the null that none of the covariates relates to the gap, with the search over 16 covariates included).' % res['iso']['n_perm'], '']
    for a, name in (('iso', 'ISS iso'), ('rule', 'ISS rule')):
        r = res[a]
        e = r['obs']['stamp_ell']
        b = r['obs'][r['best']]
        L.append('* %s minus crowdsource, cluster means: stamp ellipticity LOCO rmse %.3f against the null %.3f, signs %d / 5, Spearman %.2f (n = 5); the best of the 16 covariates is %s (rmse %.3f against %.3f, signs %d / 5). Permutation p-value that some '
                 'covariate does as well as that best: %d / %d = %.2f; for the stamp ellipticity alone %d / %d = %.2f.' % (name, e[0], e[1], e[2], e[3], r['best'], b[0], b[1], b[2], r['n_best'], r['n_perm'], r['n_best'] / r['n_perm'], r['n_ell'], r['n_perm'],
                                                                                                                     r['n_ell'] / r['n_perm']))
    L += ['', 'So: the stamp roundness is a candidate explanation, one hypothesis among several, from 5 clusters; the best covariate for the isotropic gap is itself an output of the same ISS fit (the residual scale), and the sky, noise and density covariates '
          'also beat the null (the erratum above). The report does not claim that roundness predicts or explains the dependence.', '']
    return L


def merge_section(raw, idx, pks, clusters):
    L = ['## 2. The repair examined on the spent patches (POST HOC, development): merge of units whose kernels have coherence >= rho', '',
         'The repair acts on the finished catalogue: units whose kernels have normalised inner product at least rho (`astro_clusterdep.merge_coherent`: closest pair first, flux-weighted centroid, summed flux) are one source. The table gives, per row and '
         'cluster, the mean over the cluster\'s patches of F at mag < 22 (the unmerged row first). rho* = %g is the value of the grid with the largest mean over the 19 patches of the ISS rule row; it is the value declared in `clusterdep_plan.md`.' % RHO_STAR, '']
    for a, name in (('rule', 'ISS rule'), ('iso', 'ISS iso')):
        rows = []
        for c in ['M2', 'M53', 'M92', 'M3', 'M15', 'all 19 patches']:
            ii = [i for i in range(len(idx)) if c == 'all 19 patches' or clusters[i] == c]
            cells = [fmt(np.mean([tot_f(raw, pks[i] + a)[0] for i in ii]))]
            for rho in RHO_GRID:
                cells.append(fmt(np.mean([tot_f(raw, '%s%s_m%g' % (pks[i], a, rho))[0] for i in ii])))
            cells.append(fmt(np.mean([tot_f(raw, pks[i] + 'cs')[0] for i in ii])))
            rows.append([c, str(len(ii))] + cells)
        L += ['**%s**' % name, '', X.table(['cluster', 'patches', 'unmerged'] + ['rho = %g' % r for r in RHO_GRID] + ['crowdsource (front end)'], rows), '']
    best = max(RHO_GRID, key=lambda r: (round(np.mean([tot_f(raw, '%s%s_m%g' % (pk, 'rule', r))[0] for pk in pks]), 6), -r))
    L += ['The largest mean over the 19 patches of the ISS rule row: rho = %g (declared rho* = %g: %s).' % (best, RHO_STAR, 'the same' if best == RHO_STAR else 'DIFFERENT, see the plan'), '']
    # merged gap by cluster
    L += ['**Gap to the front-end crowdsource after the merge at rho* (patches pooled, paired block bootstrap)**', '']
    rowsc = []
    for c in ['M2', 'M53', 'M92', 'M3', 'M15']:
        ii = [i for i in range(len(idx)) if clusters[i] == c]
        out = [c]
        for a in ('rule', 'rule_m%g' % RHO_STAR):
            ba, sizes = blocks_of(raw, [pks[i] + a for i in ii])
            bb, _ = blocks_of(raw, [pks[i] + 'cs' for i in ii])
            df = X.boot_diff(ba, bb, X.boot_idx_strat(sizes))[2]
            out.append('%s %s' % (X.fmt_iv(df[0], df[1]), X.verdict(*df)))
        rowsc.append(out)
    L += [X.table(['cluster', 'ISS rule minus crowdsource', 'ISS rule merged at rho* minus crowdsource'], rowsc), '',
          '(Development numbers: rho* was chosen on these patches, so these intervals are not a test. The test is section 3.)', '']
    return L


def r1_section(raw):
    cases = [c for c in R1_CASES if r1_key(*c) + 'prm' in raw]
    if not cases:
        return []
    L = ['## 2b. The first attempt at a repair, examined on seven spent cases (POST HOC, development): a noise model with a model-error floor', '',
         'Weights 1 / sqrt(s0^2 var0 + eta^2 qd^2): var0 the frozen variance, qd the smoothed absolute value of the high-passed model of the first-pass catalogue, (s0 >= 1, eta >= 0) fitted to that catalogue\'s own residual '
         '(`astro_clusterdep.calibrate`; label-free); then the whole guarded path is run again with these weights (`run_path_weighted`: the frozen program with other weights), read out at kappa 2 (the tuned value), 3 and 5. '
         'It is the principled form of "the residual scale of the round-kernel catalogue predicts the gap". It is NOT carried into the held-out test: it was dropped on the two M92 patches (F moved by 0.00 to 0.02, the number of units around bright stars did not fall) before the plan was written; the other five cases were run afterwards, with the results below (a gain of 0.04 on M53, a loss of 0.02 on the M2 control, little elsewhere). It is a lead that has no held-out evidence.', '']
    rows = []
    for c in cases:
        label, pxy, row = c
        pk = pkey(label, pxy)
        k = r1_key(*c)
        s0, eta = raw[k + 'prm']
        fz = tot_f(raw, pk + row)[0]
        cells = [fmt(tot_f(raw, '%sk%g' % (k, kap))[0]) for kap in R1_KAPPAS]
        rows.append([label, '(%d, %d)' % pxy, row, '%.2f' % s0, '%.3f' % eta, fmt(fz)] + cells + [fmt(tot_f(raw, pk + 'cs')[0]), '%.2f -> %.2f' % (raw[pk + row + '|split'][1], raw[k + 'split'][1]), '%d' % raw[k + 'counts'][6]])
    L += [X.table(['frame', 'patch', 'row', 's0', 'eta', 'F frozen', 'F calibrated kappa 2', 'kappa 3', 'kappa 5', 'F crowdsource (front end)', 'dets within 1.5 px of a bright star, frozen -> calibrated', 'guard activations'], rows), '']
    return L


def heldout_index(raw):
    out = {}
    for k in raw:
        if k.startswith('H|') and k.endswith('|kind'):
            _, c, xy = k.split('|')[:3]
            x, y = xy[1:].split('_y')
            out.setdefault(c, []).append((int(x), int(y)))
    return {c: sorted(v) for c, v in sorted(out.items())}


HELD_ROWS = [('ISS frozen (kernel rule)', 'rule'), ('ISS repaired (kernel rule, merge at rho*)', 'rule_m%g' % RHO_STAR), ('ISS frozen (isotropic kernel)', 'iso'), ('ISS repaired (isotropic kernel, merge at rho*)', 'iso_m%g' % RHO_STAR),
             ('crowdsource (front end)', 'crowdsource (front end)'), ('crowdsource (tuned)', 'crowdsource (tuned)'), ('photutils (front end)', 'photutils (front end)'), ('photutils (tuned)', 'photutils (tuned)')]
PRIMARY = [('O2', 'ISS repaired minus crowdsource (front end)', 'rule_m%g' % RHO_STAR, 'crowdsource (front end)'),
           ('O2f', 'ISS frozen minus crowdsource (front end)', 'rule', 'crowdsource (front end)'),
           ('O1', 'ISS repaired minus crowdsource (tuned)', 'rule_m%g' % RHO_STAR, 'crowdsource (tuned)'),
           ('O1f', 'ISS frozen minus crowdsource (tuned)', 'rule', 'crowdsource (tuned)'),
           ('O3', 'ISS repaired minus photutils (front end)', 'rule_m%g' % RHO_STAR, 'photutils (front end)'),
           ('O3f', 'ISS frozen minus photutils (front end)', 'rule', 'photutils (front end)'),
           ('R', 'ISS repaired minus ISS frozen', 'rule_m%g' % RHO_STAR, 'rule')]


def hblocks(raw, c, plist, key, tname):
    return blocks_of(raw, [hkey(c, p) + tname + key for p in plist])


def verdicts_heldout(raw, tname):
    """{(id, cluster): (dF, interval, verdict)} per held-out cluster (patches pooled, stratified bootstrap, seed 0)."""
    hi = heldout_index(raw)
    out = {}
    for c, pl in hi.items():
        for cid, name, a, b in PRIMARY:
            ba, sizes = hblocks(raw, c, pl, a, tname)
            bb, _ = hblocks(raw, c, pl, b, tname)
            df = X.boot_diff(ba, bb, X.boot_idx_strat(sizes))[2]
            out[(cid, c)] = (df[0], df[1], X.verdict(*df))
    return out


def pooled_heldout(raw, tname, a, b):
    """Pooled over every held-out patch: dF(a - b) with the stratified block bootstrap, and a leave-one-cluster-out view: the mean over clusters of the cluster's dF."""
    hi = heldout_index(raw)
    allp = [(c, p) for c, pl in hi.items() for p in pl]
    ba = np.concatenate([raw[hkey(c, p) + tname + a + '|blk22'] for c, p in allp])
    bb = np.concatenate([raw[hkey(c, p) + tname + b + '|blk22'] for c, p in allp])
    sizes = [len(raw[hkey(c, p) + tname + a + '|blk22']) for c, p in allp]
    return X.boot_diff(ba, bb, X.boot_idx_strat(sizes))[2]


SEED_SHARE_N = 30                                    # bootstrap seeds 0 ... 29 of the seed-sensitivity check


def failed_registration(raw):
    """{(cluster, patch)} of the held-out patches whose registration fit failed, by the rule stated in section 3.3 of the report: the corrected scoring is worse than HALF the scoring as it is for both the frozen ISS
    (kernel rule) and the front-end crowdsource (a failed fit of the offset, not a property of the methods). Post hoc: the plan had no registration gate."""
    out = set()
    for c, pl in heldout_index(raw).items():
        for p in pl:
            k = hkey(c, p)
            f = {(r, t): tot_f(raw, k + t + r)[0] for r in ('rule', 'crowdsource (front end)') for t in ('reg|', '')}
            if all(f[(r, 'reg|')] < 0.5 * f[(r, '')] for r in ('rule', 'crowdsource (front end)')):
                out.add((c, p))
    return out


def verdicts_subset(raw, tname, skip=(), seed=0):
    """{(id, cluster): (dF, interval, verdict)} as `verdicts_heldout`, but without the patches in `skip` ({(cluster, patch)}) and with the bootstrap seed `seed`; a cluster left without a patch is absent."""
    out = {}
    for c, pl in heldout_index(raw).items():
        pl = [p for p in pl if (c, p) not in skip]
        if not pl:
            continue
        for cid, name, a, b in PRIMARY:
            ba, sizes = hblocks(raw, c, pl, a, tname)
            bb, _ = hblocks(raw, c, pl, b, tname)
            df = X.boot_diff(ba, bb, X.boot_idx_strat(sizes, seed=seed))[2]
            out[(cid, c)] = (df[0], df[1], X.verdict(*df))
    return out


def count3(verdicts):
    """'w won / t tied / l lost' of a list of verdict words."""
    return '%d won / %d tied / %d lost' % (verdicts.count('win'), verdicts.count('tie'), verdicts.count('loss'))


def side_by_side_section(raw):
    """3.0: the two scorings side by side (the as-is scoring is the headline robustness check), the clusters without a failed registration, and the seed sensitivity of the cluster verdicts."""
    hi = heldout_index(raw)
    clist = list(hi)
    bad = failed_registration(raw)
    L = ['### 3.0 Headline robustness check: the two scorings side by side, the failed registrations, the bootstrap seeds', '',
         'POST HOC (the plan had no registration gate; the registration fit is described in section 3.3). The cluster verdicts of the frozen ISS (O2f) and of the repaired ISS (O2) against the front-end crowdsource under the PRIMARY scoring '
         '(registration-corrected truth) and under the SECONDARY scoring (the truth as it is), side by side. The as-is scoring is the robustness check to read first.', '']
    cvp, cva = verdicts_heldout(raw, 'reg|'), verdicts_heldout(raw, '')
    rows = []
    for cid, name, a, b in PRIMARY[:2]:
        for lab, cv in (('primary (corrected truth)', cvp), ('secondary (truth as is)', cva)):
            ver = [cv[(cid, c)][2] for c in clist]
            rows.append([cid, lab] + ['%s %s' % (X.fmt_iv(cv[(cid, c)][0], cv[(cid, c)][1]), cv[(cid, c)][2]) for c in clist] + [count3(ver)])
    L += [X.table(['outcome', 'scoring'] + ['%s (%d patch%s)' % (c, len(hi[c]), '' if len(hi[c]) == 1 else 'es') for c in clist] + ['clusters won / tied / lost'], rows), '']
    n_p = sum(len(v) for v in hi.values())
    if bad:
        left = [c for c in clist if any((c, p) not in bad for p in hi[c])]
        L += ['**Failed registration fits.** %d of %d patches: %s. Rule (post hoc): the corrected scoring is worse than half of the scoring as it is for both the frozen ISS and the front-end crowdsource (section 3.3); the fitted offset is taken to be wrong there, so the corrected '
              'scoring is treated as void there. %d of %d clusters carry information under the corrected scoring (%s); %s.' % (
                  len(bad), n_p, ', '.join('%s (%d, %d)' % ((c,) + p) for c, p in sorted(bad)), len(left), len(clist), ', '.join(left),
                  'clusters without a usable patch there: ' + ', '.join(c for c in clist if c not in left) if len(left) < len(clist) else 'every cluster keeps at least one patch'), '']
        rows = []
        for cid, name, a, b in PRIMARY[:2]:
            for lab, tname in (('primary (corrected truth), failed patches excluded', 'reg|'), ('secondary (truth as is), the same patches excluded', '')):
                cv = verdicts_subset(raw, tname, bad)
                ver = [cv[(cid, c)][2] for c in left]
                rows.append([cid, lab] + ['%s %s' % (X.fmt_iv(cv[(cid, c)][0], cv[(cid, c)][1]), cv[(cid, c)][2]) if (cid, c) in cv else 'no information' for c in clist] + [count3(ver)])
        L += [X.table(['outcome', 'scoring'] + clist + ['informative clusters won / tied / lost'], rows), '']
    # seed sensitivity
    rows, knife = [], []
    shares = {}
    for sd in range(SEED_SHARE_N):
        shares[sd] = verdicts_subset(raw, 'reg|', (), sd)
    for cid, name, a, b in PRIMARY[:2]:
        cells = []
        for c in clist:
            v = [shares[sd][(cid, c)][2] for sd in shares]
            cells.append('%d / %d / %d' % (v.count('win'), v.count('tie'), v.count('loss')))
            if cvp[(cid, c)][2] == 'tie' and (v.count('loss') or v.count('win')):
                knife.append((cid, c, v.count('win'), v.count('loss')))
        rows.append([cid, name] + cells)
    L += ['**Bootstrap seed sensitivity** (PRIMARY scoring; the cluster verdicts recomputed for the seeds 0 to %d of the stratified block bootstrap; cells: seeds giving win / tie / loss, of %d):' % (SEED_SHARE_N - 1, SEED_SHARE_N), '',
          X.table(['outcome', 'comparison'] + clist, rows), '']
    if knife:
        L += ['A tie at seed 0 that other seeds turn into a win or a loss is a knife-edge, not a tie: ' + '; '.join('%s at %s: %d of %d seeds give a win, %d a loss' % (cid, c, w, SEED_SHARE_N, l) for cid, c, w, l in knife) + '.', '']
    return L


def heldout_section(raw):
    hi = heldout_index(raw)
    if not hi:
        return ['## 3. Held-out test', '', 'No held-out result in the raw yet.', '']
    L = ['## 3. Held-out test (PRE-REGISTERED: `clusterdep_plan.md`, committed before the run)', '',
         'The clusters and patches were chosen by the rule of the plan (geometry only: `clusterdep_selection.txt`), the comparators are the frozen ones (nothing re-tuned), the repair is the merge at rho* = %g declared in the plan. '
         'PRIMARY scoring: the truth shifted by the rigid offset measured before any detection existed (the fresh protocol\'s secondary scoring); SECONDARY: the truth as it is. Intervals: 20 x 20 px block bootstrap, 2000 resamples, seed 0, '
         'stratified within patches, the same resamples for every row of a table.' % RHO_STAR, '']
    L += side_by_side_section(raw)
    for tname, tlabel in (('reg|', 'PRIMARY scoring (registration-corrected truth)'), ('', 'SECONDARY scoring (truth as it is)')):
        L += ['### 3.%d %s' % (1 if tname else 2, tlabel), '']
        cv = verdicts_heldout(raw, tname)
        clist = list(hi)
        rows = []
        for cid, name, a, b in PRIMARY:
            cells, ver = [], []
            for c in clist:
                dF, lohi, v = cv[(cid, c)]
                cells.append('%s %s' % (X.fmt_iv(dF, lohi), v))
                ver.append(v)
            rows.append([cid, name] + cells + ['%d / %d / %d' % (ver.count('win'), ver.count('tie'), ver.count('loss'))])
        L += [X.table(['outcome', 'comparison'] + ['%s (%d patch%s)' % (c, len(hi[c]), '' if len(hi[c]) == 1 else 'es') for c in clist] + ['clusters won / tied / lost'], rows), '']
        rowsp = []
        for cid, name, a, b in PRIMARY:
            df = pooled_heldout(raw, tname, a, b)
            rowsp.append([cid, name, X.fmt_iv(df[0], df[1]), X.verdict(*df)])
        L += ['Pooled over all %d held-out patches (descriptive; the patches of one cluster are not independent):' % sum(len(v) for v in hi.values()), '', X.table(['outcome', 'comparison', 'dF pooled', 'verdict'], rowsp), '']
        rowsf = []
        for c in clist:
            for label, key in HELD_ROWS:
                blk, sizes = hblocks(raw, c, hi[c], key, tname)
                idx = X.boot_idx_strat(sizes)
                (cc, cci), (pp, pci), (ff, fci) = X.boot_cpf(blk, idx)
                rowsf.append([c, label, X.fmt_iv(cc, cci), X.fmt_iv(pp, pci), X.fmt_iv(ff, fci)])
        L += [X.table(['cluster', 'catalogue', 'completeness', 'purity', 'F'], rowsf), '']
    rowsc = []
    for c, pl in hi.items():
        for p in pl:
            k = hkey(c, p)
            sh = raw[k + 'shift']
            cells = [fmt(tot_f(raw, k + t + r)[0]) for r in ('rule', 'crowdsource (front end)') for t in ('reg|', '')]
            rowsc.append([c, '(%d, %d)' % p, '(%+.2f, %+.2f)' % tuple(sh), '%.2f' % np.hypot(*sh)] + cells)
    L += ['### 3.3 The registration correction patch by patch', '',
          'The rigid offset was fitted to the pixels before any detection existed (a PSF-convolved copy of the truth list against the pixels). If the corrected scoring is much WORSE than the scoring as it is for every method, the fitted offset is '
          'wrong for that patch (a failed fit, not a property of the methods); the table lets that be read off.', '',
          X.table(['cluster', 'patch', 'offset (px)', '|offset|', 'F ISS frozen, corrected', 'F ISS frozen, as is', 'F crowdsource (front end), corrected', 'F crowdsource (front end), as is'], rowsc), '']
    return L


ROLE_OF = {'M2': 'tuning (the configuration was tuned on it)', 'M53': 'post hoc (the kernel rule was designed after its losses with the isotropic row)',
           'M92': 'post hoc (the kernel rule was designed after its losses with the isotropic row)', 'M3': 'spent (scored in an earlier round; the repaired column is post hoc)',
           'M15': 'spent (scored in an earlier round; the repaired column is post hoc)'}
HELD_ROLE = 'held out (pre-registered; geometry-selected, one amendment of the bright-star gate after a dry run)'


def iso_verdicts_heldout(raw, tname):
    """{cluster: verdict} of the originally frozen isotropic-kernel ISS row against the front-end crowdsource for every held-out cluster (patches pooled, stratified bootstrap, seed 0)."""
    out = {}
    for c, pl in heldout_index(raw).items():
        ba, sizes = hblocks(raw, c, pl, 'iso', tname)
        bb, _ = hblocks(raw, c, pl, 'crowdsource (front end)', tname)
        out[c] = X.verdict(*X.boot_diff(ba, bb, X.boot_idx_strat(sizes))[2])
    return out


def tally_section(raw):
    """The cluster tally of the ISS (kernel-rule row) against the front-end crowdsource over every cluster scored so far: the 5 spent clusters (POST HOC for the repaired row) and the held-out ones (primary scoring)."""
    spent, held = {}, {}
    if 'meta|labels' in raw:
        idx = spent_index(raw)
        pks = [pkey(l, p) for l, p in idx]
        clusters = np.array([CLUSTER_OF[l] for l, p in idx])
        for c in ['M2', 'M53', 'M92', 'M3', 'M15']:
            ii = [i for i in range(len(idx)) if clusters[i] == c]
            row = []
            for a in ('rule', 'rule_m%g' % RHO_STAR):
                ba, sizes = blocks_of(raw, [pks[i] + a for i in ii])
                bb, _ = blocks_of(raw, [pks[i] + 'cs' for i in ii])
                df = X.boot_diff(ba, bb, X.boot_idx_strat(sizes))[2]
                row.append((df[0], df[1], X.verdict(*df)))
            spent[c] = row
    cv = verdicts_heldout(raw, 'reg|') if heldout_index(raw) else {}
    for c in heldout_index(raw):
        held[c] = [cv[('O2f', c)], cv[('O2', c)]]
    if not spent and not held:
        return []
    L = ['## 4. Tally over every cluster scored so far (ISS kernel-rule row against the front-end crowdsource)', '',
         'Frozen = the ISS as frozen; repaired = after the merge at rho* = %g. For the five spent clusters the repaired column is POST HOC (rho* was chosen on their patches); for the held-out clusters both columns are the pre-registered outcomes O2f and O2 '
         '(registration-corrected truth)%s. A cluster counts as won / tied / lost by the 95 %% interval of its pooled paired dF.' % (
             RHO_STAR, '' if not held else '; %d of the %d held-out patches had a failed registration fit and carry no information under that scoring (section 3.0)' % (len(failed_registration(raw)), sum(len(v) for v in heldout_index(raw).values()))), '']
    rows = []
    for src, name in ((spent, 'spent'), (held, 'held out')):
        for c, (fz, rp) in src.items():
            rows.append([c, name, ROLE_OF[c] if name == 'spent' else HELD_ROLE, '%s %s' % (X.fmt_iv(fz[0], fz[1]), fz[2]), '%s %s' % (X.fmt_iv(rp[0], rp[1]), rp[2])])
    L += [X.table(['cluster', 'set', 'role', 'frozen ISS minus crowdsource (front end)', 'repaired ISS minus crowdsource (front end)'], rows), '']
    for label, src in (('the 5 spent clusters', spent), ('the held-out clusters', held), ('all clusters scored', {**spent, **held})):
        if not src:
            continue
        cells = []
        for j in (0, 1):
            v = [r[j][2] for r in src.values()]
            cells.append(count3(v))
        L.append('* %s (n = %d): frozen %s; repaired %s.' % (label, len(src), cells[0], cells[1]))
    groups = [('tuning cluster', [c for c in spent if ROLE_OF[c].startswith('tuning')]), ('post hoc design of the kernel rule', [c for c in spent if ROLE_OF[c].startswith('post hoc')]),
              ('other spent clusters', [c for c in spent if ROLE_OF[c].startswith('spent')]), ('held out', list(held))]
    allr = {**spent, **held}
    parts = ['%s (%s): frozen %s, repaired %s' % (lab, ', '.join(cs), count3([allr[c][0][2] for c in cs]), count3([allr[c][1][2] for c in cs])) for lab, cs in groups if cs]
    L.append('* By role: ' + '; '.join(parts) + '.  The 11 clusters are not a random sample; the repaired tally is partly post hoc (rho* was chosen on the 5 spent clusters), and for the kernel-rule row M53 and M92 are not independent of its design.')
    if spent and held:
        isp = []
        for c in spent:
            idxs = [i for i in range(len(spent_index(raw))) if CLUSTER_OF[spent_index(raw)[i][0]] == c]
            pks = [pkey(*spent_index(raw)[i]) for i in idxs]
            ba, sizes = blocks_of(raw, [pk + 'iso' for pk in pks])
            bb, _ = blocks_of(raw, [pk + 'cs' for pk in pks])
            isp.append(X.verdict(*X.boot_diff(ba, bb, X.boot_idx_strat(sizes))[2]))
        L.append('* The originally frozen isotropic-kernel row against the front-end crowdsource: spent clusters %s; held-out clusters %s (PRIMARY scoring) or %s (truth as it is).' % (
            count3(isp), count3(list(iso_verdicts_heldout(raw, 'reg|').values())), count3(list(iso_verdicts_heldout(raw, '').values()))))
    bad = failed_registration(raw) if held else set()
    if bad:
        left = [c for c in held if any((c, p) not in bad for p in heldout_index(raw)[c])]
        cvp = verdicts_subset(raw, 'reg|', bad)
        L.append('* Held-out clusters with a usable registration (%d of %d; the failed patches are excluded, post hoc): frozen %s; repaired %s (PRIMARY scoring).' % (
            len(left), len(held), count3([cvp[('O2f', c)][2] for c in left]), count3([cvp[('O2', c)][2] for c in left])))
    return L + ['']


def report(raw):
    L = ['# Why does the ISS lead over the front-end crowdsource depend on the cluster? Diagnosis, one repair, held-out test', '',
         'Generated by `experiments/exp_astro_clusterdep.py` from `clusterdep_raw.npz` alone; every number below is printed by the script. Sections 1 and 2 are POST HOC (the spent patches, read from the stored catalogues, nothing re-run); '
         'section 3 is the pre-registered held-out test of `clusterdep_plan.md`. Scoring as `exp_astro_m2`: match at 0.5 px and 0.5 mag, F-measure of completeness and purity at mag < 22, 95 % percentile intervals of a bootstrap over the '
         '20 x 20 px blocks (2000 resamples, seed 0).', '']
    if 'meta|labels' in raw:
        d, idx, pks, clusters = diag_section(raw)
        L += d + merge_section(raw, idx, pks, clusters) + r1_section(raw)
    L += heldout_section(raw) + tally_section(raw)
    return '\n'.join(L) + '\n'


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--diag', action='store_true')
    ap.add_argument('--r1', action='store_true')
    ap.add_argument('--freeze', action='store_true')
    ap.add_argument('--run', action='store_true')
    ap.add_argument('--from-raw', action='store_true')
    a = ap.parse_args()
    if a.diag:
        part_diag()
    if a.r1:
        part_r1()
    if a.freeze:
        write_freeze()
    if a.run:
        part_run()
    if a.from_raw:
        with np.load(RAW, allow_pickle=False) as z:
            raw = {k: z[k] for k in z.files}
        with open(MD, 'w') as f:
            f.write(report(raw))
        print('wrote', MD)


if __name__ == '__main__':
    main()
