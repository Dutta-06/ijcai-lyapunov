"""
The kernel of the sliding path against the PSF stamp: a POST HOC study after the frozen transfer test
====================================================================================================

What the frozen test showed (`exp_astro_m2.py`, part T2, results/astro/m2.md): the configuration chosen on the M2 development patches beat crowdsource and photutils on every M2 test
patch, but on the two clusters imaged with worse seeing (M53, M92; nothing re-tuned) it did not: M53 a small gain over crowdsource, M92 a loss.

What the images and the PSF stamps alone show (no catalogue, no truth): `imaging.fit_psf_parts` fits the frame's psField stamp by a mixture of three isotropic Gaussians. On
M2's stamps it does so to 0.1 % of the peak (rms); on the stamps of M53 and M92, which are elliptical and lumpy in the core (5 to 15 % of the peak in single pixels), only to 0.65 to 0.84 %,
whatever the start (a multi-start finds the same optimum). The path then explains the bright stars with extra units, and the standardized residual of its fit has a robust standard deviation of
1.2 to 2.3 instead of 1 (M2: 0.94 to 1.05), so that its significance levels are optimistic.

The kernel family is a cause that the algorithm can address (this script does not show it is the only one: the bright-star counts below use the truth and do not isolate it): the sliding path needs the kernel only as a sum of separable terms, and `psf_general.GeneralPixelPSF` gives each part its own
widths and offset. The rule `psf_general.fit_adequate` keeps the isotropic fit where its rms is at most 2e-3 of the peak (all six M2 stamps: nothing changes there, bit for bit) and takes the
general family fitted to the stamp otherwise (M53 and M92: 5 to 6 parts, rms 1.0e-3 to 1.7e-3).

This script runs the path with the kernel that rule chooses on the six second-cluster patches, at the SAME frozen configuration and a priori default, and reports it beside the frozen results
of the same patches. It is NOT a held-out test: the family, the rule and its goal were fixed after the transfer test had been seen (their only inputs are the stamps), on the same
patches. What it measures is how much of the loss the kernel change recovers on the same patches.

    python experiments/exp_astro_psf.py --stamps             the stamp fits of every patch (M2 development and test patches, M53, M92), label-free -> results/astro/psf_raw.npz
    python experiments/exp_astro_psf.py --run [--workers 4]  the path with the adequate kernel on the second-cluster patches
    python experiments/exp_astro_psf.py --resid              the residual scale of the frozen M2 catalogues (no path is run)
    python experiments/exp_astro_psf.py --from-raw           regenerates results/astro/psf.md (reads psf_raw.npz and the frozen results/astro/m2_raw.npz)
"""
import argparse
import os
import sys
import time

os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
from scipy.ndimage import gaussian_filter

HERE = os.path.abspath(os.path.dirname(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)
import astro_common as A                  # noqa: E402
import exp_astro_m2 as X                  # noqa: E402
from src.core import imaging as IM        # noqa: E402
from src.core import psf_general as G     # noqa: E402

OUT = os.environ.get('PSF_OUT', os.path.join(REPO, 'results', 'astro'))
RAW, MD = os.path.join(OUT, 'psf_raw.npz'), os.path.join(OUT, 'psf.md')
M2RAW = os.environ.get('PSF_M2_RAW', X.RAW)          # the frozen results, read only
CLUSTERS = ('m53', 'm92')
PATCHES = [('m2', p) for p in X.DEV + X.TESTS] + [(ds, p) for ds in CLUSTERS for p in X.SECOND[ds]['patches']]


def residual_sd(patch, P, det):
    """The robust standard deviation (1.4826 MAD about the median) of the standardized residual of the high-passed patch against a catalogue det (n, 3) = x, y, flux, over the
    interior (5 px from the border), for the kernel P (with its high-pass) and the weights of the path: 1 where the noise model and the kernel are right."""
    dn = patch['dn']
    expected = np.maximum(gaussian_filter(dn, 1.5, mode='reflect'), np.percentile(dn, 1))
    sw = 1.0 / np.sqrt(expected / patch['gain'] + patch['dark'])
    y = P.highpass(dn - np.median(dn))
    z = (y - P.synth(det[:, :2], det[:, 2])) * sw
    zi = z[5:-5, 5:-5]
    return float(1.4826 * np.median(np.abs(zi - np.median(zi))))


def run_path_psf(patch, hp, rex, fit, kappas=X.KAPPAS):
    """`exp_astro_m2.run_path` (time rule, dual order, clean-up by the units' own significance, plus the key 0.0: the BIC block) with the kernel of `fit` (`psf_general.fit_adequate`);
    the same code, the kernel built by `psf_general.make_psf`. Returns (kappa -> (x, y, flux), seconds, entrants, info), info = dict(bic_index, n_blocks, rsd: kappa -> residual_sd)."""
    ny, nx = patch['dn'].shape
    P = G.make_psf(fit, ny, nx, hp)
    L = IM.Lattice(P, h=X.LAT_H, margin=X.MARGIN)
    dn = patch['dn']
    expected = np.maximum(gaussian_filter(dn, 1.5, mode='reflect'), np.percentile(dn, 1))
    sw = 1.0 / np.sqrt(expected / patch['gain'] + patch['dark'])
    y = P.highpass(dn - np.median(dn))
    reach = P.reach()
    scale = IM.noise_scale(L, sw ** 2)
    t0 = time.time()
    R = IM.sliding_path(P, y, sw, L, max_entries=X.MAX_ENTRIES, tol=1e-6, passes=1, rtrunc=10.0 / reach, r_ex=rex, fit_radius=X.FIT_RADIUS, order='dual', tau_max=1.0 / (min(kappas) * scale))
    secs = time.time() - t0
    out = {}
    for k in kappas:
        i, K = IM.stop_time(R['hist'], R['taus'], 1.0 / (k * scale))
        h = R['hist'][i]
        keep = IM.unit_significance(P, sw, h['c'], h['w']) >= k
        out[k] = (h['c'][keep, 0].copy(), h['c'][keep, 1].copy(), h['w'][keep].copy())
    i, K = IM.stop_bic(R['hist'], ny * nx, n_par=3)
    h = R['hist'][i]
    kb = float(np.sqrt(3.0 * np.log(ny * nx)))
    keep = IM.unit_significance(P, sw, h['c'], h['w']) >= kb
    out[0.0] = (h['c'][keep, 0].copy(), h['c'][keep, 1].copy(), h['w'][keep].copy())
    rsd = {k: residual_sd(patch, P, np.c_[d[0], d[1], d[2]]) for k, d in out.items()}
    return out, secs, len(R['taus']), dict(bic_index=int(i), n_blocks=len(R['hist']), rsd=rsd)


def stamp_shape(stamp, radius=5.0):
    """The shape of a PSF stamp from its second moments inside `radius` pixels of the middle pixel: [ellipticity, orientation of the long axis in degrees from the x axis, sd along x,
    sd along y]; ellipticity = sqrt(e1^2 + e2^2) with e1 = (mxx - myy) / (mxx + myy), e2 = 2 mxy / (mxx + myy)."""
    st = np.asarray(stamp, float)
    n = st.shape[0]
    yy, xx = np.mgrid[0:n, 0:n]
    c = (n - 1) / 2.0
    w = st * (np.hypot(yy - c, xx - c) <= radius)
    w = w / w.sum()
    mxx, myy, mxy = (w * (xx - c) ** 2).sum(), (w * (yy - c) ** 2).sum(), (w * (xx - c) * (yy - c)).sum()
    e1, e2 = (mxx - myy) / (mxx + myy), 2 * mxy / (mxx + myy)
    return np.array([np.hypot(e1, e2), np.degrees(0.5 * np.arctan2(e2, e1)), np.sqrt(mxx), np.sqrt(myy)])


def stamp_arrays(patch):
    """The label-free summary of a patch's stamp fit: dict of arrays."""
    f = G.fit_adequate(patch['psf'])
    p5 = np.array(f['parts5'], float) if f['parts5'] is not None else np.zeros((0, 5))
    return dict(kind=np.array([f['kind']]), rms=np.array([f['rms_isotropic'], f['rms']]), sigma=np.array([f['sigma']]), parts5=p5, shape=stamp_shape(patch['psf'])), f


def part_stamps(out_path=RAW):
    for ds, pxy in PATCHES:
        arrays, f = stamp_arrays(X.load_patch(*pxy, ds=ds))
        X.save_part(out_path, 'S|%s|x%d_y%d' % (ds, pxy[0], pxy[1]), arrays)
        print('S %s %s: isotropic rms %.5f -> %s, rms %.5f, %d parts' % (ds, pxy, f['rms_isotropic'], f['kind'], f['rms'], len(f['parts5']) if f['parts5'] else len(f['parts'])), flush=True)


def part_resid(out_path=RAW, m2raw=M2RAW):
    """The residual scale of the FROZEN catalogues of the M2 test patches (the isotropic kernel, the stored detections at the reported significance): label-free, no path is run."""
    froz = X.load_raw(m2raw)
    for pxy in X.TESTS:
        patch = X.load_patch(*pxy, ds='m2')
        sig, parts, _ = IM.fit_psf_parts(patch['psf'], n_parts=3, sigma0=1.0)
        arrays = {}
        for label, h in (('ours', froz['T|x%d_y%d|cfg' % pxy][0]), ('default', X.DEFAULT[0])):
            P_iso = IM.PixelPSF(*patch['dn'].shape, sig, parts=parts, hp=h or None)
            arrays['%s|rsd' % label] = np.array([residual_sd(patch, P_iso, froz['T|x%d_y%d|%s|det' % (pxy[0], pxy[1], label)])])
        X.save_part(out_path, 'R|m2|x%d_y%d' % pxy, arrays)
        print('R m2 %s: residual sd %.3f (tuned configuration) %.3f (a priori default)' % (pxy, arrays['ours|rsd'][0], arrays['default|rsd'][0]), flush=True)


def _psf_patch_arrays(args):
    """The arrays of one second-cluster patch: the stamp fit, the path with the adequate kernel at the frozen configuration ('ours') and at the a priori default ('default'), scored against the
    truth as in `exp_astro_m2._test_patch_arrays`, with the residual scale of the frozen catalogue (the isotropic kernel) and of the new one at the reported significance."""
    ds, pxy, hp, rex, kap, m2raw = args
    patch = X.load_patch(*pxy, ds=ds)
    t = A.truth_for_patch(X.truth_arrays(ds), pxy[0], pxy[1])
    arrays, fit = stamp_arrays(patch)
    arrays['cfg'] = np.array([hp, rex, kap])
    froz = X.load_raw(m2raw)
    for label, (h, r), dk in (('ours', (hp, rex), kap), ('default', (X.DEFAULT[0], X.DEFAULT[1]), X.DEFAULT[2])):
        dets, secs, n, info = run_path_psf(patch, h, r, fit)
        print('P patch %s %s %s: sigma_b %g r_ex %g: kernel %s, %d entries in %.0f s, BIC block %d of %d' % (label, ds, pxy, h, r, fit['kind'], n, secs, info['bic_index'], info['n_blocks']), flush=True)
        arrays['%s|secs' % label] = np.array([secs, n])
        arrays['%s|bic' % label] = np.array([info['bic_index'], info['n_blocks'], int(n >= X.MAX_ENTRIES)])
        for k, d in dets.items():
            arrays.update(X._pack('%s|k%g' % (label, k), X.score_dets(d, t)))
        arrays['%s|det' % label] = np.c_[dets[dk][0], dets[dk][1], dets[dk][2]]
        sig, parts, _ = IM.fit_psf_parts(patch['psf'], n_parts=3, sigma0=1.0)
        P_iso = IM.PixelPSF(*patch['dn'].shape, sig, parts=parts, hp=h or None)
        arrays['%s|rsd' % label] = np.array([residual_sd(patch, P_iso, froz['T2|%s|x%d_y%d|%s|det' % (ds, pxy[0], pxy[1], label)]), info['rsd'][dk]])
    return arrays


def part_run(out_path=RAW, workers=1, m2raw=M2RAW):
    """The path with the adequate kernel on each second-cluster patch (a patch whose result exists is not run again), at the configuration frozen on the M2 development patches."""
    froz = X.load_raw(m2raw)
    X.check_freeze(froz)
    hp, rex, kap, f = X.choose(froz)
    cur = X.load_raw(out_path) if os.path.exists(out_path) else {}
    todo = [(ds, pxy) for ds in CLUSTERS for pxy in X.SECOND[ds]['patches'] if 'P|%s|x%d_y%d|cfg' % (ds, pxy[0], pxy[1]) not in cur]
    jobs = [(ds, pxy, hp, rex, kap, m2raw) for ds, pxy in todo]
    if workers > 1 and len(jobs) > 1:
        from multiprocessing import Pool
        with Pool(min(workers, len(jobs))) as pool:
            results = pool.map(_psf_patch_arrays, jobs, chunksize=1)
    else:
        results = [_psf_patch_arrays(j) for j in jobs]
    for (ds, pxy), arrays in zip(todo, results):
        X.save_part(out_path, 'P|%s|x%d_y%d' % (ds, pxy[0], pxy[1]), arrays)


# ----------------------------------------------------------------------------------------------------------------------------------
# report
# ----------------------------------------------------------------------------------------------------------------------------------
FROZEN_ROWS = [('ISS as frozen (isotropic kernel), tuned on the M2 development patches (kappa 2)', 'ours|k2'), ('ISS as frozen, a priori default (sigma_b 5, r_ex 1.2, kappa 5)', 'default|k5'),
               ('ISS as frozen, a priori default with the BIC stop', 'default|k0')]
NEW_ROWS = [('ISS with the adequate kernel, the same configuration (kappa 2)', 'ours|k2'), ('ISS with the adequate kernel, a priori default (kappa 5)', 'default|k5'),
            ('ISS with the adequate kernel, a priori default with the BIC stop', 'default|k0')]
BASE_ROWS = [('crowdsource (default)', 'crowdsource (default)'), ('crowdsource (tuned)', 'crowdsource (tuned)'), ('photutils (default)', 'photutils (default)'), ('photutils (tuned)', 'photutils (tuned)')]
PAIRS = [('adequate kernel vs frozen, tuned configuration', 'new', 'ours|k2', 'froz', 'ours|k2'), ('adequate kernel vs frozen, a priori default', 'new', 'default|k5', 'froz', 'default|k5'),
         ('adequate kernel vs frozen, BIC stop', 'new', 'default|k0', 'froz', 'default|k0'), ('adequate kernel (tuned configuration) vs crowdsource (tuned)', 'new', 'ours|k2', 'froz', 'crowdsource (tuned)'),
         ('adequate kernel (tuned configuration) vs photutils (tuned)', 'new', 'ours|k2', 'froz', 'photutils (tuned)'),
         ('adequate kernel (a priori default) vs crowdsource (default)', 'new', 'default|k5', 'froz', 'crowdsource (default)'),
         ('adequate kernel (BIC stop) vs crowdsource (default)', 'new', 'default|k0', 'froz', 'crowdsource (default)')]


NB = 9                                                  # the magnitude bins below 19.5 (`astro_common.BIN_EDGES[9]`)


def bright(raw, pre, key):
    """(matched, in the truth) of the stars brighter than 19.5 mag for the catalogue `key` of the patch with key prefix `pre`, or None if the raw has none."""
    if pre + key + '|n_true' not in raw:
        return None
    return int(raw[pre + key + '|n_matched'][:NB].sum()), int(raw[pre + key + '|n_true'][:NB].sum())


def _blk(raw, ds, pxy, which, key):
    pre = ('P|%s|x%d_y%d|' if which == 'new' else 'T2|%s|x%d_y%d|') % (ds, pxy[0], pxy[1])
    return raw[pre + key + '|blk22'] if pre + key + '|blk22' in raw else None


def report(raw, m2raw=None):
    """The markdown report from the raw of this script and the frozen raw (read only)."""
    froz = m2raw if m2raw is not None else X.load_raw(M2RAW)
    L = ['# The kernel of the sliding path against the PSF stamp: a post hoc study after the frozen transfer test', '',
         'Generated by `experiments/exp_astro_psf.py` from `psf_raw.npz` and the frozen `m2_raw.npz`; every number below is printed by the script. **POST HOC, NOT A HELD-OUT TEST.** The frozen '
         'transfer test (`m2.md`, M53 and M92 with the configuration chosen on M2) was seen first; the kernel family, the adequacy rule (the isotropic fit is kept where its rms residual against the '
         'stamp\'s peak is at most %g, which holds for all six M2 stamps) and its goal were fixed afterwards from the stamps alone (no catalogue, no truth), and the same patches are run again. '
         'The comparison with the frozen results says how much of the loss a change of kernel recovers on these patches, not how the method would do on new data, and the verdict columns below are descriptive for the same reason.' % G.RMS_GOAL, '']
    rows = []
    for ds, pxy in PATCHES:
        k = 'S|%s|x%d_y%d|' % (ds, pxy[0], pxy[1])
        if k + 'kind' not in raw:
            continue
        rms, kind, p5 = raw[k + 'rms'], str(raw[k + 'kind'][0]), raw[k + 'parts5']
        sh = raw.get(k + 'shape')
        rows.append([ds, '(%d, %d)' % pxy, '%.5f' % rms[0], kind, '%d' % (len(p5) if kind == 'general' else 3), '%.5f' % rms[1], '-' if sh is None else '%.3f' % sh[0], '-' if sh is None else '%.0f' % sh[1]])
    L += ['## The stamp fits (label-free)\n', 'The rms residual of the fit against the stamp\'s peak: of the isotropic three-part mixture (all the frozen results), the kind chosen by the rule, the number of parts and the rms of '
          'the kernel chosen; the ellipticity of the stamp from its second moments inside 5 px and the orientation of its long axis from the x axis (degrees).\n',
          X.table(['set', 'patch', 'isotropic rms', 'kernel chosen', 'parts', 'rms of the kernel chosen', 'ellipticity', 'orientation'], rows), '']
    for ds in CLUSTERS:
        plist = [p for p in X.SECOND[ds]['patches'] if 'P|%s|x%d_y%d|cfg' % (ds, p[0], p[1]) in raw]
        if not plist:
            continue
        title = X.SECOND[ds]['title']
        kap = float(raw['P|%s|x%d_y%d|cfg' % (ds, plist[0][0], plist[0][1])][2])
        sources = {'new': raw, 'froz': froz}
        L += ['## %s: the %d test patches\n' % (title, len(plist))]
        for pxy in plist:
            L += ['### Patch (%d, %d)\n' % pxy]
            idx = X.boot_idx(25)
            rows = []
            for label, key in FROZEN_ROWS + NEW_ROWS + BASE_ROWS:
                which = 'new' if (label, key) in NEW_ROWS else 'froz'
                b = _blk(sources[which], ds, pxy, which, key)
                if b is None:
                    continue
                (c, cci), (pu, pci), (ff, fci) = X.boot_cpf(b, idx)
                rows.append([label, X.fmt_iv(c, cci), X.fmt_iv(pu, pci), X.fmt_iv(ff, fci)])
            L += [X.table(['catalogue', 'completeness < 22', 'purity < 22', 'F'], rows), '']
        blk = {}
        for which, rowsdef in (('froz', FROZEN_ROWS + BASE_ROWS), ('new', NEW_ROWS)):
            for label, key in rowsdef:
                parts = [_blk(sources[which], ds, p, which, key) for p in plist]
                if all(q is not None for q in parts):
                    blk[(which, key)] = np.concatenate(parts)
        idx = X.boot_idx_strat([25] * len(plist))
        L += ['### All %d patches together\n' % len(plist), 'The blocks of the patches pooled (%d blocks); the bootstrap resamples the blocks within each patch; the same resamples for every row.\n' % (25 * len(plist))]
        rows = []
        for label, key in FROZEN_ROWS + NEW_ROWS + BASE_ROWS:
            which = 'new' if (label, key) in NEW_ROWS else 'froz'
            if (which, key) not in blk:
                continue
            tot = blk[(which, key)].sum(0)
            (c, cci), (pu, pci), (ff, fci) = X.boot_cpf(blk[(which, key)], idx)
            rows.append([label, '%d / %d' % (tot[1], tot[0]), X.fmt_iv(c, cci), '%d / %d' % (tot[3], tot[2]), X.fmt_iv(pu, pci), X.fmt_iv(ff, fci)])
        L += [X.table(['catalogue', 'matched / truth < 22', 'completeness', 'matched / detected < 22', 'purity', 'F'], rows), '', 'Paired differences on the pooled blocks:\n']
        rows = []
        for name, wa, ka, wb, kb in PAIRS:
            if (wa, ka) in blk and (wb, kb) in blk:
                dc, dp, df = X.boot_diff(blk[(wa, ka)], blk[(wb, kb)], idx)
                rows.append([name, X.fmt_iv(dc[0], dc[1]), X.fmt_iv(dp[0], dp[1]), X.fmt_iv(df[0], df[1]), X.verdict(*df)])
        L += [X.table(['comparison', 'd completeness', 'd purity', 'd F', 'verdict on F'], rows), '']
    rows = []
    for ds, plist in [('m2', X.TESTS)] + [(c, X.SECOND[c]['patches']) for c in CLUSTERS]:
        for pxy in plist:
            if ds == 'm2':
                pre = 'R|m2|x%d_y%d|' % pxy
                if pre + 'ours|rsd' in raw:
                    rows.append([ds, '(%d, %d)' % pxy, '%.2f' % raw[pre + 'ours|rsd'][0], '-', '%.2f' % raw[pre + 'default|rsd'][0], '-'])
            else:
                pre = 'P|%s|x%d_y%d|' % (ds, pxy[0], pxy[1])
                if pre + 'ours|rsd' in raw:
                    rows.append([ds, '(%d, %d)' % pxy, '%.2f' % raw[pre + 'ours|rsd'][0], '%.2f' % raw[pre + 'ours|rsd'][1], '%.2f' % raw[pre + 'default|rsd'][0], '%.2f' % raw[pre + 'default|rsd'][1]])
    if rows:
        L += ['## The residual scale of the fitted catalogues (label-free)\n', 'The robust standard deviation of the standardized residual of the fitted catalogue (1: the noise model and the kernel are right), at the reported significance, for '
              'the frozen isotropic kernel and for the adequate kernel (the M2 kernels are unchanged by the rule: no second value).\n',
              X.table(['set', 'patch', 'frozen, tuned configuration', 'adequate kernel, same', 'frozen, a priori default', 'adequate kernel, same'], rows), '']
    rows = []
    for ds, plist in [('m2', X.TESTS)] + [(c, X.SECOND[c]['patches']) for c in CLUSTERS]:
        for pxy in plist:
            pf = ('T|x%d_y%d|' if ds == 'm2' else 'T2|' + ds + '|x%d_y%d|') % pxy
            pn = 'P|%s|x%d_y%d|' % (ds, pxy[0], pxy[1])
            cells = []
            for src, pre, key in ((froz, pf, 'ours|k2'), (raw, pn, 'ours|k2'), (froz, pf, 'crowdsource (tuned)'), (froz, pf, 'photutils (tuned)')):
                b = bright(src, pre, key)
                cells.append('-' if b is None else '%d / %d' % b)
            if cells[0] != '-':
                rows.append([ds, '(%d, %d)' % pxy] + cells)
    if rows:
        L += ['## Stars brighter than 19.5 mag: matched / in the truth\n', 'Descriptive (the stars with the best signal to noise; a star split in several units is not matched; these counts use the truth and do not show that the kernel is the cause of the splitting).\n',
              X.table(['set', 'patch', 'ISS as frozen (tuned, kappa 2)', 'ISS with the adequate kernel (same)', 'crowdsource (tuned)', 'photutils (tuned)'], rows), '']
    return '\n'.join(L) + '\n'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--stamps', action='store_true')
    ap.add_argument('--run', action='store_true')
    ap.add_argument('--resid', action='store_true')
    ap.add_argument('--from-raw', action='store_true')
    ap.add_argument('--workers', type=int, default=int(os.environ.get('ASTRO_WORKERS', '3')))
    A_ = ap.parse_args()
    if A_.stamps:
        part_stamps()
    elif A_.run:
        part_run(workers=A_.workers)
    elif A_.resid:
        part_resid()
    elif A_.from_raw:
        with open(MD, 'w') as f:
            f.write(report(X.load_raw(RAW)))
        print('wrote', MD)
    else:
        ap.error('give --stamps, --run or --from-raw')


if __name__ == '__main__':
    main()
