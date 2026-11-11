"""
Two fresh globular clusters for the cluster-level test: the selection rule, the patches, the registry (SDSS r band against the HST HUGS catalogue)
=================================================================================================================================================

The astronomy evidence of the project is three clusters (M2 tuned; M53 and M92 transfer), one SDSS frame each. `results/astro/fresh_rule.md` fixes, before anything was executed, the rule that
picks two more clusters; this module is the rule (`select`), the patch builder (`build`, the patch rules of `astro_second.py` imported and applied unchanged) and the registry that the experiment
scripts read. Nothing here runs a cataloguing method and no truth is compared with a detection. It imports `astro_second.py` (hash-frozen, read only).

The rule (version 1; see `fresh_rule.md` for the reasons and for the amendments A1 to A3 that were needed), on the output of `astro_second.py survey` (stored in `results/astro/fresh_survey.txt`) and the psField of each cluster's best frame:
  1. frame = the best field of the survey for the cluster (the cluster is dropped when the survey finds none);
  2. PSF_STATUS <= 2 and PSF_NSTAR(r) >= 1;
  3. FWHM of the psField stamp at the frame centre (column 1023.5, row 744.0) in [2.3, 3.6] px;
  4. number of survey windows >= 100;
  5. order by the number of windows (descending, ties: the smaller NGC number); the first two that also yield >= 2 patches (6) are the fresh clusters;
  6. the patch rules of M53 and M92 plus the patch-centre stamp FWHM in [2.3, 3.6] px.
Universe: the clusters of the coverage survey with a downloaded HUGS catalogue, without NGC 7089, 5024 and 6341.

Amendments A1 to A3 (`fresh_rule.md`, declared when version 1 of the rule left no usable cluster): the frames of a cluster are ALL the field boxes of the survey loop (`fields`, `select-frames`; `survey_fields`,
`apply_rule_frames`), the frame is chosen among those passing items 2 to 4 (most windows first), a cluster with no patch in its frame is tried in its next passing frame, and a cluster with a single patch stays.

    python -I experiments/astro_fresh.py select --survey results/astro/fresh_survey.txt [--out FILE]     version 1 of the rule on the survey's best frames (needs the psFields; domain environment: astropy)
    python -I experiments/astro_fresh.py fields --out results/astro/fresh_fields.json                     A1: every field box of the universe's clusters
    python -I experiments/astro_fresh.py select-frames --fields results/astro/fresh_fields.json [--out FILE] [--list-psfield]    A1: the rule per frame
    python -I experiments/astro_fresh.py build --cluster ngc7078 --frame 2507:6:39 [--dry-run]            truth, patches, table (needs the frame; $ASTRO_DATA_4 for the output)
    python -I experiments/astro_fresh.py registry --clusters ngc5272 ngc7078 --out results/astro/fresh_patches.json
"""
import argparse
import json
import os
import re
import sys

import numpy as np

HERE = os.path.abspath(os.path.dirname(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
import astro_second as S2                 # noqa: E402

ROOT4 = os.environ.get('ASTRO_DATA_4', '/home/user/astro_data_4')
ROOT2 = os.environ.get('ASTRO_DATA_2', '/home/user/astro_data_2')
REGISTRY = os.environ.get('FRESH_REGISTRY', os.path.join(REPO, 'results', 'astro', 'fresh_patches.json'))
USED = ('ngc7089', 'ngc5024', 'ngc6341')                                      # M2, M53, M92
UNIVERSE = ('ngc5053', 'ngc5272', 'ngc5466', 'ngc5904', 'ngc6205', 'ngc6254', 'ngc6535', 'ngc6838', 'ngc6934', 'ngc6981', 'ngc7078')
RULE = dict(status_max=2, nstar_min=1, fwhm_lo=2.3, fwhm_hi=3.6, windows_min=100, centre=(1023.5, 744.0), min_patches=2, n_fresh=2)
SURVEY_LINE = re.compile(r'^(ngc\d+)\s+.*?\|\s+(\d+)\s+(\d)\s+(\d+)\s+\|\s+([\d.]+)\s+([\d.]+)\s+(\d+)\s+(\S+)\s+\|\s+psp_status\s+(\d+)\s+psf_nstar\(r\)\s+(\d+)\s+sat\s+(\d+)')


# ----------------------------------------------------------------------------- the rule
def parse_survey(text):
    """The rows of `astro_second.py survey` output that name a best field: list of dict(cluster, run, camcol, field, seeing, frac, windows, density, status, nstar, sat)."""
    out = []
    for line in text.splitlines():
        m = SURVEY_LINE.match(line)
        if m:
            c, run, cam, fld, see, frac, nw, dens, ps, nst, sat = m.groups()
            out.append(dict(cluster=c, run=int(run), camcol=int(cam), field=int(fld), seeing=float(see), frac=float(frac), windows=int(nw), density=dens, status=int(ps), nstar=int(nst), sat=float(sat)))
    return out


def psfield_path(cluster, run, camcol, field, root=None):
    return os.path.join(root or ROOT4, cluster, 'sdss_psfield_%06d_%d_%04d' % (run, camcol, field), 'psField-%06d-%d-%04d.fit' % (run, camcol, field))


def stamp_fwhm(kl, xpos, ypos):
    """FWHM (px) of the psField reconstruction at (column xpos, row ypos): `astro_second.psf_fwhm` of the unit-sum 51 x 51 stamp."""
    return S2.psf_fwhm(S2.unit_psf(S2.sdss_psf_recon(kl, xpos, ypos)))


def apply_rule(rows, fwhm, rule=None, universe=UNIVERSE):
    """The rule of `fresh_rule.md` items 1 to 5 on survey rows and the stamp FWHM per cluster (dict cluster -> px, absent = unknown): list of dict(row + fwhm, ok_psf, ok_seeing, ok_windows, passes) for
    every cluster of the universe that has a survey row, in rank order (passing first, by windows descending then NGC number; then the rest by NGC number). A cluster of the universe without a row is dropped
    by item 1 and appears in `missing` of the result's second value."""
    rule = rule or RULE
    byc = {r['cluster']: r for r in rows}
    out, missing = [], []
    for c in universe:
        if c not in byc:
            missing.append(c)
            continue
        r = dict(byc[c])
        f = fwhm.get(c)
        r['fwhm'] = f
        r['ok_psf'] = r['status'] <= rule['status_max'] and r['nstar'] >= rule['nstar_min']
        r['ok_seeing'] = f is not None and rule['fwhm_lo'] <= f <= rule['fwhm_hi']
        r['ok_windows'] = r['windows'] >= rule['windows_min']
        r['passes'] = bool(r['ok_psf'] and r['ok_seeing'] and r['ok_windows'])
        out.append(r)
    ngc = lambda r: int(r['cluster'][3:])
    out.sort(key=lambda r: (not r['passes'], -r['windows'] if r['passes'] else 0, ngc(r)))
    return out, missing


def choose_fresh(ranked, patch_counts, rule=None):
    """Item 5: walking the passing clusters in rank order, the first `n_fresh` with at least `min_patches` patches (patch_counts: cluster -> number of patches, None = not tried).
    Returns (chosen clusters, the clusters tried and dropped for too few patches)."""
    rule = rule or RULE
    chosen, dropped = [], []
    for r in ranked:
        if not r['passes'] or len(chosen) >= rule['n_fresh']:
            continue
        n = patch_counts.get(r['cluster'])
        if n is None:
            break                                                           # the next cluster must be built before the order can go on
        (chosen if n >= rule['min_patches'] else dropped).append(r['cluster'])
    return chosen, dropped


def cmd_select(args):
    with open(args.survey) as f:
        rows = parse_survey(f.read())
    fwhm = {}
    for r in rows:
        p = psfield_path(r['cluster'], r['run'], r['camcol'], r['field'])
        if os.path.exists(p):
            fwhm[r['cluster']] = stamp_fwhm(S2.read_psfield(p)['kl'], *RULE['centre'])
    ranked, missing = apply_rule(rows, fwhm)
    L = ['The rule of fresh_rule.md, items 1 to 5, on %s' % args.survey, 'rule constants: %s' % json.dumps(RULE),
         'universe: %s ; used already (excluded): %s' % (', '.join(UNIVERSE), ', '.join(USED)), '',
         '%-8s %-18s %7s %6s %7s %-10s %6s %6s %9s | %-6s %-8s %-9s | %s' % ('cluster', 'frame (run cam fld)', 'seeing"', 'foot', 'windows', 'density', 'status', 'nstar', 'FWHM px', 'PSF', 'seeing', 'windows', 'passes')]
    for r in ranked:
        L.append('%-8s %-18s %7.3f %6.2f %7d %-10s %6d %6d %9s | %-6s %-8s %-9s | %s' % (
            r['cluster'], '%d %d %d' % (r['run'], r['camcol'], r['field']), r['seeing'], r['frac'], r['windows'], r['density'], r['status'], r['nstar'],
            '%.3f' % r['fwhm'] if r['fwhm'] is not None else 'n/a', 'ok' if r['ok_psf'] else 'FAIL', 'ok' if r['ok_seeing'] else 'FAIL', 'ok' if r['ok_windows'] else 'FAIL', 'PASS' if r['passes'] else 'no'))
    if missing:
        L.append('dropped by item 1 (the survey names no field box for them): %s' % ', '.join(missing))
    L.append('')
    L.append('order of the passing clusters (item 5): %s' % ', '.join(r['cluster'] for r in ranked if r['passes']))
    text = '\n'.join(L) + '\n'
    print(text)
    if args.out:
        with open(args.out, 'w') as f:
            f.write(text)


# ----------------------------------------------------------------------------- amendment A1: every frame of a cluster (fresh_rule.md)
def survey_fields(cluster, centre, runs, hugs, root2=None):
    """Every field box that the loop of `astro_second.cmd_survey` finds for a cluster (at least 2 % of the F606W < 22 footprint inside the box), with its survey columns: list of dict(cluster, run, camcol, field,
    seeing (r PSF_WIDTH, arcsec), frac, windows, dmin, dmax, status, nstar, sat). `centre`: dict(ra, dec); `runs`: the photoRunAll arrays; `hugs`: `astro_second.read_hugs` of the cluster. The loop is that of
    `cmd_survey` (same functions, same cuts); it keeps all the records instead of the best."""
    from astropy.io import fits
    root2 = root2 or ROOT2
    ok = (hugs['m606'] > -90) & (hugs['m606'] < 90)
    ra, dec, m = hugs['ra'][ok], hugs['dec'][ok], hugs['m606'][ok]
    out = []
    for rc in sorted({r['run'] for r in S2.run_candidates(runs, centre['ra'], centre['dec'])}):
        i = int(np.flatnonzero(runs['RUN'] == rc)[0])
        node, incl = float(runs['NODE'][i]), float(runs['INCL'][i])
        for cam in range(1, 7):
            pfp = os.path.join(root2, 'survey', 'sdss_photofield_%06d_%d' % (rc, cam), 'photoField-%06d-%d.fits' % (rc, cam))
            if not os.path.exists(pfp):
                continue
            pf = fits.getdata(pfp, 1)
            mu0, nu0 = S2.eq_to_gc(centre['ra'], centre['dec'], node, incl)
            dm = (mu0 - pf['MU_START'] + 180) % 360 - 180
            near = np.flatnonzero((dm > -0.4) & (dm < 0.55) & (nu0 > np.minimum(pf['NU_START'], pf['NU_END']) - 0.5) & (nu0 < np.maximum(pf['NU_START'], pf['NU_END']) + 0.5))
            for k in near:
                box = {q: float(pf[q][k]) for q in ('MU_START', 'MU_END', 'NU_START', 'NU_END')}
                x, y, inside, w, h = S2.sources_in_field(ra, dec, node, incl, box)
                b = m < S2.MAG_DEPTH
                frac = float((inside & b).sum() / max(b.sum(), 1))
                if frac < 0.02:
                    continue
                tr = dict(px=x[inside], py=y[inside], m606=m[inside])
                x0s, y0s, cnt = S2.window_counts(tr['px'], tr['py'], tr['m606'], (int(h), int(w)))
                occ = S2.footprint_cells(tr['px'], tr['py'], (int(h), int(w)))
                dens = []
                for jj, yy in enumerate(y0s):
                    for ii, xx in enumerate(x0s):
                        d = cnt[jj, ii] / 1e4
                        if S2.CLASSES[0][0] <= d <= S2.CLASSES[-1][1] and xx >= 30 and yy >= 30 and xx + 100 <= w - 30 and yy + 100 <= h - 30 and S2.inside_footprint(occ, int(xx), int(yy)):
                            dens.append(d)
                out.append(dict(cluster=cluster, run=int(rc), camcol=int(cam), field=int(pf['FIELD'][k]), seeing=float(pf['PSF_WIDTH'][k][2]), frac=frac, windows=len(dens),
                                dmin=float(min(dens)) if dens else None, dmax=float(max(dens)) if dens else None, status=int(pf['PSP_STATUS'][k]), nstar=int(pf['PSF_NSTAR'][k][2]),
                                sat=float(pf['SATURATION_LEVEL'][k][2])))
    return out


def cmd_fields(args):
    """Writes the all-frames table (A1) for the clusters of the universe: JSON with every field box found, and a text table."""
    from astropy.io import fits
    sums = json.load(open(os.path.join(REPO, 'data', 'astro', 'checksums_second.json')))
    centres = sums['cluster_centres']
    t = fits.getdata(os.path.join(ROOT2, 'survey', 'photorunall', 'photoRunAll-dr17.fits'), 1)
    runs = {k: np.asarray(t[k]) for k in ('RUN', 'NODE', 'INCL', 'MU_REF', 'FIELD_REF', 'STARTFIELD', 'ENDFIELD', 'STRIPE', 'STRIP')}
    rows = []
    for c in (args.clusters or UNIVERSE):
        cat = S2.read_hugs(os.path.join(ROOT2, c, 'hugs_' + c, 'hlsp_hugs_hst_wfc3-uvis-acs-wfc_%s_multi_v1_catalog-meth1.txt' % c))
        rows += survey_fields(c, centres[c], runs, cat)
        print(c, sum(1 for r in rows if r['cluster'] == c), 'field boxes', flush=True)
    with open(args.out, 'w') as f:
        json.dump(dict(note='every field box of the survey loop with at least 2 % of the F606W < 22 footprint (astro_fresh.py fields)', fields=rows), f, indent=1)
    print('wrote', args.out, len(rows), 'fields')


def frame_id(r):
    return '%s:%d:%d:%d' % (r['cluster'], r['run'], r['camcol'], r['field'])


def apply_rule_frames(fields, fwhm, rule=None, universe=UNIVERSE):
    """Amendment A1: items 2, 3, 4 per frame (`fields`: list of `survey_fields` records; fwhm: frame id -> px, absent = not measured because items 2 and 4 failed); the frame of a cluster is the passing frame with
    the most windows (ties: smaller PSF_WIDTH, then run, camcol, field); clusters ordered by their frame's windows (descending, ties: smaller NGC number). Returns (ranked records of the chosen frames, every
    frame with its flags)."""
    rule = rule or RULE
    allf = []
    for r in fields:
        if r['cluster'] not in universe:
            continue
        r = dict(r)
        f = fwhm.get(frame_id(r))
        r['fwhm'] = f
        r['ok_psf'] = r['status'] <= rule['status_max'] and r['nstar'] >= rule['nstar_min']
        r['ok_seeing'] = f is not None and rule['fwhm_lo'] <= f <= rule['fwhm_hi']
        r['ok_windows'] = r['windows'] >= rule['windows_min']
        r['passes'] = bool(r['ok_psf'] and r['ok_seeing'] and r['ok_windows'])
        allf.append(r)
    best = {}
    for r in allf:
        if r['passes']:
            key = (-r['windows'], r['seeing'], r['run'], r['camcol'], r['field'])
            if r['cluster'] not in best or key < best[r['cluster']][0]:
                best[r['cluster']] = (key, r)
    ranked = sorted((v[1] for v in best.values()), key=lambda r: (-r['windows'], int(r['cluster'][3:])))
    return ranked, allf


def need_psfield(fields, rule=None, universe=UNIVERSE):
    """The frames of `fields` that pass items 2 and 4 (the psField is fetched only for these)."""
    rule = rule or RULE
    return [r for r in fields if r['cluster'] in universe and r['status'] <= rule['status_max'] and r['nstar'] >= rule['nstar_min'] and r['windows'] >= rule['windows_min']]


def cmd_select_frames(args):
    with open(args.fields) as f:
        fields = json.load(f)['fields']
    if args.list_psfield:
        for r in need_psfield(fields):
            print(frame_id(r))
        return
    fwhm = {}
    for r in need_psfield(fields):
        p = psfield_path(r['cluster'], r['run'], r['camcol'], r['field'])
        if os.path.exists(p):
            fwhm[frame_id(r)] = stamp_fwhm(S2.read_psfield(p)['kl'], *RULE['centre'])
    ranked, allf = apply_rule_frames(fields, fwhm)
    L = ['Amendment A1 of fresh_rule.md: items 2, 3, 4 on every frame, the frame of a cluster = the passing frame with the most windows', 'rule constants: %s' % json.dumps(RULE), '',
         '%-8s %-18s %7s %6s %7s %-10s %6s %6s %9s | %-6s %-8s %-9s | %s' % ('cluster', 'frame (run cam fld)', 'seeing"', 'foot', 'windows', 'density', 'status', 'nstar', 'FWHM px', 'PSF', 'seeing', 'windows', 'passes')]
    for r in sorted(allf, key=lambda r: (int(r['cluster'][3:]), -r['windows'], r['run'], r['camcol'], r['field'])):
        L.append('%-8s %-18s %7.3f %6.2f %7d %-10s %6d %6d %9s | %-6s %-8s %-9s | %s' % (
            r['cluster'], '%d %d %d' % (r['run'], r['camcol'], r['field']), r['seeing'], r['frac'], r['windows'], ('%.2f-%.2f' % (r['dmin'], r['dmax'])) if r['dmin'] is not None else '-', r['status'], r['nstar'],
            '%.3f' % r['fwhm'] if r['fwhm'] is not None else 'n/a (2 or 4 failed)', 'ok' if r['ok_psf'] else 'FAIL', 'ok' if r['ok_seeing'] else ('FAIL' if r['fwhm'] is not None else 'n/a'),
            'ok' if r['ok_windows'] else 'FAIL', 'PASS' if r['passes'] else 'no'))
    L += ['', 'frame of each cluster that has a passing frame, in the order of item 5: ' + ', '.join('%s (%d windows)' % (frame_id(r), r['windows']) for r in ranked)]
    text = '\n'.join(L) + '\n'
    print(text)
    if args.out:
        with open(args.out, 'w') as f:
            f.write(text)


# ----------------------------------------------------------------------------- paths and loaders (the registry)
def frame_prefix(run, camcol, field):
    return 'patch_r_%06d_%d_%04d' % (run, camcol, field)


def load_registry(path=None):
    """The registry written by `registry`: {cluster: dict(title, run, camcol, field, patches [[x0, y0], ...], root, prefix, truth, ...)}."""
    with open(path or REGISTRY) as f:
        return json.load(f)['clusters']


def patch_path(x0, y0, entry, root=None):
    return os.path.join(root or entry['root'], '%s_x%d_y%d_100x100' % (entry['prefix'], x0, y0), 'patch_arrays.npz')


def load_patch(x0, y0, entry, root=None):
    """dict(dn, psf, gain, dark) exactly like `exp_astro_m2.load_patch`."""
    z = np.load(patch_path(x0, y0, entry, root))
    return dict(dn=z['dn'].astype(float), psf=z['psf51'].astype(float), gain=float(z['gain']), dark=float(z['dark_variance']))


def truth_arrays(entry, root=None):
    z = np.load(os.path.join(root or entry['root'], entry['truth']))
    return dict(px=z['px'], py=z['py'], m606=z['m606'])


# ----------------------------------------------------------------------------- the patches
def fwhm_filter(cands, kl, rule=None):
    """The candidates (dicts with x0, y0) whose patch-centre stamp FWHM is in [fwhm_lo, fwhm_hi]; each kept candidate gets key 'fwhm'. Returns (kept, number removed)."""
    rule = rule or RULE
    kept = []
    for q in cands:
        f = stamp_fwhm(kl, q['x0'] + (S2.N_PATCH - 1) / 2.0, q['y0'] + (S2.N_PATCH - 1) / 2.0)
        if rule['fwhm_lo'] <= f <= rule['fwhm_hi']:
            kept.append(dict(q, fwhm=f))
    return kept, len(cands) - len(kept)


def cluster_paths(cluster, run, camcol, field):
    return dict(frame=os.path.join(ROOT4, cluster, 'sdss_frame_r_%06d_%d_%04d' % (run, camcol, field), 'frame-r-%06d-%d-%04d.fits.bz2' % (run, camcol, field)),
                psfield=psfield_path(cluster, run, camcol, field),
                catalog=os.path.join(ROOT2, cluster, 'hugs_' + cluster, 'hlsp_hugs_hst_wfc3-uvis-acs-wfc_%s_multi_v1_catalog-meth1.txt' % cluster),
                photofield=os.path.join(ROOT2, 'survey', 'sdss_photofield_%06d_%d' % (run, camcol), 'photoField-%06d-%d.fits' % (run, camcol)))


def cmd_build(args):
    c = args.cluster
    if args.frame:
        run, cam, fld = (int(v) for v in args.frame.split(':'))
    else:
        with open(args.survey) as f:
            r0 = {r['cluster']: r for r in parse_survey(f.read())}[c]
        run, cam, fld = r0['run'], r0['camcol'], r0['field']
    P = cluster_paths(c, run, cam, fld)
    fr, ps = S2.read_frame(P['frame']), S2.read_psfield(P['psfield'])
    cat = S2.read_hugs(P['catalog'])
    truth = S2.make_truth(cat, S2.world2pix_for(fr, 'astrans'))
    truth_h = S2.make_truth(cat, S2.world2pix_for(fr, 'header'))
    miss = S2.unmeasured_f606w(cat, S2.world2pix_for(fr, 'astrans'))
    print('cluster %s frame run %d camcol %d field %d (rerun 301); NMGY %.8f; DATE-OBS %s; psField gain %.4f dark variance %.4f' % (c, run, cam, fld, fr['header']['NMGY'], fr['header']['DATE-OBS'], ps['gain'], ps['dark_variance']))
    sky = S2.sky_image(fr['allsky'], fr['xinterp'], fr['yinterp'])
    dn_full = S2.frame_to_dn(fr['nmgy'], fr['cal'], sky)
    sat = S2._saturation_level(P['photofield'], fld)
    screened = S2.screen_patches(truth, dn_full, sat, missing=miss)
    cands = [dict(x0=q['x0'], y0=q['y0'], count=q['count'], density=q['density'], dn_max=q['dn_max']) for q in screened if all(q[k] for k in S2.RULES)]
    cands, n_removed = fwhm_filter(cands, ps['kl'])
    fun = S2.funnel(screened) + [('patch-centre stamp FWHM in [%g, %g] px' % (RULE['fwhm_lo'], RULE['fwhm_hi']), len(cands))]
    print('saturation level (photoField, r) %.0f DN; origins on the 10 px grid after each rule: %s' % (sat, '; '.join('%s %d' % f_ for f_ in fun)))
    print('candidate origins satisfying all rules: %d' % len(cands))
    for (lo, hi) in S2.CLASSES:
        print('   class [%.2f, %.2f]: %d candidates' % (lo, hi, sum(1 for q in cands if lo <= q['density'] <= hi)))
    chosen = S2.select_patches(cands)
    print('%d patches selected: %s' % (len(chosen), [(q['x0'], q['y0'], round(q['density'], 4)) for q in chosen]))
    if args.dry_run:
        return []
    os.makedirs(ROOT4, exist_ok=True)
    for suffix, t in (('', truth), ('_header', truth_h)):
        np.savez(os.path.join(ROOT4, 'truth_hugs_%s_meth1%s.npz' % (c, suffix)), **t)
    table, pad = [], 12
    for q in chosen:
        x0, y0 = q['x0'], q['y0']
        arrays, meta = S2.extract_patch(fr, ps, x0, y0, dn_full=dn_full, sky=sky)
        name = S2.patch_name(run, cam, fld, x0, y0)
        outdir = os.path.join(ROOT4, name)
        S2.write_patch(outdir, arrays, meta, fr['header'])
        t = S2.truth_for_patch(truth, x0, y0, S2.N_PATCH)
        files = ('patch_arrays.npz', 'patch_dn.fits', 'meta.json')
        fw = S2.psf_fwhm(arrays['psf51'])
        bc = S2.block_counts(truth, x0, y0)
        row = dict(name=name, x0=x0, y0=y0, n22=int((t['mag'] < 22).sum()), density22=float((t['mag'] < 22).sum()) / 1e4, n21=int((t['mag'] < 21).sum()), n23=int((t['mag'] < 23).sum()), n_all=int(len(t['mag'])),
                   n22_scored=int(((t['mag'] < 22) & t['inner']).sum()), brightest=float(t['mag'].min()), dn_max=float(arrays['dn'].max()), dn_median=float(np.median(arrays['dn'])),
                   sky_dn_median=float(np.median(arrays['sky_dn'])), fwhm_px=fw, crowding=S2.crowding(truth, x0, y0, fw),
                   peak_snr_at_22=S2.peak_snr(22.0, arrays['psf51'].max(), float(fr['header']['NMGY']), float(np.median(arrays['sky_dn'])), ps['gain'], ps['dark_variance']),
                   blocks_lt_24=dict(min=int(bc.min()), median=float(np.median(bc)), max=int(bc.max()), below_50=int((bc < 50).sum())),
                   sha256=dict(patch_arrays_content=S2.content_sha256(os.path.join(outdir, files[0])), patch_dn_fits=S2.file_sha256(os.path.join(outdir, files[1])), meta_json=S2.file_sha256(os.path.join(outdir, files[2]))))
        # DESCRIPTIVE ONLY (nothing is decided or corrected from it): the net rigid offset of the truth mapping against the image, by the PSF-convolved model of the truth (astro_second.rigid_shift)
        shifts = {}
        for label, tt in (('asTrans', truth), ('header', truth_h)):
            tp = S2.truth_for_patch(tt, x0 - pad, y0 - pad, S2.N_PATCH + 2 * pad)
            sx, sy, chi_min, chi_0 = S2.rigid_shift(arrays['dn'], ps['gain'], ps['dark_variance'], arrays['psf51'], tp['x'] - pad, tp['y'] - pad, tp['mag'])
            shifts[label] = dict(dx=sx, dy=sy, chi2_best=chi_min, chi2_zero=chi_0)
        row['image_shift_px'] = shifts
        table.append(row)
        print('patch %s x0 %d y0 %d F606W<22: %d (%.4f /px) <21: %d <23: %d all: %d scored<22: %d brightest %.2f DNmax %.0f FWHM %.2f px blocks F606W<24 min %d median %.0f max %d ; image shift asTrans (%+.3f, %+.3f) chi2 %.0f -> %.0f' % (
            name, x0, y0, row['n22'], row['density22'], row['n21'], row['n23'], row['n_all'], row['n22_scored'], row['brightest'], row['dn_max'], row['fwhm_px'], row['blocks_lt_24']['min'], row['blocks_lt_24']['median'],
            row['blocks_lt_24']['max'], shifts['asTrans']['dx'], shifts['asTrans']['dy'], shifts['asTrans']['chi2_zero'], shifts['asTrans']['chi2_best']), flush=True)
    out = os.path.join(ROOT4, 'fresh_%s_patches.json' % c)
    with open(out, 'w') as f:
        json.dump(dict(cluster=c, run=run, camcol=cam, field=fld, nmgy=float(fr['header']['NMGY']), date_obs=fr['header']['DATE-OBS'], saturation_level_dn=sat, n_candidates=len(cands), funnel=fun,
                       candidates_per_class=[sum(1 for q in cands if lo <= q['density'] <= hi) for lo, hi in S2.CLASSES], patches=table), f, indent=1)
    print('wrote', out)
    return table


def cmd_registry(args):
    """results/astro/fresh_patches.json: for each fresh cluster its frame, patches and file names (the patches' hashes are in the freeze record, not here)."""
    clusters = {}
    for c in args.clusters:
        with open(os.path.join(ROOT4, 'fresh_%s_patches.json' % c)) as f:
            d = json.load(f)
        title = {'ngc5272': 'M3 (NGC 5272)', 'ngc6205': 'M13 (NGC 6205)', 'ngc6254': 'M10 (NGC 6254)', 'ngc7078': 'M15 (NGC 7078)', 'ngc6934': 'NGC 6934', 'ngc5466': 'NGC 5466', 'ngc5053': 'NGC 5053'}.get(c, c)
        clusters[c] = dict(title='%s, SDSS run %d camcol %d field %d' % (title, d['run'], d['camcol'], d['field']), run=d['run'], camcol=d['camcol'], field=d['field'], patches=[[p['x0'], p['y0']] for p in d['patches']],
                           root=ROOT4, prefix=frame_prefix(d['run'], d['camcol'], d['field']), truth='truth_hugs_%s_meth1.npz' % c, densities=[round(p['density22'], 4) for p in d['patches']],
                           fwhm_px=[round(p['fwhm_px'], 4) for p in d['patches']], date_obs=d['date_obs'],
                           shifts=[[round(p['image_shift_px']['asTrans']['dx'], 4), round(p['image_shift_px']['asTrans']['dy'], 4)] for p in d['patches']])      # the net rigid offset of the truth mapping against the image (px to ADD to the truth), DESCRIPTIVE / the secondary scoring
    with open(args.out, 'w') as f:
        json.dump(dict(note='registry of the fresh clusters (experiments/astro_fresh.py registry); root is $ASTRO_DATA_4 at the time of writing', clusters=clusters), f, indent=1, sort_keys=True)
        f.write('\n')
    print('wrote', args.out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('select')
    p.add_argument('--survey', required=True)
    p.add_argument('--out')
    p = sub.add_parser('fields')
    p.add_argument('--clusters', nargs='*')
    p.add_argument('--out', required=True)
    p = sub.add_parser('select-frames')
    p.add_argument('--fields', required=True)
    p.add_argument('--out')
    p.add_argument('--list-psfield', action='store_true')
    p = sub.add_parser('build')
    p.add_argument('--cluster', required=True, choices=sorted(UNIVERSE))
    p.add_argument('--survey', default=os.path.join(REPO, 'results', 'astro', 'fresh_survey.txt'))
    p.add_argument('--frame', help='run:camcol:field of the frame (amendment A1: the chosen frame of the cluster); default: the survey\'s best frame')
    p.add_argument('--dry-run', action='store_true')
    p = sub.add_parser('registry')
    p.add_argument('--clusters', nargs='+', required=True)
    p.add_argument('--out', default=REGISTRY)
    a = ap.parse_args()
    {'select': cmd_select, 'fields': cmd_fields, 'select-frames': cmd_select_frames, 'build': cmd_build, 'registry': cmd_registry}[a.cmd](a)


if __name__ == '__main__':
    main()
