"""
Selection and construction of the held-out clusters of the cluster-dependence test (geometry only: no catalogue of any method is built, no detection is compared with a truth source)
================================================================================================================================================================================

The rule (with the amendment BRIGHT_LIMIT, DN_MARGIN below) is declared in `results/astro/clusterdep_plan.md` (section 3) and implemented here by reusing `astro_fresh.py` and `astro_second.py` (hash-frozen, imported, never edited) with a different rule
dictionary. Universe: the clusters of the coverage survey that no earlier test used (not M2, M53, M92, M3, M15). Frame gates (per frame): PSF status <= 3 (any), at least 0 PSF stars, FWHM of the
psField stamp at the frame centre in [2.0, 3.8] px, at least 30 survey windows (100 x 100 px origins of moderate density 0.03 to 0.15 per px inside the HUGS footprint). The clusters are ordered by the
windows of their best passing frame; walking that order, a cluster's passing frames are tried in the order of `astro_fresh.apply_rule_frames` (most windows first) and the first that leaves at least one
patch under the patch rules of M53 and M92 (`astro_second.py`, unchanged, plus the patch-centre FWHM gate [2.0, 3.8]) is the cluster's frame; the first six clusters with a patch are taken; up to two patches per
cluster (the first two of `astro_second.select_patches`: the scarcest density class first). Needs the domain environment (astropy).

    python -I experiments/clusterdep_select.py select [--out results/astro/clusterdep_selection.txt]   dry runs (the frames are fetched with data/astro/fetch_clusterdep.py first)
    python -I experiments/clusterdep_select.py build --out results/astro/clusterdep_patches.json       writes the patches of the chosen clusters and the registry
"""
import argparse
import contextlib
import io
import json
import os
import sys

import numpy as np

HERE = os.path.abspath(os.path.dirname(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
import astro_fresh as F                   # noqa: E402
import astro_second as S2                 # noqa: E402

UNIVERSE = ('ngc5053', 'ngc5466', 'ngc5904', 'ngc6205', 'ngc6254', 'ngc6535', 'ngc6838', 'ngc6934', 'ngc6981')
RULE = dict(status_max=3, nstar_min=0, fwhm_lo=2.0, fwhm_hi=3.8, windows_min=30, centre=(1023.5, 744.0), min_patches=1, n_fresh=6)
PATCHES_PER_CLUSTER = 2
BRIGHT_LIMIT, DN_MARGIN = 13.0, 6                                          # amendment of the plan (section 3.3): HST-bright F606W limit and the margin (px) of the bright-star and saturation rules; M53/M92/M3/M15: 14.5 and 10
FIELDS = os.path.join(REPO, 'results', 'astro', 'fresh_fields.json')
TITLES = {'ngc5053': 'NGC 5053', 'ngc5466': 'NGC 5466', 'ngc5904': 'M5 (NGC 5904)', 'ngc6205': 'M13 (NGC 6205)', 'ngc6254': 'M10 (NGC 6254)', 'ngc6535': 'NGC 6535', 'ngc6838': 'M71 (NGC 6838)', 'ngc6934': 'NGC 6934'}


def install_rule():
    """The relaxed constants into `astro_fresh` (its `fwhm_filter` reads `RULE` at call time)."""
    F.RULE.clear()
    F.RULE.update(RULE)
    S2.BRIGHT_LIMIT, S2.DN_MARGIN = BRIGHT_LIMIT, DN_MARGIN


def passing_frames():
    """(ranked clusters, every frame with its flags): `astro_fresh.apply_rule_frames` with this module's rule and universe, plus, per cluster, ALL its passing frames in the order of that rule."""
    with open(FIELDS) as f:
        fields = json.load(f)['fields']
    fwhm = {}
    for r in F.need_psfield(fields, RULE, UNIVERSE):
        p = F.psfield_path(r['cluster'], r['run'], r['camcol'], r['field'])
        if os.path.exists(p):
            fwhm[F.frame_id(r)] = F.stamp_fwhm(S2.read_psfield(p)['kl'], *RULE['centre'])
    ranked, allf = F.apply_rule_frames(fields, fwhm, RULE, UNIVERSE)
    per = {}
    for r in allf:
        if r['passes']:
            per.setdefault(r['cluster'], []).append(r)
    for c in per:
        per[c].sort(key=lambda r: (-r['windows'], r['seeing'], r['run'], r['camcol'], r['field']))
    return ranked, allf, per


def try_frame(r, dry=True):
    """Runs `astro_fresh.cmd_build` on the frame r (dry run or real) and returns (captured text, number of patches)."""
    ns = argparse.Namespace(cluster=r['cluster'], frame='%d:%d:%d' % (r['run'], r['camcol'], r['field']), survey=None, dry_run=dry)
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        out = F.cmd_build(ns)
    txt = buf.getvalue()
    n = 0
    for line in txt.splitlines():                                           # "<n> patches selected: [...]"
        if ' patches selected:' in line:
            n = int(line.split(' patches selected:')[0].strip())
    return txt, n


def cmd_select(args):
    install_rule()
    ranked, allf, per = passing_frames()
    L = ['Selection rule of clusterdep_plan.md section 3, geometry only (dry runs of astro_fresh.cmd_build with the relaxed rule)', 'rule constants: %s' % json.dumps(RULE), 'universe: %s' % ', '.join(UNIVERSE), '',
         '%-8s %-18s %7s %6s %7s %-10s %6s %6s %9s | %s' % ('cluster', 'frame (run cam fld)', 'seeing"', 'foot', 'windows', 'density', 'status', 'nstar', 'FWHM px', 'passes')]
    for r in sorted(allf, key=lambda r: (int(r['cluster'][3:]), -r['windows'], r['run'], r['camcol'], r['field'])):
        L.append('%-8s %-18s %7.3f %6.2f %7d %-10s %6d %6d %9s | %s' % (r['cluster'], '%d %d %d' % (r['run'], r['camcol'], r['field']), r['seeing'], r['frac'], r['windows'],
                 ('%.2f-%.2f' % (r['dmin'], r['dmax'])) if r['dmin'] is not None else '-', r['status'], r['nstar'], '%.3f' % r['fwhm'] if r['fwhm'] is not None else 'n/a', 'PASS' if r['passes'] else 'no'))
    L += ['', 'rank order of the clusters (windows of the best passing frame): ' + ', '.join('%s (%d)' % (r['cluster'], r['windows']) for r in ranked), '']
    chosen = []
    for r0 in ranked:
        if len(chosen) >= RULE['n_fresh']:
            break
        for r in per[r0['cluster']]:
            txt, n = try_frame(r, dry=True)
            L.append('--- dry run %s:%d:%d:%d (%d windows)' % (r['cluster'], r['run'], r['camcol'], r['field'], r['windows']))
            L += ['    ' + ln for ln in txt.splitlines() if ln.startswith('saturation') or ln.startswith('candidate') or 'class' in ln or 'patches selected' in ln]
            if n >= RULE['min_patches']:
                chosen.append((r, n))
                break
    L += ['', 'chosen clusters (frame, patches available): ' + '; '.join('%s %d:%d:%d (%d)' % (r['cluster'], r['run'], r['camcol'], r['field'], n) for r, n in chosen)]
    text = '\n'.join(L) + '\n'
    print(text)
    if args.out:
        with open(args.out, 'w') as f:
            f.write(text)


def cmd_build(args):
    install_rule()
    ranked, allf, per = passing_frames()
    with open(args.selection) as f:
        sel = [ln for ln in f.read().splitlines() if ln.startswith('chosen clusters')][0]
    names = [t.strip() for t in sel.split(': ', 1)[1].split(';')]
    clusters = {}
    for item in names:
        c, fr = item.split()[0], item.split()[1]
        run, cam, fld = (int(v) for v in fr.split(':'))
        r = [x for x in per[c] if (x['run'], x['camcol'], x['field']) == (run, cam, fld)][0]
        txt, n = try_frame(r, dry=False)
        with open(os.path.join(F.ROOT4, 'fresh_%s_patches.json' % c)) as f:
            d = json.load(f)
        keep = d['patches'][:PATCHES_PER_CLUSTER]
        clusters[c] = dict(title='%s, SDSS run %d camcol %d field %d' % (TITLES[c], d['run'], d['camcol'], d['field']), run=d['run'], camcol=d['camcol'], field=d['field'], patches=[[p['x0'], p['y0']] for p in keep],
                           root=F.ROOT4, prefix=F.frame_prefix(d['run'], d['camcol'], d['field']), truth='truth_hugs_%s_meth1.npz' % c, densities=[round(p['density22'], 4) for p in keep],
                           fwhm_px=[round(p['fwhm_px'], 4) for p in keep], date_obs=d['date_obs'], status=r['status'], nstar=r['nstar'], windows=r['windows'],
                           shifts=[[round(p['image_shift_px']['asTrans']['dx'], 4), round(p['image_shift_px']['asTrans']['dy'], 4)] for p in keep],
                           patches_built=[[p['x0'], p['y0']] for p in d['patches']])
        print(txt)
    with open(args.out, 'w') as f:
        json.dump(dict(note='registry of the held-out clusters of the cluster-dependence test (experiments/clusterdep_select.py build); root is $ASTRO_DATA_4 at the time of writing', clusters=clusters), f, indent=1, sort_keys=True)
        f.write('\n')
    print('wrote', args.out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('select')
    p.add_argument('--out')
    p = sub.add_parser('build')
    p.add_argument('--selection', default=os.path.join(REPO, 'results', 'astro', 'clusterdep_selection.txt'))
    p.add_argument('--out', required=True)
    a = ap.parse_args()
    {'select': cmd_select, 'build': cmd_build}[a.cmd](a)


if __name__ == '__main__':
    main()
