"""
Does the order of entry matter? Close pairs of emitters, the accumulated dual against the instantaneous correlation (WP7.2)
=========================================================================================================================

On frames of the 2016 challenge's EM-CCD model with a known Gaussian PSF (the model is the truth: no halo, no mismatch, so that only the ORDER can differ) a scene is one pair of
emitters at a separation d (in PSF widths) and a brightness ratio r, at a random place and direction, at two photon levels of the brighter emitter. The sliding inverse scale space (the
accumulated dual enters first, `order='dual'`) and matching pursuit with sliding (the instantaneous correlation, `order='greedy'`) fit the same frames with the same model, the same
blocks and the same stop (BIC, nothing to tune; the significance stop gives the same picture, see --stop). A pair counts as RESOLVED when exactly two units are returned and each lies
within 0.7 width of one emitter (one-to-one); otherwise it is merged (one unit), split (three or more) or mislocated. The table is the resolved fraction with a Wilson interval, the
paired difference of the two orders over the same scenes with a bootstrap interval, and the entry orders that differ.

    python experiments/exp_order_pairs.py --run [--workers 3] [--scenes 100]       (writes results/smlm/order_pairs_raw.npz)
    python experiments/exp_order_pairs.py --from-raw                                (regenerates results/smlm/order_pairs.md)
"""
import argparse
import os
import sys
import time

os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np

HERE = os.path.abspath(os.path.dirname(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, REPO)
from src.core import imaging as IM      # noqa: E402

OUT = os.environ.get('OP_OUT', os.path.join(REPO, 'results', 'smlm'))
RAW, MD = os.path.join(OUT, 'order_pairs_raw.npz'), os.path.join(OUT, 'order_pairs.md')
TAGGED = 'order_pairs_%s_raw.npz'                       # a second stop rule is stored beside the first (--tag)
SIG, NPIX, BG, VAR_A, VAR_B0 = 1.4, 32, 100.0, 2.236, 0.078
SEPS = (1.2, 1.6, 2.0, 2.5, 3.0)           # separations in widths
RATIOS = (1.0, 3.0, 10.0)                  # brightness ratios
PHOTONS = (2000.0, 6000.0)                 # photons of the brighter emitter
TOL = 0.7                                  # widths
ORDERS = ('dual', 'greedy')


def scene(rng, d, ratio, ph):
    """Positions (2, 2) in px and photons (2): the pair's centre uniform in [11, 21]^2, a random direction, the brighter emitter first."""
    c = rng.uniform(11.0, 21.0, 2)
    th = rng.uniform(0, 2 * np.pi)
    off = 0.5 * d * SIG * np.array([np.cos(th), np.sin(th)])
    return np.stack([c - off, c + off]), np.array([ph, ph / ratio])


def classify(c, w, pos):
    """0 mislocated, 1 resolved (two units, each within TOL widths of one emitter), 2 merged (one unit), 3 split (three or more units)."""
    n = len(w)
    if n == 1:
        return 2
    if n >= 3:
        return 3
    if n == 0:
        return 0
    d = np.linalg.norm(c[:, None, :2] - pos[None, :, :], axis=2)
    ok = (d[0, 0] <= TOL * SIG and d[1, 1] <= TOL * SIG) or (d[0, 1] <= TOL * SIG and d[1, 0] <= TOL * SIG)
    return 1 if ok else 0


_G = {}


def _init():
    P = IM.PixelPSF(NPIX, NPIX, SIG, hp=5.0)
    _G.update(P=P, Pp=IM.PixelPSF(NPIX, NPIX, SIG), L=IM.Lattice(P, h=0.5, margin=1.0))


def trial(args):
    """One scene (d, ratio, photons index, seed): the class under each order and the number of units."""
    d, ratio, ip, seed, stop, clean = args
    rng = np.random.default_rng(seed)
    pos, w = scene(rng, d, ratio, PHOTONS[ip])
    frame = IM.simulate_emccd(_G['Pp'], pos, w, BG, rng)
    out = []
    for order in ORDERS:
        r = IM.fit_frame(frame, _G['P'], _G['L'], BG, VAR_A, VAR_B0, widths=(1.0, 3.2), stop=stop, max_entries=8, order=order, clean=clean)
        out.append((classify(r['c'], r['w'], pos), len(r['w'])))
    return out


def run(workers, scenes, stop, out_path=RAW, clean=None):
    from multiprocessing import Pool
    jobs = [(d, r, ip, 1000 * (10 * i + j) + s, stop, clean) for i, d in enumerate(SEPS) for j, r in enumerate(RATIOS) for ip in range(len(PHOTONS)) for s in range(scenes)]
    t0 = time.time()
    with Pool(workers, initializer=_init) as pool:
        res = pool.map(trial, jobs, chunksize=4)
    cls = np.array([[o[0] for o in r] for r in res]).reshape(len(SEPS), len(RATIOS), len(PHOTONS), scenes, len(ORDERS))
    nu = np.array([[o[1] for o in r] for r in res]).reshape(cls.shape)
    np.savez_compressed(out_path, cls=cls, nu=nu, seps=np.array(SEPS), ratios=np.array(RATIOS), photons=np.array(PHOTONS), scenes=scenes, stop=np.array([str(stop[0]), float(stop[1])], dtype=object).astype(str), clean=np.array([np.nan if clean is None else float(clean)]))
    print('%d scenes x %d orders in %.0f s' % (len(jobs), len(ORDERS), time.time() - t0), flush=True)


def wilson(k, n, z=1.96):
    if n == 0:
        return float('nan'), float('nan')
    p = k / n
    den = 1 + z * z / n
    mid = (p + z * z / (2 * n)) / den
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return mid - half, mid + half


def report(path=RAW, boot=2000):
    z = np.load(path, allow_pickle=False)
    cls, seps, ratios, photons = z['cls'], z['seps'], z['ratios'], z['photons']
    stop = tuple(z['stop'])
    clean = float(z['clean'][0]) if 'clean' in z.files else float('nan')
    n = cls.shape[3]
    L = ['# Does the order of entry matter? Close pairs, the accumulated dual against the instantaneous correlation', '',
         'Generated by `experiments/exp_order_pairs.py` from `order_pairs_raw.npz`; every number below is printed by the script. Setting: the script\'s docstring. %d scenes per cell, '
         'stop %s %s%s; a pair is resolved when exactly two units are returned, each within %g widths of one emitter. Intervals: 95 %% Wilson (proportions) and 95 %% bootstrap over scenes '
         '(%d resamples, seed 0; paired differences).' % (n, stop[0], stop[1], '' if np.isnan(clean) else ', then the units of own significance below %g discarded' % clean, TOL, boot), '']
    rng = np.random.default_rng(0)
    idx = rng.integers(0, n, size=(boot, n))
    tot = np.zeros((len(ORDERS), 4))
    rows = []
    for i, d in enumerate(seps):
        for j, r in enumerate(ratios):
            for k, ph in enumerate(photons):
                c = cls[i, j, k]                                  # (scenes, orders)
                res = (c == 1)
                cells = []
                for o in range(len(ORDERS)):
                    lo, hi = wilson(int(res[:, o].sum()), n)
                    cells.append('%.0f %% [%.0f, %.0f]' % (100 * res[:, o].mean(), 100 * lo, 100 * hi))
                diff = res[:, 0].astype(float) - res[:, 1].astype(float)
                bs = diff[idx].mean(1)
                lo, hi = np.percentile(bs, [2.5, 97.5])
                differ = int((res[:, 0] != res[:, 1]).sum())
                rows.append(['%.1f' % d, '%g' % r, '%g' % ph] + cells + ['%+.1f [%+.1f, %+.1f]' % (100 * diff.mean(), 100 * lo, 100 * hi), str(differ)])
                for o in range(len(ORDERS)):
                    for q in range(4):
                        tot[o, q] += int((c[:, o] == q).sum())
    L.append('| separation (widths) | flux ratio | photons (brighter) | resolved, accumulated dual (ISS) | resolved, instantaneous correlation | ISS minus greedy, points | scenes where they differ |')
    L.append('|---|---|---|---|---|---|---|')
    L += ['| ' + ' | '.join(r) + ' |' for r in rows]
    L += ['', 'All %d scenes: outcome counts (mislocated, resolved, merged, split) by order.' % (cls.shape[0] * cls.shape[1] * cls.shape[2] * n), '']
    L.append('| order | mislocated | resolved | merged | split |')
    L.append('|---|---|---|---|---|')
    for o, name in enumerate(('accumulated dual (ISS)', 'instantaneous correlation (greedy)')):
        L.append('| %s | %d | %d | %d | %d |' % (name, *tot[o]))
    resolved = cls == 1
    d_all = resolved[..., 0].astype(float) - resolved[..., 1].astype(float)
    flat = d_all.reshape(-1, n)
    bs = flat[:, idx].mean(axis=(0, 2))                                       # the scenes resampled jointly over all cells
    lo, hi = np.percentile(bs, [2.5, 97.5])
    L += ['', 'Over every cell the resolved fraction of the accumulated dual minus that of the instantaneous correlation is %+.2f points [%+.2f, %+.2f] (paired over scenes, %d scenes).' % (
        100 * d_all.mean(), 100 * lo, 100 * hi, d_all.size), '']
    return '\n'.join(L) + '\n'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', action='store_true')
    ap.add_argument('--from-raw', action='store_true')
    ap.add_argument('--workers', type=int, default=int(os.environ.get('OP_WORKERS', '3')))
    ap.add_argument('--scenes', type=int, default=100)
    ap.add_argument('--stop', default='bic:4', help="the stop rule, kind:parameter (bic:4 or snr:5)")
    ap.add_argument('--clean', type=float, default=None, help='after the stop, discard the units whose own significance is below this value (imaging.unit_significance)')
    ap.add_argument('--tag', default=None, help='store the raw / md under this tag (order_pairs_<tag>_raw.npz, order_pairs_<tag>.md) instead of the default names')
    A = ap.parse_args()
    kind, par = A.stop.split(':')
    raw, md = (RAW, MD) if A.tag is None else (os.path.join(OUT, TAGGED % A.tag), os.path.join(OUT, 'order_pairs_%s.md' % A.tag))
    if A.run:
        run(A.workers, A.scenes, (kind, float(par)), out_path=raw, clean=A.clean)
    if A.from_raw:
        with open(md, 'w') as f:
            f.write(report(raw))
        print('wrote', md)


if __name__ == '__main__':
    main()
