"""
Round theory3, direction D: the per-target critical mass of the readout runs (the raw behind the mismatch analysis of theory3d.md).

The 200 benchmark targets of Part C of `exp_parametrization.py` (seed 901) are re-run by plain gradient flow (the same cell, the same path as the stored reference) over a
grid of the start eps, to find for each target the eps below which the formed set is the sliding path's: depth 3 (D3, w = u^3 - v^3) on half-decades 10^-4 ... 10^-11 for the
118 targets of runner-up margin below 0.03, the paper's unit (P0) on 10^-5 ... 10^-9 for the 55 targets of margin below 0.01. A cell records: set_equal, order_equal, n_formed,
form_order, dormant displacement and maximum, exact, status, seconds.

Usage: python experiments/theory3d_scan.py [--run] [--readouts D3,P0] [--limit N]    -> results/synthetic/theory3d_scan_raw.npz (resumable; checkpoint every 4 targets)
Environment: T3D_OUT (output directory), OMP_NUM_THREADS=1 is set.
"""
import os
import sys
import time
import warnings

os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np

HERE = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..'))
import exp_parametrization as P      # noqa: E402

OUT = os.environ.get('T3D_OUT', os.path.join(HERE, '..', 'results', 'synthetic'))
RAW = os.path.join(OUT, 'theory3d_scan_raw.npz')
PARAM_RAW = os.environ.get('T3D_PARAM_RAW', os.path.join(HERE, '..', 'results', 'chaos', 'parametrization_raw.npz'))
SPEC = {'D3': dict(eps=tuple(10.0 ** (-e / 2.0) for e in range(8, 23)), margin_max=0.03),
        'P0': dict(eps=(1e-5, 1e-6, 1e-7, 1e-8, 1e-9), margin_max=0.01)}
KEEP = ('set_equal', 'order_equal', 'n_formed', 'form_order', 'dormant_disp', 'dormant_wmax', 'exact', 'status', 'secs')


def load_targets(param_raw=PARAM_RAW):
    """the 200 targets of Part C with their stored sliding paths: [dict(i, cs, As, ref)]"""
    z = np.load(param_raw, allow_pickle=True)
    return list(z['C'][0])


def load_scan(path=RAW):
    if not os.path.exists(path):
        return {}
    z = np.load(path, allow_pickle=True)
    return {k: z[k][0] for k in z.files}


def save_scan(path, res):
    tmp = path + '.tmp.npz'
    np.savez_compressed(tmp, **{k: np.array([v], dtype=object) for k, v in res.items()})
    os.replace(tmp, path)


def run_cell(name, t, eps, setup):
    L, K, c0 = setup
    f = K.synth(t['cs'][:, None], t['As'])
    ref = t['ref']
    tau_end = 1.05 * max(ref['tau']) + 1.0                # the horizon of Part C
    c = P.cell(P.make_readout(name), K, f, c0, t['cs'], t['As'], float(eps), ref, tau_end)
    return {k: c[k] for k in KEEP}


def scan(readouts=('D3', 'P0'), limit=None, path=RAW, log=True):
    warnings.filterwarnings('ignore')
    T = load_targets()
    setup = P.setup()
    res = load_scan(path)
    t0, done_since = time.time(), 0
    for name in readouts:
        sp = SPEC[name]
        tab = res.setdefault(name, {})                # {(target index, eps): cell}
        todo = [t for t in T if t['ref']['margin'] < sp['margin_max']]
        for t in (todo if limit is None else todo[:limit]):
            for eps in sp['eps']:
                key = (int(t['i']), float(eps))
                if key not in tab:
                    tab[key] = run_cell(name, t, eps, setup)
            done_since += 1
            if done_since % 4 == 0:
                save_scan(path, res)
                if log:
                    print(f"[{time.time() - t0:.0f}s] {name} target {t['i']}", flush=True)
        save_scan(path, res)
    return res


def main():
    if '--run' not in sys.argv:
        print(__doc__)
        return
    ro = tuple(sys.argv[sys.argv.index('--readouts') + 1].split(',')) if '--readouts' in sys.argv else ('D3', 'P0')
    lim = int(sys.argv[sys.argv.index('--limit') + 1]) if '--limit' in sys.argv else None
    os.makedirs(OUT, exist_ok=True)
    scan(ro, lim)


if __name__ == '__main__':
    main()
