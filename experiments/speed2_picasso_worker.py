"""
Picasso's cost per frame, one frame at a time (run by `exp_speed2.py` in an interpreter that can import picasso; single thread)
================================================================================================================================

    python -I speed2_picasso_worker.py <frames.npy (n, 64, 64) DN> <params.json> <out.npz> <repeats>

Identification + fit (`smlm_common._picasso_localize`, the call the ladder uses) of each frame separately after a JIT warm-up on 8 other frames; per repeat, per frame, the wall
seconds (perf_counter) and the process CPU seconds. Saves `wall`, `cpu` (repeats, n) and the number of localizations per frame.
"""
import json
import os
import sys
import time
import warnings

os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('NUMBA_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
HERE = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, HERE)
import numpy as np                      # noqa: E402
import smlm_common as C                 # noqa: E402


def main(frames_path, params_path, out_path, repeats):
    import picasso  # noqa: F401
    from picasso import localize as PL
    frames = C._worker_frames(np.load(frames_path))
    p = C._merge_params(json.load(open(params_path)))
    n = len(frames)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        C._picasso_localize(frames[:8], p, PL)
        wall, cpu, nloc = np.zeros((repeats, n)), np.zeros((repeats, n)), np.zeros(n, int)
        for r in range(repeats):
            for i in range(n):
                t0, c0 = time.perf_counter(), time.process_time()
                tab = C._picasso_localize(frames[i:i + 1], p, PL)
                wall[r, i], cpu[r, i] = time.perf_counter() - t0, time.process_time() - c0
                nloc[i] = len(tab['x'])
    np.savez(out_path, wall=wall, cpu=cpu, nloc=nloc, picasso_version=getattr(picasso, '__version__', 'unknown'))


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]))
