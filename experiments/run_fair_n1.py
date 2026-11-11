"""Driver: run every fair kernel baseline for the given seeds, one process per seed in parallel."""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
METHODS = ['proposed', 'ntk', 'sparse_gp', 'exact_gp', 'deep_gp']

if __name__ == "__main__":
    seeds = [int(s) for s in sys.argv[1:]] or [0, 1, 2]
    env = dict(os.environ, PYTHONUTF8="1", OMP_NUM_THREADS="2", OPENBLAS_NUM_THREADS="2", MKL_NUM_THREADS="2")
    procs = []
    for s in seeds:
        code = (f"import subprocess,sys; [subprocess.run([sys.executable, r'{os.path.join(HERE, 'kernel_baselines_fair.py')}', m, '{s}'], check=False) "
                f"for m in {METHODS!r}]")
        log = open(os.path.join(HERE, '..', 'results', f'fair_seed{s}.log'), 'w')
        procs.append(subprocess.Popen([sys.executable, '-c', code], stdout=log, stderr=subprocess.STDOUT, env=env))
    for p in procs:
        p.wait()
    print("all seeds done")
