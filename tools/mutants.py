"""Single-point mutation testing of one experiment script: each mutant against the script's tests (pytest -x) and against its
`--from-raw` regeneration of the md (byte comparison with the repository's md). Nothing runs in the repository: each slot is a private
copy of src/, experiments/, the tests file and the raw under the work directory.

    venv/bin/python tools/mutants.py --exp exp_width_mismatch --tests tests/test_width_mismatch.py \\
        --raw results/chaos/width_mismatch_raw.npz --md results/chaos/width_mismatch.md --out-env WM_OUT --workers-env WM_WORKERS \\
        --spec tools/mutants/round84_D58_own.py --work <scratch dir> [--slot a --i0 0 --step 3] [--only ID,ID] [-k EXPR]

A spec file defines M (or MUTANTS), a list of tuples (id, old, new), (id, what, old, new) or (id, what, old, new, k): `old` replaced by
`new` at its k-th occurrence (from 0) in experiments/<exp>.py; an id whose `old` is not found is recorded as 'anchor not found'. Run
several slots side by side with --i0 0/1/2 --step 3 (distinct --slot names). Results, one JSON line per mutant, in <work>/res_<slot>.jsonl:
tests (True: the tests kill it), md (True: the md changes), fr_rc (`--from-raw`'s exit code), secs, last (pytest's last line). A mutant
that neither kills is a survivor: kill it with a test, or record why it is equivalent on the data.

The kill criterion is the tests file alone (`pytest -x` on --tests, in the slot): not any check the experiment script makes itself, and not the md
column (which is reported, not part of the criterion). Every other repository file that the tests read has to be copied into the slot with --extra
(relative paths; a .npz is also copied next to the raw into the slot's output directory): mutating an experiment script whose tests read other
files fails even the unmutated control without them. For experiments/exp_parametrization.py (tests/test_parametrization.py) these are
`--extra results/chaos/transit_certificate_raw.npz --extra results/chaos/parametrization.md` (the target sample reads the first; the report test
compares the md with the raw's). A control run is a spec with a single mutant whose `old` is '' (the source is left unmutated): it has to pass."""
import argparse
import importlib.util
import json
import os
import shutil
import subprocess
import time

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
ap = argparse.ArgumentParser()
for a in ('--exp', '--tests', '--raw', '--md', '--out-env', '--spec', '--work'):
    ap.add_argument(a, required=True)
ap.add_argument('--workers-env')
ap.add_argument('--slot', default='a')
ap.add_argument('--i0', type=int, default=0)
ap.add_argument('--step', type=int, default=1)
ap.add_argument('--only')
ap.add_argument('-k', dest='kexpr')
ap.add_argument('--extra', action='append', default=[], help='more repository files the tests read (relative paths), copied into the slot')
ap.add_argument('--timeout', type=int, default=1800)
A = ap.parse_args()

spec_ = importlib.util.spec_from_file_location('mutant_spec', A.spec)
mod = importlib.util.module_from_spec(spec_)
spec_.loader.exec_module(mod)
LIST = []
for m in getattr(mod, 'MUTANTS', None) or mod.M:
    if len(m) == 3:
        LIST.append((m[0], '', m[1], m[2], 0))
    elif len(m) == 4:
        LIST.append((m[0], m[1], m[2], m[3], 0))
    else:
        LIST.append(tuple(m))
only = set(A.only.split(',')) if A.only else None

S = os.path.join(os.path.abspath(A.work), f'slot_{A.slot}')
shutil.rmtree(S, ignore_errors=True)
shutil.copytree(os.path.join(REPO, 'src'), os.path.join(S, 'src'), ignore=shutil.ignore_patterns('__pycache__'))
os.makedirs(os.path.join(S, 'experiments'))
for f in os.listdir(os.path.join(REPO, 'experiments')):
    if f.endswith('.py'):
        shutil.copy(os.path.join(REPO, 'experiments', f), os.path.join(S, 'experiments', f))
for rel in [A.tests, A.raw] + A.extra:
    os.makedirs(os.path.dirname(os.path.join(S, rel)), exist_ok=True)
    shutil.copy(os.path.join(REPO, rel), os.path.join(S, rel))
OUTD = os.path.join(S, 'out')
os.makedirs(OUTD)
shutil.copy(os.path.join(REPO, A.raw), OUTD)
for rel in A.extra:
    if rel.endswith('.npz'):
        shutil.copy(os.path.join(REPO, rel), OUTD)
src_path = os.path.join(S, 'experiments', A.exp + '.py')
orig = open(src_path).read()
md0 = open(os.path.join(REPO, A.md)).read()
md_out = os.path.join(OUTD, os.path.basename(A.md))
env = dict(os.environ, OMP_NUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1', **{A.out_env: OUTD})
if A.workers_env:
    env[A.workers_env] = '1'
PY = os.path.join(REPO, 'venv', 'bin', 'python')
res_path = os.path.join(os.path.abspath(A.work), f'res_{A.slot}.jsonl')
for mid, what, old, new, k in LIST[A.i0::A.step]:
    if only is not None and mid not in only:
        continue
    rec = dict(id=mid, what=what)
    idx = -1
    for _ in range(k + 1):
        idx = orig.find(old, idx + 1) if old else 0
        if idx < 0:
            break
    if idx < 0:
        rec['error'] = 'anchor not found'
        open(res_path, 'a').write(json.dumps(rec) + '\n')
        print(mid, 'ANCHOR NOT FOUND', flush=True)
        continue
    open(src_path, 'w').write(orig[:idx] + new + orig[idx + len(old):] if old else orig)
    open(md_out, 'w').write(md0)                                        # the tests see the repository's md (a test may compare it with the raw)
    t0 = time.time()
    cmd = ['nice', '-n', '10', PY, '-m', 'pytest', '-x', '-q', '-p', 'no:cacheprovider', A.tests] + (['-k', A.kexpr] if A.kexpr else [])
    try:
        p = subprocess.run(cmd, cwd=S, env=env, capture_output=True, text=True, timeout=A.timeout)
        tests, last = p.returncode != 0, (p.stdout.strip().splitlines() or [''])[-1][:150]
    except subprocess.TimeoutExpired:
        tests, last = True, 'timeout'
    if os.path.exists(md_out):
        os.remove(md_out)
    q = subprocess.run(['nice', '-n', '10', PY, os.path.join('experiments', A.exp + '.py'), '--from-raw'], cwd=S, env=env,
                       capture_output=True, text=True, timeout=A.timeout)
    md = open(md_out).read() if os.path.exists(md_out) else None
    rec.update(tests=tests, md=(md != md0), fr_rc=q.returncode, secs=round(time.time() - t0, 1), last=last)
    open(res_path, 'a').write(json.dumps(rec) + '\n')
    print(mid, 'tests' if tests else '-', 'md' if rec['md'] else '-', rec['secs'], flush=True)
open(src_path, 'w').write(orig)
