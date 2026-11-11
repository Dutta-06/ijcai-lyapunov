"""Merge parts of an experiment's raw from run directories into the repository's raw (each run writes its own raw in its *_OUT dir, so
that no two processes write the same file):

    venv/bin/python tools/merge_parts.py results/chaos/<name>_raw.npz <part>=<run dir>/<name>_raw.npz [...]

Then regenerate the md with the experiment's `--from-raw`, twice, and check that the two agree byte for byte."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'experiments'))
import exp_cohesion_leak as E                                           # noqa: E402  (the raws' load/save)

raw = sys.argv[1]
res = E.load(raw)
for arg in sys.argv[2:]:
    part, src = arg.split('=', 1)
    r = E.load(src)
    if part not in r:
        raise SystemExit(f"{part} is not in {src}")
    res[part] = r[part]
    print(part, 'from', src, len(r[part]) if hasattr(r[part], '__len__') else '')
E.save(raw, **res)
print('saved', sorted(res))
