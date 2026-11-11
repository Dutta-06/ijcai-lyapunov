"""
A per-run cache for the long experiments: with the environment variable EXP_CACHE set to a directory, a job decorated with `cached` stores
its result there under a key made of its file, name, arguments and the code's version, and returns the stored result when the same job is
asked again — so an interrupted experiment restarts where it stopped. Without EXP_CACHE the jobs run as usual. The results are the same
either way (the jobs are deterministic in their arguments and the code). The code's version is a hash of every `src/core/*.py` and
`experiments/*.py` file, taken when this module is imported, so a result computed by other code is never returned (referee round 28).
"""
import functools
import glob
import hashlib
import os
import pickle

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))


def code_version():
    h = hashlib.sha1()
    for p in sorted(glob.glob(os.path.join(_ROOT, 'src', 'core', '*.py')) + glob.glob(os.path.join(_ROOT, 'experiments', '*.py'))):
        h.update(os.path.relpath(p, _ROOT).encode())
        with open(p, 'rb') as f:
            h.update(f.read())
    return h.hexdigest()


VERSION = code_version()


def cached(fn):
    @functools.wraps(fn)
    def wrapper(args):
        d = os.environ.get('EXP_CACHE')
        if not d:
            return fn(args)
        src = os.path.splitext(os.path.basename(fn.__code__.co_filename))[0]      # the same whether run as a script or imported
        key = hashlib.sha1(repr((src, fn.__name__, args, VERSION)).encode()).hexdigest()
        path = os.path.join(d, f'{src}.{fn.__name__}_{key}.pkl')
        if os.path.exists(path):
            with open(path, 'rb') as f:
                return pickle.load(f)
        out = fn(args)
        os.makedirs(d, exist_ok=True)
        tmp = f'{path}.tmp{os.getpid()}'
        with open(tmp, 'wb') as f:
            pickle.dump(out, f)
        os.replace(tmp, path)
        return out
    return wrapper
