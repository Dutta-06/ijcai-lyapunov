"""
The astronomy path with the original (guarded) or the fast evaluation, through the frozen `exp_astro_m2.run_path`, nothing frozen modified
=============================================================================================================================================

`exp_astro_m2.run_path` calls `imaging.sliding_path` (looked up at call time). As `exp_astro_cyclefix.py` does, `engine(name)` replaces that function in THIS process for the duration of a
`with` block: 'old' by `imaging_guard.sliding_path_guarded` (the corrected study's path), 'new' by `imaging_fast.sliding_path_fast` on a `FastLattice` of the same lattice (guard on). The last path's
events, entry times and barred points are kept in `TRACE`. `cap(n)` sets the entry cap of the path (`exp_astro_m2.MAX_ENTRIES`) for a `with` block (the truncated paths of the speed study).
"""
import contextlib
import os
import sys

HERE = os.path.abspath(os.path.dirname(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)
import exp_astro_m2 as X                    # noqa: E402
from src.core import imaging as IM          # noqa: E402
from src.core import imaging_fast as F      # noqa: E402
from src.core import imaging_guard as G     # noqa: E402

TRACE = {}


def _path_old(*args, **kw):
    res = G.sliding_path_guarded(*args, **kw)
    TRACE.update(events=res['events'], taus=res['taus'], barred=res['barred'])
    return res


def _path_old_tight(*args, **kw):
    """The original guarded path with the block-fit tolerance 10 times smaller (the control of the speed study: how far the original's own result moves with the iteration tolerance)."""
    return _path_old(*args, **dict(kw, tol=kw['tol'] / 10.0))


def _path_new(psf, y, sw, lat, *args, **kw):
    res = F.sliding_path_fast(psf, y, sw, F.FastLattice(lat.psf, lat.h, lat.margin), *args, **kw)
    TRACE.update(events=res['events'], taus=res['taus'], barred=res['barred'])
    return res


@contextlib.contextmanager
def engine(name):
    saved = IM.sliding_path
    IM.sliding_path = {'old': _path_old, 'new': _path_new, 'old_tight': _path_old_tight}[name]
    try:
        yield
    finally:
        IM.sliding_path = saved


@contextlib.contextmanager
def cap(n):
    saved = X.MAX_ENTRIES
    X.MAX_ENTRIES = int(n)
    try:
        yield
    finally:
        X.MAX_ENTRIES = saved
