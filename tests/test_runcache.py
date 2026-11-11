"""The per-run cache of the long experiments (src/core/runcache.py): with EXP_CACHE set a job's result is stored and returned again without
recomputation, identical; without it the job runs as usual."""
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.core.runcache import cached                                    # noqa: E402

CALLS = []


@cached
def _job(args):
    CALLS.append(args)
    return {'x': args[0] * 2, 'tag': args[1]}


def test_the_cache_stores_and_returns_without_recomputing(tmp_path, monkeypatch):
    monkeypatch.delenv('EXP_CACHE', raising=False)
    CALLS.clear()
    assert _job((3, 'a')) == {'x': 6, 'tag': 'a'} and _job((3, 'a')) == {'x': 6, 'tag': 'a'} and len(CALLS) == 2
    monkeypatch.setenv('EXP_CACHE', str(tmp_path))
    CALLS.clear()
    r1 = _job((3, 'a')); r2 = _job((3, 'a')); r3 = _job((4, 'a'))
    assert r1 == r2 == {'x': 6, 'tag': 'a'} and r3 == {'x': 8, 'tag': 'a'} and len(CALLS) == 2
    assert len([f for f in os.listdir(tmp_path) if f.endswith('.pkl')]) == 2


def test_the_key_carries_the_codes_version(tmp_path, monkeypatch):
    # referee round 28: a result computed by other code is never returned — the key includes a hash of src/core and experiments
    from src.core import runcache
    monkeypatch.setenv('EXP_CACHE', str(tmp_path))
    CALLS.clear()
    _job((5, 'b'))
    monkeypatch.setattr(runcache, 'VERSION', runcache.VERSION + '-edited')
    _job((5, 'b'))
    assert len(CALLS) == 2 and len([f for f in os.listdir(tmp_path) if f.endswith('.pkl')]) == 2
    assert runcache.code_version() == runcache.code_version() and len(runcache.code_version()) == 40
