"""The cycle guard of the sliding path (src/core/imaging_guard.py): the enter/exit limit cycle of `imaging.sliding_path` (a lattice kernel enters on the correlation at its lattice start,
a unit leaves on the correlation at its slid centre, the exit puts the dual back on the wall, so the kernel re-enters at dt = 0, at the same entry time, for ever) is planted by replacing
`imaging.block_fit` with a stub, and the guard is held to its definition: a lattice point that has exited twice at the same entry time is barred for the rest of the path; a path without a
double exit at one entry time is the unguarded path bit for bit; `cycle_guard=False` is `imaging.sliding_path` always; the source of the two functions differs only in the guard's lines."""
import collections
import difflib
import inspect
import os
import re
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from src.core import imaging as IM            # noqa: E402
from src.core import imaging_guard as G       # noqa: E402

KW = dict(max_entries=8, tol=1e-8, passes=1, rtrunc=4.0)
HOLE = (30.0, 30.0)


def _world(seed=3, n=40, sigma=1.3, stars=((12.3, 11.8, 5000.0),), noise=0.5, hole=None):
    """A frame with Gaussian-PSF stars, an optional negative blob (a region where the correlation of the residual with a kernel is negative) and a little noise; (psf, y, sw, lattice)."""
    psf = IM.PixelPSF(n, n, sigma)
    rng = np.random.default_rng(seed)
    y = psf.synth(np.array([[s[0], s[1]] for s in stars]), np.array([s[2] for s in stars]))
    if hole is not None:
        y = y + psf.synth(np.array([[hole[0], hole[1]]]), np.array([hole[2]]))
    y = y + noise * rng.standard_normal((n, n))
    return psf, y, np.full((n, n), 1.0 / noise), IM.Lattice(psf, h=0.5, margin=1.0)


def _planted_world(second=None, noise=0.5, unit_weights=False):
    """One star (and, if `second` is given, another one of nearly the same brightness 7 px away, whose lattice kernel enters shortly after the first's) and a negative blob far from them;
    the first lattice kernel to enter is the one nearest the first star (BAD); returns (psf, y, sw, lat, BAD id, BAD's lattice start)."""
    stars = ((12.3, 11.8, 5000.0),) + ((second,) if second else ())
    psf, y, sw, lat = _world(hole=HOLE + (-3000.0,), stars=stars, noise=noise)
    if unit_weights:
        sw = np.ones_like(sw)
    bad = IM.sliding_path(psf, y, sw, lat, **dict(KW, max_entries=1))['events'][0][1]
    return psf, y, sw, lat, bad, lat.pos(*bad)


REX = pytest.mark.parametrize('r_ex', [0.0, 0.2, 1.2], ids=['no_zone', 'zone_0.2', 'zone_1.2'])           # the exclusion zone adds to `blocked` in its own line: the bar has to survive that line


def _plant(monkeypatch, starts, when=None):
    """Replaces `imaging.block_fit` (both path functions look it up at call time) by a stub: whenever the block holds a unit whose lattice start is in `starts` (the k-th such call for k in `when`,
    or every one if `when` is None) the stub returns that unit's coefficient 0 with its centre at the negative blob, where the correlation of the residual is not positive, so the unit
    leaves; every other call is the real fit. Returns the list that counts the calls that held a planted unit (clear it between two runs)."""
    real = IM.block_fit
    starts = np.atleast_2d(starts)
    held = []

    def stub(psf, y, sw, c, w, lo, hi, **kw):
        c0 = 0.5 * (lo + hi)                                                  # the lattice starts of the block's units
        j = np.flatnonzero(np.any(np.all(np.isclose(c0[:, None, :], starts[None, :, :]), axis=2), axis=1))
        if len(j):
            held.append(len(held) + 1)
            if when is None or held[-1] in when:
                c2, w2 = c.copy(), w.copy()
                w2[j], c2[j] = 0.0, np.array(HOLE)
                return c2, w2, dict(nfev=0, status=0)
        return real(psf, y, sw, c, w, lo, hi, **kw)
    monkeypatch.setattr(IM, 'block_fit', stub)
    return held


def _same(a, b):
    """The two path results are equal bit for bit (the guard's extra key `barred` aside)."""
    assert a['events'] == b['events']
    assert np.array_equal(a['taus'], b['taus']) and np.array_equal(a['snr'], b['snr'])
    assert len(a['hist']) == len(b['hist'])
    for ha, hb in zip(a['hist'], b['hist']):
        assert ha['tau'] == hb['tau'] and ha['n_ent'] == hb['n_ent'] and ha['ids'] == hb['ids'] and ha['nfev'] == hb['nfev'] and ha['chi2'] == hb['chi2']
        assert np.array_equal(ha['c'], hb['c']) and np.array_equal(ha['w'], hb['w'])


# ---- the planted cycle ------------------------------------------------------------------------------------------------------------------------

@REX
def test_the_unguarded_path_loops_at_one_entry_time_until_the_cap(monkeypatch, r_ex):
    psf, y, sw, lat, bad, start = _planted_world()
    _plant(monkeypatch, start)
    R0 = IM.sliding_path(psf, y, sw, lat, **dict(KW, r_ex=r_ex))
    assert len(R0['taus']) == KW['max_entries'] and len(set(R0['taus'])) == 1                    # tau never advances
    assert [e[1:] for e in R0['events']] == [(bad, 'enter'), (bad, 'exit')] * KW['max_entries']  # one lattice point, enter and exit, over and over
    assert len({e[0] for e in R0['events']}) == 1


@REX
def test_the_guard_bars_the_point_after_its_second_exit_and_the_path_goes_on(monkeypatch, r_ex):
    psf, y, sw, lat, bad, start = _planted_world()
    _plant(monkeypatch, start)
    Rg = G.sliding_path_guarded(psf, y, sw, lat, **dict(KW, r_ex=r_ex))
    t0 = Rg['events'][0][0]
    assert Rg['events'][:4] == [(t0, bad, 'enter'), (t0, bad, 'exit'), (t0, bad, 'enter'), (t0, bad, 'exit')]       # two strikes, not one
    assert Rg['barred'] == [(t0, bad)]                                                          # the entry time and the lattice point of the second exit
    after = Rg['events'][4:]
    assert all(e[1] != bad for e in after) and len({e[1] for e in after if e[2] == 'enter'}) > 1  # never again; other kernels enter
    assert len(Rg['taus']) == KW['max_entries'] and list(Rg['taus'][:2]) == [t0, t0] and all(t > t0 for t in Rg['taus'][2:])      # tau advances again after the bar
    assert any(h['w'].size and h['w'].max() > 1000.0 for h in Rg['hist'])                       # the star is fitted by the kernels that were let in


@REX
def test_cycle_guard_false_is_the_unguarded_path_even_in_a_cycle(monkeypatch, r_ex):
    psf, y, sw, lat, bad, start = _planted_world()
    _plant(monkeypatch, start)
    R0 = IM.sliding_path(psf, y, sw, lat, **dict(KW, r_ex=r_ex))
    Ro = G.sliding_path_guarded(psf, y, sw, lat, cycle_guard=False, **dict(KW, r_ex=r_ex))
    _same(R0, Ro)
    assert Ro['barred'] == [] and len(set(Ro['taus'])) == 1


def test_the_switch_is_on_by_default(monkeypatch):
    psf, y, sw, lat, bad, start = _planted_world()
    _plant(monkeypatch, start)
    assert inspect.signature(G.sliding_path_guarded).parameters['cycle_guard'].default is True
    assert G.sliding_path_guarded(psf, y, sw, lat, **KW)['barred'] == G.sliding_path_guarded(psf, y, sw, lat, cycle_guard=True, **KW)['barred'] != []


def test_a_single_exit_followed_by_a_successful_reentry_is_untouched(monkeypatch):
    psf, y, sw, lat, bad, start = _planted_world()
    held = _plant(monkeypatch, start, when={1})                                                 # the unit is forced out once, at its first entry
    R0 = IM.sliding_path(psf, y, sw, lat, **KW)
    n_held = len(held)
    held.clear()
    Rg = G.sliding_path_guarded(psf, y, sw, lat, **KW)
    assert sum(1 for e in R0['events'] if e[2] == 'exit' and e[1] == bad) == 1                  # the plant happened: exactly one exit of the point, then it stays
    assert Rg['barred'] == [] and n_held == len(held)
    _same(R0, Rg)


def test_two_exits_of_one_point_at_different_entry_times_are_not_a_cycle(monkeypatch):
    psf, y, sw, lat, bad, start = _planted_world()
    held = _plant(monkeypatch, start, when={1, 3})                                              # forced out at its entry and again at the third block fit that holds it (a later entry time)
    R0 = IM.sliding_path(psf, y, sw, lat, **KW)
    exits = [e[0] for e in R0['events'] if e[2] == 'exit' and e[1] == bad]
    assert len(exits) == 2 and exits[0] != exits[1] and abs(exits[1] - exits[0]) > 1e-6        # the plant happened, at two entry times
    held.clear()
    Rg = G.sliding_path_guarded(psf, y, sw, lat, **KW)
    assert Rg['barred'] == []
    _same(R0, Rg)


@REX
def test_two_exits_of_one_point_a_few_millionths_apart_are_still_two_entry_times(monkeypatch, r_ex):
    psf, y, sw, lat, bad, start = _planted_world(second=(19.3, 11.8, 4990.0))
    held = _plant(monkeypatch, start, when={1, 3})                                              # forced out at its entry and again when the second star's kernel enters, 3.5e-6 later
    kw = dict(KW, r_ex=r_ex)
    R0 = IM.sliding_path(psf, y, sw, lat, **kw)
    t1, t2 = [e[0] for e in R0['events'] if e[2] == 'exit' and e[1] == bad]
    assert 1e-8 < t2 - t1 < 1e-5                                                                # far above the guard's 1e-12 and far below any slack a careless tolerance would give
    held.clear()
    Rg = G.sliding_path_guarded(psf, y, sw, lat, **kw)
    assert Rg['barred'] == []
    _same(R0, Rg)


def test_two_different_points_that_leave_once_at_the_same_entry_time_are_not_a_cycle(monkeypatch):
    """A star halfway between two lattice kernels (a lattice of 1 px): both enter together (a tie), the stub forces both out at once, both come back and stay. Two exits at one entry time, of two points."""
    psf, y, sw, _ = _world(hole=HOLE + (-3000.0,), stars=((12.5, 12.0, 5000.0),), noise=1e-6)
    lat = IM.Lattice(psf, h=1.0, margin=1.0)
    first = IM.sliding_path(psf, y, sw, lat, **dict(KW, max_entries=2))['events']
    ids = [e[1] for e in first if e[2] == 'enter']
    assert len(ids) == 2 and first[0][0] == first[1][0]                                         # a tie: two kernels at one entry time
    held = _plant(monkeypatch, np.array([lat.pos(*i) for i in ids]), when={1})
    R0 = IM.sliding_path(psf, y, sw, lat, **KW)
    exits = [(e[0], e[1]) for e in R0['events'] if e[2] == 'exit']
    assert len(exits) == 2 and exits[0][0] == exits[1][0] and exits[0][1] != exits[1][1]
    held.clear()
    Rg = G.sliding_path_guarded(psf, y, sw, lat, **KW)
    assert Rg['barred'] == []                                                                   # each point has exited once: the strikes are counted per lattice point
    _same(R0, Rg)


@REX
def test_two_exits_of_one_point_tens_of_picoseconds_apart_are_still_two_entry_times(monkeypatch, r_ex):
    """The mirror image of the first star (the pixel grid and the lattice are symmetric about x = 16) with an amplitude lower by 1e-8: its kernel enters 4.5e-11 after the first's, the second exit of
    the first point is that much later than its first; the guard's slack is 1e-12."""
    psf, y, sw, lat, bad, start = _planted_world(second=(19.7, 11.8, 5000.0 * (1 - 1e-8)), noise=1e-9, unit_weights=True)
    held = _plant(monkeypatch, start, when={1, 3})
    kw = dict(KW, r_ex=r_ex)
    R0 = IM.sliding_path(psf, y, sw, lat, **kw)
    t1, t2 = [e[0] for e in R0['events'] if e[2] == 'exit' and e[1] == bad]
    assert 1e-12 < t2 - t1 < 1e-10 and t1 > 1e-3
    held.clear()
    Rg = G.sliding_path_guarded(psf, y, sw, lat, **kw)
    assert Rg['barred'] == []
    _same(R0, Rg)


def test_each_point_of_a_cycle_is_barred_at_its_own_second_exit(monkeypatch):
    psf, y, sw, lat, bad, start = _planted_world()
    held = _plant(monkeypatch, start)
    nxt = next(e[1] for e in G.sliding_path_guarded(psf, y, sw, lat, **KW)['events'] if e[2] == 'enter' and e[1] != bad)      # the kernel that enters once BAD is barred
    held.clear()
    _plant(monkeypatch, np.array([start, lat.pos(*nxt)]))
    Rg = G.sliding_path_guarded(psf, y, sw, lat, **dict(KW, max_entries=12))
    assert [b for t, b in Rg['barred']] == [bad, nxt] and Rg['barred'][0][0] < Rg['barred'][1][0]
    exits = collections.Counter(e[1] for e in Rg['events'] if e[2] == 'exit')
    assert exits[bad] == 2 and exits[nxt] == 2                                                  # two strikes each, no more
    last = max(t for t, b in Rg['barred'])
    assert all(e[1] not in (bad, nxt) for e in Rg['events'] if e[0] > last) and len(Rg['taus']) == 12


def test_greedy_order_has_no_repeated_entry_time_so_the_guard_never_acts(monkeypatch):
    psf, y, sw, lat, bad, start = _planted_world()
    _plant(monkeypatch, start)
    kw = dict(KW, order='greedy', max_entries=5)
    R0 = IM.sliding_path(psf, y, sw, lat, **kw)
    Rg = G.sliding_path_guarded(psf, y, sw, lat, **kw)
    assert list(R0['taus']) == [1.0, 2.0, 3.0, 4.0, 5.0] and len({e[1] for e in R0['events']}) == 1     # tau is the step count: the same point cycles at distinct 'times'
    assert Rg['barred'] == []                                                                   # the documented limit of the guard
    _same(R0, Rg)


# ---- no double exit: the unguarded path, bit for bit -------------------------------------------------------------------------------------------

OPTIONS = [dict(), dict(fit_radius=4.0), dict(r_ex=1.5), dict(sig_bounds=(1.0, 2.0)), dict(snr_min=3.0, snr_extra=2), dict(tau_max=0.05), dict(local=False, joint=False), dict(order='greedy'),
           dict(sw_fn=lambda m: 1.0 / np.sqrt(0.25 + 0.01 * np.maximum(m, 0.0)), passes=2)]


@pytest.mark.parametrize('opt', OPTIONS, ids=lambda o: ','.join(sorted(o)) or 'defaults')
def test_without_a_double_exit_the_guarded_path_is_the_unguarded_path_bit_for_bit(opt):
    psf, y, sw, lat = _world(seed=5, n=26, stars=((7.4, 8.1, 3000.0), (9.0, 9.2, 1500.0), (18.6, 15.3, 2200.0)), noise=1.0)
    kw = dict(KW, **opt)
    R0 = IM.sliding_path(psf, y, sw, lat, **kw)
    Rg = G.sliding_path_guarded(psf, y, sw, lat, **kw)
    Ro = G.sliding_path_guarded(psf, y, sw, lat, cycle_guard=False, **kw)
    assert Rg['barred'] == [] and Ro['barred'] == []
    _same(R0, Rg)
    _same(R0, Ro)


def test_a_real_path_with_an_exit_and_a_reentry_but_no_double_exit_is_untouched():
    stars = ((7.4, 8.1, 600.0), (8.9, 9.0, 500.0), (14.6, 15.3, 400.0), (14.0, 20.0, 450.0), (20.0, 8.0, 300.0))
    psf, y, sw, lat = _world(seed=1, n=28, stars=stars, noise=2.0)
    kw = dict(max_entries=30, tol=1e-8, passes=1, rtrunc=4.0)
    R0 = IM.sliding_path(psf, y, sw, lat, **kw)
    Rg = G.sliding_path_guarded(psf, y, sw, lat, **kw)
    exits = collections.Counter((t, b) for t, b, k in R0['events'] if k == 'exit')
    assert sum(exits.values()) >= 1 and max(exits.values()) == 1                                # a world in which a unit leaves, but no point leaves twice at one entry time
    assert Rg['barred'] == []
    _same(R0, Rg)


# ---- the source -------------------------------------------------------------------------------------------------------------------------------

GUARD_WORDS = re.compile(r'cycle_guard|barred|last_exit|bars|sliding_path_guarded|IM\.block_fit|IM\._cluster')
ALLOWED_NEW = {'blocked = blocked | lat.near(c[:, :2], r_ex)'}
REMOVED_OK = ('def sliding_path(', 'sw_fn=None, max_nfev=200', 'blocked = taken', 'S = _cluster(', 'cS, wS, info = block_fit(', 'c, w, info = block_fit(', 'return dict(events=events',
              'deviations by which the residual')


def test_the_source_differs_from_imaging_sliding_path_only_in_the_lines_of_the_guard():
    a = [l.strip() for l in inspect.getsource(IM.sliding_path).splitlines()]
    b = [l.strip() for l in inspect.getsource(G.sliding_path_guarded).splitlines()]
    new, gone = [], []
    for tag, i1, i2, j1, j2 in difflib.SequenceMatcher(None, a, b, autojunk=False).get_opcodes():
        if tag != 'equal':
            gone += a[i1:i2]
            new += b[j1:j2]
    bad_new = [l for l in new if not (GUARD_WORDS.search(l) or l in ALLOWED_NEW)]
    assert new and not bad_new, bad_new                                                          # every new line mentions the guard (or is the one line that adds the exclusion zone to `blocked`)
    bad_gone = [l for l in gone if not l.startswith(REMOVED_OK)]
    assert not bad_gone, bad_gone                                                                # and the lines of the original that are gone are the ones those replace
    assert len(gone) <= 9 and sum(1 for l in new if 'cycle_guard' in l) >= 4                     # the def line, the docstring entry, the two uses
    assert 'barred' not in inspect.getsource(IM.sliding_path)                                    # the frozen function is untouched


def test_the_guarded_function_takes_the_arguments_of_the_original_and_one_more():
    po, pg = inspect.signature(IM.sliding_path).parameters, inspect.signature(G.sliding_path_guarded).parameters
    assert list(pg) == list(po) + ['cycle_guard'] and all(pg[k].default == po[k].default for k in po)


def test_the_guard_module_looks_up_the_helpers_of_imaging_at_call_time():
    """`block_fit` and `_cluster` are looked up on the module `imaging` when the path runs, so replacing them there (a test, an experiment) acts on both path functions."""
    src = inspect.getsource(G.sliding_path_guarded)
    assert 'IM.block_fit(' in src and 'IM._cluster(' in src and not re.search(r'(?<![.\w])block_fit\(', src.replace('IM.block_fit(', ''))
