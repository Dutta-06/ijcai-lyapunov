"""Tests of experiments/exp_mirror_branch.py, Part G (round 144): the mixed start of Part C resolved on both scales."""
import os
import sys

import numpy as np
import pytest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'experiments'))
from src.core import mirror_free as MF
import exp_mirror_branch as X

RAW = os.path.join(ROOT, 'results', 'chaos', 'mirror_branch_raw.npz')


@pytest.fixture(scope='module')
def mixed():
    m = X.mixed_C()
    fr = np.flatnonzero(~m['odd'])
    off = np.setdiff1d(np.arange(X.DC), m['S'])
    ob = np.array([j for j in off if m['odd'][j]])
    return m, fr, off, ob


def test_mixed_design_is_part_C_and_its_vertex_is_degenerate(mixed):
    m, fr, off, ob = mixed
    assert m['S'].tolist() == [1, 5, 9] and m['S'].size < X.NC          # fewer than N atoms: the truth's support
    assert np.allclose(m['Phi'] @ m['wsc'], m['y'], atol=1e-10)
    assert np.all(m['sg'][m['S']] * m['wsc'][m['S']] > 0)


def test_qp_satisfies_its_kkt_conditions(mixed):
    m, fr, off, ob = mixed
    Phi, S, sg, rho = m['Phi'], m['S'], m['sg'], m['rho']
    zs, act, viol = X.mixed_qp(Phi, S, sg, rho, ob, fr)
    assert viol < 1e-8
    assert np.allclose(Phi[:, S].T @ zs, rho[S] * sg[S], atol=1e-9)
    h = sg[ob] * (Phi[:, ob].T @ zs) + rho[ob]
    assert np.all(h <= 1e-8) and np.allclose(h[np.isin(ob, act)], 0, atol=1e-8)
    # stationarity: Φ_F Φ_Fᵀζ* lies in the span of Φ_S and the active atoms' columns, with nonnegative multipliers on the latter
    A = np.hstack([Phi[:, S]] + [sg[j] * Phi[:, [j]] for j in act])
    lam = np.linalg.lstsq(A, -Phi[:, fr] @ (Phi[:, fr].T @ zs), rcond=None)[0]
    assert np.linalg.norm(A @ lam + Phi[:, fr] @ (Phi[:, fr].T @ zs)) < 1e-8 * np.linalg.norm(Phi[:, fr] @ (Phi[:, fr].T @ zs))
    assert np.all(lam[S.size:] >= -1e-10)


def test_resolved_projection_agrees_with_projection_where_it_resolves(mixed):
    """At ε = 10⁻⁴ the off-support terms lie far above `projection`'s stopping residual: both solvers give the same point."""
    m, fr, off, ob = mixed
    L, Phi, y = m['L'], m['Phi'], m['y']
    eps = 1e-4
    P = MF.ParticleSet(L, np.where(m['odd'], m['sg'] * m['amp0'], 0.0), np.full(X.DC, eps))
    e0 = P.eta(P.x0)
    wb, it, res = MF.projection(P, Phi, y, e0, exact=True)
    st = P.state_of_w(wb, np.flatnonzero(wb != 0))
    z0 = np.linalg.lstsq(Phi.T, P.eta(st) - e0, rcond=None)[0]
    r = X.mixed_resolved(m, eps, z0 * 1.05, True)                       # from a perturbed start
    assert r['res'] < 1e-12 and r['perp'] < 1e-12
    assert abs(r["free"] / np.abs(wb[fr]).sum() - 1) < 1e-6          # projection: residual 1.7e-12


def test_resolved_projection_is_the_bregman_projection_below_the_old_floor(mixed):
    """At ε = 10⁻¹² the resolved point is the Bregman projection (η − η(0) ∈ range Φᵀ by construction; the fit holds on both scales),
    and the off-support contribution to Φw lies below `projection`'s stopping residual."""
    m, fr, off, ob = mixed
    Phi, y = m['Phi'], m['y']
    zs = X.mixed_qp(Phi, m['S'], m['sg'], m['rho'], ob, fr)[0]
    r = X.mixed_resolved(m, 1e-12, zs * np.log(1e12), True)
    assert r['res'] < 1e-12 and r['perp'] < 1e-12
    assert r['offout'] < 1e-12 * (1 + np.linalg.norm(y))
    rl = X.mixed_resolved(m, 1e-12, r['z'], False)                         # the linear core agrees to its O(ςGεη²) error
    assert abs(rl['free'] / r['free'] - 1) < 1e-4


def test_raw_part_G_trends_and_consistency_with_part_C():
    if not os.path.exists(RAW):
        pytest.skip('no raw file')
    z = np.load(RAW, allow_pickle=True)
    raw = {k: z[k][0] for k in z.files}
    assert set('ABCDEFG') <= set(raw)
    G, C = raw['G'], raw['C']
    # C's projection rows are C's mixed-start rows
    for c, mx in zip(G['cproj'], C['mx']):
        assert c['eps'] == mx['eps'] and c['free'] == mx['free'] and c['res'] == mx['res']
    rows = G['rows']
    assert rows[0]['eps'] == 1e-4 and abs(rows[0]["free"] / C["mx"][0]["free"] - 1) < 1e-6
    assert all(r['res'] < 1e-10 and r['perp'] < 1e-12 for r in rows)
    r2L = [r['free_e2L'] for r in rows]
    assert all(b > a for a, b in zip(r2L, r2L[1:])) and r2L[-1] < G['const']    # rising toward the leading order from below
    assert abs(r2L[-1] / G['const'] - 1) < 0.1
    act = G['act']
    ex = [min(r['expo'][j] for j in act) for r in rows]
    assert all(b > a for a, b in zip(ex[1:], ex[2:])) and 1.9 < ex[-1] < 2.0
    dev = [r['zdev'] for r in rows]
    assert all(b < a for a, b in zip(dev, dev[1:]))
    assert G['gordan'] and G['core_err'] < 1e-4
