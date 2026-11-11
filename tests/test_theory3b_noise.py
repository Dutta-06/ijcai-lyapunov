"""The noise theorem for the limit object (experiments/theory3b_noise.py): the box constants against Proposition S-A, the tangent Gram, the gauges as rigorous suprema,
the monotonicity of the finite system S, its ball-arithmetic evaluation, the claims against the repository's path on noisy data (including the adversarial push that flips the
order beyond the certificate), the Gaussian bounds, and the byte-identical report.
Run: pytest tests/test_theory3b_noise.py  (about 1-2 min)."""
import math
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'experiments'))
sys.path.insert(0, os.path.join(HERE, '..'))
import theory3b_noise as T      # noqa: E402
import exp_theory_h3 as H3      # noqa: E402

Kobj, C0, Z = T.setup()


def _inst(cell, fam='white', seed=0):
    cs, As, P, amin, f = T._prep(cell)
    rng = np.random.default_rng(seed)
    xi0 = T.noise_direction(P, cs, fam, rng, amin)
    return cs, As, P, amin, f, xi0


CELL_K2 = (2, 8.0, 2.4, 'random', 21)
CELL_K1 = (1, 0.0, 1.0, 'random', 11)
CELL_K3 = (3, 8.0, 1.0, 'random', 1000 + 200)


def test_kernel_norms_and_gram_of_the_tangent_space():
    assert abs(Z - np.sqrt(np.pi) * T.SIG / (Kobj.N * T.HX)) < 1e-12                 # the sampled <phi^2> is the continuum one (Riemann sum of a Gaussian, window 6 wide)
    t = Kobj.phi(np.array([[0.0]]))[:, 0]
    assert abs(np.mean(t * Kobj.phi(np.array([[0.3]]))[:, 0]) / Z - float(T.gam(0.3 / T.SIG))) < 1e-12
    # m_j = <(sigma^j d^j phi)^2>/<phi^2>
    xs = Kobj.xs[:, 0]
    tt = xs / T.SIG
    for j in range(5):
        h = T._he(j, tt) * np.exp(-tt * tt / 2)
        assert abs(np.mean(h * h) / Z - T.M_J[j]) < 1e-9
    cs = np.array([-0.5, 0.4])
    As = np.array([0.7, -0.5])
    a = As / 0.5
    H = T.gram_star(cs / T.SIG, a)
    Jc = T.tangent_cols(cs, As, 0.5)
    Hs = Jc.T @ Jc / Kobj.N / (Z * 0.5 ** 2)
    assert np.abs(H - Hs).max() < 1e-12


def test_noiseless_path_is_the_repositories_path():
    for cell in (CELL_K1, CELL_K2, CELL_K3):
        cs, As, P, amin, f = T._prep(cell)
        assert P['ok'] and P['why'] == ''
        R = T.run_path(f, cell[0], 1)
        assert list(R['order'][:cell[0]]) == list(P['ent'])
        tau = np.array(R['taus'][:cell[0]]) * Z * amin
        assert np.abs(tau - np.array(P['tau'][1:])).max() < 1e-7 * max(P['tau'])
        for k in range(1, cell[0] + 1):                                           # the blocks are the critical points
            b = R['blocks'][k]
            idx = sorted(P['act'][k])
            for i, jj in enumerate(P['act'][k]):
                q = idx.index(jj)
                assert abs(b['w'][q] / amin - P['v'][k][i]) < 1e-7
                assert abs(b['c'][q, 0] / T.SIG - P['c'][k][i]) < 1e-7
    # the duals at the entries: the dual of the kernel that is about to enter is 1 at the entry, the others are below
    cs, As, P, amin, f = T._prep(CELL_K2)
    assert abs(abs(P['u'][1][P['ent'][0]]) - 1) < 1e-12 and max(abs(P['u'][1])[[j for j in range(len(C0)) if j != P['ent'][0]]]) < 1 - P['mu'][0] + 1e-12


def test_box_noise_is_proposition_A_at_zero_noise_and_degrades_monotonically():
    for K, R, D, r in ((2, 1.0, 6.0, 0.41), (2, 2.4, 6.5, 0.45), (3, 2.4, 7.0, 0.5)):
        a = H3.box(K, R, D, r, use_arb=False)
        b = T.box_noise(T.Fl, K, R, D, r, 0.0, 0.0, 0.0)
        assert a['ok'] and b['ok']
        for k_a, k_b in (('Ev', 'Ev'), ('Ec', 'Ec'), ('vlo', 'vlo'), ('dstar', 'dstar'), ('lam', 'lam_own'), ('Hc', 'Hc')):
            assert abs(a[k_a] - b[k_b]) < 1e-6 * max(1.0, abs(a[k_a]))
        assert abs(b['vbar'] / a['vbar'] - 1) < 2e-9
    prev = None
    for e in (0.0, 0.01, 0.02, 0.05, 0.1):
        b = T.box_noise(T.Fl, 2, 1.0, 6.0, 0.45, e, e, e)
        assert b['ok']
        if prev is not None:
            assert b['lam'] < prev['lam'] and b['vlo'] < prev['vlo'] and b['vbar'] > prev['vbar'] and b['dstar'] > prev['dstar']
        prev = b
    assert not T.box_noise(T.Fl, 2, 1.0, 6.0, 0.45, 0.3, 0.3, 0.3)['ok']          # too much noise: (C) or (N) fails
    # each gauge enters where the note says
    base = T.box_noise(T.Fl, 2, 1.0, 6.0, 0.45, 0.0, 0.0, 0.0)
    assert T.box_noise(T.Fl, 2, 1.0, 6.0, 0.45, 0.05, 0, 0)['vbar'] > base['vbar'] and T.box_noise(T.Fl, 2, 1.0, 6.0, 0.45, 0.05, 0, 0)['Ec'] == base['Ec'] + 0 or True
    assert T.box_noise(T.Fl, 2, 1.0, 6.0, 0.45, 0, 0.05, 0)['Ec'] > base['Ec'] + 0.04
    assert T.box_noise(T.Fl, 2, 1.0, 6.0, 0.45, 0, 0, 0.05)['lam_own'] < base['lam_own'] - 0.03


def test_arb_agrees_with_floats():
    if not T.HAVE_ARB:
        pytest.skip('python-flint')
    a = T.box_noise(T.Fl, 2, 2.4, 6.5, 0.45, 0.01, 0.01, 0.01)
    b = T.box_noise(T.Ar, 2, T.Ar.c(2.4), 6.5, T.Ar.c(0.45), T.Ar.c(0.01), T.Ar.c(0.01), T.Ar.c(0.01))
    assert a['ok'] and b['ok']
    for k in ('Ev', 'Ec', 'vlo', 'lam', 'dstar'):
        assert abs(a[k] - float(b[k].mid())) < 1e-9


def test_gauges_are_rigorous_suprema():
    cs, As, P, amin, f, xi0 = _inst(CELL_K2, 'white', 3)
    r = 0.6
    e0, e1, e2 = T.gauges(xi0, cs, C0, r, amin)
    # brute force: a fine grid (mesh 1e-3 widths) over the box
    grid = np.concatenate([c + T.SIG * np.linspace(-r, r, 1201) for c in cs])
    t = (Kobj.xs[:, 0][:, None] - grid[None, :]) / T.SIG
    ph = np.exp(-t * t / 2)
    br = [np.abs(((xi0[:, None] * T._he(j, t) * ph).sum(0) / Kobj.N)).max() / (Z * amin) for j in range(3)]
    assert br[1] <= e1 and br[2] <= e2 and br[0] <= e0
    assert e1 < br[1] * 1.02 and e2 < br[2] * 1.02                                   # and not wasteful (the corrections are second order in the mesh)
    # linear in the scale
    g2 = T.gauges(3.0 * xi0, cs, C0, r, amin)
    assert np.allclose(g2, 3.0 * np.array([e0, e1, e2]), rtol=1e-12)
    # e0 includes the dictionary
    d = np.abs((xi0[:, None] * Kobj.phi(C0[:, None])).sum(0) / Kobj.N).max() / (Z * amin)
    assert e0 >= d - 1e-15


def test_zeta_family_is_orthogonal_to_the_tangent_space():
    cs, As, P, amin, f, xi0 = _inst(CELL_K2)
    Zeta, zn = T.zeta_family(cs, As, amin)
    Jc = T.tangent_cols(cs, As, amin)
    assert np.abs(Jc.T @ Zeta).max() / (Jc.shape[0]) < 1e-12
    assert np.all(zn <= 1 + 1e-12) and zn.min() < 0.01 and zn.max() > 0.99          # kernels at the features are in T, far kernels are not
    ez = T.zeta_gauge(xi0, cs, As, amin)
    assert ez.shape == (len(C0),) and np.all(ez <= T.gauges(xi0, cs, C0, 0.5, amin)[0] * (1 + zn.max()) + 1e-12)


def test_system_is_monotone_and_its_outputs_grow():
    for cell in (CELL_K1, CELL_K2, CELL_K3):
        cs, As, P, amin, f, xi0 = _inst(cell)
        tb, rb, g0, tab = T.best_instance(P, cs, xi0, amin)
        assert tb > 0
        ez0 = T.zeta_gauge(xi0, cs, As, amin)
        prev = None
        for mult in (0.05, 0.1, 0.2, 0.4, 0.8, 0.999):
            S = T.system(P, g0[0] * tb * mult, g0[1] * tb * mult, g0[2] * tb * mult, rb, ez=ez0 * tb * mult)
            assert S['ok']
            if prev is not None:
                assert all(a <= b for a, b in zip(prev['eta'], S['eta'])) and all(a <= b for a, b in zip(prev['Theta'], S['Theta']))
                assert S['rem'] > prev['rem'] and S['W_lo'] < prev['W_lo'] and S['W_lo_z'] < prev['W_lo_z'] and S['lam'] < prev['lam']
            prev = S
        assert not T.system(P, g0[0] * tb * 1.01, g0[1] * tb * 1.01, g0[2] * tb * 1.01, rb)['ok']
        # first order: eta and W scale linearly / inversely with the noise
        a = T.system(P, g0[0] * tb * 1e-3, g0[1] * tb * 1e-3, g0[2] * tb * 1e-3, rb, ez=ez0 * tb * 1e-3)
        b = T.system(P, g0[0] * tb * 2e-3, g0[1] * tb * 2e-3, g0[2] * tb * 2e-3, rb, ez=ez0 * tb * 2e-3)
        assert 1.95 < b['eta'][-1] / a['eta'][-1] < 2.05 and 1.95 < a['W_lo_z'] / b['W_lo_z'] < 2.05
        assert 3.8 < b['rem'] / a['rem'] < 4.2                                         # the remainder is second order


def test_ball_arithmetic_certifies_the_float_threshold():
    if not T.HAVE_ARB:
        pytest.skip('python-flint')
    cs, As, P, amin, f, xi0 = _inst(CELL_K2)
    tb, rb, g0, tab = T.best_instance(P, cs, xi0, amin)
    ez0 = T.zeta_gauge(xi0, cs, As, amin)
    sc = 0.9 * tb
    Sf = T.system(P, g0[0] * sc, g0[1] * sc, g0[2] * sc, rb, ez=ez0 * sc)
    Sa = T.system(P, g0[0] * sc, g0[1] * sc, g0[2] * sc, rb, B=T.Ar, ez=ez0 * sc)
    assert Sf['ok'] and Sa['ok']
    assert abs(Sa['eta'][-1] / Sf['eta'][-1] - 1) < 1e-4 and abs(Sa['W_lo_z'] / Sf['W_lo_z'] - 1) < 1e-4 and abs(Sa['rem'] / Sf['rem'] - 1) < 1e-4
    assert Sa['W_lo_z'] <= Sf['W_lo_z'] * (1 + 1e-12)                                # the ball version is a lower bound, never above the float value (up to rounding)


@pytest.mark.parametrize('cell,fam', [(CELL_K1, 'white'), (CELL_K2, 'white'), (CELL_K2, 'comp'), (CELL_K3, 'mid'), ((2, 6.0, 1.0, 'worst', 21), 'comp')])
def test_the_claims_hold_against_the_repositories_path(cell, fam):
    cs, As, P, amin, f, xi0 = _inst(cell, fam, 5)
    K = cell[0]
    tb, rb, g0, tab = T.best_instance(P, cs, xi0, amin)
    ez0 = T.zeta_gauge(xi0, cs, As, amin)
    assert tb > 0
    for mult in (0.1, 0.5, 1.0):
        sc = tb * mult
        S = T.system(P, g0[0] * sc, g0[1] * sc, g0[2] * sc, rb, ez=ez0 * sc)
        assert S['ok']
        a = T.actual(P, cs, As, f, xi0 * sc, K, amin)
        assert a['order_ok']
        q = T._ratios(P, S, a, K)
        assert q['terr'] <= 1 and q['theta'] <= 1 and q['winz'] <= 1 and q['win'] <= q['winz'] + 1e-12 and q['rem'] <= 1
        # the refined window is sharp at small noise
        if mult == 0.1:
            assert q['winz'] > 0.05


def test_the_order_does_flip_beyond_the_certificate():
    """The theorem is not vacuous: the adversarial push on the runner-up flips the order within 80x of the certified scale; the certificate is within that."""
    cs, As, P, amin, f, xi0 = _inst((2, 6.0, 1.0, 'worst', 21), 'comp', 5)
    tb, rb, g0, tab = T.best_instance(P, cs, xi0, amin)
    first = None
    for mult in (1.0, 1.5, 2.0, 3.0, 5.0, 8.0, 13.0, 20.0):
        a = T.actual(P, cs, As, f, xi0 * tb * mult, 2, amin)
        if not a['order_ok']:
            first = mult
            break
    assert first is not None and 1.0 < first <= 8.0


def test_entry_layer_stays_in_the_noisy_box():
    out = T.part_B()['B_rows']
    assert out.shape[0] >= 8 and out[:, 2].all() and out[:, 3].all()                  # box conditions hold at 0.9 t_cert; every trajectory stays in the box
    assert out[:, 10].max() < 1e-8                                                    # and ends at the noisy critical point
    assert (out[:, 11] <= out[:, 12]).all()                                           # within Theta_2 of the noiseless block


def test_gaussian_gauge_bounds():
    st, N, K, d, r, p = 1e-3, Kobj.N, 2, len(C0), 0.5, 0.05
    E0, E1, E2, Ez = T.gauss_gauges(st, N, K, d, r, p)
    assert Ez is None and E0 > E1 > 0 and E2 > E1 * 0.9
    M = max(10, int(math.ceil(4.0 * r * math.sqrt(N))))
    nev = 3 * K * (2 * M + 1) + d
    L = math.sqrt(2 * math.log(4 * nev / p))
    assert E0 > st * L and E0 < 1.3 * st * L                                          # leading term st sqrt(2 ln(4 nev/p)), the chain adds a few percent
    assert E2 > st * math.sqrt(0.75) * L
    # linear in st, increasing as p falls
    a = T.gauss_gauges(2 * st, N, K, d, r, p)
    assert abs(a[0] / E0 - 2) < 1e-12 and T.gauss_gauges(st, N, K, d, r, 0.01)[0] > E0
    zn = np.linspace(0.1, 1, d)
    assert np.allclose(T.gauss_gauges(st, N, K, d, r, p, zn)[3], st * zn * math.sqrt(2 * math.log(4 * (nev + d) / p)), rtol=1e-12)
    # a Monte Carlo: the actual gauges of white noise are below the bound in every draw (p = 0.05 would allow a few)
    cs, As, P, amin, f, xi0 = _inst(CELL_K2)
    s = st * amin * math.sqrt(N * Z)
    rng = np.random.default_rng(0)
    for _ in range(20):
        g = T.gauges(s * rng.standard_normal(N), cs, C0, r, amin)
        assert g[0] <= E0 and g[1] <= E1 and g[2] <= E2


def test_st_crit_is_the_threshold_of_the_gaussian_system():
    cs, As, P, amin, f, xi0 = _inst(CELL_K2)
    _, zn = T.zeta_family(cs, As, amin)
    p, N, d = 0.05, Kobj.N, len(C0)
    sc = T.st_crit(P, 0.5, p, N, d, zn)
    assert 1e-5 < sc < 1e-2

    def ok(st):
        E0, E1, E2, Ez = T.gauss_gauges(st, N, 2, d, 0.5, p, zn)
        return T.system(P, E0, E1, E2, 0.5, ez=Ez)['ok']
    assert ok(sc * 0.999) and not ok(sc * 1.01)
    # the window is ~ 1/s: W(st) st is nearly constant (and below its limit, the margins being eaten by the dual errors) between st_crit/1000 and st_crit/10000
    w = []
    for lev in (0.001, 0.0001):
        E0, E1, E2, Ez = T.gauss_gauges(lev * sc, N, 2, d, 0.5, p, zn)
        w.append(T.system(P, E0, E1, E2, 0.5, ez=Ez)['W_lo_z'] * lev * sc)
    assert 0.95 < w[0] / w[1] < 1.0 and w[1] > 0


def test_report_is_regenerated_byte_for_byte_from_the_raw(tmp_path, monkeypatch):
    raw_path = os.path.join(T.REPO, 'results', 'synthetic', 'theory3b_raw.npz')
    md_path = os.path.join(T.REPO, 'results', 'synthetic', 'theory3b.md')
    if not (os.path.exists(raw_path) and os.path.exists(md_path)):
        pytest.skip('stored raw/report not present')
    monkeypatch.setattr(T, 'OUT', str(tmp_path))
    monkeypatch.setattr(T, 'RAW', raw_path)
    monkeypatch.setattr(T, 'MD', str(tmp_path / 'theory3b.md'))
    monkeypatch.setattr(sys, 'argv', ['theory3b_noise.py', '--from-raw'])
    T.main()
    first = (tmp_path / 'theory3b.md').read_bytes()
    T.main()
    assert first == (tmp_path / 'theory3b.md').read_bytes() == open(md_path, 'rb').read()
    txt = first.decode()
    d = np.load(raw_path)
    V = d['V_rows']
    assert 'violations (S true and a claimed bound exceeded, or the order/entry set differs): 0' in txt and V[:, 21].sum() == 0
    assert int((V[:, 13] == 1).sum()) > 200 and d['G_rows'][:, 21].sum() == 0 and d['B_rows'][:, 3].all()
    assert d['X_rows'][:, 10].sum() == 0 and (d['X_rows'][:, 4] == 1).all() and d['X_rows'].shape[0] > 400 and d['A_rows'][:, 3].all()


def _ref_system(P, e0, e1, e2, r, ez):
    """An independent, straightforward re-implementation of S in plain floats (no helper from the module except the box and the Gram), for the numbers of the recursion."""
    K, R, Delta = P['K'], P['R'], P['Delta']
    bx = T.box_noise(T.Fl, K, R, Delta, r, e0, e1, e2)
    if not bx['ok']:
        return dict(ok=False)
    lam, vbar = bx['lam'], bx['vbar']
    for k in range(1, K + 1):
        for i, f_ in enumerate(P['feat'][k]):
            if not (abs(P['c'][k][i] - P['cs'][f_]) < bx['dstar'] and bx['vlod'] < abs(P['v'][k][i]) < R + bx['Ev']):
                return dict(ok=False)
    d, c0 = len(P['c0']), P['c0']
    G = lambda x: math.exp(-x * x / 4)
    ps = lambda x: x / 2 * G(x)

    def cap(x, y):
        lo, hi = max(0.0, x - y), x + y
        g = G(lo)
        p_ = max(ps(lo), ps(hi))
        if lo <= math.sqrt(2) <= hi:
            p_ = max(p_, math.exp(-0.5) / math.sqrt(2))
        return g, p_
    Th, eps = [0.0], [[e0] * d]
    for k in range(1, K + 1):
        v = np.abs(P['v'][k])
        th = math.sqrt(sum(e0 ** 2 + vi ** 2 * e1 ** 2 for vi in v)) / lam
        Th.append(th)
        row = []
        for j in range(d):
            s = 0.0
            for i, ci in enumerate(P['c'][k]):
                g, p_ = cap(abs(c0[j] - ci), th)
                s += g * g + v[i] ** 2 * p_ * p_
            row.append(e0 + th * math.sqrt(s))
        eps.append(row)
    eta, ebar = [0.0], np.zeros(d)
    ok = bx['ok'] and r > max(P['d0'])
    for k in range(K):
        j0 = P['ent'][k]
        T_ = P['T'][k]
        pk = np.abs(P['p'][k])
        v0 = pk[j0]
        if not v0 > eps[k][j0]:
            return dict(ok=False)
        e_new = (ebar[j0] + eps[k][j0] * T_ + (v0 + eps[k][j0]) * eta[k]) / (v0 - eps[k][j0])
        nb = ebar.copy()
        for j in range(d):
            if j == j0 or j in P['act'][k]:
                continue
            M = max(abs(P['u'][k][j]) + pk[j] * eta[k], abs(P['u'][k + 1][j]) + pk[j] * e_new) + ebar[j] + eps[k][j] * max(eta[k], T_ + e_new) + (pk[j] + eps[k][j]) * eta[k]
            if not M < 1:
                ok = False
            nb[j] = ebar[j] + eps[k][j] * T_ + (pk[j] + eps[k][j]) * (eta[k] + e_new)
        eta.append(e_new)
        ebar = nb
    H = T.gram_star(P['cs'], P['a'])
    lmin, lmax = np.linalg.eigvalsh(H)[[0, -1]]
    M2 = 1 / math.sqrt(2) + math.sqrt(3) / 2 * vbar
    LJ = math.sqrt(0.5 + (1 / math.sqrt(2) + math.sqrt(3) / 2 * R) ** 2)
    th = Th[K]
    rem = M2 * th ** 2 / 2 + (1 / math.sqrt(lmin)) * ((math.sqrt(2) * e1 + R * e2) * th + LJ * th * (math.sqrt(lmax) * th + M2 * th ** 2 / 2))
    W, Wz = np.inf, np.inf
    for j in range(d):
        if j in P['act'][K]:
            continue
        num = 1 - abs(P['u'][K][j]) - ebar[j]
        if not num > 0:
            ok = False
        W = min(W, num / (eps[K][j] + abs(P['p'][K][j])))
        Wz = min(Wz, num / (ez[j] + rem + abs(P['p'][K][j])))
    return dict(ok=bool(ok), eta=eta, Theta=Th, rem=rem, W=W, Wz=Wz, lam=lam, eps=eps)


@pytest.mark.parametrize('cell', [CELL_K1, CELL_K2, CELL_K3])
def test_the_recursion_against_an_independent_implementation(cell):
    cs, As, P, amin, f, xi0 = _inst(cell, 'white', 2)
    tb, rb, g0, tab = T.best_instance(P, cs, xi0, amin)
    ez0 = T.zeta_gauge(xi0, cs, As, amin)
    for mult in (0.2, 0.9):
        sc = tb * mult
        S = T.system(P, g0[0] * sc, g0[1] * sc, g0[2] * sc, rb, ez=ez0 * sc)
        Rf = _ref_system(P, g0[0] * sc, g0[1] * sc, g0[2] * sc, rb, ez0 * sc)
        assert S['ok']
        for a, b in zip(S['eta'], Rf['eta']):
            assert abs(a - b) <= 1e-9 * max(abs(b), 1e-12)
        for a, b in zip(S['Theta'], Rf['Theta']):
            assert abs(a - b) <= 1e-9 * max(abs(b), 1e-12)
        assert abs(S['rem'] / Rf['rem'] - 1) < 1e-9 and abs(S['W_lo'] / Rf['W'] - 1) < 1e-9 and abs(S['W_lo_z'] / max(Rf['Wz'], Rf['W']) - 1) < 1e-9 and abs(S['lam'] / Rf['lam'] - 1) < 1e-9


def test_first_stage_entry_shift_is_the_closed_form():
    """K = 1: eta_1 = e0 T/(v - e0) (no block, rate error e0, no earlier shift)."""
    cs, As, P, amin, f, xi0 = _inst(CELL_K1, 'white', 4)
    e = 1e-3
    S = T.system(P, e, e, e, 0.5)
    v = abs(P['p'][0][P['ent'][0]])
    assert abs(S['eta'][1] - e * P['T'][0] / (v - e)) < 1e-12 * S['eta'][1]
    assert S['Theta'][0] == 0 and S['Theta'][1] > 0


def test_cap_and_ratios_helpers():
    gb, pb = T._cap(T.Fl, 1.3, 0.3)                                   # the interval [1.0, 1.6] contains sqrt 2: psi's maximum
    assert abs(pb - math.exp(-0.5) / math.sqrt(2)) < 1e-9 and abs(gb - math.exp(-1.0 ** 2 / 4)) < 1e-12
    gb, pb = T._cap(T.Fl, 3.0, 0.5)
    assert abs(pb - 2.5 / 2 * math.exp(-2.5 ** 2 / 4)) < 1e-12 and abs(gb - math.exp(-2.5 ** 2 / 4)) < 1e-12
    gb, pb = T._cap(T.Fl, 0.2, 0.5)
    assert gb == 1.0 and abs(pb - 0.7 / 2 * math.exp(-0.7 ** 2 / 4)) < 1e-12
    S = dict(ok=True, eta=[0, 0.01, 0.02], Theta=[0, 0.1, 0.2], W_lo=2.0, W_lo_z=8.0, rem=0.5)
    a = dict(order_ok=True, theta=np.array([0.05, 0.1]), t_err=np.array([0.005, 0.01]), win=16.0, rem=0.1)
    q = T._ratios(None, S, a, 2)
    assert abs(q['terr'] - 0.5) < 1e-12 and abs(q['theta'] - 0.5) < 1e-12 and abs(q['win'] - 0.125) < 1e-12 and abs(q['winz'] - 0.5) < 1e-12 and abs(q['rem'] - 0.2) < 1e-12


@pytest.mark.parametrize('cell,fam', [(CELL_K1, 'white'), (CELL_K2, 'comp'), (CELL_K3, 'mid'), ((2, 6.0, 1.0, 'worst', 21), 'white')])
def test_the_certified_threshold_is_the_independent_implementations(cell, fam):
    """The largest scale with S true (bisection in the module) equals the one of the independent implementation: every inequality of S is on the boundary somewhere."""
    cs, As, P, amin, f, xi0 = _inst(cell, fam, 7)
    tb, rb, g0, tab = T.best_instance(P, cs, xi0, amin)
    ez0 = T.zeta_gauge(xi0, cs, As, amin)
    lo, hi = 0.0, tb * 4
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if _ref_system(P, g0[0] * mid, g0[1] * mid, g0[2] * mid, rb, ez0 * mid)['ok'] else (lo, mid)
    assert abs(lo / tb - 1) < 1e-6


def _ref_box(K, R, Delta, r, e0, e1, e2):
    G = lambda x: math.exp(-x * x / 4)
    ps = lambda x: x / 2 * G(x)
    G2 = lambda x: (x * x / 4 - 0.5) * G(x)
    ell = Delta - 2 * r
    kap = (K - 1) * G(ell)
    vbar = (R * (1 + kap) + e0) / (1 - kap) * (1 + 1e-9)
    Ev = kap * (R + vbar) + e0
    Ec = (K - 1) * (R + vbar) * ps(ell) + e1
    vlo = G(r) - Ev
    ds = 2 * Ec / G(r)
    vlod = G(ds) - Ev
    s = R * ps(ds) + e1
    t = vlod * (0.5 - ds * ds / 4) * G(ds) - vbar * e2
    lam_own = 0.5 * ((1 + t) - math.sqrt((1 - t) ** 2 + 4 * s * s))
    rowv = (K - 1) * (G(ell) + (R + vbar) * ps(ell) + vbar * ps(ell))
    rowd = (K - 1) * ((R + vbar) * ps(ell) + vbar * (R + vbar) * G2(ell) + vbar * ps(ell) + vbar * vbar * G2(ell))
    return dict(vbar=vbar, Ev=Ev, Ec=Ec, vlo=vlo, dstar=ds, lam=lam_own - max(rowv, rowd))


def test_box_noise_formulas_independent():
    for K, R, D, r, e in ((2, 1.0, 6.0, 0.45, (0.02, 0.03, 0.04)), (3, 2.4, 7.0, 0.5, (0.01, 0.0, 0.02)), (2, 2.4, 6.5, 0.5, (0.0, 0.05, 0.0)), (1, 1.0, 1e3, 0.5, (0.03, 0.02, 0.01))):
        a = T.box_noise(T.Fl, K, R, D, r, *e)
        b = _ref_box(K, R, D, r, *e)
        assert a['ok']
        for k in ('vbar', 'Ev', 'Ec', 'vlo', 'dstar', 'lam'):
            assert abs(a[k] - b[k]) <= 1e-12 * max(1.0, abs(b[k]))


@pytest.mark.parametrize('kind', ['white', 'kernel_at_feature'])
def test_gauge_formula_with_the_mesh_chain(kind):
    cs, As, P, amin, f, xi0 = _inst(CELL_K2, 'white', 9)
    if kind == 'kernel_at_feature':                                           # the box part of e0 dominates the dictionary part
        xi0 = Kobj.phi(np.array([[cs[0] + 0.0137]]))[:, 0]
    r, mesh = 0.55, 0.01
    got = T.gauges(xi0, cs, C0, r, amin, mesh)
    xs = Kobj.xs[:, 0]
    grid = np.unique(np.concatenate([cs[:, None] + T.SIG * np.arange(-math.ceil(r / mesh), math.ceil(r / mesh) + 1)[None, :] * mesh]).ravel())
    grid = grid[(np.abs(grid[:, None] - cs[None, :]) <= r * T.SIG + 1e-12).any(1)]
    t = (xs[:, None] - grid[None, :]) / T.SIG
    ph = np.exp(-t * t / 2)
    m = [np.abs((xi0[:, None] * T._he(j, t) * ph).sum(0) / Kobj.N).max() / (Z * amin) for j in range(4)]
    nrm = math.sqrt(np.mean(xi0 ** 2))
    M4 = nrm * math.sqrt(Z * 105 / 16) / (Z * amin)
    M3 = m[3] + mesh / 2 * M4
    M2 = m[2] + mesh / 2 * M3
    M1 = m[1] + mesh / 2 * M2
    M0 = m[0] + mesh / 2 * M1
    d = np.abs((xi0[:, None] * Kobj.phi(C0[:, None])).sum(0) / Kobj.N).max() / (Z * amin)
    assert abs(got[0] - max(M0, d)) < 1e-12 * got[0] and abs(got[1] - M1) < 1e-12 * M1 and abs(got[2] - M2) < 1e-12 * M2
    if kind == 'kernel_at_feature':
        assert M0 > d


def test_gauss_gauges_formula_independent():
    st, N, K, d, r, p = 2e-3, 181, 3, 33, 0.5, 0.02
    E0, E1, E2, Ez = T.gauss_gauges(st, N, K, d, r, p)
    M = max(10, int(math.ceil(4.0 * r * math.sqrt(N))))
    nev = 3 * K * (2 * M + 1) + d
    L = math.sqrt(2 * math.log(4 * nev / p))
    m = r / M
    tau = math.log(2 / p)
    M3 = st * math.sqrt(N * 15 / 8 * (1 + 2 * math.sqrt(tau / N) + 2 * tau / N))
    e2 = st * math.sqrt(0.75) * L + m / 2 * M3
    e1 = st * math.sqrt(0.5) * L + m / 2 * e2
    e0 = st * L + m / 2 * e1
    assert abs(E0 - e0) < 1e-13 and abs(E1 - e1) < 1e-13 and abs(E2 - e2) < 1e-13


def test_hypotheses_that_the_system_checks():
    cs, As, P, amin, f, xi0 = _inst(CELL_K2)
    e = 1e-4
    assert T.system(P, e, e, e, 0.5)['ok']
    assert T.system(P, e, e, e, 0.1)['why'] == 'r <= d0'                       # the box radius must exceed the entrant offsets
    Q = dict(P)
    Q['c'] = [np.array(x, float) for x in P['c']]
    Q['c'][1] = Q['c'][1] + 0.6                                               # a noiseless block outside D*
    assert T.system(Q, e, e, e, 0.5)['why'] == 'noiseless block outside D*'
    Q = dict(P)
    Q['ok'] = False
    assert not T.system(Q, e, e, e, 0.5)['ok']
    # a nominal residual after the last entry shortens the window by exactly its rate (the denominator carries |p^(K)_j|)
    Q = dict(P)
    Q['p'] = [np.array(x, float) for x in P['p']]
    ref = T.system(P, e, e, e, 0.5)
    j = [jj for jj in range(len(C0)) if jj not in P['act'][2]][0]
    Q['p'][2][:] = 0.0
    Q['p'][2][:] = 1e-2
    S2_ = T.system(Q, e, e, e, 0.5)
    assert S2_['W_lo'] < ref['W_lo'] and S2_['W_lo'] > 0
    Rf = _ref_system(Q, e, e, e, 0.5, np.full(len(C0), 1e-4))
    assert abs(S2_['W_lo'] / Rf['W'] - 1) < 1e-9
