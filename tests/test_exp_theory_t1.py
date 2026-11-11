"""Theorem 1 with an explicit threshold: the unit's closed forms (Lemma 2.1-2.4), the box constants and Hessian row sums, the limit recursion, the certificate against
the integrated network, the ingredient inequalities of the proof, and the deterministic report (`--from-raw` byte for byte)."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'experiments'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import exp_theory_h3 as H3      # noqa: E402
import exp_theory_t1 as T       # noqa: E402


def test_unit_closed_forms_match_the_repository_unit_and_lemma_2_1():
    from src.core import factorized_rbf as FR
    from src.core.mirror_free import FreeLaw
    eps = 1e-8
    U = T.Unit(eps)
    law = FreeLaw(**T.LAW)
    ro = FR.PaperLaw(law)
    for x in (0.01, 0.4, 1.3, 2.7):
        beta = eps + 2 * U.vs * np.sinh(x) ** 2                       # the conserved curve beta = eps + 2 v_s sinh^2(alpha/a_max)
        w_repo = float(ro.w(np.array([[x * T.LAW['a_max'], np.sqrt(beta)]]))[0])
        assert abs(float(U.w_of_x(x)) - w_repo) < 1e-14 * max(1.0, abs(w_repo))
    # d eta/dx = a_max/(G q) and M = dw/d eta by differences of the closed forms
    x = np.array([0.2, 1.0, 2.0])
    h = 1e-6
    eta = lambda xx: U.lam - U.s_of_x(xx)
    deta = (eta(x + h) - eta(x - h)) / (2 * h)
    assert np.allclose(deta, U.deta_dx(x), rtol=1e-6)
    dw = (U.w_of_x(x + h) - U.w_of_x(x - h)) / (2 * h)
    assert np.allclose(dw / deta, U.M_of_x(x), rtol=1e-6)
    xs = np.linspace(0.0, 4.0, 400)
    assert np.all(np.diff(U.s_of_x(xs[1:])) < 0) and np.all(np.diff(U.M_of_x(xs)) > 0) and np.all(np.diff(U.w_of_x(xs)) > 0)
    # the wall: lambda sqrt(eps) -> c_u
    cu = np.pi * T.LAW['a_max'] / (2 * np.sqrt(2.0) * T.LAW['G'] * np.sqrt(U.vs))
    assert abs(U.lam * np.sqrt(eps) / cu - 1) < 1e-4


def test_mass_identity_and_the_dormant_core_bound():
    for eps in (1e-6, 1e-10):
        U = T.Unit(eps)
        x = np.linspace(0, 2.5, 100001)
        eta = U.lam - U.s_of_x(np.maximum(x, 1e-12))
        eta[0] = 0.0
        q = np.trapezoid(U.w_of_x(x) * np.gradient(eta, x), x)
        assert abs(q - U.Psi_x(2.5)) < 1e-6 * U.Psi_x(2.5)                  # Lemma 2.4
        for mu in (0.05, 0.1, 0.3, 0.5):
            ex, cl = U.core_bound(mu)
            b = float(U.B(mu * U.lam))
            assert b <= float(ex) * (1 + 1e-6) and float(ex) <= float(cl) * (1 + 1e-9)       # Lemma 2.2: exact <= first bound <= closed form
    U = T.Unit(1e-12)
    for s in (1e3, 1e4):
        assert abs(float(U.B(s)) * s ** 3 / (4 / (T.LAW['G'] * T.LAW['gamma'] ** 2 * T.LAW['m'])) - 1) < 2e-3        # Cor. 2.3: B(s) s^3 -> 4/(G gamma^2 m)


def test_box_constants_agree_with_the_note_and_the_regime_is_slightly_smaller():
    for (K, R, D, r, Om) in [(2, 1.0, 6.0, 0.5, 0.0), (2, 2.4, 8.0, 0.6, 0.05), (3, 1.0, 8.0, 0.5, 0.02)]:
        a = T.box_t1(K, R, D, r, Om)
        b = H3.box(K, R, D, r, use_arb=False, Omega=Om, need_N=True)
        assert a['ok'] and b['ok']
        for ka, kb in (('vbar', 'vbar'), ('Ev', 'Ev'), ('Ec', 'Ec'), ('vlo', 'vlo'), ('Cmargin', 'Cmargin'), ('Hc', 'Hc')):
            assert abs(a[ka] - b[kb]) < 1e-9 * max(1.0, abs(b[kb]))
    thr = T.delta_t1(2, 1.0, 0.4167, 0.0)
    ent = H3.entry_threshold(2, 1.0, 0.4167)
    assert ent - 0.02 <= thr <= ent + 0.6                                  # Delta_T1 is Prop. A's threshold plus a small margin for the Gershgorin conditions
    assert T.delta_t1(2, 2.4, 0.4167, 0.0) > thr and T.delta_t1(3, 1.0, 0.4167, 0.0) > thr


def test_cross_hessian_row_sums_bound_the_exact_hessian():
    rng = np.random.default_rng(0)
    G, dG = (lambda d: np.exp(-d ** 2 / 4)), (lambda d: -d / 2 * np.exp(-d ** 2 / 4))

    def grad(y, cs, As):
        k = len(y) // 2
        v, d = y[:k], y[k:]
        c = cs[:k] + d
        w = v * np.sign(As[:k])
        g = G(c[:, None] - cs[None, :]) @ As - G(c[:, None] - c[None, :]) @ w
        M = dG(c[:, None] - c[None, :])
        np.fill_diagonal(M, 0.0)
        h = dG(c[:, None] - cs[None, :]) @ As - M @ w
        return np.r_[-np.sign(As[:k]) * g, -w * h]
    for (K, R, D, r) in [(2, 1.4, 6.2, 0.6), (3, 1.5, 6.0, 0.5)]:
        b = T.box_t1(K, R, D, r, 0.0)
        cs = np.arange(K) * D
        a = np.r_[1.0, rng.uniform(1, R, K - 1)]
        As = a * np.where(rng.random(K) < 0.5, -1, 1)
        wv = wd = 0.0
        for _ in range(150):
            y = np.r_[rng.uniform(b['vlo'], b['vbar'], K), rng.uniform(-r, r, K)]
            n = len(y)
            H = np.zeros((n, n))
            for j in range(n):
                e = np.zeros(n)
                e[j] = 1e-6
                H[:, j] = (grad(y + e, cs, As) - grad(y - e, cs, As)) / 2e-6
            H = 0.5 * (H + H.T)
            Hown = np.zeros_like(H)
            for i in range(K):
                ai, di, vi = abs(As[i]), y[K + i], y[i]
                Hown[i, i] = 1
                Hown[i, K + i] = Hown[K + i, i] = ai * di / 2 * G(di)
                Hown[K + i, K + i] = vi * ai * (0.5 - di ** 2 / 4) * G(di)
            C = H - Hown
            wv, wd = max(wv, np.abs(C[:K]).sum(1).max()), max(wd, np.abs(C[K:]).sum(1).max())
        assert wv <= b['rowv'] and wd <= b['rowd']


def test_limit_recursion_first_entry_and_fitted_final_block():
    cs, As, c0 = T.make_target(2, 8.0, 1.0, 13)
    prep = T.prepare(cs, As, c0)
    L = prep['L']
    p0 = T.Gam(c0[:, None] - cs[None, :]) @ As
    j1 = int(np.argmax(np.abs(p0)))
    assert L['events'][0]['kernel'] == j1 and abs(L['events'][0]['tau'] - 1 / abs(p0[j1])) < 1e-12 * L['events'][0]['tau']
    assert np.allclose(np.sort(np.abs(L['events'][-1]['w'])), np.sort(np.abs(As)), atol=1e-9)         # the final block is the target
    assert np.abs(L['P'][-1]).max() < 1e-9                                                              # and the duals stop
    assert prep['sign_ok'] and 0 < prep['mu_min'] < 1 and prep['gap_min'] > 0


def test_certificate_k1_follows_the_network_and_its_ingredients_hold():
    cs, As, c0 = T.make_target(1, 0.0, 1.0, 2)
    prep = T.prepare(cs, As, c0)
    o3, pi3 = T.best_cert(prep, 1e-3)
    assert not o3['ok'] and pi3 < 0                                        # far above the threshold
    eps = 1e-8
    o, pi = T.best_cert(prep, eps)
    assert o['ok'] and pi >= 0
    sol, U = T.simulate(prep, eps, nlocal=300)
    ok, exact = T.follows(prep, eps, sol, U)
    assert ok and exact
    row = T.diagnostics(prep, eps, T.PARS[pi])
    c = {k: row[i] for i, k in enumerate(T.D_COLS[3:])}
    assert c['E_over_Dv'] <= 1 and c['shift_ok'] == 1 and c['T1_ok'] == 1 and c['rel_ok'] == 1 and c['Ee_ok'] == 1
    assert c['mono'] == 1 and c['gslack_min'] > 0 and c['d_ok'] == 1 and c['w_ok'] == 1 and c['omega_act'] <= c['omega_bound'] and c['formed_ok'] == 1


def test_certificate_k2_has_no_violation_on_a_coarse_grid():
    cs, As, c0 = T.make_target(2, 8.0, 1.0, 13)
    prep = T.prepare(cs, As, c0)
    okc, okm, okx, Dd, Pi = T.scan_target(prep, np.array([-4.0, -6.0, -8.0, -10.0]))
    assert not np.any((okc > 0) & (okm == 0))                              # cert true => the network forms exactly the limit's units
    assert okc[-1] == 1 and okm[-1] == 1 and okc[0] == 0


def test_report_is_regenerated_byte_for_byte_from_the_raw(tmp_path):
    out = str(tmp_path)
    old = (T.OUT, T.RAW, T.MD)
    try:
        T.OUT, T.RAW, T.MD = out, os.path.join(out, 'theory_t1_raw.npz'), os.path.join(out, 'theory_t1.md')
        T.main(['--run', '--quick', '--parts', 'A,C,D,F,G,H'])
        a = open(T.MD, 'rb').read()
        raw1 = open(T.RAW, 'rb').read()
        T.main(['--from-raw'])
        b = open(T.MD, 'rb').read()
        assert a == b and len(a) > 1000
        assert open(T.RAW, 'rb').read() == raw1
        assert b'violations' in a and b'**0**' in a
    finally:
        T.OUT, T.RAW, T.MD = old
