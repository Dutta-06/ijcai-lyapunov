"""Tests of src/core/imaging.py: the pixel-integrated PSF (normalisation, derivatives, separable synthesis/correlation), the start lattice,
the block fit (non-negative, bounded centres), the sliding path on synthetic frames (one, two, three emitters; noise-free exactness;
entry order; exits) and the two stopping rules. Fast: 24 × 24 frames."""
import numpy as np
import pytest

from src.core import imaging as IM

SIG = 1.2


def psf(n=24, sig=SIG):
    return IM.PixelPSF(n, n, sig)


# ---- the PSF ----------------------------------------------------------------------------------------------------------------------

def test_psf_sums_to_one_and_peaks_at_the_centre():
    P = IM.PixelPSF(41, 41, 1.3)
    phi = P.phi(np.array([[20.3, 19.6]]))[:, 0].reshape(41, 41)
    assert abs(phi.sum() - 1.0) < 1e-12
    i, j = np.unravel_index(phi.argmax(), phi.shape)
    assert (i, j) == (20, 20)                         # row 20 (y = 19.6 → 20), column 20 (x = 20.3 → 20)
    assert phi[i, j] > phi[i, j - 1] and phi[i, j] > phi[i + 1, j]


def test_psf_matches_direct_pixel_integral():
    """Ψ is the integral of the Gaussian density over the pixel: compare with a fine midpoint quadrature."""
    P = IM.PixelPSF(9, 9, 0.9)
    c = np.array([[4.2, 3.7]])
    phi = P.phi(c)[:, 0].reshape(9, 9)
    m = 400
    t = (np.arange(m) + 0.5) / m - 0.5
    g = lambda u, c0: np.exp(-(u - c0) ** 2 / (2 * 0.9 ** 2)) / (np.sqrt(2 * np.pi) * 0.9)
    ref = np.zeros((9, 9))
    for i in range(9):
        for j in range(9):
            ref[i, j] = g(j + t, 4.2).mean() * g(i + t, 3.7).mean()
    assert np.abs(phi - ref).max() < 1e-7


def test_psf_derivative_by_finite_differences():
    P = psf()
    c = np.array([[10.3, 11.9], [5.1, 7.4]])
    h = 1e-6
    ph = P.phi(c)
    gx = (P.phi(c + [h, 0]) - P.phi(c - [h, 0])) / (2 * h)
    gy = (P.phi(c + [0, h]) - P.phi(c - [0, h])) / (2 * h)
    gy_, gx_, dgy, dgx = P.factors(c)
    dx = (gy_[:, None, :] * dgx[None, :, :]).reshape(P.N, -1)
    dy = (dgy[:, None, :] * gx_[None, :, :]).reshape(P.N, -1)
    assert np.abs(dx - gx).max() < 1e-8 and np.abs(dy - gy).max() < 1e-8
    assert np.abs(ph - (gy_[:, None, :] * gx_[None, :, :]).reshape(P.N, -1)).max() == 0


def test_synth_and_corr_equal_the_dense_products():
    rng = np.random.default_rng(1)
    P = psf()
    c = rng.uniform(2, 21, (5, 2))
    w = rng.uniform(1, 5, 5)
    Rw = rng.standard_normal((24, 24))
    Phi = P.phi(c)
    assert np.abs(P.synth(c, w).ravel() - Phi @ w).max() < 1e-12
    v, vx, vy = P.corr(c, Rw)
    assert np.abs(v - Phi.T @ Rw.ravel()).max() < 1e-11
    h = 1e-6
    for d, g in ((0, vx), (1, vy)):
        e = np.zeros(2); e[d] = h
        num = (P.corr(c + e, Rw)[0] - P.corr(c - e, Rw)[0]) / (2 * h)
        assert np.abs(num - g).max() < 1e-7
    assert P.synth(np.zeros((0, 2)), np.zeros(0)).shape == (24, 24) and not P.synth(np.zeros((0, 2)), np.zeros(0)).any()


def test_jacobian_by_finite_differences():
    rng = np.random.default_rng(2)
    P = psf()
    c = rng.uniform(3, 20, (3, 2))
    w = rng.uniform(1, 4, 3)
    sw = rng.uniform(0.5, 1.5, (24, 24))
    y = rng.standard_normal((24, 24))
    z = np.r_[w, c[:, 0], c[:, 1]]
    k = 3

    def res(z):
        return (sw * (y - P.synth(np.c_[z[k:2 * k], z[2 * k:]], z[:k]))).ravel()
    J = P.jac(c, w, sw)
    Jn = np.empty_like(J)
    for i in range(3 * k):
        e = np.zeros(3 * k); e[i] = 1e-6
        Jn[:, i] = (res(z + e) - res(z - e)) / 2e-6
    assert np.abs(J - Jn).max() < 1e-6


def test_lattice_correlation_equals_per_kernel_correlation():
    rng = np.random.default_rng(3)
    P = psf()
    L = IM.Lattice(P, h=0.5, margin=1.0)
    Rw = rng.standard_normal((24, 24))
    V = L.corr(Rw)
    assert V.shape == (L.ys.size, L.xs.size) == (51, 51)      # from −1 to 24 in steps of ½
    for b, a in ((0, 0), (7, 31), (50, 50), (20, 3)):
        assert abs(V[b, a] - P.corr(L.pos(b, a)[None, :], Rw)[0][0]) < 1e-11
    assert np.allclose(L.pos(7, 31), [L.xs[31], L.ys[7]]) and L.xs[0] == -1.0 and abs(L.xs[-1] - 24.0) < 1e-9


# ---- the block fit -----------------------------------------------------------------------------------------------------------------

def test_block_fit_recovers_noise_free_emitters_from_displaced_starts():
    P = psf()
    ct = np.array([[8.3, 9.1], [14.6, 12.2]])
    wt = np.array([800.0, 500.0])
    y = P.synth(ct, wt)
    c0 = ct + [[0.4, -0.3], [-0.5, 0.35]]
    c, w, info = IM.block_fit(P, y, 1.0, c0, np.zeros(2), c0 - 2, c0 + 2)
    assert np.abs(c - ct).max() < 1e-6 and np.abs(w - wt).max() < 1e-4
    c1, w1, _ = IM.block_fit(P, y, 1.0, ct, np.zeros(2), ct - 2, ct + 2, joint=False)
    assert np.abs(w1 - wt).max() < 1e-8 and np.abs(c1 - ct).max() == 0     # stage 1 alone: the coefficients at the centres


def test_block_fit_is_non_negative_and_respects_the_box():
    P = psf()
    y = P.synth(np.array([[12.0, 12.0]]), np.array([500.0]))
    c0 = np.array([[12.0, 12.0], [6.0, 6.0]])
    c, w, _ = IM.block_fit(P, y, 1.0, c0, np.zeros(2), c0 - 1, c0 + 1)
    assert (w >= 0).all() and w[1] < 1e-6 * w[0] and abs(w[0] - 500) < 1e-6
    c, w, _ = IM.block_fit(P, y, 1.0, np.array([[11.0, 11.0]]), np.zeros(1), np.array([[10.5, 10.5]]), np.array([[11.5, 11.5]]))
    assert np.all(np.abs(c - 11.5) < 1e-9) or np.all(c <= 11.5 + 1e-9)       # the true centre (12) is outside: the box binds
    assert c.max() <= 11.5 + 1e-9


# ---- the path ----------------------------------------------------------------------------------------------------------------------

def run(P, y, sw=1.0, h=0.5, **kw):
    L = IM.Lattice(P, h=h, margin=1.0)
    return IM.sliding_path(P, y, sw, L, **kw), L


def test_path_single_emitter_noise_free():
    P = psf()
    ct, wt = np.array([[11.7, 9.2]]), np.array([700.0])
    R, L = run(P, P.synth(ct, wt))
    H = R['hist']
    assert len(R['taus']) >= 1 and H[0]['n_ent'] == 0 and len(H[0]['w']) == 0
    # first entrant: the lattice point maximising the correlation with the data = the one nearest the emitter
    b, a = R['events'][0][1]
    assert np.abs(L.pos(b, a) - ct[0]).max() <= 0.25 + 1e-9
    h1 = H[1]
    assert len(h1['w']) == 1 and np.abs(h1['c'] - ct).max() < 1e-6 and abs(h1['w'][0] - 700) < 1e-4 and h1['chi2'] < 1e-14
    assert R['events'][0][2] == 'enter' and len(R['events']) == 1 and len(H) == 2       # exact fit: the path ends


def test_path_orders_by_amplitude_and_recovers_three_separated_emitters():
    P = psf(32)
    ct = np.array([[8.2, 8.4], [22.1, 9.3], [15.4, 23.2]])
    wt = np.array([300.0, 900.0, 1500.0])
    R, L = run(P, P.synth(ct, wt))
    H = R['hist']
    assert R['taus'].size == 3 and np.all(np.diff(R['taus']) > 0)
    # entry order: the brightest first (their correlations are the largest)
    first = [L.pos(*e[1]) for e in R['events'] if e[2] == 'enter']
    order = [int(np.argmin(np.abs(ct - f).sum(1))) for f in first]
    assert order == [2, 1, 0]
    last = H[-1]
    srt = np.lexsort((last['c'][:, 1], last['c'][:, 0]))
    ref = np.lexsort((ct[:, 1], ct[:, 0]))
    assert np.abs(last['c'][srt] - ct[ref]).max() < 1e-5 and np.abs(last['w'][srt] - wt[ref]).max() < 1e-3
    assert last['chi2'] < 1e-12 and last['n_ent'] == 3


def test_path_resolves_two_close_emitters():
    """Two emitters 1.4 σ apart: the first entrant sits between them or on the brighter one; the pair is resolved after the second."""
    P = psf()
    sep = 1.4 * SIG
    ct = np.array([[10.0, 12.0], [10.0 + sep, 12.0]])
    wt = np.array([1000.0, 600.0])
    R, L = run(P, P.synth(ct, wt))
    last = R['hist'][-1]
    srt = np.argsort(last['c'][:, 0])
    assert last['w'].size == 2 and np.abs(last['c'][srt] - ct).max() < 1e-4 and np.abs(last['w'][srt] - wt).max() < 1e-2


def test_ties_enter_together_and_the_symmetric_pair_is_found():
    """Two equal emitters mirrored about the lattice: the duals tie, both kernels enter in one event (same τ)."""
    P = psf()
    ct = np.array([[8.0, 12.0], [16.0, 12.0]])
    R, L = run(P, P.synth(ct, np.array([800.0, 800.0])))
    assert R['taus'].size == 2 and abs(R['taus'][1] - R['taus'][0]) < 1e-9 * R['taus'][0]
    assert len(R['hist']) == 2 and R['hist'][1]['w'].size == 2 and [e[2] for e in R['events']] == ['enter', 'enter']        # one event, one block of two
    last = R['hist'][-1]
    assert np.abs(np.sort(last['c'][:, 0]) - [8.0, 16.0]).max() < 1e-5 and np.abs(last['w'] - 800).max() < 1e-3


def test_entry_time_scales_inversely_with_amplitude():
    """τ_1 = ρ/max_j v_j: doubling the data halves it (the correlation is linear in y); the lattice winner is unchanged."""
    P = psf()
    y = P.synth(np.array([[11.3, 13.1]]), np.array([400.0]))
    R1, _ = run(P, y)
    R2, _ = run(P, 2 * y)
    assert abs(R2['taus'][0] * 2 / R1['taus'][0] - 1) < 1e-12 and R1['events'][0][1] == R2['events'][0][1]
    R3, _ = run(P, y, rho=3.0)
    assert abs(R3['taus'][0] / R1['taus'][0] - 3) < 1e-12


def test_a_kernel_pushed_out_by_the_others_exits():
    """Seed a path whose first entrant (the lattice point between two close emitters, correlation-wise) is later made redundant: the
    exit event is recorded with the kernel's dual back at the wall, and the final block has no zero-coefficient kernel."""
    P = psf()
    ct = np.array([[10.0, 12.0], [11.5, 12.0]])
    R, L = run(P, P.synth(ct, np.array([700.0, 700.0])), h=1.0)
    kinds = [e[2] for e in R['events']]
    last = R['hist'][-1]
    assert (last['w'] > 1e-6).all() and last['chi2'] < 1e-10
    assert kinds.count('enter') - kinds.count('exit') == last['w'].size


def test_path_stops_at_max_entries_and_on_empty_data():
    P = psf()
    rng = np.random.default_rng(4)
    y = 2.0 * rng.standard_normal((24, 24))                    # pure noise: the path never ends by itself
    R, _ = run(P, y, max_entries=6)
    assert R['taus'].size == 6 and R['hist'][-1]['n_ent'] == 6
    R0, _ = run(P, np.zeros((24, 24)))
    assert R0['taus'].size == 0 and len(R0['hist']) == 1


# ---- noise and the stopping rules ------------------------------------------------------------------------------------------------

def noisy(P, ct, wt, sd, seed):
    rng = np.random.default_rng(seed)
    return P.synth(ct, wt) + sd * rng.standard_normal((P.ny, P.nx))


def test_discrepancy_and_gap_rules_find_the_true_count_in_gaussian_noise():
    P = psf(32)
    ct = np.array([[8.3, 8.4], [22.1, 9.3], [15.4, 23.2], [24.5, 24.1]])
    wt = np.array([900.0, 1400.0, 1000.0, 700.0])
    sd = 6.0
    hits_d, hits_g = 0, 0
    for seed in range(5):
        y = noisy(P, ct, wt, sd, seed)
        R, _ = run(P, y, sw=1.0 / sd, max_entries=12)
        i = IM.stop_discrepancy(R['hist'], P.N, z=2.0)
        hits_d += len(R['hist'][i]['w']) == 4
        j, K, ratio = IM.stop_gap(R['hist'], R['taus'])
        hits_g += (K == 4 and len(R['hist'][j]['w']) == 4)
        h = R['hist'][i]
        srt = np.lexsort((h['c'][:, 1], h['c'][:, 0]))
        ref = np.lexsort((ct[:, 1], ct[:, 0]))
        if len(h['w']) == 4:
            assert np.abs(h['c'][srt] - ct[ref]).max() < 0.15 and np.abs(h['w'][srt] / wt[ref] - 1).max() < 0.15
    assert hits_d == 5 and hits_g == 5


def test_discrepancy_rule_on_pure_noise_returns_the_empty_block():
    P = psf()
    rng = np.random.default_rng(5)
    sd = 3.0
    y = sd * rng.standard_normal((24, 24))
    R, _ = run(P, y, sw=1.0 / sd, max_entries=5)
    assert IM.stop_discrepancy(R['hist'], P.N, z=3.0) == 0


def test_stop_gap_edge_cases():
    hist = [dict(n_ent=i, w=np.zeros(i), c=np.zeros((i, 2)), chi2=1.0) for i in range(5)]
    assert IM.stop_gap(hist, [1.0], k_min=1) == (0, 0, 1.0)
    idx, K, ratio = IM.stop_gap(hist, [1.0, 1.1, 1.2, 9.0], k_min=1)
    assert K == 3 and idx == 3 and abs(ratio - 7.5) < 1e-12                  # the jump between the 3rd and 4th entry: keep 3
    idx, K, ratio = IM.stop_gap(hist, [1.0, 5.0, 5.5, 6.0], k_min=1)
    assert K == 1 and idx == 1 and abs(ratio - 5.0) < 1e-12
    idx, K, ratio = IM.stop_gap(hist, [1.0, 5.0, 5.5, 6.0], k_min=2)         # k_min: the first jump (1 → 2) is not eligible
    assert K == 2 and abs(ratio - 1.1) < 1e-12
    # the discrepancy rule on a hand-made history: dof correction and the z-band
    N = 100
    hh = [dict(w=np.zeros(k), chi2=c) for k, c in ((0, 5.0), (1, 1.6), (2, 1.0), (3, 0.9))]
    assert IM.stop_discrepancy(hh, N, z=2.0) == 2
    assert IM.stop_discrepancy(hh, N, z=2.0, chi_target=0.01) == 3           # none within the band: the last block


def test_reweighting_from_the_model_runs_and_keeps_the_fit():
    P = psf(32)
    ct = np.array([[10.2, 12.4], [21.7, 20.1]])
    wt = np.array([900.0, 600.0])
    rng = np.random.default_rng(6)
    bg = 5.0
    y = rng.poisson(bg + P.synth(ct, wt)).astype(float) - bg
    L = IM.Lattice(P, h=0.5, margin=1.0)
    fn = lambda model: 1.0 / np.sqrt(bg + np.maximum(model, 0.0))
    R = IM.sliding_path(P, y, fn(0 * y), L, max_entries=8, sw_fn=fn)
    i = IM.stop_discrepancy(R['hist'], P.N, z=3.0)
    h = R['hist'][i]
    assert len(h['w']) == 2
    srt = np.lexsort((h['c'][:, 1], h['c'][:, 0]))
    ref = np.lexsort((ct[:, 1], ct[:, 0]))
    assert np.abs(h['c'][srt] - ct[ref]).max() < 0.1


def test_simulate_frame_statistics():
    P = psf(32)
    rng = np.random.default_rng(7)
    c = np.array([[16.0, 16.0]])
    fr = np.array([IM.simulate_frame(P, c, [1000.0], 3.0, rng) for _ in range(300)])
    mean = 3.0 + P.synth(c, np.array([1000.0]))
    assert abs(fr.mean(0)[16, 16] / mean[16, 16] - 1) < 0.02
    assert np.abs(fr.var(0)[14:19, 14:19] / mean[14:19, 14:19] - 1).max() < 0.35


def test_emccd_noise_model_matches_its_variance_formula():
    """The exact EM-CCD simulation (gamma register, read-out, quantisation) in photon units: mean = photons, variance =
    (2/QE)·photons + read-out² (the 'excess noise' factor 2 of the EM register, divided by the quantum efficiency)."""
    P = psf(32)
    rng = np.random.default_rng(8)
    c = np.array([[16.0, 16.0]])
    fr = np.array([IM.simulate_emccd(P, c, [2000.0], 2.0, rng) for _ in range(1500)])
    mean = 2.0 + P.synth(c, np.array([2000.0]))
    for (i, j) in ((16, 16), (16, 18), (12, 12)):
        assert abs(fr[:, i, j].mean() / mean[i, j] - 1) < 0.03
        assert abs(fr[:, i, j].var() / IM.emccd_variance(mean[i, j]) - 1) < 0.15
    assert abs(IM.emccd_variance(0.0) - (74.4 / 45 / 6.0) ** 2) < 1e-12 and abs(IM.emccd_variance(1.0) - IM.emccd_variance(0.0) - 2 / 0.9) < 1e-12


# ---- widths (moving widths, Proposition D‴.59) -------------------------------------------------------------------------------------

def test_width_derivative_and_three_column_jacobian_by_finite_differences():
    rng = np.random.default_rng(9)
    P = psf()
    c = np.c_[rng.uniform(4, 19, (3, 2)), rng.uniform(0.9, 2.2, 3)]
    w = rng.uniform(1, 4, 3)
    sw = rng.uniform(0.5, 1.5, (24, 24))
    y = rng.standard_normal((24, 24))
    k = 3
    ph = P.phi(c)
    assert np.abs(ph - P.phi(c[:, :2])).max() > 1e-3
    h = 1e-6
    f = P.factors(c, deriv_s=True)
    d = np.zeros_like(c); d[:, 2] = h
    gy_p, gx_p = P.factors(c + d)[:2]
    gy_m, gx_m = P.factors(c - d)[:2]
    assert np.abs((gx_p - gx_m) / (2 * h) - f[5]).max() < 1e-8 and np.abs((gy_p - gy_m) / (2 * h) - f[4]).max() < 1e-8
    z = np.r_[w, c.T.ravel()]

    def res(z):
        return (sw * (y - P.synth(z[k:].reshape(3, k).T, z[:k]))).ravel()
    J = P.jac(c, w, sw)
    assert J.shape == (576, 12)
    Jn = np.empty_like(J)
    for i in range(12):
        e = np.zeros(12); e[i] = 1e-6
        Jn[:, i] = (res(z + e) - res(z - e)) / 2e-6
    assert np.abs(J - Jn).max() < 1e-6


def test_block_fit_recovers_a_broad_emitter_with_its_own_width():
    P = psf(32)
    ct = np.array([[15.3, 16.1, 2.1]])
    wt = np.array([900.0])
    y = P.synth(ct, wt)
    c0 = np.array([[14.9, 16.5, SIG]])
    lo, hi = c0 - 2, c0 + 2
    lo[:, 2], hi[:, 2] = 0.8, 3.5
    c, w, info = IM.block_fit(P, y, 1.0, c0, np.zeros(1), lo, hi)
    assert np.abs(c - ct).max() < 1e-6 and abs(w[0] - 900) < 1e-3
    hi[:, 2] = 1.5                                               # the bound binds
    c, w, _ = IM.block_fit(P, y, 1.0, c0, np.zeros(1), lo, hi)
    assert c[0, 2] <= 1.5 + 1e-9


def test_path_with_sliding_widths_fits_a_defocused_and_a_sharp_emitter_exactly():
    P = psf(40)
    ct = np.array([[12.2, 12.6, 1.2], [27.4, 24.3, 2.3]])
    wt = np.array([1500.0, 2500.0])
    y = P.synth(ct, wt)
    L = IM.Lattice(P, h=0.5, margin=1.0)
    R = IM.sliding_path(P, y, 1.0, L, max_entries=6, sig_bounds=(0.8, 3.5))
    last = R['hist'][-1]
    assert last['c'].shape == (2, 3) and last['chi2'] < 1e-12
    srt = np.argsort(last['c'][:, 0])
    assert np.abs(last['c'][srt] - ct).max() < 1e-5 and np.abs(last['w'][srt] - wt).max() < 1e-2
    # a fixed width cannot fit the same frame: a residual stays (the model-class mismatch the widths remove)
    R2 = IM.sliding_path(P, y, 1.0, L, max_entries=2)
    assert R2['hist'][-1]['chi2'] > 1e-6 * R['hist'][0]['chi2']


# ---- parts (halo) and the high-pass -------------------------------------------------------------------------------------------------

def test_parts_are_a_sum_of_gaussians_of_scaled_widths():
    P = IM.PixelPSF(30, 30, 1.2, parts=((0.7, 1.0), (0.3, 3.0)))
    G1, G3 = IM.PixelPSF(30, 30, 1.2), IM.PixelPSF(30, 30, 3.6)
    c = np.array([[14.3, 15.2], [8.0, 20.5]])
    assert np.abs(P.phi(c) - (0.7 * G1.phi(c) + 0.3 * G3.phi(c))).max() < 1e-14
    cw = np.c_[c, [1.5, 2.0]]                                    # own widths: both parts scale with the kernel's width
    ref = 0.7 * IM.PixelPSF(30, 30, 1.0).phi(np.c_[c, [1.5, 2.0]]) + 0.3 * IM.PixelPSF(30, 30, 1.0).phi(np.c_[c, 3 * np.array([1.5, 2.0])])
    assert np.abs(P.phi(cw) - ref).max() < 1e-14


def test_jacobian_with_parts_and_widths_by_finite_differences():
    rng = np.random.default_rng(10)
    P = IM.PixelPSF(24, 24, 1.2, parts=((0.7, 1.0), (0.3, 3.0)))
    c = np.c_[rng.uniform(5, 18, (2, 2)), rng.uniform(1.0, 2.0, 2)]
    w = rng.uniform(1, 4, 2)
    sw = rng.uniform(0.5, 1.5, (24, 24))
    y = rng.standard_normal((24, 24))
    k = 2
    z = np.r_[w, c.T.ravel()]

    def res(z):
        return (sw * (y - P.synth(z[k:].reshape(3, k).T, z[:k]))).ravel()
    J = P.jac(c, w, sw)
    Jn = np.empty_like(J)
    for i in range(J.shape[1]):
        e = np.zeros(J.shape[1]); e[i] = 1e-6
        Jn[:, i] = (res(z + e) - res(z - e)) / 2e-6
    assert np.abs(J - Jn).max() < 1e-6


def test_highpass_model_is_the_filtered_plain_model_also_at_the_border():
    Pp = IM.PixelPSF(28, 28, 1.3)
    Ph = IM.PixelPSF(28, 28, 1.3, hp=3.0)
    c = np.array([[2.1, 3.4, 1.3], [14.0, 13.5, 1.9], [26.2, 24.7, 1.5]])      # near two borders and in the middle, own widths
    w = np.array([500.0, 900.0, 300.0])
    assert np.abs(Ph.synth(c, w) - Ph.highpass(Pp.synth(c, w))).max() < 1e-10
    # the filter kills a constant (its rows sum to one) and passes nothing without a high-pass
    assert np.abs(Ph.highpass(np.full((28, 28), 7.0))).max() < 1e-12
    y = np.random.default_rng(0).standard_normal((28, 28))
    assert Pp.highpass(y) is y and np.abs(Ph.highpass(y) - (y - Ph.Sy @ y @ Ph.Sx.T)).max() == 0
    # the smoother is row-stochastic and symmetric in the interior
    S = IM.smoother(28, 3.0)
    assert np.abs(S.sum(1) - 1).max() < 1e-14 and S[0, 0] > S[0, 5] > 0 and S[0, 13] == 0


def test_highpass_corr_jacobian_window_and_lattice_are_consistent():
    rng = np.random.default_rng(11)
    P = IM.PixelPSF(30, 30, 1.2, hp=3.0, parts=((0.8, 1.0), (0.2, 3.0)))
    c = np.c_[rng.uniform(4, 25, (3, 2)), rng.uniform(1.0, 2.0, 3)]
    w = rng.uniform(1, 4, 3)
    Rw = rng.standard_normal((30, 30))
    Phi = P.phi(c)
    assert np.abs(Phi.T @ Rw.ravel() - P.corr(c, Rw)[0]).max() < 1e-10
    assert np.abs(P.synth(c, w).ravel() - Phi @ w).max() < 1e-12
    # window: the cropped full-frame model
    pw = P.window(5, 21, 8, 27)
    assert np.abs(pw.synth(c, w) - P.synth(c, w)[5:21, 8:27]).max() < 1e-12
    pw2 = pw.window(2, 9, 3, 10)                                  # a window of a window
    assert np.abs(pw2.synth(c, w) - P.synth(c, w)[7:14, 11:18]).max() < 1e-12
    # lattice correlation at the PSF's own width
    L = IM.Lattice(P, h=0.5, margin=1.0)
    V = L.corr(Rw)
    for b, a in ((0, 0), (10, 20), (40, 55)):
        assert abs(V[b, a] - P.corr(np.r_[L.pos(b, a), 1.2][None, :], Rw)[0][0]) < 1e-10


def test_path_with_a_highpass_psf_recovers_noise_free_emitters_and_ignores_a_pedestal():
    P = IM.PixelPSF(36, 36, 1.3, hp=4.0)
    ct = np.array([[10.3, 11.1, 1.2], [25.4, 22.6, 1.8]])
    wt = np.array([1200.0, 2000.0])
    y = P.synth(ct, wt)
    L = IM.Lattice(P, h=0.5, margin=1.0)
    R = IM.sliding_path(P, y, 1.0, L, max_entries=5, sig_bounds=(0.8, 3.0))
    last = R['hist'][-1]
    srt = np.argsort(last['c'][:, 0])
    assert last['w'].size == 2 and np.abs(last['c'][srt] - ct).max() < 1e-5 and np.abs(last['w'][srt] - wt).max() < 1e-2
    # a constant pedestal on the raw frame is removed by the filter: the filtered data are unchanged
    Pp = IM.PixelPSF(36, 36, 1.3)
    raw = Pp.synth(ct, wt)
    assert np.abs(P.highpass(raw + 50.0) - P.highpass(raw)).max() < 1e-12


def test_exclusion_zone_keeps_later_entrants_away_from_active_kernels():
    """A broad emitter fitted with a too-narrow fixed PSF leaves rings and shoulders that pull extra kernels in next to the first;
    with r_ex no kernel enters within r_ex of an active centre (and without it some do: the test is not vacuous)."""
    P = psf(40)
    y = IM.PixelPSF(40, 40, 2.2).synth(np.array([[20.2, 19.7]]), np.array([3000.0]))
    L = IM.Lattice(P, h=0.5, margin=1.0)

    def min_dists(R):
        out = []
        ents = [e for e in R['events'] if e[2] == 'enter']
        for n, e in enumerate(ents[1:], start=1):
            prev = R['hist'][n]['c']                                  # the block before the (n+1)-th entry (no exits in this frame)
            out.append(np.linalg.norm(prev - L.pos(*e[1]), axis=1).min())
        return np.array(out)
    R0 = IM.sliding_path(P, y, 1.0, L, max_entries=8)
    R1 = IM.sliding_path(P, y, 1.0, L, max_entries=8, r_ex=2.5)
    assert all(e[2] == 'enter' for e in R0['events']) and all(e[2] == 'enter' for e in R1['events'])
    assert min_dists(R0).min() < 2.5
    assert min_dists(R1).min() >= 2.5 - 1e-9


# ---- the significance rule and the gated gap rule ------------------------------------------------------------------------------

def test_entrant_significance_is_the_matched_filter_snr_and_the_rules_use_it():
    P = psf(32)
    ct = np.array([[8.3, 8.4], [22.1, 9.3], [15.4, 23.2]])
    wt = np.array([2500.0, 900.0, 300.0])
    sd = 6.0
    y = noisy(P, ct, wt, sd, 3)
    L = IM.Lattice(P, h=0.5, margin=1.0)
    R = IM.sliding_path(P, y, 1.0 / sd, L, max_entries=10)
    snr = R['snr']
    assert snr.shape == R['taus'].shape
    # first entrant: snr = (Σ y φ_j / sd²)/ (√(Σ φ_j²)/sd) = Σ y φ_j / (sd √Σφ_j²): recompute directly
    b, a = R['events'][0][1]
    ph = P.phi(L.pos(b, a)[None, :])[:, 0]
    assert abs(snr[0] - (y.ravel() @ ph) / (sd * np.linalg.norm(ph))) < 1e-9
    # the strong emitters are far above 4.5, the pure-noise entrants are below it
    assert (snr[:2] > 8).all() and snr[3:].max() < 4.5
    i, K = IM.stop_snr(R['hist'], R['taus'], snr, kappa=4.5)
    assert K == 3 and R['hist'][i]['n_ent'] == 3 and len(R['hist'][i]['w']) == 3
    assert IM.stop_snr(R['hist'], R['taus'], snr, kappa=1e9)[1] == 0 and IM.stop_snr(R['hist'], R['taus'], snr, kappa=0.0)[1] == len(snr)
    # gated gap rule: an empty frame is called empty, the gap rule alone would call one
    rng = np.random.default_rng(12)
    y0 = sd * rng.standard_normal((32, 32))
    R0 = IM.sliding_path(P, y0, 1.0 / sd, L, max_entries=6)
    assert IM.stop_gap(R0['hist'], R0['taus'])[1] >= 1
    assert IM.stop_gap(R0['hist'], R0['taus'], snr=R0['snr'], gate=4.5)[1] == 0
    assert IM.stop_gap(R['hist'], R['taus'], snr=snr, gate=4.5)[1] >= 1


def test_snr_stop_cuts_the_noise_tail_after_the_first_insignificant_entrant():
    P = psf(32)
    ct = np.array([[8.3, 8.4], [22.1, 9.3]])
    wt = np.array([2500.0, 1500.0])
    sd = 6.0
    y = noisy(P, ct, wt, sd, 5)
    L = IM.Lattice(P, h=0.5, margin=1.0)
    full = IM.sliding_path(P, y, 1.0 / sd, L, max_entries=12)
    cut = IM.sliding_path(P, y, 1.0 / sd, L, max_entries=12, snr_min=4.5, snr_extra=2)
    k = int(np.flatnonzero(full['snr'] < 4.5)[0])
    assert k == 2 and len(cut['taus']) == k + 1 + 2 and np.allclose(cut['taus'], full['taus'][:k + 3], rtol=1e-6)
    assert IM.stop_snr(cut['hist'], cut['taus'], cut['snr'], 4.5)[1] == IM.stop_snr(full['hist'], full['taus'], full['snr'], 4.5)[1] == 2


# ---- the greedy order (matching pursuit with sliding) -------------------------------------------------------------------------------

def test_norm2_is_the_weighted_squared_norm_of_every_lattice_kernel():
    rng = np.random.default_rng(13)
    P = IM.PixelPSF(24, 24, 1.2, hp=3.0, parts=((0.8, 1.0), (0.2, 3.0)))
    L = IM.Lattice(P, h=1.0, margin=1.0)
    W = rng.uniform(0.5, 2.0, (24, 24))
    N2 = L.norm2(W)
    for b, a in ((0, 0), (5, 9), (25, 25), (13, 2)):
        ph = P.phi(L.pos(b, a)[None, :])[:, 0]
        assert abs(N2[b, a] - np.sum(W.ravel() * ph ** 2)) < 1e-12 * max(1.0, N2[b, a])


def test_greedy_order_picks_the_most_significant_candidate_and_recovers_separated_emitters():
    P = psf(32)
    ct = np.array([[8.3, 8.4], [22.1, 9.3], [15.4, 23.2]])
    wt = np.array([300.0, 900.0, 1500.0])
    sd = 6.0
    y = noisy(P, ct, wt, sd, 2)
    L = IM.Lattice(P, h=0.5, margin=1.0)
    Rg = IM.sliding_path(P, y, 1.0 / sd, L, max_entries=6, order='greedy')
    Rd = IM.sliding_path(P, y, 1.0 / sd, L, max_entries=6)
    assert np.allclose(Rg['taus'], np.arange(1, 7))                       # the step count
    snr = Rg['snr']
    assert snr[0] >= snr[1] >= snr[2] > 4.5 and snr[3] < 4.5            # three emitters, then noise
    # the first candidate is the maximum of v/√(ΣWφ²) over the lattice
    W = np.full((32, 32), 1.0 / sd ** 2)
    z = (L.corr(W * y) / np.sqrt(L.norm2(W))).max()
    assert abs(snr[0] - z) < 1e-9
    ig, Kg = IM.stop_snr(Rg['hist'], Rg['taus'], Rg['snr'], 4.5)
    idd, Kd = IM.stop_snr(Rd['hist'], Rd['taus'], Rd['snr'], 4.5)
    assert Kg == Kd == 3
    cg, cd = Rg['hist'][ig]['c'], Rd['hist'][idd]['c']
    assert np.abs(cg[np.lexsort((cg[:, 1], cg[:, 0]))] - cd[np.lexsort((cd[:, 1], cd[:, 0]))]).max() < 1e-3


def test_fit_radius_local_slides_agree_with_the_global_fit_in_a_crowded_field():
    """Thirty well-separated-by-1.7-widths emitters in a 48 x 48 frame: refitting only the neighbours of each entrant (4 px) gives the same final block
    as the transitive cluster (which here spans the whole frame), noise-free, after a final polish of everything."""
    P = psf(48)
    rng = np.random.default_rng(14)
    ct = []
    while len(ct) < 30:
        c = rng.uniform(3, 45, 2)
        if all(np.linalg.norm(c - q) > 2.0 for q in ct):
            ct.append(c)
    ct = np.array(ct)
    wt = rng.uniform(500, 2000, 30)
    y = P.synth(ct, wt)
    L = IM.Lattice(P, h=0.5, margin=1.0)
    Rl = IM.sliding_path(P, y, 1.0, L, max_entries=34, fit_radius=4.0, tol=1e-9)
    last = Rl['hist'][-1]
    c, w, _ = IM.refine(P, y, 1.0, last['c'], last['w'], box=0.5)
    d = np.linalg.norm(c[:, None, :] - ct[None, :, :], axis=2)
    assert len(last['w']) >= 30
    assert d.min(0).max() < 1e-3                                    # every true emitter has a unit within 1e-3 px
    assert np.mean((y - P.synth(c, w)) ** 2) < 1e-8


def test_lattice_near_mask_equals_the_dense_distance_computation():
    P = psf(30)
    L = IM.Lattice(P, h=0.5, margin=1.0)
    rng = np.random.default_rng(15)
    c = np.c_[rng.uniform(-1.5, 30.0, 25), rng.uniform(-1.5, 30.0, 25)]            # some beyond the lattice
    for r in (0.7, 2.0, 3.3):
        dx = L.xs[None, :, None] - c[:, 0][None, None, :]
        dy = L.ys[:, None, None] - c[:, 1][None, None, :]
        dense = ((dx ** 2 + dy ** 2) < r ** 2).any(axis=2)
        assert (L.near(c, r) == dense).all()
    assert not L.near(np.zeros((0, 2)), 2.0).any()


def test_fit_psf_parts_recovers_a_known_mixture_and_fits_a_moffat_stamp():
    Q = IM.PixelPSF(41, 41, 1.1, parts=((0.7, 1.0), (0.2, 2.6), (0.1, 6.0)))
    stamp = Q.phi(np.array([[20.0, 20.0]]))[:, 0].reshape(41, 41)
    sg, parts, rms = IM.fit_psf_parts(stamp, n_parts=3, sigma0=1.0)
    assert abs(sg - 1.1) < 5e-3 and rms < 1e-4
    assert abs(sum(a for a, k in parts) - 1) < 1e-12 and parts[0][1] == 1.0
    ref = Q.phi(np.array([[20.0, 20.0]]))[:, 0]
    fit = IM.PixelPSF(41, 41, sg, parts=parts).phi(np.array([[20.0, 20.0]]))[:, 0]
    assert np.abs(fit - ref).max() / ref.max() < 2e-3
    # a Moffat profile with power-law wings: three Gaussians reproduce it to a few per cent of the peak
    yy, xx = np.mgrid[0:41, 0:41]
    mof = (1 + ((xx - 20) ** 2 + (yy - 20) ** 2) / (2.2 ** 2 * 3.5)) ** -3.5
    sg2, parts2, rms2 = IM.fit_psf_parts(mof, n_parts=3, sigma0=1.0)
    assert rms2 < 0.02 and sg2 > 0.5 and all(k >= 1.0 for a, k in parts2)


def test_stopping_time_keeps_the_entrants_of_average_significance_and_stops_the_path():
    P = psf(32)
    ct = np.array([[8.3, 8.4], [22.1, 9.3], [15.4, 23.2]])
    wt = np.array([2500.0, 1500.0, 700.0])
    sd = 6.0
    y = noisy(P, ct, wt, sd, 7)
    L = IM.Lattice(P, h=0.5, margin=1.0)
    scale = IM.noise_scale(L, np.full((32, 32), 1.0 / sd ** 2))
    interior = np.sqrt(np.sum(P.phi(np.array([[16.0, 16.0]])) ** 2) / sd ** 2)                 # the squared norm of a kernel away from the border, over sd
    assert abs(scale - interior) < 1e-5 * interior and 0.03 < scale < 0.05
    full = IM.sliding_path(P, y, 1.0 / sd, L, max_entries=8)
    kappa = 5.0
    tau_stop = 1.0 / (kappa * scale)
    i, K = IM.stop_time(full['hist'], full['taus'], tau_stop)
    assert K == 3 and full['hist'][i]['n_ent'] == 3
    cut = IM.sliding_path(P, y, 1.0 / sd, L, max_entries=8, tau_max=tau_stop)
    assert len(cut['taus']) == 3 and np.allclose(cut['taus'], full['taus'][:3], rtol=1e-9)
    assert IM.stop_time(full['hist'], full['taus'], 0.0)[1] == 0 and IM.stop_time(full['hist'], full['taus'], 1e9)[1] == len(full['taus'])


def test_fit_frame_recovers_sources_of_an_emccd_frame_with_every_stop_rule():
    P = IM.PixelPSF(64, 64, 1.4, parts=((0.75, 1.0), (0.25, 3.0)), hp=5.0)
    Pp = IM.PixelPSF(64, 64, 1.4, parts=((0.75, 1.0), (0.25, 3.0)))
    L = IM.Lattice(P, h=0.5, margin=1.0)
    ct = np.array([[14.3, 20.1], [40.6, 33.2], [27.4, 50.7]])
    wt = np.array([5000.0, 3000.0, 4000.0])
    rng = np.random.default_rng(21)
    frame = IM.simulate_emccd(Pp, ct, wt, 100.0, rng)
    for stop in (('snr', 5.0), ('time', 5.0), ('disc', 3.0), ('gap', 5.0), ('bic', 4)):
        out = IM.fit_frame(frame, P, L, 100.0, 2.236, 0.078, widths=(1.0, 3.2), stop=stop, max_entries=12)
        c, w = out['c'], out['w']
        assert len(w) == 3
        d = np.linalg.norm(c[:, None, :2] - ct[None, :, :], axis=2)
        assert d.min(0).max() < 0.35 and d.min(1).max() < 0.35 and np.all(w > 1500)
    g = IM.fit_frame(frame, P, L, 100.0, 2.236, 0.078, widths=(1.0, 3.2), stop=('snr', 5.0), max_entries=12, order='greedy')      # the instantaneous correlation instead of the accumulated dual
    d = np.linalg.norm(g['c'][:, None, :2] - ct[None, :, :], axis=2)
    assert len(g['w']) == 3 and d.min(0).max() < 0.35 and d.min(1).max() < 0.35
    # a fixed variance map and no widths run as well; an empty frame gives nothing
    var = 2.236 * 100.0 * np.ones((64, 64))
    out = IM.fit_frame(frame, P, L, 100.0, var, 0.078, widths=None, stop=('snr', 5.0), max_entries=12)
    assert out['c'].shape[1] == 2 and len(out['w']) == 3
    empty = IM.simulate_emccd(Pp, np.zeros((0, 2)), np.zeros(0), 100.0, np.random.default_rng(22))
    assert len(IM.fit_frame(empty, P, L, 100.0, 2.236, 0.078, widths=(1.0, 3.2), stop=('snr', 5.0), max_entries=8)['w']) == 0
    with pytest.raises(ValueError):
        IM.fit_frame(frame, P, L, 100.0, 2.236, 0.078, stop=('nonsense', 1))


def test_unit_significance_is_amplitude_over_the_noise_it_would_have_alone():
    P = psf(32)
    c = np.array([[10.3, 12.4], [21.7, 20.1]])
    w = np.array([400.0, 0.0])
    sd = 5.0
    z = IM.unit_significance(P, 1.0 / sd, c, w)
    ph = P.phi(c)
    assert abs(z[0] - 400.0 * np.linalg.norm(ph[:, 0]) / sd) < 1e-12 and z[1] == 0.0
    assert IM.unit_significance(P, 1.0 / sd, np.zeros((0, 2)), np.zeros(0)).size == 0
    # the weights enter squared and chunks do not change the answer
    sw = np.random.default_rng(1).uniform(0.1, 0.3, (32, 32))
    cc = np.random.default_rng(2).uniform(3, 28, (7, 2))
    ww = np.random.default_rng(3).uniform(100, 900, 7)
    full = IM.unit_significance(P, sw, cc, ww)
    assert np.allclose(full, IM.unit_significance(P, sw, cc, ww, chunk=3), rtol=1e-12)
    assert np.allclose(full, ww * np.sqrt((sw.ravel() ** 2) @ P.phi(cc) ** 2), rtol=1e-12)


def test_stop_bic_keeps_a_unit_iff_it_lowers_the_weighted_sum_of_squares_by_more_than_n_par_ln_n():
    N = 4096
    pen = 4 * np.log(N)                                                       # 33.27 for a unit of four parameters
    chi2 = lambda ss: ss / N                                                   # hist stores the weighted MEAN square residual
    ss = [1000.0, 500.0, 440.0, 420.0, 415.0]                                  # sums of squares after 0..4 units: drops 500, 60, 20, 5
    hist = [dict(w=np.zeros(k), chi2=chi2(v)) for k, v in enumerate(ss)]
    assert IM.stop_bic(hist, N, n_par=4) == (2, 2)                            # the third unit (drop 20 < 33.3) is not worth its parameters
    ss2 = [1000.0, 1000.0 - pen - 1e-6, 1000.0 - 2 * pen - 2e-6]               # every unit drops the sum of squares by a hair more than the penalty
    assert IM.stop_bic([dict(w=np.zeros(k), chi2=chi2(v)) for k, v in enumerate(ss2)], N, n_par=4)[1] == 2
    ss3 = [1000.0, 1000.0 - pen + 1e-6]                                        # a hair less: not kept
    assert IM.stop_bic([dict(w=np.zeros(k), chi2=chi2(v)) for k, v in enumerate(ss3)], N, n_par=4) == (0, 0)
    ss4 = [1000.0, 1000.0 - pen]                                               # exactly equal costs: the smaller model (the first minimum)
    assert IM.stop_bic([dict(w=np.zeros(k), chi2=chi2(v)) for k, v in enumerate(ss4)], N, n_par=4)[1] == 0
    hist28 = [dict(w=np.zeros(k), chi2=chi2(v)) for k, v in enumerate([1000.0, 972.0])]
    assert IM.stop_bic(hist28, N, n_par=3)[1] == 1 and IM.stop_bic(hist28, N, n_par=4)[1] == 0      # a drop of 28 beats 3 ln N = 24.9 but not 4 ln N
    assert IM.stop_bic([dict(w=np.zeros(k), chi2=v / 100) for k, v in enumerate([1000.0, 972.0])], 100, n_par=4)[1] == 1      # the penalty grows with ln N: 4 ln 100 = 18.4 < 28
    hist_k = [dict(w=np.zeros(k), chi2=chi2(v)) for k, v in enumerate([1000.0, 900.0, 800.0])]
    assert IM.stop_bic(hist_k, N, n_par=4) == (2, 2)                          # K̂ is the number of units of the chosen block, not the index


def test_stop_bic_on_a_synthetic_frame_finds_the_sources_and_nothing_in_noise():
    P = psf(32)
    ct = np.array([[10.2, 12.4], [21.7, 20.1]])
    wt = np.array([900.0, 600.0])
    rng = np.random.default_rng(6)
    bg = 5.0
    fn = lambda model: 1.0 / np.sqrt(bg + np.maximum(model, 0.0))
    L = IM.Lattice(P, h=0.5, margin=1.0)
    y = rng.poisson(bg + P.synth(ct, wt)).astype(float) - bg
    R = IM.sliding_path(P, y, fn(0 * y), L, max_entries=8, sw_fn=fn)
    assert IM.stop_bic(R['hist'], P.N, n_par=3)[1] == 2
    sd = 3.0
    yn = sd * np.random.default_rng(5).standard_normal((32, 32))
    Rn, _ = run(P, yn, sw=1.0 / sd, max_entries=5)
    assert IM.stop_bic(Rn['hist'], P.N, n_par=3)[1] == 0                      # pure noise: no unit pays for its parameters


def test_stop_bic_returns_the_unit_count_of_the_block_not_its_index():
    N = 4096
    # blocks along a path with an exit: 0, 1, 2, 1 units; the cheapest is the last block (index 3, one unit)
    sizes = [0, 1, 2, 1]
    ss = [1000.0, 500.0, 499.0, 480.0]
    hist = [dict(w=np.zeros(k), chi2=v / N) for k, v in zip(sizes, ss)]
    assert IM.stop_bic(hist, N, n_par=3) == (3, 1)


def test_fit_frame_bic_uses_its_own_parameter_count_and_the_orders_differ_in_their_times(monkeypatch):
    P = IM.PixelPSF(64, 64, 1.4, parts=((0.75, 1.0), (0.25, 3.0)), hp=5.0)
    Pp = IM.PixelPSF(64, 64, 1.4, parts=((0.75, 1.0), (0.25, 3.0)))
    L = IM.Lattice(P, h=0.5, margin=1.0)
    ct = np.array([[14.3, 20.1], [40.6, 33.2], [27.4, 50.7]])
    wt = np.array([5000.0, 3000.0, 4000.0])
    frame = IM.simulate_emccd(Pp, ct, wt, 100.0, np.random.default_rng(21))
    kw = dict(widths=(1.0, 3.2), max_entries=12)
    assert len(IM.fit_frame(frame, P, L, 100.0, 2.236, 0.078, stop=('bic', 4), **kw)['w']) == 3
    assert len(IM.fit_frame(frame, P, L, 100.0, 2.236, 0.078, stop=('bic', 1e5), **kw)['w']) == 0      # a huge price per parameter keeps nothing
    dual = IM.fit_frame(frame, P, L, 100.0, 2.236, 0.078, stop=('snr', 5.0), **kw)
    greedy = IM.fit_frame(frame, P, L, 100.0, 2.236, 0.078, stop=('snr', 5.0), order='greedy', **kw)
    assert np.allclose(greedy['taus'], np.arange(1, len(greedy['taus']) + 1))                            # greedy: the step count
    assert not np.allclose(dual['taus'], np.arange(1, len(dual['taus']) + 1))                            # dual: entry times of the accumulated duals
    calls = []
    orig = IM.sliding_path
    monkeypatch.setattr(IM, 'sliding_path', lambda *a, **k: (calls.append(k), orig(*a, **k))[1])
    IM.fit_frame(frame, P, L, 100.0, 2.236, 0.078, stop=('bic', 4), **kw)
    IM.fit_frame(frame, P, L, 100.0, 2.236, 0.078, stop=('snr', 5.0), **kw)
    assert calls[0]['snr_min'] == 3.0 and calls[1]['snr_min'] == 3.0                                    # the path runs to the noise tail before the BIC / significance rule cuts it


def test_greedy_order_normalises_the_correlation_by_the_weights():
    P = psf(32)
    L = IM.Lattice(P, h=0.5, margin=1.0)
    ct = np.array([[8.3, 16.4], [23.6, 15.7]])
    # source A (left, weight 1) with flux 10 units, source B (right, weight 3) with flux 2: the weighted correlation is larger for B (2 x 9 > 10 x 1), the significance for A (2 x 3 < 10 x 1)
    sw = np.ones((32, 32))
    sw[:, 16:] = 3.0
    y = P.synth(ct, np.array([10.0, 2.0]))
    Rg = IM.sliding_path(P, y, sw, L, max_entries=1, order='greedy')
    z = (L.corr(sw ** 2 * y) / np.sqrt(L.norm2(sw ** 2))).max()
    assert abs(Rg['snr'][0] - z) < 1e-9
    assert np.linalg.norm(Rg['hist'][1]['c'][0] - ct[0]) < 0.6                                           # A enters first
    raw = (L.corr(sw ** 2 * y)).max()
    assert raw > 0 and (L.corr(sw ** 2 * y) / np.sqrt(L.norm2(sw ** 2))).argmax() != L.corr(sw ** 2 * y).argmax()    # the two criteria really disagree on this frame



def test_smoother_support_is_four_widths_and_its_rows_are_normalised_gaussians():
    S = IM.smoother(60, 2.0)
    assert np.allclose(S.sum(1), 1.0)
    g = np.exp(-0.5 * (np.arange(-8, 9) / 2.0) ** 2)
    assert np.allclose(S[30, 22:39], g / g.sum(), atol=1e-15) and S[30, 21] == 0 and S[30, 39] == 0                  # truncated at 4 sb = 8 pixels


def test_reach_of_a_highpassed_kernel_is_the_root_of_the_sum_of_squares():
    parts = ((0.5, 1.0), (0.5, 3.0))
    P = IM.PixelPSF(10, 10, 0.8, parts=parts, hp=3.0)
    assert abs(P.reach() - np.hypot(2.4, 3.0)) < 1e-12 and abs(P.reach(1.0) - np.hypot(3.0, 3.0)) < 1e-12
    assert abs(IM.PixelPSF(10, 10, 0.8, parts=parts).reach() - 2.4) < 1e-12


def test_window_of_a_window_is_the_cropped_kernel_also_with_a_highpass():
    P = IM.PixelPSF(30, 34, 1.1, parts=((0.7, 1.0), (0.3, 2.5)), hp=3.0)
    c = np.array([[17.2, 14.6]])
    full = P.phi(c)[:, 0].reshape(30, 34)
    Wd = P.window(8, 21, 10, 26)
    W2 = Wd.window(2, 9, 3, 12)
    assert np.allclose(Wd.phi(c)[:, 0].reshape(13, 16), full[8:21, 10:26], atol=1e-13, rtol=0)
    assert np.allclose(W2.phi(c)[:, 0].reshape(7, 9), full[10:17, 13:22], atol=1e-13, rtol=0)



def test_the_gap_jump_stops_the_path_a_few_entries_after_the_record_ratio():
    P = psf(32)
    ct = np.array([[8.0, 8.0], [22.0, 8.0], [8.0, 22.0], [22.0, 22.0]])
    y = P.synth(ct, np.array([900.0, 850.0, 60.0, 5.0]))
    R0, _ = run(P, y)
    assert R0['taus'].size == 4                                                      # without a stop the noise-free path enters all four
    R1, _ = run(P, y, jump=5.0, extra=0)
    assert R1['taus'].size == 3 and R1['taus'][2] / R1['taus'][1] > 5.0              # the record ratio (14) is at the third entrant: stop right after it
    R2, _ = run(P, y, jump=5.0, extra=1)
    assert R2['taus'].size == 4
    R3, _ = run(P, y, jump=100.0, extra=0)
    assert R3['taus'].size == 4                                                      # no ratio reaches the jump


def test_a_neighbour_within_the_fit_radius_is_refitted_with_the_entrant():
    P = psf(30)
    ct = np.array([[10.0, 12.0], [13.0, 12.0]])
    y = P.synth(ct, np.array([800.0, 600.0]))
    R, _ = run(P, y, fit_radius=4.0, tol=1e-10)
    last = R['hist'][-1]
    o = np.argsort(last['c'][:, 0])
    assert last['w'].size == 2 and np.abs(last['w'][o] - [800.0, 600.0]).max() < 1e-2 and np.abs(last['c'][o, 0] - [10.0, 13.0]).max() < 1e-4
    assert abs(R['hist'][1]['w'][0] - R['hist'][2]['w'][0]) > 10.0 if len(R['hist']) > 2 else True        # the first unit's amplitude moved when its neighbour entered


def test_the_local_window_is_wide_enough_to_agree_with_the_global_fit():
    P = psf(36)
    ct = np.array([[14.0, 16.0], [19.0, 18.0]])
    y = P.synth(ct, np.array([900.0, 700.0]))
    Rl, _ = run(P, y, local=True, tol=1e-11)
    Rg, _ = run(P, y, local=False, tol=1e-11)
    wl, wg = Rl['hist'][-1]['w'], Rg['hist'][-1]['w']
    assert wl.size == wg.size == 2 and np.abs(np.sort(wl) - np.sort(wg)).max() / 900.0 < 1e-5


def test_the_window_reach_follows_the_largest_width_the_kernels_may_slide_to():
    P = psf(44)
    y = IM.PixelPSF(44, 44, 2.6).synth(np.array([[22.0, 22.0]]), np.array([2000.0]))      # a broad emitter, wider than the PSF the lattice kernels start from
    L = IM.Lattice(P, h=0.5, margin=1.0)
    R = IM.sliding_path(P, y, 1.0, L, max_entries=3, sig_bounds=(1.0, 3.2), local=True, tol=1e-11, rtrunc=5.5)
    last = R['hist'][-1]
    assert last['w'].size == 1 and abs(last['w'][0] - 2000.0) / 2000.0 < 1e-4 and abs(last['c'][0, 2] - 2.6) < 1e-3     # a broad emitter is fitted exactly inside its local window


def test_discrepancy_rule_degrees_of_freedom_are_n_minus_three_per_unit():
    hist = [dict(w=np.zeros(k), chi2=c) for k, c in enumerate([5.0, 3.0, 2.0, 1.20, 1.0])]
    # block of 3 units: chi2 x N / (N - 9) = 1.3187 > 1 + 2 sqrt(2/91) = 1.2962 (rejected with N - 3k; it would pass with N - 2k); block of 4 units passes
    assert IM.stop_discrepancy(hist, 100) == 4
    assert IM.stop_discrepancy(hist[:4], 100) == 3                                   # the last block is returned when none passes


def test_block_fit_clips_a_start_outside_its_box():
    P = psf(20)
    y = P.synth(np.array([[7.0, 7.0]]), np.array([500.0]))
    c, w, info = IM.block_fit(P, y, 1.0, np.array([[5.0, 5.0]]), np.array([1.0]), np.array([[6.0, 6.0]]), np.array([[8.0, 8.0]]), tol=1e-11)
    assert np.abs(c - [7.0, 7.0]).max() < 1e-5 and abs(w[0] - 500.0) < 1e-3


def _local_vs_global(P, y, sw, **kw):
    """The largest relative difference of the amplitudes after two entries between the local fits (the window around the cluster) and the global ones."""
    L = IM.Lattice(P, h=0.5, margin=1.0)
    a = IM.sliding_path(P, y, sw, L, local=True, tol=1e-11, max_entries=2, **kw)['hist'][-1]
    b = IM.sliding_path(P, y, sw, L, local=False, tol=1e-11, max_entries=2, **kw)['hist'][-1]
    assert a['w'].size == b['w'].size
    ia, ib = np.lexsort(a['c'].T), np.lexsort(b['c'].T)
    return float(np.abs(a['w'][ia] - b['w'][ib]).max() / np.abs(b['w']).max())


@pytest.mark.parametrize('seed', [0, 1])
def test_the_local_window_reaches_far_enough_for_the_tails_in_noise(seed):
    # with noise the window is not harmless: the pixels beyond it would still pull the fit; at the window of 5.5 widths local and global agree to 1e-8, at half that only to 1e-6
    P = psf(40)
    y = P.synth(np.array([[14.0, 16.0], [19.0, 18.0]]), np.array([900.0, 700.0])) + np.random.default_rng(seed).standard_normal((40, 40))
    assert _local_vs_global(P, y, 1.0) < 1e-8


@pytest.mark.parametrize('seed', [0, 1])
def test_the_local_window_follows_the_widest_kernel_in_noise(seed):
    # the window reaches rtrunc x the widest kernel the fit may slide to (3.2 here), not the PSF's width: a broad emitter in noise is fitted as the global fit does it
    P = psf(44)
    y = IM.PixelPSF(44, 44, 2.6).synth(np.array([[22.0, 22.0]]), np.array([2000.0])) + np.random.default_rng(seed).standard_normal((44, 44))
    assert _local_vs_global(P, y, 1.0, sig_bounds=(1.0, 3.2)) < 1e-3


def test_a_neighbour_between_half_and_the_whole_fit_radius_is_refitted_with_the_entrant():
    P = psf(36)
    y = P.synth(np.array([[10.0, 12.0], [14.5, 12.0]]), np.array([1000.0, 150.0]))
    R, _ = run(P, y, fit_radius=6.0, tol=1e-10)            # the neighbour is 4.5 px away: inside 6, outside 3
    H = R['hist']
    assert len(H) == 3 and np.abs(np.sort(H[-1]['w']) - [150.0, 1000.0]).max() < 1e-3 and H[-1]['chi2'] < 1e-12
    R3, _ = run(P, y, fit_radius=3.0, tol=1e-10)           # a radius below the separation leaves the first unit at the amplitude it had before the entrant came (inflated by the tail)
    assert R3['hist'][2]['w'][0] > 1004.0 and len(R3['hist']) > 3


def test_refine_keeps_the_centres_inside_its_box_and_reweights_from_the_fitted_model():
    P = psf(24)
    y = P.synth(np.array([[13.0, 12.0]]), np.array([600.0]))
    c0, w0 = np.array([[11.5, 12.0]]), np.array([300.0])
    c1, w1, _ = IM.refine(P, y, 1.0, c0, w0, box=1.0, rounds=1, tol=1e-11)
    assert abs(c1[0, 0] - 12.5) < 1e-6                                                      # stops at the edge of the box (a round moves a centre by at most `box`): the emitter is 1.5 px away
    c2, w2, _ = IM.refine(P, y, 1.0, c0, w0, box=2.0, rounds=1, tol=1e-11)
    assert abs(c2[0, 0] - 13.0) < 1e-5 and abs(w2[0] - 600.0) < 1e-2
    sw_fn = lambda m: 1.0 / np.sqrt(1.0 + np.maximum(m, 0.0))
    c3, w3, sw3 = IM.refine(P, y, 1.0, c2, w2, sw_fn=sw_fn, rounds=2, tol=1e-11)
    assert np.allclose(sw3, sw_fn(P.synth(c3, w3))) and sw3.min() < 0.2 and not np.allclose(sw3, 1.0)      # the returned weights are those of the final model


def test_fit_frame_time_rule_keeps_every_entrant_before_the_stopping_time_whatever_its_significance():
    P = IM.PixelPSF(30, 30, 1.3)
    L = IM.Lattice(P, h=0.5, margin=1.0)
    rng = np.random.default_rng(4)
    y = rng.poisson(np.maximum(P.synth(np.array([[10.3, 11.2], [20.4, 18.1], [14.0, 22.5]]), np.array([140.0, 70.0, 280.0])) + 20.0, 0)).astype(float)
    ot = IM.fit_frame(y, P, L, 20.0, 1.0, 0.0, stop=('time', 3.0), max_entries=8, r_ex=2.0)
    osn = IM.fit_frame(y, P, L, 20.0, 1.0, 0.0, stop=('snr', 3.0), max_entries=8, r_ex=2.0)
    tau_stop = 1.0 / (3.0 * IM.noise_scale(L, (1.0 / np.sqrt(20.0 * np.ones((30, 30)))) ** 2))
    assert ot['taus'][-1] <= tau_stop and ot['index'] == len(ot['hist']) - 1               # the path ran to the stopping time and its last block is the answer
    assert len(ot['w']) >= len(osn['w']) + 2                                               # the significance rule cuts earlier: the fourth entrant (significance 2.2) ends it, the time rule goes on
