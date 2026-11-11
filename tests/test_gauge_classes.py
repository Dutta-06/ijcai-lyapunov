"""
Proposition D‴.42 (PARADIGM.md §5d): several particles per class. Checked here: at an exact fit with each class's particles together,
log det ∇²_ΘJ minus its one-particle-per-class value at the same class masses is a constant, whatever the split within the classes; the
class masses' drift and diffusion in Katzenberger's limit are 3/m times the one-particle-per-class ones in the slow-replicator limit, and
the fractions within a class do not diffuse; the correction per class (`fixman_classes`) is the one-particle-per-class correction at the
class masses times 3/m, given to every particle of the class, and its tracked contributions are summed per class; off M the masses form
reads the given contributions (the plug-in form the current ones), each class its own τ_k, the floor at fixman_floor·W_k near a class's
saturation; `limit_drift` leaves the eigenbasis for its Schur method at n = 3's degenerate equal split, the class law 1/3 there, the two
methods agreeing at a generic point.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'experiments')))
import exp_cohesion_leak as E                                            # noqa: E402
import exp_gauge_classes as C                                            # noqa: E402
from src.core import mass_transfer as MT                                 # noqa: E402


def test_the_determinant_does_not_see_the_split():
    net, base = E.setup()
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2)
    Wc = np.array([0.30, 0.35, 0.35])
    l1 = C.log_det(loop, E.exact_rep(net, Wc), 3 * Wc)
    d = [C.log_det(loop, *C.exact_m(net, Wc, s, 2)) - l1 for s in ((0.5, 0.5, 0.5), (0.1, 0.8, 0.35))]
    Wc2 = np.array([0.40, 0.25, 0.35])
    d.append(C.log_det(loop, *C.exact_m(net, Wc2, (0.3, 0.6, 0.5), 2)) - C.log_det(loop, E.exact_rep(net, Wc2), 3 * Wc2))
    assert np.ptp(d) < 1e-5


def test_the_class_masses_follow_the_law_at_m_and_the_fractions_do_not_diffuse():
    r = C.b_job((2, 1e-5, (0.30, 0.35), (0.3, 0.6, 0.8)))
    assert np.allclose(r['cls'] / r['cls3'], 0.5, rtol=2e-4) and np.allclose(r['Dcls'] / r['Dcls3'], 0.5, rtol=1e-6)
    assert np.abs(r['frac']).max() < 1e-3 * np.abs(r['cls']).max() and r['Dfrac'] < 1e-12 * r['Dcls'].max()


def test_the_correction_per_class():
    net, base = E.setup()
    W = np.array([0.30, 0.35, 0.35])
    loop3 = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2, fixman=2.0, fixman_tau=1e-4, fixman_A=E.A_BUMP.copy())
    c3 = loop3.fixman_correction(E.exact_rep(net, W), 3 * W)
    lab = np.repeat(np.arange(3), 2)
    loop6 = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2, fixman=2.0, fixman_tau=1e-4, fixman_A=E.A_BUMP.copy(), fixman_classes=lab)
    Th, beta = C.exact_m(net, W, (0.2, 0.6, 0.5), 2)
    c6 = loop6.fixman_correction(Th, beta)
    assert np.allclose(c6, 0.5 * c3[lab], rtol=1e-12)
    # the plug-in form per class reads the class's current contribution: on M the same
    plug = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2, fixman=2.0, fixman_tau=1e-4, fixman_classes=lab)
    assert np.allclose(plug.fixman_correction(Th, beta), c6, rtol=1e-12)
    # the tracked contributions are summed per class: on M without dither they stay at the bumps' amplitudes
    tr = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2, fixman=2.0, fixman_tau=0.0, fixman_A=E.A_BUMP + 0.01, fixman_TA=5.0,
                                   fixman_classes=lab)
    tr.run_dithered(Th, beta, 2.5, 0.25, 0.0, np.random.default_rng(0))
    assert np.allclose(tr.fixman_A - E.A_BUMP, 0.01 * (1 - 0.25 / 5.0) ** 10, rtol=1e-6)


def test_the_correction_per_class_off_the_manifold_and_at_the_floor():
    net, base = E.setup()
    lab = np.repeat(np.arange(3), 2)
    Th, beta = C.exact_m(net, np.array([0.30, 0.35, 0.35]), (0.2, 0.6, 0.5), 2)
    Th = Th + 0.05 * np.random.default_rng(3).standard_normal(Th.shape)  # off M: the current contributions are not Â
    A_hat, tau_k = E.A_BUMP.copy(), np.array([1e-4, 2e-4, 3e-4])
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2, fixman=2.0, fixman_tau=tau_k, fixman_A=A_hat, fixman_classes=lab)
    Wc = np.bincount(lab, beta, 3) / 6
    q = (np.abs(A_hat) / net.a_max / Wc) ** 2
    hand = 2.0 * (tau_k / 6) * (1 + q) / np.maximum(Wc * (1 - q), 0.01 * Wc)
    assert np.allclose(loop.fixman_correction(Th, beta), hand[lab], rtol=1e-13)
    plug = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2, fixman=2.0, fixman_tau=tau_k, fixman_classes=lab)
    assert np.abs(plug.fixman_correction(Th, beta) - hand[lab]).max() > 1e-3 * np.abs(hand).max()
    Wm = np.abs(A_hat) / net.a_max                                        # within δ_cap of class 1's saturation: the floor binds
    Wc2 = np.array([1.0005 * Wm[0], 0.5, 0.5 - 1.0005 * Wm[0]])
    beta2 = C.split_beta(Wc2, (0.3, 0.5, 0.6), 2)
    q2 = (Wm / Wc2) ** 2
    assert Wc2[0] * (1 - q2[0]) < 0.01 * Wc2[0]
    floor = 2.0 * (tau_k[0] / 6) * (1 + q2[0]) / (0.01 * Wc2[0])
    assert np.allclose(loop.fixman_correction(Th, beta2)[lab == 0], floor, rtol=1e-13)


def test_limit_drift_at_a_degenerate_equal_split():
    import exp_gauge_drift as X
    from src.core import gauge_drift as GDR
    r = C.b_job((3, 1e-2, (1 / 3, 1 / 3), (1 / 3,) * 6))                 # the eigenbasis is singular here (cond ~ 10¹⁸)
    assert np.allclose(r['cls'] / r['cls3'], 1 / 3, rtol=1e-4) and np.allclose(r['Dcls'] / r['Dcls3'], 1 / 3, rtol=1e-6)
    assert np.abs(r['frac']).max() < 1e-5 * np.abs(r['cls']).max()
    F, chart, coords, S = X.instance(1e-2)
    w = np.array([0.40, 0.25])
    d1 = GDR.limit_drift(F, chart, coords, w, S, h=X.H_ID, h2=X.H2_ID)[0]
    d2 = GDR.limit_drift(F, chart, coords, w, S, h=X.H_ID, h2=X.H2_ID, method='schur')[0]
    assert np.allclose(d1, d2, rtol=1e-5)
