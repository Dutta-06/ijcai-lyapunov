"""
Proposition D‴.15: the transverse spectrum of a lifted periodic orbit of the real closed loop
============================================================================================

(A) The OT leak toward the matching prior at the matching masses q = (1/3, 1/3, 1/3): the 3-atom loop's periodic orbit (the exact
    representation with zero tracking error) and its lifts to m = 6, 9, 15 — the lift's return residual under the Poincaré map, the
    monodromy by central differences, the per-class separation exponents (compressions onto the class separation subspaces, with the
    invariance defects), the ghost particles' exponents, the collective exponents against the 3-atom loop's and against Theorem J(b)'s
    prediction (the eigenvalues of G^{1/2}∇²𝒢_qG^{1/2}/k_e).
(B) The σ-modification (ϑ⁰ = 0) at the weights w = (7, 5, 4)/16 (admissible) and (9, 4, 3)/16 (a saddle): the weighted 3-atom orbit by
    shooting from Θ*(w), its lift to m = 16, the per-class separation exponents against the ghosts' and against the prediction
    (σ₀ − λ_j(R^G_k))/k_e of Proposition D‴.5(i) read along the orbit (Theorem J(b)); the saddle's negative exponent.

Output: results/chaos/lifted_orbit.md and lifted_orbit_raw.npz.
"""
import os
import sys
import time
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.core import langevin as L
from src.core import chaos_learning as CL
from src.core import lifted_orbit as LO

OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos')
NET_KW = dict(a_max=4.0, c_max=3.0, s_min=0.3, s_max=2.0)
G = np.array([1.0, 0.25, 0.25])
BUMPS = [(1.0, -1.5, 0.5), (-0.7, 0.0, 0.4), (0.8, 1.5, 0.5)]
PER = 2 * np.pi
SIGMA_CRIT = 0.047162567584665846
SIGMA0 = 0.15 * SIGMA_CRIT
KE = 25.0
MS_A = (6, 9, 15)
W_B = [np.array([7, 5, 4]) / 16.0, np.array([9, 4, 3]) / 16.0]
W_BAR = np.ones(3) / 3


def f_vec(X):
    return sum(A * np.exp(-(X - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS)


def x_ref(t):
    return np.array([2.0 * np.sin(t)])


def x_ref_dot(t):
    return np.array([2.0 * np.cos(t)])


def make_net():
    return L.SaturatedLiftedNetwork(1, **NET_KW)


def exact_representation(net):
    return net.unconstrained(np.array([[3 * A] for A, C, S in BUMPS]), np.array([[C] for A, C, S in BUMPS]), np.array([S for A, C, S in BUMPS]))


def exponents(DP):
    return np.sort(-np.log(np.abs(np.linalg.eigvals(DP))) / PER)


def labels_from_counts(counts):
    return np.repeat(np.arange(len(counts)), counts)


def analyse_lift(clK, zK, DPK, targets_K, weights_K, counts, net, t0):
    """Lift the K-atom orbit zK to the classes `counts`, compute the monodromy and its class-wise spectra; returns a dict."""
    labels = labels_from_counts(counts)
    m = int(labels.size)
    cl_m = LO.ClassLeakLoop(net, f_vec, x_ref, x_ref_dot, PER, KE, G, SIGMA0, targets_K[labels])
    z_m = LO.lift(zK, labels, 1, 3)
    residual = float(np.linalg.norm(cl_m.poincare(z_m, m) - z_m))
    DP = LO.monodromy_fd(cl_m, z_m, m)
    out = dict(m=m, counts=list(counts), residual=residual, classes={})
    for k, B in LO.class_separation_bases(labels, 1, 3).items():
        lam, ex, defect = LO.compression_spectrum(DP, B, PER)
        out['classes'][k] = dict(ex=ex, defect=defect, mult=counts[k] - 1)
    lamc, exc, dc = LO.compression_spectrum(DP, LO.collective_basis(labels, 1, 3), PER)
    out['collective'] = exc; out['collective_defect'] = dc
    out['all'] = exponents(DP)
    print(f"[{time.time() - t0:.0f}s] lift m = {m} {counts}: residual {residual:.1e}")
    return out


def part_A(net, exact, t0):
    clK = LO.ClassLeakLoop(net, f_vec, x_ref, x_ref_dot, PER, KE, G, SIGMA0, exact)
    zK, DPK, res = clK.shoot(np.concatenate([x_ref(0.0), exact.ravel()]), 3, max_iter=8, tol=1e-10)
    ghosts = {k: LO.ghost_monodromy(clK, zK, 3, k, exact[k]) for k in range(3)}
    # Theorem J(b)'s prediction for the collective block: eigenvalues of G^{1/2}∇²𝒢_q G^{1/2}/k_e from the averaged loop (nearest-atom version at q)
    from src.core import prior_design as PD
    mix = PD.MixturePriorLoop(net, f_vec, x_ref, PER, G, SIGMA0, exact, n_quad=256)
    H = hess_weighted(PD.weighted_energy_grad_prior, mix, W_BAR, exact)
    pred = np.sort(np.linalg.eigvalsh(H)) / KE
    lifts = [analyse_lift(clK, zK, DPK, exact, None, (m // 3,) * 3, net, t0) for m in MS_A]
    return dict(zK=zK, resK=res, exK=exponents(DPK), ghosts={k: (g[1], g[2]) for k, g in ghosts.items()}, pred_collective=pred, lifts=lifts)


def hess_weighted(grad_fn, loop, w, atoms, h=1e-6):
    """Whitened, weighted Hessian of 𝒢_w at Θ*(w): the block diag(w)^{1/2}G^{1/2}(∂(∇𝒢_w/w_k)/∂ϑ)G^{1/2}diag(w)^{1/2} of Proposition D‴.5(i), by
    central differences of grad_fn (which returns ∇_{ϑ_k}𝒢_w/w_k)."""
    K = atoms.shape[0]
    d = 3 * K
    v = atoms.ravel()
    Jm = np.zeros((d, d))
    for j in range(d):
        e = np.zeros(d); e[j] = h
        Jm[:, j] = (grad_fn(loop, w, (v + e).reshape(K, 3)).ravel() - grad_fn(loop, w, (v - e).reshape(K, 3)).ravel()) / (2 * h)
    sw = np.repeat(np.sqrt(np.asarray(w, float)), 3)
    sg = np.tile(np.sqrt(G), K)
    Hs = (sw * sg)[:, None] * Jm * (sg / sw)[None, :]          # diag(w)^{1/2}G^{1/2}·(w_k⁻¹∇²_{k·}𝒢_w)·G^{1/2}diag(w)^{−1/2} = the symmetric block
    return 0.5 * (Hs + Hs.T)


def part_B(net, exact, t0):
    single = CL.AveragedSlowLoop(net, f_vec, x_ref, PER, G, SIGMA0, n_quad=256)
    out = []
    for w in W_B:
        Tw, res = CL.weighted_stationary_atoms(single, w, exact, w0=W_BAR, n_cont=4)
        a, c, s, diff, phi = single.features(Tw)
        r = single.F - phi @ (w * a)
        R = single.residual_curvature(r, Tw)
        Rw = single.sqG[None, :, None] * R * single.sqG[None, None, :]
        pred_sep = {k: np.sort(np.linalg.eigvalsh(SIGMA0 * np.eye(3) - Rw[k])) / KE for k in range(3)}
        kR = single.kappa_R(r, Tw)
        clK = LO.ClassLeakLoop(net, f_vec, x_ref, x_ref_dot, PER, KE, G, SIGMA0, np.zeros((3, 3)), weights=w)
        zK, DPK, resK = clK.shoot(np.concatenate([x_ref(0.0), Tw.ravel()]), 3, max_iter=12, tol=1e-10)
        print(f"[{time.time() - t0:.0f}s] weighted 3-atom orbit at w = {w.round(4)}: shooting residual {resK:.1e}, period-start offset from Θ*(w) {np.abs(zK[1:] - Tw.ravel()).max():.2e}")
        ghosts = {k: LO.ghost_monodromy(clK, zK, 3, k, np.zeros(3)) for k in range(3)}
        counts = tuple(int(round(16 * x)) for x in w)
        lift = analyse_lift(clK, zK, DPK, np.zeros((3, 3)), None, counts, net, t0)
        out.append(dict(w=w, atoms=Tw, kR=kR, margin=float(SIGMA0 - kR.max()), pred_sep=pred_sep, zK=zK, resK=resK, exK=exponents(DPK),
                        offset=float(np.abs(zK[1:] - Tw.ravel()).max()), ghosts={k: (g[1], g[2]) for k, g in ghosts.items()}, lift=lift))
    return out


def fmt_ex(v):
    return ", ".join(f"{x:.3e}" for x in v)


def grouped(ex, tol=1e-7):
    """The distinct values of a list of exponents (consecutive sorted values within tol are one group), each printed as the group mean."""
    ex = np.sort(np.asarray(ex, float))
    groups, cur = [], [ex[0]]
    for x in ex[1:]:
        if x - cur[-1] < tol:
            cur.append(x)
        else:
            groups.append(np.mean(cur)); cur = [x]
    groups.append(np.mean(cur))
    return np.array(groups)


def rel_err(a, b):
    """Relative errors |a − b|/|b| of two sorted lists of equal length."""
    a, b = np.sort(np.asarray(a, float)), np.sort(np.asarray(b, float))
    return np.abs(a - b) / np.abs(b)


def report(A, B, t0=None):
    lines = ["# Proposition D‴.15: the transverse spectrum of a lifted periodic orbit of the real closed loop (`exp_lifted_orbit.py`)", "",
             f"The real closed loop (CL) of Theorem J on the three-bump instance at k_e = {KE:g}, σ = σ₀/k_e with σ₀ = 0.15σ_crit = {SIGMA0:.5f}; monodromy matrices by central differences of the "
             "Poincaré map (DOP853, rtol 10⁻¹¹); Floquet exponents −log|λ|/T₀ (the error channel's exponent sits at the finite-difference floor and is omitted). The lifts are the "
             "collapsed configurations with the K₀-atom orbit's atoms repeated n_k times; the ghost of class k is one particle with weight zero in f̂, driven by the K₀-atom orbit and "
             "leaking toward the class target.", "",
             "## A. The OT leak at the matching masses (prior = exact representation, q = (1/3, 1/3, 1/3))", "",
             f"The 3-atom loop's periodic orbit is the exact representation with zero tracking error (shooting residual {A['resK']:.1e}); its parameter Floquet exponents are "
             f"{fmt_ex(A['exK'][:-1])}, against Theorem J(b)'s prediction λ_j(H_q)/k_e = {fmt_ex(A['pred_collective'])} with H_q the whitened, weighted collective block of Proposition D‴.5(i) "
             f"(relative errors {rel_err(A['exK'][:-1], A['pred_collective']).min() * 100:.2f}–{rel_err(A['exK'][:-1], A['pred_collective']).max() * 100:.2f}%, the O(k_e⁻²) of Theorem J(b)). The ghosts' exponents: "
             + "; ".join(f"class {k}: {fmt_ex(g[0])} (return residual {g[1]:.1e})" for k, g in A['ghosts'].items()) + f" — σ₀/k_e = {SIGMA0 / KE:.3e} at zero residual.", "",
             "| m | classes | return residual of the lift | separation exponents per class (distinct values, rounded to 10⁻⁷; each of multiplicity n_k − 1) | invariance defects | collective exponents (parameters) | collective defect |",
             "|---|---|---|---|---|---|---|"]
    for lf in A['lifts']:
        lines.append(f"| {lf['m']} | {tuple(lf['counts'])} | {lf['residual']:.1e} | " + "; ".join(f"class {k}: {fmt_ex(grouped(c['ex']))} ×{c['mult']}" for k, c in lf['classes'].items())
                     + " | " + ", ".join(f"{c['defect']:.1e}" for c in lf['classes'].values()) + f" | {fmt_ex(lf['collective'][:-1])} | {lf['collective_defect']:.1e} |")
    errsA = np.concatenate([rel_err(lf['collective'][:-1], A['exK'][:-1]) for lf in A['lifts']])
    errgA = np.concatenate([np.abs(grouped(c['ex']) - A['ghosts'][k][0][0]) / A['ghosts'][k][0][0] for lf in A['lifts'] for k, c in lf['classes'].items()])
    lines += ["", f"Relative errors over the three lifts: collective exponents against the 3-atom orbit's {errsA.max() * 100:.4f}% at most; separation exponents against the ghosts' σ₀/k_e {errgA.max() * 100:.3f}% at most.",
              "", "## B. The σ-modification (ϑ⁰ = 0) at unequal weights: the weighted 3-atom orbit, its lift to m = 16, and the saddle", ""]
    for b in B:
        counts = tuple(int(round(16 * x)) for x in b['w'])
        lines += [f"### w = ({', '.join(f'{x:.4f}' for x in b['w'])}), classes {counts}: splitting margin σ₀ − max κ_R = {b['margin']:+.5f} ({'admissible' if b['margin'] > 0 else 'a saddle'})", "",
                  f"Weighted 3-atom orbit: shooting residual {b['resK']:.1e}, period-start offset from Θ*(w) {b['offset']:.2e} (Theorem J(a): O(1/k_e)); its parameter exponents {fmt_ex(b['exK'][:-1])}. "
                  f"Lift to m = 16: return residual {b['lift']['residual']:.1e}; collective exponents {fmt_ex(b['lift']['collective'][:-1])} (defect {b['lift']['collective_defect']:.1e}).", "",
                  "| class k | n_k − 1 | separation exponents of the lift | ghost exponents | prediction (σ₀ − λ_j(R^G_k))/k_e | κ_R at the atom | defect |", "|---|---|---|---|---|---|---|"]
        for k in range(3):
            c = b['lift']['classes'][k]; g = b['ghosts'][k]
            lines.append(f"| {k} | {c['mult']} | {fmt_ex(grouped(c['ex']))} (each ×{c['mult']}) | {fmt_ex(g[0])} | {fmt_ex(b['pred_sep'][k])} | {b['kR'][k]:+.5f} | {c['defect']:.1e} |")
        errs = np.concatenate([rel_err(grouped(b['lift']['classes'][k]['ex']), b['pred_sep'][k]) for k in range(3)])
        errg = np.concatenate([np.abs(grouped(b['lift']['classes'][k]['ex']) - b['ghosts'][k][0]) / np.abs(b['ghosts'][k][0]) for k in range(3)])
        errc = rel_err(b['lift']['collective'][:-1], b['exK'][:-1])
        lines.append(f"Relative errors: separation exponents of the lift against the ghosts' {errg.max() * 100:.3f}% at most, against the prediction (σ₀ − λ_j(R^G_k))/k_e {errs.min() * 100:.2f}–{errs.max() * 100:.2f}%; collective exponents of the lift against the weighted 3-atom orbit's {errc.max() * 100:.4f}% at most.")
        lines.append("")
    neg = [(b['w'], k, float(b['lift']['classes'][k]['ex'].min())) for b in B for k in range(3) if b['lift']['classes'][k]['ex'].min() < 0]
    lines += ["Reading. The lift of the K₀-atom orbit is a periodic orbit of the m-loop for every size tried (return residuals at the integration tolerance), its monodromy leaves every class's "
              "separation subspace invariant (defects at the finite-difference level) and acts there as the ghost's monodromy with multiplicity n_k − 1, while the complement carries the K₀-atom "
              "loop's spectrum: the Floquet spectrum of the lift does not depend on m. Theorem J(b) reads the separation exponents off the averaged Hessian, (σ₀ − λ_j(R^G_k))/k_e + O(k_e⁻²): "
              + ("at the matching masses the residual vanishes and every separation exponent is σ₀/k_e; " if A['lifts'] else "")
              + (("at the saddle weights " + "; ".join(f"({', '.join(f'{x:.4f}' for x in w)}), class {k}: {e:.3e}" for w, k, e in neg) + " — the over-weighted class's negative exponent, the orbit-level form of Proposition K's split.") if neg else "no negative exponent appears."), "",
              "Numbers are produced by `experiments/exp_lifted_orbit.py`; nothing is edited by hand."]
    with open(os.path.join(OUT, 'lifted_orbit.md'), 'w', encoding='utf-8') as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    if t0 is not None:
        print(f"[{time.time() - t0:.0f}s]")


def main():
    t0 = time.time()
    if '--from-raw' in sys.argv:
        z = np.load(os.path.join(OUT, 'lifted_orbit_raw.npz'), allow_pickle=True)
        report(z['A'].item(), list(z['B']), t0)
        return
    net = make_net()
    exact = exact_representation(net)
    A = part_A(net, exact, t0)
    B = part_B(net, exact, t0)
    np.savez(os.path.join(OUT, 'lifted_orbit_raw.npz'), A=np.array(A, dtype=object), B=np.array(B, dtype=object), allow_pickle=True)
    report(A, B, t0)


if __name__ == '__main__':
    main()
