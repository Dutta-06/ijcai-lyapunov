"""
Theorem J in the real closed loop: periodic orbits near the critical points of the regularised occupation loss
=============================================================================================================

For two instances and k_e ∈ {10, 20, 40, 80, 160}:
  1. find a critical point Θ* of J_m(Θ) = m L̄_m(Θ) + (σ₀/2)Σ‖ϑ_i − ϑ⁰‖²_{G⁻¹} (L-BFGS + Newton) and its Hessian;
  2. find the T₀-periodic orbit of the closed loop (σ = σ₀/k_e, no dither) by Newton shooting on the Poincaré map from
     (x_m(0), Θ*), with the monodromy matrix as a by-product;
  3. compare the period mean of the orbit with Θ* − J⁻¹(b_frozen + b_osc) (Proposition H.7 + Remark H.7a) and with the
     frozen-measure part alone;
  4. compare the parameter Floquet exponents with the eigenvalues of (1/k_e)G_m^{1/2}∇²J_m(Θ*)G_m^{1/2}.

Instance A: n = 2, m = 3, the plant and elliptic reference of the orientation experiment, σ₀ = 0.5: the minimiser is
a collapsed configuration (three identical kernels), a non-degenerate critical point with permutation symmetry.
Instance B: n = 1, m = 2, a two-bump plant on x_m = 2 sin t, σ₀ = 0.05: the minimiser has two distinct kernels.

Output: results/periodic_orbit/table.md (and raw.npz; set PO_REUSE=1 to regenerate the table from it).
"""
import os
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.core import langevin as L
from src.core import periodic as PO

OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'periodic_orbit')
PER = 2 * np.pi
KES = (10.0, 20.0, 40.0, 80.0, 160.0)


def instance_A():
    net = L.SaturatedLiftedNetwork(2, a_max=3.0, c_max=2.0, s_min=0.3, s_max=1.5)
    f = lambda X: np.stack([np.sin(X[:, 1]) + 0.3 * X[:, 0], np.cos(X[:, 0]) - 0.2 * X[:, 1]], 1)

    def J(X):
        out = np.zeros((len(X), 2, 2))
        out[:, 0, 0] = 0.3; out[:, 0, 1] = np.cos(X[:, 1]); out[:, 1, 0] = -np.sin(X[:, 0]); out[:, 1, 1] = -0.2
        return out
    xr = lambda t: np.array([1.2 * np.cos(t), 0.8 * np.sin(t)]); xd = lambda t: np.array([-1.2 * np.sin(t), 0.8 * np.cos(t)])
    rng = np.random.default_rng(11)
    m = 3
    Th0 = net.unconstrained(rng.normal(0, 0.4, (m, 2)), rng.uniform(-1, 1, (m, 2)), rng.uniform(0.5, 1.0, m))
    return dict(name='A (n = 2, m = 3, collapsed minimiser)', net=net, f=f, jac=J, xr=xr, xd=xd, G=np.array([4.0, 4.0, 1.0, 1.0, 1.0]),
                sigma0=0.5, Th0=Th0, m=m)


def instance_B():
    net = L.SaturatedLiftedNetwork(1, a_max=3.0, c_max=2.5, s_min=0.2, s_max=1.5)
    f = lambda X: 1.6 * np.exp(-(X - 1.0) ** 2 / 0.4) - 1.2 * np.exp(-(X + 1.0) ** 2 / 0.4) + 0.1 * X

    def J(X):
        d = (1.6 * np.exp(-(X - 1.0) ** 2 / 0.4) * (-2 * (X - 1.0) / 0.4) - 1.2 * np.exp(-(X + 1.0) ** 2 / 0.4) * (-2 * (X + 1.0) / 0.4) + 0.1)
        return d[:, :, None]
    xr = lambda t: np.array([2.0 * np.sin(t)]); xd = lambda t: np.array([2.0 * np.cos(t)])
    m = 2
    Th0 = net.unconstrained(np.array([[1.0], [-1.0]]), np.array([[0.8], [-0.8]]), np.array([0.6, 0.6]))
    return dict(name='B (n = 1, m = 2, two distinct kernels)', net=net, f=f, jac=J, xr=xr, xd=xd, G=np.array([2.0, 1.0, 1.0]),
                sigma0=0.05, Th0=Th0, m=m)


def _trace_jacobian(cl, t, z, m, h=1e-6):
    """tr ∂F/∂z of the closed-loop vector field at (t, z) by central differences."""
    tr = 0.0
    for k in range(len(z)):
        e = np.zeros_like(z); e[k] = h
        tr += (cl.rhs(t, z + e, m)[k] - cl.rhs(t, z - e, m)[k]) / (2 * h)
    return tr


def liouville_fast_exponent(cl, zp, DP, m, n_state, N=400):
    """
    Liouville's identity Σ_all log|λ| = ∫₀^{T₀} tr ∂F/∂z dt along the periodic orbit gives the product of the n_state
    error-channel multipliers, which the finite-difference monodromy resolves only as noise, from the resolved parameter
    multipliers.  Returns the mean error-channel Floquet exponent −(1/(n_state T₀))·Σ_fast log|λ| (positive = stable).
    """
    from scipy.integrate import solve_ivp
    ts = np.linspace(0, cl.T0, N, endpoint=False)
    sol = solve_ivp(lambda t, z: cl.rhs(t, z, m), (0.0, cl.T0), zp, method='DOP853', rtol=1e-11, atol=1e-13, t_eval=ts)
    int_tr = float(np.mean([_trace_jacobian(cl, t, sol.y[:, k], m) for k, t in enumerate(ts)]) * cl.T0)
    logs = np.sort(np.log(np.abs(np.linalg.eigvals(DP))))          # ascending: the n_state most negative are the error channel
    return -(int_tr - float(np.sum(logs[n_state:]))) / (n_state * cl.T0)


def run(args):
    i, ke = args
    inst = [instance_A, instance_B][i]()                       # rebuilt in the worker (closures do not pickle)
    net, m = inst['net'], inst['m']
    cl = PO.ClosedLoop(net, inst['f'], inst['xr'], inst['xd'], PER, ke, inst['G'], inst['sigma0'])
    orb = cl.orbit(N=2000)
    Th_star, H = cl.find_critical_point(inst['Th0'], orb)
    z0 = np.concatenate([inst['xr'](0.0), Th_star.ravel()])
    zp, DP, res = cl.shoot(z0, m)
    Thbar, Th_t, E_t, _ = cl.period_mean(zp, m)
    full, frozen = cl.offset_prediction(Th_star, orb, inst['jac'])
    fast, slow = PO.floquet_exponents(DP, PER, net.n)
    pred = PO.predicted_parameter_exponents(H, inst['G'], m, ke)
    fast_L = liouville_fast_exponent(cl, zp, DP, m, net.n)
    a_s, c_s, s_s = net.physical(Th_star)
    return dict(ke=ke, Th_star=Th_star, H=H, zp=zp, res=res, Thbar=Thbar, full=full, frozen=frozen, slow=slow, pred=pred, fast_L=fast_L,
                phys=(a_s, c_s, s_s),
                emax=float(np.abs(E_t).max()), spread=float(np.max(np.abs(Th_star - Th_star.mean(0)))),
                hess_min=float(np.linalg.eigvalsh(H)[0]), Jm=cl.J_m(Th_star, orb),
                osc_amp=float(np.max(np.abs(Th_t - Thbar[None]))))


def table(inst, rows):
    G, m = inst['G'], inst['m']
    w = 1 / np.tile(G, m)
    nrm = lambda v: float(np.sqrt(np.sum(w * np.ravel(v) ** 2)))
    cos = lambda a, b: float(np.sum(w * np.ravel(a) * np.ravel(b)) / (nrm(a) * nrm(b)))
    L = [f"### Instance {inst['name']}", "",
         f"σ₀ = {inst['sigma0']:g}, G = diag{tuple(float(g) for g in G)}; Θ* has J_m = {rows[0]['Jm']:.6f}, λ_min(∇²J_m) = {rows[0]['hess_min']:.4f}, "
         f"kernel spread {rows[0]['spread']:.2e} (0 means a collapsed configuration); its kernels have amplitudes "
         + ", ".join("(" + ", ".join(f"{v:.3f}" for v in a) + ")" for a in rows[0]['phys'][0]) + ", centres "
         + ", ".join("(" + ", ".join(f"{v:.3f}" for v in c) + ")" for c in rows[0]['phys'][1]) + " and widths "
         + ", ".join(f"{v:.3f}" for v in rows[0]['phys'][2]) + ".", "",
         "| k_e | shooting residual | max‖e‖ on the orbit | ‖Θ̄_p − Θ*‖ | full prediction: size ratio, cosine, vector error | frozen only: size ratio, cosine, vector error | parameter Floquet exponents against (1/k_e)λ(G_m^{1/2}∇²J_mG_m^{1/2}): max rel. error, median rel. error | k_e × max rel. error | mean error-channel exponent from Liouville's identity, minus k_e |",
         "|---|---|---|---|---|---|---|---|---|"]
    verr = []
    for r in rows:
        off = r['Thbar'] - r['Th_star']
        ve = nrm(off - r['full']) / nrm(r['full']); vf = nrm(off - r['frozen']) / nrm(r['frozen'])
        rel = np.abs(r['slow'] - r['pred']) / r['pred']
        verr.append(ve)
        L.append(f"| {r['ke']:g} | {r['res']:.1e} | {r['emax']:.2e} | {nrm(off):.3e} | {nrm(off) / nrm(r['full']):.3f}, {cos(off, r['full']):.4f}, {ve:.3f} | "
                 f"{nrm(off) / nrm(r['frozen']):.3f}, {cos(off, r['frozen']):.4f}, {vf:.3f} | {rel.max():.3e}, {np.median(rel):.3e} | {r['ke'] * rel.max():.3f} | {r['fast_L']:.3f}, {r['fast_L'] - r['ke']:+.3f} |")
    kes = np.array([r['ke'] for r in rows]); offs = np.array([nrm(r['Thbar'] - r['Th_star']) for r in rows])
    slope = np.polyfit(np.log(kes), np.log(offs), 1)[0]
    L += ["", f"Fitted exponent of ‖Θ̄_p − Θ*‖ against k_e: {slope:.2f} (Theorem J: −1). Vector error of the full prediction: "
          + ", ".join(f"{v:.3f}" for v in verr) + " (Theorem J(c): O(1/k_e))."]
    return L


def main():
    import multiprocessing as mp
    os.makedirs(OUT, exist_ok=True)
    insts = [instance_A(), instance_B()]
    raw = os.path.join(OUT, 'raw.npz')
    results = {}
    if os.path.exists(raw) and os.environ.get('PO_REUSE', '0') == '1':
        z = np.load(raw, allow_pickle=True)
        results = z['results'].item()
    else:
        jobs = [(i, ke) for i in range(len(insts)) for ke in KES]
        with mp.Pool(min(len(jobs), int(os.environ.get("PO_WORKERS", max(1, os.cpu_count()))))) as pool:
            rows = pool.map(run, jobs)
        for (i, ke), r in zip(jobs, rows):
            results[(i, ke)] = r
        np.savez(raw, results=np.array(results, dtype=object))
    lines = ["# Theorem J: periodic orbits of the real closed loop near the critical points of J_m (`exp_periodic_orbit.py`)", "",
             "The period mean Θ̄_p of the exact periodic orbit (Newton shooting on the Poincaré map, residual < 1e-10) is compared with "
             "the critical point Θ* of J_m and with the first-order prediction Θ* − J⁻¹(b_frozen + b_osc) of Proposition H.7 and "
             "Remark H.7a. The parameter Floquet exponents of the orbit (from the monodromy matrix) are compared with the eigenvalues "
             "of the averaged linearisation. Distances and cosines are in the labelled G⁻¹ metric. The error-channel multipliers are "
             "below the resolution of the finite-difference monodromy (exponents of order k_e); their product is recovered from Liouville's "
             "identity Σ log|λ| = ∫tr ∂F/∂z dt along the orbit, and the mean error-channel exponent is reported against k_e (Theorem J(b): "
             "k_e − L with L = O(1)).", ""]
    for i, inst in enumerate(insts):
        lines += table(inst, [results[(i, ke)] for ke in KES]) + [""]
    lines.append("Numbers are produced by `experiments/exp_periodic_orbit.py`; nothing is edited by hand.")
    with open(os.path.join(OUT, 'table.md'), 'w', encoding='utf-8') as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == '__main__':
    main()
