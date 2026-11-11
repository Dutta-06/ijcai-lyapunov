"""
Proposition H.7 in the real closed loop: the learned kernels depend on the ORIENTATION of the reference orbit
===========================================================================================================

The averaged system (𝔄) depends on the reference only through its occupation measure ρ, which is the same for an
orbit and its reversal.  Proposition H.7 says the exactly averaged closed loop has an orientation-odd correction
(1/k_e²)G∇U_odd (filter lag; for n ≥ 2 it contains the orbit-area circulation kernel), and along the real loop the
second-order averaging term b_osc of the intra-period parameter oscillation (Remark H.7a) is orientation-odd too.
Hence the equilibria of the real closed loop for the two orientations differ, by

    Θ₊ − Θ₋ ≈ −2 J⁻¹ (b_odd + b_osc),     b_odd,i = (1/k_e²) G ∇_ϑ U_odd(ϑ_i),   b_osc = (1/k_e²)⟨D_ΘH·u⟩,

where J is the Jacobian of the leading-order averaged drift at the equilibrium (linear response).  In the
scaling σ = σ₀/k_e of Conjecture H.6 this difference is O(1/k_e).

Setting: n = 2, m = 6 saturated kernels, f(x) = (sin x₂ + 0.3x₁, cos x₁ − 0.2x₂), x_m = (1.2 cos t, ±0.8 sin t),
G = diag(4, 4, 1, 1, 1), ϑ⁰ = 0, σ = σ₀/k_e (σ₀ = 2 by default), no dither (deterministic), the same initial
parameters for both orientations.  The closed loop is simulated with the project code (langevin.simulate).

Outputs: results/orbit_orientation/table.md
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
from src.core import averaging as A
from src.core.plant_models import PlantModel

OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'orbit_orientation')
PERIOD = 2 * np.pi
STEPS_PER_PERIOD = 3000
DT = PERIOD / STEPS_PER_PERIOD
REC = 100                                                   # 30 snapshots per period
M, N = 6, 2
GAINS = np.array([4.0, 4.0, 1.0, 1.0, 1.0])
SIGMA0 = float(os.environ.get('OO_SIGMA0', 2.0))
KES = (10.0, 20.0, 40.0)
SLOW_PERIODS = float(os.environ.get('OO_SLOW', 16.0))      # horizon = SLOW_PERIODS · k_e periods


class EllipsePlant(PlantModel):
    def __init__(self, orient):
        super().__init__(n=2, name='ellipse')
        self.o = float(orient)

    def f(self, x, t):
        return np.array([np.sin(x[1]) + 0.3 * x[0], np.cos(x[0]) - 0.2 * x[1]])

    def x_ref(self, t):
        return np.array([1.2 * np.cos(t), self.o * 0.8 * np.sin(t)])

    def x_ref_dot(self, t):
        return np.array([-1.2 * np.sin(t), self.o * 0.8 * np.cos(t)])


def f_vec(X):
    return np.stack([np.sin(X[:, 1]) + 0.3 * X[:, 0], np.cos(X[:, 0]) - 0.2 * X[:, 1]], 1)


def jac_vec(X):
    J = np.zeros((len(X), 2, 2))
    J[:, 0, 0] = 0.3; J[:, 0, 1] = np.cos(X[:, 1]); J[:, 1, 0] = -np.sin(X[:, 0]); J[:, 1, 1] = -0.2
    return J


def net():
    return L.SaturatedLiftedNetwork(N, a_max=3.0, c_max=2.0, s_min=0.3, s_max=1.5)


def initial():
    rng = np.random.default_rng(11)
    a = rng.normal(0, 0.4, size=(M, N)); c = rng.uniform(-1.0, 1.0, size=(M, N)); s = rng.uniform(0.5, 1.0, size=M)
    return net().unconstrained(a, c, s)


def run(args):
    ke, orient = args
    nt = net()
    ctrl = L.LangevinAdaptiveController(nt, ke=ke, G_diag=GAINS, sigma=SIGMA0 / ke, theta0=np.zeros(5), tau=0.0)
    plant = EllipsePlant(orient)
    n_per = int(round(SLOW_PERIODS * ke))
    res = L.simulate(plant, ctrl, plant.x_ref(0.0), initial(), T=n_per * PERIOD, dt=DT,
                     rng=np.random.default_rng(0), record_every=REC)
    Th = res['Th']
    k = STEPS_PER_PERIOD // REC
    last = Th[-k:].mean(0)
    earlier = Th[-11 * k:-10 * k].mean(0)
    return dict(ke=ke, orient=orient, last=last, earlier=earlier)


def lab_dist(P, Q):
    return float(np.sqrt(np.mean(np.sum((P - Q) ** 2 / GAINS[None, :], axis=1))))


def kernel_spread(Th):
    """Largest coordinate deviation of a kernel from the mean kernel (0 for a collapsed configuration)."""
    return float(np.max(np.abs(Th - Th.mean(0))))


def predicted_difference(ke, Th_eq):
    """
    Linear response at the (+)-orientation equilibrium.  Returns (frozen, full):
      frozen = −2 J⁻¹ b_odd                 (frozen-measure correction of Proposition H.7 only)
      full   = −2 J⁻¹ (b_odd + b_osc)       (plus the second-order averaging term of the parameter oscillation)
    """
    nt = net()
    sig = SIGMA0 / ke
    xr = lambda t: np.array([1.2 * np.cos(t), 0.8 * np.sin(t)])
    xd = lambda t: np.array([-1.2 * np.sin(t), 0.8 * np.cos(t)])
    orb = A.Orbit(xr, xd, PERIOD, N=2000, sub=1)
    m, p = Th_eq.shape
    h = 1e-5

    def grad_at_particles(Th, fun):
        E = np.vstack([Th + h * np.eye(p)[k][None, :] for k in range(p)] + [Th - h * np.eye(p)[k][None, :] for k in range(p)])
        v = fun(Th, E).reshape(2, p, m)
        return (v[0] - v[1]).T / (2 * h)                     # (m, p)

    def lead_drift(Th):
        fm = A.FrozenMeasure(nt, Th)
        g = grad_at_particles(Th, lambda T_, E: A.leading_potential(f_vec, fm, orb, ke, E))
        return (GAINS[None, :] * g - sig * Th).ravel()

    fm = A.FrozenMeasure(nt, Th_eq)
    g_odd = grad_at_particles(Th_eq, lambda T_, E: A.expansion_terms(f_vec, jac_vec, fm, orb, E)[1])
    b_odd = (GAINS[None, :] * g_odd / ke ** 2).ravel()
    b_osc = (A.oscillation_correction(nt, f_vec, orb, Th_eq, GAINS) / ke ** 2).ravel()
    d = m * p
    J = np.zeros((d, d))
    eps = 1e-5
    base = Th_eq.ravel()
    for j in range(d):
        e = np.zeros(d); e[j] = eps
        J[:, j] = (lead_drift((base + e).reshape(m, p)) - lead_drift((base - e).reshape(m, p))) / (2 * eps)
    frozen = (-2.0 * np.linalg.solve(J, b_odd)).reshape(m, p)
    full = (-2.0 * np.linalg.solve(J, b_odd + b_osc)).reshape(m, p)
    return frozen, full


def main():
    import multiprocessing as mp
    os.makedirs(OUT, exist_ok=True)
    raw = os.path.join(OUT, f'raw_sigma{SIGMA0:g}_slow{SLOW_PERIODS:g}.npz')
    if os.path.exists(raw) and os.environ.get('OO_REUSE', '1') == '1':
        z = np.load(raw)
        R = {(float(k), float(o)): dict(last=z[f'last_{k:g}_{o:+d}'], earlier=z[f'earlier_{k:g}_{o:+d}']) for k in KES for o in (+1, -1)}
    else:
        jobs = [(ke, o) for ke in KES for o in (+1, -1)]
        with mp.Pool(len(jobs)) as pool:
            rows = pool.map(run, jobs)
        R = {(r['ke'], float(r['orient'])): r for r in rows}
        np.savez(raw, **{f'last_{k:g}_{o:+d}': R[(k, float(o))]['last'] for k in KES for o in (+1, -1)},
                 **{f'earlier_{k:g}_{o:+d}': R[(k, float(o))]['earlier'] for k in KES for o in (+1, -1)})
    w = 1 / np.tile(GAINS, M)
    cos = lambda a, b: float(np.sum(w * a.ravel() * b.ravel()) / np.sqrt(np.sum(w * a.ravel() ** 2) * np.sum(w * b.ravel() ** 2)))
    lines = ["# Orientation of the reference orbit changes the learned kernels (Proposition H.7, real closed loop)", "",
             f"n = 2, m = {M}, G = diag{tuple(float(g) for g in GAINS)}, σ = {SIGMA0:g}/k_e, no dither, x_m = (1.2 cos t, ±0.8 sin t) "
             "(the same occupation measure for both orientations), same initial parameters, horizon "
             f"{SLOW_PERIODS:g}·k_e periods.  Distances are labelled G⁻¹ distances between period-mean parameters.", "",
             "Two linear-response predictions of Θ₊ − Θ₋ at the (+) equilibrium: *frozen* uses only the frozen-measure "
             "correction of Proposition H.7 (−2J⁻¹b_odd); *full* adds the second-order averaging term of the intra-period "
             "parameter oscillation (−2J⁻¹(b_odd + b_osc)).", "",
             "| k_e | measured ‖Θ₊ − Θ₋‖ | residual drift, last 10 periods (+, −) | frozen: size ratio measured/predicted | frozen: cosine | "
             "full: size ratio measured/predicted | full: cosine | full: vector error ‖meas − pred‖/‖pred‖ |", "|---|---|---|---|---|---|---|---|"]
    meas_n = []
    for ke in KES:
        Pp, Pm = R[(ke, 1.0)]['last'], R[(ke, -1.0)]['last']
        meas = Pp - Pm
        frozen, full = predicted_difference(ke, Pp)
        nm = lab_dist(Pp, Pm)
        meas_n.append(nm)
        rp = lab_dist(R[(ke, 1.0)]['last'], R[(ke, 1.0)]['earlier']); rm = lab_dist(R[(ke, -1.0)]['last'], R[(ke, -1.0)]['earlier'])
        nf = lab_dist(frozen, 0 * frozen); nu = lab_dist(full, 0 * full)
        verr = lab_dist(meas, full) / nu
        lines.append(f"| {ke:g} | {nm:.3e} | {rp:.1e}, {rm:.1e} | {nm / nf:.3f} | {cos(meas, frozen):.4f} | {nm / nu:.3f} | {cos(meas, full):.4f} | {verr:.3f} |")
    slope = float(np.polyfit(np.log(np.array(KES)), np.log(np.array(meas_n)), 1)[0])
    end_spread = max(kernel_spread(R[(ke, o)]['last']) for ke in KES for o in (1.0, -1.0))
    lines += ["", f"Fitted exponent of the measured difference against k_e: {slope:.2f} (prediction: −1).",
              f"Kernel spread (largest coordinate deviation of a kernel from the mean kernel): {kernel_spread(initial()):.2f} "
              f"at the start, at most {end_spread:.1e} at the end of every run. The {M} kernels converge to a common "
              "configuration, so the learned network is effectively a single kernel and the test probes the orientation "
              f"effect at a collapsed equilibrium of the m = {M} system.",
              "Reading: the occupation measure is identical for the two orientations, so (𝔄) cannot tell them apart. The "
              "difference is the orientation-odd part of the first-order correction of the real averaged loop; compare the "
              "frozen and full columns to see how much of it the frozen-measure expansion alone accounts for."]
    with open(os.path.join(OUT, f'table_sigma{SIGMA0:g}_slow{SLOW_PERIODS:g}.md'), 'w', encoding='utf-8') as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
