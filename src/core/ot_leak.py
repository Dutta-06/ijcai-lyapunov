"""
Proposition D‴.12 (PARADIGM.md §5d): the barycentric OT leak toward a prior measure — the mass constraint against the
weight fluctuation
====================================================================================================================

The nearest-atom leak of Proposition D‴.9 conserves the cluster masses of the network: with a locally constant assignment no
particle changes cell, so the cluster weights w^m of the initial cloud are carried forever, and the stationary floor of
Proposition D‴.5 — the (1/m − 1/M) law — is the memory of the initial condition that this conservation law keeps.  The
barycentric OT leak of Theorem C′ (§4b) toward ν = Σ_j q_jδ_{ϑ⁰_j} replaces the conservation by a constraint: the optimal plan γ
from μ_Θ to ν (row sums 1/m, column sums q_j) gives every prior atom its mass, and the leak points at the barycentre of the prior
mass assigned to each particle,

    dϑ_i/ds = G⟨∇_ϑ(r_Θ Φ̃(·;ϑ_i))⟩_ρ − σ₀(ϑ_i − T̄_i),        T̄_i := m Σ_j γ_{ij} ϑ⁰_j,                               (𝔄_m^{OT})

the G-gradient flow of J_m^{OT}(Θ) = mL̄(μ_Θ) + (σ₀m/2)·W²_{G⁻¹}(μ_Θ, ν) wherever the plan is unique (Danskin; Theorem C′(i)).
When mq_j ∈ ℕ for every j (m *commensurate* with ν) the vertex plans are integral — every particle is assigned to one atom and
atom j receives exactly mq_j particles — so (𝔄_m^{OT}) is (𝔄_m^ν) with a frozen assignment, and the replicated weighted
configuration Θ*_m(q) (mq_k particles at the weighted stationary atom ϑ*_k(q)) is a critical point for *every* commensurate m:
the stationary measure Σ_kq_kδ_{ϑ*_k}, the learned function and hence the stationary floors do not depend on m.  For other m a
vertex plan splits at most K₀ − 1 particles between atoms (a vertex of the transportation polytope is supported on a forest).

This module implements (𝔄_m^{OT}) as a subclass of `prior_design.MixturePriorLoop` (whose nearest-atom leak it keeps as
`self.nearest` for the comparison on the same clouds), the exact W²_{G⁻¹} between two empirical measures, the half-separation
check of the diagonal plan, and the coupled run of `experiments/exp_ot_leak.py`.
"""
import numpy as np
from scipy.optimize import linear_sum_assignment, linprog
from scipy import sparse

from .chaos_learning import _sig, cluster_weights
from .prior_design import MixturePriorLoop, weighted_stationary_atoms_prior


def _cost_matrix(A, B, G):
    """C_{ij} = ‖A_i − B_j‖²_{G⁻¹} without forming the (mA, mB, 3) difference array."""
    Ginv = 1.0 / np.asarray(G, float)
    AA = np.sum(A ** 2 * Ginv[None, :], axis=1)
    BB = np.sum(B ** 2 * Ginv[None, :], axis=1)
    C = AA[:, None] + BB[None, :] - 2.0 * (A * Ginv[None, :]) @ B.T
    return np.maximum(C, 0.0)


def transport_plan(C, a, b, tol=1e-12):
    """A vertex optimal plan of the transportation LP min Σγ_{ij}C_{ij}, row sums a, column sums b (Σa = Σb), by the HiGHS dual
    simplex (sparse constraints, one redundant equality dropped).  Returns γ ≥ 0 clipped at 0."""
    m, n = C.shape
    rows = sparse.kron(sparse.eye(m, format='csr'), np.ones((1, n)), format='csr')          # row sums
    cols = sparse.kron(np.ones((1, m)), sparse.eye(n, format='csr'), format='csr')          # column sums
    A_eq = sparse.vstack([rows, cols[:-1]], format='csr')
    b_eq = np.concatenate([a, b[:-1]])
    res = linprog(C.ravel(), A_eq=A_eq, b_eq=b_eq, bounds=(0, None), method='highs-ds')
    if not res.success:
        raise RuntimeError(res.message)
    gamma = res.x.reshape(m, n)
    gamma[np.abs(gamma) < tol] = 0.0
    return np.clip(gamma, 0.0, None)


def W2_empirical(A, B, G):
    """W²_{G⁻¹} between the uniform empirical measures of the rows of A (m rows) and B (M rows), exactly: when one size divides
    the other the smaller cloud is replicated (the measure is unchanged) and the assignment problem is solved (Birkhoff);
    otherwise the transportation LP with uniform masses."""
    A = np.asarray(A, float); B = np.asarray(B, float)
    mA, mB = A.shape[0], B.shape[0]
    if mA == mB or max(mA, mB) % min(mA, mB) == 0:
        if mA < mB:
            A = np.repeat(A, mB // mA, axis=0)
        elif mB < mA:
            B = np.repeat(B, mA // mB, axis=0)
        C = _cost_matrix(A, B, G)
        rows, cols = linear_sum_assignment(C)
        return float(C[rows, cols].mean())
    C = _cost_matrix(A, B, G)
    gamma = transport_plan(C, np.full(mA, 1.0 / mA), np.full(mB, 1.0 / mB))
    return float(np.sum(gamma * C))


class OTPriorLoop(MixturePriorLoop):
    """(𝔄_m^{OT}): the averaged slow loop with the barycentric OT leak toward the prior measure ν = Σ_j q_jδ_{prior_j}
    (prior: K₀ × 3 unconstrained coordinates; q: masses, default uniform).  `self.nearest` is the nearest-atom loop (𝔄_m^ν) of
    Proposition D‴.9 on the same instance."""

    def __init__(self, net, f_vec, x_ref, period, G, sigma0, prior, q=None, n_quad=256, plan_tol=1e-9):
        super().__init__(net, f_vec, x_ref, period, G, sigma0, prior, n_quad=n_quad)
        self.K0 = self.prior.shape[0]
        q = np.ones(self.K0) / self.K0 if q is None else np.asarray(q, float)
        assert q.shape == (self.K0,) and np.all(q > 0)
        self.q = q / q.sum()
        self.plan_tol = plan_tol
        self.nearest = MixturePriorLoop(net, f_vec, x_ref, period, G, sigma0, prior, n_quad=n_quad)
        d = np.sqrt(_cost_matrix(self.prior, self.prior, self.G))
        self.Delta_nu = float(d[~np.eye(self.K0, dtype=bool)].min()) if self.K0 > 1 else np.inf   # min inter-atom distance (G⁻¹)

    # ------------------------------------------------------------------ the plan
    def commensurate(self, m):
        """mq_j ∈ ℕ for every j."""
        counts = m * self.q
        return bool(np.allclose(counts, np.round(counts), atol=self.plan_tol))

    def cost_matrix(self, Th):
        """C_{ij} = ‖ϑ_i − ϑ⁰_j‖²_{G⁻¹}, (m, K₀)."""
        return _cost_matrix(np.asarray(Th, float), self.prior, self.G)

    def plan(self, Th):
        """A vertex optimal plan γ (m × K₀; row sums 1/m, column sums q_j) for the cost ‖·‖²_{G⁻¹}, and W²_{G⁻¹}(μ_Θ, ν) = Σγ_{ij}C_{ij}.
        Commensurate m: the assignment problem on the atoms replicated mq_j times (an integral plan, Hungarian algorithm; total
        unimodularity makes it optimal among all plans); otherwise the transportation LP (a vertex with at most K₀ − 1 split rows)."""
        Th = np.asarray(Th, float)
        m = Th.shape[0]
        C = self.cost_matrix(Th)
        if self.commensurate(m):
            n = np.round(m * self.q).astype(int)
            rep = np.repeat(np.arange(self.K0), n)
            rows, cols = linear_sum_assignment(C[:, rep])
            gamma = np.zeros((m, self.K0))
            gamma[rows, rep[cols]] = 1.0 / m
        else:
            gamma = transport_plan(C, np.full(m, 1.0 / m), self.q)
        return gamma, float(np.sum(gamma * C))

    def leak_target(self, Th, gamma=None):
        """T̄_i = mΣ_jγ_{ij}ϑ⁰_j, the barycentric projection of the plan, (m, 3)."""
        gamma = self.plan(Th)[0] if gamma is None else gamma
        return Th.shape[0] * (gamma @ self.prior)

    def assignment(self, Th, gamma=None):
        """The atom carrying the largest share of particle i's mass, (m,) (the assigned atom when the plan is integral)."""
        gamma = self.plan(Th)[0] if gamma is None else gamma
        return np.argmax(gamma, axis=1)

    def split_mask(self, Th, gamma=None):
        """Rows of the plan with two or more nonzero entries: the particles split between atoms (at most K₀ − 1 at a vertex)."""
        gamma = self.plan(Th)[0] if gamma is None else gamma
        return (gamma > 1e-6 / Th.shape[0]).sum(axis=1) >= 2

    def transit_mask(self, Th, gamma=None):
        """Unsplit particles whose assigned atom is not their nearest atom: the mass in transit between the clusters (a split row has no
        single assigned atom and is excluded — referee round 11)."""
        gamma = self.plan(Th)[0] if gamma is None else gamma
        return (self.assignment(Th, gamma) != self.nearest.assignment(Th)) & ~self.split_mask(Th, gamma)

    def W2(self, Th):
        """W²_{G⁻¹}(μ_Θ, ν)."""
        return self.plan(Th)[1]

    # ------------------------------------------------------------------ the objective and the flow
    def J(self, Th):
        """J_m^{OT}(Θ) = mL̄(μ_Θ) + (σ₀m/2)W²_{G⁻¹}(μ_Θ, ν)."""
        m = Th.shape[0]
        return m * self.L_bar(Th) + 0.5 * self.sigma0 * m * self.W2(Th)

    def drift(self, Th, T=None):
        """−G∇J_m^{OT} with the leak target T (default: the barycentric projection of the plan at Th)."""
        T = self.leak_target(Th) if T is None else T
        a, c, s, diff, phi = self.features(Th)
        m = Th.shape[0]
        r = self.F - phi @ a / m
        net = self.net
        da = 1.0 - (a / net.a_max) ** 2
        dc = 1.0 - (c / net.c_max) ** 2
        sk = _sig(Th[:, 2]); ds = (net.s_max - net.s_min) * sk * (1 - sk)
        g = np.empty_like(Th)
        g[:, 0] = (phi.T @ r) * da
        g[:, 1] = a * dc / s ** 2 * ((phi * diff).T @ r)
        g[:, 2] = a * ds / s ** 3 * ((phi * diff ** 2).T @ r)
        return self.G[None, :] * g / self.N - self.sigma0 * (Th - T)

    # ------------------------------------------------------------------ the diagonal plan and the replicated configuration
    def diagonal_plan_margin(self, atoms):
        """(ε, Δ_ν/2) with ε := max_k‖ϑ*_k − ϑ⁰_k‖_{G⁻¹}: when ε < Δ_ν/2 the diagonal plan (class k → atom k) is the unique optimal plan
        of every replicated configuration of `atoms` with masses q (c-cyclical monotonicity, Proposition D‴.12(iii))."""
        eps = float(np.sqrt(np.sum((np.asarray(atoms, float) - self.prior) ** 2 / self.G[None, :], axis=1)).max())
        return eps, 0.5 * self.Delta_nu

    def replicated(self, atoms, m):
        """Θ*_m: round(mq_k) particles at atoms[k] (m commensurate); returns (Θ (m, 3), class labels (m,))."""
        n = np.round(m * self.q).astype(int)
        assert n.sum() == m, "m is not commensurate with the masses q"
        return np.repeat(np.asarray(atoms, float), n, axis=0), np.repeat(np.arange(self.K0), n)

    def stationary_atoms(self, Th0, w0=None, n_cont=1, **kw):
        """The weighted stationary point Θ*(q) of 𝒢_q with the prior atoms (Proposition D‴.5 with the masses q): Newton on the
        nearest-atom loop (the diagonal plan makes the OT leak at the replicated configuration the nearest-atom leak toward ϑ⁰_k).
        Returns (Θ*(q), max|∇|, (ε, Δ_ν/2))."""
        atoms, res = weighted_stationary_atoms_prior(self.nearest, self.q, Th0, w0=w0, n_cont=n_cont, **kw)
        return atoms, res, self.diagonal_plan_margin(atoms)

    # ------------------------------------------------------------------ the coupled run
    def coupled_run_ot(self, Th_all, m_list, S_end, dt, rec_every=50, plan_every=5, leak='ot', label_atoms=None, w2_every=1,
                       record_final=True):
        """
        As AveragedSlowLoop.coupled_run (the m-systems are the first m rows of Th_all, the reference its M rows; midpoint rule in
        slow time), with the leak targets recomputed every `plan_every` steps (the vertex-selection solutions of §4b; the targets are
        frozen inside a step) and the records of Proposition D‴.12.  For every system, the reference included: the cluster weights
        by nearest atom of `label_atoms` (default: the prior atoms), the mass in transit (assigned ≠ nearest atom), the number of
        split particles, the fit error and W²_{G⁻¹}(μ_Θ, ν).  For every m-system against the reference: the synchronous D (mean
        square, G⁻¹ metric, the same particles), the function-level distance² ‖f̂_m − f̂_M‖²_{L²(ρ)}, and the exact W²_{G⁻¹}(μ_m, μ_M)
        (every `w2_every` records).  leak = 'nearest' runs (𝔄_m^ν) of Proposition D‴.9 with the same recorder, for the comparison on
        the same clouds.
        """
        assert leak in ('ot', 'nearest')
        M = Th_all.shape[0]
        m_list = [int(m) for m in m_list]
        sys_ = {m: Th_all[:m].copy() for m in m_list + [M]}
        atoms = self.prior if label_atoms is None else np.asarray(label_atoms, float)
        n_steps = int(round(S_end / dt))
        T = {}
        rec = dict(s=[], w={m: [] for m in sys_}, transit={m: [] for m in sys_}, n_split={m: [] for m in sys_}, fit={m: [] for m in sys_},
                   W2nu={m: [] for m in sys_}, D={m: [] for m in m_list}, hnorm2={m: [] for m in m_list}, W2={m: [] for m in m_list})
        for k in range(n_steps):
            if k % plan_every == 0:
                for m, Th in sys_.items():
                    T[m] = self.leak_target(Th) if leak == 'ot' else self.nearest.leak_target(Th)
            new = {}
            for m, Th in sys_.items():
                k1 = self.drift(Th, T[m])
                k2 = self.drift(Th + 0.5 * dt * k1, T[m])
                new[m] = Th + dt * k2
            sys_ = new
            if k % rec_every == 0 or k == n_steps - 1:
                rec['s'].append((k + 1) * dt)
                i_rec = len(rec['s']) - 1
                ref = sys_[M]
                fh_ref = self.f_hat(ref)
                for m, Th in sys_.items():
                    gamma, W2nu = self.plan(Th)
                    rec['w'][m].append(cluster_weights(Th, atoms, self.G)[1])
                    rec['transit'][m].append(float(np.mean(self.transit_mask(Th, gamma))) if leak == 'ot' else 0.0)
                    rec['n_split'][m].append(int(self.split_mask(Th, gamma).sum()) if leak == 'ot' else 0)
                    rec['fit'][m].append(self.fit_error(Th))
                    rec['W2nu'][m].append(W2nu)
                for m in m_list:
                    d = sys_[m] - ref[:m]
                    rec['D'][m].append(float(np.mean(np.sum(d ** 2 / self.G[None, :], axis=1))))
                    rec['hnorm2'][m].append(float(np.mean((self.f_hat(sys_[m]) - fh_ref) ** 2)))
                    rec['W2'][m].append(W2_empirical(sys_[m], ref, self.G) if i_rec % w2_every == 0 else float('nan'))
        rec['s'] = np.array(rec['s'])
        for key in ('w', 'transit', 'n_split', 'fit', 'W2nu', 'D', 'hnorm2', 'W2'):
            rec[key] = {m: np.array(v) for m, v in rec[key].items()}
        if record_final:
            rec['final'] = sys_
            if leak == 'ot':
                rec['final_plan'] = {m: self.plan(Th)[0] for m, Th in sys_.items()}
        return rec


# ---------------------------------------------------------------------- post-processing
def equalisation_time(s, w, q, m):
    """The first record time at which the cluster weights equal q to within half a particle (|w_k − q_k| < 1/(2m) for all k) and stay
    so until the end; nan if never."""
    ok = np.all(np.abs(np.asarray(w) - np.asarray(q)[None, :]) < 0.5 / m, axis=1)
    if not ok[-1]:
        return float('nan')
    idx = len(ok) - 1
    while idx > 0 and ok[idx - 1]:
        idx -= 1
    return float(s[idx])


def bridge_report(loop, Th, gamma):
    """For the split particles of a configuration: the plan row (m·γ_i, the barycentric weights), the distance to the nearest prior
    atom in units of Δ_ν and the distance to the barycentric target T̄_i (G⁻¹ metric)."""
    mask = loop.split_mask(Th, gamma)
    out = []
    m = Th.shape[0]
    T = loop.leak_target(Th, gamma)
    for i in np.flatnonzero(mask):
        d_atoms = np.sqrt(loop.cost_matrix(Th[i:i + 1])[0])
        out.append(dict(i=int(i), weights=m * gamma[i], dist_nearest_over_Delta=float(d_atoms.min() / loop.Delta_nu),
                        dist_target=float(np.sqrt(np.sum((Th[i] - T[i]) ** 2 / loop.G)))))
    return out


# ---------------------------------------------------------------------- dimension-agnostic leak targets (the real loop in any n)
def barycentric_targets(Th, prior, q, G):
    """T̄ = mγ·prior for a vertex optimal plan γ from the uniform empirical measure of the rows of Th (m, p) to ν = Σ_j q_jδ_{prior_j}
    (K₀, p), cost ‖·‖²_{G⁻¹}; the assignment route when mq_j ∈ ℕ, the transportation LP otherwise.  Returns (T̄ (m, p), γ (m, K₀))."""
    Th = np.asarray(Th, float); prior = np.atleast_2d(np.asarray(prior, float)); q = np.asarray(q, float) / np.sum(q)
    m, K0 = Th.shape[0], prior.shape[0]
    C = _cost_matrix(Th, prior, G)
    counts = m * q
    if np.allclose(counts, np.round(counts), atol=1e-9):
        n = np.round(counts).astype(int)
        rep = np.repeat(np.arange(K0), n)
        rows, cols = linear_sum_assignment(C[:, rep])
        gamma = np.zeros((m, K0)); gamma[rows, rep[cols]] = 1.0 / m
    else:
        gamma = transport_plan(C, np.full(m, 1.0 / m), q)
    return m * (gamma @ prior), gamma


def nearest_targets(Th, prior, G):
    """The nearest-atom targets prior[j(ϑ_i)] in the G⁻¹ metric, (m, p), and the labels."""
    C = _cost_matrix(np.asarray(Th, float), np.atleast_2d(np.asarray(prior, float)), G)
    lab = np.argmin(C, axis=1)
    return np.atleast_2d(prior)[lab], lab
