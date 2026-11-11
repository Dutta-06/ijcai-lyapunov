"""
A batch of independent dithered runs of `mass_transfer.ReplicatorCohesionLoop` (Propositions D‴.36–D‴.39), integrated together and
checkpointed per run (referee round 143: the runs of `exp_gauge_fixman.py` Part J and `exp_gauge_plateau.py` Part G).

`BatchDither(loop, **per_run)` evaluates the loop's drift for B configurations at once (arrays (B, m, 3) and (B, m)), with the formulas of
`ReplicatorCohesionLoop.drift`: the fit's gradient through the saturations, the cohesion force (compact kernel), a box's dead zone, and the
replicator's fitness with the Fixman correction in its masses form or its plug-in form, the denominator floored at `fixman_floor`·W_i. The
mass rate γ, the gain κ_F, the correction's τ and the floor may be given per run (arrays (B,)). The entropic leak, the tracked
contributions and the per-class correction are not implemented (asserted). `step` is one step of `ReplicatorCohesionLoop.run_dithered`'s
scheme: the positions by Euler–Maruyama with the noise √(2τ dt G), log β by Euler, renormalised to Σβ = m.

`run_batch` integrates a list of runs. Each run keeps its own generator and draws its noise in blocks, `rng.standard_normal((block, m,
3))`, which is the stream of `block` successive draws of shape (m, 3) — the draws of `run_dithered` — so a run's noise does not depend on
the batch it is in. The runs' states (positions, log-masses, step count, the generator's state, an exit step, a record accumulator) are
saved after every block to a checkpoint file and the integration resumes from it; `wall` stops it at a block boundary. A run ends at its
own number of steps, or (with `inside`) at the first step after which it has left a domain.

The drift agrees with `ReplicatorCohesionLoop.drift` to rounding, and a short run reproduces `run_dithered` (tests/test_gauge_batch.py).
"""
import os
import pickle
import time

import numpy as np

from src.core import cohesion_leak as CO


class BatchDither:
    def __init__(self, loop, gamma=None, fixman=None, tau_F=None, floor=None):
        assert not loop.eps and not loop.fixman_TA and loop.fixman_classes is None and loop.balanced
        assert loop.kind == 'compact'
        self.loop = loop
        base = loop.base
        self.net = base.net
        self.X, self.F, self.N = base.X, base.F, base.N
        self.G = np.asarray(loop.G, float)
        self.R, self.sigma0 = loop.R, loop.sigma0
        col = lambda v, d: np.asarray(d if v is None else v, float).reshape(-1, 1) if np.ndim(v) else float(d if v is None else v)
        self.gamma, self.fixman = col(gamma, loop.gamma), col(fixman, loop.fixman)
        self.tau_F, self.floor = col(tau_F, loop.fixman_tau), col(floor, loop.fixman_floor)
        self.A_hat = None if loop.fixman_A is None else np.asarray(loop.fixman_A, float)
        self.box, self.box_sigma = loop.box, loop.box_sigma

    def take(self, idx):
        """The batch restricted to the runs idx (per-run parameters sliced)."""
        out = object.__new__(BatchDither)
        out.__dict__.update(self.__dict__)
        for k in ('gamma', 'fixman', 'tau_F', 'floor'):
            v = getattr(self, k)
            if np.ndim(v):
                setattr(out, k, v[idx])
        return out

    def features(self, Th):
        """a, c, s (B, m) and diff, phi (B, m, N) (the transpose of AveragedSlowLoop.features' layout)."""
        net = self.net
        a = net.a_max * np.tanh(Th[:, :, 0] / net.a_max)
        c = net.c_max * np.tanh(Th[:, :, 1] / net.c_max)
        s = net.s_min + (net.s_max - net.s_min) / (1.0 + np.exp(-np.clip(Th[:, :, 2], -500.0, 500.0)))
        diff = self.X[None, None, :] - c[:, :, None]
        phi = np.exp(-diff ** 2 / (2 * s[:, :, None] ** 2))
        return a, c, s, diff, phi

    def drift(self, Th, beta):
        """(dΘ (B, m, 3), dβ (B, m)) of the loop at each of the B configurations."""
        net = self.net
        B, m, _ = Th.shape
        a, c, s, diff, phi = self.features(Th)
        r = self.F[None, :] - np.matmul((beta * a)[:, None, :], phi)[:, 0, :] / m             # (B, N)
        rc = r[:, :, None]
        pd = phi * diff
        pr = np.matmul(phi, rc)[:, :, 0]
        pdr = np.matmul(pd, rc)[:, :, 0]
        pd2r = np.matmul(pd * diff, rc)[:, :, 0]
        da = 1.0 - (a / net.a_max) ** 2
        dc = 1.0 - (c / net.c_max) ** 2
        sk = 1.0 / (1.0 + np.exp(-np.clip(Th[:, :, 2], -500.0, 500.0))); ds = (net.s_max - net.s_min) * sk * (1 - sk)
        g = np.empty_like(Th)
        g[:, :, 0] = pr * da
        g[:, :, 1] = a * dc / s ** 2 * pdr
        g[:, :, 2] = a * ds / s ** 3 * pd2r
        g /= self.N
        D = Th[:, :, None, :] - Th[:, None, :, :]                                            # (B, m, m, 3)
        d2 = np.sum(D ** 2 / self.G[None, None, None, :], axis=3)
        coh = -(self.sigma0 / m) * np.einsum('bij,bijk->bik', CO.psi(d2 / self.R ** 2, 'compact'), D)
        dTh = beta[:, :, None] * self.G[None, None, :] * g + coh
        if self.box is not None:
            dTh = dTh - self.box_sigma * (Th - np.clip(Th, -self.box[None, None, :], self.box[None, None, :]))
        if not np.any(self.gamma):
            return dTh, np.zeros_like(beta)
        u = a * pr / self.N
        if np.any(self.fixman):
            u = u + self.correction(a, beta)
        lam = np.sum(beta * u, axis=1, keepdims=True) / beta.sum(axis=1, keepdims=True)
        return dTh, self.gamma * beta * (u - lam)

    def correction(self, a, beta):
        """κ_F(τ/m)(1 + q)/max(W(1 − q), floor·W), q = (|Â|/(a_max W))² (the masses form) or (a/a_max)² (the plug-in form)."""
        m = beta.shape[1]
        W = beta / m
        amax = self.net.a_max
        q = (a / amax) ** 2 if self.A_hat is None else (np.abs(self.A_hat)[None, :] / amax / W) ** 2
        return self.fixman * (self.tau_F / m) * (1 + q) / np.maximum(W * (1 - q), self.floor * W)

    def step(self, Th, lb, dt, sq, z):
        """One step of run_dithered: Θ + dt·dΘ + sq·z (sq (3,) or (B, 3)), log β + dt·dβ/β renormalised to Σβ = m."""
        b = np.exp(lb)
        dT, db = self.drift(Th, b)
        sq = np.asarray(sq, float)
        Th = Th + dt * dT + (sq[None, None, :] if sq.ndim == 1 else sq[:, None, :]) * z
        if np.any(self.gamma):
            m = lb.shape[1]
            lb = lb + dt * db / b
            lb = lb - np.log(np.exp(lb).sum(axis=1, keepdims=True) / m)
        return Th, lb


def _dump(path, obj):
    tmp = f'{path}.tmp{os.getpid()}'
    with open(tmp, 'wb') as f:
        pickle.dump(obj, f)
    os.replace(tmp, path)


def load_states(path):
    if path is None or not os.path.exists(path):
        return {}
    with open(path, 'rb') as f:
        return pickle.load(f)


def run_batch(bd, runs, dt, block, ckpt=None, wall=np.inf, inside=None, rec_every=None, record=None, log=None):
    """Integrates the runs together. `bd` is a BatchDither whose per-run parameters follow the order of `runs`; each run is a dict with
    'key' (hashable), 'seed', 'n' (its number of steps, a multiple of `block`), 'sq' (3,) (√(2τ dt G)), 'Th0' (m, 3) and 'lb0' (m,), and
    optionally 'acc0' (its record accumulator). With `inside` (a function of Θ (B, m, 3) returning (B,) booleans) a run ends at the first
    step after which it is outside; with `record` (record(acc, k, Θ, log β) for one run, returning the updated accumulator) every
    `rec_every` steps (a divisor of `block`) and at its last step. Returns {key: state}, state = dict(Th, lb, k, rng, done, exit, acc);
    exit is the step count at the exit or None. The states are saved to `ckpt` after every block; the integration stops at a block
    boundary once `wall` seconds have passed, after at least one block (the states are then incomplete: some have done = False)."""
    t0 = time.time()
    saved = load_states(ckpt)
    states = {}
    for r in runs:
        st = saved.get(r['key'])
        if st is None:
            st = dict(Th=np.array(r['Th0'], float), lb=np.array(r['lb0'], float), k=0, rng=np.random.default_rng(r['seed']).bit_generator.state,
                      done=False, exit=None, acc=r.get('acc0'))
        assert r['n'] % block == 0 and (st['exit'] is not None or st['k'] % block == 0)
        if st['exit'] is None and st['k'] < r['n']:
            st['done'] = False                                                # a cap raised since the checkpoint
        states[r['key']] = st
    act = [i for i, r in enumerate(runs) if not states[r['key']]['done']]
    if not act:
        return states
    keys = [runs[i]['key'] for i in act]
    gens = []
    for kk in keys:
        g = np.random.default_rng(); g.bit_generator.state = states[kk]['rng']; gens.append(g)
    Th = np.array([states[kk]['Th'] for kk in keys]); lb = np.array([states[kk]['lb'] for kk in keys])
    k = np.array([states[kk]['k'] for kk in keys]); n = np.array([runs[i]['n'] for i in act])
    sq = np.array([runs[i]['sq'] for i in act]); acc = [states[kk]['acc'] for kk in keys]
    b = bd.take(np.array(act))
    m = Th.shape[1]
    nb = 0

    def sync(sel=None):
        for j, kk in enumerate(keys):
            if sel is not None and j not in sel:
                continue
            states[kk].update(Th=Th[j].copy(), lb=lb[j].copy(), k=int(k[j]), acc=acc[j])

    def drop(keep):
        nonlocal Th, lb, k, n, sq, acc, b, keys, gens, act, Z
        keep = np.asarray(keep, bool)
        ix = np.flatnonzero(keep)
        Th, lb, k, n, sq, Z = Th[ix], lb[ix], k[ix], n[ix], sq[ix], Z[ix]
        acc = [acc[j] for j in ix]; keys = [keys[j] for j in ix]; gens = [gens[j] for j in ix]; act = [act[j] for j in ix]
        b = bd.take(np.array(act, int)) if act else None

    Z = None
    while keys:
        off = int(k[0] % block)
        if off == 0:                                                          # a block boundary: checkpoint, then draw the next block
            sync()
            for j, kk in enumerate(keys):
                states[kk]['rng'] = gens[j].bit_generator.state
            if ckpt is not None:
                _dump(ckpt, states)
            if nb and time.time() - t0 > wall:                                # at least one block per call
                if log:
                    log(f"stopped at the wall ({time.time() - t0:.0f} s): {len(keys)} runs open")
                return states
            if nb and log:
                log(f"[{time.time() - t0:.0f}s] block {nb}: {len(keys)} runs open, step {int(k.min())}–{int(k.max())}")
            nb += 1
            Z = np.array([g.standard_normal((block, m, 3)) for g in gens])
        if record is not None and off % rec_every == 0:
            for j in range(len(keys)):
                acc[j] = record(acc[j], int(k[j]), Th[j], lb[j])
        Th, lb = b.step(Th, lb, dt, sq, Z[:, off])
        k += 1
        ended = k >= n
        if inside is not None:
            out = ~inside(Th)
            for j in np.flatnonzero(out):
                states[keys[j]].update(exit=int(k[j]))
            ended = ended | out
        if ended.any():
            for j in np.flatnonzero(ended):
                if record is not None and states[keys[j]]['exit'] is None:
                    acc[j] = record(acc[j], int(k[j]), Th[j], lb[j])
                if states[keys[j]]['exit'] is None:                           # at its cap, a block boundary: the generator has drawn
                    states[keys[j]]['rng'] = gens[j].bit_generator.state      # exactly k steps, so a raised cap continues the run
                states[keys[j]]['done'] = True
            sync(set(np.flatnonzero(ended).tolist()))
            drop(~ended)
            if not keys and ckpt is not None:
                _dump(ckpt, states)
    return states
