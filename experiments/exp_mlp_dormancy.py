"""
Do the paper's two practical messages hold outside RBF / KAN?  A one-hidden-layer tanh MLP student, small initialisation, readout of effective depth p in {1, 2, 3}.
==================================================================================================================================================================

Plan: results/ml/mlp_plan.md (pre-registration, committed before any test run).  Tests: tests/test_exp_mlp_dormancy.py.  numpy only (float32 forward / backward, float64 parameters).

Teacher  f*(x) = sum_{k<=K} b_k tanh(w_k . x),  x ~ N(0, I_d), |w_k| = 1 (well-separated: pairwise |cos| <= 0.5), b_k = +-U(0.5, 1.5).
Student  f(x)  = sum_{i<=H} a_i tanh(w_i . x), H = 128, no biases.  Readout parametrised at effective depth p: p = 1: a = u;  p = 2: a = u^2 - v^2;  p = 3: a = u^3 - v^3 (u = v = s at the
         start for p >= 2 so that a = 0; p = 1: u_i = s * N(0, 1)); inner weights w_i ~ N(0, s^2 / d) per entry (|w_i| about s); s^2 = eps.
Loss     0.5 * mean squared error on n = 4000 noisy training rows, full batch; early stopping on 1000 noisy validation rows (best-validation iterate kept).
Arms     gd (one rate), adam (default betas, eps_adam 1e-8), shared_v (Adam with ONE scalar second moment = EMA of mean g^2 over all parameters), sign (sign descent), adam_big (Adam with
         eps_adam = 0.1: GD-like where the gradient is small).
Outcomes number of active neurons after merging (strength |a_i| |w_i| above tau, then single-linkage components of the oriented directions sign(a_i) w_i / |w_i| with cos > c): count error
         |count - K|; test MSE against the CLEAN teacher on 20000 rows relative to the teacher variance; recovery of the teacher directions (best cosine per teacher); convergence flags.

    python experiments/exp_mlp_dormancy.py --stage dev|test --out DIR [--shard i/n] [--workers W] [--freeze results/ml/mlp_freeze.json]
    python experiments/exp_mlp_dormancy.py --merge --stage dev|test --out DIR [--raw OUT.npz]
    python experiments/exp_mlp_dormancy.py --from-raw [RAW.npz] [--md OUT.md]            (the report; byte-identical regeneration)
"""
import argparse
import itertools
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_OUT = os.environ.get('MLP_OUT', os.path.join(HERE, '..', 'results', 'ml'))

# ------------------------------------------------------------------------------------------------------------------------------------
# registries (fixed in results/ml/mlp_plan.md)
# ------------------------------------------------------------------------------------------------------------------------------------
H, N_TR, N_VAL, N_TE = 128, 4000, 1000, 20000
KS, DS, NOISES = (2, 4, 8), (10, 30), (0.05, 0.0)
EPS = (1e-2, 1e-4, 1e-6)
DEPTHS = (1, 2, 3)
ARMS = ('gd', 'adam', 'shared_v', 'sign', 'adam_big')
COMBOS = (('gd', 1), ('gd', 2), ('gd', 3), ('adam', 3), ('shared_v', 3), ('sign', 3), ('adam_big', 3))     # (arm, depth) cells that are run
EA_BIG = 0.1                                                                                           # eps_adam of the large-eps arm
COS_MAX = 0.5                                                                                          # teacher directions: pairwise |cos| <= COS_MAX
TMAX, PATIENCE, EVERY = 80000, 5000, 50                                                                # step budget, minimum patience (steps), validation period
ARM_FRAC, IMPROVE = 0.9, 1e-3                                                                          # armed once val < ARM_FRAC * val0; an improvement is a relative drop > IMPROVE
DIVERGE = 1e6                                                                                          # diverged if val > DIVERGE * val0 or not finite
# disjointness of development and test: teacher ids and data / initialisation seeds are drawn from separate ranges
DEV_TIDS, TEST_TIDS = (901,), tuple(range(1, 11))
NOISE0_TIDS, NEVER_TIDS = tuple(range(1, 4)), (1,)                                       # test design: noise 0 only for these teachers per cell; (gd, depth 3, eps 1e-6) only for these
STAGE_SEED = dict(dev=7_000_000, test=1_000_000)
LR_GRID = {('gd', 1): (0.1, 0.25, 0.6, 1.5), ('gd', 2): (0.05, 0.15, 0.4, 1.0), ('gd', 3): (0.05, 0.15, 0.4, 1.0), ('adam', 3): (3e-4, 1e-3, 3e-3, 1e-2),
           ('shared_v', 3): (3e-3, 1e-2, 3e-2, 1e-1), ('sign', 3): (1e-4, 3e-4, 1e-3, 3e-3), ('adam_big', 3): (3e-2, 1e-1, 3e-1, 1.0)}     # development grid per (arm, depth): 4 rates each
LR_EXTRA = {('sign', 3, 1e-2): (3e-5,), ('sign', 3, 1e-4): (3e-5,), ('sign', 3, 1e-6): (3e-5,), ('shared_v', 3, 1e-2): (1e-3,), ('shared_v', 3, 1e-4): (1e-3,), ('shared_v', 3, 1e-6): (0.3,),
            ('adam_big', 3, 1e-6): (3.0,)}                                                             # plan amendment 2: one more rate beyond the grid edge where the first development round chose an edge rate with a better score than its neighbour
DEV_NOISES = (0.05,)                                                                                   # development runs use the noisy level only
FIT_DEV = 0.5                                                                                          # an eps is feasible for depth-3 GD if at least this share of its development runs fits
# declared analysis settings
TAU = 0.05                                                                                             # activity threshold, strength |a| |w|
COS_MERGE = 0.95                                                                                       # merge cosine
TAUS = (0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5)                                                   # sweep of activity thresholds
COSES = (0.5, 0.7, 0.8, 0.9, 0.95, 0.98, 0.99, 1.01)                                                   # sweep of merge cosines (1.01: no merging)
REC_COS = 0.9                                                                                          # a teacher direction is recovered if some active neuron has oriented cos > REC_COS
FIT_REL = 0.2                                                                                          # a run "fits" if the relative clean test MSE < FIT_REL
N_BOOT, BOOT_SEED = 20000, 20261110
MARGIN_CE = 0.5                                                                                        # (R), (O): a count-error gain must be at least this many neurons
ABS_T = 0.01                                                                                           # (T): non-inferiority margin of the relative test MSE (difference, in teacher variances)
SWEEP_SHARE = 0.5                                                                                      # (O-sweep): share of sweep cells with the contrast above 0
HALF_H = 0.5                                                                                           # (R): depth 1 forms at least HALF_H * H neurons (before merging) at eps <= 1e-4


# ------------------------------------------------------------------------------------------------------------------------------------
# teacher and data
# ------------------------------------------------------------------------------------------------------------------------------------
def make_teacher(K, d, tid):
    """Teacher: K unit directions with pairwise |cos| <= COS_MAX (rejection) and amplitudes b_k = +-U(0.5, 1.5)."""
    rng = np.random.default_rng([11, int(K), int(d), int(tid)])
    while True:
        Wt = rng.standard_normal((K, d))
        Wt /= np.linalg.norm(Wt, axis=1, keepdims=True)
        C = np.abs(Wt @ Wt.T) - np.eye(K)
        if C.max() <= COS_MAX:
            break
    b = rng.uniform(0.5, 1.5, K) * rng.choice([-1.0, 1.0], K)
    return Wt, b


def teacher_fn(Wt, b, X):
    return np.tanh(X @ Wt.T) @ b


def make_data(stage, K, d, tid, noise):
    """Training (4000), validation (1000, noisy), clean test (20000) rows; the seed depends on the stage, the cell, the teacher and the noise level."""
    Wt, b = make_teacher(K, d, tid)
    rng = np.random.default_rng([STAGE_SEED[stage], int(K), int(d), int(tid), int(round(noise * 1000))])
    X = rng.standard_normal((N_TR, d))
    Xv = rng.standard_normal((N_VAL, d))
    Xt = rng.standard_normal((N_TE, d))
    y = teacher_fn(Wt, b, X) + noise * rng.standard_normal(N_TR)
    yv = teacher_fn(Wt, b, Xv) + noise * rng.standard_normal(N_VAL)
    yt = teacher_fn(Wt, b, Xt)
    return dict(Wt=Wt, b=b, X=X, y=y, Xv=Xv, yv=yv, Xt=Xt, yt=yt)


# ------------------------------------------------------------------------------------------------------------------------------------
# student
# ------------------------------------------------------------------------------------------------------------------------------------
def readout(P, depth):
    """Output weights a (float64) from the flat parameter vector P = [W (H x d), u (H), v (H) for depth >= 2]; products factorised (u - v)(...) against cancellation."""
    u = P[-2 * H:-H] if depth >= 2 else P[-H:]
    if depth == 1:
        return u
    v = P[-H:]
    if depth == 2:
        return (u - v) * (u + v)
    return (u - v) * (u * u + u * v + v * v)


def n_params(d, depth):
    return H * d + (H if depth == 1 else 2 * H)


def init_params(d, depth, eps, seed):
    """Inner weights w_i ~ N(0, s^2 / d) per entry (|w_i| about s, s = sqrt(eps)); u = v = s (depth >= 2) or u_i = s N(0, 1) (depth 1); the same seed gives the same inner weights
    (and, for depth 1, the same standard normal) for every depth and arm."""
    s = np.sqrt(eps)
    rng = np.random.default_rng([5, int(seed)])
    W = rng.standard_normal((H, d)) * s / np.sqrt(d)
    z = rng.standard_normal(H)
    if depth == 1:
        return np.concatenate([W.ravel(), s * z])
    return np.concatenate([W.ravel(), np.full(H, s), np.full(H, s)])


_BUF = {}


def loss_grad(P, depth, d, X, y, want_grad=True):
    """0.5 * mean squared error and the NEGATIVE gradient with respect to P (float32 matrix products, float64 parameters)."""
    W = P[:H * d].reshape(H, d)
    a = readout(P, depth)
    n = X.shape[0]
    buf = _BUF.get(n)
    if buf is None:
        buf = _BUF[n] = (np.empty((n, H), np.float32), np.empty((n, H), np.float32))
    Tn, B = buf
    np.matmul(X, W.astype(np.float32).T, out=Tn)
    np.tanh(Tn, out=Tn)
    a32 = a.astype(np.float32)
    r = Tn @ a32 - y
    L = 0.5 * float(np.mean(r.astype(np.float64) ** 2))
    if not want_grad:
        return L, None
    Ga = (Tn.T @ r).astype(np.float64) / n                                  # dL/da
    np.multiply(Tn, Tn, out=B)
    np.subtract(1.0, B, out=B)
    B *= a32
    B *= r[:, None]
    Gm = (B.T @ X).astype(np.float64) / n
    if depth == 1:
        g = np.concatenate([Gm.ravel(), Ga])
    else:
        u, v = P[-2 * H:-H], P[-H:]
        p = depth
        g = np.concatenate([Gm.ravel(), p * u ** (p - 1) * Ga, -p * v ** (p - 1) * Ga])
    return L, -g


def mse(P, depth, d, X, y):
    W = P[:H * d].reshape(H, d)
    f = np.tanh(X @ W.astype(np.float32).T) @ readout(P, depth).astype(np.float32)
    return float(np.mean((f.astype(np.float64) - y) ** 2))


class Opt:
    """One optimizer on the flat vector (callers pass the NEGATIVE gradient): gd p += lr g; sign p += lr sign(g); adam (bias-corrected, per-coordinate second moment);
    shared_v: the same with ONE scalar second moment v = EMA of mean(g^2); adam_big: adam with eps_adam = EA_BIG."""

    def __init__(self, arm, lr, beta1=0.9, beta2=0.999):
        assert arm in ARMS
        self.arm, self.lr, self.b1, self.b2 = arm, float(lr), beta1, beta2
        self.ea = EA_BIG if arm == 'adam_big' else 1e-8
        self.m = self.v = None
        self.t = 0

    def step(self, P, g):
        self.t += 1
        if self.arm == 'gd':
            return P + self.lr * g
        if self.arm == 'sign':
            return P + self.lr * np.sign(g)
        if self.m is None:
            self.m, self.v = np.zeros_like(P), np.zeros_like(P)
        self.m = self.b1 * self.m + (1 - self.b1) * g
        gg = np.mean(g ** 2) * np.ones_like(g) if self.arm == 'shared_v' else g ** 2
        self.v = self.b2 * self.v + (1 - self.b2) * gg
        mh, vh = self.m / (1 - self.b1 ** self.t), self.v / (1 - self.b2 ** self.t)
        return P + self.lr * mh / (np.sqrt(vh) + self.ea)


def fit(data, arm, depth, eps, lr, seed, tmax=TMAX, patience=PATIENCE):
    """Train from the small start.  Returns the best-validation iterate and the flags: `armed` (val fell below ARM_FRAC * val0 at some check), `converged` (armed and no relative
    improvement above IMPROVE for max(`patience`, last-improvement step) steps), `diverged`, `steps` (steps run), `best_step`.  A run that never arms runs to `tmax` (budget exhausted)."""
    X, y = data['X'].astype(np.float32), data['y'].astype(np.float32)
    Xv, yv = data['Xv'].astype(np.float32), data['yv']
    d = X.shape[1]
    P = init_params(d, depth, eps, seed)
    opt = Opt(arm, lr)
    val0 = mse(P, depth, d, Xv, yv)
    best, bestP, best_step, last_imp, armed, diverged, converged = val0, P.copy(), 0, 0, False, False, False
    t = 0
    while t < tmax:
        _, g = loss_grad(P, depth, d, X, y)
        P = opt.step(P, g)
        t += 1
        if t % EVERY == 0:
            val = mse(P, depth, d, Xv, yv)
            if not np.isfinite(val) or val > DIVERGE * val0:
                diverged = True
                break
            if val < best:                                                         # the iterate of lowest validation MSE is kept; only a relative drop above IMPROVE resets the patience
                if val < best * (1 - IMPROVE):
                    last_imp = t
                best, bestP, best_step = val, P.copy(), t
            if not armed and val < ARM_FRAC * val0:
                armed, last_imp = True, t
            if armed and t - last_imp >= max(patience, last_imp):                  # doubling rule: stop when the time since the last improvement exceeds the time it took to get there
                converged = True
                break
    return dict(P=bestP, val=best, val0=val0, armed=armed, converged=converged, diverged=diverged, steps=t, best_step=best_step)


# ------------------------------------------------------------------------------------------------------------------------------------
# outcomes: strengths, oriented directions, merged counts, recovery
# ------------------------------------------------------------------------------------------------------------------------------------
def neurons(P, depth, d):
    """Oriented unit directions (sign(a_i) w_i / |w_i|) and strengths |a_i| |w_i|."""
    W = P[:H * d].reshape(H, d)
    a = readout(P, depth)
    nw = np.linalg.norm(W, axis=1)
    V = np.sign(a)[:, None] * W / np.maximum(nw, 1e-300)[:, None]
    return V, np.abs(a) * nw


def n_components(C, act, cos):
    """Single-linkage components among the neurons `act` (index array) when two are linked by C[i, j] > cos (C: cosine matrix of the oriented directions)."""
    n = len(act)
    if n == 0:
        return 0
    A = C[np.ix_(act, act)] > cos
    lab = np.arange(n)
    while True:
        new = np.where(A, lab[None, :], n + 1).min(axis=1)
        new = np.minimum(new, lab)
        if (new == lab).all():
            break
        lab = new
    return int(len(np.unique(lab)))


def merged_count(V, strength, tau=TAU, cos=COS_MERGE):
    act = np.flatnonzero(strength > tau)
    return n_components(V @ V.T, act, cos)


def count_grid(V, strength):
    """Merged counts on the grid TAUS x COSES (rows: tau, columns: merge cosine)."""
    C = V @ V.T
    out = np.zeros((len(TAUS), len(COSES)), dtype=np.int16)
    for i, tau in enumerate(TAUS):
        act = np.flatnonzero(strength > tau)
        for j, c in enumerate(COSES):
            out[i, j] = n_components(C, act, c)
    return out


def best_cos(V, strength, Wt, b):
    """For every teacher direction (oriented sign(b_k) w_k) the largest oriented cosine to an active neuron, for each tau in TAUS (-1 if none is active); padded to 8 teachers with nan."""
    T = np.sign(b)[:, None] * Wt
    S = T @ V.T                                                              # K x H
    out = np.full((len(TAUS), 8), np.nan)
    for i, tau in enumerate(TAUS):
        act = strength > tau
        out[i, :len(b)] = S[:, act].max(axis=1) if act.any() else -1.0
    return out


def run_one(stage, K, d, tid, noise, arm, depth, eps, lr, data=None):
    """One training and its outcomes (the record stored in the raw)."""
    data = data or make_data(stage, K, d, tid, noise)
    seed = STAGE_SEED[stage] + 1000 * K + 100 * d + tid
    t0 = time.time()
    res = fit(data, arm, depth, eps, lr, seed)
    P = res['P']
    V, st = neurons(P, depth, d)
    vt = float(np.var(data['yt']))
    test = mse(P, depth, d, data['Xt'].astype(np.float32), data['yt'])
    rec = dict(stage=stage, K=K, d=d, tid=tid, noise=noise, arm=arm, depth=depth, eps=eps, lr=lr, seed=seed, val=res['val'], val0=res['val0'], test_rel=test / vt, armed=res['armed'],
               converged=res['converged'], diverged=res['diverged'], steps=res['steps'], best_step=res['best_step'], secs=time.time() - t0,
               counts=count_grid(V, st), bestcos=best_cos(V, st, data['Wt'], data['b']), n_active=int((st > TAU).sum()), strength_max=float(st.max()))
    return rec


# ------------------------------------------------------------------------------------------------------------------------------------
# jobs, running, raw
# ------------------------------------------------------------------------------------------------------------------------------------
def freeze_key(arm, depth, eps):
    return '%s|%d|%g' % (arm, depth, eps)


def test_member(tid, noise, arm, depth, eps):
    """The reduced test design (plan amendment 2): noise 0 only for the teachers NOISE0_TIDS; (gd, 3, 1e-6) only for the teachers NEVER_TIDS at noise 0.05; arm / depth / eps None: the target-level check."""
    if noise == 0.0 and tid not in NOISE0_TIDS:
        return False
    if arm == 'gd' and depth == 3 and eps == 1e-6 and (noise == 0.0 or tid not in NEVER_TIDS):
        return False
    return True


def make_jobs(stage, freeze=None, tids=None, noises=None):
    """All (K, d, tid, noise, arm, depth, eps, lr) of a stage, grouped by data set (one group = one (K, d, tid, noise)).  dev: every rate of LR_GRID; test: the frozen rate."""
    tids = tids or (DEV_TIDS if stage == 'dev' else TEST_TIDS)
    noises = noises or (DEV_NOISES if stage == 'dev' else NOISES)
    groups = []
    for K, d, tid, noise in itertools.product(KS, DS, tids, noises):
        jobs = []
        if stage == 'test' and not test_member(tid, noise, None, None, None):
            continue
        for arm, depth in COMBOS:
            for eps in EPS:
                if stage == 'test' and not test_member(tid, noise, arm, depth, eps):
                    continue
                lrs = LR_GRID[(arm, depth)] + LR_EXTRA.get((arm, depth, eps), ()) if stage == 'dev' else (freeze['lr'][freeze_key(arm, depth, eps)],)
                jobs += [(K, d, tid, noise, arm, depth, eps, lr) for lr in lrs]
        groups.append(jobs)
    return groups


def job_key(stage, j):
    return '%s|%d|%d|%d|%g|%s|%d|%g|%g' % ((stage,) + tuple(j))


def run_job(args):
    stage, j = args
    rec = run_one(stage, *j)
    rec['key'] = job_key(stage, j)
    return rec


def rec_to_json(r):
    r = dict(r)
    r['counts'] = r['counts'].tolist()
    r['bestcos'] = [[None if np.isnan(x) else float(x) for x in row] for row in r['bestcos']]
    return json.dumps(r, sort_keys=True)


def read_jsonl(path):
    recs = []
    if os.path.exists(path):
        with open(path) as f:
            for line in f:
                line = line.strip()
                if line:
                    recs.append(json.loads(line))
    return recs


def run_stage(stage, out, shard, workers, freeze_path):
    """Run the jobs i, i + n, i + 2n, ... of the flattened job list (cost-balanced by interleaving); append every record to the shard's jsonl (resumable)."""
    freeze = json.load(open(freeze_path)) if stage == 'test' else None
    jobs = [j for g in make_jobs(stage, freeze) for j in g]
    i, n = shard
    jobs = jobs[i::n]
    path = os.path.join(out, 'mlp_%s_part%d_of_%d.jsonl' % (stage, i, n))
    done = {r['key'] for fn in sorted(os.listdir(out)) if fn.startswith('mlp_%s_part' % stage) and fn.endswith('.jsonl') for r in read_jsonl(os.path.join(out, fn))}
    todo = [(stage, j) for j in jobs if job_key(stage, j) not in done]
    print('stage %s shard %d/%d: %d jobs to run, %d done already' % (stage, i, n, len(todo), len(done)), flush=True)
    t0 = time.time()
    with open(path, 'a') as f:
        if workers > 1:
            import multiprocessing as mp
            with mp.Pool(workers) as pool:
                it = pool.imap_unordered(run_job, todo)
                for k, r in enumerate(it):
                    f.write(rec_to_json(r) + '\n')
                    f.flush()
                    print('job %d/%d done, %.0f s' % (k + 1, len(todo), time.time() - t0), flush=True)
        else:
            for k, a in enumerate(todo):
                f.write(rec_to_json(run_job(a)) + '\n')
                f.flush()
                print('job %d/%d done, %.0f s' % (k + 1, len(todo), time.time() - t0), flush=True)
    return path


SCALARS = ('K', 'd', 'tid', 'noise', 'depth', 'eps', 'lr', 'seed', 'val', 'val0', 'test_rel', 'armed', 'converged', 'diverged', 'steps', 'best_step', 'secs', 'n_active', 'strength_max')


def save_raw(path, recs, freeze=None):
    """Raw npz: one array per field, rows sorted by key (so that the raw does not depend on the shard order); the frozen file (rates, eps sets) is stored inside the raw as `freeze_json`."""
    recs = sorted(recs, key=lambda r: r['key'])
    arrs = {k: np.array([r[k] for r in recs], dtype=float) for k in SCALARS}
    arrs['stage'] = np.array([r['stage'] for r in recs])
    arrs['arm'] = np.array([r['arm'] for r in recs])
    arrs['key'] = np.array([r['key'] for r in recs])
    arrs['counts'] = np.array([r['counts'] for r in recs], dtype=np.int16)
    arrs['bestcos'] = np.array([[[np.nan if x is None else x for x in row] for row in r['bestcos']] for r in recs], dtype=float)
    if freeze is not None:
        arrs['freeze_json'] = np.array(json.dumps(freeze, sort_keys=True))
    np.savez_compressed(path, **arrs)


def load_raw(path):
    z = np.load(path)
    return {k: z[k] for k in z.files}


def merge(stage, out, raw, freeze=None):
    recs, seen = [], set()
    for fn in sorted(os.listdir(out)):
        if fn.startswith('mlp_%s_part' % stage) and fn.endswith('.jsonl'):
            for r in read_jsonl(os.path.join(out, fn)):
                if r['key'] not in seen:
                    seen.add(r['key'])
                    recs.append(r)
    save_raw(raw, recs, freeze)
    print('merged', len(recs), 'records into', raw)


# ------------------------------------------------------------------------------------------------------------------------------------
# analysis of a raw (everything below is computed from the raw alone, plus the frozen file for the rates and the eps sets)
# ------------------------------------------------------------------------------------------------------------------------------------
def fmt(x, n=3):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return 'n/a'
    return ('%.' + str(n) + 'f') % x


def fmte(x):
    return 'n/a' if x is None or np.isnan(x) else '%.2e' % x


class Table:
    """The runs of a raw indexed by (arm, depth, eps) and, inside, by target (K, d, tid) and noise."""

    def __init__(self, raw):
        self.raw = raw
        self.idx = {}
        for i in range(len(raw['key'])):
            k = (str(raw['arm'][i]), int(raw['depth'][i]), float(raw['eps'][i]))
            self.idx.setdefault(k, []).append(i)
        for k in self.idx:
            self.idx[k].sort(key=lambda i: (raw['K'][i], raw['d'][i], raw['tid'][i], raw['noise'][i]))
        self.combos = sorted(self.idx, key=lambda k: (COMBOS.index((k[0], k[1])), -k[2]))

    def rows(self, combo):
        return np.array(self.idx[combo], dtype=int)

    def align(self, ca, cb):
        """Rows of combos ca and cb restricted to the (K, d, tid, noise) runs that exist in both, in the same order."""
        r = self.raw

        def keys(rows):
            return {(int(r['K'][i]), int(r['d'][i]), int(r['tid'][i]), float(r['noise'][i])): i for i in rows}
        ka, kb = keys(self.rows(ca)), keys(self.rows(cb))
        common = sorted(set(ka) & set(kb))
        return np.array([ka[k] for k in common], dtype=int), np.array([kb[k] for k in common], dtype=int)

    def tau_i(self, tau):
        return TAUS.index(tau)

    def count_err(self, rows, tau=TAU, cos=COS_MERGE):
        r = self.raw
        return np.abs(r['counts'][rows, TAUS.index(tau), COSES.index(cos)].astype(float) - r['K'][rows])

    def counts(self, rows, tau=TAU, cos=COS_MERGE):
        return self.raw['counts'][rows, TAUS.index(tau), COSES.index(cos)].astype(float)

    def recovery(self, rows, tau=TAU):
        bc = self.raw['bestcos'][rows, TAUS.index(tau), :]
        return np.nanmean(np.where(bc > REC_COS, 1.0, 0.0) + np.where(np.isnan(bc), np.nan, 0.0), axis=1)

    def fit_flag(self, rows):
        return self.raw['test_rel'][rows] < FIT_REL

    def status(self, rows):
        r = self.raw
        div, conv, armed = r['diverged'][rows] > 0, r['converged'][rows] > 0, r['armed'][rows] > 0
        return dict(n=len(rows), diverged=int(div.sum()), converged=int((conv & ~div).sum()), budget_armed=int((armed & ~conv & ~div).sum()), never_armed=int((~armed & ~div).sum()))

    def target_cells(self, rows):
        """Per-target averages over the noise levels, grouped by cell (K, d): {(K, d): (tids, per-row index lists)}; rows must be sorted by (K, d, tid, noise)."""
        r = self.raw
        out = {}
        for j, i in enumerate(rows):
            out.setdefault((int(r['K'][i]), int(r['d'][i])), {}).setdefault(int(r['tid'][i]), []).append(j)
        return {c: [v[t] for t in sorted(v)] for c, v in sorted(out.items())}


def per_target(x, cells):
    """x: values per row (aligned with `rows`); returns {cell: array of per-target means over the noise levels}."""
    return {c: np.array([np.mean(x[js]) for js in tl]) for c, tl in cells.items()}


def boot_mean(pt_by_cell, B=N_BOOT, seed=BOOT_SEED):
    """Stratified (by cell) cluster bootstrap over targets of the mean over cells of the cell mean; returns (mean, lo, hi) with the 2.5 / 97.5 percentiles."""
    rng = np.random.default_rng(seed)
    cells = [pt_by_cell[c] for c in sorted(pt_by_cell)]
    point = float(np.mean([x.mean() for x in cells]))
    tot = np.zeros(B)
    for x in cells:
        tot += x[rng.integers(0, len(x), (B, len(x)))].mean(axis=1)
    tot /= len(cells)
    lo, hi = np.percentile(tot, [2.5, 97.5])
    return point, float(lo), float(hi)


def contrast(T, ca, cb, f, B=N_BOOT, seed=BOOT_SEED):
    """Paired contrast of the per-run statistic f (callable (table, rows) -> array) between combos ca and cb: mean of f(ca) - f(cb) with its bootstrap interval; the runs of ca and cb are paired
    row by row (the same target and noise); only the (target, noise) pairs run in both combos enter."""
    ra, rb = T.align(ca, cb)
    assert len(ra) > 0
    cells = T.target_cells(ra)
    d = f(T, ra) - f(T, rb)
    return boot_mean(per_target(d, cells), B, seed)


def f_ce(tau=TAU, cos=COS_MERGE):
    return lambda T, rows: T.count_err(rows, tau, cos)


def f_test(T, rows):
    return T.raw['test_rel'][rows]


def holds_gain(m, lo, margin=MARGIN_CE):
    """A count-error gain (other minus ours) is established if its interval is above 0 and its mean is at least `margin`."""
    return bool(lo > 0 and m >= margin)


def noninferior(m, lo, hi, margin=ABS_T):
    return bool(hi <= margin)


def load_freeze(path):
    with open(path) as f:
        return json.load(f)


def report(raw_path, freeze_path=None):
    """The report of a test raw; the frozen rates and eps sets are those stored in the raw (or, if `freeze_path` is given, that file)."""
    raw = load_raw(raw_path)
    fr = load_freeze(freeze_path) if freeze_path else json.loads(str(raw['freeze_json']))
    T = Table(raw)
    L = []
    P = L.append
    nrun = len(raw['key'])
    tids = sorted(set(int(x) for x in raw['tid']))
    E_R, E_O = fr['E_R'], fr['E_O']
    P('# MLP student: does small-start gradient training of a depth-p readout leave unneeded hidden units dormant (and do Adam-type optimizers erase it)?')
    P('')
    P('Generated by `experiments/exp_mlp_dormancy.py --from-raw`; plan: `results/ml/mlp_plan.md`; frozen rates and eps sets: `results/ml/mlp_freeze.json` (stored inside the raw). Every number below is computed from the stored raw.')
    P('')
    P('## 1. Design and size')
    P('')
    P('* Teacher `sum_k b_k tanh(w_k . x)`, K in %s, d in %s; %d test targets (teacher ids %d-%d) per (K, d) at label noise sd %g, and noise 0 for the first %d of them (the same teacher and student initialisation); (gd, depth 3, eps 1e-6) only for the first %d teachers per cell (plan amendment 2): %d runs in all.' % (
        KS, DS, len(tids), tids[0], tids[-1], NOISES[0], len(NOISE0_TIDS), len(NEVER_TIDS), nrun))
    P('* Student: H = %d, n = %d, validation %d (noisy), test %d (clean teacher). Step budget %d, validation every %d steps, patience max(%d, step of the last improvement), armed once the validation MSE falls below %.1f x its initial value.' % (
        H, N_TR, N_VAL, N_TE, TMAX, EVERY, PATIENCE, ARM_FRAC))
    P('* Declared outcome: a neuron is active if |a_i| |w_i| > %g; active neurons whose oriented directions sign(a_i) w_i / |w_i| have cosine > %g are merged (single linkage); count error = |merged count - K|; relative test MSE = clean test MSE / teacher variance; a run "fits" if it is below %g; a teacher direction is recovered if some active neuron has oriented cosine > %g to it.' % (
        TAU, COS_MERGE, FIT_REL, REC_COS))
    P('* Rates (frozen on development targets, one rate per arm, depth and eps): ' + '; '.join('%s: %s' % (k, fr['lr'][k]) for k in sorted(fr['lr'])))
    P('* eps sets fixed on development: E_R (R is tested at) = %s, E_O (O and T are tested at) = %s.' % (E_R, E_O))
    P('')
    P('## 2. Convergence (every run is in every table; none is dropped)')
    P('')
    P('| arm | depth | eps | rate | runs | converged | armed, budget exhausted | never left the start | diverged |')
    P('|---|---|---|---|---|---|---|---|---|')
    for c in T.combos:
        st = T.status(T.rows(c))
        P('| %s | %d | %g | %g | %d | %d | %d | %d | %d |' % (c[0], c[1], c[2], fr['lr'][freeze_key(*c)], st['n'], st['converged'], st['budget_armed'], st['never_armed'], st['diverged']))
    P('')
    P('## 3. Count and accuracy at the declared setting (means over all runs)')
    P('')
    P('| arm | depth | eps | merged count | count error | active before merging | median rel. test MSE | fits | recovered share |')
    P('|---|---|---|---|---|---|---|---|---|')
    for c in T.combos:
        rows = T.rows(c)
        P('| %s | %d | %g | %s | %s | %s | %s | %s | %s |' % (c[0], c[1], c[2], fmt(T.counts(rows).mean(), 2), fmt(T.count_err(rows).mean(), 2), fmt(raw['n_active'][rows].mean(), 1),
                                                           fmte(np.median(raw['test_rel'][rows])), fmt(T.fit_flag(rows).mean(), 2), fmt(np.nanmean(T.recovery(rows)), 2)))
    P('')
    P('Count error by K (mean over d, tid and noise):')
    P('')
    P('| arm | depth | eps | ' + ' | '.join('K = %d' % k for k in KS) + ' |')
    P('|---|---|---|' + '---|' * len(KS))
    for c in T.combos:
        rows = T.rows(c)
        ce = T.count_err(rows)
        P('| %s | %d | %g | ' % c + ' | '.join(fmt(ce[raw['K'][rows] == k].mean(), 2) for k in KS) + ' |')
    P('')
    V = {}
    # ---- R
    P('## 4. (R) depth of the readout, plain GD')
    P('')
    P('Paired contrast (count error of depth p minus count error of depth 3), mean over targets with a 95%% cluster-bootstrap interval (%d resamples over targets within each (K, d), cells weighted equally); gain established if the interval is above 0 and the mean is at least %g neuron.' % (N_BOOT, MARGIN_CE))
    P('')
    P('| eps | contrast | mean | 95% interval | established | same, runs where both fit: n, mean, interval |')
    P('|---|---|---|---|---|---|')
    r_ok = {}
    for eps in EPS:
        for p in (2, 1):
            ca, cb = ('gd', p, eps), ('gd', 3, eps)
            m, lo, hi = contrast(T, ca, cb, f_ce())
            ok = holds_gain(m, lo)
            r_ok[(eps, p)] = ok
            ra, rb = T.align(ca, cb)
            both = T.fit_flag(ra) & T.fit_flag(rb)
            if both.sum() >= 10:
                d = (T.count_err(ra) - T.count_err(rb))[both]
                cells = T.target_cells(ra)
                dd = np.full(len(ra), np.nan)
                dd[both] = d
                pt = {}
                for c, tl in cells.items():
                    vals = [np.nanmean(dd[js]) for js in tl if not np.all(np.isnan(dd[js]))]
                    pt[c] = np.array(vals)
                if all(len(v) > 0 for v in pt.values()):
                    mc, loc, hic = boot_mean(pt)
                    cond = '%d, %s, [%s, %s]' % (both.sum(), fmt(mc, 2), fmt(loc, 2), fmt(hic, 2))
                else:
                    cond = '%d runs, a cell is empty' % both.sum()
            else:
                cond = '%d runs (too few)' % both.sum()
            P('| %g | depth %d - depth 3 | %s | [%s, %s] | %s | %s |' % (eps, p, fmt(m, 2), fmt(lo, 2), fmt(hi, 2), 'yes' if ok else 'no', cond))
    V['R'] = bool(E_R) and all(r_ok[(e, p)] for e in E_R for p in (1, 2))
    P('')
    n1 = {e: float(np.median(raw['n_active'][T.rows(('gd', 1, e))])) for e in EPS}
    P('Depth 1 (GD): median number of active neurons before merging, by eps: ' + ', '.join('%g: %s' % (e, fmt(n1[e], 1)) for e in EPS) + ' (H = %d; registered: at least %g H = %g at eps <= 1e-4).' % (H, HALF_H, HALF_H * H))
    V['R_H'] = all(n1[e] >= HALF_H * H for e in EPS if e <= 1e-4)
    P('')
    P('**(R) verdict** (every eps in E_R = %s, both contrasts established): **%s**.  Depth 1 forms about H neurons: **%s**.' % (E_R, 'HELD' if V['R'] else 'FAILED', 'HELD' if V['R_H'] else 'FAILED'))
    P('')
    # ---- O
    P('## 5. (O) optimizers at depth 3')
    P('')
    P('Paired contrast (count error of the comparator minus count error of the dormancy-keeping arm), same convention. Registered: GD and shared-v keep dormancy; Adam and sign descent erase it.')
    P('')
    P('| eps | contrast | mean | 95% interval | established |')
    P('|---|---|---|---|---|')
    pairs = (('adam', 'gd'), ('sign', 'gd'), ('adam', 'shared_v'), ('sign', 'shared_v'))
    o_ok = {}
    for eps in EPS:
        for a, b in pairs:
            m, lo, hi = contrast(T, (a, 3, eps), (b, 3, eps), f_ce())
            o_ok[(eps, a, b)] = holds_gain(m, lo)
            P('| %g | %s - %s | %s | [%s, %s] | %s |' % (eps, a, b, fmt(m, 2), fmt(lo, 2), fmt(hi, 2), 'yes' if o_ok[(eps, a, b)] else 'no'))
    V['O'] = bool(E_O) and all(o_ok[(e, a, b)] for e in E_O for a, b in pairs)
    P('')
    P('Adam with a large eps_adam (depth 3): count error %s; shared-v: %s; GD: %s (mean by eps %s).' % (
        ', '.join('%g: %s' % (e, fmt(T.count_err(T.rows(('adam_big', 3, e))).mean(), 2)) for e in EPS),
        ', '.join('%g: %s' % (e, fmt(T.count_err(T.rows(('shared_v', 3, e))).mean(), 2)) for e in EPS),
        ', '.join('%g: %s' % (e, fmt(T.count_err(T.rows(('gd', 3, e))).mean(), 2)) for e in EPS), list(EPS)))
    P('')
    P('**(O) verdict at the declared setting** (every eps in E_O = %s, all four contrasts established): **%s**.' % (E_O, 'HELD' if V['O'] else 'FAILED'))
    P('')
    P('### 5a. Sweep of the merge cosine and the activity threshold (the gap Adam - GD and sign - GD in count error)')
    P('')
    P('Cells show the mean paired difference (comparator minus GD) of the count error; `*` marks a cell whose 95%% interval (%d resamples) is above 0, `-` one whose interval is below 0, nothing an interval containing 0. Rows: activity threshold tau; columns: merge cosine (1.01 = no merging); the declared cell is tau = %g, cos = %g.' % (2000, TAU, COS_MERGE))
    sweep = {}
    for eps in EPS:
        for a, b in (('adam', 'gd'), ('sign', 'gd'), ('adam', 'shared_v'), ('sign', 'shared_v')):
            g = np.zeros((len(TAUS), len(COSES)))
            sg = np.zeros((len(TAUS), len(COSES)), dtype=int)
            for i, tau in enumerate(TAUS):
                for j, cs in enumerate(COSES):
                    m, lo, hi = contrast(T, (a, 3, eps), (b, 3, eps), f_ce(tau, cs), B=2000, seed=BOOT_SEED + 1)
                    g[i, j], sg[i, j] = m, (1 if lo > 0 else (-1 if hi < 0 else 0))
            sweep[(eps, a, b)] = (g, sg)
        for a, b in (('adam', 'gd'), ('sign', 'gd')):
            g, sg = sweep[(eps, a, b)]
            P('')
            P('eps = %g, %s - %s:' % (eps, a, b))
            P('')
            P('| tau \\ cos | ' + ' | '.join('%g' % c for c in COSES) + ' |')
            P('|---|' + '---|' * len(COSES))
            for i, tau in enumerate(TAUS):
                P('| %g | ' % tau + ' | '.join(fmt(g[i, j], 2) + ('*' if sg[i, j] > 0 else ('-' if sg[i, j] < 0 else '')) for j in range(len(COSES))) + ' |')
    P('')
    P('Share of the %d sweep cells with the comparator worse than the dormancy-keeping arm (interval above 0), `+`, and better (interval below 0), `-`:' % (len(TAUS) * len(COSES)))
    P('')
    P('| eps | contrast | worse (+) | better (-) |')
    P('|---|---|---|---|')
    for eps in EPS:
        for a, b in pairs:
            g, sg = sweep[(eps, a, b)]
            P('| %g | %s - %s | %s | %s |' % (eps, a, b, fmt((sg > 0).mean(), 2), fmt((sg < 0).mean(), 2)))
    V['O_sweep'] = bool(E_O) and all((sweep[(e, a, b)][1] > 0).mean() >= SWEEP_SHARE for e in E_O for a, b in pairs)
    P('')
    P('**(O) sweep verdict** (at every eps in E_O and for all four contrasts at least %g of the cells have the comparator worse, interval above 0): **%s**.' % (SWEEP_SHARE, 'HELD' if V['O_sweep'] else 'FAILED'))
    P('')
    # ---- T
    P('## 6. (T) accuracy of the depth-3 GD student')
    P('')
    P('Paired difference of the relative clean test MSE (depth-3 GD minus the comparator); non-inferior if the upper end of the 95%% interval is at most %g (1%% of the teacher variance).' % ABS_T)
    P('')
    P('| eps | comparator | mean difference | 95% interval | ratio of medians | non-inferior |')
    P('|---|---|---|---|---|---|')
    t_ok = {}
    comps = [('gd', 1), ('gd', 2), ('adam', 3), ('shared_v', 3), ('sign', 3), ('adam_big', 3)]
    for eps in EPS:
        for arm, dp in comps:
            m, lo, hi = contrast(T, ('gd', 3, eps), (arm, dp, eps), f_test)
            ok = noninferior(m, lo, hi)
            t_ok[(eps, arm, dp)] = ok
            ra, rb = T.align(('gd', 3, eps), (arm, dp, eps))
            P('| %g | %s, depth %d | %s | [%s, %s] | %s | %s |' % (eps, arm, dp, fmte(m), fmte(lo), fmte(hi), fmt(np.median(raw['test_rel'][ra]) / np.median(raw['test_rel'][rb]), 2), 'yes' if ok else 'no'))
    V['T'] = bool(E_O) and all(t_ok[(e, a, d)] for e in E_O for a, d in comps)
    P('')
    P('**(T) verdict** (every eps in E_O = %s, every comparator): **%s**.' % (E_O, 'HELD' if V['T'] else 'FAILED'))
    P('')
    P('## 7. Verdicts')
    P('')
    P('| prediction | verdict |')
    P('|---|---|')
    P('| (R) depth 3 GD has a smaller count error than depth 2 and depth 1 at every eps in E_R = %s | %s |' % (E_R, 'held' if V['R'] else 'failed'))
    P('| (R2) depth 1 forms at least %g H neurons at eps <= 1e-4 | %s |' % (HALF_H, 'held' if V['R_H'] else 'failed'))
    P('| (O) at depth 3 GD and shared-v have a smaller count error than Adam and sign descent at every eps in E_O = %s | %s |' % (E_O, 'held' if V['O'] else 'failed'))
    P('| (O-sweep) the same contrasts are above 0 in at least %g of the tau x cos sweep cells at every eps in E_O | %s |' % (SWEEP_SHARE, 'held' if V['O_sweep'] else 'failed'))
    P('| (T) depth-3 GD test MSE not worse than any comparator by more than %g (of the teacher variance) at every eps in E_O | %s |' % (ABS_T, 'held' if V['T'] else 'failed'))
    P('')
    return '\n'.join(L) + '\n'


def dev_choose(raw):
    """Development: for every (arm, depth, eps) the rate with the lowest mean log10(val / val0) over the development runs (ties: the smaller rate).  Returns the freeze dict."""
    T = Table(raw)
    lr, scores, edge = {}, {}, []
    for c in T.combos:
        rows = T.rows(c)
        lrs = sorted(set(float(x) for x in raw['lr'][rows]))
        sc = []
        for x in lrs:
            rr = rows[raw['lr'][rows] == x]
            sc.append(float(np.mean(np.log10(raw['val'][rr] / raw['val0'][rr]))))
        best = int(np.argmin(np.round(sc, 9)))
        lr[freeze_key(*c)] = lrs[best]
        scores[freeze_key(*c)] = dict(zip(['%g' % x for x in lrs], sc))
        if best in (0, len(lrs) - 1):
            edge.append(freeze_key(*c))
    fits = {}
    for eps in EPS:
        rows = T.rows(('gd', 3, eps))
        rr = rows[raw['lr'][rows] == lr[freeze_key('gd', 3, eps)]]
        fits[eps] = float(np.mean(raw['test_rel'][rr] < FIT_REL))
    E_O = [e for e in EPS if fits[e] >= FIT_DEV]
    E_R = [e for e in E_O if e <= 1e-4]
    return dict(lr=lr, E_O=E_O, E_R=E_R, gd3_fit_share={'%g' % e: fits[e] for e in EPS}, edge=edge, scores=scores)


def dev_report(raw_path):
    raw = load_raw(raw_path)
    fr = dev_choose(raw)
    T = Table(raw)
    L = ['# MLP dormancy: development (rate tuning), development targets only', '',
         'Generated by `experiments/exp_mlp_dormancy.py --dev-report`. Development: teacher id %s, K in %s, d in %s, noise %s, %d runs. Score = mean over the development runs of log10(validation MSE of the kept iterate / initial validation MSE) (lower is better); the rate with the lowest score is frozen (ties: the smaller rate).' % (
             DEV_TIDS, KS, DS, DEV_NOISES, len(raw['key'])), '',
         '| arm | depth | eps | rates: score (chosen rate in bold) | diverged runs at the chosen rate | never armed at the chosen rate | fit share at the chosen rate |', '|---|---|---|---|---|---|---|']
    for c in T.combos:
        k = freeze_key(*c)
        rows = T.rows(c)
        rr = rows[raw['lr'][rows] == fr['lr'][k]]
        cells = ', '.join(('**%s: %s**' if float(x) == fr['lr'][k] else '%s: %s') % (x, fmt(v, 2)) for x, v in fr['scores'][k].items())
        L.append('| %s | %d | %g | %s | %d/%d | %d/%d | %s |' % (c[0], c[1], c[2], cells, int((raw['diverged'][rr] > 0).sum()), len(rr), int((raw['armed'][rr] == 0).sum()), len(rr), fmt(np.mean(raw['test_rel'][rr] < FIT_REL), 2)))
    L += ['', 'Frozen rates at the edge of their grid: %s.' % (', '.join(fr['edge']) or 'none'), '',
          'Fit share of depth-3 GD at its frozen rate by eps: %s; eps feasible for depth-3 GD (share >= %g): E_O = %s; E_R (E_O with eps <= 1e-4) = %s.' % (fr['gd3_fit_share'], FIT_DEV, fr['E_O'], fr['E_R']), '']
    return '\n'.join(L) + '\n', fr


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--stage', choices=('dev', 'test'))
    ap.add_argument('--out', default=DEFAULT_OUT)
    ap.add_argument('--shard', default='0/1')
    ap.add_argument('--workers', type=int, default=1)
    ap.add_argument('--freeze', default=os.path.join(DEFAULT_OUT, 'mlp_freeze.json'))
    ap.add_argument('--merge', action='store_true')
    ap.add_argument('--raw')
    ap.add_argument('--from-raw', dest='from_raw', nargs='?', const='')
    ap.add_argument('--md')
    ap.add_argument('--dev-report', dest='dev_report')
    ap.add_argument('--write-freeze', dest='write_freeze')
    A = ap.parse_args(argv)
    A.freeze_given = '--freeze' in (argv if argv is not None else sys.argv)
    if A.merge:
        fr = load_freeze(A.freeze) if A.stage == 'test' else None
        merge(A.stage, A.out, A.raw or os.path.join(A.out, 'mlp_%s_raw.npz' % A.stage), fr)
        return
    if A.dev_report:
        md, fr = dev_report(A.dev_report)
        with open(A.md or os.path.join(A.out, 'mlp_dev.md'), 'w') as f:
            f.write(md)
        if A.write_freeze:
            with open(A.write_freeze, 'w') as f:
                json.dump(fr, f, indent=1, sort_keys=True)
                f.write('\n')
        return
    if A.from_raw is not None:
        raw = A.from_raw or os.path.join(A.out, 'mlp_raw.npz')
        md = A.md or os.path.join(A.out, 'mlp.md')
        with open(md, 'w') as f:
            f.write(report(raw, A.freeze if A.freeze_given else None))
        print('wrote', md)
        return
    i, n = (int(x) for x in A.shard.split('/'))
    os.makedirs(A.out, exist_ok=True)
    run_stage(A.stage, A.out, (i, n), A.workers, A.freeze)


if __name__ == '__main__':
    main()
