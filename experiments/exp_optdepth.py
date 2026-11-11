"""
Two pre-registered prediction tests of the theory: the optimizer (P1) and the depth of the readout (P2)
======================================================================================================

The theory: small-start training of a factorised Gaussian-unit model chooses its own unit count when (R) the readout has effective depth >= 3 and (O) the
optimizer does not normalise each coordinate by its own gradient history.  Plan, registries and decision rules: results/ml/optdepth_plan.md (committed before any
run on the held-out targets); the report is results/ml/optdepth.md (+ optdepth_raw.npz); this script regenerates it byte-identically from the raw (--from-raw).

P1 (optimizer).  The additive Gaussian-unit KAN [5, 1] of the KAN spike (experiments/exp_kan_spike.py: `make_target`, `sample_synth`, the activity rule), in
numpy (torch is not installed; the model is the numpy re-implementation `AddKan` of experiments/theory3f_escape.py, extended here by the 'free' control).  Arms:
plain GD (the frozen configuration of results/ml/kan_freeze.json), Adam (default eps_adam), Adam with a large eps_adam (GD-like rates), Adam with ONE shared
scalar second moment, sign descent, and the control (free weights, L1, magnitude pruning, fine-tune).  Every method's active units are merged by single linkage
at one radius (0.3, declared beforehand).  Primary metric: exact-edge share on the edges with K > 0 (cluster bootstrap over the targets).

P2 (depth).  (i) the dormancy exponents |w_dormant| ~ eps^(p/2) at a fixed margin: a scalar unit with a prescribed drive and the repository's flow (moving centres,
unequal masses); (ii) the fraction of the 200 benchmark targets of the readouts study (`exp_parametrization.sample_targets`) on which the formed set of the
depth-p readout equals the sliding path's, p = 2, 3, 4, 5, eps = 1e-6, 1e-8.

    python experiments/exp_optdepth.py --check-unused
    python experiments/exp_optdepth.py --run p1|p2 --stage dev|test --work DIR [--freeze results/ml/optdepth_freeze.json]
    python experiments/exp_optdepth.py --merge --stage dev|test --work DIR --raw RAW.npz
    python experiments/exp_optdepth.py --freeze-from RAW_DEV.npz --out results/ml/optdepth_freeze.json
    python experiments/exp_optdepth.py --from-raw RAW.npz --md OUT.md          (numpy only; byte-identical regeneration)
"""
import argparse
import json
import os
import sys
import time

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import exp_kan_spike as KS  # noqa: E402
from exp_kan_merge import merged_count  # noqa: E402

OUT = os.environ.get('OD_OUT', os.path.join(REPO, 'results', 'ml'))

# ------------------------------------------------------------------------------------------------------------------------------------
# registries (everything the plan fixes)
# ------------------------------------------------------------------------------------------------------------------------------------
DEV_TIDS, DEV_SEEDS = (401, 402, 403), (40, 41)                  # development targets / data seeds (fresh: see check_unused)
TEST_TIDS, TEST_SEEDS = tuple(range(411, 431)), (42, 43, 44)     # held-out targets / data seeds
NOISES = (0.1, 0.4)
USED_TIDS = tuple(range(1, 11)) + (101, 102, 103) + tuple(range(301, 307))     # the ids of the spike (test, development), theory3f, the referee's reruns
USED_SEEDS = (0, 1, 2, 3)                                         # the data seeds of the spike (0, 1, 2), theory3f (0), the referee's reruns (1, 2, 3)
M_UNITS, EPS0 = KS.M_UNITS, 1e-3
STEPS, PATIENCE, EVERY, FT_STEPS = KS.STEPS, KS.PATIENCE, KS.EVERY, KS.FT_STEPS
A_CAP = KS.AMP_CAP
PZ_SIG, PZ_LRC = 0.12, 12.0                                      # exp_parametrization: kernel width SIG and the rate Gm of the centres (and of the linear readout D1)
LAMS = KS.LAMS
CTRL_LR = 0.03                                                   # kan_freeze.json lr_synth.ctrl
GD_RATES = dict(a=1.0, n=0.05, c=1.0, ls=1.0, b=0.1)             # kan_freeze.json ours_synth: lr 1, nu_ratio 0.05, bias 0.1 (exp_kan_spike.run_methods)
GD_CFG = dict(lr=1.0)
MERGE_R = 0.3                                                    # the declared merge radius of the primary metric
RADII = (0.0, 0.1, 0.2, 0.3, 0.4, 0.5)
N_BOOT, BOOT_SEED = 20000, 20261110
ARMS = ('gd', 'adam', 'adam_big', 'shared_v', 'sign')
ARM_LABEL = dict(gd='plain GD (frozen configuration)', adam='Adam (default eps_adam = 1e-8)', adam_big='Adam, large eps_adam (GD-like rates)',
                 shared_v='Adam, one shared second moment', sign='sign descent', ctrl='control: L1 + prune (penalty by validation, one-SE)',
                 ctrl_or='control: L1 + prune (penalty by the truth at the radius; oracle)')
DEV_GRID = dict(adam=[dict(lr=l) for l in (0.003, 0.01, 0.03, 0.1)],                                   # amendment 2: the first grids (0.003-0.03; 0.01-1; 0.003-0.03; 0.0003-0.01) had their best
                adam_big=[dict(ea=e) for e in (0.01, 0.1, 1.0, 10.0)],                                  # value at an edge; one more value beyond each edge was added (0.1; 10; 0.001; 0.0001)
                shared_v=[dict(lr=l, pre=p) for l in (0.001, 0.003, 0.01, 0.03) for p in (False, True)],
                sign=[dict(lr=l) for l in (0.0001, 0.0003, 0.001, 0.003, 0.01)])
GATE_GAP, GATE_NEAR = 0.20, 0.10                                 # prediction (a): gap at least 0.20; (b): GD and shared-v within 0.10

# P2
EPS_I = (1e-4, 1e-6, 1e-8, 1e-10)                                # (i) the eps grid of the exponent fit
DEPTHS_I = (3, 4, 5)
MUS_SCALAR = (0.25, 0.5)
TOL_EXP = 0.05                                                   # exponent within 0.05 of p/2
TOL_SCALAR = 0.01
RAND_NET = 0.5                                                   # unequal masses: b_i = exp(0.5 N(0, 1)) (exp_parametrization part G)
MASS_SEED = 4242
MARGIN_MIN = 0.05                                                # a target enters the exponent fit only if the realised margin of the dormant maximum is at least this at every eps (amendment 1: 0.10 -> 0.05)
MIN_TARGETS = 4                                                  # fewer included targets than this: the exponent verdict of the depth is INCONCLUSIVE (amendment 1)
M_NET = 8                                                        # targets of the 200-sample used in (i) besides the design
MARGIN_PICK = 0.15                                               # (i): targets of the sample with a path runner-up margin of at least this
DEPTHS_II = (2, 3, 4, 5)
EPS_II = (1e-6, 1e-8)
EPS_II_SEC = (1e-3, 1e-4)                                        # secondary eps of (ii), where the depths are not yet saturated (amendment 1)
EPS_II_ALL = EPS_II_SEC + EPS_II
N_BOOT_II = 20000


# ------------------------------------------------------------------------------------------------------------------------------------
# unused-ids check
# ------------------------------------------------------------------------------------------------------------------------------------
def data_seed(tid, noise, seed):
    """The integer seed of KS.sample_synth's generator (a function of the target id, noise and data seed)."""
    return 10 ** 6 + 1000 * int(tid) + 10 * int(round(noise * 100)) + int(seed)


def check_unused():
    """Verify that the new ids / seeds were used by no earlier file: (1) not in USED_TIDS (the ids of kan_plan.md, kan_dev.md, kan_raw.npz, kan_dev_raw.npz, theory3f's
    TARGETS_K, the referee's reruns 101-103 / 301-306); (2) the generating seed of every new target and of every new data draw is different from every earlier one;
    (3) no new target has the bumps (or the per-edge counts and centres) of an earlier target. Returns a list of text lines and a bool."""
    lines, ok = [], True
    new_t = tuple(DEV_TIDS) + tuple(TEST_TIDS)
    clash = sorted(set(new_t) & set(USED_TIDS))
    ok &= not clash
    lines.append('new ids %s; earlier ids %s; overlap: %s' % ('%d-%d, %d-%d' % (DEV_TIDS[0], DEV_TIDS[-1], TEST_TIDS[0], TEST_TIDS[-1]), 'spike 1-10, 101-103, theory3f 101-103, referee 101-103 / 301-306', clash or 'none'))
    old_seeds = {1000 + t for t in USED_TIDS}
    new_gen = {1000 + t for t in new_t}
    ok &= not (old_seeds & new_gen)
    lines.append('target-generator seeds 1000 + id: new %d-%d, earlier %s; overlap: %s' % (min(new_gen), max(new_gen), sorted(old_seeds), sorted(old_seeds & new_gen) or 'none'))
    old_data = {data_seed(t, n, s) for t in USED_TIDS for n in NOISES for s in USED_SEEDS}
    new_data = {data_seed(t, n, s) for t in new_t for n in NOISES for s in tuple(DEV_SEEDS) + tuple(TEST_SEEDS)}
    ok &= not (old_data & new_data) and not (set(DEV_SEEDS) | set(TEST_SEEDS)) & set(USED_SEEDS)
    lines.append('data-draw seeds: %d earlier combinations (ids x noises x seeds %s), %d new combinations (seeds %s, %s); overlap: %d' % (
        len(old_data), USED_SEEDS, len(new_data), DEV_SEEDS, TEST_SEEDS, len(old_data & new_data)))
    old_b = {t: KS.make_target(t) for t in USED_TIDS}
    same = []
    for t in new_t:
        tg = KS.make_target(t)
        for u, o in old_b.items():
            if np.array_equal(tg['K'], o['K']) and all(a.shape == b.shape and np.allclose(a, b) for a, b in zip(tg['bumps'], o['bumps'])):
                same.append((t, u))
    ok &= not same
    lines.append('generated targets: %d new against %d earlier, identical bumps: %s' % (len(new_t), len(old_b), same or 'none'))
    ks = sorted(tuple(KS.make_target(t)['K'].tolist()) for t in new_t)
    lines.append('K vectors of the new targets are all different from each other: %s' % (len(set(ks)) == len(ks)))
    return lines, bool(ok)


# ------------------------------------------------------------------------------------------------------------------------------------
# P1: the model and the optimizers (numpy)
# ------------------------------------------------------------------------------------------------------------------------------------
class Opt:
    """One optimizer on the vector of ALL parameters (callers pass the NEGATIVE gradient): 'gd' p += lr*pre*g; 'sign' p += lr*sign(g);
    'adam' p += lr*mhat/(sqrt(vhat) + eps_adam) (bias corrected, per-coordinate second moment, the PyTorch rule); 'adam_g' the same with ONE scalar second
    moment v = EMA of mean((pre*g)^2); lr may be a vector."""

    def __init__(self, kind, lr, beta1=0.9, beta2=0.999, eps_adam=1e-8, pre=None):
        assert kind in ('gd', 'sign', 'adam', 'adam_g')
        self.kind, self.lr, self.b1, self.b2, self.ea, self.pre = kind, lr, float(beta1), float(beta2), float(eps_adam), pre
        self.m = self.v = None
        self.t = 0

    def step(self, p, g):
        self.t += 1
        k = self.kind
        if self.pre is not None and k in ('gd', 'adam_g'):
            g = self.pre * g
        if k == 'gd':
            return p + self.lr * g
        if k == 'sign':
            return p + self.lr * np.sign(g)
        if self.m is None:
            self.m, self.v = np.zeros_like(p), np.zeros_like(p)
        self.m = self.b1 * self.m + (1 - self.b1) * g
        gg = g ** 2 if k == 'adam' else np.mean(g ** 2) * np.ones_like(g)
        self.v = self.b2 * self.v + (1 - self.b2) * gg
        mh, vh = self.m / (1 - self.b1 ** self.t), self.v / (1 - self.b2 ** self.t)
        den = np.sqrt(vh) + self.ea
        return p + self.lr * np.divide(mh, den, out=np.zeros_like(den), where=den > 0)


class AddNet:
    """The additive KAN [d, 1]: f(x) = bias + sum_jl w_jl exp(-(x_j - c_jl)^2/(2 sigma_jl^2)), centres c and log-widths trained (clamped to [log 0.05, log 3]),
    initial centres on a uniform grid of [X_LO, X_HI], width = spacing; loss 0.5 * mean r^2 (+ l1 * sum |w|).
      mode 'fact': w = nu^2 A tanh(alpha/A), alpha(0) = 0, nu(0)^2 = eps (the factorised small start);  mode 'free': w direct, N(0, 0.1^2), prunable by `mask`."""

    def __init__(self, d, M, mode, eps=EPS0, seed=0):
        assert mode in ('fact', 'free')
        grid = np.linspace(KS.X_LO, KS.X_HI, M)
        self.mode, self.A = mode, A_CAP
        self.p = dict(c=np.tile(grid, (d, 1)), ls=np.full((d, M), np.log(grid[1] - grid[0])), b=np.zeros(1))
        if mode == 'fact':
            self.p.update(a=np.zeros((d, M)), n=np.full((d, M), np.sqrt(eps)))
            self.keys = ('a', 'n', 'c', 'ls', 'b')
        else:
            self.p.update(w=0.1 * np.random.default_rng(int(seed)).standard_normal((d, M)))
            self.keys = ('w', 'c', 'ls', 'b')
        self.mask = np.ones((d, M))

    def weights(self):
        p = self.p
        if self.mode == 'fact':
            return p['n'] ** 2 * self.A * np.tanh(p['a'] / self.A) * self.mask
        return p['w'] * self.mask

    def widths(self):
        return np.exp(np.clip(self.p['ls'], np.log(0.05), np.log(3.0)))

    def fwd(self, X):
        s = self.widths()
        u = (X[:, :, None] - self.p['c'][None]) / s[None]
        phi = np.exp(-0.5 * u * u)
        w = self.weights()
        return (phi * w[None]).sum((1, 2)) + self.p['b'][0], phi, u, s, w

    def grad(self, X, y, l1=0.0):
        """The NEGATIVE gradient of 0.5 * mean r^2 + l1 * sum |w| as a dict over the keys."""
        p = self.p
        f, phi, u, s, w = self.fwd(X)
        r = y - f
        rp = r[:, None, None] * phi
        gw = rp.mean(0)
        inside = (p['ls'] > np.log(0.05)) & (p['ls'] < np.log(3.0))
        out = dict(c=(rp * u / s[None]).mean(0) * w, ls=(rp * u * u).mean(0) * w * inside, b=np.array([r.mean()]))
        if self.mode == 'fact':
            th = np.tanh(p['a'] / self.A)
            out.update(a=gw * p['n'] ** 2 * (1 - th ** 2), n=gw * 2 * p['n'] * self.A * th)
        else:
            out['w'] = (gw - l1 * np.sign(p['w'])) * self.mask
        return out

    def vec(self, d):
        return np.concatenate([d[k].ravel() for k in self.keys])

    def setvec(self, z):
        i = 0
        for k in self.keys:
            n = self.p[k].size
            self.p[k] = z[i:i + n].reshape(self.p[k].shape)
            i += n

    def rates_vec(self, rates):
        return self.vec({k: np.full_like(self.p[k], rates[k]) for k in self.keys})

    def activity(self, X, tau=KS.TAU):
        """(rms, sd, active) of the units on the first 2000 rows of X (the spike's `layer_activity`), arrays of shape (1, d, M)."""
        return KS.layer_activity(self.weights()[None], self.p['c'][None], self.widths()[None], X[:2000], tau)

    def counts(self, X, tau=KS.TAU):
        return KS.edge_counts(self.activity(X, tau)[2])


def fit(model, data, opt, steps, l1=0.0, patience=PATIENCE):
    """The spike's `fit`: full batch, validation MSE every EVERY steps, the best iterate kept, stop when the validation MSE has not improved for `patience` steps counted
    from the first time it is below 0.9 of its initial value. Returns (steps, best validation MSE)."""
    Xtr, ytr, Xva, yva = data['Xtr'], data['ytr'], data['Xva'], data['yva']
    best, best_p, v0, armed, since = np.inf, None, None, False, 0
    t = 0
    for t in range(steps + 1):
        if t % EVERY == 0:
            v = float(np.mean((model.fwd(Xva)[0] - yva) ** 2))
            if not np.isfinite(v):
                break
            if v0 is None:
                v0 = v
            if v < best:
                best, since, best_p = v, 0, {k: x.copy() for k, x in model.p.items()}
            else:
                since += EVERY
            if v < 0.9 * v0:
                armed = True
            if armed and since >= patience:
                break
        if t == steps:
            break
        model.setvec(opt.step(model.vec(model.p), model.vec(model.grad(Xtr, ytr, l1))))
    if best_p is not None:
        model.p = best_p
    return int(t), float(best)


def make_opt(arm, cfg, model):
    rv = model.rates_vec(GD_RATES)
    if arm == 'gd':
        return Opt('gd', cfg['lr'] * rv)
    if arm == 'adam':
        return Opt('adam', cfg['lr'])
    if arm == 'adam_big':
        return Opt('adam', cfg['ea'] * rv, eps_adam=cfg['ea'])
    if arm == 'shared_v':
        return Opt('adam_g', cfg['lr'], pre=rv if cfg['pre'] else None)
    if arm == 'sign':
        return Opt('sign', cfg['lr'])
    raise ValueError(arm)


def summarise(model, data, tg, info, arm, conf, key):
    """One flat JSON-able record of a trained model: counts and centres of the active units per edge (tau = TAU), parameters, validation MSE, clean test RMSE."""
    X = data['Xtr']
    act = model.activity(X)[2][0]
    cc = model.p['c']
    cnt = KS.edge_counts(act[None])
    return dict(kind='p1', tid=int(key[0]), noise=float(key[1]), seed=int(key[2]), arm=arm, conf=conf, steps=int(info[0]), val_mse=float(info[1]),
                test_rmse=float(np.sqrt(np.mean((model.fwd(data['Xte'])[0] - data['yte']) ** 2))), n_val=int(len(data['yva'])),
                counts=[int(x) for x in cnt], centres=[[round(float(c), 4) for c, a in zip(cc[i], act[i]) if a] for i in range(act.shape[0])],
                P=int(3 * int(act.sum()) + 1), truth=[int(x) for x in tg['K']], truth_c=[[float(x) for x in b[:, 1]] for b in tg['bumps']])


def control_runs(data, tg, key, seed):
    """The control: free weights, Adam (lr CTRL_LR) with an L1 penalty, magnitude pruning at TAU, a FT_STEPS fine-tune; one record per penalty."""
    recs = []
    for lam in LAMS:
        m = AddNet(KS.D_SYN, M_UNITS, 'free', seed=seed)
        i1 = fit(m, data, Opt('adam', CTRL_LR), STEPS, l1=lam)
        act = m.activity(data['Xtr'])[2][0]
        m.mask = m.mask * act.astype(float)
        i2 = fit(m, data, Opt('adam', CTRL_LR / 3.0), FT_STEPS)
        recs.append(summarise(m, data, tg, (i1[0] + i2[0], i2[1]), 'ctrl', json.dumps(dict(lam=lam), sort_keys=True), key))
    return recs


def p1_job(spec):
    """spec = (tid, noise, seed, {arm: [cfg, ...]}, with_ctrl): all configurations of all arms (and the control) on one (target, noise, seed)."""
    tid, noise, seed, arms, with_ctrl = spec
    t0 = time.time()
    tg = KS.make_target(tid)
    data = KS.sample_synth(tg, noise, seed)
    key = (tid, noise, seed)
    recs = []
    for arm in ARMS:
        for cfg in arms.get(arm, []):
            m = AddNet(KS.D_SYN, M_UNITS, 'fact', EPS0)
            t1 = time.time()
            info = fit(m, data, make_opt(arm, cfg, m), STEPS)
            r = summarise(m, data, tg, info, arm, json.dumps(cfg, sort_keys=True), key)
            r['secs'] = round(time.time() - t1, 2)
            recs.append(r)
    if with_ctrl:
        t1 = time.time()
        cr = control_runs(data, tg, key, seed)
        for r in cr:
            r['secs'] = round((time.time() - t1) / len(cr), 2)
        recs += cr
    return recs, time.time() - t0


# ------------------------------------------------------------------------------------------------------------------------------------
# P1: analysis (numpy only)
# ------------------------------------------------------------------------------------------------------------------------------------
def group_runs(recs, arms_needed=None):
    """{(tid, noise, seed): {arm: [records]}} of the P1 records."""
    runs = {}
    for r in recs:
        if r.get('kind') == 'p1':
            runs.setdefault((r['tid'], r['noise'], r['seed']), {}).setdefault(r['arm'], []).append(r)
    return runs


def merged_counts(rec, radius):
    return np.array([merged_count(cs, radius) for cs in rec['centres']])


def exact_vec(rec, radius):
    return merged_counts(rec, radius) == np.array(rec['truth'])


def share_pos(rec, radius):
    """Exact-edge share of one record on its edges with K > 0 (nan if there are none)."""
    t = np.array(rec['truth'])
    return float(np.mean(exact_vec(rec, radius)[t > 0])) if (t > 0).any() else float('nan')


def pick_ctrl(recs, radius=None):
    """The control record of one run: by validation (one-SE rule, spike module's own `one_se_pick`) if radius is None, else the oracle: the penalty whose merged counts at
    `radius` have the highest exact share on the K > 0 edges (ties: fewer active parameters, then the smaller penalty)."""
    if radius is None:
        cands = [dict(val=r['val_mse'], P=r['P'], r=r) for r in recs]
        return KS.one_se_pick(cands, recs[0]['n_val'])['r']
    return max(recs, key=lambda r: (share_pos(r, radius), -r['P'], -json.loads(r['conf'])['lam']))


def src(arm):
    """The key of a run's record dict that holds the arm's records (both controls come from 'ctrl')."""
    return 'ctrl' if arm in ('ctrl', 'ctrl_or') else arm


def select_arm_rec(runs_arm, arm, freeze=None, radius=None):
    """The one record of `arm` in a run: the frozen configuration (or the only one), the control by the rules of pick_ctrl."""
    if arm == 'ctrl':
        return pick_ctrl(runs_arm, None)
    if arm == 'ctrl_or':
        return pick_ctrl(runs_arm, MERGE_R if radius is None else radius)
    if freeze is not None:
        conf = json.dumps(freeze[arm], sort_keys=True)
        c = [r for r in runs_arm if r['conf'] == conf]
        assert len(c) == 1, (arm, conf, len(c))
        return c[0]
    assert len(runs_arm) == 1
    return runs_arm[0]


def complete_runs(runs, arms):
    need = list(arms)
    return {k: v for k, v in sorted(runs.items()) if all(a in v for a in need)}


def target_means(runs, arm, fn, freeze=None):
    """{tid: mean over the target's runs (noise, seed) of fn(record)} for the arm, in ascending tid."""
    by = {}
    for k, v in runs.items():
        by.setdefault(k[0], []).append(fn(select_arm_rec(v[src(arm)], arm, freeze)))
    return {t: float(np.mean(x)) for t, x in sorted(by.items())}


def arm_table(runs, arm, radius, freeze=None):
    """Per-target K>0 exact shares at the radius, and pooled statistics over all runs."""
    recs = [select_arm_rec(v[src(arm)], arm, freeze, radius) for k, v in runs.items()]
    sh = target_means(runs, arm, lambda r: share_pos(r, radius), freeze) if arm != 'ctrl_or' else _oracle_target_means(runs, radius, freeze)
    T = np.array([r['truth'] for r in recs])
    C = np.array([merged_counts(r, radius) for r in recs])
    pos, emp = T > 0, T == 0
    all_exact = float(np.mean((C == T).all(1)))
    return dict(share_pos=float(np.mean(list(sh.values()))), by_target=sh, share_all=float((C == T).mean()), share_empty=float((C == T)[emp].mean()),
                surplus_pos=float((np.maximum(C - T, 0) * pos).sum() / pos.sum()), missing_pos=float((np.maximum(T - C, 0) * pos).sum() / pos.sum()),
                all5=all_exact, n_runs=len(recs), rmse=float(np.mean([r['test_rmse'] for r in recs])), steps=float(np.mean([r['steps'] for r in recs])),
                cap=int(sum(r['steps'] >= STEPS for r in recs)), secs=float(np.mean([r.get('secs', float('nan')) for r in recs])),
                P=float(np.mean([r['P'] for r in recs])))


def _oracle_target_means(runs, radius, freeze=None):
    by = {}
    for k, v in runs.items():
        by.setdefault(k[0], []).append(share_pos(pick_ctrl(v['ctrl'], radius), radius))
    return {t: float(np.mean(x)) for t, x in sorted(by.items())}


def contrast(runs, a, b, radius, freeze=None):
    """Paired difference a - b of the per-target K>0 exact shares (cluster bootstrap over the targets): (mean, lo, hi, per-target diffs)."""
    sa = arm_table(runs, a, radius, freeze)['by_target']
    sb = arm_table(runs, b, radius, freeze)['by_target']
    d = [sa[t] - sb[t] for t in sa]
    m, lo, hi = KS.boot_ci(d, n_boot=N_BOOT, seed=BOOT_SEED)
    return m, lo, hi, d


def p1_verdict(runs, radius, freeze=None):
    """The pre-registered predictions (a), (b) and the secondary (d), at one merge radius."""
    out = dict(radius=radius, a=[], b=None, d=None)
    for hi_arm in ('gd', 'shared_v'):
        for lo_arm in ('adam', 'sign'):
            m, lo, hi, _ = contrast(runs, hi_arm, lo_arm, radius, freeze)
            out['a'].append(dict(hi=hi_arm, lo=lo_arm, mean=m, lo_ci=lo, hi_ci=hi, ok=bool(m >= GATE_GAP and lo > 0.0)))
    m, lo, hi, _ = contrast(runs, 'gd', 'shared_v', radius, freeze)
    out['b'] = dict(mean=m, lo_ci=lo, hi_ci=hi, ok=bool(abs(m) <= GATE_NEAR))
    m, lo, hi, _ = contrast(runs, 'adam_big', 'gd', radius, freeze)
    out['d'] = dict(mean=m, lo_ci=lo, hi_ci=hi, ok=bool(abs(m) <= GATE_NEAR))
    out['a_ok'] = all(x['ok'] for x in out['a'])
    out['b_ok'] = out['b']['ok']
    return out


def freeze_from_dev(recs):
    """Pick, per arm, the dev configuration with the highest mean K>0 exact share at the declared radius (ties: the higher unmerged share, then the grid order).
    Returns {arm: cfg} with GD's frozen configuration, and the table of scores."""
    runs = group_runs(recs)
    frz, table = dict(gd=dict(GD_CFG)), {}
    for arm in ARMS[1:]:
        scores = []
        for cfg in DEV_GRID[arm]:
            conf = json.dumps(cfg, sort_keys=True)
            rs = [[r for r in v[arm] if r['conf'] == conf][0] for k, v in sorted(runs.items()) if arm in v]
            scores.append((float(np.mean([share_pos(r, MERGE_R) for r in rs])), float(np.mean([share_pos(r, 0.0) for r in rs])), cfg, len(rs)))
        best = max(range(len(scores)), key=lambda i: (scores[i][0], scores[i][1], -i))
        frz[arm] = scores[best][2]
        table[arm] = scores
    return frz, table


# ------------------------------------------------------------------------------------------------------------------------------------
# P2 (i): scalar unit with a prescribed drive, and the repository's flow
# ------------------------------------------------------------------------------------------------------------------------------------
def readout(p):
    import exp_parametrization as PZ
    if p == 1:
        return PZ.make_readout('D1')
    return PZ.make_readout('D%d' % p)


def scalar_h(t):
    """The drive's dual shape: eta(t) = (1 - mu) * lam * h(t), h(t) = sin(3 pi t): |h| = 1 at t = 1/6, 1/2, 5/6 (the dual reaches the margin 3 times, with both signs)."""
    return np.sin(3 * np.pi * t)


def scalar_unit(p, eps, mu, rtol=1e-12):
    """Max |w| over t in [0, 1] of one unit of depth p (the repository's readout class, w = u^p - v^p, (u, v)(0) = (s, s), s^2 = eps) driven by g(t) = d eta/dt,
    eta(t) = (1 - mu) * lam * h(t), lam = readout.scale(eps) (the wall height): the dual stays at the relative margin mu of the wall.  p = 1 has no wall: eta = h (amplitude 1).
    Integrated with the unit's own ODE p' = lr * g * dw(p)."""
    from scipy.integrate import solve_ivp
    ro = readout(p)
    lam = 1.0 if p == 1 else float(ro.scale(eps))
    amp = (1.0 if p == 1 else (1.0 - mu) * lam)
    lr = ro.lr
    p0 = ro.init(1, eps)[0]

    def rhs(t, z):
        g = amp * 3 * np.pi * np.cos(3 * np.pi * t)
        return lr * g * ro.dw(z[None, :])[0]
    ts = np.unique(np.r_[np.linspace(0, 1, 2001), 1 / 6, 0.5, 5 / 6])
    atol = np.maximum(np.abs(p0), 1e-3 if p == 1 else 1e-300) * 1e-10
    sol = solve_ivp(rhs, (0, 1), p0, method='DOP853', rtol=rtol, atol=atol, t_eval=ts)
    if sol.status != 0:
        return dict(wmax=float('nan'), status='failed')
    W = np.array([float(ro.w(sol.y[:, k][None, :])[0]) for k in range(sol.y.shape[1])])
    return dict(wmax=float(np.max(np.abs(W))), status='ok')


def closed_form_wmax(p, eps, mu):
    """The exact sup of |w| for |eta| <= (1 - mu) lam: s^p [mu^-q - (2 - mu)^-q] (p >= 3, q = p/(p-2)); p = 2: 2 s^2 sinh((1 - mu) ln(1/s^2)); p = 1: 1 * lr * 1."""
    s2 = eps
    if p == 1:
        return float(PZ_LRC)
    if p == 2:
        return float(2 * s2 * np.sinh((1 - mu) * np.log(1 / s2)))
    q = p / (p - 2.0)
    return float(s2 ** (p / 2.0) * (mu ** (-q) - (2.0 - mu) ** (-q)))


def margin_factor(p, s2, mu):
    """The part of the closed-form sup that depends on the margin: p >= 3: mu^-q - (2 - mu)^-q (so that sup |w| = s^p times it); p = 2: 2 sinh((1 - mu) ln(1/s2)) (sup |w| = s2 times it)."""
    if p == 2:
        return float(2 * np.sinh((1.0 - mu) * np.log(1.0 / s2)))
    q = p / (p - 2.0)
    return float(mu ** (-q) - (2.0 - mu) ** (-q))


def x_of_w(p, w, s2):
    """The relative dual x = |eta|/lam of a unit with |w| = w (p >= 2) from the closed form w = W(eta): bisection in x on [0, 1). p = 2: asinh(w / (2 s2)) / ln(1/s2)."""
    w = abs(float(w))
    if p == 2:
        return float(np.arcsinh(w / (2 * s2)) / np.log(1 / s2))
    q = p / (p - 2.0)
    sp = s2 ** (p / 2.0)

    def F(x):
        return sp * (np.exp(-q * np.log1p(-x)) - np.exp(-q * np.log1p(x)))
    if w <= 0:
        return 0.0
    lo, hi = 0.0, 1.0 - 1e-15
    if F(hi) < w:
        return hi
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if F(mid) < w:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def slope(eps_list, y_list):
    """Least-squares slope of log10 y against log10 eps."""
    x, y = np.log10(np.asarray(eps_list, float)), np.log10(np.asarray(y_list, float))
    return float(np.polyfit(x, y, 1)[0])


def dev_target(seed):
    """A development target for (i): three features in [-0.7, 0.7] at least 0.35 apart, amplitudes +-(0.5 ... 1.0) (the design's kind), generated from `seed`."""
    rng = np.random.default_rng(seed)
    while True:
        cs = np.sort(rng.uniform(-0.7, 0.7, 3))
        if np.min(np.diff(cs)) >= 0.35:
            break
    As = rng.choice([-1.0, 1.0], 3) * rng.uniform(0.5, 1.0, 3)
    return cs, As


def net_target_ref(cs, As):
    import exp_parametrization as PZ
    L, K, c0 = PZ.setup()
    f = K.synth(np.asarray(cs, float)[:, None], np.asarray(As, float))
    ref = PZ.reference(L, K, c0, f)
    tau_end = 1.05 * max(ref['tau']) + 1.0 if ref['tau'] else PZ.TAU_B
    return ref, tau_end


def net_cell(cs, As, p, eps, rand, tau_end):
    """The repository's flow of the depth-p readout on the target (moving centres, rate Gm; unequal masses b_i = exp(rand N(0, 1)) when rand > 0): the largest |w| over
    time of every unit, the final |w|, the masses, the formed threshold."""
    import exp_parametrization as PZ
    from src.core import factorized_rbf as FR
    from src.core import rich_transit as RT
    L, K, c0 = PZ.setup()
    f = K.synth(np.asarray(cs, float)[:, None], np.asarray(As, float))
    ro = PZ.make_readout('D%d' % p)
    kw = dict(t_scale=RT.transit_scale(L, eps)) if p == 1 else {}
    d = len(c0)
    b = ro.masses(d, 1.0, np.random.default_rng(MASS_SEED), rand) if rand > 0 else np.ones(d)
    t0 = time.time()
    try:
        res = FR.train_flow(ro, K, f, c0, eps, tau_end, PZ.LR_C, max_secs=PZ.MAX_SECS, rng=np.random.default_rng(MASS_SEED) if rand > 0 else None, rand=rand, **kw)
    except Exception as e:                                   # a failed run is a recorded result, not a crash
        return dict(status=('error: ' + type(e).__name__)[:16], secs=time.time() - t0)
    W = res['W']
    return dict(status=str(res['status'])[:12], secs=float(time.time() - t0), wmax=np.abs(W).max(0).tolist(), wend=np.abs(W[-1]).tolist(), b=[float(x) for x in b],
                thr=float(FR.formed_threshold(np.asarray(As, float))))


def p2i_net_jobs(targets, rand_arms):
    """All network cells of (i): (target index, mass arm, depth p, eps) for p in 1..5."""
    return [(i, arm, p, eps) for i in range(len(targets)) for arm in rand_arms for p in (1, 2) + DEPTHS_I for eps in EPS_I]


def pick_i_targets(sample, refs):
    """The test targets of (i): the sample targets (200 of the readouts study) with a finite path runner-up margin (margin_never) of at least MARGIN_PICK, sorted by
    that margin (ties by index); M_NET of them at the positions round(linspace(0, n-1, M_NET))."""
    c = [(refs[i]['margin_never'], i) for i in range(len(sample)) if np.isfinite(refs[i]['margin_never']) and refs[i]['margin_never'] >= MARGIN_PICK]
    c.sort()
    pos = np.unique(np.round(np.linspace(0, len(c) - 1, M_NET)).astype(int))
    return [c[k][1] for k in pos]


def fit_exponents(cells, p, arm, margin_min=None):
    """Per target: the exponent of the largest dormant coefficient against eps for depth p in the mass arm. cells: list of cell records (kind 'p2n') of ONE target and arm
    and depth, over the eps grid. The dormant population N* = units not formed (|w_end| <= thr) at every valid eps; the cell's maximum is over N* of the unit's largest
    |w| over time. A cell with status != 'ok' is dropped; fewer than 3 valid cells -> no exponent. The realised margin of the maximising unit is 1 - x(w); the target is
    excluded when that margin is below MARGIN_MIN at any eps. Returns dict(slope, slope_corrected, eps, wmax, margin, status, ...)."""
    ok = sorted([c for c in cells if c['status'] == 'ok'], key=lambda c: c['eps'])
    out = dict(n_valid=len(ok), slope=float('nan'), slope_corr=float('nan'), eps=[c['eps'] for c in ok], wmax=[], margin=[], excluded=None, idx=[])
    if len(ok) < 3:
        out['excluded'] = 'fewer than 3 valid eps'
        return out
    formed_any = np.zeros(len(ok[0]['wend']), bool)
    for c in ok:
        formed_any |= np.array(c['wend']) > c['thr']
    pop = np.flatnonzero(~formed_any)
    if pop.size == 0:
        out['excluded'] = 'no dormant unit'
        return out
    corr = []
    for c in ok:
        wm = np.array(c['wmax'])[pop]
        j = int(pop[int(np.argmax(wm))])
        w = float(np.max(wm))
        s2 = c['eps'] * c['b'][j]
        if w <= 0.0:
            out['excluded'] = 'zero coefficient'
            return out
        if p == 1:
            out['wmax'].append(w)
            out['margin'].append(float('nan'))
            out['idx'].append(j)
            continue
        x = x_of_w(p, w, s2)
        mu = 1.0 - x
        out['wmax'].append(w)
        out['margin'].append(mu)
        out['idx'].append(j)
        corr.append(w / margin_factor(p, s2, mu))
    mm = MARGIN_MIN if margin_min is None else margin_min
    if p >= 2 and min(out['margin']) < mm:
        out['excluded'] = 'realised margin %.3f below %.2f' % (min(out['margin']), mm)
    if min(out['wmax']) <= 0:
        out['excluded'] = out['excluded'] or 'zero coefficient'
        return out
    out['slope'] = slope(out['eps'], out['wmax'])
    if p >= 2 and corr:
        out['slope_corr'] = slope(out['eps'], corr)    # w with the realised-margin factor divided out: s^p = (eps b_j)^(p/2) of the maximiser (p/2 up to the unit identity changing)
    return out


# ------------------------------------------------------------------------------------------------------------------------------------
# P2 (ii): formed-set agreement on the 200 targets
# ------------------------------------------------------------------------------------------------------------------------------------
def mu_theta(p, eps, thr):
    """The formation margin: the margin below which a dormant unit can reach the formed threshold thr: thr^(-(p-2)/p) eps^((p-2)/2) (b = 1); depth 2: ln(1/thr)/ln(1/eps)."""
    if p == 2:
        return float(np.log(1 / thr) / np.log(1 / eps))
    return float(thr ** (-(p - 2.0) / p) * eps ** ((p - 2.0) / 2.0))


def p2ii_cell(K, f, c0, cs, As, ref, tau_end, p, eps):
    import exp_parametrization as PZ
    a = PZ.cell(PZ.make_readout('D%d' % p), K, f, c0, cs, As, eps, ref, tau_end)
    return dict(kind='p2c', p=int(p), eps=float(eps), set_equal=bool(a['set_equal']), order_equal=bool(a['order_equal']), n_formed=int(a['n_formed']),
                exact=bool(a['exact']), status=str(a['status']), secs=float(a['secs']), tau_rel=float(a['tau_rel']) if np.isfinite(a['tau_rel']) else None)


# ------------------------------------------------------------------------------------------------------------------------------------
# drivers
# ------------------------------------------------------------------------------------------------------------------------------------
def jsonable(x):
    if isinstance(x, dict):
        return {str(k): jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [jsonable(v) for v in x]
    if isinstance(x, np.ndarray):
        return jsonable(x.tolist())
    if isinstance(x, np.bool_):
        return bool(x)
    if isinstance(x, np.integer):
        return int(x)
    if isinstance(x, np.floating):
        return float(x)
    return x


def append_jsonl(path, recs):
    with open(path, 'a') as f:
        for r in recs:
            f.write(json.dumps(jsonable(r), sort_keys=True) + '\n')


def read_jsonl(path):
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return [json.loads(line) for line in f if line.strip()]


def run_p1(stage, work, freeze=None):
    """Run (or resume, record by record) the P1 jobs of a stage; a record is (target, noise, seed, arm, configuration)."""
    tids, seeds = (DEV_TIDS, DEV_SEEDS) if stage == 'dev' else (TEST_TIDS, TEST_SEEDS)
    path = os.path.join(work, 'p1_%s.jsonl' % stage)
    have = read_jsonl(path)
    done = {(r['tid'], r['noise'], r['seed'], r['arm'], r['conf']) for r in have}
    if stage == 'dev':
        arms = {a: DEV_GRID[a] for a in DEV_GRID}
        arms['gd'] = [dict(GD_CFG)]
    else:
        arms = {a: [freeze[a]] for a in ARMS}
    for tid in tids:
        for noise in NOISES:
            for seed in seeds:
                todo = {a: [c for c in cs if (tid, noise, seed, a, json.dumps(c, sort_keys=True)) not in done] for a, cs in arms.items()}
                with_ctrl = not any(d[:3] == (tid, noise, seed) and d[3] == 'ctrl' for d in done)
                if not any(todo.values()) and not with_ctrl:
                    continue
                recs, secs = p1_job((tid, noise, seed, todo, with_ctrl))
                append_jsonl(path, recs)
                print('p1 %s tid %d noise %g seed %d: %d records, %.0f s' % (stage, tid, noise, seed, len(recs), secs), flush=True)


def p2_sample():
    import exp_parametrization as PZ
    return PZ.sample_targets()


def p2_dev_pool(n):
    """n development targets for (ii): rows of the certificate raw (families A, B) that are NOT in the 200-target sample (first n in file order)."""
    import exp_parametrization as PZ
    import exp_cohesion_leak as E
    raw = E.load(PZ.CERT_RAW)
    smp = p2_sample()
    keys = {(tuple(np.round(t['cs'], 9)), tuple(np.round(t['As'], 9))) for t in smp}
    out = []
    for fam in ('A', 'B'):
        for r in raw[fam]:
            if r.get('skip'):
                continue
            cs, As = np.asarray(r['cs'], float), np.asarray(r['As'], float)
            if (tuple(np.round(cs, 9)), tuple(np.round(As, 9))) in keys:
                continue
            out.append(dict(fam=fam, cs=cs, As=As))
            if len(out) == n:
                return out
    return out


def run_p2(stage, work, freeze=None):
    import exp_parametrization as PZ
    from src.core import factorized_rbf as FR
    path = os.path.join(work, 'p2_%s.jsonl' % stage)
    L, K, c0 = PZ.setup()
    done = {(r['kind'], r.get('tgt'), r.get('arm'), r.get('p'), r.get('eps')) for r in read_jsonl(path)}
    # ---- scalar unit
    for p in (1, 2) + DEPTHS_I:
        for mu in (MUS_SCALAR if p > 1 else (0.0,)):
            for eps in EPS_I:
                if ('p2s', None, repr(mu), p, eps) in done:
                    continue
                r = scalar_unit(p, eps, mu)
                append_jsonl(path, [dict(kind='p2s', arm=repr(mu), p=p, eps=eps, mu=mu, wmax=r['wmax'], status=r['status'], exact=closed_form_wmax(p, eps, mu))])
    # ---- the targets of (ii) (all 200 at test, 12 development targets at dev) and of (i)
    if stage == 'dev':
        sample = p2_dev_pool(12)
        tg_i = [dict(name='design', cs=list(PZ.CSTAR), As=list(PZ.ASTAR))] + [dict(name='dev%d' % s, cs=dev_target(s)[0].tolist(), As=dev_target(s)[1].tolist()) for s in (7001, 7002, 7003)]
    else:
        sample = p2_sample()
    refs = [PZ.reference(L, K, c0, K.synth(t['cs'][:, None], t['As'])) for t in sample]
    if stage == 'test':
        pick = pick_i_targets(sample, refs)
        tg_i = [dict(name='design', cs=list(PZ.CSTAR), As=list(PZ.ASTAR))] + [dict(name='s%d' % i, cs=sample[i]['cs'].tolist(), As=sample[i]['As'].tolist(), sample_index=int(i)) for i in pick]
    # ---- (i) the network
    for ti, t in enumerate(tg_i):
        ref, tau_end = net_target_ref(t['cs'], t['As'])
        if ('p2t', ti, None, None, None) not in done:
            append_jsonl(path, [dict(kind='p2t', tgt=ti, name=t['name'], cs=t['cs'], As=t['As'], sample_index=t.get('sample_index'), margin=ref['margin'], margin_never=ref['margin_never'],
                                     order=ref['order'], tau_end=tau_end)])
        for arm, rand in (('un', RAND_NET), ('eq', 0.0)):
            for p in (1, 2) + DEPTHS_I:
                for eps in EPS_I:
                    if ('p2n', ti, arm, p, eps) in done:
                        continue
                    c = net_cell(t['cs'], t['As'], p, eps, rand, tau_end)
                    append_jsonl(path, [dict(kind='p2n', tgt=ti, arm=arm, p=p, eps=eps, **c)])
        print('p2 %s (i) target %d/%d (%s) done' % (stage, ti + 1, len(tg_i), t['name']), flush=True)
    # ---- (ii)
    for i, t in enumerate(sample):
        ref = refs[i]
        f = K.synth(t['cs'][:, None], t['As'])
        tau_end = 1.05 * max(ref['tau']) + 1.0 if ref['tau'] else PZ.TAU_B
        if ('p2r', i, None, None, None) not in done:
            append_jsonl(path, [dict(kind='p2r', tgt=i, fam=t['fam'], cs=t['cs'], As=t['As'], margin=ref['margin'], margin_never=ref['margin_never'], thr=float(FR.formed_threshold(t['As'])),
                                     order=ref['order'], set=ref['set'], relaxed=ref['relaxed'])])
        for p in DEPTHS_II:
            for eps in EPS_II_ALL:
                if ('p2c', i, None, p, eps) in done:
                    continue
                c = p2ii_cell(K, f, c0, t['cs'], t['As'], ref, tau_end, p, eps)
                c['tgt'] = i
                append_jsonl(path, [c])
        if (i + 1) % 20 == 0:
            print('p2 %s (ii) %d/%d targets' % (stage, i + 1, len(sample)), flush=True)


def merge_raw(stage, work, raw):
    recs = read_jsonl(os.path.join(work, 'p1_%s.jsonl' % stage)) + read_jsonl(os.path.join(work, 'p2_%s.jsonl' % stage))
    meta = dict(stage=stage, n_records=len(recs), dev_tids=list(DEV_TIDS), test_tids=list(TEST_TIDS), dev_seeds=list(DEV_SEEDS), test_seeds=list(TEST_SEEDS), noises=list(NOISES))
    save_raw(raw, meta, recs)
    return meta, recs


def save_raw(path, meta, recs):
    np.savez(path, meta=np.array(json.dumps(meta, sort_keys=True)), records=np.array(json.dumps(recs, sort_keys=True)))


def load_raw(path):
    z = np.load(path, allow_pickle=False)
    return json.loads(str(z['meta'])), json.loads(str(z['records']))


# ------------------------------------------------------------------------------------------------------------------------------------
# report
# ------------------------------------------------------------------------------------------------------------------------------------
def f3(x, n=3):
    return 'nan' if x is None or (isinstance(x, float) and not np.isfinite(x)) else ('%.' + str(n) + 'f') % x


def fg(x, n=3):
    return 'nan' if x is None or (isinstance(x, float) and not np.isfinite(x)) else ('%.' + str(n) + 'g') % x


def ci(m, lo, hi, n=3):
    return '%s [%s, %s]' % (f3(m, n), f3(lo, n), f3(hi, n))


def yn(b):
    return 'yes' if b else 'no'


def table(head, rows):
    return ['| ' + ' | '.join(head) + ' |', '|' + '---|' * len(head)] + ['| ' + ' | '.join(str(c) for c in r) + ' |' for r in rows]


def p1_report(recs, meta, freeze):
    """Markdown lines of the P1 section and the verdict dict."""
    runs = complete_runs(group_runs(recs), [a for a in ARMS] + ['ctrl'])
    tids = sorted({k[0] for k in runs})
    L = []
    if not runs:
        return ['No complete P1 runs.'], None
    nr, nt = len(runs), len(tids)
    arms_all = list(ARMS) + ['ctrl', 'ctrl_or']
    tabs = {r: {a: arm_table(runs, a, r, freeze) for a in arms_all} for r in RADII}
    V = {r: p1_verdict(runs, r, freeze) for r in RADII}
    L += ['## P1. Optimizers on the additive Gaussian-unit KAN', '',
          'Size: %d targets (ids %d-%d) x %d noise levels x %d data seeds = %d runs per arm (stage %s). Every arm is trained on the same data with the same early stopping (validation MSE every %d steps, '
          'the best iterate kept, patience %d from the first 10%% drop, cap %d steps) from the same factorised small start (eps = %g, M = %d units per edge). The exact-edge share of an edge is 1 when the '
          'number of active units (root-mean-square contribution at least %g of the output sd) after merging equals the number K of true bumps. **Primary metric**: exact-edge share on edges with K > 0 '
          '(per target: mean over its runs; then the mean over targets; intervals are 95%% percentile intervals of a cluster bootstrap over the targets, %d resamples). Merge radius %g (declared in the plan), '
          'single linkage on the centres of the active units, the same for every method.' % (
              nt, min(tids), max(tids), len(NOISES), nr // (nt * len(NOISES)), nr, meta.get('stage'), EVERY, PATIENCE, STEPS, EPS0, M_UNITS, KS.TAU, N_BOOT, MERGE_R), '']
    L += ['Frozen configurations (chosen on development targets only, plan and `optdepth_dev.md`): ' + '; '.join('%s: %s' % (a, json.dumps(freeze[a], sort_keys=True)) for a in ARMS) + '; control: lr %g, penalties %s, one-SE choice.' % (
        CTRL_LR, ', '.join('%g' % x for x in LAMS)), '']
    # --- primary table
    rows = []
    for a in arms_all:
        t = tabs[MERGE_R][a]
        u = tabs[0.0][a]
        lo = KS.boot_ci(list(t['by_target'].values()), n_boot=N_BOOT, seed=BOOT_SEED)
        rows.append((ARM_LABEL[a], ci(*lo), f3(u['share_pos']), f3(t['share_all']), f3(u['share_all']), f3(t['share_empty']), f3(t['surplus_pos'], 2), f3(t['missing_pos'], 2),
                     '%d of %d' % (round(t['all5'] * nr), nr), f3(t['rmse']), '%.0f' % t['steps'], '%d' % t['cap']))
    L += ['### P1.1 Exact-edge share by arm', '',
          '"K>0, merged" is the primary metric (radius %g); the other columns are the same runs unmerged, on all edges, on the empty edges, and the surplus and missing units per K>0 edge after merging.' % MERGE_R, '']
    L += table(['arm', 'K>0 exact, merged %g [95%%]' % MERGE_R, 'K>0 exact, unmerged', 'all edges, merged', 'all edges, unmerged', 'K=0 edges, merged', 'surplus units / K>0 edge', 'missing units / K>0 edge',
                'runs exact on all 5 edges', 'clean test RMSE', 'mean steps', 'runs at the cap'], rows) + ['']
    # --- contrasts
    V0 = V[MERGE_R]
    L += ['### P1.2 The pre-registered predictions (radius %g)' % MERGE_R, '',
          '* (a) GD and shared-v each exceed Adam and sign descent by at least %.2f in the K>0 exact share, with the lower end of the paired cluster-bootstrap interval above 0.' % GATE_GAP, '']
    L += table(['contrast', 'difference [95% interval]', 'at least %.2f and interval above 0' % GATE_GAP],
               [('%s - %s' % (x['hi'], x['lo']), ci(x['mean'], x['lo_ci'], x['hi_ci']), yn(x['ok'])) for x in V0['a']]) + ['']
    L += ['* (b) GD and shared-v within %.2f of each other: GD - shared_v = %s; |difference| <= %.2f: %s.' % (GATE_NEAR, ci(V0['b']['mean'], V0['b']['lo_ci'], V0['b']['hi_ci']), GATE_NEAR, yn(V0['b']['ok'])), '']
    ctr, ctro = tabs[MERGE_R]['ctrl'], tabs[MERGE_R]['ctrl_or']
    L += ['* (c) no prediction for the control; reported: K>0 exact share, merged %g: %s (penalty by validation), %s (penalty by the truth); unmerged %s. GD minus control: %s.' % (
        MERGE_R, f3(ctr['share_pos']), f3(ctro['share_pos']), f3(tabs[0.0]['ctrl']['share_pos']), ci(*contrast(runs, 'gd', 'ctrl', MERGE_R, freeze)[:3])), '']
    L += ['* (d) secondary, not part of the verdict: Adam with a large eps_adam within %.2f of GD: adam_big - GD = %s; within: %s.' % (GATE_NEAR, ci(V0['d']['mean'], V0['d']['lo_ci'], V0['d']['hi_ci']), yn(V0['d']['ok'])), '']
    L += ['**Verdict P1 (primary, radius %g): (a) %s, (b) %s.**' % (MERGE_R, 'HELD' if V0['a_ok'] else 'FAILED', 'HELD' if V0['b_ok'] else 'FAILED'), '']
    # --- unmerged and sweep
    L += ['### P1.3 The same rule at the other radii (the unmerged row is the secondary reading; the others are POST HOC)', '']
    rows = []
    for r in RADII:
        v = V[r]
        t = tabs[r]
        rows.append(('%g%s' % (r, ' (primary)' if r == MERGE_R else (' (unmerged, secondary)' if r == 0.0 else ' (post hoc)')),
                     *(f3(t[a]['share_pos']) for a in ('gd', 'shared_v', 'adam_big', 'adam', 'sign', 'ctrl')),
                     f3(min(x['mean'] for x in v['a'])), f3(min(x['lo_ci'] for x in v['a'])), yn(v['a_ok']), f3(v['b']['mean']), yn(v['b_ok'])))
    L += table(['radius', 'GD', 'shared-v', 'Adam large eps', 'Adam', 'sign', 'control (valid.)', 'smallest (a) gap', 'smallest (a) lower end', '(a) held', 'GD - shared-v', '(b) held'], rows) + ['']
    # --- per-target
    L += ['### P1.4 Per target (K>0 exact share, merged %g, mean over the target\'s %d runs)' % (MERGE_R, nr // nt), '']
    rows = []
    for t in tids:
        K = KS.make_target(t)['K']
        rows.append(('%d' % t, ' '.join(str(int(x)) for x in K), *(f3(tabs[MERGE_R][a]['by_target'][t]) for a in ('gd', 'shared_v', 'adam_big', 'adam', 'sign', 'ctrl'))))
    L += table(['target', 'K per edge', 'GD', 'shared-v', 'Adam large eps', 'Adam', 'sign', 'control'], rows) + ['']
    return L, dict(V=V, tabs=tabs, nr=nr, nt=nt)


def p2_report(recs, meta):
    L = []
    s = [r for r in recs if r.get('kind') == 'p2s']
    nets = [r for r in recs if r.get('kind') == 'p2n']
    tg = {r['tgt']: r for r in recs if r.get('kind') == 'p2t'}
    cells = [r for r in recs if r.get('kind') == 'p2c']
    refs = {r['tgt']: r for r in recs if r.get('kind') == 'p2r'}
    out = dict(i_scalar=None, i_net=None, ii=None)
    L += ['## P2. Depth of the readout', '']
    # ---------------- (i) scalar
    if s:
        L += ['### P2.1 Scalar unit under a prescribed drive (the closed-form unit, integrated with the repository\'s readout classes)', '',
              'One unit of depth p (w = u^p - v^p, (u, v)(0) = (s, s), s^2 = eps) driven by g(t) = d eta/dt with eta(t) = (1 - mu) lam sin(3 pi t): the dual reaches the relative margin mu of the wall '
              'lam three times with both signs (p = 2: the wall is ln(1/eps)/4; p = 1 has no wall: eta = sin(3 pi t)). Integrated with DOP853 (rtol 1e-12). Prediction: the largest |w| scales as eps^(p/2) '
              'for p >= 3 (exact value s^p [mu^-q - (2-mu)^-q], q = p/(p-2)), as eps^mu for p = 2, and does not depend on eps for p = 1 (w = Gm eta). Exponent of the fit over eps = %s.' % ', '.join('%g' % e for e in EPS_I), '']
        rows, ok_all = [], True
        for p in (1, 2) + DEPTHS_I:
            for mu in (MUS_SCALAR if p > 1 else (0.0,)):
                rr = sorted([r for r in s if r['p'] == p and abs(r['mu'] - mu) < 1e-12 and r['status'] == 'ok'], key=lambda r: r['eps'])
                if len(rr) < len(EPS_I):
                    rows.append((p, mu if p > 1 else '-', 'failed', '', '', '', ''))
                    ok_all = False
                    continue
                sl = slope([r['eps'] for r in rr], [r['wmax'] for r in rr])
                sx = slope([r['eps'] for r in rr], [r['exact'] for r in rr])
                pred = 0.0 if p == 1 else (mu if p == 2 else p / 2.0)
                tol = TOL_SCALAR if p != 2 else 0.02
                ok = abs(sl - pred) <= tol
                ok_all &= ok
                rel = max(abs(r['wmax'] / r['exact'] - 1) for r in rr)
                rows.append((p, '%g' % mu if p > 1 else '-', f3(sl, 4), f3(pred, 3), f3(sx, 4), '%.1e' % rel, yn(ok)))
        L += table(['p', 'mu', 'fitted exponent', 'predicted', 'exponent of the closed form', 'max rel. deviation of w_max from the closed form', 'within %g (p = 2: 0.02)' % TOL_SCALAR], rows) + ['']
        L += ['**Verdict P2.1 (scalar unit): %s.**' % ('HELD' if ok_all else 'FAILED'), '']
        out['i_scalar'] = ok_all
    # ---------------- (i) network
    if nets:
        L += ['### P2.2 The largest dormant coefficient on the repository\'s flow (moving centres; unequal masses b_i = exp(%g N(0,1)) in the arm "unequal", equal masses in the arm "equal")' % RAND_NET, '',
              'Flow: `exp_parametrization` (21 kernels of width %g on [-1, 1], centre rate Gm = %g, plain gradient flow from the small start s^2 = eps b_i) on the design target and on %d targets of the readouts study\'s '
              '200-target sample whose sliding path has a runner-up margin of at least %g; horizon tau_end = 1.05 max(path entry times) + 1. For each target, depth and arm: N* = units that are never formed '
              '(|w| <= the formed threshold) at the end at any eps of the grid; the largest dormant coefficient is the maximum over N* and over time of |w|; the exponent is the least-squares slope of '
              'log10 of it against log10 eps over eps = %s (cells with a failed integration are dropped; fewer than 3 valid cells: no exponent). A target enters the fit only if the realised margin of the '
              'maximiser, 1 - |eta|/lam, with |eta|/lam inverted from the closed form w = W(eta), is at least %g at every eps. Prediction (the plan): mean exponent over targets within %g of p/2 for p = 3, 4, 5.' % (
                  PZ_SIG, PZ_LRC, len(tg) - 1, MARGIN_PICK, ', '.join('%g' % e for e in EPS_I), MARGIN_MIN, TOL_EXP), '']
        res = {}
        res_all = {}
        for arm in ('un', 'eq'):
            for p in (1, 2) + DEPTHS_I:
                ex = []
                for ti in sorted(tg):
                    cs = [r for r in nets if r['tgt'] == ti and r['arm'] == arm and r['p'] == p]
                    ex.append((ti, fit_exponents(cs, p, arm)))
                res[(arm, p)] = ex
                res_all[(arm, p)] = [(ti, fit_exponents([r for r in nets if r['tgt'] == ti and r['arm'] == arm and r['p'] == p], p, arm, margin_min=0.0)) for ti in sorted(tg)]
        rows = []
        verdict_arm = {}
        for arm, lab in (('un', 'unequal masses (primary)'), ('eq', 'equal masses (secondary)')):
            oks = []
            for p in (1, 2) + DEPTHS_I:
                ex = res[(arm, p)]
                use = [e for ti, e in ex if np.isfinite(e['slope']) and e['excluded'] is None]
                sl = [e['slope'] for e in use]
                alls = [e['slope'] for ti, e in res_all[(arm, p)] if np.isfinite(e['slope'])]
                if not sl:
                    rows.append((lab, p, 0, len(ex), '', '', '', '', f3(float(np.mean(alls))) if alls else '', 'no target' if p >= 3 else 'contrast'))
                    if p >= 3:
                        oks.append(None)
                    continue
                m, lo, hi = KS.boot_ci(sl, n_boot=N_BOOT_II, seed=BOOT_SEED) if len(sl) > 1 else (sl[0], sl[0], sl[0])
                corr = [e['slope_corr'] for e in use if np.isfinite(e['slope_corr'])]
                if p >= 3:
                    pred = p / 2.0
                    ok = abs(m - pred) <= TOL_EXP
                    oks.append(bool(ok) if len(use) >= MIN_TARGETS else None)
                    verd = ('yes' if ok else 'no') + ('' if len(use) >= MIN_TARGETS else ' (fewer than %d targets: inconclusive)' % MIN_TARGETS)
                    pcol = f3(pred)
                elif p == 2:
                    mu_real = float(np.mean([np.mean(e['margin']) for e in use]))
                    pcol = '%s (the realised margin)' % f3(mu_real)
                    verd = 'contrast'
                else:
                    pcol = 'not 0 (the horizon scales with eps, see text)'
                    verd = 'contrast'
                rows.append((lab, p, len(use), len(ex), ci(m, lo, hi), '%s - %s' % (f3(min(sl)), f3(max(sl))), pcol, f3(float(np.mean(corr))) if corr else '', f3(float(np.mean(alls))) if alls else '', verd))
            verdict_arm[arm] = ('INCONCLUSIVE' if any(o is None for o in oks) else ('HELD' if all(oks) and len(oks) == 3 else 'FAILED'))
        L += table(['mass arm', 'p', 'targets in the fit', 'targets', 'mean exponent [95%, bootstrap over targets]', 'range over targets', 'predicted', 'mean exponent after removing the realised-margin factor (identity check)',
                    'sensitivity: mean exponent over all targets with a valid fit, no margin exclusion', 'within %g of p/2' % TOL_EXP], rows) + ['']
        # realised margins and details for primary arm
        rows = []
        for p in DEPTHS_I:
            for ti, e in res[('un', p)]:
                rows.append((p, tg[ti]['name'], f3(tg[ti]['margin_never']), f3(e['slope']), f3(e['slope_corr']), ', '.join(f3(m, 3) for m in e['margin']), e['excluded'] or ''))
        L += ['Per target, unequal masses (path margin = the equal-mass sliding path\'s runner-up margin; the realised margins are those of the maximising unit at eps = %s):' % ', '.join('%g' % x for x in EPS_I), '']
        L += table(['p', 'target', 'path margin', 'exponent', 'exponent at the realised margin', 'realised margins', 'excluded'], rows) + ['']
        d1 = [np.median([max(c['wmax']) for c in nets if c['tgt'] == ti and c['arm'] == 'un' and c['p'] == 1 and c['status'] == 'ok' and abs(c['eps'] - 1e-8) < 1e-20]) for ti in sorted(tg)]
        d3 = [np.median([max(c['wmax']) for c in nets if c['tgt'] == ti and c['arm'] == 'un' and c['p'] == 3 and c['status'] == 'ok' and abs(c['eps'] - 1e-8) < 1e-20]) for ti in sorted(tg)]
        L += ['Depth 1 (contrast): the readout has no wall and no dormancy; its dormant coefficients are of order 1e-3 to 1e-1 at every eps, with an exponent that is not 0 only because the flow of D1 is integrated to the horizon '
              '`t_scale = transit_scale(eps)`, which grows as eps^(-1/2). Median over targets of the largest dormant coefficient at eps = 1e-8, unequal masses: depth 1 %.1e, depth 3 %.1e.' % (float(np.median(d1)), float(np.median(d3))), '']
        L += ['**Verdict P2.2 (network, unequal masses, exponent within %g of p/2 for p = 3, 4 and 5, at least %d targets per depth): %s.** Equal masses (secondary): %s.' % (TOL_EXP, MIN_TARGETS, verdict_arm['un'], verdict_arm['eq']), '']
        out['i_net'] = dict(un=verdict_arm['un'], eq=verdict_arm['eq'], res=res)
    # ---------------- (ii)
    if cells:
        n = len(refs)
        L += ['### P2.3 Formed-set agreement with the sliding path on the readouts study\'s %d benchmark targets' % n, '',
              'For each target the sliding path (`exp_parametrization.reference`) and the flow of the depth-p readout (D2 ... D5; plain gradient flow, moving centres, equal masses, eps = %s). Success: the '
              'formed set (|w| > min(0.01 max|A*|, 0.1 min|A*|) at the end) equals the path\'s final active set; a failed integration counts as a failure. Prediction (the plan): f2 < f3 <= f4 <= f5 at each '
              'eps and f_p(1e-8) >= f_p(1e-6) for each p.' % ', '.join('%g' % e for e in EPS_II), '']
        F, byt = {}, {}
        for p in DEPTHS_II:
            for e in EPS_II_ALL:
                v = np.zeros(n)
                for c in cells:
                    if c['p'] == p and abs(c['eps'] - e) < 1e-30 + 1e-12 * e:
                        v[c['tgt']] = 1.0 if (c['set_equal'] and c['status'].startswith('ok')) else 0.0
                byt[(p, e)] = v
                F[(p, e)] = float(v.mean())
        rows = []
        for p in DEPTHS_II:
            row = ['D%d' % p]
            for e in EPS_II:
                v = byt[(p, e)]
                m, lo, hi = KS.boot_ci(v.tolist(), n_boot=N_BOOT_II, seed=BOOT_SEED)
                row.append('%s [%s, %s] (%d/%d)' % (f3(m), f3(lo), f3(hi), int(v.sum()), n))
            rows.append(row)
        L += table(['readout'] + ['set equal, eps = %g' % e for e in EPS_II], rows) + ['']
        L += ['POST HOC (not part of the verdict): the same fractions among the cells whose integration succeeded (status ok), i.e. with the failed integrations removed from numerator and denominator:', '']
        rows = []
        for p in DEPTHS_II:
            row = ['D%d' % p]
            for e in EPS_II_ALL:
                cc = [c for c in cells if c['p'] == p and abs(c['eps'] - e) < 1e-12 * e and c['status'].startswith('ok')]
                row.append('%s (%d/%d)' % (f3(sum(1 for c in cc if c['set_equal']) / max(len(cc), 1)), sum(1 for c in cc if c['set_equal']), len(cc)))
            rows.append(row)
        L += table(['readout'] + ['eps = %g' % e for e in EPS_II_ALL], rows) + ['']
        L += ['Same, order of entry equal to the path\'s, and failed integrations:', '']
        rows = []
        for p in DEPTHS_II:
            row = ['D%d' % p]
            for e in EPS_II:
                cc = [c for c in cells if c['p'] == p and abs(c['eps'] - e) < 1e-12 * e]
                row.append('order equal %d/%d, status not ok %d' % (sum(1 for c in cc if c['order_equal'] and c['status'].startswith('ok')), n, sum(1 for c in cc if not c['status'].startswith('ok'))))
            rows.append(row)
        L += table(['readout'] + ['eps = %g' % e for e in EPS_II], rows) + ['']
        # paired differences
        L += ['Paired differences of the success indicator (cluster bootstrap over the targets):', '']
        rows = []
        checks = []
        for e in EPS_II:
            for pa, pb in ((2, 3), (3, 4), (4, 5)):
                d = (byt[(pb, e)] - byt[(pa, e)]).tolist()
                m, lo, hi = KS.boot_ci(d, n_boot=N_BOOT_II, seed=BOOT_SEED)
                rows.append(('eps = %g' % e, 'f%d - f%d' % (pb, pa), ci(m, lo, hi)))
            checks.append(F[(2, e)] < F[(3, e)] and F[(3, e)] <= F[(4, e)] and F[(4, e)] <= F[(5, e)])
        for p in DEPTHS_II:
            d = (byt[(p, EPS_II[1])] - byt[(p, EPS_II[0])]).tolist()
            m, lo, hi = KS.boot_ci(d, n_boot=N_BOOT_II, seed=BOOT_SEED)
            rows.append(('D%d' % p, 'f(%g) - f(%g)' % (EPS_II[1], EPS_II[0]), ci(m, lo, hi)))
            checks.append(F[(p, EPS_II[1])] >= F[(p, EPS_II[0])])
        L += table(['cell', 'difference', 'mean [95% interval]'], rows) + ['']
        ok_ii = all(checks)
        sec = [F[(2, e)] <= F[(3, e)] and F[(3, e)] <= F[(4, e)] and F[(4, e)] <= F[(5, e)] for e in EPS_II_SEC]
        mono = [all(F[(p, EPS_II_ALL[k + 1])] >= F[(p, EPS_II_ALL[k])] for k in range(len(EPS_II_ALL) - 1)) for p in DEPTHS_II]
        L += ['Secondary (amendment 1; where the depths are not yet saturated): success fractions at all four eps, with the non-strict ordering f2 <= f3 <= f4 <= f5 at eps = %s and monotonicity in eps (f at a smaller eps >= f at a larger one) for every depth.' % ', '.join('%g' % e for e in EPS_II_SEC), '']
        rows = []
        for p in DEPTHS_II:
            rows.append(['D%d' % p] + ['%s (%d/%d)' % (f3(F[(p, e)]), int(byt[(p, e)].sum()), n) for e in EPS_II_ALL])
        L += table(['readout'] + ['eps = %g' % e for e in EPS_II_ALL], rows) + ['']
        L += ['Secondary readings: ordering f2 <= f3 <= f4 <= f5 at eps = %s: %s; f non-decreasing as eps decreases, for D2 / D3 / D4 / D5: %s.' % (
            ', '.join('%g' % e for e in EPS_II_SEC), ', '.join(yn(x) for x in sec), ' / '.join(yn(x) for x in mono)), '']
        # prediction of the margin law, per target
        L += ['Per-target reading of the formation-margin law (secondary): predicted success iff margin_never(path) > mu_theta(p, eps) = thr^(-(p-2)/p) eps^((p-2)/2) (p = 2: ln(1/thr)/ln(1/eps)), thr the formed threshold of the target.', '']
        rows = []
        for p in DEPTHS_II:
            for e in EPS_II:
                pred = np.array([1.0 if (np.isnan(refs[i]['margin_never']) or refs[i]['margin_never'] > mu_theta(p, e, refs[i]['thr'])) else 0.0 for i in range(n)])
                v = byt[(p, e)]
                tp, fp, fn, tn = int(((pred == 1) & (v == 1)).sum()), int(((pred == 1) & (v == 0)).sum()), int(((pred == 0) & (v == 1)).sum()), int(((pred == 0) & (v == 0)).sum())
                rows.append(('D%d' % p, '%g' % e, f3(float(pred.mean())), f3(float(v.mean())), f3((tp + tn) / n), '%d / %d / %d / %d' % (tp, fp, fn, tn)))
        L += table(['readout', 'eps', 'predicted success rate', 'observed', 'accuracy of the per-target prediction', 'pred.+obs.+ / pred.+obs.- / pred.-obs.+ / pred.-obs.-'], rows) + ['']
        sat = all(F[(p, e)] >= 0.995 for p in (3, 4, 5) for e in EPS_II)
        L += ['**Verdict P2.3 (ordering f2 < f3 <= f4 <= f5 at both eps and f(1e-8) >= f(1e-6) for every depth): %s.**%s' % (
            'HELD' if ok_ii else 'FAILED', ' (The depth-3 to depth-5 part is saturated at 100%: the ordering there is a tie.)' if sat else ''), '']
        out['ii'] = dict(ok=ok_ii, sec_order=all(sec), sec_mono=all(mono), F={'%d_%g' % k: v for k, v in F.items()})
    return L, out


def verdict_table(p1, p2):
    rows = []
    if p1:
        V = p1['V'][MERGE_R]
        rows.append(('P1 (a)', 'GD and shared-v each exceed Adam and sign descent by >= %.2f in K>0 exact share, paired intervals above 0 (radius %g)' % (GATE_GAP, MERGE_R), 'HELD' if V['a_ok'] else 'FAILED'))
        rows.append(('P1 (b)', 'GD and shared-v within %.2f of each other' % GATE_NEAR, 'HELD' if V['b_ok'] else 'FAILED'))
        rows.append(('P1 (c)', 'no claim about the control; reported', 'reported'))
        rows.append(('P1 (d), secondary', 'Adam with large eps_adam within %.2f of GD' % GATE_NEAR, 'HELD' if V['d']['ok'] else 'FAILED'))
    if p2:
        if p2['i_scalar'] is not None:
            rows.append(('P2.1 scalar unit', 'exponent p/2 (p = 3, 4, 5) within %g; p = 2: mu; p = 1: 0' % TOL_SCALAR, 'HELD' if p2['i_scalar'] else 'FAILED'))
        if p2['i_net'] is not None:
            rows.append(('P2.2 network, unequal masses', 'exponent within %g of p/2 for p = 3, 4, 5' % TOL_EXP, p2['i_net']['un']))
            rows.append(('P2.2 network, equal masses, secondary', 'same rule', p2['i_net']['eq']))
        if p2['ii'] is not None:
            rows.append(('P2.3 formed-set agreement', 'f2 < f3 <= f4 <= f5 at both eps; f(1e-8) >= f(1e-6) for each depth', 'HELD' if p2['ii']['ok'] else 'FAILED'))
            rows.append(('P2.3 secondary (eps 1e-3, 1e-4; amendment 1)', 'f2 <= f3 <= f4 <= f5 at eps = 1e-3, 1e-4; f non-decreasing as eps decreases', 'HELD' if (p2['ii']['sec_order'] and p2['ii']['sec_mono']) else 'FAILED'))
    return table(['prediction', 'rule (declared in the plan)', 'result'], rows)


def report(meta, recs, freeze=None):
    L = ['# Optimizer and depth: two pre-registered prediction tests of the theory (%s)' % meta.get('stage'), '',
         'Generated by `experiments/exp_optdepth.py --from-raw` from `results/ml/optdepth%s_raw.npz` (%d records); plan and decision rules: `results/ml/optdepth_plan.md`. Everything below is computed '
         'from the stored records (no number is typed by hand).' % ('_dev' if meta.get('stage') == 'dev' else '', len(recs)), '']
    p1l, p1 = p1_report(recs, meta, freeze)
    p2l, p2 = p2_report(recs, meta)
    L += ['## Verdicts', ''] + verdict_table(p1, p2) + ['']
    L += p1l + p2l
    return '\n'.join(L) + '\n', p1, p2


def load_freeze(path):
    with open(path) as f:
        return json.load(f)['p1']


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--check-unused', action='store_true')
    ap.add_argument('--run', choices=('p1', 'p2'))
    ap.add_argument('--stage', choices=('dev', 'test'))
    ap.add_argument('--work')
    ap.add_argument('--freeze', default=os.path.join(REPO, 'results', 'ml', 'optdepth_freeze.json'))
    ap.add_argument('--merge', action='store_true')
    ap.add_argument('--raw')
    ap.add_argument('--freeze-from')
    ap.add_argument('--out')
    ap.add_argument('--from-raw', nargs='?', const='')
    ap.add_argument('--md')
    a = ap.parse_args()
    if a.check_unused:
        lines, ok = check_unused()
        print('\n'.join(lines))
        print('all unused: %s' % ok)
        sys.exit(0 if ok else 1)
    if a.run:
        os.makedirs(a.work, exist_ok=True)
        frz = load_freeze(a.freeze) if a.stage == 'test' else None
        (run_p1 if a.run == 'p1' else run_p2)(a.stage, a.work, frz)
        return
    if a.merge:
        merge_raw(a.stage, a.work, a.raw)
        return
    if a.freeze_from:
        meta, recs = load_raw(a.freeze_from)
        frz, table_ = freeze_from_dev(recs)
        with open(a.out, 'w') as f:
            json.dump(dict(p1=frz), f, indent=1, sort_keys=True)
            f.write('\n')
        return
    if a.from_raw is not None:
        raw = a.from_raw or os.path.join(OUT, 'optdepth_raw.npz')
        meta, recs = load_raw(raw)
        frz = load_freeze(a.freeze) if meta.get('stage') == 'test' else None
        if meta.get('stage') == 'dev':
            txt = dev_report(meta, recs)
        else:
            txt = report(meta, recs, frz)[0]
        md = a.md or os.path.join(OUT, 'optdepth.md')
        with open(md, 'w') as f:
            f.write(txt)


def dev_report(meta, recs):
    frz, table_ = freeze_from_dev(recs)
    L = ['# Optimizer and depth: development report (tuning on development targets only)', '',
         'Generated by `experiments/exp_optdepth.py --from-raw` from `results/ml/optdepth_dev_raw.npz` (%d records); plan: `results/ml/optdepth_plan.md`. Development targets %s, data seeds %s, noise levels %s '
         '(%d runs per configuration). Selection rule: per arm the configuration with the highest mean K>0 exact-edge share after merging at radius %g (ties: higher unmerged share, then grid order). Plain GD is '
         'the frozen configuration of `kan_freeze.json` and is not tuned.' % (len(recs), list(DEV_TIDS), list(DEV_SEEDS), list(NOISES), len(DEV_TIDS) * len(DEV_SEEDS) * len(NOISES), MERGE_R), '']
    L += ['## P1 tuning', '']
    for arm in ARMS[1:]:
        L += ['### %s' % ARM_LABEL[arm], '']
        L += table(['configuration', 'K>0 exact, merged %g' % MERGE_R, 'K>0 exact, unmerged', 'chosen'], [(json.dumps(s[2], sort_keys=True), f3(s[0]), f3(s[1]), 'yes' if s[2] == frz[arm] else '') for s in table_[arm]]) + ['']
    L += ['Chosen: ' + '; '.join('%s: %s' % (a, json.dumps(frz[a], sort_keys=True)) for a in ARMS), '']
    runs = complete_runs(group_runs(recs), ['gd', 'ctrl'] + [a for a in ARMS[1:]])
    if runs:
        tabs = {a: arm_table(runs, a, MERGE_R, frz) for a in list(ARMS) + ['ctrl']}
        un = {a: arm_table(runs, a, 0.0, frz) for a in list(ARMS) + ['ctrl']}
        L += ['## Development runs of the chosen configurations (%d runs; for the record, not part of any verdict)' % len(runs), '']
        L += table(['arm', 'K>0 exact, merged %g' % MERGE_R, 'K>0 exact, unmerged', 'runs exact on all 5 edges', 'mean steps', 'runs at the cap'],
                   [(ARM_LABEL[a], f3(tabs[a]['share_pos']), f3(un[a]['share_pos']), '%d of %d' % (round(tabs[a]['all5'] * len(runs)), len(runs)), '%.0f' % tabs[a]['steps'], tabs[a]['cap']) for a in list(ARMS) + ['ctrl']]) + ['']
    p2l, p2 = p2_report(recs, meta)
    L += ['## P2 development'] + [''] + p2l
    return '\n'.join(L) + '\n'


if __name__ == '__main__':
    main()
