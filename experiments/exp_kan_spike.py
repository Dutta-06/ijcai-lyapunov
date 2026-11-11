"""
Does the small-start factorised Gaussian-unit readout select the number of units in a Kolmogorov-Arnold network (KAN) with Gaussian-RBF edges?
=========================================================================================================================================

SPIKE with a kill rule (plan and GO/KILL rule: results/ml/kan_plan.md, committed before any test set was scored).  A KAN edge i -> j is a sum
of Gaussian units  phi_ij(x) = sum_k w_ijk exp(-(x - c_ijk)^2 / (2 s_ijk^2))  with learnable centres c and widths s.  Four trainings of such edges:

    fact   (ours)  w = nu^2 * A * tanh(alpha / A), alpha(0) = 0, nu(0)^2 = eps  (the paper's free mass-amplitude unit, effective depth 3);
                   trained by plain gradient descent (variant 'gd', the paper's flow) or Adam (variant 'adam'), early stopping on validation
    free   (control) w direct, standard init N(0, 0.1^2), same centres and widths, L1 penalty lambda on w, magnitude pruning and a fine-tune;
                   the penalty is chosen per run by validation (one-standard-error rule), the 'oracle' choice (by the truth) is reported too
    fixed  (FastKAN-style) centres on a fixed uniform grid of G points, width = spacing, only w trained (G = 2 ... 12)
    mlp    (matched parameters) two hidden layers of width h, SiLU

Active unit: its root-mean-square contribution on the training inputs is at least TAU * (standard deviation of its output node's pre-activation).
The same rule counts units of every KAN, and active parameters are 3 per active unit (w, c, s; 1 for fixed grids) plus one bias per node.

Benchmarks.  (a) SYNTH: additive targets sum_i g_i(x_i) of d = 5 inputs, g_i a sum of K_i in {0..4} Gaussian bumps (K known, some edges K = 0), two noise
levels, three seeds; the 1-D model-order claim is testable on [d, 1] KANs.  (b) REAL: tabular OpenML regression sets chosen by a fixed rule
(`select_real`), KANs [d, 4, 1] (primary) and [d, 1].  Development targets / sets are used for choices; the test is run once.

    python experiments/exp_kan_spike.py --stage dev|test --out DIR [--freeze FILE] [--workers 4]        (needs torch; run in a Kaggle CPU notebook)
    python experiments/exp_kan_spike.py --from-raw [RAW.npz] [--md OUT.md]            (numpy only; byte-identical regeneration; default KAN_OUT/kan_raw.npz -> kan.md)
"""
import argparse
import json
import os
import sys
import time

import numpy as np

# ------------------------------------------------------------------------------------------------------------------------------------
# registries (everything the plan fixes)
# ------------------------------------------------------------------------------------------------------------------------------------
D_SYN = 5                       # inputs of a synthetic target
X_LO, X_HI = -2.5, 2.5          # inputs ~ U(X_LO, X_HI); every grid / initial centres span the same interval
TEST_TARGETS = tuple(range(1, 11))
DEV_TARGETS = (101, 102, 103)
NOISES = (0.1, 0.4)             # noise sd as a fraction of the (unit) standard deviation of the target
SEEDS = (0, 1, 2)
DEV_SEEDS = (0, 1)
N_TRAIN, N_VAL, N_TEST_SYN = 1500, 500, 2000
TAU = 0.02                      # activity threshold (fraction of the output node's standard deviation)
TAUS = (0.01, 0.02, 0.05)       # sensitivity
M_UNITS = 10                    # initial units per edge of 'fact' and 'free'
LAMS = (0.0, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2)   # L1 penalties of the control (0: prune only)
GS = (2, 3, 4, 6, 8, 12)        # FastKAN-style grid sizes
HS = (2, 3, 4, 6, 8, 12, 16, 24)  # MLP widths
STEPS, PATIENCE, EVERY = 12000, 2000, 25
PAT_ADAM = 500                  # patience of the Adam-trained models on the real sets of the test stage (compute); everything else uses PATIENCE
REAL_LAMS = (0.0, 3e-3, 3e-2)   # penalties of the control on the real sets of the test stage (compute)
FT_STEPS = 600
AMP_CAP = 3.0                   # A of the amplitude tanh
HID = 4                         # hidden nodes of the two-layer real architectures
PRIMARY = 'H'                   # primary real architecture: ours hybrid [d, 4, 1] (first layer factorised from a small start, second layer standard), against the baselines on 'B'
BASE_ARCH = {'A': 'A', 'B': 'B', 'H': 'B'}   # the architecture whose baselines a model of the given architecture is compared with
REAL_CAP = 4000                 # rows kept of a real set (random subsample, fixed seed), split 60/20/20
NI_MARGIN = 1.05                # non-inferiority margin on test RMSE
GO_SHARE_GAIN = 0.10            # 'clearly higher': paired gain in exact-edge share at least this large (and the interval above 0)
GO_REAL_SETS = 6                # of 8
N_BOOT, BOOT_SEED = 10000, 12345
EXCLUDE_REAL = ('kin8nm', 'puma8nh', 'bank8fm', 'cpu_small', 'titanic')   # simulated (Delve), a subset of a listed set, a class-label task
DEV_RANKS = (3, 6, 9)           # ranks (1-based) among the first 11 qualifying sets that are development sets
OPENML_LIST = ('https://www.openml.org/api/v1/json/data/list/status/active/number_instances/1000..20000/number_features/6..31/'
               'number_missing_values/0/limit/1000/offset/%d')

# ------------------------------------------------------------------------------------------------------------------------------------
# synthetic targets (pure numpy)
# ------------------------------------------------------------------------------------------------------------------------------------
def make_target(tid):
    """Additive target of D_SYN inputs: per input K_i in 0..4 Gaussian bumps (centre, width, signed amplitude); at least one irrelevant input (K = 0),
    at least one rich edge (K >= 3), 6 <= sum K <= 11, centres of one edge at least 1.0 apart. Scaled to unit variance under U(X_LO, X_HI)^d."""
    rng = np.random.default_rng(1000 + int(tid))
    while True:
        K = rng.choice(5, size=D_SYN, p=[0.2, 0.25, 0.25, 0.15, 0.15])
        if K.min() == 0 and K.max() >= 3 and 6 <= K.sum() <= 11:
            break
    bumps = []
    for k in K:
        while True:
            cs = rng.uniform(-1.8, 1.8, size=k)
            if k < 2 or np.min(np.abs(cs[:, None] - cs[None]) + 9 * np.eye(k)) >= 1.0:
                break
        ss = rng.uniform(0.25, 0.45, size=k)
        aa = rng.choice([-1.0, 1.0], size=k) * rng.uniform(0.8, 1.5, size=k)
        bumps.append(np.stack([aa, cs, ss], 1) if k else np.zeros((0, 3)))
    t = dict(tid=int(tid), K=K.astype(int), bumps=bumps, scale=1.0)
    ref = np.random.default_rng(7).uniform(X_LO, X_HI, size=(200000, D_SYN))
    t['scale'] = float(np.std(target_clean(t, ref)))
    return t


def target_clean(t, X):
    y = np.zeros(X.shape[0])
    for i, b in enumerate(t['bumps']):
        for a, c, s in b:
            y += a * np.exp(-0.5 * ((X[:, i] - c) / s) ** 2)
    return y / t['scale']


def sample_synth(t, noise, seed):
    """Train / val (noisy) and test (clean) data of one run."""
    rng = np.random.default_rng(10 ** 6 + 1000 * int(t['tid']) + 10 * int(round(noise * 100)) + int(seed))
    X = rng.uniform(X_LO, X_HI, size=(N_TRAIN + N_VAL + N_TEST_SYN, D_SYN))
    f = target_clean(t, X)
    y = f + noise * rng.standard_normal(f.shape)
    a, b = N_TRAIN, N_TRAIN + N_VAL
    return dict(Xtr=X[:a], ytr=y[:a], Xva=X[a:b], yva=y[a:b], Xte=X[b:], yte=f[b:])


# ------------------------------------------------------------------------------------------------------------------------------------
# the activity rule and parameter counts (pure numpy)
# ------------------------------------------------------------------------------------------------------------------------------------
def unit_contrib(w, c, s, x):
    """Contribution of every unit of one layer on its inputs. w, c, s: (out, in, M); x: (n, in). Returns (n, out, in, M)."""
    u = (x[:, None, :, None] - c[None]) / s[None]
    return w[None] * np.exp(-0.5 * u * u)


def layer_activity(w, c, s, x, tau=TAU):
    """(rms, sd_out, active): rms (out, in, M) of each unit's contribution, sd_out (out,) standard deviation of the node's pre-activation sum,
    active = rms >= tau * sd_out (and rms > 0)."""
    U = unit_contrib(w, c, s, x)
    rms = np.sqrt(np.mean(U ** 2, axis=0))
    sd = np.std(U.sum(axis=(2, 3)), axis=0)
    act = (rms >= tau * sd[:, None, None]) & (rms > 0)
    return rms, sd, act


def active_params(n_active, n_nodes, per_unit=3):
    return int(per_unit) * int(n_active) + int(n_nodes)


def edge_counts(act):
    """Active units per edge of a single-output layer: act (1, in, M) -> (in,)."""
    return act[0].sum(axis=1).astype(int)


def mlp_params(d, h):
    return d * h + h + h * h + h + h + 1


# ------------------------------------------------------------------------------------------------------------------------------------
# statistics (pure numpy)
# ------------------------------------------------------------------------------------------------------------------------------------
def boot_ci(values_by_cluster, n_boot=N_BOOT, seed=BOOT_SEED, alpha=0.05):
    """Mean over clusters and a percentile interval of the cluster bootstrap (resample clusters with replacement)."""
    v = np.asarray(values_by_cluster, float)
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, v.size, size=(n_boot, v.size))
    m = v[idx].mean(axis=1)
    return float(v.mean()), float(np.quantile(m, alpha / 2)), float(np.quantile(m, 1 - alpha / 2))


def best_at_most(points, P):
    """Test RMSE of the best member of a family using at most P active parameters: points [(P_k, rmse_k)]. If none is that small, the smallest member
    is used (flag True)."""
    pts = [(p, r) for p, r in points if p <= P]
    if pts:
        return min(r for _, r in pts), False
    p0 = min(p for p, _ in points)
    return min(r for p, r in points if p == p0), True


def one_se_pick(cands, n_val):
    """cands: list of dicts with 'val' (validation MSE) and 'P' (active parameters): the candidate with the fewest active parameters among those with
    val <= best * (1 + sqrt(2 / n_val)) (the one-standard-error rule for a mean squared error); ties: the smaller val."""
    best = min(c['val'] for c in cands)
    ok = [c for c in cands if c['val'] <= best * (1.0 + np.sqrt(2.0 / n_val))]
    return min(ok, key=lambda c: (c['P'], c['val']))


def select_real(listing, loaded):
    """The rule choosing the real sets. listing: [(did, name, n, nfeat, nnum, nsym)] in any order; loaded(did) -> (X, y) arrays or None.
    Ascending id; keep numeric-only sets; first set of each lower-cased name; names in EXCLUDE_REAL or starting 'fri_' dropped; target needs >= 50 distinct
    finite values and nonzero spread, features finite; first 11 qualifying sets: ranks DEV_RANKS are development, the other 8 test."""
    seen, out = set(), []
    for did, name, n, nf, nn, ns in sorted(listing):
        low = name.lower()
        if ns != 0 or nn != nf or low in seen or low in EXCLUDE_REAL or low.startswith('fri_'):
            continue
        seen.add(low)
        xy = loaded(did)
        if xy is None:
            continue
        X, y = xy
        if len(np.unique(y)) >= 50 and np.isfinite(X).all() and np.isfinite(y).all() and np.std(y) > 0:
            out.append((did, name))
        if len(out) == 11:
            break
    dev = [out[r - 1] for r in DEV_RANKS if r <= len(out)]
    test = [o for i, o in enumerate(out) if (i + 1) not in DEV_RANKS]
    return dev, test


def split_real(n, seed):
    """Indices of train / val / test (60/20/20) of a real set of n rows after the subsample cap; the cap's subsample is part of the split seed."""
    rng = np.random.default_rng(777 + int(seed))
    idx = rng.permutation(n)[:REAL_CAP]
    a, b = int(0.6 * idx.size), int(0.8 * idx.size)
    return idx[:a], idx[a:b], idx[b:]


# ------------------------------------------------------------------------------------------------------------------------------------
# torch part (imported lazily: the report and the tests of the rules above need only numpy)
# ------------------------------------------------------------------------------------------------------------------------------------
def build_torch():
    import torch
    import torch.nn as nn

    class Layer(nn.Module):
        def __init__(self, din, dout, M, mode, eps=1e-3, rng=None):
            super().__init__()
            self.mode, self.din, self.dout, self.M = mode, din, dout, M
            if mode == 'fixed':
                grid = torch.linspace(X_LO, X_HI, M)
                sp = float(grid[1] - grid[0]) if M > 1 else (X_HI - X_LO)
                self.register_buffer('c', grid.expand(dout, din, M).clone())
                self.register_buffer('ls', torch.full((dout, din, M), float(np.log(sp))))
            else:
                grid = torch.linspace(X_LO, X_HI, M)
                sp = float(grid[1] - grid[0])
                self.c = nn.Parameter(grid.expand(dout, din, M).clone())
                self.ls = nn.Parameter(torch.full((dout, din, M), float(np.log(sp))))
            if mode == 'fact':
                self.alpha = nn.Parameter(torch.zeros(dout, din, M))
                self.nu = nn.Parameter(torch.full((dout, din, M), float(np.sqrt(eps))))
            else:
                self.w = nn.Parameter(0.1 * torch.randn(dout, din, M, generator=rng))
            self.register_buffer('mask', torch.ones(dout, din, M))
            self.bias = nn.Parameter(torch.zeros(dout))

        def weights(self):
            if self.mode == 'fact':
                return self.nu ** 2 * AMP_CAP * torch.tanh(self.alpha / AMP_CAP) * self.mask
            return self.w * self.mask

        def widths(self):
            return torch.exp(torch.clamp(self.ls, float(np.log(0.05)), float(np.log(3.0))))

        def forward(self, x):
            u = (x[:, None, :, None] - self.c[None]) / self.widths()[None]
            return (self.weights()[None] * torch.exp(-0.5 * u * u)).sum(dim=(2, 3)) + self.bias

    class KAN(nn.Module):
        def __init__(self, widths, M, mode, eps=1e-3, seed=0, hybrid=False):
            super().__init__()
            g = torch.Generator().manual_seed(int(seed))
            modes = [mode] + [('free' if hybrid else mode)] * (len(widths) - 2)
            self.layers = nn.ModuleList([Layer(widths[i], widths[i + 1], M, modes[i], eps, g) for i in range(len(widths) - 1)])
            for i in range(len(widths) - 2):
                self.register_buffer('nm%d' % i, torch.zeros(widths[i + 1]))
                self.register_buffer('ns%d' % i, torch.ones(widths[i + 1]))

        def forward(self, x, train=False):
            for i, L in enumerate(self.layers):
                x = L(x)
                if i < len(self.layers) - 1:             # detached standardisation of the hidden pre-activations (statistics of the training batch)
                    if train:
                        getattr(self, 'nm%d' % i).copy_(x.mean(0).detach())
                        getattr(self, 'ns%d' % i).copy_(x.std(0).detach().clamp_min(0.05))
                    x = (x - getattr(self, 'nm%d' % i)) / getattr(self, 'ns%d' % i)
            return x[:, 0]

        def inputs_by_layer(self, x):
            outs = []
            for i, L in enumerate(self.layers):
                outs.append(x.detach().cpu().numpy().astype(np.float64))
                x = L(x)
                if i < len(self.layers) - 1:
                    x = (x - getattr(self, 'nm%d' % i)) / getattr(self, 'ns%d' % i)
            return outs

    class MLP(nn.Module):
        def __init__(self, d, h, seed=0):
            super().__init__()
            torch.manual_seed(int(seed))
            self.net = nn.Sequential(nn.Linear(d, h), nn.SiLU(), nn.Linear(h, h), nn.SiLU(), nn.Linear(h, 1))

        def forward(self, x, train=False):
            return self.net(x)[:, 0]

    return torch, nn, KAN, MLP


def fit(model, data, kind, lr, steps=STEPS, l1=0.0, group_lr=None, patience=None):
    """Train on (Xtr, ytr) full batch with 0.5 * MSE (+ l1 * sum |w|), keep the iterate of lowest validation MSE (patience PATIENCE steps, counted from
    the first time the validation MSE is below 0.9 of its initial value). kind: 'adam' | 'gd'. Returns the model with the best state loaded and a dict."""
    torch = sys.modules['torch']
    Xtr, ytr, Xva, yva = (data[k] for k in ('Xtr', 'ytr', 'Xva', 'yva'))
    params = [p for p in model.parameters() if p.requires_grad]
    if kind == 'adam':
        opt = torch.optim.Adam(params, lr=lr)
    else:
        groups = []
        for nme, p in model.named_parameters():
            key = nme.split('.')[-1]
            groups.append(dict(params=[p], lr=group_lr.get(nme, group_lr.get(key, lr))))
        opt = torch.optim.SGD(groups, lr=lr)
    patience = PATIENCE if patience is None else patience
    best, best_state, v0, armed, since = np.inf, None, None, False, 0
    t0 = time.time()
    step = 0
    for step in range(steps + 1):
        if step % EVERY == 0:
            with torch.no_grad():
                v = float(torch.mean((model(Xva) - yva) ** 2))
            if not np.isfinite(v):
                break
            if v0 is None:
                v0 = v
            if v < best:
                best, since = v, 0
                best_state = {k: x.clone() for k, x in model.state_dict().items()}
            else:
                since += EVERY
            if v < 0.9 * v0:
                armed = True
            if armed and since >= patience:
                break
        if step == steps:
            break
        opt.zero_grad()
        r = model(Xtr, train=True) - ytr
        loss = 0.5 * torch.mean(r * r)
        if l1 > 0:
            loss = loss + l1 * sum(L.weights().abs().sum() for L in model.layers)
        loss.backward()
        opt.step()
    if best_state is not None:
        model.load_state_dict(best_state)
    return dict(val_mse=float(best), steps=int(step), secs=time.time() - t0)


def kan_activity(model, Xtr_t, tau):
    """Per-layer activity of a trained KAN on the training inputs (first 2000 rows): list of (rms, sd, act); active units and active parameters."""
    torch = sys.modules['torch']
    with torch.no_grad():
        xs = model.inputs_by_layer(Xtr_t[:2000])
        acts = []
        for L, x in zip(model.layers, xs):
            w = L.weights().detach().numpy().astype(np.float64)
            c = L.c.detach().numpy().astype(np.float64)
            s = L.widths().detach().numpy().astype(np.float64)
            acts.append((layer_activity(w, c, s, x, tau), (w, c, s)))
    return acts


def kan_summary(model, Xtr_t, per_unit):
    out = dict()
    for tau in TAUS:
        acts = kan_activity(model, Xtr_t, tau)
        na = sum(int(a[0][2].sum()) for a in acts)
        nn_ = sum(a[0][2].shape[0] for a in acts)
        out['P_%g' % tau] = active_params(na, nn_, per_unit)
        out['n_%g' % tau] = na
        if len(acts) == 1:
            out['counts_%g' % tau] = edge_counts(acts[0][0][2]).tolist()
            if tau == TAU:
                act, cc = acts[0][0][2][0], acts[0][1][1][0]
                out['centres_%g' % tau] = [[round(float(c), 4) for c, a in zip(cc[i], act[i]) if a] for i in range(act.shape[0])]
    return out


def prune_and_finetune(model, data, lr, Xtr_t, patience=None):
    """Magnitude pruning at TAU (units below the activity rule removed) then FT_STEPS of Adam without penalty (validation-best)."""
    torch = sys.modules['torch']
    acts = kan_activity(model, Xtr_t, TAU)
    for L, a in zip(model.layers, acts):
        L.mask.mul_(torch.tensor(a[0][2].astype(np.float32)))
    return fit(model, data, 'adam', lr / 3.0, steps=FT_STEPS, patience=patience)


def rmse(model, X, y):
    torch = sys.modules['torch']
    with torch.no_grad():
        return float(torch.sqrt(torch.mean((model(X) - y) ** 2)))


def to_t(d):
    torch = sys.modules['torch']
    return {k: torch.tensor(v, dtype=torch.float32) for k, v in d.items()}


def run_methods(data, d, arch, cfg, test_name='Xte', seed=0, grids=None):
    """All methods on one dataset (dict with numpy Xtr, ytr, Xva, yva, Xte, yte). Returns records. cfg: dict of hyper-parameters (see FREEZE / DEV_GRID)."""
    torch, nn, KAN, MLP = build_torch()
    torch.set_num_threads(1)
    T = to_t(data)
    widths = [d] + ([HID] if arch in ('B', 'H') else []) + [1]
    run_ours = arch in cfg.get('ours_archs', ('A', 'B', 'H'))
    run_base = arch != 'H'
    recs = []
    nva = len(data['yva'])
    pa = cfg.get('pat_adam', PATIENCE)

    def rec(method, conf, model, info, per_unit, mlp_h=None):
        r = dict(method=method, conf=conf, val_mse=info['val_mse'], steps=info['steps'], secs=info['secs'],
                 test_rmse=rmse(model, T['Xte'], T['yte']), n_val=nva)
        if mlp_h is not None:
            r['P_0.02'] = mlp_params(d, mlp_h)
            for tau in TAUS:
                r['P_%g' % tau] = r['P_0.02']
        else:
            r.update(kan_summary(model, T['Xtr'], per_unit))
        recs.append(r)

    for o in ((cfg['ours_by_arch'][arch] if 'ours_by_arch' in cfg else cfg['ours']) if run_ours else []):
        m = KAN(widths, M_UNITS, 'fact', eps=o['eps'], seed=seed, hybrid=(arch == 'H'))
        if o['kind'] == 'gd':
            gl = dict(alpha=o['lr'], c=o['lr'], ls=o['lr'], nu=o['lr'] * o['nu_ratio'], bias=o['lr'] * 0.1)
            for nme, _ in m.named_parameters():                               # later layers (hybrid): all their parameters at l2_ratio times the rate of the first layer
                if not nme.startswith('layers.0.'):
                    gl[nme] = o['lr'] * o.get('l2_ratio', 1.0) * (0.1 if nme.endswith('bias') else 1.0)
            info = fit(m, T, 'gd', o['lr'], steps=o.get('steps', STEPS), group_lr=gl)
        else:
            info = fit(m, T, 'adam', o['lr'], steps=o.get('steps', STEPS), patience=pa)
        rec('ours', json.dumps(o, sort_keys=True), m, info, 3)
    for lr in ((cfg['lr_ctrl'] if 'lr_ctrl' in cfg else cfg['lr_free']) if run_base else ()):
        for lam in cfg.get('lams', LAMS):
            m = KAN(widths, M_UNITS, 'free', seed=seed)
            info = fit(m, T, 'adam', lr, l1=lam, patience=pa)
            i2 = prune_and_finetune(m, T, lr, T['Xtr'], patience=pa)
            info = dict(val_mse=i2['val_mse'], steps=info['steps'] + i2['steps'], secs=info['secs'] + i2['secs'])
            rec('ctrl', json.dumps(dict(lr=lr, lam=lam), sort_keys=True), m, info, 3)
    for lr in ((cfg['lr_fixed'] if 'lr_fixed' in cfg else cfg['lr_free']) if run_base else ()):
        for G in cfg.get('gs', GS):
            m = KAN(widths, G, 'fixed', seed=seed)
            info = fit(m, T, 'adam', lr, patience=pa)
            rec('fixed', json.dumps(dict(lr=lr, G=G), sort_keys=True), m, info, 1)
    for lr in (cfg['lr_mlp'] if run_base else ()):
        for h in cfg.get('hs', HS):
            m = MLP(d, h, seed=seed)
            info = fit(m, T, 'adam', lr, patience=pa)
            rec('mlp', json.dumps(dict(lr=lr, h=h), sort_keys=True), m, info, 0, mlp_h=h)
    return recs


# ------------------------------------------------------------------------------------------------------------------------------------
# configurations
# ------------------------------------------------------------------------------------------------------------------------------------
DEV_REAL_OURS = [dict(kind='gd', lr=1.0, eps=1e-3, nu_ratio=0.05), dict(kind='gd', lr=0.5, eps=1e-3, nu_ratio=0.05), dict(kind='gd', lr=1.0, eps=1e-2, nu_ratio=0.05)]
DEV_GRID = dict(
    ours=[dict(kind='gd', lr=lr, eps=eps, nu_ratio=0.05) for lr in (0.5, 1.0, 2.0) for eps in (1e-2, 1e-3)] +
         [dict(kind='adam', lr=lr, eps=eps) for lr in (0.01, 0.03) for eps in (1e-2, 1e-4)],
    lr_free=(0.003, 0.01, 0.03), lr_mlp=(0.003, 0.01))


def load_real(did):
    import warnings
    warnings.filterwarnings('ignore')
    from sklearn.datasets import fetch_openml
    X, y = fetch_openml(data_id=int(did), return_X_y=True, as_frame=False, parser='auto')
    return np.asarray(X, float), np.asarray(y, float)


def real_listing():
    import urllib.request
    rows, off = [], 0
    while True:
        ds = json.loads(urllib.request.urlopen(OPENML_LIST % off, timeout=60).read())['data']['dataset']
        rows += ds
        if len(ds) < 1000:
            break
        off += 1000
    out = []
    for r in rows:
        q = {x['name']: x['value'] for x in r.get('quality', [])}
        try:
            out.append((int(r['did']), r['name'], *(int(float(q[k])) for k in ('NumberOfInstances', 'NumberOfFeatures', 'NumberOfNumericFeatures',
                                                                                 'NumberOfSymbolicFeatures'))))
        except Exception:
            continue
    return out


def real_data(did, seed):
    """Split data of a real set: quantile-normal inputs and standardised target, both fitted on the training rows."""
    from sklearn.preprocessing import QuantileTransformer
    X, y = load_real(did)
    tr, va, te = split_real(len(y), seed)
    qt = QuantileTransformer(n_quantiles=min(1000, len(tr)), output_distribution='normal', random_state=0).fit(X[tr])
    mu, sd = y[tr].mean(), y[tr].std()
    f = lambda i: (np.clip(qt.transform(X[i]), X_LO * 1.6, X_HI * 1.6), (y[i] - mu) / sd)
    d = {}
    for k, i in (('tr', tr), ('va', va), ('te', te)):
        d['X' + k], d['y' + k] = f(i)
    return d


def job(spec):
    """One (dataset, seed) job: spec = (kind, key..., cfg)."""
    kind = spec[0]
    t0 = time.time()
    if kind == 'synth':
        _, tid, noise, seed, arch, cfg = spec
        t = make_target(tid)
        data = sample_synth(t, noise, seed)
        recs = run_methods(data, D_SYN, 'A', cfg, seed=seed)
        for r in recs:
            r.update(kind='synth', tid=tid, noise=noise, seed=seed, arch='A', truth=t['K'].tolist(), truth_c=[b[:, 1].tolist() for b in t['bumps']])
    else:
        _, did, name, seed, arch, cfg = spec
        if 'lams_real' in cfg:
            cfg = dict(cfg, lams=cfg['lams_real'], pat_adam=PAT_ADAM)
        data = real_data(did, seed)
        recs = run_methods(data, data['Xtr'].shape[1], arch, cfg, seed=seed)
        for r in recs:
            r.update(kind='real', did=did, name=name, seed=seed, arch=arch, d=int(data['Xtr'].shape[1]), n_train=int(len(data['ytr'])))
    return recs, time.time() - t0


def save_raw(path, meta, recs):
    np.savez(path, meta=np.array(json.dumps(meta, sort_keys=True)), records=np.array(json.dumps(recs, sort_keys=True)))


def load_raw(path):
    z = np.load(path, allow_pickle=False)
    return json.loads(str(z['meta'])), json.loads(str(z['records']))


# ------------------------------------------------------------------------------------------------------------------------------------
# report (numpy only)
# ------------------------------------------------------------------------------------------------------------------------------------
def pick_val(recs, method, key='conf'):
    """The configuration of `method` with the lowest validation MSE per (run) among the given records, ties by order."""
    c = [r for r in recs if r['method'] == method]
    return min(c, key=lambda r: r['val_mse']) if c else None


def run_key(r):
    return (r['tid'], r['noise'], r['seed']) if r['kind'] == 'synth' else (r['did'], r['seed'], r['arch'])


def exact_share(r, tau=TAU):
    c = np.array(r['counts_%g' % tau])
    return float(np.mean(c == np.array(r['truth'])))


def synth_summary(recs, tau=TAU):
    """Per-run exact-edge shares of ours / control (one-SE by validation) / control (oracle by truth) / fixed grid (by validation), and count errors."""
    runs = {}
    for r in recs:
        if r['kind'] == 'synth':
            runs.setdefault(run_key(r), []).append(r)
    out = []
    for k in sorted(runs):
        rs = runs[k]
        o = [r for r in rs if r['method'] == 'ours']
        ctrl = [r for r in rs if r['method'] == 'ctrl']
        if not o or not ctrl or not any(r['method'] == 'fixed' for r in rs):
            continue                                            # an incomplete run is not scored
        o = o[0]
        cands = [dict(val=r['val_mse'], P=r['P_%g' % tau], r=r) for r in ctrl]
        cv = one_se_pick(cands, ctrl[0]['n_val'])['r']
        co = max(ctrl, key=lambda r: (exact_share(r, tau), -r['P_%g' % tau]))
        fx = pick_val(rs, 'fixed')
        row = dict(key=k, truth=o['truth'])
        for nm, r in (('ours', o), ('ctrl', cv), ('ctrl_or', co), ('fixed', fx)):
            c = np.array(r['counts_%g' % tau])
            row[nm] = dict(share=exact_share(r, tau), over=int(np.sum(np.maximum(c - np.array(o['truth']), 0))),
                           under=int(np.sum(np.maximum(np.array(o['truth']) - c, 0))), P=r['P_%g' % tau], rmse=r['test_rmse'],
                           secs=r['secs'], conf=r['conf'], match=centre_match(r))
        out.append(row)
    return out


def centre_match(r, tol=0.2):
    """(recall, precision) of the active units' centres against the true bump centres of a 1-layer run: the share of true bumps with an active unit within tol, and
    the share of active units within tol of a true bump. None when the record has no centres."""
    if 'centres_%g' % TAU not in r or 'truth_c' not in r:
        return None
    hit = tot_t = act = near = 0
    for found, true in zip(r['centres_%g' % TAU], r['truth_c']):
        for c in true:
            tot_t += 1
            hit += any(abs(f - c) <= tol for f in found)
        for f in found:
            act += 1
            near += any(abs(f - c) <= tol for c in true)
    return (hit / tot_t if tot_t else 1.0, near / act if act else 1.0)


def paired_gain(rows, a, b):
    """Gain in exact-edge share a - b, cluster bootstrap over targets (each target's gain averaged over noise levels and seeds)."""
    by = {}
    for r in rows:
        by.setdefault(r['key'][0], []).append(r[a]['share'] - r[b]['share'])
    return boot_ci([np.mean(v) for _, v in sorted(by.items())])


def real_summary(recs, arch=PRIMARY):
    """Per real set: ours (architecture `arch`) vs the best baseline of BASE_ARCH[arch] using at most as many active parameters (per split seed, then averaged on the
    log scale)."""
    runs = {}
    for r in recs:
        if r['kind'] == 'real' and (r['arch'] == arch and r['method'] == 'ours' or r['arch'] == BASE_ARCH[arch] and r['method'] != 'ours'):
            runs.setdefault((r['did'], r['name']), {}).setdefault(r['seed'], []).append(r)
    out = []
    for (did, name), seeds in sorted(runs.items()):
        lr_, rows = [], []
        for s, rs in sorted(seeds.items()):
            o = [r for r in rs if r['method'] == 'ours']
            if not o:
                continue
            o = o[0]
            P = o['P_%g' % TAU]
            best, who, flag = np.inf, '', False
            for fam in ('fixed', 'ctrl', 'mlp'):
                cand = [r for r in rs if r['method'] == fam]
                if not cand:
                    continue
                # a family's members: its (lr) per member chosen by validation within each capacity setting
                caps = {}
                for r in cand:
                    cj = json.loads(r['conf'])
                    capkey = tuple(sorted((k, v) for k, v in cj.items() if k != 'lr'))
                    if capkey not in caps or r['val_mse'] < caps[capkey]['val_mse']:
                        caps[capkey] = r
                pts = [(r['P_%g' % TAU], r['test_rmse']) for r in caps.values()]
                v, fl = best_at_most(pts, P)
                if v < best:
                    best, who, flag = v, fam, fl
            lr_.append(np.log(o['test_rmse'] / best))
            rows.append(dict(seed=s, ours=o['test_rmse'], P=P, base=best, who=who, clamped=flag, secs=o['secs']))
        m = float(np.mean(lr_))
        out.append(dict(did=did, name=name, logratio=m, ratio=float(np.exp(m)), ni=bool(m <= np.log(NI_MARGIN)), rows=rows,
                        ours_mean=float(np.mean([x['ours'] for x in rows])), P_mean=float(np.mean([x['P'] for x in rows])),
                        base_mean=float(np.mean([x['base'] for x in rows])), P_cv=float(np.std([x['P'] for x in rows]) / max(np.mean([x['P'] for x in rows]), 1)),
                        ours_sd=float(np.std([x['ours'] for x in rows])), who=sorted(set(x['who'] for x in rows)),
                        clamped=int(sum(x['clamped'] for x in rows))))
    return out


def verdict(synth_rows, real_rows):
    g = paired_gain(synth_rows, 'ours', 'ctrl')
    n_ni = sum(1 for r in real_rows if r['ni'])
    go_syn = (g[1] > 0.0) and (g[0] >= GO_SHARE_GAIN)
    go_real = n_ni >= GO_REAL_SETS
    return dict(gain=g, n_ni=n_ni, n_sets=len(real_rows), go_syn=bool(go_syn), go_real=bool(go_real), go=bool(go_syn and go_real))


def fmt(x, n=3):
    return ('%.' + str(n) + 'f') % x


def report(meta, recs):
    L = ['# Gaussian-unit KAN with the factorised small-start readout: spike report', '',
         'Generated by `experiments/exp_kan_spike.py --report --from-raw`; plan and rule: `results/ml/kan_plan.md`. Stage: %s. Runs: %d records.' % (meta.get('stage'), len(recs)), '']
    rows = synth_summary(recs)
    L += ['## Synthetic additive targets: active units per edge against the truth', '',
          'Runs: %d (targets x noise x seeds); an edge is exact when the active-unit count equals the number of true bumps (tau = %g).' % (len(rows), TAU), '',
          '| method | exact-edge share | over-counted units / run | under-counted units / run | active params | clean-test RMSE | train+eval s |', '|---|---|---|---|---|---|---|']
    for nm, lab in (('ours', 'ours (factorised, small start)'), ('ctrl', 'control: L1 + prune, penalty by validation'), ('ctrl_or', 'control: penalty by the truth (oracle)'),
                    ('fixed', 'fixed grid, G by validation')):
        sh = np.mean([r[nm]['share'] for r in rows])
        L.append('| %s | %s | %s | %s | %s | %s | %s |' % (lab, fmt(sh), fmt(np.mean([r[nm]['over'] for r in rows]), 2), fmt(np.mean([r[nm]['under'] for r in rows]), 2),
                                                           fmt(np.mean([r[nm]['P'] for r in rows]), 1), fmt(np.mean([r[nm]['rmse'] for r in rows])),
                                                           fmt(np.mean([r[nm]['secs'] for r in rows]), 1)))
    L.append('')
    for a, b, lab in (('ours', 'ctrl', 'ours - control (validation penalty)'), ('ours', 'ctrl_or', 'ours - control (oracle penalty)'), ('ours', 'fixed', 'ours - fixed grid')):
        m, lo, hi = paired_gain(rows, a, b)
        L.append('* paired gain in exact-edge share, %s: %s [%s, %s] (cluster bootstrap over targets, %d resamples)' % (lab, fmt(m), fmt(lo), fmt(hi), N_BOOT))
    for tau in TAUS:
        rt = synth_summary(recs, tau)
        L.append('* threshold %g: exact share ours %s, control %s, control oracle %s' % (tau, fmt(np.mean([r['ours']['share'] for r in rt])),
                                                                                        fmt(np.mean([r['ctrl']['share'] for r in rt])), fmt(np.mean([r['ctrl_or']['share'] for r in rt]))))
    if all(r[nm]['match'] is not None for r in rows for nm in ('ours', 'ctrl')):
        L += ['', 'Do the active units sit on the true bumps? Share of true bumps with an active unit within 0.2 (recall) and share of active units within 0.2 of a true bump (precision):', '',
              '| method | recall | precision |', '|---|---|---|']
        for nm, lab in (('ours', 'ours'), ('ctrl', 'control (validation penalty)'), ('ctrl_or', 'control (oracle penalty)')):
            L.append('| %s | %s | %s |' % (lab, fmt(np.mean([r[nm]['match'][0] for r in rows])), fmt(np.mean([r[nm]['match'][1] for r in rows]))))
    L += ['', '### Per target (exact-edge share, mean over noise levels and seeds)', '', '| target | true units per edge | ours | control | control oracle | fixed grid |', '|---|---|---|---|---|---|']
    for tid in sorted(set(r['key'][0] for r in rows)):
        rr = [r for r in rows if r['key'][0] == tid]
        L.append('| %s | %s | %s | %s | %s | %s |' % (tid, rr[0]['truth'], *(fmt(np.mean([r[m]['share'] for r in rr])) for m in ('ours', 'ctrl', 'ctrl_or', 'fixed'))))
    L += ['', '### Per noise level', '', '| noise sd | ours | control | control oracle |', '|---|---|---|---|']
    for nz in sorted(set(r['key'][1] for r in rows)):
        rr = [r for r in rows if r['key'][1] == nz]
        L.append('| %g | %s | %s | %s |' % (nz, *(fmt(np.mean([r[m]['share'] for r in rr])) for m in ('ours', 'ctrl', 'ctrl_or'))))
    L += ['', '### Stability across seeds (sd of the exact-edge share of ours over seeds, mean over target x noise cells)', '']
    cells = {}
    for r in rows:
        cells.setdefault(r['key'][:2], []).append(r['ours']['share'])
    sds = [np.std(v) for v in cells.values()]
    cc = {}
    for r in rows:
        cc.setdefault(r['key'][:2], []).append(r['ctrl']['share'])
    L.append('ours %s, control %s' % (fmt(np.mean(sds)), fmt(np.mean([np.std(v) for v in cc.values()]))))
    real = real_summary(recs, PRIMARY)
    realA = real_summary(recs, 'A')
    L += ['', '## Real tabular sets: test RMSE at matched active parameters', '',
          'Ours: KAN [d, %d, 1] with the first layer factorised from a small start and the second layer standard (hybrid); baselines: KAN [d, %d, 1] of the three families. Ours is one model per split; the baselines are families (fixed grids, control over penalties, MLP widths) and the comparison is with the best baseline '
          'using at most as many active parameters as ours (a family whose smallest member is larger is compared with that member: `clamped`). Non-inferior: mean log ratio <= log(%g).' % (HID, HID, NI_MARGIN), '',
          '| set (OpenML id) | d | active params (ours) | RMSE ours | RMSE best baseline | ratio | best family | sd RMSE over splits | CV of active params | clamped | non-inferior |',
          '|---|---|---|---|---|---|---|---|---|---|---|']
    dd = {(r['did']): r['d'] for r in recs if r['kind'] == 'real'}
    for r in real:
        L.append('| %s (%d) | %d | %s | %s | %s | %s | %s | %s | %s | %d | %s |' % (r['name'], r['did'], dd[r['did']], fmt(r['P_mean'], 1), fmt(r['ours_mean']), fmt(r['base_mean']),
                                                                                  fmt(r['ratio']), '/'.join(r['who']), fmt(r['ours_sd']), fmt(r['P_cv'], 2), r['clamped'],
                                                                                  'yes' if r['ni'] else 'no'))
    if realA:
        L += ['', 'KAN [d, 1] (additive), same comparison (secondary):', '', '| set | active params (ours) | RMSE ours | RMSE best baseline | ratio | non-inferior |', '|---|---|---|---|---|---|']
        for r in realA:
            L.append('| %s | %s | %s | %s | %s | %s |' % (r['name'], fmt(r['P_mean'], 1), fmt(r['ours_mean']), fmt(r['base_mean']), fmt(r['ratio']), 'yes' if r['ni'] else 'no'))
    tm = {}
    for r in recs:
        if r['kind'] == 'real' and r['arch'] in ('B', 'H'):
            tm.setdefault(r['method'], []).append(r['secs'])
    L += ['', 'Mean training time per model, real sets, [d, %d, 1] (s): %s' % (HID, ', '.join('%s %s' % (k, fmt(np.mean(v), 1)) for k, v in sorted(tm.items()))), '']
    pure = [r for r in recs if r['kind'] == 'real' and r['arch'] == 'B' and r['method'] == 'ours']
    if pure:
        L += ['Pure small start in both layers of [d, %d, 1] (ours on architecture B), test RMSE (standardised target) and active parameters, mean over splits:' % HID, '',
              '| set | RMSE | active params |', '|---|---|---|']
        for nm in sorted(set(r['name'] for r in pure)):
            rr = [r for r in pure if r['name'] == nm]
            L.append('| %s | %s | %s |' % (nm, fmt(np.mean([r['test_rmse'] for r in rr])), fmt(np.mean([r['P_0.02'] for r in rr]), 1)))
        L.append('')
    v = verdict(rows, real)
    L += ['## Verdict by the plan\'s rule', '',
          '* synthetic: paired gain over the validation-penalised control %s [%s, %s]; needs lower bound > 0 and gain >= %g: %s' % (fmt(v['gain'][0]), fmt(v['gain'][1]), fmt(v['gain'][2]), GO_SHARE_GAIN,
                                                                                                             'met' if v['go_syn'] else 'NOT met'),
          '* real: non-inferior on %d of %d sets; needs >= %d: %s' % (v['n_ni'], v['n_sets'], GO_REAL_SETS, 'met' if v['go_real'] else 'NOT met'),
          '* **%s**' % ('GO' if v['go'] else 'KILL'), '']
    if meta.get('stage') == 'dev':
        L += dev_section(recs)
    return '\n'.join(L)


# ------------------------------------------------------------------------------------------------------------------------------------
# stages
# ------------------------------------------------------------------------------------------------------------------------------------
def choose_freeze(recs):
    """Development rule (plan amendment 4). ours_synth (also ours on [d, 1] real sets) = the configuration with the highest mean exact-edge share on the development
    targets; ours (hybrid, [d, 4, 1]) = the configuration with the lowest mean log-ratio of development test RMSE to the best member of the baselines of B on the
    development real sets. Baselines: per kind of data (synthetic, real) and family the learning rate with the lowest mean log validation MSE over the runs."""
    ours = {}
    for r in recs:
        if r['method'] == 'ours':
            ours.setdefault(r['conf'], []).append(r)
    synth, real = {}, {}
    for conf, rs in ours.items():
        sy = [exact_share(r) for r in rs if r['kind'] == 'synth']
        if sy:
            synth[conf] = float(np.mean(sy))
        ratios = []
        for r in rs:
            if r['kind'] == 'real' and r['arch'] == PRIMARY:
                same = [x for x in recs if x['kind'] == 'real' and x['did'] == r['did'] and x['seed'] == r['seed'] and x['arch'] == BASE_ARCH[PRIMARY] and x['method'] in ('fixed', 'ctrl', 'mlp')]
                ratios.append(np.log(r['test_rmse'] / min(x['test_rmse'] for x in same)))
        if ratios:
            real[conf] = float(np.mean(ratios))
    freeze = {}
    freeze['ours_synth'] = json.loads(max(synth, key=lambda c: (synth[c], c))) if synth else None
    freeze['ours'] = json.loads(min(real, key=lambda c: (real[c], c))) if real else freeze['ours_synth']
    freeze['ours_scores'] = {'synth': dict(sorted(synth.items())), 'real': dict(sorted(real.items()))}
    for kind in ('synth', 'real'):
        lrs = {}
        for fam in ('ctrl', 'fixed', 'mlp'):
            by = {}
            for r in recs:
                if r['method'] == fam and r['kind'] == kind:
                    by.setdefault(json.loads(r['conf'])['lr'], []).append(np.log(r['val_mse']))
            if by:
                lrs[fam] = min(by, key=lambda lr: (float(np.mean(by[lr])), lr))
        freeze['lr_' + kind] = lrs
    return freeze


def dev_section(recs):
    fz = choose_freeze(recs)
    L = ['## Development stage: the configurations of ours and the frozen choice', '',
         'Ours on the synthetic development targets, mean exact-edge share (the sections above use the first configuration of ours only):', '', '| configuration | exact-edge share | |', '|---|---|---|']
    for c, sc in fz['ours_scores']['synth'].items():
        L.append('| %s | %s | %s |' % (c, fmt(sc), 'chosen (synthetic and [d, 1])' if json.loads(c) == fz['ours_synth'] else ''))
    L += ['', 'Ours (hybrid [d, %d, 1]) on the development real sets, mean log-ratio of test RMSE to the best member of the baselines (no matching of parameters):' % HID, '',
          '| configuration | log RMSE ratio | |', '|---|---|---|']
    for c, sc in fz['ours_scores']['real'].items():
        L.append('| %s | %s | %s |' % (c, fmt(sc), 'chosen' if json.loads(c) == fz['ours'] else ''))
    L += ['', 'Learning rates of the baselines (lowest mean log validation MSE): synthetic %s; real %s.' % (json.dumps(fz['lr_synth'], sort_keys=True), json.dumps(fz['lr_real'], sort_keys=True)), '']
    return L


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--stage', choices=['dev', 'test', 'smoke', 'sets', 'freeze', 'merge'])
    ap.add_argument('--out', default='.')
    ap.add_argument('--freeze')
    ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--report', action='store_true')
    ap.add_argument('--from-raw', nargs='?', const='')           # with no path: KAN_OUT/kan_raw.npz -> KAN_OUT/kan.md (default KAN_OUT: results/ml)
    ap.add_argument('--md')
    ap.add_argument('--only', choices=['synth', 'real'], help='development / test: only this kind of job (partial runs)')
    ap.add_argument('--raws', help='freeze stage: comma-separated development raws to merge')
    ap.add_argument('--archs', help='real jobs: comma-separated architectures (A, B, H)')
    ap.add_argument('--dev-ours', help='development real stage: json list overriding DEV_REAL_OURS')
    ap.add_argument('--limit', type=int, help='only the first N jobs (smoke runs)')
    ap.add_argument('--tag', help='suffix of the raw file name (partial runs)')
    ap.add_argument('--ids', help='comma-separated OpenML ids: only these real sets (partial runs)')
    ap.add_argument('--sets', help='json file with {"dev": [[did, name], ...], "test": [...]}')
    a = ap.parse_args()
    if a.report or a.from_raw is not None:
        outd = os.environ.get('KAN_OUT', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results', 'ml'))
        meta, recs = load_raw(a.from_raw or os.path.join(outd, 'kan_raw.npz'))
        txt = report(meta, recs)
        md = a.md or (os.path.join(outd, 'kan.md') if not a.from_raw else None)
        if md:
            open(md, 'w').write(txt + '\n')
        else:
            print(txt)
        return
    os.makedirs(a.out, exist_ok=True)
    if a.stage == 'sets':
        lst = real_listing()
        cache = {}

        def loaded(did):
            try:
                cache[did] = load_real(did)
                return cache[did]
            except Exception as e:
                print('skip', did, repr(e)[:80], flush=True)
                return None
        dev, test = select_real(lst, loaded)
        import hashlib
        info = {}
        for did, name in dev + test:
            X, y = cache[did]
            info[str(did)] = dict(name=name, n=int(len(y)), d=int(X.shape[1]), sha=hashlib.sha256(np.round(X, 6).tobytes() + np.round(y, 6).tobytes()).hexdigest())
        json.dump(dict(dev=dev, test=test, info=info), open(os.path.join(a.out, 'sets.json'), 'w'), indent=1)
        print(json.dumps(dict(dev=dev, test=test, info=info), indent=1))
        return
    if a.stage == 'freeze':
        recs = []
        for f in a.raws.split(','):
            recs += load_raw(f)[1]
        fz = choose_freeze(recs)
        json.dump(fz, open(os.path.join(a.out, 'kan_freeze.json'), 'w'), indent=1, sort_keys=True)
        meta = dict(stage='dev', merged=[os.path.basename(f) for f in a.raws.split(',')])
        save_raw(os.path.join(a.out, 'kan_dev_raw.npz'), meta, recs)
        open(os.path.join(a.out, 'kan_dev.md'), 'w').write(report(meta, recs) + '\n')
        print(json.dumps(fz, indent=1, sort_keys=True))
        return
    if a.stage == 'merge':                       # the partial test runs into kan_raw.npz and kan.md
        recs = []
        for f in a.raws.split(','):
            recs += load_raw(f)[1]
        meta = dict(stage='test', merged=[os.path.basename(f) for f in a.raws.split(',')])
        save_raw(os.path.join(a.out, 'kan_raw.npz'), meta, recs)
        open(os.path.join(a.out, 'kan.md'), 'w').write(report(meta, recs) + '\n')
        return
    import multiprocessing as mp
    sets = json.load(open(a.sets))
    if a.stage in ('dev', 'smoke'):
        cfg = DEV_GRID
        tids, seeds = DEV_TARGETS, DEV_SEEDS
        rsets = sets['dev']
        rseeds = (0,)
    else:
        fz = json.load(open(a.freeze))
        obya = {'A': [fz['ours_synth']], 'H': [fz['ours']]}
        cfg = dict(ours_by_arch=obya, lr_ctrl=(fz['lr_synth']['ctrl'],), lr_fixed=(fz['lr_synth']['fixed'],), lr_mlp=(fz['lr_synth']['mlp'],))
        cfg_real = dict(ours_by_arch=obya, lr_ctrl=(fz['lr_real']['ctrl'],), lr_fixed=(fz['lr_real']['fixed'],), lr_mlp=(fz['lr_real']['mlp'],), lams_real=REAL_LAMS,
                        ours_archs=('A', 'H'))
        tids, seeds = TEST_TARGETS, SEEDS
        rsets = sets['test']
        rseeds = SEEDS
    jobs = []
    if a.stage == 'smoke':
        small = dict(cfg, ours=cfg['ours'][:1], lr_free=(0.01,), lr_mlp=(0.01,), lams=(0.0, 1e-3), gs=(4, 8), hs=(4, 8))
        jobs = [('synth', DEV_TARGETS[0], 0.1, 0, 'A', small), ('real', rsets[0][0], rsets[0][1], 0, 'B', small)]
    else:
        for tid in tids:
            for nz in NOISES:
                for s in seeds:
                    jobs.append(('synth', tid, nz, s, 'A', cfg))
        rcfg = dict(cfg, ours=json.loads(a.dev_ours) if a.dev_ours else DEV_REAL_OURS, lams_real=REAL_LAMS) if a.stage == 'dev' else cfg_real
        for did, name in rsets:
            for s in rseeds:
                for arch in (a.archs.split(',') if a.archs else (('H', 'A') if a.stage == 'dev' else ('B', 'H', 'A'))):
                    jobs.append(('real', did, name, s, arch, rcfg))
    if a.only:
        jobs = [j for j in jobs if j[0] == a.only]
    if a.ids:
        keep = {int(x) for x in a.ids.split(',')}
        jobs = [j for j in jobs if j[0] == 'synth' or j[1] in keep]
    if a.limit:
        jobs = jobs[:a.limit]
    t0 = time.time()
    recs = []
    with mp.get_context('fork').Pool(a.workers) as pool:
        for rs, secs in pool.imap_unordered(job, jobs):
            recs += rs
            print('job done %d/%d %.0fs total %.0fs' % (len({run_key(r) for r in recs}), len(jobs), secs, time.time() - t0), flush=True)
    meta = dict(stage=a.stage, jobs=len(jobs), secs=time.time() - t0, cfg=json.loads(json.dumps(cfg)), sets=sets)
    save_raw(os.path.join(a.out, 'kan_%s%s%s_raw.npz' % (a.stage, '_' + a.only if a.only else '', '_' + a.tag if a.tag else '')), meta, recs)
    if a.stage == 'dev' and not a.only:
        fz = choose_freeze(recs)
        json.dump(fz, open(os.path.join(a.out, 'kan_freeze.json'), 'w'), indent=1, sort_keys=True)
        print(json.dumps(fz, indent=1, sort_keys=True))
    try:
        print(report(meta, recs)[:8000])
    except Exception as e:           # the raw is saved: the report can be regenerated
        print('report failed', repr(e))


if __name__ == '__main__':
    main()
