"""
Feasibility spike: interpretable additive models whose shape functions are small RBF networks with the factorized small-initialisation
readout (the sliding inverse scale space unit), against NAM / EBM / spline GAM and the same RBF parametrisation with a standard random
start and a penalty.

    y = b + sum_i f_i(x_i),      f_i(x) = sum_k w_ik exp(-(x - c_ik)^2 / (2 s^2)),      w_ik = nu_ik^2 a_max tanh(alpha_ik / a_max)

alpha_ik(0) = 0, nu_ik(0)^2 = eps. Every feature starts with the same K candidate bumps (K is a ceiling, not a size); plain gradient
descent switches bumps on one at a time, the others stay dormant (|w| <= THETA), and the number of active bumps per feature is whatever
the validation-based stopping time leaves. The plan, the datasets, the comparators and the GO / KILL rule are in results/ml/gam_plan.md
(committed before any test dataset was scored). This file is shipped to a private Kaggle CPU notebook (torch, interpret, xgboost, sklearn
live there; nothing but numpy / scipy is needed locally for the pure parts, which the tests cover).

    stage dev    : tune every method on the four development datasets (seeds 0, 1); writes a jsonl
    stage freeze : (local) pick each method's configuration from the dev jsonl -> results/ml/gam_frozen.json
    stage test   : run the frozen configurations once on the twelve test datasets (seeds 0..4); writes jsonl + npz
    --from-raw   : regenerate results/ml/gam.md byte-identically from results/ml/gam_raw.npz

Usage (Kaggle): python exp_gam_spike.py --stage dev|test --datasets a,b,c --out DIR [--frozen JSON] [--methods m1,m2]
"""
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.environ.get('GAM_OUT', os.path.join(HERE, '..', 'results', 'ml'))

# --------------------------------------------------------------------------------------------------------------------------------
# constants fixed by the plan
# --------------------------------------------------------------------------------------------------------------------------------
CAP = 8000                     # rows per dataset (the brief allows 20 000)
SUB_SEED = 777                 # the seeded subsample rule: rows = sorted(RandomState(SUB_SEED).permutation(n)[:CAP])
SPLIT = (0.6, 0.2, 0.2)        # train / validation / test
TEST_SEEDS = (0, 1, 2, 3, 4)
DEV_SEEDS = (0, 1)
THETA = 0.01                   # a bump is active when |w| > THETA (target in units of its training sd; logits for classification)
K_UNITS = 16                   # bump ceiling per feature (the same for every feature and dataset)
A_MAX = 3.0
GRID = 100                     # grid points on [0, 1] for stored shapes
ROLE_SEED = 20261109           # dev / test assignment
MAX_LEVELS = 10                # categorical columns with more levels are dropped; D_MAX encoded columns
D_MAX = 40
N_MIN = 1000
N_NAMES = {'reg': 8, 'clf': 8}
N_DEV = {'reg': 2, 'clf': 2}

# candidate lists in the order of the plan: (name, openml data id or None, kind). The rule `choose_datasets` takes, in this order, the
# first N_NAMES[kind] candidates that load and pass the filters (n >= N_MIN rows, 3 <= encoded columns <= D_MAX).
CANDIDATES = [
    ('california_housing', None, 'reg'), ('Bike_Sharing_Demand', 42712, 'reg'), ('concrete_compressive_strength', None, 'reg'),
    ('abalone', 183, 'reg'), ('wine_quality', 287, 'reg'), ('airfoil_self_noise', None, 'reg'), ('cpu_act', 197, 'reg'),
    ('kin8nm', 189, 'reg'), ('diamonds', 42225, 'reg'), ('house_sales', 42092, 'reg'), ('elevators', 216, 'reg'), ('pol', 201, 'reg'),
    ('adult', 1590, 'clf'), ('bank-marketing', 1461, 'clf'), ('default-of-credit-card-clients', 42477, 'clf'),
    ('MagicTelescope', 1120, 'clf'), ('higgs', 23512, 'clf'), ('spambase', 44, 'clf'), ('electricity', 151, 'clf'),
    ('credit-g', 31, 'clf'), ('phoneme', 1489, 'clf'), ('churn', 40701, 'clf'), ('eeg-eye-state', 1471, 'clf'), ('kc1', 1067, 'clf'),
]
# frozen by the rule from results/ml/gam_probe.json (see tests/test_gam_spike.py::test_registry_matches_rule); filled before any run
# (name, openml data id or None, kind, role): the rule `choose_datasets` + `assign_roles` applied to results/ml/gam_probe.json (n and the
# encoded width of every candidate, measured before any model was run)
DATASETS = [('california_housing', None, 'reg', 'dev'), ('Bike_Sharing_Demand', 42712, 'reg', 'test'),
            ('concrete_compressive_strength', None, 'reg', 'test'), ('abalone', 183, 'reg', 'test'), ('wine_quality', 287, 'reg', 'test'),
            ('airfoil_self_noise', None, 'reg', 'test'), ('cpu_act', 197, 'reg', 'dev'), ('kin8nm', 189, 'reg', 'test'),
            ('adult', 1590, 'clf', 'test'), ('bank-marketing', 1461, 'clf', 'test'), ('default-of-credit-card-clients', 42477, 'clf', 'test'),
            ('MagicTelescope', 1120, 'clf', 'dev'), ('higgs', 23512, 'clf', 'test'), ('electricity', 151, 'clf', 'test'),
            ('phoneme', 1489, 'clf', 'test'), ('churn', 40701, 'clf', 'dev')]

METHODS = ['ours_gd', 'ours_sgd', 'ours_adam', 'ctrl_l1', 'ctrl_l2', 'nam', 'nam_matched', 'ebm', 'gam_spline', 'xgb', 'mlp']
# development grids (the per-run choices inside a method, e.g. the penalty of the control, are made on that run's validation split)
GRIDS = {
    'ours_gd': [dict(eps=e, sigma=s, lr=0.5) for e in (1e-2, 1e-3, 1e-4) for s in (0.05, 0.08)]       # ids 0-5: first round
    + [dict(eps=e, sigma=s, lr=1.0) for e in (1e-2, 1e-3) for s in (0.03, 0.05, 0.08)],               # ids 6-11: second round (see the plan's amendments)
    'ours_adam': [dict(eps=1e-4, sigma=s, lr=l) for s in (0.05, 0.08) for l in (3e-3, 1e-2)],
    'ours_sgd': [dict()],                                      # no tuning of its own: the frozen ours_gd configuration, mini-batches of 512
    'ctrl_l1': [dict(sigma=s, lr=1e-2) for s in (0.03, 0.05, 0.08)],
    'ctrl_l2': [dict(sigma=s, lr=1e-2) for s in (0.03, 0.05, 0.08)],
    'nam': [dict(h=h, lr=l, wd=w) for h in (16, 48) for l in (3e-3, 1e-2) for w in (0.0, 1e-4)],
    'ebm': [dict(max_bins=b, lr=l) for b in (32, 256) for l in (0.01, 0.03)],
    'gam_spline': [dict(n_knots=k) for k in (6, 12, 24)],
    'xgb': [dict(depth=dp) for dp in (3, 6)],
    'mlp': [dict(hidden=hd, alpha=a) for hd in (64, 128) for a in (1e-4, 1e-2)],
}
LAMBDAS = {'l1': (1e-4, 3e-4, 1e-3, 3e-3, 1e-2), 'l2': (1e-5, 1e-4, 1e-3, 1e-2)}
NAM1_LADDER = (1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 24, 32, 48, 64)
T_MAX = {'gd': 6000, 'sgd': 6000, 'adam': 3000, 'ctrl': 3000}
PATIENCE = {'gd': 800, 'sgd': 800, 'adam': 300, 'ctrl': 300}
MIN_DELTA = 1e-4               # a validation improvement counts when it is at least this fraction of the best loss
NAM_EPOCHS = 300
NAM_RESTARTS = {'nam': 1, 'nam_matched': 3}      # the matched NAM, tiny and prone to dead units, gets the best of three starts on validation
EVAL_EVERY = 25
TRACE_EVERY = 250

# verdict constants (plan section 5)
MARGIN = {'reg': 0.02, 'clf': 0.005}        # non-inferiority margin: relative RMSE / absolute AUC
TIE = {'reg': 0.005, 'clf': 0.002}          # tie band of the comparison with the control
PARAM_RATIO = 0.7                           # "at least 30 % fewer active parameters"
N_NEED_A, N_NEED_B, N_NEED_C = 8, 8, 7
P_ALPHA = 0.05
EPS = 1e-9                                  # float slack on the boundaries of the rules (inclusive)


# --------------------------------------------------------------------------------------------------------------------------------
# pure parts: data rules
# --------------------------------------------------------------------------------------------------------------------------------
def choose_datasets(info, n_names=None):
    """info: {name: dict(n, d_enc, kind)} for the candidates that loaded. Takes, per kind, the first n_names[kind] candidates (in the
    CANDIDATES order) with n >= N_MIN and 3 <= d_enc <= D_MAX. Returns the names in that order."""
    n_names = n_names or N_NAMES
    got = {'reg': [], 'clf': []}
    for name, _, kind in CANDIDATES:
        r = info.get(name)
        if r is None or r['kind'] != kind or len(got[kind]) >= n_names[kind]:
            continue
        if r['n'] >= N_MIN and 3 <= r['d_enc'] <= D_MAX:
            got[kind].append(name)
    return got['reg'] + got['clf']


def assign_roles(names, kinds, seed=ROLE_SEED, n_dev=None):
    """Seeded rule: within each kind the names (sorted) are permuted by RandomState(seed); the first n_dev[kind] are development
    datasets, the rest test datasets. Returns {name: 'dev' | 'test'}."""
    n_dev = n_dev or N_DEV
    rng = np.random.RandomState(seed)
    roles = {}
    for kind in ('reg', 'clf'):
        ns = sorted(n for n in names if kinds[n] == kind)
        perm = rng.permutation(len(ns))
        for rank, j in enumerate(perm):
            roles[ns[j]] = 'dev' if rank < n_dev[kind] else 'test'
    return roles


def subsample_idx(n, cap=CAP, seed=SUB_SEED):
    if n <= cap:
        return np.arange(n)
    return np.sort(np.random.RandomState(seed).permutation(n)[:cap])


def split_idx(n, seed, frac=SPLIT):
    """Train / validation / test indices of a seeded random split."""
    p = np.random.RandomState(1000 + seed).permutation(n)
    a, b = int(round(frac[0] * n)), int(round((frac[0] + frac[1]) * n))
    return p[:a], p[a:b], p[b:]


def encode_columns(cols, max_levels=MAX_LEVELS):
    """cols: list of (name, values, is_cat) -> (X float array (n, d_enc), names). Numeric: missing -> median, constant dropped.
    Categorical: more than max_levels levels (or one level) dropped; two levels -> one 0/1 column (the later level); k >= 3 levels -> k
    indicator columns; missing -> all zeros."""
    out, names = [], []
    for name, v, is_cat in cols:
        if is_cat:
            s = np.array([None if (x is None or (isinstance(x, float) and np.isnan(x))) else str(x) for x in v], dtype=object)
            levels = sorted(set(x for x in s if x is not None))
            if len(levels) < 2 or len(levels) > max_levels:
                continue
            use = levels[1:] if len(levels) == 2 else levels
            for lv in use:
                out.append(np.array([1.0 if x == lv else 0.0 for x in s]))
                names.append(f'{name}={lv}')
        else:
            a = np.asarray(v, float)
            if np.all(np.isnan(a)):
                continue
            a = np.where(np.isnan(a), np.nanmedian(a), a)
            if a.max() == a.min():
                continue
            out.append(a)
            names.append(name)
    n = len(cols[0][1]) if cols else 0
    return (np.column_stack(out) if out else np.zeros((n, 0))), names


class Quantizer:
    """Per-column map to [0, 1]: columns with at most two distinct training values go to {0, 1}; the others by the training empirical
    CDF (mean-rank for ties, linear in between, 0 and 1 outside). The same front end for every method."""

    def fit(self, X):
        self.maps = []
        for j in range(X.shape[1]):
            u = np.unique(X[:, j])
            if len(u) <= 2:
                self.maps.append(('bin', u))
                continue
            q = np.quantile(X[:, j], np.linspace(0, 1, 513))
            uq, inv = np.unique(q, return_inverse=True)
            pr = np.zeros(len(uq))
            np.add.at(pr, inv, np.linspace(0, 1, 513))
            pr /= np.bincount(inv)
            self.maps.append(('cdf', uq, pr))
        return self

    def transform(self, X):
        Z = np.empty(X.shape, float)
        for j, m in enumerate(self.maps):
            if m[0] == 'bin':
                u = m[1]
                Z[:, j] = 0.0 if len(u) == 1 else (X[:, j] > u[0]).astype(float)
            else:
                Z[:, j] = np.interp(X[:, j], m[1], m[2])
        return Z


# --------------------------------------------------------------------------------------------------------------------------------
# pure parts: metrics and counts
# --------------------------------------------------------------------------------------------------------------------------------
def auc_score(y, s):
    from scipy.stats import rankdata
    y = np.asarray(y) > 0.5
    n1, n0 = int(y.sum()), int((~y).sum())
    if n1 == 0 or n0 == 0:
        return float('nan')
    r = rankdata(s)
    return float((r[y].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def rel_rmse(yz, pz):
    """RMSE in units of the training sd of the target (yz and pz are already standardised)."""
    return float(np.sqrt(np.mean((np.asarray(yz) - np.asarray(pz)) ** 2)))


def score_of(kind, metric):
    """Lower is better: relative RMSE for regression, 1 - AUC for classification."""
    return metric if kind == 'reg' else 1.0 - metric


def active_mask(w, theta=THETA):
    return np.abs(np.asarray(w)) > theta


def ours_params(mask):
    """Stored numbers needed to evaluate the shapes: weight and centre of every active bump (the width is one global constant)."""
    return int(2 * np.sum(mask))


def segments(values, rtol=1e-9):
    """Number of runs of (numerically) equal consecutive values: the pieces of a step function."""
    v = np.asarray(values, float)
    if len(v) == 0:
        return 0
    scale = max(np.max(np.abs(v)), 1e-300)
    return int(1 + np.sum(np.abs(np.diff(v)) > rtol * scale))


def ebm_params(term_scores):
    """A step function with s pieces is s values and s - 1 cut points."""
    return int(sum(2 * segments(np.asarray(t)[1:-1]) - 1 for t in term_scores))


def nam_params(d, h):
    return int(d * (h * h + 4 * h + 1) + 1)


def nam1_params(d, h):
    """d features, each a 1 -> h -> 1 ReLU network (3h + 1 numbers), plus the global bias."""
    return int(d * (3 * h + 1) + 1)


def nam_matched_h(d, target, ladder=NAM1_LADDER):
    """Smallest hidden width of the ladder whose d single-hidden-layer per-feature networks have at least `target` parameters in all
    (conservative for ours: the comparator never has fewer parameters)."""
    for h in ladder:
        if nam1_params(d, h) >= target:
            return h
    return ladder[-1]


def n_extrema(shape, rel_tol=0.05):
    """Turning points of a sampled curve, ignoring reversals smaller than rel_tol of its range (a hysteresis filter)."""
    s = np.asarray(shape, float)
    rng = s.max() - s.min()
    if rng <= 0:
        return 0
    tol = rel_tol * rng
    ext, direction, ref = 0, 0, s[0]
    for x in s[1:]:
        if direction == 0:
            if x > ref + tol:
                direction, ref = 1, x
            elif x < ref - tol:
                direction, ref = -1, x
        elif direction == 1:
            if x > ref:
                ref = x
            elif x < ref - tol:
                ext, direction, ref = ext + 1, -1, x
        else:
            if x < ref:
                ref = x
            elif x > ref + tol:
                ext, direction, ref = ext + 1, 1, x
    return int(ext)


def shape_stability(shapes, min_sd=1e-6):
    """shapes: (S, d, G) centred shapes of one method over S splits. For every feature the mean Pearson correlation over pairs of
    splits (a pair with one flat shape counts 0, a pair of flat shapes is skipped); features are weighted by their mean shape sd.
    Returns (weighted mean correlation, number of features that matter) or (nan, 0)."""
    S, d, _ = shapes.shape
    num = den = 0.0
    used = 0
    for j in range(d):
        sd = shapes[:, j].std(1)
        cs = []
        for a in range(S):
            for b in range(a + 1, S):
                fa, fb = sd[a] > min_sd, sd[b] > min_sd
                if fa and fb:
                    cs.append(float(np.corrcoef(shapes[a, j], shapes[b, j])[0, 1]))
                elif fa != fb:
                    cs.append(0.0)
        if cs:
            wgt = float(sd.mean())
            num += wgt * float(np.mean(cs))
            den += wgt
            used += 1
    return (num / den, used) if den > 0 else (float('nan'), 0)


# --------------------------------------------------------------------------------------------------------------------------------
# pure parts: dev -> frozen configuration, and the verdict
# --------------------------------------------------------------------------------------------------------------------------------
def pick_configs(recs):
    """recs: dev records (dict with method, cfg_id, dataset, seed, kind, score). For every method the configuration with the smallest
    mean over datasets of (mean-over-seeds score / the best configuration's mean-over-seeds score on that dataset). Ties: lowest id."""
    best = {}
    for m in sorted(set(r['method'] for r in recs)):
        rm = [r for r in recs if r['method'] == m]
        cfgs = sorted(set(r['cfg_id'] for r in rm))
        dss = sorted(set(r['dataset'] for r in rm))
        mean = {(c, d): np.mean([r['score'] for r in rm if r['cfg_id'] == c and r['dataset'] == d]) for c in cfgs for d in dss
                if any(r['cfg_id'] == c and r['dataset'] == d for r in rm)}
        ratio = {}
        for c in cfgs:
            rs = []
            for d in dss:
                if (c, d) in mean:
                    mn = min(mean[(c2, d)] for c2 in cfgs if (c2, d) in mean)
                    rs.append(mean[(c, d)] / max(mn, 1e-12))
            ratio[c] = float(np.mean(rs))
        best[m] = min(cfgs, key=lambda c: (ratio[c], c))
    return best


def paired_one_sided_p(diff):
    """p-value of H0: mean(diff) <= 0 by the paired t statistic (diff > 0 favours ours)."""
    from scipy.stats import t as tdist
    d = np.asarray(diff, float)
    n = len(d)
    if n < 2:
        return 1.0
    sd = d.std(ddof=1)
    if sd == 0:
        return 0.0 if d.mean() > 0 else 1.0
    return float(tdist.sf(d.mean() / (sd / np.sqrt(n)), n - 1))


def rel_gap(kind, s_other, s_ours):
    """How much worse `ours` is than `other`: relative for regression, absolute (1 - AUC) for classification."""
    return (s_ours - s_other) / s_other if kind == 'reg' else s_ours - s_other


def criteria_for_dataset(kind, sc, par, margin_scale=1.0, tie_scale=1.0, ctrls=('ctrl_l1', 'ctrl_l2'), p_alpha=P_ALPHA, a_params=True):
    """sc[m]: per-seed scores (lower better, aligned over seeds) for the methods ours_gd, nam, nam_matched, ebm, ctrl_l1, ctrl_l2;
    par[m]: mean active parameters. Returns the dict of the three criteria of the plan. The keyword arguments are for the POST HOC sensitivity
    section of the report only (the defaults are the pre-registered rule)."""
    ours = np.asarray(sc['ours_gd'])
    mean = {m: float(np.mean(sc[m])) for m in sc}
    best = min(('nam', 'ebm'), key=lambda m: mean[m])
    gap = rel_gap(kind, mean[best], mean['ours_gd'])
    ni = gap <= MARGIN[kind] * margin_scale + EPS
    fewer = par['ours_gd'] <= PARAM_RATIO * par[best] + EPS
    A = bool(ni and (fewer or not a_params))
    diff = np.asarray(sc['nam_matched']) - ours
    p = paired_one_sided_p(diff)
    B = bool(diff.mean() > 0 and p < p_alpha and par['nam_matched'] >= par['ours_gd'])
    cb = min(ctrls, key=lambda m: mean[m])
    g = rel_gap(kind, mean[cb], mean['ours_gd'])           # > 0: ours worse
    acc = 'win' if g < -TIE[kind] * tie_scale - EPS else ('loss' if g > TIE[kind] * tie_scale + EPS else 'tie')
    par_win = par['ours_gd'] <= PARAM_RATIO * par['ctrl_l1'] + EPS
    C = 'win' if (acc == 'win' or (acc == 'tie' and par_win)) else ('loss' if acc == 'loss' else 'tie')
    return dict(best=best, gap=float(gap), noninferior=bool(ni), fewer=bool(fewer), A=A, B=B, p_B=p, ctrl_best=cb, ctrl_gap=float(g),
                acc=acc, par_win=bool(par_win), C=C)


def verdict_from(crit):
    """crit: {dataset: criteria dict}. GO iff (A on >= 8 or B on >= 8) and C 'win' on >= 7."""
    nA = sum(c['A'] for c in crit.values())
    nB = sum(c['B'] for c in crit.values())
    nC = sum(c['C'] == 'win' for c in crit.values())
    go = (nA >= N_NEED_A or nB >= N_NEED_B) and nC >= N_NEED_C
    return dict(nA=int(nA), nB=int(nB), nC=int(nC), n=len(crit), go=bool(go))


# --------------------------------------------------------------------------------------------------------------------------------
# data loading (Kaggle) and preparation
# --------------------------------------------------------------------------------------------------------------------------------
def load_openml(name, did, kind):
    from sklearn.datasets import fetch_openml
    d = fetch_openml(data_id=did, as_frame=True, parser='auto') if did else fetch_openml(name=name, version='active', as_frame=True, parser='auto')
    X, y = d.data, d.target
    keep = ~y.isna().to_numpy()
    cols = [(c, X[c].to_numpy()[keep], str(X[c].dtype) in ('category', 'object')) for c in X.columns]
    Xe, names = encode_columns(cols)
    yv = y.to_numpy()[keep]
    if kind == 'reg':
        yv = np.asarray(yv, float)
    else:
        lv = sorted(set(map(str, yv)))
        assert len(lv) == 2, (name, lv)
        yv = (np.array(list(map(str, yv))) == lv[1]).astype(float)
    return Xe, yv, names


def probe_all():
    info = {}
    for name, did, kind in CANDIDATES:
        try:
            Xe, yv, names = load_openml(name, did, kind)
            info[name] = dict(n=int(len(yv)), d_enc=int(Xe.shape[1]), kind=kind, ok=True)
        except Exception as e:
            info[name] = dict(n=0, d_enc=0, kind=kind, ok=False, err=repr(e)[:120])
        print('PROBE', name, json.dumps(info[name]), flush=True)
    return info


def prepare(name, did, kind, seed, cache={}):
    """One split of one dataset: front end fitted on the training rows, target standardised (regression)."""
    if name not in cache:
        Xe, yv, _ = load_openml(name, did, kind)
        idx = subsample_idx(len(yv))
        cache.clear()
        cache[name] = (Xe[idx], yv[idx])
    X, y = cache[name]
    tr, va, te = split_idx(len(y), seed)
    Q = Quantizer().fit(X[tr])
    D = dict(kind=kind, d=X.shape[1], n=len(y))
    for nm, ix in (('tr', tr), ('va', va), ('te', te)):
        D['X' + nm] = Q.transform(X[ix]).astype(np.float32)
    if kind == 'reg':
        mu, sd = y[tr].mean(), y[tr].std()
        for nm, ix in (('tr', tr), ('va', va), ('te', te)):
            D['y' + nm] = ((y[ix] - mu) / sd).astype(np.float32)
        D['b0'] = 0.0
    else:
        for nm, ix in (('tr', tr), ('va', va), ('te', te)):
            D['y' + nm] = y[ix].astype(np.float32)
        p = float(np.clip(y[tr].mean(), 1e-3, 1 - 1e-3))
        D['b0'] = float(np.log(p / (1 - p)))
    return D


# --------------------------------------------------------------------------------------------------------------------------------
# models (torch parts; imported lazily)
# --------------------------------------------------------------------------------------------------------------------------------
def _metric(D, pred):
    """pred: standardised prediction (regression) or logit (classification) on the test rows."""
    if D['kind'] == 'reg':
        return rel_rmse(D['yte'], pred)
    return auc_score(D['yte'], pred)


def _grid(G=GRID):
    return np.linspace(0.0, 1.0, G)


def _init_centres(d, K, seed):
    r = np.random.RandomState(seed)
    c = (np.arange(K) + 0.5) / K
    return np.tile(c, (d, 1)) + r.uniform(-0.5, 0.5, (d, K)) / K


def _rbf_shapes(w, c, sigma):
    g = _grid()
    phi = np.exp(-0.5 * ((g[None, None, :] - c[:, :, None]) / sigma) ** 2)
    f = np.einsum('dk,dkg->dg', w, phi)
    return f - f.mean(1, keepdims=True)


def _finish_rbf(D, w, c, sigma, pred, steps, t0, extra):
    mask = active_mask(w)
    shapes = _rbf_shapes(w, c, sigma)
    return dict(pred=pred, shapes=shapes.astype(np.float32), params=ours_params(mask), units=mask.sum(1).astype(int).tolist(),
                steps=int(steps), time=time.perf_counter() - t0, extra=extra)


def _rbf_grads(X, y, z, w, sigma, b0, clf):
    """Loss and its gradients with respect to the weights w and the centres z (centres in units of the width), by hand."""
    import torch
    n = X.shape[0]
    u = X[:, :, None] / sigma - z[None]
    phi = torch.exp(-0.5 * u * u)
    out = b0 + torch.einsum('ndk,dk->n', phi, w)
    if clf:
        e = (torch.sigmoid(out) - y) / n
        loss = float(torch.nn.functional.softplus(out).sub(y * out).mean())
    else:
        e = (out - y) / n
        loss = float(0.5 * ((out - y) ** 2).mean())
    gw = torch.einsum('n,ndk->dk', e, phi)
    gz = w * torch.einsum('n,ndk->dk', e, phi * u)
    return loss, gw, gz


def _rbf_out(X, z, w, sigma, b0):
    import torch
    u = X[:, :, None] / sigma - z[None]
    return b0 + torch.einsum('ndk,dk->n', torch.exp(-0.5 * u * u), w)


def _val_loss(out, y, clf):
    import torch
    return float(torch.nn.functional.softplus(out).sub(y * out).mean()) if clf else float(0.5 * ((out - y) ** 2).mean())


def train_rbf(D, cfg, seed, method, lam=None):
    """The shared training loop of the RBF additive model (gradients by hand; checked against autograd by `check_grads`).
    method 'gd' | 'sgd' | 'adam': factorised small-start readout (ours); 'l1' | 'l2': free weights from a random start, trained by Adam with an
    L1 penalty (proximal step in Adam's diagonal metric) or an L2 penalty (the control).
    Stopping (all methods): the parameters at the last validation improvement of at least MIN_DELTA (relative), after PATIENCE steps
    without one (for the factorised readout only once a bump is active; for gd / sgd the patience is max(PATIENCE, t_best): the path is
    in logarithmic time, an entry can take as long again as everything so far), or at T_MAX."""
    import torch
    torch.set_num_threads(int(os.environ.get('GAM_THREADS', '4')))
    t0 = time.perf_counter()
    d, K, sigma = D['d'], K_UNITS, cfg['sigma']
    clf = D['kind'] == 'clf'
    Xtr, ytr, Xva, yva, Xte = (torch.as_tensor(D[k]) for k in ('Xtr', 'ytr', 'Xva', 'yva', 'Xte'))
    b0 = D['b0']
    fact = method in ('gd', 'sgd', 'adam')
    lr0 = cfg['lr'] * (4.0 if (clf and method in ('gd', 'sgd')) else 1.0)
    if method in ('gd', 'sgd'):
        lr0 *= min(1.0, 12.0 / max(d, 12))                   # curvature grows with the number of overlapping shapes
    key = 'ctrl' if not fact else method
    tmax, pat = T_MAX[key], PATIENCE[key]
    ntr = Xtr.shape[0]
    for attempt in range(4):                                   # divergence guard: halve the step and restart
        lr = lr0 / 2 ** attempt
        gen = torch.Generator().manual_seed(seed)
        z = torch.tensor(_init_centres(d, K, seed) / sigma, dtype=torch.float32)
        if fact:
            al = torch.zeros(d, K)
            nu = torch.full((d, K), float(np.sqrt(cfg['eps'])))
            params = [al, nu, z]
        else:
            wf = 0.1 * torch.randn(d, K, generator=gen)
            params = [wf, z]
            mom, var = [torch.zeros_like(p) for p in params], [torch.zeros_like(p) for p in params]
        if method == 'adam':
            for p in params:
                p.requires_grad_(True)
            opt = torch.optim.Adam(params, lr=lr)

        def weights():
            return nu * nu * A_MAX * torch.tanh(al / A_MAX) if fact else wf

        best, best_t, best_snap, ok, trace = float('inf'), 0, None, True, []
        perm, pos, bs = None, 0, 512
        for t in range(tmax + 1):
            if t % EVAL_EVERY == 0:
                with torch.no_grad():
                    w = weights()
                    lv = _val_loss(_rbf_out(Xva, z, w, sigma, b0), yva, clf)
                    if not np.isfinite(lv):
                        ok = False
                        break
                    if lv < best * (1 - MIN_DELTA) or best_snap is None:
                        best, best_t = lv, t
                        best_snap = (w.detach().clone(), z.detach().clone())
                    if t % TRACE_EVERY == 0:
                        trace.append((t, int((w.abs() > THETA).sum())))
                    any_active = bool((w.abs() > THETA).any()) if fact else True
                    if any_active and t - best_t > (max(pat, best_t) if method in ('gd', 'sgd') else pat):
                        break
            if t == tmax:
                break
            if method == 'sgd':
                if perm is None or pos + bs > ntr:
                    perm, pos = torch.randperm(ntr, generator=gen), 0
                ix = perm[pos:pos + bs]
                pos += bs
                Xb, yb = Xtr[ix], ytr[ix]
            else:
                Xb, yb = Xtr, ytr
            with torch.no_grad():
                w = weights()
                _, gw, gz = _rbf_grads(Xb, yb, z, w, sigma, b0, clf)
                if method == 'l2':
                    gw = gw + lam * w
                if fact:
                    th = torch.tanh(al / A_MAX)
                    g_al = gw * nu * nu * (1 - th * th)
                    g_nu = gw * 2 * nu * A_MAX * th
                    grads = [g_al, g_nu, gz]
                else:
                    grads = [gw, gz]
                if method == 'adam':
                    for p, g in zip(params, grads):
                        p.grad = g
                    opt.step()
                elif not fact:                                  # the control: Adam, the L1 penalty as a proximal step in Adam's metric
                    for i, (p, g) in enumerate(zip(params, grads)):
                        mom[i] = 0.9 * mom[i] + 0.1 * g
                        var[i] = 0.999 * var[i] + 0.001 * g * g
                        step = lr * (mom[i] / (1 - 0.9 ** (t + 1))) / (torch.sqrt(var[i] / (1 - 0.999 ** (t + 1))) + 1e-8)
                        if method == 'l1' and i == 0:
                            u_ = p - step
                            thr = lr * lam / (torch.sqrt(var[i] / (1 - 0.999 ** (t + 1))) + 1e-8)
                            p.copy_(torch.sign(u_) * torch.clamp(u_.abs() - thr, min=0.0))
                        else:
                            p -= step
                else:
                    for p, g in zip(params, grads):
                        p -= lr * g
        if ok and best_snap is not None:
            break
    w, zz = best_snap
    with torch.no_grad():
        pred = _rbf_out(Xte, zz, w, sigma, b0).numpy()
    return _finish_rbf(D, w.numpy(), (zz * sigma).numpy(), sigma, pred, best_t, t0,
                       dict(val_loss=float(best), lam=lam, attempts=attempt, trace=trace, cfg=cfg))


def check_grads(seed=0):
    """Hand-written gradients of `_rbf_grads` against autograd (regression and classification)."""
    import torch
    rng = np.random.RandomState(seed)
    X = torch.tensor(rng.rand(50, 3), dtype=torch.float64)
    z = torch.tensor(rng.rand(3, 4) * 10, dtype=torch.float64, requires_grad=True)
    w = torch.tensor(rng.randn(3, 4), dtype=torch.float64, requires_grad=True)
    worst = 0.0
    for clf in (False, True):
        y = torch.tensor((rng.rand(50) > 0.5) if clf else rng.randn(50), dtype=torch.float64)
        sigma, b0 = 0.1, 0.3
        u = X[:, :, None] / sigma - z[None]
        out = b0 + torch.einsum('ndk,dk->n', torch.exp(-0.5 * u * u), w)
        loss = torch.nn.functional.binary_cross_entropy_with_logits(out, y) if clf else 0.5 * ((out - y) ** 2).mean()
        gw, gz = torch.autograd.grad(loss, [w, z])
        _, hw, hz = _rbf_grads(X, y, z.detach(), w.detach(), sigma, b0, clf)
        worst = max(worst, float((gw - hw).abs().max()), float((gz - hz).abs().max()))
    return worst


def run_ctrl(D, cfg, seed, pen):
    """The control: per run, the penalty is chosen on the validation split from LAMBDAS[pen] (more tuning than ours gets)."""
    import torch
    best = None
    for lam in LAMBDAS[pen]:
        r = train_rbf(D, cfg, seed, pen, lam)
        if best is None or r['extra']['val_loss'] < best['extra']['val_loss']:
            best = r
    best['extra']['lam_grid'] = list(LAMBDAS[pen])
    return best


def train_nam(D, cfg, seed, h=None, layers=2, sel=None):
    import torch
    import torch.nn.functional as F
    torch.set_num_threads(int(os.environ.get('GAM_THREADS', '4')))
    t0 = time.perf_counter()
    torch.manual_seed(seed)
    h = h or cfg['h']
    d_all = D['d']
    sel = list(range(d_all)) if sel is None else list(sel)
    d = len(sel)
    clf = D['kind'] == 'clf'
    Xtr, Xva, Xte = (torch.as_tensor(D[k][:, sel]) for k in ('Xtr', 'Xva', 'Xte'))
    ytr, yva = torch.as_tensor(D['ytr']), torch.as_tensor(D['yva'])
    W1 = torch.where(torch.rand(d, h) < 0.5, -1.0, 1.0) * (0.5 + torch.rand(d, h))       # ReLU kinks spread over [0, 1]: b1 = -W1 x0
    b1 = (-W1 * torch.rand(d, h)).requires_grad_(True)
    W1.requires_grad_(True)
    if layers == 2:
        W2 = (torch.randn(d, h, h) / np.sqrt(h)).requires_grad_(True)
        b2 = torch.zeros(d, h, requires_grad=True)
    W3 = (torch.randn(d, h) / np.sqrt(h)).requires_grad_(True)
    B0 = torch.tensor(float(D['b0']), requires_grad=True)
    params = [W1, b1, W2, b2, W3, B0] if layers == 2 else [W1, b1, W3, B0]

    def feats(X):
        h1 = torch.relu(X[:, :, None] * W1[None] + b1[None])
        if layers == 2:
            h1 = torch.relu(torch.einsum('bdh,dhk->bdk', h1, W2) + b2[None])
        return (h1 * W3[None]).sum(-1)

    def out_of(X):
        return B0 + feats(X).sum(1)

    def loss_of(o, y):
        return F.binary_cross_entropy_with_logits(o, y) if clf else 0.5 * ((o - y) ** 2).mean()

    opt = torch.optim.Adam(params, lr=cfg['lr'], weight_decay=cfg.get('wd', 0.0))
    best, best_ep, snap, ntr = float('inf'), 0, None, Xtr.shape[0]
    for ep in range(NAM_EPOCHS):
        perm = torch.randperm(ntr)
        for i in range(0, ntr, 512):
            ix = perm[i:i + 512]
            opt.zero_grad()
            loss_of(out_of(Xtr[ix]), ytr[ix]).backward()
            opt.step()
        with torch.no_grad():
            lv = float(loss_of(out_of(Xva), yva))
        if lv < best:
            best, best_ep = lv, ep
            snap = [p.detach().clone() for p in params]
        if ep - best_ep > 20:
            break
    with torch.no_grad():
        for p, s in zip(params, snap):
            p.copy_(s)
        pred = out_of(Xte).numpy()
        g = torch.as_tensor(np.tile(_grid()[:, None], (1, d)).astype(np.float32))
        fs = feats(g).numpy().T
    f = np.zeros((d_all, len(_grid())), np.float32)
    f[sel] = fs
    f = f - f.mean(1, keepdims=True)
    npar = nam_params(d, h) if layers == 2 else nam1_params(d, h)
    return dict(pred=pred, shapes=f.astype(np.float32), params=npar, units=None, steps=best_ep + 1,
                time=time.perf_counter() - t0, extra=dict(val_loss=best, h=h, layers=layers, n_sel=d, cfg=cfg))


def _tab_shapes_from_blocks(blocks):
    f = np.stack(blocks)
    return (f - f.mean(1, keepdims=True)).astype(np.float32)


def run_ebm(D, cfg, seed):
    from interpret.glassbox import ExplainableBoostingClassifier, ExplainableBoostingRegressor
    t0 = time.perf_counter()
    Xf = np.vstack([D['Xtr'], D['Xva']])
    yf = np.concatenate([D['ytr'], D['yva']])
    kw = dict(interactions=0, max_bins=cfg['max_bins'], learning_rate=cfg['lr'], random_state=seed, n_jobs=4)
    if D['kind'] == 'reg':
        m = ExplainableBoostingRegressor(**kw).fit(Xf, yf)
        pred = m.predict(D['Xte'])
    else:
        m = ExplainableBoostingClassifier(**kw).fit(Xf, yf.astype(int))
        pred = m.predict_proba(D['Xte'])[:, 1]
    d = D['d']
    g = np.tile(_grid()[:, None], (1, d))
    contrib = np.asarray(m.eval_terms(g))
    assert contrib.shape == (len(g), d)
    f = contrib.T
    return dict(pred=pred, shapes=(f - f.mean(1, keepdims=True)).astype(np.float32), params=ebm_params(m.term_scores_), units=None,
                steps=0, time=time.perf_counter() - t0, extra=dict(cfg=cfg))


def run_spline(D, cfg, seed):
    from sklearn.linear_model import LogisticRegression, Ridge
    from sklearn.preprocessing import SplineTransformer
    t0 = time.perf_counter()
    d, k = D['d'], cfg['n_knots']
    st = SplineTransformer(n_knots=k, degree=3, knots='uniform', include_bias=False, extrapolation='constant')
    st.fit(np.vstack([[np.zeros(d)], [np.ones(d)]]))
    B = {nm: st.transform(D['X' + nm]) for nm in ('tr', 'va', 'te')}
    nb = B['tr'].shape[1] // d
    best = None
    for a in (1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0):
        if D['kind'] == 'reg':
            m = Ridge(alpha=a).fit(B['tr'], D['ytr'])
            pv = m.predict(B['va'])
            v = float(np.mean((pv - D['yva']) ** 2))
        else:
            m = LogisticRegression(C=1.0 / (a * len(D['ytr']) * 10.0 + 1e-9), max_iter=300).fit(B['tr'], D['ytr'].astype(int))
            p = np.clip(m.predict_proba(B['va'])[:, 1], 1e-6, 1 - 1e-6)
            v = float(-np.mean(D['yva'] * np.log(p) + (1 - D['yva']) * np.log(1 - p)))
        if best is None or v < best[0]:
            best = (v, a, m)
    v, a, m = best
    pred = m.predict(B['te']) if D['kind'] == 'reg' else m.decision_function(B['te'])
    coef = np.ravel(m.coef_)
    G = st.transform(np.tile(_grid()[:, None], (1, d)))
    f = _tab_shapes_from_blocks([G[:, j * nb:(j + 1) * nb] @ coef[j * nb:(j + 1) * nb] for j in range(d)])
    return dict(pred=pred, shapes=f, params=int(d * nb), units=None, steps=0, time=time.perf_counter() - t0, extra=dict(alpha=a, cfg=cfg))


def run_xgb(D, cfg, seed):
    import xgboost as xgb
    t0 = time.perf_counter()
    kw = dict(n_estimators=2000, learning_rate=0.05, max_depth=cfg['depth'], subsample=0.8, colsample_bytree=0.8, early_stopping_rounds=50,
              n_jobs=4, random_state=seed)
    if D['kind'] == 'reg':
        m = xgb.XGBRegressor(**kw).fit(D['Xtr'], D['ytr'], eval_set=[(D['Xva'], D['yva'])], verbose=False)
        pred = m.predict(D['Xte'])
    else:
        m = xgb.XGBClassifier(**kw).fit(D['Xtr'], D['ytr'], eval_set=[(D['Xva'], D['yva'])], verbose=False)
        pred = m.predict_proba(D['Xte'])[:, 1]
    return dict(pred=pred, shapes=None, params=int(m.best_iteration + 1) * (2 ** cfg['depth']), units=None, steps=int(m.best_iteration + 1),
                time=time.perf_counter() - t0, extra=dict(cfg=cfg))


def run_mlp(D, cfg, seed):
    from sklearn.neural_network import MLPClassifier, MLPRegressor
    t0 = time.perf_counter()
    Xf = np.vstack([D['Xtr'], D['Xva']])
    yf = np.concatenate([D['ytr'], D['yva']])
    kw = dict(hidden_layer_sizes=(cfg['hidden'], cfg['hidden']), alpha=cfg['alpha'], early_stopping=True, n_iter_no_change=15,
              max_iter=300, random_state=seed)
    if D['kind'] == 'reg':
        m = MLPRegressor(**kw).fit(Xf, yf)
        pred = m.predict(D['Xte'])
    else:
        m = MLPClassifier(**kw).fit(Xf, yf.astype(int))
        pred = m.predict_proba(D['Xte'])[:, 1]
    npar = sum(c.size for c in m.coefs_) + sum(b.size for b in m.intercepts_)
    return dict(pred=pred, shapes=None, params=int(npar), units=None, steps=int(m.n_iter_), time=time.perf_counter() - t0, extra=dict(cfg=cfg))


def run_method(method, D, cfg, seed, ours_params_=None, ours_units_=None):
    if method == 'ours_gd':
        return train_rbf(D, cfg, seed, 'gd')
    if method == 'ours_sgd':
        return train_rbf(D, cfg, seed, 'sgd')
    if method == 'ours_adam':
        return train_rbf(D, cfg, seed, 'adam')
    if method == 'ctrl_l1':
        return run_ctrl(D, cfg, seed, 'l1')
    if method == 'ctrl_l2':
        return run_ctrl(D, cfg, seed, 'l2')
    if method == 'nam':
        return train_nam(D, cfg, seed)
    if method == 'nam_matched':
        sel = [j for j in range(D['d']) if ours_units_[j] > 0] or [0]          # the features ours uses (a free feature selection for the comparator)
        h = nam_matched_h(len(sel), ours_params_)
        rs = [train_nam(D, cfg, seed + 100 * r, h=h, layers=1, sel=sel) for r in range(NAM_RESTARTS['nam_matched'])]
        best = min(rs, key=lambda r: r['extra']['val_loss'])
        best['time'] = sum(r['time'] for r in rs)
        return best
    if method == 'ebm':
        return run_ebm(D, cfg, seed)
    if method == 'gam_spline':
        return run_spline(D, cfg, seed)
    if method == 'xgb':
        return run_xgb(D, cfg, seed)
    if method == 'mlp':
        return run_mlp(D, cfg, seed)
    raise KeyError(method)


def shape_summary(shapes, sd_floor=0.05):
    """Used features (shape sd above sd_floor in target-sd units) and mean number of turning points over them."""
    if shapes is None:
        return None, None
    sd = shapes.std(1)
    used = np.where(sd > sd_floor)[0]
    if len(used) == 0:
        return 0, 0.0
    return int(len(used)), float(np.mean([n_extrema(shapes[j]) for j in used]))


def merged_cfg(method, cfg, frozen):
    """ours_sgd inherits (eps, sigma) from the frozen ours_gd; nam_matched uses the frozen nam's lr and weight decay."""
    if method == 'ours_sgd':
        return dict(frozen['ours_gd'])
    if method == 'nam_matched':
        return dict(frozen['nam'])
    return cfg


def run_stage(stage, ds_list, out_dir, frozen=None, methods=None, min_id=0):
    os.makedirs(out_dir, exist_ok=True)
    seeds = DEV_SEEDS if stage == 'dev' else TEST_SEEDS
    methods = methods or METHODS
    jl = os.path.join(out_dir, f'gam_{stage}_{"_".join(n[:6] for n, _, _, _ in ds_list)}.jsonl')
    shapes_store = {}
    for name, did, kind, role in ds_list:
        for seed in seeds:
            D = prepare(name, did, kind, seed)
            ours_p = ours_u = None
            order = [m for m in METHODS if m in methods]
            for m in order:
                if stage == 'dev':
                    cfgs = [] if m in ('nam_matched', 'ours_sgd') else [(i, c) for i, c in enumerate(GRIDS[m]) if i >= (min_id if m == 'ours_gd' else 0)]
                else:
                    cfgs = [(frozen.get(m + '_id', 0), merged_cfg(m, frozen.get(m), frozen))]
                for cid, cfg in cfgs:
                    if m == 'nam_matched' and ours_p is None:
                        continue
                    try:
                        r = run_method(m, D, cfg, seed, ours_p, ours_u)
                    except Exception as e:
                        rec = dict(dataset=name, kind=kind, seed=seed, method=m, cfg_id=cid, error=repr(e)[:300])
                        open(jl, 'a').write(json.dumps(rec) + '\n')
                        print('ERROR', name, seed, m, repr(e)[:200], flush=True)
                        continue
                    met = _metric(D, r['pred'])
                    nused, ext = shape_summary(r['shapes'])
                    rec = dict(dataset=name, kind=kind, role=role, seed=seed, method=m, cfg_id=cid, metric=met, score=score_of(kind, met),
                               params=r['params'], units=r['units'], steps=r['steps'], time=round(r['time'], 2), d=D['d'], n_used_feats=nused,
                               extrema=ext, extra=json.loads(json.dumps(r['extra'], default=float)))
                    open(jl, 'a').write(json.dumps(rec) + '\n')
                    print('REC', name, seed, m, cid, round(met, 4), r['params'], round(r['time'], 1), flush=True)
                    if m == 'ours_gd':
                        ours_p, ours_u = r['params'], r['units']
                    if stage == 'test' and r['shapes'] is not None:
                        shapes_store.setdefault(f'shp__{name}__{m}', {})[seed] = r['shapes']
    if shapes_store:
        arrs = {k: np.stack([v[s] for s in sorted(v)]) for k, v in shapes_store.items()}
        np.savez_compressed(os.path.join(out_dir, f'gam_shapes_{"_".join(n[:6] for n, _, _, _ in ds_list)}.npz'), **arrs)
    return jl


# --------------------------------------------------------------------------------------------------------------------------------
# report
# --------------------------------------------------------------------------------------------------------------------------------
def load_records(path):
    z = np.load(path, allow_pickle=False)
    recs = [json.loads(s) for s in z['records_json']]
    shapes = {k: z[k] for k in z.files if k.startswith('shp__')}
    meta = json.loads(str(z['meta_json']))
    return recs, shapes, meta


def save_raw(path, recs, shapes, meta):
    np.savez_compressed(path, records_json=np.array([json.dumps(r, sort_keys=True) for r in recs]), meta_json=np.array(json.dumps(meta, sort_keys=True)),
                        **shapes)


def _by(recs):
    d = {}
    for r in recs:
        if 'error' in r:
            continue
        d.setdefault((r['dataset'], r['method']), {})[r['seed']] = r
    return d


def _fmt(x, nd=4):
    return 'nan' if x is None or (isinstance(x, float) and np.isnan(x)) else f'{x:.{nd}f}'


def dev_summary_lines(dev):
    """Markdown lines summarising the development (tuning) records: per method and configuration the mean over the development datasets of
    (mean-over-seeds score / the best configuration's score on that dataset) and the mean stored parameters; the chosen configuration is marked."""
    L = ['## Development (tuning) summary', '',
         'Development datasets and seeds only. "ratio": mean over datasets of (score / the best configuration\'s score on that dataset), the quantity `pick_configs` minimises; '
         '"*" marks the configuration frozen for the test run. `ctrlgd_*` is the first-round gradient-descent control (replaced by the Adam control; never a candidate).', '',
         '| method | id | configuration | ratio | mean stored params |', '|---|---|---|---|---|']
    cand = [r for r in dev if r['method'] in GRIDS]
    pick = pick_configs(cand)
    for m in METHODS:
        rm = [r for r in cand if r['method'] == m]
        if not rm:
            continue
        cfgs = sorted(set(r['cfg_id'] for r in rm))
        dss = sorted(set(r['dataset'] for r in rm))
        mean = {(c, d): float(np.mean([r['score'] for r in rm if r['cfg_id'] == c and r['dataset'] == d])) for c in cfgs for d in dss
                if any(r['cfg_id'] == c and r['dataset'] == d for r in rm)}
        for c in cfgs:
            ratios = [mean[(c, d)] / max(min(mean[(c2, d)] for c2 in cfgs if (c2, d) in mean), 1e-12) for d in dss if (c, d) in mean]
            par = float(np.mean([r['params'] for r in rm if r['cfg_id'] == c]))
            L.append(f'| {m} | {c}{"*" if pick[m] == c else ""} | {json.dumps(GRIDS[m][c], sort_keys=True)} | {np.mean(ratios):.4f} | {par:.0f} |')
    for m in sorted(set(r['method'] for r in dev if r['method'] not in GRIDS)):
        rm = [r for r in dev if r['method'] == m and 'error' not in r]
        for c in sorted(set(r['cfg_id'] for r in rm)):
            rc = [r for r in rm if r['cfg_id'] == c]
            L.append(f'| {m} | {c} | first round, plain gradient descent | - | {np.mean([r["params"] for r in rc]):.0f} |')
    return L + ['']


def render_report(recs, shapes, meta):
    by = _by(recs)
    dss = sorted(set(r['dataset'] for r in recs if 'error' not in r))
    kinds = {r['dataset']: r['kind'] for r in recs if 'error' not in r}
    seeds = sorted(set(r['seed'] for r in recs if 'error' not in r))
    L = []
    P = L.append
    P('# Additive RBF shape functions with the factorized small start: feasibility spike (test datasets)')
    P('')
    P(f'Generated by `experiments/exp_gam_spike.py --from-raw` from `gam_raw.npz` ({len(recs)} records, {len(dss)} test datasets, seeds {seeds[0]}..{seeds[-1]}). '
      'Plan and rule: `gam_plan.md`; frozen configurations: `gam_frozen.json`. Metric: relative RMSE (test RMSE in training-sd units, lower is better) for '
      'regression; AUC for classification. "params": numbers stored to evaluate the shapes (ours and the control: weight and centre of every bump with '
      f'|w| > {THETA}; NAM: all weights; EBM: values and cut points of the step functions; spline: coefficients).')
    P('')

    def val(ds, m):
        rs = by.get((ds, m), {})
        return np.array([rs[s]['metric'] for s in sorted(rs)]) if rs else None

    def sco(ds, m):
        rs = by.get((ds, m), {})
        return np.array([rs[s]['score'] for s in sorted(rs)]) if rs else None

    def pars(ds, m):
        rs = by.get((ds, m), {})
        return float(np.mean([rs[s]['params'] for s in rs])) if rs else float('nan')

    for kind, title in (('reg', 'Regression: relative RMSE (mean +- sd over seeds)'), ('clf', 'Classification: AUC (mean +- sd over seeds)')):
        ds_k = [d for d in dss if kinds[d] == kind]
        if not ds_k:
            continue
        P(f'## {title}')
        P('')
        P('| dataset | ' + ' | '.join(METHODS) + ' |')
        P('|---|' + '---|' * len(METHODS))
        for ds in ds_k:
            row = []
            for m in METHODS:
                v = val(ds, m)
                row.append('-' if v is None else f'{v.mean():.4f}+-{v.std(ddof=1) if len(v) > 1 else 0:.4f}')
            P(f'| {ds} | ' + ' | '.join(row) + ' |')
        P('')
    P('## Stored parameters (mean over seeds)')
    P('')
    P('| dataset | ' + ' | '.join(METHODS) + ' |')
    P('|---|' + '---|' * len(METHODS))
    for ds in dss:
        P(f'| {ds} | ' + ' | '.join(f'{pars(ds, m):.0f}' if not np.isnan(pars(ds, m)) else '-' for m in METHODS) + ' |')
    P('')
    P('## Bumps per feature chosen by training (ours, full-batch GD)')
    P('')
    P(f'K = {K_UNITS} candidate bumps per feature at the start, for every feature. active = |w| > {THETA}.')
    P('')
    P('| dataset | d | features with >=1 active bump | active bumps (total) | per used feature: min / median / max | Adam: active bumps | SGD: active bumps | control L1: active bumps |')
    P('|---|---|---|---|---|---|---|---|')

    def tot(ds, m):
        rs = by.get((ds, m), {})
        return float(np.mean([sum(rs[s]['units']) for s in rs])) if rs else float('nan')

    for ds in dss:
        rs = by.get((ds, 'ours_gd'), {})
        if not rs:
            continue
        allu = np.array([rs[s]['units'] for s in sorted(rs)])
        used = allu[allu > 0]
        nf = float(np.mean((allu > 0).sum(1)))
        mm = f'{used.min()} / {np.median(used):.0f} / {used.max()}' if len(used) else '-'
        P(f'| {ds} | {rs[min(rs)]["d"]} | {nf:.1f} | {tot(ds, "ours_gd"):.1f} | {mm} | {tot(ds, "ours_adam"):.1f} | {tot(ds, "ours_sgd"):.1f} | {tot(ds, "ctrl_l1"):.1f} |')
    P('')

    # criteria
    P('## Pre-registered criteria (ours = full-batch GD row)')
    P('')
    crit = {}
    need = ['ours_gd', 'nam', 'nam_matched', 'ebm', 'ctrl_l1', 'ctrl_l2']
    for ds in dss:
        if all(sco(ds, m) is not None and len(sco(ds, m)) == len(seeds) for m in need):
            crit[ds] = criteria_for_dataset(kinds[ds], {m: sco(ds, m) for m in need}, {m: pars(ds, m) for m in need})
    P('A = within the margin of the better of NAM / EBM (relative RMSE +2 %, AUC -0.005) AND stored parameters <= 0.7 x those of that comparator. '
      'B = beats the NAM with at least as many parameters (smallest ladder width), paired one-sided t-test p < 0.05. '
      'C = against the better of the two controls: win if more accurate beyond the tie band (0.5 % relative RMSE, 0.002 AUC), or tied in accuracy with '
      '<= 0.7 x the parameters of the L1 control; loss if less accurate beyond the band.')
    P('')
    P('| dataset | best of NAM/EBM | gap to it | A | params ours / comparator | B (p) | NAM-matched params | C vs control | control gap |')
    P('|---|---|---|---|---|---|---|---|---|')
    for ds, c in crit.items():
        P(f'| {ds} | {c["best"]} | {c["gap"]:+.4f} | {"yes" if c["A"] else "no"} | {pars(ds, "ours_gd"):.0f} / {pars(ds, c["best"]):.0f} | '
          f'{"yes" if c["B"] else "no"} ({c["p_B"]:.3f}) | {pars(ds, "nam_matched"):.0f} | {c["C"]} ({c["ctrl_best"]}, par {"yes" if c["par_win"] else "no"}) | {c["ctrl_gap"]:+.4f} |')
    P('')
    v = verdict_from(crit)
    P('## Verdict')
    P('')
    P(f'A on {v["nA"]} of {v["n"]} (needs {N_NEED_A}); B on {v["nB"]} of {v["n"]} (needs {N_NEED_B}); C (win vs control) on {v["nC"]} of {v["n"]} (needs {N_NEED_C}).')
    P('')
    P(f'**{"GO" if v["go"] else "KILL"}** = (A >= {N_NEED_A} or B >= {N_NEED_B}) and C >= {N_NEED_C}.')
    P('')
    if meta.get('dev'):
        L.extend(dev_summary_lines(meta['dev']))
    # POST HOC sensitivity
    P('## POST HOC sensitivity (written after the test results were read; not the pre-registered rule)')
    P('')
    P('The same per-dataset criteria with one change at a time. The A variant drops the parameter condition and widens the margin to 1.5x; the control variants change the tie band or the control set. '
      'GO would need (A >= 8 or B >= 8) and C >= 7 in every row.')
    P('')
    P('| variant | A | B | C wins | C losses | would be GO |')
    P('|---|---|---|---|---|---|')
    variants = [('pre-registered rule', {}), ('C: tie band x2 (1 % / 0.004)', dict(tie_scale=2.0)), ('C: tie band x4 (2 % / 0.008)', dict(tie_scale=4.0)),
                ('C: against ctrl_l1 only', dict(ctrls=('ctrl_l1',))), ('C: against ctrl_l2 only', dict(ctrls=('ctrl_l2',))),
                ('A: margin x1.5, no parameter condition', dict(margin_scale=1.5, a_params=False)), ('B: p < 0.10', dict(p_alpha=0.10))]
    for name, kw in variants:
        cr = {ds: criteria_for_dataset(kinds[ds], {m: sco(ds, m) for m in need}, {m: pars(ds, m) for m in need}, **kw) for ds in crit}
        vv = verdict_from(cr)
        nloss = sum(c['C'] == 'loss' for c in cr.values())
        P(f'| {name} | {vv["nA"]} | {vv["nB"]} | {vv["nC"]} | {nloss} | {"yes" if vv["go"] else "no"} |')
    P('')
    wins = {m: sum(1 for ds in dss if sco(ds, m) is not None and sco(ds, 'ours_gd') is not None
                   and rel_gap(kinds[ds], float(np.mean(sco(ds, m))), float(np.mean(sco(ds, 'ours_gd')))) <= TIE[kinds[ds]] + EPS) for m in ('ctrl_l1', 'ctrl_l2')}
    P(f'Datasets on which ours_gd is not worse than the control beyond the tie band: ctrl_l1 {wins["ctrl_l1"]} of {len(dss)}, ctrl_l2 {wins["ctrl_l2"]} of {len(dss)}. '
      f'Ratio of mean stored parameters ours_gd / ctrl_l1: {np.mean([pars(ds, "ours_gd") for ds in dss]) / np.mean([pars(ds, "ctrl_l1") for ds in dss]):.2f}.')
    P('')
    # supplementary
    P('## Supplementary (not part of the rule)')
    P('')
    P('Mean over test datasets of the score gap to the full-batch GD row (relative for regression, absolute 1 - AUC for classification; positive = worse than ours):')
    P('')
    P('| method | mean gap | datasets better than ours (beyond the tie band) | mean stored params | mean time (s) |')
    P('|---|---|---|---|---|')
    for m in METHODS:
        gaps, better, ps, ts = [], 0, [], []
        for ds in dss:
            a, b = sco(ds, m), sco(ds, 'ours_gd')
            if a is None or b is None:
                continue
            g = rel_gap(kinds[ds], float(np.mean(b)), float(np.mean(a)))
            gaps.append(g)
            better += g < -TIE[kinds[ds]]
            ps.append(pars(ds, m))
            ts.append(np.mean([by[(ds, m)][s]['time'] for s in by[(ds, m)]]))
        if gaps:
            P(f'| {m} | {np.mean(gaps):+.4f} | {better} of {len(gaps)} | {np.mean(ps):.0f} | {np.mean(ts):.1f} |')
    P('')
    P('Shape stability across the five splits (importance-weighted mean Pearson correlation of the shape of a feature between two splits, mean over datasets) '
      'and complexity (turning points per used feature, mean over datasets and seeds):')
    P('')
    P('| method | stability | turning points per used feature |')
    P('|---|---|---|')
    for m in ('ours_gd', 'ours_sgd', 'ours_adam', 'ctrl_l1', 'ctrl_l2', 'nam', 'nam_matched', 'ebm', 'gam_spline'):
        st, ex = [], []
        for ds in dss:
            arr = shapes.get(f'shp__{ds}__{m}')
            if arr is not None and len(arr) > 1:
                s_, _ = shape_stability(arr.astype(float))
                if not np.isnan(s_):
                    st.append(s_)
            for s in by.get((ds, m), {}).values():
                if s.get('extrema') is not None and s.get('n_used_feats'):
                    ex.append(s['extrema'])
        P(f'| {m} | {_fmt(float(np.mean(st)), 3) if st else "-"} | {_fmt(float(np.mean(ex)), 2) if ex else "-"} |')
    P('')
    P(f'Divergence restarts (step halved): {sum(r.get("extra", {}).get("attempts", 0) > 0 for r in recs if "error" not in r)} runs. '
      f'Errors: {sum("error" in r for r in recs)}.')
    return '\n'.join(L) + '\n'


# --------------------------------------------------------------------------------------------------------------------------------
def main(argv):
    a = {k: v for k, v in zip(argv[::2], argv[1::2])} if '--from-raw' not in argv else {}
    if '--from-raw' in argv:
        recs, shapes, meta = load_records(os.path.join(OUT, 'gam_raw.npz'))
        open(os.path.join(OUT, 'gam.md'), 'w').write(render_report(recs, shapes, meta))
        print('wrote', os.path.join(OUT, 'gam.md'))
        return
    stage = a['--stage']
    out = a.get('--out', OUT)
    if stage == 'none':                                  # import-only (the smoke driver runs the methods itself)
        return
    if stage == 'probe':
        info = probe_all()
        os.makedirs(out, exist_ok=True)
        json.dump(info, open(os.path.join(out, 'gam_probe.json'), 'w'), indent=1, sort_keys=True)
        return
    if stage == 'freeze':
        recs = []
        for p in a['--dev'].split(','):
            recs += [json.loads(l) for l in open(p) if l.strip()]
        with open(os.path.join(out, 'gam_dev_all.jsonl'), 'w') as f:
            for r in sorted(recs, key=lambda r: (r['dataset'], r['method'], r.get('cfg_id', 0), r['seed'])):
                f.write(json.dumps(r, sort_keys=True) + '\n')
        recs = [r for r in recs if 'error' not in r and r['method'] in GRIDS]        # the gradient-descent controls of the first round are not candidates
        pick = pick_configs(recs)
        frozen = {}
        for m, cid in pick.items():
            frozen[m] = GRIDS[m][cid]
            frozen[m + '_id'] = cid
        frozen['nam_matched'] = dict(frozen['nam'])
        frozen['ours_sgd'] = dict(frozen['ours_gd'])
        json.dump(frozen, open(os.path.join(out, 'gam_frozen.json'), 'w'), indent=1, sort_keys=True)
        print(json.dumps(frozen, indent=1, sort_keys=True))
        return
    if stage == 'assemble':                              # (local) the test jsonl files and shape archives of the notebooks -> gam_raw.npz
        recs = []
        for p in a['--jsonl'].split(','):
            recs += [json.loads(l) for l in open(p) if l.strip()]
        order = {n: i for i, (n, _, _, _) in enumerate(DATASETS)}
        recs.sort(key=lambda r: (order[r['dataset']], r['seed'], METHODS.index(r['method'])))
        shapes = {}
        for p in a['--shapes'].split(','):
            z = np.load(p)
            shapes.update({k: z[k] for k in z.files})
        dev = [dict(method=r['method'], cfg_id=r['cfg_id'], dataset=r['dataset'], seed=r['seed'], score=r['score'], params=r['params'])
               for r in (json.loads(l) for l in open(a['--dev-file']) if l.strip()) if 'error' not in r]
        meta = dict(frozen=json.load(open(a['--frozen-file'])), n_records=len(recs), note=a.get('--note', ''), dev=dev)
        save_raw(os.path.join(out, 'gam_raw.npz'), recs, shapes, meta)
        print('wrote gam_raw.npz', len(recs), len(shapes))
        return
    names = a['--datasets'].split(',')
    ds_list = [d for d in DATASETS if d[0] in names]
    frozen = json.loads(a['--frozen']) if '--frozen' in a else None
    methods = a['--methods'].split(',') if '--methods' in a else None
    print('RUN', stage, [d[0] for d in ds_list], flush=True)
    run_stage(stage, ds_list, out, frozen, methods, int(a.get('--min-id', 0)))


if __name__ == '__main__':
    main(sys.argv[1:])
