"""
The closed loop of Proposition D‴.49 integrated in C (`closed_loop.c`, compiled on first use with the system's C compiler into a cache
directory — $LCK_BUILD or the system's temporary directory —, called through ctypes): the scalar plant ẋ = f(x) + u,
u = ẋ_m − k_ee − f̂(x̃), e = x̃ − x_m, x̃ = x + n(t) the measured state, x_m = x₀ + a sin(ωt) + a₂ sin(ω₂t), f a sum of Gaussian bumps,
f̂ = Σ_iw_iφ_i the kernels on the paradigm's free law from α = 0 with the masses b0_i — kept on their curves, β_i = b0_i + 2c sinh²(α_i/a_max)
(Q exact) — and moving centres (the paradigm's metric, Wasserstein–Fisher–Rao, or fixed), optionally with Proposition D‴.46's leak −λη_i
toward the start. Dormand–Prince 5(4); the same equations as `exp_closed_transit.real_loop`, about two orders of magnitude faster.
"""
import ctypes
import hashlib
import os
import subprocess
import tempfile

import numpy as np

_SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'closed_loop.c')
_LIB = None
_P = ctypes.POINTER(ctypes.c_double)
_I = ctypes.POINTER(ctypes.c_int)
METRICS = {'paradigm': 0, 'wfr': 1, 'fixed': 2}


def _lib():
    global _LIB
    if _LIB is None:
        with open(_SRC, 'rb') as fh:
            tag = hashlib.sha1(fh.read()).hexdigest()[:16]
        d = os.environ.get('LCK_BUILD', os.path.join(tempfile.gettempdir(), 'lck_build'))
        os.makedirs(d, exist_ok=True)
        so = os.path.join(d, f'closed_loop_{tag}.so')
        if not os.path.exists(so):
            tmp = f'{so}.tmp{os.getpid()}'
            subprocess.check_call(['cc', '-O3', '-fPIC', '-shared', '-o', tmp, _SRC, '-lm'])
            os.replace(tmp, so)
        lib = ctypes.CDLL(so)
        lib.simulate.restype = ctypes.c_long
        lib.simulate.argtypes = [_I, _P] + [_P] * 13
        _LIB = lib
    return _LIB


def _arr(a):
    a = np.ascontiguousarray(np.atleast_1d(np.asarray(a, float)))
    return a, a.ctypes.data_as(_P)


def tones(rms, w_lo, w_hi, k=16, seed=0):
    """A measurement noise of the given rms: k tones of equal amplitude at frequencies uniform in [w_lo, w_hi], random phases."""
    rng = np.random.default_rng(seed)
    return dict(w=rng.uniform(w_lo, w_hi, k), ph=rng.uniform(0, 2 * np.pi, k), amp=np.full(k, rms * np.sqrt(2.0 / k)))


def simulate(L, b0, c0, sig, bumps, ke, ref, tout, x_init=None, metric='paradigm', regx=False, noise=None, lam=0.0, th=None,
             per=None, rtol=1e-9, atol=1e-12, hmax=0.0, max_steps=10 ** 9, gate=None, etrue=False, cbox=np.inf, al0=None):
    """Integrate the closed loop. L: the free law (m, G, γ, a_max); b0: the masses; c0: the kernels' starting centres; sig: their width;
    bumps = (centres, amplitudes, widths) of f; ref = (x₀, a, ω) or (x₀, a, ω, a₂, ω₂); tout: increasing output times (tout[0] the
    start); x_init: x(0) (default x_m(0)); regx: evaluate the adaptation's regressor at x_m instead of x̃; noise: `tones(...)` or None;
    lam: the leak; th: the formation threshold on |w| (entry times recorded; default no entries); per: the window for max|e| at the start
    and the end (default the reference's period); gate: Proposition D‴.57's dead zone, dict(mode='energy' or 'formation', Tf, lo, hi,
    theta, z0) — the error's energy filtered over Tf, the drive gated smoothly between the thresholds lo and hi on its rms, every kernel
    ('energy') or those with |w| below about theta ('formation'), z0 the filter's start (default (2·hi)²: open); etrue: the adaptation
    driven by the true error x − x_m (the control and the gate by the measured one); cbox: the centres bounded to [−cbox, cbox] (referee
    round 50's tests of the mechanism); al0: the kernels' starting α (default 0). Returns dict(t, x, al, c,
    w, entries (time of the first |w_i| > th, NaN if never), e_first, e_last (max|x̃ − x_m| over the first and the last window per),
    e_true_last, e_true_rms_last (the true error x − x_m's max and rms over the last window), x_range (the state's least and largest value
    at the accepted steps), steps[, z the filtered energy])."""
    lib = _lib()
    b0, c0 = np.broadcast_to(np.asarray(b0, float), np.shape(c0)).copy(), np.asarray(c0, float)
    d = c0.size
    bc, ba, bs = (np.atleast_1d(np.asarray(v, float)) for v in bumps)
    bs = np.broadcast_to(bs, bc.shape).copy()
    ref = tuple(ref) + (0.0, 0.0) * (len(ref) == 3)
    x0, a, om, a2, om2 = (float(v) for v in ref)
    tout = np.asarray(tout, float)
    if per is None:
        per = 2 * np.pi / om
    nz = noise or dict(w=np.zeros(0), ph=np.zeros(0), amp=np.zeros(0))
    gm = {None: 0, 'energy': 1, 'formation': 2}[None if gate is None else gate['mode']]
    n = 1 + 2 * d + (gm > 0)
    al0 = np.zeros(d) if al0 is None else np.broadcast_to(np.asarray(al0, float), (d,))
    y0 = np.r_[x0 if x_init is None else float(x_init), al0, c0]                  # x_m(0) = x₀
    if gm:
        y0 = np.r_[y0, gate.get('z0', (2 * gate['hi']) ** 2)]
    ip = np.array([d, bc.size, METRICS[metric], int(bool(regx)), len(nz['w']), tout.size, gm, int(bool(etrue))], dtype=np.int32)
    dp = np.array([L.m, L.G, L.c, L.A, sig, ke, lam, x0, a, om, a2, om2, rtol, atol, hmax, np.inf if th is None else th, per,
                   float(max_steps)] + ([gate['Tf'], gate['lo'] ** 2, gate['hi'] ** 2, gate.get('theta', 0.0)] if gm else [1.0, 0.0, 1.0, 0.0])
                  + [float(cbox)])
    yout = np.zeros((tout.size, n))
    ent = np.zeros(d)
    emax = np.zeros(6)
    work = np.zeros(10 * n)
    keep = [_arr(v) for v in (bc, ba, bs, b0, nz['w'], nz['ph'], nz['amp'], y0, tout)]
    steps = lib.simulate(ip.ctypes.data_as(_I), dp.ctypes.data_as(_P), *[k[1] for k in keep],
                         yout.ctypes.data_as(_P), ent.ctypes.data_as(_P), emax.ctypes.data_as(_P), work.ctypes.data_as(_P))
    if steps < 0:
        raise RuntimeError('closed_loop.simulate: ' + ('max_steps exceeded' if steps == -1 else 'step size underflow'))
    al, c = yout[:, 1:1 + d], yout[:, 1 + d:1 + 2 * d]
    w = (b0 + 2 * L.c * np.sinh(al / L.A) ** 2) * L.A * np.tanh(al / L.A) / L.m
    out = dict(t=tout, x=yout[:, 0], al=al, c=c, w=w, entries=ent, e_first=float(emax[0]), e_last=float(emax[1]), steps=int(steps),
               e_true_last=float(emax[2]), e_true_rms_last=float(np.sqrt(emax[3] / per)), x_range=(float(emax[4]), float(emax[5])))
    if gm:
        out['z'] = yout[:, 1 + 2 * d]
    return out
