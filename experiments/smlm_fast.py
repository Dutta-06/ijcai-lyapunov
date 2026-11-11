"""
The packaged SMLM frame fit of `exp_smlm_iss.run_frame`, computed by `src/core/imaging_fast.py` (the same algorithm, a cheaper arithmetic, an exact early exit)
================================================================================================================================================================

`run_frame_fast(args)` has the signature and the return value of `exp_smlm_iss.run_frame` and reads the same worker state (`exp_smlm_iss._init` fills it); the
only differences are `imaging_fast.sliding_path_fast` for `imaging.sliding_path` and `imaging_fast.refine_fast` for `imaging.refine`, and, for a stop rule of the form
('snr', kappa) alone, the early exit `stop_kappa` (the path ends at the first entrant below kappa: for every other rule, or several rules, the whole tail is kept).
Because the path stops there, the returned `taus` and `snr` end at that entrant (the detections of the 'snr' rule are those of the full path). `early=False` keeps the whole tail.
"""
import os
import sys
import time

HERE = os.path.abspath(os.path.dirname(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)
import numpy as np                          # noqa: E402
import exp_smlm_iss as ISS                  # noqa: E402
import smlm_common as C                     # noqa: E402
from src.core import imaging as IM          # noqa: E402
from src.core import imaging_fast as F      # noqa: E402


def run_frame_fast(args, early=True):
    k, frame = args
    t0 = time.perf_counter()
    W = ISS._W
    P, L, bg, o = W['P'], W['L'], W['bg'], W['opts']
    widths = o.get('widths', ISS.WIDTHS)
    nc = 3 if widths else 2
    y0 = frame - bg
    y = P.highpass(y0 - np.median(y0))
    sw0 = ISS._sw(np.zeros((0, nc)), np.zeros(0))
    rules = o['rules']
    stop_kappa = rules[0][1] if (early and len(rules) == 1 and rules[0][0] == 'snr') else None
    R = F.sliding_path_fast(P, y, sw0, L, max_entries=o.get('max_entries', ISS.MAX_ENTRIES), sw_fn=ISS._sw, sw_on='params', sig_bounds=widths, tol=o.get('tol', 1e-8),
                            max_nfev=o.get('max_nfev', 200), passes=1, rtrunc=4.0, r_ex=W['r_ex'], snr_min=o.get('snr_min', ISS.SNR_MIN), snr_extra=ISS.SNR_EXTRA,
                            jump=o.get('jump'), extra=6, order=o.get('order', 'dual'), cycle_guard=False, stop_kappa=stop_kappa)
    p = dict(hist=R['hist'], taus=R['taus'], snr=R['snr'])
    cache, out = {}, []
    for rule in rules:
        i = ISS._rule_block(rule, p)
        if i not in cache:
            h = R['hist'][i]
            cc, ww = h['c'], h['w']
            if len(ww):
                cc, ww, _ = F.refine_fast(P, y, sw0, cc, ww, sig_bounds=widths, box=1.0, sw_fn=ISS._sw, sw_on='params')
            sg = cc[:, 2] if nc == 3 else np.full(len(ww), P.sigma)
            cache[i] = np.c_[(cc[:, 0] + 0.5) * C.PIXEL_NM, (cc[:, 1] + 0.5) * C.PIXEL_NM, ww, sg] if len(ww) else np.zeros((0, 4))
        out.append(cache[i])
    return out, R['taus'], R['snr'], time.perf_counter() - t0
