"""Single-point mutants of src/core/imaging.py for what was added after tools/mutants/wp2_imaging_own.py: the BIC rule, its use in fit_frame and the greedy order, for tools/mutants.py:

    venv/bin/python tools/mutants.py --exp ../src/core/imaging --tests tests/test_imaging.py --raw results/smlm/calibration.md --md results/smlm/calibration.md \\
        --out-env MUT_OUT --spec tools/mutants/wp2_imaging_bic.py --work <scratch dir> -k "bic or greedy or fit_frame or unit_significance" [--slot a --i0 0 --step 2]

(only `tests` means anything for these mutants). A tuple is (id, what, old, new[, k]): `old` replaced by `new` at its k-th occurrence in src/core/imaging.py.
Result: 14 of 15 killed (tests/test_imaging.py, the original selection plus test_stop_bic_returns_the_unit_count_of_the_block_not_its_index,
test_fit_frame_bic_uses_its_own_parameter_count_and_the_orders_differ_in_their_times, test_greedy_order_normalises_the_correlation_by_the_weights).
G03 survives and is equivalent on any data without an exact tie of two lattice significances (a tie needs bit-identical floating-point sums of mirrored candidates).
"""
M = [
    # stop_bic
    ('B01', 'bic: no ln N in the penalty', "cost = [N * h['chi2'] + n_par * len(h['w']) * np.log(N) for h in hist]", "cost = [N * h['chi2'] + n_par * len(h['w']) for h in hist]"),
    ('B02', 'bic: log10 instead of the natural log', "cost = [N * h['chi2'] + n_par * len(h['w']) * np.log(N) for h in hist]", "cost = [N * h['chi2'] + n_par * len(h['w']) * np.log10(N) for h in hist]"),
    ('B03', 'bic: the parameters per unit ignored', "cost = [N * h['chi2'] + n_par * len(h['w']) * np.log(N) for h in hist]", "cost = [N * h['chi2'] + len(h['w']) * np.log(N) for h in hist]"),
    ('B04', 'bic: the largest cost chosen', "i = int(np.argmin(cost))\n    return i, len(hist[i]['w'])", "i = int(np.argmax(cost))\n    return i, len(hist[i]['w'])"),
    ('B05', 'bic: returns the index as the count', "i = int(np.argmin(cost))\n    return i, len(hist[i]['w'])", "i = int(np.argmin(cost))\n    return i, i"),
    ('B06', 'bic: chi2 not scaled by N', "cost = [N * h['chi2'] + n_par * len(h['w']) * np.log(N) for h in hist]", "cost = [h['chi2'] + n_par * len(h['w']) * np.log(N) for h in hist]"),
    ('B07', 'bic: the last minimum instead of the first', "i = int(np.argmin(cost))\n    return i, len(hist[i]['w'])", "i = int(len(cost) - 1 - np.argmin(cost[::-1]))\n    return i, len(hist[i]['w'])"),
    # fit_frame
    ('F03', 'fit_frame: BIC with a fixed three parameters per unit', "i = stop_bic(R['hist'], ny * nx, n_par=int(par))[0]", "i = stop_bic(R['hist'], ny * nx, n_par=3)[0]"),
    ('F04', 'fit_frame: BIC counts the samples of one row', "i = stop_bic(R['hist'], ny * nx, n_par=int(par))[0]", "i = stop_bic(R['hist'], nx, n_par=int(par))[0]"),
    ('F05', 'fit_frame: the order option is not passed to the path', "tol=tol, passes=1, rtrunc=rtrunc, r_ex=r_ex, fit_radius=fit_radius, order=order)", "tol=tol, passes=1, rtrunc=rtrunc, r_ex=r_ex, fit_radius=fit_radius)"),
    ('F06', 'fit_frame: the BIC path stops at the first insignificant entrant', "    if kind in ('snr', 'bic'):\n        kw.update(snr_min=3.0, snr_extra=snr_extra)", "    if kind in ('snr', 'bic'):\n        kw.update(snr_min=6.0, snr_extra=snr_extra)"),
    # greedy order
    ('G01', 'greedy: the least significant candidate enters', "ent = np.argwhere(z == z.max())[:1]", "ent = np.argwhere(z == z[np.isfinite(z)].min())[:1]"),
    ('G02', 'greedy: raw correlation instead of the significance', "z = np.where(grow, v / np.sqrt(np.maximum(lat.norm2(swv ** 2), 1e-300)), -np.inf)", "z = np.where(grow, v, -np.inf)"),
    ('G03', 'greedy: ties all enter together', "ent = np.argwhere(z == z.max())[:1]", "ent = np.argwhere(z == z.max())"),
    ('G04', 'greedy: the time does not advance', "ent = np.argwhere(z == z.max())[:1]\n            tau += 1.0", "ent = np.argwhere(z == z.max())[:1]\n            tau += 0.0"),
]
