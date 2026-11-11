"""
A guarded copy of `imaging.sliding_path`: the enter/exit limit cycle of the entry rule
====================================================================================================================================

Defect (found by the referee of the repeat-frame test, round REP, finding 2; reproduced on repeat patch m2r3 (110, 501), a priori default
(sigma_b 5, r_ex 1.2), where it stalls the path after 155 distinct entries at tau = 2.18969 in 12 s of run time): the dual of a lattice kernel
that is not in the block grows with the correlation of the residual at the kernel's *lattice start* (`v` in `sliding_path`), but a unit of the
block *exits* when its coefficient is 0 and the correlation at its *slid centre* is not positive (`g`). The block fit of a crowded cluster can
return the entrant with w = 0 (its centre pushed to the edge of its box) although the correlation at its lattice start is positive; the exit
puts the dual back on the wall, so the kernel re-enters at dt = 0 at the same entry time, is refitted to w = 0 and exits again, and so on
until `max_entries` is reached. The time rule (`tau_max`) never acts, because tau does not advance: a run of the astronomy experiments with
the stopping time of 12 to 48 stops after about 150 entries whatever its configuration says, and its `max_entries` flag (the cap) is not a
censoring of a long path but the signature of the cycle. The other units of the cluster are refitted a little further at each turn of the
cycle, so the cycle also acts as extra block-fit iterations on that cluster.

The guard: a lattice point that has exited twice at the same entry time is barred from entering for the rest of the path (`barred` in the
result lists them). A first exit followed by a successful re-entry is untouched, so a path without a double exit at one entry time is the
unguarded path bit for bit (`cycle_guard=False` is the unguarded path in every case; a test compares both with `imaging.sliding_path`).
The change is confined to the lines that mention `cycle_guard`, `barred`, `last_exit` or `bars`, the function's name, the `IM.` prefix of `block_fit` and
`_cluster`, and the one line that adds the exclusion zone to `blocked` (a test compares the source with `imaging.sliding_path`). `imaging.py` itself is hash-frozen by
the astronomy freeze records and is not modified: the experiments that re-run the affected paths call this function (`experiments/exp_astro_cyclefix.py`).

Limit: the rule keys on the entry time. In `order='greedy'` tau is the step count and advances by 1 at every entry, so no two exits share an entry time and the guard
never acts (the path is the unguarded path, `barred` is empty); the astronomy paths and the SMLM frames use `order='dual'`. Only repeated entries at one entry time are
removed: a cycle whose turns advance tau, or a point that exits and re-enters at different times, is not touched.
"""
import numpy as np

from . import imaging as IM


def sliding_path_guarded(psf, y, sw, lat, rho=1.0, max_entries=40, box=2.0, tie_tol=1e-9, vtol=1e-12, chi_stop=0.0, joint=True,
                 sw_fn=None, max_nfev=200, rtrunc=5.5, local=True, jump=None, extra=3, sig_bounds=None, tol=1e-12, passes=2, r_ex=0.0, sw_on='model', snr_min=None, snr_extra=3, order='dual', fit_radius=None, tau_max=None, cycle_guard=True):
    """The sliding inverse scale space of the frame y (weights W = sw²) on the lattice `lat`, for non-negative amplitudes. Starts with the
    empty block; each step the first lattice kernel to reach the wall enters (a tie, within tie_tol, together), the block is refitted
    (`block_fit`, centres within `box` pixels of their lattice starts), kernels at 0 with a non-positive correlation leave; stops when
    max_entries have entered, when no inactive kernel is growing (all correlations ≤ vtol times the first step's largest), or when
    the weighted mean square residual falls below chi_stop. sw_fn(model frame) → new sqrt weights, applied after each block fit
    (iteratively reweighted: variance from the fitted model); None keeps sw.
    local: refit only the cluster of kernels that interact with the entrant (centres closer than 2·rtrunc·σ, transitively), on the window
    around them (rtrunc·σ beyond the outermost centre), against the residual of the others — the block fit's exact stationary point up to
    the PSF's truncated tail (≈ 1e-6 of its peak at rtrunc = 5.5); the kernels outside the cluster are at stationary points already.
    jump, extra: a stop for the cost of the noise-fitting tail, not a rule: with ratios r_k = τ_{k+1}/τ_k (K̂ = k), once the largest ratio so
    far, at K̂ = k*, is ≥ jump the path stops `extra` entries after k* (a later, larger ratio would be missed: `stop_gap` of the stopped
    path is the gap rule restricted to the entries seen); None runs on to max_entries.
    sig_bounds = (σ_lo, σ_hi): every kernel carries its own width, entering at the PSF's σ and sliding inside the bounds with its centre
    (the moving widths of Proposition D‴.59); the block's c then has three columns (x, y, σ). None: all kernels have the PSF's width.
    snr_min, snr_extra: stop `snr_extra` entries after the first entrant whose matched-filter significance is below snr_min — the cost
    of the noise-fitting tail again, for the significance rule (`stop_snr`) and a few entries of tail beyond it (None: off).
    tau_max: stop before the first entry after this time (the inverse-scale-space stopping time: an entrant at τ has a time-averaged correlation of ρ/τ, so
    τ_max = ρ/(κ·scale) keeps the entrants of average significance ≥ κ, `noise_scale`; unlike the entrant's instantaneous significance it is not fooled by candidates
    that accumulated their dual before a neighbour was fitted).
    fit_radius: refit only the kernels within this many pixels of an entrant (not the transitive cluster): the local slides for crowded fields, where
    every kernel is connected to every other; the rest stay where they are and enter the residual, so the block is stationary only up to their
    response (None: the transitive cluster of `local`).
    order: 'dual' (the sliding inverse scale space: entry by the integral of the correlation) or 'greedy' (matching pursuit with sliding: the candidate
    of the largest matched-filter significance v/√ΣWφ² enters at each step, τ is the step count — the same block fits, exits and stop rules).
    sw_on: 'model' calls sw_fn(model frame), 'params' calls sw_fn(c, w) (the variance of a high-passed fit depends on the unfiltered signal).
    r_ex: the resolution limit — lattice kernels within r_ex pixels of an active kernel's centre do not grow (their duals stay): the
    dictionary's coherence control, no second emitter closer than the resolution is admitted (0: off).
    tol, passes: the block fit's tolerance, and the number of (fit, reweight) rounds per event (1 for data: the weights follow the fit
    one event behind and `refine` polishes the block that is finally used).
    cycle_guard: a lattice point that has exited twice at the same entry time (an enter/exit limit cycle: its entry is decided on the correlation at its lattice start, its exit on the correlation at its slid centre, so the re-entry can follow at dt = 0 for ever) is barred from entering for the rest of the path; False: the path of `imaging.sliding_path`, bit for bit.
    Returns dict(events = [(τ, (b, a), 'enter'/'exit')], hist = [dict(tau, n_ent, c, w, chi2, ids, nfev)] — the block after each event
    (hist[0] the empty block), taus = the entry times in order, ids = the lattice indices of the block's kernels, snr = each entrant's
    matched-filter significance at its entry, v/√(ΣW φ²) with v the correlation that took it to the wall: the number of noise standard
    deviations by which the residual correlates with the entrant's kernel), barred = [(τ, (b, a))], the lattice points barred by the guard, each with the entry time of its second exit
    (empty if the guard never acted or cycle_guard is False)."""
    ny, nx = psf.ny, psf.nx
    swv = np.asarray(sw, float) * np.ones((ny, nx))
    U = np.zeros(lat.shape)
    taken = np.zeros(lat.shape, bool)
    barred = np.zeros(lat.shape, bool)
    last_exit, bars = {}, []
    nc = 2 if sig_bounds is None else 3
    ids, c, w, c0 = [], np.zeros((0, nc)), np.zeros(0), np.zeros((0, nc))
    tau, n_ent, events, hist, taus, snrs = 0.0, 0, [], [], [], []

    def bounds(c0_):
        lo, hi = c0_ - box, c0_ + box
        if nc == 3:
            lo[:, 2], hi[:, 2] = sig_bounds
        return lo, hi

    def start(b, a):
        return np.r_[lat.pos(b, a), psf.sigma] if nc == 3 else lat.pos(b, a)
    reach = psf.reach(psf.sigma if sig_bounds is None else sig_bounds[1])
    dint = 2.0 * rtrunc * reach
    model = np.zeros((ny, nx))

    def chi2():
        return float(np.mean((swv * (y - model)) ** 2))
    hist.append(dict(tau=0.0, n_ent=0, c=c.copy(), w=w.copy(), chi2=chi2(), ids=[], nfev=0))
    vscale, stop_at, best, stop_snr_at = None, None, (0.0, 0), None
    for _ in range(11 * max_entries):
        Rw = swv ** 2 * (y - model)
        V = lat.corr(Rw)
        if vscale is None:
            vscale = float(np.abs(V).max())
            if vscale <= 0.0:
                break
        blocked = taken | barred if cycle_guard else taken
        if r_ex > 0.0 and len(w):
            blocked = blocked | lat.near(c[:, :2], r_ex)
        v = np.where(blocked, 0.0, V)
        grow = v > vtol * vscale
        if not grow.any():
            break
        if order == 'greedy':
            z = np.where(grow, v / np.sqrt(np.maximum(lat.norm2(swv ** 2), 1e-300)), -np.inf)
            ent = np.argwhere(z == z.max())[:1]
            tau += 1.0
        else:
            dt = np.full(lat.shape, np.inf)
            dt[grow] = (rho - U[grow]) / v[grow]
            dmin = float(dt.min())
            if tau_max is not None and tau + dmin > tau_max:
                break
            ent = np.argwhere(dt <= dmin * (1.0 + tie_tol) + tie_tol * 1e-300)
            U[~blocked] += dmin * v[~blocked]
            tau += dmin
        seeds = []
        for b, a in ent:
            U[b, a] = rho
            taken[b, a] = True
            seeds.append(len(ids))
            ids.append((int(b), int(a)))
            c0 = np.vstack([c0, start(b, a)[None, :]])
            c = np.vstack([c, start(b, a)[None, :]])
            w = np.r_[w, 0.0]
            events.append((tau, (int(b), int(a)), 'enter'))
            ph_ = psf.phi(start(b, a)[None, :])[:, 0]
            snrs.append(float(v[b, a] / np.sqrt(np.sum(swv.ravel() ** 2 * ph_ ** 2))))
            if snr_min is not None and snrs[-1] < snr_min and stop_snr_at is None:
                stop_snr_at = n_ent + 1 + snr_extra
            if taus and tau / taus[-1] > best[0]:
                best = (tau / taus[-1], n_ent)
                stop_at = None if (jump is None or best[0] < jump) else n_ent + 1 + extra
            taus.append(tau)
            n_ent += 1
        if not local:
            S = np.arange(len(w))
        elif fit_radius is not None:
            S = np.flatnonzero(np.linalg.norm(c[:, None, :2] - c[seeds][None, :, :2], axis=2).min(1) < fit_radius)
        else:
            S = IM._cluster(c, seeds, dint)
        nfev = 0
        for _pass in range(1 if sw_fn is None else passes):
            if local:
                r = rtrunc * reach
                ya, yb = max(int(np.floor(c[S, 1].min() - r)), 0), min(int(np.ceil(c[S, 1].max() + r)) + 1, ny)
                xa, xb = max(int(np.floor(c[S, 0].min() - r)), 0), min(int(np.ceil(c[S, 0].max() + r)) + 1, nx)
                oth = np.setdiff1d(np.arange(len(w)), S)
                m_oth = psf.synth(c[oth], w[oth])[ya:yb, xa:xb]
                pw = psf.window(ya, yb, xa, xb)
                lo, hi = bounds(c0[S])
                cS, wS, info = IM.block_fit(pw, y[ya:yb, xa:xb] - m_oth, swv[ya:yb, xa:xb], c[S], w[S], lo, hi,
                                         joint=joint, max_nfev=max_nfev, tol=tol)
                c[S], w[S] = cS, wS
            else:
                lo, hi = bounds(c0)
                c, w, info = IM.block_fit(psf, y, swv, c, w, lo, hi, joint=joint, max_nfev=max_nfev, tol=tol)
            nfev += info['nfev']
            model = psf.synth(c, w)
            if sw_fn is not None:
                swv = np.asarray(sw_fn(model) if sw_on == 'model' else sw_fn(c, w), float) * np.ones((ny, nx))
        # a kernel of the cluster at the bound w = 0 whose correlation does not push it back in leaves
        Rw = swv ** 2 * (y - model)
        g = psf.corr(c, Rw)[0]
        gone = (w <= 1e-12 * max(w.max(), 1e-300)) & (g <= vtol * vscale)
        for i in np.flatnonzero(gone)[::-1]:
            events.append((tau, ids[i], 'exit'))
            taken[ids[i]] = False
            U[ids[i]] = rho
            if cycle_guard:
                if abs(tau - last_exit.get(ids[i], -np.inf)) <= 1e-12 * max(1.0, abs(tau)):
                    barred[ids[i]] = True
                    bars.append((tau, ids[i]))
                last_exit[ids[i]] = tau
            del ids[i]
            c, w, c0 = np.delete(c, i, 0), np.delete(w, i), np.delete(c0, i, 0)
        if gone.any():
            model = psf.synth(c, w)
        hist.append(dict(tau=tau, n_ent=n_ent, c=c.copy(), w=w.copy(), chi2=chi2(), ids=list(ids), nfev=nfev))
        if (hist[-1]['chi2'] <= chi_stop or n_ent >= max_entries or (stop_at is not None and n_ent >= stop_at)
                or (stop_snr_at is not None and n_ent >= stop_snr_at)):
            break
    return dict(events=events, hist=hist, taus=np.array(taus), snr=np.array(snrs), barred=bars)
