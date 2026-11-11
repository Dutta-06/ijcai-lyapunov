"""
The sliding inverse scale space on an image: pixel-integrated Gaussian point-spread functions, a sub-pixel start lattice, the
sign-constrained block fit, and the two stopping rules (WP2 of the imaging study)
====================================================================================================================================

The setting of `rich_transit_nd` on a camera frame y ∈ ℝ^{ny×nx} (background removed) that is a sum of K point emitters seen through a
PSF of known width: y = Σ_i w_i φ(c_i) + noise, φ(c) the PSF *integrated over each pixel* (so the kernel is exactly what a pixel records:
`PixelPSF`), the coefficients w_i ≥ 0 photons, the centres c_i = (x_i, y_i) pixel units, a pixel (row i, column j) centred at
(x, y) = (j, i). The sliding path of Paper 17 in its limit form (the training of a factorised readout from a small initialisation, with
the kernels' centres sliding, converges to it; `rich_transit_nd.sliding_iss_nd`):

  * the inactive kernels of the *start lattice* (spacing h, `Lattice`) carry a dual u_j that grows at the rate v_j = ⟨W r φ_j⟩ of the
    weighted residual r of the current block and enters the block when it reaches the wall ρ (equal masses: ρ_j = ρ; a tie enters together);
  * the block relaxes to a stationary point of the weighted fit ½Σ_p W_p(y_p − Σ_i w_iφ_i(p))² under w_i ≥ 0, the centres sliding
    (`block_fit`: the coefficients first by non-negative least squares at the centres in place, then all together by a trust-region fit —
    the same polish `rich_transit_nd.polish_block_nd` uses, accepted unchanged: the unpolished flow reaches the same stationary point
    on the benchmark of Paper 17);
  * a kernel whose coefficient falls to 0 with a non-positive correlation leaves, its dual back on the wall (`exit`).

Every dual is piecewise linear in τ (the block's residual is constant between events), so the path is exact: no time stepping. All the
correlations of the whole lattice at once are one matrix product (the lattice is a tensor grid and the PSF is separable): V = G_yᵀ(W r)G_x.
The path's entry times are the quantities a stopping rule reads: the discrepancy rule (`stop_discrepancy`: the first block whose χ² is
within the noise) needs the noise level (the variance model of the camera), the largest-relative-gap rule (`stop_gap`: the block before the
entry whose time jumps most, relative to the entry before) needs none. Which of them works under which noise is the experiment's question
(`experiments/exp_smlm_spike.py`), not a property proved here: the proved statement is for noise-free data (Theorem D‴.53) and the
first-order noise response (D‴.56); with Poisson/EM noise, heteroscedastic weights and many emitters, both are heuristics until tested.
"""
import numpy as np
from scipy.optimize import least_squares, nnls
from scipy.special import erf

_SQ2 = np.sqrt(2.0)
_SQ2PI = np.sqrt(2.0 * np.pi)


def smoother(n, sb):
    """The n × n smoothing matrix S of the frame axis: row i the Gaussian of width sb centred at i, truncated at 4·sb and renormalised to
    sum 1 (so near the border it averages what is there): B = S_y ⊗ S_x smooths a frame as S_y Y S_xᵀ."""
    d = np.arange(n)
    S = np.exp(-0.5 * ((d[:, None] - d[None, :]) / sb) ** 2) * (np.abs(d[:, None] - d[None, :]) <= 4 * sb)
    return S / S.sum(1, keepdims=True)


class PixelPSF:
    """The PSF of standard deviation `sigma` (pixels) integrated over each pixel of an ny × nx frame, per axis
    Ψ(p; c) = ½[erf((p + ½ − c)/(√2σ)) − erf((p − ½ − c)/(√2σ))] — the probability that a photon of an emitter at c lands in pixel p —
    and φ(c) = Ψ_y ⊗ Ψ_x (it sums to 1 over the plane). Centres: (k, 2) = (x, y), x the column coordinate, or (k, 3) = (x, y, σ) with a
    width of its own for each kernel (the moving widths of Proposition D‴.59).
    parts = ((a_m, κ_m), …): a kernel that is a sum of separable Gaussians, φ = Σ_m a_m Ψ(κ_mσ) ⊗ Ψ(κ_mσ) (the default, one part (1, 1): a
    Gaussian; a wide second part models the halo of a measured PSF).
    hp: the width (pixels) of a high-pass: the model is H φ = φ − S_yφ_yS_xᵀ-smoothed φ (S = `smoother`, border-consistent) — to be
    compared with the frame H y (`highpass`), which removes pedestals and the broad halo: each part gets a partner −a_m(Sφ_m), also
    separable, so the lattice products and the Jacobian below keep their cost."""

    def __init__(self, ny, nx, sigma, y0=0, x0=0, parts=((1.0, 1.0),), hp=None, _parent=None):
        self.ny, self.nx, self.sigma = int(ny), int(nx), float(sigma)
        self.y0, self.x0 = int(y0), int(x0)
        self.N = self.ny * self.nx
        self.py = self.y0 + np.arange(self.ny, dtype=float)
        self.px = self.x0 + np.arange(self.nx, dtype=float)
        self.parts, self.hp = tuple((float(a), float(k)) for a, k in parts), hp
        self._parent = _parent
        if _parent is None and hp:
            self.Sy, self.Sx = smoother(self.ny, hp), smoother(self.nx, hp)

    def window(self, ya, yb, xa, xb):
        """The same PSF on the sub-frame of rows ya:yb and columns xa:xb (absolute pixel coordinates are kept; a high-pass is the full
        frame's, cropped)."""
        root = self._parent or self
        w = PixelPSF(yb - ya, xb - xa, self.sigma, self.y0 + ya, self.x0 + xa, self.parts, self.hp, _parent=root)
        w._crop = (ya + (self.y0 - root.y0), yb + (self.y0 - root.y0), xa + (self.x0 - root.x0), xb + (self.x0 - root.x0)) if self._parent is None \
            else (ya + self._crop[0], yb + self._crop[0], xa + self._crop[2], xb + self._crop[2])
        return w

    def reach(self, s=None):
        """The largest width among the kernel's terms at width s (default the PSF's): max_m κ_mσ, and √((κσ)² + hp²) for a high-passed part."""
        s = self.sigma if s is None else s
        r = max(k * s for a, k in self.parts)
        return float(np.hypot(r, self.hp) if self.hp else r)

    def highpass(self, y):
        """H y = y − S_y y S_xᵀ for a full frame (identity without a high-pass)."""
        root = self._parent or self
        return y if not self.hp else y - root.Sy @ y @ root.Sx.T

    def axis(self, p, c, s=None, deriv_s=False):
        """Ψ(p_a; c_i) and ∂Ψ/∂c_i, both (len(p), len(c)): ∂Ψ/∂c = −[exp(−u₊²) − exp(−u₋²)]/(√(2π)σ), u± = (p ± ½ − c)/(√2σ); the
        width σ is the PSF's own or, with s (len(c)), one per kernel; with deriv_s also ∂Ψ/∂σ = −[u₊exp(−u₊²) − u₋exp(−u₋²)]/(√πσ)."""
        sg = self.sigma if s is None else np.asarray(s, float)[None, :]
        up = (p[:, None] + 0.5 - c[None, :]) / (_SQ2 * sg)
        um = (p[:, None] - 0.5 - c[None, :]) / (_SQ2 * sg)
        eu, em = np.exp(-up ** 2), np.exp(-um ** 2)
        out = (0.5 * (erf(up) - erf(um)), -(eu - em) / (_SQ2PI * sg))
        if deriv_s:
            out += (-(up * eu - um * em) / (np.sqrt(np.pi) * sg),)
        return out

    def comps(self, c, deriv_s=False):
        """The separable terms of the kernels c: a list of (coefficient, gy, gx, dgy, dgx[, sgy, sgx]) with φ_i = Σ coefficient·gy_i ⊗ gx_i
        (columns i) — one term per part, and a partner per part with a high-pass; sg: ∂/∂σ of the kernel's own width."""
        c = np.atleast_2d(np.asarray(c, float))
        c = c.reshape(-1, c.shape[-1])
        root = self._parent or self
        py, px = (root.py, root.px) if self.hp else (self.py, self.px)
        s0 = c[:, 2] if c.shape[1] > 2 else np.full(len(c), self.sigma)
        out = []
        for a, kap in self.parts:
            ax = self.axis(px, c[:, 0], kap * s0, deriv_s)
            ay = self.axis(py, c[:, 1], kap * s0, deriv_s)
            fx = [ax[0], ax[1]] + ([kap * ax[2]] if deriv_s else [])
            fy = [ay[0], ay[1]] + ([kap * ay[2]] if deriv_s else [])
            terms = [(a, fy, fx)]
            if self.hp:
                terms.append((-a, [root.Sy @ f for f in fy], [root.Sx @ f for f in fx]))
            for coef, fy_, fx_ in terms:
                if self.hp and self._parent is not None:
                    ya, yb, xa, xb = self._crop
                    fy_, fx_ = [f[ya:yb] for f in fy_], [f[xa:xb] for f in fx_]
                out.append((coef, fy_[0], fx_[0], fy_[1], fx_[1]) + ((fy_[2], fx_[2]) if deriv_s else ()))
        return out

    def factors(self, c, deriv_s=False):
        """(Ψ_y (ny, k), Ψ_x (nx, k), ∂Ψ_y, ∂Ψ_x) of the Gaussian kernels c (a single part, no high-pass) and, with deriv_s, also
        (∂_σΨ_y, ∂_σΨ_x)."""
        assert len(self.parts) == 1 and not self.hp and self.parts[0] == (1.0, 1.0)
        t = self.comps(c, deriv_s)[0]
        return t[1:3] + t[3:5] + t[5:]

    def phi(self, c):
        """φ_i on the frame, one column per kernel: (N, k) (row-major pixels)."""
        return sum(a * (gy[:, None, :] * gx[None, :, :]).reshape(self.N, -1) for a, gy, gx, *_ in self.comps(c))

    def synth(self, c, w):
        """Σ_i w_iφ(c_i) as a frame (ny, nx)."""
        if np.size(w) == 0:
            return np.zeros((self.ny, self.nx))
        w = np.asarray(w, float)
        return sum(a * ((gy * w) @ gx.T) for a, gy, gx, *_ in self.comps(c))

    def corr(self, c, Rw):
        """v_i = Σ_p Rw_p φ_i(p) and its gradient in the centre (∂v/∂x_i, ∂v/∂y_i), for a weighted residual frame Rw = W∘r."""
        v = vx = vy = 0.0
        for a, gy, gx, dgy, dgx in self.comps(c):
            A, B = gy.T @ Rw, dgy.T @ Rw
            v, vx, vy = v + a * (A * gx.T).sum(1), vx + a * (A * dgx.T).sum(1), vy + a * (B * gx.T).sum(1)
        return v, vx, vy

    def jac(self, c, w, sw):
        """The Jacobian of the residual sw∘(y − Σ w_iφ_i) in z = (w, x, y) — and σ, if c has a third column —: (N, 3k) or (N, 4k)."""
        c = np.atleast_2d(np.asarray(c, float))
        wid = c.shape[1] > 2
        k = len(w)
        s = np.asarray(sw, float).reshape(-1, 1)
        W = np.asarray(w)[None, :]
        cw = cx = cy = cs = 0.0
        for a, gy, gx, dgy, dgx, *sg in self.comps(c, deriv_s=wid):
            cw = cw + a * (gy[:, None, :] * gx[None, :, :]).reshape(self.N, k)
            cx = cx + a * (gy[:, None, :] * dgx[None, :, :]).reshape(self.N, k)
            cy = cy + a * (dgy[:, None, :] * gx[None, :, :]).reshape(self.N, k)
            if wid:
                cs = cs + a * ((sg[0][:, None, :] * gx[None, :, :]) + (gy[:, None, :] * sg[1][None, :, :])).reshape(self.N, k)
        cols = [cw, cx * W, cy * W] + ([cs * W] if wid else [])
        return -s * np.concatenate(cols, axis=1)


class Lattice:
    """The start lattice: tensor grid with spacing h (pixels) over the frame and a margin around it (emitters just outside still
    shine in). Candidate (b, a) is at x = xs[a], y = ys[b]. `corr(Rw)` returns all correlations ⟨Rw φ_(b,a)⟩ as (nly, nlx) — the kernels
    at the PSF's own width, high-passed like the model if the PSF is."""

    def __init__(self, psf, h=0.5, margin=1.0):
        self.psf, self.h, self.margin = psf, float(h), float(margin)
        self.xs = np.arange(-margin, psf.nx - 1 + margin + 1e-9, h)
        self.ys = np.arange(-margin, psf.ny - 1 + margin + 1e-9, h)
        self.shape = (self.ys.size, self.xs.size)
        # lattice factors per axis: the kernel terms at the lattice abscissae (x) and ordinates (y)
        tx = psf.comps(np.c_[self.xs, np.full_like(self.xs, 0.0)])
        ty = psf.comps(np.c_[np.full_like(self.ys, 0.0), self.ys])
        self.terms = [(a, gy_, gx_) for (a, _, gx_, *_), (_, gy_, _, *__) in zip(tx, ty)]

    def corr(self, Rw):
        return sum(a * (gy.T @ Rw @ gx) for a, gy, gx in self.terms)

    def near(self, c, r):
        """Boolean mask (nly, nlx) of the lattice points closer than r to a centre of c (k, 2) = (x, y): the disks are stamped around the
        nearest lattice point of each centre, so the cost is k·(2r/h)² rather than k·nly·nlx."""
        m = np.zeros(self.shape, bool)
        if len(c) == 0:
            return m
        reach = int(np.ceil(r / self.h)) + 1
        off = np.arange(-reach, reach + 1)
        ox, oy = np.meshgrid(off, off)
        ox, oy = ox.ravel(), oy.ravel()
        ia = np.rint((c[:, 0] - self.xs[0]) / self.h).astype(int)
        ib = np.rint((c[:, 1] - self.ys[0]) / self.h).astype(int)
        A, Bi = ia[:, None] + ox[None, :], ib[:, None] + oy[None, :]
        ok = (A >= 0) & (A < self.shape[1]) & (Bi >= 0) & (Bi < self.shape[0])
        A, Bi = A[ok], Bi[ok]
        kk = np.broadcast_to(np.arange(len(c))[:, None], ok.shape)[ok]
        d2 = (self.xs[A] - c[kk, 0]) ** 2 + (self.ys[Bi] - c[kk, 1]) ** 2
        m[Bi[d2 < r * r], A[d2 < r * r]] = True
        return m

    def norm2(self, W):
        """Σ_p W_p φ_(b,a)(p)² for every candidate, (nly, nlx): the squared noise scale of the correlation when W is the inverse variance. The kernel
        is Σ_m a_m g_y,m ⊗ g_x,m, so its square is the sum over pairs of terms of (g_y,m g_y,m′) ⊗ (g_x,m g_x,m′)."""
        out = 0.0
        for a, gy, gx in self.terms:
            for a2, gy2, gx2 in self.terms:
                out = out + a * a2 * ((gy * gy2).T @ W @ (gx * gx2))
        return out

    def pos(self, b, a):
        return np.array([self.xs[a], self.ys[b]])


def block_fit(psf, y, sw, c, w, c_lo, c_hi, joint=True, max_nfev=200, tol=1e-12):
    """The block's fit from (c (k, 2 or 3), w (k)): min ½Σ(sw∘(y − Σ w_iφ_i))² over w ≥ 0 and the centres (and widths, third column) inside
    [c_lo, c_hi]. Stage 1: the coefficients at the centres in place (non-negative least squares, exact); stage 2 (joint): all together by
    the trust-region reflective solver with the analytic Jacobian (tolerance `tol` on the step, the cost and the gradient).
    Returns (c, w, info)."""
    k, nc = len(w), c.shape[1]
    if k == 0:
        return c, w, dict(nfev=0, status=0)
    swv = np.asarray(sw, float) * np.ones((psf.ny, psf.nx))
    yv = (swv * y).ravel()
    w1, _ = nnls(swv.reshape(-1, 1) * psf.phi(c), yv, maxiter=200 * k)
    if not joint:
        return c, w1, dict(nfev=0, status=0)
    lo = np.r_[np.zeros(k), c_lo.T.ravel()]
    hi = np.r_[np.full(k, np.inf), c_hi.T.ravel()]
    z0 = np.r_[w1, np.clip(c, c_lo, c_hi).T.ravel()]

    def res(z):
        return (swv * (y - psf.synth(z[k:].reshape(nc, k).T, z[:k]))).ravel()

    def jac(z):
        return psf.jac(z[k:].reshape(nc, k).T, z[:k], swv)
    sol = least_squares(res, z0, jac=jac, bounds=(lo, hi), method='trf', x_scale='jac', xtol=tol, ftol=tol, gtol=tol, max_nfev=max_nfev)
    return sol.x[k:].reshape(nc, k).T.copy(), sol.x[:k], dict(nfev=int(sol.nfev), status=int(sol.status))


def _cluster(c, seeds, dint):
    """Indices of the kernels connected to a seed through chains of centre distances < dint (kernels whose footprints overlap)."""
    k = len(c)
    inside = np.zeros(k, bool)
    inside[list(seeds)] = True
    frontier = list(seeds)
    while frontier:
        d = np.linalg.norm(c[:, None, :] - c[frontier][None, :, :], axis=2).min(1)
        new = (~inside) & (d < dint)
        inside |= new
        frontier = list(np.flatnonzero(new))
    return np.flatnonzero(inside)


def sliding_path(psf, y, sw, lat, rho=1.0, max_entries=40, box=2.0, tie_tol=1e-9, vtol=1e-12, chi_stop=0.0, joint=True,
                 sw_fn=None, max_nfev=200, rtrunc=5.5, local=True, jump=None, extra=3, sig_bounds=None, tol=1e-12, passes=2, r_ex=0.0, sw_on='model', snr_min=None, snr_extra=3, order='dual', fit_radius=None, tau_max=None):
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
    Returns dict(events = [(τ, (b, a), 'enter'/'exit')], hist = [dict(tau, n_ent, c, w, chi2, ids, nfev)] — the block after each event
    (hist[0] the empty block), taus = the entry times in order, ids = the lattice indices of the block's kernels, snr = each entrant's
    matched-filter significance at its entry, v/√(ΣW φ²) with v the correlation that took it to the wall: the number of noise standard
    deviations by which the residual correlates with the entrant's kernel)."""
    ny, nx = psf.ny, psf.nx
    swv = np.asarray(sw, float) * np.ones((ny, nx))
    U = np.zeros(lat.shape)
    taken = np.zeros(lat.shape, bool)
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
        blocked = taken
        if r_ex > 0.0 and len(w):
            blocked = taken | lat.near(c[:, :2], r_ex)
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
            S = _cluster(c, seeds, dint)
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
                cS, wS, info = block_fit(pw, y[ya:yb, xa:xb] - m_oth, swv[ya:yb, xa:xb], c[S], w[S], lo, hi,
                                         joint=joint, max_nfev=max_nfev, tol=tol)
                c[S], w[S] = cS, wS
            else:
                lo, hi = bounds(c0)
                c, w, info = block_fit(psf, y, swv, c, w, lo, hi, joint=joint, max_nfev=max_nfev, tol=tol)
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
            del ids[i]
            c, w, c0 = np.delete(c, i, 0), np.delete(w, i), np.delete(c0, i, 0)
        if gone.any():
            model = psf.synth(c, w)
        hist.append(dict(tau=tau, n_ent=n_ent, c=c.copy(), w=w.copy(), chi2=chi2(), ids=list(ids), nfev=nfev))
        if (hist[-1]['chi2'] <= chi_stop or n_ent >= max_entries or (stop_at is not None and n_ent >= stop_at)
                or (stop_snr_at is not None and n_ent >= stop_snr_at)):
            break
    return dict(events=events, hist=hist, taus=np.array(taus), snr=np.array(snrs))


def refine(psf, y, sw, c, w, sig_bounds=None, box=1.0, sw_fn=None, rounds=2, tol=1e-10, max_nfev=100, sw_on='model'):
    """Polish a block (c (k, 2 or 3), w) against the whole frame: centres within `box` pixels of where they are, widths inside
    sig_bounds, the weights re-estimated from the fitted model between the `rounds` (sw_fn(model) → sqrt weights; None keeps sw). The
    blocks of `sliding_path` are stationary up to the truncated tails and the one-event-old weights; this removes both. Returns (c, w, sw)."""
    swv = np.asarray(sw, float) * np.ones((psf.ny, psf.nx))
    c, w = np.array(c, float), np.array(w, float)
    if len(w) == 0:
        return c, w, swv
    for _ in range(rounds):
        lo, hi = c - box, c + box
        if c.shape[1] == 3:
            lo[:, 2], hi[:, 2] = sig_bounds
        c, w, _ = block_fit(psf, y, swv, c, w, lo, hi, max_nfev=max_nfev, tol=tol)
        if sw_fn is not None:
            swv = np.asarray(sw_fn(psf.synth(c, w)) if sw_on == 'model' else sw_fn(c, w), float) * np.ones((psf.ny, psf.nx))
    return c, w, swv


def fit_psf_parts(stamp, n_parts=3, sigma0=1.0):
    """Fit the PSF of an image stamp (a square array centred on the middle pixel, normalised to sum 1) by a sum of separable pixel-integrated Gaussians:
    φ = Σ_m a_m Ψ(κ_mσ) ⊗ Ψ(κ_mσ), κ_1 = 1, Σ a_m = 1. Returns (sigma, parts, rms) with parts = ((a_m, κ_m), …) for `PixelPSF(…, sigma, parts=parts)` and
    rms the root-mean-square residual relative to the stamp's peak."""
    stamp = np.asarray(stamp, float)
    stamp = stamp / stamp.sum()
    n = stamp.shape[0]
    P = PixelPSF(n, n, 1.0)
    ctr = np.array([[(n - 1) / 2.0, (n - 1) / 2.0]])

    def parts_of(q):
        a = np.r_[q[1:n_parts], 0.0]
        kap = np.r_[1.0, q[n_parts:]]
        a_ = np.r_[1.0 - np.sum(q[1:n_parts]), q[1:n_parts]]
        return tuple(zip(a_, kap))

    def model(q):
        Q = PixelPSF(n, n, q[0], parts=parts_of(q))
        return Q.phi(ctr)[:, 0].reshape(n, n)

    q0 = np.r_[sigma0, np.full(n_parts - 1, 0.15), np.geomspace(2.0, 6.0, n_parts - 1)]
    lo = np.r_[0.3, np.zeros(n_parts - 1), np.full(n_parts - 1, 1.2)]
    hi = np.r_[6.0, np.full(n_parts - 1, 0.9), np.full(n_parts - 1, 20.0)]
    sol = least_squares(lambda q: (model(q) - stamp).ravel() / stamp.max(), q0, bounds=(lo, hi), x_scale=np.r_[0.3, np.full(n_parts - 1, 0.1), np.full(n_parts - 1, 1.0)])
    q = sol.x
    return float(q[0]), parts_of(q), float(np.sqrt(np.mean(sol.fun ** 2)))


def stop_discrepancy(hist, N, z=2.0, chi_target=1.0):
    """The discrepancy rule: the first block whose degrees-of-freedom-corrected χ², Σ W r²/(N − 3k), is within z standard deviations
    (√(2/(N − 3k))) of chi_target (1 for a variance model that is right). Returns the index into hist (the last block if none)."""
    for i, h in enumerate(hist):
        k = len(h['w'])
        dof = max(N - 3 * k, 1)
        if h['chi2'] * N / dof <= chi_target * (1.0 + z * np.sqrt(2.0 / dof)):
            return i
    return len(hist) - 1


def stop_bic(hist, N, n_par=3):
    """The BIC rule, which has no constant to tune: K̂ = the number of units of the block that minimises N·χ² + n_par·k·ln N (χ² = hist's weighted mean square residual, k the units
    of the block, n_par the parameters of a unit: 3 for position and amplitude, 4 with a learned width; N the number of samples), i.e. a unit is kept if it lowers the weighted
    sum of squares by more than n_par·ln N. For a unit entering alone that is a matched-filter significance of √(n_par ln N): 5.0 to 5.8 for frames of 4000 to 10⁴ pixels. Returns
    (index into hist, K̂); the first minimum if several."""
    cost = [N * h['chi2'] + n_par * len(h['w']) * np.log(N) for h in hist]
    i = int(np.argmin(cost))
    return i, len(hist[i]['w'])


def noise_scale(lat, W):
    """The noise scale of the lattice correlations, the median over the candidates of √(Σ_p W_p φ_j(p)²) (W the inverse variance): a source of matched-filter
    significance κ correlates at κ·scale, so under the dual's unit wall it enters at τ ≈ 1/(κ·scale)."""
    return float(np.sqrt(np.median(lat.norm2(W))))


def unit_significance(psf, sw, c, w, chunk=256):
    """The own significance of each fitted unit, w_i·√(Σ_p sw_p² φ_i(p)²): its amplitude over the noise sd it would have if it stood alone (the high-passed kernel, the
    weights sw). A unit whose amplitude collapsed during the fits (a ghost) or a faint blend is small here whatever its entry time was."""
    c = np.atleast_2d(np.asarray(c, float))
    w = np.asarray(w, float)
    if len(w) == 0:
        return np.zeros(0)
    W = (np.asarray(sw, float) * np.ones((psf.ny, psf.nx))).ravel() ** 2
    out = np.empty(len(w))
    for a in range(0, len(w), chunk):
        ph = psf.phi(c[a:a + chunk])
        out[a:a + chunk] = w[a:a + chunk] * np.sqrt(W @ ph ** 2)
    return out


def stop_time(hist, taus, tau_stop):
    """The stopping time of the inverse scale space: the block after the last entrant with τ ≤ tau_stop. Returns (index into hist, K̂)."""
    K = int(np.sum(np.asarray(taus, float) <= tau_stop))
    return max(i for i, h in enumerate(hist) if h['n_ent'] <= K), K


def stop_snr(hist, taus, snr, kappa=4.5):
    """The significance rule: K̂ = the number of leading entries whose matched-filter significance is ≥ kappa (a detection threshold in noise
    standard deviations; ≈ 4.5 gives a few 10⁻³ false alarms per frame of ~10³ independent resolution cells); the block returned is the one
    just before the (K̂+1)-th entry. Returns (index into hist, K̂)."""
    snr = np.asarray(snr, float)
    bad = np.flatnonzero(snr < kappa)
    K = int(bad[0]) if bad.size else len(snr)
    return max(i for i, h in enumerate(hist) if h['n_ent'] <= K), K


def stop_gap(hist, taus, k_min=1, snr=None, gate=None):
    """The largest-relative-gap rule: with the entry times τ_1 ≤ τ_2 ≤ … , K̂ = argmax_{k ≥ k_min} τ_{k+1}/τ_k; the block returned is the one
    just before the (K̂+1)-th entry. Returns (index into hist, K̂, ratio); (0, 0, 1) when fewer than k_min + 1 kernels entered, or — with
    snr and gate — when the first entrant is not significant (snr[0] < gate: an empty frame, which a ratio cannot see)."""
    taus = np.asarray(taus, float)
    if taus.size <= k_min or (gate is not None and snr is not None and np.asarray(snr)[0] < gate):
        return 0, 0, 1.0
    ratio = taus[k_min:] / np.maximum(taus[k_min - 1:-1], 1e-300)
    j = int(np.argmax(ratio)) + k_min              # K̂
    idx = max(i for i, h in enumerate(hist) if h['n_ent'] <= j)
    return idx, j, float(ratio.max())


def simulate_frame(psf, c, photons, bg, rng, read_sd=0.0):
    """A frame with Poisson photon noise: counts ~ Poisson(bg + Σ_i photons_i φ(c_i)) per pixel, plus Gaussian read noise read_sd."""
    mean = bg + psf.synth(np.atleast_2d(c), np.asarray(photons, float))
    return rng.poisson(mean).astype(float) + read_sd * rng.standard_normal(mean.shape)


# the EM-CCD of the 2016 SMLM challenge's simulator (srm.epfl.ch dataset pages, "Parameters of simulation"): quantum efficiency 0.90,
# EM gain 300 (gamma-distributed multiplication), read-out noise 74.4 e⁻ (Gaussian), clock-induced charge 0.002 e⁻/pixel, 45 e⁻ per
# digital number (DN), baseline 100 DN, 16 bit — a total gain of QE·300/45 = 6 DN per photon. The excess noise of the EM register makes
# the variance of the photon-unit image x = (DN − baseline)/6 equal to (2/QE)·(photons at the pixel) + (read-out/45/6)².
EMCCD = dict(qe=0.90, em_gain=300.0, read_e=74.4, cic=0.002, e_per_adu=45.0, baseline=100.0)


def emccd_variance(mean_photons, cam=EMCCD):
    """The variance of the photon-unit image x = (DN − baseline)/gain at a pixel of mean `mean_photons` (photons, including
    background): (2/QE)·mean + (read-out e⁻/e⁻ per DN/gain)² [+ the clock-induced charge's, negligible]."""
    g = cam['qe'] * cam['em_gain'] / cam['e_per_adu']
    return (2.0 / cam['qe']) * np.asarray(mean_photons, float) + (cam['read_e'] / cam['e_per_adu'] / g) ** 2


def simulate_emccd(psf, c, photons, bg, rng, cam=EMCCD):
    """A frame of the challenge's EM-CCD in photon units: photoelectrons n ~ Poisson(QE(bg + Σ photons_iφ_i) + CIC), the EM register's
    output Gamma(shape n, scale gain) (0 for n = 0), Gaussian read-out noise, quantisation to integer DN, returned as
    (DN − baseline)/(QE·gain/e⁻-per-DN) — photons, plus noise."""
    mean = bg + psf.synth(np.atleast_2d(c), np.asarray(photons, float))
    n = rng.poisson(cam['qe'] * mean + cam['cic'])
    out = np.where(n > 0, rng.gamma(np.maximum(n, 1), cam['em_gain']), 0.0)
    dn = np.round(cam['baseline'] + (out + cam['read_e'] * rng.standard_normal(mean.shape)) / cam['e_per_adu'])
    return (dn - cam['baseline']) / (cam['qe'] * cam['em_gain'] / cam['e_per_adu'])


def fit_frame(frame, psf, lat, bg, var_a, var_b0=0.0, widths=None, stop=('snr', 5.0), max_entries=40, r_ex=2.0, fit_radius=None, rtrunc=4.0, refine_box=1.0,
              snr_extra=3, tol=1e-8, clean=None, order='dual'):
    """Localise the point sources of one frame by the sliding inverse scale space (the packaged algorithm).
    frame: (ny, nx) in photon-like units; bg: scalar or (ny, nx) background (subtracted; with a high-pass psf its level does not matter but it enters the variance);
    the variance of a pixel is var_a·(bg + model) + var_b0 (model = the fitted sources, unfiltered: iterated reweighting), or, if var_a is an array (ny, nx), that fixed
    variance map. psf: `PixelPSF` (parts, hp), lat: `Lattice` on it. widths = (σ_lo, σ_hi) lets every kernel slide its width, None keeps the PSF's.
    stop = ('snr', κ): the entrants of matched-filter significance ≥ κ (the sparse regime); ('time', κ): the entrants whose entry time is ≤ 1/(κ·noise scale) (the
    inverse-scale-space stopping time, for crowded fields); ('disc', z): the discrepancy rule; ('gap', gate): the largest relative gap gated by the first entrant's
    significance; ('bic', n_par): the block of the smallest N·chi^2 + n_par·k·ln N (`stop_bic`, nothing to tune; n_par = 3 for position and amplitude, 4 with learned widths).
    order = 'dual' (the accumulated dual enters first: the sliding inverse scale space) or 'greedy' (the instantaneous correlation: matching pursuit with sliding). Returns dict(c (k, 2 or 3), w (k), hist, taus, snr, index, n_entries): the refined detections of the chosen block and the path."""
    ny, nx = frame.shape
    fixed = np.ndim(var_a) == 2
    nc = 3 if widths else 2
    Pp = PixelPSF(ny, nx, psf.sigma, parts=psf.parts)
    bgv = np.asarray(bg, float) * np.ones((ny, nx))

    def sw_fn(c, w):
        if fixed:
            return 1.0 / np.sqrt(np.asarray(var_a, float) + var_b0)
        m = Pp.synth(c, w) if len(w) else 0.0
        return 1.0 / np.sqrt(var_a * (bgv + np.maximum(m, 0.0)) + var_b0)
    y0 = np.asarray(frame, float) - bgv
    y = psf.highpass(y0 - np.median(y0))
    sw0 = sw_fn(np.zeros((0, nc)), np.zeros(0))
    kind, par = stop
    kw = dict(max_entries=max_entries, sw_fn=None if fixed else sw_fn, sw_on='params', sig_bounds=widths, tol=tol, passes=1, rtrunc=rtrunc, r_ex=r_ex, fit_radius=fit_radius, order=order)
    if kind in ('snr', 'bic'):
        kw.update(snr_min=3.0, snr_extra=snr_extra)
    elif kind == 'time':
        kw.update(tau_max=1.0 / (par * noise_scale(lat, sw0 ** 2)))
    R = sliding_path(psf, y, sw0, lat, **kw)
    if kind == 'snr':
        i = stop_snr(R['hist'], R['taus'], R['snr'], par)[0]
    elif kind == 'bic':
        i = stop_bic(R['hist'], ny * nx, n_par=int(par))[0]
    elif kind == 'time':
        i = stop_time(R['hist'], R['taus'], 1.0 / (par * noise_scale(lat, sw0 ** 2)))[0]
    elif kind == 'disc':
        i = stop_discrepancy(R['hist'], ny * nx, z=par)
    elif kind == 'gap':
        i = stop_gap(R['hist'], R['taus'], snr=R['snr'], gate=par)[0]
    else:
        raise ValueError('unknown stop rule %r' % (kind,))
    h = R['hist'][i]
    c, w = h['c'], h['w']
    if len(w) and refine_box:
        c, w, _ = refine(psf, y, sw0, c, w, sig_bounds=widths, box=refine_box, sw_fn=None if fixed else sw_fn, sw_on='params')
    if clean is not None and len(w):
        keep = unit_significance(psf, sw0, c, w) >= clean
        c, w = c[keep], w[keep]
    return dict(c=c, w=w, hist=R['hist'], taus=R['taus'], snr=R['snr'], index=i, n_entries=len(R['taus']))

