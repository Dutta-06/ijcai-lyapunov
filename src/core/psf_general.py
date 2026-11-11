"""A kernel family for frames whose PSF is not a round mixture of Gaussians (the sliding path of `imaging.py` needs the kernel only as a sum of separable terms).

`imaging.PixelPSF` models the kernel as a mixture of isotropic Gaussians that share one centre. That fits a smooth round PSF (the M2 frame's psField stamp, to 0.1 % of its
peak) and not a measured stamp that is elliptical or lumpy in the core (the M53 and M92 frames', 0.7 to 0.8 % rms, 5 to 15 % of the peak in single core pixels). `GeneralPixelPSF`
keeps everything the path relies on (a sum of separable pixel-integrated Gaussians, the border-consistent high-pass partners, the sliding centre, the window of a local fit) and lets
every part have its own width per axis and its own offset from the kernel's centre: part = (a, kx, ky, dx, dy), the weight a, the widths kx and ky in pixels and the offset (dx, dy) in
pixels. The weights sum to 1, so the kernel has unit integral like the stamp.

`fit_general_parts` fits a stamp by such parts, starting from the isotropic fit and adding one part at a time where the residual is largest; `fit_adequate` is the rule that keeps
the isotropic fit wherever it is already adequate (its rms against the stamp's peak is at most `RMS_GOAL`), so that nothing changes for a frame like M2's, and takes the richer family
only where it is not. The criterion uses the stamp alone: no catalogue, no truth.
"""
import numpy as np
from scipy.optimize import least_squares

from . import imaging as IM

RMS_GOAL = 2e-3                 # adequate: the rms residual against the stamp's peak (M2's isotropic fits: 1.0e-3; M53's and M92's: 6.5e-3 to 8.4e-3)
K_MAX = 10                      # at most this many parts
W_LO, W_HI = 0.35, 12.0         # bounds of a part's width (pixels)
D_MAX = 2.5                     # bound of a part's offset (pixels)


class GeneralPixelPSF(IM.PixelPSF):
    """`PixelPSF` with parts (a, kx, ky, dx, dy): see the module docstring. Widths are in pixels (the base width is 1), the kernels are (k, 2) = (x, y) centres (no width of its own)."""

    def __init__(self, ny, nx, parts, hp=None, y0=0, x0=0, _parent=None):
        super().__init__(ny, nx, 1.0, y0, x0, parts=((1.0, 1.0),), hp=hp, _parent=_parent)
        self.parts5 = tuple(tuple(float(v) for v in p) for p in parts)
        assert all(len(p) == 5 for p in self.parts5) and len(self.parts5) > 0
        self.parts = None                                   # the base class's parts are not used: a base method that reads them fails instead of returning another kernel

    def window(self, ya, yb, xa, xb):
        root = self._parent or self
        w = GeneralPixelPSF(yb - ya, xb - xa, self.parts5, self.hp, self.y0 + ya, self.x0 + xa, _parent=root)
        w._crop = (ya + (self.y0 - root.y0), yb + (self.y0 - root.y0), xa + (self.x0 - root.x0), xb + (self.x0 - root.x0)) if self._parent is None \
            else (ya + self._crop[0], yb + self._crop[0], xa + self._crop[2], xb + self._crop[2])
        return w

    def reach(self, s=None):
        """The largest extent of a part from the kernel's centre, max_m (max(kx, ky) + max(|dx|, |dy|)), and with a high-pass its square root of the sum of squares with hp."""
        assert s is None or s == self.sigma, 'no per-kernel width'
        r = max(max(kx, ky) + max(abs(dx), abs(dy)) for a, kx, ky, dx, dy in self.parts5)
        return float(np.hypot(r, self.hp) if self.hp else r)

    def factors(self, c, deriv_s=False):
        raise NotImplementedError('the general kernel has several separable terms')

    def comps(self, c, deriv_s=False):
        """The separable terms of the kernels c (k, 2): a list of (coefficient, gy, gx, dgy, dgx), one per part and a partner (−a, smoothed) per part with a high-pass; the part m
        of the kernel at c is the Gaussian of widths (kx, ky) at c + (dx, dy), so the derivative with respect to c is that of the part's own centre."""
        assert not deriv_s, 'no per-kernel width'
        c = np.atleast_2d(np.asarray(c, float))
        c = c.reshape(-1, c.shape[-1])
        assert c.shape[1] == 2
        root = self._parent or self
        py, px = (root.py, root.px) if self.hp else (self.py, self.px)
        out = []
        for a, kx, ky, dx, dy in self.parts5:
            ax = self.axis(px, c[:, 0] + dx if dx else c[:, 0], np.full(len(c), kx), False)
            ay = self.axis(py, c[:, 1] + dy if dy else c[:, 1], np.full(len(c), ky), False)
            fx, fy = [ax[0], ax[1]], [ay[0], ay[1]]
            terms = [(a, fy, fx)]
            if self.hp:
                terms.append((-a, [root.Sy @ f for f in fy], [root.Sx @ f for f in fx]))
            for coef, fy_, fx_ in terms:
                if self.hp and self._parent is not None:
                    ya, yb, xa, xb = self._crop
                    fy_, fx_ = [f[ya:yb] for f in fy_], [f[xa:xb] for f in fx_]
                out.append((coef, fy_[0], fx_[0], fy_[1], fx_[1]))
        return out


def isotropic_as_general(sigma, parts):
    """The parts (a, κ) of `PixelPSF(…, sigma, parts)` as general parts (a, κσ, κσ, 0, 0)."""
    return tuple((float(a), float(k) * float(sigma), float(k) * float(sigma), 0.0, 0.0) for a, k in parts)


def stamp_model(parts, n):
    """The kernel of the general parts on an n × n stamp, centred on the middle pixel, (n, n)."""
    Q = GeneralPixelPSF(n, n, parts)
    return Q.phi(np.array([[(n - 1) / 2.0, (n - 1) / 2.0]]))[:, 0].reshape(n, n)


def _normalised(q, K):
    """The parameter vector as a (K, 5) array of parts with weights |a| normalised to 1."""
    P = np.asarray(q, float).reshape(K, 5).copy()
    w = np.abs(P[:, 0])
    P[:, 0] = w / w.sum()
    return P


def _refit(stamp, P, pk):
    """Least squares on all the parts' parameters; returns (rms against the peak, parts)."""
    K = len(P)
    n = stamp.shape[0]
    q0 = np.asarray(P, float).ravel()
    lo = np.tile([0.0, W_LO, W_LO, -D_MAX, -D_MAX], K)
    hi = np.tile([1.0, W_HI, W_HI, D_MAX, D_MAX], K)
    q0 = np.clip(q0, lo + 1e-9, hi - 1e-9)
    sol = least_squares(lambda q: (stamp_model(tuple(map(tuple, _normalised(q, K))), n) - stamp).ravel() / pk, q0, bounds=(lo, hi), x_scale=np.tile([0.1, 0.3, 0.3, 0.3, 0.3], K))
    return float(np.sqrt(np.mean(sol.fun ** 2))), tuple(tuple(float(v) for v in r) for r in _normalised(sol.x, K))


def fit_general_parts(stamp, start, rms_goal=RMS_GOAL, k_max=K_MAX):
    """Fit a unit-sum square stamp (centred on its middle pixel) by general parts, greedily. `start` is the list of parts to begin with (for example `isotropic_as_general` of the
    isotropic fit); all the parts' parameters are refitted; a part is added at the pixel of the largest absolute residual (weights 0.01, 0.03, 0.08 tried, width 0.8 px) as long
    as the rms residual against the stamp's peak exceeds `rms_goal`, there are fewer than `k_max` parts and the new part lowers it. Returns (rms, parts)."""
    stamp = np.asarray(stamp, float)
    stamp = stamp / stamp.sum()
    n = stamp.shape[0]
    pk = float(stamp.max())
    rms, P = _refit(stamp, [list(p) for p in start], pk)
    while rms > rms_goal and len(P) < k_max:
        R = (stamp_model(P, n) - stamp) / pk
        iy, ix = np.unravel_index(np.argmax(np.abs(R)), R.shape)
        c = (n - 1) / 2.0
        best = None
        for w0 in (0.01, 0.03, 0.08):
            r2, P2 = _refit(stamp, [list(p) for p in P] + [[w0, 0.8, 0.8, ix - c, iy - c]], pk)
            if best is None or r2 < best[0]:
                best = (r2, P2)
        if best[0] >= rms - 1e-6:
            break
        rms, P = best
    return rms, P


def fit_adequate(stamp, rms_goal=RMS_GOAL, k_max=K_MAX):
    """The kernel of a frame from its PSF stamp: the isotropic mixture of `imaging.fit_psf_parts` if its rms residual against the stamp's peak is at most `rms_goal`, else the general
    parts fitted from it. Returns dict(kind 'isotropic' | 'general', sigma, parts (the isotropic (a, κ) parts, or None), parts5 (the general parts, or None), rms_isotropic, rms (of the
    kernel returned))."""
    sigma, parts, rms0 = IM.fit_psf_parts(stamp, n_parts=3, sigma0=1.0)
    if rms0 <= rms_goal:
        return dict(kind='isotropic', sigma=float(sigma), parts=tuple(parts), parts5=None, rms_isotropic=float(rms0), rms=float(rms0))
    rms, P = fit_general_parts(stamp, isotropic_as_general(sigma, parts), rms_goal, k_max)
    return dict(kind='general', sigma=float(sigma), parts=None, parts5=P, rms_isotropic=float(rms0), rms=float(rms))


def make_psf(fit, ny, nx, hp=None):
    """The `PixelPSF` (isotropic) or `GeneralPixelPSF` of a `fit_adequate` result on an ny × nx frame with a high-pass of width hp (pixels) or none."""
    if fit['kind'] == 'isotropic':
        return IM.PixelPSF(ny, nx, fit['sigma'], parts=fit['parts'], hp=hp or None)
    return GeneralPixelPSF(ny, nx, fit['parts5'], hp=hp or None)
