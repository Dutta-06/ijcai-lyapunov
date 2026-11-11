"""
Corrected tail envelopes of the Gram below the table start (round theory3, direction N; finding E-2 of review/REVIEW_roundT3E.md).

`theory3e_kernels.Gram.tail(order, l)` returns the table value of the suffix maximum sup over [L_MIN, inf) for EVERY l <= L_MIN = 0.4, i.e. an underestimate of sup over [l, inf) whenever
the Gram's derivatives are larger on [l, L_MIN] (the entrants and their lattice neighbours are evaluated at X = 0 by the T1 certificate: Gamma(0) = 1 against 0.961 for the Gaussian,
sup |Gamma''| = 1/2 against 0.442).  The old file is not edited; this module defines the corrected function and binds it to the class (`install()`, run at import, idempotent):

    tail(order, l)  for l >= L_MIN:  unchanged (the table);
    tail(0, l)      for l <  L_MIN:  max( table , |Gamma(0)| )  (|Gamma(d)| <= Gamma(0) for an autocorrelation);
    tail(order, l)  for l <  L_MIN:  max( table value at L_MIN , sup over [l, L_MIN] of |Gamma^(order)| ),  order 1, 2, the second term by ball enclosures on 16 cells;
                                     order 'H' (Hessian norm, n >= 2): also the tangential eigenvalue |Gamma'|/d (ptilde) on the cells;
                                     if a cell enclosure is not finite (a Gram whose raw formula is singular at 0): max(table, 1, M1, M2) as the global fallback.
`tail_old` keeps the uncorrected function for comparison.
"""
import math
import os
import sys

HERE = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..'))
import theory3e_kernels as TK     # noqa: E402

L_MIN = TK.L_MIN
NCELL = 16
tail_old = TK.Gram.tail if not hasattr(TK.Gram, '_tail_old_t3n') else TK.Gram._tail_old_t3n


def tail_fixed(self, order, l, hess_n=None):
    """Upper bound of sup_{d >= l} |Gamma^{(order)}(d)| (order 0, 1, 2) or of the Hessian norm (order 'H'), valid for every l >= 0."""
    base = tail_old(self, order, l, hess_n)
    l = float(l)
    if l >= L_MIN:
        return base
    nn = self.n if hess_n is None else hess_n
    if order == 0:                                  # |Gamma(d)| <= Gamma(0) (an autocorrelation): the point value, no cell overestimate
        try:
            g0 = TK.up(abs(self.G(TK.arb(0))))
        except (ValueError, ZeroDivisionError, ArithmeticError, AssertionError):
            g0 = 1.0                                # Bessel kernels: the raw formula is singular at 0; the Gram is normalised, Gamma(0) = 1
        return max(base, g0)
    cells = self._cells(max(l, 0.0), L_MIN, NCELL)
    vals = [base]
    try:
        for a, b in cells:
            c = TK.cell(a, b)
            if order == 1:
                vals.append(TK.up(abs(self.G1(c))))
            else:
                vals.append(TK.up(abs(self.G2(c))))
                if order == 'H' and nn > 1:
                    vals.append(TK.up(abs(self.ptilde_ball(a, b))))
        if all(math.isfinite(v) for v in vals):
            return max(vals)
    except (ValueError, ZeroDivisionError, ArithmeticError, AssertionError):
        pass
    g = [base, 1.0, self.M1(), self.M2(hess_n)]
    return max(g)


def install():
    if not hasattr(TK.Gram, '_tail_old_t3n'):
        TK.Gram._tail_old_t3n = TK.Gram.tail
    TK.Gram.tail = tail_fixed
    global tail_old
    tail_old = TK.Gram._tail_old_t3n


install()
