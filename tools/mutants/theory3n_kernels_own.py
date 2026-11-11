# (id, what, old, new, k): single-point mutants of experiments/theory3n_kernels.py (Gram tail envelopes below the table start 0.4; finding E-2 of review/REVIEW_roundT3E.md).
# Run: tools/mutants.py --exp theory3n_kernels --tests tests/theory3n_test.py --raw results/synthetic/theory3n_raw.npz --md results/synthetic/theory3n.md --out-env T3N_OUT
#      --extra results/synthetic/theory3n.md --extra tests/theory3c_test.py --spec tools/mutants/theory3n_kernels_own.py --work <scratch dir> -k "gram_tail"
# Result (11 mutants, criterion: the tests file alone): 9 killed.  Survivors, equivalent on the data: K06 (tangential eigenvalue dropped from the Hessian norm: equal to |Gamma''| at 0 for the Gaussian and
#   dominated by the global fallback for the Bessel Grams); K10 (the global fallback without M2: the trivial bound 1.0 of the fallback exceeds M2 for every normalised Gram tested).
M = [
 ('H00', 'control (unmutated)', '', '', 0),
 ('K01', 'no correction below the table start', "    TK.Gram.tail = tail_fixed", "    pass", 0),
 ('K02', 'order 0: Gamma(0) dropped', "        return max(base, g0)", "        return base", 0),
 ('K03', 'order 0: fallback 0 for the Bessel kernels', "            g0 = 1.0                                # Bessel", "            g0 = 0.0                                # Bessel", 0),
 ('K04', 'order 1 evaluated with the second derivative', "                vals.append(TK.up(abs(self.G1(c))))", "                vals.append(TK.up(abs(self.G2(c))))", 0),
 ('K05', 'order 2 evaluated with the first derivative', "                vals.append(TK.up(abs(self.G2(c))))\n                if order == 'H'", "                vals.append(TK.up(abs(self.G1(c))))\n                if order == 'H'", 0),
 ('K06', 'Hessian: tangential eigenvalue dropped', "                    vals.append(TK.up(abs(self.ptilde_ball(a, b))))", "                    pass", 0),
 ('K07', 'cells do not reach down to l', "cells = self._cells(max(l, 0.0), L_MIN, NCELL)", "cells = self._cells(max(l, 0.2), L_MIN, NCELL)", 0),
 ('K08', 'the table value dropped', "    vals = [base]\n", "    vals = [0.0]\n", 0),
 ('K09', 'the unchanged branch starts above the table start', "    if l >= L_MIN:\n        return base", "    if l >= 2 * L_MIN:\n        return base", 0),
 ('K10', 'the global fallback loses M2', "g = [base, 1.0, self.M1(), self.M2(hess_n)]", "g = [base, 1.0, self.M1()]", 0),
 ('K11', 'a single cell', "NCELL = 16", "NCELL = 1", 0),
]
