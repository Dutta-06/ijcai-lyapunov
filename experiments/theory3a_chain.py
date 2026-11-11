"""
`theory_t1b_rigorous.verify` with the displacement accounting of Lemma 5.2 made sound (a finding of round theory3 A).

Lemma 5.2 bounds the displacement of a dormant kernel by  |c_l - c_l^0| <= mu_c ( int_0^t H^(st(t')) B(t') dt' + peak_l [ ... ] )  with the STAGE-RESOLVED H^(st(t')) (st(t') = the number of limit entries
before t').  `exp_theory_t1.cert` and `theory_t1b_rigorous.verify` use  mu_c ( H^(st(t)) I(t) + hinc )  with I(t) = int_0^t B the cumulative mass and H^(st(t)) the CURRENT stage's H only.
Since H^(st) decreases with the stage (the residual of the fit shrinks), the mass accumulated in an earlier stage is then weighted with a smaller H than in the lemma: the displacement bound dcb
(and, through  int H^(st) dcb dt, the drift term of E_j) is underestimated whenever a kernel has accumulated mass in an earlier stage.  Here the same source is patched at import (three anchored
replacements, each of which must occur exactly once) so that the displacement integral is  W = sum over the earlier sub-intervals of H^(stage) x (their mass)  +  H^(st) x (mass so far in the
current sub-interval):  dcb = mu_c (W + hinc).  Nothing else changes.  `verify_fixed` has the signature of `verify`.
"""
import inspect
import os
import sys

HERE = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..'))
import theory_t1b_rigorous as RG   # noqa: E402

PATCH = [
    ("        Icum, dsum = arb(0), arb(0)\n",
     "        Icum, dsum, Wprev = arb(0), arb(0), arb(0)\n"),
    ("                    dcb = muc * (Hst[st] * Ig + hinc)\n",
     "                    dcb = muc * (Wprev + Hst[st] * (Ig - Icum) + hinc)\n"),
    ("                Icum = up(last)\n",
     "                Wprev = up(Wprev + Hst[st] * (up(last) - Icum))\n                Icum = up(last)\n"),
]


def build():
    src = inspect.getsource(RG.verify)
    for old, new in PATCH:
        if src.count(old) != 1:
            raise RuntimeError('anchor %r found %d times' % (old, src.count(old)))
        src = src.replace(old, new)
    ns = dict(RG.__dict__)
    exec(compile(src, 'theory3a_chain:verify_fixed', 'exec'), ns)
    return ns['verify']


verify_fixed = build()
