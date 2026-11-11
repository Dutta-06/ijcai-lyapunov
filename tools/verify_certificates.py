#!/usr/bin/env python3
"""Standalone verifier of the computer-assisted certificates D and VA (paper 18's supplement, manuscript/paper18/supplement.tex, Appendix
A; the paper's Appendix A until referee round 148) from stored data.

Until referee round 148 the proof of paper 18's Theorem 2 (the transform's limit on |η| < 2) used two certificates beyond the unit disc
and off the axes — D, the domination condition on 1 < |η| ≤ 1.8927 (Proposition A.3), and VA at ε₀ = 0.11 ("VA011"), the valley
condition on 1.89 ≤ |η| < 2 (Proposition A.13). Since round 148 that theorem is proved by hand (the valley condition on 1 < |η| < 2
along a three-segment path: paper 18's Proposition 4.8, review/RESEARCH_note_annulus.md, Theorem A), and D and VA are an independent
computer-assisted cross-check, which no proof uses. They are checked here from the data in results/certificates/ alone, without
re-running the programs that proposed the boxes and paths (referee rounds 135 and 136). This file imports nothing from experiments/:
every inequality is re-implemented from the lemmas' statements in Appendix A, and every covering statement is re-derived from the stored
boxes themselves. This program, not the proposing programs' certify step, is the check of record for the cross-check.

Usage:  venv/bin/python -I tools/verify_certificates.py [--dir results/certificates] [--only D,VA] [--d-sample N] [--quiet]
Exit status 0 when every check passes (with --d-sample N > 1 the run is partial and says so), 1 on any failure, 2 on unreadable
or missing data. One line per lemma: [ok] or [FAIL], the counts and the largest rigorous bound; one [info] line that is not a check.

Labels. The output and this file use Appendix A's numbers (unchanged by its move to the supplement). The drafts (review/drafts/D61.md, Appendix A) and the review rounds up
to 135 used other names; the correspondence is:
  paper 18                                    drafts (D61)
  Theorem 2 (the transform's limit)           Theorem 3
  Lemma A.1  (the supremum over x)            Lemma D1
  Lemma A.2  (the box bound)                  Lemma D2
  Proposition A.3 (certificate D)             "certificate D"
  Lemma A.4  (the forms f_A, f_B)             Lemma VA1
  Lemma A.5  (series and tails)               Lemma VA2
  Lemma A.6  (closed forms at ϱ = 0)          Lemma VA3
  Lemma A.7  (the mean value form)            Lemma VA4
  Lemma A.8  (domains and cover)              Lemma VA5
  Lemma A.9  (the outer piece)                Lemma VA6
  Lemma A.10 (connector and handoff)          Lemma VA7
  Lemma A.11 (tiling)                         Lemma VA8
  Lemma A.12 (the end piece)                  Lemma VA9 ("the end piece")
  Proposition A.13 (certificate VA)           "certificate VA"
  equation (A.1) (the bound E)                (A.1)

What is checked (in this order; the first check's failure stops the run):
  SHA256SUMS  every data file read below is listed in SHA256SUMS and matches its hash.
  meta        the parameters are the certificates' (D: [1, 1.8927] × [0, π/2]; VA: ε₀ = 0.11, δ₀ = 0.05, END_R = 0.02, q_A = 2.5,
              p_B = 0.5, the junction's u = 4), R1 = fl(2 − ε₀), ε₀/4 < δ₀ (Lemma A.8 (i)), q_A·p_B > 1 and q_A ≥ 1 (Lemma A.8
              (ii), (iii)).
  A.2         on every stored leaf R1 ≤ R ≤ R2, f1 ≤ φ ≤ f2, Lemma A.2's bound H = max_{w ∈ W} Ḡ(w) is < 0 (Arb, 120 bits; the
              candidates w_m, the kink and the quadratics' roots enclosed by the stable formula; the top edge at the true π/2).
  A.3         D's tiling: the leaves lie in [1, 1.8927] × [0, fl(π/2)], have positive area and pairwise disjoint interiors (a sweep
              over R with the active φ-intervals kept sorted), and their areas sum to the rectangle's exactly (dyadic rationals as
              Python integers). A finite union of closed boxes in a closed rectangle with disjoint interiors and the same area is
              the rectangle; with the top edge read as the true π/2, the quadrant 1 ≤ |η| ≤ 1.8927, 0 ≤ arg η ≤ π/2 is covered.
  A.4–A.7     on every stored box of the inner pieces A-J, A-FC, B-J, B-F and on every branch the box needs (decided here, Lemma
              A.8): the mean value bound of Lemma A.7 for the form f_A or f_B (Lemma A.4) — the series S, C_v, Λ, H₁ with their
              tails (Lemma A.5; the tails computed here in Arb, so with no 1 + 10⁻⁷ factor, and floored at 10⁻³⁰⁰ as in the lemma),
              the closed forms at ϱ = 0 (Lemma A.6) in outward-rounded endpoint interval arithmetic, the derivative by forward-mode
              dual numbers — is < 0 (Arb, 128 bits).
  A.8         the branches: form A's box is evaluated on F where α > b̂ is proved on the box, on C where α < b̂, else on both;
              form B's boxes must have b̂p < α_B (the floor) proved. The junction's boxes have u = 4 exactly. The cover: the
              stored boxes cover each piece's closed parameter domain — A-FC: 0 ≤ u ≤ 4, 0 ≤ q ≤ q_A, 0 ≤ ϱ ≤ δ₀, uϱ ≤ ε₀; A-J:
              0 ≤ q ≤ q_A, 0 ≤ λ ≤ 1, 0 ≤ ϱ ≤ ε₀/4; B-F: 0 ≤ u ≤ 4, 0 ≤ p ≤ p_B, 0 ≤ ϑ ≤ (π/2)², pϑ ≤ δ₀, upϑ ≤ ε₀; B-J: the same
              at u = 4 —, proved by a recursive decomposition of the domain's bounding box at the stored boxes' faces until every
              piece lies inside one stored box or has an interior disjoint from the domain (exact rationals; (π/2)² in Arb). The
              domains are closures of their interiors (monotone constraints), so pieces whose interior misses them add nothing.
              The analytic parts of Lemma A.8 — (i) X ∈ [x_b, δ₀] gives u ≤ 4, (ii) qp = 1 sends q > q_A to p < p_B, (iii) the cap
              has p ≥ 1 and lies in form A — are proved in Appendix A and rest on the constants checked under meta. The proposer's
              root grids, clips and float bisection are not used.
  A.9         every stored segment of every η-box (the connector's pieces and the polyline's) has Lemma A.9's mean value bound
              E(M; η_c) + sup|∂_sE| + r_R sup|∂_RE| + r_θ sup|∂_θE| < 0 over the η-box (Arb, 128 bits; the top edge the true π/2).
  A.10        per η-box: lo ≤ hi finite; the connector's pieces lie on Re z = δ₀ and chain δ₀ + i·lo to δ₀ + i·hi; lo ≤ y_c ≤ hi
              and the polyline's first vertex is δ₀ + i·y_c; the polyline's pieces chain its first vertex to its last through
              every stored vertex; |1 − z_n| < END_R in balls; the stored sub-boxes tile the η-box (as in D's tiling) and on each
              Y(η) = min(κ(η)δ₀, θ) ∈ [lo, hi] is proved (κ = ω²sin 2θ/((1 − ω²)₊ + 2ω²sin²θ), (sin 2θ)₊, 0 ≤ Y ≤ θ where the
              denominator's ball meets 0, the top edge the true π/2).
  A.11        the η-boxes tile [R1, 2] × [0, fl(π/2)] (as in D's tiling), R1 = fl(2 − ε₀).
  A.12        the end piece: log ε_e + log(2 + ε_e) + (R2²/4)(1 + ε_e)² + H((R2/2)(2 + ε_e)) < 0 at R2 = 2, ε_e = END_R, with
              (R2/2)(2 + ε_e) > 1.
  seam        (under meta) 2 − ε₀ and R1 at most the real 1.89 ≤ 1.8927, in rationals: VA011's annulus 1.89 ≤ |η| < 2 starts
              inside D's disc (Proposition A.13).
Informational, not a check (it never changes the exit status): the forms f_A and f_B as implemented here, against the definition
(A.1) of E, E = log|1 − z²| + Re(η²z²/4) + Λ(ν), at sample points of the path (f_A on F, C and J; f_B on F and J with its dropped
term added back), the largest relative deviation printed. The forms are re-typed from Lemma A.4's statement (and so share its
algebra with the proposer); this line is a sanity check of the transcription, not a proof of the algebra.

Not checked here (analytic, proved in Appendix A or the paper): Lemma A.1 and the algebra of Lemmas A.4 and A.6 (that the forms and closed forms
are what they claim), Lemma A.8 (i)–(iii) beyond their constants, the real segment z ∈ [−1, x_b] of the path, the reflections
η ↦ η̄, −η, and that the valley and domination conditions give the transform's limit (Lemma 4.6 of the paper; Theorem 2 itself no
longer needs D or VA).

The trusted base: (1) Arb's ball arithmetic through python-flint (0.9.0 here): +, −, ×, ÷, √, log, sin, cos, exp, |·|, π on
real and complex balls enclose every value on the input balls, comparisons of balls are true only when true for every pair of
points, .lower()/.upper() return exact endpoints; (2) python-flint's exact conversion of an IEEE-754 double to a ball (arb(x),
acb(x, y)); (3) Python's float.fromhex (exact) for the JSON data and the struct/array module's binary64 reading for D's leaves,
CPython's binary64 comparisons (exact), fl(π)/2 (an exact halving: the stored top edge) and fl(2 − ε₀) (rounded; only compared
with the stored R1, which the seam then uses exactly), math.nextafter for a double above (π/2)² (then proved ≥ (π/2)² in Arb);
(4) Python's integers and fractions.Fraction (exact rationals) for the areas, the cover and the seams; (5) hashlib's SHA-256, gzip
and json for reading the data; (6) the CPython interpreter and this file. No floating-point arithmetic decides anything here except
through (2)–(3); every float written below is either data or a printed summary (the informational line's sample points are floats
fed exactly into balls).
"""
import argparse
import array
import bisect
import gzip
import hashlib
import heapq
import json
import math
import os
import random
import sys
import time
from fractions import Fraction as Fr

from flint import acb, arb, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_DIR = os.path.normpath(os.path.join(HERE, '..', 'results', 'certificates'))
FILES = ('meta.json', 'D_leaves.f64le.gz', 'VA011_inner.jsonl.gz', 'VA011_outer.jsonl.gz')
THEOREM = dict(D_R2=1.8927, eps0=0.11, delta0=0.05, end_r=0.02, QA=2.5, PB=0.5, UJ=4.0)   # the certificates' parameters (paper 18's Appendix A)
PI2 = math.pi / 2                       # the double fl(π)/2 (exact halving): the stored top edge, read as the true π/2
D_PREC, VA_PREC = 120, 128
TAIL_FLOOR = 1e-300                     # Lemma A.5's floor of the tails' radius
FX = float.fromhex


class Fail(Exception):
    """A box or piece outside the lemmas' domain (a raise is a failure, never a pass)."""


def hull(a, b):
    return arb(a).union(arb(b))


def pos(x):
    """x₊ on a ball."""
    if x >= 0:
        return x
    if x <= 0:
        return arb(0)
    return hull(0, x.upper())


def ball_top(f1, f2):
    """The angle interval [f1, f2] as a ball, its top edge the true π/2 where it is the double fl(π)/2."""
    return hull(f1, arb.pi() / 2) if f2 == PI2 else hull(f1, f2)


# ======================================================================================================== reading the data
def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def check_sums(d):
    """Every file in FILES listed in SHA256SUMS with a matching hash. Returns (ok, message)."""
    p = os.path.join(d, 'SHA256SUMS')
    if not os.path.exists(p):
        return False, 'SHA256SUMS missing'
    listed = {}
    for ln in open(p).read().splitlines():
        if ln.strip():
            h, name = ln.split(None, 1)
            listed[name.strip().lstrip('*')] = h.lower()
    missing = [f for f in FILES if f not in listed]
    if missing:
        return False, f'not listed in SHA256SUMS: {missing}'
    absent = [f for f in FILES if not os.path.exists(os.path.join(d, f))]
    if absent:
        return False, f'missing files: {absent}'
    bad = [f for f in FILES if sha256(os.path.join(d, f)) != listed[f]]
    if bad:
        return False, f'hash mismatch: {bad}'
    return True, f'{len(FILES)} files match'


def load(d):
    """The stored data: meta (dict), D's leaves [(R1, R2, f1, f2, H)], the inner boxes and the outer η-boxes (floats exact)."""
    meta = json.load(open(os.path.join(d, 'meta.json')))
    raw = gzip.open(os.path.join(d, 'D_leaves.f64le.gz'), 'rb').read()                # the hashed file, whatever meta says
    a = array.array('d')
    a.frombytes(raw)
    if sys.byteorder != 'little':
        a.byteswap()
    if len(a) % 5:
        raise ValueError('D_leaves: not a multiple of 5 doubles')
    leaves = [tuple(a[i:i + 5]) for i in range(0, len(a), 5)]
    inner = []
    for ln in gzip.open(os.path.join(d, 'VA011_inner.jsonl.gz'), 'rt'):
        r = json.loads(ln)
        inner.append(dict(piece=r['piece'], iv={k: (FX(v[0]), FX(v[1])) for k, v in r['iv'].items()},
                          recs=[(br, route, FX(v)) for br, route, v in r['recs']]))
    outer = []
    cz = lambda p: complex(FX(p[0]), FX(p[1]))
    for ln in gzip.open(os.path.join(d, 'VA011_outer.jsonl.gz'), 'rt'):
        r = json.loads(ln)
        outer.append(dict(box=tuple(FX(x) for x in r['box']), lo=FX(r['lo']), hi=FX(r['hi']), yc=FX(r['yc']),
                          pts=[cz(p) for p in r['pts']], conn=[(cz(s[:2]), cz(s[2:])) for s in r['conn']],
                          poly=[(cz(s[:2]), cz(s[2:])) for s in r['poly']], sub=[tuple(FX(x) for x in s) for s in r['sub']]))
    return meta, leaves, inner, outer


# ============================================================================================== tiling, in exact rationals
def tiles(rect, boxes):
    """Whether the closed boxes (x1, x2, y1, y2) (doubles) tile the closed rectangle rect: each inside it with positive area,
    pairwise disjoint interiors (a sweep over x: the boxes whose x-ranges contain the sweep position overlap in x, so their
    y-intervals must be disjoint; they are kept sorted by their lower ends, and a new interval can only meet its predecessor),
    and the areas summing to rect's exactly (dyadic rationals scaled to integers). Returns (ok, reason)."""
    X1, X2, Y1, Y2 = rect
    for i, (x1, x2, y1, y2) in enumerate(boxes):
        if not (X1 <= x1 < x2 <= X2 and Y1 <= y1 < y2 <= Y2):
            return False, f'box {i} not inside the rectangle with positive area'
    vals = [v for b in boxes for v in b[:4]] + list(rect)
    K = max(v.as_integer_ratio()[1].bit_length() - 1 for v in vals)
    sc = lambda v: (lambda n, dd: n * ((1 << K) // dd))(*v.as_integer_ratio())
    area = sum((sc(b[1]) - sc(b[0])) * (sc(b[3]) - sc(b[2])) for b in boxes)
    if area != (sc(X2) - sc(X1)) * (sc(Y2) - sc(Y1)):
        return False, f'areas sum to {float(Fr(area, 1 << (2 * K)))}, not {float(Fr(sc(X2) - sc(X1), 1 << K) * Fr(sc(Y2) - sc(Y1), 1 << K))}'
    order = sorted(range(len(boxes)), key=lambda i: boxes[i][0])
    keys, ends, heap = [], [], []
    for i in order:
        x1, x2, y1, y2 = boxes[i]
        while heap and heap[0][0] <= x1:
            _, k = heapq.heappop(heap)
            j = bisect.bisect_left(keys, k)
            assert keys[j] == k
            del keys[j], ends[j]
        j = bisect.bisect_left(keys, y2)
        if j > 0 and ends[j - 1] > y1:
            return False, f'box {i} overlaps another'
        keys.insert(j, y1)
        ends.insert(j, y2)
        heapq.heappush(heap, (x2, y1))
    return True, 'tiled'


# ============================================================================================================ D: Lemma A.2
def d_bound(R1, R2, f1, f2):
    """Lemma A.2's bound H on the box R1 ≤ R ≤ R2, f1 ≤ φ ≤ f2 ⊂ [1, 2) × [0, π/2] as an arb upper end, or None (box fails).
    Ḡ(w) = log w + (p_lo/2)(w² − ½) + M(a_hi(1 − w)) − ½ − log R1 + log 4, M(a) = a² (a ≤ 1), 2a − 1 (a > 1); W: w_m = 1 − 1/R2,
    the kink 1 − 1/a_hi and the roots in (0, w_m] of (p_lo + 2a_hi²)w² − 2a_hi²w + 1 and p_lo w² − 2a_hi w + 1."""
    half = arb(1) / 2
    Ra, Rb, fa = arb(R1), arb(R2), arb(f1)
    fb = arb.pi() / 2 if f2 == PI2 else arb(f2)
    p = None
    for R in (Ra, Rb):
        for f in (fa, fb):
            v = (R * R * (2 * f).cos()).lower()
            p = v if p is None or v < p else p
    a = (Rb * fa.cos()).upper()
    wm = (1 - 1 / Rb).upper()
    if not (wm * wm < half):
        return None
    W = [wm]
    if a > 1:
        k = 1 - 1 / a
        if not (k.upper() <= 0) and not (k.lower() > wm):
            W.append(k)
    for qa, qb in ((p + 2 * a * a, -2 * a * a), (p, -2 * a)):           # qa w² + qb w + 1
        if not (qb < 0):
            return None
        disc = qb * qb - 4 * qa
        if disc < 0:
            continue                                                    # no real root
        if not (disc >= 0):
            return None
        q = (disc.sqrt() - qb) / 2                                      # > 0; the roots 1/q and q/qa
        roots = [1 / q]
        if qa > 0 or qa < 0:
            roots.append(q / qa)
        elif not (q.lower() > abs(qa).upper() * wm):                   # qa's sign undecided: q/qa must lie beyond w_m
            return None
        W += [r for r in roots if not (r.upper() <= 0) and not (r.lower() > wm)]
    best = None
    c = -half - Ra.log() + arb(4).log()
    for w in W:
        if not (w > 0):
            return None
        s = (a * (1 - w)).upper()
        M = s * s if s <= 1 else 2 * s - 1
        v = (w.log() + p / 2 * (w * w - half) + M + c).upper()
        best = v if best is None or v > best else best
    return best


def verify_D(meta, leaves, sample=1, log=print):
    ctx.prec = D_PREC
    m = meta['D']
    rect = (FX(m['R1']), FX(m['R2']), FX(m['f1']), FX(m['f2']))
    ok = True
    t0 = time.time()
    if not (rect == (1.0, THEOREM['D_R2'], 0.0, PI2) and m['leaves'] == len(leaves) and len(leaves) > 0):
        log(f"[FAIL] D meta: rectangle {rect} or count {m['leaves']} ≠ {len(leaves)} leaves")
        return False
    nf, worst, dev, n = 0, -math.inf, 0.0, 0
    for i in range(0, len(leaves), sample):
        R1, R2, f1, f2, H = leaves[i]
        n += 1
        b = d_bound(R1, R2, f1, f2) if (1.0 <= R1 < R2 < 2.0 and 0.0 <= f1 < f2 <= PI2) else None
        if b is None or not (b < 0):
            nf += 1
            continue
        worst = max(worst, float(b))
        dev = max(dev, abs(float(b) - H))
    part = '' if sample == 1 else f' (subsample: every {sample}th leaf)'
    good = nf == 0 and n > 0
    log(f"[{'ok' if good else 'FAIL'}] Lemma A.2 (box bound < 0): {n} leaves{part}, {nf} failed, largest bound {worst:.5g} "
        f"(max |Δ| against the proposer's bound {dev:.1g}) [{time.time() - t0:.0f}s]")
    ok &= good
    t0 = time.time()
    good, why = tiles(rect, [x[:4] for x in leaves])
    log(f"[{'ok' if good else 'FAIL'}] Proposition A.3 (D's tiling): {len(leaves)} leaves, [1, {rect[1]}] × [0, π/2]: {why} (exact area, disjoint "
        f"interiors) [{time.time() - t0:.0f}s]")
    return ok and good


# ===================================================================================== VA: dual numbers, intervals, series
class Du:
    """A dual number over arb: value v, derivative d with respect to the expansion variable."""
    __slots__ = ('v', 'd')

    def __init__(s, v, d=0):
        s.v = v if isinstance(v, arb) else arb(v)
        s.d = d if isinstance(d, arb) else arb(d)

    def __add__(s, o):
        return Du(s.v + o.v, s.d + o.d) if isinstance(o, Du) else Du(s.v + o, s.d)
    __radd__ = __add__

    def __neg__(s):
        return Du(-s.v, -s.d)

    def __sub__(s, o):
        return s + (-o)

    def __rsub__(s, o):
        return (-s) + o

    def __mul__(s, o):
        return Du(s.v * o.v, s.d * o.v + s.v * o.d) if isinstance(o, Du) else Du(s.v * o, s.d * o)
    __rmul__ = __mul__

    def __truediv__(s, o):
        if isinstance(o, Du):
            nz(o.v)
            return Du(s.v / o.v, (s.d * o.v - s.v * o.d) / (o.v * o.v))
        nz(o)
        return Du(s.v / o, s.d / o)

    def __rtruediv__(s, o):
        nz(s.v)
        return Du(o / s.v, -o * s.d / (s.v * s.v))


def nz(x):
    if not (x > 0 or x < 0):
        raise Fail('division by a ball meeting 0')


def div(a, b):
    if isinstance(b, Du) or isinstance(a, Du):
        return (a if isinstance(a, Du) else Du(a)) / b
    nz(b)
    return a / b


def val(x):
    return x.v if isinstance(x, Du) else x


def sqrt_(x):
    if not (val(x) > 0):
        raise Fail('square root of a ball not > 0')
    if isinstance(x, Du):
        r = x.v.sqrt()
        return Du(r, x.d / (2 * r))
    return x.sqrt()


def pos_(x):
    """x₊; on duals with derivative 1{x > 0}x′, enclosed by the hull of 0 and x′ where x's ball meets 0 (Lemma A.7)."""
    if not isinstance(x, Du):
        return pos(x)
    if x.v >= 0:
        return x
    if x.v <= 0:
        return Du(0, 0)
    return Du(hull(0, x.v.upper()), arb(0).union(x.d))


def series(x, c, K, t0, t1):
    """Σ_{k<K} c[k] x^k by Horner in Arb plus a ball of radius t0(r) (value) / t1(r) (derivative), r = |x|'s upper end, the
    tails computed in Arb (upper ends) and floored at TAIL_FLOOR (Lemma A.5)."""
    xv = val(x)
    r = abs(xv).upper()
    if not (r >= 0):
        raise Fail('series argument not finite')
    rad = lambda t: (t if t > TAIL_FLOOR else arb(TAIL_FLOOR))
    e0, e1 = rad(t0(r).upper()), rad(t1(r).upper())
    if not (e0 >= 0 and e1 >= 0):
        raise Fail('series tail not finite')
    p = arb(0)
    for ck in reversed(c):
        p = p * xv + ck
    p = p + hull(-e0, e0)
    if not isinstance(x, Du):
        return p
    dp = arb(0)
    for k in range(K - 1, 0, -1):
        dp = dp * xv + k * c[k]
    dp = dp + hull(-e1, e1)
    return Du(p, dp * x.d)


def _ratio_ok(t):
    if not (t < 1):
        raise Fail('series tail ratio not < 1')
    return t


ctx.prec = VA_PREC                                     # the series' coefficients as balls at VA's precision (Lemma A.5)
_F = [arb(math.factorial(k)) for k in range(80)]
KS, KC, KL, KH = 16, 16, 30, 26                        # the terms of S, C_v, Λ, H₁ (Lemma A.5)
CS = [arb((-1) ** k) / _F[2 * k + 1] for k in range(KS)]
CC = [arb((-1) ** k) / _F[2 * k + 2] for k in range(KC)]
CL = [arb((-1) ** k) / (k + 2) for k in range(KL)]
CH = [arb((-1) ** (k + 1)) / (k + 3) for k in range(KH)]


def S(w):
    """sin√w/√w = Σ (−w)^k/(2k + 1)!."""
    K = KS
    return series(w, CS, K, lambda r: r ** K / _F[2 * K + 1] / (1 - _ratio_ok(r / ((2 * K + 2) * (2 * K + 3)))),
                  lambda r: K * r ** (K - 1) / _F[2 * K + 1] / (1 - _ratio_ok((K + 1) * r / (K * (2 * K + 2) * (2 * K + 3)))))


def Cv(w):
    """(1 − cos√w)/w = Σ (−w)^k/(2k + 2)!."""
    K = KC
    return series(w, CC, K, lambda r: r ** K / _F[2 * K + 2] / (1 - _ratio_ok(r / ((2 * K + 3) * (2 * K + 4)))),
                  lambda r: K * r ** (K - 1) / _F[2 * K + 2] / (1 - _ratio_ok((K + 1) * r / (K * (2 * K + 3) * (2 * K + 4)))))


def Lam(v):
    """(v − log(1 + v))/v² = Σ (−v)^k/(k + 2), |v| ≤ ½."""
    if not (abs(val(v)).upper() <= arb(1) / 2):
        raise Fail('|v| > 1/2')
    K = KL
    return series(v, CL, K, lambda r: r ** K / ((K + 2) * (1 - r)), lambda r: r ** (K - 1) / (1 - r))


def H1(d):
    """(h(1 + δ) − δ²)/δ³ = Σ (−1)^{k+1} δ^k/(k + 3), |δ| ≤ ½."""
    if not (abs(val(d)).upper() <= arb(1) / 2):
        raise Fail('|δ| > 1/2')
    K = KH
    return series(d, CH, K, lambda r: r ** K / ((K + 3) * (1 - r)), lambda r: r ** (K - 1) / (1 - r))


class Iv:
    """Endpoint interval arithmetic with arb endpoints, every result's ends rounded outward (.lower()/.upper())."""
    __slots__ = ('lo', 'hi')

    def __init__(s, lo, hi=None):
        lo = arb(lo)
        hi = lo if hi is None else arb(hi)
        s.lo, s.hi = lo.lower(), hi.upper()

    @staticmethod
    def of(o):
        return o if isinstance(o, Iv) else Iv(o)

    def __add__(s, o):
        o = Iv.of(o)
        return Iv((s.lo + o.lo).lower(), (s.hi + o.hi).upper())
    __radd__ = __add__

    def __neg__(s):
        return Iv(-s.hi, -s.lo)

    def __sub__(s, o):
        return s + (-Iv.of(o))

    def __rsub__(s, o):
        return Iv.of(o) + (-s)

    def __mul__(s, o):
        o = Iv.of(o)
        ps = [s.lo * o.lo, s.lo * o.hi, s.hi * o.lo, s.hi * o.hi]
        lo, hi = ps[0].lower(), ps[0].upper()
        for x in ps[1:]:
            lo = x.lower() if x.lower() < lo else lo
            hi = x.upper() if x.upper() > hi else hi
        return Iv(lo, hi)
    __rmul__ = __mul__

    def __truediv__(s, o):
        o = Iv.of(o)
        if not (o.lo > 0 or o.hi < 0):
            raise Fail('interval division by an interval meeting 0')
        return s * Iv((1 / o.hi).lower(), (1 / o.lo).upper())

    def __rtruediv__(s, o):
        return Iv.of(o) / s

    def sq(s):
        if s.lo >= 0:
            return Iv((s.lo * s.lo).lower(), (s.hi * s.hi).upper())
        if s.hi <= 0:
            return Iv((s.hi * s.hi).lower(), (s.lo * s.lo).upper())
        m = s.hi if s.hi > -s.lo else -s.lo
        return Iv(0, (m * m).upper())

    def pw(s, n):
        r = Iv(1)
        for _ in range(n):
            r = r * s
        return Iv(r.lo if r.lo > 0 else arb(0), r.hi) if n % 2 == 0 else r

    def pos(s):
        return Iv(s.lo if s.lo > 0 else arb(0), s.hi if s.hi > 0 else arb(0))


# ======================================================================================================= A.4: the forms
def fA(br, u, q, rho, lam=None):
    """E/X² in (u, q, ϱ) = (ε/X, θ²/X, X) on the branch br ∈ {J, F, C} (Lemma A.4's f_A; (Re σ − 1)₊ for (|Re σ| − 1)₊)."""
    w = q * rho
    eps = u * rho
    om = 1 - eps / 2
    om2 = om * om
    S1, S4, Cw = S(w), S(4 * w), Cv(w)
    alpha = u * (1 - eps / 4) + 2 * om2 * q * S1 * S1
    bh = 2 * om2 * S4
    if br == 'F':
        m = div(bh, alpha)
    elif br == 'J':
        m = div(lam * bh, alpha)
    else:
        m = Du(1) if isinstance(rho, Du) else arb(1)
    Q = q * m * (alpha * m - 2 * bh) - alpha * rho
    zz = rho + q * m * m
    vr = -2 * (rho - q * m * m) + rho * zz * zz
    L = zz * zz / 2 - vr * vr * Lam(rho * vr) / 2
    a1z = sqrt_(1 + rho * (2 + zz))
    Dh = div(om2 * (2 + zz) - u * (1 - eps / 4), om * a1z + 1)
    Ph = -u / 2 + om * (1 - (1 + rho) * q * Cw - q * m * S1)
    Dp, Pp = pos_(Dh), pos_(Ph)
    return Q + L + (Dp * Dp * (1 + rho * Dp * H1(rho * Dp)) - Pp * Pp)


def fB(br, u, p, th2, lam=None):
    """E/X² in (u, p, ϑ) = (ε/X, X/θ², θ²) on br ∈ {J, F} (Lemma A.4's f_B, the last term dropped)."""
    X = p * th2
    eps = u * X
    om = 1 - eps / 2
    om2 = om * om
    S1, S4 = S(th2), S(4 * th2)
    alpha = u * p * (1 - eps / 4) + 2 * om2 * S1 * S1
    bh = 2 * om2 * S4
    mt = div(bh, alpha) if br == 'F' else div(lam * bh, alpha)
    Q = mt * (alpha * mt - 2 * bh) - th2 * alpha
    zz = p * (th2 + mt * mt)
    vr = -2 * p * (th2 - mt * mt) + X * zz * zz
    L = zz * zz / 2 - vr * vr * Lam(X * vr) / 2
    a1z = sqrt_(1 + X * (2 + zz))
    Dh = div(om2 * (2 + zz) - u * (1 - eps / 4), om * a1z + 1)
    Dp = pos_(Dh)
    return Q + L + Dp * Dp * (1 + X * Dp * H1(X * Dp))


# ============================================================================================ A.6: the closed forms at 0
def _smaller(a, b):
    return a if a.hi < b.hi else b


def fA0(br, u, q, lam=None):
    """An upper bound of f_A at ϱ = 0 on the box (Lemma A.6), u, q, λ endpoint intervals; s = u + 2q."""
    s = u + 2 * q
    h = arb(1) / 2
    if br == 'C':
        P0 = 1 - u * h - q * (3 * h)
        if P0.lo >= 0:
            return -(u * q) - q.sq() * h                               # D₀₊² − P₀₊² = 2q(2 − u − q) exactly
        D0 = 1 - u * h + q * h
        b1 = q * (s - 4) - q.sq() * h + D0.pos().sq() - Iv(P0.pos().lo).sq()
        b2 = q * (s - 4 + 4 * D0.pos()) - q.sq() * h
        return _smaller(b1, b2)
    if br == 'F':
        D0 = 1 - u * h + 2 * q / s.sq()
        b1 = -(4 * q / s) - 8 * q.sq() / s.pw(4) + D0.pos().sq()
        b2 = -8 * q.sq() / s.pw(4) + q * (-4 / s + 2 * D0.pos() * (2 / s.sq() + h + 2 / s))
        return _smaller(b1, b2)
    D0 = 1 - u * h + 2 * q * lam.sq() / s.sq()
    return -((4 * q / s) * (lam * (2 - lam))) - 8 * lam.pw(4) * q.sq() / s.pw(4) + D0.pos().sq()


def fB0(br, u, p, lam=None):
    """An upper bound of f_B at ϑ = 0 (f_A0's floor and junction at q = 1/p, −P₀₊² dropped); s = up + 2."""
    s = u * p + 2
    h = arb(1) / 2
    if br == 'F':
        return -(4 / s) - 8 * p.sq() / s.pw(4) + (1 - u * h + 2 * p / s.sq()).pos().sq()
    return -((4 / s) * (lam * (2 - lam))) - 8 * lam.pw(4) * p.sq() / s.pw(4) + (1 - u * h + 2 * p * lam.sq() / s.sq()).pos().sq()


# ===================================================================================================== A.7: mean value
def mean_value(f, f0, br, ball, ivs, lo, hi):
    """Lemma A.7 for f(br, *ball, x[, λ]) on the expansion variable's [lo, hi]: returns (ok, route, bound or sup f′)."""
    try:
        der = f(br, *ball[:2], Du(hull(lo, hi), 1), *ball[2:]).d.upper()
        if der > 0:
            dp = der
        elif der <= 0:
            dp = arb(0)
        else:
            return False, 'x', math.inf
        v0 = f0(br, *ivs).hi
        if lo == 0:
            if v0 <= 0 and der < 0:
                return True, 'd', float(der)
            b = v0 + hi * dp
            return bool(b < 0), 'v', float(b)
        fl = f(br, *ball[:2], arb(lo), *ball[2:]).upper()
        d0 = f(br, *ball[:2], Du(hull(0, lo), 1), *ball[2:]).d.upper()
        fl2 = v0 + lo * d0
        if fl2 < fl:
            fl = fl2
        b = fl + (arb(hi) - arb(lo)) * dp
        return bool(b < 0), 'v', float(b)
    except (Fail, ZeroDivisionError):
        return False, 'x', math.inf


def branches_A(u, q, r):
    """Lemma A.8 (the branches): the branches the form-A box needs over ϱ ∈ [r_lo, r_hi]: F where α > b̂ (κX < θ), C where α < b̂, else both."""
    rho = hull(*r)
    w = q * rho
    eps = u * rho
    om = 1 - eps / 2
    d = u * (1 - eps / 4) + 2 * om * om * q * S(w) * S(w) - 2 * om * om * S(4 * w)
    return ('F',) if d > 0 else (('C',) if d < 0 else ('F', 'C'))


def floor_B(u, p, t):
    """Lemma A.8 (the branches): b̂p < α_B (the floor) over the form-B box, ϑ ∈ [t_lo, t_hi]."""
    th2 = hull(*t)
    X = p * th2
    eps = u * X
    om = 1 - eps / 2
    return bool(2 * om * om * S(4 * th2) * p < u * p * (1 - eps / 4) + 2 * om * om * S(th2) * S(th2))


PIECES = ('A-J', 'A-FC', 'B-J', 'B-F')
COORDS = {'A-J': ('q', 'l', 'r'), 'A-FC': ('u', 'q', 'r'), 'B-J': ('p', 'l', 't'), 'B-F': ('u', 'p', 't')}


def check_inner_box(name, iv, UJ):
    """Lemmas A.4–A.7 and Lemma A.8's branches on one stored box: (ok, [(branch, route, value)], reason)."""
    form, kind = name.split('-')
    sec, ev = ('q', 'r') if form == 'A' else ('p', 't')
    need = ('u', sec, ev) + (('l',) if kind == 'J' else ())
    if sorted(iv) != sorted(need):
        return False, [], 'coordinates'
    for k in need:
        if not (0.0 <= iv[k][0] <= iv[k][1]):
            return False, [], f'interval {k}'
    if kind == 'J' and iv['u'] != (UJ, UJ):
        return False, [], 'junction u ≠ 4'
    bu, bs = hull(*iv['u']), hull(*iv[sec])
    ball = [bu, bs] + ([hull(*iv['l'])] if kind == 'J' else [])
    ivs = [Iv(*iv['u']), Iv(*iv[sec])] + ([Iv(*iv['l'])] if kind == 'J' else [])
    lo, hi = iv[ev]
    if form == 'A':
        brs = ('J',) if kind == 'J' else branches_A(bu, bs, iv['r'])
        f, f0 = fA, fA0
    else:
        if kind == 'F' and not floor_B(bu, bs, iv['t']):
            return False, [], 'floor not proved'
        brs = (kind,)
        f, f0 = fB, fB0
    out = []
    for br in brs:
        ok, route, v = mean_value(f, f0, br, ball, ivs, lo, hi)
        if not ok:
            return False, out, f'branch {br}'
        out.append((br, route, v))
    return True, out, ''


def cover(rect, boxes, interior_misses):
    """Whether the closed boxes cover D ∩ rect, D a closed domain that is the closure of its interior: rect decomposed at the
    boxes' faces until each piece lies in one box or its interior misses D (interior_misses). All coordinates Fractions.
    Returns (ok, pieces in boxes, pieces outside D, a failing piece or None)."""
    dim = len(rect)
    stack = [(tuple(rect), list(range(len(boxes))))]
    nin, nout = 0, 0
    while stack:
        B, cand = stack.pop()
        if interior_misses(B):
            nout += 1
            continue
        c2 = [i for i in cand if all(boxes[i][k][0] < B[k][1] and B[k][0] < boxes[i][k][1] for k in range(dim))]
        if any(all(boxes[i][k][0] <= B[k][0] and B[k][1] <= boxes[i][k][1] for k in range(dim)) for i in c2):
            nin += 1
            continue
        if not c2:
            return False, nin, nout, B
        cut = None
        for i in c2:
            for k in range(dim):
                for v in boxes[i][k]:
                    if B[k][0] < v < B[k][1]:
                        cut = (k, v)
                        break
                if cut:
                    break
            if cut:
                break
        k, v = cut
        lo_, hi_ = list(B), list(B)
        lo_[k], hi_[k] = (B[k][0], v), (v, B[k][1])
        stack += [(tuple(hi_), c2), (tuple(lo_), c2)]
    return True, nin, nout, None


def verify_inner(meta, inner, log=print):
    """Lemmas A.4–A.8 on the stored inner boxes. Returns ok."""
    ctx.prec = VA_PREC
    m = meta['VA011']
    e0, d0, QA, PB, UJ = (FX(m[k]) for k in ('eps0', 'delta0', 'QA', 'PB', 'UJ'))
    E0, D0 = Fr(e0), Fr(d0)
    pi2sq = arb.pi() ** 2 / 4
    T = math.nextafter(float(pi2sq.upper()), math.inf)                 # a double ≥ (π/2)²: the domain's bounding box in ϑ
    ok_all = bool(arb(T) >= pi2sq)                                     # (proved, not assumed)
    by = {k: [] for k in PIECES}
    for r in inner:
        if r['piece'] not in by:
            log(f"[FAIL] inner pieces: unknown piece {r['piece']}")
            return False
        by[r['piece']].append(r)
    for name in PIECES:
        t0 = time.time()
        rs = by[name]
        nf, wv, wd, nd, nv, dev, why = 0, -math.inf, -math.inf, 0, 0, 0.0, ''
        for r in rs:
            ok, out, reason = check_inner_box(name, r['iv'], UJ)
            if not ok:
                nf += 1
                why = why or reason
                continue
            st = {br: v for br, _, v in r['recs']}
            for br, route, v in out:
                if route == 'd':
                    nd, wd = nd + 1, max(wd, v)
                else:
                    nv, wv = nv + 1, max(wv, v)
                    if br in st:
                        dev = max(dev, abs(v - st[br]))
        good = nf == 0 and len(rs) > 0 and len(rs) == m['pieces'][name]['boxes']
        log(f"[{'ok' if good else 'FAIL'}] Lemmas A.4–A.7 {name}: {len(rs)} boxes, {nf} failed{' (' + why + ')' if why else ''}; route d "
            f"(f₀ ≤ 0, sup f′ < 0) {nd}, largest sup f′ {wd:.3g}; route v {nv}, largest bound {wv:.4g} (max |Δ| against the "
            f"proposer's {dev:.1g}) [{time.time() - t0:.1f}s]")
        ok_all &= good
        # Lemma A.8: the cover of the piece's domain by its stored boxes
        t0 = time.time()
        cs = COORDS[name]
        boxes = [tuple((Fr(r['iv'][k][0]), Fr(r['iv'][k][1])) for k in cs) for r in rs]
        if name == 'A-FC':
            rect = [(Fr(0), Fr(UJ)), (Fr(0), Fr(QA)), (Fr(0), D0)]
            miss = lambda B: B[0][0] * B[2][0] >= E0                   # uϱ ≤ ε₀ (ϱ ≤ δ₀ is the bounding box)
        elif name == 'A-J':
            rect = [(Fr(0), Fr(QA)), (Fr(0), Fr(1)), (Fr(0), E0 / 4)]
            miss = lambda B: False                                     # ϱ = t = ε/4 ≤ ε₀/4 is the bounding box
        else:
            def big_t(t):                                              # ϑ ≥ (π/2)², t a face (a double), in Arb
                if Fr(float(t)) != t:
                    raise Fail('a face that is not a double')
                return bool(arb(float(t)) >= pi2sq)
            if name == 'B-F':
                rect = [(Fr(0), Fr(UJ)), (Fr(0), Fr(PB)), (Fr(0), Fr(T))]
                miss = lambda B: B[1][0] * B[2][0] >= D0 or B[0][0] * B[1][0] * B[2][0] >= E0 or big_t(B[2][0])
            else:
                uB = Fr(UJ)
                rect = [(Fr(0), Fr(PB)), (Fr(0), Fr(1)), (Fr(0), Fr(T))]
                miss = lambda B: B[0][0] * B[2][0] >= D0 or uB * B[0][0] * B[2][0] >= E0 or big_t(B[2][0])
        good, nin, nout, bad = cover(rect, boxes, miss)
        log(f"[{'ok' if good else 'FAIL'}] Lemma A.8 {name} (branches above; cover): the domain covered by its {len(boxes)} boxes "
            f"({nin} pieces inside boxes, {nout} outside the domain){'' if good else ', uncovered piece ' + str([tuple(map(float, x)) for x in bad])}"
            f" [{time.time() - t0:.1f}s]")
        ok_all &= good
    return ok_all


# ================================================================================================ A.9–A.12: the outer piece
def E_mid(z, eta):
    """E(z; η) = log|1 − z²| + Re(cz²) + C(σ) on balls, C ≤ h₊ where |σ|'s ball meets 1."""
    c = eta * eta / 4
    sg = eta / 2 * (1 + z)
    r = abs(sg)
    h = (r * r - 1) / 2 - r.log() if r > 0 else None
    if r > 1:
        a = pos(abs(sg.real) - 1)
        C = h - a * a
    elif r <= 1:
        C = arb(0)
    else:
        if h is None:
            raise Fail('|σ| meets 0')
        C = hull(0, pos(h).upper())
    w = abs(1 - z * z)
    if not (w > 0):
        raise Fail('|1 − z²| meets 0')
    return w.log() + (c * z * z).real + C


def seg_bound(P, Q, box):
    """Lemma A.9's bound of E over the segment [P, Q] × the η-box (R1, R2, f1, f2): an arb, or None when undefined."""
    R1, R2, f1, f2 = box
    Pz, Qz = acb(P.real, P.imag), acb(Q.real, Q.imag)
    M, H = (Pz + Qz) / 2, (Qz - Pz) / 2
    Rc, Tc = (arb(R1) + arb(R2)) / 2, (arb(f1) + arb(f2)) / 2
    Tb = ball_top(f1, f2)
    rR = (arb(R2) - arb(R1)) / 2
    a1, a2 = (Tb.upper() - Tc).upper(), (Tc - arb(f1)).upper()
    rT = a1 if a1 > a2 else a2
    try:
        E0 = E_mid(M, Rc * acb(0, Tc).exp())
    except Fail:
        return None
    z = M + arb(0, 1) * H                                              # s ∈ [−1, 1]
    Rb = hull(R1, R2)
    e1 = acb(0, Tb).exp()
    eta = Rb * e1
    c = eta * eta / 4
    w = 1 - z * z
    if not (abs(w) > 0):
        return None
    gp = 2 * z * (c - 1 / w)
    sg = eta / 2 * (1 + z)
    r = abs(sg)
    if not (r > 0):
        return None
    hp = pos(r - 1 / r)
    a = sg.real
    ta = pos(abs(a) - 1)
    if ta <= 0:
        G = hp * sg / r
    else:
        if a > 0:
            sgn = 1
        elif a < 0:
            sgn = -1
        else:
            return None
        G = hp * sg / r - 2 * sgn * ta
    Gc = G.conjugate()
    dS = (gp * H).real + (Gc * eta / 2 * H).real
    dR = (z * z * Rb * e1 * e1 / 2).real + (Gc * e1 * (1 + z) / 2).real
    dT = (2 * acb(0, 1) * c * z * z).real + (Gc * acb(0, 1) * sg).real
    return E0 + abs(dS) + abs(dR) * rR + abs(dT) * rT


def y_enclosed(sb, lo, hi, d0):
    """Lemma A.10: Y(η) = min(κ(η)δ₀, θ) ∈ [lo, hi] for every η of the sub-box sb = (r1, r2, t1, t2) inside the sector."""
    r1, r2, t1, t2 = sb
    Rb, Tb = hull(r1, r2), ball_top(t1, t2)
    om2 = (Rb / 2) ** 2
    sn = Tb.sin()
    den = pos(1 - om2) + 2 * om2 * sn * sn
    tlo, thi = arb(t1), ((arb.pi() / 2).upper() if t2 == PI2 else arb(t2))
    if not (den > 0):
        ylo, yhi = arb(0), thi                                          # 0 ≤ Y ≤ θ
    else:
        K = om2 * pos((2 * Tb).sin()) / den * d0
        klo, khi = K.lower(), K.upper()
        klo = klo if klo > 0 else arb(0)
        ylo = klo if klo < tlo else tlo
        yhi = khi if khi < thi else thi
    return bool(ylo >= arb(lo) and yhi <= arb(hi))


def chain(pieces, A, B):
    return bool(pieces) and pieces[0][0] == A and pieces[-1][1] == B and all(p[1] == q[0] for p, q in zip(pieces, pieces[1:]))


def verify_outer(meta, outer, log=print):
    ctx.prec = VA_PREC
    m = meta['VA011']
    e0, d0, er, R1 = (FX(m[k]) for k in ('eps0', 'delta0', 'end_r', 'R1'))
    ok_all = True
    t0 = time.time()
    nseg, nf, worst, fails = 0, 0, -math.inf, []
    for k, o in enumerate(outer):
        for P, Q in o['conn'] + o['poly']:
            nseg += 1
            b = seg_bound(P, Q, o['box'])
            if b is None or not (b < 0):
                nf += 1
                fails.append(k)
                continue
            worst = max(worst, float(b.upper()))
    good = nf == 0 and nseg == m['outer']['segments'] and len(outer) == m['outer']['boxes'] and nseg > 0
    log(f"[{'ok' if good else 'FAIL'}] Lemma A.9 (outer segments): {nseg} segments on {len(outer)} η-boxes, {nf} failed"
        f"{' (boxes ' + str(sorted(set(fails))[:5]) + ')' if nf else ''}, largest bound {worst:.4g} [{time.time() - t0:.1f}s]")
    ok_all &= good
    t0 = time.time()
    bad, nsub = [], 0
    for k, o in enumerate(outer):
        R1_, R2_, f1, f2 = o['box']
        lo, hi, yc, pts = o['lo'], o['hi'], o['yc'], o['pts']
        if not (math.isfinite(lo) and math.isfinite(hi) and lo <= hi):
            bad.append((k, 'connector [lo, hi]'))
            continue
        if not (all(P.real == d0 and Q.real == d0 for P, Q in o['conn']) and chain(o['conn'], complex(d0, lo), complex(d0, hi))):
            bad.append((k, 'connector chain'))
        if not (lo <= yc <= hi and pts and pts[0] == complex(d0, yc)):
            bad.append((k, 'y_c / polyline start'))
        verts = {x for s in o['poly'] for x in s}
        if not (pts and chain(o['poly'], pts[0], pts[-1]) and all(p in verts for p in pts)):
            bad.append((k, 'polyline chain'))
        if not (pts and abs(1 - acb(pts[-1].real, pts[-1].imag)) < er):
            bad.append((k, 'polyline end'))
        if not (0 < R1_ <= R2_ <= 2.0 and 0.0 <= f1 <= f2 <= PI2):
            bad.append((k, 'outside the sector'))
            continue
        ok_t, why = tiles(o['box'], o['sub'])
        if not ok_t:
            bad.append((k, 'sub-boxes: ' + why))
        nsub += len(o['sub'])
        if not all(y_enclosed(s, lo, hi, d0) for s in o['sub']):
            bad.append((k, 'Y not enclosed'))
    good = not bad and nsub == m['handoff_sub_boxes'] and nsub > 0
    log(f"[{'ok' if good else 'FAIL'}] Lemma A.10 (connector and handoff): {len(outer)} η-boxes, connectors and polylines chained, "
        f"{nsub} sub-boxes tiling the η-boxes with Y(η) ∈ [lo, hi] proved, {len(bad)} failures{' ' + str(bad[:5]) if bad else ''}"
        f" [{time.time() - t0:.1f}s]")
    ok_all &= good
    t0 = time.time()
    good, why = tiles((R1, 2.0, 0.0, PI2), [o['box'] for o in outer])
    log(f"[{'ok' if good else 'FAIL'}] Lemma A.11 (tiling): {len(outer)} η-boxes, [{R1}, 2] × [0, π/2]: {why} (exact area, disjoint "
        f"interiors) [{time.time() - t0:.1f}s]")
    ok_all &= good
    e, R = arb(er), arb(FX(m['end']['R2']))
    r = R / 2 * (2 + e)
    b = e.log() + (2 + e).log() + R * R / 4 * (1 + e) ** 2 + (r * r - 1) / 2 - r.log()
    good = bool(r > 1 and b < 0 and FX(m['end']['R2']) == 2.0 and FX(m['end']['eps_e']) == er)
    log(f"[{'ok' if good else 'FAIL'}] Lemma A.12 (end piece): R2 = 2, ε_e = {er}, bound {float(b.upper()):.5g}")
    return ok_all and good


def verify_meta(meta, log=print):
    """The parameters are the certificates'; Lemma A.8 (i)–(iii)'s constants; the seam between D and VA011, in rationals."""
    try:
        m, dm = meta['VA011'], meta['D']
        e0, d0, er, QA, PB, UJ, R1 = (FX(m[k]) for k in ('eps0', 'delta0', 'end_r', 'QA', 'PB', 'UJ', 'R1'))
        DR = FX(dm['R2'])
    except (KeyError, ValueError, TypeError) as x:
        log(f'[FAIL] meta: unreadable ({x})')
        return False
    why = []
    if meta.get('format') != 1:
        why.append('format')
    if (e0, d0, er, QA, PB, UJ, DR) != tuple(THEOREM[k] for k in ('eps0', 'delta0', 'end_r', 'QA', 'PB', 'UJ', 'D_R2')):
        why.append("parameters are not the certificates'")
    if R1 != 2.0 - e0:
        why.append('R1 ≠ fl(2 − ε₀)')
    if not (Fr(e0) / 4 < Fr(d0) and Fr(e0) < Fr(1, 5)):
        why.append('Lemma A.8 (i): ε₀/4 < δ₀')
    if not (Fr(QA) * Fr(PB) > 1 and Fr(QA) >= 1):
        why.append('Lemma A.8 (ii)–(iii): q_A·p_B > 1, q_A ≥ 1')
    if not (Fr(2) - Fr(e0) <= Fr(DR) and Fr(R1) <= Fr(DR)):
        why.append('seam')
    if not (Fr(2) - Fr(e0) <= Fr('1.89') and Fr(R1) <= Fr('1.89')):     # the real annulus 1.89 ≤ |η| < 2 (Proposition A.13)
        why.append('1.89')
    if not (m.get('prec') == VA_PREC and dm.get('prec') == D_PREC):
        why.append('precision')
    log(f"[{'ok' if not why else 'FAIL'}] meta and seam: ε₀ = {e0}, δ₀ = {d0}, END_R = {er}, q_A = {QA}, p_B = {PB}, u_J = {UJ}; "
        f"R1 = fl(2 − ε₀) = {R1}; ε₀/4 < δ₀; q_A·p_B > 1; seam: 2 − ε₀ and R1 at most 1.89 (real) ≤ {DR} (rationals)"
        f"{': ' + ', '.join(why) if why else ''}")
    return not why


def E_def(z, eta):
    """E(z; η) of (A.1) from its definition: log|1 − z²| + Re(η²z²/4) + Λ(ν), ν = η(1 + z)/2, Λ(ν) = H(|ν|) − (|Re ν| − 1)₊² for
    |ν| > 1 and 0 for |ν| ≤ 1, H(ϱ) = (ϱ² − 1)/2 − log ϱ (balls; None where |ν|'s ball meets 1)."""
    nu = eta * (1 + z) / 2
    r = abs(nu)
    if r > 1:
        Lm = (r * r - 1) / 2 - r.log() - pos(abs(nu.real) - 1) ** 2
    elif r <= 1:
        Lm = arb(0)
    else:
        return None
    return abs(1 - z * z).log() + (eta * eta * z * z / 4).real + Lm


def forms_sample(meta, n=1000, seed=136, log=print):
    """Informational, not a check (referee round 136, M2 (iii)): Lemma A.4's forms as implemented here (fA, fB) against the
    definition (A.1) of E at random points of the path P_η, 0 < ε ≤ ε₀, 0 ≤ θ ≤ π/2, n points in each of four groups sampled in
    the forms' own coordinates: f_A on the floor F or the cap C (u ∈ (0, 4], q ∈ [0, q_A], X ≤ δ₀, uX ≤ ε₀), f_A on the junction J
    (X = x_b = ε/4, q ∈ [0, q_A]), and f_B on F and on J (p ∈ (0, p_B], X = pϑ ≤ δ₀, uX ≤ ε₀) with its dropped term
    −((|Re ν| − 1)₊/X)² added back. Prints the largest relative deviation |f·X² − E|/|E|. Never raises and never changes the exit
    status. The sample points are doubles, taken exactly into balls."""
    try:
        ctx.prec = VA_PREC
        m = meta['VA011']
        e0, d0, QA, PB = (FX(m[k]) for k in ('eps0', 'delta0', 'QA', 'PB'))
        rnd = random.Random(seed)
        cnt = dict(AF=0, AC=0, AJ=0, BF=0, BJ=0)
        worst, skipped, meets = 0.0, 0, [0]

        def setup(eps, th):
            e, t = arb(eps), arb(th)
            om = 1 - e / 2
            eta = (2 - e) * acb(0, t).exp()
            A = om * om * acb(0, 2 * t).exp() - 1
            return e, t, eta, A.imag / abs(A.real)

        def compare(f, X, z, eta, dropped=False):
            nonlocal worst
            Ed = E_def(z, eta)
            if Ed is None or not Ed.is_finite() or Ed == 0:
                return False
            g = f * X * X
            if dropped:                                                # f_B drops −(|Re ν| − 1)₊²; add it back
                nu = eta * (1 + z) / 2
                g = g - pos(abs(nu.real) - 1) ** 2
            if not g.is_finite():
                return False
            worst = max(worst, float(abs(g.mid() - Ed.mid())) / abs(float(Ed.mid())))
            meets[0] += bool(g.overlaps(Ed))
            return True

        def logu(a, b):
            return a * (b / a) ** rnd.random()

        for group in ('A', 'AJ', 'B', 'BJ'):
            for _ in range(n):
                try:
                    if group == 'A':                                   # F or C: x_b ≤ X ≤ δ₀, θ² = qX
                        u, q = rnd.uniform(1e-3, 4.0), rnd.uniform(0.0, QA)
                        Xf = logu(1e-8, min(d0, e0 / u))
                        e, t, eta, k = setup(u * Xf, math.sqrt(q * Xf))
                        X = arb(Xf)
                        if not (4 * X >= e):
                            skipped += 1
                            continue
                        br = 'F' if k * X < t else ('C' if k * X > t else None)
                        if br is None:
                            skipped += 1
                            continue
                        z = acb(X, k * X if br == 'F' else t)
                        if compare(fA(br, e / X, t * t / X, X), X, z, eta):
                            cnt['A' + br] += 1
                    elif group == 'AJ':                                # J: X = x_b = ε/4, 0 ≤ Y ≤ min(κx_b, θ), λ = Y/(κx_b)
                        epsf, q = logu(1e-8, e0), rnd.uniform(1e-6, QA)
                        e, t, eta, k = setup(epsf, math.sqrt(q * epsf / 4))
                        X = arb(epsf / 4)
                        Y = (k * X if k * X < t else t) * rnd.random()
                        if compare(fA('J', e / X, t * t / X, X, Y / (k * X)), X, acb(X, Y), eta):
                            cnt['AJ'] += 1
                    elif group == 'B':                                 # F in (u, p, ϑ): X = pϑ ≤ δ₀, uX ≤ ε₀
                        u, pp = rnd.uniform(1e-3, 4.0), rnd.uniform(1e-4, PB)
                        th2 = logu(1e-10, min(PI2 * PI2, d0 / pp, e0 / (u * pp)))
                        thf = math.sqrt(th2)
                        e, t, eta, k = setup(u * pp * thf * thf, thf)
                        X = arb(pp * thf * thf)
                        if not (4 * X >= e and k * X < t):             # the floor (p < p_B < 1 forces it, Lemma A.8 (iii))
                            skipped += 1
                            continue
                        if compare(fB('F', e / X, X / (t * t), t * t), X, acb(X, k * X), eta, dropped=True):
                            cnt['BF'] += 1
                    else:                                              # J in (u, p, ϑ): X = x_b = ε/4, θ² = X/p ≤ (π/2)²
                        pp = rnd.uniform(1e-4, PB)
                        epsf = logu(1e-8, min(e0, 4 * pp * PI2 * PI2))
                        e, t, eta, k = setup(epsf, math.sqrt(epsf / 4 / pp))
                        X = arb(epsf / 4)
                        Y = (k * X if k * X < t else t) * rnd.random()
                        if compare(fB('J', e / X, X / (t * t), t * t, Y / (k * X)), X, acb(X, Y), eta, dropped=True):
                            cnt['BJ'] += 1
                except Fail:
                    skipped += 1
        tot = sum(cnt.values())
        log(f"[info] forms against (A.1) (informational, not a check): {tot} random path points, 0 < ε ≤ ε₀ — f_A on F {cnt['AF']}, "
            f"C {cnt['AC']}, J {cnt['AJ']}; f_B on F {cnt['BF']}, J {cnt['BJ']} (its dropped term added back); largest relative "
            f"deviation {worst:.1e}, the two balls overlapping at {meets[0]} of {tot}{f' ({skipped} draws skipped)' if skipped else ''}")
    except Exception as x:                                             # informational: never a failure of the run
        log(f'[info] forms against (A.1) (informational, not a check): not evaluated ({x!r})')


def run(d, only=('D', 'VA'), d_sample=1, log=print):
    """Every check on the data in d. Returns the exit status (0 verified, 1 a failure, 2 unreadable data)."""
    t0 = time.time()
    ok, msg = check_sums(d)
    log(f"[{'ok' if ok else 'FAIL'}] SHA256SUMS: {msg}")
    if not ok:
        log('NOT VERIFIED: the data are not the hashed ones')
        return 1
    try:
        meta, leaves, inner, outer = load(d)
    except (OSError, ValueError, KeyError, TypeError, IndexError) as x:
        log(f'[FAIL] data unreadable: {x!r}')
        log('NOT VERIFIED')
        return 2
    ok = verify_meta(meta, log)
    if 'D' in only:
        ok &= verify_D(meta, leaves, d_sample, log)
    if 'VA' in only:
        ok &= verify_inner(meta, inner, log)
        ok &= verify_outer(meta, outer, log)
        forms_sample(meta, log=log)
    gaps = ([] if set(only) == {'D', 'VA'} else ['only ' + ','.join(only)]) + ([f'Lemma A.2 on every {d_sample}th leaf'] if d_sample != 1 else [])
    if ok:
        claim = ([f"D (Proposition A.3): the domination condition on 1 < |η| ≤ {THEOREM['D_R2']}"] if 'D' in only else []) + (
            [f"VA (Proposition A.13): the valley condition along P_η on 2 − ε₀ ≤ |η| < 2, ε₀ = {THEOREM['eps0']}"] if 'VA' in only else [])
        log(f"{'PARTIALLY VERIFIED (' + '; '.join(gaps) + ')' if gaps else 'VERIFIED'}: {'; '.join(claim)}; 0 ≤ arg η ≤ π/2 "
            f"[{time.time() - t0:.0f}s]")
        return 0
    log(f'NOT VERIFIED [{time.time() - t0:.0f}s]')
    return 1


def main(argv=None):
    ap = argparse.ArgumentParser(description="Verify the stored certificates D and VA at ε₀ = 0.11 (paper 18's Appendix A; a cross-check since round 148).")
    ap.add_argument('--dir', default=DEFAULT_DIR, help='the certificates directory (default results/certificates)')
    ap.add_argument('--only', default='D,VA', help='D, VA or D,VA (default)')
    ap.add_argument('--d-sample', type=int, default=1, help="check Lemma A.2's bound on every N-th leaf only (the tiling always on all)")
    ap.add_argument('--quiet', action='store_true')
    a = ap.parse_args(argv)
    only = tuple(x for x in a.only.split(',') if x)
    if any(x not in ('D', 'VA') for x in only) or a.d_sample < 1:
        ap.error('--only takes D, VA; --d-sample ≥ 1')
    return run(a.dir, only, a.d_sample, (lambda s: None) if a.quiet else (lambda s: print(s, flush=True)))


if __name__ == '__main__':
    sys.exit(main())
