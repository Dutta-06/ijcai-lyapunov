# The stored certificates D and VA (paper 18, Appendix A: a cross-check since round 148)

Referee rounds 135 (I4, review M5) and 136 (J3–J6). This directory holds the data of the two computer-assisted certificates
that the proof of Theorem 2 of paper 18 (the transform's limit) used until referee round 148:

- **D**: the domination condition on 1 < |η| ≤ 1.8927 (Lemma A.2, Proposition A.3).
- **VA011**: the valley condition on 1.89 ≤ |η| < 2, that is VA at ε₀ = 0.11 (Lemmas A.4–A.12, Proposition A.13).

Since round 148 the theorem does not depend on them: the valley condition on 1 < |η| < 2 is proved by hand, along a
three-segment path (paper 18's Proposition 4.8; `review/RESEARCH_note_annulus.md`, Theorem A, refereed in round 148,
`review/REVIEW_round148_D3eq.md`). D and VA011 are now an independent computer-assisted cross-check of that hand proof,
and Appendix A, which states them, moved from the paper to its supplement (`manuscript/paper18/supplement.tex`) with
its numbering unchanged.

`tools/verify_certificates.py` checks these certificates from the files here alone. It does not re-run the
programs that proposed the boxes and paths, and it imports nothing from `experiments/`. Since round 136 it is the
check of record for these certificates: Appendix A and the paper's code-availability statement describe what this
verifier does.

V (1.8927 ≤ |η| ≤ 1.995) and VA at ε₀ = 0.005 are further cross-checks (V is §S5 of the supplement). They are not
exported here.

## Labels: Appendix A's numbers and the drafts' names

The verifier's output, its docstring and this file use Appendix A's numbers (the paper's until round 148, the
supplement's since; unchanged by the move). The drafts (`review/drafts/D61.md`,
Appendix A) and the review rounds up to 135 used other names:

| paper 18 | drafts (D61) | content |
|---|---|---|
| Theorem 2 | Theorem 3 | the transform's limit on \|η\| < 2 |
| Lemma A.1 | Lemma D1 | the supremum over x (analytic; not checked by the verifier) |
| Lemma A.2 | Lemma D2 | the box bound of D |
| Proposition A.3 | certificate D | D's leaves and their tiling |
| Lemma A.4 | Lemma VA1 | the forms f_A, f_B (algebra proved in Appendix A) |
| Lemma A.5 | Lemma VA2 | series and tails |
| Lemma A.6 | Lemma VA3 | closed forms at ϱ = 0 (algebra proved in Appendix A) |
| Lemma A.7 | Lemma VA4 | the mean value form |
| Lemma A.8 | Lemma VA5 | domains and cover |
| Lemma A.9 | Lemma VA6 | the outer piece |
| Lemma A.10 | Lemma VA7 | connector and handoff |
| Lemma A.11 | Lemma VA8 | tiling of the η-boxes |
| Lemma A.12 | Lemma VA9 (the end piece) | the end piece |
| Proposition A.13 | certificate VA | VA at ε₀ = 0.11 |
| (A.1) | (A.1) | the bound E of the valley condition |

## Files

| file | bytes | content |
|---|---|---|
| `meta.json` | 1930 | the parameters, the counts, the formats |
| `D_leaves.f64le.gz` | 4408325 | D's 402664 accepted leaves |
| `VA011_inner.jsonl.gz` | 66957 | the 4189 accepted parameter boxes of the inner pieces A-J, A-FC, B-J, B-F |
| `VA011_outer.jsonl.gz` | 428496 | the 220 accepted η-boxes of the outer piece, with their connectors, polylines, segments and handoff sub-boxes |
| `SHA256SUMS` | 334 | `sha256sum` format, over the four files above (`sha256sum -c SHA256SUMS`) |

All floats are exact IEEE-754 binary64 values. In the JSON files each float is a `float.hex` string, for example
`"0x1.c28f5c28f5c29p-4"` for ε₀ = 0.11; `float.fromhex` reads it back exactly. The `.gz` files are gzip with
mtime 0, so their bytes depend only on the data, and a re-export gives the same hashes.

- **`meta.json`**: `D` holds the rectangle [R1, R2] × [f1, f2] = [1, 1.8927] × [0, fl(π/2)], the leaf count, the
  precision (120 bits) and the layout. `VA011` holds ε₀, δ₀ = 0.05, END_R = 0.02, q_A = 2.5, p_B = 0.5, the
  junction's u = 4, R1 = fl(2 − ε₀) = 1.89, the precision (128 bits), the series' terms (S and C_v 16, Λ 30, H₁ 26),
  the tail floor 10⁻³⁰⁰, the per-piece counts (boxes accepted and skipped), the outer counts, the number of handoff
  sub-boxes, and the end piece's R2 = 2 and ε_e = END_R. The `proposer_*` entries are the proposing program's own
  bounds. They are informational only, and the verifier recomputes every bound itself.
- **`D_leaves.f64le.gz`**: gzip of n × 5 little-endian binary64 values, row-major, one row per leaf in the order
  accepted: R1, R2, f1, f2 (the leaf is R1 ≤ |η| ≤ R2, f1 ≤ arg η ≤ f2), then H, the proposer's rigorous bound
  (informational).
- **`VA011_inner.jsonl.gz`**: one JSON object per line, one line per accepted box. Each line has:
  - `piece`: one of A-J, A-FC, B-J, B-F.
  - `iv`: the box's intervals, as they were evaluated (after the domain clip). Form A uses `u`, `q`, `r` (= ϱ);
    form B uses `u`, `p`, `t` (= ϑ); the junction adds `l` (= λ), with `u` = [4, 4].
  - `recs`: `[branch, route, bound]` for each branch evaluated, where route `d` means f₀ ≤ 0 and sup f′ < 0, and
    route `v` means the mean value bound. This is informational.
- **`VA011_outer.jsonl.gz`**: one JSON object per line, one line per accepted η-box, with these keys:
  - `box`: [R1, R2, f1, f2].
  - `lo`, `hi`: the connector's ends at Re z = δ₀.
  - `yc`: the polyline's start δ₀ + i·y_c.
  - `pts`: the polyline's vertices [x, y].
  - `conn`: the connector's certified pieces [Px, Py, Qx, Qy], in the order certified, which is the bisected
    segments chaining δ₀ + i·lo to δ₀ + i·hi.
  - `poly`: the polyline's certified pieces, chaining the first vertex to the last.
  - `sub`: the handoff's sub-boxes [r1, r2, t1, t2], on which Y(η) ∈ [lo, hi] was proved.

## Reproduce and verify

```
# export (recomputes D and VA011, asserts equality with results/chaos/gh_peano_raw.npz, writes this directory; ~21 s)
OMP_NUM_THREADS=1 venv/bin/python experiments/exp_gh_peano.py --export-cert [DIR]
# verify (about 14 s on the shared four-core machine, one process)
OMP_NUM_THREADS=1 venv/bin/python -I tools/verify_certificates.py --dir results/certificates [--only D,VA] [--d-sample N]
```

The export flag writes nothing else. The md and the raw are untouched, and `--from-raw` still reproduces
`results/chaos/gh_peano.md` byte for byte. The only changes to the proposing programs are three optional `keep=`
hooks, which collect what was accepted:

- `dom_bb(keep=)` in `experiments/exp_gh_peano.py`;
- `bb(keep=)`, which already existed, in `experiments/gh_lastann.py`;
- `_y_enclosed(keep=)` in `experiments/gh_lastann.py`.

The verifier's output on the stored data (exit status 0; round 136):

```
[ok] SHA256SUMS: 4 files match
[ok] meta and seam: ε₀ = 0.11, δ₀ = 0.05, END_R = 0.02, q_A = 2.5, p_B = 0.5, u_J = 4.0; R1 = fl(2 − ε₀) = 1.89; ε₀/4 < δ₀; q_A·p_B > 1; seam: 2 − ε₀ and R1 at most 1.89 (real) ≤ 1.8927 (rationals)
[ok] Lemma A.2 (box bound < 0): 402664 leaves, 0 failed, largest bound -1.9432e-11 (max |Δ| against the proposer's bound 6e-17) [7s]
[ok] Proposition A.3 (D's tiling): 402664 leaves, [1, 1.8927] × [0, π/2]: tiled (exact area, disjoint interiors) [3s]
[ok] Lemmas A.4–A.7 A-J: 54 boxes, 0 failed; route d (f₀ ≤ 0, sup f′ < 0) 27, largest sup f′ -0.013; route v 27, largest bound -0.06487 (max |Δ| against the proposer's 0) [0.0s]
[ok] Lemma A.8 A-J (branches above; cover): the domain covered by its 54 boxes (72 pieces inside boxes, 0 outside the domain) [0.0s]
[ok] Lemmas A.4–A.7 A-FC: 2533 boxes, 0 failed; route d (f₀ ≤ 0, sup f′ < 0) 101, largest sup f′ -0.00458; route v 3184, largest bound -8.781e-14 (max |Δ| against the proposer's 0) [0.8s]
[ok] Lemma A.8 A-FC (branches above; cover): the domain covered by its 2533 boxes (3370 pieces inside boxes, 169 outside the domain) [0.4s]
[ok] Lemmas A.4–A.7 B-J: 413 boxes, 0 failed; route d (f₀ ≤ 0, sup f′ < 0) 13, largest sup f′ -0.0837; route v 400, largest bound -0.0004668 (max |Δ| against the proposer's 0) [0.1s]
[ok] Lemma A.8 B-J (branches above; cover): the domain covered by its 413 boxes (530 pieces inside boxes, 58 outside the domain) [0.1s]
[ok] Lemmas A.4–A.7 B-F: 1189 boxes, 0 failed; route d (f₀ ≤ 0, sup f′ < 0) 0, largest sup f′ -inf; route v 1189, largest bound -5.329e-05 (max |Δ| against the proposer's 0) [0.3s]
[ok] Lemma A.8 B-F (branches above; cover): the domain covered by its 1189 boxes (1195 pieces inside boxes, 106 outside the domain) [0.2s]
[ok] Lemma A.9 (outer segments): 16083 segments on 220 η-boxes, 0 failed, largest bound -4.228e-09 [0.3s]
[ok] Lemma A.10 (connector and handoff): 220 η-boxes, connectors and polylines chained, 1306 sub-boxes tiling the η-boxes with Y(η) ∈ [lo, hi] proved, 0 failures [0.0s]
[ok] Lemma A.11 (tiling): 220 η-boxes, [1.89, 2] × [0, π/2]: tiled (exact area, disjoint interiors) [0.0s]
[ok] Lemma A.12 (end piece): R2 = 2, ε_e = 0.02, bound -1.3314
[info] forms against (A.1) (informational, not a check): 4000 random path points, 0 < ε ≤ ε₀ — f_A on F 892, C 108, J 1000; f_B on F 1000, J 1000 (its dropped term added back); largest relative deviation 2.1e-16, the two balls overlapping at 4000 of 4000
VERIFIED: D (Proposition A.3): the domination condition on 1 < |η| ≤ 1.8927; VA (Proposition A.13): the valley condition along P_η on 2 − ε₀ ≤ |η| < 2, ε₀ = 0.11; 0 ≤ arg η ≤ π/2 [13s]
```

A-FC's route v count, 3184, exceeds its 2533 boxes because a box on which neither α > b̂ nor α < b̂ is proved is
evaluated on both branches, F and C.

## What the verifier re-implements, and what it does not check

The verifier re-implements every computational step of D and VA011 from the lemmas' statements. It shares no code
with the proposers. Its forms f_A, f_B and closed forms are re-typed from the statements of Lemmas A.4 and A.6 (with the
proposer's variable names and order of operations), so it shares their algebra with the proposer; that algebra is proved
in Appendix A, not by the verifier, and the verifier's informational line samples the forms against the definition (A.1).

- **Lemma A.2**: the box bound, at 120 bits, on every leaf. The candidates are w_m, the kink, and the two quadratics'
  roots by the stable formula; the top edge is the true π/2.
- **Proposition A.3 (D's tiling)**: the leaves lie inside [1, 1.8927] × [0, fl(π/2)] with positive area. Their areas
  sum to the rectangle's exactly, using dyadic rationals as integers. Their interiors are pairwise disjoint, checked by
  a sweep over R that keeps the active φ-intervals sorted.
- **Lemmas A.4–A.7**: on every stored inner box and branch, at 128 bits:
  - the forms f_A and f_B (Lemma A.4);
  - the series S, C_v, Λ and H₁ with their tails (Lemma A.5). The tails are computed in Arb, so they need no 1 + 10⁻⁷
    factor (that factor belongs to the proposer's double-precision tails), and are floored at 10⁻³⁰⁰;
  - the closed forms at ϱ = 0 (Lemma A.6), in outward-rounded endpoint intervals;
  - the derivative, by dual numbers;
  - the mean value form's three cases (Lemma A.7).
- **Lemma A.8**: the branches are decided by the verifier itself: F, C or both for form A, and the floor proved for
  B-F. The junction's u is exactly 4. The cover is checked for each piece: its stored boxes cover its whole closed
  parameter domain, for example 0 ≤ u ≤ 4, 0 ≤ q ≤ q_A, 0 ≤ ϱ ≤ δ₀, uϱ ≤ ε₀ for A-FC. This is proved by decomposing the
  domain's bounding box at the boxes' faces, in rationals, with (π/2)² handled in Arb. The proposer's root grids, clips,
  float-midpoint bisection and skipped boxes are not used; the paper's Lemma A.8 is stated in this form since round 136.
- **Lemma A.9**: the mean value bound, on every stored segment over its η-box.
- **Lemma A.10**: the connector [lo, hi] is finite. Its pieces lie on Re z = δ₀ and chain lo to hi. y_c ∈ [lo, hi], and
  the polyline starts at δ₀ + i·y_c. The polyline's pieces chain its first vertex to its last, through every stored
  vertex. |1 − z_n| < END_R holds in balls. The sub-boxes tile the η-box, and Y(η) ∈ [lo, hi] is proved on each.
- **Lemma A.11**: the η-boxes tile [R1, 2] × [0, fl(π/2)].
- **Lemma A.12**: the end piece's bound.
- **meta**: the parameters are the certificates'. ε₀/4 < δ₀, q_A·p_B > 1 and q_A ≥ 1 (the constants behind Lemma A.8
  (i)–(iii)). The seam is checked in rationals: 2 − ε₀ and R1 are at most the real 1.89, which is at most 1.8927.
- **Informational (not a check, never changes the exit status)**: the forms f_A (on F, C and J) and f_B (on F and J, its
  dropped term added back) against the definition (A.1) of E at 4000 random path points; the largest relative deviation
  and the number of points where the two balls overlap are printed.

The verifier does not check the analytic statements, which are proved in Appendix A or the paper:

- Lemma A.1;
- the algebra of Lemmas A.4 and A.6 (that the forms and closed forms are what they claim);
- Lemma A.8 (i)–(iii) beyond their constants;
- the real segment z ∈ [−1, x_b] of the path;
- the reflections η ↦ η̄, −η;
- the deduction of the transform's limit from the two conditions (Lemma 4.6 of the paper; since round 148 Theorem 2
  itself is proved without D or VA).

Nothing computational of D or VA011 remains outside the verifier.

**Trusted base.** The verifier relies on the following:

- Arb, through python-flint 0.9.0;
- python-flint's exact conversion of a double to a ball;
- `float.fromhex` and the binary64 reading of `array`;
- CPython's exact float comparisons and the operations listed in the verifier's docstring;
- Python integers and `fractions.Fraction`;
- hashlib, gzip and json;
- the interpreter.

## Tests

The tests are in `tests/test_gh_peano.py`, under `round135`.

- **The verifier passes on the stored data.** This test runs the full verifier and also runs it in a fresh
  interpreter, checking that nothing from `experiments/` is imported.
- **The verifier fails on planted corruptions**, and each corruption is caught by its own check:
  - a D leaf removed (area);
  - a D leaf replaced by an equal-sized copy (overlap);
  - D's worst leaf widened (its bound is no longer negative);
  - A-FC's worst box widened in q (bound +3.0·10⁻¹⁴);
  - an inner box removed (cover);
  - a B-F box widened in p (floor not provable);
  - a polyline vertex moved, a segment's end moved, or the connector's end moved (Lemma A.10);
  - a sub-box removed (Lemma A.10);
  - an η-box removed (Lemma A.11);
  - a polyline piece moved far (Lemma A.9);
  - a hash mismatch;
  - a leaf removed with SHA256SUMS rewritten (exit status 1).
- **The re-implementation against the proposers** on random boxes: the same bounds to rounding, and the
  verifier's radii no wider. The forms are also checked against E's definition at path points.
- **The export reproduces this directory byte for byte.**

Round 136 changed only the verifier's labels (to the paper's numbers) and added its informational line. The tests'
pinned output strings were updated accordingly. No check was changed or weakened. Round 148 changed only wording (the
docstring, the argparse description, one meta message: "parameters are not the certificates'"), so that nothing says
the theorem depends on the certificates; no check or printed label changed.
