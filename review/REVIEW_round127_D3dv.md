# Referee report D3dv, round 127: Proposition D‴.61 as revised after round 126

## Scope

**What I reviewed.** Commit 31487a6 on claude/vigilant-clarke-hf939b. I exported it with `git archive 31487a6 | tar -x` into W/repo, where W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref127. W/repo/venv is a symlink to the repository's venv (Python 3.13.18, python-flint 0.9.0). The repository's HEAD moved to d3db5ac during the round because of other sessions. I reviewed only 31487a6.

**Hashes (sha256 prefixes) [own]:**
- gh_peano.md 2258956259220fdd
- gh_peano_raw.npz 44a59d992d31adf2
- gh_lastann.py 7289dd4a036517aa
- exp_gh_peano.py cff57b8def76e017
- test_gh_peano.py 915cad8e4d0d96cd
- D61.md eb451c09d29d735f
- round126_D61_own.py 5893a7f453bbc51c

**How I ran things.**
- I used read-only git only.
- Everything ran under `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 10` with GP_OUT inside W, and never more than two heavy processes at once.
- I killed no process; every process of mine exited on its own.

**Read in full:**
- D61.md;
- in exp_gh_peano.py: parts D, DB, V, VC and VA, `d_lines`, `v_lines`, `vc_lines`, `va_lines`, `report`, and everything they call (`lam_plus`, `gap_closed`, `bound_sup`, `box_bound(_ball)`, `dom_bb`, `cert_*`, `c_need(_ball)`, `r_star`);
- gh_lastann.py;
- gh_peano.md;
- test_gh_peano.py from the part V tests to the end, the rest skimmed;
- REVIEW_round126_D3du.md;
- the REVIEW_LOG entries for rounds 120–126;
- `git diff 5f56b63 31487a6` for these files;
- for the status of D3du-6, the relevant parts of HANDOFF.md.

## Summary

**Design.** Round 126 made these changes:
- `handoff_check`:
  - the η-boxes' areas are summed in rationals and must equal (2 − R1)·fl(π/2) exactly, and `tiled` = exact ∧ inside ∧ disjoint;
  - the polyline's last vertex must lie inside |1 − z| < END_R in balls (the end link);
- `_y_enclosed` returns None outside 0 < R1 ≤ R2 ≤ 2, 0 ≤ f1 ≤ f2 ≤ π/2, and after 200 000 sub-boxes;
- `soundness` takes corners over the non-degenerate coordinates at both ends of the expansion variable, each distinct point once;
- new plants: T1, a box outside the sector, P8, the connector's last piece and the whole connector, a middle polyline half-piece, `_y_enclosed` → None, lo raised, the true π/2 on the top strip, `joints=False`;
- `tools/mutants/round126_D61_own.py` (R00–R10);
- D61 and md wording.

This round's main question is the proof as a whole. Conjecture P is now stated as proved, computer-assisted, on |η| < 2. It rests on (a), Duhamel's identity, the domination route (Ω) with D's branch and bound, and the valley route (Ω′) with V and VA.

**Verdicts.**
1. **Reproducibility is exact** [own]:
   - `--from-raw` regenerates the md byte for byte;
   - 42/42 tests pass (29.3 s);
   - part VA recomputes equal to the raw (4.8 s), and so does part DB: 4935 and 402 664 boxes, ball_worst −3.2813e-7 and −1.9432e-11 (15.5 s);
   - V's three bands with the ball recheck equal the raw dict for dict: 760 / 2011 / 10 658 boxes, ball_worst −2.4767e-12, −1.0258e-12, −1.0004e-12 (6.1 s, 29.3 s, 266.7 s).
2. **Round 126's changes are sound.**
   - The exact area plus containment plus disjoint interiors is a tiling of the float rectangle [1.995, 2] × [0, fl(π/2)]. The sliver up to the true π/2 is carried by the top boxes, which use the true π/2 everywhere.
   - The end link closes the chain to part V's end piece. It uses the same END_R constant as `piece_end`. On the data, |1 − end| ≤ 0.00995.
   - The budget and the preconditions only ever turn a proof into a failure, never the reverse.
   - The new plants test what they say. Two of them leave unsound survivors, though, because they sit far from the threshold or move two quantities at once (D3dv-3).
3. **The proof as a whole: I found no hole between the three certificates.**
   - The radial seams are the same floats: DOM_BB[−1] = CERT_BANDS[0][0] = 1.8927 and CERT_BANDS[−1][1] = 2 − EPS0 = 1.995.
   - The reflections hold: E is even under η ↦ η̄ (with z ↦ z̄) and under η ↦ −η (with σ ↦ −σ), and so is D's gap.
   - At arg η = π/2, D, V and VA all take the true π/2. At arg η = 0, VA's q = 0 face carries the certificate by its derivative, and the real axis is in Ω analytically anyway.
   - Uniformity on compact sets: ε ≥ ε_min, there are finitely many boxes, and on route-'d' boxes E ≤ X³·sup f′ ≤ −c(ε/4)³.
   - I re-derived the analytic steps: (a)'s representation, circle identity and prefactor; the contour bound → Λ⁺; the band near |σ| = 1; the chord and the 0-side; both axes; the Laplace constant; (iv)'s tails c² > 4 ln u/(u − 1) and mass c² > 4r. All are correct.
   - **The one load-bearing step neither proved nor verifiably cited is Duhamel's identity itself** (D3dv-1). It is true: I checked it exactly in rationals as a Beta-mixture identity. Its proof is five lines and should be written out.
   - The computer-assisted steps are exactly the three stated, plus the c_need thresholds in balls, which the md lists.
4. **The report's conclusion is not tied to the seams it relies on** (D3dv-2). `va_lines` prints "on the whole open disc |η| < 2" and "1.995 ≤ |η| < 2" without consulting V, D, VA's own R1, or the end piece. Only the tests' pins guard the seams.
5. **Wording.** The md and D61 match the code apart from cosmetic overstatements (D3dv-4, -5). The proof route and its computer-assisted parts are stated honestly.
6. **Mutants.**
   - Shipped R00–R10: 8 killed. R06 and R08 survive; I agree they are equivalent (reasons below).
   - Of my 15 single-point mutants of gh_lastann.py, 5 are killed. Of the 10 survivors, 5 are unsound in principle (Z02–Z06), 4 are sound or equivalent (Z07, Z08, Z09, Z13), and Z10 is a known sampling-tolerance survivor (S05).
   - Of the 3 on `va_lines`' guards: G01 is unsound in principle and survives; G02 and G03 are equivalent (redundant).

**No major issue.** The proof is complete modulo writing out Duhamel's identity.

**Findings:**
1. **D3dv-1 (minor).** Duhamel's identity is used without proof. It is cited [memory] to Saad 1992, a finite-matrix Krylov result, while J is infinite and unbounded, and the Status says it is "not claimed".
2. **D3dv-2 (minor).** VA's whole-disc conclusion is printed unconditionally on V, D, the 1.995 seam and the end piece.
3. **D3dv-3 (minor).** Round 126's plants leave unsound survivors:
   - Z02, the end link at 2·END_R;
   - Z03, `tiled` without `inside`;
   - Z04–Z06, the preconditions' other sides.
4. **D3dv-4 (cosmetic).** "Distinct" counts: A-FC's 24 770 are point–branch evaluations; it has 18 326 distinct points, so the total is 29 056, not HANDOFF's 35 500.
5. **D3dv-5 (cosmetic).** Wording in D61 and the md, six points (see the table).
6. **D3dv-6 (cosmetic).** HANDOFF residue.

## What I checked and how

**(1) Reproducibility** [own]
- **The md.** `GP_OUT=W/out exp_gh_peano.py --from-raw`, then `cmp` against the committed md: identical.
- **The tests.** `pytest -q -p no:cacheprovider --basetemp W/pt/base1 tests/test_gh_peano.py`: 42 passed in 29.34 s (W/logs/pytest_gh.log).
- **The parts against the raw.** W/scripts/recompute.py runs `part_VA()`, `part_DB()` and `cert_band(*CERT_BANDS[i], ball=True)` for i = 0, 1, 2. All are equal to the raw (W/logs/recompute.log). V's end pieces are −1.382677, −1.357139, −1.344304, and its ball_short counts are 220, 495, 2584, as printed.

**(2) Round 126's changes** [own, reading]
- **Exact tiling.** `Fraction` of each float is exact, and `inside` and `disjoint` are exact float comparisons. A finite union U of closed boxes inside a closed rectangle R, with disjoint interiors and Σ area = area(R), equals R: a point of R ∖ U has a neighbourhood of positive measure missing from U.
  - The rectangle's top is fl(π/2). The boxes touching it therefore have f2 == math.pi/2, and `th_ball`, `_y_enclosed` (Tb, t_hi) and `handoff_check`'s samples extend them to the true π/2. `inside` needs no tolerance.
  - My plant W/scripts/tile_plant.py shifts a bottom-row box down by 2⁻¹⁴ in |η|, which is exact in floats and preserves the area. It gives inside False, disjoint True, area_exact True, tiled False: correct, and stopped only by `inside`.
  - Round 126's own plant (a box to 2.001) also changes the area, so it cannot see whether `tiled` uses `inside` (Z03 below).
- **End link.**
  - `_chain(pl, pts[0], pts[−1])` makes the certified chain end at pts[−1], and the link tests |1 − pts[−1]| < END_R in balls.
  - `piece_end` bounds every segment from a point with |1 − z| ≤ END_R to z = 1 at R2 = 2. Its bound is monotone in R, which I checked term by term: log|1 − z²| ≤ log e + log(2 + e), Re(cz²) ≤ (R²/4)(1 + e)², C ≤ h((R/2)(2 + e)).
  - Measured: |1 − end| ∈ [0.00801, 0.00995] over the 76 boxes (W/logs/sound_count.log).
- **Budget and preconditions.** Both only return None, which `handoff_check` reports as 'Y not enclosed'. They cannot produce a false proof.
  - The clamps (pos(1 − ω²), pos(sin 2θ), k_lo ≥ 0) are valid exactly on the precondition's domain.
  - On the data, no box needs more than 50 sub-boxes: Z08 (budget 50) passes everything. The budget is therefore never exercised. The top-strip test is fast because of the depth limit, not the budget.
- **Soundness corners.** `free` and the `itertools.product` corners are as described. But in A-FC, 668 of the 2011 boxes carry both branches F and C, and each point is counted once per branch: 18 326 distinct points against 24 770 printed (D3dv-4). A-J, B-J and B-F have one branch per box, so their counts are distinct. The largest excess is −1.547e-13.
- **Shipped mutants.** Run with `tools/mutants.py --exp gh_lastann -k "part_va or md_is_the_report"` (W/logs/mut_shipped.log):
  - R01–R05, R07, R09 and R10 are killed.
  - R06 (θ's upper end at the float π/2) is equivalent: y_hi = min(k_hi, t_hi), and t_hi could bind only if sup Y = θ at the true π/2, which is impossible because κ = 0 there.
  - R08 (the fallback's lower end θ₁) is equivalent: den's ball meets 0 only if sin Tb ∋ 0, i.e. t1 = 0 under the precondition.
  - The md column is vacuous for `--exp gh_lastann`, as documented.

**(3) The proof as a whole** [own: re-derivation and exact checks]
- **Duhamel's identity**, (E − Q_K)[e^{ξz}] = (1/(K − 1)!)∫_0^ξ(ξ − u)^K e^{(ξ−u)²/2}Q_K[e^{uz}He_{K−1}]du.
  - **What the text does.** D61 introduces it as "the error representation of the K-step Lanczos (Krylov) approximation … (Saad 1992) [memory]" (ll. 55–58). The proof of (b) starts from it (ll. 129–130). The Status calls it "not claimed" (ll. 265–266), and the only check is numerical (D: quadrature against the exact transform to 8.7·10⁻¹⁵ at five points).
  - **What the citation covers.** Saad's setting is a finite matrix [snippet: secondary sources cite "Saad (1992, Thm. 5.1)" for an error formula of Krylov approximations to e^{A}v]. Here J is infinite and unbounded.
  - **An equivalent form.** In D61's s-form the identity says S =_d sT + (1 − s)Z, with s ~ Beta(K, K + 1), T the nodes' uniform Dirichlet average, Z ~ N(0, 1), all independent; the Beta density is s^{K−1}(1 − s)^K·K·C(2K, K).
  - **My exact check** (W/scripts/duhamel_exact.py). The moments of that mixture equal (ii)'s E S^{2j} = (2j)!(2K)!/(2K + 2j)!·P_j(K) exactly in rationals for K ≤ 8, j ≤ 6: 48 identities, 0 mismatches. The nodes enter through h_m(z) from e_{2i} = (−1)^i K!/(i!(K − 2i)!2^i).
  - **A proof** is in D3dv-1's fix.
- **(a).**
  - The representation He_K(√K(ζ + 1/ζ)) = K^{K/2}ζ^{−K}e^{Kζ²/2}I_K(ζ) holds by the shift t → t + iζ.
  - The Gaussian bound |·| ≤ e^{−K(1−r²)t²/2} holds.
  - The circle identity is (1 − cos ψ)(1 − |μ|²cos(2 arg μ + ψ)): I expanded both sides.
  - The prefactor (K − 1)!K^{1−K}e^K ~ (2πK)^{1/2} and the Laplace term combine to exactly 1 − μ².
  - The ellipse encloses the nodes because √K(ρ + 1/ρ) > 2√K > 2√(K − 1).
- **The contour bound.**
  - On |ζ| = ρ₀ the bound is max_θ Re f_σ(ρ₀e^{iθ}) − (1 + Re σ²/2 + log|σ|) + o(1). The o(1) is uniform on compact σ because I_K → (1 − ζ²)^{−1/2} uniformly on |ζ| ≤ ρ₀.
  - At ρ₀ = 1 the maximum is max_c(2ac − c² + ½) = M(a) + ½, which gives Λ⁺. The convergence as ρ₀ → 1 is uniform on compact σ.
  - Λ⁺ = h(|σ|) − (|Re σ| − 1)₊², with h(1) = h′(1) = 0, so F_η is C¹ across |σ| = 1 and across |Re σ| = 1. This also confirms V's gradient bound |∇C|² = h′² − 4(a − 1)(1 − a/ρ²) ≤ h′².
- **The band near |σ| = 1.** Choose r ∈ (max|η|/2, 1) on the compact set and r₀ with h(r₀) + δ(ρ₀) below a third of the margin. Then the strict margin at s = 1/|η| (Λ⁺ = 0 there and g_η(1/|η|) < g_η(½)) extends to the band. Ω is open (the supremum is continuous and tends to −∞ as |η| ↓ 1).
- **The Laplace constant.** K·C(2K, K)·2·4^{−K}(1 − η²/4)(π/(K(4 − η²)))^{1/2} → (1 − η²/4)^{1/2}, using Re G″(½) = −2(4 − Re η²) < 0.
- **The axes.**
  - Real axis: the derivative (w − 1)²/w ≥ 0. The gap and its first two derivatives vanish at 2, and the third, 2/(η − 1)³ − 4/η³, is positive for η < 2^{1/3}/(2^{1/3} − 1) ≈ 4.85.
  - Imaginary axis: concavity follows from 1/(y − 1)² ≥ 4/y² > 2/y² on (1, 2].
  - The md's −¼(2 − η)³ and −½(2 − y)² follow.
- **The domination closed form.** `gap_closed` equals g_η(s) + Λ⁺(ηs) − g_η(½) with w = 1 − s. The box bound's monotonicities hold: w² < ¼ < ½, M is increasing in |a|, and −log R is taken at R1.
- **Route (2)'s hypotheses.**
  - The chord [(1 − t)η/2, (1 + t)η/2] runs along u = η/|η|, and Re(φ″u²) = 2(Re η² − 4)/|η|² < 0.
  - On it, and on the 0-side, E = log(1 − x²) + x²Re c < 0 except at x = 0, inside |σ| ≤ (|η|/2)(1 + t) < 1.
  - VA's t = ε/4 < ε/(2 − ε). V's t = 0.9(2 − R2)/R2 ≤ 0.9(2 − R)/R.
  - The junction never overshoots the cap: κt ≤ θ·εω²/(2ε(1 − ε/4)) < θ/2.
  - Forms A and B overlap (q ≤ 2.5 or p ≤ 0.5), and the cap lies inside A, since on the cap q ≤ θκ ≤ θ cot θ ≤ 1.
  - The inner pieces cover their parameter domains by construction: the roots share edges exactly, splits use shared midpoints, and skips happen only when ρ_lo > rmax.
  - The connector carries δ₀ + iY(η) to δ₀ + iy_c. The polyline ends within END_R. The end piece goes to z = 1.
- **Uniformity on compact sets.** A compact K ⊂ {|η| < 2} splits into K ∩ {|η| ≤ 1.8927} ⊂ Ω and K ∩ {|η| ≥ 1.8927} ⊂ Ω′. Each piece is uniform. Ω′'s margins come from finitely many boxes, with E ≤ X²·bd (route 'v') or E ≤ X³·sup f′ (route 'd').
- **(iv).**
  - The tails: |ε̃_K|²e^{−Ax²−By²} ≤ 4e^{−(u−1)|ξ|²/2}.
  - The mass at 2K/γ(θ) ≤ 4Kr < c²K, using ln u > 1 − 1/u.
  - The constants: off the cusp r^{2K}(1 − r)/(a(2γx(γy − γx))^{1/2}); at the cusp (πK)^{1/2}r^{2K}₂F₁(−½, −½; 1; r²)·2r/a, through (K!)²/(2K)! ~ (πK)^{1/2}4^{−K}.
  - All re-derived; all correct.

**(4) Seams in the report** (W/scripts/seam.py, W/logs/seam.log) [own]

| plant on the raw | what V and D print | VA's conclusion "on the whole open disc \|η\| < 2" |
|---|---|---|
| V's last band failed | V: "holds on \|η\| ≤ 1.99" | still printed |
| V's last band ending at 1.994 | V: "holds on \|η\| ≤ 1.994" | still printed |
| D's branch and bound unresolved | D: "no radius beyond \|η\| ≤ 1 is verified" | still printed |
| VA's R1 = 1.996 (meta ε₀ 0.004) | — | still printed, and "holds on 1.995 ≤ \|η\| < 2" still printed (hard-coded) |

**(5) Own mutants** [own: W/scripts/spec127_la.py, spec127_gp.py; W/mut/res_la.jsonl, res_gp.jsonl]

gh_lastann.py, `-k "part_va or md_is_the_report"`:

| id | mutant | result |
|---|---|---|
| Z01 | end link < 0.5 | killed |
| Z02 | end link < 2·END_R | **survives, unsound** |
| Z03 | `tiled` without `inside` | **survives, unsound** |
| Z04 | precondition allows f1 ≥ −1 | **survives, unsound in principle** |
| Z05 | precondition allows f2 ≤ 2 | **survives, unsound in principle** |
| Z06 | precondition allows R1 > −1 | **survives, unsound in principle** |
| Z07 | budget removed | survives, sound |
| Z08 | budget 50 | survives, sound; never exercised |
| Z09 | exact area as ≥ | survives, equivalent given inside and disjoint |
| Z10 | soundness tolerance 1e-6 | survives (round 122's S05) |
| Z11 | λ's corners dropped | killed |
| Z12 | end link at pts[−2] | killed |
| Z13 | end link in floats | survives, equivalent up to rounding |
| Z14 | random points only | killed |
| Z15 | END_R = 0.03 | killed by the raw |

exp_gh_peano.py, the full test file:

| id | mutant | result |
|---|---|---|
| G01 | `va_lines` without the end piece's `ok` | **survives, unsound**: the 'end' record has no 'failed' key |
| G02 | without the soundness failures | survives, equivalent: `ok`'s first clause already requires every record's failed == 0 |
| G03 | without the handoff failures | survives, equivalent, for the same reason |

**(6) Literature.** Web searches only; arxiv.org was not reachable from this session, so for arXiv I saw snippets only (see Sources).

## Verdicts per part

| part | verdict |
|---|---|
| (i), (ii), (iii)(c), (v) | Unchanged; not re-refereed beyond what (iii) and (iv) use. (ii)'s P_j is still "a sketch with routine steps", as stated. |
| (iii)(a) | Re-derived: proved as written. |
| (iii)(b), route (1) (Ω; D) | Sound. D's branch and bound recomputes equal to the raw in balls, and the closed form and box bound are valid. The only gap is that Duhamel's identity is not proved in the text (D3dv-1). |
| (iii)(b), route (2) (Ω′; V, VA) | Sound. The hypotheses match what V and VA certify. The seams at 1.8927 and 1.995 are the same floats. The reflections, arg η = 0 and arg η = π/2 are covered, and uniformity on compact sets holds. |
| Part V | Recomputed equal to the raw in balls; consistent with D61's (2). |
| Part VA (with round 126) | Sound as shipped, and recomputed equal. The tiling is exact, and the end link and `_y_enclosed`'s guards are correct. Some plants are too loose (D3dv-3), and the conclusion is not seam-checked (D3dv-2). |
| (iv) | The Laplace step, tails and constants are re-derived and correct, given (iii)(b) on \|η\| ≤ c < 2. |
| md, D61 wording | Every number I checked is in the raw. Cosmetic overstatements remain (D3dv-4, -5). |
| Tests | 42 pass; R05, the P8 plant, T1 and W01 are now in. Unsound survivors: Z02–Z06 and G01. |

## Status of the round-126 issues

- **D3du-1 (the tiling to 10⁻¹²): resolved.** The rational area is exact, and plant T1 is in and kills R04. Residue: Z03 (`tiled` without `inside`) survives, because the box-outside plant also changes the area (D3dv-3).
- **D3du-2 (the end link): resolved.** The link is checked in balls, and P8 kills R03. Residue: P8 ends 0.5 or more away, so the threshold is unpinned and Z02 (2·END_R) survives (D3dv-3).
- **D3du-3 (failure paths unplanted): mostly resolved.**
  - In: the connector's last piece; the whole connector; a middle half-piece; `_y_enclosed` → None (76 failures); lo raised; a box outside; W01. R01 and R02 are killed.
  - Residue: Z03, and the preconditions' other sides (Z04–Z06).
- **D3du-4 (the true π/2): resolved.** The top-strip plant kills R05, and the log's C01 line is corrected. The fix text's "fast either way with the budget" is inexact: the depth limit makes it fast, and the budget is never reached (Z08).
- **D3du-5 (soundness corners): resolved in code.** Residue: the printed "distinct" count for A-FC is point–branch evaluations (D3dv-4).
- **D3du-6 (code and text): resolved.**
  - In: the preconditions, the budget, the docstring, HANDOFF §5's paragraph.
  - Residue: HANDOFF l. 62's part list lacks VC and VA; l. 317 says "Arb assumed correct"; the "35 500 distinct points" figure (D3dv-6).

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3dv-1 | minor | **Duhamel's identity, the starting point of (iii)(b), is not proved.**<br>• D61 introduces it as the Krylov error representation "(Saad 1992) [memory]" (ll. 55–58), uses it in the proof of (b) (ll. 129–130), and calls it "not claimed" (ll. 265–266) [md].<br>• Saad's result concerns finite matrices [snippet]; J is infinite and unbounded.<br>• The only check is numerical (D, five points, ≤ 8.7·10⁻¹⁵ [md]).<br>• It is true: in the s-form it is S =_d sT + (1 − s)Z with s ~ Beta(K, K + 1), and the moments agree exactly with (ii) for K ≤ 8, j ≤ 6 (48 identities, 0 mismatches) [own: duhamel_exact.py]. | **D61 Proof (iii), before (a):** "Duhamel's identity. With J₀ = J_K ⊕ J′ and V = J − J₀ (the entry √K between the levels K − 1 and K), J^m − J₀^m = Σ_{j<m}J^{m−1−j}VJ₀^j entrywise (finite sums: the matrices are banded), and J₀^je₀ = J_K^je₀ lies below level K, so (E − Q_K)[z^m] = √KΣ_{j<m}(J^{m−1−j})_{0K}(J_K^j)_{K−1,0}, (J^n)_{0K} = E[Z^nHe_K(Z)]/√K!, (J_K^j)_{K−1,0} = Q_K[z^jHe_{K−1}]/√(K−1)!; multiply by ξ^m/m! and sum (absolutely: E e^{\|ξ\|\|Z\|}\|He_K(Z)\| < ∞), with ∫_0^ξ(ξ − u)^nu^jdu = n!j!ξ^{n+j+1}/(n+j+1)! and E[e^{wZ}He_K(Z)] = w^Ke^{w²/2}; equivalently S =_d sT + (1 − s)Z, s ~ Beta(K, K + 1)."<br>**D61 Status:** "Duhamel's identity proved (the Krylov representation's infinite-matrix form)"; keep Saad as context.<br>**Tests:** pin the Beta-mixture moments for K ≤ 8. |
| D3dv-2 | minor | **VA's conclusion is not tied to the seams.** `va_lines` prints "So the valley condition holds on 1.995 ≤ \|η\| < 2 … and, with V and D's branch and bound, on the whole open disc \|η\| < 2" from VA's records alone.<br>• It stays printed when V's last band fails (V then says \|η\| ≤ 1.99), when V ends at 1.994, when D is unresolved, and when VA starts at 1.996 (its "1.995" is hard-coded) [own: seam.py].<br>• It also does not depend on the end piece's `ok`: G01 survives, because the 'end' record has no 'failed' key [own].<br>• On the shipped data it is true; only the tests' pins (CERT_BANDS[−1] == (1.99, 1.995), O['R1'] == 1.995, DOM_BB[−1] == CERT_BANDS[0][0]) guard it. | **exp_gh_peano.py `report`/`va_lines`:** pass V and D. Print the whole-disc sentence only if D's bb is resolved in balls, v_lines' contiguous R_v == 2 − meta.eps0 (exact float), and every piece incl. end passes; print the annulus from meta.eps0.<br>**Tests:** the four seam plants above and an `end` record with ok=False, each expecting no conclusion. |
| D3dv-3 | minor | **Round 126's plants leave unsound survivors** [own: res_la.jsonl].<br>• Z02 (end link at 2·END_R): P8 ends 0.5 or more from z = 1, so a polyline ending at 0.03 would pass while `piece_end` covers only 0.02.<br>• Z03 (`tiled` without `inside`): the box-outside plant changes the area, so `tiled` is False anyway. An area-preserving shift of a bottom-row box by 2⁻¹⁴ (inside False, disjoint True, area exact) is stopped only by `inside` [own: tile_plant.py].<br>• Z04–Z06: only the R2 > 2 side of `_y_enclosed`'s precondition is planted.<br>• Z07 and Z08 (the budget removed, or 50) survive; this is sound, since the budget never binds. | **test_gh_peano.py `test_part_va_…`:**<br>• P7: the last vertex and its piece's end moved to 0.97 + 0j, expecting 'polyline end'.<br>• The 2⁻¹⁴ shift plant, `assert not t['tiled']`.<br>• `_y_enclosed(1.995, 2.0, −1e-3, 0.01, …)`, `(…, 1.5, π/2 + 1e-3, …)` and `(0.0, 1.0, …)` → None.<br>• Optionally a budget plant (budget=4 on the NaN box → None).<br>**REVIEW_LOG round 126:** "fast either way" — by the depth, not the budget. |
| D3dv-4 | cosmetic | **"Distinct seeded points".** In A-FC, 668 of 2011 boxes carry both branches, and `soundness` counts each point once per branch: 24 770 evaluations, 18 326 distinct points. The four pieces' total is 29 056 distinct points, not HANDOFF's "35 500 distinct points" [own: sound_count.py; md]. | `soundness`: return `points` (distinct) and `evals`. md/D61 Numerics: "18326 distinct points (24770 point–branch evaluations)". HANDOFF: the correct total. |
| D3dv-5 | cosmetic | **Wording in D61 and the md:**<br>(a) "Λ⁺ bounds its rate (attained, numerically and as K → ∞ …: fits at K ≤ 340 approach it without reaching it)": attainment is a heuristic, and the sentence contradicts itself.<br>(b) "of exactly its area" is the float rectangle [1.995, 2] × [0, fl(π/2)]; the sliver to the true π/2 is carried by the top boxes.<br>(c) "b = (1 + ε/4)η/2 (part V's chord point" — V's b is at t = 0.9(2 − R2)/R2; VA's is another point of the same chord.<br>(d) "verified again in rounds 122 (…) and 126": round 126 checked the joint, the tiling and `_y_enclosed`, not the forms.<br>(e) "(both ends of the expansion variable)": the lower end is used only when > 0.<br>(f) The complex-phase Laplace method and Golub–Welsch are used without a reference [md]. | D61 (iii)(b) and md:<br>(a) "conjecturally attained (fits approach it)".<br>(b) "exactly its area in rationals, the top edge then taken at the true π/2".<br>(c) "a point of part V's chord".<br>(d) "rechecked (joint, tiling) in round 126".<br>(e) "and at its lower end where positive".<br>(f) Cite a Laplace-method source (e.g. O'Sullivan 2019 [snippet], or Olver's textbook [memory]) and Golub–Welsch 1969 [memory]. |
| D3dv-6 | cosmetic | **HANDOFF residue:**<br>• l. 62 lists `exp_gh_peano.py`'s parts as "P, M, L, LC, S, R, Q, D, V, DB", without VC and VA;<br>• l. 317 says "Arb assumed correct", where D61 assumes Arb's, python-flint's and CPython's float arithmetic and `math.nextafter`;<br>• l. 197, "35 500 distinct points" (see D3dv-4) [md]. | HANDOFF: the part list with VC and VA; D61's assumption sentence; the corrected count. |

## Sources

**Repository** [md][raw] at 31487a6:
- review/drafts/D61.md;
- experiments/exp_gh_peano.py and experiments/gh_lastann.py;
- results/chaos/gh_peano.md and gh_peano_raw.npz;
- tests/test_gh_peano.py;
- tools/mutants.py and tools/mutants/round126_D61_own.py;
- review/REVIEW_round126_D3du.md and REVIEW_LOG.md (rounds 120–126);
- HANDOFF.md (§§ on D‴.61).

**My scripts and logs** [own], all under W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref127:
- scripts: recompute.py, seam.py, sound_count.py, tile_plant.py, duhamel_exact.py, spec127_la.py, spec127_gp.py;
- logs: pytest_gh.log, recompute.log, seam.log, sound_count.log, tile_plant.log, duhamel_exact.log, mut_shipped.log, mut_la.log, mut_gp.log;
- mutant results: mut/res_shipped.jsonl, mut/res_la.jsonl, mut/res_gp.jsonl.

**Literature.** Web search; arXiv pages could not be fetched from this session, so for those I saw snippets only.
- Convergence analysis of Hermite approximations for analytic functions, arXiv:2312.07940 [snippet]: root-exponential Gauss–Hermite convergence, the bound K·e^{−2ρ√(2n)}, said to be sharp. It is a rate bound, not P's limit.
- Chawla–Jain, Asymptotic error estimates for the Gauss quadrature formula, Math. Comp. 22 (1968) (AMS page in results) [snippet]: contour asymptotics via Barrett's Q_n/P_n.
- P. Liu, Asymptotic Error Analysis of Gauss Quadrature for Nonsmooth Functions, arXiv:2509.23532 (2025) [snippet]: Gauss–Legendre, contour-integral remainder, the leading term in cos((2n + 1)φ).
- Jia–Lv, A posteriori error estimates of Krylov subspace approximations to matrix functions, arXiv:1307.7219 [snippet]: builds on "the error expansion of the Arnoldi approximation … due to Saad". A secondary source cites Saad (1992, Thm. 5.1), SIAM J. Numer. Anal. 29, 209–228 [snippet]. These confirm Saad's setting is finite matrices; I did not see the corner-entry statement itself.
- Bennett–Howls–Nemes–Olde Daalhuis, globally valid error bounds for steepest descent, arXiv:1710.10073 [snippet]; O'Sullivan, Perron's saddle-point method, Pacific J. Math. 298 (2019) [snippet]. Rigorous saddle-point and Laplace-method references for D3dv-5(f).

**Not found:** in five searches I found no prior statement of conjecture P's limit, ĝ_K(√Kη)e^{−Kη²/4} → (1 − η²/4)^{1/2}, nor of the Peano kernel's Beta-mixture law. I also found no computer-assisted, validated-numerics proof of a comparable Laplace-transform or saddle-point limit. The closest are analytic asymptotics of Gauss rules (Chawla–Jain; Liu) and rigorous steepest-descent remainders (Bennett et al.).

## Lesson

1. **Check the seams where the conclusion is printed.** A theorem assembled from three certificates is only as strong as the code that joins them. Here the seams at 1.8927 and 1.995 are pinned in tests, but the sentence "on the whole open disc" is printed without looking at V, D or VA's own radius. Put the seam check where the conclusion is printed.
2. **Plant at the threshold.** A plant 26 times past END_R does not pin END_R. A plant that moves the area as well as containment cannot tell whether `tiled` uses `inside`.
3. **Write out the elementary step the proof starts from.** A step called "not claimed" can still be load-bearing. Duhamel's identity for the infinite Jacobi matrix takes five lines of telescoping on banded matrices, and those five lines, not a finite-matrix citation from memory, are what make (iii)(b) a proof.
