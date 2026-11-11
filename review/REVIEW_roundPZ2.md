# Referee report PZ2: branch pairs, dictionary spacing, the Euclidean metric (WP1)

(Archived verbatim from the referee subagent's final message; the referee worked on commit a00b405 in its own copy and did not modify the repository.)

## Scope

Reviewed state: commit **a00b405** (branch claude/lyapunov-context-setup-yg51ew), exported with `git archive` into W/repo, W = scratchpad/refPZ2.
- The repository was not modified. Its venv was used read-only. A symlink to it inside W/repo lets `tools/mutants.py` find `venv/bin/python`.
- Every run used `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 10`, PZ_OUT pointing into W/out, and no network. At most two heavy processes ran at once. The other session job kept the load at 4-7 on 4 cores. Nothing was written outside W. I confirmed afterwards that no process of mine is left.
- Read first: review/REVIEW_round156_PZ1.md.
- Read in full: `src/core/factorized_rbf.py`, `experiments/exp_parametrization.py` (K, L and what they call, including `margins`; `path_duals_nd`, `BumpsND`), md sections K and L, both test files, `tools/mutants/wp1_factorized_own.py`, `tools/mutants.py`, Paper 17's Lemma (o), the "starts on branches" proposition and its supplement proof, and `git diff 5d75dc1 a00b405 -- manuscript`.
- CLAIMS.md section 5 was read last. SESSION_NOTES.md was read after that.

Not done:
- No pdflatex, so the tex was read and checked mechanically.
- No literature search.
- Parts A-J were not rerun.
- My independent engine reuses the study's reference path (`P.reference`) for horizons and reference sets. PZ1 verified that path against the stored orders, and I checked its duals indirectly (item 3).

## Summary

**K1 (PaperPairs, time scale, common early factor): PARTLY CONFIRMED.**
- The pair start, the scale ln(1/ε)/(γ a_max tanh(a0/a_max)) and the factor 1 − c/L_ε are right. My c (quadrature) is 6.2085 for a0 = 0.1 and 4.0524 for a0 = 0.3. The measured and predicted errors agree to 0.0117 on the 31 cells. All 110 K cells reproduce with an independent engine.
- "Either unit" does not flatter the pairs.
- "Every entry is early by this common factor, derived to leading order in 1/L_ε" is exact only for the first entrant. It holds for later atoms only while their duals are straight. On 14 coupled targets it fails by a factor of up to 4 (PZ2-1).

**K2 (pairs between depth 2 and 3): CONFIRMED.**
- An a0-scan on the design is monotone from P0-like (a0 → 0) to D2-like (a0 ≥ 1). The "5 of 11" and "9 of 11" figures reproduce.
- Wording: the displacement comparison is medians only, and 0.0217 is not ≤ 0.02 (PZ2-7).
- The data are pre-asymptotic. At smaller ε (own cells) the set count rises to 8-9 of 11 and exactness to 11 of 11 (PZ2-2).

**L1 (dictionary spacing): CONFIRMED.**
- All cells reproduce: 3 jobs bit for bit, 13 cells with an independent engine, 15 Radau cells included.
- `margin_never` is independently verified: P0's dormant coefficients match the creep law with μ never to 1.000 at ε = 1e-12, in all 9 configurations.
- The rule "D2 passes iff μ > ln(1/θ)/ln(1/ε)" holds out of sample. In 6 dictionaries × 43 values of ε, the observed transition is at 0.90-1.02 of the rule's L. All 36 new D2 cells are consistent with it.
- Qualification: the on-grid arm was run only where the rule says D2 passes (PZ2-5).

**E1 (Euclidean metric in Paper 17): CONFIRMED.**
- Derived by hand and checked by finite differences at 3 constant sets with G ≠ 1 and a_max ≠ 1.
- The tex edits are brace- and ref-clean.
- Defects around it: the committed PDFs still carry the withdrawn sentence (PZ2-4), and the letter ν clashes with the paper's margins (PZ2-11).

**M1 (mutation run): CONFIRMED.**
- I re-ran all 123 mutants with the repository's tool. 122 were killed by assertion failures and S14 survived. S14 is equivalent: 0 differences in 400 random records with NaN times.
- All 12 pass-1 survivors are killed. The unmutated control passes (66 tests).
- Remarks: P12 is behaviourally equivalent and is killed by an attribute assertion, N01 by a constant pin, and 21 functions have no mutant at all (PZ2-9).

Numbered findings:
- PZ2-1 (major): "every entry early by one common factor" is exact only for the first entrant.
- PZ2-2 (minor): part K tests a sliding extension of a fixed-centre proposition, at pre-asymptotic ε.
- PZ2-3 (minor): CLAIMS.md section 5 and SESSION_NOTES.md predate a00b405.
- PZ2-4 (minor): the committed PDFs are stale.
- PZ2-5 (minor): the on-grid arm covers only the spacings where D2 passes.
- PZ2-6 (minor): two margin definitions, and a NaN convention for the rule.
- PZ2-7 (minor): K2's dormant-displacement comparison uses medians.
- PZ2-8 (cosmetic): "between depth 2 and 3" rests on two a0 and on medians.
- PZ2-9 (minor): the mutation claim's wording and its coverage.
- PZ2-10 (cosmetic): the md header gives the Radau count for C only.
- PZ2-11 (cosmetic): notation clash for ν.

## What I checked and how

**1. Pair time scale, factor, "either unit".**

*Time scale and c.*
- Own derivation from β = ε + 2v_s(sinh²(α/a_max) − sinh²(a0/a_max)) and dη = dβ/(γaβ), with a constant drive g = 1. Three routes: direct DOP853 integration of the law, quadrature in α, quadrature in ln β. They agree to 8 digits [own: derive_c.py].
- Ratio of the formation time to the wall time:

| a0 | ε | ratio | 1 − ratio |
|---|---|---|---|
| 0.1 | 1e-6 | 0.55108 | 0.4489 |
| 0.1 | 1e-12 | 0.77531 | 0.2247 |
| 0.3 | 1e-6 | 0.70674 | 0.2933 |
| 0.3 | 1e-12 | 0.85334 | 0.1467 |

- The md prints 0.449, 0.225, 0.293 and 0.147 [md K]. The limit c(ε → 0) is 6.2085 (a0 = 0.1) and 4.0524 (a0 = 0.3). The finite-ε values are 6.2020 and 4.0516 at 1e-6.
- The first entrant is early by exactly this factor on every target I ran. Measured earliness is 0.4480-0.4518 against 0.4481 at 1e-6, 0.1121-0.1131 against 0.1122 at 1e-24, and 0.0561-0.0565 against 0.0561 at 1e-48 [own: coupled2.py].
- The measured error is c/L_ε all the way to ε = 1e-96 (0.0282 against 6.2/221), so there is no hidden constant offset in the scale.

*Why "every entry" fails.*
- The formed level is crossed at dual level κ_f = 1 − c/L, not at the wall. A dual that is straight from τ = 0 gives t_j = κ_f τ_j. After an earlier entry the dual's slope changes, and the earliness of atom j becomes slope-dependent at the same order 1/L.
- Test: the 14 sample targets of part C with the closest features (family Db, K = 2), pairs a0 = 0.1 [own: coupled2.py].
- Second atom, opposite-signed pairs (6 targets). Earliness equals the prediction; the ratio error is ≤ 0.011 at 1e-6.
- Second atom, same-signed pairs (8 targets, As = [1, 1]). Earliness against prediction:

| ε | measured | predicted |
|---|---|---|
| 1e-6 | 0.618-0.832 | 0.448 |
| 1e-24 | 0.353-0.507 | 0.112 |
| 1e-48 | 0.183-0.317 | 0.056 |

  The ratio error (entry-time ratio against the path's) is 0.31-0.70, 0.27-0.44 and 0.13-0.28 at the three ε.
- Mechanism, target 121: both atoms pass the formed level at τ = 3.05 and 3.35. The path enters them at 5.53 and 19.86 because the second dual is at 0.91 of its wall at the first entry and then crawls.
- The md's own part C shows the same tail for D2. The 95th percentile of the τ error is 0.654 against a median of 0.339 at 1e-6, and 0.42 against 0.17 at 1e-12 [md C].
- The 31 cells of part K are all well-separated family-A "nearest" targets, so the law looks exact there. The md says so ("the easy family"). The docstring of `form_ratio` and the md text do not.

*"Either unit".*
- `one_each` is true in 44 of 44 pair cells [raw] and in all 66 extra cells I ran (n_formed = n_atoms) [own].
- On 9 settings (3 targets × 3 cells) the partner's |w| is at most 1e-10 of the formed unit's, and at most 2.4e-8 absolute against a threshold of 0.009. "Both units" forms nothing, "pair sum" and "either" give identical sets [own: partner.py].
- Where `set` is true, `exact` is also true. So no stricter definition changes K2. The threshold-free gap K, which carries the ordering, does not depend on the definition anyway.

**2. E1 and the tex.**

*Derivation.* ∂αF = −gβa′/m, ∂βF = −ga/m and ∂cF = −w⟨r∇cφ⟩ give α̇ = Gm(−∂αF), β̇ = γmβ(−∂βF) and ċ = Gm(−∂cF). With ν = β^{1/2}: ∂νF = 2ν∂βF, ν̇ = β̇/(2ν) = (γm/4)(−∂νF), and (γmβ)⁻¹dβ² = (γm/4)⁻¹dν². So the rates are (Gm, γm/4, Gm), and lr_ν/lr_α = γ/(4G) = v_s/a_max² by v_s = γa_max²/(4G) (paper17.tex, Setting) [snippet]. Fisher-Rao on (0, ∞) is flat in √β [memory].

*Numerics* [own: e1_check.py]. Own RHS from eqs. (free) and (centres), central differences of F(α, ν, c):

| (m, G, γ, a_max) | agreement of α̇, ν̇, ċ with −lr·∂F |
|---|---|
| (12, 1, 0.1, 1) | 1.8e-9 to 2.9e-9 |
| (5, 2.3, 0.37, 0.7) | 4.7e-10 to 7.6e-10 |
| (0.8, 0.45, 3.1, 1.9) | 1.3e-10 to 2.4e-10 |

- Wrong-rate controls fail by 0.5 to 1.3 relative.
- The metric identity holds, and Q = β − v_s cosh(2α/a_max) is conserved along a driven trajectory (spread ≤ 4.4e-15).

*Tex.*
- The 4 added lines are brace-balanced with even `$` parity. The whole files are balanced. The refs used (lem:curve, eq:free, eq:centres, prop:unit, prop:branch, prop:creep, thm:sliding and the S- refs) exist. The macros are the paper's [own: tex_check.py].
- The README row D3eu-1 is consistent with the text.
- Stale PDFs and the ν clash: PZ2-4, PZ2-11.

**3. Part L.**

*Recomputation.*
- Bit for bit (every field but `secs`) with the committed code: jobs L (13, off), (9, on) and (21, off), 8 records each [own: bitident_kl.py].
- 13 cells with the independent engine, 6 of them Radau-fallback cells (on-grid P0 and D3, off-grid D3), all agree on n_formed, gap_k to 3 decimals, exactness and rms [own: l_own.py].

*`margins` is correct.*
- P0's dormant max |w| at ε = 1e-12 [raw], against a_max ε^{3/2} cosθ/(m(2v_s)^{1/2} sin³θ) with θ = πμ/2.
- With μ never the ratio is 1.000 in all 9 (nd, variant) configurations. With the later-entrant μ it is 0.004-0.48 where the two margins differ.

*Out of sample: more cells.* 36 new D2 cells, 21 predicted to fail and 13 to pass, all as predicted. The 2 remaining cells (3 kernels) have μ never = NaN; D2 passes, which is the "no never-entering kernel, so pass" convention (PZ2-6) [own: l_extra.py]. Selected results:

| dictionary | ε | result |
|---|---|---|
| 25 kernels, off-grid | 1e-10 | 8 units, fails |
| 31 kernels, off-grid | 1e-48 | 3 units, passes |
| 21 kernels, off-grid | 1e-36 | passes |
| 5, 4, 3 kernels, off-grid | 1e-6 | pass (5 kernels not exact at 1× horizon) |
| 41 kernels, P0 and D3 | 1e-6, 1e-12 | set and exact; gap 1.95 and 2.51 decades at 1e-6 |
| 21 kernels, on-grid, D2 | 1e-6 and 1e-12 | 9 and 6 units, fails |
| 16 kernels, on-grid, D2 | 1e-6 | fails |

*Out of sample: where D2 starts to pass* [own: rule_scan.py]. A scan of ε = 10^-6 … 10^-90 gives the first ε from which D2 forms the path's set. L_obs/L_rule is:

| kernels | 11 | 13 | 16 | 21 | 25 | 31 |
|---|---|---|---|---|---|---|
| L_obs/L_rule | 1.02 | 0.90 | 1.00 | 0.93 | 0.97 | 1.02 |

The rule has no fitted parameter. Its "iff" is correct to about ±10% in L.

*Part D with μ never* [own: partd_rule.py]. The rule agrees on 21, 21, 21, 20, 20, 21 of 21, against 21, 21, 20, 19, 20, 21 with Part D's margin.

*On-grid.* It is a legitimate control (it separates the margin from the sliding), but it is a different problem. The features are exactly in the dictionary, so no sliding is needed. It was run only where the rule predicts D2 passes (PZ2-5).

**4. Integration, statuses, Radau.**

*Independent engine.* myfr.py: own RHS written from the paper's formulas, Radau at rtol 1e-10, DOP853 for three cells [own].
- All 110 K cells (44 pairs, 66 controls) were recomputed. 0 mismatches in formed atoms, unit count and exactness; the largest relative difference in gap_k is 1.7e-6.
- DOP853 reproduces pairs a0 = 0.3 at 1e-6 (6 atoms), D2 1e-6 and D2 1e-12.
- The 5 K cells that used the Radau fallback (D3 at 1e-12 on targets 1, 3, 8, 10, and P0 at 1e-6 on target 3) agree.

*Statuses.* All K and L cells are `ok`; B has one recorded failure (D5 at 1e-10) [raw, md]. None of my cells needed a status the study lacks.

*Radau usage.* K: 5 of 110 cells. L: 15 of 72 (2 of them with two pieces). C: 104 of 1200 [md header]. K and L counts are in the raw only (PZ2-10).

**5. Reproducibility.**
- The md regenerated twice from the raw with `--from-raw` is byte-identical to the committed md (md5 decf4dca47bdd9cb8565cf0a8243487b) [own].
- Tests: 66 passed in 17 s (test_factorized_rbf.py) and 26 passed in 54 s (test_parametrization.py).
- Bit for bit (every field but `secs`): K jobs for targets 3, 8, 10 (the Radau cells included) and L jobs (13, off), (9, on), (21, off) [own].

**6. Mutants.**
- All 123 mutants of the spec were re-run with `tools/mutants.py` (1,488 s in total). The 12 pass-1 survivors M31, M33, M46, M48, M55, M66, M67, N05, N07, N05b, S15 and P12 are all killed.
- Count: 122 killed, S14 survives. In every kill pytest reports "1 failed, N passed", with no timeout and no collection error. The unmutated control passes (66 tests).
- S14 is equivalent: on 400 random records with NaN formation times at arbitrary positions, original and mutant give identical output, with warnings turned into errors [own: s14_equiv.py].

*Soft kills.*
- P12: PaperPairs never reads `amp0`, so the mutant is behaviourally equivalent. It is killed by `assert ro.amp0 == 0.0` (test line 462).
- N01: killed by `assert FR.STIFF_BUDGET == 150000`.
- Wording: "killed by tests + the experiment's checks" is not what the tool does. Its kill criterion is test_factorized_rbf.py only.

*Coverage* [own: mut_cover.py]. 21 functions and methods have no mutant, among them all of PowerN (D4 and D5), `_unpack`, `Unit`, `pts` and `Linear.init`.
- 10 extra mutants of mine in those functions: 6 killed, 4 survive.
- Survivors: Q04 (PowerN.init with s = ε), Q05 (assert N ≥ 1), Q09 (Linear.init random sign, unused), Q10 (equivalent).
- Only Q04 is substantive: no engine test pins PowerN's start, so the D4/D5 claims rest on the raw alone.
- 8 mutants of the new K/L analysis code (checked against tests/test_parametrization.py after adding the two data files the tests need to the slot): 6 killed, 2 survive.
- E01 (cell_atoms: the `n_atoms == K` clause) is redundant given the unit count. E04 (a fixed 0.009 for the per-target prediction threshold) is invisible on the design and within the 0.02 tolerance elsewhere.

**7. Overreach, CLAIMS.md against the md.**
- Section 5 (mtime 13:17) predates the commit (15:56). It contains no K or L numbers and says "PZ1-9 branch pairs — not done", "the dictionary-spacing sweep is the referee's, not committed", "mutation run not repeated", and that Paper 17 "needs a correction" (D8). All of that is now out of date (PZ2-3).
- The numbers it does quote match the md verbatim: part C (187/197, 192/198, 166/191, 176/195, gaps 5.4/8.4 and 4.5/7.3, displacements 0.0357σ/0.000567σ and 0.0188σ/0.000299σ, 2.57σ, 104 Radau cells), part D (21, 21, 20, 19, 20, 21 of 21; thresholds 0.33…0.014; top-K 0, 1, 10, 13, 20, 21), parts E, H, I and J.
- Is "between depth 2 and 3" an artefact of two a0? No. On the design [own: a0_scan.py], gap K decreases monotonically with a0. At ε = 1e-6: a0 = 0.003, 0.01, 0.03, 0.1 give 3 atoms and gap_k 5.66, 4.97, 3.84, 1.89; a0 = 0.2, 0.3, 0.5, 1, 3 give 5, 6, 8, 9, 9 atoms and gap_k 1.03, 0.71, 0.45, 0.12, 0.125. The reference values are P0 6.15, D3 5.27 and D2 0.144 (9 atoms). At ε = 1e-12, a0 ≥ 1 gives 6 atoms and gap_k 0.85 and 0.65, against D2's 0.68. The family interpolates from P0 to D2, and a0 is the interpolation parameter.
- Smaller ε for pairs (11 targets, own cells, not committed) [own: k_ext.py]:

| ε | a0 = 0.1: set / exact | a0 = 0.3: set / exact |
|---|---|---|
| 1e-6 [raw] | 2 / 5 of 11 | 0 / 0 of 11 |
| 1e-12 [raw] | 5 / 9 of 11 | 2 / 5 of 11 |
| 1e-24 | 8 / 11 of 11 | 5 / 9 of 11 |
| 1e-48 | 8 / 11 of 11 | 6 / 10 of 11 |
| 1e-96 | 9 / 11 of 11 | 8 / 11 of 11 |

  The median gap K goes from 0.97 to 8.54 decades (a0 = 0.1) and from 0.38 to 8.94 (a0 = 0.3). The median dormant displacement falls to 0.000σ. The persistent mismatches are the three smallest margins (0.0013, 0.0057, 0.0133). The common-factor error is ≤ 0.003 and the ratio error ≤ 0.006 there.

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| PZ2-1 | major | "The common factor by which **every** entry is early … to leading order in 1/L_ε" (`PaperPairs.form_ratio` docstring; md K "predicted"; md D for D2) is exact only for the first entrant. The argument "walls scaled by one factor give the same path in scaled time" scales the *entry* levels, but the formed level (w = θ) is crossed earlier, at dual level 1 − c/L, not at the wall. Atom j is early by c/L only if its dual is straight until its wall; otherwise its earliness is slope-dependent at the same order 1/L. [own] 14 coupled part-C targets (Db, K = 2), pairs a0 = 0.1: first entrant early by 0.4481 as predicted on 14 of 14; second atom of 6 opposite-sign pairs within 0.011; second atom of 8 same-sign pairs early by 0.618-0.832 against 0.448 at 1e-6, 0.353-0.507 against 0.112 at 1e-24, 0.183-0.317 against 0.056 at 1e-48 (its c is about 11-23, not 6.2). Target 121: both atoms formed at τ = 3.05 and 3.35, entered by the path at 5.53 and 19.86. The part K sample (11 well-separated family-A targets, ratio error ≤ 0.0306) could not show it. [md C] D2's 95th-percentile τ error is 0.654 against a median of 0.339, the same tail. | `form_ratio` docstring and md K/D text: state that the first entrant is early by exactly 1 − c/L (derivation), and later atoms only while their duals are straight (separated features). Drop "every entry". Add coupled targets (C's Db 120-125, 130, 131) to part K and print τ ratio err. per target class. K1's "agree to 0.012 / 0.031" belongs to the stated sample. |
| PZ2-2 | minor | Part K compares the pair network with the *sliding* path of the atoms. The proposition "starts on branches" (paper17.tex l.272; supplement §pf:branch) concludes "the conclusion of Theorem iss", the fixed-centre theorem; the sliding limit is proved only for α(0) = 0 starts. The md does not say that the sliding counterpart is a conjecture. ε = 1e-6 and 1e-12 are L_ε = 13.8 and 27.6, with c/L = 0.45 and 0.22: pre-asymptotic. Set 2/11 and 5/11 (a0 = 0.1) are not evidence against the limit; the own cells above show 8, 8, 9 of 11 at 1e-24, 1e-48, 1e-96, exact 11 of 11. | md K intro: one sentence that the sliding extension is tested, not proved, and that ε = 1e-12 has not reached the regime. Extend K to ε = 1e-24, 1e-48, 1e-96 (66 cells, about 4 min; my cells are in W/out/k_ext.jsonl, script W/scripts/k_ext.py). |
| PZ2-3 | minor | CLAIMS.md section 5 (mtime 13:17) and SESSION_NOTES.md "Next, in order" (item 1; the state table at its top is current) predate a00b405. Stale: "PZ1-9 branch pairs — not done"; "the dictionary-spacing sweep is the referee's, not committed"; "mutation run not repeated"; D8 "Paper 17 needs a correction"; "on-grid features … (referee)"; "random starts … not in the committed experiment". SESSION_NOTES still words D2's 4.74/ln(1/ε) as a law. | Update section 5: add K and L with the numbers in this report, the metric correction as applied, and the 123-mutant run (122 killed, S14 equivalent). Fix the "Next" list in SESSION_NOTES. |
| PZ2-4 | minor | `manuscript/paper17/paper17.pdf` (last rebuilt in fada203, 2026-11-05) still prints "of Fisher–Rao type in the mass (Lemma 2.1(o)), not a Euclidean one" [snippet: pdftotext l.274]. `supplement.pdf` lacks the added ν proof line. The README says no pdflatex was available. | Rebuild both PDFs, or say in the README that they predate the correction. |
| PZ2-5 | minor | Part L's on-grid arm is run only where overlap ≤ 0.5 (11, 9, 7 kernels), exactly where the rule says D2 passes. D2 on-grid fails at 21 kernels (9 and 6 units at 1e-6 and 1e-12; passes at 1e-24; μ never 0.152) and at 16 kernels (9 units at 1e-6; passes at 1e-12; μ never 0.262); it passes at 13 [own: l_extra.py]. On-grid is also a different problem (no sliding needed). The md's "passes on-grid" thus reads as a spacing-independent control. | Run P0, D3, D2 on-grid at all six spacings (about 1 min) and print them. Phrase: "a control that removes the sliding; D2 passes iff the same margin rule holds". |
| PZ2-6 | minor | Two margins for one rule. Part D's agreement uses the later-entrant margin (122/126 over 6 ε); with μ never it is 124/126 [own: partd_rule.py]. Part L argues for μ never (14/18 against 18/18), and `margins` returns NaN for it when every kernel enters: `NaN > thr` is False, so the rule "predicts fail" while D2 passes (3 kernels, own cells). The 18 L cells are 9 dictionaries × 2 ε, which one reading of "18/18" overstates. | Use μ never in part D. Define NaN as +∞ in `report` and in the tests. State that the rule is validated on the 9 dictionaries of L and, out of sample, on 6 more (L_obs/L_rule 0.90-1.02). |
| PZ2-7 | minor | K2: "dormant units move 0.2–1.2 σ (P0: ≤ 0.02 σ)" compares medians over 11 targets [md K, raw]. P0's median is 0.0217 at 1e-6, not ≤ 0.02. Per-target maxima: P0 0.511σ at 1e-6; pairs 2.672σ (a0 = 0.1, 1e-6) and 2.568σ (a0 = 0.3, 1e-12); D3 0.89σ. | Print the maximum next to the median, or say "median". |
| PZ2-8 | cosmetic | The ordering D2 < pairs(0.3) < pairs(0.1) < D3 is of medians (0.11 < 0.38 < 0.97 < 4.65 at 1e-6; 0.20 < 1.16 < 4.54 < 12.93 at 1e-12); per target it fails once (target 2 at 1e-6: pairs 2.70 against D3 2.08). It is supported, and is an interpolation in a0 (a0-scan above), not a two-point coincidence. | Add the a0-scan row (W/scripts/a0_scan.py) and say "a0 is the interpolation parameter". |
| PZ2-9 | minor | M1: the tool's kill criterion is test_factorized_rbf.py alone; "+ the experiment's checks" is not part of it. P12 is behaviourally equivalent (killed by an attribute assert), N01 killed by a constant pin. 21 functions have no mutant, including PowerN (D4 and D5); my Q04 (PowerN.init s = ε) survives, so D4/D5 starts are pinned by the raw only. My E01 (the `n_atoms == K` clause) and E04 (per-target θ in the prediction) survive. Mutating the experiment script with the tool needs `--extra results/chaos/transit_certificate_raw.npz --extra results/chaos/parametrization.md`, otherwise even the control fails (the tool's docstring does not say so). | Docstring of the spec: "122 killed by test_factorized_rbf.py, 1 equivalent (S14); P12 equivalent in behaviour". Add tests for PowerN.init, the `n_atoms` clause (two units of one atom matched to two features), and a K job with θ ≠ 0.009 recomputed. Note the `--extra` requirement in `tools/mutants.py`. |
| PZ2-10 | cosmetic | The md header gives the Radau count for part C only (104); K (5 of 110) and L (15 of 72) are in the raw. HANDOFF.md calls K a long run (255 s) and not L (235 s). | Extend the header sentence to K and L; fix HANDOFF. |
| PZ2-11 | cosmetic | ν_i = β_i^{1/2} (paper17.tex l.66, l.127, Lemma (o); supplement) collides with the paper's ν for certificate margins (ν_0, ν_sk, ν_p, ν_V; l.225 "the margins of Section 6 are written ν"; Table `tab:notation` l.251-252). The new symbol is used in more than one place but not listed in the table. PARADIGM.md uses ϱ_i := √β_i. | Rename to ϱ_i or list ν_i in `tab:notation`. Optionally add "(constant diagonal weights, equivalent to a linear reparametrisation)" after "group learning rates". |

## Sources

- [md] results/chaos/parametrization.md (K l.312-347, L l.348-440, header, C, D, E, H, I, J).
- [raw] results/chaos/parametrization_raw.npz (K, L, B, D; per-cell fields `radau`, `aborted`, `one_each`, `dormant_wmax`).
- [snippet] src/core/factorized_rbf.py; experiments/exp_parametrization.py; tests/test_factorized_rbf.py and test_parametrization.py; tools/mutants.py and tools/mutants/wp1_factorized_own.py; manuscript/paper17/paper17.tex (l.66, 127, 160-215, 225-252, 272-277), supplement.tex (§pf:curve, §pf:branch), README.md; `git diff 5d75dc1 a00b405`; pdftotext of the committed PDFs; HANDOFF.md diff.
- [own] The referee's scripts in refPZ2/scripts: derive_c.py, myfr.py (independent engine), k_own.py, k_all_own.py, l_own.py, a0_scan.py, coupled.py (printed "pred" and "diff" columns mislabeled; superseded), coupled2.py (corrected), coupled_detail.py, partner.py, k_ext.py, k_ext_table.py, l_extra.py, rule_scan.py, partd_rule.py, e1_check.py, tex_check.py, bitident_kl.py, mut_audit.py, mut_cover.py, mut_which.py, s14_equiv.py, extra_spec.py, exp_mut_spec.py, exp_mut_spec2.py, ctrl_spec.py. Results in refPZ2/out (k_ext.jsonl, l_extra.jsonl), logs in refPZ2/logs, mutant results in refPZ2/mut.
- CLAIMS.md (section 5, read last) and SESSION_NOTES.md.
- [memory] Fisher-Rao on (0, ∞) is flat in √β; the central-difference error is O(h²); LSODA's stiffness switch can stall (documented in `rich_transit.relax_chunk`).

## Lesson

A law derived "with no fitted parameter" earns trust only where it can fail. Here the factor 1 − c/L is exact for the first entrant (14 of 14 targets, 3 digits), and the study's sample hid that it fails for the second atom of a same-signed coupled pair by a factor of up to 4, because "formed" is not "entered". The margin rule is the opposite case. It was tested where it could fail (6 dictionaries, 43 values of ε each, transition at 0.90-1.02 of the predicted L) and held. Both checks were cheap. State a derived law together with the premise it needs, then run the one case that breaks the premise. Separately, artefacts outside the commit (the PDFs, CLAIMS.md, SESSION_NOTES.md) go stale silently; each needs a regeneration step or a date.
