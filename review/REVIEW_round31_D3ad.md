# Referee report D3ad (round 31): Propositions D‴.39 and D‴.40, and the round-30 revisions of D‴.41 and D‴.42

Object: Proposition D‴.39 as in PARADIGM.md §5d (with RESULTS, [S67], `exp_gauge_plateau.py`, `results/chaos/gauge_plateau.md`); the draft
of Proposition D‴.40 (scratchpad `aniso/D40_draft.md`, `numerics_D40.md`, `pieces_D40.md`, `lit_S68.md`, the interim render of Parts A–D;
`src/core/gauge_noise.py`, `experiments/exp_gauge_noise.py`, `src/core/mass_transfer.py`, `tests/test_gauge_noise.py`); the round-30
revisions of D‴.41 and D‴.42 (`label/D41_r30_full.md`, `pieces_D41_r30.md`, `paper15_r30.md`, `multi/D42_r30_full.md`, `pieces_D42_r30.md`,
`results/chaos/gauge_label.md`, `gauge_classes.md`), at 7adc47e with Part E of `exp_gauge_noise.py` still running. The report as delivered
(the referee's scratch paths are under the session's scratchpad, `review31/`).

I made no changes to the repository. All work was done in the scratch copy …/scratchpad/review31/repo. Scripts are in …/review31/work/ and
outputs in …/review31/out/. Every run used OMP_NUM_THREADS=1 and nice 10, with at most 2 processes, each under 15 min, and no EXP_CACHE.

**One rule broken.** I used `pgrep -f` once, in a wait loop. It matched its own command line and spun; I stopped it with TaskStop. Nothing
was killed, and the loop only read the process list.

## What I checked and how

**Read**
- D‴.39: PARADIGM §5d (statement, proof, numerics), the table row, RESULTS, [S67]; `gauge_plateau.md`.
- D‴.40: the draft, `numerics_D40`, `pieces_D40`, `lit_S68`, `NOTES_D40` and the interim A–D render.
- Code and tests: `gauge_noise.py`, `exp_gauge_noise.py`, `mass_transfer.py`, `exp_gauge_label.py` (Part G) and the tests.
- D‴.41/42: `D41_r30_full`, `pieces_D41_r30`, `paper15_r30`, `D42_r30_full`, `pieces_D42_r30`, `gauge_label.md`, `gauge_classes.md`,
  `REVIEW_round30_D3ac.md`.

**Derived by hand**
- g′ = (1 + q)/(βa′).
- The generator identity, including the box and cohesion terms.
- The leading order and the Shahshahani-gradient reading.
- D = 0 from block structure.
- (iv)(c)'s spectral formula and Ψ(∞).
- (iv)(e): X_h − BX_hBᵀ = 2hS and the pair weight.
- τ_eff, and the barrier ODE.

**Part D replicated exactly**
- `work/d_diag.py` uses the same seed and draws. The first 8 batches match the stored U and Rhs exactly (max difference 0.0) for iso
  10⁻², iso 10⁻³ and corr 10⁻³.
- It adds a decomposition of E[u_k] into −τ_kE[g′_k], the cohesion term and the box term.
- It adds a chain-exact Euler identity, E[Ψ(α_{n+1}) − Ψ(α_n)] = 0 with Ψ = (a_max²/4βG_α)cosh(2α/a_max), which gives the Euler bias
  directly, plus the exact martingale standard error.
- It records plateau and box occupation, and autocorrelations.
- Extra runs, 8·10⁴ slow units each: amplitudes only at 10⁻² and 10⁻³; centres and widths only at 10⁻²; and a start on the plateau at
  τ = 2·10⁻⁴ (`work/p_diag.py`).

**Part E (still running)**
- I read the 32 finished `e_series` pickles in the coordinator's runcache, read-only (`work/peekE*.py`).
- I computed the effective noise levels of 'correlated 1', W♯_τ, and the barrier ODE (`work/misc2.py`).

**Part G**
- I checked G1a's zeros at γ = 10⁻⁵ for stability and basin (`work/g1check2.py`).
- I continued the corner zeros in τ_c and ran flows from W₀ with a 10× finer step (`work/g1check4.py`).
- Ψ(∞) was integrated numerically (`work/misc1.py`).

**Tests**
- The four requested files: 29/29 pass.
- All seven gauge files: 57/57 pass.

**Mutation testing**
- 26 mutants of the D‴.40 code (`work/mutate40.py`, `out/mutants.txt`); results under "Code and tests" below.

**From memory, no web search**
- Rugh, PRL 78:772 (1997), and Jepps–Ayton–Evans, PRE 62:4757 (2000): configurational temperature.
- Echeverría (1982): characterising invariant measures by E[Lf] = 0.
- Stein's method for diffusions: Barbour 1990; Gorham–Mackey.
- MYULA: Durmus–Moulines–Pereyra 2018; Brosse et al., COLT 2017.
- Ioannou's switching σ-modification.

## Verdicts

### D‴.39

**(i) Correct.**

**(ii) The Arrhenius slope has no error bar.**
- I refit the three escape rates of (C) by weighted least squares: slope (6.2 ± 1.1)·10⁻⁴, i.e. 0.81 ± 0.14 of ΔJ.
- So "82%" and "the prefactor still varies" are within 1.3σ of Freidlin–Wentzell's limit.

**(iii)(a) Correct.**

**(iii)(b) True for the deterministic objects, overstated for the stationary law.**
- The restored law spends 9–19% of its time on the plateau at τ = 2·10⁻⁴ (its own Part D), 37% at 10⁻³ and 23% at 10⁻² (my runs), at the
  wall.
- On the plateau (τ = 2·10⁻⁴) the box term in the masses' fitness is 1–3× the noise's selection: (3.0, 0.69, 1.7)·10⁻³ against τE[g′] of
  (1.0, 0.77, 0.90)·10⁻³.
- "Without bias" therefore does not hold for the stationary averages that D‴.40 uses.

**(iii)(c) Sketched, as the Status says.**
- V is only C^{1,1} (its Hessian jumps by σ_B at the faces).
- Katzenberger's limit and the entropic drift (third derivatives) are undefined in a boundary layer at P_B's faces.
- The vertex claim is not reached numerically: in (E) the masses are still 0.10–0.18 from the nearest vertex.

**(iii)(d) Fine**, given D‴.38.

**Numbers.** All numbers match `gauge_plateau.md`.

**Typos.** "(C) (C)", "(D) (D)", "(E) (E)", "(ii)..".

**Naming.** "Dead zone" is a Moreau–Yosida (switching-σ) penalty. MYULA is a missing precedent (memory).

### D‴.40

**(i) False as stated.** The claim is "exactly … in any stationary law".
- Without the box there is no stationary law (D‴.39(i)), so the claim is vacuous.
- With the box the exact identity is E[u_k] = −τ_kE[g′_k] + E[g_k∂_{α_k}(V + J_coh)].
- The proof's "the box makes α_k bounded" is false: the wall is a soft quadratic, and α_k was outside the box on 0.05–38% of steps. A hard
  wall would add a local-time term instead.

The measured residual E[u_k] + τ_kE[g′_k], against (i)'s claim that it is zero:

| run | z-scores (classes 1, 2, 3) |
|---|---|
| isotropic, τ = 10⁻² | 10.7, 8.6, 3.7 |
| isotropic, τ = 10⁻³ | 4.7, 4.2, 2.6 |
| correlated 1, τ = 10⁻³ | 2.8, 8.7 (class 3 not significant) |
| amplitudes only, τ = 10⁻² | 37, 5.1, 4.0 |

- At isotropic τ = 10⁻², E[u₂] = +0.0128 ± 0.0042 (8 batches; the md's 20-batch value is +0.0159 ± 0.0031). (i) requires a negative
  value, since g′ > 0.
- The cohesion term was exactly 0 in all six runs I decomposed. The whole correction is the box's.

**"Only each amplitude's own intensity enters; the correlations and the centres' and widths' noise do not enter".** This is exact only for
the prefactor τ_k. E[g′_k(α_k)] depends on all of Σ.
- Same τ_k = τ, amplitudes only against isotropic, at τ = 10⁻²: E[−τ₂g′₂] = −0.0230 against −0.0335, and E[u₂] = −0.0157 ± 0.0021
  against +0.0128 ± 0.0042 (6σ).
- At τ = 10⁻³ the E[g′] differences are 3–3.6σ, and E[u] is 2–4× larger in magnitude without centre/width noise.

**Leading order.**
- The dropped O(τ²) term, −τ_k(g″μ_α + ½g‴C_αα), depends on the full Σ.
- The relative correction is 2–4% at τ = 10⁻³ without the plateau, and 12–55% with it.
- Parts A and B compare two leading-order computations, so they cannot test "exactly".
- "A gradient whatever Σ" holds only at leading order (Shahshahani gradient of F_τ). `gauge_noise.py`'s own docstring says "in general not
  a gradient".

**(ii)**
- D = 0 is correct: block-diagonality suffices. "g is diagonal" is inaccurate, because the mass block in the simplex chart is not diagonal.
- "Zero drift exactly" fails in the boxed law. With centres and widths dithered only at τ = 10⁻², E[u₂] = +0.0244 ± 0.0050, all of it the
  box term; α₂ is outside the box 38% of the time.
- In the quasi-stationary well it holds: Part E's 'nonlinear only' stays at 1/3 ± 4·10⁻⁴.

**(iii)**
- Holds at leading order.
- The plug-in form's "exactly" inherits (i)'s failure and also needs the denominator floor to be inactive.
- τ̂ is biased by the drift's contribution, Δt·E[b_α²]/(2G_ατ). It reads 1.0005–1.013 with the drift, all 12 values high; the noise-only
  quadratic variation reads 0.996–1.003.

**Part D's statistics.**
- The gauge identity holds pointwise, so Part D's "difference" is exactly Lψ. Its χ² therefore tests only the generic E_π[Lψ] = 0
  (stationarity plus the Euler bias), not (i).
- The batch means are adequate for the difference: its τ_int is 36–50 units, and the batch-means standard errors agree with the exact
  martingale ones within ±30%. I found no 3300-unit correlation in it.
- The E[u_k] column's standard errors are underestimated. The plateau gives a slow tail (autocorrelation 0.07–0.09 at 10³ units), and the
  standard error of u₂ grows from 2.4·10⁻³ to 4.2·10⁻³ between 10³- and 10⁴-unit batches.
- The test's power is 7–22% of τ_kg′_k.
- The Euler bias is 0.9–1.4% of τ_kg′_k (chain-exact estimator), not "at most 4·10⁻³". The draft's bound is evaluated at the fit, but the
  law sits at the wall, where the curvature is σ_B = 1.

**Part E (unfinished, run without the box).**
- 'Correlated 1' heats the centres and widths 4.0–11.7× the dither's level (amplitudes only 1.46–2.26×).
- All 8 finished uncorrected seeds escaped by 2.1–5.3·10⁴ slow units:
  - 66–97% of records have |α| > 6, with |α| up to 39;
  - the fit error is 0.067–0.083;
  - W₃ fell to 0.065, below W_min,3 = 0.2.
- The κ_F = 2 class-levels variant therefore cannot test (iii).
- W♯_τ = (0.353, 0.294, 0.353) lies only 0.028 from W♯, against an expected standard error of about 0.01.
- 'Nonlinear only + barrier' matches its ODE: (0.3758, 0.2995, 0.3247) against (0.3762, 0.2992, 0.3246).

**Numbers in the draft.**
- "Up to 3·10⁻²" should be 5.4·10⁻² (correlated 1, τ = 10⁻², class 1), and the shift is the box's, not "the box and the cohesion".
- RESULTS and the table row say "10⁻⁶–10⁻⁷"; the measured range is 2.7·10⁻⁸–1.0·10⁻⁶.
- (C) quotes only the W = 1/3 rows. The other point's values are 3.2·10⁻⁴, 7.2·10⁻⁵, 7.9·10⁻⁶ and 8.0·10⁻⁷.
- All other numbers match the interim md.

**Literature.** (i) is the stationary Stein / adjoint identity with a test function of one coordinate: the per-coordinate configurational
temperature (Rugh; Jepps–Ayton–Evans) and generalised equipartition out of equilibrium (memory). The only specific ingredient is the gauge
identity. The "New" claims are overstated.

### D‴.41 and D‴.42 after round 30

**D3ac items**
- **D3ac-1:** fixed at γ = 10⁻² (statement and basins).
- **D3ac-2:** fixed, except paper15's title.
- **D3ac-3:** fixed (Richardson extrapolation).
- **D3ac-4 to D3ac-7:** fixed.
- **D3ac-8:** fixed, but the 1 + 0.8 sin t profile at W = 1/3 still has local exponent 1.41 between 10⁻³ and 10⁻⁴. "γ ≲ 10⁻³" should be
  10⁻⁴ there.
- **D3ac-9, D3ac-10, D3ac-12, D3ac-13:** fixed.
- **D3ac-11:** partly fixed. coord='a' and Heun with dither are now tested; the Part G code (g_flow, g_zero, g1m_job) is not.
- **D3ac-14:** mostly fixed; the componentwise ratios remain.

**(iv)(c)**
- The formula f(λ) = ∫ŝω²/(λ² + ω²) is correct in the continuum: checked numerically for OU, f = 1/(1 + λτ_c), and for Gaussian
  correlations, f ≈ (λτ_c)⁻².
- f is decreasing, Ψ is concave, and a smooth noise's Ψ is bounded.
- The proof's Ψ(∞) = ∫S|ω|/(4∫S) is wrong by 2π: it should be (π/2)⟨|ω|⟩. For τ_c = 1, 1.2533 numerically against 0.1995 from the
  draft's formula.

**(iv)(e)**
- X_h − BX_hBᵀ = 2hS, with pair weight 2/(λ_i + λ_l − hλ_iλ_l), is correct (algebra checked; G2 agrees to 7·10⁻⁴).

**Multistability against G1b**
- At γ = 10⁻² it matches, but "W♯ an unstable node" should read "an unstable node 1.6–3.3·10⁻³ from W♯".
- At γ = 10⁻⁵, G1b — "the whole zero set" — misses three stable zeros at τ_c = 10⁴ and 10⁵. I confirmed them by continuation and
  fine-step flows:

| τ_c | stable zero | stability exponents |
|---|---|---|
| 1.08·10⁴ | (0.25017, 0.54967, 0.20016) | −27.3, −27.0 |
| 10⁴ | (0.62456, 0.17526, 0.20018) | −41.5, −29.1 |
| 10⁴ | (0.25018, 0.17526, 0.57456) | −41.6, −29.3 |
| 10⁵ | (0.25003, 0.54994, 0.20003) | −7.1, −7.0 |
| 10⁵ | (0.62493, 0.17504, 0.20003) | −8.3, −7.2 |
| 10⁵ | (0.25003, 0.17504, 0.57493) | −8.3, −7.0 |

- With a 2·10⁻⁴ step, the flows from W₀ and its neighbour end at (0.62456, 0.17526, 0.20018), not at G1b's "boundary".
- The causes are g_flow's 2·10⁻³ step, which overshoots the guard, and g_zero's acceptance test, which is relative to the drift at its
  own start.
- So "where from τ_c = 10⁴ on every flow ends" is an artefact, and "at a finite replicator rate the selection is multistable instead" is
  wrong: it is equally multistable at γ = 10⁻⁵.

**Numbers.** All C, C2, D, E, F, G2, B and D‴.42 numbers match the md files, except the slips in D3ad-15 and D3ad-16.

### Code and tests
- **Tests:** all pass.
- **Mutants:** 26 run, 13 killed, 13 survived.
- **Survivors that matter:**
  - NL1–NL4: `run_dithered`'s `noise_L` path is completely untested. Even ignoring noise_L survives, and Part E's correlated variants use
    this path.
  - LV1–LV3: `levels` is untested, factor-2 error included.
  - FX1–FX2: the per-class levels in the one-particle-per-class branch, which Part E's class-levels variant uses, are untested.
  - GN7: `replicator_flow` is untested.
  - CL1: the experiment's `closed` is untested.
- **Survivors that are expected:**
  - GN8 and GN9 drop Σ's off-diagonal entries, or keep only the amplitude–amplitude block. They are equivalent at leading order, as the
    claim predicts.
  - MK3 is equivalent on binary masks.
- **Untested identities:** neither the exact generator identity nor the pointwise gauge identity has a test.

## Issues

| id | severity | issue | proposed fix (where) |
|---|---|---|---|
| D3ad-1 | major | D‴.40(i) "exactly … in any stationary law" omits E[g_k∂_{α_k}V]. Without the box there is no stationary law (vacuous); with it, (i) fails at z up to 10.7 (iso τ = 10⁻²) and 4.7 (10⁻³), with E[u₂] > 0. The proof's "the box makes α bounded" is false. "(i) proved" is unfounded. | D40 (i), proof and Status: state the identity with the box and cohesion terms. Present the box-free form as holding in the exact fit's quasi-stationary well up to exit-rate terms. Report the box term's size (Part D: separate columns). |
| D3ad-2 | major | "Only each amplitude's own intensity enters; centres'/widths' noise and correlations do not" is exact only for the prefactor. E[g′_k] depends on all of Σ: amplitudes only vs isotropic at the same τ_k gives E[u₂] −0.0157 vs +0.0128 (6σ, τ = 10⁻²) and 2–4× at 10⁻³. The O(τ²) term depends on Σ. | D40 title, (i), ROW0, BULLET, README, AUDIT: restrict to the leading order. Add an amplitudes-only vs isotropic row to Part D. |
| D3ad-3 | major | Part E runs without the box. 'Correlated 1' heats the centres and widths 4–12×; 8/8 uncorrected seeds escape by 5·10⁴ units, with \|α\| up to 39 and W₃ < W_min. The class-levels variant cannot test (iii). W♯_τ and W♯ differ by only 0.028. | exp_gauge_noise E: use P.BOX, a correlated Σ whose centre/width block is at the dither's level, and enough seeds to separate W♯_τ from W♯. Report the distance to both. |
| D3ad-4 | major | Part G1b ("the whole zero set") misses three stable zeros at γ = 10⁻⁵, τ_c ≥ 10⁴ (exponents −7 to −42). W₀'s flow reaches the class-1 corner. "Every flow ends [at the saturation]" and "multistable instead" (at finite γ only) are wrong. | g_flow: step capped by the distance to the guard. g_zero: accept relative to ref. Rerun G1m. D41 (iv)(e) last sentence and the (G) numerics; pieces (G). |
| D3ad-5 | minor | D‴.40(ii)'s "zero drift exactly" fails in the boxed law (centres/widths only, τ = 10⁻²: E[u₂] = +0.0244 ± 0.005). "g is diagonal" is inaccurate. | D40 (ii) and proof: qualify as in D3ad-1; "no cross block between centres/widths and amplitudes/masses". |
| D3ad-6 | minor | τ̂ is biased by the drift (+0.05–1.3%, systematic). The plug-in form's "exactly" needs the floor inactive and inherits D3ad-1. | D40 (iii): read the noise increments Δα − Δt·b (the loop knows b), or state the O(Δt) bias. |
| D3ad-7 | minor | Part D's χ² tests E[Lψ] = 0 only (the difference is Lψ pointwise). The Euler bias is 0.9–1.4%, not ≤ 0.4%. E[u_k]'s standard errors are underestimated (slow plateau tail). Power is 7–22%. | exp_gauge_noise D: report (i)'s residual and its z, the chain-exact Euler estimator and the martingale standard errors. Fix the numerics text. |
| D3ad-8 | minor | Numbers: "up to 3·10⁻²" should be 5.4·10⁻², from the box, not the cohesion. "10⁻⁶–10⁻⁷" should be 2.7·10⁻⁸–1.0·10⁻⁶. (C) quotes only W = 1/3. `gauge_noise.py`'s docstring contradicts (i). | numerics_D40, pieces RESULTS and ROW0, gauge_noise.py docstring. |
| D3ad-9 | minor | (i) is the stationary Stein / adjoint identity (configurational temperature, generalised equipartition; memory). The blue-moon analogy is weaker. The "New" is overstated. | lit_S68 [S68], AUDIT, Status "New": novelty is the gauge identity making u_k a one-coordinate observable. |
| D3ad-10 | minor | D‴.39(iii)(b) "without bias" does not hold for the stationary averages: the plateau holds 9–37% of the law at τ = 2·10⁻⁴–10⁻³, and the box term there is 1–3× the selection. | PARADIGM D‴.39 (iii)(b) and table row: limit the claim to the deterministic objects and Laplace order; quantify the plateau's share. |
| D3ad-11 | minor | D‴.39 (C): slope (6.2 ± 1.1)·10⁻⁴ = 0.81 ± 0.14 ΔJ. "82%" and "the prefactor still varies" are within 1.3σ of the Freidlin–Wentzell limit. | gauge_plateau C, D‴.39 numerics: give error bars. |
| D3ad-12 | minor | The C^{1,1} wall makes Katzenberger's limit and the entropic drift undefined at P_B's faces. "Dead zone" is a misnomer; MYULA precedent missing (memory). | D‴.39 (iii)(c), [S67]: use a C³ penalty or flag the layer; rename. |
| D3ad-13 | minor | Mutation survivors: noise_L (NL1–4), levels (LV1–3), per-class levels in the one-per-class branch (FX1–2), replicator_flow, closed. No test of the generator identity or the pointwise gauge identity. | tests/test_gauge_noise.py. |
| D3ad-14 | minor | D41 (iv)(c) proof: Ψ(∞) = ⟨\|ω\|⟩/4 should be (π/2)⟨\|ω\|⟩ (off by 2π). | D41_r30 proof (c). |
| D3ad-15 | minor | Slips: A2 "2–12%" should be 0.2–12%. "W₀ alone off by up to 2.2%" should be 2.7%. G1a's exponents are normalised by W₀'s near-zero drift (spikes −698, −1.86·10³). D42 (B)'s range covers n = 2 only. D42 (F) "zero to the resolution": values exceed the printed resolution by up to 2000×, with a ~1% floor at γ = 10⁻⁶. | D41/D42 numerics, pieces, gauge_label G1a units, gauge_classes F. |
| D3ad-16 | cosmetic | D‴.39 typos ("(C) (C)", "(D) (D)", "(E) (E)", "(ii).."). "1.0017 against 1.0016" is not in the md (1.002/1.002 there). paper15: PCPLACE placeholder, the title omits "Markov", "for white noise linear" sits next to the convex Damian–Ma–Lee Ψ. "W♯ an unstable node". "Weighted log-determinant". | PARADIGM, D41_r30, paper15_r30, D40. |

## Lesson

D‴.40's identity is exact, but only for a system that has no stationary law. In the boxed system, where the stationary law exists, the box
adds a term, and Part D's own numbers already show that term to be the dominant one. Part D's χ² test checks only the generic E[Lψ] = 0,
so it could not catch this. The same pattern runs through this round:
- "Without bias", "every flow ends" and "the whole zero set" were each read from objects that do not hold in the regime actually used:
  deterministic fixed points instead of stationary averages, a coarse flow step at a guard, and a leading-order trace formula standing in
  for an exact claim.
- Part E is testing (iii) in the regime D‴.39 says lacks a stationary law.

Before calling a statement exact, write the identity with every term the simulation contains and report each term's size. Make sure an
experiment's hypotheses actually hold in the runs that test it. Cross-check every continuation against a global search with steps
resolved to the geometry.
