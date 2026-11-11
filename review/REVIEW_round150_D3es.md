# Referee report D3es, round 150: the Paper 14 manuscript, second review

## Scope

I reviewed the Paper 14 manuscript at commit eaea29c (branch claude/vigilant-clarke-hf939b), as revised after round 143:
- `manuscript/paper14/paper14.tex` (28 pp., §1–9 to p. 25) and `supplement.tex` (15 pp.), with both PDFs;
- `README.md`, `make_figs.py` and `numcheck14.py`;
- the code and raw files behind the new runs: `src/core/gauge_batch.py`, `exp_gauge_fixman.py` Part J, `exp_gauge_plateau.py` Part G.

I read these against:
- `review/REVIEW_round143_D3el.md`;
- REVIEW_LOG's Round 143 entry and the two later Paper 14 entries ("runs declined after round 143"; "mutation tests … shortened");
- the README sections "Revision after round 143", "Runs after round 143" and "Shortening";
- PARADIGM.md's dated block "Corrections after referee round 143".

I worked on an export at `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref150/repo`, with the venv symlinked. The repository was not modified. Venue: SIAM J. Applied Dynamical Systems or Multiscale Model. Simul. Line numbers refer to `paper14.tex` unless marked S (`supplement.tex`).

## Status of the round-143 issues

| D3el id | resolved? | evidence | what remains |
|---|---|---|---|
| 1 (leak window) | yes | Derived st. 5.5 now assumes γ²τ/ε ≪ δ² and a stiffness small against "the modes that carry mass" (l. 396). The caveat is restored with numbers (l. 405, §8.1 l. 490, S2.2). S1.13 re-derives C = O(γ²τ/ε) [own: spread/δ = 0.166, 0.096, 0.053 recomputed]. Table S4 is relabelled "at the window's edge". Prop. 3.3(b) and Derived st. 5.5 now share the regime ε/τ fixed. | The condition "modes that carry mass" is not checked: the md gives the slowest mode's mass share but not the second mode's. At ε = 10⁻³ the margin is only 2.00·10⁻³/9.00·10⁻⁴ ≈ 2.2. The paper says it is "not quantified", which is honest. |
| 2 (chart) | yes | The remark after Lemma 4.2 and S1.6 give ℓ_ψ = 2Σ log(W_k/ψ′(A_k/W_k)). I re-derived it: ∂_ψ̃ f̂ = W_kφ_k/ψ′(a_k), and the tanh case has ψ′ = 1/a′. For ψ = id the flow is dW_k ∝ 2(KW_k − 1), bounded. §1.3 is qualified (l. 119). | Abstract wording (D3es-3): a physically white dither still drives the masses to the faces of 𝒲, with a bounded force. The Jacobian explains the divergence, not the direction. |
| 3 (floor) | yes | Derived st. 5.4 (l. 388) is re-derived. I checked the cap's fraction (1 − W_min²/W²)/0.01 ≈ 200δ/W_min, δ_cap = 0.0012595 and the threshold 0.002975. New runs in Part J. | The abstract overstates Part J's agreement (D3es-2). |
| 4 (Péclet) | yes | Pe = \|b\|δ/D = O(δ/γ) (l. 319). I recomputed 0.60 = 2.055·10⁻⁷ × 0.0833/2.85·10⁻⁸ and 6.02. §8.1 now compares ensembles. | "γ of about 10⁻³ or less" (l. 325, S1.7) is generous: at 10⁻³ the free spread at 1/(γτ) is 0.070 against the distance 0.083 (D3es-10). |
| 5 (obstruction) | yes | ll. 80, 193, 341; S1.8; Open problem (3) recast as an α-chart extension. | none |
| 6 (Arrhenius) | yes | Prefactor refits (S2.3, recomputed); Prop. 6.2 now has a bounded smooth domain; separation of scales stated; Part G added. | Wording "barrier confirmed to about 7%" (D3es-1). Prop. 6.2's start point and invariance hypothesis (D3es-4). |
| 7 (stationary claims) | yes | "at fixed masses" in the abstract, l. 424 and the §6 title. Prop. 6.3 is split ((b) deterministic, before the first exit); Remark 6.4 is formal and carries the box-force limitation; "vertex endpoint not observed" is stated. | none |
| 8 (literature) | yes | Shalova–Schlichting–Peletier, Anderson 1985, Egardt, Ioannou–Kokotović 1983/84, Bornemann–Schütte and Legoll–Lelièvre are added. The Fixman analogy's inverse direction is stated correctly (l. 123). | SSP's Corollary 5.7 (SGLD gives the log pseudo-determinant drift) is the precise precedent and is not cited (D3es-8). |
| 9 (labels) | yes | Prop. 5.1 is split ((c) proved, (d) conditional). Prop. 7.2(d) gives the exact K/m ratio, which I re-derived: the mass part is Σ mπ_iδW²/(γW) = n·(KδW²/(γW)), the position part n×, so D scales by n·n⁻² = 1/n. Remark 7.3 is reworded. | Prop. 5.2's condition is italicised but sits in its fourth sentence, not the first (cosmetic). |
| 10 | yes | The remainder is continued to γ = 3·10⁻⁶ (−2.545 → −2.646 per γ²τ). | none |
| 11 | yes | "3τ" replaced; "0.1δ*" withdrawn; "as τ → 0" removed. | none |
| 12 (notation) | partly | Now resolved: s versus π, the domain 𝒪, q in place of c in the proof of Lemma 4.2. | New or remaining clashes: B (box), B_k (class mass sum) and B_i (Brownian); p (density, coordinate count p, p = 3); c as the prefactor exponent in §8.4 and S2.3/S2.5; t_{α,k} undefined (S l. 61); T_A (S l. 213). See D3es-11. |
| 13 | yes | Tested: without `gauge_drift_raw.npz` the script raises FileNotFoundError and writes no md. | none |
| 14 | yes | Fig. 3(c) shows the masses form, Part E. | none |
| 15 | yes | S2.2 states that the seeds are shared. | Part J's seeds are shared across gains *and* τ (the same seed numbers 800–815, so the same normal draws), which is not stated (D3es-5). |
| 16 | yes | The averaged loop is derived (l. 149, [SB89]); MSC 93E15 and 60J70 added. | none |

## New findings

### Major

None.

### Minor

1. **Part G: "the barrier confirmed to about 7%" (l. 525, p. 24; Table 1, l. 106).**
   - 7% is the standard error at a *fixed* prefactor exponent c. The spread over 0 ≤ c ≤ 1 is 0.91–1.17, and c = 1 is 2.4 s.e. off.
   - With the barrier fixed at ΔJ_𝒪 the data do determine c [own]: c = 0.36, with 0.11–0.61 at Δχ² = 1 and −0.14–0.86 at Δχ² = 4. So τ⁻¹ is excluded at about 2.5σ, and the full prediction (ΔJ_𝒪, τ^{−1/2}) fits with χ² = 2.56 on 5 dof (p = 0.77).
   - So "prefactor undetermined" holds only when the barrier is free. The data are better than the text says.

2. **The abstract (l. 64) overstates Part J.**
   - It says the occupation "agrees with [δ^{κ_F−1}] within errors of about 0.1–0.2".
   - At κ_F = 2 the pooled exponent is 0.652 ± 0.217, 0.35 (1.6 s.e.) below κ_F − 1 = 1, and 0.17 below the layer law that the runs are compared with.
   - §8.3's "consistent, not sharp" is the right wording; the abstract should match it.

3. **The abstract's chart sentence (l. 64).** "It tends to −∞ at every saturation, through the Jacobian of the chart …, so the dither drives the masses toward the saturated fits" reads as if the drive were the Jacobian's. By (ii) of the remark, a physically white dither also drives the masses to the faces (ℓ_phys = 2Σ log W_k, flow ∝ 2(KW_k − 1)). What is chart-specific is the divergence, the layer that does not stop the drift, the exponent, and the vanishing gain at the endpoint.

4. **Prop. 6.2 (l. 427, p. 17).**
   - The exit time's start point is not stated; it should be the exact fit, or uniform on compacts whose flow stays in 𝒪.
   - The cited theorem (e.g. Dembo–Zeitouni (A-1)) needs trajectories starting in 𝒪 to stay in 𝒪; "closure flows to the exact fit" is weaker.
   - For this gradient system the conclusion still holds. The reversed gradient path from x* to the minimiser y* of J on ∂𝒪 cannot leave 𝒪 earlier, because J decreases along the flow and y* minimises J on ∂𝒪; any exit costs at least min_{∂𝒪} J − J*. One sentence in the proof fixes this.
   - The numerical basin check (S2.5) starts the flow from the minimiser plus 63 random boundary points, all at 480–4058τ. The low-energy boundary region, where exits happen, is untested apart from y* itself; the other L-BFGS local minima on ∂𝒪 would be the informative starts. The exit points are not recorded, so no check confirms that exits go through y*, i.e. by broadening χ₂.

5. **Part J's independence assumptions (l. 514; S2.5).**
   - The three classes share runs (this is stated). The κ_F and τ groups also share seed numbers 800–815, hence the same normal draws (common random numbers); this is not stated.
   - A delete-one-seed jackknife of the pooled estimator, removing the seed from every group at once [own], gives s.e. 0.133 (κ_F = 1) and 0.242 (κ_F = 2), against the naive 0.125 and 0.217.
   - The κ_F = 2 − κ_F = 1 difference, jackknifed on the paired seeds, is 0.697 ± 0.214 (z = 3.25, against the paper's independence-based 2.8).
   - Conclusions unchanged; the robust numbers should replace the naive ones.

6. **Part J's reference law (l. 514).** Katzenberger's 1-D zero-flux exponent at γ = 0.1 along W₂ = W₃ (the Part C machinery, run at γ = 0.1, seconds of CPU) [own, using `exp_gauge_fixman.drift`]:

   | κ_F | δ = 0.0063 | δ = 0.02 | δ = 0.05 | layer law at the same δ |
   |---|---|---|---|---|
   | 0 | −0.988 | −0.943 | −0.772 | −1.012 … −1.091 |
   | 1 | −0.027 | −0.058 | −0.035 | −0.025 … −0.167 |
   | 2 | +0.934 | +0.826 | +0.703 | +0.963 … +0.758 |

   - Inside the window the two references differ by at most about 0.13 for κ_F ≥ 1, below the errors, so the verdicts stand.
   - At γ = 0.1, however, δ*_k exceeds the distance from 1/3 to every face, so all three classes are in their layers at once. Quoting the limit's own exponent over the bins, or adding it as a second reference, would close the gap between the "layer law" and the limit equation.

7. **The cap is "a quarter of the finite-τ threshold" only for class 1 (l. 514, p. 23).** For classes 2 and 3 it is 0.12 and 0.20 of the threshold (0.00053 against 0.0043; 0.00060 against 0.0030).

8. **SIAM's abstract limit.** The abstract is 361 words (math counted as one word); SIADS and MMS allow 250.

9. **Exit-rate lag (S2.5).**
   - No exit occurs before about 2100 time units at any τ, comparable to the slowest normal relaxation time 1/3.06·10⁻⁴ ≈ 3300. So n/T slightly underestimates the quasi-stationary rate at the largest τ.
   - Discarding the first 1650, 3300 or 6600 time units shifts B/ΔJ_𝒪 by +0.02 to +0.04 [own], within the s.e.
   - Worth one sentence, or a shifted-exponential fit.

10. **Extrapolations without uncertainty, and γ ≈ 10⁻³ called deterministic.**
    - "About 7·10⁶" (l. 434, l. 523) carries the barrier's s.e. of ±0.14·ΔJ/τ ≈ ±1.1 in the exponent, a factor of about 3.
    - "0.77–0.86 times 1/(γτ)" (l. 525) spans only the three c values; the barrier's own ±0.41 in the exponent adds about a factor of 1.5 either way.
    - In l. 325, "γ about 10⁻³ or less" should read ≲ 10⁻⁴ for a sharp deterministic reading (Pe = 6 at 10⁻³ is marginal).

11. **Notation and cosmetics.**
    - The clashes listed under D3el-12; the c in "for c = 0" at l. 525 is used before τ^{−c} is introduced.
    - Fig. 3(a) puts the γ = 0.1 loop exponents on a panel whose curves are the γ = 10⁻² limit equation, with no mark for their own reference: add short ticks at the layer law's window exponents (−1.04, −0.08/−0.11, +0.82/+0.86). The legend's floor and no-floor entries are hard to tell apart at print size.
    - "1.54 s.e." (l. 488) is not derived visibly (README l. 265 has it).
    - A numbered remark in place of `remark*` would make "the remark after Lemma 4.2" citable.

## Verification I did

All marked [own]. My scripts are in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref150/own/`.

- **Regeneration from raw.**
  - `GF_OUT` → a directory holding `gauge_fixman_raw.npz` and `gauge_drift_raw.npz`, `--from-raw`: `gauge_fixman.md` is byte-identical (cmp).
  - `GP_OUT` with `gauge_plateau_raw.npz`: `gauge_plateau.md` is byte-identical.
  - Without `gauge_drift_raw.npz` the script stops with FileNotFoundError before writing.
- **Part J recomputed from `J_runs`** with my own window, fit and layer-law integration (`partJ.py`).
  - All 15 rows of Table S5 are reproduced to the printed digits; pooled −0.536, −0.045 ± 0.125, +0.652 ± 0.217; layer laws −1.045, −0.093, +0.823.
  - A multinomial power-law maximum-likelihood estimator agrees: pooled −0.534, −0.037, +0.711 ± 0.225.
  - Seed-level jackknife: see D3es-5.
  - From the raw (`partJ2.py`): no corrected record inside δ_cap or beyond a saturation; closest approach 7.18·10⁻⁴ (cap 5.3·10⁻⁴); κ_F = 0's last-quarter fraction beyond the saturation 0.1426.
- **Part G recomputed from `G_runs`** (`partG.py`, `partG_lag.py`).
  - Rates, exits (41, 46, 37, 38, 34, 22) and ΔJ_𝒪/τ = 2.52–5.73 reproduced.
  - The count-weighted fits reproduce 0.909/1.039/1.169 ± 0.065, with χ² = 2.38/2.20/2.13. An independent Poisson maximum-likelihood fit agrees to 0.001.
  - The profile in c with the barrier fixed and the lag sensitivity are in D3es-1 and D3es-9.
- **Batch integrator against `run_dithered`** (`batchcheck.py`). One heterogeneous batch of four runs: (γ, κ_F, τ, floor) = (0.1, 0, 10⁻⁴, 0.006), (0.1, 1, 10⁻⁴, 0.006), (0.1, 2, 5·10⁻⁵, 0.00424) and (0.01, 2, 10⁻⁴, 0.01), each for 8000 steps.
  - The largest difference was 2.7·10⁻¹⁵ in Θ and 7.8·10⁻¹⁶ in β.
  - A wall stop after one block followed by a resume from the checkpoint is bit-identical.
  - A raised cap continues the run identically (1.3·10⁻¹⁵).
- **Katzenberger 1-D exponents at γ = 0.1** (`katzJ.py`): D3es-6.
- **Tests:** `test_gauge_batch.py`, `test_gauge_fixman.py`, `test_gauge_plateau.py`: 39 passed (30 s).
- **Build:** three alternating pdflatex passes in scratch give 28 + 15 pages. The logs have 0 warnings and no undefined or multiply defined references, and there are no "??" in the text. Ten underfull boxes in the paper and four in the supplement, no overfull ones. Acknowledgements on p. 26.
- **`numcheck14.py`:** 974 tokens; 9 not in the md (1.54, 2.4, 2.6, 2.8; 0.631, 2.14·10⁻⁴, 2.25, 7.32·10⁻⁴, 8.39·10⁻⁴), each listed in the README as derived.
- **Shortening.** Sentence-level comparison of the pre-shortening paper (eaea29c^) and supplement against the new paper plus supplement (`cmpproofs.py`, `cmpall.py`):
  - every proof and derivation block survives verbatim, up to hard-coded references;
  - the other differences are condensations, Table 1's shortened rows, and the remark's restructuring.
  - Table 1's caption ("except Lemma 2.1 and Prop. 6.2") matches what the paper contains.
  - The main statements stay self-contained, with key-step pointers.
- **Hand re-derivations:**
  - ℓ_ψ in a general chart; the ψ = id flow 2(KW_k − 1);
  - the floor's exponent fraction, δ_cap and the non-normalisability;
  - the layer law p ∝ (W a′)^{κ_F}|dα/dδ| = δ^{κ_F−1}(δ + 2W_min)^{κ_F−1}W^{−κ_F};
  - the window spread γ(τ/ε)^{1/2} against 0.0044;
  - Pe;
  - Prop. 7.2(d)'s n·n⁻²;
  - the arithmetic in §8.1, §8.3 and §8.4 (0.95/1.36/0.90; 0.57/0.74; 2.8; 5.25; 1.3, 0.6, 2.4).

## Literature

- **Shalova–Schlichting–Peletier**, arXiv:2404.12293v1, read in full [paper].
  - The paper's description ("structure of the noise decides both the limiting dynamics and its time scale") is accurate.
  - Their Corollary 5.7 (SGLD, L̂ = L + ½wᵀη, degenerate case, time scale α²σ²) gives Brownian motion on Γ plus the drift −½P∇log|∇²L|₊ on the same time scale. That is the exact precedent for Prop. 4.1 and Derived st. 4.3 with isotropic additive noise, and should be cited beside LWA22 in the paragraph at l. 119.
  - Their §6.2 (correlated noise changes the regulariser) belongs with Open problem (1).
  - SSP do not discuss chart or reference-measure dependence (they note DropConnect's rescaling invariance), so the paper's chart point stands as its own.
- **Anderson 1985, Egardt 1979, Ioannou–Kokotović 1983/84:** placed correctly (drift without persistent excitation, bursting, the σ-modification's origin) [memory, plus the paper's search notes].
- **Bornemann–Schütte 1997 and Legoll–Lelièvre 2010:** apt for the stiff-limit free energy and for effective dynamics [memory].
- **The Fixman analogy's direction:** now right (l. 123). Fixman's correction makes the rigid system sample like the stiff one; this one removes the stiff free energy from the dithered stiff system, closer to adaptive biasing force (ABF) with a closed-form bias.
- **Suggested addition** for the τ^{−1/2} prefactor of an exit through a non-critical point of ∂𝒪 (l. 523, S2.3): Helffer–Nier (Mém. SMF 105, 2006) or Di Gesù–Lelièvre–Le Peutrec–Nectoux ("Sharp asymptotics of the first exit point density", 2019) [memory; to be checked].

## Overall assessment

**Recommendation: minor revision** (accept after minor changes), for SIADS or MMS.

All eight round-143 majors are resolved correctly. The changed mathematics checks:
- the leak window in its derived form;
- the chart identity and ℓ_ψ;
- the floor's non-confinement;
- Pe;
- the non-compactness obstruction;
- the bounded-domain exit;
- the "fixed masses" scoping;
- Prop. 6.3's split;
- Prop. 7.2(d).

The two new runs are well designed and honestly labelled:
- **Part J:** a τ^{1/2} floor whose cap stays below the threshold, γ = 0.1 to open the layer, and a jackknife over seeds.
- **Part G:** a smooth bounded domain, its barrier found by optimisation, and fits with three prefactors.
- **The batch integrator** reproduces the single runs to rounding, including resume.
- **The shortening** lost nothing.

What remains is wording and robustness:
- the barrier ("7%") understates what the data show about the prefactor;
- the abstract overstates Part J and conflates the chart's divergence with the direction of the drive, and is over SIAM's length;
- Prop. 6.2 needs its start point and one sentence on positive invariance;
- the pooled errors should be jackknifed at the seed level (the result is stronger: z = 3.25);
- notation clashes and figure references.

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3es-1 | minor | "Barrier confirmed to about 7%": that is the s.e. at fixed c; over c ∈ [0, 1] the fit is 0.91–1.17. With B = ΔJ_𝒪 fixed, c = 0.36 (0.11–0.61; c = 1 excluded at about 2.5σ); (ΔJ_𝒪, τ^{−1/2}) gives χ² = 2.56 on 5 dof [own partG.py] | §8.4 l. 525, Table 1 l. 106, S2.5: "1.04 ± 0.07 at the expected c = ½; with the barrier fixed, c = 0.36 ± 0.25, consistent with ½"; drop "prefactor undetermined", or keep it only for B free |
| D3es-2 | minor | Abstract: "agrees … within errors of about 0.1–0.2" at κ_F = 1, 2, but κ_F = 2 is 0.35 (1.6 s.e.) below δ¹ and is compared with the layer law, not δ^{κ_F−1} [md J] | Abstract l. 64: "consistent with the layer law of the limit equation (exponents −0.05 ± 0.13 and 0.65 ± 0.22 against −0.09 and 0.82), not a sharp test" |
| D3es-3 | minor | Abstract attributes the drive to saturation to the Jacobian; by the remark (ii) a physically white dither also drives to the faces, with a bounded force | Abstract l. 64: "the divergence (not the direction) is the Jacobian …" |
| D3es-4 | minor | Prop. 6.2: no start point; FW/Dembo–Zeitouni need 𝒪 positively invariant. The basin check starts only from the minimiser plus 63 random boundary points at 480–4058τ; exit points are not recorded | l. 427 and proof l. 431: start at the exact fit, plus the gradient-path argument (D3es text). S2.5: flow from the L-BFGS local minima on ∂𝒪; record and report the exit points |
| D3es-5 | minor | Part J pools fits as independent although classes, gains and τ share seeds (identical draws). Seed-level jackknife: pooled s.e. 0.133 and 0.242 (naive 0.125, 0.217); κ_F = 2 − 1 difference z = 3.25 (paper 2.8) [own] | l. 514, S2.5: state the common random numbers; replace the naive pooled s.e. and the 2.8 with the jackknifed values |
| D3es-6 | minor | Reference law at γ = 0.1, where all classes are inside their layers: Katzenberger's 1-D exponent differs from the layer law by ≤ 0.13 (κ_F ≥ 1) in the window, up to 0.3 at κ_F = 0 [own katzJ.py] | §8.3 or S2.5: add the limit equation's own window exponent as a second reference (Part C's code at γ = 0.1) |
| D3es-7 | minor | Cap "a quarter of the threshold": true for class 1 only (0.12, 0.20 for classes 2, 3) | l. 514: "for the first class, and below a quarter for the others" |
| D3es-8 | minor | SSP's Corollary 5.7 (SGLD gives the log pseudo-determinant drift on the Brownian motion's time scale) and §6.2 not cited; no reference for the τ^{−1/2} exit prefactor | l. 117/119; Open problem (1); l. 523 and S2.3 (Helffer–Nier; Di Gesù et al.) |
| D3es-9 | minor | Exit-rate lag: no exits before about 2100 time units, against a relaxation time of about 3300; a lag shifts B/ΔJ_𝒪 by +0.02–0.04 [own] | S2.5: one sentence, or a shifted-exponential fit |
| D3es-10 | minor | Extrapolated exit times quoted without uncertainty (factor about 3 for 7·10⁶); "γ about 10⁻³ or less" for the deterministic reading although the spread is 0.070 against 0.083 | l. 434, l. 523, l. 525, l. 325, S1.7 |
| D3es-11 | minor | Abstract 361 words (SIAM limit 250) | l. 64 |
| D3es-12 | cosmetic | Notation: B/B_k/B_i, p, c (prefactor), t_{α,k} undefined (S l. 61), T_A; c used before τ^{−c} is defined; Prop. 5.2's condition not first; Fig. 3(a) without the layer-law marks for the open markers; legend; 1.54 not derived visibly; unnumbered remark | Table 2, §7, §8.4, S1.2, S2.2, Fig. 3 / `make_figs.py` |

## Sources

- **[read]** `paper14.tex` (all 652 lines); `supplement.tex` (all 298); `review/REVIEW_round143_D3el.md` (all); REVIEW_LOG ll. 5321–5596 (Round 143, the runs, the mutants and the shortening); README sections (headings, ll. 355–438 in part, ll. 247–265); PARADIGM.md "Corrections after referee round 143" (l. 4928 on); `src/core/gauge_batch.py`, `mass_transfer.run_dithered`, `exp_gauge_fixman.py` (header, Part J), `exp_gauge_plateau.py` (Part G); the pre-shortening `paper14.tex` and `supplement.tex` (eaea29c^) for the comparison.
- **[md]/[raw]** `gauge_drift.md` Part I; `gauge_fixman.md` (regenerated); `gauge_plateau.md` (regenerated); `gauge_fixman_raw.npz` (J_runs, I); `gauge_plateau_raw.npz` (G_runs, G, G_dom).
- **[own]** Scripts `npzkeys.py`, `peek.py`, `partJ.py`, `partJ2.py`, `partG.py`, `partG_lag.py`, `katzJ.py` (output `katzJ.out`), `batchcheck.py`, `cmpproofs.py`, `cmpall.py` and `tests.out`, all in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref150/own/`.
- **[own] outputs elsewhere** in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref150/`:
  - regenerated md: `gf/gauge_fixman.md`, `gp/gauge_plateau.md`;
  - the build and Fig. 3 renders: `build/`;
  - the pre-shortening sources: `old/`.
- **[paper]** Shalova–Schlichting–Peletier, arXiv:2404.12293v1, read in full: abstract, §1.1–1.4, Table 1, §5.1.2 (Corollary 5.7), §6.2. Downloaded to `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref150/dl/ssp/` and read only through pdftotext.
- **[memory]** Dembo–Zeitouni's exit assumptions (A-1)–(A-4); Helffer–Nier, and Di Gesù–Lelièvre–Le Peutrec–Nectoux, for the non-critical boundary prefactor; Anderson, Egardt, Ioannou–Kokotović, Bornemann–Schütte, Legoll–Lelièvre.

## Lesson

A fit that leaves both the exponent and the prefactor free reports the weaker of two facts. Here, with the barrier free, the prefactor looked undetermined. With the barrier fixed at its independently computed value, the same six rates pick out the prefactor the theory expects and reject τ⁻¹.

Likewise, a pooled standard error that "treats the fits as independent" should be checked by a jackknife over the shared random inputs (here the seeds). It costs nothing, and here it strengthened the conclusion (z = 3.25 against 2.8) instead of weakening it.

Before calling a parameter undetermined, or an error conservative, run the fit with what theory fixes, and resample what the runs share.
