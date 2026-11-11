# Referee report R3: "Implicit Model Order Selection in RBF Networks" (full draft: main text + supplement)

Stance: senior, adversarial, fair. Nothing in the repository or the paper directory was modified. My scripts are in `/tmp/claude-0/-home-user/aa355708-d4bf-50e9-ad81-220b64c70ccb/scratchpad/referee_R3_work/scripts/` (`ci_check.py`, `kappa_check.py`, `recount*.py`, `lemma1_check.py`). No network was used. Files named `*repeat*` were ignored.

## 0. What I did

- Read the brief and CLAIMS.md (LIVE block). Read every `.tex` file of the main text and supplement, all 12 `tab_*.tex`, and `refs.bib`. Rendered all 9 pages of `main.pdf` and all 15 pages of `supp.pdf` and inspected the figures.
- Read the code that decides fairness: `exp_astro_m2.py` (CS_FILTERS), `exp_astro_baselines.py` (`run_crowdsource`), `imaging.py` (refit, `r_ex`, stops), `fig1_hero.py`.
- Numbers:
  - Matched all 309 "x [a, b]" strings in the `.tex` files against `results/**/*.md`. 283 match verbatim. The other 26 are re-expressions (sign-flipped losses, integer Wilson percentages), and I traced each one to its report line.
  - Checked all 80 point estimates of Supp. Table 6 against `m2.md` (0 mismatches).
  - Recounted the 1154-target statistics from `transit_certificate_raw.npz`.
  - Re-derived Lemma 1 and the core profile by hand and integrated them numerically.
- Not run: the SMLM or astronomy pipelines, the 48/200/720-target experiments. Byte-identical regeneration of reports is not verified.

## 1. The claims as I read them (before checking)

- **C1 Mechanism.** A Gaussian-kernel network whose coefficients are "mass × saturating amplitude", with all parameters (centres too) trained by gradient flow from a tiny scale ε, acts like a sparse-recovery path.
  - Units switch on one at a time, in the order of accumulated correlation with the residual.
  - A switched-on unit slides onto a feature; the rest stay at O(ε^{3/2}) with centres moving O(ε).
  - So "the number of units formed is the model order", with no second stage.
- **C2 Limit theorem (Thm 1).** Under (H1)–(H4) the training follows the "sliding inverse scale space" (ISS). The duals deviate at ε^{1/2} with an identified first-order term. For fixed centres this is an adaptation of Berthier and Pillaud-Vivien.
- **C3 Recovery.** For separated Gaussian features of kernel width the path ends with exactly K units (Thm 2, outline; Δ0 uncomputed and above every tested separation). The numerical record is 1154 + 720 targets. A "skeleton" predicts the outcome (Prop. 2).
- **C4 Noise.** To first order, with white noise of level s, there is a window of length ∝ 1/s in which stopping gives the oracle fit.
- **C5 Readouts.** The phenomenon needs effective depth ≥ 3. Depth 1 (the ordinary linear readout) does not select. Depth 2 only on a logarithmic scale. Adam does not.
- **C6 Algorithm.** Packaging the path gives a method "with no training data and no regularization parameter".
  - It beats tuned, offset-corrected Picasso on 8 simulated sets (15/16).
  - It scores F = 0.51 vs 0.25 (crowdsource) on 4 M2 patches.
  - It fails on M53/M92 until the kernel is fixed (post hoc).
  - It costs 10²–10⁴ × Picasso.
- **C7 Stated negatives.** The dual order gives no accuracy advantage over greedy. The FRC gain is the localization count. The CPU learned baseline is weak. No DECODE.

## 2. Verified (what I checked and agrees)

- **V1. Lemma 1 (supp_theory.tex:11–21).**
  - I re-derived (o)–(iv) by hand and integrated the flow numerically (DOP853, rtol 1e-12, arbitrary drive).
  - Q = β − v_s cosh(2α/a_max) is constant to 12 digits. η̇ = g holds (η and ∫g agree to 10 digits).
  - λ₁ = Λρ(1−εb/2v_s)^{1/2} exactly, and η(α) − (α − a_max coth)/(2Gv_s) = λ₁ to 1e-8 relative.
  - The core profile W = a_max(εb)^{3/2}cosθ_μ/(m(2v_s)^{1/2} sin³θ_μ) has ratio 0.99996, 0.9997, 0.9988 at μ = 0.5, 0.2, 0.1 (ε = 1e-6), i.e. 1−O(εμ⁻²).
  - So the mass–amplitude law, the wall, the ε^{-1/2} time scale, and ε^{3/2}μ^{-3} for dormant units are right. Eqs. (1)–(2) in sec_setting.tex agree with Lemma 1(o).
- **V2. Prop. 1.** I redid the argument.
  - F = (Z/2)(A²−2wAΓ+w²) has a unique critical point with w≠0, namely (A,c*).
  - The Lemma "mono" computation q̇ = Λ0GmZ²A²w|d|²Γ²/(4σ⁴) − Zẇ is correct, and the contradiction argument at a first zero of q is sound.
  - Other duals end at Γ_j/Γ_e < 1.
- **V3. Thm 2 arithmetic and architecture.**
  - T_max ≤ 1/(Z|A|_min Γ(h√n/2)) and the choice of r0 give far kernels ≤ μ0/2.
  - The first-entry error is Γ(Δ−h√n/2), not Γ(Δ).
  - For K=2, R=2.4, μ0=1−Γ(h)=0.159 (h/σ=0.833) I get r0=4.07 and 2r0σ=8.1σ, which is the "≥ 8σ" of Remark (a).
  - The induction structure is sound. The gap is Lemma 4, not the architecture.
- **V4. Prop. 3 (supp_theory2.tex:19–21).**
  - I re-derived the entry-time recursion (i) and Σ⊥=⟨φ_j,R_Kφ_l⟩/N, and the law of sW.
  - The discrepancy-rule "late" probability is P(χ²₅₅>1.21·61)=0.046, matching the report.
  - The oracle error is 2Ks²/N for 2 parameters per unit.
  - Numerics (`noisy_certificate.md`): W_obs/W_pred = 1.0001 [0.983, 1.030] at s=0.003. W·s has median 13.5, 14.3, 13.0, 12.7 for s ≤ 0.1, so W ∝ 1/s.
  - It degrades at s=0.3 (order swapped on 81/120 draws; discrepancy rule early on 89/208).
- **V5. Numbers.**
  - Recount from raw: 1154 targets, 763 "nearest", 391 "other", 0 non-exact, max rms 3.1e-14. Linear skeleton sign right on 1149, fitted skeleton on 1154, amplitude-free certificate on 903. 1034 one-dimensional + 120 planar.
  - 720-target counts, 379 gated, 0.99 (`separated_regime.md`). 12/48 (H1) targets, slopes 1.501–1.504 / 1.001–1.002 / 1.233–1.533 (`transit_limit.md`). 98/200 and 187–198/200 (`parametrization.md` C).
  - Table 9 percentages, "16 of 16" network-follows-path, 1.3e-9 / 3.7e-8. FRC numbers. CPU times. Learned-baseline tables.
  - Summary counts reproduce: E 15/1/0, J 13/3/0, RMSE 8/5/3 from `iss.md` + `smlm2013.md` records; astronomy 6/2/2 and 8/1/1 from `m2.md`.
- **V6. Cross-references.** Supplement Tables 1, 2, 6, 7, 8, 9 and Fig. 1, the Theorem/Proposition numbers and section numbers cited in the main text match the compiled PDFs.
- **V7. Protocol.** The M2 frozen test has mechanical hashes (`freeze.json`). Dev patches (110,630), (210,630) are disjoint from the test patches, though adjacent windows of one frame. κ=2 at the grid edge is disclosed. Truth counts 923–2057 are right.
- **V8. Unflattering results are reported.**
  - Dual = greedy in SMLM (ΔE 0.0 [−0.0, 0.0]; 0.1 [−2.5, 2.1]).
  - Astronomy greedy−dual −0.009 [−0.029, 0.012].
  - FRC gain is count (53.4 vs 55.1 nm at equal count).
  - The M92 loss is reported.
  - The 18.6-point fragility to a wrong noise law is reported.

I found no wrong formula in what I could check. The mathematical problems below are gaps and statement precision, not errors.

## 3. Findings

### A. Claims versus what the paper shows

**F1 [blocker] "Answers both on its own", "no training data and no regularization parameter" is not what the paper demonstrates.**
- Location: `abstract.tex:1`; `sec_intro.tex:4` ("no penalty, threshold or training set entered"), `:11`; `sec_conclusion.tex:3`.
- Evidence:
  - With noise the network does not stop. The paper itself says "stopping needs a rule" (`sec_theory.tex:30`): discrepancy needs s; significance needs s and κ; largest gap is "observed, not proved". In iterative regularization the stopping time is the regularization parameter, and Prop. 3's window ∝ 1/s says exactly that.
  - The parameter-free options lose. The largest-gap rule is −35.9 [−45.0, −26.7] and −12.9 [−16.5, −9.9] E-points against the significance rule (`sec_exp.tex:11`). The "nothing to set" label-free BIC row of Table 1 has E = 52.9/37.2/35.1/23.4 against tuned Picasso 63.3/36.7/46.9/24.6, i.e. −10.4, +0.5, −11.8, −1.2. It beats only Picasso's defaults.
  - The "tuned" rows behind the 15/16 use 168 operating points chosen on up to 1000 labelled frames. They also use σ0 = 1.4 px measured on labelled development frames (0.68/0.997 px for 2013), a position offset estimated on labels (`supp_data.tex:5–7,15`), and a noise law calibrated on unlabelled development movies. This is development data in the sense of `sec_intro.tex:2`, which describes the competitors' second stage as "a penalty weight, a detection threshold, a stopping criterion, tuned on development data".
  - `sec_alg.tex:18` ("Inputs: a PSF model, a noise model, and the stop") and `sec_intro.tex:11` ("only a significance level") are honest. The abstract drops that qualifier.
- Fix:
  - Say what is implicit: the order of entry and the support at any stopping time. Present κ (or BIC) as the single stopping level, with its default.
  - Delete "no training data" or write "no model is trained; PSF and noise law are calibrated from unlabelled frames".
  - Scope "no threshold" in `sec_intro.tex:4` to the noiseless illustration.
  - Put the BIC-row verdict in the introduction.

**F2 [blocker] The thesis (implicit selection in RBF networks) and the real-data evidence concern different objects. Nothing connects them on imaging data.**
- Location: `sec_alg.tex:2,16–17`; `sec_exp.tex:5,11,19`.
- Evidence:
  - "Since the limit of training is the path, we compute the path." But the imaging algorithm departs from Def. 1/Thm 1 in at least these ways:
    - high-pass of frame and kernels (F 0.653→0.260 without it);
    - noise-weighted correlation;
    - a pixel-integrated PSF with learned per-unit width and halo;
    - exclusion radius `r_ex` (−0.050 [−0.083, −0.021]);
    - windowed refits by sign-constrained least squares (`imaging.py:239`, `least_squares(..., method='trf')`, a local optimiser, not the entry-layer flow of (H3));
    - a time-stop with clean-up, or BIC.
  - Agreement of the refit with the inner flow was checked on 48 one-dimensional targets only (3.1e-11, supp 1.4).
  - No factorized network is ever trained on imaging data. The network-vs-path evidence is 16 noisy 1-D trials (ε=1e-8), 200 noiseless 1-D targets, and 8 targets for the rates.
  - The authors' own ablations show the gains do not come from the accumulated-dual order:
    - SMLM greedy = dual (N1.LD E 68.2 and 68.2; N1.HD ΔE 0.1 [−2.5, 2.1]).
    - Astronomy −0.009 [−0.029, 0.012].
    - Synthetic: MP+slide 84% vs ISS 82%.
    - Close pairs: the dual is net worse, −1.3 [−2.1, −0.7] (BIC) and −7.7 [−8.4, −7.1] (κ=5).
  - Under the same (SMLM, instantaneous-significance) stop the astronomy greedy path scores F=0.644 and the dual path 0.412 (`m2.md`, "Design variants"). The dual needs its own time-stop and clean-up to reach 0.653. The "within 0.009" in `sec_exp.tex:19` compares across stop rules.
  - Sec. 5.1: the ISS limit path is the slowest accurate method (22.8 s/trial against 0.46 for MP+slide, 0.69 for SFW, 9.2 for CPGD; `noise_order.md` "Cost").
  - A reader who accepts every number still learns nothing about RBF networks choosing their order. The imaging results would read the same under the title "A calibrated sliding greedy localizer".
- Fix:
  - (a) One experiment: train the factorized network by gradient flow on the same weighted high-passed loss on a 40×40 px simulated SMLM or astronomy patch, and show it forms the algorithm's units.
  - (b) In any case, re-frame §4–5 as "an application of the path", list the departures, and state that the imaging results do not test the implicit-selection mechanism.

**F3 [major] Status and precision of the theory statements.**
- Location: `sec_intro.tex:8`; `sec_theory.tex:3,8,10,18–20,23`; `supp_theory.tex:2,28,43,47`; `supp_theory2.tex:11–13`; `supp_sep.tex:17–20,44–48`.
- Evidence:
  - (i) "Proved under stated hypotheses" (`sec_intro.tex:8`) vs "proof outlines" (abstract) vs "Proof outline (five steps)" (`supp_theory.tex:43`).
    - Step 3's near-field estimate "forces ε0 exponentially small" (Remark 1; |s_θ+σ_k| spans 2.7–9.9·10⁴ over 48 targets).
    - Lemma 4 (coupled entry) is a sketch with δ0 "not computed". (H3) is "assumed otherwise" (`sec_theory.tex:10`).
    - `sec_theory.tex:10` says the fixed-centre convergence "is proved", but `supp_theory.tex:28` says "Proof sketch" and "The rate ε^{1/2} is observed, not proved".
  - (ii) Thm 2.
    - On the proof's constants Δ0 ≥ 8σ (≈12–14σ with δ_Δ). Tested separations are ≤ 8σ (K=2) and 5.8σ (K=3). Remark (a): "The theorem therefore covers none of the tested targets".
    - Its last sentence chains two outlines into a statement about networks. The internal status was "Conjecture" (CLAIMS D1).
    - Cheap test not run: extend Table S1 to Δ = 10, 12, 14σ in a wider window (1-D, seconds per target).
  - (iii) Prop. 2 (`sec_theory.tex:23`) is a list of counts. The supplement's "Proposition (skeletons' accuracy)" has neither proof nor sketch, although `supp_theory.tex:2` promises each result is "proved or sketched".
  - (iv) Def. 1 is under-determined: "the active block … relaxes to a stationary point of F", but F is non-convex in the centres. The object is fixed only by the entry-layer flow of supp (H3), which is assumed to converge. The algorithm's refit is a different relaxation.
  - (v) Thm 1 as stated.
    - "(H1–H4)" glosses only three hypotheses. H4 (regular crossings) is never described. μ is undefined.
    - (b) "between entries" is true only on [τ_k+Cε^{1/2}log(1/ε), τ_{k+1}−ε^{1/4}] (`supp_theory.tex:41`).
    - (c) "enters" must mean "dual reaches its wall". The coefficient crossing time carries an extra s^θ_k.
    - "Converges at rate ε^{1/2}" holds for duals and entry times. In L² of the output the fixed-centre report shows slopes 0.260/0.696/0.271 ≈ ε^{1/4} (`mirror_transit.md` l.23,39,55).
    - "Observed remainder O(ε)" (`sec_theory.tex:10`): the dual-remainder slopes are 0.391–0.501 in ε^{1/2} units (ε^{0.89–1.0}); the block distance is ε^{1.23–1.53}.
- Fix: relabel ("Theorem 1 (proof outline)", "Claim 2", Prop. 2 → Observation); add the missing definitions of H4 and μ; name the norm; align intro and abstract.

**F4 [major] "A closed-form skeleton predicts which units carry them" (`abstract.tex:1`); "we can say which units and when" (`sec_conclusion.tex:3`).**
- Evidence:
  - Prop. 2 predicts one binary event: "the path is exact with the nearest kernels". The linear skeleton is right on 1149/1154 (always-nearest 763, certificate 903). The fitted skeleton's 1154 is by construction (`supp_theory2.tex:14` says so).
  - The raw records store skeleton margins (`m_skel`, `m_fit`) but no skeleton entry sequence. Which kernels the 391 "other" cases use was not scored.
  - The skeleton takes the target's true centres and amplitudes as input, and is a recursion with projections, not closed form.
- Fix: "predicts, from the target's features, whether the path uses the nearest kernels (1149/1154; always-nearest 763)". Delete "which units and when", or score the entry set and times on the 391.

**F5 [major] "15 of 16 comparisons (one tie)" (abstract; `sec_intro.tex:11`; `sec_exp.tex:9,24`).**
- Evidence:
  - 16 = 8 sets × {E-tuned point, J-tuned point}, with E evaluated at both. For the 2016 sets E at the J-tuned point appears in no table. Table 1 shows E at the E-optimum and J at the J-optimum.
  - "4 wins at the J-optimum" (`sec_exp.tex:9`) follows a sentence about J and reads as J, where Table 1 shows 2 wins and 2 ties.
  - The 16 that the paper plots (Supp. Fig. 1: E at E-opt, J at J-opt) read 13 wins / 3 ties / 0 losses (E 7/1, J 6/2).
  - The sets are not 8 independent benchmarks: the four 2016 sets are one structure at two SNRs × two densities; the two 2013 pairs share a ground truth (`supp_data.tex:15`).
  - On Bundled HD both methods have E = 5–21 (RMSE-dominated; ΔJ = +57), so its E "win" is cosmetic.
  - RMSE is 8/5/3, with three losses (up to +39 nm), all on the 2013 bundles.
- Fix: "E at the E-tuned point 7 wins/1 tie over 8 sets, J at the J-tuned point 6/2, none lost; four conditions of one simulated structure plus four 2013 sets". Drop "of 16".

**F6 [major] Transfer failure and post hoc fix: honest in §5.4 and Table S8, not in the abstract or intro. Purity is hidden.**
- Evidence:
  - `abstract.tex:1` says "fails until its PSF kernel matches the PSF stamp"; `sec_intro.tex:11` says "does not transfer … until it does". "Until" asserts a cure that exists only as a post hoc change chosen after seeing the same six patches (`psf.md`: "POST HOC, NOT A HELD-OUT TEST"). The held-out test is a `\tbd`.
  - Table 2 gives only F for M53/M92. Purity is 0.272 (M53, tuned) against 0.416 (crowdsource) and 0.558 (photutils). On M92 it is 0.116 against 0.433/0.597.
  - With the adequate kernel M92 purity is 0.253 against 0.433 (−0.180 [−0.235, −0.129]), i.e. 75% of detections are false. The lead is completeness bought with purity.
  - `sec_exp.tex:20` says "(at a lower purity; Supplement Table 8)", but Table S8 has no purity column.
- Fix: add "(post hoc)" and the purity to the abstract clause, e.g. "a post hoc kernel change gives F 0.277 vs 0.226 on M92 at purity 0.25 vs 0.43". Add completeness and purity columns for M53/M92 to Table 2 or S8.

**F7 [major] Astronomy: the baselines are less "tuned" than the text says, and the ISS's largest effect is available to no baseline.**
- Location: `sec_exp.tex:2,18`; `supp_data.tex:24`; `exp_astro_m2.py:64`; `exp_astro_baselines.py:158`.
- Evidence:
  - The 48 "crowdsource variants" are post-filters (min flux/dflux, fracflux, qf, spread bound) on one default `fit_im` run. They can only delete detections. The development search returned the unfiltered default: pooled counts are identical (911 matched / 1866 detected) and Table 2 rows 3–4 are identical.
  - "Crowdsource, tuned" is therefore "crowdsource, default". No pipeline parameter (detection threshold, blend threshold, iterations, sky model) was searched, while the ISS's tuned parameters (σ_b, r_ex, κ) are pipeline parameters. "At least the same tuning effort" is not met in kind.
  - The ISS's dominant ingredient is the 3-px border-consistent high-pass: without it F falls from 0.653 to 0.260, below crowdsource's 0.377 on that patch. No baseline was run on a high-passed patch with a high-passed PSF stamp.
  - κ=2 is the grid edge (disclosed). Table S6 shows F still rising with κ up to 8 on M92. PCAT was not run (disclosed).
- Fix: run crowdsource and photutils on the high-passed patch and stamp, and add a pipeline-parameter sweep for crowdsource. Otherwise write "default" for crowdsource and state that the search space was post-filters.

**F8 [major] Contribution 3 ("Which readouts select") is undocumented and its scope exceeds its evidence.**
- Location: `sec_intro.tex:10`; `sec_theory.tex:33`; `main.tex:23`.
- Evidence:
  - The supplement has no section on it. Searching it for depth, Adam and fixed-step hits only the U-Net paragraph. The paper's documentation is Fig. 2 and one paragraph.
  - "Plain gradient flow or fixed-step descent" rests on `parametrization.md` Part F. That is D3 at ε=1e-6 on the design and 4 targets (6 runs), and P0 on the design (1 run). P0 at ε=1e-4 is not exact on any of its 5 runs.
  - "Adam forms all units" rests on one design (Part H). P0 at step 1e-6 forms nothing (rms 0.42).
  - The depth-2 verdict uses 21 family-A targets ("the easy family") plus the 200-target set, which is "a stratified sample, harder than the population" (CLAIMS §5; not said in Fig. 2). Everything is 1-D with 21 kernels (12 planar targets).
  - "Wherever the margin exceeds a few tens of ε^{1/2}" describes a post hoc stratification at 0.02. At ε=1e-8, 0.02 is 200 ε^{1/2}, and mismatches persist to ≈57 ε^{1/2}.
  - The caption defines "formed" as |w|>0.01 max|A*|. The counts use min(0.01 max|A*|, 0.1 min|A*|) (`parametrization.md` l.3), which matters for the amplitude ratios up to 749 in family B. θ, "margin" and the Radau-fallback cells are never defined or mentioned.
- Fix: a one-page supplement section (design, threshold, margin law, Parts B/C/F/H as tables). Scope the sentence to "gradient flow (and, on 7 runs, fixed-step descent)". Delete "Adam forms all units" or show it.

**F9 [major] Title and first-page scope.**
- Evidence:
  - Fig. 2(a) shows the ordinary linear readout (depth 1) forming 12 of 21 units. The mechanism needs a factorized readout and α(0)=0, ν(0)²=εb exactly (`sec_setting.tex:10`).
  - With small random initial amplitudes the limit is a different, logarithmic-scale object. K units still form (9–10/11 exact) but not the path's units (`parametrization.md` Part G).
  - The Fig. 1 caption (`main.tex:14`) says "plain gradient flow … three features" and omits "factorized readout". The target is hand-picked and off-grid (`fig1_hero.py`).
  - A reader of page 1 concludes that ordinary RBF networks select their order.
- Fix: title "… in RBF Networks with a Factorized Readout". Fig. 1 caption: "factorized readout; hand-chosen target". One sentence on initialisation in the abstract.

**F10 [major] The SMLM evidence base.**
- Evidence:
  - One baseline: Picasso, a one-peak-per-blob localizer. The report's own recall ceiling for it on HD sets is 59–65% (`picasso_baseline.md` l.606, 1100). There is no multi-emitter fitter (DAOSTORM, ThunderSTORM multi-emitter, FALCON), which is the like-for-like comparison for a calibrated sliding fitter.
  - "Widely used tools tuned with far more effort" (`sec_conclusion.tex:3`) means Picasso. On the 2013 sets its defaults give E = −3 to −13 (Table S4).
  - All 8 sets are simulations (disclosed). The "a priori default" is "fixed in the code before its test run" but encodes what was learned on these simulations' labelled development frames (halo 0.25 at 3σ, σ0 = 1.4 from labels, `r_ex` = 2). On 2013 the Gaussian PSF lies inside the ISS kernels, and the noise slope 1.41 is "taken from the page".
  - The only real data is FRC with no truth: the gain "is mostly the count", repeat-localization precision is 7.1 vs 6.1 nm, χ²=3.6, and 26% of units sit at the width bound.
  - Learned baseline: 289k-parameter U-Net, 3–12 CPU-minutes, 997 labelled emitters at LD and 80 frames at HD. E rises 53.7→60.4 from 100 to 800 frames. It is uninformative about DECODE-class methods (acknowledged).
- Fix: add one multi-emitter baseline, or restrict the conclusion to "Picasso". Rename "a priori" (see F27) and say "four conditions of one structure".

**F11 [major] Positioning against prior art is not on page 1, and close prior art is missing.**
- Evidence:
  - `sec_intro.tex` cites only Chen 1991 and Platt 1991. The distinction from the closest work [4] appears at `sec_theory.tex:5` (p. 2) and in §6 (p. 7).
  - Missing, to my recollection (offline, please verify):
    - Osher et al. 2016 (`osher16` is in `refs.bib` but never cited);
    - Vaškevičius–Kanade–Rebeschini 2019 and Li–Ma–Zhang 2021 (Hadamard-type over-parametrisation with early stopping, the precedent for Prop. 3 and the depth contrast);
    - Poon–Peyré 2021;
    - Maennel et al. 2018 and Lyu et al. 2021 for early alignment in two-layer nets;
    - DAOSTORM, FALCON and 3B for SMLM.
  - `sec_related.tex:3` claims "on separated features all succeed and the difference is the tuning each needs". In Table 9 every method needs a stop. ISS+BIC needs exactly what MP+slide+BIC needs (82 vs 84%) and is the slowest.
- Fix: a 3-sentence positioning paragraph before the contributions: closest prior [4] (fixed centres); new = moving centres, rate and first-order terms, noise window, depth contrast for moving features. Cite the missing works. Delete or support the "tuning" sentence.

### B. Exposition

**F12 [major] The first page and the density of §3 and §5.**
- Evidence:
  - `abstract.tex:1`: 251 words, one paragraph, 12 numbers. The clause "(a perturbative argument covers generic features beyond a separation that exceeds our benchmarks'; observed on 1154 benchmark targets down to 2σ)" is unreadable. "Skeleton", "s", "oracle fit", "effective depth", "Picasso", "offset-corrected" are undefined for an IJCAI reader.
  - `sec_intro.tex:11` is a 140-word sentence with 8 numbers. `:13` ("We do not claim …") is defensive. "Reaches a wall" (`:4`) is undefined on p. 1.
  - Undefined at first use: ISS (used from `sec_exp.tex:5`, never defined), margin, R_K (`sec_theory.tex:28`), θ (`:33`), Γ and δ_Δ (`:18–20`), G, γ, m, a_max (`sec_setting.tex`), H4.
  - `sec_theory.tex:20` is one 330-word paragraph mixing proof idea, Δ0 constants and eight numerical results. The Thm 2 statement is 190 words.
  - `sec_exp.tex:9,11,12,19,20` are 330–450-word paragraphs with 20–30 numbers duplicating the tables; `:12` packs five topics.
  - Related Work is §6. Legends in Figs. 1–2 are about 5 pt at print size.
  - Fig. 2(b)'s x-axis is cumulative ("margin at least μ0"), which the caption does not make obvious. Fig. 1's "stay at O(ε^{3/2})" is contradicted by dormant units near ε=1e-8 (they are O(ε^{3/2}μ⁻³)).
- Fix: move numeric prose into tables, split `sec_theory.tex:20` into statement / idea / test, define the symbols, shorten the abstract. Suggested abstract (about 160 words):
  > A Gaussian-kernel network whose coefficients are a mass times a saturating amplitude, trained by gradient flow from a small scale ε, adds units one at a time: a unit switches on when the accumulated correlation of its kernel with the residual reaches a wall, slides onto a feature, and the others stay dormant. Under checkable hypotheses (proof outline in the supplement) the training converges, as ε→0, to an explicit path with first-order corrections; for separated features the path has one unit per feature, and under white noise of level s, stopping within a window of length ∝1/s gives the oracle fit at first order. The effect needs readouts of effective depth ≥3. Computing the path directly yields a localization/photometry algorithm with one dimensionless stopping level: it beats tuned Picasso on simulated microscopy and tuned crowdsource on M2; its failure on clusters with worse seeing is diagnosed, and a post hoc kernel change removes part of it.

### C. Consistency and hygiene (minor unless stated)

| # | Location | Problem and evidence | Fix |
|---|---|---|---|
| F13 | `abstract.tex:1`; `sec_intro.tex:9`; `sec_theory.tex:20`; `supp_theory2.tex:14` | "down to / from 2σ". The 1154 contain 40 targets at 1.5σ (family D(b), s=1.5 rows; raw recount: min 1-D separation 1.5σ). | "1.5σ" |
| F14 | `sec_theory.tex:20` | "amplitude ratios up to 2.4" is the 720-target experiment's (`supp_sep.tex:45`) and families A, C. 329 of the 1154 have ratio >2.5 (B up to 749, D(a) up to 10³). | "up to 10³" |
| F15 | `sec_exp.tex:9` | "4 wins at the J-optimum" is ambiguous (see F5). | State metric |
| F16 | `tab_smlm_rmse.tex:11` | Caption: "lower or equal … at the E-optimum on three of four". The report records 1 win, 3 ties, 0 losses (not worse on 4/4; ISS point estimate higher on 3). | Correct |
| F17 | `sec_exp.tex:19` | "On the test patches F changes by at most 0.035 over κ∈[2,4]" holds for the 4 M2 patches. Table S6 (10 patches) has 0.069 on M53 (640,680), 0.285→0.354. | "M2 test patches" |
| F18 | `sec_exp.tex:15` | "(ratio 0.82 [0.71, 0.94])" beside 40.2 vs 55.1 nm (0.73). 0.82 is the bootstrap value, which `real_frc.md` says is biased coarse (replicates hold ~63% of the blocks). | Give 0.73 and explain |
| F19 | `sec_exp.tex:20` | "at a lower purity; Supplement Table 8": Table S8 has no purity. | See F6 |
| F20 | `sec_exp.tex:5` | "the 16 of 36 targets whose entry sequences differ": 17 differ (1 same kernels/other order, 16 different kernels). | 17 |
| F21 | `sec_exp.tex:12` | "A background wrong by a factor 2 … costs 18.6": that is background ×0.5. ×2 gives +4.3 [−1.3, 10.9]. "Halved photon counts cost the network 5.4": the N1.LD interval [−15.4, 2.5] contains 0. Sweeps use 200/100 frames. | Fix wording; give CIs |
| F22 | `tab_learned.tex:22` vs `supp_data.tex:17` | Net (b) is "no labels" in the caption, but its density and photon counts are resampled from development labels. | Reword |
| F23 | `sec_exp.tex:11` | The discrepancy rule is quoted at z=2 for both sets. At HD z=3 is better (−6.2 vs −8.0; `iss.md`). | Use the best z |
| F24 | `supp_data.tex:52`; `tab_noise.tex` | "The design of the parametrization study" under "Synthetic order selection". Table 9 Wilson intervals treat 432 trials as independent (36 targets × 12). | Fix label; cluster by target |
| F25 | supplement numbering | Supp "Theorem 3 (… Theorem 2 of the paper)". "by Thm. 2" means the supplement's Theorem 2 = the paper's Thm 1. Supp Prop. 2 ≠ paper Prop. 2. | Use S-prefixed or shared numbering |
| F26 | `supp_theory.tex:38`; `supp_sep.tex:20`; `supp_theory2.tex:14`; `tab_kappa.tex:16`; `tab_2013.tex:30`; `supp_data.tex:45,55`; `supp_sep.tex:50` | Dangling internal references in a double-blind PDF: "referee round IS2", "the report", "in the repository", repo paths. Overfull hboxes up to 69 pt (supp p. 11, baselines paragraph). | Delete or point to the anonymised archive |
| F27 | throughout | Terminology: "a priori default" (frozen after development on these simulations, uses labelled σ0), "label-free" (HD sets borrow the LD camera calibration, `supp_data.tex:5`), "held-out" (test patches / repeat frames / validation), "tuned" (crowdsource). | Define once; "frozen default" |
| F28 | `sec_exp.tex:12`; `sec_intro.tex:11` | Cost. ISS times are for the a priori configuration (0.363/1.78 s); Picasso's are defaults (0.0005–0.03 s wall; tuned Picasso rungs untimed). E-tuned ISS takes 0.52–6.54 s/frame (`iss.md` l.136,297,456,599). Astronomy run times (94–522 s per 100×100 px patch; `m2.md`) are in no table. "Two to four orders" is SMLM-only. | Add an astronomy line; time tuned Picasso |
| F29 | `sec_alg.tex:17` | "BIC is the same test at κ=√(n_par ln N)" is approximate (post-refit block χ² vs pre-refit entrant z). "Does not depend on the camera or the noise level" holds only if the noise law is right (background halved at fixed slope: −18.6 E; M92 residual scale 1.8–2.3). | Qualify |
| F30 | `sec_theory.tex:28` | R_K undefined; "2Ks²/N" is for 1-D (2 parameters per unit; 3–4 in imaging). | Define |
| F31 | main text | No reproducibility statement (code, scripts, data licences, hardware, seeds). Raw data are not redistributed (fetch + SHA-256, `supp_data.tex:3`); the challenge pages state no licence (CLAIMS §6). | Add a 3-line paragraph and anonymised link |
| F32 | `preamble.tex:14,33–38` | 7 pages are met only with tightened spacing (parskip 0, display skip 4 pt, textfloatsep 6 pt, nosep lists, small captions). Expect +0.3–0.5 page under the official style. | Pre-plan cuts: §5.1 → caption; §5.3 FRC → supplement; the 2013/learned/cost paragraph → 3 sentences |

**Order of work (cost).**
- Hours: F1, F4–F6, F9, F12–F27 (text and captions).
- About a day: F3 relabelling, F8 supplement section, F6 table columns, extending Table S1 to Δ ≥ 10–14σ.
- 2–3 days: F7 baseline re-runs, F2(a) network-in-the-loop, F10 multi-emitter baseline.

## 4. Can the draft be completed without the known open items?

- **Held-out repeat frames.** Yes, if abstract, intro and Table 2 treat the kernel fix as a post hoc diagnosis (F6). It is the single most valuable astronomy addition. If it fails, the narrative flips, so fix the abstract wording only after the result.
- **GPU run (DECODE).** Yes. No claim depends on it (Scope: no superiority over deep learning). Shorten the U-Net paragraph to 3 sentences.
- **Official style.** It cannot be completed until the style file is published. Plan for overflow (F32).

## 5. Claims against evidence

| # | Claim (where) | Evidence | Supports as worded? | Note |
|---|---|---|---|---|
| 1 | Network "answers both on its own" (abstract) | Fig. 1; Thm 1; Fig. 2; 16 noisy trials | Partly: noiseless exact-K only; noise needs a stop | F1 |
| 2 | Converges "at rate ε^{1/2}" to the sliding ISS; "theorem under stated hypotheses" (abstract; Contr. 1) | supp 1.3 outline; 1.4: 8 targets (slopes 1.50/1.00/0.5); H1 on 12/48 | Yes as hypothesis-conditional outline. Intro "proved": no. Rate is for duals and entry times (ε^{1/4} in L² output) | F3 |
| 3 | One unit per feature; Δ0 > benchmarks; 1154 targets down to 2σ | Thm 2 outline; Table S1 (720/720); 1154/1154 recounted | Numerics yes; theorem vacuous at tested Δ (stated); "2σ" should be 1.5σ; "ratio ≤2.4" wrong | F3, F13, F14 |
| 4 | "Closed-form skeleton predicts which units carry them" | Prop. 2: 1149/1154 binary vs 903/783/763 | No: binary event, oracle input, fitted 1154 circular | F4 |
| 5 | Window ∝ 1/s, oracle fit at first order | Prop. 3; `noisy_certificate.md` (ratio 1.0001; W·s≈13–14) | Yes: white Gaussian, 1-D, N=61, K=3, s≤0.1; fails at s=0.3; Poisson/imaging unproven; parameter-free stop is empirical (gap rule 71%) | — |
| 6 | Depth ≥3 selects; depth 1/2 do not | Fig. 2; `parametrization.md` C | Yes for gradient flow on 200 1-D targets with a margin proviso; fixed-step and Adam thin; not in supplement | F8 |
| 7 | "No training data and no regularization parameter" | Sec. 4 Inputs; Table 1 label-free; ablations | No: κ/BIC; labelled development frames for tuned rows; label-free loses to tuned Picasso on LD by 10–12 | F1 |
| 8 | Higher efficiency than tuned offset-corrected Picasso, 15/16 (1 tie) | Table 1, S2, S4, S7; reports | Counts right (E at two operating points); comparisons not independent; plotted metrics 13/3/0; RMSE 8/5/3 | F5 |
| 9 | M2 F 0.51 vs 0.25 on held-out patches | Table 2; `m2.md`; `freeze.json` | Yes (numbers exact; 4/4 patches); windows of one frame; "tuned crowdsource" = default; baselines lack high-pass | F7 |
| 10 | M53/M92 "fails until kernel matches stamp, which we diagnose" | Table 2, S8; `psf.md` | Failure and diagnosis supported; "until" is post hoc on the same patches; purity 0.25 vs 0.43 | F6 |
| 11 | Cost 10²–10⁴ × Picasso | `speed.md`; `iss.md` | SMLM only; astronomy cost absent | F28 |
| 12 | Parity with best explicit methods "with nothing to set" | Table 9 | Accuracy yes (82 vs 83–85). "Nothing to set" holds equally for MP+slide+BIC (84%); ISS slowest | F2, F11 |
| 13 | "The network follows the path" | `noise_order.md`: 16/16, rel. rms 1.3e-9/3.7e-8 | Yes, but 16 trials, ε=1e-8, 1-D | F2 |
| 14 | FRC: gain is mostly the count | `real_frc.md` | Yes (honest) | — |
| 15 | Learned baseline below tuned ISS | Table S5; `learned.md` | Only "CPU-scale U-Net"; not DECODE-class | F10 |
| 16 | "We can say which units and when" (conclusion) | — | No | F4 |

## 6. IJCAI-style review

**Summary.** The paper shows that gradient-flow training of a Gaussian-kernel network with a factorized "mass × amplitude" readout, started from a small scale ε, builds a path. Units enter one at a time when an accumulated correlation reaches a wall, slide onto features, and the rest stay dormant. It states a limit theorem for moving centres under hypotheses (H1–H4), with outlined proofs. It gives an outlined theorem for separated features, a first-order noise theory with a window ∝ 1/s, and an empirical contrast between readout depths. It then packages the path as an imaging algorithm. That algorithm is compared with Picasso on eight simulated SMLM sets and with crowdsource/photutils on SDSS patches of M2, M53 and M92, with a reported transfer failure and a post hoc fix.

**Strengths.**
1. A clean, checkable mechanism: a conservation law and mirror coordinate, a wall ∝ ε^{-1/2}, dormant units at ε^{3/2}μ⁻³. I verified the algebra and the first-order theory; the rates match experiments to 3–4 digits on 8 targets.
2. Unusual honesty. Scope and limitation paragraphs; dual order gives no accuracy gain; the FRC gain is count; the M92 loss is reported and diagnosed; the fix is labelled post hoc in §5.4 and Table S8.
3. Experimental hygiene. Frozen, hashed protocols; paired block bootstraps; a Picasso ladder of 12–15k configurations; defaults and tuned rows; ablations; a mismatch sweep. I found no wrong number in 309 intervals, 80 κ values and a recount of the 1154-target table.
4. A neat falsifiable observation: depth ≥3 vs 2 vs 1, with the margin law and the ln(1/θ)/ln(1/ε) prediction.
5. The M2 frozen test is well controlled and the effect is large.

**Weaknesses.**
1. The headline claim (implicit selection, "on its own", no threshold, no training data) is contradicted by the paper's own §3 (stopping needs a rule), the Table 1 label-free rows, and the labelled development data behind the tuned rows (F1).
2. The thesis is about RBF networks, but the application results concern a different object: a calibrated sliding greedy fitter with a high-pass, learned widths, exclusion zones and a bespoke stop. By the authors' ablations the accumulated-dual order is not what wins, and no network is trained on imaging data (F2).
3. The theory is outlined, not proved. Thm 1 holds on a minority of random targets (12/48 at margin 0.05); Thm 2 is vacuous at every tested separation; Prop. 2 is a table of counts; labelling is inconsistent between abstract, intro and supplement (F3, F4).
4. Baselines are thin. SMLM: one single-emitter baseline. Astronomy: "tuned crowdsource" is the default and no baseline gets the high-pass that accounts for most of the ISS's gain (F7, F10). The two real-data domains are, respectively, simulation plus FRC, and one frame per cluster.
5. Positioning and clarity. Prior art is missing or not on page 1; the title overgeneralises; the abstract and §5 are very dense (F9, F11, F12).
6. The practical improvement is hard to value: slower than Picasso by 10²–10⁴, slower than MP+slide on the synthetic benchmark, and unproven outside simulations.

**Questions for the authors.**
1. Does a factorized network trained by gradient flow on a small imaging patch (weighted, high-passed loss) form the units your algorithm returns? If not, what do the imaging results say about implicit order selection?
2. With the same high-pass, learned widths, exclusion and calibration, matching pursuit with sliding gives the same E/F. What does the accumulated dual contribute in practice?
3. How do crowdsource and photutils do on the high-passed patch and stamp, and with pipeline parameters tuned? How does a multi-emitter SMLM fitter compare?
4. What is Δ0 for K=2 (even numerically)? Does Table S1 still show exactness, and the decoupled entries, at 10–14σ?
5. How many random targets satisfy (H1)–(H4) at margin 0.01? How small must ε be for the first-order statements when |s_θ+σ_k| reaches 10⁴?
6. What stop lands in the noise window without knowing s or K? The gap rule gets 71% (synthetic) and −36/−13 points (SMLM).
7. Why does M92's κ=2 produce 868 detections for 638 stars? Is F still the right objective for catalogue use at purity 0.25?
8. What are the held-out repeat-frame results, and what would you claim if they are negative?
9. How robust is the selection result to a generic random initialisation (non-zero α(0)) and to Poisson noise?

**Scores.**
- Soundness 3/5
- Significance 3/5
- Novelty 3/5
- Clarity 2/5
- Overall 5/10 (marginally below the acceptance threshold)
- Confidence 3/5

**Recommendation.** Borderline, lean weak reject in the present form. I would move to 6 if the abstract and intro are re-scoped (F1, F4–F6, F9), the theory is relabelled (F3), the readouts supplement is added (F8), purity and the baseline fairness issues are fixed (F6, F7), and one network-in-the-loop imaging experiment is added (F2). The paper's honesty is unusual and I would like to see it accepted after revision.

## 7. What I could not check

- Proofs beyond what I re-derived: Steps 2–5 of the Thm 1 outline (Gronwall and matching estimates), Lemma 4 (coupled entry), Newton–Kantorovich constants, "Prop. (skeletons' accuracy)", and (H3) beyond the 48 targets.
- Berthier and Pillaud-Vivien (arXiv:2606.23198) and all other literature (offline): the novelty boundary and the missing citations in F11 rest on my recollection.
- Re-running the SMLM or astronomy pipelines, the 48/200/720-target experiments, or byte-identical regeneration. I only compared printed numbers and recounted from one raw file.
- That the freeze and test protocols were actually followed historically (I read `freeze.json` and code comments, not git history).
- The IJCAI-27 style and its true page count; bibliographic details (journal volumes, pages), which I could not verify offline.
- Figures at print scale (rendered at 110–200 dpi); data conventions (challenge matching, HUGS alignment) and data licences.
- The held-out repeat-frame test and DECODE (excluded by the brief).