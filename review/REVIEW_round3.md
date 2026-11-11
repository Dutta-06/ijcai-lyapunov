# Adversarial review, round 3 (5 refuters + judge), on PARADIGM.md version 3

Raw reports, verbatim. The judge's consolidated verdict comes first; the refuters' reports follow. Fixes and where they were applied are listed in `REVIEW_LOG.md`.

## Judge

I judged the five round-3 refuter reports, about 75 issues. After merging duplicates there are 49 confirmed issues: 1 major, 22 minor and 26 cosmetic. Nothing is fatal, 1 issue is rejected, and 5 points are left uncertain. My scripts and outputs are in <scratchpad>\round3\judge\ (written 'judge/' below).

**The project changed while I judged.** I judged against version 3 as reviewed (01:00:39, md5 7a8afd1d…). The run was not static:
- Between 01:31 and 02:14 another process edited the tests (test_averaging, test_langevin, test_transport), the experiments (exp_slow_manifold_slope, exp_orbit_orientation, exp_langevin_tradeoff, new exp_gibbs_reference), the results (table_slope.md, the tradeoff table, new gibbs_reference.md), LITERATURE_TRANSPORT.md, RESULTS_PARADIGM.md and PARADIGM.md itself (now 1719+ lines).
- That process applied almost every refuter-proposed fix.
- I diffed the 02:14 PARADIGM.md against the reviewed snapshot and checked the applied text against my own computations. It is correct except for the residual errors below.
- Each confirmed item gives the fix in v3 terms and a STATUS line saying what the 02:14 files already contain and what still needs changing.

**Residual errors in the applied text (02:14), each with an exact correction in its item:**
1. G(a) bullet: 'sixteen times less at k_e = 80' should be about 22 times (3.7e-8 measured).
2. §6b numerics: 'about 4.8 with the actual widths' and 'about 20 in the slope rows' no longer match the redefined ω_eff (2.6–4.8 and about 16).
3. H.7(e) circulation proof, step 5: 'for every ε, hence zero' gets the logic wrong. Zero follows by letting k_e → ∞, because the cyclic sum does not depend on k_e.
4. H.7a bullet: 'against 0.7% measured by the reviewer' should read 'asserted at 3%; 0.7% measured'.
5. Dictionary σ-mod row: the λ-convexity sentence now reads as if it applied to the W-block-only object.
6. Remark H.0 table: still no note that the n = 2, k_e = 10 entry is not computed by the cited test.
7. Stale companion and code text:
   - LITERATURE_TRANSPORT.md:130, 'given by the area that the reference orbit encloses';
   - RESULTS_PARADIGM.md:178, tolerances given as measurements;
   - averaging.py docstrings, 'unique periodic solution' and 'decays like exp(−k_e T₀)';
   - exp_orbit_orientation.py module docstring, which names b_odd only;
   - tests/test_transport.py, whose `__main__` block still precedes two tests.

**Key verifications:**
- G(b′) offset: the displayed offset is 1+S times too large. Actual/displayed is 0.26–0.30, against 1/(1+S) = 0.256–0.298. The corrected formula agrees to within 0.4–0.8%, and the truncated flow ascends L*.
- b_osc relative sizes: 0.918 (n = 1) and 4.444 (n = 2), against 0.416 and 4.054 for the frozen terms. b_osc scales as G², with a ratio of exactly 4.000 when G doubles.
- H.7 table numbers: 0.537 → 0.247; 0.038 → 0.019 and 0.0225 → 0.0115; Stokes to 6.8e-15; parity to 2.2e-15; non-reciprocity −0.07107 against −0.07035 and −0.04706 against −0.04637; full matrix within 1.8%; Richardson residual 1.04%; b_osc value test 0.73%.
- f(x) = x² gives a second periodic orbit and blow-up, for k_e = 10, 40 and 160.
- Floquet rate under ω-rescaling at a fixed PE level: 4.66e-2 → 1.17e-4 for ω = 1 → 32.
- Saddle eigenvalues: −5.4e-5 for (𝔄) and −2.2e-4 exactly averaged, both in s⁻¹.
- Sampling floor of W₂ between independent 32-point samples: mean 0.670.
- L̄(μ*_τ) = 0.111183.
- The F(ii′) test still passes with the leak term deleted (worst ratio 0.308).
- The satisfying-run decomposition: (b′)/VarPro = 0.9999–1.0000; frozen-parameter ratio 0.951–0.956; weights running at frozen geometry give 0.969–0.970.

### Confirmed issues

#### [G1-offset] major (G)

**Location.** §6b Proposition G(b′), first display 'a − a*(ξ) = (δ/g_a)M⁻¹Da*Γ∇L*(ξ) + h.o.t.'; Derivation step 3 'h − a* = −(k_e/(mg_a))M⁻¹Dh ξ̇'

**Issue.** The displayed slow-manifold offset is the S → 0 truncation. In the regime that (b′) is stated for (S = O(1)) it is off by a factor of about 1 + S, and it contradicts the displayed flow. Substituted into ∇_ξL̄ it gives the first-order Neumann truncation, which can ascend L*. The (b′) flow display itself is correct.

**Evidence.** judge/g_offset.py: joint averaged flow (a) in the slope configuration, after the boundary layer.
- δ = 0.47, S = 2.35–2.90:
  - actual a − a* = −0.0373 to −0.0399;
  - displayed formula −0.124 to −0.155 (actual/displayed 0.259–0.300, against 1/(1+S) = 0.256–0.298);
  - corrected (δ/g_a)M⁻¹Da*(Γ⁻¹+(δ/g_a)Da*ᵀDa*)⁻¹∇L* agrees to within 0.4–0.8%.
- δ = 0.047: the displayed formula is 20–25% off; the corrected one agrees to within 1%.
- At ξ₀ (δ = 0.47), dL*/dt is +3.66e-6 along the truncated flow and −1.51e-6 along (b′).

**Fix.** Replace '(b′) *Three scales, under (H0), (H1), (H2a′) and (H2b).* After a boundary layer, the weights follow the slow manifold:   a − a*(ξ) = (δ/g_a)M⁻¹Da*Γ∇L*(ξ) + h.o.t.   The geometry follows, to first order,' with '(b′) *Three scales, for the averaged system (a), under (H0), (H1), (H2a′) and (H2b).* After a boundary layer, the weights follow the slow manifold:   a − a*(ξ) = −(k_e/(mg_a))M⁻¹Da*·ξ̇ + h.o.t. = (δ/g_a)M⁻¹Da*(Γ⁻¹ + (δ/g_a)Da*ᵀDa*)⁻¹∇L*(ξ) + h.o.t.   This reduces to (δ/g_a)M⁻¹Da*Γ∇L*(ξ) when S ≪ 1. The geometry follows, to leading order in 1/R and for every S,'.

In step 3, after 'h − a* = −(k_e/(mg_a))M⁻¹Dh ξ̇.', append ' This must be evaluated with the (b′) velocity, not the VarPro velocity.'

STATUS: applied verbatim in the 02:14 file and verified correct. No further change.

#### [G2-weight-block] minor (G)

**Location.** §6b (b′) 'After a boundary layer, the weights follow the slow manifold'; Numerics row `satisfying`, G(a) cell 'geometry block cosine 0.999 / 0.998 (see below)'; bullet 'The weight block. In the `satisfying` run the weights sit on the slow manifold. Their net drift is small: the predicted weight speed is 0.17–0.28 of the geometry speed'

**Issue.** Merges G #2 and G #7.
- In the closed loop the weights sit at the O(1/k_e)-shifted equilibrium of Proposition H.7, 10–100 times farther from a* than the slow-manifold offset.
- The quoted 0.17–0.28 is the predicted restoring velocity, not the net drift.
- The table cell leaves out the failing joint G(a) numbers.

**Evidence.** judge/n5_copy.py on raw_sim_satisfying.npz, t ≥ 515 s: |a_pm − a*| = 3.1e-3 to 5.9e-3, against a slow-manifold offset of 4.0e-5 to 2.9e-4.

judge/n5b_copy.py: the actual weight/geometry speed is 0.001–0.008, against a predicted 0.167–0.300.

table_satisfying.md: joint cosine 0.9848 / 0.9624, ratio 0.956 / 0.935.

**Fix.** Three replacements:
- G(a) cell → 'joint cosine 0.985 / 0.962, size 0.96 / 0.94; geometry block 0.999 / 0.998; the weight block fails (see below)'.
- After the (b′) metric-tensor sentence, add: '- In the closed loop, the weight equilibrium is further shifted by O(ω_eff/k_e) (Proposition H.7). This shift can exceed the offset displayed above: in the `satisfying` run it is 3–6·10⁻³, against an offset of 4·10⁻⁵–3·10⁻⁴ (review round 3).'
- Weight-block bullet → '**The weight block.** In the `satisfying` run the weights sit at a* + O(1/k_e), within 0.1% of a*. Their actual drift is tiny: 0.001–0.008 of the geometry speed in the G⁻¹ norm. The prediction −(m/k_e)g_a∇_aL̄ at the period-mean weights is 0.17–0.30 of the geometry speed, 40–200 times the actual drift. It is the restoring velocity of the O(1/k_e) offset, which (a) does not contain. So the weight block fails (cosine −0.21 / 0.50, carrying 92–95% of the residual), while the geometry block matches.'

STATUS: applied in the 02:14 file and correct.

#### [G3-satisfying-ratio] minor (G)

**Location.** §6b 'The size ratios': 'In the `satisfying` run, the centred ratio of 0.97 is that same correction (frozen-parameter ratio 0.956 at k_e = 20, 0.996 at k_e = 80; review round 2) plus about 0.6% from S.'

**Issue.** The attribution is wrong in size and in sign.
- With Γ = I the slope S can only lower the predicted speed, and here it changes it by less than 0.01%.
- The +1.4% comes from the weights' oscillation within each period.

**Evidence.** judge/n5_copy.py: (b′)/VarPro = 0.9999–1.0000. The frozen-parameter ratio is 0.951–0.956 along the run, while the closed loop gives 0.9693–0.9703.

judge/n5b_copy.py, weights running under (L) at frozen geometry:
- g_a = 13: 0.9691, 0.9701, 0.9702;
- g_a/4: 0.9600, 0.9587, 0.9563.
Centred windows of ±2 to ±40 periods all give 0.969–0.972.

**Fix.** Replace the sentence with: 'In the `satisfying` run the centred ratio of 0.97 is the error-channel factor (frozen-parameter ratio 0.951–0.956 along the run; 0.996 at k_e = 80) times a +1.4% averaging correction from the weights' oscillation within each period, the analogue of Remark H.7a. With the weights running under (L) at frozen geometry the ratio is 0.969–0.970, and with g_a/4 it is 0.956–0.960 (review round 3). The slope correction is below 0.01% here, because ∇L* is nearly orthogonal to the range of Da*ᵀ.'

STATUS: applied in the 02:14 file and correct.

#### [G4-slope-boundary-layer] minor (G)

**Location.** §6b Numerics table, slope rows (δ = 0.47: 'cosine 0.9994, ratio 1.10'; '(b′): cosine 0.9967, ratio 1.17'); experiments/exp_slow_manifold_slope.py start a(0) = a*(ξ₀)

**Issue.** The slope comparison included two errors that the text does not mention:
- the weight boundary layer, because the run started off the slow manifold;
- the averaging error, with ε_avg ≈ 0.08 at δ = 0.47.
Most of the 17% gap to (b′) was therefore not an error of (b′).

**Evidence.** The script was already changed (01:34) to start on the corrected slow manifold and to add a run with all gains divided by 4. The regenerated table_slope.md (01:46) gives:
- δ = 0.47: closed/(a) 1.104; closed/(b′) 1.121, cosine 0.9986; (a)/(b′) 1.015;
- all gains ÷ 4: closed/(a) 1.049; closed/(b′) 1.065;
- δ = 0.047: 1.044 / 1.047;
- δ = 0.0047: 1.044 / 1.044.
This matches the refuter's decomposition.

**Fix.** Replace the three v3 slope rows with:
'| `slope`, δ = 0.47 (m = 1) | 81–178 | 2.0–2.8 | cosine 0.9992, size 1.10 | (b) fails: ratio 0.42, cosine 0.93, and the centre moves the opposite way; (b′): cosine 0.9986, ratio 1.12 |'
'| `slope`, δ = 0.47, all gains ÷ 4 | 81–178 | 2.0–2.8 | cosine 1.0000, size 1.05 | (b) fails: ratio 0.40, cosine 0.92; (b′): cosine 0.9999, ratio 1.065 |'
'| `slope`, δ = 0.047 | 654–1779 | 0.20–0.32 | cosine 1.0000, size 1.04 | (b): ratio 0.85, cosine 0.9986; (b′): ratio 1.05, cosine 1.0000 |'
'| `slope`, δ = 0.0047 | 13162–17790 | 0.020–0.023 | cosine 1.0000, size 1.04 | (b): ratio 1.02; (b′): ratio 1.04; both cosines 1.0000 |'

Add the 'Error decomposition' bullet: ‖(a)‖/‖(b′)‖ = 1.015; ‖closed‖/‖(a)‖ = 1.10 → 1.05 at gains ÷ 4; the remaining ≈ 4.4% is the error-channel factor.

STATUS: code, table and document already applied (01:34–02:14) and correct.

#### [G5-Lbar-monotone] minor (G)

**Location.** §6b Proposition G(a), third bullet: 'Along the closed loop, its value at period-averaged parameters decreases, up to O(ε_avg) oscillation within a period.'

**Issue.** Near equilibrium, L̄ at period means can rise steadily over many periods as the weights settle at the O(1/k_e)-shifted equilibrium. This is not a within-period oscillation.

**Evidence.** judge/g5/first_case.py and n4b, the refuter's scalar closed loop. Slope configuration, δ = 1e-5:
- k_e = 20, a(0) = a* − 0.01232: L̄ rises from 2.12435e-3 to 2.12516e-3 (+8.13e-7), increasing in 35 of 46 steps; final a − a* = +0.0080.
- k_e = 80: +3.7e-8, increasing in 80 of 94 steps.
The measured ratio is 22, against 16 for an exact 1/k_e² law.

**Fix.** Replace the bullet with: 'L̄ is non-increasing along the averaged flow. Along the closed loop, L̄ at period-averaged parameters decreases while the averaged dissipation (m/k_e)‖∇L̄‖²_G dominates the O(ω_eff/k_e) corrections of Proposition H.7 and Remark H.7a. Near equilibrium it can increase by O((ω_eff/k_e)²), as the weights settle at their shifted equilibrium (review round 3: +8·10⁻⁷ over 35 periods at k_e = 20, and 3.7·10⁻⁸ at k_e = 80).'

STATUS: applied in the 02:14 file, with one wrong factor. Replace 'sixteen times less at k_e = 80' with 'and 3.7·10⁻⁸, about 22 times less, at k_e = 80'.

#### [G6-step1] minor (G)

**Location.** §6b Derivation step 1: 'and e = r(x_m)/k_e + O(ω_eff·sup‖r‖/k_e²)'

**Issue.** The v3 ω_eff contains the Lipschitz constant of f − f̂ but not the bandwidth of r along the reference. So sup‖ṙ‖/k_e² can be far larger than the stated remainder scale.

**Evidence.** judge/g_step1.py: f = 0.01 sin 100x, f̂ = 0, x_m = sin t, s_min = 0.25, k_e = 100. Then k_e/ω_eff = 25 and sup‖r‖/k_e = 1e-4, so (H2a′) holds.
- sup|e − r(x_m)/k_e| = 7.1e-5, against a claimed scale of 4.0e-6 (17.8 times larger).
- The RMS of the remainder is 54% of the RMS of r/k_e.

**Fix.** Write 'e = r(x_m)/k_e + O(sup‖ṙ‖/k_e² + sup‖∂_x r‖·sup‖r‖/k_e²)'. Add the bandwidth sup_t‖(d/dt)r(x_m(t))‖/sup‖r‖ as a fourth entry of ω_eff in (H2a′).

STATUS: applied in the 02:14 file and correct.

#### [G7-gain-condition] minor (G)

**Location.** §6b Status, 'What the designer supplies': 'g_c, g_s ≪ g_a·min(m₀/‖Γ∇²L*‖, 1/‖Da*Γ^{1/2}‖²)'; the same in LITERATURE_TRANSPORT.md §G

**Issue.** The condition counts Γ twice, because g_c = δγ_c and g_s = δγ_s. With Γ = γI it reads δγ² ≪ …, which is off by a factor γ from (H2b) and (H2c).

**Evidence.** (H2b) is equivalent to δ‖Γ^{1/2}∇²L*Γ^{1/2}‖ ≪ g_a m₀, and (H2c) to δ‖Da*Γ^{1/2}‖² ≪ g_a.

**Fix.** Replace with 'δ ≪ g_a·min(m₀/‖Γ^{1/2}∇²L*Γ^{1/2}‖, 1/‖Da*Γ^{1/2}‖²). Equivalently, ‖G_ξ^{1/2}∇²L*G_ξ^{1/2}‖ ≪ g_a m₀ and ‖Da*G_ξ^{1/2}‖² ≪ g_a, with G_ξ := δΓ.' Make the same change in LITERATURE_TRANSPORT.md.

STATUS: applied in both files and correct.

#### [G8-numerics-reporting] minor (G)

**Location.** §6b Numerics heading '(`results/variable_projection`; centred differences of period means)', column 'G(b) and the induced-metric flow (b′)', and the (H2a) columns of the tables

**Issue.** Severity raised from cosmetic. Three reporting problems:
- The slope rows use end-to-end displacements, not centred differences.
- (b′) is evaluated only in the slope rows.
- Only the second inequality of (H2a′) is reported. The first, k_e ≫ ω_eff, gives k_e/ω_eff = 20/(2.5/0.25) = 2.0 in the m = 2 and m = 4 runs, so (H2a′) holds only marginally there.

**Evidence.** exp_variable_projection.py: sup|ẋ_m| = 2.5, box s_min = 0.25, k_e = 20. judge/g_omega_eff2.py: with the actual widths and the bandwidth term now in ω_eff, k_e/ω_eff is 3.6–4.4 (`violating`) and 2.6–4.8 (`satisfying`, mostly 2.6; the residual bandwidth is about 7.7). judge/g_omega_eff.py: the slope rows give about 16 (16.1–17.1 at ξ₀).

**Fix.** Use the heading and bullets now in the file:
- the first two rows use centred differences;
- the slope rows use first-to-last displacements and start on the slow manifold;
- (b′) is evaluated only in the slope rows.

STATUS: applied, except two numbers that no longer match the redefined ω_eff. Replace 'k_e/ω_eff is about 2.0 in the first two rows (about 4.8 with the actual widths), so (H2a′) holds only marginally there. It is about 20 in the slope rows.' with 'k_e/ω_eff is 2.0 in the first two rows with the box's s_min, and 2.6–4.8 with the actual widths (the residual bandwidth, about 7.7, then dominates in the `satisfying` run), so (H2a′) holds only marginally there. It is about 16 in the slope rows.'

#### [G9-first-order-wording] cosmetic (G)

**Location.** §6b (b′) 'The geometry follows, to first order,'; Numerics '(b′) is the correct first-order slow flow'; §8 item 4 'higher-order corrections in S'

**Issue.** (b′) is leading order in 1/R and exact in S, so the remaining corrections are in 1/R, not in S.

**Evidence.** Derivation step 3: the resummed form solves the invariance equation once Dh − Da* and the quadratic term are dropped, and both of those are O(1/R).

**Fix.** Make three replacements:
- 'to first order,' → 'to leading order in 1/R and for every S,';
- 'correct first-order slow flow' → 'correct leading-order slow flow';
- §8 'higher-order corrections in S' → 'the O(1/R) corrections, which involve D²a* and ∂_ξM'.

STATUS: applied in the 02:14 file and correct.

#### [G10-RS-gain-shape] cosmetic (G)

**Location.** §6b '(H2b) does not imply (H2c). The product R·S = … is fixed by the problem.'

**Issue.** R·S does not depend on g_a, δ or k_e, but it does depend on the shape of Γ.

**Evidence.** judge/g_rs.py, at ξ₀: R·S = 165.7, 355.5 and 1149.2 for γ_s/γ_c = 1e-4, 1 and 1e4.

**Fix.** '… does not depend on g_a, δ or k_e. It is fixed by the problem and by the shape of Γ: in the slope configuration it ranges from 166 to 1149 over γ_s/γ_c ∈ [10⁻⁴, 10⁴] (review round 3).'

STATUS: applied in the 02:14 file and correct.

#### [G11-undefined-symbols] cosmetic (G)

**Location.** §6b 'O(ε_avg)', '(ω_eff/k_e, ε_avg, S, 1/R) → 0', '[0, b/ε_avg]', '𝒢′ ⋐ 𝒢', 'M ⪰ m₀I on 𝒢', 'R := g_a m₀/(δ‖Γ^{1/2}∇²L*Γ^{1/2}‖)'

**Issue.** ε_avg and 𝒢 are never defined in §6b, and §6c reuses 𝒢 for a functional. R is defined pointwise, but the limit statement needs a sup.

**Evidence.** grep: ε_avg is used 6 times and defined nowhere; 𝒢 is used in §6b and §6c with different meanings.

**Fix.** Add to the Setting of §6b:
- '𝒦 is the set of geometries ξ on which (H1) holds, and 𝒦′ ⋐ 𝒦 is compact. (𝒦 is not the functional 𝒢(μ) of §6c.)';
- 'R := g_a m₀/(δ·sup_{𝒦′}‖Γ^{1/2}∇²L*Γ^{1/2}‖)';
- 'ε_avg := (m/k_e)λ_max(G^{1/2}∇²L̄G^{1/2})/ω_eff, the second ratio in (H2a′). The tables report its inverse, with ω_ref in place of ω_eff.'
Rename 𝒢 to 𝒦 throughout §6b.

STATUS: applied in the 02:14 file and correct.

#### [H7a-global-uniqueness] minor (H7-math, H7-code)

**Location.** §6c Proposition H.7 preamble 'there is k₀, independent of μ'; (a) 'has a unique T₀-periodic solution e_μ. It attracts every solution exponentially'; proof 'Exponential attraction follows from the linearisation −k_e + Lip(r_μ).'; src/core/averaging.py docstrings

**Issue.** The hypotheses (f ∈ C², autonomous) allow superlinear f. Then the uniqueness and global attraction of the periodic error are false, and only a local statement (in a ball) is proved.

**Evidence.** judge/h7a_global.py: n = 1, f(x) = x², x_m = sin t, f̂_μ = 0. There is a second, repelling periodic orbit, periodic to 1e-12:
- k_e = 10: e ∈ [7.89, 11.91];
- k_e = 40: e ∈ [37.97, 41.98];
- k_e = 160: e ∈ [157.99, 161.99].
Starting 0.5 above its maximum, |e| reaches 1e6 at t = 0.16, 0.071 and 0.026 respectively. The code is unaffected, because it starts from e = 0.

**Fix.** Preamble → 'Then for every R > 0 there is k₀(R), independent of μ, such that for k_e ≥ k₀(R) the following hold.'

(a) → 'has exactly one T₀-periodic solution e_μ with sup_t‖e_μ‖ ≤ R. It attracts exponentially every solution with ‖e(0)‖ ≤ R, and'.

After 'uniformly in t and μ.' add: 'If f is one-sided Lipschitz, ⟨f(x) − f(y), x − y⟩ ≤ L_f‖x − y‖², then for k_e > L_f + a_max√n·e^{−1/2}/s_min the periodic solution is unique and attracts every solution. Without such a condition uniqueness is only local: for f(x) = x² there is a second, repelling periodic solution near e ≈ k_e, and solutions above it blow up.'

Proof: replace 'Exponential attraction follows from the linearisation −k_e + Lip(r_μ).' with: 'The tube {‖e‖ ≤ R} is forward invariant once k_eR > sup_{‖e‖≤R}‖r_μ(x_m + e)‖. Inside it, solutions contract at rate k_e − Lip_R(r_μ), which gives existence, uniqueness and attraction there. The global statement uses L_f in place of the local Lipschitz constant of f.'

STATUS: the document part is applied in the 02:14 file and correct. Still to do, in src/core/averaging.py:
- module docstring: 'k_e large it has a unique periodic solution e_μ' → 'k_e large it has a unique periodic solution e_μ near 0';
- periodic_error docstring: 'the transient decays like exp(−k_e T₀)' → 'the transient decays like exp(−(k_e − Lip r_μ)T₀) per period'.

#### [H7e-remainder-derivative] minor (H7-math)

**Location.** §6c H.7(e) 'Its antisymmetric part is (1/k_e²)(A_even + 2K) + O(k_e⁻³)'; proof '(e) Differentiate (b) in μ.'

**Issue.** (b) bounds the remainder uniformly in ϑ and μ and in its ϑ-gradients, but not its first variation in μ. So 'differentiate (b) in μ' is not justified as written. The conclusion is true.

**Evidence.** I re-derived the refuter's circulation argument:
- For quadratic U the Green integrand is affine, so the circulation equals −(ε²/2)Σ_cyc(S − Sᵀ) at the centroid.
- The remainder contributes at most 6ε·sup|R|, since Σ_j|dw_j| = 2ε on each edge.
The exact cyclic sums also converge at the claimed rate (judge/h7_table.py): k_e²·cyclic sum −0.07107 against −0.07035 (n = 1) and −0.04706 against −0.04637 (n = 2) at k_e = 80.

An alternative patch also works: the μ-derivative w of e_μ solves the linear periodic equation ẇ = (∇r_μ(x_m+e_μ) − k_e)w − Φ̃(x_m+e_μ;ϑ′), and the argument of (a) bounds its expansion with f ∈ C² and x_m ∈ C².

**Fix.** Use the circulation route, which is now in the file: restate (e) for the k_e⁻² coefficient, and prove the criterion through the triangle T_ε of measures μ_{a,b}.

STATUS: applied in the 02:14 file, with one logical slip in the last step. Replace '5. So if P_μ were a first variation up to O(k_e⁻³), the cyclic sum would be O(1/(εk_e)) for every ε, hence zero.' with '5. So if, for all large k_e, P_μ were a first variation up to O(k_e⁻³), the cyclic sum would be O(1/(εk_e)). It does not depend on k_e, so letting k_e → ∞ at fixed ε gives zero.'

#### [H0-gains] minor (H7-code, H7-math, H8-H0-H6)

**Location.** §6c Remark H.0: 'The first corrections are O(1/k_e) relative to the drift and do not depend on the gains. Proposition H.7 gives them exactly at frozen μ.'; experiments/exp_orbit_orientation.py module docstring 'Θ₊ − Θ₋ ≈ −2 J⁻¹ b_odd'

**Issue.** Along the real loop, the oscillation term b_osc of Remark H.7a is of the same order and quadratic in G, so the first corrections do depend on the gains. The sentence contradicts Remark H.7a. H7-code proposed major; I downgrade to minor, because:
- the correct statement is already in H.7a;
- the H.0 table is explicitly labelled frozen-μ;
- no theorem depends on the sentence.

**Evidence.** judge/h0_bosc.py, values of k_e × relative size:

| instance | frozen terms | b_osc | full first order |
|---|---|---|---|
| n = 1, Lg instance, G = diag(4,1,1) | 0.416 | 0.918 | 0.898 |
| n = 2, Av instance (seed 3), G = I | 4.054 | 4.444 | 4.133 |

Doubling G multiplies b_osc by 4.0000. judge/h0_split.py: for n = 2, U_even alone gives 4.034 and U_odd alone 1.528.

**Fix.** Replace the bullet with: '**The first corrections.** At frozen μ they are O(1/k_e) relative to the drift, do not depend on the gains, and are given exactly by Proposition H.7. Along the loop, the within-period oscillation of the parameters adds the second-order averaging term b_osc of Remark H.7a. Its size relative to the frozen terms is O(|G|T₀), so it does depend on the gains. In the formal limit below, both are O(k_e⁻²). Measured relative errors of the frozen-μ averaged drift against (𝔄), with b_osc excluded:'

After the table add: 'At the gains of these instances b_osc is not small: k_e × its relative size is 0.92 (n = 1, G = diag(4,1,1)) and 4.4 (n = 2, G = I), against 0.42 and 4.05 for the frozen terms (review round 3).'

STATUS: the document part is applied in the 02:14 file and correct (its 'U_odd alone would give 1.5' checks out at 1.528). Still to do: in experiments/exp_orbit_orientation.py, module docstring line 10, replace 'Θ₊ − Θ₋ ≈ −2 J⁻¹ b_odd,     b_odd,i = (1/k_e²) G ∇_ϑ U_odd(ϑ_i),' with 'Θ₊ − Θ₋ ≈ −2 J⁻¹ (b_odd + b_osc),     b_odd,i = (1/k_e²) G ∇_ϑ U_odd(ϑ_i), b_osc the oscillation term of PARADIGM.md Remark H.7a,'.

#### [H7-area-claims] minor (H7-math, H7-code)

**Location.** §8 item 2 '*Proved:* … Its antisymmetric interaction is generated by the area that the reference orbit encloses.'; §0 'Its kernel is the area enclosed by the reference orbit (Stokes).'; LITERATURE_TRANSPORT.md:130

**Issue.** The antisymmetric part is A_even + 2K.
- A_even needs no enclosed area, and it is the only part for n = 1.
- K is an integral over the enclosed region, not the area.
- 'Proved' overstates the status: the non-vanishing of the cyclic sums is checked in instances, not proved.

**Evidence.** For n = 1 the cyclic sum is −0.0704, with K ≡ 0. For n = 2 it splits into a state-shift part of −0.089 and a circulation part of +0.043.

Σ_k det[∇Φ̃_k(ϑ), ∇Φ̃_k(ϑ′)] = (aᵀa′)det[∇φ, ∇φ′]. So K = 0 whenever aᵀa′ = 0 or c = c′ (the two gradients are then parallel), whatever the enclosed area.

**Fix.** In §0 replace the second-part bullet with: '- The second part is an antisymmetric interaction between kernels (n ≥ 2). Its kernel is an integral of Σ_k det[∇Φ̃_k(ϑ), ∇Φ̃_k(ϑ′)] over the region that the reference orbit encloses (Stokes).'

In §8 replace the '*Proved:*' bullet with:
'- *Proved at a frozen measure:* the averaged adaptive dynamics is a semi-gradient flow plus a lag term (Proposition H.7). Along the real loop a second-order averaging term is added (Remark H.7a).
- Its interaction has two antisymmetric parts: the state-shift part A_even, which is present for every n and is the only part for n = 1; and the circulation kernel K, an integral over the region that the reference orbit encloses (n ≥ 2).
- Non-reciprocity is checked in the instances tested, not proved in general.'

STATUS: PARADIGM.md is applied in the 02:14 file and correct. Still to do: in LITERATURE_TRANSPORT.md:130, replace '- It includes an antisymmetric kernel–kernel interaction given by the area that the reference orbit encloses (Stokes).' with '- For n ≥ 2 it includes an antisymmetric kernel–kernel interaction, the circulation kernel, an integral over the region that the reference orbit encloses (Stokes). The state-shift part is non-reciprocal as well, and it is the only source for n = 1.'

#### [H7c-TD-analogy] minor (H7-math)

**Location.** §0 'ignoring that this distribution depends on the parameters, as temporal-difference learning does'; H.7(c) 'This is a semi-gradient, as in temporal-difference learning [S10], [S11].'; LITERATURE_TRANSPORT.md [S11] annotation

**Issue.** The citation is misused. In [S11] the data distribution is the stationary law of a fixed policy and does not depend on the parameters; TD is a semi-gradient because it freezes the bootstrapped target. What H.7(c) freezes is the parameter-dependent data distribution. That is the structure of performative prediction (Perdomo, Zrnic, Mendler-Dünner and Hardt, ICML 2020).

**Evidence.** Standard definitions: Sutton–Barto 2nd edition §9.3, and Tsitsiklis–Van Roy 1997 (policy evaluation for a fixed Markov chain).

The analogy itself is sound. The defect −(1/k_e²)∫(Φ̃_ϑ·∇)(½‖r_μ‖²)dρ is exactly what the semi-gradient leaves out: −(1/k_e)V_L + (1/k_e²)U_even = −(1/k_e)δJ^cl/δμ − (1/k_e²)∫(Φ̃_ϑ·∇)(½‖r_μ‖²)dρ + O(k_e⁻³). I checked analytically that the antisymmetric part of the defect's second variation equals A_even.

**Fix.** H.7(c): replace the TD bullet with: 'This is a semi-gradient: the update freezes part of its own parameter dependence. In TD learning the frozen part is the bootstrapped target, and the data distribution does not depend on the parameters [S10], [S11]. Here the frozen part is the data distribution ρ^cl_μ. This is the structure of performative prediction [S14].'

§0: replace 'as temporal-difference learning does' with 'as in performative prediction. Like semi-gradient TD, the resulting update is not the gradient of any objective'.

Also:
- add [S14] = Perdomo, Zrnic, Mendler-Dünner, Hardt, Performative prediction, ICML 2020, to LITERATURE_TRANSPORT.md;
- annotate [S11] as 'holds the bootstrapped target fixed; the data distribution is that of a fixed policy';
- add [S14] to the §8 context and the §9 reference list.

STATUS: applied in PARADIGM.md and LITERATURE_TRANSPORT.md (02:13–02:14) and correct.

#### [H7a-3pct-claim] minor (H7-code)

**Location.** §0 'The prediction is accurate to 3% at k_e = 40'; §0 table row H.7 'to 3% (cosine 0.994, k_e = 40)'; §8 item 2 'to 3%'

**Issue.** The 3% is the size ratio only. As a vector the agreement is 11% at k_e = 40.

**Evidence.** results/orbit_orientation/table_sigma2_slow16.md, regenerated at 01:35, gives the vector error ‖meas − pred‖/‖pred‖ as 0.364, 0.195 and 0.109 at k_e = 10, 20 and 40. This is consistent with ratio 0.969 and cosine 0.9944: √(r² − 2r·cos + 1) = 0.109.

**Fix.** §0 bullet: 'At k_e = 40 the prediction matches to 3% in size and 11% as a vector (cosine 0.994). The vector error halves with each doubling of k_e (`exp_orbit_orientation.py`).'

In the §0 table, §8, H.7a and §7, write '3% in size and 11% as a vector'. Add the column 'full: vector error' (0.36, 0.20, 0.11) to the H.7a table.

STATUS: applied in the 02:14 file and correct.

#### [H7-status-and-checked-table] minor (H7-math, H7-code)

**Location.** §0 table row H.7 status 'expansion, Stokes identity, parity and non-reciprocity checked to about 1%'; H.7 'Checked' header 'five tests' and rows ('both orientations', 'agreement to 1e-6', 'to 1e-10', 'reciprocal to 1e-14'); results/RESULTS_PARADIGM.md:178

**Issue.** The status line and the Checked table misreport what is measured.
- The expansion is checked through its rate: 25% for n = 1 at k_e = 80.
- Stokes and parity hold to machine precision, so 1e-6 and 1e-10 are tolerances, not measurements.
- The clockwise orbit gives different numbers from the counter-clockwise one.
- S₀ is symmetric by construction, so '1e-14' tests nothing.
- The file has seven tests, not five.

**Evidence.** judge/h7_table.py:
- expansion 0.5372 → 0.2468 (n = 1), 0.0382 → 0.0188 (n = 2, +), 0.0225 → 0.0115 (n = 2, −); Richardson residual for n = 1 is 1.04%;
- Stokes 6.8e-15; K + Kᵀ = 4e-16; parity 2.2e-15 / 7e-16;
- cyclic sums agree to 1.0% and 1.5%; full matrix within 1.78% and 1.81%;
- S₀ − S₀ᵀ is exactly 0.

**Fix.** §0 status → 'Frozen-measure part proved. Stokes identity and parity hold to machine precision. The expansion is checked through its O(1/k_e) convergence. Non-reciprocity of the exact second variation is checked to 3% (asserted; 1–1.5% measured).'

Checked-table rows →
- '0.54 → 0.25 for n = 1, 0.038 → 0.019 for n = 2 counter-clockwise, and 0.022 → 0.011 clockwise; a Richardson step at k_e = 80 and 160 brings n = 1 within 3% (asserted)';
- 'agreement to 7e-15 relative (asserted at 1e-6)';
- 'to 2e-15 relative (asserted at 1e-10)';
- 'the full matrix k_e²S + k_eS₀ is within 5% of S_even + K …'.
Header → 'seven tests: four for H.7, one for Remark H.0, two for Remark H.7a'.

STATUS: PARADIGM.md is applied in the 02:14 file and correct. Still to do: in results/RESULTS_PARADIGM.md:178, replace '| circulation kernel K | ≡ 0 (1e-12) | antisymmetric (1e-12); equals the enclosed-area integral to 1e-6 |' with '| circulation kernel K | ≡ 0 (1e-15) | antisymmetric (4e-16); equals the enclosed-area integral to 7e-15 (asserted at 1e-6) |'.

#### [H7-tests-weak] minor (H7-code, H8-H0-H6)

**Location.** tests/test_averaging.py non-reciprocity test (cyclic sums only, 10% tolerance); expansion test n = 1 threshold 0.3; tests/test_langevin.py 'errs[-1] < 0.35*errs[0]'; no value test of b_osc (§7 row 'H.7a … odd under orbit reversal (1e-6)')

**Issue.** At review time the tests did not test what the document claimed:
- the 1/k_e rate was not enforced;
- the full S_even + K matrix was not checked;
- the n = 1 expansion accepted a ±24% scaling of U_even;
- the value of b_osc was never tested, and −b_osc would have passed.

**Evidence.** Code fixes are already in tests/ (01:31–01:49):
- test_langevin asserts 0.4 < ratio < 0.62;
- the expansion test adds a Richardson step (residual 1.04%);
- the non-reciprocity test asserts the full matrix within 5% (measured 1.8%) and cyclic sums within 3% (measured 1.0% and 1.5%);
- the new test_oscillation_correction_matches_a_one_period_flow measures 0.73% (judge/osc_value.py), against 199.6% for −b_osc.
All 13 averaging and Langevin tests pass.

**Fix.** No further code change; the test fixes are applied. The document was updated in the §7 row for H.7a and the Checked table.

STATUS: one wording problem remains in the 02:14 file (Remark H.7a). Replace 'Its value matches an independent one-period near-identity flow to within 3% at k_e = 80, against 0.7% measured by the reviewer' with 'Its value matches an independent one-period near-identity flow at k_e = 80 (asserted to 3%; 0.7% measured)'.

#### [H0-table] cosmetic (H7-code, H8-H0-H6)

**Location.** §6c Remark H.0 table: 'n = 1 | 4.2% | 2.1% | …'; 'n = 2, orbit enclosing area | 43% | …' with the cited tests; results/RESULTS_PARADIGM.md

**Issue.** The table does not match its cited tests:
- The cited Lg test gives 2.047% at k_e = 20, which rounds to 2.0%.
- The n = 2 test runs only k_e = 20, 40 and 80, so the 43% is not produced by it.
- The label suggests that the enclosed area causes the larger error. It does not: the state-shift term U_even dominates.

**Evidence.** judge/h0_table.py:
- n = 1: 4.198, 2.047, 0.914 and 0.475%;
- n = 2, counter-clockwise: 42.74, 21.17, 10.40 and 5.14%;
- n = 2, reversed orbit: 35.35, 19.10, 9.83 and 4.97% (k_e × error → 3.98).
judge/h0_split.py: U_even alone gives 4.034 and U_odd alone 1.528.

**Fix.** Change '2.1%' to '2.0%' in PARADIGM.md and RESULTS_PARADIGM.md, and rename the second row 'n = 2, ellipse'.

STATUS: applied, together with the U_even attribution bullet. Still to do: after '- k_e = 10 is outside the averaging regime for the second instance.' add '(The k_e = 10 entry of the n = 2 row is computed with the same code but is not asserted; the test runs k_e = 20, 40 and 80.)'

#### [H0-period-averaging] minor (H8-H0-H6)

**Location.** §6c Remark H.0 '*Period averaging:* T₀·sup‖G∇V[μ]‖ ≪ ℓ and τλ_GT₀ ≪ ℓ²'

**Issue.** For σ > 0 the sup over ϑ ∈ ℝ^p of the confinement drift σ(ϑ − ϑ⁰) is +∞, so the condition can never hold as written. The confinement drift does not depend on time.

**Evidence.** G∇V = (1/k_e)G∇V_L + σ(ϑ − ϑ⁰). The first term is bounded by λ_G(F_ρ + K_Φ)sup‖∇_ϑΦ̃‖/k_e.

**Fix.** Replace with '*Period averaging:* T₀·sup_ϑ‖(1/k_e)G∇V_L[μ](ϑ)‖ ≪ ℓ, σT₀ ≪ 1 and τλ_GT₀ ≪ ℓ², where ℓ is the ϑ-scale of Φ̃. The sup is at most (T₀/k_e)λ_G(F_ρ + K_Φ)·sup‖∇_ϑΦ̃‖. The confinement drift does not depend on time and needs no averaging.'

STATUS: applied in the 02:14 file and correct.

#### [H0-tikhonov-vs-evidence] minor (H8-H0-H6)

**Location.** §6c Remark H.0 Tikhonov condition, and Evidence 'Limits of the evidence'

**Issue.** All Langevin runs use k_e = 10, about 3.5 times below the stated Tikhonov scale for their network. The limits list does not say so.

**Evidence.** For the network a_max = 4, s_min = 0.15 and f = sin 2x + 0.3x:
- F_ρ = 1.2469;
- (F_ρ + K_Φ)/s_min = 34.98;
- Lip(f) + K_Φ/s_min = 28.97.
The test_langevin instance uses this same network, and its frozen-μ error at k_e = 10 is 4.2%. With b_osc included, the first-order prediction is about 9% (0.898/k_e, G = diag(4,1,1); judge/h0_bosc.py).

**Fix.** Add to 'Limits of the evidence': '- All Langevin runs use k_e = 10, below the Tikhonov scale of Remark H.0 for this network ((F_ρ + K_Φ)/s_min ≈ 35). The measured frozen-μ drift error at k_e = 10 is 4.2% (n = 1), so the condition is conservative, but the agreement with (𝔄) is empirical.'

STATUS: applied in the 02:14 file and correct. Optionally append: 'Including the oscillation term of Remark H.7a, the first-order prediction is about 9%.'

#### [H2-fit-saddle] minor (H8-H0-H6)

**Location.** §6c Remark H.2 'The fit.': 'It suggests that the symmetric configuration is a *saddle* of the averaged landscape, not a barrier-protected minimum. That would explain why practice is far faster than the worst-case bound.'

**Issue.** The log(1/τ) fit is not evidence of saddle escape, for two reasons:
- It wins only on seed-averaged data.
- The saddle is so weak that a saddle-escape law would have a slope of 10³–10⁴ s, against the fitted 427 s.

**Evidence.** Regenerated results/langevin_tradeoff/table.md (01:48): per-seed RMS relative residuals are 0.274 (log law) and 0.330 (Arrhenius).

judge/saddle_exact_copy.py → saddle_out.txt, eigenvalues of G·Hess V at the Dirac fixed point:
- (𝔄): (−5.41e-5, 6.05e-3, 2.61e-2) s⁻¹;
- exactly averaged, k_e = 10: (−2.25e-4, 5.99e-3, 2.41e-2).
So 1/(2λ) is 9.2e3 s or 2.2e3 s, and over 1500 s the linear growth is at most e^{0.34}.

**Fix.** Replace the 'The fit.' bullet with:
'- **The fit.** The seed-averaged escape time fits t_esc ≈ −1693 + 427·log(1/τ) (4.5% relative residual) better than an Arrhenius law (25%). On the eight per-seed values both laws fit poorly (27% and 33%; `langevin_tradeoff/table.md`).
- **No saddle explanation.** The noiseless collapsed end point is a saddle of the averaged landscape, but a weak one: its unstable rate is about 5·10⁻⁵ s⁻¹ for (𝔄) and 2·10⁻⁴ s⁻¹ for the exactly averaged loop at k_e = 10. A saddle-escape law would therefore have a log(1/τ) slope of 10³–10⁴ s.
- **What the data show.** At the tested temperatures escape is governed by diffusion, by the approach transient and by the τ-dependent threshold. The data do not identify the escape mechanism.'

Also change §7 'saddle against Arrhenius fit' to 'log(1/τ) and Arrhenius fits, seed-averaged and per seed'.

STATUS: applied in the 02:14 file and correct. Optional wording fix: 'fit equally poorly (27% and 33%' → 'both fit poorly (27% and 33%'.

#### [H3-bakry-emery] minor (H8-H0-H6)

**Location.** §6c Remark H.3: 'The argument is Bakry–Émery with a Neumann boundary, then Holley–Stroock'

**Issue.** At σ = 0 the reference measure is uniform and its curvature is 0. The Bakry–Émery criterion needs curvature κ > 0, so it gives no LSI here. The 1/diam² LSI comes from a diameter-based theorem instead.

**Evidence.** CD(0,∞) alone implies no LSI. The diameter bound is what supplies the constant.

**Fix.** Replace the bullet with: '- The argument is a diameter-based log-Sobolev inequality for the uniform measure on a convex body, followed by Holley–Stroock. The diameter-based inequality comes from CD(0,∞) with a convex Neumann boundary plus a diameter bound (F.-Y. Wang 1997), or from Frieze–Kannan (1999) for log-concave measures on convex bodies. Together these give LSI(c·diam⁻²·e^{−2K_Φ(F_ρ+K_Φ)/(k_eτ)}).' Then add '- Bakry–Émery alone gives nothing at σ = 0.'

STATUS: applied in the 02:14 file and correct. The exact theorem numbers of the two references are unverified (see uncertain).

#### [H7-antisym-convention] cosmetic (H7-math, H7-code)

**Location.** §6c H.7(e) 'Its antisymmetric part is (1/k_e²)(A_even + 2K)'

**Issue.** The document uses S − Sᵀ, without the factor ½, as the 'antisymmetric part'.

**Evidence.** A_even = S_even − S_evenᵀ and 2K = K − Kᵀ. cyclic_defect also uses S − Sᵀ.

**Fix.** After '(e) *Non-reciprocity.*' insert: 'Write S(ϑ,ϑ′) := δP_μ(ϑ)/δμ(ϑ′) for the derivative in the direction δ_ϑ′, and call S(ϑ,ϑ′) − S(ϑ′,ϑ) its antisymmetric part (without the factor ½).'

STATUS: applied in the 02:14 file and correct.

#### [H7-definitions] cosmetic (H7-math)

**Location.** §6c H.7(b) 'is G∇_ϑP_μ(ϑ) − σ(ϑ − ϑ⁰), with P_μ(ϑ) = …'; (a) '∇r_μ(x_m) r_μ(x_m)'; design paragraph (L_τ) with g_t

**Issue.** Several objects in H.7 are implicit:
- P_μ is defined only through its expansion, so it is fixed only up to the gauge c(μ).
- ∇r must be the Jacobian. Read as the gradient, the n = 2 expansion fails.
- g_t is defined only in §1, with the unsaturated Φ.

**Evidence.** The code uses jac[t,k,l] = ∂r_k/∂x_l. A first variation is defined up to c(μ), which changes the antisymmetric part only by a coboundary.

**Fix.** Make four insertions:
- (b): '… is G∇_ϑP_μ(ϑ) − σ(ϑ − ϑ⁰), where P_μ(ϑ) := ⟨e_μᵀΦ̃(x_m + e_μ;ϑ)⟩_t satisfies'.
- (a): 'Here ∇r_μ is the Jacobian (∂r_k/∂x_l), and ∂_xΦ̃ = (∂Φ̃_k/∂x_l).'
- (e): 'Adding a function of ϑ alone, or a constant c(μ), to P_μ changes the antisymmetric part by a coboundary, which the criterion below allows.'
- (L_τ): 'g_t(ϑ) := e(t)ᵀΦ̃(x(t);ϑ)'.

STATUS: applied in the 02:14 file and correct.

#### [H7c-neglected-dependence] cosmetic (H7-math)

**Location.** §6c H.7(c) 'The neglected dependence of ρ^cl_μ on μ is O(1/k_e) relative.'

**Issue.** The sentence is true but understates the term. The neglected term has the same order as U_even, and it carries all of A_even.

**Evidence.** Analytic check: −(1/k_e)V_L + (1/k_e²)U_even = −(1/k_e)δJ^cl/δμ − (1/k_e²)∫(Φ̃_ϑ·∇)(½‖r_μ‖²)dρ + O(k_e⁻³), with J^cl(μ) := ½∫‖r_μ‖²dρ^cl_μ.

The second variation of N(ϑ) := −∫(Φ̃_ϑ·∇)(½‖r_μ‖²)dρ is ∫[rᵀ∂_xΦ̃_ϑ′Φ̃_ϑ + Φ̃_ϑ′ᵀJ_rΦ̃_ϑ]dρ. Its antisymmetric part is exactly A_even.

**Fix.** Replace with: '**The performative defect carries the even non-reciprocity.** With J^cl(μ) := ½∫‖r_μ‖²dρ^cl_μ,   −(1/k_e)V_L[μ] + (1/k_e²)U_even[μ] = −(1/k_e)δJ^cl/δμ − (1/k_e²)∫(Φ̃_ϑ·∇)(½‖r_μ‖²)dρ + O(k_e⁻³).   The first term is a true first variation, hence reciprocal. The second has the same order as U_even and carries all of A_even in (e).'

STATUS: applied in the 02:14 file and correct.

#### [S-range] cosmetic (H7-math)

**Location.** Companion files: '[S1]–[S9] point to the audit file'

**Issue.** [S10] onwards are now cited, so the stated range is out of date.

**Evidence.** §9 said [S1]–[S13]. [S14] has since been added.

**Fix.** '[S1]–[S14] point to the audit file' (and the same in §9).

STATUS: applied in the 02:14 file and correct.

#### [test-main-blocks] cosmetic (H7-math, H7-code)

**Location.** tests/test_averaging.py (fixed); tests/test_transport.py line 457 'if __name__ == "__main__":'

**Issue.** Script mode skips every test defined after the __main__ block.
- test_averaging.py is fixed: the block is now at the end.
- Judge's extension: tests/test_transport.py still has the block before test_corollary_B5_fails_along_lsigma_and_its_leak_variant_holds and test_theorem_F_leak_clause_is_tight_when_the_leak_pays_for_all_motion.

**Evidence.** grep: the __main__ block is at test_transport.py:457, and the two tests start at lines 464 and 499 (file of 520 lines, 01:49 version).

**Fix.** In tests/test_transport.py, move the block 'if __name__ == "__main__":\n    for name, fn in list(globals().items()):\n        if name.startswith("test_") and callable(fn):\n            fn(); print(name, "passed")' from before the '# ---- Corollary B.5 along (L) vs (L_σ)' section to the end of the file. No change is needed in test_averaging.py.

#### [H-evidence-W] cosmetic (H8-H0-H6)

**Location.** §6c Evidence, 'Equalisation': 'their final measures are 4.2× closer in W than without dither (0.40 against 1.69)'

**Issue.** The comparison is not what the sentence suggests:
- The W is the Euclidean W₂ over the physical (a, c, s), not W_{G⁻¹}.
- The two runs share their noise (seed 100 + seed), so the distance measures a synchronous coupling.
- 0.40 is at the sampling floor.

**Evidence.** exp_langevin_global.py: w2_phys and default_rng(100 + seed). judge/floor/out.txt: independent 32-point samples of μ*_τ are on average 0.670 apart (median 0.603; 10–90% range 0.465–0.987).

**Fix.** Replace with: 'their final 32-particle configurations, driven by the same noise realisation, are 4.2× closer in W₂ over the physical coordinates (a, c, s) than without dither (0.40 against 1.69); two independent 32-point samples of μ*_τ are typically about 0.6 apart, so 0.40 is at the sampling floor (review round 3);'

STATUS: applied in the 02:14 file and correct.

#### [H-route] cosmetic (H8-H0-H6)

**Location.** §6c Route to Conjecture H.6, bullets 1–3

**Issue.** Four points are imprecise:
- 'O(1/k_e)' is used for a lower bound.
- b is defined as the frozen part only, then b_osc is added.
- H.8 needs F_τ(μ_0) < ∞.
- The constant e^{1049} at the experiment's τ₀ is not stated.

**Evidence.** 2K_Φ(F_ρ + K_Φ)/τ₀ = 2·4·5.2469/0.04 = 1049.4. The whitened Gibbs measures depend only on (σ₀, τ₀), which gives τρ_τ = Θ(1/k_e).

**Fix.** Bullet 1: 'τρ_τ ≥ (σ₀/k_e)·e^{−2K_Φ(F_ρ+K_Φ)/τ₀}, so 1/(τρ_τ) = O(k_e). In this scaling the whitened proximal Gibbs measures do not depend on k_e, so τρ_τ = Θ(1/k_e). At the experiment's τ₀ = 0.04 the constant is e^{1049}.'

Bullet 2: 'Define b := (1/k_e²)[G∇(U_even + U_odd) + b_osc] + O(k_e⁻³).'

Bullet 3: append 'This needs F_τ(μ_0) < ∞, or a P₂ start, after the Girsanov regularisation of H(ii), which applies because b is bounded.'

STATUS: applied in the 02:14 file and correct.

#### [H6-hedge] cosmetic (H8-H0-H6)

**Location.** §6c Conjecture H.6, last bullet: 'By Proposition H.7 it is not a gradient flow.'

**Issue.** The bullet is stronger than H.7(e), which is hedged: it shows non-reciprocity in checked instances and excludes only the W_{G⁻¹} gradient-flow structure.

**Evidence.** judge/nonrecip_exp.py, on the experiments' network (a_max = 4, c_max = 3, s_min = 0.15, spread 32-atom measure):
- |cyclic sum|/max|S_even| over 20 random triples: median 0.87;
- K ≡ 0 (4e-18);
- exact k_e²·cyclic at k_e = 10: −1.8e-4, against −2.2e-4 predicted.

**Fix.** Replace with: '- By Proposition H.7(e) it is not the W_{G⁻¹}-gradient flow of any free energy. The cyclic sums of its antisymmetric second variation are non-zero in the checked instances, including the network of the experiments (review round 3).'

STATUS: applied in the 02:14 file and correct.

#### [H-notation] cosmetic (H8-H0-H6)

**Location.** §6c Design 'ϑ = (α, γ, κ)' against H(i) 'γ := N(ϑ⁰,(τ/σ)G)'; H(iv) 'α > α*' and 'F_{τ₀}', 'C₀ = σ/τ₀' against H.6 'τ := τ₀/k_e'

**Issue.** Three symbols are each used for two objects: γ, α and τ₀.

**Evidence.** H(iv) uses τ₀ for the initial annealing temperature. H.6 and the Route use τ₀ = k_eτ, and both meanings appear in exponents.

**Fix.** Rename:
- the Gaussian reference measure γ → γ_τ;
- the annealing constant α → α_a (also in Remark H.2);
- the initial annealing temperature τ₀ → τ_init, so that C₀ = σ/τ_init.
Keep τ₀ := k_eτ for H.6 and the Route.

STATUS: applied in the 02:14 file and correct.

#### [H1-reference-values] cosmetic (H8-H0-H6)

**Location.** §6c Remark H.1 'L̄(μ*_τ) = 0.111 (review round 2) … offline optimum … about 7·10⁻⁶'; Conjecture H.6 '0.1141 against … 0.1112 (review round 2)'

**Issue.** The reference values had no project provenance. 7.63e-6 is quoted as 'about 7·10⁻⁶', and it is a best-of-5 upper bound.

**Evidence.** Both results tables give 7.63e-06. judge/floor/out.txt: L̄(μ*_τ) = 0.111183. The new results/langevin_global/gibbs_reference.md (02:07) gives 0.11118, 0.00850, 0.11415, and L² offsets of 0.00681, 0.00143 and 0.00097. These agree with the round-3 reviewer's independent values.

**Fix.** Write 'about 7.6·10⁻⁶ (best of 5 L-BFGS restarts, an upper bound)'. Cite `exp_gibbs_reference.py` → `results/langevin_global/gibbs_reference.md` for 0.11118 and 0.11415, in place of 'review round 2'.

STATUS: applied in the 02:14 file and correct.

#### [H-misc] cosmetic (H8-H0-H6)

**Location.** Prop. H.8 'This is true for smooth bounded perturbations, but it is not verified here in general.'; H(ii) 'the estimates hold from t₀ on'; H.6 'satisfy the tube condition of Theorem D(b)'

**Issue.** Three side statements are imprecise:
- The H.8 sentence asserts the chain rule without proof.
- H(ii) does not say what replaces F_τ(μ_0) from t₀ on.
- D(b) is stated with a_max, but the saturated features need K_Φ.

**Evidence.** [S3] Lemma 2.2 needs integrated Fisher information. For the saturated features ‖Φ̃‖ ≤ K_Φ = a_max√n.

**Fix.** Three replacements:
- 'It is expected for smooth bounded perturbations (parabolic regularity and the entropy-dissipation estimate), but it is not verified here.'
- '… and the estimates hold with (μ_0, t) replaced by (μ_{t₀}, t − t₀).'
- '… the tube condition of Theorem D(b), with a_max replaced by K_Φ.'

STATUS: applied in the 02:14 file and correct.

#### [Cv-strict-case] minor (AF-regression)

**Location.** §4 proof of Theorem C(v), '*Strict case.* For a moving ideal whose optimal ideal–prior assignment changes, min_π sup_t b_lab(π)(t) > b.'

**Issue.** The claim is false in general. A switching optimal assignment does not force strict inequality.

**Evidence.** Analytic counterexample: m = 2, G = I, prior a ∈ {0, 1}, and a*(t) = (u, 1 − u) with u = −0.025 + 0.575 sin t.
- The costs are c_id = u² and c_swap = (1 − u)², and the swap is optimal when u > ½ (u reaches 0.55).
- b² = sup min(u², (1 − u)²) = 0.36, attained at u = −0.6.
- sup_t b_lab(id)² = 0.36, so min_π sup_t b_lab(π)(t) = 0.6 = b.

**Fix.** Replace the bullet with: '- *Strict case.* If the optimal ideal–prior assignment changes along the ideal, the inequality can be strict. For example, a*_1(t) = t and a*_2(t) = 1 − t on [0,1], with prior a ∈ {0, 1}, give b = 0.5 against 1.0.' Then add: '- *Not always strict.* For a*(t) = (u(t), 1 − u(t)) with u(t) = −0.025 + 0.575 sin t, the optimal assignment alternates between the identity and the swap, yet min_π sup_t b_lab(π)(t) = sup_t b_lab(id)(t) = 0.6 = b (review round 3).'

STATUS: applied in the 02:14 file and correct.

#### [E1-rate-not-bounded-by-PE] minor (AF-regression)

**Location.** §6 Remark E.1 'λ is bounded below by a function of the PE level αλ_min(G)/m, k_e, T and sup‖Ψ_G‖'; sentence after Theorem E 'The constants depend on m, G, k_e, T, α, sup‖Ψ‖, sup‖Ψ̇‖, …'; tests/test_transport.py E-test comment

**Issue.** Merges AF #2 and AF #8. No positive function of the PE level, k_e, T and sup‖Ψ_G‖ bounds λ below. The constants depend on Ψ(·) itself, which contradicts the uniformity that the document disclaims a few lines later.

**Evidence.** judge/e1_floquet.py: m = n = 1, ϑ* = (1.5, 0.3, 0.8), G = I, k_e = 3, x_m = 1.5 sin ωt, T = 2π. The PE level is 1.1329 and sup|Ψ| = 1.7791 for every ω, while the Floquet rate is:

| ω | 1 | 2 | 4 | 8 | 16 | 32 |
|---|---|---|---|---|---|---|
| rate | 4.66e-2 | 2.16e-2 | 6.83e-3 | 1.83e-3 | 4.66e-4 | 1.17e-4 |

This is the ω⁻² decay that averaging theory predicts, since Re H(jω) ~ k_e/ω².

**Fix.** Remark E.1: '… λ > 0 exists for each admissible Ψ, but it depends on Ψ(·), not only on its PE level. No positive function of αλ_min(G)/m, k_e, T and sup‖Ψ_G‖ alone bounds it below.' Follow with two sub-bullets: the ω-sweep example (rate from 4.7e-2 at ω = 1 to 2.9e-5 at ω = 64, like ω⁻²), and 'A lower bound must also involve sup‖Ψ̇_G‖ (step 3), and none is proved here.'

After Theorem E: 'k and λ depend on Ψ(·) itself, through the converse Lyapunov function of its linearisation, and on m, G and k_e. r₀ and c₀ depend in addition on a_max, s_min, ρ and δ. …'

STATUS: applied in PARADIGM.md (02:14) and in the test comment (01:49), and correct. The ω = 64 value follows the verified ω⁻² trend.

#### [E2-coupling-constants] minor (AF-regression)

**Location.** §6 Corollary E.2: 'Combine this with Theorem D(c) on [0, T₀] … Both solutions stay in the closed ρ-tube, which is compact.'

**Issue.** D(c) uses M_x = sup Lip_X f(·,t) over the whole union tube X. Theorem E constrains f only inside the moving tube, and under (S1) the constant is uniform only on bounded time intervals, which makes T₀ implicit. The result is true, because only the segments inside B(x_m(t), r) matter.

**Evidence.** Both solutions stay in D (proof step 4). The x-segment lies in B(x_m(t), ρ), where f = f̂*, so Lip(f̂*) and Lip(f̂) are each at most a_max J_c. The ϑ-segment lies in the convex inner set Θ_K^ε, where Proj_K is inactive.

**Fix.** Replace the two sentences with: 'Combine this with the coupling estimate of Theorem D(c) on [0, T₀], with T₀ := (Λ + λ)^{−1}log(CR₀/δ₀).' Follow with three bullets:
- 'Both solutions stay in D (proof step 4), so ‖x − x_m(t)‖ < r ≤ ρ.'
- 'The segment between x_t and x′_t lies in B(x_m(t), ρ), where f = f̂*, and Proj_K = id.'
- 'Hence (c) holds with M_x := k_e + 2a_max J_c, with L_x and L_ϑ computed for Ē = r and an inactive projection, and with L_2 = λ_G^{1/2}L₁. This Λ depends neither on T₀ nor on f outside the tube.'

STATUS: applied in the 02:14 file and correct.

#### [E-status] minor (AF-regression)

**Location.** §0 table, row E, status 'proved; linearisation checked; convergence observed'

**Issue.** The status overclaims. The only route to uniform exponential stability of the linearisation is [S8, Thm 1], whose hypotheses the proof admits were checked only against the abstract. The regularity class matters here, as the sin t² example and the ω-sweep show.

**Evidence.** §6 proof step 3: 'We could not check the exact richness definition and regularity class of [S8] against its full text.'

**Fix.** Status → 'proved, conditional on [S8, Thm 1] as stated in its abstract (its hypotheses were not checked against the full text); linearisation checked; convergence observed'.

STATUS: applied in the 02:14 file and correct.

#### [F-leak-test-power] minor (AF-regression)

**Location.** §0 row F 'proved and checked, including the leak clause and a tight instance'; §5b 'Checked at every sample … (ii′) along (L_σ)'; tests/test_transport.py::test_theorem_F_leak_clause_along_ot_sigma_mod

**Issue.** Severity raised from cosmetic. At review time the only leak-clause test could not detect an error in the leak term, so 'checked, including the leak clause' was not supported. It was also checked on a grid of pairs, not at every sample.

**Evidence.** judge/f_leak_power.py, same run as the test: worst ratio W/bound is 0.237 with the leak term and 0.308 with it deleted, so the test passes either way.

The new test test_theorem_F_leak_clause_is_tight_when_the_leak_pays_for_all_motion (01:49) passes. It has e ≡ 0, W ≈ 1.05 and ratio ≥ 0.99, and it fails if the leak term is dropped.

**Fix.** Code: already applied (the new tight test).

§5b: '- (ii′) along (L_σ), on a grid of pairs, including an instance with e ≡ 0 where it is tight (ratio ≥ 0.99) and the (L) reading fails.'

§0 row F: 'proved and checked, including the leak clause, a tight leak instance with e ≡ 0 and a tight velocity instance'.

§7 row F (ii′): also cite the new test.

STATUS: applied in the 02:14 file and correct.

#### [S5-tube-symbol] cosmetic (AF-regression)

**Location.** §1 (S5): 'take X to be the closure of ∪_t 𝒯(t) from Theorem D(b)'

**Issue.** 𝒯 is undefined in version 3, because D(b) was rewritten in terms of r(F).

**Evidence.** grep: 𝒯 occurs only on this line.

**Fix.** 'For bounded f, take X := ∪_t{x : ‖x − x_m(t)‖ ≤ r(F)} with F = sup‖f‖ (Theorem D(b)).'

STATUS: applied in the 02:14 file and correct.

#### [B5-variant-limit] cosmetic (AF-regression)

**Location.** §3 after the B.5 variant: 'As η̄ → 0 the variant tends to …'; proof 'replace b² by b² + 4η̄R'

**Issue.** The limit statement has two gaps, and the replacement is not the best available:
- The limit also needs T₂ → ∞.
- Along (L_σ) with b > 0, η̄ cannot tend to 0.
- The R-free replacement (b + 2η̄)² is sharper.

**Evidence.** judge/b5_w.py: W/b stays in [0.813, 0.884] on [20, 60] for b = 0.002, 0.01 and 0.05, and the RMS is linear in b (3.63e-4, 1.82e-3, 9.08e-3).

**Fix.** 'As η̄ → 0 and T₂ → ∞ the variant tends formally to …. Along (L_σ) with b > 0, the measured η̄ stays of order b (about 0.8–0.9b in the test instance), so the bound is first order in b.' In the proof: 'replace b² by (b + 2η̄)² (from √C ≤ b + 2η in C(ii)), or by b² + 4η̄R.'

STATUS: applied in the 02:14 file (it says 'about 0.9b') and correct.

#### [A-projection-row] cosmetic (AF-regression)

**Location.** §2 (iii) Projection row: 'In the example m = 1, n = 1, G = I, x = 0.7, e = 1 and ϑ = (0.93, 0.2, 1.5), with a in the layer: … −0.060'

**Issue.** The value −0.060 depends on (a_max, ε_K), and the example does not state them.

**Evidence.** judge/a_projection_row_copy.py, Jacobian entry by box:

| (a_max, ε_K) | (1, 0.15) | (1, 0.1) | (1, 0.3) | (0.95, 0.1) | (2, 0.15) |
|---|---|---|---|---|---|
| entry | −0.0598 | −0.0189 | −0.1236 | −0.1345 | 0 (not in the layer) |

h(0.93) = ((0.93 − 0.85)/0.15)² = 0.284.

**Fix.** Insert 'a_max = 1, ε_K = 0.15 (so h(0.93) = 0.284) and c_max = 2' into the example.

STATUS: applied in the 02:14 file and correct.

#### [A-only-if-activity] cosmetic (AF-regression)

**Location.** §2 proof of (i) 'on five active joint balls'; §7 row 'five active joint balls are each distinguished'

**Issue.** Only four of the five joint balls are active.

**Evidence.** In the test, (W_max, ε_W) = (3.0, 0.5) gives ρ = 2.5 ≥ √2·0.85. That is the second case of the construction: the joint ball is inactive and the per-kernel field is active.

**Fix.** 'on five joint balls. Four are active against an inactive per-kernel field, and one is inactive against an active per-kernel field.' Make the same change in §7.

STATUS: applied in the 02:14 file and correct.

#### [map-rows] cosmetic (AF-regression)

**Location.** §7 rows 'B.5 along (L) and (L_σ)', 'C(i), Lemma C.3', 'C(iv)'; list 'Not checked in the repository'

**Issue.** The verification map misdescribes four things:
- B.5 is not tested along (L).
- The C(i) test integrates the unprojected raw_velocity.
- At review time the C(iv) test asserted the version-2 gap condition; it now asserts version 3's.
- The 'not checked' list was incomplete.

**Evidence.** tests/test_transport.py: only the (L_σ) B.5 test exists; raw_velocity is used at line 259; the 01:49 version asserts 4η̄(b+η̄) < Δ_gap. No test covers B.6/B.7, the Filippov form of C(i), the moving-ideal part of C(v), Proposition C.5's D⁺V_OT formula, D(b), Corollary E.2 or Remark E.3.

**Fix.** Rows:
- 'B.5: failure of the plain bound along (L_σ); the (L_σ) variant holds';
- 'C(i), Lemma C.3 (unprojected pure leak)';
- C(iv) '…, including version 3's gap condition 4η̄(b + η̄) < Δ_gap'.

'Not checked' list: add 'B.5 along (L), and B.6 and B.7; the Filippov form of C(i), the inequality of C(v), and Proposition C.5's D⁺V_OT = σ_m(δ − δ²); D(b); Corollary E.2 and Remark E.3'.

STATUS: applied in the 02:14 file and correct.

#### [E-sint2-initial-state] cosmetic (AF-regression)

**Location.** §6 proof step 3: 'the ratio ‖w(t₀+200)‖/‖w(t₀)‖ is 0.9715, 0.9942 and 0.9980'

**Issue.** The ratios depend on the initial state, and the text does not give it.

**Evidence.** judge/e_sint2.py, ratios at t₀ = 5, 20, 50:
- from (0, 1): 0.9715, 0.9942, 0.9980;
- from (1, 0): 0.0935, 0.0136, 0.0075;
- from (1, 1)/√2: 0.753, 0.693, 0.711.

**Fix.** '… ζ̇ = Ψe, from w(t₀) = (e, ζ) = (0, 1), the ratio …'

STATUS: applied in the 02:14 file and correct.

#### [Db-strict-inequality] cosmetic (AF-regression)

**Location.** §5 proof of D(b): '‖e(t)‖ ≤ max(‖e(0)‖, (F + a_max)/k_e) < r(F)'

**Issue.** For e(0) = 0 the maximum equals r(F), so the strict inequality fails.

**Evidence.** r(F) = ‖e(0)‖ + (F + a_max)/k_e.

**Fix.** '‖e(t)‖ ≤ e^{−k_e t}‖e(0)‖ + (1 − e^{−k_e t})(F + a_max)/k_e < r(F) for every finite t, since a_max > 0.'

STATUS: applied in the 02:14 file and correct.

#### [Dc-L2-cost] cosmetic (AF-regression)

**Location.** §5 D(c) L_2 bullet: '≤ L_2 times the coupling's L²(G⁻¹) cost'

**Issue.** The bound holds with the square root of the transport cost, not with the cost itself.

**Evidence.** The two sides have different homogeneity.

**Fix.** 'For any coupling γ, ‖∫Φ(x;·)d(μ − μ′)‖ ≤ L_2(∫‖ϑ − ϑ′‖²_{G⁻¹}dγ)^{1/2}.'

STATUS: applied in the 02:14 file and correct.

#### [A-sigma-mod-row] cosmetic (AF-regression)

**Location.** §2 dictionary row 'classical σ-modification Γ(y − σθ) (leak inside the gain, `src/baselines/sigma_mod.py`)'

**Issue.** src/baselines/sigma_mod.py adapts only W (c_dot = 0, sigma_dot = 0). The D-weighted object is the one for a full-vector leak.

**Evidence.** src/baselines/sigma_mod.py lines 52–53 and 89–90.

**Fix.** Row label: '(leak inside the gain; `src/baselines/sigma_mod.py` implements its W-block, with c and σ frozen)'.

STATUS: applied, but the applied sentence order now attaches the λ-convexity to the W-block object. In the 02:14 row, replace 'For the W-block alone, (σ/2)∫(a − a⁰)ᵀ(I_n/m)(a − a⁰)dμ. In the G⁻¹ metric it is λ-convex with λ = σ·min(g_a/m, m g_c, m g_s), which depends on m.' with 'In the G⁻¹ metric the full-vector object is λ-convex with λ = σ·min(g_a/m, m g_c, m g_s), which depends on m. For the W-block alone the object is (σ/2)∫(a − a⁰)ᵀ(I_n/m)(a − a⁰)dμ, which is convex only in a (λ = σg_a/m on the a-block, 0 in (c, s)).'

#### [Ci-filippov-tau] cosmetic (AF-regression)

**Location.** §4 Theorem C(i): 'the bound weakens to D⁺Z ≤ −2σ_mZ + (σ_m/2)Σ_kθ_kA_{ττ_k} + (2Z)^{1/2}Γ'

**Issue.** τ is not specified. The proof holds for any optimal τ, so the minimum over optimal τ may be taken.

**Evidence.** The proof uses Z ≤ Z_τ for an arbitrary optimal τ.

**Fix.** '… ≤ −2σ_mZ + (σ_m/2)min_τ Σ_kθ_kA_{ττ_k} + (2Z)^{1/2}Γ, with the minimum over assignments τ optimal at t.'

STATUS: applied in the 02:14 file and correct.

#### [B4-sqrt-m] cosmetic (AF-regression)

**Location.** §3 Remark B.4: 'a set √m times smaller than the per-kernel product of Theorem A'

**Issue.** It is the radius that is √m times smaller, not the set.

**Evidence.** The joint radius is c_max. The ball circumscribing the per-kernel product has radius √m·c_max. I recomputed 137.8 (m = 6), 101.15 (m = 64) and 214.37.

**Fix.** '… used a joint centre ball of radius c_max. That radius is √m times smaller than the radius √m·c_max of the ball that circumscribes the per-kernel product of Theorem A.'

STATUS: applied in the 02:14 file and correct.

### Rejected

- (AF-regression) review/REVIEW_LOG.md has no round-2 section, so the header 'lists every confirmed issue and its fix' is false. The log parts of round-2 fixes AB.2 and F.1 were not applied. — The issue is already resolved. The current REVIEW_LOG.md (00:37:53, after REVIEW_round2.md at 23:41) has a full 'Round 2' section with verdict counts and rows for H.1–H.13, G.1–G.12, C.1–C.10, AB.1–AB.14, D.1–D.9, E.1–E.11 and F.1–F.9. The B.6 row carries the AB.2 erratum ('214 for every m … OT is not smaller globally'), and an 'Errata to round 1' list logs the F.1 contraction-identity erratum. The 'Pending: re-run round 2' line is gone; the only 'Pending' line is under the round-3 heading, which is still in progress. The header is accurate for rounds 1–2.
- (H7-code) The Remark H.0 gain-independence issue is major. — The issue itself is confirmed (item H0-gains), but at minor severity. Remark H.7a already states the correct fact ('not negligible as k_e → ∞ at fixed G'), the H.0 table is explicitly labelled frozen-μ, and no theorem or proof step depends on the sentence. It is an internal inconsistency in a remark labelled formal, not a stated result that is false.

### Uncertain

- Remark H.3 references. The mathematical point is certain: the Bakry–Émery criterion gives nothing at curvature 0, and a 1/diam² log-Sobolev inequality holds for log-concave measures on convex bodies. The citations are from memory, though: F.-Y. Wang 1997 (PTRF) and Frieze–Kannan, Ann. Appl. Probab. 9 (1999). I did not check the theorem numbers or exact constants.
- Theorem E, step 3. AF suggested an alternative self-contained proof: uniform complete observability via output injection, with (PE) plus the Ψ̇ bound implying the integrated observability Gramian bound. I did not verify it or adopt it. The conditional-status fix (E-status) does not need it.
- H.7(e) circulation proof, now in the document. I re-derived steps 1–4: Green's theorem on an affine integrand gives −(ε²/2)Σ_cyc at the centroid, and Σ_j|dw_j| = 2ε per edge gives the 6ε bound. I did not re-run the claim that 'the circulation equals the centroid formula to all printed digits' (ε = 0.05, n = 2). For quadratic U it is an exact identity, so equality to rounding is expected.
- Concurrent edits. PARADIGM.md, LITERATURE_TRANSPORT.md, RESULTS_PARADIGM.md, the tests and the experiments were being edited during judging. My STATUS notes and the residual-error list reflect the files as of 02:14:12. I did not review edits after that.
- Two G-block cosmetic numbers were not reproduced. The refuter reported (b′) matching the joint flow to within 1.5%, 0.6% and 0.2% at S = 2, 8.5 and 17. No fix text depends on these numbers.

## Refuter: H7-math

### Holds

- Scope: reviewed the current PARADIGM.md (saved 01:00:39, which already restricts H.7(e) to W_{G⁻¹} and adds Remark H.7a), the current averaging.py and test_averaging.py (6/6 pass), and LITERATURE_TRANSPORT.md (01:00:16, now defines [S10]–[S13]). Scripts are in <scratchpad>\round3\H7_math\.
- (a) The expansion e_μ = r/k_e + (∇r r − ṙ)/k_e² + O(k_e⁻³) is correct when ∇r is the Jacobian. It is uniform over μ ∈ P(ℝ^p): ‖f̂_μ‖ ≤ K_Φ and every x-derivative of Φ̃ is bounded by the saturation, e.g. ‖∂_xΦ̃‖ ≤ a_max√n·e^{−1/2}/s_min. The hypotheses f, x_m ∈ C² are enough for the O(k_e⁻³) remainder.
- (b) The expansion was checked against the exact periodic solution (reproduce.py). k_e·(k_e²(P_μ − P₀) − U) converges: for n = 1 to (−0.344, −0.252, 0.484, −0.383) at k_e = 320, and likewise for n = 2 in both orientations. The reported relative errors reproduce: 0.537 → 0.247 (n = 1), 0.038 → 0.019 (n = 2, +) and 0.023 → 0.012 (n = 2, −). Extending to k_e = 160 and 320 gives 0.118 and 0.058 for n = 1. The identity (∇r r)ᵀΦ̃ + rᵀ(∂_xΦ̃)r = (r·∇)(rᵀΦ̃) holds.
- Chain rule through the saturations holds. The period average of langevin.SaturatedLiftedNetwork.grad_g along (x_m + e_μ, e_μ) equals ∇_ϑP_μ by finite differences to 1.5e-11 (n = 1) and 2.5e-11 (n = 2), with |a|/a_max up to 0.75 and |c|/c_max up to 0.94 (chainrule.py). So the averaged drift of (L_τ) is G∇_ϑP_μ − σ(ϑ − ϑ⁰).
- (c) The semi-gradient identity holds to O(k_e⁻³): k_e³·(lhs − rhs) converges, n = 1 ≈ (0.021, 0.041, −0.024, 0.028) and n = 2 ≈ (−0.097, 0.026, −0.088, 0.043). The neglected dependence of ρ^cl_μ on μ is 0.92/k_e (n = 1) and 1.0/k_e (n = 2) relative, as stated (secondvar.py).
- (d) K is antisymmetric (|K + Kᵀ| = 4e-16). The Stokes/Green formula, with the counterclockwise orientation of (1.2 cos t, 0.8 sin t), matches the loop integral to 7e-15 relative, far better than the stated 1e-6. K ≡ 0 for n = 1 because 1-forms on ℝ are exact. Reversing the orbit keeps ρ, leaves U_even unchanged and flips U_odd.
- (e) The S_even formula equals finite differences of U_even in μ to 4e-12, and K equals the μ-variation of U_odd to 2e-12, so A_even = S_even − S_evenᵀ is correct. The cyclic sums reproduce the table: n = 1, −0.0704 predicted against −0.0711 exact (k_e = 80); n = 2, −0.0894 + 0.0430 = −0.0464 predicted against −0.0471 exact. They are non-zero for all 8 seeds tried, for n = 1 and n = 2 (secondvar.py).
- The coboundary criterion is the right gauge-invariant test. A first variation is defined up to a constant c(μ), and a second variation up to a function a(ϑ) of ϑ alone. Each changes the antisymmetric part only by a coboundary, which the criterion already allows. The cyclic-sum characterisation of coboundaries is correct.
- The current conclusion, 'not the W_{G⁻¹}-gradient flow of any free energy', is justified. For τ > 0, the diffusion τ∇·(G∇μ) together with an entropy term forces the mobility G among Otto-type metrics with local mobility. So the hedge 'a different, μ-dependent metric is not excluded' is accurate.
- Non-reciprocity persists at the experiments' k_e = 10 (ke10.py). Exact k_e²·cyclic sum: −0.0618 against −0.0704 predicted for the test instance, and 0.0127 against 0.0092 for a 32-atom spread measure on exp_langevin_global's network. The exact/predicted ratio is 1.008–1.078 for k_e = 20–80. So the Evidence sentence 'already non-reciprocal for n = 1, through the state shift' is supported at k_e = 10.
- Remark H.7a's mathematics holds. The second-order averaging formula is correct: the −⟨Du·⟨H⟩⟩ term vanishes because u has zero mean. b_osc is odd under reversal, since u_rev(t) = −u(−t). The 'non-reciprocity survives' argument is valid, because the average of two W_{G⁻¹}-gradient flows is again one.
- Round-2 fixes H.1 (Remark H.0 marked formal; an n = 2 averaging test) and H.2 (H.6 restated; the fixed-k_e counterexample; H.7 and H.8) are applied. Version 3 also corrects round 2's claim that the non-gradient term vanishes for n = 1, since A_even ≠ 0 there.

### Issues

#### minor: §6c Proposition H.7(a): 'has a unique T₀-periodic solution e_μ. It attracts every solution exponentially'; proof of (a): 'Exponential attraction follows from the linearisation −k_e + Lip(r_μ).'

**Claim.** For autonomous f ∈ C² and k_e ≥ k₀, e_μ is the only T₀-periodic solution and attracts every solution exponentially.

**Argument.** The hypotheses allow f(x) = x². Take n = 1, x_m = sin t, μ = δ at α = 0 (so a = 0 and f̂ = 0), and k_e = 10, giving ė = (sin t + e)² − 10e. Besides the attracting periodic solution (e ∈ [0.0009, 0.1248]) there is a repelling one with e ∈ [7.891, 11.910], found by backward integration and periodic to 3e-12. The solution from e(0) = 12 blows up: |e| reaches 1e6 at t = 0.176 (<scratchpad>\round3\H7_math\part_a_global.py and part_a_global2.py). For any k_e, every e(0) > k_e + 2 blows up, and by the intermediate value theorem there is a second periodic orbit, so no k₀ works. The contraction and the linearisation in the proof are local only.

**Fix.** Replace the two sentences with: 'The frozen-μ error equation ė = r_μ(x_m + e) − k_e e has a unique T₀-periodic solution e_μ with sup_t‖e_μ‖ ≤ 1. The ball {‖e‖ ≤ 1} is forward invariant, and e_μ attracts every solution starting in it at exponential rate at least k_e − L₁, where L₁ := Lip(f on the 1-neighbourhood of 𝒪) + a_max√n·e^{−1/2}/s_min. If f is also one-sided Lipschitz, ⟨f(x) − f(y), x − y⟩ ≤ L_f‖x − y‖², then for k_e > L_f + a_max√n·e^{−1/2}/s_min the periodic solution is unique and attracts every solution.' Replace the proof bullet with: 'On ‖e‖ = 1, d/dt ½‖e‖² ≤ sup_{dist(x,𝒪)≤1}‖f‖ + K_Φ − k_e < 0 for k_e ≥ k₀. In the ball, d/dt ½‖e − e′‖² ≤ −(k_e − L₁)‖e − e′‖². The global statement uses L_f in place of the local Lipschitz constant of f.'

#### minor: §6c Proposition H.7(e): 'Its antisymmetric part is (1/k_e²)(A_even + 2K) + O(k_e⁻³)'; proof of (e): 'Differentiate (b) in μ.'

**Claim.** The O(k_e⁻³) remainder of (b) may be differentiated in μ, so the antisymmetric part of the exact second variation is (1/k_e²)(A_even + 2K) + O(k_e⁻³).

**Argument.** (b) bounds the remainder uniformly in ϑ and μ and in its ϑ-gradients, but not its first variation in μ. A remainder that is small in sup norm can have a large μ-derivative, so neither the '+ O(k_e⁻³)' nor the step 'Differentiate (b) in μ' follows from what is stated. The conclusion itself is true, and it can be proved with no μ-derivative of the remainder. U_even + U_odd is quadratic in μ, so its circulation around a small triangle of probability measures equals exactly −(ε²/2) times the cyclic sum at the centroid: 5.952830e-05 against 5.952830e-05 for ε = 0.05, n = 2 (loopcheck.py in the scratch directory above). A first variation plus a constant c(μ) has zero circulation, and the remainder contributes at most 6ε·sup|R|. The exact cyclic sums also converge at the rate claimed: k_e·(exact − predicted) = −0.011, −0.045, −0.058 at k_e = 20, 40, 80 for n = 1 (ke10.py).

**Fix.** In (e), replace 'Its antisymmetric part is (1/k_e²)(A_even + 2K) + O(k_e⁻³), with' by 'The second variation of its k_e⁻² coefficient U_even + U_odd has antisymmetric part A_even + 2K, with'. Replace 'Differentiate (b) in μ.' by: 'Fix ϑ₀, ϑ₁, ϑ₂ and ε ∈ (0,1). For (a,b) in T_ε := {a, b ≥ 0, a + b ≤ ε}, let μ_{a,b} := (1−ε)μ + (ε−a−b)δ_{ϑ₀} + aδ_{ϑ₁} + bδ_{ϑ₂}. The circulation ∮_{∂T_ε}⟨Q_{μ_{a,b}}, dμ_{a,b}⟩ vanishes when Q_μ is a first variation plus a constant c(μ), since dμ_{a,b} has zero mass; in particular it vanishes for −(1/k_e)V_L. For U_even + U_odd, which is quadratic in μ, it equals −(ε²/2)·Σ_cyc(A_even + 2K)(ϑ_i, ϑ_{i+1}) at the centroid of T_ε. The remainder of (b) contributes at most 6ε·sup|R| = O(εk_e⁻³). So if P_μ were a first variation up to O(k_e⁻³), the cyclic sum would be O(1/(εk_e)) for every ε, hence zero.' Keep the formulas for S_even and K and the coboundary lemma.

#### minor: §0 contributions table, row H.7, status: 'expansion, Stokes identity, parity and non-reciprocity checked to about 1%'

**Claim.** The expansion of P_μ is checked to about 1%.

**Argument.** The H.7 'Checked' table itself gives a relative error of k_e²(P_μ − P₀) against U_even + U_odd of 0.54 and 0.25 at k_e = 40 and 80 for n = 1. I reproduce 0.537 and 0.247, and get 0.118 and 0.058 at k_e = 160 and 320 (reproduce.py). For n = 2 it is 0.038/0.019 (+ orientation) and 0.023/0.012 (− orientation). The expansion is correct: k_e³ times the remainder converges to (−0.344, −0.252, 0.484, −0.383), against max|U_even| = 0.026. But it is checked through its rate, not to 1%. Stokes is checked to 7e-15, parity to 1e-10, and non-reciprocity to about 1%.

**Fix.** Replace with: 'frozen-measure part proved; the expansion is checked through its convergence rate (relative error of the k_e⁻² coefficient 25% for n = 1 and 1.9% for n = 2 at k_e = 80, halving with each doubling of k_e); the Stokes identity, parity and non-reciprocity are checked to 1e-6, 1e-10 and about 1%.'

#### minor: §8 item 2: '*Proved:* at finite gain, the averaged adaptive dynamics is a semi-gradient flow. Its antisymmetric interaction is generated by the area that the reference orbit encloses.'; §0: 'Its kernel is the area enclosed by the reference orbit (Stokes).'

**Claim.** The antisymmetric interaction comes from the area the orbit encloses, and the kernel is that area.

**Argument.** By H.7(e), the antisymmetric part of the interaction is A_even + 2K. A_even needs no enclosed area and is the only source for n = 1 (cyclic sum −0.0704 with K ≡ 0). In the n = 2 instance it is also the larger part (−0.089 against +0.043). K is not the area either. For Gaussian features, Σ_k det[∇Φ̃_k(ϑ), ∇Φ̃_k(ϑ′)] = (aᵀa′)·det[∇φ, ∇φ′], so K = 0 whenever aᵀa′ = 0 or c = c′, whatever the enclosed area. On the same ellipse, K = 0.0 for a ⟂ a′ and 6.1e-3 for a ∥ a′ (loopcheck.py). 'Proved' also overstates: the semi-gradient form is proved only at a frozen measure, Remark H.7a adds a further term along the real loop, and non-reciprocity is checked in two instances, not proved.

**Fix.** §8: '*Proved at a frozen measure:* the averaged adaptive dynamics is a semi-gradient flow plus a lag term (Proposition H.7). Its interaction has two antisymmetric parts. The first is A_even, from the state shift; it is present for every n and is the only one for n = 1. The second is the circulation kernel K, an integral over the region that the reference orbit encloses (n ≥ 2). Non-reciprocity is checked in two instances.' §0: 'Its kernel is an integral of Σ_k det[∇Φ̃_k(ϑ), ∇Φ̃_k(ϑ′)] over the region that the reference orbit encloses (Stokes).'

#### minor: §0: 'ignoring that this distribution depends on the parameters, as temporal-difference learning does'; §6c H.7(c): 'This is a semi-gradient, as in temporal-difference learning [S10], [S11].'; LITERATURE_TRANSPORT.md [S11]: 'The semi-gradient (data distribution held fixed) structure of TD.'

**Claim.** TD learning is a semi-gradient because it ignores the dependence of the data distribution on the parameters.

**Argument.** In [S11], states are sampled from the stationary law of a fixed Markov chain (policy evaluation), which does not depend on the parameters. TD is a semi-gradient because it does not differentiate the bootstrapped target γV_θ(s′) (Sutton–Barto, 2nd ed., §9.3). [S10] also samples from a fixed distribution. What H.7(c) leaves out is the dependence of the data distribution ρ^cl_μ on μ. That is the structure of performative prediction (Perdomo, Zrnic, Mendler-Dünner, Hardt, ICML 2020), where gradient steps taken with the distribution held at its current value converge to performatively stable points, not optimal ones. The analogy with TD's non-symmetric matrix Φᵀ D(γP − I)Φ fits the lag term better: K correlates features at time t with features at the filter-delayed time.

**Fix.** H.7(c): 'This is a semi-gradient: the update freezes part of its own parameter dependence, as TD learning does [S10], [S11]. In TD the frozen part is the bootstrapped target, and the data distribution does not depend on the parameters. Here the frozen part is the data distribution ρ^cl_μ, as in performative prediction (Perdomo et al., ICML 2020).' §0: replace 'as temporal-difference learning does' with 'as in performative prediction; like semi-gradient TD, the resulting update is not the gradient of any objective'. [S11] annotation: 'The semi-gradient (bootstrapped target held fixed) structure of TD; the data distribution is that of a fixed policy.' Add the Perdomo et al. entry to the audit file.

#### minor: §6c Remark H.0: 'The first corrections are O(1/k_e) relative to the drift and do not depend on the gains. Proposition H.7 gives them exactly at frozen μ.'

**Claim.** The first-order corrections to (𝔄) do not depend on the gains.

**Argument.** This is true only of the frozen-measure corrections. The oscillation term of Remark H.7a, b_osc = ⟨D_ΘH·u⟩, is quadratic in G (H ∝ G and u ∝ G·T₀), and the document itself says its size relative to the H.7 terms grows like |G|T₀. So in the formal limit stated in Remark H.0 (k_e → ∞ with G fixed) it has the same order and depends on G. In the orbit-orientation runs (G = diag(4,4,1,1,1)), the H.7 terms alone overpredict by 2.29–2.67×, while the full theory matches to 3%.

**Fix.** Replace with: 'The frozen-measure corrections are O(1/k_e) relative to the drift and do not depend on the gains; Proposition H.7 gives them. Along the real loop, the oscillation of the parameters within each period adds a correction of the same order, whose size relative to them grows like |G|T₀ (Remark H.7a). The table measures the frozen-measure part:'

#### cosmetic: §6c H.7(e): 'Its antisymmetric part is (1/k_e²)(A_even + 2K)'; 'Checked' row: 'k_e²·cyclic sum of the antisymmetric part'

**Claim.** The antisymmetric part of the second variation is (1/k_e²)(A_even + 2K).

**Argument.** A_even = S_even − S_evenᵀ and 2K = K − Kᵀ, so 'antisymmetric part' here means S − Sᵀ, as in averaging.cyclic_defect. The usual meaning is ½(S − Sᵀ). Under that convention the coefficient is ½A_even + K and the reported cyclic sums halve: −0.0352 and −0.0232 for the predictions (secondvar.py).

**Fix.** After '(e) *Non-reciprocity.*' insert: 'Write S(ϑ,ϑ′) := δP_μ(ϑ)/δμ(ϑ′) for the derivative in the direction δ_ϑ′, and call S(ϑ,ϑ′) − S(ϑ′,ϑ) its antisymmetric part.'

#### cosmetic: §6c H.7(b): 'is G∇_ϑP_μ(ϑ) − σ(ϑ − ϑ⁰), with P_μ(ϑ) = …'; (a): '∇r_μ(x_m) r_μ(x_m)'; proof of (b): '(∇r r)ᵀΦ̃ + rᵀ(∂_xΦ̃)r = (r·∇)(rᵀΦ̃)'; design paragraph: '(L_τ)' with g_t

**Claim.** The objects in H.7 are defined implicitly.

**Argument.** (1) P_μ is fixed only through its ϑ-gradient, so it is defined up to a constant c(μ), and 'the second variation of P_μ' is then ambiguous. The criterion is invariant under this gauge and under adding a(ϑ), but the text should say so. 'Symmetric at order 1/k_e' holds for directions δ_ϑ′; for directions δ_ϑ′ − μ it holds only up to a coboundary. (2) The identity needs ∇r to be the Jacobian (∂r_k/∂x_l). Read as the gradient, ∇r r = ∇(½‖r‖²) and the identity fails for a non-symmetric Jacobian: the n = 2 expansion error becomes 1.56 instead of 0.019 at k_e = 80 (extras.py). (3) g_t is defined only in §1, with the unsaturated Φ.

**Fix.** In (b): '… is G∇_ϑP_μ(ϑ) − σ(ϑ − ϑ⁰), where P_μ(ϑ) := ⟨e_μᵀΦ̃(x_m + e_μ;ϑ)⟩_t satisfies'. In (a): 'Here ∇r_μ is the Jacobian (∂r_k/∂x_l), and ∂_xΦ̃ = (∂Φ̃_k/∂x_l).' In (e): 'Second variations are taken in the directions δ_ϑ′. Adding a function of ϑ alone, or a constant c(μ) to P_μ, changes the antisymmetric part by a coboundary, which the criterion allows.' In the §6c design paragraph, after (L_τ): 'with g_t(ϑ) := e(t)ᵀΦ̃(x(t);ϑ).'

#### cosmetic: §6c H.7(c): 'The neglected dependence of ρ^cl_μ on μ is O(1/k_e) relative.'

**Claim.** The neglected dependence is O(1/k_e) relative.

**Argument.** The statement is true: 0.92/k_e for n = 1 and 1.0/k_e for n = 2, from finite differences of J^cl(μ) = ½∫‖r_μ‖²dρ^cl_μ (secondvar.py). But the neglected term has the same order as U_even, and it is exactly the even non-reciprocity. It equals −(1/k_e²)∫(Φ̃_ϑ·∇)(½‖r_μ‖²)dρ + O(k_e⁻³): for n = 2 at k_e = 80, k_e·(full − semi) = (−0.1064, 0.0640, 0.0116, −0.1989) against (−0.1064, 0.0639, 0.0114, −0.1961) from the formula (extras.py). The antisymmetric part of the second variation of this term is exactly A_even. As written, the sentence suggests the omission is negligible at the order H.7 resolves.

**Fix.** Replace with: 'The neglected dependence has the same order as U_even: −(1/k_e)V_L + (1/k_e²)U_even = −(1/k_e)δJ^cl/δμ − (1/k_e²)∫(Φ̃_ϑ·∇)(½‖r_μ‖²)dρ + O(k_e⁻³), with J^cl(μ) := ½∫‖r_μ‖²dρ^cl_μ. The last term carries all of A_even.'

#### cosmetic: PARADIGM.md 'Companion files': 'Reference numbers [1]–[89] point to review/LITERATURE_refs.md, and [S1]–[S9] point to the audit file.'

**Claim.** The audit file holds [S1]–[S9].

**Argument.** [S10]–[S13] are now cited (H.7(c), the §0 table, §8, Corollary B.6) and are defined in LITERATURE_TRANSPORT.md. §9 already says [S1]–[S13].

**Fix.** Change '[S1]–[S9] point to the audit file' to '[S1]–[S13] point to the audit file'.

#### cosmetic: tests/test_averaging.py: test_oscillation_correction_is_odd_under_orbit_reversal is defined after the 'if __name__ == "__main__":' block

**Claim.** The file's script mode runs every test.

**Argument.** pytest collects the test (6 passed). But 'python tests/test_averaging.py' runs the main block before this function is defined, so the H.7a test is skipped without any message.

**Fix.** Move the 'if __name__ == "__main__":' block to the end of tests/test_averaging.py.

## Refuter: H7-code

### Holds

- Versions reviewed (they changed during the review): PARADIGM.md 01:00:39, src/core/averaging.py 00:51:36 (adds oscillation_correction), tests/test_averaging.py 00:59:07 (6 tests), tests/test_langevin.py 23:44. Under pytest, all 6 averaging tests and Lg::test_averaging_error_of_the_formal_reduction_decreases_like_one_over_ke pass. Nothing was written under the project root. Scripts: <scratchpad>/round3/averaging-impl/ (common.py, indep.py = independent Fourier-collocation implementation).
- periodic_error returns the periodic solution. Against an independent Fourier-collocation solve (N = 2000), the RK4 output (N = 2000, sub = 2) has max|Δe| = 4.9e-11/1.0e-10 (n=1) and 5.5e-12/1.1e-11 (n=2) at k_e = 40/80, and k_e²|ΔP| ≤ 1.8e-8 against measured quantities of 0.02–0.2. The output samples are time-aligned with orbit.xm. The transient is irrelevant after one period (check_misc.py).
- averaged_potential computes ⟨e_μᵀΦ̃(x_m+e_μ;ϑ)⟩ and leading_potential computes −V_L/k_e. Both use the same quadrature, so k_e²(P−P₀) carries no quadrature mismatch.
- expansion_terms implements H.7(b) exactly, with the convention [t,j,k,l] = ∂_lΦ̃_k, the correct signs and (∇r r)_k = Σ_l ∂_l r_k r_l. An independent implementation (finite differences in x along r, spectral d/dt) agrees to ≤ 2e-12. Richardson extrapolation of the exact k_e²(P−P₀) over k_e = 40…1280 (collocation) reproduces U_even+U_odd to 5.9e-8 (n=1) and 4.4e-10 / 8.0e-11 (n=2, both orientations) (check_expansion.py). So the large n = 1 test error (0.54 → 0.25) is a large k_e⁻³ coefficient, not a formula error.
- H.7(a) expansion (locally): k_e³·max_t|e_μ − r/k_e − (∇r r − ṙ)/k_e²| stays bounded (→ 26.0, 3.28, 1.71 for the three test cases, k_e = 40…640). H.7(c) semi-gradient form converges at O(1/k_e) (n=1: 0.029 → 0.0027 for k_e = 40 → 640) (check_secondvar.py).
- H.7(d) holds and is implemented correctly:
- U_odd = −⟨ḟᵀΦ̃⟩ + ∫K dμ to 6e-16.
- circulation_kernel equals a spectral-d/dt version to 1e-16.
- The Green test integrates the right 2-form Σ_k dΦ̃_k(ϑ)∧dΦ̃_k(ϑ') = det[∇Φ̃_k(ϑ),∇Φ̃_k(ϑ')]dA, with the right (counter-clockwise) orientation and Jacobian 0.96ρ. Agreement is 6.8e-15.
- |K+Kᵀ| = 4e-16, K ≡ 0 for n = 1 (9e-16), and orbit-reversal parity holds to 2e-15.
- second_variation_formula is right as a full matrix, not only in its antisymmetric part:
- S_even + K equals the Richardson limit of k_e²S_exact + k_eS₀ (k_e = 80…1280) to 8.6e-9, 7.1e-9 and 4.2e-9 (n=1, n=2 ±).
- S_even equals the weight-derivative of an independent U_even to 2e-10.
- The leading second variation is −(1/k_e)⟨Φ̃_ϑᵀΦ̃_ϑ'⟩ as stated.
- exact_second_variation (central differences in the weight of an added signed atom) is legitimate. The change of total mass adds a function of ϑ only (a coboundary), which cancels in cyclic sums. With h = 1e-3 the truncation and round-off errors are negligible. An independent collocation computation at k_e = 80 gives the same cyclic sums (−0.071066, −0.047056).
- All quoted numbers reproduce:
- H.7 table: 0.5372→0.2468 (n=1), 0.0382→0.0188 (n=2, +); −0.07107 vs −0.07035; −0.04706 vs −0.04637; shift −0.08935, circulation +0.04298.
- Remark H.0 rows: 4.198/2.047/0.914/0.475 % and 42.7/21.2/10.4/5.14 %, with k_e×err 0.37–0.42 and 4.11–4.27.
- The Euler stepping of the n=1 Lg test agrees with RK4 (4.199/2.054/0.919/0.472 %) (repro_tables.py, repro_H0.py).
- H.7 predicts the Remark H.0 constants: k_e×error 0.416 for the n=1 instance (measured ≈ 0.4) and 4.05 for the n=2 instance (measured 4.11–4.27) (predict_H0_n1.py, split_H0.py).
- The expansion and Stokes tests discriminate. Running the real test functions on mutated code, all 14 formula mutations fail them: U_odd sign, (∇r)ᵀr in place of ∇r r, dropped or halved rᵀ∂Φ̃r, dropped ∇r term, transposed Jacobian in ṙ, ±5/10/20 % scalings, and transposed or sign-flipped ∂_xΦ̃. The non-reciprocity test fails for Kᵀ, 2K, K/2, K = 0, a transposed (Φ̃_j·∇) term and ±8 % S_even (mutate.py, check_misc.py).
- The new oscillation_correction is correct. An independent near-identity Poincaré-map construction (cumulative-trapezoid antiderivative, RK4 one-period flow of Θ̇ = H/k_e, Ψ(Θ) = Θ − u(Θ,0)/k_e) converges to it: 3.3 % → 0.38 % (n=2) and 12 % → 1.9 % (n=1) for k_e = 20 → 160. It is exactly odd under reversal and scales as G² (×4.000 for 2G). So the H.7a 'Size' and 'Parity' bullets hold, and the 'Non-reciprocity survives' argument is valid (check_osc.py).
- The round-2 fixes in this block were applied:
- Remark H.0 is labelled formal.
- The H(i) row says 'does not test averaging'.
- New 1/k_e tests exist for n = 1 and for an area-enclosing n = 2 orbit.
- v3 correctly overturns round 2's statement that the non-gradient part vanishes for n = 1: the extrapolated cyclic sum is −0.0704 ≠ 0, all from A_even.

### Issues

#### major: §6c Remark H.0, bullet 'The first corrections': 'The first corrections are O(1/k_e) relative to the drift and do not depend on the gains. Proposition H.7 gives them exactly at frozen μ.' (and the table after it, presented as the measured error of the reduction)

**Claim.** The O(1/k_e) correction of the closed loop relative to (𝔄) is gain-independent and is given by Proposition H.7.

**Argument.** Along the real loop the parameters move within each period. Second-order averaging then adds b_osc = ⟨D_ΘH·u⟩/k_e², which is of the same order in 1/k_e and ∝ G² (∝ |G|T₀ relative to the frozen terms).

The document's own Remark H.7a says so: 'It is not negligible as k_e → ∞ at fixed G'. Remark H.0 therefore contradicts H.7a, and in H.0's own 'exact formal limit' (k_e → ∞ with G fixed) b_osc is a first-order term.

On the very instances of the H.0 table, at their gains:
- k_e × relative size of b_osc is 0.92 against 0.42 for the frozen terms (n=1, G=diag(4,1,1)), and 4.44 against 4.05 (n=2, G=I).
- k_e × the full first-order relative error is 0.90 and 4.13. For n = 1 this is twice the quoted ≈ 0.4.
- Replacing G by 2G multiplies b_osc by 4.000.

The repo's new experiment confirms it: the frozen terms account for less than half of the orientation-odd response (size ratio 2.3–2.7).

Scripts: <scratchpad>/round3/averaging-impl/check_osc.py and H0_full.py.

**Fix.** Replace the bullet by: '- **The first corrections** are O(1/k_e) relative to the drift. At frozen μ they do not depend on the gains, and Proposition H.7 gives them exactly. Along the real loop the parameters also move within each period, and second-order averaging adds the oscillation term b_osc of Remark H.7a. It is of the same order in 1/k_e, and relative to the frozen terms it scales like |G|T₀. Measured relative errors of the frozen-μ averaged drift against (𝔄), b_osc excluded: [table unchanged]. At the gains of these instances b_osc is not small: k_e × its relative size is 0.92 (n = 1, G = diag(4,1,1)) and 4.4 (n = 2, G = I), against 0.42 and 4.05 for the frozen terms; k_e × the full first-order relative error is 0.90 and 4.1.'

Append to the 'exact formal limit' bullet: 'In this limit both the frozen terms and b_osc are first-order corrections.'

Also correct the docstring of experiments/exp_orbit_orientation.py: 'Θ₊ − Θ₋ ≈ −2 J⁻¹ b_odd' → 'Θ₊ − Θ₋ ≈ −2 J⁻¹ (b_odd + b_osc)'.

#### minor: §6c Proposition H.7(a): 'has a unique T₀-periodic solution e_μ. It attracts every solution exponentially'. Also the proof bullet 'Exponential attraction follows from the linearisation −k_e + Lip(r_μ)', and the averaging.py docstrings 'k_e large it has a unique periodic solution e_μ' and 'the transient decays like exp(−k_e T₀)'.

**Claim.** Under 'f ∈ C²(ℝⁿ) autonomous' there is k₀ such that for k_e ≥ k₀ the T₀-periodic solution is globally unique and globally attracting.

**Argument.** Counterexample: n = 1, f(x) = x³ (C², autonomous), x_m = 0.5 sin t, μ = the test's 5 saturated atoms. For every k_e there is a second T₀-periodic solution near e ≈ √k_e:
- k_e = 40: e ∈ [5.56, 7.06], collocation residual 4e-12.
- k_e = 160: e ∈ [11.9, 13.4].
- k_e = 640: e ∈ [24.5, 26.0].

The solution starting at e₂(0) + 0.5 blows up in finite time: |e| > 1e6 at t = 0.025, 0.008 and 0.0026 respectively. The proof only uses bounds 'near the orbit', so only a local statement is proved. The code is unaffected: it starts from e = 0 and converges to the small solution.

Script: <scratchpad>/round3/averaging-impl/counterexample_a.py

**Fix.** Replace (a) by: '(a) *Periodic error.* Fix R > 0. There is k₀(R), independent of μ, such that for k_e ≥ k₀(R) the frozen-μ error equation ė = r_μ(x_m + e) − k_e e has exactly one T₀-periodic solution e_μ with sup_t‖e_μ‖ ≤ R. It attracts exponentially every solution with ‖e(0)‖ ≤ R, and e_μ = … (unchanged) uniformly in t and μ. If f is globally Lipschitz, uniqueness and attraction are global.'

Replace the proof bullet by: 'The tube {‖e‖ ≤ R} is forward invariant once k_eR > sup_{‖e‖≤R}‖r_μ(x_m+e)‖, and solutions in it contract at rate k_e − Lip_R(r_μ).'

In averaging.py:
- 'unique periodic solution e_μ' → 'unique periodic solution e_μ near 0';
- 'decays like exp(−k_e T₀)' → 'decays like exp(−(k_e − Lip r_μ)T₀)'.

#### minor: §0 table, row H.7, status column: 'expansion, Stokes identity, parity and non-reciprocity checked to about 1%'

**Claim.** The four H.7 checks agree to about 1 %.

**Argument.** At k_e = 80 the expansion test's relative errors are 0.247 (n=1), 0.0188 (n=2, +) and 0.0115 (n=2, −). The expansion is checked only through its O(1/k_e) convergence, not to 1 %.

Stokes holds to 6.8e-15 and parity to 2e-15; parity holds exactly by construction. Only non-reciprocity is 'about 1 %': 1.0 % (n=1) and 1.5 % (n=2).

The expansion does hold to 6e-8 after Richardson extrapolation (check_expansion.py), but no repository test does this.

Scripts: <scratchpad>/round3/averaging-impl/repro_tables.py and check_expansion.py.

**Fix.** Replace with: 'frozen-measure part proved; Stokes identity and orbit-reversal parity hold to machine precision (7e-15, 2e-15); the expansion is checked through its O(1/k_e) convergence (relative error 0.54 → 0.25 for n = 1 and 0.038 → 0.019 for n = 2 at k_e = 40 → 80); the non-reciprocity of the exact second variation is checked to 1.5 % at k_e = 80.'

#### minor: §8 item 2: '*Proved:* at finite gain, the averaged adaptive dynamics is a semi-gradient flow. Its antisymmetric interaction is generated by the area that the reference orbit encloses.'

**Claim.** The non-reciprocity comes from the enclosed area, and this is proved.

**Argument.** This contradicts H.7(e), §0 ('This already holds for n = 1, so it is not a two-dimensional effect') and the H.7 table. The antisymmetric part is A_even + 2K. A_even needs no enclosed area and is the only source for n = 1 (cyclic sum −0.0704 with K ≡ 0).

In the checked n = 2 instance, the state-shift cyclic sum (−0.0894) is twice the circulation part (+0.0430). With the orbit reversed the total is −0.1323.

'Proved' also overstates: (e) only shows non-vanishing cyclic sums in the checked instances, at a frozen measure.

Scripts: <scratchpad>/round3/averaging-impl/repro_tables.py and check_secondvar.py.

**Fix.** Replace with: '*Proved (at a frozen measure; non-reciprocity in the checked instances):* at finite gain, the averaged adaptive dynamics is a semi-gradient flow plus an orientation-odd lag term. Its non-reciprocity has two sources: the state shift (A_even), which needs no enclosed area and is the only source for n = 1, and the orbit-area circulation kernel 2K for n ≥ 2. In the checked n = 2 instance the state-shift part of the cyclic sum (−0.089) is about twice the circulation part (+0.043).'

#### minor: §0 bullet 'The prediction is accurate to 3% at k_e = 40 (`exp_orbit_orientation.py`)'; §0 table row H.7 'predicts the orientation dependence of the learned kernels in the real closed loop to 3% (cosine 0.994, k_e = 40)'; §8 item 2 'predicts the O(1/k_e) difference of the learned kernels to 3%'

**Claim.** The full first-order prediction of Θ₊ − Θ₋ is accurate to 3 %.

**Argument.** The 3 % is the size ratio only (0.969). Using results/orbit_orientation/raw_sigma2_slow16.npz and the experiment's own predicted_difference, the vector error ‖meas − full‖_{G⁻¹}/‖full‖_{G⁻¹} is 0.364, 0.195 and 0.109 at k_e = 10, 20 and 40. This is consistent with cosine 0.994 at ratio 0.969.

The theory is supported, since the vector error halves with each doubling of k_e. But '3 %' overstates the agreement about 3.6×. Remark H.7a and §7 correctly say 'in size'.

Script: <scratchpad>/round3/averaging-impl/orient_vec_err.py

**Fix.** - §0 bullet: 'At k_e = 40 the prediction matches the measured difference to 3 % in size and 11 % as a vector (cosine 0.994); the vector error halves with each doubling of k_e (36 %, 20 %, 11 % at k_e = 10, 20, 40).'
- §0 table: 'to 3 % in size and 11 % as a vector (cosine 0.994, k_e = 40)'.
- §8: 'to 3 % in size and 11 % as a vector at k_e = 40'.
- Add a vector-error column to the experiment's table.

#### minor: tests/test_averaging.py::test_exactly_averaged_dynamics_is_non_reciprocal_for_n1_and_n2, docstring 'times k_e² matches the predicted coefficient S_even + K to within 10%'; the formula is PARADIGM.md H.7 proof of (e), 'S_even(ϑ,ϑ') = −∫[(Φ̃_ϑ'·∇)(r_μᵀΦ̃_ϑ) + (r_μ·∇)(Φ̃_ϑ'ᵀΦ̃_ϑ)]dρ'

**Claim.** The test checks the full predicted second-variation coefficient S_even + K.

**Argument.** Only the cyclic sum of S − Sᵀ is compared. Dropping the symmetric term (r·∇)(Φ̃_jᵀΦ̃_i) from second_variation_formula still passes the real test (mutate.py svf: 'S_drop_term2 PASS'). The formula is in fact right: the full matrix is verified to 9e-9 by Richardson extrapolation (check_secondvar.py). But no repository test checks its symmetric half, and the 10 % tolerance is loose against the measured 1.0 % and 1.5 %.

Scripts: <scratchpad>/round3/averaging-impl/mutate.py and proposed_asserts.py.

**Fix.** Docstring: 'The cyclic sum of the antisymmetric part of k_e² times the exact second variation (k_e = 80) is within 3 % of that of S_even + K, and the full matrix k_e²S + k_eS₀ is within 5 % of S_even + K.'

Code, replacing the 'exact = …' and 'assert abs(exact - pred) < 0.1 …' lines:
Sx = A.exact_second_variation(f, fm, orb, ke, TE[:3])
Phi, _ = A.features(net, orb.xm, TE[:3]); S0 = np.einsum('tik,tjk->ij', Phi, Phi) / orb.N
assert np.max(np.abs(ke**2 * Sx + ke * S0 - (S_e + K))) < 0.05 * np.max(np.abs(S_e + K))
exact = ke ** 2 * A.cyclic_defect(Sx)
assert abs(exact - pred) < 0.03 * abs(pred), (n, exact, pred)

Measured: full matrix 0.0178 (n=1) and 0.0181 (n=2); dropping the symmetric term gives 1.86 and 0.27.

#### minor: tests/test_langevin.py::test_averaging_error_of_the_formal_reduction_decreases_like_one_over_ke: 'assert errs[-1] < 0.35 * errs[0], errs   # ~1/k_e'. tests/test_averaging.py::test_frozen_measure_expansion_converges_at_rate_one_over_ke: n = 1 threshold 'errs[1] < (0.3 if n == 1 else 0.03)'. §7 row 'Remark H.0 | … decays like 1/k_e, n = 1 and n = 2'.

**Claim.** The tests establish 1/k_e decay (n = 1) and the n = 1 expansion.

**Argument.** A factor 0.35 over an 8× range of k_e accepts any decay faster than k_e^(−0.505), so the test does not check '1/k_e'. The measured per-doubling ratios are 0.49, 0.45 and 0.52.

The n = 1 sub-case of the expansion test accepts U_even scaled by any λ ∈ [0.76, 1.24], and 72 % of random 10 % perturbations (thresholds.py). The test as a whole stays sharp only because the n = 2 sub-cases share the code path.

A Richardson check at k_e = 80/160 is sharp for n = 1: residual 0.0104 with the correct U, 0.045/0.056 for U_even ×1.05/0.95, run time 4 s (rich_n1.py).

Scripts: <scratchpad>/round3/averaging-impl/thresholds.py and rich_n1.py.

**Fix.** In test_langevin, replace the last assertion by:
r = errs[1:] / errs[:-1]
assert np.all(r > 0.4) and np.all(r < 0.62), errs   # halves when k_e doubles (measured 0.49, 0.45, 0.52)
assert errs[0] > 1e-3, errs

In the expansion test, for n = 1, add a k_e = 160 solve (Orbit(..., N=2000, sub=4)) and assert:
rich = 2*C160 - C80
assert np.max(np.abs(rich - (Ue + Uo))) < 0.03 * np.max(np.abs(Ue + Uo))

#### minor: §6c Remark H.7a 'Parity. b_osc is odd under reversal of the orbit (`test_oscillation_correction_is_odd_under_orbit_reversal`)'; §7 row 'H.7a, oscillation term | odd under orbit reversal (1e-6)'

**Claim.** b_osc (averaging.oscillation_correction) is checked.

**Argument.** Oddness under reversal holds exactly for any expression of the form ⟨B(H,u)⟩ with u the antiderivative. −b_osc or 2·b_osc pass this test equally, so the test never checks the value. The value is supported only indirectly, by the experiment's size ratio.

I verified the value independently (check_osc.py). A 2-second value test would guard it: a near-identity one-period flow at k_e = 80 gives relative error 0.0073, against 2.0 for −b_osc and 0.095 for 1.1·b_osc.

Scripts: <scratchpad>/round3/averaging-impl/check_osc.py and proposed_asserts.py.

**Fix.** Add to tests/test_averaging.py, above the __main__ block:
def test_oscillation_correction_matches_a_one_period_flow():
    net, Th, _ = _setup(2, seed=3); f, J, xr, xd = _plant_2d(+1); G = np.ones(5); ke = 80.0
    def H(T, t):
        x = xr(t); r = f(x[None, :])[0] - net.f_hat(x, T); return G[None, :] * net.grad_g(x, r, T)
    def u0m(T, Nu=2000):
        Hs = np.array([H(T, t) for t in np.arange(Nu) * PER / Nu]); hs = Hs - Hs.mean(0)
        I = np.concatenate([np.zeros((1,) + hs.shape[1:]), np.cumsum(0.5 * (hs[1:] + hs[:-1]), 0) * PER / Nu])
        return -I.mean(0), Hs.mean(0)
    Th0 = Th + u0m(Th)[0] / ke; X = Th0.copy(); dt = PER / 2000
    for k in range(2000):
        t = k * dt; k1 = H(X, t) / ke; k2 = H(X + dt / 2 * k1, t + dt / 2) / ke; k3 = H(X + dt / 2 * k2, t + dt / 2) / ke; k4 = H(X + dt * k3, t + dt) / ke
        X = X + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    P0 = Th0 - u0m(Th0)[0] / ke; P1 = X - u0m(X)[0] / ke; Ym = 0.5 * (P0 + P1)
    D = ke ** 2 * ((P1 - P0) / PER - u0m(Ym)[1] / ke)
    osc = A.oscillation_correction(net, f, A.Orbit(xr, xd, PER, N=2000, sub=1), Ym, G)
    assert np.max(np.abs(D - osc)) < 0.03 * np.max(np.abs(osc))

In §7: 'H.7a, oscillation term | value against a one-period near-identity flow (0.7 % at k_e = 80); odd under orbit reversal | Av::test_oscillation_correction_matches_a_one_period_flow, Av::test_oscillation_correction_is_odd_under_orbit_reversal'.

#### cosmetic: §6c Proposition H.7(e): 'Its antisymmetric part is (1/k_e²)(A_even + 2K) + O(k_e⁻³)'

**Claim.** The antisymmetric part of the second variation equals (1/k_e²)(A_even + 2K).

**Argument.** A_even is defined as the full difference S_even − S_evenᵀ, and cyclic_defect uses S − Sᵀ. The stated identity is right for S − Sᵀ (cyclic sums −0.0704 = cyclic(A_even + 2K)). With the standard antisymmetric part ½(S − Sᵀ), the right-hand side would be (1/k_e²)(½A_even + K).

**Fix.** Write: 'Its antisymmetric part, written S − Sᵀ (without the factor ½), is (1/k_e²)(A_even + 2K) + O(k_e⁻³)'.

#### cosmetic: §6c H.7 'Checked' table, row non-reciprocity: 'the leading-order interaction is reciprocal to 1e-14'; tests/test_averaging.py 'assert abs(A.cyclic_defect(S0)) < 1e-14'

**Claim.** A check shows that the leading-order interaction is reciprocal.

**Argument.** S0 = einsum('tik,tjk->ij', Phi, Phi)/N is symmetric by construction. S0 − S0ᵀ is exactly 0.0 (repro_tables.py), so the assertion cannot fail and tests nothing about the exact second variation.

**Fix.** Text: 'the leading-order interaction −(1/k_e)⟨Φ̃_ϑᵀΦ̃_ϑ'⟩_ρ is symmetric by construction; that k_e·S_exact + S₀ = O(1/k_e) is checked by the full-matrix assertion'. Replace the 1e-14 assertion by the full-matrix assertion of the previous test-docstring issue.

#### cosmetic: §6c H.7 'Checked' table: '(… 0.038 → 0.019 for n = 2, at k_e = 40 → 80; both orientations)', 'agreement to 1e-6 relative', 'to 1e-10 relative'; the same in results/RESULTS_PARADIGM.md

**Claim.** Both n = 2 orientations give 0.038 → 0.019. Stokes agrees to 1e-6 and parity holds to 1e-10.

**Argument.** The clockwise orbit gives 0.0225 → 0.0115. The actual agreements are 6.8e-15 (Stokes) and 2.2e-15 (parity). The quoted values are the assertion tolerances, not the measurements (repro_tables.py).

**Fix.** '(0.54 → 0.25 for n = 1; 0.038 → 0.019 (counter-clockwise) and 0.022 → 0.011 (clockwise) for n = 2, at k_e = 40 → 80)'; 'agreement to 7e-15 relative (asserted at 1e-6)'; 'to 2e-15 relative (asserted at 1e-10)'.

#### cosmetic: §6c Remark H.0 table: '| n = 1 | 4.2% | 2.1% | 0.91% | 0.48% |' and '| n = 2, orbit enclosing area | 43% | … | `test_averaged_drift_error_decays_like_one_over_ke_for_an_orbit_enclosing_area` |'; the same in results/RESULTS_PARADIGM.md

**Claim.** The cited tests produce these entries, and the n = 2 error is associated with the enclosed area.

**Argument.** - The cited Lg test gives 2.047 % at k_e = 20, which rounds to 2.0 %; 2.1 % is the RK4 value (2.054 %).
- The n = 2 test runs only k_e = 20, 40 and 80. The 43 % entry (42.7 %, reproduced with the same code) is not produced by the cited test.
- The label 'orbit enclosing area' suggests the area causes the 10× larger error. H.7 attributes it to the orientation-even U_even: k_e×error 4.05, against 1.53 from U_odd alone, and the reversed orbit gives 4.0 (split_H0.py, repro_H0.py).

**Fix.** Change the n = 1 row to '| n = 1 | 4.2% | 2.0% | 0.91% | 0.48% |' and rename the second row 'n = 2, ellipse'. Add below the table: 'The k_e = 10 entry of the n = 2 row uses the same code but is not asserted (the test runs k_e = 20, 40, 80). The n = 2 error is dominated by the orientation-even state-shift term: H.7 predicts k_e × error 4.05, of which U_odd alone would give 1.5; the reversed orbit gives 4.0.'

#### cosmetic: tests/test_averaging.py: 'def test_oscillation_correction_is_odd_under_orbit_reversal' is defined after the 'if __name__ == "__main__":' block. PARADIGM.md §6c: '*Checked* (`tests/test_averaging.py`, five tests, using `src/core/averaging.py`)'.

**Claim.** Running the file executes all its tests, and the file has five tests.

**Argument.** Running `python tests/test_averaging.py` executes only five tests; the new sixth one is defined after the loop runs. I verified this: script mode prints 5 'passed' lines, while pytest collects 6.

**Fix.** Move test_oscillation_correction_is_odd_under_orbit_reversal (and any new test) above the __main__ block. Change the text to '*Checked* (`tests/test_averaging.py`, six tests: four for H.7, one for Remark H.0, one for Remark H.7a, using `src/core/averaging.py`)'.

## Refuter: H8-H0-H6

### Holds

- Version note: PARADIGM.md was rewritten during this review (mtime 01:00:39; it grew from 1609 to 1678 lines and gained Remark H.7a and the b_osc sentence in the Route). Everything below refers to that current version. A copy is at <scratchpad>\round3\H\PARADIGM_snapshot_0127.md. All my scripts are in that same folder.
- H(i) holds. The first variation is exact. Writing F_τ = (1/k_e)L̄ + τH(·|γ) + const with γ = N(ϑ⁰,(τ/σ)G) is correct, and the constant is τ[(p/2)log 2π + ½log det((τ/σ)G)]. 'Finite exactly when H(μ|γ) < ∞' is true: 0 ≤ L̄ ≤ ½(F_ρ+K_Φ)² for every μ ∈ P(ℝ^p), and H(μ|γ) < ∞ forces μ ∈ P₂ (Donsker–Varadhan with ψ = ε‖ϑ‖²).
- H(ii) holds. Whitening, Assumptions 1–3 of [S3], Bakry–Émery (σ/τ) and the Holley–Stroock oscillation 2K_Φ(F_ρ+K_Φ)/(k_eτ) all check. Thm 3.2, Cor 3.3, Lemma 3.4 and Prop 5.1 are transcribed correctly; I read them in mfld.txt.
- The P₂-start argument in H(ii) holds, constant included. Freezing the flow of marginals turns the McKean–Vlasov particle into a linear SDE with a bounded drift difference, so Girsanov (Novikov holds trivially) plus data processing give KL ≤ (1/(4τ))E∫‖b‖² ≤ B²t₀/(4τ). The constant is sharp: for constant b and σ = 0, KL(N(bt₀, 2τt₀) | N(0, 2τt₀)) = b²t₀/(4τ) exactly. The OU density bound (2π(τ/σ)(1−e^{−2σt₀}))^{−p/2} and the bound B are correct.
- H(iii) holds. The continuity constant (F_ρ+K_Φ)L_Φ is right, with L_Φ the dual-norm Lipschitz constant in the G⁻¹ metric. The bias bound L̄(μ*_τ) ≤ inf_ν{L̄(ν) + k_eτH(ν|γ)} is right.
- H(iv) holds. In [S3] Thm 4.1, ρ_τ is used only along τ_t ≤ τ₀ (Step 1, t ≥ t₀). Lemma 4.2 needs only some LSI at τ = 1, which holds. 𝒢 ≥ 0, and α* = 2K_Φ(F_ρ+K_Φ)/k_e matches the exponent.
- Proposition H.8 holds, given the chain-rule assumption. Young's inequality with ε = 1 is the optimal split. Then ∫‖ξ‖² = τ²I(μ|μ̂) ≥ 2ρ_ττ²H(μ|μ̂) ≥ 2ρ_ττ(F−F*), by the LSI and the upper entropy sandwich. Grönwall gives the first two lines. With Talagrand for μ*_τ (itself a proximal Gibbs measure) and the lower sandwich, limsup W² ≤ (2/(τρ_τ))·β̄²/(2τρ_τ), so the factor β̄/(τρ_τ) is right. It is sharp for OU with a constant drift (b/σ). The chain-rule caveat is stated and is listed as missing item 3.
- The Route's scaling claims hold. In the H.6 scaling the whitened proximal Gibbs measures do not depend on k_e, so τρ_τ = Θ(1/k_e). β̄ = O(k_e⁻²) holds including the oscillation term. On the n = 1 test instance, k_e²|b_frozen| ≈ 0.29–0.32 and k_e²|b_osc| ≈ 0.072 are constant over k_e = 10, 20, 40, 80, with |b_osc|/|b| ≈ 0.25 (second_order.py). With dither, the extra Itô and diffusion terms of the averaging transformation are zero-mean at O(k_e⁻²), so drift-only perturbations suffice at this order. The radius O(1/k_e)e^{2K_Φ(F_ρ+K_Φ)/τ₀}/σ₀ follows.
- Every quoted number reproduces (consts.py, gibbs_dual.py, exact_avg2.py).
- Constants: F_ρ = 1.24688 and α* = 4.1975. log(2τρ_τ) ≥ −8400.53 and −1052.60, so 'e^{−8401}' and 'e^{−1053}' are valid lower bounds. L̄(0) = 0.3113.
- Grid minimiser: L̄(μ*_τ) = 0.111183, by an independent dual solver on two grids.
- Exactly averaged stationary state at k_e = 10: L̄ = 0.11415 (residual 6e-11), against 0.11118.
- Remark H.0 table: n = 1 errors 4.198%, 2.047%, 0.914%, 0.475%; n = 2 errors 42.7%, 21.2%, 10.4%, 5.14%, with k_e·error = 4.27 → 4.08.
- The tables and fits reproduce.
- The τ table matches results/langevin_tradeoff/table.md.
- That run used G = diag(4,1,1): the deterministic τ = 0 trajectories reproduce exactly with this G and not with diag(2,.5,.5) (tradeoff_G.py).
- Fit coefficients −1693 + 427·log(1/τ), residuals 0.045 and 0.253, and thresholds 0.022 → 0.13 reproduce.
- Annealing α values are 1.56e-3, 1.25e-2 and 2.03e-2.
- Loss ratios 1.03 and 5.2, RMS 0.0465/0.0452, W 0.403/1.691 (4.2×), and 0.166 → 0.048 match the tables.
- Further checks are consistent with the text.
- Finite m: the m = 32 averaged-system Gibbs measure has E[L̄(μ^32)] = 0.113 by MALA, with instantaneous sd 0.011 (finite_m_gibbs.py). So 'about 0.11 matches 0.111' is fair.
- At σ = 0.002, τ = 5e-4, L̄(μ*_τ) ≈ 0.0085 (gibbs_cold.py), against spread 0.009 and collapsed 0.048. This supports the metastability reading.
- In the H.6 scaling, ‖f̂_avg − f̂*‖_{L²(ρ)} = 0.0068, 0.0014, 0.0010, 0.0006 at k_e = 10, 20, 40, 80, which supports r(k_e) → 0. k_e = 10 is pre-asymptotic: the L̄ gap goes +0.0030, +0.0003, −0.0002, −0.0002.
- Non-reciprocity at the experiment's own network (a_max=4, s_min=0.15): the cyclic defect of S_even is O(1) relative to max|S_even| over 20 random triples, and K ≡ 0 to 1e-15 (nonrecip_exp.py).
- All round-2 fixes H.1–H.13 are applied: Remark H.0, the H.6 restatement, P₂/Girsanov, autonomous f, the new (iii), the new (iv), the conditional H.5, the Remark H.2 corrections, titles, preamble, tests, and bookkeeping including expit and the τ ≤ 0 guard. tests/test_langevin.py and tests/test_averaging.py: 13/13 pass.

### Issues

#### minor: §6c Remark H.2, bullet 'The fit.': 'It suggests that the symmetric configuration is a *saddle* of the averaged landscape, not a barrier-protected minimum. That would explain why practice is far faster than the worst-case bound.'

**Claim.** The log(1/τ) fit (4.5% residual) beating Arrhenius (25%) is evidence for saddle escape from the collapsed configuration.

**Argument.** 1. The contrast exists only after averaging the two seeds, whose escape times differ by up to 2.3× (970 vs 430 s at τ = 4e-3). On the eight per-seed values the log(1/τ) law has a 27% relative residual and the Arrhenius law 29% (fit.py).
2. I linearised the splitting dynamics at the noiseless collapsed end point (saddle_exact.py, saddle_check.py; σ = 0.002, k_e = 10, G = diag(4,1,1)).
   - For (𝔄), the Dirac fixed point is ϑ̄ = (1.199, 0.939, −1.431), with L̄ ≈ 0.166. G·Hess V[δ_ϑ̄](ϑ̄) has eigenvalues (−5.4e-5, 6.0e-3, 2.6e-2) s⁻¹.
   - For the exactly averaged loop at k_e = 10 they are (−2.2e-4, 6.0e-3, 2.4e-2) s⁻¹.
   - So the configuration is a saddle, but a very weak one. A saddle-escape law t ≈ (1/(2λ))log(ℓ²λ/τ) needs τ ≪ ℓ²λ ≈ 1e-5 (ℓ ≈ 0.5), and its slope would be 1/(2λ) ≈ 2·10³–10⁴ s, not 427 s.
   - At the tested τ ≥ 1e-3, motion along the splitting direction is diffusion-dominated.
3. The threshold depends on τ, and the escape time includes a transient of about 500 s.
The conclusion 'saddle' happens to be right, but the fit is not evidence for it.

**Fix.** Replace the bullet with: "**The fit.** The seed-averaged escape time fits t_esc ≈ −1693 + 427·log(1/τ) (4.5% relative residual) better than an Arrhenius law (25%). On the eight per-seed values the two laws fit equally poorly (27% and 29%). The noiseless collapsed end point is a saddle of the averaged landscape, but a weak one: its unstable rate is ≈ 5·10⁻⁵ s⁻¹ for (𝔄) and ≈ 2·10⁻⁴ s⁻¹ for the exactly averaged loop at k_e = 10. A saddle-escape law would therefore need τ ≲ 10⁻⁵ and would have a log(1/τ) slope of 10³–10⁴ s. At the tested temperatures escape is governed by diffusion, by the approach transient and by the τ-dependent threshold. The data do not identify the escape mechanism."

#### minor: §6c Remark H.0, bullet 'The first corrections': 'are O(1/k_e) relative to the drift and do not depend on the gains. Proposition H.7 gives them exactly at frozen μ.'

**Claim.** The first-order corrections to (𝔄) are independent of the gains.

**Argument.** Remark H.7a, added in v3, says the oscillation term b_osc has size, relative to the frozen-measure terms, scaling like |G|T₀, and that it 'is not negligible as k_e → ∞ at fixed G'. That is exactly the formal limit Remark H.0 names in its next bullet. Numerically, on the n = 1 test instance with G = diag(4,1,1), k_e²|b_osc| ≈ 0.072 is constant over k_e = 10–80, about 25% of the frozen term (second_order.py), and b_osc ∝ G². So the full first correction along the loop does depend on the gains. The sentence is true only for the frozen-(μ, ϑ) correction that the table measures.

**Fix.** Replace with: "**The first corrections.** At frozen μ (and frozen ϑ) they are O(1/k_e) relative to the drift and do not depend on the gains; Proposition H.7 gives them exactly. Along the loop, the within-period oscillation of the parameters adds the second-order averaging term b_osc of Remark H.7a. Its size relative to the frozen terms is O(|G|T₀), so it does depend on the gains, and in the formal limit below both are O(k_e⁻²). Measured relative errors of the frozen-μ averaged drift against (𝔄):"

#### minor: §6c Remark H.0, 'Period averaging: T₀·sup‖G∇V[μ]‖ ≪ ℓ and τλ_GT₀ ≪ ℓ²'

**Claim.** This is a sufficient condition for the period-averaging step.

**Argument.** G∇V[μ] = (1/k_e)G∇V_L[μ] + σ(ϑ − ϑ⁰). For σ > 0 its supremum over ϑ ∈ ℝ^p (the support of μ_t when τ > 0) is +∞, so the condition can never hold as written. The confinement drift does not depend on time and needs no averaging. Only the oscillating, loss-driven part needs to be small over a period, and that part is bounded.

**Fix.** Replace with: "*Period averaging:* T₀·sup_ϑ‖(1/k_e)G∇V_L[μ](ϑ)‖ ≪ ℓ (this sup is at most (T₀/k_e)λ_G(F_ρ+K_Φ)·sup‖∇_ϑΦ̃‖), σT₀ ≪ 1 and τλ_GT₀ ≪ ℓ², where ℓ is the ϑ-scale of Φ̃. The confinement drift does not depend on time and needs no averaging."

#### minor: §6c Remark H.0 'Tikhonov: k_e ≫ max(ω_ref, Lip(f) + K_Φ/s_min, (F_ρ + K_Φ)/s_min)', together with the Evidence 'Limits of the evidence' list

**Claim.** The Evidence supports the averaged theory and lists its limits.

**Argument.** Every Langevin run uses k_e = 10 with a_max = 4, s_min = 0.15 and f = sin 2x + 0.3x. For this network:
- (F_ρ+K_Φ)/s_min = 5.247/0.15 ≈ 35;
- Lip(f) + K_Φ/s_min = 2.3 + 26.7 ≈ 29.
So the runs sit about 3.5× below the stated Tikhonov scale, not ≫ above it. The limits list mentions n = 1 and fixed k_e but not this. The measured frozen-μ error at k_e = 10 is 4.2% for this very network (h0_table.py). So the condition is conservative, and the agreement is empirical rather than covered by the formal regime.

**Fix.** Add to 'Limits of the evidence': "All Langevin runs use k_e = 10, below the Tikhonov scale of Remark H.0 for this network ((F_ρ+K_Φ)/s_min ≈ 35). The measured frozen-μ drift error at k_e = 10 is 4.2% (n = 1), so the condition is conservative, but the agreement with (𝔄) is empirical."

#### minor: §6c Remark H.3: 'Such dynamics satisfy a log-Sobolev inequality with constant of order 1/diam², even for σ = 0. The argument is Bakry–Émery with a Neumann boundary, then Holley–Stroock'

**Claim.** Bakry–Émery with a Neumann boundary gives an LSI of order 1/diam² at σ = 0.

**Argument.** At σ = 0 the reference measure is uniform on Θ_K, and Hess V = 0. Bakry–Émery needs curvature ≥ κ > 0; with a convex Neumann boundary it yields only CD(0,∞), which gives no LSI by itself. The 1/diam² LSI for convex bodies is a different theorem: either a diameter-based LSI for log-concave measures on convex bodies (Frieze–Kannan 1999; Bobkov 1999), or CD(0,∞) plus a diameter bound (F.-Y. Wang 1997). The claimed constant is right; the cited argument does not give it.

**Fix.** Replace the second bullet with: "The argument is a diameter-based log-Sobolev inequality for the uniform measure on a convex body (CD(0,∞) with a convex Neumann boundary plus a diameter bound, e.g. F.-Y. Wang 1997; or Frieze–Kannan 1999 for log-concave measures on convex bodies), then Holley–Stroock. This gives LSI(c·diam⁻²·e^{−2K_Φ(F_ρ+K_Φ)/(k_eτ)}). Bakry–Émery alone gives nothing at σ = 0."

#### cosmetic: §6c Evidence, 'Equalisation' bullet: 'their final measures are 4.2× closer in W than without dither (0.40 against 1.69)'

**Claim.** This W comparison shows that the two starts equalise.

**Argument.** The W here is not the theorem's metric, and the two runs are coupled.
- exp_langevin_global.w2_phys computes W₂ over the physical coordinates (a, c, s) between the two 32-particle final configurations. It is not W_{G⁻¹} over ϑ.
- run() seeds the noise with default_rng(100 + seed) whatever the start, so collapsed and spread runs of one seed share their Brownian increments. 0.40 therefore measures a synchronous coupling.
- Two independent 32-point samples of μ*_τ are W₂ = 0.67 apart on average (median 0.60, 10–90% range 0.47–0.99; w2_floor.py). So 0.40 is at or below the sampling floor: consistent with equal laws, but no finer test of them.

**Fix.** Replace with: "their final 32-particle configurations, driven by the same noise realisation, are 4.2× closer in W₂ over the physical coordinates (a, c, s) than without dither (0.40 against 1.69). Two independent 32-point samples of μ*_τ are typically about 0.6 apart, so 0.40 is at the sampling floor;"

#### cosmetic: §6c 'Route to Conjecture H.6', bullets 1 and 3: 'This is O(1/k_e), with a constant independent of k_e.' and 'Proposition H.8 then confines every solution of "(𝔄) + b"'

**Claim.** These state the scaling of τρ_τ and what Proposition H.8 then gives.

**Argument.** 1. The bullet states a lower bound, but 'O(1/k_e)' is an upper-bound notation. What the ISS gain needs is 1/(τρ_τ) = O(k_e). Both directions do hold, because in this scaling the whitened proximal Gibbs measures do not depend on k_e.
2. The symbol b in bullet 3 is defined in bullet 2 as the frozen-μ term only, while bullet 2 adds b_osc as 'another' term.
3. H.8 needs F_τ(μ_0) < ∞.
4. The constant at the experiment's τ₀ = k_eτ = 0.04 is e^{1049}, which the reader should see.

**Fix.** Bullet 1: "τρ_τ ≥ (σ₀/k_e)·e^{−2K_Φ(F_ρ+K_Φ)/τ₀}, so 1/(τρ_τ) = O(k_e). In this scaling the whitened proximal Gibbs measures do not depend on k_e, so τρ_τ = Θ(1/k_e). At the experiment's τ₀ = 0.04 the constant is e^{1049}." In bullet 2, define b := (1/k_e²)[G∇(U_even + U_odd) + b_osc] + O(k_e⁻³). Bullet 3: "Proposition H.8 then confines every solution of (𝔄) + b with F_τ(μ_0) < ∞ (or from P₂, after the Girsanov regularisation of H(ii), which applies because b is bounded) to …"

#### cosmetic: §6c Conjecture H.6, last bullet: 'By Proposition H.7 it is not a gradient flow.'

**Claim.** The exactly averaged system at fixed k_e is not a gradient flow.

**Argument.** H.7(e) is now hedged: 'So in these instances … It is not the W_{G⁻¹}-gradient flow of any free energy. A different, μ-dependent metric is not excluded.' Its checked instances use a_max = 3, c_max = 2, s_min = 0.3 and random 5-atom measures. H.6 asserts the unhedged version. The claim does hold at the experiment's own network: the cyclic defect of S_even is O(1) relative to max|S_even| for 20 random triples, with K ≡ 0 to 1e-15 (nonrecip_exp.py).

**Fix.** Replace with: "By Proposition H.7(e) it is not the W_{G⁻¹}-gradient flow of any free energy: the cyclic sums of its antisymmetric second variation are nonzero in the checked instances, including the network of the experiments."

#### cosmetic: §6c notation: Design 'ϑ = (α, γ, κ)'; H(i) 'γ := N(ϑ⁰, (τ/σ)G)'; H(iv) 'α > α*', 'F_{τ₀}(μ_0)', 'C₀ = σ/τ₀'; H.6 and the Route 'τ := τ₀/k_e', 'e^{−2K_Φ(F_ρ+K_Φ)/τ₀}'

**Claim.** The symbols are used consistently.

**Argument.** Three symbols each name two things:
- γ is both the unconstrained centre coordinate and the Gaussian reference measure.
- α is both the weight coordinate and the annealing constant.
- τ₀ is the initial annealing temperature in (iv), and k_eτ in H.6 and the Route.
The τ₀ clash matters most, because both meanings sit in exponents a few lines apart: C₀ = σ/τ₀ against e^{−…/τ₀}.

**Fix.** Rename the Gaussian reference measure to γ_τ (or 𝒩_τ), the annealing constant to α_a, and the initial annealing temperature to τ_init (so C₀ = σ/τ_init). Keep τ₀ := k_eτ for H.6 and the Route.

#### cosmetic: §6c Remark H.0 table: rows 'n = 1 | 4.2% | 2.1% | …' and 'n = 2, orbit enclosing area | 43% | …', with their test names

**Claim.** The cited tests produce these numbers.

**Argument.** 1. test_averaged_drift_error_decays_like_one_over_ke_for_an_orbit_enclosing_area evaluates only k_e = 20, 40, 80, so the 43% at k_e = 10 is not computed by the cited test. The value itself is right (42.7%, h0_table.py).
2. The n = 1 value at k_e = 20 is 2.047%, which rounds to 2.0%, not 2.1%.
3. The n = 1 test asserts only monotone decrease and errs[-1] < 0.35·errs[0]. That allows a decay as slow as k_e^{−0.5}, weaker than its name 'decreases_like_one_over_ke'.

**Fix.** Change '2.1%' to '2.0%'. Either add k_e = 10 to the n = 2 test or footnote the 43% as 'computed separately'. Tighten the n = 1 assertion, e.g. assert 0.4 < errs[i+1]/errs[i] < 0.62 (the observed ratios are 0.49, 0.45, 0.52).

#### cosmetic: §6c Remark H.1: 'a grid computation gives L̄(μ*_τ) = 0.111 (review round 2). The unregularised offline optimum over 32 kernels is about 7·10⁻⁶'; Conjecture H.6: 'L̄ = 0.1141 against L̄(μ*_τ) = 0.1112 at k_e = 10 (review round 2)'

**Claim.** These are the quoted reference values.

**Argument.** The values are correct: I get 0.11118 by an independent dual solver and 0.11415 for the exactly averaged stationary state (gibbs_dual.py, exact_avg2.py). But no script or table in the project computes them; they come from a review scratchpad. The offline optimum in both results tables is 7.63e-06, the best of 5 L-BFGS restarts and hence an upper bound. That rounds to 8·10⁻⁶.

**Fix.** Write "about 7.6·10⁻⁶ (best of 5 L-BFGS restarts)". Add a project script (e.g. experiments/exp_gibbs_reference.py, writing results/langevin_global/gibbs_reference.md) that computes L̄(μ*_τ) and the exactly averaged stationary state, and cite it instead of 'review round 2'.

#### cosmetic: §6c Proposition H.8: 'This is true for smooth bounded perturbations, but it is not verified here in general.'; Theorem H(ii): 'and the estimates hold from t₀ on'; Conjecture H.6: 'satisfy the tube condition of Theorem D(b)'

**Claim.** Three side statements are precise as written.

**Argument.** 1. The first clause of the H.8 sentence is asserted without proof or reference. [S3] Lemma 2.2 needs ∫‖∇log μ_t‖²dμ_t dt < ∞, and the AGS bypass used in [S3] Thm 3.2 does not apply to a non-gradient perturbed flow.
2. 'Hold from t₀ on' does not say what replaces F_τ(μ_0).
3. Theorem D(b) is stated with a_max as the bound on ‖∫Φdμ‖. For the saturated features that bound is K_Φ = a_max√n.

**Fix.** 1. H.8: "It is expected for smooth bounded perturbations (parabolic regularity and the entropy-dissipation estimate), but it is not verified here."
2. H(ii): "…and the estimates hold with (μ_0, t) replaced by (μ_{t₀}, t − t₀)."
3. H.6: "satisfy the tube condition of Theorem D(b) with a_max replaced by K_Φ".

## Refuter: G

### Holds

- (b′) flow formula and its derivation hold. From the exact identity L̄ = L* + ½(a−a*)ᵀM(a−a*), ∂²_{ξa}L̄ = −Da*ᵀM at a = a*(ξ), Da* = −M⁻¹∂²_{aξ}L̄ and h − a* = −(k_e/(mg_a))M⁻¹Dh ξ̇, I get (Γ⁻¹ + (δ/g_a)Da*ᵀDh)ξ̇ = −(mδ/k_e)∇L*. The m/k_e factors cancel correctly. Replacing Dh by Da* costs only O(S/(R(1+S))) relative, so there is no missing term at leading order in 1/R.
- (b′) checked numerically against the joint averaged flow started on the slow manifold. m = 1 slope configuration: agreement within 1.5% (S = 2–2.9), 0.3% (S = 0.2–0.3) and 0.2–0.6% (S = 8.5–17). m = 2 with anisotropic Γ, S = 1–3, R = 73–377: size within 0.1–0.3%, cosine ≥ 0.99999. This confirms the kernel coupling through Da* and where Γ sits. Scripts: scratchpad\round3\propG\n1_slow_manifold.py, n6_m2_bprime.py, n13_large_S.py.
- The induced-metric identification holds. With the unnormalised G⁻¹ metric, the tensor induced on the graph of a* is (δΓ)⁻¹ + g_a⁻¹Da*ᵀDa*, and −(m/k_e)·(that tensor)⁻¹∇L* equals (b′).
- (a) holds: ϑ̇_i = −(m/k_e)G∇_{ϑ_i}L̄ = −(1/k_e)G∇(δL̄/δμ)(ϑ_i), and this matches the (𝔄) drift at σ = τ = 0. The relation to H(i) is worded correctly.
- Remarks on the hypotheses hold: m₀ ≤ 1/m², so the slowest weight rate is at most g_a/(k_e m); the (H2c) sufficient condition (evaluated at a = a*) is valid; R and S contain neither k_e nor ω_ref. (H2b) does not imply (H2c): in the slope configuration R·S ≥ 166 for every gain shape Γ tested.
- (b) holds as the S → 0 limit of (b′), and the plain-VarPro error scales like S. In the slope rows, (b) is off from (a) by 22% at S = 0.2–0.3 and by 2% at S = 0.02.
- 'No mean-field form' holds. (b) is the gradient flow of L*/k_e in (1/m)Σ‖ξ̇_i‖²_{(δΓ)⁻¹}. [S2] (arXiv:2504.18208, eq. 10) defines ℒ^λ(μ) = min_u (1/2λ)‖Φ_μu − Y‖²_{L²(ρ)} + ‖u‖²_{L²(μ)}, so λℒ^λ → L*. ℒ⁰ is min-norm interpolation and its flow is a weighted ultra-fast diffusion, as the document says. The support dependence and the W-discontinuity argument are correct.
- Every number in the Numerics table matches table_violating.md, table_satisfying.md and table_slope.md. My own scalar re-implementation (scratchpad\round3\propG\cl1.py) reproduces the δ = 0.47 slope row: closed-loop Δξ = (−0.00298, −0.02904), ratios 1.098 / 1.166 / 0.439, max|a| = 2.27.
- The t_first/t_last alignment is consistent. Closed-loop period means are centred at 0.475T and (n_per − 0.525)T, while the flows are evaluated at 0.5T and (n_per − 0.5)T. The offset is the same at both ends, so the durations match exactly and the effect is negligible.
- 'At small δ the 4% is the O(1/k_e) error-channel correction; it does not depend on δ' holds. At δ = 0.0047, closed/(a) is 1.0439 at k_e = 20 and 1.0109 at k_e = 80; at δ = 0.047 and 0.0047 it is 1.043 and 1.044 (scratchpad\round3\propG\n3_decompose.py).
- 'In all runs, L̄ at period-mean parameters is non-increasing' holds, including the three slope runs (0 increasing steps). Note that exp_slow_manifold_slope.py does not record L̄ itself.
- Step 1's exact identity ẇ = −k_e w − ṙ/k_e and the bound ‖w‖ ≤ e^{−k_e t}‖w(0)‖ + sup‖ṙ‖/k_e² are correct.
- The round-2 fixes for G.1–G.10 and G.12 were applied: (H0), (H2a′), Γ-weighted (H2b), (H2c), the induced-metric flow, the no-mean-field paragraph, the exact step-1 identity, the composition caveat, the o(1) limit with its initial condition, the S and joint-rate columns, centred differences, and the relabelled rows. Three fixes were applied but are themselves inaccurate (issues 1, 5, 6), and the G.11 attribution text is wrong (issue 3).

### Issues

#### major: §6b, Proposition G(b′), first display: 'After a boundary layer, the weights follow the slow manifold: a − a*(ξ) = (δ/g_a)M⁻¹Da*Γ∇L*(ξ) + h.o.t.'

**Claim.** Under (H0), (H1), (H2a′) and (H2b), without (H2c), the slow-manifold offset is (δ/g_a)M⁻¹Da*Γ∇L*.

**Argument.** Step 3's own identity h − a* = −(k_e/(mg_a))M⁻¹Da*ξ̇ must be evaluated with the (b′) velocity, not the VarPro velocity. That gives (δ/g_a)M⁻¹Da*(Γ⁻¹ + (δ/g_a)Da*ᵀDa*)⁻¹∇L*. The displayed formula is the S → 0 truncation of this, and its 'h.o.t.' are O(S) relative, which is O(1) in exactly the regime where (b′) is needed. For m = 1 the displayed offset is exactly (1 + S_loc) times the correct one.

Evidence:
- Slope configuration, δ = 0.47, along the joint averaged flow after the boundary layer (t = 80–230 s): actual a − a* = −0.0372 to −0.0399; displayed formula −0.124 to −0.154 (actual/displayed = 0.26–0.30 = 1/(1+S)); corrected formula within 0.4–0.8%. At δ = 0.047 the displayed formula is 20–25% off. Script: scratchpad\round3\propG\n1_slow_manifold.py.
- m = 2, S = 3: actual offset (−0.0340, −0.0002); displayed (−0.311, +0.088), so even a sign is wrong; corrected (−0.0337, −0.0005). Script: n6_m2_bprime.py.
- The display contradicts the displayed flow. Substituting it into ∇_ξL̄ = ∇L* − Da*ᵀM(a − a*) gives ξ̇ = −(mδ/k_e)(Γ − (δ/g_a)ΓDa*ᵀDa*Γ)∇L*. At ξ0, δ = 0.47, this has dL*/dt = +3.7e-6 (ascent), against −1.5e-6 for (b′). Over 232 s it raises L* by 3.9e-4 and has cosine −0.44 with the joint flow (n13_large_S.py).
- Origin: the round-2 judge proposed this formula for step 3 of (b), which assumes S ≪ 1. Version 3 moved it into (b′), where S need not be small.

**Fix.** Replace the display with: 'a − a*(ξ) = −(k_e/(mg_a))M⁻¹Da*ξ̇ + h.o.t. = (δ/g_a)M⁻¹Da*(Γ⁻¹ + (δ/g_a)Da*ᵀDa*)⁻¹∇L*(ξ) + h.o.t., which reduces to (δ/g_a)M⁻¹Da*Γ∇L*(ξ) when S ≪ 1.' In step 3, add after 'h − a* = −(k_e/(mg_a))M⁻¹Dh ξ̇': '(evaluated with the (b′) velocity, not the VarPro velocity)'.

#### minor: §6b, Proposition G(b′) ('After a boundary layer, the weights follow the slow manifold … + h.o.t.') and Numerics bullet 'The weight block. In the `satisfying` run the weights sit on the slow manifold.'

**Claim.** In the closed loop the weights sit on the slow manifold, at the displayed offset plus higher-order terms.

**Argument.** (b′) is derived for the averaged flow (a). In the closed loop, the weight equilibrium is shifted by O(ω_eff/k_e) (Proposition H.7). This shift is not of higher order than the O(δ/g_a) offset.
- Stored satisfying run: |a_pm − a*| = 3.1e-3 to 5.9e-3 for t = 515–2864 s, against a slow-manifold offset |(δ/g_a)M⁻¹Da*H⁻¹∇L*| of 4e-5 to 2.9e-4. The shift is 10–100 times larger (scratchpad\round3\propG\n5_satisfying.py).
- The shift scales like 1/k_e: in the slope configuration with frozen geometry it is +0.00803 at k_e = 20 and +0.00197 at k_e = 80 (n4_Lbar_increase.py).

**Fix.** Begin (b′) with 'For the averaged system (a), that is, to leading order in ω_eff/k_e, …'. After the slow-manifold display, add: 'In the closed loop the weight equilibrium is further shifted by O(ω_eff/k_e) (Proposition H.7). This shift can exceed the displayed offset: 3–6e-3 against 4e-5–3e-4 in the `satisfying` run.' In the weight-block bullet, replace 'the weights sit on the slow manifold' with 'the weights sit at the closed loop's shifted equilibrium a* + O(1/k_e), within 0.1% of a*'.

#### minor: §6b Numerics, bullet 'The size ratios': 'In the `satisfying` run, the centred ratio of 0.97 is that same correction (frozen-parameter ratio 0.956 at k_e = 20, 0.996 at k_e = 80; review round 2) plus about 0.6% from S.'

**Claim.** 0.97 equals the error-channel ratio 0.956 plus a 0.6% contribution from S.

**Argument.** This attribution is wrong in both size and sign.
- The S term is negligible. Along the run, the ratio of (b′) to VarPro is 0.9999–1.0000, so S changes the predicted speed by at most 0.01%. With Γ = I the induced metric can only slow the flow, so it cannot raise 0.956 toward 0.97.
- The error-channel ratio does not explain the gap. The frozen-parameter ratio stays at 0.951–0.956 along the whole run, whether evaluated at a* or at the run's period-mean weights. Yet the centred closed-loop ratio is 0.9693–0.9703 (n5_satisfying.py).
- The +1.4% is an averaging effect from the weights oscillating within each period. With the weights running under (L) at frozen geometry (g_a = 13), the averaged geometry velocity ratio is 0.969–0.970, which matches the run. With g_a/4 it falls to 0.956–0.960, i.e. the excess shrinks about 4× (n5b_satisfying_weights_free.py). Here the fastest weight rate divided by ω_ref is about 0.14.
- The window is not the cause: centred windows of ±2 to ±40 periods all give 0.969–0.972.

**Fix.** Replace the sentence with: 'In the `satisfying` run the centred ratio of 0.97 is the error-channel factor (frozen-parameter ratio 0.951–0.956 along the run; 0.996 at k_e = 80) times a +1.4% averaging correction from the weights' oscillation within each period (fastest weight rate/ω_ref ≈ 0.14). With the weights running at frozen geometry the ratio is 0.969–0.970; with g_a/4 it is 0.956–0.960. The slope correction is below 0.01% here, because ∇L* is nearly orthogonal to the range of Da*ᵀ.'

#### minor: §6b Numerics table, row '`slope`, δ = 0.47' ('cosine 0.9994, ratio 1.10'; '(b′): cosine 0.9967, ratio 1.17'), and the bullets '(b′) tracks the closed loop at every S' and 'The size ratios. At small δ the 4% …'

**Claim.** The slope experiment compares (b′) like with like, and the size ratios are explained by the O(1/k_e) correction.

**Argument.** (b′) is stated 'after a boundary layer'. But exp_slow_manifold_slope.py starts the closed loop and (a) at a = a*(ξ0), off the slow manifold, and integrates (b′) from t = 0. The weights take 40–80 s of the 232 s horizon to reach the manifold. The δ = 0.47 row also has ε_avg ≈ 1/12.5 = 0.08. Neither effect is mentioned.

Decomposition (scratchpad\round3\propG\n3_decompose.py, which reproduces the table row):
- (a)/(b′) is 1.0615 when started at a*, but 1.0148 when started on the slow manifold. About 4.6% is the boundary layer.
- closed/(a) is 1.098 at full gains. Dividing all gains by 4 (R and S unchanged, ε_avg/4) gives 1.045–1.049; dividing by 16 gives 1.033; k_e = 80 gives 1.028.
- closed/(b′) goes from 1.166 (table) to 1.121 (start on the manifold), 1.065 (gains/4) and 1.048 (gains/16).
- So the 17% is roughly: 4.5% boundary layer, 6–7% averaging error, 3–4% error channel, and 1.5% intrinsic O(1/R) error of (b′).

**Fix.** In exp_slow_manifold_slope.py, start the closed loop and (a) at a(0) = h(ξ0) from the corrected (b′) offset, or compare centred velocities after the boundary layer; optionally also divide all gains by 4. Add to 'The size ratios': 'At δ = 0.47 the larger ratios (1.10 for (a), 1.17 for (b′)) come from the averaging error (ε_avg ≈ 0.08; dividing all gains by 4 at fixed R and S brings closed/(a) to 1.05). For (b′) they also come from the boundary layer included in the comparison: started on the slow manifold, (a)/(b′) = 1.015.'

#### minor: §6b, Proposition G(a), third bullet: 'Along the closed loop, its value at period-averaged parameters decreases, up to O(ε_avg) oscillation within a period.'

**Claim.** L̄ at period-averaged parameters decreases along the closed loop, apart from within-period oscillation.

**Argument.** At finite k_e the closed loop converges to an O(ω_eff/k_e)-shifted weight equilibrium, not to argmin_a L̄ (Proposition H.7). So near equilibrium L̄ rises systematically over many periods, which is not a within-period oscillation.

Counterexample (scratchpad\round3\propG\n4b_Lbar_increase.py): slope configuration, δ = 1e-5 (G ≻ 0; (H0) and (H2a′) hold). Start with a(0) = a*(ξ0) − 0.01232, so that the first period mean equals a*(ξ0). L̄ at period means then rises in 35 consecutive periods, from 2.12435e-3 to 2.12516e-3 (+8.1e-7 = ½m₀·0.008²). At k_e = 80 the rise is 3.7e-8, about 16× smaller, i.e. O(1/k_e²).

**Fix.** Replace with: 'Along the closed loop, L̄ at period-averaged parameters decreases while the averaged dissipation (m/k_e)‖∇L̄‖²_G dominates the O(ω_eff/k_e) non-gradient correction. Near equilibrium it can increase by O((ω_eff/k_e)²) as the weights settle at their shifted equilibrium (Proposition H.7).'

#### minor: §6b Derivation, step 1: 'and e = r(x_m)/k_e + O(ω_eff·sup‖r‖/k_e²)'

**Claim.** Under (H2a′), the error channel satisfies e = r(x_m)/k_e up to O(ω_eff·sup‖r‖/k_e²).

**Argument.** ω_eff contains the Lipschitz constant of f − f̂ but not its spatial bandwidth. So ṙ ≈ ∂_x r·ẋ_m can be much larger than ω_eff·sup‖r‖.

Counterexample (scratchpad\round3\propG\n7_step1.py): n = 1, f̂ = 0, f(x) = 0.01 sin(100x), x_m = sin t, s_min = 0.25. Then ω_eff = 4, k_e = 100 (k_e/ω_eff = 25) and sup‖r‖/k_e = 1e-4 ≪ s_min, so (H2a′) holds. Yet sup|e − r(x_m)/k_e| = 7.1e-5, against a claimed scale of 4e-6 (17.8× larger), and the RMS remainder is 54% of RMS r/k_e.

Conclusion (a) is not affected here, because the remainder is high-frequency and averages out against the band-limited ∇_ϑΦ(x_m(t)).

**Fix.** Write: 'So ‖w(t)‖ ≤ e^{−k_e t}‖w(0)‖ + sup‖ṙ‖/k_e², and e = r(x_m)/k_e + O(sup‖ṙ‖/k_e² + sup‖∂_x r‖·sup‖r‖/k_e²).' Alternatively, add the bandwidth sup_t‖(d/dt)r(x_m(t))‖/sup‖r‖ to ω_eff in (H2a′).

#### minor: §6b Numerics table, row '`satisfying`', G(a) cell 'geometry block cosine 0.999 / 0.998 (see below)'; bullet 'The weight block … Their net drift is small: the predicted weight speed is 0.17–0.28 of the geometry speed in the G⁻¹ norm.'

**Claim.** The weights' net drift is small, as measured by the predicted weight speed of 0.17–0.28; G(a) is reported through the geometry block only.

**Argument.** - The table cell is selective. table_satisfying.md gives a joint G(a) cosine of 0.9848 / 0.9624 and ratio 0.956 / 0.935, and 0% of samples have cosine > 0.99. None of this appears in §6b.
- The predicted weight speed does not measure the net drift. Measured with centred differences in the G⁻¹ norm, the actual weight speed is 0.001–0.008 of the geometry speed, while the predicted one is 0.167–0.300: the prediction is 40–200× the actual drift (scratchpad\round3\propG\n5b_satisfying_weights_free.py).
- The prediction −(m/k_e)g_a∇_aL̄ at the period-mean weights is almost entirely the restoring velocity of the O(1/k_e) offset, which the closed loop balances with its own O(1/k_e) correction. The bullet's conclusion is right, but its supporting sentence misreads what the 0.17–0.28 measures.

**Fix.** Table cell: 'joint cosine 0.985 / 0.962, size ratio 0.96 / 0.94; geometry block cosine 0.999 / 0.998; weight block fails (see below)'. Bullet: 'The weights sit at a* + O(1/k_e). Their actual drift is tiny: 0.001–0.008 of the geometry speed in the G⁻¹ norm. The prediction −(m/k_e)g_a∇_aL̄ at the period-mean weights is 0.17–0.30 of the geometry speed, i.e. 40–200 times the actual drift. It is the restoring velocity of the O(1/k_e) offset, which (a) does not contain. So the weight block fails (cosine −0.21 / 0.50, 92–95% of the residual), while the geometry block matches.'

#### minor: §6b Status, 'What the designer supplies': 'g_c, g_s ≪ g_a·min(m₀/‖Γ∇²L*‖, 1/‖Da*Γ^{1/2}‖²)'. The same text is in LITERATURE_TRANSPORT.md §G.

**Claim.** This is the gain condition that combines (H2b) and (H2c).

**Argument.** Since g_c = δγ_c and g_s = δγ_s, the condition counts Γ twice. (H2b) and (H2c) give δ‖Γ^{1/2}∇²L*Γ^{1/2}‖ ≪ g_a m₀ and δ‖Da*Γ^{1/2}‖² ≪ g_a. With Γ = γI the Status line reads δγ² ≪ g_a·min(m₀/‖∇²L*‖, 1/‖Da*‖²), which is off by a factor γ. It is dimensionally consistent only when Γ = I.

**Fix.** Replace with: 'δ ≪ g_a·min(m₀/‖Γ^{1/2}∇²L*Γ^{1/2}‖, 1/‖Da*Γ^{1/2}‖²); equivalently, with G_ξ := δΓ = diag(g_c I_n, g_s): ‖G_ξ^{1/2}∇²L*G_ξ^{1/2}‖ ≪ g_a m₀ and ‖Da*G_ξ^{1/2}‖² ≪ g_a.' Make the same change in LITERATURE_TRANSPORT.md §G.

#### cosmetic: §6b (b′): 'The geometry follows, to first order,'; Numerics: '(b′) is the correct first-order slow flow'; §8 item 4: 'Open: … higher-order corrections in S'

**Claim.** (b′) is a first-order slow flow, and corrections of higher order in S remain open.

**Argument.** (b′) is leading order in 1/R and exact in S. Dropping Dη and the quadratic term, the invariance equation is solved exactly by the resummed form. Numerically, (b′) matches the joint flow started on the slow manifold within 1.5%, 0.6% and 0.2% at S = 2, 8.5 and 17 (R = 178, 42, 21). The strictly first-order-in-S truncation has cosine −0.44, −0.25 and −0.10 with the joint flow and increases L* (scratchpad\round3\propG\n13_large_S.py). So 'first order' invites the wrong truncation, and the remaining corrections are in 1/R, not in S.

**Fix.** Write 'The geometry follows, to leading order in 1/R and for every S,' and '(b′) is the correct leading-order slow flow'. In §8, replace 'higher-order corrections in S' with 'the O(1/R) corrections (they involve D²a* and ∂_ξM)'.

#### cosmetic: §6b remark '(H2b) does not imply (H2c). The product R·S = m₀‖Da*Γ^{1/2}‖²/‖Γ^{1/2}∇²L*Γ^{1/2}‖ is fixed by the problem.'

**Claim.** R·S is fixed by the problem.

**Argument.** R·S does not depend on g_a, δ or k_e, but it does depend on the gain shape Γ. In the slope configuration at ξ0 it is 166, 355 and 1149 for γ_s/γ_c = 1e-4, 1 and 1e4 (scratchpad\round3\propG\n8_RS_gamma.py). The conclusion still holds, since R·S ≥ 166 for every Γ tested.

**Fix.** Write: 'The product R·S = … does not depend on g_a, δ or k_e. It is fixed by the problem and the gain shape Γ: in the slope configuration it ranges from 166 to 1149 over γ_s/γ_c ∈ [10⁻⁴, 10⁴].'

#### cosmetic: §6b: 'O(ε_avg)' in (a); '(ω_eff/k_e, ε_avg, S, 1/R) → 0' in (b); '[0, b/ε_avg]' in step 2; '𝒢′ ⋐ 𝒢' in (H0); 'M ⪰ m₀I on 𝒢' in (H1); 'R := g_a m₀/(δ‖Γ^{1/2}∇²L*Γ^{1/2}‖)'

**Claim.** These symbols are defined.

**Argument.** ε_avg is not defined anywhere in the document, yet (b)'s limit treats it as an independent parameter. 𝒢 is not defined in §6b, and §6c uses 𝒢(μ) for a functional. R is defined pointwise, but its use in the limit needs a sup over 𝒢′.

**Fix.** Add to the Setting: 'ε_avg := (m/k_e)λ_max(G^{1/2}∇²L̄G^{1/2})/ω_eff, the second ratio in (H2a′); the tables report its inverse with ω_ref in place of ω_eff.' Also add: '𝒢 denotes the set of geometries ξ on which (H1) holds.' Rename one of the two 𝒢 symbols. Define R with sup_{𝒢′}‖Γ^{1/2}∇²L*Γ^{1/2}‖.

#### cosmetic: §6b Numerics heading '(`results/variable_projection`; centred differences of period means)', the column header 'G(b) and the induced-metric flow (b′)', and the (H2a) columns of the tables

**Claim.** All rows use centred differences, (b′) is reported where the column says so, and (H2a′) is checked.

**Argument.** - The slope rows use end-to-end displacements between the first and last period means, not centred differences.
- (b′) is computed only in the slope rows, not in `violating` or `satisfying`.
- The tables report only the second inequality of (H2a′) (ω_ref over the joint rate). The first, k_e ≫ ω_eff, gives k_e/ω_eff = 20/(2.5/0.25) = 2.0 in the m = 2 and m = 4 runs, so it holds only marginally there; with the actual widths the ratio is about 4.8 (scratchpad\round3\propG\n12_omega_eff.py).

**Fix.** Change the heading to: '(centred differences of period means in the first two rows; displacement between the first and last period means in the slope rows)'. Write '(b′) was evaluated only in the slope rows'. Report k_e/ω_eff: 2.0 for `violating` and `satisfying`, about 20 for `slope`.

## Refuter: AF-regression

### Holds

- Tests: all 16 tests in tests/test_transport.py pass (306 s, project venv, OPENBLAS_NUM_THREADS=1, no cache). Every test name cited in the §7 rows for A–F exists.
- Corollary B.5 and its (L_σ) variant. The derivation from C(iii) is correct. Under 4η̄(b+η̄) < Δ_gap, C(ii) gives C = b(t)² ≤ b². Then max_η(−σ_mη² + σ_mηb) = σ_m b²/4, and integrating and completing the square gives the displayed RMS bound. The 'without the gap condition' clause is valid: C ≤ b² + 4η̄R gives −σ_mη² + σ_mη√C ≤ σ_mC/4. Reproduced test instance (b = 0.01): RMS 1.816e-3, plain B.5 bound 5.20e-4 (exceeded 3.49×), variant 2.59e-3. RMS is linear in b: 3.63e-4, 1.82e-3, 9.08e-3 for b = 0.002, 0.01, 0.05. The plain-(L) limit ε_max/k_e is also correct.
- Remark B.4 numbers recomputed by hand (J = 3.1784, L₁ = 25.447, Δ = 14.373): classical constant on the per-kernel set 214.37 for every m, crude OT form 533.25, block split 382.03, per particle 214.37. With a joint centre radius c_max, v2's figures reproduce: 137.8 at m = 6 and 101.15 at m = 64.
- Lemma P.1: the exact identity matches TransportController.project to 7e-16 on 80,000 particle-states (anisotropic G, both s-faces), and the largest ratio is 1.0000000. The Gaussian suprema 0.9165/s_min and 1.573/s_min² are confirmed. The (P) argument holds for both the ball factors and the box factor.
- Theorem A: the scaling identity, the gain map, the two-case only-if construction (each case checked) and the m = 1 coincidence hold. The classical-leak algebra D = diag(I/m, mI, m) with λ = σ·min(g_a/m, m g_c, m g_s) holds. The projection-row Jacobian asymmetry is reproduced: −0.0598 at a_max = 1, ε_K = 0.15.
- Theorem B, Lemma B.1/B.2 and Corollary B.3 re-derived. Corollary B.6 holds: the coupling bound, the block constant, the centroid second-order bound, the fact that a best equal-weight approximation is a centroid quantiser, the Maurey form via Hilbert-space concentration, and the cited results within their hypotheses (Chevallier for p ≥ 3, Fournier–Guillin for p > 4, Zador lower bound). The Remark B.7 coupling sketch is consistent.
- Theorem C. The stationary Filippov example at the tie holds (e ≡ 0 is consistent), as do the polarisation identity for the Filippov bound and its tightness (0 = −σ_m/4 + σ_m/4). Also verified: (i) for selection solutions; both bounds of (ii) and the gap clause; both displays of (iii); the Young splits, α, β and invariance in (iv). Scope: the antisymmetry argument gives Δ_gap ≤ 2R²/m, hence η̄ < √Δ_gap/2 ≤ R/√(2m) and b < η̄/√2, and switching forces Δ_gap = 0. The statement of (v) and its 0.5 vs 1.0 example hold. Lemma C.3: the affine-in-λ argument and the strictness at ties hold. Proposition C.5: the swap is optimal iff δ < ½, D⁺V_OT = σ_m(δ−δ²), and a₁+a₂ is conserved. Transient numbers: 2.2585/0.3917 = 5.8×, and 0.0771/0.0743 is +3.8%. Remark C.6 wording is correct.
- Theorem D. In (b), the first-exit argument and the linear-growth formula for F hold. For the Van der Pol plant, min_F S(F) − F is attained at F = 0 with S(0) ≈ 51.56. On the experiment's box L_2 = λ_G^{1/2}L₁ = 33.077 and L_2·15 = 496.1. Also holding: Grönwall on S = A + D in (c), the bootstrap in (d), the integrated form in (e), and the lemma of Remark D.1, re-proved (case h' = 0 gives ≥ 0; when both are active, √(rr') ≥ θ_max − ε_K). The numerics (1.68/1.64/1.69/1.68; 4.6e-3 down to below 1e-4) match results/mean_field_m/table.md.
- Theorem F: (i)–(iii), (i′)–(iii′), d_lab ≥ W, L_G ≤ λ_G^{1/2}L₁, the W_∞ bound and the Reading all hold. The e ≡ 0 (L_σ) run reproduces: sup‖e‖ = 7e-14, W(μ_10, μ_0) = 1.0535, leak-clause ratio 0.999998.
- Theorem E. Verified: the Jacobian A(t) and QA + AᵀQ = diag(−k_eI, 0); the ζ change of variables; the PE transfer; the richness bound; the closed form of d² and the m = 1 convention; the c₀ and C norm conversions and the strictness; uniqueness via the two Khalil balls; the limits. Remark E.3 holds (V_OT = d²/8, two optimal assignments, W ≥ d/2). Corollary E.2: the merging bound and the Hölder algebra hold, and both pieces equal (CR₀)^{Λ/(Λ+λ)}δ₀^{λ/(Λ+λ)}. The sin t² ratios reproduce from w(t₀) = (0,1). The Remark E.1 table ranges and the rank correlation 0.60 match results/pe_convergence/table.md.
- These round-2 fixes are applied correctly: AB.1, AB.3–AB.14, C.1–C.3, C.5–C.10, D.1–D.9, E.1, E.3–E.5 and E.7–E.11, plus the PARADIGM.md parts of AB.2, F.1 and F.2–F.9. C.4 is applied except for the 'Strict case' sentence. E.2 and E.6 are applied but introduce the contradictions reported below. The REVIEW_LOG parts of AB.2 and F.1 are not applied.

### Issues

#### minor: §4, proof of Theorem C(v), bullet '*Strict case.* For a moving ideal whose optimal ideal–prior assignment changes, min_π sup_t b_lab(π)(t) > b.'

**Claim.** Whenever the optimal ideal–prior assignment changes along a moving ideal, the best fixed-label bias strictly exceeds b.

**Argument.** This is false in general. The sentence turns the round-2 judge's phrase 'strictly whenever the optimal assignment changes' into a theorem. The statement of (v) itself, '≥ b, with equality when one assignment is optimal at every t', is correct.

Counterexample (C^∞ and admissible):
- Setting: m = 2, n = 1, G = I, prior a ∈ {0,1}, other coordinates equal, a_max − ε_K ≥ 1.6.
- Ideal: a*(t) = (u(t), 1 − u(t)) with u(t) = −0.025 + 0.575 sin t.
- Costs: c_id = u² and c_swap = (1 − u)².
- The swap is strictly optimal when u > ½ (e.g. t = π/2, u = 0.55), and the identity when u < ½. So the optimal assignment switches.
- b² = sup_t min(u², (1−u)²) = 0.36, attained at u = −0.6.
- sup_t b_lab(id)(t)² = sup u² = 0.36.
- Hence min_π sup_t b_lab(π)(t) = 0.6 = b.

Scripts: scratchpad/round3/fixes_1to6/cv_strict_case_smooth.py and cv_strict_case.py. The latter also reproduces the document's 0.5 vs 1.0.

**Fix.** Replace the bullet by: '*Strict case.* If the optimal ideal–prior assignment changes along the ideal, the inequality can be strict: a*_1(t) = t and a*_2(t) = 1 − t on [0,1], with prior a ∈ {0, 1}, give b = 0.5 against 1.0. It need not be strict: for a*(t) = (u(t), 1 − u(t)) with u(t) = −0.025 + 0.575 sin t the optimal assignment alternates between the identity and the swap, yet min_π sup_t b_lab(π)(t) = sup_t b_lab(id)(t) = 0.6 = b.'

#### minor: §6 Remark E.1, bullet 'The rate is not set by the PE level.', sentence 'λ is bounded below by a function of the PE level αλ_min(G)/m, k_e, T and sup‖Ψ_G‖, but the PE level does not determine it.' The same claim is in the comment of tests/test_transport.py::test_theorem_E_matched_error_converges_labelled_does_not, line 435: 'The rate is bounded below by a function of the PE level'.

**Claim.** Some positive function of (αλ_min(G)/m, k_e, T, sup‖Ψ_G‖) bounds the rate λ of Theorem E from below.

**Argument.** The claim is false, and it contradicts the document twice:
- Proof step 3 says 'The bound on Ψ̇ cannot be dropped'.
- The paragraph after Theorem E says that uniformity of k and λ over the class with given (T, α, sup‖Ψ‖, sup‖Ψ̇‖) is not proved. The remark asserts this uniformity, and without sup‖Ψ̇‖.

Counterexample inside Theorem E's own setting (Gaussian feature, sinusoidal reference):
- Setting: m = n = 1, ϑ* = (1.5, 0.3, 0.8), G = I, k_e = 3, x_m(t) = 1.5 sin(ωt) with integer ω, T = 2π.
- Every window of length T contains ω whole periods. So λ_min(∫_t^{t+T}Ψ_GᵀΨ_G) = 1.1329 and sup‖Ψ_G‖ = 1.779 for every ω, and k_e and T are fixed.
- The exact Floquet rate of the step-3 linearisation is 4.66e-2, 2.16e-2, 6.83e-3, 1.83e-3, 4.66e-4, 1.17e-4 and 2.93e-5 for ω = 1, 2, 4, 8, 16, 32 and 64. It falls like ω⁻².

So no positive function of these four quantities bounds λ below; only sup‖Ψ̇‖ (∝ ω) changes. Script: e1_rate_vs_omega.py. It also reproduces round 2's rate of 4.66e-2 at ω = 1.

**Fix.** Replace the sentence by: 'λ > 0 exists for each admissible Ψ, but it depends on Ψ(·) and not only on its PE level: no positive function of αλ_min(G)/m, k_e, T and sup‖Ψ_G‖ alone bounds it below. For m = n = 1, ϑ* = (1.5, 0.3, 0.8), G = I, k_e = 3 and x_m = 1.5 sin ωt (T = 2π) these four quantities do not depend on ω, while the exact Floquet rate falls from 4.7e-2 (ω = 1) to 2.9e-5 (ω = 64), like ω⁻². A lower bound must also involve sup‖Ψ̇_G‖ (step 3), and none is proved here.'

In the test comment, replace 'The rate is bounded below by a function of the PE level (it is not set by it, and' with 'The rate depends on Ψ(·), not only on the PE level (it is not set by it, and'.

#### minor: §6 Corollary E.2: 'Combine this with Theorem D(c) on [0, T₀], with T₀ := (Λ + λ)^{−1}log(CR₀/δ₀). Both solutions stay in the closed ρ-tube, which is compact. This gives uniform-in-time Hölder continuity in the basin'

**Claim.** Theorem D(c), applied with X equal to the closed ρ-tube, supplies a Λ for which T₀ and the exponent λ/(Λ+λ) are well defined and uniform in time.

**Argument.** D(c) sets M_x := sup_{t≤T} Lip_X f(·,t), the Lipschitz constant of f(·,t) on all of X = ∪_s B̄(x_m(s), ρ). This creates two problems.

First, Theorem E assumes nothing about f outside the moving tube {‖x − x_m(t)‖ < ρ}; (S1) is not among its assumptions. X contains points at distance > ρ from x_m(t) whenever x_m moves. For example, ε(x,t) = ((‖x − x_m(t)‖ − ρ)₊)^{1/2}·v satisfies E's tube hypothesis but gives Lip_X f(·,t) = ∞.

Second, even under (S1) the constant is uniform only on bounded time intervals. Then Λ = Λ(T₀), the definition T₀ := (Λ+λ)^{−1}log(CR₀/δ₀) is implicit, and the 'uniform-in-time' exponent can degenerate as δ₀ → 0.

The result is true. The coupling estimate only uses the Lipschitz constant on the segment [x_t, x'_t] ⊂ B(x_m(t), r), where f = f̂* and the projection is inactive.

**Fix.** Replace the two sentences by: 'Combine this with the coupling estimate of Theorem D(c) on [0, T₀], T₀ := (Λ + λ)^{−1}log(CR₀/δ₀). Both solutions stay in D (proof step 4), so ‖x − x_m(t)‖ < r ≤ ρ, the segment between x_t and x'_t lies in B(x_m(t), ρ) where f = f̂*, and Proj_K = id. Hence (c) holds with M_x := k_e + 2a_max J_c, with L_x and L_ϑ computed for Ē = r and an inactive projection, and with L_2 = λ_G^{1/2}L₁; this Λ depends neither on T₀ nor on f outside the tube.'

#### minor: Header, 'Companion files': 'Review log. review/REVIEW_LOG.md lists every confirmed issue and its fix.' Also round-2 fixes AB.2 ('Also correct the 3.9–5.3× row in review/REVIEW_LOG.md') and F.1 ('Log the erratum in review/REVIEW_LOG.md').

**Claim.** The review log lists every confirmed issue and its fix.

**Argument.** review/REVIEW_LOG.md is dated 14 Oct 21:56, before REVIEW_round2.md (23:41), and has no round-2 section:
- It still says round 2 'produced no verdicts' and 'Pending: re-run round 2 on C.3, C.5, B.6, F, G, E.2 and H'.
- Row B.6 still says 'OT constant is 3.9–5.3× the classical one'. v3's Remark B.4 shows both constants are 214 on the same per-kernel set; the crude OT form is 533/214 = 2.5×.
- The erratum for v2's contraction identity, which F.1 asked to log (missing positive part; the identity is per factor), is absent.

So the log parts of round-2 fixes AB.2 and F.1 were not applied, and the header statement is false.

**Fix.** Add a 'Round 2' section to review/REVIEW_LOG.md, with one row per confirmed issue (AB.1–AB.14, C.1–C.10, D.1–D.9, E.1–E.11, F.1–F.9, and the G and H issues) and its fix.
- In row B.6, replace '3.9–5.3× the classical one' by 'on the same per-kernel set the classical constant is 214 for every m; the OT linear branch is 533 (crude), 382 (block split through W) and 214 (per particle)'.
- Add the erratum: 'the v2 identity ‖Proj(ϑ,Gy)‖² = ‖y‖²_G − h(2−h)(ϑᵀGy)²/(ϑᵀGϑ) lacks the positive part and holds only per factor; replaced by Lemma P.1'.
- Delete the 'Pending' line.

Alternatively, change the header to 'lists the round-1 issues; round-2 issues and fixes are in review/REVIEW_round2.md'.

#### minor: §0 table, row E, status 'proved; linearisation checked; convergence observed'. Compare §6 proof step 3, last bullet: 'We could not check the exact richness definition and regularity class of [S8] against its full text.'

**Claim.** Theorem E is proved.

**Argument.** The only route to uniform exponential stability (UES) of the linearisation is [S8, Thm 1]. The document admits it matched the hypotheses of that theorem against the abstract only.

The regularity class is exactly what decides UES here:
- Step 3's own sin t² example is PE and bounded but not UES.
- The ω-sweep of the Remark E.1 issue shows the rate is governed by Ψ̇.

An unqualified 'proved' therefore overstates the status. The proof text itself is transparent about the gap.

**Fix.** Change the status to 'proved, conditional on [S8, Thm 1] as stated in its abstract (hypotheses not checked against the full text); linearisation checked; convergence observed'.

Alternatively, give a self-contained step 3:
- V = ½‖w‖² gives V̇ = −k_e‖e‖².
- UES then follows from uniform complete observability of ([I 0], A(t)).
- By output injection with the bounded K = [k_eI; −Ψ_Gᵀ], this is equivalent to observability of ([I 0], [[0, −Ψ_G], [0, 0]]), i.e. uniform positivity of ∫_t^{t+T}[I, −∫_t^sΨ_G]ᵀ[I, −∫_t^sΨ_G]ds.
- (PE) together with the bound on Ψ̇_G implies this positivity.

Cite the observability-implies-UES lemma used, with its theorem number, after checking it.

#### cosmetic: §5b 'Checked at every sample … (ii′) along (L_σ)'; §0 row F 'proved and checked, including the leak clause'; §7 row 'F (ii′) along (L_σ)'; tests/test_transport.py::test_theorem_F_leak_clause_along_ot_sigma_mod

**Claim.** The leak clause (ii′) is checked by the test.

**Argument.** The test cannot detect an error in the leak term:
- In the test's own run, the worst ratio W(μ_{t₁},μ_{t₀})/bound is 0.237.
- With the leak term σ_mW(μ_r,μ⁰) deleted from the bound, the worst ratio is 0.308, so the test still passes.
- The leak contributes 9% of the integral, and the assignment never switches.
(Script: f_leak_test_power.py.)

A discriminating instance exists: round 2's run with e ≡ 0.
- Setting: n = 1, m = 2, the test box and gains, σ_m = 0.5; ideal two copies of (1,0,1); prior a⁰ = (1.5, 0.5); start a = (3, −1) with (c,s) = (0,1); x(0) = x_m(0).
- Result: sup‖e‖ = 7e-14, W(μ_10,μ_0) = 1.0535, worst ratio 0.999998.
(Script: f_tight_leak_instance.py.)

Also, 'at every sample' is inaccurate for (ii′), which the test checks on a grid of pairs.

**Fix.** Add the e ≡ 0 instance to the test, asserting:
- sup‖e‖ < 1e-10;
- W(μ_10, μ_0) > 1.05;
- 0.99 ≤ max W(μ_{t₁},μ_{t₀}) / ∫(λ_G^{1/2}L₁‖e‖ + σ_mW(μ_r,μ⁰))dr ≤ 1 + 1e-3.

In §5b, write '(ii′) along (L_σ), on a grid of pairs, including an instance with e ≡ 0 where it is tight (ratio ≥ 0.99)'.

#### cosmetic: §1 (S5): 'For bounded f, take X to be the closure of ∪_t 𝒯(t) from Theorem D(b).'

**Claim.** (S5) can be discharged with the tube 𝒯(t) of Theorem D(b).

**Argument.** 𝒯 was the tube symbol of v2's D(b). The rewritten D(b) (round-2 fix D.2) defines only r(F) and S(F). The symbol 𝒯 is now undefined anywhere in v3 (grep finds only this line).

**Fix.** Write: 'For bounded f, take X := ∪_t{x : ‖x − x_m(t)‖ ≤ r(F)} with F = sup‖f‖ (Theorem D(b)). Then ε_max is the supremum over this X, and it depends on ‖e(0)‖.'

#### cosmetic: §6, after Theorem E: 'The constants depend on m, G, k_e, T, α, sup‖Ψ‖, sup‖Ψ̇‖, a_max, s_min, ρ and δ. The cited theorems are of existence type, so it is not proved here that k and λ are uniform over the class of Ψ with given (T, α, sup‖Ψ‖, sup‖Ψ̇‖).'

**Claim.** The constants are functions of the listed quantities.

**Argument.** Read as 'depend only on', the first sentence is exactly the uniformity that the second sentence disclaims. k and λ depend on the function Ψ(·) itself, through the converse Lyapunov function of the linearisation. The ω-sweep in the Remark E.1 issue shows this dependence is real.

**Fix.** 'k and λ depend on Ψ(·) itself (through the converse Lyapunov function of its linearisation), on m, G and k_e; r₀ and c₀ depend in addition on a_max, s_min, ρ and δ. It is not proved here that k and λ can be chosen uniformly over the class of Ψ with given (T, α, sup‖Ψ‖, sup‖Ψ̇‖).'

#### cosmetic: §3, after the Corollary B.5 variant: 'As η̄ → 0 the variant tends to ε_max/(2k_e) + (ε_max²/(4k_e²) + σ_m b²/(4k_e))^{1/2}'. Also the proof: 'Without the gap condition, replace b² by b² + 4η̄R.'

**Claim.** This is the limit of the (L_σ) bound, and the no-gap replacement is b² + 4η̄R.

**Argument.** (a) The limit also needs T₂ → ∞, because of the V_OT(T₁)/(k_e(T₂−T₁)) term.

(b) Along (L_σ) with b > 0, η̄ cannot tend to 0, because the leaky equilibrium is not ϑ*. In the test instance, W stays in [0.81b, 0.88b] on [20, 60] for b = 0.002, 0.01 and 0.05, and C(iv) forces η̄² > 2β/α ≥ 2b². So the limit is formal. The conclusion 'first order in b' still holds, since D̄ = O(η̄²).

(c) The no-gap clause is valid but depends on R. The R-free bound √C ≤ b + 2η of C(ii) gives σ_mC/4 ≤ σ_m(b+2η̄)²/4. In the test instance this bound is 9.8e-3 against 0.129, with measured RMS 1.8e-3. Script: b5_numbers.py.

**Fix.** Write: 'As η̄ → 0 and T₂ → ∞ the variant tends formally to …; along (L_σ) with b > 0 the measured η̄ stays of order b (about 0.9b in the test instance), so the bound is first order in b.'

In the proof, write: 'Without the gap condition, replace b² by (b + 2η̄)² (from √C ≤ b + 2η in C(ii)), or by b² + 4η̄R.'

#### cosmetic: §2, proof of (iii), 'Projection row': 'In the example m = 1, n = 1, G = I, x = 0.7, e = 1 and ϑ = (0.93, 0.2, 1.5), with a in the layer: the c-derivative of the a-component is −0.060'

**Claim.** This data determines the Jacobian entry −0.060.

**Argument.** The value depends on (a_max, ε_K), which are not given. With the project code it is:
- −0.0598 for (a_max, ε_K) = (1, 0.15);
- −0.0189 for (1, 0.1);
- −0.1236 for (1, 0.3);
- −0.1345 for (0.95, 0.1);
- 0 when a is not in the layer (a_max = 2).
The qualitative conclusion, an asymmetric Jacobian, holds for every box with a in the layer. Script: a_projection_row.py.

**Fix.** Insert 'a_max = 1, ε_K = 0.15 (so h(0.93) = 0.284) and c_max = 2' after 'ϑ = (0.93, 0.2, 1.5)'.

#### cosmetic: §2, proof of (i): 'the only-if construction against ProjectionOperator.project_matrix on five active joint balls'. §7 row 'A(i) projected rates, only-if': 'five active joint balls are each distinguished'.

**Claim.** The test distinguishes five active joint balls.

**Argument.** In the test, (W_max, ε_W) = (3.0, 0.5) falls in the second case: ρ = 2.5 ≥ √2·r_in = 1.20. There the joint ball is inactive and the per-kernel field is active. Only four of the five joint balls are active (script a_only_if_activity.py).

**Fix.** Write 'on five joint balls (four active against an inactive per-kernel field, one inactive against an active per-kernel field)', in both places.

#### cosmetic: §7 verification map: row 'B.5 along (L) and (L_σ)'; row 'C(i), Lemma C.3'; row 'C(iv)' ('an instance satisfying every hypothesis, asserted'); the list 'Not checked in the repository'.

**Claim.** The map states accurately what each test checks and what is not checked.

**Argument.** (a) No test checks B.5 along (L). The cited test runs (L_σ) only.

(b) test_theorem_C_i_exact_dissipation_pure_leak integrates raw_velocity, i.e. the unprojected flow of Lemma C.3. The claim that Proj_K is inactive along the pure leak is not tested.

(c) test_theorem_C_iv sets η̄ = 0.9·Δ_gap/(4R), the v2 condition. It does not assert v3's condition 4η̄(b+η̄) < Δ_gap. That condition does hold: 4·0.267·(0.048+0.267) = 0.34 < 5.44.

(d) The 'Not checked' list omits several items:
- B.5 along (L);
- B.6 and B.7;
- the Filippov form of C(i);
- the inequality of C(v);
- Proposition C.5's D⁺V_OT = σ_m(δ − δ²);
- D(b);
- Corollary E.2 and Remark E.3.

**Fix.** - Rename the B.5 row 'B.5: failure of the plain bound along (L_σ); the (L_σ) variant holds'.
- In the C(i) row, write 'Lemma C.3 (unprojected pure leak)'.
- Assert 4η̄(b+η̄) < Δ_gap in the C(iv) test.
- Extend 'Not checked in the repository' with the items in (d).

#### cosmetic: §6, proof step 3: 'the ratio ‖w(t₀+200)‖/‖w(t₀)‖ is 0.9715, 0.9942 and 0.9980 for t₀ = 5, 20 and 50'

**Claim.** The ratio takes these values.

**Argument.** The ratio depends on the initial state w(t₀):
- from w(t₀) = (0, 1) I reproduce 0.9715, 0.9942 and 0.9980;
- from (1, 0) it is 0.094, 0.014 and 0.0075;
- from (1, 1)/√2 it is 0.753, 0.693 and 0.711.
Showing non-UES needs the supremum over initial states, which (0, 1) supplies. Script: e1_rate_vs_omega.py.

**Fix.** Write '… from w(t₀) = (e, ζ) = (0, 1), the ratio ‖w(t₀+200)‖/‖w(t₀)‖ is …'.

#### cosmetic: §5, proof of Theorem D(b): 'Variation of constants gives ‖e(t)‖ ≤ max(‖e(0)‖, (F + a_max)/k_e) < r(F).'

**Claim.** max(‖e(0)‖, (F + a_max)/k_e) < r(F).

**Argument.** If e(0) = 0, the maximum equals (F + a_max)/k_e = r(F), so the strict inequality is false. The first-exit argument still works with the convex-combination bound e^{−k_e t}‖e(0)‖ + (1 − e^{−k_e t})(F + a_max)/k_e, which is < r(F) for every finite t because a_max > 0.

**Fix.** 'Variation of constants gives ‖e(t)‖ ≤ e^{−k_e t}‖e(0)‖ + (1 − e^{−k_e t})(F + a_max)/k_e < r(F) for every finite t, since a_max > 0.'

#### cosmetic: §5, Theorem D(c), bullet L_2: 'For any coupling, ‖∫Φ(x;·)d(μ − μ')‖ ≤ L_2 times the coupling's L²(G⁻¹) cost.'

**Claim.** The bound is L_2 times the transport cost of the coupling.

**Argument.** The transport cost is ∫‖ϑ−ϑ'‖²_{G⁻¹}dγ. The inequality holds with its square root, the L² norm, not with the cost itself: the two sides have different homogeneity.

**Fix.** 'For any coupling γ, ‖∫Φ(x;·)d(μ − μ')‖ ≤ L_2(∫‖ϑ − ϑ'‖²_{G⁻¹}dγ)^{1/2}.'

#### cosmetic: §2, A(iii), dictionary row 'classical σ-modification Γ(y − σθ) (leak inside the gain, src/baselines/sigma_mod.py) | … D = diag(I_n/m, mI_n, m)'

**Claim.** The cited file implements the leak whose measure-level object is the D-weighted energy.

**Argument.** src/baselines/sigma_mod.py adapts only W, with Ẇ = Γ_W(φeᵀ − σW). Its own comment says 'c and sigma are not adapted', and sigma_dot is zero. The object for that file is (σ/2)∫(a − a⁰)ᵀ(I_n/m)(a − a⁰)dμ with (c,s) frozen. The D = diag(I_n/m, mI_n, m) object describes a leak on all three blocks, which no file in the repository implements.

**Fix.** Write 'classical σ-modification Γ(y − σθ) (leak inside the gain; src/baselines/sigma_mod.py implements its W-block with c and σ frozen)'. Keep the D-weighted object, stated for the full-vector leak.

#### cosmetic: §4, Theorem C(i): 'For a Filippov solution with weights θ_k on tied optimal assignments τ_k, the bound weakens to D⁺Z ≤ −2σ_mZ + (σ_m/2)Σ_kθ_kA_{ττ_k} + (2Z)^{1/2}Γ.'

**Claim.** This is the Filippov bound.

**Argument.** τ is not specified in the statement. The proof yields the bound for any τ optimal at t, which enters through Z ≤ Z_τ. So the minimum over optimal τ may be taken.

**Fix.** '… ≤ −2σ_mZ + (σ_m/2)min_τ Σ_kθ_kA_{ττ_k} + (2Z)^{1/2}Γ, the minimum over assignments τ optimal at t.'

#### cosmetic: §3, Remark B.4: 'Version 2's comparison figure of "101–138" used a joint centre radius c_max, a set √m times smaller than the per-kernel product of Theorem A.'

**Claim.** The joint centre set is √m times smaller.

**Argument.** It is the radius that is √m times smaller: c_max against the radius √m·c_max of the ball that circumscribes the per-kernel product. The volume ratio is not √m. Recomputed: 137.8 (m = 6) and 101.15 (m = 64) with joint radius c_max, and 214.37 for every m with radius √m·c_max.

**Fix.** '… used a joint centre ball of radius c_max, √m times smaller than the radius √m·c_max of the ball that circumscribes the per-kernel product of Theorem A.'

