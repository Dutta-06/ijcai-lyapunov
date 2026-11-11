# Referee D3az, round 53: Proposition D‴.54 (the transit in ℝⁿ; (v) rewritten as a race)

**Summary.** The fixes to (i)–(iv) and to the bookkeeping hold:
- the md regenerates byte for byte from the raw;
- all 22 tests pass;
- every number in the numerics paragraph matches the md, except the ones computed by hand.

The projected second-moment law in (v) is right, and the race is real. An ε-free reduced problem (C4's two kernels after feature 0's entry, with features 1 and 2 left in the residual, integrated in s) reproduces C4's departure times, and predicts all 15 of C3's outcomes and 28 of C4's 30. But (v) states the hold as a general law, and it holds only on one side:

- **The far side.** The proof says E*(w) increases with w "along every branch computed (V)". V computes only the branch toward the unfitted feature.
  - On the mirror branch (the runner-up displaced away from the unfitted feature), E*(w) falls with w at every step and the rest never folds, so a free coefficient grows.
  - At the second-moment order E* does not depend on w at all; its w-dependence enters at odd (third-moment) order, so the two branches are mirror images.
  - I mirrored C2's near tie across the face, keeping feature 2 at 6.08σ but on the entrant's side. A tail that C4 sees captured is then held on the far side and grows. At x = 7, ε = 10⁻¹³ the loop ends with two kernels at feature 0, one more than (iv)'s count.
- **Rounding picks the side for formed runner-ups.**
  - In the inner layer the runner-up and the entrant collapse onto one point, to 1.5·10⁻¹⁴σ in the loop.
  - For every formed runner-up (gap x ≲ 6 here), the integration error at that point decides the side. Toward the unfitted feature, the runner-up holds, departs and races; away from it, it is held and its coefficient grows.
  - In the ε-free problem, the wrong side or no exit occurs in 12, 11, 7, 2 and 0 of 16 runs at rtol 10⁻⁹, 10⁻¹⁰, 10⁻¹¹, 10⁻¹² and 10⁻¹³.
  - C4b's one disagreement is this effect, not "the saddle δ = 0". The rtol-10⁻¹² run left to the far side, 0.270σ from feature 0, and its coefficient grew from 0.02150 to 0.02320 by τ₂.

Other findings:
- C4's "stranded" tails never left their own feature (they stayed within 0.44σ of it).
- V3's "entrant to A′" cell is an artefact of how v_job starts its continuation; the fold there is c = 3.135.
- The fold that C3 and C4 divide by is not the configuration's own fold (0.01668, not 0.0172).
- Several numbers in the draft were computed by hand and have wrong ranges.
- 13 of 25 mutants survive, 6 of them changing the md.

## 1. What I checked and how

**Environment and repository state**
- I read the repository at 886793f and wrote nothing to it. In it I ran only `git log`, `git show`, `git status` and `cmp`.
- During the review, other sessions moved HEAD to 265a133 (commits for D‴.52, D‴.57 and D‴.58). None of them touches D54's files: `cmp` confirms the seven D54 files equal the copy I reviewed, and the round-49 log entry is unchanged.
- The `.pyc` files dated 10:32–10:33 in the repository's `experiments/__pycache__` and `tests/__pycache__` (gitignored) are not mine. Every one of my runs imported from the scratch copy, and nothing under `venv/` or `src/` changed.
- Scratch tree: /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref53/ (scripts in …/ref53/scripts, outputs in …/ref53/out).
- I used the repository's venv python with OMP_NUM_THREADS=1 and nice -n 15, at most two processes at a time. The longest run took 145 s. I killed no process.

**Regeneration and tests**
- `--from-raw` on the copy reproduces the md byte for byte (md5 0564f14b…).
- 22 of 22 tests pass in 20.6 s (PYTHONDONTWRITEBYTECODE=1, `-p no:cacheprovider`, run in the copy).

**Numerics paragraph.** I checked every number against the md, and the hand-quoted T3 and C4 numbers against the raw (`t3check.py`, `raw_c4.py`).

**(v), static checks**
- The projected law and its next order, by hand.
- The A′ → 0 limit at d = 5, 6, 7σ (`vlimit.py`).
- The branch away from the unfitted feature, in 1D (`away.py`).
- V3's (5σ, A′/A = 3.5) cell (`v3check.py`, `v3fold.py`).
- The static fold in C4's own planar geometry, on the samples and with feature 1 (`fold_c4.py`).

**(v), the ε-free reduced problem** (`reduced.py`, `side.py`, `delta.py`)
- The model: C2's near tie at μ → 0, with the entrant (0, 2H) and the runner-up (H, 2H). It starts as `inner_k` does (xi0 behind the wall, the runner-up x behind), runs Λ₀ times the free law at b0 = 0, keeps features 1 and 2 in the residual, and integrates in s to 3·10⁸. This is C4 with s₂ = ∞.
- A scan of 110 runs: x = 2.5 … 10, xi0 ∈ {400, 600}, rtol 10⁻⁹ … 10⁻¹³.
- The C3 gaps.
- A mirrored target: features 1 and 2 reflected across x = H/2, 30 runs.

**Loop reruns** (c4_job's settings; `c4traj.py`, `c4mirror.py`)
- (x, ε, rtol) = (4.9, 10⁻¹¹, 10⁻¹²) with the full trajectory.
- (7, 10⁻⁸, 10⁻¹²) followed past τ₂.
- The mirrored target at (7, 10⁻¹³) and (7, 10⁻¹⁰).

**Mutants** (`mutate.py`, `mutate2.py`)
- 25 single replacements in `exp_transit_nd.py`: v_job and its energy, fold_level, the V3 cells, c4_job's departure and hold, the C4b agreement count, c3_outcome, the D2 shift, T3's product, v_law, V's secants and C3's fold column.
- Each was run against the tests (`pytest -x`) and the md regeneration (TN_OUT pointed at a scratch directory), with the md diffed against the original.
- Round 49's three surviving `inner_k` mutants, re-run against the current tests.

**Literature.** Web searches; arxiv.org and ams.org were not reachable.

## 2. Verdicts

**(v), the static law: holds.**
- *Re-derived:*
  - (∂²_eφ)_⊥ = σ⁻²(u_e² − ½)φ;
  - B = ‖(∂²_eφ)_⊥‖² = π^{n/2}σ^{n−4}/2 (Var u² = ½ for u ~ N(0, ½));
  - ∫(u_e² − ½)e^{−|u−m|²} = π^{n/2}m_e² with m_e = d/(2σ), which gives P = A′π^{n/2}σ^{n−2}(d/σ)²e^{−d²/(4σ²)}/4;
  - M₂ = 2P/B = σ²ℓ(d), and |δ|² = σ²ℓw_e/(w(w_e + w)).
- *No transverse pull.* For a direction at angle θ to e, P(θ) = P cos²θ, so the rest lies on e.
- *A′ → 0 limit (w = 0.3):*

  | d | A′ = 10⁻² | 10⁻³ | 10⁻⁴ |
  |---|---|---|---|
  | 5σ | 1.0125 | 1.0038 | 1.0012 |
  | 6σ | 1.0039 | 1.0012 | 1.0004 |
  | 7σ | 1.0010 | 1.0003 | 1.0001 |

  These are rest/law. The unprojected law's ratio tends to 1.044, 1.029, 1.021 instead.
- *What the limit shows.* The error falls like A′^{1/2}, that is, linearly in |δ|. The first correction is therefore odd in δ (the third moment), not "the binding softening", which would be even (D3az-10).
- *Fold in C4's own geometry.* On the samples, with feature 1 (−0.6 at 6.25σ) present, the toward fold is w_f = 0.01668, with the rest 0.842σ from feature 0. Without feature 1 it is 0.01720, V's line value. Feature 1's opposite-sign binding lowers the fold by 3.0% (D3az-5).

**(v), the hold: right on the toward branch, wrong as stated (D3az-1).**
- *At second order E* has no w.* min over M₂ of [(1/8)M₂²B − ½M₂P] = −P²/(2B), independent of w. The w-dependence first enters through the third moment M₃ = w(1 − w²/w_e²)δ³, whose sign is the side's.
- *Measured in 1D, w from 0.6 down to 10⁻⁴* (secant dE*/dw):

  | d | toward branch | far branch |
  |---|---|---|
  | 5σ | +3.0·10⁻⁵ … +3.8·10⁻⁴; folds at 0.154 | −5.8·10⁻⁴ … −3.0·10⁻⁵ |
  | 6σ | +9.92·10⁻⁸ … +4.35·10⁻⁵; folds at 0.0208 | −9.91·10⁻⁸ … −3.77·10⁻⁵ |
  | 7σ | +7.5·10⁻¹¹ … +8.8·10⁻⁷; folds at 0.0017 | −4.4·10⁻⁷ … −7.5·10⁻¹¹ |

  - The far secant is negative at 52 of 52 steps for each d, and there are far rests at every w down to 10⁻⁴: −1.59σ, −1.45σ, −0.74σ at d = 5, 6, 7σ.
  - The two columns are mirror images; the even part is at most 7%, at the fold end.
- *In the plane* (C4's geometry) the far rest persists down to w = 0.003, at 0.59σ, with no fold.
- *The hold ends; no stall.* On the toward side the reduced problem departs at:

  | x | 2.5 | 3 | 3.5 | 4 | 4.26 | 4.5 | 4.9 |
  |---|---|---|---|---|---|---|---|
  | w at entry | 0.110 | 0.080 | 0.057 | 0.040 | 0.034 | 0.029 | 0.0215 |
  | s_d | 4.8·10⁷ | 2.52·10⁷ | 1.06·10⁷ | 4.09·10⁶ | 2.43·10⁶ | 1.22·10⁶ | 3.32·10⁵ |

  The local exponents of s_d against w are 2.0 to 4.6.
- *x = 3.* C4's "bound at every mass" is the race's finite horizon, not a stall. The loop's coefficient at s₂ = 9.7·10⁶ (ε = 10⁻¹²) is 0.0637, against the reduced problem's 0.064 there, and the hold ends at 2.52·10⁷. Capture is predicted for ε ≲ (9.7/2.5·10⁷)² ≈ 1.5·10⁻¹³ on the toward side.

**(v), the race on the toward side: confirmed; s_d is ε-free.**
- *Reduced problem against C4's departures:*

  | x | ε-free s_d | C4 (10⁻⁸ … 10⁻¹²) |
  |---|---|---|
  | 4.9 | 3.32·10⁵ | 3.3–3.5·10⁵ |
  | 7 | 9.86·10⁴ | 0.98–1.04·10⁵ |
  | 10 | 1.52·10⁵ | 1.52–1.67·10⁵ |
  | 17.3 | 5.84·10⁵ | 5.84–6.28·10⁵ |
  | 24.5 | 1.73·10⁶ | 1.73–1.85·10⁶ |

  C4's output grid resolves s to a factor 1.107 per output.
- *Further departures below the fold:* 6.4·10⁴, 7.8·10⁴, 1.17·10⁵, 3.59·10⁵, 5.62·10⁵ and 1.02·10⁶ at x = 5.48, 6.04, 8, 14.5, 17.0 and 20.6. s_d is smallest near x ≈ 5.5. The capture set is therefore an interval of x that widens on both sides as ε falls.
- *The rule "captured iff s_d(x) < s₂", with s₂ = 9.70ε^{−1/2}:*
  - It predicts all 15 of C3's outcomes. The capture tie at x = 17.0 has s₂ = 5.48·10⁵ against s_d = 5.62·10⁵, within 2.5%.
  - It predicts 28 of C4's 30. The two exceptions:
    - (7, 10⁻⁸): s_d exceeds s₂ by 1.4%, yet the loop had already arrived. At τ₂ the runner-up holds 0.6995 at feature 2, while that feature's path entrant has 8.3·10⁻⁵.
    - (4.9, 10⁻¹¹, rtol 10⁻¹²): the run that went to the far side.
- *The reach grows only slowly.* The largest captured coefficient grows like s₂^{1/4.6 … 1/2}, that is, like ε^{−0.11 … −0.25}.

**(v), the far side, dynamically (D3az-1).**
- *Reduced problem, mirrored target.* At x = 7, 8, 10 and 17.3 the runner-up goes to the far side in 6 of 6 settings for each x, rests 0.23–0.9σ from feature 0, and never departs up to s = 3·10⁸. Its coefficient grows:
  - x = 7: 0.00494 → 0.0105 at s = 10⁷ → 0.153 at s = 3·10⁸; it crosses θ|A| = 0.008 at s ≈ 5.6·10⁶;
  - x = 10: 0.00080 → 0.0105 at s = 2·10⁷.
- *Loop, mirrored target:*
  - (7, 10⁻¹³): w goes 0.00494 → 0.02235 at τ₂, on the far side (0.27σ at τ₂, at most 0.49σ). The final counts are [2, 1, 1] ("bound"); (iv) gives 1 (x = 7 > x_θ = 6.30).
  - (7, 10⁻¹⁰): w goes 0.00494 → 0.00542, below θ|A|. The counts are [1, 1, 1], and c3_outcome labels it "stranded" although the tail never left (0.016σ from feature 0 at the end).
  - On C4's own target, x = 7 is captured at every ε.
- *Which side, in exact arithmetic.* The pair's collapse keeps its direction: the e-share of c_r − c_e stays 0.49–0.50 while |c_r − c_e| falls from 6·10⁻⁴ to 1.1·10⁻⁹σ. The side is therefore sign((c_r⁰ − c_e⁰)·e): +0.48 in C4 (toward), −0.48 in the mirror (far).

**C4b and rounding (D3az-2).**
- *The rtol-10⁻¹² run at (4.9, 10⁻¹¹), traced:*
  - |c_r − c_e| goes 6.1·10⁻⁴ (s = 110) → 1.1·10⁻⁹ (s = 302) → 1.5·10⁻¹⁴σ (s = 502). There the e-component flips sign (−3.7·10⁻¹⁴ at s = 833).
  - It then grows at λ ≈ 8.2·10⁻⁴ per unit s, from 10⁻¹² (s = 6.3·10³) to 1.2·10⁻⁴ (s = 2.9·10⁴), and reaches the far rest at 0.270σ by s = 4.8·10⁴.
  - On the far rest w grows 0.02150 → 0.02200 (s = 9.1·10⁵) → 0.02320 at τ₂.
  - After feature 2's own entrant fits that feature, the runner-up returns to 0.0021σ from feature 0. That is the md's "0.0021σ at the end".
- *Why tightening the tolerance helps.* In the reduced scan at x = 2.5–6 (80 runs), 21 runs leave to the far side and 11 never leave the saddle by s = 3·10⁷. By rtol (10⁻⁹ … 10⁻¹³) the wrong-side-or-stuck counts are 12, 11, 7, 2, 0 of 16. At x = 7–10 all 30 runs go toward, with identical departures.
- *What this means for the reported runs.* C4's formed runner-ups at rtol 10⁻¹² are thus about 1-in-8 draws. C4b at 10⁻¹³ went toward in 10 of 10. C3 ran at 10⁻¹² and C2 at 10⁻⁹.
- *The departure time itself.* For formed runner-ups, the saddle exit starts from the integrator's error floor. The exact inner problem's offset is smaller (its collapse continues below 10⁻¹⁴σ), which adds λ⁻¹ ln(floor/exact) to s_d: 2.8·10³ per decade of offset at x = 4.9.

**C3, C4 and C2's "stranded" (D3az-3).**
- C4's six "stranded" rows (x = 10 at 10⁻⁸; 17.3 at 10⁻⁸, 10⁻⁹; 24.5 at 10⁻⁸ … 10⁻¹⁰) stayed within 0.17–0.44σ of feature 0 throughout (d_held), and ended 0.08–0.43σ from it. Each had s_d(x) > s₂, so they had not departed. They are not "away from both features".
- C2's ε = 10⁻¹⁰ tail (gap 30) records no distance. The ε-free problem puts it 0.39σ from feature 0 at s₂ = 9.7·10⁵ (it departs at 3.26·10⁶).

**V and V3.**
- All md values reproduce from the raw.
- The (5σ, 3.5) cell is an artefact (D3az-4). A branch of rests runs from w = 0.79 (w_e = 0.032 at −1.20σ) down to a fold at w_f = 0.42369, so c = 3.135, matching D3av's 3.14. v_job's cold start at 0.99A puts the free kernel at 1% of A, and LM sends it to A′.
- "267 steps" counts 51 twice: V3's A′/A = 0.875 column repeats V's d = 5, 6, 7 branches, so there are 216 distinct steps (D3az-7).

**(i)–(iv) and the other numerics.**
- All checked numbers match the md: A, T and T2, B, B4 and B4x, C, D, D2, H, C2, W1–W5, V, V3, C3, C4, and T3's counts (9 of 10 and 6 of 10).
- D2's shifts give 128.7, 130 and 130·ε^{1/2}.
- Draft numbers that are not in the md (D3az-7):
  - "0.0008 of the fold" (the md prints "(0.00)");
  - "(and 10⁻¹¹)": no rtol-10⁻¹¹ run exists in the raw or the md;
  - T3's percentage ranges;
  - "about 1.5%" for the B4 grid error.

**Round-49 fixes.**

| issue | status |
|---|---|
| D3av-1 | partly fixed. (v) is a race; C4 and C4b were added; "1.1 times" was dropped. But the far side is missing (D3az-1), C4b is misdiagnosed (D3az-2), and (iv)'s sentence is garbled (D3az-8). |
| D3av-2 | fixed: statement, proof, v_law and test; verified to 1.0001–1.0012 at A′ = 10⁻⁴ |
| D3av-3 | fixed (0.94–0.76σ); the stated reason, "the binding softens", is not supported (D3az-10) |
| D3av-4 | fixed, except the (5σ, 3.5) cell (D3az-4) |
| D3av-5 | reworded, but "away from both features" contradicts C4's own d_end (D3az-3) |
| D3av-6 | fixed; w_unbind remains as dead code (D3az-12) |
| D3av-7 | fixed: 213 continued, 211 settled, 2 stopped by the clock, counts unchanged, shares moved ≤ 1.6·10⁻⁶ |
| D3av-8 | fixed (D2) |
| D3av-9 | fixed; the hand-quoted ranges are partly wrong (D3az-7) |
| D3av-10, -11, -13, -14 | fixed |
| D3av-12 | partly fixed: round 49's inner_k survivors still survive (D3az-9) |

**Tests and mutants (D3az-9).**
- 12 of 25 mutants are killed.
- Survivors that change the md:

  | mutant | md lines changed |
  |---|---|
  | D2's x_θ taken at t = 0 | 3 |
  | D2's shift sign | 6 |
  | T3's product without the (1 − p) factors | 10 |
  | T3's planar product at t = 0 | 10 |
  | V's secants over V only | 1 |
  | C3's fold column with d for d² | 30 |

- Survivors that leave the md unchanged:
  - v_job with 4 fold bisections, and along_fold from the unbisected state;
  - V3's lower bracket (equivalent to 3 decimals);
  - c4_job's departure at 1σ, its held window including the departure, and its w at the departure output. No test runs c4_job, and `--from-raw` cannot see it.
  - C4b's agreement counted by n_fin (equivalent on the data).
- inner_k: the post-chunk clock check and the zero-drift branch still survive. `max_secs = 0` trips the in-integration check, never the post-chunk one. "≤ vs <" survives too, and is probably equivalent.

**Literature.**
- I found no precedent for:
  - the hold and the race of a partly formed particle against a second target's activation;
  - the gap-dependent k-fold shares;
  - the codimension law for ties at Voronoi faces.
- Duval–Peyré's thin-grid Lasso (two neighbouring spikes per spike), Denoyelle–Duval–Peyré–Soubies's sliding Frank–Wolfe and Aurenhammer's power diagrams are correctly cited.
- Missing neighbours (D3az-11):
  - Chizat's conic particle gradient descent (Math. Program. 2022): the Wasserstein–Fisher–Rao particle flow for sparse deconvolution, with positions and multiplicative weights and several particles per spike — the closest dynamical relative of the free law;
  - Chizat–Bach (2018), on few-particle failure;
  - Boursier–Flammarion (JMLR 2025): early alignment and spurious stationary points from small initialisation, also in the mean-field limit — an analogue of bound or stranded runner-ups;
  - Pesme–Flammarion (NeurIPS 2023) and Kunin et al. (alternating gradient flows): dormant units, an analogue of the frozen tails;
  - Poon–Peyré (multi-dimensional super-resolution): for pairs of spikes in d ≥ 2 the relative position matters.
- Sources:
  - https://arxiv.org/pdf/1811.06416
  - https://arxiv.org/pdf/1503.08577
  - https://arxiv.org/pdf/1907.10300
  - https://arxiv.org/abs/1805.09545
  - https://jmlr.org/papers/v26/24-1523.html
  - https://proceedings.neurips.cc/paper_files/paper/2023/hash/17a9ab4190289f0e1504bbb98d1d111a-Abstract.html
  - https://arxiv.org/html/2506.06489v4
  - https://ar5iv.arxiv.org/html/1709.03157
  - https://link.springer.com/doi/10.1007/pl00009187

## 3. Issues

| id | severity | issue | proposed fix (where) |
|---|---|---|---|
| D3az-1 | major | **The hold is one-sided; the race misses the far side.**<br>- D54.md:139–141 claims E* increases with w "along every branch computed (V)". V computes only the toward branch.<br>- On the far branch the secant is < 0 at 52 of 52 steps (d = 5, 6, 7σ), with mirror-image magnitudes (6σ: −9.91·10⁻⁸ … −3.77·10⁻⁵ against +9.92·10⁻⁸ … +4.35·10⁻⁵). Far rests exist down to w = 10⁻⁴ (1D) and 0.003 (plane), with no fold.<br>- At second order E* − E₀ = −P²/(2B) does not depend on w; the hold is odd (third) order.<br>- Mirrored near tie, reduced problem: x = 7–17.3 are held on the far side and grow (x = 7: 0.00494 → 0.0105 at s = 10⁷, crossing θ at s ≈ 5.6·10⁶).<br>- Mirrored near tie, loop at (7, 10⁻¹³): [2, 1, 1], one more than (iv)'s count. At 10⁻¹⁰: [1, 1, 1].<br>- The side is sign((c_r⁰ − c_e⁰)·e); the collapse keeps the direction to 10⁻⁹σ. | (v), statement (D54.md:84–107): the runner-up displaced toward the unfitted feature decays and races; one displaced away is held and grows, and can form, giving (iv)'s count + 1, with a reach that grows as ε → 0. Name the side rule.<br>Proof: E*'s w-dependence is third-order and odd.<br>(iv): the exception can be ±1.<br>V: run v_job with the tail started at −10⁻³σe.<br>C4: add the mirrored target.<br>Qualify D‴.51(ii)(c) likewise. |
| D3az-2 | major | **Rounding decides the side for formed runner-ups; C4b is misdiagnosed.**<br>- The pair collapses to 1.5·10⁻¹⁴σ, and the sign flips there.<br>- The run at (4.9, 10⁻¹¹, 10⁻¹²) went to the far rest (0.270σ), with w growing 0.02150 → 0.02320 by τ₂. It was not "on the saddle δ = 0 … to the end" (D54.md:219, and the log's D3av-1 row).<br>- Reduced problem, x = 2.5–6: wrong side or stuck in 12, 11, 7, 2, 0 of 16 runs at rtol 10⁻⁹ … 10⁻¹³.<br>- So "29 of 30 agree" is not evidence of robustness; C3 (10⁻¹²) and C2 (10⁻⁹) are exposed.<br>- s_d contains λ⁻¹ ln(1/offset) (2.8·10³ per decade at x = 4.9). | c4_job and c2_job: record \|c_r − c_e\|, its e-component after the layer, and the side at τ₂.<br>Run formed runner-ups at rtol ≤ 10⁻¹³, or carry the relative offset as a state variable, or seed it along the approach direction and report s_d's dependence on the seed.<br>Correct D54.md:218–220 and the log. |
| D3az-3 | minor | **"Stranded … away from both features" is wrong for C4.**<br>- C4's six stranded tails stayed within 0.44σ of feature 0 (d_held 0.17–0.44σ, d_end 0.08–0.43σ; md:302–314), with s_d > s₂: they never departed.<br>- C2's ε = 10⁻¹⁰ claim (D54.md:199) has no measurement; the ε-free problem gives 0.39σ at s₂.<br>- c3_outcome's residual cut labels a tail by w·offset, not by its position (the mirrored 10⁻¹⁰ run: "stranded" at 0.016σ). | (v) (D54.md:102–103), C3's caption, c3_outcome: classify by the recorded distance, as "not departed by s₂ (recalled)"; record the distance in c2_job. |
| D3az-4 | minor | **V3's (5σ, A′/A = 3.5) cell is an artefact.** The branch exists from w = 0.79 to w_f = 0.42369, c = 3.135. "Off the branch" (D54.md:207) is wrong. test_transit_nd.py:330–331 enshrines the artefact. | v_job (exp_transit_nd.py:784–787): continue from a bound state found by fresh starts; correct the cell, the draft and the test. |
| D3az-5 | minor | **The fold is unlike-for-like.** C3 and C4 divide by V's 1D single-feature fold, 0.0172. C4's own geometry gives 0.01668 (−3.0%, from feature 1). Against it the departures are at 0.97, 1.00 (C4) and 1.03 (C4b), i.e. at the fold, not "0.94–0.97 of the fold" (D54.md:215). w_dep is the last output before the 2σ crossing. | fold_level, C3/C4 (exp_transit_nd.py:850–858, 1429–1432): use the configuration's own fold (a continuation of about 1 min). |
| D3az-6 | minor | **The race's statement is underspecified.**<br>- "The fold is the capture threshold only in the limit of a zero hold" (D54.md:104–105) reads as a hold of zero duration, the opposite limit. No s₂ makes the fold an edge of the capture set.<br>- s_d(x) is given only at C4's six gaps, quantised to a factor 1.107. At (7, 10⁻⁸) the outcome is decided within 1.4%. | (v): state the capture set {x: s_d(x) < s₂} and tabulate s_d(x) from the ε-free problem (seconds per x, §2); say "no decay before s₂" for the limit. |
| D3az-7 | minor | **Hand-computed or unrecorded numbers.**<br>- "0.0008 of the fold": the md prints (0.00).<br>- "(and 10⁻¹¹)": no such run is recorded.<br>- T3: "exact to 0.2%" (it is −0.24% at the gap-9 corners); "0.3–7%" (−0.01 … −7.2%); "20–28% at ≈ 10⁻³" (20–37%; −50 … −76% at 2–6·10⁻⁴); "H/2's 0.011" (0.0110, 0.0106, 0.0102).<br>- "267 steps": 216 are distinct.<br>- W4's range omits x_r = 4 (−1.8·10⁻⁴).<br>- "About 1.5%" for the B4 grid is 1.2% (15.078 against 15.261), and it is the product rule's grid error, not the 4-kernel problem's. | report(): print T3's per-corner relative errors, C4b's w_dep and the rtol-10⁻¹¹ run; quote ranges from the md (D54.md:179, 205, 216, 220–226, 235). |
| D3az-8 | minor | **(iv)'s sentence has no main clause** (D54.md:65–67). The destination's count is not covered either: C3's capture tie (x = 17.0) gives feature 2 two kernels, where (iv) gives 1, and feature 2 has no same-sign unfitted feature at its own entry. | Rewrite as two sentences, and name the destination of a migrant. |
| D3az-9 | minor | **Tests miss the round-49 code.** 13 of 25 mutants survive; 6 of them change the md (D2 ×2, T3 ×2, V secants, C3 fold column). c4_job is untested. inner_k's post-chunk clock check and zero-drift branch are still untested (D3av-12 not met). | Add report tests on synthetic D2, T3 and C3 raws, a c4_job test on a cheap target, and inner_k with max_secs > 0 and a zero drift. |
| D3az-10 | minor | **"The binding softens as \|δ\| grows"** (D54.md:95, 139) is unsupported. The first correction to the law is linear in \|δ\| (rest/law − 1 falls like A′^{1/2}: 1.25%, 0.38%, 0.12% at d = 5σ), an odd effect. | Proof: attribute it to the third moment. |
| D3az-11 | minor | **The novelty claim needs its neighbours.** Missing: Chizat's conic particle GD (WFR), Chizat–Bach, Boursier–Flammarion, Pesme–Flammarion / Kunin et al., Poon–Peyré. | Status (D54.md:153–162). |
| D3az-12 | cosmetic | w_unbind is dead code, with a stale comment (exp_transit_nd.py:729–738, 1044–1047). "A′/A = 0.44, 0.88" should read 0.437, 0.875. The md prints "(0.00)" for the tails. "max w held" on the far side is the grown coefficient. | code; md formats |

## 4. Lesson

A sign computed on one branch does not settle a drift whose sign comes from a broken symmetry. At the second-moment order the held runner-up's energy does not depend on its coefficient at all. The hold lives at odd order, where the toward and far branches are mirror images. The round-49 revision computed the branch it expected and generalised its sign. Mirroring the geometry is the cheap test (2 s per gap in the ε-free problem), and it reverses the outcome: a tail C4 sees captured is held and forms. The same symmetry sets a numerical trap. The inner layer collapses the pair below double precision before an instability chooses the side, so for every formed runner-up the side depends on the integration error, and comparing two tolerances near that floor shuffles the error rather than testing convergence. Record the quantity that decides — here the pair's relative offset after the layer — and solve the ε-free problem, which reproduces the loop's departures in seconds, before scanning masses.
