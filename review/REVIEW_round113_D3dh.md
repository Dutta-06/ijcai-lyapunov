# Referee report D3dh, round 113: Proposition D‴.58 as revised after round 107

## Scope

- **Subject.** Commit 5be79c1 on claude/vigilant-clarke-hf939b, which applies round 107 (D3db), diffed against ce7cfac.
  - Exported with `git archive 5be79c1 | tar -x -C W/repo`, W = `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref113`.
  - WM_OUT = W/out and GP_OUT = W/gpout, each with a copy of its raw.
  - My scripts are in W/own, logs in W/logs, mutant work in W/mutk and W/mutf, and my own spec in W/mut113/spec113.py.
- **Read in full.**
  - D58.md, plus a word diff against ce7cfac.
  - Every hunk of `exp_width_mismatch.py` changed since ce7cfac:
    - `D61_STATUS` and a5's label;
    - the merged threshold line;
    - the distinct counts (`dist`, `dst_`, `dst`);
    - the outward exponent range;
    - the A6, A6E and D headings;
    - `a6x_lines` line by line: `fit4`, `ratio_hp`, the windows, `frame`/`fr`, `sig3`/`rng3`, the later window, the last window over the proved rule's bias, the four-term fits, `nlab`, `nfast` and the verdict.
  - width_mismatch.md.
  - `tests/test_width_mismatch.py`: both `test_round107_*` tests line by line, and the round-107 edits to the round-88, 92, 96 and 101 tests.
  - `tools/mutants/round107_D58_own.py`.
  - REVIEW_round107_D3db.md and the Round 107 entry of REVIEW_LOG.md.
  - From D‴.61:
    - D61.md (iii)(b), (iv) and its Status;
    - gh_peano.md parts D and V;
    - the raw's `DB` record. DB has no md heading of its own: `report` substitutes it for D's branch and bound.
  - HANDOFF's changed rows.
- **Run settings.**
  - Every run used the repository's venv with `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 10`. pytest always ran with `--basetemp` inside W.
  - My runner `W/own/mutants113.py` differs from `tools/mutants.py` in three lines: REPO = W/repo, PY = the repository's venv, and a per-slot `--basetemp`.
  - Two deviations, both mine:
    - Once, three short exact-error jobs ran at the same time (≈11 s each).
    - One foreground `sleep` was moved to the background by the harness.
  - Otherwise at most two heavy processes of mine ran at a time.
  - I stopped three of my own mutant-runner processes (two in W/mut107, one in W/mutf) and their children, after checking /proc/PID/cmdline and cwd, to restart them faster. All my processes have exited.
- **The repository was not modified.** I ran only `git archive/log/diff/show/status`.
- **The repository moved during the review.**
  - HEAD moved from 5be79c1 to eb22ec4 (two commits).
  - Neither commit touches D58.md, `exp_width_mismatch.py`, width_mismatch.md, its raw or its tests.
  - eb22ec4 edits D61.md (part VC; the Status's tier wording is reflowed but unchanged) and gh_peano_raw.npz. I checked the new raw: its V and DB records are identical to 5be79c1's.
- **Tags.**
  - [md]: the committed md or draft.
  - [raw]: the committed raw.
  - [own]: my code, or the repository's functions on inputs I changed.
  - [snippet]: web-search results only.
  - [memory]: recalled facts.

## Summary

**Design.**
- Regenerated both mds from their raws and ran the 63 tests.
- Reran the repository's `a6_job` on one A6 and one A6X rule against the raw. Computed the exact errors at every K = 19 … 61 for three rules.
- Recomputed every new A6X number three ways: the float solve, exact rationals from the doubles, and exact rationals from the stored 25-digit errors (t(v) exact).
- Computed the shrink factors that a pure power remainder gives for these windows.
- Ran four-, five- and six-term fits.
- Planted limit offsets at each A6X rule, and at the calibrator itself, through the repository's `a6x_lines`.
- Ran `a6_job` at K = 200, 240, 320 for the two other proved fast-branch rules at τ/σ = 2 (Mercer frame, 1.25·v_hs). With them I tested the verdict under another calibrator, and on a proved rule relabelled as unproved.
- Checked D‴.61's tier against its raw.
- Mutants: round 107's spec, plus 33 of my own and a control, each first against the A6/round tests with the md compared, and the survivors then against the other 54 tests.
- Five web searches.

**Verdicts.**
- **Reproducibility is exact.**
  - Both mds are byte-identical.
  - 63/63 tests pass.
  - The A6 (1, moment matching) rerun and the A6X (2, 2·v_hs) rerun are identical to the raw in every field.
  - My consecutive-K errors reproduce the stored e₀ and e₁ at K = 20, 40 and 60 to all 25 digits.
- **Every new number in D58 and the md is right**, with one inward bound: D58's Numerics "within … 5.85·10⁻⁹" (the true value is 5.8504·10⁻⁹).
- **The tier "proved by D‴.61 (iv) for t < 0.9900 … both steps rigorous in ball arithmetic" is earned.**
  - D61's Status says it.
  - D‴.61's raw records both rechecks as negative balls at 120 bits (python-flint 0.9.0), with no failed or open box, covering 1 < |η| ≤ 1.8927 and then 1.8927 ≤ |η| ≤ 1.995 contiguously.
  - σ*/σ < 199 is exact.
  - The raw holds aggregates, not the boxes, so the rigour itself rests on D‴.61's code and its referees (rounds 104, 108, 110).
- **A6X's "calibrated by it" is not a calibration.**
  - The threshold 3 was set post hoc, with 2.70 in hand, on one calibrator.
  - At the same τ/σ = 2 and in the same window, the three proved fast-branch rules' biases differ by 29× (my reruns).
  - Relabel one of them as unproved and the criterion rejects its proved limit ("up to 29 times: NOT calibrated").
- **The four-term fits show less than the md implies.** The outlier's −6.5·10⁻¹³ is a transient of the formula's second branch, with rate (h(τ²/v)/t(v))² = 0.850630 per kernel, measured at 0.8503–0.8508. It is not a fit bias. Modelled, the outlier's limit lies within 6.5·10⁻¹⁶ of 1/t(v).
  - Six-term fits put the three moment-matching limits within 2·10⁻¹⁸ … 2·10⁻¹⁷ of 1/t(v), at the proved rule's level of 7·10⁻¹⁸.
  - So the conclusion stands, much more strongly, by an argument the draft does not make.
- **Mutants.**
  - Round 107's spec reproduces the log: 19 killed (3 of them only by the md regression), and Y17 equivalent on the data.
  - Of my 33, 15 are killed by unit tests and 9 only by the md regression (7 of the 12 such kills over both specs are in the four-term path). The 9 survivors are all equivalent on the data; 3 of them would matter with more calibrators.
- **No major issue.**

**Findings.**
1. **D3dh-1 (minor).** "Calibrated by the proved rule" is not a calibration. The threshold is post hoc, on a single sample of a bias that varies 29× among proved rules. It rejects a proved limit. Its acceptance window is asymmetric and flips with a 1.0·10⁻¹² offset. The "later window" comparison is window geometry. At τ/σ = 1/4 no calibrator of this class can exist.
2. **D3dh-2 (minor).** What the four-term fits show. The outlier's 6.5·10⁻¹³ is its second-axis transient (predicted −6.47·10⁻¹³, found −6.52·10⁻¹³). Other proved rules' four-term fits are 1.1–1.2·10⁻¹², above the outlier's. The evidence that does hold is the geometric convergence of successive orders, plus the transient model.
3. **D3dh-3 (minor).** The verdict's logic:
   - it prints "calibrated" while its own necessary condition fails;
   - it never checks the calibrator's own convergence;
   - it ignores the later windows and the four-term fits that D58 calls the calibration;
   - the plants never touch the four-term path;
   - "(A6X extends 2 of the 1 …)" in the unproved-frame case.
4. **D3dh-4 (cosmetic).** Test coverage: 12 md-changing mutants are caught only by the md regression; survivors Z06, Z08 and Z09 are unconstrained; V's contiguity with DB is unpinned; one weak assertion.
5. **D3dh-5 (cosmetic).** Residues: "within 5.85·10⁻⁹" inward; HANDOFF §5 stale; "the 2 rules unproved"; round 107's R19/R21 not recorded.

## What I checked and how

**1. Reproduction** [md, raw, own]
- `--from-raw` into W/out and W/gpout, then `cmp`: both identical.
  - width_mismatch.md md5 d4a793cd…; its raw f1079628… (unchanged since ce7cfac).
  - gh_peano.md 0ea434f5…; its raw 515fe6a9….
- Full test file: 63 passed in 54.5 s.
- `W/own/job.py` runs the repository's `a6_job` and compares every field except `secs` by `==`:
  - A6 (1, moment matching): identical, 18.8 s.
  - A6X (2, 2·v_hs) at 200, 240, 320: identical, 97.9 s.
- `W/own/parity.py` calls `a6_error` at K = 19 … 61 in the job's precision for (¼, 2·v_hs), (2, 2·v_hs) and (1, moment matching). It reproduces e₀ and e₁ at K = 20, 40, 60 to all 25 digits (asserted).

**2. D‴.61's tier against its raw** (`W/own/insp.py`, `gpcheck.py`) [raw, md]

| record | range | ball bound (120 bits, flint 0.9.0) | failed / open |
|---|---|---|---|
| DB | 1 < \|η\| ≤ 1.78 | −3.2813·10⁻⁷ | 0 open |
| DB | 1 < \|η\| ≤ 1.8927 (402664 boxes) | −1.9432·10⁻¹¹ | 0 open |
| V | 1.8927–1.98 | −2.4767·10⁻¹²; end piece −1.38268 | 0 failed |
| V | 1.98–1.99 | −1.0258·10⁻¹²; end piece −1.35714 | 0 failed |
| V | 1.99–1.995 | −1.0004·10⁻¹²; end piece −1.34430 | 0 failed |

- D's own `dom.bb` carries no ball fields, so `report` takes DB, as round 110 arranged.
- The ranges are contiguous. r*(1.995) = 0.990046 rounds down to 0.9900. (1 + 0.99)/(1 − 0.99) = 199 exactly.
- c_need(0.99) ≤ 1.99497695 (gh_peano.md V) is computed at report time, not stored.
- D61's (iv) and its Status read "for r < 0.9900 … with V, both steps rigorous in ball arithmetic" at 5be79c1, and the same at eb22ec4.
- So D58's statements at (i) :14–16, (ii) :25–26, Status :171–174, :183–193, the sketch :273–277 and New :295–300 match their source.
- D58's "without them for t < 0.0966" matches D61's "(for r < 0.0966 without it)".

**3. A6X recomputed** (`W/own/a6x_check.py`; fits from the 25-digit errors, t(v) exact) [own, raw]

| rule | fit 1 at 80–160 / 160–240 / 200–320 | fit 2 (same windows) | shrink fit 1, overall / later | shrink fit 2, overall / later | last window ÷ proved (fit 1 / fit 2) |
|---|---|---|---|---|---|
| (¼, 2·v_hs) | +8.7052·10⁻⁸ / +4.9173·10⁻¹² / +1.5465·10⁻¹² | +3.5976·10⁻⁷ / −6.6644·10⁻⁹ / −3.3378·10⁻⁹ | 5.629·10⁴ / 3.18 | 107.8 / 1.997 | 0.462 / 0.749 |
| (½, MM) | +1.9291·10⁻¹⁰ / +2.2930·10⁻¹¹ / +9.0281·10⁻¹² | +1.9118·10⁻⁸ / +3.6612·10⁻⁹ / +1.8091·10⁻⁹ | 21.37 / 2.54 | 10.57 / 2.024 | 2.699 / 0.406 |
| (1, MM) | −7.5169·10⁻¹¹ / −8.9683·10⁻¹² / −3.5349·10⁻¹² | −5.8904·10⁻⁸ / −1.1718·10⁻⁸ / −5.8504·10⁻⁹ | 21.26 / 2.537 | 10.07 / 2.003 | 1.057 / 1.313 |
| (2, MM) | −1.6811·10⁻¹¹ / −2.2259·10⁻¹² / −9.0022·10⁻¹³ | −3.5212·10⁻⁸ / −7.0284·10⁻⁹ / −3.5122·10⁻⁹ | 18.67 / 2.473 | 10.03 / 2.001 | 0.269 / 0.788 |
| (2, 2·v_hs), proved | −6.9382·10⁻¹¹ / −8.4436·10⁻¹² / −3.3450·10⁻¹² | −4.4933·10⁻⁸ / −8.9285·10⁻⁹ / −4.4564·10⁻⁹ | 20.74 / 2.524 | 10.08 / 2.004 | 1 / 1 |

- **Every printed cell matches.**
  - Overall shrinks: 18.6–5.63·10⁴ and 10–108; the proved rule's 20.7 and 10.
  - Later window: 2.47–3.18 and 1.99–2.03; the proved rule's 2.52 and 2.
  - Last window: 0.268–2.7.
  - Four-term fits: within 6.6·10⁻¹³ (true −6.52·10⁻¹³); the proved rule's 2.0·10⁻¹⁴ (true −1.98·10⁻¹⁴).
  - Coverage: "4 of the 9"; the outlier's +1.55·10⁻¹² and −3.34·10⁻⁹.
- **Directed rounding.**
  - The 0.268 rests on the doubles' resolution: 0.2689 from the doubles, 0.2691 from the 25-digit errors. It is still a valid lower end.
  - D58's Numerics "within 9.03·10⁻¹² and 5.85·10⁻⁹": the second is inward (5.8504·10⁻⁹).
- **Window geometry** (`geom.py`, exact rationals). For a pure remainder K⁻ᵖ the extrapolated limit's error shrinks by a factor that depends only on p and the K values, not on the rule:

  | fit | remainder | 80–160 → 200–320 | 160–240 → 200–320 |
  |---|---|---|---|
  | fit 1 | K⁻⁴ | 21.11 | 2.533 |
  | fit 1 | K⁻⁵ | 46.5 | 3.177 |
  | fit 2 | K⁻³ | 10 exactly | 2 exactly |

  The proved rule's 20.74/2.524 and 10.08/2.004 are those numbers.

**4. Is it a calibration?** [own, raw]
- **The criterion.** The verdict is "calibrated" iff the maximum of |d_last|/|d_last(proved)| over the 4 unproved rules and both fits is ≤ 3, with one proved rule. On the data the maximum is 2.699 (MM ½, fit 1). The log records that the threshold was chosen in round 107 with the referee's 2.70 in hand.
- **Planted offsets** (`plants.py`, `plant2.py`: the repository's `a6x_lines`, every ratio of one rule in A6 and A6X raised by δ·(1/t(v))):
  - The outlier is "calibrated" for δ = −1·10⁻¹¹ … +5·10⁻¹², but not at −1.9·10⁻¹¹ (5.22) or +1·10⁻¹¹ (3.45).
  - MM ½ is "calibrated" for −1.9·10⁻¹¹ … +1.0·10⁻¹², and NOT at +1.01·10⁻¹². A 10⁻¹² perturbation thus flips the published verdict, and the printed text then reads "deviate by up to 3 times …: NOT calibrated".
  - MM 2 at +1·10⁻¹¹ prints: "NOT every deviation shrinks as the window moves out, and the extended rules without a frame deviate by at most 3 times the proved rule's bias in the same window: calibrated by it".
  - At the calibrator:
    - +5·10⁻¹², which cancels its bias, gives NOT (5.46);
    - +3·10⁻¹¹ gives "calibrated … at most 2 times", while the calibrator's own later-window factor is 0.808 (its deviation grows) and its overall factor 1.47.
- **The proved rules' biases differ.**
  - At A6's own window (80–160) the four proved fast-branch rules far from the cusp have fit-1 deviations:
    - −2.49·10⁻¹¹ at (1, 1.25·v_hs);
    - −6.94·10⁻¹¹ at (2, 2·v_hs);
    - −1.29·10⁻⁹ at (2, 1.25·v_hs);
    - −2.10·10⁻⁹ at (2, Mercer frame).

    Their fit-2 deviations are −5.6·10⁻⁸ … −1.74·10⁻⁷ (`a6cells.py`).
  - At 200–320 (my reruns, `altcal.py`), the three proved rules at τ/σ = 2:

    | proved rule | fit 1 | fit 2 |
    |---|---|---|
    | 2·v_hs (the calibrator) | −3.345·10⁻¹² | −4.456·10⁻⁹ |
    | 1.25·v_hs | −5.781·10⁻¹¹ | −1.372·10⁻⁸ |
    | Mercer frame | −9.694·10⁻¹¹ | −1.623·10⁻⁸ |

    The published calibrator has the smallest bias, 17–29× below the other two.
  - With either of the others as calibrator, the four unproved rules sit at 0.0093–0.427 of its bias: "calibrated" by a threshold 17–29× looser.
- **Relabelled proved rule** (`relabel.py`).
  - Setup: my (2, Mercer frame) A6X run added to A6X. Its v is moved so that ρ > 1, and its ratios and t(v) are left untouched. Its limit 7/3 is proved by D‴.61 and confirmed by A5.
  - Result: the repository's code prints "deviate by up to 29 times the proved rule's bias: NOT calibrated by it".
  - So the criterion rejects a true limit.
- **Same-τ/σ calibrators.**
  - At τ/σ = ¼ the frame interval is v/v_hs ∈ (0.9701, 1.0308), inside A6_NEAR. No far fast-branch frame rule can exist there.
  - At τ/σ = ½ one can, for v/v_hs ∈ [1.0513, 1.1180), but none is extended.

**5. What the four-term fits show: the outlier's second branch** (`transient.py`, `fit5_all.py`, `fit5_extra.py`) [own, raw]
- **Its residual is geometric.** I fitted the four-term power law through 160–320 and measured the outlier's ratio against it at K = 20 … 120:

  | K | 20 | 30 | 40 | 50 | 60 | 80 | 120 |
  |---|---|---|---|---|---|---|---|
  | residual | +3.13·10⁻³ | +6.40·10⁻⁴ | +1.28·10⁻⁴ | +2.54·10⁻⁵ | +5.03·10⁻⁶ | +1.98·10⁻⁷ | +3.01·10⁻¹⁰ |

  - It falls by 0.8506–0.8508 per kernel over K = 40–61, and by 0.8503 from K = 80 to 120.
  - The prediction is (h(τ²/v)/t(v))² = (−0.3197938/0.3467367)² = 0.8506295: the second argument of t(v)'s max, squared, as if e_K² summed the two axes' contributions.
  - The amplitude is B = 0.2383, the same in three fits.
- **The controls have no such term.**
  - At (2, 2·v_hs), (h₁/t)² = 0.0077; at moment matching, h₁ = 0.
  - Their residuals are 3–4 orders smaller (≈10⁻⁸ at K = 20, 6.5·10⁻¹⁴ at K = 120) and power-law: their per-kernel factor rises toward 1.
  - This is evidence for the two-axis structure beyond D‴.61's frames (ρ = 16.5 here).
- **What the transient accounts for.**
  - A6's outlier (8.7·10⁻⁸ at 80–160).
  - Its 5.6·10⁴ overall shrink, and its later-window 3.18 against the K⁻⁴ value of 2.53.
  - Its four-term deviation. The fit's weight on K = 160 is −1.3913, and the transient there is 4.65·10⁻¹³, so the predicted deviation is −6.47·10⁻¹³; the fit gives −6.52·10⁻¹³.
  - Its power fits worsening as terms are added: −6.5·10⁻¹³, then +9.3·10⁻¹¹, then −2.9·10⁻⁹.
- **With the transient modelled** (q fixed at the predicted value): +6.5·10⁻¹⁶ through 120–320, and +3.1·10⁻¹⁵ through 80–320 with a K⁻⁵ term.
- **Successive orders elsewhere** (four-term 160–320 / five-term 120–320 / six-term 80–320):

  | rule | four-term | five-term | six-term |
  |---|---|---|---|
  | MM ½ | −3.9·10⁻¹⁴ | −2.1·10⁻¹⁶ | +1.9·10⁻¹⁷ |
  | MM 1 | +8.6·10⁻¹⁵ | −3.3·10⁻¹⁶ | +2.1·10⁻¹⁸ |
  | MM 2 | −3.6·10⁻¹⁴ | +5.2·10⁻¹⁶ | −3.6·10⁻¹⁸ |
  | proved (2, 2·v_hs) | −2.0·10⁻¹⁴ | +4.9·10⁻¹⁶ | −7.4·10⁻¹⁸ |
  | proved (2, 1.25·v_hs) | +1.1·10⁻¹² | −3.5·10⁻¹⁴ | +1.9·10⁻¹⁵ |
  | proved (2, Mercer frame) | +1.2·10⁻¹² | +2.5·10⁻¹³ | −1.3·10⁻¹³ |

  - So the contrast "6.6·10⁻¹³ (the proved 2.0·10⁻¹⁴)" depends on which proved rule is used.

**6. Tests, line by line** [md, own]
- `test_round107_d61_status_pinned_…`:
  - It pins `D61_RSTAR_V` to `floor(GP.r_star(1.995))` and checks V's and DB's ball fields, `R == DOM_BB` and max R2 = 1.995.
  - It does not check that V's first R1 equals DB's last R, or that the bands are contiguous. A gap would pass.
  - It checks "both" and "ball arithmetic" as strings.
  - The verdict and coverage are asserted on `report()` ✓, and the 5·10⁻⁸ plant ✓.
- `test_round107_a6x_verdict_plants`:
  - The 5·10⁻¹¹ plant ✓, though its "deviate by up to 1" also matches 15.4 or 100.
  - The K⁻³ plant kills fit-1-only ✓; the t(v) ≥ r* mock ✓.
  - The mocks' e₀ = 1e-10 and e₁ = 1e-11 at every K make the four-term fit meaningless in every mock: the four-term path is pinned only by the md regression.
  - The same mock with a frame rule at t(v) ≥ r* prints "(A6X extends 2 of the 1 fast-branch rules without a frame)", which the test does not assert against.
- The edits to the round-88, 92, 96 and 101 tests (distinct counts, the outward exponent, the new threshold wording, `D61_RSTAR_V` in a5) are correct.

**7. Mutants** (`tools/mutants/round107_D58_own.py`; mine in W/mut113/spec113.py, Z01–Z33 plus control Z00) [own]
- **Method.** Pass 1 used `-k "(a5 or a6 or round or report or quoted) and not round88_precision_and_the_md"` (17 tests), with the md compared. Pass 2 ran the survivors against the other 54 tests. Both controls pass.
- **Round 107's spec.**
  - Killed by unit tests: 16 (Y01, Y02, Y04–Y08, Y11–Y16, Y18–Y20).
  - Killed only by the md regression: Y03 (frame census by name), Y09 and Y10 (the four-term form and source).
  - Y17 survives, equivalent on the data.
  - This matches the log's "19 killed … Y17 equivalent".
- **My 33.**
  - Killed by unit tests: 15 (Z05, Z10, Z13, Z16, Z17, Z19–Z23, Z25–Z27, Z30, Z31).
  - Changing the md but killed only by its regression: 9.
    - Z01: `ratio_hp` inverted.
    - Z02, Z03: the four-term values or K from 120.
    - Z14: the proved rule's single values rounded up.
    - Z15, Z33: the later window replaced by the far one, on either fit.
    - Z24: the frame census' distinct count to 2 decimals.
    - Z28: the proved rule's four-term maximum taken over the unproved.
    - Z32: the four-term rate perturbed by 10⁻¹².
  - Survivors, all 63 tests passing, all equivalent on the data: 9.
    - Z04 (`fit4` without its pivot division): equivalent by construction, since A[0][0] stays 1.
    - Z06: "proved" at D61_RSTAR = 0.8052.
    - Z07: ≤ r*.
    - Z08: the necessary condition on fit 1 only.
    - Z09: the smallest proved bias as the denominator.
    - Z11: cal < 3.
    - Z12: `round(cal)`.
    - Z18: a5's ≤.
    - Z29: the four-term guard `all`.
- **Which change published results.** None of the md-only mutants flips the verdict, because the verdict reads none of those numbers (D3dh-3). But they change the published md (later-window factors, four-term values, the census) while every unit test passes. Z06, Z08 and Z09 would change the verdict with data the tests do not plant: a frame rule with 0.8052 ≤ t(v) < 0.99, a fit-2-only growth, or two calibrators (as §4 shows).

**8. Literature** (five searches) [snippet]
- Karvonen–Oates–Girolami (Math. Comp. 2021) and Karvonen–Särkkä (BIT 2019): consistent with D58's descriptions; their rates are not sharp in the sense D58 claims.
- Arb (Johansson) and the FLINT docs; python-flint 0.9.0 with Arb balls in a recent computer-assisted bound (arXiv:2609.20074). This is the same toolchain as D‴.61's raw, and supports "computer-assisted, rigorous in ball arithmetic" as a recognised tier.
- Mikkelsen–López-Villellas 2025, on when Richardson's error estimate can be trusted.
- Sidi's GREP, whose general basis functions φᵢ allow adding a known-rate exponential term to the extrapolation basis. This is the precedent for §5's transient model.
- D58's literature passages are unchanged since ce7cfac.

## Verdicts per part

| claim (D58 lines) | verdict | evidence |
|---|---|---|
| (i) :14–16: "proved by D‴.61 (iv) for t < 0.9900 (σ*/σ < 199; computer-assisted, both … rigorous in ball arithmetic, and without them for t < 0.0966)" | Earned | [raw] §2; [md] D61 (iv), Status |
| (ii) :25–26: proved for frame rules with t(v) < 0.9900, computer-assisted | Earned | §2 |
| Status (i) :171–174 (199 exact; one tier) and (ii) :183–193 (19.924, 140.712; 4.187, 6.514; "(17 distinct)", "(21 of 35 distinct)") | Right | [md], round 107's closed forms |
| A5 :232–233: "0.4286 < 0.9900, computer-assisted" | Right | [md] |
| A6 :240–249: "32 rules (31 distinct)", exponent "1.99 … 2.04", "18 … (17 distinct)" | Right and outward | [md], [own] |
| :255–257: "an outlier, pre-asymptotic" | Right, and now explained: a second-axis transient, rate 0.850630 | §5; D3dh-2 |
| :257–260: A6X's numbers; "necessary, not a calibration" | Right | §3 |
| :260–263: "The same window calibrates: … as the proved rule's do …; 0.268–2.7 times …; four-term fits … within 6.6·10⁻¹³, the proved rule's within 2.0·10⁻¹⁴" | **Numbers right, the reading not** | §3–§5; D3dh-1, -2 |
| Sketch :273–277; New :295–300 | Right | [md] |
| Numerics (A6X) :328–329 | Right, except "5.85·10⁻⁹" (5.8504·10⁻⁹) | D3dh-5 |
| md: A6/A6E/D headings, a5, threshold line, distinct counts | Right | [md] |
| md :183 (A6X line), verdict "calibrated by it" | Numbers right; the verdict is not a calibration and its logic is defective | D3dh-1, -3 |
| (iii)–(v), their proofs, B, C, D | Unchanged since ce7cfac; the tests pass | [md] |

**What is proved, and how, at 5be79c1.**
- **Proved:**
  - (i)'s identity, its bounds and the rate t as a K-th root;
  - (v);
  - the order of λ, λ/(2 − λ) and q against t.
- **Proved, computer-assisted** (D‴.61, both steps in ball arithmetic):
  - the sharp constant for t < 0.99, i.e. σ*/σ < 199 (without computer for t < 0.0966);
  - the rate t(v) and the sharp asymptotics at every frame rule with t(v) < 0.99: A6's 22 (21 distinct), and the Fock and Mercer-zero rules below their thresholds.
- **Numerical:**
  - moment matching's λ. By my six-term fits, A6X's three limits agree with λ to ≲ 2·10⁻¹⁷;
  - the 14 frameless rules: the slow branch by bracket; 4 of the 9 fast-branch rules by A6X (the outlier within 6.5·10⁻¹⁶ once its transient is modelled);
  - (iii) and (iv).
- **Conjectural:**
  - t(v) without a frame (beyond these numbers) or with t(v) ≥ 0.99;
  - the sharp constant beyond σ*/σ = 199.

## Status of the round-107 issues

| D3db | status | note |
|---|---|---|
| 1 | fixed | One tier everywhere; both false statements removed; 199. `D61_STATUS` is pinned to D‴.61's raw (V's bands and DB's radii negative in balls, R2 = 1.995, `D61_RSTAR_V` = floor r*(1.995)). It is a constant, not derived from the raw: acceptable. Residue: V–DB contiguity is not pinned (D3dh-4) |
| 2 | partly fixed | "4 of the 9" ✓; "necessary, not a calibration" ✓; the later window, the last window and the four-term fits are printed ✓. But the new verdict, "calibrated by it", is itself not a calibration (D3dh-1), the four-term fits are misread at the outlier (D3dh-2), and the verdict's logic is defective (D3dh-3). The optional extension was not done; I ran two of the calibrators it named |
| 3 | fixed | Exponent 1.99 … 2.04; A6X ranges directed at 3 significant digits; 199. New residue: D58's "5.85·10⁻⁹" (D3dh-5) |
| 4 | fixed, residues | Distinct counts; singular forms; "proved" needs t(v) < r*; A6X in the Numerics; HANDOFF row. Residues: HANDOFF §5 still lists round 107 as running and asks for the fix it already made; "at the 2 rules unproved"; "(… 2 of the 1 …)" (D3dh-3) |
| 5 | fixed, residues | Assertions on `report()`; `D61_RSTAR_V` pinned; planted verdicts. Residues: the four-term path and the later windows are pinned only by the md regression; the referee's R19/R21 equivalences are not recorded in the log's fix column (D3dh-4, -5) |

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3dh-1 | minor | **"Calibrated by the proved rule" is not a calibration.**<br>• The verdict is "max over 4 rules × 2 fits of the last-window deviation ÷ one proved rule's ≤ 3". The threshold was set with the data's 2.70 in hand (margin 1.11).<br>• At MM ½ an offset of +1.01·10⁻¹² flips it. The acceptance windows are asymmetric: −1.9·10⁻¹¹ … +1.0·10⁻¹² at MM ½, about −1·10⁻¹¹ … +5·10⁻¹² at the outlier [own, repository's `a6x_lines`].<br>• The proved rules' biases at the same τ/σ = 2 and window 200–320 are −3.3·10⁻¹², −5.8·10⁻¹¹ and −9.7·10⁻¹¹ (fit 1). The calibrator is the smallest; at A6's 80–160 window the four proved fast-branch rules span 84× [own reruns, raw].<br>• (2, Mercer frame), whose limit is proved, relabelled as unproved: the criterion prints "up to 29 times: NOT calibrated".<br>• "As the proved rule's do (2.52 and 2)" is window geometry: K⁻⁴ and K⁻³ remainders give 2.533 and exactly 2 for any rule [own, exact].<br>• At τ/σ = ¼ no far fast-branch frame rule exists. At ½ one could (v/v_hs ∈ [1.051, 1.118)) but none is extended. | **D58 :260–263 and the md's verdict:** replace "calibrated by it" with what the data show:<br>• successive orders converge geometrically at the three MM rules as at (2, 2·v_hs) (4/5/6-term: ≲ 4·10⁻¹⁴, 5·10⁻¹⁶, 2·10⁻¹⁷);<br>• the outlier within 6.5·10⁻¹⁶ with its transient (D3dh-2);<br>• "under an expansion in integer powers of 1/K, numerical, not a proof".<br>If a threshold is kept, take it from several proved rules (add (2, Mercer frame), (2, 1.25·v_hs), (1, 1.25·v_hs) and a (½, 1.1·v_hs) rule to `A6X_RULES`; ≈95 s each at τ/σ = 2) and call it a heuristic. Say in the text that no calibrator can exist at τ/σ = ¼ |
| D3dh-2 | minor | **The four-term fits at the outlier measure a transient, not a bias.**<br>• The ratio's residual against the 160–320 power law falls by 0.8503–0.8508 per kernel (K = 40–120), = (h(τ²/v)/t(v))² = 0.850630, amplitude 0.2383. It is absent at the calibrator ((h₁/t)² = 0.0077) and at moment matching (h₁ = 0).<br>• It explains A6's 8.7·10⁻⁸, the 5.6·10⁴ shrink, the later window's 3.18 and the four-term −6.52·10⁻¹³ (predicted −6.47·10⁻¹³ from the weight −1.391 at K = 160). The outlier's five- and six-term fits worsen (+9.3·10⁻¹¹, −2.9·10⁻⁹).<br>• With the transient: +6.5·10⁻¹⁶ (120–320) [own, raw; my consecutive-K errors].<br>• The proved (2, 1.25·v_hs) and (2, Mercer frame) have four-term deviations of 1.1–1.2·10⁻¹², above the outlier's. | **D58 :255–257:** "pre-asymptotic: a geometric transient at the formula's second branch, (h(τ²/v)/t(v))² per kernel" — also evidence for the two-axis picture without a frame.<br>**`a6x_lines`:** where (h₁/t)² > 0.5, fit the transient at that fixed q; print the five- and six-term fits; drop the "(the proved 2.0e-14)" contrast or print all proved rules' values |
| D3dh-3 | minor | **The verdict's logic.**<br>• It prints "calibrated" while its necessary condition fails (MM 2 + 10⁻¹¹: "NOT every deviation shrinks …, and … calibrated by it").<br>• It never checks the calibrator itself: proved + 3·10⁻¹¹ (later factor 0.808) is still "calibrated"; proved + 5·10⁻¹², which cancels its bias, makes everything NOT.<br>• It ignores the later windows and the four-term fits, which D58 presents as the calibration.<br>• The four-term fits come from e₀/e₁ while every plant (the repository's and mine) changes `ratio`, so that path is never planted. Under a planted 10⁻¹¹ offset it still reads 6.6·10⁻¹³.<br>• An extended frame rule with t(v) ≥ r* gives "(A6X extends 2 of the 1 …)" (`len(nf)` against a frameless count).<br>• At cal ∈ (3, 3.005) it prints "up to 3 times … NOT calibrated".<br>No published number changes. | `a6x_lines`:<br>• make the verdict conditional on the shrink conditions, including the calibrator's own later-window factor near its geometric value;<br>• plant through e₀/e₁ as well as the ratio;<br>• count `nfast` over the same population as `nf` (fast branch, not proved);<br>• print cal to one more digit.<br>Tests for each case |
| D3dh-4 | cosmetic | **Test coverage.**<br>• 12 of 34 md-changing mutants are caught only by the md regression: the four-term path (Y09, Y10, Z01–Z03, Z28, Z32), the later window (Z15, Z33), the proved rule's single values (Z14), and the frame census' distinct count (Y03, Z24).<br>• Survivors Z06 (proved at 0.8052), Z08 (fit 2 dropped from the necessary condition) and Z09 (min over calibrators) are equivalent on the data but unconstrained.<br>• The D‴.61 pin does not check that V's R1 equals DB's last R or that the bands are contiguous.<br>• "deviate by up to 1" also matches 15.4. | Mock rules with consistent e₀/e₁ (from an exact model sequence) asserting the four-term value and the later-window factors. A mock with two calibrators (kills Z09), a frame rule with 0.8052 ≤ t(v) < 0.99 (Z06), and a fit-2-only growth (Z08). Assert `bs[0]['R1'] == DOM_BB[-1]` and `R2[i] == R1[i+1]`. Assert the full cal |
| D3dh-5 | cosmetic | **Residues.**<br>• D58 :329 "within 9.03·10⁻¹² and 5.85·10⁻⁹": the second is inward (5.8504·10⁻⁹).<br>• HANDOFF §5 still lists round 107 as in flight and asks for the D58 wording change that 5be79c1 made.<br>• md "at the 2 rules unproved" (in the unproved-frame case).<br>• Round 107's R19/R21 equivalences are not recorded in the log's fix column. | "5.86·10⁻⁹" (or `up_e`); HANDOFF §5; "at the 2 unproved rules"; log line |

## Sources

- [Karvonen–Oates–Girolami, Integration in RKHS of Gaussian kernels, arXiv:2004.12654; Math. Comp. 90 (2021)](https://arxiv.org/pdf/2004.12654) and the [AMS page](https://www.ams.org/journals/mcom/2021-90-331/S0025-5718-2021-03659-5/) [snippet]. Scaled Gauss–Hermite rules without a length-scale restriction; non-matching upper and lower bounds (the exponent form differs between versions).
- [Karvonen–Särkkä, Gaussian kernel quadrature at scaled Gauss–Hermite nodes, arXiv:1803.09532](https://arxiv.org/pdf/1803.09532) [snippet]. Nodes from a truncated Mercer expansion.
- [Johansson, Arb, arXiv:1611.02831](https://ar5iv.arxiv.org/html/1611.02831) and the [FLINT feature overview](https://flintlib.org/doc/overview.html) [snippet]. Midpoint-radius ball arithmetic; enclosures can overestimate.
- [A computer-assisted upper bound of 1.7813 for the real Grothendieck constant, arXiv:2609.20074](https://arxiv.org/pdf/2609.20074) [snippet]. Outward-rounded Arb balls via python-flint 0.9.0 for every accepted inequality: the same toolchain as D‴.61's raw.
- [Mikkelsen–López-Villellas, How Accurate is Richardson's Error Estimate? (2025)](https://zaguan.unizar.es/record/168390) [snippet]. When the estimate can be trusted, in terms of Richardson's fraction.
- [Sidi, GREP^(m), Math. Comp. (P107)](https://asidi.cswp.cs.technion.ac.il/wp-content/uploads/sites/142/2023/03/P107_Math.Comp_.GREP^m.pdf) and [GREP, JIMA 1979 (P006)](https://asidi.cswp.cs.technion.ac.il/wp-content/uploads/sites/142/2023/03/P006_JIMA.GREP_.pdf) [snippet]. Richardson generalised to arbitrary asymptotic basis functions.
- Repository [md, raw]:
  - D58.md, D61.md, width_mismatch.md and gh_peano.md at 5be79c1, and D61.md and gh_peano_raw.npz at eb22ec4;
  - REVIEW_round107_D3db.md; REVIEW_LOG.md Round 107;
  - `git diff ce7cfac 5be79c1`.
- [memory]: Richardson extrapolation's error is set by the first neglected term. A known case calibrates an unknown one only if the two share the expansion's structure.

## Lesson

**A calibration by a known case needs the known case to share the unknown one's error structure, and needs more than one known case.** Here the one proved rule has the smallest bias of the three at its τ/σ, by up to 29×, and none of the outlier's second-branch transient. The cheap test is to relabel a proved rule as unproved and see whether the criterion accepts it. It does not: "up to 29 times: NOT calibrated".

**When adding terms to a fit makes its limit worse, look for a term outside the basis before calling the data pre-asymptotic.** Here the formula's own second argument predicted the missing term's rate to four digits. With it included, the outlier's limit improves from 6.5·10⁻¹³ to 6.5·10⁻¹⁶.
