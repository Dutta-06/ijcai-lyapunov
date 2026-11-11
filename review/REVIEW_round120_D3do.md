# Referee report D3do, round 120: Proposition D‴.61 as revised after round 117 — certifying the last annulus

## Scope

**Commit and setup.**
- Reviewed: commit a2f3b77 on claude/vigilant-clarke-hf939b, exported with `git archive a2f3b77 | tar -x -C W/repo`. W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref120.
- Hashes match `git show a2f3b77:…` (sha256 prefixes) [own]: exp_gh_peano.py 435ce048…, test_gh_peano.py 773cc869…, D61.md 30dc054b…, gh_peano_raw.npz ad579462…, gh_peano.md 576fefbe….
- W/repo/venv is a symlink to the repository's venv (python-flint 0.9.0). GP_OUT pointed at W/out, which held copies of the raw, the md and width_mismatch_raw.npz.
- Everything ran under `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 10`, with at most two heavy processes of mine at a time. Pytest always ran with `--basetemp` inside W.
- I used only read-only git commands and wrote nothing into the repository. I killed only my own processes, each verified through /proc/PID/cmdline.
- During the round, other sessions moved the repository's HEAD to b2088e8 and left uncommitted noise_resonance changes and an untracked REVIEW_round119_D3dn.md. None of that is mine; this review is of the a2f3b77 export.

**Read in full.**
- review/drafts/D61.md, especially (iii)(b)'s Status with Lemma VC and the path family Im z = min(κRe z, θ), and Numerics (V) and (VC).
- In experiments/exp_gh_peano.py: `cert_E`, `cert_grad`, `cert_flow`, `cert_seg_*`, `cert_box*`, `cert_band`, `cert_end*`, `_cert_E_ball`, `corner_*`, `part_VC`, `vc_lines`.
- results/chaos/gh_peano.md (V, VC).
- review/REVIEW_round117_D3dl.md, all of it including §§6–7; review/REVIEW_round115_D3dj.md §3.
- tests/test_gh_peano.py (the V and VC tests).
- tools/mutants.py and tools/mutants/round117_D61_own.py.
- HANDOFF §§5, 7, 10 and the Round 117 entry of REVIEW_LOG.md.
- `git diff eb22ec4 a2f3b77` of these files.

**Deliverable.** My certificate is in W/own: lastann.py, selfcheck.py, coverage.py, vs_repo.py and a README. Logs are in W/logs and mutant runs in W/mut. The certificate is self-contained (python-flint only) and deterministic.

## Summary

**Design (the revision).** Round 117's fixes are wording and test changes to part VC:
- the junction from b;
- "the four paths are one family";
- a VC_MID point with Re σ > 1;
- the path grid's range stored from the grid used;
- Ψ_φ computed in balls, and the ε lists stored.

No proved claim changes. (iii)(b) is still a conjecture on 1.995 < |η| < 2.

**Design (my certificate).** I certify the valley condition along the family Im z = min(κRe z, θ). The path starts at b = (1 + ε/4)η/2 with a vertical junction, runs out to Re z = δ₀ = 0.05, and from there a connector plus a descent polyline per η-box goes to |1 − z| ≤ 0.02, followed by part V's end piece.

The key choice is to blow up in the **path point's** weight-0 coordinates rather than in the parameters:
- **Form A**: (u, q, ρ) = (ε/X, θ²/X, X).
- **Form B**: (u, p, ϑ) = (ε/X, X/θ², θ²).

Each form is an exact identity for E/X², with every division by ε, θ or X done analytically. The four regimes of part VC become faces of one compact (u, q)-square, and the matching region becomes its corner. The box bound is the mean value form f(ρ) ≤ f(0) + ρ·sup f′:
- f(0) is the weight-4 limit (Ψ̃ in (u, q) coordinates), in closed form and computed in endpoint intervals;
- sup f′ is obtained by forward differentiation in Arb at 128 bits;
- where f(0) vanishes (the junction's λ = 0, the face q = 0, and the corner u = q = 0), the certificate rests on sup f′ < 0.

**Verdicts.**
- **The whole last annulus is certified** (D3do-1). Every η with 1.995 ≤ |η| < 2 and 0 ≤ arg η ≤ π/2 has an explicit path with E < 0, proved in ball arithmetic:
  - 3140 parameter boxes, 75 η-boxes and 5126 outer segments, none failed, in about 3 s;
  - with ε₀ = 0.11 the same code covers 1.89 ≤ |η| < 2, overlapping the domination disc |η| ≤ 1.8927.

  With the reflections and D's branch and bound, the valley condition then holds on the whole open disc |η| < 2. Consequently (iii)(b) holds there (computer-assisted), (iv) holds for every frame with r < 1, and D‴.58's sharp constant holds for every σ*/σ.
- **The certificate checks itself.**
  - Its forms agree with `cert_E` and with F_η's definition to ≤ 7.6e-23 at 1479 points, and with the repository's `_cert_E_ball` to 2.3e-23 at 339 points.
  - The true bound lies below the certified bound at 9401 points sampled inside every accepted box.
  - All 28000 sampled path points lie in accepted boxes.
  - Wrong paths are rejected.
- **Reproducibility is exact.**
  - `--from-raw` regenerates the md byte for byte (sha256 576fefbe…4b2c).
  - All 39 tests pass (18.5 s).
  - The 7 shipped mutants are killed and the control survives.
  - Of my 10 mutants, 7 are killed; one non-equivalent survivor is cosmetic.
- **Round 117's fixes are in.** What remains is residue, all cosmetic (D3do-4 to D3do-6).

**Findings.**
1. **D3do-1 (major, constructive).** 1.995 ≤ |η| < 2 is certified in ball arithmetic (§ The certificate, § Coverage). D61's Status ("a conjecture on 1.995 < |η| < 2"), (iv)'s "conditionally on conjecture P" and the claims built on them are out of date once this is adopted.
2. **D3do-2 (minor).** Two charts suffice. "Four parameter charts, plus inner, matching and outer charts in z" overstates the remaining work. In (u, q) = (ε/X, θ²/X):
   - the ψ chart is the rays q/u = ψ²;
   - the φ and (ε/θ, θ/√ε) charts are the face q = 0;
   - matching is the corner u = q = 0, where f(0) = −uq − q²/2 and sup f′ ≈ −1/3, that is E ≈ −εθ² − θ⁴/2 − X³/3;
   - the (√ε/θ, θ) chart is q → ∞, handled by form B with p = 1/q and ϑ up to (π/2)².
3. **D3do-3 (minor).** Part V is subsumed: at ε₀ = 0.11, 4189 parameter boxes and 220 η-boxes certify 1.89 ≤ |η| < 2 in about 4 s, against V's 13429 η-boxes [md].
4. **D3do-4 (cosmetic).** D61 keeps the stale sentence "A proof needs, per regime, its path … the regimes' overlaps …" right after the corrected plan, which says no path matching is needed.
5. **D3do-5 (cosmetic).** "In Ψ̃'s coordinates exactly y = min(sx, ψ)" (D61, md) is exact only as ε → 0. At ε > 0 the floor is y = (A_i/|A_r|)x, with A_i/|A_r| = 2ψω²sinc 2θ/(1 − ε/4 + 2ω²ψ²sinc²θ). The cap y = ψ is exact.
6. **D3do-6 (cosmetic).** Residue of round 117:
   - the coef's ε labels: the md text uses a literal, and the stored `eps` is a second literal, not the loop's;
   - mutant Y02 (mid's labels taken from the constant) survives;
   - no VC_MID point lies on the cap's far piece;
   - the label "(0.2, 3, 0.1, Re σ > 1)" sits inside the tuple;
   - D61's Numerics (VC) still lacks rounds 115 and 117's numbers;
   - HANDOFF §7 still says "Write it as D‴.61 (iii)(b) on the whole disc, add the connectivity check to part D".

## What I checked and how

**(1) Reproducibility** [own]
- `--from-raw` into W/out: `cmp` identical to the committed md, sha256 576fefbe8d4f7f981a1117ab0bbf475ad47b990e7848bbba3ca2093b40a64b2c.
- `pytest -q -p no:cacheprovider --basetemp=W/pt/base1 tests/test_gh_peano.py`: 39 passed in 18.52 s (W/logs/pytest_gh.log). HANDOFF §10 also says 39.
- The tests' `assert GP.part_VC() == raw['VC']` holds, so part VC's recomputation equals the raw.

**(2) Round 117's fixes against the code and the md** [md][raw][own]
- **Junction and one-family wording.** Present in D61's Status, `vc_lines` and HANDOFF §5. Residue in D3do-4 and D3do-5.
- **The Re σ > 1 point (ψ, x, y) = (0.2, 3, 0.1).** I recomputed it from F_η's definition (not from `cert_E`) at 300 bits: −0.007273, −0.005839, −0.005825 at ε = 10⁻⁴, 10⁻⁶, 10⁻⁸, with Re σ > 1 at all three, which is the md's line. The md's other three points also reproduce verbatim (−0.890684…, 0.467173…, 31.513697…).
- **Stored grids.** `path` holds xmin = 10⁻⁴, xmax = 1602, xmax0 = 200, n_psi = 145, taken from the grid used.
- **Ψ_φ in balls.** Done; y is still the float `corner_floor`, which is harmless.
- **Coefficient ε labels.** Stored as a literal ['1e-6', '1e-8'] separate from the loop's literal; the md text "at ε = 10⁻⁶, 10⁻⁸" is a third literal (D3do-6).

**(3) Tests and mutants** [own] (W/tools_ref/mutants_bt.py is tools/mutants.py with REPO = W/repo and a per-slot `--basetemp`)
- **Shipped round117_D61_own.py.** Q00, the control, survives (39 passed). Q01–Q07 are all killed, and Q07 also changes the md. This matches REVIEW_LOG.
- **My spec (W/mut/spec120.py).**
  - Killed (7):
    - Y01: the coef loop's ε changed to 10⁻⁷ while the label stays 10⁻⁸;
    - Y03: the sign of ψy in P;
    - Y04: Ψ_φ back in float;
    - Y05: xmax0 taken from the last ψ;
    - Y06: the far check moved to 1 + 4ψ²;
    - Y09: the Re σ > 1 point moved onto y = ψ;
    - Y10: the normalisation ψ instead of ψ².

    Most fall to `part_VC() == raw['VC']`.
  - Survive (3): Y02 (`vc_lines` takes mid's ε labels from the constant) is non-equivalent but equivalent on the data. Y07 (label at P ≥ 0) and Y08 (|Re σ| → Re σ) are equivalent on the data.

**(4) The certificate's mathematics** [own]
- **The bound.** E = log|1 − z²| + Re(cz²) + C(σ) is exactly F_η(σ) − F_η(η/2). Checked: `E_centred` against `E_def` (D61's F_η with Λ⁺) agree to 9.4e-35.
- **Route (2)'s hypotheses.**
  - b = (1 + t)η/2 with t = ε/4 < (2 − |η|)/|η|; part V's chosen chord point is at 0.9 of the disc's edge, and D61's (2) allows any such t.
  - The chord crosses the saddle along the real segment, with Re A < 0.
  - The 0-side and the chord are V's exact argument.
- **The identities behind the forms** (my derivation; validated in (5)):
  - log|1 − z²| + Re z² = |z|⁴/2 − v²Λ(v)/2, with v = |1 − z²|² − 1 = −2Re z² + |z|⁴;
  - h(1 + δ) = δ²(1 + δH₁(δ));
  - |σ| − 1 = (ω²|1 + z|² − 1)/(ω|1 + z| + 1);
  - Re σ − 1 = −ε/2 + ω[X − (1 + X)θ²Cv(θ²) − Y·θ·sinc θ];
  - on the floor, Re(Az²) = −X²|A|²/|Re A|.
- **The closed forms at ρ = 0** reproduce Lemma VC's function in (u, q) coordinates; for example the cap's far value is −uq − q²/2 = (−ψ² − ψ⁴/2)/x².
- **The ε = 0 face.** At u = 0 the limit along the path is:
  - −q²/2 for q ≤ 2/3;
  - 7q²/4 − 3q + 1 for 2/3 ≤ q ≤ 1;
  - −(1 − 1/(2q))² for q ≥ 1 (the floor).

  It is negative for every q > 0, so the family is a valid valley path uniformly up to |η| = 2.

**(5) Validation** [own]
- 200 random η: ε log-uniform in [10⁻⁹, 0.005], θ drawn from the φ, ψ and fixed-θ regimes; 1479 path points (J at λ = 0.1, 0.5, 1; F; C).
- Form A (947 points): relative deviation from `E_centred` and from F_η's definition ≤ 9.3e-25.
- Form B (574 points): ≤ 7.6e-23.
- The closed forms against the forms at ρ = 0 (398 and 400 points): 1.5e-59 and 2.1e-59.
- Against the repository's own `_cert_E_ball` (vs_repo.py, 339 points): ≤ 2.3e-23.

**(6) Self-checks** [own] (W/logs/selfcheck.log, coverage.log)
- **Soundness on every accepted box.** At 9401 random points across all 3808 accepted box-branches, the true bound (from the definition, 200 bits) minus the certified bound is ≤ −7.4e-6, and the true bound is always negative.
- **Soundness on random boxes.** 60 random accepted boxes (1500 points): ≤ −1.25e-3. 60 outer segments (1200 points): ≤ −5.3e-4.
- **Coverage.** 4000 random η (ε down to 10⁻¹², θ = 0, θ down to 10⁻¹², θ within 10⁻¹² of π/2) give 28000 path points. Every one lies in an accepted closed box of its branch. At θ = 0 floor and cap coincide, and the certificate's α > b̂ rule applies.
- **Negative controls: wrong caps.** Caps at Im z = 2θ, θ/2 and 0 are rejected. Their failure is genuine: at ψ = 1, ε = 10⁻⁹, x = 1000, the bound over ε² is +997, +248 and +999, against −1.83 on the true cap; a cap at 1.1θ gives +8.2.
- **Negative controls: the straight segment.** In arg η ∈ [0.30, 0.302], 1.999 ≤ |η| ≤ 2, 128 of 400 outer segments of the straight segment fail. Again genuine: at X = 0.05 the bound is +1.96e-3.
- **Bugs caught in my own code during development, before any run was used:**
  - a NaN derivative read as 0;
  - a NaN preferred in a minimum;
  - the closed form's corner at (u, q) = (2, 0), the disc's edge at θ = 0, where floor/cap switch and D = P = 0 coincide. It needs the convexity bound D₊² − P₊² ≤ 2D₊(D − P) = 4qD₊; without it two boxes failed after 181075.

**(7) Robustness** [own]
- With ε₀ = 0.0005, 0.02 and 0.11, and with δ₀ = 0.02 and 0.1, everything passes.
- Example counts at δ₀ = 0.02: 1685/320/754 boxes and 170 η-boxes. At δ₀ = 0.1: 3286/320/811 and 67.

## The certificate (per piece)

**Settings.** W/own/lastann.py, sha256 78de9795…; run `python lastann.py all --json logs/all.json`; 128 bits; ε₀ = 0.005, δ₀ = 0.05, QA = 2.5, PB = 0.5, t = ε/4. Times are on a machine with load average about 6.

**Form A** (X = ρ, ε = uρ, θ² = w = qρ, ω = 1 − ε/2).
- Auxiliary series: S(w) = sinc √w, Cv(w) = (1 − cos √w)/w, Λ(v) = (v − log(1 + v))/v², H₁(δ) = (h(1 + δ) − δ²)/δ³, all with rigorous tails.
- Coefficients: α = u(1 − ε/4) + 2ω²qS(w)² and b̂ = 2ω²S(4w), so that |Re A| = ρα and Im A = θb̂.
- The path: Y = θm, with m = b̂/α on the floor, m = 1 on the cap, and m = λb̂/α on the junction (at u = 4).
- E/X² is the sum of:
  - Q = qm(αm − 2b̂) − αρ;
  - L = ζ²/2 − (v/X)²Λ(v)/2, where ζ = ρ + qm² and v/X = −2(ρ − qm²) + ρζ²;
  - C = D̂₊²(1 + ρD̂₊H₁(ρD̂₊)) − P̂₊², where D̂ = [ω²(2 + ζ) − u(1 − ε/4)]/(ω√(1 + ρ(2 + ζ)) + 1) and P̂ = −u/2 + ω[1 − (1 + ρ)qCv(w) − qmS(w)].
- Closed form at ρ = 0, with s = u + 2q:
  - floor: −4q/s − 8q²/s⁴ + C₀;
  - cap: q(s − 4) − q²/2 + C₀, which equals −uq − q²/2 where 1 − u/2 − 3q/2 ≥ 0;
  - junction: −(4q/s)λ(2 − λ) − 8λ⁴q²/s⁴;
  - here C₀ ≤ min(D₀₊² − P₀₊², 2D₀₊(D₀ − P₀)).

**Form B** (X = pϑ, θ² = ϑ, ε = uX).
- Coefficients: α_B = up(1 − ε/4) + 2ω²S(ϑ)², b̂ = 2ω²S(4ϑ); the path Y = θpm̃.
- E/X² is the sum of:
  - Q = m̃(α_Bm̃ − 2b̂) − ϑα_B;
  - L = ζ²/2 − (v/X)²Λ/2, where ζ = p(ϑ + m̃²) and v/X = −2p(ϑ − m̃²) + Xζ²;
  - C ≤ D̂₊²(1 + XD̂₊H₁(XD̂₊)).
- The closed form is form A's with q = 1/p.

| piece | ranges | boxes (skipped) | failed | worst bound (value route) | worst sup f′ (derivative route) | smallest sides | time |
|---|---|---|---|---|---|---|---|
| A-J (junction) | u = 4, q ∈ [0, 2.5], λ ∈ [0, 1], ρ ≤ ε₀/4 | 40 | 0 | −0.124 | −0.618 (14 boxes) | q 0.25, λ 0.25 | 0.0 s |
| A-FC (floor, cap, matching) | u ∈ [0, 4], q ∈ [0, 2.5], ρ ≤ min(δ₀, ε₀/u) | 2011 | 0 | −8.8e-14 | −1.1e-3 (86) | u, q 1/128; ρ 5.8e-13 | 1.6 s |
| B-J | u = 4, p ∈ [0, 0.5], λ, ϑ ≤ (π/2)², X ≤ δ₀, ε ≤ ε₀ | 320 (634) | 0 | −1.3e-3 | −0.028 (12) | p 1/64, λ 1/8, ϑ 3.6e-4 | 0.2 s |
| B-F (floor, θ up to π/2) | u ∈ [0, 4], p ∈ [0, 0.5], ϑ ≤ (π/2)² | 769 (1319) | 0 | −5.3e-5 | – | u 1/8, p 1/64, ϑ 3.6e-4 | 0.5 s |
| O (outer) | 1.995 ≤ \|η\| ≤ 2, 0 ≤ arg η ≤ π/2 (true π/2): connector at Re z = δ₀ plus descent polyline to \|1 − z\| < 0.01 | 75 η-boxes, 5126 segments | 0 | −1.5e-7 | mean value form in (s, \|η\|, arg η) | Δ\|η\| 0.005, Δθ 0.0061 | 0.2 s |
| E (end) | \|1 − z\| ≤ 0.02 to z = 1, R2 = 2 | 1 | 0 | −1.3314 | – | – | – |

- The "skipped" boxes lie outside ε ≤ ε₀ or X ≤ δ₀.
- The value-route worst bounds are acceptance margins of the branch and bound, not the true margins. The true minimum along the path at the handoff X = 0.05 is about −0.03·X² near θ = 0.05 and about −0.7·X² near θ = 0.4, from F_η's definition.
- At ε₀ = 0.11 (logs/all_eps011.json) the counts are 54, 2533, 413 and 1189 parameter boxes and 220 η-boxes (16082 segments, worst −4.2e-9), nothing failed, about 4 s.

## Coverage

For every η = (2 − ε)e^{iθ} with 0 < ε ≤ 0.005 and 0 ≤ θ ≤ π/2, the path is certified piece by piece. The lower edge is the float 1.995, which is exactly part V's last edge Rb; ε₀ − (2 − 1.995) > 0. The top edge is the true π/2 (ϑ ≤ up(π²/4) in form B; `th_ball` in the outer piece).

- **0-side and chord.** V's exact argument applies on the real segment up to b = (1 + ε/4)η/2.
- **J.** From b up to the floor point (t, κt), which lies below the cap:
  - in form A, α ≥ 4 − ε > 2 ≥ b̂;
  - in form B, `Bfloor_ok` checks it per box.
- **F, then C, up to X = δ₀.** Every path point with X ≤ δ₀ lies in form A (q = θ²/X ≤ 2.5) or in form B (p = X/θ² ≤ 0.5), because QA·PB = 1.25 > 1. The branch at each point is the floor if α > b̂ (⇔ κX < θ) and the cap otherwise; boxes where the branch is undecided certify both. Form B's floors never meet the cap (tan θ/θ ≥ 1), and this is checked per box.
- **O.** From the handoff δ₀ + i·min(κδ₀, θ), whose range over the η-box is enclosed rigorously, the connector and the box centre's polyline lead to |1 − z| < 0.02, checked rigorously.
- **E.** V's end piece completes the path.
- **The rest of the annulus.** The reflections η ↦ η̄ (z ↦ z̄) and η ↦ −η (C depends on |σ| and |Re σ| only) give the full annulus.
- **The disc.** With V (|η| ≤ 1.995) and D's branch and bound (|η| ≤ 1.8927), or with my run at ε₀ = 0.11 in place of V, the valley condition holds on the whole open disc |η| < 2.

**So the pieces cover 1.995 ≤ |η| < 2 completely.** What is not certified, or remains to be done:
- |η| = 2 itself, which is not needed;
- an independent second derivation of the forms' algebra (mine, validated numerically as in (5));
- adoption into the repository, with tests.

## Status of the round-117 issues

| id | status |
|---|---|
| D3dl-1 (junction) | **Resolved** in wording (D61's Status, `vc_lines`, HANDOFF §5) [md]. Now certified: A-J and B-J, through f(0) ∝ λ and sup f′ ≈ −(4 + 2q) (D3do-1). |
| D3dl-2 (one family) | **Resolved.** Residue: the stale "per regime … overlaps" sentence (D3do-4) and "exactly y = min(sx, ψ)" (D3do-5). |
| D3dl-3 (ψ chart certificate) | **Deferred** at a2f3b77 ("not in the repository"). **Superseded** by D3do-1, which covers the whole annulus rather than the ψ chart's inner part. |
| D3dl-4 (test gaps) | **Resolved.** Q01–Q03 are killed; the grid is stored and pinned; the Re σ > 1 point is checked and confirmed from the definition. Residue: no VC_MID point on the cap's far piece (y = ψ, x ≥ ½ + 3ψ²/2), the second half of round 117's proposal. |
| D3dl-5 (md wording) | **Resolved:** "0.0001 ≤ x", "145 ψ", "−0.0043", "D₊²", the lemma's place, Ψ_φ in balls, mid's ε from the raw. Residue: coef's ε labels are literals; Y02 survives (D3do-6). |
| D3dl-6 (residue) | **Partly resolved.** Fixed: "enclosed sample boxes" in the md and HANDOFF, the matching condition, and HANDOFF §7's φ₀/φ₁ sentence. Still open: HANDOFF §7's "Write it as D‴.61 (iii)(b) on the whole disc, add the connectivity check to part D (η up to 1.99, the margin)" with its scipy.ndimage plan, and D61's Numerics (VC), which were not extended (deferred, as logged). |

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3do-1 | major (constructive) | **The last annulus is certified.** Every η with 1.995 ≤ \|η\| < 2, 0 ≤ arg η ≤ π/2 has an explicit path (b → J → F → C → O → E) with E < 0 in Arb: 3140 parameter boxes, 75 η-boxes, 5126 segments, none failed, about 3 s [own]. Validated against `cert_E`, F_η's definition and `_cert_E_ball` (≤ 7.6e-23). Sound at 9401 points sampled in every accepted box (true − bound ≤ −7.4e-6). Coverage complete at 28000 sampled path points. Wrong paths rejected. D61 still says "a conjecture on 1.995 < \|η\| < 2" and makes (iv) beyond r = 0.99 conditional on P [md]. | Adopt W/own/lastann.py as a part (say VA) of exp_gh_peano.py. Store the per-piece counts and bounds in the raw and add an md line. Tests: the forms against `cert_E` at thin points; the closed forms against ρ = 0; a full rerun (about 3 s); the soundness sampling; two negative controls (the cap at 2θ; the straight segment in the sector). Then D61: (iii)(b) on \|η\| < 2 (computer-assisted); (iv) for every frame r < 1; the half-squeeze rule's sharp constant for every σ*/σ; Status and Numerics. Update PAPERS' claim and HANDOFF §§5, 7, 8. Have one referee re-derive the identities listed above. |
| D3do-2 | minor | **Two charts suffice.** "Four parameter charts … and in z an inner, a matching and an outer chart" (D61, md, HANDOFF §5) overstates the work. In (u, q) = (ε/X, θ²/X) the regimes are faces of one square: q = 0 for φ and (ε/θ, θ/√ε); the corner for matching; q → ∞, via form B, for (√ε/θ, θ). The next orders (Ψ_φ, −X³/3, V's ε³) appear as sup ∂f/∂X on the degenerate faces [own]. | Replace the plan in D61's "What remains", `vc_lines`' last sentence and HANDOFF §5 by the two-form description, or by "certified (part VA)". |
| D3do-3 | minor | **Part V is subsumed.** At ε₀ = 0.11 the same code certifies 1.89 ≤ \|η\| < 2 (4189 parameter boxes, 220 η-boxes, about 4 s), against V's 13429 η-boxes and 10.7·10⁶ segments [own][md]. | Keep V as an independent cross-check, or replace it, stating the overlap with \|η\| ≤ 1.8927. |
| D3do-4 | cosmetic | D61 line 218 still reads "A proof needs, per regime, its path certified …, the regimes' overlaps, and inner, matching and outer charts in z", contradicting the preceding "no matching of paths needed" [md]. | Delete it in D61's Status. |
| D3do-5 | cosmetic | "In Ψ̃'s coordinates exactly y = min(sx, ψ)" (D61, `vc_lines`) holds only as ε → 0. At ε > 0 the floor's slope is A_i/\|A_r\| [own]. | Write "to leading order (the cap y = ψ exactly)". |
| D3do-6 | cosmetic | **Round-117 residue** [md][own]: (i) the coef's ε labels: the md text literal "at ε = 10⁻⁶, 10⁻⁸" and the stored list are not the loop's ε; (ii) Y02 survives (no test feeds `vc_lines` a raw whose ε differs from the constant); (iii) no VC_MID point lies on y = ψ, x ≥ x₂; (iv) "(0.2, 3, 0.1, Re σ > 1)" puts the label inside the tuple; (v) D61's Numerics (VC) lack rounds 115 and 117's numbers; (vi) HANDOFF §7's whole-disc/connectivity plan is stale. | (i) one ε tuple for the loop, the raw and the md; (ii) a test with altered ε lists; (iii) add e.g. (0.2, 3, 0.2); (iv) "…: (0.2, 3, 0.1), Re σ > 1"; (v) extend Numerics (VC); (vi) drop the HANDOFF §7 plan. |

## Sources

- review/drafts/D61.md, results/chaos/gh_peano.md, review/REVIEW_round117_D3dl.md §§6–7, review/REVIEW_round115_D3dj.md §3, HANDOFF.md, REVIEW_LOG.md [md]; gh_peano_raw.npz [raw].
- Johansson, Arb: efficient arbitrary-precision midpoint-radius interval arithmetic, IEEE Trans. Comput. 66 (2017), arXiv:1611.02831 [memory; as cited in D61].
- R. E. Moore, Interval Analysis (1966): centred and mean value forms [memory].
- L. B. Rall, Automatic Differentiation (1981): forward-mode differentiation in interval arithmetic [memory].
- W. Tucker, Validated Numerics (Princeton 2011) [memory].
- Weighted blow-up and distinguished limits (Kruskal–Newton diagram; arXiv:1807.03107), as cited by round 117 [snippet, via round 117].
- No web search this round.

## Lesson

1. **Blow up in the path point's weight-0 coordinates, not in the parameters.** With (u, q) = (ε/X, θ²/X), four parameter regimes and the inner and matching charts become faces and a corner of one compact square. The radius X then becomes the variable of a mean value form whose f(0) is the limit in closed form.
2. **On a face where the limit vanishes, certify the derivative, not the value.** f(0) ≤ 0 with sup f′ < 0 gives f < 0 for X > 0. The blow-up's next order is just ∂f/∂X there, enclosed automatically by forward differentiation in balls.
3. **Three guards made this sound.**
   - Evaluate closed forms near such faces in endpoint intervals, because midpoint-radius balls lose the sign of products of intervals at 0.
   - Use the convexity bound a₊² − b₊² ≤ 2a₊(a − b) where two positive parts switch together.
   - Never let a NaN pass as 0.

   Then check the certificate itself by sampling the true bound inside its own accepted boxes, by sampling coverage, and by negative controls on wrong paths.

**Files.** All under /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref120:
- own/: lastann.py, selfcheck.py, coverage.py, vs_repo.py, README.md;
- logs/: all.json, all_eps011.json, run_all.log, run_eps011.log, selfcheck.log, coverage.log, vs_repo.log, pytest_gh.log;
- mut/: spec120.py, log_ship.txt, log_own.txt, ship/res_a.jsonl, own/res_b.jsonl;
- tools_ref/mutants_bt.py;
- out/ (the regenerated md and the raw copies).
