# Handoff to the next session (2026-11-01, second session)

This note lets a new session carry the research on without the previous sessions' context. The scratch areas are gone; everything needed
is in the repository. Read this file first, then `review/drafts/README.md`, then the last entries of `review/REVIEW_LOG.md` (Rounds 84,
90, 85, for the format), then the referee reports named in §5–§7.

## 0. First commands

```bash
cd /home/user/Lyapunov-Convergence-Kernels-Showcase
git status -sb                      # your designated branch; start it from origin/main if it is new (§2)
git log --oneline -5                # 5e42eee (round 85 closed), 4bc7b07 (rounds 84 and 90), …
python3 -m venv venv && venv/bin/pip install -r tools/requirements-pinned.txt     # venv/ is not in git; the exact versions of the raws
OMP_NUM_THREADS=1 venv/bin/python -m pytest -q -p no:cacheprovider tests/test_width_mismatch.py tests/test_gh_peano.py   # 63 + 42 pass
OMP_NUM_THREADS=1 venv/bin/python -m pytest -q -p no:cacheprovider tests/test_noise_resonance.py                         # 71 pass
```

Make a scratch directory (the session's scratchpad) and call it `$S` below. Put every run directory, patch script and log there, never in
the repository.

## 1. The project in one page

The repository started from one idea: kernels that adapt their own centres, widths and weights inside a control loop, with a Lyapunov
function proving convergence. It has grown into a research program in seven threads (status: proved / partly proved / numerical /
conjectural):

| thread | results | status |
|---|---|---|
| Convergence foundations | E′ (explicit strict Lyapunov function under persistent excitation), E.0′, E.0″ (sharp rates, non-periodic excitation), I (non-reciprocity), J, J′ (finite-m averaging, stationary law under dither) | proved |
| Mean-field limit | D′, D″, D‴: propagation of chaos with a plant in the loop, finite horizon and uniform in time | proved |
| Leaks and priors | C′, D‴.9–D‴.31: optimal-transport leaks, the mass constraint, a PL basin uniform in m (D‴.21), adaptive priors, measurement noise | mostly proved |
| Noise selects kernels | D‴.32–D‴.47: the gauge under dither (Katzenberger drift, Fixman correction), the colour of noise as a sparsity prior, the free law as a mirror flow | proofs + numerics |
| Sparse kernel formation ("rich transit") | D‴.48–D‴.60: one kernel per feature from small masses; sliding limit (Theorem D‴.53), exact recovery (Theorem D‴.55); close features, ℝⁿ, noise, dead zones, width mismatch, moving widths, noise resonance | D‴.53, D‴.55 proved; the rest partly proved + numerics |
| Gaussian quadrature of features (new) | D‴.58 (a wider feature as a sequence of Gauss–Hermite rules; the half-squeeze rule's exact series; one formula t(v) for every scaled rule's rate, proved for every rule with a Mehler frame and, by D‴.62's Bochner representation on one axis, for the rules without one) and D‴.61 (the Gauss–Hermite Peano kernel as a law; conjecture P proved on the whole open disc |η| < 2, three computer-assisted steps rigorous in ball arithmetic — D's branch and bound, V, VA (the last annulus, rounds 120–121) —, which gives D‴.58's sharp constants for every frame; part VC the corner's blow-up limit; since round 148 also proved by hand, the three certificates now cross-checks), and D‴.62 (the rules without a Mehler frame: the rate and sharp constant, c₂ = s₁/2 proved since round 139) | proved (D‴.61 by hand and computer-assisted; D‴.62 by hand); D‴.58's other parts as its Status states |
| Applications | cogging-torque servo, second-order plants, NN tracking control, mixture-prior compression | simulations |

Where things live:

- `PARADIGM.md` — every result through D‴.48. `results/RESULTS_PARADIGM.md`, `README.md` (summary table), `PAPERS.md` (18 paper plans:
  7 core, 10 applied, and Paper 18 on the scaled Gauss–Hermite rules), `LITERATURE_TRANSPORT.md`, `ASSESSMENT.md`.
- `manuscript/` — LaTeX (`lyapunov.tex`, `sections/*.tex`), built to `lyapunov.pdf`; it covers every result through D‴.62 (D‴.40, D‴.44, D‴.49–D‴.62 added in manuscript part 5, 2026-11-06, with the
  dated corrections to D‴.36–D‴.48; the drafts transcribed from `review/drafts/` at 8322a5d, so a later change to a draft must be ported to its
  `sections/thmDppp_NN.tex`). Build: pdflatex three times in a copy of `manuscript/` beside a `results` link; needs the `cm-super` fonts.
- `review/drafts/Dxx.md` — the newer propositions, refereed but not yet moved into `PARADIGM.md`: D‴.44, 49, 50, 52, 54, 56, 57, 58, 59,
  60, 61. `review/drafts/README.md` has one status row per draft.
- `review/REVIEW_roundNN_<id>.md` — the referee reports; `review/REVIEW_LOG.md` — the append-only log, one entry per applied round.
- `src/core/` — library code; `experiments/exp_<name>.py` — each writes `results/chaos/<name>.md` and `<name>_raw.npz`;
  `tests/test_<name>.py`; `tools/` — the process tools (§4), `tools/mutants/` the mutant specs.

The drafts and their experiments (output and worker env vars):

| draft | experiment | env vars |
|---|---|---|
| D‴.44 | `exp_mass_dictionary.py` | MD_OUT, MD_WORKERS |
| D‴.49 | `exp_closed_transit.py` | CT_OUT, CT_WORKERS |
| D‴.50 | `exp_transit_noise.py` | TN_OUT, TN_WORKERS |
| D‴.52 | `exp_close_features.py` | CF_OUT, CF_WORKERS |
| D‴.54 | `exp_transit_nd.py` | TN_OUT |
| D‴.56 | `exp_noisy_certificate.py` | NC_OUT, NC_WORKERS |
| D‴.57 | `exp_dead_zone.py` | DZ_OUT, DZ_WORKERS |
| D‴.58 | `exp_width_mismatch.py` | WM_OUT, WM_WORKERS |
| D‴.59 | `exp_moving_width.py` | MW_OUT, MW_WORKERS |
| D‴.60 | `exp_noise_resonance.py` | NR_OUT, NR_WORKERS |
| D‴.61 | `exp_gh_peano.py` (parts P, M, L, LC, S, R, Q, D, V, DB, VC, VA — VA runs `gh_lastann.py`; R reads D‴.58's raw) | GP_OUT, GP_WORKERS |
| the transit across readouts (depth 1, 2, 3; fixed-step descent; Adam; branch pairs; dictionary spacing) | `exp_parametrization.py` (parts A–L; C, J, K and L are long runs: with `PZ_WORKERS=2` on the shared 4-core machine K took 11 min and L 5.3 min of wall-clock time (the cells' own times add up to 1,266 s and 608 s: the 255 s and 235 s quoted before the PZ2 round were such sums for the smaller K and L, not wall-clock times), D 1.3 min; engine `src/core/factorized_rbf.py`) | PZ_OUT, PZ_WORKERS |

Every experiment takes `--only PART,PART` (run those parts, keep the rest of the raw) and `--from-raw` (regenerate the md from the raw
without computing). Raws are saved and loaded with `exp_cohesion_leak.save` / `exp_cohesion_leak.load`; loaded lists are numpy object
arrays.

## 2. Conventions (all in force)

Branch and commits:
- Develop, commit and push only on the branch your session's instructions designate (this session's: `claude/vigilant-clarke-hf939b`;
  the earlier `claude/festive-bardeen-eopej2` was merged into `main` as PR #1 and is finished): `git push -u origin <branch>`, retried up
  to 4 times with 2/4/8/16 s backoff on network failure. If the designated branch is new, start it from `origin/main`, which holds all the
  research. No pull request unless the user asks.
- Commit with `git commit -F <msgfile>`, adding paths explicitly (never `git add -A`). End the message with the attribution lines your
  session's system prompt gives. Never rewrite history. No model identifiers in repository files.
- A container can restart and kill background jobs, or be reclaimed with uncommitted work. Commit at the end of each round; if pausing
  mid-round, commit a clearly labelled work-in-progress commit (as a4bfe1e).
- Optional: while a round is mid-edit, its WIP files may be hidden with `git update-index --assume-unchanged` (log each flag and unflag
  in `$S/LOCAL_GIT_FLAGS.md`), so that the tree is clean between turns; unflag before committing.

Results:
- Every number comes from a script. A draft quotes the md's printed numbers verbatim. No arithmetic on md numbers by hand, including
  roundings, ratios, and min/max over printed lists: make the script print it. `tools/numcheck.py` lists the numbers a draft's added
  lines quote that the md does not print.
- `--from-raw` must regenerate the md byte-identically: run it twice and `cmp`.
- Every new or changed result line or function gets tests; pin case lists and constants against their stated contents and the raw, and use
  `np.isclose(..., atol=0)` on small quantities (the default atol 1e-8 makes assertions on ~1e-13 vacuous). Each round runs single-point
  mutants of its new code (`tools/mutants.py`): kill each survivor with a test, or record why it is equivalent on the data.
- Every new result or revision goes to an adversarial referee subagent (§4). It must not modify the repository and returns its report as
  its final message. Archive it as `review/REVIEW_roundNN_<id>.md` with `tools/extract_report.py`.
- `review/REVIEW_LOG.md` is append-only. In-place bracketed errata are allowed only when a referee asks for them.

Processes:
- Python as `OMP_NUM_THREADS=1 nice venv/bin/python ...`; pytest as `venv/bin/python -m pytest -q -p no:cacheprovider`.
- 4 cores, 16 GB: at most about three heavy processes at once.
- Two processes must never write the same raw: each run gets its own output dir (`WM_OUT=$S/runX`, …) holding a copy of the repo's raw;
  merge afterwards with `tools/merge_parts.py`.
- Background jobs as separate statements `( cmd; echo "<label> finished rc=$? $(date)" >> $S/queue.log ) > /dev/null 2>&1 &`. Never end an
  `&&` chain with `&`. Watch the queue with a Monitor:
  `tail -n0 -F $S/queue.log | grep -E --line-buffered "finished|done|rc=|Traceback|Error|Killed"` (30 min; re-arm on expiry).
- No foreground `sleep`. Never `pkill -f`/`pgrep -f`; kill only explicit PIDs verified through `/proc/PID/cmdline`.

Editing:
- Edit atomically: write `.tmp`, `py_compile`, `mv`. Read a file before editing it. Python patch scripts (exact `old → new` replacements,
  each asserted to occur once) in `$S` are the usual way to change the large experiment files.
- Raw-loaded lists are numpy object arrays: never `x or []`, never `if not arr`.
- Quoted heredocs (`<<'EOF'`). In Python string patches beware `'''` adjacency and `\\n` double-unescaping (put test snippets in files).
- Escape `|` in markdown table cells as `\|` (`\\|` inside Python string sources). Avoid `rev` on multibyte text.
- Drafts are wrapped at width 140. For a changed paragraph: join it, replace, then
  `textwrap.fill(width=140, break_long_words=False, break_on_hyphens=False)`; rewrap only the changed blocks. Never use the file's max line
  length as the width. (`tools/rewrap.py` rewraps a whole file's prose; prefer per-block rewrapping.)

Other:
- arxiv.org, ams.org and ar5iv.labs.arxiv.org are unreachable from these containers; tag literature seen only in search results
  `[snippet]`. When citing web sources to the user, give markdown links.
- mpmath and sympy are not installed; exact arithmetic uses `decimal` and `fractions`.
- Use they/them for people.

## 3. What the user asked for

The standing request (verbatim, abridged): "Analyze the repo and weigh its integral principles, contributions and its impact and
novelty. I want to continue the research done here, making this work more and more impactful and novel, while ensuring that the work you
do does not deviate far at any point from the central topic. Keep yourself aware at all times of how novel the direction you're going in
is, and at all times ensure that whatever you are coding is also mathematically integral and the foundational rigour of it is strong.
Keep pushing on and on … if there are hurdles or … negative results … analyze the faults and come up with explanations and strong
solutions … make this research foundation something on which a very high number of different papers can be written in future, both
application based and core contribution papers … keep me updated as you go about your work." Give the user short status lines while
working.

## 4. How a round works

1. Read the report (`review/REVIEW_roundNN_<id>.md`): its Summary, its Issues table (id, severity, evidence, proposed fix) and its
   Lesson.
2. Change the code (new parts, new report lines), each with tests.
3. Run the new parts in private output dirs (§2), merge them, regenerate the md twice and `cmp`.
4. Run the experiment's full test file. Run mutants on the round's new code:
   ```bash
   venv/bin/python tools/mutants.py --exp exp_width_mismatch --tests tests/test_width_mismatch.py \
     --raw results/chaos/width_mismatch_raw.npz --md results/chaos/width_mismatch.md --out-env WM_OUT --workers-env WM_WORKERS \
     --spec tools/mutants/round84_D58_own.py --work $S/mut --slot a --i0 0 --step 3      # and slots b, c with --i0 1, 2
   ```
   Each slot is a private copy; results land in `$S/mut/res_<slot>.jsonl`. A spec is a list `M` of `(id, what, old, new, k)` tuples. Each
   run gets the repository's md, so a test may compare md and raw. To rerun survivors against new tests: `--only ID,ID,...` in a new slot.
5. Revise the draft: numbers verbatim from the md (`tools/numcheck.py <draft> <md> <rev>`), rewrapped at 140; update the draft's row in
   `review/drafts/README.md`.
6. Append the round's entry to `review/REVIEW_LOG.md`, in the format of the Round 83 entry: a heading
   `## Round NN (adversarial review of Proposition D‴.XX as revised after round MM, <date>)`, then **Design.** (what the referee did),
   **Verdicts (<id>).**, the issues table `| id | sev. | issue | fix (where) |` with the fix as applied, and **Lesson.**
7. Commit, push, and launch the next referee round on the revision; archive its report when it arrives, and commit that.

Referee ids run … D3ck (90, D‴.61), D3cl (91, D‴.61), D3cm (92, D‴.58), D3cn (93, D‴.61), D3co (94, D‴.60); the next is **D3cp = round 95**,
and so on.

A referee brief (a general-purpose subagent; adapt the names):

> You are an adversarial referee, report D3ci, round 88, for Proposition D‴.58 of the repository at
> /home/user/Lyapunov-Convergence-Kernels-Showcase (branch <your branch>, commit <hash>). Do not modify the repository in
> any way: export it with `git archive <hash> | tar -x -C W/repo` into your own work directory W under <scratchpad>/ref88, and run
> everything there (WM_OUT pointing into W; `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 10`, the repo's venv, at most two heavy
> processes). Read in full: review/drafts/D58.md, experiments/exp_width_mismatch.py, results/chaos/width_mismatch.md,
> tests/test_width_mismatch.py, review/REVIEW_round84_D3ce.md and the Round 84 entry of review/REVIEW_LOG.md, and
> `git diff f710a0d <hash>` for these files. Check: reproducibility (the md byte for byte, the tests, a sample of jobs bit-identically
> against the raw); every claim of the draft against the md and the raw (numbers verbatim, no overreach); the mathematics (what is proved,
> what is numerical, what is conjectural; hypotheses); the code, recomputing key quantities with your own code; the tests, with
> single-point mutants of this round's new or changed code (report survivors, and which change published results); the literature (web
> search; cite only what you saw); the status of round 84's issues. Tag evidence [md], [raw], [own], [snippet], [memory]. Start the
> report with "# Referee report D3ci, round 88: Proposition D‴.58 as revised after round 84" and use the sections Scope, Summary (design,
> verdicts, numbered findings D3ci-1, … with severity major/minor/cosmetic), What I checked and how, Verdicts per part, Status of the
> round-84 issues, Issues (a table: id | sev. | issue and evidence | proposed fix and where), Sources, Lesson. Your final message must be
> the full report and nothing else.

Archive it: `venv/bin/python tools/extract_report.py <the agent's output file> review/REVIEW_round88_D3ci.md "# Referee report D3ci, round 88"`.
If the agent's final message is available directly, writing it to the file is equivalent.


## 5. In flight at handoff

**Latest (2026-11-02, end of session):** D‴.59 round 133 applied (a80663c: the scaling as an interval, a real prefactor
difference between the arithmetics named); D‴.60 round 134 applied (10e4244: C8F, the full 2⁵ separation — C8A closed); Theorem 3's
certificates checkable by `tools/verify_certificates.py` from `results/certificates/` (b5c6fa3, 13 s); the Paper 18 manuscript revised
after round 135 (4a005cb, 085250b). Since then: KOG21, KS19, GRR24 and others read in full (a0b250a, `review/LITERATURE_gauss_hermite_rkhs.md`; one novelty claim
corrected: the classical rule's rate λ was observed numerically in KOG21); TeX Live installed (apt) and the manuscript compiles cleanly
(f4825c8, 35 pages, `manuscript/paper18/paper18.pdf`). A second journal-style round (136, D3ee) recommended minor revision; applied in 638fb1b (paper 35 pp. + supplement 4 pp.; the verifier the certificate of record). Open: KW12 and KSW17
are not on arXiv and their pages are unreachable (allow web.maths.unsw.edu.au); referee rounds on the latest D‴.59 and D‴.60 revisions.

Reports received and being applied by subagents (each then needs: verify, log entry, README row, commit, next referee round); if an
application did not reach the repository, redo it from the archived report; the next referee id is D3ee = 136:
- **Round 128 on D‴.60** (D3dw) — applied (a15cc7e): the brackets rounded outward, the matched arrival at every window, 'cold_c' (the
  centre necessary and sufficient at the 11 alarms), the kill plants, the labels; a referee round on this revision is due.
- **Round 131 on D‴.61** (D3dz) — applied (f3b22eb): VA011 (VA at ε₀ = 0.11) makes the proof's chain D + VA, V a cross-check;
  Appendix A states every certified inequality as a lemma (it found the series tails' underflow, now floored); seam completeness;
  citations. D3dz-2 (split into Theorems 1–3 and Corollary 4, the statement under 400 words, no referee rounds in the mathematics,
  one path lemma) is for the Paper 18 manuscript; D58's assumption sentence (D3dz-8) for D58's next revision.
- **Round 129 on D‴.58** (D3dx) — applied (df069e9): A6C (exact errors at every K = 159 … 400) shows the Mercer frames' stall is the
  window; the transient carried only above an order's remainder; c₂ derived (formal), the check with it fixed passes at all three proved
  rules. Pending for D58's next revision: cite D‴.62 for the rules without a frame (after round 132's application) and round 131's
  D3dz-8 (D58's assumption sentence aligned with D61's list).
- **Round 130 on D‴.59** (D3dy) — applied (e851e9b): A2P on 800 fresh draws, Λ/(11 ln 2) = 1.013 ± 0.009; A2H, the prefactor's
  step dependence (shared in part); computed verdicts; regeneration tests. A referee round on this revision is due.
- **D‴.61's older residue** — round 131 found the unit-test targets of round 121's G01, G03, G04, G08 correct (0 failures in its own
  checks) and VC outside the proof's chain; they remain regression guards only.
- **D‴.62 — the rules without a frame, proved and refereed** (8ed0f31; round 132 applied in 846c023). A Bochner representation
  e_K² = σ*(v/τ²)E[φ_K(a + b)φ_K(a − b)] puts the representing Gaussian on one axis exactly when |ρ| ≥ 1, so D‴.61 (iii)(b) is needed
  only on the axes' intervals (proved by hand): D‴.58's rate t(v) and a sharp constant C(v) for moment matching and every |ρ| ≥ 1 rule,
  no computer assistance (moment matching: e_K ~ 2^{−1/4}(1 + λ)^{3/4}λ^K); two ridges that add; c₂ formally (proved since round 139:
  the c₂ note below). Round 132's independent
  attack found no gap; it found that the result contradicts Kuo–Woźniakowski 2012's threshold 2αγ² = 1 for the classical rule (moment
  matching converges at every length-scale), stated conditionally (papers seen through abstracts). Propagated into D58 (i), (ii) and
  its Status (with round 131's assumption list), D61 (iv)'s last sentence, PAPERS.md Paper 18 and width_mismatch.md's labels.

**D‴.61's conjecture P is proved on the whole open disc |η| < 2 — by hand since round 148** (the research note
`review/RESEARCH_note_annulus.md`: the valley condition on 1 < |η| < 2 along the three-segment path Π_η, chord to the unit circle, then
vertical, then horizontal, on which the bound beyond the circle is monotone; refereed in round 148, D3eq: "the proof stands"; applied to
D‴.61 with a dated upgrade note and propagated to Paper 18 — no computer-assisted step left, its Appendix A moved to the supplement, 28 pp.
—, D58, D62, the c₂ note, PAPERS.md and the scripts' report text; D, V, VA and VA011 now independent cross-checks). Until then it was
computer-assisted: round 120's referee certified the valley condition
on the last annulus 1.995 ≤ |η| < 2 along part VC's path family Im z = min(κRe z, arg η) (two scaled forms of E/X², exact on the path up to two one-sided simplifications, in the path
point's weight-0 coordinates (u, q, ρ) = (ε/X, θ²/X, X) and (u, p, ϑ) = (ε/X, X/θ², θ²), mean value forms in Arb, an outer polyline per
η-box, V's end piece), adopted as experiments/gh_lastann.py (part VA). Round 121's referee verified it independently and found one
soundness defect — in the η-box nearest η = 2 Re A's ball contains 0, κ's ball is NaN, and the connector was silently dropped (round 120's
certificate was not a proof as submitted) — and fixed it; with the fix: 3140 parameter boxes, 76 η-boxes, 5193 segments, none failed;
rounds 122 and 126 then made the checks check what was certified: the joint rebuilt from the certified pieces (Y(η) proved inside every
connector on 484 sub-boxes, the pieces chaining the connector and the polyline to the end piece), the η-boxes disjoint and of exactly the
sector's area (in rationals), the forms' value below the bound at 29 056 distinct points (35 500 point–branch evaluations); round 127 reread the
whole proof (no gap between the certificates) and had Duhamel's identity proved in the text (telescoping on the banded Jacobi matrices). Hence (iii)(b) on |η| < 2, uniformly on compact subsets;
(iv) for every frame with r < 1; D‴.58's t(v) for every rule with a Mehler frame and the half-squeeze sharp constant for every σ*/σ.
Still not proved (then): the rules without a frame (|ρ| ≥ 1, moment matching; D‴.62 since), and an analytic (not computer-assisted) proof
of the corner (the annulus by hand since round 148, along Π_η, which avoids the corner).
New in round 113 (worth a paper remark): A6X's outlier (τ/σ = ¼, 2·v_hs, no frame, ρ = 16.5) carries a geometric transient at
q = (h(τ²/v)/t(v))², the formula's other branch, as if the two arguments' contributions to e_K² added — numerical evidence for
the two-axis (Fock-integral) picture beyond the Mehler frames; a proof would extend D‴.61 (iv) to |ρ| > 1.

- **Paper 18 manuscript** — `manuscript/paper18/` (835d9a0; revised after the journal-style referee round 135, D3ed, major revision:
  4a005cb, with `response_round135.md`). The mathematics was checked correct; remaining: (I1) read Kuo–Woźniakowski 2012,
  Kuo–Sloan–Woźniakowski 2017, Karvonen–Särkkä 2019 and Karvonen–Oates–Girolami 2021 in full — arxiv.org and the publishers' pages
  are refused by this environment's network policy (the user must allow them) —, and (I4) a standalone certificate verifier reading
  stored data (being built by a subagent: tools/verify_certificates.py, results/certificates/). Not compiled (no pdflatex here).

- **Optional tracks started (2026-11-02, the user's go-ahead):** referee rounds 137 (D‴.59) and 138 (D‴.60 as a whole); a research
  note proving (F) to every order and c₂ = s₁/2 for the frameless rules (`review/RESEARCH_note_c2.md`, 9398242; referee round 139,
  D3eh: accepted with minor revisions, applied and propagated — c₂ = s₁/2 proved for every rule without a frame (Proposition A at n = 3,
  Proposition B, Corollary B′), F₂'s closed form numerical and not needed, the frame rules (Remark F, a sketch) and the cusp still formal —
  to D‴.62 (d) with a dated correction, D‴.58 (ii), Paper 18 (open problem (4), Table 1, §5, supplement §S1) and PAPERS.md §18); the Paper 17 draft (`manuscript/paper17/`, 9398242; round 140 (D3ei, major revision) applied in d6994d3 — 39 pp., still to shorten, mutants on its new Parts U, W, G, H not run; D‴.48/51/53/55 corrected in PARADIGM.md, a dated block before D‴.55); the Paper 16 draft (`manuscript/paper16/`, a833d10; round 141 (D3ej, major revision) applied in 58e68c6 — 25 + 7 pp.; D‴.45 (iv)(b) and D‴.47 (v) corrected inline); the Paper 15 draft (`manuscript/paper15/`, e82680b; round 142 (D3ek, major revision) applied in 81976da — Theorem 3.4 (b) proved for Gaussian noise of rational spectrum by state augmentation, 26 + 11 pp., main text 23 pp. against a 20-page target; D‴.41, D‴.43 corrected); the Paper 14 draft (`manuscript/paper14/`, 1de112e; round 143 (D3el, major revision) applied in f501b9f — 32 + 8 pp., to shorten; D‴.36–D‴.39, D‴.42 corrected in a dated block before D‴.40). Since then (2026-11-05): second referee rounds on all four — 144 (Paper 16), 145 (Paper 15), 149 (Paper 17), 150 (Paper 14), each minor revision, all applied; Paper 17 positioned as a proof of AGF's open gradient-flow correspondence (Kunin et al., round 149), with its main proofs in Appendix A (34 pp.); Paper 15's C² assumption proved (Lemma 2.6); the declined runs done (Paper 14 Parts J, G with src/core/gauge_batch.py; Paper 17 Parts I, J, K); every new md part mutation-tested; D‴.44 (round 146) and D‴.40 (round 147: white noise in the unconstrained chart; colour, Stratonovich and the physical chart as a remark) re-reviewed and applied; the valley condition proved by hand on the whole disc (review/RESEARCH_note_annulus.md, round 148; certificates D, V, VA now cross-checks). Then (2026-11-05): D‴.40's early-phase gap resolved as a fluctuation (Part I, 192 fresh seeds); Paper 14 gained Section 7 (white noise of any covariance in the unconstrained chart) and a third round (153: Part J, the Euler bias halves with the step; Part K, the long runs' offsets decomposed); Paper 17's plane computations (Parts L, M: the least certificate at truncation 10⁻¹³, the plane BLASSO) and a third round (152: the AGF claim scoped — an AGF-type correspondence beyond AGF's setting); Paper 15's third round (151: Lemma 2.6 correct). The O(γ²τ)-drift question then answered (review/RESEARCH_note_gamma2_drift.md, gauge_noise.md Part M): the loop carries Katzenberger's drift at about 0.8 of the prediction (Richardson over the step); Part K's shortfall was mainly its first-order offset (the fit's motion with the masses, an (I − M)⁻¹ factor); the remaining ~20% is real (fresh seeds, Part N) and a finite-τ correction is favoured but not established (slope z ≈ −2, τ → 0 limit consistent with 1) — settling it needs several times Part N's seeds at τ = 3·10⁻⁴, a lower-variance estimator, or the loop at smaller τ; Papers 14–17 could go to a fourth round or to submission.

## 6. What this session did (2026-11-01)

- **D‴.58, round 84 closed.** A6 run to K = 160: the conjectured rate t(v) = tanh((r + |ln(v/v_hs)|)/2) holds within 6.1·10⁻⁵ at all 32
  rules off the cusp; the τ²/v branch converges more slowly (excess as K^{−1.52…−1.98} at K = 160, steepening) but keeps falling; the cusp
  within 3.1·10⁻⁷ with e_K ∝ K^{1/4}t^K. Two failing tests fixed; seven vacuous np.isclose assertions made exact; mutation testing (157
  applicable; 14 survived; 4 killed since; 10 recorded).
- **D‴.60, round 85 closed.** C7L and C8N run; the C9A line's mispairing of leads with runs found and fixed; the draft revised finding by
  finding (the ejecting set window-dependent; C8 a replay with nulls showing a leave needs both the drifted coefficient and the coming
  noise); mutation testing (109 mutants; 41 survived; 29 killed since; 12 recorded).
- **D‴.61, new, refereed (round 90) and revised.** From the attempt to prove D‴.58's conjecture: the Gauss–Hermite Peano kernel is the law
  of a Dirichlet average of the doubled nodes and a squared-Hermite point (proved); its moments are weighted Dyck paths with Var =
  (K² + K + 1)/((K + 1)(2K + 1)) (proved); conjecture P, ĝ_K(√Kη)e^{−Kη²/4} → (1 − η²/4)^{1/2} for complex |η| < 2 (the exact cumulants to
  K = 5120 match −(2j − 1)!/4^j for j ≤ 6, to 2.5·10⁻⁷ … 1.9·10⁻³ after extrapolation; the transform on and off the axes to K = 160); and,
  proved conditionally on P, the rate t(v) of every scaled rule with a Mehler frame (22 of A6's 36 rules) and its sharp constant —
  including D‴.58 (i)'s conjectured half-squeeze constant with ₂F₁(−½, −½; 1; t²) —, within 0.01 of the exact errors at 16 of 22 rules at
  K = 160 (the others on the real axis, still converging). It explains D‴.58's slow branch: the τ²/v rules have their mass on the real axis,
  where P converges slowest. PAPERS.md Paper 18.
- **Process.** The branch moved (§2). `tools/mutants.py` now gives each run the repository's md.

- **Later the same day (second context).** Rounds 88 (D‴.58), 89 (D‴.60) and 91 (D‴.61) applied, each with a log entry and mutation
  testing. The main advance: **conjecture Q is proved and conjecture P on most of the disc** (D‴.61 (iii)). The steps:
  - an exact representation He_K(√K(ζ + 1/ζ)) = K^{K/2}ζ^{−K}e^{Kζ²/2}I_K(ζ) with a Gaussian bound (the exterior Plancherel–Rotach
    asymptotics with their proof);
  - the circle identity Re f(μ) − Re f(μe^{iψ}) = (1 − cos ψ)(1 − |μ|²cos(2 arg μ + ψ)) and the maximum principle give c_K(√Kμ) →
    1 − μ² on the disc;
  - Laplace's method on Duhamel's identity, under a domination condition (proved on both axes up to 2, grid-checked on |η| < 1.8927),
    gives P there and the cumulant limits.
  - Hence D‴.61 (iv)'s sharp asymptotics hold at every frame with r < 0.8053: all 22 of A6's rules with a frame, and D‴.58 (i)'s
    half-squeeze sharp constant for σ*/σ < 9.27.
  - The round-91 referee found the identity and the condition independently.
  - Also: the corner ratio's Dirichlet form, and its continuation on the imaginary axis to y = 1.8637, with the heuristic alias rate
    (11/18)y² − log y − 3/2 beyond.
  - Numerically, a **deformed Duhamel path** connects 0 → η/2 → η inside the bound's sublevel set at every tested η up to |η| = 1.98,
    sector included (scratch notes, not yet in the repository): P on the whole disc looks within reach (§7).
  - Later the same day (rounds 99–121): the deformed path certified to |η| = 1.995 (part V), the corner's four regimes as one path family
    (part VC, Lemma VC), the last annulus certified by round 120's referee (part VA) and verified by round 121's, which found a NaN that
    dropped one connector and fixed it: **conjecture P proved on the whole open disc, computer-assisted** (by hand since round 148); D‴.60's alarms all the cold
    restart's (round 119; round 114's "neither" withdrawn); D‴.59's long-double exit on the stored data (round 116).

## 7. Queue

Referee reports received and not yet applied (each needs: apply → log entry → commit → next referee round):
- **D‴.59, round 87 (D3ch)**, `review/REVIEW_round87_D3ch.md`: D3ch-1 (major) A2's 3A row at σ*/σ₀ = 0.25 is held on the saddle by
  LSODA (other solvers leave it; the plateau then an exact end); minors: A2E's damping (BDF's order-2 steps amplify, order-1 damp; record
  the order), the split criterion's classes (ends that fit none; Σw's sign flips within one structure; the null test's norms), "approaching
  a null cluster" (say "undecided at 10⁹" and print the ratio to the residual), tests (a full-md test; `sq_args` field by field; `sq_job` with
  two amplitudes; GL's floor), cosmetics.
- **D‴.52, round 86 (D3cg)**, `review/REVIEW_round86_D3cg.md`: all minor (G against the whole band, not like for like; the null-as-unit
  tests count runs; the plug-in's loss below the band; tests miss the central lines; the double-bootstrap remedy untested; wording).
- **D‴.54, round 79 (D3bz)**: minors (test gaps; "deciding x ≤ 5.48 in any arithmetic" overreaches; the coordinates' mechanism; cosmetics).
- Older: **D‴.49, round 51 (D3ax)** and **D‴.50, round 52 (D3ay)** (check what later rounds superseded).
- Re-reviews due: D‴.44 (after round 34); D‴.57 (after round 50).

Documentation and the manuscript:
- Move the finished drafts into `PARADIGM.md` / `results/RESULTS_PARADIGM.md` / `README.md` / `PAPERS.md` once their rounds converge
  (pending since long ago: D‴.40, D‴.44, D‴.50, D‴.55).
- The LaTeX manuscript: done (2026-11-06) — D‴.40, D‴.44, D‴.49–D‴.62 and the corrections added, `lyapunov.pdf` rebuilt (541 pages).

Open mathematics (highest value, §8):
- **The corner at η = 2: done, computer-assisted** (part VA, rounds 120–121; see §5), **and by hand since round 148** (the annulus note's
  path Π_η avoids the corner; the hand proof along VA's own path below is no longer needed). What would have added value: an analytic proof of the
  last annulus — the four limits (Ψ_φ, Ψ_∞, Ψ̃_ψ, the fixed-θ limit) each have a proved negative path (Lemma VC), and in (u, q) =
  (ε/X, θ²/X) they are faces of one square, so a hand proof needs a remainder bound uniform on that square (VA certifies it in boxes);
  and D61's Numerics (VC) still lacks rounds 115 and 117's numbers, and VC_MID a point on the cap's far piece (round 120, deferred).
- The corner ratio beyond the disc on the imaginary axis: my threshold 1.8637 and rate were wrong (round 93): the first alias is the
  saddle of I_K's integrand at |ζ| = 1/√2, rate (5/8)y² − log y − ½ − (3/2)log 2, zero 1.8589. Fix part D (`alias_rate`, a Prony fit
  over consecutive K from the exact integer series) and the text.
- The rules outside the frame picture (moment matching, |ρ| ≥ 1): their rates (moment matching's λ is numerical, A3).
- The half-squeeze rule's tail constant ₂F₁(−½, −½; 1; t²) − 1 (D‴.58 (i)) now follows from P; a direct proof would be independent.

## 8. Where the impact is (for choosing what to do next)

Ranked by novelty × impact:
1. **The scaled Gauss–Hermite rules' sharp rates (Paper 18; D‴.58, D‴.61).** One formula t(v) for every scaled rule's worst-case error in
   the Gaussian kernel's space, minimised at the half-squeeze scaling with a cusp, and sharp constants — where the literature (Kuo–
   Woźniakowski; Karvonen–Oates–Girolami) has non-matching bounds. Now proved for every frame (r < 1; by hand since round 148, the
   certificates D, V and VA in ball arithmetic an independent cross-check), through the Peano kernel's transform and the corner ratio of e^{uJ_K}:
   conjecture P holds on the whole open disc. A standalone numerical-analysis paper (SIAM J. Numer. Anal., Math. Comp., FoCM) is ready
   to draft; the open ends were the rules without a frame (|ρ| ≥ 1; D‴.62 since) and a hand proof of the corner (round 148 since).
2. The rich transit (PAPERS.md, Paper 17): an adaptive controller grows a sparse kernel network one feature at a time; D‴.53, D‴.55 proved;
   the fitted skeleton predicts the outcome on all 1154 test targets.
3. Noise selecting among equally good fits (Papers 14–15).
4. Uniform-in-time propagation of chaos with a plant in the loop (Paper 3); adaptive priors removing a leak's bias (Paper 7).

Recent rounds have spent much effort making numerics referee-proof; necessary, but not new science. After the in-flight rounds converge,
the highest-value work is: draft Paper 18 (conjecture P is now proved, by hand since round 148); extend D‴.61 (iv) beyond the frames (|ρ| > 1,
round 113's two-axis transient); consolidate the drafts into PARADIGM.md and the manuscript; draft Papers 14–17 as papers.

## 9. Environment

- Python 3.13 in `venv/` (not in git). The raws were produced with numpy 2.5.4 and scipy 1.18.1 (`tools/requirements-pinned.txt`). Tests
  that rerun a job and compare it bit-identically with the raw may depend on these versions and on the BLAS; if they fail after a
  reinstall, check the versions before suspecting the code.
- `pdflatex` was installed in the previous container (`apt-get install texlive-latex-extra` or similar); a new container may lack it.
- Unreachable hosts: arxiv.org, ams.org, ar5iv.labs.arxiv.org. Web search snippets are reachable.
- Test-suite times at handoff: `tests/test_width_mismatch.py` about 55 s, `tests/test_noise_resonance.py` about 100 s, `tests/test_gh_peano.py` about 17 s.


## 10. State at handoff

- Last commits: eb6e7d3 (D‴.61 round 127: the proof read as a whole, Duhamel proved), 6465011 (D‴.59 round 125), d3db5ac (D‴.58
  round 123), b361d6d (D‴.60 round 124), 31487a6 (D‴.61 round 126), 5f56b63 (D‴.61 round 122).
- Tests (committed): `test_width_mismatch.py` 66, `test_gh_peano.py` 43, `test_noise_resonance.py` 83 (more in the round-128
  application), `test_moving_width.py` 56.
- Every md regenerates byte-identically from its raw.
- Referee ids: D3db = 107 … D3dw = 128 applied; D3dx = 129 (D‴.58) applied; D3dy = 130 (D‴.59) applied; D3dz = 131 (D‴.61)
  applied; D3ea = 132 (D‴.62) applied; D3eb = 133 (D‴.59), D3ec = 134 (D‴.60, C8A closed) and D3ed = 135 (the manuscript) applied; D3ee = 136 (the revised manuscript, minor revision) applied (638fb1b); D3ef = 137 (D‴.59; ulp-scale traps at 3A, part A2T) applied (75469ae); D3eg = 138 (D‴.60 as a whole) applied; D3eh = 139 (the c₂ note) applied; D3ei = 140 (the Paper 17 draft) applied; D3ej = 141 (the Paper 16 draft) applied; D3ek = 142 (the Paper 15 draft) applied; D3el = 143 (the Paper 14 draft) applied; D3em = 144 … D3es = 150 applied (Papers 16, 15, 14, 17 second rounds; D‴.44; D‴.40; the annulus note); D3et = 151 (Paper 15, third), D3eu = 152 (Paper 17, third), D3ev = 153 (Paper 14, third), D3ew = 154 (Paper 18 after the by-hand annulus), D3ex = 155 (acceptance check of Papers 16, 14) applied; next D3ey = 156.
- **D‴.61's state (the main advance of this session):** Q proved; the transform's limit proved on |η| < 1 and on both axes; on |η| ≤
  1.8927 by the branch and bound (rechecked in ball arithmetic, part DB); on 1.8927 ≤ |η| ≤ 1.995 by part V (the deformed path's bound
  in the saddle's coordinates, log|1 − z²| + Re(η²z²/4) + C(σ), certified per η-box, every bound in ball arithmetic); on 1.995 ≤ |η| < 2
  by part VA (round 120's certificate, round 121's verification and fix). Hence (iv) for every frame with r < 1 and D‴.58 (i)'s
  half-squeeze sharp constant for every σ*/σ; all three computer-assisted steps rigorous in ball arithmetic (Arb's, python-flint's and CPython's float arithmetic
  and `math.nextafter` assumed correct). Since round 148 the whole disc is by hand (`review/RESEARCH_note_annulus.md`, Theorem A; D61's
  upgrade note), and D, V and VA are independent cross-checks.
- D‴.59: round 97's "invariant for the problem as posed" was wrong (round 100): Lemma D‴.59a; the held runs are BDF's fixed point.

## State at the close of this stretch (2026-11-05)

Papers 14–18 (`manuscript/paper14` … `paper18`, each with a supplement, README, number check and figures) have been through referee
rounds to a minor-revision or accept verdict, every verdict applied:

| Paper | rounds | last verdict |
|---|---|---|
| 14 (gauge noise: Katzenberger's drift, a Fixman correction, white noise of any covariance) | 143, 150, 153, 155 | accept after editorial fixes (applied) |
| 15 (the colour of label noise as a sparsity prior) | 142, 145, 151 | minor, applied (Lemma 2.6 refereed correct) |
| 16 (the free mass–amplitude law) | 141, 144, 155 | accept after a one-line fix (applied) |
| 17 (the small-mass transit, the sliding limit, an AGF-type correspondence) | 140, 149, 152 | minor, applied |
| 18 (sharp rates of scaled Gauss–Hermite rules; conjecture P by hand on the whole disc) | 135, 136, 154 | minor, applied |

Left for the authors: placeholders (author, affiliation, DOI, code-release URL), the bibliographic `% CHECK` items, venue formatting.
Open research: D‴.40's finite-τ question (PARADIGM §5d, D‴.40's Part N); Paper 17's exact plane certificate beyond the truncation; the open problems listed in
each paper's last section.

## The IJCAI-27 stretch (2026-11-07 and 08): imaging experiments, baselines, astronomy

Work on the paper "Implicit model order selection in RBF networks" (Paper 17's sliding limit applied to localization microscopy and crowded-field astronomy). The paper text is NOT in this
repository (it must stay private until the decision); this section indexes the code and the results that its numbers come from. Every report regenerates byte-identically from its raw
file with the command in the last column (set the output-directory variable to a scratch directory to compare); every script has tests, and the new ones have mutation
specifications under `tools/mutants/` (`wp_*_own.py`; the docstring of each states how to run it and what survives and why).

| experiment | script | report (raw) | regenerate / output variable |
|---|---|---|---|
| separated regime: the path against its decoupled limit (720 targets; part C: 180 targets at 10, 12 and 14 sigma) | `experiments/exp_separated_regime.py` | `results/chaos/separated_regime.md` | `--from-raw` / `SEP_OUT` |
| synthetic order selection against the baselines, order strata, the network against the path (noisy) | `experiments/exp_noise_order.py`, `src/core/baselines.py` | `results/chaos/noise_order.md` | `--from-raw` / `NO_OUT` |
| order of entry on controlled close pairs | `experiments/exp_order_pairs.py` | `results/smlm/order_pairs{,_snr5,_bicc}.md` | `--from-raw` |
| 2016 localization sets: the ISS, the BIC stop, the label-free width, the design variants, against the Picasso ladder | `experiments/exp_smlm_iss.py`, `exp_smlm_picasso.py`, `smlm_common.py`, `src/core/imaging.py` | `results/smlm/iss.md`, `picasso_baseline.md`, `calibration.md` | `--from-raw` / `ISS_OUT` |
| the 2013 simulated sets | `experiments/exp_smlm2013.py`, `smlm2013.py`, `data/smlm/fetch_2013.py` | `results/smlm/smlm2013.md` | `--from-raw` / `SMLM2013_OUT` (`SMLM_DATA_2013=/nonexistent`) |
| a real sequence by Fourier ring correlation | `experiments/exp_smlm_real.py`, `frc.py` | `results/smlm/real_frc.md` | `--from-raw` |
| a CPU learned localizer, in distribution and under mismatch (needs torch for the training stages only) | `experiments/exp_smlm_learned.py`, `src/core/learned_localizer.py` | `results/smlm/learned.md` | `--from-raw` / `LEARNED_OUT` |
| cost of the frame fit in CPU time, two tolerances | `experiments/exp_smlm_speed.py` | `results/smlm/speed.md` | `--from-raw` / `SPEED_OUT` |
| crowded-field astronomy: M2 against HST, the frozen test on four patches, the transfer test on M53 and M92, the design variants | `experiments/exp_astro_m2.py`, `exp_astro_baselines.py`, `astro_common.py`, `astro_second.py`, `astro_patch_survey.py` | `results/astro/m2.md` (`freeze.json`, `m2_raw.npz`, `baseline_catalogs*.npz`) | `--from-raw` / `ASTRO_OUT` |
| the kernel of the path against the PSF stamp (POST HOC, labelled so) | `experiments/exp_astro_psf.py`, `src/core/psf_general.py` | `results/astro/psf.md` | `--from-raw` / `PSF_OUT` |
| the astronomy baselines given the ISS's front end (3 or 5 px high-pass of the frame and the PSF stamp, the ISS's weights) and a wider pipeline search; family declared before scoring (`baselines_hp_plan.md`), selected on the two M2 development patches | `experiments/exp_astro_baselines_hp.py`, `exp_astro_hp.py`, `tools/mutants_py.py` | `results/astro/baselines_hp.md` (`baselines_hp_raw.npz`, `baseline_catalogs_hp.npz`, `baselines_hp_plan.md`) | `--from-raw` / `HP_OUT` |
| the baselines with the ISS's front end, round 2 (declared after referee BHP: bounded photutils fits, crowdsource passes and a 2 px high-pass) | `experiments/exp_astro_baselines_hp2.py`, `exp_astro_hp2.py` | `results/astro/baselines_hp2.md` (`baselines_hp2_plan.md`, `baselines_hp2_raw.npz`, `baseline_catalogs_hp2.npz`) | `--from-raw` / `HP2_OUT` |
| post hoc sensitivity of the held-out repeat test (truth offset, block size, tolerances, depth, scaled lengths, camera column) | `experiments/exp_astro_repeat_sens.py` | `results/astro/repeat_sens.md` (`repeat_sens_raw.npz`) | `--from-raw` / `REPEAT_SENS_OUT` (domain environment) |
| the PSF stamp against stacks of stars, exact-stamp synthetic worlds, splitting and flux bias of bright stars (post hoc, label-using) | `experiments/exp_astro_diag.py` | `results/astro/diag.md` (`diag_raw.npz`) | `--from-raw` / `DIAG_OUT` |
| guarded re-runs of the astronomy paths hit by the enter/exit limit cycle of `imaging.sliding_path` (12 stage-1 paths, check mode on 65 other published paths, corrected development table and corrected M92/repeat rows; post hoc), plus an SMLM spot check | `src/core/imaging_guard.py` (guarded copy; `imaging.py` stays frozen), `experiments/exp_astro_cyclefix.py`, `exp_cycle_scan.py` | `results/astro/cyclefix.md` (`cyclefix_plan.md`, `cyclefix_raw.npz`), `results/smlm/cycle_scan.md` | `--from-raw` / `CYCLEFIX_OUT` |
| the factorized network trained by gradient flow on imaging frames against the path (24 declared frames) | `experiments/exp_network_imaging.py` | `results/synthetic/network_imaging.md` | `--from-raw` / `NI_OUT` |
| held-out repeat imaging of M2: the four test windows in three other SDSS frames (stamp FWHM 2.3, 3.0-3.1, 3.5 px), six patches; analysis declared and frozen before the run (`repeat_plan.md`, `freeze_repeat.json`) | `experiments/exp_astro_repeat.py`, `astro_repeat.py`, `exp_astro_baselines_repeat.py`, `data/astro/fetch_repeat.py` | `results/astro/repeat.md` (`repeat_raw.npz`, `repeat_spec.md`, `baseline_catalogs_repeat.npz`) | `--from-raw` / `REPEAT_OUT` |

Conventions in force: development data choose, a mechanical freeze (`results/astro/freeze.json`, hashes of the code, constants, catalogues and choices) is checked by the astronomy test parts,
test patches and frames are run once; nothing is tuned on a test set; a result that was seen before a design choice is labelled post hoc. Raw third-party data are never in the repository
(fetch scripts with checksums: `data/smlm/`, `data/astro/`). Referee reports are archived verbatim in `review/` (`REVIEW_roundSM1`, `IS1a`, `IS2`, ...).
