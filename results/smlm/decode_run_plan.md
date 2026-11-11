# Plan of the DECODE run on the four 2D challenge sets (written before the final run)

**Written 2026-11-09 at about 09:22 UTC (file saved before the final run was submitted at 09:23:25 UTC; no detection on a test frame existed at either time).** The run is made on a free cloud GPU session (two Tesla T4, 15 GB each; the session's shell has no other GPU) instead of
the authors' machine of `decode_run_spec.md`. Everything below was fixed before any DECODE detection on a TEST frame existed. What had been run before this plan: a pipeline test on the CPU of a
cloud session, and two short GPU timing runs (3-4 epochs on MT0.N1.LD, 860 s and 640 s of GPU time) that inferred and scored the DEVELOPMENT frames only (best development E after 3 epochs: 57.4; after 2
epochs: 54.0) and no test frame; they fixed the speed figures that set the budget below, nothing else was tuned on them.

## What is run

* **DECODE**, source of commit `4c0f38c33681d1c3005bc18bd99a11f46d30f1c1` (`decode.__version__ = 0.10.1dev1`, the clone the spec was written from; the notebook downloads that commit's archive from
  GitHub and checks the sha256 of five of its files against the clone's before installing it with `pip --no-deps`). Environment: conda (micromamba 2.10) with python 3.8.18, **pytorch 1.7.1
  (cuda 11.0.221, cudnn 8.0.5)**, numpy 1.19.5, scipy 1.8.0, **spline 0.10.0 (the turagalab build with CUDA kernels)**, mkl 2024.0 (pytorch 1.7.1 does not import with a newer mkl), the rest as
  `environment.yaml` of the commit; these are the versions of the commit's own `environment.yaml` (the code does `from torch._six import container_abcs`, which torch >= 1.9 removed: the first environment
  with torch 1.12 failed at import, nothing was patched). The exact package versions are returned by the run (`versions.txt`) and quoted in `decode.md`.
* **Where the data come from.** The notebook (private, internet on) downloads the eight public files (four stacks, four `activations.csv`) from `srm.epfl.ch` and checks their sha256 against
  `data/smlm/checksums.json`. Nothing of this repository, no result, no ISS code, no title or author information is on the notebook. The ground-truth file contains every frame; the worker keeps
  the rows of the development frames the moment it has read the file (`load_gt_dev`) and never uses the others.
* **One model per set** (MT0.N1.LD, MT0.N1.HD, MT0.N2.LD, MT0.N2.HD), the first training run of each is the only one: **one configuration per set** (the ISS tried 8 model configurations x 7 thresholds x 3
  post-filters and Picasso ~10^2 configurations on the same development frames; the cost of ours is stated next to the result). Seed 0 (python, numpy, torch, cuda), no restart unless DECODE's own
  convergence check (`auto_restart_param`) triggers one.

## What DECODE is given (all from the information that the ISS and Picasso also receive: the camera constants of the dataset pages, and the development frames and their ground truth)

| item | value | source |
|---|---|---|
| camera | baseline 100 DN, e/ADU 45, EM gain 300, read sigma 74.4 e, spurious 0.002, **QE 0.90**, pixel 100 nm, 64 x 64 | dataset pages (as the spec) |
| background (photons) | N1 99.93, N2 9.91, simulated uniformly in +-10 % | mean of the quiet development frames (no emitter in the frame or its neighbours) of MT0.N1.LD / MT0.N2.LD; the HD sets use their noise class's value (they have no quiet frames) |
| photons per emitter | normal(mean, sd) of the set's development ground truth, clamped at its 1 % quantile (N1.LD 4693 / 3097 / 140; N1.HD 4658 / 3116 / 163; N2.LD 507 / 335 / 15; N2.HD 503 / 337 / 18); emitters below the same photon count are not labelled (`emitter_label_photon_min`) | development ground truth of the set |
| axial positions | uniform on [0, 600] nm (symmetric PSF, so the sign of z is not identifiable: 0..z_hi); 600 = the 97 % quantile of the development |z| rounded up to 50 nm | development ground truth |
| emitters per frame (simulated) | LD 3 (denser than the data's 1.25, so that overlaps are seen), HD 14 (DECODE draws ~0.9 of the nominal number: ~12.7 realised against the data's 12.9); lifetime 2.3 frames | development ground truth (density, mean on-time from the whole-set runs in `calibration.md`) |
| **PSF** | a **cubic-spline calibration built from the radial single-emitter profile**: for the single-emitter development frames of MT0.N1.LD (frames 0-999) the pixel value above the quiet-frame background per ground-truth photon as a function of the distance of the pixel centre to the ground-truth position, in 5 bins of |z| (86 / 144 / 75 / 23 / 32 emitters in 0-100 / 100-200 / 200-350 / 350-500 / 500-800 nm); smoothed beyond 6 px, PCHIP in r, linear in |z|; zero beyond 26 px. 51 x 51 px ROI, 41 z cells of 50 nm. The profile integrates to 1.42 photons per ground-truth photon (halo included: `calibration.md` section (b)), so DECODE's photon count is in ground-truth units. Option 1 of the spec (empirical z-dependent PSF from the development frames), with a radially symmetric model: the asymmetry of the real PSF (the +11 / -9 nm offset of a symmetric estimator) is not represented and is removed afterwards by the development-estimated constant offset. No bead stack, no SMAP: the .mat file is written by the notebook in SMAP's layout (`SXY.cspline.{coeff, x0, z0, dz}`) and its coefficients checked against the compiled spline library (single emitters rendered with DECODE's own `CubicSplinePSF` agree with the profile to a correlation of 0.994, brightest pixel at the right row and column). One PSF for all four sets (the profile does not depend on the noise class: `calibration.md`), measured on the development frames of MT0.N1.LD. | development frames and ground truth of MT0.N1.LD |

Known handicaps, to be quoted next to the numbers: (i) the PSF is a symmetric radial model measured on 360 single-emitter frames, not a bead calibration; (ii) the simulated axial prior is DECODE's
uniform one, the data's is concentrated near the focus (median |z| 143 nm); (iii) the training budget is a small fraction of the authors' (below); (iv) the HD sets reuse the PSF and the
background of the LD sets.

## Training

* DECODE's defaults (`reference.yaml`) except: `Camera.qe` 0.9, `Simulation` as above, `max_number_targets` 40 (LD) / 60 (HD) instead of 250 (DECODE raises an error if a simulated frame has more emitters;
  an iteration took ~1.0 s per batch of 64 on a T4 with 250 padded targets; the loss tensors scale with this number, how much of the second it accounts for was not measured), `TestSet.test_size` 256, `Hardware.num_worker_train` 2 per process, `HyperParameter.epochs`
  a large number (training is stopped by the clock, below), `learning_rate_scheduler_param.gamma` 1 (DECODE's step schedule switched off). Network `SigmaMUNet` (48 features, default), loss
  `GaussianMMLoss`, AdamW (lr 2e-4, weight decay 0.1), batch 64, 10,000 freshly simulated frames per epoch (`pseudo_ds_size`), 3 input channels (frame window).
* **Budget.** The sets are trained two at a time, one per GPU, in two waves (MT0.N1.LD with MT0.N1.HD, then MT0.N2.LD with MT0.N2.HD). Wall-clock budget per set: at most 5,000 s of training (the
  driver reduces it if the session bound of 10,800 s would be exceeded); the epoch loop stops when another epoch would overrun the budget. The learning rate decays exponentially by a factor 10 over
  the budget (lr(t) = 2e-4 x 0.1^(t / budget), set after every epoch), in place of DECODE's step schedule. The GPU time of the whole session is bounded by 10,800 s (the account's GPU allowance
  for the week has ~57,000 s left; the whole job is held to 14,400 s, of which 1,500 s are used by the tests above). DECODE's authors train for ~1000 epochs of this size (hours to days); this run gets
  at most a few tens to ~150 epochs depending on the speed reached: the number of epochs actually run is reported per set.

## Selection (development frames only)

1. Three snapshots of the weights per set: the first epoch end at or after 1/3 and at or after 2/3 of the budget, and the last epoch.
2. Each snapshot is run with DECODE's own inference pipeline (`SigmaMUNet`, camera back-conversion, non-maximum suppression and probability aggregation of `SpatialIntegration`) with the
   detection threshold lowered to 0.05 on the 1000 (LD) / 100 (HD) development frames; every candidate (frame, x, y, photons, probability) is kept.
3. **Axis convention**: DECODE's first coordinate runs along the image rows (checked by rendering an emitter at a known coordinate: brightest pixel at that row), the ground truth's x along the columns.
   The notebook converts both ways (x_nm = (y_px + 0.5) 100 and the plain one), scores each on the development frames and keeps the one with the higher best development E (expected: the swap; the
   alternative is recorded).
4. **Snapshot**: the one with the highest development E (offset-corrected, best of the 12 thresholds 0.05..0.98 of `smlm_decode_ingest.THRESHOLDS`; ties to the later snapshot), E computed with a copy of
   the harness's matching (Hungarian, 250 nm, pooled J and RMSE). `experiments/exp_smlm_decode.py` recomputes the development scores with the project's own `smlm_decode_ingest.evaluate` and records
   whether it agrees.
5. Only then are the TEST frames (LD 1000-1499, HD 100-499) run, once, with the selected snapshot and convention, and the candidates saved.
6. The probability threshold (twice: by E and by J) and the constant position offset are chosen on the development frames by `smlm_decode_ingest.evaluate` and applied to the test candidates. DECODE-E
   (E-chosen) is the primary row, compared with the E-tuned ISS (offset removed) and Picasso R2E+c, which use the development labels in the same way.

Metrics, matching and the paired block bootstrap (blocks of 10 frames, 1000 resamples, seed 0) are those of `experiments/smlm_common.py`, as for the other rows; the comparison rows are the stored ones of
`learned_raw.npz` (ISS a priori default, E-tuned ISS, Picasso R0 and R2E+c, and the CPU-scale CNN of `learned.md`).

## What will be reported whatever the outcome

The numbers of all four sets with their intervals, including losses; the number of epochs and seconds per set; the development scores of the three snapshots; the axis convention found; every deviation
from this plan (written into `decode.md`), and the GPU seconds used. The notebook is private; the Kaggle notebooks created for the work are listed in the hand-over message so they can be deleted.

## Appendix, added after the run (2026-11-09, 12:22 UTC)

The run followed the plan; no step of the protocol was changed. The notebook text that was submitted at 09:23:25 UTC is `results/smlm/decode_gpu_notebook.txt` (sha256 `9a2baad8a7b62adea10412f0d458c54fa9f2f725c090d4ece583f2ee5c2c91b6`); it runs only inside such a session
(it installs the environment, fetches the public data, builds the calibration, trains, selects and returns the result files through its log). The number of epochs reached, 36-37 per set (4,942-4,994 s of training each), is an outcome
of the speed, not a choice: 0.79-0.85 s per iteration of batch 64 on a T4, only 15-20 % less than the ~1.0 s measured with 250 padded targets, so the loss over padded targets was not the main cost (the plan's "a few tens to ~150 epochs"
bounded it correctly). The returned files (candidates of the selected snapshots on the development and the test frames, the selection records, the training summaries, versions) are what `decode.md` and `decode_raw.npz` hold;
the weights (8.8 MB per set) and the notebook's logs stayed on the cloud session. The account's GPU counter went from 104,287 s to 116,830 s over the whole job (12,543 s: 860 s and 640 s for the two timing runs, 11,043 s for this run).
