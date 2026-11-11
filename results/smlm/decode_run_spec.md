# Running DECODE on the 2D single-molecule benchmark sets — what to run on the GPU machine and what to return

Written for the authors; the session that wrote it has no GPU and has not run DECODE. Everything below that is a *measured fact about the data* comes from the
development frames of each set (`calibration.md`, `iss.md`); everything about DECODE comes from reading its repository at the pinned commit (the clone's
`decode/utils/reference_files/*.yaml`, `param_friendly.yaml` and the two command-line entry points), not from running it. **Check each DECODE key against your installed
version before the first training run.**

## 1. What is wanted

For each of the four sets `MT0.N1.LD`, `MT0.N1.HD`, `MT0.N2.LD`, `MT0.N2.HD` (the 2D microtubule training sets of the 2016 localisation challenge, public ground truth),
a DECODE model trained for that set's imaging conditions, applied to the development and the test frames, returning *all candidates with their detection probability* (not
only those above the default threshold). The session then chooses the probability threshold on the development frames and scores the test frames with the same
harness, matching, offset correction and bootstrap as Picasso and the sliding inverse scale space (`experiments/smlm_decode_ingest.py`). Nothing about the test frames may
enter any choice of yours (training data, hyper-parameters, PSF, stopping): only the development frames and their ground truth (exported) may.

## 2. The frames

On any machine that has this repository and the data (`python data/smlm/fetch.py --sets MT0.N1.LD MT0.N1.HD MT0.N2.LD MT0.N2.HD`, ~230 MB, sha256 in
`data/smlm/checksums.json`):

```bash
python experiments/smlm_decode_export.py --out /path/to/decode_run
```

writes `<out>/<set>/dev.tif`, `test.tif` (16-bit TIFF stacks, DN, baseline 100 DN as in the challenge), `frames.json` (indices, camera) and `dev_gt.csv` (the ground truth of the
development frames only: `frame,x_nm,y_nm,z_nm,photons`, frame 0-based within `dev.tif`; positions in nm in the **corner convention**, the centre of pixel j at (j + 0.5)·100 nm, x
along the columns). Protocol: LD sets dev = stack frames 0–999, test = 1000–1499; HD sets dev = 0–99, test = 100–499 (500 and 400 test frames).

## 3. Software

DECODE at commit `4c0f38c33681d1c3005bc18bd99a11f46d30f1c1` (the version of the repository clone examined; `git log -1` there is 2023-06-22), its `environment.yaml`, a CUDA GPU.
Record in your return: the DECODE commit, torch and CUDA versions, the GPU model, the random seeds, the wall-clock training time per set, and the number of training
runs/configurations tried per set (the other methods were tuned on the same development frames with comparable effort — Picasso over ~10² configurations, the inverse scale
space over 8 model configurations × 7 thresholds × 3 post-filters; say what yours was).

## 4. Parameters that are known (from the challenge's simulator and the development frames)

| item | value | source |
|---|---|---|
| frame size | 64 × 64 px, pixel 100 nm | dataset page |
| `Camera.baseline` / `e_per_adu` / `em_gain` / `read_sigma` / `spur_noise` | 100 / 45 / 300 / 74.4 / 0.002 | dataset page; the same numbers are the example in DECODE's `param_friendly.yaml` ("values of the camera used in the challenge") |
| `Camera.qe` | 0.90 (DECODE's default is 1.0: set it) | dataset page ("photon converter factor 0.90 e⁻/photon") |
| `Camera.px_size` | [100.0, 100.0] | |
| background (photons per pixel) | N1: 99.9 (sd 15 per pixel, flat); N2: 9.5–9.9 | `calibration.md`, quiet frames |
| noise | var = 2.23 x + 0.1 (photon units; excess-noise factor 2/QE) | `calibration.md` |
| emitters per frame | LD 1.25; HD 12.9 (development frames) | |
| photons per emitter | N1: mean 4690, sd 3100, median 4260, 5/95 % = 390 / 10100; N2: mean 507, sd 335, median 460, 5/95 % = 42 / 1090 | development frames' ground truth |
| z range of the emitters | −750 … 750 nm (sd 240 nm) | ground truth |
| on-time | an emitter stays on for 2.3 frames on average | whole-set ground truth |
| matching tolerance for the evaluation | 250 nm lateral (`Evaluation.dist_lat: 250`) | the challenge |

Suggested starting values for the training simulation (yours to adjust on the *development* frames): `Simulation.img_size: [64, 64]`, `emitter_extent` x, y = `[[-0.5, 63.5], [-0.5, 63.5]]`
(DECODE's convention: pixel centres at the integers), `bg_uniform` = ±10 % of the set's background, `intensity_mu_sig` = the set's (mean, sd) above, `emitter_av` = the set's density
(for LD train with a few emitters per frame at least, so that the network sees overlaps; the test density is 1.25), `lifetime_avg: 2.3`, `Hardware.device: cuda:0`.
`Evaluation.match_dims: 2` (the lateral metrics only).

## 5. The one decision that is not ours: the PSF

DECODE trains on simulated frames and needs the PSF of the microscope as a cubic-spline calibration (`InOut.calibration_file`), normally measured from a bead stack. **No bead stack
is available for the 2D modality of these sets.** What the data say about the PSF (`calibration.md`, referee report `review/REVIEW_roundSM1.md`):

* it is not the optics' Gaussian of 0.93 px: the core is ~1.3 px wide in focus, broadens with |z| (1.28 px at |z| < 100 nm, ~2.4 px at 350–700 nm: the 2D modality images
  the same z-range, with a defocus broadening and no astigmatism), and it is not Gaussian — it has a halo: about 44 % of the emitted photons arrive beyond 5 px from the centre;
* the centre of a symmetric Gaussian fit is displaced from the simulated position by an estimator- and z-dependent offset of about (+11, −9) nm in focus.

Options, in the order we would try them (record which one you used; the dev frames decide between them if you have time):
1. **An empirical z-dependent PSF built from the development frames and their ground truth** (isolated bright emitters, aligned on their known position, averaged by |z| bin and
   oversampled), turned into a bead-stack-like TIFF and calibrated with DECODE's SMAP workflow. We can provide the stack builder if you want to go this way (it is not written yet:
   say so and it will be, from `dev.tif`/`dev_gt.csv`).
2. **A synthetic bead stack** from the measured radial profile by |z| (core width and halo fraction per |z| bin, in `calibration.md`), fed to the same calibration step.
3. Any calibration you already have for a comparable symmetric-defocus PSF; state its origin.
Because the PSF is a modelling input that DECODE needs and the sliding inverse scale space does not need (it fits the width of each emitter), a mismatch here is a *DECODE-side* handicap
that must be stated next to its numbers.

## 6. Train, infer, return

```bash
# training: edit a copy of param_friendly.yaml (section 4), then
python -m decode.neuralfitter.train.train -p param_<set>.yaml -l runs -c <set>
# inference on the exported TIFFs: one fit meta file per set and split (the template is decode/utils/reference_files/fit.yaml):
#   Frames.path: dev.tif | test.tif;  Model.path / param_path: your trained model; Output.path: <set>_dev.csv | <set>_test.csv
python -m decode.neuralfitter.inference.infer -p fit_<set>_dev.yaml
```

Set the post-processing threshold of the model used for inference low (`PostProcessingParam.raw_th: 0.05`, default 0.5) so that low-probability candidates are kept; do not apply
any other filter. Convert DECODE's output (`frame_ix` 0-based within the TIFF, `x`, `y` in px, `phot`, `prob`) to **exactly** these columns, one CSV per set and split:

```
frame,x_nm,y_nm,photons,prob          frame = frame_ix;  x_nm = (x_px + 0.5) * 100;  y_nm = (y_px + 0.5) * 100;  photons = phot;  prob = prob
```

(DECODE's pixel centres are at the integers; the +0.5 px converts to the corner convention of the ground truth. A slip here shows up as a mean signed error of ±50 nm on the
development frames, which the ingest script prints as a warning.) Return: `<set>_dev.csv`, `<set>_test.csv` for the four sets, a text file with the items of section 3, and the
trained models if you can. Then, anywhere with this repository:

```bash
python experiments/smlm_decode_ingest.py --set MT0.N1.LD --dev MT0.N1.LD_dev.csv --test MT0.N1.LD_test.csv --save
```

which chooses the probability threshold on the development frames twice (by the efficiency E and by the Jaccard index J), removes the development-estimated constant offset, scores the test
frames, and stores the detections in `results/smlm/decode_raw.npz` for the paired comparison with Picasso and the sliding inverse scale space.

## 7. What to say next to the numbers

* DECODE's own leaderboard results on the challenge's *test* datasets (hidden ground truth, 3D and 2D) are published by its authors; they are **not comparable** with numbers on the MT0 training sets
  and cannot be used here.
* The MT0 sets are simulated and share one structure (three microtubules): they are not independent datasets. A network trained on simulated emitters scattered over the field is not
  trained on the microtubule structure, but state it.
* The sliding inverse scale space is a model-based, single-frame method with no training and a per-frame cost of ~1 s (LD) to ~5 s (HD) on one CPU core; DECODE trains for hours on a GPU and then
  fits a frame in milliseconds. Report both.
