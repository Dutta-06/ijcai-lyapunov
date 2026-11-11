# Peak / pulse counting on real 1-D signals: dataset scouting (workstream J, step 1)

Question: is there a public dataset of real 1-D signals in which the NUMBER (ideally also the positions) of peak or pulse components has a
ground truth that is independent of any fit, with at least 100 signals and roughly Gaussian-shaped components, so that the paper's
mechanism (a small-start factorized Gaussian network whose training chooses the number of units) can be tested on its model order?

Scouting date: 2026-11-10. Hosts tried (from this machine, through the allowlist proxy): reachable were Kaggle (datasets through the MCP
tools, files through signed `storage.googleapis.com` links), Zenodo (API and files), GitHub (API and raw), PyPI. Refused by the proxy
(HTTP 403 on CONNECT, not worked around): physionet.org, api.figshare.com, huggingface.co, archive.ics.uci.edu, www.openml.org, rruff.info,
data.mendeley.com, osf.io, www.nist.gov, dataverse.harvard.edu, export.arxiv.org, www.ncbi.nlm.nih.gov. So the classical spectroscopy
repositories (RRUFF, NIST, PhysioNet, UCI, OpenML) were not available directly; only what was mirrored on Kaggle or Zenodo.

## Candidates examined

| # | dataset (host, id) | size, licence | what the ground truth is | independent of a fit? | peak shape | verdict |
|---|---|---|---|---|---|---|
| 1 | **CPSC 2019 training set mirror** (Zenodo 23179391; official challenge data, ICBEB 2019) | 2000 single-lead ECG records, 10 s at 500 Hz, 58 MB; record states CC BY 4.0, the data stay under the challenge terms | R-peak (QRS) positions in every record, annotated by experts; 6 to 38 beats per record (mean 14.7); RR intervals down to 0.2 s (5 % of intervals below 0.32 s) | yes: human annotation of the waveform, no model involved; count is the number of annotated R peaks | QRS is a narrow biphasic spike, not Gaussian; after band-pass, rectification and smoothing (the Pan-Tompkins front end) each beat is a smooth bump of about 0.1 s, roughly Gaussian | **used for step 2** (qualifies on the letter; Gaussianity only after a front end) |
| 2 | MIT-BIH Noise Stress Test (Kaggle elieli94/mit-nst; WFDB files) and MIT-BIH arrhythmia on Kaggle (mondejar/mitbih-database etc.) | 70 MB / 500 MB, ODC-By / CC0 | beat annotations of two clean records (118, 119) re-used under six calibrated noise levels (24 to -6 dB) | yes (cardiologist annotation; noise calibration known) | as in 1 | interesting for a noise-ladder follow-up, but only two underlying records (little diversity), so not chosen |
| 3 | ISP ECG delineation (Zenodo 14679837) | 475 records, 41 MB, CC BY-ND 4.0 | P, QRS and T waves segmented by two cardiologists | yes | P and T are Gaussian-like (the McSharry ECG model is a sum of Gaussians) but the QRS is a sharper multi-lobed complex; the natural fitted count of a QRS is 2 to 3 Gaussians, the label is 1 | rejected: label count and model count are different objects, so a "count error" would measure the convention, not the order selection |
| 4 | RPT15 demonstration subset (Zenodo 21652811) | 5 subjects, 6 MB, ODbL | R, P, T peaks, manual | yes | as 3 | rejected: 5 subjects only |
| 5 | HPGe gamma-ray spectra of certified reference materials (Zenodo 23093533) | 255 spectra, 5 MB, CC BY 4.0 | certified activity concentrations of the 226Ra series, 228Ra, 40K, 137Cs; no peak list | the expected line list could be rebuilt from nuclear data (branching ratios, an efficiency curve) but nuclear data hosts (NNDC, NIST) are blocked and the expected detectable-line count depends on an efficiency model | photopeaks are genuinely Gaussian (best physical fit of all candidates) but counts are low (Poisson), lines overlap with continuum and each other | not usable in this spike: no peak-count label; reconstructing one is a project of its own. Best physics fit if someone builds the line list |
| 6 | Gamma spectra, IAEA inspired (Kaggle marufmasum/gamma-spectra-iaea) | 58 MB, other | isotopic concentrations of U, Pu, MOX samples, 16384 channels | no peak list | Gaussian photopeaks, hundreds of lines | rejected: no peak count label |
| 7 | Liquid-chromatography Raman of four amino acids (Kaggle sergioalejandrod/raman-spectroscopy; Hiramatsu) | 62 MB, CC BY 4.0 | one run, four analytes known (glycine, leucine, phenylalanine, tryptophan) | elution peaks correspond to the known analytes, but it is a single run (a few records) | Gaussian-like elution peaks | rejected: far fewer than 100 independent signals |
| 8 | RRUFF powder XRD (Kaggle brandonpian/rruff-powder; DIF peak lists and calculated profiles) | 234 MB, MIT | peak lists from the crystal structure | independent of any fit, but the profiles are CALCULATED from the structure, not measured; RRUFF itself (experimental patterns) is blocked | pseudo-Voigt | rejected as real data (it would be a synthetic benchmark with realistic peak lists); the experimental XY_RAW archive was not reachable |
| 9 | Peak Detection Dataset (Kaggle milankalkenings/peak-detection-dataset) | 29 MB, CC0 | 30000 length-50 signals with a peak count in `peaks(TARGET)` (2 to 9) | the count is by construction | synthetic (inspected: values in [0, 1], generated) | rejected: synthetic, not real data |
| 10 | Synthetic multi-phase XRD (Kaggle tarikegemenliksiz/...), simulated NaI spectra (sebas02/...), synthetic urban gamma spectra (Zenodo 7987216), chromatogram annotation-study data (GitHub catac0mb/...) | various | by construction | n/a | synthetic or simulated | rejected: not real |
| 11 | Ground-truth recordings for spike sorting (Zenodo 1205233); Spikefinder (GitHub codeneuro) | 36 GB; data hosted elsewhere | spike times from simultaneous intracellular recording | yes | spikes are not Gaussian bumps; calcium fluorescence is an exponential-decay kernel | rejected: shape mismatch (and size) |
| 12 | Cell / microplastic / bacterial Raman sets, FTIR sets, NMR sets, MALDI/CASMI mass-spectrometry sets, MIT-BIH variants (Kaggle) | various | class labels or molecule identities, not peak counts | no | n/a | rejected: no component-count label |
| 13 | Dig4Bio Raman competition (Kaggle) | competition data | concentrations of glucose, sodium acetate, magnesium sulfate | the label is a concentration, not a count of peaks | n/a | rejected: no peak-count target |

The Kaggle search returned only the first page of each query (12 to 20 items) and the Kaggle dataset search is keyword based, so this is a
scouting survey, not an exhaustive catalogue. Several searches (gamma, chromatogram, photoelectron waveform, spike sorting, XPS) returned
nothing relevant.

## Honest assessment of the fit to the mechanism

* No candidate offers real-measurement spectra (Raman, IR, XRD, MS, gamma) with an independent peak-count label. The nearest physics fit
  (gamma photopeaks of certified materials, #5) has no peak list, and the reference archives that would supply one (NIST, RRUFF, NNDC) are
  not reachable from here. This is the main negative finding of the scouting.
* The only dataset that has (a) real signals, (b) an independent human count label and (c) at least 100 records is the ECG family, and within it
  CPSC 2019 is the largest and has the cleanest label (one R peak per beat). Its weaknesses for this paper are stated plainly:
  1. The QRS complex is not a Gaussian; the test needs a fixed front end (band-pass, rectify, smooth) that turns each beat into a bump. The
     front end is a design choice and is applied to the signal before any method sees it. All methods get the same envelope.
  2. Beats are well separated for most records (median RR 0.64 s against a bump width of about 0.1 s), so most records have no overlap and a
     peak finder is close to optimal. The mechanism's claimed strength is resolving overlapping components; only records with fast rhythms
     (RR near 0.25 to 0.35 s) test that, and the noisy or artefact-laden records test robustness. The prior expectation is a tie with a tuned
     peak finder; a win would need to come from the artefact records.
  3. The label is a beat count, which equals the number of bumps only if the front end suppresses P and T waves and noise. Errors of a method
     therefore mix model-order error with front-end error, identically for every method that shares the envelope.
* Baselines available: `scipy.signal.find_peaks` on the same envelope with tuned prominence/distance, an envelope matched filter, adaptive
  threshold (Pan-Tompkins type), forward selection with BIC (`baselines.stepwise_bic`) and matching pursuit with sliding
  (`baselines.omp_slide`) in this repository. `lmfit` is not installed here and not needed (forward selection with BIC is the same model
  family). Published deep models for CPSC 2019 exist but are not run here (no GPU, no time); they are mentioned as a ceiling that this spike
  does not compete with.

Decision: step 2 is run on CPSC 2019 as a narrow pulse-counting test, with the front end and every method tuned on a development split only
and the verdict decided by `results/ml/peakcount_plan.md`. The dataset is borderline for "Gaussian-ish peaks" (the bump is Gaussian-like only
after the front end), and this is stated in the plan and will be restated in the report whatever the verdict.

Data handling: the raw archive is not in git. `experiments/exp_peakcount.py --fetch` re-downloads it from Zenodo (record 23179391, file
`CPSC2019-train.zip`) and checks SHA-256 `1a838bcab24d5b93de5dfebd0c9f5b43b9a621415a6472f2056f6b79c7d5c7fc` (value measured on 2026-11-10; also
the MD5 6bf346bde7962f28f7464dca9c35a252 stated by the record).

## Outcome of step 2 (added after the test run)
The spike on CPSC 2019 was run as pre-registered (`results/ml/peakcount_plan.md`, `results/ml/peakcount.md`): verdict KILL. On 500 held-out records the sliding inverse scale space path
counts R peaks with a mean absolute error of 3.72 beats, the best tuned classical counter (matched filter on the same envelope) with 1.17; the path is statistically tied with greedy sliding matching
pursuit (3.67) and with forward selection (3.44). The scouting caveat 2 above (beats rarely overlap, so a peak finder is near optimal) is consistent with this, but the loss is larger than a tie: the
path families have no refractory (minimum distance) prior and their residual-based stopping rules do not separate beats from T waves and noise bumps.
