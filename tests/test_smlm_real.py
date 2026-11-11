"""Tests of the pure functions of experiments/exp_smlm_real.py (the real-data experiment): the loader of per-frame TIFF zips, the label-free calibration
(noise law from a synthetic noisy 128 x 128 movie with a planted law, the PSF-width rule), the ISS pipeline on noisy frames with emitters of known position,
the chunk cache and the worker pool, the repeat-localization precision, the densest window, the agreement, the packing and the report (deterministic).
No downloads and no Picasso needed. A few seconds."""
import io
import json
import os
import sys
import zipfile

import numpy as np
import pytest
from PIL import Image

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'experiments'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import exp_smlm_real as R      # noqa: E402
import frc as FRC              # noqa: E402
from src.core import imaging as IM  # noqa: E402

NY = NX = 128
A_TRUE, V0_TRUE, X0 = 0.8, 120.0, 200.0         # planted noise law: var = A (x - X0) + V0 for the pixel mean x


def _bg_surface():
    yy, xx = np.mgrid[0:NY, 0:NX]
    return X0 + 60.0 + 30.0 * xx / NX + 20.0 * yy / NY


HALO = ((0.75, 1.0), (0.25, 3.0))                  # the kernel family of the ISS configuration 'hp5h' (core + a halo of 3x the width)


def synth_frame(rng, centres, photons, sigma=1.5, bg=None, parts=((1.0, 1.0),)):
    """A noisy frame: background surface + emitters (pixel-integrated PSF of core width sigma, ``parts``, integrated counts ``photons``), Gaussian noise with the planted
    heteroscedastic variance A (mean - X0) + V0 (the temporal behaviour of a camera: independent pixels)."""
    bg = _bg_surface() if bg is None else bg
    P = IM.PixelPSF(NY, NX, sigma, parts=parts)
    mean = bg + (P.synth(np.atleast_2d(centres), np.asarray(photons, float)) if len(photons) else 0.0)
    return mean + np.sqrt(A_TRUE * (mean - X0) + V0_TRUE) * rng.standard_normal((NY, NX)), mean


def synth_movie(n, seed=0, n_emit=12):
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n):
        c = rng.uniform(8, 120, (n_emit, 2))
        out.append(synth_frame(rng, c, rng.uniform(8000, 30000, n_emit))[0])
    return np.stack(out)


# ------------------------------------------------------------------ the loader
def _write_zip(path, frames, names):
    with zipfile.ZipFile(path, 'w') as z:
        z.writestr(zipfile.ZipInfo('sequence/'), b'')
        for name, f in zip(names, frames):
            buf = io.BytesIO()
            Image.fromarray(np.asarray(f, dtype=np.uint16)).save(buf, format='TIFF')
            z.writestr('sequence/' + name, buf.getvalue())


def test_loader_reads_per_frame_tiffs_in_numeric_order_and_refuses_gaps_and_overruns(tmp_path):
    rng = np.random.default_rng(0)
    frames = rng.integers(0, 4000, (5, 12, 16))
    path = tmp_path / 'seq.zip'
    order = [3, 1, 4, 2, 5]                                       # stored out of order
    _write_zip(path, [frames[k - 1] for k in order], ['%05d.tif' % k for k in order])
    got = R.load_real_frames(str(path))
    assert got.shape == (5, 12, 16) and got.dtype == np.float64 and np.array_equal(got, frames)
    assert np.array_equal(R.load_real_frames(str(path), 1, 3), frames[1:3])
    assert R.load_real_frames(str(path), 2, 2).shape[0] == 0
    with pytest.raises(IndexError):
        R.load_real_frames(str(path), 0, 6)
    with pytest.raises(ValueError):
        R.load_real_frames(str(path), 3, 2)
    gap = tmp_path / 'gap.zip'
    _write_zip(gap, [frames[0], frames[1], frames[3]], ['00001.tif', '00002.tif', '00004.tif'])
    with pytest.raises(ValueError):
        R.load_real_frames(str(gap))
    empty = tmp_path / 'empty.zip'
    with zipfile.ZipFile(empty, 'w') as z:
        z.writestr('sequence/', b'')
    with pytest.raises(ValueError):
        R.load_real_frames(str(empty))


# ------------------------------------------------------------------ the noise law
def test_noise_law_recovers_a_planted_heteroscedastic_law_from_noisy_frames():
    frames = synth_movie(60, seed=1)
    law = R.calibrate_noise(frames)
    assert law['scale'] > 0 and law['n_frames'] == 12                                     # every 5th of 60 frames
    for x in (260.0, 400.0, 700.0, 1000.0):                                               # the planted variance at the mean level x
        planted = A_TRUE * (x - X0) + V0_TRUE
        est = law['a'] * (x - law['x_ref']) + law['b0']
        assert abs(est / planted - 1) < 0.15, (x, est, planted)
    assert abs(law['scale'] - 1) < 0.25                                                   # white noise: the temporal anchor agrees with the single-frame law
    assert abs(law['c'] - (law['b0'] - law['a'] * law['x_ref'])) < 1e-9


def test_noise_law_is_blind_to_a_pedestal_and_scales_with_the_noise():
    frames = synth_movie(60, seed=2)
    base = R.noise_law(frames)
    shifted = R.noise_law(frames + 500.0)                                                 # a constant added to every pixel moves the abscissa, not the law
    assert abs(shifted['a'] / base['a'] - 1) < 1e-6 and abs(shifted['x_ref'] - base['x_ref'] - 500.0) < 1e-6 and abs(shifted['b0'] / base['b0'] - 1) < 1e-6
    rng = np.random.default_rng(3)
    noisier = frames + 10.0 * rng.standard_normal(frames.shape)                           # + 100 ADU^2 of white noise
    more = R.noise_law(noisier)
    assert 80 < more['b0'] - base['b0'] < 120
    with pytest.raises(ValueError):
        R.noise_law(np.zeros((4, 12, 12)))


def test_temporal_dark_variance_ignores_a_blinking_minority_and_a_per_frame_offset():
    rng = np.random.default_rng(4)
    n = 200
    frames = 300.0 + 10.0 * rng.standard_normal((n, 16, 16)) + 25.0 * rng.standard_normal(n)[:, None, None]       # white 100 ADU^2 + a per-frame offset
    frames[:, :4, :4] += 800.0 * (rng.random((n, 4, 4)) < 0.3)                              # blinking pixels, far brighter than the dark ones
    td = R.temporal_dark_variance(frames, 0.5)
    assert abs(td['v'] / 100.0 - 1) < 0.2 and td['n'] == 128 and 295 < td['x'] < 305
    by = R.temporal_variance_by_level(frames, [250, 350, 700, 2000])
    assert by[0][1] < 150 and all(len(row) == 3 for row in by)


def test_scale_by_segment_is_the_scale_of_the_rule_on_consecutive_segments():
    frames = synth_movie(120, seed=9)
    rows = R.scale_by_segment(frames, 60)
    assert rows.shape == (2, 5) and list(rows[:, 0]) == [0, 60] and list(rows[:, 1]) == [60, 120]
    assert np.all(np.abs(rows[:, 2] - 1) < 0.3) and np.all(rows[:, 3] > 0)                       # white noise: the scale is about 1 on every segment
    assert abs(rows[0, 2] - R.calibrate_noise(frames[:60])['scale']) < 1e-12
    assert R.scale_by_segment(frames[:50], 60).shape == (0, 5)


def test_calibrate_noise_takes_the_scale_from_the_temporal_variance_only_with_enough_frames():
    frames = synth_movie(20, seed=5)
    assert R.calibrate_noise(frames)['scale'] == 1.0                                       # fewer than TD_MIN_FRAMES: the single-frame law as it is
    frames = synth_movie(60, seed=5)
    noisy = frames + np.sqrt(150.0) * np.random.default_rng(6).standard_normal(frames.shape)
    assert R.calibrate_noise(noisy)['scale'] != 1.0


# ------------------------------------------------------------------ sigma0
def _dets(rows):
    """Per-frame (m, 6) arrays (x_nm, y_nm, photons, sx_nm, sy_nm, net gradient) from lists of (x_px, y_px, photons, sx_px, sy_px)."""
    return [np.array([[x * 100, y * 100, ph, sx * 100, sy * 100, 9000.0] for x, y, ph, sx, sy in r]).reshape(-1, 6) for r in rows]


def test_sigma0_rule_uses_bright_isolated_localizations_only():
    rng = np.random.default_rng(0)
    rows = []
    for f in range(30):
        r = [(10 + 20 * k, 10, rng.uniform(6000, 7000), 1.7 + 0.1 * rng.random(), 1.8 + 0.2 * rng.random()) for k in range(5)]      # five isolated bright ones, sqrt(sx sy) ~ 1.8
        r += [(10, 100, 5000.0, 1.0, 1.0), (13, 100, 5200.0, 1.0, 1.0)]                                                              # a bright pair within 8 px: not isolated, narrow
        r += [(100, 60, 100.0, 3.0, 3.0)]                                                                                             # isolated but dim: below the median photon count
        rows.append(r)
    sig = R.sigma0_from_picasso(_dets(rows))
    # brute-force recomputation of the rule
    allr = [(f, x, y, ph, sx, sy) for f, r in enumerate(rows) for x, y, ph, sx, sy in r]
    thr = np.median([a[3] for a in allr])
    sel = []
    for f, x, y, ph, sx, sy in allr:
        iso = all(np.hypot(x - x2, y - y2) > 8.0 for f2, x2, y2, *_ in allr if f2 == f and (x2, y2) != (x, y))
        if iso and ph >= thr:
            sel.append(np.sqrt(sx * sy))
    assert not sig['fallback'] and sig['n_sel'] == len(sel) > 100 and abs(sig['sigma0'] - np.median(sel)) < 1e-9 and 1.7 < sig['sigma0'] < 2.0
    assert sig['n_all'] == 30 * 8 and sig['n_iso'] == 30 * 6 and abs(sig['photon_threshold'] - thr) < 1e-6
    assert abs(sig['median_all'] - np.median([np.sqrt(a[4] * a[5]) for a in allr])) < 1e-9
    few = R.sigma0_from_picasso(_dets([[(10, 10, 5000.0, 2.2, 2.2)]]))                                  # too few qualifying localizations: the simulated sets' width, flagged
    assert few['fallback'] and few['sigma0'] == R.ISS.SIGMA0
    assert R.sigma0_from_picasso([np.zeros((0, 6))])['fallback']
    nan = _dets([[(10, 10, 5000.0, 1.5, 1.5), (50, 10, np.nan, 1.5, 1.5)]] * 60)
    assert R.sigma0_from_picasso(nan)['n_all'] == 60                                                    # a non-finite fit is not used


def test_iss_config_is_the_simulated_sets_configuration_at_their_width_and_scales_with_sigma0():
    c = R.iss_config(R.ISS.SIGMA0)
    assert c['widths'] == pytest.approx(list(R.ISS.WIDTHS)) and c['r_ex'] == pytest.approx(R.ISS.R_EX) and c['kappa'] == 5.0 and c['cfg'] == 'hp5h'
    assert c['hp'] == 5.0 and [tuple(p) for p in c['parts']] == [(0.75, 1.0), (0.25, 3.0)] and c['lat_h'] == 0.5
    d = R.iss_config(1.4 * 1.2)
    assert d['widths'] == pytest.approx([1.2 * w for w in R.ISS.WIDTHS]) and d['r_ex'] == pytest.approx(1.2 * R.ISS.R_EX) and d['sigma0'] == pytest.approx(1.68)
    assert d['hp'] == c['hp'] and d['kappa'] == c['kappa']


# ------------------------------------------------------------------ the ISS on noisy frames
def _calibration():
    return dict(bg=_bg_surface(), a=A_TRUE, b0=V0_TRUE, x_ref=X0)


def test_pipeline_finds_planted_emitters_in_noisy_128x128_frames_and_nothing_in_an_empty_one():
    rng = np.random.default_rng(7)
    centres = np.array([[30.3, 40.6], [64.8, 70.2], [100.4, 25.7]])
    photons = np.array([9000.0, 7000.0, 11000.0])
    cal = _calibration()
    cfg = R.iss_config(1.5)
    ctx = R.make_context(cal['bg'], cal['a'], cal['b0'], cal['x_ref'], cfg, (NY, NX))
    for _ in range(2):
        frame, _ = synth_frame(rng, centres, photons, parts=HALO)
        out = R.process_frame(frame, ctx)
        det = out['det']
        assert det.shape == (3, 5) and not out['hit_cap'] and out['secs'] > 0 and out['n_path'] >= 3
        d = np.linalg.norm(det[:, None, :2] - centres[None, :, :], axis=2)
        assert d.min(axis=0).max() < 0.35 and d.min(axis=1).max() < 0.35                  # every emitter has a unit within 0.35 px (35 nm) and vice versa
        match = d.argmin(axis=0)                                                            # the detection nearest to each emitter
        assert np.allclose(det[match, 3], photons, rtol=0.15)                              # fitted integrated counts
        assert (det[:, 2] >= cfg['widths'][0] - 1e-9).all() and (det[:, 2] <= cfg['widths'][1] + 1e-9).all()
        assert (det[:, 4] > 20).all() and 0.5 < out['chi2'] < 2.0
        assert abs(out['delta']) < 30
    empty, _ = synth_frame(rng, np.zeros((0, 2)), [], parts=HALO)
    assert len(R.process_frame(empty, ctx)['det']) == 0


def test_run_iss_gives_the_same_detections_with_one_process_two_processes_and_the_cache(tmp_path):
    rng = np.random.default_rng(8)
    frames = np.stack([synth_frame(rng, np.array([[20.0 + 10.0 * k, 30.4], [90.3, 100.0 - 10.0 * k]]), [9000.0, 7000.0], parts=HALO)[0] for k in range(3)])
    cal, cfg = _calibration(), R.iss_config(1.5)
    logs = []
    one = R.run_iss(frames, cal, cfg, workers=1, log=logs.append)
    two = R.run_iss(frames, cal, cfg, workers=2, cache_dir=str(tmp_path), chunk=2, log=logs.append)
    again = R.run_iss(frames, cal, cfg, workers=1, cache_dir=str(tmp_path), chunk=2, log=logs.append)
    assert any('read from the cache' in s for s in logs) and len(list(tmp_path.glob('chunk_*.npz'))) == 2
    for a, b, c in zip(one, two, again):
        assert np.allclose(a['det'], b['det'], atol=1e-9) and np.allclose(a['det'], c['det'], atol=1e-9)
        assert a['n_path'] == b['n_path'] == c['n_path'] and len(a['det']) == 2
    # a different configuration or different frames is a different cache entry
    cfg2 = dict(cfg, kappa=6.0)
    assert R._chunk_signature(frames, (cal['bg'], [1.0]), cfg) != R._chunk_signature(frames, (cal['bg'], [1.0]), cfg2)
    assert R._chunk_signature(frames, (cal['bg'], [1.0]), cfg) != R._chunk_signature(frames + 1.0, (cal['bg'], [1.0]), cfg)


def test_bench_reports_wall_and_cpu_seconds_per_frame():
    rng = np.random.default_rng(13)
    frames = np.stack([synth_frame(rng, np.array([[20.3, 24.6]]), [9000.0], parts=HALO)[0][:48, :48] for _ in range(2)])
    cal = dict(bg=_bg_surface()[:48, :48], a=A_TRUE, b0=V0_TRUE, x_ref=X0)
    rows = R.bench_iss(frames, cal, R.iss_config(1.5), np.array([0, 1]))
    assert rows.shape == (2, 4) and list(rows[:, 0]) == [0, 1] and (rows[:, 1] > 0).all() and (rows[:, 2] > 0).all() and (rows[:, 3] == 1).all()
    assert (rows[:, 2] <= rows[:, 1] * 1.5 + 0.05).all()                                       # one process: the CPU time cannot much exceed the wall time


def test_packing_round_trips_and_converts_to_nm_in_the_corner_convention():
    res = [dict(det=np.array([[10.0, 20.0, 1.5, 1000.0, 30.0], [0.0, 0.0, 1.2, 500.0, 25.0]]), n_path=3, hit_cap=False, secs=1.0, delta=0.5, chi2=1.1),
           dict(det=np.zeros((0, 5)), n_path=0, hit_cap=True, secs=0.2, delta=-1.0, chi2=1.0),
           dict(det=np.array([[5.0, 6.0, 1.4, 800.0, 40.0]]), n_path=1, hit_cap=False, secs=0.5, delta=0.0, chi2=0.9)]
    off, flat, st = R._pack_results(res)
    back = R._unpack_results(off, flat, st)
    assert list(off) == [0, 2, 2, 3] and st.shape == (3, 5) and back[1]['hit_cap'] and back[1]['det'].shape == (0, 5)
    assert np.array_equal(back[0]['det'], res[0]['det']) and back[2]['chi2'] == 0.9
    o2, f2 = R.iss_table(res, pixel_nm=100.0)
    assert list(o2) == [0, 2, 2, 3] and np.allclose(f2[0], [1050.0, 2050.0, 1000.0, 1.5, 30.0]) and np.allclose(f2[1][:2], [50.0, 50.0])      # nm = (p + 0.5) * pixel
    assert list(R.frames_of(o2)) == [0, 0, 2]
    dets = [np.array([[1050.0, 2050.0, 100.0, 140.0, 170.0, 9000.0]]), np.zeros((0, 6))]
    o3, f3 = R.pic_table(dets)
    assert f3.shape == (1, 6) and f3[0, 2] == pytest.approx(100.0 * R.C.DN_PER_PHOTON) and f3[0, 3] == pytest.approx(1.4) and f3[0, 4] == pytest.approx(1.7)


# ------------------------------------------------------------------ the precision proxy
def _repeat_data(sigma=18.0, n_frames=400, per_frame=25, seed=0, field=6400.0):
    """Molecules that stay on for 2 consecutive frames at random positions of the field, localized with independent per-axis errors sigma."""
    rng = np.random.default_rng(seed)
    xy, fr = [], []
    for t in range(n_frames - 1):
        pos = rng.uniform(0, field, (per_frame, 2))
        for k in (t, t + 1):
            xy.append(pos + rng.normal(0, sigma, pos.shape))
            fr.append(np.full(per_frame, k))
    return np.concatenate(xy), np.concatenate(fr).astype(np.int64), n_frames


def test_repeat_precision_recovers_the_planted_localization_precision():
    xy, fr, nf = _repeat_data(18.0)
    rp = R.repeat_precision(xy, fr, nf, n_boot=0)
    assert abs(rp['fit']['sigma_nm'] / 18.0 - 1) < 0.08 and 0.8 < rp['fit']['w'] <= 1.0 and rp['n_pairs'] > 3000
    xy2, fr2, _ = _repeat_data(30.0, seed=1)
    assert R.repeat_precision(xy2, fr2, nf)['fit']['sigma_nm'] > rp['fit']['sigma_nm'] * 1.4      # monotone in the planted error
    rb = R.repeat_precision(xy, fr, nf, n_boot=8, seed=0, block=10)
    assert len(rb['boot']) == 8 and np.isfinite(rb['boot']).all() and np.all(np.abs(rb['boot'] / 18.0 - 1) < 0.2)
    rb2 = R.repeat_precision(xy, fr, nf, n_boot=8, seed=0, block=10)
    assert np.array_equal(rb['boot'], rb2['boot'])                                                # a function of the seed alone


def test_next_frame_distances_and_the_mixture_fit_on_known_inputs():
    xy = np.array([[0.0, 0.0], [100.0, 0.0], [0.0, 0.0], [3.0, 4.0], [500.0, 500.0]])
    fr = np.array([0, 0, 1, 1, 3])
    d = R.next_frame_distances(xy, fr, 5)
    assert d[0] == pytest.approx(0.0) and d[1] == pytest.approx(np.hypot(97.0, 4.0)) and np.isnan(d[2]) and np.isnan(d[3]) and np.isnan(d[4])      # frame 1 -> frame 2 is empty
    f = R.fit_repeat_mixture(np.array([10.0, 20.0]))
    assert np.isnan(f['sigma_nm'])                                                                # too few distances
    rng = np.random.default_rng(1)
    pure = rng.rayleigh(30.0, 20000)
    pure = pure[pure <= 150.0]
    fit = R.fit_repeat_mixture(pure)
    assert abs(fit['s_nm'] / 30.0 - 1) < 0.03 and fit['w'] > 0.9 and fit['sigma_nm'] == pytest.approx(fit['s_nm'] / np.sqrt(2))


# ------------------------------------------------------------------ the window, the agreement
def test_densest_window_finds_a_planted_cluster_and_breaks_ties_by_row_major_order():
    rng = np.random.default_rng(9)
    uniform = rng.uniform(0, 12800, (3000, 2))
    cluster = np.array([9000.0, 4000.0]) + rng.normal(0, 300.0, (2000, 2))
    x0, y0, x1, y1 = R.densest_window([uniform, cluster], 3000.0, 12800.0)
    assert x1 - x0 == 3000 and y1 - y0 == 3000 and x0 % 100 == 0 and y0 % 100 == 0
    assert x0 < 9000 - 600 and x1 > 9000 + 600 and y0 < 4000 - 600 and y1 > 4000 + 600
    assert R.densest_window([np.zeros((0, 2))], 3000.0, 12800.0) == (0.0, 0.0, 3000.0, 3000.0)       # no localizations at all: the first window
    with pytest.raises(ValueError):
        R.densest_window([uniform], 20000.0, 12800.0)


def test_agreement_counts_the_same_frame_partners_within_the_radius():
    a = np.array([[0.0, 0.0], [1000.0, 0.0], [0.0, 1000.0]])
    fa = np.array([0, 0, 1])
    b = np.array([[100.0, 0.0], [1000.0, 400.0], [0.0, 1050.0]])
    fb = np.array([0, 0, 0])                                                                       # the partner of a[2] is in another frame
    frac, n = R.agreement(a, fa, b, fb, radius=150.0)
    assert n == 3 and frac == pytest.approx(1.0 / 3.0)
    assert R.agreement(a, fa, b, fb, radius=500.0)[0] == pytest.approx(2.0 / 3.0)


# ------------------------------------------------------------------ the report
def _tiny_raw(nf=40, seed=0):
    """A raw dictionary of the form ``run`` writes, on synthetic localizations (two noisy versions of a filament structure), for the report."""
    rng = np.random.default_rng(seed)
    field = 12800.0

    def sites(n):
        t = rng.uniform(0, 1, n)
        return np.c_[1500 + 9000 * t, 2000 + 3000 * np.sin(3 * t) + 4000 * t]

    def table(n_per, noise, ncol, extra):
        per = []
        for f in range(nf):
            s = sites(n_per)
            cols = [s[:, 0] + rng.normal(0, noise, n_per), s[:, 1] + rng.normal(0, noise, n_per), rng.uniform(3000, 20000, n_per)] + extra(n_per)
            per.append(np.c_[tuple(cols)])
        return R.pack_rows(per, ncol)

    raw = {'README': np.array('x'), 'meta|set': np.array('Real_High_Density'), 'meta|shape': np.array([nf, 128, 128]), 'meta|pixel_nm': np.array(100.0),
           'meta|zip_sha256': np.array('0' * 64),
           'meta|versions': np.array('{"python": "3", "numpy": "2", "picasso": "0.11.3", "cpu_count": 4, "loadavg_at_end": [1.0, 1.0, 1.0], "workers": 2}')}
    raw['cal|bg'] = _bg_surface()
    raw['cal|mean_image'] = np.zeros((128, 128), np.float32)
    raw['cal|frame_delta'] = np.linspace(20, -5, nf)
    raw['cal|law'] = np.array([0.74, -45.0, 243.0, 134.0, 100, 1440000, 1.23, 0.60, 109.0, 140.0, 251.0, 892])
    raw['cal|law_x'], raw['cal|law_v'], raw['cal|law_n'] = np.array([243.0, 300.0, 500.0]), np.array([116.0, 150.0, 265.0]), np.array([100, 200, 300])
    raw['cal|temporal_levels'] = np.array([[250.0, 140.0, 686.0], [300.0, 1558.0, 2282.0]])
    raw['cal|sigma0'] = np.array('{"sigma0": 1.63, "n_sel": 7690, "n_all": 42099, "n_iso": 16387, "n_bright": 21050, "photon_threshold": 1606.0, "median_sx": 1.42, '
                                  '"median_sy": 1.76, "median_all": 1.44, "fallback": false}')
    cfg = R.iss_config(1.63)
    raw['iss|config'] = np.array(json.dumps(cfg))
    raw['iss|off'], raw['iss|flat'] = table(30, 25.0, 5, lambda n: [rng.uniform(1.17, 2.0, n), rng.uniform(10, 200, n)])
    raw['iss|stats'] = np.c_[np.full(nf, 40.0), np.zeros(nf), np.full(nf, 5.0), np.zeros(nf), np.full(nf, 3.0)]
    for key, n_per, noise in (('pic', 20, 30.0), ('pic_ng1000', 24, 32.0)):
        raw[key + '|off'], raw[key + '|flat'] = table(n_per, noise, 6, lambda n: [rng.uniform(1.3, 1.6, n), rng.uniform(1.4, 1.9, n), rng.uniform(6000, 20000, n)])
        raw[key + '|run_s'] = np.array([1.0, float(nf)])
        raw[key + '|config'] = np.array('{"min_net_gradient": 5000.0}')
    win = R.densest_window([raw['iss|flat'][:, :2], raw['pic|flat'][:, :2]], R.ZOOM_NM, field)
    raw['zoom|window_nm'] = np.array(win)
    for key in ('iss', 'pic'):
        flat = raw[key + '|flat']
        m = (flat[:, 0] >= win[0]) & (flat[:, 0] < win[2]) & (flat[:, 1] >= win[1]) & (flat[:, 1] < win[3])
        raw['zoom|' + key] = np.c_[flat[m, 0], flat[m, 1], flat[m, 2], R.frames_of(raw[key + '|off'])[m]]
    return raw


def test_report_is_a_deterministic_function_of_the_raw_and_says_what_it_cannot_measure(tmp_path, monkeypatch):
    monkeypatch.setattr(R, 'FRC_PIXEL_NM', 40.0)                      # a coarse grid keeps the test fast
    raw = _tiny_raw()
    raw['cal|scale_segments'] = np.array([[0.0, 100.0, 1.8, 250.0, 255.0], [100.0, 200.0, 1.3, 180.0, 252.0]])
    raw['iss|bench'] = np.array([[0.0, 2.0, 1.5, 30.0], [20.0, 3.0, 2.5, 31.0]])
    raw['meta|bench'] = np.array('{"n": 2, "loadavg_before": [1.0, 1.0, 1.0], "loadavg_after": [2.0, 2.0, 2.0], "cpu_count": 4}')
    t1 = R.build_md(raw, n_boot=3, n_splits=2)
    t2 = R.build_md(raw, n_boot=3, n_splits=2)
    assert t1 == t2
    assert 'Stability of that scale' in t1 and 'Timing benchmark' in t1 and '2.00 s CPU per frame' in t1 and 'recentred interval' in t1 and 'Paired comparison' in t1
    assert 'does not' in t1 and 'detection accuracy' in t1 and 'ground truth' in t1.lower() and 'ASSUMED' in t1 and '100 nm' in t1
    for needle in ('Picasso, default pipeline (R0)', 'sliding ISS, a priori configuration', 'variation:', 'alternating blocks', 'bootstrap', 'sigma_rep', 'Seconds per frame', 'sigma0 = 1.630'):
        assert needle in t1, needle
    path = tmp_path / 'raw.npz'
    R.C.save_npz_deterministic(str(path), raw)
    md1, md2 = tmp_path / 'a.md', tmp_path / 'b.md'
    x1 = R.write_report(str(path), str(md1), 3, 2, log=lambda s: None)
    x2 = R.write_report(str(path), str(md2), 3, 2, log=lambda s: None)
    assert md1.read_bytes() == md2.read_bytes() and x1 == x2                                         # byte for byte from the raw alone
    assert str(R.file_sha256(str(path)))[:16] in x1


def test_get_set_filters_non_finite_positions_and_the_width_variation():
    raw = _tiny_raw(nf=6)
    off, flat = raw['iss|off'], raw['iss|flat'].copy()
    flat[0, 0] = np.nan
    raw['iss|flat'] = flat
    xy, fr, fl = R.get_set(raw, 'iss')
    assert len(xy) == len(flat) - 1 and len(fr) == len(xy) and np.isfinite(xy).all()
    lo, hi = json.loads(str(raw['iss|config']))['widths']
    xy2, fr2, fl2 = R.get_set(raw, 'iss_interior')
    assert (fl2[:, 3] > lo + R.WIDTH_TOL).all() and len(xy2) <= len(xy)
    pxy, pfr, pfl = R.get_set(raw, 'pic')
    assert pfl.shape[1] == 6 and len(pfr) == len(pxy)
