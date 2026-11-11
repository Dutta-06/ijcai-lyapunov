"""Fast tests of src/core/learned_localizer.py: the target encoding and the decoder on synthetic maps with known answers (numpy only), the symmetry
transforms of frames and maps, the input normalisation, and, when torch is installed, a one-batch overfit, determinism of the seeded training, the weight
round trip, test-time augmentation on an exactly equivariant stand-in, and a tiny end-to-end example (train on simulated frames with 1-3 emitters, recover them
within 100 nm). The torch tests are skipped without torch (the repository's own venv has none; the learned baseline runs in a separate environment, see
`experiments/exp_smlm_learned.py`). Whole file: well under a minute on one thread (one thread on purpose: on an oversubscribed machine two OpenMP threads can be 30 times slower on kernels this small). LEARNED_LONG_TEST=1 makes the end-to-end test train ten times longer."""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from src.core import learned_localizer as LL      # noqa: E402

needs_torch = pytest.mark.skipif(LL.torch is None, reason='torch is not installed')


def _perfect_maps(emitters, shape=(64, 64)):
    """The maps a perfect network would output for the emitters: the target probability, offsets and log photons of the voting neighbourhood (zero elsewhere)."""
    t = LL.encode_targets([emitters], shape)
    return t['prob'][0].astype(float), t['off'][0] * t['mask'][0], t['logph'][0] * t['mask'][0]


# ------------------------------------------------------------------ encoding
def test_encode_marks_the_containing_pixel_and_the_voting_offsets():
    t = LL.encode_targets([np.array([[10.3, 5.8, 4000.0]])], (16, 16))
    assert t['n_emit'][0] == 1
    assert t['prob'].sum() == 1 and t['prob'][0, 6, 10] == 1                       # row round(5.8) = 6, column round(10.3) = 10
    assert np.allclose(t['off'][0, :, 6, 10], [0.3, -0.2], atol=1e-6)                # (x - j, y - i)
    assert np.allclose(t['off'][0, :, 7, 11], [10.3 - 11, 5.8 - 7], atol=1e-6)        # a neighbour votes for the same emitter
    assert t['mask'].sum() == 9 and t['mask'][0, 5:8, 9:12].all()
    assert np.allclose(t['logph'][0, 6, 10], np.log(4.0), atol=1e-6)


def test_encode_skips_emitters_outside_the_field_and_clips_neighbourhoods_at_the_border():
    t = LL.encode_targets([np.array([[-0.8, 3.0, 1000.0], [16.2, 3.0, 1000.0]])], (16, 16))      # round(-0.8) = -1 and round(16.2) = 16: outside a 16-pixel field
    assert t['n_emit'][0] == 0 and t['prob'].sum() == 0
    t = LL.encode_targets([np.array([[15.2, 3.0, 1000.0]])], (16, 16))                          # round(15.2) = 15: the last column
    assert t['n_emit'][0] == 1 and t['prob'][0, 3, 15] == 1
    t = LL.encode_targets([np.array([[0.1, 0.2, 500.0]])], (16, 16))
    assert t['prob'][0, 0, 0] == 1 and t['mask'][0].sum() == 4                                  # corner: a 2 x 2 neighbourhood


def test_encode_overlapping_neighbourhoods_go_to_the_nearer_emitter():
    t = LL.encode_targets([np.array([[8.0, 8.0, 1000.0], [10.0, 8.0, 3000.0]])], (16, 16))
    assert np.allclose(t['off'][0, :, 8, 9], [-1.0, 0.0])                            # pixel 9 is 1 px from each emitter: a tie, kept by the first
    assert np.allclose(t['off'][0, :, 8, 8], [0.0, 0.0]) and np.allclose(t['off'][0, :, 8, 10], [0.0, 0.0])


# ------------------------------------------------------------------ decoder
def test_decoder_recovers_the_encoded_positions_exactly():
    em = np.array([[10.3, 5.8, 4000.0], [40.6, 50.1, 800.0], [0.2, 63.0, 2500.0]])
    p, o, l = _perfect_maps(em)
    d = LL.decode_frame(p, o, l, threshold=0.5)
    assert d.shape == (3, 4)
    got = d[np.lexsort((d[:, 1], d[:, 0]))]
    want = em[np.lexsort((em[:, 1], em[:, 0]))]
    assert np.allclose(got[:, 0], (want[:, 0] + 0.5) * 100.0, atol=1e-6) and np.allclose(got[:, 1], (want[:, 1] + 0.5) * 100.0, atol=1e-6)
    assert np.allclose(got[:, 2], want[:, 2], rtol=1e-5) and np.allclose(got[:, 3], 1.0)


def test_decoder_pools_probability_split_over_neighbouring_pixels():
    em = np.array([[20.5, 30.0, 2000.0]])                                               # on the boundary of pixels 20 and 21
    p, o, l = _perfect_maps(em, (64, 64))
    p[30, 20], p[30, 21] = 0.6, 0.4                                                     # the network is unsure which pixel
    d = LL.decode_frame(p, o, l, threshold=0.5, window=1)
    assert len(d) == 1 and abs(d[0, 3] - 1.0) < 1e-12                                   # the window sum is the detection probability
    assert abs(d[0, 0] - (20.5 + 0.5) * 100.0) < 1e-6 and abs(d[0, 1] - 3050.0) < 1e-6   # the votes agree on the true position
    d0 = LL.decode_frame(p, o, l, threshold=0.5)                                        # the default window 0: the single pixel scores only 0.6
    assert len(d0) == 1 and abs(d0[0, 3] - 0.6) < 1e-12
    dp = LL.decode_frame(p, o, l, threshold=0.5, window=1, pos_window=1)               # pooling the votes of the window gives the same position here (the votes agree)
    assert abs(dp[0, 0] - d[0, 0]) < 1e-9 and abs(dp[0, 1] - d[0, 1]) < 1e-9


def test_decoder_threshold_equals_a_lower_threshold_filtered_afterwards():
    rng = np.random.default_rng(0)
    p = np.zeros((64, 64))
    for k in range(12):                                                                  # blobs of varied strength
        i, j = rng.integers(2, 62, 2)
        p[i - 1:i + 2, j - 1:j + 2] += rng.uniform(0.02, 0.3) * np.array([[0.1, 0.1, 0.1], [0.1, 0.6, 0.1], [0.1, 0.1, 0.1]])
        p[i, j] += rng.uniform(0.0, 0.5)
    off = rng.uniform(-0.5, 0.5, (2, 64, 64))
    lph = rng.normal(0, 1, (64, 64))
    low = LL.decode_frame(p, off, lph, threshold=0.05)
    for thr in (0.1, 0.3, 0.6):
        hi = LL.decode_frame(p, off, lph, threshold=thr)
        assert np.array_equal(hi, low[low[:, 3] >= thr])
    assert np.all(np.diff(low[:, 3]) <= 0)                                               # sorted by decreasing score


def test_decoder_nms_and_empty_maps():
    assert LL.decode_frame(np.zeros((64, 64)), np.zeros((2, 64, 64)), np.zeros((64, 64)), 0.5).shape == (0, 4)
    p = np.zeros((16, 16))
    off = np.zeros((2, 16, 16))
    p[8, 5], p[8, 7] = 0.9, 0.6                                                          # two local maxima 2 px apart
    lph = np.zeros((16, 16))
    assert len(LL.decode_frame(p, off, lph, 0.5)) == 2                                   # kept: 2 px >= 1.5 px
    off[0, 8, 7] = -0.9                                                                  # the weaker one votes for a position 1.1 px from the stronger
    d = LL.decode_frame(p, off, lph, 0.5)
    assert len(d) == 1 and abs(d[0, 3] - 0.9) < 1e-12                                    # suppressed by the stronger candidate within 1.5 px
    p = np.zeros((16, 16))
    p[8, 5] = p[8, 6] = 1.0                                                              # an exact tie between neighbours: one detection
    assert len(LL.decode_frame(p, off * 0, lph, 0.5)) == 1


# ------------------------------------------------------------------ symmetries, normalisation
def test_d4_augmentation_and_map_transforms_are_consistent():
    rng = np.random.default_rng(1)
    em = np.c_[rng.uniform(0, 31, 3), rng.uniform(0, 31, 3), rng.uniform(500, 5000, 3)]
    f = rng.normal(size=(32, 32))
    base = LL.encode_targets([em], (32, 32))
    for k in range(8):
        fk, ek = LL.augment_d4(f, em, k)
        assert np.allclose(fk, LL.d4_maps(f, k))                                          # the frame transform is the map transform
        tk = LL.encode_targets([ek], (32, 32))
        assert np.array_equal(LL.d4_maps_inverse(tk['prob'][0], k), base['prob'][0])      # the targets of the augmented emitters are the transformed targets
        assert np.array_equal(LL.d4_maps_inverse(tk['mask'][0], k), base['mask'][0])
        assert np.allclose(LL.d4_offsets_inverse(tk['off'][0], k), base['off'][0], atol=1e-5)   # offset vectors included
        assert np.array_equal(LL.d4_maps_inverse(LL.d4_maps(f, k), k), f)


def test_normalize_removes_surface_and_median_and_scales():
    rng = np.random.default_rng(2)
    bg = 100.0 + np.linspace(0, 2, 64)[None, :] * np.ones((64, 1))
    fr = bg + 7.0 + 4.0 * rng.standard_normal((5, 64, 64))
    x = LL.normalize(fr, bg, 4.0)
    assert x.dtype == np.float32 and np.allclose(np.median(x.reshape(5, -1), axis=1), 0, atol=1e-6) and abs(x.std() - 1.0) < 0.03
    assert np.allclose(LL.normalize(fr[0], bg, 4.0), x[0])                               # a single frame is accepted
    assert abs(LL.noise_scale(100.0, 2.3) - np.sqrt(2.3 * 100.0 + LL.B0)) < 1e-12


# ------------------------------------------------------------------ simulation
def test_simulator_density_units_and_labels():
    rng = np.random.default_rng(3)
    fr, ems = LL.simulate_frames(300, rng, bg=100.0, density=2.0, photons=np.array([3000.0]), sigma=1.4, shape=(32, 32))
    assert fr.shape == (300, 32, 32) and abs(np.mean([len(e) for e in ems]) - 2.0) < 0.25
    assert all(np.all((e[:, 0] >= -0.5) & (e[:, 0] < 31.5)) and np.all(e[:, 2] == 3000.0) for e in ems if len(e))
    fr2, ems2 = LL.simulate_frames(200, np.random.default_rng(4), bg=100.0, density=1.0, photons=np.array([3000.0]), sigma=1.4, shape=(32, 32), photon_scale=0.5,
                                   psf_scale=1.25)
    assert all(np.all(e[:, 2] == 1500.0) for e in ems2 if len(e))                         # the scaled photon counts are the labels
    empty = [f for f, e in zip(fr2, ems2) if len(e) == 0]
    assert abs(np.mean(empty) - 100.0) < 0.2 and abs(np.std(empty) - np.sqrt(2 / 0.9 * 100 + 0.078)) < 1.0      # the camera's noise law at a flat background


# ------------------------------------------------------------------ torch
def _tiny_batch(rng, n=8, size=24, em_max=3):
    from src.core import imaging as IM
    psf = IM.PixelPSF(size, size, 1.4)
    frames, ems = np.empty((n, size, size)), []
    for i in range(n):
        k = int(rng.integers(1, em_max + 1))
        c = rng.uniform(2, size - 3, (k, 2))
        ph = rng.uniform(2500, 6000, k)
        frames[i] = rng.poisson(100.0 + psf.synth(c, ph)).astype(float)
        ems.append(np.c_[c, ph])
    return LL.normalize(frames, 100.0, 10.0), ems


@needs_torch
def test_one_batch_overfit_decreases_the_loss():
    LL.torch.set_num_threads(1)
    x, ems = _tiny_batch(np.random.default_rng(5), 4)
    tgt = LL.encode_targets(ems, (24, 24))
    model = LL.make_model(0, c1=8, c2=16, c3=16)
    _, hist = LL.train(model, lambda step: (x, tgt), 120, lr=3e-3, warmup=10)
    first, last = np.mean([h[1] for h in hist[:5]]), np.mean([h[1] for h in hist[-5:]])
    assert last < 0.4 * first and all(np.isfinite(h[1]) for h in hist)


@needs_torch
def test_training_is_deterministic_given_the_seed():
    LL.torch.set_num_threads(1)
    x, ems = _tiny_batch(np.random.default_rng(6), 4)
    tgt = LL.encode_targets(ems, (24, 24))
    ws = []
    for _ in range(2):
        m = LL.make_model(7, c1=8, c2=16, c3=16)
        LL.train(m, lambda step: (x, tgt), 6, lr=3e-3, warmup=2)
        ws.append([p.detach().numpy().copy() for p in m.parameters()])
    assert all(np.array_equal(a, b) for a, b in zip(*ws))
    other = LL.make_model(8, c1=8, c2=16, c3=16)
    assert not all(np.array_equal(a, p.detach().numpy()) for a, p in zip(ws[0], other.parameters()))        # a different seed gives a different network


@needs_torch
def test_save_load_roundtrip_and_parameter_count(tmp_path):
    m = LL.make_model(1)
    assert LL.n_params(m) == 289040 and LL.n_params(m) <= 300000                        # the default architecture (c = 32, 64, 64) is within the budget
    m = LL.make_model(1, c1=8, c2=16, c3=16)
    x = _tiny_batch(np.random.default_rng(7), 3)[0]
    path = str(tmp_path / 'w.npz')
    LL.save(m, path, dict(seed=1, note='x'))
    m2, meta = LL.load(path)
    assert meta == dict(seed=1, note='x') and m2.arch == m.arch
    for a, b in zip(LL.predict(m, x), LL.predict(m2, x)):
        assert np.array_equal(a, b)


@needs_torch
def test_test_time_augmentation_on_an_equivariant_stand_in(monkeypatch):
    rng = np.random.default_rng(8)
    x = rng.normal(size=(3, 24, 24)).astype(np.float32)

    def fake(model, xx, batch):                                                           # maps built from the frame by symmetric operations only
        xx = xx.astype(np.float64)
        gy, gx = np.gradient(xx, axis=(1, 2))
        return 1.0 / (1.0 + np.exp(-xx)), np.stack([gx, gy], axis=1), xx ** 2
    monkeypatch.setattr(LL, '_predict_maps', fake)
    m = LL.make_model(0, c1=8, c2=16, c3=16)
    a, b = LL.predict(m, x, tta=False), LL.predict(m, x, tta=True)
    for u, v in zip(a, b):
        assert np.allclose(u, v, atol=1e-12)                                              # equivariant maps: averaging over the 8 orientations changes nothing


@needs_torch
def test_end_to_end_synthetic_example_recovers_emitters_within_100_nm():
    LL.torch.set_num_threads(1)
    long = os.environ.get('LEARNED_LONG_TEST') == '1'
    steps = 5000 if long else 500
    model = LL.make_model(3, c1=12, c2=24, c3=24)

    def batch_fn(step):
        x, ems = _tiny_batch(np.random.default_rng([11, step]), 12)
        return x, LL.encode_targets(ems, (24, 24))
    LL.train(model, batch_fn, steps, lr=3e-3, warmup=20)
    x, ems = _tiny_batch(np.random.default_rng(99), 40)
    dets = LL.detect(model, x, threshold=0.3, window=1)                              # (a tiny network spreads an emitter's probability over neighbouring pixels: the window sum collects it)
    tp = fp = fn = 0
    errs = []
    for d, e in zip(dets, ems):
        truth = (e[:, :2] + 0.5) * 100.0
        used = set()
        for r in d:
            dist = np.hypot(truth[:, 0] - r[0], truth[:, 1] - r[1])
            j = int(np.argmin(dist))
            if dist[j] <= 100.0 and j not in used:
                used.add(j)
                tp += 1
                errs.append(dist[j])
            else:
                fp += 1
        fn += len(truth) - len(used)
    assert tp / (tp + fn) >= 0.88 and tp / max(tp + fp, 1) >= 0.9, (tp, fp, fn)
    assert np.median(errs) < 40.0


# ------------------------------------------------------------------ the experiment script's own logic (numpy only)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'experiments'))
import exp_smlm_learned as E      # noqa: E402
import exp_smlm_iss as X          # noqa: E402
import smlm_common as C           # noqa: E402


def test_the_simulator_uses_the_iss_psf_model_and_its_width_bounds():
    import inspect
    default_parts = inspect.signature(LL.simulate_frames).parameters['parts'].default
    assert tuple(default_parts) == tuple(X.CONFIGS['hp5h']['parts'])                      # the halo model of the a priori default of the ISS
    w = E.width_sampler()(np.random.default_rng(0), 20000)
    assert 1.0 <= w.min() and w.max() <= 3.2 and 1.30 < np.median(w) < 1.45 and 0.0 < np.mean(w > 2.0) < 0.2


def test_mismatch_tags_and_the_isss_wrong_inputs():
    assert E.tag_cfg('sim') == E.SimCfg(1.0, 1.0, 1.0) and E.tag_cfg('sup') == E.SimCfg(1.0, 1.0, 1.0) and E.tag_cfg('sup200') == E.SimCfg(1.0, 1.0, 1.0)
    assert E.tag_cfg('psf0.8') == E.SimCfg(0.8, 1.0, 1.0) and E.tag_cfg('bg2') == E.SimCfg(1.0, 2.0, 1.0) and E.tag_cfg('ph0.5') == E.SimCfg(1.0, 1.0, 0.5)
    with pytest.raises(ValueError):
        E.tag_cfg('nonsense')
    cal = dict(bg=np.full((4, 4), 100.0), a=2.3)
    c, o = E.iss_inputs(cal, 'psf1.25')
    assert abs(o['sigma0'] - 1.75) < 1e-12 and np.all(c['bg'] == 100.0) and c['a'] == 2.3 and o['rules'] == [('snr', 5.0)]
    c, o = E.iss_inputs(cal, 'bg0.5')
    assert np.all(c['bg'] == 50.0) and c['a'] == 2.3 and 'sigma0' not in o                  # the camera law applied to the believed background
    c, o = E.iss_inputs(cal, 'bg2/a')
    assert np.all(c['bg'] == 200.0) and abs(c['a'] - 1.15) < 1e-12                            # the slope recomputed for the wrong surface: the noise at a background pixel stays right
    assert cal['a'] == 2.3 and np.all(cal['bg'] == 100.0)                                     # the caller's calibration is not modified
    assert len(E.model_keys()) == 8 + 10 and E.model_keys()[0] == ('MT0.N1.LD', 'sup')


def test_labels_of_simulated_emitters_round_trip_through_the_harness_conventions():
    ems = [np.array([[10.0, 20.0, 1000.0], [30.5, 5.25, 2000.0]]), np.zeros((0, 3)), np.array([[0.0, 63.0, 500.0]])]
    gt = E.gt_from_emitters(ems)
    assert list(gt['frame']) == [0, 0, 2] and abs(gt['x_nm'][0] - 1050.0) < 1e-9 and abs(gt['y_nm'][2] - 6350.0) < 1e-9        # nm = (px + 0.5) * 100: the corner convention
    dets = [np.c_[(e[:, 0] + 0.5) * 100.0, (e[:, 1] + 0.5) * 100.0, e[:, 2], np.ones(len(e))] if len(e) else np.zeros((0, 4)) for e in ems]
    r = C.score(dets, gt, np.arange(3))
    assert r['jaccard'] == 1.0 and r['rmse_nm'] < 1e-9
    back = E.gt_emitters(np.sort(gt, order=['frame', 'id']), 0, 3)
    assert np.allclose(back[0], ems[0]) and len(back[1]) == 0 and np.allclose(back[2], ems[2])


def test_worker_claims_are_exclusive_reentrant_and_stale_claims_are_taken_over(tmp_path):
    d = str(tmp_path)
    assert E.claim(d) and E.claim(d)                                      # ours, and again by the same process
    os.remove(os.path.join(d, 'claim.pid'))
    open(os.path.join(d, 'claim.pid'), 'w').write(str(os.getppid()))      # a live process that is not us
    assert not E.claim(d)
    open(os.path.join(d, 'claim.pid'), 'w').write('999999999')            # no such process: stale
    assert E.claim(d) and open(os.path.join(d, 'claim.pid')).read() == str(os.getpid())
