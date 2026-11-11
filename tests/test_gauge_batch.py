"""
src/core/gauge_batch.py (referee round 143): the batched drift is ReplicatorCohesionLoop.drift to rounding (with per-run γ, κ_F, τ and
floor, the masses and plug-in forms, a box), a batch of runs reproduces run_dithered for the same seeds, a run stopped at a block boundary
and resumed from its checkpoint is the uninterrupted run, a raised cap continues a finished run, and an exit is the first step outside.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'experiments')))
import exp_cohesion_leak as E                                            # noqa: E402
from src.core import mass_transfer as MT                                 # noqa: E402
from src.core.gauge_batch import BatchDither, run_batch                  # noqa: E402

TAU, DT = 1e-4, 0.25


def _configs(B, seed=0):
    net, _ = E.setup()
    rng = np.random.default_rng(seed)
    Th = np.array([E.exact_rep(net, np.ones(3) / 3) + 0.8 * rng.standard_normal((3, 3)) for _ in range(B)])
    return Th, 3 * rng.dirichlet(np.full(3, 20.0), B)


def test_the_batched_drift_is_the_loops():
    net, base = E.setup()
    Th, beta = _configs(6)
    pars = [(0.0, 0.0, 0.01, None, None), (1e-2, 2.0, 0.01, E.A_BUMP, None), (0.1, 1.0, 0.0, E.A_BUMP, None), (3e-2, 2.0, 0.006, None, None),
            (1e-2, 0.0, 0.01, None, (5.0, 3.0, 4.0)), (0.1, 2.0, 0.003, E.A_BUMP, (5.0, 3.0, 4.0))]
    for g, k, fl, A, box in pars:
        loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=g, fixman=k, fixman_tau=TAU, fixman_floor=fl,
                                         fixman_A=(None if A is None else A.copy()), box=box)
        dT, db = BatchDither(loop).drift(Th, beta)
        for i in range(len(Th)):
            ref = loop.drift(Th[i], beta[i])
            assert np.allclose(dT[i], ref[0], rtol=0, atol=1e-14) and np.allclose(db[i], ref[1], rtol=0, atol=1e-14)
    # per-run parameters
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1.0, fixman=1.0, fixman_tau=1.0, fixman_A=E.A_BUMP.copy())
    gs, ks, ts, fs = [1e-2, 0.1, 0.0, 3e-2, 0.1, 1e-2], [0.0, 1.0, 2.0, 2.0, 1.0, 0.5], [1e-4, 5e-5, 1e-4, 2e-4, 1e-4, 1e-4], [0.01, 0.006, 0.0, 0.004, 0.01, 0.0]
    bd = BatchDither(loop, gamma=gs, fixman=ks, tau_F=ts, floor=fs)
    dT, db = bd.drift(Th, beta)
    sub = bd.take(np.array([4, 1]))
    dT2, db2 = sub.drift(Th[[4, 1]], beta[[4, 1]])
    for i in range(len(Th)):
        ref = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=gs[i], fixman=ks[i], fixman_tau=ts[i], fixman_floor=fs[i],
                                        fixman_A=E.A_BUMP.copy()).drift(Th[i], beta[i])
        assert np.allclose(dT[i], ref[0], rtol=0, atol=1e-14) and np.allclose(db[i], ref[1], rtol=0, atol=1e-14)
    assert np.array_equal(dT2, dT[[4, 1]]) and np.array_equal(db2, db[[4, 1]])


def test_a_batch_reproduces_run_dithered_and_resumes_from_its_checkpoint(tmp_path):
    net, base = E.setup()
    Th0 = E.exact_rep(net, np.ones(3) / 3)
    loops = [MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=g, fixman=k, fixman_tau=TAU, fixman_A=E.A_BUMP.copy(), fixman_floor=0.006)
             for g, k in ((0.1, 2.0), (1e-2, 1.0), (0.0, 0.0))]
    ref = [l.run_dithered(Th0, np.ones(3), 400.0, DT, TAU, np.random.default_rng(31 + i), rec_every=10 ** 9) for i, l in enumerate(loops)]
    bd = BatchDither(loops[0], gamma=[0.1, 1e-2, 0.0], fixman=[2.0, 1.0, 0.0], tau_F=[TAU] * 3, floor=[0.006] * 3)
    runs = [dict(key=i, seed=31 + i, n=1600, sq=np.sqrt(2 * TAU * DT * E.G), Th0=Th0, lb0=np.zeros(3)) for i in range(3)]
    ck = str(tmp_path / 'ck.pkl')
    st = run_batch(bd, runs, DT, 400, ckpt=ck, wall=-1.0)                    # stops after one block
    assert not any(v['done'] for v in st.values()) and all(v['k'] == 400 for v in st.values())
    st = run_batch(bd, runs, DT, 400, ckpt=ck)
    for i in range(3):
        assert st[i]['done'] and st[i]['k'] == 1600 and st[i]['exit'] is None
        assert np.allclose(st[i]['Th'], ref[i][0], rtol=0, atol=1e-12) and np.allclose(np.exp(st[i]['lb']), ref[i][1], rtol=0, atol=1e-12)
    # a raised cap continues the finished runs from the checkpoint: the same as one longer run
    long = loops[1].run_dithered(Th0, np.ones(3), 600.0, DT, TAU, np.random.default_rng(32), rec_every=10 ** 9)
    runs2 = [dict(r, n=2400) for r in runs]
    st2 = run_batch(bd, runs2, DT, 400, ckpt=ck)
    assert st2[1]['k'] == 2400 and np.allclose(st2[1]['Th'], long[0], rtol=0, atol=1e-12)


def test_records_and_exits():
    net, base = E.setup()
    Th0 = E.exact_rep(net, np.ones(3) / 3)
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=0.0)
    runs = [dict(key=i, seed=50 + i, n=800, sq=np.sqrt(2 * 3e-3 * DT * E.G), Th0=Th0, lb0=np.zeros(3), acc0=[]) for i in range(4)]
    inside = lambda T: np.abs(T[:, 1, 2] - Th0[1, 2]) < 0.15                  # a slab in the second width
    rec = lambda acc, k, T, lb: acc + [(k, T[1, 2])]
    st = run_batch(BatchDither(loop), runs, DT, 200, inside=inside, rec_every=100, record=rec)
    for i in range(4):
        path = []
        loop.run_dithered(Th0, np.ones(3), 800 * DT, DT, 3e-3, np.random.default_rng(50 + i), rec_every=1, callback=lambda s_, T, b: path.append(T[1, 2]))
        path = np.array(path)
        out = np.flatnonzero(np.abs(path - Th0[1, 2]) >= 0.15)
        if out.size:
            assert st[i]['exit'] == out[0] and st[i]['done']
            assert [k for k, _ in st[i]['acc']] == list(range(0, out[0], 100))
        else:
            assert st[i]['exit'] is None and [k for k, _ in st[i]['acc']] == list(range(0, 801, 100))
        assert np.allclose([v for _, v in st[i]['acc']], path[[k for k, _ in st[i]['acc']]], rtol=0, atol=1e-12)
    assert sum(st[i]['exit'] is not None for i in range(4)) >= 1


def test_exits_survive_a_checkpoint_and_resume(tmp_path):
    # runs that exited before the stop keep their exit; the others resume where they stopped (the same exits as one uninterrupted call)
    net, base = E.setup()
    Th0 = E.exact_rep(net, np.ones(3) / 3)
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=0.0)
    inside = lambda T: np.abs(T[:, 1, 2] - Th0[1, 2]) < 0.4
    runs = [dict(key=i, seed=60 + i, n=1200, sq=np.sqrt(2 * 3e-3 * DT * E.G), Th0=Th0, lb0=np.zeros(3)) for i in range(6)]
    one = run_batch(BatchDither(loop), runs, DT, 200, inside=inside)
    ck = str(tmp_path / 'ck.pkl')
    import pickle
    st = run_batch(BatchDither(loop), runs, DT, 200, ckpt=ck, inside=inside, wall=0.0)
    for _ in range(20):
        if all(v['done'] for v in st.values()):
            break
        st = run_batch(BatchDither(loop), runs, DT, 200, ckpt=ck, inside=inside, wall=0.0)
    assert all(v['done'] for v in st.values())
    assert [st[i]['exit'] for i in range(6)] == [one[i]['exit'] for i in range(6)]
    assert any(one[i]['exit'] is not None for i in range(6)) and any(one[i]['exit'] is None for i in range(6))
    assert pickle.load(open(ck, 'rb'))[2]['exit'] == one[2]['exit']


def test_the_batched_drift_where_the_cohesion_the_box_and_the_floor_act():
    # configurations where the stored runs' terms switch on: particles within the cohesion's range of each other (in the metric G⁻¹),
    # positions outside the box (the dead zone pulls them back), and masses near a saturation, where the floored denominator binds
    net, base = E.setup()
    rng = np.random.default_rng(11)
    T0 = E.exact_rep(net, np.ones(3) / 3)
    close = np.array([np.tile(T0[1], (3, 1)) + 0.15 * rng.standard_normal((3, 3)) for _ in range(3)])
    out = np.array([T0 + np.array([[6.5, 0.2, -0.1], [-7.0, 3.6, 0.3], [0.5, -0.4, 4.8]]) * s for s in (1.0, 1.2, 0.9)])
    Th = np.r_[close, out]
    beta = 3 * np.array([[0.2505, 0.37, 0.3795], [0.33, 0.176, 0.494], [0.34, 0.36, 0.30],
                         [0.2503, 0.30, 0.4497], [0.30, 0.1755, 0.5245], [0.33, 0.33, 0.34]])
    d2 = np.sum((close[:, :, None, :] - close[:, None, :, :]) ** 2 / np.asarray(E.G)[None, None, None, :], axis=3)
    assert np.any((d2 > 0) & (d2 < 1.0))                                     # within the cohesion's range R = 1
    box = np.array([5.0, 3.0, 4.0])
    assert np.all(np.any(np.abs(out) > box, axis=(1, 2)))
    for g, k, fl, A, bx in ((1e-2, 2.0, 0.01, E.A_BUMP, box), (0.1, 1.0, 0.006, E.A_BUMP, None), (3e-2, 2.0, 0.01, None, box)):
        loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=g, fixman=k, fixman_tau=TAU, fixman_floor=fl,
                                         fixman_A=(None if A is None else A.copy()), box=(None if bx is None else tuple(bx)))
        bd = BatchDither(loop)
        dT, db = bd.drift(Th, beta)
        for i in range(len(Th)):
            ref = loop.drift(Th[i], beta[i])
            assert np.allclose(dT[i], ref[0], rtol=0, atol=1e-14) and np.allclose(db[i], ref[1], rtol=0, atol=1e-14)
        if A is not None:                                                    # the floor binds at the masses next to the first face
            W = beta / 3
            q = (np.abs(A)[None, :] / net.a_max / W) ** 2
            assert np.any(W * (1 - q) < fl * W)
            c = bd.correction(np.zeros_like(W), beta)
            assert np.allclose(c, k * (TAU / 3) * (1 + q) / np.maximum(W * (1 - q), fl * W), rtol=1e-13)
