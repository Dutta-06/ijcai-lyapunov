"""The astronomy report end to end on a synthetic raw (the keys part_T and part_V write), so that a bug in the report cannot hide until the real runs are done."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'experiments'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import astro_common as A      # noqa: E402
import exp_astro_m2 as X      # noqa: E402


def _synthetic_catalogs(rng):
    """A truth list of 60 stars in a 100 x 100 patch and four catalogues of different quality (each (x, y, flux))."""
    pos = rng.uniform(5, 95, (60, 2))
    mag = rng.uniform(17.0, 22.6, 60)
    truth = dict(x=pos[:, 0], y=pos[:, 1], mag=mag, inner=np.ones(60, bool))
    cats = {}
    for name, miss, fake, jit in (('good', 0.1, 3, 0.1), ('mid', 0.3, 8, 0.2), ('poor', 0.5, 15, 0.3), ('bic', 0.15, 4, 0.1)):
        keep = rng.uniform(size=60) > miss
        xy = pos[keep] + rng.normal(0, jit, (keep.sum(), 2))
        fl = A.mag_to_counts(mag[keep] + rng.normal(0, 0.1, keep.sum()))
        fx = rng.uniform(5, 95, (fake, 2))
        ff = A.mag_to_counts(rng.uniform(20.0, 23.0, fake))
        cats[name] = (np.r_[xy[:, 0], fx[:, 0]], np.r_[xy[:, 1], fx[:, 1]], np.r_[fl, ff])
    return truth, cats


def test_report_runs_on_a_synthetic_raw_and_is_deterministic(monkeypatch):
    rng = np.random.default_rng(5)
    truth, cats = _synthetic_catalogs(rng)
    pxy = (310, 630)
    monkeypatch.setattr(X, 'TESTS', [pxy, (410, 630)])
    kap = 5.0
    raw = {}
    monkeypatch.setattr(X, 'SECOND', {'m53': dict(title='M53 (NGC 5024)', root='/nonexistent', prefix='p', truth='t.npz', patches=[(1, 2)])})
    groups = [('T|x%d_y%d|', p) for p in X.TESTS] + [('T2|m53|x%d_y%d|', (1, 2))]
    for pi, (kfmt, p) in enumerate(groups):
        pre = kfmt % p
        raw[pre + 'cfg'] = np.array([5.0, 1.2, kap])
        for label, key in (('ours', 'good'), ('default', 'mid')):
            raw[pre + label + '|secs'] = np.array([12.0, 90.0])
            for k in X.KAPPAS:
                for kk, v in X._pack('%s|k%g' % (label, k), A.score(np.c_[cats[key][0], cats[key][1]], cats[key][2], truth)).items():
                    raw[pre + kk] = v
        for kk, v in X._pack('ours|k0', A.score(np.c_[cats['bic'][0], cats['bic'][1]], cats['bic'][2], truth)).items():
            raw[pre + kk] = v
        for kk, v in X._pack('default|k0', A.score(np.c_[cats['bic'][0], cats['bic'][1]], cats['bic'][2], truth)).items():
            raw[pre + kk] = v
        for label, key in (('photutils (default)', 'poor'), ('photutils (tuned)', 'mid'), ('crowdsource (default)', 'poor'), ('crowdsource (tuned)', 'mid')):
            for kk, v in X._pack(label, A.score(np.c_[cats[key][0], cats[key][1]], cats[key][2], truth)).items():
                raw[pre + kk] = v
            if '(tuned)' in label:
                raw[pre + label.split(' ')[0] + ' (tuned)|name'] = np.array(['variant'])
        raw[pre + 'ours|bic'] = np.array([10, 20, 0])
        raw[pre + 'default|bic'] = np.array([19, 20, 1])
        if pi == 0:
            for kk, v in X._pack('daophot (as shipped)', A.score(np.c_[cats['poor'][0], cats['poor'][1]], cats['poor'][2], truth)).items():
                raw[pre + kk] = v
            for kk, v in X._pack('daophot (published)', A.score(np.c_[cats['poor'][0], cats['poor'][1]], cats['poor'][2], truth)).items():
                raw[pre + kk] = v
            for kk, v in X._pack('greedy_snr|k5', A.score(np.c_[cats['good'][0], cats['good'][1]], cats['good'][2], truth)).items():
                raw['V|x%d_y%d|%s' % (p[0], p[1], kk)] = v
            raw['V|x%d_y%d|greedy_snr|secs' % p] = np.array([10.0, 80.0])
            for k in X.KAPPAS:
                for kk, v in X._pack('greedy_snr|k%g' % k, A.score(np.c_[cats['good'][0], cats['good'][1]], cats['good'][2], truth)).items():
                    raw['V|x%d_y%d|%s' % (p[0], p[1], kk)] = v
    md = X.report(raw)
    assert md == X.report(raw)
    for needle in ('M2 (NGC 7089), SDSS run 2583: test patch (310, 630)', 'M2 (NGC 7089), SDSS run 2583: test patch (410, 630)', 'M2 (NGC 7089), SDSS run 2583: all 2 test patches together', 'M53 (NGC 5024): test patch (1, 2)', 'a priori default with the BIC stop', 'y shifted by -1 px (the shift was fitted to the HUGS truth of this patch)', 'DAOPHOT catalogue of the benchmark as shipped (no shift)', 'unmatched (all mags)', 'BIC readout: block 10 of 20 (interior)', 'blocks of 20 x 20 px without such a source', 'stratified', 
                   'Paired differences', 'Design variants of the path', 'Completeness by magnitude', 'verdict on F', 'wins / ties / losses'):
        assert needle in md, needle
    assert 'nan' not in md.lower().replace('nan,', '')                                    # no NaN cell
    pooled_md = md.split('SDSS run 2583: all 2 test patches together')[1].split('M53 (NGC 5024)')[0]
    assert 'y shifted by -1 px' not in pooled_md                                         # a catalogue of one patch is not pooled
    assert '| photutils (default) |' in pooled_md and '| sliding ISS' in pooled_md        # the others are
    assert '(50 blocks)' in pooled_md.replace('pooled (50 blocks)', '(50 blocks)')       # 2 patches x 25 blocks
