"""
Export the frames of the SMLM benchmark for a run of DECODE (or any other GPU method) on the authors' machine
==============================================================================================================

Writes, for each set, under <out>/<set>/:  dev.tif and test.tif (16-bit TIFF stacks of the development and test frames, in DN, exactly as the challenge's stack: baseline 100 DN),
frames.json (the 0-based stack indices of each file, the camera and the protocol), dev_gt.csv (the ground truth of the DEVELOPMENT frames only: frame, x_nm, y_nm, z_nm, photons;
frame 0-based, the test frames' ground truth is never exported). The raw data stay out of the repository: <out> is an ordinary data directory.

    python experiments/smlm_decode_export.py --out /path/to/decode_run [--sets MT0.N1.LD MT0.N1.HD MT0.N2.LD MT0.N2.HD]

The protocol (shared with `exp_smlm_iss.py` and `exp_smlm_picasso.py`): LD sets, development frames 0-999 and test frames 1000-1499; HD sets, development 0-99, test 100-499.
`results/smlm/decode_run_spec.md` says what to run on these files and what to return.
"""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, HERE)
import smlm_common as C     # noqa: E402

PROTOCOL = {'LD': dict(dev=(0, 1000), test=(1000, 1500)), 'HD': dict(dev=(0, 100), test=(100, 500))}
CAMERA = dict(baseline_dn=C.BASELINE_DN, em_gain=C.EM_GAIN, e_per_adu=C.E_PER_ADU, qe=C.QE, read_sigma_e=74.4, spur_noise=0.002, pixel_nm=C.PIXEL_NM, size_px=64)


def protocol(name):
    return PROTOCOL['LD' if name.endswith('LD') else 'HD']


def write_tiff(path, frames):
    """A multi-page 16-bit TIFF of frames (n, 64, 64) with integer DN values, written with Pillow."""
    from PIL import Image
    a = np.asarray(frames)
    if not np.all(a == np.round(a)) or a.min() < 0 or a.max() > 65535:
        raise ValueError('frames are not 16-bit integers')
    ims = [Image.fromarray(f.astype(np.uint16)) for f in a]
    ims[0].save(path, save_all=True, append_images=ims[1:], compression=None)


def export_set(name, out, root=None):
    d = os.path.join(out, name)
    os.makedirs(d, exist_ok=True)
    stack = C.data_path(name, 'stack', root)
    gt = C.load_gt(C.data_path(name, 'gt', root))
    meta = dict(set=name, camera=CAMERA, indexing='0-based indices into the original stack', files={})
    for which in ('dev', 'test'):
        a, b = protocol(name)[which]
        fr = C.load_frames(stack, a, b)
        write_tiff(os.path.join(d, which + '.tif'), fr)
        meta['files'][which] = dict(file=which + '.tif', first=a, last=b - 1, n=b - a)
    a, b = protocol(name)['dev']
    g = gt[(gt['frame'] >= a) & (gt['frame'] < b)]
    with open(os.path.join(d, 'dev_gt.csv'), 'w') as f:
        f.write('frame,x_nm,y_nm,z_nm,photons\n')
        for r in g:
            f.write('%d,%.3f,%.3f,%.3f,%.3f\n' % (r['frame'] - a, r['x_nm'], r['y_nm'], r['z_nm'], r['photons']))
    meta['dev_gt'] = dict(file='dev_gt.csv', frame_index='0-based within dev.tif', rows=int(len(g)))
    json.dump(meta, open(os.path.join(d, 'frames.json'), 'w'), indent=1)
    return meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True)
    ap.add_argument('--sets', nargs='*', default=list(C.SETS))
    A = ap.parse_args()
    for name in A.sets:
        m = export_set(name, A.out)
        print('%s: dev %d frames, test %d frames, %d ground-truth rows of the development frames' % (name, m['files']['dev']['n'], m['files']['test']['n'], m['dev_gt']['rows']))


if __name__ == '__main__':
    main()
