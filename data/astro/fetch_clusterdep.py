#!/usr/bin/env python3
"""Data access for the cluster-dependence held-out test (`results/astro/clusterdep_plan.md`): psField (0.3 MB) and r-band frame (3 MB) of the frames of further HUGS-covered globular clusters.

A thin companion of `data/astro/fetch_fresh.py` (whose helpers it imports; `fetch_fresh.py` and `checksums_fresh.json` are not touched): only the checksum file (`data/astro/checksums_clusterdep.json`)
differs; the download root is `$ASTRO_DATA_4` (default /home/user/astro_data_4, never inside the repository) so that the frames lie beside the earlier fresh ones. Standard library only; the host is data.sdss.org
(SDSS SAS dr17) and nothing else; a refusal (401/403/404) is reported and never worked around. The HUGS catalogues and the photoField files are the ones that `fetch_second.py` recorded in `checksums_second.json`.

    python -I data/astro/fetch_clusterdep.py psfield --frames ngc6934:1738:1:29 ...
    python -I data/astro/fetch_clusterdep.py frame   --frames ngc6934:1738:1:29 ...
    python -I data/astro/fetch_clusterdep.py verify  --frames ...
"""
import argparse
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location('fetch_fresh', os.path.join(HERE, 'fetch_fresh.py'))
FF = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(FF)

CHECKSUMS = os.path.join(HERE, 'checksums_clusterdep.json')
FF.CHECKSUMS = CHECKSUMS
FF.FS.CHECKSUMS = CHECKSUMS


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('kind', choices=('psfield', 'frame', 'verify'))
    ap.add_argument('--frames', nargs='+', required=True, help='cluster:run:camcol:field')
    a = ap.parse_args()
    frames = FF.parse_frames(a.frames)
    if a.kind == 'verify':
        return FF.verify(FF.items_for(frames, ('psfield', 'frame')), sums_path=CHECKSUMS)
    FF.fetch(FF.items_for(frames, (a.kind,)))


if __name__ == '__main__':
    sys.exit(main() or 0)
