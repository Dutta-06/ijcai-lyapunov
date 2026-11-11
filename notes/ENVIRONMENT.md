# Environment

- Python 3.13.18 (CPython), Linux x86-64, CPU only.
- Exact package pins: `tools/requirements-pinned.txt` (numpy, scipy, matplotlib, python-flint for ball arithmetic, scikit-learn, pytest, ...). Astronomy extras: `tools/requirements-astro.txt`. Loose list: `requirements.txt`.
- Scripts are run from the repository root as `python experiments/<name>.py` (most take `--from-raw` to regenerate the report byte-identically from the stored `*_raw.npz`).
- The ball-arithmetic verification (theory3n, theory3y) needs python-flint (arb).
- Tests: `pytest tests/` (the mechreal tests take about 25 s, the theory3y tests about 26 s, theory3n about a minute).
