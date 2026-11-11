"""tools/mutants.py with another interpreter for the mutants' tests (the repository's own runner hard-codes `venv/bin/python`): for experiment scripts whose tests need packages that are not in the repository's
environment (the astronomy baselines need astropy, photutils and crowdsource). Same arguments as tools/mutants.py, plus the environment variable MUTANT_PYTHON, the interpreter (a script that execs it is fine;
it must be able to run `-m pytest`).

    MUTANT_PYTHON=/path/to/interpreter venv/bin/python tools/mutants_py.py --exp exp_astro_baselines_hp --tests tests/test_exp_astro_baselines_hp.py ... --spec tools/mutants/wp_astrobaseshp_own.py --work <dir>

The file tools/mutants.py is read, its one line `PY = ...` replaced, and the result executed as that file (so REPO and the rest are unchanged); nothing is written.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PATH = os.path.join(HERE, 'mutants.py')
OLD = "PY = os.path.join(REPO, 'venv', 'bin', 'python')"
src = open(PATH).read()
if src.count(OLD) != 1:
    raise SystemExit('tools/mutants.py has changed: the line %r is not there exactly once' % OLD)
if not os.environ.get('MUTANT_PYTHON'):
    raise SystemExit('set MUTANT_PYTHON to the interpreter that runs the tests')
exec(compile(src.replace(OLD, "PY = os.environ['MUTANT_PYTHON']"), PATH, 'exec'), {'__name__': '__main__', '__file__': PATH})
