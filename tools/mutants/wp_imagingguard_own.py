"""Single-point mutants of src/core/imaging_guard.py (the cycle guard of the sliding path), for tools/mutants.py:

    venv/bin/python tools/mutants.py --exp ../src/core/imaging_guard --tests tests/test_imaging_guard.py --raw results/smlm/calibration.md --md results/smlm/calibration.md \\
        --out-env MUT_OUT --spec tools/mutants/wp_imagingguard_own.py --work <scratch dir> [--slot a --i0 0 --step 2]

(the tool's md / `--from-raw` columns mean nothing for these mutants: only `tests` does, True: tests/test_imaging_guard.py kills the mutant). A tuple is (id, what, old, new[, k]): `old`
replaced by `new` at its k-th occurrence (0-based) in src/core/imaging_guard.py.

The mutants are the guard's own lines (the two-strike rule and its comparison, the tolerance, the record of the last exit, the application of the bar to `blocked`, the append to `barred`, the
switch `cycle_guard`, the permanent bar) and three lines of the cycle itself that the guard copies from `imaging.sliding_path` (the exit resets the dual and releases the point).

Result (tests/test_imaging_guard.py, 34 tests, 39 mutants): 34 are killed, 5 survive, all five equivalent on every path the code can produce:
  G04  `<` instead of `<=` in the comparison of the two exit times: the slack 1e-12 x max(1, tau) is positive, and the repeats of a cycle are at exactly the same tau (the re-entry is at dt = 0), so
       the difference is 0 and 0 < slack holds as 0 <= slack does;
  G10  slack 0: the same exact repeats (difference 0 <= 0); the slack only matters for two exits that differ by less than 1e-12, which the dual dynamics do not produce at a re-entry (dt is exactly 0);
  G11  slack 1e-12 (absolute, not scaled with tau) and G12 slack 1e-12 x tau (purely relative): the same for exact repeats; they differ from the code only for two different exit times closer than
       1e-12 x max(1, tau), and the test worlds have none (the closest pair planted is 4.5e-11 apart, `test_two_exits_of_one_point_tens_of_picoseconds_apart_are_still_two_entry_times`);
  G21  the bar applied to `blocked` even with the switch off: `barred` is only ever set under the switch, so with the switch off it is all False and `taken | barred` is `taken`.
  Also worth saying: G24 (the bar recorded in `taken`, which blocks the point just the same) is behaviourally equivalent and is killed only by the test that compares the source with
  `imaging.sliding_path` (its new line does not mention the guard); G13 (a first exit compared with entry time 0) is killed by the tie test, whose entry times are 5e-15 (below the slack).
  The mutants whose old text is a line copied from `imaging.sliding_path` (H01 to H03) are killed by the planted-cycle tests: without the reset of the dual, the release of the point or the
  exit condition there is no cycle to guard.
"""
TOL = "abs(tau - last_exit.get(ids[i], -np.inf)) <= 1e-12 * max(1.0, abs(tau))"
M = [
    ('G01', 'one strike bars the point (the double-exit condition always true)', "if " + TOL + ":", "if True:"),
    ('G02', 'the comparison is == instead of <= (the slack is never met exactly)', TOL, "abs(tau - last_exit.get(ids[i], -np.inf)) == 1e-12 * max(1.0, abs(tau))"),
    ('G03', 'the comparison is >= (bars exits far apart in time)', TOL, "abs(tau - last_exit.get(ids[i], -np.inf)) >= 1e-12 * max(1.0, abs(tau))"),
    ('G04', 'the comparison is < (strict)', TOL, "abs(tau - last_exit.get(ids[i], -np.inf)) < 1e-12 * max(1.0, abs(tau))"),
    ('G05', 'tolerance 1.0 (an absolute slack of 1 below tau 1)', "<= 1e-12 * max(1.0, abs(tau))", "<= 1.0 * max(1.0, abs(tau))"),
    ('G06', 'tolerance 1e-3', "<= 1e-12 * max(1.0, abs(tau))", "<= 1e-3 * max(1.0, abs(tau))"),
    ('G07', 'tolerance 1e-5', "<= 1e-12 * max(1.0, abs(tau))", "<= 1e-5 * max(1.0, abs(tau))"),
    ('G08', 'tolerance 1e-7', "<= 1e-12 * max(1.0, abs(tau))", "<= 1e-7 * max(1.0, abs(tau))"),
    ('G09', 'tolerance 1e-9', "<= 1e-12 * max(1.0, abs(tau))", "<= 1e-9 * max(1.0, abs(tau))"),
    ('G10', 'tolerance 0 (exact equality of the two entry times)', "<= 1e-12 * max(1.0, abs(tau))", "<= 0.0 * max(1.0, abs(tau))"),
    ('G11', 'the slack is absolute (no scaling with tau)', "<= 1e-12 * max(1.0, abs(tau))", "<= 1e-12 * 1.0"),
    ('G12', 'the slack is purely relative to tau', "<= 1e-12 * max(1.0, abs(tau))", "<= 1e-12 * abs(tau)"),
    ('G13', 'the first exit of a point is compared with entry time 0', "last_exit.get(ids[i], -np.inf)", "last_exit.get(ids[i], 0.0)"),
    ('G14', 'the last exit is never recorded', "                last_exit[ids[i]] = tau\n", "                pass\n"),
    ('G15', 'the last exit is recorded only when the point is barred', "                    bars.append((tau, ids[i]))\n                last_exit[ids[i]] = tau\n", "                    bars.append((tau, ids[i]))\n                    last_exit[ids[i]] = tau\n"),
    ('G16', 'the exit is recorded before it is compared (a point is compared with itself)', "            if cycle_guard:\n                if abs(", "            if cycle_guard:\n                last_exit[ids[i]] = tau\n                if abs("),
    ('G17', 'the strikes are counted for all points together (one global last exit)', "last_exit.get(ids[i], -np.inf)", "last_exit.get(0, -np.inf)", 0),
    ('G18', 'the strikes are counted for all points together (the record)', "                last_exit[ids[i]] = tau\n", "                last_exit[0] = tau\n"),
    ('G19', 'the bar is not applied to the entry rule', "        blocked = taken | barred if cycle_guard else taken\n", "        blocked = taken\n"),
    ('G20', 'the bar is combined with the taken points by and, not or', "blocked = taken | barred if cycle_guard else taken", "blocked = taken & barred if cycle_guard else taken"),
    ('G21', 'the bar is applied even with the switch off', "blocked = taken | barred if cycle_guard else taken", "blocked = taken | barred"),
    ('G22', 'the bar is lost when the exclusion zone is added to `blocked`', "            blocked = blocked | lat.near(c[:, :2], r_ex)", "            blocked = taken | lat.near(c[:, :2], r_ex)"),
    ('G23', 'the bar is never set', "                    barred[ids[i]] = True\n", "                    pass\n"),
    ('G24', 'the bar is set on the array of the taken points (the point stays blocked: equivalent)', "                    barred[ids[i]] = True\n", "                    taken[ids[i]] = True\n"),
    ('G25', 'the bar is set on the wrong lattice point (the neighbour in a)', "                    barred[ids[i]] = True\n", "                    barred[ids[i][0], max(ids[i][1] - 1, 0)] = True\n"),
    ('G26', 'the bar is lifted again at the next pass (it is not permanent)', "        blocked = taken | barred if cycle_guard else taken\n", "        blocked = taken | barred if cycle_guard else taken\n        barred[:] = False\n"),
    ('G27', 'the point is not appended to the list of the barred', "                    bars.append((tau, ids[i]))\n", "                    pass\n"),
    ('G28', 'the barred point is listed without its time', "bars.append((tau, ids[i]))", "bars.append(ids[i])"),
    ('G29', 'the barred point is listed with entry time 0', "bars.append((tau, ids[i]))", "bars.append((0.0, ids[i]))"),
    ('G30', 'the barred point is listed with its indices swapped', "bars.append((tau, ids[i]))", "bars.append((tau, ids[i][::-1]))"),
    ('G31', 'the list of the barred is not returned', "snr=np.array(snrs), barred=bars)", "snr=np.array(snrs), barred=[])"),
    ('G32', 'the list of the barred is not a key of the result', ", barred=bars)", ")"),
    ('G33', 'the strikes are counted and the bars set even when the switch is off', "            if cycle_guard:\n                if abs(", "            if True:\n                if abs("),
    ('G34', 'the strikes are counted only when the switch is off', "            if cycle_guard:\n                if abs(", "            if not cycle_guard:\n                if abs("),
    ('G35', 'the default of the switch is off', "tau_max=None, cycle_guard=True):", "tau_max=None, cycle_guard=False):"),
    ('G36', 'every lattice point starts barred', "    barred = np.zeros(lat.shape, bool)\n", "    barred = np.ones(lat.shape, bool)\n"),
    ('H01', 'the exit does not put the dual back on the wall (copied line)', "            U[ids[i]] = rho\n", ""),
    ('H02', 'the exit does not release the lattice point (copied line)', "            taken[ids[i]] = False\n", ""),
    ('H03', 'the exit test reads the correlation with the other sign (copied line)', "gone = (w <= 1e-12 * max(w.max(), 1e-300)) & (g <= vtol * vscale)", "gone = (w <= 1e-12 * max(w.max(), 1e-300)) & (g >= vtol * vscale)"),
]
