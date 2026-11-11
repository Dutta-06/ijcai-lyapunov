# (id, what, old, new, k): single-point mutants of exp_gh_peano.py — round 110's revision (the ball bound's two quadratic pieces against a
# dense sup, part DB's stored bound, the reporting guards: open boxes, a missing or positive ball bound, no fallback to the grid's radius,
# D's own bb kept, c_need at the floors in ball arithmetic).
M = [
 ('Z00', 'control', "VC_EPS = (1e-2, 1e-3)", "VC_EPS = (1e-2, 1e-3)", 0),
 ('Z01', "ball bound: the a ≤ 1 piece's roots dropped", "    for qa, qb in ((p + 2 * B * B, -2 * B * B), (p, -2 * B)):\n        if not (qb < 0):", "    for qa, qb in ((p, -2 * B),):\n        if not (qb < 0):", 0),
 ('Z02', "ball bound: the a ≤ 1 piece's quadratic wrong", "    for qa, qb in ((p + 2 * B * B, -2 * B * B), (p, -2 * B)):\n        if not (qb < 0):", "    for qa, qb in ((p + B * B, -2 * B * B), (p, -2 * B)):\n        if not (qb < 0):", 0),
 ('Z03', 'ball bound: the kink dropped', "        if not (wk.upper() <= 0) and not (wk.lower() > wm):\n            cand.append(wk)", "        if False:\n            cand.append(wk)", 0),
 ('Z04', 'part_DB stores the float bound', "    return [dict(R=R, n=n, worst=w, open=o, ball_worst=bw) for R in DOM_BB", "    return [dict(R=R, n=n, worst=w, open=o, ball_worst=w) for R in DOM_BB", 0),
 ('Z05', "d_lines: open boxes ignored", "    bbv = [x['R'] for x in dm.get('bb', []) if x['open'] == 0]", "    bbv = [x['R'] for x in dm.get('bb', [])]", 0),
 ('Z06', "d_lines: the grid's radius as a fallback", "    R_ver = min(max(bbv), R_dom) if bbv else None", "    R_ver = min(max(bbv or [R_dom]), R_dom)", 0),
 ('Z07', "v_lines: open boxes ignored", "    R_bb = max([x['R'] for x in dm.get('bb', []) if x['open'] == 0] or [0.0])", "    R_bb = max([x['R'] for x in dm.get('bb', [])] or [0.0])", 0),
 ('Z08', "v_lines: no guard without a verified branch and bound", "    if R_bb == 0.0:", "    if False:", 0),
 ('Z09', "v_lines: 'both' without the ball bound's sign", "x.get('ball_worst', 0) < 0 and x['open'] == 0 for x in dm.get('bb', [])", "x['open'] == 0 for x in dm.get('bb', [])", 0),
 ('Z10', "report: DB always replaces D's bb", "    if DB is not None and isinstance(D, dict) and D.get('dom') and not (", "    if DB is not None and isinstance(D, dict) and D.get('dom') and not (False and ", 0),
 ('Z11', "d_lines: 'inside' unconditional", "{'inside' if c_need_ball(rs_) < R_ver else 'NOT inside'}", "inside", 0),
 ('Z12', "v_lines: 'inside' unconditional", "{'inside' if c_need_ball(rv) < R_v else 'NOT inside'}", "inside", 0),
 ('Z13', "c_need_ball: the lower end", "return math.ceil(float((4 * u.log() / (u - 1)).sqrt().upper()) * 1e8) / 1e8", "return math.floor(float((4 * u.log() / (u - 1)).sqrt().lower()) * 1e8) / 1e8", 0),
]
