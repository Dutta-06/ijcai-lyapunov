# (id, what, old, new, k): single-point mutants of exp_transit_certificate.py — Part K after round 149 (D3er-6, D3er-11): `k_classify`
# (Part I's failures by type from the stored solutions: lost, spurious, displaced, wrong sign), `part_K` (the classes, the displacement
# radius K_RD·σ, Part G's entry sequences) and the K section of `report`.
# Run with -k "k_classify or part_K or report_K" and --extra results/chaos/transit_certificate.md (the tests recount printed rows).
# Result: 9 mutants and a control (the control survives: 4 tests pass, the md regenerates). All 9 killed by the tests; Q08 and Q09 also
# change the md. Q09 survived the first run (the md's 640 was checked against the stored md only) and is killed by
# test_report_K_on_planted_rows as extended (the regenerated section's 640). Q05, Q06, Q07 are equivalent on the data (no wrong sign, every
# sequence and class equal) and are killed by planted rows. No survivor.
M = [
 ('Q00', 'control', "def k_classify(sol, cs, As, tol, rd):", "def k_classify(sol, cs, As, tol, rd):", 0),
 ('Q01', "displaced: every nonzero beyond the tolerance", "disp, spur = dm[(dm > tol) & (dm <= rd)], dm[dm > rd]",
  "disp, spur = dm[dm > tol], dm[dm > rd]", 0),
 ('Q02', "lost: no nonzero within the tolerance (not the radius)", "if not np.any((near == k) & (dm <= rd))]",
  "if not np.any((near == k) & (dm <= tol))]", 0),
 ('Q03', "displaced before spurious", "'spurious' if len(spur) else 'displaced' if len(disp)", "'displaced' if len(disp) else 'spurious' if len(spur)", 0),
 ('Q04', "the radius in grid spacings of the plane, not in σ", "tol, rd = I_TOL * H_BP[0], K_RD * TM.setup()[1].sig",
  "tol, rd = I_TOL * H_BP[0], K_RD * H_BP[1]", 0),
 ('Q05', "the wrong-sign check dropped", "sign = any(c > 0 and np.sign(x[(near == k) & (dm <= tol)].sum()) != np.sign(A)",
  "sign = False and any(c > 0 and np.sign(x[(near == k) & (dm <= tol)].sum()) != np.sign(A)", 0),
 ('Q06', "the sequences compared as sets", "seq=[dict(src=s, same=o == f) for s, o, f in seq]", "seq=[dict(src=s, same=sorted(o) == sorted(f)) for s, o, f in seq]", 0),
 ('Q07', "the plane's classes not compared", "same=g['out'] == g['out_flow']) for g in G['plane']", "same=True) for g in G['plane']", 0),
 ('Q08', "the report's failures without the spurious class", "bad = [t for t in ts if t['cls'] != 'rec']",
  "bad = [t for t in ts if t['cls'] in ('lost', 'displaced', 'sign')]", 0),
 ('Q09', "the sequences' own-row count with D(a) for D(b)", "src['D'][x['src'][1]]['tag'][0] == 'b') for x in sq)",
  "src['D'][x['src'][1]]['tag'][0] == 'a') for x in sq)", 0),
]
