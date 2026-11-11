# (id, what, old, new, k): single-point mutants of gh_lastann.py — round 122's checks (the joint rebuilt from the pieces certified, Y's
# enclosure proved on sub-boxes, the tiling by containment and disjointness, the soundness corners). Run with --exp gh_lastann against
# tests/test_gh_peano.py -k part_va (gh_lastann writes no md: only the tests decide).
M = [
 ('P00', 'control', "def _chain(pieces, A, B):", "def _chain(pieces, A, B):", 0),
 ('P01', "the connector's cover unchecked", "            if not _chain(pc[:ncon], complex(delta0, lo), complex(delta0, hi)):", "            if False:", 0),
 ('P02', "the polyline's start unchecked", "            if not (J['pts'] and J['pts'][0] == complex(delta0, J['yc'])):", "            if False:", 0),
 ('P03', 'disjointness unchecked', "    disjoint = all(min(p[1], q[1]) <= max(p[0], q[0])", "    disjoint = True or all(min(p[1], q[1]) <= max(p[0], q[0])", 0),
 ('P04', "Y's enclosure not computed", "        if ok(r1, r2, t1, t2):", "        if True:", 0),
 ('P05', 'the soundness corners dropped', "            for j in range(per + 4):", "            for j in range(per):", 0),
 ('P06', "κ's clamp at 0 dropped", "            k_lo = k_lo if k_lo > 0 else arb(0)", "            k_lo = k_lo", 0),
 ('P07', 'the tiling by the area alone', "tiled=bool(rel < 1e-12 and inside and disjoint))", "tiled=bool(rel < 1e-12))", 0),
 ('P08', "the polyline's vertices unchecked", " and all(p_ in verts for p_ in J['pts'])):", "):", 0),
 ('P09', 'y_c unchecked', "            if not (lo <= J['yc'] <= hi):", "            if False:", 0),
]
