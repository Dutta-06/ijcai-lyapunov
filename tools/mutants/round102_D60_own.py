# (id, what, old, new, k): single-point mutants of exp_noise_resonance.py — round 102's new code (state_start, c8a_job's rerun with the
# metric, the one-ingredient variants and the noise-free probe, c8a_lines' computed conclusions, the anchored zero raise, the C8R clause).
M = [
 ('X00', 'control', 'C8R_BACK = (20, 40, 80) ', 'C8R_BACK = (20, 40, 80) ', 0),
 ('U01', 'state_start: inst coefficient from the mean', "al = raised(L, al, idx, k, abs(float(W[p - 1][k] if cw == 'mean' else Wi[k])) - w0)", "al = raised(L, al, idx, k, abs(float(W[p - 1][k])) - w0)", 0),
 ('U02', 'state_start: inst centre from the mean', "c[idx[k]] = float(C[p - 1][k] if cc == 'mean' else Ci[k])", "c[idx[k]] = float(C[p - 1][k])", 0),
 ('U03', 'state_start: the mean of period p', "abs(float(W[p - 1][k] if cw == 'mean' else Wi[k]))", "abs(float(W[p][k] if cw == 'mean' else Wi[k]))", 0),
 ('U04', 'c8a_job: the rerun without the metric', "r, L, al0, c0, idx, Wp, Cp = c_run(m, ke, s, dict(opts, metric=mt))   # round 102", "r, L, al0, c0, idx, Wp, Cp = c_run(m, ke, s, dict(opts))   # round 102", 0),
 ('U05', 'c8a_job: the instantaneous state one output late', "Wi, Ci = [float(v) for v in r['w'][p * kk, idx]], [float(v) for v in r['c'][p * kk, idx]]", "Wi, Ci = [float(v) for v in r['w'][p * kk + 1, idx]], [float(v) for v in r['c'][p * kk + 1, idx]]", 0),
 ('U06', 'c8a_job: c1 of kernel 2', "c1=float(C[p - 1][0] - CT.CSTAR[0])", "c1=float(C[p - 1][1] - CT.CSTAR[1])", 0),
 ('U07', 'the quiet probe under the noise', "res=list(edge_run(m, ke, mt, 0, 0.0, NP=NP_C7, start=offsets_start(W, C, pp, False, sub))))", "res=list(win_run(m, ke, mt, s, opts, 0.0, pp, start=offsets_start(W, C, pp, False, sub))))", 0),
 ('U08', 'variants: w+c1 with kernel 2', "C8A_VARIANTS = (('w+c1', (('mean', 'mean'), ('mean', 'fit'), ('mean', 'fit')))", "C8A_VARIANTS = (('w+c1', (('mean', 'fit'), ('mean', 'mean'), ('mean', 'fit')))", 0),
 ('U09', 'lines: no-alarm clause without the check', "    if all(c == 'kept' for c in oc['w']):", "    if True:", 0),
 ('U10', 'lines: the centre clause half-checked', "if 'w+c1' in oc and oc['w+c1'] == oc['wc'] and oc['w+c23'] == oc['w']:", "if 'w+c1' in oc and oc['w+c1'] == oc['wc']:", 0),
 ('U11', 'lines: the instantaneous clause unchecked', "if ki and all(oc['i'][k_] == 'removed' for k_ in ki):", "if ki:", 0),
 ('U12', 'lines: probe clause unchecked', "if all(l_ == (2 in o) for o, l_ in has3) else '')", "if True else '')", 0),
 ('U13', 'lines: the collision clause unchecked', "and any(x['only'] == [2] and x['res'][0] >= 0 for x in pb['rows']) else ''))", "or True else ''))", 0),
 ('U14', 'lines: kept as delayed', "('kept' if lp == red else ('delayed' if lp > red else 'earlier'))", "('kept' if lp >= red else 'earlier')", 0),
 ('U15', 'anchored zero raise not counted', "                ancr0 += 1\n", "                pass\n", 0),
 ('U16', 'C8R clause whatever the passes', "if cov and all(z['c8r'] for z in cov) else ''))", "if cov else ''))", 0),
]
