# (id, what, old, new, k): single-point mutants of src/core/gauge_batch.py — the batch integrator added with the runs declined after round 143
# (`BatchDither`: the per-run parameters, `take`, `features`, `drift`, `correction`, `step`; `run_batch`: the per-run generators and their
# blocks of draws, the noise's indexing, the records, the exits, the checkpoints). tools/mutants.py mutates experiments/<exp>.py, so this
# spec is run with `--exp ../src/core/gauge_batch` (the slot's src/core/gauge_batch.py) and `--tests tests/test_gauge_batch.py`; no md
# reads the integrator (the md's regeneration then runs the module itself, which writes nothing: the md column is void, kills are the tests').
# Run: tools/mutants.py --exp ../src/core/gauge_batch --tests tests/test_gauge_batch.py --raw results/chaos/gauge_plateau_raw.npz
# --md results/chaos/gauge_plateau.md --out-env GP_OUT (one worker; ~1–4 s per mutant; the md column is void, the control included).
# Result: 24 mutants and a control. Against the tests as committed: 21 killed by the tests; the 3 others (B10: the floor unfloored; B14: the
# cohesion without the metric; B15: the box's force reversed) killed by test_the_batched_drift_where_the_cohesion_the_box_and_the_floor_act
# (the committed configurations never put two particles within the cohesion's range, a position outside the box, or a mass within the
# floor's cap). No survivor.
M = [
 ('B00', 'control', "def take(self, idx):", "def take(self, idx):", 0),
 ('B01', "the noise one draw behind", "Th, lb = b.step(Th, lb, dt, sq, Z[:, off])", "Th, lb = b.step(Th, lb, dt, sq, Z[:, off - 1])", 0),
 ('B02', "the block's draws in (m, block, 3)", "Z = np.array([g.standard_normal((block, m, 3)) for g in gens])",
  "Z = np.array([g.standard_normal((m, block, 3)).transpose(1, 0, 2) for g in gens])", 0),
 ('B03', "one generator for the batch", "g = np.random.default_rng(); g.bit_generator.state = states[kk]['rng']; gens.append(g)",
  "g = gens[0] if gens else np.random.default_rng(); g.bit_generator.state = states[kk]['rng']; gens.append(g)", 0),
 ('B04', "the noise with dt, not as given", "Th = Th + dt * dT + (sq[None, None, :] if sq.ndim == 1 else sq[:, None, :]) * z",
  "Th = Th + dt * dT + dt * (sq[None, None, :] if sq.ndim == 1 else sq[:, None, :]) * z", 0),
 ('B05', "the per-run amplitude by class", "else sq[:, None, :]) * z", "else sq[:, :, None]) * z", 0),
 ('B06', "log β without the renormalisation", "lb = lb - np.log(np.exp(lb).sum(axis=1, keepdims=True) / m)", "lb = lb - 0.0 * np.log(np.exp(lb).sum(axis=1, keepdims=True) / m)", 0),
 ('B07', "log β by dβ, not dβ/β", "lb = lb + dt * db / b", "lb = lb + dt * db", 0),
 ('B08', "the replicator's mean fitness unweighted", "lam = np.sum(beta * u, axis=1, keepdims=True) / beta.sum(axis=1, keepdims=True)",
  "lam = np.mean(u, axis=1, keepdims=True)", 0),
 ('B09', "the fitness without 1/N", "u = a * pr / self.N", "u = a * pr", 0),
 ('B10', "the correction's floor on W(1 − q) unfloored", "/ np.maximum(W * (1 - q), self.floor * W)", "/ (W * (1 - q))", 0),
 ('B11', "the masses form's q without 1/W", "(np.abs(self.A_hat)[None, :] / amax / W) ** 2", "(np.abs(self.A_hat)[None, :] / amax) ** 2", 0),
 ('B12', "the correction's τ over m²", "return self.fixman * (self.tau_F / m) * (1 + q)", "return self.fixman * (self.tau_F / m ** 2) * (1 + q)", 0),
 ('B13', "the width's chain rule with s²", "g[:, :, 2] = a * ds / s ** 3 * pd2r", "g[:, :, 2] = a * ds / s ** 2 * pd2r", 0),
 ('B14', "the cohesion without the metric", "d2 = np.sum(D ** 2 / self.G[None, None, None, :], axis=3)", "d2 = np.sum(D ** 2, axis=3)", 0),
 ('B15', "the box's dead zone pulls outward", "dTh = dTh - self.box_sigma * (Th - np.clip(", "dTh = dTh + self.box_sigma * (Th - np.clip(", 0),
 ('B16', "take ignores the per-run γ", "for k in ('gamma', 'fixman', 'tau_F', 'floor'):", "for k in ('fixman', 'tau_F', 'floor'):", 0),
 ('B17', "the record before the step at every step", "if record is not None and off % rec_every == 0:", "if record is not None and off % rec_every == 1:", 0),
 ('B18', "no record at the cap", "if record is not None and states[keys[j]]['exit'] is None:", "if record is not None and states[keys[j]]['exit'] is not None:", 0),
 ('B19', "the exit one step late", "states[keys[j]].update(exit=int(k[j]))", "states[keys[j]].update(exit=int(k[j]) + 1)", 0),
 ('B20', "the generator's state not saved at the cap", "states[keys[j]]['rng'] = gens[j].bit_generator.state      # exactly",
  "pass      # exactly", 0),
 ('B21', "no checkpoint at the block boundary", "            if ckpt is not None:\n                _dump(ckpt, states)\n            if nb",
  "            if False:\n                _dump(ckpt, states)\n            if nb", 0),
 ('B22', "a raised cap not reopened", "            st['done'] = False                                                # a cap raised",
  "            pass                                                # a cap raised", 0),
 ('B23', "the per-run amplitudes not dropped with their runs", "Th, lb, k, n, sq, Z = Th[ix], lb[ix], k[ix], n[ix], sq[ix], Z[ix]",
  "Th, lb, k, n, Z = Th[ix], lb[ix], k[ix], n[ix], Z[ix]", 0),
 ('B24', "the plug-in form's q = |a|/a_max", "q = (a / amax) ** 2 if self.A_hat is None", "q = np.abs(a / amax) if self.A_hat is None", 0),
]
