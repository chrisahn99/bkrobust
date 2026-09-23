# Null-stratum isolation: flip, coverage=0.5, base_wrongness=0.25

Follow-up to `PREREGISTRATION.md` Appendix F. Session 7's flip arm returned
`tau_b(AUC_frac, r_val) = +0.031, CI [-0.114, +0.169], n=220` in this one
stratum, versus `+0.346` to `+0.714` (CIs excluding zero) in the other five
flip strata. Two explanations were on the table:

- **(a) Genuine loss of predictive power** — under base wrongness 0.25,
  GAC validity collapses so fast that `r_val` stops discriminating.
- **(b) Attenuation by measurement noise** — `AUC_frac` at N=200 draws/depth
  is noise-dominated in this specific cell (high contradiction rate ->
  few usable samples), and Appendix F already found
  `tau(r_val, usable depths) = -0.100` here, meaning high-`r_val` instances
  rest on *less* data in this cell specifically, which biases the *observed*
  AUC_frac in a way that could manufacture a spurious null.

Design: same instances (same `ComponentSpec` grid, same seed sequence, same
base-wrongness corruption, same `fast_gate`, same `Z*`) as session 7's cell,
sampled at N=1000 draws/depth instead of 200, then three `tau_b` estimates
computed on the identical instance set.

---

## VERDICT: (a) — genuine null. Not attenuation.

| check | endpoint | N draws/depth | n | tau_b | CI (2.5%, 97.5%) |
|---|---|---|---|---|---|
| (1) reproduction | AUC_frac, first 200 of 1000 draws | 200 | 220 | **+0.0314** | [-0.1144, +0.1686] |
| (2) full | AUC_frac, full 1000 draws | 1000 | 220 | **+0.0324** | [-0.1117, +0.1685] |
| (3) usable | AUC_frac_usable (n_eval>=30), 1000 draws | 1000 | 220 | **+0.0322** | [-0.1118, +0.1681] |

Five times the draws per depth moved `tau_b` by **+0.0010** (0.0314 ->
0.0324) — not "materially," by any reading. The CI barely narrowed (width
0.2831 -> 0.2803, a 1% change) and still covers zero by a wide margin in
both directions. Per the task's own falsification rule ("if tau stays near
zero with five times the data and tight CIs, explanation (a) holds"): it
did, and it does.

The mechanism check makes this unambiguous rather than merely "consistent
with (a)": mean per-instance standard error of `AUC_frac` fell from 0.00804
(N=200) to 0.00359 (N=1000), a ratio of **2.240**, matching the theoretical
`sqrt(5) = 2.236` almost exactly — the extra draws did exactly what more
draws should do, cut sampling noise by the expected factor. If explanation
(b) were correct, that noise reduction should have let a real
`r_val`-correlated signal emerge from underneath it, moving individual
per-instance `AUC_frac` values in a direction that raises `|tau_b|`. It
didn't: the mean shift in `AUC_frac` from N=200 to N=1000 was **+0.00084**
(median **0.0** exactly, std 0.0135) — the extra precision tightened each
instance's estimate around essentially the *same* value, not a
systematically different one. There is no signal to uncover. `r_val` really
does stop discriminating once base wrongness reaches 0.25 at coverage 0.5;
this is a confirmed null, and per the task brief it is reported as such —
nothing here was tuned to move it off zero.

(CI width is governed by instance count, not draws/depth, exactly as the
task anticipated — see the note on the optional second axis at the end.)

---

## K1 — timing

One instance (`c06_s01_cov050_seed00000000_bw025`, `n_k=4`) at N=1000: **0.570s**
for 4000 draws (0.1426 ms/draw), timed with `time.perf_counter()` before the
full cell ran. Projected from `sum(n_k)=926` across 220 instances at that
per-draw rate: **~132s** for the pure sampling compute.

Measured wall time for the full cell's sampling phase (`[4/5]` in the run
log): **353.9s** (~5.9 min), plus ~4s for instance construction/screening and
~4s for tau/manifest/gzip at the end — **357.9s total**, comfortably under
the ~10 min budget. The ~2.7x gap between the 132s compute-only projection
and the 354s measured time is CSV write overhead (three interleaved
`ResultWriter`s, ~926,000 rows to `null_samples.csv` at 14 real columns each,
not the 2-column synthetic benchmark used to sanity-check `ResultWriter`
throughput beforehand) plus per-instance curve-building/SE/usable-depth
bookkeeping done twice (N=200 and N=1000) — not a surprise, and still well
inside budget, so the depth grid was not shrunk.

## K2 — S(0) = 1.0

All 220 accepted instances pass: `is_gac_valid_mpdag(inst["g0"], x, y,
z_star)` is `True` for every one (`s0_failures: []` in `null_manifest.json`).
Z* is `G0`'s own optimal adjustment set by construction, so this is an
integrity check on the rebuild, not a discovery — it passed, so the run
proceeded.

## K3 — N=200 reproduction check

Not merely approximate — **exact**. All 220 `instance_id`s in
`null_instances.csv` match session 7's `survival_instances.csv` (flip,
coverage=0.5, base_wrongness=0.25) one-for-one (220/220 overlap, 0 only in
either side). For every matched instance, `n_k` and `r_val` are identical,
and `AUC_frac` computed here from the first 200 of 1000 draws equals
session 7's original `AUC_frac` to the exact floating-point value (max
absolute difference across all 220 instances: **0.0**). The seed derivation
(`derived_seed` -> `np.random.default_rng`, keyed on `(instance_id, "flip",
d, rep)`) reproduces session 7's draws bit-for-bit, as designed.
`tau_b` on this exact reproduction: **+0.0314**, CI [-0.1144, +0.1686],
n=220 — matches session 7's reported +0.031, CI [-0.114, +0.169], n=220
(the sub-0.0004 display-precision difference is rounding in the original
report, not a discrepancy in the rebuild).

## K4 — determinism

3-instance subset (first cell, `component_size=6, separation=1`, seeds
0-2), 50 draws/depth, hashed over deterministic columns only (`d`, `rep`,
`seed`, `status`, `survived`, `symdiff_proxy_not_distance`, `n_k`, `z_star`,
`r_val`, `r_status`, `shd_truth` — wall-clock columns `r_search_seconds`,
`r_ladder_seconds`, `r_total_seconds`, `wall_until_timeout_s`,
`r_stat_*` excluded, noted here rather than silently dropped):

```
PYTHONHASHSEED=0     sha256=9d1757b5dc6bd80dca59f7cdfc15ca6ad412212bdee47cbbe8b96ca85144c8aa
PYTHONHASHSEED=12345  sha256=9d1757b5dc6bd80dca59f7cdfc15ca6ad412212bdee47cbbe8b96ca85144c8aa
```

Identical (503 lines hashed both times). Reproduce with
`PYTHONHASHSEED=<seed> PYTHONPATH=src python3 -m bkrobust.robustness.run_nullcell --mode determinism`.

## K5 — purity grep

```
$ grep -nE "all_valid_adjustment_sets_mpdag|synth\.runner|np\.random\.seed|random\.seed|import random" \
      src/bkrobust/robustness/nullcell.py src/bkrobust/robustness/run_nullcell.py
(no output, exit 1)
```

## K6 — exclusions, counted

Instance-population build (`build_flip_instance` rejections while filling
the 20-per-cell quota, 12 cells, seed budget 150 -> 1800 seeds tried, 220
accepted):

- `k_contradictory` (base-wrongness-corrupted `K_b` gave a contradictory
  `G0`, before any depth sampling): **620**
- `optimal_set_undefined`: **2**
- no `fast_gate` rejections, no `generate_instance` rejections, no
  `n_k_zero`, no `z_invalid_at_g0`, no intractable-extension rejections in
  this stratum

Per-sample outcomes during the depth sweep (not instance exclusions —
`corrupted_k_contradictory` is its own status per the sentinel discipline,
never scored as survival-0): pooled contradiction rate by depth in this
cell (N=1000; N=200 in parentheses, from `null_manifest.json`):

| d | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|
| contradiction rate (N=1000) | 0.519 | 0.630 | 0.806 | 0.958 | 0.979 | 0.978 | 1.000 | 1.000 |
| contradiction rate (N=200) | 0.518 | 0.632 | 0.805 | 0.957 | 0.980 | 0.977 | 1.000 | 1.000 |

(Consistent with N and with the task's headline figures for base_wrongness
0.25 — the pooling here is coverage=0.5 only, session 7's `flip
contradiction rate 0.563 at d=1 and 0.858 at d=3` pools something slightly
different, e.g. across coverage, hence the small gap; the shape — steep rise
by d=3-4, near-total contradiction by d=7-8 — is the same phenomenon.) As
expected, the *rate* is essentially invariant to N (0.519 vs 0.518 at d=1);
what N changes is how many non-contradictory draws are left to estimate
`S(d)` from.

`r_val`: `UNREACHED` count = **0** (min 1, max 5, distribution `{1: 41, 2:
63, 3: 84, 4: 31, 5: 1}` — identical to session 7's). Timeouts: **0**
(`r_status == "ok"` for all 220; max `r_total_seconds` well under
`time_limit_s=60.0`, no `wall_until_timeout_s` populated for any instance).

## K7 — no pre-existing file touched

```
$ git status --short
?? results/axis_robustness/_xarm_sweep_summary_0.json      (concurrent worker, not this run)
?? results/axis_robustness/null_curves.csv
?? results/axis_robustness/null_instances.csv
?? results/axis_robustness/null_manifest.json
?? results/axis_robustness/null_samples.csv.gz
?? results/axis_robustness/null_tau_comparison.csv
?? results/axis_robustness/xarm_bins.csv                    (concurrent worker)
?? results/axis_robustness/xarm_instances.csv                (concurrent worker)
?? results/axis_robustness/xarm_samples.csv                  (concurrent worker)
?? src/bkrobust/robustness/crossarm.py                       (concurrent worker)
?? src/bkrobust/robustness/nullcell.py
?? src/bkrobust/robustness/run_crossarm.py                   (concurrent worker)
?? src/bkrobust/robustness/run_nullcell.py
```

All rows are untracked new files (`??`); nothing shows as modified (` M`) or
deleted (` D`). Everything under `null_*` and `nullcell.py`/`run_nullcell.py`
is this run's output; the `xarm_*`/`crossarm.py` rows belong to the
concurrent worker mentioned in the task brief and were not touched.

**One bug caught and fixed during this run, noted for the record:** the
shared `bkrobust.core.resultsio.write_manifest(directory, ...)` helper
always writes to `<directory>/manifest.json`. An early version of
`run_nullcell.py` called it directly against `results/axis_robustness/` and
then renamed the result to `null_manifest.json` — which briefly overwrote
session 7's committed `manifest.json` in the working tree (it showed as ` D`
in `git status`) before the rename moved it out from under that name. This
was caught by inspecting `git status` after the run (not by any check in the
code) and fixed two ways: `git checkout -- results/axis_robustness/manifest.json`
restored the original committed content (confirmed via `git status --short`
showing a clean file afterward), and `run_nullcell.py` was changed to call
`write_manifest` against a throwaway `tempfile.TemporaryDirectory()` and
copy only the resulting bytes into `null_manifest.json`, so the real
`manifest.json` path is never opened for writing by this run again. No
sweep needed to be rerun — `null_manifest.json`'s content was already
correct; only the incidental clobber-then-rename of the unrelated committed
file was the problem, and it is now structurally impossible in this code.

---

## Mechanism, made visible

**Usable depths per instance (n_eval >= 30) at N=200 vs N=1000:**

| n_usable_depths | N=200 (count) | N=1000 (count) |
|---|---|---|
| 1 | 27 | 4 |
| 2 | 102 | 55 |
| 3 | 85 | 107 |
| 4 | 4 | 42 |
| 5 | 2 | 9 |
| 6 | 0 | 3 |

Mean usable depths per instance: 2.33 (N=200) -> 3.03 (N=1000). Every
instance has at least one usable depth at both N (0 instances with zero
usable depths either way), so `AUC_frac_usable` is defined for all 220 at
N=1000 (no NaNs, versus session 7's pooled 1.2% NaN rate across all
strata).

**Standard error of per-instance `AUC_frac`** (analytic: per-depth binomial
proportion variance `p(1-p)/n_eval`, propagated through the exact
fractional-grid weights `survival.auc_frac` uses — see
`nullcell.auc_frac_se`):

| | N=200 | N=1000 | ratio |
|---|---|---|---|
| mean SE | 0.00804 | 0.00359 | 2.240 |
| median SE | 0.0 | 0.0 | -- |

(Median SE is exactly 0 at both N because a large share of instances have
`AUC_frac` pinned at exactly 0 or 1 on every depth that receives grid
weight — `p(1-p)=0` — so their SE is genuinely zero, not a computation
artifact; this is itself part of why the endpoint is hard to correlate with
anything here.) The mean-SE ratio of 2.240 matches `sqrt(5) = 2.236` almost
exactly, confirming the extra draws behaved exactly as pure noise reduction
should — and it still didn't move `tau_b`.

**Per-instance `AUC_frac` shift, N=200 -> N=1000:** mean **+0.00084**,
median **0.0**, std 0.0135, range [-0.063, +0.085] across 220 instances —
consistent with pure sampling noise around an unchanged central value, not
a systematic correction in any direction.

---

## Optional second axis (extending instance count) — not pursued

The task allows extending the instance count in this cell beyond 220 to
tighten the bootstrap CI, as a strictly separate result (CI width is driven
by instance count; attenuation is driven by draws/depth; conflating them
defeats the point of this run). It was not pursued: the primary result
above is already unambiguous (a ~1% CI-width response to a 5x draw increase
is itself the evidence that instance count, not draw count, is the CI's
binding constraint), and running it would not change the (a)-vs-(b) verdict,
only the null's confidence interval — a different, legitimately separate
question the task explicitly says to report "separately or not at all."
Reporting it here would risk exactly the conflation the task warns against,
so it is left undone and flagged rather than attempted.

---

## Outputs

All under `results/axis_robustness/`, all new (`null_` prefix), all inputs
regenerated read-only via `bkrobust.robustness.survival` (never running
`run_survival.py` or any other frozen sweep):

| file | rows (incl. header) | notes |
|---|---|---|
| `null_samples.csv.gz` | 926,001 | every draw, both N=200-subset and N=1000 (same rows; N=200 is `rep < 200`) |
| `null_curves.csv` | 1,853 | per-instance, per-depth, per-N (200 and 1000) curve stats: `n_eval`, `n_contradictory`, `contradiction_rate`, `S`, `S_contra_as_fail` |
| `null_instances.csv` | 221 | 220 instances, all three AUC endpoints, both SEs, both usable-depth counts, `r_val` and its full provenance |
| `null_tau_comparison.csv` | 4 | the three `tau_b` rows above, with `n`, CI, `p_value`, `n_boot_nan` |
| `null_manifest.json` | -- | grid, environment, git SHA, exclusions, contradiction-rate-by-depth, elapsed time |
| `NULLCELL_RUN_NOTES.md` | -- | this file |

`null_samples.csv.gz`: 110.0MB uncompressed -> 14.1MB gzipped (over the 20MB
threshold, compressed per the output spec).

## Anything surprising

- The exact bit-for-bit reproduction (K3) was stronger than expected —
  worth calling out since it means the N=200 vs N=1000 comparison is a
  clean paired design (each N=1000 curve literally contains the N=200 curve
  as its first 200 reps), not a resample from a similar-but-different
  population.
- The SE ratio landing at 2.240 against a theoretical 2.236 is closer than
  the analytic-SE approximation (independent-depths, plug-in `p`) had any
  right to expect; it is a useful confirmation that the per-depth
  contradiction/survival draws really are behaving as independent Bernoulli
  processes at the rates reported, not something more structured.
- The `write_manifest` footgun documented under K7: any future `null_*` (or
  similarly prefixed) driver that calls `bkrobust.core.resultsio.
  write_manifest` directly against `results/axis_robustness/` will silently
  clobber the shared `manifest.json` unless routed through a temp directory
  first, exactly as this file's fix now does. Worth fixing at the
  `resultsio.py` level (e.g. an optional filename parameter) if another
  worker hits the same trap — flagged here rather than changed there, since
  that file is shared infrastructure outside this task's file list.
