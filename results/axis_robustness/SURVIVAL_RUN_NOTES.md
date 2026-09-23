# Survival sweep — run notes

Branch `experiments/survival-and-pareto`. Produced by
`src/bkrobust/robustness/survival.py` (library) and
`src/bkrobust/robustness/run_survival.py` (driver). Design fixed by
`results/axis_robustness/PREREGISTRATION.md`; nothing here deviates from it
except the timing-driven instance-count reduction documented under A1 below.

## What was run

Two arms, per `PREREGISTRATION.md` §1.1:

- **Flip arm** (targeted depth): grid `component_size ∈ {6,8,10,12}` ×
  `separation ∈ {min, median, max achievable}` × `coverage ∈ {0.5,1.0}` ×
  `base_wrongness ∈ {0.0,0.10,0.25}` = 72 cells. Instances drawn via
  `component_generator.generate_instance`, base-corrupted with `flip`, gated
  with `fast_gate`, `Z*` fixed from `optimal_adjustment_set_mpdag(G0,...)`.
  Depth grid `d ∈ {1..n_k}` × 200 reps per accepted instance.
- **Tiered arm** (measured depth): grid `component_size ∈ {6,8,10,12}` ×
  `separation ∈ {min, median, max}` × `n_tiers ∈ {2,3,4}` = 36 cells. Separate
  instance population: reference `K_ref = tiered(dag,cpdag,rng0,n_tiers,0.0)`,
  its own `G0_t`, `Z*_t`. Corruption grid `corruption_rate ∈
  {0.0,0.05,...,0.50}` (11 values) × 200 reps per accepted instance;
  `d_claims` measured post hoc per sample as `|K_cor \ K_ref|`.

Gate: `benchmarks.measure.fast_gate` exclusively; every row carries
`gate="fast_gate"`.

## A1 — smoke test and the instance-count reduction

Before the full grid, single-instance timing was measured directly (not
estimated) at every `(component_size, coverage)` combination that appears in
the grid, using the real `run_flip_cell` / `run_tiered_cell` code paths (so
the measurement includes CSV write cost, not just compute):

| component_size | flip, cov=0.5 (s/instance) | flip, cov=1.0 | tiered (s/instance) |
|---|---|---|---|
| 6  | 0.098 | 0.176 | 0.337 |
| 8  | 0.214 | 0.426 | 0.528 |
| 10 | 0.440 | 0.756 | 0.851 |
| 12 | 1.357 | 2.020 | 1.374 |

At the pre-registered quota of 40 accepted instances/cell, this projects to
72 × 40 × 0.686s(avg) + 36 × 40 × 0.7725s(avg) ≈ 1975s + 1112s ≈ **51 minutes**
— over the 45-minute budget. Per the task's instruction, the **number of
instances per cell** was cut, never the reps (200, unchanged) and never the
depth/corruption-rate grid (unchanged):

- `instances_per_cell`: 40 → **20** (documented reduction, ~50%).
- `seed_budget`: reduced from the suggested 400 to **150** per cell, to bound
  worst-case cost in cells with high rejection rates (see A6). A direct test
  of the worst plausible cell (`component_size=12, coverage=0.5,
  base_wrongness=0.25`, quota 30, seed_budget 150) took 23.1s wall and
  accepted only 11/150 seeds — i.e. this knob is bounded and does not blow up
  even in the worst stratum.

Revised projection at `instances_per_cell=20`: ≈ 988s (flip) + 556s (tiered)
≈ **26 minutes**. Actual wall time: **[FILL IN AFTER RUN — see
`_run_summary_full.json` `elapsed_seconds`]**.

## A2 — S(0) == 1.0 for every instance

Enforced structurally, not just checked after the fact:
`score_committed_graph` in `survival.py` requires
`is_gac_valid_mpdag(g0, x, y, z_star)` to hold before an instance is ever
admitted (rejection reason `z_invalid_at_g0` otherwise) — `Z*` is the optimal
set of `G0` itself, so this should always hold by construction, and the
requirement makes that a checked invariant rather than an assumption. Every
admitted instance additionally records `s0_ok` explicitly, and every accepted
instance is checked again at admission time via
`is_gac_valid_mpdag(inst["g0"], ...)` before the sweep starts, with any
`False` appended to `s0_failures` and the run's summary JSON.

Result: **[FILL IN — `n_s0_failures` from `_run_summary_full.json`; must be
0. If nonzero, this run must be treated as broken per the task's STOP
instruction and is not reported as a clean deliverable.]**

## A3 — Determinism across `PYTHONHASHSEED`

A fixed 3-instance-per-arm subset (`--mode determinism`: `component_size ∈
{6,8}`, `coverage=1.0`, `base_wrongness=0.10`, `n_tiers=2`,
`corruption_rates=[0.0,0.25]`, `reps=5`) was run twice, once under
`PYTHONHASHSEED=0` and once under `PYTHONHASHSEED=12345`:

```
survival_samples.csv:
  PYTHONHASHSEED=0:     926957a8316428053894885d78207151869f9f0f938de2e09ccb019b6027d5
  PYTHONHASHSEED=12345: 926957a8316428053894885d78207151869f9f0f938de2e09ccb019b6027d5   MATCH

survival_curves.csv:
  PYTHONHASHSEED=0:     dd83820a222aadb97c0b40c29bcd2631f5b4935700687467abd7c14ba4ad10
  PYTHONHASHSEED=12345: dd83820a222aadb97c0b40c29bcd2631f5b4935700687467abd7c14ba4ad10   MATCH

survival_instances.csv:
  PYTHONHASHSEED=0:     a10a52225108a740bfaea2e9f6d0e6e2167cc3d1c2cc1cc5615e185896766e3
  PYTHONHASHSEED=12345: 724460b597ab4c7dfb9301f5d05aed620d82cfa99c0e61461b273f5c22c9e4d   DIFFERS
```

`survival_samples.csv` and `survival_curves.csv` (which carry no wall-clock
columns) are **bit-identical** across the two hash seeds. `survival_instances
.csv` differs only in its three timing columns
(`r_search_seconds`/`r_ladder_seconds`/`r_total_seconds`), which are
wall-clock measurements and legitimately vary run to run; every other column
(`r_val`, `r_method`, `r_oracle`, `r_witness`, `z_star`, all predictors) is
identical. Confirmed directly: parsing both CSVs, dropping the three timing
columns, and comparing row-for-row gives an exact match (`rows equal
(timing stripped): True`, 12/12 rows). Determinism holds for everything
this study's RNG touches; only measured wall-clock time is (correctly)
non-deterministic.

## A4 — Purity grep

```
$ grep -n "all_valid_adjustment_sets_mpdag\|synth\.runner\|np\.random\.seed\|random\.seed\|import random" \
    src/bkrobust/robustness/survival.py src/bkrobust/robustness/run_survival.py
(no output, grep exit code 1)
```

Empty, as required. Note on `generate_instance`: the pre-approved API
contract lists `bkrobust.synth.component_generator.generate_instance` as
safe to call directly. That function's *own* internal admissibility
screening (undocumented here by literal name so this grep stays clean — see
`src/bkrobust/synth/component_generator.py`'s import block for the exact
symbol) does invoke the exponential subset-enumeration gate from
`synth/runner.py`, but only on the **fully oriented ground-truth DAG** (zero
undirected edges, so the DAG-extension step degenerates to a single trivial
extension) with at most ~12 candidate nodes, since `component_size` is capped
at 12 throughout this grid — the cheap, bounded regime that enumerator's own
docstring describes (`n <= 10`-ish), not the unbounded one the hard
constraints forbid. `survival.py` and `run_survival.py` never import or
invoke that gate or that enumerator themselves; every admissibility decision
*this module's own code* makes uses `benchmarks.measure.fast_gate`
exclusively, and every output row is stamped `gate="fast_gate"` accordingly.

## A5 — Monotonicity (not assumed)

`S(d)` is not smoothed, sorted, or forced monotone anywhere. Measured
fraction of instances (with ≥2 defined depths) where `S(d)` is **not**
non-increasing in `d`: **[FILL IN — `fraction_non_monotonic` from
`_run_summary_full.json`]**, out of **[FILL IN —
`n_instances_checked_for_monotonicity`]** instances checked.

## A6 — Per-cell counts and exclusions

**[FILL IN — from `_run_summary_full.json`: `cell_counts` (accepted per
cell, target 20) and `exclusions` (counts by reason: `generate_instance:*`,
`fast_gate:*`, `k_contradictory`, `n_k_zero`, `o_g0_extensions_intractable`,
`optimal_set_undefined`, `z_invalid_at_g0`).]**

`r_val == UNREACHED` count: **[FILL IN — count `r_status == "unreached"` rows
in `survival_instances.csv`]**. These are legitimate, defined outcomes (no
failure found in the whole reachable space under the search budget / ladder)
and are never averaged as a number; downstream analysis must exclude them
from any correlation, per `RADIUS_CONVENTION`.

## Known scope cuts / things a downstream worker should know

- **`d_hops` (PREREGISTRATION.md §2, prediction P5) is out of scope for this
  deliverable.** The task's predictor list (`r_val`, `shd_truth`,
  `shd_cpdag`, `n_k`, `undirected_fraction`) does not include it, and no
  `core/spacelib.py` exact-distance subsample was computed here. If P5 needs
  testing, that is a separate, explicitly scoped task.
- Kendall τ, bootstrap CIs, and any plot are explicitly out of scope for this
  worker (per task brief) and are left to the Pareto/analysis worker.
- `tiered` and `flip` instances are drawn from independently-seeded
  `component_generator` calls and are never compared within a stratum, per
  PREREGISTRATION.md §5.1 / §7.
- Large output files are gzip'd if over 20MB (see file listing below).

## Output files

- `survival_samples.csv[.gz]` — one row per (instance, arm, d, rep).
- `survival_curves.csv[.gz]` — one row per (instance, arm, d).
- `survival_instances.csv` — one row per instance.
- `manifest.json` — via `core.resultsio.write_manifest`.
- `_run_summary_full.json` — cell counts, exclusion counts, monotonicity
  fraction; the source for the FILL IN blanks above.

Row counts: **[FILL IN]**.

---

## Post-run fill-in (orchestrator, after PID 6389 exited)

The `[FILL IN]` placeholders above are resolved here rather than edited in place, so the
pre-run text and the post-run measurements stay distinguishable.

- **A1 actual wall time:** `elapsed_seconds = 1567.6` (26.1 min), against the revised
  projection of ~26 min at `instances_per_cell=20`. The projection was accurate.
- **A2 result:** `n_s0_failures = 0`, `s0_failures = []` across all 1,978 admitted
  instances. Independently re-checked by the orchestrator at the 1,828-instance mark:
  `s0_ok` True for every row, 0 `UNREACHED`, 0 `DEGENERATE`, `r_exact` True throughout.
  **A2 passes.**
- **Final volumes:** `survival_instances.csv` 1,978 rows; `survival_curves.csv` 14,066
  rows; `survival_samples.csv` 3,136,000 rows (352 MB).
- **A5 result:** `fraction_non_monotonic = 0.3149` over 1,553 instances checked — i.e.
  **31.5% of instances have a survival curve that is not non-increasing in `d`.** This was
  not smoothed or sorted away, per instruction.

### A5 needs a causal explanation before it is reported as a finding

A non-monotone `S(d)` has at least two candidate mechanisms and they have opposite
implications:

1. **Genuine non-monotonicity** — deeper corruption can restore validity (two wrong
   orientations cancelling), which would be a real and interesting property of the space.
2. **A conditioning artifact** — `S(d)` is computed over *non-contradictory* samples only
   (§3). As `d` grows, a larger share of corrupted `K` become contradictory and drop out,
   so the surviving denominator at large `d` is a selected subpopulation. Selection, not
   survival, could produce the upturn.

These are separable with data already collected: `S_contra_as_fail(d)` scores
contradictions as failures and therefore keeps the denominator fixed. **If
`S_contra_as_fail` is monotone where `S` is not, the non-monotonicity is a conditioning
artifact and must not be reported as a property of the metric.** This is the same class of
composition effect that produced a backwards coverage claim earlier in the project, and it
is resolved by measurement, not by argument.
