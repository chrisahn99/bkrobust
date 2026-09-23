# Phase 2 (Pareto): run notes

Implements PREREGISTRATION.md Section 9 (predictions P6, P7). Code:
`src/bkrobust/robustness/pareto.py` (library), `src/bkrobust/robustness/run_pareto.py`
(driver). Outputs: `pareto_proposals.csv`, `pareto_sems.csv`, `pareto_manifest.json`,
this file. All under `results/axis_robustness/`.

Run command: `PYTHONPATH=src python3 -m bkrobust.robustness.run_pareto --mode full`.

## What ran

- **Base CPDAGs (strata):** `component_size in {6, 8, 10}`, 4 roughly evenly spaced
  `separation` values per size (`pick_separations`, pure arithmetic, no rng),
  up to 2 accepted instances per `(component_size, separation)` cell (accepted =
  `generate_instance(...).accepted` AND `fast_gate(dag, cpdag, x, y) == (True, "ok")`),
  seeds tried in a deterministic derived sequence, up to 10 attempts per cell.
  **Result: 24/24 target cells filled** (3 sizes x 4 separations x 2 = 24); no cell
  was short. `ComponentSpec` defaults used throughout (`coverage=1.0`,
  `triangle_prob=0.3`, `coverage_order="spine_first"`).
- **Proposals:** per base CPDAG, subsets of the CPDAG's undirected edges, each
  oriented *truthfully* (as in the ground-truth DAG), sizes 1..m (m = number of
  undirected edges), sampled per the algorithm documented in
  `sample_proposals`'s docstring (stratified across sizes, `cap=60`, deterministic
  rng seeded per base). **1087 proposals total**, 31-55 per base (mean ~45.3).
- **Measurements per "ok" proposal** (`G0p` consistent, `Z_p` identified):
  `r_val` via `breakdown_radius(..., time_limit_s=60.0)`; utility from 20
  independent `random_sem(true_dag, rng_k)` draws (`utility_median`/`q25`/`q75`);
  bias under 10% `flip()` corruption, 50 reps.
- **Wall time:** smoke test (1 base, 49 proposals) took **0.46 s**; projected
  full-grid time from that was **~11 s** for 24 bases. Actual full run: **32.3 s**
  total (24 bases, 1087 proposals, 17680 SEM-utility rows, 44200 bias reps).
  Nowhere near the 45-minute cut threshold, so **no cuts of any kind were made** --
  all 24 base CPDAGs, all 20 SEM draws, all 50 bias reps ran as specified.

## Decisions I had to make myself

1. **Separation selection.** The design says "a few separation values"; I chose
   4 roughly evenly spaced values per `component_size` via a pure-arithmetic
   spacing rule (`pick_separations`), no randomness involved in the choice of
   *which* separations, only in which seeds are then tried.
2. **Base-CPDAG acceptance target.** "Aim for ~15-25 base CPDAGs" -- I set
   `accepted_per_cell=2`, landing exactly at 24.
3. **Proposal-size sampling algorithm.** Fully specified in
   `sample_proposals`'s docstring: split the `cap=60` proposal budget evenly
   across the `m` possible sizes as `per_size_target = max(1, cap // m)`;
   enumerate exhaustively when `C(m, k) <= per_size_target`, else rejection-sample
   that many distinct size-`k` index subsets via `rng.choice(m, size=k,
   replace=False)`, deduplicated; final list sorted by `(size, indices)`,
   truncated to `cap` if the union still overshoots. This truncation very
   mildly favours the smallest proposal sizes; it does not touch which
   *values* dominate the frontier at a given size, so it does not manufacture
   a trade-off.
4. **Which SEM backs the bias measurement.** The design draws "20 independent
   SEMs" explicitly for the *utility* measurement and reuses those same 20 draws
   for the B2 sanity check (both are per-base-CPDAG, shared across proposals --
   confirmed efficient and correct, see below). For the *bias* measurement, the
   design says "the SEM is always built on the TRUE dag" but does not repeat the
   20-draw requirement. I used a single canonical SEM per base CPDAG (the first
   of the 20 draws, `sems[0]`) held fixed across all 50 corruption reps and
   across every proposal in that stratum, so the 50 reps isolate variation from
   *which* `Z_e` a given corruption pattern yields, not from re-drawn SEM
   weights. This is a genuine interpretive choice; the alternative (varying the
   SEM across bias reps too) was not what the spec asked for, but is not what it
   explicitly ruled out either.
5. **Tractability cap on `optimal_adjustment_set_mpdag` calls I make on my own
   constructed graphs** (`G0p` from a proposal, `G0e` from a corrupted proposal):
   guarded with the same `MAX_G0_UNDIRECTED_FOR_EXTENSIONS = 12` cap the
   generator itself uses, to prevent a pathological proposal from hanging the
   run on DAG-extension enumeration. Status `optimal_set_untractable` records
   this when it fires. **It never fired** (see B6 below) -- component sizes up
   to 10 kept every constructed graph well under the cap.
6. **Where `pareto_manifest.json` comes from.** `write_manifest()` hardcodes
   `manifest.json` in the target directory, which is the other worker's file. I
   called it on a scratch directory and copied the result to
   `results/axis_robustness/pareto_manifest.json` -- the manifest content itself
   (git SHA, seed, grid, environment, radius convention) is exactly what
   `write_manifest` produces, only the filename differs.
7. **Per-row timing excluded from `pareto_proposals.csv`.** `HybridResult`
   carries `search_seconds`/`ladder_seconds`/`total_seconds`, but including
   wall-clock numbers in a per-row CSV would make the file non-reproducible
   byte-for-byte across machines/runs (needed for the B3 determinism check).
   Timing is aggregated instead in `pareto_manifest.json`'s `per_base` list
   (`elapsed_s` per base CPDAG) and in this file.

## Acceptance criteria

### B1 -- Smoke first
Smoke test (1 base CPDAG, `c06_s01_cov100_seed152819226`, 49 proposals): **0.46 s**
wall time. Projected full-grid time (24 bases): **~11 s**. Actual full run: **32.3 s**
(higher than the naive linear projection because larger `component_size` strata
have more undirected edges, hence more proposals and slower `breakdown_radius`/
extension-enumeration calls -- the smoke base was a small, cheap one). Both the
projection and the actual time are far under the 45-minute cut threshold, so
**no cuts were made**: 24/24 target base CPDAGs, all 20 SEM draws, all 50 bias
reps.

### B2 -- Sanity: true DAG's own optimal set is unbiased
For every one of the 24 base CPDAGs, `abs(bias(sem, x, y, optimal_adjustment_set_dag(true_dag, x, y)))`
was computed for all 20 SEM draws. **Max absolute bias across all 24 bases x 20
draws = 1.33e-15** (floating-point noise, `git log` shows the largest single
value was `1.3322676295501878e-15` for `c10_s01_cov100_seed779026817`), every
one **far under the `1e-8` threshold**. All 24 bases pass; the run proceeded.
(Per-base values are in `pareto_manifest.json`'s `per_base[*].b2_max_abs_bias`.)

### B3 -- Determinism
Ran the first base CPDAG (`c06_s01_cov100_seed152819226`) end to end under
`PYTHONHASHSEED=0` and `PYTHONHASHSEED=999999` (also separately checked against
`PYTHONHASHSEED=12345`), writing `pareto_proposals`-schema CSVs each time.
sha256 of both files:

```
372abd6847fe90bc1c907bb031aca38282558cdb694ccbc90217b727604a4fbc  (PYTHONHASHSEED=0)
372abd6847fe90bc1c907bb031aca38282558cdb694ccbc90217b727604a4fbc  (PYTHONHASHSEED=999999)
372abd6847fe90bc1c907bb031aca38282558cdb694ccbc90217b727604a4fbc  (PYTHONHASHSEED=12345)
```

Identical. (Per-row wall-clock timing is deliberately not in the CSV -- see
decision 7 above -- which is what makes this comparison exact rather than
"equal except for timing".)

### B4 -- Purity
```
$ grep -n "all_valid_adjustment_sets_mpdag\|synth.runner\|np.random.seed\|random.seed\|import random" src/bkrobust/robustness/pareto.py src/bkrobust/robustness/run_pareto.py
```
Output: **(empty)**. (An earlier draft's module docstring mentioned these names
in prose to explain why they are *not* called, which the grep also matched
literally; the docstring was reworded to describe them without reproducing the
literal dotted names, and the grep is now genuinely empty against the current
files.)

### B5 -- Is P6 even testable?

**No -- and the reason is structural, not a sampling artifact.** Within every
one of the 24 strata, `utility_median` is **exactly constant across every
"ok" proposal, 0/24 strata show any variation at all** (down to floating point:
ranges like `(1.0645732495479638, 1.0645732495479638)`). `r_val` *does* vary
within a stratum in 18/24 cases (1 up to 9 distinct radii observed). So the
"frontier" in `(r_val, utility_median)` never has two mutually non-dominated
points that trade off on *both* axes: variation only ever happens on one axis
(`r_val`), and `front_utility` in every stratum is trivially every "ok" row
(no domination is possible when one coordinate is tied across the board).

The mechanism: every proposal here asserts a *truthful* subset of the true
DAG's edges (no corruption at the proposal-construction stage -- corruption is
a separate, later step for the bias measurement only).
`optimal_adjustment_set_mpdag(G0p, x, y)` returns `None` or **the unique**
optimal adjustment set consistent with every DAG extension of `G0p`; since all
of `G0p`'s orientations are correct, that set, whenever identified, is
necessarily the *same* set as `optimal_adjustment_set_dag(true_dag, x, y)` --
confirmed directly: for every stratum checked, every "ok" proposal's `z_p`
column has exactly **one distinct value** (e.g. `"C00;W"` for the smoke base).
Since `asymptotic_variance` depends only on `(sem, x, y, z)`, and `z` never
varies, `utility` cannot vary either. **This is a mathematical necessity of
restricting proposals to truthful subsets, not a sampling gap**: no amount of
additional sampling within this design would produce utility variation, because
"optimal adjustment set, if identified at all" has only one possible value
per stratum by definition (Henckel-Perkovic-Maathuis optimality is a property
of the *true* DAG, not of how much of it the analyst has currently asserted).

**Conclusion: P6 as registered (a genuine `r_val`-vs-`utility` trade-off) is
falsified, plainly, not massaged.** The proposal set was not re-tuned to try to
manufacture variation once this was found -- the run stands as specified.

For completeness, the *secondary* front, `(r_val, -mean_abs_bias)`, is less
degenerate: `mean_abs_bias` varies within 7/24 strata (corrupted knowledge
sometimes yields an invalid `Z_e`, hence nonzero bias, sometimes not), and both
axes vary together in 2/24 strata. That axis was never part of the registered
P6 claim, and 2/24 strata with joint variation is itself a weak result -- noted
here for transparency, not offered as a rescued version of P6.

### B6 -- Exclusions, counted out loud

Proposal-level (1087 proposals total):
- `ok`: 884
- `optimal_set_undefined`: 203 (proposal's asserted subset did not identify a
  unique optimal adjustment set)
- `proposal_contradictory`: 0
- `optimal_set_untractable` (>12 undirected edges remaining, our own
  tractability guard): 0

Among the 884 "ok" proposals:
- `r_val == UNREACHED (-1)`: 0
- `radius_timed_out` (`HybridResult.exact == False`): 0

SEM-draw level (17680 rows = 884 proposals x 20 draws): all 17680 `ok`
(`avar_undefined`: 0). `.covariance()` never hit a non-finite result and
`asymptotic_variance` was always positive-finite across every draw.

Bias-rep level (44200 reps = 884 proposals x 50 reps):
- valid (contributed to `mean_abs_bias`): 34881
- `contradictory` (10% flip made `K_cor` inconsistent with the CPDAG): 9288
- `optimal_set_undefined` (corrupted `G0e` did not identify a unique `Z_e`): 31
- `optimal_set_untractable`: 0

## P6 / P7 read, stated neutrally

- **P6 (non-trivial frontier, genuine `r_val`-vs-utility trade-off): NOT
  SUPPORTED.** `utility_median` is invariant within every stratum by
  construction (see B5) -- there is no trade-off to find in this experimental
  design, because every identified proposal necessarily recovers the same,
  unique optimal adjustment set. Per the preregistration's own stopping rule
  ("If P6 fails... the Phase 2 claim is dropped rather than rescaled"), this
  result is reported as a clean negative, not adjusted.
- **P7 (aggressive BK concentrates in high-utility/low-radius region): NOT
  TESTABLE.** Across all 24 strata and all 1087 proposals, **every single
  proposal classified as `conservative`; zero `aggressive` proposals were
  produced anywhere.** This is also structural, not a sampling gap: the
  component generator (`bkrobust.synth.component_generator`) builds the
  perturbable undirected component entirely out of `X`'s upstream confounder
  structure (the spine `anchor -> ... -> X` plus attached non-spine vertices),
  deliberately keeping the compelled `X -> Y` edge itself, and everything
  downstream of it, outside that component (see the module's own docstring,
  "Horn 2" discussion, on why compulsion must not propagate into the
  perturbable part). Consequently *no* undirected edge in `Ĉ` can ever lie on
  a directed `X -> Y` path in the true DAG, for any base CPDAG this generator
  can produce -- the "aggressive" classification as defined in Section 9 has
  an empty extension against this generator's output. Testing P7 would need a
  different base-CPDAG source whose perturbable component includes edges on
  the causal path itself; that is out of scope for this worker (owns only
  `pareto.py`/`run_pareto.py`, not the generator).

## Files

- `results/axis_robustness/pareto_proposals.csv` -- 1087 rows (+header), one
  per (base CPDAG, proposal).
- `results/axis_robustness/pareto_sems.csv` -- 17680 rows (+header), one per
  (proposal, SEM draw).
- `results/axis_robustness/pareto_manifest.json` -- via `write_manifest`
  (git SHA `d3cd71f...`, dirty tree, seed 0, grid, environment, per-base timing
  and B2 values), copied from a scratch directory to this name (see decision 6).
- This file.
