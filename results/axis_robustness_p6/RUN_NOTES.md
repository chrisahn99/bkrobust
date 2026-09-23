# Point 6 survival rerun — run notes

Branch `revision/experiments_points67`. Produced by
`src/bkrobust/robustness/survival_p6.py` (library) and
`src/bkrobust/robustness/run_survival_p6.py` (driver), both new files that
import from -- and never modify -- `src/bkrobust/robustness/survival.py` /
`run_survival.py`, per this repo's standing precedent (commit 748f2e1) that a
module which produced a committed result is never edited after the fact.

This rerun answers reviewer Point 6: the committed
`results/axis_robustness/` survival table was run at `N=200` draws/grid-point
and `instances_per_cell=20`, with `table_tau_comparisons.md` explicitly
flagging two strata as unresolved/unstable at that budget (footnote
`[^unstable]`):

- flip, coverage=0.5, base_wrongness=0.25 -- already partly addressed by an
  ad hoc `N=1000` re-run (`results/axis_robustness/null_*` -- "Appendix H"),
  which found a **confirmed null** (`τ_b(r_val) = +0.032`, CI includes 0).
- flip, coverage=1.0, base_wrongness=0.00 -- **no** `N=1000` re-run existed;
  flagged **unresolved**.

This rerun does not itself compute τ_b or any paired bootstrap (Track D
owns that analysis); it produces the underlying `N=1000`,
`instances_per_cell=40` sample/curve/instance tables, with corruption
denominators and new baseline predictors, across the **full** pre-registered
grid (not just the two flagged strata), so Track D's analysis is not
confined to a hand-picked subset.

## What was run

Same two arms and same grid as `results/axis_robustness/SURVIVAL_RUN_NOTES.md`,
at the originally pre-registered budget the committed run had to cut for
time:

- `reps`: 200 -> **1000**
- `instances_per_cell`: 20 -> **40** (the original pre-registration quota)
- `seed_budget`: 150 -> **600** (scaled up so the low-acceptance strata --
  see "Cells below quota" below -- get a fair chance at 40/cell, not just the
  easy ones)
- Grid unchanged: `component_size ∈ {6,8,10,12}` × `separation ∈ {min,
  median, max achievable}` × `coverage ∈ {0.5,1.0}` × `base_wrongness ∈
  {0.0,0.10,0.25}` (flip, 72 cells) and `component_size` × `separation` ×
  `n_tiers ∈ {2,3,4}` (tiered, 36 cells); `corruption_rate` grid unchanged
  (11 points, 0.0 to 0.5 step 0.05).
- Gate: `benchmarks.measure.fast_gate` exclusively, as before; every row
  carries `gate="fast_gate"`.

### Parallelisation (disjoint slices, merged deterministically)

Run as 9 disjoint OS processes (split by `component_size`, and by
`base_wrongness` within `component_size=12` to balance load -- that stratum
is the most expensive per instance), **at most 5 concurrent** as instructed
(the machine has 10 cores shared with three other tracks):

| shard | arm | component_size | other restriction |
|---|---|---|---|
| `flip_c06` | flip | 6 | -- |
| `flip_c08` | flip | 8 | -- |
| `flip_c10` | flip | 10 | -- |
| `flip_c12a` | flip | 12 | `base_wrongness ∈ {0.0, 0.10}` |
| `flip_c12b` | flip | 12 | `base_wrongness = 0.25` |
| `tiered_c06` | tiered | 6 | -- |
| `tiered_c08` | tiered | 8 | -- |
| `tiered_c10` | tiered | 10 | -- |
| `tiered_c12` | tiered | 12 | -- |

Each shard writes its own `survival_samples.csv` / `survival_curves.csv` /
`survival_instances.csv` / `manifest.json` / `_run_summary_<shard>.json`
under `results/axis_robustness_p6/shards/<shard>/`, and logs verbosely
(every cell, wall time, running row totals) to
`results/axis_robustness_p6/_logs/<shard>.log`. `run_survival_p6.py --merge`
then concatenates all nine shards' CSVs and sorts them deterministically
(`instance_id`, then `d`/`corruption_rate`, then `rep`) into the top-level
`results/axis_robustness_p6/{survival_samples.csv.gz,survival_curves.csv,
survival_instances.csv}`, and writes a combined `manifest.json` recording
every shard's own grid, seeds and elapsed time.

**STATUS: [FILL IN AFTER MERGE]** -- wall time, cell counts, and the
sections below are completed once all nine shards finish and are merged.
See `_logs/` for the live per-shard audit trail in the meantime.

## Task 4 — the nesting check (done *before* the full run, as instructed)

Verified on two committed cells, both flip and tiered arm, using the
*unmodified* `survival.sample_flip_state` / `sample_tiered_state` (not the
`_p6` instrumented versions, to test the nesting property against the
literal originally-committed code path):

1. **Flip**, `c06_s01_cov050_seed00000000_bw000` (component_size=6,
   separation=1, coverage=0.5, base_wrongness=0.0, seed=0, `n_k=4`): drew
   1000 reps per depth via `sample_flip_state`, built the survival curve
   from only `rep < 200`, and compared to a fresh, independent `reps=200`
   draw and to `results/axis_robustness/survival_curves.csv`'s own row for
   this instance. All three agree **exactly** at every depth `d=1..4`
   (`n_eval`, `n_contradictory`, `S`, `S_contra_as_fail` bit-for-bit equal,
   including the `S=""` sentinel at `d=3,4` where every draw was
   contradictory).
2. **Tiered**, `c06_s01_cov100_seed00000002_tiered_nt2` (component_size=6,
   separation=1, `n_tiers=2`, seed=2): same check across the full
   `corruption_rate` grid (measured depth `d=0..5` realised for this
   instance). Exact match at every depth against
   `results/axis_robustness/survival_curves.csv`.

**Verdict: PASS.** This is expected, not just hoped for: `derived_seed`
hashes `(instance_id, "flip"|"tiered", d|corruption_rate, rep)` -- the seed
for rep `r` never depends on the *total* number of reps requested, only on
`r` itself, and each `(d, rep)` draw gets its own freshly constructed
`np.random.Generator`. There is no shared mutable state across reps in
`survival.py` that could make a longer run perturb an earlier rep's draw.
The empirical check confirms no accidental violation of that guarantee (e.g.
via hidden global RNG use), not just the guarantee itself. Consequence: the
`N=200` (committed) and `N=1000` (this rerun) survival numbers for the *same*
instance are a paired comparison on nested draws, not two independent
samples -- Track D can treat rep `r<200` of this rerun as literally the same
draw the committed run made.

Separately (Task 1 integrity check, also done before the full run): the
`_p6` instrumented samplers (`sample_flip_state_p6` / `sample_tiered_state_p6`)
were checked against the unmodified originals over 200 flip draws and 200
tiered draws each (same instances as above) -- `status`, `survived`,
`symdiff_proxy_not_distance` and `seed` match bit-for-bit in all 400 draws,
and the new denominator columns satisfy the expected identities
(`corruption_rejected + corruption_accepted == 1`;
`n_closure_orientations_changed == symdiff_proxy_not_distance` whenever
accepted; `closure_inert == 1` iff accepted and that quantity is 0) on every
row.

## Task 1 — corruption denominators (per draw)

Four denominators, reported separately as instructed, never merged into one
"survival rate":

- `corruption_rejected` / `corruption_accepted`: partitions every sampled
  draw (`rejected + accepted = 1`) by whether Meek closure refused the
  corrupted state (`apply_orientations` returned `None`).
- `n_claims_attempted` / `n_claims_reversed`: the claim-level target vs. the
  verified actual claim difference from the reference `K`.
- `n_closure_orientations_changed` / `closure_inert`: the directed-edge
  symmetric difference between the corrupted closure and `G0`, and whether
  an *accepted* corruption changed nothing (0 when rejected -- there is no
  closed graph to compare, not a hidden zero).

**Finding, reported per the "report it loudly" instruction, not buried:**
`n_claims_attempted == n_claims_reversed` on every single row checked so far
(400/400 in the pre-run integrity check above; **[FILL IN: full-run total
from the merged denominator totals]**). This is not a coding accident where
one column was aliased to the other -- both are computed independently (see
`survival_p6.py`'s `sample_flip_state_p6` / `sample_tiered_state_p6`) -- it
is a fact about the corruption operators themselves:
- `flip(k, rng, rate)` always reverses *exactly* `round(rate * n) = d`
  distinct, previously-unreversed claims (drawn without replacement from
  `k`'s own entries), so the flip arm's "targeted" `d` and the verified
  post-hoc claim difference coincide by construction; there is no partial-
  failure mode at the operator level.
- `tiered(...)`'s corruption knob is a *node*-relocation rate, not a
  claim-level target at all (a single misplaced node can flip several
  cross-tier claims at once, or none, depending on its edges) -- so the
  tiered arm has no independent "attempted claims" quantity to compare
  against; both columns are the same measured `|K_cor \ K_ref|` there by
  necessity, per the task brief's own parenthetical ("`d` in flip arm;
  measured in tiered arm").

The only place these four denominators diverge in an informative way is
therefore `corruption_rejected/accepted` and `closure_inert`, not
`attempted/reversed` -- i.e. the corruption *operators* never partially
fail; only the downstream Meek closure can refuse a state, and only among
accepted states can the closure end up unchanged (inert). **[FILL IN: final
per-stratum rejected/accepted/inert rates from the merged denominator
totals, flagging any stratum where `closure_inert` is a large share of
`corruption_accepted` -- an inert-but-accepted corruption trivially inflates
naive survival.]**

## Task 2 — new baseline predictors

- `separation`: already a column in the existing schema (it is the swept
  grid parameter); carried through here as its own predictor column
  together with `separation_status`. **`separation_status` is `"measured"`
  for every admitted row, by construction, not by omission**:
  `component_generator.generate_instance`'s own admission gate rejects
  (`reason = REASON_PARAMS_NOT_REALISED`) any draw whose
  `realised_separation` is not defined and exactly equal to the requested
  `spec.separation` -- so an admitted instance's separation is *never*
  undefined; this was verified by reading `component_generator.py`'s gate
  logic directly (lines ~665-671), not assumed.
- `k_g0`: number of closure orientations (directed edges in `G0` not
  directed in the CPDAG). **This is definitionally identical to the
  already-existing `shd_cpdag` column** (`survival.shd_cpdag_count(g0,
  cpdag)`), verified equal on every row of every smoke/pilot run performed
  so far. It is logged under its own name because that is the name
  `experiments/stage0_claim_radius.py`'s own output schema uses, not
  because it is an independent baseline -- **it carries zero information
  beyond `shd_cpdag`** and should not be treated as a fourth predictor
  alongside `phi_1`/`r_claim`/`separation` in any downstream analysis.
- `phi_1`: single-claim deletion fragility, matching
  `experiments/stage0_claim_radius.py` lines 140-151 exactly (drop one
  claim, re-close, test GAC-validity of `Z*`; `phi_1 = n_broken / n_k`,
  with claim-retractions that make the closure itself contradictory
  counted in `phi_1_n_inconsistent`, excluded from the numerator but not
  from `n_k`).
- `r_claim`: minimum number of claim *retractions* (deletions) that
  invalidates `Z*`, found by brute-force search capped at depth 3
  (`CLAIM_RADIUS_MAX_DEPTH` in `survival_p6.py`; `r_claim_censored=1` /
  `r_claim=""` if no hit within the cap). Depth 1 is exhaustive by
  construction (`phi_1_n_broken > 0` iff `r_claim == 1`), so the cap only
  bites at depths 2-3, where the search space (`C(n_k,2) + C(n_k,3)`) stays
  in the low hundreds for every `n_k` this grid produces (`n_k <=
  component_size - 1 <= 11`).

### The stated "known trap" — checked, not assumed, and it is narrower than stated

The task brief warned: *"at coverage 1.0 the admission gate is itself the
single-retraction test, so `phi_1 > 0` and `r_claim == 1` may hold by
construction on every admitted row."* A pilot run
(`component_size ∈ {6,8}`, `coverage=1.0`, `base_wrongness=0.0`, 8
instances/cell across all three achievable separations) shows this is
**real but strictly confined to the minimum-separation stratum**, not
universal across coverage=1.0:

| separation (of 3 achievable) | `phi_1` | `r_claim` |
|---|---|---|
| minimum (s=1) | `> 0` on 16/16 admitted instances | `= 1` on 16/16 |
| median (s=3 at c=6, s=4 at c=8) | `= 0.0` on 16/16 | `= 3` (found at the depth-3 cap) on all 16 |
| maximum (s=5 at c=6, s=7 at c=8) | `= 0.0` on 16/16 | **censored** (not found within depth 3) on all 16 |

This is coherent with the generator's own design intent
(`component_generator.py`'s docstring: separation is *designed* to be the
quantity "radius tracks"), not a bug: `r_claim` tracks `separation` almost
exactly in this pilot (`r_claim ≈ min(separation, 3)`), and the
"phi_1-degenerate" trap is a minimum-separation-stratum artifact, not a
coverage=1.0-wide one. **[FILL IN: full-grid confirmation/exceptions once
all coverage=1.0 cells across all four component sizes and three
base_wrongness values are in -- report any component_size/base_wrongness
combination where this pattern does NOT hold, exactly as it does not hold
uniformly for phi_1/r_claim degeneracy outside minimum separation.]**

`shd_truth` (an existing, not new, baseline) is separately known to be
constant (`≡ 0`) in the `coverage=1.0, base_wrongness=0.0` stratum
specifically (footnote `[^shdconst]` in `table_tau_comparisons.md`) -- that
degeneracy is unrelated to the phi_1/r_claim one above (different stratum
shape: it is bw=0.00-wide, not separation-narrow) and pre-dates this rerun;
noted here only so the two "predictor is constant somewhere" facts are not
conflated.

## Cells below quota (acceptance-rate, not budget, limited)

**[FILL IN AFTER MERGE]**: the committed `results/axis_robustness` run
(`seed_budget=150`, quota=20) had several structurally thin cells that did
not reach quota even at 150 seed attempts -- most severely
`tiered:c12_s1_nt2` (0/150) and `flip:c12_s11_cov1.0_b0.25` (1/150), both in
the high-base-wrongness / minimum-tiered-separation corner. This rerun uses
`seed_budget=600` (4x) to give the doubled quota (40) a fair chance, but a
structurally near-empty cell can still fall short -- report the final
accepted count for every cell here, flagging any still below 40/cell, with
its `exclusions` breakdown (`k_contradictory` vs `optimal_set_undefined` vs
`generate_instance:*` vs `fast_gate:*`), not just its raw count.

## s0 checks and monotonicity

**[FILL IN AFTER MERGE]**: `n_s0_failures` (must be 0 -- `Z*` is `G0`'s own
optimal set, so `S(0) == 1.0` is a checked invariant, not an assumption; a
nonzero count here means this run is broken and must be reported as such,
not smoothed over) and `fraction_non_monotonic` (measured, never enforced,
per `survival.fraction_non_monotonic`'s own docstring).

## Wall time

**[FILL IN AFTER MERGE]**: per-shard `elapsed_seconds` from each shard's
`_run_summary_<shard>.json`, their sum (total CPU-time), and the actual
wall-clock time from launch to the last shard's `SHARD_EXIT` in
`_logs/_launcher.log` (wall time < sum of shard times, since up to 5 ran
concurrently; the shared 10-core machine also had three other tracks
running, so wall time is not simply sum/5).

## What this rerun does *not* do

No τ_b, no paired bootstrap, no cross-predictor comparison against the
committed `N=200` table -- that is Track D's analysis, over these files.
This rerun's job ends at a clean, fully logged, nesting-verified set of
result files with the requested denominators and baselines attached.
