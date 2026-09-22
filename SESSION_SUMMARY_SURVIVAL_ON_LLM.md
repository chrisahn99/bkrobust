# Session summary: survival and ranking on elicited LLM knowledge

Branch `revision/survival_on_llm`. Night session, 2026-09-22.

## What the session set out to do

Two things, both about where the analyst's background knowledge `K` comes from.

**Point 1.** The survival track built `K` with
`benchmarks.measure.select_knowledge(dag, cpdag, coverage)`, which reads the
**ground-truth DAG** and returns claims that are correct by construction. The
only wrongness in that design was the synthetic `base_wrongness` dial applied
on top. The claims an LLM actually asserted, sitting in
`results/elicit/knowledge.json` since the elicitation sessions, were never on
that path. Put them on it.

**Point 2.** Use the model panel as the variation axis. Within one
`(network, coverage)` cell of the committed corpus, `|K|` is fixed and
`shd_truth` is identically zero at `base_wrongness = 0`, so there is nothing
to rank. Ten elicitation conditions give ten different knowledge states per
network instead.

Tasks 3 (an empirically-grounded error generator to replace uniform flip) and
4 (upstream query-timeout and denominator fixes) were explicitly out of scope
and were not attempted.

## What was written

Three new modules and one runner. All lint-clean under the repo's ruff config.

| file | what it does |
| --- | --- |
| `src/bkrobust/robustness/llm_survival.py` | Reads `K` from `knowledge.json` for a `(condition, network)` pair. Measures `k_accuracy`, `assert_rate` and the off-skeleton audit. Builds `G0` by Meek closure of the elicited claims. |
| `src/bkrobust/robustness/run_llm_survival.py` | Sharded, resumable driver. One shard per `(condition, network)`; 240 shards. `list` / `run` / `pool` / `status` / `manifest`. |
| `src/bkrobust/robustness/llm_analyse.py` | Ranking tables. Four stratifications; two endpoints; Kendall tau-b with a cluster bootstrap over networks. |
| `scripts/run_llm_survival_panel.sh` | The four-stage runner: manifest, status, pool, analysis. Verbose, timestamped, append-only logging to `results/axis_robustness_llm/_logs/`. |

### `real_survival.py` was not modified, and that is deliberate

The brief said to modify `real_survival.py` "or the equivalent survival/ranking
script". Its own module docstring records the repository's standing rule: the
synthetic counterpart is left unmodified so every number in
`report_fragility_and_pareto.md` stays reproducible from the code that produced
it, and the committed 831-row corpus in `results/axis_robustness_real/` depends
on `real_survival.py` the same way. Editing it in place would have put those
frozen results at risk for no scientific gain.

Instead `llm_survival.py` is its elicited-knowledge sibling and **imports**
`real_survival`'s grid, seeding, closure cache, radius machinery and cell
statistics wholesale — `depth_grid`, `derived_seed`, `flip_draw`,
`ClosureCache`, `radius_columns`, `cell_statistics`, `auc_over_grid`,
`auc_usable`. The two arms therefore cannot drift apart on what a grid point,
a seed or an endpoint means, which is the property that mattered. Likewise
`llm_analyse.py` imports `_build_unit` and `tau_for_stratum` from
`real_analyse.py` rather than reimplementing them.

`knowledge.json` was **read only**, never modified.

## What was run

```
scripts/run_llm_survival_panel.sh --workers 8
```

- 240 shards = 10 conditions x 24 networks. **0 failures.**
- Sweep wall clock **2284 s** (38 min), 8 workers. Analysis a further ~200 s.
- 1000 draws per grid point, the pre-registered count. No wall-cap censoring
  anywhere (`n_cells_censored_wall_cap == 0` on all 240 markers).
- Long tail: the three `diabetes` shards (413 nodes, 89 pairs) took ~2200 s
  each, run in parallel; everything else finished inside 212 s.
- 5400 analysis units, 5400 instance rows, 19 870 cell rows.

Outputs under `results/axis_robustness_llm/`: `manifest.json`, `_done/`,
`shards/`, `_logs/`, `analysis_units.csv`, `analysis_tau.csv`,
`analysis_panel.csv`, `analysis_per_condition.csv`, `analysis_summary.json`,
and four verification artefacts described below.

## Integrity checks

`ORCHESTRATOR_VERIFY.txt`, produced by `orchestrator_verify.py`, which imports
nothing from the analysis harness and recomputes everything from the shard
JSONL and the bundle:

| check | result |
| --- | --- |
| Each shard's `K` equals `knowledge.json`'s `k` for its cell | 240/240, 0 mismatches |
| `k_b_sha256` recomputed from the bundle, per instance row | 5400/5400 exact |
| `k_accuracy` recomputed against the true DAG | 0 mismatches |
| Claims outside the CPDAG skeleton | 0 |
| `k_source` on every instance row | `knowledge.json` x 5400 |
| `S` defined iff `n_eval > 0` | 0 violations |
| Forbidden `seconds` key on any row | 0 |
| Wall-cap censored cells | 0 |
| Networks with >1 distinct `\|K\|` across conditions | 24/24 |

**`select_knowledge` is not on the execution path.** This was traced rather
than asserted. The name appears in `llm_survival`'s dependency graph because
`real_survival.py` binds it at module scope, but the only call site in the
whole graph is inside `run_real_survival.run_flip_shard`, which this arm never
enters. `run_llm_survival` imports exactly four names from that module —
`JsonlWriter`, `_censored_cell`, `load_frame`, `sha_file` — all pure I/O or
hashing, none of which reference it.

**`pathfinder`** is in the frozen frame but absent from every condition of the
bundle: it was never elicited. Its 3 frame rows produce no shard and the
exclusion is stamped in `manifest.json` under `networks_absent_from_bundle`,
so the network count cannot silently drop from 25 to 24 unremarked.

## The panel

Per-network `|K|` varies across conditions in **all 24** networks — `ecoli70`
spans 3 to 23 claims over ten distinct values, `sachs` 1 to 14 over eight. The
two-dimensional variation the design needed is present.

Accuracy is reported two ways because they rank the conditions differently and
neither is privileged. Pooled weights a network by how many claims the model
made there; the network-mean gives each network one vote.

| condition | model | naming | Σ\|K\| | med \|K\| | acc (pooled) | acc (net-mean) |
| --- | --- | --- | --- | --- | --- | --- |
| `D_LLM` | qwen2.5:7b-instruct | real | 130 | 5.0 | 0.677 | 0.638 |
| `D_LLM_INSTR` | qwen2.5:7b-instruct | real | 109 | 3.0 | 0.606 | 0.641 |
| `D_LLM_GEMMA_27B` | gemma-3-27b-it | real | 146 | 6.0 | 0.651 | 0.669 |
| `D_LLM_32B` | Qwen2.5-32B-Instruct | real | 158 | 6.0 | 0.715 | 0.692 |
| `D_LLM_72B` | Qwen2.5-72B-Instruct | real | 191 | 7.0 | 0.660 | 0.716 |
| `D_LLM_72B_INSTR` | Qwen2.5-72B-Instruct | real | 174 | 6.0 | 0.655 | 0.735 |
| `D_LLM_SRC` | llama3.1:8b | real | 160 | 4.5 | 0.675 | 0.684 |
| `D_LLM_SRC_70B` | Llama-3.3-70B-Instruct | real | 185 | 7.5 | 0.697 | 0.699 |
| `D_SCRAMBLED` | qwen2.5:7b-instruct | scrambled | 108 | 4.0 | 0.630 | 0.694 |
| `D_SCRAMBLED_72B` | Qwen2.5-72B-Instruct | scrambled | 140 | 4.5 | **0.743** | 0.721 |

### Scale buys claims, not accuracy

Paired per-network sign tests along each family ladder:

- **`|K|` rises with scale, significantly.** 7B → 72B: 18/20 and 21/23
  networks up, sign-test p ≈ 4e-4 and 7e-5.
- **Accuracy does not.** Every ladder step is non-significant
  (p from 0.21 to 0.97), and 32B is not monotonic between 7B and 72B.

So the panel's `|K|` axis is a genuine model-scale effect; its accuracy axis is
not. The brief's framing of the panel as supplying "varying baseline
accuracies" from model scale is not supported — accuracy varies, but not with
size.

### The scrambled control fails again

`D_SCRAMBLED*` answer about relabelled variables and should be a semantic
floor. Paired per-network, they are not:

| comparison | real > scrambled | scrambled > real | ties | sign-test p |
| --- | --- | --- | --- | --- |
| `D_LLM` vs `D_SCRAMBLED` | 7 | 12 | 4 | 0.359 |
| `D_LLM_INSTR` vs `D_SCRAMBLED` | 6 | 11 | 6 | 0.332 |
| `D_LLM_72B` vs `D_SCRAMBLED_72B` | 11 | 10 | 1 | 1.000 |
| `D_LLM_72B_INSTR` vs `D_SCRAMBLED_72B` | 11 | 8 | 3 | 0.648 |

At 7B the scrambled arm is *better* more often than worse. This reproduces the
earlier finding that the LLM attribution fails its matched control: the claims
are not demonstrably driven by real domain semantics rather than skeleton
structure and base rates. It is a property of the elicitation, not of anything
introduced this session, but it bounds what the panel can be said to measure.

## The ranking result

Pooled over the 8 real-naming conditions. n = 2422 units, 24 network clusters,
Kendall tau-b, 10 000-resample cluster bootstrap over **networks** (a network's
rows are dependent twice over here: one corrupted state is shared by all its
pairs within a condition, and the same structure recurs across conditions).

| predictor | tau_b | 95% CI | LOO flips |
| --- | --- | --- | --- |
| **`radius` (r_val)** | **+0.319** | **[0.122, 0.485]** | no |
| `k_g0` | +0.214 | [0.079, 0.339] | no |
| `n_k` | +0.185 | [0.071, 0.327] | no |
| `shd_truth` | +0.054 | [−0.077, 0.200] | no |
| `k_accuracy` | +0.001 | [−0.132, 0.135] | **yes** |
| `assert_rate` | −0.061 | [−0.173, 0.063] | no |

**`r_val` has the largest point estimate but is not shown to beat the
strongest baselines.** Paired network-cluster bootstrap on the *difference* in
tau, computed twice independently:

| difference | mean | 95% CI | verdict |
| --- | --- | --- | --- |
| `radius − n_k` | +0.132 | [−0.082, +0.321] | includes zero — **not shown to beat** |
| `radius − k_g0` | +0.112 | [−0.087, +0.300] | includes zero — **not shown to beat** |
| `radius − shd_truth` | +0.261 | [+0.038, +0.463] | excludes zero — **beats** |

The defensible claim is: *on knowledge an analyst actually supplied, `r_val`
ranks queries by survival, and it beats the SHD-to-truth baseline. It is not
shown to beat `k_g0` or `|K|`.* This is a weaker claim than the committed
corpus supports, and it is the one the data supports here.

### The radius concentration is gone

The motivation for this work was that the committed evaluations piled up at
radius 1 with no variance to rank. On elicited `G0`, `r_val` spreads from 1 to
13 with roughly half the mass at 1:

| condition | n | median | max | share(r=1) |
| --- | --- | --- | --- | --- |
| `D_LLM_SRC_70B` | 305 | 1 | 13 | 0.505 |
| `D_LLM_SRC` | 283 | 1 | 13 | 0.534 |
| `D_LLM_72B_INSTR` | 315 | 1 | 13 | 0.556 |
| `D_SCRAMBLED` | 223 | 2 | 11 | 0.453 |

(`UNREACHED` excluded throughout, never averaged; 0 rows carried it.)

## Three caveats that limit the result

### 1. The endpoint does not survive rescoring contradictory draws

Elicited `K` is dense on an almost-determined CPDAG and Meek-consistent by
construction, so reversing claims very often yields a claim set admitting no
MPDAG. Pooled contradiction rate by depth fraction:

```
d/|K|  0.1    0.2    0.3    0.4    0.5    0.6    0.7    0.8    0.9    1.0
       0.478  0.555  0.606  0.759  0.699  0.819  0.714  0.840  0.846  0.593
```

`S` is correctly undefined at saturated grid points, and `auc_over_grid` then
falls back to the nearest defined point — so `AUC_frac` on this arm is carried
by the low-depth end. The second endpoint `AUC_frac_contra_as_fail` scores a
contradictory draw as a failure instead of dropping it, computed through the
*same* `auc_over_grid` from the `S_contra_as_fail` column every cell already
carries. Under it **every predictor's CI includes zero and `radius` flips
sign** (−0.120, [−0.285, 0.129]).

So the ranking result is conditional on how contradictory corruptions are
scored, and at these rates that is not a minor modelling choice. Note also the
**non-monotonicity**: reversing *every* claim (0.593) is less often
contradictory than reversing 90% of them (0.846). The uniform flip operator is
a poor corruption model for elicited knowledge. This is the strongest argument
the session produced for Task 3, which was out of scope.

### 2. Only 54% of units are scorable, and the sample is selected on `|K|`

| status | n | share |
| --- | --- | --- |
| `ok` | 2905 | 0.538 |
| `optimal_set_undefined` | 1542 | 0.286 |
| `o_g0_extensions_intractable` | 743 | 0.138 |
| `n_k_zero` | 210 | 0.039 |

Partial elicited `K` frequently fails to determine an optimal adjustment set at
all — a real finding about what LLM-supplied knowledge delivers. But it means
entry into the ranking sample depends on `|K|`, the predictor being ranked.
The selection is real and runs in the inflating direction: ok-share against
mean `|K|` across conditions gives Spearman rho = 0.81 (p = 0.015), and
**24/24 networks individually** show a positive tau between `n_k` and being
scorable (median 0.47). `optimal_set_undefined` concentrates at small `|K|`
(0.507 of units at `|K|` 1–2, falling to 0.213 at `|K|` ≥ 11).

The `panel_balanced` stratification restricts to triples scorable in *every*
real-naming condition, which removes the channel entirely — but there are only
31 such triples across 7 networks:

| predictor | tau_b (balanced, 7 clusters) | 95% CI |
| --- | --- | --- |
| `radius` | +0.499 | [0.240, 0.624] |
| `k_g0` | +0.380 | [−0.074, 0.562] |
| `n_k` | +0.333 | [−0.125, 0.484] |
| `shd_truth` | +0.285 | [−0.167, 0.456] |

**This does not demonstrate that the pooled `n_k` effect is a selection
artefact**, and it would be wrong to read it that way. An intermediate panel
(scorable in ≥6 of 8 conditions; 205 triples, 20 networks) gives `n_k` tau =
+0.255 with a CI that *excludes* zero — a **higher** point estimate than
pooled. Across balancing levels 8/8 → 7/8 → 6/8 → pooled, `n_k` runs
0.333 → 0.262 → 0.255 → 0.185, with the CI crossing zero only in the two
smallest panels. That pattern tracks cluster count (7 → 16 → 20 → 24) at least
as well as it tracks selection strength. The honest statement is that the
selection exists and points the right way to matter, and this data cannot
separate it from a power effect.

`radius` survives every balancing level with a CI excluding zero, which is the
one thing the control does establish.

### 3. `k_accuracy` ranks nothing

tau = +0.001 with an LOO verdict flip. How *correct* an elicited knowledge
state is carries no information about how well its queries survive corruption,
on this corpus. Given caveat 1 — the corruption operator does not resemble the
errors models actually make — this is not yet evidence that accuracy is
irrelevant, only that it is irrelevant to *uniform reversal* of the claims.

## Verification

Four artefacts, three of them produced by independent workers that were given
the definitions and told to reimplement from scratch rather than import the
harness:

| file | what it checks |
| --- | --- |
| `ORCHESTRATOR_VERIFY.txt` | K provenance, accuracy recomputation, sentinel discipline, radius spread, censoring — nothing imported from the analysis harness |
| `VERIFY_PANEL.txt` | `select_knowledge` call-graph trace, `k_b_sha256` fidelity per row, panel description, scale ladders, scrambled control, degenerate census |
| `VERIFY_SELECTION.txt` | the selection audit, balanced-panel re-derivation, intermediate balancing, paired radius-vs-baseline bootstrap |
| `VERIFY_ENDPOINT.txt` | independent recomputation of `AUC_frac`, `AUC_frac_contra_as_fail` and every pooled tau, plus file-digest verification |

Results of the independent recomputations:

- **Per-unit endpoints.** All 5400 units recomputed from the shard JSONL with
  the grid, nearest-point and tie rules reimplemented from the spec:
  `AUC_frac` and `AUC_frac_contra_as_fail` both **0 mismatches** against
  `analysis_units.csv` at 1e-9.
- **File digests.** 240/240 shards, both JSONL files each, recomputed sha256
  matches the completion marker. 0 mismatches.
- **Cell arithmetic.** 19 870 cells: 0 with `S` defined at `n_eval == 0`, 0
  with `n_eval + n_contradictory != n_draws`, 0 with `n_survived > n_eval`, 0
  carrying a `seconds` key, and `n_draws` uniformly 1000.
- **Pooled taus.** The five `AUC_frac` taus reproduce to <=1e-6. The five
  `AUC_frac_contra_as_fail` taus were flagged as disagreeing, but every diff
  is between 2.5e-6 and 3.1e-5 on identical row counts. The harness
  reproduces its own `analysis_tau.csv` from `analysis_units.csv` to nine
  decimal places, and the worker's per-unit values for that endpoint matched
  the same CSV to 1e-9, so the residual sits in the worker's tau step rather
  than in the pipeline. It is four orders of magnitude below the CI widths and
  changes no verdict; it was not chased further and is recorded here as
  unresolved on the verification side, not as a pipeline defect.

The paired bootstrap on `radius − baseline` and the scrambled-control sign
tests were each recomputed a second time by the orchestrator and agree to
three decimal places.

## What this changes, and what is still open

The survival/ranking track now runs on knowledge an analyst actually supplied.
The radius-1 concentration that made ranking vacuous is gone, and `r_val`
ranks. But the claim it supports is narrower than the committed corpus's: it
beats SHD-to-truth, not `k_g0` or `|K|`, and the whole result is conditional
on a corruption operator that visibly does not fit elicited knowledge.

Open, in priority order:

1. **Task 3 — the empirically-grounded error generator.** Contradiction rates
   of 0.48–0.85 and their non-monotonicity in depth say the uniform flip is
   wrong here. The elicitation data needed to fit a real error model (which
   pairs each model gets wrong, and how) is already in `knowledge.json`.
2. **The `optimal_set_undefined` population.** 28.6% of units. Whether an
   analyst with partial LLM knowledge can form an adjustment set at all is
   arguably a more interesting question than how that set survives corruption.
3. **`o_g0_extensions_intractable`** at 13.8%, concentrated in `arth150` and
   `diabetes` — a measurement limit of this machine, not a structural fact,
   and it removes two of the largest networks from most conditions.
