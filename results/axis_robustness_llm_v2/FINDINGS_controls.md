# Controls and statuses on the LLM-elicited panel (v2)

Script: `experiments/llm_controls_v2.py`; output: `controls_k_accuracy.json`.
Nothing under `results/matched_control/` or `results/axis_robustness_llm/` was modified.

## Matched random-knowledge control (review point 4.7)

`results/matched_control/arms/D_RANDOM_MATCHED` holds, per network, |K|-matched
random orientations of the primary model's own CPDAG pairs (Meek-checked, 20
seeds). Claim accuracy k_acc against the data-generating graph
(network-cluster bootstrap, 10,000 resamples):

| knowledge | k_acc [95% CI] |
|---|---|
| primary model (Qwen2.5 7B) | 63.8% [49.0, 78.1] |
| eight real-name conditions pooled | 68.4% [57.6, 78.5] |
| matched random | 56.1% [51.5, 60.5] |

Paired per network: primary − random = +0.037 [−0.103, +0.172];
pooled real-name − random = +0.085 [−0.020, +0.183]. Neither excludes zero.

The survival panel was not run for the random arm: its radius computation alone
took about 8.75 h (`results/matched_control/_logs/random_rand{A,B,Huge}.log`).

## Status taxonomy (review point 6)

On `analysis_units.csv` (5,400 units):

| status | units |
|---|---|
| finite r_val (scored) | 2,905 |
| `optimal_set_undefined` (no valid set, OR empty valid set; see below) | 1,542 |
| `o_g0_extensions_intractable` | 743 |
| `n_k_zero` (no claims asserted) | 210 |
| Z invalid at G0 (r = 0) | 0 (pre-filtered) |
| r = inf | 0 observed |
| timeout | 0 observed |

Knowledge repaired for consistency: 74 of 320 (condition, network) states, 474
claims dropped (upstream, `results/elicit/knowledge.json`).

`real_survival.commit_z_star` tests `if not z:`, which is true both for an
unidentified effect (None) and for an identified, valid, empty optimal set, so
the two are merged in `optimal_set_undefined`; `status_split.json` separates them.

The 33 historical `r_val = 0` rows are in `results/final_table/instances.jsonl`
(the older `experiments/final_table.py` pipeline): there `Z` is the empty set
and already invalid at G0, because that pipeline does not check validity at G0
before calling `breakdown_radius`; `final_table.py` reports them as
`n_degenerate` and excludes them from averages.
