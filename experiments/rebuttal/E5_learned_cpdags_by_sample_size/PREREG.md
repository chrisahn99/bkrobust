# E5 — validity at truth and radius on learned CPDAGs, by sample size (prereg, 2026-11-06)

Written 2026-11-06, before any run. Not edited after the first run; changes go in a dated addendum below.

**Question.** When the equivalence class is learned from data, how often is the set the analyst commits to invalid in the true DAG, and does that share fall as n grows from 200 to 5000? Serves reviewer points h8pC-W2 and C2Bd-W2.

**Phase A (reproduction gate).** Recompute all 560 rows of `results/pest/rows.csv` (244 certifiable) on the reference learned CPDAG (PC, n = 5000, seed 0) with the submitted code. Every deterministic field (status, committed set size, claim radius, validity at truth, buckets B2/B4/B5/B7, knowledge hash) must match exactly. The stored hop radius was computed with the back-door predicate, so it is compared with the submitted `breakdown_radius(criterion="backdoor")`: equal where both are exact, and never below the other's anytime lower bound. The reported radius is the GAC one (`criterion="gac"`). One disagreement makes the run exit non-zero and stops interpretation.

**Phase B.** Every stored learned graph (PC and GES, n in {200, 1000, 5000}, all seeds in `graphs.jsonl` and `graphs_retry.jsonl`, 214 graphs) plus the true CPDAG of each network as the no-discovery-error anchor. The 80 frame queries (x, y) are held fixed. Arms are as in the P-EST panel, restricted to each graph: D_LLM, D_RAND (3 reps), D_DEGEN, A_TRUE, A_ORACLE. A sixth arm, A_ORACLE_R, takes the truthful orientation of every undirected edge of the graph that exists in the truth, adding them one by one in sorted order and skipping any that the learned arrows make impossible to close. The stored A_ORACLE, by contrast, becomes `knowledge_inconsistent` on the whole graph as soon as one such edge conflicts (all 20 arth150 rows on the reference graph). For every row we record status, the GAC hop radius, the claim radius, validity at truth (GAC on the true DAG, networkx code), and whether the truth extends the analyst's graph.

**Registered statistic.** Arm A_ORACLE_R, method PC. At each n, the pooled share of certifiable rows (status ok) whose committed set is invalid at the truth. The statistic is Δ = share(n = 5000) − share(n = 200). Its 95% interval comes from a percentile bootstrap (2000 draws, seed 20261106) that resamples graphs (seeds) with replacement within each (network, n) cell.

**Prediction.** Invalidity falls with n: Δ < 0 with the upper end of the interval below 0. The true-CPDAG anchor gives 0 invalid on A_ORACLE_R and D_DEGEN, which holds by construction.

**Against us.** If the interval of Δ reaches 0 or above, invalidity does not fall with n, and discovery error is systematic here. We will then say that the radius is conditional on the class and that discovery error needs its own control (a skeleton-stability check beside the radius). We will not claim that the error shrinks with data.

**Also reported, not registered.** The same Δ for every arm and for GES. The share computed with the denominator held fixed (invalid and certifiable over all frame rows). The certifiable share. The distribution of GAC radii (finite, unbounded, lower bound only). Counts of B6/B7 dangerous rows. Each per network.

**Consistency checks (exit non-zero).**
1. Any Phase A disagreement.
2. A row where the truth is a consistent extension of the analyst's graph, yet the committed set is invalid at the truth.
3. A disagreement between the GAC radius by criterion and by enumeration, both exact, on graphs with at most 8 undirected edges.
4. Phase B rows on the reference graph that differ from the Phase A rows.

**Populations, seeds, exclusions.** Networks: ecoli70, arth150, magic-niab, magic-irri. D_RAND seeds: [20260913, stable_seed(network), rep], as in the stored sweep. No row is dropped. Knowledge that cannot be closed is a counted `knowledge_inconsistent` row. A radius that hits the per-row wall (120 s) is recorded as `wall_timeout`. The committed set and the validity at truth are still computed for that row. A cell with more than 5% timeouts is flagged in the table.

**Stopping rule.** One launch of the full grid, resumed only after a crash. No reruns with changed settings. If a consistency check fails, we fix and rerun, and report both runs.

## Deviations from the plan
- D_LLM knowledge on each learned graph is the stored supplier's claims restricted to edges left undirected in that graph. Claims that cannot be closed are dropped greedily, in stored order. There is no new elicitation, so edges the supplier was never asked about stay open.
- On graphs above 60 nodes (arth150, magic-irri) the SAT ladder is skipped, as in the stored sweep. The search answers alone, and a search that runs out of budget gives the lower bound depth + 1. The submitted `breakdown_radius` has no such switch, so this dispatch is a small wrapper in `run.py` that calls the same submitted functions.
- The radius uses the GAC (reported). Phase A also computes the back-door radius, only to compare with the stored numbers.
