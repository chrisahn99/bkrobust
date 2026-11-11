# E24 — real far-away discovery errors transplanted into the true graph (prereg, 2026-11-06)

Written 2026-11-06, before any run of this code. Not edited after the first run; changes go in DEVIATIONS.md.

**Question.** Reviewers C2Bd (W2) and h8pC (W2, W4a) ask whether the radius, which assumes a correct class, can be defended by checking the class only near the query. If only the discovery errors that PC and GES actually made away from the query are kept, and the truth is restored near the query, do committed sets still become invalid in the true DAG?

**Population.** Every learned graph stored for the P-EST panel (PC and GES, n in {200, 1000, 5000}, all seeds in `graphs.jsonl` and `graphs_retry.jsonl`; networks ecoli70, arth150, magic-niab, magic-irri; 214 graphs) and the 80 stored frame queries (x, y). A unit is a (learned graph, query) pair that is certifiable under the repaired oracle arm A_ORACLE_R of E5: G0 = the learned CPDAG with the truthful orientation of every undirected edge that exists in the truth, added greedily in sorted order, skipping any that do not close; the committed set of G0 exists. E5 found 987 such units (PC 408, GES 579). This run recomputes them and must agree with E5's rows on status and validity at truth.

**Hybrid DAG.** N_d is the set of nodes within d hops of {x, y} in the true skeleton. One consistent DAG extension E of G0 is taken, deterministically (Dor and Tarsi's algorithm, sinks taken in sorted order). The hybrid H_d has, on every pair with both ends in N_d, the true DAG's edge or non-edge, and on every other pair the edge or non-edge of E. Levels d = 1, 2, 3 are the experiment. Two control levels are also run: d = none (H = E, no transplant) and d = all (N covers every node, so H must equal the true DAG). A cyclic H_d is discarded and counted; a unit whose G0 has no consistent extension is discarded and counted at every level.

**Arm on the hybrid.** C_d = CPDAG(H_d). Knowledge = truthful orientations of the undirected edges of C_d that exist in the truth, repaired as in E5 (added in sorted order, skipping any that do not close). The committed set is computed on the resulting graph by the same closed-form rule as E5 and tested in the true DAG with the generalised adjustment criterion (criterion "gac", networkx implementation).

**Registered statistic.** At d = 2, pooled over PC and GES: the share of hybrids with a committed set whose set is invalid in the true DAG. 95% Wilson interval with the sample size divided by the design effect of clustering by learned graph (ratio-estimator variance over graphs; design effect floored at 1). Set against the true-class anchor (E5: 0 of 25 invalid on the true CPDAGs, recomputed here) and the full learned graphs (E5, PC: 0.54 at n = 200 to 0.66 at n = 5000). Reported with it: per d, the number of cyclic hybrids, of units without a consistent extension, and of hybrids with no committed set.

**Prediction.** Restoring the truth within two hops removes most invalidity: the d = 2 share is below 0.05.

**Against the authors.** A d = 2 share of 0.05 or more (point estimate). Errors far from the query then break committed sets through orientation propagation, no local check of the class suffices, and the Limitations paragraph says so in those words with the measured share.

**Also reported, not registered.** The share at d = none, 1, 3, per method and per network; the share restricted to hybrids whose remaining errors all lie at two or more hops (min endpoint distance of every wrong pair at least 2), since a learned extra edge from x to a far node is kept by the construction above; counts of remaining skeleton and orientation errors.

**Consistency checks (exit non-zero).**
1. A unit's recomputed status or validity on G0 differs from E5's A_ORACLE_R row (when E5's rows are present).
2. E is not a consistent extension of G0 by `is_consistent_extension`.
3. H_d differs from the truth on a pair inside N_d.
4. At d = all, H differs from the true DAG, or a committed set is invalid in the truth.
5. The truth is a consistent extension of the hybrid's knowledge graph, yet the committed set is invalid in the truth.
6. The incremental knowledge repair on a hybrid differs from the E5 repair (checked on the first query of every graph).
7. The true-class anchor is not 25 certifiable rows with 0 invalid.

**Seeds.** None needed: every step is deterministic (sorted orders throughout).

**Stopping rule.** One launch of the full grid, resumed only after a crash. No reruns with changed settings. If a consistency check fails, we fix the code and rerun from scratch, and report both runs.
