# E23 preregistration: PC and GES agreement as a check for discovery error

Dated 2026-11-06. Written before any run. Not edited after the first full run; later changes go in DEVIATIONS.md with a date.

**Question.** When PC and GES are run on the same sample, does disagreement between the two learned states about the committed adjustment set flag the sets that are invalid in the generating DAG? Bootstrap stability did not (E6: AUC 0.245). Serves reviewer points C2Bd-W2, C2Bd-Q3 and h8pC-W2.

**Data.** Stored panel only, no discovery is rerun. `experiments/rebuttal/_data/pest/graphs.jsonl` and `graphs_retry.jsonl` hold the PC and GES graphs learned from the same sample (`sample_for(network, n, seed)`) for each (network, n, seed). A pair exists when both methods have a graph for that triple. arth150 has GES graphs only for seeds at n = 1000 and 5000 (2 each), so it contributes few pairs; the per-network row counts are reported.

**Population.** Every certifiable row (status `ok`) of arm A_ORACLE_R on a learned graph that has a partner graph from the other method, over the 80 frame queries, n in {200, 1000, 5000}, both directions (PC row checked against the GES graph, GES row checked against the PC graph). The source state is rebuilt with E5's code (`plans_transfer`, `apply_orientations`, `committed_adjustment_set`). If E5's rows are present, the rebuilt status, set size, knowledge hash and validity at truth must equal E5's row for the same unit, and the set of certifiable units must be the same.

**Procedure.** For each row (graph G, query x, y, committed set Z):
1. Build the partner state: the partner graph G' with its own A_ORACLE_R knowledge (`plans_transfer` on G'), closed with `apply_orientations`. This is the same knowledge rule applied to the other method's graph.
2. `z_valid_other` = Z passes the generalised adjustment criterion (GAC) in that state (`is_gac_valid_mpdag`). A state that cannot be closed counts as not valid.
3. `z_same_other` = the committed set read off the partner state equals Z.
4. Registered flag: `flag = 1 - z_valid_other`. Union flag (secondary): `flag_or_diff = flag or not z_same_other`.
5. Also recorded, secondary: the same check when the source knowledge itself is moved to G' (claims on edges undirected in G', dropped greedily if they do not close).

**Registered statistic.** AUC (Mann-Whitney, ties one half) of `flag` for "Z invalid at truth" (GAC on the true DAG), pooled over the three sample sizes and both directions. 95% interval: two-stage cluster bootstrap as in E6, 2000 draws, seed 20261106: resample networks with replacement, then query pairs (x, y) within each drawn network, keeping every row of a drawn pair (all n, seeds and directions). Draws without both classes are dropped and counted.

**Secondary, descriptive.** The same AUC by n and by direction; the AUC of the union flag and of the moved-knowledge flag; per-network row counts, invalid counts and AUC; the shared-error share, the share of invalid-at-truth rows on which the partner state also finds Z valid (flag 0); the false-alarm share, valid-at-truth rows with flag 1.

**Prediction.** PC and GES err differently (E5: invalidity rises with n for PC and falls for GES), so disagreement carries some signal. We predict AUC between 0.6 and 0.75 and a shared-error share of at least 0.3.

**Bar.** 0.65.

**Against the authors.** An AUC of 0.65 or below, or a lower bound under 0.5. The two algorithms then make the same mistakes near the query, and agreement between them does not detect discovery error, just as resampling did not. The response then says that neither resampling nor a second algorithm certifies the equivalence class, and that the radius is conditional on the class. If the AUC is above 0.65 with a lower bound at or above 0.5, the response offers cross-algorithm agreement as a companion check and prints the shared-error share beside it.

**Consistency checks (exit 2).**
1. A rebuilt source row that differs from E5's row (status, set size, knowledge hash, validity at truth), or a different set of certifiable units, when E5's rows are present.
2. Z not GAC-valid in its own source state (the committed set must pass on the graph it was read from).
3. A pair whose two graphs are identical but whose flag is 1.

**Seeds, stopping rule, exclusions.** Deterministic: no random draws except the interval seed. One full launch, resumed only after a crash. No rerun with changed settings. Rows without a partner graph are counted and excluded. No other row is dropped.

**Criterion.** GAC throughout (state: `bkrobust.gac.mpdag_level`; truth: the networkx GAC audit copied from E5).
