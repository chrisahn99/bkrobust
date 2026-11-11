# E7 — one hidden variable, near or far from the query (prereg, 2026-11-06)

Written 2026-11-06, before any run. Not edited after the first run; changes go in a dated addendum. Serves C2Bd-W2 (latent confounders) and h8pC-W2 (correct class assumed).

**Question.** If one variable of the true DAG is hidden before PC runs, how often does the committed adjustment set become invalid in the true DAG, and does that depend on how far the hidden variable sits from the query (x, y)?

**Populations.** The four Gaussian networks of the learned-CPDAG panel (ecoli70, arth150, magic-niab, magic-irri), the panel's 80 frame queries (20 per network, `frame.jsonl`), n = 5000, sample seeds 0-4: the panel's own samples (`sample_for`, seed `[20260912, stable_seed(net), 5000, s]`). PC as in the panel: Fisher-z, alpha 0.01, stable, then the Meek closure, or the raw graph where the closure fails.

**Hidden variables.** For query i of a network, with `rng = default_rng([20261106, stable_seed(net), i])`, draw one node among the true-skeleton neighbours of x or y (hop 1), and one among the nodes at skeleton distance 2 or more from both (hop 2+), each when the pool is not empty. The network's hidden set H is the union, the same for every seed. A unit is (network, seed, h) for h in H, plus the unit with nothing hidden. In a hidden unit, column h is dropped and PC runs on the rest. Every (unit, query) pair with h not in {x, y} is scored, and its stratum is set by that pair's own distance d(h, {x, y}) in the true skeleton: near if d = 1, far if d >= 2 (or unreachable).

**Arms.** A_ORACLE (registered): the truthful orientation of every undirected edge of the learned graph that exists in the truth; an edge absent from the truth gets no claim. D_DEGEN (reported): no knowledge. The analyst graph is `apply_orientations(learned, K)`; the committed set Z is `committed_adjustment_set` on it. Z is invalid at truth when Z exists and fails the generalised adjustment criterion (GAC) in the full true DAG, hidden node included.

**Registered statistic.** For each stratum s, Δ_s = mean over its pairs of [Z invalid at truth with h hidden] − [Z invalid at truth with nothing hidden, same network, seed and query], arm A_ORACLE. A pair with no committed set counts 0. 95% interval: percentile bootstrap, 2000 draws, seed 20261106, resampling hidden variables with replacement within each network (a drawn h carries all its seeds and queries).

**Prediction.** Δ_near > Δ_far, and Δ_far < 0.03 with the upper end of its interval below 0.10.

**Against us.** If Δ_far > 0.10, a hidden variable far from the query also breaks the committed set, and no local check covers it. We then restate that the method assumes causal sufficiency (Section 3) and report the result as a limitation, without claiming robustness to latents anywhere.

**Also reported, not registered.** Δ per network, per seed, for D_DEGEN, and the transition share P(invalid with h | valid without). The commit rate per stratum. The GAC radius (A_ORACLE) on hidden and full graphs, and among pairs that became invalid, the share whose radius was 2 or more or unreached (the radius gave no warning). Pairs where h is a common ancestor of x and y.

**Radius.** `criterion="gac"`, search budget 3. On graphs above 60 nodes (arth150, magic-irri) the SAT ladder is skipped, as in the panel sweep; a search that runs out of budget gives the lower bound 4. Ladder rungs 20 s; 120 s wall per row (`wall_timeout`).

**Consistency checks (exit 3).** (1) Every unit regenerates the panel sample hash. (2) The nothing-hidden graph equals the stored panel graph (`graphs.jsonl`). (3) GAC at truth by the submitted code and by the networkx audit agree. (4) Z is GAC-valid on its own analyst graph and its radius is not 0. (5) On about one radius row in ten (fixed by a hash of unit and query) of a graph of at most 60 nodes, the submitted `breakdown_radius` agrees with the dispatch when both are exact, and neither is below the other's lower bound. (6) Nothing hidden: truth extends the analyst graph implies Z valid at truth.

**Exclusions and stopping rule.** Pairs with h in {x, y} are not formed. A PC run past 900 s is a `timeout` unit; its pairs are dropped and counted, and a network with more than 5% timed-out units is flagged. One launch of the full grid, resumed only after a crash; no reruns with changed settings. A failed consistency check is fixed and rerun, and both runs are reported.

## Deviations from the plan
- The plan says "hide one variable at hop 1, or at hop 2 or more". Hidden variables are drawn per query but each PC run is scored on every query of its network, so far pairs outnumber near pairs; the strata are defined per pair.
- The plan cites `ledger_sweep.py`; its `instrument_row` passes options the submitted `breakdown_radius` lacks, so the radius uses a small dispatch in `run.py` over the same submitted functions.
