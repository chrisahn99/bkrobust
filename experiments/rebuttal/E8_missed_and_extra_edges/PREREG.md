# E8 — Missed and extra edges: how local is the damage (preregistration, 2026-11-06)

Written before any run of this code. Not to be edited after the first run; later notes go in a separate file.

**Question.** Reviewers C2Bd (W2, Q3) and h8pC (W2) object that the radius assumes a correct equivalence class. (a) When discovery misses one edge, is the damage to the committed adjustment set confined near the query? (b) Can one extra edge, oriented as the true order allows, make the committed set invalid in the true DAG?

**Population.** The generator of the earlier skeleton study (X3, 2026-08-31), with fresh seeds: SCM seeds 1,000,000 + s (X3 used 0 to 15,999), master seed 20261106 for the size and density draws. Ensembles original (p 5 to 8), licensed (p 5 to 10), large (p 15 to 25); 12,000, 16,000 and 4,800 draws. Linear-Gaussian iSCM, population covariance, so no sampling error. Knowledge K = true orientations of 4 undirected CPDAG edges.

**Eligibility (as X3 after its amendment).** G0 = Meek(CPDAG, K) exists, O*(G0) is defined, and it is valid in the true DAG D. Dropped draws are counted by reason. Validity is the generalized adjustment criterion (gac), decided by two independent implementations; any disagreement stops the run.

**Deletion arm.** Every edge of D deleted in turn; the analyst graph is Meek(CPDAG(D minus the edge), the surviving K). O* is recomputed and classified against D as unchanged, changed but valid, silent bias (invalid in D), not amenable (loud), or path blowup. Hop = distance of the deleted edge to {x, y} on the true CPDAG skeleton.

**Registered statistic 1.** r2 = silent-bias count / (unchanged + changed valid + silent bias) over deletions at hop 2 or more, arm with knowledge, pooled over the three ensembles. Interval: 95% SCM-level cluster bootstrap (B = 4000, seed 20261106), Wilson shown beside it.

**Prediction 1.** r2 near the X3 value 0.000423, upper bound below 0.002. Verdict: SUPPORTED if the bootstrap upper bound is below 0.002; AGAINST US if the point estimate is 0.002 or more; INCONCLUSIVE otherwise. Against us, the authors say that a local skeleton-stability check is not enough and that discovery must be checked beyond two hops; the Limitations paragraph is rewritten to that.

**Addition arm.** Every pair of nodes absent in D, with every orientation that keeps D plus the edge acyclic (one orientation when one node is an ancestor of the other, both when neither is). The analyst graph is Meek(CPDAG(D plus the edge), K). O* is classified against the true D as above.

**Registered statistic 2.** F = number of additions, arm with knowledge, where O* is defined and invalid in D.

**Prediction 2.** F = 0. The argument (not yet a proof): D plus the edge lies in the analyst's class, so O* is valid there, and deleting an edge cannot open a path the criterion needs blocked. One failure, F of 1 or more, refutes the claim that extra edges with true orientations are harmless; the authors then withdraw that sentence and report the failing graph.

**Secondary, reported, not tested.** Hop 0 and hop 1 rates; the arm without knowledge (where the CPDAG alone is amenable); per-ensemble rates; for additions, the loud and efficiency-only rates and the share of committed sets that remain valid in D plus the edge.

**Exclusions.** Path blowup is counted and excluded from denominators. Rows with an oracle disagreement or a theory anomaly (Meek failure or O* invalid in D plus the edge in the addition arm) are reported; any such row means the run is fixed and redone from scratch before any statistic is read.

**Stopping rule.** Every drawn unit is run once. No interim look at either statistic, no extra draws after looking. A crash is resumed, not restarted. If fewer than 5,000 deletions fall at hop 2 or more, statistic 1 is reported as underpowered.

## Deviations from the plan
- The plan says "inserts one absent edge along a topological order". We test every absent pair and every acyclic orientation, which covers every topological order of D. The smoke test caps additions at 12 per SCM; the full run caps nothing.
- Second validity oracle: `bkrobust.gac.is_gac_valid_dag` from the submitted code, beside the vendored X3 criterion; on every 25th SCM, `is_gac_valid_mpdag` also checks that O* is valid on the analyst graph.
