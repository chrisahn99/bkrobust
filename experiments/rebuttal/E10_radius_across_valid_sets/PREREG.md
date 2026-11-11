# E10 preregistration: the radius across valid adjustment sets (2026-11-06, frozen before the first run)

**Question.** For one query, how much does the breakdown radius vary across the adjustment sets valid in G0, and how often is the committed optimal set O strictly beaten by another valid set? (Reviewer qQYf, Q2 and minor comment 2.)

**Populations.** (1) Published networks: every released row of `per_query_records.csv` with `status = ok` in the ten elicited arms and `D_RAND` (three replicates), state rebuilt from the released files and checked against the row's hashes and counts. `A_TRUE`, `A_ORACLE` (reference arms) and `D_DEGEN` (no knowledge, so every radius is infinite) are excluded. (2) Census: every CPDAG on 4 nodes with 1 to 6 undirected edges and on 5 nodes with 1 to 4, every knowledge state G0 of its corrected space, every ordered pair (x, y).

**Sets.** All GAC-valid sets in G0, exhaustively, when at most 12 nodes lie outside the forbidden set (at most 512 kept, O always among them); otherwise O, the canonical set and valid sets found by a seeded local walk, at most 20 per query. Census: exhaustive.

**Radius.** Criterion `gac` throughout. Networks of at most 60 nodes: the submitted `breakdown_radius` (search budget = retractable edges if at most 12, else 3; ladder 60 s per rung). Larger networks: the same top-state check, then the submitted cover generator in an upward search with no depth cap, 4,000 closures and 20 s per set; a stopped search gives a lower bound (censored). Validity at the true DAG for every set, by two independent implementations. Census: hybrid radius and the exact BFS radius on the enumerated space.

**Registered statistic.** Among network queries with at least one non-empty valid set besides O, the fraction where some listed set is certified strictly more robust than O (infinite beats every finite radius; a lower bound counts only when it exceeds O's exact radius), over queries whose comparison is decided. 95% Clopper-Pearson interval, plus a network-clustered bootstrap interval (2,000 resamples). Also reported: share of queries whose radii differ across sets (spread), the largest spread, the share of censored radii and of undecided queries, the census counts, and the worst case with undecided queries counted as beaten.

**Prediction.** O is never beaten in the census under gac (as under back-door, 0 of 249,732), and on networks the fraction is at most 5% with its upper bound below 5%. The radius differs across sets on a fifth to a third of queries.

**Against the authors.** If O is beaten on more than 5% of network queries with two or more valid sets (point estimate), or on any census instance, we report the most robust set beside O with its variance cost (asymptotic variance ratio under five random linear SEMs on the true DAG) and drop any wording that O is the robust choice. If more than half of the comparisons are undecided, we state that the evidence on large networks is a lower bound.

**Seeds and exclusions.** D_RAND knowledge seed 20260913 (released). Set sampling and SEMs seeded by a hash of the query key. Rows whose rebuild does not match the release are reported, not analysed, and make the run exit non-zero. Empty sets are listed but not counted as alternatives in the registered statistic (secondary line includes them).

**Stopping rule.** One pass over all units; a resume only reruns units that raised. No budget changes after the first full run; any later rerun is reported beside this one.

**Consistency failures (exit 1).** Hybrid and enumerated-space radius differ (census); ladder and search both exact and different; the two truth oracles differ; O's exact radius or truth validity differs from the released row; a rebuild mismatch.

## Deviations from the plan

- Above 60 nodes the submitted `breakdown_radius` would enter its SAT ladder, whose quartic model build does not return at that size; the released sweep skipped the ladder there too. A wrapper keeps the top-state check and replaces the ladder by a budgeted search over the submitted cover generator (`vendor_search.py`).
- "At most 20 sets" is met by a seeded local walk, so on large networks "O never beaten" is about listed sets, not all valid sets.
- The census uses the corrected space (`build_space_correct`); the back-door census used the older enumerator, which misses states at 5 or more undirected edges.
