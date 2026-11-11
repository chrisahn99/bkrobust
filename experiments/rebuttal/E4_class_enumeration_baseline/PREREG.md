# E4 preregistration: enumerating the class as a timing baseline (2026-11-06)

Written before any run. Not to be edited after the first run; later notes go in a separate file.

**Question.** Corollary 1(c) gives the radius as min |K_G0 \ dir(D)| over DAGs D of the class [Ĉ] in which Z is not a valid adjustment set. For which class sizes |[Ĉ]| is listing [Ĉ] and taking that minimum faster than the hybrid (`bkrobust.hybrid.breakdown_radius`)? Serves reviewer points qQYf-W1 and h8pC-W1.

**Populations.** (1) four-node: every CPDAG on 4 nodes with an undirected edge, every G0 of its corrected space, every ordered (X, Y) with a non-empty optimal set Z valid in G0 (back-door filter, the selection behind the 1,044 instances of App. D); the count is reported, not assumed. (2) highk: the 144 rows of `results/axisb3/highk.jsonl`, rebuilt with `dense_instance` logic, Z = optimal set of G0. (3) timing: the 256 rows of `results/e2e_speedup_gated_v2/e2e_speedup.jsonl`, rebuilt from (n, effective_seed, edge_prob) with the stored Z. No random draws; every instance is fixed by its stored parameters.

**Methods, per instance, in one process, sequentially.** hybrid: `breakdown_radius(..., criterion="gac")`, default budget b = 3. enum: depth-first listing of [Ĉ] by orienting one undirected edge at a time and closing with Meek's rules, GAC checked in every DAG, cap 10^5 DAGs and 120 s. bnb (secondary): the same listing with a branch-and-bound prune on the number of K_G0 orientations already reversed. lib (check only, k ≤ 14): the brute-force `enumerate_dag_extensions`. Class size is counted exactly by root-picking (He, Jia and Yu 2015), independently of the listing.

**Registered statistic.** For each population and each class-size bin (1, 2-10, 11-100, 101-10^3, 10^3-10^4, 10^4-10^5, >10^5), the fraction of instances on which enum finishes (not capped, not timed out) with a wall time below the hybrid's, with a Wilson 95% interval. Pooled over populations: the crossover bin, the lowest bin whose fraction is below 0.5. Secondary: the same for bnb, and median time ratio enum/hybrid per bin with a bootstrap 95% interval (2,000 resamples, seed 0).

**Prediction.** enum is faster only on very small classes: fraction above 0.5 for class size 1 to 10, below 0.5 from 11-100 upward; below 0.5 overall on highk and timing; the crossover bin is 11-100 or lower.

**Against the authors.** enum faster on a majority of the timing population (the realistic one). The authors then say so in the response, present the retraction search and the SAT stage as needed only for large classes, and report the crossover class size.

**Consistency (any failure is a bug, the run exits non-zero).** Complete enum radius equals the hybrid radius (∞ when no DAG fails). Complete enum count equals the counted class size. A capped enum that found a failing DAG gives an upper bound: hybrid ≤ it. bnb when complete equals the hybrid. lib radius and DAG set equal enum's. The first 200 DAGs listed are consistent extensions of Ĉ. A hybrid answer with exact = False is flagged and left out of equality checks.

**Exclusions.** Instances whose rebuilt graph or Z does not match the stored fields (size of Z, number of undirected edges) are excluded and counted. Instances whose class could not be counted within 120 s go to an "unknown" bin. Nothing else is excluded; capped or timed-out enum counts as "not faster".

**Stopping rule.** One pass over every instance of the three populations. No re-runs to change a number; a crash is resumed, never restarted from scratch.

**Timing caveat.** Times are wall clock inside a worker; up to 8 workers run concurrently on the workstation, so absolute times carry contention. Both methods of an instance run in the same process back to back, so the comparison is within instance.

## Deviations from the plan

- The plan names `enumerate_dag_extensions` (meek.py line 362). It tries all 2^k orientations whatever its `limit`, which cannot run at k up to 51 (highk). The registered enum leg is a depth-first listing with Meek closure; the library function is run as a cross-check where k ≤ 14.
- The plan says "the 1,044 four-node instances". The script that selected them is not in the repository; the selection is reconstructed from the App. D description and the hybrid test, and the count it gives is reported.
- The highk rows store |Z| but not Z; Z is recomputed as the optimal set of G0 and checked against the stored size.
