# E3 preregistration: depth budget b and SAT mode (2026-11-06; frozen at first run)

**Question.** Which depth budget b of the upward search minimises total time of the radius computation, and is one CP-SAT minimise call faster than the ladder of feasibility calls? (Reviewer qQYf, W1.)

**Population.** The 144 rows of `results/axisb3/highk.jsonl` (dense Erdos-Renyi, n = 8..18, edge probability 0.6/0.75/0.9). Each instance is rebuilt from (n, edge_prob, seed): `erdos_renyi_dag` with `default_rng(seed)`, nodes relabelled to the unpadded names the file uses, CPDAG, G0 = Meek(CPDAG, knowledge_to_recover). Z = O-set of G0 at the recorded (x, y). Criterion: GAC everywhere. The CPDAG (top-state) check is on in every arm.

**Arms.** b in {0, 1, 2, 3, 4, 6} crossed with SAT mode {ladder (submitted `radius_e1`), minimise (one CP-SAT call minimising the retraction count over the same E1 model)}: 12 arms per instance, plus one reference call to the submitted `bkrobust.hybrid.breakdown_radius` at b = 3. Arm order rotates per instance; caches are cleared before every arm. CP-SAT: 1 worker, seed 0, 300 s per call. The search leg of each arm has a 600 s wall-clock cap; the SAT leg is bounded by the 300 s per-call limit.

**Registered statistics.**
1. Consistency: the number of instances on which two completed exact arms (or the reference) return different radii, or a SAT witness does not fail on replay. Prediction: 0.
2. Budget: T(b) = summed total time over instances, ladder mode. Report b* = argmin T(b) and the ratio T(b)/T(3) for every b with a paired bootstrap 95% interval (10,000 instance resamples, seed 0). Prediction: b* in {1, 2, 3}, and every ratio for b in {1, 2, 4} within [0.5, 2].
3. SAT mode: on instances where the SAT stage runs at b = 0 (finite radius, Z invalid at the CPDAG), the geometric mean over instances of (minimise SAT time / ladder SAT time), with a bootstrap 95% interval. Prediction: interval contains 1 or lies above it (minimise is not faster), because most radii are 1 or 2 and the ladder stops on the first or second rung.

**What would go against the authors, and what they say then.** Any disagreement in (1) is a bug, reported as such, and the radius claims in Sections 4-5 are re-checked before the response is sent. If T(b) barely moves with b, we say b affects running time only, as Section 5 states. If minimise is faster with an interval wholly below 1, we report it, make minimise the default, and note that an OPTIMAL status certifies shells 0..r-1 clean just as the UNSAT rungs do.

**Exclusions.** An instance whose recorded (x, y) gives no defined, non-empty, GAC-valid O-set at G0 uses the first sorted-permutation pair that does (marked `query_source = fallback`). An instance with no valid pair at all is written as an error row and excluded from (2) and (3). Capped arms are kept at the cap as lower bounds and counted; they are never dropped.

**Stopping rule.** One pass over all 144 instances. No reruns on the basis of results; a crashed shard is resumed, which skips completed instances.

**Deviations from the plan.**
- The script that generated `highk.jsonl` is not in any repository, and its node names are unpadded while the current generator pads them. The rebuilt graphs may differ from the App. D ones; each row records whether m_undirected, k_g0, largest component and |Z| match the file. The statistics are computed on the rebuilt instances regardless.
- The App. D timings used the back-door criterion; this run uses GAC, so radii are compared across arms, not with the file.
- The minimise mode is written in this folder (`sat_minimise.py`) on top of the submitted E1 model builder, not in a src/ worktree.
- The 600 s cap on the search leg is added so that b = 6 cannot hang the run on an instance of large radius.
