# E2 preregistration: running time of the final implementation by graph size (2026-11-06)

Written before any run. Not to be edited after the first run; later notes go in a separate file.

**Question.** Under the generalized adjustment criterion (criterion "gac"), how do the CPDAG check, a single validity call, the depth-3 retraction search, the SAT build and the end-to-end hybrid scale with the number of nodes n (8 to 400) and the number k of undirected CPDAG edges? Serves reviewers h8pC-W1 and qQYf-W1.

**Populations.** Family "appD": the Appendix D protocol, Erdos-Renyi n in {8,10,12,14,16,18,20,24}, p in {0.3,0.5,0.7,0.85}, seeds 0-7 (256 units). Family "sparse": Erdos-Renyi with expected degree d in {2,4}, p = d/(n-1), n in {16,24,32,48,64,100,150,200,300,400}, seeds 0-7 (160 units). G0 is the data-generating DAG (full knowledge). An instance whose CPDAG has no undirected edge, or with no admissible query, is redrawn with seed + 10000*a, a = 1..100, as in the submitted e2e_speedup.py gate. Query: the first ordered pair (X, Y) in sorted order whose optimal adjustment set in G0 is non-empty and GAC-valid in G0.

**What is timed, each in a fresh process.** CPDAG construction; the top-state check (one cold GAC call on the CPDAG) and one cold GAC call on G0 (median of 5 cold calls each); the bounded search (depth 3, only when the top check fails, as in the hybrid); the SAT build (always) and the SAT ladder alone (only when the top check fails); the hybrid breakdown_radius(criterion="gac", search_budget=3) end to end; and, when the hybrid reaches the SAT stage, the unbounded search as an independent check of the radius. Cap: 120 s per leg (Appendix D), 4 GB of memory per leg. A leg over the cap or the memory limit is censored, never dropped.

**Registered statistic.** For each family and n: the median end-to-end hybrid time with a distribution-free 95% interval from order statistics (censored runs count as +infinity), the maximum, and the share censored. Headline: the largest n in the sparse family up to which every hybrid run finished within 60 s ("practical size"). Secondary: the same median and maximum per stage, and per bin of k.

**Prediction.** The top check and the single validity call stay below 0.1 s median at every n up to 400. On appD the hybrid median stays below 10 ms, as in the legacy back-door timings. In the sparse family the bounded search and the SAT build grow fast with k; the practical size is at least 64 and below 400.

**Outcome against the authors.** If some hybrid run exceeds 60 s at n <= 48, or the practical size is below 64, exact computation is not practical at moderate size. We then state the size up to which it is practical, and say that beyond it only the anytime lower bounds (no failure within the searched depth, or within the SAT rungs certified unsatisfiable) are available.

**Consistency (run exits non-zero on any failure).** The hybrid radius equals the bounded search radius when that search finds a failure, equals the SAT ladder radius when both are exact, and equals the unbounded search when both finish; the top-check verdict matches the hybrid's top_state dispatch; for n <= 24 the Appendix D G0 construction equals the DAG and the fast optimal-set routine equals the repository's.

**Exclusions and stopping.** Excluded from the statistics, and counted: units whose gate is exhausted, and units that fail for an infrastructure reason. Per series (family and parameter), a leg is skipped at larger n once at least 75% of its attempts (at least 4) were censored at some smaller n; skipped legs are reported as such and counted as censored in the hybrid statistic. The run stops when every unit of the fixed grid has a row. Units with status "error" may be rerun once (--retry-errors); nothing else is rerun, and the cap is not changed after results are seen.

**Deviations from the plan.**
1. Appendix D chose the query with the back-door enumeration oracle; here validity is GAC (polynomial), and the optimal set is computed directly on G0, which is a DAG. For n <= 24 the result is checked against the repository function.
2. Appendix D used p up to 0.85 for every n; for n above 24 a dense Erdos-Renyi graph has almost no undirected edges and tens of thousands of edges, so the large-n family is sparse (constant expected degree). The two families overlap at n = 16 and 24.
3. For n > 24, G0 is taken as the DAG directly instead of Meek(CPDAG, greedy knowledge); both are the same graph (checked at n <= 24), and the greedy knowledge routine is quadratic in k.
4. The space-building leg of Appendix D is not rerun: it is exponential in k and was already shown to exceed 120 s beyond k = 10.
5. Runs are serial (one worker) so that timings are not perturbed by concurrent jobs.
