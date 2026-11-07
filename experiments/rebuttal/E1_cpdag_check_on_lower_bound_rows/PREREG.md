# E1 preregistration: the validity check at the CPDAG on the 372 lower-bound rows (2026-11-06)

Written before any run; not to be edited after the first full run.

**Question.** Of the 372 released radii reported only as lower bounds (`release/per_query_records.csv`, `hop_method = local_up_fast_budget`), how many does the validity check at the CPDAG (Algorithm 1, line 2; `method = top_state` in the submitted `breakdown_radius`) settle as infinite?

**Registered statistic.** S = number of rows whose committed set Z is valid at the CPDAG under the generalized adjustment criterion (`is_gac_valid_mpdag(cpdag, x, y, Z)`), as S / N with a Clopper-Pearson 95% interval, N = rows rebuilt exactly (hashes of CPDAG, DAG and knowledge, |Z|, verdict and edge counts equal to the released row). Rows that fail the rebuild are excluded and counted; 372 is printed beside N.

**Secondary.** The same count under the back-door criterion; rows where the two criteria disagree at the CPDAG; on rows not settled, the outcome of a search with ten times the released closure budget (finite radius found, bound raised, unchanged); whether the released search settings reproduce the released `hop_lower_bound`, and under which criterion.

**Prediction.** More than half of the rows settle as infinite under the GAC; zero hard consistency failures.

**Against the authors.** (a) A hard failure: a row settled as infinite on which a search (released or ours) found a failure; a back-door-valid set that is GAC-invalid; the submitted `breakdown_radius` disagreeing with the direct predicate; an exhaustive search finding no failure on a row that fails at the CPDAG. The run exits non-zero and the authors report the bug and its fix before any number. (b) Few rows settle (S < N/2): the rest are reported as lower bounds with their values, raised where the larger search raises them, and the check is said not to close this population.

**Population.** Exactly these 372 rows, all arms; D_RAND knowledge rebuilt from seed 20260913 and the replicate index as in the ledger sweep. No other seeds. The CPDAG is the true DAG's, as released.

**Stopping rule.** One pass over the 372 rows with the registered budgets: release settings (4000 closures; depth = number of retractable edges if at most 12, else 3) and the larger search (40000 closures, 300 s per search). No other budget enters the registered statistic.

## Deviations from the plan

1. The submitted `breakdown_radius` is called only on rows its top-state check answers: past it, the E1 ladder builds an O(n^4) model that does not return on these 64-413-node graphs (the released sweep turned it off above 60 nodes). The larger-budget rerun is an upward BFS over the submitted `local_up_covers_fast`, with closure and wall-clock budgets (`search.py`).
2. `committed_adjustment_set` is absent from the submitted tree; it is vendored unchanged (`vendor_optimal.py`).
3. `search.py` certifies r >= d + 1 with d the smallest depth left unexpanded. The working-tree search behind the release used the largest popped depth, which overstates the bound by one when the closure budget stops mid-shell; both are recorded per row.
4. The working-tree `breakdown_radius` that produced the release decides validity with the back-door predicate, not the GAC; the released bound is compared with both.
