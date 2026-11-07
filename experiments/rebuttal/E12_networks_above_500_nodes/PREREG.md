# E12 preregistration: radius and audit on the five networks above 500 nodes (2026-11-06)

Written before any run; not to be edited after the first full run.

**Question.** Can the radius be computed on link, munin, munin2, munin3 and munin4 (724 to 1041 nodes), which the released ledger excluded with its 500-node cap (App. I, line 2118), and what does the audit at the true DAG show there?

**Population.** The 100 frozen frame rows of these five networks (`results/frame/frame.jsonl`, 20 per network), crossed with the released conditions: the ten LLM knowledge states of `release/knowledge_states.json`, D_RAND (3 replicates, seed 20260913, rebuilt as in `ledger_sweep.py`), D_DEGEN, A_TRUE and A_ORACLE. About 1,600 units. The CPDAG is the true DAG's, as released. No other seeds, no other networks.

**Procedure per unit.** Build G0 = Meek(CPDAG, K) and the committed set Z on G0 (released code, vendored). Then, observables only: (1) the validity check of Z at the CPDAG under the GAC; (2) if Z is valid there, the radius is infinite and the submitted `breakdown_radius` is called to confirm `top_state`; (3) if not, a budgeted copy of the hybrid's upward search (same covers, same order; depth rule of the ledger; 40,000 closures, 120 s), then the submitted `breakdown_radius` under a hard wall clock (60 s, or twice the search time plus 30 s when the search found the radius) and a 2 GB data-memory cap. Then the claim radius and phi_1 under 120 s. Only after that, the audit at the true DAG (validity of Z by the networkx GAC, the four audit buckets). Time per call is recorded for every leg.

**Registered statistic.** Two numbers. N_net = how many of the 5 networks have at least one certifiable unit (a committed set exists) whose hop radius is exact (infinite by the CPDAG check, or a finite radius found). C = censored fraction: certifiable units with no exact hop radius, divided by certifiable units, with a Clopper-Pearson 95% interval. Reported per network and per condition beside it.

**Secondary.** Fraction of certifiable units whose Z is valid at the truth, by condition group; bucket counts (dangerous, held, safe, slack) for the hop and claim instruments; median and maximum time of the CPDAG check, the search and the hybrid call, by network; how many hybrid calls completed, timed out or hit the memory cap.

**Prediction.** N_net >= 4 and C < 0.10, because most components in these networks are small (link rows have components of size 2) and the CPDAG check settles most rows in one call. Zero consistency failures.

**Against the authors.** (a) N_net < 2 or C > 0.5: the reply states that the evidence covers networks up to 413 nodes and that the radius was not computable above 500 nodes within these budgets, with the censoring table. (b) A consistency failure (CPDAG check and hybrid disagree; search and hybrid return different exact radii or witnesses; a search exhausts the up-set without failure while Z fails at the CPDAG; Z invalid at G0): the run exits non-zero and the bug is reported and fixed before any number is used.

**Exclusions.** None by result. Units whose knowledge is inconsistent with the CPDAG, or with no committed set, are counted and kept out of the denominator of C, as in the ledger. A_ORACLE knowledge is built by the released greedy routine under a 900 s limit per network; a network where it does not finish has its A_ORACLE units recorded as `oracle_knowledge_timeout` and counted.

**Stopping rule.** One pass over all units with the budgets above. No rerun with larger budgets enters the registered statistic; a larger-budget pass, if any, is reported separately as exploratory.

## Deviations from the plan

1. The plan says to run `ledger_sweep.py` with the cap lifted. That script imports two modules absent from the submitted tree and calls `breakdown_radius` with options (`use_ladder`, `search_closure_budget`) the submitted signature does not have. E12 reimplements its per-row logic on the submitted tree, with its scoring functions copied unchanged (`vendor_ledger.py`, `vendor_optimal.py`, `vendor_audit.py`).
2. The submitted hybrid has no closure or time budget, and its SAT ladder builds an O(n^4) model; it is therefore run under an external wall clock and memory cap, and the budgeted search supplies the anytime lower bound. The released sweep had a 5 s ladder limit and no ladder above 60 nodes.
3. The claim radius has a 120 s wall clock that the ledger did not have; on timeout it certifies only r_claim >= 1.
