# E9 deviations from PREREG.md

## 2026-11-06: arth150 is excluded, because its state has no DAG to simulate from

**What happened.** The first full run on the workstation stopped with exit 3 before any PC run: 38 rows reported "frozen certifiable row has no rebuilt plan", all of them arth150 (D_LLM 7, D_RAND 5 + 6 + 10, D_DEGEN 5, A_TRUE 5). Its output is kept on the workstation in `results/rebuttal/E9_failed_run1/`.

**Cause.** The rows themselves rebuild exactly as in E6. What failed is the step E6 does not have: the simulated truth. The arth150 reference CPDAG of the panel (`_data/pest/chat_e.json`, sha256 88b0c645ffbeee7e, 86 arcs, 41 undirected edges) is Meek-closed and acyclic, but one of its undirected components is not chordal: the four undirected edges 422-714, 714-47, 47-736, 736-422 form a chordless 4-cycle. PC on a finite sample can return such a graph; it is not a valid CPDAG, and it has no consistent DAG extension, since every orientation of a chordless 4-cycle adds a v-structure or a cycle. `dor_tarsi` correctly returns None on it. None of the knowledge sets of the panel touches those four edges, so every arth150 state (D_LLM, D_RAND reps 0-2, D_DEGEN, A_TRUE) keeps the 4-cycle and none has an extension either (A_ORACLE has no certifiable arth150 row).

In the old `run.py` the network-level check on D_0 did `continue` before the network's rows were visited, so the rows were never marked as seen and came out as "no rebuilt plan". The real reason, "the reference CPDAG has no consistent DAG extension", went into `model_problems`, which the driver only printed after the exit-3 branch, so it never appeared. The earlier smoke run passed because it uses ecoli70 and magic-niab only.

**What changed.**

1. `build_reference` now rebuilds and checks every frozen row (K hash, committed set and verdict, GAC at the truth) whether or not the row can be simulated. A row whose network has no D_0, or whose own state has no consistent extension, is no longer silently skipped: it goes into an `excluded` list with its reason, and its group is not simulated. A network with no simulated group gets no units, the identity unit included.
2. `reference.jsonl` still has one line per frozen certifiable row (244); excluded rows carry an `excluded` field with the reason, the others `excluded: null`. A new `exclusions.json` gives the counts and the excluded row ids.
3. `model_problems.txt` is written before the exit-3 branch, so a fitted-model problem can no longer be hidden by a rebuild problem.
4. `analyse.py` scores only non-excluded rows and reports `rows_excluded`, `rows_excluded_invalid`, `excluded_reasons` and `networks_excluded`; an excluded network is not listed as "under 160 simulations".
5. `run.py --check-reference` runs only the rebuild and the model fit (about 2 s), for checking the population before a launch.

**Effect on the registered population.** 244 rows are rebuilt and match the frozen panel; 206 are simulated and 38 (all arth150, 35 of them invalid at the truth) are excluded. The registered AUC(E9) and Delta are therefore computed on 206 rows, not 244 (Delta was already defined on the rows present in both experiments). The response must say this, and E6's AUC on the same 206 rows (`auc_e6_paired`) is the comparison, not E6's AUC on 244.

**Alternative not taken.** Simulating arth150 from a DAG that breaks the 4-cycle (a non-consistent extension that adds one v-structure) would keep the 38 rows, but that DAG is not a member of the analyst's state, its fitted law differs from one choice to another, and PREREG.md's model (a consistent extension of G0, the same law for every group) does not cover it. That would be a new design decision for the authors, not a repair.
