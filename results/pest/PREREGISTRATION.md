# Pre-registration — the estimated-CPDAG panel (P-EST, Table 7)

**Written 12 September, after the ledger of `results/ledger/` was complete and read,
after the discovery backend was installed and timed on two networks for wall-clock
only, and before any discovery run of this panel, any structure statistic, any
frame row, any supplier answer or any sweep row of it had been produced or looked
at.** Nothing below is edited after `results/pest/structure.csv` is written. Outcomes
go in Appendix A.

---

## 1. What the panel is

Every row of the ledger was computed on the CPDAG of the true DAG and stamped
`via_Chat_only`. This panel lifts that idealisation on the four networks that carry
real Gaussian parameters (`ecoli70`, `arth150`, `magic-niab`, `magic-irri`): the
CPDAG is estimated from simulated data, and the frame, the questionnaire, the
knowledge and every instrument are computed from that estimate. The true DAG draws
the samples, as nature would, and audits the frozen rows, and is on no other path.
The deployment rows therefore carry `true_dag_on_path = "false"`, because
`chat_source_sha256` differs from `generating_dag_sha256`; the two audit arms,
whose knowledge is read off the truth, keep `via_K`.

## 2. Fixed before the run

- **Data.** Ancestral samples from the linear-Gaussian SEM with the fitted
  coefficients and residual variances of each network's file, at
  n in {200, 1000, 5000}, seeds 0..9, generator seeded by (panel seed, a stable
  hash of the network name, n, seed).
- **Discovery.** PC with the Fisher-z test at alpha = 0.01, stable skeleton, the
  package's default collider handling; GES with the BIC score. Package
  `causal-learn`, version recorded in `structure_summary.json`. A learned graph is
  converted to the repository's MPDAG; a conflict edge (both endpoints arrowheads)
  is recorded and carried as undirected; whether the graph arrives Meek-closed is
  recorded, and the Meek closure is what every later stage uses.
- **Structure statistics per run.** Undirected fraction, chain-component partition,
  skeleton precision and recall against the truth, orientation accuracy on the
  learned directed edges that exist in the truth, and, per (network, n, method),
  the mean pairwise adjusted Rand index of the chain-component partitions across
  seeds and against the oracle CPDAG's partition.
- **Reference run.** PC at n = 5000, seed 0, is C-hat_E for each network. If its
  closure fails the next seed is taken and the fallback is printed.
- **Compute bound, declared here.** Every discovery run has a 900 s wall; a run
  that hits it is a `timeout` row and is counted, not dropped. GES on `arth150`
  (107 nodes) is run on seeds 0..2 only at every n; every other cell runs the ten
  seeds. PC, which supplies the reference, is not restricted. The bound was set
  from one timing: PC on `arth150` at n = 5000 took 208 s and GES had not
  returned after four minutes when this was written.
- **Frame.** `frame_build`'s rules on C-hat_E: treatment in or beside a chain
  component of size at least two, outcome a possible descendant, cap 20 per
  network by seeded stride, inclusion probability recorded.
- **Questionnaire.** `questionnaire_build`'s renderer and domain dictionary on
  C-hat_E's components, same cap of twenty nodes, same control block, per-network
  seed from a stable hash so the prompts reproduce across processes.
- **Supplier.** The primary local supplier only, `qwen2.5:7b-instruct`, real
  names, order seed 0, temperature zero, cached; Meek repair as `elicit_run` does.
  If the host is unreachable the prompts are exported and the panel stops there.
- **Arms.** `D_LLM`, `D_RAND` (3 reps on the model's edges), `D_DEGEN`, `A_TRUE`
  (truthful orientation of the model's edges that exist in the truth; a claim on an
  edge absent from the truth has no truthful orientation and is dropped, counted),
  `A_ORACLE` (truthful orientation of every undirected edge of C-hat_E that exists
  in the truth).
- **Instruments.** Exactly `ledger_sweep`'s: the committed set, `r_claim` to depth 3
  under the 300-subset budget, the hop radius with its dispatch leg, `phi_1`, the
  free rule, closure size and claim count. Instrument rows for every arm are written
  and hashed before the audit runs.
- **Audit.** Validity at the truth by the generalised adjustment criterion on the
  generating DAG (`benchmarks.audit`); the four buckets per instrument; commission
  and omission. Two statuses the oracle panel cannot have, recorded on every row
  and never used as filters: `claim_on_absent_edge` when an asserted claim names an
  edge the true DAG does not contain, and `y_not_true_descendant` when the outcome
  is not a descendant of the treatment in the truth.
- **A-CAL-2.** On every certifiable row, the OLS coefficient of the treatment in
  the regression of the outcome on the treatment and the committed set, computed
  on the reference sample (n = 5000, seed 0, the sample C-hat_E was estimated
  from), against the true total effect of the SEM; the realised absolute bias, its
  Frisch-Waugh-Lovell standard error and the 95 percent band, reported
  conditional on the claim certificate's certified floor.

## 3. Predictions

**P1. Estimation leaves more open than the oracle, and less as n grows.** At
n = 200 the mean undirected fraction of the estimated CPDAG exceeds the oracle
CPDAG's on every network, for both methods; the mean undirected fraction is
non-increasing in n on every network, for both methods. Falsified by one network
where the n = 200 fraction is at or below the oracle's, or by one where the
n = 5000 fraction exceeds the n = 200 fraction.

**P2. The certifiable share falls.** On `D_LLM`, the fraction of frame rows with a
certificate (`status = ok`) on this panel is below the fraction on the oracle panel
restricted to the same four networks (`results/ledger/rows.csv`, 26 of 80 rows,
0.325). Falsified if the panel's share is at or above 0.325.

**P3. The claim certificate stays dangerous on no row where the truth extends the
analyst's graph.** On every arm, among certifiable rows whose asserted claims all
name edges the truth contains and all agree with it in direction where they
exist as claims of the audit (that is, rows whose `truth_status` is `ok`), the
`B7_claim` bucket is `dangerous` on zero rows. Falsified by one such row. Rows
with a wrong skeleton are reported beside, on their own line, with whatever
bucket they get: the structural argument of the ledger's L6 assumed the truth
extends the CPDAG, and this panel is the first place that assumption can fail.

**P4. Wrong-skeleton rows are counted, not dropped.** Every frame row appears in
`rows.csv` for every arm, the two truth statuses are non-empty on at least one
network, and the table prints their counts. Falsified by a missing row or a
silent drop.

**P5. Realised bias is larger at the smallest certified floor.** Pooled over the
deployment arms, the mean realised absolute bias on certifiable rows with a
certified floor of 1 exceeds the mean on rows with a certified floor of at least
2, unbounded included. Falsified if the order is reversed or the two are equal.

## 4. Decision rules

- P3 falsified: halt and inspect the row before anything else is read; on this
  panel a dangerous certificate on a correct-skeleton row is either a bug or a
  result, exactly as on the ledger.
- P2 falsified: the estimated CPDAG did not cost the analyst certificates on these
  networks; report it and do not generalise from four networks.
- P5 falsified: the certified floor does not order bias on this corpus; the
  calibration is reported as failed and the certificate is not described as a
  bias proxy anywhere.

## 5. What this panel cannot say

Four networks, one reference run each, one seven-billion-parameter supplier. The
network cluster bootstrap on four clusters is indicative and is labelled so. The
reference CPDAG is one draw at one n; the ten seeds measure how much it would move.
No human or harvested knowledge, as on the ledger.

---

## Appendix A — outcomes

Scored 12 September from `results/pest/TABLE7.md` and `results/pest/summary.json`, which
`experiments/pest_table.py` builds from the frozen rows. Intervals are a four-network cluster
bootstrap and are indicative only; counts are exact.

### A.0 What changed between registration and the table, before any prediction was scored

- **Two references fell back, as the registered rule says.** PC at n = 5000, seed 0, failed
  its Meek closure on `arth150` and on `magic-niab`. The references are seed 1 and seed 5.
  Across the panel 36 of 120 PC runs failed the closure, 24 of them because the returned
  graph already held a directed cycle. No GES run failed it.
- **The GES bound on `arth150` bit harder than registered.** Run alone, one GES fit at
  n = 5000 took 348 s. With six workers sharing the machine, 7 of the 9 registered cells
  hit the 900 s wall: all three seeds at n = 200, and two of three at n = 1000 and at
  n = 5000. Those seven cells were rerun once, on three workers, under the same wall
  (`structure_retry.csv`). Two finished: n = 1000 seed 1 and n = 5000 seed 1. The
  first-pass rows are unchanged and the retry is printed as its own block. GES on
  `arth150` at n = 200 has no completed run.
- **The supplier host answered.** 23 items, 3 of them from cache, 813 s. Nothing was
  exported.
- **The sweep ran twice, and its frozen instrument rows are byte-identical across both runs**
  (`frozen_rows.sha256` = `a4ff90e6…`). The first table fired P3's decision rule. The second
  run added one audit column, `truth_extends_g0`, and per-network error counts for C-hat_E.
  Both were defined **after** that first table was read and after the inspection below, and
  are post hoc. No instrument number moved.
- **A label was corrected.** The first table called every row whose `truth_status` is not
  `ok` a "wrong skeleton" row. That is false: the flag looks at the claims and at the query,
  and the outcome can fail to descend from the treatment on the oracle frame too. The
  label is now `flagged`.
- **The pipeline was checked before the real sweep.** The oracle CPDAG was fed in place of
  C-hat_E, together with the ledger's own knowledge and frame, on these four networks. The
  sweep reproduced all 420 comparable rows of `results/ledger/rows.csv` with no difference
  in any field. The refactored `questionnaire_build.py` reproduced the frozen questionnaire's
  hash under a fixed hash seed. The refactored `elicit_run.py` reproduced the ledger's
  `D_LLM` answers and knowledge from the cache.

### A.1 Verdicts

**P1 NOT SUPPORTED as registered: it holds for PC and reverses for GES.** Mean undirected
fraction at n = 200, 1000, 5000 against the oracle's:

| network | oracle | PC | GES |
|---|---|---|---|
| ecoli70 | 0.357 | 0.854, 0.545, 0.221 | 0.039, 0.182, 0.222 |
| arth150 | 0.300 | 0.441, 0.302, 0.284 | no run, 0.014, 0.030 |
| magic-niab | 0.152 | 0.211, 0.050, 0.060 | 0.040, 0.072, 0.106 |
| magic-irri | 0.118 | 0.383, 0.190, 0.124 | 0.029, 0.052, 0.063 |

- PC is above the oracle at n = 200 on all four networks. It is non-increasing on three;
  on `magic-niab` it rises from 0.050 to 0.060.
- GES is below the oracle at n = 200 on every network with a run, and rises with n on every
  one of them. At n = 200, BIC returns a sparse graph that is nearly fully oriented, with
  skeleton precision 0.47 to 0.67. A low undirected fraction there is a wrong graph
  oriented, not information.
- The reference sample size does not rescue the partition. For PC at n = 5000, the adjusted
  Rand index of chain-component partitions is 0.373, 0.728, 0.605 and 0.744 across seeds,
  and 0.116, 0.332, 0.470 and 0.164 against the oracle. Skeleton precision is 0.94 to 0.97;
  recall is 0.81 to 1.00. The directed edges that exist in the truth point the right way on
  0.55 to 0.72 of them.

**P2 FALSIFIED.** `D_LLM` certifies 38 of 80 frame rows (0.475), against 26 of 80 (0.325)
on the oracle panel.

- Every arm certifies more: `D_DEGEN` 19 against 11, `A_TRUE` 34 against 18, `D_RAND` 116
  of 240 against 66 of 180. `A_ORACLE` ties at 37.
- The estimate orients more than the truth does. It leaves 8, 41, 5 and 11 edges open where
  the oracle leaves 25, 45, 10 and 12, and it compels 12 of 52, 40 of 86, 14 of 62 and 25 of
  90 edges the wrong way (ecoli70, arth150, magic-niab, magic-irri).
- Identifiability here is bought with wrong orientations. The committed set is invalid at
  the truth on 0.737 of `D_LLM`'s certifiable rows against 0.269 on the oracle panel, and on
  0.737 of `D_DEGEN`'s against 0.000.
- As the decision rule says, this is reported and not generalised beyond four networks. On
  an estimated CPDAG, the certifiable share is not a measure of how good the analyst's
  position is.

**P3 FALSIFIED as operationalised, and the premise in its title held on no row.**

- On the registered operationalisation, certifiable rows with `truth_status = ok`, the claim
  certificate is dangerous on 3 `D_LLM`, 9 `D_RAND`, 3 `D_DEGEN`, 8 `A_TRUE` and 10
  `A_ORACLE` rows.
- **The decision rule fired, and the rows were inspected before anything else was read.**
  There were six rows on `D_LLM` and `D_DEGEN`, covering three queries on two networks. No
  bug was found. Each committed set is valid on the analyst's graph and contains a node the
  truth forbids.
  - On `magic-irri`, for the effect on `FT`, C-hat_E compels `HT -> FT` and `YLD -> FT`
    against the truth, so both variables, descendants of the outcome in the truth, enter the
    set.
  - On `magic-irri`, for the effect of `GTEMP`, C-hat_E lacks the edge `GTEMP -- GW`, so
    `GW`, which the truth forbids, enters the set.
  - On `arth150`, C-hat_E lacks three of the treatment's edges, and 11 members of the set
    are forbidden in the truth.
- No retraction of claims reaches either error, so no certificate computed on the CPDAG can
  see it.
- The registered proxy was the wrong operationalisation of the title's premise. It looks
  only at the claims and at the query, and an error in the skeleton or in a compelled edge
  passes it.
- Measured exactly (post hoc), the generating DAG is a consistent extension of the
  analyst's graph on **0 of 244** certifiable rows on every arm. The reason is that every
  reference skeleton differs from the truth, with 12/2, 29/6, 0/1 and 5/4 edges missing/extra (ecoli70, arth150, magic-niab, magic-irri).
  The ledger's L6 argument therefore has no row to apply to on this panel. "Dangerous on no
  row where the truth extends g0" is vacuously true and says nothing.
- What is measured instead:
  - The claim certificate is dangerous on 14 of 38 `D_LLM` rows, 43 of 116 `D_RAND` rows and
    14 of 19 `D_DEGEN` rows. On the same four networks of the oracle panel it is dangerous on
    none.
  - The hop radius is dangerous on 14, 44 and 14. The two instruments differ by one row on
    the whole panel.
  - Discovery error, not claim error, drives the count: `A_TRUE`, which asserts no false
    claim, is dangerous on 23 of 34.
- Per network, `D_LLM` is dangerous on 6 of 8 rows on ecoli70, 5 of 7 on arth150, 0 of 14 on
  magic-niab and 3 of 9 on magic-irri.
- `magic-niab`, whose reference recovers the whole skeleton plus one extra edge, is the only
  network where neither the model arm nor the chance arm is ever dangerous: its 6 and 14
  invalid sets are all detected. Even there, the audit arms are dangerous on 3 and 4 rows,
  because a compelled edge is wrong.

**P4 PARTLY SUPPORTED.**

- Nothing was dropped. There are 80 rows per arm on the 80 frame rows, and 240 for the
  three chance replicates. The 20 `A_ORACLE` rows on `arth150`, whose truthful orientation
  of every open edge is inconsistent with the estimate's compelled edges, are counted as
  `knowledge_inconsistent`.
- `y_not_true_descendant` fires on 55 of 80 frame rows, on all four networks.
- `claim_on_absent_edge` fires on no row. The supplier's 18 claims after repair all name
  edges the truth contains, and the other arms reuse those edges or filter them on the
  truth by construction. The registered "both statuses non-empty" fails for this one.
- Neither status is what the prediction's phrase "rows where the skeleton is wrong" needed.
  The descendant status also fires on the oracle frame of these networks (43 of 80), and the
  claim status cannot see an error away from the claims. The status that answers the phrase
  is `truth_extends_g0`, and it is post hoc.

**P5 SUPPORTED as registered, and not robust.**

- Pooled over the three deployment arms, the mean realised absolute bias is **0.110** on the
  71 certifiable rows with a certified floor of 1, against **0.063** on the 102 rows with a
  floor of at least 2, unbounded included.
- Three reads undo the order:
  - The medians reverse: 0.021 against 0.028.
  - Without `magic-irri` the means reverse: 0.047 on 52 rows against 0.050 on 83. The
    pooled order comes from that network's chance arm, whose 14 rows at floor 1 average
    0.344.
  - Of the eight network-by-arm cells with rows at both floors, four go each way.
- The FWL 95 percent band covers the true effect on 0.634 of floor-1 rows and on 0.431 of
  the rest. On an invalid set the bias is structural, and a sampling band says nothing
  about it. On valid sets at floor 2 or more, the mean bias is 0.004.
- The registered decision rule for a falsified P5 does not formally apply. Its conclusion
  is the right reading anyway: on this panel the certified floor does not order bias, and
  the certificate is not to be described as a bias proxy.

### A.2 Not registered, read off the same table

- **Commission on the estimated CPDAG.** The primary supplier's 18 post-repair claims are
  false at the truth on 11: 2 of 4 on ecoli70, 3 of 5 on arth150, 3 of 3 on magic-niab, 3 of
  6 on magic-irri. On `arth150`, 11 of the 18 claims before repair were false. Meek repair
  dropped 13 claims and kept 3 false ones of 5. A questionnaire drawn from an estimated
  component asks about edges whose neighbourhood is itself wrong, and this is not comparable
  with the ledger's commission on the oracle components.
- **The sentence this panel supports.** The claim certificate is a statement conditional on
  the CPDAG. On an estimated CPDAG, that condition fails on every row of these four networks.
  The certificate then over-certifies at the rate the committed set is invalid for reasons
  outside the claims: 0.37 of certifiable rows for the model and chance arms, 0.74 with no
  knowledge. The hop radius over-certifies on the same rows. A deployment claim for either
  instrument needs a certificate that retracts the discovery output as well as the
  knowledge. This panel measures the size of the gap; it does not close it.
