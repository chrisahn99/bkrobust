# Pre-registration — the knowledge interval and the null radius (measurement M1)

**Written 12 September, after the ledger of `results/ledger/`, the frozen frame and
the elicited knowledge were complete and read, after the effect-set enumeration
passed its test against brute-force DAG extensions
(`tests/estimation/test_knowledge_interval.py`), and before
`experiments/interval_panel.py` had produced any row.** The primitives were timed
on four networks for wall-clock only (closure, committed set, one claim-radius
search). Nothing in sections 1 to 6 is edited after `results/interval/rows.csv`
is written. Outcomes go in Appendix A.

---

## 1. The object

The analyst's claims `K` orient the CPDAG `C` of the true DAG into
`G0 = Meek(C, K)`. The retraction ball `B_r(K)` is every Meek re-closure of `K`
minus `S` for `|S| <= r`. For each graph in the ball, the set of total effects of
`X` on `Y` it leaves possible is computed IDA-style: every parent set of `X` that
some DAG extension of the graph realises, found semi-locally (orient the
undirected edges incident to `X` in every way `apply_orientations` accepts, trying
only cliques of the undirected neighbours), and for each parent set `P` the
coefficient of `X` in the regression of `Y` on `X` and `P` (zero when `Y` is in
`P`). The knowledge interval `I_r` is the hull over the ball; from a sample
covariance the band is the hull of the per-effect 95 % confidence intervals
(classical OLS standard error with an intercept, `n - |P| - 2` degrees of freedom,
normal quantile). The class-wide interval is the same over `C`.

Three exact reductions are used and were tested against brute force before this
file was written (random DAGs on 5 to 8 nodes, claims partly reversed, 7,770
treatment and outcome pairs, 2,035 with more than one parent set):

- only the treatment's chain component of `C` enters the parent sets, so the ball
  is enumerated over the claims inside that component (`K_loc`), and claims
  elsewhere leave every interval unchanged;
- a local closure is inconsistent iff the full closure is;
- the parent sets reachable within budget `r` contain those within `r - 1`
  (retracting more never removes an extension), so the hull at budget `r` is
  non-decreasing in `r`.

**The null radius `r0`** is the smallest `r` in 0..3 at which the interval
contains zero, with the zero test `lo <= tol and hi >= -tol`,
`tol = 1e-6 * sd(Y) / sd(X)` from the covariance in use (the fix to the
`has_zero` knife-edge found in the earlier real-coefficient run). It is computed
three ways: from the population interval, from the band, and from the plug-in
hull of the sample point estimates. Status: `exact` when found; `never` when the
ball was enumerated to `|K_loc| <= 3` and no budget contains zero (retracting all
local claims is the class itself); `>3` when `|K_loc| > 3` and nothing up to 3
contains zero; `censored_at_depth_d` when the subset budget stopped the ball
first. Censored and `never` values enter comparisons as 4.

**The validity radius `r_val`** is the ledger's claim radius (`r_claim`, depth 3,
300-subset budget, generalised adjustment criterion on each re-closure, over the
network-level claim set). It is read from `results/ledger/rows.csv` for `D_LLM`,
`D_LLM_72B`, `A_ORACLE` and the truthful arm (`A_TRUE`, which is the controlled
condition at `b = 0`), and computed with `ledger_sweep.breaking_depths` on the
same inputs for `b = 1, 2`.

## 2. Populations, parameters and data

- **Panel B**, real coefficients: `ecoli70`, `arth150`, `magic-niab`,
  `magic-irri`, with the fitted coefficients and residual variances of their
  files. The edge set is asserted equal to the parsed network's.
- **Panel A**, semi-synthetic: the other 23 networks of `results/ledger/rows.csv`
  (all at most 500 nodes). Coefficients uniform in magnitude on [0.5, 1.5] with a
  random sign, drawn in sorted edge order; unit noise variance. The generator is
  `numpy.random.default_rng` seeded with the first 16 hex digits of
  `sha256("interval-sem-v1:<network>")`. Labelled semi-synthetic wherever shown.
- **Queries**: the frozen frame rows of those networks (528 rows; 448 in A, 80 in
  B). Table 2 adds the author-declared exposure and outcome of the applied DAGs
  (`results/frame/declared_rows.csv`, status `ok`), which are not frame rows and do
  not enter Table 1.
- **Truth**: the total effect by path sum, `(I - B)^-1`.
- **Data**: population covariance (no band), and ancestral samples of every node
  at `n` in {1000, 20000}, seeds 0..19, generator seeded by
  `sha256("interval-data-v1:<network>:<n>:<seed>")`; sample covariance with
  divisor `n - 1`.
- **Bounds**: at most 4,096 clique candidates per graph (a censored graph is a
  censored row, recorded); at most 300 non-empty retraction subsets per row and
  condition, cumulated over depths 1..3 as in the ledger.

## 3. Knowledge conditions

- **Controlled, `CTRL_b` for `b = 0, 1, 2`.** The truthful orientations of the
  primary supplier's own edges (the `A_TRUE` construction). **Design choice,
  written before running:** the `b` reversals are drawn per row among that row's
  local claims `K_loc`, uniformly among the `b`-subsets whose reversal the local
  Meek closure accepts, seeded by `sha256("interval-reverse-v1:<network>:<x>:<y>:<b>")`.
  Reason: by the component reduction a reversal outside the treatment's component
  changes no interval, so drawing among all of the network's claims would label
  as `b = 1` rows whose interval saw no wrong claim. Claims outside the component
  stay truthful, so the committed set and `r_val` are computed on the full graph
  with them. A row where `|K_loc| < b` or where no `b`-subset survives the closure
  is recorded as `b_unreachable` with that reason and leaves the `b` column.
- **Elicited, `D_LLM` and `D_LLM_72B`**, as frozen. Each row records the false
  claims inside the treatment's component (`n_false_local`, the stratifier) and
  over the whole network (`n_false_total`, reported beside it). The local count is
  the stratifier for the same reason as above.
- **`A_ORACLE`** (the recovering set) is computed on every row for Table 2 and
  Figure 1 only.

## 4. Procedures and metric

Same estimator in every row of Table 1.

- **PLAIN**: the 95 % interval of the coefficient on the committed set
  `Z = O*(G0)` (`committed_adjustment_set` on the full analyst graph).
  **Design choice:** where `G0` commits to no set (`not_amenable`, or
  `y_not_possible_descendant`, where the graph says the effect is zero), PLAIN is
  the band over `G0`'s own possible parent sets, i.e. the interval at budget zero.
  Reason: excluding those rows would condition Table 1 on the analyst's graph
  being amenable, which drops precisely the rows where a wrong claim orients the
  query away, and a zero-width report of zero cannot be calibrated. A secondary
  table restricted to rows with a committed set is also printed.
- **FAR1**: the band over `B_0` together with the single retractions of local claims
  incident to neither `X` nor `Y`. **NEAR1** (for the locality section, not a
  Table 1 row): the same over claims incident to `X` or `Y`.
- **I1, I2**: the band over `B_1`, `B_2`. **BLANKET**: the class-wide band.
- **CH_ORACLE** (feasibility M2): the Cinelli and Hazlett omitted-variable bound
  around PLAIN on rows with a committed set, `|bias| <= sqrt(R2_Y R2_X / (1 - R2_X))
  * sd(Y|X,Z) / sd(X|Z)`, interval `beta_hat -/+ (bound + 1.96 se_adj)` with
  `se_adj = se * sqrt((1 - R2_Y) / (1 - R2_X))`. The oracle parameter is the pair
  of population partial R-squared values of the omitted block
  `U = V_true` minus `Z`, where `V_true` is the true optimal set when `Y` descends
  from `X`, the true parents of `X` when it does not and `Y` is not a parent of
  `X`, and undefined otherwise. **The parameter is defined** on a row iff `Z` plus
  `U` is a valid adjustment set at the truth: only then is the bias of `Z` the
  omitted-variable bias the bound models. On other rows (a descendant or a
  collider put into `Z`, or no valid set) it is reported as ill-defined, and the
  formula applied there is shown only as the P9 diagnostic.
- **Informative stratum**: rows whose class-wide population interval has width
  above `tol`. It does not depend on `K`. Table 1 conditions on it and prints its
  size.
- **Direction-dispute share.** **Definition fixed here:** the brief's parenthesis
  (the width restricted to parent sets under which `Y` is a possible descendant of
  `X`, over the full class width) is the share that is *not* in dispute; the
  earlier synthetic runs measured it at 0.31 to 0.36 and called the complement,
  about two thirds, the dispute. Both are recorded: `restricted_share` and
  `dispute_share = 1 - restricted_share`. P5 is scored on `dispute_share`.
- **Metric.** Per replicate `i` with interval midpoint `m_i` and half-width
  `h_i`, `lambda_i = |tau - m_i| / h_i` (0 when `h_i = 0` and `tau = m_i`,
  infinite when `h_i = 0` otherwise). Per cell (panel, `b`, `n`) and procedure,
  `lambda*` is the smallest value with at least 95 % of `lambda_i <= lambda*`;
  the calibrated width is `mean(2 lambda* h_i)`, and
  `L95 = calibrated width / calibrated BLANKET width` on the same replicates. Raw
  coverage and the ratio of mean raw widths are printed beneath. A replicate is a
  (row, seed) pair. **L95 is computed at `n = 1000` and `n = 20000` only**: at the
  population most procedures are points and no scaling can calibrate them, so the
  population cells print raw coverage and width ratio only.
- **Dispersion.** Network cluster bootstrap, `B = 2000`, generator seed 20260912,
  percentile 95 % intervals for each L95 and for the ratio `L95(I1) / L95(PLAIN)`.
  Panel B has four clusters and its intervals are printed with that caveat.
- **Column populations.** Primary: every informative row where that column's `b`
  is reachable. Secondary: the rows where `b = 2` is reachable, used for all three
  columns, so the columns share one population.

## 5. Predictions

Each is scored in every (panel, `n` in {1000, 20000}) cell unless stated. A
prediction is **supported** if it holds in all four cells, **partial** if in some,
**falsified** if in none. Orderings are on point estimates; the bootstrap
intervals are reported beside them and do not change the verdict.

**P1. Insurance costs when nothing is wrong.** At `b = 0`, PLAIN has the smallest
L95 of the Table 1 rows and `L95(I1) > L95(PLAIN)`. Falsified in a cell if any
other row's L95 is at most PLAIN's.

**P2. One wrong claim is bought back.** At `b = 1`, PLAIN's raw coverage is below
0.95, `L95(PLAIN) > L95(I1)` and `L95(I1) < L95(BLANKET) = 1`. Falsified in a cell
if any of the three fails.

**P3. Far claims buy nothing.** At `b = 1` and `b = 2`, FAR1's raw coverage is
within two points of PLAIN's. Falsified in a cell if the absolute difference
exceeds 0.02 at either `b`.

**P4. Two wrong claims need two retractions.** At `b = 2`, I2's raw coverage is at
least 0.95 and I1's is below 0.95. Falsified in a cell if either fails.

**P5. The width is mostly a dispute over direction.** The median `dispute_share`
over informative rows, pooled over both panels, is at least 0.5. Falsified below
0.5. Reported per panel as well.

**P6. The informative stratum is a minority.** It holds between 15 % and 40 % of
the 528 frame rows, pooled. Falsified outside that band.

**P7. The plug-in hull is anticonservative and the band is not.** At `n = 1000`,
over every frame row and replicate of the three controlled conditions where `b` is
reachable, pooled over panels, the plug-in `r0` exceeds the population `r0` on at
least 15 % of replicates and the band's `r0` exceeds it on at most 5 %. Falsified
if either fails.

**P8. Elicited knowledge reproduces the ordering.** Pooling `D_LLM` and
`D_LLM_72B` rows, stratified by `n_false_local`: at 0 false claims the P1 ordering
holds, and at 1 false claim the three P2 statements hold. Cells with fewer than 20
informative rows are printed as underpowered and not scored; the prediction is
scored on the remaining cells.

**P9. A confounding bound tuned to its own oracle misses an orientation error.**
On rows of the controlled `b = 1, 2` and elicited conditions whose committed set
contains a true descendant of `X` and is invalid at the truth, the CH_ORACLE
formula's raw coverage at `n = 20000`, pooled over panels, is below 0.95.
Falsified at 0.95 or above; underpowered below 20 such rows. As a check on the
implementation, its population coverage on rows where the parameter is defined
must be 1.

## 6. Table 2 and Figure 1 rules

- **Knowledge for Table 2 and Figure 1**: the recovering set `A_ORACLE`, the one
  the declared-query radius in `results/frame/declared_rows.csv` was computed on,
  so `r_val` there is comparable.
- **Declared queries**: every applied DAG with status `ok` there, with the
  semi-synthetic parameters of Panel A, labelled. Excluded ones are listed with
  their status.
- **Real-coefficient queries**: per Panel B network, among frame rows whose
  `A_ORACLE` graph commits to a set (so there is an estimate to report), finite
  population `r0` first, then frame order; the first five. The brief says "finite
  r0 first, then frame order"; the committed-set filter is added because a row with
  no set has no estimate column.
- **Per row**: `|K|` and `|K_loc|`; the estimate and its 95 % interval at
  `n = 20000`, seed 0; `I1` and `I2` as bands at `n = 20000`, seed 0, with the
  population intervals beside; `r0` from the population, from the band at seed 0,
  and its range over the 20 seeds; `r_val`; the claims whose retraction first
  makes the band at seed 0 contain zero, in the network's variable names: the
  retraction set at budget `r0` whose endpoint crossed zero.
- **Figure 1 candidates**: `A_ORACLE` first, then `CTRL_b0`; real-coefficient
  networks first, then applied DAGs (declared queries, then frame rows). A
  candidate has population `r0` and seed-0 band `r0` both in {2, 3}, a band that
  excludes zero at every budget below `r0`, and a plateau: the band width at
  budget `r0 - 1` at most 1.5 times the width at budget 0. Ranked by the jump
  `width(r0) / width(r0 - 1)`. Up to six. If fewer than six qualify, the nearest
  misses are listed with the failed condition named.

## 7. What this measurement cannot say

The CPDAG is the true one on every row. Estimation is linear-Gaussian with a
correct model class. Panel A's parameters are drawn, not fitted. The frame is
capped at twenty pairs per network, and the informative stratum is where every
width statement lives. Panel B has four networks, so its bootstrap is coarse.

---

## Appendix A — outcomes, 13 September, after the run was read

### A.0 What happened between this file and the tables

- **Before the run.** A crash-only run on two networks (`asia`, `Polzer_2012`) into a
  scratch directory, read for status counts and for the implementation invariants of
  A.1 as violation counts. No metric of section 4 was computed on it.
- **After the first full run was read**, three changes that touch no row and no
  prediction rule. The per-row wall-clock column was dropped so the rows file is
  deterministic. Per-budget profiles, first kept for `A_ORACLE` and `CTRL_b0` only,
  were extended to every condition on the applied DAGs and the real-coefficient
  networks, because the Figure 1 rule of section 6 returned no candidate (A.4). In
  those profiles an empty committed set was first written as missing, which dropped
  rows with an empty set from the Table 2 selection; that was a bug in the profile
  writer, fixed, and Table 2 is from the fixed run. In the table script, bootstrap
  intervals are order statistics, so an infinite upper end is printed as `inf)`
  instead of being dropped, and JSON outputs write non-finite values as strings.
- **Reproducibility.** Five full runs after the timing column was dropped, with 4, 3,
  6, 7 and 7 worker processes, gave a byte-identical `rows.csv` (sha256 `ef7a98d0…`)
  and identical uncompressed replicates (sha256 `c0a9e783…`), the latter also equal
  to the first run's; a run takes about 90 s on 7 processes.

### A.1 Implementation checks (all passed)

On every row with status `ok` (2,749 rows): the class-wide population interval
covers the truth; under truthful knowledge (`CTRL_b0`, `A_ORACLE`) the budget-zero
interval covers the truth; `I0 ⊆ I1 ⊆ I2 ⊆ I3 ⊆ BLANKET` and `FAR1, NEAR1 ⊆ I1`; a
committed set valid at the truth has population coefficient equal to the truth; the
committed verdict equals the ledger's on all 2,112 rows whose `r_val` is read from
the ledger; the realised number of false local claims equals `b` on every controlled
row. Zero violations. The confounding bound at its oracle parameter covers the truth
in the population on every row where the parameter is defined (coverage 1.000). An
independent simulation of twelve real-coefficient rows under the recovering set
(100 datasets of 5,000 each, ordinary least squares on the committed set) gave
coverage 0.90 to 1.00 and absolute bias below 0.004.

### A.2 The predictions

**P1 SUPPORTED** (4 of 4 cells). At `b = 0`, L95 of PLAIN / I1: Panel A 0.464 / 0.957
(`n = 1000`), 0.407 / 0.922 (`n = 20000`); Panel B 0.840 / 0.979, 0.764 / 0.957. The
bootstrap interval of `L95(I1) / L95(PLAIN)` is above 1 in all four cells (Panel A
[1.69, 2.69] and [1.76, 3.52]; Panel B [1.12, 1.28] and [1.16, 1.49]).

**P2 SUPPORTED** (4 of 4). At `b = 1`, PLAIN coverage / L95 against I1's L95: Panel A
0.597 / 5.12 vs 0.851 and 0.573 / 14.37 vs 0.758; Panel B 0.885 / 3.32 vs 0.974 and
0.787 / 8.74 vs 0.952. PLAIN's bootstrap intervals have an infinite upper end in Panel
A, because 3.6 % of its replicates are zero-width reports of a zero effect (the
analyst graph makes the outcome a parent of the treatment) on rows where the effect is
not zero, and some resamples push that share past 5 %.

**P3 PARTIAL** (holds in Panel B, falsified in Panel A). FAR1 minus PLAIN coverage:
Panel A +0.044 and +0.060 (`b = 1, 2`, `n = 1000`), +0.042 and +0.052 (`n = 20000`);
Panel B between -0.006 and +0.005. In the population, where PLAIN equals the budget-zero
interval, the Panel A gap is +0.019 and +0.031, so about half of the finite-sample gap
is far retractions moving the parent sets through Meek propagation inside the component
and the other half is FAR1 containing the parent-set band over `G0`, which is wider than
the committed set's interval. Far claims are not inert on the semi-synthetic panel;
near claims carry almost all of it (NEAR1 coverage 0.960 against FAR1 0.640 at `b = 1`,
`n = 1000`).

**P4 SUPPORTED** (4 of 4). At `b = 2`, I2 / I1 coverage: Panel A 0.984 / 0.781 and
0.981 / 0.757; Panel B 0.979 / 0.947 and 0.960 / 0.909. The Panel B cell at
`n = 1000` holds by 0.003.

**P5 SUPPORTED.** Median `dispute_share` 1.000 over the 460 informative rows (Panel A
1.000, Panel B 1.000); mean 0.630. The restricted width is exactly zero on 58.0 % of
informative rows and its mean share is 0.370, which reproduces the synthetic runs
(60.7 % and 0.314).

**P6 FALSIFIED.** The informative stratum is 460 of 528 frame rows, 87.1 % (Panel A
391 of 448, Panel B 69 of 80), not 15 to 40 %. The earlier fifth was measured on
randomly drawn problems; this frame draws treatments in or beside a chain component
with the outcome a possible descendant, so the class leaves the effect open on almost
every row. The non-informative rows are separated cleanly from the tolerance: their
class width is at most 1.6e-9 times `tol`, the smallest informative width is 9.4
times `tol`.

**P7 SUPPORTED.** At `n = 1000`, over 23,120 replicates of the controlled conditions,
the plug-in hull states a null radius above the population one on 39.7 % and the band
on 1.7 % (Panel A 36.9 % and 1.6 %, Panel B 55.9 % and 2.4 %). At `n = 20000`: 39.0 %
and 1.8 %. On the elicited rows at `n = 1000`: 38.7 % and 1.7 %.

**P8 SUPPORTED** (4 of 4, all cells powered: false = 0 has 450 rows in A and 85 in B,
false = 1 has 208 and 23). At one false local claim, PLAIN coverage / L95 against I1:
Panel A 0.600 / 6.92 vs 0.867 (`n = 1000`) and 0.552 / 15.64 vs 0.854; Panel B 0.696 /
2.41 vs 0.901 and 0.574 / 4.28 vs 0.825.

**P9 SUPPORTED.** On 102 rows over 17 networks whose committed set contains a true
descendant of the treatment and is invalid at the truth, the formula at its oracle
parameter covers 13.4 % of replicates at `n = 20000` and 6.9 % in the population.

### A.3 Feasibility M2, the verdict

The oracle parameter is **well defined on 666 of the 808** committed-set rows under
possibly wrong knowledge (`CTRL_b1`, `CTRL_b2`, both elicited conditions). It is
**ill-defined on 142**: on 101 the committed set contains a true descendant of the
treatment, so no omitted block makes the set valid and the bias is not an
omitted-variable bias; on 41 the truth admits no valid adjustment set for the pair
(the outcome is a parent of the treatment, so the effect is zero and no set blocks the
edge). No row failed because of a collider alone. Where it is defined the bound, at its
oracle, covers and is narrow (Panel A, `b = 1`, `n = 20000`: coverage 0.964, L95 0.445
against I1's 0.632 on the same rows), which is what an oracle buys. A Table 1 row for it
can be printed only on the defined rows, with this count beside it; on the descendant
rows the same formula covers 13.4 %.

### A.4 Table 2 and Figure 1

Table 2 has the nine declared queries with status `ok` and 17 real-coefficient queries:
five each on `ecoli70`, `magic-niab` and `magic-irri`, and on `arth150` the only two
frame rows whose recovering-set graph commits to a set.
Seven of the nine declared queries have no claim inside the treatment's component, so
their interval is a point and `r0` is `never`; `mediator` is the one with a finite
`r0` (2). The rule of section 6 returned **no Figure 1 candidate**: no row under the
recovering set or the truthful controlled knowledge has both radii in {2, 3} on a
real-coefficient network, and the declared and applied-DAG rows that do have them widen
by a factor of 27 to 34 before the break, so the plateau condition fails. A second list, **outside the
pre-registered rule** and labelled so in `fig1_candidates.json`, keeps the radii
condition, drops the plateau condition, admits every knowledge condition and puts rows
whose local claims are all true first; 46 rows qualify. Its first row is a
real-coefficient query whose two local claims are both true: `magic-niab`,
`G311 -> G43`, true effect -0.255, estimate -0.251 [-0.259, -0.244], interval flat at
budget 1 (band width ratio 1.06) and reaching zero at budget 2 when both claims are
retracted, `r0 = 2` on all twenty seeds, `r_val = 2`. The zero it reaches is the
graph in which the outcome becomes a parent of the treatment. The claims come from the
larger elicited supplier; the recovering set asserts only one of them.

### A.5 Censoring and bounds

- Parent-set candidates: no graph reached the 4,096 cap.
- Retraction subsets (300 per row and condition): the ball stopped at depth 3 on 20
  rows of each controlled condition (`hailfinder`, 13 local claims), 20 rows of
  `D_LLM` (`hailfinder`), and 58 rows of `D_LLM_72B` (`insurance` 20, `paths` 20,
  `arth150` 12, `ecoli70` 6). The population `r0` is censored at depth 3 on 11 of
  them; I1 and I2 are exact everywhere.
- Unreachable `b`: `b = 1` on 72 rows (fewer local claims than `b`); `b = 2` on 356 rows
  (288 with too few local claims, 68 with no Meek-consistent pair of reversals).
- `r_val` was read from the ledger for 2,112 rows and computed for 628; `arth150`
  took 95 s including sampling, with no censoring beyond the 12 rows above.
- Panel B has four networks; its bootstrap intervals rest on four clusters.

---

## M1b predictions — policy rows: the radius stop rule and the leave-one-out screen

**Written 13 September, after M1 was run, read and scored (Appendix A), and before any
M1b quantity was computed.** What the author of this section had already seen and
could have anchored on: every M1 table, including Table 2, where `r_val = 1` on most
real-coefficient queries under the recovering set, and the M1 secondary table, where
Panel B has 7, 17 and 14 informative rows with a committed set at `b = 0, 1, 2`. The
refusal rates, the screen's triggers, the ceiling and every policy cell were not
computed from the M1 rows before this section was written. Only the ledger's
censoring statuses for `A_TRUE`, `D_LLM` and `D_LLM_72B` were checked: 16 rows are
censored at depth 3 and none at depth 1 or 2, so the certificate at `r <= 2` is exact
on every Table 1 row. Nothing in this section is edited after the M1b rows are
written; outcomes go in Appendix B. M1's rows, replicates and cells must not move: the
M1 columns of the new files are checked byte-identical against the M1 hashes of
Appendix A.0.

### B.1 The procedures, fixed here

- **STOP1, STOP2.** Certified at budget `r` iff `G0` commits to a set and `r_val > r`,
  read from the same `r_val` fields as M1 (the ledger's claim radius for `CTRL_b0` and
  the elicited conditions, `breaking_depths` for `b = 1, 2`): an exact value is compared
  directly; `unreached`, `no_retractable_edges`, `gt_d` with `d >= r` and
  `censored_at_depth_d` with `d > r` certify; anything else refuses and is counted.
  Certified: report PLAIN. Otherwise report the band `I_r`.
- **SCREEN.** For each single claim of the network-level claim set (the set `r_val`
  retracts), retract it, re-close under Meek and recompute the committed result with
  `committed_adjustment_set`. Triggered iff some result differs from `G0`'s. Triggered:
  report `I_1`; otherwise report PLAIN.
- **SCREEN2.** The same over every retraction of one or two claims; triggered: report
  `I_2`. It is defined when the retractions fit the 300-subset budget
  (`k + C(k, 2) <= 300`, which holds for every claim set of the Table 1 conditions,
  at most 21 claims); a row over budget would be counted as censored and reported
  as triggered.
- **The committed result compared by the screens** is the set when there is one, and
  otherwise the verdict class (`y_not_possible_descendant`, or no set). Closed-form and
  canonical sets are compared as sets.
- **Rows where `G0` commits to no set.** The stop rule has no `Z` to certify and, read
  literally, widens to `I_r`; the screen compares verdicts and reports PLAIN (M1's
  budget-zero fallback) unless a retraction changes the verdict or produces a set. This
  is the literal reading of both rules and it makes them differ on such rows by
  definition. **Choice written now:** the policy comparisons (Q1, Q2, Q4, the ceiling)
  are scored on informative rows **with a committed set**, the population on which the
  design run defined both rules; the literal all-row numbers are printed beside them and
  not scored. Table 1 cells (Q3, Q5) use M1's population, all informative rows, with the
  literal reading.
- **Closures for the screens** are assembled from per-component closures (`G0` with
  each chain component replaced by the Meek closure of its own claims). This equals the
  full closure by the component reduction; a test asserting equality with
  `apply_orientations` on the full graph is added before the run.

### B.2 What is computed

- The four policies as Table 1 rows at the population and at `n = 1000, 20000`, same
  replicates, same metric, same bootstrap (M1's cells are computed first with M1's
  generator so their bootstrap draws are unchanged; M1b cells use a separate generator
  seeded 20260913).
- **The starred comparison**: at `b = 2`, per panel and `n`, STOP2 minus SCREEN in raw
  coverage and in L95 over the same replicates, with a network cluster bootstrap
  (`B = 2000`, one resample drives both policies), and the per-network differences.
- **The identity check**: at `b = 1`, the number of rows where STOP1 and SCREEN take
  different decisions, with each difference classified (the set changes while `Z`
  stays valid; the screen stays silent while `Z` breaks, which the argument says cannot
  happen; no committed set).
- **The price of the guarantee**: at `b = 0`, the refusal rate of STOP1, STOP2, SCREEN
  and SCREEN2, the share of informative rows with a committed set on which the policy
  widens, per panel, with a network bootstrap interval; the all-row literal rate beside.
- **The ceiling**: the share of informative rows with a committed set in the
  interaction stratum, `r_val = 2` or some pair retraction changes the committed result
  while no single one does, per panel and condition.
- **An ill-posed point, written before the run.** STOP2 and SCREEN do not differ only
  where two retractions interact. Where a single retraction breaks `Z` both refuse, but
  STOP2 reports `I_2` and SCREEN `I_1`; where the set changes at depth one while `Z`
  stays valid to depth two, STOP2 reports PLAIN and SCREEN `I_1`. A coverage gap at
  `b = 2` can therefore come from the fallback's budget rather than the certificate's
  depth. The depth-matched comparison is STOP2 against SCREEN2, which by the depth-one
  argument can differ only where the set changes while `Z` stays valid. The gap is
  therefore decomposed by stratum (interaction stratum; `r_val = 1`; the rest), the
  literal share of rows where the two policies' reports differ is printed beside the
  interaction ceiling, and STOP2 against SCREEN2 is reported. Q2 is scored as written.

### B.3 Predictions

**Q1. The depth-one identity.** At `b = 1`, on informative rows with a committed set
pooled over both panels, STOP1 and SCREEN take different decisions on at most 1 % of
rows. Falsified above 1 %. Differences are listed with their class either way.

**Q2. The starred comparison separates, inside a small ceiling.** At `b = 2`, Panel A,
informative rows with a committed set, STOP2's raw coverage exceeds SCREEN's with a
network-bootstrap interval of the difference above zero at `n = 1000` and at
`n = 20000`; and the interaction ceiling at `b = 2` is at most 0.31 in each panel.
Supported if all hold, partial if some, falsified if none.

**Q3. The stop rule keeps the narrow report where errors cannot matter.** At `b = 0`,
Table 1 population, `L95(STOP1) < L95(I1)` in both panels at both `n`. Scored per cell,
supported if 4 of 4.

**Q4. The guarantee is affordable.** At `b = 0`, STOP1's refusal rate on informative
rows with a committed set is below 0.5 in Panel A and in Panel B. Falsified in a panel
at 0.5 or above. Panel B rests on 7 rows and is flagged as such.

**Q5. The stop rule covers after one wrong claim.** At `b = 1`, Table 1 population,
STOP1's raw coverage is at least 0.95 in both panels at both `n`. Scored per cell,
supported if 4 of 4.

## Appendix B — M1b outcomes, 13 September, after the run was read

### B.0 What was run and what did not move

- One additional test before the run: per-component closures equal the full closure on
  300 random claim sets (`test_closure_by_components_equals_full_closure`). Then a
  crash-only run on `asia`, `Polzer_2012` and `diabetes` into a scratch directory, read
  for decision counts. It showed rows where the stop rule certifies and the screen
  triggers, each with the set changing while `Z` stays valid, as B.1 expects.
- The full panel, 7 processes, 50 s. **M1 did not move**: projected onto M1's
  columns, the new `rows.csv` hashes to `ef7a98d0…` and the uncompressed replicates to
  `c0a9e783…`, the M1 hashes of A.0; `profiles.jsonl`, `table2.json`,
  `fig1_candidates.json` and `TABLE2_REPORTS.md` are byte-identical; all 3,287 M1
  values of `table1.json` and every M1 section of `summary.json` are unchanged.
- No row needed the pair budget (SCREEN2 exhaustive on every row); the stop rule's
  certificate was exact, `unreached`, `gt_3` or `censored_at_depth_3` everywhere, so
  no refusal came from censoring. On all 2,212 policy rows, a stop rule that refuses at
  budget `r` is matched by a screen that triggers at depth at most `r`: zero exceptions.

### B.1 The predictions

**Q1 FALSIFIED, narrowly.** At `b = 1`, informative rows with a committed set, STOP1
and SCREEN decide differently on 2 of 162 rows (1.23 %, one row over the 1 % bar). Both
are the expected class, a single retraction changing the committed set while `Z` stays
valid: `Kampen_2014`, ALN on PER (`Z` = {AFF, SAN}, `r_val = 3`), and `Schipf_2010`,
A on S (`Z` empty, `r_val = 2`), where retracting one claim makes the closed form return
a non-empty optimal set that the empty set is still valid for. The same class accounts
for 2 of 73 rows at `b = 2`, 0 of 141 at `b = 0`, 2 of 161 and 1 of 217 under the two
elicited conditions. On all informative rows, under the literal reading, 112 of 403
differ at `b = 1`, and 110 of those are rows with no committed set.

**Q2 SUPPORTED as written, and the separation is not the one the design run
wanted.** At `b = 2`, Panel A, rows with a committed set (59 rows, 8 networks): STOP2
minus SCREEN coverage +0.256 [+0.149, +0.346] at `n = 1000` (0.986 against 0.731) and
+0.269 [+0.163, +0.361] at `n = 20000` (0.971 against 0.702); L95 0.740 against 10.27
and 0.733 against 43.66. The interaction ceiling at `b = 2` is 0.237 [0.10, 0.33] in
Panel A and 0.286 [0.00, 0.53] in Panel B (14 rows). Two facts qualify it:

- About a third of the gap is the fallback's budget, not the certificate's depth. On the
  34 rows with `r_val = 1` both rules refuse, STOP2 reports `I_2` and SCREEN `I_1`; that
  stratum contributes +0.092 of +0.256 at `n = 1000` (+0.101 of +0.269 at `n = 20000`).
  The interaction stratum (14 rows) contributes +0.164 (+0.169), the rest +0.001.
- **A free screen at depth two matches the stop rule exactly.** STOP2 minus SCREEN2:
  coverage +0.000 [+0.000, +0.000] at `n = 1000` and -0.002 [-0.003, +0.000] at
  `n = 20000`, L95 -0.002 and -0.000, and zero in Panel B. The two decide differently on
  5 of the 754 committed informative rows over the five conditions, every one a set change
  with `Z` valid. The comparison fixed in advance separates only because its free
  baseline stops at depth one; the depth-matched baseline costs the same kind of graph
  operation and ties.

Panel B at `b = 2` (14 rows, 4 networks): +0.075 [+0.000, +0.175] and +0.121
[+0.000, +0.283], all of it from the 4 interaction rows; its intervals touch zero.

**Q3 SUPPORTED on the letter, and empty.** At `b = 0`, all informative rows,
`L95(STOP1) < L95(I1)` in all four cells, by 0.0011, 0.0004, 0.0007 and 0.0003. Under
the literal reading STOP1 widens on 97 % (Panel A) and 99 % (Panel B) of those rows, so
it is `I_1` almost everywhere; the prediction's reason, keeping the narrow report where
errors cannot matter, is not what the margin shows.

**Q4 FALSIFIED.** At `b = 0`, informative rows with a committed set, STOP1 refuses on
0.910 [0.76, 1.00] of 134 rows in Panel A (20 networks) and 0.857 [0.80, 1.00] of 7 rows
in Panel B (3 networks). STOP2 0.955 and 1.000; SCREEN is identical to STOP1 (0.910,
0.857); SCREEN2 0.963 and 1.000. Context computed after the fact: on committed rows
outside the informative stratum the stop rule refuses on none (57 rows in A, 11 in B),
so over every committed row the rate is 0.639 (A) and 0.333 (B). The refusals sit
exactly where the knowledge carries identification, and there they are nearly
universal: a single true claim, retracted, breaks the committed set on nine informative
rows in ten. The literal all-row rates are STOP1 0.969 and 0.986, SCREEN 0.517 and
0.406.

**Q5 SUPPORTED.** At `b = 1`, all informative rows, STOP1 covers 0.979 and 0.979 in
Panel A and 0.980 and 0.966 in Panel B, with L95 0.715 and 0.702 (Panel A) against
I1's 0.851 and 0.758.

### B.2 Ill-posed points, as found

- **The rules are undefined without a committed set.** Read literally, the stop rule
  widens and the screen compares verdicts, which makes them differ on those rows by
  construction (110 of the 112 all-row differences at `b = 1`). The policy comparisons
  are therefore on committed rows, where Panel B has 7, 17 and 14 informative rows.
- **STOP2 against SCREEN compares depth and fallback budget at once.** The depth-only
  comparison, STOP2 against SCREEN2, is a second identity and ties.
- **The refusal rate's denominator decides its size**: 0.91 on informative committed
  rows, 0.64 on all committed rows, 0 on committed rows the class already identifies.
