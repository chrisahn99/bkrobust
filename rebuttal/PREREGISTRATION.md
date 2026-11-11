# Rebuttal experiments: plan and expected outcomes

Branch `iclr_2027_rebuttal`, cut from `iclr_2027` at `5421cab` (the submitted
code). Experiment ids follow `rebuttal/experiments.pdf` in the Overleaf project
(José's plan). Every expectation below was written **before** the corresponding
run; deviations are reported, not edited away.

Machine: Apple M5 laptop, 10 cores, 32 GB (not betelgeuse). Wall-clock times
are therefore laptop times and are labelled as such.

Data available locally: the submitted repository (results/, frame, elicitation
bundle, networks) and the release folder of the paper
(`release/per_query_records.csv`). Not available: José's `results/pest`,
`results/ledger`, `bkrobust-ze`, the closed-form optimal set. Experiments that
needed them are re-implemented from the submitted code and say so.

No new package was installed. PC discovery (E5-E7) is implemented in-house
(PC-stable skeleton with Fisher-z tests, v-structures, Meek closure from
`bkrobust.demo.meek`). AID (E15) is implemented in-house from its definition.

## Priorities

| tier | experiments | points |
|---|---|---|
| A, cheap, must | E21, E13, E16, E17, W4 tables, E1 | C2Bd-W1/Q2, qQYf-M2/Q1, C2Bd-W3/W4, qQYf-S2, h8pC-W1 |
| B, main compute | E2, E3, E10, E15, E18, E12, E11 | h8pC-W1, qQYf-W1, qQYf-Q2, h8pC-W3, C2Bd-Q1, h8pC-W4a/b |
| C, discovery error | E5, E6, E7, E8 (insertion arm) | C2Bd-W2/Q3, h8pC-W2 |
| D, if time | E19, E4, E9, E14, E20 | optional |

## Expected outcomes (written before running)

- **E21 toys.** Walkthrough toy (U-X, X-Y, X-W, W-Y; K = U->X, X->Y, X->W):
  expect a small space (< 20 states), r_val in {1, 2}. Figure 2(a) toy
  (A-W, W-X, W-Y, X-Y; K = A->W, X->Y): Figure 2(b) prints r_val = 3 next to a
  schematic; expect the true value to be smaller (1 or 2), since |K_G0| = 4.
- **E13.** (a) Units with no optimal set: expect most non-amenable or with
  extensions disagreeing; some amenable units (valid sets exist) are possible.
  (c) Corollary 1(c) on rows invalid at truth: expect 0 violations (any
  violation is a bug).
- **E16.** Radius on units where separation is undefined: expect positive
  tau_b, interval possibly including zero (fewer networks).
- **E17.** Expect at least one applied network with r = 1 and r >= 3 queries in
  one state; the share invalid at truth by radius bin likely similar for LLM
  and random states (the gap describes CPDAGs more than knowledge).
- **W4.** Corrected-cohort Table 13 per condition: point estimates move by
  < 0.05 for most conditions.
- **E1.** Of the 372 lower-bound rows, expect a majority settled by the top
  check as infinite (Z valid in C_hat) and the rest finite at depth > 3;
  0 rows where a found failure coexists with an infinite answer.
- **E2.** GAC hybrid: median time grows polynomially with n while the top
  check settles most queries; expect tail cases (SAT) of seconds to a minute
  for n <= 100; censoring possible beyond.
- **E3.** Same radius in every (b, mode) arm. Expect b in {2, 3} near the
  optimum; the single minimise call faster than the ladder on deep radii.
- **E10.** Expect the radius to vary across valid sets on a minority of
  queries, and O to be beaten on < 5% of queries with >= 2 valid sets.
- **E15.** Expect oracle AID (G0 vs truth) to rank survival weakly (like SHD);
  analyst-side AID (G0 vs C_hat) to rank positively but below r_val on tiered
  strata.
- **E18.** Post hoc standardisation: profiles identical to 1e-9 (else bug).
  Internal standardisation: ceilings and signatures change; radius unchanged.
- **E12.** Expect the top check to settle many rows on the five > 500-node
  networks within seconds; finite radii computable for >= 2 networks.
- **E11.** Expect a positive within-state tau_b on the CPDAG-only frame,
  smaller than 0.43; interval may include zero.
- **E5.** Expect invalidity at truth on learned CPDAGs to fall with n but stay
  well above zero at n = 5000 (systematic orientation/skeleton error).
- **E6.** Expect bootstrap stability of Z to separate invalid from valid sets
  modestly (AUC 0.6-0.75).
- **E7.** Expect a hidden variable adjacent to the query to raise invalidity
  more than a distant one.
- **E8 (insertion arm).** Expect at least one case where adding a true-oriented
  extra edge breaks a set valid in the true DAG? No: an extra edge is absent
  from the truth, so validity at truth is unchanged by definition; the
  meaningful test is whether Z committed at the perturbed G0 stays valid at
  truth. Expect damage concentrated within hop 1 of the query.

---

# E4L (10 November 2026): listing the class, at scale

Written before any run of this experiment. José's E4 (branch `ze/rebuttal`,
`experiments/rebuttal/E4_class_enumeration_baseline`) compared listing the class
with the hybrid on three small populations (four-node; 144 regenerated dense;
256 timing instances, n <= 24; caps 10^5 DAGs, 120 s) and found listing faster
on 8% of the timing instances and 27% of the four-node ones. E4L does not repeat
it; it extends the same comparison to the populations the paper and the
response actually use, and measures how long listing would take where it cannot
finish.

**Populations** (all rebuilt from stored parameters with the submitted code;
Z = closed-form O(G0), checked valid at G0):
- P-E2: the 2,676 queries of E2 (n = 8..400, dense App. D grid and sparse
  degree 2-4, full and half recovering knowledge);
- P-REL: the 2,365 committed release queries of the ten elicited conditions;
- P-BIG: the committed queries on the five networks above 500 nodes (E12);
- P-SYN: the 4,584 synthetic survival instances of Section 6.

**Legs, per instance, in one process:** class size |[C^]| counted exactly by
root-picking (He, Jia & Yu 2015; cap 60 s); hybrid (final, GAC, b = 3);
LIST = depth-first listing of [C^] with Meek closure and DAG-level GAC in every
DAG (Corollary 1(c)); BNB = the same listing pruned on reversed K_G0
orientations. LIST runs to completion when |[C^]| <= 10^6 (cap 300 s);
otherwise for 20 s only, to measure its throughput. BNB cap 120 s. The listing
and pruning rules follow José's E4 so the two runs are comparable.

**Registered statistics.** (1) Per population and class-size bin
(1, 2-10, 11-100, 101-10^3, 10^3-10^4, 10^4-10^6, 10^6-10^9, >10^9, unknown):
fraction of instances where LIST finishes faster than the hybrid (Wilson 95%).
(2) Projected LIST time = |[C^]| / measured throughput, by population, against
the hybrid time; share of instances where the projection exceeds 1 hour and
1 year. (3) Same as (1) for BNB. (4) Validation of the projection on a
long-cap subsample: 60 instances with |[C^]| in 10^5..10^8, stratified by
decade, LIST given 1,800 s; ratio of actual DAGs listed to projected.

**Prediction.** Listing is faster only on classes of size <= 100; on P-E2 at
n >= 100 and on P-BIG its projected time exceeds 1 hour on most instances with a
finite radius, while the hybrid takes under a second; BNB beats the hybrid on
few instances outside small classes. Throughput is roughly constant per
instance, so the projection is accurate within a factor of 2.

**Consistency (any failure is a bug).** A complete LIST or BNB radius equals the
hybrid radius; a complete LIST count equals |[C^]|; a capped LIST that found a
failing DAG gives an upper bound >= hybrid radius.

**Against us.** If BNB matches or beats the hybrid on most large-class
instances, the response says a pruned listing is a competitive exact method and
the hybrid's advantage is engineering, not asymptotics.

**E4L amendment after the smoke run (24 units), before the full run.** A class-size
census (counting only, 15 s cap) gives about 280 units with |[C^]| in 10^4..10^6 and
about 830 above 10^6; at the caps above the main run needs about 8 hours on 8
workers. The long-cap validation subsample is therefore cut from 60 to 30
instances (10 per decade 10^5, 10^6, 10^7), still with 1,800 s each. Nothing else
changes. Smoke: no consistency failure; listing and pruned listing agree with the
hybrid on every complete unit.

# E3S (10 November 2026): the SAT stage against other exact solvers

Written before any run of this experiment (branch `iclr_2027_rebuttal_solvers`,
cut from `iclr_2027_rebuttal` at the E4L results). Addresses qQYf-W1 ("the
choice of SAT over other optimization approaches ... is not thoroughly
motivated").

**Question.** The SAT stage of Algorithm 1 answers "is there a Meek-closed state
within k retractions of G0 where Z fails (GAC)?" with CP-SAT, as a ladder over k.
Do other exact formulations of the same problem do better?

**One model, several solvers.** The submitted E1 encoding (`bkrobust.sat.e1._build`,
criterion "gac") is built once per query. Its constraint types are bool_or,
bool_and (enforced), exactly_one, lin_max over Booleans, linear (cardinalities
and enforced position differences for acyclicity). It is translated mechanically,
with no reformulation, to
- a MILP (big-M on enforced linear constraints with M from the variable bounds;
  lin_max of Booleans as out >= each, out <= sum), minimising the number of
  retractions;
- a weighted CNF (positions in order encoding; cardinalities by a sequential
  counter), each retraction a unit soft clause of weight 1 (MaxSAT).
Translators are checked on every instance: every returned witness is decoded to
a state, re-closed, and checked to make Z fail with the claimed number of
retractions; every exact answer must equal the hybrid radius.

**Arms** (no search stage, start at k = 0, single thread, cap 300 s per arm per query):
A1 CP-SAT ladder (submitted); A2 CP-SAT one minimise call; A3 CP-SAT minimise with
core-based search only (`optimize_with_core`, a MaxSAT-style algorithm inside
CP-SAT); A4 MILP SCIP; A5 MILP HiGHS; A6 MILP CBC (all three bundled with OR-Tools
9.15); A7 MaxSAT RC2 and A8 plain SAT ladder with CaDiCaL (PySAT), only if PySAT
is installed (decision recorded below before the run); A9 explicit retraction
search with unbounded depth (`radius_local_up_fast`, max_depth = |K_G0|), the
"explicit search" of the question.

**Populations.**
- P-E3: the 132 E3 queries (dense ER, n = 10..24, radii 1..13; 63 with r >= 4).
- P-SATR: the E4L queries the hybrid sent to its SAT stage (1,911; 1,847
  synthetic, 56 release, 8 E2); all if the budget allows, otherwise a seeded
  sample of 400 stratified by population and radius.
- P-SCALE: the deepest finite query per (n, family) for n in {24, 32, 48, 64}
  from the E2 sparse and dense generators, 10 per n, to see where each
  formulation stops building or solving.

**Registered statistics.** Per arm and radius class (1-3, 4-6, >= 7): share
solved exactly within the cap (Wilson 95%); median time (build + solve) with
order-statistic interval; paired time ratio to A1 (median, cluster bootstrap
over graphs); share of queries on which the arm is fastest; model size
(variables, constraints, clauses). The ladder's certificate time (UNSAT rungs)
reported separately for A1 and A8.

**Prediction.** CP-SAT (A1 or A2) is fastest on most queries in every radius
class. MILP solves fewer deep queries within the cap and is slower by a median
factor of at least 10 on r >= 4 (weak LP relaxation of logical and big-M
constraints). MaxSAT (A7, A3) is within a factor of 3 of CP-SAT and may win on
some deep queries. Explicit search (A9) wins on r <= 2 and loses from r >= 4.

**Against us.** If any non-CP-SAT arm is fastest on the majority of r >= 4
queries, the response says so and that CP-SAT was an engineering choice, not an
optimal one; App. D.3 reports it.

**Consistency (any failure is a bug, to be fixed before analysis).** Every exact
arm returns the hybrid radius; every witness re-checks.

**E3S decision before the run.** PySAT is installed (python-sat 1.9.dev21, user approval
10 November 2026), so arms A7 (RC2) and A8 (CaDiCaL ladder) are included.

**E3S amendment after the smoke run (12 queries) and a 16-query census of
P-SATR, before the full run (10 November 2026, 19:40).** Smoke: no radius
mismatch, every witness re-checks; the three MILP arms hit the cap on most
r >= 4 queries. P-SATR is cheap (median under 0.3 s per arm except CBC), so
all 1,911 queries are run and no sample is drawn. The P-SCALE query id now
includes the edge probability (two smoke queries shared an id). Ladder rung
times are recorded for A1 and A8. Nothing else changes; cap 300 s, 8 workers.

**E3S amendment 2, after the full run (10 November 2026, 22:05).** The full run
(18,747 tasks, 19:35-22:01) produced 378 rows with status
`error:AssertionError:ok`: all 9 arms on 42 P-SATR queries. Cause: a driver bug,
not a solver failure. `build()` in `rebuttal/e3s_solvers.py` asserted a
non-empty optimal set, but on these 42 queries O(G0) is the empty set (a valid
adjustment set, and included in E4L). Fix: the assertion now checks only the
closed-form status. The 378 error rows were moved to
`results/rebuttal/e3s/errors_emptyZ.jsonl` and the 378 tasks rerun with the same
settings (cap 300 s, 8 workers; log `run_amend2.log`). No other row was touched.
After the rerun: 0 errors, 0 radius mismatches, 0 witness failures. The
descriptive statistics in `rebuttal/e3s_explore.py` (fastest arm among the
solver arms only, solve-only ratios, explicit search by radius) were not
registered and are reported as exploratory.

---

# E30-E34 (11 November 2026): finite-sample discovery and the radius

Written before the main run of any of these experiments (branch
`iclr_2027_rebuttal_discovery`). Addresses C2Bd-W2 ("finite-sample discovery errors
... are not captured"), C2Bd-Q3 and h8pC-W2 ("the metric neglects the bias caused by
the recovery of the CPDAG"). Laptop (Apple M5, 10 cores, 32 GB): all times are laptop
times. Package installed with the user's approval (11 November 2026): causal-learn
0.1.4.8 (pip, into `.venv`). No Java, no pytetrad.

**Framing (the user's, kept).** The radius is a sensitivity analysis of the orientation
knowledge *conditional on the equivalence class*. Judging "discovery + radius" as one
pipeline mixes two error sources. So the registered statistics separate "class locally
correct near the query" from "class wrong near the query", and ask how fast the former
becomes the norm as n grows, per learner.

## What was seen before writing this (disclosed, not hidden)

- Population edge strength (no data): partial correlation of each true edge given the
  child's other parents; 50%-power sample size of a Fisher-z test at alpha = 0.01 is
  <= 5,000 for 99% of edges in every network, <= 7,651 except one arth150 edge
  (about 750,000).
- **Oracle check of E5's PC at its depth cap.** With oracle d-separation tests and
  conditioning sets of size <= 3 (E5's setting), PC returns the true CPDAG on ecoli70
  only; it keeps 73 / 48 / 40 extra adjacencies on magic-irri / magic-niab / arth150
  (true edges 102 / 66 / 150). With no depth cap it returns the true CPDAG
  (`pc_oracle_check.py` used the maximum degree, which is why E5's validation passed).
  So E5's learner is not consistent on three of the four networks; part of the E5
  plateau may be this cap, not the sample size. E5's learner is kept as arm `pc01d3`;
  the other PC arms have no cap.
- causal-learn timings on seed 0, n = 5,000 (edge counts only): BOSS 0.6 / 4.6 s,
  GRaSP 0.7 / 4.4 s, GES 1.5 / 118 s on magic-niab / arth150.
- Structure errors (missing / extra / mark differences against the true CPDAG) of
  pc01d3, pc01, pcan, cpc, hc on magic-irri and arth150, seed 0, n = 100 and 100,000.
- Five smoke tasks of the driver (committed / invalid counts, all arms pooled):
  oracle magic-niab 82 / 3 (LLM arm); pc01d3 magic-niab n = 1,000 seed 0: 1 / 0;
  boss ecoli70 5,000 s0: 60 / 7; hc ecoli70 5,000 s0: 38 / 9; cpc ecoli70 5,000 s0: 7 / 0.
  Theorem violations: 0. These tasks are part of the registered population and are rerun.

## Common design

- **Networks.** ecoli70, magic-irri, magic-niab, arth150 (fitted linear-Gaussian
  bnlearn networks, `discovery.GAUSSIAN`), sampled ancestrally with `discovery.sample`.
- **Ladder.** n in {100, 200, 500, 1k, 2k, 5k, 10k, 20k, 50k, 100k}; seeds 0-9. Seeds
  0-4 at n in {200, 1k, 5k} are E5's data sets (same sampler).
- **Learners** (`rebuttal/learners.py`; all on the same samples):
  `pc01d3` E5's PC-stable (Fisher-z, alpha 0.01, depth <= 3, sepset collider rule);
  `pc01`, `pc001` the same without depth cap at alpha 0.01, 0.001;
  `pcan` without cap, alpha_n = 1/sqrt(n);
  `cpc`, `mpc` as pcan with the conservative / majority collider rule (ambiguous
  triples left unoriented; Meek may propagate through them, unlike Tetrad's CPC);
  `hc` greedy hill climbing over DAGs with Gaussian BIC (lambda 0.5), CPDAG of the result;
  `ges`, `boss`, `grasp` causal-learn, `local_score_BIC_from_cov`, lambda 0.5 (standard
  BIC; BOSS/GRaSP's package default is 2), GRaSP depth 3.
- **Frame and arms.** E5's 40 pairs per network (`e5_discovery.frame`: y a true
  descendant of x, x on an undirected edge of the true CPDAG) and E5's knowledge arms
  (`e5_discovery.arms`): TRUE (truthful orientation of every undirected edge of C_hat
  that is a true edge, added in sorted order, skipping those that do not close; José's
  A_ORACLE_R), LLM (released D_LLM claims), NONE. **TRUE is the registered arm.**
  Reference: the same pipeline on the true CPDAG ("oracle").
- **Per committed set** (status ok: amenable, O(G0) exists and is valid at G0):
  validity in the true DAG (GAC); radius on C_hat (hybrid, 30 s per rung, exactness
  recorded); local check on A = An_true({x, y} u Z) (inclusive):
  `class_ok` = no true edge inside A missing from C_hat, no directed C_hat edge inside A
  reversing a true edge, and the true orientations of the undirected C_hat edges inside
  A Meek-close; `g0_ok` = the same on G0; `k_err` = |K_G0 \ dir(D)| for one D in
  [C_hat] built from the class_ok orientations (closure, then G0's closure orientations
  that still close, sorted, then Dor-Tarsi).
  Theorem used (Remark 8(i) localised to the ancestral set A, with Corollary 1(c)):
  `class_ok` and (r = inf or k_err < r) imply Z valid in the true DAG.
  Exploratory: agreement of C_hat with the true CPDAG on all pairs within 1 or 2 hops
  of {x, y} in the true skeleton.
- **Per graph:** skeleton precision / recall, v-structure precision / recall, SHD
  between CPDAGs, number of undirected edges, whether C_hat is a CPDAG (equals the
  CPDAG of its Dor-Tarsi extension), learning time.
- **Intervals.** Percentile bootstrap, 2,000 draws, seed 20261111, resampling seeds
  with replacement within each (network, n) cell (graphs are the clusters); pooled
  statistics weight networks by their committed rows, as in E5.
- **Timeouts and errors.** A task that raises is recorded with its traceback and
  counted; no task is dropped silently. Radius rungs that time out are recorded as
  non-exact.

## E30: the sample-size ladder and the coherence threshold

**Registered statistics** (TRUE arm, per learner, pooled over networks; also per network):
I(n) = share of committed sets invalid in the true DAG (oracle: 0 by construction);
Cm(n) = share of the 160 frame queries with a committed set (oracle: 1.00 in E5);
V(n) = share of frame queries with a committed set that is valid (fixed denominator).
**Coherence threshold** n*(learner) = the smallest ladder n such that at n and at every
larger ladder n, I(n) <= 0.05 and Cm(n) >= 0.95 (point estimates); "not reached" if none.
Also reported at (0.10, 0.90).

**Predictions.** (a) For ges, boss, grasp and pcan, I(n) falls with n: I(100k) < I(1k)
with the bootstrap interval of the difference below 0. (b) At least one of ges, boss,
grasp reaches n*(0.05, 0.95) by 100k. (c) pc01d3 does not reach it (depth cap), and at
n >= 20k its I(n) exceeds that of pc01. (d) seeds 0-4 of pc01d3 at n in {200, 1k, 5k}
reproduce E5's TRUE-arm rows (status, Z, validity) exactly.

**Against us.** No learner reaches n*(0.05, 0.95) by n = 100k: then we report that on
these networks the learned class is not coherent near the query even at 100k samples,
for any learner tried, and the answer says so with the measured I(100k).

## E31: learners compared

The same statistics side by side on the same samples, plus graph-level accuracy and
learning time. **Prediction:** the score-based learners (ges, boss, grasp) have lower
I(n) than PC at alpha 0.01 at n >= 5k. Reported, not registered: cpc and mpc against
pcan (do careful collider rules lower invalidity at the cost of fewer committed sets?).

## E32: the conditional guarantee (main)

**Registered statistics** (TRUE arm, per learner and n): P_LC(n) = share of committed
sets with `class_ok`; the decomposition
I(n) = (1 - P_LC) * P(invalid | not class_ok) + P_LC * P(invalid | class_ok);
the number of **theorem violations** (class_ok, r = inf or k_err < r, Z invalid in the
truth), in every arm (TRUE, LLM, NONE), learner and n.

**Predictions.** (a) 0 theorem violations (any violation is a bug or a false theorem
and stops the analysis until found). (b) P(invalid | class_ok) is near 0 in the TRUE
arm (it is 0 whenever k_err < r or r = inf). (c) For ges, boss, grasp and pcan,
P_LC rises with n and reaches >= 0.9 at n = 100k. (d) P(invalid | not class_ok) stays
>= 0.3 at every n with at least 30 such rows, i.e. invalidity is carried by local
class error.

**Against us.** P_LC stays below 0.5 at n = 100k for every learner: local class errors
do not become rare on these networks, and the answer says so.

## E33: does the radius on the learned class converge to the oracle radius?

**Registered statistics** (TRUE arm; queries committed both on C_hat and on the true
CPDAG): P(r_hat = r_oracle), Kendall tau_b(r_hat, r_oracle) (inf ranked highest),
and P(Z_hat = Z_oracle), per learner and n. Also for the LLM arm (secondary).

**Prediction.** For ges, boss, grasp and pcan, P(r_hat = r_oracle) rises with n and is
>= 0.8 at 100k. **Against us:** below 0.5 at 100k for every learner.

## E34: does the radius still rank survival on learned classes?

**Design** (`rebuttal/e34_survival.py`, to be written; reads the learned CPDAGs stored by
E30, no rerun of discovery). Learners pcan, ges, boss; n in {500, 5k, 50k}; seeds 0-4;
plus the true CPDAG. Units: committed sets of the TRUE arm (registered) and the LLM
arm (secondary). Claims K = the orientations the arm actually applied. The paper's
flip arm: at each depth of `real_survival.depth_grid(|K|)`, 200 draws (the paper uses
1,000; reduced for time) of `flip_draw` (seeded by the unit id), Meek-close on C_hat,
survival = share of non-contradictory draws in which Z stays valid;
survival AUC = `auc_over_grid` over the depth grid.

**Registered statistic.** Kendall tau_b between the radius on C_hat (inf highest) and
the survival AUC, pooled over units per (learner, n), with the graph-cluster bootstrap;
the same on the true CPDAG as reference.

**Prediction.** tau_b > 0 with the interval excluding 0 for every learner at n >= 5k,
and within 0.15 of the oracle-class value at n = 50k. **Against us:** an interval that
includes 0, or a negative tau_b, at n >= 5k.

## Exploratory (labelled so wherever reported)

Radius bins and the AUC of 1/r for invalidity at truth on learned classes, by n;
edge strength near the query against local error; cpc/mpc against pcan; agreement
within 1 and 2 hops; the LLM and NONE arms beyond what is registered above.

## Stopping rule

One launch of the full grid (E30-E33 share one run), resumed only after a crash.
No rerun with changed settings. A failed consistency check is fixed, rerun, and both
runs are reported as amendments below.

**E30-E34 amendment 1 (11 November 2026, during the E30 run, before any E30 row was
analysed): E35, local agreement between learners as a data-only flag (C2Bd-Q3).**
E30 stores every learned CPDAG, all learners on the same samples, so this costs no new
run. For each (network, n, seed) and each TRUE-arm committed set of the reference
learner `boss`, the region R = {x, y} u Z plus their neighbours in C_hat_boss (data
only, no truth). Flag = the CPDAGs of boss, ges, grasp and pcan on the same sample
disagree on some pair inside R (adjacency or mark). **Registered statistics:** the
invalid share among unflagged sets, per n, with the graph-cluster bootstrap; AUC of the
flag for invalidity at truth. **Prediction:** unflagged invalid share <= 0.05 at
n >= 10k; AUC > 0.6. **Against us:** the AUC interval includes 0.5, or the unflagged
invalid share is above 0.10 at n >= 10k. Reported with the share of sets flagged
(a flag that flags everything is useless).

**E30-E34 amendment 2 (11 November 2026, before any E34 run): E34 cost.** One Meek
closure of the corrupted claims takes 0.8 / 3.3 / 7.8 / 54 ms on the true CPDAGs of
magic-niab / magic-irri / ecoli70 / arth150 (laptop), so the registered E34 (200 draws
per depth, every committed set) would take more than 40 laptop-hours. Changed before
running: 100 draws per depth (the paper uses 1,000), and at most the first 10 committed
sets per (graph, arm) in frame order (fixed before seeing any outcome). Everything else
as registered.

**E30-E34 amendment 3 (11 November 2026, 05:57, after 100 of 4,004 E30 tasks, before
any analysis): learning time cap.** causal-learn's GES took 3,141 s on one magic-irri
graph (n = 100k, seed 3), and 7 of 8 workers had been inside GES fits on arth150 /
magic-irri for over an hour. The run was stopped and resumed (the 100 finished tasks are
kept, unchanged) with a cap of 900 laptop-seconds on the learning step (SIGALRM). A
capped task is recorded as `learn_timeout` and counted per (learner, network, n); it
has no rows. The statistics are computed on the graphs that finished, and the counts
of capped graphs are reported beside every pooled number for that learner. The one
3,141 s graph, which finished before the cap existed, is kept.

**E30-E34 amendment 4 (11 November 2026, 11:25): interruption.** The E30 run was
stopped on the user's request (laptop to be closed) after 3,667 of 4,004 tasks, all
written lines valid. Tasks in flight were discarded, not written. The run is resumed
with the identical command; finished tasks are skipped, so the only effect is that the
in-flight tasks are recomputed from scratch on the same samples. No setting changes.

**E30-E34 note 5 (11 November 2026, 12:00): interim analysis on 3,667 of 4,004 tasks,
at the user's request (deadline).** Missing: 337 arth150 tasks (boss, grasp, hc, pc01,
pc001, pc01d3, pcan at n >= 5k). Because the pooled four-network numbers then change
composition with n, the interim report pools the three complete networks (ecoli70,
magic-irri, magic-niab) and gives arth150 separately. The registered four-network
statistics are computed after the run completes and supersede the interim ones. E34 has
not been run.
