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
