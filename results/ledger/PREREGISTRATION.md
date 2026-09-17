# Pre-registration — the certificate ledger

**Written 13 September, after the questionnaire was frozen and hashed, after a
two-network pilot of the primary supplier (asia, Schipf_2010: nine pairs, eight
asserted, one declined, no Meek inconsistency, one control edge wrong), and
before any other supplier answer or any ledger row had been looked at.** The
supplier run was in progress when this was written; its outputs were not read.

Nothing below is edited after `results/ledger/rows.csv` is written. Outcomes go
in Appendix A.

---

## 1. What the ledger is

One row per knowledge arm, on the frozen frame of `results/frame/frame.jsonl`
restricted to the networks the questionnaire covers. Deployment arms assert
nothing read off the truth: a language model on real names (`D_LLM`), the same
model on a second presentation order (`D_LLM_INSTR`), a second model family
(`D_LLM_SRC`), the same model on scrambled names (`D_SCRAMBLED`), a chance-level
control on the model's own edges (`D_RAND`), and no knowledge (`D_DEGEN`). Audit
arms are read off the truth and sit below the rule: truthful orientations on the
model's edges (`A_TRUE`) and the greedy recovering set (`A_ORACLE`).

The CPDAG is the true DAG's, so every deployment row is stamped
`via_Chat_only`. That is the idealisation panel. It is not a deployment claim.

## 2. Fixed before the sweep

- The committed set is the polynomial closed form with the canonical fallback
  (`mpdag_criterion.optimal`), validated against the enumerating definition on
  822 rows with 818 exact agreements and no invalid set. Rows where the
  analyst's graph makes the outcome not a possible descendant of the treatment
  are a status, `y_not_possible_descendant`, not a certificate.
- `r_claim` by exhaustive retraction to depth 3 under a 300-subset budget;
  `r_hop` from the hybrid with a 5 s ladder wall and its exactness carried;
  `phi_1` exhaustively; `NO_RETRACTABLE_EDGES` when the knowledge orients
  nothing the CPDAG left open.
- Commission is the number of asserted claims false at the truth; omission is
  the number of orientable edges in the treatment's component the arm did not
  assert. They are separate columns and are never summed.
- Every rate carries a network cluster bootstrap. The five networks above 500
  nodes are excluded from the sweep for cost and the exclusion is printed.

## 3. Predictions

**L1. The primary supplier is wrong, and measurably.** Its false-claim fraction
over asserted claims lies in [0.10, 0.45] pooled over networks. Falsified below
0.05, in which case a seven-billion-parameter model reproduces these benchmarks
and the scrambled arm becomes the headline; falsified above 0.50, in which case
it is no better than chance and the arm is renamed a mechanical source.

**L2. Order elicitation keeps most networks Meek-consistent.** The protocol's A3
check predicted a strictly positive rejection rate on a majority of networks
under elicitation. This registers the opposite for order-based elicitation:
fewer than half the networks need any claim dropped for consistency. A cycle
cannot come from an order, so inconsistency can only come from an orientation
that creates a new v-structure against the CPDAG. Falsified if half or more
networks need repair.

**L3. Names carry information.** The scrambled arm's false-claim fraction
exceeds the real-name arm's. If they are within each other's intervals, the
names bought nothing, and whether the model was reciting or reasoning cannot be
decided on this corpus; the paper says so.

**L4. Chance is worse than the model, and this is what makes the model a
knowledge source rather than a cardinality.** `D_RAND`'s false-claim fraction is
near one half, and its committed set is invalid at the truth more often than
`D_LLM`'s, with intervals that separate. Falsified if they overlap, in which
case the model carries no knowledge on this substrate and the paper leads with
the identifiability step and the audit alone.

**L5. Knowledge buys identifiability first.** `D_DEGEN` certifies fewer than
one in ten frame rows. Every knowledge arm certifies more. This is the step
the other branch found from the efficiency side, measured on the frame.

**L6. The certificate never over-certifies on any deployment arm.** The
retraction certificate's dangerous rate is zero on every arm; the free rule's
is at least as high on every arm. The certificate detects every truth-invalid
committed set. Falsified by one dangerous row for the certificate.

**L7. Certified moves order by truth.** On the same edges, the truthful arm
certifies at least as many safe moves as the model arm, which certifies at
least as many as the chance arm. If the chance arm certifies more than the
model, the certificate is rewarding well-placed wrong knowledge and the paper
must say what it is a certificate of.

**L8. Instrument variance is smaller than source variance.** The disagreement
between the model and itself on a second order is smaller than between the two
families. Falsified if reversed; the prompt would then be doing more work than
the knowledge.

**L9. The compelled-edge control agrees with the audit.** The primary
supplier's error rate on edges the data already orient, which needs no truth,
lies within the interval of its false-claim fraction on open edges, which does.
If it does, the truth-free bridge the protocol wants for the no-truth tier is
supported on this corpus; if not, it is reported as failed.

**L10. Leverage stays heavy-tailed under elicited knowledge.** Median 1,
maximum above 5.

## 4. Decision rules

- L6 falsified: halt and inspect the row before anything else is read; a
  dangerous certificate row is either a bug or the result.
- L4 falsified: the elicited arm is uninformative on this substrate; the
  ledger is reported with that verdict and the abstract leads with the step.
- L1 low side: the scrambled arm and the control block decide whether the
  benchmark was recited; the paper reports whichever they say.

## 5. What this ledger cannot say

The CPDAG is the oracle's on every row. The supplier is a seven-billion-parameter
model and an eight-billion one; neither is a domain expert and the paper does
not call them one. Two of the four named suppliers in the protocol, the source
publications and human analysts, are not here. Five networks are excluded for
size. The frame is capped at twenty pairs per network.

---

## Appendix A — outcomes

### A.0 One instrument setting changed before any row was read

The first full sweep was stopped after eleven networks. On the fourth frame
row of the 223-node network the hop instrument fell through to the ladder,
whose model is built over every ordered node triple; the build had not
returned after twelve minutes and would not have on the 413-node network
either. The wall in section 2 bounds the solver, not the build. The rows
written up to that point were discarded unread and the sweep was rerun with
the hop instrument set as follows: the upward walk is exhaustive when the
analyst's graph has at most twelve retractable edges, otherwise depth three;
the ladder runs only on graphs of at most sixty nodes; the walk spends at most
four thousand closures. A walk that stops on a budget reports the depth it did
check, so the row carries the anytime statement `r_hop >= depth + 1` and is
marked inexact; the audit bucket for that instrument uses that lower bound,
which is what the claim instrument already did with its own censoring. No
prediction was touched. Two smaller things fixed in the same rerun: the
chance arm's seed now derives from a stable hash of the network name rather
than the interpreter's salted one, so a rerun draws the same orientations; and
the sweep can be split across processes by network and joined back.

### A.0b The halt rule fired on the first read, and found two bugs and one result

The first complete table showed the claim certificate dangerous on seven rows,
six under the order replicate and one under the second family, and the hop
certificate dangerous on twenty-nine. L6's decision rule says halt and inspect
before anything else is read, and the inspection sorted the seven into two
bugs and the twenty-nine into a result.

Six rows on the seventeen-node textbook chain, where the replicate asserted
thirteen claims and every one was false: the exhaustive claim search was
censored by its subset budget at depth three, and the bucket treated a
censored radius as an unbounded one. The certificate had said `r_claim >= 3`;
the bucket scored it as infinite. Fixed by scoring every censored radius at
the floor it certifies, which is what the claim column already did for a
search that finished its depth. A search that tries every subset and finds no
failure is now stamped `unreached` rather than `gt_d`, so the two cases are
told apart.

One row on the gene-expression network: the certificate said no retraction of
up to three of four claims breaks the set, the truth said the set was
invalid, and three of the four claims were false, which cannot both hold. An
independent implementation of the generalised adjustment criterion on
networkx sided with the certificate. The repository's DAG oracle is Pearl's
back-door criterion, which forbids every descendant of the treatment, and the
committed set contained one that lies on no causal path to the outcome. Across
the 1,169 committed sets of the elicited and oracle arms the back-door oracle
rejected four such sets and agreed everywhere else. The audit now scores
validity with the criterion the instrument uses, written separately in
`benchmarks.audit`, and keeps the back-door verdict as a column.

The twenty-nine hop rows are the unit gap and are left as the result. On
twenty-two of them the claim radius is exactly one, the hop radius exactly
two, one claim is false and the set is invalid at the truth: retracting the
false claim from the CPDAG removes its propagations with it and breaks the
set in one claim, while the hop walk un-orients one edge at a time and needs
two. The hop instrument certifies a safe move that does not exist. That is
the reason the paper's certificate is in claims.

The sweep was rerun with the two fixes and nothing else changed; the numbers
in A.2 are from that rerun.

### A.1 Resolved by the knowledge-level audit, 13 September, before the sweep was read

`experiments/knowledge_audit.py` scores the frozen knowledge against the true
DAGs and reads nothing from the sweep. Rates are network means with a cluster
bootstrap over the 32 networks the questionnaire covers.

**L1 SUPPORTED.** The primary supplier's commission is **0.325 [0.210, 0.449]**
by network mean and 53 of 171 pooled. Inside the registered band. A
seven-billion-parameter model is wrong on about a third of what it asserts,
and it asserts 171 of 369 pairs, declining 137. By tier: 0.27 on the benchmark
networks, 0.28 on the applied papers, 0.44 on the gene-expression networks,
0.75 on the textbook diagrams with abstract variable names, where there is
nothing to know.

**L2 SUPPORTED, against the protocol's own expectation.** Six of 32 networks
needed a consistency repair under the primary supplier, seven under its
replicate, eight under the second family, four under scrambled names. A
minority every time. Order elicitation does what it was chosen for. The repair
is not small where it fires: the primary supplier made 214 claims, 43 of them
were dropped to restore consistency, and 171 remained. The consistency check
practitioners already run therefore removes one in five of what this source
says, and does nothing about the rest. The second family lost 41 of 264, the
order replicate 57 of 196, the scrambled run 12 of 163.

**L3 NOT SUPPORTED.** Scrambled names give commission **0.343 [0.260, 0.431]**
against 0.325 [0.210, 0.449] on real names; on the compelled-edge control,
0.289 against 0.191, intervals overlapping. The model is about as often wrong
when it answers under either naming; what changes is that it declines more
under scrambled names (omission 0.568 against 0.450). Names did not buy
accuracy. The reading is that this source is not reciting the benchmarks; a
reciting source would collapse under scrambled names and it does not. That is
the outcome the contamination probe exists to be able to return, and it is
reported as such: the recitation worry is not detected on this supplier, and
neither is much domain knowledge.

**L8 DIRECTION SUPPORTED, NOT RESOLVED.** Instrument disagreement 0.160
[0.053, 0.286] on 87 shared edges against source disagreement 0.220
[0.118, 0.333] on 119. The order right, the intervals overlapping. The
scrambled replicate disagrees with the real-name run on 0.252 of 74 shared
edges, which is the largest of the three and is the same statement as L3 from
the other side.

**L9 NOT SUPPORTED, in the direction the protocol pre-declared.** The
compelled-edge control, which needs no truth, gives an error rate of **0.191
[0.094, 0.293]**; the open-edge commission, which does, is 0.325
[0.210, 0.449]. The truth-free rate sits below the truth-requiring interval. The
protocol's fourth realism rule said the bridge would be biased this way, because
compelled edges sit in v-structures and cascades and are not exchangeable with
the edges an analyst must assert; the measured correction is about 0.13 on this
corpus, and it varies by network. The bridge to a no-truth tier is therefore
reported as failed as an unbiased estimator and usable only with that
correction carried.

### A.2 Resolved on the ledger, 12 September, after the rerun of A.0b

528 frame rows on 27 networks, eight arms, 5,140 evaluations. Five networks
above 500 nodes were excluded as registered; a sixth in the frame, the
109-node network, carries no row because every one of its chain components
is above the questionnaire's twenty-node cap, so no supplier was asked about
it. The primary supplier can certify 229 of its 528 rows; on 155 its own
graph says the treatment cannot cause the outcome, on 144 the query is not
amenable. Those are statuses on the frozen frame, not filters.

**L4 NOT SUPPORTED as registered.** The chance arm's false-claim fraction is
**0.479 [0.398, 0.559]**, near one half as predicted. Against the primary
supplier's 0.329 [0.184, 0.494] the intervals overlap, and on the outcome
that matters, the committed set invalid at the truth, 0.378 [0.295, 0.460]
against 0.288 [0.155, 0.435] overlap as well. The registered test asked for
separation and did not get it. A paired contrast on the 24 shared networks,
which the registration did not name and is reported as post hoc, puts the
chance arm 0.143 [0.006, 0.287] above the model in false claims and 0.082
[-0.020, 0.176] above it in invalid sets. The model asserts something better
than a coin on what it says; whether that buys a valid set more often is not
resolved on 24 networks. The decision rule applies as written: the elicited
arm is not shown to be informative on this substrate, the ledger says so, and
the abstract leads with the identifiability step and the audit.

**L5 HALF SUPPORTED.** No knowledge certifies **68 of 528 rows, 0.129**,
above the registered one in ten. Every knowledge arm certifies more: 229,
206, 243, 197 for the four elicited conditions, 649 of 1,444 for chance, 209
and 294 for the two audit arms. The step is three to four times the empty
arm, which is what the prediction was about; the threshold was set too low
by three points and is reported as missed.

**L6 SUPPORTED, and structural.** After the two fixes of A.0b the claim
certificate is dangerous on **no row of any arm**, 1,524 certifiable
deployment rows and 503 audit rows, and detects every truth-invalid
committed set on every deployment arm. The free rule's dangerous rate is at
least as high everywhere: equal at zero on the primary supplier, 0.009 to
0.015 on the other four. This is not a statistical finding and the paper
must not present it as one. If the committed set is invalid at the truth and
the arm made `w` false claims, retracting those `w` claims leaves a graph the
truth extends, on which the criterion fails, so `r_claim <= w` whenever the
search reaches depth `w`, and a search censored below `w` certifies a floor
below `w`. The zero is what a sound criterion and an audit that asks the same
question must return; it was returned once the audit did. The informative
contrast is the hop instrument on the same rows: dangerous on 29 of 1,524
deployment rows, 0.015 to 0.036 by arm, all above the rule-of-three floor
where the count allows it, for the reason A.0b gives. The naive rungs are far
worse: counting the closure over-certifies on 0.22 to 0.34 of rows, counting
the claims on 0.17 to 0.34.

**L7 NOT RESOLVED.** The measure registered as "safe moves" was the mean of
`r_claim - 1` on rows with an exact radius; on the first table it gave 0.07
for the truthful arm against 0.15 for the model, the opposite order to the
one predicted, and it discards every censored row, which favours arms whose
sets break early. The table now reports the moves the certificate licenses
on every certifiable row, the exact radius less one or the censored floor,
capped at three; this was defined after the first read and is said so. On
it the truthful arm certifies **1.086 [0.767, 1.469]** moves per row, the
model 1.044 [0.752, 1.388], chance 1.057 [0.769, 1.378]; paired, truthful
minus model is +0.018 [-0.099, +0.168] and model minus chance +0.031
[-0.150, +0.214]. The three are the same number within noise under either
measure. The registered alarm, chance certifying more than the model, does
not fire; neither does the order. What this says is what the protocol's
question asked: the certificate counts how many claims stand between the
committed set and its failure, which is a property of the asserted set's
shape around the query, not of its truth. Truth is what the audit bucket
measures, and there the arms separate; the move count does not. The paper
states the certificate as that and not as a proxy for accuracy.

**L10 SUPPORTED.** On the 1,390 rows where both radii are exact and finite,
leverage `r_hop / r_claim` has **median 1**, is above 1 on 0.138 of rows, and
reaches **12** on the oracle arm and 9 on the deployment arms; above 5 on
0.006 of rows. Heavy-tailed in the registered sense, with a thin tail: eight
rows.

**Not registered: the paired instrument ladder.** The protocol's B6 asks for
the claim certificate against the other instruments on the same rows, paired
by network, scored on validity and sharpness. Validity is absolute: the claim
certificate over-certifies on no row on any arm; the free rule does on three
to six rows on four of the five deployment arms and on none for the primary
supplier; the hop radius on three to ten rows on every arm. Sharpness is the
false-alarm rate on truth-valid rows, claim minus the other instrument. Against
the free rule the claim certificate is sharper by **0.206 [0.053, 0.380]** on
the primary supplier, 0.188 on its order replicate, 0.246 on scrambled names,
0.308 on chance, each interval below zero; on the second family 0.129
[-0.038, 0.294], underpowered. Against the hop radius the claim certificate is
less sharp by 0.144 [0.025, 0.297] on the primary supplier and by 0.12 to 0.19
on the others, every interval above zero. The three instruments therefore sit
where the protocol said they would: the free rule is safe and blunt, the hop
radius is sharp and unsafe, the claim certificate is safe and sharper than
the safe alternative, on every arm, with the textbook tier in or out. That is
the sentence the abstract is written from, and the table prints the verdict
per arm rather than pooling them.

**Not registered, read off the same table.** The false-claim fraction of the
elicited arms on the certifiable rows, 0.29 to 0.39, matches the
knowledge-level audit of A.1 on the whole questionnaire. The committed set is
invalid at the truth on 0.23 to 0.34 of certifiable rows under elicited
knowledge, on 0.38 under chance, on none under truthful or oracle knowledge
and on none with no knowledge. The instrument replicate certifies 206 rows
against the primary's 229 on the same networks; the second family 243; the
scrambled run 197. The oracle arm, whose graph is fully oriented, has the
highest share of `r_claim = 1`, 0.769, which is the E4 result again: the
committed set sits one claim from failure most often when the knowledge is
complete and true.

### Errata, 13 September (appended; nothing above edited)

Found while writing the apostila of 13 September, against the files named:
- The header's «Written 13 September» should read 12 September: the sweep, the tables and
  Appendix A.2 are dated 12 September and the file's own timestamps agree. Nothing in the
  predictions depends on the date.
- A.1 (L1, commission by tier) calls T3 «the gene-expression networks». Only ecoli70 and arth150
  are gene networks; magic-niab and magic-irri are plant-genetics networks from MAGIC populations
  (`results/substrate/TABLE2_SUBSTRATE.md`). Read «the four Gaussian genomics networks».

### Errata to the errata, 13 September (appended)

- The first erratum above is itself wrong about the date. The session transcript places the
  writing of this file on 11 September at 15:18 UTC, five minutes after the primary supplier's
  elicitation started and 37 minutes before the first sweep; the header's «13 September» should
  read 11 September. Appendix A.2 was written on 12 September. The ordering the header asserts
  (written before any other supplier answer or ledger row was read) holds on that timeline.
- A.0b does not record the severity-ordering fix that `docs/PROTOCOL_RUN.md` says is dated here.
  Recorded now: on 12 September the audit's severity label was decided validity first, because
  the DAG-level optimal-set formula returns the empty set where the truth orients the query away,
  so an invalid committed set equal to it had been labelled benign. The sweep was rerun and diffed:
  5,140 rows, only the severity column changed, on 113 rows; no certificate, bucket or verdict
  moved. On the current 14-condition ledger the same situation covers 206 rows.
- The A.2 paragraph on the paired instrument ladder calls the free rule «safe». It is dangerous on
  no row of the primary supplier and on 3 to 6 rows in four of the other five deployment arms, 15
  rows in all. Read «as safe as the claim radius on the primary supplier, and nearly safe
  elsewhere».
- §2 states the committed set was «validated against the enumerating definition on 822 rows with
  818 exact agreements and no invalid set». That differential seeded its two half-oriented
  populations with Python's salted `hash()` and could not be reproduced. Rerun on 13 September
  with stable seeds: 834 of 845 agreements, 11 rows where the closed form strictly contains the
  enumerated set, no invalid set, 131 rows recovered. The property §2 relies on, that the closed
  form never returns an invalid set, holds on both runs.
