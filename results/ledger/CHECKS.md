# The protocol's falsifiable checks, as run

Each entry names the check in the protocol, the pre-registered direction, and
what the ledger returned. Numbers are printed by the script named beside each.

## A3 (a) — the cascade test (`experiments/cascade_test.py`)

Pre-registered: the greedy recovering set compels more edges per assertion than
an elicited set, paired by network, interval excluding zero; if not, the claim
that the committed generator "cascades unrealistically" is withdrawn as a
justification for the redesign, which then rests on the invariant alone.

Result: closure per assertion, oracle minus the primary 7B supplier, **+1.13
[−0.26, +2.86]** on 27 networks, larger on 14 of 27. **Not separated.** The
claim is withdrawn as written. Against the 72B model the difference is +1.28
[+0.43, +2.50] and against the 70B second family +1.62 [+0.32, +3.27], both
separated: the larger suppliers assert many claims that compel little, the
oracle set and the 7B assert few that compel much. Truthful orientations on the
7B's own edges compel less than its false ones, −0.63 [−1.52, −0.02].

## A1, second rule — the illustration gate (`experiments/funnel_tables.py`)

Report the maximum radius and the count at r ≥ 3 twice, with and without the
textbook tier. Under elicited knowledge the maximum claim radius is **3** on the
primary supplier and most arms, with or without the textbook tier, and **2** on the 72B and the 70B
arms (corrected 2026-09-13 from "3 on every arm", against `TABLE0_FUNNEL.md`); r ≥ 3 occurs on 8 of the primary
supplier's 229 certifiable rows, 6 outside the textbook tier. The certificate's
measured range on real structure under elicited knowledge is **{1, 2, 3}** and
every sentence about larger radii is a sentence about a hop count.

## A5 (a) — the Lever 0 admission count (`experiments/funnel_tables.py`)

Pre-registered threshold 384 certifiable rows whose analyst graph keeps an
undirected edge, per group. On the committed corpus the count was 0 of 831.
Under the primary supplier it is **112 of 229** certifiable rows; 79 of 285
under the 72B; 388 of 649 under chance; 0 of 294 under the oracle set. On the
elicited arms it stays under the threshold (the chance arm, at 388, is above it; corrected
2026-09-13 from "per arm"), and the abstract's scope sentence is therefore
kept: on this frame about half the certifiable queries are answered on a graph
the knowledge only partly orients.

## Table 8 — the decision table (`experiments/decision_table.py`)

The first-failing retraction set is a single claim on **141 of 157** exact
rows of the primary supplier. In the audit, the named claim was false on 68 of
157; when it was false the committed set was invalid at the truth on 58 of 68,
and when it was true the set was valid on 85 of 89. The sentence "the estimate
holds unless this claim is retracted" points at the claim that in fact broke
the set on 58 of the 62 invalid rows.

## A4 (b) — the monotonicity audit (`experiments/monotonicity_audit.py`)

Pre-registered: failure at a retraction set implies failure at every superset;
one counterexample halts. **5,000 sampled pairs on 469 ledger rows across the
elicited, truthful and oracle arms, zero violations**; the one-sided 95 percent
upper bound on the violation rate is 0.06 percent. Theorem 1 as implemented
stands, and with it the minimal-correction-set reading of the certificate.

## A1 — sortability, the M-bias fixture, the selection node, the Sachs anchor (`experiments/sortability_probe.py`, `results/substrate/DECISIONS.md`)

Var-sortability on the four Gaussian networks, n = 5,000, ten seeds: 0.62
(arth150), 0.71 (ecoli70), 0.59 (magic-irri), 0.33 (magic-niab); R-squared
sortability 0.25, 0.38, 0.61, 0.76. None is near one and the two criteria do
not track each other, so neither ordering shortcut explains what a discovery
method finds on this tier. Every other network is written as undefined with its
reason. M-bias is included as an explicit five-node fixture with its latent
common causes measured, with the written argument that this makes it a
calibration fixture rather than the M-bias trap; the selection indicator of
Didelez 2010 is recovered as ordinary and its 48 rows named as a stratum; the
Sachs file equals the 17-arc network distributed with bnlearn and differs from
the 20-edge Tetrad ground truth by one reversal and three additions, so the
paper names the object it uses and never calls it a consensus.

## A3 — the ancestral granularity (`experiments/ancestral_block.py`, `results/elicit/PREREGISTRATION_ANCESTRAL.md`)

Order elicitation yields tiers and required edges by construction, and a
forbidden edge on a skeleton edge is the required edge the other way; the one
granularity left is ancestral claims on non-adjacent pairs. Asked on 156 pairs
at skeleton distance two across 32 networks, the primary supplier asserted
ancestry on **6**, answered "neither" on 45 and declined 105. All six reduce
soundly to orientations (the reduction fires only when every possibly directed
path shares its first or last undirected edge), giving ten orientation claims
of which **four are false**. Merged into the order-derived set on the four
networks concerned, the condition certifies the same rows with the same
buckets and differs by one true claim. All three registered predictions were
falsified: this supplier does not volunteer ancestral knowledge, and the
granularity adds nothing measurable on this frame. Reported as a negative
result with its count.

## A6-0c — the source-overlap probe (`experiments/source_overlap_probe.py`, `results/elicit/PREREGISTRATION_SOURCE_OVERLAP.md`)

Asked cold, from the variable list alone, to name the publication behind each
network and its declared exposure, outcome and adjustment set: the 7B primary
supplier names the source on **0 of 29** networks (it declines on 27 and
fabricates twice, at high confidence); the 8B second family names it on **0 of
29** and fabricates on all 29. The declared exposure is matched on at most one
of ten scorable networks, the outcome on four, at about chance. **No network is
recited by either model.** With the scrambled-name replicate this closes the
contamination screen on these suppliers from both sides: no recitation is
detected, and no domain knowledge beyond what the variable names give away.
Two instrument defects were found and dated in that file's appendix: an
invalid JSON grammar for the decline token, and a process-salted permutation
seed, replaced by a stable one.
