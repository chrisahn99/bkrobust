# Pre-registration — the ancestral block

**Written 12 September, after `results/elicit/questionnaire_ancestral.json` was
built and hashed (sha256 `fc66f983a1d986f611ec8d30176ebaf5a762646da99fd4256f337c48a6dfe2fd`:
45 items over 32 networks, 156 non-adjacent pairs at skeleton distance two,
drawn from 638 candidates under a cap of eight per network) and before any
model answer to it was requested or read.** The reduction rule and its unit
tests (`src/bkrobust/knowledge/ancestral.py`, `tests/knowledge/`) were written
and passing before this file; they are not edited after the answers arrive.

The block asks the primary supplier of the ledger, the seven-billion
instruction-tuned model on real names at order seed 0, "Is A a cause of B,
possibly through other variables?" for pairs the causal-order questionnaire
could not ask, because they are not adjacent. Each answer is `a causes b`,
`b causes a`, `neither` or `DECLINE` with a quarantined confidence. An
assertion is a directional answer; `neither` and `DECLINE` are not assertions.
An assertion *reduces* when the sound rule converts it into at least one
orientation claim on an undirected edge of the CPDAG (unique undirected first
edge of every possibly directed path away from the ancestor, or unique
undirected last edge into the descendant); otherwise it is `not_reducible`.
Nothing is guessed.

Nothing below is edited after `results/elicit/knowledge_ancestral.json` is
written. Outcomes go in the appendix.

## Predictions

**N1. The model asserts ancestry on more than half the pairs.** Assertions
(`a causes b` plus `b causes a`) exceed half of the 156 pairs asked. The order
questionnaire drew an assertion on 171 of 369 adjacent pairs (46%); a question
that permits mediation is easier to say yes to, so the rate rises above one
half. Falsified if assertions are at or below 78.

**N2. Fewer than half the assertions reduce soundly.** Distance-two pairs
inside chordal components usually have more than one possibly directed path
between them, so the end edges diverge and the rule stays silent. Falsified if
at least half the assertions carry a reduced orientation.

**N3. The reduced claims are no worse than the order-derived ones.** The
commission of the reduced orientation claims, false at the truth over reduced
claims pooled across networks, is at or below **0.325**, the network-mean
commission of the same supplier's order-derived claims recorded in
`results/elicit/knowledge_audit.json`. The rule fires only where the class
leaves one path end, so a wrong reduced claim needs a wrong ancestral claim,
and ancestry through a mediator is the easier fact. Falsified if the pooled
commission exceeds 0.325. The per-network mean is printed beside it and is
secondary, because most networks contribute one or two reduced claims.

## Not predicted

How many ledger rows the merged arm `D_LLM_PLUS_ANC` certifies, and whether
its committed set is invalid at the truth more or less often than `D_LLM`'s.
Those are reported, not predicted: the merged arm differs from `D_LLM` only on
the networks where at least one claim reduced, and the sweep is run only there.

---

## Appendix A — outcomes, written after `knowledge_ancestral.json`

The 45 items were answered by the primary supplier at temperature zero and
seed zero, every answer cached under the hash of model and prompt; every item
parsed. Counts, from `results/elicit/knowledge_ancestral.json`:

| asked | asserted | reduced | not reducible | neither | declined | not reached | invented |
|---|---|---|---|---|---|---|---|
| 156 | 6 | 6 | 0 | 45 | 105 | 0 | 2 |

**N1 falsified.** Six assertions on 156 pairs, 3.8%, against a predicted
majority. The model declined 105 pairs and answered `neither` on 45. A
question that permits mediation was not easier to say yes to; the supplier
that asserted 46% of the adjacent pairs asserted 4% of the distance-two pairs.
Whole networks were declined outright (andes, ecoli70, munin3, paths, Kampen)
and others answered `neither` throughout (hailfinder, water).

**N2 falsified.** All six assertions reduced, four at both ends and two at the
first edge only, ten orientation claims in all; no assertion was
`not_reducible`. The prediction assumed divergent paths; the pairs the model
was willing to assert sit on components small enough that the end edges are
unique. Whether the rule stays silent on the larger components is not tested
by this run, because the model declined them.

**N3 falsified.** Four of the ten reduced orientations are false at the
truth, pooled commission 0.400, above the 0.325 reference; the per-network
mean over the five networks with a reduced claim is 0.333, also above it.
Four of the six ancestral assertions themselves are false at the truth
(arth150 twice, asia, munin), so the reduced claims were wrong because the
ancestry was wrong, exactly as the rule guarantees and not in the direction
the prediction hoped. Ten claims is too few to separate 0.400 from 0.325; the
sign is recorded, the interval is not claimed.

**The merged arm, `D_LLM_PLUS_ANC`.** Of the ten reduced orientations, three
were already asserted by the order, one contradicted it (asia, where the
order was right and the ancestry wrong) and was dropped, and six were new. Of
the six new, four were dropped for Meek consistency against the order's own
frozen set: the three in arth150, and the one in child, where the reduced
claim was true at the truth and the order's single claim it contradicted was
false, so the direct-question-wins policy kept the wrong one. Two claims were
added, one each in Shrier_2008 and munin; munin lies above the sweep's
500-node cap, so the merged arm differs from `D_LLM` on exactly one ledger
network.

**Not predicted, reported.** `ledger_sweep.py --arms D_LLM_PLUS_ANC
--knowledge-extra results/elicit/knowledge_ancestral.json` on the four
networks under the size cap with a reduced claim (Shrier_2008, arth150, asia,
child), written to `results/ledger/parts_anc/`: 80 frame rows, 32 certifiable,
41 query reversed, 7 not amenable, the same three statuses row for row as
`D_LLM` on the same networks. The merged arm differs from `D_LLM` by one claim
on Shrier_2008, true at the truth: median claims 3 against 2, certified moves
0.5625 against 0.5313 on the 32 rows, false claims over asserted 0.052
against 0.060 pooled and 0.250 against 0.250 by network mean, committed set
invalid at the truth on 3 rows against 3 (all in child, all from the order's
own false claim), hop bucket 29 held, 2 safe, 1 dangerous under both arms,
claim bucket 27 held, 3 safe, 2 slack under both. The ancestral block changed
the ledger by one true orientation on one network and no certificate.
