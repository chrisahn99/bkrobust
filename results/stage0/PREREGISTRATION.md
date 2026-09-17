# Pre-registration — Stage 0

**Written 12 September, before the calibration sweep and before the depth-3
claim radius were run.** The three checks that were already run when this was
written are marked as such and are not predicted here: the free rule against
the committed radius (`stage0_ladder.py`), the single-claim retraction at
coverage 1.0 (`stage0_claim_radius.py`), and the reproduction of a sample of
committed radii through `hybrid.breakdown_radius` (`stage0_reproduce.py`, 30
rows over both dispatch legs, 0 disagreements).

Nothing below may be edited after the first new result file lands in
`results/stage0/`. Corrections go in a dated appendix at the foot.

---

## 1. What Stage 0 is for

The certificate is currently reported in covering hops. The unit an analyst
acts on is the claim they asserted. Stage 0 asks three things, on data already
committed, before any knowledge is elicited:

1. How far apart are the two units, and does the gap have structure?
2. When the knowledge is actually wrong, does the radius over-certify less
   often than the free statistics an analyst already has?
3. Is there a third radius, the distance to losing identifiability, and is the
   left tail the current pipeline deletes empty or structured?

## 2. Definitions, fixed here

- `r_hop` — what `hybrid.breakdown_radius` returns. Covering steps on the
  fixed-skeleton sub-poset. Carries the Conjecture 2 assumption; the error is
  one-sided, so it can only be too large.
- `r_claim` — the size of the smallest subset of the **asserted** set `K` whose
  retraction, with Meek re-closure, makes the committed `Z` invalid at the
  resulting graph. Computed by exhaustive enumeration to depth 3 with a `>3`
  sentinel. No assumption.
- `r_id` — the size of the smallest subset of `K` whose retraction makes the
  query non-amenable, that is, makes the effect not identifiable by adjustment
  at all. Same enumeration, same sentinel. `r_id = 0` where the query is
  already non-amenable.
- `leverage` — `r_hop / r_claim`, the number of covering steps one asserted
  claim carries.
- `phi_1` — the fraction of single retractions from `K` that make `Z` invalid.
- Free rule (`B2`) — the hop distance from `X` to the nearest member of `Z`
  inside `X`'s undirected component, with "no member in the component"
  predicting 1. Computable from the CPDAG and `K`.
- `w_c` — commission count. The number of claims in the asserted set that are
  false in the true DAG. Omission is a separate axis and is not folded in,
  because retraction cannot repair a claim that was never made.

## 3. The calibration sweep

Population: the admissible pairs of `results/axisa3/instances.jsonl`, minus
`pathfinder`. Knowledge is the greedy recovering set, selected at coverage 1.0
and 0.5. **Coverage is a prefix of one seeded permutation per network**, not
the committed even stride, because the stride is not nested: `K(0.25)` is not a
subset of `K(0.5)` on diabetes, ecoli70, magic-irri and munin1, so a committed
coverage sweep is not a retraction sequence. The committed coverage-1.0 rows
are unaffected by the change and are the differential test.

Corruption: `synth.knowledge.flip` at rates 0, 0.1, 0.25 and 0.5, seeded,
several replicates per rate above 0 and one at rate 0 where the draw is
deterministic. Every corrupted set is Meek-closed; a set admitting no
consistent MPDAG is status `knowledge_inconsistent` and is **not** scored as a
failure, because the analyst finds out.

Instruments, all computed from observables only, on the corrupted knowledge:
`B0` certify always, `B1` component size, `B2` the free rule, `B3`
`min(s, |K_G0|)`, `B4` `|K_G0|`, `B5` `|K|`, `B6` `r_hop`, `B7` `r_claim`.

The truth is revealed once, after every instrument has answered, and only to
score. Each instrument's prediction `r` is bucketed against the truth:

| | `Z` invalid at the truth | `Z` valid at the truth |
|---|---|---|
| `w_c <= r - 1` | **dangerous** | held |
| `w_c >= r` | safe | slack |

**Validity** is the dangerous rate. **Sharpness** is the mean certified safe
moves, `r - 1`. An instrument that always answers 1 is perfectly valid and has
zero sharpness, so neither axis alone decides anything.

Reporting: the network is the unit. The intraclass correlation of the `r = 1`
indicator is 0.40 over 25 networks, giving a design effect near 14 and an
effective sample near 60, so every rate carries a network cluster bootstrap and
no per-row interval is printed. Rates are reported over the full denominator
with the undefined-separation stratum included, never conditional on the
measurable subset.

## 4. Predictions

Written before the sweep. Each says what would falsify it.

**P1. `r_claim <= r_hop` on every row.** Certain by the definitions if both are
computed correctly; it is here as an implementation check. A single row with
`r_claim > r_hop` is a bug in one of the two and halts the sweep.

**P2. Re-denominating in claims can only raise the share at 1.** The share of
rows at `r_hop = 1` is 61.1% pair-weighted on the committed data. The share at
`r_claim = 1` will be higher. Falsified if it is lower.

**P3. At coverage 1.0 and corruption rate 0, `r_claim = 1` on every row.** Not
a discovery: `fast_gate` admits a pair only if some single retraction from the
recovering set invalidates the optimal set, which at full coverage is the
definition. Measured at 540 of 540 before this file was written. It is
registered so that the *departure* from it under corruption is the result.

**P4. Under corruption, `r_claim >= 2` appears.** The gate selected on the
uncorrupted recovering set, so nothing forces a corrupted set to fail at one
retraction. Falsified if `r_claim = 1` on every corrupted row too, in which
case the claim unit carries no information anywhere on this substrate and the
paper says so.

**P5. The free rule is the more dangerous instrument.** On the committed rows
`B2` exceeds the committed radius on 16 of 831 and falls below it on 179, so it
is a conservative bound there. Under corruption the prediction is the opposite
direction: because `B2` ignores the knowledge entirely, its dangerous rate will
**exceed** `B7`'s at every corruption rate above 0. Falsified if `B2`'s
dangerous rate is within the cluster-bootstrap interval of `B7`'s, in which
case the certificate buys no validity over a hop count and the contribution
narrows to `phi_1`, the witness and the identifiability step.

**P6. `B4` and `B5` are the least valid.** The closure size and the asserted
count are unrelated to where the failure is. Their dangerous rates will exceed
both `B2`'s and `B7`'s. This is the direction the survival work on the other
branch already found for `|K|`, where the association with survival is negative
rather than absent, so a confirmation here is not a surprise and is reported as
a replication.

**P7. Sharpness orders the other way.** `B0` is maximally sharp and maximally
dangerous; `B7` is the least sharp of the informative rungs. The claim the
paper can make is a validity-sharpness frontier, not a single winner, and the
question is whether `B7` sits on it.

**P8. `r_id > r_claim` wherever both are finite.** Losing validity is easier
than losing identifiability, since an invalid set is still a set. Falsified by
one row with `r_id < r_claim`, which would mean the amenability test and the
validity test disagree about what a retraction does.

**P9. The left tail is not empty.** Among the 518 rejected rows, all
non-amenable at their own `G0`, `r_id = 0` by definition. Among admissible
rows, `r_id` will take at least three distinct values. Falsified if `r_id` is
constant, in which case it carries nothing and is dropped.

**P10. Leverage is heavy-tailed.** The ratio `r_hop / r_claim` will have a
median near 1 and a maximum in the double digits, driven by networks whose
recovering set is a single claim. `paths` asserts one claim and its committed
radius is 14. Falsified if the ratio is within `[1, 2]` everywhere.

## 5. Decision rules

- **P1 or P8 falsified** — halt. These are implementation invariants, not
  hypotheses, and a violation means one of the three radii is wrong.
- **P4 falsified** — the claim unit is uninformative on this substrate. Report
  it, and the elicited arm becomes the only place the certificate can be
  measured at all.
- **P5 falsified** — the certificate does not beat the free rule on validity.
  The paper reports the equivalence with its interval and leads with the
  identifiability step and the audit instead. This is the outcome the sweep is
  built to be able to reach, and it is reportable rather than fatal.
- **P9 falsified** — `r_id` is dropped from the paper.

## 6. What Stage 0 cannot decide

The knowledge here is the true recovering set, corrupted mechanically. It is
not elicited, so nothing in Stage 0 licenses a claim about a real knowledge
source's error process, and every row is stamped as an audit-arm row. The
corruption is uniform and independent, which the matched cross-arm work on the
other branch shows is not the more detectable error process; a correlated
operator is available there and is not run here. Sample sizes are effective
sizes near 60, so an equivalence verdict will be reported as an interval and
labelled underpowered rather than as a null.

---

## Appendix A — outcomes, 12 September

The calibration sweep ran to completion after this file was written: 8,553
scored rows over 24 networks from 16,496 draws, 3,122 s. Of the corrupted rows,
**963 have a committed adjustment set that is invalid at the truth**, which is
the outcome everything here is scored against. Verdicts below, prediction by
prediction, in the order they were registered.

**P1 `r_claim <= r_hop`. HELD.** Zero violations in 8,553 rows. The invariant
stands and the sweep was not halted.

**P2 the share at 1 rises under re-denomination. SUPPORTED.** 61.9% of rows
have `r_hop = 1`; **93.2%** have `r_claim = 1`. Re-denominating in claims moves
almost the whole population onto the smallest certificate.

**P3 `r_claim = 1` everywhere at full coverage and zero corruption. HELD**, as
it had to be: 667 of 667 uncorrupted rows.

**P4 `r_claim >= 2` appears under corruption. SUPPORTED.** 560 rows reach 2, 37
reach 3, across every corruption rate (72 at 0.1, 267 at 0.25, 221 at 0.5). The
claim unit is not degenerate once the knowledge is not the truth, which is the
reason the elicited arm can measure anything at all.

**P5 the free rule is the more dangerous instrument. DIRECTION CONFIRMED,
MAGNITUDE UNRESOLVED.** The free rule over-certifies on 8 of 963 truth-invalid
rows (0.0010) against 0 of 963 for both radii. The direction is as predicted,
but the network cluster-bootstrap interval on the free rule is [0, 0.0038] and
contains zero, so the miss-rate gap alone is not resolved at 24 networks. All
eight misses are on one network, `Kampen_2014`, which is what a clustered
interval is built to notice.

What **is** resolved is the rest of the comparison, and it is the stronger
result. Against the free rule the hop radius detects every truly invalid row
rather than 99.2% of them, raises fewer false alarms (0.404 against 0.498), and
certifies **2.6 times more safe moves** (0.792 against 0.305). It is better on
all four axes at once, not traded off. The worry this sweep was built to test,
that the certificate is an expensive re-derivation of a hop count, is answered:
that holds on truthful knowledge, where the gate pins everything at 1, and it
fails once the knowledge is wrong.

**P6 the count baselines are the least valid. SUPPORTED, and more starkly than
predicted.** The closure size, the asserted count and certify-always all have a
detection rate of **exactly zero**: on 963 rows where the set was invalid at the
truth, they warned on none. Their over-certification rate is 0.122 against 0.001
for the free rule and 0 for the radii. This replicates, by a different route and
on a different endpoint, the finding on the other branch that the claim count is
not an inert baseline.

**P7 sharpness orders against validity. PARTIALLY SUPPORTED, and the stated
form was wrong.** The prediction named certify-always as maximally sharp; its
sharpness is not a number, since it certifies everything, and the sharpest
finite instrument is the closure size at 11.5 moves, followed by the asserted
count at 6.3. Those are the two that never warn, so the frontier is real: the
instruments that certify the most moves are the ones that catch nothing. Among
the instruments that do detect, the ordering is the hop radius 0.792, the free
rule 0.305, the claim radius 0.083. The claim radius is the least sharp of the
informative rungs, as predicted.

**P8 `r_id >= r_claim`. HELD.** Zero violations. Losing validity is never
harder than losing identifiability on this corpus.

**P9 the identifiability radius is not constant. SUPPORTED**, at the floor of
what was predicted: it takes three values, 1 on 7,974 rows, 2 on 523, 3 on 37.
It is informative but it tracks the claim radius closely, and it does not earn a
separate column in the paper on this evidence.

**P10 leverage is heavy-tailed. SUPPORTED.** Median 1.0, maximum **14.0**, and
812 rows above 2. One asserted claim carries up to fourteen covering steps.

### Two things that were not predicted and are worth carrying

**Most corrupted knowledge announces itself.** 1,366 of 16,496 draws produce a
knowledge set with no consistent MPDAG at all, so Meek closure fails and the
analyst finds out before computing anything. These are recorded as a status and
are never scored as failures. The direction agrees with the self-revealing-error
result on the other branch.

**The enumeration bound bit where it was expected to.** 2,694 draws were
recorded `extensions_intractable` at a cap of eight undirected edges in the
analyst graph, and only 19 rows carried a censored radius. The cost is in
enumerating the extensions of a corrupted analyst graph, not in the retraction
search: the median row settles after closing a single subset.

### Scope, unchanged

The knowledge here is the true recovering set corrupted mechanically at declared
rates. It is not elicited, every row is an audit-arm row, and nothing above
licenses a claim about a real knowledge source. The corruption is uniform and
independent, which the matched cross-arm work on the other branch shows is the
*more* silently damaging process and the less detectable one; a correlated
operator exists there and was not run here.
