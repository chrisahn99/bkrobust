# The protocol as run

What each stage of the evaluation protocol became in code, what it measured,
and where it departs from the protocol as written. Numbers are in
`results/ledger/TABLE4.md` and `results/stage0/STAGE0.md`, both generated;
this file says what produced them.

## The invariant, and how it is enforced

The radius is a function of the CPDAG, the knowledge, the query, and nothing
else. The true DAG enters twice: to draw the CPDAG the frame is built on, and,
after every row is frozen, to score it. Each ledger row carries
`true_dag_on_path` with one of the four values in `core.conventions`:
`via_Chat_only` on every deployment arm here, because the CPDAG is the true
DAG's and no discovery run replaced it; `via_K` on the audit arms, whose
knowledge was read off the truth. Rows of the estimated-CPDAG panel, whose
CPDAG was learned from a sample, carry `false`. A row's
`chat_source_sha256` and `generating_dag_sha256` are printed so the stamp can
be checked rather than trusted.

## Stage 0 — the ladder and the unit

`experiments/stage0_*.py`, reported in `results/stage0/STAGE0.md`. The free
hop-count rule against the committed radius, the claim-unit radius on the
committed rows, the reproduction gate, the corrupted-knowledge calibration, and
the differential for the five patches to the measurement module.

## Stage 0B — prior art

`docs/PRIOR_ART.md`. The retraction operator is a minimal correction set. Cited,
not claimed.

## Stage 1 — the committed set in polynomial time

`mpdag_criterion.optimal.committed_adjustment_set`. The optimal-set formula of
Henckel, Perković and Maathuis lifted to a maximal PDAG, with the canonical set
as fallback and a verdict saying which was committed to. Validated by
`experiments/lever0_differential.py` against the enumerating definition on three
populations: fully oriented, half-oriented truthful, half-oriented corrupted.
Agreement 834 of 845 where both exist, no invalid set anywhere, and 131 rows
recovered that enumeration could not reach. Where the enumerating set is a
strict subset of the closed form's, 11 rows, the closed form is valid but not
the minimal optimal set, so any efficiency claim drops there and the verdict
column says so. (Rerun on 13 September with stable seeds: the first run seeded
its two half-oriented populations with Python's salted `hash()`, so its 818 of
822, 110 and 4 could not be reproduced; the conclusion is unchanged.)

A query the analyst's graph orients away, so that the outcome is no longer a
possible descendant of the treatment, is a verdict of its own rather than a
certificate. The frame drew the pair from a CPDAG in which the effect was
possible; the knowledge said otherwise. That is information about the
knowledge and it is counted.

## Stage 2 — the frame

`experiments/frame_build.py`. Treatments in or beside a chain component of at
least two nodes, outcomes among their possible descendants, nothing else
filtered, per-network cap of twenty by seeded stride, inclusion probability
recorded. 80,132 candidate pairs over 33 networks; all 543 committed pairs
inside it. Every condition the old gate used as a filter is a status column on
this frozen list.

## Stage 3 — the questionnaire

`experiments/questionnaire_build.py`. One question per chain component that
contains a frame treatment, capped at twenty nodes, asked as a causal order and
read back per pair as `a->b`, `b->a` or `DECLINE` with a quarantined confidence.
Two naming conditions: real names with the data dictionary an analyst would
have, meaning a domain sentence that never names the benchmark and the
variable's states where the file has them; and scrambled tokens with a neutral
domain sentence. Two presentation orders under recorded seeds, with a startup
assertion that the two render distinct prompts. A compelled-edge control block
per network, rendered exactly like an open question. 120 components, 369 open
pairs, 540 rendered items, frozen and hashed before any supplier saw them.

### Knowledge granularity

The protocol's questionnaire was to elicit four granularities: required edges,
forbidden edges, tiers and ancestral constraints. Three of them are covered by
the causal order by construction, under the orientation-only scope. A causal
order over a component is a tier assignment with one variable per tier, and a
tier assignment forbids exactly the edges pointing from a later tier to an
earlier one, so on the fixed skeleton its whole content is the orientation of
each adjacent pair from earlier to later, which is what is read back. A
required edge `a -> b` on a skeleton edge is the orientation claim `a->b`; on a
non-adjacent pair it would assert an adjacency the data rejected, which is
outside the scope since the skeleton is the CPDAG's. A forbidden edge `a -> b`
on a skeleton edge is the same claim as the required edge `b -> a`, because the
pair is adjacent and the edge has exactly two orientations; on a non-adjacent
pair it forbids nothing. Ancestral constraints are the one granularity the
order does not express, since they name a directed path rather than an edge,
and they are elicited by a separate block, below.

`experiments/ancestral_block.py`. For every chain component of the frozen
questionnaire, the pairs of component nodes that are not adjacent and lie at
distance two in the component's undirected skeleton, capped at eight per
network by a seeded draw, asked as "Is A a cause of B, possibly through other
variables?" with the answer grammar `a causes b`, `b causes a`, `neither`,
`DECLINE` and a quarantined confidence, rendered like the open items (real
names, the domain sentence, the variable states, order seed 0). The frozen
questionnaire is untouched; the block is `results/elicit/questionnaire_ancestral.json`
with its own hash. An ancestral assertion becomes orientation claims only when
it entails them on the CPDAG: when every possibly directed path from A to B
shares its first edge and that edge is undirected it is oriented away from A,
and when every such path shares its last edge and that edge is undirected it is
oriented into B; every other assertion is `not_reducible` and counted, never
converted. The rule is in `bkrobust.knowledge.ancestral`, sound by
construction and incomplete by design, with a unit test where it fires and one
where it must not. The reduced claims are frozen and hashed before the truth
scores them, and the merged arm `D_LLM_PLUS_ANC`, the primary supplier's
orientation set plus the reduced claims, is read by the sweep from the second
knowledge file through `--knowledge-extra`. Predictions were registered before
any answer was read, in `results/elicit/PREREGISTRATION_ANCESTRAL.md`.

As run: 156 pairs over 32 networks. The primary supplier asserted ancestry on
6, answered `neither` on 45 and declined 105, whole networks among them. All
6 assertions reduced, to 10 orientation claims, of which 4 are false at the
truth (pooled commission 0.400, against 0.325 for the same supplier's
order-derived claims; four of the six ancestral assertions were themselves
false). Merged into `D_LLM_PLUS_ANC`: 3 of the 10 were already asserted by
the order, 1 contradicted it and the order kept the edge, 4 were dropped for
Meek consistency against the frozen order set (one of them true, where the
order's claim it contradicted was false), 2 were added. The sweep on the four
networks under the size cap with a reduced claim, `results/ledger/parts_anc/`:
80 rows, 32 certifiable, unchanged from `D_LLM`; the merged arm differs by one
true claim on one network, certified moves 0.5625 against 0.5313 on the same
rows, pooled commission 0.052 against 0.060, invalid sets 3 against 3, every
audit bucket identical. All three pre-registered predictions were falsified,
each recorded in the file's appendix. The granularity is now elicited and
reducible; on this supplier it is mostly declined.

## Stage 4 — the suppliers

`experiments/elicit_run.py`. Every call at temperature zero and seed zero,
cached against the hash of model and prompt, never repeated once cached. The
models run on the local GPU host through its HTTP API. The primary supplier is
a seven-billion-parameter instruction-tuned model; the second family is an
eight-billion-parameter model from a different lineage. Neither is called an
expert anywhere.

Four conditions, each a ledger column: the primary supplier on real names; the
same model on the second presentation order, which isolates instrument variance;
the second family on real names, which isolates source variance; and the primary
supplier on scrambled names, which is the recitation probe. The order the model
returns is the primary instrument and the per-pair directions are read against
it; a stated direction that contradicts the model's own order is counted and
the order wins. `DECLINE` is an answer; a pair never mentioned is `not_reached`.

Per network the orientations are Meek-closed against the CPDAG. A set with no
consistent closure is recorded as such, that is the rejection rate the protocol
asks for, and is then repaired as a practitioner would: drop claims in
increasing confidence until the closure succeeds, and record the count. The
post-repair set is frozen and hashed before any radius is computed.

Four larger checkpoints answer the same frozen questionnaire on the cluster
through vLLM, one user message per item through each model's own chat
template, temperature zero, seed zero, the same token budget:
`Qwen2.5-72B-Instruct` at both presentation orders and on scrambled names,
`Qwen2.5-32B-Instruct`, `Llama-3.3-70B-Instruct` and `gemma-3-27b-it` on real
names. `experiments/elicit_export.py` writes the prompts, `experiments/remote/`
holds the runner and the batch script, `experiments/elicit_import.py` puts the
answers into the same cache under the same key, and `elicit_run.py` refuses to
call a cluster model itself. Their predictions were registered before any
answer was read, in `results/elicit/PREREGISTRATION_SCALE.md`, and scored
in its appendix. The primary supplier of the ledger remains the seven-billion
model; the chance and truthful audit arms stay matched to its edges. The short
version: the 72B model asserts 281 pairs where the 7B asserted 171, is wrong
on a smaller fraction of them with a paired interval that does not exclude
zero, and its committed set is invalid at the truth on the same fraction of
certifiable rows, 0.298 against 0.288. Scale bought assertions and accuracy
per assertion, not validity. The claim certificate is dangerous on none of
the 3,014 certifiable deployment rows; the hop instrument on one row of the
72B arm and on none of the gemma arm. The truth-free compelled-edge bridge,
biased at 7B, is unbiased at 72B.

Two of the protocol's suppliers are not here: knowledge harvested from source
publications, and human analysts. Both are human hours, and both are deferred
with that reason.

## Stage 5 — the sweep

`experiments/ledger_sweep.py`. For every frozen frame row and every arm: the
analyst graph, the committed set with its verdict, and then, from observables
alone, the retraction certificate to depth three under a subset budget, the hop
radius with its dispatch leg, the fraction of single retractions that break the
set, the free rule, the closure size and the claim count.

The hop radius is the upward walk of the hybrid, exhaustive when the analyst's
graph has at most twelve retractable edges and to depth three otherwise, under
a budget of four thousand closures; the ladder, whose model is built over every
ordered node triple, runs only on graphs of at most sixty nodes and with a
five-second wall per rung. A walk that stops on a budget carries the depth it
checked, so the row states `r_hop >= depth + 1` and is marked inexact, and the
audit scores that instrument on the lower bound. The first sweep hung in the
ladder's build on the 223-node network and was discarded unread; the
pre-registration's Appendix A.0 records the change. The row is frozen there. Two
controls that assert nothing from the truth sit beside the suppliers: a
chance-level arm on the primary supplier's own edges with orientations drawn
uniformly among Meek-consistent choices, and an empty knowledge set. Five
networks above 500 nodes are excluded for cost and printed as excluded.

## Stage 6 — the estimated-CPDAG panel

`experiments/pest_panel.py` and `experiments/pest_table.py`, reported in
`results/pest/TABLE7.md`, with its predictions and their outcomes in
`results/pest/PREREGISTRATION.md`. On the four networks with fitted Gaussian
parameters, samples at three sizes and ten seeds, PC with a Fisher-z test at
level 0.01 and GES with BIC; one reference run per network at the largest
sample becomes the CPDAG the frame, the questionnaire, the primary supplier and
the sweep are built on. Every row carries `false`.

What it says, and it bounds the whole certificate. The truth extends the
analyst's graph on none of the 244 certifiable rows: every learned skeleton
misses or adds edges, and a quarter of the compelled orientations point the
wrong way. The committed set is invalid at the truth on 28 of the primary
supplier's 38 certifiable rows, against about a quarter on the oracle panel.
The claim radius is dangerous on 14 of them and the hop radius on the same 14;
with no knowledge at all the certificate is dangerous on 14 of 19. The
guarantee that the claim radius cannot over-certify holds because retracting
the false claims lands in a class that contains the truth, and a learned CPDAG
whose class does not contain the truth voids that premise. The certificate
prices error in the analyst's claims given a correct equivalence class. It
says nothing about discovery error, which on these networks at five thousand
samples invalidates more committed sets than the claims do. The one network
whose learned skeleton is nearly exact is the one where neither instrument is
dangerous on any deployment row (the two audit arms are dangerous on 3 and 4 rows there;
qualified 2026-09-13).

## Stage 7 — the audit

Same script, after the freeze. Commission is the number of asserted claims
false at the truth; omission is the number of orientable edges in the
treatment's component the arm never asserted. Whether the committed set is
valid at the truth, its severity against the true optimal set, and the
four-way bucket for each instrument's certified radius: dangerous when the
certificate covered more errors than were made on an invalid set, held, safe,
slack.

Validity at the truth is the generalised adjustment criterion on the
generating DAG, implemented in `benchmarks.audit` over networkx and sharing no
code with the criterion package the instrument uses. The repository's DAG
oracle is Pearl's back-door criterion, which forbids every descendant of the
treatment; it is sufficient and not necessary, and on the ledger it rejected
four committed sets the adjustment formula accepts, each containing a
descendant of the treatment on no causal path to the outcome. Its verdict is
kept as a column beside the audit's. Severity is decided validity first: the
DAG-level optimal-set formula returns the empty set where the truth orients the
query away, and a committed set equal to it can still be invalid, so equality
with that set is read as benign only once the set is valid. A censored claim
radius enters the bucket
as the floor it certifies, `r_claim >= d`, not as an unbounded one; both of
these were found when the pre-registered halt rule fired on the first read of
the table, and both are dated in the pre-registration's appendix. Two audit arms below the rule: truthful orientations on the supplier's
own edges, and the greedy recovering set.

## Departures from the protocol as written, each with its reason

- **Four experiment scripts seeded with Python's salted `hash()`** (found 2026-09-12 and
  2026-09-13): the chance arm of the sweep, the Lever 0 differential, the corrupted-knowledge
  calibration, the superseded frame-status run, and the questionnaire builder. All now use a
  SHA-256 seed. The sweep and the differential were rerun and their numbers updated; the
  calibration is being rerun. The questionnaire frozen on 12 September used the salted seed for
  its scrambled-name permutation and its compelled-edge sample, so it cannot be regenerated byte
  for byte; it stays the authoritative input the suppliers answered, and the builder now refuses
  to overwrite a frozen questionnaire whose hash differs, writing a rebuild beside it instead.

- **The CPDAG is the oracle's on every row of the main ledger.** The
  estimated panel exists on the four networks with fitted parameters, and on
  24 of the 39 networks there is no covariance to run discovery on. That panel
  shows the certificate's guarantee does not survive discovery error, so no
  deployment claim is made beyond the statement that the certificate is
  computable from the learned CPDAG and the claims.
- **The second supplier family is small.** An eight-billion-parameter model,
  because the twenty-billion reasoning model on the same host spends its whole
  budget thinking and returns nothing within its wall. A larger family is a
  matter of downloading one.
- **Harvested and human knowledge are absent**, as above.
- **One network in the frame has no elicited knowledge.** The questionnaire
  asks one causal order per chain component and caps a component at twenty
  nodes; every component of the 109-node network with 122 undirected edges is
  above the cap, so no supplier was asked about it and it carries no ledger
  row. The cap is the questionnaire's, not the frame's, and is printed.
- **The chance-level control is matched to the primary supplier's edges**, as
  the protocol says, and therefore inherits which edges were spoken about from
  a source that itself chose them without the truth. That is the intended
  matching, not a leak.
- **The commission count treats a false claim as one wrong claim.** The
  protocol's operator prices a flip at two atomic moves. The audit bucket
  follows the protocol's own A-CAL-1 definition, which counts false asserted
  claims; the two-move reading would make every instrument look safer and is
  not used.
