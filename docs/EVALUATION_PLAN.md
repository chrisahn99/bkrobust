# Evaluation plan for the ICLR 2027 submission

Working note, 11 September. What I checked on the committed data, what the
evaluation looks like, and what runs before the 25 September deadline. The
abstract and the author list close on 18 September.

Every number below is printed by `experiments/stage0_ladder.py` or
`experiments/stage0_claim_radius.py`. Both read `results/axisa3/` and write
`results/stage0/`. The second one skips `pathfinder` (85-node component, three
censored rows).

## What the committed data say

**The gate fixes the radius in claim units.** `fast_gate` admits a pair only
if retracting one claim from the recovering set, with Meek re-closure, makes
the optimal set invalid. At coverage 1.0 the asserted set is the recovering set
and `G0` is the true DAG, so that condition is the definition of a radius of 1
in asserted claims. Measured on the 540 coverage-1.0 rows: 540 of 540 fail
under one retraction. Among them the committed hop radius is 1 in 377, 2 in
102, 3 in 43 and 4 to 14 in 18. The network with radius 14, `paths`, asserts
one claim; retract it and the 15-node chain un-orients. The hop count measures
the length of a Meek cascade, not a number of errors. The median asserted set
has 5 claims. The closure size `k_g0`, which the reports call the knowledge
size, has median 12.

**The free rule is a conservative bound.** Predict the separation where it is
measured and 1 where it is undefined. Against the committed radius it is exact
on 636 of 831 rows and within one on 813. It over-certifies on 16 rows, all
with radius 1 and separation 2 (barley, win95pts, Schipf_2010, magic-niab,
asia), and under-certifies on 179. Of those 179, 69 sit in the undefined
stratum, where 69 of 368 rows have radius 2 and the rule has nothing to say.
The session-5 law transfers as `r >= min(s, k_g0)`: equality on 337 of the 463
measurable rows, radius larger on 110, smaller on 16. `r` never exceeds
`k_g0`.

**The coverage sweep is not nested.** The stride in `select_knowledge` picks
different indices at 0.25 and at 0.5. `K(0.25)` is not a subset of `K(0.5)` on
diabetes, ecoli70, magic-irri and munin1. One pair is admissible at coverage
1.0 and 0.25 and not at 0.5. On the 105 pairs admissible at all three
coverages the radius is unchanged in 73, lower at lower coverage in 32 and
higher in none.

**Two counts.** The 831 rows are 543 distinct pairs. The `o_g0_not_identified`
rejections are non-amenable at `G0`, so no adjustment set can be certified on
them under any definition of the optimal set. The instances that die at partial
coverage die because of which claims the stride keeps.

**phi_1.** The fraction of single retractions that break the optimal set runs
from 0.03 to 1.0 on the coverage-1.0 rows, median 0.30. It varies where the
claim radius does not.

## What this changes

The certificate is priced in the analyst's asserted claims (`r_claim`). The
hop radius is reported beside it as the geometry, and the ratio
`r_hop / r_claim` is the leverage of a claim: how many data-compelled
orientations hang on it. Safe moves counted in hops are not a headline.

A claim radius above 1 can only be measured on knowledge that did not come
from the truth. The elicited arm is therefore where the paper's numbers come
from, and the committed rows become the oracle control below the line.

The free rule is reported as a lower bound with its 16 exceptions, not as a
competitor to beat. Instruments are compared on two axes, validity (rate of
over-certification once the truth is revealed) and sharpness (certified safe
moves given validity), paired within instance and bootstrapped over networks.
The intraclass correlation of the `r = 1` indicator is 0.40 over 25 networks,
so the effective sample is about 60 and no per-row interval is printed.

The `survival-and-pareto` branch already settles three things I had planned to
run. Efficiency is invariant under truthful knowledge, so knowledge buys
identifiability and then robustness, never precision. Correlated errors are
more self-revealing than uniform ones, and uniform errors fail silently about
1.7 times more often. Strict dominance over SHD and `|K|` holds on 7 of 9
strata on the conservative endpoint and is not claimed. The plan builds on
those results instead of repeating them.

## The design

Invariant. The radius is a function of `(C_hat, K, X, Y, Sigma_hat, n)`. The
true DAG draws the samples and audits the certificate afterwards. It is never
on the computation path. Every row carries a hash of the CPDAG source and of
the generating DAG; where the two are equal the row is in the audit arm.

Frame. One candidate set per network, frozen before any knowledge is
collected. `X` lies in or next to a chain component of size at least 2; `Y` is
a possible descendant of `X`. Every knowledge-dependent condition is a status
column on that frozen list, never a filter, so every knowledge source reports
the same denominator. Two panels: the oracle CPDAG on the 25 yielding
networks, labelled as such, and an estimated CPDAG (PC and GES at three sample
sizes, ten seeds) on the four Gaussian networks, which is the only panel that
carries a deployment claim.

Knowledge. A questionnaire over each chain component that contains a candidate
treatment, asking for a causal order over the component rather than for edge
directions, since pairwise questions produce cycles. Sources above the line: a
language model at temperature 0, prompt hashed and cached, in two model
families and two presentation orders so that instrument variance and source
variance are separate columns, and with real and with scrambled variable names
so that recitation shows up as invariance; a random control matched to the
model's answers on size, granularity and distance from the treatment; the
empty set as a status line. Below the line: the truthful recovering set, its
corrupted versions, and the committed rows.

Radius. `r_claim` by exhaustive retraction subsets to depth 3 with a `>3`
sentinel; `r_hop` from the existing search with the dispatch leg printed on
every row; `phi_1`; and an identifiability radius, the number of retractions
until the effect stops being identifiable by adjustment, which is monotone and
is defined on the rows the current pipeline rejects. Polynomial GAC gate only.
Wall cap per instance, censoring recorded as a status.

Audit, opened only after every deployment row is frozen and hashed. Commission
and omission counted separately; the consistent-but-false fraction as the
reporting object; the four buckets (dangerous, held, safe, slack) with the
rule-of-three floor printed beside the dangerous rate; clustering of wrong
assertions tested in graph space and in prompt position; invariance under
scrambled names.

## What the nine pages hold

Figure 1: one real instance, the first-failing retraction set named in domain
terms, and the hop count that missed it. Table 1: the applicability funnel and
the ladder. Table 2: the certificate ledger, deployment sources above the
rule and controls below it. Table 3: the decision table. Table 4: worst-case
radius against average-case survival on the same real instances. Everything
else goes to the appendix.

## Order of work

12 September. Environment on Python 3.12 with ortools. Stage 0 on the
committed rows (done, `results/stage0/`). Five small patches, each with a
differential test against the 831 rows: nested coverage by prefix of one
seeded permutation; a three-valued verdict for the optimal set; the
`g0_undirected_edges` field set in every branch; the dispatch leg on every row;
`not_amenable` split from `optimal_not_identified`.

13 to 16 September. The frame and its population table; the questionnaire;
elicitation; the matched control; the radius sweep; the ladder on the elicited
knowledge; survival and cross-arm corruption on the real instances; the audit;
the sampling band on the four Gaussian networks.

17 September. Title, abstract, authors.

19 to 24 September. Writing. Content freeze on the 23rd.

What I would like from the other branch: the N = 1000 rerun of the nine
strata; the survival and cross-arm code run on the 831 real rows; merging
`survival-and-pareto` into `main` so the paper reads from one branch.

Deferred past the deadline: knowledge harvested from source publications, a
human analyst round, the estimated-CPDAG panel on the discrete networks, a
second corpus, sortability probes, the M-bias expansion, the RoCA port.

## Methods paragraph

`paper/sections/evaluation-protocol.tex` holds the paragraph as drafted. One
sentence needs the two-panel clause: the frame is drawn from an estimated
CPDAG only where data can be simulated.
