# Evaluation plan for the ICLR 2027 submission

Working note, 11 September. What I checked on the committed data, what the
evaluation looks like, and what runs before the 25 September deadline. The
abstract and the author list close on 18 September.

Every number below is printed by a script in `experiments/`: `stage0_ladder`,
`stage0_claim_radius`, `stage0_reproduce`, `stage0_calibration` and
`stage0_differential`. All of them read `results/axisa3/` and write
`results/stage0/`, and `stage0_report` assembles the canonical
`results/stage0/numbers.json` that `STAGE0.md` is generated from, so no figure
is typed twice. They skip `pathfinder` (85-node component, three censored
rows). Predictions were written first, in `results/stage0/PREREGISTRATION.md`.

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

## What has run, as of 12 September

Environment, and two things that do not install on their own: `ortools` is
imported at module level by `hybrid.py` and is declared nowhere, so a clean
install cannot compute a radius; and the suite needs 3.11 or later to collect at
all. Both are in `results/stage0/ENV.md` with the test counts.

The gate first: thirty committed rows recomputed through `breakdown_radius`
across both dispatch legs, zero disagreements. Nothing else was run until that
passed.

The five patches landed in `benchmarks/measure.py`, all additive, all with the
committed behaviour as the default: a `nested` selection mode beside the stride,
a three-valued verdict for the optimal set, `g0_undirected_edges` and `k_g0`
assigned on every branch, the dispatch leg and the search budget as columns, and
`not_amenable` split out of `o_g0_not_identified`. Guarded by
`tests/benchmarks/test_measure_selection.py` and by
`experiments/stage0_differential.py`.

The differential says exactly what moved, and it found a second silent default.
Every committed field on every admissible row reproduces. The only differences
are `k_g0` on rejected rows, where the committed value was a dataclass default
of zero that had never been assigned, the same defect as `g0_undirected_edges`
and in the same record. Every sampled row carrying the old label comes back
`not_amenable`, which confirms through the pipeline what a standalone diagnostic
had found. Under the nested selection the coverage-1.0 rows are untouched, and
at partial coverage about a fifth of rows change status in both directions:
instances the stride killed as non-amenable become admissible and a smaller
number go the other way. Which claims the selection keeps decides whether the
instance survives, and the stride was choosing them alphabetically.

The prior-art gate ran and it retires a claim. `docs/PRIOR_ART.md` has the
eight verdicts. The short version: the smallest subset of asserted constraints
whose retraction breaks a property is a **minimal correction set**, named and
studied in the satisfiability literature since the 1990s, and we should cite it
rather than name it. Nothing was found that applies it to causal graphs,
background knowledge or adjustment-set validity, so the instantiation is where
the contribution has to live. Continuous-parameter relatives, stability radius,
radius of robust feasibility and distance to ill-posedness, are not our object.

## Phase B, what has run

**The frame is built and it is drawn from the CPDAG alone.** Treatment in or
beside a chain component, outcome among its possible descendants, nothing else
filtered: **80,132 candidate pairs over 33 networks**, against 543 admissible
pairs over 25 under the true-DAG gate, an expansion of 148 times. Every one of
the 543 committed pairs lies inside the frame and none outside it, so the frame
strictly contains the population every previous result was measured on and no
committed number is orphaned. Eight networks yield a frame and yielded no
committed instance at all; those are the ones the gate removed. The frozen
sample is 648 rows at a cap of twenty per network with the inclusion
probability recorded. `experiments/frame_build.py`, `results/frame/`.

**The decision table exists, and it is the column no baseline can produce.** On
the 540 committed rows at full coverage, the first-failing retraction set is
named in the analyst's own variables: the estimate for occupation on hormone
therapy holds unless age to occupation is retracted; for alcohol on diabetes
unless alcohol to smoking is; for bronchitis on dyspnoea unless smoking to
bronchitis is. Two things about that set are worth the paper.

It is **unique on 411 of 540 rows**, median one and at most three, so the
sentence is well posed rather than one of many.

And on all 540 rows the retraction costs **identifiability**, not merely
validity: the claim whose removal invalidates the committed set also makes the
query non-amenable, so there is no stratum here where the set breaks and another
valid set survives. Knowledge is buying identifiability rather than a better
set, which is the same conclusion Phase 2 of the other branch reached from the
efficiency side.

**No instrument ever reads below the claim depth.** Against the depth an analyst
can act on, the free rule reads high on 108 rows and exact on 432, the hop
radius high on 163 and exact on 377, the closure size high on 530 and exact on
10, with widest gaps of 6, 13 and 44. None of them ever reads low. That is not
an error any of them makes; it is the unit mismatch, and it is the reason the
certificate has to be re-denominated before a number goes in front of an
analyst.

**The elicited arm ran, on local models, and the ledger exists.** The
questionnaire asks one causal order per chain component of at most twenty
nodes, 540 rendered items over 32 networks, frozen and hashed first. Four
supplier conditions on the local GPU host at temperature zero with a cache: a
seven-billion-parameter model on real names, the same model on a second
presentation order and on scrambled names, and an eight-billion model of a
second family. The Meek repair dropped 43 of the primary supplier's 214 claims
and fired on six of 32 networks. The sweep put eight arms on 528 frame rows
over 27 networks, 5,140 evaluations, and the audit opened only after the
rows were frozen. `results/ledger/TABLE4.md` is the protocol's Table 4 and
the paper's Table 2; the ten registered predictions are scored in
`results/ledger/PREREGISTRATION.md`.

What it says. The primary supplier can certify 229 of its 528 rows; on 155 its
own graph orients the query away and on 144 the query is not amenable. It is
wrong on a third of what it asserts, 0.329 [0.184, 0.494], and its committed
set is invalid at the truth on 0.288 [0.155, 0.435] of the certifiable rows;
chance on the same edges sits at 0.479 and 0.378. Paired over 24 networks the
model beats chance on commission, +0.14 [+0.01, +0.29], and is not shown to
beat it on validity, +0.08 [-0.02, +0.18]. Knowledge buys identifiability
first: the empty arm certifies 68 rows, every knowledge arm three to four
times that.

The certificate in claims is dangerous on no row of any arm, 1,524
certifiable deployment rows, and that is a property of the criterion rather
than a finding: retracting the false claims leaves a graph the truth extends,
on which the set fails, so the radius is at most the number of false claims.
The hop instrument on the same rows is dangerous on 29, between 0.015 and
0.036 by arm, the unit gap measured under elicited knowledge; counting the
closure over-certifies on 0.22 to 0.34 of rows and counting the claims on
0.17 to 0.34. The certificate's move count does not order truth, model and
chance, 1.09, 1.04 and 1.06 per row: it measures where the committed set sits
in the asserted set, not whether the asserted set is true. The paper says so.

The pre-registered halt rule fired on the first read and found two things in
the repository. The DAG-level validity oracle is Pearl's back-door criterion,
which forbids every descendant of the treatment; the certificate decides with
the generalised adjustment criterion, and the two disagree on four of 1,169
committed sets. The audit now uses the generalised criterion written
separately in `benchmarks/audit.py`. And the ladder's constraint model is
built over every ordered node triple, so it never returns above about 150
nodes; the hop instrument walks exhaustively where the retractable set is
small and reports a censored floor where it is not. Both are dated in the
pre-registration's appendix.

**Four larger suppliers ran on the cluster the same day**, on the same
frozen questionnaire, with ten predictions registered first and scored in
`results/elicit/PREREGISTRATION_SCALE.md`. The 72B model of the primary
family answers what the 7B declined, 281 assertions against 171, and is wrong
on fewer of them, 0.233 against 0.325 by network mean, paired −0.10 [−0.22,
+0.02]; its committed set is invalid at the truth just as often, 0.298 against
0.288, because more claims give a wrong one more chances to sit on the query.
The 7B stays the primary row; the 72B is the scale replicate. The claim
certificate stays dangerous on no row, 3,014 certifiable deployment rows; the
hop instrument's over-certification falls from four rows to one and does not
close. Scrambled names cost the 72B accuracy in direction only. The
compelled-edge bridge is unbiased at 72B and biased at 7B. The third family
is the outlier on the high side.

Not run: the estimated-CPDAG panel, so every deployment row is stamped
`via_Chat_only`; the harvested and human arms.

## Order of work

12 September. Environment on Python 3.12 with ortools. Stage 0 on the
committed rows. Five small patches, each with a differential test against the
831 rows: nested coverage by prefix of one seeded permutation; a three-valued
verdict for the optimal set; the `g0_undirected_edges` field set in every
branch; the dispatch leg on every row; `not_amenable` split from
`optimal_not_identified`.

12 September, continued. The frame and its population table; the
questionnaire; elicitation on local models; the matched control; the radius
sweep; the audit; the ledger. Done, as above.

13 to 16 September. The paired instrument ladder on the elicited knowledge;
survival and cross-arm corruption on the real instances; the sampling band on
the four Gaussian networks; a larger supplier family on the cluster.

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
