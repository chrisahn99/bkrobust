# LOG — audit trail

Append-only. What was launched, what came back, what was integrated, what was
discarded and why.

## Session start

- Branch gate: `git branch --show-current` -> `experiments/synth_graphs_and_heuristics`. Passed.
- Read existing work: `src/bkrobust/demo/` (14 modules), `tests/demo/` (7 files,
  133 tests), `report.md`, `results/breakdown_radius_demo/`.
- Starting state: 4d0de90, clean tree, no `synth/` or `search/` yet.

## T0 — conventions and frozen core

- **T0.1 radius convention resolved.** Chose `r = min distance from G0 to a
  failing graph` (shells `0..r-1` certified clean; safe-moves = `r-1`).
  Rationale in `src/bkrobust/core/conventions.py`. Prior numbers (3,3,2) carry
  over unchanged; the other convention would have shifted every published
  number by one. Asserted element-by-element in
  `tests/core/test_conventions.py::test_shells_below_radius_are_genuinely_clean`.
- **T0.2 core frozen**: `core/{conventions,oracle,spacelib,instance,resultsio}.py`.
  Reuses the demo's already-cross-checked primitives rather than reimplementing
  them. Integration check: the new core reproduces the published radii 3/3/2.
- **T0.3 determinism** extended to the core: 3 hash-seed-invariance tests.
- 152 tests passing, ruff clean.
- Inherited-code fix: `src/bkrobust/data/loaders.py` carried a dead
  `import pandas as pd` under TYPE_CHECKING (scaffold commit c201cfb, never
  used). Removed so the repo stays lint-clean. No behavioural change; that
  module is stubbed.

## B1 — the two structural conjectures

- **Conjecture 1 is a theorem, not a conjecture.** Validity is a for-all over
  represented DAGs, so a failing `D in [G]` is still in `[G']` whenever
  `[G] <= [G']`. One-line proof, recorded in the module docstring. One caveat:
  this implementation defines `is_valid` as False for an EMPTY extension set, so
  an empty-extension graph would break the implication — it does not bite
  because `enumerate_space` admits no such element. Checked, not assumed.
- **Conjecture 2 exhaustively tested, no counterexample so far.**
  - n=3: 150 radius comparisons, 0 violations (exhaustive).
  - n=4: 13,344 radius comparisons over 126 spaces, 0 violations (exhaustive).
  - n=5: running, scope = all CPDAGs with <=6 undirected edges (99% of the 8782).
- **Free validation:** the enumerator reproduces the published counts of
  labelled DAGs (25, 543, 29281) AND of Markov equivalence classes
  (11, 185, 8782) on 3/4/5 nodes. Both pinned in a test.

## B2/B3 — accelerated exact search and differential validation

- Implemented `search/exact.py`: `radius_local_up` does upward BFS with covers
  generated LOCALLY, so the space is never enumerated. Exploits upward-closure
  by never expanding above a known failure. Carries an anytime mode that returns
  a labelled bound rather than a wrong exact answer.
- **Differential validation: 5,304 cases (n=3 and n=4 exhaustive), 0
  disagreements** with full-space BFS. Committed to `results/search/differential.csv`.
- Critical correctness check: locally generated upper covers were compared
  against the space-derived covers element by element and match exactly
  (`test_local_covers_match_space_derived_covers`). Had local generation missed
  a cover, distances would silently come out too large.
- **Honest speedup accounting.** Comparing per-query times with a pre-built
  space flatters BFS, because `local_up` never builds one. On the single-query
  comparison (BFS pays space construction; local_up does not) the speedup grows
  with the undirected-edge count k: 25x at k=2, 81x at k=6, 133x at k=7 on n=5
  CPDAGs. Where many queries share one CPDAG, BFS amortises construction and the
  advantage shrinks — both regimes must be reported.

## Bug found in my own work: results directories were gitignored

`.gitignore` line `/results/*` (added during the earlier demo session, with an
exception only for `results/breakdown_radius_demo/`) silently excluded
`results/synth/` and `results/search/`. Consequence: the pre-registration and
the Axis B differential rows were NOT committed by the commits that claimed
them. Caught when `git commit` reported "nothing to commit" for the
pre-registration, which should have been a new file.

**Correction to the record:** commit 133ab44's message says the differential
rows were "committed to results/search/differential.csv". They were not; they
are committed here instead. The numbers in that message are unaffected — only
the claim about where the file landed was wrong.

Exceptions added for both directories. Nothing was lost; the files were on disk
throughout.

## B2.3/B2.4 — bounds, and the L == U certificate rate

Exhaustive over all n=4 CPDAGs, 792 instances with a finite exact radius:

- `descendant_lower_bound` (L) covers ONE failure mode: some `z in Z` becoming a
  possible descendant of X. **Admissible in every case where it is defined —
  732/792, zero violations** — and **tight (L == r) in 660**.
- `guided_upper_bound` (U) was a valid upper bound in all 792 cases.
- **L == U == r in 660/792 = 83.3%**: the exact radius is certified with no
  enumeration at all in five out of six instances.

Important limitation, stated rather than buried: L is undefined in 60/792
instances (7.6%), where no `z` can be made a descendant of X yet the radius is
still finite. In those the binding failure mode is a back-door path becoming
unblocked, for which no bound is implemented. So L is NOT a complete lower bound
on its own, and `L = UNREACHED` must be read as "this mode gives no bound",
never as infinity. The 83.3% figure is the rate at which the pair certifies
exactly, not a claim that L alone certifies.

## Bug in my own harness: O(|space|^2) check run per combination

The n=5 conjecture sweep stalled for minutes on a single dense CPDAG, twice, at
the same index. Diagnosis: `check_conjecture1` is O(|space|^2) and I was calling
it once per (G0, X, Y, Z) combination. For a 120-element space with ~19k
combinations that is ~276M pair comparisons on one CPDAG.

Two fixes, in order of what actually mattered:
1. (minor) The space-size guard ran after `build_space`, so the cubic
   covering-relation computation happened before the graph could be rejected.
   Now elements are enumerated first and covers built only if the size passes.
2. (the real one) Conjecture 1 is a THEOREM, so checking it guards the
   implementation rather than testing mathematics. It is now sampled at a
   recorded rate instead of run exhaustively.

n=4 re-verified after both changes: 13,344 comparisons, 0 violations, unchanged,
and runtime fell from ~10s to 3.9s.

## B1.2 — Conjecture 2 exhaustive result (n=5 complete)

| n | CPDAGs | spaces | radius comparisons | C2 counterexamples |
|---|---|---|---|---|
| 3 | 11 | 7 | 150 | 0 |
| 4 | 185 | 126 | 13,344 | 0 |
| 5 | 8,782 | 6,030 | **1,971,820** | **0** |

Total ~1.985M comparisons, zero counterexamples. Conjecture 1 sampled alongside:
9,859 checks at n=5, zero violations.

Scope, stated rather than implied: at n=5, 2,616 CPDAGs have no undirected edge
(nothing for knowledge to orient), 131 have more than 6 undirected edges, and 5
had a space above 200 elements. So 136 of the 6,166 CPDAGs that *could* carry
knowledge were skipped (2.2%), all of them the densest. The claim is therefore
"no counterexample over all CPDAGs on at most 5 nodes with at most 6 undirected
edges and a space of at most 200 elements", not "over all CPDAGs on 5 nodes".

## Finding: Meek(cpdag, K) can fall outside `enumerate_space` (0.09% of cases)

The main Axis A run crashed with `source is not an element of the space`. Root
cause is a genuine inconsistency in the inherited validity definition, not a
coding slip.

Minimal captured case:

    cpdag : V1->V4 V2->V0 V2->V4 V3->V0 V4->V0 V5->V0 V5->V4 V6->V0
            V1-V2 V1-V3 V1-V6 V2-V3 V2-V5 V2-V6 V3-V5
    g0    : V1->V4 V1->V6 V2->V0 V2->V4 V2->V6 V3->V0 V3->V2 V4->V0
            V5->V0 V5->V4 V6->V0 V6->V1  V1-V2 V1-V3 V2-V5 V3-V5

`g0` is Meek-closed, keeps every compelled edge, satisfies `[g0] subset [cpdag]`,
and represents **5** consistent DAGs -- yet `is_valid_mpdag(g0)` is False, so
`enumerate_space` omits it and BFS from it raises.

Why: `is_chordal_components` inspects the UNDIRECTED subgraph only. There the
component {V1,V2,V3,V5} is the 4-cycle V1-V2-V5-V3-V1 with no undirected chord.
But V2 and V3 *are* adjacent, via the DIRECTED edge V3->V2. So the cycle is
chordless only if directed edges are ignored, and a legitimate knowledge state
is excluded from the space.

Frequency: **2 of 2235 (0.09%)** Meek closures across four generators.

Decision: gate it, do not silently drop it, and do not change the validity
definition in this session. Relaxing `is_valid_mpdag` would invalidate the
space enumeration that every result so far rests on (including S2's
chain-space-is-6 result and the 1.985M-comparison conjecture sweep) and would
require re-running all of it. The instances are recorded as a rejection reason
so the rate stays visible, and the question is carried to NEXT.md.

## Axis A first run — two gate bugs found, partial results kept

Main grid (864-point pilot, then 8,640-point run) over the three healthy
generators. The run crashed at grid index 4934 but incremental writing preserved
377 rows, which is what that design is for.

Two defects found by inspecting the output rather than trusting it:

1. **`g0_not_in_space` crash** (0.09% of Meek closures) — recorded above.
2. **34/377 rows (9%) had `r_val == 0`**, which the frozen convention forbids.
   All had `Z = {}`: the runner proposed target pairs with no causal path from X
   to Y, so `cn(X,Y)` is empty, `O*` is empty, and the empty set is not valid
   because a back-door path is open. The convention document already anticipated
   this ("r = 0 ... gated out at generation time"); the runner's gate was
   incomplete. Sent back with two new required gate reasons (`no_causal_path`,
   `z_invalid_at_g0`).

Preliminary H1 on the surviving rows, **after excluding the 34 invalid ones**
(340 finite instances, Erdos-Renyi, n=6..9):

    r_val = 1 : 272  (80.0%)
    r_val = 2 :  58  (17.1%)
    r_val = 3 :   8  ( 2.4%)
    r_val = 4 :   2  ( 0.6%)
    UNREACHED :   3

Against the exhaustive n<=4 census (90.9% at r=1) this is a lower saturation
rate, in the direction H1 predicted: mass at r >= 2 grows with graph size. Not
yet a claim -- the run must be redone once the gates are fixed.

Preliminary H5 (middle regime, r_eps > r_val) is strongly epsilon-dependent:
0.0% at eps=0.02 and 0.05, 0.8% at 0.1, 16.3% at 0.2.

## Correction: my transcription of the failing graph was faulty

When sending the `g0_not_in_space` diagnosis to the generator subagent I pasted a
`g0` edge string containing BOTH `V1->V6` and `V6->V1`. The MPDAG constructor
provably rejects mutually directed pairs (verified), so no such graph can exist —
the paste was garbled by me, not produced by the code. The subagent spotted the
contradiction, declined to trust it, and independently reconstructed the mechanism
as a minimal 4-node example. That is the better artifact and the right response.

The underlying phenomenon was real and is unaffected: I had independently measured
`Meek(cpdag,K)` falling outside `enumerate_space` at 2/2235 (0.09%) before sending,
and the subagent's `g0_not_in_space` gate now covers it.

## Axis A main results (clean run, fixed gates)

Grid: 12,960 points over erdos_renyi / scale_free / block, n=6..9, four
corruptions, three rates. **1,064 accepted (8.2%)**, zero run errors, zero rows
with `r_val == 0`.

Rejections by reason (the rate is itself a result):
  empty_set_trivially_valid                6,889
  treatment_not_in_or_adjacent_to_component 3,107
  no_atomic_perturbation_changes_validity    696
  z_not_identified                           620
  k_assumed_inconsistent                     341
  cpdag_too_large_for_bfs                    160
  z_invalid_at_g0                             83

## H1 and H2 headline numbers

Exhaustive census, all CPDAGs on 5 nodes (134,140 instances, 88,700 with a
finite radius):
    r=1  71,840 (81.0%)   r=2  14,100 (15.9%)   r=3  2,760 (3.1%)
    UNREACHED 45,440
Ensembles, n=6..9 (1,050 finite): r=1 78.6%, r=2 17.9%, r=3 2.8%, r=4 0.5%.
Stable across generators (76.9-80.4%) and across n (77.4-81.2%).

**H2: 0 / 134,140 instances at n=5, and 0 / 1,044 at n=4, have a valid
adjustment set strictly more robust than the optimal one.** The
efficiency-robustness frontier is flat in every instance examined.

## CORRECTION: my H2 statistic was wrong, and the frontier is NOT flat

The earlier claim "0 of 134,140 instances have a strictly more robust valid
adjustment set" is **wrong**. It was an artefact of my own analysis code, in two
compounding ways:

1. `h2_frontier` excluded every comparison in which either radius was UNREACHED,
   on the grounds that UNREACHED is not a number. It is not a number, but it is
   not missing data either: a set that never fails anywhere in the fully
   enumerated space is *strictly more robust* than one failing at distance 1 —
   the largest gap there is. Discarding those comparisons discarded exactly the
   cases where the frontier is non-flat.
2. The census capped candidate adjustment sets at 10 sorted by size, which can
   hide a more robust larger set.

Corrected statistic (n=4, all CPDAGs with <=6 undirected edges, no cap,
UNREACHED ordered above every finite radius):

    instances with >=2 valid sets and O* defined : 2,652
    frontier NON-FLAT                            :   816  (30.77%)
    O* STRICTLY BEATEN                           :     0  ( 0.00%)
    O* never fails anywhere (UNREACHED)          : 1,584  (59.7%)

Excluding the empty adjustment set (trivially robust, and gated out of the
ensembles anyway) the non-flat rate is 10.96% of 2,628 instances; O* is still
never beaten.

**The scientific conclusion changes.** The frontier is genuinely non-flat — valid
adjustment sets do differ in robustness — but there is no efficiency-robustness
TRADEOFF, because the optimal (most efficient) set is never strictly less robust
than any alternative. That is a different and more informative statement than
either "the frontier is flat" or "a tradeoff exists".

This also reframes H3: designed adversarial families are not needed to exhibit a
non-flat frontier, since ~31% of ordinary small instances already have one.

## H6 calibration and session close

H6 over 89 instances across the three generators:
  coverage = 1.0000 in EVERY instance (the correctness gate, not a finding)
  conservativeness mean 0.181, range 0.000-0.625
    block 0.236 | erdos_renyi 0.189 | scale_free 0.149
  Spearman(undirected-component size, conservativeness) = 0.564 -- the
  knowledge-intersected component predicts conservativeness, as H6 predicted.
Conservativeness is much lower than the prior report's single-example 0.43-0.47.

Deliverables written: report_synth_and_search.md, NEXT.md, PLAN.md revisions.
277 tests passing, ruff clean.

## Frontier sweep complete (n=4 and n=5)

| n | instances | non-flat | O* strictly beaten | O* never fails |
|---|---|---|---|---|
| 4 |   2,652 |    816 (30.8%) | 0 | 1,584 (59.7%) |
| 5 | 247,080 | 53,880 (21.8%) | 0 | 147,740 (59.8%) |
| total | **249,732** | 54,696 (21.9%) | **0 (0.00%)** | 149,324 (59.8%) |

Across a quarter of a million exhaustively enumerated instances the optimal
adjustment set is never strictly less robust than any alternative valid set,
while the frontier itself is non-flat in over a fifth of them.

## Session PDF report

Built `SESSION_REPORT.pdf` (23 pages) at the repository root as the single
centralised record: Task 0, both axes, every hypothesis verdict including the
failed predictions and the hypothesis that was not run, all nine bugs, the
orchestration history, scope limits, and reproduction commands.

Two gaps closed on the way:

1. **No figures had been produced this session.** The 20 files in `figures/`
   were all from the prior worked-example session. Nine session figures were
   generated from the committed results (`analysis/session_figures.py`), plus
   the speedup measurement, which had been computed ad hoc and never saved --
   so it could not have been plotted or reproduced. It is now in
   `results/search/speedup.csv` (165 queries, all agreeing with BFS).
2. **Two rendering bugs, both caught by actually looking at the output.**
   reportlab does not wrap bare strings in table cells -- it runs them off the
   page and clips them -- so the headline table was silently truncated in the
   first build; cells are now Paragraphs. And the base-14 fonts are WinAnsi
   only, so 9 distinct characters used in the prose (arrow, subset, >=, C-hat,
   QED box, the o-double-acute in Erdos) rendered as empty boxes. DejaVu is now
   registered from matplotlib's bundled fonts, which keeps the document
   portable rather than depending on macOS system fonts.

The PDF's commit table is generated by `git log` at build time, so it lists the
commits up to but not including the commit that adds the PDF itself.

# ==================== SESSION 2 (Axis B deep) ====================

## Task 0 RESOLVED — the chordality test was wrong, not the validity notion

Decisive experiment: compare `enumerate_space(C)` against the set of *reachable
knowledge states* `{Meek(C,K) : K consistent}`, which is what the space is
supposed to model.

Result over all CPDAGs on 3 and 4 nodes with an undirected edge (133 CPDAGs):
7 mismatch, and **every mismatch is one-way — reachable states EXCLUDED from the
space, never the reverse**.

Every excluded graph is:
  - reachable (an explicit K produces it),
  - Meek-closed under R1-R4,
  - non-empty in DAG extensions,
  - and **maximally oriented** — every undirected edge genuinely undecided
    across its extensions, which is the defining property of an MPDAG.

The only test it fails is chordality of the undirected subgraph. That condition
characterises CPDAGs (essential graphs are chain graphs with chordal chain
components); it is not a property of an MPDAG built by adding knowledge, because
knowledge can place a directed edge inside an otherwise-undirected component.
The chord that would make the component chordal is then present but DIRECTED, so
the undirected subgraph reads as chordless.

**Minimal witness** (as small as it gets — one assertion):
    cpdag : V0-V1 V0-V2 V0-V3 V1-V2 V1-V3   (K4 minus V2-V3)
    K     : V0->V1
    state : V0->V1 V0-V2 V0-V3 V1-V2 V1-V3
    meek_closed=True  extensions=5  maximally_oriented=True
    old is_valid_mpdag=False  -> excluded

**Corrected membership test** (fixpoint, verified equivalent to reachability on
1,588 states): `apply_orientations(C, G.directed \ C.directed) == G`.

**The effect is much larger than previously reported, and grows with density.**
The prior session's 0.09% was per Meek-closure. Measured per element of the
space (n=5):
    k<=4 : 0.00% missing,  0/12 CPDAGs affected
    k=5  : 3.85% missing, 12/12 affected
    k=6  : 2.14% missing,  9/12 affected
    k=7  : 8.25% missing, 12/12 affected
    k=8  : 10.19% missing, 12/12 affected
So at k=5 EVERY CPDAG loses elements, and the previous session's sweeps ran to
k=6. Its conclusions are more fragile than they looked; this must be said.

Mitigation chosen: the old `enumerate_space` is left UNMODIFIED so prior
committed results stay reproducible from the code that produced them. New work
uses `search/space_fixed.py`. Everything dependent is re-run side by side.

Good news for the inherited headline: the worked example's CPDAG has no excluded
states (48 elements either way), so the prior report's radii 3/3/2 are unaffected.

## Bug (MINE, RECURRING): results/axisb2/ was gitignored

Exactly last session's bug 1, repeated. `.gitignore` line `/results/*` whitelists
only `results/breakdown_radius_demo/`, `results/synth/` and `results/search/`, so
everything written to `results/axisb2/` this session was untracked — including
`task0_space_membership.json`, which the Task 0 commit message claimed to add.

**Caught by the back-door-bound subagent**, which noticed its own output was
untracked and said so in its report rather than assuming the orchestrator had it
in hand. I had not checked, despite having written the identical mistake up in
last session's report as bug 1. Recording it again, with the aggravating factor
that it was a known failure mode.

Correction to the record: commit 748f2e1 ("Task 0 ...") did NOT commit
`results/axisb2/task0_space_membership.json`. The file existed on disk
throughout; it is committed here. No numbers are affected.

## Session 2 close

Delivered: THEOREMS.md (Conjecture 2 proved modulo Anti-Exchange Case B),
report_axisb_deep.md + REPORT_AXISB_DEEP.pdf (9 pages, 5 figures), updated
NEXT.md, results under results/axisb2/ with manifests, and 5 new figures.

Subagents: two launched, both reviewed rather than trusted.
  - conjecture2 re-run on the corrected space: n=3 and n=4 exhaustive (0
    counterexamples), n=5 partial at 62.4% when its parent stopped. Its
    incremental checkpointing preserved 814,736 comparisons.
  - back-door lower bound: L now defined in 100% of finite-radius instances, 0
    admissibility violations. I re-verified admissibility independently on the
    old-space scope (792 instances, 0 violations, tight 90.9%) before accepting.
    It also caught my gitignore recurrence.

Tests 310 passing, ruff clean across the repo.

## Correction: the n=5 corrected sweep completed after I reported it as partial

I reported the n=5 corrected sweep as 62.4% complete (5,480 of 8,782 CPDAGs,
814,736 radius comparisons), because its incremental checkpoint file was the only
evidence available at the time its parent agent paused. The sweep subsequently
ran to completion:

    all 8,782 CPDAGs, 6,030 spaces, 1,977,820 radius comparisons, 0 counterexamples
    720 space elements added by the Task 0 correction
    13,860 comparisons involved a state the OLD space did not contain

The earlier number was accurate when written and is superseded, not wrong. The
report, THEOREMS.md and the PDF are updated to the completed figures.

The density gap remains open on the same terms: 131 CPDAGs exceed k=6 and 5 have
a corrected space above 200 elements. A closure run over those 136 was under way
at session end; partial results (9 CPDAGs at k=7..10, 18,177 comparisons, 0
counterexamples) are checkpointed in n5_dense.jsonl.

## n=5 density gap CLOSED (135 by sweep, 1 by the proof chain)

The conjecture2 subagent completed the density closure: 136 CPDAGs attempted,
135 run to completion (299,520 further radius comparisons, 0 counterexamples),
1 recorded as infeasible rather than dropped — the complete K5 skeleton, whose
corrected space has 4,231 elements and whose covering relation costs
O(N^3) ~ 7.6e10 operations.

I closed that last one a different way. Anti-Exchange needs only closure
computations, not the covering relation, so it is cheap exactly where enumeration
is not; and the chain AE => Lemma R => Property S => Lemma L => Conjecture 2 is
proved. On the complete K5 CPDAG: 4,231 closed sets, 60,140 closures, 30,840
applicable triples, 0 violations, 10 seconds.

This is a CONDITIONAL closure and is labelled as such: it settles Conjecture 2
there by inheriting the Anti-Exchange dependency, rather than testing Conjecture 2
directly. The distinction matters and is stated in the report.

Net effect: the previous report's most pointed caveat -- "all CPDAGs on <=5 nodes
with <=6 undirected edges" rather than "all CPDAGs on <=5 nodes" -- is retired.

Total corrected-space evidence: 2,290,834 exhaustive radius comparisons across
n=3,4,5, zero counterexamples.

## Bug (MINE): two silent no-op string replaces

While integrating the density-gap closure I updated `report_axisb_deep.md` but the
PDF builder kept its own hard-coded text, so I committed and pushed a PDF that
contradicted the markdown ("the density gap is still open"). Two of my
`str.replace` calls had also silently no-opped earlier because the pattern used a
hyphen where the source had an em dash, and I had not asserted on them.

Caught by verifying the built PDF's extracted text against what I had just
claimed, rather than trusting that the rebuild picked the change up. Fixed by
replacing every block with an explicit `assert` on the anchor before writing, and
by adding a post-build check that greps the rendered PDF for the claims it is
supposed to make. Lesson recorded: a rebuild is not a verification.

---

## Session 3 — the declarative encoding

### Bug (MINE): the closure encoding was too strong, and lost half the spaces

Stage 1 reproduced only **66 of 133** CPDAG spaces. I had encoded each Meek rule
as an implication and dropped R3's and R4's *undirectedness* premises, reasoning
from R1 and R2 — which do tolerate the omission, because the alternatives to
their consequent are a new v-structure or a cycle, both independently forbidden.
R3 and R4 have no such backstop.

Caught by comparing the solution set against the enumerated corrected space
rather than by reading the encoding. Fixed by encoding every rule as a
**forbidden firing configuration**: one clause per premise-matching tuple ruling
out (premises ∧ ¬consequent). Lesson: a "simplification" that is sound for two
members of a rule family is not thereby sound for the family.

### Bug (MINE): the collider clause erred toward under-reporting robustness

Stage 2 disagreed with the oracle on 144 of 26,304 cases. My non-collider clause
required `into_l ∨ into_r`; a collider needs **both** arrows pointing in. So a
chain through a member of `Z` was accepted as a collider, blocked paths read as
open, and non-failures read as failures.

The direction is what matters. This makes radii **too small** — and a spuriously
small E3 radius is exactly the shape of a Conjecture 2 counterexample. Had it
survived into the stress test it would have manufactured a false headline result:
"Conjecture 2 refuted", with a witness that does not exist.

It was caught only because Stage 2 was run as its own layer against a reference
oracle *before any radius was computed*. Had I gone straight from "the encoding
compiles" to "here are the radii", the numbers would have looked plausible.

### Bug (MINE): `UNREACHED` conflated with "budget exhausted"

E1's ladder has a natural top; exhausting it proves no failure exists in the
up-set. E3's `max_k` is a pure budget, so exhausting it proves only
`radius > max_k`. My first version returned `UNREACHED` there, which reads as
"no failure exists". Now records `exhausted_to` and sets `exact=False`. Found
while writing the comparison logic for the n = 4 sweep, not by a failing test —
which is to say it would have shipped if I had only run the tests.

### Not a bug: the hang that turned out to be the result

The first high-`k` sweep stalled with no output. Cause: `local_up` on a
degenerate instance at k = 15 must exhaust an up-set of up to 2¹⁵ states, each
requiring a Meek closure and a DAG-extension enumeration. I killed it and re-ran
with `local_up` in a subprocess under a hard 60 s cap.

The instinct to treat this as an obstacle would have been wrong. Recording the
timeout as a datum — 26 of 49 degenerate instances unfinished at 60 s, against
E1 answering all 49, the worst in 8.745 s — is the session's central measurement.

### Mechanical, but worth recording

A `sed` invocation used `|` as its delimiter while the replacement text contained
`|` in a type annotation, corrupting a dataclass field into `None = None| None =
None`. Separately, a blanket `str.replace` of an import line also rewrote that
line *inside an embedded source string* used by a subprocess test, producing an
unterminated literal. Both were caught immediately by the linter and the test
run, and both were edits made without first reading the line being changed.

### Determinism

Radii, witnesses, walks **and solver conflict/branch counters** hash identically
across `PYTHONHASHSEED` 0, 1 and 12345 with a single worker and a fixed seed.
The counters are in the hash on purpose: matching radii alone would not detect a
search that explored a different tree and landed on the same answer.

---

## Session 4 — the oracle, and the cost of enumeration

### Not a bug, but the session's pivot: the brief's premise was wrong

The brief predicted `local_up` would be oracle-bound, so that replacing the
oracle would move its ceiling substantially. Profiling said 96.1%
enumeration-bound but only **25.6%** oracle — **70.4%** was cover-minimality
enumeration. By Amdahl that caps the criterion at 1.3× on its own.

The reason this mattered is that the measurement was *front-loaded*, exactly as
the brief instructed. Had it been run after building the criterion, the criterion
would have been built to do a job it cannot do, and the 70.4% would plausibly
have gone unnoticed because nothing would have been looking at it.

My own recorded prediction was also wrong, in the same direction but less so: I
predicted a roughly even split and a ~2× cap. The correct move was to measure,
not to reason from the call graph.

### Bug (MINE, in delegation): correctness gated, performance not

I delegated the criterion with an exhaustive correctness acceptance criterion and
**no performance acceptance criterion**. It came back verified on 74,568
exhaustive cases and ~2.4M sampled, 0 disagreements — and **exponential in the
vertex count**, because condition (c) enumerated every simple path between `X`
and `Y`. Up to 2.35 s per call; over a 60 s cap on 12 of 81 instances at n = 12
and 14, on instances where the enumeration oracle it is meant to *replace* never
timed out once.

The component was thoroughly verified and useless for its purpose, and only one
of those two properties was checked. Caught by running the envelope sweep, which
recorded criterion timeouts where there should have been none — and only because
the sweep timed *both* oracles rather than assuming the new one was faster.

Lesson: an acceptance criterion that does not mention the property the component
exists for is not an acceptance criterion for that component.

### Bug (MINE, in the harness): timeouts masquerading as measurements

My subprocess drivers wrote `{"timed_out_at_s": CAP, "seconds": <wall clock>}`
for a censored run, reusing the same `seconds` key that a real measurement uses.
Any analysis filtering on the presence of `seconds` would silently treat a
timeout as a measurement, and my first analysis pass did exactly that — it
reported "criterion completed = 80, criterion timeouts = 11" for the same 80
rows, which is what exposed it. Renamed to `wall_until_timeout_s`.

The same class of error as session 3's `UNREACHED`-versus-budget-exhausted bug:
a censored outcome given the same representation as a real one.

### Two textbook definitions that fail on knowledge-carrying MPDAGs

Both surfaced through differential testing during the criterion's construction,
and both are worth recording because they are the same trap twice.

`possde(X,G) = ⋃_D de(X,D)` is **false** under naive possibly-causal paths: in
`V0→V1, V0−V2, V1−V2` the path `⟨V1,V2,V0⟩` is possibly causal, yet every
extension keeps `V0→V1` because orienting it forward closes a cycle. Restricting
to unshielded paths repairs it (570 → 0 mismatches at n ≤ 4).

Amenability must **not** use that same reduction, because short-circuiting
`x→v→w→y` deletes the undirected first edge that is the whole question. Decided
instead by imposing `x→v` and Meek-closing (702 / 1,380 / 0 mismatches for the
three formulations).

A lemma that holds on CPDAGs need not survive the addition of background
knowledge. This repository has now been bitten by that three times: the
chordality filter in session 2, and both of these.

### The good news, recorded at equal weight

The sharp test the brief designed for E1 **passes**. E1 must witness every
invalidity including the non-amenable kind, which is structurally different from
the two failure modes it encodes; failure would have made radii too *large*,
overstating robustness. 74,568 cases, 0 disagreements against both oracles, with
the non-amenable stratum isolated at 28,464 cases (38.2%) and every one of them
witnessed. Session 3 tested that layer but never isolated the stratum.

---

## Session 5 — correct definitions, real components, saturation

### Notes taken at the start, so no prior report needs re-reading

Session 1 baselines this session is measured against: `r_val = 1` in **81.0%**
of the exhaustive n = 5 census (88,700 finite; r=2 15.9%, r=3 3.1%) and 78.6%
across n = 6…9 ensembles, remarkably stable across generators (76.9–80.4%).
Frontier: `O*` never strictly beaten, **0 of 249,732** instances with ≥2 valid
sets. H5: 0.0% at small ε, 14.1% at ε = 0.2. H6: coverage 1.0000,
conservativeness 0.181. Degeneracy gate rejected **91.8%**, mostly "empty set
trivially valid". H4 never run — now four sessions unrun.

Session 4 envelope: largest undirected component **6** at n = 24; 159 of 256
instances at components 2–4. This is the structural reason the saturation
question could not be answered before.

### Phase structure

Phase 1 (definitions) and Phase 2 (generator) are disjoint and were launched in
parallel. Phase 3 (census) depends on both. Phase 4 is analysis and report.

### Decision: reversing session 4's reconciliation, and why it was right then

Session 4 found the repo's oracle is Pearl's back-door criterion, not the GAC,
and reconciled the MPDAG criterion *toward* back-door by enlarging the forbidden
set to `possde(X,G) ∪ {X,Y}`. That kept four sessions of results mutually
comparable, which was the right call at the time. This session makes the GAC the
definition at both layers and reports old versus new side by side. The back-door
implementations stay in the tree unmodified as second oracles, so nothing prior
becomes unreproducible.

### Caveat recorded on sight, before the generator's own report

The component generator reports **acceptance rate 1.0** (1,584 of 1,584), where
session 1's random Erdős–Rényi draws were rejected 91.8% of the time. That is
not a bug — the generator calls `synth.runner.gate` and the instances pass —
but it is a change in kind that the report must state plainly: **this is a
designed family, not a random sample of realistic CPDAGs.** Session 1 made the
same caveat about generator-specificity for its own families, and it applies
with more force here because the construction is aimed at the hypothesis.

What it buys is a grid that is not confounded: 66 `(c, s)` cells covering
`s = 1…c−1` for `c = 2…12`, so separation varies within a fixed component size
and component size varies within a fixed separation. That is exactly what the
pre-registered joint-table analysis of H7 versus H8 needs.

### Timing hygiene note, recorded when the decision was made

The random-ensemble control and the frontier sweep were run **concurrently**,
which contaminates wall-clock fields in both. That was a deliberate trade: the
headline of each is a distribution of radii and candidate counts, not a cost,
and the night is finite. **The clean runtime-versus-radius numbers are the
census's**, which ran alone. Any cost claim in the report cites the census, and
`results/axisa2/manifest.json` records which files' `seconds` fields are
load-contaminated so nobody later mistakes them for measurements.

### Bug (MINE, twice over) — the sweep was bottlenecked on the gate, not the radius

The random-ensemble control crawled: minutes per instance at n = 15–20. I first
blamed the radius search and re-scoped the sweep twice around it (600 s cap →
120 s cap; then a GAC spot-check instead of a second full radius computation).
Both were reasonable but both were treating a symptom.

Timing a single instance end-to-end showed the truth: **187 seconds, and the
instance was then REJECTED by the gate** — the radius was never computed at all.
`synth.runner.gate` calls `all_valid_adjustment_sets_mpdag`, which enumerates
every subset of V. At n = 20 that is 2²⁰ per candidate pair.

Two separate errors of mine, and the mechanism that caught each:

1. **I set `search_budget=64` on the hybrid**, which means the bounded search
   explores to depth 64 and the E1 ladder — built in session 4 precisely to
   answer the UNSAT side fast — never fires. Caught by reading my own worker
   while looking for the real bottleneck.
2. **I diagnosed by inference rather than measurement**, twice, and re-scoped a
   sweep on each wrong diagnosis. Caught by finally timing one instance instead
   of reasoning about which part was slow. The same lesson session 4 recorded
   about front-loading a measurement, and I did not apply it here until the third
   attempt.

Fix: a `fast_gate` with the same verdict vocabulary, replacing the three
enumeration-dependent checks with equivalents that cost one back-door test each.
**187 s → 0.2 s**, roughly 900×.

**A third error, caught by the differential test rather than by me.** I
documented the replacement's "no atomic perturbation" check as strictly WEAKER
than the original. It is strictly STRICTER: the original sets `sanity = True` if
*any* valid set is perturbable, so restricting to `O` can only reject more.
Measured: 16 disagreements in 46,800 cases (0.034%), every one in that direction.
The excluded instances are ones where `O` is robust to every atomic perturbation
— i.e. `r_val(O) = UNREACHED` — so dropping them removes maximally-robust
instances and biases **against** this session's hypothesis. Conservative, and
recorded rather than buried.

### The sharpest result of the session was found by a subagent rendering a figure

While building the two-way table, the figures subagent noticed and verified that
`r = min(s, |K_{G₀}|)` holds in **780/780** rows at coverage 0.5. I re-derived it
independently across both coverage levels: **1,572 of 1,572, 100.00%**.

That is a better statement than the `r = s` I had. It unifies the two coverage
regimes in one formula and names both binding constraints — how far the
adjustment set sits from the treatment, and how much the analyst claimed to know
— where I had described the coverage-0.5 case only qualitatively. It is a law of
the designed family rather than a theorem, and the report says so, but it holds
with no exceptions across 1,572 instances.

Worth recording as a process note: the instruction to review deliverables rather
than trust them is usually about catching errors. This time the review found
something better than what was asked for.

### Two pre-registered predictions of mine failed, and both are reported

1. **H10.** I predicted `O*` would be strictly beaten on a non-zero but small
   fraction, under 5%. It is beaten **zero** times, across 235,534 GAC-admissible
   candidate sets this session.
2. **H5-revisited (`r_ε`).** I predicted the middle regime would appear once
   radii were large. It does not — 0.0% at every ε and every radius stratum out
   to 8, with mean absolute bias at `G₀` exactly 0 in all 560 instances.

Both are recorded because the pre-registration named them in advance and a
prediction that fails is only useful if it is reported as having failed. The
second is the more interesting: it means there is no intermediate
"slightly biased" regime for the radius to be measuring past, which strengthens
`r_val` as the right cut-point rather than undermining it.

---

## Session 6 — real graphs

### Acquisition, resolved before anything was planned

Egress is partial and the obvious route is dead: **`bnlearn.com` is unreachable**
(`TLSV1_ALERT_PROTOCOL_VERSION`), so `pgmpy.utils.get_example_model`, which
downloads from there, cannot be used. Guessed GitHub raw URLs for `.bif` files
all returned 404 — a reminder that guessing URLs is not acquisition.

What worked: downloading the **pgmpy sdist from PyPI** and reading its package
data directly. It ships 24 `.bif.gz` benchmark networks and 9 applied-paper
`.txt` DAGs. `pgmpy` itself is **not installed** — recent versions pull heavy
dependencies — so the archive is read and parsed here, which also means the
parser is ours and can be validated against an independent count rather than
trusted.

Everything therefore traces to one downloaded artefact with a checksum, which is
what the hard rule requires.

### Ordering as a control

The descriptive table is committed **before** the first radius, so it cannot be
shaped by what the radii turn out to be. And because the full pre-registration
has to wait for Phase 2 (it must fix the network list), the **directional
predictions were committed separately, before Phase 2 finished** — see
`results/axisa3/prediction_before_looking.md`. A pre-registration written after
the descriptive table could not have protected those.

### Phase 3 first look — the factual question resolves positive

Early rows (small networks only; the large ones run last because the driver
sorts smallest-first so results arrive early):

- **Separation ≥ 2 in 46.4% of measurable pairs**, against Erdős–Rényi's 14.2%,
  with a maximum of 7. That is **above** my pre-recorded 15–40% prediction, so the
  answer to the session's factual question is positive and my prediction was too
  pessimistic.
- **`r_val = 1` in 50.2%** of admissible instances, against ~80% in random
  ensembles, with a maximum radius of 14.
- **The law holds in 73.1%**, inside my predicted 70–100% band. The exceptions
  are the informative part and they run **mostly one way**: 48 cases of radius 2
  where `min(s, |K_{G₀}|) = 1`, 40 of radius 3 where the law says 2. So on real
  structure the true radius is often **larger** than the designed family's law
  predicts — the spine construction was conservative, because breaking the single
  path between `X` and the adjustment set does not always invalidate it when real
  graphs carry redundant blocking routes. Five cases run the other way and need
  separate characterisation.

### The tractability ceiling is in the GATE, not the radius — measured, not inferred

The measurement sweep stalled twice, and both times I timed one instance end to
end rather than reasoning about it (session 5's lesson, applied first this time
rather than third).

**Stall 1 — the gate's validity predicate.** `fast_gate`'s perturbation loop
called `is_valid_adjustment_set_mpdag`, which enumerates DAG extensions. By
Theorem 14 the polynomial GAC predicate coincides with it exactly on
`Z = O(...)`, which is what the gate passes. Verified rather than assumed: 2,240
pairs across child, insurance, asia, sachs and water, **identical verdicts,
8.2× faster** (`results/axisa3/gate_predicate_swap.json`).

**Stall 2 — `optimal_adjustment_set_mpdag` at low coverage.** It enumerates the
extensions of `G0`, and at coverage 0.25 on a large real component `G0` still
carries many undirected edges. Bounded by the same constant session 5 used
(`MAX_G0_UNDIRECTED_FOR_EXTENSIONS = 12`) and given **its own status**,
`o_g0_extensions_intractable` — a measurement limit is not a structural
rejection, and conflating them would misstate how often the framework applies.

**Stall 3 — pathfinder, and an optimisation I did not make.** The gate's
perturbation loop calls `apply_orientations` once per knowledge edge: on
pathfinder that is 122 from-scratch Meek closures on a 109-node graph, measured
at **29 s per pair**, so its 11,772 pairs would take ~95 hours.

The obvious fix is to test the one-level upper covers of `G0` instead, which is
**746× faster** on pathfinder. **I tested it and it is wrong.** Retracting one
knowledge edge and re-closing can land strictly above a cover, so the loop tests
a superset of what the covers test; the two disagreed on child and pathfinder.
Recorded because the speedup was tempting and only the differential test stopped
it — the same discipline that caught session 5's `fast_gate` direction error.

So pathfinder and the largest components are recorded as **censored with a
measured reason**, which is a better result than session 4's ceiling: this time
the ceiling is where the method actually breaks, and the break is in the
degeneracy gate rather than in the radius computation.

### The verifier caught a miscount before publication

The acquisition table said the sdist ships **11** dagitty applied-paper DAGs. It
ships **12**; I had silently folded the M-bias exclusion into the count of what
was shipped, which conflates "what the archive contains" with "what was usable".
Corrected to "12 (11 usable)" with M-bias named in the same row. Caught by
`session6_verify`, not by rereading.

Small, but exactly the class of error the verifier exists for: a number that is
defensible in isolation and wrong in its column.

### Three of five pre-recorded predictions were wrong

| | predicted | actual |
|---|---|---|
| median undirected fraction | 15–45% | **10.7%** |
| networks with max component > 6 | more than half | **11 of 39** |
| separation ≥ 2 | 15–40% | **45.4%** |
| law agreement | 70–100% | 72.8% ✓ |
| hybrid fails at component 15–40 | — | **never failed; the gate did** |

The descriptive pair were wrong in the same direction and the headline one in the
other, so the session's answer is stronger than I expected while the framework's
applicability is narrower. The H14 miss is the most useful: I predicted the wrong
component entirely. The radius handled an 85-vertex component in 0.16 s; the
degeneracy gate — a screening convenience — is what does not scale.

### A figure review that caught a substantive error, not a cosmetic one

The first render of the separation figure labelled its x-axis *"separation
between the optimal set and the nearest admissible rival"*. That is not what
separation is — it is the graph distance **inside `X`'s undirected component,
from the treatment `X` to the nearest member of `O(G₀)`**, and nothing to do with
rival sets. The figure would have carried a wrong definition of the session's
central quantity into the report.

Fixed, along with a legend sitting on top of two annotations and per-bar labels
colliding on the small bars. Recorded because the reviewing instruction is
usually framed as catching bugs in code, and this was a bug in *exposition* —
the numbers were right and the caption was wrong, which is the harder kind to
notice.

### Bug (MINE, in the report's reasoning) — a "confirmation" that pointed the wrong way

I wrote that the coverage sweep was a "genuine out-of-sample confirmation" of
session 5's law, on the grounds that a less-informed analyst should have a
smaller radius, and cited `r = 1` falling 69.8% → 51.1% → 34.0%.

**That is backwards.** A falling share at `r = 1` means radii moving *up*, which
is the opposite of what a shrinking `|K_{G₀}|` predicts. I had taken a number
that moved and assumed it moved my way.

The resolution, once measured properly: the three coverage levels are **not the
same population**. Admissible instances collapse 543 → 182 → 106, because a
less-informed analyst more often cannot identify `O(G₀)` at all, and the
survivors have systematically larger separation (median 1 → 2). The `r = 1` share
is a composition effect and is evidence for nothing.

What *does* support the law is the ceiling, which I had not looked at: **max `r`
falls 14 → 4 → 4** as median `|K_{G₀}|` falls 12 → 6 → 4, and the inequality
`r_val ≤ min(s, |K_{G₀}|)` holds within every coverage stratum.

**Caught by the figures subagent**, which refused to render the monotone story
the brief and my text implied and said so explicitly. Worth recording twice over:
the review instruction is usually framed as catching bugs in delegated code, and
this is the second time this session it caught a bug in *my* exposition instead —
the first being a figure axis carrying a wrong definition of separation. Both
were cases where every individual number was right.

### A limitation found on review: the ball is centred on the truth

`select_knowledge` draws the analyst's claims from `knowledge_to_recover`, which
reads the *true* orientation off the ground-truth DAG. So the coverage sweep
varies **how much** the analyst knows and never **whether they are right**. At
coverage 1.0 — 543 of 831 admissible instances — `G₀` **is the ground-truth DAG
exactly**: zero undirected edges left, zero orientations contradicting the truth.
At 0.5 and 0.25 it is under-determined but still never contradicts the truth.

Session 1 had exactly this machinery and used it: the worked example's scenario B
is a *flipped* assertion. Session 6 imports none of `flip`, `omit`, `compound` or
`tiered`. The gap was mine and the main report did not flag it.

**Measured rather than conceded.** Same admissible pairs on twelve networks with
a fraction of claims reversed, three draws per rate: as the wrong fraction rises
from 0% to 51.6%, the `r = 1` share **falls** 57.8% → 42.1% and the median radius
**rises** 1 → 2, with the maxima for radius (14) and separation (7) unchanged.
**Wrong knowledge makes radii larger.** Centring on true knowledge was therefore
the choice least favourable to the finding, and the finding survives the
realistic case with room to spare.

Worth stating precisely why the framework itself is unaffected: `breakdown_radius`
takes only `(Ĉ, K, X, Y, Z)` and never consults the true DAG, so a wrong `G₀` is
just another element of the same space. Only the *centre* of the ball moves, and
only the sampling of centres was restricted.

Still uncovered: the corruption is uniform random reversal, whereas a real
analyst's errors are likely correlated. `tiered` is the vehicle for that and was
not used.
