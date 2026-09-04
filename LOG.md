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
