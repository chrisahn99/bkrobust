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
