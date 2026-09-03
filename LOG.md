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
