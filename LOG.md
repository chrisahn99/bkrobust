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
