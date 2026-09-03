# LOG — audit trail

Append-only. What was launched, what came back, what was integrated, what was
discarded and why.

## Session start

- Branch gate: `git branch --show-current` -> `experiments/synth_graphs_and_heuristics`. Passed.
- Read existing work: `src/bkrobust/demo/` (14 modules), `tests/demo/` (7 files,
  133 tests), `report.md`, `results/breakdown_radius_demo/`.
- Starting state: 4d0de90, clean tree, no `synth/` or `search/` yet.
