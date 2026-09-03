# PLAN — synthetic-graph generalisation and search heuristics

Living document. Re-prioritised as results arrive; changes are recorded at the
bottom with a reason, never silently rewritten.

Branch: `experiments/synth_graphs_and_heuristics`. Never commit to `main`.

## Standing constraints

- No global RNG; every random draw descends from an explicit root seed.
- No dependence on set/dict iteration order — not in RNG consumption, not in
  float summation. Output bit-identical across runs and `PYTHONHASHSEED`.
- Cross-check every new primitive against an independent implementation.
- Calibration coverage < 1.00 is a **bug**: halt and investigate.
- Existing `results/breakdown_radius_demo/` is never modified or regenerated.

## Task 0 — conventions and frozen core (BLOCKING, orchestrator only)

- [x] T0.0 Branch gate.
- [ ] T0.1 Resolve the radius off-by-one. Normative definition in a module
      docstring, asserted in a test.
- [ ] T0.2 Freeze the shared core: validity oracle, space/distance, instance
      schema, results contract. Committed before any subagent is spawned.
- [ ] T0.3 Determinism regression test extended to the new core.

## Axis B first (before large Axis A runs)

Conjectures 1 and 2 are the highest-value items in the session: either could
collapse the search space, and Axis A's scale ladder depends on the outcome.
They are cheap to test exhaustively on small graphs, so they come first.

- [ ] B1.1 Conjecture 1 (failure is upward-closed) — exhaustive test.
- [ ] B1.2 Conjecture 2 (retraction-optimal witnesses) — exhaustive
      counterexample search over all small CPDAGs and knowledge states.
- [ ] B2.3 Admissible lower bound `r >= L` from structural failure events.
- [ ] B2.4 Guided upper bound `r <= U`; report frequency of `L == U`.
- [ ] B2.1 Relevance pruning / quotient, with exactness proof obligation.
- [ ] B2.2 Component decomposition.
- [ ] B2.5 Best-first over shells using L as the heuristic.
- [ ] B2.6-9 Symmetry reduction, incremental closure, declarative encoding,
      anytime degradation. Lower priority; take as time allows.
- [ ] B3 Differential testing against BFS on every feasible instance.

## Axis A

- [ ] A2.0 Minimal generator + knowledge simulation + instance runner, small
      scale, so Axis B has instances to race on.
- [ ] A1 Pre-registration written and committed BEFORE the first large run.
- [ ] A2 Generators: Erdos-Renyi, scale-free, block, and the designed
      adversarial families for H3 with an explicit coupling parameter.
- [ ] A3 Knowledge simulation: omission, flip, compound, tiered.
- [ ] A4 Scale ladder, cost reported against the largest chordal component the
      knowledge intersects (not against n).

## Deliverables

- [ ] `results/synth/preregistration.md` (unedited after first commit)
- [ ] `results/synth/`, `results/search/` with manifests
- [ ] `report_synth_and_search.md`
- [ ] `NEXT.md`
- [x] `PLAN.md`, `LOG.md`

## Sequencing rationale

Axis B's two conjectures gate everything: if failure is upward-closed and
witnesses are retraction-optimal, the search space drops from ~3^m to ~2^k and
Axis A can reach materially larger graphs. So B1 runs before the large A runs,
against a small A pipeline built purely to supply instances.

## Revisions

(append here, with date and reason)
