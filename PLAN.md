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
- [x] T0.1 Resolve the radius off-by-one. Normative definition in a module
      docstring, asserted in a test.
- [x] T0.2 Freeze the shared core: validity oracle, space/distance, instance
      schema, results contract. Committed before any subagent is spawned.
- [x] T0.3 Determinism regression test extended to the new core.

## Axis B first (before large Axis A runs)

Conjectures 1 and 2 are the highest-value items in the session: either could
collapse the search space, and Axis A's scale ladder depends on the outcome.
They are cheap to test exhaustively on small graphs, so they come first.

- [x] B1.1 Conjecture 1 (failure is upward-closed) — exhaustive test.
- [x] B1.2 Conjecture 2 (retraction-optimal witnesses) — exhaustive
      counterexample search over all small CPDAGs and knowledge states.
- [x] B2.3 Admissible lower bound `r >= L` from structural failure events.
- [x] B2.4 Guided upper bound `r <= U`; report frequency of `L == U`.
- [ ] B2.1 Relevance pruning / quotient, with exactness proof obligation.
- [ ] B2.2 Component decomposition.
- [ ] B2.5 Best-first over shells using L as the heuristic.
- [ ] B2.6-9 Symmetry reduction, incremental closure, declarative encoding,
      anytime degradation. Lower priority; take as time allows.
- [x] B3 Differential testing against BFS on every feasible instance.

## Axis A

- [x] A2.0 Minimal generator + knowledge simulation + instance runner, small
      scale, so Axis B has instances to race on.
- [x] A1 Pre-registration written and committed BEFORE the first large run.
- [x] A2 Generators: Erdos-Renyi, scale-free, block, and the designed
      adversarial families for H3 with an explicit coupling parameter.
- [x] A3 Knowledge simulation: omission, flip, compound, tiered.
- [ ] A4 Scale ladder, cost reported against the largest chordal component the
      knowledge intersects (not against n).

## Deliverables

- [x] `results/synth/preregistration.md` (unedited after first commit)
- [x] `results/synth/`, `results/search/` with manifests
- [x] `report_synth_and_search.md`
- [x] `NEXT.md`
- [x] `PLAN.md`, `LOG.md`

## Sequencing rationale

Axis B's two conjectures gate everything: if failure is upward-closed and
witnesses are retraction-optimal, the search space drops from ~3^m to ~2^k and
Axis A can reach materially larger graphs. So B1 runs before the large A runs,
against a small A pipeline built purely to supply instances.

## Revisions

1. **Axis B moved ahead of the large Axis A runs.** The two conjectures were
   cheap to test exhaustively and gate everything downstream. Conjecture 1 turned
   out to be a theorem, so the effort went to Conjecture 2 and the accelerated
   search.
2. **Frozen-core change: instance params serialise to one `params_json` column**
   rather than exploding into `param_<key>`. Exploding them makes the CSV header
   depend on which generator produced the row, so a heterogeneous grid drifts
   schema mid-run. Caught by the results contract's own strict guard on the first
   Axis A run; no completed work needed re-running.
3. **Conjecture-1 checking switched from exhaustive to sampled.** It is
   O(|space|^2) per call and a theorem, so running it on every combination was
   guarding the implementation at the cost of dominating the sweep.
4. **H3's designed adversarial families abandoned as unachievable.** Identification
   requires compelled edges around the target; perturbability requires undirected
   ones. No construction satisfying both was found. Superseded by an exhaustive
   frontier sweep, which is stronger evidence anyway.
5. **B2.6-B2.8 (symmetry reduction, incremental closure, declarative encoding)
   not implemented.** The L/U certificate looked more valuable per unit of effort
   and was done instead. Recorded in NEXT.md.
6. **H4 not run.** Effort went to correcting the H2 statistic instead. Recorded
   as the cheapest remaining item.
