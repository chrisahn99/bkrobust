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
- [x] B2.1 Relevance pruning / quotient, with exactness proof obligation.
- [x] B2.2 Component decomposition.
- [x] B2.5 Best-first over shells using L as the heuristic.
- [x] B2.6-9 Symmetry reduction, incremental closure, declarative encoding,
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

---

# SESSION 2 — Axis B: proof, scale, and bounds

Branch unchanged: `experiments/synth_graphs_and_heuristics`.

## Task 0 (BLOCKING) — space membership
- [x] Determine whether the validity definition or the chordality test is wrong.
- [x] Measure how the exclusion rate behaves as n and density grow.
- [x] Re-run everything that depends on the enumeration; report old vs new.

## Conjecture 2
- [x] C2.0 Verify the semimodularity refutation computationally (record the dead end).
- [x] C2.1 Do joins exist? Is the space a join-semilattice?
- [x] C2.2 Lemma L: d(G0, G v G0) <= d(G0, G). Exhaustive, with slack distribution.
- [x] C2.3 Attempt a proof of Lemma L; identify the exact local property needed.
- [ ] C2.4 Close the n=5 dense gap (the 136 skipped CPDAGs).
- [ ] C2.5 n=6 stratified search, steered by non-gradedness.

## Speedup at scale
- [x] S1 Push k; measure up-set size vs space size; empirical scaling exponents.
- [x] S2 Machine-independent counters + peak memory.
- [x] S3 Crossover query count, as a practitioner rule.

## Bounds
- [x] B1 Back-door lower bound; verify admissibility exhaustively.
- [x] B2 Tighten U by local search.
- [ ] B3 Characterise when L = U.
- [x] B4 Scale the certification rate to n=5 exhaustive, n=6 sampled.

## Session 2 revisions
(append with reasons)

## Session 2 revisions

1. **Semimodularity was NOT avoided, contrary to the brief.** The brief's
   refutation conflates one knowledge assertion with one covering step; the poset
   is graded on all 203 spaces tested. Pursuing it produced the proof chain, so
   this was the single most consequential deviation of the session.
2. **The old `enumerate_space` was left unmodified** rather than fixed in place,
   so prior committed results stay reproducible from the code that produced them.
   Corrected logic lives in `search/space_fixed.py`.
3. **C2.4 (n=5 density gap) and C2.5 (n=6 stratified search) not completed.** The
   corrected n=5 sweep reached 62.4% before its agent stopped, and the 80 densest
   CPDAGs remain unexamined. Recorded as open in NEXT.md rather than glossed.
4. **B3 (characterise when L = U) not done**; effort went to the proof chain,
   which turned out to be the higher-value target.
5. **§6 fallback items not triggered** — §§3-5 did not stall.

---

# Session 3 — tractability (declarative encoding)

**Scope taken from the brief.** In: the SAT/CP encoding (item B2.8), scaling, and
stress-testing Conjecture 2 beyond brute force. Explicitly parked and not
touched: the Anti-Exchange Case B proof, the L/U bounds work, and Axis A.

## What was executed

1. **Three encodings, not one** (`src/bkrobust/sat/`). E1 retraction-only, E2
   join-distance, E3 assumption-free path unrolling. Kept separate because they
   assume different things and the differences *are* the experiment.
2. **Layered validation before any measurement.** Closure encoding vs the
   corrected space; failure predicate vs the validity oracle; radii vs BFS and
   `local_up`; then cross-encoding. Each layer found a bug before the next ran.
3. **Certificates over solver claims** (`sat/verify.py`). Every E3 witness walk
   is replayed through the ordinary graph code and re-checked without consulting
   the encoding.
4. **Scaling study** to n = 20 across six generators, and a **targeted high-`k`
   hunt** on dense CPDAGs reaching `|K_{G₀}| = 51`.
5. **Report numbers are derived, not transcribed.** The PDF computes every figure
   it quotes at build time, and `analysis/session3_verify.py` independently
   asserts the markdown agrees with the same files.

## Deviations from the plan, and why

1. **E2 needed no extra copy of the orientation variables.** The brief budgeted
   one for the join. Theorem 2 plus closure-under-intersection makes the join
   determined rather than searched for, collapsing the objective to a symmetric
   difference over one copy. Assumptions unchanged; only the cost changed.
2. **`local_up` had to be run under a hard wall-clock cap, in a subprocess.**
   The first high-`k` sweep hung, because `local_up` on a degenerate instance at
   k = 15 must exhaust an up-set of up to 2¹⁵ states. Rather than lower the
   ambition, the timeout was made a *recorded datum*. That change is what
   produced the session's central result (§3.3 of the report); the hang was the
   measurement.
3. **E3 was budgeted by radius, not by instance size.** It unrolls `r+1` copies
   of the closure encoding, so it ran only where `r ≤ 4` (general) or `r ≤ 3`
   (hunt). Instances at radius 5, 6 and 8 in the hunt are therefore **not
   covered** by the Conjecture 2 test. Recorded as a scope limit in the report,
   not glossed.
4. **The §3.6 fallback was not triggered as a fallback, but one of its items is
   now the top recommendation.** The encoding was not a dead end, so no switch
   was recorded. Component decomposition nonetheless became the highest-value
   next change, because the measured profile points at it: build cost dominates
   and is `O(n⁴)` in the vertex count while the answer depends only on the
   knowledge-intersected component.
5. **No crossover against `local_up` was found on the SAT side, and none is
   expected.** The brief asked for the crossover; the honest answer is that it
   exists only on the UNSAT side. Reported as the headline rather than buried.
