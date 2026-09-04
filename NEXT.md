# NEXT — status after session 2 (Axis B deep)

Supersedes the session-1 version. Session-1 items that are now done are marked.

## What is now safe to put in a paper

1. **Failure is upward-closed** — a theorem, one-line proof. *(unchanged)*
2. **The space of knowledge states is a join-semilattice**, with the explicit join
   `G ∨ H = Meek(Ĉ, K_G ∩ K_H)`. **Proved** this session; verified on 260,716
   pairs. This is new and quotable on its own.
3. **Conjecture 2 (retraction-optimal witnesses) is proved**, modulo a single
   local property — Anti-Exchange Case B — whose Case A is proved and whose Case
   B has zero violations in 5,254 applicable triples. State it as: *proved
   conditional on Anti-Exchange*, not as an empirical regularity. This is a large
   upgrade on session 1's "~2M comparisons, no counterexample".
4. **`radius_local_up`'s exactness is exactly Anti-Exchange Case B**, and any
   failure is **one-sided**: radii come out too large, i.e. optimistic about
   robustness. Say this wherever exactness is claimed.
5. **The space-free method's cost tracks the up-set, not the space** — measured
   1.59^k against the up-set's 1.57^k, versus 1.96^k for the space and 3.97^k for
   BFS construction. Confirms the predicted `2^(#knowledge-oriented edges)`
   behaviour.
6. **The crossover rule**: building the space pays off only above roughly 175
   queries per CPDAG at k = 7, and ~50,000 at k = 10. For one-off certification,
   never build the space.
7. **The combined lower bound is defined in 100%** of finite-radius instances at
   n = 4 (was 22.8%), with zero admissibility violations. The certification rate
   is 85.4%, but quote it with the denominator caveat in §4 of the report.
8. **The breakdown radius does not saturate at 1** *(session 1, unchanged)*.
9. **No efficiency–robustness tradeoff** *(session 1, unchanged)*.

**Not yet safe.**

- **Anti-Exchange Case B as proved.** It is not.
- **Conjecture 2 unconditionally.** It rests on the above.
- **Full n = 5 completeness.** The corrected sweep is complete for k ≤ 6
  (all 8,782 CPDAGs, 1,977,820 comparisons, 0 counterexamples), but the 136
  densest CPDAGs — k > 6, or a corrected space above 200 elements — are still
  unexamined, and that is precisely where a counterexample would live.
- **Certification rates beyond n = 4.**
- **Anything at realistic graph sizes.** The ladder reaches k = 10, n = 7.
- **The naive `K`-count baseline (H4).** Still not run, two sessions on.

## The experiments that should come next

**1. Prove Anti-Exchange Case B.** Everything now hangs on one statement with a
clean reading: *no two distinct orientations are perfectly correlated across the
represented DAGs*. The obstruction is identified — the natural argument uses
chordality of chain components to reorient freely, and under background knowledge
a component need not be chordal (Task 0). Two routes worth trying: (a) prove it
for the cross-component case first, where orientations are nearly independent,
and characterise the same-component case separately; (b) attack it through the
convex-geometry structure, which is now established empirically (anti-exchange
holds, covers add exactly one element, the poset is graded).

**2. Close the n = 5 density gap.** The k ≤ 6 sweep is now complete (1,977,820
comparisons, 0 counterexamples); the 136 densest CPDAGs remain. A closure run was
under way at session end and had cleared 9 of them at k = 7…10 with 0
counterexamples. Finishing it is now much cheaper than it looks, because the
space-free method is
32,667× faster than BFS-with-construction at k = 10 and could do the sweep
without building spaces at all.

**3. Push the scale ladder using the space-free method.** The scaling data says
BFS construction is the binding cost (3.97^k) and `local_up` is not (1.59^k). Any
result that only needs radii — not the whole space — can now go materially
further than n = 7. The obstacle is that verification currently requires BFS as
ground truth; with Conjecture 2 proved conditional on Anti-Exchange, that
dependency is weaker than it was.

**Cheap and still unrun:** H4, the naive `K`-count baseline. Two sessions of
being the last unrun pre-registered hypothesis.

## Threads not pursued

- §3.6 fallback items (component decomposition, SAT/ILP/ASP encoding, incremental
  Meek closure, symmetry reduction) were not reached — §§3–5 did not stall.
  The SAT encoding remains the largest unexplored direction.
- `improved_upper_bound` exists but had nothing to close at n = 4, because the
  guided construction is already exact there. Its value, if any, is at larger k.
- Axis A remains paused.
