# NEXT — status after session 2 (Axis B deep)

Supersedes the session-1 version. Session-1 items that are now done are marked.

## What is now safe to put in a paper

1. **Failure is upward-closed** — a theorem, one-line proof. *(unchanged)*
2. **The space of knowledge states is a join-semilattice**, with the explicit join
   `G ∨ H = Meek(Ĉ, K_G ∩ K_H)`. **Proved** this session; verified on 260,716
   pairs. This is new and quotable on its own.
3. **Conjecture 2 (retraction-optimal witnesses) is proved**, modulo a single
   local property — Anti-Exchange Case B — whose Case A is proved and whose Case
   B has zero violations in 36,094 applicable triples, including the complete K₅
   skeleton. State it as: *proved
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
- **n = 5 is now complete**, with one qualification worth keeping: 135 of the
  136 densest CPDAGs were closed by direct sweep, and the last (the complete K₅
  skeleton, 4,231 space elements) was closed *via the proof chain* by verifying
  Anti-Exchange there, not by enumerating Conjecture 2. So that one inherits the
  Anti-Exchange dependency.
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

**2. Push the scale ladder to n = 6.** The n = 5 density gap is closed, so the
next frontier is n = 6, where exhaustive MEC enumeration is out of reach and a
stratified design is needed. This is now much cheaper than it looks, because the
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

---

# Addendum — status after session 3 (Axis B tractability)

Session-2 items are unchanged unless noted. Nothing in this session touched
Anti-Exchange Case B, the L/U bounds, or Axis A.

## What is now safe to put in a paper, added this session

1. **The breakdown radius is computable without constructing the space.** Three
   CP-SAT encodings, agreeing with brute-force BFS and `local_up` on every
   instance where either could run (1,044 exhaustive at n = 4; 346 of 360
   generated), zero disagreements.
2. **An UNSAT rung is a certificate that shells 0…k are clean**, produced without
   enumerating those shells. This is the form a robustness claim actually wants.
3. **The tractability picture is a split, not a speedup.** `local_up` wins
   decisively when a failure is near; the encoding wins — unboundedly — when it
   must be proved that no failure exists. The governing parameter is the
   **radius**, not `k`.
4. **Conjecture 2 survived an assumption-free test that could have refuted it**,
   over 1,007 certified witness walks including a deliberate hunt on dense,
   large-component, high-`|K_{G₀}|` instances.

## Carried forward, sharpened

- **Anti-Exchange Case B** is still the single open property, and E1's exactness
  still rests on it. This session gives it *no* new support: E3's soundness uses
  only the easy direction of Lemma R, so agreement between E1 and E3 bears on
  Conjecture 2 without touching Lemma R at all.
- **H4**, the naive `K`-count baseline, is still the last unrun pre-registered
  hypothesis — now for a third session.
- **Axis A** remains paused.

## New, and ordered by value

1. **Build the hybrid.** Run `local_up` under a small depth budget; if it finds a
   failure it answers in microseconds. If it exhausts the budget, hand the
   instance to E1. The two are strong in disjoint regimes and the discriminator —
   has a failure been found yet — is free at runtime. This is the practical
   deliverable and the data specifies it completely.
2. **Restrict the encoding to the knowledge-intersected component before emitting
   clauses.** Build cost dominates and is `O(n⁴)` in the vertex count, while the
   answer depends only on that component. Highest-value change to the encoding.
3. **Extend E3 past small radii** by incremental unrolling that reuses solver
   state across rungs. Right now E3 is capped at `r ≤ 3–4`, so the hunt's
   radius-5, -6 and -8 instances are untested for Conjecture 2.
4. **Explain the zero margins.** Every `r_E3 − r_E1` is exactly 0 across 1,007
   instances — never violated, but never with slack either. That is consistent
   with Conjecture 2 holding, and equally with the two encodings searching the
   same object for a structural reason not yet identified. Worth understanding
   *before* the tie is presented as evidence.
5. **Find E1's actual breaking point.** No E1 run hit its time limit here, so the
   reported ceiling is where I stopped, not where the method fails.

## Environment debt

23 tests fail and 3 files fail to collect on this machine, all `ImportError` on
`TypeAlias` / `StrEnum` in the `REPO_INIT` scaffold stubs, which need Python
3.11; the machine runs 3.9.6. Not caused by any session's work, but it means
`pytest tests` is not currently a clean signal and someone will eventually
mistake it for one.
