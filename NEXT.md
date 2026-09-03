# NEXT — what this session advances, what is safe to claim, what comes next

## What is now safe to put in a paper

**Safe, with the stated scope.**

1. **Failure is upward-closed — as a theorem, with a one-line proof.** Validity is
   a for-all over represented DAGs, so the DAG witnessing failure survives upward.
   State the empty-extension convention caveat alongside it. This licenses the
   pruning and memoisation the search relies on.
2. **The breakdown radius does not saturate at 1.** `r_val = 1` in ~80% of
   instances, but ~20% have radius ≥ 2, stable across three generators, across
   n = 6…9, and in an exhaustive census over every CPDAG on ≤ 5 nodes. This
   retires the project's principal risk.
3. **There is no efficiency–robustness tradeoff.** The frontier is genuinely
   non-flat (30.8% of instances), yet the optimal adjustment set is *never*
   strictly beaten (0 of 2,652). An analyst choosing `O*` gives up nothing in
   robustness. This is the most useful single claim the session produced, and it
   corrects the prior report's flat-frontier reading in both directions.
4. **An exact radius can be certified without enumerating the space in 83.3% of
   instances**, via matching lower and upper bounds — with the caveat that `L`
   covers only one of the two failure modes.
5. **The degeneracy rate is high and must be reported with any ensemble number.**
   91.8% of randomly generated instances are degenerate for this study, mostly
   because there is no confounding. Any rate quoted without this denominator is
   misleading.

**Not yet safe.**

- **Conjecture 2 as a theorem.** ~2M comparisons without a counterexample is
  strong evidence, not a proof, and the densest 2.2% of n=5 CPDAGs were skipped.
  Quote it as an empirically-supported conjecture with the enumeration scope.
- **Any claim at realistic graph sizes.** Everything here is n ≤ 9.
- **`r_ε` as an independent quantity.** It is nearly redundant with `r_val` unless
  ε is a substantial fraction of the effect (0.0% middle regime at ε = 0.02–0.05;
  14.1% at ε = 0.2).
- **Anything about the naive `K`-count baseline.** H4 was not run.

## The two or three experiments that should come next

**1. Prove or refute Conjecture 2.** Highest value by a distance. If the nearest
failure always lies in the up-set, the search space drops from ~`3^m` to ~`2^k`
*with a guarantee*, and the accelerated method stops being conditional. The proof
sketch to attack: for a failing `G`, does the join of `G0` and `G` (which fails,
by upward-closure) always lie at distance ≤ `d(G0, G)`? A counterexample would
need a failure incomparable to `G0` strictly nearer than any failure above it —
search there specifically, rather than uniformly. Failing a proof, close the 2.2%
scope gap at n=5 and push to n=6 with the space-free method.

**2. Complete the lower bound to cover the back-door failure mode.** `L` is
undefined in 7.6% of instances because the binding failure is a back-door path
becoming unblocked, not a descendant appearing. A bound for that mode would make
`L` a genuine certificate rather than half of one, and would likely push the
83.3% exact-certification rate substantially higher. This is the single change
that most improves the practical method.

**3. Resolve the MPDAG validity question (bug 6).** A Meek-closed graph
representing 5 DAGs, keeping every compelled edge and satisfying `[G] ⊆ [Ĉ]`, is
nonetheless excluded from the space, because chordality is checked on the
undirected subgraph alone and a component whose only chord is *directed* reads as
chordless. It affects 0.09% of Meek closures and is currently gated. Decide
whether the validity definition or the chordality test is wrong, then re-run
whatever depends on it. Do this before scaling up, not after: every result rests
on the current enumeration.

**Cheap and worth doing:** run H4. The runner needs one extra column and the
naive-radius machinery already exists in `bkrobust.demo.baseline`. It would
settle whether the prior report's off-by-one was systematic or a one-example
artefact.

## Work packages advanced

- *Structural theory*: Conjecture 1 promoted to theorem; Conjecture 2 supported at
  ~2M comparisons; the flat-frontier claim corrected to "non-flat, but `O*` never
  beaten"; the identification-vs-perturbability tension identified as a
  structural obstruction to a designed decoupled family.
- *Method*: an exact space-free search, differentially validated on 5,304 cases,
  with measured speedup and an honest statement of where it does not help; plus
  bounds certifying 83.3% of instances without enumeration.
- *Infrastructure*: a frozen core with a settled radius convention, an
  append-only results contract with self-describing manifests, and determinism
  guarantees enforced across `PYTHONHASHSEED`.

## Open threads not pursued

- Symmetry reduction via automorphisms, incremental Meek closure, and the
  SAT/ILP/ASP encoding (B2.6–B2.8) were specified but not implemented; the L/U
  certificate looked more promising per unit of effort and was done instead.
- Component decomposition (B2.2) was not implemented; the space-free search made
  it less urgent, but it is the natural route to larger graphs.
- The scale ladder above n = 9 was not attempted, because resolving item 3 above
  should come first.
