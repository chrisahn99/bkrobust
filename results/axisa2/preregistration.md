# Pre-registration — session 5 (Axis A revisited: definitions and components)

**Written before the first large run and not edited afterwards.** Anything
learned later goes in the report, marked as post-hoc.

Git SHA at writing: recorded in `results/axisa2/manifest.json`.

---

## Background and the question

Session 1 measured `r_val = 1` in **81.0%** of the exhaustive n = 5 census
(88,700 finite radii; r=2 15.9%, r=3 3.1%) and **78.6%** across ensembles at
n = 6…9. A method that returns 1 almost always is not worth publishing, so the
question is whether that saturation is real or an artefact of two things:

1. **Definitional.** The repository's oracle is Pearl's back-door criterion,
   which is sufficient but not necessary for adjustment. The GAC is sound *and*
   complete and accepts strictly more sets, so fewer graphs in the ball are
   failures and the nearest failure can only move further away.
2. **Structural.** Undirected components were tiny — session 4's largest at
   n = 24 was 6 vertices, with 159 of 256 instances at 2–4. A component of
   three vertices cannot exhibit a radius of 3 whatever `n` is.

---

## Hypotheses

### H7 — separation (primary)

**Statement.** `r_val` increases with the **separation** `s`: the graph distance,
within the undirected component, from `X` to the nearest member of `O(G₀)`.

**Mechanism.** Failure requires the treatment to become an ancestor of a member
of `Z`, or a back-door/non-causal path to open. Separation is what stands in the
way: each intervening undirected edge is one more orientation that must be
retracted before the failure becomes reachable.

**Prediction.** Strong and monotone. Median `r_val` rises with `s`; the fraction
at `r_val = 1` falls monotonically in `s`.

**Falsified if** the radius distribution is flat in `s` — specifically, if the
fraction at `r_val = 1` differs by less than 10 percentage points between the
lowest and highest achieved separation stratum, at `n ≥ 30` instances per
stratum.

**Pre-committed quantitative guess.** At `s ≥ 3`, fewer than 50% of instances at
`r_val = 1`.

### H8 — component size (primary, with a genuine counter-mechanism)

**Statement.** Mass at `r_val ≥ 2` increases with component size `c`.

**The counter-mechanism, stated in advance.** A larger component gives more
separation (raising the radius) but also admits **more DAG extensions per
release** — so any single retraction has more chances to produce an extension in
which `Z` fails, which *lowers* the radius. These pull in opposite directions and
the net sign is genuinely open.

**Prediction.** Weakly positive net, dominated by separation.

**Analysis commitment, made before seeing data.** `c` and `s` are correlated by
construction (large `s` requires large `c`). They will therefore be reported
**jointly, as a two-way table of the radius distribution over `(c, s)`**, not as
two marginals. The claim "separation dominates" will be judged by whether the
radius still rises with `s` **within** a fixed `c` stratum, and whether it rises
with `c` **within** a fixed `s` stratum. If the cells are too thin to support
that, the conclusion is "underpowered", not a marginal comparison.

**Falsified if** at fixed `s`, larger `c` gives *lower* mass at `r_val ≥ 2` — in
which case the extension-count mechanism dominates, which is a publishable
finding in its own right and must be reported as such.

### H9 — definitional (direction forced, magnitude open)

**Statement.** GAC radii exceed back-door radii, and the fraction at
`r_val = 1` falls.

**The direction is forced**, not predicted: GAC accepts a superset of sets, so
`r_val(GAC) ≥ r_val(back-door)` on every instance. This is asserted as an
**invariant on every instance in every sweep**; a violation is a bug, halts the
run, and is not a finding.

**The open question is magnitude.** Pre-committed guess: **fewer than 15%** of
back-door `r_val = 1` instances get a strictly larger GAC radius, because the
optimal adjustment set `O*` should never contain a forbidden node, so the two
criteria ought to coincide exactly where most prior results were measured. If
that guess is badly wrong in either direction it is informative.

**Headline statistic of the definitional phase:** of the instances with
`r_val = 1` under back-door, the fraction with a strictly larger GAC radius.

### H10 — the frontier at realistic sizes

**Statement.** With larger components and GAC-admissible candidates, `O*` is
strictly beaten on robustness by some other valid adjustment set.

**Why it must be re-run rather than carried over.** Session 1's "`O*` is never
strictly beaten, 0 of 249,732" is a statement about **back-door-valid candidates
only**. GAC admits sets that were never in that pool.

**Prediction.** `O*` is beaten on a non-zero but small fraction — under 5%.

**Falsified if** it is never beaten across the new sweep, which would strengthen
session 1's negative rather than overturn it.

### H5-revisited — the `r_ε` middle regime (secondary)

Session 1 found 0.0% at small ε and 14.1% at ε = 0.2, against a pre-registered
20–50%. **Prediction:** the middle regime appears at small ε once radii are
larger, because it needs room between "valid" and "badly biased" to exist.

---

## Analysis commitments

- **Stratified reporting, never a single headline percentage.** Session 1's mean
  over an easy-dominated sample was misleading; the radius distribution is
  reported by `(c, s)` cell with per-cell counts visible.
- **Both definitions computed on every instance**, in one sweep, so the H9
  invariant is checked per instance rather than in aggregate.
- **Realised, not intended, parameters** are the ones analysed. Intended values
  are recorded alongside for the acceptance-rate accounting.
- **Sentinels.** `UNREACHED` is not a number: never averaged, never plotted
  numerically. Budget exhaustion is a distinct outcome from "no failure exists"
  and gets its own key.
- **Runtime against radius is a result**, recorded in its own right. If H7/H8
  hold, radii rise, the ladder goes deeper, and the expensive UNSAT-side work
  dominates. Cost growing is therefore *weak evidence the hypothesis is right*,
  and must not be presented as a mere bookkeeping note.
- **Assumption status on every number.** Both search legs are upward searches, so
  radii are exact **under Conjecture 2**, which rests on Anti-Exchange Case B —
  verified, not proved. Its error direction makes radii **too large**, which is
  **favourable to the hypothesis being tested here**. Any refutation of
  saturation therefore rests on an unproven assumption pointing in the same
  direction as the desired conclusion, and that caveat must be visible wherever
  the refutation is claimed.

## What would make this session a clean negative

If, under GAC definitions, on components of size 6–12 with separation swept
across its achievable range, the radius distribution still puts ~80% of mass at
`r_val = 1`, then saturation is real and structural. That goes in the first
paragraph of the report, undressed, and redirects the project toward the
`r_ε`/diagnostic framing rather than invalidating the work.
