# Pre-registration — Axis A, synthetic graph ensembles

**Committed before the first large run. Not edited afterwards.** Deviations,
additions and anything learned that changes the analysis are recorded in
`report_synth_and_search.md`, not here.

Git SHA at writing: recorded in the commit that adds this file.

## Scope and standing assumptions

Fixed skeleton, causal sufficiency, faithfulness, oracle conditional-independence
testing. Finite-sample skeleton error is outside this analysis. Every claim is
about the ensembles named below and is not a claim about real data.

## Radius convention

`r = min distance from G0 to a graph where the property fails`, so shells
`0 .. r-1` are certified clean and the practitioner's safe-move count is `r-1`.
Fixed in `src/bkrobust/core/conventions.py` before any run.

## Statistics fixed in advance

- Radii are keyed off **mean** absolute bias, never the sampled maximum. The
  maximum is a sampled proxy, is unstable, and is reported as a diagnostic only.
- The governing cost parameter is the **size and density of the largest chordal
  component the knowledge intersects**, not `n`. It goes on the abscissa of
  every scaling plot.
- Every result row carries the method that produced it (`bfs_exact`,
  `local_up_exact`, or a labelled bound). No Axis A claim may rest on an
  unvalidated heuristic.
- `UNREACHED` is never averaged, plotted on a numeric axis, or compared as a
  number.

---

## H1 — Saturation

**Question.** Is `r_val = 1` the common case? The project's principal risk was
that the radius saturates at 1 and carries no information.

**Prediction.** `r_val = 1` is common but not dominant: I predict the mode is 1
and the distribution has non-trivial mass at 2 and above, with the mass at
`>= 2` increasing as the knowledge-intersected component grows.

**Statistic.** The **full distribution** of `r_val` per generator and per
knowledge corruption, reported as counts, not a mean.

**Falsified if.** More than 90% of non-degenerate instances have `r_val = 1`
across all generators. That would make the radius near-vacuous and would be a
serious negative result for the project.

## H2 — Frontier

**Question.** Does any valid adjustment set beat the optimal one on robustness,
i.e. is `max_Z r_val(Z) - r_val(O) > 0` ever positive?

**Prediction.** Yes, at a low but non-zero rate overall, and concentrated in the
decoupled families of H3.

**Statistic.** The distribution of `max_Z r_val(Z) - r_val(O)`, and the fraction
of instances where it exceeds 0.

**Falsified if.** The difference is identically 0 across every generator
including the designed decoupled families. That would say the flat frontier of
the prior report is a general phenomenon rather than an artefact of one graph.

## H3 — The flat-frontier mechanism (the key designed test)

**Background.** The prior report explained its flat frontier structurally: at
the witness, one DAG extension made the treatment an ancestor of the *entire*
confounder block at once, and since the empty set was confounded, every valid
adjustment set had to contain one of those covariates. So all sets failed
together.

**Design.** `decoupled_backdoor_dag(coupling)` builds two back-door routes whose
blocking sets are disjoint and, at `coupling = 0`, lie in **different chordal
components**. A perturbation confined to one component then cannot invalidate
the sets that block via the other, so the mechanism above cannot operate.

**Prediction.** At `coupling = 0` the frontier is **not** flat and `Z* != O` at a
materially higher rate than at `coupling = 1`. The rate should decrease
monotonically in `coupling`.

**Statistic.** Fraction of instances with `max_Z r_val(Z) > r_val(O)`, as a
function of `coupling`, swept over the full `[0, 1]` range.

**Falsified if.** The frontier stays flat at `coupling = 0`. **This outcome is
more important than a confirmation**: it would mean the prior report's
explanation is wrong or incomplete, and it must be reported as the headline
rather than buried.

## H4 — Naive baseline

**Question.** Is the naive `K`-count radius always at most the model radius? Is
the off-by-one seen in the prior report systematic or an artefact of one example?

**Prediction.** The naive count is `<=` the model radius in the large majority of
instances but **not always**; the strict off-by-one will not survive as a
universal rule.

**Statistic.** The joint distribution of `(naive, model)` radii; the fraction
with `naive > model` (which, if non-zero, refutes a universal inequality); and
the distribution of the inconsistent-single-flip fraction.

**Falsified if.** `naive > model` occurs in more than a negligible share, which
would refute the ordering entirely rather than merely the off-by-one.

## H5 — Regimes

**Question.** How often does the middle regime occur — `r_val` low but `r_eps`
high, i.e. fragile identification with a stable estimate? This determines
whether `r_eps` deserves to carry part of the contribution.

**Prediction.** The middle regime occurs in a substantial minority (I predict
20-50%) of non-degenerate instances.

**Statistic.** Joint distribution of `(r_val, r_eps)` at the pre-set epsilon
grid; fraction with `r_eps > r_val`.

**Falsified if.** `r_eps == r_val` almost always, which would make `r_eps`
redundant and the middle regime an artefact of the single prior example.

## H6 — Calibration

**Question.** How conservative is the certificate, and which graph descriptors
predict it?

**Prediction.** Coverage is exactly 1.00 everywhere **by construction** —
anything else is a bug and halts the run, not a finding. Conservativeness varies
materially across generators and increases with the size of the
knowledge-intersected component.

**Statistic.** Distribution of conservativeness; its rank correlation with each
graph descriptor.

**Falsified if.** Conservativeness is essentially constant across generators and
uncorrelated with every descriptor, meaning no descriptor predicts it.

---

## Degeneracy gate, fixed in advance

Instances are rejected, with the reason recorded and rejection counts reported,
when: the treatment is neither in nor adjacent to an undirected component; the
empty set is already a valid adjustment set; no valid adjustment set exists; or
no single atomic perturbation changes any validity status. **The rejection rate
is itself a reported result**, not a filtering detail.

## Analysis rules fixed in advance

- No generator, parameter, or epsilon will be changed after seeing results in
  order to obtain a nicer outcome. If a design changes after results are seen,
  the change and its reason are stated explicitly and the affected work re-run.
- Negative and null results are reported with the same prominence as positive
  ones.
- Every figure traces to a committed results file.
