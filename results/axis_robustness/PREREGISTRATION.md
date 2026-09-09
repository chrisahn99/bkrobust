# Pre-registration — Axis Robustness (survival curves + Pareto frontier)

**Branch:** `experiments/survival-and-pareto`
**Written:** 2026-09-09, after a read-only API scoping pass, before any survival curve,
AUC, or rank correlation was computed.
**Status at time of writing:** no experiment code written, no data generated.

This file follows the project convention of prediction-before-looking. It exists so the
Phase 1 rank correlations cannot be read as post-hoc. Nothing below may be edited after
the first result file lands in `results/axis_robustness/`; corrections go in a dated
appendix at the foot of this file.

---

## 1. What is being measured

An analyst holds a CPDAG `Ĉ`, asserts background knowledge `K`, Meek-closes to `G₀`, and
reads an adjustment set off `G₀`. We corrupt `K` and ask whether the set they already
committed to is still valid.

**Operational definitions, fixed in advance:**

- `Z* := optimal_adjustment_set_mpdag(G₀, x, y)` (`demo/evaluate.py:288`) — computed once,
  at depth 0, then **held fixed**. Survival concerns the analyst's committed set, not
  whether *some* valid set still exists. Recomputing `Z` after corruption would measure a
  different and much weaker claim. Instances where this returns `None` (extensions
  disagree ⇒ non-identification) are excluded at generation with status
  `optimal_set_undefined`, not scored.
- **Corruption depth `d_claims`** = number of asserted orientations that differ from the
  truthful knowledge set.
- **Survival** at a sampled corrupted state `G`: `is_gac_valid_mpdag(G, x, y, Z*)` is True
  (`gac/mpdag_level.py:298`, polynomial, does not enumerate `[G]`).
- `S(d)` = fraction of sampled states at depth `d` that survive, `N = 200` samples.

### 1.1 The two corruption processes are not symmetric, and are not treated as such

Scoping established that `synth/knowledge.py::flip` and `::tiered` have **different
signatures and different semantics**:

- `flip(k, rng, rate)` — corrupts an *existing* knowledge list; reverses
  `round(rate·|k|)` entries drawn without replacement. Independent, per-edge.
- `tiered(dag, cpdag, rng, n_tiers, corruption_rate=0.0)` — **generative, not corruptive.**
  It does not take a `k`. It ranks nodes by longest-path depth, splits them into
  `n_tiers` contiguous blocks, and asserts every cross-tier CPDAG edge earlier→later.
  `corruption_rate` relocates that fraction of *nodes* to a random different tier, so a
  single misplaced node can flip several assertions at once.

So `tiered` cannot be driven to a target `d_claims` directly. Protocol: generate
`K_ref := tiered(dag, cpdag, rng, n_tiers, 0.0)` and
`K_cor := tiered(dag, cpdag, rng, n_tiers, rate)`; **measure** `d_claims` post hoc as the
number of assertions in `K_cor` that contradict or are absent from `K_ref`. The tiered arm
is therefore a *measured-depth* design, the flip arm a *targeted-depth* design. They are
analysed as separate strata and never merged into one curve.

---

## 2. `d_hops` — what is measurable and what is not

`r_val` is defined in BFS hops on the covering-relation neighbour graph; the corruption
process is defined in claims. These are different units, and the mapping is not identity.

Scoping established a hard constraint: **no cheap single-pair distance function exists.**
`search/exact_fast.py:91 radius_local_up_fast` searches *upward* only, and a flip-corrupted
`G` is generally **not above** `G₀` in the refinement order — neither is a refinement of
the other — so the upward search can never reach it. `radius_local_up_fast` is therefore
**not** a distance oracle for this design, and will not be used as one. (This is also why
the radius error is one-sided: an upward-only search can miss a nearer non-upward failure,
which inflates `r`, never deflates it.)

Consequences, fixed in advance:

- **Primary depth axis is `d_claims`.** Every headline curve is in claim units.
- **`d_hops` is computed exactly**, via `core/spacelib.py:53 build_space` +
  `:74 distances_from`, **only on a small-instance subsample** where full space
  enumeration is tractable (`component_size ≤ 8`, budget-checked before the sweep).
  P5 is tested on that subsample and its size is reported.
- **Symmetric difference `|dir(G₀) Δ dir(G)|` is recorded on every state as a cheap
  proxy, and is labelled a proxy everywhere it appears.** `THEOREMS.md` §10–§11 give only
  `r ≤ r_E2 ≤ r_E1` — the symmetric-difference objective is a documented **upper-bound
  surrogate**, and the converse (Lemma R) is *verified, not proved*. It is never presented
  as a distance.

---

## 3. Sentinel discipline (fixed in advance)

- Corrupted `K` admitting no consistent MPDAG (`demo/meek.py:207 apply_orientations`
  returns `None`) is status `corrupted_k_contradictory`. It is **not** a survival-0
  outcome. `S(d)` is computed over non-contradictory samples only; the contradiction rate
  is reported per depth as its own series.
- Pre-specified sensitivity: recompute everything with contradictions scored as failures
  (`S_contra_as_fail`). Both appear in the report. If they disagree, that is reported
  prominently rather than resolved by choosing one.
- `UNREACHED = -1` is never averaged, plotted numerically, or fed to a correlation.
  Instances with `r_val == UNREACHED` are excluded from τ and counted out loud.
- `DEGENERATE = 0` instances are gated out at generation, not filtered post hoc.

---

## 4. The AUC confound, and how it is neutralized

Unnormalized area `Σ_{d=1..|K|} S(d)` grows mechanically with `|K|`, which would
manufacture a correlation with one of the very baselines under test.

- **Primary endpoint: `AUC_frac`** — `S` on the fixed fractional grid
  `d/|K| ∈ {0.1,…,1.0}` (nearest achievable integer `d`), then averaged. Dimensionless,
  comparable across instances with different `|K|`.
- **Secondary: `AUC_abs`** — the raw sum, reported only next to `|K|` so the confound is
  visible rather than hidden.

Primary endpoint is `AUC_frac`. Fixed now.

---

## 5. Predictors under test

| Predictor | Definition | Unit |
|---|---|---|
| `r_val` | `hybrid.py:102 breakdown_radius(...).radius` | hops |
| `shd_truth` | SHD(`G₀`, true DAG) | edges |
| `shd_cpdag` | SHD(`Ĉ`, `G₀`) — orientations committed incl. Meek propagation | edges |
| `n_k` | `|K|`, the naive claim count (project hypothesis **H4**, unrun for six sessions) | claims |
| `undirected_fraction` | of `Ĉ`, instance-level | ratio |

**`shd_truth` is degenerate on true-knowledge instances.** At coverage 1.0 `G₀` *is* the
ground-truth DAG, so `shd_truth ≡ 0` and τ against it is undefined. The Phase 1 population
must therefore span corruption > 0 at the *instance* level. Any stratum in which a
predictor is constant is reported as **"undefined (predictor constant)"** — never as
τ = 0, and never as evidence that the baseline "fails".

### 5.1 Instance population (fixed in advance)

To make `shd_truth` non-degenerate, the analyst's *starting* knowledge must itself vary in
correctness. Instances are generated as
`component_generator.generate_instance(ComponentSpec(component_size, separation, coverage),
seed)`, then the truthful `K` is pre-corrupted at a **base wrongness rate**
`b ∈ {0.0, 0.10, 0.25}` via `flip`, and `G₀ := apply_orientations(Ĉ, K_b)`. The survival
sweep then corrupts *further*, from that `G₀`.

This mirrors the `results/axisa3/wrong_knowledge.jsonl` design and yields `shd_truth`
spread across instances without altering the endpoint. Note the endpoint remains
*structural survival* — `is_gac_valid_mpdag(G, x, y, Z*)`, validity of the committed set
across states in the space — not whether `Z*` is correct in the true DAG. Those are
different questions and only the first is what `r_val` is defined to measure.

Because `tiered` generates its own knowledge set, the tiered arm's instances have their own
`K_ref`, `G₀`, `Z*` and therefore their own `r_val`. The two arms are **different instance
populations** and are never compared within a stratum.

---

## 6. Predictions (the part that must not move)

- **P1.** τ_b(`AUC_frac`, `r_val`) > 0 in the majority of strata. Direction: positive.
- **P2.** τ_b(`AUC_frac`, `r_val`) exceeds τ_b(`AUC_frac`, `n_k`) in the majority of
  strata where both are defined.
- **P3.** τ_b(`AUC_frac`, `n_k`) is not distinguishable from 0 (bootstrap CI covers 0) in
  the majority of strata. This is H4 finally being run.
- **P4.** No directional prediction for `tiered` vs `flip` mean AUC. Directional
  prediction on spread: **`tiered` produces higher across-instance variance in `AUC_frac`
  than `flip`**, because whether a correlated block lands on the X–Y causal path is a coin
  flip that uniform corruption averages over.
- **P5.** τ on `d_hops` and on `d_claims` agree in sign, on the small-instance subsample
  of §2. If they do not, the headline claim is withdrawn pending investigation.

**Pre-committed falsification.** If τ(`AUC_frac`, `r_val`) is not positive in the majority
of strata, or if `n_k` matches or beats `r_val`, that is the result and it is reported as
the result. The word "perfectly" appears nowhere in this analysis plan: the endpoint is a
rank correlation with a bootstrap interval, not a claim of perfect ordering. A finding
that `r_val` merely *ranks better than* the baselines is the realistic success case and is
what will be claimed if observed.

---

## 7. Analysis plan

- **Kendall τ-b.** Ties are expected and material (`r_val` is small-integer valued).
  Never τ-a.
- **Stratified, never pooled.** Strata = (corruption process) × (coverage) × (family or
  network). No pooled τ as a headline. Admissible instances collapse as coverage falls and
  survivors have larger separation, so pooling across coverage measures composition, not
  the law.
- **Fixed instance set** within any comparison across coverage.
- **Bootstrap CI:** 10,000 resamples over *instances*, not over Monte-Carlo draws, seed 0.
- **Gate: `benchmarks.measure.fast_gate` exclusively**, recorded in every row.
  `synth.runner.gate` is **not** used: scoping confirmed it calls
  `all_valid_adjustment_sets_mpdag` (`demo/evaluate.py:197`, exponential) internally. Note
  `fast_gate` is strictly stricter (16/46,800 disagreements, one direction, biasing against
  the hypothesis). A number without its gate label is not quotable.
- Every radius row carries `HybridResult.assumes` verbatim: *"Conjecture 2 (hence
  Anti-Exchange Case B, verified not proved)"*. All upward searches are exact iff
  Conjecture 2 holds; the error is one-sided (radii can only be too large).
- `all_valid_adjustment_sets_mpdag` is never on a decision path.

## 8. Determinism

`num_workers: 1`, `random_seed: 0`, CP-SAT single-threaded. Scoping confirmed the live
tree contains **zero** uses of `random.seed`, `np.random.seed`, or `import random` — all
randomness flows through explicit `numpy.random.Generator` arguments. That property is
preserved: every sampler takes an explicit `Generator`, seeded by a documented derivation
from `(instance_id, process, depth, rep)`. Bit-identical output across two values of
`PYTHONHASHSEED` is **checked, not assumed**.

Result IO uses `core/resultsio.py`: `ResultWriter` (CSV, schema-drift detection, flush per
row — satisfies the incremental-write constraint) and `write_manifest` (captures git SHA,
dirty flag, interpreter, platform, library versions and `RADIUS_CONVENTION` automatically).
No jsonl writer exists in the live tree; CSV is used rather than hand-rolling one.

## 9. Phase 2 (Pareto) — registered in outline

Endpoint: for each admissible BK proposal on a fixed CPDAG, the pair (`r_val`, utility),
with utility from `demo/evaluate.py:503 asymptotic_variance` (n-free, Frisch–Waugh form;
caller multiplies by `1/n`) and, separately, `:498 bias` under a fixed 10% domain-error
injection. `LinearSEM` (`:330`) requires a **fully oriented DAG**, so SEMs are built on the
ground-truth DAG via `random_sem(dag, rng, ...)` (`:423`); the MPDAG state supplies only
the adjustment set, never the SEM.

- **P6.** The non-dominated frontier is non-trivial: there exist proposals strictly better
  on utility and strictly worse on `r_val`. The objectives genuinely trade off.
- **P7.** Aggressive BK (claims on edges lying on X→Y paths) concentrates in the
  high-utility / low-`r_val` region; conservative BK (peripheral edges) in the
  low-utility / high-`r_val` region.

If P6 fails — frontier is a single point, or the objectives are monotonically aligned —
there is no "fragility of precision" trade-off and the Phase 2 claim is dropped rather
than rescaled.

---

# Appendix A — 2026-09-09, after the first Phase 2 run

Per §0, the body above is frozen. This appendix records a falsification and registers a
revised question. Everything below is **post-hoc motivated and labelled as such**; it does
not inherit the pre-registered status of §6 or §9.

## A.1 P6 is falsified, and the reason is structural

Run: 24 base CPDAGs, 1,087 proposals, 17,680 SEM draws (`pareto_proposals.csv`,
`pareto_sems.csv`, `pareto_manifest.json`). Independently re-checked by the orchestrator:

- strata with more than one distinct `Z_p`: **0 of 24**
- strata with more than one distinct `utility_median`: **0 of 24**
- strata with more than one distinct `r_val`: **18 of 24**

**Proposition (why this had to happen).** Every proposal in this design asserts *truthful*
orientations, so the true DAG is always a member of `[G₀]`, the extension set of `G₀`.
`optimal_adjustment_set_mpdag` returns non-`None` precisely when all extensions agree on
the optimal set. If they agree, the common value is in particular the true DAG's optimal
set. Therefore **every identifying truthful proposal yields the same `Z`, hence the same
asymptotic variance.** Utility is invariant by construction; one axis of the intended
Pareto plane is a constant, so no trade-off can exist to detect.

This is not a null result to be worked around. It is a fact about the estimand:
**with truthful background knowledge there is no efficiency price for robustness, because
the optimal adjustment set is an invariant of the MPDAG whenever it is identified at all.**
The "fragility of precision" cannot live in the efficiency of the optimal set. If it exists
anywhere, it lives in **identifiability** — 203 of 1,087 proposals (18.7%) failed to
identify any optimal set.

**P6 as stated in §9 is dropped.** It is not restated in weakened form, and the secondary
bias-based front (both axes jointly varying in only 2 of 24 strata) is not offered as a
rescue.

## A.2 P7 is untestable on this generator, and why

**Zero** of 1,087 proposals classified as "aggressive". `synth/component_generator.py`
builds the perturbable undirected component entirely out of `X`'s upstream confounder
structure, holding `X→Y` and everything downstream of `X` outside that component by
construction. No undirected edge in any CPDAG this generator emits can lie on a directed
`X→Y` path, so the aggressive/conservative contrast has an empty cell. This is a property
of the instance family, not of the metric.

## A.3 Revised questions (post-hoc, newly registered, not pre-registered)

- **P6′ (achievable efficiency).** Utility is redefined as the efficiency of the best set
  the analyst can actually *justify from `G₀`*, over a fixed candidate menu of at most 12
  sets constructed once per base CPDAG and filtered by `is_gac_valid_mpdag`. Weaker
  knowledge should leave fewer menu members valid and force a less efficient fallback.
  Prediction: within a stratum, `r_val` and achievable utility are negatively associated
  in a majority of strata where both vary. **If achievable utility is also invariant, the
  Phase 2 claim is abandoned entirely rather than redefined a third time.**
- **P7′ (aggressive BK).** Tested on real benchmark CPDAGs from
  `results/axisa3/networks/`, where undirected edges *can* lie on `X→Y` paths, rather than
  on the synthetic generator. If no benchmark network admits such a pair either, the
  aggressive/conservative contrast is reported as structurally unavailable and dropped.

The menu is explicitly capped and is **not** `all_valid_adjustment_sets_mpdag`; the 2^V
enumerator remains off every decision path.

---

# Appendix B — 2026-09-09, after the Phase 2 round-2 run

## B.1 P6′ is falsified on the same mechanism. Phase 2 is abandoned.

Run: `pareto2_proposals.csv` (1,431 rows), `pareto2_menus.csv` (243),
`pareto2_aggressive_census.csv` (39), 97.6 s wall. Orchestrator re-verified:

- strata with >1 distinct `achievable_utility_median`: **0 of 32**
- strata with >1 distinct `n_valid_menu_members`: **0 of 32**
- strata with >1 distinct `winning_member_index`: **0 of 32**
- strata with >1 distinct `r_val`: **23 of 32**

Redefining utility as *achievable* efficiency over a bounded 12-set menu did not
introduce variation. Appendix A.3 pre-committed the response to exactly this outcome:
**Phase 2 is abandoned. No third definition of utility is attempted.**

## B.2 What was actually learned (this is a finding, not a consolation)

The two independent falsifications converge on one statement:

> **For truthful background knowledge, statistical efficiency is an invariant of
> `(Ĉ, x, y)` — not a function of what the analyst asserts. Knowledge moves the
> robustness axis and only the robustness axis.**

Round 1 showed the *optimal* set is shared by all identifying proposals. Round 2 showed
the stronger property: the entire GAC-valid subset of a fixed menu is identical across
proposals within a stratum, so even the fallback the analyst would be forced onto does not
move. Efficiency is constant while `r_val` varies in 23 of 32 strata.

The consequence for practice is the reverse of the "fragility of precision" hypothesis.
Precision is not traded against robustness. Precision is **fixed**, and marginal truthful
background knowledge purchases **only fragility**. What knowledge does buy is
**identifiability**: 385 of 1,431 proposals (26.9%) could not justify any valid set from
the menu at all, and 203 of 1,087 in round 1 identified no optimal set.

So the operative trade-off is **identifiability against robustness, and it is a step, not
a frontier**: assert enough to identify, and every further claim is pure `r_val` cost at
zero precision gain. This is a descriptive reading of data already collected, **not** a
third redefinition — no new utility axis is proposed and no new experiment is run for it
here. Testing it properly (does `r_val` fall monotonically with proposal size beyond the
minimum identifying set?) is a separate question for a later session and is not claimed now.

## B.3 P7′ — structurally available, but the contrast is untestable

The census is the useful part and stands on its own: sampling ≤200 descendant-ordered
pairs per network across the 39-network corpus (M-bias excluded as an ADMG),
**28 networks admit at least one (x, y) pair with an undirected CPDAG edge lying on a
directed x→y path**, 396 such pairs in total — densest in `paths` (96), `Polzer_2012` (64),
`Kampen_2014` (33), `child` (22), `sachs` (22), `insurance` (21). Round 1's zero was an
artifact of `component_generator`, not a property of real structure.

*(The round-2 worker's own summary said 23 of 39 networks; the committed CSV gives 28. The
file is authoritative; the 396-pair total agrees.)*

The follow-up gated 8 real-network bases (`mediator` ×3, `sachs` ×2, `Kampen_2014` ×3),
yielding 215 aggressive and 129 conservative proposals — so the empty cell is genuinely
filled. But utility is invariant there too (0 of 8 strata vary), so the intended
aggressive-vs-utility contrast cannot be run. The residual `r_val`-only contrast is
non-directional (2 bases lower for aggressive, 2 tied, 1 reversed) and is reported as
**inconclusive**, not as support in either direction.

**P7′ is not supported and not refuted; it is untestable given B.2.**

---

# Appendix C — 2026-09-09, after the Phase 1 sweep completed

## C.1 A measured fact that forces a change of primary endpoint

The sweep completed: 1,978 instances, 3,136,000 samples, 1,567 s. `n_s0_failures = 0`.
Then the contradiction rate turned out to be far higher than the design assumed.

Pooled contradiction rate by corruption depth (`survival_curves.csv`, orchestrator-verified):

| `d` | 1 | 2 | 3 | 4 | 5 | 6 | 8 | 10 | 12 |
|---|---|---|---|---|---|---|---|---|---|
| contradiction rate | 0.490 | 0.773 | 0.896 | 0.893 | 0.930 | 0.947 | 0.985 | 0.984 | 1.000 |

Split by arm and base wrongness at low depth, it is stable and not an artifact of the
pre-corruption:

| arm | base wrongness | d=1 | d=2 | d=3 |
|---|---|---|---|---|
| flip | 0.00 | 0.651 | 0.799 | 0.884 |
| flip | 0.10 | 0.645 | 0.795 | 0.881 |
| flip | 0.25 | 0.564 | 0.718 | 0.859 |
| tiered | — | 0.354 | 0.769 | 0.917 |

**Consequence.** `S(d)` was pre-registered (§3) over non-contradictory samples only. At
`d ≥ 3` that denominator is under 15% of draws and at `d ≥ 12` it is zero. Median depths
per instance = 9, but **median depths retaining ≥30 usable samples = 3**; 18.7% of
instances have fewer than 2 usable depths and 22 have none beyond `d = 0`. `AUC_frac` as
defined in §4 therefore averages a grid that is mostly Monte-Carlo noise.

**Change, and it is a deviation from §4 driven by measurement, not by preference:**
the primary endpoint becomes the AUC of **`S_contra_as_fail`**, whose denominator is fixed
at 200 for every depth. `S` (contradictions excluded) becomes the sensitivity analysis —
the two swap roles. §7's Task-0 instruction to the analysis worker pre-authorized exactly
this swap conditional on the measurement, so it is not a free choice made after seeing a
correlation; no τ had been computed when this was decided.

Semantically the swapped endpoint is also the better one, which is worth stating plainly
rather than treating as luck: if the analyst's corrupted beliefs admit no consistent MPDAG,
their committed adjustment set has failed them. "No graph at all" is a failure of the
analyst's position, not an undefined outcome.

## C.2 Monotonicity: mixed, mostly artifact

Of 1,557 instances with curves, 629 (40.4%) are non-monotone under `S` and 592 (38.0%)
under `S_contra_as_fail`; **381 — 60.6% of the non-monotone set — become monotone once
contradictions are scored as failures.** So the majority of the non-monotonicity is the
conditioning artifact anticipated in the run notes, but a residue survives.

**Caveat the analysis must respect:** the residue is measured with a zero-tolerance test
(any uptick at all counts). With 200 reps the standard error of `S` is ≈0.035, so noise
alone produces upticks. The residual non-monotonicity is **not** established as genuine
until it is retested against Monte-Carlo error. It is not claimed here.

## C.3 The finding this actually produced

> **Most orientation errors are self-revealing.** Reversing a single knowledge claim makes
> the knowledge set inconsistent with the CPDAG 56–65% of the time under uniform flips and
> 35% under correlated tiered corruption; by depth 3 it is ~88%. The majority of wrong
> background knowledge cannot be expressed as an MPDAG at all — Meek closure fails and the
> analyst finds out.

This reframes the robustness question rather than answering it differently. The dangerous
errors are the *minority that stay consistent*, because those are the ones an analyst
cannot detect from the graph alone. The space `r_val` is defined over consists precisely of
consistent MPDAGs, so `r_val` is measuring the residual risk that survives self-consistency
filtering — which is the risk that matters. This holds independently of how the rank
correlations come out.

It also explains why uniform `flip` and correlated `tiered` differ most at `d = 1`
(0.65 vs 0.35): block-correlated errors move a whole tier coherently and are far more
likely to stay internally consistent, hence **harder to detect**. That is a point in favour
of taking correlated analyst error seriously, and it was the open question left by
`report_real_graphs.md` §6.

## C.4 The swapped endpoint is not degenerate (checked before any τ was computed)

Scoring contradictions as failures raises an obvious objection: if `S_contra_as_fail` were
essentially `1 − contradiction_rate`, the endpoint would measure structural rigidity rather
than adjustment-set validity, and any correlation with `r_val` would be circular. Measured:

- `mean(S_contra_as_fail − (1 − contradiction_rate)) = −0.172`, sd 0.228; the two agree to
  within 0.02 in only 38.3% of depths and within 0.10 in 50.1%.
- Among **non-contradictory** samples at depths with ≥30 usable draws (213 depths),
  mean `S = 0.397`, median 0.487; only 4.7% of those depths have `S = 1.0` and 54.0% have
  `S ≤ 0.5`.

So validity failure contributes substantially beyond contradiction: among corrupted
knowledge sets that *remain consistent*, the committed adjustment set still fails about half
the time. The endpoint carries genuine GAC-validity signal. **Not degenerate.**

---

# Appendix D — 2026-09-09, correcting Appendix C

## D.0 Why C needed correcting

The sweep that produced Appendix C was destroyed mid-analysis when its worker relaunched
and overwrote its own outputs in place. It was regenerated (run 2, 1,978 instances,
14,066 curves, 3,136,000 samples, 1,656 s, `n_s0_failures = 0`, exclusions
`k_contradictory = 2869` / `optimal_set_undefined = 777` — identical to run 1).

Run 2 also carries a new `base_instance_id` column, which was the point of the relaunch.
Investigating the difference showed the relaunch was **not** cosmetic, contrary to the
orchestrator's first reading: in run 1, `instance_id` was not unique, and 1,978 instances
collapsed into ~1,553 groups when keyed by it. Every Appendix C statistic computed by
*grouping on `instance_id`* was therefore computed over merged instances. Statistics
computed by *pooling rows* were unaffected.

Run 2's `instance_id` is unique (1,978 distinct of 1,978). Run 2 is authoritative.
All numbers below supersede Appendix C.

## D.1 CONFIRMED — the contradiction rates and the C.3 finding

Pooled rates are reproduced to the third decimal: `d=1` 0.491, `d=2` 0.773, `d=3` 0.896,
`d=12` 1.000. By arm and base wrongness at `d=1`: flip 0.653 / 0.647 / 0.563 at
`b = 0.0 / 0.10 / 0.25`, tiered 0.354. Usable-depth structure holds: median 7 depths per
instance, **median 3 with ≥30 usable samples**, 20.7% of instances with fewer than 2.

**C.3 stands unchanged.** Most orientation errors are self-revealing, and correlated
(tiered) errors are markedly harder to detect at `d=1` (0.354 vs 0.653) than uniform ones.

## D.2 CORRECTED — monotonicity

| | Appendix C (run 1, collided) | run 2 (authoritative) |
|---|---|---|
| non-monotone under `S` | 40.4% | **28.3%** |
| non-monotone under `S_contra_as_fail` | 38.0% | **31.4%** |
| of the non-monotone set, fixed by scoring contradictions | 60.6% | **78.9%** |

The ordering reverses: `S_contra_as_fail` is *more* often non-monotone than `S`, not less.
Both are true simultaneously — 441 of 559 instances non-monotone under `S` become monotone
under `S_contra_as_fail`, while 503 instances monotone under `S` become non-monotone under
it. The cause is visible in D.1: **the contradiction rate is itself not monotone in `d`**
(0.986 at `d=8`, 0.937 at `d=9`, 0.984 at `d=10`, 0.887 at `d=11`), and
`S_contra_as_fail` inherits that.

Under a 2-SE noise tolerance (0.0707 at 200 reps) the figures fall to 14.7% (`S`) and
12.0% (`S_contra_as_fail`), so roughly half of the raw non-monotonicity is Monte-Carlo
noise and roughly half is not. Neither series is cleanly monotone.

## D.3 RETRACTED — C.4's "not degenerate" verdict, and the endpoint swap it justified

C.4 concluded that `S_contra_as_fail` carries validity signal beyond contradiction
propensity. On run 2 that is **false**:

| | Appendix C (run 1) | run 2 (authoritative) |
|---|---|---|
| `mean(S_contra_as_fail − (1 − contradiction rate))` | −0.172 | **−0.056** |
| agreement within 0.02 | 38.3% of depths | **78.1% of depths** |
| conditional `S` among non-contradictory, mean / median | 0.397 / 0.487 | **0.717 / 0.952** |
| depths with conditional `S ≤ 0.5` | 54.0% | **23.6%** |

`S_contra_as_fail` coincides with `1 − contradiction_rate` to within 0.02 in 78% of depths.
It is largely a re-expression of detectability, not of adjustment-set validity.

**Therefore C.1's endpoint swap is withdrawn and §4's pre-registered primary is restored:
`AUC_frac` over `S` (contradictions excluded).** The pre-registration was right the first
time — contradiction is not a validity failure — and the swap was made on collided data.
No deviation from §4 stands.

To keep both mechanisms visible rather than choosing between them, the analysis reports the
**exact decomposition**

> `S_contra_as_fail(d) = (1 − contradiction_rate(d)) · S(d)`

whose two factors are *detectability* and *residual validity given consistency*. They
answer different questions and neither is a nuisance parameter of the other.

## D.4 What the corrected numbers actually say

Median conditional survival is **0.952**: among corrupted knowledge sets that stay
internally consistent, the committed adjustment set usually stays valid. Combined with D.1,
the picture is a two-stage filter — most errors are self-revealing, and most of the
survivors are harmless. The dangerous case is **consistent *and* invalidating**, and it is
rare.

This makes the Phase 1 question sharper rather than weaker, and it is worth stating why the
two quantities can disagree without contradiction: `r_val` is a **worst-case** distance to
the nearest failing graph, while `S(d)` is an **average-case** survival probability over a
shell. `r_val = 1` is common in this population while average survival at `d = 1` is high —
both are true, because a failure existing at distance 1 says nothing about how much of the
shell fails. Whether the worst-case radius predicts average-case survival is exactly the
open question, and these numbers confirm it is non-trivial rather than answering it.

---

# Appendix E — 2026-09-09, a measurement defect in the tiered arm

## E.1 The defect

Found while checking why `mean_S` at `d = 0` was 0.973 rather than 1.0.

§1.1 specified the tiered arm as a *measured-depth* design: `d_claims` = the number of
assertions in `K_cor` "that contradict or are absent from `K_ref`". As implemented this
counts `|K_cor \ K_ref|` — assertions **added or reversed** — and therefore **misses
assertions removed**. `tiered` corrupts by relocating nodes between tiers; when relocation
puts two nodes in the *same* tier, their cross-tier assertion simply disappears. `K_cor` is
then a strict subset of `K_ref`, so `|K_cor \ K_ref| = 0` while the resulting MPDAG is
strictly less oriented than `G₀` and can invalidate `Z*`.

Evidence (only the tiered arm has `d = 0` rows; the flip arm's depth is targeted and exact):

- 408 of 668 tiered `d = 0` rows have `S(0) < 0.999`, with `n_eval` up to 994 — far above
  the 200 reps a single depth should hold.
- The `d = 0` bin is contaminated by every corruption rate: 133,600 samples at rate 0.0
  (all surviving, correct) plus 86,557 at 0.05, 42,185 at 0.10, 22,947 at 0.15, and so on.
- Those contaminants include genuine failures — 1,385 at rate 0.20, 1,371 at 0.25.

`s0_ok` is `True` for all 1,978 instances and `n_s0_failures = 0`, so **`G₀` itself is
sound and A2 is not violated**. The defect is in the depth *metric*, not in the graphs, the
oracle, or the gate.

## E.2 The fix — no regeneration required

The tiered arm is re-binned on its **native experimental knob, `corruption_rate`**, which
is exactly recorded per sample in `survival_samples.csv`. The result is a clean, balanced,
monotone design:

| corruption_rate | 0.0 | 0.05 | 0.10 | 0.15 | 0.20 | 0.25 | 0.30 | 0.35 | 0.40 | 0.45 | 0.50 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| contradiction rate | 0.000 | 0.229 | 0.404 | 0.574 | 0.648 | 0.692 | 0.742 | 0.801 | 0.813 | 0.833 | 0.852 |
| `S` given consistent | 1.000 | 0.994 | 0.963 | 0.937 | 0.898 | 0.879 | 0.848 | 0.753 | 0.740 | 0.671 | 0.618 |
| n | 133,600 at every rate |

*(Corrected 2026-09-09: this row originally printed 0.798 at rate 0.35, a transcription
error by the orchestrator. Independent recomputation during the re-binning gives 0.7526,
confirming 0.753. Verified against `analysis_tiered_rebinned.csv`.)*

`S = 1.000` at rate 0.0 exactly as A2 requires, both series monotone, n balanced. This axis
is better than the one pre-registered: it is the knob actually turned, needs no post-hoc
measurement, and cannot be corrupted by a counting rule.

**Scope of the correction.** The **flip arm is unaffected** — its depth is targeted
(`rate = d/|K|`, asserted exact at generation) and it has no `d = 0` rows. All flip-arm
numbers in Appendices C, D and the run notes stand.

## E.3 RETRACTED — the cross-arm detectability comparison

D.1 and C.3 stated that correlated (tiered) errors are markedly harder to detect than
uniform ones, quoting 0.354 versus 0.653 at `d = 1`. **The tiered half of that comparison
came from a contaminated bin and is withdrawn.**

The within-arm gradient survives and is solid: tiered contradiction rises monotonically
from 0.000 to 0.852 across the corruption grid. But the arms cannot be matched — one
relocates a fraction of *nodes*, the other reverses a count of *claims*, and §1.1 already
declares them separate instance populations. **No cross-arm detectability claim is
supported by this design**, and none is made. Establishing one needs a corruption process
parameterised identically across both, which is a design change, not a reanalysis.

**The C.3 finding itself stands on the flip arm alone**, where it is unaffected: reversing
a single truthful claim renders the knowledge set inconsistent 65.3% / 64.7% / 56.3% of the
time at base wrongness 0.0 / 0.10 / 0.25. Most orientation errors are self-revealing.

## E.4 Consequence for P4

P4 predicted higher across-instance variance in `AUC_frac` for tiered than flip. Tiered
`AUC_frac` as originally computed rests on the defective axis, so P4 must be evaluated on
the re-binned rate axis — and the two arms' AUCs are then in different units (fraction of
nodes relocated vs fraction of claims reversed). P4 is therefore reported as
**not comparable as pre-registered**, with both arms' variances given separately and the
unit mismatch stated, rather than forced into a verdict.

---

# Appendix F — 2026-09-09, effective sample size behind `AUC_frac`

## F.1 One of the two discordance spotlights does not survive verification

The analysis verified that each spotlight instance's `AUC_frac` was *recomputed correctly*
from the curves. It did not check how much data each grid point rested on. Doing so:

**`c06_s01_cov050_seed00000007_bw000` — SOLID, and it is the headline case.**
`shd_truth = 0` (so `G₀` *is* the ground-truth DAG — the analyst's knowledge is perfectly
correct), `n_k = 2`, `r_val = 1`, `AUC_frac = 0.0`. `S(1) = 0.0` on 107 non-contradictory
samples and `S(2) = 0.0` on 200. Every single corruption invalidates the committed set.
SHD says "perfect"; `r_val = 1` says "one move from invalid"; the survival curve agrees with
`r_val`. Zero SHD is compatible with maximal fragility.

**`c12_s11_cov100_seed00000134_bw025` — WEAK, and is withdrawn as evidence.**
`AUC_frac = 1.0` is computed from depths with `n_eval` of 35, 10, **1** and **1**, with
several depths at 0. Four surviving draws are not a survival curve. It is not quoted.

## F.2 The systematic version of the problem

Effective sample size behind `AUC_frac` varies across instances, and correlates with
`r_val` **in different directions in different strata**:

| stratum | τ_b(`r_val`, usable depths) |
|---|---|
| flip cov=0.5 bw=0.00 | **+0.318** |
| flip cov=0.5 bw=0.10 | +0.192 |
| flip cov=0.5 bw=0.25 | −0.100 |
| flip cov=1.0 bw=0.00 | −0.245 |
| flip cov=1.0 bw=0.10 | −0.252 |
| flip cov=1.0 bw=0.25 | −0.308 |
| tiered (pooled) | −0.438 |

And `AUC_frac` rises with the number of usable depths (median 0.60 at 1 usable depth,
0.90 at 3; n = 386 and 667).

So the bias is **not uniform**:

- Where τ(`r_val`, usable) < 0 — all coverage-1.0 flip strata and tiered — high-`r_val`
  instances rest on *less* data, which pushes their measured AUC *down*. The bias runs
  **against** the hypothesis; those strata are conservative.
- Where it is positive — flip cov=0.5 at base wrongness 0.00 and 0.10 — the bias runs
  **for** the hypothesis. These are the strata carrying the two largest τ values
  (+0.714, +0.380), so those two specifically may be inflated.

## F.3 Control and consequence

The pre-specified sensitivity `AUC_frac_usable` (depths with `n_eval ≥ 30` only) is exactly
the right control. Under it: **P1 falls from 9/9 to 8/9 strata and P2 from 9/9 to 7/9.**
The verdicts weaken but hold; 23 instances (1.2%) have no usable depth and become NaN.

**Reporting rule adopted:** both endpoints are reported side by side wherever P1 or P2 is
quoted, and the `AUC_frac_usable` figure is the one quoted in any abstract or summary claim,
because it is the conservative one. `AUC_frac` alone is not quoted as a headline.

## F.4 Housekeeping

`run_analyse.py` sits at the repository **root** rather than in
`src/bkrobust/robustness/`. On inspection it is root-anchored *by construction* — it
resolves both `src/` and `results/axis_robustness/` via `Path(__file__).parent` — so it
works as placed and must not be moved without editing those paths. It differs from the
repo's convention of putting drivers inside the package (`benchmarks/measure.py`,
`synth/runner.py`); noted rather than "fixed", since relocating a working driver to satisfy
a convention is how this repo has broken things before.
