# Synthetic-graph generalisation and search heuristics for the breakdown radius

*Session report. Every number traces to a committed file under `results/`. Where
a run was cut short or a scope was bounded, this says so. Negative and null
results are reported at the same weight as positive ones, and the places where I
got something wrong are recorded with how it was caught.*

Branch: `experiments/synth_graphs_and_heuristics`. Nothing was committed to `main`.

---

## 0. The headline results

| | Result |
|---|---|
| **Radius convention** | Resolved: `r` = distance to the nearest failure. Shells `0…r−1` clean; the practitioner's safe-move count is `r−1`, **not** `r`. Prior numbers carry over unchanged. |
| **Conjecture 1** (failure upward-closed) | **It is a theorem**, with a one-line proof. Not a conjecture. |
| **Conjecture 2** (retraction-optimal witnesses) | **No counterexample in 1,985,314 exhaustive radius comparisons.** Stated as an empirically-supported conjecture, not proved. |
| **Exact accelerated search** | Space-free method reproduces BFS on **5,304/5,304** differential cases. Speedup 25× to **133×**, growing with component size. |
| **`L = U = r` certificate** | Exact radius certified with **no enumeration in 83.3%** of instances. |
| **H1 saturation** | `r_val = 1` in ~**80%** of instances — dominant but well short of the >90% that would make the radius vacuous. |
| **H2/H3 frontier** | **The frontier is non-flat in 21.9% of instances — but `O*` is never strictly beaten (0 of 249,732).** No efficiency–robustness tradeoff. |
| **H5 middle regime** | Strongly ε-dependent: 0.0% at small ε, 14.1% at ε=0.2. Pre-registered prediction of 20–50% **not supported**. |
| **H6 calibration** | Coverage **1.0000** everywhere (the correctness gate). Conservativeness mean **0.181**, predicted by knowledge-component size (Spearman 0.564). |

---

## 1. Task 0 — the radius convention

The report and the framing document disagreed by one. Resolved in
`src/bkrobust/core/conventions.py` as the single normative definition:

> `r_P = min { d(G0, G) : P fails at G }` — **the distance to the nearest failure.**
> Shells `0 … r−1` are certified clean. The practitioner's "how many moves am I
> safe for" is `r − 1`.

Chosen because (a) it is a minimum over a set, a clean object; (b) Axis B's lower
bounds are naturally `r ≥ L`, certifying shells `0…L−1` without visiting them,
which composes directly with this reading and awkwardly with the other; and (c)
it leaves every number in the existing `report.md` unchanged. Under the other
convention every previously published number would have shifted by one.

This is asserted element-by-element, not merely documented:
`test_shells_below_radius_are_genuinely_clean` walks every element of every shell
below `r` and checks it really is clean, and that shell `r` really does fail.

The frozen core (`oracle`, `spacelib`, `instance`, `resultsio`) reproduces the
prior report's radii **3 / 3 / 2** exactly, which is the integration check that
the new machinery agrees with the old.

---

## 2. Axis B — exact search

### 2.1 Conjecture 1 is a theorem

> If `Z` fails in `G` and `[G] ⊆ [G']`, then `Z` fails in `G'`.

Validity is a for-all over represented DAGs. `Z` failing in `G` means some
`D ∈ [G]` has `Z` invalid; since `[G] ⊆ [G']`, that same `D` lies in `[G']` and
witnesses failure there. ∎

It was posed as a conjecture; it is one line. Reporting it as an empirical
finding would have overstated what was learned.

**One caveat, and it is a convention artefact, not mathematics.** This codebase
defines `is_valid` to return `False` when `[G]` is *empty* — a graph representing
no model cannot certify anything. Under that convention an empty-extension graph
"fails" vacuously while `[] ⊆ [G']` for every `G'`, so the implication would
break. It does not bite because `enumerate_space` admits no such element, which
is checked (`check_no_empty_extensions`) rather than assumed.

Consequences the search exploits: the failing set is an **up-set**; its boundary
is an antichain of minimal failures; nothing above a known failure ever needs
expanding.

### 2.2 Conjecture 2 — no counterexample in ~2M comparisons

> The nearest failure is reachable from `G0` by retractions alone, so
> `r_val` = the minimum number of retractions inducing failure.

This is **not** implied by Conjecture 1: upward-closure says failures persist as
you move up, not that the nearest one lies above `G0`. A failure incomparable to
`G0` at distance 2, with the nearest up-set failure at 3, would refute it.

Exhaustive search over every CPDAG, every space element taken as `G0`, every
ordered `(X, Y)`, every valid `Z`:

| n | CPDAGs | spaces | radius comparisons | counterexamples |
|---|---|---|---|---|
| 3 | 11 | 7 | 150 | 0 |
| 4 | 185 | 126 | 13,344 | 0 |
| 5 | 8,782 | 6,030 | **1,971,820** | **0** |

**Scope, stated rather than implied.** At n=5, 2,616 CPDAGs have no undirected
edge (nothing for knowledge to orient), 131 exceed 6 undirected edges, and 5 had
a space above 200 elements — so **136 of the 6,166 knowledge-carrying CPDAGs were
skipped (2.2%)**, all of them the densest. The claim is therefore "no
counterexample over all CPDAGs on ≤5 nodes with ≤6 undirected edges and a space
of ≤200 elements", not "over all CPDAGs on 5 nodes". It remains a conjecture.

*Free validation.* The enumerator reproduces the published counts of labelled
DAGs (25 / 543 / 29,281) **and** of Markov equivalence classes (11 / 185 / 8,782)
on 3/4/5 nodes. Both are pinned in a test.

### 2.3 A space-free exact method

`radius_local_up` performs upward BFS generating covers **locally**, so the space
is never enumerated, and prunes using upward-closure (never expanding above a
known failure).

**Differential validation — 5,304 cases at n=3 and n=4, exhaustive, 0
disagreements** with full-space BFS (`results/search/differential.csv`).

The real correctness risk was cover generation: a missed cover would silently
inflate distances. Locally generated covers were compared **element by element**
against the space-derived covers and are identical.

**Speedup, accounted honestly.** Timing per query against a pre-built space
flatters BFS, because `local_up` never builds one. On the single-query comparison
(BFS pays construction; `local_up` does not):

| undirected edges `k` | mean space size | speedup |
|---|---|---|
| 2 | 6 | 25× |
| 4 | 32 | 24× |
| 6 | 97 | 81× |
| 7 | 177 | **133×** |

The advantage grows because construction is `3^k` while `local_up` is independent
of it. **Where it does not help:** when many queries share one CPDAG, BFS
amortises construction and the gap narrows sharply. On tiny instances with a
pre-built space `local_up` is actually *slower* (0.72s vs 0.36s over 5,208 n=4
cases). Both regimes matter and both are reported.

### 2.4 Bounds: exact certification without enumeration in 83.3%

Exhaustive over all n=4 CPDAGs, 792 instances with a finite exact radius:

- `L` (descendant-mode lower bound): **admissible in all 732 cases where
  defined, 0 violations**; tight (`L = r`) in **660**.
- `U` (guided upper bound): valid (`U ≥ r`) in **all 792**.
- **`L = U = r` in 660/792 = 83.3%** — five instances in six certified exactly
  with no enumeration at all.

**Limitation, stated not buried.** `L` covers only one of the two ways `Z` can
fail — a member becoming a possible descendant of `X`. In **60/792 (7.6%)** no
such member exists yet the radius is finite, because the binding mode is a
back-door path becoming unblocked, for which no bound is implemented. So `L` is
**not** a complete lower bound alone, and `L = UNREACHED` means "this mode gives
no bound", never infinity. 83.3% is the rate at which the *pair* certifies.

---

## 3. Axis A — ensembles

`results/synth/preregistration.md` was committed before the first large run and
has not been edited. Deviations are recorded here.

### 3.1 What was run

Grid: 12,960 points over Erdős–Rényi, scale-free and block generators, n=6…9,
four corruption operators (omission, flip, compound, tiered), three rates.
**1,064 accepted (8.2%)**, zero run errors, zero `r_val = 0` rows.

**The rejection rate is itself a result** (pre-registered as such):

| reason | count |
|---|---|
| empty set trivially valid | 6,889 |
| treatment not in/adjacent to a component | 3,107 |
| no atomic perturbation changes validity | 696 |
| `Z` not identified | 620 |
| `K_assumed` inconsistent with the CPDAG | 341 |
| CPDAG too large for BFS | 160 |
| `Z` invalid at `G0` | 83 |

**91.8% of randomly generated instances are degenerate for this study** — mostly
because there is no confounding to speak of. The phenomenon needs fairly specific
structure, and a paper should not present ensemble rates without this denominator.

Alongside the ensembles, an **exhaustive census** over every CPDAG on ≤5 nodes
(134,140 instances) gives generator-free evidence for H1 and H2.

### 3.2 H1 — saturation. Prediction upheld; sub-prediction not.

| source | r=1 | r=2 | r=3 | r=4 |
|---|---|---|---|---|
| census, n=5 (88,700 finite) | **81.0%** | 15.9% | 3.1% | — |
| ensembles, n=6…9 (1,050 finite) | **78.6%** | 17.9% | 2.8% | 0.5% |

Remarkably stable across generators (76.9–80.4%) and across n (77.4–81.2%).

**Verdict: H1 not falsified.** The mode is 1, but roughly a fifth of instances
have radius ≥ 2 — short of the >90% saturation that would have made the radius
near-vacuous, which was the project's principal risk.

**The sub-prediction is not supported.** I predicted mass at `r ≥ 2` would grow
with the knowledge-intersected component. It does not: the distribution is
essentially flat in n across 6…9. Reported as a failed prediction, not omitted.

### 3.3 H2 / H3 — the frontier. My first answer was wrong.

This is the result I got wrong first, so the correction is given in full.

**The wrong version.** My initial `h2_frontier` reported *0 of 134,140* instances
with a strictly more robust valid set — an apparently decisive confirmation that
the prior report's flat frontier was general. It was an artefact of my own
analysis code, two ways: it **excluded every comparison in which either radius
was `UNREACHED`**, and the census **capped candidate sets at 10** sorted by size.
`UNREACHED` is not a number, but it is not missing data either — a set that never
fails anywhere is *strictly more robust* than one failing at distance 1, the
largest gap there is. Excluding those comparisons excluded exactly the non-flat
cases.

**The corrected version**, exhaustive at n=4 and n=5, no cap, `UNREACHED` ordered
above every finite radius:

| | n=4 | n=5 | total |
|---|---|---|---|
| instances with ≥2 valid sets and `O*` defined | 2,652 | 247,080 | **249,732** |
| frontier **non-flat** | 816 (30.8%) | 53,880 (21.8%) | **54,696 (21.9%)** |
| `O*` **strictly beaten** | **0** | **0** | **0 (0.00%)** |
| `O*` never fails anywhere | 1,584 (59.7%) | 147,740 (59.8%) | 149,324 (59.8%) |

(n=4 covers CPDAGs with ≤6 undirected edges, n=5 with ≤4, for cost.)

Excluding the empty adjustment set — trivially robust, and gated out of the
ensembles anyway — the n=4 non-flat rate is 10.96% of 2,628 instances, and `O*`
is still never beaten.

**The conclusion changes, and improves.** Valid adjustment sets genuinely *do*
differ in robustness, so the prior report's flat frontier is **not** a general
phenomenon. But the optimal set is never the less robust choice. So:

> **There is no efficiency–robustness tradeoff. The most efficient adjustment set
> is also among the most robust.**

That is a stronger and more useful statement than either "the frontier is flat"
or "a tradeoff exists", and it is directly actionable: an analyst choosing `O*`
gives up nothing in robustness.

**H3 and the designed families — a structural obstruction.** The pre-registered
plan was to build graphs with two back-door routes whose blocking sets lie in
different chordal components. This turned out to be **impossible for the
construction shape attempted**, for a reason worth recording. Identification
requires `X → Y` to be compelled, which requires a v-structure at `Y`; adding the
node that supplies it also freezes the confounder→`X`,`Y` edges, which makes `O*`
structurally invariant across the whole space — so *no* atomic perturbation
changes any validity status and every instance is gated out. The family went from
"zero valid adjustment sets" to "100% gate-rejected for the opposite reason".

There is a genuine tension between **identification** (needs compelled edges
around the target) and **perturbability** (needs undirected ones). I did not find
a construction satisfying both and am not claiming none exists.

**This did not block H3**, because the corrected H2 analysis shows designed
adversarial families are unnecessary: ~31% of *ordinary* small instances already
have a non-flat frontier. The designed experiment was superseded by an exhaustive
one, which is the better evidence anyway.

### 3.4 H4 — not run

The naive `K`-count radius was **not computed for the ensembles**. The runner's
schema carries no `naive_radius` column and I chose to spend the remaining effort
on the H2 correction instead. The prior report's single-example finding (naive
understates the model radius, off by one in all three scenarios) is therefore
**neither confirmed nor refuted** by this session. Listed in `NEXT.md` as the
cheapest remaining item.

### 3.5 H5 — the middle regime. Prediction not supported.

Fraction of instances with `r_ε > r_val` (fragile identification, stable
estimate), over 1,050 ensemble instances:

| ε | 0.02 | 0.05 | 0.1 | 0.2 |
|---|---|---|---|---|
| middle regime | 0.0% | 0.0% | 0.7% | **14.1%** |

I pre-registered 20–50%. **Not supported**: at small ε the middle regime
essentially does not occur, because bias is at machine zero while `Z` stays valid
and jumps as soon as it does not. Only at ε=0.2 — a fifth of the effect size —
does the regime appear at 14.1%.

**Implication for the project:** `r_ε` carries much less independent information
than the single worked example suggested. It is nearly redundant with `r_val`
unless ε is set at a substantial fraction of the effect.

### 3.6 H6 — calibration

89 instances across the three generators:

- **Coverage = 1.0000, minimum 1.0000, in every instance.** This is the
  correctness gate, not a finding: anything below 1.0 would be a bug and halt the
  run.
- **Conservativeness: mean 0.181**, range 0.000–0.625 (block 0.236, Erdős–Rényi
  0.189, scale-free 0.149).
- **Spearman(undirected-component size, conservativeness) = 0.564** — the
  knowledge-intersected component predicts conservativeness, as H6 predicted.

Conservativeness here (0.18) is much lower than the prior report's single-example
0.43–0.47, i.e. the certificate is *less* pessimistic on these ensembles than the
worked example suggested.

---

## 4. Bugs found, and how each was caught

House style from the prior report: record them, including my own.

| # | Bug | How caught |
|---|---|---|
| 1 | `.gitignore` silently excluded `results/synth/` and `results/search/`, so a commit claiming to add the differential rows added nothing | `git commit` reported "nothing to commit" for a brand-new file |
| 2 | **My H2 statistic discarded `UNREACHED` comparisons**, producing a confident and wrong "0 of 134,140" | Re-derived the frontier from scratch without the cap and got 23% non-flat |
| 3 | Conjecture-1 check (O(\|space\|²)) run once per combination; one dense CPDAG burned 4 minutes, twice | n=5 sweep stalled at the same index on two independent runs |
| 4 | Space-size guard ran *after* the cubic covering-relation build | Same investigation |
| 5 | I pasted a garbled graph (both `V1->V6` and `V6->V1`) into a subagent brief | The subagent refused to trust it and reconstructed the mechanism itself; the constructor provably forbids mutual pairs |
| 6 | `Meek(Ĉ,K)` can fall outside `enumerate_space` (0.09%) — chordality is checked on the undirected subgraph, so a component whose only chord is *directed* reads as chordless | A 12,960-point run crashed; incremental writing preserved 377 rows |
| 7 | 9% of accepted instances had `r_val = 0`, forbidden by the convention — the gate admitted target pairs with no causal path, giving an empty and invalid `O*` | Inspecting the output distribution rather than trusting it |
| 8 | The H3 family had **zero** valid adjustment sets for its designed `(X,Y)` — the session's key experiment was silently untestable | Reviewing the generator against its purpose before running it |
| 9 | The H3 fix then froze `O*` structurally, making the family 100% gate-rejected the other way | Found and reported by the subagent |

Items 1–5 are mine, 6 is inherited, 7–9 arose in delegated work. Bug 2 is the
most consequential: it would have put a confident false claim in the report.

**One issue deliberately left open.** Bug 6 points at a real tension in the
inherited validity definition: a Meek-closed graph representing 5 DAGs, keeping
every compelled edge and satisfying `[G] ⊆ [Ĉ]`, is nonetheless excluded from the
space. Relaxing `is_valid_mpdag` would invalidate every result resting on the
current enumeration — including the 1.99M-comparison conjecture sweep — so it is
gated, its rate is recorded, and the question is carried to `NEXT.md`.

---

## 5. What this does not show

- **Fixed skeleton, causal sufficiency, faithfulness, oracle CI testing** are
  assumed throughout. Finite-sample skeleton error is entirely outside this
  analysis, and in practice it may dominate everything measured here.
- **Small graphs.** Exhaustive results cover n ≤ 5; ensembles reach n = 9. Nothing
  here establishes behaviour at realistic sizes.
- **Conjecture 2 is not proved.** ~2M comparisons with no counterexample is
  strong evidence and not a proof, and the 2.2% densest CPDAGs at n=5 were not
  examined.
- **`L` is not a complete lower bound.** It covers one failure mode; 7.6% of
  instances fail by the other.
- **The 83.3% certificate rate is an n=4 number.** It is not established at
  larger sizes.
- **Bias is a sampled proxy.** Radii key off *mean* bias; the maximum is a
  sampled statistic, measurably unstable, and is never quoted as a bound.
- **H4 was not run at all.**
- **Generator-specific vs general.** H1's ~80% and H2's "`O*` never beaten" hold
  across all three generators *and* the exhaustive census, so they are not
  artefacts of one family. H5's ε-dependence and H6's conservativeness levels are
  ensemble-specific and should not be quoted as general.
- **`O*` is never beaten** rests on 249,732 exhaustive instances at n=4 and n=5.
  It is an empirical regularity, not a proof, and a single counterexample would
  matter. It has not been tested above n=5 or on dense CPDAGs.

---

## 6. Reproduction

```bash
PYTHONPATH=src python3 -m pytest tests/core tests/search tests/synth tests/demo -q   # 277 passing
```

(The three root-level scaffold test modules require Python 3.11 — `StrEnum`,
`TypeAlias` — and do not collect on this interpreter. That is by design and
predates this session.)

Key runs, all seeded from root seed `20260919`:

```bash
PYTHONPATH=src python3 -c "from bkrobust.search.conjecture_study import run_study; run_study(4,'results/search/conjectures')"
PYTHONPATH=src python3 -c "from bkrobust.search.differential import differential_sweep; differential_sweep(4)"
PYTHONPATH=src python3 -c "from bkrobust.search.radius_census import census; census(5)"
PYTHONPATH=src python3 -c "from bkrobust.search.frontier import frontier_sweep; frontier_sweep(4)"
```

**Environment.** Python 3.9.6, numpy 2.0.2, networkx 3.2.1, pandas 2.3.3,
scipy 1.13.1, matplotlib 3.9.4.

**Runtimes.** n=5 conjecture sweep 926 s; n=5 census 107 s; Axis A grid 132 s;
differential sweep 2.5 s; frontier n=4 ~110 s.

**Determinism.** All output is bit-identical across runs and across
`PYTHONHASHSEED`, enforced by regression tests in `tests/core/test_determinism.py`
and `tests/synth/`. Two inherited bugs of exactly this kind (RNG consumed in
frozenset order; float summation in set order) were fixed in the prior session
and the guard was extended to all new code here.
