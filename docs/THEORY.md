# Theory

Definitions, notation, and the three theorem targets. This document is the
**project's original planning note**; it states targets, not results.

> **Status and precedence (updated).** Since this note was written, the
> order-theoretic programme has been settled and now lives in
> [`THEOREMS.md`](../THEOREMS.md), with the bias radius in
> [`R_EPSILON_THEORY.md`](R_EPSILON_THEORY.md). Where the two disagree,
> **`THEOREMS.md` wins**: it carries the Main Theorem hierarchy (`THEOREMS.md`
> §0b), the proof-structure diagram, and the parametric assumption the bias
> results require (`THEOREMS.md` §16b). This file is retained for the targets
> that remain genuinely open — T1 and T3 — and for the framing in §2, §5 and §6.

**Notation is standardised as follows.** Early drafts wrote radii as `δ`; the
settled notation is `r`, and the two are the same objects:

| This file (old) | Standard (use this) | Meaning |
|---|---|---|
| `δ_valid` | `r_val` | validity / breakdown radius |
| `δ_opt` | `r_opt` | optimality radius |
| — | `r_ε` | ε-bias radius (`R_EPSILON_THEORY.md`) |
| `M = M(C, K)` | `G`, with `G₀` the analyst's state | a knowledge state (MPDAG) |
| `C` | `Ĉ` | the estimated CPDAG |
| `d(·,·)` | `d(·,·)` on the covering graph of `⪯` | the poset distance; `shd` is a baseline, not the metric |

The radius convention is the one fixed in `THEOREMS.md` §0:
`r = min{ d(G₀,G) : the property fails at G }`, so the closed ball of radius
`r − 1` is certified clean and the practitioner's safe-move count is `r − 1`.
A **ball** is `{G : d ≤ k}`; a **sphere** (loosely, a "shell") is `{G : d = k}`.
Guarantees are stated on balls, never on shells.

---

## 1. Notation

| Symbol | Reads as | Notes |
|---|---|---|
| `G` | the true DAG | Unknown in practice; known in simulation. |
| `C` | the CPDAG | The Markov equivalence class of `G`, as output by discovery. |
| `K` | background knowledge | Required/forbidden edges, tiers, ancestral constraints. |
| `M = M(C, K)` | the MPDAG | Meek closure of `C` under `K`. Undefined when Algorithm 1 FAILs. |
| `X`, `Y` | treatment, outcome | The target pair. |
| `cn(X, Y)` | causal nodes | Nodes on a proper possibly-causal path from `X` to `Y`. |
| `forb(X, Y)` | forbidden set | Possible descendants of `cn(X, Y)`, plus `X`. |
| `Z` | an adjustment set | Any covariate set. |
| `O*(M)` | optimal adjustment set | Henckel–Perković–Maathuis, read off `M`. |
| `asVar(Z; θ)` | asymptotic variance | Of the `Z`-adjusted estimator under SCM `θ`. |
| `d(·,·)` | a distance | `shd`, `orientation_flips`, or `poset`. See `graphs/distances.py`. |
| `δ` | perturbation radius | Distance from true knowledge, in the metric `d`. |
| `K_δ(C, G)` | the knowledge ball | `K` consistent with `C`, false w.r.t. `G`, at distance `δ`. |
| `δ_opt` | optimality radius | Smallest `δ` at which `O*` stops being optimal. |
| `δ_valid` | validity radius | Smallest `δ` at which `O*` stops being valid. |
| `τ`, `τ̂` | true effect, estimate | Total effect of `X` on `Y`. |

---

## 2. The gap this project addresses

Background knowledge is checked for **consistency** and never for **truth**.

- **Consistency.** `K` is consistent with `C` if some DAG in the equivalence
  class of `C` satisfies every constraint in `K` — equivalently, if Meek's
  Algorithm 1 terminates without FAIL. This is computable from observed data.
- **Truth.** `K` is true if `G` satisfies every constraint in `K`. This requires
  `G`, which is exactly what discovery was run to find.

The set of **consistent-but-false** knowledge is therefore non-empty for any
`C` with at least one undirected edge, and by construction its members are
indistinguishable from true knowledge given the data. No amount of care with
existing tooling detects them.

The three-node chain is the whole problem in miniature. `A - B - C` is the
CPDAG of both `A → B → C` and `C → B → A`. An expert asserting the reverse of
the truth passes every check, and is wrong about both edges.

---

## 3. The two breakdown radii

Fix `C`, `G`, `X`, `Y`, a distance `d`, and (for optimality) an SCM `θ`.

```
δ_valid = min { δ : ∃ K ∈ K_δ(C, G) such that O*(M(C, K)) is not valid for (X, Y) in G }
δ_opt   = min { δ : ∃ K ∈ K_δ(C, G) such that O*(M(C, K)) is not optimal for (X, Y) in G under θ }
```

Both are **worst-case over the ball**: the smallest radius at which *something*
can go wrong. That is the right shape for a guarantee ("below this radius you
are safe regardless of which wrong claims were made") and the wrong shape for a
forecast, so the experiments report the distribution over the ball alongside the
minima — see `radius_distribution` in `theory/radius.py`.

Three regimes follow:

| Regime | Condition | Consequence | Detectable by the analyst? |
|---|---|---|---|
| Benign | `δ < δ_opt` | `O*` unchanged | n/a — nothing happened |
| Silent degradation | `δ_opt ≤ δ < δ_valid` | unbiased, larger variance | **No** |
| Bias | `δ ≥ δ_valid` | biased at any `n` | **No** |

The middle regime is the one worth naming. The estimate is consistent, the
diagnostics are clean, and the analyst is paying in effective sample size with
no signal that they are doing so.

---

## 4. Theorem targets

### T1 — Ordering of the radii — **OPEN**

> For any `C`, `G` in its equivalence class, target pair `(X, Y)`, and SCM `θ`:
> `δ_opt ≤ δ_valid`.

**Intuition.** Optimality is a property of the graph *and* the SCM; validity is
a property of the graph alone. A perturbation that changes `O*` at all
generically changes its variance, whereas only perturbations that change `O*` in
particular ways — introducing a forbidden node, or dropping a node needed to
block a back-door path — break validity. The optimality-breaking set should
therefore be the larger one, and its minimum the smaller radius.

**Where it could fail.** The argument above is generic, not universal. There may
be perturbations that break validity while leaving the variance exactly
unchanged — ties in `asVar` are not measure-zero on structured SCMs. If so, T1
holds only under a genericity condition on `θ`, and stating that condition
precisely is part of the target.

**Status.** Unproved. `verify_ordering` in `theory/radius.py` checks it
empirically, and `tests/test_radius.py::test_ordering_conjecture_on_canonical_graphs`
checks it on the canonical graphs. An exact-method violation is a counterexample
and is more valuable than a proof of the conjecture; the sweeps persist the full
witness for any such case.

### T2 — Bias bound beyond `r_val` — **SUPERSEDED (settled in a different form)**

> **Resolved, but not as stated.** The programme this target describes was
> carried out in `R_EPSILON_THEORY.md`, which supplies an *exact* worst-case
> bias functional `B(G)` and the radius `r_ε` at which it first exceeds a
> tolerance, together with an attained band on the realised error
> (`THEOREMS.md` §17–§21b). It does **not** take the Cinelli–Hazlett partial-`R²`
> form below; the certificate is an exact supremum over the ambiguity the
> knowledge state leaves open, not a two-parameter bound.
>
> **It carries a parametric assumption that this target left implicit**: a
> nonsingular jointly Gaussian observational law from a linear-Gaussian SCM
> (`THEOREMS.md` §16b). The original phrasing "in the linear-Gaussian case" was
> right to name the class; later write-ups dropped it and claimed the bias
> results needed no assumption beyond `r_val`'s. That claim is withdrawn.
>
> The partial-`R²` reparameterisation remains an open and worthwhile target, for
> the reason given below — analyst-readable units. Original statement:

> For `δ ≥ δ_valid` in the linear-Gaussian case, the bias of the `O*(M(C,K))`-adjusted
> estimator is bounded by a function of two partial `R²` quantities — the
> mishandled variables' association with the treatment and with the outcome —
> and the bound is attained.

**Why this form.** It is the Cinelli–Hazlett omitted-variable-bias
parameterisation, and it is chosen because both parameters live on `[0, 1)` and
can be reasoned about by an analyst with no access to `G`. A bound stated in
terms of edge coefficients would be correct and useless to the person who needs
it.

**Where it could fail.** An invalid adjustment set can err by *inclusion* — a
collider, or a descendant of the treatment — not only by omission, and the OVB
factoring is derived for the omission case. The inclusion case may need its own
bound, in which case T2 splits.

**Status.** Superseded as above; the partial-`R²` form is unproved. See
`theory/bias_bound.py`.

### T3 — Efficiency gap within the band — **OPEN**

> For `δ_opt ≤ δ < δ_valid`, the variance ratio
> `asVar(O*(M(C,K))) / asVar(O*(M(C,K_true)))` is bounded above by a function of
> the cascade size of `K` and the graph's local structure at `cn(X, Y)`.

**Why cascade size.** The perturbation itself is not what costs variance — what
costs variance is how many orientations near the causal nodes the perturbation
*forces*. A large cascade far from `cn(X, Y)` should be free. If that is right,
the bound depends on the cascade's intersection with the neighbourhood of
`cn(X, Y)` rather than on its total size, and the statement above is the coarse
version of a sharper one.

**Status.** Unproved. See `theory/efficiency_gap.py`.

---

## 5. The certificate

T1–T3 all require `G`. A practitioner has `C`, `K`, and data.

The certificate (`theory/certificate.py`) is the usable form. It reports the
radii **relative to the analyst's own knowledge as the reference point** rather
than relative to an unknown truth — which is well-defined, because `K_δ` is
built from the observed `C` and nothing in its definition needs to know which
member of the ball is true.

It answers: *how wrong would my knowledge have to be before my conclusion
changes?* Not *is my knowledge right?* The second question is unanswerable; the
first is the one sensitivity analysis has always asked in place of it.

---

## 6. Open questions

1. Does `δ_valid` concentrate as graphs grow, or does the worst case stay
   pathological? If it concentrates, the worst-case radius is a useful summary;
   if not, the distribution is the only honest report.
2. Are hub-heavy graphs systematically more fragile, as the cascade argument
   suggests? Scale-free versus Erdős–Rényi at matched density is the test.
3. Do amortized models inherit these radii, or fail differently? They have no
   `O*` to break, so any breakdown they show must arise by another route — and
   finding that route is the interesting outcome of the audit.
4. Is the tier kind reliably the most damaging per constraint asserted? It is
   the cheapest knowledge to state and forbids the most edges at once, which is
   a bad combination if true.
5. Can `δ` be estimated from `C` and `K` alone, without `G`? This is what the
   conditional representation-learning component is really asking, and a
   negative answer is worth as much as a positive one.
