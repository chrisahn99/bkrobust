# The ε-bias radius `r_ε` — statements, proofs, and what each one assumes

Session 8, Axis B. This document extends `THEOREMS.md` from a *validity*
certificate to a *bias* certificate. It is written to the same discipline: each
item states its hypotheses and is either **proved** (a proof a referee can
check), **inherited** (proved elsewhere and cited), or **verified** (with the
verification scope stated). Nothing is labelled proved on the strength of
evidence alone.

The reader is assumed to have `THEOREMS.md` §0–§6 to hand. Notation is
unchanged: `Ĉ` is a CPDAG, `𝔊_Ĉ` its space of knowledge states, `[G]` the DAGs
an element represents, `K_G = dir(G) \ dir(Ĉ)` its knowledge set,
`G ⪯ H ⟺ [G] ⊆ [H]`, `d` the hop count on the covering graph, `G₀` the
analyst's state, `Z` the adjustment set read off `G₀`, and

```
r_val  =  min { d(G₀,G) : Z is not a valid adjustment set for (X,Y) in G } .
```

---

## 0. What is being added, and why the existing `r_eps` does not do it

`bkrobust.synth.runner._radius_eps` already computes a quantity called `r_eps`.
It thresholds `bias_stats(...).mean_abs_bias`: for each element `G`, a linear
SEM is **redrawn at random on every DAG of `[G]`**, the absolute bias of the
`Z`-adjusted estimand is averaged over extensions × draws, and the radius is
the distance to the first element whose average exceeds `ε`.

That quantity carries no guarantee, for three independent reasons.

1. **It is a mean, and a mean over a larger set can be smaller.** Retracting
   knowledge enlarges `[G]`, so the average runs over more extensions; adding a
   low-bias extension lowers the mean. The map `G ↦ mean_abs_bias(G)` is
   therefore **not monotone** under `⪯`, the failure set `{mean > ε}` is not
   upward-closed, and none of the order-theoretic machinery of `THEOREMS.md`
   applies to it. In particular the upward search is not exact for it. §7
   exhibits explicit counterexamples.
2. **The SEM is redrawn per element.** Two elements are then compared under
   *different* data-generating parameters, so the difference between their bias
   numbers confounds the knowledge state with the draw. The paper already flags
   the related point for `max_abs_bias` ("a sampled proxy for a worst case and
   not a supremum", App. C).
3. **It is not a bound on anything the analyst holds.** The analyst has one
   dataset and one number. An average over hypothetical redraws is not an upper
   bound on the error in that number.

This document replaces the mean by an exact **supremum over the ambiguity the
knowledge state actually leaves open, at the analyst's own fixed observational
law**. That object is monotone by construction, is computable in closed form
without enumerating DAGs, and upper-bounds the analyst's realised error. Every
structural result proved for `r_val` in `THEOREMS.md` then transfers to it
verbatim, because — as §3 makes precise — those results never used anything
about failure except that it is **upward-closed**.

---

## 1. Setting: the estimand, the ambiguity set, and the bias functional

Fix a query `(X,Y)` and a set `Z` read off `G₀`, with `X, Y ∉ Z`.

**Standing model.** The observational law `P` is that of a linear-Gaussian SCM
whose DAG lies in `[Ĉ]`; write `Σ` for its population covariance. (§6 states
what survives beyond linear-Gaussian.) For a set `S ⊆ V \ {X,Y}` write

```
β(S; Σ)  :=  the coefficient on X in the population least-squares
             regression of Y on {X} ∪ S ,
```

obtained by solving the normal equations on the corresponding submatrix of `Σ`.
This is the large-sample limit of OLS, not a sample quantity.

**The analyst's estimand.**

```
θ_Z  :=  β(Z; Σ) .
```

This is **one number**. It depends on `Z` and on the data, and on nothing else —
in particular not on which DAG is the truth. It is what the pipeline reports.

**The effect a candidate truth would assign.** For a DAG `D ∈ [Ĉ]`, the parent
set `pa_D(X)` is a valid back-door adjustment set in `D`, so under `D` the total
effect of `X` on `Y` is

```
τ_D  :=  β(pa_D(X); Σ) ,      with the convention τ_D := 0 when Y ∈ pa_D(X)
```

(the convention is the usual IDA one: if `D` makes `Y` a parent of `X` then `X`
has no causal effect on `Y`). The multiset `{τ_D : D ∈ [Ĉ]}` is exactly the IDA
multiset of Maathuis–Kalisch–Bühlmann; `{τ_D : D ∈ [G]}` is its MPDAG refinement
(Perković 2020).

**Definition 1 (effect ambiguity set).** For `G ∈ 𝔊_Ĉ`,

```
T(G)  :=  { τ_D : D ∈ [G] }  ⊆ ℝ .
```

**Definition 2 (worst-case bias functional).**

```
B(G)  :=  max { |θ_Z − τ| : τ ∈ T(G) }
       =  max ( |θ_Z − min T(G)| , |θ_Z − max T(G)| ) .
```

`B(G)` is the largest error the reported number `θ_Z` can carry if the truth is
any of the DAGs the knowledge state `G` still permits. `T(G)` is finite and
non-empty for every `G ∈ 𝔊_Ĉ` (no element has an empty extension set —
`THEOREMS.md` §2), so `B` is well defined and finite everywhere.

**Definition 3 (normalisation).** Because `ε` is wanted as a percentage, define

```
B_rel(G) := B(G) / |θ_Z|          (relative to the reported estimate; |θ_Z| > 0)
B_std(G) := B(G) · sd(X) / sd(Y)  (standardised units)
```

**Remark 1 (the normaliser must not depend on `G`).** Both normalisers above
are constants of the instance. This is not cosmetic: every monotonicity result
below is stated for `B` and inherited by `c·B` for any constant `c > 0`, and is
**false** for a normalisation that varies with `G` — e.g. dividing by
`mean{|τ| : τ ∈ T(G)}`, which is exactly the kind of "relative bias" that looks
natural and destroys the certificate. Normalise by the estimate or by the data;
never by the ambiguity set.

---

## 2. The two facts about `B` that everything rests on

### Theorem A (exact zero-bias inside the certified shell) — **PROVED**

**Statement.** Let `G ∈ 𝔊_Ĉ` with `d(G₀, G) ≤ r_val − 1`. Then

```
T(G) = { θ_Z }      and therefore      B(G) = 0 .
```

Consequently `sup { B(G) : d(G₀,G) ≤ r_val − 1 } = 0`: on the whole closed ball
of radius `r_val − 1`, including `G₀` itself, the `Z`-adjusted estimand carries
**exactly zero** identification bias — not "no significant bias", but zero, at
the population level, for every observational law compatible with `Ĉ`.

**Proof.** By the definition of `r_val` as a minimum, no element at distance
`≤ r_val − 1` fails, so `Z` is a valid adjustment set for `(X,Y)` in `G`;
validity in an MPDAG is a for-all over `[G]`, so `Z` is a valid adjustment set
in every `D ∈ [G]`. Fix such a `D`. In a linear SCM a valid adjustment set
identifies the total effect: `β(Z; Σ) = τ_D`. Hence `τ_D = θ_Z` for every
`D ∈ [G]`, so `T(G) = {θ_Z}` and `B(G) = max |θ_Z − θ_Z| = 0`. ∎

**Remark 2 (this is the strongest form the claim can take, and its limits).**
Theorem A is an *identification* statement. It says the estimand is exactly
right, so there is no bias term to trade off. It says nothing about sampling
error: the analyst's finite-sample estimate still has variance of order
`n^{-1/2}`, and inside the certified shell that variance is the *only* error.
That is precisely the practical content — inside the shell the analyst may stop
worrying about the graph and reason about standard errors alone, and, as the
objective for this session put it, there is no reason to recompute estimates at
each state, since they all coincide.

**Remark 3 (`r_val` is conservative for bias, and that is the safe direction).**
The converse of Theorem A is false in general. At distance exactly `r_val` some
element has `Z` invalid in some extension `D`, but it can still happen that
`β(Z;Σ) = τ_D` for the particular `Σ` at hand — invalidity is a statement about
all laws compatible with `D`, not about one. Such coincidences are non-generic
(they confine `Σ` to a proper algebraic subset), but they are possible. So

```
r_val  ≤  min { d(G₀,G) : B(G) > 0 } ,
```

with equality for generic `Σ`. Every `ε`-radius defined below therefore
satisfies `r_ε ≥ r_val`: the validity radius is a **universal lower bound** on
every bias radius, and reporting it understates robustness rather than
overstating it.

### Theorem B (monotonicity of worst-case bias) — **PROVED**

**Statement.** If `G ⪯ H` then `T(G) ⊆ T(H)` and hence

```
B(G)  ≤  B(H) .
```

In particular `B(G ∨ H) ≥ max(B(G), B(H))`, and `B(Ĉ) ≥ B(G)` for every `G`.

**Proof.** `G ⪯ H` means `[G] ⊆ [H]` by definition of the order, so
`T(G) = {τ_D : D ∈ [G]} ⊆ {τ_D : D ∈ [H]} = T(H)`: note `τ_D` depends on `D`
and `Σ` only, never on the ambient MPDAG, which is what makes the inclusion of
index sets an inclusion of value sets. A maximum over a subset of a finite set
is no larger. For the join, `G ⪯ G ∨ H` and `H ⪯ G ∨ H` (Theorem 2 of
`THEOREMS.md`); for `Ĉ`, it is the maximum of the order. ∎

**Corollary B1 (the ε-failure set is upward-closed).** For every `ε ≥ 0`,

```
F_ε  :=  { G ∈ 𝔊_Ĉ : B(G) > ε }
```

is an up-set: `G ∈ F_ε` and `G ⪯ H` imply `H ∈ F_ε`. ∎

**Corollary B2 (nesting in ε).** `ε ≤ ε'` implies `F_{ε'} ⊆ F_ε`. ∎

**Remark 4 (what is *not* monotone, and why the distinction is the whole
point).** Theorem B is monotonicity **along the order**, i.e. along retraction.
It does **not** say that the bias of an individual state grows with its distance
from `G₀`. Two states at the same distance are typically incomparable and have
unrelated biases, and a state at distance `d+1` can have strictly smaller bias
than one at distance `d`. §7 gives explicit counterexamples found by exhaustive
search. What is monotone in the distance is the **worst case over the ball**,
which is the object defined next — and it is the right object anyway, since a
certificate has to hold for whichever wrong knowledge the analyst actually
holds, not for an average one.

---

## 3. Lifting the `r_val` machinery: it only ever used upward-closure

The following theorem is the reason this session is cheap. `THEOREMS.md` §6
proves that the nearest *failing* state is reachable by pure retraction. Reading
that proof back, the only property of "failing" it consumes is Theorem 1 —
upward-closure. Nothing about adjustment sets, d-separation or validity enters.
So the result is really a statement about arbitrary up-sets.

### Theorem U (universal retraction theorem) — **PROVED** (from Lemma L, `THEOREMS.md` §5)

**Statement.** Let `Φ ⊆ 𝔊_Ĉ` be **any** up-set (`G ∈ Φ`, `G ⪯ H` ⟹ `H ∈ Φ`),
and suppose `Φ ≠ ∅`. Put `r_Φ := min { d(G₀,G) : G ∈ Φ }`. Then the minimum is
attained at a state comparable to `G₀`:

```
r_Φ  =  min { d(G₀,G) : G ∈ Φ and G₀ ⪯ G } ,
```

and `r_Φ` equals the least number of retractions from `G₀` that lands in `Φ`.

**Proof.** Let `G ∈ Φ` attain `d(G₀,G) = r_Φ` and put `J = G ∨ G₀`. Then
`G ⪯ J`, so `J ∈ Φ` because `Φ` is an up-set; and `G₀ ⪯ J` by definition of the
join. By **Lemma L** (`THEOREMS.md` §5), `d(G₀,J) ≤ d(G₀,G) = r_Φ`; by
minimality of `r_Φ` over all of `Φ`, `d(G₀,J) ≥ r_Φ`. Hence `d(G₀,J) = r_Φ` with
`J` comparable to `G₀`. The identification of `d(G₀,J)` with a count of
retractions is Proposition R below. ∎

**Proposition R (retraction parameterisation of the up-set) — PROVED.**
For `G₀ ∈ 𝔊_Ĉ`,

```
↑G₀ := { G ∈ 𝔊_Ĉ : G₀ ⪯ G }  =  { Meek(Ĉ, K_{G₀} \ S) : S ⊆ K_{G₀} } ,
```

and for `G ∈ ↑G₀`,  `d(G₀,G) = |K_{G₀}| − |K_G|`.

**Proof.** *(⊇)* For `S ⊆ K_{G₀}`, the set `K_{G₀} \ S` is consistent with `Ĉ`
(any DAG witnessing consistency of `K_{G₀}` also satisfies the subset), so
`H := Meek(Ĉ, K_{G₀} \ S) ∈ 𝔊_Ĉ`; and every `D ⊨ K_{G₀}` satisfies
`K_{G₀} \ S`, so `[G₀] ⊆ [H]`, i.e. `G₀ ⪯ H`.
*(⊆)* If `G₀ ⪯ G` then by **Lemma O** (`THEOREMS.md` §13) `dir(G) ⊆ dir(G₀)`,
hence `K_G ⊆ K_{G₀}`; taking `S = K_{G₀} \ K_G` and using that `G` is a fixpoint
of the closure gives `Meek(Ĉ, K_{G₀} \ S) = Meek(Ĉ, K_G) = G`.
*Distance.* By **Lemma R** (`THEOREMS.md` §4b) every covering step changes
`ρ(W) = |K_W|` by exactly one, so every path in the covering graph from `G₀` to
`G` has length at least `ρ(G₀) − ρ(G)`. By **Property S** (§4c) the poset is
upper semimodular, hence satisfies the Jordan–Dedekind chain condition, so the
interval `[G, G₀]` has a maximal chain of length exactly `ρ(G₀) − ρ(G)`, which
is a path of that length. ∎

**Corollary R1 (shell parameterisation).** For `0 ≤ d ≤ |K_{G₀}|`,

```
{ G ∈ ↑G₀ : d(G₀,G) = d }  ⊆  { Meek(Ĉ, K_{G₀} \ S) : S ⊆ K_{G₀}, |S| = d }
                            ⊆  { G ∈ ↑G₀ : d(G₀,G) ≤ d } .
```

**Proof.** The first inclusion is Proposition R with `S = K_{G₀} \ K_G`, which
has `|S| = d`. The second holds because `K_{Meek(Ĉ,K₀\S)} ⊇ K_{G₀} \ S` gives
`ρ ≥ |K₀| − |S|`, i.e. distance `≤ |S|`. ∎

The point of Corollary R1 is operational: **shell `d` can be generated directly
from the `d`-subsets of `K_{G₀}`, with no traversal of shells `< d`.** §5 uses
this.

---

## 4. The bias profile and the ε-radius

**Definition 4 (bias profile).** For `d ≥ 0`,

```
β(d)   :=  max { B(G) : G ∈ 𝔊_Ĉ,  d(G₀,G) ≤ d }        (full ball)
β↑(d)  :=  max { B(G) : G ∈ ↑G₀,   d(G₀,G) ≤ d }        (retraction ball)
```

**Definition 5 (ε-bias radius).** For `ε ≥ 0`,

```
r_ε  :=  min { d(G₀,G) : B(G) > ε }  =  min { d : β(d) > ε } ,
```

with `r_ε := ∞` when `B(G) ≤ ε` everywhere.

### Theorem C (the profile is a monotone staircase, computable upward) — **PROVED**

**Statement.**

1. `β` and `β↑` are non-decreasing, and `β↑ ≤ β`.
2. `β(d) = 0` for every `d ≤ r_val − 1`.
3. `β↑(d) = max { B(G) : G ∈ ↑G₀, d(G₀,G) = d }` for `0 ≤ d ≤ |K_{G₀}|`
   (**sphere sufficiency**): only the outermost shell need be evaluated.
4. `β↑(d) = max { B(Meek(Ĉ, K_{G₀} \ S)) : S ⊆ K_{G₀}, |S| = d }`
   (**direct shell evaluation**): shell `d` needs no shell below it.
5. `r_ε = min { d : β↑(d) > ε }` — the first crossing of `β` and of `β↑`
   coincide, even though the two profiles differ above it.
6. `ε ↦ r_ε` is non-decreasing, and `r_val ≤ r_0 ≤ r_ε` for every `ε ≥ 0`.

**Proof.**

*(1)* Balls grow with `d`, and `↑G₀ ∩ ball ⊆ ball`.

*(2)* Theorem A.

*(3)* `≥` is immediate since the sphere is contained in the ball. For `≤`, take
`G ∈ ↑G₀` with `d(G₀,G) = d' < d`. By Proposition R, `|K_G| = |K_{G₀}| − d' > 0`
whenever `d ≤ |K_{G₀}|`, so `G ≠ Ĉ`; the finite interval `[G, Ĉ]` therefore
contains an upper cover `G'` of `G`, which by Lemma R satisfies
`|K_{G'}| = |K_G| − 1`, i.e. `d(G₀,G') = d' + 1`, and `G ⪯ G'` gives
`B(G) ≤ B(G')` by Theorem B. Iterating `d − d'` times produces `G'' ∈ ↑G₀` at
distance exactly `d` with `B(G) ≤ B(G'')`. Hence the ball maximum is attained on
the sphere.

*(4)* By Corollary R1 the family `{Meek(Ĉ, K₀\S) : |S| = d}` is sandwiched
between the sphere at `d` and the ball at `d`; its maximum of `B` is therefore
between `max` over the sphere and `max` over the ball, which are equal by (3).

*(5)* Put `d* = min{d : β(d) > ε}`. Pick `G` with `d(G₀,G) = d*`, `B(G) > ε`.
Let `J = G ∨ G₀ ∈ ↑G₀`. Theorem B gives `B(J) ≥ B(G) > ε`, and Lemma L gives
`d(G₀,J) ≤ d*`, so `β↑(d*) > ε`. Conversely `β↑ ≤ β`, so `β↑(d) > ε` forces
`β(d) > ε`. The two first-crossing indices therefore agree. (This is exactly
Theorem U applied to `Φ = F_ε`, which is an up-set by Corollary B1.)

*(6)* Corollary B2 gives `F_{ε'} ⊆ F_ε` for `ε ≤ ε'`, so the minimum distance
into `F_{ε'}` is at least that into `F_ε`. The chain `r_val ≤ r_0` is Remark 3,
and `r_0 ≤ r_ε` is the case `ε' = ε, ε = 0`. ∎

### Theorem D (the practitioner's certificate) — **PROVED**

**Statement.** Let `D*` be the true DAG, `τ* = τ_{D*}` the true total effect,
and suppose the analyst's asserted orientations are wrong in at most `k` of
them — formally, that there is `S ⊆ K_{G₀}` with `|S| ≤ k` and `D* ⊨ K_{G₀} \ S`.
Then the realised identification error of the reported estimand obeys

```
|θ_Z − τ*|  ≤  β↑(k)  ≤  β(k) .
```

Consequently:

```
k < r_val        ⟹   θ_Z = τ*  exactly            (zero bias)
k < r_ε          ⟹   |θ_Z − τ*| ≤ ε
r_{ε_x} ≤ k < r_{ε_y}  ⟹   ε_x  <  β↑(k)  ≤  ε_y
```

so for an error budget of `k` claims the worst-case bias is pinned into the
half-open band `(ε_x, ε_y]`: **no larger than `ε_y`, and demonstrably able to
exceed `ε_x`.**

**Proof.** Put `H = Meek(Ĉ, K_{G₀} \ S)`. By Proposition R, `H ∈ ↑G₀` and
`d(G₀,H) ≤ |S| ≤ k`. Since `D* ⊨ K_{G₀} \ S` and `D* ∈ [Ĉ]`, Meek's theorem
gives `D* ∈ [H]`, hence `τ* ∈ T(H)` and `|θ_Z − τ*| ≤ B(H) ≤ β↑(k)` by
Definition 4. The three consequences follow from Theorem C(2) (which gives
`β↑(k) = 0` for `k < r_val`), from `k < r_ε ⟹ β(k) ≤ ε ⟹ β↑(k) ≤ ε`, and from
`k ≥ r_{ε_x} ⟹ β(k) > ε_x`, combined with Theorem C(5) which lets `β` be
replaced by `β↑` at the crossing. ∎

### Theorem D' (the band is tight, not merely an upper bound) — **PROVED**

**Statement.** For every `k ≥ 0` there exists a DAG `D ∈ [Ĉ]` such that

1. `D` is compatible with the observed law (it is Markov equivalent to the truth,
   so it fits `Σ` exactly);
2. `D` satisfies all but at most `k` of the analyst's orientations `K_{G₀}`; and
3. `|θ_Z − τ_D| = β↑(k)`.

So `β↑(k)` is not a bound that happens to hold: it is *attained* by a world the
analyst cannot rule out on the evidence they have, and it therefore cannot be
improved without an assumption beyond the CPDAG, the data, and the error budget.

**Proof.** Let `H ∈ ↑G₀` with `d(G₀,H) ≤ k` attain `β↑(k)`, and let `D ∈ [H]`
attain the maximum in `B(H) = max_{D' ∈ [H]} |θ_Z − τ_{D'}|`. Then (3) holds by
construction. `D ∈ [H] ⊆ [Ĉ]`, and every element of `[Ĉ]` is Markov equivalent
to the truth and hence fits `Σ`, giving (1). For (2): `K_H ⊆ K_{G₀}` by
Proposition R, and `|K_{G₀}| − |K_H| = d(G₀,H) ≤ k`, so `D ⊨ K_H` fails at most
`k` of `K_{G₀}`. ∎

**Corollary (the band is exact).** With `r_{ε_x} ≤ k < r_{ε_y}`, the quantity
`β↑(k)` reported in Theorem D is the exact worst case over the budget, and the
band `(ε_x, ε_y]` is the tightest one the supplied grid can express. Refining
the grid refines the band; nothing else does.

**Remark 5 (why the certificate is sound and not merely worst-case).** The
truth `D*` is *inside* the ambiguity set `T(H)`, by the displayed argument. So
`B(H)` is a genuine upper bound on the analyst's realised error, not a bound on
a hypothetical one. This is what distinguishes `B` from the sampled
`mean_abs_bias`: the latter is an average over draws none of which need be the
truth.

**Remark 6 (units).** `k` counts orientations in the **closure** `K_{G₀}`, the
same unit as `r_val` in the paper. The claim-level analogue `r_claim` (paper
§6.3) needs the same statement with `S` ranging over sets of *asserted claims*
and `K_{G₀}` replaced by their closure; every step above goes through unchanged
because Proposition R only needs `S` to index a subset of a consistent
orientation set, so `r_ε^claim ≥ r_claim` holds by the same argument.

---

## 5. Computing `r_ε` — the algorithm and its bounds

Theorem C converts "compute `r_ε` for arbitrary `ε`" from a search problem per
`ε` into a single monotone-staircase computation. Five facts do the work.

**(i) The search space is the up-set, not the space.** By Theorem U, only
`↑G₀ = {Meek(Ĉ, K₀\S)}` need be searched: at most `2^{|K_{G₀}|}` closures rather
than `|𝔊_Ĉ|` elements, and far fewer after deduplication, since many `S` share a
closure. This is the same reduction the paper already exploits for `r_val`, now
licensed for `r_ε` by the same theorem.

**(ii) Shells below `r_val` are free.** By Theorem A, `B ≡ 0` there. The
expensive part of an `r_ε` computation — evaluating `B` — is **never run inside
the certified shell**. The incremental cost of `r_ε` over `r_val` is therefore
only the shells `r_val … r_ε`. `r_val` itself comes from the existing
accelerated hybrid (`bkrobust.hybrid.breakdown_radius`), unchanged.

**(iii) Only the frontier shell is evaluated.** Theorem C(3): `β↑(d)` is the max
over the sphere, not the ball. Interior states need no `B` evaluation at all.

**(iv) One traversal answers every `ε`.** The staircase `β↑(r_val), β↑(r_val+1),
…` is computed once; `r_ε = min{d : β↑(d) > ε}` is then a lookup for **any**
`ε`, including `ε` chosen after the fact. There is no per-`ε` search, and
consequently no cost in supporting a whole grid `ε_x < ε_y < …` — which is what
the interval certificate of Theorem D needs.

**(v) `B` is evaluated semi-locally, without enumerating `[G]`.** See
Theorem E.

### Theorem E (semi-local evaluation of the ambiguity set) — **PROVED**

**Statement.** For `G ∈ 𝔊_Ĉ`, let `pa_G(X)` be the parents of `X` in `G` and
`nb_G(X)` its undirected neighbours. Then

```
{ pa_D(X) : D ∈ [G] }  =  { pa_G(X) ∪ S :  S ⊆ nb_G(X),
                            Meek( Ĉ,  K_G ∪ {s→X : s ∈ S} ∪ {X→t : t ∈ nb_G(X)\S} )
                            ≠ FAIL } ,
```

and hence `T(G) = { β(pa_G(X) ∪ S; Σ) : S admissible }`.

**Proof.** *(⊆)* Let `D ∈ [G]` and put `S = pa_D(X) \ pa_G(X)`. Every parent of
`X` in `G` is a parent in `D` and every non-neighbour relation is inherited, so
`S ⊆ nb_G(X)` and `pa_D(X) = pa_G(X) ∪ S`. `D` orients `s→X` for `s ∈ S` and
`X→t` for `t ∈ nb_G(X)\S`, and `D ⊨ K_G`, so `D` witnesses the consistency of
the displayed orientation set with `Ĉ`; by Meek's theorem the closure does not
FAIL.
*(⊇)* If the closure does not FAIL, then by Meek's theorem it has a non-empty
extension set, and any `D` in it satisfies `K_G` (so `D ∈ [G]`) and has
`pa_D(X) = pa_G(X) ∪ S` exactly. ∎

**Cost.** `2^{|nb_G(X)|}` Meek closures, each polynomial, with no DAG
enumeration anywhere. `|nb_G(X)|` is the number of *undirected* edges at the
treatment, which is small in every corpus this project measures. The exact
enumeration route (`enumerate_dag_extensions`) is retained and used as a
differential test, never on the hot path — the same discipline as Lemma O in
`THEOREMS.md` §13.

### Two-sided bounds without exploring the shells above

The objective asked for a way to **bound** `r_ε` without enumerating all graphs
of superior radius. Three bounds, each cheap and each provable.

**Bound L1 (free lower bounds).** `r_ε ≥ r_val` for all `ε ≥ 0` (Remark 3), and
`r_ε ≥ r_{ε'}` for `ε ≥ ε'` (Theorem C(6)). So the already-computed validity
radius, and any smaller-`ε` radius, are lower certificates at zero extra cost.

**Bound L2 (shell-probe lower bound).** Evaluating a **single** shell `d` by
Theorem C(4) — without touching any shell below or above it — gives
`β↑(d)`. If `β↑(d) ≤ ε` then `r_ε > d`, by monotonicity of `β↑`. A lower bound
is thus purchasable at the price of one shell.

**Bound U1 (anytime chain upper bound).** Any single monotone retraction chain
`G₀ ⋖ G₁ ⋖ … ⋖ G_m` with `B(G_m) > ε` certifies `r_ε ≤ m`, by Definition 5. A
greedy chain that at each step retracts the orientation maximising `B` yields
such a witness in `O(|K_{G₀}|²)` closures — no shell is enumerated at all. The
witness is replayable and checkable independently, exactly as
`sat.verify.verify_walk` does for E3.

**Bound U2 (the ∞ test).** `B(Ĉ) ≤ ε ⟹ r_ε = ∞`, by Theorem B (`Ĉ` is the
maximum). One evaluation of `B` at the top of the space settles, for every `ε`
above `B(Ĉ)`, that no perturbation whatsoever can push the bias past `ε`. This
is the bias analogue of the paper's "no failure anywhere" case and, like it,
is the case where an exhaustive search would be most expensive.

**Bracketing.** L1/L2 and U1 together return an interval `[lb, ub] ∋ r_ε`
after work that is linear in `|K_{G₀}|`, and the interval can be tightened
shell by shell in either direction. Because `β↑` is monotone (Theorem C(1)),
**binary search on `d` is valid**: `r_ε` is located in `⌈log₂|K_{G₀}|⌉` shell
evaluations, using Theorem C(4) to build each probed shell directly.

### The dispatch

Two exact strategies, with complementary costs — mirroring the paper's existing
hybrid:

* **Incremental** — evaluate shells `r_val, r_val+1, …` in order, stopping at
  the first crossing. Cost `Σ_{d ≤ r_ε} |shell d|`. Best when `r_ε` is small,
  which the paper's measurements say is the common case.
* **Bisection** — binary search on `d ∈ [r_val, |K_{G₀}|]`, each probe built
  directly by Theorem C(4). Cost `O(log |K_{G₀}|)` shells, but a probe near
  `|K_{G₀}|/2` is the largest shell. Best when `r_ε` is large or unknown.

The discriminator is free at runtime, as in `hybrid.py`: run incremental under a
shell budget; if it has not crossed, switch to bisection. Both are exact, so the
dispatch is purely a performance decision and the returned radius is identical.

---

## 6. Scope, and what each statement assumes

| # | Statement | Status | Rests on |
|---|---|---|---|
| A | `B ≡ 0` on the ball of radius `r_val − 1` | **PROVED** | definition of `r_val`; validity ⟹ identification in linear SCMs |
| B | `G ⪯ H ⟹ B(G) ≤ B(H)` | **PROVED** | definition of the order only |
| B1 | `F_ε` is upward-closed | **PROVED** | B |
| B2 | `F_{ε'} ⊆ F_ε` for `ε ≤ ε'` | **PROVED** | — |
| U | nearest state of **any** up-set is reachable by retraction | **PROVED** | Lemma L (`THEOREMS.md` §5) |
| R | `↑G₀ = {Meek(Ĉ,K₀\S)}`, `d = ρ(G₀) − ρ(G)` | **PROVED** | Lemma O, Lemma R, Property S |
| C | monotone staircase; sphere sufficiency; direct shells; one pass for all `ε` | **PROVED** | A, B, U, R |
| D | practitioner certificate, band `(ε_x, ε_y]` | **PROVED** | C, Meek's theorem |
| E | semi-local `T(G)`, no DAG enumeration | **PROVED** | Meek's theorem |

**Inherited assumptions.** Lemma L, Lemma R and Property S are taken from
`THEOREMS.md` §4b–§5 exactly as the paper takes them for `r_val`. This is
deliberate: **`r_ε` inherits the assumption chain of `r_val` and adds nothing to
it.** Any future strengthening or weakening of that chain moves both radii
together. In particular the one-sided error direction is unchanged — an upward
search can only return a radius that is too large, i.e. can only overstate
robustness, and every result object carries that caveat.

**A gap in the supplied Anti-Exchange proof, found and closed.** The
`THEOREMS.md` supplied for this session upgrades Anti-Exchange Case B from
*verified* to *proved* via Chickering's covered-edge-reversal theorem. Step 2 of
that argument asserts that the reversal sequence between two extensions
`D_A, D_B ∈ [G]` stays inside `[G]` on the grounds that it "operates exclusively
on uncompelled edges". That does not follow: Chickering's theorem gives a path
inside the Markov equivalence class `[Ĉ]`, and the orientations of `K_G` are
precisely edges **uncompelled in `[Ĉ]`**, so they are exactly what the theorem is
free to reverse. The classical substitute — reorienting freely inside a chordal
chain component — is also unavailable, because an MPDAG carrying background
knowledge can have a non-chordal undirected part, which is the phenomenon
`THEOREMS.md` §1 documents and which we re-confirmed computationally this
session.

The conclusion is nonetheless true, for a different reason, and `THEOREMS.md`
§24 now carries the corrected proof. What is needed is not the *existence* of a
Chickering path but its **length**: the strong form of the theorem produces
*exactly `m`* distinct covered-edge reversals, where `m` is the number of
differing edges. A walk of `m` unit steps taking the difference count from `m`
to `0` must decrease at every step, so every reversal flips a *currently
differing* edge and no edge on which `D_A` and `D_B` agree is ever touched. Since
both extend `G`, they agree on all of `dir(G)`, so every intermediate DAG
contains `dir(G)` and lies in `[G]`. Case B, and with it Anti-Exchange, Lemma R,
Property S, Lemma L and Conjecture 2, are therefore **unconditional**.

Two consequences for this document. First, the assumption chain above is
discharged: `r_val` and `r_ε` are both exact, not exact-modulo-a-conjecture, and
the one-sided caveat Sessions 3-7 attached to every radius can be dropped.
Second, nothing in the results here changes, because they were stated relative
to Lemma L as a hypothesis and that hypothesis is now a theorem.

**Beyond linear-Gaussian.** Definitions 1–2 need only that (a) the analyst's
estimand is a functional of `P` and `Z` alone, and (b) each `D ∈ [Ĉ]` assigns a
number `τ_D` identified from `P`. Both hold for the g-formula under any SCM with
positivity, with `τ_D = E[Y | do(X)]` contrasts computed by adjusting for
`pa_D(X)`. Theorems A, B, B1, B2, U, R, C, D and E go through verbatim in that
generality — none of their proofs uses linearity. Linearity is used only to make
`β(S;Σ)` a closed-form linear-algebra evaluation, i.e. it is a statement about
the *cost* of `B`, not about its *validity*. The empirical study of §7 is
linear-Gaussian because that is where exactness is cheap.

**Finite samples.** `Σ` is estimated. Write `B̂` for the plug-in built from
`Σ̂`. Since `T(G)` is a max over a family of regression coefficients indexed by
the finitely many possible parent sets, and each is a smooth function of `Σ` at
any `Σ ≻ 0`, the delta method gives `sup_G |B̂(G) − B(G)| = O_p(n^{-1/2})`
uniformly over the (finite) space. A conservative radius is obtained by
thresholding an **upper** confidence bound on `B` instead of `B̂`: inflating `B`
enlarges `F_ε`, which can only *decrease* `r_ε`, i.e. can only understate
robustness. §7 measures this.

---

## 7. What the empirical study checks

The claims above are theorems, so the study is a *falsification* harness, not
evidence for them; it is there to catch implementation error, and to measure the
quantities the theorems do not determine.

**Falsifiable predictions (any violation is a bug or a refutation).**

* **P1** `B(G) = 0` for every `G` with `d(G₀,G) < r_val`. (Theorem A)
* **P2** `B(G) ≤ B(H)` for every ordered pair `G ⪯ H`. (Theorem B)
* **P3** `β↑(d)` non-decreasing; `β↑(d) = max` over sphere `d`. (Theorem C 1,3)
* **P4** `r_ε` from the staircase equals `r_ε` from an independent brute-force
  BFS over the fully enumerated space with predicate `B > ε`. (Theorem C 4,5)
* **P5** `r_val ≤ r_0 ≤ r_{ε}` and `ε ↦ r_ε` non-decreasing. (Theorem C 6)
* **P6** For every `S ⊆ K_{G₀}`, the realised error `|θ_Z − τ_{D*}|` of the
  state `Meek(Ĉ,K₀\S)` that contains the truth is `≤ β↑(|S|)`. (Theorem D)
* **P7** The semi-local `T(G)` of Theorem E equals the enumeration-based one.

**Measured, not predicted.**

* The shape of the staircase past `r_val` — step or ramp. The paper's App. C
  table shows a ramp on the running example; whether that generalises is open.
* Conservativeness: the gap between `β↑(k)` and the realised error, the bias
  analogue of the paper's 0.43–0.47 coverage-conservativeness figure.
* The cost of `r_ε` relative to `r_val`, and the incremental-vs-bisection
  crossover.
* **Counterexamples to the naive readings**, which are the negative results this
  session owes the reader: (a) individual-state bias is not monotone in
  distance; (b) the incumbent `mean_abs_bias` is not monotone under `⪯`, so its
  `r_eps` is not computable by upward search and carries no certificate.
