# Session 8 — the ε-bias radius, and an unconditional `r_val`

Branch `experiments/r_epsilon_tests`. Everything below traces to a file under
`results/epsilon/`, `THEOREMS.md` or `docs/R_EPSILON_THEORY.md`. No number here
was written from memory.

---

## 0. The headline

Two results, one expected and one not.

**Expected.** `r_val` extends to a bias certificate. There is a worst-case bias
functional `B` that is *exactly zero* on the ball `r_val` already certifies,
monotone under knowledge retraction, computable without enumerating graphs, and
an upper bound on the error the analyst actually carries. The radius `r_ε` at
which it first exceeds a threshold is computed by the *same* machinery as
`r_val`, because the retraction reduction was never about validity — it is a
statement about arbitrary up-sets. Pairing two thresholds gives the band the
session asked for: *at most `ε_y`, and demonstrably able to exceed `ε_x`*.

**Not expected.** The Anti-Exchange proof supplied for this session has a false
premise, and the gap is closable. The supplied Case B argues that Chickering's
reversal path between two extensions of `G` stays inside `[G]` because it
"operates exclusively on uncompelled edges" — but the orientations of `K_G` are
*precisely* the edges uncompelled in `[Ĉ]`, so that is exactly what the theorem
is free to reverse. The classical repair (reorient freely inside a chordal chain
component) is unavailable, because an MPDAG carrying background knowledge can
have a non-chordal one — the phenomenon `THEOREMS.md` §1 already documents.

The conclusion is nonetheless true, for a different reason. The **strong** form
of Chickering's theorem produces *exactly `m`* reversals for `m` differing
edges. A walk of `m` unit steps taking the difference count from `m` to `0` must
decrease at every step, so every reversal flips a *currently differing* edge and
never touches one on which the two graphs agree. Two extensions of `G` agree on
all of `dir(G)`, so every intermediate retains it and stays in `[G]`.

**Consequence.** Anti-Exchange, Lemma R, Property S, Lemma L and Conjecture 2
are now unconditional. The one-sided caveat that Sessions 3–7 attached to every
reported radius can be dropped: `radius_local_up` is exact, not
exact-if-a-conjecture-holds. No previously reported number changes.

---

## 1. What was added

| | |
|---|---|
| `THEOREMS.md` | supplied file saved verbatim, then §§16–24 appended: the bias functional, Theorems A/B/U/R/C/D/D′/E, and §24 — the gap and its repair |
| `docs/R_EPSILON_THEORY.md` | the full statements and proofs, with a per-item status table and the scope of every verification |
| `src/bkrobust/epsilon/` | `bias.py`, `profile.py`, `certify.py`, `finite_sample.py`, `figures.py`, plus the three harnesses `verify.py`, `chickering.py`, `counterexamples.py` |
| `tests/epsilon/` | 26 tests, all passing |
| `experiments/` | `run_r_epsilon.py` (the sweep), `run_r_epsilon_worked.py`, `run_r_epsilon_finite_sample.py`, `run_c2_rate_mc_aware.py` |
| `paper/sections/epsilon-radius.tex` | the drafted section, no `TODO` numbers left |
| `figures/` | `fige1_staircase`, `fige2_what_epsilon_buys`, `fige3_cost` |

---

## 2. The construction, in one page

The analyst holds one observational law with covariance `Σ` and reports one
number, `θ_Z = β(Z;Σ)`, which does not depend on which DAG is the truth. Each
candidate truth `D` assigns the effect its own value `τ_D = β(pa_D(X);Σ)` — the
IDA quantity. So a knowledge state leaves a finite set of possible effects,
`T(G) = {τ_D : D ∈ [G]}`, and

```
B(G) = max { |θ_Z − τ| : τ ∈ T(G) } .
```

* **Theorem A.** `d(G₀,G) < r_val ⟹ T(G) = {θ_Z} ⟹ B(G) = 0` — exactly zero, at
  the population level, for every law compatible with `Ĉ`. Inside the certified
  ball there is nothing to recompute: every state gives the same estimand.
* **Theorem B.** `T` grows with `[G]`, so `B` is monotone under the order and
  `{B > ε}` is an up-set.
* **Theorem U.** The §6 retraction argument consumed *only* upward-closure of
  failure. It is therefore a theorem about arbitrary up-sets, and `r_ε` inherits
  it free.
* **Theorem C.** The staircase `β↑` is non-decreasing, zero below `r_val`,
  attained on the sphere, and each shell is generated directly from the
  `d`-subsets of `K_{G₀}` with no shell below it. One traversal answers every
  `ε`; monotonicity licenses bisection.
* **Theorem D.** If at most `k` of `K_{G₀}` are false of the truth, the realised
  error is `≤ β↑(k)`. Hence the band.
* **Theorem D′.** The band is *attained* by a world the analyst cannot rule out,
  so it cannot be tightened without a further assumption.
* **Theorem E.** `T(G)` is read off the possible parent sets of `X`, obtained by
  orienting the undirected edges at `X` each way and Meek-closing —
  `2^{deg(X)}` polynomial closures, no DAG enumeration.

**`r_ε` adds no assumption to `r_val`.** It rests on Lemma L, Lemma R, Property
S and Lemma O exactly as the existing radius does — and, after §24, those are
unconditional.

---

## 3. Verification

### 3.1 The repair (`results/epsilon/chickering/`) — exhaustive at n ≤ 5

| check | objects | violations |
|---|---|---|
| L1 — a covered reversal preserves the equivalence class | 50,758 | **0** |
| L2 — a difference-reducing sequence of exactly `m` reversals exists | 206,878 | **0** |
| L3 — every intermediate DAG stays inside `[G]` | 2,907,242 | **0** |
| L4 — `[G]` is connected under reversals that stay in `[G]` | 118,643 | **0** |
| AE-B (semantic) — no two undirected edges induce the same bipartition of `[G]` | 293,814 | **0** |

Scope: **all** DAGs, CPDAGs and knowledge states at n = 3, 4, 5 — all 8,782
CPDAGs and all 116,992 knowledge states at n = 5, with no sampling anywhere;
plus an n = 6 spot check for L1. L2 and L3 are the literal computational content
of the two repair lemmas; AE-B is an independent check that does not pass
through Chickering at all. **3,584,369 checks, 0 violations.** Verification is
not proof, and §24 carries the proof; the sweep is there because a proof with a
slip in it and a harness with a bug in it rarely agree to 3.5 million cases.

### 3.2 The construction (`results/epsilon/verification/`)

400 screened instances × 3 independent parameter draws. Tolerance `1e-9`.

| | checked | violations |
|---|---|---|
| P1 zero bias inside the shell | 1,758 | 0 |
| P2 monotone over ordered pairs | 372,099 | 0 |
| P3 staircase / sphere sufficiency | 17,160 | 0 |
| P4 incremental = bisection = brute-force BFS | 5,844 | 0 |
| P5 nesting in `ε`, and `≥ r_val` | 5,844 | 0 |
| P6 certificate audit | 22,680 | 0 |
| P7 semi-local = enumerated parent sets | 11,413 | 0 |

**436,798 checks, 0 violations.**

---

## 4. The simulation study (`results/epsilon/study/`)

900 accepted instances; five generators, n = 6–12, three knowledge coverages,
four corruption regimes.

**What `ε` buys is a ceiling, not a longer runway.** `r_ε > r_val` on 1.0% of
instances at `ε = 1%`, rising to 24.5% at 25% and 75.6% at 100% — and when it
does, the gain is almost always exactly one shell. What rises far faster is the
fraction where the threshold is *never reached*: 1.4% → **70.1%** at `ε = 100%`.
On seven instances in ten, no perturbation of the analyst's knowledge, however
large, can move the estimate by its own magnitude. `r_val` cannot say that at
all, and it comes from the same traversal.

**`r_0 = r_val` on all 900 instances, 0 exceptions** — the converse to Theorem A
fails only non-generically, and the non-generic case was never hit.

**The staircase is a step 84.1% of the time**, a ramp 15.1%, flat 0.8%. The
gradual rise the paper's Appendix C shows is, in the large majority of
instances, an artefact of averaging over a shell whose contaminated fraction is
growing — not of the worst case rising.

**The certificate held on every one of the 884 instances where `β↑(k)` was
actually computed — 0 violations.** Sixteen further instances had the traversal
stop before depth `k`, so the recorded `β↑(k)` is a lower bound and the
comparison is not a test of Theorem D; three of those had a realised error above
that lower bound and were originally reported as violations. Recomputing their
full staircases, the certificate holds in all three — **with exact equality**,
which is Theorem D′ showing up. Among instances with a non-degenerate bound, the
conservativeness ratio has median `1.000` and maximum `1.000`: the bound is
routinely attained, not merely valid.

**Cost.** Both exact strategies agree on every instance. Incremental evaluates a
median of 1 state (mean 2.09, max 65); bisection a median of 1 (mean 5.14, max
172). Bisection is the wrong default *for this regime* — `r_ε` is small almost
everywhere, so a mid-lattice probe is wasted — which is why the dispatch tries
incremental first. The anytime greedy chain, which enumerates no shell at all,
returns the **exact** radius on 98.4% of instances.

**Against the incumbent.** On a matched subsample the new radius disagrees with
the mean-based `_radius_eps` on at least one threshold in 84% of instances, in
both directions. On the same instances the mean was non-monotone under the order
in **49 of 50** cases, so the upward search used to compute its radius is not
licensed there at all.

---

## 5. The negative results (`results/epsilon/counterexamples/`)

Both naive readings of "bias grows past `r_val`" are false, and both refutations
are minimal and replayable (`replay.py` re-derives every headline number from
the stored witness).

**C1 — bias is not monotone in distance.** At four nodes: a state at distance 1
with `B = 1.3405`, and a state at distance 2 with `B = 0.0000`. They are
incomparable in the order, so Theorem B is untouched; what fails is the reading
that confuses distance with knowledge. This is why `β↑` maximises over the shell
rather than following a path.

**C2 — the incumbent mean is not monotone, and its radius is wrong.** The
flagship witness has `[G] ⊊ [H]` with mean bias `0.96–1.03` at `G` against
`0.43–0.46` at `H`, reproduced in 8 of 8 independent repetitions at a minimum of
**11.4 combined standard errors**. Rate, scored on identical pairs with a
Monte-Carlo-aware rule (3 combined s.e. in each of two repetitions, 300 draws —
so a *lower bound*): **34.7%** of comparable pairs for the mean, **0.0%** for
`B`. Consequence: the retraction-only search disagrees with brute force on
**269 of 644** threshold probes, and in the recorded cases by reporting that no
perturbation reaches the threshold anywhere when one does at distance two — an
error in the direction that overstates robustness.

The first pass at that rate used an `1e-6` threshold at 80 draws, which sits far
below the Monte-Carlo error and conflated "the estimand is non-monotone" with
"this estimate of it was". `experiments/run_c2_rate_mc_aware.py` separates them;
the `34.7%` above is the separated number.

---

## 6. Finite samples (`results/epsilon/finite_sample/`)

`r_val` reads no numbers and has no sampling distribution. `r_ε` reads a
covariance and does, and a plug-in `r_ε` reported as exact is a confidence
statement with the confidence removed. Over 120 instances, n = 100 … 10,000:

* the **plug-in** radius is too large — the direction that overstates
  robustness — in up to **5.0%** of cases at n = 100, falling to 0–0.8% at
  n = 10,000;
* thresholding a bootstrap upper bound instead, and reporting the 5% quantile of
  the radius, gives **0.0% in every cell**, at the cost of understating the
  radius on 1–17%.

The conservative radius is the one to report.

---

## 7. What is not established

* **Verification is not proof.** §24's repair is a proof; the n ≤ 5 sweep
  corroborates it and does not extend it. L1's n = 6 check is a spot check only.
* **The `34.7%` non-monotonicity rate is a lower bound** by construction. The
  true population rate is higher; how much higher is not established.
* **Linear-Gaussian is where the *cost* claims hold.** Theorems A–E need only
  that the estimand is a functional of `P` and `Z`, and that each DAG identifies
  a number from `P`; the g-formula under positivity satisfies both. Linearity
  makes `B` a closed-form linear-algebra evaluation, nothing more. Nothing
  beyond linear-Gaussian was measured.
* **Scope of the sweep**: n ≤ 12, `|undirected| ≤ 10`, synthetic generators
  only. The real-graph corpus was not touched; `r_ε` on the 85-vertex component
  is untested.
* **`Z` is the optimal adjustment set throughout.** The behaviour of `r_ε`
  across the other valid sets, which the paper reports for `r_val`, was not
  measured.
* The **staircase-shape split** (84% step) is a property of this corpus, not a
  theorem. Nothing forbids a long ramp.

---

## 8. One thing to decide

`hybrid.py`'s `HybridResult.assumes` still reads *"Conjecture 2 (hence
Anti-Exchange Case B, verified not proved)"*, and the module docstring says the
search is exact *iff* Conjecture 2 holds. After §24 that is no longer the
caveat the result carries. Updating it changes no computed number, but it does
change the string attached to every radius this repository has produced, and the
paper's §A.3 and the reproducibility statement would have to move with it. Left
unchanged pending that decision.
