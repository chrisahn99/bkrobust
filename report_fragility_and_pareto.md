# Robustness survival and the Pareto frontier of causal assumptions

**Branch:** `experiments/survival-and-pareto`
**Pre-registration:** `results/axis_robustness/PREREGISTRATION.md` (§1–§9 written before any
survival curve or rank correlation existed; Appendices A–F record every falsification,
correction and retraction, in order)
**Gate:** `benchmarks.measure.fast_gate` exclusively, recorded in every row. `synth.runner.gate`
was not used — it calls the exponential `all_valid_adjustment_sets_mpdag` internally. `fast_gate`
is strictly stricter (16/46,800 disagreements, one direction), which biases **against** the
hypothesis under test.
**Assumption carried on every radius:** *Conjecture 2 (hence Anti-Exchange Case B, verified not
proved).* All upward searches are exact iff it holds; the error is one-sided, so radii can only
be too **large**.

---

## 1. What was asked, and what the data actually supports

The brief asked for two demonstrations: that `r_val` ranks average-case structural survival
where SHD and `|K|` fail, and that precision trades off against robustness along a Pareto
frontier ("the fragility of precision").

**The first is supported.** **The second is false**, and the reason it is false is the more
interesting result of the two.

| Prediction | Verdict | Evidence |
|---|---|---|
| **P1** `r_val` positively ranks survival | **SUPPORTED** | τ_b > 0 in 9/9 strata; 8/9 under the conservative endpoint |
| **P2** `r_val` outranks `\|K\|` | **SUPPORTED** | 9/9 strata; 7/9 under the conservative endpoint |
| **P3** `\|K\|` is a null predictor (project hypothesis **H4**) | **NOT SUPPORTED** | CI covers 0 in only 3/9; `\|K\|` is *negatively* associated |
| **P4** tiered AUC variance > flip | **NOT COMPARABLE** | unit mismatch; no verdict forced (Appendix E.4) |
| **P5** claim-axis and hop-axis agree in sign | **PASS** | 6/6 strata, independent subsample |
| **P6** a non-trivial efficiency/robustness frontier exists | **FALSIFIED** | utility exactly invariant in 32/32 strata |
| **P6′** frontier exists for *achievable* efficiency | **FALSIFIED** | invariant again; Phase 2 abandoned per pre-registered trigger |
| **P7 / P7′** aggressive BK trades utility for radius | **UNTESTABLE** | one axis is constant; census reported instead |

---

## 2. Phase 1 — survival under corrupted background knowledge

### 2.1 Design

An analyst holds CPDAG `Ĉ`, asserts knowledge `K`, Meek-closes to `G₀`, and reads
`Z* = optimal_adjustment_set_mpdag(G₀, x, y)` **once, at depth 0, then holds it fixed**. `K` is
then corrupted and we ask whether the set they already committed to survives —
`is_gac_valid_mpdag(G, x, y, Z*)`, polynomial, no enumeration of `[G]`.

1,978 instances; 3,136,000 sampled states; 200 reps per depth; 1,656 s. Instances span
component size 6–12, three separations, coverage ∈ {0.5, 1.0}, and a **base wrongness**
`b ∈ {0, 0.10, 0.25}` applied before the sweep so that `SHD(G₀, truth)` is not constant — at
coverage 1.0 with truthful `K`, `G₀` *is* the ground-truth DAG and SHD ≡ 0, which would have
made the baseline comparison rigged by construction.

Two arms, never merged: **flip** (uniform, independent, targeted depth) and **tiered**
(block-correlated, re-binned on `corruption_rate`; see §2.5).

### 2.2 P1 — `r_val` ranks survival

τ_b(AUC, `r_val`), 10,000-resample bootstrap over instances, seed 0:

| stratum | τ_b | 95% CI | n |
|---|---|---|---|
| flip cov=0.5 bw=0.00 | +0.714 | [+0.625, +0.798] | 240 |
| flip cov=1.0 bw=0.00 | +0.519 | [+0.442, +0.591] | 240 |
| flip cov=1.0 bw=0.10 | +0.388 | [+0.296, +0.478] | 236 |
| flip cov=0.5 bw=0.10 | +0.380 | [+0.265, +0.488] | 240 |
| flip cov=1.0 bw=0.25 | +0.346 | [+0.209, +0.472] | 134 |
| flip cov=0.5 bw=0.25 | +0.031 | [−0.114, +0.169] | 220 |
| tiered n_tiers=3 | +0.811 | [+0.787, +0.833] | 238 |
| tiered n_tiers=2 | +0.799 | [+0.761, +0.834] | 190 |
| tiered n_tiers=4 | +0.777 | [+0.753, +0.798] | 240 |

Positive in 9/9; CI excludes zero in 8/9. The exception is the highest-corruption,
lowest-coverage flip cell, where the association is indistinguishable from zero — reported, not
dropped.

**Effect sizes** move monotonically: for tiered `n_tiers=3`, median AUC rises 0.67 → 0.94 → 0.99
→ 0.98 → 1.00 across `r_val` buckets 1→5+.

### 2.3 P2/P3 — the baselines fail, but not the way the brief assumed

`r_val` outranks `|K|` in **9/9** strata. But `|K|` is **not inert**, so H4's framing is wrong:
its association with survival is real and **negative** — flip −0.194/−0.205 at high corruption,
tiered −0.425 to −0.470, all with CIs excluding zero.

> More asserted claims means *worse* survival, not better.

`SHD(G₀, truth)` fails for a different reason again: it is `undefined (predictor constant)` at
`b = 0` **by construction**, and elsewhere its sign flips across strata (+0.240, +0.027, −0.175,
−0.171, −0.134). It has no consistent direction. That is a fairer indictment than "SHD tangles":
a third of the design cannot score it at all.

### 2.4 The discordance case

**`c06_s01_cov050_seed00000007_bw000`.** `shd_truth = 0` — `G₀` *is* the ground-truth DAG, the
analyst's knowledge is perfectly correct — with `n_k = 2`. Yet `AUC_frac = 0.0`: `S(1) = 0.0` on
107 non-contradictory samples, `S(2) = 0.0` on 200. Every corruption invalidates the committed
set, and `r_val = 1` says exactly that.

**Zero SHD is compatible with maximal fragility.** SHD and `|K|` describe how much was asserted;
`r_val` describes how much can go wrong before the answer breaks.

A second, opposite spotlight (`AUC_frac = 1.0` at `r_val = 8`) was **withdrawn on verification**:
its AUC rested on depths with `n_eval` of 35, 10, 1 and 1. Four surviving draws are not a survival
curve.

### 2.5 Two defects found and corrected mid-analysis

**The tiered depth metric was broken (Appendix E).** It counted assertions added or reversed but
not *removed*; `tiered` relocates nodes, and relocation into a shared tier silently deletes a
cross-tier assertion. Its `d = 0` bin was contaminated by every corruption rate and contained
genuine failures. Re-binning on `corruption_rate` — the knob actually turned — gives a clean
design: contradiction rises monotonically 0.000 → 0.852, `S` falls 1.000 → 0.618, `S = 1.000`
exactly at rate 0 for all 668 instances individually, n = 133,600 per rate. Tiered τ values
**survive** the correction (0.777–0.811). The flip arm was never affected.

**Effective sample size varies with `r_val` (Appendix F).** τ(`r_val`, usable depths) runs from
+0.318 to −0.438 depending on stratum, so the bias runs *against* the hypothesis in all
coverage-1.0 and tiered strata and *for* it in the two strata carrying the largest τ. Under the
pre-specified control `AUC_frac_usable` (depths with `n_eval ≥ 30`), **P1 falls 9/9 → 8/9 and P2
9/9 → 7/9**. Weakened, not overturned. Per Appendix F.3 the conservative figure is the one quoted
in any summary claim.

### 2.6 P5 — the units question

`r_val` counts BFS hops; corruption counts claims. These are not the same axis, and
`radius_local_up_fast` cannot measure the distance because it searches *upward* only while a
flip-corrupted state is not a refinement of `G₀`. Exact hop distance therefore required full
space enumeration on a small-component subsample (1,256 instances, component size 3–7; spaces of
32–607 elements, 0.007–12.2 s to build).

One corrupted claim ≈ **2.9 hops** (sd 1.01), converging to exactly 2·d at d ≥ 5. **P5 passes in
6/6 strata**, signs agreeing, τ positive on both axes. Cross-checked against the main dataset:
6/6 sign agreement, 4/6 rough magnitude agreement; the two divergent cells are plausibly a
component-size effect (main spans 6–12, subsample 3–7) — reported, not resolved.

### 2.7 The finding that was not predicted

Contradiction rates by corruption depth, flip arm:

| depth | 1 | 2 | 3 | 5 | 8 | 12 |
|---|---|---|---|---|---|---|
| rate | 0.491 | 0.773 | 0.896 | 0.930 | 0.986 | 1.000 |

At `d = 1` and base wrongness 0.00/0.10/0.25 the rates are 0.653/0.647/0.563.

> **Most orientation errors are self-revealing.** Reversing a single truthful claim renders the
> knowledge set inconsistent with the CPDAG about 6 times in 10. Meek closure fails and the
> analyst finds out.

And among corrupted sets that *stay* consistent, median survival is **0.952** — most survivors
are harmless. The picture is a two-stage filter: most errors announce themselves, most of the
rest are benign, and the dangerous case is **consistent *and* invalidating**. That is precisely
the case `r_val` is defined over, since its space contains only consistent MPDAGs.

This also dissolves an apparent tension. `r_val = 1` is common while average survival at `d = 1`
is high. Both are true: `r_val` is a **worst-case** distance to the nearest failing graph, `S(d)`
an **average-case** probability over a shell. A failure at distance 1 says nothing about how much
of the shell fails. That the worst case predicts the average case is the finding, not a tautology.

---

## 3. Phase 2 — there is no Pareto frontier

### 3.1 Falsified twice

Round 1 (24 base CPDAGs, 1,087 truthful proposals, 17,680 SEM draws): strata with more than one
distinct `utility_median`: **0 of 24**, while `r_val` varied in 18. Round 2 redefined utility as
*achievable* efficiency over a capped 12-set menu: **0 of 32** strata varied in achievable
utility, in the count of valid menu members, or in the winning member — while `r_val` varied in
23. The pre-registration's trigger was explicit, so Phase 2 was **abandoned rather than
redefined a third time**.

### 3.2 Why — and this is the substantive result

> Every truthful proposal leaves the true DAG inside `[G₀]`. `optimal_adjustment_set_mpdag`
> returns non-`None` exactly when all extensions agree; if they agree, the common answer *is* the
> true DAG's optimal set. Every identifying truthful proposal therefore yields the same `Z`, hence
> the same asymptotic variance.

**For truthful background knowledge, statistical efficiency is an invariant of `(Ĉ, x, y)` — not
a function of what the analyst asserts. Knowledge moves the robustness axis and only the
robustness axis.**

So the trade-off is not precision against robustness. Precision is *fixed*. What knowledge
purchases is **identifiability** — 385 of 1,431 proposals (26.9%) could justify no valid set at
all. The operative trade-off is **identifiability against robustness, and it is a step, not a
frontier**: assert enough to identify, and every further claim is pure `r_val` cost at zero
precision gain.

This converges with Phase 1's P3 from the opposite direction: there, more claims measurably
*reduced* survival. Two independent designs, same conclusion — marginal truthful background
knowledge is not free, and its price is paid entirely in robustness.

Testing the step hypothesis properly — does `r_val` fall monotonically with proposal size beyond
the minimum identifying set? — is a separate question and is **not claimed here**.

### 3.3 P7′ — the census stands even though the contrast cannot be run

Round 1 found **zero** "aggressive" proposals, an artifact of `component_generator`, which builds
the undirected component purely from `X`'s upstream confounders and keeps `X→Y` and its
descendants outside it by construction.

On real structure the picture is different. Sampling ≤200 descendant-ordered pairs per network
across the 39-network corpus (M-bias excluded as an ADMG), **28 networks admit at least one
(x, y) pair with an undirected CPDAG edge on a directed causal path — 396 pairs**, densest in
`paths` (96), `Polzer_2012` (64), `Kampen_2014` (33), `child` (22), `sachs` (22), `insurance`
(21). Eight real-network bases were gated and yielded 215 aggressive and 129 conservative
proposals, so the empty cell is genuinely filled — but utility is invariant there too (0 of 8
strata vary), so the intended contrast cannot be run. The residual `r_val`-only contrast is
non-directional (2 bases lower for aggressive, 2 tied, 1 reversed) and is **inconclusive**.

---

## 4. What is not claimed

- **No cross-arm detectability comparison.** An earlier draft claimed correlated errors are
  harder to detect than uniform ones (0.354 vs 0.653 at `d = 1`). The tiered half came from the
  contaminated bin, and the comparison was never sound in principle: one arm relocates a fraction
  of *nodes*, the other reverses a count of *claims*. Retracted in Appendix E.3. Establishing it
  needs a corruption process parameterised identically across both arms — a design change, not a
  reanalysis.
- **No claim of perfect ranking.** The endpoint is a rank correlation with a bootstrap interval.
  One stratum is null.
- **Synthetic structure only** for Phase 1. Real networks appear only in the P7′ census.
- **Oracle CI throughout.** No finite-sample discovery; causal sufficiency assumed. These remain
  the project's largest scope limits and this work does not touch them.
- **`AUC_frac` is not quoted alone.** The conservative `AUC_frac_usable` figure governs summary
  claims (Appendix F.3).

## 5. Provenance

All results under `results/axis_robustness/`: `survival_instances.csv` (1,978),
`survival_curves.csv` (14,066), `survival_samples.csv.gz` (3,136,000 rows, 48 MB compressed),
`hops_samples.csv.gz` (6 MB), `analysis_tau_primary.csv`,
`analysis_tau_secondary.csv`, `analysis_tiered_rebinned.csv` (668), `analysis_decomposition.csv`,
`analysis_discordance.csv`, `analysis_effect_sizes.csv`, `pareto_proposals.csv` (1,087),
`pareto_sems.csv` (17,680), `pareto2_proposals.csv` (1,431), `pareto2_menus.csv` (243),
`pareto2_aggressive_census.csv` (39), `hops_p5.csv`, `hops_instances.csv`,
`hops_cost_census.csv`, plus manifests capturing git SHA, interpreter, platform, library
versions and `RADIUS_CONVENTION`.

`num_workers: 1`, `random_seed: 0`, CP-SAT single-threaded, no global RNG anywhere in the live
tree. Determinism verified across `PYTHONHASHSEED` 0 / 12345 for every stage; bootstrap CIs
bit-identical across independent runs. The decomposition identity
`S_contra_as_fail(d) = (1 − contradiction_rate(d)) · S(d)` holds to 2.2e-16 on all 8,841 rows
where `S` is defined.

**Known incident.** The first survival sweep was destroyed when its worker relaunched and
overwrote its own outputs in place; it was regenerated (run 2) and reproduces the original's
pooled statistics exactly. The relaunch also fixed a real defect — run 1's `instance_id` was not
unique — so several Appendix C statistics computed by grouping on it were wrong and are corrected
in Appendix D, including a retracted endpoint swap. Full account in Appendices C–D.
