# ρ-breakdown-knowledge — PILOT RESULTS

**Idea #1 of the latent-causal campaign** (`output/2026-07-17_latent-causal-chris/reports/03-RANKED-IDEAS.md`)
*"How Wrong Would the Expert Have to Be? A Breakdown Radius for Background Knowledge in Causal Adjustment"*

**Date:** 2026-07-17 · **Compute:** betelgeuse, 16 CPU cores, ~20 min wall total (GPU unused)
**Scale:** 600 primary SCMs + 1 650 sweep SCMs + 60 nonlinear SCMs; 7 688 brute-force-validated Meek closures
**Verdict: SUPPORTED. The falsifier did not fire.**

---

## 1. The question

Background knowledge `K` about edge orientations cannot be tested against the observational
distribution (BRIEF item 10: `testable(K) ⊥ value(K)`; power = size at every *n*). The response
to an untestable assumption is to **price** it, not test it — Rosenbaum's Γ prices unmeasured
confounding; **ρ prices wrong orientations**.

> Does a **breakdown radius** for background knowledge exist, and is it informative?

Two numbers decide it, and neither was entailed:

1. **What fraction of single-orientation errors change `O*` at all?**
2. **Is the ignorance interval finite and informative at ρ = 1, 2?**

**Pre-registered falsifier:** *unbounded interval at ρ = 1* **AND** *near-universal `O*` disruption*
⇒ breakdown radius is 0 ⇒ **no paper**.

## 2. Setup

- **Graphs.** Random DAGs, `p ∈ {5,…,8}`, ER over a random topological order, expected degree
  ∈ {1.5, 2.0, 2.5}. CPDAG `C` = skeleton + v-structures + Meek closure.
- **Parameterisation.** **iSCM** (Ormaniec et al., **arXiv:2406.11601**) — each variable
  internally standardized during ancestral sampling against an `N = 100 000` reference sample; the
  standardization constants `(μⱼ, sdⱼ)` are then **frozen into the SCM**. For the linear arm the
  frozen SCM is again a linear SEM with coefficients `wᵢⱼ/sdⱼ` and noise variance `σⱼ²/sdⱼ²`, so
  **Σ, the true total effect, and every adjusted estimate are available in closed form**. The
  primary arm therefore carries **zero Monte-Carlo error**: each number is exact for its SCM.
- **Estimand / estimator.** True total effect `τ = ((I−A)⁻¹)[X,Y]`. The estimate under adjustment
  set `Z` is the **population** OLS coefficient of `X` in `Y ~ X + Z`, from Σ. Valid `Z` ⇒ exactly
  `τ` (verified to 5e-15, test T5).
- **Knowledge `K`.** Orientations of undirected edges of `C`, **true by construction** (hence
  Meek-consistent). `|K| ≤ 4`. Two regimes:
  - **generic** — a random subset of `C`'s undirected edges;
  - **tiered** — a random tiering consistent with `D*`'s topological order (Bang & Didelez,
    **arXiv:2306.01638**); `K` = the cross-tier orientations.
- **The ρ-ball.** All `C(|K|, ρ)` ways the expert could have got **exactly ρ of their claims
  backwards**, ρ = 1, 2, 3. Each perturbed `K'` is checked for Meek-consistency (Perković et al.
  Algorithm 1); consistent ones give an MPDAG `G_{K'}`, from which `O*(X,Y,G_{K'}) =
  pa(cn(X,Y,G)) \ forb(X,Y,G)` (**Henckel/Perković/Maathuis, arXiv:1907.02435**) and an estimate.
- **Ignorance interval at ρ** = [min, max] of those estimates over all consistent, amenable members
  at distance ≤ ρ (including the ρ = 0 point, which is `τ`).
- **Scale.** 600 linear SCMs (579 with a generic `K`, 180 with a tiered `K`), master seed
  `20260717`; per-SCM seeds `0…1799` recorded in `results/linear_raw.json`.

### 2.1 The machinery is tested, not trusted

`code/test_machinery.py` — **all pass**, locally and on betelgeuse:

| test | what it proves |
|---|---|
| T0 | Meek R1–R4 fire on the four textbook configurations; collider/chain CPDAGs correct |
| T1 | **CPDAG == the edges common to every DAG in the MEC**, by brute-force enumeration, 200 graphs |
| T2 | **Meek-consistency verdict AND the MPDAG match brute force on 7 688 knowledge sets** — the closure is validated against the *definition* (common orientations over all MEC members agreeing with `K`), not against a re-implementation. This is what a wrong R4 would have failed. |
| T3 | `O*` under **true** `K` is a valid adjustment set in `D*` (187 amenable cases, 0 failures) |
| T4 | **Monte-Carlo check of the analytic quantities**: Σ, the population OLS coefficient, and `τ` all match a 200 000-sample simulation of the frozen iSCM (max err 0.009, 0.003, 0.002) |
| T5 | every valid `Z` reproduces `τ` exactly (545 sets, max err 5.3e-15) |

### 2.2 Benchmark-leak diagnostics (mandatory, from run one)

Reisach/**Tami** var-sortability (**arXiv:2303.18211**), over all 600 SCMs, iSCM vs the naive
linear-SEM sampler on the *same graphs and weights*:

| diagnostic | naive SEM (control) | **iSCM (used)** |
|---|---|---|
| var-sortability | **0.937** [0.929, 0.945] | **0.508** [0.488, 0.530] |
| R²-sortability | 0.567 [0.546, 0.590] | 0.388 [0.365, 0.408] |

The control reproduces the known pathology (var-sortability ≈ 0.94 — `sortnregress` would read the
causal order straight off the marginal variances). **iSCM removes it (0.508 ≈ chance).** Honest
caveat: iSCM does **not** neutralise R²-sortability, it mildly **inverts** it (0.388 < 0.5) — a
known property of the parameterisation. Since nothing in this pilot *learns* a graph from data, no
result here is exposed to either leak; the diagnostics are reported because the campaign's
downstream claims will be.

---

## 3. Results — the two numbers

### NUMBER 1 — the fate of ONE wrong orientation (generic `K`, 2 058 ρ=1 perturbations)

| outcome | fraction of **all** single errors | fraction of **Meek-consistent** ones |
|---|---|---|
| **Meek-INCONSISTENT** — caught for free by the graphical check | **33.0 %** | — |
| `O*` **unchanged** — zero bias | 42.3 % | 63.1 % |
| **not amenable** — no valid adjustment set; the analyst gets a **loud failure**, not a wrong number | 15.3 % | 22.9 % |
| `O*` **changed and INVALID** — **silent bias** | **9.4 %** | **14.1 %** |
| `O*` changed but still valid | **0.0 %** (0 / 754) | 0.0 % |

> **NUMBER 1 = 36.9 %** of Meek-consistent single-orientation errors change `O*` at all
> (tiered: 41.0 %). **Only 14.1 % produce silent bias** (tiered: 17.0 %).
> Per-SCM mean 0.146 [0.124, 0.168] (generic), 0.175 [0.137, 0.218] (tiered).

**Near-universal disruption was required for the falsifier. It is not there.** Two thirds of
consistent single errors leave the optimal adjustment set *untouched*, and a third of all single
errors never even reach the analyst — Meek's own consistency check rejects them.

**This is the candidate explanation of b-LOAD's own unexplained empirical robustness under
"moderate structural noise"** (arXiv:2607.04447; BRIEF item 11 — *"Zé & Chris's own paper's hole
to close"*): most orientation noise never touches `pa(cn(X,Y)) \ forb(X,Y)`, because `O*` is a
*local* functional of the graph and most edges are not local to `(X,Y)`.

### NUMBER 2 — the ignorance interval (generic `K`)

| ρ | width/\|τ\| median | width/\|τ\| mean | q90 | interval **excludes 0** | ≥1 non-identified member |
|---|---|---|---|---|---|
| 1 | **0.000** | 0.351 [0.269, 0.439] | 0.869 | **94.3 %** [92.4, 96.2] | 50.9 % |
| 2 | 0.000 | 0.480 [0.376, 0.586] | 1.255 | 91.5 % [89.3, 93.8] | 67.2 % |
| 3 | 0.000 | 0.514 [0.419, 0.623] | 1.495 | 90.3 % [87.9, 92.7] | 69.4 % |

> **NUMBER 2 = yes.** The interval is **finite** (it is a max/min over a finite enumerated ball),
> and it is **informative**: for **73 % of SCMs its width at ρ = 1 is exactly 0** — one wrong edge
> changes the reported effect *not at all* — and it **excludes zero for 94.3 % of SCMs**. The sign
> of the effect survives a single expert error in 19 of 20 problems.

### The breakdown radius ρ\* itself

Smallest ρ at which *some* Meek-consistent member overturns the conclusion:

| criterion | ρ\* = 1 | ρ\* ≤ 2 | **never overturned in the ball** |
|---|---|---|---|
| **sign of τ reversed** | 5.7 % [4.0, 7.8] | 8.5 % | **90.3 %** [87.7, 92.6] |
| \|bias\| > 50 % of \|τ\| | 16.8 % | 23.3 % | 75.5 % |
| effect **not identified** (some member non-amenable) | 50.9 % | 67.2 % | 30.6 % |

**A breakdown radius exists, it is > 0 for ~90 % of problems, and it is reportable.**

### The finding that decides §3 of the paper: **the theorem cannot be a continuity bound**

Conditional on `O*` moving at all (ρ = 1, n = 274):

| \|bias\|/\|τ\| | median | mean | q90 | **max** |
|---|---|---|---|---|
| | **0.591** | 1.167 | 3.015 | **13.137** |

BRIEF item 11 warned that *"δ small ⇒ bias small" is FALSE in general* and that bias under an
invalid set is unbounded. **Confirmed, sharply.** When the damage lands, it is not small — the
median hit is 59 % of the true effect and the tail reaches 13×.

**So the breakdown radius survives on *probability*, not on *magnitude*.** ρ = 1 is safe because
the perturbation usually misses `O*` entirely, **not** because a near-miss is nearly right. This
**vindicates horn (a)** (a breakdown-radius result) and **kills horn (b)-as-continuity** — exactly
the §3 decision the ranked idea said the pilot must make. It also means the paper's theorem must
be **combinatorial** (which perturbations reach `pa(cn) \ forb`), never a Lipschitz-style bound.

### The consistency–correctness gap, quantified

`100.0 %` of Meek-consistent perturbed `K'` expel the true DAG from the MPDAG — **but this number
is entailed by the design, not a finding**: `K` is true, so reversing any of its statements makes
`K'` false, and no MPDAG built from a false orientation can contain `D*`. **The non-trivial
quantity is how often Meek's check catches it:**

| | Meek CATCHES the false `K'` | passes silently |
|---|---|---|
| generic `K` | **33.0 %** | 67.0 % |
| tiered `K` | 23.4 % | 76.6 % |

**Two thirds of false knowledge sails through the consistency check** — Perković et al.'s
*"consistent"* is a graphical check that never verifies truth, and the gap is large. **But** of
that 67 %, only ~14 % goes on to bias the estimate; the rest is absorbed by `O*`'s locality (63 %)
or converted into a visible non-identification (23 %). *The gap is real and wide; the damage that
crosses it is narrow.* That asymmetry is the paper.

### A structural finding: **knowledge is dangerous exactly where it is useful**

| stratum | n | ρ=1 width/\|τ\| mean | ρ=1 frac `O*` invalid |
|---|---|---|---|
| CPDAG **already amenable** rel. (X,Y) — `K` is only needed for *efficiency* | 75 (12.5 %) | **0.000 [0.000, 0.000]** | **0.000 [0.000, 0.000]** |
| CPDAG **not** amenable — `K` is load-bearing for *identification* | 504 | 0.403 [0.303, 0.505] | 0.209 [0.179, 0.237] |

**Every unit of damage, in all 600 SCMs, lives in the stratum where the knowledge was doing
identification work.** Where the CPDAG alone already identifies the effect, wrong knowledge was
harmless in **75/75** SCMs — sharper: **0 of the 430 Meek-consistent perturbations across those
75 SCMs, at any ρ, moved the estimate by more than 1e-9.** This is a third duality alongside item
10's utility–testability one:
**utility ⊥ safety**. It is a conjecture-shaped observation (see Limitations) and an obvious
Theorem-1 candidate for the paper.

### Tiered vs generic: **the enabling lemma did not help here — reported honestly**

The hypothesis under test was that **tiered / R1-closed knowledge controls the Meek cascade**
(Bang & Didelez, arXiv:2306.01638). In this design it did the opposite:

| | cascade (extra orientations beyond `K`) | frac of flips that cascade at all | `O*` disrupted | `O*` invalid |
|---|---|---|---|---|
| generic `K` | 0.278 [0.246, 0.310] | 19.9 % | 36.9 % | 14.1 % |
| **tiered `K`** | **0.760** [0.677, 0.843] | **50.1 %** | **41.0 %** | **17.0 %** |

**Caveat that matters (see Limitations §5.2):** flipping a statement of a tiered `K` produces a
`K'` that is *no longer a tiering*, so R2–R4 become free to fire — this arm tests *"is knowledge
elicited as a tiering more robust to misstatement?"* (answer here: **no, mildly worse**), **not**
Bang & Didelez's lemma, whose proper test is a perturbation of the *tier assignment* that keeps
`K'` tiered. That experiment is the #1 follow-up and it is not run here. **The R1-closed
sublattice finiteness result therefore remains unsupported by this pilot — neither confirmed nor
refuted.**

---

## 4. Verdict

### **SUPPORTED** — the falsifier did not fire.

The falsifier required **both** conditions. Neither holds:

| falsifier condition | required | **observed** |
|---|---|---|
| ignorance interval **unbounded at ρ = 1** | yes | **NO** — finite; median width **0**; excludes zero for **94.3 %** of SCMs |
| **near-universal** `O*` disruption | yes | **NO** — **36.9 %** disrupted, **14.1 %** harmful |

**A breakdown radius for background knowledge exists, is > 0 for ~90 % of problems, and is
informative enough to report.** The paper is alive and the pilot has also settled its §3: it must
be a **breakdown-radius / combinatorial** theorem (horn a), because the bias, when it lands, is
large and heavy-tailed (median 0.59 |τ|, max 13 |τ|) — a continuity bound is not available.

**What the pilot hands the paper, beyond a green light:**
1. **NUMBER 1 (36.9 % / 14.1 %) as the explanation of b-LOAD's own unexplained robustness** —
   closing arXiv:2607.04447's `[TO DO]` with a mechanism (`O*` is local; most noise is not).
2. **The utility ⊥ safety stratification** (0/75 vs 504) — a Theorem-1 candidate.
3. **The four-way outcome decomposition** (33 % caught free / 42 % inert / 15 % loud / 9 % silent) —
   a reviewer-legible frame, and the honest reason the answer is not "one wrong edge ⇒ unbounded bias".
4. **A negative on the tiered lever**, which redirects §3 rather than flattering it.

---

## 5. Limitations — stated in full

**5.1 The DGP is self-authored — partially mitigated, not cured.** *The campaign has already been
burned by a self-authored DGP passing as a general finding.* Random ER DAGs are **not** the space
of real causal problems, and every headline fraction is a number **about this graph ensemble**.
The §6b sweep honours the Herman et al. (arXiv:2503.17037) dense-graph caveat and confirms the
mechanism's *predicted* direction (more protective as `p` grows; less as density grows), which is
real evidence against artefact — but it **does not** cover: scale-free / hub graphs (where a
high-degree `X` would put many edges *inside* `pa(cn)` — plausibly the worst case and untested),
degree > 6, `p` > 10, non-Gaussian noise, or any real dataset. `NUM1 invalid` is still rising at
the edge of the grid. **Claim licensed: ER, `p ≤ 10`, degree ≤ 6. Nothing wider.**

**5.2 The tiered arm does not test Bang & Didelez.** Flipping statements breaks the tiering, so
the R1-only property is destroyed by the perturbation itself. The tiered results answer a
different (still useful) question. **The R1-closed finiteness result — the paper's stated
mathematical core — is untested.** Highest-priority follow-up: perturb the *tier assignment*
(move one node between tiers), keeping `K'` tiered, and re-measure the cascade.

**5.3 The perturbation model is "the expert reverses ρ of their claims."** It does not cover:
the expert asserting an edge that does not exist in the skeleton; *omitting* knowledge; knowledge
about non-adjacent pairs; or errors correlated across statements (a domain expert with one wrong
mental model gets *many* correlated statements wrong — arguably the realistic case, and it is
exactly the case where a ρ-ball with independent flips is optimistic).

**5.4 `|K| ≤ 4`.** The ball at ρ = 3 is then nearly the whole space of orientations of `K`, so
ρ = 3 is close to "the expert knows nothing." Larger `|K|` with the same ρ would give a
*relatively* smaller perturbation and, plausibly, a more favourable picture. **The ρ = 3 column is
the least trustworthy in this report.**

**5.5 `changed_still_valid = 0 / 754` is an empirical zero, not a theorem.** The observation that
`O*` *never* changed while remaining valid is striking and is a natural lemma candidate — but it
was observed, not proved, and only on this ensemble. Do not print it as a theorem without a proof
or a counterexample search.

**5.6 The `100 %` expulsion figure is entailed** by starting from a true `K`. It is reported for
completeness and must **not** be presented as a finding. The informative gap is 67 % (uncaught),
not 100 %.

**5.7 Population, not finite-sample.** The linear arm has no estimation error *by design* — this
buys decisiveness (no noise can mask the effect of a perturbation) and costs realism. The
ignorance interval reported is the **pure knowledge-induced** spread; a real analyst's interval is
this convolved with sampling error, and the two must be combined before any coverage claim. **The
conformal / finite-sample half of the paper — Zé's half — is not touched by this pilot.**

**5.8 Single X, single Y, no hidden confounding, causal sufficiency, faithfulness assumed.**
Sufficiency is load-bearing (it is also what makes the idea immune to D'Amour 1902.10286 —
there is no latent U to bite). Relaxing it (ADMGs) changes the machinery and is where
Taeb/Guo/Henckel 2511.10625 already lives.

**5.9 The nonlinear arm is secondary and finite-sample** (§6). Its estimator error (0.165 |τ|) is
*larger than the ρ=1 signal it measures* (0.116 |τ|), so it can support the qualitative
replication and nothing more. Its tiered cell (n = 16) is uninterpretable.

**5.10 The race is not addressed by this pilot.** Henckel is one author on all three nearest
competitors. Nothing here checks whether **Asiaee arXiv:2603.02204**'s Thm 1 transfers by
relabelling — the ranked idea flags that as a **free, do-it-first** gate. It remains open.

---

## 6. Nonlinear arm — the picture is not a linear-algebra artefact

60 nonlinear iSCMs (mechanisms `Σᵢ wᵢⱼ gᵢⱼ(xᵢ) + εⱼ`, `g ∈ {id, tanh, sin, x²−1, softplus}`,
internally standardized). `τ = E[Y|do(X=+1)] − E[Y|do(X=−1)]` by 100 000-sample interventional
Monte Carlo on the frozen SCM; estimates by gradient-boosted g-formula on `n = 4 000`.
**This arm is finite-sample and noisy — read it qualitatively only.**

| | generic `K` (59 SCMs) | tiered `K` (16 SCMs) |
|---|---|---|
| **estimator-noise baseline** (ρ=0, *correct* `O*`): `\|est−τ\|/\|τ\|` | **0.165** [0.118, 0.217] | 0.190 [0.077, 0.328] |
| **NUMBER 1** — consistent flips changing `O*` | **31.0 %** | 45.2 % |
| — of which making `O*` invalid | **10.6 %** | 19.0 % |
| `O*` changed but still valid | **0.0 %** | 0.0 % |
| ρ=1 width/\|τ\| median / mean | **0.000** / 0.116 [0.046, 0.207] | 0.000 / 0.219 |
| ρ=1 interval excludes 0 | **98.3 %** [94.9, 100] | 93.8 % |

**Both headline numbers replicate** (31.0 % vs 36.9 % disrupted; 10.6 % vs 14.1 % harmful; median
width 0; sign informative). The `changed-but-valid = 0` observation replicates too.

**One honest wrinkle that matters:** the ρ=1 knowledge-induced spread (**0.116 |τ|**) is *smaller
than the estimator's own error* (**0.165 |τ|**). At `n = 4 000` with a nonparametric estimator, a
single orientation error moves the answer **less than the regression does**. That is good news for
the thesis and bad news for this arm's resolution: it cannot separate small width differences from
noise, and it should not be used for any quantitative claim. `n_scm = 16` for tiered is too small
to interpret at all.

## 6b. Density & size sweep — Limitation 5.1, honoured

Herman et al. (**arXiv:2503.17037**) warn that dense graphs break this kind of result. **Tested
directly**: 150 linear iSCMs per cell, `p ∈ {6, 8, 10}` × expected degree ∈ {1.5, 2.5, 4.0, 6.0}
(generic `K`; degree 6 at `p=6` omitted as near-complete). Figure `fig3_density_sweep.pdf`.

| p | deg | \|U\| | Meek-inconsistent | **NUM1 changed** | **NUM1 invalid** | width med | width mean | **excludes 0** |
|---|---|---|---|---|---|---|---|---|
| 6 | 1.5 | 3.7 | 0.331 | 0.425 | 0.115 | **0.000** | 0.267 | 0.959 |
| 6 | 2.5 | 4.3 | 0.345 | 0.410 | 0.208 | **0.000** | 0.389 | 0.959 |
| 6 | 4.0 | 4.9 | 0.322 | 0.466 | 0.256 | **0.000** | 0.755 | 0.882 |
| 8 | 1.5 | 3.9 | 0.293 | 0.300 | 0.057 | **0.000** | 0.249 | 0.980 |
| 8 | 2.5 | 4.1 | 0.338 | 0.321 | 0.138 | **0.000** | 0.287 | 0.932 |
| 8 | 4.0 | 4.8 | 0.322 | 0.350 | 0.219 | **0.000** | 0.840 | 0.911 |
| 8 | 6.0 | 5.8 | 0.277 | 0.368 | 0.224 | **0.000** | 1.416 | 0.915 |
| 10 | 1.5 | 4.0 | 0.244 | 0.242 | 0.077 | **0.000** | 0.177 | 0.973 |
| 10 | 2.5 | 4.2 | 0.277 | 0.237 | 0.111 | **0.000** | 0.524 | 0.938 |
| 10 | 4.0 | 4.3 | 0.317 | 0.270 | 0.196 | **0.000** | 0.821 | 0.933 |
| 10 | 6.0 | 4.7 | 0.333 | 0.337 | **0.259** | **0.000** | 1.380 | **0.896** |

Three things, and the third is the one to watch:

1. **The `O*`-locality mechanism is confirmed as the explanation.** At fixed degree, growing `p`
   makes the knowledge *more* robust (deg 1.5: 0.425 → 0.300 → 0.242 as p = 6 → 8 → 10; deg 2.5:
   0.410 → 0.321 → 0.237). This was **predicted in advance** from the mechanism — more nodes at
   fixed degree ⇒ more edges far from `(X,Y)` ⇒ fewer perturbations reach `pa(cn) \ forb`. The
   prediction held in 2/2 series. This is the strongest evidence in the pilot that the effect is
   *the mechanism* and not an artefact of `p ≤ 8`.
2. **Density erodes it gracefully, and never breaks it.** Silent-bias rate rises monotonically with
   degree (0.057 → 0.259) and mean width reaches 1.4 |τ| at degree 6 — but **the median width is
   exactly 0.000 in all 11 cells**, and the interval still excludes zero for **89.6–98.0 %** of
   SCMs in every cell. **The falsifier does not fire in any cell**, including the densest.
3. **⚠️ The extrapolation is not safe.** `NUM1 invalid` is still *rising* at degree 6 and shows no
   plateau. Nothing here says it stays below the falisifier at degree 10, or on scale-free /
   hub-structured graphs, which were not tested. **The claim licensed by this sweep is "holds for
   `p ≤ 10`, ER, degree ≤ 6", not "holds generally."**

---

## 7. Reproducibility

```bash
# on betelgeuse (16 CPU cores, no GPU needed, no SLURM)
ssh costaj@100.110.205.73
cd ~/latent-causal/rho-breakdown-knowledge/code
PY=~/miniforge3/envs/tfmgpu/bin/python

$PY test_machinery.py                                    # ~3 min, MUST print ALL TESTS PASSED
$PY run_linear.py 600 ../results/linear_raw.json         # ~40 s on 16 cores
$PY analyse.py ../results/linear_raw.json ../results     # ~10 s
# NOTE: pin threads or sklearn's internal pool fights the process pool (load avg hit 161/16)
OMP_NUM_THREADS=1 $PY run_nonlinear.py 60 ../results/nonlinear_raw.json   # ~8 min on 16 cores
$PY analyse_nonlinear.py ../results/nonlinear_raw.json ../results
OMP_NUM_THREADS=1 $PY sweep_density.py 150 ../results/density_sweep.json  # ~6 min, 11 cells
$PY analyse_sweep.py ../results/density_sweep.json
$PY figures.py ../results                                # -> fig1/fig2 .pdf + .png
$PY figure3.py ../results/summary_sweep.json             # -> fig3
```

**Seeds.** Sweep: base `3000000 + cell_index*100000 + k`. Linear: master `20260717`; job seeds `0…1799` (`p`, `deg` drawn from it); each SCM's
seed is stored in `linear_raw.json`. Diagnostics use seed + 1e7. Nonlinear: master `20260718`, job
seeds `500000…`. Bootstrap CIs: 2 000 resamples, seed 7 (linear) / 11 (nonlinear); all intervals
are percentile bootstrap over **SCMs**.

**Files.**
- `code/graphs.py` — DAG/CPDAG/Meek R1–R4/MPDAG/consistency/brute-force MEC
- `code/adjust.py` — cn / forb / possDe / amenability / `O*` / generalized adjustment criterion / population OLS
- `code/scm.py` — iSCM (linear + nonlinear) + naive control + var-/R²-sortability
- `code/test_machinery.py` — **the validation suite (T0–T5)**
- `code/run_linear.py`, `code/run_nonlinear.py`, `code/analyse.py`, `code/analyse_nonlinear.py`, `code/figures.py`
- `results/linear_raw.json` (every member of every ball), `results/summary.json`, `results/plotdata.json`
- `code/sweep_density.py`, `code/analyse_sweep.py`, `code/figure3.py` — the §6b robustness sweep
- `results/fig1_breakdown_radius.pdf` — **the headline figure**
- `results/fig2_diagnostics_and_bias.pdf` — sortability + the bias-is-large result
- `results/fig3_density_sweep.pdf` — the density/size robustness grid
- `results/summary_nonlinear.json`, `results/summary_sweep.json`, `results/density_sweep.json`

**Not reused:** Chris's `/Users/josecosta/LOAD_in_MPDAG` (causal-learn + rpy2 dependency; the Meek
and `O*` machinery here is independent and brute-force-validated, which is a feature — it means
b-LOAD's numbers and these were not produced by the same possibly-wrong closure).
