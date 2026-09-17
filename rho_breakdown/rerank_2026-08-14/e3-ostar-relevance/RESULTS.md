# E3 — O\*-relevance vs edit distance: RESULTS

**Run 2026-08-14 on Zé's Mac, base env `/Users/josecosta/miniforge3/bin/python3` (3.12.7, numpy
2.3.3, scipy 1.16.3, networkx 3.5). Total compute ≈ 30 s. No GPU, no betelgeuse, no new simulation.**

Pre-registration: [`PREREGISTRATION.md`](PREREGISTRATION.md), written and saved **before** any number
was computed.

```
cd /Users/josecosta/mugango/output/2026-08-14_latent-causal-iclr2027-rerank/e3-ostar-relevance
python3 e3_features.py     # regenerate (D,C), emit 8,085 feature rows      (0.7 s)
python3 e3_analyse.py      # pre-registered head-to-head, cluster bootstrap (22 s)
python3 e3_decisive.py     # matched-granularity comparators, counterexample (4 s)
```

---

## 0. VERDICT IN ONE PARAGRAPH

**The pre-registered rule as literally written returned SUPPORTED. I am overriding it to
PARTIAL / ABSORBED, because my own honest-risk comparator was mis-operationalised.** The O\*-relevance
metric beats edit distance enormously and unambiguously (**AUC 0.842 vs 0.637, ΔAUC +0.205
[+0.187, +0.224]**) — so "edit distance is the wrong yardstick" is *measured, not asserted*. But
against a **matched-granularity** locality feature that uses only hop-distance to the query and no
O\* machinery at all, the gain collapses to **+0.029 [+0.015, +0.042]** — below the 0.05 bar I set in
advance for "this is a section". And a post-hoc decomposition shows the specifically-O\* part of the
metric (`pa(cn) \ forb` membership) is **anti-predictive, AUC 0.409**. The predictive content is
**causal-path locality**, not optimal-adjustment-set membership. **Contribution: one sentence and
one figure in #1's §4, not a section, and the campaign's #4 abstract is over-claimed.**

---

## 1. Integrity gates — all passed

`D` and `C` are not stored in `linear_raw.json`; they are regenerated from the stored
`(seed, p, deg)` triple. Every gate in the pre-registration passed on **600/600 SCMs**:

| gate | check | result |
|---|---|---|
| G1 | regenerated `n_undirected` == stored | 600/600 |
| G2 | regenerated `n_edges` == stored | 600/600 |
| G3 | regenerated `cpdag_amenable` == stored | 600/600 |
| G4 | recomputed `O0` from stored `K` == stored `O0` | 759/759 arms |
| G5 | flip-set enumeration order reproduces every member's `rho` | 8,085/8,085 |

`gate_failures = {}`. **759 arms, 8,085 perturbation records, 754 silent (9.33%)** — reproducing the
ground-truth recomputation of 754/8085 = 9.3% exactly.

**C13 capability check (the vault rule: a null must show its instrument could produce a non-null).**
The oracle feature `ostar_changed` scores **AUC = 1.0000** on the primary population through the
identical code path, so the AUC machinery has full range on this outcome. ⚠️ Note this oracle is a
*tautology*, and that is itself a finding worth stating in the paper: on the undetectable population,
`silent ⟺ ostar_changed` **exactly** — because `changed_still_valid = 0/754`. There is no such thing
as an O\* that moved and stayed valid in this ensemble.

---

## 2. The metric, written down (the campaign never did)

`O* = pa(cn(X,Y,G0), G0) \ forb(X,Y,G0)`. Define the **O\*-defining region**
`R = {X,Y} ∪ cn ∪ pa(cn) ∪ forb` — every node that appears anywhere in that expression. For a
statement `k = (u→v) ∈ K`:

```
r(k)  = 1 if {u,v} ∩ R ≠ ∅ else 0
r2(k) = 2 if {u,v} ∩ (cn ∪ {X,Y}) ≠ ∅      # causal-path core: can destroy AMENABILITY
        1 elif {u,v} ∩ (pa ∪ forb) ≠ ∅     # adjustment machinery: can change the ESTIMATE
        0 otherwise                         # appears in no term of the formula
d_ostar(F) = Σ_{k∈F} r(k)     d_ostar_graded(F) = Σ_{k∈F} r2(k)     d_edit(F) = |F| = ρ
```

No coefficient is fitted anywhere. The two tiers are the two syntactic roles a node can play in
`pa(cn) \ forb`, ordered a priori by "destroys the estimand" > "changes the estimate" > "absent".
Everything is computable from `(C, K, X, Y)` — **no true DAG, no perturbation enumeration.**

**Marginal distribution (generic arm, 2,058 statements):** `r=1` for **1,682 (81.7%)**; `r2=2` for
1,398; `r2=1` for 284; `r2=0` for 376. Statements touching X or Y: 1,264. **Crucially, of the 794
statements NOT incident to X or Y, 418 (52.6%) are still O\*-relevant** — so the metric is genuinely
not identical to query-incidence, and `d_ostar ≢ d_edit` (no degeneracy).

---

## 3. HEAD-TO-HEAD — primary population

**P3 "undetectable" = Meek-consistent ∧ amenable, generic arm: n = 2,052 members, 569 silent, 561
SCMs.** This is the only population where the metric has a decision to inform: the analyst ran the
pipeline, got a number, and nothing looked wrong. Event = `ostar_valid == False` (silent bias).
All CIs are **cluster bootstraps resampling the 561 SCMs** (1,000 reps; 2,000 for §5).

| feature | what it knows | AUC | 95% CI |
|---|---|---|---|
| `ORACLE ostar_changed` | *post*-perturbation (no analyst has it) | **1.0000** | — |
| **`d_ostar_graded`** | **O\* region, graded — PRE-REGISTERED PRIMARY** | **0.8417** | [0.8252, 0.8591] |
| `n_xy_inc` | count of wrong statements touching X or Y | 0.8131 | [0.7930, 0.8338] |
| `d_dist_graded` | hop-distance, graded 2/1/0 — matched granularity | 0.8089 | [0.7881, 0.8296] |
| `d_ostar` | O\* region, binary relevance count | 0.7650 | [0.7429, 0.7856] |
| `d_gdist` | −(min hop-distance over F) | 0.7639 | [0.7438, 0.7841] |
| `d_xy_inc` | binary: any wrong statement touches X or Y | 0.7638 | [0.7437, 0.7838] |
| `n_dist_le1` | count of statements within 1 hop | 0.7159 | [0.6920, 0.7393] |
| **`d_edit`** | **ρ = number of wrong orientations — BASELINE** | **0.6367** | [0.6150, 0.6576] |
| `d_ostar_any` | binary: any relevant statement wrong | 0.6200 | [0.6049, 0.6344] |

### Pre-registered decision, condition by condition

| # | condition | measured | verdict |
|---|---|---|---|
| 1 | `A_ost ≥ 0.70` | 0.8417 | ✅ PASS |
| 2 | `A_ost − A_edit ≥ 0.10`, CI>0 | **+0.2051 [+0.1869, +0.2239]** | ✅ PASS |
| 3 | `A_ost − A_xy_inc ≥ 0.05`, CI>0 | **+0.0780 [+0.0599, +0.0942]** | ✅ PASS (see §4) |

**Literal pre-registered verdict: SUPPORTED.**

### Continuous outcome — Spearman ρ with |bias|/|τ| (same population, n = 2,052)

`d_ostar_graded` **+0.4759** [+0.425, +0.520] · `d_dist_graded` +0.4426 · `n_xy_inc` +0.4222 ·
`d_ostar` +0.4137 · `d_gdist` +0.3941 · `d_xy_inc` +0.3755 · **`d_edit` +0.2035** [+0.159, +0.245].
Same ordering, same story: a large win over edit distance, a small one over locality.

---

## 4. 🔴 THE HONEST RISK FIRED — my condition (3) comparator was handicapped

The pre-registration named condition (3) *"the honest-risk comparator … if `d_ostar` only re-encodes
'is the statement incident to X or Y', the contribution is one sentence, not a section."* That was
the right question **operationalised with the wrong instrument**: `d_xy_inc` is **binary** while
`d_ostar_graded` ranges 0–6. Any graded score beats a binary one on AUC whenever signal exists,
because ties inside the two bins score 0.5. The +0.078 is partly a granularity artefact, not
information.

Rebuilt with **matched-granularity** comparators carrying *only* hop-distance — no `cn`, no `pa`, no
`forb`, no Meek closure beyond what produced `G0`:

| paired ΔAUC (cluster bootstrap) | point | 95% CI | vs the 0.05 bar |
|---|---|---|---|
| `d_ostar_graded` − `d_edit` | **+0.2051** | [+0.1869, +0.2239] | far above |
| `d_ostar_graded` − `d_xy_inc` (pre-registered, **binary**) | +0.0780 | [+0.0599, +0.0942] | above |
| `d_ostar_graded` − `n_xy_inc` (**count**, matched) | **+0.0287** | [+0.0154, +0.0416] | ❌ **below** |
| `d_ostar_graded` − `d_dist_graded` (**graded 2/1/0**, matched) | **+0.0329** | [+0.0259, +0.0399] | ❌ **below** |
| `d_ostar` − `n_xy_inc` | **−0.0480** | [−0.0737, −0.0259] | ❌ **worse than trivial** |

**Under the correctly-operationalised version of my own stated intent, condition (3) FAILS.** The
gain is real (CIs exclude 0) but is roughly **0.03 AUC**, not the ≥0.05 I fixed in advance as the
threshold for "the estimand induces a metric the graph does not". The last row is the sharpest: the
*binary* relevance metric `d_ostar` is **significantly worse** than simply counting how many wrong
statements touch the query.

**Nested logistic** (`silent ~ d_dist_graded + d_edit` vs `+ d_ostar_graded`): LR = 296.1, df 1,
β_ostar = **+2.539**. Direction and magnitude are real. ⚠️ The nominal p (2.4e-66) is
anti-conservative — members are clustered in SCMs and this test ignores that. Read the sign, not the p.

---

## 5. POST-HOC DECOMPOSITION — what actually carries the signal

**Labelled post-hoc: these features were built after seeing §3–4 and are diagnostic, not confirmatory.**
Splitting `d_ostar_graded = 2·n_core + 1·n_mach`:

| component | meaning | AUC |
|---|---|---|
| `n_core` | wrong statements touching `cn ∪ {X,Y}` (a causal-path node) | **0.8502** |
| `n_xy_inc` | wrong statements touching `{X,Y}` | 0.8131 |
| `n_core_notxy` | touching the `cn` *interior* only (mean 0.081/member) | 0.6063 |
| **`n_mach`** | **touching `pa ∪ forb` only — the specifically-O\* tier** | **0.4092** |
| `2·n_core` (drop the machinery term entirely) | | **0.8502** |
| `2·n_core + n_mach` = `d_ostar_graded` (as pre-registered) | | 0.8417 |

Paired cluster bootstrap, 2,000 reps: `n_core − n_xy_inc` = **+0.0371 [+0.0283, +0.0467]**;
`n_core − d_ostar_graded` = +0.0084 [−0.0015, +0.0183]. `corr(n_core, n_xy_inc) = 0.936`.

**Three consequences, and they are the scientific content of this experiment.**

1. **The `pa(cn) \ forb` tier — the part that makes the metric about the *optimal adjustment set*
   rather than about causal paths — is anti-predictive (AUC 0.409, below chance).** Including it
   *lowers* the metric's AUC from 0.850 to 0.842. My own graded weighting is suboptimal, and the
   reason is mechanical: `forb` is large (in the seed-508 case it is all 7 nodes), so tier-1
   membership is close to "not tier 2", which is negatively associated with damage.
2. **The whole useful content is "does the wrong statement touch a node on a possibly-causal path
   from X to Y".** That is 93.6% correlated with "does it touch X or Y" and buys +0.037 over it.
3. **Therefore the honest name for the mechanism is CAUSAL-PATH LOCALITY, not O\*-relevance.** The
   campaign's #4 title — *"O\*-Relevance as the Right Metric on Background Knowledge"* — attributes
   the effect to the wrong object.

---

## 6. 🔴 A COUNTEREXAMPLE TO THE CAMPAIGN'S HEADLINE, VERIFIED FROM SCRATCH

The campaign's #4 abstract states: *"an error outside O\*'s neighbourhood costs exactly zero"*, and
the vault repeats *"silent bias exactly 0 at graph-distance ≥1 (0/791 statements)"*.

**I reproduce that at ρ = 1 and I break it at ρ ≥ 2.**

At **ρ = 1**, generic arm, every statement, by hop-distance to `{X,Y}`:

| dist | 0 | 1 | 2 | 3 | ≥4 / unreachable |
|---|---|---|---|---|---|
| n | 1,264 | 540 | 93 | 19 | 2 |
| silent | **194** | **0** | **0** | **0** | **0** |
| rate | 15.35% | 0 | 0 | 0 | 0 |

**0/654 at distance ≥1 — the locality result holds exactly for single errors.**

At **ρ ≥ 2 it does not.** Pooling all ρ on P3, by *minimum* distance over the flip-set:

| min dist | 0 | 1 | 2 | 3 | ≥4 |
|---|---|---|---|---|---|
| n | 1,266 | 529 | 79 | 17 | 2 |
| silent | 568 | **1** | 0 | 0 | 0 |

**The one case, re-derived from scratch by regenerating the SCM independently of the stored JSON:**

```
seed 508, p=7, deg=2.0, x=4, y=0, tau=-0.138462
K = [(6,5), (1,3), (4,6), (5,2)]          O0 = {}  (matches stored)
node hop-distances to {x,y}: {0:0, 4:0, 1:1, 5:1, 6:1, 2:2, 3:2}
flip {0,1}: reverse (6->5) and (1->3)  — NEITHER is incident to X or Y, both at distance 1
perturbed K = [(5,6), (3,1), (4,6), (5,2)]
perturbed O* = {1,5}     valid in the true DAG? FALSE     |bias|/|tau| = 0.318
```

Two statements, each individually harmless and each at distance ≥1, **jointly** move O\* from `{}` to
an invalid `{1,5}` and produce 32% silent relative bias. This is a **Meek-cascade interaction**: R1
chains along chordless undirected paths of arbitrary length, so the influence cone is not hop-bounded.

**Consequences for the paper, both binding.**
- The sentence *"an error outside O\*'s neighbourhood costs exactly zero"* **may not be written.** It
  is true for ρ = 1 and false for ρ = 2 on our own data. It must become *"single orientation errors at
  graph-distance ≥1 produced zero silent bias (0/654); joint errors did not — 1/529 at min-distance 1"*.
  This is the same discipline the vault already imposes on `changed_still_valid = 0` (empirical zero,
  not theorem) and on the no-continuity-bound rule.
- **This single case is the strongest argument FOR the O\*-relevance metric**, and it should be the
  figure: `d_ostar_graded = 4` flags it (both statements touch `cn = {0,1,5,6}`), while every
  distance-based feature misses it. The metric buys exactly the cases where the cascade reaches in.
  It buys them at a rate of ~1 in 529, which is why it is worth a sentence and not a section.

### ⚠️ The campaign's "0/791" does not reproduce as a denominator

Under my definition (min over endpoints of skeleton hop-distance to `{X,Y}`) the ρ=1 distance-≥1
denominator is **654** (generic arm) or **1,002** (both arms). Neither is 791. The *qualitative*
result — exactly zero — reproduces in both. The number 791 came from the verification agent's own
re-run and is **not recoverable from the on-disk file**; do not print it. **165 of 2,673 statements
(6.2%) have an endpoint unreachable from `{X,Y}` in the skeleton** and land in the "≥4" bucket — a
labelling detail the campaign never stated.

---

## 7. Where the discrimination actually lives

Stratifying P3 (generic) on whether any wrong statement touches X or Y:

| stratum | n | silent | rate | AUC `d_ostar_graded` | AUC `d_edit` |
|---|---|---|---|---|---|
| `xy_inc = 0` | 786 | **1** | **0.13%** | 1.000 (1 positive — meaningless) | 0.796 |
| `xy_inc = 1` | 1,266 | 568 | **44.87%** | 0.709 | 0.596 |

**A 350× rate difference.** Query-incidence alone is a near-perfect screen; inside the dangerous
stratum the O\* grading buys 0.709 vs 0.596 for edit distance. The `xy_inc = 0` AUC of 1.000 rests on
a **single** positive event and must not be quoted.

---

## 8. Limitations

1. **Scope is the pilot's scope, unwidened:** ER graphs, p ∈ {5..8}, expected degree ∈ {1.5, 2.0, 2.5},
   |K| ∈ {3,4}, ρ ≤ 3, linear-Gaussian iSCM, single X, single Y, causal sufficiency, faithfulness,
   population (no sampling error). Nothing here is licensed on scale-free graphs, p > 10, degree > 6,
   or finite samples.
2. **|K| ≤ 4 caps ρ ≤ 3**, so "joint errors break locality" is measured on flip-sets of size ≤ 3 only.
   The single counterexample is at ρ = 2; the rate at larger |K| is unmeasured and could be much higher.
   This is the same instrument degeneracy that E1 exists to fix, and it bounds this result too.
3. **The metric is one of many possible operationalisations.** `r2`'s tier weights were fixed a priori
   and turned out suboptimal (§5). A fitted weighting would score higher and would no longer be
   assumption-free — and would need a train/test split this analysis does not have.
4. **§5 is post-hoc.** `n_core` was constructed after seeing §3–4. It is a diagnostic decomposition,
   not a confirmatory test, and must be re-derived on fresh data before being printed as a result.
5. **The tiered arm is reported but not used for the verdict** (`P3_undetectable_tiered` is in
   `results/e3_analysis.json`); the pre-registration named `generic` as primary.
6. **AUC treats all silent-bias events alike**, but the pilot's own §3 shows damage is heavy-tailed
   (median 0.591|τ|, max 13.137|τ|). The Spearman column partly covers this; a cost-weighted metric
   does not exist here.

---

## 9. What this changes for the ICLR paper

- **Novelty leg (ii) — "O\*-relevance reweighting" — cannot be a headline or a section.** It is
  measured at +0.03 AUC over a feature any reader can compute from the query alone. The refinement
  phase had already demoted it from headline to positioning on prior-art grounds (gadjid 2402.08616,
  Fang 2207.05067); **this experiment independently demotes it on evidence.** Two independent
  demotions of the same leg.
- **What survives, and it is worth writing:** *edit distance is measurably the wrong yardstick*
  (ΔAUC +0.205 over ρ, a large and tight effect), and *the right yardstick is causal-path locality*.
  That is one sentence, one stacked-bar figure (the four-way decomposition by distance), and the
  seed-508 counterexample as the honest caveat.
- **It strengthens the "utility ⊥ safety" abstract sentence** the ICLR-anatomy constraints already
  mandate: damage is concentrated 350× in the stratum where knowledge touches the query.
- **It forces a phrasing correction** in any b-LOAD camera-ready text and in #1: the zero-at-distance
  result is a **ρ = 1** statement with a **known counterexample at ρ = 2**.
