# X2 RESULTS — locality: can a statement far from the query move `O*`?

**Run 2026-08-19 on betelgeuse.** Pre-registration: `PREREG.md` (frozen before implementation).
Binding design audit: `AUDIT.md`. Every mitigation marked *fatal* or *major* was honoured; the
three places where I disagree with the auditor are recorded in §11 and were implemented anyway.

---

## 0. Headline

> **The falsifier fired. Two of the three sentences are refuted outright; the third survives only
> with a `ρ = 1` qualifier it does not currently carry.**
>
> * **S5 as written — REFUTED.** 70 knowledge sets in E1′'s own archive change `O*` while **every**
>   misstated statement sits at finite skeleton distance ≥ 1 from the query. All 70 are silent
>   **bias** (the new `O*` is invalid in the true DAG). 61 at ρ = 2, 9 at ρ = 3. Smallest instance
>   **p = 5**, exhibited in full with its Meek rule chain (§4.4). The number "0 of 791" is not
>   repaired by this run; it is replaced.
> * **The radius is sharp, and that is the useful result.** Every one of the 70, and every one of
>   the 980 identification losses, involves a misstatement at distance **exactly 1**. Over the whole
>   ball (ρ ≤ 4), knowledge sets whose *closest* misstatement is at distance ≥ 2 produced
>   **0 silent bias and 0 identification loss in 6,726 amenable perturbations**. The elicitation
>   radius for the reversal operator is **1** — not 0, not unbounded.
> * **S5 restricted to ρ = 1 — NOT refuted, and now backed by something better than a sample.**
>   Exhaustive single-statement sweeps found **0** `O*`-changes at `dmin ≥ 1` over **54,079**
>   misstatement trials across five ensembles (`SUPPORTED` / `SUPPORTED-WEAK` / `INCONCLUSIVE` per
>   ensemble by the M7 thresholds — see §4.2 for which is which, no ensemble was enlarged after the
>   fact), **and** over the *complete* enumeration of every CPDAG, every Meek-consistent knowledge
>   state, every query and every statement at **p ≤ 5** — 0 of 136,200 far trials against
>   184,332 of 307,584 near ones, with no sampling error in that zero.
> * **S8 — NOT SUPPORTED AT ANY RADIUS.** Refuted at `r = 0`: `Q₀` = 0.0143–0.5231, with the
>   Clopper–Pearson lower bound strictly above the 0.01 REFUTED bar in **4 of 5** usable cells
>   (`large` misses it by 10⁻⁵: CI [0.009999, 0.019869]). At `r = 1` the pruner is safe but prunes
>   nothing (mean `Cost₁` = 0.28–0.97; `Q₁` CI spans the 0.001 SUPPORTED bar on
>   `original`/`licensed`). No radius satisfies both halves of the registered bar.
> * **S9 — REFUTED.** S9 asserts `r*_elicit = 0`. Measured: **1** under the reversal operator and
>   **≥ 4** under the addition operator. The reason is not that far knowledge is dangerous but that
>   it is **valuable** — a single statement 1 to ≥ 4 hops out turned a *non-identified* query into
>   an identified one **29,364** times (`E_gain`, never measured by any prior experiment).
>   The paper tells the practitioner to skip exactly that knowledge.
> * The apparent E1′-vs-S5 contradiction (0.935 censoring at distance 1) is **not** a refutation.
>   It is 37 of 38 accounted for by `ρ*_any` folding in **loud** identification failure, and 1 of 38
>   by a genuine ρ ≥ 2 `O*`-change. Resolved in code, §3.

Every fraction below carries `N`, `D_s5`, `D_am`, `D_con`, `D_all`. Every stratum carries its `n`.
`disconnected` is never merged into any "≥ k" bucket.

---

## 1. What was pre-registered, what the audit changed, and what actually ran

`AUDIT.md` returned **FLAWED** with six *fatal* and seven *major* findings. The design that ran is
the pre-registration **as amended by M1–M19**, not the pre-registration as written.

| # | mitigation | how it was honoured |
|---|---|---|
| M1 | per-ensemble G4 bands, G4 relabelled a reproduction check | gate **G4a**, §2 |
| M2 | ARM A demoted to mechanism/`E_gain`; out of the S5 denominator; exhaustive `p ∈ {4,5,6}` lemma | ARM A kept and reported; §7 lemma at p = 4, 5 exhaustive, p = 6 **sampled** (deviation, §10) |
| M3 | new primary S5 arm: base `K` + **exhaustive** statement sweep | **ARM E**, §4.2 |
| M4 | ARM D moved to disclosure; §4.2's prediction paragraph deleted | §5 states ARM D is recovered from archived data; no prediction is claimed |
| M5 | `Cost_r` aggregator = mean over SCMs + per-SCM fraction; Pareto frontier | §5 |
| M6 | S9 bars replaced by `r*_gain` / `r*_ident` / `r*_elicit` + rates with own denominators | §6 |
| M7 | S5 thresholds per ensemble from realised `D_s5` (SUPPORTED ≥ 20,000 · WEAK ≥ 3,000 · below = INCONCLUSIVE), budgeted before the run | §4.2; **the thresholds were never moved and `n_SCM` was never raised after seeing a result** — see §10.a |
| M8 | `D_s5 = D_am ∧ ¬true_in_D ∧ finite δ`, printed as a fifth column; `Ω` on finite `dmin` only | every S5 table |
| M9 | naive **and** cluster rule-of-three bounds | §4.2 |
| M10 | `dmax` promoted to a headline row | §4.3 |
| M11 | §1.2's superset/subset inversion corrected; skeleton `r*_elicit` = a **lower** bound | §6, §10 |
| M12 | S3/S4/S5 as one contingency table with the decomposition | §8 |
| M13 | within-SCM matched near/far contrast | §4.2 |
| M14 | composite `P(E_set ∪ E_ident ∪ E_gain \| consistent)` as co-primary | §6 |
| M15 | `E_ident`/`E_gain` split into `nopath`/`unamenable` at record time | §6 — and it is a **qualitative** split, not a cosmetic one |
| M16 | ARM C candidate set on `C`, τ floor fixed (1e-6 draw / 1e-3 sensitivity) | `x2lib.draw_scm(query_rule="argmax_delta_C")`; ensembles `deep`, `deep8` |
| M17 | G1/G2/G3/G5 relabelled unit tests; **G10** placebo added | §2 |
| M18 | ambiguities resolved; `P(cascade > 0)` registered | §4.2, §5 |
| M19 | Clopper–Pearson on `Q_r`; cells with `D_am < 300` excluded from verdicts | §5 |

**Machinery tests were written and run before any experiment** (`code/test_machinery.py`),
in the spirit of the pilot's own `test_machinery.py`. **17 of 17 pass, in 29.8 s.** The two that
carry the most weight:

* **T5** — the one-pass `O*` used everywhere in X2 was checked against **both**
  `adjust.optimal_adjustment_set` and an independently written brute-force implementation
  (BFS-over-path-prefixes, transitive-closure `poss_de`), **exhaustively over every CPDAG and every
  reachable single-statement MPDAG at p ≤ 5 and every ordered query: 0 mismatches over 756,740
  (graph, query) triples.**
* **T6** — `is_valid_adjustment_set`, on which the whole `E_bias` notion rests, was checked against
  the only independent ground truth available: a valid adjustment set must give population OLS
  `β = τ` exactly. **0 violations over 3,272 valid sets**; 0 of 4,286 invalid sets coincidentally
  equal.

Full list, with counts, in `results/machinery_tests.json`.

---

## 2. Gates — all of them, pass or fail, before any verdict

| gate | what it is | result |
|---|---|---|
| machinery tests | 17 brute-force checks (§1) | **17/17 PASS** |
| **G3** replication (M17: a unit test, cannot fail on the science) | my independent implementation re-derives every structural field of the E1′ archive from `(seed, p, deg)` alone — `n_edges`, `n_undirected`, `cpdag_amenable`, `stmt_hopdist`, `O0`, and per member `consistent`, `amenable`, `ostar_changed`, `ostar_valid`, `cascade`, `expels_true_dag` | **PASS — 0 mismatches over 176,592 archived members in 6 cells** |
| **G4a** positive control, reproduction (M1) | `L_set(0, ρ=1)` per statement, `D_am` denominator, against the bands the AUDIT fixed from archived data before implementation | **PASS 5/5** — original 0.3701 [0.32,0.42] n=2,721 · licensed 0.4660 [0.42,0.52] n=1,944 · large 0.2027 [0.15,0.25] n=1,026 · k8/K4 0.7582 [0.71,0.81] n=91 · k8/K8 0.5552 [0.51,0.61] n=1,077 |
| **G4b** ARM E dynamic range (floor 0.05, fixed before ARM E ran) | the new primary arm's detector must fire where events are known to be | **FAIL on `large`** (see below) |
| **G7** sentinel | no `1<<20` leaks into a "≥ k" bucket | **PASS** — and T2 proves it is decidable: a skeleton edge is *wholly* inside or *wholly* outside the query's component, so `disconnected` is a property of the statement, not of an endpoint |
| **G9** abort accounting | path-enumeration guard | **PASS — 0 aborts in every ensemble** (`PATH_CAP` 200,000, `SCM_CAP` 20 s never tripped) |
| **G10** placebo ball (M17, the one gate that could fail) | rerun the whole ARM D partition on E1′'s `members_frozen` ball; `Q_r` must be 0 for every `r` | **PASS — 0 of 39,400 (SCM, r) cells** |
| **G11** null-trial placebo (new) | in ARM E a swept *true* statement on an edge already in the base rebuilds the base exactly (T7c); it must fire nothing | **PASS — 0 events over 144,740 null trials** |
| G1, G2, G5 | Meek idempotence, skeleton invariance, dead-switch | relabelled **unit tests** per M17; run and passed (T4a, T4b, T9); they cannot fail |

🔴 **G4b failed on `large` and I am not smoothing it.** ARM E's `L_set(0, ρ=1)` is 0.1016
(licensed), **0.0363 (large — below the 0.05 floor)**, 0.0868 (k8), 0.0879 (deep), 0.0689 (deep8).

The cause is measured, not argued. ARM E's near stratum mixes two operators and one of them is
**structurally incapable** of moving `O*`. Splitting them (misstatements only, `dmin = 0`):

| ensemble | `rep` (reverse a statement in the base — the campaign's own operator) | `add` (assert a new required edge on an untouched edge) | ARM B's archived `L_set(0)` for comparison |
|---|---|---|---|
| `licensed` | **3,355 / 7,305 = 0.4593** | **0 / 2,858 = 0.0000** | 0.4660 |
| `large` | **391 / 1,922 = 0.2034** | **0 / 1,920 = 0.0000** | 0.2027 |
| `k8` | **218 / 297 = 0.7340** | **0 / 328 = 0.0000** | 0.7582 (`k8`/K4) |
| `deep` | **388 / 498 = 0.7791** | **0 / 101 = 0.0000** | — (ARM C ensemble) |
| `deep8` | **71 / 83 = 0.8554** | **0 / 36 = 0.0000** | — (ARM C ensemble) |

The `rep` sub-arm reproduces the archived ARM B rate to within 0.007, 0.001 and 0.024 on the three
boxes where the comparison exists — and it does so through **X2's own, independently written code
path**, on **different SCM draws** (ARM E uses seed 20260820/20260821, ARM B used 20260819), which
is a stronger statement about the detector than the archive-based G4a. The blended rate is low because half the sweep is the `add` operator, whose zero
§7 shows to be a **structural fact** (0 of 498,792 at p ≤ 5, exhaustively, at every distance
including 0) rather than a dead detector.

So G4b's floor was set on the wrong object — a blend of two operators with different physics — and
**the zero at `dmin ≥ 1` remains readable**: G4a passes 5/5 on the archive-backed operator that
contains every known event, the `rep` sub-arm reproduces it inside ARM E's own code path, the
within-SCM matched contrast (§4.2) reproduces the near/far gap on the *same* graphs, and the p ≤ 5
result is an enumeration rather than a sample. I record G4b as a **defect in the gate's
specification**, not as a licence to ignore a failed gate.

---

## 3. Step 0 — the E1′ / S5 contradiction, dissolved in code

E1′ reported `ρ*_any` censoring of 0.935 at "distance 1" (`licensed`, n = 229) and 0.977 (`large`,
n = 638), i.e. 6.5 % and 2.3 % **not** censored, against S5's "0 of 791". X2 reproduces E1′'s table
from raw and then decomposes it.

**Reproduction (n = 20,000, `TAU_FLOOR` = 1e-3, per-SCM `HOPMIN` strata — E1′'s own unit):**

| cell | n analysed | dist 0 | dist 1 | dist 2 | dist 3 | dist ≥4 | **disconnected** |
|---|---|---|---|---|---|---|---|
| `original`/K4 | 2,637 | 0.0374 (n=2,539) | 1.000 (n=93) | 1.000 (n=1) | — | — | 1.000 (n=4) |
| `licensed`/K4 | 2,087 | 0.0697 (n=1,852) | 0.9345 (n=229) | 1.000 (n=4) | — | — | 1.000 (n=2) |
| `large`/K4 | 2,443 | 0.3302 (n=1,472) | 0.9765 (n=638) | 1.000 (n=217) | 1.000 (n=39) | 1.000 (n=20) | 1.000 (n=57) |
| `k8`/K4 | 127 | 0.0300 (n=100) | 0.7778 (n=27) | — | — | — | — |
| `k8`/K6 | 329 | 0.0195 (n=307) | 0.9545 (n=22) | — | — | — | — |
| `k8`/K8 | 585 | 0.0160 (n=563) | 0.9545 (n=22) | — | — | — | — |

Every published E1′ cell reproduces exactly. Two cells E1′'s `RESULTS.md` omitted are printed here
(`k8`/K4 0.7778 at n=27 — the weakest in the corpus — and `k8`/K6), and the `large` "dist ≥ 2
(n=331)" row is broken into 217 + 39 + 20 **+ 57 disconnected**, which is the 17 % of that
denominator that is at *infinite* distance, not two hops.

**Decomposition of every non-censored `HOPMIN ≥ 1` SCM** (identical at n = 20,000 and n = ∞):

| cell | n at `HOPMIN ≥ 1` | not censored | **ident-only** | `ρ*_se = 1` | `ρ*_se ≥ 2` |
|---|---|---|---|---|---|
| `original` | 98 | 0 | 0 | 0 | 0 |
| `licensed` | 235 | 15 | **15** | 0 | 0 |
| `large` | 971 | 15 | **14** | 0 | **1** |
| `k8`/K4 | 27 | 6 | **6** | 0 | 0 |
| `k8`/K6 | 22 | 1 | **1** | 0 | 0 |
| `k8`/K8 | 22 | 1 | **1** | 0 | 0 |
| **total** | **1,375** | **38** | **37** | **0** | **1** |

Against the four hypotheses the caller named:

* **(c) `ρ*_any` folds in the identification criterion — 37 of 38.** `ρ*_any = min(ρ*_se, ρ*_ident)`
  (`analyse_arm1.py:140`). Those 37 SCMs lost *amenability*, which is **loud**: the analyst cannot
  report a number at all. They are not silent bias and they do not bear on S5, which is about
  `O*` moving. They bear on **S9**, and they refute it (§6).
* **(a) a genuine `O*`-change — 1 of 38, at ρ ≥ 2**, on `large`. This is not an artefact; it is the
  same phenomenon §4.4 exhibits 70 times.
* **(b) the min-vs-max endpoint artefact — not a bug, but the *whole* question.** `dmin(s) ≥ 1` is
  *identical* to "the statement is not incident to X or Y" (machinery test T3, 0 mismatches). Under
  the other reading (`dmax ≥ 1`) the claim flips from 0 % to 100 % true. §4.3.
* **(d) a bug — no.** Gate G3 re-derived every archived member independently: 0 mismatches over
  176,592.

**Nothing about `ρ*_se` fires at `HOPMIN ≥ 1` in any cell at ρ = 1.** The 0.935 and the "0 of 791"
were never in contradiction; they were never comparable.

---

## 4. S5

### 4.1 The registered question

> "no misstatement at graph-distance at least 1 from the query moved it (0 of 791)"

### 4.2 ρ = 1, ARM E (exhaustive statement sweep) — the primary arm after M3

Operator: draw an SCM; draw a base `K` of `k_base ∈ {1,…,4}` **true** statements exactly as
`run_arm1.py:139`; then sweep **every** undirected edge of `C` in **both** orientations, replacing
the base statement on that edge if present, appending otherwise. This removes E1′'s `|K| = 4`
sampling of *which* edges may be misstated — rare far edges are no longer under-sampled — while
keeping the base context that gives the detector its dynamic range.

| ensemble | SCMs | `N_set` | **`D_s5`** | `D_am` | `D_con` | `D_all` | clusters | 3/`D_s5` | 3/n_SCM | design eff. | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `licensed` | 28,000 | **0** | 18,327 | 18,327 | 73,303 | 84,095 | 8,643 | 1.64e-04 | 3.47e-04 | 2.12 | **SUPPORTED-WEAK** |
| `large` | 10,000 | **0** | 23,056 | 23,056 | 31,922 | 36,479 | 6,581 | 1.30e-04 | 4.56e-04 | 3.50 | **SUPPORTED** |
| `k8` | 12,000 | **0** | 2,476 | 2,476 | 61,141 | 66,096 | 532 | 1.21e-03 | 5.64e-03 | 4.65 | **INCONCLUSIVE** |
| `deep` | 7,000 | **0** | 7,959 | 7,959 | 17,052 | 21,365 | 3,213 | 3.77e-04 | 9.34e-04 | 2.48 | **SUPPORTED-WEAK** |
| `deep8` | 2,400 | **0** | 2,261 | 2,261 | 10,542 | 13,050 | 567 | 1.33e-03 | 5.29e-03 | 3.99 | **INCONCLUSIVE** |
| **pooled, uniform-draw only** | | **0** | **43,859** | | | | 15,756 | 6.84e-05 | 1.90e-04 | | |
| **ARM C (`argmax δ`), not pooled with the above** | | **0** | **10,220** | | | | 3,780 | 2.94e-04 | 7.94e-04 | | |

PREREG 7.9: ARM C's `argmax δ(X,Y)` query rule makes `deep`/`deep8` non-comparable to the uniform-draw ensembles **by design**, so no rate and no rule-of-three bound is pooled across the two groups. The all-ensemble figure quoted in the headline is a **count of trials in which no counterexample was found**, not a pooled rate.

* **Positive control on the same code path** (`dmin = 0`, ARM E, misstatements and true statements pooled): `licensed` 3,355/33,015 = 0.1016 · `large` 391/10,761 = 0.0363 · `k8` 218/2,513 = 0.0867 · `deep` 388/4,413 = 0.0879 · `deep8` 71/1,031 = 0.0689.

* **Within-SCM matched contrast (M13)** — restricted to the SCMs that contribute *both* a near and a far amenable trial, so the comparison is on the same graphs: `licensed` n=7,864 SCMs, near 2,971/25,207, far **0**/27,684 · `large` n=3,808 SCMs, near 376/9,783, far **0**/26,009 · `k8` n=452 SCMs, near 218/2,495, far **0**/3,802 · `deep` n=2,772 SCMs, near 377/4,026, far **0**/13,953 · `deep8` n=595 SCMs, near 69/1,004, far **0**/5,078. The near/far gap is not a comparison between two disjoint sets of graphs.

* **A *correct* far statement never moved `O*` either** (the split that would refute S8 and S9 without refuting S5): `N_set(far, true_in_D) = 0` in every arm and ensemble, over 106,406 amenable true-statement far trials.

* **Disconnected statements** (`δ = ∞`, kept out of every `≥ k` bucket per G7): `N_set = 0` over 29,464 amenable ARM E trials. This is a definitional zero — Meek R1–R4 act on adjacent triples and `pcp`/`cn`/`forb`/`pa` live in the query's component — and it is reported separately rather than folded into the far denominator.

* **ARM A** (base `K` = ∅, `G0 = C`, demoted to a mechanism arm per M2): `N_set = 0` and `N_ident = 0` at *every* distance including `dmin = 0`, over 82,548 amenable trials. Its zero carries no information about S5 and is excluded from the pooled denominator, exactly as the audit required. Its content is §7.

* **Operator split (`rep` = reverse a statement in the base — the campaign's own operator; `add` = assert a new required edge on an untouched edge — never run by any prior experiment).** Misstatements only: `licensed` — `add` near 0/2,858 = 0.0000, far **0**/8,425, `rep` near 3,355/7,305 = 0.4593, far **0**/9,902 · `large` — `add` near 0/1,920 = 0.0000, far **0**/13,126, `rep` near 391/1,922 = 0.2034, far **0**/9,930 · `k8` — `add` near 0/328 = 0.0000, far **0**/1,853, `rep` near 218/297 = 0.7340, far **0**/623 · `deep` — `add` near 0/101 = 0.0000, far **0**/3,816, `rep` near 388/498 = 0.7791, far **0**/4,143 · `deep8` — `add` near 0/36 = 0.0000, far **0**/1,602, `rep` near 71/83 = 0.8554, far **0**/659. Both operators give zero at `dmin ≥ 1`; the `add` operator also gives zero at `dmin = 0`, which §7 shows is structural rather than a detector failure.

### 4.3 The `dmax` reading, and the joint distribution (M10)

Two joints were asked for and both are reported. The **(distance × event)** joint is the body of
every table in this file — each stratum carries all six events with their own denominators, never a
single marginal rate. The **(`dmin` × `dmax`)** joint has a closed form, and X2 verified it rather
than assuming it:

> **Lemma (checked, 0 violations over 36,616 skeleton edges at p ∈ [4,13], deg ∈ [1,6]).** For any
> edge `{u,v}` of the skeleton, `|δ(u) − δ(v)| ≤ 1`, and `δ(u)` is finite iff `δ(v)` is (machinery
> test T2). Hence **`dmax ∈ {dmin, dmin + 1}` always**, and the (`dmin`, `dmax`) joint is supported
> on exactly four classes:
>
> | class | meaning | observed share (3,000 random CPDAGs) |
> |---|---|---|
> | `(0, 0)` | the statement **is** the X–Y edge | 1,229 |
> | `(0, 1)` | exactly one endpoint is X or Y | 14,287 |
> | `(d, d)` or `(d, d+1)`, `d ≥ 1` | neither endpoint is the query | 19,937 |
> | `(∞, ∞)` | the edge is in another component | 1,063 |
>
> So the `dmin` marginal plus the `dmax` marginal **is** the joint, up to one bit, and no
> cross-tabulation can hide anything. That is why the two conventions can be reported side by side.

"The graph-distance of a statement" is nevertheless undefined for an edge. Under `dmax` (a
statement is *near* only if **both** endpoints are the query, i.e. it *is* the X–Y edge) the
arithmetic inverts completely:

| ensemble (ARM E) | `N_set` at `dmax = 0` | `D_am` at `dmax = 0` | `N_set` at `dmax ≥ 1` | `D_am` at `dmax ≥ 1` |
|---|---|---|---|---|
| `licensed` | 0 | 4,086 | **3,355** | 72,119 |
| `large` | 0 | 995 | **391** | 59,444 |
| `k8` | 0 | 209 | **218** | 8,304 |
| `deep` | 0 | 18 | **388** | 23,693 |
| `deep8` | 0 | 0 | **71** | 6,916 |

On the archive the same inversion is exact: `N_set(dmax ≥ 1)` equals the **total** `N_set` in every
one of the six cells (1,007 / 906 / 208 / 69 / 280 / 598), and `N_set(dmax = 0)` is 0 with
`D_am(dmax = 0) = 0` — misstating the X–Y edge itself **never** moves `O*`; it destroys amenability
every time (`N_ident = D_con` in that row).

**So the defensible sentence is not about distance at all.** It is:

> *no misstatement on an edge **not incident to X or Y** changed `O*` at ρ = 1; every observed
> change involved an edge incident to the query, and misstating the X–Y edge itself never moved
> `O*` — it destroyed identification.*

### 4.4 ρ ≥ 2 — 🔴 THE FALSIFIER FIRED, 70 TIMES

Restricted to members where **every** flipped statement is at finite `dmin ≥ 1` (M18's "all flips
far" stratum), from the E1′ archive, recomputed and independently verified by gate G3.

⚠️ Read the operator before reading the table. The base `K` is built from `D`'s own orientations
(`run_arm1.py:139`), so **every unflipped statement is TRUE**. A member in this stratum therefore
contains *only* misstatements, *all* of them at distance ≥ 1 from the query. There is no near
misstatement anywhere in the knowledge set doing the work. This is precisely the configuration S5
says cannot move `O*`.

| cell | ρ=2 `N_set`/`D_am` | ρ=3 `N_set`/`D_am` | ρ=4 `N_set`/`D_am` |
|---|---|---|---|
| `original`/K4 | **14** / 842 | **1** / 266 | 0 / 26 |
| `licensed`/K4 | **7** / 1,305 | 0 / 479 | 0 / 70 |
| `large`/K4 | **1** / 3,958 | 0 / 1,498 | 0 / 221 |
| `k8`/K4 | **2** / 127 | 0 / 61 | — |
| `k8`/K6 | **12** / 328 | **2** / 219 | — |
| `k8`/K8 | **25** / 698 | **6** / 517 | — |
| **pooled** | **61** / 7,258 | **9** / 3,040 | **0** / 317 |

**70 of 10,615** amenable all-far perturbations changed `O*`. **All 70 are `E_bias`** — the new
`O*` is *invalid* in the true DAG, so all 70 are silent bias, not benign relabelling. Distribution
by `p`: 5→20, 6→14, 7→12, 8→14, 9→7, 10→2, 20→1. Every one of the 70 has **all flipped statements at
`dmin` exactly 1**.

**And the radius is sharp.** Pooling the *entire* ball (ρ = 1…4) over all six cells, restricted to
knowledge sets in which **every** misstatement is at finite `dmin ≥ 1`, then stratified by the
*closest* misstatement:

| closest misstatement at `dmin` = | `D_all` | `D_con` | `D_am` | `N_set` | `N_bias` | `N_ident` |
|---|---|---|---|---|---|---|
| **1** | 31,195 | 15,259 | 14,279 | **70** | **70** | **980** |
| **2** | 6,673 | 4,271 | 4,271 | **0** | **0** | **0** |
| 3 | 2,071 | 1,482 | 1,482 | 0 | 0 | 0 |
| 4 | 855 | 650 | 650 | 0 | 0 | 0 |
| 5 | 253 | 191 | 191 | 0 | 0 | 0 |
| 6 | 116 | 86 | 86 | 0 | 0 | 0 |
| 7 | 42 | 31 | 31 | 0 | 0 | 0 |
| 8 | 14 | 10 | 10 | 0 | 0 | 0 |
| 9 | 4 | 4 | 4 | 0 | 0 | 0 |
| 10 | 1 | 1 | 1 | 0 | 0 | 0 |
| **≥ 2 pooled** | **10,029** | **6,726** | **6,726** | **0** | **0** | **0** |
| **all far** | 41,224 | 21,985 | 21,005 | 70 | 70 | 980 |

**Every event in the entire corpus — silent bias and loud identification loss alike — involves at
least one misstatement at distance exactly 1.** Beyond one hop, over 6,726 amenable perturbations
of arbitrary order, nothing happened at all: `N_set = N_bias = N_ident = 0`. For the *reversal*
operator the elicitation radius is therefore **exactly 1**, not 0 and not unbounded — which is a
sharper and far more useful statement than either S5 or S9 currently makes.

**Minimality (PREREG 4.5.7).** The smallest instance is **p = 5**. An exhaustive ρ = 2 far-pair
search over every CPDAG and every reachable MPDAG at **p = 4** returns **0 hits of 0 eligible
trials** — the configuration (an amenable base with two undirected edges both at distance ≥ 1)
*cannot exist* at p = 4. So p = 5 is minimal, not merely smallest-found.

**Exhibit 1 of 3** (all three in `results/exhibit.json`, regeneration verified exact —
`regen_ok = {xy: True, tau: True, K: True}`):

```
ensemble k8 / arm K6 · seed 393 · p = 5 · deg = 6.0
D : 1->0 2->0 3->0 4->0 1->2 3->1 1->4 3->2 4->2 3->4
C : 0--1 0--2 0--3 0--4 1--2 1--3 1--4 2--3 2--4 3--4      (fully undirected)
X = 1 , Y = 2 , tau = -0.684396
delta : {0:1, 1:0, 2:0, 3:1, 4:1}                           (skeleton BFS from {X,Y})
K  = [3->4, 1->4, 3->1, 3->0, 4->2, 4->0]                   (all TRUE in D)
flipped: statements 3 and 5, i.e. 3->0 becomes 0->3 and 4->0 becomes 0->4
        dmin = [1, 1]   dmax = [1, 1]     — neither endpoint is X or Y
O*(G0) = {3}      beta = -0.684396   valid in D = True
O*(G') = {0,3}    beta = -0.538350   valid in D = FALSE    |bias|/|tau| = 0.2134

Meek chain that carried the misstatement into pa(cn)\forb:
  assert 3->4   propagated: (none)
  assert 1->4   propagated: (none)
  assert 3->1   propagated: (none)
  assert 0->3   propagated: R2([3]) => 0->1 ;  R2([1]) => 0->4
  assert 4->2   propagated: R2([4]) => 0->2 ; R2([4]) => 1->2 ; R2([1]) => 3->2
  assert 0->4   propagated: (none)
```

The mechanism is exactly the one R1 (theorist) named as the unmet proof obligation: **Meek R2 fired
twice from a statement whose two endpoints are both at distance 1, orienting `0->1` and `0->4`, and
that put node 0 into `pa(cn(X,Y))`.** Closure is the propagation channel, and it reaches the query
neighbourhood from outside it.

Two further exhibits (`k8`/K8 seed 4 and seed 59, both p = 5) are in `results/exhibit.json` with the
same structure; seed 59 has `|bias|/|τ| = 1.2723` and flips the **sign** of the estimate
(β: −0.691 → +0.188).

### 4.5 Verdict on S5

| statement | verdict |
|---|---|
| S5 **as printed** in `C-triage.tex:13` (no ρ qualifier, "graph-distance ≥ 1", "0 of 791") | 🔴 **REFUTED.** 70 counterexamples, all silently biasing, minimal p = 5, exhibited with the Meek chain. The number 791 is also not reproducible from any artefact and was ordered retired on 2026-08-14. |
| S5 **restricted to ρ = 1** and to *"not incident to X or Y"* | **SUPPORTED** on `large`; **SUPPORTED-WEAK** on `licensed`, `deep`; **INCONCLUSIVE** (denominator too small) on `k8`, `deep8`. Strengthened out of all proportion to the sampled arms by the **complete enumeration at p ≤ 5** (§7), which has no sampling error at all. |
| S5 as a **property** ("`O*` is a local functional of the graph") | **not licensed by any measurement.** It requires the theorem, which X2 does not discharge (§10.1). |

**Pre-committed mandated wording (PREREG 4.1), now instantiated:**

> *A single orientation error on an edge not incident to X or Y never changed `O*` — 0 of 54,079
> misstatement trials across five ensembles, and 0 of 136,200 in the complete enumeration of every
> CPDAG, every knowledge state and every query at p ≤ 5 — despite Meek closure orienting additional
> edges in 42–67 % of those trials. **Joint errors did**: 70 of 10,615 knowledge sets whose
> misstatements all lie at distance ≥ 1 moved `O*`, and all 70 produced an invalid adjustment set.
> The radius is sharp: in all 70 the closest misstatement is at distance exactly 1, and knowledge
> sets whose closest misstatement is at distance ≥ 2 produced no change of any kind — 0 silent bias
> and 0 identification loss in 6,726 amenable perturbations of any order.*

---

## 5. S8 — the pruning claim

ARM D is **not an unrun experiment**; per M4 it is a re-partition of an already-enumerated ball and
is reported as a **measurement recovered from archived data**, with no prediction claimed. It was
computed twice — once from E1′'s stored estimates (interval version `Q_r`) and once from `O*` sets
re-derived by gate G3 (set version `Qset_r`, new in X2).

`Q_r` = P(the ignorance interval computed on the pruned ball differs from the full ball).
`Qset_r` = P(the *set of reachable `O*`* differs). `Cost_r` = mean over SCMs of `|B^r_ρ| / |B_ρ|`,
with the per-SCM fraction below 0.5 also printed (M5). Clopper–Pearson 95 % intervals (M19).

| cell | n | `Q₀` [CI] | `Qset₀` | `Q₁` [CI] | `Qset₁` | `Q₂` | mean `Cost₀` | mean `Cost₁` | mean `Cost₂` |
|---|---|---|---|---|---|---|---|---|---|
| `original`/K4 | 2,637 | 0.1244 [0.1120,0.1376] | 0.1566 | 0.0011 [0.00023,0.00332] | 0.0011 | 0.0000 | 0.453 | 0.827 | 0.898 |
| `licensed`/K4 | 2,087 | 0.1452 [0.1303,0.1610] | 0.1979 | 0.0014 [0.00030,0.00420] | 0.0014 | 0.0000 | 0.356 | 0.808 | 0.899 |
| `large`/K4 | 2,443 | 0.0143 [0.0100,0.0199] | 0.0156 | 0.0000 [0,0.00151] | 0.0000 | 0.0000 | 0.140 | 0.281 | 0.490 |
| `k8`/K4 † | 127 | 0.1969 [0.1316,0.2767] | 0.2441 | 0.0000 [0,0.02863] | 0.0000 | 0.0000 | 0.385 | 0.921 | 0.980 |
| `k8`/K6 | 329 | 0.4043 [0.3508,0.4595] | 0.5653 | 0.0000 [0,0.01115] | 0.0000 | 0.0000 | 0.348 | 0.948 | 0.986 |
| `k8`/K8 | 585 | 0.5231 [0.4817,0.5642] | 0.7829 | 0.0000 [0,0.00629] | 0.0000 | 0.0000 | 0.308 | 0.966 | 0.991 |

† `k8`/K4 has n = 127 < 300 and is **excluded from every verdict** per M19; it is printed, not used.
`Cost₋₁ = 1/16 = 0.0625` exactly (arithmetic, not a measurement) and `Q₋₁` = 0.0970–0.9009.

**The (Q, Cost) Pareto frontier is the reportable object, and it has no good point:**

* `r = 0` — cheap (`Cost₀` 0.14–0.45) and **unsafe**: the Clopper–Pearson lower bound on `Q₀` is
  strictly above the 0.01 REFUTED bar on `original` (0.1120), `licensed` (0.1303), `k8`/K6 (0.3508)
  and `k8`/K8 (0.4817) — **4 of 5** usable cells. `large` is INCONCLUSIVE by a hair (`Q₀` = 0.0143,
  CI [0.009999, 0.019869]; the point estimate is above the bar, the interval straddles it by 10⁻⁵),
  and I report it as INCONCLUSIVE rather than rounding it into the refutation. PREREG 4.2's rule is
  *"REFUTED iff `Q_r ≥ 0.01` on **any** ensemble"*, so **S8 is REFUTED at r = 0.**
* `r = 1` — safe on `large`/`k8` but the CI on `original`/`licensed` (`Q₁` = 0.0011/0.0014) **spans**
  the 0.001 SUPPORTED bar, so **INCONCLUSIVE** by M19; and the cost half fails outright on
  `original` (0.827), `licensed` (0.808), `k8`/K6 (0.948), `k8`/K8 (0.966). A radius that is safe
  and prunes 3 % of the ball is not a saving.
* `r ≥ 2` — `Q` and `Qset` are 0 everywhere, but `Cost₂` = 0.49–0.99: pruning at radius 2 discards
  almost nothing.

`Qset_r` is **strictly larger** than `Q_r` at r = 0 in every cell (0.157 vs 0.124; 0.198 vs 0.145;
0.565 vs 0.404; 0.783 vs 0.523). Two members can reach different `O*` and still bracket the same
interval, so the interval version *understates* how much pruning changes the analysis.

**Verdict: S8 NOT SUPPORTED AT ANY RADIUS** — REFUTED at `r = 0`, INCONCLUSIVE at `r ≥ 1`. The
sentence cannot be repaired by choosing a radius. Separately, and independently of these numbers,
X2 tests only the **validity** of pruning; no pruned enumerator exists and no wall-clock was
measured, so *"cost does not scale with |K|"* remains an **unimplemented, unbenchmarked** claim
(§10.2).

---

## 6. S9 — the elicitation radius

S9 is not S5. A statement that destroys or creates identification has changed the analysis
materially without ever moving `O*`. Per M6, the registered statistic is the **radius**, not a
threshold count.

| cell | `r*_set` | `r*_ident` | `r*_gain` | **`r*_elicit`** | `N_gain` by `dmin` | `N_ident` by `dmin` | composite near/far |
|---|---|---|---|---|---|---|---|
| `licensed` / ARM A | -1 | -1 | 4 | **4** | 1:9,298 2:563 3:21 ge4:2 disc:0 | 1:0 2:0 3:0 ge4:0 disc:0 | 0.0974 / 0.0588 |
| `licensed` / ARM E | 0 | 2 | 3 | **3** | 1:6,702 2:130 3:1 ge4:0 disc:0 | 1:637 2:3 3:0 ge4:0 disc:0 | 0.1294 / 0.0475 |
| `large` / ARM A | -1 | -1 | 4 | **4** | 1:2,356 2:595 3:66 ge4:2 disc:0 | 1:0 2:0 3:0 ge4:0 disc:0 | 0.3397 / 0.0414 |
| `large` / ARM E | 0 | 2 | 3 | **3** | 1:1,193 2:223 3:26 ge4:0 disc:0 | 1:68 2:3 3:0 ge4:0 disc:0 | 0.2940 / 0.0221 |
| `k8` / ARM A | -1 | -1 | 4 | **4** | 1:2,590 2:107 3:8 ge4:1 disc:0 | 1:0 2:0 3:0 ge4:0 disc:0 | 0.0063 / 0.0205 |
| `k8` / ARM E | 0 | 2 | 4 | **4** | 1:3,763 2:56 3:2 ge4:1 disc:0 | 1:227 2:1 3:0 ge4:0 disc:0 | 0.0353 / 0.0318 |
| `deep` / ARM A | -1 | -1 | 2 | **2** | 1:506 2:51 3:0 ge4:0 disc:0 | 1:0 2:0 3:0 ge4:0 disc:0 | 0.3963 / 0.0130 |
| `deep` / ARM E | 0 | 1 | 2 | **2** | 1:401 2:14 3:0 ge4:0 disc:0 | 1:44 2:0 3:0 ge4:0 disc:0 | 0.4351 / 0.0119 |
| `deep8` / ARM A | -1 | -1 | 3 | **3** | 1:261 2:55 3:4 ge4:0 disc:0 | 1:0 2:0 3:0 ge4:0 disc:0 | 0.2891 / 0.0123 |
| `deep8` / ARM E | 0 | 1 | 3 | **3** | 1:341 2:24 3:1 ge4:0 disc:0 | 1:17 2:0 3:0 ge4:0 disc:0 | 0.3646 / 0.0162 |

`r* = -1` means the event never fired. `ge4` is the `dmin ≥ 4` bucket. The composite is `P(E_set ∪ E_ident ∪ E_gain | consistent)` on `D_con` (M14), the quantity S9 is actually about and the one immune to the competing-risk conditioning.

**`E_ident` — the M15 split is qualitative, not cosmetic.** At `dmin = 0` the majority of
identification losses are **`nopath`**: the closure destroys every proper possibly-causal path, so
the analyst does not read "not identified" — they read *"there is no causal effect"*, which is a
**wrong answer**, not a lost one (ARM E, `dmin = 0`: licensed 4031 `nopath` vs 1229 `unamenable`; large 1432 vs 214;
k8 113 vs 251; deep 2595 vs 62; deep8 539 vs 22). **At `dmin ≥ 1` the `nopath` count is 0 in every ensemble and every
stratum** — far misstatements can un-identify a query but they never fabricate a null effect. That
distinction changes what the paper's replacement sentence has to say, and it did not exist before.

**`E_gain` — never measured by any prior experiment, and it is what actually kills S9.**
A single statement at finite distance ≥ 1 turned a *non-identified* query into an identified one
**29,364** times across the five ensembles and both arms. Every single gain is `gain_orient`
(the possibly-causal paths already existed and the far statement oriented the first edge);
`gain_paths` is **0** everywhere, so no far statement ever *created* a causal path. And the
registered integrity surprise — `N_gain(disconnected) ≥ 1`, which would be a bug — is
**0 in every ensemble and every arm**, as it must be. Far `E_ident` fires
**1,000** times, of which **0** are `nopath`.

**The radius depends on the operator, and X2 reports both rather than one.**

| operator | arm | `r*_set` (silent bias) | `r*_ident` (identification lost) | `r*_gain` (identification created) | **`r*_elicit`** |
|---|---|---|---|---|---|
| **reverse** an asserted orientation, whole ball ρ ≤ 4 | ARM B archive (§4.4) | **1** | **1** | n/a (base always amenable) | **1** |
| reverse, ρ = 1 only | ARM B archive (§3) | **0** | 1 | n/a | 1 |
| **add** a required edge on an untouched edge, plus reverse | ARM A + ARM E | **0** | **2** | **≥ 4** | **≥ 4** |

**Verdict: S9 REFUTED.** Under the reversal operator the elicitation radius is **exactly 1** —
sharp, and far more useful than S5/S9's "0". Under the *addition* operator, which no prior
experiment ran, it is **≥ 4**, driven entirely by `E_gain`. Either way S9-as-written (`r*_elicit = 0`)
is false. All radii are skeleton distances and are therefore **lower** bounds on the
orientation-respecting radius (M11: `d_skel ≤ d_dir` makes the skeleton far-set the *smaller* one).

The replacement sentence cannot be *"knowledge outside that neighbourhood need not be elicited"*
under any reading. What is true is narrower and more useful:

> *A misstatement more than one hop from the query never silently biased the estimate — 0 of 6,726
> perturbations of any order — but knowledge that far out, and further, routinely decides whether
> the effect is identified at all, in both directions. Far knowledge need not be elicited to avoid
> bias; it must be elicited to obtain identification.*

---

## 7. The lemma — exhaustive existence search at p ≤ 5, plus a bounded p = 6 sample

M2 asks for an exhaustive small-`p` search as the honest substitute for the theorem. Because
`E_set`, `E_ident` and `E_gain` are functions of `(G0, X, Y, s)` alone — no SCM, no `Σ`, no truth —
the ρ = 1 question can be **enumerated** rather than sampled:

1. every labelled DAG on `p` nodes → CPDAG → dedupe;
2. from each CPDAG, BFS over **every reachable MPDAG** (orient one undirected edge, Meek-close,
   reject a cycle or a new v-structure) — this is exactly `{ M(C,K) : K Meek-consistent }`, i.e.
   *every knowledge state*;
3. every ordered query, every undirected edge of `G0`, both orientations.

Two operators are scored: `add` (assert `s` on an edge the base left open) and `reverse` (compare
the two orientations of the same edge against each other — the pilot's and E1′'s operator,
conditional on the rest of `K`).

| p | mode | CPDAGs | MPDAGs | stratum | `N_set_rev` / `D_am_rev` | `N_set_add` / `D_am_add` |
|---|---|---|---|---|---|---|
| 4 | **exhaustive** | 185 | 1,601 | `dmin = 0` | **1,272 / 2,304 = 0.552** | 0 / 2,016 |
| | | | | `dmin ≥ 1` | **0 / 480** | 0 / 912 |
| | | | | disconnected | 0 / 12 | 0 / 24 |
| 5 | **exhaustive** | 8,782 | 116,992 | `dmin = 0` | **183,060 / 305,280 = 0.600** | 0 / 233,280 |
| | | | | `dmin = 1` | **0 / 132,120** | 0 / 252,720 |
| | | | | `dmin = 2` | **0 / 3,600** | 0 / 7,200 |
| | | | | disconnected | 0 / 1,320 | 0 / 2,640 |

**p = 6 is a SAMPLE, not an enumeration** (§10.3): 676 distinct CPDAGs drawn at random, of which
**639 were scanned and 37 skipped whole** by a registered cap of 400 reachable MPDAGs per CPDAG
(a skipped CPDAG is never partially scanned — a partial scan would drop knowledge states
non-uniformly). 16,555 MPDAGs. `reverse`: **34,427 / 59,984 = 0.574** at `dmin = 0`;
**0 / 42,876** at `dmin = 1`; **0 / 2,648** at `dmin = 2`; **0 / 121** at `dmin = 3`;
**0 / 3,056** disconnected. `add`: **0 / 88,912** at every distance.

⚠️ The first p = 6 attempt (2,335 CPDAGs, no cap) was **killed at 20.9 minutes** for breaching the
PREREG 3.1 wall-clock budget, and produced nothing. The capped re-run above finished in ~2 minutes.
The cap is a **selection rule on the ensemble**, so it is stated here rather than hidden in the
code: it excludes the densest 5.5 % of sampled CPDAGs, which are exactly the ones with the most
knowledge states — the p = 6 row is therefore biased *against* finding a counterexample and is
reported as a bounded sample, not as evidence of absence.

**At p ≤ 5 this is a complete enumeration**: 0 of **136,200** far reversal trials moved `O*`,
against **184,332 of 307,584** near ones. There is no sampling error in that zero.

Two structural facts fall out and are worth more than the zero:

* **`E_meek ≡ 0` when the statement is on an edge the MPDAG still leaves undirected.** 0 MeekFails
  in 9,300,976 such single-statement additions at p ≤ 5. Every undirected edge of an MPDAG is
  orientable both ways inside its restricted equivalence class, so such a statement can never
  produce a directed cycle or a new v-structure. The campaign's ~35 % Meek-rejection rate therefore
  does **not** measure "expert statements are often caught"; it measures *contradicting something
  the rest of the knowledge set had already compelled* — directly, by reversing an asserted
  orientation, or indirectly, through the closure of the others. That is a property of the
  perturbation operator, not of expert error.
* **`add` from an amenable base never moves `O*`, at any distance** (0 of 498,792 at p ≤ 5). Every
  reversal event therefore lives in the region where the *common base is not amenable* and the two
  orientations give different identified `O*`. That is the shape of the theorem R1 asked for, and it
  is why ARM A (base = bare CPDAG) is null. The mechanism, measured on every ensemble:

  | ensemble | `amen(C)` rate | mean \|`cn`\| given `amen(C)` | mean undirected edges **incident to `cn`** | share of amenable CPDAGs with ≥ 1 such edge |
  |---|---|---|---|---|
  | `licensed` | 0.1305 | 1.470 | **0.2480** | **9.2 %** |
  | `large` | 0.5016 | 1.470 | **0.0068** | **0.7 %** |
  | `k8` | 0.0132 | 1.799 | **1.1698** | **15.1 %** |
  | `deep` | 0.1233 | 3.294 | **0.0359** | **3.5 %** |
  | `deep8` | 0.0242 | 3.517 | **0.0000** | **0.0 %** |

  Amenability plus Meek closure has already compelled essentially every edge touching the causal-node
  set, so `pa(cn)` is frozen and one further orientation cannot move it — at any distance, near
  included. Note `deep`/`deep8` have the largest `cn` (3.3–3.5 nodes, from the `argmax δ(X,Y)` query
  rule that deliberately lengthens the causal paths) and the *fewest* undirected edges touching it,
  which is the opposite of what PREREG §2.3's "longer causal paths ⇒ locality breaks" mechanism
  predicted. That prediction is refuted, and it is the one place X2's ARM C had a registered
  direction to be wrong about.

---

## 8. S3, S4 and S5 are one contingency table, not three findings (M12)

Per statement, ρ = 1, on the ARM B archive:

| cell | S3 `D_con/D_all` | S4 `N_bias/D_con` | = `P(dmin=0 \| consistent)` × `L_bias(dmin=0)` | `N_bias` at `dmin ≥ 1` |
|---|---|---|---|---|
| `original`/K4 | 0.6424 (6,776/10,548) | **0.1486** (1,007/6,776) | 0.6234 × 0.2384 = 0.1486 | **0** |
| `licensed`/K4 | 0.6763 (5,646/8,348) | 0.1605 (906/5,646) | 0.5177 × 0.3100 = 0.1605 | **0** |
| `large`/K4 | 0.8278 (8,089/9,772) | **0.0257** (208/8,089) | 0.2120 × 0.1213 = 0.0257 | **0** |
| `k8`/K4 | 0.7854 (399/508) | 0.1729 (69/399) | 0.5063 × 0.3416 = 0.1729 | **0** |
| `k8`/K6 | 0.6373 (1,258/1,974) | 0.2226 (280/1,258) | 0.5882 × 0.3784 = 0.2226 | **0** |
| `k8`/K8 | 0.5128 (2,400/4,680) | **0.2492** (598/2,400) | 0.6067 × 0.4107 = 0.2492 | **0** |

The decomposition is **exact to four decimals in all six cells**. S4's "14.1 %" is not a property of
expert error; it is `P(the misstatement was incident to the query)` times `P(it biased | incident)`,
and both factors move with `p` and density — the same quantity is **0.0257** on `large` and
**0.2492** on `k8`/K8, a 9.7× spread. S4's numerator and S5's zero are the two cells of **one**
partition of the same trials. **Citing them as independent support double-counts one measurement.**

Two identities, recomputed rather than assumed (PREREG 5.5, 5.6):

* **`E_bias ≡ E_set`** — 0 discrepancies in all six cells. Every `O*` change was invalid in the true
  DAG. `changed_still_valid = 0` replicates.
* **`E_est(∞) ≡ E_set`** — 0 discrepancies. So *"`O*` changed"*, *"the estimate moved"* and *"silent
  bias"* are **one** measurement with three names, and must not be cited as three pieces of
  evidence.

---

## 9. What the closure was doing while nothing happened (M18)

A bare zero invites "the perturbation did nothing". It did plenty. `P(cascade > 0 | consistent)` —
the fraction of far trials in which Meek closure oriented **additional** edges beyond the asserted
one:

| ensemble (ARM E, consistent trials) | `P(cascade>0)` at `dmin = 0` | `P(cascade>0)` at `dmin ≥ 1` |
|---|---|---|
| `licensed` | 0.5379 | **0.5395** |
| `large` | 0.5061 | **0.4657** |
| `k8` | 0.5902 | **0.6036** |
| `deep` | 0.3781 | **0.4220** |
| `deep8` | 0.6464 | **0.6730** |

On the archive, at `dmin ≥ 2`: `original` 74/285 = 0.26, `licensed` 103/334 = 0.31, `large`
750/1,624 = 0.46, `k8`/K4 10/13 = 0.77. So the honest sentence is *"0 of 54,079 **despite the
closure orienting additional edges in a quarter to two-thirds of those trials**"* — which is a
substantive claim about where propagation stops, not an observation that nothing propagated.

---

## 10. Limitations — in writing, so a gap cannot later read as coverage

1. 🔴 **This is not a proof.** An empirical zero, even an exhaustive one at p ≤ 5, is not the theorem
   R1 asked for (*"no distant statement can enter the neighbourhood through Meek closure"*). And at
   ρ ≥ 2 the theorem is **false** — §4.4 exhibits the counterexample. What X2 licenses is
   *"checked at ρ = 1; exhaustively true at p ≤ 5; refuted at ρ ≥ 2 from p = 5 upward."*
2. 🔴 **S8 is only half-tested.** X2 measures the **validity** of pruning (`Q_r`, `Qset_r`), never
   its cost. No pruned enumerator was implemented; `Cost_r` is the combinatorial ratio
   `|B^r_ρ|/|B_ρ|`, not a timing. "Cost does not scale with |K|" is unbenchmarked even where
   pruning is safe.
3. <a name="s103"></a>**The p = 6 lemma is a SAMPLE, not an enumeration, and a *capped* one** — a
   deviation from M2's wording. There are 3.78 M labelled DAGs on 6 nodes and step 3 would be ~10⁸
   `O*` evaluations. The uncapped 2,335-CPDAG attempt was killed at 20.9 min for breaching the
   PREREG 3.1 budget; the reported run caps reachable MPDAGs at 400 per CPDAG and skips 37 of 676
   whole. Those 37 are the densest CPDAGs — the ones with the most knowledge states — so the p = 6
   row is biased *against* finding a counterexample. Labelled as a bounded sample everywhere it
   appears, and it carries no weight the p ≤ 5 enumeration does not already carry.
4. 🔴 **The X1 class is excluded and can invalidate every distance in this file.** `graphs.py:192`
   (**not** `:179-180`, a citation slip in the PREREG that the AUDIT caught) refuses a statement on
   a non-adjacent pair, so "required edges asserted between non-adjacent variables" (S7) cannot be
   perturbed here. This is not merely uncovered: **a spurious edge adds to the skeleton and
   therefore changes `δ(·)` itself.** Every X2 distance assumes the skeleton is correct. Registered
   dependency on X1, not a footnote.
5. **The CPDAG is the oracle's.** `C = dag_to_cpdag(D)`. An analyst runs PC/GES and gets `Ĉ ≠ C`.
   Distances, admissibility and `U(C)` are all defined against the true CPDAG.
6. **Oracle leakage in the query draw.** `run_arm1.py:127-129` selects `(X,Y)` from pairs with a
   possibly-causal path *in `D`*, and the uniform-draw ensembles inherit it for comparability with
   E1′. ARM C (`deep`, `deep8`) defines its candidate set on `C` instead, per M16. The measured
   overlap: the `D`-based set is a **strict subset** of the `C`-based one in every SCM
   (`n_both = n_D` exactly), and the `C`-based set is ~1.5–1.9× larger (mean |pairs_D| 29.6 vs
   |pairs_C| 44.3 on `deep`; 33.0 vs 61.5 on `deep8`). So M16's fix is not cosmetic: ARM C actually
   selects its query from a wider, non-oracle set. `τ` floor: 1e-6 at draw, 1e-3 as the reported
   sensitivity (M16).
7. **`r*_elicit` is measured on the skeleton and is a LOWER bound** (M11). `d_skel ≤ d_dir`, so
   `{d_skel ≥ 1} ⊆ {d_dir ≥ 1}`: the skeleton far-set is the smaller one and the skeleton claim is
   the **weaker** one. The PREREG asserted the opposite; it is wrong and it is corrected here.
8. **`ρ > 4` is not explored.** The ball is capped at the pilot's `MAX_K = 4` (3 on `k8`). A pruning
   failure first appearing at ρ = 5 would not be seen. Given that ρ = 2 already refutes S5, this
   matters only for the *rate*, not the existence.
9. **The competing-risk conditioning is real** (AUDIT A12). `D_am` excludes trials that destroyed
   amenability, and that hazard is itself a function of distance (`D_am/D_con` = 0.60–0.74 at
   `dmin = 0`, **1.000** at `dmin ≥ 2` in every cell). `L_set` alone cannot separate "far statements
   do nothing" from "far statements are filtered differently", which is why the composite
   `P(E_set ∪ E_ident ∪ E_gain | consistent)` is reported beside it as co-primary (§6).
10. **Scope.** Linear-Gaussian iSCM, ER DAGs `p ∈ [5,25]` (plus complete enumeration at `p ≤ 5`),
    expected degree ≤ 6, single-node `X`, single-node `Y`, orientation statements only, population
    `Σ` (`n = ∞` decisive, `n = 2×10⁴` reported). Not licensed: nonlinear mechanisms, latent
    confounding, multi-node interventions, tiered knowledge, real or semi-synthetic data.
11. **No p-value.** `Q_r` and `L_set` are worst-case-over-a-ball sensitivity quantities in the
    Rosenbaum-Γ tradition. Read as tests they are uncorrected multiplicity. They are not tests.
12. **ARM C is not poolable** with the uniform-draw ensembles by design (`argmax δ(X,Y)` query rule)
    and is not pooled anywhere.

### 10.a Deviations from the pre-registration, disclosed

* **No ensemble was enlarged after seeing its result.** The first pass realised a lower per-SCM
  `D_s5` yield than the pilot projected, which leaves `licensed` short of M7's 20,000 SUPPORTED
  threshold. I considered topping it up — M7 explicitly budgets for 20,000 — and **decided against
  it**, because "choose `N` after seeing `n`" is the objection AUDIT A5 raises and the exhaustive
  `p ≤ 5` enumeration (§7) already provides a stronger zero than any sampled SUPPORTED could. The
  threshold stands where it was written and `licensed` is reported at whatever verdict its realised
  denominator earns.
* **The `op` (`add` / `rep`) split was added after the first pass** and required re-running the
  scans at the **identical seed and identical `n_SCM`**. Pass 1 is kept verbatim in
  `results/pass1/` on betelgeuse.
* 🔴 **The scan is not bit-reproducible, and re-running it measured by how much.** `run_scan.py`
  collects results with `Pool.imap_unordered` and stops at the first `n_SCM` *keepers*, so **which**
  keepers make the cut depends on worker scheduling. (E1′'s `run_arm1.py` uses ordered `imap` and is
  deterministic; this is my regression, not an inherited one.) Re-running `licensed` and `large` at
  the same seed and the same `n_SCM` gave `D_s5` = 18,327 vs 18,330 and 23,056 vs 23,062 — a drift
  of **3 and 6 trials, 0.02 % and 0.03 %**. The S5 numerator was 0 in both passes. That drift is the
  measured reproducibility of every scan number in this file, and no conclusion here turns on a
  quantity smaller than it.
* **G4b failed** (§2). Reported, not repaired.
* **p = 6 lemma sampled and capped, not enumerated; the first attempt was killed for breaching the 20-minute budget** (§10.3).

---

## 11. Where I disagree with the auditor — implemented anyway, recorded here

1. **A2's demotion of ARM A is right, but its reason understates the finding.** The auditor calls
   ARM A "null by construction". The exhaustive p ≤ 5 enumeration shows the null is *general*: the
   `add` operator moves `O*` **never** from an amenable base, at any distance, including 0. That is
   not a property of ARM A's ensemble; it is a property of the operator. So ARM A is not a wasted
   arm — it is the only arm that isolates the operator whose zero is a lemma rather than a sample.
2. **A5 treats "the SUPPORTED bar is unreachable" as a fatal pre-registration defect.** I implemented
   M7 as written. But raising `n` against a **fixed** threshold with a **zero** numerator is not the
   mirror image of a free bar: a free bar cannot fail, whereas this one can (one counterexample
   flips it to REFUTED at any `n`). The asymmetry of an existence question makes "choose `N` after
   seeing `n`" much less dangerous here than the audit implies. Having said that, I did **not** use
   the freedom: no ensemble was enlarged after its result was known (§10.a).
3. **A15 says G5 "cannot fail" and is therefore not a placebo.** True, and G10 is the better gate.
   But G10 as specified also cannot fail for an *interval* reason — the frozen ball has every `est`
   equal to `est0` by construction, so `Q_r = 0` is arithmetic once the partition indexing is
   correct. Its real content is that it catches a **bookkeeping** bug in the prune-set construction,
   which is exactly what it is for; I ran it (0 of 39,400) and I state its scope rather than
   claiming more.

---

## 12. Vault items this run bears on — flagged, not acted on

* 🔴 `wiki/activities/active-claims.md:570` records b-LOAD's graph-distance falsifier as *"currently
  NOT firing"* on the strength of `0/791`. **It fires.** X2 produces **70** instances, all silently
  biasing, minimum p = 5, fully regenerable — not one borderline case. The b-LOAD CIKM camera-ready
  is **2026-08-23**. Under Operating Rule #1 this entry must be corrected before any locality
  sentence is pasted into that camera-ready. **Zé's decision and Zé's wording; no agent edits a live
  claim.**
* The AUDIT's point stands and is sharpened: the entry **cannot be repaired by substituting a new
  denominator**, because the replacement number is 0 under `dmin` and equal to the *total* event
  count under `dmax`. The claim needs new wording, not a new figure.
* `wiki/projects/rho-breakdown-adjustment.md:204-205, :224` still carries `0/791` twice, once with
  the word *"provably"*, which no measurement of any size can support — and which is now known to be
  false at ρ ≥ 2.
* E1′ `AUDIT.md` A10 remains open: there is no ρ-breakdown entry in `active-claims.md` and
  `grep -rl 'relevance:.*rho-breakdown' wiki/literature/` returns zero notes. X2 produces evidence
  for a claim that is not registered.

---

## 13. Files

```
output/2026-08-19_x1-x2-locality-spurious/x2/
  PREREG.md   AUDIT.md   RESULTS.md (this file)
  code/    x2lib.py  test_machinery.py  step0_archive.py  run_scan.py
           verify_armB.py  run_lemma.py  exhibit.py  analyse_x2.py  run_all3.sh
           graphs.py  adjust.py  scm.py  se.py        (copied UNMODIFIED from e1prime-se)
  results/ machinery_tests.json  step0_archive.json  verify_armB.json
           lemma_p{4,5,6}.json    exhibit.json       x2_analysis.json
  logs/    scan.log  scan3.log  scan4.log  verify_armB.log  lemma_p5.log
           lemma_p6.log  exhibit.log  analyse.log
```

The raw per-SCM scan files (`scan_{licensed,large,k8,deep,deep8}.json`, ~100 MB, plus
`results/pass1/`) stay on betelgeuse at `~/latent-causal/x2-locality/results/` and are **not**
copied into the vault — same convention as the E1′ mirror, which carries only the `*_analysis_*`
files. `x2_analysis.json` contains every aggregate any table in this file uses.

Work dir on betelgeuse: `~/latent-causal/x2-locality/{code,results,logs}`.
**Nothing under `~/latent-causal/e1prime-se/` or `~/latent-causal/rho-breakdown-knowledge/` was
modified. No git commit was made by any agent.**
