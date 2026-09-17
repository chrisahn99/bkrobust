# E3 — O\*-relevance vs edit distance: PRE-REGISTRATION

**Written 2026-08-14 BEFORE any number was computed.** Nothing below was chosen after seeing a
result. The analysis script (`e3_ostar_relevance.py`) was written after this file and its thresholds
are read from it.

---

## 1. The metric, defined (the campaign never wrote it down)

### 1.1 Objects the analyst actually has

An analyst holds a CPDAG `C` (estimated from data), a query `(X, Y)`, and a stated set of orientation
constraints `K` elicited from an expert. Meek closure gives the MPDAG `G0 = closure(C, K)` and the
optimal adjustment set

```
O*(G0, X, Y) = pa( cn(X,Y,G0), G0 ) \ forb(X,Y,G0)
```

(Henckel/Perković/Maathuis 1907.02435). The analyst does **not** know which statements of `K` are
wrong. A useful metric must therefore be computable from `(C, K, X, Y)` alone — **no true DAG, and
no enumeration of the perturbation ball.** A "metric" that needs the perturbed graph is not a
cheaper alternative to ρ\*; it *is* ρ\*.

### 1.2 The O\*-defining region

`O*` is a function of exactly three node sets computed on `G0`:

- `cn  = cn(X,Y,G0)`   — nodes on proper possibly-causal paths X→…→Y, excluding X
- `pa  = pa(cn, G0)`   — their parents (the set O\* is drawn *from*)
- `frb = forb(X,Y,G0)` — possible descendants of `cn`, plus X (the set O\* is *reduced by*)

Define the **O\*-defining region**

```
R(G0,X,Y) := {X, Y} ∪ cn ∪ pa ∪ frb
```

Every node whose incident-edge orientations can change the *value* of the expression
`pa(cn) \ forb` lies in `R`. Nodes outside `R` enter the formula nowhere.

### 1.3 Relevance of a single knowledge statement

For `k = (u → v) ∈ K`:

```
r(k)  :=  1  if {u,v} ∩ R ≠ ∅   else 0                          (binary relevance)

r2(k) :=  2  if {u,v} ∩ (cn ∪ {X})        ≠ ∅                    (graded relevance)
          1  elif {u,v} ∩ (pa ∪ frb)      ≠ ∅
          0  otherwise
```

**Justification of each term, and why there are no free weights.** The two tiers of `r2` are not
tuned: they are the two syntactic roles a node can play in `pa(cn) \ forb`. Tier 2 is the
*causal-path core* — a node in `cn ∪ {X}` sits on the possibly-causal paths themselves, so
mis-orienting an edge at it can destroy **amenability**, i.e. remove the estimand entirely (a
categorical failure). Tier 1 is the *adjustment machinery* — a node in `pa ∪ frb` can enter or
leave `O*` without touching amenability (a quantitative failure). Tier 0 nodes appear in no term of
the formula. The ordinal ranking 2 > 1 > 0 is fixed a priori by "destroys the estimand" >
"changes the estimate" > "appears nowhere". No coefficient is fitted anywhere in this experiment.

### 1.4 The metric on a knowledge error

For a set `F ⊆ K` of statements the expert got backwards:

| name | definition | what it says |
|---|---|---|
| `d_edit(F)` | `|F|` (= ρ) | **BASELINE.** Count the wrong orientations. The currency the campaign and every located competitor uses. |
| `d_ostar(F)` | `Σ_{k∈F} r(k)` | **CANDIDATE.** Count the *relevant* wrong orientations. |
| `d_ostar_graded(F)` | `Σ_{k∈F} r2(k)` | **CANDIDATE.** Weight them by structural role. |
| `d_ostar_any(F)` | `1[∃k∈F: r(k)=1]` | **CANDIDATE (binary).** Did the expert get *any* relevant statement wrong? |

### 1.5 The three comparators that decide whether this is a section or a sentence

| name | definition | why it is here |
|---|---|---|
| `d_edit` | `|F|` | the metric to beat |
| `d_xy_inc` | `1[∃k∈F : {u,v} ∩ {X,Y} ≠ ∅]` | **the honest-risk comparator.** The pilot already showed all silent bias sits at graph-distance 0. If `d_ostar` only re-encodes "is the statement incident to X or Y", the contribution is one sentence, not a section. |
| `d_gdist` | `min_{k∈F} dist(k)`, `dist(k) = min hop-distance in skeleton(C) from {u,v} to {X,Y}` | the continuous version of the same risk |

---

## 2. Data, population, outcome

**Data.** `/Users/josecosta/research-pilots/latent-causal-rho-breakdown-knowledge/results/linear_raw.json`
— 600 SCMs, 8,085 perturbation records, read-only (copied, never modified). `D` and `C` are not
stored and are regenerated bitwise from the stored `(seed, p, deg)` triple. **Integrity gate:** the
run aborts unless the regenerated `n_undirected`, `n_edges` and `O0` match the stored values for
every SCM, and unless the flip-set enumeration order reproduces every member's stored `rho`.

**Populations** (reported in this order; **P3 is primary**):

- **P1 — the whole ball.** Every member. Event = silent bias.
- **P2 — consistent.** Members surviving the free Meek-consistency check.
- **P3 — undetectable (PRIMARY).** Members that are Meek-consistent *and* amenable: the analyst
  runs the pipeline, gets a number, and nothing looks wrong. Event = `ostar_valid == False`
  ("silent bias"). This is the only population where the metric has a decision to inform.

**Outcomes.**

- **Binary:** `silent = consistent ∧ amenable ∧ ¬ostar_valid`. Scored by **AUC** (Mann-Whitney).
- **Continuous:** `|bias| / |tau|` on P3. Scored by **Spearman ρ**.

**Uncertainty.** Members are clustered within SCM, so all CIs come from a **cluster bootstrap
resampling the 600 SCMs** (1,000 replicates, percentile CIs). ΔAUC CIs are computed **paired** on
the same resample. The primary arm is `generic` (n = 579 SCMs); `tiered` is reported separately and
is not used for the verdict.

---

## 3. PRE-REGISTERED DECISION RULE

Let `A_ost` = AUC of `d_ostar_graded` on P3, `A_edit` = AUC of `d_edit`, `A_xy` = AUC of
`d_xy_inc`. (`d_ostar_graded` is nominated as the primary candidate *before* running, because it is
the metric the campaign's abstract describes — "an error outside O\*'s neighbourhood costs exactly
zero" — in its most informative form. `d_ostar` and `d_ostar_any` are reported as secondaries.)

**SUPPORTED — worth a paper section (§4 of #1).** All three must hold:
1. `A_ost ≥ 0.70` (conventional floor for decision-usable discrimination, Hosmer–Lemeshow);
2. `A_ost − A_edit ≥ 0.10`, lower bound of the paired cluster-bootstrap 95% CI > 0;
3. `A_ost − A_xy ≥ 0.05`, lower bound of the paired cluster-bootstrap 95% CI > 0.

**PARTIAL / ABSORBED — worth one sentence, not a section.** (1) and (2) hold but (3) fails.
Reported verdict: *the mechanism is graph-locality; the O\*-relevance metric is a repackaging of
"is this statement incident to the query", and the contribution is one sentence in #1's §4.*

**REFUTED as a contribution.** (1) fails, or (2) fails.

**Why 0.10 and 0.05.** 0.10 is the smallest AUC gain over simple error-counting that could support
the sentence "edit distance is the wrong yardstick" — a 0.02–0.05 gain is a tweak, not a change of
metric. 0.05 is roughly 2.5–5 cluster-bootstrap standard errors at n = 579 SCMs (expected SE
0.01–0.02), i.e. the smallest gap that would still be visible after the CI is drawn. Both are fixed
here, before any number exists.

**C13 compliance (the vault's standing rule: a null must show its instrument could have produced a
non-null).** The run reports the AUC of the **oracle feature** `ostar_changed` (a post-perturbation
quantity no analyst has) on the same population with the same code path. If that oracle does not
clear 0.90, the AUC machinery lacks range on this outcome and **no null below is interpretable**.
The run also reports the marginal distribution of `r(k)` over all statements: if `r(k) = 1` for
(essentially) every statement, `d_ostar ≡ d_edit` by construction and the comparison is vacuous —
that too is a reportable degeneracy, not a null.

---

## 4. What each result would mean

- **SUPPORTED** → novelty claim (iii) survives as #1's §4 with a measured head-to-head, and the
  paper can state that the estimand induces a metric that dominates edit distance.
- **PARTIAL** → the campaign's #4 abstract is over-claimed: "we define an O\*-relevance metric" must
  become "silent bias is confined to statements incident to the query", one sentence and one figure.
- **REFUTED** → the O\*-relevance metric is dropped; #1 rests on legs (1) perturbed-object and
  (3) sample-size-indexed ρ\* only, and the ICLR abstract must not promise a new metric.
