# X3 — skeleton perturbation: does the locality law survive when the DISCOVERY is wrong?

**Written 2026-08-31, before the first number.** Prompted by the second research question posted to
the Teams channel for the Thursday meeting: *« l'assumption que la data driven method est correcte
est très forte (et très contestable auprès des chercheurs) »*.

## 1. The question, and why it is not answered yet

Every measurement in this line uses the **oracle CPDAG**. The locality law
([[../2026-08-19_e1prime-se-criterion/RESULTS.md]]) says the breakdown radius is informative only
when an elicited statement touches the query neighbourhood: censoring 0.037–0.070 at hop 0 against
0.935–1.000 at hop ≥ 1, a 13–27× separation invariant in `n`.

That was measured for **orientation** statements. The claim that would answer the research question
is about **adjacency** errors, and it has never been measured:

> **H1 (the inference under test).** Damage from a skeleton error obeys the same locality law, because
> `O* = pa(cn(x,y)) \ forb(x,y)` is a local functional of the graph and does not read structure far
> from `{x, y}`.

If H1 holds, the answer to the research question becomes quantitative and much weaker than "the
discovery method must be correct": **the discovery method only has to be correct within `k` hops of
the query, and `k` is measured.** If H1 fails, the framework inherits the full strength of the
assumption and that must be written in the paper.

## 2. Model of a skeleton error — stated as an assumption, not derived

At finite samples, particularly `n ≲ d`, CI tests fail chiefly by **missing weak edges**, and a
missing edge is exactly what breaks back-door blocking. We model one such failure as:

> The analyst's graph is `Ĉ = dag_to_cpdag(D')` where `D' = D` with one edge deleted.

⚠️ **This is a modelling choice with a known gap.** A CI-test failure does not in general return the
CPDAG of a sub-DAG; it returns whatever the algorithm's constraint set implies, which may not be the
CPDAG of any DAG missing exactly that edge. The choice buys a always-valid CPDAG and an
interpretable single-edge perturbation. It is stated here so that no result is read as more general
than it is. **Edge ADDITION is not tested** (the failure direction at `n ≲ d` is omission); that is a
declared scope limit, not an oversight.

## 3. Design

For each SCM: draw `D = random_dag(p, deg)`, `C = dag_to_cpdag(D)`, pick a query `(x, y)` with at
least one possibly-causal path, build `make_linear_iscm(D)` and the population `Σ`, require
`|τ| ≥ 1e-6`. Stream and grids imported unchanged from `run_x1.py` so the ensembles are comparable.

**Eligibility (applied before any outcome is computed).** `C` must be amenable for `(x, y)` and
`O*(C)` must be a valid adjustment set in `D`. SCMs failing either are dropped and counted.
This doubles as the **placebo**: with no deletion the silent-bias rate is 0 by construction.

**Perturbation.** Exhaustive over every edge of `D`, one at a time. No sampling, therefore no
sampling noise in the stratification.

**Two arms.**

| arm | analyst's graph | what it isolates |
|---|---|---|
| **S0** (primary) | `G₀ = Ĉ` | the discovery error alone, no background knowledge. This is the research question as asked. |
| **SK** (secondary) | `G₀ = Meek(Ĉ, K↾)` | the same error inside our actual setting. `K↾` keeps the statements of `K` whose edge still exists and is still undirected in `Ĉ`; the surviving count is recorded. `MeekFail` is counted as *caught free*. |

`K` is built exactly as in `run_x1.py`: the true orientations of up to `MAX_K = 4` undirected edges
of the TRUE CPDAG `C`, in the same permuted order.

**Outcome per (SCM, deleted edge), evaluated against the TRUE DAG `D`:**

| outcome | meaning |
|---|---|
| `meekfail` | (SK only) `K↾` is inconsistent with `Ĉ` — caught free by the graphical check |
| `not_amenable` | `O*(G₀)` undefined — **loud** failure, the analyst sees it |
| `unchanged` | `O*` equals the unperturbed `O*` |
| `changed_valid` | `O*` moved but is still valid in `D` — efficiency only |
| **`silent_bias`** | `O*` moved and is **not** valid in `D` — the only outcome with no alarm |

Bias magnitude is `|ols_coefficient(Σ, x, y, O*) − τ|` on the **population** `Σ`, so Monte-Carlo
error is exactly zero.

**Stratifier.** `hop = stmt_hop(distC, a, b)` for the deleted edge `(a, b)`, with
`distC = hop_dist_from(skeleton(C), [x, y], p)` — minimum over both endpoints, on the TRUE CPDAG
skeleton. Identical convention to `x1_ops.py:271`, so the numbers are comparable to the locality
table.

**Ensembles.** `original` (p 5–8), `licensed` (p 5–10), `large` (p 15–25), grids unchanged.
`large` is the only ensemble carrying real mass at hop ≥ 2.

## 4. Pre-registered decision rule

Primary statistic: `r(h)` = silent-bias rate among perturbations at hop `h`, arm S0, ensemble
`licensed`.

| verdict | condition |
|---|---|
| **SUPPORTED** | `r(0) / r(≥1) ≥ 5` **and** `r(≥1) < 0.02` |
| **INCONCLUSIVE** | anything between |
| **REFUTED** | `r(0) / r(≥1) < 2` **or** `r(≥1) ≥ 0.05` |

A REFUTED verdict is a real outcome and gets written into the paper as a limitation, not buried.
It would mean the framework needs the discovery output to be globally correct, which is the strong
form of the assumption the research question objects to.

Confidence intervals: Wilson, on the perturbation count. ⚠️ **Perturbations from one SCM are not
independent** (same `D`, same query). The rates are therefore reported with a **paired SCM-level
cluster bootstrap** as the primary interval, with the Wilson interval shown beside it, and the design
effect printed. This is the correction X1 already applies and the one the 22/08 assessment was right
to ask for.

## 5. Gates, declared now

| gate | condition | why |
|---|---|---|
| **G1** placebo | silent-bias rate with **no** deletion = 0.0000 | eligibility makes this true by construction; if it is not, the harness is wrong |
| **G2** comparability | `licensed` reproduces `run_x1.py`'s SCM stream (same seeds → same `x, y, K, O*(C)`) | the ensembles must be the same population as the locality table |
| **G3** dynamic range | the overall silent-bias rate is neither 0 nor 1 | an instrument stuck at either end measures nothing |
| **G4** hop coverage | at least 200 perturbations at hop ≥ 1 in `licensed` | the contrast needs mass on both sides |
| **G5** amenability | the `not_amenable` rate is reported, never folded into `silent_bias` | a loud failure is a different event, and pricing it as a breakdown inflates fragility |

## 6. What this experiment cannot answer

- It does not model **edge addition**, nor a discovery algorithm's actual error process.
- It deletes **one** edge. Nothing here speaks to several simultaneous skeleton errors.
- The SCMs are **linear-Gaussian iSCM**. The var-sortability leak is removed by the parameterisation
  and nothing else about nonlinearity is claimed.
- `Σ` is the **population** covariance. This is a structural statement, not a finite-sample one.

---

## 7. AMENDMENT, 2026-08-31, before any outcome was computed

A 60-draw machinery smoke test was run to check the harness. **No outcome was computed from it**;
only the eligibility status distribution was read, which is a property of the design, not a result:

```
drop_C_nonamenable   18
ok                    1
```

**The TRUE CPDAG is not amenable for the query in 18 of 19 returned draws.** With no background
knowledge the optimal adjustment set is not defined at all, so arm S0 as pre-registered would have
run on a ~5 % subset selected precisely on amenability. That is a bad primary arm: the contrast
would be confounded with the selection.

It is also, on its own, a substantive answer to the research question, and it will be reported as
one: **without background knowledge the question does not arise**, because the discovery output
alone identifies nothing here. That is why `K` is in the framework.

### The change

- **Arm SK becomes PRIMARY.** Eligibility becomes exactly `run_x1.py`'s: `G₀ = Meek(C, K)` must exist
  (no `MeekFail`), `O*(G₀)` must be defined, and it must be valid in `D`. This is our actual setting
  and the one the research question is about.
- **Arm S0 becomes SECONDARY**, computed only on the sub-population where `C` is itself amenable, and
  reported with that selection stated.
- `cpdag_amenable` is recorded for every eligible SCM so the rate above is measured on the full run
  rather than quoted from a smoke test.

The decision rule of §4 is unchanged and now reads on **arm SK, ensemble `licensed`**. Gates
unchanged, except **G1's placebo** now reads on the SK baseline, and **G4's 200-perturbation floor**
now reads on arm SK.

Nothing else was inspected before this amendment was written.
