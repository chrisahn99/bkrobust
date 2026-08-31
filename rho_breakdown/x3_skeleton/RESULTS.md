# X3 — skeleton perturbation: results

**Run 2026-08-31.** Pre-registered in `PREREG.md` before the first number, amended once (§7) after a
machinery smoke test and before any outcome was computed.
**9 721 analysed SCMs, 130 120 single-edge deletions**, population `Σ`, zero Monte-Carlo error.
Reproduce: `code/run_x3.py` then `code/analyse_x3.py`.

---

## The verdict the rule returns: INCONCLUSIVE

| pre-registered statistic (arm SK, `licensed`) | value | bar | |
|---|---|---|---|
| `r(hop 0) / r(hop >= 1)` | **7.10** | >= 5 | ✅ |
| `r(hop >= 1)` | **0.0231** | < 0.02 | ❌ |

**SUPPORTED needs both. It gets one.** The rule therefore returns **INCONCLUSIVE**, and that is the
result. H1 as written is not established: a missed edge one hop from the query still breaks the
adjustment set in **2.3 %** of cases, which is above the bar I set in advance.

It is not REFUTED either (that needed ratio < 2 or `r(hop >= 1)` >= 0.05), so the locality intuition
survives to be re-registered on a sharper cut. Which brings the next section, and its label.

---

## 🔴 The sharp cut is at hop >= 2, and it is POST HOC

I did not pre-register the hop-1 / hop-2 split; I pre-registered hop 0 against hop >= 1. Lumping
hop 1 in with hop >= 2 is what produced the INCONCLUSIVE, because **hop 1 is where the residual
damage lives and hop >= 2 is where it stops**. Pooled over all three ensembles, arm SK:

| distance from the deleted edge to {x, y} | silent bias | rate |
|---|---|---|
| hop 0 | 5 066 / 36 829 | **0.1376** |
| hop >= 2 | **15 / 35 501** | **0.000423** |

**A 326-fold separation**, on 35 501 perturbations at hop >= 2. Per ensemble:

| ensemble | hop 0 | hop 1 | hop >= 2 |
|---|---|---|---|
| `original` (p 5–8) | 0.1599 | 0.0458 | **0.0000** (0 / 1 879) |
| `licensed` (p 5–10) | 0.1638 | 0.0269 | **0.000286** (1 / 3 493) |
| `large` (p 15–25) | 0.0718 | 0.0176 | **0.000465** (14 / 30 129) |

⚠️ **This cut was chosen after seeing the data. It is a hypothesis, not a result, and it must be
re-registered and re-run before it is written into a paper.** The honest statement today is: *the
pre-registered test is inconclusive, and the data suggest the threshold sits at two hops rather than
one.* Writing it the other way round is the failure mode this vault has a register for.

---

## What the run says regardless of the verdict

### 1. Without background knowledge, the question does not arise

`O*` is undefined on the raw discovery output for most problems, because the CPDAG is not amenable
for the query:

| ensemble | CPDAG amenable without `K` |
|---|---|
| `original` | **224 / 2 748 = 8.2 %** |
| `licensed` | **652 / 3 853 = 16.9 %** |
| `large` | **1 740 / 3 120 = 55.8 %** |

This is why `K` is in the framework at all, and it reframes the research question: the objection
*"the assumption that the data-driven method is correct is very strong"* presumes the analyst could
proceed without knowledge. On 83 % of `licensed` problems there is nothing to proceed with.

### 2. 🔑 The efficiency-only band is NOT empty here, and that repairs a collision

The band where `O*` moves but stays valid — `Severity.EFFICIENCY` in Chris's `THEORY.md`, *"the one
worth naming"*, on which T3 and all of `efficiency_gap.py` rest — was measured **empty** by our
orientation operator: **0 / 754**. Under skeleton perturbation it is **not** empty:

| ensemble | `changed_valid` |
|---|---|
| `original` | 463 / 19 740 = **2.35 %** |
| `licensed` | 1 258 / 47 395 = **2.65 %** |
| `large` | 3 485 / 62 985 = **5.53 %** |

**The band exists; our operator was the wrong instrument for it.** That is a correction to our side
of the 2026-08-20 analysis, not to his, and it should be handed over as such.

### 3. The full outcome taxonomy, arm SK

| outcome | `original` | `licensed` | `large` |
|---|---|---|---|
| `unchanged` | 0.7212 | 0.7883 | 0.8432 |
| `not_amenable` (**loud**, the analyst sees it) | 0.1728 | 0.1106 | 0.0848 |
| `silent_bias` | 0.0825 | 0.0746 | 0.0166 |
| `changed_valid` (efficiency only) | 0.0235 | 0.0265 | 0.0553 |

The `not_amenable` column is the fourth exit that `cl(G, a -> b)` does not currently have: it returns
bottom for a cycle or a forbidden structure, and nothing for *`O*` undefined*. At 8–17 % it is not a
corner case, and pricing it as a breakdown would inflate fragility.

### 4. ❌ "The edges a CI test misses are weak, so they do not matter" is not supported

Median absolute weight of the deleted edge, silent-bias events against everything else:

| ensemble | silent | not silent |
|---|---|---|
| `original` | 0.6825 | 0.6600 |
| `licensed` | 0.5726 | 0.5583 |
| `large` | 0.5675 | 0.6028 |

Indistinguishable, and inverted in `large`. **Edge strength does not predict damage. Distance
does.** The comfortable version of the finite-sample story does not survive contact with the
measurement.

### 5. Locality cuts the FREQUENCY of damage, not its SEVERITY

Relative absolute bias among silent-bias events, arm SK:

| ensemble | hop 0 median / q90 / max | hop >= 1 median / q90 / max |
|---|---|---|
| `original` | 0.61 / 2.86 / 228 | 0.58 / 2.54 / 49 |
| `licensed` | 0.47 / 2.13 / 1 496 | 0.39 / 1.75 / 1 204 |
| `large` | 0.48 / 2.15 / 92 | 0.43 / 1.65 / 77 |

Heavy-tailed at every distance, and barely lighter far away. So the claim licensed by this run is
*a distant skeleton error rarely does damage*, *never* *a distant skeleton error does little damage*.
Same shape as the pilot's orientation finding, and the same reason a continuity bound is unavailable.

### 6. Clustering is real and mild, and it is corrected

Design effects run **0.73 to 1.81** across every cell. The paired SCM-level cluster-bootstrap
intervals sit essentially on top of the Wilson intervals (`licensed` hop 0: cluster
[0.1586, 0.1692] against Wilson [0.1585, 0.1693]). The clustering objection raised on 2026-08-22 is
answered here by construction rather than by argument.

---

## Gates

| gate | result |
|---|---|
| **G1** placebo | 0.0000 by construction: eligibility requires the unperturbed `O*` to be valid in `D` |
| **G2** comparability | SCM stream imported unchanged from `run_x1.py` through `K`; grids identical |
| **G3** dynamic range | overall silent rate 0.0839 on `licensed`, neither 0 nor 1 ✅ |
| **G4** hop mass | 23 940 perturbations at hop >= 1 on `licensed`, floor was 200 ✅ |
| **G5** amenability | reported separately at 8–17 %, never folded into `silent_bias` ✅ |

## What this cannot answer

- **Edge addition is not tested.** Declared in the prereg: the failure direction at `n <= d` is
  omission. A discovery method that hallucinates an edge is outside this run.
- **One edge at a time.** Nothing here speaks to several simultaneous skeleton errors, and they are
  the realistic case.
- **The error model is an assumption.** `Chat = dag_to_cpdag(D minus one edge)` is not what a CI-test
  failure actually returns. It buys an always-valid CPDAG and an interpretable perturbation, and it
  is stated rather than derived.
- **Linear-Gaussian iSCM, population `Σ`.** This is a structural statement, not a finite-sample one,
  and the regime the b-LOAD reviewers called too narrow is still the regime.

## Related

- `PREREG.md` — written first, amended once before any outcome
- [[../2026-08-19_e1prime-se-criterion/RESULTS.md]] — the orientation locality law this tests against
- [[../../wiki/activities/reunion-bk-robustesse-2026-09-03]] — the research question this answers
- [[../../wiki/projects/rho-breakdown-adjustment]]
