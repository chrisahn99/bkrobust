# X3 — skeleton perturbation, and a concrete use of the framework

Built 2026-08-31 for the Thursday meeting on the robustness of background knowledge. Two
deliverables, one for each half of what the Teams message asked for.

| what the message asked | what is here |
|---|---|
| *« préparer un cas concret d'utilisation du framework »* | `results/worked_example.txt` (+ `code/worked_example.py`) |
| *« l'assumption que la data driven method est correcte est très forte »* | `PREREG.md` → `RESULTS.md` (+ `code/run_x3.py`, `code/analyse_x3.py`) |

## 1. The experiment: does the locality law survive a wrong SKELETON?

Every measurement in this line uses the oracle CPDAG, so the research question had no answer.
X3 gives it one, on **9 721 SCMs and 130 120 single-edge deletions**, population Σ, zero Monte-Carlo
error.

🔴 **The pre-registered rule returns INCONCLUSIVE**, and that is the headline. It required
`r(hop 0)/r(hop ≥ 1) ≥ 5` **and** `r(hop ≥ 1) < 0.02`. Measured: ratio **7.10** ✅, but
`r(hop ≥ 1) = 0.0231` ❌. A missed edge one hop from the query still breaks the adjustment set
2.3 % of the time, above the bar set in advance.

⚠️ **The sharp cut is at hop ≥ 2, and it is POST HOC.** Pooled: **0.1376** at hop 0 against
**0.000423** at hop ≥ 2 (15 / 35 501), a **326×** separation. This must be re-registered and re-run
before it goes in a paper. Today's honest sentence is: *the pre-registered test is inconclusive, and
the data suggest the threshold is two hops rather than one.*

Four things the run says regardless of the verdict:

- **Without background knowledge the question does not arise.** The CPDAG alone is amenable for the
  query in **8.2 % / 16.9 % / 55.8 %** of problems. On 83 % of `licensed` there is nothing to
  proceed with, which is *why* `K` is in the framework.
- 🔑 **The efficiency-only band is NOT empty here.** `changed_valid` = **2.35 / 2.65 / 5.53 %**,
  against **0 / 754** for our orientation operator. Chris's `Severity.EFFICIENCY` band exists; our
  instrument was the wrong one. That corrects **our** side of the 2026-08-20 analysis.
- ❌ **"The edges a CI test misses are weak, so they don't matter" is not supported.** Median |w| of
  the deleted edge is indistinguishable between damaging and harmless deletions, and inverted on
  `large`. Distance predicts damage; strength does not.
- **Locality cuts the frequency of damage, not its severity.** Bias is heavy-tailed at every hop
  (licensed hop ≥ 1: median 0.39, q90 1.75, max 1204).

Design effects run 0.73–1.81 and the cluster-bootstrap intervals sit on top of the Wilson ones, so
the clustering objection of 2026-08-22 is answered by construction.

## 2. The worked example: one problem, end to end

`results/worked_example.txt`. Synthetic, from the same generator as every other run; the variable
names are attached **by topological position** so the skin is coherent, and they carry no empirical
claim. Query: effect of `Statin` on `CardiacEvent`.

It walks the chain the framing asks for — CPDAG → `K` → `G₀` → `O*` → `r_val` + witness → `r_ε` —
and lands on the case that justifies the whole design:

```
r_val = 2      the adjustment set stops being provably valid here
r_eps > 2      the estimate does not move at all until r = 3
```

**Fragile identification, stable estimate** — the middle row of the framing's own §5 table, and the
row it says a validity-only analysis raises a false alarm in. Collapsing the two radii into one would
have called this analysis fragile at `r = 2`, while the number it is alarmed about had not moved by
1e-9. At `r = 3` it moves by **0.4520 against an effect of −0.4191**, i.e. 108 %.

And the witness, which is the output that matters:

> *"The estimate holds unless `BMI → BloodPressure` is retracted and `Exercise → BMI` is retracted."*

A domain expert can argue with that sentence without any metric literacy.

## 3. 🆕 A structural remark that fell out of building it

Across every case searched, the last move of the minimal breaking path was **always a retraction,
never an assertion**. That is forced, not luck:

> if `G'` refines `G` then `[G'] ⊆ [G]`, so a `Z` valid in every DAG of `[G]` is valid in every DAG
> of `[G']`.

**Orienting one more edge can never invalidate a fixed valid `Z`.** For `r_val`, where `Z` is fixed
at `O*(G₀)`, only the **upward** moves bind. The downward move earns its place through `r_ε`, where
`O*` is re-derived and the estimate does move. This removes an entire direction from the `r_val`
search, which matters for §8's cost statement.

## Files

```
PREREG.md                  written first, amended once (§7) before any outcome
RESULTS.md                 the full run, gates, and what it cannot answer
code/run_x3.py             the experiment; stream imported unchanged from run_x1.py
code/analyse_x3.py         exactly the pre-registered statistics
code/worked_example.py     the end-to-end case (--search to re-pick a seed)
results/x3_{ens}.json      raw run registers
results/x3_analysis.json   every number in RESULTS.md
results/worked_example.txt the rendered case
```

## Related

- [[../../wiki/activities/reunion-bk-robustesse-2026-09-03|The meeting this was built for]]
- [[../../wiki/projects/rho-breakdown-adjustment|ρ-breakdown]]
- [[../2026-08-19_e1prime-se-criterion/RESULTS.md|E1′ — the orientation locality law this tests against]]
- [[../2026-08-31_chris-rho-breakdown-package/README|The package for Chris]]
