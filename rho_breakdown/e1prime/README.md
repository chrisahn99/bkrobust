---
title: "E1′ — the `se` criterion: ρ* stops being an oracle statistic"
type: synthesis
created: 2026-08-19
updated: 2026-08-19
tags: [research, causal-inference, iclr-2027, rho-breakdown, prereg, betelgeuse]
sources:
  - output/2026-08-14_latent-causal-iclr2027-rerank/RESULTS.md
  - output/2026-08-14_latent-causal-iclr2027-rerank/RANKING.md
  - output/2026-08-14_latent-causal-iclr2027-rerank/CHRIS-MEETING.md
  - output/2026-08-14_apostila-rho-breakdown/07-experimentos.tex
  - wiki/projects/rho-breakdown-adjustment.md
---

# E1′ — the `se` criterion

**Run 2026-08-19 on betelgeuse** (16 cores, no GPU, no network), `~/latent-causal/e1prime-se`.
Pre-registration: [`PREREG.md`](PREREG.md) + Addendum 1. Design audit: [`AUDIT.md`](AUDIT.md).

## What this fixes

E1 (14/08) answered the reviewer's *"ρ\* is 4-valued and censored at the top, 90 % of the mass in
one bin"* by raising |K| from 4 to 8 and switching from `sign` to `rel50`. SUPPORTED on 3 of 4
ensembles; **3 of 3 adversarial verifiers refuted the interpretation.** Six defects were on record.
This run attacks the one nobody had touched — and it is the deepest:

> 🔴 **Both of E1's criteria compare the perturbed estimate to `tau`, the true total effect.**
> An analyst does not have `tau`. The paper's product is therefore not computable by the person it
> is for. `RANKING.md:325-331` says the *only* thing separating ρ\* from gadjid / SID / the
> Taeb–Guo–Henckel radius is that those are computed around a **known** G₀ — and then concedes the
> pilot never exercises the analyst-facing quantity. `RANKING.md:839` lists `crit='se'` under
> **"Experiments not run"**, calls it *"the single highest-value remaining experiment"*, and prices
> it at ~1 day.

The replacement yardstick is the analyst's **own standard error**:

```
rho*_se(n) = min { rho : max_m |est(K') - est0| > z * SE_n(est0) }
```

*How few expert mistakes would move my estimate by more than my own sampling error.* No `tau`.

## What is new relative to E1, beyond the criterion

| | E1 (14/08) | E1′ |
|---|---|---|
| yardstick | `tau` (oracle) | `z·SE_n` (analyst's own) |
| |K| | **8** — definable on 1.5 % of the licensed box | **4**, the pilot's own cap, ensemble unconditioned |
| ball | ρ ≤ 3, **budget-censored** | **full**, ρ = 1…\|K\| — censoring is measured, not imposed |
| the bar | censored ≤ 0.60 — *cleared for free* by `ident` | a **location** (`n50`, `nND` ≤ 20 000), which nothing clears for free |
| censored bin | one bin | **decomposed**: caught-free / loud / genuinely robust (Addendum 1) |
| expert | **true by construction** | ARM 2: genuinely wrong, `r ∈ {0,1,2}` reversals, finite samples |
| negative control | none | placebo P1 (frozen `O0` → must censor 100 %) + floor P2 |

## Arms

- **ARM 1** — population Σ, exact, zero Monte Carlo. 4 ensembles (`licensed`, `original`, `large`,
  `k8`), n-grid = the registered falsifier's own {200 … 2e4, ∞} plus diagnostics, z ∈ {1, 1.96, 2.58}.
- **ARM 2** — finite-sample analyst simulation. The analyst gets only `(C, K_asserted, Σ̂, n)`.
  Coverage of `tau` and interval **width** for `naive` vs `robust(ρ)`. `tau` enters once, to
  adjudicate — never as input.

## Reproduce

```bash
ssh betelgeuse
cd ~/latent-causal/e1prime-se/code
PY=$HOME/e1venv/bin/python NW=12 ./run_all.sh
```

## Results

_(filled by the run — `RESULTS.md`)_

## Related Pages

- [[wiki/projects/rho-breakdown-adjustment|ρ-Breakdown Radius — the project page]]
- [[output/2026-08-14_latent-causal-iclr2027-rerank/RESULTS|E1 (14/08), the run this fixes]]
- [[output/2026-08-14_latent-causal-iclr2027-rerank/RANKING|The ranking + the un-run experiment list]]
- [[output/2026-08-14_apostila-rho-breakdown/README|A apostila ρ-breakdown (cap. 7 = E1)]]
- [[wiki/entities/chris|Chris (Seong Woo Ahn)]]
