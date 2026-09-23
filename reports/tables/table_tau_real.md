# Real-structure anchor table — τ_b(predictor, survival AUC) by stratum

The real-structure counterpart of `table_tau_comparisons.md`, on the **committed real-network corpus**: the 831 admissible rows on 25 networks in `results/axisa3/instances.jsonl`. Same schema, same footnoting discipline. **The two tables are never pooled and no row here may be compared with a synthetic row** — the corpora differ in structure, in how knowledge is obtained, and in what an instance is.

Three things about this corpus change how the table must be read, and all three were measured before any τ existed (`results/axis_robustness_real/PREREGISTRATION.md` §1.5, §2.1, Appendix A):

1. **The unit of analysis is the network, not the pair.** Every interval here is a cluster bootstrap over the 25 networks, 10,000 resamples, seed 0 — not over rows. The intraclass correlation of the `r = 1` indicator is 0.40, a design effect of about 12.7, so a row-level interval would be roughly three and a half times too narrow. `networks` is printed beside `n` for that reason, and no per-row interval appears anywhere.
2. **`|K|` and `k_g0` are constants within a `(network, coverage)` cell** on this corpus: `select_knowledge` returns one claim set per network and coverage, and `flip` preserves list length. Those two baselines can therefore rank *networks* but not *queries*, which is a weaker thing than they could do synthetically. A win over them here is worth correspondingly less, and that is stated rather than banked.
3. **`SHD(G₀, truth)` needs base wrongness to exist at all.** The corpus draws the analyst's claims from the true orientations, so `G₀` *is* the ground-truth DAG at coverage 1.0 and `SHD ≡ 0`. Base wrongness is applied before the sweep, exactly as the synthetic design did. Appendix A records that the rate-based mechanism reverses nothing at all on most of this corpus because `|K|` is 1 or 2 there, and adds an absolute one-claim level, reported in the supplementary table.

## The nine pre-registered strata

| stratum | n | networks | endpoint used | τ_b(r_hop) [95% CI] | τ_b(shd_truth) [95% CI] | τ_b(|K|) [95% CI] | τ_b(k_g0) [95% CI] |
|---|---|---|---|---|---|---|---|
| flip, coverage=0.5, base_wrongness=0.00 | 182 | 17 | `AUC_frac_usable` | +0.191 [-0.441, +0.646][^loo] | +0.426 [-0.095, +0.777] | +0.179 [-0.400, +0.791][^loo] | +0.269 [-0.382, +0.724][^loo] |
| flip, coverage=0.5, base_wrongness=0.10 | 546 | 17 | `AUC_frac_usable` | +0.140 [-0.442, +0.561][^loo] | +0.432 [-0.060, +0.740] | +0.154 [-0.394, +0.770][^loo] | +0.241 [-0.399, +0.663][^loo] |
| flip, coverage=0.5, base_wrongness=0.25 | 515 | 17 | `AUC_frac_usable` | +0.112 [-0.468, +0.567][^loo] | +0.540 [+0.190, +0.753] | +0.277 [-0.283, +0.839] | +0.309 [-0.343, +0.649] |
| flip, coverage=1.0, base_wrongness=0.00 | 540 | 24 | `AUC_frac_usable` | +0.064 [-0.351, +0.381][^loo] | undefined (predictor constant)[^const] | +0.303 [+0.114, +0.565] | +0.150 [-0.144, +0.375] |
| flip, coverage=1.0, base_wrongness=0.10 | 1579 | 24 | `AUC_frac_usable` | +0.126 [-0.299, +0.414][^loo] | +0.260 [+0.043, +0.522] | +0.321 [+0.121, +0.584] | +0.163 [-0.141, +0.377] |
| flip, coverage=1.0, base_wrongness=0.25 | 1386 | 23 | `AUC_frac_usable` | +0.292 [-0.251, +0.608] | +0.212 [-0.026, +0.512] | +0.223 [-0.032, +0.589] | +0.163 [-0.159, +0.355] |
| tiered, n_tiers=2 | 60 | 11 | `AUC_rate_usable` | +0.580 [+0.120, +0.854] | +0.616 [+0.009, +0.781] | -0.173 [-0.694, +0.515][^loo] | -0.092 [-0.519, +0.523][^loo] |
| tiered, n_tiers=3 | 176 | 18 | `AUC_rate_usable` | +0.378 [+0.061, +0.687] | +0.308 [-0.163, +0.553] | -0.066 [-0.469, +0.349][^loo] | +0.227 [-0.133, +0.515] |
| tiered, n_tiers=4 | 265 | 19 | `AUC_rate_usable` | +0.483 [+0.278, +0.703] | +0.255 [-0.193, +0.491] | +0.124 [-0.343, +0.446][^loo] | +0.271 [-0.097, +0.489] |

## Supplementary strata

Not among the nine, and not comparable with the synthetic table. Coverage 0.25 has no synthetic counterpart; the `bw_abs=1` strata exist to give `SHD(G₀, truth)` a non-degenerate comparison on every network (Appendix A). **No rate is compared across coverage levels**: admissible instances collapse 543 → 182 → 106 as coverage falls and the survivors have larger separation, so such a comparison measures composition ([RE-4], [RE-6]).

| stratum | n | networks | endpoint used | τ_b(r_hop) [95% CI] | τ_b(shd_truth) [95% CI] | τ_b(|K|) [95% CI] | τ_b(k_g0) [95% CI] |
|---|---|---|---|---|---|---|---|
| flip, coverage=0.25, base_wrongness=0.00 | 106 | 8 | `AUC_frac_usable` | -0.527 [-0.922, -0.090] | -0.322 [-0.792, +0.670][^loo] | +0.824 [+0.411, +0.990] | -0.618 [-0.899, +0.091] |
| flip, coverage=0.25, base_wrongness=0.10 | 318 | 8 | `AUC_frac_usable` | -0.513 [-0.918, -0.097] | -0.337 [-0.760, +0.629][^loo] | +0.793 [+0.365, +0.987] | -0.565 [-0.866, +0.116] |
| flip, coverage=0.25, base_wrongness=0.25 | 318 | 8 | `AUC_frac_usable` | -0.521 [-0.917, -0.097] | -0.386 [-0.775, +0.653][^loo] | +0.802 [+0.377, +0.987] | -0.524 [-0.869, +0.194] |
| flip, coverage=0.25, one claim reversed | 112 | 5 | `AUC_frac_usable` | +0.079 [-0.213, +0.721][^loo] | -0.264 [-0.377, +0.288][^loo] | undefined (predictor constant)[^const] | +0.005 [-0.084, +0.056][^loo] |
| flip, coverage=0.5, one claim reversed | 310 | 15 | `AUC_frac_usable` | +0.444 [+0.107, +0.775] | -0.302 [-0.621, +0.202] | -0.042 [-0.505, +0.590][^loo] | +0.311 [-0.117, +0.786] |
| flip, coverage=1.0, one claim reversed | 1146 | 21 | `AUC_frac_usable` | +0.300 [+0.045, +0.543] | +0.020 [-0.258, +0.246][^loo] | +0.179 [-0.035, +0.414] | +0.072 [-0.195, +0.356][^loo] |

### Footnotes

[^const]: **Undefined by construction, not zero.** The predictor takes a single value across the whole stratum, so no rank correlation exists. Printed as `undefined (predictor constant)` rather than `τ = 0`, and never counted as evidence that the baseline "fails".
[^loo]: **⚠ single-network leverage — do not quote this cell on its own.** The cell's verdict does not survive removing one network: either some leave-one-network-out τ has the opposite sign to the full-sample τ, or the full-sample interval excludes zero while the leave-one-out range spans it. Within-network variation in `r_hop` is thin on this corpus — at coverage 1.0, 10 of the 25 networks have a single distinct radius across all their rows, while `paths` alone spans `r = 1…14` in 20 rows — so a single network can carry a stratum. `analysis_tau.csv` names the network at the minimum and the maximum.

**Assumption (every `r_hop` row):** Conjecture 2 (hence Anti-Exchange Case B, verified not proved); the error is one-sided, so radii can only be too large, never too small (`THEOREMS.md` §4c, §6, §8).

**Gate:** `benchmarks.measure.fast_gate` exclusively, on every row.

**Draws:** N = 1000 per grid point. **Endpoint:** the conservative `_usable` variant is tabulated; the raw variant is in `analysis_tau.csv` and every disagreement is marked.

**Provenance:** `results/axis_robustness_real/analysis_tau.csv`, built by `bkrobust.robustness.real_analyse` from the sharded sweep, which was built on the frozen frame `results/axis_robustness_real/frame.jsonl` (sha256 `73881798319206be…`). Regenerate with `PYTHONPATH=src python -m bkrobust.robustness.build_table_tau_real`.
