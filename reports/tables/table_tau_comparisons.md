# Publication anchor table -- τ_b(predictor, survival AUC) by stratum

This table reports τ_b between three predictors (`r_val`, `shd_truth`, `n_k`) and the conservative-control survival AUC (grid points with `n_eval ≥ 30`) across the nine pre-registered strata, using `AUC_frac_usable` for the flip arm (targeted-depth `d_claims` axis) and the newly computed `AUC_rate_usable` for the tiered arm (`corruption_rate` axis, per Appendix E) -- the two arms' endpoints are on different axes and are never compared row-to-row. `r_val`'s CI excludes zero in 8 of 9 strata and is the sole predictor to do so (bolded) only in flip cov=0.5 bw=0.10; elsewhere it shares significance with, or is beaten by, the SHD/`n_k` baselines, or (flip cov=1.0 bw=0.00) all three are indistinguishable from zero. `n_k`'s association is negative rather than absent: in every tiered stratum (n_tiers 2/3/4, τ ≈ −0.42/−0.38/−0.38) and in the higher-precision null-cell re-run below, its 95% CI excludes zero on the negative side, so H4's "inert baseline" framing is wrong. One stratum, flip cov=0.5 bw=0.25, carries a confirmed null for `r_val` at high precision: the N=1000 re-run (Appendix H) gives τ = +0.032, essentially unchanged from τ = +0.031 at N=200 unfiltered (five times the draws moved it +0.001), even though that same stratum's N=200 `AUC_frac_usable` value used in the main table row above is a markedly higher +0.279 -- a usable-depth-filter artifact at low draw count, not a draw-count attenuation effect (see the source-data note in the accompanying report). `shd_truth` is undefined by construction (predictor constant) in exactly one stratum, flip cov=1.0 bw=0.00, where full coverage and zero base wrongness make `G₀` the ground truth for every instance.

| stratum | n | endpoint used | τ_b(r_val) [95% CI] | τ_b(shd_truth) [95% CI] | τ_b(n_k) [95% CI] |
|---|---|---|---|---|---|
| flip, coverage=0.5, base_wrongness=0.00 | 239 | `AUC_frac_usable` | 0.634 [0.533, 0.728] | 0.204 [0.105, 0.300] | 0.184 [0.076, 0.289] |
| flip, coverage=0.5, base_wrongness=0.10 | 240 | `AUC_frac_usable` | **0.334 [0.206, 0.453]** | -0.000 [-0.111, 0.107] | -0.004 [-0.114, 0.105] |
| flip, coverage=0.5, base_wrongness=0.25 | 220 | `AUC_frac_usable` | 0.279 [0.153, 0.397] ⚠[^unstable] | 0.032 [-0.063, 0.127] | 0.131 [0.000, 0.256] |
| flip, coverage=1.0, base_wrongness=0.00 | 229 | `AUC_frac_usable` | -0.091 [-0.204, 0.030] ⚠[^unstable] | undefined (predictor constant)[^shdconst] | 0.021 [-0.062, 0.103] |
| flip, coverage=1.0, base_wrongness=0.10 | 226 | `AUC_frac_usable` | 0.474 [0.393, 0.548] | 0.263 [0.160, 0.367] | 0.261 [0.171, 0.349] |
| flip, coverage=1.0, base_wrongness=0.25 | 133 | `AUC_frac_usable` | 0.355 [0.223, 0.476] | 0.582 [0.488, 0.666] | 0.557 [0.473, 0.635] |
| tiered, n_tiers=2 | 190 | `AUC_rate_usable` | 0.791 [0.763, 0.814] | 0.361 [0.287, 0.434] | -0.423 [-0.516, -0.327] |
| tiered, n_tiers=3 | 238 | `AUC_rate_usable` | 0.756 [0.728, 0.781] | 0.169 [0.079, 0.257] | -0.382 [-0.459, -0.303] |
| tiered, n_tiers=4 | 240 | `AUC_rate_usable` | 0.722 [0.690, 0.750] | 0.195 [0.103, 0.286] | -0.377 [-0.453, -0.300] |

**Supplementary -- null-cell high-precision re-measurement**

| stratum | n | endpoint used | τ_b(r_val) [95% CI] | τ_b(shd_truth) [95% CI] | τ_b(n_k) [95% CI] |
|---|---|---|---|---|---|
| flip, coverage=0.5, base_wrongness=0.25 -- null-cell re-run, N=1000 | 220 | `AUC_frac_usable` | 0.032 [-0.112, 0.168] | -0.170 [-0.260, -0.080] | -0.200 [-0.334, -0.068] |

[^unstable]: **Unresolved at N = 200 — do not quote this cell on its own.** The
`n_eval >= 30` filter is not a conservative control: it conditions on `n_eval`, which is
itself correlated with `r_val` (tau from +0.318 to -0.438 across strata, Appendix F.2), so
it selects different grid points for high- and low-`r_val` instances. In these two strata it
changes the verdict, in opposite directions. On the unfiltered `AUC_frac`, coverage=0.5
bw=0.25 is +0.031 (CI includes 0) against +0.279 here, and coverage=1.0 bw=0.00 is +0.519
(CI excludes 0) against -0.091 here -- a sign flip of -0.610. The diagnostic is decisive:
at N = 1000 draws on the same 220 instances, the two endpoints agree to 0.0002 (Appendix H),
so divergence signals too few draws rather than a safer estimate. For coverage=0.5 bw=0.25
the supplementary N = 1000 row below supersedes this cell and the stratum is a **confirmed
null**. For coverage=1.0 bw=0.00 no N = 1000 run exists and the stratum is **unresolved**.
See Appendix J.

[^shdconst]: `shd_truth ≡ 0` by construction there (full coverage, zero base wrongness => `G₀` *is* the ground-truth DAG for every instance in this stratum).

**Assumption (every `r_val` row):** Conjecture 2 (hence Anti-Exchange Case B, verified not proved); the error is one-sided, so radii can only be too large.

**Gate:** benchmarks.measure.fast_gate exclusively.

## Notes

### Tiered `AUC_rate_usable` filter statistics (the gap this table closes)

Rates surviving the `n_eval_rate_* >= 30` filter, per instance, out of the 11-point `corruption_rate` grid (668 tiered instances total):

- min: 3
- median: 8.0
- max: 11
- instances excluded for zero usable rates: 0 (of 668)

### Source-data note: the N=200 vs. N=1000 `AUC_frac_usable` discrepancy at flip cov=0.5 bw=0.25

The main table row above uses session 7's original N=200-draw survival data, filtered to `n_eval >= 30`, per the task's instruction to use the conservative usable control everywhere. For this specific stratum that gives τ_b(`r_val`, `AUC_frac_usable`) = +0.279 (CI excludes zero). The Appendix H null-cell re-run, same 220 instances, N=1000 draws, same `n_eval >= 30` filter, gives τ_b = +0.032 (CI includes zero) -- essentially the *unfiltered* N=200 value (+0.031, PREREGISTRATION.md Appendix H.1) rather than the filtered N=200 value. This is not a contradiction in the source files -- both numbers are reproduced exactly from their respective committed CSVs (see L1 in the accompanying report) -- but it is worth flagging: at N=200 draws the `n_eval >= 30` filter is aggressive in this thin, high-contradiction-rate stratum, and which depths it excludes correlates with `r_val` strongly enough to move τ from +0.03 to +0.28. At N=1000 draws the same threshold is far less exclusionary (raw counts are 5x higher), the selection effect attenuates, and τ reverts to the unfiltered value. The N=1000 result is the more trustworthy one for this stratum precisely because it is less exposed to that selection effect -- which is the whole point of the higher-precision re-run -- but the main table still reports the N=200 usable-filtered figure per the task's fixed rule of using `AUC_*_usable` for every main-table cell, with this note attached rather than silently overriding that rule.
