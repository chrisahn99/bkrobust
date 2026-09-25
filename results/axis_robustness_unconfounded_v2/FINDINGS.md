# Tiered survival with queries drawn uniformly (v2)

Generator: `src/bkrobust/robustness/unconfounded_generator.py` (rule in its docstring).
Driver: `src/bkrobust/robustness/run_unconfounded_v2.py`; analysis:
`experiments/unconfounded_analysis_v2.py`. Same tiered protocol as
`results/axis_robustness_p6/` (component sizes 6-12, 2-4 tiers, 1,000 draws per
grid point, fast_gate, AUC_frac), except that the treatment is drawn uniformly
among admissible vertices instead of at a target separation. 480 instances
(40 per cell), wall clock ~11 min. Per-sample shards and
`survival_samples.csv.gz` are regenerable and not committed.

r_val = 1 on 316/480 (66%). Separation measured on all 480; r_val != separation on 62%.

Marginal Kendall tau_b with AUC_frac (instance bootstrap 95% CI):

| stratum | r_val | separation | n_k | k_g0 | shd_truth |
|---|---|---|---|---|---|
| 2 tiers | 0.45 [0.36, 0.52] | 0.38 | -0.09 | 0.00 | 0.10 |
| 3 tiers | 0.35 [0.22, 0.47] | 0.48 | 0.11 | 0.16 | 0.19 |
| 4 tiers | 0.20 [0.07, 0.32] | 0.42 | 0.20 | 0.25 | 0.28 |
| pooled | 0.29 [0.22, 0.36] | 0.34 [0.28, 0.40] | 0.17 | 0.22 | -0.01 |

Paired r_val minus separation: +0.06 (CI contains 0), -0.13 [-0.24, -0.02],
-0.23 [-0.33, -0.13], pooled -0.05 [-0.13, 0.02]. Pooled, r_val leads SHD
(+0.31 [0.20, 0.41]), n_k (+0.12 [0.03, 0.21]) and r_claim (+0.10 [0.06, 0.13])
and ties k_g0 (+0.07, CI contains 0). Full tables: `marginal_tau.csv`,
`paired_delta.csv`.

Verdict: without separation targeting, the radius does not beat separation as a
predictor of average-case survival; separation wins at 3 and 4 tiers.
