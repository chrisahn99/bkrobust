# Phase 1 Analysis Notes

Endpoint per PREREGISTRATION.md Appendix D.3 (restoring section 4): **primary = `AUC_frac` over `S`, contradictions excluded**, for the flip arm. Per **Appendix E**, the tiered arm's `AUC_frac` (gridded on the defective `d_claims` metric, which undercounts corruption depth because it misses assertions *removed* by node relocation) is re-binned on `corruption_rate` -- the tiered arm's native, exactly-recorded experimental knob -- giving **`AUC_frac_rate`**, the tiered arm's primary endpoint for every table below. `AUC_frac_usable` (restricted to depths with `n_eval >= 30`, still on the d_claims axis for both arms) is the required sensitivity analysis, not a redefinition; for the tiered arm it inherits the d_claims defect and is reported only for continuity.

## Appendix E -- tiered re-binning verification (E1)

Pooled tiered-arm samples (/Users/ahn/Documents/Research/iclr27/bkrobust/results/axis_robustness/survival_samples.csv, tiered rows only, 1469600 rows across 11 rates) checked against PREREGISTRATION.md Appendix E.2's target table:
 rate      n  n_matches_target  contradiction_rate  target_contradiction_rate  contradiction_abs_dev  S_given_consistent  target_S_given_consistent  S_abs_dev
 0.00 133600              True            0.000000                      0.000               0.000000            1.000000                      1.000   0.000000
 0.05 133600              True            0.228862                      0.229               0.000138            0.994477                      0.994   0.000477
 0.10 133600              True            0.404326                      0.404               0.000326            0.962969                      0.963   0.000031
 0.15 133600              True            0.574394                      0.574               0.000394            0.937004                      0.937   0.000004
 0.20 133600              True            0.647582                      0.648               0.000418            0.897946                      0.898   0.000054
 0.25 133600              True            0.691699                      0.692               0.000301            0.878681                      0.879   0.000319
 0.30 133600              True            0.742216                      0.742               0.000216            0.847706                      0.848   0.000294
 0.35 133600              True            0.800584                      0.801               0.000416            0.752571                      0.753   0.000429
 0.40 133600              True            0.813121                      0.813               0.000121            0.739616                      0.740   0.000384
 0.45 133600              True            0.832590                      0.833               0.000410            0.670974                      0.671   0.000026
 0.50 133600              True            0.851729                      0.852               0.000271            0.617800                      0.618   0.000200
- n matches target (133,600) at every rate: True
- S(rate=0.0) == 1.000 exactly: True
- **max_abs_deviation from the target table: 0.000477** (target table itself is printed to 3 decimals, so deviations at this scale are rounding, not disagreement) -- **MATCHES** (within tol=0.002).
- Per-instance re-binning written to `analysis_tiered_rebinned.csv` (668 tiered instances, `AUC_frac_rate` undefined for 0 of them).
- Flip arm is unaffected (targeted depth, exact at generation, no d=0 rows) and never touches survival_samples.csv, per Appendix E.2 scope.

## Verdicts (primary endpoint: flip=AUC_frac, tiered=AUC_frac_rate)

- **P1** (SUPPORTED), endpoints used: AUC_frac / AUC_frac_rate: tau_b(endpoint, r_val) > 0 in 9/9 defined primary strata (of 9 total).
   arm  base_wrongness  coverage  n_tiers      endpoint    tau_b  ci_lo_2p5  ci_hi_97p5   n  underpowered
  flip            0.00       0.5      NaN      AUC_frac 0.713720   0.625286    0.797966 240         False
  flip            0.00       1.0      NaN      AUC_frac 0.519242   0.441501    0.590725 240         False
  flip            0.10       0.5      NaN      AUC_frac 0.379724   0.264552    0.488149 240         False
  flip            0.10       1.0      NaN      AUC_frac 0.387687   0.295891    0.477832 236         False
  flip            0.25       0.5      NaN      AUC_frac 0.031381  -0.114416    0.168647 220         False
  flip            0.25       1.0      NaN      AUC_frac 0.346081   0.209074    0.472125 134         False
tiered             NaN       NaN      2.0 AUC_frac_rate 0.799283   0.761308    0.833541 190         False
tiered             NaN       NaN      3.0 AUC_frac_rate 0.810814   0.786762    0.832792 238         False
tiered             NaN       NaN      4.0 AUC_frac_rate 0.777266   0.752649    0.798170 240         False

- **P2** (SUPPORTED), endpoints used: AUC_frac / AUC_frac_rate: r_val beats n_k in 9/9 strata where both defined.
   arm  base_wrongness  coverage  n_tiers endpoint_rval  tau_b_rval  tau_b_nk
  flip            0.00       0.5      NaN      AUC_frac    0.713720  0.233918
  flip            0.00       1.0      NaN      AUC_frac    0.519242 -0.043271
  flip            0.10       0.5      NaN      AUC_frac    0.379724 -0.040990
  flip            0.10       1.0      NaN      AUC_frac    0.387687 -0.082396
  flip            0.25       0.5      NaN      AUC_frac    0.031381 -0.204544
  flip            0.25       1.0      NaN      AUC_frac    0.346081 -0.193695
tiered             NaN       NaN      2.0 AUC_frac_rate    0.799283 -0.433353
tiered             NaN       NaN      3.0 AUC_frac_rate    0.810814 -0.469739
tiered             NaN       NaN      4.0 AUC_frac_rate    0.777266 -0.425099

- **P3 / H4** (NOT SUPPORTED), endpoints used: AUC_frac / AUC_frac_rate: n_k's CI covers 0 in 3/9 defined strata.
   arm  base_wrongness  coverage  n_tiers      endpoint     tau_b  ci_lo_2p5  ci_hi_97p5   n
  flip            0.00       0.5      NaN      AUC_frac  0.233918   0.117972    0.354242 240
  flip            0.00       1.0      NaN      AUC_frac -0.043271  -0.146295    0.058993 240
  flip            0.10       0.5      NaN      AUC_frac -0.040990  -0.142665    0.061697 240
  flip            0.10       1.0      NaN      AUC_frac -0.082396  -0.188848    0.025398 236
  flip            0.25       0.5      NaN      AUC_frac -0.204544  -0.338929   -0.071216 220
  flip            0.25       1.0      NaN      AUC_frac -0.193695  -0.332259   -0.053190 134
tiered             NaN       NaN      2.0 AUC_frac_rate -0.433353  -0.515114   -0.347388 190
tiered             NaN       NaN      3.0 AUC_frac_rate -0.469739  -0.541185   -0.395815 238
tiered             NaN       NaN      4.0 AUC_frac_rate -0.425099  -0.502142   -0.346591 240

- **P4** (NOT COMPARABLE AS PRE-REGISTERED (unit mismatch, Appendix E.4)).
  flip AUC_frac is gridded on d/n_k (fraction of claims reversed); tiered AUC_frac_rate is gridded on corruption_rate (fraction of nodes relocated). Different axes, different units -- per Appendix E.4 the raw variance numbers above are reported for completeness only and are NOT a same-units comparison.
  var(flip AUC_frac) = 0.0640; var(tiered AUC_frac_rate) = 0.0139 (raw, non-comparable, tiered numerically higher: False). Matched (component_size, separation) cells: 0/12 numerically favor tiered (reported for completeness; not a same-units comparison; no verdict is forced).
 component_size  separation  n_flip  n_tiered  var_flip_AUC_frac  var_tiered_AUC_frac_rate  tiered_higher
              6           1     120        60           0.058112                  0.009954          False
              6           3     120        60           0.020541                  0.010412          False
              6           5     120        60           0.008964                  0.001098          False
              8           1     120        48           0.070736                  0.006413          False
              8           4     112        60           0.025341                  0.011335          False
              8           7     106        60           0.012040                  0.000547          False
             10           1     109        42           0.077784                  0.004804          False
             10           5     110        60           0.037938                  0.008330          False
             10           9     104        60           0.009805                  0.000324          False
             12           1     103        38           0.104364                  0.004948          False
             12           6      96        60           0.046279                  0.012073          False
             12          11      90        60           0.016360                  0.000635          False

## E6 -- sensitivity endpoint (AUC_frac_usable): do verdicts change?

- P1: primary=SUPPORTED, usable=SUPPORTED (8/9) -- UNCHANGED
- P2: primary=SUPPORTED, usable=SUPPORTED (7/9) -- UNCHANGED
- P3: primary=NOT SUPPORTED, usable=NOT SUPPORTED (2/9) -- UNCHANGED
- P4: no sensitivity variant post-Appendix-E (see above) -- the primary comparison is already NOT COMPARABLE by unit mismatch, independent of any usable-depth threshold.
- n instances with 0 usable (n_eval>=30) depths: 23 (1.2%); AUC_frac_usable undefined (NaN) for 23 instances. **Caveat:** for the tiered arm, AUC_frac_usable remains on the defective d_claims axis (no rate-binned "usable" variant was computed), so its P1-P3 rows above should be read as a d_claims-axis sensitivity check only, not as a check on AUC_frac_rate itself.

## Cross-check against hops_p5.csv (independent small-component subsample)

 coverage  base_wrongness  n_main  tau_main  ci_lo_2p5  ci_hi_97p5  n_hops  tau_hops  sign_agree  abs_diff  rough_magnitude_agree
      0.5            0.00     240  0.713720   0.625286    0.797966     210  0.863783        True  0.150063                   True
      0.5            0.10     240  0.379724   0.264552    0.488149     210  0.852911        True  0.473187                  False
      0.5            0.25     220  0.031381  -0.114416    0.168647     210  0.639799        True  0.608418                  False
      1.0            0.00     240  0.519242   0.441501    0.590725     210  0.687700        True  0.168458                   True
      1.0            0.10     236  0.387687   0.295891    0.477832     210  0.578777        True  0.191090                   True
      1.0            0.25     134  0.346081   0.209074    0.472125     206  0.375766        True  0.029685                   True

Sign agreement: 6/6 strata. Rough magnitude agreement (|main - hops| <= 0.35): 4/6.

## E1 -- decomposition identity

S_contra_as_fail(d) = (1 - contradiction_rate(d)) * S(d), checked on 8841 of 14066 curve rows (S undefined on 5225, all-contradictory depths, where the identity does not apply). max_abs_deviation=2.220e-16, mean_abs_deviation=2.325e-17, n_rows with deviation > 1e-6: 0.

## E2 -- determinism

Bootstrap CIs identical across two independent runs: True.

## E3 -- purity grep

```
/Users/ahn/Documents/Research/iclr27/bkrobust/src/bkrobust/robustness/analyse.py:7: ``all_valid_adjustment_sets_mpdag`` or ``synth.runner``. Reads
/Users/ahn/Documents/Research/iclr27/bkrobust/src/bkrobust/robustness/analyse.py:38: All randomness is explicit ``numpy.random.Generator`` -- no ``np.random.seed``,
/Users/ahn/Documents/Research/iclr27/bkrobust/src/bkrobust/robustness/analyse.py:39: no ``random.seed``, no bare ``np.random.*``.
/Users/ahn/Documents/Research/iclr27/bkrobust/run_analyse.py:30: r"all_valid_adjustment_sets_mpdag|synth\.runner|np\.random\.seed|random\.seed|import random"
```

## E4 -- stratum n and underpowered flags (primary strata)

   arm  base_wrongness  coverage  n_tiers        endpoint   n  underpowered
  flip            0.00       0.5      NaN        AUC_frac 240         False
  flip            0.00       0.5      NaN AUC_frac_usable 239         False
  flip            0.00       1.0      NaN        AUC_frac 240         False
  flip            0.00       1.0      NaN AUC_frac_usable 229         False
  flip            0.10       0.5      NaN        AUC_frac 240         False
  flip            0.10       0.5      NaN AUC_frac_usable 240         False
  flip            0.10       1.0      NaN        AUC_frac 236         False
  flip            0.10       1.0      NaN AUC_frac_usable 226         False
  flip            0.25       0.5      NaN        AUC_frac 220         False
  flip            0.25       0.5      NaN AUC_frac_usable 220         False
  flip            0.25       1.0      NaN        AUC_frac 134         False
  flip            0.25       1.0      NaN AUC_frac_usable 133         False
tiered             NaN       NaN      2.0   AUC_frac_rate 190         False
tiered             NaN       NaN      2.0 AUC_frac_usable 190         False
tiered             NaN       NaN      3.0   AUC_frac_rate 238         False
tiered             NaN       NaN      3.0 AUC_frac_usable 238         False
tiered             NaN       NaN      4.0   AUC_frac_rate 240         False
tiered             NaN       NaN      4.0 AUC_frac_usable 240         False

## E5 -- discordance spotlight verification

AUC_frac recomputed directly from survival_curves.csv matches the recorded value (within 1e-6) for 30/30 spotlighted instances (this check is on the d_claims axis AUC_frac column; it applies to both arms but for tiered instances is a plumbing check on AUC_frac, not on the AUC_frac_rate that actually drove their ranking -- see the next line).
                           instance_id                                           discordance_type  recorded_AUC_frac  recomputed_AUC_frac  auc_matches  first_S  last_S  n_curve_rows
     c08_s04_cov100_seed00000101_bw025 high_high (looks far from truth, survives deep corruption)           0.935185             0.935185         True 1.000000     1.0             8
     c12_s01_cov050_seed00000036_bw025 high_high (looks far from truth, survives deep corruption)           0.942105             0.942105         True 1.000000     1.0             6
c12_s06_cov100_seed00000010_tiered_nt4 high_high (looks far from truth, survives deep corruption)           1.000000             1.000000         True 1.000000     1.0            11
     c08_s04_cov100_seed00000128_bw025 high_high (looks far from truth, survives deep corruption)           0.964010             0.964010         True 1.000000     1.0             8
     c12_s01_cov050_seed00000023_bw025 high_high (looks far from truth, survives deep corruption)           0.942857             0.942857         True 1.000000     1.0             7
     c12_s01_cov050_seed00000031_bw010 high_high (looks far from truth, survives deep corruption)           0.905577             0.905577         True 0.760000     1.0             8
     c10_s05_cov100_seed00000036_bw025 high_high (looks far from truth, survives deep corruption)           1.000000             1.000000         True 1.000000     1.0             9
     c10_s05_cov100_seed00000140_bw025 high_high (looks far from truth, survives deep corruption)           1.000000             1.000000         True 1.000000     1.0             9
     c12_s01_cov050_seed00000020_bw010 high_high (looks far from truth, survives deep corruption)           0.915120             0.915120         True 0.656716     1.0             8
     c10_s05_cov100_seed00000009_bw025 high_high (looks far from truth, survives deep corruption)           0.975000             0.975000         True 1.000000     1.0            10
     c12_s11_cov100_seed00000134_bw025 high_high (looks far from truth, survives deep corruption)           1.000000             1.000000         True 1.000000     1.0            11
     c10_s05_cov100_seed00000059_bw025 high_high (looks far from truth, survives deep corruption)           1.000000             1.000000         True 1.000000     1.0            12
     c12_s06_cov100_seed00000094_bw025 high_high (looks far from truth, survives deep corruption)           1.000000             1.000000         True 1.000000     1.0            12
     c12_s06_cov100_seed00000127_bw025 high_high (looks far from truth, survives deep corruption)           1.000000             1.000000         True 1.000000     1.0            12
     c12_s01_cov100_seed00000040_bw025 high_high (looks far from truth, survives deep corruption)           0.925000             0.925000         True 1.000000     1.0            14
     c06_s01_cov050_seed00000007_bw000                 low_low (looks near truth, collapses fast)           0.000000             0.000000         True 0.000000     0.0             2
     c06_s01_cov050_seed00000007_bw010                 low_low (looks near truth, collapses fast)           0.000000             0.000000         True 0.000000     0.0             2
     c06_s01_cov050_seed00000007_bw025                 low_low (looks near truth, collapses fast)           0.000000             0.000000         True 0.000000     0.0             2
     c06_s01_cov050_seed00000013_bw000                 low_low (looks near truth, collapses fast)           0.000000             0.000000         True 0.000000     0.0             3
     c06_s01_cov050_seed00000013_bw010                 low_low (looks near truth, collapses fast)           0.000000             0.000000         True 0.000000     0.0             3
     c08_s01_cov050_seed00000007_bw000                 low_low (looks near truth, collapses fast)           0.128155             0.128155         True 0.427184     0.0             4
     c08_s01_cov050_seed00000013_bw010                 low_low (looks near truth, collapses fast)           0.129032             0.129032         True 0.430108     0.0             4
     c06_s01_cov050_seed00000016_bw010                 low_low (looks near truth, collapses fast)           0.227273             0.227273         True 0.454545     0.0             3
     c08_s01_cov050_seed00000013_bw000                 low_low (looks near truth, collapses fast)           0.135789             0.135789         True 0.452632     0.0             4
     c06_s01_cov050_seed00000005_bw000                 low_low (looks near truth, collapses fast)           0.244094             0.244094         True 0.488189     0.0             3
     c08_s01_cov050_seed00000007_bw010                 low_low (looks near truth, collapses fast)           0.162000             0.162000         True 0.540000     0.0             4
     c06_s01_cov050_seed00000016_bw000                 low_low (looks near truth, collapses fast)           0.257246             0.257246         True 0.514493     0.0             3
     c06_s01_cov050_seed00000005_bw010                 low_low (looks near truth, collapses fast)           0.288194             0.288194         True 0.576389     0.0             3
     c06_s01_cov100_seed00000007_bw000                 low_low (looks near truth, collapses fast)           0.000000             0.000000         True 0.000000     0.0             5
     c06_s01_cov100_seed00000007_bw010                 low_low (looks near truth, collapses fast)           0.000000             0.000000         True 0.000000     0.0             5
Supplementary check for the 1 spotlighted tiered instances: their AUC_frac_rate in the spotlight table matches analysis_tiered_rebinned.csv's own value for the same instance_id for 1/1.

## Sanity checks

- r_val == UNREACHED count: 0 (task expects 0): CONFIRMED
- r_assumes constant across dataset: "Conjecture 2 (hence Anti-Exchange Case B, verified not proved)"

## Effect sizes (median of each arm's primary endpoint by predictor bucket, primary strata; `endpoint` column states which axis: AUC_frac for flip, AUC_frac_rate for tiered)

   arm  base_wrongness  coverage  n_tiers predictor bucket      endpoint   n  median_endpoint_value  mean_endpoint_value
  flip            0.00       0.5      NaN     r_val      1      AUC_frac  80               0.375138             0.434704
  flip            0.00       0.5      NaN     r_val      2      AUC_frac  29               0.700000             0.700000
  flip            0.00       0.5      NaN     r_val      3      AUC_frac  11               0.800000             0.791466
  flip            0.00       0.5      NaN     r_val      4      AUC_frac  65               0.800000             0.799838
  flip            0.00       0.5      NaN     r_val     5+      AUC_frac  55               0.900000             0.879358
  flip            0.10       0.5      NaN     r_val      1      AUC_frac  65               0.411884             0.429974
  flip            0.10       0.5      NaN     r_val      2      AUC_frac  44               0.700000             0.655944
  flip            0.10       0.5      NaN     r_val      3      AUC_frac  11               0.800000             0.789880
  flip            0.10       0.5      NaN     r_val      4      AUC_frac  68               0.800000             0.793087
  flip            0.10       0.5      NaN     r_val     5+      AUC_frac  52               0.700000             0.746154
  flip            0.25       0.5      NaN     r_val      1      AUC_frac  41               0.254918             0.321890
  flip            0.25       0.5      NaN     r_val      2      AUC_frac  63               0.700000             0.660664
  flip            0.25       0.5      NaN     r_val      3      AUC_frac  84               0.600000             0.594482
  flip            0.25       0.5      NaN     r_val      4      AUC_frac  31               0.500000             0.561890
  flip            0.25       0.5      NaN     r_val     5+      AUC_frac   1               0.633117             0.633117
  flip            0.00       1.0      NaN     r_val      1      AUC_frac  80               0.260360             0.361691
  flip            0.00       1.0      NaN     r_val      3      AUC_frac  20               0.572486             0.633927
  flip            0.00       1.0      NaN     r_val      4      AUC_frac  20               0.535390             0.568454
  flip            0.00       1.0      NaN     r_val     5+      AUC_frac 120               0.900000             0.820097
  flip            0.10       1.0      NaN     r_val      1      AUC_frac  47               0.204009             0.382638
  flip            0.10       1.0      NaN     r_val      2      AUC_frac  37               0.346571             0.458278
  flip            0.10       1.0      NaN     r_val      3      AUC_frac  18               0.500000             0.518877
  flip            0.10       1.0      NaN     r_val      4      AUC_frac  17               0.500000             0.549364
  flip            0.10       1.0      NaN     r_val     5+      AUC_frac 117               0.800000             0.746912
  flip            0.25       1.0      NaN     r_val      1      AUC_frac  18               0.195398             0.288492
  flip            0.25       1.0      NaN     r_val      2      AUC_frac  37               0.354902             0.474268
  flip            0.25       1.0      NaN     r_val      3      AUC_frac  21               0.400866             0.521620
  flip            0.25       1.0      NaN     r_val      4      AUC_frac  33               0.700000             0.630239
  flip            0.25       1.0      NaN     r_val     5+      AUC_frac  25               0.600000             0.664381
tiered             NaN       NaN      2.0     r_val      1 AUC_frac_rate  58               0.701295             0.676229
tiered             NaN       NaN      2.0     r_val      2 AUC_frac_rate  39               0.819560             0.804019
tiered             NaN       NaN      2.0     r_val      3 AUC_frac_rate  31               0.913733             0.903904
tiered             NaN       NaN      2.0     r_val      4 AUC_frac_rate  22               0.951542             0.949903
tiered             NaN       NaN      2.0     r_val     5+ AUC_frac_rate  40               0.979040             0.971223
tiered             NaN       NaN      3.0     r_val      1 AUC_frac_rate  78               0.778168             0.777478
tiered             NaN       NaN      3.0     r_val      2 AUC_frac_rate  33               0.921626             0.919176
tiered             NaN       NaN      3.0     r_val      3 AUC_frac_rate  28               0.975025             0.965042
tiered             NaN       NaN      3.0     r_val      4 AUC_frac_rate  37               0.971528             0.972689
tiered             NaN       NaN      3.0     r_val     5+ AUC_frac_rate  62               1.000000             0.998479
tiered             NaN       NaN      4.0     r_val      1 AUC_frac_rate  80               0.812048             0.794068
tiered             NaN       NaN      4.0     r_val      3 AUC_frac_rate  40               0.960954             0.957705
tiered             NaN       NaN      4.0     r_val      4 AUC_frac_rate  33               1.000000             0.990844
tiered             NaN       NaN      4.0     r_val     5+ AUC_frac_rate  87               1.000000             0.998896
  flip            0.00       0.5      NaN shd_truth      0      AUC_frac  25               0.800000             0.624338
  flip            0.00       0.5      NaN shd_truth      1      AUC_frac  55               0.700000             0.570224
  flip            0.00       0.5      NaN shd_truth      2      AUC_frac  51               0.800000             0.687851
  flip            0.00       0.5      NaN shd_truth      3      AUC_frac  59               0.700000             0.690652
  flip            0.00       0.5      NaN shd_truth      4      AUC_frac   9               0.778909             0.712441
  flip            0.00       0.5      NaN shd_truth     5+      AUC_frac  41               0.900000             0.851829
  flip            0.10       0.5      NaN shd_truth      0      AUC_frac  20               0.555996             0.558943
  flip            0.10       0.5      NaN shd_truth      1      AUC_frac  41               0.787059             0.637420
  flip            0.10       0.5      NaN shd_truth      2      AUC_frac  35               0.700000             0.625946
  flip            0.10       0.5      NaN shd_truth      3      AUC_frac  57               0.700000             0.733436
  flip            0.10       0.5      NaN shd_truth      4      AUC_frac  17               0.589210             0.511268
  flip            0.10       0.5      NaN shd_truth     5+      AUC_frac  70               0.700000             0.692998
  flip            0.25       0.5      NaN shd_truth      0      AUC_frac   5               0.480000             0.401000
  flip            0.25       0.5      NaN shd_truth      1      AUC_frac   7               0.700000             0.700000
  flip            0.25       0.5      NaN shd_truth      2      AUC_frac  18               0.650000             0.624482
  flip            0.25       0.5      NaN shd_truth      3      AUC_frac  51               0.700000             0.613106
  flip            0.25       0.5      NaN shd_truth      4      AUC_frac  25               0.536901             0.482977
  flip            0.25       0.5      NaN shd_truth     5+      AUC_frac 114               0.600000             0.537887
  flip            0.00       1.0      NaN shd_truth      0      AUC_frac 240               0.637503             0.630811
  flip            0.10       1.0      NaN shd_truth      0      AUC_frac  39               0.900000             0.706406
  flip            0.10       1.0      NaN shd_truth      2      AUC_frac 197               0.600000             0.575930
  flip            0.25       1.0      NaN shd_truth      2      AUC_frac  40               0.700000             0.592894
  flip            0.25       1.0      NaN shd_truth      4      AUC_frac  80               0.449643             0.491424
  flip            0.25       1.0      NaN shd_truth     5+      AUC_frac  14               0.450433             0.576609
tiered             NaN       NaN      2.0 shd_truth      2 AUC_frac_rate  47               0.756631             0.786774
tiered             NaN       NaN      2.0 shd_truth      3 AUC_frac_rate  46               0.821843             0.815377
tiered             NaN       NaN      2.0 shd_truth      4 AUC_frac_rate  41               0.920094             0.873824
tiered             NaN       NaN      2.0 shd_truth     5+ AUC_frac_rate  56               0.876017             0.857738
tiered             NaN       NaN      3.0 shd_truth      1 AUC_frac_rate  48               0.957111             0.905897
tiered             NaN       NaN      3.0 shd_truth      2 AUC_frac_rate  98               0.973261             0.916217
tiered             NaN       NaN      3.0 shd_truth      3 AUC_frac_rate  66               0.973533             0.912842
tiered             NaN       NaN      3.0 shd_truth      4 AUC_frac_rate  21               0.831461             0.848790
tiered             NaN       NaN      3.0 shd_truth     5+ AUC_frac_rate   5               0.954545             0.909602
tiered             NaN       NaN      4.0 shd_truth      0 AUC_frac_rate  52               0.941601             0.900616
tiered             NaN       NaN      4.0 shd_truth      1 AUC_frac_rate  98               0.997043             0.930279
tiered             NaN       NaN      4.0 shd_truth      2 AUC_frac_rate  74               0.990570             0.928655
tiered             NaN       NaN      4.0 shd_truth      3 AUC_frac_rate  13               1.000000             0.928211
tiered             NaN       NaN      4.0 shd_truth      4 AUC_frac_rate   3               0.858921             0.882946
  flip            0.00       0.5      NaN       n_k      2      AUC_frac  34               0.700000             0.658235
  flip            0.00       0.5      NaN       n_k      3      AUC_frac  20               0.800000             0.603591
  flip            0.00       0.5      NaN       n_k      4      AUC_frac  87               0.800000             0.708907
  flip            0.00       0.5      NaN       n_k     5+      AUC_frac  99               0.900000             0.686968
  flip            0.10       0.5      NaN       n_k      2      AUC_frac  34               0.700000             0.658824
  flip            0.10       0.5      NaN       n_k      3      AUC_frac  20               0.800000             0.604333
  flip            0.10       0.5      NaN       n_k      4      AUC_frac  88               0.800000             0.708909
  flip            0.10       0.5      NaN       n_k     5+      AUC_frac  98               0.700000             0.626100
  flip            0.25       0.5      NaN       n_k      2      AUC_frac  42               0.700000             0.664405
  flip            0.25       0.5      NaN       n_k      3      AUC_frac  11               0.246154             0.235506
  flip            0.25       0.5      NaN       n_k      4      AUC_frac  95               0.600000             0.552763
  flip            0.25       0.5      NaN       n_k     5+      AUC_frac  72               0.500000             0.552771
  flip            0.00       1.0      NaN       n_k     5+      AUC_frac 240               0.637503             0.630811
  flip            0.10       1.0      NaN       n_k     5+      AUC_frac 236               0.660604             0.597491
  flip            0.25       1.0      NaN       n_k     5+      AUC_frac 134               0.598889             0.530613
tiered             NaN       NaN      2.0       n_k      1 AUC_frac_rate 111               0.934331             0.886669
tiered             NaN       NaN      2.0       n_k      2 AUC_frac_rate  29               0.795807             0.774755
tiered             NaN       NaN      2.0       n_k      3 AUC_frac_rate  30               0.743297             0.751816
tiered             NaN       NaN      2.0       n_k      4 AUC_frac_rate  11               0.739945             0.750708
tiered             NaN       NaN      2.0       n_k     5+ AUC_frac_rate   9               0.723471             0.738376
tiered             NaN       NaN      3.0       n_k      2 AUC_frac_rate  92               1.000000             0.976818
tiered             NaN       NaN      3.0       n_k      3 AUC_frac_rate  29               0.922944             0.878326
tiered             NaN       NaN      3.0       n_k      4 AUC_frac_rate  33               0.928305             0.886298
tiered             NaN       NaN      3.0       n_k     5+ AUC_frac_rate  84               0.845927             0.848881
tiered             NaN       NaN      4.0       n_k      3 AUC_frac_rate  84               1.000000             0.993864
tiered             NaN       NaN      4.0       n_k      4 AUC_frac_rate  17               0.943279             0.888837
tiered             NaN       NaN      4.0       n_k     5+ AUC_frac_rate 139               0.887017             0.883745
