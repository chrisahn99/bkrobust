# E6 preregistration: bootstrap stability as a companion check for discovery error

Dated 2026-11-06. Not to be edited after the first full run; later changes go under "Deviations" with a date.

**Question.** On the learned-CPDAG panel (PC, Fisher-z, alpha 0.01, n = 5000, four Gaussian networks), does the share of bootstrap resamples in which the committed adjustment set Z stays valid flag the rows whose Z is invalid in the true DAG?

**Population.** The 244 certifiable rows of `results/pest/rows.csv` (status `ok`; arms D_LLM, D_RAND reps 0-2, D_DEGEN, A_TRUE, A_ORACLE; networks ecoli70, arth150, magic-niab, magic-irri). Each row is rebuilt from the reference CPDAG and its arm's knowledge K; the rebuilt K hash, Z size and validity at truth must equal the frozen row, or the run stops.

**Procedure.** Per network, B = 100 nonparametric resamples (rows drawn with replacement) of the reference n = 5000 sample, seed `[20261106, stable_seed(network), b]`. On each: PC as in the panel, conversion, Meek closure (the raw graph if closure fails, as the panel did), then the same K imposed with `apply_orientations`. Z is valid on a replicate when that state exists and Z passes the generalised adjustment criterion (GAC) on it (`is_gac_valid_mpdag`). An inconsistent K or a failed check counts as not valid. A PC run past 900 s is a timeout and leaves the denominator. One extra identity replicate (no resampling) per network must give back the reference CPDAG hash exactly.

**Registered statistic.** Per row, stability s = valid replicates / replicates with a finished PC. Score = 1 - s. AUC of the score for "Z invalid at truth" (Mann-Whitney, ties count one half), over the 244 rows. 95% interval: two-stage cluster bootstrap, 2000 draws, seed 20261106: resample networks with replacement, then query pairs (x, y) within each drawn network, keeping every arm's row of a drawn pair. Draws without both classes are dropped and counted.

**Secondary, descriptive only.** AUC per network and per arm; AUC of neighbourhood instability (the skeleton adjacencies of x and y, and the edge marks around x and y in the state, unchanged from the reference); AUC of "committed set on the replicate differs from Z"; share of invalid rows with s >= 0.95 ("stable but wrong"); share of valid rows with s < 0.5.

**Prediction.** Discovery errors at n = 5000 are mostly systematic: a weak edge PC misses on the sample it also misses on most resamples. We predict AUC between 0.55 and 0.75, and that at least a third of invalid rows have s >= 0.95.

**Against the authors.** AUC <= 0.6 (point estimate). The response then says bootstrap stability is a variance check that does not certify the equivalence class, and the paper does not offer it as a fix for discovery error. AUC >= 0.8 with an interval above 0.6 is reported as a useful companion check, with the stable-but-wrong share printed beside it.

**Stopping rule.** B = 100 per network, fixed in advance; no extension or rerun because of the result. Timeouts are not rerun. A network with fewer than 80 finished replicates is reported as such and kept.

**Exclusions.** None beyond the non-certifiable rows. A rebuild mismatch stops the run; it never drops a row.

**Criterion.** GAC throughout (state level: `gac.mpdag_level`; truth: `gac.dag_level`, checked against the frozen networkx audit).

## Deviations from the plan

- None yet. Planned as "rerun PC, rebuild the state with the same knowledge"; the radius is not recomputed on replicates (not in the specification).
