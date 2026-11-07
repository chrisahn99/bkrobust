# E9 preregistration: parametric bootstrap of validity under a fitted model

Dated 2026-11-06. Not to be edited after the first full run; later changes go under "Deviations" with a date.

**Question.** On the learned-CPDAG panel (PC, Fisher-z, alpha 0.01, n = 5000, four Gaussian networks), does simulating from a linear-Gaussian model fitted to one DAG of the analyst's state flag the rows whose committed set Z is invalid in the true DAG better than the plain row bootstrap of E6?

**Population.** The 244 certifiable rows of `results/pest/rows.csv`, rebuilt exactly as in E6 (arms D_LLM, D_RAND reps 0-2, D_DEGEN, A_TRUE, A_ORACLE); a rebuild that disagrees with the frozen row stops the run.

**Fitted model.** Per row group (network, arm, rep) the analyst's state is G0 = the reference CPDAG with the arm's K imposed. Its simulated truth D_g is one consistent DAG extension of G0 (Dor-Tarsi, sorted tie-break). Coefficients and residual variances are fitted by least squares without intercept on the reference n = 5000 sample, variance = RSS/n (Gaussian maximum likelihood). All D_g of a network extend the same CPDAG, so their fitted covariances must coincide; the run checks this (relative difference < 1e-6) and simulates once per network from D_0, the extension of the D_DEGEN state. Z must be GAC-valid in every D_g, or the run stops.

**Procedure.** Per network, S = 200 datasets of n = 5000 drawn from the fitted D_0, seed `[20261109, stable_seed(network), b]`. On each: PC as in the panel, conversion, Meek closure (the raw graph if closure fails, as the panel did), the same K imposed with `apply_orientations`, then the committed set Z_b recomputed on that state with the panel's rule. A PC run past 900 s is a timeout and leaves the denominator. One identity unit per network runs PC on the reference sample and must give back the reference CPDAG.

**Registered statistic.** Per row, parametric validity v = share of finished simulations in which Z_b exists and is GAC-valid in D_g (no set, an inconsistent K or a failed check count as not valid). Score = 1 - v. AUC of the score for "Z invalid at truth" (Mann-Whitney, ties one half) over the rows, and the paired difference Delta = AUC(E9) - AUC(E6 primary score) over the rows present in both. 95% intervals: two-stage cluster bootstrap, 2000 draws, seed 20261106 (networks, then query pairs within a network, every arm's row of a drawn pair kept), the same draws for both AUCs.

**Secondary, descriptive only.** AUC of the parametric analogue of E6 (fixed Z GAC-valid on the simulated state); AUC of "Z_b differs from Z"; AUC of the conditional rate (Z_b invalid in D_g among simulations with a Z_b); AUC of the score averaged over M = 5 random consistent extensions of G0 (same simulated data, different simulated truths); per network and per arm AUCs; share of invalid rows with v >= 0.95.

**Prediction.** The fitted model treats the analyst's state as the truth, so it can only see errors that sampling noise around that state would cause, not the edges PC missed or misoriented on the real sample. We predict AUC(E9) between 0.5 and 0.7 and Delta between -0.05 and +0.05.

**Against the authors.** Delta < 0.05 (point estimate): the parametric bootstrap adds cost without information, it is left out of the response, and the response says that a bootstrap of either kind measures sampling variability around the learned state and does not certify the equivalence class. Delta >= 0.05 with an interval above 0 is reported beside E6 as the better companion check.

**When it runs.** Only after E6's full run has a results.json with a non-null AUC (the launcher checks). If E6 is not finished, Delta is not computed and only AUC(E9) is printed, marked provisional.

**Stopping rule.** S = 200 per network, fixed in advance; no extension or rerun because of the result. Timeouts are not rerun. A network with fewer than 160 finished simulations is reported as such and kept.

**Exclusions.** None beyond the non-certifiable rows. **Criterion.** GAC throughout (`gac.mpdag_level` on states, `gac.dag_level` on DAGs).

## Deviations from the plan

- Plan: "fit one DAG of the analyst's state, simulate 200 datasets". Fitting one DAG per row group would cost 28 x 200 PC runs; because Markov-equivalent Gaussian DAGs have the same maximum-likelihood fit, one simulation per network serves every group, checked numerically, and the cost stays at 2x E6.
