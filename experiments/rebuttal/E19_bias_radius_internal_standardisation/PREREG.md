# E19 preregistration: the 900-instance bias-radius study under internal standardisation

Dated 2026-11-06. Written before any run; not to be edited after the first run (later notes go in a separate file).

Question. Reviewer C2Bd (Q1) asks whether standardising variables changes the Section 7 / Appendix E.4 bias results. Do the App. E.4 shares (threshold never reached on 70.1% of 900 instances at epsilon = 100%; staircase a step on 84.1%) hold when the linear SEM is internally standardised (iSCM, Ormaniec et al.)?

Units. The 900 instances with status "ok" in results/epsilon/study/results.csv (the App. E.4 population), rebuilt from their stored index with experiments/run_r_epsilon.py's grid and seed (20260916 + index). Graph, treatment, outcome, knowledge, G0 and Z are unchanged; the check that they match the stored row is a consistency gate.

Arms, all from the same coefficient draw (submitted random_sem, |w| in U(0.3, 1.5), noise variance in U(0.5, 1.5)):
- naive: the draw as is (control; must reproduce the stored row: r_val, r_eps at the seven tolerances, r_0, staircase to 1e-9 relative).
- posthoc: every variable rescaled to unit variance after the fact. Relative bias is scale-free, so every r_eps and the shape must equal naive; any difference is a bug.
- internal: iSCM, each variable standardised in topological order before its children are generated, in closed form (population variances, no reference sample). Every implied variance must be 1 to 1e-9.

Criterion: generalised adjustment criterion ("gac", the certify default); recorded in every row.

Registered statistic (internal arm, instances "ok" in all three arms):
(1) share with r_eps = UNREACHED at epsilon = 100%; (2) share of staircases classified "step" among shape != "n/a", with the classifier of run_r_epsilon.py. Each with a Wilson 95% interval, plus the paired difference internal minus naive with a 95% bootstrap interval (2000 resamples of instances, seed 0).
Secondary, per arm: the never-reached share per tolerance, the share with r_eps > r_val per tolerance (App. E.4 definition), ramp share, share with r_0 = r_val, certificate-holds share, population var-sortability and R2-sortability per generator.

Prediction. r_val is identical across arms by construction. The never-reached share at 100% stays at or above 50% and the step share at or above 70% under internal standardisation; posthoc equals naive exactly.

Against the authors. If the never-reached share at 100% falls below 40% or the step share below 60%, App. E.4 will print these sentences for each generator (naive and internally standardised side by side) and say the ceiling is generator-dependent. A posthoc difference is a bug, reported as such before any other number.

Exclusions. An instance that times out (default 60 s per arm, raised from the study's 8 s because the population is now fixed by the stored rows rather than selected by the cap) or errors in any arm is excluded from the paired statistic and counted. If more than 5% are excluded, the statistic is reported on the reduced set with that fraction stated.

Stopping rule. One full pass over the 900 instances; no reruns with other seeds or tolerances to change the result. A consistency failure stops the reporting until explained.

## Deviations from the plan

- The plan names "the --sem switch from E18". E18's code did not exist when this was written, so the three SEM variants are implemented here (sem_variants.py), in closed form. The pilot iSCM constructor (pilot tree, rho_breakdown/pilot/code/scm.py) estimates standard deviations from a 100,000-row reference sample and draws |w| in U(0.5, 2.0); here they are exact population values and the submitted coefficient range is kept so the arms differ only by the standardisation.
- The study's side analyses that do not bear on the question (incumbent mean-bias comparison, shell-skip saving, greedy bound) are not rerun.
- The incremental and bisection strategies are both run at the largest tolerance as a consistency check, as in the study.
