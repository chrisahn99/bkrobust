# E18 preregistration: Section 7 profiles under standardised models (2026-11-06)

Written before the first run; not to be edited after it. Serves reviewer point C2Bd-Q1.

**Question.** Do the mean-bias profiles r_mu and the ceilings B(Chat)/|theta_Z| of Tables 2 and 6 change when the linear-Gaussian model is standardised, post hoc or internally (Ormaniec et al. 2025), or when the coefficients are drawn again?

**Population.** The 49 queries of Table 6 (39 rows of results/mean_table_relevant_v2/x_component plus the 10 rows of extended/ with finite r_val and non-zero bias), D_LLM knowledge, Z the optimal set of G0. No other query is added or removed.

**Arms, paired within a draw.** naive: random_sem of the submitted code (|w| in U(0.3, 1.5), random sign, noise variance U(0.5, 1.5)). posthoc: the covariance of the same draw rescaled to unit variances. internal: the same weights and noises, each variable divided by its exact population sd in topological order (the population limit of the iSCM constructor; unit variances by construction). Ten draws, d = 0..9. Draw 0 uses the stored seed derive_seed(20260917, net, x, y); draw d >= 1 uses derive_seed(20260917, "net|draw<d>", x, y).

**Unit.** One (query, draw, arm): its profile mu(d)/|theta_Z| over the states agreeing with G0 outside X's component, r_mu at the six tolerances 2, 10, 25, 40, 60, 90 %, the ceiling B(Chat)/|theta_Z|, the GAC radius, var- and R2-sortability. 1,470 units.

**Registered statistic.** For each (arm, draw): the number of distinct r_mu signatures over the 49 queries (signature = the six r_mu values, UNR kept as a value), and with the ceiling rounded to 3 decimals. Primary: the paired difference internal minus naive in distinct signatures, mean over the 10 draws with a 95 % t-interval. Secondary: the same for signature+ceiling; per-query agreement of signatures with the naive arm of the same draw; the two Table 2 contrasts per draw: (C1) barley dg25->s2225 crosses 2-60 % at one distance and 90 % later or never; (C2) Schipf_2010 A->U and Kampen_2014 AFF->CDR both reach 90 % later than dg25->s2225 (UNR counts as later).

**Prediction.** posthoc equals naive on every unit (signatures, rounded ceiling, profile and ceiling within 1e-9 relative, plus 1e-8 on the relative scale for bias_at's absolute 1e-9 zero clamp). The radius is identical across arms and draws (it reads no coefficients). Under internal standardisation and under redraws signatures move on some queries, but the count stays well above the radius's 3 distinct values; naive draw 0 reproduces 22 and 32.

**Against us.** Any posthoc difference is a bug: run.py exits non-zero and nothing is reported until it is fixed. If the median internal count over draws is at most half the median naive count, we say r_mu is generator-dependent, report it so in Appendix E.5, and keep the radius as the only generator-free result. If C1 or C2 fails in most draws of an arm, the Table 2 illustration is reported as draw-specific.

**Consistency checks (exit non-zero).** Naive draw 0 must reproduce every stored profile field to 1e-9, the stored b_top of the extended rows to 1e-9 and the printed Table 6 ceilings to 5e-4; the GAC radius must equal the stored r_val; the internal arm must have unit variances to 1e-9.

**Exclusions.** A query whose slice enumeration exceeds the per-query budget (default 3600 s) is recorded as timeout and listed; theta_Z = 0 in an arm excludes that unit, listed. No other exclusion.

**Stopping rule.** One full run of all 1,470 units. No reruns with other seeds, tolerances or draws after seeing results; a rerun is allowed only to fix a consistency failure, and is reported.

## Deviations from the plan
- The plan named experiments/mean_table_relevant_v2.py and mean_table_xcomp_extend_v2.py "with --sem, --draws". Those files are not modified; run.py imports their functions and caches the SEM-free part (slice, distances, possible parent sets of X) once per query.
- The internal arm uses the exact population sd instead of the 100,000-sample reference of scm.py, and the submitted weight generator instead of scm.py's U(0.5, 2.0), so the three arms differ only by standardisation.
