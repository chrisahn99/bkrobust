# E20 preregistration: Section 7 profiles under the UUMC generator (2026-11-06)

Written before any run. Not to be edited after the first run; later notes go in a separate file.

**Question.** Reviewer C2Bd-Q1 asks whether the mean-bias profiles of Section 7 depend on how coefficients are drawn. Do the profiles still separate queries that share a radius when the coefficients come from UUMC (Herman, Wahl, Ninad, Runge, CLeaR 2025, Algorithm 1), which removes var-sortability and does not show iSCM's reversed R2-sortability?

**Units.** The 49 queries of Table 6: the 39 `ok` rows of results/mean_table_relevant_v2/x_component plus the 10 `ok` rows of results/mean_table_relevant_v2/extended with finite r_val and nonzero bias. D_LLM knowledge, Z is the optimal set O(G0), and the radius uses the GAC. There are three arms and draws 0 to 9. The naive arm is the submitted random_sem with seed derive_seed(20260917 + k, net, x, y), so draw 0 is the stored instance. The internally standardised arm (iSCM) standardises naive draw k exactly, at the population level. The UUMC arm is Algorithm 1 in population form, seeded from SeedSequence([derive_seed(20260917 + k, net, x, y), 1]).

**Registered statistic.** For each arm and draw, S = the number of distinct r_mu signatures (the tuple of r_mu at 2, 10, 25, 40, 60 and 90 %) over the queries computed in that draw. Also C = distinct (signature, B_rel(Chat) rounded to 3 decimals), R = distinct (r_val, B_rel(Chat) rounded), and added = C - R. The headline is the median of S over the ten UUMC draws, with a 95 % percentile bootstrap interval over draws (4000 resamples, seed 0) and the min and max. It is set beside the iSCM median and the naive draw-0 value, which must be 22 (C 32, R 22, B alone 21).

**Secondary.** The Table 2 contrast holds in a draw when barley dg25->s2225 has a finite r_mu(60 %) that is smaller than r_mu(60 %) of both Schipf_2010 A->U and Kampen_2014 AFF->CDR, with UNR counted as infinity. We report the share of draws in which it holds, per arm. We also report pair-counted R2-sortability and var-sortability (Herman et al. App. B, on the exact covariance) per arm, and the Reisach path-length version.

**Prediction.** Under UUMC, median S >= 10 and added > 0 in at least 8 of 10 draws. UUMC var-sortability is 0.5 by construction, and its R2-sortability sits closer to 0.5 than the naive arm's.

**Against the authors.** UUMC reverses the iSCM result if either of these holds:
(a) the Table 2 contrast holds in at least 6 of 10 draws under one of the two arms and in at most 4 of 10 under the other;
(b) median added is > 0 under one of the two arms and 0 under the other.
If so, every generator is reported side by side and the paper says the profile depends on the generator. If median S under UUMC is <= 3, the profile adds nothing to the radius under UUMC, and Section 7 says so.

**Consistency checks (any failure makes the run exit non-zero, and no number is reported):**
- naive draw 0 reproduces the stored theta_Z, profile mu, r_mu and B(Chat) to 1e-9;
- post hoc standardisation of naive draw 0 leaves every relative B unchanged to 1e-9;
- the enumerated and semilocal possible-parent sets agree on every state;
- cached B equals bias_at(method="enumerate") on sampled states;
- B = 0 at every distance below r_val (Theorem 1), in every arm and draw;
- the recomputed GAC radius equals the stored r_val;
- every UUMC and iSCM variance equals 1.

**Exclusions.** A query is excluded from a draw when |theta_Z| < 1e-12. A query whose state cache exceeds the budget (3600 s) is excluded from every arm. Both are counted and printed.

**Stopping rule.** One full run of 49 queries x 3 arms x 10 draws. No draws are added after the results are seen. After a consistency failure the bug is fixed and the whole run is redone, and the failure is reported.

## Deviations from the plan
- The plan says to add UUMC as a fourth arm of E18 through E18's --sem switch. That switch does not exist on the submitted branch. This run is therefore self-contained: it carries its own naive and iSCM comparators, so UUMC is set against an iSCM computed in exactly the same pipeline.
- iSCM here is the exact population version applied to the same naive weights (|w| in U(0.3, 1.5), noise variance in U(0.5, 1.5)). It is not the finite-sample reference standardisation with weights in U(0.5, 2) of the pilot tree's rho_breakdown/pilot/code/scm.py. If E18's iSCM uses that constructor, the two iSCM arms are not identical.
- UUMC is ported from the paper's Algorithm 1 (arXiv:2503.17037v3), not from the UUMCdata package, which is not installed. Here r is drawn as u^(1/d), and the noise is unit for roots.
