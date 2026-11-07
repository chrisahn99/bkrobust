# E14 preregistration: can a data test across sets valid in G0 flag a committed set invalid at truth? (2026-11-06)

Written before any run; not to be edited after the first full run.

**Question.** Reviewer qQYf (Q1) asks whether an invalid committed set can be detected. Su and Henckel (2022) test whether the OLS coefficients of X from several adjustment sets valid in a candidate graph agree. Applied to the knowledge-induced graph G0 of each released query, how often does it reject when the committed set Z is invalid in the true DAG (detection), and when Z is valid (false alarm)?

**Population.** Released rows of `per_query_records.csv` with `status = ok`, deduplicated on (network, x, y, k_sha256) so a state shared by several arms counts once. The state is rebuilt as the released sweep built it and checked against the row's hashes, |Z| and verdict; mismatches are excluded and counted. Candidate sets valid in G0 under the GAC: Z, the canonical set, pa_G0(x), minimal sets found by greedy removal from the canonical set (3 seeded orders), and Z with one element removed or one canonical element added; kept if GAC-valid in G0, at most 20, Z first. Units with fewer than two sets are excluded and counted.

**Data.** For each unit and replicate r (50 replicates), a linear SEM on the true DAG restricted to the ancestors of x, y and every candidate set: weights and noise variances from the submitted `random_sem` (|w| in [0.3, 1.5], random sign), rescaled node by node to unit marginal variance; each node's noise is one of uniform, Laplace or centred exponential, chosen by the seed, scaled to its variance. n = 1,000 and 10,000, independent samples, same SEM per replicate. Seeds: numpy `default_rng([20261106, unit hash, r])` for the SEM, `[20261106, unit hash, r, n]` for the sample.

**Test.** For each set, the OLS coefficient of x in y ~ 1 + x + Z_k and its influence function (Frisch-Waugh residual of x times the regression residual). The contrasts beta_k - beta_Z have a joint sandwich covariance from the influence functions; its rank is estimated by keeping eigenvalues above n^(-1/2) times the variance of Z's influence function (Su and Henckel's rank-estimation strategy). Wald statistic on the kept directions, chi-square with that many degrees of freedom, reject at alpha = 0.05. Sets with more than n/10 members are dropped at that n. A unit whose rank is 0 does not reject.

**Registered statistic.** Delta = D - F at n = 10,000, where D is the mean over units with Z GAC-invalid in the true DAG of the per-unit rejection rate over replicates, and F the same over units with Z GAC-valid in the true DAG. 95% interval: cluster bootstrap over units, 2,000 resamples, seed 20261106. Reported by released radius r_hop of the committed set (1, 2, 3 or more, no failure found) and for n = 1,000.

**Secondary.** Size on units where every compared set is valid in the true DAG (should be near 0.05); false alarms on units where Z is valid but some other set is not; the population ceiling, the share of invalid-Z replicates whose population coefficients differ across sets; a Bonferroni pairwise version of the test.

**Prediction.** D well above F at n = 10,000 (Delta >= 0.3), with D falling as the radius grows; a population ceiling below 1, because some wrong orientations bias every set equally.

**Against the authors.** Delta < 0.1 at n = 10,000. The authors then say that detection needs checking the named orientations or interventional data, as the Limitations will state. Hard failures exit non-zero and come before any number: a candidate set not GAC-valid in G0; the DAG-level and MPDAG-level GAC disagreeing on the true DAG; population coefficients that differ although every set is valid at truth, or a truth-valid Z whose population coefficient differs from the total effect.

**Stopping rule.** One pass over the population at 50 replicates and both n. No other n, alpha, replicate count or rank rule enters the registered statistic.

## Deviations from the plan

1. The plan takes the queries and their valid sets from E10, which does not exist yet; the candidate sets are generated here by the rule above, and the radius is the released r_hop of the committed set (computed by the released sweep, back-door criterion per E1), not recomputed per set.
2. The plan asks for a non-Gaussian sampler in `src/bkrobust/demo/evaluate.py`; it is written in this folder (`sem_sampler.py`) and reuses the submitted `random_sem` and `LinearSEM`.
3. No reference code of the test is used (no network on the compute machine); `coef_test.py` implements its rank-estimation strategy from the paper's description, so it may differ from the authors' RoCA package in details of the rank rule.
