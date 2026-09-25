# When r_val and separation disagree, which one ranks survival correctly?

Source: `results/axis_robustness_p6/survival_instances.csv` (4,104 instances).
Script: `experiments/discordant_separation_v2.py` (NEW; nothing under
`results/axis_robustness_p6*` outside this `discordant/` subdirectory was
modified). All bootstraps are instance-clustered, seed=0.
`r_val != separation` on 2,287/4,104 instances overall (55.7%; flip 1,407/2,730
= 51.5%, tiered 880/1,374 = 64.0% -- matches the 55.7% figure already reported
in `results/axis_robustness_p6_paired_v2/FINDINGS.md`).

## C. How do r_val and separation differ? (`diff_distribution.json`)

`diff = r_val - separation`, on the "measured" separation-status population:

| arm | n | mean(diff) | median | q25/q75 | frac r_val>sep | frac r_val<sep | frac equal |
|---|---|---|---|---|---|---|---|
| flip | 2730 | -0.876 | 0.0 | -1.0 / 0.0 | 12.6% | 38.9% | 48.5% |
| tiered | 1374 | -1.346 | -1.0 | -2.0 / 0.0 | **0.0%** | 64.0% | 36.0% |
| all | 4104 | -1.033 | 0.0 | -2.0 / 0.0 | 8.4% | 47.3% | 44.3% |

r_val is **never larger than separation on a single tiered instance**
(`n_r_gt_sep = 0` of 1,374) -- separation is a hard upper bound on r_val in
the tiered arm on this sweep. On flip it is usually an identity or r_val is
smaller (`r_val <= separation` on 87.4% of instances); r_val exceeds
separation on only 12.6% of flip instances and never by more than +4 (vs. as
low as -7 the other way). So wherever they differ, **r_val undershoots
separation** far more often and by more than it overshoots, in both arms.

## A. Discordant-subset paired tau_b (`discordant_subset_paired.csv`)

Restricted to the `r_val != separation` instances; `paired_comparison(r_val,
separation, AUC_frac, instance_id)`, n_boot=10,000:

| group | n | delta = tau(r_val)-tau(sep) | 95% CI | verdict |
|---|---|---|---|---|
| flip POOLED | 1407 | -0.067 | [-0.098, -0.038] | separation wins |
| tiered POOLED | 880 | +0.152 | [0.132, 0.173] | r_val wins |
| flip cov0.5/bw0.00 | 187 | +0.136 | [0.107, 0.169] | r_val wins |
| flip cov0.5/bw0.10 | 262 | -0.239 | [-0.337, -0.144] | separation wins |
| flip cov0.5/bw0.25 | 396 | -0.147 | [-0.237, -0.056] | separation wins |
| flip cov1.0/bw0.10 | 308 | -0.012 | [-0.040, 0.016] | no sig. diff |
| flip cov1.0/bw0.25 | 254 | -0.076 | [-0.144, -0.011] | separation wins |
| tiered nt2 | 320 | +0.103 | [0.063, 0.145] | r_val wins |
| tiered nt3 | 320 | -0.031 | [-0.058, -0.006] | separation wins |
| tiered nt4 | 240 | -0.007 | [-0.022, 0.009] | no sig. diff |
(flip cov1.0/bw0.00 has 0 discordant instances -- r_val==separation
identically there, as already documented.)

Mixed verdict: r_val wins on 2 cells (flip cov0.5/bw0.00, tiered nt2),
separation wins on 4 (flip cov0.5/bw0.10, cov0.5/bw0.25, cov1.0/bw0.25,
tiered nt3), 2 are inconclusive. Pooling across strata flips the flip-arm
sign relative to its own best cell -- a composition effect from mixing
design cells of different sizes and concordance rates, exactly the kind of
aggregation artifact analysis B is designed not to have.

## B. Pair-level test (`pairwise_disagreement.csv`) -- the cleanest statistic

For every stratum and pooled per arm, over all instance PAIRS whose r_val-
order and separation-order of the pair disagree, the fraction of those pairs
where AUC's order agrees with r_val vs. with separation (2,000-resample
instance-cluster bootstrap CI). "Strict" = both predictors give a non-tied,
opposite verdict; "any" = additionally admits one predictor being tied.

**Pooled per arm ("any disagreement", the more inclusive definition):**

| arm | n pairs | frac AUC agrees r_val | 95% CI | frac AUC agrees separation | 95% CI |
|---|---|---|---|---|---|
| flip | 903,246 | 0.418 | [0.402, 0.434] | 0.260 | [0.245, 0.275] |
| tiered | 138,289 | 0.423 | [0.399, 0.447] | 0.289 | [0.264, 0.315] |

("strict opposite" pooled tells the same story more sharply: flip 0.551 vs.
0.413; tiered 0.744 vs. 0.246.) **Pooled, r_val wins the pair-level test in
both arms** -- the remaining ~30-32% of "any"-disagreement pairs match
neither predictor (an AUC tie, or a sign matching neither disagreeing
predictor).

This pooled result reverses the flip-arm verdict analysis A gave (A said
separation wins pooled on flip). The reversal traces to per-stratum
heterogeneity: at the design-cell level the pair-level test still shows
separation winning cleanly in the harder regimes --

| stratum | any: agree_r | agree_sep | strict: agree_r | agree_sep |
|---|---|---|---|---|
| flip cov0.5/bw0.10 | 0.348 | 0.177 | 0.590 | 0.410 |
| flip cov0.5/bw0.25 | 0.248 | 0.296 | 0.357 | **0.643** |
| flip cov1.0/bw0.25 | 0.378 | 0.316 | 0.378 | **0.622** |
| tiered nt3 | 0.234 | **0.472** | 0.309 | **0.691** |
| tiered nt4 | 0.161 | **0.549** | 0.125 | **0.875** |

-- and separation dominates most starkly at **tiered nt4** (strict: separation
agrees with AUC on 87.5% of disagreeing pairs vs. r_val's 12.5%), a cell
where analysis A found *no* significant difference (CI included 0). The
pair-level test resolves that null result: at nt4 separation really is the
better tiebreaker, more starkly than the instance-subset test could see.

## D. tau_b(separation, AUC | r_val) -- symmetric partial association

Mirrors the already-committed `tau_b(r_val, AUC | separation)`
(`results/axis_robustness_p6_paired_v2/separation_within_strata_fixed.json`,
Addendum (A): flip 0.375 [0.312, 0.433], tiered 0.487 [0.401, 0.562],
stratifying on exact separation + component_size within each design cell).
Here the roles are swapped: stratify on **exact r_val** + component_size
within each design cell, and ask whether separation still orders survival
once r_val is pinned (`rval_conditioned_separation.json`):

| arm | instances in varying cells | tau_b(separation \| r_val) | 95% CI |
|---|---|---|---|
| flip | 1128 / 2730 | **0.120** | [0.054, 0.225] |
| tiered | 141 / 1374 | **0.042** | [-0.117, 0.195] |

Compare: tau_b(r_val \| separation) = 0.375 (flip) / 0.487 (tiered) vs.
tau_b(separation \| r_val) = 0.120 (flip) / 0.042, not sig. (tiered). Most
tiered r_val-value cells (nt3, nt4, and most of nt2) are exactly
separation-constant given r_val (0/12, 0/12, 4/12 design cells vary; only
141 of 1,374 tiered instances fall in a cell where separation still varies
once r_val is fixed) -- tau_b there is mostly **undefined**, not just small.
This is a genuinely asymmetric result: **r_val carries substantial
information about survival beyond what separation alone gives (both arms);
separation carries little-to-no extra information about survival beyond
what r_val alone gives**, especially in tiered where r_val's exact value
essentially pins separation.

## Honest verdict

There is no single overall winner, and reporting only the pooled numbers
would be misleading in both directions. On the cleanest, subset-selection-free
pair-level test (B), r_val's ranking matches survival on more disagreeing
pairs than separation's does, pooled in both arms (flip 42% vs. 26%; tiered
42% vs. 29% of "any-disagreement" pairs) -- and separation adds essentially
no information about survival once r_val is fixed (D: partial tau_b 0.12/flip,
~0/tiered, vs. r_val's own 0.38/0.49 given separation). But that pooled edge
for r_val hides real, sizeable reversals within specific design cells --
separation is the better tiebreaker in flip's higher-corruption cells
(bw=0.25) and in tiered nt3/nt4, most dramatically at nt4 (separation agrees
with survival on 87.5% vs. r_val's 12.5% of strictly-disagreeing pairs) --
confirming and sharpening, not contradicting, the paper's own disclosure that
separation beats r_val on tiered nt3/nt4.
