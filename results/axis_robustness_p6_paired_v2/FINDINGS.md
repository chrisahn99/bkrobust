# Follow-up analysis on the committed synthetic sweep (v2)

Source: `results/axis_robustness_p6/survival_instances.csv` (4,104 instances),
`survival_curves.csv`, `survival_samples.csv.gz`. Committed paired bootstrap:
`results/axis_robustness_p6_paired/paired_diff_synth.csv`. All outputs here are
NEW, written by `experiments/paired_synth_v2.py`; nothing under
`results/axis_robustness_p6*` was modified. `n_boot=10000`, `seed=0` throughout,
matching the committed run.

## 0. Reproduction check

Rerunning `paired_comparison` for the five UNCHANGED baselines (shd_truth, n_k,
k_g0, separation, r_claim) against the committed `paired_diff_synth.csv`:
**0 discrepancies** (`reproduction_check.json` is `[]`; every tau_a_point,
tau_b_point, delta_point, CI bound, frac_delta_gt_0 and status matches to
1e-9). r_claim censoring is handled identically to the committed run: censored
rows carry a blank `r_claim` cell (never a value), so they are dropped by the
same missing-value path `paired_resample._clean_matched` already uses; no
explicit `censoring_columns` mapping was needed in either run. Per-stratum n
used for r_claim (`manifest.json: r_claim_n_used_by_stratum`):

| stratum | stratum_n | n_r_claim_defined | n_censored(blank) |
|---|---|---|---|
| flip_cov0.5_bw0.0 | 480 | 240 | 240 |
| flip_cov0.5_bw0.1 | 480 | 240 | 240 |
| flip_cov0.5_bw0.25 | 476 | 379 | 97 |
| flip_cov1.0_bw0.0 | 480 | 200 | 280 |
| flip_cov1.0_bw0.1 | 480 | 212 | 268 |
| flip_cov1.0_bw0.25 | 334 | 196 | 138 |
| tiered_nt2 | 414 | 414 | 0 |
| tiered_nt3 | 480 | 480 | 0 |
| tiered_nt4 | 480 | 480 | 0 |

## 1. Direction fix: Delta = taub(r_val) - taub(neg_phi_1), neg_phi_1 = -phi_1

`paired_diff_synth_v2.csv`, n=4104 total (per-stratum n above). Positive Delta
means r_val still outranks (-phi_1) once phi_1's fragility sign is corrected.

| stratum | n | taub(r_val) | taub(-phi_1) | Delta | 95% CI |
|---|---|---|---|---|---|
| flip cov0.5 bw0.00 | 480 | 0.614 | 0.379 | +0.235 | [0.197, 0.275] |
| flip cov0.5 bw0.10 | 480 | 0.328 | 0.346 | -0.017 | [-0.077, 0.040] (CI includes 0) |
| flip cov0.5 bw0.25 | 476 | -0.027 | 0.342 | -0.369 | [-0.444, -0.296] |
| flip cov1.0 bw0.00 | 480 | 0.469 | 0.362 | +0.107 | [0.078, 0.137] |
| flip cov1.0 bw0.10 | 480 | 0.418 | 0.304 | +0.114 | [0.055, 0.174] |
| flip cov1.0 bw0.25 | 334 | 0.280 | 0.300 | -0.020 | [-0.096, 0.060] (CI includes 0) |
| tiered nt2 | 414 | 0.819 | -0.568 | +1.387 | [1.340, 1.431] |
| tiered nt3 | 480 | 0.807 | 0.534 | +0.274 | [0.239, 0.309] |
| tiered nt4 | 480 | 0.748 | 0.616 | +0.132 | [0.104, 0.161] |

Once the direction bug is fixed, r_val beats -phi_1 with a CI excluding 0 in
6 of 9 strata (all defined tiered strata, and 3 of 6 flip strata); 3 flip
strata are inconclusive (CI straddles 0), and r_val never significantly loses.
This replaces the committed table's (wrongly-signed) `taub(r_val) -
taub(phi_1)` comparison, which read as always-positive for the wrong reason.

## Delta vs. separation (disclosed baseline)

| stratum | n | taub(r_val) | taub(separation) | Delta | 95% CI |
|---|---|---|---|---|---|
| flip cov0.5 bw0.00 | 480 | 0.614 | 0.461 | +0.153 | [0.118, 0.190] |
| flip cov0.5 bw0.10 | 480 | 0.328 | 0.245 | +0.083 | [0.040, 0.128] |
| flip cov0.5 bw0.25 | 476 | -0.027 | 0.005 | -0.032 | [-0.097, 0.032] (CI includes 0) |
| flip cov1.0 bw0.00 | 480 | 0.469 | 0.469 | 0.000 | identity (r_val == separation exactly here) |
| flip cov1.0 bw0.10 | 480 | 0.418 | 0.419 | -0.000 | [-0.019, 0.018] (CI includes 0) |
| flip cov1.0 bw0.25 | 334 | 0.280 | 0.312 | -0.033 | [-0.085, 0.018] (CI includes 0) |
| tiered nt2 | 414 | 0.819 | 0.734 | +0.085 | [0.056, 0.115] |
| tiered nt3 | 480 | 0.807 | 0.820 | -0.013 | [-0.025, -0.002] |
| tiered nt4 | 480 | 0.748 | 0.778 | -0.031 | [-0.044, -0.019] |

Confirms the review point: **separation beats r_val on tiered nt3 and nt4**
(CI entirely negative, -0.013 and -0.031), same direction the paper already
flags for the old sweep. On flip it is a wash (identity at cov1.0/bw0.00,
CI includes 0 elsewhere) or favours r_val (cov0.5 strata).

## Delta vs. r_claim

| stratum | n(matched) | taub(r_val, matched) | taub(r_claim) | Delta | 95% CI |
|---|---|---|---|---|---|
| flip cov0.5 bw0.00 | 240 | - | 0.369 | 0.000 | identity (r_val==r_claim on uncensored rows) |
| flip cov0.5 bw0.10 | 240 | - | 0.370 | 0.000 | identity |
| flip cov0.5 bw0.25 | 379 | - | 0.133 | -0.030 | [-0.072, 0.011] (CI includes 0) |
| flip cov1.0 bw0.00 | 200 | - | 0.257 | 0.000 | identity |
| flip cov1.0 bw0.10 | 212 | - | 0.326 | 0.000 | identity |
| flip cov1.0 bw0.25 | 196 | - | 0.293 | 0.000 | identity |
| tiered nt2 | 414 | 0.819 | undefined (r_claim constant) | n/a | status=undefined |
| tiered nt3 | 480 | 0.807 | 0.659 | +0.148 | [0.124, 0.173] |
| tiered nt4 | 480 | 0.748 | 0.714 | +0.034 | [0.012, 0.058] |

## 2. Separation facts on the current 4,104-instance sweep

`separation_status == "measured"` on all 4,104 instances (0 undefined).

**r_val != separation:**

| scope | n | n(r_val!=sep) | frac |
|---|---|---|---|
| overall | 4104 | 2287 | 55.7% |
| flip | 2730 | 1407 | 51.5% |
| tiered | 1374 | 880 | 64.0% |

By stratum: flip cov1.0/bw0.00 is an exact identity (0% differ, r_val==sep
by construction there); all other strata differ on 39%-83% of instances
(flip cov0.5/bw0.25: 83.2%, flip cov1.0/bw0.25: 76.0%, tiered nt2: 77.3%,
nt3: 66.7%, nt4: 50.0%; see `separation_facts.json:neq_facts.by_stratum` for
the full table). This 55.7% overall figure is the current-sweep analogue of
the paper's stale "56% of the 1,978 instances" -- coincidentally close in
magnitude but computed on a different (larger, current) corpus and not
otherwise validated against the old number.

**Within-separation-strata taub(r_val, AUC_frac).** Strata = design cell that
holds everything but separation fixed: (component_size, coverage,
base_wrongness) for flip [24 strata, all used], (component_size, n_tiers) for
tiered [12 strata, all used]. No script implementing the paper's original
combination method was found in this repo (grepped `src/bkrobust/robustness`
and `experiments` for "1978"/"56%"/"separation strata"/"stratified tau" --
nothing; the old 1,978-instance sweep and its analysis script are not present
here), so per the task brief's fallback this uses a pair-count-weighted
average of per-stratum taub (weight = n*(n-1)/2), with an instance-cluster
bootstrap CI recomputed the same way on each resample:

| arm | n | strata used | taub (weighted) | 95% CI |
|---|---|---|---|---|
| tiered | 1374 | 12/12 | 0.781 | [0.767, 0.794] |
| flip | 2730 | 24/24 | 0.494 | [0.455, 0.535] |

Both intervals exclude 0 comfortably; r_val's ranking of survival is not an
artefact of the separation design, on the current sweep. (These numbers are
NOT directly comparable to the paper's 0.55/0.24 on the old 1,978-instance
sweep -- different corpus, different stratum definition/count -- so this
section should be read as fresh current-sweep numbers, not a reproduction.)

Same computation for the other predictors (comparison only):

| arm | predictor | taub (weighted) | 95% CI |
|---|---|---|---|
| tiered | n_k | -0.595 | [-0.620, -0.570] |
| tiered | k_g0 | -0.418 | [-0.452, -0.384] |
| tiered | shd_truth | -0.271 | [-0.309, -0.231] |
| tiered | neg_phi_1 | +0.260 | [0.199, 0.318] |
| flip | n_k | -0.219 | [-0.260, -0.178] |
| flip | k_g0 | -0.247 | [-0.284, -0.211] |
| flip | shd_truth | +0.021 | [-0.020, 0.062] (CI includes 0) |
| flip | neg_phi_1 | +0.378 | [0.336, 0.418] |

r_val's within-strata taub (0.781 tiered, 0.494 flip) exceeds every other
predictor's including neg_phi_1, within separation strata specifically
(unlike the raw per-stratum paired table above, where separation itself
edges out r_val on tiered nt3/nt4).

## 3. AUC support (reviewer 4.4)

Full support = `n_depths_all_contradictory == 0` (no grid point where every
draw was contradictory).

| arm | n | full support | some undefined gridpoint | frac undefined | defined gridpoints (mean/median/min/max) |
|---|---|---|---|---|---|
| flip | 2730 | 853 (31.2%) | 1877 (68.7%) | 68.7% | 4.65 / 4 / 1 / 11 |
| tiered | 1374 | 0 (0.0%) | 1374 (100%) | 100% | 5.01 / 5 / 3 / 10 |

**Every tiered instance has at least one grid point where all draws were
contradictory** -- full-support restriction eliminates the tiered arm
entirely (`n_full_support_total = 853`, all from flip). Per-stratum
full-support fraction for flip collapses sharply with intensity:
cov0.5/bw0.00 66.9% down to cov1.0/bw0.25 0.9%.

Re-running the paired comparison restricted to the 853 full-support (flip
only) instances (`full_support_paired.csv`):

| stratum | n | Delta vs separation | Delta vs neg_phi_1 |
|---|---|---|---|
| flip cov0.5 bw0.00 | 321 | +0.444 [0.381, 0.511] | +0.573 [0.492, 0.672] |
| flip cov0.5 bw0.10 | 250 | +0.601 [0.516, 0.686] | +0.506 [0.408, 0.620] |
| flip cov0.5 bw0.25 | 86 | +0.032 [-0.392, 0.386] (CI includes 0) | -0.082 [-0.372, 0.135] (CI includes 0) |
| flip cov1.0 bw0.00 | 150 | 0.000 (identity) | undefined (neg_phi_1 constant) |
| flip cov1.0 bw0.10 | 43 | +0.004 [0.000, 0.013] | undefined (neg_phi_1 constant) |
| flip cov1.0 bw0.25 | 3 | 0.000 (n too small) | undefined |

**Conclusion: cannot be checked for tiered** (0 full-support instances there,
so the separation-beats-r_val result on tiered nt3/nt4 has no full-support
analogue to test). On flip, where full support exists, r_val's edge over
separation and neg_phi_1 survives at low intensity (cov0.5 strata) and is
either an identity or too small-n to resolve at high intensity -- consistent
with (not contradicting) the whole-sample conclusions, but the flip-only,
small-n subset (down to n=3 at cov1.0/bw0.25) limits how much weight this
robustness check can bear.

Accepted (consistent, non-contradictory) draws per intensity, from
`survival_curves.csv`'s own `n_eval` (median / min across instances), full
table in `auc_support.json:accepted_draws_per_intensity`; representative
points: flip d=1 median 1000 (max support), d=11 median ~40; tiered d=0
median ~904 declining to single digits by d=8-9.

## 4. Closure-orientation-change conditioning

`survival_samples.csv.gz` (32.79M per-draw rows) DOES carry the signal:
`closure_inert` (1 iff `n_closure_orientations_changed == 0`), logged only
for `corruption_accepted == 1` (consistent) draws -- confirmed against
`survival_p6.py`'s own docstring and column list. No corruption sweep was
rerun; this is a streaming aggregation (one `awk` pass, 64.6s) of the
existing file.

- **p0 = P(inert | consistent)**: 0.0 at every intensity for the flip arm
  -- every consistent flip draw changes at least one closure orientation.
  For tiered it is intensity-dependent and drops fast:
  d=0: 0.857, d=1: 0.114, d=2: 0.016, d=3: 0.0025, d=4+: approx 0.

- **Survival conditional on a changed state** (`survival_given_changed`):
  flip declines from 0.873 (d=1) to ~0.00002 (d=11), roughly tracking the
  unconditional survival curve since p0 approx 0 there means "changed" is
  nearly all consistent draws. Tiered similarly declines from 0.713 (d=0) to
  0.0 by d=9 (n=3).

- **taub(r_val, AUC of changed-only survival)**: flip 0.390 (n=2730), tiered
  0.542 (n=1374), pooled 0.363 (n=4104) -- positive and comparable in sign
  and rough magnitude to the whole-sample taub(r_val, AUC_frac) reported in
  the committed paired table, so conditioning on "closure state actually
  changed" does not overturn r_val's ranking.

Full per-arm/intensity tables in `changed_state_analysis.json`.

## Addendum: corrected within-separation strata (A) and AUC common-support check (B)

New files: `separation_within_strata_fixed.json`, `auc_grid_diagnostics.json`,
`auc_common_support_marginal.json`, `auc_common_support_paired.csv`. Script:
`experiments/paired_synth_v2_addendum.py`.

### (A) Within-separation strata, separation now held fixed

Corrected stratum = (paired-table stratum, component_size, separation) -- the
generator's own ~40-instance admission cell drawn at one exact target
separation (Section 2 above had this backwards: its "strata" let separation
vary). Only cells where r_val itself varies contribute a defined tau_b.

| arm | instances | cells total | cells with r_val variation | instances in those cells | taub(r_val, AUC), pair-weighted | 95% CI |
|---|---|---|---|---|---|---|
| flip | 2730 | 72 | 30 | 1132 | 0.375 | [0.312, 0.433] |
| tiered | 1374 | 36 | 8 | 320 | 0.487 | [0.401, 0.562] |

Per paired-table stratum (tiered only, as requested):

| stratum | cells varying / total | instances contributing | taub(r_val) | 95% CI |
|---|---|---|---|---|
| tiered nt2 | 4/12 | 160 | 0.608 | [0.519, 0.678] |
| tiered nt3 | 3/12 | 120 | 0.412 | [0.276, 0.535] |
| tiered nt4 | 1/12 | 40 | 0.224 | [-0.073, 0.499] (CI includes 0) |

Other predictors, pair-weighted, separation truly fixed:

| arm | predictor | taub | 95% CI |
|---|---|---|---|
| flip | n_k | 0.036 | [-0.034, 0.108] (CI incl. 0) |
| flip | k_g0 | 0.020 | [-0.024, 0.065] (CI incl. 0) |
| flip | shd_truth | 0.002 | [-0.051, 0.054] (CI incl. 0) |
| flip | neg_phi_1 | 0.179 | [0.123, 0.232] |
| tiered | n_k | -0.024 | [-0.087, 0.037] (CI incl. 0) |
| tiered | k_g0 | -0.037 | [-0.096, 0.021] (CI incl. 0) |
| tiered | shd_truth | -0.056 | [-0.112, 0.003] (CI incl. 0) |
| tiered | neg_phi_1 | -0.039 | [-0.107, 0.031] (CI incl. 0) |

With separation genuinely pinned, r_val is the only predictor here with a
CI clearly excluding 0 in both arms (n_k/k_g0/shd_truth/neg_phi_1 are all CI
incl. 0 except flip neg_phi_1, weaker at 0.179 vs r_val's 0.375). Caveat:
only 30/72 flip cells and 8/36 tiered cells have any r_val variation at all
(most ~40-instance cells are r_val-constant), so this result rests on a
minority of cells/instances (1132/2730 flip, 320/1374 tiered) -- reported
directly above, not hidden.

### (B) How AUC_frac is computed, and which grid points are undefined

`survival.auc_frac` (survival.py:559-583) is **not** a trapezoid and **not**
a mean over defined raw depths. It evaluates a fixed 10-point fractional grid
`FRAC_GRID = 0.1, 0.2, ..., 1.0` (fractions of the instance's own `n_k`), and
for each fraction takes `S` at the *nearest depth with a defined S* (ties
broken toward the smaller depth) -- so it is a per-instance normalized,
nearest-neighbour-imputed 10-point average. The raw intensity grid itself is
not comparable across instances: flip's `d` is a dense target `1..n_k`
(n_k varies by component_size/coverage); tiered's `d` is a *measured*
post-hoc claim-delta (irregular, not a swept knob -- the actual knob is
`corruption_rate` in `{0.00, 0.05, ..., 0.50}`). This normalization is why
`auc_frac` is defined for nearly every instance even though many raw grid
points are contradictory.

Undefined S is concentrated at high intensity, as expected: by depth rank
(1st, 2nd, ... smallest measured depth), `frac_undefined` rises from 0% at
rank 1 to 66% (flip) / 96% (tiered) by rank 8 (`auc_grid_diagnostics.json:
by_depth_rank`). Separately, the *exact* target depth for a given
`FRAC_GRID` fraction is itself rarely the one with defined S (needs the
nearest-neighbour fallback 45-100% of the time across fractions, both arms)
-- expected, since `round(frac*n_k)` is a specific integer and depths are
sparse; this is the imputation mechanism working as designed, not a support
failure per se.

**Common-support recompute.** Per paired-table stratum, took the intersection
of `{d : S(d) defined}` across every instance in the stratum (the 90%-prefix
fallback was never needed -- every stratum had a non-empty exact
intersection, 1-3 common depths), recomputed AUC via the same nearest-
neighbour algorithm restricted to that depth set, and reran
`tau_b(r_val, AUC_common)` and the paired deltas vs. separation, neg_phi_1,
n_k, shd_truth (**n_boot=2000, not 10,000** -- rebuilding the paired tau on
a resampled index set with a Python nearest-neighbour search per row made
10,000 too slow within budget; all `auc_common_support_paired.csv` rows carry
`n_boot=2000`). All 4,104 instances were retained (support depths existed for
every instance; no drops).

| stratum | taub(r_val, AUC_common) | Δ vs separation (AUC_common) | 95% CI |
|---|---|---|---|
| tiered nt3 | 0.805 | **-0.019** | [-0.032, -0.007] |
| tiered nt4 | 0.751 | **-0.042** | [-0.056, -0.029] |
| tiered nt2 | 0.825 | +0.087 | [0.058, 0.118] |
| flip (all 6 strata) | 0.66-0.74 | 0 to +0.20 (CI incl. 0 only at cov1.0/bw0.10) | mostly r_val wins or ties |

**Conclusion: unchanged.** Restricting to common-support depths, separation
still beats r_val on tiered nt3 (Δ=-0.019) and nt4 (Δ=-0.042), same sign and
similar magnitude to the whole-sample result in Section 1 (-0.013, -0.031);
tiered nt2 and every flip stratum still favour r_val or tie. Full tables
(incl. Δ vs neg_phi_1/n_k/shd_truth) in `auc_common_support_paired.csv`.
