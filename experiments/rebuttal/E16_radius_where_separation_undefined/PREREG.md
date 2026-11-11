# E16 preregistration: radius ranking where separation is undefined (2026-11-06)

Serves reviewer C2Bd, weakness 3 ("separation performs comparably; the radius is defined where s is absent").

**Question.** On the published-network units of Table 1 where the separation proxy s is undefined, does the breakdown radius still rank survival?

**Registered statistic.** Kendall tau_b of the radius (GAC criterion) with the survival AUC (`AUC_frac`) on the pooled real-name units of the submitted cohort that have a radius, an AUC, and no separation value (expected n = 1,033, 14 networks). 95% percentile interval from a network-cluster bootstrap, 10,000 resamples, seed 0, computed by the submitted `real_analyse.tau_for_stratum`.

**Secondary statistics.** (a) The same tau_b on the corrected cohort (submitted + the 94 empty-valid units). (b) Within-state tau_b of the radius on the undefined units of radius-varying (network, condition) states, both cohorts. (c) On the corrected cohort, the within-state paired lead tau_b(radius) - tau_b(separation) on units where separation is defined, same network draws for both, 10,000 resamples. Reference cells that must reproduce published numbers: pooled radius 0.319 (n 2,422), within-state radius 0.419 (1,973 units), pooled separation 0.284 and paired lead +0.039 (1,389 units), within-state paired lead +0.187 (1,213 units, 79 states).

**Radius.** The stored back-door radius is replaced by the GAC radius of `results/radius_gac_recompute/` (joined on shard_id, frame_row_id). Every scorable unit must have an exact GAC radius equal to the stored one; any missing or different value stops the run with a non-zero exit.

**Prediction.** tau_b > 0 with the lower bound of the interval above zero on the 1,033 undefined units; within-state paired lead on the corrected cohort positive with interval above zero (submitted cohort: +0.19 [0.07, 0.90]).

**Outcome against the authors.** If the registered interval contains zero, the rebuttal says only that the radius returns a value where separation does not, and makes no ranking claim on those units. If the corrected-cohort paired lead interval contains zero, the rebuttal drops "leads separation within state" and reports both tau_b side by side.

**Populations, seeds, exclusions.** Real-name naming only; units with missing radius, UNREACHED radius or missing AUC dropped (the Table 1 rule). Network is the bootstrap cluster. Seed 0 everywhere, as in Table 1. No other exclusion; no unit is dropped after seeing results.

**Stopping rule.** One full run with 10,000 resamples. No re-run with other seeds, endpoints or subsets unless a consistency check fails, in which case the defect is fixed and the run repeated in full, with both runs reported.

## Deviations from the plan

- The undefined subset is the same 1,033 units in both cohorts: all 94 added units have separation measured. (a) can therefore differ from the registered value only through the within-state state selection, and is reported anyway.
- No script in the repository produced the published within-state paired lead (+0.19 [0.07, 0.90], `task5_followup.json`). It is reimplemented here with the draws of `llm_paired_v2._within_state_bootstrap`; its point must match to 1e-6 (hard check), its interval is compared at tolerance 0.02 and a mismatch is reported, not fatal.
- Input data are the files committed in the submitted repository, not copies under `_data/`.
