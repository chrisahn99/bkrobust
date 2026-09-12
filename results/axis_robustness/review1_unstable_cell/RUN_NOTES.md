# Re-run of the second unstable stratum at N = 1000

`flip, coverage = 1.0, base_wrongness = 0.00`. The published anchor table carried an
"unstable at N = 200" marker on two cells. Only the null cell (coverage 0.5, base
wrongness 0.25) had been re-run at N = 1000; this run retires the other.

Design, seeds, endpoint and bootstrap are those of the committed null-cell run
(`bkrobust.robustness.run_nullcell.run_full_mode`); only the stratum constants differ.
Driver: `experiments/review1_unstable_stratum_n1000.py`.

| check | endpoint | N | n | tau_b | 95% CI |
|---|---|---|---:|---:|---|
| reproduction | `AUC_frac`, first 200 of 1000 draws | 200 | 240 | +0.519 | [+0.442, +0.591] |
| full | `AUC_frac` | 1000 | 240 | **+0.574** | [+0.502, +0.641] |
| restricted | `AUC_frac_usable`, `n_eval >= 30` | 1000 | 240 | **+0.537** | [+0.484, +0.585] |

**Guard.** The reproduction check returns +0.519, matching the committed unfiltered
N = 200 value for this stratum exactly, so the instance population and the endpoint are
the published ones.

**Verdict.** The published restricted-endpoint value for this cell was **-0.091**, with
an interval covering zero, on 229 of 240 instances (11 had no usable depth at N = 200).
At N = 1000 every instance has a usable depth and the same endpoint gives **+0.537**,
with the interval clear of zero. The negative sign was an artefact of the
`n_eval >= 30` filter at low draw count, not a property of the stratum: the filter
selects different depths for high- and low-radius instances, and that selection
disappears once the raw counts are five times larger. The cell is **resolved positive**
and the "do not quote alone" marker no longer applies to it.

The other marked cell (coverage 0.5, base wrongness 0.25) remains a confirmed null at
N = 1000 (tau = +0.032).
