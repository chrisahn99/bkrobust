# Every number in REPORT.md, with the file it comes from

Paths are relative to this directory. Where a number was corrected or retired after first being
reported, the correction is recorded here rather than quietly applied.

| number | meaning | source |
|---|---|---|
| flip = 2 units | our operator reverses an existing skeleton edge, so it is a retraction plus an assertion in the framing's units | `pilot/code/run_linear.py:132` |
| `MAX_K = 4` | global budget of the pilot, i.e. radius 8 in the framing's units | `pilot/code/run_linear.py:34` |
| 0.59 median, up to 128× | damage conditional on `O*` moving, heavy-tailed, scale-free graphs | `pilot/RESULTS.md` |
| 7 688 | knowledge sets on which Meek closure was brute-force validated against its definition | `pilot/RESULTS.md`, `pilot/code/test_machinery.py` |
| 0.937 → 0.508 | var-sortability control, naive parameterisation → iSCM (the Reisach leak, removed) | `pilot/RESULTS.md` |
| 33.0 % | single-orientation errors that are Meek-INCONSISTENT, caught free by the graphical check | `pilot/RESULTS.md`, outcome table under NUMBER 1 |
| 15.3 % / 22.9 % | loud non-amenable failures, as a share of all / of Meek-consistent single errors | same table |
| 9.4 % / 14.1 % | silent bias, as a share of all / of Meek-consistent single errors | same table |
| 0.0 % (0 / 754) | `O*` changed but still valid, i.e. the efficiency-only band, measured EMPTY | same table |
| 0.037 / 0.070 / 0.016 / 0.330 | censoring of `rho*_any` at hop distance 0, by ensemble | `e1prime/RESULTS.md` (primary), cross-checked in `locality/x1/RESULTS.md` |
| 1.000 / 0.935 / 0.955 / 0.977 | censoring at hop distance 1, by ensemble | `e1prime/RESULTS.md` |
| 1.000 (n = 331) | censoring at hop distance >= 2, `large` | `e1prime/RESULTS.md` |
| 13-27× | the separation between distance 0 and distance >= 1 | derived in `e1prime/RESULTS.md` |
| 0.9345 | the same censoring at n = 2 000, 20 000 and infinity, to four decimals | `e1prime/RESULTS.md` |
| 60.3 / 88.7 / 96.3 % | share of SCMs with at least one statement at distance 0, `large` / `licensed` / `original` | `e1prime/RESULTS.md`, `locality/VERDICT.md` |
| hop anchor | `hop_dist_from(skeleton(C), [x, y])`, minimum over both endpoints of the statement | `locality/x1/code/x1_ops.py:271` |
| 69.4 / 30.6 / 12.9 % | b-LOAD pool composition on `licensed` (5 555 / 2 445 / 1 030 of 8 000). `original` 73.1 / 27.0, `large` 91.4 / 8.6. ⚠️ ORACLE labels: not a triage rule | `locality/x1/results/x1_analysis_licensed.json` → `P3_A26.composition_Sloc_pool` |
| 32.4 % vs 33.0 % | E1-prime's Meek-inconsistency rate against the pilot's, two independent runs agreeing to half a point | `e1prime/RESULTS.md`, `pilot/RESULTS.md` |
| 0.554 → 0.515 | censoring of `rho*_any`, n = 200 → infinity, against a hard floor | `e1prime/RESULTS.md` |
| 64.6 % | licensed problems with `d(1) = 0` exactly, at every sample size | `e1prime/RESULTS.md` |
| 0.166 and 69.0 % | composite `min(rho*_se, rho*_ident)`: censoring, and mass at `rho* = 1` | `e1prime/RESULTS.md` |
| 32.4 % | Meek-inconsistency rate, from the G1 denominator correction | `e1prime/RESULTS.md`, `e1prime/code/analyse_arm1.py` |
| 23.6 % / 38.6 % | single / double errors caught free by Meek in ARM 2 | `e1prime/RESULTS.md` |
| 0.9509 | naive coverage at r = 0, n = 20 000, against 0.95 nominal; zero regeneration-gate failures | `e1prime/RESULTS.md`, `e1prime/code/analyse_arm2.py` |
| 0.951 / 0.720 / 0.568 | naive CI coverage at r = 0 / 1 / 2, n = 20 000 | `e1prime/RESULTS.md` |
| 1.61 → 6.84 | mean robust/naive width ratio, n = 200 → 20 000 | `e1prime/RESULTS.md` |
| 2 500 × 20 | SCMs × datasets in ARM 2 | `e1prime/PREREG.md` |
| complete at `p <= 5` | the enumeration in X2 is exhaustive there; `p = 6` is a sample | `locality/x2/RESULTS.md`, `locality/x2/code/run_lemma.py` |

## Corrections we applied to ourselves

| what | first reported | corrected | why |
|---|---|---|---|
| `0/791` silent bias at distance >= 1 | pilot, 17 July | **retired** | no provenance found in the audit of 19 August. Use the censoring table instead |
| silent-bias rate at `rho = 1` | 0.141 | 0.161 on the corrected denominator | the first cut divided by all `rho = 1` members; 0.141 uses consistent members only. Both denominators are now reported |
| width ratio | median 1.00 | mean 6.84× reported beside it | the ball is inert for 69 % of problems, so the median is 1.00 by arithmetic. The pre-registered statistic is honoured and the mean is printed next to it |
| G2 constant | "within 0.03 of 0.336" | comparator is 0.271 | 0.336 is the `original` grid; the `k8` ensemble uses the `licensed` grid |

## What has never been run

- Anything with an **estimated** CPDAG. Every run here uses the oracle CPDAG, which is the regime all
  four b-LOAD reviewers called too narrow. The licensed phrasing is "computable from
  `(C, K, Sigma-hat, n)`", never "robust under an estimated CPDAG".
- Any **nonlinear or non-Gaussian** mechanism beyond the 60 nonlinear SCMs in the pilot.
- The **skeleton perturbation** restricted to borderline or subsampling-unstable adjacencies, which is
  the natural repair for the oracle-CI assumption in section 2 of the framing.
