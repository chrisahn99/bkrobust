# E19 deviations from PREREG.md (written 2026-11-06 after the first full run, before the rerun)

## 1. Threshold ties at the 100 % tolerance

**What happened.** The first full run (6 Nov 18:54, 900 instances, 34 s) stopped with exit 3: the naive arm
differed from the stored study on 16 of 900 instances. On every one of them a relative bias equals the 100 %
tolerance to within a few units in the last place (for example 1.0000000000000018 against 1.0000000000000005).
The library counts a tolerance as reached when the bias is strictly above it, so the radius at 100 %, the length
of the stored staircase and its shape flip with the rounding of the last digit. Nothing else differed: graph,
query, knowledge, adjustment set, validity radius and r_0 matched on all 900. The first run's output is kept in
`results/rebuttal/E19_failed_run1/`.

**Change (decided by the authors before the analysis was read).**
- A relative bias within 1e-9 (relative) of a tolerance is a tie.
- Every arm now also reports, for every tolerance, the radius under both conventions read off the same
  staircase: `r_eps_ge` (bias counted as exceeding when it is at least the tolerance minus 1e-9) and
  `r_eps_gt` (when it is above the tolerance plus 1e-9), and the staircase shape under each convention.
  To read both, the staircase is traversed one step past the largest tolerance (an extra tolerance of
  1 + 1e-6 is passed to `certify`). The stored staircase, `r_eps`, the shape and the certificate audit are
  still computed on the staircase truncated by the library's own rule, so the registered statistic keeps
  the study's convention.
- The reproduction check no longer counts a difference in `r_eps`, shape or staircase length as a failure
  when the naive or stored staircase holds such a tie, the differing tolerances are all tied, and nothing
  else differs. Such rows are listed in `tie_only_differences` and counted. Every other check is unchanged.

**What is unchanged.** Units, arms, seeds, coefficient draws, criterion, registered statistic, thresholds of the
outcome against the authors, exclusions and stopping rule.

**How to undo.** Restore `run.py` and `analyse.py` from the first run's copies (checksums in
`results/rebuttal/E19_failed_run1/run_config.shard0of1.json` name the commit; the pre-change files are in the
laptop scratch copy `run.py.orig`, `analyse.py.orig`) and read `E19_failed_run1/`.
