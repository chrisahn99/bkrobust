# E1 PRE-REGISTRATION — rho* instrument rebuild (rel50 @ |K|=8 vs sign @ |K|<=4)

Written 2026-08-14 BEFORE any |K|=8 simulation was run. Baseline numbers below were
already public in the campaign record (RESULTS.md / summary.json); the |K|=8 outcome
was not observed at the time of writing.

## Quantity under test

The 4-bin distribution of the breakdown radius
rho* in {1, 2, 3, >3}, generic arm, ball radius capped at rho <= 3,
under `crit='rel50'` (|est - tau| > 0.5|tau|), as computed by
`breakdown_radius()` at code/analyse.py:42.

`rho* = 4` in the code means "never overturned inside the ball" = CENSORED.

## Baseline (already on record, |K| <= 4, 600 SCMs)

| crit  | rho*=1 | never (censored) |
|-------|--------|------------------|
| sign  | 0.057  | **0.903**        |
| rel50 | 0.168  | **0.755**        |
| ident | 0.509  | 0.306            |

The verifier's objection is that under `sign` at |K|<=4 the diagnostic returns
"no breakdown" for 90% of problems. The campaign's declared fix is rel50 at |K|=8,
claimed at 26/14/4/56% across rho*=1/2/3/>4 with mean interval width 1.0|tau|.

## Pre-registered decision rule

Let `c8` = censored fraction (rho* > 3) on the |K|=8 arm under rel50, and
`maxbin8` = the largest of the four bin masses.

- **REFUTED** if `c8 >= 0.755`.
  Grounding: 0.755 is the measured rel50 censoring at |K|<=4 on data already on disk,
  and is the threshold in the vault's own registered falsifier for this claim. If
  raising |K| from 4 to 8 does not beat the |K|<=4 rel50 baseline, the declared fix
  buys nothing and the instrument-degeneracy objection stands unanswered.

- **SUPPORTED** if `c8 <= 0.60` AND `maxbin8 <= 0.60` AND at least 3 of the 4 bins
  each hold >= 0.05 of the mass.
  Grounding: 0.60 is the campaign's own claimed |K|=8 figure (56% censored) rounded
  to a bar it must clear. The >=5%-in-3-bins clause is the instrument-quality
  requirement: a radius that is effectively 2-valued is not a radius.

- **INCONCLUSIVE** otherwise (0.60 < c8 < 0.755, or c8 <= 0.60 with a degenerate
  bin profile).

## Mandatory controls

1. **Matched |K|=4 control on the SAME SCMs.** MAX_K (run_linear.py:34) sits upstream
   of a conditional `rng.choice`, so simply re-running with a different MAX_K yields a
   different ensemble. Control here is stronger than seed-matching: one uniform random
   ordering of the orientable statements is drawn per SCM and the |K|=4 arm is its
   first 4 elements, so K4 is a subset of K8 on an identical graph, identical (x,y),
   identical Sigma, identical tau. K4 remains a uniform random 4-subset, i.e.
   distributionally identical to the original protocol.

2. **C13 compliance** (vault registry rule: a null must show the instrument was capable
   of producing a non-null). Demonstrated on data in hand: the same
   `breakdown_radius()` over the same 600 SCMs returns 30.6% censoring under
   `crit='ident'` versus 90.3% under `crit='sign'`. The instrument has a measured
   dynamic range of 60 pp across criteria, so a high censoring rate at |K|=8 would be
   a property of the estimand, not of a dead switch.

3. **Scope must be reported, not assumed.** The |K|=8 arm requires CPDAGs with >= 8
   undirected edges. Whether that is reachable inside the claim's licensed scope
   (ER, p <= 10, expected degree <= 6) is measured and reported as part of the result,
   not silently worked around.

## Secondary (reported, not decisive)

- Mean and median ignorance-interval width / |tau| at rho = 1, 2, 3, both arms.
- Fraction of SCMs with rho* > 0 (trivially 1.0 by construction; reported as the
  fraction whose interval at rho<=3 is non-degenerate).
- Head-to-head table of all three criteria (sign / rel50 / ident) x both |K| arms.
