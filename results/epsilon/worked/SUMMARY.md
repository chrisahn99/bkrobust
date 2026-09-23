# The epsilon-bias radius on the running example

Companion to Appendix C of the paper. The bias column there is a mean over
drawn coefficients; here it is the exact worst case `beta_up(d)` at one fixed
documented SEM draw (seed 20260916), with the distribution over
400 further draws reported below it.

Two things to read carefully.

**Shells past `|K_G0|` carry a dash, and that is not a gap.** The staircase is
defined on the retraction up-set of `G0`, which by Proposition R reaches depth
`|K_G0|` and no further. States beyond that depth are reached by *orienting*
edges the analyst declined to orient, not by retracting claims they made, so no
error budget denominated in their own claims can reach them: an analyst with
`|K_G0|` orientations cannot have more than `|K_G0|` of them wrong. They are
therefore outside Theorem D by construction. The radii are nonetheless
cross-checked against a brute-force BFS over the *whole* space, which does visit
them, and agree -- which is Theorem C.5 holding in practice.

**The worst case is a step here, not a ramp.** Appendix C's mean rises shell by
shell; the exact worst case does not. At the headline draw it jumps from zero to
`B(Chat)` -- the global maximum over the entire space -- at exactly `r_val` and
stays there. The Appendix C ramp is an averaging artefact: what rises with depth
is the *proportion* of states at which `Z` fails (the `Z invalid` column), not how
badly the worst one fails. Over the draw distribution below, roughly half the
draws saturate immediately and roughly half rise once, so neither shape is
universal; but on this instance the practitioner-facing consequence is that all of
the bias risk arrives at `r_val`.

**`r_eps = r_val` is not a null result.** Even when the epsilon-radius buys no
extra certified shell, it supplies what `r_val` cannot: a magnitude. `r_val` says
the third revision breaks the set; the band says what breaking it costs.

## Scenario A

- `Z` = ['Age', 'Smoke'], `r_val` = 3 (via local_up_fast), `|K_G0|` = 4
- reported estimate `theta_Z` = 0.326211, true effect = 0.326211, realised relative error = 0.0
- up-set radii vs brute-force-over-the-whole-space radii agree: **True**

| shell | states | Z invalid | worst-case bias | relative |
|---|---|---|---|---|
| 0 | 1 | 0 | 0.0000 | 0.0% |
| 1 | 3 | 0 | 0.0000 | 0.0% |
| 2 | 4 | 0 | 0.0000 | 0.0% |
| 3 | 7 | 1 | 0.1629 | 49.9% |
| 4 | 8 | 3 | 0.1629 | 49.9% |
| 5 | 9 | 5 | - | - |
| 6 | 5 | 4 | - | - |
| 7 | 6 | 5 | - | - |
| 8 | 3 | 3 | - | - |
| 9 | 2 | 2 | - | - |

`r_eps` (relative thresholds): {'0.01': 3, '0.05': 3, '0.1': 3, '0.25': 3, '0.5': -1, '1': -1}

## Scenario B

- `Z` = ['Age', 'Smoke'], `r_val` = 3 (via local_up_fast), `|K_G0|` = 5
- reported estimate `theta_Z` = 0.326211, true effect = 0.326211, realised relative error = 0.0
- up-set radii vs brute-force-over-the-whole-space radii agree: **True**

| shell | states | Z invalid | worst-case bias | relative |
|---|---|---|---|---|
| 0 | 1 | 0 | 0.0000 | 0.0% |
| 1 | 2 | 0 | 0.0000 | 0.0% |
| 2 | 4 | 0 | 0.0000 | 0.0% |
| 3 | 5 | 1 | 0.1629 | 49.9% |
| 4 | 8 | 3 | 0.1629 | 49.9% |
| 5 | 8 | 4 | 0.1629 | 49.9% |
| 6 | 8 | 5 | - | - |
| 7 | 4 | 3 | - | - |
| 8 | 5 | 4 | - | - |
| 9 | 2 | 2 | - | - |
| 10 | 1 | 1 | - | - |

`r_eps` (relative thresholds): {'0.01': 3, '0.05': 3, '0.1': 3, '0.25': 3, '0.5': -1, '1': -1}

## Scenario C

- `Z` = ['Age', 'Smoke'], `r_val` = 2 (via local_up_fast), `|K_G0|` = 4
- reported estimate `theta_Z` = 0.326211, true effect = 0.326211, realised relative error = 0.0
- up-set radii vs brute-force-over-the-whole-space radii agree: **True**

| shell | states | Z invalid | worst-case bias | relative |
|---|---|---|---|---|
| 0 | 1 | 0 | 0.0000 | 0.0% |
| 1 | 4 | 0 | 0.0000 | 0.0% |
| 2 | 6 | 1 | 0.1629 | 49.9% |
| 3 | 9 | 3 | 0.1629 | 49.9% |
| 4 | 8 | 4 | 0.1629 | 49.9% |
| 5 | 8 | 5 | - | - |
| 6 | 4 | 3 | - | - |
| 7 | 5 | 4 | - | - |
| 8 | 2 | 2 | - | - |
| 9 | 1 | 1 | - | - |

`r_eps` (relative thresholds): {'0.01': 2, '0.05': 2, '0.1': 2, '0.25': 2, '0.5': -1, '1': -1}

## Parameter dependence

### Scenario A (`r_val` = 3, fixed by the graph alone)

- staircase shape over 400 draws: {'flat': 221, 'step': 179}
- relative rise past `r_val`: {'mean': 0.062241, 'median': 0.0, 'min': 0.0, 'max': 3.086348}
- `r_eps` distribution: {'0.01': {'-1': 6, '3': 378, '4': 16}, '0.05': {'-1': 92, '3': 268, '4': 40}, '0.1': {'-1': 183, '3': 179, '4': 38}, '0.25': {'-1': 290, '3': 93, '4': 17}, '0.5': {'-1': 334, '3': 56, '4': 10}, '1': {'-1': 370, '3': 26, '4': 4}}

### Scenario B (`r_val` = 3, fixed by the graph alone)

- staircase shape over 400 draws: {'flat': 221, 'step': 179}
- relative rise past `r_val`: {'mean': 0.062241, 'median': 0.0, 'min': 0.0, 'max': 3.086348}
- `r_eps` distribution: {'0.01': {'-1': 6, '3': 378, '5': 16}, '0.05': {'-1': 92, '3': 268, '5': 40}, '0.1': {'-1': 183, '3': 179, '5': 38}, '0.25': {'-1': 290, '3': 93, '5': 17}, '0.5': {'-1': 334, '3': 56, '5': 10}, '1': {'-1': 370, '3': 26, '5': 4}}

### Scenario C (`r_val` = 2, fixed by the graph alone)

- staircase shape over 400 draws: {'flat': 221, 'step': 179}
- relative rise past `r_val`: {'mean': 0.062241, 'median': 0.0, 'min': 0.0, 'max': 3.086348}
- `r_eps` distribution: {'0.01': {'-1': 6, '2': 378, '4': 16}, '0.05': {'-1': 92, '2': 268, '4': 40}, '0.1': {'-1': 183, '2': 179, '4': 38}, '0.25': {'-1': 290, '2': 93, '4': 17}, '0.5': {'-1': 334, '2': 56, '4': 10}, '1': {'-1': 370, '2': 26, '4': 4}}
