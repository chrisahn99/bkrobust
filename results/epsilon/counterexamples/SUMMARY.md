# Two naive readings of "bias grows past `r_val`", and why both are false

This directory holds the negative results for the ε-bias radius. Both are
needed: the first says why `B` has to be maximised over a shell rather than read
off one perturbation, and the second says why the incumbent mean-based `r_eps`
cannot be computed the way it currently is.

Everything here is replayable. `replay.py` re-derives each headline number from
the stored witness — the SEM weights and noise variances are persisted, so no
search has to be rerun to check a claim.

---

## C1 — the bias of an individual state is *not* monotone in its distance from `G0`

**Found, minimally, at four nodes** (355 witnesses at `n = 4`, from 761 screened
candidates; the search started at `n = 4` and never needed to go higher).

```
CPDAG      V0->V1  V2->V1  V3->V1   V0-V3   V2-V3
true DAG   V0->V1  V2->V1  V3->V0   V3->V1  V3->V2
G0         V0->V1  V2->V1  V3->V0   V3->V1  V3->V2      (|K_G0| = 2)
query      X = V2,  Y = V1,  Z = {V0, V3}
```

| state | reached by | `d(G0, ·)` | `B` |
|---|---|---|---|
| `G` = `V0->V1 V2->V1 V3->V0 V3->V1 V2-V3` | retracting `V3->V2` | **1** | **1.3405** |
| `H` = `V0->V1 V0->V3 V2->V1 V3->V1 V3->V2` | retracting and re-orienting `V0—V3` | **2** | **0.0000** |

One shell *further out*, the bias is exactly zero. The margin is the whole
1.3405, so this is not a numerical artefact.

**Why this is consistent with Theorem B, and what it actually shows.** Theorem B
is monotonicity along the **order**, not along the **distance**, and `G` and `H`
are incomparable: `H` flips `V3->V0` to `V0->V3`, so it is not above `G` and
nothing relates their biases. Distance counts elementary revisions in either
direction; the order counts only knowledge. The two coincide on the retraction
up-set and nowhere else.

The consequence for the method is exactly the design of `beta_up`: a certificate
cannot be built by walking one perturbation path and watching the bias, because
the bias along such a path is not monotone. It has to be the **maximum over the
whole shell**, which is what Theorem C.3 makes affordable by showing the ball
maximum is attained on the sphere.

---

## C2 — the incumbent `mean_abs_bias` is not monotone under the order, and its radius is wrong because of it

`bkrobust.synth.runner._radius_eps` thresholds
`bkrobust.core.oracle.bias_stats(...).mean_abs_bias`, which redraws a SEM on
every DAG of `[G]` and averages. Retracting knowledge enlarges `[G]`, and an
average over a larger set can be smaller, so the quantity has no reason to be
monotone — and is not.

### C2a — the flagship witness

`[G]` is a **proper** subset of `[H]`, yet the mean bias at `G` is more than
twice that at `H`:

```
G   V0->V1  V0->V2  V0->V3  V2->V1  V3->V2
H   V0->V1  V0-V2   V0-V3   V1-V2   V2-V3
```

| | `mean_abs_bias(G)` | `mean_abs_bias(H)` | deficit | deficit / combined MC s.e. |
|---|---|---|---|---|
| across 8 independent repetitions, 400 draws each | 0.964 – 1.033 | 0.433 – 0.461 | 0.517 – 0.583 | **11.4 – 13.4** |

All 8 repetitions reproduce the sign (reproduction rate 1.00), at a minimum of
**11.4 combined standard errors**. This is a property of the population mean, not
of a draw.

### C2b — how often, after removing Monte-Carlo noise

The first pass counted a violation whenever the deficit exceeded `1e-6` at 80
draws. That threshold sits far below the Monte-Carlo error of an 80-draw mean, so
the resulting rate conflated "the estimand is non-monotone" with "this estimate
of it happened to be". `experiments/run_c2_rate_mc_aware.py` separates them: a
pair counts only if the deficit exceeds **3 combined standard errors in both of
two independent repetitions**, at 300 draws per state. That rule is conservative
and will miss real violations with small margins, so the rate is a **lower
bound**.

Scored on **identical pairs**, over 45 instances and 801 comparable pairs:

| quantity | violations | rate |
|---|---|---|
| `mean_abs_bias`, naive `1e-6` rule | 351 / 801 | 43.8% |
| `mean_abs_bias`, MC-aware rule (**lower bound**) | 278 / 801 | **34.7%** |
| `B` (closed form, tolerance `1e-9`) | **0** / 801 | **0.0%** |

Margins of the MC-aware violations, in combined standard errors: min 3.0,
median **9.6**, max 30.0 — so these are not borderline.

Running both quantities on the same pairs is the point. Two separate sweeps
would leave open that the instances differed.

### C2c — the consequence: the incumbent radius is not merely unjustified, it is wrong

Upward (retraction-only) search is exact for `B` because `{B > ε}` is an up-set
(Theorem U). It is **not** licensed for `mean_abs_bias`, and the lack of a
licence is not theoretical. Comparing the retraction-only radius against a
genuine brute-force BFS over the whole enumerated corrected space, with the same
mean-based predicate:

- 40 instances, **644 (instance, ε) probes**, **269 disagreements** (41.8%).

The disagreements run in the dangerous direction. In the first recorded case the
retraction-only search returns `UNREACHED` — *no perturbation reaches this bias
anywhere* — while brute force finds a violating state at distance 2. A search
that reports infinite robustness where a two-step perturbation already exceeds
the threshold is not a conservative approximation; it is an incorrect answer that
overstates robustness.

---

## What was checked and found *not* to be a problem

On the same 801 pairs, `B` had **zero** monotonicity violations. Had it had any,
that would have refuted Theorem B and would be the most important line in this
file. It is recorded here as a control precisely so that the C2 rates cannot be
read as a property of the harness rather than of the quantity.

## Scope, stated so it is not overread

- Node counts 4 and 5 only, with `|undirected edges| <= 6`, so that the full
  corrected space is enumerable and the "brute force" comparison is genuinely
  exhaustive rather than itself a search.
- `B` is evaluated at one SEM draw per instance; `mean_abs_bias` at 300 draws per
  state per repetition. The asymmetry is not a bias in the comparison: `B` is a
  closed-form maximum at a fixed covariance and has no draw-to-draw variation to
  average out, which is itself one of the differences at issue.
- The C2b rate is a lower bound by construction. The true population rate is
  higher; how much higher is not established here.
