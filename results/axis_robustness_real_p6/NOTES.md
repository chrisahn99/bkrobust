# NOTES — reviewer Point 6 on the real corpus (TRACK D)

Scope: two reviewer gaps against the committed, verified real-structure survival sweep
(`results/axis_robustness_real/`, 831-row frame, 725/725 completed shards, `make
re11-verify` 26/26). **Nothing under `results/axis_robustness_real/` is modified or
re-hashed.** Everything here is additional, read-only against that tree.

- `src/bkrobust/robustness/real_baselines_p6.py` — computes `phi_1` and `r_claim` on
  the frozen 831-row frame (`separation` and `k_g0` already existed and are carried
  through, not reimplemented), and an extended tau table ranking seven predictors
  instead of the committed four.
- `src/bkrobust/robustness/paired_resample.py` — a generic paired-cluster-bootstrap
  module (predictor columns + cluster column + endpoint column in, delta + CI +
  P(delta>0) + leave-one-cluster-out out), applied here to `radius` vs. each of the
  six baselines.

Run order: `real_baselines_p6.py` first (writes `baselines_per_instance.csv`, needed
by the second script), then `paired_resample.py`. Both took under 12 minutes total on
this machine (`baselines_manifest.json`: 666.7s; `paired_diff_manifest.json`: 190.8s).

---

## 1. Which baselines are degenerate by construction

### `r_claim` — completely degenerate, everywhere, not just at coverage 1.0

**`r_claim == 1` on all 831 of 831 frame rows.** Every single admitted row's `Z*` is
invalidated by retracting exactly one of the analyst's asserted claims. The task brief
predicted this *at coverage 1.0* (543 rows), because `benchmarks.measure.fast_gate`'s
admission test literally is that check there (its `k_true` equals the row's own `K`
exactly when coverage is 1.0, and its `o` equals the row's own `Z*` there too). That
prediction is confirmed empirically (`baselines_manifest.json:degeneracy_check_coverage_1`:
`543/543` resolved at depth 1). **What was not predicted, and is the actual finding, is
that the same thing holds at coverage 0.5 (182/182) and coverage 0.25 (106/106) as well**
— every admissible pair on this corpus, at every coverage level tested, has *some*
single claim whose retraction breaks `Z*`. `r_claim` therefore carries **zero ranking
information** anywhere in this corpus: every `tau_b` and every paired comparison against
it reports `"undefined (predictor constant)"` (`analysis_tau_extended.csv`,
`paired_diff_real.csv`). This is reported as the finding it is, not deepened into a
combinatorial search that could not change a constant into a variable.

### `phi_1` — not constant, but sign-degenerate by the same construction

`phi_1` (the *fraction* of `K` whose individual retraction breaks `Z*`) is never zero
on this corpus, for the identical reason `r_claim` is never above 1: the admission gate
guarantees at least one such claim exists. But **`phi_1`'s magnitude does vary** —
0.013–1.0 at coverage 1.0 (mean 0.343), 0.2–1.0 at coverage 0.5 (mean 0.452), 0.5–1.0
at coverage 0.25 (mean 0.717) — so, unlike `r_claim`, it is not useless as a ranking
predictor. Its mean rising as coverage falls is expected (fewer claims left to
retract, so each one matters more) and is not itself informative about the radius axis.

### `k_g0` and `separation` — already known, carried through, not new findings here

`k_g0` was already one of the four ranked predictors in the committed `analysis_tau.csv`
and is carried into `analysis_tau_extended.csv` unchanged (imported from
`real_analyse.tau_for_stratum`, not recomputed). `separation` was already a column on
every analysis unit but had never been *ranked* — only used to stratify. It is not
degenerate in the same sense as `r_claim`/`phi_1`, but it is **structurally undefined on
44% of the corpus** (368 of 831 frame rows carry `separation_status =
"no_z_member_in_component"`, i.e. no member of `Z*` lies in `X`'s undirected component,
so "distance to `Z*`" has no referent). Every tau and paired computation here drops
those rows — never imputes zero or any other number — and the drop count is on every row
of `analysis_tau_extended.csv` (`n_excluded_nan`) and `paired_diff_real.csv`
(`n_excluded`). In the primary flip-bw0 strata this removes 44–56% of the units before
any correlation is computed (`n_excluded_nan` at `flip_cov100_bw000`: 237 of 543;
`flip_cov050_bw000`: 85 of 182; `flip_cov025_bw000`: 48 of 106).

---

## 2. Scope decision: where `phi_1`/`r_claim` can honestly be compared

`phi_1` and `r_claim` are computed from the frame's own `K` (`select_knowledge(dag,
cpdag, coverage)`, the analyst's **uncorrupted** claim set) and the frame's own `Z*`.
That is exactly the flip arm's `base_wrongness = 0.00` construction — and *only* that
one. Every other unit in the committed sweep — flip at `base_wrongness` 0.10/0.25/bw_abs=1,
and every tiered unit (`K_ref` built from temporal tiers, not `select_knowledge`) — has a
**different, corrupted or independently-generated** operative `K`, which the frozen
frame does not capture and which recomputing was out of scope for "the frozen 831-row
frame." Attaching the frame's `phi_1`/`r_claim` value to those units would silently
compare a stale predictor to an endpoint it was never computed against.

So `real_baselines_p6.build_extended_units` attaches `phi_1`/`r_claim` **only** to
flip-arm units with `base_wrongness == 0.0` — three strata: `flip_cov100_bw000`,
`flip_cov050_bw000`, `flip_cov025_bw000`. In every other stratum (9 more flip strata,
3 tiered strata) both columns are `None` for every unit, which correctly surfaces as
`status = "insufficient_n"` (paired) or `"insufficient_n"`/dropped-to-zero (tau) rather
than a fabricated result. This is stated here once rather than annotated on all 12
affected rows of the table. `separation`, `k_g0`, `shd_truth`, `n_k` carry no such
restriction — they are already correctly computed per-arm, per-unit by the existing
pipeline, so they are ranked across all 15 strata.

---

## 3. Paired `r_val` (radius) vs. each baseline

All bootstraps: 10,000 resamples, seed 0, cluster = network, percentile CI at
2.5/97.5%, matched sample (a row missing *either* predictor is excluded from **both**
tau computations in every resample — never an independent per-predictor drop, which
would let the pair be scored on different populations). Endpoint is the `*_usable`
variant per arm (`AUC_frac_usable` flip, `AUC_rate_usable` tiered), per Appendix E.
Full table: `paired_diff_real.csv` (90 rows = 15 strata x 6 baselines); per-network LOO
detail: `paired_diff_loo_detail.json`.

### r_val vs. `r_claim` — no comparison possible

`r_claim` is constant everywhere (Sec. 1), so every one of the 15 rows reports
`"undefined (predictor or endpoint constant on matched sample)"`. This is not a null
result about radius; it is the absence of a baseline to compare against.

### r_val vs. `phi_1` — radius wins, robustly, everywhere it can be tested

The only three strata where the comparison is valid (Sec. 2) all show radius
significantly and robustly outranking `phi_1`:

| stratum | delta (tau_radius − tau_phi1) | 95% CI | P(delta>0) | LOO range | LOO flips? |
|---|---|---|---|---|---|
| `flip_cov100_bw000` | **+0.491** | [0.168, 0.858] | 0.9996 | [0.332, 0.583] | No |
| `flip_cov050_bw000` | **+0.392** | [0.120, 0.904] | 0.997 | [0.284, 0.573] | No |
| `flip_cov025_bw000` | **+0.298** | [0.035, 0.365] | 0.977 | [0.176, 0.322] | No |

All three CIs exclude zero, all three LOO ranges stay entirely positive (no single
network's removal reverses the sign), and the effect is consistent in direction and
roughly monotonic in coverage. This is the strongest, cleanest positive result in this
track's output: **radius outranks the single-claim-fragility baseline at every coverage
level where the comparison is well-posed.**

### r_val vs. `shd_truth` — one robust win, otherwise not significant

Ten of eleven applicable strata (`shd_truth` is defined for both arms) show a CI
spanning zero, several also flagged `loo_verdict_flips = True` (the point estimate
itself is not stable under leave-one-network-out). The exception:

| stratum | delta | 95% CI | P(delta>0) | LOO range | LOO flips? |
|---|---|---|---|---|---|
| `flip_cov100_bwa1` (supplementary, |K|·base_wrongness=1 absolute-error variant) | +0.279 | [0.075, 0.483] | 0.994 | [0.227, 0.346] | No |

One genuine, LOO-stable win, in a supplementary (not primary) stratum. Everywhere else
— including all three primary bw0 flip strata and all three tiered strata — there is
**no significant paired difference** between radius and `shd_truth`. (`flip_cov100_bw000`
itself reports `shd_truth` as constant, exactly as `PREREGISTRATION.md` Sec. 2.1
predicted for the `b=0`/coverage=1.0 stratum, so no paired delta exists there at all.)

### r_val vs. `n_k` — radius **loses**, robustly, at coverage 0.25

This is the null/adverse result the task brief asked to be reported at the same weight
as a positive one, and it is real:

| stratum | delta | 95% CI | P(delta>0) | LOO range | LOO flips? |
|---|---|---|---|---|---|
| `flip_cov025_bw000` | **−1.351** | [−1.903, −0.557] | 0.000 | [−1.708, −0.947] | No |
| `flip_cov025_bw010` | **−1.306** | [−1.895, −0.507] | 0.000 | [−1.677, −0.894] | No |
| `flip_cov025_bw025` | **−1.323** | [−1.895, −0.525] | 0.000 | [−1.698, −0.915] | No |

Replicated across all three base-wrongness levels at coverage 0.25, CI entirely
negative in each, LOO range entirely negative in each (`n_k` beats radius even with any
one network dropped). At coverage 0.25, `n_k` (`tau_b = +0.824`, `analysis_tau_extended.csv`)
correlates strongly and positively with survival while `radius` (`tau_b = −0.527`)
correlates negatively, on a thin sample (8 networks, 106 admissible pairs, `flip_cov025_bw000`)
— exactly the corpus's known thin-strata risk (`PREREGISTRATION.md` Sec. 6.2's
leave-one-network-out warning). **This is not softened, re-stratified, or explained
away here: on this corpus, at the sparsest coverage level, the naive count of asserted
claims outranks the breakdown radius, and the finding survives removing any single
network.** At coverage 0.5 and 1.0, and on the tiered arm, no significant paired
difference against `n_k` is detected in either direction (all CIs span zero).

### r_val vs. `k_g0` — no significant paired difference detected, anywhere

All 15 strata have a CI spanning zero. About half are additionally flagged
`loo_verdict_flips = True`, meaning even the (non-significant) point estimate changes
sign or crosses zero under some single-network removal. No claim of superiority or
inferiority is supportable here in either direction.

### r_val vs. `separation` — no reliable significant difference

Fourteen of fifteen strata have a CI spanning zero. The apparent exception,
`flip_cov025_bwa1` (delta = −0.023, CI = [−0.023, −0.017]), **rests on `n_clusters = 2`**
— only two networks contribute to that matched subsample. A percentile bootstrap over
two clusters with replacement has at most three distinct resample outcomes, so this
"CI" is not a meaningful interval and the apparent significance is an artifact of too
few clusters, not evidence. **This is flagged explicitly and excluded from any claim of
a real finding**; nowhere else does `separation` show a significant paired difference.

---

## 4. Summary, stated plainly

- `r_claim` is uninformative everywhere on this corpus (constant at 1) — a genuine,
  stronger-than-predicted degeneracy, not limited to coverage 1.0.
- `phi_1` is never zero (construction artifact of the admission gate) but varies enough
  in magnitude to be a real ranking predictor, and **radius robustly, significantly
  outranks it at every coverage level where the comparison is valid** (the three flip
  bw=0.00 strata).
- Radius has **one** robust win against `shd_truth` (a supplementary stratum) and **no**
  significant difference against it in every primary stratum.
- Radius **loses**, robustly and reproducibly across three base-wrongness replicates, to
  the naive `n_k` baseline at coverage 0.25 — reported here at full weight, not
  rescued by re-stratifying or swapping the endpoint.
- Radius shows **no** significant paired difference against `k_g0` or `separation`
  anywhere reliable (the one nominal exception for `separation` is an n_clusters=2
  artifact, named as such).
- Roughly a third of all 90 paired cells (33/90) are flagged `loo_verdict_flips = True`
  — a reminder, consistent with `PREREGISTRATION.md` Appendix H, that most of these
  strata are thin enough for a single network to change the verdict, and any of these
  numbers quoted alone without its LOO column overstates its own certainty.

## 5. Reproduction

```
PYTHONPATH=src /usr/bin/python3 src/bkrobust/robustness/real_baselines_p6.py
PYTHONPATH=src /usr/bin/python3 src/bkrobust/robustness/paired_resample.py
```

Both scripts refuse to run if `results/axis_robustness_real/frame.jsonl`'s SHA-256
does not match `frame_hash.json` (imported check, `real_analyse.verify_frame`).
`real_baselines_p6.py` additionally asserts that its recomputed `select_knowledge(dag,
cpdag, coverage)` reproduces the frame's own `n_k` for every `(network, coverage)`
cell before trusting anything downstream of it — a mismatch is a crash, not a silent
row, matching this corpus's existing convention for hard invariants.

`paired_resample.py`'s core API (`paired_cluster_bootstrap`, `leave_one_cluster_out_delta`,
`paired_comparison`) takes plain `list[dict]` rows and column names; it does not import
`real_analyse` or read anything under `results/axis_robustness_real/` at module scope,
so the same functions apply unchanged to the synthetic N=1000 dataset another track is
producing at `results/axis_robustness_p6/`, whenever that directory exists.
