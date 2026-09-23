# Matched Cross-Arm Corruption Experiment — Run Notes

2026-09-10. Fixes the two defects that made session 7 retract its cross-arm
detectability claim (PREREGISTRATION.md Appendix E.3): different units
(fraction of claims reversed vs. fraction of nodes relocated) and different
instance populations (`tiered` is generative, `flip` corrupts an existing
`K`). This run builds every instance **once**, the tiered arm's own way
(`K_ref = tiered(dag, cpdag, rng, n_tiers, 0.0)`, `G0`, `Z*` read off it),
then corrupts that *same* `(cpdag, K_ref, G0, Z*)` **both ways**, scored on a
shared, unified intensity axis: `|dir(G0) Δ dir(G)| / |dir(G0)|`.

Code: `src/bkrobust/robustness/crossarm.py` (corruption/scoring logic),
`src/bkrobust/robustness/run_crossarm.py` (CLI/orchestration/I/O).

---

## VERDICT — matched-instance, matched-intensity comparison

> **Do correlated (tiered) errors degrade GAC validity faster than uniform
> (flip) errors at identical corruption intensity?**

**Yes, in the intensity range where both arms have real data (≈0.15–0.55),
and the answer is *not* a null.** Restricting to (instance, bin) cells where
*both* arms have ≥30 non-contradictory samples for that instance (per the
brief), and running a per-bin sign test plus paired bootstrap over matched
instances (`xarm_paired_test.csv`):

| bin | n instances matched | mean(S_tiered − S_flip) | sign-test p |
|---|---|---|---|
| 0.15 | 19 | **+0.034** | 0.77 (n.s.) |
| 0.20 | 47 | **−0.021** | 1.7e-4 |
| 0.25 | 35 | **−0.016** | 0.043 |
| 0.30 | 44 | **−0.072** | 9.0e-6 |
| 0.35 | 68 | **−0.053** | 3.3e-6 |
| 0.40 | 18 | **−0.088** | 0.031 |
| 0.55 | 16 | **−0.213** | 0.0042 |

(Full table, including thinly-supported bins with n_instances_matched < 15,
in `xarm_paired_test.csv`; those are reported but not leaned on.)

Six of the seven bins with reasonable instance support (n ≥ 15) show
`S_tiered < S_flip` — correlated corruption survives *less* often than
uniform corruption at the same measured intensity — and five of those six
are individually significant by sign test (α = 0.05). Only bin 0.15 (n=19)
goes the other way, and is not significant. This is a **consistent,
mostly-significant, non-null result in one direction**: at matched
intensity, tiered errors are *more* dangerous to `Z*`'s validity than flip
errors, not less.

This is compounded by the contradiction-rate side of the question (pooled
across all instances *within each shared intensity bin*, from
`xarm_bins.csv`, controlling for intensity — not restricted to the
paired-test's ≥30 subset): tiered's contradiction rate is **higher** than
flip's at every shared intensity bin in [0.10, 0.85), often by 20–30 points
(e.g. bin 0.2: tiered 0.754 vs. flip 0.447; bin 0.3: tiered 0.715 vs. flip
0.439; bin 0.5: tiered 0.786 vs. flip 0.652). Session 7 found correlated
errors "appeared less self-revealing" but explicitly could not separate
that from the population confound (Appendix E.3). This matched design can,
and it says the opposite: **at matched intensity, correlated errors are
actually more self-revealing (more likely to hit a logical contradiction),
not less** — and, among the errors that don't self-reveal, still more
likely to silently invalidate `Z*`. Both halves of "tiered looks safer"
from session 7 do not survive the matched design.

(Caution on aggregation: this bin-controlled comparison is the one that
matters and is unambiguous. An *un*-controlled, whole-arm pooled
contradiction rate is genuinely ambiguous depending on how it's weighted —
per-instance-averaged it's tiered 0.633 vs. flip 0.524 [flip looks lower],
but sample-weighted-pooled it's tiered 0.633 vs. flip 0.670 [flip looks
higher] — because flip's grid density (`len(K_ref)` points) varies by
instance while tiered's is fixed at 10 points/instance, so pooling and
per-instance-averaging weight instances differently. Both unconditional
numbers are reported in J6 for completeness, but neither is intensity-
matched, and they disagree with each other, which is itself evidence that
an unconditional aggregate is not a safe substitute for the binned
comparison above.)

Reported as found, not as expected: this is *not* the null the brief
explicitly allowed for. It is a real effect in a design built specifically
to be capable of returning one.

**Caveat, stated plainly (J5):** this verdict covers intensity ≈0.10–0.85.
Below ≈0.10 only tiered has (heavy) coverage — flip's minimum achievable
intensity is bounded below by its "raw symdiff never below 2" floor, so
sub-0.10 bins are tiered-only and uncomparable, exactly as anticipated in
the brief. Above ≈0.85 both arms thin out fast (few instances reach that
much disagreement without every draw contradicting). The bins actually
carrying the verdict (0.20–0.55) sit in the well-populated middle of both
arms' achieved range.

---

## J1 — Timing and projection

`time-one` on isolated instances before committing to the full grid
(`run_crossarm.py time-one`):

| component_size | separation | n_k | n_dir_g0 | build | tiered sampling (2000) | flip sampling | r_val | total |
|---|---|---|---|---|---|---|---|---|
| 6 | 3 (median) | 3 | 6 | 3.1 ms | 0.31 s | 0.088 s (600) | 0.3 ms | **0.41 s** |
| 10 | 1 | 4 | 10 | 27.9 ms | 0.70 s | 0.24 s (800) | 0.5 ms | **0.97 s** |
| 10 | 5 (median) | 3 | 10 | 15.6 ms | 0.76 s | 0.21 s (600) | 61.9 ms | **1.04 s** |
| 10 | 9 (max) | 2 | 10 | 40.3 ms | 0.53 s | 0.09 s (400) | 39.6 ms | **0.70 s** |

Per-sample cost ≈0.15–0.35 ms regardless of component size (the dominant
cost is Meek closure + `is_gac_valid_mpdag`, both polynomial in graph size,
not exponential). `r_val` (`breakdown_radius`, `time_limit_s=60`) was under
70 ms in every timed case — nowhere near its budget.

**Projection at the time:** ~1260 target instances (63 cells × 20) ×
~0.5–1.0 s/instance ≈ 10–21 minutes of pure compute. A full-cell dry run
(`component_size=8, separation=4, n_tiers=3`, 20/20 accepted) measured
**12.1 s for 20 instances** (0.6 s/instance), confirming the projection.
**Measured total** (sum of the 9 chunk `elapsed_seconds`, run sequentially):
**944.2 s ≈ 15.7 minutes** of actual compute across the whole grid — within
the projected range, no cut to `instances_per_cell` (kept at 20) or to the
separations grid (kept at the *full* `achievable_separations(c)`, not the
3-point min/median/max subset originally considered for time-safety) was
needed. **N (200 reps) was never touched.**

One structural (not time) shortfall: cell `component_size=10, separation=1,
n_tiers=2` reached only 11/20 accepted instances after exhausting a
400-seed budget (389 of that cell's seeds rejected as
`optimal_set_undefined` — with only 2 tiers, low separation leaves too much
of the CPDAG undirected in `K_ref` for `optimal_adjustment_set_mpdag` to
identify a unique set across extensions). This is a low-yield *region of the
generative process*, not a time cut, and is reported rather than padded.
Every other one of the 63 cells reached the full 20/20.

## J2 — S = 1.0 at intensity 0

Both arms share one `G0`/`Z*` per instance, and `Z*` is checked valid at
`G0` (`is_gac_valid_mpdag(g0, x, y, z_star)`) as an **admission gate** at
instance-build time (`z_invalid_at_g0` rejection reason) — so every
accepted instance satisfies this by construction, for both arms
simultaneously (there is only one `G0` in this design). Verified, not
merely assumed: `s0_ok` is recomputed independently after admission and
written to `xarm_instances.csv`. **Result: `s0_ok = True` for all 1251/1251
accepted instances, 0 failures.** (The swept grids don't sample
`corruption_rate = 0` / `d = 0` directly — tiered starts at 0.05, flip at
`d = 1` — since that state is `G0` itself, checked once per instance rather
than resampled 200 times; sampling it would trivially return `S = 1.0`
by the same admission-gate logic.)

## J3 — Determinism

One small cell (`component_size=6, separation=3, n_tiers=2,
instances_per_cell=3, reps=200`) run twice, isolated from the main output
directory, under `PYTHONHASHSEED=0` and `PYTHONHASHSEED=12345`:

```
sha256(xarm_bins.csv), PYTHONHASHSEED=0:     7322197792891fab6bb81a4ee27815e1366d26c2afa88313c3098888b1c6bdd4
sha256(xarm_bins.csv), PYTHONHASHSEED=12345: 7322197792891fab6bb81a4ee27815e1366d26c2afa88313c3098888b1c6bdd4
```

**Identical.** No columns needed exclusion — `xarm_bins.csv` carries no
wall-clock columns at all (those live only in `xarm_instances.csv`, via
`r_search_seconds` etc., which are not part of this determinism check).
`xarm_samples.csv` was also byte-identical between the two runs, as a bonus
check beyond what was asked.

## J4 — Purity

```
$ grep -n "all_valid_adjustment_sets_mpdag\|synth\.runner\|np\.random\.seed\|random\.seed\|import random" \
    src/bkrobust/robustness/crossarm.py src/bkrobust/robustness/run_crossarm.py
(no output — clean)
```

Neither file imports `bkrobust.synth.runner`, calls
`all_valid_adjustment_sets_mpdag`, touches `np.random.seed`/`random.seed`,
or imports the `random` module. `benchmarks.measure.fast_gate` is the only
gate used, stamped as the literal string `"fast_gate"` on every row. Every
draw goes through `np.random.default_rng(derived_seed(...))` — `derived_seed`
uses `hashlib.sha256`, never Python's salted `hash()`.

## J5 — Achieved intensity histogram (non-contradictory samples only)

Pooled `n_noncontra` per bin, both arms, from `xarm_bins.csv` (full detail
in `_xarm_finalize_summary.json`; only bins with material mass shown):

| bin | tiered | flip |
|---|---:|---:|
| 0.0 | 287,557 | 0 |
| 0.05 | 54,197 | 0 |
| 0.10 | 230,510 | 168 |
| 0.15 | 62,132 | 2,200 |
| 0.20 | 47,217 | 5,859 |
| 0.25 | 86,435 | 9,409 |
| 0.30 | 43,990 | 10,772 |
| 0.35 | 18,480 | 11,619 |
| 0.40 | 33,229 | 16,741 |
| 0.45 | 838 | 4,904 |
| 0.50 | 28,682 | 8,989 |
| 0.55 | 13,521 | 8,818 |
| 0.60 | 5,596 | 13,275 |
| 0.65 | 3,100 | 10,956 |
| 0.70 | 1,766 | 19,425 |
| 0.75 | 1,192 | 18,269 |
| 0.80 | 81 | 7,273 |
| 0.85 | 77 | 15,499 |
| 0.90 | 0 | 3,697 |
| 1.00+ | 29 | ~146,000 |

**As anticipated in the brief:** flip's raw symdiff floor (reversing one
claim changes ≥2 orientations) pushes essentially all of its mass to
intensity ≥ 0.10, with a very long tail past 1.0 (flip can, and often does,
disagree with `G0` on *more* directed edges than `G0` had to begin with —
1.6 was the observed max bin). Tiered is the opposite: 74% of its
non-contradictory mass sits at intensity ≤ 0.10, including a large exact-0
spike (a relocated node's crossing edges can re-derive the same
orientations). **Bins below ≈0.10 support tiered only; flip is essentially
absent there. Bins above ≈0.85 thin out on both sides.** The matched
comparison in the verdict above is restricted to where both arms actually
have mass — the middle of the axis, 0.15–0.55 — exactly per this histogram,
not by post-hoc bin selection to manufacture a result.

## J6 — Exclusions, counted

**Instance-admission exclusions** (across the whole grid, summed from each
chunk's printed summary — see note below on why these live here rather than
in one aggregated JSON):

| reason | count |
|---|---|
| `optimal_set_undefined` | 903 |
| `generate_instance:*`, `fast_gate:*`, `n_k_zero`, `k_ref_contradictory`, `o_g0_extensions_intractable`, `z_invalid_at_g0`, `dir_g0_empty` | 0 (none fired) |

All 903 exclusions were `optimal_set_undefined` (no DAG extension of `G0`
agreed on one optimal adjustment set — expected at low `n_tiers`/low
separation, see J1). Every other guarded rejection reason in
`build_matched_instance` fired zero times across 1251 + 903 = 2154 attempted
seeds.

**Per-sample outcomes:**
- Contradictions (`apply_orientations` → `None`): summed from
  `xarm_bins.csv`'s `n_contradictory` column — **1,583,371 / 2,502,000**
  tiered-arm samples (pooled rate 0.6328) and **637,544 / 951,600** flip-arm
  samples (pooled rate 0.6700) were contradictory, *unconditionally* (not
  matched to intensity). Averaged per-instance instead of pooled by sample
  count, the same two arms give tiered 0.6328 / flip 0.5238 — the two
  unconditional statistics disagree with each other (see the caution note
  in the verdict section), which is exactly why the matched, per-bin table
  there (not either of these two numbers) is the one to trust; it shows
  tiered's contradiction rate higher than flip's at every shared bin in
  [0.10, 0.85).
- Grid points with **zero** non-contradictory draws out of 200 reps
  (`"unresolved_all_contradictory"`, filed as their own bin rather than
  guessed at): 6 (tiered) / 757 (flip) grid-point×instance cells — flip's
  count is large because many instances have small `n_k` and high-`d` grid
  points there reverse most/all claims, contradicting on every one of 200
  reps.
- `r_val == UNREACHED` (`-1`): **0** instances.
- Timeouts (`r_status == "timeout"`, `time_limit_s = 60`): **0** instances.
  Every one of 1251 `r_val` computations returned `"ok"`; max observed
  `r_total_seconds` was well under a second in the J1 timing spot-checks.

**Note on the summary JSON:** a bug in the original `tag = int(t_start)`
(monotonic `perf_counter`, not a wall-clock epoch — arbitrary and can
collide across separate process invocations) meant `_xarm_sweep_summary_*.json`
files from earlier chunks were overwritten by later ones sharing the same
integer tag; only the last chunk's summary survives on disk. The exclusion
counts above were reconstructed from each chunk's own stdout (captured live
during the run, not re-derived after the fact) and cross-checked against
`xarm_instances.csv`'s row count (1251 accepted + 903 rejected = 2154
`build_matched_instance` calls made, seed-budget-bounded per cell). The bug
is fixed in `run_crossarm.py` (tag now `f"{int(time.time()*1000)}_{pid}"`)
for any future run; it affected only a non-contractual diagnostic file, not
any of the seven required outputs, and no data was lost.

## J7 — Touched files

```
$ git status --short
 M results/axis_robustness/PREREGISTRATION.md
?? results/axis_robustness/NULLCELL_RUN_NOTES.md
?? results/axis_robustness/_xarm_finalize_summary.json
?? results/axis_robustness/_xarm_sweep_summary_0.json
?? results/axis_robustness/null_curves.csv
?? results/axis_robustness/null_instances.csv
?? results/axis_robustness/null_manifest.json
?? results/axis_robustness/null_samples.csv.gz
?? results/axis_robustness/null_tau_comparison.csv
?? results/axis_robustness/xarm_bins.csv
?? results/axis_robustness/xarm_instances.csv
?? results/axis_robustness/xarm_manifest.json
?? results/axis_robustness/xarm_paired_test.csv
?? results/axis_robustness/xarm_samples.csv.gz
?? src/bkrobust/robustness/crossarm.py
?? src/bkrobust/robustness/nullcell.py
?? src/bkrobust/robustness/run_crossarm.py
?? src/bkrobust/robustness/run_nullcell.py
```

**The modified file (`PREREGISTRATION.md`) and the `null_*`/`nullcell.py`
files are not mine** — they belong to a different, concurrently-running
task in this same checkout (visible mid-run in `git status`; its diff to
`PREREGISTRATION.md` is a new "Appendix H — the null stratum is real",
about a session-7 τ_b stratum, unrelated to this cross-arm design). This
session never opened, read for editing, or wrote to `PREREGISTRATION.md` or
any `null_*`/`nullcell*` path. Confirmed via `git diff --stat` (only
`PREREGISTRATION.md` shows as modified-not-created, and its diff is
entirely Appendix H content) and by this session's own tool history writing
only `crossarm.py`, `run_crossarm.py`, and `xarm_`/`_xarm_`-prefixed files
under `results/axis_robustness/`.

**Created by this task, all new files, nothing pre-existing touched:**
`src/bkrobust/robustness/crossarm.py`, `src/bkrobust/robustness/run_crossarm.py`,
`results/axis_robustness/xarm_samples.csv.gz`, `xarm_bins.csv`,
`xarm_instances.csv`, `xarm_paired_test.csv`, `xarm_manifest.json`,
`_xarm_finalize_summary.json`, `_xarm_sweep_summary_0.json`, and this file.

---

## Output files and row counts

| file | rows (excl. header) | notes |
|---|---:|---|
| `xarm_samples.csv.gz` | 3,453,600 | gzipped (uncompressed ≈447 MB, well over the 20 MB threshold; compressed 51.7 MB) |
| `xarm_bins.csv` | 9,799 | per (instance × arm × intensity bin, or the `unresolved_all_contradictory` sentinel) |
| `xarm_instances.csv` | 1,251 | one row per accepted matched instance |
| `xarm_paired_test.csv` | 15 | one row per shared intensity bin with ≥1 matched instance pair |
| `xarm_manifest.json` | — | grid, git SHA, environment, gate |

Grid actually run: `component_size ∈ {6, 8, 10}` × `separation ∈
achievable_separations(component_size)` (**full set**, 5/7/9 values — not
cut to a 3-point subset) × `n_tiers ∈ {2, 3, 4}` = **63 cells**, target 20
instances/cell (1260), **1251 accepted** (one cell short at 11/20, see J1/J6).
`N = 200` reps per (instance, arm, grid point) throughout, never reduced.
Tiered grid: `corruption_rate ∈ {0.05, 0.10, ..., 0.50}` (10 points, per the
brief). Flip grid: `d/len(K_ref)` for `d = 1..len(K_ref)`, per instance
(`n_k` ranged 1–11, mean 3.80, across the accepted instances).

`AUC_intensity_usable` (mean `S` over intensity bins in `[0, 1)` with ≥30
non-contradictory samples, per instance, per arm): defined for **1251/1251**
instances on the tiered arm (mean 0.883) but only **1026/1251** on the flip
arm (mean 0.751 where defined) — the other 225 instances' flip sweep never
put ≥30 non-contradictory draws in any bin below 1.0 at all (consistent with
J5: small-`n_k` instances' flip grid can sit entirely above intensity 1.0).
This is itself informative, not a gap to paper over: flip's usable range is
narrower than tiered's for a meaningful fraction of instances, precisely
because flip cannot express low intensity when `n_k` is small (its coarsest
possible corruption, `d=1`, already crosses into high intensity).

---

## What surprised me

1. **The direction of the effect is the opposite of the session-7
   impression, on both fronts.** Session 7's (confounded, later-retracted)
   read was that correlated errors were less self-revealing. This design,
   built specifically to test that cleanly, finds correlated errors *more*
   self-revealing (higher contradiction rate at matched intensity) *and*
   more dangerous when they don't self-reveal (lower `S`). Worth flagging
   for anyone tempted to still informally trust session 7's phrasing even
   after its retraction — the corrected picture isn't just "unsupported,
   no claim," it's "actively the other way" once measured properly.
2. **Speed.** Given the brief's warning about 10-minute budgets and cutting
   instance counts, I expected to have to cut the grid. Per-sample cost
   (~0.15–0.35 ms) made the whole 63-cell, 1260-target-instance,
   200-reps-per-point grid finish in under 16 minutes of compute — no cut
   to `instances_per_cell` or to the separations grid was needed anywhere
   except the one structurally-limited cell.
3. **`n_tiers` matters more than separation for instance yield.** The
   `optimal_set_undefined` rejection is almost entirely an `n_tiers=2`
   phenomenon (406 of 903 total exclusions in the `n_tiers=2` chunks for
   `component_size=10` alone) — coarser tiering leaves too much of `K_ref`
   unoriented to identify a unique optimal set, independent of how far `X`
   is from the adjustment set. `n_tiers=4` at the same low separation
   accepted 20/23 seeds tried versus 11/400 at `n_tiers=2`.
4. **A large fraction of flip's grid is "unresolved" (all-contradictory).**
   757 (instance × grid-point) cells never produced a single non-contradictory
   draw across 200 reps — almost all at high `d` on small-`n_k` instances,
   where reversing most/all asserted claims collides with Meek closure
   essentially every time. This is the flip-side confirmation of the brief's
   warning that the two processes cover the axis unevenly.
