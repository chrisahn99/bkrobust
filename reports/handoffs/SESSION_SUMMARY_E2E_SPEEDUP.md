# Session summary: end-to-end speedup for Section 5.3

Branch `revision/e2e_speedup`, merged into `main`. 2026-09-23.
Repository: <https://github.com/chrisahn99/bkrobust>

This file is a handoff for a session that will edit the paper (the LaTeX
repository, which this session could not touch). It says what was found, what the
paper should change, and where each number lives. **Every number here is in a
committed file listed in §6. Do not quote a number that is not.**

---

## 1. What the session set out to do

Section 5.3 (`\subsection{A hybrid, and what it costs}`, label `sec:hybrid`) ends
with this TODO:

> A single end-to-end speedup figure combining the two components —
> search-versus-space and hybrid-versus-search — on one instance family. The
> components are measured on different sweeps and cannot be composed post hoc.

The two components came from sweeps that share no instances:

| component | source | family |
|---|---|---|
| search vs space (`fig:scaling`, `figures/scaling.pdf`) | `results/axisb2/scaling/` | n = 5–7, k ≤ 10 |
| hybrid vs search (the "256 instances, max 21.6 s" paragraph) | `results/axisb4/hybrid_envelope.jsonl` | ER, n = 8–24 |

**Done:** one sweep measures all legs on the same instances and the same query,
so `e2e = (space/search) × (search/hybrid)` holds exactly for each instance.

## 2. What was built and run

| commit | content |
|---|---|
| `279b011` | `experiments/e2e_speedup.py` (driver) + first, ungated run `results/e2e_speedup/` |
| `e6ff160` | **Dispatch fix** in `src/bkrobust/search/exact_fast.py` |
| `7e46433` | `--gate` option, `experiments/plot_e2e.py`, **final run `results/e2e_speedup_gated/`** + figures |

**Use `results/e2e_speedup_gated/` for the paper.** `results/e2e_speedup/` is the
before-fix record only.

**Design of the final sweep.**
- **Instances.** 256 ER instances: n ∈ {8,10,12,14,16,18,20,24} ×
  p ∈ {0.3, 0.5, 0.7, 0.85} × seeds 0–7, built with
  `bkrobust.gac.sweep.dense_instance`. G0 is the true DAG, with full knowledge.
- **k ≥ 1 gate.** If a draw's CPDAG has no undirected edge (k = 0), the cell is
  redrawn deterministically with seed + 10,000·attempt. 242 cells used attempt 0,
  13 used attempt 1 and 1 used attempt 2.
- **Query.** The first sorted (x, y) pair whose optimal adjustment set is
  non-empty and valid; this is the rule in `src/bkrobust/search/scaling.py`.
- **Four legs per instance**, each in its own subprocess, run one at a time
  (serial), with a 120 s cap per leg:

| leg | what | oracle |
|---|---|---|
| `space` | build the full space (`build_space_correct`) + BFS | criterion (enumeration also timed) |
| `search_frozen` | `radius_local_up`, unbounded; the "search" of the old figures | enumeration |
| `search_fast` | `radius_local_up_fast`, unbounded | criterion |
| `hybrid` | `breakdown_radius`, defaults (budget 3) | criterion |

A capped run is recorded as censored and never carries a time. A ratio whose
numerator was capped is a **lower bound** and sits in a separate `*_lb` column.

**The dispatch fix.** `radius_local_up_fast` marked the depth budget as hit when
it reached any state at that depth, even a state with nowhere further to go. So
a search that had finished within the budget still reported `exact=False`, and
the hybrid ran the SAT ladder for nothing. The fix sets `budget_hit` only when a
state at the budget has an unvisited cover.
- Radii cannot change, because states at the budget are still not expanded.
- Checked on 238 instances × budgets 0–3 with 0 violations; the tests in
  `tests/hybrid` and `tests/search` pass (66/66).
- `src/bkrobust/search/exact.py` has the same pattern and was deliberately left
  alone; it is the frozen reference implementation.

## 3. Results

### Integrity checks (recounted independently of the driver)

- 256 rows, 256 unique cells, **0 rows with k = 0**.
- Radius agreement across every completed leg: **256 / 256**.
- Capped runs: space 22, search_frozen 6, search_fast 2, **hybrid 0**.
- The identity e2e = product of the two component ratios holds on all 234 rows
  where every leg completed.
- **Every instance where the search alone hit the cap has no reachable failure.**
  All 8 (6 frozen, 2 fast) are UNREACHED, and the hybrid answered each one via the
  E1 ladder in 3.4–45.6 s.

### Per-k medians (criterion chain unless marked)

Times are in seconds. Ratios are medians over rows where both legs finished.

| k | rows | space | search_fast | hybrid | space/search | search/hybrid | **e2e** | frozen/hybrid (enum) |
|---|---|---|---|---|---|---|---|---|
| 1 | 32 | 0.0006 | 0.0003 | 0.0003 | 2.3 | 0.99 | **2.3** | 1.98 |
| 2 | 43 | 0.0015 | 0.0004 | 0.0004 | 3.6 | 1.00 | **3.5** | 2.61 |
| 3 | 35 | 0.0070 | 0.0011 | 0.0012 | 5.7 | 1.00 | **5.6** | 3.98 |
| 4 | 35 | 0.020 | 0.0021 | 0.152 | 10.3 | 0.013 | **0.093** | 0.074 |
| 5 | 24 | 0.069 | 0.0024 | 0.153 | 19.2 | 0.019 | **0.37** | 0.15 |
| 6 | 21 | 0.186 | 0.0041 | 0.081 | 44.9 | 0.077 | **2.3** | 0.70 |
| 7 | 12 | 1.00 | 0.0051 | 0.205 | 757 | 0.49 | **690** | 1.02 |
| 8 | 19 | 4.35 | 0.0022 | 0.0022 | 1,054 | 0.89 | **1,043** | 1.38 |
| 9 | 7 | 23.7 | 0.073 | 2.11 | 323 | 0.059 | **13.0** | 1.49 |
| 10 | 8 | 53.9 (6 of 8 finished) | 0.0012 | 0.0013 | 78,971 | 0.74 | **72,959** | 2.07 |
| 11–26 | 20 | capped on every row | — | — | — | 0.70–2.68 | lower bounds only | 1.36–7.13 |

- Hybrid over all 256 instances: median **0.0011 s**, max **45.6 s** (n = 24,
  p = 0.85, k = 26).
- Frozen search / hybrid over all 250 rows where both finished: median **1.92×**.

### What the numbers say

1. **Building the space scales exponentially in k.** The space leg's median time
   grows from 0.0006 s at k = 1 to 53.9 s at k = 10, and it hits the cap on 2 of 8
   instances at k = 10 and on all 20 at k ≥ 11. The searches never grow like that.
2. **The end-to-end speedup is carried by search-versus-space.** It reaches a
   median of about 73,000× at k = 10.
3. **Space/search is not monotone in k** (k = 8: 1,054×; k = 9: 323×;
   k = 10: 78,971×). The space cost rises steadily, but the search cost depends on
   how close the nearest failure is. The current `fig:scaling` caption, "The
   single-query speedup grows exponentially with k", is true of the space cost and
   of the overall trend, but not of every k.
4. **The hybrid's value is completion, not speed.**
   - It finished all 256 instances. The unbounded fast search hit the cap on 2
     and the frozen search on 6, all of them failure-free.
   - Against the fast search, the hybrid is not faster in general. At k = 1–3 they
     tie; at k = 4–6 and 9 the hybrid is slower.
   - Against the frozen search, which is the "search" in the paper's existing
     sentence, the hybrid is faster at k ≤ 3 and at k ≥ 7, and slower at k = 4–6.
5. **The hybrid is slower at k = 4–6 and 9 by design, not because of a bug.**
   There the up-set of G0 is taller than `DEFAULT_SEARCH_BUDGET = 3` (in
   `src/bkrobust/hybrid.py`). The budget is genuinely hit, and the hybrid calls the
   SAT ladder as designed. I re-checked all 78 ladder calls at k ≤ 14 and none was
   wasted. On this family the SAT ladder is simply slower than letting the fast
   search finish. The k = 3 dip in the first run *was* a bug, and is fixed (§2).
6. **These are not the paper's cited instances** (see §5, item 1).

## 4. What needs to change in the paper (Section 5.3)

These are proposals for the paper session; wording is the author's call.

### 4.1 Figure

The figures are in `results/e2e_speedup_gated/figures/`, as PDF (use this),
SVG and PNG.

| figure | shows | suggested use |
|---|---|---|
| `e2e_divergence.pdf` | median seconds vs k for space, search_fast and hybrid, with IQR; individual runs where fewer than 3 completed; triangles where a run hit the cap | **This is the TODO's figure.** Replace `fig:scaling`, or put it beside that figure. |
| `e2e_rescue_scatter.pdf` | search_fast vs hybrid for each instance, coloured by the hybrid's method | Shows that the hybrid ≈ search when the search finishes, and rescues only the capped runs. Appendix, or a second panel. |
| `e2e_censoring_bars.pdf` | capped runs per leg: 22 / 6 / 2 / 0 | Low information; a sentence or a table row may do. |

Suggested caption for `e2e_divergence`: "Median time per query against the number
of undirected edges k, on the same 256 instances for all three methods. Shaded
bands are IQRs, drawn where at least three runs completed; dots are individual
runs; triangles mark k where at least one run hit the 120 s cap. Building the space
grows exponentially in k and no longer completes beyond k = 10; the search and the
hybrid stay below a second on most instances."

### 4.2 Text to update

1. **Remove the `\todo{...}`** at the end of §5.3; it is satisfied.
2. **`fig:scaling` caption** ("The single-query speedup grows exponentially with
   k"): soften to the space cost, or to the trend (§3, item 3). If `fig:scaling`
   is replaced, the "illustrates the first performance advantage" sentence should
   cite the new figure.
3. **The envelope paragraph** ("We evaluated this hybrid method on 256 generated
   instances ranging from 8 to 24 variables … median runtime well under 0.1 s and a
   maximum of 21.6 s"). Two options:
   - **Keep it, citing the old run** (`results/axisb4/hybrid_envelope.jsonl`). It
     is still correct for that run. The new figure would then be on a different
     draw, and the text must say so.
   - **Switch to this sweep (recommended)**, so the text and the figure share one
     family: "The hybrid completed all 256 instances (median 0.0011 s, maximum
     45.6 s). The upward search alone hit a 120 s cap on 6 instances with the
     enumeration oracle and on 2 with the criterion; every one of them has no
     reachable failure."
4. **"The second advantage" paragraph.** It is consistent with the data as long
   as it stays about *completion* on failure-free instances. **Do not** add a
   claim that the hybrid is faster than the search in general: against the fast
   search it is not (§3, items 4–5). If a median speedup is wanted, the defensible
   one is frozen search / hybrid = 1.92× (median, 250 rows), stated with its
   k-dependence.
5. **Optional, honest note:** the budget of 3 makes the hybrid slower than the
   search on mid-sized up-sets (k = 4–6) in this family. That is a performance
   choice, not a correctness issue.

### 4.3 Not to be claimed

- That these are the session-4 envelope instances.
- A single "end-to-end speedup" number without its k. The median over all rows
  mixes regimes and hides the lower-bound rows.
- Any ratio at k ≥ 11 as a measurement. They are lower bounds (min 2.6×, from the
  120 s cap).

## 5. Missing information and open items

1. **Session-4 envelope driver: missing.** The script that produced
   `results/axisb4/hybrid_envelope.jsonl` was never committed. Neither
   `dense_instance` (0/256 match) nor seven seeding variants reproduce its graphs,
   and a k ≥ 1 gate does not either, because the graphs themselves differ. Without
   that script the cited "256 instances, max 21.6 s" cannot be regenerated or
   measured with the space leg. If a copy exists outside the repo, it closes this
   gap.
2. **Search budget: open decision.** The hybrid-slower-than-search dip at k = 4–6
   is governed by `DEFAULT_SEARCH_BUDGET = 3` in `src/bkrobust/hybrid.py`. Raising
   it, or choosing it adaptively from the up-set size, would be a method change and
   was not authorised this session.
3. **No regression test for the dispatch fix.** It was verified by hand (§2), but
   `tests/` was outside this session's file scope. Suggested test: an instance
   whose up-set height equals the budget returns `exact=True`, and radii are
   identical across budgets 0–3.
4. **Lower bounds at k ≥ 11 are weak.** With a 120 s cap the space leg gives e2e
   lower bounds of only 2.6–40,386×. A longer space cap (for example 600 s, a few
   more hours) would tighten them.
5. **The environment differs from earlier sweeps.** `.venv` is Python 3.14.7 with
   numpy 2.5.3; `axisb2`/`axisb4` ran on 3.9.6. Ratios within this sweep are
   sound, but absolute seconds are not comparable across sweeps.
6. **One measurement per leg.** Sub-millisecond times are noisy. This matters for
   the k ≤ 3 ratios (all ≈ 1–6×), not for the space blow-up.
7. **Unrelated, carried over:** `hybrid.py`'s `assumes` string still says
   Anti-Exchange Case B is "verified not proved", while the paper states it is
   proved. This is an open decision from session 8, not touched here.
8. **`results/*` is gitignored.** The session's results were force-added
   (`git add -f`). A future session committing new results must do the same.

## 6. Where everything is (GitHub, `main`)

Base: <https://github.com/chrisahn99/bkrobust/tree/main>

| what | path |
|---|---|
| This report | [`SESSION_SUMMARY_E2E_SPEEDUP.md`](https://github.com/chrisahn99/bkrobust/blob/main/SESSION_SUMMARY_E2E_SPEEDUP.md) |
| **Final data (use this)** | [`results/e2e_speedup_gated/`](https://github.com/chrisahn99/bkrobust/tree/main/results/e2e_speedup_gated) |
| Per-instance flat table with all ratios | [`results/e2e_speedup_gated/e2e_speedup.csv`](https://github.com/chrisahn99/bkrobust/blob/main/results/e2e_speedup_gated/e2e_speedup.csv) |
| Per-instance full record (every leg, counters) | [`results/e2e_speedup_gated/e2e_speedup.jsonl`](https://github.com/chrisahn99/bkrobust/blob/main/results/e2e_speedup_gated/e2e_speedup.jsonl) |
| Per-k aggregates, gate histogram | [`results/e2e_speedup_gated/summary.json`](https://github.com/chrisahn99/bkrobust/blob/main/results/e2e_speedup_gated/summary.json) |
| Provenance (sha, env, grid, cap, command) | [`results/e2e_speedup_gated/manifest.json`](https://github.com/chrisahn99/bkrobust/blob/main/results/e2e_speedup_gated/manifest.json) |
| Findings for the final run | [`results/e2e_speedup_gated/FINDINGS.md`](https://github.com/chrisahn99/bkrobust/blob/main/results/e2e_speedup_gated/FINDINGS.md) |
| **Figures (PDF/SVG/PNG)** | [`results/e2e_speedup_gated/figures/`](https://github.com/chrisahn99/bkrobust/tree/main/results/e2e_speedup_gated/figures) |
| Proof the envelope does not reproduce | [`results/e2e_speedup_gated/query_rule_verification.json`](https://github.com/chrisahn99/bkrobust/blob/main/results/e2e_speedup_gated/query_rule_verification.json) |
| Before-fix record (do not cite) | [`results/e2e_speedup/`](https://github.com/chrisahn99/bkrobust/tree/main/results/e2e_speedup) |
| Sweep driver | [`experiments/e2e_speedup.py`](https://github.com/chrisahn99/bkrobust/blob/main/experiments/e2e_speedup.py) |
| Plot script | [`experiments/plot_e2e.py`](https://github.com/chrisahn99/bkrobust/blob/main/experiments/plot_e2e.py) |
| Dispatch fix | [`src/bkrobust/search/exact_fast.py`](https://github.com/chrisahn99/bkrobust/blob/main/src/bkrobust/search/exact_fast.py), commit `e6ff160` |
| Old component sources (for comparison) | `results/axisb2/scaling/`, `results/axisb4/hybrid_envelope.jsonl` |

**CSV columns the paper session will need:** `k_undirected`, `space_total_crit_s`,
`search_fast_s`, `search_frozen_s`, `hybrid_total_s`, `hybrid_method`,
`*_censored`, `speedup_space_vs_search`, `speedup_search_vs_hybrid`,
`speedup_e2e`, `speedup_e2e_lb` (lower bounds), `speedup_search_vs_hybrid_enum`.

**To regenerate:**

```
PYTHONPATH=src .venv/bin/python experiments/e2e_speedup.py --gate --out results/e2e_speedup_gated   # ~79 min, serial
PYTHONPATH=src .venv/bin/python experiments/plot_e2e.py results/e2e_speedup_gated/e2e_speedup.csv
```
