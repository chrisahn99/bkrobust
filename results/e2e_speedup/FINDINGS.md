# End-to-end speedup sweep (Sec. 5.3 TODO) — findings

Session of 2026-09-23, branch `revision/e2e_speedup`. Driver:
`experiments/e2e_speedup.py`. Every number below is recomputed from
`e2e_speedup.jsonl` / `e2e_speedup.csv` in this directory.

## What was measured

One sweep, one instance family, four legs per instance on the **same** graph and the
**same** query `(x, y, Z)`:

| leg | what | oracle |
|---|---|---|
| `space` | `build_space_correct` + BFS radius (single-query cost = build + query) | enumeration and criterion (both timed) |
| `search_frozen` | `radius_local_up`, unbounded (the "search" of both prior sweeps) | enumeration |
| `search_fast` | `radius_local_up_fast`, unbounded (the hybrid's own search leg) | criterion |
| `hybrid` | `breakdown_radius`, defaults (budget 3, criterion) | criterion |

Family: `dense_instance(n, seed, edge_prob)`, n ∈ {8,10,12,14,16,18,20,24} ×
p ∈ {0.3,0.5,0.7,0.85} × seed 0..7 = 256 instances. Query rule: first `(x, y)` in
sorted-node permutation order with a non-empty, valid `optimal_adjustment_set_mpdag`
(the rule of `src/bkrobust/search/scaling.py`). Each leg in a fresh subprocess, legs
run serially, 120 s wall cap per leg. Total wall time 4618 s.

Composition: `e2e = space / hybrid = (space / search) × (search / hybrid)`, an exact
per-instance identity (verified on all 235 rows where every leg completed).

## Integrity checks (recounted independently of the driver)

- 256 rows, 256 unique instances.
- Radius agreement across every completed leg: **256 / 256**, 0 disagreements.
- Censored legs: space 21, search_frozen 6, search_fast 2, hybrid **0**. No censored
  leg carries a time; ratios with a censored numerator are in `*_lb` columns only.

## Per-k medians (k = undirected edges of the CPDAG; criterion chain)

| k | space / search_fast | search_fast / hybrid | e2e |
|---|---|---|---|
| 0–2 | 1.6–3.6 | ≈ 1.0 | 1.6–3.6 |
| 3–5 | 5.3–19.1 | 0.008–0.019 | 0.044–0.37 |
| 6 | 51.4 | 0.50 | 29.0 |
| 7 | 751 | 0.48 | 666 |
| 8 | 1,051 | 0.83 | 1,044 |
| 9 | 341 | 0.059 | 12.9 |
| 10 | 78,724 | 0.86–0.90 | 71,225 (6 of 8 rows; 2 are lower bounds) |
| ≥ 11 | space censored on every row | — | lower bounds only (min 2.7) |

Against the frozen search (enumeration chain), search_frozen / hybrid is 1.3–2.5 at
k ≤ 2, 0.035–0.15 at k = 3–5, and 1.07–1.97 at k ≥ 6.

Hybrid overall: median 0.0012 s, max 44.4 s.

## Findings

1. **The end-to-end speedup is carried by space-versus-search.** The hybrid is
   decisive over the unbounded fast search on only 2 instances (k = 21, 26), where
   the search hit the cap and the hybrid finished in 27.7 s and 44.4 s.
2. **Hybrid < search at k = 3–5 is a dispatch defect, not SAT overhead.** In
   `src/bkrobust/search/exact_fast.py` (~line 130), any state popped at
   depth == `max_depth` sets `budget_hit`, even when it has no covers. When the
   up-set has height ≤ the budget, the search was in fact exhaustive, yet it returns
   `exact=False` and the hybrid runs the E1 ladder. 105 of the 109 ladder dispatches
   were slower than letting the search finish. Radii are unaffected. **Do not quote a
   search-versus-hybrid number from this sweep until the fix lands and the hybrid
   leg is re-run.**
3. **This is not the session-4 envelope.** `dense_instance` with the envelope's
   `(n, p, seed)` reproduces 0 / 256 of `results/axisb4/hybrid_envelope.jsonl`
   (`query_rule_verification.json`); the envelope has no k = 0 graphs, so its
   (never committed) driver had a resample/gate step. The paper's "256 instances,
   max 21.6 s" still describes the old run.

## Caveats

- Environment: Python 3.14.7 / numpy 2.5.3, unlike the 3.9.6 of `axisb2`/`axisb4`.
  Ratios within this sweep are comparable; absolute seconds across sweeps are not.
- `manifest.json` records `git_sha 05d3680, git_dirty true`: the driver was
  uncommitted when it ran and is committed alongside these results.
- The 120 s cap makes the k ≥ 11 e2e lower bounds weak.
- Single cold measurement per leg; sub-millisecond timings are noisy.

## Open items

1. Session-4 envelope driver, if a copy exists, to run on the cited instances.
2. Dispatch fix in `exact_fast.py` + re-run of the hybrid leg only.
3. Figure script for the Sec. 5.3 figure (data: `e2e_speedup.csv`).
