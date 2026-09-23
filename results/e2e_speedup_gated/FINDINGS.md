# End-to-end speedup, gated sweep after the dispatch fix — findings

Session of 2026-09-23, branch `revision/e2e_speedup`. Supersedes
`results/e2e_speedup/` for the Sec. 5.3 figure; that directory stays as the
before-fix record. Every number below was recomputed from `e2e_speedup.jsonl` /
`e2e_speedup.csv` in this directory, independently of the driver's `summary.json`.

## What changed since `results/e2e_speedup/`

1. **Dispatch fix** in `src/bkrobust/search/exact_fast.py`: a state popped at
   `depth == max_depth` now sets `budget_hit` only if it has an unseen cover. An
   up-set exhausted within the budget returns `exact=True`, so the hybrid no longer
   hands it to the E1 ladder. Radii cannot change (states at the budget are still
   not expanded); verified on 238 instances × budgets 0–3, 0 violations, and
   `tests/hybrid` + `tests/search` 66/66 pass. `src/bkrobust/search/exact.py` (the
   frozen reference) has the same pattern and was deliberately left untouched.
2. **k ≥ 1 gate** in `experiments/e2e_speedup.py --gate`: each grid cell
   `(n, p, seed)` tries attempt 0 = `dense_instance(n, seed, p)`, then effective
   seed `seed + 10_000·a`, and accepts the first draw with `k_undirected ≥ 1` and a
   valid query. Attempts used: 242 cells at 0, 13 at 1, 1 at 2, none exhausted.
   **This is not the session-4 envelope** (`results/axisb4/hybrid_envelope.jsonl`):
   the graphs themselves differ, not only the k = 0 cells, and the envelope's
   driver was never committed. It is a k ≥ 1-gated draw of the same grid.

Same legs, cap (120 s), serial subprocess timing and censoring representation as
before. Total wall time 4733 s.

## Integrity checks

- 256 rows, 256 unique cells, **0 rows with k = 0** (min k = 1).
- Radius agreement across every completed leg: 256 / 256.
- Censored: space 22, search_frozen 6, search_fast 2, hybrid 0. No censored leg
  carries a time. `e2e = (space/search) × (search/hybrid)` exact on all 234 rows
  where every leg completed.
- The run log shows one clean pass 1 → 256 (the worker restarted once, from
  scratch, at 19/256 to add a column; no rows mix the two starts).
- Every E1-ladder dispatch at k ≤ 14 (78 of 85) re-checked: the budget-3 search is
  genuinely non-exhaustive on all 78, **0 wasted dispatches** remain. The 7 at
  k > 14 were not re-checked (cost).

## Per-k medians (criterion chain)

| k | rows | search_fast / hybrid (before fix → after) | e2e | ladder share |
|---|---|---|---|---|
| 1 | 32 | 0.996 → 0.992 | 2.29 | 0 / 32 |
| 2 | 43 | 0.995 → 0.996 | 3.53 | 0 / 43 |
| 3 | 35 | **0.008 → 0.996** | 5.56 | 0 / 35 |
| 4 | 35 | 0.014 → 0.013 | 0.093 | 22 / 35 |
| 5 | 24 | 0.019 → 0.018 | 0.37 | 18 / 24 |
| 6 | 21 | 0.50 → 0.077 | 2.29 | 11 / 21 |
| 7 | 12 | 0.48 → 0.49 | 690 | 6 / 12 |
| 8 | 19 | 0.83 → 0.89 | 1,043 | 9 / 19 |
| 9 | 7 | 0.059 → 0.059 | 13.0 | 6 / 7 |
| 10 | 8 | 0.86 → 0.74 | 72,959 | 3 / 8 |

"Before" medians are from `results/e2e_speedup/`, a different (ungated) draw, so
only the k = 3 change is attributable to the fix alone. Hybrid overall: median
0.0011 s, max 45.6 s.

## Findings

1. **The k = 3 artefact is gone.** It was the dispatch defect.
2. **The k = 4–5 (and 6, 9) dip is not an artefact.** There the up-set is taller
   than the default budget of 3, so the ladder dispatch is the hybrid working as
   designed; on this family the unbounded fast search is simply cheaper than the
   SAT ladder. It is a property of `DEFAULT_SEARCH_BUDGET = 3` in
   `src/bkrobust/hybrid.py`, a performance choice, not a bug. Changing it is a
   method change and was out of scope.
3. **The e2e speedup is carried by space-versus-search.** The hybrid rescues the
   fast search on 2 instances, where it hit the cap (see `e2e_rescue_scatter`).
   The hybrid's guarantee here is 0 / 256 censored, not speed.

## Figures (`figures/`, PDF + SVG + PNG, `experiments/plot_e2e.py`)

- `e2e_divergence`: per-k median seconds with IQR for space, search_fast and
  hybrid; drawn only where ≥ 3 runs completed and a majority completed. Sparser k
  shows the individual runs; triangles mark k with ≥ 1 run censored at the cap.
- `e2e_rescue_scatter`: search_fast vs hybrid per instance, coloured by the
  hybrid's method; censored search_fast drawn at the cap as lower bounds.
- `e2e_censoring_bars`: censored runs per leg out of 256.

## Caveats

Python 3.14.7 / numpy 2.5.3 (absolute seconds not comparable to `axisb2`/`axisb4`);
single cold measurement per leg; 120 s cap makes the k ≥ 11 e2e lower bounds weak.
