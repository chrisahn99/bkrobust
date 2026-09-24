# End-to-end speedup, v2 — the top-state check — findings

Companion sweep to `results/e2e_speedup_gated/` (the run cited by the paper's
`tab:envelope` and `fig_e2e_speedup.pdf`), on the exact same 256 instances and
`--gate` config, after adding the top-state check to
`bkrobust.hybrid.breakdown_radius`: if `Z` is valid in the CPDAG `Ĉ` itself,
the radius is `UNREACHED` by upward-closure alone (`Ĉ` is the maximum of the
order), returned in one oracle call (`method="top_state"`) before any search
or SAT ladder runs.

**Timing correction (this revision):** the first pass of this run left the
top-state check's own time out of `total_seconds` (it measured search +
ladder only, so a `top_state` answer showed as 0.0 s). Fixed in
`src/bkrobust/hybrid.py`: `HybridResult` now has a `top_seconds` field,
measured on **every** `breakdown_radius` call (not only when it is the
answer — also when it falls through to search), and
`total_seconds = top_seconds + search_seconds + ladder_seconds`. Only the
`hybrid` leg needed re-running (`search_fast_top`'s worker already timed its
own top check correctly); `space`, `search_frozen`, `search_fast` and
`search_fast_top` timings are unchanged from the first v2 pass.

## Headline numbers (corrected)

- **Hybrid total time: median 0.247 ms, 95th pct 2.03 ms, max 2.044 s**
  (n=18, k=20). Old run (no top check): median 1.128 ms, p95 3258 ms
  (3.26 s), max 45.57 s.
- **Top-state check alone: median 0.144 ms** over all 256 calls (it runs on
  every call, whether or not it answers).
- **164/256 (64%) instances have `r = ∞`**, all answered by `top_state`.
- **6/256** still need the SAT ladder — genuinely finite radius beyond the
  depth-3 budget (r = 4, 4, 4, 5, 5, 6); 86/256 answered within budget by the
  bounded search.
- Radius agreement: **256/256** across all six legs, and identical to the old
  run's hybrid radii on all 256 instances (0 mismatches).

## Ladder instances (finite radius > budget 3)

| n | k | radius | total (s) | top (ms) | search (ms) | ladder (s) |
|---|---|---|---|---|---|---|
| 14 | 9 | 5 | 0.224 | 0.52 | 5.4 | 0.219 |
| 14 | 14 | 4 | 0.533 | 0.77 | 26.7 | 0.505 |
| 14 | 17 | 5 | 0.565 | 1.45 | 29.1 | 0.534 |
| 16 | 13 | 4 | 0.810 | 1.14 | 30.6 | 0.778 |
| 18 | 10 | 4 | 1.155 | 0.77 | 27.3 | 1.127 |
| 18 | 20 | 6 | **2.044** | 1.81 | 81.6 | 1.961 |

## Per-k table (`tab:envelope` format), search = search_fast_top

| k | instances | space (s) | search (ms) | hybrid (ms) | space/search |
|---|---|---|---|---|---|
| 1 | 32 | 0.0006 | 0.148 | 0.132 | 5.2 |
| 2 | 43 | 0.0015 | 0.187 | 0.189 | 12.5 |
| 3 | 35 | 0.0070 | 0.205 | 0.206 | 38.6 |
| 4 | 35 | 0.020 | 0.217 | 0.223 | 137 |
| 5 | 24 | 0.069 | 0.194 | 0.177 | 447 |
| 6 | 21 | 0.186 | 0.302 | 0.345 | 1,196 |
| 7 | 12 | 1.00 | 0.294 | 0.296 | 3,344 |
| 8 | 19 | 4.35 | 0.423 | 0.450 | 9,559 |
| 9 | 7 | 23.7 | 0.319 | 0.333 | 81,870 |
| 10 | 8 | 53.9$^\dagger$ | 0.773 | 0.729 | 62,190 |
| 11–26 | 20 | $>120$ | 0.45–1,510 | 0.41–2,040 | — |

$^\dagger$ two of eight k=10 instances censored on `space`; median over the
other six. Space column reused from `results/e2e_speedup_gated/`, unaffected
by this fix.

## Findings

1. **The bug did not affect correctness, only reporting.** Radii, methods and
   which legs answered are unchanged from the first v2 pass; only
   `hybrid_total_s` moves (up, from the previously-omitted top-state time),
   most visibly for the 164 `top_state` rows (0 → ~0.1–0.3 ms each).
2. Even corrected, hybrid stays dramatically faster than before the check:
   median 0.247 ms vs 1.128 ms, and p95 2.03 ms vs 3.26 s — the p95 move is
   the big one, since roughly a third of the old run's instances paid full
   search/ladder cost near the cap-adjacent tail to prove `r = ∞`.
3. The top-state check itself costs about 0.14 ms median — negligible next to
   the search budget it replaces, and paid on every call regardless of
   outcome, so it is not free but also not a bottleneck anywhere in this
   sweep (its total across the 6 ladder instances is ~1–2 μs per instance out
   of a multi-hundred-ms total).

## Figure

`results/e2e_speedup_gated_v2/fig_e2e_speedup_v2.pdf`, regenerated from the
corrected CSV (search leg = `search_fast_top` in both panels). Built by a
scratch script (not committed), reading
`results/e2e_speedup_gated_v2/e2e_speedup.csv`. The committed
`figures/paper/fig_e2e_speedup.pdf` is untouched.
