# Zero-radius audit (reviewer Point 7)

**Track C. Owns `experiments/zero_radius_audit.py` and `results/zero_radius/` only.**
Does not modify `experiments/final_table.py`, `src/bkrobust/hybrid.py`, or anything
under `results/final_table/`; those are read-only inputs here.

## Interpreters used (two different ones; see below)

- **This audit's own recomputation** (`experiments/zero_radius_audit.py`, and every
  ad-hoc check quoted in this file) ran under `PYTHONPATH=src /usr/bin/python3`
  (Python 3.9.6), per this track's environment instruction. It imports only
  `bkrobust.demo.*` and `bkrobust.mpdag_criterion.*` -- graph/GAC recomputation
  at G0 -- and never `bkrobust.epsilon.certify`, which uses `zip(..., strict=True)`
  (Python 3.10+ only, `certify.py:261`) and cannot even be imported under 3.9.6.
  `results/axisa3/instances.jsonl` (the 831-row real-corpus baseline re-verified
  below) was likewise produced by 3.9-compatible code.
- **`results/final_table/` and `results/final_table_eps_gt1/`** (read here, not
  recomputed) were produced by `experiments/final_table.py`, which *does* call
  `certify()` and therefore needed `.venv/bin/python` (3.10+) to run at all. This
  audit does not re-run that script and does not need the venv anywhere.

Both interpreters agree on everything this audit checked: every recomputed fact
below (Z, its GAC-validity, amenability) is a pure graph computation shared by
both environments, and where it can be cross-checked against a `certify()`-derived
number (the committed `r_val`), it matches exactly (see the crosstab below).

## What was already established, and re-verified here

1. **The 831-row real corpus has none of this problem.**
   `results/axisa3/instances.jsonl` has 116,094 raw rows (most `admissible: false`,
   i.e. no causal path -- not part of the radius corpus at all). Restricting to
   `admissible: true` gives exactly **831** rows, all with a finite `radius`
   (`radius != -1`), **zero** with `radius == 0`, **zero** with `z_size == 0`, and
   a minimum finite radius of **1**. Re-verified directly by this script's
   `check_real_corpus_baseline()`, not assumed. This 831-row corpus is unrelated
   to the reviewer's 33 rows: those live only in the *LLM-elicited-knowledge*
   table, `results/final_table/`.

2. **The 33 rows are already labelled `degenerate` and already excluded from
   every median.** `experiments/final_table.py::summarise_network` splits
   `degenerate = [r for r in solved if r["r_val"] == 0]` out from `ok` *before*
   any `_cell`/median is computed (final_table.py:356-361), and
   `docs/PAPER_INTEGRATION_FINAL_TABLE.md` Sec 3.3/5 says so explicitly ("Degenerate
   queries are not fragile queries... Excluded from every median"). So the
   reviewer's literal worry -- that these are silently averaged into a robustness
   statistic -- is not what is happening. What this audit adds is a finer
   diagnosis of *why* `r_val = 0`, because "already invalid at G0" turns out to
   bundle two structurally different failures, one of which is a real defect.

## Audit scope and method

The audit re-examines every **non-informative** instance in
`results/final_table/instances.jsonl`:

- **33** rows with `status == "ok"` and `r_val == 0` ("degenerate")
- **27** rows with `status == "error"` (all 27 carry the identical message
  `"no adjustment set is identified at G0; nothing to certify"` -- verified;
  no other error message occurs)

= **60** non-informative instances, out of 160 total. (The remaining 100 are 66
informative + 34 timeout; see the reconciliation table below.)

`results/final_table_eps_gt1/instances.jsonl` was checked for schema and
instance-set consistency: it is the **same 160 (network, x, y, condition)
instances** (same seed, same `D_LLM` knowledge, a wider epsilon grid only). 158/160
rows have an identical `(status, r_val)`; the other 2 (`munin4`, two queries)
flip from `("ok", 0)` to `("timeout", None)` -- a timing-boundary flake near the
120s cap (recomputing one of those same two rows here took 130-144s across two
runs of this very script, so the flake is unsurprising and not a classification
disagreement). The audit uses `results/final_table/` as canonical and does not
double-count the `eps_gt1` copies.

For each of the 60 instances, the network's (DAG, CPDAG) is rebuilt exactly as
`experiments/final_table.py::load_networks` does, `K` is loaded from
`results/elicit/knowledge.json["D_LLM"]`, and `G0 = apply_orientations(cpdag, K)`
is rebuilt exactly as `certify()` does. Then, independently of any committed
label or `r_val`:

- `is_amenable(G0, x, y)` (`bkrobust.mpdag_criterion`)
- `Z = optimal_adjustment_set_mpdag(G0, x, y)` (`bkrobust.demo.evaluate`, the same
  function `certify()` calls to derive its default `z`)
- if `Z` is not `None`: `is_valid_mpdag(G0, x, y, frozenset(Z))` and, if that is
  `False`, `why_invalid(G0, x, y, frozenset(Z))` for the specific reason

Each instance is placed into exactly one class from those facts alone:

| Class | Condition |
|---|---|
| **VALID_EMPTY** | `Z` identified, `|Z| == 0`, and `is_valid_mpdag(...) == True` |
| **NO_SET_EXISTS** | `optimal_adjustment_set_mpdag` returns `None` (no consistent DAG extension, or extensions disagree on the optimal set) |
| **INVALID_RETURNED** | `Z` identified but `is_valid_mpdag(...) == False` -- a precondition violation |

## Results: the exact counts

```
VALID_EMPTY      =  0
NO_SET_EXISTS    = 27   (matches the committed "no identified Z" count exactly)
INVALID_RETURNED = 33   (matches the committed "degenerate" count exactly)
```

Crosstab of committed label vs. recomputed class (perfect agreement -- every
committed `degenerate` row recomputes to `INVALID_RETURNED`, every committed
`no_identified_z` row recomputes to `NO_SET_EXISTS`, with zero rows landing
anywhere else):

```
degenerate      -> INVALID_RETURNED : 33
no_identified_z -> NO_SET_EXISTS    : 27
```

**There is no `VALID_EMPTY` instance among the 60 non-informative rows.** Every
`r_val == 0` row really is a precondition violation on recomputation, not a
legitimate empty-and-valid query mislabelled as degenerate.

Context (not part of the 60-row audit scope, computed as a cross-check): among
the **66 informative** rows (`status == "ok"`, `r_val != 0`), **27/66** have a
committed `Z` that is itself empty (`|Z| == 0`) *and* GAC-valid, with a defined,
nonzero radius. This demonstrates directly that "empty adjustment set" and
"degenerate query" are different properties in this corpus: `VALID_EMPTY` is a
real, populated, and unremarkable category -- it simply never occurs among the
33+27 flagged rows.

### The UNREACHED sentinel is orthogonal to this split (flagged during review; addressed here)

`r_val == -1` is `bkrobust.core.conventions.UNREACHED`, a sentinel meaning "no
reachable perturbation invalidates Z" -- never a radius, never averaged
(`PAPER_INTEGRATION_FINAL_TABLE.md` Sec 5). It occurs *inside* the 66-row
"informative" bucket (`status == "ok"` and `r_val != 0` is true for `r_val == -1`
too), not inside either of the two buckets this audit classifies:

```
informative (66) = 53 finite-radius rows + 13 UNREACHED-sentinel rows
```

Re-verified directly from `results/final_table/instances.jsonl`: of 99
`status == "ok"` rows, 13 have `r_val == -1`, 33 have `r_val == 0`, and 53 have
`r_val > 0` (distribution `{1: 46, 2: 2, 3: 3, 4: 1, 11: 1}`, matching
`PAPER_INTEGRATION_FINAL_TABLE.md` Sec 3.2's table exactly). The 33
`r_val == 0` rows and the 27 `status == "error"` rows can never be `-1` by
construction (`degenerate` is specifically `r_val == 0`, distinct from
`UNREACHED == -1`; `error` rows have `r_val = None`), so no UNREACHED row is or
could be folded into `VALID_EMPTY` / `NO_SET_EXISTS` / `INVALID_RETURNED`. This
is recorded explicitly in `summary.json` under `unreached_sentinel` so a reader
never mistakes "66 informative" for 66 homogeneous finite radii.

## Root cause of the 33 `INVALID_RETURNED` rows (the actual finding)

This is a genuine, reproducible defect, and it decomposes into exactly two
mechanisms, both invisible to `optimal_adjustment_set_mpdag`'s own agreement
check:

**Mechanism 1 -- non-amenable G0 (8/33).** `is_amenable(G0, x, y) == False` for
these 8 (`alarm`, `ecoli70` x2, `hailfinder` (1 of its 5), `hepar2`, `munin1` x2,
`munin4` (1 of 3), `win95pts`). Amenability is an *MPDAG-level* property about
whether undirected edges leave the causal direction out of `x` ambiguous. Every
individual fully-oriented DAG extension is trivially amenable (there are no
undirected edges left to be ambiguous about), so per-extension analysis is
structurally blind to this precondition. `optimal_adjustment_set_mpdag` only
checks that the per-extension optimal sets *agree*, never that `G0` itself is
amenable -- so it can (and here, does) return an agreed-upon `Z` for a query
`is_valid_mpdag` correctly rejects outright via `why_invalid == "not_amenable"`
(checked before `Z` is even examined, per `criterion.py`'s priority order).

**Mechanism 2 -- no causal path from x to y at all (25/33).** For the remaining
25, `G0` *is* amenable, `Z = {}` is returned, and `why_invalid == "open_noncausal_path"`.
Recomputing `possibly_causal_paths(G0, x, y)` (`bkrobust.mpdag_criterion.paths`)
for all 25 gives **zero possibly-causal paths in every one of them** -- `X` does
not (possibly) cause `Y` at all under this MPDAG and this knowledge. The
Henckel-Perkovic-Maathuis optimal-adjustment-set formula implemented in
`optimal_adjustment_set_dag` (`src/bkrobust/demo/evaluate.py:255-285`) is
`O = pa(cn(x,y)) \ (cn(x,y) union {x})`; when `cn(x,y)` (the causal nodes) is
empty -- as it is, trivially, in *every* DAG extension whenever there is no
causal path at all -- the formula degenerates to `O = {}` in every extension.
Agreement across extensions is then automatic, and `optimal_adjustment_set_mpdag`
returns `{}` -- without ever checking whether the empty set actually blocks the
(here, necessarily entirely non-causal, since none of the paths are causal)
confounding paths that remain between `X` and `Y`. It generally does not: that
is exactly what `why_invalid == "open_noncausal_path"` reports on independent
recomputation.

**Both mechanisms share one root cause.** `optimal_adjustment_set_mpdag`
(`src/bkrobust/demo/evaluate.py:288-321`) determines `Z` purely by checking that
`optimal_adjustment_set_dag`'s per-extension answer agrees across every DAG
extension of `G0`. It never calls `is_valid_mpdag` or `is_amenable`
(`bkrobust.mpdag_criterion`) to confirm the agreed-upon set is actually
MPDAG-level GAC-valid before returning it. `certify()`
(`src/bkrobust/epsilon/certify.py:218-222`) then trusts any non-`None` return
value unconditionally:

```python
if z is None:
    derived = optimal_adjustment_set_mpdag(g0, x, y)
    if derived is None:
        raise ValueError("no adjustment set is identified at G0; nothing to certify")
    z = frozenset(derived)
```

`hybrid.breakdown_radius` (`src/bkrobust/hybrid.py` around line 224) is **not**
the bug -- it is the safety net that catches the mismatch, by running the same
`is_valid_mpdag` check on `G0` at search depth 0 and reporting `method =
"degenerate"` / `radius = 0` when it fails. That detection is correct and
working as designed. The defect is upstream: `optimal_adjustment_set_mpdag` can
hand `certify()` (and hence `breakdown_radius`) a `Z` that a one-line check
against the repository's own validated GAC criterion (`is_valid_mpdag`, proven
to agree with the enumeration oracle exhaustively -- see
`mpdag_criterion/criterion.py`'s own docstring) would have rejected. The 27
"no identified Z" rows show the codebase already has the right behaviour for
"nothing usable was found" (raise / report unidentified); the 33 degenerate
rows are the same situation reached by a different code path that fails to
raise.

**Important scope limitation, stated precisely so this is not overstated.**
This audit establishes that the *specific* `Z` returned by
`optimal_adjustment_set_mpdag` fails GAC at `G0` for all 33 rows. It does
**not** establish that no valid adjustment set exists at all for these queries
(a different, nonempty `Z` might satisfy GAC even where the "optimal" one does
not, in the 25 open-path cases in particular -- though not in the 8
not-amenable cases, where `is_valid_mpdag`'s "not_amenable" verdict is
unconditional on `Z` and rules out every possible set). So while the report
below recommends merging the *symptom* handling of these 33 with the existing
27 "unidentified" bucket, it should not be silently reclassified as
`NO_SET_EXISTS` without that check; per the task's own class (c) definition,
this is exactly a "genuine implementation defect" and must be surfaced as one,
not quietly filed as "unidentified."

## Reconciliation against the committed counts

| Committed | n | Recomputed here |
|---|---|---|
| degenerate | 33 | 33/33 recompute to `INVALID_RETURNED` (0 to any other class) |
| no identified Z | 27 | 27/27 recompute to `NO_SET_EXISTS` (0 to any other class) |
| timeout | 34 | out of audit scope (unchanged; not touched) |
| informative | 66 | out of audit scope for the 3-way split; 53 finite + 13 UNREACHED (both re-verified) |
| **total** | **160** | **160** (33+27+34+66, exactly reproduced) |

`summary.json`'s `reconciliation_against_committed_counts.matches_committed_counts`
is `true`.

## What the final output arrays should do differently

1. **`VALID_EMPTY` (n = 0 among the flagged rows, but n = 27 among the
   informative rows already in the table).** These belong in the analysis with
   a defined radius exactly as they are today. Nothing changes here; this
   audit's contribution is confirming that "empty Z" is not itself a red flag.

2. **`NO_SET_EXISTS` (n = 27, unchanged).** Correctly excluded from every
   median already, and correctly reported as "unidentified" (the `error`
   status, with `ValueError: no adjustment set is identified at G0`). Belongs
   in the denominator as "unidentified," never as `radius = 0`. No change
   needed to how these are already handled.

3. **`INVALID_RETURNED` (n = 33).** These are currently reported as
   `method = "degenerate"`, `r_val = 0`, and are already excluded from medians
   -- so *numerically* nothing here is wrong in the paper's committed medians.
   But the *label* is imprecise in a way worth fixing before this is written up
   for the paper: "the optimal adjustment set is already empty and invalid" is
   true but reads as if the pipeline computed a real, if unlucky, answer. What
   actually happened is that `optimal_adjustment_set_mpdag` handed `certify()`
   a `Z` its own validated criterion rejects, for two specific, now-understood
   reasons (non-amenable `G0`, or `X` has no causal path to `Y` at all). The
   recommended fix -- **not implemented here, out of this track's file
   ownership** (only `src/bkrobust/demo/evaluate.py` would need a one-line
   `is_valid_mpdag` check added to `optimal_adjustment_set_mpdag` before it
   returns a non-`None` `Z`) -- would turn these 33 into the same `ValueError`
   /"unidentified" path the other 27 already take, collapsing the "degenerate"
   method in `hybrid.py` to a pure defensive fallback that should then never
   actually fire from `certify()`'s default-`z` path. Until that fix lands,
   these 33 must be surfaced explicitly as what they are -- a known, diagnosed
   implementation gap in optimal-adjustment-set derivation, not a property of
   the data or a mere "ill-posed query" -- rather than filtered silently. A
   follow-up task has been flagged for this (see the spawned suggestion).

## Bottom line

- **A genuine class-(c) precondition violation does exist**, and it accounts
  for all 33 of the reviewer's flagged rows: `INVALID_RETURNED = 33`,
  `NO_SET_EXISTS = 27`, `VALID_EMPTY = 0` among the 60 non-informative rows.
- It is **not** a numerical error in the paper's committed medians (the 33 were
  already excluded, as documented). It **is** a mislabelling: "degenerate" /
  "already invalid at G0" describes the symptom `hybrid.py` correctly detects,
  not the cause, which is a missing validity check in
  `optimal_adjustment_set_mpdag` (`src/bkrobust/demo/evaluate.py`) triggered by
  two identified, disjoint mechanisms (8 non-amenable, 25 no-causal-path).
- The 831-row real corpus (`results/axisa3/`) has none of this: 0 zero-radius,
  0 empty-Z, min radius 1 -- re-confirmed here, not merely assumed. The defect
  is specific to the LLM-elicited-knowledge table's particular (network, x, y, K)
  combinations, where the LLM's assertions happen to produce non-amenable or
  causally-disconnected `(x, y)` pairs at a nontrivial rate (33/160 = 21% of
  this table's rows).
