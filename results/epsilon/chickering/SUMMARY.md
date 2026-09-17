# Chickering strong-form check -- L1, L2, L3, L4, AE-B

Computational check of the repair proposed for the gap flagged in `THEOREMS.md` §24 / `docs/R_EPSILON_THEORY.md` §6: whether the covered-edge-reversal path between two DAG extensions of an MPDAG `G` can be forced to stay inside `[G]` (not merely inside `[C-hat]`) by using the *strong*, exact-`m`-step form of Chickering's theorem, restricted at each step to reversing a currently-differing edge.

- interpreter: `/Users/ahn/Documents/Research/iclr27/bkrobust/.venv/bin/python`
- python: `3.14.7 (main, Aug  5 2026, 10:29:49) [Clang 21.0.0 (clang-2100.1.1.101)]`
- git HEAD: `5b90a0cccd73a975d270952c4da0a69db7f723ba`
- seed: `20260916`
- total elapsed: 240.9s
- total violations across all checks: **0**

## Results

| check | checked | violations | scope |
|---|---|---|---|
| L1_equivalence_preservation | 50758 | 0 | n=3: exhaustive (all_dags); n=4: exhaustive (all_dags); n=5: exhaustive (all_dags) |
| L1_equivalence_preservation_n6_sample | 7034 | 0 | n=6: random sample, n_samples=4000, seed=20260916 |
| L2_strong_form_exact_m_path | 206878 | 0 | n=3: exhaustive (all equivalence classes, all ordered pairs); n=4: exhaustive (all equivalence classes, all ordered pairs); n=5: exhaustive (all equivalence classes, all ordered pairs) |
| L3_path_stays_in_G | 2907242 | 0 | n=3: exhaustive; n=4: exhaustive; n=5: exhaustive (all 8,782 CPDAGs on 5 nodes) |
| L4_extension_connectivity | 118643 | 0 | n=3: exhaustive; n=4: exhaustive; n=5: exhaustive (all 8,782 CPDAGs on 5 nodes) |
| AE_B_semantic | 293814 | 0 | n=3: exhaustive; n=4: exhaustive; n=5: exhaustive (all 8,782 CPDAGs on 5 nodes) |

## Per-check scope detail

### L1_equivalence_preservation

```json
{
  "3": {
    "n_covered_edge_reversals_checked": 30,
    "n_dags": 25,
    "n_violations_this_n": 0,
    "note": "exhaustive (all_dags)"
  },
  "4": {
    "n_covered_edge_reversals_checked": 828,
    "n_dags": 543,
    "n_violations_this_n": 0,
    "note": "exhaustive (all_dags)"
  },
  "5": {
    "n_covered_edge_reversals_checked": 49900,
    "n_dags": 29281,
    "n_violations_this_n": 0,
    "note": "exhaustive (all_dags)"
  }
}
```

No violations found within the stated scope.

### L1_equivalence_preservation_n6_sample

```json
{
  "6": {
    "n_covered_edge_reversals_checked": 7034,
    "n_dags": 4000,
    "n_violations_this_n": 0,
    "note": "random sample, n_samples=4000, seed=20260916"
  }
}
```

No violations found within the stated scope.

### L2_strong_form_exact_m_path

```json
{
  "3": {
    "n_classes_nontrivial": 7,
    "n_classes_total": 11,
    "n_ordered_pairs_checked": 54,
    "n_violations_this_n": 0,
    "note": "exhaustive (all equivalence classes, all ordered pairs)"
  },
  "4": {
    "n_classes_nontrivial": 126,
    "n_classes_total": 185,
    "n_ordered_pairs_checked": 2424,
    "n_violations_this_n": 0,
    "note": "exhaustive (all equivalence classes, all ordered pairs)"
  },
  "5": {
    "n_classes_nontrivial": 6166,
    "n_classes_total": 8782,
    "n_ordered_pairs_checked": 204400,
    "n_violations_this_n": 0,
    "note": "exhaustive (all equivalence classes, all ordered pairs)"
  }
}
```

No violations found within the stated scope.

### L3_path_stays_in_G

```json
{
  "3": {
    "G_note": "exhaustive (all elements of enumerate_space_correct)",
    "cpdag_note": "exhaustive",
    "cpdags_available": 11,
    "cpdags_examined": 11,
    "n_ordered_pairs_checked": 114,
    "n_violations_this_n": 0,
    "pair_note": "exhaustive (all ordered pairs of extensions)",
    "total_G_available": 50
  },
  "4": {
    "G_note": "exhaustive (all elements of enumerate_space_correct)",
    "cpdag_note": "exhaustive",
    "cpdags_available": 185,
    "cpdags_examined": 185,
    "n_ordered_pairs_checked": 12708,
    "n_violations_this_n": 0,
    "pair_note": "exhaustive (all ordered pairs of extensions)",
    "total_G_available": 1601
  },
  "5": {
    "G_note": "exhaustive (all elements of enumerate_space_correct)",
    "cpdag_note": "exhaustive (all 8,782 CPDAGs on 5 nodes)",
    "cpdags_available": 8782,
    "cpdags_examined": 8782,
    "n_ordered_pairs_checked": 2894420,
    "n_violations_this_n": 0,
    "pair_note": "exhaustive (all ordered pairs of extensions)",
    "total_G_available": 116992
  }
}
```

No violations found within the stated scope.

### L4_extension_connectivity

```json
{
  "3": {
    "G_note": "exhaustive (all elements of enumerate_space_correct)",
    "cpdag_note": "exhaustive",
    "cpdags_examined": 11,
    "n_G_graphs_checked": 50,
    "n_violations_this_n": 0,
    "total_G_available": 50
  },
  "4": {
    "G_note": "exhaustive (all elements of enumerate_space_correct)",
    "cpdag_note": "exhaustive",
    "cpdags_examined": 185,
    "n_G_graphs_checked": 1601,
    "n_violations_this_n": 0,
    "total_G_available": 1601
  },
  "5": {
    "G_note": "exhaustive (all elements of enumerate_space_correct)",
    "cpdag_note": "exhaustive (all 8,782 CPDAGs on 5 nodes)",
    "cpdags_examined": 8782,
    "n_G_graphs_checked": 116992,
    "n_violations_this_n": 0,
    "total_G_available": 116992
  }
}
```

No violations found within the stated scope.

### AE_B_semantic

```json
{
  "3": {
    "G_note": "exhaustive (all elements of enumerate_space_correct)",
    "cpdag_note": "exhaustive",
    "cpdags_examined": 11,
    "n_triples_checked": 12,
    "n_violations_this_n": 0,
    "total_G_available": 50
  },
  "4": {
    "G_note": "exhaustive (all elements of enumerate_space_correct)",
    "cpdag_note": "exhaustive",
    "cpdags_examined": 185,
    "n_triples_checked": 1632,
    "n_violations_this_n": 0,
    "total_G_available": 1601
  },
  "5": {
    "G_note": "exhaustive (all elements of enumerate_space_correct)",
    "cpdag_note": "exhaustive (all 8,782 CPDAGs on 5 nodes)",
    "cpdags_examined": 8782,
    "n_triples_checked": 292170,
    "n_violations_this_n": 0,
    "total_G_available": 116992
  }
}
```

No violations found within the stated scope.

## Reading these results

- Verification is not proof. Every check above (L1, L2, L3, L4, AE-B) is run **exhaustively** for n=3, 4 and 5 -- all DAGs / all CPDAGs / all elements of `enumerate_space_correct` / all ordered pairs of extensions, with no per-CPDAG or per-G sampling caps anywhere in this run (the full n=5 sweep, including all 8,782 CPDAGs and all 116,992 knowledge states, completed in about four minutes, well inside the 30-minute budget). L1 additionally spot-checks a random sample of 4,000 DAGs on n=6 nodes, documented per-check above. No CPDAG or knowledge state on n<=5 was skipped or sampled.
- Any violation of L2, L3, L4 or AE-B is the single most valuable output of this run: it is a counterexample to the repair argument, not a bug to be smoothed over. It is persisted in full in the corresponding JSON file, with the CPDAG, `G`, both DAGs, the attempted path, and (for L3) exactly where it broke.
- L1 violations would indicate a bug in `covered_edges`/`reverse`/`dag_to_cpdag`, not a mathematical discovery -- L1 is the textbook fact the whole repair leans on.