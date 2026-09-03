# S2 contract — space construction (written by the orchestrator, implemented by S2)

Module: `src/bkrobust/demo/space.py`  Tests: `tests/demo/test_space.py`

Everything is brute force. No conjectured characterisation of the covering
relation may be used: the covering relation is computed from model inclusion by
enumerating represented DAGs, exactly as the brief requires.

## Definitions

- The **space** `𝔊_Ĉ` is the set of all valid MPDAGs obtained by orienting some
  subset of `Ĉ`'s undirected edges and closing under Meek. Formally: all `G`
  that are valid MPDAGs (`meek.is_valid_mpdag`) AND satisfy `[G] ⊆ [Ĉ]`.

  > CORRECTION (orchestrator, after S2 reported). The first version of this
  > spec defined the space as "same skeleton as `Ĉ`, contains every directed
  > edge of `Ĉ`, and is a valid MPDAG". That is UNSOUND: `is_valid_mpdag(G)`
  > checks `G` against itself, so on the chain `A-B-C` the collider `A->B<-C`
  > passes it trivially -- a full DAG is vacuously Meek-closed and extends
  > itself -- while lying in a DIFFERENT Markov equivalence class. The model
  > inclusion condition `[G] ⊆ [Ĉ]` is required and is what the covering
  > relation is defined by anyway. S2 caught this and fixed it.
- `[G]` is the set of consistent DAG extensions of `G`
  (`meek.enumerate_dag_extensions`).
- **Model inclusion**: `G1 ⪯ G2` iff `[G1] ⊆ [G2]`. More orientations means a
  smaller model, so `Ĉ` itself is the unique maximum.
- **Covering**: `G1 ⋖ G2` iff `G1 ≺ G2` and no `H` in the space has
  `G1 ≺ H ≺ G2`.
- The **neighbour graph** is the undirected graph whose edges are exactly the
  covering pairs. Distance is BFS hop count on it.

## Required API

```python
enumerate_space(cpdag: MPDAG) -> list[MPDAG]
    # deterministic, sorted by (n_directed_edges, edge_string)

represented_dags(space: list[MPDAG]) -> dict[MPDAG, frozenset[MPDAG]]
    # [G] for each G, computed once and cached -- everything else reads this

covering_pairs(space, reps) -> set[tuple[MPDAG, MPDAG]]
    # (lower, upper) with lower ≺ upper, strictly-between check by brute force

neighbour_graph(space, covers) -> dict[MPDAG, set[MPDAG]]
    # symmetric adjacency

bfs_distances(nbrs, source) -> dict[MPDAG, int]
    # unreachable elements omitted from the dict, NOT recorded as infinity

all_pairs_distances(nbrs) -> dict[tuple[MPDAG, MPDAG], int]

check_metric_axioms(dists, space) -> dict[str, object]
    # keys: identity_of_indiscernibles (bool), symmetry (bool),
    #       triangle_inequality (bool), n_triples_checked (int),
    #       violations (list of readable strings, empty if none)
    # Check ALL triples. State the count in the return.

atomic_moves(g_from: MPDAG, g_to: MPDAG) -> list[str]
    # human-readable move sequence for adjacent-or-not pairs:
    # "orient A->B" (an assertion) / "un-orient A-B" (a retraction).
    # A FLIP must render as TWO moves: un-orient then orient the other way.
```

## Tests S2 must write

- On a 3-node chain CPDAG `A-B-C`: the space is exactly 6 elements.

  > CORRECTION (orchestrator). This spec originally said 4 (the CPDAG plus the
  > 3 DAGs). That was wrong. There are two additional non-maximal elements,
  > `B->A` with `B-C` still undirected and `B->C` with `A-B` still undirected.
  > Meek's R1 only propagates from an edge oriented INTO the shared node, so
  > orienting AWAY from B forces nothing and leaves a genuine partially
  > oriented element. Verified independently by the orchestrator.
- On a 4-clique CPDAG: the space's maximal elements (the DAGs) number `4! = 24`,
  since a clique's v-structure-free acyclic orientations are its linear orders.
- `Ĉ` is the unique maximum: `[Ĉ] ⊇ [G]` for every `G`.
- Every DAG in the space is minimal (its `[G]` is a singleton).
- Covering pairs differ by at least one orientation and nothing lies strictly
  between (re-verify by direct search, independently of how they were built).
- Metric axioms hold over ALL triples; report how many were checked.
- `bfs_distances(nbrs, cpdag)[cpdag] == 0`.
- A flip yields exactly 2 moves from `atomic_moves`; a single assertion yields 1.
- Cross-check `enumerate_space` against an independent brute-force enumeration
  over all `3^k` orientation-state assignments to the `k` undirected edges,
  filtering by `is_valid_mpdag`. The two sets must be identical. This is the
  brief's "every valid MPDAG refining Ĉ is enumerated" check.
