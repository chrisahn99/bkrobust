# E22: which discovery error breaks each committed set (prereg, 2026-11-06)

Written 2026-11-06, before any run. Not edited after the first run; changes go in `DEVIATIONS.md`.

**Question.** For every committed set that is invalid in the generating DAG, which discovery errors does the learned
class contain near the query, and which adjustment condition fails at the truth? Serves reviewer points C2Bd-W2,
h8pC-W2 and C2Bd-Q3. No discovery is run: the learned graphs and rows are those stored by E5.

**Population.** Every row of `results/rebuttal/E5/rows.shard0of1.jsonl` with phase B, `graph_kind` learned, arm
A_ORACLE_R and status ok. At planning time this is 987 rows on 178 learned graphs: 534 invalid at the truth and 453
valid (the plan said 593 valid; the stored file holds 453). Networks ecoli70, arth150, magic-niab, magic-irri.

**Reconstruction check (exit non-zero, run stops).** For each learned graph, the A_ORACLE_R knowledge is rebuilt with
E5's functions (`truthful` on the sorted undirected edges, then `greedy_consistent`), G0 = `apply_orientations`, and
Z = `committed_adjustment_set(G0, x, y)`. Every row must reproduce the stored `k_sha256`, `n_k`, `committed_verdict`,
`z_size`, `y_true_descendant` and `z_valid_at_truth`. One mismatch stops the run with exit 2 and the row is reported.

**Census (learned CPDAG against the true CPDAG, net difference).** Error types: missing adjacency, extra adjacency,
false v-structure (in the learned graph, not in the true CPDAG), missed v-structure (the reverse), reversed compelled
edge (directed both ways in the two graphs), under-oriented and over-oriented edges (directed in one, undirected in the
other), and, per row, query direction wrong (y is a possible descendant of x in G0 but not a descendant in the true
DAG). Each error records its nodes, its hop distance from {x, y} in the true skeleton (minimum over its nodes; null when
unreachable), and whether all its nodes lie in A = An_true({x, y} ∪ Z) (ancestors taken inclusively). The full list is
written for valid and invalid rows alike.

**Failing condition at the truth (GAC).** Checked in this order, the first that holds is the first failing condition:
(1) y is not a descendant of x in the true DAG; (2) Z contains a forbidden node; (3) a non-causal path is open given Z
in the proper back-door graph. All three flags are recorded. A set is invalid exactly when (2) or (3) holds; that must
equal the stored validity (consistency check).

**Registered statistic.** The number of invalid rows whose differences from the truth inside A are only extra
adjacencies with orientations that agree with the truth. Operationally, inside A: no true edge is missing from G0; no
directed edge of G0 contradicts the true orientation; and imposing the true orientation on every undirected G0 edge
that is a true edge succeeds under the Meek closure. This is the hypothesis of statement S3 restricted to A, which is
where validity in the true DAG is decided (every causal path, every forbidden node in Z and every d-connecting path
lies in A). Under it some DAG D in [G0] contains the true DAG on A with the same orientations, Z is valid in D, and S3
gives Z valid in the true DAG. Reported with its share of invalid rows and a 95 % percentile interval from a bootstrap
that resamples learned graphs (2000 draws, seed 20261106). For each row passing the hypothesis, the run also builds one
such D and checks that Z is valid in it and that D contains the true DAG on A.

**Prediction.** 0 rows.

**Against us.** One or more such rows. Then S3, or the code, disagrees with the data; we say so in the answer, and
find which by reading the rows. Also against us: more than 20 % of invalid rows have no census error of any type within
two hops of x or y. Then the failure is not local, and the answer must not say that discovery errors near the query
explain the invalid sets.

**Reported with it, not registered.** The share of invalid rows by first failing condition and by nearest error type
(ties broken by a fixed order), with the same graph-clustered intervals. The literal variant of the statistic on the
census (invalid rows whose census errors inside A are all extra adjacencies). Per error type, the share of valid and of
invalid rows with that error within two hops. Each split by network.

**Seeds, stopping.** Deterministic: no random draw except the bootstrap (seed 20261106). One launch of the full set,
resumed only after a crash. No rerun with changed settings; if a consistency check fails, fix, rerun and report both.
