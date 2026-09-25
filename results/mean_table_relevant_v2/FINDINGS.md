# Mean-bias profile restricted to the query's relevant components (v2)

Script: `experiments/mean_table_relevant_v2.py` (reuses the per-instance objects of
`experiments/mean_table_fullspace.py`: same panel, seeds, SEM and space, so B(G)
values are identical; only the averaging set changes).

Two rules for which undirected chain components of the CPDAG count as relevant:

- `conservative` (top-level files here): components containing X, a node of
  Z ∪ {Y}, or a possible ancestor of X or Y. 6 of 39 queries change.
- `x_component` (`x_component/`, used in the v2 paper): only X's own component.
  Justification: B(G) = max_D |theta_Z − beta(Pa_D(X))| depends only on Pa_D(X),
  and chain components of a CPDAG are oriented independently. Checked on the data:
  `n_invariance_violations = 0` on all 39 queries (B constant across states with the
  same restriction to X's component). 16 of 39 queries have other components, and on
  all 16 the whole-space average delays the r_mu onsets.

Distinct signatures over the 39 queries (x_component): r_mu profile 17 (whole space
23); profile with B(Chat) 24 (whole space 29); B(Chat) alone 16; (r_val, B(Chat)) 17.
