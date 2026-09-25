# X-component mean-bias profile on the 18 larger instances (v2)

Script: `experiments/mean_table_xcomp_extend_v2.py`. Enumerates only the states
that agree with G0 outside X's undirected component (distance by the closed
form |K_G Δ K_H|), instead of the whole space. Validation: reproduces all 39
queries of `../x_component/` exactly (`validation.json`: all_match = true).

Of the 18 instances skipped by the whole-space run: 14 computed (900 s budget),
4 timed out (ecoli70 cspG→hupB, paths 1→10, 1→8, 11→14; X's component has 13-14
undirected edges; cost is DAG-extension enumeration inside `bias_at`). Of the
14: 4 have no failure (r = inf), 10 have finite r_val and nonzero bias and are
added to the paper's table (49 rows in total).

Enlarged set (49): non-monotone 9; non-persistent cells 16 of 271 finite, on 5
instances; distinct r_mu signatures 22, with B(Chat) 32; B(Chat) alone 21;
(r_val, B(Chat)) 22; Pearson(#distinct finite r_mu, D) = 0.75; D-bin means
1.00 (D<=3, n=21), 1.25 (4-6, n=12), 2.82 (7-9, n=11), 3.60 (>=10, n=5).
