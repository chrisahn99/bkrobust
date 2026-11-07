# E21: compute the walkthrough toy exactly before drawing it

Computes, with the submitted code (GAC), every number the walkthrough figure will show for the
four-node toy (CPDAG U-X, X-Y, X-W, W-Y; K = U->X, X->Y, X->W), the radius of the Figure 2(a)
toy, and the running example (scenario B): Z = O(G0), |space|, the up-set of G0 with validity,
nearest failing states, DAGs below G0, non-states, claim radius. Each radius is computed six ways
and must agree (exit 2 otherwise). Predictions and deviations: PREREG.md.

    .venv/bin/python experiments/rebuttal/E21_walkthrough_toy/run.py --smoke --out results/rebuttal/E21_smoke
    bash experiments/rebuttal/E21_walkthrough_toy/launch.sh

Duration: seconds (smoke < 30 s, full well under 5 min on one core).
Outputs in results/rebuttal/E21/: rows.jsonl (one row per toy, with the full space and covers
for drawing), config.json, results.json and the printed table from analyse.py, run.log, DONE.
