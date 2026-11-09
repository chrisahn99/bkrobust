# E28: joint radius over orientation retractions and edge additions near the query

For the 984 scored E22 rows (learned PC and GES CPDAGs, A_ORACLE_R knowledge; 534 invalid at the truth, 452 of them with a true edge missing inside An_true({x, y} u Z)), rebuilds G, K and the committed set Z with E5's functions (exit 2 on any mismatch), then computes the joint radius: the fewest moves, each the retraction of one orientation of K or the addition of one directed edge between non-adjacent nodes within two hops of {x, y} u Z, after which Z fails the GAC. Depth 1 is exhaustive and lists every breaking move; depth 2 runs only when depth 1 has none. Beside it, the orientation-only radius (the ledger's `breaking_depths`) per row, and whether a true missing edge is among the distance-1 witnesses. Plan and registered statistics: `PREREG.md`.

Needs the finished E22 run (`results/rebuttal/E22/rows.shard0of1.jsonl`, `run.done`) and `experiments/rebuttal/_data/pest/`. No discovery is run.

    .venv/bin/python experiments/rebuttal/E28_joint_radius_edges/run.py --smoke && \
    .venv/bin/python experiments/rebuttal/E28_joint_radius_edges/analyse.py --out results/rebuttal/E28_smoke --n-boot 200
    bash experiments/rebuttal/E28_joint_radius_edges/launch.sh 2,3,4,5,6,7     # six shards, one core each
    # Slurm: run.py takes --shard I --n-shards K like the other experiments

Duration: smoke 6 rows in 34 s on one core (median 2.3 s per row, at most 15 s on arth150). Full run: 984 rows, a few core-hours, so under an hour on 6 cores unless depth 2 is needed often (600 s cap per row).
Outputs in `results/rebuttal/E28/`: `rows.shardIof6.jsonl` (one line per row: r_orient, r_joint, witness lists, coverage fields; resumable), `config.*.json`, `run.shard*.log`, `shard*.done`, `analyse.txt`, `results.json`, `run.done`.
Fidelity: `python ../fidelity_gate.py FIDELITY.yaml --no-judge` after the run.
`vendor_bkrobust.py` and `vendor_e5.py` are E22's copies, unchanged (header names their sha256).
