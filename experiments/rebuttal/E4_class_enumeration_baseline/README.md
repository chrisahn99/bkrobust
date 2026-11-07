# E4: listing the class as a timing baseline

Measures for which class sizes |[C^]| computing the radius by Corollary 1(c) (list every DAG of the class, take the
fewest reversed G0 orientations over DAGs where Z fails) is faster than `bkrobust.hybrid.breakdown_radius`, on the
four-node, highk (144) and timing (256) instances, under GAC. Every complete listing must return the hybrid's radius.
Registered statistic and prediction: `PREREG.md`. Run checklist: `FIDELITY.yaml`.

    .venv/bin/python experiments/rebuttal/E4_class_enumeration_baseline/run.py --smoke          # < 60 s, results/rebuttal/E4_smoke
    bash experiments/rebuttal/E4_class_enumeration_baseline/launch.sh 6              # full run, nohup
    .venv/bin/python experiments/rebuttal/E4_class_enumeration_baseline/analyse.py results/rebuttal/E4

Expected duration: 2 to 6 h on 6 workers (dominated by capped listings, at most 2 x 120 s per highk/timing instance).
Outputs in results/rebuttal/E4: rows_shard*.jsonl (one row per instance), manifest_shard*.json, logs/,
analyse.txt, results.json, E4.done. Array jobs: --shard I --n-shards K, then analyse the shared directory.
