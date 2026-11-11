# E20: Section 7 profiles under the UUMC generator

Recomputes the mean-bias profiles r_mu and the ceiling B(Chat) of the 49 Table 6 queries under three coefficient generators, with ten draws each: the naive submitted random_sem, its exact internal standardisation (iSCM), and UUMC (Herman et al. 2025, Algorithm 1). The registered statistic is the number of distinct r_mu signatures per draw. The run also counts signatures added beyond (radius, ceiling), checks the Table 2 contrast, and reports pair-counted var- and R2-sortability. See PREREG.md.

    PYTHONPATH=src .venv/bin/python experiments/rebuttal/E20_uumc_generator/run.py --smoke --workers 2 && .venv/bin/python experiments/rebuttal/E20_uumc_generator/analyse.py --out results/rebuttal/E20_smoke
    bash experiments/rebuttal/E20_uumc_generator/launch.sh 6

The smoke run covers 4 queries, 3 arms and 2 draws, in well under a minute. The full run covers 49 x 3 x 10 = 1470 units. Its cost is the per-query graph cache, done once: about 10 min of CPU, mostly child CO2->DuctFlow at about 5 min and the four Kampen_2014 queries at about 1 min each. The estimate is 15 to 40 min with 6 workers; it has not been measured.

Outputs, in results/rebuttal/E20/:
- units.shard*.jsonl: one row per unit, with criterion and commit;
- config.shard*.json, run.log and DONE;
- results.json, the registered statistic and the verdict.

Both scripts exit 2 on any consistency failure.
