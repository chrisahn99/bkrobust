# E3: depth budget b and SAT mode on the 144 high-k instances (reviewer qQYf, W1)

Measures, on the instances of `results/axisb3/highk.jsonl` (GAC, CPDAG check on), the total time of the radius computation for b in {0,1,2,3,4,6} x SAT mode {ladder, one CP-SAT minimise call}, and checks that all 12 arms and the library `breakdown_radius` return the same radius. Registered statistics and deviations: `PREREG.md`.

Commands, from the repository root:

    .venv/bin/python experiments/rebuttal/E3_depth_budget_and_sat_mode/run.py --smoke && .venv/bin/python experiments/rebuttal/E3_depth_budget_and_sat_mode/analyse.py --smoke
    bash experiments/rebuttal/E3_depth_budget_and_sat_mode/launch.sh   # then analyse.py without --smoke

Smoke: 3 instances at n = 8, all arms, plus a cap probe; under a minute. Full: serial, expected 1 to 3 h (App. D summed E1 time over the file is 134 s; the b = 6 search and the 13 calls per instance dominate).

Outputs in `results/rebuttal/E3/{smoke,full}/`: `rows_shard*.jsonl` (one line per instance), `config_shard*.json` (read by `FIDELITY.yaml`), `results.json` and the printed table from `analyse.py`; the full run also writes `run.log`, `run.pid` and `DONE` (holds the exit code; 2 = a consistency failure). Delete `DONE` only to resume a crashed run.
