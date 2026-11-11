# Rebuttal experiments (ICLR 2027, submission 46323)

This folder holds the experiments run during the discussion period, E1 to E29. Each one answers one
reviewer point. Every experiment folder is self-contained: a pre-registration written before the run
(`PREREG.md`), the code (`run.py`, `analyse.py`, and `vendor_*.py` copies of older helpers with their
source commit in the header), a launcher (`launch.sh`), a fidelity rubric (`FIDELITY.yaml`), a short
`README.md` with the exact smoke and full commands, and, when the run did not go as planned,
`DEVIATIONS.md`. The code imports the submitted package from `src/` of this repository; nothing here
changes `src/`.

Outputs are written to `results/rebuttal/<E>/` (or `results/rebuttal/<E>_smoke/`), which git ignores.

## What each experiment measures

| id | folder | what it measures |
|---|---|---|
| E1 | `E1_cpdag_check_on_lower_bound_rows` | how many of the 372 released lower-bound radii become exact (infinite) once the committed set is checked at the CPDAG |
| E2 | `E2_running_time_by_graph_size` | running time of the submitted hybrid under the GAC, from 8 to 400 nodes |
| E3 | `E3_depth_budget_and_sat_mode` | effect of the depth budget b and of the SAT mode on the 144 high-k instances |
| E4 | `E4_class_enumeration_baseline` | listing the Markov equivalence class as a timing baseline against the hybrid, with agreement of the radii |
| E5 | `E5_learned_cpdags_by_sample_size` | validity at the true graph and the radius when the CPDAG is learned from data, by sample size |
| E6 | `E6_bootstrap_stability_check` | bootstrap stability of the committed set as a companion check for discovery error |
| E7 | `E7_one_hidden_variable` | one hidden variable, placed near or far from the query |
| E8 | `E8_missed_and_extra_edges` | how local the damage of a single missed or extra edge is |
| E9 | `E9_parametric_bootstrap_validity` | validity of the committed set under a parametric bootstrap from a fitted linear Gaussian model |
| E10 | `E10_radius_across_valid_sets` | the radius across all valid adjustment sets of a query |
| E11 | `E11_survival_on_cpdag_only_queries` | survival ranking on queries chosen from the CPDAG alone |
| E12 | `E12_networks_above_500_nodes` | radius and audit on the five networks above 500 nodes |
| E13 | `E13_committed_set_consistency_checks` | consistency checks on committed sets and the error budget |
| E14 | `E14_data_test_for_invalid_sets` | whether a data test across sets valid in G0 flags an invalid committed set |
| E15 | `E15_aid_as_survival_baseline` | the adjustment identification distance as a survival baseline |
| E16 | `E16_radius_where_separation_undefined` | radius ranking on queries where separation is undefined |
| E17 | `E17_worked_case_and_validity_by_radius` | a worked applied case, and validity at the true graph by radius |
| E18 | `E18_section7_standardised_models` | the Section 7 bias profiles under standardised models, ten coefficient draws |
| E19 | `E19_bias_radius_internal_standardisation` | the Appendix E.4 bias-radius study under internal standardisation |
| E20 | `E20_uumc_generator` | the Section 7 bias profiles under the UUMC generator |
| E21 | `E21_walkthrough_toy` | exact enumeration of the walkthrough toy example drawn in the paper |
| E22 | `E22_e22_discovery_error_census` | which discovery error breaks each committed set |
| E23 | `E23_e23_pc_ges_disagreement_check` | PC and GES agreement as a check for discovery error |
| E24 | `E24_e24_far_error_transplant` | real discovery errors far from the query, transplanted into the true graph |
| E25 | `E25_amenability_check` | whether the nearest failing states are exactly the non-amenable states |
| E26 | `E26_taeb_edge_addition_tolerance` | an edge-addition tolerance in the sense of Taeb et al., compared with the radius |
| E27 | `E27_amenability_identity_new_populations` | the radius versus the amenability radius on query populations not used before |
| E28 | `E28_joint_radius_edges` | a joint radius over orientation retractions and edge additions near the query, against the discovery errors of the learned graphs |
| E29 | `E29_radius_vs_degree_sep_undefined` | the radius against the undirected degree of the treatment where the separation score is undefined |

## Setting up

From the repository root, with Python 3.11 or later:

    python -m venv .venv
    .venv/bin/pip install -e .
    .venv/bin/pip install pyyaml          # only for fidelity_gate.py

A few experiments list extra packages in their own `requirements.txt` (empty when there are none).
The launchers expect the environment at `.venv/bin/python` and are run from the repository root.
Single-threaded BLAS is assumed: `export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1`.

## Running one experiment: smoke first, then the full run

Each experiment `README.md` gives the two commands. The pattern is:

    # 1. smoke: a few units, small budgets, usually under a minute
    .venv/bin/python experiments/rebuttal/E1_cpdag_check_on_lower_bound_rows/run.py --smoke
    .venv/bin/python experiments/rebuttal/E1_cpdag_check_on_lower_bound_rows/analyse.py --out results/rebuttal/E1_smoke

    # 2. full run, detached, with a worker count; writes run.log, run.pid and a DONE file
    bash experiments/rebuttal/E1_cpdag_check_on_lower_bound_rows/launch.sh 8

The launchers refuse to start when a run is already live or finished in the output folder, and never
delete anything. Runs are resumable: a restarted run skips the units already written.

`_smoke/run_all_smoke.sh [E1 E2 ...]` runs the smoke test of every (or of the listed) experiments one
after the other with a ten-minute limit each, and writes a summary to `_smoke/SUMMARY.txt`.

### The two machines

The runs used two machines, and either one is enough for any experiment:

* **a 16-core workstation** (15 GB of memory). The launchers cap the worker count (8 for most
  experiments) to stay within memory. `_queue/run_lane.sh <lane> <cpu-list> <workers> E5 E6 ...`
  runs several experiments one after the other on a fixed set of cores, waiting for each DONE file.
* **a Slurm cluster**, for the long experiments: every `run.py` except E25's accepts
  `--shard I --n-shards K`. `_slurm/job_smoke.slurm` runs one smoke test on a compute node and
  `_slurm/job_array.slurm` runs one array, for example

      mkdir -p logs
      sbatch --job-name=rb_E26 --array=0-7 --export=ALL,EXP=E26,NSHARDS=8,EXTRA="" experiments/rebuttal/_slurm/job_array.slurm

  Set `BKR_ROOT` to the checkout on the cluster (default: the submit directory) and `BKR_VENV` to a
  Python environment (default: `$BKR_ROOT/venv.tar.gz`, unpacked on the node). Partition, time and
  memory in the two files are examples; adjust them to your cluster. The analysis step
  (`analyse.py --out results/rebuttal/<E>`) runs afterwards on any machine.

## Data files

The released inputs used by most experiments are committed in `_data/` (sha256):

| file | sha256 |
|---|---|
| `_data/per_query_records.csv` | `13541300e559b7aae464d35121b1f5f8a9678209dac99063ba1e31768412452e` |
| `_data/knowledge_states.json` | `45119393225fee57ed1dd4c002130f805d6a655c2ee7648006470918c2dfc7bb` |
| `_data/chris_e10/amenability_radius.csv` | `55b19b1e6510d8ef2caeb2df39b2c3eb11547f78523de90808a687c105c4dd7e` |
| `_data/pest/chat_e.json` | `7f2d22ca4b778f25376d4071ed3f6ddae2e2a085f308a1e8940f21637318e30a` |
| `_data/pest/frame.jsonl` | `530bbdd0619b0e23606a8513ae5d643ce6ee9cc471d9add43360632a3609a530` |
| `_data/pest/graphs.jsonl` | `708871335eff33eac4db09754ab9c561e8fa7ba552cbce0b8cf679375115b402` |
| `_data/pest/graphs_retry.jsonl` | `71f2bffada8d467294766eb26c5f6e472fcdd4973aa4ead750285cc5d1f2b093` |
| `_data/pest/knowledge.json` | `fbe8741137df8b670114e72b6fdc6c2ac9e79d11a1cc83bc48550b021604040c` |
| `_data/pest/rows.csv` | `0763a714f002e3698351ff6fd39770aace23f1e983a2f6b21d808338c94bcd7c` |
| `_data/pest/structure.csv` | `25c3f88099c54e777a642e6e100d79b93d2bd93d5b0dc7e6255ed9f319f89c41` |

`per_query_records.csv` and `knowledge_states.json` are the paper's released files.
`chris_e10/amenability_radius.csv` is the first author's amenability-radius table on the released
queries. `pest/` is the learned-CPDAG panel (bkrobust main branch, commit ba7ddc12, `results/pest/`).
`E7_one_hidden_variable/panel_samples/` holds the two panel samples that do not regenerate bit for bit
on x86_64 (see that folder's `DEVIATIONS.md`). Networks come from `results/axisa3/networks/` of this
repository. Check the files with `shasum -a 256 -c` against the table, and many `run.py` scripts
also record or verify the hashes of their inputs.

## Pre-registration

Each `PREREG.md` was written and frozen before the first full run of its experiment. It states the
question, the population, the registered statistic, the decision rule, and what result would count
against the paper. `_queue/PREREG_SHA256.txt` holds the hash of each pre-registration and fidelity
rubric at the time it was frozen. Anything changed afterwards is written in that experiment's
`DEVIATIONS.md`, with the cause and the alternative not taken. Some runs also record the hash of
`PREREG.md` in their config, and their fidelity rubric checks that the run started under that text.
`_queue/REDACTION.md` lists the files whose text changed only to remove local paths and machine
names before sharing, with the registered and the current hash.

## Fidelity gate

A run with plausible numbers can still be the wrong experiment. Each `FIDELITY.yaml` lists, before the
run, what a run must have done to count as that experiment: values read from the run's saved config
(`config` items), read-only shell checks run from the repository root (`check` items), and lines that
a person decides by reading the listed code (`judge` items). After a full run:

    python experiments/rebuttal/fidelity_gate.py experiments/rebuttal/E1_cpdag_check_on_lower_bound_rows/FIDELITY.yaml

prints one verdict per item and writes `FIDELITY.verdict.json` beside the rubric (ignored by git).
Exit 0 means every decided item passes, 1 a failure, 2 something unclear. The verdict is decided
apart from the result: a run that fails its rubric is reported as not being that experiment,
whatever it measured.
