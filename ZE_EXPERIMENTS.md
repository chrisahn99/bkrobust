# Zé's experiments on the breakdown radius, 17 July to 16 September 2026

Everything I ran on our line, on one branch (`ze/experiments`), on top of your history so that it
merges. The explanation, in order, with every number tied to a file here, is the polycopié:
**`rho_breakdown/polycopie_2026-09-16/polycopie.pdf`**. This file is the map.

Abstract and author lock **Friday 18 September (AoE)**, paper **Friday 25 September (AoE)**.

## What is on the branch

```
iclr2027-eval-plan (ba7ddc1, already on GitHub)       evaluation plan, protocol, stage-0 claims
+ merge of rho-breakdown-evidence (never pushed)       the August package, rho_breakdown/
+ commit 1  src/, tests/, .gitignore                   library changes of the protocol run
+ commit 2  experiments/, docs/                        protocol-run scripts and docs
+ commit 3  results/                                   protocol-run results and pre-registrations
+ commit 4  rho_breakdown/<dated folders>              the experiment folders of 14 Aug to 16 Sep
+ commit 5  this file + rho_breakdown/polycopie_2026-09-16/
```

`git diff d3cd71f ze/experiments --stat` shows exactly what I added to your code. Your `main` has
moved since `d3cd71f`; merge rather than rebase.

## 1. The protocol run of 11 to 13 September (inside your package layout)

Read `docs/PROTOCOL_RUN.md` first: it lists the stages in the order they ran (stage 0 the ladder
and the unit, 0B prior art, 1 the committed set in polynomial time, 2 the frame, 3 the questionnaire,
4 the suppliers, 5 the sweep, 6 the estimated-CPDAG panel, 7 the audit) and every departure from the
protocol as written, with its reason.

| folder | what it holds | headline file |
|---|---|---|
| `results/stage0/` | the free-rule ladder, the claim-unit radius on the committed rows, the corrupted-knowledge calibration, the Lever 0 differential | `STAGE0.md`, `PREREGISTRATION.md` |
| `results/frame/` | the observables-only frame (80,132 candidate pairs, 33 networks), witnesses, the declared queries; `status_*` is the first status run, superseded by the ledger's statuses (`results/ledger/sweep_summary.json`) | `frame_summary.json`, `per_network.csv` |
| `results/elicit/` | the frozen, hashed questionnaire (540 items), every supplier's answers (7B to 72B), the knowledge they became, the source-overlap probe, the ancestral block | `PREREGISTRATION_SCALE.md`, `SOURCE_OVERLAP.md` |
| `results/ledger/` | the certificate ledger (8,308 rows, 14 knowledge conditions), Table 4, the checks, the funnel, the decision table | `TABLE4.md`, `CHECKS.md`, `PREREGISTRATION.md` |
| `results/pest/` | the estimated-CPDAG panel (PC, GES; Table 7) | `TABLE7.md`, `PREREGISTRATION.md` |
| `results/interval/` | M1 and M1b: the knowledge interval, the null radius, bands, Table 1 (L95), Table 2 (sensitivity reports) | `TABLE1_L95.md`, `TABLE2_REPORTS.md`, `PREREGISTRATION.md` |
| `results/substrate/` | sortability probes, the per-network structure table, fixture decisions | `TABLE2_SUBSTRATE.md`, `DECISIONS.md` |

New library code: `src/bkrobust/estimation/knowledge_interval.py`, `knowledge/ancestral.py`,
`mpdag_criterion/optimal.py` (closed-form O*), `benchmarks/discovery.py` (PC/GES arm),
`benchmarks/audit.py` (the generalised adjustment criterion re-implemented in networkx, sharing no code
with the instrument). New tests: `tests/benchmarks/test_{audit_oracle,discovery,measure_selection,sortability_probe}.py`,
`tests/criterion/test_optimal_closed_form.py`, `tests/estimation/`, `tests/knowledge/`: **39 passed**
(16 September, `pytest -q` on those paths).

Large models ran on Jean Zay with vLLM: `experiments/remote/elicit_vllm.py`, `elicit_jz.slurm`
(account line `#SBATCH --account=zkh@h100` is my project; change it). Their raw answers are in
`results/elicit/remote/*.jsonl` and enter the cache through `experiments/elicit_import.py`;
`elicit_run.py` refuses to call a cluster model directly.

### Caveats about these files

- **Stage-0 calibration**: `results/stage0/calibration_{rows.csv,summary.json}` are the rerun with
  SHA-256 seeds (8,911 scored of 16,496 draws). `STAGE0.md` and `PREREGISTRATION.md` still quote the
  first run (8,553 rows), which seeded with Python's salted `hash()` and is not reproducible across
  processes. The rerun changes one reading: the hop radius goes from 0 to 3 dangerous rows on the
  corrupted rows (all at corruption rate 0.5), the free rule from 8 to 9, the claim radius stays at 0.
  Salted `hash()` seeded six things in all: the chance arm of the ledger sweep, the Lever 0 differential,
  the calibration, the first frame-status run, the questionnaire builder and the source-overlap probe
  (`docs/PROTOCOL_RUN.md` heads the list «Four experiment scripts» and names five). All now use SHA-256
  seeds; the sweep, the differential, the calibration and the probe were rerun. The frozen questionnaire does not
  regenerate byte for byte and stays the authoritative input; its builder now refuses to overwrite it.
- **T3 selection** (`experiments/interval_tables.py`): the printed T3 takes five rows per network
  sorted with a finite null radius first, which over-represents exactly the outcome the table
  illustrates (76 % finite among printed rows, 43 % among the 37 eligible). Fix: drop that sort key,
  or print all 37. See `rho_breakdown/panelB_T3_2026-09-16/`.
- **Validity oracle**: `is_valid_adjustment_set_dag` is Pearl's back-door criterion; the certificate
  uses the generalised adjustment criterion. Audits must use `bkrobust.benchmarks.audit`.
  `tests/benchmarks/test_audit_oracle.py` guards it. Two places still differ and are not yet changed:
  `experiments/declared_queries.py` scores `z_valid_at_truth` with the back-door function (the committed
  set judged at the truth is not the truth's optimal set, so the criteria can disagree, as they did on 4
  of 1,169 ledger rows; rerun it with the GAC audit before quoting that column), and the hop
  radius, through `hybrid.py`, tests validity with `mpdag_criterion.is_valid_mpdag` while the claim
  radius uses `gac.is_gac_valid_mpdag`. The stricter test can only make `r_hop` smaller, so the hop
  unit's over-certification counts are, if anything, conservative.
- **August `rho_breakdown/NUMBERS.md`, two corrections found on 16 September**: the row
  «0.554 → 0.515, censoring of `rho*_any`» is the `rho*_se` series (`e1prime/RESULTS.md` table,
  n = 200 → ∞); `rho*_any` goes 0.183 → 0.166. And «13-27×» separation by hop distance holds for the
  `original` (26.7×) and `licensed` (13.4×) ensembles only; on `large` it is about 3×.
- **`rho_breakdown/eval_protocol_2026-09-11/02_PROTOCOL.md`** was written before the code diagnosis of
  the same day and still says G_0 has no undirected edge on all 831 rows and makes Lever 0 the blocking
  step; both were withdrawn in `06_CODE_FIXES.md` (288 rows have a free edge; Lever 0 recovers 0 of 518).
  `stage0/rclaim_cov1.py` writes to a hard-coded scratch directory; edit `SCR` before running it.
- Not on the branch on purpose: `figures/*` and `results/axisa3/networks/acquisition_manifest.json`
  (regenerated locally by a reproduction run), the pgmpy sdist tarball.

## 2. `rho_breakdown/`: the experiment folders, in date order

Everything here is CPU-only unless stated. Most scripts keep the absolute paths of the machine they
ran on. The mapping is:

| path inside scripts | here |
|---|---|
| `/Users/josecosta/bkrobust/` | the repository root |
| `/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/` | `rho_breakdown/locality/` |
| `/Users/josecosta/mugango/output/2026-08-19_e1prime-se-criterion/` | `rho_breakdown/e1prime/` |
| `/Users/josecosta/mugango/output/2026-08-31_x3-skeleton-perturbation/` | `rho_breakdown/x3_skeleton/` |
| `/Users/josecosta/mugango/output/2026-09-07_radius-redef/` | `rho_breakdown/radius_redef_2026-09-07/` |
| `/Users/josecosta/research-pilots/latent-causal-rho-breakdown-knowledge/` | `rho_breakdown/pilot/` |
| `/home/costaj/latent-causal/…` | the Linux workstation some August sweeps ran on; the outputs the reports use are in each folder's `results/` |

| folder | date | question | read |
|---|---|---|---|
| `pilot/` | 17 Jul | does knowledge-induced breakdown exist at all on SCMs? 600 + 1,650 + 60 SCMs, iSCM, closed-form linear SEMs | `RESULTS.md` |
| `rerank_2026-08-14/` | 14 Aug | three experiments that decided the direction: e2 is the radius finite, e3 is O* relevant, e4 what Guo and Perković actually cover (the paper PDFs are not included, see their arXiv ids) | `RESULTS.md`, `e*/RESULTS.md` |
| `e1prime/` | 19 Aug | the `se` criterion that removes the oracle statistic ρ*, the censoring floor, calibration and coverage | `PREREG.md`, `RESULTS.md`, `AUDIT.md` |
| `locality/x1/` | 19 Aug | hop-distance stratification, operator arms, pool composition; `audit/` holds the eligibility audit and the refuter | `PREREG.md`, `RESULTS.md`, `AUDIT.md` |
| `locality/x2/` | 19–23 Aug | exhaustive enumeration with no SCM and no truth (complete at p ≤ 5); `cfirewall/` the C-FIREWALL census; `chat/` the estimated-Ĉ arm; `firewall/`, `falsegain/`, `refutation/` the arms that tried to break it | `RESULTS.md`, `AUDIT.md`, `cfirewall/STATUS.md` |
| `locality/VERDICT.md` | 19 Aug | what survived of the locality law | |
| `x3_skeleton/` | 31 Aug | does the locality law survive a wrong skeleton, and one worked example end to end | `README.md`, `RESULTS.md` |
| `REPORT.md`, `NUMBERS.md`, `deck.pdf` | 31 Aug | the August report (every number mapped to its file, four self-corrections) | |
| `radius_redef_2026-09-07/` | 7–8 Sep | what should the radius measure: the candidate redefinitions, their audits (`audit_*.py`), the retraction hull, `certified-core/`, the null radius r_0; `runs_2026-09-08/` the five-node worked example, the stratum and locality runs and their JSON (`resultados.json`, `exemplo.json`) | `EXPERIMENTS.md`, `IMPROVEMENT.md` |
| `rnull_real_2026-09-09/` | 9 Sep | r_0 on four Gaussian networks with fitted coefficients, from Σ̂ | `README.md` |
| `eval_protocol_2026-09-11/` | 10–11 Sep | where to evaluate: the literature ledger, the protocol «Beat the BFS First», the structured judgements, the code diagnosis of the old corpus, the stage-0 probes | `02_PROTOCOL.md`, `06_CODE_FIXES.md`, `10_PROTOCOL_RUN_TABLE4.md` |
| `panelB_T3_2026-09-16/` | 16 Sep | what a Table 1-B cell measures, T3 in full, the one-claim cliff, the T3 selection audit; `make_figs.py` reproduces the 27 panel-A cells of `TABLE1_L95.md` and asserts it | `PANELB-T3.pdf`, `README.md` |
| `table1_dossier_2026-09-16/` | 16 Sep | thirteen Table 1 candidates, claims C1–C11, recommendation T11 | `TABLE1-DOSSIER.pdf`, `README.md` |
| `theory_checks_2026-09-16/` | 16 Sep | enumeration of the four-node classes used to check the draft's theory | the script |
| `polycopie_2026-09-16/` | 16 Sep | the explanation of all of the above, the review of the 16 September draft, and the merge | `polycopie.pdf` |

### The five notes written in Portuguese, in one paragraph each

- `radius_redef_2026-09-07/IMPROVEMENT.md` (8 Sep). Two rounds. (1) The null radius r_0, the smallest
  retraction budget at which the hull of possible effects contains zero, looked worse than r_val on
  all problems (84 % «never»), but only 21.2 % of sampled problems are informative (the class leaves
  some uncertainty about the effect). Conditioning on that stratum (1,800 problems, 381 informative),
  r_val collapses (89.5 % at 1, 1.51 effective values) and r_0 opens (3.46 effective values, 0.60 →
  1.79 bits). So r_0 becomes the main quantity, r_val is reported beside it, and the informative
  stratum is a declared scope condition. (2) Locality measured in coverage with a negative control, at
  a matched budget of one retraction: retracting claims that do not touch the query moves coverage by
  +0.2 pp with one error and +0.7 pp with two; retracting the ones that do recovers almost all of it
  (0.998 vs 1.000) at 0.621 vs 0.630 of blanket width. Locality pays in cost (28 % fewer enumerations)
  and in what to elicit, not in width: it is lemma L2 seen from the other side.
- `rnull_real_2026-09-09/README.md` (9 Sep). r_0 computed from Σ (population or sample) by a routine
  that never receives the true DAG, on the Gaussian networks with fitted coefficients. `load_real.py`
  loads a network as a `LinearSEM`, `rnull_real.py` computes the estimand, its standard error, the hull
  and the radii. Two qualifications the README omits: `arth150` hit the 20-minute wall at coverage 0.5,
  so the summary covers ecoli70, magic-niab and magic-irri; and the run used the oracle CPDAG with
  knowledge from `select_knowledge`, which reads the true DAG, so it tests the estimator, not a
  deployment. It runs, but only with a confidence band: the plug-in hull is anti-conservative.
- `locality/x2/chat/RESULTS.md` (23 Aug). The firewall with an estimated Ĉ (PC written from scratch,
  gated on recovering the CPDAG exactly under an oracle CI test: 150/150 and 60/60). Conditioning on
  the base report already being valid, the rate of newly invalid sets that pass the check falls from
  10.87 % at n = 200 to 2.44 % at n = 50,000, against 17–30 % for unchecked sets: the firewall's
  efficacy is a monotone function of discovery quality.
- `locality/x2/cfirewall/STATUS.md` (23 Aug). C-FIREWALL, not refuted: 0 violations in 70,375,325
  trials (census at p = 3, 4, 5; samples at p = 6, 7, 8), rule-of-three bound < 4.3 × 10⁻⁸, detector
  power 21.8 % invalid in the control arm. A census at p = 6 is intractable (≈ 160 core-hours) and
  mostly wasted, hence sampling, labelled as such.
- `eval_protocol_2026-09-11/06_CODE_FIXES.md` (11 Sep). Two of my earlier claims withdrawn:
  `g0_undirected_edges` is only assigned in the `O_INTRACTABLE` branch of `measure.py`, so it is not
  evidence that G_0 was fully oriented (543 of 831 were, 288 were not); and Lever 0 recovers nothing:
  all 518 `o_g0_not_identified` rejections are non-amenable, 0 amenable. The deletion is right, the
  label is wrong. Where the defect really is: `select_knowledge` thins the recovering set in
  alphabetical order with a constant step; on 40 non-degenerate cases in the eight smallest networks,
  29 (72 %) have another subset of the same size that keeps amenability. Nine fixes in three layers,
  the last being that `select_knowledge` reads the true DAG, which only elicited knowledge changes.

## 3. Reproducing

```bash
python -m venv .venv && . .venv/bin/activate && pip install -e .   # plus ortools for the CP-SAT ladder
pytest -q tests/benchmarks tests/criterion tests/estimation tests/knowledge
```

Then follow `docs/PROTOCOL_RUN.md` stage by stage; each `experiments/*.py` states in its docstring
what it writes. The supplier stage needs either ollama with qwen2.5:7b and llama3.1:8b, or the cached
answers already in `results/elicit/`. For the dated folders, each `README.md`/`RESULTS.md` gives the command; mind the absolute paths
above.

## 4. What is not here

- The paper drafts (yours of 16 September, and my outline of 13 September); ask me for the outline.
- Third-party PDFs (Guo and Perković 2010.08611; Taeb, Guo and Henckel 2511.10625): fetch from arXiv.
