# E3S night session: brief

Prepared 10 November 2026. Branch `iclr_2027_rebuttal_solvers` (from
`iclr_2027_rebuttal` at `c0bb2db`, E4L results). Preregistration: section E3S of
`rebuttal/PREREGISTRATION.md`, committed before any code (`8fc9677`).

## The review point

qQYf-W1: "Algorithm 1 contains two parts: the upward search with depth b and the
SAT certificate. The choice of SAT over other optimization approaches, when it
becomes preferable to explicit search, and how to select the appropriate
parameter is not thoroughly motivated."

Draft answer (author's): three parts. (1) Why a constraint solver: placeholder
"[We also tested ILP and MaxSAT formulations of the same query (new App. D.3): ...]"
-- this session fills it. (2) When explicit search is preferable: regime table,
already supported by E3 and E4L (App. D.2). (3) How to select b: b affects time
only (Remark 6); b = 3 default (App. D, E3).

## What is already known (do not redo)

- E3 (`rebuttal/e3_budget_sat.py`, `results/rebuttal/e3/budget_sat.csv`): 132
  dense queries, 12 arms (b in 0..6 x ladder/minimise), all agree; minimise
  cuts worst case 28 s -> 6.5 s at b = 3.
- E3b (`results/rebuttal/e3/large_sat_cost.csv`): the E1 encoding is O(n^4)
  constraints (145k at n = 48, 3.7M at n = 150); building it exceeds 600 s at
  n >= 300. So solver comparisons are meaningful only up to about n = 64.
- E4L (`results/rebuttal/e4l/rows.jsonl`): 1,911 queries reach the SAT stage
  (`hyb_method == "e1_ladder"`), almost all synthetic with tiny classes; listing
  beats the hybrid on 1,910 of them.
- The model's constraint types (checked on a 14-node instance): bool_or 24,664,
  bool_and (enforced) 798, bool_or (enforced) 508, linear (enforced) 276,
  linear 184, exactly_one 62, lin_max 56; linear coefficients only {1} and
  {-1, 1}; 14 integer position variables. Export with
  `model.ExportToFile(path)` then parse with `text_format` into
  `cp_model_pb2.CpModelProto` (the Python `Proto()` object of OR-Tools 9.15 has
  no `SerializeToString`).

## Solvers available

- OR-Tools 9.15 in `.venv`: CP-SAT; MILP via `pywraplp` with SCIP, HiGHS, CBC
  (checked: all three create). No Gurobi/CPLEX/Xpress.
- PySAT (RC2 MaxSAT, CaDiCaL 1.9.5): installed with the user's approval (decision log).

## Plan for the session

1. `rebuttal/e3s_translate.py`: proto -> MILP (pywraplp) and proto -> WCNF
   (PySAT, if available). Unit test on 20 small queries: every arm's radius =
   hybrid radius; every witness decodes to a state where Z fails with that many
   retractions (use `bkrobust.demo.meek.apply_orientations` on the retained
   orientations, `common.is_gac_valid_mpdag`).
2. `rebuttal/e3s_solvers.py`: driver, arms A1-A9 of the preregistration,
   resumable JSONL output to `results/rebuttal/e3s/rows.jsonl`, one process per
   query (rebuild model per arm so build time is counted per arm), SIGALRM kill
   at cap + 60 s, 8 workers, `OMP_NUM_THREADS=1`, `caffeinate`.
3. Smoke (10 queries), then a time census; if the full budget exceeds about
   8 hours, cut P-SATR to the stratified 400 sample and record that as an
   amendment in PREREGISTRATION.md before the full run.
4. `rebuttal/e3s_analyse.py`: registered statistics -> `analysis.json`;
   `rebuttal/e3s_figure.py`: time per arm against radius (log scale).
5. Commit and push results on `iclr_2027_rebuttal_solvers` (never on
   `iclr_2027`; the anonymous repo serves it).
6. Overleaf (`6ab5a4476c152c4289d9bb25`): `git pull` first.
   - New files only in `rebuttal/`: `session_report_e3s.tex/.pdf` (same style as
     `session_report_e4l.tex`): summary, design, results tables, prediction
     check, against-us outcomes, draft text for the bracket in the answer.
   - App. D.3 in `sections/appendix/d_computation.tex` only if the decision log
     says so: new `\subsection` after D.2, red, unnumbered tables (numbered
     floats would shift Tables 7-13 cited in the reviews and answers).
   - Verify with `git diff --name-status origin/main HEAD` before pushing.
7. Send the report PDF to the user.

## Constraints

- Do not modify `src/`; drivers in `rebuttal/`, outputs in `results/rebuttal/e3s/`.
- Report every outcome, including against-us ones; no post hoc arm selection.
- Laptop timings (10-core M5, 8 workers): say so; ratios are the robust part.
- Never touch José's files on Overleaf except App. D.3 if approved.

## Decision log (filled before the run)

- PySAT install: approved by the user (10 Nov 2026); `python-sat 1.9.dev21` installed in `.venv` and checked (RC2 and CaDiCaL 1.9.5 run). Arms A7 and A8 are in. Record the version in the report.
- App. D.3 in the paper: approved by the user (10 Nov 2026). Write it into `sections/appendix/d_computation.tex` after D.2 (red, unnumbered tables) and push to Overleaf `main`, together with the session report as new files in `rebuttal/`.

## Prompt to start the night session

> Run the E3S experiment prepared in `rebuttal/E3S_SESSION_BRIEF.md` on branch
> `iclr_2027_rebuttal_solvers` (preregistration E3S in
> `rebuttal/PREREGISTRATION.md`). Follow the plan in the brief, run it
> overnight, then commit and push the results, write the session report into
> the Overleaf `rebuttal/` folder as new files, and draft App. D.3 as recorded
> in the brief's decision log. Ask nothing; record any deviation as an
> amendment.
