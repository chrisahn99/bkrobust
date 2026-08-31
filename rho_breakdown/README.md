# rho_breakdown

The measurements behind the breakdown-radius direction, and the code that produced them. Shared with
Chris on 2026-08-31, in answer to the request of 21 August.

Read `REPORT.md` first, or `deck.pdf` (16 slides) if you would rather have it as a talk. `deck.tex`
is the source.

`REPORT.md` and the deck carry the same content. `NUMBERS.md` maps every figure in the report to the file it came from, and
records the four corrections we applied to ourselves and the three things that have never been run.

```
x3_skeleton/ 2026-08-31  skeleton perturbation + the worked example (see its own README)
pilot/      2026-07-17   600 + 1 650 + 60 SCMs, iSCM, closed-form linear SEMs
e1prime/    2026-08-19   the `se` criterion, the censoring floor, calibration and coverage (ARM 2)
locality/x1 2026-08-19   hop-distance stratification, operator arms, pool composition
locality/x2 2026-08-19   exhaustive enumeration, no SCM and no truth, complete at p <= 5
```

Each of `e1prime/` and `locality/x*/` carries a `PREREG.md` written before the first number and an
`AUDIT.md` written against the design. Read those before the results if you want to know what was
declared in advance.

Everything is CPU-only. Requirements are numpy, scipy, pandas and networkx; nothing here needs a GPU
or a private dataset.

```bash
python pilot/code/run_linear.py          # then pilot/code/analyse.py
bash   e1prime/code/run_all.sh           # ARM 1 and ARM 2
python locality/x1/code/run_x1.py        # then analyse_x1.py, stats_x1.py
python locality/x2/code/run_lemma.py     # exhaustive, no sampling
```
