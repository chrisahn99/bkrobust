# Matched control for reviewer point 7 — does LLM knowledge explain radius-1 concentration?

Written by the orchestrator after the worker that built this arm was killed
mid-session by a rate limit. Everything below is recomputed from the committed
files in this directory, not from the worker's report.

## The question

`results/final_table/` observes that `r_val` concentrates at 1 on real networks
at LLM-elicited background knowledge — **46 of 53** finite informative instances,
86.8%. The paper attributes that concentration to the model's knowledge. The
reviewer objects that the attribution is unjustified without a control:
concentration at 1 could equally be a property of the CPDAGs, or of *how many*
claims were asserted, rather than of what the model knows.

## The arms

All arms share the questionnaire (`questionnaire_sha256` recorded in
`manifest.json`), the 32-network universe, the query seed (20260917),
`--queries-per-network 5` and `--skip pathfinder`, so the `(network, x, y)`
queries are **matched by construction**. All arms ran under one interpreter,
`.venv/bin/python` 3.14.7 with numpy 2.5.3 / scipy 1.18.1 — the interpreter that
produced the committed `D_LLM` table, established from the fact that
`certify.py:261` uses `zip(..., strict=True)` (3.10+) and the committed table
carries `r_eps` values. No shim is used anywhere.

| arm | what it controls | n | informative | share at `r_val == 1` |
|---|---|---|---|---|
| `D_LLM` (treatment, committed, reused verbatim) | — | 160 | 53 | **86.8%** (46/53) |
| `D_SCRAMBLED` | knowledge **content** (names scrambled) | 160 | 42 | 76.2% (32/42) |
| `D_SCRAMBLED_72B` | same, larger model | 160 | 31 | 83.9% (26/31) |
| `D_RANDOM_MATCHED` | knowledge **quantity** (random orientations at matched \|K\|) | 2585 | 700 | **74.9%** (524/700) |

The quantity control exists because the semantic control does not hold \|K\|
fixed: `D_LLM` asserts 171 claims across the corpus and `D_SCRAMBLED` only 151,
and `r_val` depends on \|K\|. `D_RANDOM_MATCHED` draws, per network, exactly the
`D_LLM` count of orientations uniformly from the CPDAG's undirected edges,
redrawing on Meek-inconsistency, over 20 seeds per network (2 seeds on the six
networks — `arth150`, `diabetes`, `link`, `munin`, `munin2`, `munin3` — that
time out on 100% of `D_LLM`'s queries regardless of knowledge and so contribute
no informative instance to match against; that reduction is deliberate and is
recorded here rather than made silently).

## Result — the attribution is not supported

Paired over the matched queries, clustered by network, 10,000 resamples:

| comparison | n pairs | networks | paired Δ share-at-1 | 95% CI |
|---|---|---|---|---|
| `D_LLM` − `D_RANDOM_MATCHED` | 53 | 20 | **+0.038** | **[−0.079, +0.134]** |
| `D_LLM` − `D_SCRAMBLED` | 29 | 13 | +0.103 | [0.000, +0.333] |

**Against the quantity control the difference is +0.038 with an interval that
comfortably contains zero.** Random orientations, at the same knowledge budget,
concentrate at radius 1 essentially as hard as elicited LLM knowledge does. On
this evidence the radius-1 concentration is a property of the real CPDAGs and
the size of the asserted knowledge set, **not** of the model's knowledge, and
the paper's causal attribution should be withdrawn or restated as an
observation about the corpus.

The semantic contrast is larger (+0.103) but its interval's lower bound sits
exactly at zero, so it does not establish a difference either.

## The caveat that cuts both ways

The scrambled arm is **not measurably less accurate** than the LLM arm.
`compelled_control` accuracy is 77.4% (`D_LLM`), 67.3% (`D_SCRAMBLED`) and 75.0%
(`D_SCRAMBLED_72B`), and the paired per-network difference over the 22 networks
carrying the data is **+0.092, CI [−0.058, +0.259]** — it includes zero.

That weakens the semantic control's interpretability in both directions: if
scrambling variable names did not reliably remove the model's ability to answer
correctly, then a *null* semantic contrast cannot be read as "knowledge does not
matter", because the manipulation may not have removed the knowledge. The
quantity control does not have this problem — a uniformly random draw carries no
knowledge by construction — which is why it, not the scrambled arm, carries the
conclusion above.

## Denominators, which differ by arm and are never dropped

Outcome classes per arm are in `comparison/summary.json` under
`arm_denominators`. They are not comparable as raw counts: `D_RANDOM_MATCHED`
has 2,585 instances (many seeds per query) against 160 for the elicited arms,
and its 865 degenerate / 539 error / 221 timeout rows are carried explicitly.
`UNREACHED` is `−1`, a sentinel, and is never averaged as a radius.

## Two defects found and fixed here

1. `cmd_compare` built the paired scrambled statistic with `float("nan")` for
   any instance that was informative but **not** at radius 1, then included it,
   so the whole paired mean came back `NaN`. Every non-radius-1 instance poisoned
   it. Fixed to contribute `0.0`.
2. `control_accuracy` and `per_network_accuracy_diffs` assumed every network
   carries a `compelled_control` block. Networks whose elicitation never reached
   a compelled question carry `None`, which raised `TypeError`. Those networks
   are now skipped and counted, never imputed as agreement.

## Reproduce

```
PYTHONPATH=src:experiments .venv/bin/python experiments/matched_control.py verify
PYTHONPATH=src:experiments .venv/bin/python experiments/matched_control.py table-arm --condition D_SCRAMBLED
PYTHONPATH=src:experiments .venv/bin/python experiments/matched_control.py random --networks <list> --seeds 20 --shard-tag randA
PYTHONPATH=src:experiments .venv/bin/python experiments/matched_control.py merge-random
PYTHONPATH=src:experiments .venv/bin/python experiments/matched_control.py compare
```

`r_eps` on real networks is SEM-conditional (Σ comes from a seeded
linear-Gaussian SEM attached to the ground-truth DAG); `r_val`, which carries
every number above, is purely graph-theoretic and carries no such assumption.
