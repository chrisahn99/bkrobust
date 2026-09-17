# Paper integration brief — the real-network radius table at LLM-elicited knowledge

**Who this is for.** A session sitting down to update `ICLR_2027___bkrobust_v0` with the
final-table work on branch `experiments/final_table`. You do not need to read any code:
§3 below maps every claim to the file it comes from, and §4 says where each edit goes.

**Three house rules this brief inherits.**

1. **Never invent a number.** Every figure below traces to `results/final_table/`. If a
   number you want is not in §3, it is a `TODO`, not a guess.
2. **The `r` vs `r−1` convention** (`src/bkrobust/core/conventions.py`): `r` is the distance
   to the nearest failure; the number of *safe* moves is `r−1`. This matters more here than
   in any previous session, because the headline value is `r = 1`, which means **zero** safe
   retractions — see §3.2. Do not paraphrase `r = 1` as "tolerates one error".
3. **The radius is called `r_val`, everywhere.** `r_val` is the name in the code, in every
   committed result, in `epsilon-radius.tex` and in the paper, and it has no alias — earlier
   drafts of this brief used one and were wrong. When the *unit* rather than the radius is
   meant, write "hop units" and "claim units" in words: the unit axis and the property axis
   are different distinctions, and only the latter names radii (`r_val`, `r_opt`, `r_eps`).

---

## 1. The two things that change the paper

### 1.1 The ε-radius section's own stated limitation is closed

`PAPER_INTEGRATION_R_EPSILON.md` §5 lists, among what must not be claimed, that **the
real-graph corpus was not touched** — the ε-radius construction was evaluated on synthetic
instances only. It has now been run on 32 real benchmark networks, at background knowledge
elicited from language models rather than assumed.

No previously reported number changes — what changes is that the sentence "evaluated on
synthetic instances" can be replaced by a real-network panel, and the limitation bullet can
be struck.

### 1.2 On real networks, at LLM knowledge, the radius is 1 — and tolerating bias does not raise it

Two findings, and the second is the uncomfortable one.

**`r_val = 1` almost everywhere.** Of the informative instances, 46/53 have `r_val = 1`, and
16 of the 20 networks with a defined median sit at 1. By the `r−1` convention that is *zero*
safe retractions: the committed adjustment set stops being valid at the very first wrong
orientation the model asserted.

**Widening the question from validity to bias buys almost nothing here.** `r_eps` equals
`r_val` at every ε in `{0.01, 0.05, 0.10, 0.25, 0.50, 1.00}` for all but 2 of 53 informative
instances, and at network level only `paths` gains (1 → 11, at ε = 100%). On synthetic
instances the ε-refinement bought a great deal; on this corpus, at this knowledge source, it
does not. That **qualifies** the ε-radius story rather than extending it, and §5 says how far
that qualification may be pushed.

---

## 2. Files to read, in priority order

### Tier 1 — read these

| File | What it is |
|---|---|
| `docs/PAPER_INTEGRATION_FINAL_TABLE.md` | This brief. |
| `results/final_table/summary.json` | Per-network aggregates, the ε grid, the CLI used, and the `assumptions` block that must travel with every `r_eps` number. |

### Tier 2 — evidence, when you need to cite or check a number

| File | What it holds |
|---|---|
| `results/final_table/instances.jsonl` | One record per (network, x, y): `r_val`, `r_eps` per ε, `theta_z`, `n_knowledge`, `status`, `seconds`. 160 rows. Every number in §3 is recomputable from this file alone. |
| `results/elicit/knowledge.json` | The elicited knowledge sets K themselves, per condition per network, with `k`, `k_sha256` and the repair counts. |
| `experiments/final_table.py` | The generating script. Module docstring carries the assumption statement. |

### Tier 3 — how to regenerate

| Purpose | Command |
|---|---|
| The sweep (≈74 min) | `PYTHONPATH=src .venv/bin/python experiments/final_table.py --queries-per-network 5 --time-limit 120 --out results/final_table` |
| Re-render the table only | `PYTHONPATH=src .venv/bin/python experiments/final_table.py --report-only --out results/final_table` |
| A different knowledge source | add `--condition D_LLM_72B` (ten conditions exist; see `results/elicit/knowledge.json`) |

---

## 3. The findings, with their evidence

### 3.1 What was actually run

32 real networks from `results/axisa3/networks/example_models/`, each at its LLM-elicited
knowledge set K for condition `D_LLM`, on the first 5 (x, y) queries per network from the
frozen frame `results/frame/frame.jsonl`. 160 instances. One `epsilon.certify` call per
instance returns both radii from a single traversal.
→ `results/final_table/summary.json` (`args` block)

### 3.2 The radius is 1, meaning zero safe moves

Among the 53 informative instances with a finite radius, the distribution of `r_val` is:

| `r_val` | 1 | 2 | 3 | 4 | 11 |
|---|---|---|---|---|---|
| instances | 46 | 2 | 3 | 1 | 1 |

**46/53 = 87% sit at `r_val = 1`.** Taking networks as the unit of analysis, as
`PAPER_NARRATIVE.md` requires, 16 of the 20 networks with a defined median are at 1, and the
median over networks is 1.
→ `results/final_table/instances.jsonl`, `results/final_table/summary.json`

Stated in the paper's own units: **the analyst has zero safe moves.** One wrong orientation,
out of a set the model asserted with no access to the truth, invalidates the adjustment set.

### 3.3 The denominators, which are not flattering and must be reported

160 instances do not yield 160 numbers:

| Outcome | n | Meaning |
|---|---|---|
| informative | 66 | A radius was computed. |
| degenerate | 33 | `r_val = 0`, `method = "degenerate"`: the optimal adjustment set at G₀ is **empty and already invalid**, so the query is ill-posed, not fragile. Excluded from every median. |
| timeout | 34 | Exceeded the 120 s per-instance budget. |
| no identified Z | 27 | `certify` raised: no adjustment set is identified at G₀, so there is nothing to certify. |

Only **21 of 32 networks** yield at least one informative query. Eleven are blank:
`arth150`, `diabetes`, `hailfinder`, `insurance`, `link`, `munin`, `munin1`, `munin2`,
`munin3`, `munin4`, `sachs`.
→ `results/final_table/summary.json` (`n_ok`, `n_degenerate`, `n_timeout`, `n_error` per network)

The degenerate class is the one most likely to be misread. `r_val = 0` is `hybrid`'s verdict
that the question is degenerate (`hybrid.py:224`), verified directly: on `child`, three zero
instances all return `method='degenerate'` with `|Z| = 0`. Averaging those zeros into a
robustness statistic would have understated the radius by a third of the sample; an earlier
draft of the aggregation did exactly that, and it is why the `ok/dg/to/er` column exists.

### 3.4 The ε-refinement does not pay on this corpus

`r_eps` exceeds `r_val` at some ε in only **2 of 53** informative instances (`Kampen_2014`,
`Schipf_2010`). At network level the median `r_eps` equals the median `r_val` at every ε in
the grid; the single visible gain is `paths`, 1 → 11 at ε = 100%.
→ `results/final_table/instances.jsonl`

The honest reading is that on real graphs, at LLM-elicited knowledge, the bias crosses any
tolerance essentially at the validity boundary. The "ramp, not a step" behaviour that the
synthetic sweep demonstrates is not contradicted — it is simply not *reached* here, because
the radius is already 1 and there is no room between 0 and 1 for a ramp to be visible.

### 3.5 Why eleven networks are blank, and what it is not

The timeouts are **not** shell enumeration. `link` has |K| = 2, i.e. four retraction shells,
and still exhausts its budget. Profiling shows the cost is upstream of any radius:
`optimal_adjustment_set_mpdag(g0, x, y)`, which `certify` calls whenever `z` is not supplied,
alone exceeds 90 s on `link`'s 724 nodes. `arth150` (107 nodes, |K| = 14) is the genuinely
shell-bound case. Two different causes, and only the second would be helped by a bounds-based
strategy.

---

## 4. Concrete edits to the existing paper

| Where | Edit |
|---|---|
| `paper/sections/epsilon-radius.tex`, §`sec:epsilon-eval` | Add a new `\paragraph{On real graphs}` after the existing evaluation paragraphs (end of file, after line 269), carrying §3.2 and §3.4. Keep it prose: this section contains no `table`/`tabular` environments today — see §6.1 before introducing one. |
| `PAPER_INTEGRATION_R_EPSILON.md` §5 | Strike the "the real-graph corpus was not touched" bullet; replace with the narrower scope statement in §5 below. |
| `PAPER_NARRATIVE.md` §5, Result B | Add a sentence noting the real-network panel now exists but does **not** close `[RE-11]` — see §5. |
| `PAPER_NARRATIVE.md` §10 (provenance map) | Add `results/final_table/instances.jsonl` and `summary.json`. |
| Abstract / intro, only if §6.2 is resolved in favour of foregrounding | The "zero safe moves at LLM-elicited knowledge" claim is the strongest sentence this work supports. It is also a negative result about the knowledge source, not about the instrument. |
| Reproducibility statement | Add the two commands in Tier 3. Note that the sweep needs no network access and no model API. |

---

## 5. What must NOT be claimed

- **This table is not pre-registered.** `docs/EVALUATION_PLAN.md` names `r_claim`, `r_val`,
  `phi_1` and an identifiability radius; `r_eps` appears nowhere in its radius list or table
  budget. Present this as a new result, not as the delivery of a planned one.
- **It does not close `[RE-11]`.** That item is specifically worst-case radius against
  *average-case survival* on real instances. No survival AUC is computed here. Citing this
  table as `[RE-11]`'s completion would be wrong.
- **`r_eps` on real networks is SEM-conditional; `r_val` is not.** The corpus ships graph
  structure only, so Σ comes from a seeded linear-Gaussian SEM attached to the ground-truth
  DAG (`demo.evaluate.random_sem`). `r_val` is purely graph-theoretic and carries no such
  assumption. Any sentence reporting an `r_eps` number on real data must carry the assumption
  with it; the `assumptions` block in `summary.json` is the canonical wording.
- **Eleven blank networks are unknown, not robust.** They timed out or were degenerate. They
  must never be reported as `UNREACHED`, and must not be silently dropped from the
  denominator: 21/32 is the honest network count.
- **`UNREACHED` is never averaged.** It is `−1` in the JSON and a sentinel, not a radius.
  13 informative instances and the whole of `win95pts` are `UNREACHED`.
- **Degenerate queries are not fragile queries.** The 33 `r_val = 0` rows mean the adjustment
  set was already invalid at G₀ with an empty Z. Do not fold them into a "radius is low"
  claim; they say the query was ill-posed under that K.
- **The ε-null result is about this corpus and this knowledge source.** It does not show the
  ε-refinement is useless in general — the synthetic sweep shows otherwise. Say "does not
  separate from `r_val` here", not "does not separate".
- **One condition only.** Everything above is `D_LLM`. Nine other conditions (larger models,
  scrambled names, a second presentation order) sit unused in `results/elicit/knowledge.json`.

---

## 6. Open decisions, deliberately not taken

**6.1 Table or prose.** `paper/sections/epsilon-radius.tex` reports every evaluation result as
prose under `\paragraph` headers and contains no tabular environment at all. A 32-row × 6-ε
table would be a structural departure from that section's convention. The full table is
rendered by `--report-only` if you want it; the alternative is three sentences carrying §3.2
and §3.4 and a pointer to `summary.json`. **Decide before drafting the paragraph.**

**6.2 How far to foreground a negative result.** §3.2 and §3.4 together say the instrument
finds almost no robustness to measure at LLM-elicited knowledge. That is a real finding about
the knowledge source and arguably the most quotable sentence available — but it sits awkwardly
beside a paper whose contribution is the instrument. This is an authorial call.

**6.3 Whether to fill the eleven blanks, and at what cost.** Passing `z=` explicitly to
`certify` sidesteps the `optimal_adjustment_set_mpdag` bottleneck (§3.5) and would likely
recover the large-|V| networks. But both radii are *defined relative to the analyst's
committed adjustment set*, so substituting a cheaper Z changes the question being asked rather
than computing the same answer faster. It is a scientific choice, not an optimisation, and it
is not taken here.

**6.4 Nothing outstanding on naming.** The repo previously carried a conflict here: §11 of
`PAPER_NARRATIVE.md` mandated a separate name for the geometric radius and forbade `r_val`,
while the code, `epsilon-radius.tex` and `PAPER_INTEGRATION_R_EPSILON.md` all used `r_val`.
That conflict is resolved in favour of `r_val`, which is now the single name used repo-wide.
Nothing is left to decide.
