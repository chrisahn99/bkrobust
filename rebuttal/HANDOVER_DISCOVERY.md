# Handover: discovery error (E5 and beyond) for the ICLR 2027 rebuttal

Written 10 November 2026 (evening), at the end of a long session. It is for a
fresh Claude Code session started in this repository. Read it fully before
doing anything. Then read the files it points to under "Read first".

---

## 0. One-paragraph context

The paper is ICLR 2027 submission 46323, "breakdown radius" of an adjustment
set. Given a CPDAG Ĉ, the analyst's orientation knowledge (MPDAG G0) and a
committed adjustment set Z (the optimal set O, criterion GAC), the radius
`r_val` is the fewest retractions of the analyst's closure orientations
K_G0 after which Z is no longer valid. The radius is a **sensitivity analysis
of the orientation knowledge, conditional on the equivalence class being
correct**. The current reviews are qQYf 3, C2Bd 3 and h8pC 2 (weak reject).
The discussion closes **18 November 2026, 23:59 AoE** (19 Nov 12:59 Paris).
The first round of answers is already posted on OpenReview.

---

## 1. Where things are

| what | path |
|---|---|
| Code repo (this one) | `/Users/ahn/Documents/Research/iclr27/bkrobust` (GitHub remote `origin`) |
| Python env | `.venv/` (run everything as `PYTHONPATH=src:rebuttal .venv/bin/python ...`) |
| Overleaf project (git clone) | `/Users/ahn/Documents/Research/iclr27/6ab5a4476c152c4289d9bb25` (`git pull` before any edit; push to `main`) |
| Rebuttal drivers and outputs | `rebuttal/*.py`, `results/rebuttal/<exp>/` |
| Preregistrations (ours) | `rebuttal/PREREGISTRATION.md` (one section per experiment, amendments appended, never edited away) |
| Overleaf rebuttal folder | `rebuttal/` in the Overleaf clone: `reviews.tex/.pdf`, `responses.tex/.pdf`, `results.tex/.pdf`, `session_report*.tex/.pdf`, `README.md`, `changes_note.md` |
| Paper appendix on discovery error | Overleaf `sections/appendix/g_experiments.tex`, subsection "Learned equivalence classes", `\label{app:discovery}`, Remark `rem:discovery` (printed as App. H.4, Remark 8) |

### Branches (GitHub, code repo)

| branch | role | push? |
|---|---|---|
| `iclr_2027` | **the submitted code; the anonymous repository serves it** | **NEVER push, never commit** |
| `iclr_2027_rebuttal` | first rebuttal session (E1, E2, E3, E5/E6/E7, E8, E10–E18, E21, E4L, Figure 3 toy generator); cut from `iclr_2027` at `5421cab` | yes |
| `iclr_2027_rebuttal_solvers` | E3S (SAT stage vs MILP/MaxSAT/SAT/explicit), merged from `iclr_2027_rebuttal`; head `e0ae89e` | yes |
| **`iclr_2027_rebuttal_discovery`** | **this work.** Cut from `iclr_2027_rebuttal_solvers`, so it contains everything above. This handover is its first commit. | **yes, work here** |
| `ze/rebuttal` | José's (co-author) rebuttal experiments E1–E29, under `experiments/rebuttal/<E..>/` with PREREG.md, run.py, analyse.py and FIDELITY.yaml. His learned-graph data is in `experiments/rebuttal/_data/pest/` (graphs.jsonl, graphs_retry.jsonl, frame.jsonl, knowledge.json, rows.csv). His **result files are not in git**; the numbers are only in Overleaf `rebuttal/results.tex`. | **read only, never push** |
| `main`, `revision/*`, `experiments/*`, `iclr_2027_v2/v3` | older history | do not touch |

Read a file from José's branch without checking it out:
`git show origin/ze/rebuttal:experiments/rebuttal/E5_learned_cpdags_by_sample_size/PREREG.md`

---

## 2. Standing rules (from the user, all still in force)

- **GitHub:** commit and push only on rebuttal branches (`iclr_2027_rebuttal*`). Never on `iclr_2027`. Do not modify `src/`; drivers go in `rebuttal/`, outputs in `results/rebuttal/<exp>/`.
- **Overleaf:**
  - New material goes in **new files** in `rebuttal/`.
  - Paper changes go in red (`\rev{...}` or `{\color{red} ...}`).
  - Use **no numbered floats in the appendix.** Use unnumbered `center` tables or figures, because a new numbered float would shift Tables 7–13, which the reviews and answers cite.
  - **Do not edit José's files**: `results.tex`, `README.md`, `changes_note.md`, and his sections, except where the user approves a specific edit.
  - Before pushing, verify with `git diff --name-status origin/main HEAD`. Build with `latexmk` (the paper builds into `build/`) and check for errors and new overfull boxes.
- **`responses.tex`:** this is the text **as posted**. Red text means additions not yet posted. Its header says which additions the red covers; update that line when you add red.
- **Science norms:**
  - Write the preregistration in `rebuttal/PREREGISTRATION.md` **before any run**: question, populations, registered statistic, prediction, an "against us" clause, consistency checks.
  - Record every deviation as a dated amendment.
  - Report every outcome, including against-us ones.
  - No post hoc selection of arms or subsets presented as registered. Exploratory analyses are labelled exploratory.
- **Machine:** a 10-core Apple M5 laptop, 32 GB. Use 8 workers, one thread each (`OMP_NUM_THREADS=1`), and run long jobs under `caffeinate`. **Label all timings as laptop timings.**
- **Installing packages:** ask the user before installing anything. PySAT was installed with approval. `causal-learn`, `pytetrad` and `gcastle` are **not** installed.
- **Privacy:** the user's email is private. Never put it in a commit, URL or request.
- **Commit trailer:** `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.

---

## 3. Assessment of the three reviews and the state of every point

Reviews, word for word with point labels: Overleaf `rebuttal/reviews.tex`.
Answers as posted: `rebuttal/responses.tex`.

### qQYf (3, weak accept, confidence 3)
- **S2 (practical motivation):** answered.
- **W1 (why SAT; when explicit search; how to pick b):** answered.
  - New since posting (red in responses.tex, not yet posted): the E4L paragraph (listing the class) and the E3S paragraph (MILP/MaxSAT/SAT comparison, new App. D.3).
  - CP-SAT is the fastest solver. Explicit search wins while the retraction space is small and fails at radius ≥ 11.
- **M2, Q1 (uniqueness of the optimal set; invalid set with a larger radius):** answered. Every set valid in G0 shares one radius (Prop. C.7.1). App. C.7 was restructured with a summary table, a diagram and a 3-node example.
- **M3, M4:** wording, answered.

### C2Bd (3, weak accept, confidence 2) — the "second reviewer"
- **W1, Q2 (visual intuition):**
  - New Figure 2, the lattice walkthrough.
  - The Figure 3 toy was replaced on 10 Nov by one whose radius really is 2 (it had been printed as 3; the true value was 1).
- **W2 (THE TARGET OF THIS SESSION):**
  - Reviewer's words: "The method assumes a correct skeleton and v-structures; finite-sample discovery errors or unmeasured confounders outside CPDAG are not captured. It would strengthen the paper if it included some experimentation on how sensitive the method is to these."
  - Posted answer: concedes that the radius is conditional on the class, and points to App. H.4. It gives:
    - Remark 8: an extra adjacency cannot break Z within the budget, while missing adjacencies and false or missed v-structures can;
    - an extra edge never broke a committed set (0 of 1,229,413 cases);
    - a missed edge two or more hops away rarely did (0.00062 of 35,488).
  - **No positive finite-sample result is in the paper or the answer yet.**
- **W3 (separation does as well):** conceded, with nuance.
- **W4 (cohort tables):** answered.
- **Q1 (standardisation):** answered (E18).
- **Q3 (couple with skeleton-stability or bootstrap checks?):** posted answer:
  - no check of that kind comes with a certificate;
  - we tried bootstrap stability, parametric bootstrap and PC–GES agreement, and none is in the paper;
  - every invalid set traces to a missing adjacency or reversed arrow among the ancestors of X, Y, Z.

### h8pC (2, weak reject, confidence 3) — the score to move
- **W1 (complexity):**
  - App. C.8: NP membership, an FPT result, NP-hardness open.
  - App. D.2: timings up to 400 nodes.
  - New red text: E4L listing comparison.
- **W2 ("The paper assumes that the Markov equivalent class is correct ... the metric neglects the bias caused by the recovery of the CPDAG"):**
  - Same family as C2Bd-W2.
  - Posted answer: agree; the radius is conditional; App. H.4.
  - **A positive finite-sample result would help here too.**
- **W3 (novelty vs graph-space distances, Taeb et al. 2026):** answered.
  - The poset is graded.
  - The radius equals the distance to the nearest non-amenable state (Cor. C.7.1).
  - AID baseline: +0.19.
- **W4(1), W4(2) (small networks; pair selection):** answered (E12, E11; τ_b 0.63 on CPDAG-only queries).

---

## 4. Everything already measured on discovery error (do not redo blindly)

### Ours (branch `iclr_2027_rebuttal`; files `rebuttal/discovery.py`, `rebuttal/e5_discovery.py`, `rebuttal/e5_analyse.py`; outputs `results/rebuttal/e5/`)

**Setup**
- In-house PC-stable (`rebuttal/discovery.py`, numpy/scipy only): Fisher-z tests, α = 0.01, max depth 3; v-structures strongest first; one Meek closure.
- Validated against an oracle: it returns the true CPDAG with oracle tests (`rebuttal/pc_oracle_check.py`).
- Two earlier buggy runs are kept in `results/rebuttal/e5/superseded_*` and are **not used**:
  - Meek closure after every v-structure;
  - Python `hash` seeding.

**Networks and data**
- Four bnlearn Gaussian networks with fitted linear-Gaussian CPDs: ecoli70, magic-irri, magic-niab, arth150 (44 to 107 nodes; `D.GAUSSIAN`).
- Data are sampled ancestrally.
- n ∈ {200, 1000, 5000}, 5 seeds each, plus the oracle CPDAG.
- 40 frame pairs per network: y is a descendant of x in the truth, and x touches an undirected edge of the oracle CPDAG.

**Knowledge arms on the learned CPDAG**
- TRUE: true orientation of each undirected edge.
- LLM: the released D_LLM claims.
- NONE.

**Results (`results/rebuttal/e5/analysis.json`)**
- TRUE arm, share of committed sets invalid in the true DAG:
  - n=200: 0.336 (116 sets); n=1000: 0.406 (283); n=5000: 0.373 (260); oracle: 0.
  - The share does **not fall** with n.
- The certifiable share rises with n: 0.145 → 0.354 → 0.325, against 1.0 on the oracle.
- Invalidity is concentrated near skeleton errors: 0.41 with an error near the query against 0.235 without (TRUE); AUC 0.71.
- By radius bin on learned graphs (TRUE): r=1 0.216, r=2 0.15, r≥3 0.19, **r=∞ 0.587**.
  - Infinite radii on learned classes are often artefacts of missing edges.
  - AUC of 1/r for invalidity: 0.34, i.e. the wrong direction.
- E6, bootstrap stability as a flag:
  - AUC 0.68;
  - the 7 of 284 sets stable in all 30 replicates were all valid;
  - at stability ≥ 0.9, 14% are invalid, against 49% among unstable sets.
- E7, one hidden variable: no average effect (0.34 at hop 1, 0.32 further away, 0.30 baseline). A hidden common cause of X and Y gives 0.45 (5 of 11).
- E8 (`rebuttal/e8_edges.py`), population level:
  - an extra edge never breaks a set (0 of 23,895);
  - a missed edge breaks one only within one hop of {X, Y} ∪ Z;
  - Proposition: deleting edges from a DAG keeps a valid set valid (proved in `session_report.tex` and as Remark 8(i) in the paper).

### José's (branch `ze/rebuttal`; numbers from Overleaf `rebuttal/results.tex`)

**E5 setup**
- PC **and GES** at n ∈ {200, 1000, 5000}, every stored seed: 214 learned graphs, plus the true CPDAG as anchor.
- Graphs come from the stored P-EST panel (`_data/pest/graphs*.jsonl`). Each graph record carries `network`, `n`, `seed` and `method` (`pc` or `ges`); the learner settings are not in this file, so read the P-EST code on `main` (`results/pest`, commit `ba7ddc12`) or ask the user.
- Arms: D_LLM, D_RAND, D_DEGEN, A_TRUE, A_ORACLE, A_ORACLE_R.
- Registered statistic: Δ = invalid share(n=5000) − share(n=200), A_ORACLE_R, PC.

**E5 results**
- PC: Δ = **+0.121 [+0.014, +0.221]**, the wrong direction (against us). Without arth150: +0.105 [−0.005, 0.217].
- **GES: Δ = −0.295 [−0.396, −0.194]: invalidity FALLS with n under GES.** This is the most promising existing positive finding; it is unregistered for GES (reported "also").

**Other rows**
- E6 (bootstrap flag): AUC pooled 0.245, stratified 0.442; against us.
- E9 (parametric bootstrap): AUC 0.59; no claim.
- **E22, census of what broke each invalid set (for us):**
  - 0 of 534 invalid sets have extra edges only near the query;
  - 0 of 534 have no error within two hops;
  - every invalid set traces to a missing or reversed true edge near the query.
- E23, PC–GES agreement as a flag: AUC 0.616 [0.522, 0.725]; flags 51% of valid sets; against us.
- **E24, far discovery errors kept and truth restored near the query (for us):**
  - invalid share 0.021 [0.012, 0.038] at d=2, pooled;
  - PC alone 0.038 [0.021, 0.068].
  - **Discovery error far from the query is nearly harmless.**
- E28: joint radius over retractions plus edge additions; read its README.

### Our reading so far
- The **global** invalid share on learned classes is dominated by local skeleton or orientation errors around X, Y, Z.
- Under PC with a fixed α that local error does not shrink at n ≤ 5000 on these networks; under GES it does.
- Nobody has yet measured:
  - sample sizes above 5000;
  - other CPDAG-producing learners;
  - whether the radius computed on the learned class **converges to the oracle radius**;
  - whether the **within-class** claims of the paper (the radius ranks survival under knowledge corruption) **survive** on learned classes.

---

## 5. The objective of this session

The user's framing (keep it):
- The radius is a sensitivity analysis **on the knowledge**, conditional on the class.
- Assessing "discovery + knowledge" as a whole pipeline is partly ill-posed.
- But reviewer C2Bd asked (and h8pC-W2 echoes it), so: **find as many solid quantitative results as possible that support the work and its relevance**, and turn the negative E5 finding into a positive one where the data honestly allow.
- In particular:
  - study in depth how finite-sample causal discovery limits the framework;
  - find **at what sample size results become coherent**;
  - try **other CPDAG-producing discovery methods** besides PC and GES.

Honesty constraint (project norm, non-negotiable): preregister, report against-us outcomes, and label exploratory analyses. "Turning negative into positive" means finding the **right questions** (decompositions, convergence, conditional guarantees), not hiding results.

---

## 6. Suggested plan (a starting point; refine, then preregister as E30+)

Use ids **E30 onwards**: José used E1–E29 and we used E1–E21, E3S and E4L.

**E30 — sample-size ladder and the coherence threshold.**
- Setup:
  - Same four Gaussian networks, plus the synthetic generators if useful.
  - n ∈ {100, 200, 500, 1k, 2k, 5k, 10k, 20k, 50k, 100k}, at least 10 seeds.
  - Knowledge arm TRUE (truthful orientations on the learned class; A_ORACLE_R in José's terms) as the registered arm; LLM and NONE as secondary arms.
- Per graph, record:
  - SHD to the true CPDAG;
  - skeleton precision and recall;
  - v-structure precision and recall;
  - **local** error, defined as an error among the ancestors of {X, Y} ∪ Z or within one hop.
- Per query, record:
  - committed or not;
  - Z valid at truth;
  - r on the learned class and r on the oracle class.
- Registered statistics: the invalid share by n with bootstrap over seeds, and the **first n at which the invalid share is within δ of the oracle** (define δ beforehand, e.g. ≤ 0.05), per method.
- Prediction: falls with n for consistent methods (GES with BIC; PC with α decreasing in n); may plateau for PC with fixed α. Faithfulness violations or near-violations in the fitted networks may cap it. Report the plateau if so.

**E31 — other CPDAG learners.**
- Candidates (all output CPDAGs or patterns):
  - PC with α schedules (0.05, 0.01, 0.001, α_n = 1/n or 0.05/√n);
  - conservative PC and majority-rule PC (more careful v-structures, the error type Remark 8 says matters);
  - GES/fGES (BIC, penalty discount);
  - **GRaSP** and **BOSS** (permutation-based; strong finite-sample accuracy in the Tetrad benchmarks);
  - MMHC/hill-climbing DAG converted to a CPDAG;
  - NOTEARS or DAGMA DAG converted to a CPDAG (beware varsortability; fine as an extra arm).
- Implementation options:
  - `causal-learn` (PC, GES, GRaSP, BOSS);
  - `pytetrad` (needs Java);
  - in-house code.
- **Ask the user before installing.**
- Validate every learner with oracle tests or the oracle CPDAG where possible.
- Compare all learners on the same samples.

**E32 — the conditional guarantee (likely the strongest positive result).**
- Among queries whose **local** structure is correct, the radius's guarantee holds exactly: Z is valid at truth whenever the truth agrees with G0 within the budget (Remark 8(i) / Cor. 1(c)). Concretely: no missing or reversed edge among the ancestors of {X, Y} ∪ Z, and no false or missed v-structure there.
- Measure the share of queries that are locally correct as a function of n and learner. This share should rise toward 1 for consistent learners even when global SHD stays high.
- Then the headline becomes: "the radius's certificate fails only through local discovery errors, and the probability of a local error falls to X% by n = N".
- E22 and E24 (José) already point this way. A per-query local-correctness check is cheap and gives the user a practical diagnostic.

**E33 — convergence of the radius itself.**
- Measure P(r_learned = r_oracle) and Kendall τ between r_learned and r_oracle, by n and learner, on the frame queries.
- Prediction: rises with n for consistent learners. A "the radius is a consistent estimand" result would be strongly positive.

**E34 — the paper's within-class claim on learned classes.**
- Run the Section 6 survival experiment (corrupt the analyst's claims, re-close, check Z's validity in the corrupted state; survival AUC; Kendall τ_b between radius and survival) with Ĉ learned at each n.
- The radius should still rank survival, because survival is defined within the class. This shows the method's purpose (sensitivity to knowledge) is intact under discovery error.
- Read `rebuttal/e11_*.py` and `rebuttal/e15_*.py` for the frame and the survival machinery, and the paper's Section 6.

**Optional:**
- A latent-variable arm with FCI (or RFCI) → PAG, measuring how often the query is still amenable in the PAG. This is outside the paper's scope: the method assumes causal sufficiency, and a PAG extension is future work. Only do it if cheap.
- Sweeping the PC α for E6-style flags.

**Deliverables at the end:**
1. Results committed and pushed on `iclr_2027_rebuttal_discovery`.
2. A session report `session_report_discovery.tex/.pdf` (+ figures) as **new files** in Overleaf `rebuttal/`. Use the style of `session_report_e4l.tex` or `session_report_e3s.tex`:
   - summary;
   - design;
   - tables;
   - a preregistered-prediction table with holds / partly / fails;
   - things to be careful about;
   - **draft answer text** for C2Bd-W2, C2Bd-Q3 and h8pC-W2.
3. If the user approves, a red addition to App. H.4 (`g_experiments.tex`, after the existing H.4 paragraphs; unnumbered table or figure only), and red paragraphs in `responses.tex` under C2Bd-W2 and h8pC-W2. Update the red-text header note.
4. Send the report PDF with SendUserFile.

---

## 7. Pitfalls already hit (save yourself the time)

- **PC bugs:**
  - Meek closure must run once, after all v-structures. Closing after each one lets R1 pre-empt true v-structures.
  - Never seed with Python's `hash()`, which is per-process. Use `hashlib`.
- **Empty optimal set.** O can be empty (∅ is a valid adjustment set). Do not assert Z non-empty. This cost a rerun in E3S.
- **Compute limits.**
  - On graphs above about 60 nodes, José skipped the SAT ladder and used search with a lower-bound fallback.
  - `bkrobust.search.exact_fast.radius_local_up_fast` with no budget explodes beyond radius about 7 on dense graphs.
  - The E1 SAT model is O(n^4) and too slow to build at n ≥ 300.
  - Use the hybrid (`common.hybrid`) or CP-SAT minimisation for deep radii.
- **Helpers live in `rebuttal/common.py` and `rebuttal/oset.py`:**
  - `MPDAG`, `apply_orientations`, `is_closed`, `k_of`, `is_gac_valid_mpdag`, `hybrid`, `OUT`, `RELEASE`;
  - `optimal_set_closed_form(g0, x, y)` returns `(z, status)`.
- **Network files:**
  - The Gaussian network files are loaded by `discovery.py` (see `D.GAUSSIAN` and the loader there).
  - `results/axisa3/networks/example_models` also holds small textbook DAGs.
- **José's pipeline** reads `_data/pest/` on his branch. His E5 rows are not in git. If you need them, ask the user (they are on betelgeuse). The Overleaf `results.tex` numbers are authoritative for what he found.
- **Overleaf builds:**
  - `latexmk` writes the paper to `build/main.pdf`.
  - Rebuttal docs build in place. Run `latexmk -c` and do not commit `.aux`, `.log` or similar files.
- The paper appendix letters as printed: D = Computation (D.1 encodings, D.2 running time, D.3 solvers); H.4 = learned equivalence classes (it lives in `g_experiments.tex`); C.7 = all valid sets share one radius; C.8 = complexity.

---

## 8. Read first (in this order)

1. Overleaf `rebuttal/reviews.tex` (C2Bd and h8pC sections) and `rebuttal/responses.tex` (C2Bd-W2, C2Bd-Q3, h8pC-W2 answers).
2. Overleaf `sections/appendix/g_experiments.tex` from `\subsection{\rev{Learned equivalence classes}}` to the end of that subsection (Remark 8 and its proof).
3. Overleaf `rebuttal/results.tex`, the "Discovery error" block (E5–E9, E22–E24).
4. Overleaf `rebuttal/session_report.tex`, `\section{Discovery error}` (our E5–E8 write-up).
5. This repo: `rebuttal/PREREGISTRATION.md` (top: plan and E5 expectation; bottom: the E3S section as a template), `rebuttal/discovery.py`, `rebuttal/e5_discovery.py`, `rebuttal/e5_analyse.py`, `results/rebuttal/e5/analysis.json`.
6. José's `git show origin/ze/rebuttal:experiments/rebuttal/E5_learned_cpdags_by_sample_size/PREREG.md`, plus the E22 and E24 READMEs and PREREGs.

---

## 9. Prompt to start the new session

The session is started from the parent folder `/Users/ahn/Documents/Research/iclr27`
(not a git repo; it holds this repo as `bkrobust/` and the Overleaf clone as
`6ab5a4476c152c4289d9bb25/`, plus older snapshots to ignore). The full start prompt,
with the folder map, branch checkout and rules, is in
`/Users/ahn/Documents/Research/iclr27/PROMPT_discovery.md`. All paths in this
handover are relative to `bkrobust/`; "Overleaf" means `../6ab5a4476c152c4289d9bb25/`.
