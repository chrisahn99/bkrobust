# Table 0 — the applicability funnel

Per substrate tier: what the frame contains and what the primary supplier can certify. The intraclass correlation, design effect and effective sample size are for the certifiable indicator with the network as the cluster. Networks above 500 nodes and the one whose chain components all exceed the questionnaire cap are in the frame but not swept.

| tier | networks parsed / with frame / swept | nodes | undirected edges | components ≥ 2 | candidate treatments | candidate pairs | sampled rows | committed pairs (old gate) | certifiable (primary) | ICC / design effect / n_eff |
|---|---|---|---|---|---|---|---|---|---|---|
| T0 | 3 / 2 / 2 | 26 | 19 | 2 | 19 | 252 | 32 | 23 | 13 of 32 | 0.0 / 1.0 / 32.0 |
| T1 | 24 / 19 / 13 | 6655 | 615 | 232 | 792 | 76565 | 380 | 343 | 91 of 260 | 0.2158 / 5.1 / 51.0 |
| T2 | 8 / 8 / 8 | 120 | 41 | 14 | 45 | 482 | 156 | 117 | 99 of 156 | 0.1487 / 3.75 / 41.6 |
| T3 | 4 / 4 / 4 | 261 | 92 | 28 | 113 | 2833 | 80 | 60 | 26 of 80 | 0.1636 / 4.11 / 19.5 |

## The radius printed twice, and the Lever 0 admission count

| arm | certifiable | analyst graph keeps an undirected edge | max r_claim all / without T0 | r_claim ≥ 3 all / without T0 | max r_hop all / without T0 |
|---|---|---|---|---|---|
| `D_DEGEN` | 68 | 68 | None / None | 0 / 0 | None / None |
| `D_LLM` | 229 | 112 | 3 / 3 | 8 / 6 | 11 / 4 |
| `D_LLM_32B` | 204 | 98 | 3 / 3 | 3 / 3 | 11 / 3 |
| `D_LLM_72B` | 285 | 79 | 2 / 2 | 0 / 0 | 11 / 5 |
| `D_LLM_72B_INSTR` | 279 | 100 | 3 / 3 | 6 / 6 | 11 / 5 |
| `D_LLM_GEMMA_27B` | 263 | 71 | 3 / 3 | 1 / 1 | 4 / 4 |
| `D_LLM_INSTR` | 206 | 128 | 3 / 3 | 1 / 1 | 11 / 7 |
| `D_LLM_SRC` | 243 | 100 | 3 / 2 | 2 / 0 | 11 / 3 |
| `D_LLM_SRC_70B` | 267 | 103 | 2 / 2 | 0 / 0 | 11 / 6 |
| `D_RAND` | 649 | 388 | 3 / 3 | 12 / 11 | 6 / 6 |
| `D_SCRAMBLED` | 197 | 176 | 3 / 3 | 3 / 3 | 9 / 6 |
| `D_SCRAMBLED_72B` | 192 | 116 | 3 / 3 | 10 / 9 | 7 / 3 |
| `A_ORACLE` | 294 | 0 | 1 / 1 | 0 / 0 | 12 / 4 |
| `A_TRUE` | 209 | 99 | 3 / 2 | 1 / 0 | 12 / 2 |
