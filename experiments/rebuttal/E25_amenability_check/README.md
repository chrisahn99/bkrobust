# E25: the nearest failing states are the non-amenable states (amenability check)

Tests, query by query, the consequence of Proposition C.7.1(iv) with Perković (UAI 2020, Prop. 4.2): every nearest failing state of the committed set is non-amenable, and the radius equals the amenability radius. Populations: the 3,014 released rows of the E10 frame, the two E11 frames under the eight real-name conditions, and the 4,584 stored synthetic instances; the first author's 6,949 queries (`_data/chris_e10/amenability_radius.csv`) are joined in. Registered plan: PREREG.md.

    .venv/bin/python experiments/rebuttal/E25_amenability_check/run.py --smoke && \
      .venv/bin/python experiments/rebuttal/E25_amenability_check/analyse.py --out results/rebuttal/E25_smoke
    bash experiments/rebuttal/E25_amenability_check/launch.sh 6

Files: `run.py` (level scan, submitted radius, amenability search, timings; resumable), `analyse.py` (registered shares with Clopper-Pearson intervals, comparison with the first author, exceptions with their graphs), `b_amenability.py` (literal Definitions 3.1 and 4.3 of Perković et al. 2017, a second implementation), `vendor_states.py` (E11 and E15 state builders), `vendor_ledger.py` and `vendor_search.py` (from E10).

Outputs in `results/rebuttal/E25/`: `rows.jsonl` (one line per query), `config.json`, `run.log`, `analyse.log`, `results.json`, `DONE`.
Smoke: 18 queries (6 released, 6 E11, 6 synthetic), about one second. Full run: estimated 1 to 4 h on 6 workers.
