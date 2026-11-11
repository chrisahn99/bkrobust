# E27 deviations from PREREG.md

1. (2026-11-07, about 02:05 CET, during the first full pass.) Items on the munin family (first
   D_LLM_SRC|munin and D_LLM_SRC_70B|munin, 26 items by 02:13, among them `paths` items handled by a worker whose
   network cache already held a munin graph; the final list is errors.shard-0-of-1.jsonl) raised
   MemoryError in `build_state` or in network loading (munin family, 1,041 nodes) under the 1.5 GB per-worker data cap set by `launch.sh`. The error was raised
   outside the per-call guard, so the whole item wrote no row (errors.shard-0-of-1.jsonl). The cap is
   a machine-safety limit, not a budget of the statistic; the first pass is kept as is. A second pass,
   `resume_after_done.sh`, waits for DONE and reruns only the units missing from disk (those
   items) with 2 workers and a 4 GB cap, then reanalyses into results.json and writes DONE2. Both
   numbers are reported: the first-pass analysis is kept as analyse_pass1.txt / results_pass1.json.
   No time limit, frame, condition or statistic changed.
2. (2026-11-07, 02:20 CET, launch check; a citation note, not a change.) PREREG.md cites amenability as
   "Perković et al. (2017, 2018)" without a definition number. The predicate the run computes is the submitted
   `is_amenable`, the first author's e10b rule; the definition it implements is b-amenability, Perković, Kalisch
   and Maathuis (UAI 2017), Definition 4.3 (i), with Theorem 4.4 for soundness and completeness, as E25's
   PREREG.md states and checks literally on the nearest failing states. Any sentence that uses E27 cites those
   numbers. The frozen PREREG.md is unchanged (sha256 4108c5a1…d39a9504 still matches the queue).
