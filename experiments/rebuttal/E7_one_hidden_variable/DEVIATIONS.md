# E7 deviations

## 2026-11-06: two panel samples are loaded from a file instead of regenerated

The first full launch on the workstation (results/rebuttal/E7_failed_run1) stopped at consistency check (1) after seven units: the nothing-hidden unit of magic-niab seed 2 regenerated a sample with hash 059f97046f5ded15, where structure.csv stores 603fc76330509e1f (structure.csv, row `magic-niab,5000,2,pc`).

The retry file graphs_retry.jsonl is not the cause. It holds only two arth150 GES graphs (n = 1000 and n = 5000, seed 1) re-run after a timeout, and the panel's sample code (pest_panel.py, `sample_for` and `structure_task`) has no retry seed: the sample hash is taken before discovery runs, from the seeds [20260912, stable_seed(net), n, seed] alone.

The cause is the platform. The panel was made on the first author's laptop (arm64, numpy 2.5.3); the workstation is x86_64 with numpy 2.5.4. On the laptop the same code gives the stored hash for magic-niab seed 2, magic-irri seed 3 and the two controls tried. On the workstation, regenerating all 120 panel samples (4 networks, n in {200, 1000, 5000}, seeds 0-9) gives 106 matching hashes and 14 that differ. For magic-niab seed 2, exactly one of the 220,000 raw standard-normal draws differs, by one ulp (4.4e-16, the tail value -3.7988 at index 83333), which most likely comes from the platform's math library in numpy's normal sampler; this reaches 5 of the 220,000 sample entries (largest difference 8.9e-16). For magic-irri seed 3, 6 of 320,000 entries differ, again by at most 8.9e-16. PC on the the workstation sample of magic-niab seed 2 still learns the stored panel graph (hash 3dfc6c3b437447c0).

Within E7's grid (n = 5000, seeds 0-4), two cells do not regenerate on the workstation: magic-niab seed 2 and magic-irri seed 3. The other 18 do.

What changed. make_panel_samples.py writes the panel's own sample for a cell to panel_samples/, and refuses unless the machine running it reproduces the stored hash; it was run on the laptop for those two cells. run.py first regenerates each sample as before; only when the hash differs does it load the file from panel_samples/, and it keeps the file only when its hash equals the stored one. Before any PC runs, run.py now checks every (network, seed) cell of the grid this way and stops with exit 3 if a cell matches neither; `--check-samples` runs that check alone. Each unit line records `sample_source` (regenerated or vendored) and the regenerated hash, and config.json records the source of every cell and the sha256 of each file used. No sample is excluded. The population, seeds, method and statistic of PREREG.md are unchanged, and check (1) is kept in its strict form: every unit uses a sample whose hash equals the panel's.

Both runs are reported, as PREREG.md requires: the failed launch is kept in results/rebuttal/E7_failed_run1 and was not deleted.
