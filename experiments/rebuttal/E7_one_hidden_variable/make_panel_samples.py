"""Write the panel's own n = 5000 samples for the cells that do not regenerate bit for bit elsewhere.

Run it on a machine where ``sample_for`` reproduces the hash stored in structure.csv (the
panel was made on the first author's laptop: arm64, numpy 2.5.3). It refuses to write a
sample whose hash differs from the stored one, so a file in panel_samples/ is always the
panel's sample. run.py loads a file only when regeneration fails and the file's hash matches.

    python experiments/rebuttal/E7_one_hidden_variable/make_panel_samples.py magic-niab:2 magic-irri:3
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
for _p in (ROOT / "src", HERE):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import numpy as np  # noqa: E402

from run import PANEL_N, PANEL_SAMPLES, panel_sample_path  # noqa: E402
from vendor_pest import data_sha, load_gaussian, sample_for  # noqa: E402


def main(argv: list[str]) -> int:
    pest = ROOT / "experiments" / "rebuttal" / "_data" / "pest" / "structure.csv"
    models = ROOT / "results" / "axisa3" / "networks" / "example_models"
    stored = {
        (r["network"], int(r["seed"])): r["sample_sha256"]
        for r in csv.DictReader(pest.open())
        if r["method"] == "pc" and int(r["n"]) == PANEL_N and r["status"] == "ok"
    }
    if not argv:
        print(__doc__, file=sys.stderr)
        return 2
    PANEL_SAMPLES.mkdir(exist_ok=True)
    bad = 0
    for cell in argv:
        net, seed = cell.split(":")
        dag, weights, noise, _ = load_gaussian(models, net)
        data = sample_for(dag, weights, noise, net, PANEL_N, int(seed))
        h, want = data_sha(data), stored[(net, int(seed))]
        if h != want:
            print(f"{net} seed {seed}: this machine gives {h}, panel {want}; not written", file=sys.stderr)
            bad += 1
            continue
        np.save(panel_sample_path(net, int(seed)), np.ascontiguousarray(data), allow_pickle=False)
        print(f"{net} seed {seed}: {h} written to {panel_sample_path(net, int(seed))}")
    return 3 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
