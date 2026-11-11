"""E11: adapt the frozen CPDAG-only frame (results/frame/frame.jsonl) to the row
format of the submitted survival runner. Networks above 500 nodes (link, munin,
munin2-4; see E12) and pathfinder (absent from the elicitation bundle) are left
out, which leaves the 528 queries per condition of the paper's release.
Separation and the programmatic-K radius are not part of this frame (None)."""
import hashlib, json, sys
from pathlib import Path
from common import ROOT
from bkrobust.robustness import real_survival as rs
from bkrobust.demo.graph import undirected_components

SKIP = {"link", "munin", "munin2", "munin3", "munin4", "pathfinder"}
out = ROOT / "results/rebuttal/e11_frame"
out.mkdir(parents=True, exist_ok=True)
rows, lcs = [], {}
for line in open(ROOT / "results/frame/frame.jsonl"):
    r = json.loads(line)
    if r["network"] in SKIP:
        continue
    if r["network"] not in lcs:
        _, cp = rs.load_network(r["network"])
        lcs[r["network"]] = max((len(c) for c in undirected_components(cp)), default=0)
    rid = hashlib.sha256(f"{r['network']}|{r['x']}|{r['y']}".encode()).hexdigest()[:16]
    rows.append({"row_id": rid, "network": r["network"], "x": r["x"], "y": r["y"],
                 "largest_component_size": lcs[r["network"]],
                 "component_size_recomputed": r["component_size"],
                 "separation_recomputed": None, "separation_status_recomputed": "not_computed",
                 "radius_recomputed": None, "frame_source": r["frame_source"],
                 "in_committed_admissible": r["in_committed_admissible"]})
txt = "".join(json.dumps(r) + "\n" for r in rows)
(out / "frame.jsonl").write_text(txt)
(out / "frame_hash.json").write_text(json.dumps({"frame_sha256": hashlib.sha256(txt.encode()).hexdigest()}))
print(len(rows), "rows,", len(lcs), "networks")
