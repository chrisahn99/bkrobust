"""E10 summary: radius spread across valid sets, how often O is beaten."""
import json
import pandas as pd
from common import OUT
d = pd.read_csv(OUT / "e10" / "valid_sets.csv", low_memory=False)
out = {"status": d.status.value_counts().to_dict()}
ok = d[d.status == "ok"].copy()
for c in ("varies", "O_beaten", "O_is_max", "all_exact", "O_valid_truth"):
    ok[c] = ok[c].astype(str) == "True"
multi = ok[ok.n_sets >= 2]
out.update({
    "n_queries": int(len(ok)), "n_queries_2plus_sets": int(len(multi)),
    "mode": ok["mode"].value_counts().to_dict(),
    "median_sets": float(ok.n_sets.median()), "max_sets": int(ok.n_sets.max()),
    "share_radius_varies_2plus": float(multi.varies.mean()),
    "share_O_beaten_2plus": float(multi.O_beaten.mean()),
    "n_O_beaten": int(multi.O_beaten.sum()),
    "share_O_is_max_2plus": float(multi.O_is_max.mean()),
    "all_exact": bool(ok.all_exact.all()),
    "O_beaten_by_rO": multi.groupby(multi.r_O.astype(str)).O_beaten.agg(["sum", "count"]).to_dict(),
    "varies_by_rO": multi.groupby(multi.r_O.astype(str)).varies.agg(["sum", "count"]).to_dict(),
    "O_beaten_exhaustive": int(multi[multi["mode"] == "exhaustive"].O_beaten.sum()),
    "n_exhaustive_2plus": int((multi["mode"] == "exhaustive").sum()),
    "O_beaten_examples": multi[multi.O_beaten][["arm", "network", "x", "y", "mode", "n_sets", "r_O", "r_max", "O_valid_truth", "maxr_set_valid_truth"]].head(15).to_dict(orient="records"),
    "O_beaten_by_arm": multi.groupby("arm").O_beaten.sum().to_dict(),
    "O_beaten_maxset_valid_truth": multi[multi.O_beaten].maxr_set_valid_truth.astype(str).value_counts().to_dict(),
    "O_beaten_O_valid_truth": multi[multi.O_beaten].O_valid_truth.value_counts().to_dict(),
})
(OUT / "e10" / "summary.json").write_text(json.dumps(out, indent=1, default=str))
print(json.dumps(out, indent=1, default=str))
