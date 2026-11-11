"""E2 summary: per n (pooled over sparse families) and query class, the median,
95th percentile and maximum seconds, the leg that answered, and the size of the
search problem (undirected edges, |K_G0|)."""
import json
import pandas as pd
from common import OUT
d = pd.read_csv(OUT / "e2" / "timing.csv")
d = d[d.method != "no_query"].copy()
d["fam"] = d.family.where(d.family == "dense_appD", "sparse")
g = d.groupby(["fam", "n", "query_class"])
s = g.agg(N=("total_s", "size"), median_s=("total_s", "median"),
          p95_s=("total_s", lambda x: x.quantile(0.95)), max_s=("total_s", "max"),
          top=("method", lambda x: (x == "top_state").sum()),
          search=("method", lambda x: (x == "local_up_fast").sum()),
          sat=("method", lambda x: (x == "e1_ladder").sum()),
          k_und_med=("k_undirected", "median"), k_g0_med=("k_g0", "median"),
          validity_call_med=("validity_call_s", "median"), build_med=("build_s", "median"))
print(s.round(5).to_string())
tot = {"n_queries": int(len(d)), "censored": int(d.censored.astype(str).eq("True").sum()),
       "method_counts": d.method.value_counts().to_dict(),
       "share_under_1s": float((d.total_s < 1).mean()),
       "share_under_0.1s": float((d.total_s < 0.1).mean()),
       "max_s": float(d.total_s.max()),
       "sat_rows": d[d.method == "e1_ladder"][["family", "n", "knowledge", "k_g0", "r", "total_s", "ladder_s"]].to_dict(orient="records"),
       "slowest": d.nlargest(8, "total_s")[["family", "n", "knowledge", "query_class", "k_g0", "r", "method", "total_s", "search_s", "ladder_s"]].to_dict(orient="records")}
(OUT / "e2" / "summary.json").write_text(json.dumps({"by_cell": s.reset_index().to_dict(orient="records"), **tot}, indent=1, default=str))
print(json.dumps({k: v for k, v in tot.items() if k != "by_cell"}, indent=1, default=str))
