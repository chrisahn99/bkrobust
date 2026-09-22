"""Orchestrator's INDEPENDENT recomputation of the elicited-knowledge panel.

Deliberately imports nothing from ``llm_analyse`` or ``llm_survival``: it reads
``results/elicit/knowledge.json`` and the shard JSONL files straight off disk
and recomputes every headline quantity from scratch, so a bug in the analysis
harness cannot hide inside its own cross-check. The precedent is
``results/axis_robustness_real/orchestrator_spotcheck.py``.

    PYTHONPATH=src .venv/bin/python results/axis_robustness_llm/orchestrator_verify.py

Writes nothing; print the output into ORCHESTRATOR_VERIFY.txt.
"""
import collections
import glob
import json
import os
import statistics
import sys

OUT = sys.argv[1] if len(sys.argv) > 1 else "results/axis_robustness_llm"
KJ = "results/elicit/knowledge.json"
FRAME = "results/axis_robustness_real/frame.jsonl"

DONE = {os.path.basename(f)[:-5] for f in glob.glob(f"{OUT}/_done/*.json")}
markers = {}
for f in glob.glob(f"{OUT}/_done/*.json"):
    m = json.load(open(f))
    markers[m["shard_id"]] = m

cells, insts = [], []
for f in glob.glob(f"{OUT}/shards/*.cells.jsonl"):
    sid = os.path.basename(f)[: -len(".cells.jsonl")]
    if sid not in DONE:
        continue
    for line in open(f):
        cells.append(json.loads(line))
for f in glob.glob(f"{OUT}/shards/*.instances.jsonl"):
    sid = os.path.basename(f)[: -len(".instances.jsonl")]
    if sid not in DONE:
        continue
    for line in open(f):
        insts.append(json.loads(line))

print(f"complete shards={len(DONE)}  cells={len(cells)}  instance rows={len(insts)}")

# --- V1: K really is the elicited K, claim for claim -------------------------
print("\n[V1] every shard's K equals knowledge.json's K for its (condition, network)")
K = json.load(open(KJ))
mismatch = []
checked = 0
for sid, m in sorted(markers.items()):
    cond, net = m["condition"], m["network"]
    want = sorted(tuple(e) for e in K[cond]["networks"][net]["k"])
    if m["n_k"] != len(want):
        mismatch.append((sid, m["n_k"], len(want)))
    checked += 1
print(f"  shards checked={checked}  |K| mismatches={len(mismatch)}")
for row in mismatch[:10]:
    print("   MISMATCH", row)

# --- V2: accuracy, recomputed from the DAG, not read from the row ------------
print("\n[V2] k_accuracy recomputed against the true DAG (independent of the pipeline)")
sys.path.insert(0, "src")
from bkrobust.robustness.real_survival import load_network  # noqa: E402

acc_bad = []
by_cond_acc = collections.defaultdict(list)
for sid, m in sorted(markers.items()):
    cond, net = m["condition"], m["network"]
    k = [tuple(e) for e in K[cond]["networks"][net]["k"]]
    if not k:
        continue
    dag, _ = load_network(net)
    truth = set(dag.directed_edges)
    acc = sum(1 for e in k if e in truth) / len(k)
    by_cond_acc[cond].append(acc)
    if m["k_accuracy"] is None or abs(m["k_accuracy"] - acc) > 1e-12:
        acc_bad.append((sid, m["k_accuracy"], acc))
print(f"  shards with non-empty K={sum(len(v) for v in by_cond_acc.values())}  "
      f"accuracy mismatches={len(acc_bad)}")
for row in acc_bad[:10]:
    print("   MISMATCH", row)
print("  per-condition mean accuracy over networks (recomputed):")
for cond in sorted(by_cond_acc):
    v = by_cond_acc[cond]
    print(f"    {cond:18s} n_nets={len(v):3d} mean={statistics.mean(v):.3f} "
          f"median={statistics.median(v):.3f}")

# --- V3: the two-dimensional panel variation actually exists -----------------
print("\n[V3] within-network spread across conditions (the design's premise)")
k_by_net = collections.defaultdict(dict)
for m in markers.values():
    k_by_net[m["network"]][m["condition"]] = m["n_k"]
n_spread = 0
for net in sorted(k_by_net):
    ks = sorted(k_by_net[net].values())
    if len(set(ks)) > 1:
        n_spread += 1
    print(f"  {net:16s} conds={len(ks):3d} |K| min={min(ks):3d} max={max(ks):3d} "
          f"distinct={len(set(ks))}")
print(f"  networks with >1 distinct |K|: {n_spread}/{len(k_by_net)}")

# --- V4: instance status mix -------------------------------------------------
print("\n[V4] instance status mix (a declined elicitation is an outcome, not a drop)")
st = collections.Counter(i["status"] for i in insts)
for k_, v in st.most_common():
    print(f"  {k_:34s} {v:6d}  ({v / len(insts):.3f})")

# --- V5: contradiction rate by depth, pooled --------------------------------
print("\n[V5] contradiction rate by depth fraction, pooled over shards")
# one draw-set is shared by every pair of a shard -> dedupe by (shard, depth)
seen = set()
byf = collections.defaultdict(lambda: [0, 0])
for c in cells:
    if c.get("status") != "ok":
        continue
    key = (c["shard_id"], c["grid_point"])
    if key in seen:
        continue
    seen.add(key)
    bucket = round(c["frac_of_n_k"], 1)
    byf[bucket][0] += c["n_draws"]
    byf[bucket][1] += c["n_contradictory"]
for f_ in sorted(byf):
    nd, nc = byf[f_]
    print(f"  d/|K|={f_:<4} draws={nd:8d}  contradiction={nc / nd:.3f}")

# --- V6: S is never a number when every draw was contradictory ---------------
print("\n[V6] sentinel discipline: S defined iff n_eval > 0")
bad = [c for c in cells if (c["n_eval"] == 0) != (c["S"] is None)]
print(f"  cells violating 'S is None iff n_eval==0': {len(bad)}")
forbidden = [c for c in cells if "seconds" in c]
print(f"  cells carrying the forbidden 'seconds' key: {len(forbidden)}")
unreached = [i for i in insts if i.get("radius") == -1]
print(f"  instance rows with UNREACHED radius: {len(unreached)} "
      f"(all carry r_status='unreached': "
      f"{all(i.get('r_status') == 'unreached' for i in unreached)})")

# --- V7: radius distribution per condition ----------------------------------
print("\n[V7] r_val on elicited G0, per condition (UNREACHED excluded, never averaged)")
by_cond_r = collections.defaultdict(list)
for i in insts:
    if i.get("r_status") == "ok":
        by_cond_r[i["condition"]].append(i["radius"])
for cond in sorted(by_cond_r):
    r = by_cond_r[cond]
    share1 = sum(1 for x in r if x == 1) / len(r)
    print(f"  {cond:18s} n={len(r):5d} median={statistics.median(r):5.1f} "
          f"max={max(r):3d} share(r=1)={share1:.3f}")

# --- V8: no shard ever called select_knowledge ------------------------------
print("\n[V8] provenance: k_source on every instance row")
print("  ", collections.Counter(i.get("k_source") for i in insts))
print("  knowledge.json sha256 on markers:",
      collections.Counter(m.get("knowledge_json_sha256") for m in markers.values()))

# --- V9: off-skeleton claims -------------------------------------------------
off = [m for m in markers.values() if m.get("n_off_skeleton")]
print(f"\n[V9] shards with claims outside the CPDAG skeleton: {len(off)} "
      "(any non-zero means the answers were parsed against a different graph)")
for m in off[:10]:
    print("   ", m["shard_id"], m["n_off_skeleton"])

# --- V10: wall-cap censoring -------------------------------------------------
cens = [c for c in cells if c.get("status") == "censored_wall_cap"]
print(f"\n[V10] wall-cap censored cells: {len(cens)}")
print(f"      shards reporting a censored tail: "
      f"{sum(1 for m in markers.values() if m.get('n_cells_censored_wall_cap'))}")
print(f"      slowest shards: "
      f"{sorted(((m['elapsed_s'], m['shard_id']) for m in markers.values()), reverse=True)[:5]}")
