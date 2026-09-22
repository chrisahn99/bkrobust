"""Orchestrator's INDEPENDENT check of falsification trigger T5 and of the
achieved intensity coverage of each cross-arm process. Imports nothing from the
analysis harness, on purpose: it exists to be compared against it."""
import json, glob, os, collections, statistics
DONE={os.path.basename(f)[:-5] for f in glob.glob('results/axis_robustness_real/_done/*.json')}
cells=[]
for f in glob.glob('results/axis_robustness_real/shards/xarm__*.cells.jsonl'):
    sid=os.path.basename(f)[:-len('.cells.jsonl')]
    if sid in DONE:
        for l in open(f): cells.append(json.loads(l))
print("xarm cells read:", len(cells))
sha=collections.defaultdict(set)
for c in cells:
    if c.get('instance_sha256'): sha[(c['network'],c['x'],c['y'],c['n_tiers'])].add(c['instance_sha256'])
bad=[k for k,v in sha.items() if len(v)>1]
print(f"instances whose two arms disagree on instance_sha256: {len(bad)} of {len(sha)}")
pool=collections.defaultdict(lambda: [0,0,0,0])
for c in cells:
    if c['status']!='ok' or c.get('intensity_bin') is None: continue
    k=(c['network'],c['x'],c['y'],c['n_tiers'],c['arm'],c['intensity_bin'])
    p=pool[k]; p[0]+=c['n_draws']; p[1]+=c['n_contradictory']; p[2]+=c['n_survived']; p[3]+=c['n_eval']
hist=collections.defaultdict(collections.Counter)
for (net,x,y,nt,arm,b),p in pool.items(): hist[arm][b]+=1
bins=sorted({b for h in hist.values() for b in h})
print("\nachieved intensity coverage ((instance, bin) cells per arm)")
print(f"{'bin':>6s} {'xarm_tiered':>12s} {'xarm_flip':>10s}")
for b in bins: print(f"{b:6.2f} {hist['xarm_tiered'][b]:12d} {hist['xarm_flip'][b]:10d}")
byinst=collections.defaultdict(dict)
for (net,x,y,nt,arm,b),p in pool.items(): byinst[(net,x,y,nt,b)][arm]=p
print("\nT5: matched instances per intensity bin (both arms, n_eval >= 30)")
print(f"{'bin':>6s} {'matched':>8s} {'networks':>9s}  {'mean S_t - S_f':>15s}")
n_ok=0
for b in bins:
    rows=[]; nets=set()
    for (net,x,y,nt,bb),d in byinst.items():
        if bb!=b: continue
        t,f=d.get('xarm_tiered'),d.get('xarm_flip')
        if not t or not f or t[3]<30 or f[3]<30: continue
        rows.append(t[2]/t[3]-f[2]/f[3]); nets.add(net)
    if rows:
        if len(rows)>=15: n_ok+=1
        print(f"{b:6.2f} {len(rows):8d} {len(nets):9d}  {statistics.mean(rows):+15.4f}"
              + ("  <-- >=15 matched" if len(rows)>=15 else ""))
print(f"\nbins with >= 15 matched instances in both arms: {n_ok}")
print("T5 verdict:", "pairing workable (>=3 bins)" if n_ok>=3 else "pairing NOT workable (<3 bins)")
