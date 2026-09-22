"""Orchestrator's INDEPENDENT recomputation of key quantities, straight from the
shard files. Deliberately imports nothing from real_analyse. Used to cross-check
the analysis harness rather than to trust it."""
import json, glob, os, collections, statistics

DONE = {os.path.basename(f)[:-5] for f in glob.glob('results/axis_robustness_real/_done/*.json')}
cells, insts = [], []
for f in glob.glob('results/axis_robustness_real/shards/*.cells.jsonl'):
    sid = os.path.basename(f)[:-len('.cells.jsonl')]
    if sid not in DONE: continue
    for l in open(f): cells.append(json.loads(l))
for f in glob.glob('results/axis_robustness_real/shards/*.instances.jsonl'):
    sid = os.path.basename(f)[:-len('.instances.jsonl')]
    if sid not in DONE: continue
    for l in open(f): insts.append(json.loads(l))
print(f"complete shards={len(DONE)}  cells={len(cells)}  instance rows={len(insts)}")

# --- P5: flip-arm contradiction rate at depth 1, pooled, by base wrongness
print("\n[P5] flip-arm contradiction rate at depth d=1 (pooled draws, not a mean of rates)")
for bw in (0.0, 0.10, 0.25, None):
    sel = [c for c in cells if c['arm']=='flip' and c['grid_point']==1
           and c['base_wrongness']==bw and c['status']=='ok']
    # one draw-set is shared by all pairs of a network -> dedupe by (shard_id, grid_point)
    seen=set(); nd=nc=0
    for c in sel:
        k=(c['shard_id'], c['grid_point'])
        if k in seen: continue
        seen.add(k); nd+=c['n_draws']; nc+=c['n_contradictory']
    label = f"bw={bw}" if bw is not None else "bw_abs=1"
    print(f"  {label:12s} networks-cells={len(seen):3d}  draws={nd:7d}  contradiction={nc/nd if nd else float('nan'):.3f}")

print("\n[P5b] flip-arm contradiction rate by depth, bw=0.00 only (pooled)")
byd=collections.defaultdict(lambda:[0,0])
seen=set()
for c in cells:
    if c['arm']!='flip' or c['base_wrongness']!=0.0 or c['status']!='ok': continue
    k=(c['shard_id'],c['grid_point'])
    if k in seen: continue
    seen.add(k)
    byd[c['grid_point']][0]+=c['n_draws']; byd[c['grid_point']][1]+=c['n_contradictory']
for d in sorted(byd)[:14]:
    nd,nc=byd[d]; print(f"  d={d:3d}  draws={nd:7d}  contradiction={nc/nd:.3f}")

# --- tiered contradiction and S by rate (pooled)
print("\n[tiered] by corruption_rate, pooled over all n_tiers")
byr=collections.defaultdict(lambda:[0,0,0,0])
seen=set()
for c in cells:
    if c['arm']!='tiered' or c['status']!='ok': continue
    k=(c['shard_id'],c['grid_point'])
    if k not in seen:
        seen.add(k); byr[c['grid_point']][0]+=c['n_draws']; byr[c['grid_point']][1]+=c['n_contradictory']
    byr[c['grid_point']][2]+=c['n_survived']; byr[c['grid_point']][3]+=c['n_eval']
for r in sorted(byr):
    nd,nc,ns,ne=byr[r]
    print(f"  rate={r:<5} draws={nd:8d} contra={nc/nd if nd else 0:.3f}  S={ns/ne if ne else float('nan'):.3f}")

# --- T4 again over the full data
t0=[c for c in cells if c['arm']=='tiered' and c['grid_point']==0.0 and c['status']=='ok']
bad=[c for c in t0 if c['S']!=1.0]
print(f"\n[T4] tiered rate-0 cells={len(t0)}  S!=1.000 exactly: {len(bad)}")

# --- T6 over the full data
worst=0.0; n=0
for c in cells:
    if c['S'] is None or c['contradiction_rate'] is None: continue
    worst=max(worst, abs(c['S_contra_as_fail']-(1-c['contradiction_rate'])*c['S'])); n+=1
print(f"[T6] decomposition identity on {n} cells: max dev {worst:.3e}")

# --- statuses and n_draws
print("\n[status] instance:", collections.Counter(i['status'] for i in insts).most_common())
print("[status] cell    :", collections.Counter(c['status'] for c in cells).most_common())
print("[N] n_draws values on non-censored cells:", collections.Counter(c['n_draws'] for c in cells if c['status']=='ok'))
print("[sentinel] rows with a 'seconds' key:", sum(1 for r in insts+cells if 'seconds' in r))
print("[sentinel] UNREACHED instances:", sum(1 for i in insts if i.get('radius')==-1))
print("[sentinel] censored radii:", sum(1 for i in insts if i.get('r_status')=='censored_wall_cap'))
