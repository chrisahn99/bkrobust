import collections, hashlib, json, os, sys, time
from pathlib import Path
BK = Path(os.path.expanduser("~/bkrobust")); sys.path.insert(0, str(BK/"src"))
from bkrobust.benchmarks.describe import parse_file
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.example import dag_to_cpdag, knowledge_to_recover
from bkrobust.demo.meek import apply_orientations
from bkrobust.demo.evaluate import optimal_adjustment_set_dag
from bkrobust.gac.mpdag_level import is_gac_valid_mpdag
SCR = Path("/private/tmp/claude-501/-Users-josecosta-mugango/0c3a820e-3bf3-499f-9297-016689024c33/scratchpad")
MODELS = BK/"results/axisa3/networks/example_models"
SKIP = {"pathfinder"}
def sk(k_true, coverage):
    if coverage >= 1.0: return list(k_true)
    keep = round(len(k_true)*coverage)
    if keep <= 0: return []
    step = len(k_true)/keep
    return [k_true[min(len(k_true)-1, int(i*step))] for i in range(keep)]
def carrega(name):
    p = next((q for q in sorted(MODELS.iterdir()) if q.name.startswith(name)), None)
    pn = parse_file(p, hashlib.sha256(p.read_bytes()).hexdigest())
    dag = MPDAG(pn.nodes, directed=pn.edges); return dag, dag_to_cpdag(dag)
rows = [json.loads(l) for l in open(BK/"results/axisa3/instances.jsonl")]
adm = [r for r in rows if r.get("admissible") and r["network"] not in SKIP]
nets = sorted({r["network"] for r in adm})
G = {}; KT = {}
for n in nets:
    t=time.time(); G[n]=carrega(n); KT[n]=sorted(knowledge_to_recover(*G[n])); print("loaded", n, "|k_true|", len(KT[n]), round(time.time()-t,1), "s", flush=True)
nest = collections.Counter(); viol=[]
for n in nets:
    K = {c: set(sk(KT[n], c)) for c in (1.0, 0.5, 0.25)}
    nest["n"]+=1; nest["K.25<=K.5"] += K[0.25] <= K[0.5]; nest["K.5<=K1"] += K[0.5] <= K[1.0]
    if not K[0.25] <= K[0.5]: viol.append((n, len(K[1.0]), len(K[0.5]), len(K[0.25])))
print("NESTING:", dict(nest), "violators:", viol, flush=True)
cov1 = [r for r in adm if r["coverage"]==1.0]
t0=time.time(); res=[]
for i,r in enumerate(cov1):
    dag, cp = G[r["network"]]; K = KT[r["network"]]
    z = frozenset(optimal_adjustment_set_dag(dag, r["x"], r["y"]))
    nsf=0; nnone=0
    for k in K:
        g = apply_orientations(cp, [e for e in K if e != k])
        if g is None: nnone+=1; continue
        if not is_gac_valid_mpdag(g, r["x"], r["y"], z): nsf+=1
    res.append((r["network"], r["x"], r["y"], r["radius"], len(K), r["k_g0"], nsf, nnone))
    if i % 50 == 0: print("  ", i, "/", len(cov1), round(time.time()-t0,1), "s", flush=True)
print("coverage-1.0 instances (pathfinder excluded):", len(res))
print("r_claim = 1:", sum(1 for t in res if t[6]>0), " no single retraction invalidates:", sum(1 for t in res if t[6]==0), flush=True)
print("r_hop among r_claim=1:", sorted(collections.Counter(t[3] for t in res if t[6]>0).items()))
print("r_hop among r_claim>1:", sorted(collections.Counter(t[3] for t in res if t[6]==0).items()))
phi = sorted(t[6]/t[4] for t in res if t[4])
print("phi_1: min", round(phi[0],3), "median", round(phi[len(phi)//2],3), "max", round(phi[-1],3))
print("median |K| asserted", sorted(t[4] for t in res)[len(res)//2], "median k_g0", sorted(t[5] for t in res)[len(res)//2])
json.dump(res, open(SCR/"rclaim_cov1.json","w"))
print("DONE", flush=True)
