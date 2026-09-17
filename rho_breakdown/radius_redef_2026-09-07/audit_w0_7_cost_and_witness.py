import sys; sys.path.insert(0,'.')
from audit_w0_rk import *
from bkrobust.demo.scenario import true_dag

# --- (a) the cost claim, on the 8-variable demo -------------------------------
d = true_dag(); cp = dag_to_cpdag(d)
und = sorted(tuple(sorted(e)) for e in cp.undirected_edges)
K = [ (a,b) if (a,b) in d.directed_edges else (b,a) for (a,b) in und[:3] ]
g0 = apply_orientations(cp, K); nodes=sorted(d.nodes); x,y = nodes[0], nodes[-1]
for xx in nodes:
    for yy in nodes:
        if xx!=yy and optimal_adjustment_set_mpdag(g0,xx,yy) is not None:
            x,y = xx,yy; break
    else: continue
    break
z = frozenset(optimal_adjustment_set_mpdag(g0,x,y))
fib = enumerate_dag_extensions(cp)
t0=time.perf_counter()
for _ in range(20): r_K(fib,K,x,y,z)
t_rk=(time.perf_counter()-t0)/20
t0=time.perf_counter()
for _ in range(3): r_val(cp,g0,x,y,z)
t_rv=(time.perf_counter()-t0)/3
print(f"demo graph: |fibre|={len(fib)} |space|={len(enumerate_space(cp))}")
print(f"  r_K   {t_rk*1000:8.3f} ms   (excludes the enumerate_dag_extensions call itself)")
t0=time.perf_counter(); enumerate_dag_extensions(cp); t_enum=time.perf_counter()-t0
print(f"  + fibre enumeration {t_enum*1000:.3f} ms  -> honest total {(t_rk+t_enum)*1000:.3f} ms")
print(f"  r_val {t_rv*1000:8.3f} ms   -> speedup {(t_rv)/(t_rk+t_enum):.1f}x  (claimed ~200x)")

# --- (b) the witness set: direction (4) for free ------------------------------
rng = np.random.default_rng(5150)
probs,_ = gen(rng, 400, n=7, p=0.40, k_claims=4, n_false=0)
sizes=[]; unions=[]; nK=[]
for pr in probs:
    K=pr['K']; bad=[(m_of(d,K),d) for d in pr['fibre']
                    if not is_valid_adjustment_set_dag(d,pr['x'],pr['y'],pr['z'])]
    if not bad: continue
    r=min(b[0] for b in bad)
    wit=[frozenset((a,b) for (a,b) in K if d.is_directed_edge(b,a)) for m,d in bad if m==r]
    u=set().union(*wit)
    sizes.append(len(set(wit))); unions.append(len(u)); nK.append(len(K))
print(f"\nfinite-r_K problems: {len(sizes)}")
print(f"  distinct minimum-cardinality witness claim-sets : mean {np.mean(sizes):.2f}  max {max(sizes)}")
print(f"  |union of witnesses| / |K|                      : mean {np.mean(np.array(unions)/np.array(nK)):.3f}")
print(f"  |K| mean {np.mean(nK):.2f}, flagged claims mean {np.mean(unions):.2f}")
print(f"  problems where the witness union is a STRICT subset of K: {np.mean(np.array(unions)<np.array(nK)):.3f}")
