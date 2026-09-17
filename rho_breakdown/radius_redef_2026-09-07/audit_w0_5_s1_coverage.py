import sys; sys.path.insert(0,'.')
from audit_w0_rk import *
from audit_w0_4_widths import ols, Z95
rng = np.random.default_rng(778)
probs,_ = gen(rng, 300, n=7, p=0.40, k_claims=3, n_false=1)
rng2 = np.random.default_rng(1000)
in_fibre=0; mstar=[]; hull_cov=[]; own_ok=[]
for pr in probs:
    dag,x,y,K = pr['dag'],pr['x'],pr['y'],pr['K']
    fib = pr['fibre']
    hit = [d for d in fib if d.directed_edges==dag.directed_edges]
    in_fibre += len(hit)>0
    ms = m_of(dag, K); mstar.append(ms)
    sem = random_sem(dag, rng2); truth = sem.true_total_effect(x,y)
    sigma=sem.covariance(); nodes=list(dag.nodes); idx={v:i for i,v in enumerate(nodes)}
    data = rng2.multivariate_normal(np.zeros(len(nodes)), sigma, size=2000)
    S=[d for d in fib if m_of(d,K)<=1]
    vals=[ols(data,idx,x,optimal_adjustment_set_dag(d,x,y),y) for d in S]
    lo=min(v[0] for v in vals); hi=max(v[0] for v in vals)
    hull_cov.append(lo<=truth<=hi)
    e,s = ols(data,idx,x,optimal_adjustment_set_dag(dag,x,y),y)
    own_ok.append(abs(e-truth) <= Z95*s)
print(f"true DAG present in enumerated fibre : {in_fibre}/{len(probs)}")
print(f"m(d*) distribution                   : {Counter(mstar)}")
print(f"truth inside RAW hull of S_1         : {np.mean(hull_cov):.3f}")
print(f"true DAG's own O* CI covers truth    : {np.mean(own_ok):.3f}  (should be ~0.95)")
