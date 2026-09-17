import sys; sys.path.insert(0,'.')
from audit_w0_rk import *
from audit_w0_4_widths import ols, Z95
rng = np.random.default_rng(31415)
probs,_ = gen(rng, 300, n=7, p=0.40, k_claims=3, n_false=1)
bad_opt=0; tot=0; nonamen=0; probs_nonamen=0
for pr in probs:
    x,y = pr['x'],pr['y']
    amen_all = True
    for d in pr['fibre']:
        tot += 1
        o = optimal_adjustment_set_dag(d,x,y)
        if not is_valid_adjustment_set_dag(d,x,y,o):
            bad_opt += 1
        if y not in d.descendants(x):
            nonamen += 1; amen_all = False
    if not amen_all: probs_nonamen += 1
print(f"DAGs in fibres examined                     : {tot}")
print(f"  O*(d) is NOT a valid adjustment set in d  : {bad_opt}  ({bad_opt/tot:.3f})")
print(f"  Y not a descendant of X in d              : {nonamen}  ({nonamen/tot:.3f})")
print(f"problems where SOME d in the fibre is non-amenable: {probs_nonamen}/{len(probs)}")
# repeat the coverage check restricted to fully-amenable problems
rng2 = np.random.default_rng(2718); ok=[]; okA=[]
for pr in probs:
    dag,x,y = pr['dag'],pr['x'],pr['y']
    sem = random_sem(dag,rng2); truth = sem.true_total_effect(x,y)
    sigma=sem.covariance(); nodes=list(dag.nodes); idx={v:i for i,v in enumerate(nodes)}
    data = rng2.multivariate_normal(np.zeros(len(nodes)), sigma, size=2000)
    o = optimal_adjustment_set_dag(dag,x,y)
    e,s = ols(data,idx,x,o,y)
    c = abs(e-truth) <= Z95*s
    ok.append(c)
    if is_valid_adjustment_set_dag(dag,x,y,o): okA.append(c)
print(f"\ntrue DAG's own O* CI covers truth, ALL problems      : {np.mean(ok):.3f} (n={len(ok)})")
print(f"true DAG's own O* CI covers truth, O* actually VALID : {np.mean(okA):.3f} (n={len(okA)})")
