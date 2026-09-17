import sys, json, numpy as np
sys.path.insert(0,'/private/tmp/claude-501/-Users-josecosta-mugango/eee5422c-2e7e-4a55-8d52-0664260abb06/scratchpad/ref')
sys.path.insert(0,'/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/code')
import indep as I
from graphs import random_dag
from scm import make_linear_iscm

ex=json.load(open('/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/results/exhibit.json'))
for k,dos in enumerate(ex['dossiers']):
    p=dos['p']; seed=dos['seed']; deg=dos['deg']
    rng=np.random.default_rng(seed)
    Dm=random_dag(p,deg,rng)
    dag_edges=[(i,j) for i in range(p) for j in range(p) if Dm[i,j]==1 and Dm[j,i]==0]
    claimed=[tuple(int(t) for t in e.split('->')) for e in dos['D_edges']]
    print('=== dossier',k,'ens',dos['ens'],dos['arm'],'seed',seed,'p',p)
    print(' D regen matches claim:', sorted(dag_edges)==sorted(claimed))
    Dg=I.G(p); Dg.d=set(dag_edges)
    C=I.dag_to_cpdag(dag_edges,p)
    print(' C (mine):',C.edges())
    print(' C matches claim:', sorted(C.edges())==sorted(dos['C_edges']))
    x,y=dos['x'],dos['y']
    dd=I.dist(C,[x,y]); print(' delta (mine):',dd, 'claimed',dos['delta'])
    K=[tuple(e) for e in dos['K']]; Kp=[tuple(e) for e in dos['K_perturbed']]
    # check all base K statements are TRUE in D
    print(' all base K true in D:', all((i,j) in Dg.d for (i,j) in K))
    # check flipped set
    flipped=[tuple(e) for e in dos['flipped']]
    print(' flipped:',flipped,' dmin each:',[min(dd[a],dd[b]) for a,b in flipped],
          ' dmax each:',[max(dd[a],dd[b]) for a,b in flipped])
    try:
        G0=I.apply_K(C,K); ok0=True
    except I.Inconsistent as e: ok0=False; print(' BASE INCONSISTENT',e)
    try:
        G1=I.apply_K(C,Kp); ok1=True
    except I.Inconsistent as e: ok1=False; print(' PERT INCONSISTENT',e)
    if not (ok0 and ok1): continue
    print(' G0 (mine):',G0.edges()); print(' G0 matches claim:', sorted(G0.edges())==sorted(dos['G0_edges']))
    print(' G1 (mine):',G1.edges()); print(' G1 matches claim:', sorted(G1.edges())==sorted(dos['G1_edges']))
    print(' amenable G0/G1:',I.amenable(G0,x,y), I.amenable(G1,x,y))
    O0=I.Ostar(G0,x,y); O1=I.Ostar(G1,x,y)
    print(' O*(G0)=',sorted(O0) if O0 is not None else None,'claim',dos['O0'])
    print(' O*(G1)=',sorted(O1) if O1 is not None else None,'claim',dos['O1'])
    print(' O* CHANGED:', O0!=O1)
    # SCM regen: draw_scm consumes rng: random_dag, then pairs, then integers, then make_linear_iscm
    from adjust import possibly_causal_paths
    pairs_D=[(a,b) for a in range(p) for b in range(p) if a!=b and len(possibly_causal_paths(Dm,a,b))>0]
    xi,yi=pairs_D[int(rng.integers(len(pairs_D)))]
    scm=make_linear_iscm(Dm,rng)
    A=scm['A']; S=I.sigma_of(A,scm['omega'])
    tau=I.total_effect(A,x,y)
    print(' regen query (x,y)=',(xi,yi),'claimed',(x,y),' match:',(xi,yi)==(x,y))
    print(' tau mine %.6f claimed %.6f'%(tau,dos['tau']))
    b0=I.ols(S,x,y,O0); b1=I.ols(S,x,y,O1)
    print(' beta0 mine %.6f claimed %.6f | beta1 mine %.6f claimed %.6f'%(b0,dos['beta0'],b1,dos['beta1']))
    print(' |beta1-tau|/|tau| = %.4f  (claimed rel_bias %.4f)'%(abs(b1-tau)/abs(tau),dos['rel_bias']))
    # independent validity: brute force over all subsets? use OLS==tau as ground truth
    print(' O0 recovers tau:',abs(b0-tau)<1e-9,' O1 recovers tau:',abs(b1-tau)<1e-9)
