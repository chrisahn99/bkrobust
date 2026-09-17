import sys, numpy as np
sys.path.insert(0,'/private/tmp/claude-501/-Users-josecosta-mugango/eee5422c-2e7e-4a55-8d52-0664260abb06/scratchpad/ref')
sys.path.insert(0,'/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/code')
import indep as I
import x2lib as X
CASES=[("licensed",10442,7,2.0,3,(0,2)),("large",1847,23,2.0,3,(6,22)),
       ("licensed",2650,9,2.0,2,(8,3)),("large",3610,20,2.5,1,(17,4)),
       ("large",4802,20,2.0,3,(11,6)),("licensed",11746,9,2.0,4,(2,6))]
NEED={"licensed":3,"large":3,"k8":8}
for ens,seed,p,deg,kb_claim,stmt in CASES:
    sd=X.draw_scm(seed,p,deg,NEED[ens],"uniform")
    if sd is None: print(ens,seed,"draw None"); continue
    rng=sd["rng"]; kb=int(rng.integers(1,min(4,len(sd["K_all"]))+1))
    x,y=sd["x"],sd["y"]; Dm=sd["D"]
    K_base=sd["K_all"][:kb]
    dag=[(i,j) for i in range(p) for j in range(p) if Dm[i,j]==1 and Dm[j,i]==0]
    C=I.dag_to_cpdag(dag,p); Dg=I.G(p); Dg.d=set(dag)
    dd=I.dist(C,[x,y]); u,v=stmt
    s_true=(u,v) if Dm[u,v]==1 else (v,u); s_false=(s_true[1],s_true[0])
    Kp=[t for t in K_base if frozenset(t)!=frozenset((u,v))]+[s_false]
    G0=I.apply_K(C,[tuple(t) for t in K_base])
    print('=== %s seed=%d p=%d query (%d,%d) k_base=%d(claim %d)'%(ens,seed,p,x,y,kb,kb_claim))
    print('   stmt asserted (the MISSTATEMENT): %d->%d ; delta=(%s,%s) dmin=%s dmax=%s'%(
          s_false[0],s_false[1],dd[s_false[0]],dd[s_false[1]],min(dd[u],dd[v]),max(dd[u],dd[v])))
    print('   base K =',[tuple(t) for t in K_base],' all true in D:',all((a,b) in Dg.d for a,b in K_base))
    try:
        G1=I.apply_K(C,[tuple(t) for t in Kp])
    except I.Inconsistent as e:
        print('   perturbed inconsistent:',e); continue
    a0,a1=I.amenable(G0,x,y),I.amenable(G1,x,y)
    n0,n1=len(I.pcp(G0,x,y)),len(I.pcp(G1,x,y))
    print('   G0 amenable=%s (paths %d) O*=%s | G1 amenable=%s (paths %d)'%(
        a0,n0,sorted(I.Ostar(G0,x,y)) if a0 else None,a1,n1))
    if a0 and not a1 and n1>0:
        bad=[q for q in I.pcp(G1,x,y) if not G1.isdir(q[0],q[1])]
        print('   >>> IDENTIFICATION DESTROYED by a single dmin=%s misstatement. offending path %s'%(min(dd[u],dd[v]),bad[0]))
        print('   >>> edges newly oriented (statement + Meek cascade):',sorted(set(G1.d)-set(G0.d)),
              ' un-oriented:',sorted(set(G0.d)-set(G1.d)))
        print('   >>> distance of each newly-oriented edge endpoint:',[(e,(dd[e[0]],dd[e[1]])) for e in sorted(set(G1.d)-set(G0.d))])
