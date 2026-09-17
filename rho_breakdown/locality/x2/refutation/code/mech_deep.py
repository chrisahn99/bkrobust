"""WHY is the far zero true?  Measure whether the Meek cascade from a FAR
statement ever reaches an edge incident to cn(X,Y) -- the only edges O* can
depend on.  If it never does, the zero is a geometric fact about cn, not an
accident of sample size."""
import sys, itertools, collections, numpy as np
sys.path.insert(0,'/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/code')
from graphs import dag_to_cpdag, undirected_edges, apply_background_knowledge, MeekFail, skeleton, directed_edges
from adjust import optimal_adjustment_set, possibly_causal_paths, causal_nodes, is_amenable
INF=1<<20
def hop(C,srcs,p):
    S=skeleton(C); d=np.full(p,INF,int)
    fr=list(srcs)
    for v in fr: d[v]=0
    while fr:
        nx=[]
        for v in fr:
            for w in np.flatnonzero(S[v]):
                if d[w]>d[v]+1: d[w]=d[v]+1; nx.append(int(w))
        fr=nx
    return d
def clique_chain(rng,ncl,s):
    p=ncl*s; S=np.zeros((p,p),np.int8)
    for c in range(ncl):
        idx=list(range(c*s,(c+1)*s))
        for a,b in itertools.combinations(idx,2): S[a,b]=S[b,a]=1
        if c:
            u=rng.integers((c-1)*s,c*s); v=rng.integers(c*s,(c+1)*s); S[u,v]=S[v,u]=1
    return S
def dag_from(S,rng):
    p=S.shape[0]; o=rng.permutation(p); pos=np.empty(p,int)
    for i,v in enumerate(o): pos[v]=i
    D=np.zeros((p,p),np.int8)
    for a in range(p):
        for b in range(a+1,p):
            if S[a,b]:
                if pos[a]<pos[b]: D[a,b]=1
                else: D[b,a]=1
    return D
tot=collections.Counter()
rng0=np.random.default_rng(12345)
NG=2500
for g in range(NG):
    rng=np.random.default_rng(int(rng0.integers(1<<40)))
    S=clique_chain(rng,int(rng.integers(3,5)),int(rng.integers(3,5))); p=S.shape[0]
    D=dag_from(S,rng); C=dag_to_cpdag(D); U=undirected_edges(C)
    if len(U)<2: continue
    pairs=[(a,b) for a in range(p) for b in range(p) if a!=b and len(possibly_causal_paths(D,a,b))>0]
    if not pairs: continue
    best,bd=[],-1;
    for (aa,bb) in pairs:
        dq=hop(C,[aa],p)[bb]
        if dq<INF and dq>bd: best,bd=[(aa,bb)],dq
        elif dq==bd: best.append((aa,bb))
    x,y=best[int(rng.integers(len(best)))]
    dist=hop(C,[x,y],p)
    Kall=[(i,j) if D[i,j]==1 else (j,i) for (i,j) in U]
    idx=rng.permutation(len(Kall)); kb=int(rng.integers(0,len(Kall)+1))
    Kb=[Kall[i] for i in idx[:kb]]
    try: G0=apply_background_knowledge(C,Kb)
    except MeekFail: continue
    if not is_amenable(G0,x,y): continue
    cn0=causal_nodes(G0,x,y); O0=optimal_adjustment_set(G0,x,y)
    # undirected edges of G0 incident to cn
    U0=undirected_edges(G0)
    tot['graphs']+=1
    tot['u_inc_cn']+= sum(1 for (a,b) in U0 if a in cn0 or b in cn0)
    tot['has_u_inc_cn']+= int(any(a in cn0 or b in cn0 for (a,b) in U0))
    tot['cn_size']+=len(cn0)
    tot['max_cn_dist']+= max((int(dist[c]) for c in cn0 if dist[c]<INF), default=0)
    for (u,v) in U0:
        dm=int(min(dist[u],dist[v]))
        far = 2<=dm<INF
        if not far: continue
        for (a,b) in ((u,v),(v,u)):
            try: G=apply_background_knowledge(G0,[(a,b)])
            except MeekFail: continue
            tot['far_trials']+=1
            newdir=set(map(tuple,directed_edges(G)))-set(map(tuple,directed_edges(G0)))
            touch = any((i in cn0 or j in cn0) for (i,j) in newdir)
            tot['far_cascade_pos']+= int(len(newdir)>1)
            tot['far_touch_cn']+= int(touch)
            cn1=causal_nodes(G,x,y)
            tot['far_cn_changed']+= int(cn1!=cn0)
            O=optimal_adjustment_set(G,x,y)
            tot['far_amen']+= int(O is not None)
            if O is not None: tot['far_Ostar_changed']+= int(O!=O0)
print(f"graphs kept                             : {tot['graphs']}")
print(f"mean |cn|                               : {tot['cn_size']/tot['graphs']:.3f}")
print(f"mean max skeleton distance of a cn node : {tot['max_cn_dist']/tot['graphs']:.3f}")
print(f"mean undirected edges incident to cn    : {tot['u_inc_cn']/tot['graphs']:.4f}")
print(f"P(>=1 undirected edge incident to cn)   : {tot['has_u_inc_cn']/tot['graphs']:.4f}")
print(f"far trials (dmin>=2)                    : {tot['far_trials']}")
print(f"  of which Meek cascaded (>1 new edge)  : {tot['far_cascade_pos']} ({tot['far_cascade_pos']/max(tot['far_trials'],1):.3f})")
print(f"  cascade oriented an edge INCIDENT to cn: {tot['far_touch_cn']} ({tot['far_touch_cn']/max(tot['far_trials'],1):.5f})")
print(f"  cn(X,Y) itself changed                 : {tot['far_cn_changed']} ({tot['far_cn_changed']/max(tot['far_trials'],1):.5f})")
print(f"  still amenable                         : {tot['far_amen']}")
print(f"  O* changed                             : {tot['far_Ostar_changed']}")
