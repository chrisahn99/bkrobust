"""INDEPENDENT re-implementation. Nothing imported from x2/e1prime code except
scm.make_linear_iscm / graphs.random_dag for regenerating the SCM draw."""
import itertools, numpy as np, networkx as nx

# graph repr: dict of sets. dir[(i,j)] means i->j ; und = set of frozensets
class G:
    def __init__(self, p):
        self.p=p; self.d=set(); self.u=set()
    def copy(self):
        g=G(self.p); g.d=set(self.d); g.u=set(self.u); return g
    def adj(self,i,j): return (i,j) in self.d or (j,i) in self.d or frozenset((i,j)) in self.u
    def isdir(self,i,j): return (i,j) in self.d
    def isund(self,i,j): return frozenset((i,j)) in self.u
    def orient(self,i,j):
        self.u.discard(frozenset((i,j))); self.d.add((i,j))
    def pa(self,j): return {i for i in range(self.p) if (i,j) in self.d}
    def ch(self,i): return {j for j in range(self.p) if (i,j) in self.d}
    def nb(self,i): return {j for j in range(self.p) if frozenset((i,j)) in self.u}
    def pdn(self,i):  # possibly-directed neighbours
        return self.ch(i)|self.nb(i)
    def edges(self):
        return sorted([('%d->%d'%e) for e in self.d]+['%d--%d'%tuple(sorted(e)) for e in self.u])

def meek(g):
    g=g.copy(); ch=True
    while ch:
        ch=False
        for e in list(g.u):
            a,b=sorted(e)
            for (a,b) in ((a,b),(b,a)):
                if not g.isund(a,b): continue
                fire=False
                # R1 c->a, a--b, c not adj b
                for c in range(g.p):
                    if c in (a,b): continue
                    if g.isdir(c,a) and not g.adj(c,b): fire=True;break
                # R2 a->c->b
                if not fire:
                    for c in range(g.p):
                        if c in (a,b): continue
                        if g.isdir(a,c) and g.isdir(c,b): fire=True;break
                # R3 a--c,a--d,c->b,d->b, c,d nonadj
                if not fire:
                    cand=[c for c in range(g.p) if c not in (a,b) and g.isund(a,c) and g.isdir(c,b)]
                    for c,d in itertools.combinations(cand,2):
                        if not g.adj(c,d): fire=True;break
                # R4 a--d, d->c, c->b, a adj c, b,d nonadj
                if not fire:
                    for d in range(g.p):
                        if d in (a,b) or not g.isund(a,d) or g.adj(b,d): continue
                        for c in range(g.p):
                            if c in (a,b,d): continue
                            if g.isdir(d,c) and g.isdir(c,b) and g.adj(a,c): fire=True;break
                        if fire: break
                if fire:
                    g.orient(a,b); ch=True; break
    return g

def vstructs(g):
    out=set()
    for b in range(g.p):
        P=sorted(g.pa(b))
        for a,c in itertools.combinations(P,2):
            if not g.adj(a,c): out.add((a,b,c))
    return out

def has_cycle(g):
    G2=nx.DiGraph(); G2.add_nodes_from(range(g.p)); G2.add_edges_from(g.d)
    return not nx.is_directed_acyclic_graph(G2)

def dag_to_cpdag(dag_edges,p):
    g=G(p)
    for (i,j) in dag_edges: g.u.add(frozenset((i,j)))
    D=G(p); D.d=set(dag_edges)
    for (a,b,c) in vstructs(D):
        g.u.discard(frozenset((a,b))); g.u.discard(frozenset((c,b)))
        g.d.add((a,b)); g.d.add((c,b))
    return meek(g)

class Inconsistent(Exception): pass

def apply_K(C,K):
    g=C.copy()
    ref=vstructs(C)
    for (i,j) in K:
        if g.isdir(j,i): raise Inconsistent('conflict')
        if not g.adj(i,j): raise Inconsistent('no edge')
        if g.isund(i,j): g.orient(i,j)
        g=meek(g)
    if has_cycle(g): raise Inconsistent('cycle')
    if vstructs(g)!=ref: raise Inconsistent('new v')
    return g

def pcp(g,x,y):
    """all proper possibly-causal paths x..y (simple, no v_{i+1}->v_i)"""
    res=[]
    def dfs(path,vis):
        v=path[-1]
        if v==y: res.append(list(path)); return
        for w in sorted(g.pdn(v)):
            if w in vis: continue
            vis.add(w); path.append(w); dfs(path,vis); path.pop(); vis.remove(w)
    dfs([x],{x}); return res

def cn(g,x,y):
    s=set()
    for pth in pcp(g,x,y): s|=set(pth[1:])
    return s

def possde(g,S):
    seen=set(S); st=list(S)
    while st:
        v=st.pop()
        for w in g.pdn(v):
            if w not in seen: seen.add(w); st.append(w)
    return seen

def forb(g,x,y):
    c=cn(g,x,y)
    return ({x} if not c else possde(g,c)|{x})

def amenable(g,x,y):
    P=pcp(g,x,y)
    if not P: return False
    return all(g.isdir(p[0],p[1]) for p in P)

def Ostar(g,x,y):
    if not amenable(g,x,y): return None
    c=cn(g,x,y)
    pa={a for a in range(g.p) for b in c if g.isdir(a,b)}
    return frozenset(pa-forb(g,x,y))

def dist(g,srcs):
    """BFS on skeleton"""
    import collections
    d={v:float('inf') for v in range(g.p)}
    q=collections.deque()
    for s in srcs: d[s]=0; q.append(s)
    while q:
        v=q.popleft()
        for w in range(g.p):
            if w!=v and g.adj(v,w) and d[w]==float('inf'):
                d[w]=d[v]+1; q.append(w)
    return d

# ---- linear SEM population quantities (independent of se.py)
def sigma_of(A,omega):
    p=A.shape[0]; M=np.linalg.inv(np.eye(p)-A.T)
    return M@np.diag(omega)@M.T
def total_effect(A,x,y):
    p=A.shape[0]; return float(np.linalg.inv(np.eye(p)-A)[x,y])
def ols(S,x,y,Z):
    idx=[x]+sorted(Z)
    return float(np.linalg.solve(S[np.ix_(idx,idx)], S[np.ix_(idx,[y])])[0,0])
