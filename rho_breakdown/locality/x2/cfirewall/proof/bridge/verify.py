from common import *
from itertools import combinations
def mk(p,ds,us):
    G=np.zeros((p,p),dtype=np.int8)
    for i,j in ds: G[i,j]=1
    for i,j in us: G[i,j]=1;G[j,i]=1
    return G
def fmt(G): return " ".join([f"{i}->{j}" for (i,j) in directed_edges(G)]+[f"{i}--{j}" for (i,j) in undirected_edges(G)])
def uc(W):
    p=W.shape[0]; out=set()
    for v in range(p):
        pars=[u for u in range(p) if is_directed(W,u,v)]
        for u,w in combinations(sorted(pars),2):
            if not adjacent(W,u,w): out.add((u,v,w))
    return out

print("### W1: acyclicity of H is necessary (control-arm witness)")
G0=mk(4,[(0,3),(1,2),(3,1),(3,2)],[])
C=mk(4,[],[(0,3),(1,2),(1,3),(2,3)])
D=G0.copy()
print(" C =",fmt(C)," G0 =",fmt(G0)," D =",fmt(D))
print(" D in [G0]:", dag_agrees_with(D,G0) and v_structures(D)==v_structures(C))
H,info=bk_assert(G0,[(1,0)])
print(" H =",fmt(H)," cycle:",has_directed_cycle(H)," extendable:",pdag_extendable(H))
O0,_,a0=X.ostar_and_paths(G0,1,2); O1,_,a1=X.ostar_and_paths(H,1,2)
print(f" amen(G0)={a0} O*(1,2,G0)={sorted(O0)} valid={adjust.is_valid_adjustment_set(D,1,2,set(O0))}")
print(f" amen(H)={a1} O*(1,2,H)={sorted(O1)} valid={adjust.is_valid_adjustment_set(D,1,2,set(O1))}")

print()
print("### W2: G0-amenability is necessary -- the spurious edge MANUFACTURES amenability")
G0=mk(4,[],[(2,3)]); C=G0.copy(); D=mk(4,[(3,2)],[])
print(" C = G0 =",fmt(G0)," D =",fmt(D))
print(" D in [G0]:", dag_agrees_with(D,G0) and v_structures(D)==v_structures(C))
H,info=bk_assert(G0,[(0,2)])
print(" H =",fmt(H)," cycle:",has_directed_cycle(H)," extendable:",pdag_extendable(H),
      " conflict:",info["conflict"])
O0,n0,a0=X.ostar_and_paths(G0,2,3); O1,n1,a1=X.ostar_and_paths(H,2,3)
print(f" amen(G0,2,3)={a0} (npaths={n0})   amen(H,2,3)={a1} O*={sorted(O1)}")
print(" valid in D:",adjust.is_valid_adjustment_set(D,2,3,set(O1)))
Dp=D.copy(); Dp[0,2]=1; P=G0.copy(); P[0,2]=1
print(" C1 (D+ acyclic):",not has_directed_cycle(Dp)," C2 (uc(D+)<=uc(P)):",uc(Dp)<=uc(P))

print()
print("### W3: (C2) is sharp FOR THE BRIDGE -- D+ is not in [H]")
C=mk(4,[],[(1,2),(1,3),(2,3)]); G0=mk(4,[(1,2),(1,3)],[(2,3)]); D=mk(4,[(1,2),(1,3),(3,2)],[])
print(" C =",fmt(C)," G0 =",fmt(G0)," D =",fmt(D))
print(" D in [G0]:", dag_agrees_with(D,G0) and v_structures(D)==v_structures(C))
H,info=bk_assert(G0,[(0,2)]); P=G0.copy(); P[0,2]=1
Dp=D.copy(); Dp[0,2]=1
print(" H =",fmt(H)," D+ =",fmt(Dp))
print(" C1:",not has_directed_cycle(Dp)," C2:",uc(Dp)<=uc(P),
      " uc(D+)=",sorted(uc(Dp))," uc(P)=",sorted(uc(P)))
print(" D+ preserves H? ", all(is_directed(Dp,i,j) for (i,j) in directed_edges(H)))
O0,_,a0=X.ostar_and_paths(G0,1,2); O1,_,a1=X.ostar_and_paths(H,1,2)
print(f" O*(1,2,G0)={sorted(O0)}  O*(1,2,H)={sorted(O1)}  valid in D:",
      adjust.is_valid_adjustment_set(D,1,2,set(O1)))
