import itertools, networkx as nx
vertices=range(4)
pairs=list(itertools.combinations(vertices,2))
classes={}
counts={}
for code in itertools.product(range(3), repeat=len(pairs)):
    edges=[(a,b) if s==1 else (b,a) for (a,b),s in zip(pairs,code) if s]
    g=nx.DiGraph()
    g.add_nodes_from(vertices)
    g.add_edges_from(edges)
    if not nx.is_directed_acyclic_graph(g): continue
    skeleton=tuple(sorted(tuple(sorted(e)) for e in edges))
    colliders=tuple(sorted((min(a,b),c,max(a,b)) for c in vertices for a,b in itertools.combinations(g.predecessors(c),2) if not g.has_edge(a,b) and not g.has_edge(b,a)))
    key=(skeleton,colliders)
    classes[key]=nx.is_weakly_connected(g)
    counts[key]=counts.get(key,0)+1
print("Labeled four-node Markov equivalence classes:",len(classes))
print("Connected classes:",sum(classes.values()))

print("Classes with at least one undirected edge:",sum(n>1 for n in counts.values()))
