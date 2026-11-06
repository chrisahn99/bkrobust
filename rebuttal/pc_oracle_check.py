"""Validate the in-house PC: with oracle d-separation tests on the true DAG it
must return the true CPDAG exactly (skeleton, v-structures, Meek closure)."""
import sys
import networkx as nx
import discovery as D
from bkrobust.demo.example import dag_to_cpdag

class Oracle:
    def __init__(self, dag, nodes):
        self.G = nx.DiGraph(); self.G.add_nodes_from(nodes); self.G.add_edges_from(dag.directed_edges); self.nodes = nodes
    def pval(self, i, j, S):
        return 1.0 if nx.is_d_separator(self.G, {self.nodes[i]}, {self.nodes[j]}, {self.nodes[k] for k in S}) else 0.0

for net in sys.argv[1:]:
    dag = D.true_dag(net); nodes, X = D.sample(net, 20, 0)
    deg = max(len(dag.parents(v)) + len(dag.children(v)) for v in dag.nodes)
    D.FisherZ = lambda X_, dag=dag, nodes=nodes: Oracle(dag, nodes)
    g, info = D.pc(nodes, X, alpha=0.5, max_depth=deg)
    cp = dag_to_cpdag(dag)
    print(net, "max degree", deg, "equal CPDAG:", g == cp, info, "skel err", D.skeleton_errors(g, dag),
          "dir diff", len(set(g.directed_edges) ^ set(cp.directed_edges)),
          "und diff", len(set(g.undirected_edges) ^ set(cp.undirected_edges)), flush=True)
