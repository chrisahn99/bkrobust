#!/usr/bin/env python3
"""O exemplo trabalhado da apostila, calculado do zero.

Um unico problema, pequeno o bastante para caber num desenho, sobre o qual TODAS as
definicoes da Parte I sao exibidas com numeros: CPDAG, extensoes, conjunto de ajuste,
Theta, casco, retracao, reversao, portante, tela, orcamento, r_val.

Grava exemplo.json. Nenhum numero da Parte I e digitado a mao.
"""
import json, sys, itertools
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
import numpy as np

from bkrobust.demo.graph import MPDAG
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions
from bkrobust.demo.evaluate import (random_sem, adjusted_estimand,
                                    optimal_adjustment_set_mpdag,
                                    is_valid_adjustment_set_mpdag)
from bkrobust.demo.space import (enumerate_space, represented_dags, covering_pairs,
                                 neighbour_graph, bfs_distances)

# ---------------------------------------------------------------- o mundo verdadeiro
# C confunde X e Y; X age sobre Y atraves de M; W e um vizinho de C fora da consulta.
NODES = ["C", "M", "W", "X", "Y"]
TRUE = [("C", "X"), ("C", "Y"), ("C", "W"), ("X", "M"), ("M", "Y")]
X, Y = "X", "Y"
SEED = 4                      # escolhido para os numeros do casco caberem no texto

dag = MPDAG(nodes=NODES, directed=set(TRUE), undirected=[])
cpdag = dag_to_cpdag(dag)
und = sorted(tuple(sorted(e)) for e in cpdag.undirected_edges)

# as tres afirmacoes do especialista, todas VERDADEIRAS
K = []
for a, b in und:
    K.append((a, b) if (a, b) in dag.directed_edges else (b, a))
g0 = apply_orientations(cpdag, K)
Z = optimal_adjustment_set_mpdag(g0, X, Y)

sem = random_sem(dag, np.random.default_rng(SEED))
tau = sem.true_total_effect(X, Y)


def parents(d, v):
    return {a for a, b in d.directed_edges if b == v}


def theta(g):
    """IDA: o conjunto de efeitos totais possiveis sob g."""
    out = []
    for d in sorted(enumerate_dag_extensions(g), key=lambda d: sorted(d.directed_edges)):
        out.append(dict(edges=sorted("->".join(e) for e in d.directed_edges),
                        pais_de_X=sorted(parents(d, X)),
                        efeito=round(adjusted_estimand(sem, X, Y, parents(d, X)), 4)))
    return out


def hull(rows):
    v = [r["efeito"] for r in rows]
    return [round(min(v), 4), round(max(v), 4)] if v else None


def cobre(h, t, tol=5e-4):
    """O casco arredondado contem o efeito verdadeiro?"""
    return bool(h and h[0] - tol <= t <= h[1] + tol)


def width(h):
    return round(h[1] - h[0], 4) if h else 0.0


def apply(claims):
    """Meek(Chat, claims); None quando o fechamento recusa."""
    return apply_orientations(cpdag, list(claims)) if claims else cpdag


def edges_of(g):
    return dict(dirigidas=sorted("->".join(e) for e in g.directed_edges),
                nao_dirigidas=sorted("--".join(sorted(e)) for e in g.undirected_edges))


out = {"seed": SEED}
out["verdadeiro"] = dict(nos=NODES, arestas=sorted("->".join(e) for e in TRUE),
                         consulta=[X, Y], efeito_verdadeiro=round(tau, 4),
                         pesos={f"{a}->{b}": round(w, 4) for (a, b), w in sorted(sem.weights.items())})
out["cpdag"] = edges_of(cpdag)
out["K"] = ["->".join(e) for e in K]
out["G0"] = edges_of(g0)
out["Z"] = sorted(Z) if Z else []
out["Z_valido_em_G0"] = bool(is_valid_adjustment_set_mpdag(g0, X, Y, Z))

# ---- Theta e cascos --------------------------------------------------------------
out["theta_G0"] = theta(g0)
out["casco_G0"] = hull(out["theta_G0"])
out["theta_blanket"] = theta(cpdag)
out["casco_blanket"] = hull(out["theta_blanket"])
out["largura_blanket"] = width(out["casco_blanket"])

# ---- retracao de cada afirmacao, uma a uma ---------------------------------------
retr = []
for i in range(len(K)):
    Kp = [K[j] for j in range(len(K)) if j != i]
    g = apply(Kp)
    th = theta(g) if g is not None else []
    h = hull(th)
    retr.append(dict(
        afirmacao="->".join(K[i]), resto=["->".join(e) for e in Kp],
        grafo=edges_of(g) if g is not None else None,
        grafo_igual_a_G0=bool(g is not None and g == g0),
        casco=h, largura=width(h), n_extensoes=len(th),
        portante=bool(h is not None and (width(h) > 1e-9 or abs(h[0] - out["casco_G0"][0]) > 1e-9))))
out["retracoes"] = retr
out["portantes"] = [r["afirmacao"] for r in retr if r["portante"]]
out["nao_portantes"] = [r["afirmacao"] for r in retr if not r["portante"]]

# ---- reversao de cada afirmacao --------------------------------------------------
rev = []
for i in range(len(K)):
    Kp = [K[j] if j != i else (K[i][1], K[i][0]) for j in range(len(K))]
    g = apply(Kp)
    th = theta(g) if g is not None else []
    rev.append(dict(afirmacao="->".join(K[i]), invertida="->".join((K[i][1], K[i][0])),
                    meek_recusa=bool(g is None),
                    casco=hull(th), n_extensoes=len(th),
                    efeitos=sorted({r["efeito"] for r in th})))
out["reversoes"] = rev

# ---- L1 sobre este exemplo: retract({k}) = reverse(0) u reverse({k}) -------------
l1 = []
for i in range(len(K)):
    Kp = [K[j] for j in range(len(K)) if j != i]
    g = apply(Kp)
    ext_retr = {tuple(sorted("->".join(e) for e in d.directed_edges))
                for d in enumerate_dag_extensions(g)} if g is not None else set()
    ext_rev = {tuple(sorted("->".join(e) for e in d.directed_edges))
               for d in enumerate_dag_extensions(g0)}
    gi = apply([K[j] if j != i else (K[i][1], K[i][0]) for j in range(len(K))])
    if gi is not None:
        ext_rev |= {tuple(sorted("->".join(e) for e in d.directed_edges))
                    for d in enumerate_dag_extensions(gi)}
    l1.append(dict(afirmacao="->".join(K[i]), n_retracao=len(ext_retr),
                   n_uniao_reversoes=len(ext_rev), iguais=bool(ext_retr == ext_rev)))
out["L1_no_exemplo"] = l1

# ---- a escada RH_j ---------------------------------------------------------------
escada = []
for j in range(len(K) + 1):
    vals, subsets = [], 0
    for r in range(j + 1):
        for S in itertools.combinations(range(len(K)), r):
            Kp = [K[i] for i in range(len(K)) if i not in S]
            g = apply(Kp)
            if g is None:
                continue
            subsets += 1
            vals += theta(g)
    h = hull(vals)
    escada.append(dict(j=j, subconjuntos=subsets, casco=h, largura=width(h),
                       razao_blanket=round(width(h) / out["largura_blanket"], 4)
                       if out["largura_blanket"] > 0 else None,
                       cobre_o_verdadeiro=cobre(h, tau)))
out["escada_RH"] = escada

# ---- L2 sobre este exemplo: casco sobre as portantes = casco sobre todas ---------
def hull_over(indices):
    vals = list(out["theta_G0"])
    for i in indices:
        Kp = [K[j] for j in range(len(K)) if j != i]
        g = apply(Kp)
        if g is not None:
            vals += theta(g)
    return hull(vals)

idx_port = [i for i, r in enumerate(retr) if r["portante"]]
out["L2_no_exemplo"] = dict(
    casco_sobre_todas=hull_over(range(len(K))),
    casco_sobre_portantes=hull_over(idx_port),
    iguais=bool(hull_over(range(len(K))) == hull_over(idx_port)),
    enumeracoes_todas=len(K) + 1, enumeracoes_portantes=len(idx_port) + 1)

# ---- o braco em que o especialista INVERTE a afirmacao portante ------------------
i_port = idx_port[0]
Kerr = [K[j] if j != i_port else (K[i_port][1], K[i_port][0]) for j in range(len(K))]
g0e = apply(Kerr)
Ze = optimal_adjustment_set_mpdag(g0e, X, Y) if g0e is not None else None


def theta_err(g, z):
    """Theta calculado com o conjunto de ajuste que o analista LEU do grafo errado."""
    return theta(g)


escada_err = []
if g0e is not None:
    for j in range(len(Kerr) + 1):
        vals = []
        for r in range(j + 1):
            for S in itertools.combinations(range(len(Kerr)), r):
                Kp = [Kerr[i] for i in range(len(Kerr)) if i not in S]
                g = apply(Kp)
                if g is None:
                    continue
                vals += theta(g)
        h = hull(vals)
        escada_err.append(dict(j=j, casco=h, largura=width(h),
                               razao_blanket=round(width(h) / out["largura_blanket"], 4)
                               if out["largura_blanket"] > 0 else None,
                               cobre_o_verdadeiro=cobre(h, tau)))
out["braco_com_erro"] = dict(
    afirmacao_invertida="->".join((K[i_port][1], K[i_port][0])),
    K=["->".join(e) for e in Kerr],
    meek_aceita=bool(g0e is not None),
    G0=edges_of(g0e) if g0e is not None else None,
    Z=sorted(Ze) if Ze else [],
    Z_valido_no_mundo_verdadeiro=bool(is_valid_adjustment_set_mpdag(dag, X, Y, Ze))
    if Ze is not None else None,
    escada=escada_err)

# ---- r_val sobre este exemplo ----------------------------------------------------
space = enumerate_space(cpdag)
reps = represented_dags(space)
nbrs = neighbour_graph(space, covering_pairs(space, reps))
g0s = next((g for g in space if g == g0), None)
dist = bfs_distances(nbrs, g0s)
camadas = {}
for g in space:
    d = dist.get(g)
    if d is None:
        continue
    ok = bool(is_valid_adjustment_set_mpdag(g, X, Y, Z))
    camadas.setdefault(d, []).append(ok)
bad = sorted(d for d, oks in camadas.items() if d > 0 and not all(oks))
out["r_val"] = dict(
    tamanho_do_espaco=len(space),
    camadas={str(d): dict(n=len(oks), invalidos=int(sum(1 for o in oks if not o)))
             for d, oks in sorted(camadas.items())},
    valor=(bad[0] if bad else None))

json.dump(out, open("exemplo.json", "w"), indent=1, ensure_ascii=False)

print(f"CPDAG nao orientadas : {out['cpdag']['nao_dirigidas']}")
print(f"K                    : {out['K']}")
print(f"Z                    : {out['Z']}   valido em G0: {out['Z_valido_em_G0']}")
print(f"efeito verdadeiro    : {out['verdadeiro']['efeito_verdadeiro']}")
print(f"Theta(G0)            : {sorted({r['efeito'] for r in out['theta_G0']})}")
print(f"casco blanket        : {out['casco_blanket']}  largura {out['largura_blanket']}")
print(f"portantes            : {out['portantes']}")
print(f"nao portantes        : {out['nao_portantes']}")
for r in retr:
    print(f"   retrair {r['afirmacao']:6s} -> casco {r['casco']} "
          f"({r['n_extensoes']} ext, = G0? {r['grafo_igual_a_G0']}, portante {r['portante']})")
print(f"L1 no exemplo        : {[(d['afirmacao'], d['iguais']) for d in l1]}")
print(f"L2 iguais            : {out['L2_no_exemplo']['iguais']}  "
      f"({out['L2_no_exemplo']['enumeracoes_portantes']} vs "
      f"{out['L2_no_exemplo']['enumeracoes_todas']} enumeracoes)")
for e in escada:
    print(f"   RH_{e['j']}: casco {e['casco']} largura {e['largura']} "
          f"razao {e['razao_blanket']} cobre {e['cobre_o_verdadeiro']}")
b = out["braco_com_erro"]
print(f"BRACO COM ERRO: K = {b['K']}  Meek aceita: {b['meek_aceita']}  Z = {b['Z']}")
for e in b["escada"]:
    print(f"   RH_{e['j']}: casco {e['casco']} razao {e['razao_blanket']} cobre {e['cobre_o_verdadeiro']}")
print(f"n extensoes do CPDAG : {len(out['theta_blanket'])}")
print(f"r_val                : {out['r_val']['valor']}   camadas {out['r_val']['camadas']}")
