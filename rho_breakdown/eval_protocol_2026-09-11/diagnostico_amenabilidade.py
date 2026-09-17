"""Por que as 518 instancias sao apagadas, e se o desbaste e o culpado.

Roda contra ~/bkrobust. Read-only. Nao importa measure.py (ele puxa ortools);
select_knowledge e copiado verbatim de measure.py:103-125.

    python3 diagnostico_amenabilidade.py            # os dois diagnosticos
"""
from __future__ import annotations
import collections, hashlib, itertools, json, os, random, sys
from pathlib import Path

BK = Path(os.path.expanduser("~/bkrobust"))
sys.path.insert(0, str(BK / "src"))

from bkrobust.benchmarks.describe import parse_file
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.example import dag_to_cpdag, knowledge_to_recover
from bkrobust.demo.meek import apply_orientations
from bkrobust.mpdag_criterion.criterion import causal_nodes, forbidden_set, is_amenable
from bkrobust.gac.mpdag_level import is_gac_valid_mpdag

MODELS = BK / "results/axisa3/networks/example_models"
INSTANCES = BK / "results/axisa3/instances.jsonl"


def select_knowledge(dag, cpdag, coverage):
    """Verbatim de benchmarks/measure.py:103-125."""
    k_true = sorted(knowledge_to_recover(dag, cpdag))
    if coverage >= 1.0:
        return k_true
    keep = round(len(k_true) * coverage)
    if keep <= 0:
        return []
    step = len(k_true) / keep
    return [k_true[min(len(k_true) - 1, int(i * step))] for i in range(keep)]


def candidato(g, x, y):
    """pa_g(cn(x,y)) menos o conjunto proibido: a forma fechada, polinomial."""
    cn = causal_nodes(g, x, y)
    pa: set = set()
    for n in cn:
        pa |= g.parents(n)
    return frozenset(pa - forbidden_set(g, x, y) - cn - {x})


def carrega(name):
    p = next((q for q in sorted(MODELS.iterdir()) if q.name.startswith(name)), None)
    if p is None:
        return None
    pn = parse_file(p, hashlib.sha256(p.read_bytes()).hexdigest())
    dag = MPDAG(pn.nodes, directed=pn.edges)
    return dag, dag_to_cpdag(dag)


def diagnostico_518():
    """As rejeicoes o_g0_not_identified sao recuperaveis por Lever 0?"""
    rows = [json.loads(l) for l in open(INSTANCES)]
    rej = [r for r in rows if r.get("reject_reason") == "o_g0_not_identified"]
    cnt = collections.Counter()
    for name in {r["network"] for r in rej}:
        carga = carrega(name)
        if carga is None:
            cnt["sem_arquivo"] += 1
            continue
        dag, cp = carga
        for r in (x for x in rej if x["network"] == name):
            g0 = apply_orientations(cp, select_knowledge(dag, cp, r["coverage"]))
            if not is_amenable(g0, r["x"], r["y"]):
                cnt["nao_amenavel"] += 1
            elif is_gac_valid_mpdag(g0, r["x"], r["y"], candidato(g0, r["x"], r["y"])):
                cnt["recuperavel_por_lever0"] += 1
            else:
                cnt["amenavel_candidato_invalido"] += 1
    return {"total": len(rej), **cnt}


def diagnostico_desbaste(redes=None, cap=120, seed=20260911):
    """A nao-amenabilidade e intrinseca ao tamanho de K, ou artefato de sorted()?"""
    redes = redes or ["mediator", "asia", "Didelez_2010", "sachs", "child",
                      "Sebastiani_2005", "Acid_1996", "Shrier_2008"]
    rows = [json.loads(l) for l in open(INSTANCES)]
    rej = [r for r in rows if r.get("reject_reason") == "o_g0_not_identified"]
    rng = random.Random(seed)
    cnt = collections.Counter()
    for name in redes:
        sel = [r for r in rej if r["network"] == name]
        carga = carrega(name)
        if not sel or carga is None:
            continue
        dag, cp = carga
        k_true = sorted(knowledge_to_recover(dag, cp))
        for r in sel:
            keep = max(0, round(len(k_true) * r["coverage"]))
            if keep == 0 or keep >= len(k_true):
                cnt["degenerado"] += 1          # round(1*0.5) == 0: K fica vazio
                continue
            subs = list(itertools.combinations(k_true, keep))
            if len(subs) > cap:
                subs = rng.sample(subs, cap)
            ok = any(is_amenable(apply_orientations(cp, list(s)), r["x"], r["y"]) for s in subs)
            cnt["existe_subconjunto_amenavel" if ok else "nenhum_subconjunto_amenavel"] += 1
            cnt["_subconjuntos"] += len(subs)
    return dict(cnt)


if __name__ == "__main__":
    a = diagnostico_518()
    b = diagnostico_desbaste()
    print(json.dumps({"rejeicoes_518": a, "artefato_do_desbaste": b}, indent=1, ensure_ascii=False))
    json.dump({"rejeicoes_518": a, "artefato_do_desbaste": b},
              open("diagnostico_amenabilidade.json", "w"), indent=1, ensure_ascii=False)
