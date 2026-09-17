import json,sys
from collections import defaultdict
f=sys.argv[1]
d=json.load(open(f))
C=d["counts"]
cells=defaultdict(int)
for k,v in C.items():
    p=k.split("|")
    if p[0]!="cell": continue
    _,b,truth,s0,s1,snd = p
    cells[(b,truth,s0,s1,snd)]+=v
def rate(k,n): return f"{k}/{n} = {100*k/n:.2f}%" if n else "n=0"
print("=== ",f)
for truth in sorted(set(t for (_,t,_,_,_) in cells)):
  for b in ["0","1","2","3","disc"]:
    sub={(s0,s1,snd):v for (bb,t,s0,s1,snd),v in cells.items() if bb==b and t==truth}
    if not sub: continue
    Dbot=sum(v for (s0,_,_),v in sub.items() if s0=="REFUSE")
    fg_adj=sum(v for (s0,s1,snd),v in sub.items() if s0=="REFUSE" and s1=="POS" and snd=="False")
    fg_zero=sum(v for (s0,s1,snd),v in sub.items() if s0=="REFUSE" and s1=="ZERO" and snd=="False")
    g_adj=sum(v for (s0,s1,snd),v in sub.items() if s0=="REFUSE" and s1=="POS")
    g_any=sum(v for (s0,s1,snd),v in sub.items() if s0=="REFUSE" and s1!="REFUSE")
    Did=sum(v for (s0,_,_),v in sub.items() if s0!="REFUSE")
    sc=sum(v for (s0,s1,snd),v in sub.items() if s0!="REFUSE" and snd=="False")
    print(f" truth={truth} dmin={b}: D_bot={Dbot} FGR={rate(fg_adj+fg_zero,Dbot)} "
          f"[FG-adjust {rate(fg_adj,Dbot)} | FG-zero {rate(fg_zero,Dbot)}]  "
          f"P(wrong|adjust-gain)={rate(fg_adj,g_adj)}  P(wrong|any gain)={rate(fg_adj+fg_zero,g_any)}  "
          f"|| D_id={Did} SCR={rate(sc,Did)}")
tot_id=sum(v for (b,t,s0,s1,snd),v in cells.items() if s0!="REFUSE")
tot_id_bad=sum(v for (b,t,s0,s1,snd),v in cells.items() if s0!="REFUSE" and snd=="False")
print(f" TOTAL identified-base trials {tot_id}, unsound {tot_id_bad}")
