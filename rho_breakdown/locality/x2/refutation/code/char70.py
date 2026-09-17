import sys, json, numpy as np, collections
sys.path.insert(0,'/private/tmp/claude-501/-Users-josecosta-mugango/eee5422c-2e7e-4a55-8d52-0664260abb06/scratchpad/ref')
sys.path.insert(0,'/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/code')
import indep as I
from graphs import random_dag
ex=json.load(open('/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/results/exhibit.json'))
rows=[]
for a in ex['all_found']:
    p=a['p']; seed=a['seed']; deg=a['deg']; x=a['x']; y=a['y']
    rng=np.random.default_rng(seed); Dm=random_dag(p,deg,rng)
    dag=[(i,j) for i in range(p) for j in range(p) if Dm[i,j]==1 and Dm[j,i]==0]
    C=I.dag_to_cpdag(dag,p)
    dd=I.dist(C,[x,y])
    ecc=max(dd.values())
    nE=len(dag); maxE=p*(p-1)//2
    nU=len(C.u)
    xy_adj=C.adj(x,y)
    rows.append(dict(ens=a['ens'],arm=a['arm'],p=p,deg=deg,seed=seed,rho=a['rho'],
                     dmins=a['dmins'],ecc=ecc,nE=nE,dens=nE/maxE,nU=nU,xy_adj=bool(xy_adj),
                     frac_at_1=sum(1 for v in dd.values() if v==1)/(p-2 if p>2 else 1)))
print('%-9s %-4s %3s %5s %5s %5s %5s %6s %5s %6s'%('ens','arm','p','rho','ecc','dens','nU','xy_adj','maxdmin','deg'))
for r in rows:
    print('%-9s %-4s %3d %5d %5s %5.2f %5d %6s %5.1f %6.1f'%(r['ens'],r['arm'],r['p'],r['rho'],str(r['ecc']),r['dens'],r['nU'],r['xy_adj'],max(r['dmins']),r['deg']))
print()
print('eccentricity of {X,Y} over the 70:',collections.Counter(r['ecc'] for r in rows))
print('max dmin over flipped, over the 70:',collections.Counter(max(r['dmins']) for r in rows))
print('X adjacent to Y:',collections.Counter(r['xy_adj'] for r in rows))
print('density:',['%.2f'%r['dens'] for r in rows])
print('mean density %.3f  min %.3f  max %.3f'%(np.mean([r['dens'] for r in rows]),min(r['dens'] for r in rows),max(r['dens'] for r in rows)))
print('p distribution:',collections.Counter(r['p'] for r in rows))
print('ecc==1 count:',sum(1 for r in rows if r['ecc']==1),'of',len(rows))
json.dump(rows,open('/tmp/char70.json','w'))
