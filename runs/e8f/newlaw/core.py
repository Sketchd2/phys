import numpy as np, sys, subprocess
from scipy_free import nnls
# usage: DAMP=d core.py DIR CLOUD C6 POLAR K tags B : a Born-Mayer core A_ab exp(-B r) per atom-type pair, A >= 0,
# fitted to what the law leaves over on the pairs above 5 kcal/mol (relative error), the law itself held.
H=627.509474; d,cloud,c6,polar,K,tags,B=sys.argv[1:8]; B=float(B)
law=subprocess.run(['python3','build2.py',cloud,c6,polar,K],capture_output=True,text=True).stdout
open(d+'/law-tmp.txt','w').write(law); types=[(0,0),(0,1),(1,1)]; rows=[]
for t in tags.split(','):
    r=subprocess.run(['../../../../target/release/phys-lawenergy','law-tmp.txt',t],capture_output=True,text=True,cwd=d).stdout
    e={int(l.split()[0]):(float(l.split()[2])*H,float(l.split()[3])*H) for l in r.splitlines()}
    for l in open(f'{d}/pairs-water{t}.txt'):
        if l.startswith('#') or '|' not in l: continue
        h,x=l.split('|'); x=x.split(); at=[(np.array(list(map(float,x[i:i+3]))),int(x[i+3])) for i in range(0,len(x),4)]; n=len(at)//2
        v=dict.fromkeys(types,0.0)
        for pa,ta in at[:n]:
            for pb,tb in at[n:]: v[tuple(sorted((ta,tb)))]+=np.exp(-B*np.linalg.norm(pa-pb))*H
        rows.append(list(e[int(h.split()[0])])+[v[k] for k in types])
a=np.array(rows); ref,lw,C=a[:,0],a[:,1],a[:,2:]
m=ref>5; w=1/ref[m]; A=nnls(C[m]*w[:,None],(ref[m]-lw[m])*w,200000); fit=lw+C@A
print('core '+','.join('%d-%d=%.6g:%g'%(x,y,v,B) for (x,y),v in zip(types,A)))
for lo,hi in ((-9,0),(0,5),(5,20),(20,100),(100,1e9)):
    q=(ref>=lo)&(ref<hi)
    if q.any(): print('  ref %4g..%4g n %3d  before: mean %+.2f rms %.2f   with core: mean %+.2f rms %.2f'%(lo,hi,q.sum(),(lw-ref)[q].mean(),(lw-ref)[q].std(),(fit-ref)[q].mean(),(fit-ref)[q].std()))
