import numpy as np, sys, subprocess
from scipy_free import nnls
# usage: fitk3.py DIR CLOUD C6 POLAR tags B [cap]  -> K and a Born-Mayer core A_ab exp(-B r) per atom-type pair,
# all non-negative, by weighted least squares on every pair under `cap` kcal/mol (weight 1 / (1 + max(E, 0) / 5)).
H=627.509474; d,cloud,c6,polar,tags,B=sys.argv[1:7]; B=float(B); cap=float(sys.argv[7]) if len(sys.argv)>7 else 100.0
def energies(K):
    law=subprocess.run(['python3','build2.py',cloud,c6,polar,str(K)],capture_output=True,text=True).stdout
    open(d+'/law-tmp.txt','w').write(law); out={}
    for t in tags.split(','):
        r=subprocess.run(['../../../../target/release/phys-lawenergy','law-tmp.txt',t],capture_output=True,text=True,cwd=d).stdout
        for l in r.splitlines():
            w=l.split(); out[(t,int(w[0]))]=(float(w[2])*H,float(w[3])*H)
    return out
e0,e1=energies(0),energies(1)
keys=sorted(e0); types=[(0,0),(0,1),(1,1)]
core={}
for t in tags.split(','):
    for l in open(f'{d}/pairs-water{t}.txt'):
        if l.startswith('#') or '|' not in l: continue
        h,x=l.split('|'); x=x.split(); at=[(np.array(list(map(float,x[i:i+3]))),int(x[i+3])) for i in range(0,len(x),4)]; n=len(at)//2
        v=dict.fromkeys(types,0.0)
        for pa,ta in at[:n]:
            for pb,tb in at[n:]: v[tuple(sorted((ta,tb)))]+=np.exp(-B*np.linalg.norm(pa-pb))
        core[(t,int(h.split()[0]))]=[v[k]*H for k in types]
ref=np.array([e0[k][0] for k in keys]); base=np.array([e0[k][1] for k in keys]); s=np.array([e1[k][1]-e0[k][1] for k in keys]); C=np.array([core[k] for k in keys])
m=ref<cap; w=1/(1+np.maximum(ref[m],0)/5)
X=np.column_stack([s[m],C[m]])*w[:,None]; y=(ref[m]-base[m])*w
c=nnls(X,y); K=c[0]; A=c[1:]
fit=base+K*s+C@A
print('K = %.4f; core exp(-%g r): '%(K,B)+', '.join('%d-%d=%.6g:%g'%(a,b,x,B) for (a,b),x in zip(types,A)))
for lo,hi in ((-9,0),(0,5),(5,20),(20,100)):
    q=(ref>=lo)&(ref<hi)
    if q.any(): print('  ref %4g..%4g n %3d  rms %.3f  mean %+.3f'%(lo,hi,q.sum(),np.sqrt(((fit-ref)[q]**2).mean()),(fit-ref)[q].mean()))
i=ref.argmin(); print('  deepest ref %.2f law %.2f'%(ref[i],fit[i]))
