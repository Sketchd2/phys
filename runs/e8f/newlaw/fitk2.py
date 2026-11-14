import numpy as np, sys, subprocess
# usage: fitk2.py DIR CLOUD C6 POLAR tag[,tag..] [core]   -> prints K and the errors; run from newlaw
H=627.509474; d,cloud,c6,polar,tags=sys.argv[1:6]; core=sys.argv[6:7]
def energies(K):
    law=subprocess.run(['python3','build2.py',cloud,c6,polar,str(K)]+(['core',core[0]] if core else []),capture_output=True,text=True).stdout
    open(d+'/law-tmp.txt','w').write(law); out={}
    for t in tags.split(','):
        r=subprocess.run(['../../../../target/release/phys-lawenergy','law-tmp.txt',t],capture_output=True,text=True,cwd=d).stdout
        out[t]=np.array([[float(x) for x in l.split()[:4]] for l in r.splitlines() if l.strip()])
    return out
e0,e1=energies(0),energies(1)
def get(ts,cap=5.0):
    ref=np.concatenate([e0[t][:,2] for t in ts])*H; b=np.concatenate([e0[t][:,3] for t in ts])*H; s=np.concatenate([e1[t][:,3] for t in ts])*H-b
    m=ref<cap; return ref[m],b[m],s[m]
ts=tags.split(','); ref,b,s=get(ts); K=((ref-b)*s).sum()/(s*s).sum()
print('K = %.4f from %d pairs under 5 kcal/mol'%(K,len(ref)))
for t in ts:
    ref,b,s=get([t]); dd=b+K*s-ref; m=ref<0
    print('  %-10s n %3d rms %.3f mean %+.3f | bound n %3d mean %+.3f rms %.3f | deepest ref %.2f law %.2f'%(t,len(ref),np.sqrt((dd*dd).mean()),dd.mean(),m.sum(),dd[m].mean() if m.any() else 0,dd[m].std() if m.any() else 0,ref.min(),(b+K*s)[ref.argmin()]))
