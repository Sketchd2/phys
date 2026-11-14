import numpy as np, sys
H=627.509474
def ld(f): return np.array([[float(x) for x in l.split()[:4]] for l in open(f) if l.strip()])
sets=['-gpu-bias-mp2','-gpu-mp2','-liq-gpu-mp2','-gpu-wall-mp2']
D={t:(ld('e0%s.txt'%t),ld('e1%s.txt'%t)) for t in sets}
def get(ts,cap=5.0):
    ref=np.concatenate([D[t][0][:,2] for t in ts])*H; e0=np.concatenate([D[t][0][:,3] for t in ts])*H; s=np.concatenate([D[t][1][:,3] for t in ts])*H-e0
    m=ref<cap; return ref[m],e0[m],s[m]
ref,e0,s=get(sets[:3]); K=((ref-e0)*s).sum()/(s*s).sum()
print('K fitted to all non-wall pairs under 5 kcal/mol: %.3f (%d pairs)'%(K,len(ref)))
for ts in ([sets[0]],[sets[1]],[sets[2]],[sets[3]]):
    ref,e0,s=get(ts); d=e0+K*s-ref; b=ref<0
    print('  %-14s n %3d rms %.3f mean %+.3f | bound n %3d mean %+.3f rms %.3f'%(ts[0],len(ref),np.sqrt((d*d).mean()),d.mean(),b.sum(),d[b].mean(),d[b].std()))
# how few pairs are enough: K from random subsets of 10 and 20
rng=np.random.default_rng(1); ref,e0,s=get(sets[:3])
for n in (5,10,20,40):
    ks=[]
    for _ in range(400):
        i=rng.choice(len(ref),n,replace=False); ks.append(((ref[i]-e0[i])*s[i]).sum()/(s[i]*s[i]).sum())
    print('  K from %2d pairs: %.3f +- %.3f'%(n,np.mean(ks),np.std(ks)))
