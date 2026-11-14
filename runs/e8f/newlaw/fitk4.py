import numpy as np, sys, subprocess
# usage: fitk4.py DIR CLOUD C6 POLAR tags [core]  (run from newlaw)
# One constant K, each pair weighted by how much a liquid sees of it: 1 for a bound pair,
# exp(-E / e) for a repulsive one, e the depth of the species' own deepest pair.
H=627.509474; d,cloud,c6,polar,tags=sys.argv[1:6]; core=sys.argv[6:7]
def energies(K):
    law=subprocess.run(['python3','build2.py',cloud,c6,polar,str(K)]+(['core',core[0]] if core else []),capture_output=True,text=True).stdout
    open(d+'/law-tmp.txt','w').write(law); out=[]
    for t in tags.split(','):
        r=subprocess.run(['../../../../target/release/phys-lawenergy','law-tmp.txt',t],capture_output=True,text=True,cwd=d).stdout
        out+=[[float(x) for x in l.split()[:4]] for l in r.splitlines() if l.strip()]
    return np.array(out)
e0,e1=energies(0),energies(1)
sep=e0[:,1]; ref=e0[:,2]*H; b=e0[:,3]*H; s=e1[:,3]*H-b
eps=-ref.min(); w=np.where(ref<0,1.0,np.exp(-np.maximum(ref,0)/eps))
K=(w*(ref-b)*s).sum()/(w*s*s).sum(); Ku=((ref-b)*s)[ref<5].sum()/((s*s)[ref<5]).sum()
print('%s: deepest pair %.3f kcal/mol; K weighted %.4f (unweighted under 5 kcal/mol %.4f)'%(d,-eps,K,Ku))
for name,k in (('weighted',K),('unweighted',Ku)):
    dd=b+k*s-ref; q=np.quantile(sep,[0,0.15,0.3,0.5,0.75,1.0])
    print('  %-10s weighted rms %.3f; bound pairs mean %+.3f rms %.3f; deepest fifth of bound pairs mean %+.3f'%(name,np.sqrt((w*dd*dd).sum()/w.sum()),dd[ref<0].mean(),dd[ref<0].std(),dd[ref<np.quantile(ref[ref<0],0.2)].mean()))
    print('      by separation (A): '+'  '.join('%.1f-%.1f: %+.3f (ref %+.2f)'%(q[i],q[i+1],dd[(sep>=q[i])&(sep<=q[i+1])&(ref<2*eps)].mean(),ref[(sep>=q[i])&(sep<=q[i+1])&(ref<2*eps)].mean()) for i in range(5)))
