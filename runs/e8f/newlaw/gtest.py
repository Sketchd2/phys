#!/usr/bin/env python3
"""How well a Gaussian model of a molecule's density (gaussfit.py, the density without the 1s
cores, nuclei as points of charge Z - core) gives the exact first-order electrostatic energy
and the density overlap of pairs.  usage: gtest.py NAME PAIRS OVERLAP  (run from runs/e8f after gaussfit)"""
import sys, numpy as np
from math import erf, sqrt, pi
H=627.509474; name,pf,of=sys.argv[1:4]
atoms=[(int(l.split()[1]),np.array(list(map(float,l.split()[2:5])))) for l in open(f'density-{name}.txt') if l.startswith('atom')]
ref=np.array([r for _,r in atoms]); zeff=[z-2 if z>2 else z for z,_ in atoms]
model=np.load(f'gauss-{name}.npy')
def place(at):
    at=np.array(at); ca,cr=at.mean(0),ref.mean(0)
    U,_,Vt=np.linalg.svd((ref-cr).T@(at-ca)); d=np.sign(np.linalg.det(U@Vt)); R=U@np.diag([1,1,d])@Vt
    return [((m[:3]-cr)@R+ca,m[3],-m[4]) for m in model]+[(p,0.0,z) for p,z in zip(at,zeff)]
geo={}
for l in open(pf):
    if l.startswith('#') or '|' not in l: continue
    h,t=l.split('|'); t=t.split(); geo[int(h.split()[0])]=[np.array(list(map(float,t[q:q+3]))) for q in range(0,len(t),4)]
out=[]; n=len(atoms)
for l in open(of):
    v=l.split()
    if len(v)!=9: continue
    try: k=int(v[0])
    except: continue
    g=geo[k]; A,B=place(g[:n]),place(g[n:]); E=S=0.0
    for pa,sa,qa in A:
        for pb,sb,qb in B:
            r=np.linalg.norm(pa-pb); s2=sa*sa+sb*sb
            E+=qa*qb*(erf(r/sqrt(s2))/r if s2>0 else 1/r)
            if sa>0 and sb>0: S+=qa*qb*(pi*s2)**-1.5*np.exp(-r*r/s2)
    out.append((float(v[2])*H,float(v[4])*H,E*H,float(v[6]),S,float(v[8])*H))
o=np.array(out); m=o[:,0]<5
print(f"{name}: {len(model)} Gaussians, {m.sum()} pairs under 5 kcal/mol: electrostatics model - exact rms {np.sqrt(((o[m,2]-o[m,1])**2).mean()):.3f} mean {(o[m,2]-o[m,1]).mean():+.3f} kcal/mol (exact rms {np.sqrt((o[m,1]**2).mean()):.2f}); overlap model/exact {np.mean(o[m,4]/o[m,3]):.3f} spread {np.std(o[m,4]/o[m,3])/np.mean(o[m,4]/o[m,3])*100:.1f}%")
ex=o[m,5]-o[m,1]; K=(ex*o[m,4]*H).sum()/((o[m,4]*H)**2).sum(); print(f"   exchange = K x model overlap: K {K:.3f}, rms {np.sqrt(((K*o[m,4]*H-ex)**2).mean()):.3f} kcal/mol")
