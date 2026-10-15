import glob, re, math
import numpy as np
K=627.509474
hf={}   # (c,i,j,k) -> HF three-body kcal/mol
for f in sorted(glob.glob('runs/many/mb-hf-c*.txt')):
    c=int(re.search(r'c(\d+)\.txt',f).group(1)); pair={}; trip={}
    for l in open(f):
        w=l.split()
        if not w: continue
        if w[0]=='pair' and len(w)==7: pair[(int(w[1]),int(w[2]))]=float(w[5])
        if w[0]=='triple' and len(w)==9: trip[(int(w[1]),int(w[2]),int(w[3]))]=float(w[7])
    for (i,j,k),v in trip.items():
        hf[(c,i,j,k)]=(v-pair[(i,j)]-pair[(i,k)]-pair[(j,k)])*K
rows=[]
for l in open('runs/many/law-triples.txt'):
    w=l.split()
    key=(int(w[1]),int(w[2]),int(w[3]),int(w[4]))
    if key in hf: rows.append((hf[key],float(w[5]),float(w[6]),float(w[7]),float(w[8]),key))
a=np.array([r[:5] for r in rows]); print('triples matched',len(a))
h,l,d1,d2,d3=a.T
r=h-l
bohr=0.529177
d=np.sort(np.stack([d1,d2,d3]),axis=0)*bohr   # angstrom, sorted
print('HF3B   mean %+.4f sd %.4f ; law3B mean %+.4f sd %.4f ; resid mean %+.4f sd %.4f'%(h.mean(),h.std(),l.mean(),l.std(),r.mean(),r.std()))
print('corr(HF3B, law3B) = %.3f   slope %.3f'%(np.corrcoef(h,l)[0,1], np.polyfit(l,h,1)[0]))
# binning by the largest O-O distance and by the smallest
for name,x in (('smallest OO',d[0]),('middle OO',d[1]),('largest OO',d[2]),('sum OO',d.sum(0))):
    qs=np.quantile(x,[0,.25,.5,.75,1.0]); print('--',name)
    for lo,hi in zip(qs[:-1],qs[1:]):
        m=(x>=lo)&(x<=hi); print('  %5.2f-%5.2f A  n=%3d  HF3B %+.3f  law3B %+.3f  resid %+.3f  (resid sd %.3f)'%(lo,hi,m.sum(),h[m].mean(),l[m].mean(),r[m].mean(),r[m].std()))
# overlap-product descriptors  P = prod exp(-(dij)/rho)
print('-- product-of-exponentials descriptor P=prod_pairs exp(-d/rho): fit resid = kappa*P and HF3B = kappa*P')
for rho in (0.3,0.4,0.5,0.7,1.0):
    P=np.exp(-(d.sum(0))/rho)
    P=P/P.max()
    for tgt,name in ((r,'resid'),(h,'HF3B')):
        kap=(P*tgt).sum()/(P*P).sum(); fit=kap*P
        ss=1-((tgt-fit)**2).sum()/((tgt-tgt.mean())**2).sum()
        print('  rho %.1f  %-6s kappa(scaled) %+8.3f  R2 %.3f  corr %+.3f'%(rho,name,kap,ss,np.corrcoef(P,tgt)[0,1]))
