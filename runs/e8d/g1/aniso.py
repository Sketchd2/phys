import numpy as np
B=0.529177210903; K=627.5095
le={int(l.split()[0]):l.split() for l in open('le.txt')}
rows=[]
for l in open('pairs-water-liq-gpu-mp2.txt'):
    if l.startswith('#') or '|' not in l: continue
    h,t=l.split('|'); k=int(h.split()[0]); t=t.split()
    a=np.array([[float(t[4*i+j]) for j in range(3)] for i in range(6)])*B
    ref,law=float(le[k][2])*K,float(le[k][3])*K
    best=None
    for d,c in ((0,3),(3,0)):   # donor molecule start d, acceptor start c
        Od,Oa=a[d],a[c]
        for h_ in (1,2):
            H=a[d+h_]; r=np.linalg.norm(H-Oa)
            if best is None or r<best[0]:
                v1=Od-H; v2=Oa-H
                ang=np.degrees(np.arccos(v1@v2/np.linalg.norm(v1)/np.linalg.norm(v2)))
                bis=(a[c+1]-Oa)+(a[c+2]-Oa); bis/=np.linalg.norm(bis)
                u=H-Oa; u/=np.linalg.norm(u)
                th=np.degrees(np.arccos(u@bis))   # 0 = along bisector (away from the lone pairs? H atoms side)
                best=(r,ang,th)
    rows.append((k,float(le[k][1]),ref,law,ref-law)+best)
rows=np.array(rows)
sel=rows[(rows[:,1]>2.7)&(rows[:,1]<4.3)]
print('pairs 2.7-4.3 A:',len(sel))
print('  OO    ref    law    res   rHO  angOHO  theta_bis')
for r in sel[np.argsort(sel[:,6])]: print('%5.2f %6.2f %6.2f %6.2f %5.2f %6.1f %6.1f'%tuple(r[[1,2,3,4,5,6,7]]))
for nm,col in (('rHO',5),('ang O-H..O',6),('theta bisector',7)):
    print(nm,'corr with residual %.2f'%np.corrcoef(sel[:,col],sel[:,4])[0,1])
print()
m=(sel[:,1]>3.2)&(sel[:,1]<3.9)
for lo,hi in ((90,115),(115,135),(135,180)):
    s=sel[m&(sel[:,6]>=lo)&(sel[:,6]<hi)]
    if len(s): print('O-O 3.2-3.9, angle %d-%d: n=%d mean residual %+.3f sd %.3f  mean law %.2f'%(lo,hi,len(s),s[:,4].mean(),s[:,4].std(),s[:,3].mean()))
print('all 47: rms residual %.3f kcal/mol, rms ref %.2f'%(np.sqrt((sel[:,4]**2).mean()),np.sqrt((sel[:,2]**2).mean())))
