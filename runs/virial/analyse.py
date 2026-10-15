import re, math, subprocess
K=627.509474; BAR=69477.0; RHO=0.03343; n=5
def load(path):
    pair={};trip={}
    for l in open(path):
        w=l.split()
        if not w: continue
        if w[0]=='pair' and len(w)==7: pair[(int(w[1]),int(w[2]))]=float(w[5])
        if w[0]=='triple' and len(w)==9: trip[(int(w[1]),int(w[2]),int(w[3]))]=float(w[7])
    return pair,trip
def e3(pair,trip):
    out={}
    for (i,j,k),v in trip.items():
        if k<n and all(p in pair for p in ((i,j),(i,k),(j,k))):
            out[(i,j,k)]=(v-pair[(i,j)]-pair[(i,k)]-pair[(j,k)])*K
    return out
p1,t1=load('runs/many/mb-hf-c0.txt'); p0,t0=load('runs/virial/vir-c0-n5-L0.99.txt')
h1=e3(p1,t1); h0=e3(p0,t0)
common=sorted(set(h1)&set(h0))
print('triples with both scalings:',len(common))
dl=-math.log(0.99)
def law(scale):
    o=subprocess.run(['E:/Projects/phys/target/release/phys-manybody.exe','water','bulk-water-298.snap','runs/e8c/avg/refit1/law-opt1.txt','--max',str(n),'--centres','0','--triples','--scale',str(scale)],capture_output=True,text=True,cwd='E:/Projects/phys').stdout
    d={}
    for l in o.splitlines():
        w=l.split()
        if w and w[0]=='T': d[(int(w[2]),int(w[3]),int(w[4]))]=float(w[5])
    return d
l1=law(1.0); l0=law(0.99)
sh=sum(h1[t]-h0[t] for t in common)/dl; sl=sum(l1[t]-l0[t] for t in common)/dl
print('per triple dE3/dln(lambda) [kcal/mol]: HF vs law induction (+ = energy rises on expansion)')
for t in common: print('  ',t,'E3(1.00) HF %+.4f law %+.4f | E3(0.99) HF %+.4f law %+.4f | dE/dlnL HF %+.3f law %+.3f'%(h1[t],l1[t],h0[t],l0[t],(h1[t]-h0[t])/dl,(l1[t]-l0[t])/dl))
print('sum over the triples: dE3/dlnL  HF %+.3f   law %+.3f kcal/mol'%(sh,sl))
for name,s in (('HF',sh),('law induction',sl)):
    P=-(RHO/(3*n))*s*BAR
    print('P3 (one cluster of %d per molecule, rho=%.4f/A^3): %s %+.0f bar'%(n,RHO,name,P))
