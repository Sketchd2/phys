import glob, re, subprocess, math
K=627.509474
hf={}
for f in sorted(glob.glob('runs/many/mb-hf-c*.txt')):
    c=int(re.search(r'c(\d+)\.txt',f).group(1)); pair={}; trip={}
    for l in open(f):
        w=l.split()
        if not w: continue
        if w[0]=='pair' and len(w)==7: pair[(int(w[1]),int(w[2]))]=float(w[5])
        if w[0]=='triple' and len(w)==9: trip[(int(w[1]),int(w[2]),int(w[3]))]=float(w[7])
    for (i,j,k),v in trip.items(): hf[(c,i,j,k)]=(v-pair[(i,j)]-pair[(i,k)]-pair[(j,k)])*K
def law3(path):
    d={}
    for c in (0,20,45,70,95,120,145,170,195,210):
        o=subprocess.run(['E:/Projects/phys/target/release/phys-manybody.exe','water','bulk-water-298.snap',path,'--max','10' if c==0 else '6','--centres',str(c),'--triples'],capture_output=True,text=True,cwd='E:/Projects/phys').stdout
        for l in o.splitlines():
            w=l.split()
            if w and w[0]=='T': d[(int(w[1]),int(w[2]),int(w[3]),int(w[4]))]=float(w[5])
    return d
print('alpha scale | rms(HF3B - law3B) | slope HF on law | mean law3B | mean HF3B   [300 triples, kcal/mol]')
for f in ('1.00','0.95','0.90','0.85','0.80','0.70'):
    import os
    p='runs/e8d/ind/law-a%s.txt'%f
    if not os.path.exists(p):
        out=[]
        for l in open('runs/e8c/avg/refit1/law-opt1.txt'):
            w=l.split()
            out.append('alpha %s %r\n'%(w[1],float(w[2])*float(f)) if w and w[0]=='alpha' else l)
        open(p,'w').writelines(out)
    d=law3(p); keys=[k for k in hf if k in d]
    h=[hf[k] for k in keys]; l=[d[k] for k in keys]
    rms=math.sqrt(sum((a-b)**2 for a,b in zip(h,l))/len(h))
    sxx=sum(b*b for b in l); sxy=sum(a*b for a,b in zip(h,l))
    print('   %s      |      %.4f       |     %.3f       |  %+.4f   | %+.4f'%(f,rms,sxy/sxx,sum(l)/len(l),sum(h)/len(h)))
