import sys, math
def load(path):
    box=None; mols=[]
    for l in open(path):
        w=l.split()
        if not w: continue
        if w[0]=='box': box=[float(x) for x in w[1:4]]
        elif w[0]=='mol':
            v=[float(x) for x in w[2:]]
            mols.append([v[0:3],v[3:6],v[6:9]])
    return box,mols
def mi(d,b): return [x-b[i]*round(x/b[i]) for i,x in enumerate(d)]
def sub(a,b): return [a[i]-b[i] for i in range(3)]
def norm(a): return math.sqrt(sum(x*x for x in a))
def dot(a,b): return sum(x*y for x,y in zip(a,b))
def analyse(path):
    box,mols=load(path)
    n=len(mols); B=0.529177
    rc=3.3/B
    nb=[[] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i==j: continue
            d=mi(sub(mols[j][0],mols[i][0]),box)
            if norm(d)<rc: nb[i].append(d)
    # tetrahedral order parameter q
    qs=[]
    for i in range(n):
        v=nb[i]
        if len(v)<4: continue
        # take 4 nearest
        v=sorted(v,key=norm)[:4]
        s=0
        for a in range(3):
            for b in range(a+1,4):
                c=dot(v[a],v[b])/(norm(v[a])*norm(v[b]))
                s+=(c+1/3)**2
        qs.append(1-3/8*s)
    # OOO angle hist
    hist=[0]*18
    cnt=0
    for i in range(n):
        v=nb[i]
        for a in range(len(v)):
            for b in range(a+1,len(v)):
                c=dot(v[a],v[b])/(norm(v[a])*norm(v[b]))
                th=math.degrees(math.acos(max(-1,min(1,c))))
                hist[min(17,int(th/10))]+=1; cnt+=1
    # H-bonds: O-O<3.5 A, angle H-O..O <30 deg (donor H on i, acceptor O on j)
    hb=[0]*n
    for i in range(n):
        for j in range(n):
            if i==j: continue
            d=mi(sub(mols[j][0],mols[i][0]),box)
            if norm(d)*B>3.5: continue
            for h in (1,2):
                oh=mi(sub(mols[i][h],mols[i][0]),box)
                c=dot(oh,d)/(norm(oh)*norm(d))
                if math.degrees(math.acos(max(-1,min(1,c))))<30:
                    hb[i]+=1; hb[j]+=1
    hh=[0]*8
    for x in hb: hh[min(7,x)]+=1
    print(path)
    print("  mean q (4 nearest within 3.3 A): %.3f over %d molecules (ideal tetrahedral 1, random 0.0 ice 0.9)"%(sum(qs)/len(qs),len(qs)))
    print("  O-O-O angle histogram (10 deg bins, fraction):",[round(h/cnt,3) for h in hist])
    print("  H-bonds per molecule: mean %.2f, histogram 0..7:"%(sum(hb)/n),hh)
for p in sys.argv[1:]: analyse(p)
