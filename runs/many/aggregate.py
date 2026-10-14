import glob, re, subprocess, math, sys
K = 627.509474
law = 'runs/e8c/law-bias.txt'
res = {}   # centre -> {n: (pairs, three, cumulative flags)}
for f in sorted(glob.glob('runs/many/mb-hf-c*.txt')):
    c = int(re.search(r'c(\d+)\.txt', f).group(1))
    pair, trip = {}, {}
    for l in open(f):
        w = l.split()
        if not w: continue
        if w[0] == 'pair' and len(w) == 7:
            pair[(int(w[1]), int(w[2]))] = float(w[5])
        elif w[0] == 'triple' and len(w) == 9:
            trip[(int(w[1]), int(w[2]), int(w[3]))] = float(w[7])
    out = {}
    nmax = max([k for (_, k) in pair] + [0]) + 1
    for n in range(3, nmax + 1):
        need_p = [(i, k) for k in range(1, n) for i in range(k)]
        need_t = [(i, j, k) for k in range(2, n) for j in range(1, k) for i in range(j)]
        if all(p in pair for p in need_p) and all(t in trip for t in need_t):
            ps = sum(pair[p] for p in need_p)
            th = sum(trip[(i, j, k)] - pair[(i, j)] - pair[(i, k)] - pair[(j, k)] for (i, j, k) in need_t)
            out[n] = (ps * K, th * K)
    res[c] = out
classical = {}
for c in res:
    o = subprocess.run(['E:/Projects/phys/target/release/phys-manybody.exe', 'water', 'bulk-water-298.snap', law, '--max', '10', '--centres', str(c)], capture_output=True, text=True, cwd='E:/Projects/phys').stdout
    for l in o.splitlines():
        m = re.match(r'\s*(\d+)\s+([+-][\d.]+) \+- [\d.]+\s+[\d.nNaA]+\s+([+-][\d.]+)\s+([+-][\d.]+)\s+([+-][\d.]+)', l)
        if m: classical[(c, int(m.group(1)))] = (float(m.group(2)), float(m.group(3)), float(m.group(5)))
print('centre n | HF pair sum | HF 3-body (sum of triples) | law many-body (induction) | law three-body | law pair-sum')
for c in sorted(res):
    for n in sorted(res[c]):
        cl = classical.get((c, n), (float('nan'),) * 3)
        print(f'{c:4d} {n:2d} | {res[c][n][0]:+9.3f} | {res[c][n][1]:+8.3f} | {cl[0]:+8.3f} | {cl[1]:+8.3f} | {cl[2]:+9.2f}')
print('\nmeans over centres complete at each n:  n  centres  HF3body mean +- se   law many-body mean   law 3-body mean')
for n in range(3, 11):
    q = [res[c][n][1] for c in res if n in res[c]]
    l = [classical[(c, n)][0] for c in res if n in res[c] and (c, n) in classical]
    l3 = [classical[(c, n)][1] for c in res if n in res[c] and (c, n) in classical]
    if not q: continue
    m = sum(q) / len(q)
    se = math.sqrt(sum((x - m) ** 2 for x in q) / (len(q) - 1) / len(q)) if len(q) > 1 else float('nan')
    print(f'{n:3d} {len(q):4d}   {m:+.3f} +- {se:.3f}   {sum(l)/len(l):+.3f}   {sum(l3)/len(l3):+.3f}')
