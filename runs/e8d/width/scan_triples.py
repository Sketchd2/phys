import glob, re, subprocess, math, os
K = 627.509474
hf = {}
for f in sorted(glob.glob(r'E:\Projects\phys\runs\many\mb-hf-c*.txt')):
    c = int(re.search(r'c(\d+)\.txt', f).group(1)); pair = {}; trip = {}
    for l in open(f):
        w = l.split()
        if not w: continue
        if w[0] == 'pair' and len(w) == 7: pair[(int(w[1]), int(w[2]))] = float(w[5])
        if w[0] == 'triple' and len(w) == 9: trip[(int(w[1]), int(w[2]), int(w[3]))] = float(w[7])
    for (i, j, k), v in trip.items():
        hf[(c, i, j, k)] = (v - pair[(i, j)] - pair[(i, k)] - pair[(j, k)]) * K

template = open(r'E:\Projects\phys\runs\e8c\avg\refit1\law-opt1.txt').read().splitlines()
BIN = r'E:\Projects\phys\target\dev9\release\phys-manybody.exe'


def law3(path, width):
    d = {}
    env = dict(os.environ, PHYS_IND_WIDTH=str(width))
    for c in (0, 20, 45, 70, 95, 120, 145, 170, 195, 210):
        o = subprocess.run([BIN, 'water', 'bulk-water-298.snap', path, '--max', '10' if c == 0 else '6', '--centres', str(c), '--triples'],
                           capture_output=True, text=True, cwd=r'E:\Projects\phys', env=env).stdout
        for l in o.splitlines():
            w = l.split()
            if w and w[0] == 'T':
                d[(int(w[1]), int(w[2]), int(w[3]), int(w[4]))] = float(w[5])
    return d


print('width scale | alpha_O alpha_H | rms(HF3B - law3B) | slope HF on law | mean law | mean HF   [300 triples, kcal/mol]')
rows = [(1.0, 9.234446, 0.854324), (1.2, 8.247910, 1.465385), (1.4, 6.344199, 2.472424), (1.6, 3.218981, 3.940734)]
if os.path.exists(r'E:\Projects\phys\runs\e8d\width\extra.txt'):
    for l in open(r'E:\Projects\phys\runs\e8d\width\extra.txt'):
        w = l.split()
        if len(w) == 3: rows.append(tuple(float(x) for x in w))
for f, ao, ah in rows:
    out = []
    for l in template:
        w = l.split()
        if w and w[0] == 'alpha' and w[1] == '0': out.append('alpha 0 %r' % ao)
        elif w and w[0] == 'alpha' and w[1] == '1': out.append('alpha 1 %r' % ah)
        else: out.append(l)
    p = r'E:\Projects\phys\runs\e8d\width\law-w%s.txt' % f
    open(p, 'w').write('\n'.join(out) + '\n')
    d = law3(p, f)
    keys = [k for k in hf if k in d]
    h = [hf[k] for k in keys]; l3 = [d[k] for k in keys]
    rms = math.sqrt(sum((a - b) ** 2 for a, b in zip(h, l3)) / len(h))
    sxx = sum(b * b for b in l3); sxy = sum(a * b for a, b in zip(h, l3))
    corr = sxy / math.sqrt(sxx * sum(a * a for a in h))
    print('   %.1f      | %.3f  %.3f |      %.4f       |     %.3f       | %+.4f | %+.4f   (corr %.3f)' % (f, ao, ah, rms, sxy / sxx, sum(l3) / len(l3), sum(h) / len(h), corr))
