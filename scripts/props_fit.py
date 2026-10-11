import math, sys

A = 0.529177210903
rows = []
for l in [x for x in open(sys.argv[1]) if x.count("|") == 2]:
    head, mid, tail = l.split("|")[:3]
    h = head.split()
    r1, r2, th = float(h[1]), float(h[2]), float(h[3])
    m = mid.split()
    t = tail.split()
    rows.append(dict(r1=max(r1, r2), r2=min(r1, r2), th=th, hf=float(m[0]), mu=float(m[1]), alpha=float(t[0]), axx=float(t[1]), ayy=float(t[2]), azz=float(t[3])))
print(len(rows), "geometries")
r0 = 0.97894 / A
th0 = math.radians(104.078)

def solve(a, b):
    n = len(b)
    a = [r[:] for r in a]
    b = b[:]
    for c in range(n):
        p = max(range(c, n), key=lambda i: abs(a[i][c]))
        a[c], a[p] = a[p], a[c]
        b[c], b[p] = b[p], b[c]
        for r in range(c + 1, n):
            f = a[r][c] / a[c][c]
            for k in range(c, n):
                a[r][k] -= f * a[c][k]
            b[r] -= f * b[c]
    x = [0.0] * n
    for r in range(n - 1, -1, -1):
        x[r] = (b[r] - sum(a[r][k] * x[k] for k in range(r + 1, n))) / a[r][r]
    return x

def feats(row, names):
    s = (row["r1"] + row["r2"]) / 2 - r0
    d = row["r1"] - row["r2"]
    u = row["th"] - th0
    table = dict(one=1.0, s=s, u=u, d2=d * d, s2=s * s, u2=u * u, su=s * u)
    return [table[n] for n in names]

def fit(key, names):
    X = [feats(r, names) for r in rows]
    y = [r[key] for r in rows]
    n = len(names)
    ata = [[sum(X[i][a] * X[i][b] for i in range(len(rows))) for b in range(n)] for a in range(n)]
    atb = [sum(X[i][a] * y[i] for i in range(len(rows))) for a in range(n)]
    c = solve(ata, atb)
    res = [y[i] - sum(c[k] * X[i][k] for k in range(n)) for i in range(len(rows))]
    rms = math.sqrt(sum(e * e for e in res) / len(rows))
    return c, rms

sd = lambda key: math.sqrt(sum((r[key] - sum(q[key] for q in rows) / len(rows)) ** 2 for r in rows) / len(rows))
for key, label in (("mu", "MP2 dipole (au)"), ("alpha", "MP2 mean polarisability (bohr^3)")):
    print(f"\n{label}: spread of the data {sd(key):.4f}")
    for names in (["one"], ["one", "s", "u"], ["one", "s", "u", "su", "s2", "u2"], ["one", "s", "u", "su", "s2", "u2", "d2"]):
        c, rms = fit(key, names)
        print(f"  {' '.join(names):<24} rms {rms:.5f}   " + "  ".join(f"{n}={v:+.4f}" for n, v in zip(names, c)))
print("\nsensitivity per 0.01 A of the mean bond and per 5 degrees of angle, and the asymmetric part:")
c, _ = fit("mu", ["one", "s", "u", "su", "s2", "u2", "d2"])
print(f"  dipole: d mu / d s * (0.01 A) = {c[1] * 0.01 / A:+.5f}, d mu / d theta * 5 deg = {c[2] * math.radians(5):+.5f}, d2 coefficient * (0.1 A)^2 = {c[6] * (0.1 / A) ** 2:+.5f}")
c, _ = fit("alpha", ["one", "s", "u", "su", "s2", "u2", "d2"])
print(f"  alpha : d a / d s * (0.01 A) = {c[1] * 0.01 / A:+.5f}, d a / d theta * 5 deg = {c[2] * math.radians(5):+.5f}, d2 coefficient * (0.1 A)^2 = {c[6] * (0.1 / A) ** 2:+.5f}")
