"""Merge the E8c slope data with the G1 liquid-sampled slopes, holding out the odd G1 pairs."""
import os, shutil
base = r'E:\Projects\phys\runs\e8c\avg\refit1'
g1 = r'E:\Projects\phys\runs\e8d\g1'
out = r'E:\Projects\phys\runs\e8d\fit2'
os.makedirs(out, exist_ok=True)
for f in ['pairs-water-gpu-bias-mp2.txt', 'pairs-water-gpu-mp2.txt', 'pairs-water-gpu-wall-mp2.txt', 'rotate-water.txt', 'polar-water.txt',
          'esp-water-cloud-mp2-fit.txt', 'grow-water.state']:
    shutil.copy(os.path.join(base, f), out)
shutil.copy(os.path.join(g1, 'pairs-water-liq-gpu-mp2.txt'), out)

# merged pair-geometry file for the slopes: the bias pairs as they are, the G1 pairs at index + 1000
lines = [l for l in open(os.path.join(base, 'pairs-water-gpu-bias-mp2.txt')) if l.strip()]
g1_lines = []
for l in open(os.path.join(g1, 'pairs-water-liq-gpu-mp2.txt')):
    if l.startswith('#') or '|' not in l:
        continue
    head, tail = l.split('|', 1)
    w = head.split()
    w[0] = str(int(w[0]) + 1000)
    g1_lines.append(' '.join(w) + ' |' + tail)
open(os.path.join(out, 'pairs-water-slp.txt'), 'w').writelines(lines + g1_lines)

# slopes: the E8c ones, plus the G1 even-indexed ones for training; the odd ones kept for a test
train = [l for l in open(os.path.join(base, 'deriv-water.txt')) if l.strip()]
test = []
for l in open(os.path.join(g1, 'deriv-water.txt')):
    w = l.split()
    if len(w) < 4:
        continue
    k = int(w[0])
    w[0] = str(k + 1000)
    if k % 2 == 0:
        train.append(' '.join(w) + '\n')
    else:
        test.append(' '.join(w) + '\n')
open(os.path.join(out, 'deriv-train.txt'), 'w').writelines(train)
open(os.path.join(out, 'deriv-test.txt'), 'w').writelines(test)
open(os.path.join(out, 'deriv-water.txt'), 'w').writelines(train)
print('E8c slopes', len(train) - sum(1 for l in train if int(l.split()[0]) >= 1000), '+ G1 train', sum(1 for l in train if int(l.split()[0]) >= 1000), '; G1 test', len(test))
