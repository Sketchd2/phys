import re, subprocess
log = open(r'E:\Projects\phys\runs\e8d\full-test2.log', encoding='utf-8', errors='replace').read()
tracked = {l.strip().replace('tests/', '') for l in subprocess.run(['git', 'ls-files', 'tests'], capture_output=True, text=True, cwd=r'E:\Projects\phys').stdout.splitlines() if l.strip()}
res = {}
cur = None
for m in re.finditer(r'Running (?:tests[\\/](\w+)\.rs|(unittests)[^\n]*)|test result: (\w+)\. (\d+) passed; (\d+) failed; (\d+) ignored', log):
    if m.group(1):
        cur = m.group(1)
    elif m.group(2):
        cur = '(unit tests)'
    else:
        res[cur] = (m.group(3), int(m.group(4)), int(m.group(5)), int(m.group(6)))
tp = tf = ti = 0
ran = 0
for name, (st, p, f, i) in res.items():
    if name == '(unit tests)' or (name + '.rs') in tracked:
        tp += p; tf += f; ti += i; ran += 1
failed = [n for n, (st, p, f, i) in res.items() if f and (n == '(unit tests)' or (n + '.rs') in tracked)]
tracked_rs = sorted(t for t in tracked if t.endswith('.rs'))
missing = [t for t in tracked_rs if t[:-3] not in res]
print('binaries with results (tracked + unit):', ran, '; tracked test files:', len(tracked_rs))
print('passed', tp, 'failed', tf, 'ignored', ti)
print('binaries with failures:', failed)
print('tracked test files with no result in the log:', missing)
print('budget:', res.get('budget'))
