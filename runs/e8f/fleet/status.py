#!/usr/bin/env python3
"""Write fleet/status.html: where every species and every worker of the fleet stands.
usage, from anywhere in the repository:   python3 runs/e8f/fleet/status.py
Reads fleet/species.tsv, the *.ok files and liquid states under species/, the queues, and the
workers' branches as last fetched (run `git fetch origin 'refs/heads/claude/e8f-w-*:refs/remotes/origin/claude/e8f-w-*'`
first). fleet/state.json, if the orchestrator keeps one, adds its own notes:
{"children": 20, "window": "allowed, resets 14:30 UTC", "note": "..."}.
The page is then published with the Artifact tool from the same path, so its link stays the same."""
import os, json, subprocess, html, time, glob, re
F = os.path.dirname(os.path.abspath(__file__)); E = os.path.dirname(F); ROOT = os.path.dirname(os.path.dirname(E))
STAGES = ['grow', 'mono', 'pairs', 'full', 'law', 'plan']
def sh(*a):
    return subprocess.run(a, capture_output=True, text=True, cwd=ROOT).stdout
species = [l.rstrip('\n').split('\t') for l in open(F + '/species.tsv') if not l.startswith('#')]
state = json.load(open(F + '/state.json')) if os.path.exists(F + '/state.json') else {}
# fleet/held.tsv: species<TAB>stage<TAB>reason, for a stage the orchestrator chose not to carry on with (not an error)
held = {}
if os.path.exists(F + '/held.tsv'):
    for l in open(F + '/held.tsv'):
        if l.strip() and not l.startswith('#'):
            n, st, why = (l.rstrip('\n').split('\t') + ['', ''])[:3]; held[(n, st)] = why
rows = []; totals = {s: 0 for s in STAGES + ['liquid', 'boil']}; nheld = 0
for name, heavy, spec in species:
    d = f'{E}/species/{name}'; cells = []
    for s in STAGES:
        ok = os.path.exists(f'{d}/{s}.ok'); hole = s == 'law' and ok and open(f'{d}/law.ok').read().startswith('HOLE')
        started = os.path.isdir(d) and any(glob.glob(f'{d}/{p}') for p in {'grow': ['grow.log', 'grow.status'], 'mono': ['density.log'], 'pairs': ['geom.log'], 'full': ['full.log'], 'law': ['law.log'], 'plan': ['liquid']}[s])
        stopped = not ok and (name, s) in held
        cells.append('hole' if hole else 'ok' if ok else 'held' if stopped else 'run' if started else 'none'); totals[s] += ok and not hole
    want = sum(1 for _ in open(f'{d}/liquid/states.txt')) if os.path.exists(f'{d}/liquid/states.txt') else 0
    done = len(glob.glob(f'{d}/liquid/s-*/done')); totals['liquid'] += want > 0 and done >= want
    boil = open(f'{d}/boil.ok').read().strip() if os.path.exists(f'{d}/boil.ok') else ''; totals['boil'] += bool(boil)
    why = '; '.join(dict.fromkeys(held[(name, st)] for st in STAGES if (name, st) in held and cells[STAGES.index(st)] == 'held')); nheld += bool(why)
    if 'hole' in cells:  # a law with no liquid: say why, from law.ok
        t = open(f'{d}/law.ok').read(); bits = []
        b2 = re.search(r'300 K: ([^;\n]+)', t)
        if b2: bits.append(f'No liquid: B2(300 K) = {b2.group(1).strip()} cm3/mol, so the law has a short-range hole')
        rms = re.search(r'weighted rms ([0-9.]+)', t)
        if rms: bits.append(f'Fit rms {rms.group(1)} kcal/mol')
        dp = re.search(r'deepest: reference ([-0-9.e+]+), law ([-0-9.e+]+)', t)
        if dp: bits.append(f'deepest pair {dp.group(1)} by reference, {float(dp.group(2)):+.3f} by the law')
        why = '; '.join(x for x in ['. '.join([bits[0], ', '.join(bits[1:])]) if len(bits) > 1 else ''.join(bits), why] if x)
    rows.append((name, heavy, cells, done, want, boil, why))
workers = []
for q in sorted(glob.glob(F + '/queue/*.txt'), key=lambda p: int(os.path.basename(p)[:-4])):
    w = os.path.basename(q)[:-4]; jobs = [l.split()[0] for l in open(q) if l.strip() and not l.startswith('STOP')]
    ref = f'origin/claude/e8f-w-{w}'
    last = sh('git', 'log', '-1', '--format=%ct|%s', ref).strip()
    if not last: workers.append((w, len(jobs), 0, 0, None, 'never pushed')); continue
    ct, subj = last.split('|', 1)
    ls = lambda p: [x.split('/')[-1] for x in sh('git', 'ls-tree', '--name-only', ref, f'runs/e8f/fleet/{p}/').split()]
    done, failed = [j for j in ls('done') if j in jobs], [j for j in ls('failed') if j in jobs]
    workers.append((w, len(jobs), len(done), len(failed), int(ct), subj.replace(f'fleet worker {w}: ', '')))
now = int(time.time())
def age(t):
    if t is None: return '—'
    m = (now - t) // 60; return f'{m} min' if m < 120 else f'{m // 60} h'
live = sum(1 for w in workers if w[4] and now - w[4] < 2700 and w[2] + w[3] < w[1])
e = html.escape
out = ['''<title>Boiling Point Fleet</title>
<meta name="claude:viewport-fit-sides" content="cover">
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600&family=Space+Grotesk:wght@500;600&family=JetBrains+Mono:wght@400;500&display=swap">
<style>
/* layout: a status line, then the species grid (one row a species, one cell a stage), then the workers */
:root { --bg:#f6f7f9; --panel:#ffffff; --fg:#161b22; --muted:#5d6876; --line:#dfe3e8; --ok:#1f8a5b; --run:#b7791f; --bad:#c2413b; --none:#c9cfd6; --held:#7a5cc7; --accent:#2f5fb3;
  --display:"Space Grotesk",system-ui,sans-serif; --body:"Inter",system-ui,sans-serif; --mono:"JetBrains Mono",ui-monospace,monospace }
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) { --bg:#0f1318; --panel:#171c23; --fg:#e6eaef; --muted:#93a0af; --line:#2a323c; --ok:#4cc38a; --run:#e0a94a; --bad:#ef6f69; --none:#3a434e; --held:#a98cf0; --accent:#7ba7ee; color-scheme:dark } }
:root[data-theme="dark"] { --bg:#0f1318; --panel:#171c23; --fg:#e6eaef; --muted:#93a0af; --line:#2a323c; --ok:#4cc38a; --run:#e0a94a; --bad:#ef6f69; --none:#3a434e; --held:#a98cf0; --accent:#7ba7ee; color-scheme:dark }
body { background:var(--bg); color:var(--fg); font:14px/1.5 var(--body); padding-block:24px; padding-left:max(16px,env(safe-area-inset-left,0px)); padding-right:max(16px,env(safe-area-inset-right,0px)) }
main { max-width:1040px; margin:0 auto; display:flex; flex-direction:column; gap:24px }
h1 { font:600 26px/1.2 var(--display); margin:0; text-wrap:balance } h2 { font:600 16px/1.3 var(--display); margin:0 0 10px }
.sub { color:var(--muted); margin:4px 0 0 } .mono, td.n { font-family:var(--mono); font-variant-numeric:tabular-nums }
.tiles { display:grid; grid-template-columns:repeat(auto-fit,minmax(130px,1fr)); gap:10px }
.tile { background:var(--panel); border:1px solid var(--line); border-radius:8px; padding:10px 12px } .tile b { display:block; font:600 22px/1.2 var(--display); font-variant-numeric:tabular-nums } .tile span { color:var(--muted); font-size:12px }
.wrap { overflow-x:auto; background:var(--panel); border:1px solid var(--line); border-radius:8px }
table { border-collapse:collapse; width:100%; min-width:560px } th,td { padding:6px 10px; text-align:left; border-bottom:1px solid var(--line); white-space:nowrap }
th { font:500 11px/1.4 var(--body); letter-spacing:.05em; text-transform:uppercase; color:var(--muted) } tr:last-child td { border-bottom:0 }
.dot { display:inline-block; width:14px; height:14px; border-radius:3px; background:var(--none); vertical-align:-2px } .dot.ok { background:var(--ok) } .dot.run { background:var(--run) } .dot.hole { background:var(--bad) } .dot.held { background:var(--held) }
.pill { font:500 12px var(--body); padding:1px 8px; border-radius:99px; border:1px solid var(--line) } .pill.live { color:var(--ok); border-color:var(--ok) } .pill.quiet { color:var(--bad); border-color:var(--bad) } .pill.done { color:var(--muted) }
td.why { white-space:normal; min-width:260px; color:var(--muted); font-size:13px }
.key { display:flex; flex-wrap:wrap; gap:14px; color:var(--muted); font-size:12px; margin-top:8px } .note { color:var(--muted) }
</style>
<main>''']
out.append(f'<header><h1>Boiling Point Fleet</h1><p class="sub">Updated {time.strftime("%Y-%m-%d %H:%M UTC", time.gmtime(now))}. {e(str(state.get("note", "")))}</p></header>')
tiles = [(len(species), 'species'), (totals['grow'], 'bases grown'), (totals['pairs'], 'pair sets done'), (totals['law'], 'laws built'), (totals['boil'], 'boiling points'), (nheld, 'species stopped'), (live, 'workers pushing'), (state.get('children', '—'), 'children allowed'), (state.get('window', '—'), 'five-hour window')]
out.append('<section class="tiles">' + ''.join(f'<div class="tile"><b>{e(str(v))}</b><span>{e(l)}</span></div>' for v, l in tiles) + '</section>')
out.append('<section><h2>Species</h2><div class="wrap"><table><tr><th>Species</th><th>Heavy atoms</th>' + ''.join(f'<th>{s}</th>' for s in STAGES) + '<th>Liquid states</th><th>Boiling point</th><th>Note</th></tr>')
for name, heavy, cells, done, want, boil, why in rows:
    out.append(f'<tr><td>{e(name.replace("_", " "))}</td><td class="n">{heavy}</td>' + ''.join(f'<td><span class="dot {c}" title="{c}"></span></td>' for c in cells) + f'<td class="n">{f"{done} / {want}" if want else "—"}</td><td class="n">{e(boil) or "—"}</td><td class="why">{e(why)}</td></tr>')
out.append('</table></div><div class="key"><span><span class="dot ok"></span> finished</span><span><span class="dot run"></span> started</span><span><span class="dot hole"></span> law has a hole</span><span><span class="dot held"></span> stopped by the orchestrator (reason in Note)</span><span><span class="dot"></span> not started</span></div></section>')
out.append('<section><h2>Workers</h2><div class="wrap"><table><tr><th>Worker</th><th>State</th><th>Jobs done</th><th>Failed</th><th>Queued</th><th>Last push</th><th>Last report</th></tr>')
for w, n, d, f, t, subj in workers:
    st = ('done', 'finished') if n and d + f >= n else ('live', 'pushing') if t and now - t < 2700 else ('quiet', 'quiet')
    out.append(f'<tr><td class="n">{w}</td><td><span class="pill {st[0]}">{st[1]}</span></td><td class="n">{d}</td><td class="n">{f}</td><td class="n">{n}</td><td class="n">{age(t)}</td><td>{e(subj[:60])}</td></tr>')
out.append('</table></div><p class="note">A worker is quiet when it has unfinished jobs and has not pushed for 45 minutes; the orchestrator restarts or replaces it.</p></section></main>')
open(F + '/status.html', 'w').write('\n'.join(out))
print(f'{F}/status.html: {len(species)} species, {len(workers)} workers, {live} pushing')
