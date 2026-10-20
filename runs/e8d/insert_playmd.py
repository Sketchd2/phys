import re
src = open(r'E:\Projects\phys\runs\e8d\playmd_e8d.txt', encoding='utf-8').read()

# two statements made exact against the data
old = """induction tracks the reference (correlation 0.75, slope 0.64: it is
       about 1.5 times too strong per triple) and triples without a hydrogen-bonded
       pair carry nothing (|mean| <= 0.003, sd 0.006); a three-body exchange form
       fitted to them has R^2 0.01-0.05 and the attractive sign."""
new = """induction tracks the reference (correlation 0.75, slope 0.64: it is
       about 1.5 times too strong per triple) and triples whose closest O-O
       is over 3.06 A carry nothing (means -0.002 and +0.003 kcal/mol); a
       three-body exchange form (a product of pair overlaps) fitted to them has
       R^2 0.01-0.05 and the attractive sign."""
assert old in src
src = src.replace(old, new, 1)
old = """0.166-0.170 (E8c 0.304), inner and first-shell slope errors halved, a"""
new = """0.166-0.170 (E8c 0.304), the first-shell slope error cut to a third
       (-0.32 to -0.11) and the inner one from -0.99 to +0.10, a"""
assert old in src
src = src.replace(old, new, 1)

p = r'E:\Projects\phys\docs\PLAY.md'
s = open(p, encoding='utf-8', newline='').read()
nl = '\r\n' if '\r\n' in s else '\n'
marker = "   - **E8z, remove what was not chosen.**"
assert marker in s and 'E8d, where the liquid' not in s
block = src.replace('\r\n', '\n').replace('\n', nl)
s = s.replace(marker, block + marker, 1)
open(p, 'w', encoding='utf-8', newline='').write(s)
print('inserted', len(block), 'chars; newline', repr(nl))
