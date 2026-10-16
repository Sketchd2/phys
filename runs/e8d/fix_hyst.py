p='docs/PLAY.md'
s=open(p,encoding='utf-8',newline='').read()
old=s[s.index("Its structure is water-like at 0.95\n"):s.index("     - **A second molecule: methane.**")]
new="""Its structure is the same at fixed density from
       either start: a hysteresis test (0.95 compressed to 1.0, and 1.0 expanded to
       0.95, 30 ps each) gives g(OO) agreeing to about 0.03 at every radius, so
       there is no second structure, and the earlier 4.65 against 7.56 "neighbours
       inside the first minimum" was an artefact of locating a very shallow minimum
       (g about 0.9-1.0, measured about 0.8). At fixed radius the first shell is
       close to measured (4.5 inside 3.3 A at 1.0 against about 4.2); what the law
       lacks is the gap after it (g(3.4) 0.98 against 0.8) and the second shell
       (g(4.4) 0.98 against 1.1). The pressure is history-dependent only through
       sampling noise: -171 +- 75 from the open start, -928 +- 115 from the dense
       one, energies within 0.06-0.09 kcal/mol.
"""
s=s.replace(old,new,1)
open(p,'w',encoding='utf-8',newline='').write(s)
p='runs/many/E8d-notes.md'
s=open(p,encoding='utf-8',newline='').read()
a=s.index("**The structure is density-dependent")
b=s.index("0.95 configuration compressed, and from the 1.0 one expanded).")+len("0.95 configuration compressed, and from the 1.0 one expanded).")
s=s[:a]+"""**The structure is not hysteretic; the real deficiency is the shell gap.** Counting neighbours inside the first
g(OO) minimum made 0.95 and 1.0 g/cm3 look like two structures (4.65 against 7.6), but the minimum is so shallow
that its location jumps. Compared at fixed radius, runs started from the open and from the dense configuration agree
to about 0.03 in g(r) (g(3.4) 0.91/0.93 at 0.95; 0.98/0.97 at 1.0), neighbours inside 3.3 A are 4.47 from either
start, and the energies differ by 0.06-0.09 kcal/mol. The law's first shell is about right; its minimum is too
shallow (g(3.4) ~0.98 against ~0.8) and its second shell too weak (g(4.4) ~0.98 against ~1.1), the signature of the
bent-contact over-attraction at 3.1-3.8 A (section 10.3)."""+s[b:]
open(p,'w',encoding='utf-8',newline='').write(s)
