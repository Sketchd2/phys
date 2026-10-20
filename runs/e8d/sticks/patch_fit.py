p = r'E:\Projects\phys\runs\e8d\sticks\fit_sticks.py'
s = open(p, encoding='utf-8').read()

# allow C6 held at the engine's derived values (E6: O-O 17.48, O-H 6.331, H-H 2.376 hartree bohr^6)
s = s.replace("def fit(variant, qh_fixed=None):", "def fit(variant, qh_fixed=None, c6_fixed=None):")
s = s.replace("""    def unpack(p):
        qh = qh_fixed if qh_fixed is not None else p[9]
        return p[:9], qh
""", """    def unpack(p):
        qh = qh_fixed if qh_fixed is not None else p[-1]
        pp = np.concatenate([p[:6], np.array(c6_fixed)]) if c6_fixed is not None else p[:9]
        return pp, qh
""")
s = s.replace("""    p0 = np.array([math.log(5.0), math.log(1.0), math.log(1.0), 1.4, 1.6, 1.8, 20.0, 5.0, 2.0] + ([] if qh_fixed is not None else [0.35]))""",
"""    p0 = np.array([math.log(5.0), math.log(1.0), math.log(1.0), 1.4, 1.6, 1.8] + ([] if c6_fixed is not None else [20.0, 5.0, 2.0]) + ([] if qh_fixed is not None else [0.35]))""")
s = s.replace("for trial in range(12):", "for trial in range(6):")
s = s.replace("""fit('gascharge', qh_gas)
fit('freecharge', None)""", """E6 = [17.48, 6.331, 2.376]
if len(sys.argv) > 1 and sys.argv[1] == 'e6':
    fit('e6disp', qh_gas, E6)
else:
    fit('gascharge', qh_gas)
    fit('freecharge', None)""")
open(p, 'w', encoding='utf-8').write(s)
print('patched')
