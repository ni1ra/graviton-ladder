"""Post hoc diagnostic (written AFTER d2_killing.py ran and its plant control failed as registered).
Question: on which data does the plant (coefficient 2 -> 3) go uncaught in the wedge, and why?
Same seed and data as d2_killing.py."""
import random
import sympy as sp
import d2_killing as d

rnd = random.Random(72)
N_wedge = d.X[0]
rows = []
for deg in (2, 3, 4, 5):
    for trial in range(2):
        phi = d.random_harmonic(deg, rnd)
        Ks = d.K_scope(phi)
        for name, K in [('scope', Ks), ('hess', d.hess(phi)), ('curl', d.curl(Ks))]:
            K = K.applyfunc(sp.expand)
            if K == sp.zeros(3, 3):
                continue
            G, M, _ = d.rates(N_wedge, K)
            G3, M3, _ = d.rates(N_wedge, K, coeff=3)
            rows.append((f'{name} deg={deg} #{trial}', d.zero(M), d.zero(G), not d.zero(G3 - M3)))
            # keep the random stream aligned with d2_killing.py (it draws a random quadratic lapse per datum)
            _ = [rnd.randint(-3, 3) for _ in range(9 + 3 + 1)]
n = len(rows)
m0 = [r for r in rows if r[1]]
uncaught = [r for r in rows if not r[3]]
print(f'{n} data; wedge modular rate M = 0 on {len(m0)}; plant uncaught on {len(uncaught)}')
for r in uncaught:
    print('  uncaught:', r[0], '| M = 0:', r[1], '| G = 0:', r[2])
ok = all(r[1] and r[2] for r in uncaught) and all(r[3] for r in rows if not r[1])
print('plant caught on every datum with M != 0, and every uncaught datum has M = G = 0:', ok)
