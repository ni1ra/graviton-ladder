"""Killing-limit test with the plant control stated correctly. Written AFTER d2_killing.py ran.

d2_killing.py is kept unchanged as registered. Its plant control ("the plant must break the wedge
equality on every datum") fails as written, because on 16 of its 24 data both rates vanish on the
wedge (curl K = 0 or curl curl K = 0), so no coefficient can show there. This script reruns the same
data and the same identities, and states the control as it should have been registered:
the plant must be caught on every datum whose wedge modular rate is non-zero, and at least one
such datum must exist per degree.
Run: python3 d2b_killing_fixed_control.py   (exit code 1 on any FAIL)
"""
import sys, random
import sympy as sp
import d2_killing as d

rnd = random.Random(72)
N_wedge = d.X[0]
N_ball = (1 - sum(v**2 for v in d.X)) / 2
wedge_ok = ball_ok = gen_ok = plant_ok = True
live = {}
n = 0
for deg in (2, 3, 4, 5):
    for trial in range(2):
        phi = d.random_harmonic(deg, rnd)
        Ks = d.K_scope(phi)
        for name, K in [('scope', Ks), ('hess', d.hess(phi)), ('curl', d.curl(Ks))]:
            K = K.applyfunc(sp.expand)
            if K == sp.zeros(3, 3):
                continue
            assert d.is_tt_harmonic(K)
            n += 1
            G, M, rhs = d.rates(N_wedge, K)
            wedge_ok &= d.zero(G - M)
            if not d.zero(M):
                live[deg] = live.get(deg, 0) + 1
                G3, M3, _ = d.rates(N_wedge, K, coeff=3)
                plant_ok &= not d.zero(G3 - M3)
            G, M, rhs = d.rates(N_ball, K)
            ball_ok &= d.zero(G - M - K) and d.zero(G - M - rhs)
            Nq = sum(rnd.randint(-3, 3) * a * b for a in d.X for b in d.X) + sum(rnd.randint(-3, 3) * a for a in d.X) + rnd.randint(-3, 3)
            G, M, rhs = d.rates(Nq, K)
            gen_ok &= d.zero(G - M - rhs)

print(f'{n} data; data with non-zero wedge modular rate, by degree: {live}')
d.FAILS.clear()
d.check('wedge (Killing lapse): geometric rate = modular rate on every datum', wedge_ok)
d.check('ball: geometric rate - modular rate = K on every datum', ball_ok)
d.check('general lapse: G - M = 2 sym(K.HessN) - (Lap N) K', gen_ok)
d.check('every degree 2..5 has a datum with non-zero wedge modular rate', all(live.get(k, 0) > 0 for k in (2, 3, 4, 5)))
d.check('plant (coefficient 2 -> 3) caught on every datum with non-zero wedge modular rate', plant_ok)
print('failed:', d.FAILS)
sys.exit(1 if d.FAILS else 0)
