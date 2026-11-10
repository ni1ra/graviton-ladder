"""Green identity (5) on explicit non-degenerate q."""
import sympy as sp
import random
x, y, z = sp.symbols('x y z', real=True)
X = [x, y, z]
r2 = x**2 + y**2 + z**2
N = (1 - r2) / 2
random.seed(3)

def monos(deg):
    return [x**a * y**b * z**(deg - a - b) for a in range(deg + 1) for b in range(deg + 1 - a)]

def mono_sphere(a, b, c):
    if a % 2 or b % 2 or c % 2:
        return 0
    n = a + b + c
    return 2 * sp.gamma(sp.Rational(a + 1, 2)) * sp.gamma(sp.Rational(b + 1, 2)) * sp.gamma(sp.Rational(c + 1, 2)) / sp.gamma(sp.Rational(n + 3, 2))

def sphere_int(f):
    f = sp.Poly(sp.expand(f), x, y, z)
    return sp.simplify(sum(co * mono_sphere(*m) for m, co in f.terms()))

def ball_int(f):
    f = sp.Poly(sp.expand(f), x, y, z)
    return sp.simplify(sum(co * mono_sphere(*m) / (sum(m) + 3) for m, co in f.terms()))

def lap(f):
    return sum(sp.diff(f, v, 2) for v in X)

for trial in range(3):
    q = sp.zeros(3, 3)
    for i in range(3):
        for j in range(i, 3):
            p = sum(random.randint(-3, 3) * m for d in range(0, 5) for m in monos(d))
            q[i, j] = p
            q[j, i] = p
    xv = sp.Matrix(X)
    Rl = sum(sp.diff(q[i, j], X[i], X[j]) for i in range(3) for j in range(3)) - lap(q.trace())
    dVol = sp.Rational(1, 2) * ball_int(q.trace())
    dArea = sp.Rational(1, 2) * sphere_int(q.trace() - (xv.T * q * xv)[0])
    lhs = sp.simplify(dArea - 2 * dVol)
    rhs = sp.simplify(-sp.Rational(1, 2) * ball_int(N * Rl))
    print('trial', trial, ' dVol =', dVol, ' dArea =', dArea, ' lhs =', lhs, ' rhs =', rhs, ' equal:', sp.simplify(lhs - rhs) == 0)
