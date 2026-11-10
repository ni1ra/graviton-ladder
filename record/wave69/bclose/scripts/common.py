"""Shared exact helpers for w69-bclose.

Monomial integration over the unit sphere and ball follows wave68/q0edge/qcommon.py (closed
double-factorial formula, no nsimplify); the magnetic family follows wave64/btheorem/scripts/s9_scope.py.
All exact integrals are returned as rationals: the value is (returned number) * pi.
"""
import json, os
from fractions import Fraction
import sympy as sp

t, x, y, z = sp.symbols('t x y z', real=True)
X = (x, y, z)
X4 = (t, x, y, z)
r2 = x**2 + y**2 + z**2
N = (1 - r2)/2
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), 'out')


def _df(n):
    r = 1
    while n > 1:
        r *= n
        n -= 2
    return r


def mono_sphere(a, b, c):
    """oint x^a y^b z^c dA / pi over the unit sphere."""
    if a % 2 or b % 2 or c % 2:
        return Fraction(0)
    return Fraction(4*_df(a - 1)*_df(b - 1)*_df(c - 1), _df(a + b + c + 1))


def mono_ball(a, b, c):
    return mono_sphere(a, b, c)/(a + b + c + 3)


# ---- dict polynomials {(i,j,k): Fraction} in x, y, z ------------------------------------------
def pd(expr):
    """sympy polynomial in x,y,z (rational coefficients) -> dict."""
    p = sp.Poly(sp.expand(expr), x, y, z)
    return {m: Fraction(int(sp.numer(c)), int(sp.denom(c))) for m, c in zip(p.monoms(), p.coeffs()) if c != 0}


def pd_diff(p, d):
    out = {}
    for (i, j, k), c in p.items():
        e = [i, j, k]
        f = c
        ok = True
        for ax in range(3):
            for _ in range(d[ax]):
                if e[ax] == 0:
                    ok = False
                    break
                f *= e[ax]
                e[ax] -= 1
            if not ok:
                break
        if ok:
            out[tuple(e)] = out.get(tuple(e), 0) + f
    return out


def pd_int_prod(p, q, dom):
    """Integral of p*q over the sphere ('S') or ball ('B'), divided by pi."""
    f = mono_sphere if dom == 'S' else mono_ball
    tot = Fraction(0)
    for (a, b, c), u in p.items():
        for (d, e, g), v in q.items():
            tot += u*v*f(a + d, b + e, c + g)
    return tot


def sphere(expr):
    return sum((c*mono_sphere(*m) for m, c in pd(expr).items()), Fraction(0))


def ball(expr):
    return sum((c*mono_ball(*m) for m, c in pd(expr).items()), Fraction(0))


# ---- the scope family ---------------------------------------------------------------------------
def zonal(l):
    u = sp.Symbol('u')
    P = sp.Poly(sp.legendre(l, u), u)
    return sp.expand(sum(co*z**k*r2**((l - k)//2) for (k,), co in zip(P.monoms(), P.coeffs())))


def sectoral(l):
    return sp.expand(sp.re(sp.expand((x + sp.I*y)**l)))


def magnetic_K(phi):
    g = [sp.diff(phi, s) for s in X]
    V = [y*g[2] - z*g[1], z*g[0] - x*g[2], x*g[1] - y*g[0]]
    return [[sp.expand((sp.diff(V[j], X[i]) + sp.diff(V[i], X[j]))/2) for j in range(3)] for i in range(3)]


def members():
    out = []
    for l in (2, 3, 4):
        out.append((f'l={l} m={l}', l, sectoral(l)))
        out.append((f'l={l} m=0', l, zonal(l)))
    return out


def K2(K):
    return sp.expand(sum(K[i][j]**2 for i in range(3) for j in range(3)))


def lap(f):
    return sp.expand(sum(sp.diff(f, s, 2) for s in X))


def inv_lap(p):
    """A polynomial u with lap(u) = p (exact), by the homogeneous-component formula
    u = sum_k (-1)^k r^(2k+2) lap^k(p_n) / prod_{j=0..k} (2j+2)(2n-2j+3)  for p_n homogeneous of degree n."""
    P = sp.Poly(sp.expand(p), x, y, z)
    comps = {}
    for m, c in zip(P.monoms(), P.coeffs()):
        comps.setdefault(sum(m), 0)
        comps[sum(m)] += c*x**m[0]*y**m[1]*z**m[2]
    u = 0
    for n, pn in comps.items():
        cur = pn
        den = sp.Integer(1)
        k = 0
        while cur != 0:
            den *= (2*k + 2)*(2*n - 2*k + 3)
            u += (-1)**k*r2**(k + 1)*cur/den
            cur = lap(cur)
            k += 1
    u = sp.expand(u)
    assert sp.expand(lap(u) - p) == 0
    return u


def R_lin(q):
    tr = q[0][0] + q[1][1] + q[2][2]
    s = sum(sp.diff(q[i][j], X[i], X[j]) for i in range(3) for j in range(3))
    return sp.expand(s - lap(tr))


def symgrad(xi):
    return [[sp.expand(sp.diff(xi[j], X[i]) + sp.diff(xi[i], X[j])) for j in range(3)] for i in range(3)]


def dump(name, out):
    os.makedirs(OUT, exist_ok=True)
    checks = {k: v for k, v in out.items() if k.startswith('CHECK')}
    out['ALL_PASS'] = all(bool(v) for v in checks.values())
    with open(os.path.join(OUT, name + '.json'), 'w') as fh:
        json.dump(out, fh, indent=1, default=str)
    for k, v in checks.items():
        if not v:
            print('FAILED', k)
    print(name, 'checks:', len(checks), 'ALL_PASS:', out['ALL_PASS'])
