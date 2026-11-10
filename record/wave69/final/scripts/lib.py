"""Shared exact helpers for w69-final. Polynomials in Cartesian x, y, z over the rationals (sympy).

Conventions (Paper A r18, Section 2):
  N = (1 - r^2)/2,  (curl T)_ij = eps_ikl d_k T_lj,  (x cross T)_ij = eps_ikl x_k T_lj,
  E = Ric_3(h),  B = - sym curl K,  flow: dE = N d_t E + 2 sym(x cross B),  dB = N d_t B - 2 sym(x cross E).
"""
import json, os, sys
import sympy as sp

x, y, z = sp.symbols('x y z', real=True)
X = (x, y, z)
r2 = x**2 + y**2 + z**2
N = (1 - r2) / 2
eps = sp.LeviCivita
I3 = sp.eye(3)


def dfact(n):
    """Double factorial with (-1)!! = 1."""
    return sp.Integer(1) if n <= 0 else sp.factorial2(n)


def sphere_mono(a, b, c):
    """Integral of x^a y^b z^c over the unit sphere."""
    if a % 2 or b % 2 or c % 2:
        return sp.Integer(0)
    return 4 * sp.pi * dfact(a - 1) * dfact(b - 1) * dfact(c - 1) / dfact(a + b + c + 1)


def _terms(p):
    p = sp.Poly(sp.expand(p), x, y, z)
    return p.terms()


def sphere_int(p):
    """Exact integral over the unit sphere of a polynomial."""
    return sp.nsimplify(sum(co * sphere_mono(*mon) for mon, co in _terms(p)))


def ball_int(p):
    """Exact integral over the unit ball of a polynomial."""
    return sum(co * sphere_mono(*mon) / (sum(mon) + 3) for mon, co in _terms(p))


def grad(f):
    return [sp.diff(f, v) for v in X]


def lap(f):
    return sum(sp.diff(f, v, 2) for v in X)


def mat(fn):
    return sp.Matrix(3, 3, fn)


def mexp(M):
    return M.applyfunc(sp.expand)


def is_zero(M):
    M = sp.Matrix(M)
    return all(sp.expand(e) == 0 for e in M)


def sym(M):
    return (M + M.T) / 2


def curl(T):
    return mat(lambda i, j: sum(eps(i, k, l) * sp.diff(T[l, j], X[k]) for k in range(3) for l in range(3)))


def xcross(T):
    return mat(lambda i, j: sum(eps(i, k, l) * X[k] * T[l, j] for k in range(3) for l in range(3)))


def div(T):
    return sp.Matrix([sum(sp.diff(T[i, j], X[i]) for i in range(3)) for j in range(3)])


def tlap(T):
    return T.applyfunc(lap)


def contract(A, B):
    return sum(A[i, j] * B[i, j] for i in range(3) for j in range(3))


def ric3(h):
    """Linearised Ricci tensor of the 3-metric delta + h."""
    tr = h.trace()
    d = div(h)
    return mat(lambda i, j: sp.Rational(1, 2) * (-lap(h[i, j]) + sp.diff(d[j], X[i]) + sp.diff(d[i], X[j])
                                                  - sp.diff(tr, X[i], X[j])))


def rlin(q):
    """Linearised scalar curvature d_i d_j q_ij - Lap q_kk."""
    return sum(sp.diff(q[i, j], X[i], X[j]) for i in range(3) for j in range(3)) - lap(q.trace())


def inc(f):
    return mat(lambda i, j: sum(eps(i, k, l) * eps(j, m, n) * sp.diff(f[l, n], X[k], X[m])
                                for k in range(3) for l in range(3) for m in range(3) for n in range(3)))


def cotton(f):
    """wave69's Cotton-York operator C f = curl(inc f - (1/2) delta tr inc f)."""
    g = inc(f)
    return curl(g - sp.Rational(1, 2) * I3 * g.trace())


def harmonics(l):
    """Three real harmonic polynomials of degree l: sectoral cos, sectoral sin, zonal; and one mixture."""
    c = sp.expand(sp.re(sp.expand((x + sp.I * y)**l)))
    s = sp.expand(sp.im(sp.expand((x + sp.I * y)**l)))
    # zonal: r^l P_l(z/r) as a polynomial
    t = sp.symbols('t')
    P = sp.legendre(l, t)
    zon = sp.expand(sum(co * z**k * r2**sp.Rational(l - k, 2) for (k,), co in sp.Poly(P, t).terms()))
    # tesseral m = 1 by differentiating: (x d_z - z d_x) applied to the zonal one stays harmonic of degree l
    tes = sp.expand(x * sp.diff(zon, z) - z * sp.diff(zon, x))
    mix = sp.expand(c + 2 * zon - 3 * tes)
    return {'m=l cos': c, 'm=l sin': s, 'm=0': zon, 'mixed': mix}


def magnetic(phi):
    """v = x cross grad(phi), K = sym grad v."""
    g = grad(phi)
    v = sp.Matrix([sum(eps(i, k, l) * X[k] * g[l] for k in range(3) for l in range(3)) for i in range(3)])
    K = mat(lambda i, j: sp.Rational(1, 2) * (sp.diff(v[j], X[i]) + sp.diff(v[i], X[j])))
    return v, mexp(K)


def xvec():
    return sp.Matrix(X)


class Report:
    def __init__(self, name):
        self.name = name
        self.checks = []
        self.plants = []
        self.values = {}

    def check(self, label, ok, detail=None):
        self.checks.append({'label': label, 'pass': bool(ok), 'detail': None if detail is None else str(detail)})
        if not ok:
            print('FAIL', label, detail)

    def plant(self, label, caught, detail=None):
        """A planted error: 'caught' must be True (the check that guards it fails on the plant)."""
        self.plants.append({'label': label, 'caught': bool(caught), 'detail': None if detail is None else str(detail)})
        if not caught:
            print('PLANT NOT CAUGHT', label, detail)

    def value(self, key, val):
        self.values[key] = str(val)

    def write(self):
        out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'out', self.name + '.json')
        ok = all(c['pass'] for c in self.checks) and all(p['caught'] for p in self.plants)
        json.dump({'script': self.name, 'ALL_PASS': ok, 'n_checks': len(self.checks), 'n_plants': len(self.plants),
                   'checks': self.checks, 'plants': self.plants, 'values': self.values},
                  open(out, 'w'), indent=1)
        print(self.name, 'ALL_PASS' if ok else 'FAILED', len(self.checks), 'checks,', len(self.plants), 'plants')
        if not ok:
            sys.exit(1)
