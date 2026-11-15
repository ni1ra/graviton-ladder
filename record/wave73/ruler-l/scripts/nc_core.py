"""Shared exact machinery for lane NULL-CONE (wave 73).

Flat ball of radius 1, diamond tips at t = -1 and t = +1.  The cone is the future
light cone of the past tip: points (t, x) = (v - 1, v n), |n| = 1, 0 <= v <= 1,
generators k = d_t + d_r, affine parameter v, background expansion 2/v.

Everything is a polynomial with rational coefficients in (v, n1, n2, n3); the sphere
is handled by reducing modulo n1^2 + n2^2 + n3^2 - 1 and integrating monomials exactly.
Units: the common factor 1/(8G) is dropped, as in D-DEP's report.
"""
import itertools

import sympy as sp

t, x, y, z = sp.symbols("t x y z", real=True)
X4 = (t, x, y, z)
XS = (x, y, z)
v = sp.Symbol("v", positive=True)
n1, n2, n3 = sp.symbols("n1 n2 n3", real=True)
NV = (n1, n2, n3)
ETA = sp.diag(-1, 1, 1, 1)
SPHERE = sp.Poly(n1**2 + n2**2 + n3**2 - 1, n3)


def dfact(k):
    """(k)!! with (-1)!! = 1."""
    return sp.factorial2(k) if k > 0 else sp.Integer(1)


def red(e):
    """Normal form of a polynomial on the unit sphere (degree <= 1 in n3)."""
    e = sp.expand(e)
    if e == 0:
        return sp.Integer(0)
    return sp.expand(sp.Poly(e, n3).rem(SPHERE).as_expr())


def sphere_int(e):
    """Exact integral over the unit sphere of a polynomial in n1, n2, n3."""
    e = sp.expand(e)
    if e == 0:
        return sp.Integer(0)
    tot = sp.Integer(0)
    for (a, b, c), co in sp.Poly(e, n1, n2, n3).terms():
        if a % 2 or b % 2 or c % 2:
            continue
        tot += co * 4 * sp.pi * dfact(a - 1) * dfact(b - 1) * dfact(c - 1) / dfact(a + b + c + 1)
    return sp.expand(tot)


def ball_int(e):
    """Exact integral over the unit ball of a polynomial in x, y, z."""
    e = sp.expand(e)
    if e == 0:
        return sp.Integer(0)
    tot = sp.Integer(0)
    for (a, b, c), co in sp.Poly(e, x, y, z).terms():
        if a % 2 or b % 2 or c % 2:
            continue
        tot += co * 4 * sp.pi * dfact(a - 1) * dfact(b - 1) * dfact(c - 1) / dfact(a + b + c + 1) / (a + b + c + 3)
    return sp.expand(tot)


def surf_int(e):
    """Exact integral over the unit sphere r = 1 of a polynomial in x, y, z."""
    return sphere_int(sp.expand(e).subs({x: n1, y: n2, z: n3}, simultaneous=True))


# ---------------------------------------------------------------- harmonics and metrics
def harmonic(l):
    """A generic real harmonic polynomial of degree l (three m-values mixed)."""
    w = sp.expand((x + sp.I * y) ** l)
    w1 = sp.expand((x + sp.I * y) ** (l - 1))
    r2 = x**2 + y**2 + z**2
    zonal = sp.expand(sum(sp.Rational((-1) ** k) * sp.binomial(l, k) * sp.binomial(2 * l - 2 * k, l)
                          * z ** (l - 2 * k) * r2**k for k in range(l // 2 + 1)) / 2**l)
    phi = sp.re(w) + 2 * z * sp.im(w1) + 3 * zonal
    phi = sp.expand(phi)
    assert sp.expand(sum(sp.diff(phi, c, 2) for c in XS)) == 0
    return phi


def lap(e):
    return sp.expand(sum(sp.diff(e, c, 2) for c in XS))


def grad(e):
    return [sp.diff(e, c) for c in XS]


def hess(e):
    return sp.Matrix(3, 3, lambda i, j: sp.diff(e, XS[i], XS[j]))


def metric_electric(phi):
    """Static electric tide: h_tt = -2 phi, h_ij = -2 phi delta_ij."""
    h = sp.zeros(4, 4)
    h[0, 0] = -2 * phi
    for i in range(1, 4):
        h[i, i] = -2 * phi
    return h


def vec_potential(phi):
    """V = x cross grad(phi): tangential, divergence-free, harmonic of the same degree."""
    g = grad(phi)
    return [sp.expand(sum(sp.LeviCivita(i, j, k) * XS[j] * g[k] for j in range(3) for k in range(3)))
            for i in range(3)]


def sym_grad(V):
    return sp.Matrix(3, 3, lambda i, j: sp.expand((sp.diff(V[i], XS[j]) + sp.diff(V[j], XS[i])) / 2))


def tt_series(H, parity):
    """Transverse-traceless gauge solution of the wave equation with polynomial data H:
    parity 0: h(0) = H, dh/dt(0) = 0;  parity 1: h(0) = 0, dh/dt(0) = 2 H  (so K = H)."""
    out = sp.zeros(3, 3)
    term, k = H, 0
    while any(e != 0 for e in term):
        if parity == 0:
            out += term * t ** (2 * k) / sp.factorial(2 * k)
        else:
            out += 2 * term * t ** (2 * k + 1) / sp.factorial(2 * k + 1)
        term = term.applyfunc(lap)
        k += 1
    h = sp.zeros(4, 4)
    for i in range(3):
        for j in range(3):
            h[i + 1, j + 1] = sp.expand(out[i, j])
    return h


def rung_scalar(phi, l):
    """Scalar-type transverse traceless polynomial tensor of degree l (one rung above Hess phi)."""
    r2 = x**2 + y**2 + z**2
    g = grad(phi)
    b = -sp.Rational(6 * (l - 1), l + 9)
    c = sp.Rational(4 * l * (l - 1), l + 9)
    H = sp.Matrix(3, 3, lambda i, j: r2 * sp.diff(phi, XS[i], XS[j]) + b * (XS[i] * g[j] + XS[j] * g[i])
                  + (c * phi if i == j else 0))
    return H.applyfunc(sp.expand)


def rung_vector(phi, l):
    """Vector-type transverse traceless polynomial tensor of degree l+1 (one rung above sym grad V)."""
    r2 = x**2 + y**2 + z**2
    V = vec_potential(phi)
    K = sym_grad(V)
    b = -sp.Rational(l - 1, l + 4)
    H = sp.Matrix(3, 3, lambda i, j: r2 * K[i, j] + b * (XS[i] * V[j] + XS[j] * V[i]))
    return H.applyfunc(sp.expand)


def is_tt(H):
    tr = sp.expand(H.trace())
    div = [sp.expand(sum(sp.diff(H[i, j], XS[i]) for i in range(3))) for j in range(3)]
    return tr == 0 and all(d == 0 for d in div)


# ---------------------------------------------------------------- linearised curvature
class Curv:
    """Linearised Riemann tensor R_abcd of a metric perturbation h (4x4, polynomial)."""

    def __init__(self, h):
        self.h = h
        self.d2 = {}
        for a in range(4):
            for b in range(a, 4):
                for c in range(4):
                    for d in range(c, 4):
                        e = sp.diff(h[a, b], X4[c], X4[d]) if h[a, b] != 0 else sp.Integer(0)
                        for (p, q) in ((a, b), (b, a)):
                            for (r, s) in ((c, d), (d, c)):
                                self.d2[(p, q, r, s)] = e

    def R(self, a, b, c, d):
        D = self.d2
        return (D[(a, d, b, c)] + D[(b, c, a, d)] - D[(a, c, b, d)] - D[(b, d, a, c)]) / 2

    def ricci_zero(self):
        for b in range(4):
            for d in range(b, 4):
                e = sum(ETA[a, a] * self.R(a, b, a, d) for a in range(4))
                if sp.expand(e) != 0:
                    return False
        return True

    def E(self):
        return sp.Matrix(3, 3, lambda i, j: sp.expand(self.R(0, i + 1, 0, j + 1)))

    def B(self):
        return sp.Matrix(3, 3, lambda i, j: sp.expand(sum(
            sp.LeviCivita(i, k, m) * self.R(k + 1, m + 1, 0, j + 1) for k in range(3) for m in range(3)) / 2))

    def cone_weyl(self, nomag=False):
        """C_ij = R(k, e_i, k, e_j) on the cone, projected transverse and trace-free, as a 3x3
        matrix of polynomials in (v, n).  nomag (plant) drops the mixed time-space part."""
        kv = (sp.Integer(1),) + NV
        on = {t: v - 1, x: v * n1, y: v * n2, z: v * n3}
        M = sp.zeros(3, 3)
        for i in range(3):
            for j in range(i, 3):
                e = sp.Integer(0)
                for a in range(4):
                    for b in range(4):
                        if nomag and ((a == 0) != (b == 0)):
                            continue
                        r = self.R(a, i + 1, b, j + 1)
                        if r != 0:
                            e += kv[a] * kv[b] * r
                e = red(sp.expand(e).subs(on, simultaneous=True))
                M[i, j] = M[j, i] = e
        return tt_project(M)


def tt_project(M):
    nv = sp.Matrix(NV)
    P = sp.eye(3) - nv * nv.T
    MT = (P * M * P).applyfunc(red)
    tr = red(MT.trace())
    return (MT - P * tr / 2).applyfunc(red)


def moment(C, p):
    """sigma_p(v) = -v^(-p) int_0^v u^p C(u) du, entrywise on polynomials in v.
    p = 2 is the shear of the cone (Sachs equation with expansion 2/v)."""
    def one(e):
        e = sp.expand(e)
        if e == 0:
            return sp.Integer(0)
        out = sp.Integer(0)
        for (k,), co in sp.Poly(e, v).terms():
            out += -co * v ** (k + 1) / (k + p + 1)
        return sp.expand(out)
    return C.applyfunc(one)


def norm2(A, B=None):
    """Pointwise contraction A_ij B_ij reduced on the sphere."""
    B = A if B is None else B
    return red(sum(A[i, j] * B[i, j] for i in range(3) for j in range(3)))


def cone_int(e, weight):
    """Integral over the sphere and over 0 <= v <= 1 of weight(v) * e(v, n)."""
    e = sp.expand(weight * e)
    if e == 0:
        return sp.Integer(0)
    tot = sp.Integer(0)
    for (k,), co in sp.Poly(e, v).terms():
        tot += co / (k + 1)
    return sphere_int(tot)


def ang_int_poly(e):
    """Sphere integral keeping the dependence on v."""
    e = sp.expand(e)
    if e == 0:
        return sp.Integer(0)
    tot = sp.Integer(0)
    for (k,), co in sp.Poly(e, v).terms():
        tot += sphere_int(co) * v**k
    return sp.expand(tot)


F_CONE = v * (1 - v)            # zeta = f(v) k on the cone, R = 1


def chain_fields_on_cone(cv):
    """chi_E = v * (x.E.x), chi_B = v * (x.B.x) restricted to the cone."""
    on = {t: v - 1, x: v * n1, y: v * n2, z: v * n3}
    E, B = cv.E(), cv.B()
    pe = sp.expand(sum(XS[i] * XS[j] * E[i, j] for i in range(3) for j in range(3)))
    pb = sp.expand(sum(XS[i] * XS[j] * B[i, j] for i in range(3) for j in range(3)))
    return red(v * pe.subs(on, simultaneous=True)), red(v * pb.subs(on, simultaneous=True))


def hislop_longo_E1(psi):
    """E_1[psi] = (1/2) int_B N (psi_t^2 + |grad psi|^2) + (1/2) int_B psi^2 at t = 0, N = (1 - r^2)/2."""
    N = (1 - x**2 - y**2 - z**2) / 2
    pt = sp.diff(psi, t).subs(t, 0)
    p0 = psi.subs(t, 0)
    dens = pt**2 + sum(sp.diff(p0, c) ** 2 for c in XS)
    return ball_int(N * dens) / 2 + ball_int(p0**2) / 2


def entropy_slice(cv, l):
    """The programme's entropy: c_l (E_1[x.E.x] + E_1[x.B.x]), c_l = 4 (l-2)!/(l+2)!  (1/8G dropped)."""
    E, B = cv.E(), cv.B()
    pe = sp.expand(sum(XS[i] * XS[j] * E[i, j] for i in range(3) for j in range(3)))
    pb = sp.expand(sum(XS[i] * XS[j] * B[i, j] for i in range(3) for j in range(3)))
    cl = 4 * sp.factorial(l - 2) / sp.factorial(l + 2)
    return sp.expand(cl * (hislop_longo_E1(pe) + hislop_longo_E1(pb)))


def spatial_K(h):
    """K_ij = (1/2) dh_ij/dt at t = 0 and h_ij at t = 0."""
    K = sp.Matrix(3, 3, lambda i, j: sp.expand(sp.diff(h[i + 1, j + 1], t).subs(t, 0) / 2))
    h0 = sp.Matrix(3, 3, lambda i, j: sp.expand(h[i + 1, j + 1].subs(t, 0)))
    return K, h0


def rstr(e):
    """Exact string of a rational multiple of pi (no numerical simplification anywhere)."""
    return str(sp.factor(sp.expand(e)))
