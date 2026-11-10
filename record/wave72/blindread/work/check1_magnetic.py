"""Blind check 1: magnetic family, rates, identities, integrals. Written from the paper's text only."""
import sympy as sp
import random
from itertools import product

x, y, z, t, eps = sp.symbols('x y z t epsilon', real=True)
X = [x, y, z]
r2 = x**2 + y**2 + z**2
Nball = (1 - r2) / 2
LC = sp.LeviCivita
random.seed(7)


def monos(deg):
    return [x**a * y**b * z**(deg - a - b) for a in range(deg + 1) for b in range(deg + 1 - a)]


def lap(f):
    return sum(sp.diff(f, v, 2) for v in X)


def random_harmonic(l):
    ms = monos(l)
    cs = sp.symbols('c0:%d' % len(ms))
    f = sum(c * m for c, m in zip(cs, ms))
    eqs = sp.Poly(lap(f), x, y, z).coeffs() if lap(f) != 0 else []
    sol = sp.solve(eqs, cs, dict=True)[0] if eqs else {}
    f = f.subs(sol)
    free = [c for c in cs if f.has(c)]
    f = f.subs({c: random.randint(-5, 5) or 1 for c in free})
    assert sp.expand(lap(f)) == 0
    return sp.expand(f)


def ball_int(f):
    """exact integral of polynomial over the unit ball"""
    f = sp.Poly(sp.expand(f), x, y, z)
    tot = 0
    for (a, b, c), co in f.terms():
        if a % 2 or b % 2 or c % 2:
            continue
        n = a + b + c
        # int over sphere of x^a y^b z^c = 2 G((a+1)/2)G((b+1)/2)G((c+1)/2)/G((n+3)/2)
        ang = 2 * sp.gamma(sp.Rational(a + 1, 2)) * sp.gamma(sp.Rational(b + 1, 2)) * sp.gamma(sp.Rational(c + 1, 2)) / sp.gamma(sp.Rational(n + 3, 2))
        tot += co * ang / (n + 3)
    return sp.simplify(tot)


def symm(M):
    return (M + M.T) / 2


def grad(v):  # (grad v)_ij = d_i v_j
    return sp.Matrix(3, 3, lambda i, j: sp.diff(v[j], X[i]))


def curlT(K):  # (curl K)_ij = eps_ikl d_k K_lj
    return sp.Matrix(3, 3, lambda i, j: sum(LC(i, k, l) * sp.diff(K[l, j], X[k]) for k in range(3) for l in range(3)))


def crossT(a, T):  # (a x T)_ij = eps_ikl a_k T_lj
    return sp.Matrix(3, 3, lambda i, j: sum(LC(i, k, l) * a[k] * T[l, j] for k in range(3) for l in range(3)))


def ric_lin(q):
    tr = q.trace()
    return sp.Matrix(3, 3, lambda i, j: sp.Rational(1, 2) * (
        sum(sp.diff(q[j, k], X[k], X[i]) + sp.diff(q[i, k], X[k], X[j]) for k in range(3))
        - lap(q[i, j]) - sp.diff(tr, X[i], X[j])))


def R_lin(q):
    return sum(sp.diff(q[i, j], X[i], X[j]) for i in range(3) for j in range(3)) - lap(q.trace())


def zero(M):
    return all(sp.expand(e) == 0 for e in M)


def riemann_lin(h):
    """4D linearised Riemann R_abcd, signature -+++, coordinates (t,x,y,z); h is 4x4 lower-index."""
    C = [t, x, y, z]
    def R(a, b, c, d):
        return sp.Rational(1, 2) * (sp.diff(h[a, d], C[b], C[c]) + sp.diff(h[b, c], C[a], C[d])
                                    - sp.diff(h[a, c], C[b], C[d]) - sp.diff(h[b, d], C[a], C[c]))
    return R


def EB_from_R(R):
    E = sp.Matrix(3, 3, lambda i, j: R(0, i + 1, 0, j + 1))
    B = sp.Matrix(3, 3, lambda i, j: sp.Rational(1, 2) * sum(LC(i, k, l) * R(0, j + 1, k + 1, l + 1) for k in range(3) for l in range(3)))
    return E, B


def lie_R(R, zeta):
    C = [t, x, y, z]
    def LR(a, b, c, d):
        s = sum(zeta[e] * sp.diff(R(a, b, c, d), C[e]) for e in range(4))
        s += sum(R(e, b, c, d) * sp.diff(zeta[e], C[a]) for e in range(4))
        s += sum(R(a, e, c, d) * sp.diff(zeta[e], C[b]) for e in range(4))
        s += sum(R(a, b, e, d) * sp.diff(zeta[e], C[c]) for e in range(4))
        s += sum(R(a, b, c, e) * sp.diff(zeta[e], C[d]) for e in range(4))
        return s
    return LR


zeta = [sp.Rational(1, 2) * (1 - t**2 - r2), -t * x, -t * y, -t * z]

for l in [2, 3, 4]:
    print('==== l =', l)
    phi = random_harmonic(l)
    gphi = [sp.diff(phi, v) for v in X]
    v = sp.Matrix([sum(LC(i, k, m) * X[k] * gphi[m] for k in range(3) for m in range(3)) for i in range(3)])
    K = symm(grad(v)).applyfunc(sp.expand)
    xv = sp.Matrix(X)
    print(' traceless:', sp.expand(K.trace()) == 0)
    print(' transverse:', zero(sp.Matrix([sum(sp.diff(K[i, j], X[i]) for i in range(3)) for j in range(3)])))
    print(' harmonic:', zero(K.applyfunc(lap)))
    print(' x.K.x = 0:', sp.expand((xv.T * K * xv)[0]) == 0)
    cK = curlT(K)
    print(' curl K symmetric:', zero(cK - cK.T))
    # 4D check of E, B from h_ij = 2 eps t K_ij
    h = sp.zeros(4, 4)
    for i in range(3):
        for j in range(3):
            h[i + 1, j + 1] = 2 * eps * t * K[i, j]
    R = riemann_lin(h)
    E4, B4 = EB_from_R(R)
    print(' E = 0 from 4D Riemann:', zero(E4))
    print(' B = -eps curl K from 4D Riemann (paper definition of B):', zero(B4 + eps * cK))
    B = (-eps * cK).applyfunc(sp.expand)
    hess = sp.Matrix(3, 3, lambda i, j: sp.diff(phi, X[i], X[j]))
    print(' B = eps (l+1)/2 dd phi:', zero(B - eps * sp.Rational(l + 1, 2) * hess))
    psiB = sp.expand((xv.T * B * xv)[0])
    print(' psi_B = (1/2)(l-1)l(l+1) eps phi:', sp.expand(psiB - sp.Rational((l - 1) * l * (l + 1), 2) * eps * phi) == 0)
    rate = 2 * symm(crossT(X, B))
    print(' 2 sym(x x B) = (l+1) eps K:', zero(rate - (l + 1) * eps * K))
    print(' Ric_lin(2 N K) = (l+2) K:', zero(ric_lin(2 * Nball * K) - (l + 2) * K))
    print(' R_lin(2 N K) = 0 in the bulk:', sp.expand(R_lin(2 * Nball * K)) == 0)
    # Lie derivative check of eq (4) on this solution at t=0
    LR = lie_R(R, zeta)
    dE, dB = EB_from_R(LR)
    dE0 = dE.subs(t, 0); dB0 = dB.subs(t, 0)
    print(' Lie_zeta E at t=0 equals 2 sym(x x B):', zero(dE0 - rate))
    print(' Lie_zeta B at t=0 equals 0 (static, E=0):', zero(dB0))
    # A and (2)
    A = ball_int(Nball * sum(K[i, j]**2 for i in range(3) for j in range(3)))
    Iphi = ball_int(phi**2)
    print(' A / int phi^2 =', sp.nsimplify(A / Iphi), ' expected (l-1)l(l+1)/2 =', sp.Rational((l - 1) * l * (l + 1), 2))
    cl = sp.Rational(4) * sp.factorial(l - 2) / sp.factorial(l + 2)
    lhs = cl * sp.Rational(l + 2, 2) * ball_int((psiB / eps)**2)
    print(' c_l (l+2)/2 int psi_B^2 == eps^2 A :', sp.simplify(lhs - A) == 0)
    # int N |grad psi|^2 = l int psi^2
    print(' int N|grad phi|^2 = l int phi^2:', sp.simplify(ball_int(Nball * sum(g**2 for g in gphi)) - l * Iphi) == 0)

    # electric solution: Newtonian metric, check E and Lie derivative for B
    U = random_harmonic(l)
    h = sp.zeros(4, 4)
    h[0, 0] = -2 * eps * U
    for i in range(3):
        h[i + 1, i + 1] = -2 * eps * U
    R = riemann_lin(h)
    E4, B4 = EB_from_R(R)
    hessU = sp.Matrix(3, 3, lambda i, j: sp.diff(U, X[i], X[j]))
    print(' electric: E = dd Phi:', zero(E4 - eps * hessU), ' B = 0:', zero(B4))
    LR = lie_R(R, zeta)
    dE, dB = EB_from_R(LR)
    print(' electric: Lie_zeta E at t=0 = 0:', zero(dE.subs(t, 0)),
          '; Lie_zeta B at t=0 = -2 sym(x x E):', zero(dB.subs(t, 0) + 2 * symm(crossT(X, E4))))
    # -2 sym(x x E) = -2 K[U]
    gU = [sp.diff(U, vv) for vv in X]
    vU = sp.Matrix([sum(LC(i, k, m) * X[k] * gU[m] for k in range(3) for m in range(3)) for i in range(3)])
    KU = symm(grad(vU))
    print(' electric: -2 sym(x x E) = -2 eps K[U]:', zero(-2 * symm(crossT(X, E4)) + 2 * eps * KU))

# ---- identity (6) for a general lapse and general TT harmonic K
print('==== identity (6), random lapse, random TT harmonic K')
def random_TT_harmonic(deg):
    ms = monos(deg)
    idx = [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]
    cs = sp.symbols('k0:%d' % (6 * len(ms)))
    comp = {}
    n = 0
    for (i, j) in idx:
        comp[(i, j)] = sum(cs[n + a] * ms[a] for a in range(len(ms)))
        n += len(ms)
    Kg = sp.Matrix(3, 3, lambda i, j: comp[(min(i, j), max(i, j))])
    eqs = []
    eqs += sp.Poly(Kg.trace(), x, y, z).coeffs()
    for j in range(3):
        d = sp.expand(sum(sp.diff(Kg[i, j], X[i]) for i in range(3)))
        if d != 0:
            eqs += sp.Poly(d, x, y, z).coeffs()
    for (i, j) in idx:
        d = sp.expand(lap(Kg[i, j]))
        if d != 0:
            eqs += sp.Poly(d, x, y, z).coeffs()
    sol = sp.solve(eqs, cs, dict=True)[0]
    Kg = Kg.subs(sol)
    free = sorted(Kg.free_symbols - {x, y, z}, key=str)
    Kg = Kg.subs({c: random.randint(-4, 4) or 1 for c in free})
    return Kg.applyfunc(sp.expand), len(free)

for deg in [1, 2, 3]:
    Kg, nfree = random_TT_harmonic(deg)
    Ngen = sum(random.randint(-3, 3) * m for d in range(0, 4) for m in monos(d))
    gN = [sp.diff(Ngen, vv) for vv in X]
    hN = sp.Matrix(3, 3, lambda i, j: sp.diff(Ngen, X[i], X[j]))
    lhs = ric_lin(2 * Ngen * Kg) - 2 * symm(crossT(gN, curlT(Kg)))
    rhs = 2 * symm(Kg * hN) - lap(Ngen) * Kg
    print(' deg', deg, ' (TT harmonic space dim %d): identity (6) holds:' % nfree, zero(lhs - rhs))
    # Rindler: N = x
    lhsR = ric_lin(2 * x * Kg) - 2 * symm(crossT([1, 0, 0], curlT(Kg)))
    print('    Rindler N=x: Ric(2NK) - 2 sym(grad N x curl K) = 0:', zero(lhsR))
    # ball
    lhsB = ric_lin(2 * Nball * Kg) - 2 * symm(crossT([-x, -y, -z], curlT(Kg)))
    print('    ball: difference equals K:', zero(lhsB - Kg))

# ---- Green identity (5) for random polynomial q
print('==== Green identity (5)')
q = sp.Matrix(3, 3, lambda i, j: 0)
for i in range(3):
    for j in range(i, 3):
        p = sum(random.randint(-3, 3) * m for d in range(0, 4) for m in monos(d))
        q[i, j] = p; q[j, i] = p
# delta Vol = 1/2 int tr q ; delta Area = 1/2 oint (tr q - q_nn)
def sphere_int(f):
    f = sp.Poly(sp.expand(f), x, y, z)
    tot = 0
    for (a, b, c), co in f.terms():
        if a % 2 or b % 2 or c % 2:
            continue
        n = a + b + c
        tot += co * 2 * sp.gamma(sp.Rational(a + 1, 2)) * sp.gamma(sp.Rational(b + 1, 2)) * sp.gamma(sp.Rational(c + 1, 2)) / sp.gamma(sp.Rational(n + 3, 2))
    return sp.simplify(tot)
xv = sp.Matrix(X)
dVol = sp.Rational(1, 2) * ball_int(q.trace())
dArea = sp.Rational(1, 2) * sphere_int(q.trace() - (xv.T * q * xv)[0])  # on r=1, n=x
rhs = -sp.Rational(1, 2) * ball_int(Nball * R_lin(q))
print(' dArea - 2 dVol =', dArea - 2 * dVol, '  -1/2 int N R_lin =', rhs, ' equal:', sp.simplify(dArea - 2 * dVol - rhs) == 0)

# ---- l=2 numbers
print('==== l=2 specialisation')
Es = sp.Matrix(3, 3, lambda i, j: sp.Symbol('E%d%d' % (min(i, j), max(i, j))))
Es[2, 2] = -Es[0, 0] - Es[1, 1]
psi = sp.expand((xv.T * Es * xv)[0])
I = ball_int(psi**2)
EE = sp.expand(sum(Es[i, j]**2 for i in range(3) for j in range(3)))
print(' int_B psi^2 / (E_ij E_ij) =', sp.simplify(I / EE), ' (paper: 8 pi/105)')
c2 = sp.Rational(4) * sp.factorial(0) / sp.factorial(4)
print(' S coefficient: c_2*(3/2)*(8pi/105)/8 =', c2 * sp.Rational(3, 2) * 8 * sp.pi / 105 / 8, ' (paper: pi/420)')
print(' area coefficient on B^2: c_2*(4/2)*(8pi/105)/8 =', c2 * 2 * 8 * sp.pi / 105 / 8, ' = (4/3)*pi/420 =', sp.Rational(4, 3) * sp.pi / 420)

# ---- eq (1) vs the electric table
print('==== eq (1) evaluated for U = r^l P_l: 8G S/eps^2')
for l in range(2, 7):
    cl = sp.Rational(4) * sp.factorial(l - 2) / sp.factorial(l + 2)
    intU2 = 4 * sp.pi / ((2 * l + 1) * (2 * l + 3))
    val = cl * sp.Rational(l + 1, 2) * (l * (l - 1))**2 * intU2
    print(' l=%d: %s ; in units (c_l/16G) int psi^2: %s' % (l, sp.nsimplify(val), l + 1))
