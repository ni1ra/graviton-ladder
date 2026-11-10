"""Blind check 3: electric entropy by the paper's DIRECT route (section 3 a-c), not via duality.

Electric solution k: h_ij = -2 eps U delta_ij, momentum 0.  Modular velocity on the slice (eq. 4, rate 2 pi):
  delta E = 0,   delta B = -2 sym(x x E) = -2 eps K[U]   (times 2 pi).
Cut velocity: data (delta h = 0, delta pi) supported in the closed ball, obeying the momentum constraint
distributionally on R^3 (d_j pi_ij = 0 inside, n_j pi_ij = 0 on the sphere), whose magnetic curvature inside is
-2 K[U].  Here pi = p - delta tr p,  B_ij = -eps_ikl d_k p_lj  (same convention as the paper's B = -curl K for p = K).
Pairing (paper's normalisation extended to non-TT data: h.pi with pi = p - delta tr p):
  S = (2 pi /(32 pi G)) | int h_k . delta pi | = (eps^2/8G) |T|,   T = int_B U tr(delta pi)   [target normalised to -2K[U]].
Paper's eq. (1) predicts |T| = 2 l (l-1)/(l+2) int_B U^2.
"""
import sympy as sp
import sys

x, y, z = sp.symbols('x y z', real=True)
X = [x, y, z]
r2 = x**2 + y**2 + z**2
LC = sp.LeviCivita

def monos(deg):
    return [x**a * y**b * z**(deg - a - b) for a in range(deg + 1) for b in range(deg + 1 - a)]

def mono_sphere(a, b, c):
    if a % 2 or b % 2 or c % 2:
        return 0
    n = a + b + c
    return 2 * sp.gamma(sp.Rational(a + 1, 2)) * sp.gamma(sp.Rational(b + 1, 2)) * sp.gamma(sp.Rational(c + 1, 2)) / sp.gamma(sp.Rational(n + 3, 2))

def ball_int(f, gens):
    f = sp.Poly(sp.expand(f), *gens)
    return sp.expand(sum(co * mono_sphere(*m) / (sum(m) + 3) for m, co in f.terms()))

def run(l, U):
    assert sp.expand(sum(sp.diff(U, v, 2) for v in X)) == 0
    gU = [sp.diff(U, v) for v in X]
    vU = [sum(LC(i, k, m) * X[k] * gU[m] for k in range(3) for m in range(3)) for i in range(3)]
    KU = sp.Matrix(3, 3, lambda i, j: sp.Rational(1, 2) * (sp.diff(vU[j], X[i]) + sp.diff(vU[i], X[j])))
    degs_pi = [d for d in range(l + 2, -1, -2)]
    degs_q = [d for d in range(l + 1, -1, -2)]
    unknowns = []
    idx = [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]
    comp = {}
    for (i, j) in idx:
        e = 0
        for d in degs_pi:
            for m in monos(d):
                s = sp.Symbol('p_%d%d_%d' % (i, j, len(unknowns)))
                unknowns.append(s)
                e += s * m
        comp[(i, j)] = e
    PI = sp.Matrix(3, 3, lambda i, j: comp[(min(i, j), max(i, j))])
    q = []
    for i in range(3):
        e = 0
        for d in degs_q:
            for m in monos(d):
                s = sp.Symbol('q_%d_%d' % (i, len(unknowns)))
                unknowns.append(s)
                e += s * m
        q.append(e)
    P = PI - sp.Rational(1, 2) * PI.trace() * sp.eye(3)        # p = pi - (1/2) delta tr pi
    eqs = []
    def add(expr):
        expr = sp.expand(expr)
        if expr != 0:
            eqs.extend(sp.Poly(expr, x, y, z).coeffs())
    for i in range(3):
        add(sum(sp.diff(PI[i, j], X[j]) for j in range(3)))                         # momentum constraint
        add(sum(X[j] * PI[i, j] for j in range(3)) - (1 - r2) * q[i])               # n.pi = 0 on r=1
    for i in range(3):
        for j in range(3):
            Bij = -sum(LC(i, k, m) * sp.diff(P[m, j], X[k]) for k in range(3) for m in range(3))
            add(Bij + 2 * KU[i, j])                                                 # B[delta p] = -2 K[U]
    A, b = sp.linear_eq_to_matrix(eqs, unknowns)
    sol = sp.linsolve((A, b), unknowns)
    if sol == sp.EmptySet:
        print(' l=%d: NO polynomial solution with this ansatz' % l)
        return
    sol = list(sol)[0]
    subs = dict(zip(unknowns, sol))
    PIs = PI.subs(subs).applyfunc(sp.expand)
    free = sorted(PIs.free_symbols - {x, y, z}, key=str)
    T = ball_int(U * PIs.trace(), (x, y, z))
    T = sp.simplify(T)
    IU = ball_int(U**2, (x, y, z))
    print(' l=%d  U=%s' % (l, U))
    print('   free parameters left in delta pi:', len(free))
    print('   T = int_B U tr(delta pi) =', T, '  depends on free parameters:', bool(T.free_symbols))
    print('   |T| / int U^2 =', sp.simplify(sp.Abs(T) / IU) if not T.free_symbols else 'n/a',
          '   paper eq.(1) predicts 2l(l-1)/(l+2) =', sp.Rational(2 * l * (l - 1), l + 2))
    # show one representative solution with free parameters set to zero
    rep = PIs.subs({f: 0 for f in free})
    print('   tr(delta pi) (free params = 0) =', sp.factor(rep.trace()))

run(2, 2 * z**2 - x**2 - y**2)
run(2, x * y + 3 * y * z)
run(3, 2 * z**3 - 3 * z * (x**2 + y**2))
run(3, x * y * z + x**3 - 3 * x * y**2)
if len(sys.argv) > 1:
    run(4, 8 * z**4 - 24 * z**2 * (x**2 + y**2) + 3 * (x**2 + y**2)**2)
