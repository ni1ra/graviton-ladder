"""Blind check 5: is the 'cut velocity' of section 3(c) unique?

Find ALL polynomial metric perturbations h on the ball with
  (i)  Ric_lin[h] = (l+1) K inside the ball (the flowed electric curvature, 2 pi eps stripped),
  (ii) theta(1-r) h satisfies the linearised Hamiltonian constraint as a distribution on R^3:
       int_B h_ij (d_i d_j f - delta_ij Lap f) = 0 for all test polynomials f,
and evaluate the pairing  int_B h.K .  Paper: h = (l+1)/(l+2) * 2 N K, pairing = 2 (l+1)/(l+2) * A.
"""
import sympy as sp
import random

x, y, z = sp.symbols('x y z', real=True)
X = [x, y, z]
r2 = x**2 + y**2 + z**2
N = (1 - r2) / 2
LC = sp.LeviCivita
random.seed(11)

def monos(deg):
    return [x**a * y**b * z**(deg - a - b) for a in range(deg + 1) for b in range(deg + 1 - a)]

def mono_sphere(a, b, c):
    if a % 2 or b % 2 or c % 2:
        return 0
    n = a + b + c
    return 2 * sp.gamma(sp.Rational(a + 1, 2)) * sp.gamma(sp.Rational(b + 1, 2)) * sp.gamma(sp.Rational(c + 1, 2)) / sp.gamma(sp.Rational(n + 3, 2))

def ball_int(f):
    f = sp.Poly(sp.expand(f), x, y, z)
    return sp.expand(sum(co * mono_sphere(*m) / (sum(m) + 3) for m, co in f.terms()))

def lap(f):
    return sum(sp.diff(f, v, 2) for v in X)

def ric_lin(q):
    tr = q.trace()
    return sp.Matrix(3, 3, lambda i, j: sp.Rational(1, 2) * (
        sum(sp.diff(q[j, k], X[k], X[i]) + sp.diff(q[i, k], X[k], X[j]) for k in range(3))
        - lap(q[i, j]) - sp.diff(tr, X[i], X[j])))

def run(l, phi):
    assert sp.expand(lap(phi)) == 0
    g = [sp.diff(phi, v) for v in X]
    v = [sum(LC(i, k, m) * X[k] * g[m] for k in range(3) for m in range(3)) for i in range(3)]
    K = sp.Matrix(3, 3, lambda i, j: sp.Rational(1, 2) * (sp.diff(v[j], X[i]) + sp.diff(v[i], X[j])))
    unknowns = []
    idx = [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]
    comp = {}
    for (i, j) in idx:
        e = 0
        for d in range(l + 1, -1, -2):
            for m in monos(d):
                s = sp.Symbol('h_%d%d_%d' % (i, j, len(unknowns)))
                unknowns.append(s)
                e += s * m
        comp[(i, j)] = e
    H = sp.Matrix(3, 3, lambda i, j: comp[(min(i, j), max(i, j))])
    if CONTINUOUS:
        # finite one-particle norm: h must be continuous across the sphere, i.e. vanish on r=1
        H = (1 - r2) * H
    eqs = []
    R = ric_lin(H) - (l + 1) * K
    for i in range(3):
        for j in range(i, 3):
            e = sp.expand(R[i, j])
            if e != 0:
                eqs.extend(sp.Poly(e, x, y, z).coeffs())
    # distributional constraint, tested on all monomials f up to degree l+5
    for d in range(0, l + 6):
        for f in monos(d):
            Hf = sp.Matrix(3, 3, lambda i, j: sp.diff(f, X[i], X[j]))
            integrand = sum(H[i, j] * Hf[i, j] for i in range(3) for j in range(3)) - H.trace() * lap(f)
            e = ball_int(integrand)
            if e != 0:
                eqs.append(e)
    A_, b_ = sp.linear_eq_to_matrix(eqs, unknowns)
    sol = sp.linsolve((A_, b_), unknowns)
    if sol == sp.EmptySet:
        print(' l=%d: no solution in ansatz' % l)
        return
    sol = list(sol)[0]
    Hs = H.subs(dict(zip(unknowns, sol))).applyfunc(sp.expand)
    free = sorted(Hs.free_symbols - {x, y, z}, key=str)
    pairing = sp.simplify(ball_int(sum(Hs[i, j] * K[i, j] for i in range(3) for j in range(3))))
    A = ball_int(N * sum(K[i, j]**2 for i in range(3) for j in range(3)))
    print(' l=%d phi=%s' % (l, phi))
    print('   free (gauge) parameters left:', len(free), '  pairing depends on them:', bool(pairing.free_symbols))
    print('   pairing / A =', sp.simplify(pairing / A), '   paper: 2(l+1)/(l+2) =', sp.Rational(2 * (l + 1), l + 2))
    # is the paper's profile in the solution family?
    target = sp.Rational(l + 1, l + 2) * 2 * N * K
    diff = (Hs - target).applyfunc(sp.expand)
    eqd = []
    for i in range(3):
        for j in range(i, 3):
            if diff[i, j] != 0:
                eqd.extend(sp.Poly(diff[i, j], x, y, z).coeffs())
    ok = True
    if eqd:
        s2 = sp.linsolve(eqd, free) if free else sp.EmptySet
        ok = (s2 != sp.EmptySet)
    print('   paper profile (l+1)/(l+2)*2NK is a member of the family:', ok)

for CONTINUOUS in (False, True):
    print('=== require h = 0 on the sphere (continuity across the sphere):', CONTINUOUS)
    run(2, x * y + 2 * y * z - z**2 + x**2)
    run(3, 2 * z**3 - 3 * z * (x**2 + y**2) + x * y * z)
