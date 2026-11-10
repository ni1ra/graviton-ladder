"""s2: the area side on the scope family (THEOREM.md, Part B). Exact.

  B1  Green identity: int_B N R_lin[q] = -2 dA[q] + 4 dV[q] for every symmetric q
  B2  d(A - kV) is invariant under q -> q + sym grad xi with xi free on S iff k = 2
  B3  second-order Hamiltonian constraint on g = eta + 2 eps t K + eps^2 h2: R_lin[q2] = |K|^2
      (scalar curvature of delta + eps^2 q is eps^2 R_lin[q] + O(eps^4): checked on polynomial q)
  B4  explicit q2 = u delta with -2 Lap u = |K|^2: -2 d(A - 2V)[q2] = A; unchanged by gauge, harmonic
      conformal and source-free TT additions, while -2 dA alone moves
  B5  no first-order drive of the corner or the slice: tr K = 0 and K_nn = 0 on S (so the bilinear term
      eps * (tau, s) in the area and volume vanishes), and the induced metric has no O(eps) term
  B6  the quadratic form of the displacement is nondegenerate on degrees j >= 2
Plants: k = 2.2; the mixed reading Z1L1 = -6 dA + 14 dV; a wrong lapse; a source scaled by 11/10.
"""
import sympy as sp
from lib import *

R = Report('s2_area_side')
xv = xvec()
half = sp.Rational(1, 2)


def dA(q):
    return half * sphere_int(sp.expand(q.trace() - (xv.T * q * xv)[0]))   # n = x on S


def dV(q):
    return half * ball_int(sp.expand(q.trace()))


def symgrad2(xi):
    return mat(lambda i, j: sp.diff(xi[j], X[i]) + sp.diff(xi[i], X[j]))


def rand_sym(seed):
    import random
    rnd = random.Random(seed)
    mons = [x**a * y**b * z**c for a in range(4) for b in range(4) for c in range(4) if a + b + c <= 4]
    def p():
        return sum(rnd.randint(-3, 3) * m for m in rnd.sample(mons, 6))
    M = sp.zeros(3, 3)
    for i in range(3):
        for j in range(i, 3):
            M[i, j] = M[j, i] = p()
    return M


def rand_vec(seed):
    import random
    rnd = random.Random(seed)
    mons = [x**a * y**b * z**c for a in range(4) for b in range(4) for c in range(4) if a + b + c <= 3]
    return sp.Matrix([sum(rnd.randint(-3, 3) * m for m in rnd.sample(mons, 5)) for _ in range(3)])


# B1
for seed in (1, 2, 3):
    q = rand_sym(seed)
    lhs = ball_int(sp.expand(N * rlin(q)))
    R.check(f'B1 Green identity, random q #{seed}', sp.simplify(lhs - (-2 * dA(q) + 4 * dV(q))) == 0)
    Nw = N + (1 - r2)**2
    R.plant(f'P lapse N + (1-r^2)^2 breaks the Green identity, q #{seed}',
            sp.simplify(ball_int(sp.expand(Nw * rlin(q))) - (-2 * dA(q) + 4 * dV(q))) != 0)

# B2
for seed in (11, 12, 13):
    xi = rand_vec(seed)
    g = symgrad2(xi)
    R.check(f'B2 d(A-2V) kills sym grad xi, xi free on S, #{seed}', sp.simplify(dA(g) - 2 * dV(g)) == 0)
    R.check(f'B2 dA alone is not invariant (or xi_n has zero mean), #{seed}',
            sp.simplify(dA(g) - 2 * sphere_int(sp.expand((xv.T * xi)[0]))) == 0)
dil = symgrad2(sp.Rational(1, 10) * xv)        # a 10 per cent dilation, pure gauge
R.check('B2 dilation: dA = 8 pi /10, dV = 4 pi /10', sp.simplify(dA(dil) - sp.pi * 8 / 10) == 0 and sp.simplify(dV(dil) - sp.pi * 4 / 10) == 0)
R.plant('P k = 2.2: d(A - 2.2V) of a pure-gauge dilation is nonzero', sp.simplify(dA(dil) - sp.Rational(22, 10) * dV(dil)) != 0)
z1l1 = sp.simplify(-6 * dA(dil) + 14 * dV(dil))
R.plant('P Z1L1 (fixed zeta on a moved corner) = -6 dA + 14 dV: nonzero on a pure-gauge dilation', z1l1 != 0, z1l1)
R.value('Z1L1 on a 10 per cent dilation', z1l1)

# B3: scalar curvature of delta + e q to first order in e is e R_lin[q]
e = sp.symbols('e')
q = rand_sym(5)
g3 = sp.eye(3) + e * q
# first-order inverse and Christoffels are enough for the O(e) scalar curvature
ginv = sp.eye(3) - e * q
def chris(a, b, c):
    return sum(ginv[a, d] * (sp.diff(g3[d, b], X[c]) + sp.diff(g3[d, c], X[b]) - sp.diff(g3[b, c], X[d])) for d in range(3)) / 2
Gam = [[[chris(a, b, c) for c in range(3)] for b in range(3)] for a in range(3)]
def ricci(b, c):
    t = sum(sp.diff(Gam[a][b][c], X[a]) - sp.diff(Gam[a][b][a], X[c]) for a in range(3))
    t += sum(Gam[a][a][d] * Gam[d][b][c] - Gam[a][c][d] * Gam[d][b][a] for a in range(3) for d in range(3))
    return t
Rs = sum(ginv[b, c] * ricci(b, c) for b in range(3) for c in range(3))
R1 = sp.expand(sp.diff(Rs, e).subs(e, 0))
R.check('B3 scalar curvature of delta + e q is e R_lin[q] + O(e^2)', sp.expand(R1 - rlin(q)) == 0 and sp.expand(Rs.subs(e, 0)) == 0)


def inv_lap(p):
    """A polynomial u with Lap u = p (p polynomial): term by term on homogeneous harmonic pieces is overkill;
    solve the linear system on monomials of degree deg p + 2."""
    p = sp.expand(p)
    deg = sp.Poly(p, x, y, z).total_degree()
    mons = [x**a * y**b * z**c for a in range(deg + 3) for b in range(deg + 3) for c in range(deg + 3) if a + b + c <= deg + 2]
    cs = sp.symbols(f'c0:{len(mons)}')
    u = sum(ci * m for ci, m in zip(cs, mons))
    eqs = sp.Poly(sp.expand(lap(u) - p), x, y, z).coeffs()
    sol = sp.solve(eqs, cs, dict=True)[0]
    return sp.expand(u.subs(sol).subs({ci: 0 for ci in cs}))


for l in (2, 3, 4):
    for name in ('m=l cos', 'm=0', 'mixed'):
        phi = harmonics(l)[name]
        tag = f'l={l} {name}: '
        v, K = magnetic(phi)
        K2 = sp.expand(contract(K, K))
        A = sp.simplify(ball_int(N * K2))
        # B5
        R.check(tag + 'B5 tr K = 0 and x.K.x = 0: no eps-linear drive of corner time shift or slice', sp.expand(K.trace()) == 0 and sp.expand((xv.T * K * xv)[0]) == 0)
        # B4
        u = inv_lap(-K2 / 2)
        q2 = u * sp.eye(3)
        R.check(tag + 'B3 R_lin[q2] = |K|^2 for q2 = u delta', sp.expand(rlin(q2) - K2) == 0)
        tot = sp.simplify(-2 * (dA(q2) - 2 * dV(q2)))
        R.check(tag + 'B4 -2 d(A-2V)[q2] = A', sp.simplify(tot - A) == 0, (tot, A))
        w = harmonics(l + 1)['mixed'] + harmonics(2)['m=0']            # harmonic conformal term
        _, Ktt = magnetic(harmonics(l + 1)['m=l cos'])                  # a source-free TT term
        xi = rand_vec(7 + l) + sp.Rational(2, 7) * (1 + x**2) * xv   # a gauge vector with nonzero flux through S
        R.check(tag + 'gauge probe has nonzero flux through S', sphere_int(sp.expand((xv.T * xi)[0])) != 0)
        q2b = q2 + symgrad2(xi) + w * sp.eye(3) + 3 * Ktt
        R.check(tag + 'B4 other second-order gauges still solve the constraint', sp.expand(rlin(q2b) - K2) == 0)
        totb = sp.simplify(-2 * (dA(q2b) - 2 * dV(q2b)))
        R.check(tag + 'B4 total unchanged by gauge + harmonic conformal + TT additions', sp.simplify(totb - A) == 0)
        R.check(tag + 'B4 -2 dA alone does move (the split is gauge dependent)', sp.simplify(dA(q2b) - dA(q2)) != 0)
        R.plant(tag + 'P source x 11/10: total is 1.1 A', sp.simplify(-2 * (dA(q2 * sp.Rational(11, 10)) - 2 * dV(q2 * sp.Rational(11, 10))) - A) != 0)
        R.plant(tag + 'P k = 2.2 on the on-shell q2b differs from q2', sp.simplify((dA(q2b) - sp.Rational(22, 10) * dV(q2b)) - (dA(q2) - sp.Rational(22, 10) * dV(q2))) != 0)
        zz = sp.simplify((-6 * dA(q2b) + 14 * dV(q2b)) - (-6 * dA(q2) + 14 * dV(q2)))
        R.plant(tag + 'P Z1L1 total is gauge dependent on shell', zz != 0, zz)
        R.value(tag + 'A', A)

# B6: the quadratic form of the displacement, derived here from the flat area and volume elements
a, sg, gs, gt = sp.symbols('a sigma gsq gtq')          # gsq = |grad_S sigma|^2, gtq = |grad_S tau|^2
area_sigma = sp.series((1 + a * sg)**2 * sp.sqrt(1 + a**2 * gs / (1 + a * sg)**2), a, 0, 3).removeO()
vol_sigma = sp.series((1 + a * sg)**3 / 3, a, 0, 3).removeO()
quad_sigma = sp.expand(area_sigma - 2 * vol_sigma).coeff(a, 2)
R.check('B6 radial shift: second-order density of A - 2V is |grad sigma|^2/2 - sigma^2', sp.expand(quad_sigma - (gs / 2 - sg**2)) == 0)
R.check('B6 radial shift: first-order density of A - 2V vanishes (k = 2 is the mean curvature)', sp.expand(area_sigma - 2 * vol_sigma).coeff(a, 1) == 0)
area_tau = sp.series(sp.sqrt(1 - a**2 * gt), a, 0, 3).removeO()
R.check('B6 corner time shift: second-order area density is -|grad tau|^2/2', area_tau.coeff(a, 2) == -gt / 2)
for jj in (2, 3, 4):
    Yj = harmonics(jj)['mixed']
    en = ball_int(sp.expand(sum(sp.diff(Yj, v)**2 for v in X)))
    R.check(f'B6 harmonic extension energy int_B |grad(r^j Y)|^2 = j oint Y^2, j={jj}', sp.simplify(en - jj * sphere_int(sp.expand(Yj**2))) == 0)
# so on degree j: sigma-part (j(j+1) - 2)/2 = (j-1)(j+2)/2, tau-part -j(j+1)/2 + j = -j(j-1)/2
jv = sp.symbols('j')
R.check('B6 Hessians (j-1)(j+2)/2 and -j(j-1)/2: zero only at j = 1 (sigma) and j = 0, 1 (tau): 3 + 4 = 7 Poincare modes',
        sp.factor(jv * (jv + 1) / 2 - 1) == (jv - 1) * (jv + 2) / 2 and sp.factor(-jv * (jv + 1) / 2 + jv) == -jv * (jv - 1) / 2)
R.write()
