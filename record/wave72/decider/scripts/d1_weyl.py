"""w72-decider. Electric area side from the exact static Weyl solution, exact rational arithmetic.

F = Area - 2*Volume in the slice t = 0 of
  ds^2 = -e^{2U}dt^2 + e^{-2U}[e^{2k}(drho^2+dz^2) + rho^2 dphi^2],  U = eps*r^l*P_l(cos theta),
on the surface r = 1 + eps*a*P_l, expanded to eps^2, stationary in a.
Reported: X_l = area side / (c_l/8G * 1/2 * int_B psi^2), psi = l(l-1)U. Modular side is l+1 in these units.
Run: python3 d1_weyl.py   (prints PASS/FAIL lines; exit code 1 on any FAIL)
"""
import sys, random, itertools
import sympy as sp
from sympy import Rational as Q, pi, legendre, factorial

eps, a, r, x = sp.symbols('epsilon a r x', real=True)
rho, z = sp.symbols('rho z', real=True)
FAILS = []


def check(name, ok, detail=''):
    print(('PASS ' if ok else 'FAIL ') + name + (('  ' + str(detail)) if detail else ''))
    if not ok:
        FAILS.append(name)


def harmonic_u(l):
    """r^l P_l(z/r) as a polynomial in rho, z."""
    t = sp.symbols('t')
    P = sp.Poly(legendre(l, t), t)
    u = 0
    for (k,), c in P.terms():
        u += c * z**k * (rho**2 + z**2)**((l - k) // 2)
    return sp.expand(u)


def weyl_k(u):
    s = sp.symbols('s', positive=True)
    integrand = (s * (sp.diff(u, rho)**2 - sp.diff(u, z)**2)).subs(rho, s)
    kap = sp.expand(sp.integrate(sp.expand(integrand), (s, 0, rho)))
    assert sp.expand(sp.diff(kap, z) - 2 * rho * sp.diff(u, rho) * sp.diff(u, z)) == 0, 'k not integrable'
    return kap


def to_rx(expr):
    return sp.expand(expr.subs({rho: r * sp.sqrt(1 - x**2), z: r * x}))


def trunc(expr, n=2):
    expr = sp.expand(expr)
    return sum(expr.coeff(eps, i) * eps**i for i in range(n + 1))


def F_series(l, plant=None):
    """Coefficients (F0, F1, F2(a)) of Area - 2 Volume."""
    u = harmonic_u(l)
    kap = weyl_k(u)
    if plant == 'nok':
        kap = sp.Integer(0)
    U = to_rx(u)          # r^l P_l(x)
    Kp = to_rx(kap)       # r^{2l} * poly(x)
    w = Q(1, 2) if plant == 'conf' else Q(1)      # conformal factor e^{-2wU}
    # area factor e^{-2wU+k}, volume factor e^{-3wU+2k}, to eps^2
    areaf = 1 - 2 * w * eps * U + eps**2 * (2 * w**2 * U**2 + Kp)
    volf = 1 - 3 * w * eps * U + eps**2 * (Q(9, 2) * w**2 * U**2 + 2 * Kp)
    Pl = legendre(l, x)
    aa = 0 if plant == 'noshift' else a
    Rs = 1 + eps * aa * Pl
    dR2 = eps**2 * aa**2 * (1 - x**2) * sp.diff(Pl, x)**2      # (dR/dtheta)^2
    area_int = trunc(areaf.subs(r, Rs) * Rs * (Rs + dR2 / 2))
    vol_inner = sp.integrate(sp.expand(volf * r**2), (r, 0, r))  # antiderivative in r, zero at 0
    vol_int = trunc(vol_inner.subs(r, Rs))
    F = 2 * pi * sp.integrate(sp.expand(area_int - 2 * vol_int), (x, -1, 1))
    F = sp.expand(F)
    return [sp.expand(F.coeff(eps, i)) for i in range(3)]


def stationary(F2):
    al, be, ga = F2.coeff(a, 2), F2.coeff(a, 1), F2.coeff(a, 0)
    if al == 0:
        return None, be, ga, ga
    astar = -be / (2 * al)
    return astar, be, ga, sp.nsimplify(sp.expand(ga - be**2 / (4 * al)))


def c_l(l):
    return Q(4) * factorial(l - 2) / factorial(l + 2)


def X_of(l, plant=None):
    F0, F1, F2 = F_series(l, plant)
    astar, be, ga, F2s = stationary(sp.expand(F2))
    intu2 = 2 * pi * Q(1, 2 * l + 3) * Q(2, 2 * l + 1)
    X = sp.nsimplify(-4 * F2s / (c_l(l) * l**2 * (l - 1)**2 * intu2))
    return F0, F1, F2s, astar, X


# ---------- C1: first-law identity on the flat unit ball ----------
X3 = sp.symbols('x1 x2 x3', real=True)


def sph_mono(e):
    if any(k % 2 for k in e):
        return sp.Integer(0)
    g = sp.gamma
    return 2 * g(Q(e[0] + 1, 2)) * g(Q(e[1] + 1, 2)) * g(Q(e[2] + 1, 2)) / g(Q(sum(e) + 3, 2))


def integrate_poly(p, kind):
    """kind: 'S' sphere, 'B' ball, 'NB' ball with weight N=(1-r^2)/2."""
    p = sp.Poly(sp.expand(p), *X3)
    tot = 0
    for e, c in p.terms():
        n = sum(e)
        s = sph_mono(e)
        if kind == 'S':
            tot += c * s
        elif kind == 'B':
            tot += c * s / (n + 3)
        else:
            tot += c * s * Q(1, 2) * (Q(1, n + 3) - Q(1, n + 5))
    return sp.simplify(tot)


def c1_identity(seed, deg=4):
    rnd = random.Random(seed)
    monos = [X3[0]**i * X3[1]**j * X3[2]**k for i in range(deg + 1) for j in range(deg + 1 - i)
             for k in range(deg + 1 - i - j)]
    h = sp.zeros(3, 3)
    for i in range(3):
        for j in range(i, 3):
            h[i, j] = h[j, i] = sum(rnd.randint(-3, 3) * m for m in rnd.sample(monos, 6))
    tr = h.trace()
    r2 = sum(v**2 for v in X3)
    hnn_r2 = sum(h[i, j] * X3[i] * X3[j] for i in range(3) for j in range(3))  # n.h.n on r=1
    dA = Q(1, 2) * integrate_poly(tr - hnn_r2, 'S')
    dV = Q(1, 2) * integrate_poly(tr, 'B')
    dR = sum(sp.diff(h[i, j], X3[i], X3[j]) for i in range(3) for j in range(3)) - sum(sp.diff(tr, v, 2) for v in X3)
    rhs = -Q(1, 2) * integrate_poly(dR, 'NB')
    return sp.simplify(dA - 2 * dV - rhs) == 0, (dA - 2 * dV)


# ---------- C2: normalisation of the scope family ----------
def c2_norm(l, phi):
    xs = sp.Matrix(X3)
    grad = sp.Matrix([sp.diff(phi, v) for v in X3])
    v = xs.cross(grad)
    K = sp.Matrix(3, 3, lambda i, j: Q(1, 2) * (sp.diff(v[i], X3[j]) + sp.diff(v[j], X3[i])))
    curl = sp.Matrix(3, 3, lambda i, j: sum(sp.LeviCivita(i, k, m) * sp.diff(K[m, j], X3[k])
                                            for k in range(3) for m in range(3)))
    psi = sp.expand(sum(curl[i, j] * X3[i] * X3[j] for i in range(3) for j in range(3)))
    ratio = sp.simplify(psi / phi)
    K2 = sp.expand(sum(K[i, j]**2 for i in range(3) for j in range(3)))
    A = integrate_poly(K2, 'NB')
    E2 = Q(l + 2, 2) * integrate_poly(psi**2, 'B')
    return ratio, sp.simplify(A - c_l(l) * E2) == 0


def harmonic_xyz(l, kind):
    if kind == 'zonal':
        return sp.expand(harmonic_u(l).subs({rho: sp.sqrt(X3[0]**2 + X3[1]**2), z: X3[2]}))
    return sp.expand(sp.re(sp.expand((X3[0] + sp.I * X3[1])**(l - 1))) * X3[2])  # sectoral-ish, harmonic


# ---------- C3: scalar curvature of the Weyl slice through eps^2 ----------
def c3_scalar_curvature(l):
    u = harmonic_u(l)
    kap = weyl_k(u)
    ph = sp.symbols('phi_', real=True)
    co = [rho, z, ph]
    U_, k_ = eps * u, eps**2 * kap
    g = sp.diag(sp.exp(-2 * U_ + 2 * k_), sp.exp(-2 * U_ + 2 * k_), rho**2 * sp.exp(-2 * U_))
    gi = g.inv()
    Gam = [[[sum(gi[i, m] * (sp.diff(g[m, j], co[k]) + sp.diff(g[m, k], co[j]) - sp.diff(g[j, k], co[m]))
                 for m in range(3)) / 2 for k in range(3)] for j in range(3)] for i in range(3)]
    Rs = 0
    for j in range(3):
        for k in range(3):
            ric = sum(sp.diff(Gam[i][j][k], co[i]) for i in range(3)) - sum(sp.diff(Gam[i][j][i], co[k]) for i in range(3)) \
                + sum(Gam[i][i][m] * Gam[m][j][k] for i in range(3) for m in range(3)) \
                - sum(Gam[i][k][m] * Gam[m][j][i] for i in range(3) for m in range(3))
            Rs += gi[j, k] * ric
    ser = sp.series(sp.simplify(Rs), eps, 0, 3).removeO()
    return sp.simplify(sp.expand(ser)) == 0


if __name__ == '__main__':
    print('== controls ==')
    for seed in (1, 2, 3):
        ok, val = c1_identity(seed)
        check(f'C1 first-law identity, seed {seed}', ok, f'd(A-2V) = {val}')
    for l in (2, 3, 4):
        for kind in ('zonal', 'tess'):
            phi = harmonic_xyz(l, kind)
            assert sp.expand(sum(sp.diff(phi, v, 2) for v in X3)) == 0
            ratio, okA = c2_norm(l, phi)
            want = Q(1, 2) * (l - 1) * l * (l + 1)
            check(f'C2 psi/phi and A = c_l E_2, l={l} {kind}', (abs(ratio) == want) and okA, f'psi/phi = {ratio}')
    for l in (2, 3):
        check(f'C3 Weyl slice R = 0 through eps^2, l={l}', c3_scalar_curvature(l))
    F0, F1, F2 = F_series(1)
    astar, be, ga, F2s = stationary(sp.expand(F2))
    check('C4 l=1 (Rindler) second-order change is 0', F1 == 0 and be == 0 and F2s == 0, f'F2* = {F2s}, beta = {be}')

    print('== electric area side ==')
    rows = []
    for l in (2, 3, 4, 5, 6):
        F0, F1, F2s, astar, X = X_of(l)
        check(f'l={l}: F0 = 4pi/3 and F1 = 0', sp.simplify(F0 - 4 * pi / 3) == 0 and F1 == 0)
        bin_ = {l + 1: 'MOD', l + 2: 'CAN', l: 'SCAL'}.get(X, 'OTHER')
        rows.append((l, sp.nsimplify(-2 * F2s), astar, X, bin_))
        print(f'  l={l}: area side (1/8G) = {sp.nsimplify(-2*F2s)}   a* = {astar}   X_l = {X}   bin = {bin_}')

    print('== C5 plants (each must change X_l) ==')
    for l in (2, 3):
        base = X_of(l)[4]
        for plant in ('nok', 'noshift', 'conf'):
            Xp = X_of(l, plant)[4]
            check(f'C5 plant {plant}, l={l} changes X', Xp != base, f'X = {Xp} (base {base})')

    print('== summary ==')
    bins = {b for *_, b in rows}
    print('bins:', bins, '| controls failed:', FAILS)
    sys.exit(1 if FAILS else 0)
