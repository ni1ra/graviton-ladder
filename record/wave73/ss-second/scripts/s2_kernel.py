"""P1, D = 5/6, C1 (off-coincidence), C2, C3(b)(c). Exact rational arithmetic.
Graviton:  <Psi(1) Psi(2)> = 8 pi G  T^{ab}T^{cd} d^4 W_D,   Maxwell: <phi(1) phi(2)> = - T^{ab} d_a d_b W_D,
W_D = Gamma(D/2-1)/(4 pi^{D/2}) sigma^{-(D-2)/2},  x = u1 l1 - u2 l2 - i eps t.
Printed quantities are the kernels with W_D replaced by sigma^{-p}, p = (D-2)/2, times sigma^p (pure rationals).
"""
import sympy as sp
from common import direction, conj, dot, Ttensor, apply_T, power_derivs

R = sp.Rational

def kernels(D, par1, par2, u1, u2, eps, same_helicity=False):
    l1, m1 = direction(par1)
    l2, m2 = direction(par2)
    m2p = m2 if same_helicity else conj(m2)
    xs = sp.symbols('x0:%d' % D)
    T = Ttensor(l1, m1, l2, m2p)
    f2 = apply_T(T, {0: sp.Integer(1)}, xs)
    f4 = apply_T(T, f2, xs)
    t = [1] + [0]*(D - 1)
    xv = [u1*l1[a] - u2*l2[a] - sp.I*eps*t[a] for a in range(D)]
    sub = dict(zip(xs, xv))
    sigma = sp.expand(dot(xv, xv))
    p = R(D - 2, 2)
    g = power_derivs(p, sigma, 4)
    grav = sp.nsimplify(sp.expand(sum(P.subs(sub)*g[k] for k, P in f4.items())))
    maxw = sp.nsimplify(sp.expand(-sum(P.subs(sub)*g[k] for k, P in f2.items())))
    c = dot(m1, m2p)
    return sp.simplify(grav), sp.simplify(maxw), sigma, c

pts4 = [((R(1, 3), R(-1, 5)), (R(2, 7), R(3, 4)), R(2, 5), R(7, 10)),
        ((R(-2, 3), R(1, 2)), (R(1, 9), R(-4, 5)), R(3, 10), R(1, 4)),
        ((R(1, 7), R(1, 11)), (R(1, 6), R(1, 13)), R(9, 10), R(1, 2))]
ok = True
print('== D = 4, eps = 0, non-coincident generators (P1, C2) ==')
for par1, par2, u1, u2 in pts4:
    gk, mk, sig, c = kernels(4, par1, par2, u1, u2, 0)
    print(' graviton kernel:', gk, '  Maxwell kernel:', mk, '  sigma =', sig)
    ok &= (gk == 0 and mk == 0)
print('== D = 4, eps > 0: closed forms  grav = 384 eps^4 c^2/sigma^4 [x sigma^-1],  maxw = 8 eps^2 c/sigma^2 [x sigma^-1] ==')
for (par1, par2, u1, u2), eps in zip(pts4, [R(1, 10), R(1, 3), R(2, 1)]):
    gk, mk, sig, c = kernels(4, par1, par2, u1, u2, eps)
    dg = sp.simplify(gk - 384*eps**4*c**2/sig**4)
    dm = sp.simplify(mk - 8*eps**2*c/sig**2)
    print(' eps =', eps, ' grav diff:', dg, ' maxw diff:', dm)
    ok &= (dg == 0 and dm == 0)
print('== C3(c): same helicity <Psi_0 Psi_0>, D = 4, eps = 1/3 ==')
gk, mk, sig, c = kernels(4, pts4[0][0], pts4[0][1], pts4[0][2], pts4[0][3], R(1, 3), same_helicity=True)
print(' graviton:', gk, ' Maxwell:', mk)
ok &= (gk == 0 and mk == 0)

print('== D = 6 and D = 5, eps = 0, non-coincident generators ==')
pts6 = [((R(1, 3), R(-1, 5), R(1, 4), R(2, 3)), (R(2, 7), R(3, 4), R(-1, 2), R(1, 5)), R(2, 5), R(7, 10)),
        ((R(1, 3), R(-1, 5), R(1, 4), R(2, 3)), (R(1, 3), R(-1, 5), R(1, 4) + R(1, 20), R(2, 3)), R(2, 5), R(7, 10))]
nz6 = True
for par1, par2, u1, u2 in pts6:
    gk, mk, sig, c = kernels(6, par1, par2, u1, u2, 0)
    print(' D=6 graviton:', sp.N(gk, 8), ' Maxwell:', sp.N(mk, 8), ' sigma =', sp.N(sig, 6))
    nz6 &= (gk != 0 and mk != 0)
pts5 = [((R(1, 3), R(-1, 5), R(1, 4)), (R(2, 7), R(3, 4), R(-1, 2)), R(2, 5), R(7, 10))]
nz5 = True
for par1, par2, u1, u2 in pts5:
    gk, mk, sig, c = kernels(5, par1, par2, u1, u2, 0)
    print(' D=5 graviton:', sp.N(gk, 8), ' Maxwell:', sp.N(mk, 8), ' sigma =', sp.N(sig, 6))
    nz5 &= (gk != 0 and mk != 0)
print(' D=6 non-ultralocal (graviton and Maxwell):', nz6, '  D=5:', nz5)
print('== C3(b) plant: D = 5 kernel fed to the D = 4 ultralocality test ==')
print(' flagged as non-zero:', nz5)

print('== C1: scalar off-coincidence, d_u1 d_u2 [(u1 u2)^p sigma^-p] at eps = 0, D = 4, 5, 6 ==')
a1, a2, th = sp.symbols('u1 u2 theta', positive=True)
for D in (4, 5, 6):
    p = R(D - 2, 2)
    expr = (a1*a2)**p*(2*a1*a2*(1 - sp.cos(th)))**(-p)
    val = sp.simplify(sp.diff(expr, a1, a2))
    print(' D =', D, ':', val)
    ok &= (val == 0)
print('S2 PASS' if (ok and nz6 and nz5) else 'S2 FAIL')
