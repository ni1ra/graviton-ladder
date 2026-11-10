"""Blind check 2: electric area side from the Weyl metric, written from the paper's description.
ds^2 = -e^{2U}dt^2 + e^{-2U}[e^{2k}(drho^2+dz^2) + rho^2 dphi^2],  U = eps r^l P_l(cos th).
F = Area - 2 Vol on r = 1 + eps a P_l, expanded to eps^2, stationary in a.
"""
import sympy as sp

rho, zz, r, c, eps, a = sp.symbols('rho z r c epsilon a', real=True)   # c = cos(theta)

def solid(l):
    # r^l P_l(cos th) as polynomial in rho, z
    R = sp.sqrt(rho**2 + zz**2)
    u = sp.expand(sp.simplify(R**l * sp.legendre(l, zz / R)))
    return sp.simplify(u)

def kappa(u):
    # k = eps^2 kappa, kappa_rho = rho(u_rho^2 - u_z^2), kappa_z = 2 rho u_rho u_z, kappa=0 on axis
    ur, uz = sp.diff(u, rho), sp.diff(u, zz)
    krho = sp.expand(rho * (ur**2 - uz**2))
    kz = sp.expand(2 * rho * ur * uz)
    # integrability
    assert sp.expand(sp.diff(krho, zz) - sp.diff(kz, rho)) == 0
    k = sp.integrate(krho, (rho, 0, rho))   # k(rho=0,z)=0 on the axis (regular axis)
    assert sp.expand(sp.diff(k, zz) - kz) == 0
    return sp.expand(k)

def run(l):
    u = solid(l)
    assert sp.expand(sp.diff(u, rho, 2) + sp.diff(u, rho) / rho + sp.diff(u, zz, 2)) == 0
    k2 = kappa(u)
    s = sp.sqrt(1 - c**2)
    sub = {rho: r * s, zz: r * c}
    u_rc = sp.expand(sp.simplify(u.subs(sub)))        # = r^l P_l(c)
    k_rc = sp.expand(sp.simplify(k2.subs(sub)))
    Pl = sp.legendre(l, c)
    assert sp.expand(u_rc - r**l * Pl) == 0
    # Volume: integrand e^{-3U+2k} r^2 dr dc dphi
    volint = sp.series(sp.exp(-3 * eps * u_rc + 2 * eps**2 * k_rc), eps, 0, 3).removeO() * r**2
    Rs = 1 + eps * a * Pl
    inner = sp.integrate(sp.expand(volint), (r, 0, r))           # polynomial in r
    inner = inner.subs(r, Rs)
    inner = sp.series(sp.expand(inner), eps, 0, 3).removeO()
    Vol = 2 * sp.pi * sp.integrate(sp.expand(inner), (c, -1, 1))
    # Area: 2 pi int dth  [rho e^{-U}] * [e^{-U+k} sqrt(R^2 + R_th^2)]
    #  with dth = -dc/s, R_th = dR/dth = -s dR/dc ; rho = R s
    dRdc = sp.diff(Rs, c)
    integrand = Rs * sp.exp(-2 * eps * u_rc.subs(r, Rs) + eps**2 * k_rc.subs(r, Rs)) * sp.sqrt(Rs**2 + (1 - c**2) * dRdc**2)
    integrand = sp.series(integrand, eps, 0, 3).removeO()
    Area = 2 * sp.pi * sp.integrate(sp.expand(integrand), (c, -1, 1))
    F = sp.expand(Area - 2 * Vol)
    F0 = F.coeff(eps, 0); F1 = F.coeff(eps, 1); F2 = sp.expand(F.coeff(eps, 2))
    dF = sp.expand(sp.diff(F2, a))
    if dF == 0:
        astar = ['any (F2 independent of a)']
        F2s = sp.simplify(F2)
    else:
        astar = sp.solve(dF, a)
        F2s = sp.simplify(F2.subs(a, astar[0]))
    return F0, F1, F2, astar, F2s

for l in range(1, 7):
    F0, F1, F2, astar, F2s = run(l)
    cl = sp.Rational(4) * sp.factorial(max(l - 2, 0)) / sp.factorial(l + 2)
    pred = (cl * sp.Rational(l + 1, 2) * (l * (l - 1))**2 * 4 * sp.pi / ((2 * l + 1) * (2 * l + 3))) if l >= 2 else 0
    print('l=%d: F0=%s  F1=%s  F2(a)=%s' % (l, F0, sp.simplify(F1), sp.factor(F2)))
    print('      a*=%s (paper 2/(l+2)=%s)   area side in units eps^2/8G = -2 F2* = %s   eq.(1) gives %s   equal: %s'
          % (astar, sp.Rational(2, l + 2), sp.nsimplify(-2 * F2s), pred, sp.simplify(-2 * F2s - pred) == 0))
    if l >= 2:
        unit = cl / 2 * (l * (l - 1))**2 * 4 * sp.pi / ((2 * l + 1) * (2 * l + 3))  # (c_l/16G) int psi^2 in units eps^2/8G
        print('      in units of the section-1 table:', sp.simplify(-2 * F2s / unit))
