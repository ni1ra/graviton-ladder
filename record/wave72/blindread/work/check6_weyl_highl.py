"""Blind check 6: electric area side beyond the paper's table, via the decomposition
F2(a) = alpha a^2 - beta a + gamma  with
  alpha = 2 pi (l-1)(l+2)/(2l+1)                  (flat second variation)
  beta  = 8 pi (l-1)/(2l+1)
  gamma = 4 pi (4l-3)/((2l+1)(2l+3)) + oint kappa - 4 int_B kappa      (coordinate sphere, a=0)
where k = eps^2 kappa is the Weyl function.  Checks the conjectured closed form
  oint kappa = -4 pi l /((2l-1)(2l+1))   and hence  -2 F2* = 8 pi l (l-1)/((l+2)(2l+1)(2l+3)) = eq.(1).
"""
import sympy as sp
rho, zz, r, c = sp.symbols('rho z r c', real=True)

def kappa_on_sphere_integral(l):
    R = sp.sqrt(rho**2 + zz**2)
    u = sp.expand(sp.simplify(R**l * sp.legendre(l, zz / R)))
    ur, uz = sp.diff(u, rho), sp.diff(u, zz)
    krho = sp.expand(rho * (ur**2 - uz**2))
    kz = sp.expand(2 * rho * ur * uz)
    assert sp.expand(sp.diff(krho, zz) - sp.diff(kz, rho)) == 0
    k = sp.expand(sp.integrate(krho, (rho, 0, rho)))
    assert sp.expand(sp.diff(k, zz) - kz) == 0
    # on unit sphere: rho^2 = 1-c^2, z = c ; k is even in rho
    ks = sp.expand(k.subs(rho, sp.sqrt(1 - c**2)).subs(zz, c))
    return 2 * sp.pi * sp.integrate(ks, (c, -1, 1))

for l in range(1, 13):
    ok = kappa_on_sphere_integral(l)
    conj = -4 * sp.pi * l / ((2 * l - 1) * (2 * l + 1))
    gamma = 4 * sp.pi * (4 * l - 3) / ((2 * l + 1) * (2 * l + 3)) + ok * (1 - sp.Rational(4, 2 * l + 3))
    alpha = 2 * sp.pi * (l - 1) * (l + 2) / (2 * l + 1)
    beta = 8 * sp.pi * (l - 1) / (2 * l + 1)
    F2s = gamma - beta**2 / (4 * alpha) if l > 1 else gamma
    pred = 8 * sp.pi * l * (l - 1) / ((l + 2) * (2 * l + 1) * (2 * l + 3))
    print('l=%2d  oint kappa = %s  (conjecture %s, match %s)   -2F2* = %s   eq.(1): %s   equal: %s'
          % (l, ok, conj, sp.simplify(ok - conj) == 0, sp.nsimplify(-2 * F2s), pred, sp.simplify(-2 * F2s - pred) == 0))
