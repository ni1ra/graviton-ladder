"""Lane SUBSPACE, check C2: the multiplier in Hislop-Longo 1982, Theorem 2 (CMP 84, p. 76) is that of a
dimension-1 field, and the flow is the diamond's conformal Killing flow.
  v_lambda:  z_i(lambda) = [(1+z_i) - e^lambda (1-z_i)] / [(1+z_i) + e^lambda (1-z_i)],  z_1 = x0 + x3, z_2 = x0 - x3
  f_lambda(x) = gamma(x; lambda) f(v_{-lambda}(x)),
  gamma = 2^6 [(1+z1) + e^{-lambda}(1-z1)]^{-3} [(1-z2) + e^{lambda}(1+z2)]^{-3}
  Delta^{-i lambda} = Gamma[V(2 pi lambda)],  V(lambda) phi[f] Omega = phi[f_lambda] Omega.
Usage: python hl_gamma.py"""
import sympy as sp

lam, z1, z2 = sp.symbols("lambda z1 z2", real=True)
t, z = (z1 + z2) / 2, (z1 - z2) / 2          # x0 and x3
fails = []


def check(name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + name + ("  " + str(info) if info != "" else ""), flush=True)
    if not ok:
        fails.append(name)


def zmap(zz, l):
    return ((1 + zz) - sp.exp(l) * (1 - zz)) / ((1 + zz) + sp.exp(l) * (1 - zz))


gamma = 64 * ((1 + z1) + sp.exp(-lam) * (1 - z1)) ** -3 * ((1 - z2) + sp.exp(lam) * (1 + z2)) ** -3
check("gamma(lambda = 0) = 1", sp.simplify(gamma.subs(lam, 0) - 1) == 0)
dg = sp.simplify(sp.diff(gamma, lam).subs(lam, 0))
check("d gamma / d lambda at 0 = -3 x0", sp.simplify(dg + 3 * t) == 0, dg)

# generator of v_lambda in the (t, z) plane and the diamond field zeta = (1/2)(1 - t^2 - z^2) d_t - t z d_z (on x=y=0)
dz1 = sp.simplify(sp.diff(zmap(z1, lam), lam).subs(lam, 0))
dz2 = sp.simplify(sp.diff(zmap(z2, lam), lam).subs(lam, 0))
vt, vz = sp.simplify((dz1 + dz2) / 2), sp.simplify((dz1 - dz2) / 2)
zeta_t, zeta_z = (1 - t**2 - z**2) / 2, -t * z
check("generator of v_lambda = -zeta on the x0-x3 plane", sp.simplify(vt + zeta_t) == 0 and sp.simplify(vz + zeta_z) == 0,
      (vt, vz))
# group property of the point maps (one-parameter group)
mu = sp.symbols("mu", real=True)
check("v_lambda is a one-parameter group", sp.simplify(zmap(zmap(z1, lam), mu) - zmap(z1, lam + mu)) == 0)
# the Jacobian of v_{-lambda} in the 2-plane, and gamma against |det|^{(4 - Delta)/4} in 4d:
# a conformal map with scale factor Omega has 4d Jacobian Omega^4; here Omega^2 = (dz1'/dz1)(dz2'/dz2).
j1 = sp.diff(zmap(z1, -lam), z1)
j2 = sp.diff(zmap(z2, -lam), z2)
Om2 = sp.simplify(j1 * j2)
for Delta in (1, 2):
    # both sides are positive on the double cone; compare squares so that only rational functions appear
    # (first run compared the fractional power directly and sympy did not reduce a true zero: see logs)
    ok = sp.simplify(sp.together(gamma ** 2 - Om2 ** (4 - Delta))) == 0
    num = [sp.N((gamma - Om2 ** sp.Rational(4 - Delta, 2)).subs({lam: sp.Rational(3, 10), z1: sp.Rational(23, 63), z2: sp.Rational(-5, 63)}))]
    if Delta == 1:
        check("gamma = (conformal factor of v_{-lambda})^{4 - Delta} with Delta = 1 (test function of a dimension-1 field)", ok, num)
    else:
        print(("CAUGHT " if not ok else "MISSED ") + "plant: Delta = 2 does not reproduce gamma  " + str(num), flush=True)
        if ok:
            fails.append("plant Delta=2")
# induced infinitesimal action on the field: phi[f_lambda] = int phi f_lambda, f_lambda ~ f + lambda (dg f + zeta.grad f);
# int phi zeta.grad f = -int f (zeta.grad phi + (div zeta) phi), so
#   delta phi = -lambda (zeta.grad phi + ((div zeta) - dg) phi);  div zeta = -4 t in 4d.
# (first run had the sign of dg wrong in this scripted formula: (div zeta) + dg; kept in logs/hl_gamma_first_run.log)
div_zeta = -4 * t
check("induced field generator = zeta.d phi + (1/4)(d.zeta) phi (dimension 1), up to the sign of lambda",
      sp.simplify((div_zeta - dg) - div_zeta / 4) == 0, f"(d.zeta) - dgamma = {sp.simplify(div_zeta - dg)}, (d.zeta)/4 = {sp.simplify(div_zeta/4)}")
for Delta in (2, 3):
    # a field of dimension Delta would need (div zeta) - dg = Delta (div zeta)/4
    ok = sp.simplify((div_zeta - dg) - Delta * div_zeta / 4) == 0
    print(("CAUGHT " if not ok else "MISSED ") + f"plant: dimension {Delta} generator", flush=True)
    if ok:
        fails.append(f"plant generator Delta={Delta}")
print("DONE hl_gamma: fails =", fails)
