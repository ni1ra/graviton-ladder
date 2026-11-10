"""S11: exact l = 0 classification closure with arbitrary profiles (sympy Functions, no ansatz).
h = a(r) delta + b(r) nn, K = c(r) delta + d(r) nn.  Gauge: h -> h + (d_i xi_j + d_j xi_i), K -> K + Hess f.
Usage: python3 s11_l0.py OUT.json"""
import json, sys
import sympy as sp

x, y, z = X = sp.symbols('x y z', real=True)
r = sp.Symbol('r', positive=True)
rr = sp.sqrt(x**2 + y**2 + z**2)
n = [s/rr for s in X]
a, b, c, dd, beta, gam, R = [sp.Function(s) for s in ('a', 'b', 'c', 'd', 'beta', 'gamma', 'Rho')]
out = {}

def Dr(F, k=1):
    return sp.Subs(sp.Derivative(F(r), (r, k)), r, rr)

def to_r(e):
    """Replace sqrt(x^2+y^2+z^2) by r after differentiation; evaluate at a generic point on the ray (x,y,z)=(r,0,0)."""
    e = e.subs({x: r, y: 0, z: 0})
    return sp.simplify(e.doit())

def rad(F):  # F(r) composed with |x|
    return F(rr)

def tensor_rad(A, B):  # A(r) delta + B(r) nn
    return [[(rad(A) if i == j else 0) + rad(B)*n[i]*n[j] for j in range(3)] for i in range(3)]

def gauge_h(xi):
    return [[sp.diff(xi[j], X[i]) + sp.diff(xi[i], X[j]) for j in range(3)] for i in range(3)]

def eq_tensor(T1, T2):
    # compare all components at a generic point (x,y,z) = (r s1, r s2, r s3) on two directions (rotation covariance)
    ok = True
    for pt in [{x: r, y: 0, z: 0}, {x: r*sp.Rational(3, 5), y: r*sp.Rational(4, 5), z: 0},
               {x: r*sp.Rational(2, 7), y: r*sp.Rational(3, 7), z: r*sp.Rational(6, 7)}]:
        for i in range(3):
            for j in range(3):
                dlt = sp.simplify((T1[i][j] - T2[i][j]).subs(pt).doit())
                ok = ok and (dlt == 0)
    return ok

# (1) spatial gauge xi_i = -(1/2) a(r) x_i  removes a and maps b -> b - r a'
xi = [-sp.Rational(1, 2)*rad(a)*s for s in X]
lhs = gauge_h(xi)
target = [[(-rad(a) if i == j else 0) - rr*Dr(a)*n[i]*n[j] for j in range(3)] for i in range(3)]
out["spatial_gauge_removes_a"] = eq_tensor(lhs, target)

# (2) lapse gauge f with f'(r) = -r c(r): Hess f = -c delta - r c' nn.  Use f = F(r) with F' = -r c.
def hess_direct():
    # Hess F(|x|) = (F'/r) delta + (F'' - F'/r) nn with F' = -r c, F'' = -(c + r c')
    Fp = -rr*c(rr); Fpp = -(c(rr) + rr*Dr(c))
    return [[(Fp/rr if i == j else 0) + (Fpp - Fp/rr)*n[i]*n[j] for j in range(3)] for i in range(3)]
Hd = hess_direct()
targetK = [[(-rad(c) if i == j else 0) - rr*Dr(c)*n[i]*n[j] for j in range(3)] for i in range(3)]
out["lapse_gauge_removes_c"] = eq_tensor(Hd, targetK)
# cross-check the Hess formula itself on a concrete F
Ft = rr**4*sp.exp(-rr**2)
Hc = [[sp.diff(Ft, X[i], X[j]) for j in range(3)] for i in range(3)]
Fp = sp.diff(sp.Symbol('s')**4*sp.exp(-sp.Symbol('s')**2), sp.Symbol('s')).subs(sp.Symbol('s'), rr)
Fpp = sp.diff(sp.Symbol('s')**4*sp.exp(-sp.Symbol('s')**2), sp.Symbol('s'), 2).subs(sp.Symbol('s'), rr)
Hf = [[(Fp/rr if i == j else 0) + (Fpp - Fp/rr)*n[i]*n[j] for j in range(3)] for i in range(3)]
out["hess_formula_check"] = eq_tensor(Hc, Hf)

# (3) invariants: under the general gauge, beta = b - r a', gamma = d - r c' are invariant (symbolic xi(r), f(r))
s_ = sp.Function('s'); fl = sp.Function('l')
# general spatial gauge xi(r) n: delta a = 2 xi/r, delta b = 2 xi' - 2 xi/r
da = 2*s_(r)/r; db = 2*sp.diff(s_(r), r) - 2*s_(r)/r
out["beta_invariant"] = bool(sp.simplify(db - r*sp.diff(da, r)) == 0)
dc = sp.diff(fl(r), r)/r; ddd = sp.diff(fl(r), r, 2) - sp.diff(fl(r), r)/r
out["gamma_invariant"] = bool(sp.simplify(ddd - r*sp.diff(dc, r)) == 0)
# and the Cartesian form of the general spatial gauge (xi(r) n) matches (da, db)
xi_g = [s_(rr)*s/rr for s in X]
tg = [[(2*s_(rr)/rr if i == j else 0) + (2*Dr(s_) - 2*s_(rr)/rr)*n[i]*n[j] for j in range(3)] for i in range(3)]
out["general_spatial_gauge_form"] = eq_tensor(gauge_h(xi_g), tg)

# (4) constraints on the reduced data: q = beta nn, p = gamma nn
q = tensor_rad(lambda s: 0, beta)
tr = sum(q[k][k] for k in range(3))
Rl = sum(sp.diff(q[i][j], X[i], X[j]) for i in range(3) for j in range(3)) - sum(sp.diff(tr, s, 2) for s in X)
Rl_r = sp.simplify(Rl.subs({x: r, y: 0, z: 0}).doit())
out["Rlin_beta"] = str(Rl_r)
out["Rlin_equals_2_over_r2_(r beta)'"] = bool(sp.simplify(Rl_r - 2/r**2*sp.diff(r*beta(r), r)) == 0)
p = tensor_rad(lambda s: 0, gam)
trp = sum(p[k][k] for k in range(3))
pi = [[p[i][j] - (trp if i == j else 0) for j in range(3)] for i in range(3)]
divpi = [sum(sp.diff(pi[i][j], X[i]) for i in range(3)) for j in range(3)]
ndiv = sp.simplify(sum(n[j]*divpi[j] for j in range(3)).subs({x: r, y: 0, z: 0}).doit())
out["n_div_pi"] = str(ndiv)
out["n_div_pi_equals_2gamma_over_r"] = bool(sp.simplify(ndiv - 2*gam(r)/r) == 0)
# tangential components of div pi vanish (l = 0)
out["div_pi_tangential_zero"] = bool(sp.simplify(divpi[1].subs({x: r, y: 0, z: 0}).doit()) == 0)
# (5) regular inversion: beta = (1/(2 r)) \int_0^r s^2 Rho(s) ds solves (2/r^2)(r beta)' = Rho
s = sp.Symbol('s', positive=True)
beta_inv = sp.Integral(s**2*R(s), (s, 0, r))/(2*r)
out["inversion"] = bool(sp.simplify(2/r**2*sp.diff(r*beta_inv, r) - R(r)) == 0)
# homogeneous solutions of (r beta)' = 0: beta = C/r, singular at 0 (so regular data are determined by Rho)
C1 = sp.Symbol('C1')
out["homogeneous_is_C_over_r"] = bool(sp.simplify(sp.diff(r*(C1/r), r)) == 0)
# (6) gauge maps are smooth: xi = -(1/2) a(r) x and f = -\int_0^r s c(s) ds are smooth for smooth even a, c (by construction)
out["ALL_PASS"] = all(v for k, v in out.items() if isinstance(v, bool))
js = json.dumps(out, indent=1); print(js)
open(sys.argv[1], "w").write(js + "\n")
