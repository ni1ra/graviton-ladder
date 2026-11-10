r"""S10: uniqueness in action.  Units 1/(4G) = 1, so 16 pi G = 4 pi and FL1 reads  D(q, 0) = 2 pi \int_B N rho.
The candidate F0 = -d(A - 2V) passes every configuration; each rival violates FL1 (or T) on an explicit one.
Usage: python3 s10_power.py OUT.json"""
import json, sys
import sympy as sp
from common import x, y, z, X, R_lin, dA, dV, symgrad, ball_int, sphere_int

LC = sp.LeviCivita
r2 = x**2 + y**2 + z**2
N = (1 - r2)/2
pi = sp.pi
out = {"configs": {}, "rivals": {}}

def poisson_u(rho):
    """polynomial u with -2 lap u = 4 pi rho, so q = u delta has R_lin[q] = 4 pi rho."""
    deg = sp.Poly(sp.expand(rho), *X).total_degree() + 2
    cs = []; u = 0
    for i in range(deg + 1):
        for j in range(deg + 1 - i):
            for k in range(deg + 1 - i - j):
                cc = sp.Symbol(f'c_{i}_{j}_{k}'); cs.append(cc); u += cc*x**i*y**j*z**k
    eqs = sp.Poly(sp.expand(-2*sum(sp.diff(u, s, 2) for s in X) - 4*pi*rho), *X).coeffs()
    sol = sp.solve(eqs, cs, dict=True)[0]
    return sp.expand(u.subs(sol).subs({cc: 0 for cc in cs}))

def diag(u):
    return [[u if i == j else 0 for j in range(3)] for i in range(3)]

def add(*qs):
    return [[sum(q[i][j] for q in qs) for j in range(3)] for i in range(3)]

configs = {}
# M1: matter g = 1 + x^2 + y z (does not vanish on S); rho = g^2/2
g1 = 1 + x**2 + y*z
rho1 = sp.expand(g1**2/2)
configs["M1"] = {"q": diag(poisson_u(rho1)), "fl1": 2*pi*ball_int(N*rho1), "rho": rho1}
# M2: matter g = z (angle-dependent profile), rho = z^2/2
rho2 = sp.expand(z**2/2)
configs["M2"] = {"q": diag(poisson_u(rho2)), "fl1": 2*pi*ball_int(N*rho2), "rho": rho2}
# M3: M1 plus a source-free electric TT datum and an arbitrary gauge term (same source, same FL1 value)
phi = x**2 - y**2
V = [sum(LC(i, j, k)*X[j]*sp.diff(phi, X[k]) for j in range(3) for k in range(3)) for i in range(3)]
Kmag = [[sp.expand((sp.diff(V[j], X[i]) + sp.diff(V[i], X[j]))/2) for j in range(3)] for i in range(3)]
curlK = [[sum(LC(i, k, m)*sp.diff(Kmag[m][j], X[k]) for k in range(3) for m in range(3)) for j in range(3)] for i in range(3)]
Ee = [[sp.expand((curlK[i][j] + curlK[j][i])/2) for j in range(3)] for i in range(3)]   # even-parity TT datum
configs["E1"] = {"q": Ee, "fl1": sp.Integer(0), "rho": 0}
configs["G1_dilation"] = {"q": symgrad([x, y, z]), "fl1": sp.Integer(0), "rho": 0}
configs["G2_l2_corner"] = {"q": symgrad([sp.diff(3*z**2 - r2, s) for s in X]), "fl1": sp.Integer(0), "rho": 0}
configs["M3"] = {"q": add(configs["M1"]["q"], Ee, symgrad([y*z**2, x - z, x*y*z])), "fl1": configs["M1"]["fl1"], "rho": rho1}
# JSS slice (fixed volume): M1 plus the dilation that sets dV = 0
s = sp.Symbol('s')
qs = add(configs["M1"]["q"], [[s*e for e in row] for row in symgrad([x, y, z])])
s0 = sp.solve(sp.Eq(dV(qs), 0), s)[0]
qs = [[sp.expand(e.subs(s, s0)) for e in row] for row in qs]
configs["M1_JSS_slice"] = {"q": qs, "fl1": configs["M1"]["fl1"], "rho": rho1}

for name, c in configs.items():
    q = c["q"]
    rl = sp.simplify(R_lin(q) - 4*pi*c["rho"])
    c["dA"] = dA(q); c["dV"] = dV(q)
    out["configs"][name] = {"source_check_Rlin=4pi*rho": bool(rl == 0), "FL1_value": str(sp.simplify(c["fl1"])),
                            "dA": str(c["dA"]), "dV": str(c["dV"])}

F0 = lambda q, c: -(c["dA"] - 2*c["dV"])
Y20 = 3*z**2 - r2
rivals = {
    "F0 = -d(A-2V) (BMSV functional)": F0,
    "R1 angle-dependent area charge F0 + oint Y20 q_AA": lambda q, c: F0(q, c) + sphere_int(Y20*(sum(q[k][k] for k in range(3)) - sum(X[i]*X[j]*q[i][j] for i in range(3) for j in range(3)))),
    "R2 corner monopole F0 + dA": lambda q, c: F0(q, c) + c["dA"],
    "R3 weight N + (1-r^2)^2": lambda q, c: ball_int((N + (1 - r2)**2)*R_lin(q))/2,
    "R4 corner jet F0 + oint R_lin|_S": lambda q, c: F0(q, c) + sphere_int(R_lin(q)),
    "R6 normalisation 11/10 F0": lambda q, c: sp.Rational(11, 10)*F0(q, c),
    "JSS functional -dA (valid on its slice only)": lambda q, c: -c["dA"],
}
for rn, F in rivals.items():
    res = {}
    for name, c in configs.items():
        val = sp.simplify(F(c["q"], c))
        res[name] = {"value": str(val), "FL1_holds": bool(sp.simplify(val - c["fl1"]) == 0)}
    out["rivals"][rn] = res
# R5: a T-odd term alpha \int_B p_kk is invisible to FL1 configurations (p = 0) and is excluded by T alone
out["rivals"]["R5 T-odd F0 + int p_kk"] = "not seen by FL1 (p = 0 on every FL1 configuration); excluded by (T): D(q,p) = D(q,-p)"

ok_F0 = all(v["FL1_holds"] for v in out["rivals"]["F0 = -d(A-2V) (BMSV functional)"].values())
caught = {rn: any(not v["FL1_holds"] for v in res.values()) for rn, res in out["rivals"].items() if isinstance(res, dict) and not rn.startswith("F0") and not rn.startswith("JSS")}
jss_on_slice = out["rivals"]["JSS functional -dA (valid on its slice only)"]["M1_JSS_slice"]["FL1_holds"]
out["F0_passes_all"] = ok_F0
out["rivals_caught"] = caught
out["JSS_functional_passes_on_its_slice"] = jss_on_slice
out["ALL_PASS"] = ok_F0 and all(caught.values()) and jss_on_slice and all(v["source_check_Rlin=4pi*rho"] for v in out["configs"].values())
js = json.dumps(out, indent=1); print(js)
open(sys.argv[1], "w").write(js + "\n")
