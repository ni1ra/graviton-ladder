"""S1 (Green identity) and S2 (corner-gauge invariance of delta(A - 2V)), exact.
Usage: python3 s1s2_green_gauge.py OUT.json"""
import json, random, sys
import sympy as sp
from common import x, y, z, X, R_lin, dA, dV, symgrad, rand_poly, ball_int

rng = random.Random(64)
r2 = x**2 + y**2 + z**2
N = (1 - r2)/2
out = {"S1": [], "S1_general_w": [], "S1_plant_k22": [], "S2": [], "S2_plants": []}

for trial in range(3):
    q = [[None]*3 for _ in range(3)]
    for i in range(3):
        for j in range(i, 3):
            q[i][j] = q[j][i] = rand_poly(rng, 4)
    lhs = ball_int(N*R_lin(q))
    a, v = dA(q), dV(q)
    out["S1"].append({"trial": trial, "int_N_Rlin": str(lhs), "minus2_dA_plus4_dV": str(-2*a + 4*v),
                      "pass": bool(sp.simplify(lhs - (-2*a + 4*v)) == 0)})
    # plant: k = 2.2 in place of 2
    plant = -2*a + 2*sp.Rational(22, 10)*v
    out["S1_plant_k22"].append({"trial": trial, "caught": bool(sp.simplify(lhs - plant) != 0)})
    # general weight with w|_S = 0, d_r w|_S = -kappa (kappa = 1): w = N + c (1 - r^2)^2
    w = N + sp.Rational(3, 7)*(1 - r2)**2
    bulk = sum(q[i][j]*(sp.diff(w, X[i], X[j]) - (1 if i == j else 0)*sum(sp.diff(w, X[k], 2) for k in range(3)))
               for i in range(3) for j in range(3))
    lw = ball_int(w*R_lin(q)); rw = ball_int(bulk) - 2*a
    out["S1_general_w"].append({"trial": trial, "pass": bool(sp.simplify(lw - rw) == 0)})

# S2: corner gauge, xi not vanishing on S
for trial in range(4):
    xi = [rand_poly(rng, 3) for _ in range(3)]
    q = symgrad(xi)
    a, v = dA(q), dV(q)
    out["S2"].append({"trial": trial, "dA": str(a), "dV": str(v), "d(A-2V)": str(sp.simplify(a - 2*v)),
                      "pass_invariant": bool(sp.simplify(a - 2*v) == 0),
                      "dA_alone_moves": bool(a != 0), "dV_alone_moves": bool(v != 0),
                      "Rlin_zero": bool(sp.simplify(R_lin(q)) == 0)})
    # plant: k = 2.2
    out["S2_plants"].append({"trial": trial, "k=2.2_caught": bool(sp.simplify(a - sp.Rational(22, 10)*v) != 0)})

# radius-R control: on B_R the invariant combination is A - (2/R) V; k = 2 at R = 2 must fail
xi = [rand_poly(rng, 3) for _ in range(3)]
q = symgrad(xi)
R = sp.Integer(2)
aR, vR = dA(q, R), dV(q, R)
out["S2_plants"].append({"radius_2": True, "A-1*V_invariant": bool(sp.simplify(aR - vR) == 0),
                         "A-2V_at_R2_caught": bool(sp.simplify(aR - 2*vR) != 0)})
# dilation (the l = 0 corner displacement) and an l = 2 corner displacement, explicitly
for name, xi in [("dilation", [x, y, z]), ("grad_r2Y20", [sp.diff((2*z**2 - x**2 - y**2), s) for s in X])]:
    q = symgrad(xi)
    a, v = dA(q), dV(q)
    out["S2"].append({"explicit": name, "dA": str(a), "dV": str(v), "pass_invariant": bool(sp.simplify(a - 2*v) == 0)})

out["ALL_PASS"] = (all(t["pass"] for t in out["S1"]) and all(t["caught"] for t in out["S1_plant_k22"])
                   and all(t["pass"] for t in out["S1_general_w"])
                   and all(t["pass_invariant"] for t in out["S2"])
                   and all(t.get("Rlin_zero", True) for t in out["S2"])
                   and all(t.get("k=2.2_caught", True) and t.get("A-2V_at_R2_caught", True) and t.get("A-1*V_invariant", True)
                           for t in out["S2_plants"]))
s = json.dumps(out, indent=1); print(s)
open(sys.argv[1], "w").write(s + "\n")
