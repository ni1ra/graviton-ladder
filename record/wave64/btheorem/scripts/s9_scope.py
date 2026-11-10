"""S9: the scope family's second-order area side, exact.  h1 = 2 t K (magnetic l = 2, 3), h1|_Sigma = 0.
Usage: python3 s9_scope.py OUT.json"""
import json, sys
import sympy as sp
from common import x, y, z, X, R_lin, dA, dV, symgrad, ball_int, sphere_int

eps = sp.Symbol('epsilon')
out = {}
r2 = x**2 + y**2 + z**2
N = (1 - r2)/2
LC = sp.LeviCivita

# (a) R[delta + eps^2 q] = eps^2 R_lin[q] + O(eps^4): exact Ricci scalar for a polynomial q
q = [[x*y + z**2, 2*x*z, y**2], [2*x*z, x**2 - y*z, x*y], [y**2, x*y, 3*z**2 + x]]
q = [[q[min(i, j)][max(i, j)] for j in range(3)] for i in range(3)]
g = sp.Matrix(3, 3, lambda i, j: (1 if i == j else 0) + eps**2*q[i][j])
ginv = sp.eye(3) - eps**2*sp.Matrix(q)      # exact to O(eps^4)
Gam = [[[sum(ginv[a, b]*(sp.diff(g[b, i], X[j]) + sp.diff(g[b, j], X[i]) - sp.diff(g[i, j], X[b])) for b in range(3))/2
         for j in range(3)] for i in range(3)] for a in range(3)]
Ric = [[sum(sp.diff(Gam[a][i][j], X[a]) - sp.diff(Gam[a][i][a], X[j]) for a in range(3))
        + sum(Gam[a][a][b]*Gam[b][i][j] - Gam[a][j][b]*Gam[b][i][a] for a in range(3) for b in range(3))
        for j in range(3)] for i in range(3)]
Rs = sp.expand(sum(ginv[i, j]*Ric[i][j] for i in range(3) for j in range(3)))
R2 = sp.expand(Rs).coeff(eps, 2)
out["(a)_R_order_eps2_equals_Rlin"] = bool(sp.simplify(R2 - R_lin(q)) == 0)
out["(a)_no_order_eps1"] = bool(sp.expand(Rs).coeff(eps, 1) == 0)

def magnetic_K(phi):
    V = [sum(LC(i, j, k)*X[j]*sp.diff(phi, X[k]) for j in range(3) for k in range(3)) for i in range(3)]
    return [[sp.expand((sp.diff(V[j], X[i]) + sp.diff(V[i], X[j]))/2) for j in range(3)] for i in range(3)]

for l, phi in [(2, x**2 - y**2), (3, sp.expand(sp.re((x + sp.I*y)**3)))]:
    K = magnetic_K(phi)
    rec = {}
    rec["traceless"] = sp.simplify(sum(K[i][i] for i in range(3))) == 0
    rec["transverse"] = all(sp.simplify(sum(sp.diff(K[i][j], X[i]) for i in range(3))) == 0 for j in range(3))
    rec["harmonic"] = all(sp.simplify(sum(sp.diff(K[i][j], s, 2) for s in X)) == 0 for i in range(3) for j in range(3))
    K2 = sp.expand(sum(K[i][j]**2 for i in range(3) for j in range(3)))
    # (b) Hamiltonian constraint at O(eps^2) with gamma = delta + eps^2 q2, K_ij = eps K: R_lin[q2] = |K|^2 (tr K = 0)
    A = ball_int(N*K2); flux = sphere_int(K2)
    rec["A=int N|K|^2"] = str(A); rec["oint|K|^2"] = str(flux)
    rec["A_over_flux"] = str(sp.simplify(A/flux))
    rec["A_over_flux_equals_1/((2l+1)(2l+3))"] = bool(sp.simplify(A/flux - sp.Rational(1, (2*l + 1)*(2*l + 3))) == 0)
    # explicit q2 = u delta with -2 lap u = |K|^2 (polynomial particular solution)
    deg = sp.Poly(K2, *X).total_degree() + 2
    cs = []; u = 0
    for i in range(deg + 1):
        for j in range(deg + 1 - i):
            for k in range(deg + 1 - i - j):
                cc = sp.Symbol(f'c_{i}_{j}_{k}'); cs.append(cc); u += cc*x**i*y**j*z**k
    eqs = sp.Poly(sp.expand(-2*sum(sp.diff(u, s, 2) for s in X) - K2), *X).coeffs()
    sol = sp.solve(eqs, cs, dict=True)[0]
    u = sp.expand(u.subs(sol).subs({cc: 0 for cc in cs}))
    q2 = [[u if i == j else 0 for j in range(3)] for i in range(3)]
    rec["Rlin_q2_equals_K2"] = bool(sp.simplify(R_lin(q2) - K2) == 0)
    area_side = sp.simplify(-2*(dA(q2) - 2*dV(q2)))
    rec["-2 d(A-2V)[q2]"] = str(area_side)
    rec["area_side_equals_A"] = bool(sp.simplify(area_side - A) == 0)
    # homogeneous and gauge additions leave it unchanged; dA alone moves (plant)
    uh = x**2 - y**2 + 3*x*z + z      # harmonic: R_lin[uh delta] = 0
    q2b = [[q2[i][j] + (uh if i == j else 0) for j in range(3)] for i in range(3)]
    gq = symgrad([x + y**2, z*x, y - z**3])
    q2c = [[q2b[i][j] + gq[i][j] for j in range(3)] for i in range(3)]
    rec["Rlin_unchanged"] = bool(sp.simplify(R_lin(q2c) - K2) == 0)
    rec["area_side_unchanged_by_homog+gauge"] = bool(sp.simplify(-2*(dA(q2c) - 2*dV(q2c)) - A) == 0)
    rec["plant_dA_alone_moves"] = bool(sp.simplify(-2*dA(q2c) - (-2*dA(q2))) != 0)
    rec = {k: (bool(v) if isinstance(v, (bool, sp.logic.boolalg.BooleanAtom)) else v) for k, v in rec.items()}
    out[f"l={l}"] = rec

out["paper_A_unit_flux_gate_l2_equals_1/35"] = out["l=2"]["A_over_flux"] == "1/35"
out["ALL_PASS"] = (out["(a)_R_order_eps2_equals_Rlin"] and out["(a)_no_order_eps1"] and out["paper_A_unit_flux_gate_l2_equals_1/35"]
                   and all(v for l in (2, 3) for k, v in out[f"l={l}"].items() if isinstance(v, bool) and not k.startswith("plant"))
                   and all(out[f"l={l}"]["plant_dA_alone_moves"] for l in (2, 3)))
s = json.dumps(out, indent=1); print(s)
open(sys.argv[1], "w").write(s + "\n")
