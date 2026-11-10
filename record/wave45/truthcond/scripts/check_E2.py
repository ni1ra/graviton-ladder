#!/usr/bin/env python3
"""w45-truthcond: exact check of the clauses the truth condition quotes from Paper A r7 (Sec. 2, 3, 5).

Stationary magnetic family at d = 4: phi = r^l Y (here Re/Im (x+iy)^l and a mixed harmonic),
v = x cross grad phi, K = sym grad v (TT, harmonic), h = 2tK so E = 0 and B = -sym curl K at t = 0.
Checks, exact in sympy, l = 2..6:
  (1) delta E_mod = 2 sym(x cross B) = (l+1) K                    [eq. (3) at t = 0]
  (2) Ric3(2 N K) = (l+2) K,  N = (1 - r^2)/2                      [Prop. 2 proof]
  (3) so delta E_geom - delta E_mod = K (one unit of the datum's momentum), ratio (l+1)/(l+2)
  (4) A/E_chain = (l+2)/(l+1) from Sec. 3 closed forms (chain target and A only, no gluing)
Planted: N -> (1 - r^2)/2 * 1.1 must break (2); B -> 1.1 B must break (1).
Run: python3 check_E2.py > ../out/check_E2.json
"""
import json
import sympy as sp

x, y, z = X = sp.symbols('x y z', real=True)
r2 = x**2 + y**2 + z**2
eps = sp.LeviCivita

def grad(f): return [sp.diff(f, c) for c in X]
def sym(T): return [[sp.Rational(1, 2) * (T[i][j] + T[j][i]) for j in range(3)] for i in range(3)]
def cross_x(T):   # (x x T)_ij = eps_ikl x_k T_lj
    return [[sum(eps(i, k, l) * X[k] * T[l][j] for k in range(3) for l in range(3)) for j in range(3)] for i in range(3)]
def curl(T):      # (curl T)_ij = eps_ikl d_k T_lj
    return [[sum(eps(i, k, l) * sp.diff(T[l][j], X[k]) for k in range(3) for l in range(3)) for j in range(3)] for i in range(3)]
def ric3(h):      # linearised Ricci of a spatial perturbation
    tr = sum(h[k][k] for k in range(3))
    lap = lambda f: sum(sp.diff(f, c, 2) for c in X)
    return [[sp.Rational(1, 2) * (sum(sp.diff(h[k][j], X[k], X[i]) + sp.diff(h[k][i], X[k], X[j]) for k in range(3))
                                   - lap(h[i][j]) - sp.diff(tr, X[i], X[j])) for j in range(3)] for i in range(3)]
def ratio(A, Bm):
    """Return c if A == c*Bm exactly (c constant), else None."""
    c = None
    for i in range(3):
        for j in range(3):
            b = sp.expand(Bm[i][j]); a = sp.expand(A[i][j])
            if b == 0:
                if a != 0: return None
                continue
            q = sp.simplify(a / b)
            if q.free_symbols: return None
            if c is None: c = q
            elif sp.simplify(q - c) != 0: return None
    return c

def harmonics(l):
    w = sp.expand((x + sp.I * y)**l)
    re, im = sp.re(w), sp.im(w)
    mixed = sp.expand(re + 3 * im)
    if l >= 2:   # add an m = 0 piece too: r^l P_l(z/r)
        P = sp.legendre(l, sp.Symbol('u'))
        rl0 = sp.expand(sp.simplify((sp.sqrt(r2))**l * P.subs(sp.Symbol('u'), z / sp.sqrt(r2))))
        mixed = sp.expand(mixed + rl0)
    return {"Re(x+iy)^l": re, "Im(x+iy)^l": im, "mixed incl. m=0": mixed}

out = {"cases": [], "plants": {}}
N = (1 - r2) / 2
for l in range(2, 7):
    for name, phi in harmonics(l).items():
        assert sp.expand(sum(sp.diff(phi, c, 2) for c in X)) == 0
        g = grad(phi)
        v = [sp.expand(X[(i + 1) % 3] * g[(i + 2) % 3] - X[(i + 2) % 3] * g[(i + 1) % 3]) for i in range(3)]
        K = sym([[sp.diff(v[j], X[i]) for j in range(3)] for i in range(3)])
        tt = all(sp.expand(sum(sp.diff(K[i][j], X[i]) for i in range(3))) == 0 for j in range(3)) and sp.expand(sum(K[i][i] for i in range(3))) == 0
        Bm = [[-e for e in row] for row in sym(curl(K))]
        dEmod = [[2 * e for e in row] for row in sym(cross_x(Bm))]
        c_mod = ratio(dEmod, K)
        c_geom = ratio(ric3([[2 * N * K[i][j] for j in range(3)] for i in range(3)]), K)
        rec = {"l": l, "harmonic": name, "K_TT": tt, "dEmod_over_K": str(c_mod), "Ric3_2NK_over_K": str(c_geom),
               "difference": str(sp.simplify(c_geom - c_mod)) if None not in (c_mod, c_geom) else None,
               "pass": c_mod == l + 1 and c_geom == l + 2}
        out["cases"].append(rec)
    # (4) closed forms of Sec. 3 (per N_Y)
    Q = sp.Rational((l + 1)**3 * l**2 * (l - 1)**2, 8 * (2 * l + 3))
    A = sp.Rational(l * (l + 1) * (l - 1), 2 * (2 * l + 3))
    c_l = sp.Rational(4, (l - 1) * l * (l + 1) * (l + 2))
    out.setdefault("chain_ratio", {})[l] = {"A_over_Echain": str(A / (c_l * Q)), "pass": A / (c_l * Q) == sp.Rational(l + 2, l + 1)}

# plants at l = 2
phi = sp.expand(sp.re((x + sp.I * y)**2)); g = grad(phi)
v = [sp.expand(X[(i + 1) % 3] * g[(i + 2) % 3] - X[(i + 2) % 3] * g[(i + 1) % 3]) for i in range(3)]
K = sym([[sp.diff(v[j], X[i]) for j in range(3)] for i in range(3)])
cg = ratio(ric3([[2 * sp.Rational(11, 10) * N * K[i][j] for j in range(3)] for i in range(3)]), K)
Bm = [[-sp.Rational(11, 10) * e for e in row] for row in sym(curl(K))]
cm = ratio([[2 * e for e in row] for row in sym(cross_x(Bm))], K)
out["plants"] = {"N_x1.1_Ric3_over_K": str(cg), "caught_N": cg != 4, "B_x1.1_dEmod_over_K": str(cm), "caught_B": cm != 3}
out["ALL_PASS"] = all(c["pass"] for c in out["cases"]) and all(v["pass"] for v in out["chain_ratio"].values()) \
    and out["plants"]["caught_N"] and out["plants"]["caught_B"]
print(json.dumps(out, indent=1, default=str))
