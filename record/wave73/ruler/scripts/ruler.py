"""Lane RULER (wave 73): is the metric-versus-extrinsic-curvature asymmetry a property of the
fixed-volume ruler, or does the cone ruler have it too?   l = 2, both parities, exact.

Delta = Q_2 - Q_3 with null-cone's forms Q_p = 2 oint int_0^1 v^2 f |sigma_p|^2, f = v(1 - v),
as a bilinear form on Cauchy data at t = 0 (TT gauge):  metric block (h = H, K = 0),
extrinsic-curvature block (h = 0, K = H), cross block.  Gram matrices are printed in the chain
basis psi = r^(2j) phi, psi_dot = r^(2j) phi (j = 0..3), divided by N_Y = oint phi^2.

Usage: python ruler.py OUT.json [--plant w3]
"""
import json
import sys
import time

import sympy as sp

import nc_core as nc
from nc_core import t, v, x, y, z, XS

plant = sys.argv[sys.argv.index("--plant") + 1] if "--plant" in sys.argv else None
L = 2
NB = 4
t0 = time.time()
res = {"script": "ruler.py", "plant": plant, "checks": [], "l": L}
ok = True


def check(name, cond, detail=""):
    global ok
    cond = bool(cond)
    ok &= cond
    res["checks"].append({"name": name, "pass": cond, "detail": str(detail)})
    print(("PASS " if cond else "FAIL ") + name + (f" | {detail}" if detail != "" else ""), flush=True)


phi = nc.harmonic(L)
N_Y = nc.surf_int(phi**2)
s = x**2 + y**2 + z**2
gphi = nc.grad(phi)
Hphi = nc.hess(phi)
W2 = v**3 * nc.F_CONE if plant == "w3" else v**2 * nc.F_CONE   # weight of Q_2
W3 = v**2 * nc.F_CONE                                            # weight of Q_3 (never planted)


# ------------------------------------------------------------------ TT tensors
def even_tt(k, deg=5):
    """Even-type TT tensor with x.T.x = s^k phi (free gauge part Hess phi set to zero)."""
    cs = sp.symbols(f"c0:{4 * (deg + 1)}")
    f = [sum(cs[i * (deg + 1) + m] * s**m for m in range(deg + 1)) for i in range(4)]
    T = sp.Matrix(3, 3, lambda i, j: f[0] * Hphi[i, j] + f[1] * (XS[i] * gphi[j] + XS[j] * gphi[i])
                  + (f[2] * phi if i == j else 0) + f[3] * phi * XS[i] * XS[j])
    eqs = [sp.expand(T.trace())]
    eqs += [sp.expand(sum(sp.diff(T[i, j], XS[i]) for i in range(3))) for j in range(3)]
    eqs += [sp.expand(sum(XS[i] * XS[j] * T[i, j] for i in range(3) for j in range(3)) - s**k * phi)]
    lin = []
    for e in eqs:
        lin += [c for c in sp.Poly(e, x, y, z).coeffs()]
    sol = sp.solve(lin, cs, dict=True)
    assert sol, "no even TT tensor found"
    sol = sol[0]
    full = {c: sol.get(c, c) for c in cs}
    free = set().union(*[sp.sympify(e).free_symbols for e in full.values()]) & set(cs)
    T = T.subs(full).subs({c: 0 for c in free}).applyfunc(sp.expand)
    return T, len(free)


V = nc.vec_potential(phi)
K0 = nc.sym_grad(V)


def odd_tt(k):
    """O_0 = sym grad V;  O_(k) = a s^k K0 + s^(k-1) (xV + Vx), a = -(2k + l + 2)/((l - 1) k)."""
    if k == 0:
        return K0
    a = -sp.Rational(2 * (k - 1) + L + 4, (L - 1) * k)
    return sp.Matrix(3, 3, lambda i, j: sp.expand(a * s**k * K0[i, j] + s ** (k - 1) * (XS[i] * V[j] + XS[j] * V[i])))


# ------------------------------------------------------------------ data
def profile(e):
    """Coefficients c_j with e = sum_j c_j s^j phi (exact, asserted)."""
    e = sp.expand(e)
    if e == 0:
        return [sp.Integer(0)] * NB
    p = sp.cancel(e / phi)
    p1 = sp.Poly(sp.expand(p.subs({y: 0, z: 0})), x)
    co = [p1.coeff_monomial(x ** (2 * j)) for j in range(NB)]
    assert sp.expand(e - sum(co[j] * s**j for j in range(NB)) * phi) == 0, "chain field is not a profile times phi"
    return co


def build(parity):
    data = []          # (name, kind, cone Weyl C, chain row of length 8, cv)
    for kind in ("metric", "extrinsic"):
        for i in range(NB):
            if parity == "even":
                H, nfree = even_tt(i + 1)
                name = f"even {kind} T_{i + 1}"
                check(f"{name}: the TT tensor with x.T.x = s^k phi is unique (the gauge mode has x.T.x != 0)", nfree == 0, nfree)
            else:
                k = i + 1 if kind == "metric" else i
                H = odd_tt(k)
                name = f"odd {kind} O_{k}"
            check(f"{name}: transverse and traceless", nc.is_tt(H))
            h = nc.tt_series(H, 0 if kind == "metric" else 1)
            cv = nc.Curv(h)
            check(f"{name}: linearised vacuum equations hold exactly", cv.ricci_zero())
            W = cv.E() if parity == "even" else cv.B()
            Wo = cv.B() if parity == "even" else cv.E()
            chain = sp.expand(sum(XS[a] * XS[b] * W[a, b] for a in range(3) for b in range(3)))
            other = sp.expand(sum(XS[a] * XS[b] * Wo[a, b] for a in range(3) for b in range(3)))
            check(f"{name}: the other parity's chain field vanishes identically", other == 0)
            row = profile(chain.subs(t, 0)) + profile(sp.diff(chain, t).subs(t, 0))
            C = cv.cone_weyl()
            data.append((name, kind, C, row, cv))
            print(f"   built {name} at {time.time() - t0:.0f} s", flush=True)
    return data


def gram(mats, weight):
    n = len(mats)
    G = sp.zeros(n, n)
    for a in range(n):
        for b in range(a, n):
            G[a, b] = G[b, a] = sp.nsimplify(2 * nc.cone_int(nc.norm2(mats[a], mats[b]), weight) / N_Y, rational=True)
    return G


def e1_form():
    """c_l E_1 on the chain basis: psi block c_l (1/2 int N grad.grad + 1/2 int psi psi), psi_dot block c_l/2 int N."""
    cl = 4 * sp.factorial(L - 2) / sp.factorial(L + 2)
    N = (1 - s) / 2
    G = sp.zeros(2 * NB, 2 * NB)
    for a in range(NB):
        for b in range(NB):
            pa, pb = s**a * phi, s**b * phi
            G[a, b] = cl * (nc.ball_int(N * sum(sp.diff(pa, c) * sp.diff(pb, c) for c in XS)) / 2
                            + nc.ball_int(pa * pb) / 2) / N_Y
            G[NB + a, NB + b] = cl * nc.ball_int(N * pa * pb) / 2 / N_Y
    return G.applyfunc(lambda e: sp.nsimplify(e, rational=True))


def definiteness(G):
    """Exact: leading principal minors (Sylvester), plus float eigenvalues for the record."""
    G = sp.Matrix(G)
    if G.is_zero_matrix:
        return "zero", [], []
    minors = [G[:k, :k].det() for k in range(1, G.shape[0] + 1)]
    lam_ = sp.Symbol("lam_")
    ev = sorted(float(sp.re(r)) for r in sp.Poly(G.charpoly(lam_).as_expr(), lam_).nroots(n=30, maxsteps=500))
    if all(m > 0 for m in minors):
        verdict = "positive definite"
    elif all((-1) ** (k + 1) * m > 0 for k, m in enumerate(minors)):
        verdict = "negative definite"
    else:
        verdict = "not definite (see eigenvalues)"
    return verdict, [str(m) for m in minors], ev


def rel_spectrum(Q2, S):
    """Eigenvalues of S^-1 Q2 on a block (the possible ratios Q_2/S), as floats from exact matrices."""
    lam_ = sp.Symbol("lam_")
    pol = sp.Poly((sp.Matrix(Q2) - lam_ * sp.Matrix(S)).det(), lam_)
    return sorted(float(sp.re(r)) for r in pol.nroots(n=30, maxsteps=500))


def mstr(G):
    return [[str(e) for e in G.row(i)] for i in range(G.shape[0])]


S_FORM = e1_form()
out = {}
for parity in ("even", "odd"):
    data = build(parity)
    mats2 = [nc.moment(d[2], 2) for d in data]
    mats3 = [nc.moment(d[2], 3) for d in data]
    M = sp.Matrix([d[3] for d in data])                       # datum -> chain coordinates
    check(f"{parity}: the eight data span the chain basis (psi, psi_dot; j = 0..3)", M.det() != 0, M.det())
    Mi = M.inv()
    G2 = (Mi * gram(mats2, W2) * Mi.T).applyfunc(sp.nsimplify)
    G3 = (Mi * gram(mats3, W3) * Mi.T).applyfunc(sp.nsimplify)
    print(f"   {parity} Gram matrices at {time.time() - t0:.0f} s", flush=True)
    # which t = 0 datum feeds which chain block
    kinds = ["metric" if any(M[a, j] != 0 for a in range(NB)) else "extrinsic" for j in (0, NB)]
    blockname = {"psi": kinds[0], "psi_dot": kinds[1]}
    for a in range(2 * NB):
        lo = all(M[a, j] == 0 for j in range(NB))
        hi = all(M[a, j] == 0 for j in range(NB, 2 * NB))
        assert lo != hi, "a datum mixes psi and psi_dot"
    check(f"{parity}: Q_3 Gram equals c_l E_1 on the chain basis (slice entropy, block diagonal)",
          (G3 - S_FORM).is_zero_matrix)
    cv0 = data[0][4]
    check(f"{parity}: Q_3 of the first datum equals nc.entropy_slice",
          sp.simplify(2 * nc.cone_int(nc.norm2(mats3[0]), W3) - nc.entropy_slice(cv0, L)) == 0)
    D = G2 - G3
    blocks = {"psi": D[:NB, :NB], "psi_dot": D[NB:, NB:], "cross": D[:NB, NB:], "full": D}
    o = {"psi_is": blockname["psi"], "psi_dot_is": blockname["psi_dot"],
         "Q2": mstr(G2), "S": mstr(G3), "Delta": {}}
    for bn, B in blocks.items():
        if bn == "cross":
            o["Delta"][bn] = {"matrix": mstr(B), "is_zero": bool(B.is_zero_matrix)}
            print(f"\n[{parity}] Delta cross block (psi x psi_dot): zero = {B.is_zero_matrix}")
            sp.pprint(B)
            continue
        verdict, minors, ev = definiteness(B)
        o["Delta"][bn] = {"matrix": mstr(B), "is_zero": bool(B.is_zero_matrix), "verdict": verdict,
                          "leading_minors": minors, "eigenvalues_float": ev}
        label = f"{bn} = {blockname[bn]} data" if bn in blockname else "full 8x8"
        print(f"\n[{parity}] Delta on {label}: {verdict}; eigenvalues {['%.6g' % e for e in ev]}")
        if bn != "full":
            sp.pprint(B)
            rs = rel_spectrum(G2[:NB, :NB] if bn == "psi" else G2[NB:, NB:],
                              G3[:NB, :NB] if bn == "psi" else G3[NB:, NB:])
            o["Delta"][bn]["Q2_over_S_spectrum_float"] = rs
            print(f"   spectrum of Q_2 relative to S on this block: {['%.6g' % e for e in rs]}")
    o["static_ratio_psi00"] = str(G2[0, 0] / G3[0, 0])
    o["first_psi_dot_ratio"] = str(G2[NB, NB] / G3[NB, NB])
    print(f"\n[{parity}] Q_2/S on psi = phi ({blockname['psi']} datum): {o['static_ratio_psi00']};"
          f"  on psi_dot = phi ({blockname['psi_dot']} datum): {o['first_psi_dot_ratio']}")
    out[parity] = o
    out[parity]["_G"] = (G2, G3)

# ------------------------------------------------------------------ static reproductions, null-cone's own data
cvE = nc.Curv(nc.metric_electric(phi))
B1 = nc.Curv(nc.tt_series(K0, 1)).B()
lam = sp.cancel(B1[0, 0] / Hphi[0, 0])
cvB = nc.Curv(nc.tt_series(K0 / lam, 1))
stat = {}
for name, cv in (("static electric tide (Newtonian gauge, K = 0)", cvE), ("stationary magnetic tide h = 2tK", cvB)):
    C = cv.cone_weyl()
    q2 = 2 * nc.cone_int(nc.norm2(nc.moment(C, 2)), W2)
    q3 = 2 * nc.cone_int(nc.norm2(nc.moment(C, 3)), W3)
    stat[name] = (q2, q3)
    r = sp.cancel(q2 / q3)
    res.setdefault("static", {})[name] = str(r)
    print(f"{name}: Q_2/S = {r}")
    if plant is None:
        check(f"{name}: Q_2/S = 16/9", r == sp.Rational(16, 9), r)
    else:
        check(f"PLANT {plant} caught: {name} no longer gives 16/9", r != sp.Rational(16, 9), r)
G2e, G3e = out["even"].pop("_G")
G2o, G3o = out["odd"].pop("_G")
# chain field of both static tides is l(l-1) phi = 2 phi, so the Gram entry [0,0] times 4 N_Y must reproduce them
q2, q3 = stat["static electric tide (Newtonian gauge, K = 0)"]
check("even Gram [psi0, psi0] reproduces the static electric tide (gauge invariance; pure metric datum)",
      sp.simplify(q2 - 4 * N_Y * G2e[0, 0]) == 0 and sp.simplify(q3 - 4 * N_Y * G3e[0, 0]) == 0
      and out["even"]["psi_is"] == "metric")
q2, q3 = stat["stationary magnetic tide h = 2tK"]
check("odd Gram [psi0, psi0] reproduces the stationary magnetic tide (pure extrinsic-curvature datum)",
      sp.simplify(q2 - 4 * N_Y * G2o[0, 0]) == 0 and sp.simplify(q3 - 4 * N_Y * G3o[0, 0]) == 0
      and out["odd"]["psi_is"] == "extrinsic")
check("duality: even and odd chain-basis Gram matrices of Q_2 are identical", (G2e - G2o).is_zero_matrix)

# ------------------------------------------------------------------ profile law on monomial cone data (symbolic n)
n = sp.Symbol("n", positive=True)
k = n - 2                                             # gamma ~ v^n  <->  C ~ v^(n-2)
I = sp.integrate(v**2 * v * (1 - v) * v ** (2 * k + 2), (v, 0, 1))
law = sp.simplify((I / (k + 3) ** 2) / (I / (k + 4) ** 2) - ((n + 2) / (n + 1)) ** 2)
check("profile law: for C ~ v^(n-2), Q_2/Q_3 = ((n+2)/(n+1))^2 (symbolic n)", law == 0)

res["results"] = out
res["all_pass"] = ok
with open(sys.argv[1], "w", encoding="utf8", newline="\n") as fh:
    json.dump(res, fh, indent=1, sort_keys=True)
nf = sum(1 for c in res["checks"] if not c["pass"])
print(f"\nplant={plant} checks={len(res['checks'])} fail={nf} time={time.time() - t0:.0f} s")
