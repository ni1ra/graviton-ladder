# CHARGE-REF: the second-order area formula, own route (sympy).  Usage: python ref_area.py OUT.json
# Area side A_un[h] = - int N R2[h] - 2 (d2Area - 2 d2Vol), per 4 pi, 1/8G dropped, K = 0.
#   R2: e^2 coefficient of the scalar curvature of g = delta + e h, from Gamma = (1/2) g^-1 (dg + dg - dg) with
#       g^-1 = 1 - e h + e^2 h.h inserted as a truncated series and everything expanded by sympy (no hand-split orders).
#   Area of the coordinate sphere by the coarea form dA = sqrt(det g) sqrt(nu.g^-1.nu) dOmega; volume sqrt(det g).
import sys, json, time
import sympy as sp

OUT = sys.argv[1]
T0 = time.time()
x, y, z, e = sp.symbols("x y z e")
XS = (x, y, z)
I3 = range(3)
R = sp.Rational
RES = {"lines": []}


def say(tag, ok, text):
    s = "PASS" if ok else "FAIL"
    print(f"{s} {tag} {text}", flush=True)
    RES["lines"].append([tag, s, text])


def dfact(n):
    out = 1
    while n > 1:
        out *= n
        n -= 2
    return out


def sphere(expr):      # average over the unit sphere of a polynomial
    tot = 0
    for (a, b, c), co in sp.Poly(sp.expand(expr), x, y, z).terms():
        if a % 2 or b % 2 or c % 2:
            continue
        tot += co * R(dfact(a - 1) * dfact(b - 1) * dfact(c - 1), dfact(a + b + c + 1))
    return tot


def ball(expr, withN=False):
    tot = 0
    for (a, b, c), co in sp.Poly(sp.expand(expr), x, y, z).terms():
        if a % 2 or b % 2 or c % 2:
            continue
        n = a + b + c
        w = R(dfact(a - 1) * dfact(b - 1) * dfact(c - 1), dfact(n + 1))
        tot += co * w * (R(1, (n + 3) * (n + 5)) if withN else R(1, n + 3))
    return tot


def scalar_curv_orders(h):
    """returns (R1, R2): e and e^2 coefficients of the scalar curvature of delta + e h."""
    H = sp.Matrix(3, 3, lambda i, j: h[i][j])
    g = sp.eye(3) + e * H
    gi = sp.eye(3) - e * H + e**2 * H * H
    dg = [[[sp.diff(g[i, j], XS[k]) for k in I3] for j in I3] for i in I3]     # dg[i][j][k] = d_k g_ij
    Gam = [[[sp.expand(sum(gi[k, l] * (dg[l][j][i] + dg[l][i][j] - dg[i][j][l]) for l in I3) / 2) for j in I3] for i in I3] for k in I3]
    trunc = lambda ex: sum(sp.Poly(ex, e).coeff_monomial(e**n) * e**n for n in (1, 2)) if ex != 0 else 0
    Gam = [[[trunc(Gam[k][i][j]) for j in I3] for i in I3] for k in I3]
    Ric = [[0] * 3 for _ in I3]
    for i in I3:
        for j in I3:
            t = 0
            for k in I3:
                t += sp.diff(Gam[k][i][j], XS[k]) - sp.diff(Gam[k][i][k], XS[j])
                for m in I3:
                    t += Gam[k][k][m] * Gam[m][i][j] - Gam[k][j][m] * Gam[m][i][k]
            Ric[i][j] = t
    Rs = sp.expand(sum(gi[i, j] * Ric[i][j] for i in I3 for j in I3))
    P = sp.Poly(Rs, e)
    return P.coeff_monomial(e), P.coeff_monomial(e**2)


def second_orders(h):
    """e^2 coefficients of sqrt(det g) (volume element) and of the area element on r = 1."""
    H = sp.Matrix(3, 3, lambda i, j: h[i][j])
    nu = sp.Matrix(XS)
    tr = H.trace()
    trH2 = (H * H).trace()
    a = (nu.T * H * nu)[0]
    b = (nu.T * H * H * nu)[0]
    # generic series, done by sympy on scalars: det(1+eH) = 1 + e tr + e^2 (tr^2 - trH2)/2 + e^3 det H
    t1, t2, t3, aa, bb = sp.symbols("t1 t2 t3 aa bb")
    vol_series = sp.series(sp.sqrt(1 + e * t1 + e**2 * (t1**2 - t2) / 2 + e**3 * t3), e, 0, 3).removeO()
    nrm_series = sp.series(sp.sqrt(1 - e * aa + e**2 * bb), e, 0, 3).removeO()       # nu.g^-1.nu = 1 - e a + e^2 b + O(e^3)
    v2 = sp.Poly(vol_series, e).coeff_monomial(e**2)
    ar2 = sp.Poly(sp.expand(vol_series * nrm_series), e).coeff_monomial(e**2)
    sub = {t1: tr, t2: trH2, aa: a, bb: b}
    return sp.expand(v2.subs(sub)), sp.expand(ar2.subs(sub))


def A_area(h):
    _, R2 = scalar_curv_orders(h)
    v2, a2 = second_orders(h)
    return -ball(R2, True) - 2 * (sphere(a2) - 2 * ball(v2))


def A_formula(h, parts=False):
    dh2 = sum(sp.diff(h[i][j], XS[k])**2 for i in I3 for j in I3 for k in I3)
    hx = [sum(h[i][j] * XS[j] for j in I3) for i in I3]
    xhx = sum(XS[i] * hx[i] for i in I3)
    p = (R(1, 4) * ball(dh2, True), -R(1, 2) * sphere(sum(c**2 for c in hx)), R(1, 4) * sphere(xhx**2))
    return p if parts else sum(p)


def symgrad2(v):
    return [[sp.diff(v[j], XS[i]) + sp.diff(v[i], XS[j]) for j in I3] for i in I3]


def tadd(A, B):
    return [[sp.expand(A[i][j] + B[i][j]) for j in I3] for i in I3]


def tt_family(phi, parity, kmax=1):
    s = x**2 + y**2 + z**2
    g = [sp.diff(phi, c) for c in XS]
    if parity == "even":
        bases = [[[sp.diff(g[i], XS[j]) for j in I3] for i in I3],
                 [[XS[i] * g[j] + g[i] * XS[j] for j in I3] for i in I3],
                 [[phi if i == j else 0 for j in I3] for i in I3],
                 [[phi * XS[i] * XS[j] for j in I3] for i in I3]]
    else:
        v = [XS[1] * g[2] - XS[2] * g[1], XS[2] * g[0] - XS[0] * g[2], XS[0] * g[1] - XS[1] * g[0]]
        bases = [[[XS[i] * v[j] + v[i] * XS[j] for j in I3] for i in I3],
                 [[sp.diff(v[j], XS[i]) + sp.diff(v[i], XS[j]) for j in I3] for i in I3]]
    cs = sp.symbols(f"c0:{len(bases) * (kmax + 1)}")
    T = [[sp.expand(sum(cs[b * (kmax + 1) + k] * s**k * bases[b][i][j] for b in range(len(bases)) for k in range(kmax + 1))) for j in I3] for i in I3]
    eqs = []
    cond = [sum(T[i][i] for i in I3)] + [sum(sp.diff(T[i][j], XS[i]) for i in I3) for j in I3]
    for c in cond:
        eqs += [co for co in sp.Poly(sp.expand(c), x, y, z).coeffs()]
    sol = sp.solve(eqs, cs, dict=True)[0]
    free = [c for c in cs if c not in sol]
    out = []
    for n, f in enumerate(free):
        vals = {ff: (R(n + 2, 3) if ff == f else R(1, 5 + n)) for ff in free}     # a generic mixture, one per free parameter
        out.append([[sp.expand(T[i][j].subs(sol).subs(vals)) for j in I3] for i in I3])
    return out


def is_tt(h):
    return sp.expand(sum(h[i][i] for i in I3)) == 0 and all(sp.expand(sum(sp.diff(h[i][j], XS[i]) for i in I3)) == 0 for j in I3)


for l, phi in ((2, x * y + z**2 - x**2), (3, x * y * z + x**3 - 3 * x * y**2)):
    assert sp.expand(sum(sp.diff(phi, c, 2) for c in XS)) == 0
    NY = sphere(phi**2)
    ev = tt_family(phi, "even")
    od = tt_family(phi, "odd")
    data = {"even0": ev[0], "even1": ev[-1], "odd0": od[0], "odd1": od[-1], "even+odd": tadd(ev[-1], od[-1])}
    assert all(is_tt(h) for h in data.values())
    ok = True
    recs = []
    for nm, h in data.items():
        a1, a2 = A_area(h), A_formula(h)
        ok = ok and a1 == a2 and a1 != 0
        recs.append(f"{nm}: {a1 / NY} vs {a2 / NY}")
    say(f"A01[l={l}]", ok, "area side = four-term TT formula; " + "; ".join(recs))
    p = A_formula(ev[-1], parts=True)
    print(f"   A02[l={l}] even1 terms / N_Y: (1/4) int N |dh|^2 = {p[0] / NY}, -(1/2) oint |h nu|^2 = {p[1] / NY}, (1/4) oint (nu h nu)^2 = {p[2] / NY}; "
          f"formula without the last term misses the area side by {(A_area(ev[-1]) - p[0] - p[1]) / NY}", flush=True)

    # A03 off TT gauge
    T = ev[-1]
    hg = tadd(T, symgrad2([phi * c for c in XS]))
    hs = [[-2 * phi if i == j else 0 for j in I3] for i in I3]
    r1 = (A_area(hg) / NY, A_formula(hg) / NY)
    r2 = (A_area(hs) / NY, A_formula(hs) / NY)
    val = R(2 * l * (l - 1), (l + 2) * (2 * l + 3)) - R(4 * (l - 1), l + 2)
    say(f"A03[l={l}]", r1[0] != r1[1] and r2[0] != r2[1],
        f"off TT: T + 2 sym grad(phi x): area side {r1[0]}, four-term formula {r1[1]}; h = -2 phi delta: area side {r2[0]} (lane's X03 value {val}: {r2[0] == val}), four-term formula {r2[1]}")

    # A04 Green identity on a generic non-TT h'
    hp = [[0] * 3 for _ in I3]
    seeds = {(0, 0): x**2 * y + z, (0, 1): x * z**2 - y**3 + 2, (0, 2): x * y * z + x**2, (1, 1): y**4 - x * z, (1, 2): z**3 + x * y, (2, 2): x**2 * z**2 + y}
    for (i, j), v in seeds.items():
        hp[i][j] = hp[j][i] = v
    R1, _ = scalar_curv_orders(hp)
    tr = sum(hp[i][i] for i in I3)
    nhn = sum(XS[i] * hp[i][j] * XS[j] for i in I3 for j in I3)
    dF = R(1, 2) * sphere(tr - nhn) - ball(tr)
    if l == 2:
        say("A04", ball(R1, True) == -2 * dF and dF != 0, f"int N R_lin[h'] = {ball(R1, True)} = -2 dF[h'] = {-2 * dF}")
        # A00 (unregistered): my curvature code on an exactly flat metric, the pull-back of delta by x -> x + e xi
        xi = [x * y - z**2, y * z + x**2, x * z**2 - y]
        Jm = sp.Matrix(3, 3, lambda i, j: sp.diff(xi[i], XS[j]))
        h1 = Jm + Jm.T
        h2 = Jm.T * Jm
        hh = [[h1[i, j] + e * h2[i, j] for j in I3] for i in I3]
        H = sp.Matrix(3, 3, lambda i, j: hh[i][j])
        # scalar_curv_orders assumes h independent of e in g^-1; redo with the exact second-order inverse
        g = sp.eye(3) + e * H
        gi = sp.eye(3) - e * h1 + e**2 * (h1 * h1 - h2)
        dg = [[[sp.diff(g[i, j], XS[k]) for k in I3] for j in I3] for i in I3]
        Gam = [[[sp.expand(sum(gi[k, m] * (dg[m][j][i] + dg[m][i][j] - dg[i][j][m]) for m in I3) / 2) for j in I3] for i in I3] for k in I3]
        Rs = 0
        for i in I3:
            for j in I3:
                t = 0
                for k in I3:
                    t += sp.diff(Gam[k][i][j], XS[k]) - sp.diff(Gam[k][i][k], XS[j])
                    for m in I3:
                        t += Gam[k][k][m] * Gam[m][i][j] - Gam[k][j][m] * Gam[m][i][k]
                Rs += gi[i, j] * t
        P = sp.Poly(sp.expand(Rs), e)
        say("A00", P.coeff_monomial(e) == 0 and P.coeff_monomial(e**2) == 0, "unregistered: curvature code gives R = 0 through e^2 on an exactly flat pulled-back metric")

    # A05 edge: gauge route against the direct drive
    s = x**2 + y**2 + z**2
    gph = [sp.diff(phi, c) for c in XS]
    xi = [R(2, 5) * phi * XS[i] - R(1, 3) * gph[i] + R(1, 7) * s * gph[i] for i in I3]
    sig = sphere(sum(XS[i] * xi[i] for i in I3) * phi) / NY
    ok = sig != 0
    recs = []
    for nm in ("even0", "even1"):
        T = data[nm]
        psi = sum(XS[i] * T[i][j] * XS[j] for i in I3 for j in I3)
        Gf = (sum(XS[i] * sp.diff(psi, XS[i]) for i in I3) + 2 * psi) / 2
        Gc = sphere(Gf * phi) / NY
        lhs = A_area(tadd(T, symgrad2(xi))) - A_area(T)
        rhs = (-2 * sig * Gc - (l - 1) * (l + 2) * sig**2) * NY
        ok = ok and lhs == rhs and Gc != 0
        recs.append(f"{nm}: G = {Gc} Y, change/N_Y = {lhs / NY}")
    say(f"A05[l={l}]", ok, f"sigma = {sig} Y; " + "; ".join(recs))
    print(f"   l={l} done ({time.time() - T0:.1f}s)", flush=True)

npass = sum(1 for r_ in RES["lines"] if r_[1] == "PASS")
nfail = sum(1 for r_ in RES["lines"] if r_[1] == "FAIL")
print(f"SUMMARY pass {npass} fail {nfail} time {time.time() - T0:.1f}s")
RES["summary"] = {"pass": npass, "fail": nfail}
with open(OUT, "w", encoding="utf-8") as fh:
    json.dump(RES, fh, indent=1)
