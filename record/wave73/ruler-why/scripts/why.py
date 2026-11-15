"""Lane RULER-WHY (wave 73): which ingredient of the fixed-volume rule makes the metric-versus-
extrinsic-curvature asymmetry?  Exact rational arithmetic (sympy).

Units 1/(8G): area side = -2 x (eps^2 coefficient of the edge area at fixed size).  Unit ball.
Electric datum: exact static Weyl metric, U = eps r^l P_l(cos theta), K = 0   (pure metric data).
Magnetic datum: stationary tide h = 2tK, K = sym grad(x cross grad phi)        (pure K data).

Usage: python why.py OUT.json [--plant k1] [--lmax 3]
"""
import json
import os
import sys
import time

import sympy as sp

import nc_core as nc
from nc_core import v, x, y, z, XS

plant = sys.argv[sys.argv.index("--plant") + 1] if "--plant" in sys.argv else None
LS = (2, 3, 4)
KVOL = 1 if plant == "k1" else 2          # the "2" of Area - 2 Volume
t0 = time.time()
res = {"script": "why.py", "plant": plant, "checks": [], "table": {}, "v3f": {}}
ok = True
R_ = sp.Rational


def check(name, cond, detail=""):
    global ok
    cond = bool(cond)
    ok &= cond
    res["checks"].append({"name": name, "pass": cond, "detail": str(detail)})
    print(("PASS " if cond else "FAIL ") + name + (f" | {detail}" if detail != "" else ""), flush=True)


eps, a, c, r, X, rho, zz = sp.symbols("eps a c r X rho zz", real=True)


def trunc(e, order=2):
    e = sp.expand(e)
    return sum(e.coeff(eps, k) * eps**k for k in range(order + 1))


def xint(e):
    return sp.integrate(sp.expand(e), (X, -1, 1))


# ====================================================================== electric, brute force
def electric(l):
    P = sp.expand(sp.legendre(l, X))
    dP = sp.diff(P, X)
    NY = 4 * sp.pi / (2 * l + 1)
    # Weyl quadrature for k2 (k = eps^2 k2): k_rho = rho (phi_rho^2 - phi_z^2), k_z = 2 rho phi_rho phi_z
    leg = sp.Poly(sp.legendre(l, X), X)
    phi_c = sp.expand(sum(co * zz**j * (rho**2 + zz**2) ** ((l - j) // 2) for (j,), co in leg.terms()))
    assert sp.expand(sp.diff(phi_c, rho, 2) + sp.diff(phi_c, rho) / rho + sp.diff(phi_c, zz, 2)) == 0
    k_rho = sp.expand(rho * (sp.diff(phi_c, rho) ** 2 - sp.diff(phi_c, zz) ** 2))
    k2c = sp.integrate(k_rho, (rho, 0, rho))
    check(f"l={l}: Weyl quadrature is integrable and k vanishes on the axis",
          sp.expand(sp.diff(k2c, zz) - 2 * rho * sp.diff(phi_c, rho) * sp.diff(phi_c, zz)) == 0 and k2c.subs(rho, 0) == 0)
    # to (r, X): rho^2 = r^2 (1 - X^2), z = r X   (k2c is even in rho)
    k2 = sp.expand(sp.Poly(k2c, rho, zz).as_expr())
    k2 = sp.expand(sum(co * (r**2 * (1 - X**2)) ** (i // 2) * (r * X) ** j
                       for (i, j), co in sp.Poly(k2c, rho, zz).terms()))
    rs = 1 + eps * a * P + eps**2 * c
    rsx = eps * a * dP
    # area element: rs e^{-2U+k} sqrt(rs^2 + (1-X^2) rs_X^2) dX dphi at r = rs
    conf = trunc(1 - 2 * eps * rs**l * P + 2 * eps**2 * P**2 + eps**2 * k2.subs(r, 1))
    root = trunc(rs + (1 - X**2) * rsx**2 / 2)                 # sqrt(rs^2 + q) = rs + q/(2 rs) + ...
    A = trunc(rs * conf * root)
    # volume element r^2 e^{-3U+2k}
    dens = trunc(1 - 3 * eps * r**l * P + R_(9, 2) * eps**2 * r ** (2 * l) * P**2 + 2 * eps**2 * k2)
    Vr = sp.integrate(sp.expand(r**2 * dens), r)
    V = trunc(Vr.subs(r, rs))
    # geodesic distance from the centre to the edge point at angle X: radial length minus bending
    lin = trunc(1 - eps * r**l * P + eps**2 * r ** (2 * l) * P**2 / 2 + eps**2 * k2)
    Dr = sp.integrate(sp.expand(lin), r)
    bend = -eps**2 * (1 - X**2) * dP**2 / (2 * (l + 1) ** 2 * (2 * l + 1))
    D = trunc(Dr.subs(r, rs) + bend)
    # mean over initial directions: Jacobian of the first-order deflection Delta = -eps P_theta/(l(l+1))
    D1 = D.coeff(eps, 1)
    jac = eps**2 * (1 - X**2) * (-dP / (l * (l + 1))) * sp.diff(D1, X)
    out = {}
    A0, A1, A2 = [2 * sp.pi * xint(A.coeff(eps, k)) for k in range(3)]
    V0, V1, V2 = [2 * sp.pi * xint(V.coeff(eps, k)) for k in range(3)]
    Zs = [xint(D.coeff(eps, k)) / 2 for k in range(3)]
    Z2 = sp.expand(Zs[2] + xint(jac.coeff(eps, 2)) / 2)
    Z2_nobend = sp.expand(Z2 - xint(bend.coeff(eps, 2)) / 2)
    check(f"l={l}: zeroth and first order: A = 4 pi, V = 4 pi/3, radius 1, no eps^1 terms",
          A0 == 4 * sp.pi and V0 == 4 * sp.pi / 3 and Zs[0] == 1 and A1 == 0 and V1 == 0 and Zs[1] == 0)
    check(f"l={l}: size shift c enters as 8 pi c, 4 pi c, c",
          sp.diff(A2, c) == 8 * sp.pi and sp.diff(V2, c) == 4 * sp.pi and sp.diff(Z2, c) == 1)
    A2, V2, Z2, Z2_nobend = [e.subs(c, 0) for e in (A2, V2, Z2, Z2_nobend)]
    S = R_(2 * l * (l - 1), (l + 2) * (2 * l + 3)) * NY
    # --- D-DEP's expansion (inherited), per N_Y, evaluated on the generic harmonic
    phi = nc.harmonic(l)
    ny = nc.surf_int(phi**2)
    Nl = (1 - x**2 - y**2 - z**2) / 2
    F2_ddep = ((-a**2 + R_(l * (l + 1), 2) * a**2 - 2 * a * (l - 1) + 2)
               - 9 * nc.ball_int(phi**2) / ny - nc.ball_int(Nl * sum(sp.diff(phi, q) ** 2 for q in XS)) / ny)
    if plant is None:
        check(f"l={l}: brute-force Area - 2 Volume from the exact Weyl metric equals D-DEP's F2(a) for every a",
              sp.expand((A2 - 2 * V2) / NY - F2_ddep) == 0, sp.factor((A2 - 2 * V2) / NY))
        S_nc = nc.entropy_slice(nc.Curv(nc.metric_electric(phi)), l) / ny
        check(f"l={l}: electric entropy from nc.entropy_slice equals 2l(l-1)/((l+2)(2l+3)) N_Y",
              sp.simplify(S_nc - S / NY) == 0, S_nc)
    Fvol = sp.expand(A2 - KVOL * V2)
    Frad = sp.expand(A2 - 8 * sp.pi * Z2)
    a_vol = sp.solve(sp.diff(Fvol, a), a)[0]
    a_rad = sp.solve(sp.diff(Frad, a), a)[0]

    def ratio(F, aa):
        return sp.nsimplify(sp.simplify(-2 * F.subs(a, aa) / S), rational=True)

    out["(a) maximal slice, fixed volume, edge extremised"] = ratio(Fvol, a_vol)
    out["(b) fixed volume, geodesic ball a=1/(l+1)"] = ratio(Fvol, R_(1, l + 1))
    out["(b') fixed volume, diamond edge a=2/(l+1)"] = ratio(Fvol, R_(2, l + 1))
    out["(b'') fixed volume, isotropic-coordinate sphere a=0"] = ratio(Fvol, 0)
    out["(c1) geodesic ball, bare area at fixed geodesic radius"] = ratio(Frad, R_(1, l + 1))
    out["(c2) fixed mean geodesic radius, edge extremised"] = ratio(Frad, a_rad)
    extra = {"a_star fixed volume": str(a_vol), "a_star fixed mean radius": str(a_rad),
             "stiffness of a^2 (fixed volume, fixed radius), per N_Y": [str(sp.expand(Fvol / NY).coeff(a, 2)),
                                                                         str(sp.expand(Frad / NY).coeff(a, 2))],
             "(c1) without the geodesic-bending term": str(sp.nsimplify(sp.simplify(
                 -2 * (A2 - 8 * sp.pi * Z2_nobend).subs(a, R_(1, l + 1)) / S), rational=True)),
             "-2 Fvol(a)/S": str(sp.factor(sp.simplify(-2 * Fvol / S))),
             "-2 Frad(a)/S": str(sp.factor(sp.simplify(-2 * Frad / S)))}
    if plant is None:
        check(f"l={l}: stationary edge at fixed volume is a = 2/(l+2)", a_vol == R_(2, l + 2), a_vol)
    return out, extra


# ====================================================================== magnetic, l = 0 constraint
def magnetic(l):
    phi = nc.harmonic(l)
    K = nc.sym_grad(nc.vec_potential(phi))
    K2 = sp.expand(sum(K[i, j] ** 2 for i in range(3) for j in range(3)))
    xKx = sp.expand(sum(XS[i] * XS[j] * K[i, j] for i in range(3) for j in range(3)))
    check(f"l={l}: magnetic K is TT with x.K.x = 0: no scalar on the slice can drive an edge shape or a time shift",
          nc.is_tt(K) and xKx == 0)
    kap = nc.surf_int(K2) / (4 * sp.pi)                  # <|K|^2>(r) = kap r^(2l-2)
    # areal gauge g2 = A(r) dr^2, A = 2m/r:  R1[A x x / r^2] = (2/r^2) (r A)'   (checked on A = r^2, r^4)
    for p in (1, 2):
        s_ = x**2 + y**2 + z**2
        g = sp.Matrix(3, 3, lambda i, j: s_ ** (p - 1) * XS[i] * XS[j])          # A = r^(2p)
        R1 = sp.expand(sum(sp.diff(g[i, j], XS[i], XS[j]) for i in range(3) for j in range(3)) - nc.lap(g.trace()))
        check(f"l={l}: linearised scalar curvature of A(r) dr^2 is (2/r^2)(rA)' for A = r^{2 * p}",
              sp.expand(R1 - 2 * (2 * p + 1) * s_ ** (p - 1)) == 0)
    m = kap * r ** (2 * l + 1) / (4 * (2 * l + 1))        # 4 m'/r^2 = kap r^(2l-2)
    V2 = 4 * sp.pi * sp.integrate(r * m, (r, 0, 1))       # volume inside areal radius 1
    Z2 = sp.integrate(m / r, (r, 0, 1))                   # geodesic radius of areal radius 1
    A2 = 0
    Nl = (1 - x**2 - y**2 - z**2) / 2
    if plant is None:
        check(f"l={l}: mass-function value of -2(A - 2V) equals int_B N |K|^2 (D-DEP 2.1)",
              sp.simplify(-2 * (A2 - 2 * V2) - nc.ball_int(Nl * K2)) == 0)
    S = nc.entropy_slice(nc.Curv(nc.tt_series(K, 1)), l)
    vol = sp.nsimplify(sp.simplify(-2 * (A2 - KVOL * V2) / S), rational=True)
    rad = sp.nsimplify(sp.simplify(-2 * (A2 - 8 * sp.pi * Z2) / S), rational=True)
    return {"(a) maximal slice, fixed volume, edge extremised": vol,
            "(b) fixed volume, geodesic ball a=1/(l+1)": vol,
            "(b') fixed volume, diamond edge a=2/(l+1)": vol,
            "(b'') fixed volume, isotropic-coordinate sphere a=0": vol,
            "(c1) geodesic ball, bare area at fixed geodesic radius": rad,
            "(c2) fixed mean geodesic radius, edge extremised": rad}


for l in LS:
    e, extra = electric(l)
    m_ = magnetic(l)
    res["table"][f"l={l}"] = {k: {"electric (metric)": str(e[k]), "magnetic (K)": str(m_[k])} for k in e}
    res["table"][f"l={l}"]["_electric detail"] = extra
    print(f"   l={l} done at {time.time() - t0:.0f} s", flush=True)
    A_ = "(a) maximal slice, fixed volume, edge extremised"
    if plant is None:
        check(f"l={l}: case (a) reproduces 1 and (l+2)/(l+1)", e[A_] == 1 and m_[A_] == R_(l + 2, l + 1), (e[A_], m_[A_]))
        check(f"l={l}: (b) geodesic ball at fixed volume gives (l+2)/(2(l+1)^2) electric (D-DEP/NULL-CONE: 2/9, 5/32)",
              e["(b) fixed volume, geodesic ball a=1/(l+1)"] == R_(l + 2, 2 * (l + 1) ** 2))
        if l == 2:
            check("l=2: (c1) reproduces Wang's printed geodesic-ball bare area 14/45 (electric), 28/15 (magnetic)",
                  e["(c1) geodesic ball, bare area at fixed geodesic radius"] == R_(14, 45)
                  and m_["(c1) geodesic ball, bare area at fixed geodesic radius"] == R_(28, 15),
                  (e["(c1) geodesic ball, bare area at fixed geodesic radius"],
                   m_["(c1) geodesic ball, bare area at fixed geodesic radius"]))
    else:
        check(f"PLANT {plant} caught at l={l}: case (a) no longer gives 1 and (l+2)/(l+1)",
              e[A_] != 1 and m_[A_] != R_(l + 2, l + 1), (e[A_], m_[A_]))

# ====================================================================== v^3 f on the cone, l = 2, 3, 4
if plant is None:
    W3 = v**2 * nc.F_CONE
    W4 = v**3 * nc.F_CONE
    for l in (2, 3, 4, 5):
        phi = nc.harmonic(l)
        K0 = nc.sym_grad(nc.vec_potential(phi))
        data = {"static electric (metric)": nc.metric_electric(phi),
                "stationary magnetic h=2tK (K)": nc.tt_series(K0, 1),
                "odd metric rung, psi_dot = phi": nc.tt_series(nc.rung_vector(phi, l), 0),
                "even K rung, psi_dot = phi (fixed volume: 1 + 1/(2l+3))": nc.tt_series(nc.rung_scalar(phi, l), 1)}
        row = {}
        for name, h in data.items():
            cv = nc.Curv(h)
            assert cv.ricci_zero()
            C = cv.cone_weyl()
            s2 = nc.norm2(nc.moment(C, 2))
            q3 = nc.cone_int(nc.norm2(nc.moment(C, 3)), W3)
            row[name] = {"Q2[v^2 f]/S": str(sp.cancel(nc.cone_int(s2, W3) / q3)),
                         "Q2[v^3 f]/S": str(sp.cancel(nc.cone_int(s2, W4) / q3))}
            if name.startswith("odd"):
                B = cv.B()
                pb = sp.expand(sum(XS[i] * XS[j] * B[i, j] for i in range(3) for j in range(3)))
                check(f"l={l}: odd metric rung has psi_B = 0 and psi_B_dot proportional to phi at t = 0",
                      pb.subs(nc.t, 0) == 0 and sp.cancel(sp.diff(pb, nc.t).subs(nc.t, 0) / phi).is_number)
            print(f"   v^3 f l={l} {name}: {row[name]} at {time.time() - t0:.0f} s", flush=True)
        res["v3f"][f"l={l}"] = row
        check(f"l={l}: v^3 f on both static tides gives (l+2)/(l+1)",
              all(row[k]["Q2[v^3 f]/S"] == str(R_(l + 2, l + 1)) for k in list(data)[:2]))
        res["v3f"][f"l={l}"]["rung equals fixed-volume value 1"] = row["odd metric rung, psi_dot = phi"]["Q2[v^3 f]/S"] == "1"
    n = sp.Symbol("n", positive=True)
    k = n - 2
    num = sp.integrate(v**3 * v * (1 - v) * v ** (2 * k + 2), (v, 0, 1)) / (k + 3) ** 2
    den = sp.integrate(v**2 * v * (1 - v) * v ** (2 * k + 2), (v, 0, 1)) / (k + 4) ** 2
    check("profile law: for cone data gamma ~ v^n, Q2[v^3 f]/Q3[v^2 f] = (n+2)/(n+1) (symbolic n)",
          sp.simplify(num / den - (n + 2) / (n + 1)) == 0)
    # form test at l = 2: RULER's own plant run (read only) against the fixed-volume odd form c_l E_2
    pj = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "ruler", "out", "ruler_plant_w3.json")
    if os.path.exists(pj):
        d = json.load(open(pj, encoding="utf8"))["results"]["odd"]
        Q = sp.Matrix([[R_(e_) for e_ in rw] for rw in d["Q2"]])
        S = sp.Matrix([[R_(e_) for e_ in rw] for rw in d["S"]])
        phi = nc.harmonic(2)
        ny = nc.surf_int(phi**2)
        s_ = x**2 + y**2 + z**2
        cl = R_(1, 6)
        E2x = sp.zeros(8, 8)                               # A_fixed-volume(odd) - S = (c_l/2) int psi_B^2
        for i in range(4):
            for j in range(4):
                E2x[i, j] = cl / 2 * nc.ball_int(s_ ** (i + j) * phi**2) / ny
        diff = Q - S - E2x
        res["v3f"]["l=2 form test"] = {
            "source": "ruler/out/ruler_plant_w3.json (inherited) minus c_l E_2 (computed here)",
            "entries where Q2[v^3 f] equals the fixed-volume odd form": [[i, j] for i in range(8) for j in range(i, 8) if diff[i, j] == 0],
            "is identity of forms": bool(diff.is_zero_matrix),
            "psi block row 0 of the difference": [str(e_) for e_ in diff[0, :4]]}
        print("   l=2 form test:", res["v3f"]["l=2 form test"], flush=True)

# ====================================================================== print
print("\nArea side / S.  Rows = switch settings; columns = electric (metric data), magnetic (K data)")
for l in LS:
    print(f"\n  l = {l}")
    for k_, val in res["table"][f"l={l}"].items():
        if k_.startswith("_"):
            continue
        print(f"    {k_:58s} {val['electric (metric)']:>10s} {val['magnetic (K)']:>10s}")
    for k_, val in res["table"][f"l={l}"]["_electric detail"].items():
        print(f"      [{k_}] {val}")
for l_, row in res["v3f"].items():
    print(f"  v^3 f {l_}: {row}")
res["all_pass"] = ok
with open(sys.argv[1], "w", encoding="utf8", newline="\n") as fh:
    json.dump(res, fh, indent=1, sort_keys=True)
nf = sum(1 for c_ in res["checks"] if not c_["pass"])
print(f"\nplant={plant} checks={len(res['checks'])} fail={nf} time={time.time() - t0:.0f} s")
