"""ODD-SIGN: uniform-in-l lower bound for lam ||beta1hat^(l) e_l|| / ||e_l|| in the ODD tower.

D0(q, pi) = (0, -N q), N = 1/(1 + cosh rho); omega(k, D0 k) = -<q, N q>.
Family: u(y) = (1 - ((y-2)/4)^2)^4 on [-2, 6], y = rho - ln(2 lam); e_l = (q_l, 0), q_l = exp(-eps h_l^2) h_l^2 u, eps = 1.5.
Pass lines: ..\\logs\\01_prefixed_pass_lines.log (fixed before the first run).
Part C: ball arithmetic (python-flint) for the constants of the proved bound.
Part N: float64 grid (even-passage's discretisation of h_l) for the actual quantity at l = 16 ... 4096 and the limit model.
Usage: ..\\..\\even-nonzero\\venv\\Scripts\\python.exe odd_sign.py ..\\out\\odd_sign.json
"""
import json
import math
import sys
import time

import numpy as np
from flint import acb, arb, ctx

ctx.prec = 100
TAU = math.pi
EPS = 1.5
Y0, B, K = 2.0, 4.0, 4
L0 = 16
LS = [16, 64, 256, 1024, 4096]
YBOT = -4.8


def lam_of(l):
    return math.sqrt(l * (l + 1))


# ------------------------------------------------------------------ Part C
def lo(b):
    return float(b.mid()) - float(b.rad()) * (1 + 1e-9) - 1e-15 * abs(float(b.mid()))


def hi(b):
    return float(b.mid()) + float(b.rad()) * (1 + 1e-9) + 1e-15 * abs(float(b.mid()))


def integ(f):
    return acb.integral(lambda y, _: f(y), arb(Y0 - B), arb(Y0 + B)).real


def profile(y):
    s = (y - Y0) / B
    w = 1 - s * s
    return w ** 4, -2 * s * w ** 3, 3 * s * s * w * w - w ** 3 / 2     # u, u', u''


def Vl(y, l):
    e = (-2 * y).exp()
    if l is None:
        return e
    return e / (1 - e / (4 * l * (l + 1))) ** 2


def inv_lapse(y, l):
    """1/(lam N_l) = (1 + cosh rho)/lam = e^y (1 + e^-y/(2 lam))^2;  limit e^y."""
    if l is None:
        return y.exp()
    return y.exp() * (1 + (-y).exp() / (2 * arb(l * (l + 1)).sqrt())) ** 2


def AGQ(l):
    A = integ(lambda y: profile(y)[1] ** 2 + Vl(y, l) * profile(y)[0] ** 2)
    G = integ(lambda y: (-profile(y)[2] + Vl(y, l) * profile(y)[0]) ** 2).sqrt()
    Q = integ(lambda y: inv_lapse(y, l) * profile(y)[0] ** 2)           # <u, (lam N_l)^-1 u>
    return A, G, Q


def sup_cells(f, P, n, tail):
    """Upper bound for sup_{p>=0} f: ball evaluation on n cells of [0, P], and a proven tail bound beyond P."""
    best = tail
    h = P / n
    for i in range(n):
        v = f(arb((i + 0.5) * h, 0.5 * h * (1 + 1e-12)))
        best = max(best, hi(v))
    return best


def chain(A_lo, A_hi, G_hi, Q_hi, Sp, C):
    J0 = A_lo * (-(EPS * G_hi * G_hi / A_lo)).exp()
    J1 = A_hi / (arb(1).exp() * EPS)
    ss = (2 * J0 / J1).sqrt()
    assert hi(ss) < TAU, "closed form needs s* < tau"
    M = 8 * J0 * J0 * ss / (15 * Q_hi)
    Den = (Sp * A_hi).sqrt()
    E = (C * A_hi / 2).sqrt()
    return dict(J0_lo=J0, J1_hi=J1, sstar=ss, M_lo=M, Den_hi=Den, E_hi=E, c=M / (2 * Den * E))


def certified():
    PI = arb.pi()
    # S' = sup p^3 e^(-2 eps p^2) sinh(2 pi p).  Tail p >= 6: f <= (p^3/2) e^(2 pi p - 3 p^2), decreasing for p >= 2.
    Sp = sup_cells(lambda p: p ** 3 * (-2 * EPS * p * p).exp() * (2 * PI * p).sinh(), 6.0, 6144,
                   108.0 * math.exp(2 * math.pi * 6 - 2 * EPS * 36) * 1.01)
    # C = sup (p^3 + p^2/pi) e^(-2 eps p^2) >= sup p^3 coth(pi p) e^(-2 eps p^2).  Tail p >= 6: <= 228 e^(-108).
    C = sup_cells(lambda p: (p ** 3 + p * p / PI) * (-2 * EPS * p * p).exp(), 6.0, 6144, 228 * math.exp(-108.0))
    Sp, C = arb(Sp), arb(C)
    Ai, Gi, Qi = AGQ(None)
    A0, G0, Q0 = AGQ(L0)
    dV = integ(lambda y: ((Vl(y, L0) - Vl(y, None)) * profile(y)[0]) ** 2).sqrt()
    G_hi = Gi + dV
    out = dict(S_prime=Sp, C=C, A_inf=Ai, G_inf=Gi, Q_inf=Qi, A_16=A0, G_16=G0, Q_16=Q0, dV_16=dV, G_hi=G_hi)
    unif = chain(Ai, A0, G_hi, Q0, Sp, C)
    per_l = {}
    for l in LS:
        A, G, Q = AGQ(l)
        per_l[l] = dict(A=A, G=G, Q=Q, **chain(A, A, G, Q, Sp, C))
    lim = chain(Ai, Ai, Gi, Qi, Sp, C)
    return out, unif, per_l, lim


# ------------------------------------------------------------------ Part N
def grid(dy, ytop, ybot):
    n = int(round((ytop - ybot) / dy))
    return ybot + dy * np.arange(1, n + 1)


def model(l, dy=0.02, ytop=40.0, even_sign_plant=False):
    """Grid quantities for sector l (None = limit model). The lapse used is lam*N_l, so beta is lam*beta1hat."""
    if l is None:
        y = grid(dy, ytop, YBOT)
        pot, lapse = np.exp(-2 * y), np.exp(-y)
    else:
        lam = lam_of(l)
        y = grid(dy, ytop, max(YBOT, -math.log(2 * lam)))          # rho > 0
        rho = y + math.log(2 * lam)
        pot, lapse = lam ** 2 / np.sinh(rho) ** 2, lam / (1.0 + np.cosh(rho))
    n = len(y)
    Hm = np.zeros((n, n))
    i = np.arange(n)
    Hm[i, i] = 2.0 / dy ** 2 + pot
    Hm[i[:-1], i[:-1] + 1] = -1.0 / dy ** 2
    Hm[i[:-1] + 1, i[:-1]] = -1.0 / dy ** 2
    Hm[-1, -1] -= 1.0 / dy ** 2                         # Neumann top wall, Dirichlet at the bottom
    p2, U = np.linalg.eigh(Hm)
    assert p2.min() > 0
    p = np.sqrt(p2)
    del Hm
    s_ = (y - Y0) / B
    u = np.where(np.abs(s_) < 1, (1 - s_ * s_) ** K, 0.0)
    cu = U.T @ u
    damp = np.exp(-EPS * p2)
    a = p2 * damp * cu                                   # q = exp(-eps h^2) h^2 u in the eigenbasis
    res = {"n": n}
    res["A"] = dy * float(np.sum(p2 * cu * cu))
    res["G"] = math.sqrt(dy * float(np.sum(p2 * p2 * cu * cu)))
    res["J0"] = dy * float(np.sum(p2 * damp * cu * cu))
    res["J1"] = dy * float(np.sum(p2 * p2 * damp * cu * cu))
    res["Q"] = dy * float(np.sum(u * u / lapse))
    Nt = (U.T * lapse[None, :]) @ U                      # lam N in the eigenbasis (multiplication operator)
    f = lambda s: dy * float((np.cos(s * p) * a) @ (Nt @ (np.cos(s * p) * a)))    # <q_s, lam N q_s>
    res["integrand"] = [f(s) for s in np.linspace(0, TAU, 33)]
    g, wq = np.polynomial.legendre.leggauss(64)
    res["m"] = float(sum(0.5 * TAU * ww * f(0.5 * TAU * (uu + 1)) for uu, ww in zip(g, wq)))
    # z = h^(1/2) q / sqrt 2;  Den^2 = <z, 2 sinh(2 pi h) z> = <q, h sinh(2 pi h) q>, written so that nothing overflows
    res["Den"] = math.sqrt(dy * float(np.sum(np.exp(2 * math.pi * p - 2 * EPS * p2) * (-np.expm1(-4 * math.pi * p))
                                             / 2 * p ** 5 * cu * cu)))
    z = (np.sqrt(p) * a / math.sqrt(2)).astype(complex)
    zb = np.conj(z)
    isq = 1.0 / np.sqrt(p)
    M = isq[:, None] * Nt * isq[None, :]                 # h^(-1/2) lam N h^(-1/2)
    del Nt
    c2 = 1.0 / (-np.expm1(-2 * math.pi * p))
    s2 = np.exp(-2 * math.pi * p) * c2
    coth = c2 + s2

    def F(theta):                                        # int_0^tau e^(-i s theta) ds
        th = TAU * theta
        return TAU * np.exp(-0.5j * th) * np.sinc(th / (2 * math.pi))

    # D0 z = -(i/2) M (z + zbar);  e^{sX0} z = e^{-ish} z.  T z = A z + B zbar.
    Am = -0.5j * M * F(p[:, None] - p[None, :])
    Az, Ac, As = Am @ z, Am @ (c2 * z), Am @ (s2 * z)
    del Am
    Bm = (0.5j if even_sign_plant else -0.5j) * M * F(p[:, None] + p[None, :])
    Bz, Bc, Bs = Bm @ zb, Bm @ (c2 * zb), Bm @ (s2 * zb)
    del Bm, M
    Tz = Az + Bz
    Pv, Qv = Bc - As, Ac - Bs
    res["Im_z_Tz"] = dy * float(np.imag(np.vdot(z, Tz)))
    res["beta"] = math.sqrt(dy * float(np.sum(c2 * np.abs(Pv) ** 2 + s2 * np.abs(Qv) ** 2)))
    res["e"] = math.sqrt(dy * float(np.sum(coth * np.abs(z) ** 2)))
    res["r"] = res["beta"] / res["e"]
    res["chain"] = res["m"] / (2 * res["Den"] * res["e"])
    return res


def main():
    t0 = time.time()
    out, unif, per_l, lim = certified()
    print("== Part C: ball arithmetic (prec %d bits) ==" % ctx.prec)
    for k, v in out.items():
        print(f"  {k:10s} {v}")
    print("  uniform chain (all l >= 16):")
    for k, v in unif.items():
        print(f"    {k:8s} {v}")
    c_unif = lo(unif["c"])
    J0_lo, J1_hi, M_lo = lo(unif["J0_lo"]), hi(unif["J1_hi"]), lo(unif["M_lo"])
    Den_hi, E_hi = hi(unif["Den_hi"]), hi(unif["E_hi"])
    A_inf, A_16, Q_16 = lo(out["A_inf"]), hi(out["A_16"]), hi(out["Q_16"])
    print(f"  PROVED UNIFORM FLOOR  lam ||beta1hat e_l||/||e_l|| >= c_unif = {c_unif:.6e}  for all l >= {L0}")
    print(f"  limit-model bound (l = inf constants): {lo(lim['c']):.6e}")
    lines = {}
    lines["C1"] = c_unif > 0 and hi(unif["sstar"]) < TAU and all(math.isfinite(hi(v)) for v in unif.values())
    lines["C2"] = 1.2e-3 <= c_unif <= 2.6e-3
    lines["C3"] = (25 <= lo(out["Q_16"]) and Q_16 <= 40 and 24.5 <= lo(out["S_prime"]) and hi(out["S_prime"]) <= 26
                   and 0.113 <= lo(out["C"]) and hi(out["C"]) <= 0.118)

    print("== Part N: float64 grid, dy = 0.02 ==", flush=True)
    R = {}
    for l in LS + [None]:
        R[l] = model(l)
        r = dict(R[l])
        r.pop("integrand")
        print(f"  l={'inf' if l is None else l:>5} :", json.dumps(r), flush=True)
    sgrid = np.linspace(0, TAU, 33)
    floor = np.maximum(J0_lo - sgrid ** 2 * J1_hi / 2, 0.0) ** 2 / Q_16
    fl = {k: [] for k in ("P1", "P2", "P3", "P4", "P5", "P6")}
    for l in LS:
        r = R[l]
        fl["P1"].append(A_inf * (1 - 2e-3) <= r["A"] <= A_16 * (1 + 2e-3) and r["J0"] >= J0_lo and r["J1"] <= J1_hi)
        fl["P2"].append(bool(np.all(np.array(r["integrand"]) >= floor)))
        fl["P3"].append(r["m"] >= M_lo)
        fl["P4"].append(r["Den"] <= Den_hi and r["e"] <= E_hi)
        fl["P5"].append(r["r"] >= c_unif and r["r"] >= r["chain"] * (1 - 1e-9))
        fl["P6"].append(abs(r["Im_z_Tz"] + r["m"] / 2) <= 1e-6 * r["m"])
    for k, v in fl.items():
        lines[k] = all(v)
        if not all(v):
            print(f"  {k} per l: {v}")
    rinf = R[None]["r"]
    lines["N1"] = R[4096]["r"] / R[16]["r"] >= 0.8
    d = {l: lam_of(l) * abs(R[l]["r"] / rinf - 1) for l in LS}
    dd = [d[l] for l in (64, 256, 1024, 4096)]
    lines["N2"] = max(dd) / min(dd) <= 2
    fine_l, fine_inf, big_inf = model(1024, dy=0.01), model(None, dy=0.01), model(None, ytop=80.0)
    n3 = {"dy_l1024": abs(fine_l["r"] / R[1024]["r"] - 1), "dy_inf": abs(fine_inf["r"] / rinf - 1),
          "box_inf": abs(big_inf["r"] / rinf - 1)}
    lines["N3"] = n3["dy_l1024"] < 1e-3 and n3["dy_inf"] < 1e-3 and n3["box_inf"] < 1e-6
    lines["N4"] = all(5 <= R[l]["r"] / c_unif <= 400 for l in LS)

    # plants
    pA = {l: R[l]["r"] / lam_of(l) for l in LS}                      # factor lam dropped
    plantA_P5_fail = [pA[l] < c_unif for l in LS]
    plantA_N1_fail = pA[4096] / pA[16] < 0.8
    plantA_caught = plantA_P5_fail[3] and plantA_P5_fail[4] and plantA_N1_fail
    pB = model(64, even_sign_plant=True)                            # wrong slot: even-parity sign on the B part
    pB_p6 = abs(pB["Im_z_Tz"] + pB["m"] / 2) / pB["m"]
    pB_dr = abs(pB["r"] / R[64]["r"] - 1)
    plantB_caught = pB_p6 > 1e-2 and pB_dr > 1e-2

    print("== table ==")
    print(f"  c_unif (proved, all l >= 16) = {c_unif:.6e};  M_lo = {M_lo:.6f}; Den_hi = {Den_hi:.5f}; E_hi = {E_hi:.5f}")
    print("      l   r_l=lam||b1 e||/||e||   proved c_unif   proved c(l)   lam*m_l    Den_l     ||e_l||   chain(m/2Den e)  "
          "lam|r/rinf-1|  -2Im<z,Tz>/m")
    for l in LS:
        r = R[l]
        print(f"  {l:5d}   {r['r']:.8f}            {c_unif:.6e}    {lo(per_l[l]['c']):.6e}  {r['m']:.6f}  {r['Den']:.5f}  "
              f"{r['e']:.5f}   {r['chain']:.6f}         {d[l]:.4f}       {-2 * r['Im_z_Tz'] / r['m']:.9f}")
    print(f"    inf   {rinf:.8f}            (limit model)   {lo(lim['c']):.6e}  {R[None]['m']:.6f}  {R[None]['Den']:.5f}  "
          f"{R[None]['e']:.5f}   {R[None]['chain']:.6f}")
    print(f"  N3 relative changes: {n3}")
    print(f"  min over l, s of integrand/floor (where floor > 0): "
          f"{min(min(a / b for a, b in zip(R[l]['integrand'], floor) if b > 0) for l in LS):.3f}")
    print(f"  grid Q_l: {[f'{R[l]['Q']:.5f}' for l in LS]} (certified Q_16 <= {Q_16:.5f})")
    print(f"  PLANT A (lam dropped): values {[f'{pA[l]:.3e}' for l in LS]}; below c_unif: {plantA_P5_fail}; "
          f"ratio 4096/16 = {pA[4096] / pA[16]:.4f}; caught = {plantA_caught}")
    print(f"  PLANT B (even-parity sign on B, l = 64): |Im<z,Tz> + m/2|/m = {pB_p6:.4f} (line P6 allows 1e-6); "
          f"r = {pB['r']:.6f} vs {R[64]['r']:.6f} (shift {pB_dr:.4f}); caught = {plantB_caught}")
    lines["PLANT_A_caught"] = plantA_caught
    lines["PLANT_B_caught"] = plantB_caught
    print("== pass lines ==")
    for k, v in lines.items():
        print(f"  {k}: {'PASS' if v else 'FAIL'}")
    ok = all(lines.values())
    print("ALL LINES HELD" if ok else "SOME LINES FAILED", f"({time.time() - t0:.0f} s)")
    if len(sys.argv) > 1:
        ball = lambda v: {"ball": str(v), "lo": lo(v), "hi": hi(v)}
        with open(sys.argv[1], "w") as f:
            json.dump({"constants": {k: ball(v) for k, v in out.items()},
                       "uniform": {k: ball(v) for k, v in unif.items()},
                       "per_l": {str(l): {k: ball(v) for k, v in per_l[l].items()} for l in LS},
                       "limit": {k: ball(v) for k, v in lim.items()},
                       "grid": {("inf" if l is None else str(l)): R[l] for l in R},
                       "d": {str(l): d[l] for l in LS}, "n3": n3,
                       "plantA": {str(l): pA[l] for l in LS},
                       "plantB": {"p6": pB_p6, "dr": pB_dr, "r": pB["r"]},
                       "lines": {k: bool(v) for k, v in lines.items()}}, f, indent=1)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
