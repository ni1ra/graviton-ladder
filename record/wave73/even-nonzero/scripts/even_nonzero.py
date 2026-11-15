"""EVEN-NONZERO: beta1hat e != 0 in the Liouville limit model of the even-parity flow.

Pass lines: ..\\logs\\01_prefixed_pass_lines.log (fixed before the first run).
Part C: certified enclosures (python-flint ball arithmetic, hand-bounded head and tail).
Part N: float64 finite differences on a common y-grid: limit model, finite l = 64, 256, 1024, two plants.
Usage: ..\\venv\\Scripts\\python.exe even_nonzero.py ..\\out\\even_nonzero.json
"""
import json
import math
import sys
import time

import numpy as np
from flint import acb, arb, ctx

ctx.prec = 120
TAU = math.pi            # tau = 2 pi t, t = 0.5
EPS4 = 0.25              # phi(p) = exp(-p^2/4)
LS = [64, 256, 1024]
YBOT = -4.8              # Dirichlet wall; x = e^4.8 = 121, all vectors below e^-100 there


# ------------------------------------------------------------------ Part C
def lo(b):
    return float((arb(b.mid()) - arb(b.rad())).mid()) if not hasattr(b, "lower") else float(b.lower())


def hi(b):
    return float((arb(b.mid()) + arb(b.rad())).mid()) if not hasattr(b, "upper") else float(b.upper())


def enc(f, a, b, head, tail):
    """Enclosure of int_0^inf f = [0, head] + int_a^b f + [0, tail]; head, tail are proven upper bounds (floats)."""
    core = acb.integral(lambda p, _: f(p), arb(a), arb(b)).real
    extra = arb(0.5 * (head + tail) * (1 + 1e-12), 0.5 * (head + tail) * (1 + 1e-12))
    return core + extra


def certified():
    PI = arb.pi()
    d, P = 2.0 ** -10, 40.0
    sinh_ = lambda p: (PI * p).sinh()
    cosh_ = lambda p: (PI * p).cosh()
    phi = lambda p: (-(p * p) / 4).exp()
    # head: p^k/sinh(pi p) <= p^(k-1)/pi on [0,d];  tail: sinh(pi p) >= e^(pi p)/4 for p >= 1 and
    # int_P^inf p^k e^(-pi p) dp <= P^k e^(-pi P)/(pi - k/P).  phi <= 1 everywhere.
    head = lambda k: d ** k / (k * math.pi)
    tail = lambda k: 4 * P ** k * math.exp(-math.pi * P) / (math.pi - k / P)
    out = {}
    for k, exact in ((3, 0.25), (5, 0.5), (7, 17 / 8)):
        b = 2 * enc(lambda p, k=k: p ** k / sinh_(p), d, P, head(k), tail(k))
        out[f"moment_{k}"] = b
        out[f"moment_{k}_contains_exact"] = bool(b.contains(arb(exact)))
    I0 = 2 * enc(lambda p: phi(p) * p ** 3 / sinh_(p), d, P, head(3), tail(3))                # <pi1, pi0>
    N1sq = 2 * enc(lambda p: phi(p) ** 2 * p ** 3 / sinh_(p), d, P, head(3), tail(3))         # ||pi1||^2
    H1sq = 2 * enc(lambda p: phi(p) ** 2 * p ** 7 / sinh_(p), d, P, head(7), tail(7))         # ||h^2 pi1||^2
    # ||e||_H^2 = int phi^2 p^2 cosh/sinh^2;  head <= d cosh(pi d)/pi^2;  tail: coth <= 2, as above with k = 2
    Esq = enc(lambda p: phi(p) ** 2 * p * p * cosh_(p) / sinh_(p) ** 2, d, P,
              d * math.cosh(math.pi * d) / math.pi ** 2, 2 * tail(2))
    # Den^2 = ||(2 sinh 2 pi h)^(1/2) z||^2 = 4 int phi^2 p^2 cosh(pi p) dp, entire integrand, no head;
    # tail: 4 int_P p^2 e^(-p^2/2 + pi p) <= 4 P^2 e^(-(P/2 - pi) P)/((P/2 - pi) - 2/P)
    a = P / 2 - math.pi
    Dsq = enc(lambda p: 4 * phi(p) ** 2 * p * p * cosh_(p), 0, P, 0.0, 4 * P * P * math.exp(-a * P) / (a - 2 / P))
    w_lo = 2 * I0 * I0                      # <pi1, W pi1> >= <pi1,pi0>^2 / <pi0, W^-1 pi0> = 2 <pi1,pi0>^2
    kappa = (H1sq * N1sq).sqrt() / 2        # |<a_s,W a_s> - <a,W a>| <= kappa s^2
    sstar = (w_lo / kappa).sqrt()
    assert hi(sstar) < TAU, "s* must lie below tau for the closed form used"
    m_lo = 2 * w_lo * sstar / 3             # int_0^s* (w_lo - kappa s^2) ds; integrand >= 0 beyond s*
    m_hi = PI * N1sq / 2                    # <a_s, W a_s> <= ||W|| ||a||^2 = ||a||^2/2
    ratio_lo = (m_lo / 2) / (Dsq * Esq).sqrt()
    out.update(I0=I0, N1sq=N1sq, H1sq=H1sq, Esq=Esq, Dsq=Dsq, w_lo=w_lo, kappa=kappa, sstar=sstar,
               m_lo=m_lo, m_hi=m_hi, ratio_lo=ratio_lo)
    return out


# ------------------------------------------------------------------ Part N
def logt(rho):
    return -2.0 * np.arctanh(np.exp(-rho))


def kern_l(l, y):
    """lam * W_l(y_i, y_j), exact kernel, rho = y + ln(2 lam)."""
    lam = math.sqrt(l * (l + 1))
    a = logt(y + math.log(2 * lam))
    lo_, hi_ = np.minimum(a[:, None], a[None, :]), np.maximum(a[:, None], a[None, :])
    return lam * np.exp((l + 1) * lo_ - l * hi_) * (-np.expm1((2 * l + 1) * hi_)) / (2 * l + 1)


def kern_inf(y):
    x = np.exp(-y)
    lo_, hi_ = np.minimum(x[:, None], x[None, :]), np.maximum(x[:, None], x[None, :])
    return 0.5 * np.exp(-(hi_ - lo_)) * (-np.expm1(-2.0 * lo_))


def model(y, dy, pot, Wmat, plant_coth=False, checks=False):
    n = len(y)
    Hm = np.zeros((n, n))
    i = np.arange(n)
    Hm[i, i] = 2.0 / dy ** 2 + pot
    Hm[i[:-1], i[:-1] + 1] = -1.0 / dy ** 2
    Hm[i[:-1] + 1, i[:-1]] = -1.0 / dy ** 2
    Hm[-1, -1] -= 1.0 / dy ** 2                         # Neumann top wall, Dirichlet at YBOT
    p2, U = np.linalg.eigh(Hm)
    assert p2.min() > 0
    p = np.sqrt(p2)
    x = np.exp(-y)
    pi0 = x * np.exp(-x)
    res = {}
    if checks:
        w = Wmat @ pi0
        hp = Hm @ pi0
        res["grid_norm_pi0"] = dy * pi0 @ pi0
        res["grid_pi0_W_pi0"] = dy * pi0 @ w
        res["grid_pi0_h2_pi0"] = dy * pi0 @ hp
        res["grid_norm_h2pi0"] = dy * hp @ hp
    a = np.exp(-EPS4 * p2) * (U.T @ pi0)                # pi1 in the eigenbasis
    Wt = U.T @ Wmat @ U
    del Hm
    # m(tau) by Gauss-Legendre in s
    u, wq = np.polynomial.legendre.leggauss(64)
    m = 0.0
    for uu, ww in zip(u, wq):
        s = 0.5 * TAU * (uu + 1)
        v = np.cos(s * p) * a
        m += 0.5 * TAU * ww * dy * (v @ (Wt @ v))
    z = 1j * a / np.sqrt(2 * p)
    zb = np.conj(z)
    sq = np.sqrt(p)
    M = sq[:, None] * Wt * sq[None, :]
    if plant_coth:
        c2, s2, coth = np.ones(n), np.zeros(n), np.ones(n)
    else:
        c2 = 1.0 / (-np.expm1(-2 * math.pi * p))
        s2 = np.exp(-2 * math.pi * p) * c2
        coth = c2 + s2

    def F(theta):                                        # int_0^tau e^(-i s theta) ds, theta = matrix
        th = TAU * theta
        return TAU * np.exp(-0.5j * th) * np.sinc(th / (2 * math.pi))

    A = -0.5j * M * F(p[:, None] - p[None, :])
    Az, Ac, As = A @ z, A @ (c2 * z), A @ (s2 * z)
    del A
    B = 0.5j * M * F(p[:, None] + p[None, :])
    Bz, Bc, Bs = B @ zb, B @ (c2 * zb), B @ (s2 * zb)
    del B
    Tz = Az + Bz
    Pv, Qv = Bc - As, Ac - Bs
    res["m"] = m
    res["Im_z_Tz"] = dy * float(np.imag(np.vdot(z, Tz)))
    res["beta_sq"] = dy * float(np.sum(c2 * np.abs(Pv) ** 2 + s2 * np.abs(Qv) ** 2))
    res["e_sq"] = dy * float(np.sum(coth * np.abs(z) ** 2))
    res["R"] = res["beta_sq"] / res["e_sq"]
    res["p_min"] = float(p[0])
    return res


def grid(dy, ytop):
    n = int(round((ytop - YBOT) / dy))
    return YBOT + dy * np.arange(1, n + 1)


def main():
    t0 = time.time()
    C = certified()
    cert = {k: (v if isinstance(v, bool) else {"ball": str(v), "lo": lo(v), "hi": hi(v)}) for k, v in C.items()}
    print("== Part C: certified enclosures (ball arithmetic, prec %d bits) ==" % ctx.prec)
    for k, v in cert.items():
        print(f"  {k:28s} {v if isinstance(v, bool) else v['ball']}")
    m_lo, m_hi, r_lo = lo(C["m_lo"]), hi(C["m_hi"]), lo(C["ratio_lo"])
    print(f"  ENCLOSURE  m(pi) in [{m_lo:.6f}, {m_hi:.6f}]  =>  matrix element mu = -m/2 in [{-m_hi/2:.6f}, {-m_lo/2:.6f}]")
    print(f"  CERTIFIED  ||beta1hat e||/||e||_H >= {r_lo:.6e}   (squared: {r_lo**2:.6e})")
    lines = {}
    lines["C1"] = all(C[f"moment_{k}_contains_exact"] for k in (3, 5, 7))
    lines["C2"] = m_lo > 0
    lines["C3"] = r_lo > 0

    print("== Part N: float64 grid ==", flush=True)
    dy, ytop = 0.02, 40.0
    y = grid(dy, ytop)
    base = model(y, dy, np.exp(-2 * y), kern_inf(y) * dy, checks=True)
    print("  limit  dy=0.02 top=40 :", json.dumps(base), flush=True)
    ex = {"grid_norm_pi0": 0.25, "grid_pi0_W_pi0": 0.125, "grid_pi0_h2_pi0": 0.5, "grid_norm_h2pi0": 17 / 8}
    n1 = {k: abs(base[k] / v - 1) for k, v in ex.items()}
    lines["N1"] = all(v <= 1e-3 for v in n1.values())
    lines["N2"] = m_lo <= base["m"] <= m_hi
    lines["N3"] = abs(base["Im_z_Tz"] + base["m"] / 2) <= 1e-6 * base["m"]
    lines["N4"] = base["R"] > 0 and base["R"] >= r_lo ** 2
    fine = model(grid(0.01, 40.0), 0.01, np.exp(-2 * grid(0.01, 40.0)), kern_inf(grid(0.01, 40.0)) * 0.01)
    print("  limit  dy=0.01 top=40 :", json.dumps(fine), flush=True)
    yb = grid(0.02, 80.0)
    big = model(yb, 0.02, np.exp(-2 * yb), kern_inf(yb) * 0.02)
    print("  limit  dy=0.02 top=80 :", json.dumps(big), flush=True)
    n5 = {"dy": abs(fine["R"] / base["R"] - 1), "box": abs(big["R"] / base["R"] - 1)}
    lines["N5"] = all(v < 0.01 for v in n5.values())
    fin = {}
    for l in LS:
        lam = math.sqrt(l * (l + 1))
        r = model(y, dy, lam ** 2 / np.sinh(y + math.log(2 * lam)) ** 2, kern_l(l, y) * dy)
        r["d"] = abs(r["R"] / base["R"] - 1)
        r["lam_d"] = lam * r["d"]
        fin[l] = r
        print(f"  l={l:5d} (lam W_l, h_l)  :", json.dumps(r), flush=True)
    lines["N6"] = fin[1024]["d"] < fin[256]["d"] < fin[64]["d"] and fin[1024]["d"] < 0.01
    pa = model(y, dy, np.exp(-2 * y), np.diag(np.exp(-y)))
    pb = model(y, dy, np.exp(-2 * y), kern_inf(y) * dy, plant_coth=True)
    print("  plant A (e^-y for W)   :", json.dumps(pa), flush=True)
    print("  plant B (coth dropped) :", json.dumps(pb), flush=True)
    n7 = {"A": abs(pa["R"] / base["R"] - 1), "B": abs(pb["R"] / base["R"] - 1)}
    lines["N7"] = all(v > 0.20 for v in n7.values())

    print("== table ==")
    print(f"  R_inf = ||beta1hat e||^2/||e||^2 = {base['R']:.8f}   m(pi) = {base['m']:.8f}   mu = {base['Im_z_Tz']:.8f}")
    print(f"  N1 rel. errors {n1}")
    print(f"  N5 rel. changes {n5}")
    for l in LS:
        print(f"  l={l:5d}  R_l={fin[l]['R']:.8f}  m_l={fin[l]['m']:.8f}  |R_l/R_inf-1|={fin[l]['d']:.3e}  lam*d={fin[l]['lam_d']:.4f}")
    print(f"  plants: R_A={pa['R']:.8f} (rel {n7['A']:.3f})  R_B={pb['R']:.8f} (rel {n7['B']:.3f})")
    print("== pass lines ==")
    for k, v in lines.items():
        print(f"  {k}: {'PASS' if v else 'FAIL'}")
    ok = all(lines.values())
    print("ALL LINES HELD" if ok else "SOME LINES FAILED", f"({time.time() - t0:.0f} s)")
    if len(sys.argv) > 1:
        with open(sys.argv[1], "w") as f:
            json.dump({"certified": cert, "limit": base, "limit_fine": fine, "limit_bigbox": big,
                       "finite_l": fin, "plantA": pa, "plantB": pb, "n1": {k: float(v) for k, v in n1.items()}, "n5": n5, "n7": n7,
                       "lines": {k: bool(v) for k, v in lines.items()}}, f, indent=1)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
