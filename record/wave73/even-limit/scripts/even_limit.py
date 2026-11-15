"""EVEN-LIMIT: lam*W_l -> W_inf on L^2(dy), y = rho - ln(2 lam), x = e^-y.

W_l(rho,rho') = t_<^(l+1) (t_>^-l - t_>^(l+1)) / (2l+1), t = tanh(rho/2), lam = sqrt(l(l+1)).
W_inf(y,y')   = sinh(x_<) e^(-x_>)   (kernel w.r.t. dy dy').
Usage: python even_limit.py <out.json>
"""
import json
import math
import sys

import mpmath as mp
import numpy as np

LS = [16, 64, 256, 1024]
YTOP = 40.0
YBOT = -12.0


# ---------- test functions (smooth, compact support in y) ----------
def f1(y):
    out = np.zeros_like(y)
    m = np.abs(y) < 1
    out[m] = np.exp(-1.0 / (1.0 - y[m] ** 2))
    return out


def f2(y):
    s = (y - 1.0) / 1.5
    out = np.zeros_like(y)
    m = np.abs(s) < 1
    out[m] = np.exp(-1.0 / (1.0 - s[m] ** 2)) * np.cos(3.0 * y[m])
    return out


TESTS = {"f1": (f1, -1.0, 1.0), "f2": (f2, -0.5, 2.5)}


# ---------- kernels ----------
def logt(rho):
    # ln tanh(rho/2) = -2 artanh(e^-rho), accurate for large rho
    return -2.0 * np.arctanh(np.exp(-rho))


def kern_l(l, shift):
    """lam * W_l as a function of (y, y'), rho = y + shift. Zero where rho <= 0."""
    lam = math.sqrt(l * (l + 1))

    def K(y, yp):
        r, rp = y + shift, yp + shift
        ok = (r > 0) & (rp > 0)
        r = np.where(r > 0, r, 1.0)
        rp = np.where(rp > 0, rp, 1.0)
        a, b = logt(r), logt(rp)
        lo, hi = np.minimum(a, b), np.maximum(a, b)      # ln t_<, ln t_>
        # t_<^(l+1) t_>^-l (1 - t_>^(2l+1)) / (2l+1)
        val = np.exp((l + 1) * lo - l * hi) * (-np.expm1((2 * l + 1) * hi)) / (2 * l + 1)
        return np.where(ok, lam * val, 0.0)

    return K


def G(x, xp):
    lo, hi = np.minimum(x, xp), np.maximum(x, xp)
    return 0.5 * np.exp(-(hi - lo)) * (-np.expm1(-2.0 * lo))   # sinh(lo) e^-hi


def kern_inf(y, yp):
    return G(np.exp(-y), np.exp(-yp))


def kern_first(y, yp):
    x, xp = np.exp(-y), np.exp(-yp)
    return -0.5 * (x + xp) * G(x, xp)


# ---------- quadrature ----------
def outer_grid(ylo, yhi, step, n):
    npan = max(1, int(math.ceil((yhi - ylo) / step)))
    edges = np.linspace(ylo, yhi, npan + 1)
    u, w = np.polynomial.legendre.leggauss(n)
    h = np.diff(edges)
    y = (edges[:-1, None] + 0.5 * h[:, None] * (u[None, :] + 1)).ravel()
    wt = (0.5 * h[:, None] * w[None, :]).ravel()
    return y, wt


def apply(K, f, a, b, y, m, n):
    """(K f)(y) = int_a^b K(y,y') f(y') dy', split at the kink y'=y, m GL panels each side."""
    u, w = np.polynomial.legendre.leggauss(n)
    c = np.clip(y, a, b)
    total = np.zeros_like(y)
    for lo, hi in ((np.full_like(y, a), c), (c, np.full_like(y, b))):
        hp = (hi - lo) / m
        for j in range(m):
            left = lo + j * hp
            yp = left[:, None] + 0.5 * hp[:, None] * (u[None, :] + 1)
            wt = 0.5 * hp[:, None] * w[None, :]
            total += np.sum(K(y[:, None], yp) * f(yp) * wt, axis=1)
    return total


def run(res):
    m, n, step, no = res
    out = {}
    for name, (f, a, b) in TESTS.items():
        rows = []
        for l in LS:
            lam = math.sqrt(l * (l + 1))
            row = {"l": l, "lam": lam}
            for tag, shift in (("true", math.log(2 * lam)), ("plant", math.log(lam))):
                ycut = -shift
                y, wt = outer_grid(ycut, YTOP, step, no)
                d = apply(kern_l(l, shift), f, a, b, y, m, n) - apply(kern_inf, f, a, b, y, m, n)
                e2 = np.sum(wt * d * d)
                # below the cut W_l f = 0, so the difference is -W_inf f
                yb, wb = outer_grid(YBOT, ycut, step, no)
                e2 += np.sum(wb * apply(kern_inf, f, a, b, yb, m, n) ** 2)
                row["err_" + tag] = math.sqrt(e2)
                if tag == "true":
                    F = apply(kern_first, f, a, b, y, m, n)
                    row["F_norm"] = math.sqrt(np.sum(wt * F * F))
                    r = lam * d - F
                    row["resid"] = math.sqrt(np.sum(wt * r * r))
            rows.append(row)
        out[name] = rows
    return out


def literal_check():
    """numpy kernel vs the literal formula in mpmath at 60 digits."""
    mp.mp.dps = 60
    rng = np.random.default_rng(73)
    worst = 0.0
    for l in LS:
        lam = math.sqrt(l * (l + 1))
        sh = math.log(2 * lam)
        K = kern_l(l, sh)
        ys = rng.uniform(-sh + 0.05, 3.0, size=(20, 2))
        for y, yp in ys:
            r, rp = mp.mpf(float(y)) + mp.mpf(sh), mp.mpf(float(yp)) + mp.mpf(sh)
            t1, t2 = mp.tanh(min(r, rp) / 2), mp.tanh(max(r, rp) / 2)
            ex = mp.mpf(lam) * t1 ** (l + 1) * (t2 ** (-l) - t2 ** (l + 1)) / (2 * l + 1)
            got = float(K(np.array([y]), np.array([yp]))[0])
            if ex > mp.mpf("1e-250"):
                worst = max(worst, float(abs(got - ex) / ex))
    return worst


def domination_check():
    rng = np.random.default_rng(7)
    worst = -1e300
    neg = 0.0
    for l in LS:
        lam = math.sqrt(l * (l + 1))
        sh = math.log(2 * lam)
        # log-uniform x in (1e-6, 2 lam)
        x = np.exp(rng.uniform(math.log(1e-6), math.log(2 * lam * (1 - 1e-9)), size=(20000, 2)))
        y, yp = -np.log(x[:, 0]), -np.log(x[:, 1])
        k = kern_l(l, sh)(y, yp)
        D = np.minimum(2 * np.minimum(x[:, 0], x[:, 1]), 0.5 * np.exp(-np.abs(x[:, 0] - x[:, 1])))
        worst = max(worst, float(np.max(k - D * (1 + 1e-12))))
        neg = min(neg, float(np.min(k)))
    return worst, neg


def main():
    lines = {}
    lw = literal_check()
    lines["L1"] = lw < 1e-10
    print(f"L1 literal-formula check, worst rel diff = {lw:.2e}")
    dw, neg = domination_check()
    lines["L6"] = (dw <= 0.0) and (neg >= 0.0)
    print(f"L6 domination: max(lamW - D) = {dw:.3e} (must be <= 0), min(lamW) = {neg:.3e}")

    base = run((8, 40, 0.25, 40))
    fine = run((12, 48, 0.125, 48))

    ok2 = ok3 = ok4 = ok5 = ok7 = True
    for name in TESTS:
        print(f"\n{name}:  err = ||lam W_l f - W_inf f||_L2(dy)")
        print("    l      err          l*err     lam*err/||F||  ||lam*d-F||/||F||  err_plant(ln lam)  quad rel chg")
        rows, rf = base[name], fine[name]
        for r, q in zip(rows, rf):
            chg = abs(r["err_true"] - q["err_true"]) / q["err_true"]
            ok5 &= chg < 1e-8
            print(f"{r['l']:5d}  {q['err_true']:.6e}  {r['l']*q['err_true']:.6f}  "
                  f"{q['lam']*q['err_true']/q['F_norm']:.5f}        {q['resid']/q['F_norm']:.5f}"
                  f"            {q['err_plant']:.6e}      {chg:.1e}")
        e = [q["err_true"] for q in rf]
        p = [q["err_plant"] for q in rf]
        ok2 &= all(e[i + 1] < e[i] for i in range(3))
        slope = math.log(e[3] / e[2]) / math.log(LS[3] / LS[2])
        ok3 &= -1.05 <= slope <= -0.95
        last = rf[3]
        ok4 &= abs(last["lam"] * e[3] / last["F_norm"] - 1) < 0.03 and last["resid"] / last["F_norm"] < 0.03
        ok7 &= (p[3] > 0.5 * p[0]) and (p[3] > 100 * e[3])
        print(f"    slope(256->1024) = {slope:.4f};  ||F|| = {last['F_norm']:.6e};  ||W_inf f|| scale: plant err(1024)/err(1024) = {p[3]/e[3]:.1f}")
    lines.update({"L2": ok2, "L3": ok3, "L4": ok4, "L5": ok5, "L7_plant_caught": ok7})
    print()
    for k in sorted(lines):
        print(f"{k}: {'PASS' if lines[k] else 'FAIL'}")
    allok = all(lines.values())
    print("ALL LINES HELD" if allok else "SOME LINES FAILED")
    if len(sys.argv) > 1:
        with open(sys.argv[1], "w") as fh:
            json.dump({"lines": lines, "literal_worst": lw, "base": base, "fine": fine}, fh, indent=1)
    sys.exit(0 if allok else 1)


if __name__ == "__main__":
    main()
