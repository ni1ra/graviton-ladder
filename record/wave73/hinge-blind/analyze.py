"""HINGE-BLIND analysis: tables and the pre-registered extrapolations.

  python analyze.py [results_dir]  > results/tables.txt

Q_l = l(l+1)(1 - A_l)/g^2. All numbers come from results/raw/*.json written by run.py.
"""
import glob, json, math, os, sys
from fractions import Fraction
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "results"), "raw")


def load():
    rows = {}
    for p in sorted(glob.glob(os.path.join(RAW, "*.json"))):
        r = json.load(open(p))
        for c in r["cases"]:
            key = (r["disc"], r["n"], r["N"] // r["n"], r["l"], r["prec"], c["set"], Fraction(c["g"]))
            rows[key] = float(c["Q"]) if "Q" in c else float("nan")
    return rows


ROWS = load()


def Q(disc, n, mult, l, st, g, prec=None):
    """Q for one case; default precision = the lowest one present (the production run)."""
    ks = [k for k in ROWS if k[:4] == (disc, n, mult, l) and k[5] == st and k[6] == Fraction(g) and (prec is None or k[4] == prec)]
    if not ks:
        return None
    return ROWS[min(ks, key=lambda k: k[4])]


def Qmean(disc, n, mult, l, st, g):
    a, b = Q(disc, n, mult, l, st, g), Q(disc, n, mult, l, st, -Fraction(g))
    return None if a is None or b is None else 0.5 * (a + b)


def fmt(x, w=11, d=7):
    return (" " * w) if x is None else ("%*.*f" % (w, d, x))


def table(disc, mult, st, which, ns, ls, title):
    print("\n" + title)
    print("   l \\ n " + "".join("%11d" % n for n in ns))
    for l in ls:
        vals = []
        for n in ns:
            if which == "mean":
                vals.append(Qmean(disc, n, mult, l, st, 1))
            elif which == "mean025":
                vals.append(Qmean(disc, n, mult, l, st, Fraction(1, 4)))
            else:
                vals.append(Q(disc, n, mult, l, st, which))
        if any(v is not None for v in vals):
            print("   %5d " % l + "".join(fmt(v) for v in vals))


# ----------------------------------------------------------------------------- extrapolation in n (PREREG)
def xof(n, st):
    return 1.0 / (n + 0.5) if st != "site" else 1.0 / n


def polyfit_x(ns, qs, st, powers):
    x = np.array([xof(n, st) for n in ns])
    A = np.stack([x ** p for p in powers], axis=1)
    coef = np.linalg.lstsq(A, np.array(qs), rcond=None)[0]
    return coef[0]


def powerfit3(ns, qs, st):
    """Exact three-point fit q0 + c x^p; returns (q0, p) or None."""
    x = [xof(n, st) for n in ns]
    d1, d2 = qs[0] - qs[1], qs[1] - qs[2]
    if d1 == 0 or d2 == 0 or d1 / d2 <= 0:
        return None

    def f(p):
        return (x[0] ** p - x[1] ** p) / (x[1] ** p - x[2] ** p) - d1 / d2
    lo, hi = 0.01, 12.0
    if f(lo) * f(hi) > 0:
        return None
    for _ in range(200):
        m = 0.5 * (lo + hi)
        if f(lo) * f(m) <= 0:
            hi = m
        else:
            lo = m
    p = 0.5 * (lo + hi)
    c = d2 / (x[1] ** p - x[2] ** p)
    return qs[2] - c * x[2] ** p, p


def extrap_n(get, st, nall=(24, 32, 48, 64, 96, 128, 192)):
    """PREREG rule. get(n) -> Q or None. Returns dict(central, bracket, rule, p, alts)."""
    have = [n for n in nall if get(n) is not None]
    q = {n: get(n) for n in have}
    if not all(n in q for n in (48, 96, 192)):
        return None
    d1, d2 = q[48] - q[96], q[96] - q[192]
    p = math.log2(d1 / d2) if d1 / d2 > 0 else float("nan")
    top4 = have[-4:]
    top4d = have[-5:-1]

    def rule1(ns):
        return polyfit_x(ns, [q[n] for n in ns], st, (0, 1, 2))

    def rule2(ns):
        return polyfit_x(ns, [q[n] for n in ns], st, (0, 2, 4))

    def rule3(ns3):
        r = powerfit3(ns3, [q[n] for n in ns3], st)
        return None if r is None else r[0]
    alts = dict(lin_quad=rule1(top4), even=rule2(top4), power3=rule3((48, 96, 192)),
                lin_quad_drop=rule1(top4d), even_drop=rule2(top4d), power3_drop=rule3((32, 64, 128)) if all(n in q for n in (32, 64, 128)) else None)
    if not (p > 0.4):
        return dict(central=None, bracket=None, rule="no limit claimed (p=%.2f)" % p, p=p, alts=alts,
                    slope=(q[192] - q[96]) / math.log(2))
    if abs(p - 1) < 0.3:
        rule, central, others = "lin+quad", alts["lin_quad"], [alts["power3"], alts["lin_quad_drop"]]
    elif abs(p - 2) < 0.5:
        rule, central, others = "even", alts["even"], [alts["power3"], alts["even_drop"]]
    else:
        rule, central = "power3", alts["power3"]
        others = [alts["lin_quad"] if abs(p - 1) < abs(p - 2) else alts["even"], alts["power3_drop"]]
    others = [o for o in others if o is not None]
    br = max(abs(central - o) for o in others) if others else float("nan")
    return dict(central=central, bracket=br, rule=rule, p=p, alts=alts)


# ----------------------------------------------------------------------------- extrapolation in l (PREREG)
def extrap_l(ls, qs, lmin=4):
    ls_, qs_ = zip(*[(l, q) for l, q in zip(ls, qs) if l >= lmin])
    L = np.array([l * (l + 1.0) for l in ls_])
    lv = np.array(ls_, dtype=float)
    y = np.array(qs_)
    out = {}
    for name, cols in (("1/L,1/L^2", [np.ones_like(L), 1 / L, 1 / L ** 2]),
                       ("1/l,1/l^2", [np.ones_like(L), 1 / lv, 1 / lv ** 2]),
                       ("(b+b'ln l)/L", [np.ones_like(L), 1 / L, np.log(lv) / L])):
        A = np.stack(cols, axis=1)
        coef, res, rk, sv = np.linalg.lstsq(A, y, rcond=None)
        out[name] = (coef[0], float(np.sqrt(np.mean((A @ coef - y) ** 2))))
    return out


def main():
    ns = (24, 32, 48, 64, 96, 128, 192)
    ls = sorted({k[3] for k in ROWS if k[0] == "D2" and k[2] == 4})
    for st, nm in (("mid", "primary D2-mid"), ("site", "boundary-site variant D2-site")):
        for which, wn in ((1, "g=+1"), (-1, "g=-1"), ("mean", "mean of g=+-1")):
            table("D2", 4, st, which, ns, ls, "Q_l(n), %s, %s, N_tot=4n" % (nm, wn))
    table("D2", 4, "mid", "mean025", ns, ls, "Q_l(n), primary D2-mid, mean of g=+-1/4, N_tot=4n")
    for st in ("mid", "site"):
        table("Sred", 4, st, "mean", ns, ls, "Q_l(n), Sred-%s, mean of g=+-1, N_tot=4n" % st)

    print("\n\nExtrapolation in n at fixed l (PREREG rule); x = a")
    cont = {}
    for disc, st, which in (("D2", "mid", "mean"), ("D2", "mid", 1), ("D2", "mid", -1), ("D2", "mid", "mean025"),
                            ("D2", "site", "mean"), ("Sred", "mid", "mean")):
        print("\n  %s-%s, %s" % (disc, st, which))
        print("     l    p_eff  rule       central     bracket    | lin+quad     even      power3   | drop-largest: lin+quad   even    power3")
        for l in ls:
            if which == "mean":
                get = lambda n: Qmean(disc, n, 4, l, st, 1)
            elif which == "mean025":
                get = lambda n: Qmean(disc, n, 4, l, st, Fraction(1, 4))
            else:
                get = lambda n: Q(disc, n, 4, l, st, which)
            r = extrap_n(get, st)
            if r is None:
                continue
            a = r["alts"]
            print("   %4d  %6.2f  %-9s %s %s  | %s %s %s | %s %s %s" % (
                l, r["p"], r["rule"][:9], fmt(r["central"]), fmt(r["bracket"], 10, 7),
                fmt(a["lin_quad"], 10), fmt(a["even"], 10), fmt(a["power3"], 10),
                fmt(a["lin_quad_drop"], 10), fmt(a["even_drop"], 10), fmt(a["power3_drop"], 10)))
            cont[(disc, st, which, l)] = r

    print("\n\nExtrapolation in l (PREREG rule) of the continuum values")
    for disc, st, which in (("D2", "mid", "mean"), ("D2", "mid", 1), ("D2", "mid", -1), ("D2", "mid", "mean025"),
                            ("D2", "site", "mean"), ("Sred", "mid", "mean")):
        have = [(l, cont[(disc, st, which, l)]) for l in ls if (disc, st, which, l) in cont and cont[(disc, st, which, l)]["central"] is not None]
        if len(have) < 4:
            continue
        lv = [l for l, _ in have]
        qv = [r["central"] for _, r in have]
        bv = [r["bracket"] for _, r in have]
        print("\n  %s-%s, %s:  l = %s" % (disc, st, which, lv))
        print("    continuum Q_l: " + " ".join("%.6f" % q for q in qv))
        print("    n-brackets   : " + " ".join("%.6f" % b for b in bv))
        qd = dict(zip(lv, qv))
        for l in (2, 4, 8):
            if all(x in qd for x in (l, 2 * l, 4 * l)):
                print("    doubling ratio rho_%d = (Q_%d - Q_%d)/(Q_%d - Q_%d) = %.3f" % (
                    l, 4 * l, 2 * l, 2 * l, l, (qd[4 * l] - qd[2 * l]) / (qd[2 * l] - qd[l])))
        for lmin in (4, 6, 8):
            if len([l for l in lv if l >= lmin]) >= 4:
                fits = extrap_l(lv, qv, lmin)
                print("    l >= %d: " % lmin + "   ".join("%s -> %.6f (rms %.1e)" % (k, v[0], v[1]) for k, v in fits.items()))


if __name__ == "__main__":
    main()
