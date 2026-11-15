"""Convergence studies and planted-error controls, tabulated from results/raw.

  python controls.py [results_dir] > results/controls.txt
"""
import glob, json, os, sys
from fractions import Fraction
from mpmath import mp, mpf

mp.dps = 34
HERE = os.path.dirname(os.path.abspath(__file__))
RES = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "results")

runs = [json.load(open(p)) for p in sorted(glob.glob(os.path.join(RES, "raw", "*.json")))]


def get(disc, n, mult, l, st, g, prec=None, pick=min):
    c = []
    for r in runs:
        if (r["disc"], r["n"], r["N"], r["l"]) != (disc, n, mult * n, l):
            continue
        if prec is not None and r["prec"] != prec:
            continue
        for cs in r["cases"]:
            if cs["set"] == st and Fraction(cs["g"]) == Fraction(g) and "Q" in cs:
                c.append((r["prec"], mpf(cs["Q"])))
    if not c:
        return None
    return pick(c)[1]


def precs(disc, n, mult, l):
    return sorted({r["prec"] for r in runs if (r["disc"], r["n"], r["N"], r["l"]) == (disc, n, mult * n, l)})


def f(x, d=9):
    return "      -     " if x is None else "%.*f" % (d, float(x))


print("A. Precision: the same case at the production precision and at about 1.4x the bits (primary set, g=+1 and g=-1)")
print("   n   l   bits(prod) bits(high)   |Q_high - Q_prod| (g=+1)   (g=-1)     digits of margin at production")
for n in (24, 48, 96):
    for l in (2, 8, 32):
        ps = precs("D2", n, 4, l)
        if len(ps) < 2:
            continue
        a1, b1 = get("D2", n, 4, l, "mid", 1, ps[0]), get("D2", n, 4, l, "mid", 1, ps[-1])
        a2, b2 = get("D2", n, 4, l, "mid", -1, ps[0]), get("D2", n, 4, l, "mid", -1, ps[-1])
        marg = [r["margin_digits"] for r in runs if (r["disc"], r["n"], r["N"], r["l"], r["prec"]) == ("D2", n, 4 * n, l, ps[0])][0]
        print("  %3d %3d   %6d     %6d        %.1e                 %.1e      %.0f" % (n, l, ps[0], ps[-1], float(abs(a1 - b1)), float(abs(a2 - b2)), marg))

print("\nB. Outside size: Q(mean of g=+-1, primary set) against N_tot/n")
print("   n   l      N=2n          N=4n          N=8n          N=16n        (8n - 4n)     (16n - 8n)")
for n in (24, 48, 96):
    for l in (2, 4, 8, 32):
        v = {}
        for m in (2, 4, 8, 16):
            a, b = get("D2", n, m, l, "mid", 1), get("D2", n, m, l, "mid", -1)
            v[m] = None if a is None or b is None else (a + b) / 2
        d1 = None if v[8] is None or v[4] is None else v[8] - v[4]
        d2 = None if v[16] is None or v[8] is None else v[16] - v[8]
        print("  %3d %3d  %s  %s  %s  %s   %s   %s" % (n, l, f(v[2], 10), f(v[4], 10), f(v[8], 10), f(v[16], 10),
                                                    "    -    " if d1 is None else "%+.1e" % float(d1), "    -    " if d2 is None else "%+.1e" % float(d2)))

print("\nC. Radial discretisation: Srednicki form against the direct second difference (primary boundary, mean of g=+-1)")
print("   n   l      D2            Sred          Sred - D2")
for n in (24, 48, 96, 128, 192):
    for l in (2, 4, 8, 16, 32):
        a = [get(d, n, 4, l, "mid", g) for d in ("D2", "Sred") for g in (1, -1)]
        if None in a:
            continue
        d2, sr = (a[0] + a[1]) / 2, (a[2] + a[3]) / 2
        print("  %3d %3d  %s  %s  %+.2e" % (n, l, f(d2, 10), f(sr, 10), float(sr - d2)))

print("\nD. Boundary site of the phi^2 integral (mean of g=+-1): primary 'mid' (a=1/(n+1/2), w_n=1) against")
print("   'site' (a=1/n, w_n=1/2) and the planted errors w_n=0, w_n=2 (a=1/(n+1/2)); ratio to primary in brackets")
print("   n   l      mid           site                 w_n=0                w_n=2")
for n in (24, 48, 96):
    for l in (2, 8, 32):
        row = []
        for st in ("mid", "site", "w0", "w2"):
            a, b = get("D2", n, 4, l, st, 1), get("D2", n, 4, l, st, -1)
            row.append(None if a is None or b is None else (a + b) / 2)
        if row[0] is None:
            continue
        s = "  %3d %3d  %s" % (n, l, f(row[0], 8))
        for v in row[1:]:
            s += "   " + ("        -         " if v is None else "%s (%.3f)" % (f(v, 8), float(v / row[0])))
        print(s)

print("\nE. Sign and normalisation of the perturbation (primary set): g=+1, g=-1, and Z with the factor 1/2 dropped")
print("   (reported as l(l+1)(1-A) for nominal g=+-1, i.e. the number a wrong-normalisation run would have quoted)")
print("   n   l      g=+1          g=-1        (g=-1)/(g=+1)   no-1/2, g=+1   no-1/2, g=-1   ratio of means")
for n in (24, 48, 96):
    for l in (2, 8, 32):
        a, b = get("D2", n, 4, l, "mid", 1), get("D2", n, 4, l, "mid", -1)
        c, d = get("D2", n, 4, l, "half", 1), get("D2", n, 4, l, "half", -1)
        if None in (a, b, c, d):
            continue
        print("  %3d %3d  %s  %s    %.4f       %s    %s     %.4f" % (n, l, f(a, 8), f(b, 8), float(b / a), f(c, 8), f(d, 8), float((c + d) / (a + b))))

print("\nF. Scaling in g (primary set): mean over +-g of l(l+1)(1-A)/g^2, and the odd part [Q(-g)-Q(+g)]/(2g)")
print("   n   l    Q2 (2nd order)   |g|=1/16      |g|=1/4       |g|=1/2       |g|=1        mean(1/4)/mean(1)   odd/g: 1/16     1/4       1")
for n in (24, 48, 96):
    for l in (2, 8, 32):
        q2 = [r.get("Q2_mid") for r in runs if (r["disc"], r["n"], r["N"], r["l"]) == ("D2", n, 4 * n, l) and r.get("Q2_mid")]
        if not q2:
            for p in glob.glob(os.path.join(RES, "pt2", "D2_n%d_N%d_l%d_*.json" % (n, 4 * n, l))):
                q2 = [json.load(open(p))["Q2_mid"]]
        m, o = {}, {}
        for g in (Fraction(1, 16), Fraction(1, 4), Fraction(1, 2), Fraction(1)):
            a, b = get("D2", n, 4, l, "mid", g), get("D2", n, 4, l, "mid", -g)
            if a is None or b is None:
                continue
            m[g], o[g] = (a + b) / 2, (b - a) / 2 / mpf(g.numerator) * g.denominator
        if Fraction(1) not in m or Fraction(1, 4) not in m:
            continue
        print("  %3d %3d   %s   %s  %s  %s  %s      %.5f        %s  %s  %s" % (
            n, l, f(mpf(q2[0]), 9) if q2 else "     -     ", f(m.get(Fraction(1, 16)), 9), f(m.get(Fraction(1, 4)), 9), f(m.get(Fraction(1, 2)), 9),
            f(m.get(Fraction(1)), 9), float(m[Fraction(1, 4)] / m[Fraction(1)]),
            f(o.get(Fraction(1, 16)), 6), f(o.get(Fraction(1, 4)), 6), f(o.get(Fraction(1)), 6)))
