"""Post hoc (NOT pre-registered) continuum and large-l analysis. Kept separate from analyze.py,
which applies the PREREG rules. Motivation: the lattice numbers are exact to ~30 digits, so the
continuum limit at fixed l can be taken by exact-fit (Richardson) extrapolation through several
lattice sizes; the PREREG rule (three-parameter fits) leaves most of that accuracy unused.

Step 1, n -> infinity at fixed l. Basis B = {1, x, x^2 ln x, x^2, x^3 ln x, x^3, x^4, ...} in x = a,
exact fit through the K largest lattices. Central value: K = all lattices with n >= NMIN.
Bracket: largest deviation among (same basis, K-1 points dropping the smallest n), (same basis with
the largest n dropped), (pure polynomial basis, same points), (basis led by x ln x, same points).

Step 2, l -> infinity. Least-squares / exact polynomial fits in 1/(l + c) for c in {0, 1/2, 1},
several orders and windows, on four series that must share one limit: g = +1, g = -1, their mean,
and the second-order coefficient Q2.

  python continuum.py [results_dir] [nmin]
"""
import glob, json, math, os, sys
from fractions import Fraction
from mpmath import mp, mpf, matrix, lu_solve, log, qr_solve

mp.dps = 40
HERE = os.path.dirname(os.path.abspath(__file__))
RES = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "results")
NMIN = int(sys.argv[2]) if len(sys.argv) > 2 else 24


def load(st="mid", disc="D2", mult=4):
    out = {"p1": {}, "m1": {}, "mean": {}, "q2": {}}
    best = {}
    for p in glob.glob(os.path.join(RES, "raw", "%s_*.json" % disc)):
        r = json.load(open(p))
        if r["N"] != mult * r["n"]:
            continue
        k = (r["l"], r["n"])
        d = {Fraction(c["g"]): mpf(c["Q"]) for c in r["cases"] if c["set"] == st and "Q" in c}
        if Fraction(1) in d and Fraction(-1) in d and (k not in best or r["prec"] < best[k]):
            best[k] = r["prec"]
            out["p1"][k], out["m1"][k], out["mean"][k] = d[1], d[-1], (d[1] + d[-1]) / 2
        if "Q2_" + st in r:
            out["q2"][k] = mpf(r["Q2_" + st])
    for p in glob.glob(os.path.join(RES, "pt2", "%s_*.json" % disc)):
        r = json.load(open(p))
        if r["N"] == mult * r["n"]:
            out["q2"][(r["l"], r["n"])] = mpf(r["Q2_" + st])
    return out


def xof(n, st="mid"):
    return mpf(1) / (n + mpf(1) / 2) if st != "site" else mpf(1) / n


BASES = {
    "log2": lambda j, x: [mpf(1), x, x * x * log(x), x * x, x ** 3 * log(x), x ** 3, x ** 4, x ** 4 * log(x), x ** 5, x ** 6][j],
    "poly": lambda j, x: x ** j,
    "log1": lambda j, x: [mpf(1), x * log(x), x, x * x * log(x), x * x, x ** 3 * log(x), x ** 3, x ** 4, x ** 5, x ** 6][j],
}


def exact_fit(ns, qs, basis, st="mid"):
    k = len(ns)
    A = matrix(k, k)
    for i, n in enumerate(ns):
        for j in range(k):
            A[i, j] = BASES[basis](j, xof(n, st))
    return lu_solve(A, matrix(qs))[0]


def continuum(series, l, st="mid", kmax=7):
    ns = sorted(n for (ll, n) in series if ll == l and n >= NMIN)
    if len(ns) < 4:
        return None
    ns = ns[-kmax:]
    qs = [series[(l, n)] for n in ns]
    c = exact_fit(ns, qs, "log2", st)
    alts = [exact_fit(ns[1:], qs[1:], "log2", st), exact_fit(ns[:-1], qs[:-1], "log2", st),
            exact_fit(ns, qs, "poly", st), exact_fit(ns, qs, "log1", st)]
    return float(c), float(max(abs(c - a) for a in alts)), ns, [float(a) for a in alts]


def lfit(ls, qs, c, order, weights=None):
    """Least squares polynomial of given order in u = 1/(l+c); returns constant term."""
    A = matrix(len(ls), order + 1)
    b = matrix(len(ls), 1)
    for i, l in enumerate(ls):
        u = mpf(1) / (l + mpf(c))
        w = mpf(1) if weights is None else mpf(1) / weights[i]
        for j in range(order + 1):
            A[i, j] = w * u ** j
        b[i] = w * mpf(qs[i])
    sol, res = qr_solve(A, b)
    return float(sol[0]), float(res)


def main():
    data = load()
    ls = sorted({l for (l, n) in data["m1"]} | {l for (l, n) in data["q2"]})
    cont = {}
    print("Step 1: continuum values Q_l (n -> infinity), basis log2, n >= %d" % NMIN)
    for name in ("m1", "p1", "mean", "q2"):
        print("\n  series %s" % {"m1": "g = -1", "p1": "g = +1", "mean": "mean of g = +-1", "q2": "Q2 (g -> 0)"}[name])
        print("     l   continuum     bracket    lattices                         alternates: drop-smallest, drop-largest, poly, log1")
        for l in ls:
            r = continuum(data[name], l)
            if r is None:
                continue
            cont[(name, l)] = r
            print("   %4d  %.8f  %.1e   %-32s %s" % (l, r[0], r[1], r[2], " ".join("%.8f" % a for a in r[3])))

    print("\n\nStep 2: l -> infinity. Linear least squares in u = 1/lam, lam = l + 1/2 (lam^2 = l(l+1) + 1/4).")
    print("  forms:  poly-k      : Q_inf + sum_{j<=k} c_j u^j")
    print("          logfree-k   : poly-k + b' u ln(lam)                 (b' fitted; reported)")
    print("          logfix-k    : poly-k - (1/64) u ln(lam)             (far-region coefficient imposed)")
    print("          logfix2-k   : logfix-k + c' u^2 ln(lam)")
    forms = []
    for k in (2, 3, 4):
        forms.append(("poly-%d" % k, k, None, False))
        forms.append(("logfree-%d" % k, k, "free", False))
        forms.append(("logfix-%d" % k, k, "fix", False))
        forms.append(("logfix2-%d" % k, k, "fix", True))
    summary = {}
    for name in ("m1", "p1", "mean", "q2"):
        have = [l for l in ls if (name, l) in cont]
        if len(have) < 5:
            continue
        print("\n  series %s: l = %s" % (name, have))
        for lmin in (2, 3, 4, 6, 8):
            for lmax in sorted({max(have), 32, 24, 16}, reverse=True):
                if lmax > max(have):
                    continue
                sel = [l for l in have if lmin <= l <= lmax]
                if len(sel) < 5:
                    continue
                row = []
                for fname, k, lg, lg2 in forms:
                    npar = k + 1 + (1 if lg == "free" else 0) + (1 if lg2 else 0)
                    if len(sel) < npar + 1:
                        continue
                    A = matrix(len(sel), npar)
                    b = matrix(len(sel), 1)
                    for i, l in enumerate(sel):
                        lam = mpf(l) + mpf(1) / 2
                        u = 1 / lam
                        cols = [u ** j for j in range(k + 1)]
                        if lg == "free":
                            cols.append(u * log(lam))
                        if lg2:
                            cols.append(u * u * log(lam))
                        for j, v in enumerate(cols):
                            A[i, j] = v
                        y = mpf(cont[(name, l)][0])
                        if lg == "fix":
                            y += u * log(lam) / 64
                        b[i] = y
                    sol, res = qr_solve(A, b)
                    extra = " b'=%+.4f" % float(sol[k + 1]) if lg == "free" else ""
                    row.append("%s: %.6f (%.0e)%s" % (fname, float(sol[0]), float(res), extra))
                    summary.setdefault(fname, []).append((name, lmin, lmax, float(sol[0]), float(res)))
                print("    l in [%d,%d] (%d pts): " % (lmin, lmax, len(sel)) + " | ".join(row))


if __name__ == "__main__":
    main()
