"""Post hoc (not pre-registered) continuum extrapolation diagnostics: exact-fit polynomials in
x = a through the k largest lattices, for k = 2..7, and the effective exponent sequence.
The data are exact to ~30 digits, so the only question is the functional form in x.

  python richardson.py [results_dir] [set] [which]      which: mean | p1 | m1 | q2
"""
import glob, json, math, os, sys
from fractions import Fraction
from mpmath import mp, mpf, matrix, lu_solve, log

mp.dps = 40
HERE = os.path.dirname(os.path.abspath(__file__))
RES = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "results")
ST = sys.argv[2] if len(sys.argv) > 2 else "mid"
WHICH = sys.argv[3] if len(sys.argv) > 3 else "mean"
DISC = sys.argv[4] if len(sys.argv) > 4 else "D2"


def load():
    out = {}
    if WHICH == "q2":
        for p in glob.glob(os.path.join(RES, "pt2", "%s_*.json" % DISC)):
            r = json.load(open(p))
            if r["N"] == 4 * r["n"]:
                out[(r["l"], r["n"])] = mpf(r["Q2_" + ST])
        return out
    best = {}
    for p in glob.glob(os.path.join(RES, "raw", "%s_*.json" % DISC)):
        r = json.load(open(p))
        if r["N"] != 4 * r["n"]:
            continue
        d = {Fraction(c["g"]): mpf(c["Q"]) for c in r["cases"] if c["set"] == ST and "Q" in c}
        if Fraction(1) not in d or Fraction(-1) not in d:
            continue
        v = {"mean": (d[1] + d[-1]) / 2, "p1": d[1], "m1": d[-1]}[WHICH]
        k = (r["l"], r["n"])
        if k not in best or r["prec"] < best[k][0]:
            best[k] = (r["prec"], v)
    return {k: v[1] for k, v in best.items()}


def x_of(n):
    return mpf(1) / (n + mpf(1) / 2) if ST != "site" else mpf(1) / n


def exact_fit(ns, qs, basis):
    """Solve for coefficients of the basis functions through the points; return the constant."""
    k = len(ns)
    A = matrix(k, k)
    for i, n in enumerate(ns):
        x = x_of(n)
        for j in range(k):
            A[i, j] = basis(j, x)
    c = lu_solve(A, matrix(qs))
    return c[0]


def main():
    data = load()
    ls = sorted({l for l, n in data})
    print("set=%s which=%s disc=%s" % (ST, WHICH, DISC))
    for l in ls:
        ns = sorted(n for (ll, n) in data if ll == l)
        if len(ns) < 3:
            continue
        qs = [data[(l, n)] for n in ns]
        print("\n l=%d  n=%s" % (l, ns))
        print("   raw      : " + " ".join("%.7f" % q for q in qs))
        # local exponent from consecutive triples
        ps = []
        for i in range(len(ns) - 2):
            x0, x1, x2 = [x_of(n) for n in ns[i:i + 3]]
            d1, d2 = qs[i] - qs[i + 1], qs[i + 1] - qs[i + 2]
            if d1 / d2 > 0:
                # solve (x0^p - x1^p)/(x1^p - x2^p) = d1/d2 by bisection
                f = lambda p: (x0 ** p - x1 ** p) / (x1 ** p - x2 ** p) - d1 / d2
                lo, hi = mpf("0.01"), mpf(12)
                if f(lo) * f(hi) < 0:
                    for _ in range(80):
                        m = (lo + hi) / 2
                        if f(lo) * f(m) <= 0:
                            hi = m
                        else:
                            lo = m
                    ps.append(float((lo + hi) / 2))
                else:
                    ps.append(float("nan"))
            else:
                ps.append(float("nan"))
        print("   local p  : " + " ".join("%.3f" % p for p in ps))
        for name, basis in (("poly x^j        ", lambda j, x: x ** j),
                            ("1, x, x^2 ln x..", lambda j, x: [mpf(1), x, x * x * log(x), x * x, x ** 3 * log(x), x ** 3, x ** 4][j]),
                            ("1, x ln x, x, ..", lambda j, x: [mpf(1), x * log(x), x, x * x * log(x), x * x, x ** 3 * log(x), x ** 3][j])):
            row = []
            for k in range(2, len(ns) + 1):
                row.append(float(exact_fit(ns[-k:], qs[-k:], basis)))
            print("   %s top k=2..%d: " % (name, len(ns)) + " ".join("%.7f" % v for v in row))
            if len(ns) >= 4:
                row = []
                for k in range(2, len(ns)):
                    row.append(float(exact_fit(ns[-k - 1:-1], qs[-k - 1:-1], basis)))
                print("   %s   (largest n dropped): " % name[:14] + " ".join("%.7f" % v for v in row))


if __name__ == "__main__":
    main()
