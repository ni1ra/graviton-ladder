"""EVEN-BLIND analysis, PREREG.md rules as written.

  python eanalyze.py [nmax=0 (all)] [relcut=1e-3]

Step 1 per l (n >= max(24, 2l), eight largest, exact fit, bracket from four alternates).
Step 2: competing forms, out-of-sample splits S1, S2, S3, gates.
"""
import glob, json, os, statistics, sys
from fractions import Fraction
import numpy as np
from mpmath import mp, mpf, matrix, log, qr_solve

argv = sys.argv[:]
sys.argv = sys.argv[:1]
import continuum as C
mp.dps = 40
HERE = os.path.dirname(os.path.abspath(__file__))
NMAX = int(argv[1]) if len(argv) > 1 else 0
RELCUT = float(argv[2]) if len(argv) > 2 else 1e-3
WEIGHTED = len(argv) > 3 and argv[3] == "1"      # addendum A3
EXTRA = len(argv) > 4 and argv[4] == "1"         # addendum A1: q2only rows (n = 224, 256)


def load():
    out = {k: {} for k in ("q2", "q2_trunc", "q2_site", "p1", "m1", "mean", "p1L")}
    if EXTRA:
        for p in glob.glob(os.path.join(HERE, "results", "raw", "E_*_q2only.json")):
            r = json.load(open(p))
            if r["N"] == 4 * r["n"] and not (NMAX and r["n"] > NMAX) and "ghost" in r["q2"]:
                k = (r["l"], r["n"])
                out["q2"][k] = mpf(r["q2"]["ghost"]); out["q2_trunc"][k] = mpf(r["q2"]["trunc"]); out["q2_site"][k] = mpf(r["q2"]["site"])
    for p in glob.glob(os.path.join(HERE, "results", "raw", "E_*_main.json")):
        r = json.load(open(p))
        if r["N"] != 4 * r["n"] or (NMAX and r["n"] > NMAX):
            continue
        k = (r["l"], r["n"])
        out["q2"][k] = mpf(r["q2"]["ghost"]); out["q2_trunc"][k] = mpf(r["q2"]["trunc"]); out["q2_site"][k] = mpf(r["q2"]["site"])
        d = {Fraction(c["g"]): mpf(c["Q"]) for c in r["cases"] if c["set"] == "ghost" and "Q" in c}
        if 1 in d:
            out["p1"][k] = d[1]; out["p1L"][k] = d[1] * r["l"] / (r["l"] + 1)
        if -1 in d:
            out["m1"][k] = d[-1]
        if 1 in d and -1 in d:
            out["mean"][k] = (d[1] + d[-1]) / 2
    return out


def step1(series, l, st):
    ns = sorted(n for (ll, n) in series if ll == l and n >= max(24, 2 * l))[-8:]
    if len(ns) < 4:
        return None
    qs = [series[(l, n)] for n in ns]
    c = C.exact_fit(ns, qs, "log2", st)
    alts = [C.exact_fit(ns[1:], qs[1:], "log2", st), C.exact_fit(ns[:-1], qs[:-1], "log2", st),
            C.exact_fit(ns, qs, "poly", st), C.exact_fit(ns, qs, "log1", st)]
    return float(c), float(max(abs(c - a) for a in alts)), ns, float(qs[-1])


# ---- forms: name -> (family, list of basis functions of lam) ; power forms handled separately
def B(name):
    u = lambda lam: 1 / lam
    one = lambda lam: mpf(1)
    ul = lambda lam: log(lam) / lam
    u2 = lambda lam: 1 / lam ** 2
    u3 = lambda lam: 1 / lam ** 3
    return {
        "C1": [one, u], "C2": [one, u, u2], "C3": [one, u, u2, u3],
        "CL1": [one, u, ul], "CL2": [one, u, ul, u2], "CL3": [one, u, ul, u2, u3],
        "ZL1": [u, ul], "ZL2": [u, ul, u2], "ZL3": [u, ul, u2, u3],
        "ZG": [lambda lam: 1 / log(lam), lambda lam: 1 / log(lam) ** 2],
        "GL": [one, lambda lam: log(lam), u],
        "XF2": [one, u, u2], "XF3": [one, u, u2, u3],
    }[name]


def offset(name, lam):
    """Imposed term of the exploratory fixed-log forms (addendum 2, X1)."""
    return -log(lam) / lam / 64 if name in ("XF2", "XF3") else mpf(0)


FAM = {"C1": "const", "C2": "const", "C3": "const", "CL1": "const", "CL2": "const", "CL3": "const",
       "ZL1": "zero", "ZL2": "zero", "ZL3": "zero", "ZG": "zero", "GL": "grow", "ZP": "power", "ZP2": "power",
       "XF2": "explore", "XF3": "explore"}
ORDER = ["C1", "C2", "C3", "CL1", "CL2", "CL3", "ZP", "ZP2", "ZL1", "ZL2", "ZL3", "ZG", "GL", "XF2", "XF3"]


def npar(name):
    return {"ZP": 2, "ZP2": 3}.get(name) or len(B(name))


def linfit(name, ls, qs, sg=None):
    bs = B(name)
    A = matrix(len(ls), len(bs)); b = matrix(len(ls), 1)
    for i, l in enumerate(ls):
        lam = mpf(l) + mpf(1) / 2
        w = 1 / mpf(sg[i]) if sg else mpf(1)
        for j, fn in enumerate(bs):
            A[i, j] = w * fn(lam)
        b[i] = w * (mpf(qs[i]) - offset(name, lam))
    sol, res = qr_solve(A, b)
    return (lambda l: float(offset(name, mpf(l) + mpf(1) / 2) + sum(sol[j] * bs[j](mpf(l) + mpf(1) / 2) for j in range(len(bs))))), [float(x) for x in sol]


def powfit(name, ls, qs, sg=None):
    lam = np.array(ls, float) + 0.5
    w = 1 / np.array(sg, float) if sg else np.ones(len(ls))
    q = np.array(qs, float) * w
    best = None
    for s in np.arange(-4.0, 4.0001, 0.001):
        cols = [w * lam ** s] + ([w * lam ** (s - 1)] if name == "ZP2" else [])
        Am = np.stack(cols, 1)
        sol, *_ = np.linalg.lstsq(Am, q, rcond=None)
        r = float(((Am @ sol - q) ** 2).sum())
        if best is None or r < best[0]:
            best = (r, s, sol)
    _, s, sol = best
    return (lambda l: float(sol[0] * (l + 0.5) ** s + (sol[1] * (l + 0.5) ** (s - 1) if name == "ZP2" else 0))), [float(s)] + [float(x) for x in sol]


CONT = {}


def fit(name, ls, qs, s=None):
    if len(ls) < npar(name) + 1:
        return None
    sg = None
    if WEIGHTED and s is not None:
        sg = [max(CONT[(s, l)][1], 1e-4 * abs(CONT[(s, l)][0])) for l in ls]
    return powfit(name, ls, qs, sg) if name in ("ZP", "ZP2") else linfit(name, ls, qs, sg)


def main():
    data = load()
    cont = CONT
    print("lattices nmax=%s relcut=%g weighted=%s q2only-rows=%s" % (NMAX or "all", RELCUT, WEIGHTED, EXTRA))
    for s in ("q2", "q2_trunc", "q2_site", "p1", "m1", "mean", "p1L"):
        st = "site" if s == "q2_site" else "mid"
        print("\nStep 1, series %s:   l   continuum    bracket   rel     largest-n value  lattices" % s)
        for l in sorted({l for (l, n) in data[s]}):
            r = step1(data[s], l, st)
            if r:
                cont[(s, l)] = r
                rel = r[1] / abs(r[0])
                print("   %3d  %.8f  %.1e  %.1e  %.8f  %s%s" % (l, r[0], r[1], rel, r[3], r[2], "" if rel <= RELCUT else "  EXCLUDED"))
    use = {s: [l for l in sorted({l for (ss, l) in cont if ss == s}) if cont[(s, l)][1] <= RELCUT * abs(cont[(s, l)][0])] for s in data}
    print("\nConvention check (Q2 continuum): l  ghost  trunc-ghost  site-ghost  (brackets ghost/trunc/site)")
    for l in sorted({l for (ss, l) in cont if ss == "q2"}):
        if ("q2_trunc", l) in cont and ("q2_site", l) in cont:
            g, t, si = cont[("q2", l)], cont[("q2_trunc", l)], cont[("q2_site", l)]
            print("   %3d  %.8f  %+.1e  %+.1e  (%.1e / %.1e / %.1e)" % (l, g[0], t[0] - g[0], si[0] - g[0], g[1], t[1], si[1]))
    print("\nDoubling differences (q2 continuum): l  D_l = Q_2l - Q_l  ratio D_2l/D_l")
    for l in (2, 4, 8, 16, 3, 6, 12, 24, 5, 10, 20):
        if all(("q2", m) in cont for m in (l, 2 * l)):
            d1 = cont[("q2", 2 * l)][0] - cont[("q2", l)][0]
            rt = ""
            if ("q2", 4 * l) in cont:
                rt = "%.3f" % ((cont[("q2", 4 * l)][0] - cont[("q2", 2 * l)][0]) / d1)
            print("   %3d  %.6f  %s" % (l, d1, rt))
    SPLITS = (("S1", 3, 12, 16), ("S2", 3, 20, 24), ("S3", 4, 16, 20))
    passed = {}
    for s in ("q2", "p1", "mean", "q2_site", "q2_trunc"):
        if len(use[s]) < 5:
            print("\nseries %s: too few usable l (%s)" % (s, use[s])); continue
        print("\nStep 2 out-of-sample, series %s, usable l = %s" % (s, use[s]))
        for sp, lo, hi, tlo in SPLITS:
            tr = [l for l in use[s] if lo <= l <= hi]; te = [l for l in use[s] if l >= tlo]
            print("  split %s train %s test %s%s" % (sp, tr, te, "  NOT EVALUABLE (<2 test points)" if len(te) < 2 else ""))
            if len(te) < 2:
                continue
            preds = {}
            for name in ORDER:
                r = fit(name, tr, [cont[(s, l)][0] for l in tr], s)
                if not r:
                    continue
                pred, par = r
                preds[name] = [pred(l) for l in te]
                ratios = [abs(pred(l) - cont[(s, l)][0]) / cont[(s, l)][1] for l in te]
                rel = [abs(pred(l) - cont[(s, l)][0]) / abs(cont[(s, l)][0]) for l in te]
                ok = max(ratios) <= 2
                fam = FAM[name]
                if fam == "power":
                    fam = "zero" if par[0] < 0 else "grow"
                passed[(s, sp, name)] = (ok, fam)
                print("    %-4s %-5s %s  max miss/bracket %8.2f  max rel miss %.1e  lead par %+.6f  misses/bracket: %s" % (
                    name, fam, "PASS" if ok else "fail", max(ratios), max(rel), par[0], " ".join("%.1f" % x for x in ratios)))
            if preds:
                sp_ = [(max(p[i] for p in preds.values()) - min(p[i] for p in preds.values())) / cont[(s, l)][1] for i, l in enumerate(te)]
                print("    spread of all forms' predictions / bracket per test l: %s" % " ".join("%.1f" % x for x in sp_))
                names = list(preds)
                sep = ["%s|%s" % (a, b) for ia, a in enumerate(names) for b in names[ia + 1:]
                       if max(abs(preds[a][i] - preds[b][i]) / cont[(s, l)][1] for i, l in enumerate(te)) <= 4]
                print("    pairs NOT separated (<= 4 brackets everywhere): %s" % (" ".join(sep) or "none"))
    # gates on q2
    print("\nGATES (deciding series q2)")
    both = {name: passed.get(("q2", "S1", name), (False, ""))[0] and passed.get(("q2", "S2", name), (False, ""))[0] for name in ORDER}
    fam2 = {name: passed.get(("q2", "S2", name), (None, FAM[name]))[1] for name in ORDER}
    pc = [n for n in ORDER if both[n] and fam2[n] == "const"]
    pz = [n for n in ORDER if both[n] and fam2[n] == "zero"]
    pg = [n for n in ORDER if both[n] and fam2[n] == "grow"]
    print("  pass S1 and S2: const %s | zero %s | grow %s" % (pc, pz, pg))
    print("  (S3 passes on q2: %s)" % [n for n in ORDER if passed.get(("q2", "S3", n), (False,))[0]])

    def limits(forms, label, sers=("q2", "p1", "mean")):
        vals = []
        for s in sers:
            for lmin in (3, 4, 6):
                ls = [l for l in use[s] if l >= lmin]
                for name in forms:
                    r = fit(name, ls, [cont[(s, l)][0] for l in ls], s)
                    if r:
                        vals.append((s, lmin, name, r[1][0]))
        if vals:
            v = [x[3] for x in vals]
            print("  %s: %d fits, median %.6f range [%.6f, %.6f] half-range/median %.1f%%" % (
                label, len(v), statistics.median(v), min(v), max(v), 50 * (max(v) - min(v)) / abs(statistics.median(v))))
            for s in sers:
                w = [x[3] for x in vals if x[0] == s]
                if w:
                    print("      series %-4s median %.6f range [%.6f, %.6f]" % (s, statistics.median(w), min(w), max(w)))
            for name in forms:
                w = [x[3] for x in vals if x[2] == name]
                if w:
                    print("      form %-4s median %.6f range [%.6f, %.6f] (%d)" % (name, statistics.median(w), min(w), max(w), len(w)))
        return vals

    limits(pc, "LIMIT from passing constant forms (gate c)")
    print("  for reference, every constant-family form, all windows (spread, not gate):")
    limits(["C1", "C2", "C3", "CL1", "CL2", "CL3"], "all constant forms")
    print("  exploratory (addendum 2): X2 per-series limits with the passing constant forms; X1 fixed-log forms")
    limits(pc, "X2 passing forms on p1L (l^2 (1-A), g=1), p1, m1", ("p1L", "p1", "m1"))
    limits(["XF2", "XF3"], "X1 fixed-log forms (b = -1/64 imposed)", ("q2", "p1", "mean", "p1L"))
    print("  X1 fixed-log passes S1 and S2 on q2: %s" % [n for n in ("XF2", "XF3") if both[n]])
    print("  fitted log coefficients b (CL forms, q2, window l>=3):")
    ls = [l for l in use["q2"] if l >= 3]
    for name in ("CL1", "CL2", "CL3", "ZP", "ZP2", "GL"):
        r = fit(name, ls, [cont[("q2", l)][0] for l in ls], "q2")
        if r:
            print("      %-4s params %s" % (name, " ".join("%+.6f" % x for x in r[1])))


if __name__ == "__main__":
    main()
