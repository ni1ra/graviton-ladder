"""HINGE-BLIND-2: the PREREG-ADDENDUM rule, implemented literally (written before the n = 192 rows
existed on this machine). continuum.py prints the fits; this script applies the rule's selection,
median and range.

  python addendum_rule.py [set=mid] [disc=D2] [prereg | all | nmax] [lmax=32]

Step 1: per l, lattices n >= max(24, 2l), the eight largest; exact fit in
{1, x, x^2 ln x, x^2, x^3 ln x, x^3, x^4, x^4 ln x}; bracket = largest deviation among
drop-smallest, drop-largest, pure polynomial, basis led by x ln x.
Step 2: l with step-1 bracket < 5e-5; series Q2, mean(+-1), g=-1, g=+1; lam = l + 1/2, u = 1/lam;
fits logfix-3, logfix-4, logfix2-3, logfix2-4, logfree-3 on windows l >= 3, 4, 6 (kept when
points >= parameters + 1). Central = median of kept fits, bracket = full range.
Reported beside it (not folded in): poly-2,3,4 and logfree-4 fits, per-series medians.
"""
import sys, statistics
from mpmath import mpf, matrix, log, qr_solve
sys_argv = sys.argv[:]
sys.argv = sys.argv[:1]
import continuum as C
ST = sys_argv[1] if len(sys_argv) > 1 else "mid"
DISC = sys_argv[2] if len(sys_argv) > 2 else "D2"
NSET = sys_argv[3] if len(sys_argv) > 3 else "prereg"   # prereg | all | <nmax>
NMAX = int(NSET) if NSET.isdigit() else 0
LMAX = int(sys_argv[4]) if len(sys_argv) > 4 else 32
data = C.load(ST, DISC)
PRE = (24, 32, 48, 64, 96, 128, 192)
data = {s: {k: v for k, v in d.items() if k[0] <= LMAX and k[0] >= 2 and (NSET == "all" or (NMAX and k[1] <= NMAX and k[1] in PRE) or (NSET == "prereg" and k[1] in PRE))}
        for s, d in data.items()}
SER = ("q2", "mean", "m1", "p1")


def step1(series, l):
    ns = sorted(n for (ll, n) in series if ll == l and n >= max(24, 2 * l))[-8:]
    if len(ns) < 4:
        return None
    qs = [series[(l, n)] for n in ns]
    c = C.exact_fit(ns, qs, "log2", ST)
    alts = [C.exact_fit(ns[1:], qs[1:], "log2", ST), C.exact_fit(ns[:-1], qs[:-1], "log2", ST),
            C.exact_fit(ns, qs, "poly", ST), C.exact_fit(ns, qs, "log1", ST)]
    return float(c), float(max(abs(c - a) for a in alts)), ns


def fit(ls, qs, k, lg, lg2):
    npar = k + 1 + (lg == "free") + lg2
    if len(ls) < npar + 1:
        return None
    A = matrix(len(ls), npar); b = matrix(len(ls), 1)
    for i, l in enumerate(ls):
        lam = mpf(l) + mpf(1) / 2; u = 1 / lam
        cols = [u ** j for j in range(k + 1)]
        if lg == "free":
            cols.append(u * log(lam))
        if lg2:
            cols.append(u * u * log(lam))
        for j, v in enumerate(cols):
            A[i, j] = v
        b[i] = mpf(qs[i]) + (u * log(lam) / 64 if lg == "fix" else 0)
    sol, res = qr_solve(A, b)
    return float(sol[0]), float(res), (float(sol[k + 1]) if lg == "free" else None)


RULE = (("logfix-3", 3, "fix", False), ("logfix-4", 4, "fix", False), ("logfix2-3", 3, "fix", True),
        ("logfix2-4", 4, "fix", True), ("logfree-3", 3, "free", False))
BESIDE = (("poly-2", 2, None, False), ("poly-3", 3, None, False), ("poly-4", 4, None, False),
          ("logfree-2", 2, "free", False), ("logfree-4", 4, "free", False))

print("set=%s disc=%s lattices=%s lmax=%d" % (ST, DISC, NSET, LMAX))
cont = {}
for s in SER:
    print("\nStep 1, series %s:   l  continuum  bracket  lattices" % s)
    for l in sorted({l for (l, n) in data[s]}):
        r = step1(data[s], l)
        if r:
            cont[(s, l)] = r
            print("   %3d  %.8f  %.1e  %s%s" % (l, r[0], r[1], r[2], "" if r[1] < 5e-5 else "   EXCLUDED (bracket >= 5e-5)"))
kept, beside, per = [], [], {s: [] for s in SER}
print("\nStep 2 (rule fits):")
for s in SER:
    ls_all = [l for l in sorted({l for (ss, l) in cont if ss == s}) if cont[(s, l)][1] < 5e-5]
    for lmin in (3, 4, 6):
        ls = [l for l in ls_all if l >= lmin]
        row = []
        for name, k, lg, lg2 in RULE:
            r = fit(ls, [cont[(s, l)][0] for l in ls], k, lg, lg2)
            if r:
                kept.append(r[0]); per[s].append(r[0])
                if name == 'logfree-3':
                    beside.append((name, r[0]))
                row.append("%s %.6f%s" % (name, r[0], "" if r[2] is None else " (b'=%+.5f)" % r[2]))
        print("  %-4s l>=%d (%d pts, l=%s): %s" % (s, lmin, len(ls), ls, " | ".join(row)))
        row = []
        for name, k, lg, lg2 in BESIDE:
            r = fit(ls, [cont[(s, l)][0] for l in ls], k, lg, lg2)
            if r:
                beside.append((name, r[0]))
                row.append("%s %.6f%s" % (name, r[0], "" if r[2] is None else " (b'=%+.5f)" % r[2]))
        print("        beside: %s" % " | ".join(row))
if kept:
    print("\nADDENDUM RULE: %d fits kept. median = %.6f   range = [%.6f, %.6f]   half-range = %.6f" % (
        len(kept), statistics.median(kept), min(kept), max(kept), (max(kept) - min(kept)) / 2))
    for s in SER:
        if per[s]:
            print("   series %-4s: median %.6f  range [%.6f, %.6f]  (%d fits)" % (s, statistics.median(per[s]), min(per[s]), max(per[s]), len(per[s])))
    for nm in sorted({n for n, _ in beside}):
        v = [x for n, x in beside if n == nm]
        print("   beside %-9s: median %.6f  range [%.6f, %.6f]  (%d fits)" % (nm, statistics.median(v), min(v), max(v), len(v)))
    fl = [x for n, x in beside if n in ("logfree-2", "logfree-3")]
    if fl:
        print("FREE-LOG ESTIMATE (ADDENDUM-2 item 4): %d fits. median = %.6f   range = [%.6f, %.6f]" % (len(fl), statistics.median(fl), min(fl), max(fl)))
