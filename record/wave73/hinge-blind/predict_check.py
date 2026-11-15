"""HINGE-BLIND-2: out-of-sample test of the large-l forms. Fit on l <= LFIT (step-1 continuum values,
all lattices), predict the continuum value at larger l, compare with the step-1 value and bracket
there (which the fit never saw).   python predict_check.py [LFIT=16]
"""
import sys
a = sys.argv[:]; sys.argv = a[:1] + ["mid", "D2", "all", "64"]
import addendum_rule as R
from mpmath import mpf, matrix, log, qr_solve
LFIT = int(a[1]) if len(a) > 1 else 16
def coef(ls, qs, k, lg, lg2):
    A = matrix(len(ls), k + 1 + (lg == "free") + lg2); b = matrix(len(ls), 1)
    for i, l in enumerate(ls):
        lam = mpf(l) + mpf(1) / 2; u = 1 / lam
        cols = [u ** j for j in range(k + 1)] + ([u * log(lam)] if lg == "free" else []) + ([u * u * log(lam)] if lg2 else [])
        for j, v in enumerate(cols): A[i, j] = v
        b[i] = mpf(qs[i]) + (u * log(lam) / 64 if lg == "fix" else 0)
    return qr_solve(A, b)[0]
def ev(c, l, k, lg, lg2):
    lam = mpf(l) + mpf(1) / 2; u = 1 / lam
    cols = [u ** j for j in range(k + 1)] + ([u * log(lam)] if lg == "free" else []) + ([u * u * log(lam)] if lg2 else [])
    return float(sum(c[j] * v for j, v in enumerate(cols)) - (u * log(lam) / 64 if lg == "fix" else 0))
print("\nOUT-OF-SAMPLE: fit on 3 <= l <= %d, predict larger l (step-1 value +- bracket in brackets)" % LFIT)
for s in R.SER:
    ls = [l for l in sorted({l for (ss, l) in R.cont if ss == s}) if 3 <= l <= LFIT]
    tg = [l for l in sorted({l for (ss, l) in R.cont if ss == s}) if l > LFIT]
    print(" series %s, fit l=%s" % (s, ls))
    for name, k, lg, lg2 in (("logfix-3", 3, "fix", False), ("logfree-2", 2, "free", False), ("logfree-3", 3, "free", False), ("poly-3", 3, None, False), ("poly-4", 4, None, False)):
        if len(ls) < k + 2 + (lg == "free"): continue
        c = coef(ls, [R.cont[(s, l)][0] for l in ls], k, lg, lg2)
        print("   %-9s Q_inf=%.6f | " % (name, float(c[0])) + "  ".join("l=%d: pred %.6f, miss %+.1e (meas %.6f +- %.0e)" % (
            l, ev(c, l, k, lg, lg2), ev(c, l, k, lg, lg2) - R.cont[(s, l)][0], R.cont[(s, l)][0], R.cont[(s, l)][1]) for l in tg))
