"""ODD-SECOND second way for the sups: independent branch-and-bound cell cover.

Different library (mpmath interval arithmetic, C:\\Python314) and different argument from the
primary: no derivative, no critical-point bracket.  Every cell [a,b] of [0,P] gets an upper bound
from interval evaluation; cells whose upper bound is below the best proven lower bound (a point
value) are discarded, the rest are bisected.  Tail p >= P handled by an explicit monotone bound.
Result: sup in [L, U].  Writes out/sup_cover.json and logs/02_sup_cover.log (UTF-8).
"""
import json, os, sys, time
from mpmath import iv, mpf, mp

iv.prec = 110
mp.prec = 110
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
LOG = open(os.path.join(ROOT, "logs", "02_sup_cover.log"), "w", encoding="utf-8")
PI = iv.pi
RELGAP = mpf(10) ** (-9) * 2


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    LOG.write(s + "\n")
    LOG.flush()


def lo_(x):
    return mpf(x._mpi_[0])


def hi_(x):
    return mpf(x._mpi_[1])


def sinh(x):
    e = iv.exp(x)
    return (e - 1 / e) / 2


def coth(x):
    e = iv.exp(2 * x)
    return (e + 1) / (e - 1)


def f_S(p):
    return p ** 3 * iv.exp(-3 * p * p) * sinh(2 * PI * p)


def f_C(p):
    return (p ** 3 + p * p / PI) * iv.exp(-3 * p * p)


def f_h(p):
    # p coth(pi p) e^{-p^2}; on cells touching [0, 1/16] use coth x <= 1/x + x/3 (and >= 1/x)
    if lo_(p) < mpf(1) / 16:
        pc = iv.mpf([lo_(1 / PI), hi_(1 / PI + PI * p * p / 3)])
        return pc * iv.exp(-p * p)
    return p * coth(PI * p) * iv.exp(-p * p)


def bnb(name, f, P, tail_ub, cube=False):
    t0 = time.time()
    cells = [(mpf(i) * P / 64, mpf(i + 1) * P / 64) for i in range(64)]
    L = mpf(0)
    neval = 0
    while True:
        ub = []
        for a, b in cells:
            m = (a + b) / 2
            lo = lo_(f(iv.mpf([m, m])))         # rigorous lower bound of a point value
            if lo > L:
                L = lo
            ub.append(hi_(f(iv.mpf([a, b]))))
            neval += 2
        keep = [(c, u_) for c, u_ in zip(cells, ub) if u_ >= L]
        U = max(u_ for _, u_ in keep)
        if (U - L) <= RELGAP * L:
            break
        cells = []
        for (a, b), _ in keep:
            m = (a + b) / 2
            cells.append((a, m))
            cells.append((m, b))
    assert hi_(tail_ub) < L, "tail not excluded"
    pa = min(c[0] for c, _ in keep)
    pb = max(c[1] for c, _ in keep)
    if cube:
        L, U = lo_(iv.mpf([L, L]) ** 3), hi_(iv.mpf([U, U]) ** 3)
    L, U = L * (1 - mpf(10) ** (-22)), U * (1 + mpf(10) ** (-22))   # absorb decimal print rounding
    say(f"{name}: sup in [{mp.nstr(L, 22)}, {mp.nstr(U, 22)}]  maximiser in [{mp.nstr(pa, 12)}, "
        f"{mp.nstr(pb, 12)}]  evals={neval}  {time.time() - t0:.1f}s")
    return {"lo": mp.nstr(L, 25), "hi": mp.nstr(U, 25), "evals": neval}


res = {}
four = iv.mpf(4)
# Sodd tail, p >= 4: f <= p^3 e^{-3p^2 + 2 pi p}/2, whose log-derivative 3/p - 6p + 2 pi < 0 for p >= 4
res["Sodd"] = bnb("Sodd", f_S, mpf(4), four ** 3 * iv.exp(-3 * four * four + 2 * PI * four) / 2)
# Codd tail, p >= 4: log-derivative 2/p + 1/(p+1/pi) - 6p < 0 for p >= 1, so f <= f(4)
res["Codd"] = bnb("Codd", f_C, mpf(4), f_C(four))
# h tail, p >= 4: h <= p coth(4 pi) e^{-p^2} <= 4 coth(4 pi) e^{-16}  (p e^{-p^2} decreasing, p >= 1)
res["s_3"] = bnb("s_3 (=sup h cubed)", f_h, mpf(4), four * coth(4 * PI) * iv.exp(-four * four), cube=True)
with open(os.path.join(ROOT, "out", "sup_cover.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, indent=1)
say("wrote out/sup_cover.json")
