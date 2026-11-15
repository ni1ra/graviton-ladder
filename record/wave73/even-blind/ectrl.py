"""EVEN-BLIND controls summary (ii)-(v) + precision, from results/raw JSON."""
import glob, json, os
from fractions import Fraction as F
from mpmath import mp, mpf, log
mp.dps = 40
HERE = os.path.dirname(os.path.abspath(__file__))


def ld(n, N, l, tag):
    ps = glob.glob(os.path.join(HERE, "results", "raw", "E_n%d_N%d_l%d_p*_%s.json" % (n, N, l, tag)))
    return json.load(open(ps[0])) if ps else None


def cs(r, st):
    return {F(c["g"]): c for c in r["cases"] if c["set"] == st}


def Q(c):
    return mpf(c["Q"]) if "Q" in c else None


def s(x, d=9):
    return "none" if x is None else mp.nstr(x, d)


C3 = [(n, l) for n in (24, 48, 96) for l in (2, 8, 32)]
print("== (ii) g -> 0 (ghost): 1-A and Q = l(l+1)(1-A)/g^2 against Q2")
for n, l in C3:
    m, c = ld(n, 4 * n, l, "main"), ld(n, 4 * n, l, "ctrl")
    if not (m and c):
        continue
    g1, gs = cs(m, "ghost"), cs(c, "ghost")
    q2 = mpf(m["q2"]["ghost"])
    print("  n=%d l=%d Q2=%s" % (n, l, s(q2, 12)))
    for g in (F(1), F(1, 4), F(1, 16)):
        a, b = (g1 if g == 1 else gs)[g], (g1 if g == 1 else gs)[-g]
        qa, qb = Q(a), Q(b)
        mean = (qa + qb) / 2 if qa is not None and qb is not None else None
        print("     |g|=%-5s 1-A(+)=%s 1-A(-)=%s  Q(+)=%s Q(-)=%s mean/Q2-1=%s odd part/Q2=%s" % (
            g, a.get("one_minus_A", "ERR")[:12], b.get("one_minus_A", "ERR")[:12], s(qa), s(qb),
            s(mean / q2 - 1, 3) if mean is not None else "none", s((qb - qa) / 2 / q2, 3) if mean is not None else "none"))

print("\n== boundary conventions at finite n, exact g=+1 / g=-1 and Q2 (continuum comparison is in eanalyze step 1)")
for n, l in C3:
    m, c = ld(n, 4 * n, l, "main"), ld(n, 4 * n, l, "ctrl")
    if not (m and c):
        continue
    print("  n=%d l=%d  ghost %s/%s  trunc %s/%s  site %s/%s   Q2 ghost %s trunc %s site %s" % (
        n, l, s(Q(cs(m, "ghost")[1]), 7), s(Q(cs(m, "ghost")[-1]), 7), s(Q(cs(c, "trunc")[1]), 7), s(Q(cs(c, "trunc")[-1]), 7),
        s(Q(cs(c, "site")[1]), 7), s(Q(cs(c, "site")[-1]), 7), m["q2"]["ghost"][:10], m["q2"]["trunc"][:10], m["q2"]["site"][:10]))

print("\n== (iii) local momentum contrast Z = (g/2) Int N pi^2, N=(1-r^2)/2: 1-A (not multiplied by l(l+1))")
rows = [(n, l, ld(n, 4 * n, l, "ctrl")) for n, l in C3] + [(192, l, ld(192, 768, l, "local")) for l in (2, 8, 32)]
for n, l, c in rows:
    if not c:
        continue
    d = cs(c, "local")
    q2 = mpf(c["q2"]["local"]) / (l * (l + 1))
    print("  n=%3d l=%2d  (1-A)/g^2: g=+1 %s  g=-1 %s  g=+1/4 %s  g=-1/4 %s   second-order coeff %s   -lnA(g=1) %s" % (
        n, l, s(Q(d[F(1)]) / (l * (l + 1)), 7) if "Q" in d[F(1)] else "ERR", s(Q(d[F(-1)]) / (l * (l + 1)), 7) if "Q" in d[F(-1)] else "ERR:" + d[F(-1)].get("error", "")[:40],
        s(Q(d[F(1, 4)]) / (l * (l + 1)), 7) if F(1, 4) in d and "Q" in d[F(1, 4)] else "-",
        s(Q(d[F(-1, 4)]) / (l * (l + 1)), 7) if F(-1, 4) in d and "Q" in d[F(-1, 4)] else "-", s(q2, 7),
        s(-log(1 - mpf(d[F(1)]["one_minus_A"])), 7) if "Q" in d[F(1)] else "ERR"))

print("\n== (v) planted errors: ratio to the correct (ghost) value, Q2 and exact g=+1")
for n, l in C3:
    m, c = ld(n, 4 * n, l, "main"), ld(n, 4 * n, l, "ctrl")
    if not (m and c):
        continue
    q2, p1 = mpf(m["q2"]["ghost"]), Q(cs(m, "ghost")[1])
    out = []
    for st in ("neumann", "nocent", "ghost_x2"):
        e = Q(cs(c, st)[F(1)])
        out.append("%s: Q2 ratio %s, g=+1 ratio %s" % (st, s(mpf(c["q2"][st]) / q2, 5) if st in c["q2"] else "(4 by construction)", s(e / p1, 5) if e is not None else "ERR"))
    print("  n=%d l=%d  %s" % (n, l, " | ".join(out)))

print("\n== (iv) outside size: Q(N) - Q(4n), ghost")
for n in (24, 48, 96):
    for l in (2, 4, 8, 32):
        m = ld(n, 4 * n, l, "main")
        if not m:
            continue
        out = []
        for mult in (2, 8):
            o = ld(n, mult * n, l, "out")
            if o:
                out.append("N=%dn: dQ2 %s d(g=+1) %s" % (mult, s(mpf(o["q2"]["ghost"]) - mpf(m["q2"]["ghost"]), 3), s(Q(cs(o, "ghost")[1]) - Q(cs(m, "ghost")[1]), 3)))
        print("  n=%d l=%d  Q2(4n)=%s  %s" % (n, l, m["q2"]["ghost"][:12], " | ".join(out)))

print("\n== precision reruns at 1.4x bits: differences to the default-precision row")
for n, l in ((24, 8), (48, 16), (96, 32)):
    m, p = ld(n, 4 * n, l, "main"), ld(n, 4 * n, l, "prec")
    if m and p:
        print("  n=%d l=%d prec %d -> %d  dQ2=%s d(g=+1)=%s d(g=-1)=%s (30 digits stored)" % (
            n, l, m["prec"], p["prec"], s(mpf(p["q2"]["ghost"]) - mpf(m["q2"]["ghost"]), 3),
            s(Q(cs(p, "ghost")[1]) - Q(cs(m, "ghost")[1]), 3), s(Q(cs(p, "ghost")[-1]) - Q(cs(m, "ghost")[-1]), 3)))
