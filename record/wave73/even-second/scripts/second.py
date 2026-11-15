"""EVEN-SECOND: blind recomputation, primary route.

Integrals: exact rational closed forms (u is a polynomial) times interval
exponentials; V_16 handled by the geometric-type series
(1-q)^(-2) = sum (n+1) q^n, q = e^(-2y)/1088 <= e^4/1088 on the support,
with an explicit positive tail bound.
Sups: unique critical point (monotone log-derivative) enclosed by a verified
sign change; s_3 additionally needs a cell cover outside a window J.
Run: C:\\Python314\\python.exe second.py
"""
import json, os, sys, time
from fractions import Fraction as Fr
from math import comb
from mpmath import iv, mp, mpf

iv.dps = 50
mp.dps = 50
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "out")
NTERMS = 40
LAM2 = 272
FOURLAM2 = 4 * LAM2  # 1088
YA, YB = -2, 6

# ---------- exact polynomial helpers (coefficients in y, low degree first)
def pmul(a, b):
    r = [Fr(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            r[i + j] += x * y
    return r

def pder(a):
    return [i * a[i] for i in range(1, len(a))] or [Fr(0)]

def peval(a, x):
    r = Fr(0)
    for c in reversed(a):
        r = r * x + c
    return r

def pint(a):
    """exact integral over [YA, YB]"""
    prim = [Fr(0)] + [c / (i + 1) for i, c in enumerate(a)]
    return peval(prim, Fr(YB)) - peval(prim, Fr(YA))

def fiv(q):
    q = Fr(q)
    return iv.mpf(q.numerator) / iv.mpf(q.denominator)

def I(P, k):
    """interval enclosure of int_{YA}^{YB} P(y) e^(-k y) dy, integer k"""
    if k == 0:
        return fiv(pint(P))
    sa = sb = Fr(0)
    D = P
    j = 0
    while any(c != 0 for c in D):
        sa += peval(D, Fr(YA)) / Fr(k) ** (j + 1)
        sb += peval(D, Fr(YB)) / Fr(k) ** (j + 1)
        D = pder(D)
        j += 1
    return iv.exp(iv.mpf(-k * YA)) * fiv(sa) - iv.exp(iv.mpf(-k * YB)) * fiv(sb)

t = [Fr(-1, 2), Fr(1, 4)]                  # (y-2)/4
w = [Fr(1)] + [Fr(0)]
one_minus_t2 = [a - b for a, b in zip([Fr(1), Fr(0), Fr(0)], pmul(t, t))]
u = [Fr(1)]
for _ in range(4):
    u = pmul(u, one_minus_t2)
du, ddu = pder(u), pder(pder(u))
u2, du2, ddu2, uddu = pmul(u, u), pmul(du, du), pmul(ddu, ddu), pmul(u, ddu)
assert peval(u, Fr(2)) == 1 and peval(u, Fr(-2)) == 0 and peval(u, Fr(6)) == 0

qmax = iv.exp(iv.mpf(4)) / FOURLAM2       # sup of q on the support
lam = iv.sqrt(iv.mpf(LAM2))

def tail_T(N):
    """sum_{n>=N} (n+1) q^n <= (N+1) q^N/(1-q)^2 for 0<=q<=qmax"""
    return iv.mpf([0, ((N + 1) * qmax ** N / (1 - qmax) ** 2).b])

def series_Vu2(shift_terms, weight_fn=lambda n: n + 1):
    """sum_n weight(n)/1088^n * sum_{(s,coef)} coef * I(u^2, 2n+2+s) + tail.
    shift_terms: list of (s, coef interval) for an extra factor sum coef e^(-s y)."""
    tot = iv.mpf(0)
    for n in range(NTERMS):
        c = fiv(Fr(weight_fn(n), FOURLAM2 ** n))
        for s, coef in shift_terms:
            tot += c * coef * I(u2, 2 * n + 2 + s)
    lead = sum((coef * I(u2, 2 + s) for s, coef in shift_terms), iv.mpf(0))
    return tot + tail_T(NTERMS) * lead

R = {}
R["norm_u_sq"] = fiv(pint(u2))
R["A_inf"] = fiv(pint(du2)) + I(u2, 2)
R["A_16"] = fiv(pint(du2)) + series_Vu2([(0, iv.mpf(1))])
# planted error: bracket exponent -1 instead of -2, i.e. weights 1 instead of n+1
R["PLANT_A_16"] = fiv(pint(du2)) + series_Vu2([(0, iv.mpf(1))], weight_fn=lambda n: 1)

Ginf_sq = fiv(pint(ddu2)) - 2 * I(uddu, 2) + I(u2, 4)
R["G_inf_sq"] = Ginf_sq
R["G_inf"] = iv.sqrt(Ginf_sq)
# D^2 = int e^(-4y) u^2 ((1-q)^-2 - 1)^2,  ((1-q)^-2-1)^2 = sum_{n>=2} c_n q^n,
# c_n = C(n+3,3) - 2(n+1) >= 0
Dsq = iv.mpf(0)
for n in range(2, NTERMS):
    cn = comb(n + 3, 3) - 2 * (n + 1)
    Dsq += fiv(Fr(cn, FOURLAM2 ** n)) * I(u2, 4 + 2 * n)
ratio = qmax * Fr(NTERMS + 4, NTERMS + 1).numerator / Fr(NTERMS + 4, NTERMS + 1).denominator
assert ratio.b < 1
tailD = comb(NTERMS + 3, 3) * qmax ** NTERMS / (1 - ratio)
Dsq += iv.mpf([0, tailD.b]) * I(u2, 4)
R["D_sq"] = Dsq
R["D"] = iv.sqrt(Dsq)
R["G_hi"] = R["G_inf"] + R["D"]

# Q_16: e^y (1+e^(-y)/(2 lam))^2 = e^y + 1/lam + e^(-y)/1088
Qk = I(du2, -1) + I(du2, 0) / lam + I(du2, 1) / FOURLAM2
Qv = series_Vu2([(-1, iv.mpf(1)), (0, 1 / lam), (1, iv.mpf(1) / FOURLAM2)])
R["Q_16"] = Qk + Qv

# ---------- sups
pi = iv.pi
def sinh(x): return (iv.exp(x) - iv.exp(-x)) / 2
def cosh(x): return (iv.exp(x) + iv.exp(-x)) / 2
def coth(x): return cosh(x) / sinh(x)

def bracket_root(g_mp, g_iv, guess, delta="1e-42"):
    r = mp.findroot(g_mp, guess)
    lo, hi = r - mpf(delta), r + mpf(delta)
    assert g_iv(iv.mpf(lo)).a > 0 and g_iv(iv.mpf(hi)).b < 0, "sign change not verified"
    return iv.mpf([lo, hi])

# S' : f = p e^(-3p^2) sinh(2 pi p); (log f)' = 1/p - 6p + 2pi coth(2pi p), strictly
# decreasing on p>0 (each term is), +inf at 0+, -inf at inf: unique maximiser.
gS_mp = lambda p: 1 / p - 6 * p + 2 * mp.pi * mp.coth(2 * mp.pi * p)
gS_iv = lambda p: 1 / p - 6 * p + 2 * pi * coth(2 * pi * p)
fS = lambda p: p * iv.exp(-3 * p * p) * sinh(2 * pi * p)
pS = bracket_root(gS_mp, gS_iv, 1.2)
R["p_Sprime"] = pS
R["Sprime"] = fS(pS)

# C : f = (p+1/pi) e^(-3p^2); (log f)' = 1/(p+1/pi) - 6p strictly decreasing, >0 at 0.
# root of 6p^2 + 6p/pi - 1 = 0
pC = (-6 / pi + iv.sqrt(36 / pi ** 2 + 24)) / 12
gC_iv = lambda p: 1 / (p + 1 / pi) - 6 * p
assert 0 in gC_iv(pC)
R["p_C"] = pC
R["C"] = (pC + 1 / pi) * iv.exp(-3 * pC * pC)

# s_3 = (sup h)^3, h = p coth(pi p) e^(-p^2).
# k = (log h)' = 1/p - 2pi/sinh(2pi p) - 2p ; k' = -1/p^2 + 4pi^2 cosh(2pi p)/sinh(2pi p)^2 - 2.
# On J=[Ja,Jb]: k' <= -1/Jb^2 + 4pi^2 cosh(2pi Ja)/sinh(2pi Ja)^2 - 2 (cosh/sinh^2 decreasing).
Ja, Jb = iv.mpf("0.6"), iv.mpf("0.75")
kprime_up = -1 / Jb ** 2 + 4 * pi ** 2 * cosh(2 * pi * Ja) / sinh(2 * pi * Ja) ** 2 - 2
assert kprime_up.b < 0, "k not decreasing on J"
k_mp = lambda p: 1 / p - 2 * mp.pi / mp.sinh(2 * mp.pi * p) - 2 * p
k_iv = lambda p: 1 / p - 2 * pi / sinh(2 * pi * p) - 2 * p
ph = bracket_root(k_mp, k_iv, 0.67)
assert Ja.b < ph.a and ph.b < Jb.a
h = lambda p: p * coth(pi * p) * iv.exp(-p * p)
hstar = h(ph)
# cover [0,Ja] and [Jb,3] by cells: h <= b coth(pi b) e^(-a^2) (x coth x increasing)
def cover(a, b, inc, dec, target, depth=0):
    """count of cells proving sup_{[a,b]} inc*dec < target"""
    U = inc(iv.mpf(b)) * dec(iv.mpf(a))
    if U.b < target:
        return 1
    assert depth < 60, "cover failed"
    m = (a + b) / 2
    return cover(a, m, inc, dec, target, depth + 1) + cover(m, b, inc, dec, target, depth + 1)
inc_h = lambda p: p * coth(pi * p)
dec_h = lambda p: iv.exp(-p * p)
iv.dps = 20
ncell = cover(mpf(0), mpf("0.6"), inc_h, dec_h, hstar.a) + cover(mpf("0.75"), mpf(3), inc_h, dec_h, hstar.a)
iv.dps = 50
# tail p>=3: p coth(pi p) <= p + 1/pi, and (p+1/pi) e^(-p^2) is decreasing for p>=1
tail3 = (3 + 1 / pi) * iv.exp(iv.mpf(-9))
assert tail3.b < hstar.a
R["p_s3"] = ph
R["sup_h"] = hstar
R["s_3"] = hstar ** 3
R["s3_cover_cells"] = ncell

# ---------- derived
eps = iv.mpf(3) / 2
e = iv.exp(iv.mpf(1))
Ainf, A16, Ghi, Q16 = R["A_inf"], R["A_16"], R["G_hi"], R["Q_16"]
J0 = Ainf * iv.exp(-eps * Ghi ** 2 / Ainf)
J1 = A16 / (e * eps)
sstar = iv.sqrt(2 * J0 / J1)
M = iv.mpf(8) / 15 * J0 ** 2 * sstar / Q16
Den = iv.sqrt(R["Sprime"] * A16)
Ehi = iv.sqrt(R["C"] * A16 / 2)
c = M / (2 * Den * Ehi)
Elo = iv.sqrt(Ainf * iv.exp(-2 * eps * Ghi ** 2 / Ainf) / (2 * pi))
N3 = iv.sqrt(R["s_3"] * R["norm_u_sq"] / 2)
C0 = (1 + N3 / Elo) * pi ** 2 / 2
R.update(J0_lo=J0, J1_hi=J1, s_star=sstar, M_lo=M, Den_hi=Den, E_hi=Ehi, c=c,
         E_lo=Elo, N3_hi=N3, C0=C0, ratio_c_over_2C0=c / (2 * C0))

SIDE = dict(A_inf="lower", E_lo="lower", c="lower", J0_lo="lower", M_lo="lower",
            s_star="lower", ratio_c_over_2C0="lower",
            A_16="upper", G_hi="upper", Q_16="upper", Sprime="upper", C="upper", s_3="upper",
            N3_hi="upper", C0="upper", J1_hi="upper", Den_hi="upper", E_hi="upper")

def LO(x): return mp.make_mpf(x._mpi_[0])
def HI(x): return mp.make_mpf(x._mpi_[1])

def safe(x, side, digits=15):
    """decimal with `digits` significant digits, rounded outward on the given side"""
    from decimal import Decimal, ROUND_FLOOR, ROUND_CEILING, getcontext
    getcontext().prec = 60
    if side == "lower":
        d = Decimal(mp.nstr(LO(x), 45)) * (1 - Decimal("1e-40"))
        rnd = ROUND_FLOOR
    else:
        d = Decimal(mp.nstr(HI(x), 45)) * (1 + Decimal("1e-40"))
        rnd = ROUND_CEILING
    q = Decimal(10) ** (d.adjusted() - digits + 1)
    return str(d.quantize(q, rounding=rnd))

out = {}
print(f"{'quantity':18s} {'lower (rounded down)':24s} {'upper (rounded up)':24s} width     side")
for name, x in R.items():
    if not hasattr(x, "a"):
        print(f"{name:18s} {x}")
        out[name] = x
        continue
    lo, hi = safe(x, "lower"), safe(x, "upper")
    width = mp.nstr(HI(x) - LO(x), 3)
    print(f"{name:18s} {lo:24s} {hi:24s} {width:9s} {SIDE.get(name, '')}")
    # endpoints pushed outward by 1e-38 relative so the 45-digit decimal strings stay rigorous
    out[name] = dict(lo=mp.nstr(LO(x) - abs(LO(x)) * mpf("1e-38"), 45), hi=mp.nstr(HI(x) + abs(HI(x)) * mpf("1e-38"), 45), lo15=lo, hi15=hi,
                     side=SIDE.get(name, ""))
with open(os.path.join(OUT, "values.json"), "w") as f:
    json.dump(out, f, indent=1)

# ---------- second route for the sups: tolerance cell cover of the whole half-line
# (independent of the critical-point argument): [L, L(1+tol)] encloses the sup.
def tolcover(inc, dec, P, guess, tol=mpf("1e-9")):
    iv.dps = 20
    guess = LO(guess)
    L = LO(inc(iv.mpf(guess)) * dec(iv.mpf(guess)))
    bound = L * (1 + tol)
    stack, n = [(mpf(0), mpf(P))], 0
    while stack:
        a, b = stack.pop()
        if HI(inc(iv.mpf(b)) * dec(iv.mpf(a))) <= bound:
            n += 1
            continue
        assert b - a > mpf("1e-16"), "tolcover stuck"
        m = (a + b) / 2
        stack += [(a, m), (m, b)]
    iv.dps = 50
    return L, bound, n

print("\nSECOND ROUTE FOR SUPS (cell cover of [0,P], tail by monotone decay beyond P):")
t0 = time.time()
checks = [
    ("Sprime", lambda p: p * sinh(2 * pi * p), lambda p: iv.exp(-3 * p * p), 4, R["p_Sprime"], R["Sprime"]),
    ("C", lambda p: p + 1 / pi, lambda p: iv.exp(-3 * p * p), 4, R["p_C"], R["C"]),
    ("sup_h", inc_h, dec_h, 3, R["p_s3"], R["sup_h"]),
]
ok = True
for name, inc, dec, P, guess, prim in checks:
    L, U, n = tolcover(inc, dec, P, guess)
    agree = (L <= HI(prim)) and (LO(prim) <= U)
    ok &= agree
    print(f"  {name:8s} cover enclosure [{mp.nstr(L, 15)}, {mp.nstr(U, 15)}] cells={n} "
          f"contains primary: {agree}")
print(f"  ({time.time() - t0:.1f}s)  all sup checks agree: {ok}")
sys.exit(0 if ok else 1)
