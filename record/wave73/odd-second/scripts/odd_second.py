"""ODD-SECOND primary: rigorous enclosures (python-flint / Arb ball arithmetic).

Run with the venv interpreter that has python-flint.  Writes out/primary.json and
logs/01_primary.log (UTF-8).  Everything is computed from the definitions in the brief.
"""
import json, os, sys, time
from flint import arb, acb, fmpq, fmpq_poly, fmpz, ctx

ctx.prec = 700
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
LOG = open(os.path.join(ROOT, "logs", "01_primary.log"), "w", encoding="utf-8")


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    LOG.write(s + "\n")
    LOG.flush()


D = 30  # decimal digits printed


def lohi(x, d=D):
    """Decimal strings lo <= x <= hi, outward rounded (one extra unit of slack)."""
    t = x * arb(10) ** d
    lo = (t.mid() - t.rad()).floor().mid().floor().unique_fmpz() - 1
    hi = (t.mid() + t.rad()).ceil().mid().ceil().unique_fmpz() + 1

    def fmt(n):
        n = int(n)
        sgn = "-" if n < 0 else ""
        s = str(abs(n)).rjust(d + 1, "0")
        return sgn + s[:-d] + "." + s[-d:]
    return fmt(lo), fmt(hi)


# ---------------------------------------------------------------- profile
y = fmpq_poly([0, 1])
t = (y - 2) / 4
u = (1 - t * t) ** 4          # on [-2, 6]; zero of order 4 at both ends
du = u.derivative()
ddu = du.derivative()
YA, YB = fmpq(-2), fmpq(6)


def poly_int(P):
    Q = P.integral()
    return Q(YB) - Q(YA)          # exact fmpq


def exp_int(P, a):
    """Exact-form integral of P(y) e^{a y} over [-2,6], a a nonzero integer.
    = [ e^{ay} sum_k (-1)^k P^(k)(y) / a^(k+1) ]."""
    a = fmpq(a)
    sa = fmpq(0)
    sb = fmpq(0)
    Pk = P
    k = 0
    while Pk.degree() >= 0:
        c = fmpq((-1) ** k) / a ** (k + 1)
        sa += c * Pk(YA)
        sb += c * Pk(YB)
        Pk = Pk.derivative()
        k += 1
    return (arb(a * YB).exp() * arb(sb)) - (arb(a * YA).exp() * arb(sa))


def P_at(P, x):
    """Evaluate fmpq_poly at an acb by Horner."""
    r = acb(0)
    for c in reversed(P.coeffs()):
        r = r * x + acb(arb(c))
    return r


def quad(f):
    """Rigorous Gauss-Legendre quadrature with Arb error bounds (acb.integral)."""
    old = ctx.prec
    ctx.prec = 200
    try:
        v = acb.integral(lambda x, _: f(x), acb(-2), acb(6), rel_tol=arb(2) ** (-150),
                         abs_tol=arb(2) ** (-170), eval_limit=10 ** 7, depth_limit=10 ** 4)
    finally:
        ctx.prec = old
    assert v.imag.contains(0)
    return v.real


LAM2 = 272
FOURL2 = 4 * LAM2            # 1088
u2 = u * u

# ---------------------------------------------------------------- closed forms
T0 = time.time()
I_du2 = poly_int(du * du)             # exact
I_ddu2 = poly_int(ddu * ddu)          # exact
norm2 = poly_int(u2)                  # exact
say("exact rationals: int u'^2 =", I_du2, "  int u''^2 =", I_ddu2, "  ||u||^2 =", norm2)

E = {}                                 # cache of int u^2 e^{a y}


def Eu2(a):
    if a not in E:
        E[a] = exp_int(u2, a)
    return E[a]


A_inf = arb(I_du2) + Eu2(-2)

# V_16 = e^{-2y} (1-x)^{-2}, x = e^{-2y}/1088, 0 < x <= r := e^4/1088 on the support.
r = arb(4).exp() / FOURL2
assert r < arb("0.051")
N = 90
S = arb(0)
for n in range(N):
    S += arb(fmpq(n + 1, FOURL2 ** n)) * Eu2(-2 * (n + 1))
# tail: sum_{n>=N} (n+1) x^n <= (N+1) r^N / (1-q), q = r (N+2)/(N+1)  (term ratio bound)
q = r * fmpq(N + 2, N + 1)
tailA = (N + 1) * r ** N / (1 - q) * Eu2(-2)
tailA = arb(0).union(tailA)
A_16 = arb(I_du2) + S + tailA
say("A_16 series: N =", N, " tail bound <=", tailA.upper().str(5))

# G_inf^2 = int (u'')^2 - 2 u'' u e^{-2y} + u^2 e^{-4y}
Ginf2 = arb(I_ddu2) - 2 * exp_int(ddu * u, -2) + Eu2(-4)
G_inf = Ginf2.sqrt()

# D^2 = || (V16 - Vinf) u ||^2 ;  ((1-x)^-2 - 1)^2 = sum_n d_n x^n, d_n = C(n+3,3) - 2(n+1), n>=1
S = arb(0)
for n in range(1, N):
    dn = (n + 3) * (n + 2) * (n + 1) // 6 - 2 * (n + 1)
    assert dn >= 0
    S += arb(fmpq(dn, FOURL2 ** n)) * Eu2(-2 * (n + 2))
cN = (N + 3) * (N + 2) * (N + 1) // 6
q = r * fmpq(N + 4, N + 1)
tailD = cN * r ** N / (1 - q) * Eu2(-4)
tailD = arb(0).union(tailD)
D2 = S + tailD
Dn = D2.sqrt()
G_hi = G_inf + Dn

lam = arb(LAM2).sqrt()
Qodd_inf = Eu2(1)
Qodd_16 = Eu2(1) + arb(norm2) / lam + Eu2(-1) / (4 * LAM2)

# ---------------------------------------------------------------- quadrature (second way)
say("\n-- self-check 1: closed form vs rigorous quadrature (acb.integral) --")


def V16c(x):
    e = (-2 * x).exp()
    return e / ((1 - e / FOURL2) ** 2)


checks = []


def check(name, closed, quadv, expect_fail=False):
    ok = bool(closed.overlaps(quadv))
    diff = (closed - quadv)
    rel = (diff.abs_upper() / closed.abs_lower())
    verdict = "AGREE" if ok else "MISMATCH"
    say(f"  {name:12s} closed={closed.str(25, radius=False)}  quad={quadv.str(25)}  "
        f"|diff|<={rel.str(3, radius=False)} rel  -> {verdict}")
    checks.append({"name": name, "agree": ok, "rel_diff_upper": float(rel), "planted": expect_fail})
    return ok


q_Ainf = quad(lambda x: P_at(du, x) ** 2 + (-2 * x).exp() * P_at(u, x) ** 2)
q_A16 = quad(lambda x: P_at(du, x) ** 2 + V16c(x) * P_at(u, x) ** 2)
q_Ginf2 = quad(lambda x: (-P_at(ddu, x) + (-2 * x).exp() * P_at(u, x)) ** 2)
q_D2 = quad(lambda x: ((V16c(x) - (-2 * x).exp()) * P_at(u, x)) ** 2)
q_Q16 = quad(lambda x: x.exp() * (1 + (-x).exp() / (2 * acb(lam))) ** 2 * P_at(u, x) ** 2)
q_Qinf = quad(lambda x: x.exp() * P_at(u, x) ** 2)
q_n2 = quad(lambda x: P_at(u, x) ** 2)

allok = True
allok &= check("A_inf", A_inf, q_Ainf)
allok &= check("A_16", A_16, q_A16)
allok &= check("G_inf^2", Ginf2, q_Ginf2)
allok &= check("D^2", D2, q_D2)
allok &= check("Qodd_16", Qodd_16, q_Q16)
allok &= check("Qodd_inf", Qodd_inf, q_Qinf)
allok &= check("||u||^2", arb(norm2), q_n2)
say("  all integral checks agree:", allok)

say("\n-- planted error: closed form of Qodd_16 built with weight e^(2y) instead of e^y --")
# e^{2y}(1+e^{-y}/(2 lam))^2 = e^{2y} + e^{y}/lam + 1/(4 lam^2)
Q_plant = Eu2(2) + Eu2(1) / lam + arb(norm2) / (4 * LAM2)
plant_caught = not check("Qodd_16 PLANT", Q_plant, q_Q16, expect_fail=True)
say("  plant caught by self-check:", plant_caught)

# ---------------------------------------------------------------- sups (primary: log-derivative)
say("\n-- sups: monotone log-derivative, bracket of the unique critical point --")
PI = arb.pi()
W = arb(10) ** (-60)


def bisect(g, a, b, it=260):
    a, b = arb(a), arb(b)
    for _ in range(it):
        m = ((a + b) / 2).mid()
        gm = g(m)
        if gm > 0:
            a = m
        elif gm < 0:
            b = m
        else:
            break
    return ((a + b) / 2).mid()


def sup_by_bracket(name, f, g, a0, b0):
    p0 = bisect(g, a0, b0)
    a, b = p0 - W, p0 + W
    assert g(a) > 0 and g(b) < 0, "bracket sign check failed"
    X = a.union(b)
    lo = f(p0)                # sup >= f(p0)
    hi = f(X)                 # sup = f(p*) with p* in X, so sup in f(X)
    val = lo.union(hi)
    say(f"  {name}: maximiser in [{lohi(a, 40)[0]}, {lohi(b, 40)[1]}]")
    say(f"     value in [{lohi(val)[0]}, {lohi(val)[1]}]")
    return val, p0, X


# Sodd: f = p^3 e^{-3p^2} sinh(2 pi p);  (log f)' = 3/p - 6p + 2 pi coth(2 pi p): each term strictly
# decreasing on (0,inf) => unique zero, f increasing before and decreasing after => global max.
fS = lambda p: p ** 3 * (-3 * p * p).exp() * (2 * PI * p).sinh()
gS = lambda p: 3 / p - 6 * p + 2 * PI * (2 * PI * p).coth()
Sodd, pS, XS = sup_by_bracket("Sodd", fS, gS, "0.5", "3")

# Codd: f = (p^3 + p^2/pi) e^{-3p^2}; (log f)' = 2/p + 1/(p + 1/pi) - 6p, strictly decreasing.
fC = lambda p: (p ** 3 + p * p / PI) * (-3 * p * p).exp()
gC = lambda p: 2 / p + 1 / (p + 1 / PI) - 6 * p
Codd, pC, XC = sup_by_bracket("Codd", fC, gC, "0.2", "3")

# s_3 = (sup h)^3, h = p coth(pi p) e^{-p^2}, h(0) = 1/pi.
# (log h)' = 1/p - 2 pi / sinh(2 pi p) - 2p is NOT monotone on (0,inf) (it starts at 0 and rises),
# so: (i) on R = [0.5, 0.9] prove (log h)'' < 0 by interval subdivision => unique crit. point in R;
#     (ii) outside R exclude by cells: h < h(p0) on [0,0.05] (analytic bound), [0.05,0.5], [0.9,3];
#     (iii) p >= 3: h <= p coth(3 pi) e^{-p^2} <= 3 coth(3 pi) e^{-9} (p e^{-p^2} decreasing for p>=1).
fh = lambda p: p * (PI * p).coth() * (-p * p).exp()
gh = lambda p: 1 / p - 2 * PI / (2 * PI * p).sinh() - 2 * p
dgh = lambda p: -1 / (p * p) + 4 * PI * PI * (2 * PI * p).cosh() / (2 * PI * p).sinh() ** 2 - 2
hval, ph, Xh = sup_by_bracket("sup h (s_3 = h^3)", fh, gh, "0.5", "0.9")
ctx.prec = 80
ncell = 0
for i in range(400):                     # (i) concavity of log h on R
    c = arb(fmpq(500 + i, 1000)).union(arb(fmpq(501 + i, 1000)))
    assert dgh(c) < 0
    ncell += 1
Lh = arb(fh(ph).lower())
assert arb(0.5) < Xh and Xh < arb(0.9)
# (ii) [0, 0.05]: coth x <= 1/x + x/3  => p coth(pi p) <= 1/pi + pi p^2/3, e^{-p^2} <= 1
assert 1 / PI + PI * arb("0.05") ** 2 / 3 < Lh


def cover(a, b, n):
    global ncell
    for i in range(n):
        lo = a + (b - a) * i / n
        hi = a + (b - a) * (i + 1) / n
        c = arb(lo).union(arb(hi))
        assert fh(c) < Lh, (lo, hi)
        ncell += 1


cover(fmpq(1, 20), fmpq(1, 2), 4500)
cover(fmpq(9, 10), fmpq(3), 4200)
assert 3 * (3 * PI).coth() * arb(-9).exp() < Lh      # (iii)
ctx.prec = 700
say(f"  s_3 exclusion cover: {ncell} cells verified, concavity of log h on [0.5,0.9] verified")
s3 = hval ** 3

# ---------------------------------------------------------------- derived
eps = arb(3) / 2
J0_lo = A_inf * (-eps * G_hi ** 2 / A_inf).exp()
J1_hi = A_16 / (arb(1).exp() * eps)
sstar = (2 * J0_lo / J1_hi).sqrt()
M_lo = arb(8) / 15 * J0_lo ** 2 * sstar / Qodd_16
Den_hi = (Sodd * A_16).sqrt()
E_hi = (Codd * A_16 / 2).sqrt()
c_odd = M_lo / (2 * Den_hi * E_hi)
unorm = arb(norm2).sqrt()
E_lo = J0_lo / ((2 * PI).sqrt() * unorm)
N3_hi = (s3 * A_16 / 2).sqrt()
K = (1 + N3_hi / E_lo) / 2
C0 = K * PI * PI
ratio = c_odd / (2 * C0)
l1 = 1 / ratio

rows = [
    ("||u||^2", arb(norm2), "exact", "exact rational 262144/109395 (polynomial integral)"),
    ("||u||", unorm, "exact", "sqrt of exact rational"),
    ("int u'^2", arb(I_du2), "exact", "exact rational"),
    ("A_inf", A_inf, "LOWER", "closed form (int by parts, exact rationals x e^{4}, e^{-12})"),
    ("A_16", A_16, "UPPER", "geometric-type series in e^{-2y}/1088, 90 closed-form terms + tail"),
    ("G_inf", G_inf, "UPPER (via G_hi)", "sqrt of closed form"),
    ("||(V16-Vinf)u||", Dn, "UPPER (via G_hi)", "sqrt of series of closed forms + tail"),
    ("G_hi", G_hi, "UPPER", "sum of the two"),
    ("Qodd_inf", Qodd_inf, "(info)", "closed form"),
    ("Qodd_16", Qodd_16, "UPPER", "closed form, lambda = sqrt(272)"),
    ("Sodd", Sodd, "UPPER", "monotone log-derivative + 1e-60 bracket"),
    ("Codd", Codd, "UPPER", "monotone log-derivative + 1e-60 bracket"),
    ("s_3", s3, "UPPER", "concavity on [0.5,0.9] + bracket + exclusion cover + tail"),
    ("J0_lo", J0_lo, "LOWER", "derived"),
    ("J1_hi", J1_hi, "UPPER", "derived"),
    ("s*", sstar, "LOWER", "derived"),
    ("M_lo", M_lo, "LOWER", "derived"),
    ("Den_hi", Den_hi, "UPPER", "derived"),
    ("E_hi", E_hi, "UPPER", "derived"),
    ("c_odd", c_odd, "LOWER", "derived"),
    ("E_lo", E_lo, "LOWER", "derived"),
    ("N3_hi", N3_hi, "UPPER", "derived"),
    ("K", K, "UPPER", "derived"),
    ("C0", C0, "UPPER", "derived"),
    ("c_odd/(2 C0)", ratio, "LOWER", "derived"),
    ("l1", l1, "UPPER", "derived"),
]
say("\n-- results (outward-rounded decimal enclosures) --")
out = {"quantities": {}, "checks": checks, "integral_checks_all_agree": bool(allok),
       "plant_caught": bool(plant_caught),
       "maximisers": {"Sodd": lohi(XS, 40), "Codd": lohi(XC, 40), "s_3": lohi(Xh, 40)}}
for name, v, side, method in rows:
    lo, hi = lohi(v)
    relw = float((v.rad() * 2) / v.abs_lower())
    say(f"  {name:16s} [{lo}, {hi}]  side={side}  relwidth~{relw:.1e}")
    out["quantities"][name] = {"lo": lo, "hi": hi, "side": side, "method": method, "rel_width": relw}
with open(os.path.join(ROOT, "out", "primary.json"), "w", encoding="utf-8") as fh_:
    json.dump(out, fh_, indent=1)
say(f"\nelapsed {time.time() - T0:.1f} s; wrote out/primary.json")
sys.exit(0 if (allok and plant_caught) else 1)
