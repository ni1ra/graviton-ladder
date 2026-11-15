# EVEN-CONST (wave 73): the per-sector affinity constants of the Liouville limit model.
# Pass lines: ..\logs\01_prefixed_pass_lines.log (fixed before the first run).
# Usage: C:\Python314\python.exe even_const.py ..\out\even_const.json
import sys, json, time
import mpmath as mp
import numpy as np

T0 = time.time()
RES = {}
FAILED = []


def line(tag, ok, msg):
    print(("PASS " if ok else "FAIL ") + tag + " | " + msg, flush=True)
    RES[tag] = {"ok": bool(ok), "msg": msg}
    if not ok:
        FAILED.append(tag)


def info(tag, msg):
    print("INFO " + tag + " | " + msg, flush=True)
    RES[tag] = {"ok": None, "msg": msg}


# ---------------------------------------------------------------- framework pieces
def mu_true(p):
    return 2 * p * mp.sinh(mp.pi * p) / mp.pi**2


def mu_plant(p):
    return 2 * p * mp.tanh(mp.pi * p) / mp.pi**2


def I_odd(p, pp):
    return mp.pi**2 / (4 * mp.cosh(mp.pi * (p + pp) / 2) * mp.cosh(mp.pi * (p - pp) / 2))


def J_even(p, pp):
    return (mp.pi**2 * (mp.sech(mp.pi * (p - pp) / 2) - mp.sech(mp.pi * (p + pp) / 2))
            / (8 * p * pp * mp.cosh(mp.pi * p / 2) * mp.cosh(mp.pi * pp / 2)))


def dens_framework(p, pp, even, mu):
    P, q = p + pp, p - pp
    if even:
        M = mp.pi * mp.sqrt(mu(p) * mu(pp)) * J_even(p, pp) * mp.sqrt(p * pp) / 2
    else:
        M = mp.pi * mp.sqrt(mu(p) * mu(pp)) * I_odd(p, pp) / (2 * mp.sqrt(p * pp))
    return (M**2 / (mp.sinh(mp.pi * p) * mp.sinh(mp.pi * pp))
            * ((mp.sinh(mp.pi * P / 2) / P)**2 + (mp.sinh(mp.pi * q / 2) / q)**2) / (2 * mp.pi**2))


def dens_reduced(p, pp, even, mu):
    P, q = p + pp, p - pp
    a, b = mp.pi * P / 2, mp.pi * q / 2
    r, rp = mu(p) / (p * mp.sinh(mp.pi * p)), mu(pp) / (pp * mp.sinh(mp.pi * pp))
    T2 = (mp.tanh(mp.pi * p / 2) * mp.tanh(mp.pi * pp / 2))**2 if even else 1
    return (mp.pi**4 / 128 * r * rp * T2
            * (mp.tanh(a)**2 / (P**2 * mp.cosh(b)**2) + mp.tanh(b)**2 / (q**2 * mp.cosh(a)**2)))


def g_mp(a, b, even, plant=False):
    # reduced integrand in a = pi P/2, b = pi q/2 ; constant = (1/32) int_0^inf da int_0^a db g
    T2 = (mp.tanh((a + b) / 2) * mp.tanh((a - b) / 2))**2 if even else 1
    tb = mp.mpf(1) if b == 0 else (mp.tanh(b) / b)**2
    v = T2 * ((mp.tanh(a) / a)**2 * mp.sech(b)**2 + tb * mp.sech(a)**2)
    if plant:
        v *= mp.sech(a + b) * mp.sech(a - b)
    return v


def setting_A(even):
    def inner(a):
        if a <= 15:
            pts = [0, a]
        elif a <= 30:
            pts = [0, 15, a]
        elif a <= 60:
            pts = [0, 15, 30, a]
        else:
            pts = [0, 15, 30, 60]          # neglects < e^-120
        return mp.quad(lambda b: g_mp(a, b, even), pts)
    return mp.quad(inner, [0, 1, 3, 8, 20, 60, mp.inf]) / 32


def gl_panels(lo, hi, npan, n):
    x, w = np.polynomial.legendre.leggauss(n)
    e = np.linspace(lo, hi, npan + 1)
    mid, half = (e[1:] + e[:-1]) / 2, (e[1:] - e[:-1]) / 2
    return (mid[:, None] + half[:, None] * x).ravel(), (half[:, None] * w).ravel()


def setting_B(even, plant=False, n=10, wpan=0.5):
    b, wb = gl_panels(0.0, 40.0, int(round(40 / wpan)), n)
    s, ws = gl_panels(0.0, 60.0, int(round(60 / wpan)), n)
    out = np.empty_like(b)
    for i0 in range(0, len(b), 100):
        B = b[i0:i0 + 100, None]
        A = B + s[None, :]
        T2 = (np.tanh((A + B) / 2) * np.tanh((A - B) / 2))**2 if even else 1.0
        g = T2 * (np.tanh(A)**2 / (A**2 * np.cosh(B)**2) + np.tanh(B)**2 / (B**2 * np.cosh(A)**2))
        if plant:
            g = g / (np.cosh(A + B) * np.cosh(A - B))
        out[i0:i0 + 100] = g @ ws
    if not plant:
        out += 1.0 / ((b + 60.0) * np.cosh(b)**2)     # exact tail of 1/(a^2 cosh^2 b), a > b+60
    return float(out @ wb) / 32


def xcmsh(x):
    # x cosh x - sinh x without cancellation
    if x < 1:
        s, k = mp.mpf(0), 1
        while True:
            t = x**(2 * k + 1) * 2 * k / mp.factorial(2 * k + 1)
            s += t
            if t < mp.eps * s * mp.mpf(10)**(-5):
                return s
            k += 1
    return x * mp.cosh(x) - mp.sinh(x)


def Phi(x):
    if x == 0:
        return mp.mpf(38) / 3 - 4 * mp.pi
    return (2 - 4 * mp.pi * mp.sech(x) + 8 * x / (mp.sinh(x) * mp.cosh(x))
            + 8 * xcmsh(x) / mp.sinh(x)**3)


def th2(x):
    return mp.mpf(1) if x == 0 else (mp.tanh(x) / x)**2


def beta4():
    return (mp.zeta(4, mp.mpf(1) / 4) - mp.zeta(4, mp.mpf(3) / 4)) / 256


def a_odd_closed():
    return 7 * mp.zeta(3) / (16 * mp.pi**2)


def a_even_closed():
    return 49 * mp.zeta(3) / (16 * mp.pi**2) - 3 * beta4() / mp.pi**2 + mp.catalan / 8 - mp.log(2) / 4


# ---------------------------------------------------------------- Bessel quadratures (y = -ln x)
def Kip(p, x):
    return mp.re(mp.besselk(1j * p, x))


YPTS = [-4.2, -2, 0, 4, 8, 14, 22, 32, 45, 60]


def bessel_int(f):
    # int_0^inf f(x) dx/x  over x in (e^-60, e^4.2)
    return mp.quad(lambda y: f(mp.exp(-y)), YPTS)


# ================================================================= run
mp.mp.dps = 25
PTS3 = [(0.5, 1.2), (1.0, 1.0), (0.3, 2.0)]

# L1
errs = []
for p in (mp.mpf("0.7"), mp.mpf("1.9")):
    num = bessel_int(lambda x: Kip(p, x) * mp.exp(-x) * x)
    errs.append(abs(num - mp.pi * p / mp.sinh(mp.pi * p)))
pa = mp.quad(lambda p: mu_true(p) * (mp.pi * p / mp.sinh(mp.pi * p))**2, [0, 1, 4, 12, 40])
F2 = lambda p: (1 + p**2) * mp.pi * p / (3 * mp.sinh(mp.pi * p))
pc = mp.quad(lambda p: mu_true(p) * F2(p)**2, [0, 1, 4, 12, 40])
p0 = mp.mpf("0.7")
eF2 = abs(bessel_int(lambda x: Kip(p0, x) * mp.exp(-x) * x**2) - F2(p0))
line("L1", max(errs) < 1e-12 and abs(pa - 0.25) < 1e-18 and abs(pc - mp.mpf(3) / 8) < 1e-18 and eF2 < 1e-12,
     "K e^-x transform err %s; Plancherel x e^-x: %s (1/4), err %s; x^2 e^-x: err %s, transform err %s"
     % (mp.nstr(max(errs), 3), mp.nstr(pa, 20), mp.nstr(abs(pa - 0.25), 3), mp.nstr(abs(pc - mp.mpf(3) / 8), 3), mp.nstr(eF2, 3)))

# L2
worst = 0
for (p, pp) in PTS3:
    p, pp = mp.mpf(str(p)), mp.mpf(str(pp))
    num = bessel_int(lambda x: Kip(p, x) * Kip(pp, x) * x)
    worst = max(worst, abs(num / I_odd(p, pp) - 1))
line("L2", worst < 1e-12, "int K_ip K_ip' dx vs pi^2/(4 cosh cosh): worst rel err %s" % mp.nstr(worst, 3))

# L3
mp.mp.dps = 45
worst = 0
for (p, pp) in [(0.5, 1.2), (1.0, 1.3), (0.3, 2.0), (2.5, 0.1), (4.0, 3.0), (0.05, 0.07)]:
    p, pp = mp.mpf(str(p)), mp.mpf(str(pp))
    for even in (False, True):
        for mu in (mu_true, mu_plant):
            worst = max(worst, abs(dens_framework(p, pp, even, mu) / dens_reduced(p, pp, even, mu) - 1))
        a, b = mp.pi * (p + pp) / 2, mp.pi * (p - pp) / 2
        # reduced form in (a,b): density dp dp' = (4/pi^2) da db-measure; g/32 * (pi^2/4)... check constant
        worst = max(worst, abs(dens_reduced(p, pp, even, mu_true) * (4 / mp.pi**2) / (g_mp(a, b, even) / 32) - 1))
        worst = max(worst, abs(dens_reduced(p, pp, even, mu_plant) * (4 / mp.pi**2) / (g_mp(a, b, even, True) / 32) - 1))
line("L3", worst < 1e-25, "framework density = reduced form = g/32 in (a,b): worst rel err %s" % mp.nstr(worst, 3))
mp.mp.dps = 25

# L4
worst = 0
for (p, k) in [(0.7, 0.9), (1.6, 2.5)]:
    p, k = mp.mpf(str(p)), mp.mpf(str(k))
    num = mp.quad(lambda y: Kip(p, mp.exp(-y)) * mp.sin(k * mp.exp(-y)), [-4.2, -3.5, -3, -2.5, -2, -1.5, -1, 0, 2, 5, 10, 20, 40, 70])
    worst = max(worst, abs(num / (mp.pi * mp.sin(p * mp.asinh(k)) / (2 * p * mp.cosh(mp.pi * p / 2))) - 1))
wth = 0
for (p, pp) in PTS3:
    num = mp.quad(lambda t: mp.sin(p * t) * mp.sin(pp * t) / mp.cosh(t), mp.linspace(0, 80, 41))
    wth = max(wth, abs(num - mp.pi / 4 * (mp.sech(mp.pi * (p - pp) / 2) - mp.sech(mp.pi * (p + pp) / 2))))
line("L4", worst < 1e-10 and wth < 1e-15, "sine transform worst rel err %s; theta integral worst abs err %s" % (mp.nstr(worst, 3), mp.nstr(wth, 3)))

# L9, L10 (cheap, before the 2-d work)
mp.mp.dps = 40
worst = 0
for x in ("0.3", "1", "2.5", "6"):
    x = mp.mpf(x)
    Tf = lambda a, b: ((mp.cosh(a) - mp.cosh(b)) / (mp.cosh(a) + mp.cosh(b)))**2
    F1 = mp.quad(lambda b: Tf(x, b) * mp.sech(b)**2, [-x, 0, x])
    F2n = mp.quad(lambda a: Tf(a, x) * mp.sech(a)**2, [x, x + 2, x + 10, x + 40, x + 120])
    worst = max(worst, abs(F1 + 2 * F2n - Phi(x)))
line("L9", worst < 1e-20, "F1 + 2 F2 vs Phi: worst abs err %s" % mp.nstr(worst, 3))

e1 = abs(mp.quad(th2, [0, 1, 5, 20, mp.inf]) - 14 * mp.zeta(3) / mp.pi**2)
e2 = abs(mp.quad(lambda x: mp.tanh(x) / x * mp.sech(x)**2, [0, 1, 5, 20, 60, 120]) - 7 * mp.zeta(3) / mp.pi**2)
e3 = abs(mp.quad(lambda x: (x / mp.sinh(x) - 1) / x**2 if x > mp.mpf("1e-12") else -mp.mpf(1) / 6, [0, 1, 5, 20, mp.inf]) + mp.log(2))
e4 = abs(mp.quad(lambda x: th2(x) * mp.sech(x), [0, 1, 5, 20, 60, 140]) - (48 * beta4() / mp.pi**3 - 2 * mp.catalan / mp.pi))
line("L10", max(e1, e2, e3, e4) < 1e-20, "abs errs: tanh^2/x^2 %s; tanh sech^2/x %s; 1/(x sinh x)-1/x^2 %s; tanh^2 sech/x^2 %s"
     % tuple(mp.nstr(e, 3) for e in (e1, e2, e3, e4)))

# 1-d reduction at high precision
mp.mp.dps = 120
a1d = mp.quad(lambda x: th2(x) * Phi(x), [0, mp.mpf(1) / 2, 1, 2, 5, 12, 30, 80, mp.inf]) / 64
aE = a_even_closed()
aO = a_odd_closed()
d1 = abs(a1d - aE)
print("a_even closed = %s" % mp.nstr(aE, 50))
print("a_even 1-d    = %s" % mp.nstr(a1d, 50))
RES["a_even_closed"] = mp.nstr(aE, 60)
RES["a_even_1d"] = mp.nstr(a1d, 60)
RES["a_odd_closed"] = mp.nstr(aO, 60)

# L12 PSLQ
mp.mp.dps = 100
z3, z5, pi2, l2, G, b4 = mp.zeta(3), mp.zeta(5), mp.pi**2, mp.log(2), mp.catalan, beta4()
B0 = [a1d, 1, z3, z5, pi2, l2, G, 1 / pi2, z3 / pi2, z5 / pi2, l2 / pi2, G / pi2]
r0 = mp.pslq(B0, maxcoeff=10**4, maxsteps=200000)
r1 = mp.pslq(B0 + [b4 / pi2], maxcoeff=10**4, maxsteps=200000)
exp1 = [16, 0, 0, 0, 0, 4, -2, 0, -49, 0, 0, 0, 48]
ok1 = r1 is not None and (list(r1) == exp1 or list(r1) == [-c for c in exp1])
info("L12a", "PSLQ basis B0 (12 elements, |coeff| <= 1e4, 100 digits): %s" % ("NULL (no relation)" if r0 is None else str(r0)))
line("L12b", ok1, "PSLQ basis B0 + beta(4)/pi^2: %s" % str(r1))

# I1
mp.mp.dps = 30
bo = (2 * mp.pi + mp.quad(lambda x: ((1 + mp.tanh(x)**2) * 2 - 4) * (mp.sin(x) / x)**2 if x > 0 else 0, mp.linspace(0, 60, 31))) / 16
be = (2 * mp.pi + mp.quad(lambda x: ((1 + mp.tanh(x)**2) * Phi(x) - 4) * (mp.sin(x) / x)**2 if x > 0 else Phi(0) - 4, mp.linspace(0, 60, 31))) / 16
line("I1", abs(bo - mp.mpf("0.29166")) < 1e-5,
     "||beta_1-hat||^2_HS(tau=pi): odd %s (HINGE J(0.5)/2 = 0.29166), even %s; neither is the affinity constant (%s, %s)"
     % (mp.nstr(bo, 10), mp.nstr(be, 10), mp.nstr(aO, 8), mp.nstr(aE, 8)))
RES["beta_norm_odd"] = mp.nstr(bo, 15)
RES["beta_norm_even"] = mp.nstr(be, 15)

# I2
line("I2", mp.mpf("0.013555") <= aE <= mp.mpf("0.013625"),
     "a_even = %s inside lattice range [0.013555, 0.013625]; a_even - 0.013586 = %s (%s %%)"
     % (mp.nstr(aE, 12), mp.nstr(aE - mp.mpf("0.013586"), 3), mp.nstr((aE / mp.mpf("0.013586") - 1) * 100, 3)))

# L11 kernel cross-check by a finite-difference solve
def Kip_np(p, x, dt=0.05, tmax=40.0):
    t = np.arange(0.0, tmax + dt / 2, dt)
    wt = np.full_like(t, dt); wt[0] = dt / 2
    ct, cw = np.cosh(t), np.cos(p * t) * wt
    out = np.empty_like(x)
    for i0 in range(0, len(x), 2000):
        out[i0:i0 + 2000] = np.exp(-x[i0:i0 + 2000, None] * ct[None, :]) @ cw
    return out


def thomas(lo, di, up, rhs):
    n = len(di)
    c, d = np.empty(n), np.empty(n)
    c[0], d[0] = up[0] / di[0], rhs[0] / di[0]
    for i in range(1, n):
        m = di[i] - lo[i] * c[i - 1]
        c[i] = up[i] / m
        d[i] = (rhs[i] - lo[i] * d[i - 1]) / m
    w = np.empty(n)
    w[-1] = d[-1]
    for i in range(n - 2, -1, -1):
        w[i] = d[i] - c[i] * w[i + 1]
    return w


y0, y1, hf = -np.log(45.0), 30.0, 0.002
nf = int(round((y1 - y0) / hf)); nf += nf % 2
yf = np.linspace(y0, y1, nf + 1); hf = yf[1] - yf[0]
xf = np.exp(-yf)
orders = sorted({p for pr in PTS3 for p in pr})
KG = {p: Kip_np(p, xf) for p in orders}
kchk = max(abs(KG[p][i] - float(Kip(mp.mpf(str(p)), mp.mpf(float(xf[i]))))) for p in orders for i in (10, 2000, 5000, 9000, 15000))


def w_fd(p, pp, step):
    y, h = yf[::step], hf * step
    Kp, Kpp = KG[p][::step], KG[pp][::step]
    yi = y[1:-1]
    n = len(yi)
    lo = np.full(n, -1 / h**2 + 1 / (2 * h)); up = np.full(n, -1 / h**2 - 1 / (2 * h))
    di = 2 / h**2 + np.exp(-2 * yi)
    w = thomas(lo, di, up, np.exp(-yi) * Kpp[1:-1])
    N2 = 2 * np.sqrt(p * pp * np.sinh(np.pi * p) * np.sinh(np.pi * pp)) / np.pi**2
    return N2 * h * float(Kp[1:-1] @ w)


worst, rows = 0, []
for (p, pp) in PTS3:
    wc, wfine = w_fd(p, pp, 2), w_fd(p, pp, 1)
    wr = (4 * wfine - wc) / 3
    v = np.sqrt(p * pp * np.sinh(np.pi * p) * np.sinh(np.pi * pp)) / (np.cosh(np.pi * p) + np.cosh(np.pi * pp))
    wcl = v * np.tanh(np.pi * p / 2) * np.tanh(np.pi * pp / 2) / (p * pp)
    rows.append((p, pp, wr, wcl, wfine))
    worst = max(worst, abs(wr / wcl - 1))
line("L11", worst < 1e-4 and kchk < 1e-10,
     "finite-difference W_inf kernel vs closed form: worst rel err %.2e (Bessel grid vs mpmath %.1e); rows %s"
     % (worst, kchk, "; ".join("(%g,%g) fd %.8f closed %.8f" % r[:4] for r in rows)))

# 2-d integrals
oB, oB16 = setting_B(False), setting_B(False, n=16, wpan=0.4)
eB, eB16 = setting_B(True), setting_B(True, n=16, wpan=0.4)
oPl = setting_B(False, plant=True)
print("setting B done at %.0f s" % (time.time() - T0), flush=True)
mp.mp.dps = 20
tA = time.time()
oA = setting_A(False)
print("setting A odd  = %s  (%.0f s)" % (mp.nstr(oA, 18), time.time() - tA), flush=True)
tA = time.time()
eA = setting_A(True)
print("setting A even = %s  (%.0f s)" % (mp.nstr(eA, 18), time.time() - tA), flush=True)
mp.mp.dps = 30
rOA, rOB = abs(oA / aO - 1), abs(mp.mpf(oB) / aO - 1)
line("L5", rOA < 1e-10 and rOB < 1e-10,
     "odd: A = %s, B = %.15f, closed 7 zeta(3)/(16 pi^2) = %s; rel err A %s, B %s"
     % (mp.nstr(oA, 16), oB, mp.nstr(aO, 16), mp.nstr(rOA, 3), mp.nstr(rOB, 3)))
rPl = abs(mp.mpf(oPl) / aO - 1)
line("L6", rPl > 1e-3, "PLANT mu = p tanh: odd = %.12f, rel deviation %s (caught if > 1e-3)" % (oPl, mp.nstr(rPl, 4)))
dAB = abs(eA - mp.mpf(eB))
line("L7", dAB / aE < 1e-8 and abs(eB16 / eB - 1) < 1e-10 and abs(oB16 / oB - 1) < 1e-10,
     "even: A = %s, B = %.15f, |A-B| = %s (rel %s); B refinement rel change even %.1e, odd %.1e"
     % (mp.nstr(eA, 16), eB, mp.nstr(dAB, 3), mp.nstr(dAB / aE, 3), abs(eB16 / eB - 1), abs(oB16 / oB - 1)))
line("L8", abs(eA / a1d - 1) < 1e-8 and abs(mp.mpf(eB) / a1d - 1) < 1e-8 and d1 < mp.mpf("1e-40"),
     "even 2-d vs 1-d: rel A %s, B %s; |1-d - closed form| = %s"
     % (mp.nstr(abs(eA / a1d - 1), 3), mp.nstr(abs(mp.mpf(eB) / a1d - 1), 3), mp.nstr(d1, 3)))
RES.update({"odd_A": mp.nstr(oA, 18), "odd_B": oB, "even_A": mp.nstr(eA, 18), "even_B": eB, "odd_plant": oPl})

print("FAILED LINES: %s" % (", ".join(FAILED) if FAILED else "none"))
print("total %.0f s" % (time.time() - T0))
if len(sys.argv) > 1:
    with open(sys.argv[1], "w", encoding="utf-8") as f:
        json.dump(RES, f, indent=1)
sys.exit(1 if FAILED else 0)
