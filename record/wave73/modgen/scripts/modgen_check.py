"""MODGEN check (pass lines in ..\\logs\\00_passlines.log, fixed before the first run).

Vacuum covariance of the massless scalar on the unit ball, sector l, against the beta = 2 pi thermal form of
h^2 = -d^2/drho^2 + l(l+1)/sinh^2 rho, in both slots, plus the flow-derivative (cross) test.
Usage: C:\\Python314\\python.exe modgen_check.py ..\\out\\modgen_check.json
"""
import sys, json, time
from fractions import Fraction as Fr
from mpmath import mp, mpf, sqrt, pi, gamma, hyp2f1, hyp3f2, coth, quad, sin, cosh, tanh, log, besselj, legendre

OUT = sys.argv[1] if len(sys.argv) > 1 else None
mp.dps = 30
RES = {"lines": [], "plants": [], "tables": {}}
FAILS = []
T0 = time.time()


def line(name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + name + ("  | " + info if info else ""), flush=True)
    RES["lines"].append({"name": name, "ok": bool(ok), "info": info})
    if not ok:
        FAILS.append(name)


def plant(name, caught, info=""):
    print(("CAUGHT " if caught else "MISSED ") + name + ("  | " + info if info else ""), flush=True)
    RES["plants"].append({"name": name, "caught": bool(caught), "info": info})
    if not caught:
        FAILS.append(name)


# ---------------- polynomials in r: dict {exponent: Fraction} ----------------
def padd(a, b, s=1):
    c = dict(a)
    for e, v in b.items():
        c[e] = c.get(e, 0) + s * v
    return {e: v for e, v in c.items() if v != 0}


def pmul(a, b):
    c = {}
    for e1, v1 in a.items():
        for e2, v2 in b.items():
            c[e1 + e2] = c.get(e1 + e2, 0) + v1 * v2
    return {e: v for e, v in c.items() if v != 0}


def pder(a):
    return {e - 1: v * e for e, v in a.items() if e != 0}


def pshift(a, k):
    assert all(e + k >= 0 for e in a), "negative power"
    return {e + k: v for e, v in a.items()}


def pscale(a, s):
    return {e: v * s for e, v in a.items()}


NPOLY = {0: Fr(1, 2), 2: Fr(-1, 2)}            # N = (1 - r^2)/2
ONE_PLUS = {0: Fr(1), 2: Fr(1)}                # 1 + r^2


def lap(a, l):                                 # sector Laplacian: r^n -> (n(n+1) - l(l+1)) r^(n-2)
    out = {}
    for e, v in a.items():
        c = e * (e + 1) - l * (l + 1)
        if c != 0:
            assert e >= 2
            out[e - 2] = out.get(e - 2, 0) + v * c
    return out


def Lop(a, l, weight=True):                    # L f = -div(N grad f) + f ; div(N grad f) = N lap f - r f'
    d = padd(pmul(NPOLY, lap(a, l)), pshift(pder(a), 1), -1)
    out = pscale(d, -1)
    return padd(out, a) if weight else out


def adag(P, m, j):
    """A_j^dagger (N^m P) = N^m [N P' - m r P + j (1+r^2) P/(2r)],  A_j^dagger = d/drho + j coth rho."""
    t = padd(pmul(NPOLY, pder(P)), pscale(pshift(P, 1), m), -1)
    t = padd(t, pscale(pshift(pmul(ONE_PLUS, P), -1), Fr(j, 2)))
    return t


def peval(P, r):
    return sum(mpf(v.numerator) / v.denominator * r ** e for e, v in P.items())


# ---------------- flat side ----------------
def Q_heine(l, t):
    return sqrt(pi) * gamma(l + 1) / gamma(mpf(l) + mpf(3) / 2) * t ** (l + 1) * hyp2f1(mpf(1) / 2, l + 1, mpf(l) + mpf(3) / 2, t * t)


def Q_log(l, t):
    chi = (1 + t * t) / (2 * t)
    W = sum(legendre(k - 1, chi) * legendre(l - k, chi) / k for k in range(1, l + 1))
    return legendre(l, chi) * log((1 + t) / (1 - t)) - W


MCACHE = {}


def M(l, n):
    if (l, n) not in MCACHE:
        MCACHE[(l, n)] = quad(lambda t: t ** n * Q_heine(l, t), [0, mpf(1) / 2, mpf(9) / 10, mpf(99) / 100, 1])
    return MCACHE[(l, n)]


def B(l, a, b):
    """<a, omega^-1 b>_{R^3}/N_Y for polynomial sector profiles supported in the unit ball."""
    s = mpf(0)
    for e1, v1 in a.items():
        for e2, v2 in b.items():
            s += mpf((v1 * v2).numerator) / (v1 * v2).denominator * (M(l, e1 + 1) + M(l, e2 + 1)) / (e1 + e2 + 4)
    return s / pi


# ---------------- thermal side ----------------
def gl_nodes(n):
    xs, ws = [], []
    for i in range(1, n + 1):
        x = mp.cos(mp.pi * (i - mpf(1) / 4) / (n + mpf(1) / 2))
        for _ in range(100):
            p0, p1 = mpf(1), x
            for k in range(2, n + 1):
                p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
            dp = n * (x * p1 - p0) / (x * x - 1)
            dx = p1 / dp
            x -= dx
            if abs(dx) < mpf(10) ** (-mp.dps + 3):
                break
        p0, p1 = mpf(1), x
        for k in range(2, n + 1):
            p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
        dp = n * (x * p1 - p0) / (x * x - 1)
        xs.append(x)
        ws.append(2 / ((1 - x * x) * dp * dp))
    return xs, ws


def composite(panels, n):
    xs, ws = gl_nodes(n)
    out = []
    for a, b in zip(panels[:-1], panels[1:]):
        a, b = mpf(a), mpf(b)
        h, m = (b - a) / 2, (a + b) / 2
        out += [(m + h * x, h * w) for x, w in zip(xs, ws)]
    return out


RHO_PANELS = [mpf(i) for i in range(0, 71)]            # e^(-70) truncation, below 1e-30
K_PANELS = [mpf(i) for i in range(0, 17)]              # transforms decay like e^(-pi k)


def thermal_setup(n):
    rn = composite(RHO_PANELS, n)
    kn = composite(K_PANELS, n)
    S = [[sin(k * rho) for rho, w in rn] for k, wk in kn]
    return rn, kn, S


def data_for(l):
    """returns list of (name, slot, m, P) with g = N^m P(r), and flat polys f_i, p_j."""
    lam2 = l * (l + 1)
    fs, ps, g = [], [], []
    for i in range(3):
        f = pmul({l + 2 * i: Fr(1)}, {0: Fr(1), 2: Fr(-2), 4: Fr(1)})   # r^(l+2i) (1-r^2)^2
        fs.append(f)
        g.append(("q%d" % i, 2, {l + 2 * i + 1: Fr(4)}))                # q = r f = 4 N^2 r^(l+2i+1)
    for j in range(3):
        ps.append({l + 2 * j: Fr(1)})
        g.append(("pi%d" % j, 1, {l + 2 * j + 1: Fr(1)}))               # pi = N r p
    return fs, ps, g, lam2


def thermal(l, grid):
    rn, kn, S = grid
    fs, ps, g, lam2 = data_for(l)
    Nv = [1 / (1 + cosh(rho)) for rho, w in rn]
    rv = [tanh(rho / 2) for rho, w in rn]
    gvals, Gvals, dvals = [], [], []
    for name, m, P in g:
        gvals.append([Nn ** m * peval(P, r) for Nn, r in zip(Nv, rv)])
        # g' = dg/drho = N^m (N P' - m r P)
        dP = padd(pmul(NPOLY, pder(P)), pscale(pshift(P, 1), m), -1)
        dvals.append([Nn ** m * peval(dP, r) for Nn, r in zip(Nv, rv)])
        Pj = P
        for j in range(l, 0, -1):
            Pj = adag(Pj, m, j)
        assert min(Pj) >= 1
        Gvals.append([Nn ** m * peval(Pj, r) for Nn, r in zip(Nv, rv)])
    wr = [w for rho, w in rn]
    hat = []
    for (k, wk), Srow in zip(kn, S):
        nrm = sqrt(2 / pi)
        for j in range(1, l + 1):
            nrm /= sqrt(k * k + j * j)
        sw = [s * w for s, w in zip(Srow, wr)]
        hat.append([nrm * sum(a * b for a, b in zip(sw, Gv)) for Gv in Gvals])

    def kform(wfun, I, J):
        return [[sum(wk * wfun(k) * hv[i] * hv[j] for (k, wk), hv in zip(kn, hat)) for j in J] for i in I]

    def rform(A, Bv, I, J, extra=None):
        return [[sum(w * (A[i][n] * Bv[j][n]) * (extra[n] if extra else 1) for n, w in enumerate(wr)) for j in J] for i in I]

    Iq, Ip, Iall = [0, 1, 2], [3, 4, 5], list(range(6))
    out = {
        "pars": kform(lambda k: mpf(1), Iall, Iall),
        "direct": rform(gvals, gvals, Iall, Iall),
        "en_k": kform(lambda k: k * k, Iq, Iq),
        "en_r": [[a + lam2 * b for a, b in zip(ra, rb)] for ra, rb in zip(
            rform(dvals, dvals, Iq, Iq),
            rform(gvals, gvals, Iq, Iq, extra=[1 / mp.sinh(rho) ** 2 for rho, w in rn]))],
        "pp": kform(lambda k: coth(pi * k) / k, Ip, Ip),
        "qq": kform(lambda k: k * coth(pi * k), Iq, Iq),
        "qp": kform(lambda k: k * coth(pi * k), Iq, Ip),          # <q_i, h coth(pi h) pi_j>
        "pp_wrong": kform(lambda k: coth(pi * k / 2) / k, Ip, Ip),
        "qq_wrong": kform(lambda k: k * coth(pi * k / 2), Iq, Iq),
    }
    return out


def rel(A, Bm):
    sc = max(abs(v) for r in Bm for v in r)
    return max(abs(a - b) for ra, rb in zip(A, Bm) for a, b in zip(ra, rb)) / sc


def s(x, d=3):
    return mp.nstr(x, d)


def tab(Mx, d=17):
    return [[mp.nstr(v, d) for v in r] for r in Mx]


# ================= L0: generator identity (sympy, generic f) =================
import sympy as sp
r_, lam2_ = sp.symbols("r lambda2", positive=True)
f_ = sp.Function("f")(r_)
N_ = (1 - r_ ** 2) / 2
q_ = r_ * f_
lhs = N_ * r_ * (sp.diff(r_ ** 2 * N_ * sp.diff(f_, r_), r_) / r_ ** 2 - N_ * lam2_ * f_ / r_ ** 2 - f_)
rhs = N_ * sp.diff(N_ * sp.diff(q_, r_), r_) - lam2_ * q_ * N_ ** 2 / r_ ** 2     # d/drho = N d/dr, sinh rho = r/N
d0 = sp.simplify(sp.expand(lhs - rhs))
rho_t = mpf("1.37")
chk = max(abs(tanh(rho_t / 2) - mp.sinh(rho_t) / (1 + cosh(rho_t))),
          abs(mp.diff(lambda x: tanh(x / 2), rho_t) - 1 / (1 + cosh(rho_t))),
          abs((1 - tanh(rho_t / 2) ** 2) / 2 - 1 / (1 + cosh(rho_t))))
line("L0 generator identity N r[div(N grad f) - f] = q_rhorho - lambda^2 q/sinh^2 rho (sympy, generic f)", d0 == 0 and chk < mpf("1e-25"),
     f"difference simplifies to {d0}; r = N sinh rho, dr/drho = N, N = 1/(1+cosh rho) hold to {s(chk)}")

# ================= L2 flat guards =================
e2a = max(abs(Q_heine(l, t) - Q_log(l, t)) / abs(Q_log(l, t)) for l in (2, 3) for t in (mpf("0.6"), mpf("0.85"), mpf("0.97")))
e2b = mpf(0)
for l, n in ((2, 5), (3, 8)):
    al = mpf(n + l + 2) / 2
    ser = sqrt(pi) * gamma(l + 1) / gamma(mpf(l) + mpf(3) / 2) / (n + l + 2) * hyp3f2(mpf(1) / 2, l + 1, al, mpf(l) + mpf(3) / 2, al + 1, 1)
    e2b = max(e2b, abs(M(l, n) - ser) / abs(ser))
# Bessel-transform guard, l = 2, f = r^2 (1-r^2)^2
mp.dps = 15
l = 2
fB = pmul({l: Fr(1)}, {0: Fr(1), 2: Fr(-2), 4: Fr(1)})
rB = composite([0, mpf(1) / 4, mpf(1) / 2, mpf(3) / 4, 1], 24)
fBv = [peval(fB, r) * r * r * w for r, w in rB]
kB = composite([mpf(2) * i for i in range(0, 41)], 14)
m3 = m1 = mpf(0)
for k, wk in kB:
    ft = sum(v * sqrt(pi / (2 * k * r)) * besselj(l + mpf(1) / 2, k * r) for (r, w), v in zip(rB, fBv))
    m3 += wk * k ** 3 * ft * ft
    m1 += wk * k * ft * ft
m3, m1 = 2 * m3 / pi, 2 * m1 / pi
mp.dps = 30
p3, p1 = B(2, fB, pscale(lap(fB, 2), -1)), B(2, fB, fB)
e2c = max(abs(m3 - p3) / abs(p3), abs(m1 - p1) / abs(p1))
line("L2 flat guards: (a) Heine Q_l = log form; (b) M(n) quadrature = 3F2(1); (c) Bessel-transform = position-space, l=2",
     e2a <= mpf("1e-20") and e2b <= mpf("1e-10") and e2c <= mpf("1e-4"),
     f"(a) {s(e2a)}; (b) {s(e2b)}; (c) {s(e2c)} [k^3: {s(m3,10)} vs {s(p3,10)}; k^1: {s(m1,10)} vs {s(p1,10)}]")
print(f"   guards done, {time.time()-T0:.0f} s", flush=True)

# ================= thermal grids =================
fine = thermal_setup(30)
print(f"   fine sine table done, {time.time()-T0:.0f} s", flush=True)
coarse = thermal_setup(20)
TH, THC, FLAT = {}, {}, {}
for l in (2, 3):
    TH[l] = thermal(l, fine)
    THC[l] = thermal(l, coarse)
    fs, ps, g, lam2 = data_for(l)
    FLAT[l] = {
        "pp": [[B(l, ps[i], ps[j]) for j in range(3)] for i in range(3)],
        "qq": [[B(l, fs[i], pscale(lap(fs[j], l), -1)) for j in range(3)] for i in range(3)],
        "C1": [[B(l, pscale(lap(fs[i], l), -1), pmul(NPOLY, ps[j])) for j in range(3)] for i in range(3)],
        "C2": [[B(l, Lop(fs[i], l), ps[j]) for j in range(3)] for i in range(3)],
        "C2_noweight": [[B(l, Lop(fs[i], l, weight=False), ps[j]) for j in range(3)] for i in range(3)],
    }
    print(f"   l={l} done, {time.time()-T0:.0f} s", flush=True)

e1 = max(max(rel(TH[l]["pars"], TH[l]["direct"]), rel(TH[l]["en_k"], TH[l]["en_r"])) for l in (2, 3))
line("L1 spectral guard: Parseval (6x6) and int k^2 qhat qhat = <q, h^2 q> (3x3), l = 2, 3", e1 <= mpf("1e-12"), f"max rel err {s(e1)}")

for nm, l, key, lab in (("L3", 2, "pp", "momentum slot <p,omega^-1 p> = <pi,h^-1 coth(pi h) pi>"),
                        ("L4", 3, "pp", "momentum slot <p,omega^-1 p> = <pi,h^-1 coth(pi h) pi>"),
                        ("L5", 2, "qq", "field slot <f,omega f> = <q,h coth(pi h) q>"),
                        ("L6", 3, "qq", "field slot <f,omega f> = <q,h coth(pi h) q>")):
    e = rel(TH[l][key], FLAT[l][key])
    line(f"{nm} l={l} {lab}, 3x3 Gram", e <= mpf("1e-10"),
         f"rel err {s(e)}; flat diag {[s(FLAT[l][key][i][i],15) for i in range(3)]}; thermal diag {[s(TH[l][key][i][i],15) for i in range(3)]}")
    RES["tables"][f"{key}_l{l}"] = {"flat": tab(FLAT[l][key]), "thermal": tab(TH[l][key]), "rel_err": s(e)}

for nm, l in (("L7", 2), ("L8", 3)):
    ea, eb = rel(FLAT[l]["C2"], FLAT[l]["C1"]), rel(TH[l]["qp"], FLAT[l]["C1"])
    line(f"{nm} l={l} cross/flow test: C1 = <N p_j, omega f_i> = C2 = <L f_i, omega^-1 p_j> = C3 = <pi_j, h coth(pi h) q_i>",
         ea <= mpf("1e-10") and eb <= mpf("1e-10"),
         f"C2 vs C1 {s(ea)}; C3 vs C1 {s(eb)}; C1 diag {[s(FLAT[l]['C1'][i][i],15) for i in range(3)]}")
    RES["tables"][f"cross_l{l}"] = {"C1": tab(FLAT[l]["C1"]), "C2": tab(FLAT[l]["C2"]), "C3": tab(TH[l]["qp"]),
                                    "rel_C2_C1": s(ea), "rel_C3_C1": s(eb)}

e9 = max(rel(THC[l][k], TH[l][k]) for l in (2, 3) for k in ("pp", "qq", "qp"))
line("L9 convergence: thermal Grams, 20 against 30 nodes per panel", e9 <= mpf("1e-12"), f"max rel diff {s(e9)}")

w = {(l, k): rel(TH[l][k + "_wrong"], FLAT[l][k]) for l in (2, 3) for k in ("pp", "qq")}
plant("P1 wrong temperature coth(pi h/2), both slots, l = 2, 3", all(v > mpf("1e-3") for v in w.values()),
      "; ".join(f"l={l} {'momentum' if k=='pp' else 'field'} {s(v)}" for (l, k), v in w.items()))
w2 = {l: rel(FLAT[l]["C2_noweight"], FLAT[l]["C1"]) for l in (2, 3)}
plant("P2 weight-1 term dropped from the flow generator", all(v > mpf("1e-3") for v in w2.values()),
      "; ".join(f"l={l} {s(v)}" for l, v in w2.items()))

print("FAILED LINES: " + (", ".join(FAILS) if FAILS else "none"))
print(f"total {time.time()-T0:.0f} s")
RES["failed"] = FAILS
if OUT:
    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump(RES, fh, indent=1)
sys.exit(1 if FAILS else 0)
