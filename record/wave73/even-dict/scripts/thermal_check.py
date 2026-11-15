"""EVEN-DICT secondary check (lines L9, L10, plant P4 of ..\\logs\\00_passlines.log): momentum slot of the thermal norm.

flat vacuum:  <p_i, omega^-1 p_j>_{R^3}/N_Y = int int p_i(r) p_j(r') r r' Q_l(chi)/pi dr dr', chi = (r^2+r'^2)/(2 r r')
thermal:      <pi_i, h^-1 coth(pi h) pi_j>_{L^2(drho)} by the spectral resolution of h^2, pi = N r p.
Usage: python thermal_check.py ..\\out\\thermal_check.json
"""
import sys, json, time
from mpmath import mp, mpf, mpc, cosh, sinh, exp, log, sqrt, pi, gamma, hyp2f1, coth, quad

OUT = sys.argv[1] if len(sys.argv) > 1 else None
L = 2
KS = (0, 1)                      # p = r^(L+2k)
RES = {"lines": [], "plants": []}
FAILS = []


def line(name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + name + ("  | " + info if info else ""), flush=True)
    RES["lines"].append({"name": name, "ok": bool(ok), "info": info})
    if not ok:
        FAILS.append(name)


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


def unorm(k):
    """delta(k-k') normalisation: sinh^(l+1) 2F1(a,b;c;-sinh^2) -> 2|C| cos(k rho + phase)."""
    b = mpc(L + 1, -k) / 2
    cma = mpc(L + 2, -k) / 2
    C = gamma(mpf(L) + mpf(3) / 2) * gamma(mpc(0, -k)) / (gamma(b) * gamma(cma))
    return sqrt(2 / pi) / (2 * abs(C))


def uraw(k, rho):
    a = mpc(L + 1, k) / 2
    b = mpc(L + 1, -k) / 2
    s = sinh(rho)
    return (s ** (L + 1) * hyp2f1(a, b, mpf(L) + mpf(3) / 2, -s * s)).real


def pi_of(kk, rho):
    r = 1 - 2 / (exp(rho) + 1)
    N = 1 / (1 + cosh(rho))
    return N * r * r ** (L + 2 * kk)


t0 = time.time()
# ---- L10a: ODE guard at 30 digits ----
mp.dps = 30
worst = mpf(0)
for k, rho in ((mpf("0.7"), mpf("1.3")), (mpf(3), mpf(6)), (mpf("5.5"), mpf("0.4"))):
    h = mpf("1e-6")
    u0, up, um = uraw(k, rho), uraw(k, rho + h), uraw(k, rho - h)
    res = -(up - 2 * u0 + um) / h ** 2 + L * (L + 1) * u0 / sinh(rho) ** 2 - k * k * u0
    worst = max(worst, abs(res) / abs(u0))
okode = worst <= mpf("1e-8")

mp.dps = 20
# ---- thermal side ----
rnodes = composite("0 0.5 1 1.5 2 3 4 5 6 8 10 12 14 17 20 23 26 30 34 38 42 46".split(), 40)
knodes = composite("0 0.25 0.5 1 1.5 2 3 4 5 6 7 8.5".split(), 16)
pis = [[pi_of(kk, rho) for rho, w in rnodes] for kk in KS]
hat = []
for k, wk in knodes:
    nk = unorm(k)
    u = [nk * uraw(k, rho) for rho, w in rnodes]
    hat.append([sum(w * pv * uv for (rho, w), pv, uv in zip(rnodes, pis[i], u)) for i in range(len(KS))])
print(f"   transforms done, {time.time()-t0:.0f} s", flush=True)
n = len(KS)
therm = [[sum(wk * coth(pi * k) / k * hv[i] * hv[j] for (k, wk), hv in zip(knodes, hat)) for j in range(n)] for i in range(n)]
wrongT = [[sum(wk * coth(pi * k / 2) / k * hv[i] * hv[j] for (k, wk), hv in zip(knodes, hat)) for j in range(n)] for i in range(n)]
pars = [[sum(wk * hv[i] * hv[j] for (k, wk), hv in zip(knodes, hat)) for j in range(n)] for i in range(n)]
direct = [[sum(w * pis[i][m] * pis[j][m] for m, (rho, w) in enumerate(rnodes)) for j in range(n)] for i in range(n)]


# ---- flat side ----
def Q2(chi, r, rp):
    # run 02 crashed here on chi - 1 rounding to 0; (chi+1)/(chi-1) = ((r+rp)/(r-rp))^2 exactly
    if chi > 3:
        return sqrt(pi) * gamma(L + 1) / (gamma(mpf(L) + mpf(3) / 2) * (2 * chi) ** (L + 1)) * \
            hyp2f1(mpf(L + 1) / 2, mpf(L + 2) / 2, mpf(L) + mpf(3) / 2, 1 / chi ** 2)
    return (3 * chi * chi - 1) / 4 * 2 * log((r + rp) / (r - rp)) - 3 * chi / 2


def flat(i, j):
    ei, ej = L + 2 * KS[i] + 1, L + 2 * KS[j] + 1          # p(r) r

    def outer(r):
        def inner(rp):
            if rp <= 0 or rp >= r:
                return mpf(0)
            q = Q2((r * r + rp * rp) / (2 * r * rp), r, rp)
            return (r ** ei * rp ** ej + r ** ej * rp ** ei) * q
        return quad(inner, [0, r / 2, r])
    return quad(outer, [0, mpf(1) / 2, 1]) / pi


flatG = [[flat(i, j) for j in range(n)] for i in range(n)]
print(f"   flat side done, {time.time()-t0:.0f} s", flush=True)


def rel(A, B):
    sc = max(abs(v) for r in B for v in r)
    return max(abs(A[i][j] - B[i][j]) for i in range(n) for j in range(n)) / sc


e9, ep, e4 = rel(therm, flatG), rel(pars, direct), rel(wrongT, flatG)
line("L10 guard: u_k solves h^2 u = k^2 u (finite differences) and Parseval int pihat_i pihat_j dk = int pi_i pi_j drho",
     okode and ep <= mpf("1e-7"), f"ODE residual {mp.nstr(worst,3)}; Parseval rel err {mp.nstr(ep,3)}")
line("L9 l=2 momentum slot: <p, omega^-1 p>_flat = <pi, h^-1 coth(pi h) pi>_rho, 2x2 Gram", e9 <= mpf("1e-7"),
     f"rel err {mp.nstr(e9,3)}; flat {[[mp.nstr(v,12) for v in r] for r in flatG]}; thermal {[[mp.nstr(v,12) for v in r] for r in therm]}")
caught = e4 > mpf("1e-3")
print(("CAUGHT " if caught else "MISSED ") + f"P4 coth(pi k/2) in place of coth(pi k)  | rel err {mp.nstr(e4,3)}")
RES["plants"].append({"name": "P4 wrong temperature", "caught": bool(caught), "info": mp.nstr(e4, 3)})
if not caught:
    FAILS.append("P4")
print("FAILED LINES: " + (", ".join(FAILS) if FAILS else "none"))
RES["failed"] = FAILS
RES["flat"] = [[mp.nstr(v, 15) for v in r] for r in flatG]
RES["thermal"] = [[mp.nstr(v, 15) for v in r] for r in therm]
if OUT:
    with open(OUT, "w") as fh:
        json.dump(RES, fh, indent=1)
sys.exit(1 if FAILS else 0)
