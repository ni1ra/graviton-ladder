"""even-proof: Lemma A for the even perturbation D_E (wave73).  Native Python, numpy + mpmath + fractions.

Setting (hinge-proof Section 3): sector l, rho in (0, inf), r = t = tanh(rho/2), N = 1/(1+cosh rho) = (1-r^2)/2,
dr = N drho, h^2 = -d^2/drho^2 + lam^2/sinh^2 rho, lam^2 = l(l+1).  Canonical data (q, pi) in L^2(drho):
    psi = q/(N sinh rho),  psidot = pi/(N^2 sinh rho),  (1/2) int N psidot^2 r^2 dr = (1/2) int pi^2 drho.
Z_E = (1/2) <psidot, G_D psidot> = (1/2) <pi, W pi>_{L^2(drho)},
    W(rho, rho2) = t_<^(l+1) (t_>^(-l) - t_>^(l+1)) / (2l+1),      D_E (q, pi) = (g W pi, 0).
Claim proved by hand (see final report):  ||h W|| <= 1/lam,  hence  ||D_E||_H <= sqrt(2)|g|/lam.

Lines fixed before the first run:
  V1  exact rationals: the generalised problem  -d2w/dr2 + lam^2 w/r^2 = mu w/N  on w_k = r^(l+1+2k)(1-r^2) is
      triangular with diagonal mu_k = (k+1)(2l+2k+3)  (even-tower's spectrum 1 + 1/mu_k), l = 2..6, k = 0..5.
  V2  mpmath 30 digits, l = 2, 3, j = 0, 1, 2: int W(rho, rho2) w_j(rho2) drho2 = w_j(rho)/mu_j at 3 points, rel err < 1e-20.
  V3  grid: top three eigenvalues of the W matrix within 1e-3 relative of 1/mu_j, l = 2, 3.
  V4  grid: kernel W equals N^(1/2) (h^2 + 1/4)^-1 N^(1/2) (hinge-e's form) to 2e-3 in operator norm relative, l = 2, 3.
  C1  ||D_E||_H <= sqrt(2)/lam on every grid and every l.
  C2  lam ||h W|| <= 1 on every grid and every l.
  C3  log-log slope of ||D_E||_H against lam over l = 16..128 lies in [-1.08, -0.92].
  C4  convergence: finest two grids agree to 1% on ||D_E||_H for every l.
  P1  plant: bound with the wrong power, sqrt(2)/lam^2, must be VIOLATED for some l.
  P2  plant: Neumann Green kernel must MISS 1/(2l+3) as top eigenvalue by more than 5% at l = 2.
"""
import json, math, sys, time
from fractions import Fraction
import numpy as np
import mpmath as mp

LS = [2, 3, 4, 8, 16, 32, 64, 128]
GRIDS = [(40.0, 0.04), (40.0, 0.02), (60.0, 0.02), (40.0, 0.01)]   # (RHO_MAX, H); the last two are the finest two
ok_all = True
res = {"checks": [], "rows": []}


def check(name, ok, detail=""):
    global ok_all
    ok_all &= bool(ok)
    print(f"{'PASS' if ok else 'FAIL'}  {name}  {detail}")
    res["checks"].append({"name": name, "ok": bool(ok), "detail": detail})


def plant(name, caught, detail=""):
    global ok_all
    ok_all &= bool(caught)
    print(f"{'PLANT CAUGHT' if caught else 'PLANT MISSED'}  {name}  {detail}")
    res["checks"].append({"name": "plant: " + name, "ok": bool(caught), "detail": detail})


# ---------------- V1 exact triangularity
def mu(l, k):
    return (k + 1)*(2*l + 2*k + 3)


bad = 0
for l in range(2, 7):
    for k in range(6):
        # A_r r^m = (-m(m-1) + l(l+1)) r^(m-2);  w_k = r^a - r^(a+2), a = l+1+2k
        a = l + 1 + 2*k
        lo = Fraction(-a*(a - 1) + l*(l + 1))            # coefficient of r^(a-2)
        hi = -Fraction(-(a + 2)*(a + 1) + l*(l + 1))     # coefficient of r^a
        # rhs mu w/N = 2 mu r^a: diagonal entry hi/2, the lower power couples to w_(k-1) only
        if hi/2 != mu(l, k) or lo != -2*k*(2*l + 2*k + 1):
            bad += 1
check("V1 triangular generalised problem, diagonal (k+1)(2l+2k+3), exact, l=2..6, k=0..5", bad == 0, f"{bad} mismatches")


def eig_coeffs(l, j):
    """w_j = sum_k a_k r^(l+1+2k)(1-r^2); a_(k+1) = a_k (c_k - 2 mu_j)/d_(k+1)."""
    a = [Fraction(1)]
    for k in range(j):
        ck = 2*mu(l, k)
        dk1 = 2*(k + 1)*(2*l + 2*(k + 1) + 1)
        a.append(a[-1]*(ck - 2*mu(l, j))/dk1)
    return a


# ---------------- V2 mpmath: kernel in rho reproduces the spectrum
mp.mp.dps = 30
worst = mp.mpf(0)
for l in (2, 3):
    for j in (0, 1, 2):
        a = eig_coeffs(l, j)

        def w(rho, a=a, l=l):
            t = mp.tanh(rho/2)
            return sum(mp.mpf(ak.numerator)/ak.denominator*t**(l + 1 + 2*k) for k, ak in enumerate(a))*(1 - t*t)

        def W(rho, rp, l=l):
            t, s = mp.tanh(rho/2), mp.tanh(rp/2)
            lo, hi = (t, s) if t < s else (s, t)
            return lo**(l + 1)*(hi**(-l) - hi**(l + 1))/(2*l + 1)

        for rho in (mp.mpf("0.7"), mp.mpf("2.1"), mp.mpf("4.3")):
            val = mp.quad(lambda x: W(rho, x)*w(x), [0, rho, rho + 5, 40, 90])
            worst = max(worst, abs(val*mu(l, j)/w(rho) - 1))
check("V2 W in rho variables has eigenvalues 1/((j+1)(2l+2j+3)), l=2,3, j=0,1,2 (mpmath)", worst < mp.mpf("1e-20"), f"worst rel err {mp.nstr(worst, 3)}")


# ---------------- grid objects
def build(l, RHO_MAX, H, kernel="dirichlet"):
    n = int(round(RHO_MAX/H))
    rho = H*np.arange(1, n + 1)
    M = np.zeros((n, n))
    i = np.arange(n)
    M[i, i] = 2.0/H**2
    M[i[:-1], i[:-1] + 1] = -1.0/H**2
    M[i[:-1] + 1, i[:-1]] = -1.0/H**2
    M[-1, -1] = 1.0/H**2                                   # Neumann outer wall, as hinge-proof
    M[i, i] += l*(l + 1)/np.sinh(rho)**2
    t = np.tanh(rho/2)
    lo = np.minimum.outer(t, t)
    hi = np.maximum.outer(t, t)
    if kernel == "dirichlet":
        Wm = lo*(lo/hi)**l*(1.0 - hi**(2*l + 1))/(2*l + 1)*H
    else:                                                  # planted: Neumann Green function of the ball
        Wm = lo*(lo/hi)**l*(1.0 + l/(l + 1.0)*hi**(2*l + 1))/(2*l + 1)*H
    return rho, M, Wm


def norms(l, RHO_MAX, H):
    rho, M, Wm = build(l, RHO_MAX, H)
    e, U = np.linalg.eigh(M)
    assert e.min() > 0
    p = np.sqrt(e)
    Wp = U.T @ Wm @ U
    rP = np.sqrt(p/np.tanh(math.pi*p))                     # sqrt(h coth(pi h))
    rQ = np.sqrt(p*np.tanh(math.pi*p))                     # 1/sqrt(h^-1 coth(pi h))
    nD = float(np.linalg.norm((rP[:, None]*Wp)*rQ[None, :], 2))
    nhW = float(np.linalg.norm(p[:, None]*Wp, 2))
    return nD, nhW


# ---------------- V3, V4, P2
for l in (2, 3):
    rho, M, Wm = build(l, 40.0, 0.01)
    ev = np.sort(np.linalg.eigvalsh(Wm))[::-1][:3]
    ex = np.array([1.0/mu(l, j) for j in range(3)])
    check(f"V3 grid W top eigenvalues = 1/mu_j, l={l}", np.max(np.abs(ev/ex - 1)) < 1e-3, f"{ev} vs {ex}")
    Nh = np.sqrt(1.0/(1.0 + np.cosh(rho)))
    W2 = Nh[:, None]*np.linalg.inv(M + 0.25*np.eye(len(rho)))*Nh[None, :]
    d = np.linalg.norm(W2 - Wm, 2)/np.linalg.norm(Wm, 2)
    check(f"V4 kernel = N^1/2 (h^2+1/4)^-1 N^1/2, l={l}", d < 2e-3, f"rel op-norm diff {d:.2e}")
rho, M, WN = build(2, 40.0, 0.02, kernel="neumann")
evN = np.max(np.linalg.eigvalsh(WN))
plant("P2 Neumann kernel misses 1/(2l+3) at l=2", abs(evN*7 - 1) > 0.05, f"top eigenvalue * 7 = {evN*7:.4f}")

# ---------------- the decisive table
print("\n  l      lam  " + "  ".join(f"|D_E|({R:g},{H:g})".rjust(14) for R, H in GRIDS) + "   sqrt2/lam     ratio  lam|hW|  lam|D_E|")
table = {}
for l in LS:
    lam = math.sqrt(l*(l + 1))
    vals, hws = [], []
    t0 = time.time()
    for R, H in GRIDS:
        nD, nhW = norms(l, R, H)
        vals.append(nD)
        hws.append(nhW)
    bound = math.sqrt(2)/lam
    table[l] = (lam, vals, hws, bound)
    res["rows"].append({"l": l, "lam": lam, "normD": vals, "norm_hW": hws, "bound": bound, "grids": GRIDS})
    print(f"{l:3d}  {lam:8.3f}  " + "  ".join(f"{v:14.6e}" for v in vals) + f"   {bound:10.4e}  {vals[-1]/bound:6.4f}  {lam*hws[-1]:7.4f}  {lam*vals[-1]:7.4f}   ({time.time()-t0:.0f}s)", flush=True)

check("C1 ||D_E||_H <= sqrt(2)/lam, all l, all grids", all(v <= b for (_, vs, _, b) in table.values() for v in vs),
      f"max ratio {max(v/b for (_, vs, _, b) in table.values() for v in vs):.4f}")
check("C2 lam ||h W|| <= 1, all l, all grids", all(lam*x <= 1 for (lam, _, hw, _) in table.values() for x in hw),
      f"max {max(lam*x for (lam, _, hw, _) in table.values() for x in hw):.4f}")
xs = np.log([table[l][0] for l in LS if l >= 16])
ys = np.log([table[l][1][-1] for l in LS if l >= 16])
slope = float(np.polyfit(xs, ys, 1)[0])
check("C3 slope of log||D_E|| against log lam, l = 16..128, in [-1.08, -0.92]", -1.08 <= slope <= -0.92, f"slope {slope:.4f}")
conv = max(abs(table[l][1][-1]/table[l][1][-2] - 1) for l in LS)
check("C4 finest two grids agree to 1%", conv < 0.01, f"worst rel diff {conv:.2e}")
viol = [l for l in LS if table[l][1][-1] > math.sqrt(2)/table[l][0]**2]
plant("P1 wrong power: bound sqrt(2)/lam^2 is violated", len(viol) > 0, f"violated at l = {viol}")
res["slope"] = slope
print("\nALL LINES HELD" if ok_all else "\nSOME LINE FAILED")
if len(sys.argv) > 1:
    json.dump(res, open(sys.argv[1], "w"), indent=1, default=float)
sys.exit(0 if ok_all else 1)
