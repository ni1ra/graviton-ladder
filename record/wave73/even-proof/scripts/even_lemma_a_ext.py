"""even-proof addendum: is the l-scaling of ||D_E||_H exactly 1/lam?  Written after run 01, before this run.

Run 01 showed lam*||D_E|| still rising at l = 32 (0.40 -> 0.83), so a fixed-window slope test is pre-asymptotic.
Lines fixed before this run (grid RHO_MAX = 40, H = 0.02 and 0.01, l = 16 .. 1024):
  E1  lam*||D_E|| is increasing in l and stays below sqrt(2); lam*||hW|| is increasing and stays below 1.
  E2  local slopes d log||D_E|| / d log lam between consecutive l are monotone and approach -1; the last one
      (512 -> 1024) lies in [-1.00, -0.95].
  E3  the two grids agree to 1% at every l.
  E4  the exact-kernel matrix and the resolvent form agree to 0.5% on ||D_E|| at l = 16 .. 128, H = 0.01.
  Diagnosis of run 01's one outlier (l = 128 on the H = 0.04 grid, lam||hW|| = 1.2845), fixed here before this run:
  the exact kernel sampled on a grid is unresolved near rho = 0 once l*H is of order 1 and then gives the spurious
  value lam||hW|| ~ l*H/4 (= 1.28 there; top singular vectors sit on single grid points rho = H, 2H, 3H).
  So this run uses W = N^(1/2) (h^2 + 1/4)^-1 N^(1/2) built from the same discrete h^2 (identity V4 of run 01,
  hinge-e T2), which has no such artefact, and E4 ties it back to the exact kernel where l*H <= 1.28.
  P3  plant: exponent 1/2 (the bound ||D_E|| <= |g| (2l+3)^(-1/2) from ||W|| alone is true but not sharp):
      lam^(1/2)*||D_E|| must FALL by more than 30% from l = 64 to l = 1024, i.e. a 1/2 exponent is rejected as the scaling.
"""
import json, math, sys, time
import numpy as np
sys.path.insert(0, ".")

LS = [16, 32, 64, 128, 256, 512, 1024]
GRIDS = [(40.0, 0.02), (40.0, 0.01)]
ok_all = True
out = {"rows": [], "checks": []}


def line(kind, name, ok, detail=""):
    global ok_all
    ok_all &= bool(ok)
    tag = {"check": ("PASS", "FAIL"), "plant": ("PLANT CAUGHT", "PLANT MISSED")}[kind][0 if ok else 1]
    print(f"{tag}  {name}  {detail}")
    out["checks"].append({"name": name, "ok": bool(ok), "detail": detail})


def norms(l, RHO_MAX, H, form="resolvent"):
    n = int(round(RHO_MAX/H))
    rho = H*np.arange(1, n + 1)
    M = np.zeros((n, n))
    i = np.arange(n)
    M[i, i] = 2.0/H**2 + l*(l + 1)/np.sinh(rho)**2
    M[i[:-1], i[:-1] + 1] = -1.0/H**2
    M[i[:-1] + 1, i[:-1]] = -1.0/H**2
    M[-1, -1] -= 1.0/H**2
    t = np.tanh(rho/2)
    lo = np.minimum.outer(t, t)
    hi = np.maximum.outer(t, t)
    if form == "kernel":
        Wm = lo*(lo/hi)**l*(1.0 - hi**(2*l + 1))/(2*l + 1)*H
    else:
        Nh = np.sqrt(1.0/(1.0 + np.cosh(rho)))
        Wm = Nh[:, None]*np.linalg.inv(M + 0.25*np.eye(n))*Nh[None, :]
    e, U = np.linalg.eigh(M)
    p = np.sqrt(e)
    Wp = U.T @ Wm @ U
    rP = np.sqrt(p/np.tanh(math.pi*p))
    rQ = np.sqrt(p*np.tanh(math.pi*p))
    return float(np.linalg.norm((rP[:, None]*Wp)*rQ[None, :], 2)), float(np.linalg.norm(p[:, None]*Wp, 2))


print("   l       lam   |D_E|(H=.02)   |D_E|(H=.01)  lam|D_E|  lam|hW|  sqrt(lam)|D_E|  local slope")
tab = []
for l in LS:
    lam = math.sqrt(l*(l + 1))
    t0 = time.time()
    v = [norms(l, R, H) for R, H in GRIDS]
    sl = (math.log(v[1][0]/tab[-1][2])/math.log(lam/tab[-1][1])) if tab else float("nan")
    tab.append((l, lam, v[1][0], v[1][1], v[0][0], sl))
    out["rows"].append({"l": l, "lam": lam, "normD_coarse": v[0][0], "normD": v[1][0], "norm_hW": v[1][1], "local_slope": sl})
    print(f"{l:5d} {lam:9.3f}  {v[0][0]:13.6e}  {v[1][0]:13.6e}  {lam*v[1][0]:8.4f}  {lam*v[1][1]:7.4f}  {math.sqrt(lam)*v[1][0]:10.4f}     {sl:8.4f}   ({time.time()-t0:.0f}s)", flush=True)

x = [r[1]*r[2] for r in tab]
y = [r[1]*r[3] for r in tab]
line("check", "E1 lam||D_E|| increasing and < sqrt(2); lam||hW|| increasing and < 1",
     all(b > a for a, b in zip(x, x[1:])) and max(x) < math.sqrt(2) and all(b > a for a, b in zip(y, y[1:])) and max(y) < 1,
     f"lam||D_E||: {x[0]:.4f} .. {x[-1]:.4f}; lam||hW||: {y[0]:.4f} .. {y[-1]:.4f}")
s = [r[5] for r in tab[1:]]
line("check", "E2 local slopes monotone towards -1, last in [-1.00, -0.95]",
     all(b < a for a, b in zip(s, s[1:])) and -1.0 <= s[-1] <= -0.95, "slopes " + ", ".join(f"{q:.4f}" for q in s))
conv = max(abs(r[4]/r[2] - 1) for r in tab)
line("check", "E3 grids agree to 1%", conv < 0.01, f"worst rel diff {conv:.2e}")
worst = 0.0
for r in tab:
    if r[0] <= 128:
        k = norms(r[0], 40.0, 0.01, form="kernel")[0]
        worst = max(worst, abs(k/r[2] - 1))
line("check", "E4 exact kernel against resolvent form, l = 16..128, H = 0.01, within 0.5%", worst < 5e-3, f"worst rel diff {worst:.2e}")
h64 = [math.sqrt(r[1])*r[2] for r in tab if r[0] == 64][0]
h1024 = math.sqrt(tab[-1][1])*tab[-1][2]
line("plant", "P3 exponent 1/2 rejected: sqrt(lam)||D_E|| falls > 30% from l=64 to 1024", h1024 < 0.7*h64, f"{h64:.4f} -> {h1024:.4f}")
print("\nALL LINES HELD" if ok_all else "\nSOME LINE FAILED")
if len(sys.argv) > 1:
    json.dump(out, open(sys.argv[1], "w"), indent=1)
sys.exit(0 if ok_all else 1)
