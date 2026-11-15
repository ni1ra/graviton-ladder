"""Lane EVEN-TOWER, Addendum A: independent Gaussian affinity of the even tower's area-side generator against the
vacuum, sector by sector, hyperbolic frame.  Own code (covariance-matrix route).  numpy/scipy, double precision.
Usage: python shale_check.py [dirichlet|neumann|dirichlet2]   -> ../out/shale_check[_wall].json"""
import json, os, sys, time
import numpy as np
from scipy.linalg import eigh

HERE = os.path.dirname(os.path.abspath(__file__))
WALL = sys.argv[1] if len(sys.argv) > 1 else "dirichlet"      # dirichlet | neumann | dirichlet2 (box doubled)
RHO, H = (100.0, 0.05) if WALL == "dirichlet2" else (50.0, 0.025)
n = int(round(RHO / H)) - (0 if WALL == "neumann" else 1)
rho = H * np.arange(1, n + 1)
TAG = "" if WALL == "dirichlet" else "_" + WALL
Nl = 1.0 / (1.0 + np.cosh(rho))
AE = 0.0135994054803
REF = json.load(open(os.path.join(HERE, "..", "..", "hinge-e", "out", "e3_hyperbolic.json")))


def find(obj, l):
    if isinstance(obj, dict):
        if obj.get("l") == l and "1-A_g=1.0" in obj:
            return obj
        for v in obj.values():
            r = find(v, l)
            if r:
                return r
    elif isinstance(obj, list):
        for v in obj:
            r = find(v, l)
            if r:
                return r
    return None


def coth(x):
    return 1.0 / np.tanh(x)


def psqrt(M, power):
    w, Q = eigh(M)
    return (Q * w ** power) @ Q.T


def logdet(M):
    s, ld = np.linalg.slogdet(M)
    assert s > 0
    return ld


def affinity(l, g, mode):
    kin = np.full(n, 2.0 / H**2)
    if WALL == "neumann":
        kin[-1] = 1.0 / H**2
    M0 = np.diag(kin + l * (l + 1) / np.sinh(rho) ** 2) - (np.eye(n, k=1) + np.eye(n, k=-1)) / H**2
    p2, U = eigh(M0)
    p = np.sqrt(p2)
    A, B = np.eye(n), M0.copy()
    if mode in ("even", "plantN"):
        if mode == "even":
            R = (U / (p2 + 0.25)) @ U.T
            W = np.sqrt(Nl)[:, None] * R * np.sqrt(Nl)[None, :]
        else:
            W = np.diag(Nl)
        A = A + g * W
    elif mode == "odd":
        B = B + g * np.diag(Nl)
    Ah, Aih = psqrt(A, 0.5), psqrt(A, -0.5)
    w2, U1 = eigh(Ah @ B @ Ah)
    w = np.sqrt(w2)
    c0, c1 = coth(np.pi * p / 2), coth(np.pi * w / 2)
    Vf = 0.5 * (U * (c0 / p)) @ U.T + 0.5 * Ah @ (U1 * (c1 / w)) @ U1.T @ Ah
    Vp = 0.5 * (U * (c0 * p)) @ U.T + 0.5 * Aih @ (U1 * (c1 * w)) @ U1.T @ Aih
    logA = 0.5 * np.sum(np.log(c0)) + 0.5 * np.sum(np.log(c1)) - 0.5 * logdet(Vf) - 0.5 * logdet(Vp)
    return -np.expm1(logA), logA


out = {"checks": [], "rows": []}
def check(name, ok, info=""):
    out["checks"].append({"name": name, "pass": bool(ok), "info": str(info)})
    print(("PASS " if ok else "FAIL ") + name + "  " + str(info), flush=True)

t0 = time.time()
_, la = affinity(16, 0.0, "even")
check("S1 control g = 0: log A = 0", abs(la) < 1e-8, f"{la:.2e}")
r16 = find(REF, 16)
d, _ = affinity(16, 1.0, "odd")
check("S2 control (known answer): odd tower potential at l = 16 reproduces HINGE-E's stored odd-control number within 2 %",
      abs(d / r16["odd_control_1-A_g=1.0"] - 1) < 0.02, f"{d:.6e} vs {r16['odd_control_1-A_g=1.0']:.6e}; l^2(1-A)={256*d:.5f}")
out["odd_control_l16"] = float(d)
if WALL == "dirichlet2":
    json.dump(out, open(os.path.join(HERE, "..", "out", f"shale_check{TAG}.json"), "w"), indent=1)
    sys.exit(0)
if not out["checks"][-1]["pass"]:
    json.dump(out, open(os.path.join(HERE, "..", "out", f"shale_check{TAG}.json"), "w"), indent=1)
    sys.exit("S2 failed: routine wrong, nothing claimed")
ok3, means = True, []
for l in (16, 32, 64, 128):
    dp, _ = affinity(l, 1.0, "even")
    dm, _ = affinity(l, -1.0, "even")
    ref = find(REF, l)["1-A_g=1.0"]
    ok3 &= abs(dp / ref - 1) < 0.02
    means.append(float(l * l * (dp + dm) / 2))
    out["rows"].append({"l": l, "1-A_g=+1": float(dp), "1-A_g=-1": float(dm), "hinge_e_g=+1": ref, "l2_mean": means[-1], "l2_mean_over_aE": means[-1] / AE})
    print(f"l={l:4d}  1-A(g=1)={dp:.6e} (HINGE-E {ref:.6e}, ratio {dp/ref:.4f})  1-A(g=-1)={dm:.6e}  l^2 mean={means[-1]:.6f} = {means[-1]/AE:.4f} a_E  t={time.time()-t0:.0f}s", flush=True)
check("S3 even tower: 1 - A_l at g = 1 within 2 % of HINGE-E's stored values, l = 16..128", ok3)
check("S4 l^2 (1 - A_l) increasing over l = 16..128 and inside [0.6, 1.05] a_E: sum (2l+1)(1 - A_l) diverges logarithmically",
      all(b > a for a, b in zip(means, means[1:])) and all(0.6 * AE <= m <= 1.05 * AE for m in means), [round(float(m), 6) for m in means])
dpl, _ = affinity(16, 1.0, "plantN")
out["plant_W_to_N"] = {"value": float(dpl), "caught": bool(abs(dpl / r16["1-A_g=1.0"] - 1) >= 0.02)}
print(("CAUGHT " if out["plant_W_to_N"]["caught"] else "MISSED ") + f"plant W -> N in the kinetic slot: {dpl:.6e} vs {r16['1-A_g=1.0']:.6e}")
out["seconds"] = time.time() - t0
json.dump(out, open(os.path.join(HERE, "..", "out", f"shale_check{TAG}.json"), "w"), indent=1)
print(f"DONE {sum(c['pass'] for c in out['checks'])}/{len(out['checks'])} checks, t={time.time()-t0:.0f}s")
