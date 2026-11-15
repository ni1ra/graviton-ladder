"""hinge-proof U1..U5, C1 (see ../PREREG.md). Independent of hinge/scripts. float64.

Sector l of the conformally coupled field on R x H^3, u = sinh(rho) * radial function:
    K0 = -d^2/drho^2 + lam^2/sinh^2(rho),  lam^2 = l(l+1)   (Dirichlet at 0, Neumann outer wall)
    Kg = K0 + g * Npot(rho),  Npot = N = 1/(1+cosh rho)  (sharp)  or  N^2 (lapse-weighted control)
One-particle real space of the ball: data (q, pi) with norm
    ||k||^2 = 1/2 <q, p coth(pi p) q> + 1/2 <pi, coth(pi p)/p pi>,   p = sqrt(K0)   (beta = 2 pi).
S_t = exp(2 pi t X_g) exp(-2 pi t X_0),  X = [[0, 1], [-K, 0]].
"""
import json, math, sys, time
import numpy as np

RHO_MAX, H, T = 50.0, 0.04, 0.5
LS = [4, 8, 16, 32, 64, 128]
G0 = 1.0/64
n = int(round(RHO_MAX/H))
rho = H*np.arange(1, n + 1)
Nlap = 1.0/(1.0 + np.cosh(rho))
tau = 2*math.pi*T


def k0(l):
    M = np.zeros((n, n))
    i = np.arange(n)
    M[i, i] = 2.0/H**2
    M[i[:-1], i[:-1] + 1] = -1.0/H**2
    M[i[:-1] + 1, i[:-1]] = -1.0/H**2
    M[-1, -1] = 1.0/H**2                      # Neumann outer wall
    M[i, i] += l*(l + 1)/np.sinh(rho)**2
    return M


def smat(p, U0, Kg):
    """S_t in the K0 eigenbasis, as a 2n x 2n real matrix acting on (q, pi)."""
    w2, U = np.linalg.eigh(Kg)
    assert w2.min() > 0, w2.min()
    w = np.sqrt(w2)
    R = U0.T @ U
    C = (R*np.cos(tau*w)) @ R.T
    Si = (R*(np.sin(tau*w)/w)) @ R.T
    So = (R*(w*np.sin(tau*w))) @ R.T
    E = np.block([[C, Si], [-So, C]])
    c0, s0 = np.cos(tau*p), np.sin(tau*p)
    F = np.block([[np.diag(c0), np.diag(-s0/p)], [np.diag(p*s0), np.diag(c0)]])
    return E @ F


def hs_beta(S, p):
    """Exact squared Hilbert-Schmidt norm of the antilinear part of S (+) 1 on the thermofield space."""
    W = np.sqrt(p)
    Sqq, Sqp, Spq, Spp = S[:n, :n], S[:n, n:], S[n:, :n], S[n:, n:]
    a = (W[:, None]*Sqq)/W[None, :]
    d = (Spp*W[None, :])/W[:, None]
    b = Spq/(W[:, None]*W[None, :])
    c = Sqp*(W[:, None]*W[None, :])
    A = 0.5*((a + d) + 1j*(b - c))
    B = 0.5*((a - d) + 1j*(b + c))
    cth = 1.0/np.sqrt(-np.expm1(-2*math.pi*p))       # cosh(theta)
    sth = np.exp(-math.pi*p)*cth                     # sinh(theta)
    A1 = A - np.eye(n)
    t1 = np.sum(np.abs(cth[:, None]*A1*sth[None, :])**2)
    t2 = np.sum(np.abs(cth[:, None]*B*cth[None, :])**2)
    t3 = np.sum(np.abs(sth[:, None]*A1*cth[None, :])**2)
    t4 = np.sum(np.abs(sth[:, None]*B*sth[None, :])**2)
    return float(t1 + t2 + t3 + t4)


def opnorm(Mx, p):
    """Operator norm of a real 2n x 2n matrix in the one-particle norm on H."""
    coth = 1.0/np.tanh(math.pi*p)
    gdiag = np.concatenate([0.5*p*coth, 0.5*coth/p])
    r = np.sqrt(gdiag)
    return float(np.linalg.norm((r[:, None]*Mx)/r[None, :], 2))


rows = []
for l in LS:
    t0 = time.time()
    lam = math.sqrt(l*(l + 1))
    M0 = k0(l)
    p2, U0 = np.linalg.eigh(M0)
    p = np.sqrt(p2)
    row = {"l": l, "lam": lam, "p_min": float(p[0])}
    row["U1_mu"] = lam*float(np.linalg.norm((Nlap[:, None]*U0)/p[None, :], 2))
    I2 = np.eye(2*n)
    for tag, pot in (("sharp", Nlap), ("lapse", Nlap**2)):
        S = {}
        for g in (1.0, -1.0, G0, -G0):
            S[g] = smat(p, U0, M0 + g*np.diag(pot))
        E = {g: hs_beta(S[g], p) for g in S}
        E1 = (E[G0] + E[-G0])/(2*G0**2)
        row[f"{tag}_E1"] = E1
        row[f"{tag}_lam2_E1"] = lam**2*E1
        for g in (1.0, -1.0):
            row[f"{tag}_E_g={g:+.0f}"] = E[g]
            row[f"{tag}_dev_g={g:+.0f}"] = E[g]/E1 - 1.0
        row[f"{tag}_l2_Emean"] = l*l*0.5*(E[1.0] + E[-1.0])
        if tag == "sharp":
            B1 = (S[G0] - S[-G0])/(2*G0)
            for g in (1.0, -1.0):
                bound_D = 2*math.sqrt(2)*math.pi*abs(g)/lam
                nu = opnorm(S[g] - I2, p)
                r = opnorm(S[g] - I2 - g*B1, p)
                row[f"U2_nu_g={g:+.0f}"] = nu
                row[f"U2_bound_g={g:+.0f}"] = math.expm1(bound_D*T)
                row[f"U3_lam2_r_g={g:+.0f}"] = lam**2*r
                row[f"U3_bound_g={g:+.0f}"] = 4*math.pi**2*g*g*T*T*math.exp(bound_D*T)*1.02
    row["secs"] = round(time.time() - t0, 1)
    rows.append(row)
    print(json.dumps(row), flush=True)

by = {r["l"]: r for r in rows}
checks = {}
checks["U1"] = all(r["U1_mu"] <= 1 + 1e-9 for r in rows)
checks["U2"] = all(r[f"U2_nu_g={g}"] <= r[f"U2_bound_g={g}"] for r in rows for g in ("+1", "-1"))
checks["U3"] = all(r[f"U3_lam2_r_g={g}"] <= r[f"U3_bound_g={g}"] for r in rows for g in ("+1", "-1"))
ratios = {g: abs(by[128][f"sharp_dev_g={g}"])/abs(by[16][f"sharp_dev_g={g}"]) for g in ("+1", "-1")}
checks["U4"] = all(v <= 0.25 for v in ratios.values()) and by[128]["sharp_lam2_E1"] >= by[16]["sharp_lam2_E1"] > 0
hinge = {16: 0.244, 64: 0.278}
u5 = {l: by[l]["sharp_l2_Emean"]/hinge[l] for l in hinge}
checks["U5"] = all(abs(v - 1) <= 0.03 for v in u5.values())
c_sharp = by[64]["sharp_l2_Emean"]/by[16]["sharp_l2_Emean"]
c_lapse = by[64]["lapse_l2_Emean"]/by[16]["lapse_l2_Emean"]
checks["C1"] = (c_lapse <= 0.5) and (c_sharp >= 1.0)
summary = {"checks": checks, "U4_dev_ratio_128_over_16": ratios, "U5_ratio_to_hinge": u5,
           "C1_ratio_64_over_16": {"sharp": c_sharp, "lapse": c_lapse},
           "params": {"rho_max": RHO_MAX, "h": H, "n": n, "t": T, "g0": G0}}
print("SUMMARY", json.dumps(summary), flush=True)
out = sys.argv[1] if len(sys.argv) > 1 else None
if out:
    with open(out, "w") as f:
        json.dump({"rows": rows, "summary": summary}, f, indent=1)
