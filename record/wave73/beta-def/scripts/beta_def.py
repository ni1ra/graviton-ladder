"""BETA-DEF: antilinear part of S (+) 1 from the definition vs the closed formula. Pass lines: ../logs/01_*.log

Definition side (uses no c, s, z): K = C^{2n} with inner product <v, w> = v^dagger G w, G = M + (i/2) Omega,
M = diag(h coth(pi h)/2, coth(pi h)/(2h)); kappa(k) = k (real vectors); H' = {v : Im<k, v> = 0 for all real k}.
Formula side: L (+) Lbar, k(z) = (c z, s zbar), z = (h^1/2 q + i h^-1/2 pi)/sqrt2.
Real-linear maps on C^{2n} are 4n x 4n real matrices in coordinates [Re v; Im v].
"""
import json, math, sys
import numpy as np

n = 6
h = np.array([0.25, 0.45, 0.7, 0.95, 1.25, 1.6])
tau = math.pi
rng = np.random.default_rng(73)
I_n, Z_n = np.eye(n), np.zeros((n, n))
Om = np.block([[Z_n, I_n], [-I_n, Z_n]])
X0 = np.block([[Z_n, I_n], [-np.diag(h**2), Z_n]])
Sym = rng.standard_normal((2*n, 2*n)); Sym = 0.5*(Sym + Sym.T)
D = Om @ Sym                                   # D^T Om + Om D = 0
coth = 1/np.tanh(math.pi*h)
c = 1/np.sqrt(-np.expm1(-2*math.pi*h)); s = np.exp(-math.pi*h)*c
W = np.sqrt(h)


def expm(A):
    k = max(0, int(math.ceil(math.log2(max(np.linalg.norm(A, 1), 1e-300)))) + 4)
    B = A/2**k
    E = np.eye(A.shape[0]); term = np.eye(A.shape[0])
    for j in range(1, 30):
        term = term @ B/j
        E = E + term
    for _ in range(k):
        E = E @ E
    return E


def rl(Mc):                                    # complex-linear map as a real matrix
    return np.block([[Mc.real, -Mc.imag], [Mc.imag, Mc.real]])


m2 = 2*n
J = rl(1j*np.eye(m2))
CONJ = np.block([[np.eye(m2), np.zeros((m2, m2))], [np.zeros((m2, m2)), -np.eye(m2)]])
fro = np.linalg.norm

# ---------- definition side ----------
M = np.diag(np.concatenate([0.5*h*coth, 0.5*coth/h]))
G = M + 0.5j*Om
Lc = np.linalg.cholesky(G)                     # G = Lc Lc^dagger ; orthonormal coords y = Lc^dagger v
Hb = np.vstack([np.eye(m2), np.zeros((m2, m2))])
_, sv, Vt = np.linalg.svd(np.hstack([G.imag, G.real]))   # Im<k, v> = k^T (ImG Re v + ReG Im v)
Hp = Vt[m2:].T                                 # null space: basis of H'
Bas = np.hstack([Hb, Hp])
Bas_inv = np.linalg.inv(Bas)


def ext(S, Sp=None):                           # S (+) Sp on H (+) H' (Sp given in the Hp basis; default 1)
    blk = np.eye(2*m2); blk[:m2, :m2] = S
    if Sp is not None:
        blk[m2:, m2:] = Sp
    return Bas @ blk @ Bas_inv


def anti(R):
    return 0.5*(R + J @ R @ J)


# ---------- formula side ----------
Kc = np.block([[np.diag(c*W), 1j*np.diag(c/W)], [np.diag(s*W), -1j*np.diag(s/W)]])/math.sqrt(2)  # k -> (c z, s zbar)
U, U_inv = rl(Kc), np.linalg.inv(rl(Kc))


def AB(R):                                     # real 2n x 2n on (q, pi)  ->  R z = A z + B zbar
    a = (W[:, None]*R[:n, :n])/W[None, :]
    d = (R[n:, n:]*W[None, :])/W[:, None]
    b = R[n:, :n]/(W[:, None]*W[None, :])
    cc = R[:n, n:]*(W[:, None]*W[None, :])
    return 0.5*((a + d) + 1j*(b - cc)), 0.5*((a - d) + 1j*(b + cc))


def beta_formula(R, cc=c, ss=s):               # antilinear part of (R (+) 0)^ on (u, v), R = S - 1 or T
    A, B = AB(R)
    M1 = np.block([[cc[:, None]*B*cc[None, :], -cc[:, None]*A*ss[None, :]],
                   [ss[:, None]*A.conj()*cc[None, :], -ss[:, None]*B.conj()*ss[None, :]]])
    return rl(M1) @ CONJ, float(np.sum(np.abs(M1)**2))


def en_form(R, z, cc=c, ss=s, drop_conj=False):   # EVEN-NONZERO: (c P, s Qbar)
    A, B = AB(R)
    P = B @ (cc**2*z.conj()) - A @ (ss**2*z)
    Q = A @ (cc**2*z) - B @ (ss**2*z.conj())
    return np.concatenate([cc*P, ss*(Q if drop_conj else Q.conj())]), P, Q


def zvec(k):
    return (W*k[:n] + 1j*k[n:]/W)/math.sqrt(2)


def cvec(x):                                   # real 4n coords -> complex 2n
    return x[:m2] + 1j*x[m2:]


def rvec(v):
    return np.concatenate([v.real, v.imag])


E0m = expm(-tau*X0)
big = expm(tau*np.block([[X0, D], [np.zeros((m2, m2)), X0]]))
T = big[:m2, m2:] @ E0m                        # int_0^tau e^{sX0} D e^{-sX0} ds
Sg = lambda g: expm(tau*(X0 + g*D)) @ E0m

res, ok = {}, {}
# L1
S1 = Sg(0.1)
sing = np.linalg.svd(Bas, compute_uv=False)
res["L1"] = {"G_min_eig": float(np.linalg.eigvalsh(G).min()), "null_sv_max": float(sv[m2:].max()) if len(sv) > m2 else 0.0,
             "rank_sv_min": float(sv[:m2].min()), "basis_sv_min": float(sing.min()),
             "sympl_err": float(fro(S1.T @ Om @ S1 - Om))}
ok["L1"] = res["L1"]["G_min_eig"] > 0 and Hp.shape[1] == m2 and res["L1"]["basis_sv_min"] > 1e-6 and res["L1"]["sympl_err"] < 1e-9
# L2
res["L2"] = float(fro(Kc.conj().T @ Kc - G)/fro(G)); ok["L2"] = res["L2"] <= 1e-10
# L3: U Hp columns must be of the form (s w, c wbar): from first slot w = u/s, then second slot = c conj(w)
r3 = 0.0
for j in range(m2):
    v = cvec(U @ Hp[:, j]); w = v[:n]/s
    r3 = max(r3, float(fro(v[n:] - c*w.conj())/fro(v)))
res["L3"] = r3; ok["L3"] = r3 <= 1e-9
# L4, L5, L7
gs = [1e-1, 1e-2, 1e-3, 1e-4]
b1_def = U @ anti(ext(T + np.eye(m2)) - np.eye(2*m2)) @ U_inv      # antilinear part of T (+) 0, by definition
b1_for, hs1 = beta_formula(T)
Ly = rl(Lc.conj().T); Ly_inv = np.linalg.inv(Ly)
rows = []
for g in gs:
    S = Sg(g)
    bd_abs = anti(ext(S))
    bd = U @ bd_abs @ U_inv
    bf, hs_f = beta_formula(S - np.eye(m2))
    hs_def = 0.5*fro(Ly @ bd_abs @ Ly_inv)**2
    rows.append({"g": g, "L4_rel": float(fro(bd - bf)/fro(bf)), "hs2_def": hs_def, "hs2_fourterm": hs_f,
                 "L5_rel": abs(hs_def/hs_f - 1), "L7_r": float(fro(bd - g*b1_def)/(g*fro(b1_def)))})
res["rows"] = rows
ok["L4"] = all(r["L4_rel"] <= 1e-7 for r in rows[:3])
ok["L5"] = all(r["L5_rel"] <= 1e-7 for r in rows[:3])
rat7 = [rows[i]["L7_r"]/rows[i + 1]["L7_r"] for i in range(3)]
res["L7_decade_ratios"] = rat7; ok["L7"] = all(8 <= x <= 12 for x in rat7)
# L6, L9, L10, plants A, B
l6 = l9a = l10 = pa = pb = 0.0; l9b = True; pa_min = pb_min = 1e9
for _ in range(5):
    k = rng.standard_normal(m2); z = zvec(k)
    e_model = Kc @ k                                    # = (c z, s zbar)
    lhs = cvec(b1_def @ rvec(e_model))
    f, P, Q = en_form(T, z)
    A, B = AB(T); Tz = A @ z + B @ z.conj()
    l6 = max(l6, float(fro(lhs - f)/fro(f)))
    l9a = max(l9a, float(fro(P + Q - Tz)/fro(Tz)))
    l9b = l9b and (fro(lhs) >= fro(Tz/np.sqrt(2*np.sinh(2*math.pi*h)))*(1 - 1e-12))
    eta = np.concatenate([1j*z/c, (1j*z).conj()/s])
    x, y = np.vdot(eta, lhs).real, np.vdot(z, Tz).imag
    l10 = max(l10, abs(x - y)/abs(y))
    fa, _, _ = en_form(T, z, cc=s, ss=c); fb, _, _ = en_form(T, z, drop_conj=True)
    pa_min = min(pa_min, float(fro(lhs - fa)/fro(lhs))); pb_min = min(pb_min, float(fro(lhs - fb)/fro(lhs)))
res.update({"L6": l6, "L9_sum": l9a, "L9_bound_holds": bool(l9b), "L10": l10, "L12_plantA_min_mismatch": pa_min,
            "L13_plantB_min_mismatch": pb_min})
ok["L6"] = l6 <= 1e-8; ok["L9"] = l9a <= 1e-12 and l9b; ok["L10"] = l10 <= 1e-9
ok["L12_plantA_caught"] = pa_min >= 1e-2; ok["L13_plantB_caught"] = pb_min >= 1e-2
# L8
k = rng.standard_normal(m2); z = zvec(k); f, _, _ = en_form(T, z)
d8 = [float(fro(cvec(U @ anti(ext(Sg(g))) @ U_inv @ rvec(Kc @ k)) - g*f)/(g*fro(f))) for g in gs]
res["L8_r"] = d8; res["L8_decade_ratios"] = [d8[i]/d8[i + 1] for i in range(3)]
ok["L8"] = all(8 <= x <= 12 for x in res["L8_decade_ratios"])
# L11: i*kappa(k(psi)) = k(i coth psi) + k'(-i csch psi), tested on the definition side
psi = zvec(rng.standard_normal(m2))
def k_of_z(zz):                                         # z -> real (q, pi)
    return np.concatenate([math.sqrt(2)*zz.real/W, math.sqrt(2)*zz.imag*W])
kpsi = k_of_z(psi)
coef = Bas_inv @ (J @ np.concatenate([kpsi, np.zeros(m2)]))
hpart, hprime = coef[:m2], cvec(U @ (Hp @ coef[m2:]))
wv = -1j*psi/np.sinh(math.pi*h)
r11 = max(float(fro(hpart - k_of_z(1j*coth*psi))/fro(hpart)),
          float(fro(hprime - np.concatenate([s*wv, c*wv.conj()]))/fro(hprime)))
res["L11"] = r11; ok["L11"] = r11 <= 1e-9
# L14 plant C: S (+) S' with k'(w) -> k'(S w)
g = 1e-2; S = Sg(g)
Kp = np.block([[np.diag(s*W), 1j*np.diag(s/W)], [np.diag(c*W), -1j*np.diag(c/W)]])/math.sqrt(2)   # k -> k'(z(k))
Hp_model = U_inv @ np.vstack([Kp.real, Kp.imag])        # columns: abstract coords of k'(z(e_j))
Cm = np.linalg.lstsq(Hp, Hp_model, rcond=None)[0]       # Hp @ Cm = Hp_model
Sp_in_Hp = Cm @ S @ np.linalg.inv(Cm)
hsC = 0.5*fro(Ly @ anti(ext(S, Sp_in_Hp)) @ Ly_inv)**2
res["L14_hs2_S+S'"] = hsC; res["L14_hs2_S+1"] = rows[1]["hs2_def"]
ok["L14_plantC_differs"] = abs(hsC/rows[1]["hs2_def"] - 1) >= 0.01
# INFO
res["INFO_hs2_expX0_plus_1"] = 0.5*fro(Ly @ anti(ext(expm(tau*X0))) @ Ly_inv)**2
res["hs2_beta1"] = hs1

print("g        L4 exact-formula rel   hs2_def          hs2_fourterm     L5 rel     L7 r(g)")
for r in rows:
    print(f"{r['g']:<8.0e} {r['L4_rel']:<22.3e} {r['hs2_def']:<16.9e} {r['hs2_fourterm']:<16.9e} {r['L5_rel']:<10.2e} {r['L7_r']:.6e}")
for kk, v in res.items():
    if kk != "rows":
        print(kk, "=", v)
for kk, v in ok.items():
    print(("PASS " if v else "FAIL ") + kk)
if len(sys.argv) > 1:
    with open(sys.argv[1], "w") as fh:
        json.dump({"res": res, "ok": {a: bool(b) for a, b in ok.items()}}, fh, indent=1)
sys.exit(0 if all(ok.values()) else 1)
