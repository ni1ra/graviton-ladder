# ASSEMBLY lane, wave 73. Pass lines: logs\01_prefixed_pass_lines.log (fixed before the first run).
# Usage: python assembly_check.py <out.json>
import sys, json, math
import numpy as np

rng = np.random.default_rng(73)
fails, rows, out = [], [], {}


def line(tag, ok, msg):
    print(("PASS " if ok else "FAIL ") + tag + "  " + msg)
    if not ok:
        fails.append(tag)


def expm(A):
    # scaling and squaring with a Taylor series; A real or complex, moderate norm
    nrm = np.linalg.norm(A, 2)
    s = max(0, int(math.ceil(math.log2(max(nrm, 1e-16)))) + 3)
    B = A / (2.0 ** s)
    E = np.eye(A.shape[0], dtype=A.dtype)
    term = np.eye(A.shape[0], dtype=A.dtype)
    for j in range(1, 30):
        term = term @ B / j
        E = E + term
    for _ in range(s):
        E = E @ E
    return E


# ---------------- V1-V3, PA: vacuum expectation in a truncated Fock space
def ladder(d):
    return np.diag(np.sqrt(np.arange(1, d)), 1).astype(complex)


d1 = 160
a = ladder(d1)
ev, U = np.linalg.eigh(a + a.conj().T)
om = np.zeros(d1, complex); om[0] = 1
worst = 0.0
plantA = None
for r in (0.1, 0.5, 1.0, 2.0):
    Wom = U @ (np.exp(1j * r * ev) * (U.conj().T @ om))
    lhs = np.linalg.norm(Wom - om) ** 2
    worst = max(worst, abs(lhs - 2 * (1 - math.exp(-r * r / 2))))
    if r == 1.0:
        plantA = abs(lhs - 2 * (1 - math.exp(-r * r / 4)))
line("V1", worst <= 1e-9, "||W(r)Om - Om||^2 = 2(1 - exp(-r^2/2)), worst abs err %.2e" % worst)

d2 = 40
a2 = ladder(d2); I2 = np.eye(d2)
A1 = np.kron(a2, I2); A2 = np.kron(I2, a2)
om2 = np.zeros(d2 * d2, complex); om2[0] = 1
worst = 0.0
for k in ((0.7 + 0.2j, -0.3 + 0.5j), (0.1j, 1.1), (-0.9 - 0.4j, 0.6j)):
    ak = np.conj(k[0]) * A1 + np.conj(k[1]) * A2          # a(k), antilinear in k
    ev2, U2 = np.linalg.eigh(ak + ak.conj().T)
    val = om2.conj() @ (U2 @ (np.exp(1j * ev2) * (U2.conj().T @ om2)))
    worst = max(worst, abs(val - math.exp(-(abs(k[0]) ** 2 + abs(k[1]) ** 2) / 2)))
line("V2", worst <= 1e-9, "<Om, W(k) Om> = exp(-||k||^2/2), two modes, complex k, worst abs err %.2e" % worst)
rr = np.linspace(0, 6, 6001)
f = 2 * (1 - np.exp(-rr ** 2 / 2))
line("V3", bool(np.all(np.diff(f) > 0)), "r -> 2(1 - exp(-r^2/2)) strictly increasing on [0, 6]")
caughtA = plantA > 1e-3
print(("CAUGHT" if caughtA else "MISSED") + " PA  exponent 1/4 at r = 1: mismatch %.4f" % plantA)

# ---------------- the n = 6 thermal model
n = 6
tau, g, lam = math.pi, 1.0, 4.0
R0 = rng.standard_normal((n, n))
P0 = R0 @ R0.T + 0.5 * np.eye(n)
P0 = P0 / np.linalg.norm(P0, 2)            # fixed positive matrix, norm 1
N_generic = P0 / lam
M0 = P0 / 3.0
J = np.block([[np.zeros((2 * n, 2 * n)), -np.eye(2 * n)], [np.eye(2 * n), np.zeros((2 * n, 2 * n))]])
Z = np.zeros((n, n))


def model(hmin):
    h = np.array([hmin, 0.6, 1.0, 1.5, 2.2, 3.0])
    coth = 1.0 / np.tanh(math.pi * h)
    M = np.block([[np.diag(0.5 * h * coth), Z], [Z, np.diag(0.5 * coth / h)]])
    Om = np.block([[Z, np.eye(n)], [-np.eye(n), Z]])
    G = M + 0.5j * Om
    L = np.linalg.cholesky(G).conj().T       # G = L^H L, kappa x = L x
    Kap = np.vstack([L.real, L.imag])        # 24 x 12
    return h, coth, M, Om, Kap


def symp_complement(Kap):
    # v with Im<v, kappa x> = -v^T J kappa x = 0 for all x
    u, s, vt = np.linalg.svd((J @ Kap).T)
    return vt[2 * n:].T, s


def cases(h):
    X0 = np.block([[Z, np.eye(n)], [-np.diag(h * h), Z]])
    Dg = np.block([[Z, Z], [-N_generic, Z]])
    Ns = np.diag(h) @ M0 @ np.diag(h) / lam
    Ds = np.block([[Z, Z], [-Ns, Z]])
    De = np.block([[Z, M0 / lam], [Z, Z]])
    res = {}
    for name, D in (("a odd generic", Dg), ("b odd structured", Ds), ("c even", De)):
        res[name] = expm(tau * (X0 + g * D)) @ expm(-tau * X0)
    return res


def four_terms(S, h):
    # BETA-DEF section 1: (S - 1) z = A z + B zbar in z = (h^1/2 q + i h^-1/2 pi)/sqrt 2
    P = np.diag(np.concatenate([np.sqrt(h), 1 / np.sqrt(h)]) / math.sqrt(2))
    R = P @ (S - np.eye(2 * n)) @ np.linalg.inv(P)
    R11, R12, R21, R22 = R[:n, :n], R[:n, n:], R[n:, :n], R[n:, n:]
    A = 0.5 * ((R11 + R22) + 1j * (R21 - R12))
    B = 0.5 * ((R11 - R22) + 1j * (R21 + R12))
    c = np.diag(1 / np.sqrt(1 - np.exp(-2 * math.pi * h)))
    s = np.diag(np.sqrt(np.diag(c) ** 2 - 1))
    fro = lambda X: np.linalg.norm(X, "fro") ** 2
    return fro(c @ B @ c) + fro(c @ A @ s) + fro(s @ A @ c) + fro(s @ B @ s)


def build(Kap, Hp, S):
    Bas = np.hstack([Kap, Hp])
    Sh = np.hstack([Kap @ S, Hp]) @ np.linalg.inv(Bas)
    beta = 0.5 * (Sh + J @ Sh @ J)
    return Sh, beta


def rel(x, y):
    return abs(x - y) / max(abs(x), abs(y), 1e-300)


print()
print("%-8s %-18s %10s %12s %12s %12s %12s %10s" % ("h_min", "case", "cos(H,H')", "||S||_H", "||S-hat||_K", "||beta||_op", "||beta||_HS", "S2 rel"))
okG1 = okG2 = okS1 = okS2 = okS3 = okS4 = okS5 = okS6 = True
wG1 = wS1 = wS2 = wS6 = 0.0
for hmin in (0.25, 0.05, 0.01, 0.002):
    h, coth, M, Om, Kap = model(hmin)
    for _ in range(20):
        x, y = rng.standard_normal(2 * n), rng.standard_normal(2 * n)
        kx, ky = Kap @ x, Kap @ y
        e1 = rel(kx @ kx, x @ M @ x)
        e2 = rel(-kx @ J @ ky, 0.5 * (x @ Om @ y))
        wG1 = max(wG1, e1, e2)
    Hp, sv = symp_complement(Kap)
    smin = np.linalg.svd(np.hstack([Kap @ np.linalg.inv(np.linalg.cholesky(M).T), Hp]), compute_uv=False)[-1]
    QH = np.linalg.qr(Kap)[0]
    cosang = np.linalg.svd(QH.T @ Hp, compute_uv=False)[0]
    okG2 = okG2 and Hp.shape[1] == 2 * n and smin > 0
    Lc = np.linalg.cholesky(M).T                      # thermal norm ||x||_H = |Lc x|
    for name, S in cases(h).items():
        Sh, beta = build(Kap, Hp, S)
        nS = np.linalg.norm(Sh, 2)
        nb = np.linalg.norm(beta, 2)
        hs = math.sqrt(0.5) * np.linalg.norm(beta, "fro")
        nSH = np.linalg.norm(Lc @ S @ np.linalg.inv(Lc), 2)
        eS1 = np.linalg.norm(Sh.T @ J @ Sh - J, 2) / max(1.0, nS * nS)
        eS2 = rel(nS, nb + math.sqrt(1 + nb * nb))
        wS1, wS2 = max(wS1, eS1), max(wS2, eS2)
        okS3 = okS3 and nS <= (hs + math.sqrt(1 + hs * hs)) * (1 + 1e-9)
        mx = 0.0
        for _ in range(50):
            e = Kap @ rng.standard_normal(2 * n)
            be = 0.5 * (Sh @ e + J @ (Sh @ (J @ e)))
            mx = max(mx, np.linalg.norm(be) / np.linalg.norm(e))
        okS4 = okS4 and mx <= nb * (1 + 1e-9) and nb <= hs * (1 + 1e-9)
        if name[0] in "bc":
            okS5 = okS5 and nSH <= math.exp(math.sqrt(2) * g * tau / lam) * (1 + 1e-9)
        eS6 = rel(hs * hs, four_terms(S, h))
        wS6 = max(wS6, eS6)
        print("%-8g %-18s %10.6f %12.6f %12.6f %12.6f %12.6f %10.1e" % (hmin, name, cosang, nSH, nS, nb, hs, eS2))
        rows.append(dict(hmin=hmin, case=name, cos_angle=cosang, S_H=nSH, Shat_K=nS, beta_op=nb, beta_HS=hs,
                         S1=eS1, S2=eS2, S6=eS6, max_ratio_on_H=mx))
    if hmin == 0.25:
        # plant B: orthogonal complement in place of the symplectic complement
        u, s_, vt = np.linalg.svd(Kap.T)
        Hbad = vt[2 * n:].T
        pb = []
        for name, S in cases(h).items():
            Sh, beta = build(Kap, Hbad, S)
            nS, nb = np.linalg.norm(Sh, 2), np.linalg.norm(beta, 2)
            pb.append(max(np.linalg.norm(Sh.T @ J @ Sh - J, 2) / max(1.0, nS * nS), rel(nS, nb + math.sqrt(1 + nb * nb))))
        caughtB = all(v > 1e-3 for v in pb)
print()
line("G1", wG1 <= 1e-9, "kappa reproduces M and Omega/2, worst rel err %.2e" % wG1)
line("G2", okG2, "dim H' = 12 and H + H' spans K at every h_min")
line("S1", wS1 <= 1e-6, "S-hat symplectic on K, worst rel err %.2e" % wS1)
line("S2", wS2 <= 1e-6, "||S-hat|| = ||beta||_op + sqrt(1 + ||beta||_op^2), worst rel err %.2e" % wS2)
line("S3", okS3, "||S-hat|| <= ||beta||_HS + sqrt(1 + ||beta||_HS^2)")
line("S4", okS4, "||beta e||/||e|| <= ||beta||_op <= ||beta||_HS on 50 random e in H per row")
line("S5", okS5, "||S||_H <= exp(sqrt2 g tau/lambda) = %.4f in cases b, c" % math.exp(math.sqrt(2) * g * tau / lam))
line("S6", wS6 <= 1e-6, "(1/2)||beta||_F^2 = BETA-DEF's four terms, worst rel err %.2e" % wS6)
print(("CAUGHT" if caughtB else "MISSED") + " PB  orthogonal complement used for H': violations " + ", ".join("%.3f" % v for v in pb))

# ---------------- T1: the tower sum
c, l1 = 2.110e-3, 14714
pref = (c * g / 2) ** 2
okT = True
prev = None
print()
for L in (10 ** 5, 10 ** 6, 10 ** 7, 10 ** 8):
    l = np.arange(l1, L + 1, dtype=np.float64)
    Pm = pref * float(np.sum(1.0 / l + 1.0 / (l + 1.0)))        # (2l+1)/(l(l+1)) = 1/l + 1/(l+1)
    lower = pref * (math.log((L + 1) / l1) + math.log((L + 2) / (l1 + 1)))
    okT = okT and Pm >= lower
    inc = None if prev is None else Pm - prev
    print("T1 L = %.0e  partial sum %.6e  lower bound %.6e  decade increment %s" % (L, Pm, lower, "-" if inc is None else "%.6e" % inc))
    if L >= 10 ** 7:
        okT = okT and rel(inc, pref * 2 * math.log(10)) <= 1e-3
    prev = Pm
    out["T1_%d" % L] = Pm
line("T1", okT, "partial sums above the log lower bound; decade increment -> (c g/2)^2 2 ln 10 = %.6e" % (pref * 2 * math.log(10)))

# ---------------- O1: odd remainder constants (numeric, not certified)
p = np.arange(1e-6, 6.0, 1e-5)
fo = (p / np.tanh(math.pi * p)) ** 3 * p ** 2 * np.exp(-3 * p ** 2)
s5 = float(fo.max())
line("O1", s5 <= 0.14, "sup (p coth pi p)^3 p^2 exp(-3p^2) = %.6f at p = %.4f" % (s5, p[fo.argmax()]))
u2 = 4 * 65536 / 109395
J0lo, c_odd = 0.363568, 2.038e-3
E_lo = J0lo / math.sqrt(u2) / math.sqrt(2 * math.pi)
N3_hi = math.sqrt(0.5 * s5 * u2)
fac = 0.5 * (1 + N3_hi / E_lo)
C0 = fac * math.pi ** 2
thr = c_odd / (2 * C0)
print("O1 derived [N]: ||e_l|| >= %.5f, ||coth^{3/2} psi_l|| <= %.5f, factor %.4f, C0 = %.3f, g/lambda <= %.3e (before the exp(x) factor), l1_odd ~ %.0f |g|"
      % (E_lo, N3_hi, fac, C0, thr, 1 / thr))
out.update(odd_s5=s5, odd_E_lo=E_lo, odd_N3_hi=N3_hi, odd_C0=C0, odd_threshold=thr)

out.update(rows=rows, fails=fails, plantA=bool(caughtA), plantB=bool(caughtB))
with open(sys.argv[1], "w") as fh:
    json.dump(out, fh, indent=1)
print()
print("FAILED LINES: " + (", ".join(fails) if fails else "none") + "; plants caught: A %s, B %s" % (caughtA, caughtB))
sys.exit(0 if (not fails and caughtA and caughtB) else 1)
