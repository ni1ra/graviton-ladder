"""Lane JSS-ID (wave 73): which graviton charge does Jensen-Sorce-Speranza's A4 refer to, on odd-parity
data of a flat unit ball?  Exact rational arithmetic.  Definitions and pass/fail lines: ../PREREG.md.
Usage: python jss_id.py L [KMAX] [--no4d]      (L = 2, 3, 4; KMAX default 2)"""
import sys, os, json, time, random
sys.dont_write_bytecode = True          # never write __pycache__ into charge-id (outside this lane's fence)
from fractions import Fraction as Fr
from math import factorial
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..", "charge-id", "scripts")))
import polyalg as P
from polyalg import (add, sub, mul, scale, const, pdiff, diff, sumall, ball, sph, gauss2, dot, trace, div,
                     tlap, curl, xTx, Tx, ric_lin, rscal_lin, rscal_2, tadd, tscale, tmulp, tzero, I3, XS,
                     R2, NLAPSE, eps, harmonic, powp, symgrad, lap)

ARGS = [a for a in sys.argv[1:] if not a.startswith("--")]
L = int(ARGS[0])
KMAX = int(ARGS[1]) if len(ARGS) > 1 else 2
DO4D = "--no4d" not in sys.argv
OUT = {"l": L, "kmax": KMAX, "checks": {}, "info": {}}
FAILS = []
T0 = time.time()
R4 = range(4)


def cl(l):
    return Fr(4 * factorial(l - 2), factorial(l + 2))


CL = cl(L)


def log(*a):
    print(*a, flush=True)


def check(name, ok, info="", kind="PASS"):
    """kind = PASS for an identity that must hold, CAUGHT for a plant or contrast that must be detected."""
    OUT["checks"][name] = {"ok": bool(ok), "kind": kind, "info": str(info)}
    tag = (kind if ok else ("MISSED" if kind == "CAUGHT" else "FAIL"))
    log(f"{tag:6s} {name}" + ("  | " + str(info) if info != "" else ""))
    if not ok:
        FAILS.append(name)


def info(name, val):
    OUT["info"][name] = str(val)
    log(f"INFO   {name}: {val}")


# ------------------------------------------------------------------ helpers (re-implemented from charge_id.py)
def spoly(coefs):
    return sumall([scale(powp(R2, k), c) for k, c in enumerate(coefs)])


def dcoefs(coefs):
    return [k * c for k, c in enumerate(coefs)][1:] or [Fr(0)]


def padd(a, b):
    n = max(len(a), len(b))
    a = list(a) + [Fr(0)] * (n - len(a))
    b = list(b) + [Fr(0)] * (n - len(b))
    return [x + y for x, y in zip(a, b)]


def KP(phi):
    g = [pdiff(phi, i + 1) for i in I3]
    v = [sumall([scale(mul(XS[j + 1], g[k]), eps(i, j, k)) for j in I3 for k in I3 if eps(i, j, k)]) for i in I3]
    K = [[scale(add(pdiff(v[j], i + 1), pdiff(v[i], j + 1)), Fr(1, 2)) for j in I3] for i in I3]
    Pt = [[add(mul(XS[i + 1], v[j]), mul(v[i], XS[j + 1])) for j in I3] for i in I3]
    return v, K, Pt


def fg_from_eta(coefs, l):
    coefs = [Fr(c) for c in coefs]
    d = dcoefs(coefs)
    if P.GAUSS[0]:
        d = padd(d, [-c for c in coefs])
    g = [Fr(-2, l + 2) * c for c in d]
    f = [Fr(2, l - 1) * c for c in padd(coefs, [-c for c in [Fr(0)] + g])]
    return f, g


def build_T(phi, coefs, l):
    v, K, Pt = KP(phi)
    f, g = fg_from_eta(coefs, l)
    return tadd(tmulp(K, spoly(f)), tmulp(Pt, spoly(g)))


def zvec(w):
    return all(not c for c in w)


def zt(A):
    return all(not A[i][j] for i in I3 for j in I3)


def tsub(A, B):
    return tadd(A, tscale(B, -1))


def transp(A):
    return [[A[j][i] for j in I3] for i in I3]


def symcurl(A):
    C = curl(A)
    return tscale(tadd(C, transp(C)), Fr(1, 2))


def crossv(w, A):
    """(w x A)_ij = eps_ikl w_k A_lj"""
    C = tzero()
    for i in I3:
        for j in I3:
            C[i][j] = sumall([scale(mul(w[k], A[l][j]), eps(i, k, l)) for k in I3 for l in I3 if eps(i, k, l)])
    return C


def symcross(w, A):
    C = crossv(w, A)
    return tscale(tadd(C, transp(C)), Fr(1, 2))


def grad(p):
    return [diff(p, i + 1) for i in I3]


def psi_of_K(K):
    return scale(xTx(curl(K)), -1)


def psidot_of_h(h):
    return xTx(curl(ric_lin(h)))


def Epsi(pa, pb, w, INT=ball, N=NLAPSE):
    gg = sumall([mul(x, y) for x, y in zip(grad(pa), grad(pb))])
    return Fr(1, 2) * INT(mul(N, gg)) + Fr(w, 2) * INT(mul(pa, pb))


def Edot(pa, pb, INT=ball, N=NLAPSE):
    return Fr(1, 2) * INT(mul(N, mul(pa, pb)))


# ------------------------------------------------------------------ 4d flows on Cauchy data
def d4(p, m):
    return pdiff(p, 0) if m == 0 else diff(p, m)


def sync_H(h, K):
    """Synchronous-gauge solution to O(t^3): H_0mu = 0, H_ij = h + 2tK - t^2 Ric_lin[h] + (t^3/3) Lap K.
    Valid in any spatial gauge for h (K is TT). On TT h it is h + 2tK + (t^2/2) Lap h + ..."""
    t = P.T_
    Rh, LK = ric_lin(h), tlap(K)
    H = [[{} for _ in R4] for _ in R4]
    for i in I3:
        for j in I3:
            H[i + 1][j + 1] = sumall([h[i][j], scale(mul(t, K[i][j]), 2), scale(mul(powp(t, 2), Rh[i][j]), -1),
                                      scale(mul(powp(t, 3), LK[i][j]), Fr(1, 3))])
    return H


def lie2(H, xi, rescale=False):
    dxi = [[pdiff(xi[r], m) for m in R4] for r in R4]
    D = [[{} for _ in R4] for _ in R4]
    for m in R4:
        for n in R4:
            acc = {}
            for r in R4:
                acc = add(acc, mul(xi[r], d4(H[m][n], r)))
                acc = add(acc, mul(H[r][n], dxi[r][m]))
                acc = add(acc, mul(H[m][r], dxi[r][n]))
            if rescale:
                acc = add(acc, scale(mul(P.T_, H[m][n]), 2))
            D[m][n] = acc
    return D


def readback(D):
    dh = [[P.at_t0(D[i + 1][j + 1]) for j in I3] for i in I3]
    dK = [[sub(scale(P.tcoef(D[i + 1][j + 1], 1), Fr(1, 2)),
               scale(add(P.at_t0(d4(D[j + 1][0], i + 1)), P.at_t0(d4(D[i + 1][0], j + 1))), Fr(1, 2)))
           for j in I3] for i in I3]
    return dh, dK


def flow(h, K, xi, rescale=False):
    return readback(lie2(sync_H(h, K), xi, rescale))


ZETA = [sub(NLAPSE, scale(powp(P.T_, 2), Fr(1, 2)))] + [scale(mul(P.T_, XS[i + 1]), -1) for i in I3]
BOOST = [P.X_, P.T_, {}, {}]
NBOOST = P.X_
BOOST_BAD = [add(P.X_, scale(R2, Fr(1, 10))), P.T_, {}, {}]


def omega(a, b, INT=ball):
    return INT(sub(dot(a[1], b[0]), dot(b[1], a[0])))


def gram_flow(basis, dl, INT=ball):
    n = len(basis)
    raw = [[Fr(1, 2) * omega(basis[i], dl[j], INT) for j in range(n)] for i in range(n)]
    sym = [[(raw[i][j] + raw[j][i]) / 2 for j in range(n)] for i in range(n)]
    asym = max(abs(raw[i][j] - raw[j][i]) for i in range(n) for j in range(n))
    return sym, asym


def polar(Q, ts):
    n = len(ts)
    d = [Q(t) for t in ts]
    M = [[Fr(0)] * n for _ in range(n)]
    for i in range(n):
        M[i][i] = d[i]
        for j in range(i + 1, n):
            M[i][j] = M[j][i] = (Q(tadd(ts[i], ts[j])) - d[i] - d[j]) / 2
    return M


def mdiff(A, B):
    return max(abs(A[i][j] - B[i][j]) for i in range(len(A)) for j in range(len(A[0])))


def mlin(A, c, B):
    return [[A[i][j] + c * B[i][j] for j in range(len(A[0]))] for i in range(len(A))]


def nonzero(A):
    return any(x != 0 for r in A for x in r)


# ------------------------------------------------------------------ geometric second variations (ball)
def vol2(h):
    tr = trace(h)
    return ball(sub(scale(mul(tr, tr), Fr(1, 8)), scale(dot(h, h), Fr(1, 4))))


def area2(h):
    tr, hh, hx, hnn = trace(h), dot(h, h), Tx(h), xTx(h)
    hx2 = sumall([mul(c, c) for c in hx])
    trt = sub(tr, hnn)
    htan2 = add(sub(hh, scale(hx2, 2)), mul(hnn, hnn))
    return sph(sub(scale(mul(trt, trt), Fr(1, 8)), scale(htan2, Fr(1, 4))))


def dVol(q):
    return Fr(1, 2) * ball(trace(q))


def dArea(q):
    return Fr(1, 2) * sph(sub(trace(q), xTx(q)))


def rscal2_prod(h):
    """O(eps^2) scalar curvature of delta + eps h, written so that only single fields are differentiated.
    polyalg.rscal_2 differentiates products, which is wrong in Gaussian mode (a product carries exp(-2r^2));
    found on the first run of this lane (logs/run1_l2.log), where the R^3 control failed for that reason."""
    d1 = [[[diff(h[i][j], k + 1) for k in I3] for j in I3] for i in I3]                     # d1[i][j][k] = d_k h_ij
    d2 = [[[[diff(d1[i][j][k], m + 1) for m in I3] for k in I3] for j in I3] for i in I3]   # d_m d_k h_ij
    G1 = [[[scale(sumall([d1[k][j][i], d1[k][i][j], scale(d1[i][j][k], -1)]), Fr(1, 2)) for j in I3] for i in I3]
          for k in I3]                                                                      # G1[k][i][j]
    dG1 = lambda m, l, i, j: scale(sumall([d2[l][j][i][m], d2[l][i][j][m], scale(d2[i][j][l][m], -1)]), Fr(1, 2))
    dG2 = lambda m, k, i, j: scale(sumall([add(mul(d1[k][l][m], G1[l][i][j]), mul(h[k][l], dG1(m, l, i, j)))
                                           for l in I3]), -1)                               # d_m G2^k_ij
    G1c = [sumall([G1[k][k][l] for k in I3]) for l in I3]
    acc = {}
    for i in I3:
        for k in I3:
            acc = add(acc, dG2(k, k, i, i))
            acc = sub(acc, dG2(i, k, i, k))
        for l in I3:
            acc = add(acc, mul(G1c[l], G1[l][i][i]))
            for k in I3:
                acc = sub(acc, mul(G1[k][i][l], G1[l][i][k]))
    return sub(acc, dot(h, ric_lin(h)))


def Bh(h, INT=ball, N=NLAPSE):
    return -INT(mul(N, rscal2_prod(h)))


def BK(Ka, Kb, INT=ball, N=NLAPSE):
    return INT(mul(N, sub(dot(Ka, Kb), mul(trace(Ka), trace(Kb)))))


P2h = lambda h: Bh(h)
P3h = lambda h: Bh(h) - 2 * area2(h)
P4h = lambda h: Bh(h) - 2 * (area2(h) - 2 * vol2(h))


def P1h(h, rescale=False, xi=ZETA):
    Z = tzero()
    dh, dK = flow(h, Z, xi, rescale)
    return Fr(1, 2) * omega((h, Z), (dh, dK))


# ================================================================== data
phi1 = harmonic(L, Fr(3, 5), Fr(4, 5), (1, 2, 3), "re")
Z3 = tzero()
nb = KMAX + 1
unit = lambda k, n=nb: [Fr(1) if i == k else Fr(0) for i in range(n)]
Tk = [build_T(phi1, unit(k), L) for k in range(nb)]
check("J0a basis tensors are TT (div = 0, tr = 0, x.T.x = 0)",
      all(zvec(div(t)) and not trace(t) and not xTx(t) for t in Tk))
if FAILS:
    sys.exit("TT construction failed: stop")
basis = [(Tk[k], Z3) for k in range(nb)] + [(Z3, Tk[k]) for k in range(nb)]
n = len(basis)
HH = [[ball(dot(Tk[i], Tk[j])) for j in range(nb)] for i in range(nb)]
ps = [psi_of_K(t) for t in Tk]
pd = [psidot_of_h(t) for t in Tk]
E2psi = [[CL * Epsi(ps[i], ps[j], 2) for j in range(nb)] for i in range(nb)]
E1psi = [[CL * Epsi(ps[i], ps[j], 1) for j in range(nb)] for i in range(nb)]
Edd = [[CL * Edot(pd[i], pd[j]) for j in range(nb)] for i in range(nb)]
PSI2 = [[CL * Fr(1, 2) * ball(mul(ps[i], ps[j])) for j in range(nb)] for i in range(nb)]
log(f"l={L} c_l={CL} basis {n}, t={time.time()-T0:.1f}s")

# ================================================================== J1 / J2: Gram blocks on the TT representative
dl_pure = [flow(*b, ZETA, False) for b in basis]
dl_w = [flow(*b, ZETA, True) for b in basis]
G1p, as1 = gram_flow(basis, dl_pure)
G1w, as1w = gram_flow(basis, dl_w)
blk = lambda M, a, b: [[M[a * nb + i][b * nb + j] for j in range(nb)] for i in range(nb)]
NKK = [[ball(mul(NLAPSE, dot(Tk[i], Tk[j]))) for j in range(nb)] for i in range(nb)]
BKK = [[BK(Tk[i], Tk[j]) for j in range(nb)] for i in range(nb)]
check("J0c own second-order scalar curvature = polyalg.rscal_2 on polynomial data (TT and non-TT)",
      all(not sub(rscal2_prod(t), rscal_2(t)) for t in Tk)
      and not sub(rscal2_prod(tadd(Tk[0], tscale(symgrad([mul(R2, XS[2]), {}, mul(XS[1], XS[3])]), 2))),
                  rscal_2(tadd(Tk[0], tscale(symgrad([mul(R2, XS[2]), {}, mul(XS[1], XS[3])]), 2)))))
check("J1a K-K block: P1 (pure Lie) = P1w (rescaled) = P2 = P3 = P4 = int_B N K.K",
      mdiff(blk(G1p, 1, 1), NKK) == 0 and mdiff(blk(G1w, 1, 1), NKK) == 0 and mdiff(BKK, NKK) == 0)
check("J1b K-K block = c_l E_2 (psi part), exactly", mdiff(NKK, E2psi) == 0)
d1 = mlin(NKK, -1, E1psi)
check("J1c K-K block - c_l E_1 (psi part) = (1/2) c_l int psi psi, non-zero (weight 1 excluded in every prescription)",
      mdiff(d1, PSI2) == 0 and nonzero(PSI2), f"[0][0] = {PSI2[0][0]}")
check("J1d no h-K cross block in P1, P1w; both flows Hamiltonian on the ball (Omega(a, delta b) symmetric)",
      not nonzero(blk(G1p, 0, 1)) and not nonzero(blk(G1w, 0, 1)) and as1 == 0 and as1w == 0)
G2h = polar(P2h, Tk)
G3h = polar(P3h, Tk)
G4h = polar(P4h, Tk)
A2 = polar(area2, Tk)
check("J2a h-h block on TT: P1w = P4 = c_l (1/2) int N psi_dot psi_dot",
      mdiff(blk(G1w, 0, 0), Edd) == 0 and mdiff(G4h, Edd) == 0)
lam = {}
for nm, M in (("P1", blk(G1p, 0, 0)), ("P3", G3h), ("P2", G2h)):
    D = mlin(M, -1, Edd)
    ratios = {D[i][j] / HH[i][j] for i in range(nb) for j in range(nb) if HH[i][j] != 0}
    lam[nm] = sorted(ratios)
check("J2b h-h block on TT: P1 - P4 = (1/2) int h.h", lam["P1"] == [Fr(1, 2)], f"lambda set {[str(x) for x in lam['P1']]}")
check("J2c h-h block on TT: P3 - P4 = 1 * int h.h", lam["P3"] == [Fr(1)], f"lambda set {[str(x) for x in lam['P3']]}")
check("J2d h-h block on TT: P2 = P3 + 2 delta^2 Area[h], and delta^2 Area is not a multiple of int h.h",
      mdiff(G2h, mlin(G3h, 2, A2)) == 0 and len(lam["P2"]) > 1, f"P2 lambda set has {len(lam['P2'])} values")
check("PL4 plant caught: claim 'P1 (literal 3.28) = P4 (fixed volume) on TT' is false",
      mdiff(blk(G1p, 0, 0), G4h) != 0, kind="CAUGHT")

# ================================================================== J3 / J4: gauge table and Section 3.3
random.seed(11 + L)
rr = lambda: [Fr(random.randint(-9, 9), random.randint(1, 5)) for _ in range(nb)]
eta = rr()
while all(c == 0 for c in eta):
    eta = rr()
v1, K1, P1t = KP(phi1)
f_, g_ = fg_from_eta(eta, L)
hTT = build_T(phi1, eta, L)
hTTp = tadd(hTT, tscale(K1, Fr(7, 3)))
wrw = [mul(spoly([-x / 2 for x in f_]), c) for c in v1]
hRW = tadd(hTT, tscale(symgrad(wrw), 2))
grw = padd(g_, [-x for x in dcoefs(f_)])
reps = {"TT": hTT, "TT'": hTTp, "RW": hRW}
check("J3a representatives: TT' is TT and differs from TT; RW = alpha(r) P; psi_dot identical on all three",
      zvec(div(hTTp)) and not trace(hTTp) and not zt(tsub(hTTp, hTT))
      and zt(tsub(hRW, tmulp(P1t, spoly(grw))))
      and not sub(psidot_of_h(hTTp), psidot_of_h(hTT)) and not sub(psidot_of_h(hRW), psidot_of_h(hTT)))
pdT = psidot_of_h(hTT)
e2 = CL * Edot(pdT, pdT)
table = {}
for nm, h in reps.items():
    table[nm] = {"P1": P1h(h), "P1w": P1h(h, True), "P2": P2h(h), "P3": P3h(h), "P4": P4h(h),
                 "int|h|^2": ball(dot(h, h)), "d2Area": area2(h), "d2Vol": vol2(h)}
    log(f"TABLE  {nm:3s} " + "  ".join(f"{k}={float(v):.6g}" for k, v in table[nm].items()))
OUT["gauge_table"] = {k: {a: str(b) for a, b in v.items()} for k, v in table.items()}
info("c_l E_2 (h part) of the datum", f"{e2} = {float(e2):.6g}")
check("J3b P4 (fixed volume) is the same on TT, TT', RW and equals c_l E_2",
      table["TT"]["P4"] == table["TT'"]["P4"] == table["RW"]["P4"] == e2 and e2 != 0)
for nm in ("P1", "P2", "P3"):
    vals = [table[g][nm] for g in ("TT", "TT'", "RW")]
    check(f"J3c {nm} is gauge dependent (values differ between TT, TT', RW)", len(set(vals)) == 3, kind="CAUGHT",
          info=", ".join(f"{float(x):.6g}" for x in vals))
vals = [table[g]["P1w"] for g in ("TT", "TT'", "RW")]
info("P1w (canonical with +2tH) on TT, TT', RW", ", ".join(f"{float(x):.6g}" for x in vals)
     + ("  (gauge dependent off TT)" if len(set(vals)) > 1 else "  (invariant)"))
fs = [const(1), harmonic(L, Fr(3, 5), Fr(4, 5), (1, 2, 3), "re"), harmonic(2, Fr(1), Fr(0), (2, 3, 1), "re"),
      harmonic(L + 1, Fr(5, 13), Fr(12, 13), (3, 1, 2), "im")]
F2 = lambda q: area2(q) - 2 * vol2(q)
dF = lambda q: dArea(q) - 2 * dVol(q)


def lie3(w, h):
    return [[sumall([mul(w[k], pdiff(h[i][j], k + 1)) for k in I3]
                    + [mul(h[k][j], pdiff(w[k], i + 1)) for k in I3]
                    + [mul(h[i][k], pdiff(w[k], j + 1)) for k in I3]) for j in I3] for i in I3]


ok_first = ok_cross = True
cross_ctrl = None
for nm, h in reps.items():
    ok_first &= dVol(h) == 0 and dArea(h) == 0
    for f in fs:
        w = [mul(f, XS[i + 1]) for i in I3]
        lw = tscale(symgrad(w), 2)
        cr = (F2(tadd(h, lw)) - F2(h) - F2(lw)) + dF(lie3(w, h))
        ok_cross &= cr == 0
# power control for the cross term: an even-type perturbation h = phi * delta must couple to f = phi
he = [[mul(phi1, const(1 if i == j else 0)) for j in I3] for i in I3]
w = [mul(phi1, XS[i + 1]) for i in I3]
lw = tscale(symgrad(w), 2)
cross_ctrl = (F2(tadd(he, lw)) - F2(he) - F2(lw)) + dF(lie3(w, he))
check("J4a Section 3.3 first-order conditions hold identically on odd data: delta Vol[h] = delta Area[h] = 0", ok_first)
check("J4b Section 3.3 extremality: cross term between h and a normal displacement f vanishes (4 choices of f, 3 gauges)",
      ok_cross)
check("J4c power control: the same cross term is non-zero for an even perturbation phi*delta with f = phi",
      cross_ctrl != 0, f"{cross_ctrl}", kind="CAUGHT")
for nm, h in reps.items():
    r0 = P3h(h) - e2
    r1 = P4h(h) - e2
    info(f"P(mu) - c_l E_2 on {nm}", f"mu=0: {float(r0):.6g}, mu=1: {float(r1):.6g}  (linear in mu, slope 4 d2Vol = {float(4*vol2(h)):.6g})")
check("J4d the region-fixing freedom mu acts on the h-h block only; no mu gives c_l E_1 (K-K block is mu-independent, J1c)",
      nonzero(PSI2))

# ================================================================== J5: flows on curvature data (pointwise identities)
Xv = [XS[i + 1] for i in I3]
ok = {k: True for k in ("psiB", "psiE", "litE", "wB", "pureB", "pl2", "Eham")}
for (h, K), (dhp, dKp), (dhw, dKw) in zip(basis, dl_pure, dl_w):
    E = ric_lin(h)
    B = tscale(symcurl(K), -1)
    dE5 = tadd(tscale(tmulp(symcurl(B), NLAPSE), -1), tscale(symcross(Xv, B), 2))
    dB5 = tadd(tmulp(symcurl(E), NLAPSE), tscale(symcross(Xv, E), -2))
    psi, psd = psi_of_K(K), psidot_of_h(h)
    dv = sumall([diff(mul(NLAPSE, g), i + 1) for i, g in enumerate(grad(psi))])
    ok["psiB"] &= not sub(xTx(dB5), mul(NLAPSE, psd))
    ok["psiE"] &= not sub(xTx(curl(dE5)), sub(dv, psi))
    ok["litE"] &= zt(tsub(tsub(ric_lin(dhp), dE5), K)) and zt(tsub(dhp, dhw))
    ok["wB"] &= zt(tsub(tscale(symcurl(dKw), -1), dB5))
    ok["pureB"] &= zt(tsub(tsub(tscale(symcurl(dKp), -1), dB5), symcurl(h)))
    dE5bad = tadd(tscale(tmulp(symcurl(B), NLAPSE), -1), tscale(symcross(Xv, B), 1))
    ok["pl2"] &= not sub(xTx(curl(dE5bad)), sub(dv, psi))
check("J5a P5 (Lie derivative of the Weyl tensor): x.dB.x = N psi_dot, identically", ok["psiB"])
check("J5b P5: x.curl(dE).x = div(N grad psi) - 1*psi, identically (weight 1, the vacuum modular flow)", ok["psiE"])
check("J5c flow induced by (3.28) on E: Ric_lin[2NK] = dE_P5 + K, identically (and delta h is the same with or without +2tH)",
      ok["litE"])
check("J5d B half with the rescaling term: - sym curl(delta^w K) = dB_P5, identically", ok["wB"])
check("J5e B half of literal (3.28) on TT: - sym curl(delta_pure K) = dB_P5 + sym curl h, identically", ok["pureB"])
check("PL2 plant caught: coefficient 1 instead of 2 on x cross B in P5 breaks the weight-1 identity", not ok["pl2"],
      kind="CAUGHT")
# generator of P5 in the transfer variables: (1/2) c_l Omega_s(a, delta_5 a) = c_l E_1 on the ball, no edge term
G5 = [[Fr(0)] * n for _ in range(n)]
psA = [psi_of_K(b[1]) for b in basis]
pdA = [psidot_of_h(b[0]) for b in basis]
d5 = []
for (h, K) in basis:
    E = ric_lin(h)
    B = tscale(symcurl(K), -1)
    dE5 = tadd(tscale(tmulp(symcurl(B), NLAPSE), -1), tscale(symcross(Xv, B), 2))
    dB5 = tadd(tmulp(symcurl(E), NLAPSE), tscale(symcross(Xv, E), -2))
    d5.append((xTx(dB5), xTx(curl(dE5))))
for i in range(n):
    for j in range(n):
        G5[i][j] = CL * Fr(1, 2) * ball(sub(mul(pdA[i], d5[j][0]), mul(psA[i], d5[j][1])))
E1full = [[CL * (Epsi(psA[i], psA[j], 1) + Edot(pdA[i], pdA[j])) for j in range(n)] for i in range(n)]
E2full = [[CL * (Epsi(psA[i], psA[j], 2) + Edot(pdA[i], pdA[j])) for j in range(n)] for i in range(n)]
check("J5f generator of P5 in (psi, psi_dot): c_l E_1 on all Gram entries; and P1w = P4 = c_l E_2 on all Gram entries",
      mdiff(G5, E1full) == 0 and mdiff(G1w, E2full) == 0 and mdiff(E1full, E2full) != 0)

# ================================================================== C1(a): Killing control, pointwise
okc = {"E": True, "B": True, "plE": True}
gN = [pdiff(NBOOST, i + 1) for i in I3]
for (h, K) in basis:
    dh, dK = flow(h, K, BOOST, False)
    E = ric_lin(h)
    B = tscale(symcurl(K), -1)
    dEc = tadd(tscale(tmulp(symcurl(B), NBOOST), -1), tscale(symcross(gN, B), -2))
    dBc = tadd(tmulp(symcurl(E), NBOOST), tscale(symcross(gN, E), 2))
    okc["E"] &= zt(tsub(ric_lin(dh), dEc))
    okc["B"] &= zt(tsub(tscale(symcurl(dK), -1), dBc))
    dhb, dKb = flow(h, K, BOOST_BAD, False)
    okc["plE"] &= zt(tsub(ric_lin(dhb), dEc)) and zt(tsub(tscale(symcurl(dKb), -1), dBc))
check("C1a Killing control (boost, lapse x_1): flow of (3.28) on (E, B) = Lie derivative of the Weyl tensor, identically",
      okc["E"] and okc["B"])
check("PL1 plant caught: lapse x_1 + r^2/10 (not Killing) breaks the control identity", not okc["plE"], kind="CAUGHT")
if not (okc["E"] and okc["B"]):
    json.dump(OUT, open(os.path.join(HERE, "..", "out", f"jss_id_l{L}.json"), "w"), indent=1)
    sys.exit("CONTROL FAILED: stop, a definition is misread")

# ================================================================== J6: 4d check that the diffeomorphism acts by P5
ETA4 = [-1, 1, 1, 1]


def Dq(q):
    return [[[sub(add(pdiff(q[d][c], b), pdiff(q[d][b], c)), pdiff(q[b][c], d)) for c in R4] for b in R4] for d in R4]


def riem_lin(q):
    memo = {}

    def R(a, b, c, d):
        k = (a, b, c, d)
        if k not in memo:
            memo[k] = scale(sumall([pdiff(pdiff(q[a][d], b), c), pdiff(pdiff(q[b][c], a), d),
                                    scale(pdiff(pdiff(q[b][d], a), c), -1), scale(pdiff(pdiff(q[a][c], b), d), -1)]),
                            Fr(1, 2))
        return memo[k]
    return R


def run4d(h, K, label):
    H = sync_H(h, K)
    A = [[scale(P.T_, -2 * ETA4[a]) if a == b else {} for b in R4] for a in R4]      # Lie_zeta eta = -2 t eta
    C = lie2(H, ZETA, False)
    DA, DB, DC = Dq(A), Dq(H), Dq(C)
    g1 = lambda D: [[[scale(D[a][b][c], Fr(ETA4[a], 2)) for c in R4] for b in R4] for a in R4]
    Gs, Gk, Gc = g1(DA), g1(DB), g1(DC)
    Gsk = [[[sub(Gc[a][b][c], scale(sumall([mul(scale(H[a][d], ETA4[a] * ETA4[d]), DA[d][b][c]) for d in R4]
                                           + [mul(scale(A[a][d], ETA4[a] * ETA4[d]), DB[d][b][c]) for d in R4]),
                                    Fr(1, 2))) for c in R4] for b in R4] for a in R4]

    def Rsk(a, b, c, d, bil=True):   # O(s kappa) coefficient of R^a_bcd at t = 0
        G = Gsk if bil else Gc
        acc = sub(pdiff(G[a][d][b], c), pdiff(G[a][c][b], d))
        if bil:
            for e in R4:
                acc = add(acc, mul(Gs[a][c][e], Gk[e][d][b]))
                acc = add(acc, mul(Gk[a][c][e], Gs[e][d][b]))
                acc = sub(acc, mul(Gs[a][d][e], Gk[e][c][b]))
                acc = sub(acc, mul(Gk[a][d][e], Gs[e][c][b]))
        return P.at_t0(acc)
    RL = riem_lin(H)
    E = ric_lin(h)
    B = tscale(symcurl(K), -1)
    dE5 = tadd(tscale(tmulp(symcurl(B), NLAPSE), -1), tscale(symcross(Xv, B), 2))
    dB5 = tadd(tmulp(symcurl(E), NLAPSE), tscale(symcross(Xv, E), -2))
    # conventions at linear order
    E0 = [[P.at_t0(RL(0, i + 1, 0, j + 1)) for j in I3] for i in I3]
    B0 = [[scale(sumall([scale(P.at_t0(RL(k + 1, l + 1, 0, j + 1)), eps(i, k, l)) for k in I3 for l in I3 if eps(i, k, l)]),
                 Fr(1, 2)) for j in I3] for i in I3]
    conv = zt(tsub(E0, E)) and zt(tsub(B0, B))
    # nonlinear diffeomorphism, O(s kappa)
    En = [[scale(Rsk(0, i + 1, 0, j + 1), -1) for j in I3] for i in I3]
    Bn = [[scale(sumall([scale(Rsk(k + 1, l + 1, 0, j + 1), eps(i, k, l)) for k in I3 for l in I3 if eps(i, k, l)]),
                 Fr(1, 2)) for j in I3] for i in I3]
    ric = all(not sumall([Rsk(a, b, a, d) for a in R4]) for b in R4 for d in R4)
    # Lie derivative of the linearised Riemann tensor field, 4d
    dz = [[pdiff(ZETA[r], m) for m in R4] for r in R4]

    def lieR(a, b, c, d):
        acc = sumall([mul(ZETA[e], pdiff(RL(a, b, c, d), e)) for e in R4])
        for e in R4:
            acc = add(acc, mul(RL(e, b, c, d), dz[e][a]))
            acc = add(acc, mul(RL(a, e, c, d), dz[e][b]))
            acc = add(acc, mul(RL(a, b, e, d), dz[e][c]))
            acc = add(acc, mul(RL(a, b, c, e), dz[e][d]))
        return P.at_t0(acc)
    El = [[lieR(0, i + 1, 0, j + 1) for j in I3] for i in I3]
    Bl = [[scale(sumall([scale(lieR(k + 1, l + 1, 0, j + 1), eps(i, k, l)) for k in I3 for l in I3 if eps(i, k, l)]),
                 Fr(1, 2)) for j in I3] for i in I3]
    # plant: bilinear term dropped = linearised Riemann of Lie_zeta H alone
    Ep = [[scale(Rsk(0, i + 1, 0, j + 1, bil=False), -1) for j in I3] for i in I3]
    return {"conv": conv, "E_nl": zt(tsub(En, dE5)), "B_nl": zt(tsub(Bn, dB5)), "ricci": ric,
            "E_lie": zt(tsub(El, dE5)), "B_lie": zt(tsub(Bl, dB5)), "plant_differs": not zt(tsub(Ep, dE5)),
            "plant_diff_is_K_plus": zt(tsub(tsub(Ep, dE5), K)), "nontrivial": not zt(dE5) or not zt(dB5)}


if DO4D:
    t4 = time.time()
    hx = build_T(phi1, [Fr(2), Fr(-1, 3)], L)
    Kx = build_T(harmonic(L, Fr(5, 13), Fr(12, 13), (3, 1, 2), "im"), [Fr(1), Fr(3, 2)], L)
    r = run4d(hx, Kx, "mixed")
    OUT["run4d"] = r
    check("J0b conventions: R^lin_0i0j[H] = Ric_lin[h] and (1/2) eps_ikl R^lin_kl0j[H] = - sym curl K at t = 0", r["conv"])
    check("J6a O(s kappa) Riemann of the dragged metric: R_0i0j = dE_P5 and (1/2) eps R_kl0j = dB_P5 at t = 0 (datum non-trivial)",
          r["E_nl"] and r["B_nl"] and r["nontrivial"])
    check("J6b control: the dragged metric is Ricci flat at O(s kappa), all 16 components at t = 0", r["ricci"])
    check("J6c covariance: (Lie_zeta R^lin[H]) has the same E and B components as P5", r["E_lie"] and r["B_lie"])
    check("PL3 plant caught: without the bilinear (Lie_zeta eta, H) term the E component is not dE_P5", r["plant_differs"],
          kind="CAUGHT", info="(4d component of a field that is not a solution; the Cauchy-data statement is J5c)")
    log(f"4d block t={time.time()-t4:.1f}s")

# ================================================================== J7 / C1(b,c) / C2: R^3 with Gaussian profiles
P.GAUSS[0] = True
ng = 2
phiA = harmonic(L, Fr(3, 5), Fr(4, 5), (1, 2, 3), "re")
phiB = harmonic(L + 1, Fr(3, 5), Fr(4, 5), (1, 2, 3), "re")
TA = [build_T(phiA, unit(k, ng), L) for k in range(ng)]
TB = [build_T(phiB, unit(k, ng), L + 1) for k in range(ng)]
check("J7a Gaussian-profile tensors are TT on R^3 (degrees l and l+1)",
      all(zvec(div(t)) and not trace(t) for t in TA + TB))
gb = [(TA[k], Z3) for k in range(ng)] + [(Z3, TA[k]) for k in range(ng)]
m = len(gb)
dgw = [flow(*b, ZETA, True) for b in gb]
dgp = [flow(*b, ZETA, False) for b in gb]
Gw, asw = gram_flow(gb, dgw, gauss2)
Gp, asp = gram_flow(gb, dgp, gauss2)
psG = [psi_of_K(b[1]) for b in gb]
pdG = [psidot_of_h(b[0]) for b in gb]
E2G = [[CL * (Epsi(psG[i], psG[j], 2, gauss2) + Edot(pdG[i], pdG[j], gauss2)) for j in range(m)] for i in range(m)]
E1G = [[CL * (Epsi(psG[i], psG[j], 1, gauss2) + Edot(pdG[i], pdG[j], gauss2)) for j in range(m)] for i in range(m)]
check("J7b two-sided, zeta used on all of R^3: (1/2) Omega(a, delta^w a) = c_l E_2 over R^3 on all Gram entries; flow Hamiltonian",
      mdiff(Gw, E2G) == 0 and asw == 0 and nonzero(E2G))
HHG = [[gauss2(dot(TA[i], TA[j])) for j in range(ng)] for i in range(ng)]
blg = lambda M, a, b: [[M[a * ng + i][b * ng + j] for j in range(ng)] for i in range(ng)]
D = mlin(blg(Gp, 0, 0), -1, blg(Gw, 0, 0))
check("J7c two-sided literal (3.28): pure-Lie charge - c_l E_2 = (1/2) int_R3 h.h on TT", mdiff(D, [[x / 2 for x in r] for r in HHG]) == 0)
BG = polar(lambda h: Bh(h, gauss2), TA)
D = mlin(BG, -1, blg(E2G, 0, 0))
lamG = sorted({D[i][j] / HHG[i][j] for i in range(ng) for j in range(ng)})
check("J7d two-sided stress-tensor integral: -int_R3 N R^(2)[h] - c_l (1/2) int N psi_dot^2 = 1 * int_R3 h.h on TT",
      lamG == [Fr(1)], f"lambda set {[str(x) for x in lamG]}")
BKG = [[BK(TA[i], TA[j], gauss2) for j in range(ng)] for i in range(ng)]
check("J7e two-sided K-K block: int_R3 N K.K = c_l E_2 (psi part) != c_l E_1 (psi part)",
      mdiff(BKG, blg(E2G, 1, 1)) == 0 and mdiff(BKG, blg(E1G, 1, 1)) != 0)
# gauge move with a Gaussian vector field (stands in for a compactly supported one)
vA = KP(phiA)[0]
wg = [mul(spoly([Fr(2), Fr(-1, 3)]), c) for c in vA]
hg0 = TA[0]
hg1 = tadd(hg0, tscale(symgrad(wg), 2))
check("J7f gauge-moved Gaussian h is not TT but has the same psi_dot", not zvec(div(hg1)) and not sub(psidot_of_h(hg1), psidot_of_h(hg0)))
bz0, bz1 = Bh(hg0, gauss2), Bh(hg1, gauss2)
check("C2 contrast: for zeta the stress-tensor integral -int N R^(2)[h] over R^3 CHANGES under the gauge move",
      bz0 != bz1, kind="CAUGHT", info=f"{float(bz0):.6g} -> {float(bz1):.6g}")
# Killing control on R^3
cb = [(TA[0], Z3), (TB[0], Z3), (Z3, TA[0]), (Z3, TB[0]), (TA[1], Z3), (Z3, TB[1])]
dcb = [flow(*b, BOOST, False) for b in cb]
Gc_, asc = gram_flow(cb, dcb, gauss2)
hs = [b[0] for b in cb]
Ks = [b[1] for b in cb]
Bc = polar(lambda h: Bh(h, gauss2, NBOOST), hs)
Bc = [[Bc[i][j] + BK(Ks[i], Ks[j], gauss2, NBOOST) for j in range(len(cb))] for i in range(len(cb))]
check("C1b Killing control on R^3 (boost): canonical charge of (3.28) = stress-tensor integral on all Gram entries, non-zero, Hamiltonian",
      mdiff(Gc_, Bc) == 0 and asc == 0 and nonzero(Bc), f"max entry {float(max(abs(x) for r in Bc for x in r)):.6g}")
hmix0 = tadd(TA[0], TB[0])
hmix1 = tadd(hmix0, tscale(symgrad(wg), 2))
b0, b1 = Bh(hmix0, gauss2, NBOOST), Bh(hmix1, gauss2, NBOOST)
check("C1c Killing control: the stress-tensor integral is invariant under the same gauge move, and non-zero",
      b0 == b1 and b0 != 0, f"{float(b0):.6g}")
P.GAUSS[0] = False

OUT["fails"] = FAILS
OUT["seconds"] = time.time() - T0
os.makedirs(os.path.join(HERE, "..", "out"), exist_ok=True)
json.dump(OUT, open(os.path.join(HERE, "..", "out", f"jss_id_l{L}.json"), "w"), indent=1)
npl = sum(1 for c in OUT["checks"].values() if c["kind"] == "CAUGHT")
log(f"DONE l={L}: {len(OUT['checks'])-len(FAILS)}/{len(OUT['checks'])} as expected ({npl} plants/contrasts), "
    f"fails={FAILS}, t={time.time()-T0:.1f}s")
