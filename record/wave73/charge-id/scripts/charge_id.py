"""Lane CHARGE-ID (wave 73): is the area-side charge on general odd-parity graviton
data equal to c_l E_2 in the transfer variables?  Exact rational arithmetic.
Usage: python charge_id.py L [KMAX]      (L = 2, 3, 4; KMAX default 3)
Definitions and pass/fail lines: ../PREREG.md."""
import sys, json, time, random
from fractions import Fraction as Fr
from math import factorial
import polyalg as P
from polyalg import (add, sub, mul, scale, const, pdiff, diff, sumall, ball, sph, gauss2,
                     dot, trace, div, tlap, curl, xTx, Tx, ric_lin, rscal_lin, rscal_2,
                     tadd, tscale, tmulp, tzero, I3, XS, R2, NLAPSE, eps, harmonic, powp)

L = int(sys.argv[1])
KMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 3
CL = Fr(4 * factorial(L - 2), factorial(L + 2))
M_ = L + 2
OUT = {"l": L, "kmax": KMAX, "checks": {}}
FAILS = []
T0 = time.time()


def log(*a):
    print(*a, flush=True)


def check(name, ok, info=""):
    OUT["checks"][name] = {"pass": bool(ok), "info": str(info)}
    log(("PASS " if ok else "FAIL ") + name + ("  " + str(info) if info != "" else ""))
    if not ok:
        FAILS.append(name)


def spoly(coefs):  # polynomial in s = r^2 from coefficient list
    return sumall([scale(powp(R2, k), c) for k, c in enumerate(coefs)])


def dcoefs(coefs):  # d/ds
    return [k * c for k, c in enumerate(coefs)][1:] or [Fr(0)]


def padd(a, b):
    n = max(len(a), len(b))
    a = list(a) + [Fr(0)] * (n - len(a))
    b = list(b) + [Fr(0)] * (n - len(b))
    return [x + y for x, y in zip(a, b)]


def KP(phi):
    """v = x cross grad phi, K_l = sym grad v, P = x(x)v + v(x)x  (plain derivatives)."""
    g = [pdiff(phi, i + 1) for i in I3]
    v = [sumall([scale(mul(XS[j + 1], g[k]), eps(i, j, k)) for j in I3 for k in I3 if eps(i, j, k)])
         for i in I3]
    K = [[scale(add(pdiff(v[j], i + 1), pdiff(v[i], j + 1)), Fr(1, 2)) for j in I3] for i in I3]
    Pt = [[add(mul(XS[i + 1], v[j]), mul(v[i], XS[j + 1])) for j in I3] for i in I3]
    return v, K, Pt


def fg_from_eta(coefs):
    """TT profile pair (f, g) as coefficient lists in s=r^2 from eta (plain) or q (Gaussian mode)."""
    coefs = [Fr(c) for c in coefs]
    d = dcoefs(coefs)
    if P.GAUSS[0]:
        d = padd(d, [-c for c in coefs])          # (q' - q) e^{-s}
    g = [Fr(-2, M_) * c for c in d]               # g = -eta'/(m r) = -(2/m) d eta/ds
    sg = [Fr(0)] + g
    f = [Fr(2, L - 1) * c for c in padd(coefs, [-c for c in sg])]
    return f, g


def build_T(phi, coefs):
    v, K, Pt = KP(phi)
    f, g = fg_from_eta(coefs)
    return tadd(tmulp(K, spoly(f)), tmulp(Pt, spoly(g)))


def is_zero_vec(w):
    return all(not c for c in w)


def is_zero_t(A):
    return all(not A[i][j] for i in I3 for j in I3)


def psi_of_K(K):
    return scale(xTx(curl(K)), -1)


def psidot_of_h(h):
    return scale(xTx(curl(tscale(ric_lin(h), -1))), -1)   # K_dot = -Ric_lin[h]


def grad(p):
    return [diff(p, i + 1) for i in I3]


def Ew(pa, pda, pb, pdb, w, N=NLAPSE):
    gg = sumall([mul(x, y) for x, y in zip(grad(pa), grad(pb))])
    return Fr(1, 2) * ball(mul(N, add(mul(pda, pdb), gg))) + Fr(w, 2) * ball(mul(pa, pb))


def delta_w(h, K, lapse=NLAPSE, drop_rescale=False):
    """delta^w H = Lie_zeta H + 2 t H on H = h + 2tK + t^2/2 Lap h + t^3/3 Lap K, H_0mu = 0,
    by the 4d component formula; read back as Cauchy data at t = 0."""
    t = P.T_
    Lh, LK = tlap(h), tlap(K)
    H = [[{} for _ in range(4)] for _ in range(4)]
    for i in I3:
        for j in I3:
            H[i + 1][j + 1] = sumall([h[i][j], scale(mul(t, K[i][j]), 2),
                                      scale(mul(powp(t, 2), Lh[i][j]), Fr(1, 2)),
                                      scale(mul(powp(t, 3), LK[i][j]), Fr(1, 3))])
    zeta = [sub(lapse, scale(powp(t, 2), Fr(1, 2)))] + [scale(mul(t, XS[i + 1]), -1) for i in I3]
    dz = [[pdiff(zeta[r], m) for m in range(4)] for r in range(4)]   # dz[r][m] = d_m zeta^r
    D = [[{} for _ in range(4)] for _ in range(4)]
    for m in range(4):
        for n in range(4):
            acc = {}
            for r in range(4):
                acc = add(acc, mul(zeta[r], pdiff(H[m][n], r)))
                acc = add(acc, mul(H[r][n], dz[r][m]))
                acc = add(acc, mul(H[m][r], dz[r][n]))
            if not drop_rescale:
                acc = add(acc, scale(mul(t, H[m][n]), 2))
            D[m][n] = acc
    dh = [[P.at_t0(D[i + 1][j + 1]) for j in I3] for i in I3]
    dK = [[sub(scale(P.tcoef(D[i + 1][j + 1], 1), Fr(1, 2)),
               scale(add(P.at_t0(pdiff(D[j + 1][0], i + 1)), P.at_t0(pdiff(D[i + 1][0], j + 1))), Fr(1, 2)))
           for j in I3] for i in I3]
    assert not D[0][0]
    return dh, dK


def omegaB(a, b):   # Omega(1,2) = int_B (K1.h2 - K2.h1); a = (h, K)
    return ball(sub(dot(a[1], b[0]), dot(b[1], a[0])))


def d2F(h):
    tr = trace(h)
    hh = dot(h, h)
    vol2 = ball(sub(scale(mul(tr, tr), Fr(1, 8)), scale(hh, Fr(1, 4))))
    hx = Tx(h)
    hnn = xTx(h)
    hx2 = sumall([mul(c, c) for c in hx])
    trt = sub(tr, hnn)
    htan2 = add(sub(hh, scale(hx2, 2)), mul(hnn, hnn))
    area2 = sph(sub(scale(mul(trt, trt), Fr(1, 8)), scale(htan2, Fr(1, 4))))
    return area2 - 2 * vol2


def Qgeo_h(h, N=NLAPSE):
    return -ball(mul(N, rscal_2(h))) - 2 * d2F(h)


def Qgeo_K(Ka, Kb, N=NLAPSE):
    return ball(mul(N, sub(dot(Ka, Kb), mul(trace(Ka), trace(Kb)))))


def closed4(ha, hb, Ka, Kb, N=NLAPSE):
    gh = sumall([mul(diff(ha[i][j], k + 1), diff(hb[i][j], k + 1)) for i in I3 for j in I3 for k in I3])
    hxa, hxb = Tx(ha), Tx(hb)
    return (ball(mul(N, dot(Ka, Kb))) + Fr(1, 4) * ball(mul(N, gh))
            - Fr(1, 2) * sph(sumall([mul(x, y) for x, y in zip(hxa, hxb)])))


def matdiff(A, B):
    return max(abs(A[i][j] - B[i][j]) for i in range(len(A)) for j in range(len(A)))


# ------------------------------------------------------------------ data
phi1 = harmonic(L, Fr(3, 5), Fr(4, 5), (1, 2, 3), "re")
phi2 = harmonic(L, Fr(5, 13), Fr(12, 13), (3, 1, 2), "im")
phi3 = harmonic(L, Fr(1), Fr(0), (2, 3, 1), "re")
Z3 = tzero()
nb = KMAX + 1
unit = lambda k: [Fr(1) if i == k else Fr(0) for i in range(nb)]
Tk = [build_T(phi1, unit(k)) for k in range(nb)]
for k in range(nb):
    check(f"D0 TT of basis tensor k={k}", is_zero_vec(div(Tk[k])) and not trace(Tk[k]) and not xTx(Tk[k]))
if FAILS:
    sys.exit("TT construction failed: stop")
basis = [(Tk[k], Z3) for k in range(nb)] + [(Z3, Tk[k]) for k in range(nb)]   # (h, K)
n = len(basis)
log(f"l={L} c_l={CL} basis {n} built, t={time.time()-T0:.1f}s")

# Green identity on a random symmetric polynomial q
random.seed(7 + L)
def randpoly(deg, nterm=6):
    p = {}
    for _ in range(nterm):
        e = [0, random.randint(0, deg), random.randint(0, deg), random.randint(0, deg)]
        p = add(p, {tuple(e): Fr(random.randint(-5, 5), random.randint(1, 4))})
    return p
q = tzero()
for i in I3:
    for j in range(i, 3):
        q[i][j] = q[j][i] = randpoly(3)
lhs = ball(mul(NLAPSE, rscal_lin(q)))
rhs = -sph(sub(trace(q), xTx(q))) + 2 * ball(trace(q))
check("D1 Green identity int N R_lin[q] = -2 dF[q] on random q", lhs == rhs, f"{lhs}")

# ------------------------------------------------------------------ Gram matrices
dl = [delta_w(*b) for b in basis]
Mom = [[Fr(1, 2) * omegaB(basis[i], dl[j]) for j in range(n)] for i in range(n)]
check("T3 Omega_B(a, delta^w b) symmetric on TT data (delta^w Hamiltonian on the ball)",
      all(Mom[i][j] == Mom[j][i] for i in range(n) for j in range(n)))
ps = [psi_of_K(b[1]) for b in basis]
pd = [psidot_of_h(b[0]) for b in basis]
ME = {w: [[CL * Ew(ps[i], pd[i], ps[j], pd[j], w) for j in range(n)] for i in range(n)] for w in (1, 2)}
log(f"canonical and scalar Grams done, t={time.time()-T0:.1f}s")
M4 = [[closed4(basis[i][0], basis[j][0], basis[i][1], basis[j][1]) for j in range(n)] for i in range(n)]
# geometric: K-K block bilinear, h-h block by polarisation, no cross term by construction
Mg = [[Fr(0)] * n for _ in range(n)]
qd = [Qgeo_h(Tk[k]) for k in range(nb)]
for i in range(nb):
    Mg[i][i] = qd[i]
    for j in range(i + 1, nb):
        Mg[i][j] = Mg[j][i] = (Qgeo_h(tadd(Tk[i], Tk[j])) - qd[i] - qd[j]) / 2
    for j in range(nb):
        Mg[nb + i][nb + j] = Qgeo_K(Tk[i], Tk[j])
log(f"geometric Gram done, t={time.time()-T0:.1f}s")

d_can_E2 = matdiff(Mom, ME[2])
d_geo_E2 = matdiff(Mg, ME[2])
d_4_can = matdiff(M4, Mom)
check("T1a A_can == c_l E_2 on all Gram entries (h-h, K-K, h-K)", d_can_E2 == 0, f"max|diff|={d_can_E2}")
check("T1b A_geo == c_l E_2 on all Gram entries", d_geo_E2 == 0, f"max|diff|={d_geo_E2}")
check("T1c Lemma-3 closed form (4) == A_can", d_4_can == 0, f"max|diff|={d_4_can}")
cross = max(abs(Mom[i][nb + j]) for i in range(nb) for j in range(nb))
check("T1d no h-K cross term in A_can", cross == 0)
OUT["gram_can_4pi"] = [[str(x) for x in r] for r in Mom]

# ------------------------------------------------------------------ controls
K0 = basis[nb][1]
A0 = Mom[nb][nb]
k5 = sph(dot(K0, K0)) / ((2 * L + 1) * (2 * L + 3))
check("C1 stationary control: A = oint|K|^2/((2l+1)(2l+3)) = c_l E_2, c_l E_1/A = (l+1)/(l+2)",
      A0 == k5 and A0 == ME[2][nb][nb] and ME[1][nb][nb] / A0 == Fr(L + 1, L + 2) and A0 != 0,
      f"A/4pi={A0}, E1/E2={ME[1][nb][nb]/A0}")
if FAILS:
    json.dump(OUT, open(f"../out/charge_id_l{L}.json", "w"), indent=1)
    sys.exit("control or primary failed; see log")
if L == 2 and KMAX >= 1:
    c = [Fr(0)] * n
    c[nb], c[nb + 1] = Fr(1), Fr(-3, 5)
    qf = lambda M: sum(c[i] * c[j] * M[i][j] for i in range(n) for j in range(n))
    r = qf(ME[1]) / qf(ME[2])
    check("C2 record value E_1/E_2 = 204/265 on momentum profile 1 - 3r^2/5, l=2", r == Fr(204, 265), f"got {r}")

# ------------------------------------------------------------------ plants
Nbad = scale(sub(const(1), scale(R2, Fr(11, 10))), Fr(1, 2))
dlb = [delta_w(*b, lapse=Nbad) for b in basis]
Mb = [[Fr(1, 4) * (omegaB(basis[i], dlb[j]) + omegaB(basis[j], dlb[i])) for j in range(n)] for i in range(n)]
check("PL1 plant caught: lapse (1 - 11/10 r^2)/2 in delta^w", matdiff(Mb, ME[2]) != 0, f"max|diff|={float(matdiff(Mb, ME[2])):.3e}")
check("PL2 plant caught: weight 1 instead of 2", matdiff(Mom, ME[1]) != 0, f"max|diff|={float(matdiff(Mom, ME[1])):.3e}")
bad = Fr(4 * factorial(L - 1), factorial(L + 2)) / CL
Mc = [[bad * x for x in r] for r in ME[2]]
if L == 2:
    # preregistered plant is degenerate at l = 2 because 0! = 1!; found on the first run, recorded, not hidden
    log("NOTE PL3 as preregistered is vacuous at l=2 (4(l-1)!/(l+2)! = c_2); PL3b added after the first run")
    OUT["PL3_vacuous_at_l2"] = matdiff(Mom, Mc) == 0
else:
    check("PL3 plant caught: c_l -> 4(l-1)!/(l+2)!", matdiff(Mom, Mc) != 0, f"max|diff|={float(matdiff(Mom, Mc)):.3e}")
Mc2 = [[Fr(L + 3, L + 2) * x for x in r] for r in ME[2]]
check("PL3b plant caught (added post-run): c_l -> c_l (l+3)/(l+2)", matdiff(Mom, Mc2) != 0, f"max|diff|={float(matdiff(Mom, Mc2)):.3e}")
dld = [delta_w(*b, drop_rescale=True) for b in basis]
Md = [[Fr(1, 4) * (omegaB(basis[i], dld[j]) + omegaB(basis[j], dld[i])) for j in range(n)] for i in range(n)]
dd = matdiff(Md, ME[2])
hh = [[Fr(1, 2) * ball(dot(basis[i][0], basis[j][0])) for j in range(n)] for i in range(n)]
Mdm = [[Md[i][j] - hh[i][j] for j in range(n)] for i in range(n)]
check("PL4 plant caught: pure Lie derivative (no +2tH) differs from c_l E_2", dd != 0, f"max|diff|={float(dd):.3e}")
check("PL4b pure-Lie charge = A + (1/2) int_B |h|^2 exactly", matdiff(Mdm, Mom) == 0)

# ------------------------------------------------------------------ T2 flow, pointwise
ok1 = ok2 = ok2w1 = True
okc = True
for b, (dh, dK) in zip(basis, dl):
    r1 = sub(psi_of_K(dK), mul(NLAPSE, psidot_of_h(b[0])))
    p = psi_of_K(b[1])
    dv = sumall([diff(mul(NLAPSE, g), i + 1) for i, g in enumerate(grad(p))])
    pdh = psidot_of_h(dh)
    r2 = sub(pdh, sub(dv, scale(p, 2)))
    r2w1 = sub(pdh, sub(dv, scale(p, 1)))
    ok1 &= not r1
    ok2 &= not r2
    ok2w1 &= not r2w1
    okc &= is_zero_vec(div(dK)) and not trace(dK) and not rscal_lin(dh)
check("T2a transferred flow: psi[delta K] = N psi_dot, identically", ok1)
check("T2b transferred flow: psi_dot[delta h] = div(N grad psi) - 2 psi, identically", ok2)
check("T2c plant caught: weight 1 in the flow identity fails", not ok2w1)
check("T2d delta^w a satisfies the linearised constraints (div dK = 0, tr dK = 0, R_lin[dh] = 0)", okc)

# ------------------------------------------------------------------ random mixed-m spot check
rr = lambda: [Fr(random.randint(-9, 9), random.randint(1, 5)) for _ in range(nb)]
h = tadd(build_T(phi1, rr()), build_T(phi2, rr()))
K = tadd(build_T(phi2, rr()), build_T(phi3, rr()))
dh, dK = delta_w(h, K)
Acan = Fr(1, 2) * omegaB((h, K), (dh, dK))
Ageo = Qgeo_h(h) + Qgeo_K(K, K)
p, pdt = psi_of_K(K), psidot_of_h(h)
E2 = CL * Ew(p, pdt, p, pdt, 2)
check("S1 random mixed-harmonic datum: A_can = A_geo = c_l E_2 exactly", Acan == Ageo == E2 and Acan != 0, f"A/4pi={float(Acan):.6g}")
OUT["spot"] = {"A_can": str(Acan), "A_geo": str(Ageo), "cE2": str(E2)}

# floating-point quadrature of the same datum, independent of the monomial integral formula
try:
    import numpy as np
    def ev(pol, X, Y, Z):
        s = np.zeros_like(X)
        for (t, a, b, c), v in pol.items():
            s += float(v) * X**a * Y**b * Z**c
        return s
    nr, nt, nf = 24, 24, 48
    xr, wr = np.polynomial.legendre.leggauss(nr); xr = (xr + 1) / 2; wr = wr / 2
    xc, wc = np.polynomial.legendre.leggauss(nt)
    ph = (np.arange(nf) + 0.5) * 2 * np.pi / nf
    Rg, Cg, Pg = np.meshgrid(xr, xc, ph, indexing="ij")
    Wg = (wr[:, None, None] * Rg**2) * wc[None, :, None] * (2 * np.pi / nf)
    Sg = np.sqrt(1 - Cg**2)
    X, Y, Z = Rg * Sg * np.cos(Pg), Rg * Sg * np.sin(Pg), Rg * Cg
    nb_int = lambda pol: float(np.sum(ev(pol, X, Y, Z) * Wg))
    a_num = 0.5 * nb_int(sub(dot(K, dh), dot(dK, h)))
    gg = sumall([mul(g, g) for g in grad(p)])
    e_num = float(CL) * (0.5 * nb_int(mul(NLAPSE, add(mul(pdt, pdt), gg))) + nb_int(mul(p, p)))
    ex = float(Acan) * 4 * np.pi
    rel = max(abs(a_num - ex), abs(e_num - ex)) / abs(ex)
    check("S2 float quadrature of A_can and c_l E_2 agree with exact value (rel <= 1e-9)", rel <= 1e-9, f"rel={rel:.2e}")
except ImportError:
    log("numpy missing: S2 skipped")

# ------------------------------------------------------------------ gauge
v1, K1, P1 = KP(phi1)
ac = [Fr(random.randint(-5, 5), random.randint(1, 3)) for _ in range(3)]
wv = [mul(spoly(ac), c) for c in v1]
h1 = build_T(phi1, rr())
from polyalg import symgrad
hg = tadd(h1, tscale(symgrad(wv), 2))
check("G0 gauge-moved h is not TT (test is not vacuous)", not is_zero_vec(div(hg)))
check("G1 psi_dot invariant under odd gauge h -> h + 2 sym grad(a(r) v)", not sub(psidot_of_h(hg), psidot_of_h(h1)))
q0, q1 = Qgeo_h(h1), Qgeo_h(hg)
check("G2 A_geo invariant under random odd gauge move", q0 == q1, f"{float(q0):.6g}")
c0, c1 = closed4(h1, h1, Z3, Z3), closed4(hg, hg, Z3, Z3)
OUT["G3_closed_form_moves"] = {"TT": str(c0), "moved": str(c1)}
log(f"INFO G3 closed form (4) applied off TT gauge: TT value {float(c0):.6g}, moved value {float(c1):.6g}, "
    + ("MOVES (formula is TT-only)" if c0 != c1 else "does not move"))
# Regge-Wheeler gauge: a = -f/2 kills the K_l component
eta = rr()
f, g = fg_from_eta(eta)
hT = build_T(phi1, eta)
wrw = [mul(spoly([-x / 2 for x in f]), c) for c in v1]
hrw = tadd(hT, tscale(symgrad(wrw), 2))
grw = padd(g, [-x for x in dcoefs(f)])       # g - f'/(2r) = g - df/ds
check("G4 Regge-Wheeler representative is g_RW(r) P (no K_l, h_AB = 0)", is_zero_t(tadd(hrw, tscale(tmulp(P1, spoly(grw)), -1))))
check("G5 A_geo in Regge-Wheeler gauge = A_geo in TT gauge", Qgeo_h(hrw) == Qgeo_h(hT), f"{float(Qgeo_h(hT)):.6g}")
check("G6 psi_dot in Regge-Wheeler gauge = TT gauge", not sub(psidot_of_h(hrw), psidot_of_h(hT)))

# ------------------------------------------------------------------ T5 (added after the first run, not preregistered)
# the transfer intertwines the dynamics: psi built from Lap K is Lap psi; same for psi_dot
ok5 = True
for k in range(nb):
    ok5 &= not sub(psi_of_K(tlap(Tk[k])), P.lap(psi_of_K(Tk[k])))
    ok5 &= not sub(psidot_of_h(tlap(Tk[k])), P.lap(psidot_of_h(Tk[k])))
check("T5 (post-hoc) psi[Lap K] = Lap psi[K] and psi_dot[Lap h] = Lap psi_dot[h]: psi obeys the scalar wave equation", ok5)

# ------------------------------------------------------------------ T4 global symplectic normalisation
P.GAUSS[0] = True
Tg = [build_T(phi1, unit(k)) for k in range(nb)]
okTT = all(is_zero_vec(div(t)) and not trace(t) for t in Tg)
check("T4a Gaussian-profile tensors are TT on R^3", okTT)
ratios = []
for i in range(nb):
    pi_ = psi_of_K(Tg[i])
    for j in range(nb):
        og = gauss2(dot(Tg[i], Tg[j]))                 # Omega_grav((0,T_i),(T_j,0)) = int K_a.h_b
        os_ = -gauss2(mul(pi_, psidot_of_h(Tg[j])))    # Omega_s(a,b) = int (psidot_a psi_b - psidot_b psi_a)
        ratios.append(og / (CL * os_) if os_ != 0 else None)
P.GAUSS[0] = False
kap = set(ratios)
check("T4b Omega_grav = kappa c_l Omega_scalar on R^3 with one constant kappa", len(kap) == 1 and None not in kap, f"kappa values: {sorted(str(x) for x in kap)}")
OUT["kappa"] = sorted(str(x) for x in kap)

OUT["fails"] = FAILS
OUT["seconds"] = time.time() - T0
json.dump(OUT, open(f"../out/charge_id_l{L}.json", "w"), indent=1)
log(f"DONE l={L}: {len(OUT['checks'])-len(FAILS)}/{len(OUT['checks'])} pass, fails={FAILS}, t={time.time()-T0:.1f}s")
