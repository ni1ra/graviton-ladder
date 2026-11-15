"""Lane EVEN-TOWER (wave 73): the fixed-volume second-order area-side charge on EVEN-parity graviton data of the
flat unit ball, built from the geometric definition with the extremisation over the edge, the corner time shift and
the slice; its gauge invariance; its form in the chain variables; the vacuum side.  Exact rational arithmetic.
Usage: python even_tower.py L [KMAX]      (L >= 2; KMAX default 3)
Definitions, hand predictions and pass/fail lines: ../PREREG.md.
polyalg.py is imported read-only from charge-id.  d2F, KP, the odd build and the Fourier helpers are copied from
charge-id/scripts/charge_id.py and subspace/scripts/subspace.py (both run at import, so they cannot be imported)."""
import sys, os, json, time, random
from fractions import Fraction as Fr
from math import factorial

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..", "charge-id", "scripts")))
import polyalg as P
from polyalg import (add, sub, mul, scale, const, pdiff, diff, sumall, ball, sph, gauss2, dot, trace, div, lap,
                     tlap, curl, xTx, Tx, ric_lin, rscal_lin, rscal_2, tadd, tscale, tmulp, tzero, I3, XS, R2,
                     NLAPSE, eps, harmonic, powp, symgrad, dfact)

L = int(sys.argv[1])
KMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 3
CL = Fr(4 * factorial(L - 2), factorial(L + 2))
OUT = {"l": L, "kmax": KMAX, "checks": {}}
FAILS, MISSED = [], []
T0 = time.time()
OUTP = os.path.join(HERE, "..", "out", f"even_tower_l{L}.json")


def log(*a):
    print(*a, flush=True)


def check(name, ok, info=""):
    OUT["checks"][name] = {"pass": bool(ok), "info": str(info)}
    log(("PASS " if ok else "FAIL ") + name + ("  " + str(info) if info != "" else ""))
    if not ok:
        FAILS.append(name)


def plant(name, caught, info=""):
    OUT["checks"][name] = {"pass": bool(caught), "plant": True, "info": str(info)}
    log(("CAUGHT " if caught else "MISSED ") + name + ("  " + str(info) if info != "" else ""))
    if not caught:
        MISSED.append(name)


def dump():
    OUT["fails"], OUT["plants_missed"], OUT["seconds"] = FAILS, MISSED, time.time() - T0
    json.dump(OUT, open(OUTP, "w"), indent=1)


def stop(msg):
    log("STOP: " + msg)
    OUT["stopped"] = msg
    dump()
    sys.exit(2)


# ------------------------------------------------------------------ small helpers
def spoly(coefs):
    return sumall([scale(powp(R2, k), c) for k, c in enumerate(coefs)])


def padd(a, b):
    n = max(len(a), len(b))
    a = list(a) + [Fr(0)] * (n - len(a))
    b = list(b) + [Fr(0)] * (n - len(b))
    return [x + y for x, y in zip(a, b)]


def pscale(a, c):
    return [Fr(c) * x for x in a]


def ds(coefs):
    """d/ds of a profile; in Gaussian mode the profile is coefs(s) e^{-s}."""
    d = [k * c for k, c in enumerate(coefs)][1:] or [Fr(0)]
    if P.GAUSS[0]:
        d = padd(d, pscale(coefs, -1))
    return d


def unit(k):
    return [Fr(0)] * k + [Fr(1)]


def grad(p):
    return [diff(p, i + 1) for i in I3]


def gdot(pa, pb):
    return sumall([mul(x, y) for x, y in zip(grad(pa), grad(pb))])


def xgrad(p):
    return sumall([mul(XS[i + 1], diff(p, i + 1)) for i in I3])


def hess(p):
    return [[diff(diff(p, i + 1), j + 1) for j in I3] for i in I3]


def is_zero_vec(w):
    return all(not c for c in w)


def is_zero_t(A):
    return all(not A[i][j] for i in I3 for j in I3)


def sym(A):
    return [[scale(add(A[i][j], A[j][i]), Fr(1, 2)) for j in I3] for i in I3]


def symcurl(A):
    return sym(curl(A))


def symxcross(Mt):
    C = tzero()
    for i in I3:
        for j in I3:
            acc = {}
            for k in I3:
                for l in I3:
                    e = eps(i, k, l)
                    if e:
                        acc = add(acc, scale(mul(XS[k + 1], Mt[l][j]), e))
            C[i][j] = acc
    return sym(C)


def tcomb(cs, Ts):
    A = tzero()
    for c, T in zip(cs, Ts):
        if c != 0:
            A = tadd(A, tscale(T, c))
    return A


def solve(M, rhs):
    """exact solution of M c = rhs for a possibly singular symmetric M; free variables set to 0; None if inconsistent."""
    n = len(M)
    A = [list(M[i]) + [rhs[i]] for i in range(n)]
    piv, r = [], 0
    for c in range(n):
        p = next((i for i in range(r, n) if A[i][c] != 0), None)
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        pv = A[r][c]
        A[r] = [x / pv for x in A[r]]
        for i in range(n):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        piv.append(c)
        r += 1
    if any(A[i][n] != 0 for i in range(r, n)):
        return None, r
    sol = [Fr(0)] * n
    for i, c in enumerate(piv):
        sol[c] = A[i][n]
    return sol, r


def det(M):
    M = [list(r) for r in M]
    n, d = len(M), Fr(1)
    for c in range(n):
        p = next((i for i in range(c, n) if M[i][c] != 0), None)
        if p is None:
            return Fr(0)
        if p != c:
            M[c], M[p] = M[p], M[c]
            d = -d
        d *= M[c][c]
        for i in range(c + 1, n):
            f = M[i][c] / M[c][c]
            M[i] = [x - f * y for x, y in zip(M[i], M[c])]
    return d


def prof(p, phi):
    """coefficients F_k with p = sum_k F_k s^k phi exactly (asserted)."""
    if not p:
        return [Fr(0)]
    degs = sorted({a + b + c for (t, a, b, c) in p})
    out = []
    kmax = (degs[-1] - L) // 2
    for k in range(kmax + 1):
        ref = mul(powp(R2, k), phi)
        key = next(iter(ref))
        out.append(Fr(p.get(key, 0)) / ref[key])
    assert not sub(p, mul(spoly(out), phi)), "not of the form F(s) phi"
    return out


def GD(F, neumann=False):
    """U with -Lap(U phi) = F phi and U(1) = 0 (Dirichlet); Lap(s^k phi) = 2k(2k+2l+1) s^(k-1) phi."""
    U = [Fr(0)] + [-F[k] / (2 * (k + 1) * (2 * k + 2 * L + 3)) for k in range(len(F))]
    if neumann:
        U[0] = -sum(U[k] * (2 * k + L) for k in range(1, len(U))) / L
    else:
        U[0] = -sum(U[1:])
    return U


# ------------------------------------------------------------------ even TT tensors (companion A8 Lemma 7, n = 3)
def ebasis(phi):
    g = [pdiff(phi, i + 1) for i in I3]
    Hs = [[pdiff(g[i], j + 1) for j in I3] for i in I3]
    Pp = [[add(mul(XS[i + 1], g[j]), mul(g[i], XS[j + 1])) for j in I3] for i in I3]
    pd = [[phi if i == j else {} for j in I3] for i in I3]
    pX = [[mul(phi, mul(XS[i + 1], XS[j + 1])) for j in I3] for i in I3]
    return g, Hs, Pp, pd, pX


def build_even(phi, q, lo=None):
    lo = L if lo is None else lo
    q = [Fr(c) for c in q]
    qs = ds(q)
    beta = pscale(qs, Fr(-2, lo + 1))
    alpha = pscale(padd(q, [Fr(0)] + pscale(qs, Fr(2, lo + 1))), Fr(1, lo))
    als = ds(alpha)
    f2 = pscale(padd(beta, pscale(als, -4)), Fr(1, lo + 2))
    f3 = padd(pscale(als, -2), pscale(f2, -(lo + 1)))
    num = padd(pscale(f2, -2 * lo), pscale(f3, -3))
    assert num[0] == 0, "f4 not smooth"
    f4 = num[1:] or [Fr(0)]
    f1 = pscale(padd(alpha, [Fr(0)] + pscale(f2, -1)), Fr(1, lo - 1))
    g, Hs, Pp, pd, pX = ebasis(phi)
    return tcomb([Fr(1)] * 4, [tmulp(Hs, spoly(f1)), tmulp(Pp, spoly(f2)), tmulp(pd, spoly(f3)), tmulp(pX, spoly(f4))])


def gauge_h(phi, a, b):
    """2 sym grad xi, xi = a(s) grad phi + b(s) phi x."""
    g = [pdiff(phi, i + 1) for i in I3]
    xi = [add(mul(spoly(a), g[i]), mul(mul(spoly(b), phi), XS[i + 1])) for i in I3]
    return tscale(symgrad(xi), 2)


def gauge_K(phi, c):
    return hess(mul(spoly(c), phi))


# ------------------------------------------------------------------ odd tensors (copied from charge-id)
def KP(phi):
    g = [pdiff(phi, i + 1) for i in I3]
    v = [sumall([scale(mul(XS[j + 1], g[k]), eps(i, j, k)) for j in I3 for k in I3 if eps(i, j, k)]) for i in I3]
    K = [[scale(add(pdiff(v[j], i + 1), pdiff(v[i], j + 1)), Fr(1, 2)) for j in I3] for i in I3]
    Pt = [[add(mul(XS[i + 1], v[j]), mul(v[i], XS[j + 1])) for j in I3] for i in I3]
    return v, K, Pt


def build_odd(phi, eta, lo=None):
    lo = L if lo is None else lo
    eta = [Fr(c) for c in eta]
    g = pscale(ds(eta), Fr(-2, lo + 2))
    f = pscale(padd(eta, [Fr(0)] + pscale(g, -1)), Fr(2, lo - 1))
    v, K, Pt = KP(phi)
    return tadd(tmulp(K, spoly(f)), tmulp(Pt, spoly(g)))


def psiB_of_K(K):
    return scale(xTx(curl(K)), -1)


def psiBdot_of_h(h):
    return xTx(curl(ric_lin(h)))


# ------------------------------------------------------------------ the geometric charge at the coordinate ball
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


NQ = [0]


def Qh(h, N=NLAPSE):
    NQ[0] += 1
    return -ball(mul(N, rscal_2(h))) - 2 * d2F(h)


def BK(Ka, Kb, N=NLAPSE):
    return ball(mul(N, sub(dot(Ka, Kb), mul(trace(Ka), trace(Kb)))))


def gram_h(vs, N=NLAPSE):
    n = len(vs)
    d = [Qh(v, N) for v in vs]
    M = [[Fr(0)] * n for _ in range(n)]
    for i in range(n):
        M[i][i] = d[i]
        for j in range(i + 1, n):
            M[i][j] = M[j][i] = (Qh(tadd(vs[i], vs[j]), N) - d[i] - d[j]) / 2
    return M


def Bh(h1, h2, q1=None, q2=None):
    q1 = Qh(h1) if q1 is None else q1
    q2 = Qh(h2) if q2 is None else q2
    return (Qh(tadd(h1, h2)) - q1 - q2) / 2


def schur(Maa, Vag, Mgg):
    """A_ext on the data block: Maa - Vag Mgg^+ Vga, with consistency of every stationarity system asserted."""
    na = len(Maa)
    sols = []
    for a in range(na):
        c, rk = solve(Mgg, [-x for x in Vag[a]])
        if c is None:
            return None, None, rk
        sols.append(c)
    Ax = [[Maa[a][b] + sum(Vag[a][g] * sols[b][g] for g in range(len(Mgg))) for b in range(na)] for a in range(na)]
    return Ax, sols, rk


def matdiff(A, B):
    return max(abs(A[i][j] - B[i][j]) for i in range(len(A)) for j in range(len(A[0])))


# ================================================================== data
random.seed(11 + L)
rr = lambda n: [Fr(random.randint(-9, 9), random.randint(1, 5)) for _ in range(n)]
phi = harmonic(L, Fr(3, 5), Fr(4, 5), (1, 2, 3), "re")
phi2 = harmonic(L, Fr(5, 13), Fr(12, 13), (3, 1, 2), "im")
NY = sph(mul(phi, phi))
nb = KMAX + 1
Tk = [build_even(phi, unit(k)) for k in range(nb)]
okTT = all(is_zero_vec(div(T)) and not trace(T) and not sub(xTx(T), mul(powp(R2, k), phi)) for k, T in enumerate(Tk))
check("D0 even basis T[s^k], k=0..KMAX: div = 0, tr = 0, x.T.x = s^k phi", okTT)
g_, Hs_, Pp_, pd_, pX_ = ebasis(phi)
check("D0b T[1] = Hess(phi)/(l(l-1)) (the TT gauge mode)", is_zero_t(tadd(Tk[0], tscale(Hs_, Fr(-1, L * (L - 1))))))
if FAILS:
    stop("even TT construction failed")
log(f"l={L} c_l={CL} N_Y/4pi={NY} nb={nb}")

# C3 Green identity
def randpoly(deg, nterm=6):
    p = {}
    for _ in range(nterm):
        e = [0, random.randint(0, deg), random.randint(0, deg), random.randint(0, deg)]
        p = add(p, {tuple(e): Fr(random.randint(-5, 5), random.randint(1, 4))})
    return p
qq = tzero()
for i in I3:
    for j in range(i, 3):
        qq[i][j] = qq[j][i] = randpoly(3)
lhs = ball(mul(NLAPSE, rscal_lin(qq)))
check("C3 Green identity int N R_lin[q] = -2 dF[q] on random symmetric q", lhs == -sph(sub(trace(qq), xTx(qq))) + 2 * ball(trace(qq)), f"{lhs}")

# ================================================================== FIELD (h) BLOCK
half, zero, one = Fr(1, 2), Fr(0), Fr(1)
FH = [gauge_h(phi, [Fr(1, L)], [zero]),          # xi = grad phi / l          sigma = Y
      gauge_h(phi, [zero], [one]),               # xi = phi x                 sigma = Y
      gauge_h(phi, [zero, one], [zero])]         # xi = s grad phi            sigma = l Y
sigF = [one, one, Fr(L)]
NH = [gauge_h(phi, [half, -half], [zero]),       # xi = N grad phi            sigma = 0, interior move
      gauge_h(phi, [Fr(-1, L)], [one]),          # xi = phi x - grad phi / l  tangent on the edge
      gauge_h(phi, [zero, one, -one], [one, -one])]   # another interior move
vs = Tk + FH + NH
t1 = time.time()
Mh = gram_h(vs)
log(f"field Gram ({len(vs)} vectors, {NQ[0]} R^(2) evaluations) t={time.time()-t1:.1f}s")
iF = list(range(nb, nb + 3))
iN = list(range(nb + 3, nb + 6))
check("H1 A_un unchanged by edge-tangent and interior xi: their rows vanish against all TT data and all gauge data",
      all(Mh[i][j] == 0 for i in iN for j in range(len(vs))))
Mgg = [[Mh[i][j] for j in iF] for i in iF]
check("H1b A_un depends on xi only through sigma: gauge-family form = -(l-1)(l+2) sigma_i sigma_j oint Y^2 (rank 1)",
      all(Mgg[i][j] == -(L - 1) * (L + 2) * sigF[i] * sigF[j] * NY for i in range(3) for j in range(3)),
      f"M(sigma=Y)={Mgg[0][0]} ; flat-space value -(l-1)(l+2)N_Y={-(L-1)*(L+2)*NY}")
Maa = [[Mh[i][j] for j in range(nb)] for i in range(nb)]
Vag = [[Mh[i][j] for j in iF] for i in range(nb)]
AextH, solH, rkH = schur(Maa, Vag, Mgg)
check("Q1a(h) stationarity in sigma solvable for every TT datum, rank of the gauge form = 1", AextH is not None and rkH == 1, f"rank={rkH}")
if AextH is None:
    stop("field-block stationarity inconsistent")
psiE = [xTx(ric_lin(T)) for T in Tk]
check("D4 psi_E[T[s^k]] = -(1/2) Lap(s^k phi) (companion A8 Lemma 7)",
      all(not add(psiE[k], scale(lap(mul(powp(R2, k), phi)), half)) for k in range(nb)))


def Efield(w, N=NLAPSE):
    return [[CL * (half * ball(mul(N, gdot(psiE[i], psiE[j]))) + Fr(w, 2) * ball(mul(psiE[i], psiE[j])))
             for j in range(nb)] for i in range(nb)]


E1f, E2f = Efield(1), Efield(2)
check("Q2a FIELD block: A_ext^h = c_l E_1[psi_E] on the whole Gram (k=0..KMAX), exact", matdiff(AextH, E1f) == 0 and any(x != 0 for r in E1f for x in r),
      f"max|diff|={matdiff(AextH, E1f)} ; A_ext(T1,T1)={AextH[1][1]}")
plant("PL4 weight 2 on the even field block", matdiff(AextH, E2f) != 0, f"max|diff|={float(matdiff(AextH, E2f)):.3e}")
check("Q2a' gauge mode T[1] has a zero row in A_ext^h", all(AextH[0][j] == 0 for j in range(nb)))
# Lemma 8 closed form and Corollary 2(ii)
Gk = [scale(add(xgrad(xTx(T)), scale(xTx(T), 2)), half) for T in Tk]
L8 = [[Maa[i][j] + sph(mul(Gk[i], Gk[j])) / ((L - 1) * (L + 2)) for j in range(nb)] for i in range(nb)]
check("H3(h) companion A8 Lemma 8: A_ext^h = A_un^h + oint G_a G_b/((l-1)(l+2)), G = (1/2)(d_r + 2)(x.h.x)", matdiff(L8, AextH) == 0)
c4 = all(L * (L - 1) * Maa[0][k] == -Fr(L, 4) * (2 * k + (L + 2)) * NY for k in range(nb))
check("C4 companion A8 Cor. 2(ii): unextremised A_un(Hess phi, T[s^k]) = -(l/4)(2k + l + 2) N_Y, non-zero", c4 and Maa[0][1] != 0,
      f"A_un(Hess phi,Hess phi)={L*L*(L-1)*(L-1)*Maa[0][0]}")
check("Q1d(h) non-vacuity: A_un^h differs from A_ext^h on TT data (the extremisation is doing something)", matdiff(Maa, AextH) != 0)
# PL3 wrong lapse, reduced set
Nbad = scale(sub(const(1), scale(R2, Fr(11, 10))), half)
Mb = gram_h([Tk[1], FH[0]], Nbad)
plant("PL3 lapse (1 - 1.1 r^2)/2 in A_un: extremised value on T[s] differs from c_l E_1",
      Mb[1][1] == 0 or Mb[0][0] - Mb[0][1] ** 2 / Mb[1][1] != E1f[1][1])
dump()

# ================================================================== MOMENTUM (K) BLOCK
DK = KMAX + 1
FK = [gauge_K(phi, unit(k)) for k in range(DK + 1)]
XK = [gauge_K(phi, unit(k)) for k in (DK + 1, DK + 2)]
X0 = [mul(powp(R2, k), phi) for k in range(DK + 3)]
vk = Tk + FK + XK
nk = len(vk)
Mk = [[BK(vk[i], vk[j]) for j in range(nk)] for i in range(nk)]
jF = list(range(nb, nb + DK + 1))
jX = list(range(nb + DK + 1, nk))


def tanflat(i, j):
    gi, gj = grad(X0[i]), grad(X0[j])
    full = sumall([mul(a, b) for a, b in zip(gi, gj)])
    return sph(sub(full, mul(xgrad(X0[i]), xgrad(X0[j])))) - 2 * ball(full)


check("H2(K) pure-gauge K = Hess X0 returns flat-space geometry: A_un = oint|grad_tan tau|^2 - 2 int|grad X0|^2 (all pairs)",
      all(Mk[nb + i][nb + j] == tanflat(i, j) for i in range(DK + 3) for j in range(DK + 3)))
kap = [xTx(T) for T in Tk]
check("H3(K) drive: A_un(T_a, Hess X0) = oint kappa_a tau (depends on X0 only through the corner time shift)",
      all(Mk[a][nb + j] == sph(mul(kap[a], X0[j])) for a in range(nb) for j in range(DK + 3)))
Mgk = [[Mk[i][j] for j in jF] for i in jF]
Kaa = [[Mk[i][j] for j in range(nb)] for i in range(nb)]
Vak = [[Mk[i][j] for j in jF] for i in range(nb)]
AextK, solK, rkK = schur(Kaa, Vak, Mgk)
check("Q1a(K) stationarity in X0 solvable, gauge form non-degenerate on the polynomial family", AextK is not None and rkK == DK + 1, f"rank={rkK} of {DK+1}")
if AextK is None:
    stop("momentum-block stationarity inconsistent")
okx = all(Mk[a][x] + sum(solK[a][g] * Mk[jF[g]][x] for g in range(DK + 1)) == 0 for a in range(nb) for x in jX)
check("Q1c(K) the stationary point is stationary against directions outside the family (s^(D+1) phi, s^(D+2) phi)", okx)
okh = all(solK[a][0] == -Fr(1, L * (L - 1)) and all(c == 0 for c in solK[a][1:]) for a in range(nb))
check("H3'(K) stationary X0 is harmonic, X0 = -kappa(1)/(l(l-1)) phi (maximal slice; moved K has zero normal-normal part on the edge)", okh)
CK = [[Kaa[a][b] - sph(mul(kap[a], kap[b])) / (L * (L - 1)) for b in range(nb)] for a in range(nb)]
check("H3''(K) closed form: A_ext^K = int N K_a.K_b - oint kappa_a kappa_b/(l(l-1)) on TT K", matdiff(CK, AextK) == 0)
check("Q1d(K) non-vacuity: A_un^K differs from A_ext^K", matdiff(Kaa, AextK) != 0)

psd = [xTx(ric_lin(tscale(T, 2))) for T in Tk]
Fp = [prof(p, phi) for p in psd]
Up = [GD(F) for F in Fp]
uD = [mul(spoly(U), phi) for U in Up]
check("D4b G_D: -Lap u = psidot_E and u = 0 on the edge (u = F(s) phi with F(1) = 0)",
      all(not add(lap(uD[k]), psd[k]) and sum(Up[k]) == 0 for k in range(nb)))
SK = [[CL * half * ball(mul(NLAPSE, mul(psd[i], psd[j]))) for j in range(nb)] for i in range(nb)]
ZK = [[CL * half * ball(mul(psd[i], uD[j])) for j in range(nb)] for i in range(nb)]
SZ = [[SK[i][j] + ZK[i][j] for j in range(nb)] for i in range(nb)]
check("Q2b MOMENTUM block: A_ext^K = c_l[(1/2) int N psidot^2 + (1/2) int psidot G_D psidot] on the whole Gram, exact",
      matdiff(AextK, SZ) == 0 and any(x != 0 for r in ZK for x in r), f"max|diff|={matdiff(AextK, SZ)} ; A_ext(T1,T1)={AextK[1][1]} ; Z_E(T1,T1)={ZK[1][1]}")
check("Q2b' momentum block is NOT c_l E_w for any w (E_w has no w on this block): A_ext^K - c_l(1/2)int N psidot^2 = Z_E != 0",
      matdiff(AextK, SK) != 0, f"Z_E/S on T[s] = {ZK[1][1]/SK[1][1]}")
check("Q2b'' Z_E positive on the basis diagonal (k >= 1) and zero on the gauge mode", all(ZK[k][k] > 0 for k in range(1, nb)) and all(AextK[0][j] == 0 for j in range(nb)))
uN = [mul(spoly(GD(F, neumann=True)), phi) for F in Fp]
ZN = [[CL * half * ball(mul(psd[i], uN[j])) for j in range(nb)] for i in range(nb)]
plant("PL5a Neumann-type inverse in place of G_D", matdiff(AextK, [[SK[i][j] + ZN[i][j] for j in range(nb)] for i in range(nb)]) != 0)
plant("PL5b coefficient 2 on the G_D term", matdiff(AextK, [[SK[i][j] + 2 * ZK[i][j] for j in range(nb)] for i in range(nb)]) != 0)
# spectrum
sub3 = lambda M: [[M[i][j] for j in range(1, nb)] for i in range(1, nb)]
A3, S3 = sub3(AextK), sub3(SK)
eigs = [1 + Fr(1, (j + 1) * (2 * L + 2 * j + 3)) for j in range(nb - 1)]
okE = all(det([[A3[i][j] - lam * S3[i][j] for j in range(nb - 1)] for i in range(nb - 1)]) == 0 for lam in eigs)
d1 = det([[A3[i][j] - S3[i][j] for j in range(nb - 1)] for i in range(nb - 1)])
check("Q2d spectrum of A_ext^K against c_l(1/2)int N psidot^2 = {1 + 1/((j+1)(2l+2j+3))}, j=0..KMAX-1 (exact roots); 1 is not a root",
      okE and d1 != 0 and det(S3) != 0, f"eigenvalues {[str(e) for e in eigs]}")
OUT["even_momentum_eigs"] = [str(e) for e in eigs]
# flow level
okfl = True
for a in range(1, nb):
    Kt = tadd(Tk[a], tscale(FK[0], solK[a][0]))
    okfl &= not sub(xTx(ric_lin(tscale(tmulp(Kt, NLAPSE), 2))), add(mul(NLAPSE, psd[a]), uD[a]))
    okfl &= not sph(mul(xTx(Kt), xTx(Kt)))
check("Q2e flow of A_ext^K on the curvature: x.Ric_lin[2 N Kt].x = N psidot_E + G_D psidot_E (non-local in psidot_E)", okfl)
okloc = all(sub(xTx(ric_lin(tscale(tmulp(Tk[a], NLAPSE), 2))), mul(NLAPSE, psd[a])) for a in range(1, nb))
check("Q2e' without the corner shift the same expression is not N psidot_E either (TT representative is not the local one)", okloc)
# PL2: drop the bulk slice term from the stiffness
Mg2 = [[Mk[jF[i]][jF[j]] + 2 * ball(gdot(X0[i], X0[j])) for j in range(DK + 1)] for i in range(DK + 1)]
A2, _, _ = schur(Kaa, Vak, Mg2)
plant("PL2 (extremisation) bulk slice term -2 int|grad X0|^2 dropped from the stiffness", A2 is None or matdiff(A2, SZ) != 0)
# PL1: corner time shift frozen (X0 = N s^k phi only)
FK1 = [gauge_K(phi, padd(pscale(unit(k), half), pscale(unit(k + 1), -half))) for k in range(DK + 1)]


def aext_K_family(K, fam):
    Mg = [[BK(a, b) for b in fam] for a in fam]
    v = [BK(K, g) for g in fam]
    c, rk = solve(Mg, [-x for x in v])
    return None if c is None else BK(K, K) + sum(x * y for x, y in zip(v, c))


Kr = tcomb(rr(nb - 1), Tk[1:])
mv = tscale(FK[0], Fr(7, 3))
p1a, p1b = aext_K_family(Kr, FK1), aext_K_family(tadd(Kr, mv), FK1)
f1a, f1b = aext_K_family(Kr, FK), aext_K_family(tadd(Kr, mv), FK)
plant("PL1 (extremisation) corner time shift frozen: value changes under the harmonic X0 move, which keeps K TT",
      p1a != p1b and f1a == f1b, f"frozen: {float(p1a):.6g} -> {float(p1b):.6g} ; full: {float(f1a):.6g} -> {float(f1b):.6g}")
dump()
log(f"momentum block done t={time.time()-T0:.1f}s")

# ================================================================== Q1b: explicit gauge invariance, h and K independent
h0 = tcomb(rr(nb - 1), Tk[1:])
K0 = tcomb(rr(nb - 1), Tk[1:])
ga, gb, gc = rr(3), rr(3), rr(3)
hm = tadd(h0, gauge_h(phi, ga, gb))
Km = tadd(K0, gauge_K(phi, gc))
check("G0 moved data are not TT and the edge and corner really move (sigma != 0, tau != 0)",
      not is_zero_vec(div(hm)) and bool(trace(Km)) and (L * sum(ga) + sum(gb)) != 0 and sum(gc) != 0,
      f"sigma={L*sum(ga)+sum(gb)} Y, tau={sum(gc)} Y")
check("G1 psi_E and psidot_E are unchanged by the even gauge move", not sub(xTx(ric_lin(hm)), xTx(ric_lin(h0))) and not sub(xTx(ric_lin(Km)), xTx(ric_lin(K0))))
check("G1b moved data still satisfy the linearised constraints (R_lin[h] = 0, div K = grad tr K)",
      not rscal_lin(hm) and is_zero_vec([sub(a, b) for a, b in zip(div(Km), grad(trace(Km)))]))
qFH = [Mgg[i][i] for i in range(3)]


def aext_h(h):
    q = Qh(h)
    v = [Bh(h, FH[i], q, qFH[i]) for i in range(3)]
    c, rk = solve(Mgg, [-x for x in v])
    return q, (None if c is None else q + sum(x * y for x, y in zip(v, c)))


q0, a0 = aext_h(h0)
qm, am = aext_h(hm)
check("Q1b(h) A_ext^h(h + 2 sym grad xi) = A_ext^h(h), random xi with a, b of degree 2 (pipeline run on the non-TT datum)",
      a0 is not None and am is not None and a0 == am and a0 != 0, f"{float(a0):.8g} = {float(am):.8g}")
check("Q1d(h) the unextremised A_un^h moves under the same xi", q0 != qm, f"{float(q0):.6g} -> {float(qm):.6g}")
k0, km = aext_K_family(K0, FK), aext_K_family(Km, FK)
check("Q1b(K) A_ext^K(K + Hess X0) = A_ext^K(K), random X0 = c(s) phi of degree 2", k0 is not None and k0 == km and k0 != 0, f"{float(k0):.8g} = {float(km):.8g}")
check("Q1d(K) the unextremised A_un^K moves under the same X0", BK(K0, K0) != BK(Km, Km), f"{float(BK(K0,K0)):.6g} -> {float(BK(Km,Km)):.6g}")
# stationarity of the moved K datum against extra directions
Mg = [[BK(a, b) for b in FK] for a in FK]
v = [BK(Km, g) for g in FK]
cm, _ = solve(Mg, [-x for x in v])
Kst = tadd(Km, tcomb(cm, FK))
check("Q1c(K) moved datum: stationary against the extra directions, and the stationary slice is maximal (tr K = 0)",
      all(BK(Kst, x) == 0 for x in XK) and not trace(Kst))
psi0, pd0 = xTx(ric_lin(h0)), xTx(ric_lin(tscale(K0, 2)))
u0 = mul(spoly(GD(prof(pd0, phi))), phi)
chain = CL * (half * ball(mul(NLAPSE, add(mul(pd0, pd0), gdot(psi0, psi0)))) + half * ball(mul(psi0, psi0)) + half * ball(mul(pd0, u0)))
check("Q1/Q2 random datum, h and K independent, after the gauge move: A_ext = c_l E_1[psi_E, psidot_E] + (c_l/2)<psidot_E, G_D psidot_E>",
      am + km == chain, f"A_ext/4pi={float(am+km):.8g}")
OUT["spot"] = {"A_ext_h": str(am), "A_ext_K": str(km), "chain": str(chain)}
# mixed-m datum (two harmonics), closed forms
Tm = build_even(phi2, [Fr(0), Fr(2), Fr(-1, 3)])
hx_, Kx_ = tadd(h0, Tm), tadd(K0, Tm)
G2 = scale(add(xgrad(xTx(hx_)), scale(xTx(hx_), 2)), half)
ah = Qh(hx_) + sph(mul(G2, G2)) / ((L - 1) * (L + 2))
ak = BK(Kx_, Kx_) - sph(mul(xTx(Kx_), xTx(Kx_))) / (L * (L - 1))
pm_, pdm_ = xTx(ric_lin(Tm)), xTx(ric_lin(tscale(Tm, 2)))
um_ = mul(spoly(GD(prof(pdm_, phi2))), phi2)
px, pdx, ux = add(psi0, pm_), add(pd0, pdm_), add(u0, um_)
chx = CL * (half * ball(mul(NLAPSE, add(mul(pdx, pdx), gdot(px, px)))) + half * ball(mul(px, px)) + half * ball(mul(pdx, ux)))
check("S1 mixed-m datum (two harmonics), closed forms of H3: A_ext = chain form, exact", ah + ak == chx, f"{float(ah+ak):.8g}")
dump()
log(f"gauge tests done t={time.time()-T0:.1f}s")

# ================================================================== C1: ODD data through the same code path
nodd = 3
To = [build_odd(phi, unit(k)) for k in range(nodd)]
phiU = harmonic(L + 1, Fr(3, 5), Fr(4, 5), (1, 2, 3), "re")
Tu = build_odd(phiU, [Fr(1), Fr(-2, 3), Fr(1, 5)], L + 1)
check("C1-0 odd tensors are TT with x.T.x = 0", all(is_zero_vec(div(T)) and not trace(T) and not xTx(T) for T in To + [Tu]))
MoK = [[BK(a, b) for b in To] for a in To]
drK = all(BK(T, g) == 0 for T in To + [Tu] for g in FK + XK)
qo = [Qh(T) for T in To]
qu = Qh(Tu)
drH = all(Bh(T, FH[i], q, qFH[i]) == 0 for T, q in zip(To + [Tu], qo + [qu]) for i in range(3))
check("C1a ODD data, degree l and degree l+1 (same polynomial parity as the gauge family): zero drive against every even gauge direction, h and K slots",
      drK and drH)
MoH = [[Fr(0)] * nodd for _ in range(nodd)]
for i in range(nodd):
    MoH[i][i] = qo[i]
    for j in range(i + 1, nodd):
        MoH[i][j] = MoH[j][i] = Bh(To[i], To[j], qo[i], qo[j])
pB = [psiB_of_K(T) for T in To]
pBd = [psiBdot_of_h(T) for T in To]
E2K = [[CL * (half * ball(mul(NLAPSE, gdot(pB[i], pB[j]))) + ball(mul(pB[i], pB[j]))) for j in range(nodd)] for i in range(nodd)]
E1K = [[CL * (half * ball(mul(NLAPSE, gdot(pB[i], pB[j]))) + half * ball(mul(pB[i], pB[j]))) for j in range(nodd)] for i in range(nodd)]
E2H = [[CL * half * ball(mul(NLAPSE, mul(pBd[i], pBd[j]))) for j in range(nodd)] for i in range(nodd)]
check("C1b ODD data through the same geometric code: A_ext = A_un = c_l E_2[psi_B] on the Gram, K slot and h slot",
      matdiff(MoK, E2K) == 0 and matdiff(MoH, E2H) == 0 and MoK[0][0] != 0 and MoH[1][1] != 0, f"A(K_l)/4pi={MoK[0][0]}")
k5 = sph(dot(To[0], To[0])) / ((2 * L + 1) * (2 * L + 3))
check("C1c stationary magnetic datum: A = oint|K|^2/((2l+1)(2l+3)), c_l E_1/A = (l+1)/(l+2)", MoK[0][0] == k5 and E1K[0][0] / MoK[0][0] == Fr(L + 1, L + 2))
eo = [1 + Fr(1, (2 * j + 1) * (L + j + 1)) for j in range(nodd)]
check("INFO odd momentum spectrum against c_l E_1: 1 + 1/((2j+1)(l+j+1)); with the even one both are 1 + 2/(n(2l+n+1)), n odd / n even",
      all(det([[MoK[i][j] - lam * E1K[i][j] for j in range(nodd)] for i in range(nodd)]) == 0 for lam in eo)
      and all(eo[j] == 1 + Fr(2, (2 * j + 1) * (2 * L + 2 * j + 2)) for j in range(nodd))
      and all(eigs[j] == 1 + Fr(2, (2 * j + 2) * (2 * L + 2 * j + 3)) for j in range(nb - 1)))
if FAILS:
    stop("control C1 or an earlier line failed")

# ================================================================== C2: stationary electric datum
hel = [[scale(phi, -2) if i == j else {} for j in I3] for i in I3]
check("C2-0 h = -2 phi delta satisfies the linearised Hamiltonian constraint and is not TT", not rscal_lin(hel) and bool(trace(hel)))
qe, ae = aext_h(hel)
pe = xTx(ric_lin(hel))
e1 = CL * (half * ball(mul(NLAPSE, gdot(pe, pe))) + half * ball(mul(pe, pe)))
rec = Fr(2 * L * (L - 1), (L + 2) * (2 * L + 3)) * NY
check("C2 stationary electric datum through the full pipeline: A_ext = 2l(l-1)/((l+2)(2l+3)) N_Y = c_l E_1, psi_E = l(l-1) phi; ratio information/area = 1",
      ae == rec == e1 and not sub(pe, scale(phi, L * (L - 1))), f"A_ext/N_Y={ae/NY} ; unextremised A_un/N_Y={qe/NY}")
OUT["electric_stationary"] = {"A_ext_over_NY": str(ae / NY), "A_un_over_NY": str(qe / NY)}
if FAILS:
    stop("control C2 failed")
dump()
log(f"controls done t={time.time()-T0:.1f}s")

# ================================================================== Q3: vacuum side
# --- Q3c, Q3d on the ball (polynomial mode)
eh, eK = rr(nb - 1), rr(nb - 1)
hd, Kd = tcomb(eh, Tk[1:]), tcomb(eK, Tk[1:])
Et = ric_lin(hd)
Bt = tscale(symcurl(Kd), -1)
Edot = tscale(symcurl(Bt), -1)
Bdot = symcurl(Et)
psi, psid = xTx(Et), xTx(Edot)
check("Q3c-0 E_dot = -sym curl B equals Ric_lin[2K]; B is TT with x.B.x = 0; E is TT", is_zero_t(tadd(Edot, tscale(ric_lin(tscale(Kd, 2)), -1)))
      and is_zero_vec(div(Bt)) and not trace(Bt) and not xTx(Bt) and is_zero_vec(div(Et)) and not trace(Et) and not is_zero_t(Bt))
dE = tadd(tmulp(Edot, NLAPSE), tscale(symxcross(Bt), 2))
dB = tadd(tmulp(Bdot, NLAPSE), tscale(symxcross(Et), -2))
check("Q3c-1 x.dE.x = N psidot_E", not sub(xTx(dE), mul(NLAPSE, psid)))
dpsid = scale(xTx(symcurl(dB)), -1)
dNg = sumall([diff(mul(NLAPSE, diff(psi, i + 1)), i + 1) for i in I3])
c3 = lambda w: not sub(dpsid, sub(dNg, scale(psi, w)))
check("Q3c-2 -x.(sym curl dB).x = div(N grad psi_E) - 1*psi_E (weight 1: the conformal flow on the curvature)", c3(1))
plant("PL6a weight 2 in Q3c-2", not c3(2))
check("Q3c-3 dE, dB are transverse and traceless", is_zero_vec(div(dE)) and is_zero_vec(div(dB)) and not trace(dE) and not trace(dB))
qE = prof(psi, phi)
Fd = prof(psid, phi)
check("Q3d-1 E = T[q_E] with q_E phi = psi_E: E <-> psi_E algebraic, local both ways", is_zero_t(tadd(Et, tscale(build_even(phi, qE), -1))))
check("Q3d-2 B = odd-type TT tensor with profile eta = F/(l(l+1)), psidot_E = F phi: B <-> psidot_E algebraic, local both ways",
      is_zero_t(tadd(Bt, tscale(build_odd(phi, pscale(Fd, Fr(1, L * (L + 1)))), -1))))
check("Q3d-3 psi_E[Lap h] = Lap psi_E[h], psidot_E[Lap K] = Lap psidot_E[K]: the chain field obeys the scalar wave equation",
      all(not sub(xTx(ric_lin(tlap(T))), lap(xTx(ric_lin(T)))) for T in Tk))

# --- Fourier helpers (copied from subspace.py)
_H1 = [[Fr(1)]]


def h1(n):
    while len(_H1) <= n:
        h = _H1[-1]
        new = [Fr(0)] * (len(h) + 1)
        for i, c in enumerate(h):
            if i:
                new[i - 1] += i * c
            new[i + 1] -= c / 2
        _H1.append(new)
    return _H1[n]


def ft(p):
    if not p:
        return {}, None
    par = {(a + b + c) % 2 for (t, a, b, c) in p}
    assert len(par) == 1, "mixed parity"
    n0 = par.pop()
    out = {}
    for (t, a, b, c), v in p.items():
        sg = (-1) ** ((a + b + c - n0) // 2)
        ha, hb, hc = h1(a), h1(b), h1(c)
        for i, ci in enumerate(ha):
            if not ci:
                continue
            for j, cj in enumerate(hb):
                if not cj:
                    continue
                for k, ck in enumerate(hc):
                    if not ck:
                        continue
                    key = (0, i, j, k)
                    out[key] = out.get(key, 0) + v * sg * ci * cj * ck
    return {k: v for k, v in out.items() if v != 0}, n0


def ftt(A):
    return [[ft(A[i][j])[0] for j in I3] for i in I3]


def Iw(Q, s):
    tot = Fr(0)
    for (t, a, b, c), v in Q.items():
        if a % 2 or b % 2 or c % 2:
            continue
        D = a + b + c
        ang = Fr(dfact(a - 1) * dfact(b - 1) * dfact(c - 1), dfact(D + 1))
        m = D + 2 + s
        rad = dfact(m - 1) if s == 0 else 2 ** ((m - 1) // 2) * factorial((m - 1) // 2)
        tot += v * ang * rad
    return tot


def kdot(A, B):
    return sumall([mul(A[i][j], B[i][j]) for i in I3 for j in I3])


# --- Q3a, Q3b on R^3 (Gaussian mode)
P.GAUSS[0] = True
ng = min(nb, 3)
Tg = [build_even(phi, unit(k)) for k in range(ng)]
check("Q3-TT Gaussian-profile even tensors are TT on R^3", all(is_zero_vec(div(t)) and not trace(t) for t in Tg))
ps = [xTx(ric_lin(t)) for t in Tg]                 # psi_E of h = T
pd = [xTx(ric_lin(tscale(t, 2))) for t in Tg]      # psidot_E of K = T
FT_T = [ftt(t) for t in Tg]
FT_ps = [ft(p)[0] for p in ps]
FT_pd = [ft(p)[0] for p in pd]
R_ = range(ng)
okP = all(gauss2(dot(Tg[a], Tg[b])) == Iw(kdot(FT_T[a], FT_T[b]), 0) and gauss2(mul(ps[a], ps[b])) == Iw(mul(FT_ps[a], FT_ps[b]), 0)
          and gauss2(mul(pd[a], pd[b])) == Iw(mul(FT_pd[a], FT_pd[b]), 0) for a in R_ for b in R_)
check("Q3-B0 control: Parseval, position space = momentum space, tensors and scalars", okP)
if not okP:
    stop("Parseval control failed")
kset = {gauss2(dot(Tg[a], Tg[b])) / (CL * gauss2(mul(pd[a], ps[b]))) for a in R_ for b in R_}
check("Q3a Omega_grav = kappa c_l Omega_chain on R^3 with one constant kappa = +1: int K_a.h_b = c_l int psidot_E,a psi_E,b",
      kset == {Fr(1)}, f"kappa={sorted(str(x) for x in kset)}")
Gh = [[half * Iw(kdot(FT_T[a], FT_T[b]), 1) for b in R_] for a in R_]
Sp = [[Iw(mul(FT_ps[a], FT_ps[b]), 1) for b in R_] for a in R_]
GKm = [[2 * Iw(kdot(FT_T[a], FT_T[b]), -1) for b in R_] for a in R_]
Sd = [[Iw(mul(FT_pd[a], FT_pd[b]), -1) for b in R_] for a in R_]
Spm = [[Iw(mul(FT_ps[a], FT_ps[b]), -1) for b in R_] for a in R_]
Sdp = [[Iw(mul(FT_pd[a], FT_pd[b]), 1) for b in R_] for a in R_]
eqm = lambda A, B, c: all(A[a][b] == c * B[a][b] for a in R_ for b in R_)
nz = all(any(x != 0 for r in Mx for x in r) for Mx in (Gh, Sp, GKm, Sd))
check("Q3b-h vacuum two-point function, h-h block: int (1/2) omega h_a.h_b = c_l int omega psi_E,a psi_E,b, whole Gram",
      eqm(Gh, Sp, CL) and nz, f"{Gh[0][0]} = {CL*Sp[0][0]}")
check("Q3b-K vacuum two-point function, K-K block: int 2 K_a.K_b/omega = c_l int psidot_E,a psidot_E,b/omega, whole Gram",
      eqm(GKm, Sd, CL) and nz, f"{GKm[0][0]} = {CL*Sd[0][0]}")
plant("PL6b transfer (psi/2, 2 psidot): symplectic identity kept, two-point function broken",
      not eqm(Gh, Sp, CL / 4) and not eqm(GKm, Sd, CL * 4))
plant("PL6c omega and 1/omega exchanged on the scalar side", not eqm(Gh, Spm, CL) and not eqm(GKm, Sdp, CL))
plant("PL6d c_l scaled by (l+3)/(l+2)", not eqm(Gh, Sp, CL * Fr(L + 3, L + 2)))
P.GAUSS[0] = False

dump()
n = len(OUT["checks"])
npl = sum(1 for c in OUT["checks"].values() if c.get("plant"))
log(f"DONE l={L}: {n - npl - len(FAILS)}/{n - npl} checks pass, {npl - len(MISSED)}/{npl} plants caught, fails={FAILS}, missed={MISSED}, "
    f"R^(2) evaluations={NQ[0]}, t={time.time()-T0:.1f}s")
