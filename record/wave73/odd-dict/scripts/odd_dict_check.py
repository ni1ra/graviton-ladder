"""ODD-DICT check (wave 73, card Q08).

Odd sector l, psi = f(r) Y, psidot = p(r) Y.  q = r f, pi = N r p, r = tanh(rho/2).
Pass lines: ..\\logs\\00_passlines.log (fixed before the first run).
Usage: python odd_dict_check.py <out.json>
r side: exact rationals (3D Cartesian and radial).  rho side: Gauss-Legendre at 45 digits.
"""
import sys, json, time
from fractions import Fraction as Fr
from math import factorial
import mpmath as mp

mp.mp.dps = 45
T0 = time.time()
FAILS = []
RES = {}


def line(tag, ok, msg):
    s = "PASS" if ok else "FAIL"
    print(f"{tag:<10} {s}  {msg}", flush=True)
    if not ok:
        FAILS.append(f"{tag}: {msg}")


def plant(tag, caught, msg):
    s = "CAUGHT" if caught else "MISSED"
    print(f"{tag:<10} {s}  {msg}", flush=True)
    if not caught:
        FAILS.append(f"{tag} MISSED: {msg}")


# ---------------------------------------------------------------- Cartesian polynomials, exact
def padd(a, b, cb=1):
    r = dict(a)
    for k, v in b.items():
        nv = r.get(k, 0) + cb * v
        if nv == 0:
            r.pop(k, None)
        else:
            r[k] = nv
    return r


def pscale(a, c):
    return {k: c * v for k, v in a.items()} if c != 0 else {}


def pmul(a, b):
    r = {}
    for (a1, a2, a3), va in a.items():
        for (b1, b2, b3), vb in b.items():
            k = (a1 + b1, a2 + b2, a3 + b3)
            r[k] = r.get(k, 0) + va * vb
    return {k: v for k, v in r.items() if v != 0}


def pdiff(a, i):
    r = {}
    for k, v in a.items():
        if k[i] > 0:
            kk = list(k)
            kk[i] -= 1
            r[tuple(kk)] = r.get(tuple(kk), 0) + v * k[i]
    return r


ONE = {(0, 0, 0): Fr(1)}
XV = [{(1, 0, 0): Fr(1)}, {(0, 1, 0): Fr(1)}, {(0, 0, 1): Fr(1)}]
S2 = {(2, 0, 0): Fr(1), (0, 2, 0): Fr(1), (0, 0, 2): Fr(1)}
NLAPSE = padd(pscale(ONE, Fr(1, 2)), pscale(S2, Fr(-1, 2)))


def dfact(n):
    r = 1
    while n > 1:
        r *= n
        n -= 2
    return r


def ball_int(p):
    """int over the unit ball, in units of 4 pi."""
    t = Fr(0)
    for (a, b, c), v in p.items():
        if a % 2 or b % 2 or c % 2:
            continue
        n = a + b + c
        t += v * Fr(dfact(a - 1) * dfact(b - 1) * dfact(c - 1), dfact(n + 1) * (n + 3))
    return t


def sph_int(p):
    """int over the unit sphere, in units of 4 pi."""
    t = Fr(0)
    for (a, b, c), v in p.items():
        if a % 2 or b % 2 or c % 2:
            continue
        n = a + b + c
        t += v * Fr(dfact(a - 1) * dfact(b - 1) * dfact(c - 1), dfact(n + 1))
    return t


def harmonic(l):
    re, im = dict(ONE), {}
    for _ in range(l - 1):
        re, im = padd(pmul(re, XV[0]), pmul(im, XV[1]), -1), padd(pmul(re, XV[1]), pmul(im, XV[0]))
    return pmul(re, XV[2])


def poly_in_s(coef):
    r, pw = {}, dict(ONE)
    for c in coef:
        r = padd(r, pw, Fr(c))
        pw = pmul(pw, S2)
    return r


def grad(p):
    return [pdiff(p, i) for i in range(3)]


def lap(p):
    r = {}
    for i in range(3):
        r = padd(r, pdiff(pdiff(p, i), i))
    return r


def dot(u, v):
    r = {}
    for i in range(3):
        r = padd(r, pmul(u[i], v[i]))
    return r


# ---------------------------------------------------------------- radial polynomials, exact
def radd(a, b, cb=1):
    r = dict(a)
    for k, v in b.items():
        nv = r.get(k, 0) + cb * v
        if nv == 0:
            r.pop(k, None)
        else:
            r[k] = nv
    return r


def rmul(a, b):
    r = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            r[ka + kb] = r.get(ka + kb, 0) + va * vb
    return {k: v for k, v in r.items() if v != 0}


def rdiff(a):
    return {k - 1: v * k for k, v in a.items() if k != 0}


def rshift(a, n):
    return {k + n: v for k, v in a.items()}


def rint(a):
    return sum(v / (k + 1) for k, v in a.items())


NR = {0: Fr(1, 2), 2: Fr(-1, 2)}


def radial(l, coef):
    return {l + 2 * j: Fr(c) for j, c in enumerate(coef) if c != 0}


def reval(a, t):
    return mp.fsum(mp.mpf(v.numerator) / v.denominator * t ** k for k, v in a.items())


# ---------------------------------------------------------------- Gauss-Legendre
def gl_nodes(n):
    xs, ws = [], []
    tol = mp.mpf(10) ** (-(mp.mp.dps - 3))
    for i in range(1, n + 1):
        x = mp.cos(mp.pi * (i - mp.mpf(1) / 4) / (n + mp.mpf(1) / 2))
        for _ in range(200):
            p0, p1 = mp.mpf(1), x
            for k in range(2, n + 1):
                p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
            dp = n * (x * p1 - p0) / (x * x - 1)
            dx = p1 / dp
            x -= dx
            if abs(dx) < tol:
                break
        p0, p1 = mp.mpf(1), x
        for k in range(2, n + 1):
            p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
        dp = n * (x * p1 - p0) / (x * x - 1)
        xs.append(x)
        ws.append(2 / ((1 - x * x) * dp * dp))
    return xs, ws


PANELS = [0, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24, 32, 40, 48, 64, 80, 96, 112]
GLC = {}


def grid(n):
    if n not in GLC:
        xs, ws = gl_nodes(n)
        R, W = [], []
        for a, b in zip(PANELS[:-1], PANELS[1:]):
            a, b = mp.mpf(a), mp.mpf(b)
            for x, w in zip(xs, ws):
                R.append((a + b) / 2 + (b - a) / 2 * x)
                W.append((b - a) / 2 * w)
        GLC[n] = (R, W)
    return GLC[n]


def rho_grams(l, fr, pr, n, qpow=0, pipow=1, zpow=1):
    """rho-side Grams from the rho-side definitions.
    q = N^qpow r f ; pi = N^pipow r p ; Z weight N^zpow.  True dictionary: qpow=0, pipow=1, zpow=1."""
    R, W = grid(n)
    lam2 = l * (l + 1)
    fd = [rdiff(f) for f in fr]
    Q, DQ, PI, NN, SH = [], [], [], [], []
    for rho in R:
        t = mp.tanh(rho / 2)
        N = 1 / (1 + mp.cosh(rho))
        NN.append(N)
        SH.append(mp.sinh(rho))
        qrow, dqrow = [], []
        for f, df in zip(fr, fd):
            fv, dfv = reval(f, t), reval(df, t)
            if qpow == 0:
                qrow.append(t * fv)
                dqrow.append(N * (fv + t * dfv))          # d/drho = N d/dr
            else:                                         # q = N r f, N' = -r
                qrow.append(N * t * fv)
                dqrow.append(N * (-t * t * fv + N * (fv + t * dfv)))
        Q.append(qrow)
        DQ.append(dqrow)
        PI.append([N ** pipow * t * reval(p, t) for p in pr])
    nf, npp = len(fr), len(pr)
    mom = [[mp.fsum(W[k] * PI[k][i] * PI[k][j] for k in range(len(R))) / 2 for j in range(npp)] for i in range(npp)]
    fld = [[mp.fsum(W[k] * (DQ[k][i] * DQ[k][j] + lam2 * Q[k][i] * Q[k][j] / SH[k] ** 2) for k in range(len(R))) / 2
            for j in range(nf)] for i in range(nf)]
    sym = [[mp.fsum(W[k] * PI[k][i] * Q[k][j] for k in range(len(R))) for j in range(nf)] for i in range(npp)]
    zz = [[mp.fsum(W[k] * NN[k] ** zpow * Q[k][i] * Q[k][j] for k in range(len(R))) / 2 for j in range(nf)]
          for i in range(nf)]
    return mom, fld, sym, zz


def relerr(A, B):
    """max |A - B| / max |B| ; A mp matrix, B exact Fractions."""
    mx = max(abs(mp.mpf(b.numerator) / b.denominator) for row in B for b in row)
    d = max(abs(a - mp.mpf(b.numerator) / b.denominator) for ra, rb in zip(A, B) for a, b in zip(ra, rb))
    return d / mx


def relerr_mp(A, B):
    mx = max(abs(b) for row in B for b in row)
    return max(abs(a - b) for ra, rb in zip(A, B) for a, b in zip(ra, rb)) / mx


def sci(x):
    return mp.nstr(x, 3)


FCOEF = [[1], [0, 1], [0, 0, 1], [0, 0, 0, 1], [3, Fr(-7, 2), Fr(5, 3), Fr(2, 7), Fr(-11, 13)]]
PCOEF = [[1], [0, 1], [0, 0, 1], [0, 0, 0, 1], [1, 2, Fr(-9, 5), Fr(1, 3)]]
RHOS = ["0.3", "0.7", "1.1", "1.9", "2.7", "4.3", "6.1", "9.0"]


def run_l(l):
    out = {}
    lam2 = l * (l + 1)
    cl = Fr(4 * factorial(l - 2), factorial(l + 2))
    H = harmonic(l)
    assert lap(H) == {}, "H not harmonic"
    NY = sph_int(pmul(H, H))
    psi = [pmul(poly_in_s(c), H) for c in FCOEF]
    psd = [pmul(poly_in_s(c), H) for c in PCOEF]
    fr = [radial(l, c) for c in FCOEF]
    pr = [radial(l, c) for c in PCOEF]
    nf = len(fr)

    # ---- r side (A): genuine 3D integrals
    gp = [grad(p) for p in psi]
    momA = [[ball_int(pmul(NLAPSE, pmul(a, b))) / 2 / NY for b in psd] for a in psd]
    zA = [[ball_int(pmul(a, b)) / 2 / NY for b in psi] for a in psi]
    fldA = [[ball_int(pmul(NLAPSE, dot(gp[i], gp[j]))) / 2 / NY + zA[i][j] for j in range(nf)] for i in range(nf)]
    symA = [[ball_int(pmul(a, b)) / NY for b in psi] for a in psd]
    # ---- r side (B): radial reduction
    r2 = {2: Fr(1)}
    momB = [[rint(rmul(NR, rmul(rmul(a, b), r2))) / 2 for b in pr] for a in pr]
    zB = [[rint(rmul(rmul(a, b), r2)) / 2 for b in fr] for a in fr]
    fldB = [[rint(rmul(NR, radd(rmul(rmul(rdiff(a), rdiff(b)), r2), rmul(a, b), lam2))) / 2 + zB[i][j]
             for j, b in enumerate(fr)] for i, a in enumerate(fr)]
    symB = [[rint(rmul(rmul(a, b), r2)) for b in fr] for a in pr]
    ok0 = momA == momB and zA == zB and fldA == fldB and symA == symB
    line(f"L0  l={l}", ok0, f"3D Cartesian Grams == radial Grams exactly; N_Y/4pi = {NY}; Z[0][0] = {zA[0][0]}, E1fld[0][0] = {fldA[0][0]}")

    # ---- rho side
    mom, fld, sym, zz = rho_grams(l, fr, pr, 60)
    e1, e2, e3, e4 = relerr(mom, momA), relerr(fld, fldA), relerr(sym, symA), relerr(zz, zA)
    thr = mp.mpf("1e-30")
    line(f"L1  l={l}", e1 <= thr, f"(a1) momentum rel err {sci(e1)}")
    line(f"L2  l={l}", e2 <= thr, f"(a2) field rel err {sci(e2)}")
    line(f"L3  l={l}", e3 <= thr, f"(a3) symplectic, 25 entries, rel err {sci(e3)}")
    line(f"L4  l={l}", e4 <= thr, f"(b) Z with g = 1 rel err {sci(e4)}")
    gdev = max(abs(zz[i][j] / (mp.mpf(zA[i][j].numerator) / zA[i][j].denominator) - 1)
               for i in range(nf) for j in range(nf))
    line(f"L5  l={l}", gdev <= mp.mpf("1e-28"), f"g read back from Z: max |ratio - 1| = {sci(gdev)}")

    # full 10x10 E_2 Gram and symplectic Gram
    n10 = 2 * nf
    zero = Fr(0)
    E2r = [[(fldA[i][j] + zA[i][j]) if (i < nf and j < nf) else (momA[i - nf][j - nf] if (i >= nf and j >= nf) else zero)
            for j in range(n10)] for i in range(n10)]
    E2h = [[(fld[i][j] + zz[i][j]) if (i < nf and j < nf) else (mom[i - nf][j - nf] if (i >= nf and j >= nf) else mp.mpf(0))
            for j in range(n10)] for i in range(n10)]
    # Omega(a,b) = int(psidot_a psi_b - psidot_b psi_a); datum i<nf is (f_i,0), i>=nf is (0,p_i)
    def om(S, i, j, z):
        if i >= nf and j < nf:
            return S[i - nf][j]
        if i < nf and j >= nf:
            return -S[j - nf][i]
        return z
    Omr = [[om(symA, i, j, zero) for j in range(n10)] for i in range(n10)]
    Omh = [[om(sym, i, j, mp.mpf(0)) for j in range(n10)] for i in range(n10)]
    e6a, e6b = relerr(E2h, E2r), relerr(Omh, Omr)
    line(f"L6  l={l}", e6a <= thr and e6b <= thr, f"10x10 E_2 Gram rel err {sci(e6a)}; 10x10 symplectic Gram rel err {sci(e6b)}")

    g80 = rho_grams(l, fr, pr, 80)
    e7 = max(relerr_mp(a, b) for a, b in zip((mom, fld, sym, zz), g80))
    line(f"L7  l={l}", e7 <= mp.mpf("1e-32"), f"60 vs 80 nodes per panel: {sci(e7)}")

    # ---- (c) exact 3D: Hamiltonian field
    def flow3d(p, pd, w):
        dpsi = pmul(NLAPSE, pd) if pd else {}
        if p:
            g = grad(p)
            dv = {}
            for i in range(3):
                dv = padd(dv, pdiff(pmul(NLAPSE, g[i]), i))
            dpd = padd(dv, p, -w)
        else:
            dpd = {}
        return dpsi, dpd
    data = [(p, {}) for p in psi] + [({}, pd) for pd in psd]
    ok8 = True
    for a in range(n10):
        da = flow3d(data[a][0], data[a][1], 2)
        for b in range(n10):
            lhs = (ball_int(pmul(data[b][1], da[0])) - ball_int(pmul(da[1], data[b][0]))) / NY
            if lhs != 2 * E2r[a][b]:
                ok8 = False
    line(f"L8  l={l}", ok8, "Omega(b, X a) = 2 E_2(a,b) exactly on all 100 pairs, X = (N psidot, div(N grad psi) - 2 psi)")

    def radflow(f, w):
        # (r^2 N f')'/r^2 - N lam^2 f/r^2 - w f
        t1 = rshift(rdiff(rmul(rmul(r2, NR), rdiff(f))), -2)
        t2 = rshift(rmul(NR, f), -2)
        return radd(radd(t1, t2, -lam2), f, -w)
    ok9 = True
    dps = []
    for j, f in enumerate(fr):
        dp = radflow(f, 2)
        dps.append(dp)
        if any((k - l) % 2 or k < l for k in dp):
            ok9 = False
            continue
        mx = max(dp)
        coef = [dp.get(l + 2 * i, Fr(0)) for i in range((mx - l) // 2 + 1)]
        if flow3d(psi[j], {}, 2)[1] != pmul(poly_in_s(coef), H):
            ok9 = False
    line(f"L9  l={l}", ok9, "Cartesian div(N grad psi) - 2 psi == radial formula x Y, every f datum")

    # ---- (c) pointwise on the rho side
    def pointwise(w):
        worst, gworst, gvals = mp.mpf(0), mp.mpf(0), []
        for j, f in enumerate(fr):
            dp = radflow(f, w)
            qf = lambda x, f=f: mp.tanh(x / 2) * reval(f, mp.tanh(x / 2))
            for rs in RHOS:
                rho = mp.mpf(rs)
                t, N, sh = mp.tanh(rho / 2), 1 / (1 + mp.cosh(rho)), mp.sinh(rho)
                q = qf(rho)
                q2 = mp.diff(qf, rho, 2)
                lhs = q2 - lam2 * q / sh ** 2 - N * q
                rhs = N * t * reval(dp, t)
                sc = max(abs(q2), abs(lam2 * q / sh ** 2), abs(N * q))
                worst = max(worst, abs(lhs - rhs) / sc)
                g = (q2 - lam2 * q / sh ** 2 - rhs) / (N * q)
                gvals.append(g)
                gworst = max(gworst, abs(g - 1))
        return worst, gworst, gvals
    w10, g10, _ = pointwise(2)
    line(f"L10 l={l}", w10 <= mp.mpf("1e-28") and g10 <= mp.mpf("1e-25"),
         f"pointwise flow rel err {sci(w10)}; g read back, max |g - 1| = {sci(g10)} over {nf * len(RHOS)} points")

    # ---- (d) tensors
    v = [padd(pmul(XV[(i + 1) % 3], pdiff(H, (i + 2) % 3)), pmul(XV[(i + 2) % 3], pdiff(H, (i + 1) % 3)), -1)
         for i in range(3)]
    gv = [[pdiff(v[j], i) for j in range(3)] for i in range(3)]       # gv[i][j] = d_i v_j

    def sdiff(c):
        return [c[k] * k for k in range(1, len(c))] or [Fr(0)]

    def cadd(a, b, cb=1):
        n = max(len(a), len(b))
        return [(a[k] if k < len(a) else 0) + cb * (b[k] if k < len(b) else 0) for k in range(n)]

    def cmul(a, b):
        r = [Fr(0)] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                r[i + j] += x * y
        return r

    def Stens(eta):
        Deta = [2 * c for c in sdiff(eta)]                         # eta'/r as a polynomial in s
        g1 = [-c / (l + 2) for c in Deta]
        g2 = [c / (l - 1) for c in cadd(eta, [Fr(0)] + [c / (l + 2) for c in Deta])]
        G1, G2 = poly_in_s(g1), poly_in_s(g2)
        return [[padd(pmul(G1, padd(pmul(XV[i], v[j]), pmul(v[i], XV[j]))), pmul(G2, padd(gv[i][j], gv[j][i])))
                 for j in range(3)] for i in range(3)]

    def curl(T):                                                   # eps_ikl d_k T_lj
        return [[padd(pdiff(T[(i + 2) % 3][j], (i + 1) % 3), pdiff(T[(i + 1) % 3][j], (i + 2) % 3), -1)
                 for j in range(3)] for i in range(3)]

    def xTx(T):
        r = {}
        for i in range(3):
            for j in range(3):
                r = padd(r, pmul(pmul(XV[i], XV[j]), T[i][j]))
        return r

    def psiB(K):
        return pscale(xTx(curl(K)), Fr(-1))

    def psidotB(h):
        E = [[pscale(lap(h[i][j]), Fr(-1, 2)) for j in range(3)] for i in range(3)]
        return xTx(curl(E))

    def contract(A, B):
        r = {}
        for i in range(3):
            for j in range(3):
                r = padd(r, pmul(A[i][j], B[i][j]))
        return r

    one_s = [Fr(1), Fr(-1)]
    c3 = cmul(cmul(one_s, one_s), one_s)
    etas = [c3, cmul(c3, [Fr(0), Fr(1)]), cmul(cmul(c3, one_s), [Fr(1), Fr(2)])]
    tens = [Stens(e) for e in etas]
    ok11 = True
    for T in tens:
        tr = {}
        for i in range(3):
            tr = padd(tr, T[i][i])
        if tr != {}:
            ok11 = False
        for j in range(3):
            dv = {}
            for i in range(3):
                dv = padd(dv, pdiff(T[i][j], i))
            if dv != {}:
                ok11 = False
    line(f"L11 l={l}", ok11, "S[eta] traceless and divergence-free, three profiles")
    psK = [psiB(T) for T in tens]
    pdh = [psidotB(T) for T in tens]
    ok12 = True
    for e, a, b in zip(etas, psK, pdh):
        es = sdiff(e)
        ess = sdiff(es)
        et = cadd([Fr(0)] + [4 * c for c in ess], [(4 * l + 6) * c for c in es])
        if a != pmul(poly_in_s([lam2 * c for c in e]), H):
            ok12 = False
        if b != pmul(poly_in_s([Fr(lam2, 2) * c for c in et]), H):
            ok12 = False
    line(f"L12 l={l}", ok12, "psi_B[K=S[eta]] = l(l+1) eta H;  psidot_B[h=S[eta]] = (1/2) l(l+1) eta~ H")
    KH = [[ball_int(contract(tens[a], tens[b])) for b in range(3)] for a in range(3)]
    PP = [[ball_int(pmul(psK[a], pdh[b])) for b in range(3)] for a in range(3)]
    ok13 = all(KH[a][b] == -cl * PP[a][b] and KH[a][b] != 0 for a in range(3) for b in range(3))
    ratios = sorted({str(KH[a][b] / PP[a][b]) for a in range(3) for b in range(3) if PP[a][b] != 0})
    line(f"L13 l={l}", ok13, f"int_B K_a.h_b = -c_l int_B psi_B[K_a] psidot_B[h_b] on 3x3; ratios found {ratios}, c_l = {cl}; "
                             f"entry [0][0] = {KH[0][0]} (units 4pi)")

    # ---- plants
    _, _, _, zP = rho_grams(l, fr, pr, 60, zpow=2)
    d = relerr(zP, zA)
    plant(f"P1  l={l}", d > mp.mpf("1e-3"), f"Z with N^2: L4 deviation {sci(d)}")
    mP, _, sP, _ = rho_grams(l, fr, pr, 60, pipow=2)
    d1, d3 = relerr(mP, momA), relerr(sP, symA)
    plant(f"P2  l={l}", d1 > mp.mpf("1e-3") and d3 > mp.mpf("1e-3"), f"pi = N^2 r p: L1 deviation {sci(d1)}, L3 deviation {sci(d3)}")
    _, fP, _, zP3 = rho_grams(l, fr, pr, 60, qpow=1)
    d2, d4 = relerr(fP, fldA), relerr(zP3, zA)
    plant(f"P3  l={l}", d2 > mp.mpf("1e-3") and d4 > mp.mpf("1e-3"), f"q = N r f: L2 deviation {sci(d2)}, L4 deviation {sci(d4)}")
    _, gP, gv4 = pointwise(3)
    plant(f"P4  l={l}", gP > mp.mpf("1e-3"), f"weight 3: g read back = {mp.nstr(gv4[0], 12)} (max |g - 1| = {sci(gP)})")
    okP5 = all(KH[a][b] == cl * PP[a][b] for a in range(3) for b in range(3))
    plant(f"P5  l={l}", not okP5, "+c_l in L13 (role swap ignored)")

    out.update(dict(c_l=str(cl), NY_over_4pi=str(NY), a1=sci(e1), a2=sci(e2), a3=sci(e3), b=sci(e4), g_dev=sci(gdev),
                    E2_10x10=sci(e6a), Omega_10x10=sci(e6b), quad_60_80=sci(e7), flow_pointwise=sci(w10),
                    g_flow_dev=sci(g10), odd_bilinear_ratios=ratios,
                    Z_r=[[str(x) for x in row] for row in zA], Z_rho=[[mp.nstr(x, 35) for x in row] for row in zz],
                    plants=dict(P1=sci(d), P2=[sci(d1), sci(d3)], P3=[sci(d2), sci(d4)], P4=mp.nstr(gv4[0], 12))))
    return out


if __name__ == "__main__":
    outp = sys.argv[1] if len(sys.argv) > 1 else None
    print("ODD-DICT check; mp.dps =", mp.mp.dps)
    for l in (2, 3, 4):
        RES[f"l={l}"] = run_l(l)
        print(f"   [l={l} done, {time.time() - T0:.0f} s]", flush=True)
    print("FAILED LINES:", "none" if not FAILS else "")
    for f in FAILS:
        print("   ", f)
    RES["failed"] = FAILS
    RES["seconds"] = round(time.time() - T0, 1)
    if outp:
        with open(outp, "w", encoding="utf-8") as fh:
            json.dump(RES, fh, indent=1)
    sys.exit(1 if FAILS else 0)
