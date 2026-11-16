# CHARGE-ALL-L: the sector identity A_ext = c_l E_1 + Z for every l, both parities.
# Exact arithmetic only (sympy rational functions, Python Fractions). No nsimplify.
# Usage: python charge_all_l.py LMIN LMAX OUT.json [--nosym]
import sys, json, time
from fractions import Fraction as Fr
from itertools import product

LMIN, LMAX = int(sys.argv[1]), int(sys.argv[2])
OUTF = sys.argv[3]
NOSYM = "--nosym" in sys.argv
T0 = time.time()
RES = {"lines": [], "fails": [], "plants_missed": []}


def line(tag, ok, text=""):
    s = "PASS" if ok else "FAIL"
    print(f"{s} {tag} {text}", flush=True)
    RES["lines"].append([tag, s, text])
    if not ok:
        RES["fails"].append(tag)


def plant(tag, caught, text=""):
    s = "CAUGHT" if caught else "MISSED"
    print(f"PLANT {s} {tag} {text}", flush=True)
    RES["lines"].append([tag, s, text])
    if not caught:
        RES["plants_missed"].append(tag)


# =====================================================================
# PART S: symbolic in l
# =====================================================================
def symbolic_part():
    import sympy as sp
    l, r, a, b = sp.symbols("l r a b", positive=True)
    R = sp.Rational
    N = (1 - r**2) / 2
    cl = 4 / ((l - 1) * l * (l + 1) * (l + 2))
    Ra, Rb = sp.symbols("Ra Rb", positive=True)     # stand for r^(2a), r^(2b)
    d = lambda f: sp.diff(f, r) + (2 * a / r) * Ra * sp.diff(f, Ra) + (2 * b / r) * Rb * sp.diff(f, Rb)
    D = lambda f: d(f) / r
    til = lambda q: d(d(q)) + (2 * l + 2) * d(q) / r

    def evenf(q):
        be = -D(q) / (l + 1)
        al = (q + r**2 * D(q) / (l + 1)) / l
        f2 = (be - 2 * D(al)) / (l + 2)
        f3 = -D(al) - (l + 1) * f2
        f4 = -(2 * l * f2 + 3 * f3) / r**2
        f1 = (al - r**2 * f2) / (l - 1)
        return [f1, f2, f3, f4], al, be

    def oddg(eta):
        g1 = -D(eta) / (l + 2)
        g2 = (eta + r**2 * D(eta) / (l + 2)) / (l - 1)
        return [g1, g2]

    GE = {(0, 0): l * (l - 1) * (2 * l + 1) * (2 * l - 1), (0, 1): 2 * l * (l - 1) * (2 * l + 1), (0, 2): 0,
          (0, 3): l * (l - 1), (1, 1): 2 * l * (2 * l + 1) + 2 * l**2, (1, 2): 2 * l, (1, 3): 2 * l,
          (2, 2): 3, (2, 3): 1, (3, 3): 1}
    oe = (-2, 0, 0, 2)

    def Bev(F1, F2):   # oint_{S_r} T1.T2 = N_Y r^(2l) Bev
        tot = 0
        for i in range(4):
            for j in range(4):
                g = GE[(min(i, j), max(i, j))]
                tot += g * r**(oe[i] + oe[j]) * F1[i] * F2[j]
        return tot

    def Bod(g, h):     # oint_{S_r} S1.S2 = N_Y r^(2l) Bod   (N_v = l(l+1) N_Y)
        return l * (l + 1) * (2 * r**2 * g[0] * h[0] + 2 * (l - 1) * (g[0] * h[1] + g[1] * h[0])
                              + 2 * (l - 1) * (2 * l + 1) * g[1] * h[1] / r**2)

    def I01(expr):     # exact integral over (0,1) of r^(2l) * expr, expr a Laurent polynomial in r, Ra = r^(2a), Rb = r^(2b)
        # (run 01 bug: the first integrator split coefficient and power after sympy's expand, which had
        #  multiplied r-powers into l-denominators; the split failed silently. Now: one fraction, Poly in r.)
        num, den = sp.fraction(sp.together(expr))
        dp = sp.Poly(den, r)
        assert len(dp.terms()) == 1, den
        (m,), dc = dp.terms()[0]
        assert not dc.has(Ra) and not dc.has(Rb)
        tot = 0
        for (k, i, j), c in sp.Poly(num, r, Ra, Rb).terms():
            assert not c.has(r)
            tot += c / dc / (2 * l + k - m + 2 * a * i + 2 * b * j + 1)
        return tot

    at1 = lambda e: e.subs({r: 1, Ra: 1, Rb: 1})
    zero = lambda e: sp.cancel(sp.together(e)) == 0

    qf = sp.Function("q")(r)
    QS = sp.symbols("Q0:7")

    def jet(expr):
        expr = expr.doit()
        for k in range(6, 0, -1):
            expr = expr.subs(sp.Derivative(qf, (r, k)), QS[k])
        expr = expr.subs(qf, QS[0])
        assert not expr.has(sp.Derivative) and not expr.has(qf)
        return expr

    jz = lambda e: zero(jet(e))

    # ---- S01
    f, al, be = evenf(qf)
    rows = [2 * l * f[1] + 3 * f[2] + r**2 * f[3],
            (l - 1) * d(f[0]) / r + r * d(f[1]) + (l + 3) * f[1] + f[2],
            (l * d(f[1]) + d(f[2])) / r + r * d(f[3]) + (l + 4) * f[3],
            l * al + r**2 * be - qf]
    g = oddg(qf)
    rowo = (l - 1) * d(g[1]) / r + r * d(g[0]) + (l + 4) * g[0]
    line("S01", all(jz(x) for x in rows) and jz(rowo), "radial TT rows, even (3 rows + x.T.x) and odd (1 row)")

    # ---- S02
    fh_direct = [d(d(f[0])) + (2 * l - 2) * d(f[0]) / r + 4 * f[1],
                 d(d(f[1])) + (2 * l + 2) * d(f[1]) / r + 2 * f[3],
                 d(d(f[2])) + (2 * l + 2) * d(f[2]) / r + 2 * f[3],
                 d(d(f[3])) + (2 * l + 6) * d(f[3]) / r]
    fh = evenf(til(qf))[0]
    # odd: S = g1 P + g2 M, M = grad v + grad v^T (degree l-1, harmonic), Lap P = 2 M
    gh_direct = [d(d(g[0])) + (2 * l + 4) * d(g[0]) / r,
                 d(d(g[1])) + (2 * l) * d(g[1]) / r + 2 * g[0]]
    gh = oddg(til(qf))
    line("S02", all(jz(x - y) for x, y in zip(fh_direct, fh)) and all(jz(x - y) for x, y in zip(gh_direct, gh)),
         "Lap T[q] = T[qt], Lap S[eta] = S[etat]")

    # ---- S03 (b) pointwise
    q1 = d(qf)
    qt = til(qf)
    U = ((2 * l + 1) * qf**2 / r**3 + 2 * qf * q1 / r**2 + 2 * q1**2 / (r * (l + 1) * (l + 2))) / (l * (l - 1))
    line("S03", jz(Bev(f, f) - cl / 2 * qt**2 - ((2 * l + 2) * U / r + d(U))), "(b) pointwise potential U, q arbitrary")

    # ---- monomial Grams
    qa, qb = Ra, Rb
    w = 1          # the factor r^(2l) is supplied by I01

    def gram_b(qA, qB, corner=None):
        if corner is None:
            corner = 1 / (l * (l - 1))
        fA, fB = evenf(qA)[0], evenf(qB)[0]
        lhs = I01(N * Bev(fA, fB) * w * r**2) - at1(qA * qB) * corner
        rhs = cl / 2 * (I01(N * til(qA) * til(qB) * w * r**2) + I01(d(qA) * d(qB) * w * r**2))
        return lhs, rhs

    def gradgrad(FA, FB):   # oint grad(FA phi).grad(FB phi) / N_Y, times r^2 measure included by caller
        return (d(FA) * d(FB) * r**2 + l * r * d(FA * FB) + l * (2 * l + 1) * FA * FB) / r**2

    def gram_c(qA, qB, t=1, cfac=1):
        fA, alA, beA = evenf(qA)
        fB, alB, beB = evenf(qB)
        ftA, ftB = evenf(til(qA))[0], evenf(til(qB))[0]
        bulk = R(1, 4) * (-I01(N * (Bev(fA, ftB) + Bev(fB, ftA)) / 2 * w * r**2) - R(3, 2) * I01(Bev(fA, fB) * w * r**2))
        bd = at1(Bev(fA, fB) / 8 - (l * (2 * l + 1) * alA * alB + l * (alA * beB + alB * beA) + beA * beB) / 2 + qA * qB / 4)
        GA = at1(d(qA) + (l + 2) * qA) / 2
        GB = at1(d(qB) + (l + 2) * qB) / 2
        edge = GA * GB / ((l - 1) * (l + 2)) * (2 * t - t**2)
        FA, FB = -til(qA) / 2, -til(qB) / 2
        rhs = cfac * cl / 2 * (I01(N * gradgrad(FA, FB) * r**2) + I01(FA * FB * w * r**2))
        return bulk + bd + edge, rhs, (bulk + bd, GA, GB)

    def gram_oK(eA, eB, wt=2):
        lhs = I01(N * Bod(oddg(eA), oddg(eB)) * w * r**2)
        k = l * (l + 1)
        rhs = cl * k**2 * (R(1, 2) * I01(N * gradgrad(eA, eB) * r**2) + R(wt, 2) * I01(eA * eB * w * r**2))
        return lhs, rhs

    def gram_oh(eA, eB):
        gA, gB = oddg(eA), oddg(eB)
        gtA, gtB = oddg(til(eA)), oddg(til(eB))
        bulk = R(1, 4) * (-I01(N * (Bod(gA, gtB) + Bod(gB, gtA)) / 2 * w * r**2) - R(3, 2) * I01(Bod(gA, gB) * w * r**2))
        bd = at1(Bod(gA, gB) / 8 - l * (l + 1) * eA * eB / 2)
        k = l * (l + 1)
        rhs = cl * k**2 / 4 * R(1, 2) * I01(N * til(eA) * til(eB) * w * r**2)
        return bulk + bd, rhs

    lb, rb = gram_b(qa, qb)
    line("S04", zero(lb - rb), "(b) momentum Gram on r^(2a), r^(2b), symbolic (l,a,b)")

    # ---- S05 (c) potential
    kb = 2 / ((l + 2) * (l + 1))
    p = 2 * l + 1
    qt1 = d(qt)
    W = (-(3 * p * qf**2 / r**4 + 6 * qf * q1 / r**3 + 3 * kb * q1**2 / r**2 + 2 * qf * qt / r**2 + 2 * kb * q1 * qt / r + kb * qt**2)
         - 2 * N * (p * qf * qt / r**4 + (qf * qt1 + q1 * qt) / r**3 + kb * q1 * qt1 / r**2 + kb * (l * qt**2 / r**2 + qt * qt1 / r))) / (8 * l * (l - 1))
    dA = -R(1, 4) * N * Bev(f, fh) - R(3, 8) * Bev(f, f)
    F = -qt / 2
    dE = cl / 2 * (N * (d(F) + l * F / r)**2 + l * (l + 1) * N * F**2 / r**2 + F**2)
    bb = Bev(f, f) / 8 - (l * (2 * l + 1) * al**2 + 2 * l * al * be + be**2) / 2 + qf**2 / 4 + (q1 + (l + 2) * qf)**2 / (4 * (l - 1) * (l + 2))
    okW = jz(dA - dE - ((2 * l + 3) * W + r * d(W)))
    okb = zero(jet(W + bb).subs(r, 1))
    line("S05", okW and okb, f"(c) potential W: bulk {okW}, boundary W(1) = -b {okb}, q arbitrary")

    lc, rc, _ = gram_c(qa, qb)
    line("S06", zero(lc - rc), "(c) field Gram on r^(2a), r^(2b), symbolic (l,a,b)")
    lk, rk = gram_oK(qa, qb)
    line("S07", zero(lk - rk), "(a) odd K Gram, weight 2, symbolic (l,a,b)")
    lh, rh = gram_oh(qa, qb)
    line("S08", zero(lh - rh), "(a) odd h Gram, symbolic (l,a,b)")

    # ---- S09 odd potentials, eta arbitrary (eta = qf)
    k = l * (l + 1)
    densL = N * Bod(g, g) * r**2
    densR = cl * k**2 * (R(1, 2) * N * (q1**2 * r**2 + 2 * l * r * qf * q1 + l * (2 * l + 1) * qf**2) + qf**2 * r**2)
    pot = cl * k**2 * N * r * qf**2      # times r^(2l)
    line("S09a", jz(densL - densR - (2 * l * pot / r + d(pot))), "(a) odd K: difference = c_l l^2(l+1)^2 (N r^(2l+1) eta^2)'")
    Bo = Bod(g, g)
    mine = R(1, 4) * (-N * Bod(g, gh) - R(3, 2) * Bo) * r**2 - cl * k**2 / 8 * N * qt**2 * r**2      # times r^(2l)
    V = k * qf * (R(1, 4) * ((2 * l + 1) - (2 * l - 1) * r**2) * qf + R(1, 2) * r * (1 - r**2) * q1) / r   # times r^(2l)
    P = V - R(1, 8) * N * r**2 * (2 * l * Bo / r + d(Bo)) - R(1, 8) * r**3 * Bo
    okV = jz(mine - (2 * l * P / r + d(P)))
    okV1 = zero(jet(V - k * qf**2 / 2).subs(r, 1))
    line("S09b", okV and okV1, f"(a) odd h: companion A8 potential V at n = 3: bulk {okV}, V(1) = eta(1)^2 N_v/2 {okV1}")

    # ---- S10, S11 spectra
    def sdet(M):   # explicit 3x3 cofactor expansion (sympy's .det() on rational functions was too slow)
        M = M.applyfunc(sp.cancel)
        return (M[0, 0] * (M[1, 1] * M[2, 2] - M[1, 2] * M[2, 1]) - M[0, 1] * (M[1, 0] * M[2, 2] - M[1, 2] * M[2, 0])
                + M[0, 2] * (M[1, 0] * M[2, 1] - M[1, 1] * M[2, 0]))

    def gramnum(fn, exps):
        A = sp.zeros(3, 3)
        V_ = sp.zeros(3, 3)
        for i, ea in enumerate(exps):
            for j, eb in enumerate(exps):
                x, y = fn(sp.Integer(1) if ea == 0 else r**(2 * ea), sp.Integer(1) if eb == 0 else r**(2 * eb))
                A[i, j] = sp.cancel(x)
                V_[i, j] = sp.cancel(y)
        return A, V_

    def b_pair(qA, qB):
        lhs, _ = gram_b(qA, qB)
        vac = cl / 2 * I01(N * til(qA) * til(qB) * w * r**2)
        return lhs, vac

    A, Vv = gramnum(b_pair, [1, 2, 3])
    ok = all(zero(sdet(A - (1 + 1 / ((j + 1) * (2 * l + 2 * j + 3))) * Vv)) for j in range(3))
    ok = ok and not zero(sdet(Vv))
    line("S10", ok, "(b) spectrum 1 + 1/((j+1)(2l+2j+3)), j = 0,1,2, symbolic l")

    def k_pair(eA, eB):
        lhs, _ = gram_oK(eA, eB)
        _, vac = gram_oK(eA, eB, wt=1)
        return lhs, vac

    A, Vv = gramnum(k_pair, [0, 1, 2])
    ok = all(zero(sdet(A - (1 + 1 / ((2 * j + 1) * (l + j + 1))) * Vv)) for j in range(3))
    ok = ok and not zero(sdet(Vv))
    line("S11", ok, "(a) spectrum 1 + 1/((2k+1)(l+k+1)), k = 0,1,2, symbolic l")

    # ---- S12 static datum
    qs = -l * (l - 1) * r**2 / (2 * l + 3)
    lc_, rc_, (aun, Gs, _) = gram_c(qs, qs)
    val = 2 * l * (l - 1) / ((l + 2) * (2 * l + 3))
    Bc = (7 * l + 3) / ((2 * l + 3) * (l + 1))
    Ac = -(1 + 2 * Bc) / (4 * l + 2)
    g1 = l * Ac + Bc                       # xi.nu at r = 1 (xi = Ac r^2 grad phi + Bc phi x)
    qchk = -2 + 2 * (l + 1) * g1           # coefficient of r^2 in q* from the gauge vector
    sig = -Gs / ((l - 1) * (l + 2))
    shape = sp.cancel(sig + g1)
    unext = sp.cancel(val - (l - 1) * (l + 2) * shape**2)
    ok = zero(lc_ - val) and zero(rc_ - val) and zero(shape - 2 / (l + 2)) and zero(qchk + l * (l - 1) / (2 * l + 3))
    line("S12", ok, f"static datum: A_ext/N_Y = {sp.factor(sp.cancel(lc_))}, shape = {shape}, unextremised = {sp.factor(unext)}; at l=2: {sp.cancel(lc_).subs(l, 2)}, {shape.subs(l, 2)}, {unext.subs(l, 2)}")

    # ---- S13 null rows
    one = sp.Integer(1)
    x1, y1 = gram_b(one, qb)
    x2, y2, _ = gram_c(one, qb)
    x3, y3 = gram_oh(one, qb)
    line("S13", all(zero(t) for t in (x1, y1, x2, y2, x3, y3)), "gauge null rows q = 1 (b),(c) and eta = 1 (odd h), all b")

    # ---- plants
    t = (l + 2) / (l + 1)
    lcp, rcp, _ = gram_c(qa, qb, t=t)
    plant("SP1", not zero(lcp - rcp), "(c) shape 2/(l+1) for 2/(l+2); on the static datum at l=2 A/N_Y = %s (true 1/7)" % sp.cancel(gram_c(qs, qs, t=t)[0]).subs(l, 2))
    lbp, rbp = gram_b(qa, qb, corner=1 / (l * (l + 1)))
    plant("SP2", not zero(lbp - rbp), "(b) corner coefficient 1/(l(l+1))")
    lkp, rkp = gram_oK(qa, qb, wt=1)
    plant("SP3", not zero(lkp - rkp), "(a) weight 1 on the odd K block")
    lcq, rcq, _ = gram_c(qa, qb, cfac=(l + 3) / (l + 2))
    plant("SP4", not zero(lcq - rcq), "c_l scaled by (l+3)/(l+2) in (c)")


# =====================================================================
# PART C: Cartesian exact polynomial arithmetic
# =====================================================================
I3 = range(3)


def padd(*ps):
    out = {}
    for p_ in ps:
        for m, c in p_.items():
            v = out.get(m, 0) + c
            if v:
                out[m] = v
            elif m in out:
                del out[m]
    return out


def psc(p_, c):
    return {m: v * c for m, v in p_.items()} if c else {}


def psub(p_, q_):
    return padd(p_, psc(q_, -1))


def pmul(p_, q_):
    out = {}
    if len(p_) > len(q_):
        p_, q_ = q_, p_
    for (a1, b1, c1), u in p_.items():
        for (a2, b2, c2), v in q_.items():
            m = (a1 + a2, b1 + b2, c1 + c2)
            out[m] = out.get(m, 0) + u * v
    return {m: v for m, v in out.items() if v}


def pdiff(p_, i):
    out = {}
    for m, c in p_.items():
        if m[i]:
            mm = list(m)
            mm[i] -= 1
            out[tuple(mm)] = c * m[i]
    return out


ONE = {(0, 0, 0): Fr(1)}
XV = [{(1, 0, 0): Fr(1)}, {(0, 1, 0): Fr(1)}, {(0, 0, 1): Fr(1)}]
SS = padd(*[pmul(x, x) for x in XV])


def dfact(n):
    out = 1
    while n > 1:
        out *= n
        n -= 2
    return out


_W = {}


def wsph(m):
    if m not in _W:
        if any(e % 2 for e in m):
            _W[m] = Fr(0)
        else:
            _W[m] = Fr(dfact(m[0] - 1) * dfact(m[1] - 1) * dfact(m[2] - 1), dfact(sum(m) + 1))
    return _W[m]


def sph(p_):
    return sum((c * wsph(m) for m, c in p_.items()), Fr(0))


def ball(p_):
    return sum((c * wsph(m) / (sum(m) + 3) for m, c in p_.items()), Fr(0))


def ballN(p_):
    return sum((c * wsph(m) / ((sum(m) + 3) * (sum(m) + 5)) for m, c in p_.items()), Fr(0))


def lap(p_):
    return padd(*[pdiff(pdiff(p_, i), i) for i in I3])


def grad(p_):
    return [pdiff(p_, i) for i in I3]


def xdot(v):
    return padd(*[pmul(XV[i], v[i]) for i in I3])


def dotv(u, v):
    return padd(*[pmul(u[i], v[i]) for i in I3])


def tdot(A, B):
    return padd(*[pmul(A[i][j], B[i][j]) for i in I3 for j in I3])


def tadd(A, B, c=1):
    return [[padd(A[i][j], psc(B[i][j], c)) for j in I3] for i in I3]


def tzero(A):
    return all(not A[i][j] for i in I3 for j in I3)


def ttr(A):
    return padd(*[A[i][i] for i in I3])


def tdiv(A):
    return [padd(*[pdiff(A[i][j], i) for i in I3]) for j in I3]


def tx(A):
    return [padd(*[pmul(A[i][j], XV[j]) for j in I3]) for i in I3]


def xtx(A):
    return xdot(tx(A))


def dTdT(A, B):
    return padd(*[pmul(pdiff(A[i][j], k), pdiff(B[i][j], k)) for i in I3 for j in I3 for k in I3])


def symgrad2(v):   # d_i v_j + d_j v_i
    return [[padd(pdiff(v[j], i), pdiff(v[i], j)) for j in I3] for i in I3]


def ricci(h):
    tr = ttr(h)
    dv = tdiv(h)   # dv_j = d_k h_kj
    return [[psc(padd(pdiff(dv[j], i), pdiff(dv[i], j), psc(lap(h[i][j]), -1), psc(pdiff(pdiff(tr, i), j), -1)), Fr(1, 2))
             for j in I3] for i in I3]


EPS = {(0, 1, 2): 1, (1, 2, 0): 1, (2, 0, 1): 1, (0, 2, 1): -1, (2, 1, 0): -1, (1, 0, 2): -1}


def curl(T):   # eps_ikl d_k T_lj
    return [[padd(*[psc(pdiff(T[l_][j], k), e) for (i2, k, l_), e in EPS.items() if i2 == i]) for j in I3] for i in I3]


def cross(u, v):
    return [padd(*[psc(pmul(u[j], v[k]), e) for (i2, j, k), e in EPS.items() if i2 == i]) for i in I3]


# radial polynomials in s = r^2: lists of Fractions
def rD(q):   # (1/r) d/dr = 2 d/ds
    return [2 * k * q[k] for k in range(1, len(q))] or [Fr(0)]


def radd(p_, q_, c=1):
    n = max(len(p_), len(q_))
    return [(p_[i] if i < len(p_) else 0) + c * (q_[i] if i < len(q_) else 0) for i in range(n)]


def rsc(p_, c):
    return [c * v for v in p_]


def rs(p_):   # multiply by s
    return [Fr(0)] + list(p_)


def rdivs(p_):
    assert p_[0] == 0, "not divisible by s"
    return list(p_[1:]) or [Fr(0)]


def rtil(q, L):   # q'' + (2l+2) q'/r = 4 s q_ss + (4l+6) q_s
    qs_ = [k * q[k] for k in range(1, len(q))] or [Fr(0)]
    qss = [k * qs_[k] for k in range(1, len(qs_))] or [Fr(0)]
    return radd(rsc(rs(qss), 4), rsc(qs_, 4 * L + 6))


def r1(q):
    return sum(q, Fr(0))


def rcart(q):
    out, pw = {}, ONE
    for c in q:
        out = padd(out, psc(pw, c))
        pw = pmul(pw, SS)
    return out


def harmonic(L):
    re, im = ONE, {}
    for _ in range(L - 1):
        re, im = psub(pmul(re, XV[0]), pmul(im, XV[1])), padd(pmul(re, XV[1]), pmul(im, XV[0]))
    return pmul(re, XV[2])


def build_even(phi, q, L):
    be = rsc(rD(q), Fr(-1, L + 1))
    al = rsc(radd(q, rsc(rs(rD(q)), Fr(1, L + 1))), Fr(1, L))
    f2 = rsc(radd(be, rsc(rD(al), -2)), Fr(1, L + 2))
    f3 = radd(rsc(rD(al), -1), rsc(f2, -(L + 1)))
    f4 = rsc(rdivs(radd(rsc(f2, 2 * L), rsc(f3, 3))), -1)
    f1 = rsc(radd(al, rsc(rs(f2), -1)), Fr(1, L - 1))
    F = [rcart(x) for x in (f1, f2, f3, f4)]
    g = grad(phi)
    H = [[pdiff(g[i], j) for j in I3] for i in I3]
    T = [[padd(pmul(F[0], H[i][j]),
               pmul(F[1], padd(pmul(XV[i], g[j]), pmul(g[i], XV[j]))),
               pmul(F[2], phi) if i == j else {},
               pmul(F[3], pmul(phi, pmul(XV[i], XV[j])))) for j in I3] for i in I3]
    return T


def build_odd(phi, eta, L):
    g1 = rsc(rD(eta), Fr(-1, L + 2))
    g2 = rsc(radd(eta, rsc(rs(rD(eta)), Fr(1, L + 2))), Fr(1, L - 1))
    v = cross(XV, grad(phi))
    G1, G2 = rcart(g1), rcart(g2)
    M = symgrad2(v)
    S = [[padd(pmul(G1, padd(pmul(XV[i], v[j]), pmul(v[i], XV[j]))), pmul(G2, M[i][j])) for j in I3] for i in I3]
    return S, v


def basis_gram(phi):
    g = grad(phi)
    H = [[pdiff(g[i], j) for j in I3] for i in I3]
    P = [[padd(pmul(XV[i], g[j]), pmul(g[i], XV[j])) for j in I3] for i in I3]
    Dl = [[phi if i == j else {} for j in I3] for i in I3]
    Xx = [[pmul(phi, pmul(XV[i], XV[j])) for j in I3] for i in I3]
    Bs = [H, P, Dl, Xx]
    return {(i, j): sph(tdot(Bs[i], Bs[j])) for i in range(4) for j in range(i, 4)}


def det3(M):
    return (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
            + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))


def R2scalar(h):
    dh = [[[pdiff(h[j][k], i) for k in I3] for j in I3] for i in I3]      # dh[i][j][k] = d_i h_jk
    G1 = [[[psc(padd(dh[i][j][k], dh[j][i][k], psc(dh[k][i][j], -1)), Fr(1, 2)) for j in I3] for i in I3] for k in I3]  # G1[k][i][j]
    G2 = [[[psc(padd(*[pmul(h[k][m], G1[m][i][j]) for m in I3]), -1) for j in I3] for i in I3] for k in I3]
    tot = {}
    for i in I3:
        for j in I3:
            ric1 = padd(*[pdiff(G1[k][i][j], k) for k in I3], psc(pdiff(padd(*[G1[k][i][k] for k in I3]), j), -1))
            tot = padd(tot, psc(pmul(h[i][j], ric1), -1))
    for i in I3:
        tot = padd(tot, *[pdiff(G2[k][i][i], k) for k in I3], psc(pdiff(padd(*[G2[k][i][k] for k in I3]), i), -1))
        for k in I3:
            for m in I3:
                tot = padd(tot, pmul(G1[k][k][m], G1[m][i][i]), psc(pmul(G1[k][i][m], G1[m][i][k]), -1))
    return tot


def Aun_chr(h):
    tr = ttr(h)
    hx = tx(h)
    nhn = xdot(hx)
    hh = tdot(h, h)
    trt = psub(tr, nhn)
    htan2 = padd(hh, psc(dotv(hx, hx), -2), pmul(nhn, nhn))
    d2A = sph(padd(psc(pmul(trt, trt), Fr(1, 8)), psc(htan2, Fr(-1, 4))))
    d2V = ball(padd(psc(pmul(tr, tr), Fr(1, 8)), psc(hh, Fr(-1, 4))))
    return -ballN(R2scalar(h)) - 2 * (d2A - 2 * d2V)


def A_TT(h1, h2):
    x1, x2 = tx(h1), tx(h2)
    return Fr(1, 4) * ballN(dTdT(h1, h2)) - Fr(1, 2) * sph(dotv(x1, x2)) + Fr(1, 4) * sph(pmul(xdot(x1), xdot(x2)))


def cart_part(L):
    t0 = time.time()
    phi = harmonic(L)
    assert not lap(phi)
    NY = sph(pmul(phi, phi))
    CL = Fr(4, (L - 1) * L * (L + 1) * (L + 2))
    profs = [[Fr(1)], [Fr(0), Fr(1)], [Fr(0), Fr(0), Fr(1)], [Fr(0), Fr(0), Fr(0), Fr(1)]]
    nb = len(profs)
    tag = lambda s: f"{s}[l={L}]"

    # ----- even
    Ts = [build_even(phi, q, L) for q in profs]
    ok = all(not ttr(T) and all(not c for c in tdiv(T)) and not psub(xtx(T), pmul(rcart(q), phi)) for T, q in zip(Ts, profs))
    line(tag("C01"), ok, "T[q] TT, x.T.x = q phi")
    psi = [xtx(ricci(T)) for T in Ts]
    ok = all(not padd(ps, psc(pmul(rcart(rtil(q, L)), phi), Fr(1, 2))) for ps, q in zip(psi, profs))
    line(tag("C02"), ok, "x.Ric[T].x = -(1/2) qt phi")
    bg = basis_gram(phi)
    form = {(0, 0): L * (L - 1) * (2 * L + 1) * (2 * L - 1), (0, 1): 2 * L * (L - 1) * (2 * L + 1), (0, 2): 0, (0, 3): L * (L - 1),
            (1, 1): 2 * L * (2 * L + 1) + 2 * L * L, (1, 2): 2 * L, (1, 3): 2 * L, (2, 2): 3, (2, 3): 1, (3, 3): 1}
    line(tag("C03"), all(bg[k] == form[k] * NY for k in form), "sphere Gram of (Hess phi, P_phi, phi delta, phi X)")

    kap = [xtx(T) for T in Ts]
    Gp = [psc(padd(xdot(grad(k_)), psc(k_, 2)), Fr(1, 2)) for k_ in kap]
    Aun = [[A_TT(Ts[i], Ts[j]) for j in range(nb)] for i in range(nb)]
    GG = [[sph(pmul(Gp[i], Gp[j])) for j in range(nb)] for i in range(nb)]
    gpsi = [grad(p_) for p_ in psi]
    E1f = [[CL * (Fr(1, 2) * ballN(dotv(gpsi[i], gpsi[j])) + Fr(1, 2) * ball(pmul(psi[i], psi[j]))) for j in range(nb)] for i in range(nb)]
    K_ = (L - 1) * (L + 2)
    dif = [[Aun[i][j] + GG[i][j] / K_ - E1f[i][j] for j in range(nb)] for i in range(nb)]
    ok = all(v == 0 for row in dif for v in row) and all(Aun[0][j] + GG[0][j] / K_ == 0 for j in range(nb))
    line(tag("C04"), ok, f"(c) field Gram {nb}x{nb}; A_ext[s,s]/N_Y = {(Aun[1][1] + GG[1][1] / K_) / NY}")

    psid = [psc(p_, 2) for p_ in psi]
    us = [pmul(rcart(radd(q, [r1(q)], -1)), phi) for q in profs]
    oku = all(not padd(lap(u), pd) and sph(pmul(u, u)) == 0 for u, pd in zip(us, psid))
    gu = [grad(u) for u in us]
    Am = [[ballN(tdot(Ts[i], Ts[j])) - sph(pmul(kap[i], kap[j])) / (L * (L - 1)) for j in range(nb)] for i in range(nb)]
    Vm = [[CL * Fr(1, 2) * ballN(pmul(psid[i], psid[j])) for j in range(nb)] for i in range(nb)]
    Zm = [[CL * Fr(1, 2) * ball(dotv(gu[i], gu[j])) for j in range(nb)] for i in range(nb)]
    ok = oku and all(Am[i][j] - Vm[i][j] - Zm[i][j] == 0 for i in range(nb) for j in range(nb))
    line(tag("C05"), ok, f"(b) momentum Gram {nb}x{nb}, G_D verified {oku}; T[s]: A_ext/N_Y = {Am[1][1] / NY}, Z_E/N_Y = {Zm[1][1] / NY}, Z/vac = {Zm[1][1] / Vm[1][1]}")
    lam = [1 + Fr(1, (j + 1) * (2 * L + 2 * j + 3)) for j in range(3)]
    sub = lambda M, c, M2: [[M[i][j] - c * M2[i][j] for j in (1, 2, 3)] for i in (1, 2, 3)]
    ok = all(det3(sub(Am, lm, Vm)) == 0 for lm in lam) and det3(sub(Vm, 0, Vm)) != 0
    line(tag("C06"), ok, "(b) spectrum " + ", ".join(str(x) for x in lam))

    # static datum
    Bc = Fr(7 * L + 3, (2 * L + 3) * (L + 1))
    Ac = -(1 + 2 * Bc) / (4 * L + 2)
    gphi = grad(phi)
    xi = [padd(psc(pmul(SS, gphi[i]), Ac), psc(pmul(phi, XV[i]), Bc)) for i in I3]
    hs = [[psc(phi, -2) if i == j else {} for j in I3] for i in I3]
    qst = [Fr(0), Fr(-L * (L - 1), 2 * L + 3)]
    Tst = build_even(phi, qst, L)
    okt = tzero(tadd(tadd(hs, symgrad2(xi)), Tst, -1))
    kst = xtx(Tst)
    Gst = psc(padd(xdot(grad(kst)), psc(kst, 2)), Fr(1, 2))
    Gc = sph(pmul(Gst, phi)) / NY
    aext = A_TT(Tst, Tst) + sph(pmul(Gst, Gst)) / K_
    shape = -Gc / K_ + sph(pmul(xdot(xi), phi)) / NY
    val = Fr(2 * L * (L - 1), (L + 2) * (2 * L + 3))
    ok = okt and aext == val * NY and shape == Fr(2, L + 2)
    line(tag("C07"), ok, f"static datum: tensor identity {okt}; A_ext/N_Y = {aext / NY}; shape = {shape}; unextremised = {aext / NY - K_ * shape**2}")

    # ----- odd
    Ss = []
    for e in profs:
        S, v = build_odd(phi, e, L)
        Ss.append(S)
    ok = all(not ttr(S) and all(not c for c in tdiv(S)) and not xtx(S) and all(not psub(a_, pmul(rcart(e), b_)) for a_, b_ in zip(tx(S), v))
             for S, e in zip(Ss, profs))
    kk = L * (L + 1)
    psB = [psc(xtx(curl(S)), -1) for S in Ss]
    psBd = [xtx(curl(ricci(S))) for S in Ss]
    ok2 = all(not psub(pb, psc(pmul(rcart(e), phi), kk)) for pb, e in zip(psB, profs))
    ok3 = all(not psub(pb, psc(pmul(rcart(rtil(e, L)), phi), Fr(kk, 2))) for pb, e in zip(psBd, profs))
    line(tag("C08"), ok and ok2 and ok3, f"S[eta] TT, S.x = eta v: {ok}; psi_B = l(l+1) eta phi: {ok2}; psidot_B = (1/2) l(l+1) etat phi: {ok3}")
    gB = [grad(p_) for p_ in psB]
    AK = [[ballN(tdot(Ss[i], Ss[j])) for j in range(nb)] for i in range(nb)]
    EKg = [[CL * Fr(1, 2) * ballN(dotv(gB[i], gB[j])) for j in range(nb)] for i in range(nb)]
    EKm = [[CL * Fr(1, 2) * ball(pmul(psB[i], psB[j])) for j in range(nb)] for i in range(nb)]
    ok = all(AK[i][j] - EKg[i][j] - 2 * EKm[i][j] == 0 for i in range(nb) for j in range(nb))
    line(tag("C09"), ok, f"(a) odd K Gram {nb}x{nb} = c_l (E_1 + (1/2) int psi_B^2)")
    Ah = [[A_TT(Ss[i], Ss[j]) for j in range(nb)] for i in range(nb)]
    Eh = [[CL * Fr(1, 2) * ballN(pmul(psBd[i], psBd[j])) for j in range(nb)] for i in range(nb)]
    ok = all(Ah[i][j] - Eh[i][j] == 0 for i in range(nb) for j in range(nb)) and all(Ah[0][j] == 0 for j in range(nb))
    line(tag("C10"), ok, f"(a) odd h Gram {nb}x{nb}; eta = 1 row zero")

    # K gauge move
    c0, c1 = Fr(2, 3), Fr(-5, 7)
    X = pmul(rcart([c0, c1]), phi)
    gX = grad(X)
    HX = [[pdiff(gX[i], j) for j in I3] for i in I3]
    ok = True
    for i in (1, 2):
        Kp = tadd(Ts[i], HX)
        trK = ttr(Kp)
        lhs = ballN(psub(tdot(Kp, Kp), pmul(trK, trK)))
        xg = xdot(gX)
        rhs = ballN(tdot(Ts[i], Ts[i])) + 2 * sph(pmul(kap[i], X)) + sph(psub(dotv(gX, gX), pmul(xg, xg))) - 2 * ball(dotv(gX, gX))
        ok = ok and lhs == rhs
    # stationary point: X harmonic with tau = -kappa/(l(l-1))
    i = 2
    Xh = psc(phi, -r1(profs[i]) / (L * (L - 1)))
    gX = grad(Xh)
    HX = [[pdiff(gX[a_], b_) for b_ in I3] for a_ in I3]
    Kp = tadd(Ts[i], HX)
    ok = ok and ballN(tdot(Kp, Kp)) == Am[i][i] and not xtx(Kp) == {} and sph(pmul(xtx(Kp), xtx(Kp))) == 0
    line(tag("C11"), ok, "K gauge move formula; stationary representative has kappa = 0 on the edge and int N |K~|^2 = A_ext")
    Vk = [[EKg[i][j] + EKm[i][j] for j in range(nb)] for i in range(nb)]
    lamo = [1 + Fr(1, (2 * k + 1) * (L + k + 1)) for k in range(3)]
    sub0 = lambda M, c, M2: [[M[i][j] - c * M2[i][j] for j in (0, 1, 2)] for i in (0, 1, 2)]
    ok = all(det3(sub0(AK, lm, Vk)) == 0 for lm in lamo) and Vk[0][0] / AK[0][0] == Fr(L + 1, L + 2)
    line(tag("C12"), ok, "odd spectrum " + ", ".join(str(x) for x in lamo) + f"; stationary magnetic c_l E_1/A = {Vk[0][0] / AK[0][0]}")

    # plants
    t = Fr(L + 2, L + 1)
    bad = [[Aun[i][j] + GG[i][j] / K_ * (2 * t - t * t) - E1f[i][j] for j in range(nb)] for i in range(nb)]
    plant(tag("CP1"), any(v != 0 for row in bad for v in row) and aext + sph(pmul(Gst, Gst)) / K_ * (2 * t - t * t - 1) != val * NY,
          f"(c) shape 2/(l+1): static datum A/N_Y = {(aext + sph(pmul(Gst, Gst)) / K_ * (2 * t - t * t - 1)) / NY} against {val}")
    bad = [[ballN(tdot(Ts[i], Ts[j])) - sph(pmul(kap[i], kap[j])) / (L * (L + 1)) - Vm[i][j] - Zm[i][j] for j in range(nb)] for i in range(nb)]
    plant(tag("CP2"), any(v != 0 for row in bad for v in row), "(b) corner coefficient 1/(l(l+1))")
    plant(tag("CP3"), any(AK[i][j] - EKg[i][j] - EKm[i][j] != 0 for i in range(nb) for j in range(nb)), "(a) weight 1 on odd K")
    f_ = Fr(L + 3, L + 2)
    plant(tag("CP4"), any(Aun[i][j] + GG[i][j] / K_ - f_ * E1f[i][j] != 0 for i in range(nb) for j in range(nb)), "c_l (l+3)/(l+2) in (c)")

    # reproduction
    if L == 2:
        line("R1", aext / NY == Fr(1, 7) and aext / NY - K_ * shape**2 == Fr(-6, 7), f"static: A_ext/N_Y = {aext / NY}, unextremised {aext / NY - K_ * shape**2}")
        line("R2", Am[1][1] / NY == Fr(8, 27) and Zm[1][1] / NY == Fr(1, 27) and Am[1][1] / NY * Fr(4, 15) == Fr(32, 405) and Zm[1][1] / NY * Fr(4, 15) == Fr(4, 405)
             and Zm[1][1] / Vm[1][1] == Fr(1, 7), f"T[s] momentum: {Am[1][1] / NY} x 4/15 = {Am[1][1] / NY * Fr(4, 15)}; Z: {Zm[1][1] / NY} x 4/15 = {Zm[1][1] / NY * Fr(4, 15)}")
    if L in (2, 3, 4):
        exp = {2: "8/7, 19/18, 34/33", 3: "10/9, 23/22, 40/39", 4: "12/11, 27/26, 46/45"}[L]
        line(tag("R3"), ", ".join(str(x) for x in lam) == exp, f"spectrum {exp}")

    # ----- low l: Christoffel
    if L <= 4:
        ok = True
        for h in (Ts[1], Ts[2], tadd(Ts[1], Ts[2]), tadd(Ts[0], Ts[2]), Ss[0], Ss[1], Ss[2], tadd(Ss[1], Ss[2]), tadd(Ts[1], Ss[1])):
            ok = ok and Aun_chr(h) == A_TT(h, h)
        line(tag("X01"), ok, "Christoffel second-order area side = TT formula (9 data)")
        cs = (Fr(1, 3), Fr(-1, 2), Fr(2, 5))
        xi2 = [padd(psc(gphi[i], cs[0]), psc(pmul(phi, XV[i]), cs[1]), psc(pmul(SS, gphi[i]), cs[2])) for i in I3]
        sg = sph(pmul(xdot(xi2), phi)) / NY
        ok = sg != 0
        for i in (1, 2):
            Gi = sph(pmul(Gp[i], phi)) / NY
            ok = ok and Aun_chr(tadd(Ts[i], symgrad2(xi2))) == Aun[i][i] - 2 * sg * Gi * NY - K_ * sg * sg * NY
        line(tag("X02"), ok, f"even gauge move, sigma = {sg} Y")
        v3 = Aun_chr(hs) / NY
        line(tag("X03"), v3 == val - Fr(4 * (L - 1), L + 2), f"static datum unextremised by Christoffel: {v3}")
        v = cross(XV, gphi)
        aa = rcart([Fr(1), Fr(1, 2)])
        xo = [pmul(aa, v[i]) for i in I3]
        line(tag("X04"), Aun_chr(tadd(Ss[1], symgrad2(xo))) == Ah[1][1] and Ah[1][1] != 0, "odd gauge move leaves A_un unchanged")
    RES.setdefault("table", {})[str(L)] = {
        "c_l": str(CL), "N_Y": str(NY),
        "field_s_s_over_NY": str((Aun[1][1] + GG[1][1] / K_) / NY),
        "mom_s_s_over_NY": str(Am[1][1] / NY), "Z_s_s_over_NY": str(Zm[1][1] / NY),
        "static_over_NY": str(aext / NY), "shape": str(shape),
        "oddK_1_1_over_NY": str(AK[0][0] / NY), "oddh_s_s_over_NY": str(Ah[1][1] / NY),
        "spectrum_even": [str(x) for x in lam], "spectrum_odd": [str(x) for x in lamo], "seconds": round(time.time() - t0, 1)}
    print(f"   l={L} c_l={CL} N_Y={NY} done in {time.time() - t0:.1f}s", flush=True)


if not NOSYM:
    symbolic_part()
    print(f"   symbolic part done in {time.time() - T0:.1f}s", flush=True)
for L in range(LMIN, LMAX + 1):
    cart_part(L)
npass = sum(1 for x in RES["lines"] if x[1] == "PASS")
nfail = len(RES["fails"])
ncaught = sum(1 for x in RES["lines"] if x[1] == "CAUGHT")
print(f"SUMMARY pass {npass} fail {nfail} plants caught {ncaught} missed {len(RES['plants_missed'])} FAILED LINES: {RES['fails'] or 'none'} MISSED PLANTS: {RES['plants_missed'] or 'none'} time {time.time() - T0:.1f}s")
RES["summary"] = {"pass": npass, "fail": nfail, "caught": ncaught, "missed": len(RES["plants_missed"])}
with open(OUTF, "w", encoding="utf-8") as fh_:
    json.dump(RES, fh_, indent=1)
