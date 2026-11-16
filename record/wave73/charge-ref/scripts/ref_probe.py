# CHARGE-REF probe: the sector identity at one l, own code, from the definitions.
# Usage: python ref_probe.py L MODE OUT.json      MODE = poly (weight 1, exact) | exp (weight exp(-s), 60 digits)
# Every field is P(x) * w(s), w = exp(-LAM s), P a rational polynomial; derivative D_i = d_i - 2 LAM x_i acts on P.
import sys, json, time
from fractions import Fraction as Fr
from math import gcd
import mpmath as mp

L = int(sys.argv[1]); MODE = sys.argv[2]; OUT = sys.argv[3]
assert MODE in ("poly", "exp")
LAM = 1 if MODE == "exp" else 0
KMAX = 5
mp.mp.dps = 60
T0 = time.time()
RES = {"l": L, "mode": MODE, "lines": []}
I3 = range(3)


def say(tag, ok, text):
    s = "PASS" if ok else "FAIL"
    print(f"{s} {tag}[l={L},{MODE}] {text}", flush=True)
    RES["lines"].append([tag, s, text])


# ---------------- polynomials: dict {(a,b,c): Fraction}
def padd(*ps):
    out = {}
    for p in ps:
        for m, c in p.items():
            v = out.get(m, 0) + c
            if v:
                out[m] = v
            else:
                out.pop(m, None)
    return out


def psc(p, c):
    return {m: v * c for m, v in p.items()} if c else {}


def pmul(p, q):
    out = {}
    for (a1, b1, c1), u in p.items():
        for (a2, b2, c2), v in q.items():
            m = (a1 + a2, b1 + b2, c1 + c2)
            out[m] = out.get(m, 0) + u * v
    return {m: v for m, v in out.items() if v}


X = [{(1, 0, 0): Fr(1)}, {(0, 1, 0): Fr(1)}, {(0, 0, 1): Fr(1)}]
ONE = {(0, 0, 0): Fr(1)}
S = padd(*[pmul(x, x) for x in X])


def d0(p, i):
    out = {}
    for m, c in p.items():
        if m[i]:
            mm = list(m); mm[i] -= 1
            out[tuple(mm)] = c * m[i]
    return out


def D(p, i):          # derivative of P w
    if not LAM:
        return d0(p, i)
    return padd(d0(p, i), psc(pmul(X[i], p), -2 * LAM))


def lapD(p):
    return padd(*[D(D(p, i), i) for i in I3])


def gradD(p):
    return [D(p, i) for i in I3]


def xdot(v):
    return padd(*[pmul(X[i], v[i]) for i in I3])


def tx(A):
    return [padd(*[pmul(A[i][j], X[j]) for j in I3]) for i in I3]


def xtx(A):
    return xdot(tx(A))


def ttr(A):
    return padd(*[A[i][i] for i in I3])


def tdiv(A):
    return [padd(*[D(A[i][j], i) for i in I3]) for j in I3]


def tlin(coefs, tens):
    return [[padd(*[psc(t[i][j], c) for c, t in zip(coefs, tens)]) for j in I3] for i in I3]


def ricci(h):
    tr = ttr(h)
    dv = tdiv(h)
    return [[psc(padd(D(dv[j], i), D(dv[i], j), psc(lapD(h[i][j]), -1), psc(D(D(tr, i), j), -1)), Fr(1, 2)) for j in I3] for i in I3]


EPS = {(0, 1, 2): 1, (1, 2, 0): 1, (2, 0, 1): 1, (0, 2, 1): -1, (2, 1, 0): -1, (1, 0, 2): -1}


def curl(T):          # eps_ikl D_k T_lj
    return [[padd(*[psc(D(T[l_][j], k), e) for (i2, k, l_), e in EPS.items() if i2 == i]) for j in I3] for i in I3]


# ---------------- integrals.  pair(p,q) -> {degree: sum c1 c2 <x^m1 x^m2>_sphere}, sphere average of monomials
DF = [1] * 400
for n_ in range(2, 400):
    DF[n_] = DF[n_ - 2] * (n_ - 1) if n_ % 2 == 0 else 0     # DF[n] = (n-1)!! for even n
ODD = [1] * 400
for n_ in range(1, 400):
    ODD[n_] = ODD[n_ - 2] * n_ if n_ >= 2 else 1               # ODD[n] = n!! (used for odd n)


def intform(p):
    den = 1
    for c in p.values():
        den = den * c.denominator // gcd(den, c.denominator)
    cls = {}
    for m, c in p.items():
        cls.setdefault((m[0] % 2, m[1] % 2, m[2] % 2), []).append((m[0], m[1], m[2], int(c * den)))
    return den, cls


def pair(p, q):
    if not p or not q:
        return {}
    dp, cp = intform(p)
    dq, cq = intform(q)
    acc = {}
    for par, lp in cp.items():
        lq = cq.get(par)
        if not lq:
            continue
        for a1, b1, c1, u in lp:
            for a2, b2, c2, v in lq:
                n = a1 + b1 + c1 + a2 + b2 + c2
                acc[n] = acc.get(n, 0) + u * v * DF[a1 + a2] * DF[b1 + b2] * DF[c1 + c2]
    return {n: Fr(v, dp * dq * ODD[n + 1]) for n, v in acc.items() if v}


def dsum(*ds):
    out = {}
    for d_ in ds:
        for n, v in d_.items():
            out[n] = out.get(n, 0) + v
    return out


def dsc(d_, c):
    return {n: v * c for n, v in d_.items()}


def tpair(A, B):
    return dsum(*[pair(A[i][j], B[i][j]) for i in I3 for j in I3])


def vpair(u, v):
    return dsum(*[pair(u[i], v[i]) for i in I3])


def dpair(A, B):
    return dsum(*[pair(D(A[i][j], k), D(B[i][j], k)) for i in I3 for j in I3 for k in I3])


_J = {}


def J(n):             # int_0^1 r^n w^2 dr
    if n not in _J:
        if LAM:
            a = mp.mpf(n + 1) / 2
            _J[n] = mp.gammainc(a, 0, 2 * LAM) / (2 * mp.mpf(2 * LAM) ** a)
        else:
            _J[n] = Fr(1, n + 1)
    return _J[n]


W1 = mp.exp(-2 * LAM) if LAM else Fr(1)     # w(1)^2


def num(v):
    return (mp.mpf(v.numerator) / v.denominator) if LAM else v


def ball(d_):
    return sum((num(v) * J(n + 2) for n, v in d_.items()), 0 * W1)


def ballN(d_):
    return sum((num(v) * (J(n + 2) - J(n + 4)) / 2 for n, v in d_.items()), 0 * W1)


def sph(d_):
    return sum((num(v) for v in d_.values()), 0 * W1) * W1


def sph_plain(d_):
    return sum(d_.values(), Fr(0))


def close(a, b, rel):
    if not LAM:
        return a == b and a != 0
    return a != 0 and abs(a - b) <= rel * abs(a)


def show(v):
    return str(v) if not LAM else mp.nstr(v, 25)


# ---------------- harmonic, different from the lane's
re, im = ONE, {}
for _ in range(L - 2):
    re, im = padd(pmul(re, X[0]), psc(pmul(im, X[1]), -1)), padd(pmul(re, X[1]), pmul(im, X[0]))
rho2 = padd(pmul(X[0], X[0]), pmul(X[1], X[1]))
phi = pmul(im, padd(pmul(X[2], X[2]), psc(rho2, Fr(-1, 2 * (L - 1)))))
assert phi and not padd(*[d0(d0(phi, i), i) for i in I3]), "phi not harmonic"
NY = sph_plain(pair(phi, phi))
CL = Fr(4, (L - 1) * L * (L + 1) * (L + 2))
gphi = [d0(phi, i) for i in I3]
HESS = [[d0(gphi[i], j) for j in I3] for i in I3]
PP = [[padd(pmul(X[i], gphi[j]), pmul(gphi[i], X[j])) for j in I3] for i in I3]
DEL = [[phi if i == j else {} for j in I3] for i in I3]
XX = [[pmul(phi, pmul(X[i], X[j])) for j in I3] for i in I3]
vv = [padd(*[psc(pmul(X[j], gphi[k]), e) for (i2, j, k), e in EPS.items() if i2 == i]) for i in I3]
XV = [[padd(pmul(X[i], vv[j]), pmul(vv[i], X[j])) for j in I3] for i in I3]
MV = [[padd(d0(vv[j], i), d0(vv[i], j)) for j in I3] for i in I3]
SP = [ONE]
for _ in range(KMAX):
    SP.append(pmul(SP[-1], S))


def tmulp(A, p):
    return [[pmul(A[i][j], p) for j in I3] for i in I3]


def nullspace(cols):
    n = len(cols)
    keys = set()
    for c in cols:
        keys.update(c)
    piv = {}                                  # pivot column -> row (list of Fr), reduced
    for k in keys:
        row = [c.get(k, Fr(0)) for c in cols]
        for pc, pr in piv.items():
            if row[pc]:
                f = row[pc]
                row = [a - f * b for a, b in zip(row, pr)]
        j = next((t for t in range(n) if row[t]), None)
        if j is None:
            continue
        row = [a / row[j] for a in row]
        for pc in list(piv):
            if piv[pc][j]:
                f = piv[pc][j]
                piv[pc] = [a - f * b for a, b in zip(piv[pc], row)]
        piv[j] = row
    free = [t for t in range(n) if t not in piv]
    out = []
    for f in free:
        vec = [Fr(0)] * n
        vec[f] = Fr(1)
        for pc, pr in piv.items():
            vec[pc] = -pr[f]
        out.append(vec)
    return out


def constraint(E):
    c = {}
    for m, v in ttr(E).items():
        c[(3, m)] = v
    for j, p in enumerate(tdiv(E)):
        for m, v in p.items():
            c[(j, m)] = v
    return c


def solve_TT(bases):
    els = [tmulp(B, SP[k]) for B in bases for k in range(KMAX + 1)]
    ns = nullspace([constraint(E) for E in els])
    return els, ns


WA = [Fr(1), Fr(-2), Fr(3), Fr(1, 2), Fr(-1), Fr(2, 3), Fr(1, 5), Fr(-4, 3)]
WB = [Fr(1, 3), Fr(1), Fr(-1), Fr(2), Fr(5, 7), Fr(-3), Fr(1), Fr(1, 2)]


def combo(els, ns, w):
    coef = [sum((w[i] * ns[i][c] for i in range(len(ns))), Fr(0)) for c in range(len(els))]
    return tlin(coef, els)


NDIR = (Fr(1, 3), Fr(2, 3), Fr(2, 3))
PHIN = sum(c * NDIR[0] ** m[0] * NDIR[1] ** m[1] * NDIR[2] ** m[2] for m, c in phi.items())
assert PHIN != 0


def radial(p):        # p = F(s) phi : read F on the ray, verify the identity exactly
    F = {}
    for m, c in p.items():
        k2 = sum(m) - L
        assert k2 >= 0 and k2 % 2 == 0, "not of the form F(s) phi"
        F[k2 // 2] = F.get(k2 // 2, Fr(0)) + c * NDIR[0] ** m[0] * NDIR[1] ** m[1] * NDIR[2] ** m[2] / PHIN
    Fl = [F.get(k, Fr(0)) for k in range(max(F) + 1)] if F else [Fr(0)]
    chk = {}
    for k, c in enumerate(Fl):
        sp = ONE
        for _ in range(k):
            sp = pmul(sp, S)
        chk = padd(chk, psc(pmul(sp, phi), c))
    ok = not padd(p, psc(chk, -1))
    return Fl, ok


def green(Fa, Fb, dirichlet=True):
    # <f_a Y, G g_b Y>/N_Y with f = F(r^2) r^l w(r), kernel r_<^l (r_>^(-l-1) - r_>^l)/(2l+1)
    if not LAM:
        tot = Fr(0)
        for k, a in enumerate(Fa):
            for j, b in enumerate(Fb):
                tot += a * b * (Fr(1, (2 * j + 2 * L + 3) * (2 * k + 2 * j + 2 * L + 5)) + Fr(1, (2 * k + 2 * L + 3) * (2 * k + 2 * j + 2 * L + 5)))
        if dirichlet:
            tot -= sum(a / (2 * k + 2 * L + 3) for k, a in enumerate(Fa)) * sum(b / (2 * j + 2 * L + 3) for j, b in enumerate(Fb))
        return tot / (2 * L + 1)
    fa = [mp.mpf(c.numerator) / c.denominator for c in Fa]
    fb = [mp.mpf(c.numerator) / c.denominator for c in Fb]

    def f(F, r):
        return sum(c * r ** (2 * k) for k, c in enumerate(F)) * r ** L * mp.exp(-LAM * r * r)

    def I1(F, r):
        return sum(c * mp.gammainc(mp.mpf(2 * k + 2 * L + 3) / 2, 0, LAM * r * r) / (2 * mp.mpf(LAM) ** (mp.mpf(2 * k + 2 * L + 3) / 2)) for k, c in enumerate(F))

    def integrand(r):
        if r == 0:
            return mp.mpf(0)
        return (f(fa, r) * I1(fb, r) + f(fb, r) * I1(fa, r)) * r ** (1 - L)

    old = mp.mp.dps
    mp.mp.dps = 45
    tot = mp.quad(integrand, [0, mp.mpf(1) / 4, mp.mpf(1) / 2, mp.mpf(3) / 4, 1])
    if dirichlet:
        tot -= I1(fa, mp.mpf(1)) * I1(fb, mp.mpf(1))
    tot = tot / (2 * L + 1)
    mp.mp.dps = old
    return +tot


REL = mp.mpf(10) ** -30
RELG = mp.mpf(10) ** -20
PAIRS = (("AA", 0, 0), ("AB", 0, 1), ("BB", 1, 1))

# =================== ODD
els, ns = solve_TT([XV, MV])
print(f"odd null space dimension {len(ns)} of {len(els)} unknowns ({time.time() - T0:.1f}s)", flush=True)
Sd = [combo(els, ns, WA), combo(els, ns, WB)]
ok = all(not ttr(Sx) and all(not c for c in tdiv(Sx)) and not xtx(Sx) and any(Sx[i][j] for i in I3 for j in I3) for Sx in Sd)
# profile: S.x = eta v
say("P00o", ok, f"odd A, B trace-free, divergence-free, x.S.x = 0; null space dim {len(ns)}")
psB = [psc(xtx(curl(Sx)), -1) for Sx in Sd]
psBd = [xtx(curl(ricci(Sx))) for Sx in Sd]
etas = [radial(p) for p in psB]
RES["odd_profile_psiB_over_phi"] = [[str(c) for c in e[0]] for e in etas]
print("   psi_B/(w phi) for A:", [str(c) for c in etas[0][0]], "ok", etas[0][1], flush=True)
gB = [gradD(p) for p in psB]
Sxs = [tx(Sx) for Sx in Sd]
tab = {}
for nm, a, b in PAIRS:
    lhs = ballN(tpair(Sd[a], Sd[b]))
    rhs = CL * (ballN(vpair(gB[a], gB[b])) / 2 + ball(pair(psB[a], psB[b])))
    Zp = CL * ball(pair(psB[a], psB[b])) / 2
    say("P01", close(lhs, rhs, REL), f"odd K {nm}: lhs/N_Y = {show(lhs / NY)}  lhs - rhs = {show(lhs - rhs)}  Z/N_Y = {show(Zp / NY)}")
    tab["oddK_" + nm] = [show(lhs / NY), show(rhs / NY), show(Zp / NY)]
    lhs = ballN(dpair(Sd[a], Sd[b])) / 4 - sph(vpair(Sxs[a], Sxs[b])) / 2 + sph(pair(xdot(Sxs[a]), xdot(Sxs[b]))) / 4
    rhs = CL * ballN(pair(psBd[a], psBd[b])) / 2
    say("P02", close(lhs, rhs, REL), f"odd h {nm}: lhs/N_Y = {show(lhs / NY)}  lhs - rhs = {show(lhs - rhs)}")
    tab["oddh_" + nm] = [show(lhs / NY), show(rhs / NY)]
print(f"   odd done ({time.time() - T0:.1f}s)", flush=True)

# =================== EVEN
els, ns = solve_TT([HESS, PP, DEL, XX])
print(f"even null space dimension {len(ns)} of {len(els)} unknowns ({time.time() - T0:.1f}s)", flush=True)
Td = [combo(els, ns, WA), combo(els, ns, WB)]
ok = all(not ttr(Tx) and all(not c for c in tdiv(Tx)) and any(Tx[i][j] for i in I3 for j in I3) for Tx in Td)
say("P00e", ok, f"even A, B trace-free, divergence-free; null space dim {len(ns)}")
kap = [xtx(Tx) for Tx in Td]
qs = [radial(k) for k in kap]
RES["even_profile_q"] = [[str(c) for c in q[0]] for q in qs]
print("   q for A (x.T.x = q w phi):", [str(c) for c in qs[0][0]], "ok", qs[0][1], flush=True)
print("   q for B:", [str(c) for c in qs[1][0]], "ok", qs[1][1], flush=True)
psE = [xtx(ricci(Tx)) for Tx in Td]
psEd = [psc(p, 2) for p in psE]
Fd = [radial(p) for p in psEd]
okF = all(f[1] for f in Fd) and all(q[1] for q in qs)
Txs = [tx(Tx) for Tx in Td]
Gs = [psc(padd(xdot(gradD(k)), psc(k, 2)), Fr(1, 2)) for k in kap]
gE = [gradD(p) for p in psE]
planted = []
for nm, a, b in PAIRS:
    lhs = ballN(tpair(Td[a], Td[b])) - sph(pair(kap[a], kap[b])) / (L * (L - 1))
    vac = CL * ballN(pair(psEd[a], psEd[b])) / 2
    gr = green(Fd[a][0], Fd[b][0]) * num(NY)
    Z = CL * gr / 2
    rhs = vac + Z
    say("P03", okF and close(lhs, rhs, RELG if LAM else 0), f"even K {nm}: lhs/N_Y = {show(lhs / NY)}  lhs - rhs = {show(lhs - rhs)}  Z/N_Y = {show(Z / NY)}  radial form verified {okF}")
    tab["evenK_" + nm] = [show(lhs / NY), show(rhs / NY), show(Z / NY)]
    if a == b == 0:
        ratio = Z / vac
        bound = Fr(1, 2 * L + 3)
        say("P06", 0 < ratio < (num(bound) if LAM else bound), f"Z/vac on A = {show(ratio)}  (0 < . < 1/(2l+3) = {bound})")
    grf = green(Fd[a][0], Fd[b][0], dirichlet=False) * num(NY)
    dif = lhs - vac - CL * grf / 2
    planted.append((dif != 0) if not LAM else (abs(dif) > mp.mpf(10) ** -6 * abs(lhs)))
    print(f"   plant {nm}: free-space Green function gives lhs - rhs = {show(dif)}", flush=True)
    lhs = (ballN(dpair(Td[a], Td[b])) / 4 - sph(vpair(Txs[a], Txs[b])) / 2 + sph(pair(kap[a], kap[b])) / 4
           + sph(pair(Gs[a], Gs[b])) / ((L - 1) * (L + 2)))
    rhs = CL * (ballN(vpair(gE[a], gE[b])) / 2 + ball(pair(psE[a], psE[b])) / 2)
    say("P04", close(lhs, rhs, REL), f"even h {nm}: lhs/N_Y = {show(lhs / NY)}  lhs - rhs = {show(lhs - rhs)}")
    tab["evenh_" + nm] = [show(lhs / NY), show(rhs / NY)]
s_ = "CAUGHT" if all(planted) else "MISSED"
print(f"PLANT {s_} PL1[l={L},{MODE}] free-space Green function in P03", flush=True)
RES["lines"].append(["PL1", s_, "free-space Green function in P03"])

if not LAM:
    mono = lambda a, b, c, v=1: {(a, b, c): Fr(v)}
    Xp = padd(mono(2, 1, 0), mono(0, 0, 3), mono(1, 1, 1, -2), pmul(S, X[0]), mono(0, 4, 0), psc(pmul(padd(ONE, S), phi), Fr(1, 3)))
    gX = [d0(Xp, i) for i in I3]
    HX = [[d0(gX[i], j) for j in I3] for i in I3]
    Kp = [[padd(Td[0][i][j], HX[i][j]) for j in I3] for i in I3]
    trK = ttr(Kp)
    xg = xdot(gX)
    lhs = ballN(dsum(tpair(Kp, Kp), dsc(pair(trK, trK), -1)))
    cross = sph(pair(kap[0], Xp))
    rhs = ballN(tpair(Td[0], Td[0])) + 2 * cross + sph(dsum(vpair(gX, gX), dsc(pair(xg, xg), -1))) - 2 * ball(vpair(gX, gX))
    say("P05", lhs == rhs and cross != 0 and bool(padd(*[d0(gX[i], i) for i in I3])), f"corner identity, generic non-harmonic X: lhs - rhs = {lhs - rhs}; 2 oint kappa X / N_Y = {2 * cross / NY}")

RES["table"] = tab
RES["N_Y"] = str(NY)
RES["c_l"] = str(CL)
npass = sum(1 for x in RES["lines"] if x[1] == "PASS")
nfail = sum(1 for x in RES["lines"] if x[1] == "FAIL")
print(f"SUMMARY l={L} mode={MODE} pass {npass} fail {nfail} plant {s_} time {time.time() - T0:.1f}s", flush=True)
RES["summary"] = {"pass": npass, "fail": nfail, "plant": s_, "seconds": round(time.time() - T0, 1)}
with open(OUT, "w", encoding="utf-8") as fh:
    json.dump(RES, fh, indent=1)
