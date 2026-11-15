"""Lane CHAIN-ALL-L (wave 73): the chain fields psi_E = x.E.x, psi_B = x.B.x of linearised gravity, every l.
Exact rational arithmetic, self-contained (no import from other lanes).
Usage: python chain_all_l.py LMIN LMAX
Pass lines: ../logs/00_passlines.log (fixed before the first run)."""
import sys, os, json, time
from math import factorial, comb
from sympy.polys.rings import ring
from sympy import QQ
import sympy as sp

R, x, y, z = ring("x,y,z", QQ)
X = (x, y, z)
I3 = range(3)
S2 = x**2 + y**2 + z**2
ZERO = R.zero
GAUSS = [False]          # True: every field is (polynomial) * exp(-r^2/2)
FAILS, MISSED = [], []
OUT = {}
T0 = time.time()


def log(*a):
    print(*a, flush=True)


def check(l, name, ok, info=""):
    OUT.setdefault(str(l), {})[name] = bool(ok)
    log(f"{'PASS' if ok else 'FAIL'} l={l} {name}" + (f"  {info}" if info != "" else ""))
    if not ok:
        FAILS.append((l, name))


def plant(l, name, caught, info=""):
    OUT.setdefault(str(l), {})[name] = bool(caught)
    log(f"{'CAUGHT' if caught else 'MISSED'} l={l} {name}" + (f"  {info}" if info != "" else ""))
    if not caught:
        MISSED.append((l, name))


def eps(i, j, k):
    return (i - j) * (j - k) * (k - i) // 2


# ---------------------------------------------------------------- calculus (Gaussian aware)
def d(P, i):
    out = P.diff(X[i])
    if GAUSS[0]:
        out = out - X[i] * P
    return out


def lap(P):
    return sum((d(d(P, i), i) for i in I3), ZERO)


def tlap(T):
    return [[lap(T[i][j]) for j in I3] for i in I3]


def trace(T):
    return T[0][0] + T[1][1] + T[2][2]


def div(T):
    return [sum((d(T[i][j], i) for i in I3), ZERO) for j in I3]


def curl(T):
    C = [[ZERO] * 3 for _ in I3]
    for i in I3:
        for k in I3:
            for m in I3:
                e = eps(i, k, m)
                if e:
                    for j in I3:
                        C[i][j] = C[i][j] + e * d(T[m][j], k)
    return C


def xTx(T):
    return sum((X[i] * X[j] * T[i][j] for i in I3 for j in I3), ZERO)


def tadd(A, B):
    return [[A[i][j] + B[i][j] for j in I3] for i in I3]


def tscale(A, c):
    return [[c * A[i][j] for j in I3] for i in I3]


def tzero(A):
    return all(not A[i][j] for i in I3 for j in I3)


def tsym(A):
    return all(not (A[i][j] - A[j][i]) for i in I3 for j in I3)


def isTT(A):
    return tsym(A) and not trace(A) and all(not c for c in div(A))


# ---------------------------------------------------------------- radial profiles: coefficient lists in s = r^2
def q_(a, b=1):
    return QQ(a, b)


def padd(a, b):
    n = max(len(a), len(b))
    a = list(a) + [q_(0)] * (n - len(a))
    b = list(b) + [q_(0)] * (n - len(b))
    return [u + v for u, v in zip(a, b)]


def pscale(a, c):
    return [c * u for u in a]


def pshift(a):
    return [q_(0)] + list(a)


def pdivs(a):
    assert a[0] == 0, "profile not divisible by s"
    return list(a[1:]) or [q_(0)]


def Dp(c):
    """D q = q'/r as an operator on s-coefficients (Gaussian mode: q = c(s) exp(-s/2))."""
    dd = [2 * k * ck for k, ck in enumerate(c)][1:] or [q_(0)]
    if GAUSS[0]:
        dd = padd(dd, pscale(c, q_(-1)))
    return dd


_S2P = [R.one]


def spoly(c):
    while len(_S2P) < len(c):
        _S2P.append(_S2P[-1] * S2)
    return sum((ck * _S2P[k] for k, ck in enumerate(c) if ck != 0), ZERO)


def unit(k):
    return [q_(0)] * k + [q_(1)]


def qtilde(c, l):
    D1 = Dp(c)
    return padd(pshift(Dp(D1)), pscale(D1, q_(2 * l + 3)))


# ---------------------------------------------------------------- sector objects
def harmonic(l):
    n = l - 1
    re = ZERO
    for j in range(0, n + 1, 2):
        re = re + comb(n, j) * (-1) ** (j // 2) * x ** (n - j) * y ** j
    phi = re * z
    assert not sum((phi.diff(v).diff(v) for v in X), ZERO)
    return phi


def ebasis(phi):
    g = [phi.diff(v) for v in X]
    H = [[g[i].diff(X[j]) for j in I3] for i in I3]
    Pp = [[X[i] * g[j] + X[j] * g[i] for j in I3] for i in I3]
    Dl = [[phi if i == j else ZERO for j in I3] for i in I3]
    XX = [[phi * X[i] * X[j] for j in I3] for i in I3]
    return H, Pp, Dl, XX


def comb_t(Fs, Bs):
    return [[sum((F * B[i][j] for F, B in zip(Fs, Bs)), ZERO) for j in I3] for i in I3]


def even_f(Q, l):
    DQ = Dp(Q)
    beta = pscale(DQ, q_(-1, l + 1))
    alpha = pscale(padd(Q, pscale(pshift(DQ), q_(1, l + 1))), q_(1, l))
    Da = Dp(alpha)
    f2 = pscale(padd(beta, pscale(Da, q_(-2))), q_(1, l + 2))
    f3 = padd(pscale(Da, q_(-1)), pscale(f2, q_(-(l + 1))))
    f4 = pdivs(pscale(padd(pscale(f2, q_(2 * l)), pscale(f3, q_(3))), q_(-1)))
    f1 = pscale(padd(alpha, pscale(pshift(f2), q_(-1))), q_(1, l - 1))
    return f1, f2, f3, f4


def build_even(phi, Q, l):
    return comb_t([spoly(f) for f in even_f(Q, l)], ebasis(phi))


def obasis(phi):
    g = [phi.diff(v) for v in X]
    v = [sum((eps(i, j, k) * X[j] * g[k] for j in I3 for k in I3), ZERO) for i in I3]
    A1 = [[X[i] * v[j] + X[j] * v[i] for j in I3] for i in I3]
    A2 = [[v[j].diff(X[i]) + v[i].diff(X[j]) for j in I3] for i in I3]
    return A1, A2


def build_odd(phi, eta, l, planted=False):
    De = Dp(eta)
    g1 = pscale(De, q_(-1, l + 1 if planted else l + 2))
    g2 = pscale(padd(eta, pscale(pshift(De), q_(1, l + 2))), q_(1, l - 1))
    return comb_t([spoly(g1), spoly(g2)], obasis(phi))


# ---------------------------------------------------------------- linear algebra for the dimension checks
def nullity(cols):
    """cols: list of dicts key -> rational; nullity of the matrix with these columns."""
    keys = sorted({k for c in cols for k in c}, key=str)
    M = sp.Matrix([[sp.Rational(int(c[k].numerator), int(c[k].denominator)) if k in c else 0 for c in cols] for k in keys])
    return len(cols) - M.rank()


def flat(tag, P):
    return {(tag, m): c for m, c in P.terms()}


# ---------------------------------------------------------------- Fourier transform of polynomial * exp(-r^2/2)
_HE = [[q_(1)], [q_(0), q_(1)]]


def he(n):
    while len(_HE) <= n:
        k = len(_HE) - 1
        a, b = _HE[-1], _HE[-2]
        new = [q_(0)] + list(a)
        for i, c in enumerate(b):
            new[i] -= k * c
        _HE.append(new)
    return _HE[n]


_HEP = {}


def hep(n, iv):
    if (n, iv) not in _HEP:
        _HEP[n, iv] = sum((c * X[iv] ** i for i, c in enumerate(he(n)) if c != 0), ZERO)
    return _HEP[n, iv]


def ft(P):
    """Returns (real polynomial Pt, n0): FT = (2 pi)^(3/2) (-i)^n0 Pt(k) exp(-k^2/2)."""
    if not P:
        return ZERO, None
    par = {sum(m) % 2 for m, _ in P.terms()}
    assert len(par) == 1, "mixed parity"
    n0 = par.pop()
    out = ZERO
    for (a, b, c), v in P.terms():
        sg = (-1) ** ((a + b + c - n0) // 2)
        out = out + (v * sg) * (hep(a, 0) * hep(b, 1) * hep(c, 2))
    return out, n0


def ftt(A):
    r = [[ft(A[i][j]) for j in I3] for i in I3]
    par = {r[i][j][1] for i in I3 for j in I3 if r[i][j][1] is not None}
    assert len(par) == 1
    return [[r[i][j][0] for j in I3] for i in I3], par.pop()


def dfact(n):
    r = 1
    while n > 1:
        r *= n
        n -= 2
    return r


def Iw(Qp, s):
    """int k^s Q(k) exp(-k^2) d^3k with the common constant dropped (2 pi for s = +-1; 2 pi sqrt(pi) for s = 0)."""
    tot = q_(0)
    for (a, b, c), v in Qp.terms():
        if a % 2 or b % 2 or c % 2:
            continue
        Dg = a + b + c
        ang = q_(dfact(a - 1) * dfact(b - 1) * dfact(c - 1), dfact(Dg + 1))
        m = Dg + 2 + s
        rad = q_(dfact(m - 1), 2 ** (m // 2)) if s == 0 else q_(factorial((m - 1) // 2))
        tot += v * ang * rad
    return tot


def tpair(A, B, s):
    (a, na), (b, nb) = A, B
    assert na == nb
    return Iw(sum((a[i][j] * b[i][j] for i in I3 for j in I3), ZERO), s)


def spair(A, B, s):
    (a, na), (b, nb) = A, B
    assert na == nb
    return Iw(a * b, s)


def gram(vs, pair, s, c=1):
    return [[c * pair(vs[a], vs[b], s) for b in range(len(vs))] for a in range(len(vs))]


def eqm(A, B, c):
    return all(A[a][b] == c * B[a][b] for a in range(len(A)) for b in range(len(A)))


def nzm(A):
    return all(v != 0 for r in A for v in r)


# ================================================================== L01: the wave-equation identity, generic M
def generic_identity():
    xs = sp.symbols('x1 x2 x3')
    M = [[sp.Function('M%d%d' % (min(i, j), max(i, j)))(*xs) for j in I3] for i in I3]
    L = lambda f: sum(sp.diff(f, v, 2) for v in xs)
    lhs = L(sum(xs[i] * xs[j] * M[i][j] for i in I3 for j in I3))
    rhs = (sum(xs[i] * xs[j] * L(M[i][j]) for i in I3 for j in I3)
           + 4 * sum(xs[i] * sp.diff(M[i][j], xs[j]) for i in I3 for j in I3) + 2 * sum(M[i][i] for i in I3))
    return sp.simplify(sp.expand(lhs - rhs)) == 0


# ================================================================== per-l run
def run(l):
    t0 = time.time()
    CL = q_(4 * factorial(l - 2), factorial(l + 2))
    CLBAD = q_(4 * factorial(l - 1), factorial(l + 2))
    phi = harmonic(l)
    LL = q_(l * (l + 1))
    Qs = [unit(k) for k in range(3)]

    okTT = okq = okOTT = okOe = okX = True
    for g in (False, True):
        GAUSS[0] = g
        Ts = [build_even(phi, Q, l) for Q in Qs]
        Ss = [build_odd(phi, Q, l) for Q in Qs]
        okTT &= all(isTT(T) and not tzero(T) for T in Ts)
        okq &= all(not (xTx(T) - spoly(Q) * phi) for T, Q in zip(Ts, Qs))
        okOTT &= all(isTT(S) and not tzero(S) for S in Ss)
        okOe &= all(not xTx(S) and not (xTx(curl(S)) + LL * spoly(Q) * phi) for S, Q in zip(Ss, Qs))
        for Q, T, S in zip(Qs, Ts, Ss):
            qt = qtilde(Q, l)
            okX &= tzero(tadd(curl(S), tscale(build_even(phi, pscale(Q, -LL), l), q_(-1))))
            okX &= tzero(tadd(curl(T), tscale(build_odd(phi, pscale(qt, 1 / LL), l), q_(-1))))
            okX &= tzero(tadd(tlap(T), tscale(build_even(phi, qt, l), q_(-1))))
            okX &= tzero(tadd(tlap(S), tscale(build_odd(phi, qt, l), q_(-1))))
        if g:
            bad = build_odd(phi, Qs[1], l, planted=True)
            plant(l, "P3 odd build with g1 = -D eta/(l+1): divergence nonzero", any(c for c in div(bad)))
    check(l, "L02 E-TT", okTT)
    check(l, "L03 E-q", okq)
    check(l, "L05 O-TT", okOTT)
    check(l, "L06 O-eta", okOe)
    check(l, "L08 X-curl", okX)

    # dimension checks, polynomial mode
    GAUSS[0] = False
    H, Pp, Dl, XX = ebasis(phi)
    sc = lambda B, k: [[_S2P[k] * B[i][j] for j in I3] for i in I3]
    spoly(unit(2))
    eb = [sc(H, 0), sc(H, 1), sc(H, 2), sc(Pp, 0), sc(Pp, 1), sc(Dl, 0), sc(Dl, 1), sc(XX, 0)]

    def cons(T, extra):
        c = flat('tr', trace(T))
        for j, w in enumerate(div(T)):
            c.update(flat('d%d' % j, w))
        if extra is not None:
            c.update(flat('x', extra(T)))
        return c
    n1 = nullity([cons(T, None) for T in eb])
    n2 = nullity([cons(T, xTx) for T in eb])
    check(l, "L04 E-dim", n1 == 3 and n2 == 0, f"TT dim {n1}, with x.T.x=0: {n2}")
    A1, A2 = obasis(phi)
    ob = [sc(A1, 0), sc(A1, 1), sc(A2, 0), sc(A2, 1), sc(A2, 2)]
    n1 = nullity([cons(T, None) for T in ob])
    n2 = nullity([cons(T, lambda S: xTx(curl(S))) for T in ob])
    check(l, "L07 O-dim", n1 == 3 and n2 == 0, f"TT dim {n1}, with x.curlS.x=0: {n2}")

    # Gaussian mode from here on
    GAUSS[0] = True
    Ts = [build_even(phi, Q, l) for Q in Qs]
    Ss = [build_odd(phi, Q, l) for Q in Qs]
    half = q_(1, 2)

    def evolve(h, K):
        hn = [h, tscale(K, q_(2))]
        hn.append(tlap(hn[0]))
        hn.append(tlap(hn[1]))
        E = [tscale(tlap(hn[n]), -half) for n in (0, 1)]
        B = [tscale(curl(hn[n + 1]), -half) for n in (0, 1)]
        return E, B

    def wave_ok(M):
        return isTT(M) and not (xTx(tlap(M)) - lap(xTx(M)))
    E, B = evolve(tadd(Ts[0], Ts[1]), tadd(Ts[1], Ts[2]))
    check(l, "L09 W-ev", all(wave_ok(M) for M in E + B) and all(not xTx(M) for M in B) and all(xTx(M) for M in E))
    E, B = evolve(tadd(Ss[0], Ss[1]), tadd(Ss[1], Ss[2]))
    check(l, "L10 W-od", all(wave_ok(M) for M in E + B) and all(not xTx(M) for M in E) and all(xTx(M) for M in B))

    # chain fields of the basis data
    psE = [xTx(tscale(tlap(T), -half)) for T in Ts]            # psi_E of h = T
    pdE = [xTx(tscale(tlap(T), q_(-1))) for T in Ts]           # psidot_E of K = T
    psB = [xTx(tscale(curl(S), q_(-1))) for S in Ss]           # psi_B of K = S
    pdB = [xTx(curl(tscale(tlap(S), -half))) for S in Ss]      # psidot_B of h = S
    k2 = S2
    a, _ = ft(lap(psE[1]))
    b, _ = ft(psE[1])
    check(l, "L11 FT", not (a + k2 * b) and bool(b))

    FT_T = [ftt(T) for T in Ts]
    FT_S = [ftt(S) for S in Ss]
    FpsE, FpdE, FpsB, FpdB = ([ft(p) for p in v] for v in (psE, pdE, psB, pdB))

    Gh = gram(FT_T, tpair, 1, half)
    GK = gram(FT_T, tpair, -1, 2)
    Gs = gram(FT_T, tpair, 0)
    Sp, Sm = gram(FpsE, spair, 1), gram(FpsE, spair, -1)
    Dp_, Dm = gram(FpdE, spair, 1), gram(FpdE, spair, -1)
    Ss0 = [[spair(FpdE[a_], FpsE[b_], 0) for b_ in range(3)] for a_ in range(3)]
    check(l, "L12 G-Eh", eqm(Gh, Sp, CL) and nzm(Gh), f"[0][0]: {Gh[0][0]} = {CL * Sp[0][0]}")
    check(l, "L13 G-EK", eqm(GK, Dm, CL) and nzm(GK), f"[0][0]: {GK[0][0]} = {CL * Dm[0][0]}")
    check(l, "L14 G-Es", eqm(Gs, Ss0, CL) and nzm(Gs))

    Oh = gram(FT_S, tpair, 1, half)
    OK = gram(FT_S, tpair, -1, 2)
    Os = gram(FT_S, tpair, 0)
    Bp, Bm = gram(FpsB, spair, 1), gram(FpsB, spair, -1)
    Bdp, Bdm = gram(FpdB, spair, 1), gram(FpdB, spair, -1)
    Bs0 = [[spair(FpsB[a_], FpdB[b_], 0) for b_ in range(3)] for a_ in range(3)]
    check(l, "L15 G-OK", eqm(OK, Bp, CL) and nzm(OK), f"[0][0]: {OK[0][0]} = {CL * Bp[0][0]}")
    check(l, "L16 G-Oh", eqm(Oh, Bdm, CL) and nzm(Oh), f"[0][0]: {Oh[0][0]} = {CL * Bdm[0][0]}")
    check(l, "L17 G-Os", eqm(Os, Bs0, -CL) and nzm(Os))

    plant(l, "P1 c_l with (l-1)! in place of (l-2)!", not eqm(Gh, Sp, CLBAD) and not eqm(OK, Bp, CLBAD),
          f"ratio seen/planted = {Gh[0][0] / (CLBAD * Sp[0][0])}")
    # P1b added after run 01: P1 is vacuous at l = 2 because 1! = 0! (see logs/02_addendum.log)
    CLB2 = CL * q_(l + 3, l + 2)
    plant(l, "P1b c_l scaled by (l+3)/(l+2)", not eqm(Gh, Sp, CLB2) and not eqm(OK, Bp, CLB2))
    plant(l, "P2 omega and 1/omega exchanged on the scalar side", not eqm(Gh, Sm, CL) and not eqm(OK, Bm, CL))
    OUT[str(l)]["seconds"] = round(time.time() - t0, 1)
    OUT[str(l)]["c_l"] = str(CL)
    log(f"-- l={l} done in {time.time() - t0:.1f} s, c_l = {CL}")


if __name__ == "__main__":
    lmin, lmax = int(sys.argv[1]), int(sys.argv[2])
    ok = generic_identity()
    check(0, "L01 W0 generic identity Lap(x.M.x) = x.LapM.x + 4 x.divM + 2 trM", ok)
    for l in range(lmin, lmax + 1):
        run(l)
    OUT["fails"] = [list(f) for f in FAILS]
    OUT["plants_missed"] = [list(f) for f in MISSED]
    OUT["seconds"] = round(time.time() - T0, 1)
    outp = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "out", f"chain_all_l_{lmin}_{lmax}.json")
    os.makedirs(os.path.dirname(outp), exist_ok=True)
    json.dump(OUT, open(outp, "w"), indent=1)
    log(f"DONE l={lmin}..{lmax}: fails={FAILS} plants_missed={MISSED} total {time.time() - T0:.1f} s")
