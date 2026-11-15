"""Exact sparse polynomial algebra in (t, x, y, z) over the rationals, plus
tensor helpers and exact integrals over the unit ball / unit sphere / R^3 with
a Gaussian weight. No external dependencies. Lane CHARGE-ID, wave 73."""
from fractions import Fraction as Fr
from math import comb

NV = 4  # variables: t, x, y, z  (index 0..3)
ZERO = {}


def const(c):
    c = Fr(c)
    return {(0, 0, 0, 0): c} if c != 0 else {}


def var(i):
    e = [0, 0, 0, 0]
    e[i] = 1
    return {tuple(e): Fr(1)}


T_, X_, Y_, Z_ = var(0), var(1), var(2), var(3)
XS = [None, X_, Y_, Z_]


def add(p, q):
    r = dict(p)
    for k, v in q.items():
        s = r.get(k, 0) + v
        if s == 0:
            r.pop(k, None)
        else:
            r[k] = s
    return r


def scale(p, c):
    c = Fr(c)
    if c == 0:
        return {}
    return {k: v * c for k, v in p.items()}


def sub(p, q):
    return add(p, scale(q, -1))


def mul(p, q):
    if not p or not q:
        return {}
    if len(p) > len(q):
        p, q = q, p
    r = {}
    for k1, v1 in p.items():
        a0, a1, a2, a3 = k1
        for k2, v2 in q.items():
            k = (a0 + k2[0], a1 + k2[1], a2 + k2[2], a3 + k2[3])
            r[k] = r.get(k, 0) + v1 * v2
    return {k: v for k, v in r.items() if v != 0}


def powp(p, n):
    r = const(1)
    for _ in range(n):
        r = mul(r, p)
    return r


def pdiff(p, i):
    r = {}
    for k, v in p.items():
        if k[i]:
            e = list(k)
            c = e[i]
            e[i] -= 1
            r[tuple(e)] = v * c
    return r


# Gaussian mode: fields are (polynomial) * exp(-r^2); spatial derivative of the
# polynomial prefactor is d_i p - 2 x_i p.  t-derivatives are unaffected.
GAUSS = [False]


def diff(p, i):
    r = pdiff(p, i)
    if GAUSS[0] and i > 0:
        r = sub(r, scale(mul(XS[i], p), 2))
    return r


def sumall(ps):
    r = {}
    for p in ps:
        r = add(r, p)
    return r


R2 = sumall([mul(XS[i], XS[i]) for i in (1, 2, 3)])
NLAPSE = scale(sub(const(1), R2), Fr(1, 2))


def at_t0(p):
    return {k: v for k, v in p.items() if k[0] == 0}


def tcoef(p, n):
    return {(0,) + k[1:]: v for k, v in p.items() if k[0] == n}


def dfact(n):  # double factorial, (-1)!! = 1
    r = 1
    while n > 1:
        r *= n
        n -= 2
    return r


def sph(p):
    """oint over unit sphere, in units of 4*pi. p must be t-free."""
    s = Fr(0)
    for (t, a, b, c), v in p.items():
        assert t == 0
        if a % 2 or b % 2 or c % 2:
            continue
        s += v * Fr(dfact(a - 1) * dfact(b - 1) * dfact(c - 1), dfact(a + b + c + 1))
    return s


def ball(p):
    """int over unit ball, in units of 4*pi."""
    s = Fr(0)
    for (t, a, b, c), v in p.items():
        assert t == 0
        if a % 2 or b % 2 or c % 2:
            continue
        s += v * Fr(dfact(a - 1) * dfact(b - 1) * dfact(c - 1),
                    dfact(a + b + c + 1) * (a + b + c + 3))
    return s


def gauss2(p):
    """int over R^3 of p * exp(-2 r^2), in units of (pi/2)^(3/2)."""
    s = Fr(0)
    for (t, a, b, c), v in p.items():
        assert t == 0
        if a % 2 or b % 2 or c % 2:
            continue
        s += v * Fr(dfact(a - 1) * dfact(b - 1) * dfact(c - 1), 2 ** (a + b + c))
    return s


# ---------- tensors: 3x3 lists indexed 0..2 for x,y,z ----------
I3 = (0, 1, 2)


def eps(i, j, k):
    return ((i - j) * (j - k) * (k - i)) // 2


def tzero():
    return [[{} for _ in I3] for _ in I3]


def tadd(A, B):
    return [[add(A[i][j], B[i][j]) for j in I3] for i in I3]


def tscale(A, c):
    return [[scale(A[i][j], c) for j in I3] for i in I3]


def tmulp(A, p):
    return [[mul(A[i][j], p) for j in I3] for i in I3]


def dot(A, B):
    return sumall([mul(A[i][j], B[i][j]) for i in I3 for j in I3])


def trace(A):
    return sumall([A[i][i] for i in I3])


def div(A):
    return [sumall([diff(A[i][j], i + 1) for i in I3]) for j in I3]


def lap(p):
    return sumall([diff(diff(p, i + 1), i + 1) for i in I3])


def tlap(A):
    return [[lap(A[i][j]) for j in I3] for i in I3]


def symgrad(w):
    return [[scale(add(diff(w[j], i + 1), diff(w[i], j + 1)), Fr(1, 2)) for j in I3] for i in I3]


def curl(A):
    """(curl T)_ij = eps_ikl d_k T_lj"""
    C = tzero()
    for i in I3:
        for j in I3:
            acc = {}
            for k in I3:
                for l in I3:
                    e = eps(i, k, l)
                    if e:
                        acc = add(acc, scale(diff(A[l][j], k + 1), e))
            C[i][j] = acc
    return C


def xTx(A):
    return sumall([mul(mul(XS[i + 1], XS[j + 1]), A[i][j]) for i in I3 for j in I3])


def Tx(A):
    return [sumall([mul(A[i][j], XS[j + 1]) for j in I3]) for i in I3]


def ric_lin(h):
    """Ric_ij = 1/2(-Lap h_ij + d_i d_k h_kj + d_j d_k h_ki - d_i d_j tr h)"""
    dv = div(h)
    tr = trace(h)
    R = tzero()
    for i in I3:
        for j in I3:
            R[i][j] = scale(sumall([scale(lap(h[i][j]), -1), diff(dv[j], i + 1),
                                    diff(dv[i], j + 1),
                                    scale(diff(diff(tr, i + 1), j + 1), -1)]), Fr(1, 2))
    return R


def rscal_lin(q):
    dv = div(q)
    return sub(sumall([diff(dv[j], j + 1) for j in I3]), lap(trace(q)))


def rscal_2(h):
    """O(eps^2) coefficient of the scalar curvature of delta + eps h. No gauge assumed."""
    dh = [[[diff(h[i][j], k + 1) for k in I3] for j in I3] for i in I3]  # dh[i][j][k]=d_k h_ij
    G1 = [[[scale(sumall([dh[j][k][i], dh[i][k][j], scale(dh[i][j][k], -1)]), Fr(1, 2))
            for j in I3] for i in I3] for k in I3]  # G1[k][i][j] = Gamma1^k_ij
    G2 = [[[scale(sumall([mul(h[k][l], G1[l][i][j]) for l in I3]), -1)
            for j in I3] for i in I3] for k in I3]
    # R1_ij
    R1 = ric_lin(h)
    # trace of R2_ij = d_k G2^k_ii - d_i G2^k_ik + G1^k_kl G1^l_ii - G1^k_il G1^l_ik
    G1c = [sumall([G1[k][k][l] for k in I3]) for l in I3]  # G1^k_kl
    acc = {}
    for i in I3:
        for k in I3:
            acc = add(acc, diff(G2[k][i][i], k + 1))
            acc = sub(acc, diff(G2[k][i][k], i + 1))
        for l in I3:
            acc = add(acc, mul(G1c[l], G1[l][i][i]))
            for k in I3:
                acc = sub(acc, mul(G1[k][i][l], G1[l][i][k]))
    return sub(acc, dot(h, R1))


def harmonic(l, c, s, perm=(1, 2, 3), part="re"):
    """Re or Im of (x_a + i(c x_b + s x_c))^l with c^2+s^2=1 (rational): harmonic, degree l."""
    assert Fr(c) ** 2 + Fr(s) ** 2 == 1
    u = XS[perm[0]]
    w = add(scale(XS[perm[1]], c), scale(XS[perm[2]], s))
    out = {}
    for k in range(l + 1):
        if part == "re" and k % 2 == 0:
            sg = (-1) ** (k // 2)
        elif part == "im" and k % 2 == 1:
            sg = (-1) ** ((k - 1) // 2)
        else:
            continue
        out = add(out, scale(mul(powp(u, l - k), powp(w, k)), sg * comb(l, k)))
    assert not lap(out) or GAUSS[0]
    return out
