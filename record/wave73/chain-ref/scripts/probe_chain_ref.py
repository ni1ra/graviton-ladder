"""CHAIN-REF probe (wave 73). Independent exact test of the ratio in step (iii) of chain-all-l.
Written from the definitions; imports nothing from chain_all_l.py and uses none of its radial formulas.
TT tensors are made with the local operator  Pi = Lap^2 P_TT  applied to a seed tensor:
  Pi(M)_ij = A_ik A_jl M_kl - (1/2) A_ij A_kl M_kl,   A_ik = Lap delta_ik - d_i d_k.
All fields are (polynomial) * exp(-r^2/2).  FT convention f^(k) = int f e^{-ik.x}: x_j -> i d/dk_j,
done by the recursion g_{a+1} = g_a' - k g_a (so FT of x^a e^{-x^2/2} is i^a g_a(k) e^{-k^2/2}, constants dropped).
Graviton forms: (1/2) int w |h^|^2, 2 int |K^|^2 / w, int K.h.   Scalar: int w |psi^|^2, int |psidot^|^2 / w, int psidot psi.
"""
import sys
from math import factorial
from sympy.polys.rings import ring
from sympy import QQ, Rational, gamma, sqrt, pi, simplify

R, x, y, z = ring("x,y,z", QQ)
X = (x, y, z)
I3 = range(3)
S = x**2 + y**2 + z**2
Z = R.zero
NFAIL = [0]


def say(ok, name, info=""):
    print(("PASS " if ok else "FAIL ") + name + ("  " + str(info) if info != "" else ""), flush=True)
    if not ok:
        NFAIL[0] += 1


def eps(i, j, k):
    return (i - j) * (j - k) * (k - i) // 2


def d(P, i):                       # derivative of P * exp(-r^2/2), Gaussian stripped
    return P.diff(X[i]) - X[i] * P


def lap(P):
    return sum((d(d(P, i), i) for i in I3), Z)


def A(M):                          # (A M)_ij = A_ik M_kj  (acts on the first index)
    return [[lap(M[i][j]) - sum((d(d(M[k][j], k), i) for k in I3), Z) for j in I3] for i in I3]


def transpose(M):
    return [[M[j][i] for j in I3] for i in I3]


def Pi(M):
    AM = A(M)                                  # A_ik M_kj
    AMA = transpose(A(transpose(AM)))          # A_ik A_jl M_kl
    tr = sum((AM[k][k] for k in I3), Z)        # A_kl M_kl (M symmetric)
    Atr = [[(lap(tr) if i == j else Z) - d(d(tr, i), j) for j in I3] for i in I3]
    return [[AMA[i][j] - QQ(1, 2) * Atr[i][j] for j in I3] for i in I3]


def tlap(T):
    return [[lap(T[i][j]) for j in I3] for i in I3]


def curl(T):                       # curl(T)_ij = eps_ikl d_k T_lj
    C = [[Z] * 3 for _ in I3]
    for i in I3:
        for k in I3:
            for l in I3:
                e = eps(i, k, l)
                if e:
                    for j in I3:
                        C[i][j] = C[i][j] + e * d(T[l][j], k)
    return C


def sc(T, c):
    return [[c * T[i][j] for j in I3] for i in I3]


def xTx(T):
    return sum((X[i] * X[j] * T[i][j] for i in I3 for j in I3), Z)


def isTT(T):
    sym = all(not (T[i][j] - T[j][i]) for i in I3 for j in I3)
    tr = not (T[0][0] + T[1][1] + T[2][2])
    dv = all(not sum((d(T[i][j], i) for i in I3), Z) for j in I3)
    nz = any(T[i][j] for i in I3 for j in I3)
    return sym and tr and dv and nz


# ---------------- Fourier transform: returns (Re, Im) polynomials in k (same ring, variables reused as k)
_G = [[R.one] * 1 for _ in I3]


def g(a, iv):
    L = _G[iv]
    while len(L) <= a:
        L.append(L[-1].diff(X[iv]) - X[iv] * L[-1])
    return L[a]


def ft(P):
    re, im = Z, Z
    for (a, b, c), v in P.terms():
        term = v * g(a, 0) * g(b, 1) * g(c, 2)
        n = (a + b + c) % 4                     # i^n
        if n == 0:
            re = re + term
        elif n == 1:
            im = im + term
        elif n == 2:
            re = re - term
        else:
            im = im - term
    return re, im


def herm(F, G):                    # Re of F conj(G)
    return F[0] * G[0] + F[1] * G[1]


def dfact(n):
    r = 1
    while n > 1:
        r *= n
        n -= 2
    return r


def I(Q, s):
    """int k^s Q(k) exp(-k^2) d^3k, exact (sympy number)."""
    tot = 0
    for (a, b, c), v in Q.terms():
        if a % 2 or b % 2 or c % 2:
            continue
        D = a + b + c
        ang = Rational(dfact(a - 1) * dfact(b - 1) * dfact(c - 1), dfact(D + 1))
        tot += Rational(int(v.numerator), int(v.denominator)) * ang * 4 * pi * gamma(Rational(D + 3 + s, 2)) / 2
    return tot


def tpair(TA, TB, s):
    return I(sum((herm(TA[i][j], TB[i][j]) for i in I3 for j in I3), Z), s)


def spair(a, b, s):
    return I(herm(a, b), s)


def ftt(T):
    return [[ft(T[i][j]) for j in I3] for i in I3]


def ratio(a, b):
    if b == 0:
        return None
    return simplify(a / b)


def run(name, phi, l, profs):
    assert not sum((phi.diff(v).diff(v) for v in X), Z), "not harmonic"
    c = Rational(4 * factorial(l - 2), factorial(l + 2))
    gph = [phi.diff(v) for v in X]
    v = [sum((eps(i, j, k) * X[j] * gph[k] for j in I3 for k in I3), Z) for i in I3]
    Ts = [Pi([[X[i] * X[j] * phi * p for j in I3] for i in I3]) for p in profs]
    Ss = [Pi([[(X[i] * v[j] + X[j] * v[i]) * p for j in I3] for i in I3]) for p in profs]
    say(all(isTT(T) for T in Ts + Ss), f"{name} seeds give nonzero TT tensors")
    n = len(profs)
    # chain fields
    psE = [xTx(sc(tlap(T), QQ(-1, 2))) for T in Ts]       # psi_E of h = T
    pdE = [xTx(sc(tlap(T), QQ(-1))) for T in Ts]          # psidot_E of K = T
    psB = [xTx(sc(curl(S_), QQ(-1))) for S_ in Ss]        # psi_B of K = S
    pdB = [xTx(curl(sc(tlap(S_), QQ(-1, 2)))) for S_ in Ss]  # psidot_B of h = S
    say(all(bool(p) for p in psE + pdE + psB + pdB), f"{name} chain fields nonzero")
    say(all(not xTx(curl(T)) for T in Ts) and all(not xTx(tlap(S_)) for S_ in Ss),
        f"{name} psi_B of even data = 0 and psi_E of odd data = 0")
    # profile of psi_E: degree in r^2 of q = psi_E/phi (to show it is outside {1,s,s^2})
    print(f"   {name}: total degrees of psi_E polynomials: {[max(sum(m) for m, _ in p.terms()) for p in psE]} (phi has degree {l})")
    FT, FS = [ftt(T) for T in Ts], [ftt(S_) for S_ in Ss]
    FpsE, FpdE, FpsB, FpdB = ([ft(p) for p in L] for L in (psE, pdE, psB, pdB))
    seen = set()
    ok = True
    for a in range(n):
        for b in range(n):
            rr = [ratio(tpair(FT[a], FT[b], 1) / 2, spair(FpsE[a], FpsE[b], 1)),
                  ratio(2 * tpair(FT[a], FT[b], -1), spair(FpdE[a], FpdE[b], -1)),
                  ratio(tpair(FT[a], FT[b], 0), spair(FpdE[a], FpsE[b], 0)),     # int K_a.h_b / int psidot[K_a] psi[h_b]
                  ratio(2 * tpair(FS[a], FS[b], -1), spair(FpsB[a], FpsB[b], 1)),
                  ratio(tpair(FS[a], FS[b], 1) / 2, spair(FpdB[a], FpdB[b], -1))]
            ro = ratio(tpair(FS[a], FS[b], 0), spair(FpsB[a], FpdB[b], 0))       # int K_a.h_b / int psi_B[K_a] psidot_B[h_b]
            seen.update(str(r) for r in rr)
            ok &= all(r == c for r in rr) and ro == -c
            if a == 0 and b == n - 1:
                print(f"   {name}: [0][{b}] ratios even(h,K,symp) odd(K,h) = {rr}; odd bilinear K.h / psi_B psidot_B = {ro}")
    say(ok, f"{name} l={l}: all six ratios on the full {n}x{n} Gram = c_l = {c} (odd bilinear = -c_l)", f"distinct values seen: {sorted(seen)}")
    return FT, FS


if __name__ == "__main__":
    r2 = S
    profs = [R.one, S**3, S + 7 * S**2]         # seed profiles (times exp(-r^2/2)); Pi raises the degree further
    out = {}
    out["l3m0"] = run("l=3 m=0  phi=z(5z^2-3r^2)", z * (5 * z**2 - 3 * r2), 3, profs)
    out["l3m3"] = run("l=3 m=3  phi=Re(x+iy)^3", x**3 - 3 * x * y**2, 3, profs)
    out["l2"] = run("l=2 m=2  phi=x^2-y^2", x**2 - y**2, 2, profs[:2])
    out["l4"] = run("l=4 m=0  phi=35z^4-30z^2r^2+3r^4", 35 * z**4 - 30 * z**2 * r2 + 3 * r2**2, 4, profs[:2])
    # orthogonality between different m at the same l (needed to pass from one m to the whole sector space)
    a, b = out["l3m0"], out["l3m3"]
    zero = all(tpair(a[0][i], b[0][j], s) == 0 and tpair(a[1][i], b[1][j], s) == 0
               for i in range(3) for j in range(3) for s in (1, -1, 0))
    say(zero, "l=3: graviton forms vanish between m=0 and m=3 (even and odd)")
    print("DONE fails =", NFAIL[0])
    sys.exit(1 if NFAIL[0] else 0)
