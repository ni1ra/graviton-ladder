"""Known-answer validation of hb.py. Run:  python validate.py > logs/validate.log
Every case prints a deviation that should sit at rounding level for its arithmetic."""
import sys
from fractions import Fraction
import numpy as np
from flint import arb, arb_mat, ctx
import hb

np.set_printoptions(linewidth=200)
PREC = 400
hb.setprec(PREC)


def amat(a):
    a = np.atleast_2d(np.asarray(a, dtype=float))
    return arb_mat(a.shape[0], a.shape[1], [arb(float(x)) for x in a.flatten()])


def nmat(M):
    return np.array([[float(M[i, j]) for j in range(M.ncols())] for i in range(M.nrows())])


def f(x):
    return float(arb(x).mid())


print("== V1 thermal chain: modular Hamiltonian must be beta*H")
N, l, beta = 6, 1, 0.7
d, e = hb.lattice_coeffs(N, l, "D2")
K = np.diag([float(x) for x in d]) + np.diag([float(x) for x in e], 1) + np.diag([float(x) for x in e], -1)
w, V = np.linalg.eigh(K)
om = np.sqrt(w)
X = (V * (0.5 / om / np.tanh(beta * om / 2))) @ V.T
P = (V * (0.5 * om / np.tanh(beta * om / 2))) @ V.T
mo = hb.modular(amat(X), amat(P))
print("  max|M - beta K| = %.2e   max|N - beta 1| = %.2e" % (abs(nmat(mo["M"]) - beta * K).max(),
                                                          abs(nmat(mo["N"]) - beta * np.eye(N)).max()))

print("== V2 two-mode squeezed vacuum, reduced to one mode: eps = 2 log coth r")
for r in (0.3, 1.1):
    c = np.cosh(2 * r) / 2
    mo = hb.modular(amat([[c]]), amat([[c]]))
    print("  r=%.1f eps=%.15f exact=%.15f  M=%.15f N=%.15f" % (r, f(mo["eps"][0]), 2 * np.log(1 / np.tanh(r)),
                                                           f(mo["M"][0, 0]), f(mo["N"][0, 0])))
print("   two squeezed pairs, reduced modes mixed by a point transformation G")
r1, r2 = 0.4, 0.9
nu = np.array([np.cosh(2 * r1) / 2, np.cosh(2 * r2) / 2])
ep = np.array([2 * np.log(1 / np.tanh(r1)), 2 * np.log(1 / np.tanh(r2))])
G = np.array([[1.3, 0.4], [-0.2, 0.8]])
Gi = np.linalg.inv(G)
X = G @ np.diag(nu) @ G.T
P = Gi.T @ np.diag(nu) @ Gi
mo = hb.modular(amat(X), amat(P))
print("  max|M - exact| = %.2e  max|N - exact| = %.2e" % (abs(nmat(mo["M"]) - Gi.T @ np.diag(ep) @ Gi).max(),
                                                       abs(nmat(mo["N"]) - G @ np.diag(ep) @ G.T).max()))

print("== V3 one mode, two thermal states of the same oscillator: closed form")
for e1, e2 in ((0.5, 0.9), (2.0, 7.0)):
    E1, A1, A1i, _ = hb.williamson(amat([[e1]]), amat([[e1]]))
    E2, A2, A2i, _ = hb.williamson(amat([[e2]]), amat([[e2]]))
    lI = hb.affinity_I(E1, A1, A1i, E2, A2, A2i)
    lII, rel = hb.affinity_II(E1, A1, A1i, E2, A2, A2i)
    ex = np.sqrt(2 * np.sinh(e1 / 2) * 2 * np.sinh(e2 / 2)) / (2 * np.sinh((e1 + e2) / 4))
    print("  eps=(%.1f,%.1f) A_I=%.15f A_II=%.15f exact=%.15f" % (e1, e2, np.exp(f(lI)), np.exp(f(lII)), ex))


def fock_ops(dim, modes):
    a = np.diag(np.sqrt(np.arange(1, dim)), 1)
    x1 = (a + a.T) / np.sqrt(2)
    p1 = (a - a.T) / (1j * np.sqrt(2))
    I = np.eye(dim)
    xs, ps = [], []
    for m in range(modes):
        ox, op = np.array([[1.0]]), np.array([[1.0 + 0j]])
        for k in range(modes):
            ox = np.kron(ox, x1 if k == m else I)
            op = np.kron(op, p1 if k == m else I)
        xs.append(ox)
        ps.append(op)
    return xs, ps


def fock_K(M, Nn, xs, ps):
    n = len(xs)
    K = np.zeros(xs[0].shape, dtype=complex)
    for i in range(n):
        for j in range(n):
            K += 0.5 * M[i, j] * xs[i] @ xs[j] + 0.5 * Nn[i, j] * ps[i] @ ps[j]
    return (K + K.conj().T) / 2


def fock_sqrt_state(K):
    w, V = np.linalg.eigh(K)
    w = w - w.min()
    pr = np.exp(-w)
    pr /= pr.sum()
    return (V * np.sqrt(pr)) @ V.conj().T, (V * pr) @ V.conj().T


def gauss_affinity(M, Nn, Mp):
    E1, A1, A1i, _ = hb.williamson(amat(M), amat(Nn))
    E2, A2, A2i, _ = hb.williamson(amat(Mp), amat(Nn))
    lI = hb.affinity_I(E1, A1, A1i, E2, A2, A2i)
    lII, rel = hb.affinity_II(E1, A1, A1i, E2, A2, A2i)
    return np.exp(f(lI)), np.exp(f(lII))


print("== V4 brute-force Fock space, one mode: K = (N p^2 + M x^2)/2, K' = K + c x^2/2")
xs, ps = fock_ops(600, 1)
for M, Nn, c in ((1.3, 0.8, 2 * np.pi * 0.2), (1.3, 0.8, -2 * np.pi * 0.1), (3.0, 2.0, 2 * np.pi), (9.0, 6.0, -2 * np.pi)):
    s1, _ = fock_sqrt_state(fock_K(np.array([[M]]), np.array([[Nn]]), xs, ps))
    s2, _ = fock_sqrt_state(fock_K(np.array([[M + c]]), np.array([[Nn]]), xs, ps))
    Af = np.trace(s1 @ s2).real
    AI, AII = gauss_affinity([[M]], [[Nn]], [[M + c]])
    print("  M=%.1f N=%.1f c=%+.3f  Fock=%.13f  I=%.13f  II=%.13f  |I-Fock|=%.1e" % (M, Nn, c, Af, AI, AII, abs(AI - Af)))

print("== V5 brute-force Fock space, two coupled modes, diagonal phi^2 perturbation")
xs, ps = fock_ops(48, 2)
M = np.array([[4.0, -1.1], [-1.1, 3.0]])
Nn = np.array([[2.5, 0.6], [0.6, 3.5]])
for g, wv in ((1.0, (0.3, 0.3)), (-1.0, (0.3, 0.15)), (0.25, (0.3, 0.3))):
    Mp = M + 2 * np.pi * g * np.diag(wv)
    s1, _ = fock_sqrt_state(fock_K(M, Nn, xs, ps))
    s2, _ = fock_sqrt_state(fock_K(Mp, Nn, xs, ps))
    Af = np.trace(s1 @ s2).real
    AI, AII = gauss_affinity(M, Nn, Mp)
    print("  g=%+.2f w=%s  Fock=%.12f  I=%.12f  II=%.12f  |I-Fock|=%.1e" % (g, wv, Af, AI, AII, abs(AI - Af)))

print("== V6 lattice vacuum, n=2 inside N=10 (D2, l=1): exp(-K) must reproduce X_A, P_A in Fock space")
hb.setprec(600)
XA, PA, dg = hb.vacuum_blocks(2, 10, 1, "D2", 600)
mo = hb.modular(XA, PA)
print("  eps =", [f(x) for x in mo["eps"]])
Mn, Nnn = nmat(mo["M"]), nmat(mo["N"])
xs, ps = fock_ops(48, 2)
_, rho = fock_sqrt_state(fock_K(Mn, Nnn, xs, ps))
Xf = np.array([[np.trace(rho @ xs[i] @ xs[j]).real for j in range(2)] for i in range(2)])
Pf = np.array([[np.trace(rho @ ps[i] @ ps[j]).real for j in range(2)] for i in range(2)])
print("  max|X_fock - X_A| = %.2e   max|P_fock - P_A| = %.2e" % (abs(Xf - nmat(XA)).max(), abs(Pf - nmat(PA)).max()))
a = Fraction(2, 5)
for g in (1, -1):
    Mp = hb.perturbed_M(mo["M"], a, [1, 1], g)
    E1, A1, A1i, _ = hb.williamson(mo["M"], mo["N"])
    E2, A2, A2i, _ = hb.williamson(Mp, mo["N"])
    lI = hb.affinity_I(E1, A1, A1i, E2, A2, A2i)
    lII, rel = hb.affinity_II(E1, A1, A1i, E2, A2, A2i)
    s1, _ = fock_sqrt_state(fock_K(Mn, Nnn, xs, ps))
    s2, _ = fock_sqrt_state(fock_K(nmat(Mp), Nnn, xs, ps))
    Af = np.trace(s1 @ s2).real
    print("  g=%+d  1-A: Fock=%.12e  I=%.12e  II=%.12e" % (g, 1 - Af, f(hb.one_minus_A(lI)), f(hb.one_minus_A(lII))))

print("== V7 ball lattice: g=0 gives A=1; williamson(M,N) returns the entanglement spectrum; formula I = formula II")
for n, l, prec in ((12, 2, 1600), (16, 4, 2400), (24, 8, 4000)):
    hb.setprec(prec)
    XA, PA, dg = hb.vacuum_blocks(n, 4 * n, l, "D2", prec)
    mo = hb.modular(XA, PA)
    E1, A1, A1i, d1 = hb.williamson(mo["M"], mo["N"])
    deps = max(hb.l10((E1[k] - mo["eps"][n - 1 - k]) / E1[k]) for k in range(n))
    l0 = hb.affinity_I(mo["eps"], hb.mid(hb.scale_cols(mo["W"], [1 / x.sqrt() for x in mo["nu"]])),
                       hb.mid(hb.diagmat([x.sqrt() for x in mo["nu"]]) * mo["Winv"]), E1, A1, A1i)
    a = Fraction(2, 2 * n + 1)
    out = []
    for g in (1, -1):
        Mp = hb.perturbed_M(mo["M"], a, [1] * n, g)
        E2, A2, A2i, d2 = hb.williamson(Mp, mo["N"])
        lI = hb.affinity_I(E1, A1, A1i, E2, A2, A2i)
        lII, rel = hb.affinity_II(E1, A1, A1i, E2, A2, A2i)
        out.append((g, hb.one_minus_A(lI).str(25), hb.l10((lI - lII) / lI), rel))
    print("  n=%d l=%d prec=%d  log10 min(nu^2-1/4)=%.1f  log10 max rel|eps_W - eps|=%.1f  log10|1-A(g=0)|=%.1f"
          % (n, l, prec, mo["diag"]["log10_min_delta"], deps, hb.l10(hb.one_minus_A(l0))))
    for g, v, dd, rel in out:
        print("     g=%+d  1-A(I)=%s  log10 rel|I-II|=%.1f  (log10 rel radius of det in II: %.1f)" % (g, v, dd, rel))

print("== V8 Srednicki large-l entropy: S ~ xi (1 - log xi), xi = n(n+1)(2n+1)^2 / (64 l^2 (l+1)^2)")
for n, l in ((6, 150), (6, 300), (10, 400)):
    hb.setprec(3000)
    XA, PA, dg = hb.vacuum_blocks(n, 8 * n, l, "Sred", 3000)
    mo = hb.modular(XA, PA)
    xi = n * (n + 1) * (2 * n + 1) ** 2 / (64.0 * l ** 2 * (l + 1) ** 2)
    print("  n=%d l=%d  S=%.6e  xi(1-log xi)=%.6e  ratio=%.5f" % (n, l, f(mo["entropy"]), xi * (1 - np.log(xi)),
                                                              f(mo["entropy"]) / (xi * (1 - np.log(xi)))))

print("== V9 information: row sums of the pi-pi kernel N against the continuum 2 pi (1-r^2)/(2a)")
n, l, prec = 32, 2, 2600
hb.setprec(prec)
XA, PA, dg = hb.vacuum_blocks(n, 4 * n, l, "D2", prec)
mo = hb.modular(XA, PA)
a = 1.0 / (n + 0.5)
Nn = nmat(mo["N"])
for j in (4, 8, 16, 24, 28, 31, 32):
    r = j * a
    print("  j=%2d r=%.3f  rowsum N=%.4f  diag N=%.4f  2pi(1-r^2)/(2a)=%.4f" % (j, r, Nn[j - 1].sum(), Nn[j - 1, j - 1],
                                                                          2 * np.pi * (1 - r * r) / (2 * a)))
