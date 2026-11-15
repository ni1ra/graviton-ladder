"""EVEN-BLIND known-answer checks for the momentum-block perturbation.
E1: brute-force Fock space, two coupled modes, K' = K + (c/2) p^T G p, against affinity_I/II.
E2: second-order formula (even.second_order_mat with Gt = A^-1 G A^-T) against exact affinity at small g.
E3: lattice G (ghost, trunc, site) against the continuum Green function of -d^2/dr^2 + l(l+1)/r^2.
E4: odd second-order through second_order_mat equals hb.second_order (same matrix path)."""
from fractions import Fraction
import numpy as np
from flint import arb, arb_mat
import hb, even

hb.setprec(600)


def amat(a):
    a = np.atleast_2d(np.asarray(a, dtype=float))
    return arb_mat(a.shape[0], a.shape[1], [arb(float(x)) for x in a.flatten()])


def nmat(M):
    return np.array([[float(M[i, j]) for j in range(M.ncols())] for i in range(M.nrows())])


def f(x):
    return float(arb(x).mid())


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
        xs.append(ox); ps.append(op)
    return xs, ps


def fock_K(M, Nn, xs, ps):
    n = len(xs)
    K = np.zeros(xs[0].shape, dtype=complex)
    for i in range(n):
        for j in range(n):
            K += 0.5 * M[i, j] * xs[i] @ xs[j] + 0.5 * Nn[i, j] * ps[i] @ ps[j]
    return (K + K.conj().T) / 2


def fock_sqrt(K):
    w, V = np.linalg.eigh(K)
    pr = np.exp(-(w - w.min())); pr /= pr.sum()
    return (V * np.sqrt(pr)) @ V.conj().T


print("== E1 Fock space, two modes, momentum-block perturbation N' = N + 2 pi g a G")
xs, ps = fock_ops(48, 2)
M = np.array([[4.0, -1.1], [-1.1, 3.0]])
Nn = np.array([[2.5, 0.6], [0.6, 3.5]])
G = np.array([[0.9, 0.35], [0.35, 0.6]])
a = 0.3
E1, A1, A1i, _ = hb.williamson(amat(M), amat(Nn))
for g in (1.0, -1.0, 0.25):
    Np = Nn + 2 * np.pi * g * a * G
    E2, A2, A2i, _ = hb.williamson(amat(M), amat(Np))
    AI = np.exp(f(hb.affinity_I(E1, A1, A1i, E2, A2, A2i)))
    AII = np.exp(f(hb.affinity_II(E1, A1, A1i, E2, A2, A2i)[0]))
    Af = np.trace(fock_sqrt(fock_K(M, Nn, xs, ps)) @ fock_sqrt(fock_K(M, Np, xs, ps))).real
    print("  g=%+.2f Fock=%.12f I=%.12f II=%.12f |I-Fock|=%.1e" % (g, Af, AI, AII, abs(AI - Af)))

print("== E2 second-order formula, momentum perturbation, against exact affinity at small g")
q2 = even.q2_momentum(E1, A1i, Fraction(3, 10), amat(G))
for ge in (Fraction(1, 1000), Fraction(-1, 1000), Fraction(1, 100000)):
    Np = even.perturbed_N(amat(Nn), Fraction(3, 10), amat(G), ge)
    E2, A2, A2i, _ = hb.williamson(amat(M), Np)
    oma = hb.one_minus_A(hb.affinity_I(E1, A1, A1i, E2, A2, A2i))
    print("  g=%s (1-A)/g^2=%.12e  q2=%.12e  rel diff=%.1e" % (ge, f(oma) / float(ge) ** 2, f(q2), f(oma) / float(ge) ** 2 / f(q2) - 1))

print("== E3 lattice G against the continuum Green function, G_c = (r<^(l+1) r>^(-l) - (r r')^(l+1))/(2l+1)")
for n in (16, 64):
    for l in (2, 8):
        row = []
        for conv in ("ghost", "trunc", "site"):
            aa, Gm = even.kb_matrix(n, l, conv)
            aa = float(aa)
            worst = 0.0
            for (i, j) in ((n // 2, n // 2), (n // 2, 3 * n // 4), (3 * n // 4, 7 * n // 8), (n - 2, n - 2)):
                r1, r2 = i * aa, j * aa
                lo, hi = min(r1, r2), max(r1, r2)
                gc = (lo ** (l + 1) * hi ** (-l) - (r1 * r2) ** (l + 1)) / (2 * l + 1)
                gl = aa * f(Gm[i - 1, j - 1])      # U = a K^-1 pi/a ... kernel G_c(r,r') dr' <-> a * (K^-1) entries / a * a
                worst = max(worst, abs(gl - gc) / abs(gc))
            row.append("%s %.2e" % (conv, worst))
        print("  n=%d l=%d  max rel |a K_B^-1 - G_c|: %s" % (n, l, "  ".join(row)))

print("== E4 odd second order: second_order_mat(A^T W A) against hb.second_order")
hb.setprec(1400)
XA, PA, dg = hb.vacuum_blocks(12, 48, 2, "D2", 1400)
mo = hb.modular(XA, PA)
A1 = hb.mid(hb.scale_cols(mo["W"], [1 / x.sqrt() for x in mo["nu"]]))
aa = Fraction(2, 25)
s1 = hb.second_order(mo["eps"], A1, aa, [1] * 12)
s2 = even.second_order_mat(mo["eps"], hb.mid(A1.transpose() * A1), aa)
print("  log10 rel diff = %.1f" % hb.l10((s1 - s2) / s1))
