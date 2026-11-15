"""EVEN-BLIND additions to hb.py (hb.py itself is an unmodified copy of hinge-blind's).

Z = (g/2) a pi^T G pi  ->  N' = N + 2 pi g a G, M' = M.
G conventions (see PREREG.md): ghost, trunc, site; plants: neumann, nocent; local: diag((1-r^2)/2)/a^2.
"""
from fractions import Fraction
from flint import arb, arb_mat
import hb


def kb_matrix(n, l, conv):
    """(a, G) with G the n x n arb_mat of the lattice Dirichlet Green operator (midpoints)."""
    ll = l * (l + 1)
    if conv == "local":
        a = Fraction(2, 2 * n + 1)
        G = arb_mat(n, n)
        for j in range(1, n + 1):
            r = j * a
            G[j - 1, j - 1] = hb.A_((1 - r * r) / 2 / (a * a))
        return a, G
    if conv == "site":
        a, m = Fraction(1, n), n - 1
    else:
        a, m = Fraction(2, 2 * n + 1), n
    cent = 0 if conv == "nocent" else ll
    K = arb_mat(m, m)
    for j in range(1, m + 1):
        K[j - 1, j - 1] = hb.A_(Fraction(2) + Fraction(cent, j * j))
        if j < m:
            K[j - 1, j] = -1
            K[j, j - 1] = -1
    if conv in ("ghost", "nocent"):
        K[m - 1, m - 1] = K[m - 1, m - 1] + 1
    elif conv == "neumann":
        K[m - 1, m - 1] = K[m - 1, m - 1] - 1
    elif conv not in ("trunc", "site"):
        raise ValueError(conv)
    Gi = hb.mid(K.inv())
    if m == n:
        return a, Gi
    G = arb_mat(n, n)
    for i in range(m):
        for j in range(m):
            G[i, j] = Gi[i, j]
    return a, G


def perturbed_N(Nm, a, G, g):
    c = 2 * arb.pi() * hb.A_(g) * hb.A_(a)
    Np = hb.mid(Nm + G * c)
    return hb.mid((Np + Np.transpose()) / 2)


def second_order_mat(eps, Gt, a):
    """hb.second_order with the mode-basis matrix supplied: field perturbation Gt = A^T W A,
    momentum perturbation Gt = A^{-1} G A^{-T}."""
    n = len(eps)
    ex = [(-x).exp() for x in eps]
    exh = [(-x / 2).exp() for x in eps]
    den = [1 - x for x in ex]
    s = arb(0)
    for k in range(n):
        for l in range(k, n):
            g2 = Gt[k, l] * Gt[k, l]
            sm = eps[k] + eps[l]
            t1 = (1 - exh[k] * exh[l]) ** 2 / (sm * sm)
            if k == l:
                t2 = ex[k] / 4
            else:
                df = eps[k] - eps[l]
                t2 = (exh[l] - exh[k]) ** 2 / (df * df)
            term = g2 * (t1 + t2) / (den[k] * den[l])
            s += term if k == l else 2 * term
    c = (2 * arb.pi() * hb.A_(a)) ** 2 / 8
    return (c * s).mid()


def q2_momentum(eps, A1i, a, G):
    Gt = hb.mid(A1i * G * A1i.transpose())
    return second_order_mat(eps, Gt, a)
