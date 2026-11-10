"""Shared exact objects for w70-modflow: the helicity-2 two-point function in (E, B)
form and the conformal action delta_X on (E, B).

Conventions (REPORT section 1):
  xi = x - y = (T, X1, X2, X3), s = xi^2 = T^2 - |X|^2 (signature + - - -).
  <phi^a_ij(x) phi^b_kl(y)> = W^{ab}_{ij,kl}(xi) = c0 * Q^{ab}(omega -> i d_T, k_j -> -i d_j) (1/s),
  c0 = -1/(4 pi^2) (dropped: every statement checked is homogeneous in W),
  boundary value T -> T - i eps.
  Q^{EE} = Q^{BB} = (1/2) omega^4 P_TT(khat), Q^{EB} = -(1/2) omega^4 (khat x P_TT), Q^{BE} = -Q^{EB}.
"""
import itertools
import sympy as sp

T, X1, X2, X3 = sp.symbols('T X1 X2 X3')
XI = (T, X1, X2, X3)
S = T**2 - X1**2 - X2**2 - X3**2
om, k1, k2, k3 = sp.symbols('omega k1 k2 k3')
KV = (k1, k2, k3)
R3 = range(3)
eps = sp.LeviCivita


def momentum_polys():
    """Q^{ab}_{ij,kl} as polynomials in (omega, k), built from the TT mode sum."""
    ksq = k1**2 + k2**2 + k3**2
    pi_ = [[ksq * (1 if i == j else 0) - KV[i] * KV[j] for j in R3] for i in R3]  # omega^2 P_ij on shell
    Pi4 = {}
    for i, j, k, l in itertools.product(R3, repeat=4):
        Pi4[i, j, k, l] = sp.expand(sp.Rational(1, 2) * (pi_[i][k] * pi_[j][l] + pi_[i][l] * pi_[j][k] - pi_[i][j] * pi_[k][l]))
    QEE, QEB, QBB, QBE = {}, {}, {}, {}
    for i, j, k, l in itertools.product(R3, repeat=4):
        # omega^4 P with omega^2 = k^2: degree 4 in k.  Write it with omega explicitly: replace k^2 by omega^2
        QEE[i, j, k, l] = sp.Rational(1, 2) * Pi4[i, j, k, l]
        # (khat x P)_{ij,kl} = eps_{iab} khat_a P_{bj,kl};  omega^4 (khat x P) = eps k_a Pi4 / omega
        e1 = sp.expand(sum(eps(i, a, b) * KV[a] * Pi4[b, j, k, l] for a in R3 for b in R3))
        q, r = sp.div(e1, ksq, *KV)
        assert r == 0, 'k x Pi4 not divisible by k^2'
        QEB[i, j, k, l] = -sp.Rational(1, 2) * om * sp.expand(q)      # omega * (degree 3 in k)
        QBE[i, j, k, l] = +sp.Rational(1, 2) * om * sp.expand(q)
    # Q^{BB} = (1/2) omega^4 (khat x)^T (khat x) P, computed explicitly (must equal Q^{EE})
    for i, j, k, l in itertools.product(R3, repeat=4):
        e2 = sp.expand(sum(eps(m, a, i) * KV[a] * eps(m, c, d) * KV[c] * Pi4[d, j, k, l]
                           for m in R3 for a in R3 for c in R3 for d in R3))
        q, r = sp.div(e2, ksq, *KV)
        assert r == 0
        QBB[i, j, k, l] = sp.Rational(1, 2) * sp.expand(q)
    return QEE, QEB, QBE, QBB


_D4 = {}


def d4_num(alpha):
    """s^5 * d^alpha (1/s) for a multi-index alpha (tuple of 4 counts), as a polynomial."""
    alpha = tuple(alpha)
    if alpha not in _D4:
        e = 1 / S
        for v, n in zip(XI, alpha):
            if n:
                e = sp.diff(e, v, n)
        order = sum(alpha)
        _D4[alpha] = sp.expand(sp.cancel(e * S**(order + 1))), order
    return _D4[alpha]


def to_xspace(Q):
    """Numerator N with W = N / s^5, from a degree-4 polynomial Q(omega, k)."""
    P = sp.Poly(sp.expand(Q), om, k1, k2, k3)
    N = 0
    for (a, b1, b2, b3), c in P.terms():
        nk = b1 + b2 + b3
        if c == 0:
            continue
        if a + nk == 4:
            pass
        else:
            # EE-type entries are written with k only (omega^2 = k^2): degree 4 in k
            assert a == 0 and nk == 4, (a, nk)
        phase = sp.I**a * (-sp.I)**nk
        num, order = d4_num((a, b1, b2, b3))
        assert order == 4
        N += c * phase * num
    N = sp.expand(N)
    assert not N.has(sp.I)
    return N


def build_W():
    """Wn[a][b][i,j,k,l]: numerators over s^5; a, b in {0: E, 1: B}."""
    QEE, QEB, QBE, QBB = momentum_polys()
    Wn = [[{}, {}], [{}, {}]]
    for key in QEE:
        Wn[0][0][key] = to_xspace(QEE[key])
        Wn[0][1][key] = to_xspace(QEB[key])
        Wn[1][0][key] = to_xspace(QBE[key])
        Wn[1][1][key] = to_xspace(QBB[key])
    return Wn


# ---------------- conformal Killing fields and the action on (E, B) ----------------

def killing_fields(xs):
    """15 conformal Killing fields X^mu as polynomials in xs = (t, x1, x2, x3)."""
    t, a, b, c = xs
    sp_ = (a, b, c)
    x2 = t**2 - a**2 - b**2 - c**2
    F = {}
    for m in range(4):
        F['P%d' % m] = tuple(1 if n == m else 0 for n in range(4))
    for (i, j) in ((0, 1), (0, 2), (1, 2)):
        v = [0, 0, 0, 0]
        v[1 + j] = sp_[i]
        v[1 + i] = -sp_[j]
        F['J%d%d' % (i + 1, j + 1)] = tuple(v)          # x_i d_j - x_j d_i
    for i in R3:
        v = [sp_[i], 0, 0, 0]
        v[1 + i] = t
        F['N%d' % (i + 1)] = tuple(v)                   # x_i d_t + t d_i
    F['D'] = tuple(xs)
    for m in range(4):
        bx = t if m == 0 else -sp_[m - 1]
        F['K%d' % m] = tuple(sp.expand(2 * bx * xs[n] - x2 * (1 if n == m else 0)) for n in range(4))
    return F


def delta_matrix_part(Xf, xs, E, B, w=3, s=+1):
    """Zeroth-order part of delta_X on (E, B): returns (dE, dB) as 3x3 lists.
    dE_ij = A_im E_mj + A_jm E_im + w lam E_ij + s (eps_ikl beta_k B_lj + eps_jkl beta_k B_li)
    dB_ij = A_im B_mj + A_jm B_im + w lam B_ij - s (eps_ikl beta_k E_lj + eps_jkl beta_k E_li)
    A_im = d_i X^m - lam delta_im, beta_k = d_k X^0, lam = d.X / 4.
    The PREREG registered s = -1 (Paper A's eq. (flow), whose B is minus wave69's B). The first run
    (out/s1_out_as_registered.json) showed s = -1 fails and s = +1 holds; default corrected post hoc."""
    lam = sum(sp.diff(Xf[m], xs[m]) for m in range(4)) / 4
    A = [[sp.diff(Xf[1 + m], xs[1 + i]) - (lam if i == m else 0) for m in R3] for i in R3]
    beta = [sp.diff(Xf[0], xs[1 + k]) for k in R3]

    def rot(M):
        return [[sum(A[i][m] * M[m][j] + A[j][m] * M[i][m] for m in R3) + w * lam * M[i][j] for j in R3] for i in R3]

    def bst(M):
        return [[sum(eps(i, k, l) * beta[k] * M[l][j] + eps(j, k, l) * beta[k] * M[l][i] for k in R3 for l in R3)
                 for j in R3] for i in R3]

    rE, rB, bE, bB = rot(E), rot(B), bst(B), bst(E)
    dE = [[rE[i][j] + s * bE[i][j] for j in R3] for i in R3]
    dB = [[rB[i][j] - s * bB[i][j] for j in R3] for i in R3]
    return dE, dB
