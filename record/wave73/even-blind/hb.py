"""HINGE-BLIND core: radial-lattice vacuum covariances, modular Hamiltonian of a Gaussian state,
and the exact affinity tr(rho^{1/2} rho'^{1/2}) between exp(-K) and exp(-K - 2 pi Z).

Arbitrary precision through python-flint (arb). Ball radii are discarded after every stage
(midpoints only); accuracy is established by residual checks and by rerunning at higher
precision, not by the enclosures.

Conventions (see PREREG.md): lattice variables phi_j = phi(r_j), pi_j = a pi(r_j),
H = (1/2a)(pi.pi + phi.K.phi), X = <phi phi> = K^{-1/2}/2, P = <pi pi> = K^{1/2}/2,
K_mod = 1/2 (pi.N.pi + phi.M.phi), Z = (g/2) a sum_j w_j phi_j^2.
"""
import time
from fractions import Fraction

import numpy as np
from flint import arb, arb_mat, acb_mat, ctx
from scipy.linalg import eigh_tridiagonal


# ----------------------------------------------------------------------------- helpers
def setprec(bits):
    ctx.prec = int(bits)


def mid(M):
    return M.mid()


def A_(x):
    """arb from Fraction / int / arb, exact up to the working precision."""
    if isinstance(x, Fraction):
        return arb(x.numerator) / arb(x.denominator)
    return arb(x)


def eye(n):
    M = arb_mat(n, n)
    for i in range(n):
        M[i, i] = 1
    return M


def diagmat(v):
    n = len(v)
    M = arb_mat(n, n)
    for i in range(n):
        M[i, i] = v[i]
    return M


def scale_cols(M, v):
    """M @ diag(v)"""
    return M * diagmat(v)


def l10(x):
    """log10 |x| as a float (works far below float range); -inf for exact zero."""
    x = abs(arb(x).mid())
    if x == 0:
        return float("-inf")
    return float((x.log() / arb(10).log()).mid())


def maxabs(M):
    """log10 of the largest |entry| (midpoints)."""
    m = arb(0)
    for x in M.entries():
        a = abs(x.mid())
        if a > m:
            m = a
    return l10(m)


def maxabs_arb(M):
    m = arb(0)
    for x in M.entries():
        a = abs(x.mid())
        if a > m:
            m = a
    return m


def logdet_spd(M):
    """log det of a symmetric positive definite arb_mat via Cholesky."""
    L = cholesky(M)
    s = arb(0)
    for i in range(M.nrows()):
        s += L[i][i].log()
    return 2 * s


def cholesky(M):
    """Lower Cholesky factor as list of lists of arb (midpoints)."""
    n = M.nrows()
    L = [[arb(0)] * n for _ in range(n)]
    for i in range(n):
        Li = L[i]
        for j in range(i + 1):
            Lj = L[j]
            s = M[i, j]
            for k in range(j):
                s -= Li[k] * Lj[k]
            if i == j:
                s = s.mid()
                if not (s > 0):
                    raise ValueError("cholesky: not positive definite at row %d (pivot %s)" % (i, s.str(10)))
                Li[j] = s.sqrt().mid()
            else:
                Li[j] = (s / Lj[j]).mid()
    return L


def tri_inv_lower(L):
    """Inverse of lower-triangular list-of-lists, as arb_mat."""
    n = len(L)
    Lm = arb_mat(n, n, [x for row in L for x in row])
    return mid(Lm.inv())


def to_mat(L):
    n = len(L)
    return arb_mat(n, len(L[0]), [x for row in L for x in row])


# ----------------------------------------------------------------------------- lattice
def lattice_coeffs(N, l, disc):
    """Diagonal d_j and off-diagonal e_j (j=1..N) of the dimensionless tridiagonal K, as Fractions."""
    ll = l * (l + 1)
    d, e = [], []
    for j in range(1, N + 1):
        if disc == "D2":
            d.append(Fraction(2) + Fraction(ll, j * j))
            e.append(Fraction(-1))
        elif disc == "Sred":
            jp = Fraction(2 * j + 1, 2) ** 2
            jm = Fraction(2 * j - 1, 2) ** 2 if j >= 2 else Fraction(0)
            d.append((jp + jm + ll) / (j * j))
            e.append(-jp / (j * (j + 1)))
        else:
            raise ValueError(disc)
    return d, e[:-1]


def tridiag_eigs(dF, eF, prec, log=None):
    """All eigenvalues of the symmetric tridiagonal matrix to `prec` bits (Newton on the
    characteristic recurrence, precision ladder from float64 starting values)."""
    N = len(dF)
    lam0 = eigh_tridiagonal(np.array([float(x) for x in dF]), np.array([float(x) for x in eF]),
                            eigvals_only=True)
    lam = [arb(float(x)) for x in lam0]
    p = 64
    worst = float("-inf")
    while True:
        p = min(2 * p, prec)
        setprec(p + 30)
        d = [A_(x) for x in dF]
        e2 = [A_(x) ** 2 for x in eF]
        nsteps = 2 if p <= 128 or p == prec else 1
        for _ in range(nsteps):
            worst = float("-inf")
            for i in range(N):
                x = lam[i]
                p0, dp0 = arb(1), arb(0)
                p1, dp1 = d[0] - x, arb(-1)
                for k in range(1, N):
                    t = d[k] - x
                    p2 = (t * p1 - e2[k - 1] * p0).mid()
                    dp2 = (t * dp1 - p1 - e2[k - 1] * dp0).mid()
                    p0, dp0, p1, dp1 = p1, dp1, p2, dp2
                corr = (p1 / dp1).mid()
                lam[i] = (x - corr).mid()
                c = l10(corr / x)
                if c > worst:
                    worst = c
        if log:
            log("  eig ladder prec=%d log10 last rel corr=%.1f" % (p, worst))
        if p >= prec:
            break
    setprec(prec)
    return lam, worst


def vacuum_blocks(n, N, l, disc, prec, log=None):
    """X_A, P_A (n x n arb_mat) of the vacuum of H=(1/2)(pi.pi+phi.K.phi) restricted to sites 1..n.
    Eigenvectors by two-sided three-term recurrence glued at the float64 eigenvector's largest
    component (forward recurrence alone is unstable for the states localised at the centrifugal
    wall, whose eigenvalues lie above the band). Returns (X_A, P_A, diagnostics in log10)."""
    dF, eF = lattice_coeffs(N, l, disc)
    guard = 64
    lam, lastcorr = tridiag_eigs(dF, eF, prec + guard, log)
    dd = np.array([float(x) for x in dF]); ee = np.array([float(x) for x in eF])
    _, Vd = eigh_tridiagonal(dd, ee)
    piv = np.argmax(np.abs(Vd), axis=0)
    del Vd
    setprec(prec + guard)
    d = [A_(x) for x in dF]
    e = [A_(x) for x in eF]
    einv = [1 / x for x in e]
    cols = []
    worst_res = float("-inf")
    for i in range(N):
        x = lam[i]
        m = int(piv[i])
        v = [None] * N
        # forward 0..m
        v[0] = arb(1)
        if m >= 1:
            v[1] = ((x - d[0]) * einv[0]).mid()
        for j in range(1, m):
            v[j + 1] = (((x - d[j]) * v[j] - e[j - 1] * v[j - 1]) * einv[j]).mid()
        # backward N-1..m, then rescale to glue at m
        if m < N - 1:
            b = [None] * N
            b[N - 1] = arb(1)
            b[N - 2] = ((x - d[N - 1]) * einv[N - 2]).mid()
            for j in range(N - 2, m, -1):
                b[j - 1] = (((x - d[j]) * b[j] - e[j] * b[j + 1]) * einv[j - 1]).mid()
            s = (v[m] / b[m]).mid()
            for j in range(m + 1, N):
                v[j] = (b[j] * s).mid()
        nrm = arb(0)
        for j in range(N):
            nrm += v[j] * v[j]
        nrm = nrm.mid()
        # glue-row residual, relative to the norm
        res = (d[m] - x) * v[m]
        if m >= 1:
            res += e[m - 1] * v[m - 1]
        if m < N - 1:
            res += e[m] * v[m + 1]
        r = 0.5 * l10(res * res / nrm)
        if r > worst_res:
            worst_res = r
        s = (1 / nrm.sqrt()).mid()
        cols.append([(v[j] * s).mid() for j in range(n)])
    VA = arb_mat(n, N, [cols[k][j] for j in range(n) for k in range(N)])
    del cols
    q = [x.sqrt().sqrt().mid() for x in lam]
    B = scale_cols(VA, [1 / t for t in q])
    C = scale_cols(VA, q)
    XA = mid((B * B.transpose()) / 2)
    PA = mid((C * C.transpose()) / 2)
    I = VA * VA.transpose()
    for i in range(n):
        I[i, i] -= 1
    diag = dict(lastcorr=lastcorr, glue_res=worst_res, completeness=maxabs(I),
                lam_min=float(lam[0]), lam_max=float(lam[-1]))
    setprec(prec)
    return XA, PA, diag


# ----------------------------------------------------------------------------- symmetric eigensolver
def eigh_arb(S, check=True):
    """Eigen-decomposition of a real symmetric arb_mat: returns (list of eigenvalues, U arb_mat
    with orthonormal columns, residual diagnostics). Uses acb_mat.eig (approximate, uncertified),
    then re-orthonormalises nothing: residuals are measured and returned."""
    n = S.nrows()
    E, R = acb_mat(S).eig(right=True, nonstop=True, algorithm="approx")
    # real parts, normalise columns
    cols = []
    for k in range(n):
        col = [R[j, k] for j in range(n)]
        # pick phase: divide by the largest-modulus entry to make the vector real
        big = max(col, key=lambda z: float(abs(z).mid()))
        colr = [(z / big).real.mid() for z in col]
        nr = arb(0)
        for z in colr:
            nr += z * z
        s = (1 / nr.sqrt()).mid()
        cols.append([(z * s).mid() for z in colr])
    ev = [z.real.mid() for z in E]
    order = sorted(range(n), key=lambda k: ev[k])
    ev = [ev[k] for k in order]
    U = arb_mat(n, n, [cols[order[k]][j] for j in range(n) for k in range(n)])
    diag = {}
    if check:
        O = U.transpose() * U
        for i in range(n):
            O[i, i] -= 1
        diag["orth"] = maxabs(O)
        Rm = S * U - U * diagmat(ev)
        diag["res"] = maxabs(Rm)
    return ev, U, diag


# ----------------------------------------------------------------------------- modular Hamiltonian
def modular(XA, PA):
    """Modular Hamiltonian K = 1/2 (pi.N.pi + phi.M.phi) of the Gaussian state with <phi phi>=XA,
    <pi pi>=PA and no phi-pi correlation, from the covariances alone.
    X = W W^T, P = W^{-T} nu^2 W^{-1}; eps = log((nu+1/2)/(nu-1/2));
    N = W (eps/nu) W^T, M = W^{-T} (eps nu) W^{-1}."""
    n = XA.nrows()
    L = cholesky(XA)
    Lm = to_mat(L)
    Linv = tri_inv_lower(L)
    S = mid(Lm.transpose() * PA * Lm)
    S = mid((S + S.transpose()) / 2)
    ev, U, dg = eigh_arb(S)
    q = arb(1) / 4
    h = arb(1) / 2
    delta = [(x - q).mid() for x in ev]          # nu^2 - 1/4
    if not (delta[0] > 0):
        raise ValueError("modular: nu^2 - 1/4 not positive (%s): precision too low" % delta[0].str(10))
    nu = [x.sqrt().mid() for x in ev]
    num = [(delta[k] / (nu[k] + h)).mid() for k in range(n)]   # nu - 1/2 without cancellation
    eps = [(((nu[k] + h) / num[k]).log()).mid() for k in range(n)]
    W = mid(Lm * U)
    Winv = mid(U.transpose() * Linv)
    Nm = mid(scale_cols(W, [eps[k] / nu[k] for k in range(n)]) * W.transpose())
    Mm = mid(scale_cols(Winv.transpose(), [eps[k] * nu[k] for k in range(n)]) * Winv)
    ent = arb(0)
    for k in range(n):
        ent += (nu[k] + h) * (nu[k] + h).log() - num[k] * num[k].log()
    dg["log10_min_delta"] = l10(delta[0])
    dg["log10_max_delta"] = l10(delta[-1])
    I = W * Winv
    for i in range(n):
        I[i, i] -= 1
    dg["WWinv"] = maxabs(I)
    return dict(eps=eps, nu=nu, W=W, Winv=Winv, N=Nm, M=Mm, entropy=ent.mid(), diag=dg)


def williamson(Mm, Nm):
    """Normal modes of K = 1/2 (pi.N.pi + phi.M.phi): phi = A Q, pi = A^{-T} Pi with
    A^T M A = eps = A^{-1} N A^{-T}. Returns (eps, A, Ainv, diag). Raises if M is not positive."""
    n = Mm.nrows()
    LN = cholesky(Nm)
    LNm = to_mat(LN)
    LNinv = tri_inv_lower(LN)
    T = mid(LNm.transpose() * Mm * LNm)
    T = mid((T + T.transpose()) / 2)
    ev, U, dg = eigh_arb(T)
    if not (ev[0] > 0):
        raise ArithmeticError("williamson: M not positive on the lattice, smallest eps^2 = %s" % ev[0].str(15))
    eps = [x.sqrt().mid() for x in ev]
    rt = [x.sqrt().mid() for x in eps]
    A = mid(scale_cols(LNm * U, [1 / x for x in rt]))
    Ainv = mid(diagmat(rt) * U.transpose() * LNinv)
    I = A * Ainv
    for i in range(n):
        I[i, i] -= 1
    dg["AAinv"] = maxabs(I)
    return eps, A, Ainv, dg


def _lognorm_half(eps):
    """log tr rho^{1/2} for rho = exp(-K)/tr exp(-K): sum 1/2 log(2 sinh(eps/2)) - log(2 sinh(eps/4))."""
    s = arb(0)
    for x in eps:
        s += (2 * (x / 2).sinh()).log() / 2 - (2 * (x / 4).sinh()).log()
    return s.mid()


def _half_cov(eps, A, Ainv):
    """Covariance blocks of the normalised rho^{1/2}: Xt = A nut A^T, Pt = A^{-T} nut A^{-1}."""
    nut = [((x / 4).cosh() / (x / 4).sinh() / 2).mid() for x in eps]
    Xt = mid(scale_cols(A, nut) * A.transpose())
    Pt = mid(scale_cols(Ainv.transpose(), nut) * Ainv)
    return Xt, Pt


def affinity_I(eps1, A1, A1inv, eps2, A2, A2inv):
    """log tr(rho1^{1/2} rho2^{1/2}) by overlap of the two normalised square-root states."""
    X1, P1 = _half_cov(eps1, A1, A1inv)
    X2, P2 = _half_cov(eps2, A2, A2inv)
    ld = logdet_spd(mid(X1 + X2)) + logdet_spd(mid(P1 + P2))
    return (_lognorm_half(eps1) + _lognorm_half(eps2) - ld / 2).mid()


def _halfT(eps, A, Ainv):
    """Real form of the symplectic matrix of exp(-K/2): S [[ch, sh],[sh, ch]] S^{-1}, S=diag(A, A^{-T})."""
    n = len(eps)
    ch = [(x / 2).cosh() for x in eps]
    sh = [(x / 2).sinh() for x in eps]
    AT = A.transpose()
    AinvT = Ainv.transpose()
    b11 = scale_cols(A, ch) * Ainv
    b12 = scale_cols(A, sh) * AT
    b21 = scale_cols(AinvT, sh) * Ainv
    b22 = scale_cols(AinvT, ch) * AT
    R = arb_mat(2 * n, 2 * n)
    for i in range(n):
        for j in range(n):
            R[i, j] = b11[i, j]
            R[i, n + j] = b12[i, j]
            R[n + i, j] = b21[i, j]
            R[n + i, n + j] = b22[i, j]
    return mid(R)


def affinity_II(eps1, A1, A1inv, eps2, A2, A2inv):
    """log tr(rho1^{1/2} rho2^{1/2}) from tr(e^{-K1/2} e^{-K2/2}) = [(-1)^n det(1 - R1 R2)]^{-1/2}.
    Returns (value, relative radius of the determinant as log10)."""
    n = len(eps1)
    R = _halfT(eps1, A1, A1inv) * _halfT(eps2, A2, A2inv)
    D = -R
    for i in range(2 * n):
        D[i, i] += 1
    det = D.det()
    if n % 2:
        det = -det
    rel = l10(det.rad()) - l10(det.mid()) if det.rad() != 0 else float("-inf")
    s = arb(0)
    for x in list(eps1) + list(eps2):
        s += (2 * (x / 2).sinh()).log() / 2
    return (s - det.mid().log() / 2).mid(), rel


def perturbed_M(Mm, a, w, g):
    """M' = M + 2 pi g a diag(w)   (Z = (g/2) a sum_j w_j phi_j^2)."""
    n = Mm.nrows()
    c = 2 * arb.pi() * A_(g) * A_(a)
    Mp = arb_mat(n, n, list(Mm.entries()))
    for j in range(n):
        Mp[j, j] = (Mp[j, j] + c * A_(w[j])).mid()
    return Mp


def one_minus_A(logA):
    """1 - A from log A without cancellation."""
    return (-(logA.expm1())).mid() if hasattr(logA, "expm1") else (1 - logA.exp()).mid()


# ----------------------------------------------------------------------------- second order in g
def second_order(eps, A1, a, w):
    """c2 = lim_{g->0} (1 - A)/g^2 for rho' ∝ exp(-K - 2 pi Z), Z = (g/2) a sum_j w_j phi_j^2,
    from second-order perturbation theory of sqrt(rho) in the normal modes of K
    (phi = A1 Q, K = sum eps_k (b^+ b + 1/2)):

      1 - A = (1/2) tr[(d sqrt(rho))^2],
      c2 = (2 pi a)^2 / 8 * sum_{k,l} G_kl^2 [ (1 - e^{-(e_k+e_l)/2})^2 / (e_k+e_l)^2
                                              + (e^{-e_l/2} - e^{-e_k/2})^2 / (e_k-e_l)^2 ] / ((1-e^{-e_k})(1-e^{-e_l})),
      G = A1^T diag(w) A1, the k = l limit of the second term being e^{-e_k}/4.
    Independent of the exact-affinity code path (no second Williamson form, no determinants)."""
    n = len(eps)
    G = A1.transpose() * diagmat([A_(x) for x in w]) * A1
    ex = [(-x).exp() for x in eps]
    exh = [(-x / 2).exp() for x in eps]
    den = [1 - x for x in ex]
    s = arb(0)
    for k in range(n):
        for l in range(k, n):
            g2 = G[k, l] * G[k, l]
            sm = eps[k] + eps[l]
            t1 = (1 - exh[k] * exh[l]) ** 2 / (sm * sm)
            if k == l:
                t2 = ex[k] / 4
            else:
                df = eps[k] - eps[l]
                t2 = (exh[l] - exh[k]) ** 2 / (df * df)
            term = g2 * (t1 + t2) / (den[k] * den[l])
            s += term if k == l else 2 * term
    c = (2 * arb.pi() * A_(a)) ** 2 / 8
    return (c * s).mid()
