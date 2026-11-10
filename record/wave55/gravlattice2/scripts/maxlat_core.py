"""w47-maxlattice core: exact infinite-volume cubic-lattice vacuum correlators, ball region
algebras, Gaussian (Peschel / Casini-Huerta) relative entropy of coherent states, with centre.

No radial reduction and no chain map anywhere: the modular kernel comes only from the
field's own two-point functions on a 3D cubic lattice."""
import numpy as np
from scipy.special import ive
from scipy.integrate import quad
from itertools import product

# ---------------------------------------------------------------- propagator
_GCACHE = {}


def G_triple(d):
    """<phi_x phi_{x+d}> for the massless cubic-lattice scalar, H = 1/2 sum pi^2 + 1/2 sum (grad phi)^2.
    FT of 1/(2|khat|), |khat|^2 = sum 2(1-cos k).  1/|k| = pi^{-1/2} int t^{-1/2} e^{-t k^2},
    e^{-2t(1-cos k)} -> e^{-2t} I_d(2t).  With t = s^2: G = pi^{-1/2} int_0^inf prod ive(d_i, 2 s^2) ds."""
    d = tuple(sorted(abs(int(v)) for v in d))
    if d in _GCACHE:
        return _GCACHE[d]
    f = lambda s: ive(d[0], 2 * s * s) * ive(d[1], 2 * s * s) * ive(d[2], 2 * s * s)
    # split: [0,S1] finite, tail via substitution s = 1/w
    S1 = 20.0
    a, _ = quad(f, 0, S1, limit=400, epsabs=1e-15, epsrel=1e-13)
    g = lambda w: f(1.0 / w) / (w * w)
    b, _ = quad(g, 0, 1.0 / S1, limit=400, epsabs=1e-15, epsrel=1e-13)
    val = (a + b) / np.sqrt(np.pi)
    _GCACHE[d] = val
    return val


def G_array(Dmax):
    """Garr[i,j,k] = G(i,j,k) for 0 <= i,j,k <= Dmax."""
    Garr = np.zeros((Dmax + 1,) * 3)
    for t in product(range(Dmax + 1), repeat=3):
        if t[0] <= t[1] <= t[2]:
            v = G_triple(t)
            for p in set(__import__('itertools').permutations(t)):
                Garr[p] = v
    return Garr


def Glook(Garr, D):
    """D: (...,3) integer differences."""
    D = np.abs(D)
    return Garr[D[..., 0], D[..., 1], D[..., 2]]


E3 = np.eye(3, dtype=int)

# ---------------------------------------------------------------- geometry


def ball_sites(R, centre=(0.0, 0.0, 0.0)):
    c = np.array(centre)
    rng = np.arange(-R - 2, R + 3)
    X = np.array(np.meshgrid(rng, rng, rng, indexing='ij')).reshape(3, -1).T
    r2 = ((X - c) ** 2).sum(1)
    return X[r2 < R * R]


def key(x):
    return (int(x[0]), int(x[1]), int(x[2]))


def build_maxwell(R, centre=(0.0, 0.0, 0.0)):
    sites = ball_sites(R, centre)
    sset = {key(x) for x in sites}
    links = []  # (x, i): from x to x+e_i
    for x in sites:
        for i in range(3):
            if key(x + E3[i]) in sset:
                links.append((key(x), i))
    lidx = {l: n for n, l in enumerate(links)}

    def plaq_links(x, i, j):
        """plaquette at base x in plane (i,j), i<j, oriented by e_i x e_j = +e_k.
        boundary: (x,i)+ , (x+e_i,j)+ , (x+e_j,i)- , (x,j)-"""
        x = np.array(x)
        return [(key(x), i, +1), (key(x + E3[i]), j, +1), (key(x + E3[j]), i, -1), (key(x), j, -1)]

    planes = [(0, 1), (0, 2), (1, 2)]
    # plaquettes with all 4 links in V
    plaqs_in = []
    for x in sites:
        for (i, j) in planes:
            pl = plaq_links(x, i, j)
            if all((y, d) in lidx for (y, d, s) in pl):
                plaqs_in.append((key(x), i, j))
    # all plaquettes touching a V link (for the E correlator)
    ext = set()
    for (x, i) in links:
        for j in range(3):
            if j == i:
                continue
            a, b = min(i, j), max(i, j)
            ext.add((x, a, b))
            ext.add((key(np.array(x) - E3[j]), a, b))
    plaqs_ext = sorted(ext)
    pidx_ext = {p: n for n, p in enumerate(plaqs_ext)}
    nl = len(links)
    # curl restricted: plaqs_in x links
    from scipy.sparse import lil_matrix
    Cin = lil_matrix((len(plaqs_in), nl))
    for n, (x, i, j) in enumerate(plaqs_in):
        for (y, d, s) in plaq_links(x, i, j):
            Cin[n, lidx[(y, d)]] += s
    # curl ext: plaqs_ext x links (only V links kept: E_l = sum_p curl_pl C_p)
    Cext = lil_matrix((len(plaqs_ext), nl))
    for n, (x, i, j) in enumerate(plaqs_ext):
        for (y, d, s) in plaq_links(x, i, j):
            if (y, d) in lidx:
                Cext[n, lidx[(y, d)]] += s
    Cin, Cext = Cin.tocsr(), Cext.tocsr()
    return dict(sites=sites, links=links, plaqs_in=plaqs_in, plaqs_ext=plaqs_ext, Cin=Cin, Cext=Cext,
                plaq_links=plaq_links)


def link_pos(links):
    return np.array([np.array(x) + 0.5 * E3[i] for (x, i) in links])


def plaq_pos(plaqs):
    return np.array([np.array(x) + 0.5 * (E3[i] + E3[j]) for (x, i, j) in plaqs])


def plaq_normal(plaqs):
    return np.array([3 - i - j for (x, i, j) in plaqs])


def same_dir_G(Garr, pos, dirs):
    """matrix delta_{dir} G(pos_a - pos_b); positions of same-dir objects differ by integers."""
    n = len(pos)
    M = np.zeros((n, n))
    for d in np.unique(dirs):
        idx = np.where(dirs == d)[0]
        P = np.rint(pos[idx] * 2).astype(np.int32)
        D = np.abs(P[:, None, :] - P[None, :, :]) // 2
        M[np.ix_(idx, idx)] = Garr[D[..., 0], D[..., 1], D[..., 2]]
        del D
    return M


# ---------------------------------------------------------------- Gaussian relative entropy

def _sqrtm_psd(M):
    w, V = np.linalg.eigh(M)
    w = np.clip(w, 0, None)
    s = np.sqrt(w)
    return (V * s) @ V.T, (V / s) @ V.T


def gauss_rel(X, P, dq, dp, caps=(25.0, 30.0, 35.0)):
    """canonical pairs [q_a,p_b] = i delta_ab, vacuum covariances X=<qq>, P=<pp>, <{q,p}>=0.
    Returns S = 1/2 sum eps_k (dq'_k^2 + dp'_k^2) for each eps cap, plus nu_min."""
    Xh, Xhi = _sqrtm_psd(X)
    M = Xh @ P @ Xh
    M = 0.5 * (M + M.T)
    w, U = np.linalg.eigh(M)
    nu = np.sqrt(np.clip(w, 0, None))
    nu_min = nu.min()
    sn = np.sqrt(nu) if dq.ndim == 1 else np.sqrt(nu)[:, None]
    dq1 = sn * (U.T @ (Xhi @ dq))
    dp1 = (U.T @ (Xh @ dp)) / sn
    out = {}
    x = np.clip(2 * nu - 1, 1e-300, None)
    eps_raw = np.log((2 * nu + 1) / x)
    for c in caps:
        eps = np.minimum(eps_raw, c)
        if dq1.ndim == 1:
            out[c] = 0.5 * float(np.sum(eps * (dq1 ** 2 + dp1 ** 2)))
        else:
            out[c] = 0.5 * np.sum(eps[:, None] * (dq1 ** 2 + dp1 ** 2), axis=0)
    return out, float(nu_min)


def null_project(C, tol=1e-10):
    """orthonormal basis of the range of PSD C (drop zero-variance directions)."""
    w, V = np.linalg.eigh(0.5 * (C + C.T))
    keep = w > tol * w.max()
    return V[:, keep], V[:, ~keep]


def rel_entropy_with_centre(X, P, W, dq, dp, caps=(25.0, 30.0, 35.0), tol=1e-10):
    """q-type ops (commuting), p-type ops (commuting), [q_a, p_b] = i W_ab.
    Returns dict with S_full (classical KL of centre + conditioned quantum part) and S_drop."""
    info = {}
    Qx, Nx = null_project(X, tol)
    Qp, Np = null_project(P, tol)
    info['dq_null'] = float(np.abs(Nx.T @ dq).max()) if Nx.shape[1] else 0.0
    info['dp_null'] = float(np.abs(Np.T @ dp).max()) if Np.shape[1] else 0.0
    info['n_q'], info['n_p'] = Qx.shape[1], Qp.shape[1]
    Xt = Qx.T @ X @ Qx
    Pt = Qp.T @ P @ Qp
    Wt = Qx.T @ W @ Qp
    dqt = Qx.T @ dq
    dpt = Qp.T @ dp
    U, s, Vt = np.linalg.svd(Wt)
    r = int((s > 1e-9 * s.max()).sum())
    info['rank'] = r
    info['left_null'] = Wt.shape[0] - r
    info['right_null'] = Wt.shape[1] - r
    info['sigma_min'] = float(s[:r].min())
    Ur, sr = U[:, :r], s[:r]
    V = Vt.T
    Vr, V0 = V[:, :r], V[:, r:]
    # canonical non-central pairs
    q1 = (Ur / sr).T  # q' = Sigma^-1 U^T q
    Xq = q1 @ Xt @ q1.T
    dq1 = q1 @ dqt
    Ppp = Vr.T @ Pt @ Vr
    dp1 = Vr.T @ dpt
    if info['left_null'] > 0:
        # central q-type directions: should not occur for the electric algebra
        UL = U[:, r:]
        info['dq_leftnull'] = float(np.abs(UL.T @ dqt).max())
    res_drop, numin = gauss_rel(Xq, Ppp, dq1, dp1, caps)
    info['nu_min'] = numin
    if V0.shape[1] > 0:
        Pzz = V0.T @ Pt @ V0
        Ppz = Vr.T @ Pt @ V0
        dz = V0.T @ dpt
        Pzzi_dz = np.linalg.solve(Pzz, dz)
        KL = 0.5 * np.sum(dz * Pzzi_dz, axis=0)
        Pcond = Ppp - Ppz @ np.linalg.solve(Pzz, Ppz.T)
        dpc = dp1 - Ppz @ Pzzi_dz
        res_q, numin2 = gauss_rel(Xq, Pcond, dq1, dpc, caps)
        info['eig_Pzz_min'] = float(np.linalg.eigvalsh(Pzz).min())
    else:
        KL = 0.0
        res_q = res_drop
    info['KL'] = KL
    full = {c: KL + res_q[c] for c in caps}
    return full, res_drop, info


# ---------------------------------------------------------------- continuum quadrature on unit ball

def ball_quad(nr=24, nt=24, nphi=48):
    xr, wr = np.polynomial.legendre.leggauss(nr)
    r = 0.5 * (xr + 1)
    wr = 0.5 * wr * r ** 2
    ct, wt = np.polynomial.legendre.leggauss(nt)
    ph = 2 * np.pi * np.arange(nphi) / nphi
    wph = np.full(nphi, 2 * np.pi / nphi)
    st = np.sqrt(1 - ct ** 2)
    R_, C_, P_ = np.meshgrid(r, ct, ph, indexing='ij')
    W_ = (wr[:, None, None] * wt[None, :, None] * wph[None, None, :])
    S_ = np.sqrt(1 - C_ ** 2)
    pts = np.stack([R_ * S_ * np.cos(P_), R_ * S_ * np.sin(P_), R_ * C_], -1).reshape(-1, 3)
    # sphere
    Cs, Ps = np.meshgrid(ct, ph, indexing='ij')
    Ss = np.sqrt(1 - Cs ** 2)
    spts = np.stack([Ss * np.cos(Ps), Ss * np.sin(Ps), Cs], -1).reshape(-1, 3)
    sw = (wt[:, None] * wph[None, :]).reshape(-1)
    return pts, W_.reshape(-1), spts, sw


GL4 = np.polynomial.legendre.leggauss(6)


def line_int(F, x0, d):
    """int_0^1 F(x0 + s d) . d ds, F maps (n,3)->(n,3); x0 (m,3), d (m,3)."""
    xs, ws = GL4
    s = 0.5 * (xs + 1)
    w = 0.5 * ws
    tot = np.zeros(len(x0))
    for sk, wk in zip(s, w):
        tot += wk * np.einsum('ij,ij->i', F(x0 + sk * d), d)
    return tot
