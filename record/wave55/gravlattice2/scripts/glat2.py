"""w55-gravlattice2 solver: the wave52 computation (glat.py, unchanged) block-diagonalised by the reflection group Z2^3 of the
cubic ball (x_k - c_k -> -(x_k - c_k)), a subgroup of the octahedral group O_h.  Not a radial reduction: each of the eight parity
sectors mixes all l and both chain parities.  Only rows X[orbit representatives, :] are built; W stays sparse.
Each reflection g acts on the fields h, pi (true tensors on primal cells) by T_ab(y) -> s_a s_b T_ab(g y); on the generators
e = inc(h), b = sc(pi) it is a signed permutation e_t(p) -> sig(t, g) e_t(g p), with the sign found from the stencils."""
import numpy as np
from itertools import product
from scipy.sparse import coo_matrix, csr_matrix
from glat import INC, SC, LAM, OFF2, TCOMP, gauss_rel

GROUP = [np.array(s) for s in product((1, -1), repeat=3)]
CHARS = [np.array(k) for k in product((0, 1), repeat=3)]


def chi(k, s):
    return int(np.prod(np.where(np.array(k) == 1, s, 1)))


def stencil_sign(st, t, s):
    """sig with  st_t(p') applied to (g.field) = sig * st_t(p) applied to field,  p' = g p (centre at 0; translation invariant)."""
    ref = {(tin, tuple(int(v) for v in sh)): c for (tin, sh), c in st[t].items()}
    mapped = {}
    for (tin, sh), c in ref.items():
        a, b = TCOMP[tin]
        mapped[(tin, tuple(int(v) for v in np.array(sh) * s))] = c * s[a] * s[b]
    # e_t(gp)(h) = sum_k ref[k] h(gp + sh);  e_t(p)(U_g h) = sum ref[k] s_a s_b h(g(p + sh)) = sum mapped[k'] h(gp + sh')
    assert set(mapped) == set(ref), 'stencil not reflection-closed'
    r = {k: mapped[k] / ref[k] for k in ref}
    vals = set(np.round(list(r.values()), 12))
    assert len(vals) == 1 and abs(abs(list(vals)[0]) - 1) < 1e-12, vals
    return int(round(list(vals)[0]))


def orbits(T, P, centre2, st):
    """images under each g: perm[g][n] = index of g-image, sig[g][n] = sign.  rep[n] = smallest index in the orbit."""
    idx = {(int(t),) + tuple(int(v) for v in p): n for n, (t, p) in enumerate(zip(T, P))}
    perm = np.zeros((8, len(T)), int); sig = np.zeros((8, len(T)), int)
    for gi, s in enumerate(GROUP):
        Pg = centre2 + s * (P - centre2)
        sgn_t = [stencil_sign(st, t, s) for t in range(6)]
        for n in range(len(T)):
            perm[gi, n] = idx[(int(T[n]),) + tuple(int(v) for v in Pg[n])]
            sig[gi, n] = sgn_t[T[n]]
    rep = perm.min(0)
    return perm, sig, rep


def sector_bases(perm, sig, rep):
    """for each character k: list of rep indices kept, and sparse V (n x m) with orthonormal columns v = sum_g chi(g) U_g e_r / N."""
    reps = np.unique(rep)
    out = []
    for k in CHARS:
        rows, cols, vals, keep, norms = [], [], [], [], []
        for r in reps:
            v = {}
            for gi, s in enumerate(GROUP):
                j = perm[gi, r]
                v[j] = v.get(j, 0.0) + chi(k, s) * sig[gi, r]
            v = {j: c for j, c in v.items() if c != 0}
            if not v:
                continue
            N = np.sqrt(sum(c * c for c in v.values()))
            m = len(keep)
            for j, c in v.items():
                rows.append(j); cols.append(m); vals.append(c / N)
            keep.append(r); norms.append(N)
        V = coo_matrix((vals, (rows, cols)), shape=(perm.shape[1], len(keep))).tocsc()
        out.append((np.array(keep, int), np.array(norms), V))
    return out, reps


def corr_rows(Ti, Pi, Tj, Pj, Garr, tabs):
    """<e_i e_j> for row cells (Ti, Pi) and column cells (Tj, Pj): glat.corr_matrix without the full square."""
    M = np.zeros((len(Ti), len(Tj)))
    for t, tp in product(range(6), repeat=2):
        I = np.where(Ti == t)[0]; J = np.where(Tj == tp)[0]
        if len(I) == 0 or len(J) == 0:
            continue
        par, w = tabs[(t, tp)]
        D2 = Pi[I][:, None, :] - Pj[J][None, :, :]
        blk = np.zeros((len(I), len(J)))
        for d, v in w.items():
            A = np.abs(D2 + np.array(d))
            assert not (A % 2).any()
            A //= 2
            blk += v * Garr[A[..., 0], A[..., 1], A[..., 2]]
        M[np.ix_(I, J)] = blk
        del D2
    return M


def sector_matrix(T, P, reps, sectors, Garr, tabs):
    """X_k = V_k^T X V_k = diag(|G| / N_a) * X[rep_a, :] @ V_k  (X V_k is chi_k-isotypic)."""
    Xr = corr_rows(T[reps], P[reps], T, P, Garr, tabs)
    pos = {int(r): a for a, r in enumerate(reps)}
    out = []
    for keep, norms, V in sectors:
        rows = np.array([pos[int(r)] for r in keep], int)
        Mk = (V.T @ Xr[rows].T).T * (8.0 / norms)[:, None]
        out.append(0.5 * (Mk + Mk.T))
    del Xr
    return out


def commutator_sparse(Te, Pe, Tb, Pb):
    """W = [e, b]/i as a sparse matrix (same construction as glat.commutator, never densified)."""
    keys = {}

    def build(T, P, st):
        rows, cols, vals = [], [], []
        for n in range(len(T)):
            for (tin, s), c in st[T[n]].items():
                kk = (tin,) + tuple(P[n] + np.array(s))
                if kk not in keys:
                    keys[kk] = len(keys)
                rows.append(n); cols.append(keys[kk]); vals.append(c)
        return rows, cols, vals
    re_, ce, ve = build(Te, Pe, INC)
    rb, cb, vb = build(Tb, Pb, SC)
    na = len(keys)
    lam = np.zeros(na)
    for kk, a in keys.items():
        lam[a] = LAM[kk[0]]
    Ie = coo_matrix((ve, (re_, ce)), shape=(len(Te), na)).tocsr()
    Sb = coo_matrix((vb, (rb, cb)), shape=(len(Tb), na)).tocsr()
    return csr_matrix(Ie.multiply(lam[None, :]) @ Sb.T)


def _reduce(X, P, W, dq, dp, xmax, pmax, tol):
    w, V = np.linalg.eigh(X); Qx = V[:, w > tol * xmax]; Nx = V[:, w <= tol * xmax]
    w, V = np.linalg.eigh(P); Qp = V[:, w > tol * pmax]; Np = V[:, w <= tol * pmax]
    Wt = Qx.T @ W @ Qp
    s = np.linalg.svd(Wt, compute_uv=False) if min(Wt.shape) else np.zeros(0)
    return Qx, Nx, Qp, Np, Wt, s


def rel_entropy_sectors(Xs, Ps, Ws, dqs, dps, caps=(25.0, 30.0, 35.0), tol=1e-10):
    """glat.rel_entropy applied sector by sector, with the null and rank thresholds relative to the global maxima."""
    xmax = max(np.linalg.eigvalsh(X).max() for X in Xs if X.size)
    pmax = max(np.linalg.eigvalsh(P).max() for P in Ps if P.size)
    smax = 0.0
    for X, P, W, dq, dp in zip(Xs, Ps, Ws, dqs, dps):
        if X.size and P.size:
            s = _reduce(X, P, W, dq, dp, xmax, pmax, tol)[5]
            if s.size:
                smax = max(smax, s.max())
    nd = dqs[0].shape[1]
    tot = {c: np.zeros(nd) for c in caps}
    info = dict(n_q=0, n_p=0, rank=0, cq=0, cp=0, sigma_min=np.inf, nu_min=np.inf, KLq=np.zeros(nd), KLp=np.zeros(nd),
                dq_null_abs=np.zeros(nd), dp_null_abs=np.zeros(nd), sectors=[])
    for X, P, W, dq, dp in zip(Xs, Ps, Ws, dqs, dps):
        if not (X.size or P.size):
            continue
        if X.size == 0 or P.size == 0:
            raise RuntimeError('one-sided sector')
        Qx, Nx, Qp, Np, Wt, s = _reduce(X, P, W, dq, dp, xmax, pmax, tol)
        if Nx.shape[1]:
            info['dq_null_abs'] = np.maximum(info['dq_null_abs'], np.abs(Nx.T @ dq).max(0))
        if Np.shape[1]:
            info['dp_null_abs'] = np.maximum(info['dp_null_abs'], np.abs(Np.T @ dp).max(0))
        Xt = Qx.T @ X @ Qx; Pt = Qp.T @ P @ Qp
        dqt = Qx.T @ dq; dpt = Qp.T @ dp
        U, s, Vt = np.linalg.svd(Wt)
        r = int((s > 1e-9 * smax).sum())
        Ur, sr, U0 = U[:, :r], s[:r], U[:, r:]
        Vv = Vt.T; Vr, V0 = Vv[:, :r], Vv[:, r:]
        q1 = (Ur / sr).T
        Xqq = q1 @ Xt @ q1.T; dq1 = q1 @ dqt
        Ppp = Vr.T @ Pt @ Vr; dp1 = Vr.T @ dpt
        KLq = np.zeros(nd); KLp = np.zeros(nd)
        if U0.shape[1]:
            Czz = U0.T @ Xt @ U0; Cqz = q1 @ Xt @ U0; dz = U0.T @ dqt
            sol = np.linalg.solve(Czz, dz)
            KLq = 0.5 * np.sum(dz * sol, axis=0)
            Xqq = Xqq - Cqz @ np.linalg.solve(Czz, Cqz.T); dq1 = dq1 - Cqz @ sol
        if V0.shape[1]:
            Czz = V0.T @ Pt @ V0; Cpz = Vr.T @ Pt @ V0; dz = V0.T @ dpt
            sol = np.linalg.solve(Czz, dz)
            KLp = 0.5 * np.sum(dz * sol, axis=0)
            Ppp = Ppp - Cpz @ np.linalg.solve(Czz, Cpz.T); dp1 = dp1 - Cpz @ sol
        if r:
            res, numin = gauss_rel(0.5 * (Xqq + Xqq.T), 0.5 * (Ppp + Ppp.T), dq1, dp1, caps)
        else:
            res, numin = {c: np.zeros(nd) for c in caps}, np.inf
        for c in caps:
            tot[c] += KLq + KLp + res[c]
        info['KLq'] += KLq; info['KLp'] += KLp
        info['n_q'] += Xt.shape[0]; info['n_p'] += Pt.shape[0]; info['rank'] += r
        info['cq'] += Wt.shape[0] - r; info['cp'] += Wt.shape[1] - r
        if r:
            info['sigma_min'] = min(info['sigma_min'], float(sr.min()))
        info['nu_min'] = min(info['nu_min'], numin)
        info['sectors'].append(dict(n_e=X.shape[0], n_b=P.shape[0], n_q=Xt.shape[0], n_p=Pt.shape[0], rank=r))
    return tot, info
