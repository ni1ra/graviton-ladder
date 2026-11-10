"""w52-gravlattice core: staggered 3D cubic lattice linearised gravity, curvature-net ball algebra, exact
infinite-volume correlators from the scalar lattice propagator, Gaussian relative entropy with e- and b-type centres.
Cells, doubled coordinates: component t = (a,b) in TCOMP; primal offset off2 = 0 (a=b) or e_a+e_b (a!=b); dual = primal + (1,1,1).
Derivative d_k: f(y + e_k/2) - f(y - e_k/2) (symbol i khat_k).  No radial reduction and no chain map anywhere."""
import numpy as np
from itertools import product
from maxlat_core import G_array  # exact lattice propagator <phi phi>, FT of 1/(2|khat|) (wave47/maxlattice, unchanged)

TCOMP = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]
TIDX = {}
for n, (a, b) in enumerate(TCOMP):
    TIDX[(a, b)] = n; TIDX[(b, a)] = n
E3 = np.eye(3, dtype=int)
OFF2 = np.array([[0, 0, 0] if a == b else list(E3[a] + E3[b]) for (a, b) in TCOMP])
EPS = np.zeros((3, 3, 3))
EPS[0, 1, 2] = EPS[1, 2, 0] = EPS[2, 0, 1] = 1
EPS[0, 2, 1] = EPS[2, 1, 0] = EPS[1, 0, 2] = -1
LAM = np.array([1.0, 1.0, 1.0, 0.5, 0.5, 0.5])  # [q_t, p_t] = i LAM_t (matrix-entry pairing, symmetric variables)
KMAT = np.array([[0.5 * ((a == c) * (b == d) + (a == d) * (b == c)) - 0.5 * (a == b) * (c == d) for (c, d) in TCOMP] for (a, b) in TCOMP])


def _addterm(st, tout, tin, shift, coef):
    k = (tin, tuple(int(v) for v in shift))
    st[tout][k] = st[tout].get(k, 0.0) + coef


def stencil_inc():
    """inc(q)_ij = eps_ikl eps_jmn d_k d_m q_ln ; returns st[tout] = {(tin, shift2): coef}."""
    st = [dict() for _ in range(6)]
    for tout, (i, j) in enumerate(TCOMP):
        for k, l, m, n in product(range(3), repeat=4):
            c = EPS[i, k, l] * EPS[j, m, n]
            if c == 0:
                continue
            tin = TIDX[(l, n)]
            for sk, sm in product((1, -1), repeat=2):
                _addterm(st, tout, tin, sk * E3[k] + sm * E3[m], c * sk * sm)
    return [{k: v for k, v in d.items() if abs(v) > 1e-14} for d in st]


def stencil_sc():
    """sc(p)_ij = 1/2 (eps_ikl d_k p_lj + eps_jkl d_k p_li)."""
    st = [dict() for _ in range(6)]
    for tout, (i, j) in enumerate(TCOMP):
        for k, l in product(range(3), repeat=2):
            for (a, b) in ((i, j), (j, i)):
                c = 0.5 * EPS[a, k, l]
                if c == 0:
                    continue
                tin = TIDX[(l, b)]
                for sk in (1, -1):
                    _addterm(st, tout, tin, sk * E3[k], c * sk)
    return [{k: v for k, v in d.items() if abs(v) > 1e-14} for d in st]


INC, SC = stencil_inc(), stencil_sc()


def cells(R, centre2, dual, Rpad=0.0):
    """all cells (comp, pos2) with |pos - centre| < R + Rpad; pos2 doubled integer coords."""
    out_t, out_p = [], []
    rng = np.arange(-2 * int(R + Rpad) - 6, 2 * int(R + Rpad) + 7)
    G = np.array(np.meshgrid(rng, rng, rng, indexing='ij')).reshape(3, -1).T
    for t in range(6):
        o = OFF2[t] + (1 if dual else 0)
        P = G[np.all((G - o) % 2 == 0, axis=1)]
        d2 = ((P - centre2) ** 2).sum(1)
        P = P[d2 < (2 * (R + Rpad)) ** 2]
        out_t.append(np.full(len(P), t)); out_p.append(P)
    return np.concatenate(out_t), np.concatenate(out_p)


def stencil_inputs_inside(T, P, st, inset):
    keep = np.ones(len(T), bool)
    for n in range(len(T)):
        for (tin, s) in st[T[n]]:
            if (tin, tuple(P[n] + np.array(s))) not in inset:
                keep[n] = False
                break
    return keep


def kernel_tables(Dm):
    """kappa[(t,t')] : function of m (integer half-difference grid) for the correlator of inc-type generators.
    <e_c e_d> = sum_{(ti,s),(ti',s')} c c' K_{ti ti'} G((D2 + s - s')/2),  D2 = pos2_c - pos2_d."""
    Garr = G_array(Dm + 4)
    tabs = {}
    for t, tp in product(range(6), repeat=2):
        par = (OFF2[t] - OFF2[tp]) % 2
        w = {}
        for (ti, s), c in INC[t].items():
            for (tj, sp_), cp in INC[tp].items():
                k = KMAT[ti, tj]
                if k == 0:
                    continue
                d = tuple(np.array(s) - np.array(sp_))
                w[d] = w.get(d, 0.0) + c * cp * k
        tabs[(t, tp)] = (par, {d: v for d, v in w.items() if abs(v) > 1e-14})
    return Garr, tabs


def corr_matrix(T, P, Garr, tabs):
    n = len(T)
    M = np.zeros((n, n))
    for t, tp in product(range(6), repeat=2):
        I = np.where(T == t)[0]; J = np.where(T == tp)[0]
        if len(I) == 0 or len(J) == 0:
            continue
        par, w = tabs[(t, tp)]
        D2 = P[I][:, None, :] - P[J][None, :, :]
        blk = np.zeros((len(I), len(J)))
        for d, v in w.items():
            A = np.abs(D2 + np.array(d))
            assert not (A % 2).any()
            A //= 2
            blk += v * Garr[A[..., 0], A[..., 1], A[..., 2]]
        M[np.ix_(I, J)] = blk
        del D2
    return 0.5 * (M + M.T)


def apply_stencil_fn(T, P, st, fn_in):
    """value at output cells: sum coef * fn_in(tin, pos2_in) ; fn_in vectorised over positions."""
    out = np.zeros(len(T))
    for t in range(6):
        I = np.where(T == t)[0]
        if len(I) == 0:
            continue
        for (tin, s), c in st[t].items():
            out[I] += c * fn_in(tin, P[I] + np.array(s))
    return out


def commutator(Te, Pe, Tb, Pb):
    """W_cd = [e_c, b_d]/i = sum_a INC_{c,a} LAM_a SC_{d,a} (input cells a: primal)."""
    from scipy.sparse import coo_matrix
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
    return (Ie.multiply(lam[None, :]) @ Sb.T).toarray()


# ------------------------------------------------------------ Gaussian relative entropy with centres on both sides

def null_project(C, tol=1e-10):
    w, V = np.linalg.eigh(0.5 * (C + C.T))
    keep = w > tol * w.max()
    return V[:, keep], V[:, ~keep]


def gauss_rel(X, P, dq, dp, caps=(25.0, 30.0, 35.0)):
    """dq, dp: (n, ndata). Returns {cap: array(ndata)}, nu_min."""
    w, V = np.linalg.eigh(X)
    w = np.clip(w, 0, None)
    Xh = (V * np.sqrt(w)) @ V.T; Xhi = (V / np.sqrt(w)) @ V.T
    M = Xh @ P @ Xh; M = 0.5 * (M + M.T)
    del V
    w2, U = np.linalg.eigh(M)
    nu = np.sqrt(np.clip(w2, 0, None))
    sn = np.sqrt(nu)[:, None]
    dq1 = sn * (U.T @ (Xhi @ dq)); dp1 = (U.T @ (Xh @ dp)) / sn
    x = np.clip(2 * nu - 1, 1e-300, None)
    er = np.log((2 * nu + 1) / x)
    return {c: 0.5 * np.sum(np.minimum(er, c)[:, None] * (dq1 ** 2 + dp1 ** 2), axis=0) for c in caps}, float(nu.min())


def rel_entropy(X, P, W, dq, dp, caps=(25.0, 30.0, 35.0), tol=1e-10):
    """q-type (e), p-type (b) generators with [q_a, p_b] = i W_ab; dq (n_e, ndata), dp (n_b, ndata).
    Null directions projected; central directions on either side give a classical KL; quantum part conditioned."""
    info = {}
    Qx, Nx = null_project(X, tol); Qp, Np = null_project(P, tol)
    nq = np.abs(dq).max(0) + 1e-300; npp = np.abs(dp).max(0) + 1e-300
    info['dq_null_rel'] = (np.abs(Nx.T @ dq).max(0) / nq * (np.abs(dq).max(0) > 0)).tolist() if Nx.shape[1] else [0.0] * dq.shape[1]
    info['dp_null_rel'] = (np.abs(Np.T @ dp).max(0) / npp * (np.abs(dp).max(0) > 0)).tolist() if Np.shape[1] else [0.0] * dp.shape[1]
    del Nx, Np
    Xt = Qx.T @ X @ Qx; Pt = Qp.T @ P @ Qp; Wt = Qx.T @ W @ Qp
    dqt = Qx.T @ dq; dpt = Qp.T @ dp
    del Qx, Qp
    U, s, Vt = np.linalg.svd(Wt)
    r = int((s > 1e-9 * s.max()).sum())
    info.update(n_q=Xt.shape[0], n_p=Pt.shape[0], rank=r, cq=Wt.shape[0] - r, cp=Wt.shape[1] - r, sigma_min=float(s[:r].min()))
    Ur, sr, U0 = U[:, :r], s[:r], U[:, r:]
    V = Vt.T; Vr, V0 = V[:, :r], V[:, r:]
    q1 = (Ur / sr).T
    Xqq = q1 @ Xt @ q1.T; dq1 = q1 @ dqt
    Ppp = Vr.T @ Pt @ Vr; dp1 = Vr.T @ dpt
    nd = dq.shape[1]
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
    resq, numin = gauss_rel(0.5 * (Xqq + Xqq.T), 0.5 * (Ppp + Ppp.T), dq1, dp1, caps)
    info.update(KLq=KLq.tolist(), KLp=KLp.tolist(), nu_min=numin)
    return {c: KLq + KLp + resq[c] for c in caps}, info
