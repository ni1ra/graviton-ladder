"""S2: lattice rank count for Lemma B (PREREG section 3, S2).

D_w(V) = dim Z(V) - dim B_w(V) on Z^3 with forward differences, where
  Z(V)   = constrained fields T supported on V (div T = 0, plus the algebraic constraints built into the components)
  B_w(V) = { C f : f supported in the w-thickening V_w, and C f supported in V }.
dim Z(V)   = n_T |V| - rank(Div|_V)
dim B_w(V) = rank(C I_w) - rank(R_{outside V} C I_w)

All operators are polynomials in the commuting forward differences d_i = S_i - 1, so div o C = 0 exactly.
Ranks: exact Gaussian elimination mod p when small; otherwise eigenvalues of the exact integer Gram matrix
(float64), with the spectral gap reported. Where both run, they are compared.
"""
import itertools
import json
import sys
import time
import numpy as np
import scipy.linalg as sla
import scipy.sparse as sps
import sympy as sp
from sympy import LeviCivita as eps

P_MOD = 2147483629  # prime < 2^31
EXACT_LIMIT = int(sys.argv[2]) if len(sys.argv) > 2 else 2600  # max min(m,n) for the mod-p elimination

d = sp.symbols("d1 d2 d3")
t = sp.symbols("t1 t2 t3")  # shifts, d_i = t_i - 1
R3 = range(3)
I3 = sp.eye(3)


def curl1(A):
    return sp.Matrix(3, 3, lambda i, j: sum(eps(i, k, l) * d[k] * A[l, j] for k in R3 for l in R3))


def inc(s):
    return sp.Matrix(3, 3, lambda i, j: sum(
        eps(i, k, l) * eps(j, m, n) * d[k] * d[m] * s[l, n] for k in R3 for l in R3 for m in R3 for n in R3))


def cotton2(s):  # 2 * Cotton-York, integer coefficients
    G = inc(s)
    return curl1(2 * G - I3 * G.trace())


def stencil(expr, var):
    """coefficient of var in expr as {offset: int} with d_i -> t_i - 1"""
    c = sp.expand(sp.diff(expr, var).subs({d[i]: t[i] - 1 for i in R3}))
    if c == 0:
        return {}
    poly = sp.Poly(c, *t)
    return {tuple(int(e) for e in mon): int(co) for mon, co in poly.terms()}


def build_complex(name):
    """returns (n_pot, n_fld, Cst, Dst): stencils Cst[(out,in)], Dst[(eq,in)]"""
    if name == "maxwell":
        a = sp.symbols("a0:3")
        T = [sum(eps(i, k, l) * d[k] * a[l] for k in R3 for l in R3) for i in R3]
        Tv = sp.symbols("T0:3")
        div = [sum(d[i] * Tv[i] for i in R3)]
        pot, fld, out = list(a), list(Tv), T
    else:
        if name == "calabi":
            q = sp.symbols("q0:6")
            mk = lambda u: sp.Matrix([[u[0], u[3], u[4]], [u[3], u[1], u[5]], [u[4], u[5], u[2]]])
            read = lambda M: [M[0, 0], M[1, 1], M[2, 2], M[0, 1], M[0, 2], M[1, 2]]
            op = inc
        elif name == "cotton":
            q = sp.symbols("q0:5")
            mk = lambda u: sp.Matrix([[u[0], u[2], u[3]], [u[2], u[1], u[4]], [u[3], u[4], -u[0] - u[1]]])
            read = lambda M: [M[0, 0], M[1, 1], M[0, 1], M[0, 2], M[1, 2]]
            op = cotton2
        else:
            raise ValueError(name)
        S = mk(q)
        Tm = op(S).expand()
        assert (Tm - Tm.T).expand() == sp.zeros(3, 3)
        if name == "cotton":
            assert sp.expand(Tm.trace()) == 0
        out = read(Tm)
        Tv = sp.symbols(f"T0:{len(q)}")
        TM = mk(Tv)
        div = [sum(d[i] * TM[i, j] for i in R3) for j in R3]
        # exact check div o C = 0
        chk = [sp.expand(sum(d[i] * Tm[i, j] for i in R3)) for j in R3]
        assert all(c == 0 for c in chk)
        pot, fld = list(q), list(Tv)
    Cst = {(o, i): stencil(out[o], pot[i]) for o in range(len(out)) for i in range(len(pot))}
    Dst = {(e, i): stencil(div[e], fld[i]) for e in range(len(div)) for i in range(len(fld))}
    return len(pot), len(fld), Cst, Dst


def assemble(st, n_in, n_out, cols_sites):
    """sparse integer matrix of the stencil operator restricted to input fields supported on cols_sites.
    (Op f)(x) = sum_off coef * f(x + off). Returns (matrix, row_sites array)."""
    cols_sites = np.asarray(cols_sites, dtype=np.int64)
    offs = sorted({o for s_ in st.values() for o in s_})
    rows = set()
    for o in offs:
        rows.update(map(tuple, cols_sites - np.array(o)))
    rows = np.array(sorted(rows), dtype=np.int64)
    ridx = {tuple(r): k for k, r in enumerate(rows)}
    ri, ci, vv = [], [], []
    for (o_c, i_c), s_ in st.items():
        for off, co in s_.items():
            rs = cols_sites - np.array(off)
            rk = np.fromiter((ridx[tuple(r)] for r in rs), dtype=np.int64, count=len(rs))
            ri.append(rk * n_out + o_c)
            ci.append(np.arange(len(cols_sites)) * n_in + i_c)
            vv.append(np.full(len(rs), co, dtype=np.int64))
    M = sps.csr_matrix((np.concatenate(vv), (np.concatenate(ri), np.concatenate(ci))),
                       shape=(len(rows) * n_out, len(cols_sites) * n_in))
    return M, rows


def rank_modp(M):
    A = np.array(M.todense(), dtype=np.int64) % P_MOD
    m, n = A.shape
    if m < n:
        A = A.T.copy()
        m, n = n, m
    r = 0
    for c in range(n):
        if r == m:
            break
        nz = np.nonzero(A[r:, c])[0]
        if len(nz) == 0:
            continue
        p = r + nz[0]
        if p != r:
            A[[r, p]] = A[[p, r]]
        inv = pow(int(A[r, c]), P_MOD - 2, P_MOD)
        A[r, c:] = (A[r, c:] * inv) % P_MOD
        tgt = r + 1 + np.nonzero(A[r + 1:, c])[0]
        if len(tgt):
            A[tgt, c:] = (A[tgt, c:] - np.outer(A[tgt, c], A[r, c:]) % P_MOD) % P_MOD
        r += 1
    return r


def rank_gram(M):
    m, n = M.shape
    G = (M.T @ M) if n <= m else (M @ M.T)
    G = np.array(G.todense(), dtype=np.float64)
    ev = sla.eigvalsh(G, overwrite_a=True, check_finite=False)
    ev = np.abs(ev)
    top = ev.max() if len(ev) else 0.0
    if top == 0:
        return 0, None
    s = np.sort(ev)[::-1] / top
    # largest log-gap in the spectrum below 1e-3 decides; report it
    tol = 1e-10
    r = int((s > tol).sum())
    above = s[r - 1] if r > 0 else None
    below = s[r] if r < len(s) else 0.0
    return r, {"smallest_kept": float(above), "largest_dropped": float(below), "tol": tol}


def rank(M, label, log):
    m, n = M.shape
    if min(m, n) == 0 or M.nnz == 0:
        return 0
    t0 = time.time()
    r_e = rank_modp(M) if min(m, n) <= EXACT_LIMIT else None
    t1 = time.time()
    r_g, gap = rank_gram(M)
    t2 = time.time()
    log.append({"what": label, "shape": [m, n], "rank_modp": r_e, "rank_gram": r_g, "gram_gap": gap,
                "sec_modp": round(t1 - t0, 1), "sec_gram": round(t2 - t1, 1)})
    if r_e is not None and r_e != r_g:
        log[-1]["MISMATCH"] = True
    return r_e if r_e is not None else r_g


# ---------- regions ----------
def box(lo, hi):
    return np.array(list(itertools.product(*[range(l, h + 1) for l, h in zip(lo, hi)])), dtype=np.int64)


def region(name, w):
    """returns (V, V_w) as site arrays"""
    if name == "ball":
        g = box((-7,) * 3, (7,) * 3)
        r2 = (g ** 2).sum(1)
        return g[r2 <= 9], g[r2 <= (3 + w) ** 2]
    if name == "cube":
        return box((0,) * 3, (4,) * 3), box((-w,) * 3, (4 + w,) * 3)
    if name == "shell":
        g = box((-7,) * 3, (7,) * 3)
        r2 = (g ** 2).sum(1)
        return g[(r2 > 4) & (r2 <= 16)], g[(r2 > (2 - w / 2) ** 2) & (r2 <= (4 + w) ** 2)]
    if name in ("ring", "filled"):
        g = box((-10, -10, -4), (10, 10, 4))
        mx = np.abs(g[:, :2]).max(1)
        V = g[(mx >= 5) & (mx <= 6) & (g[:, 2] == 0)]
        inner = (5 - w) if name == "ring" else 0
        Vw = g[(mx >= inner) & (mx <= 6 + w) & (np.abs(g[:, 2]) <= w)]
        return V, Vw
    # --- fix round (post hoc, after the registered ring proved too thin to carry the charges) ---
    # ringA: cross-section 3 x 3 sites; ringB: 4 x 4 sites. "fillA/fillB": same V, hole of V_w filled.
    if name in ("ringA", "fillA", "ringB", "fillB"):
        g = box((-10, -10, -5), (10, 10, 5))
        mx = np.abs(g[:, :2]).max(1)
        z = g[:, 2]
        if name[-1] == "A":
            outer, zlo, zhi = 6, -1, 1
        else:
            outer, zlo, zhi = 7, -1, 2
        V = g[(mx >= 4) & (mx <= outer) & (z >= zlo) & (z <= zhi)]
        inner = (4 - w) if name.startswith("ring") else 0
        Vw = g[(mx >= inner) & (mx <= outer + w) & (z >= zlo - w) & (z <= zhi + w)]
        return V, Vw
    # larger ball and shell, also post hoc, to give the null cells more constrained fields
    if name == "ballL":
        g = box((-8,) * 3, (8,) * 3)
        r2 = (g ** 2).sum(1)
        return g[r2 <= 25], g[r2 <= (5 + w) ** 2]
    if name == "shellL":
        g = box((-8,) * 3, (8,) * 3)
        r2 = (g ** 2).sum(1)
        return g[(r2 > 9) & (r2 <= 36)], g[(r2 > (3 - w / 2) ** 2) & (r2 <= (6 + w) ** 2)]
    raise ValueError(name)


def D_w(cx, V, Vw, log, tag):
    n_pot, n_fld, Cst, Dst = cx
    Vset = set(map(tuple, V))
    # dim Z(V)
    Dm, _ = assemble(Dst, n_fld, max(e for e, _ in Dst) + 1, V)
    dimZ = n_fld * len(V) - rank(Dm, tag + " div|V", log)
    # C restricted to potentials on V_w
    Cm, rows = assemble(Cst, n_pot, n_fld, Vw)
    outside = np.array([tuple(r) not in Vset for r in rows])
    rowmask = np.repeat(outside, n_fld)
    r_all = rank(Cm, tag + " C I_w", log)
    r_out = rank(Cm[rowmask], tag + " R_out C I_w", log)
    dimB = r_all - r_out
    return {"sites_V": int(len(V)), "sites_Vw": int(len(Vw)), "dimZ": int(dimZ), "dimB_w": int(dimB),
            "D_w": int(dimZ - dimB)}


def main():
    out_path = sys.argv[1]
    plan = json.loads(sys.argv[3]) if len(sys.argv) > 3 else None
    default = {"ball": [0, 1, 2, 3], "cube": [0, 1, 2, 3], "shell": [0, 1, 2], "ring": [0, 1, 2, 3],
               "filled": [0, 1, 2, 3]}
    plan = plan or {c: default for c in ("maxwell", "calabi", "cotton")}
    res, log = {}, []
    for cname, regs in plan.items():
        cx = build_complex(cname)
        res[cname] = {}
        for rname, ws in regs.items():
            res[cname][rname] = {}
            for w in ws:
                V, Vw = region(rname, w)
                t0 = time.time()
                r = D_w(cx, V, Vw, log, f"{cname}/{rname}/w{w}")
                r["sec"] = round(time.time() - t0, 1)
                res[cname][rname][str(w)] = r
                print(cname, rname, "w=", w, r, flush=True)
                json.dump({"results": res, "rank_log": log, "prime": P_MOD, "exact_limit": EXACT_LIMIT},
                          open(out_path, "w"), indent=1)
    mism = [l for l in log if l.get("MISMATCH")]
    print("rank method mismatches:", len(mism))


if __name__ == "__main__":
    main()
