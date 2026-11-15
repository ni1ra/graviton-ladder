"""EVEN-REMAINDER: uniform-in-l bound on the second-order remainder of the cocycle for EVEN-PASSAGE's family.

Pass lines: ..\\logs\\01_prefixed_pass_lines.log (fixed before the first run).
Part C: ball arithmetic (python-flint) for E_lo, N3_hi, C0 and (via even_passage.certified) c.
Part N: float64 grid (even-passage's discretisation). Full cocycle S = exp(tau(X0+gD)) exp(-tau X0) by exact
        functional calculus of Omega^2 = h^2 + g h W h (cross-checked against a Taylor matrix exponential).
Usage: ..\\..\\even-nonzero\\venv\\Scripts\\python.exe -B even_remainder.py ..\\out\\even_remainder.json
"""
import json
import math
import os
import sys
import time

sys.dont_write_bytecode = True          # nothing may be written into even-passage\
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "even-passage", "scripts"))

import numpy as np
from flint import arb, ctx

import even_passage as ep

ctx.prec = 100
TAU, EPS = ep.TAU, ep.EPS
LS = [16, 64, 256, 1024, 16384]
GS = [1.0, 8.0]
SQ2 = math.sqrt(2.0)
lam_of, lo, hi = ep.lam_of, ep.lo, ep.hi


# ------------------------------------------------------------------ Part C
def certified():
    PI = arb.pi()
    out, unif, _, _ = ep.certified()
    c = lo(unif["c"])
    A_inf, G_hi = out["A_inf"], out["G_hi"]
    # (c) ||e_l||^2 >= (1/2pi) <hu, e^(-2 eps h^2) hu> >= (1/2pi) A_l exp(-2 eps G_l^2/A_l), monotone in (A up, G down)
    E_lo = (A_inf * (-(2 * EPS * G_hi * G_hi / A_inf)).exp() / (2 * PI)).sqrt()
    # (b) sup_p (p coth(pi p))^3 e^(-2 eps p^2).  p coth(pi p) is increasing, e^(-2 eps p^2) decreasing:
    #     on a cell [a, b] the product is <= (b coth(pi b))^3 e^(-2 eps a^2).  Tail p >= 6: <= (6 + 1/pi)^3 e^(-108).
    P, n = 6.0, 6144
    best = (6 + 1 / math.pi) ** 3 * math.exp(-108.0) * 1.01
    arg = 0.0
    for i in range(n):
        a_, b_ = arb(i) * arb(P) / n, arb(i + 1) * arb(P) / n
        cothb = (PI * b_).cosh() / (PI * b_).sinh()
        v = hi((b_ * cothb) ** 3 * (-2 * EPS * a_ * a_).exp())
        if v > best:
            best, arg = v, (i + 0.5) * P / n
    sup3 = arb(best)
    u2 = ep.integ(lambda y: ep.profile(y)[0] ** 2)
    u2_exact = arb(4) * arb(65536) / arb(109395)                 # 4 * 2^17 (8!)^2 / 17!
    assert (u2 - u2_exact).abs_upper() < 1e-20
    N3_hi = (sup3 * u2 / 2).sqrt()
    Rb = N3_hi / E_lo
    C0 = (1 + Rb) * PI * PI / 2
    return dict(c=c, E_lo=lo(E_lo), N3_hi=hi(N3_hi), Rb=hi(Rb), C0=hi(C0), sup3=best, sup3_argmax=arg,
                u2=hi(u2), A_inf_lo=lo(A_inf), G_hi=hi(G_hi), E_hi=hi(unif["E_hi"]))


def em1x(x):
    """Upper bound for e^x - 1 - x (x >= 0), float with outward fudge."""
    return (math.expm1(x) - x) * (1 + 1e-12) + 1e-300 if x > 1e-4 else 0.5 * x * x * math.exp(x) * (1 + 1e-12)


# ------------------------------------------------------------------ Part N
def setup(l, dy=0.02, ytop=40.0, ybot=None):
    lam = lam_of(l)
    y = ep.grid(dy, ytop, max(ep.YBOT if ybot is None else ybot, -math.log(2 * lam)))
    pot, Wmat = lam ** 2 / np.sinh(y + math.log(2 * lam)) ** 2, ep.kern_l(l, y) * dy      # Wmat = lam * W_l
    n = len(y)
    Hm = np.zeros((n, n))
    i = np.arange(n)
    Hm[i, i] = 2.0 / dy ** 2 + pot
    Hm[i[:-1], i[:-1] + 1] = -1.0 / dy ** 2
    Hm[i[:-1] + 1, i[:-1]] = -1.0 / dy ** 2
    Hm[-1, -1] -= 1.0 / dy ** 2
    p2, U = np.linalg.eigh(Hm)
    assert p2.min() > 0
    p = np.sqrt(p2)
    s_ = (y - ep.Y0) / ep.B
    u = np.where(np.abs(s_) < 1, (1 - s_ * s_) ** ep.K, 0.0)
    cu = U.T @ u
    a = p2 * np.exp(-EPS * p2) * cu
    W = (U.T @ Wmat @ U) / lam                                    # true W_l in the eigenbasis of h
    W = 0.5 * (W + W.T)
    c2 = 1.0 / (-np.expm1(-2 * math.pi * p))
    s2 = np.exp(-2 * math.pi * p) * c2
    sq = np.sqrt(p)
    M = sq[:, None] * W * sq[None, :]                             # h^(1/2) W h^(1/2)

    def F(theta):
        th = TAU * theta
        return TAU * np.exp(-0.5j * th) * np.sinc(th / (2 * math.pi))

    A = -0.5j * M * F(p[:, None] - p[None, :])                    # T z = A z + B zbar, T = int_0^tau e^(sX0) D_1 e^(-sX0) ds
    Bm = 0.5j * M * F(p[:, None] + p[None, :])
    return dict(l=l, lam=lam, dy=dy, n=n, p=p, a=a, W=W, M=M, A=A, B=Bm, c2=c2, s2=s2, coth=c2 + s2)


def pert(st, g):
    p = st["p"]
    Om2 = g * (p[:, None] * st["W"] * p[None, :])
    Om2[np.diag_indices_from(Om2)] += p * p
    w2, V = np.linalg.eigh(Om2)
    assert w2.min() > 0
    return V, np.sqrt(w2)


def fwd(st, V, w, X, Y):
    """exp(tau(X0 + g D_1)) on real coordinates x = h^(1/2) q, y = h^(-1/2) pi (columns)."""
    sp = np.sqrt(st["p"])[:, None]
    xi, eta = V.T @ (sp * X), V.T @ (Y / sp)
    cw, sw, ww = np.cos(TAU * w)[:, None], np.sin(TAU * w)[:, None], w[:, None]
    return (V @ (cw * xi + ww * sw * eta)) / sp, sp * (V @ (cw * eta - sw / ww * xi))


def back0(st, X, Y):
    """exp(-tau X0)."""
    c0, s0 = np.cos(TAU * st["p"])[:, None], np.sin(TAU * st["p"])[:, None]
    return c0 * X - s0 * Y, s0 * X + c0 * Y


def S_apply(st, V, w, X, Y):
    return fwd(st, V, w, *back0(st, X, Y))


def T_apply(st, X, Y):
    Z = (X + 1j * Y) / SQ2
    t = st["A"] @ Z + st["B"] @ np.conj(Z)
    return SQ2 * t.real, SQ2 * t.imag


def beta_norm(st, w1, w2):
    """|| (1/2)[k(w1) + i k(w2)] ||_K = || (cP, s Qbar) ||, P = (w1 + i w2)/2, |Q| = |w1 - i w2|/2."""
    P, Q = 0.5 * (w1 + 1j * w2), 0.5 * (w1 - 1j * w2)
    return math.sqrt(st["dy"] * float(np.sum(st["c2"] * np.abs(P) ** 2 + st["s2"] * np.abs(Q) ** 2)))


def vectors(st, gs):
    p, a, coth, dy = st["p"], st["a"], st["coth"], st["dy"]
    yv = a / np.sqrt(p)                                           # psi = i a/sqrt(2p):  x = 0, y = a/sqrt(p)
    n = st["n"]
    X = np.stack([np.zeros(n), -coth * yv], axis=1)               # columns: psi, i coth psi
    Y = np.stack([yv, np.zeros(n)], axis=1)
    res = {"e": math.sqrt(dy * float(np.sum(coth * yv * yv)) / 2),
           "n3": math.sqrt(dy * float(np.sum(coth ** 3 * yv * yv)) / 2)}
    tx, ty = T_apply(st, X, Y)
    tw = (tx + 1j * ty) / SQ2
    res["beta1"] = beta_norm(st, tw[:, 0], tw[:, 1])              # ||beta1hat e||, unscaled
    for g in gs:
        V, w = pert(st, g)
        sx, sy = S_apply(st, V, w, X, Y)
        sw = ((sx - X) + 1j * (sy - Y)) / SQ2
        rw = sw - g * tw
        res[g] = {"beta": beta_norm(st, sw[:, 0], sw[:, 1]), "rem": beta_norm(st, rw[:, 0], rw[:, 1])}
    return res


def hnorm(st, Mat):
    cd = np.sqrt(np.concatenate([st["coth"], st["coth"]]))
    return float(np.linalg.norm(cd[:, None] * Mat / cd[None, :], 2))


def opnorms(st, gs):
    n = st["n"]
    I, Z = np.eye(n), np.zeros((n, n))
    X, Y = np.hstack([I, Z]), np.hstack([Z, I])
    res = {"D": hnorm(st, np.block([[Z, st["M"]], [Z, Z]]))}
    Tm = np.vstack(T_apply(st, X, Y))
    for g in gs:
        V, w = pert(st, g)
        Sm = np.vstack(S_apply(st, V, w, X, Y))
        res[g] = {"S": hnorm(st, Sm), "R2": hnorm(st, Sm - np.eye(2 * n) - g * Tm)}
    return res


def expm_taylor(G):
    k = max(0, int(math.ceil(math.log2(np.linalg.norm(G, 1)))) + 3)
    Gs = G / 2.0 ** k
    E, term = np.eye(len(G)), np.eye(len(G))
    for j in range(1, 25):
        term = term @ Gs / j
        E = E + term
    for _ in range(k):
        E = E @ E
    return E


def expm_check(l, g):
    st = setup(l, dy=0.2)
    n, p = st["n"], st["p"]
    I, Z = np.eye(n), np.zeros((n, n))
    V, w = pert(st, g)
    Efc = np.vstack(fwd(st, V, w, np.hstack([I, Z]), np.hstack([Z, I])))
    G = np.block([[Z, np.diag(p) + g * st["M"]], [-np.diag(p), Z]])
    Eex = expm_taylor(TAU * G)
    return float(np.linalg.norm(Efc - Eex) / np.linalg.norm(Eex))


def main():
    t0 = time.time()
    C = certified()
    c, E_lo, N3_hi, Rb, C0 = C["c"], C["E_lo"], C["N3_hi"], C["Rb"], C["C0"]
    print("== Part C: certified constants (ball arithmetic, %d bits; values are outward-rounded floats) ==" % ctx.prec)
    for k, v in C.items():
        print(f"  {k:12s} {v!r}")
    l1 = lambda g: 2 * C0 * 1.001 * abs(g) / c
    print(f"  l_1(1) = 2*C0*1.001/c = {l1(1):.1f}   (C = 1.001*C0 = {1.001 * C0:.4f}, valid while sqrt(2) pi g/lam <= 9.9e-4)")
    lines = {}
    lines["C1"] = (all(math.isfinite(v) for v in C.values()) and 0.145 <= E_lo <= 0.155 and 0.30 <= N3_hi <= 0.34
                   and 14 <= C0 <= 17)
    lines["C2"] = abs(c / 2.110e-3 - 1) < 5e-4

    rem_hi = lambda lam, g: 0.5 * (1 + Rb) * em1x(TAU * SQ2 * abs(g) / lam)
    low = lambda lam, g: c * abs(g) / lam - rem_hi(lam, g)

    print("== Part N: fine grid dy = 0.02 ==", flush=True)
    R, ST = {}, {}
    for l in LS:
        st = setup(l)
        R[l] = vectors(st, GS)
        print(f"  l={l:6d} n={st['n']} ||e||={R[l]['e']:.6f} ||coth^1.5 psi||={R[l]['n3']:.6f} "
              f"lam*b1/e={st['lam'] * R[l]['beta1'] / R[l]['e']:.8f}", flush=True)
        del st
    OO, RCO, O, RC = {}, {}, {}, {}
    for title, yb, Od, Rd in (("as pre-fixed: box [-4.8, 40] (kernel unresolved at the bottom, see logs 03, 04)", None, OO, RCO),
                              ("box [-3, 40] (added after run 1)", -3.0, O, RC)):
        print(f"== coarse grid dy = 0.04, operator norms on H, {title} ==", flush=True)
        for l in LS:
            st = setup(l, dy=0.04, ybot=yb)
            Od[l] = opnorms(st, GS)
            Rd[l] = vectors(st, GS)
            print(f"  l={l:6d} n={st['n']} lam||D||={st['lam'] * Od[l]['D']:.5f} " + " ".join(
                f"g={g:g}: ||S||={Od[l][g]['S']:.8f} lam^2||R2||={st['lam'] ** 2 * Od[l][g]['R2']:.5f}" for g in GS), flush=True)
            del st

    P = {k: [] for k in ("P1", "P2", "P3", "P4", "P5", "P6", "P5b", "P6b")}
    rows = []
    for l in LS:
        lam, r = lam_of(l), R[l]
        P["P1"].append(E_lo <= r["e"] <= C["E_hi"])
        P["P2"].append(r["n3"] <= N3_hi)
        P["P5"].append(OO[l]["D"] <= SQ2 / lam)
        P["P5b"].append(O[l]["D"] <= SQ2 / lam)
        for g in GS:
            x = TAU * SQ2 * g / lam
            act, first, rem = r[g]["beta"] / r["e"], g * r["beta1"] / r["e"], r[g]["rem"] / r["e"]
            rh, lw = rem_hi(lam, g), low(lam, g)
            rows.append(dict(l=l, g=g, lam=lam, actual=act, first=first, rem=rem, rem_hi=rh, low=lw, half=c * g / (2 * lam),
                             lam2rem=lam * lam * rem, opR2=O[l][g]["R2"], opS=O[l][g]["S"]))
            P["P3"].append(rem <= rh)
            P["P4"].append(act >= lw)
            for kO, kR, Od, Rd in (("P5", "P6", OO, RCO), ("P5b", "P6b", O, RC)):
                P[kO].append(Od[l][g]["S"] <= math.exp(x) * (1 + 1e-12) and Od[l][g]["R2"] <= em1x(x))
                P[kR].append(Rd[l][g]["rem"] <= 0.5 * Od[l][g]["R2"] * (Rd[l]["e"] + Rd[l]["n3"]) * (1 + 1e-9))
    for k, v in P.items():
        lines[k] = all(v)
        if not all(v):
            print(f"  {k} detail: {v}")
    row = {(r["l"], r["g"]): r for r in rows}
    lines["N1"] = (all(40 <= row[(l, 8.0)]["rem"] / row[(l, 1.0)]["rem"] <= 90 for l in LS)
                   and abs(row[(1024, 1.0)]["lam2rem"] / row[(256, 1.0)]["lam2rem"] - 1) <= 0.15)
    ref = {16: 0.121788, 64: 0.123952, 256: 0.124501, 1024: 0.124639}
    n2 = {l: abs(lam_of(l) * R[l]["beta1"] / R[l]["e"] / ref[l] - 1) for l in ref}
    lines["N2"] = max(n2.values()) < 1e-5
    n3 = {f"l{l}_g{g:g}": expm_check(l, g) for l in (16, 64) for g in GS}
    lines["N3"] = max(n3.values()) < 1e-8
    slack = {l: row[(l, 1.0)]["rem_hi"] / row[(l, 1.0)]["rem"] for l in LS}
    lines["N4"] = all(3 <= s <= 300 for s in slack.values())
    r16k = row[(16384, 1.0)]
    lines["N5"] = (r16k["low"] > r16k["half"] and r16k["actual"] >= r16k["low"]
                   and all(row[(l, 1.0)]["low"] < 0 for l in (16, 64, 256, 1024)))
    lines["N6"] = all(row[(l, 8.0)]["low"] < 0 and row[(l, 8.0)]["actual"] >= row[(l, 8.0)]["half"] for l in LS)
    stf = setup(1024, dy=0.01)
    rf = vectors(stf, [1.0])
    del stf
    n7 = {"lam2rem": abs((rf[1.0]["rem"] / rf["e"]) / row[(1024, 1.0)]["rem"] - 1),
          "ratio": abs((rf[1.0]["beta"] / rf["e"]) / row[(1024, 1.0)]["actual"] - 1)}
    lines["N7"] = n7["lam2rem"] < 0.02 and n7["ratio"] < 1e-3
    n8 = {l: abs((RC[l][1.0]["rem"] / RC[l]["e"]) / row[(l, 1.0)]["rem"] - 1) for l in LS}
    lines["N8"] = max(n8.values()) < 0.05

    # plants
    pA = {l: C0 * math.exp(TAU * SQ2 / lam_of(l)) / lam_of(l) ** 3 for l in LS}
    pA_viol = {l: row[(l, 1.0)]["rem"] > pA[l] for l in LS}
    lines["PLANT_A_caught"] = pA_viol[1024] and pA_viol[16384]
    pB = {(l, g): em1x(0.5 * SQ2 * g / lam_of(l)) for l in LS for g in GS}
    pB_viol = {k: O[k[0]][k[1]]["R2"] > v for k, v in pB.items()}
    lines["PLANT_B_caught"] = any(pB_viol.values())

    print("== table ==  (c = %.4e, C0 = %.4f, Rb = %.4f, E_lo = %.5f, N3_hi = %.5f)" % (c, C0, Rb, E_lo, N3_hi))
    print("      l    g   actual ||b e||/||e||   g||b1 e||/||e||   actual rem     proved rem_hi   proved low      c g/(2 lam)   "
          "lam^2*rem   rem_hi/rem  lam^2||R2||_H  bound lam^2(e^x-1-x)")
    for r in rows:
        x = TAU * SQ2 * r["g"] / r["lam"]
        print(f"  {r['l']:5d}  {r['g']:3g}   {r['actual']:.9e}       {r['first']:.9e}   {r['rem']:.6e}   {r['rem_hi']:.6e}   "
              f"{r['low']:+.6e}   {r['half']:.6e}   {r['lam2rem']:.5f}    {r['rem_hi'] / r['rem']:8.2f}   "
              f"{r['lam'] ** 2 * r['opR2']:.5f}        {r['lam'] ** 2 * em1x(x):.5f}")
    print("  grid norms: ||e_l|| " + ", ".join(f"{R[l]['e']:.5f}" for l in LS) + f"  (proved range [{E_lo:.5f}, {C['E_hi']:.5f}])")
    print("              ||coth^1.5 psi_l|| " + ", ".join(f"{R[l]['n3']:.5f}" for l in LS) + f"  (proved <= {N3_hi:.5f})")
    print("  lam||D_1||_H (box [-3, 40]): " + ", ".join(f"{lam_of(l) * O[l]['D']:.5f}" for l in LS) + "  (proved <= 1.41421)")
    print(f"  N2 rel. diff to EVEN-PASSAGE table: {n2}")
    print(f"  N3 expm cross-check: {n3}")
    print(f"  N7 grid stability at l=1024: {n7}")
    print(f"  N8 box [-3, 40], dy 0.04 vs fine grid, rel. diff of remainder: {n8}")
    print("  lam||D_1||_H on the pre-fixed box [-4.8, 40]: " + ", ".join(f"{lam_of(l) * OO[l]['D']:.5f}" for l in LS))
    print(f"  PLANT A (g^2/lam^3): planted {[f'{pA[l]:.3e}' for l in LS]} vs rem {[f'{row[(l, 1.0)]['rem']:.3e}' for l in LS]}; "
          f"violated {list(pA_viol.values())}; caught = {lines['PLANT_A_caught']}")
    print(f"  PLANT B (tau = 0.5): planted lam^2 bound {[f'{lam_of(l) ** 2 * pB[(l, 1.0)]:.4f}' for l in LS]} (g=1) vs lam^2||R2|| "
          f"{[f'{lam_of(l) ** 2 * O[l][1.0]['R2']:.4f}' for l in LS]}; violated {pB_viol}; caught = {lines['PLANT_B_caught']}")
    print("== pass lines ==")
    for k, v in lines.items():
        print(f"  {k}: {'PASS' if v else 'FAIL'}")
    ok = all(lines.values())
    print("ALL LINES HELD" if ok else "SOME LINES FAILED", f"({time.time() - t0:.0f} s)")
    if len(sys.argv) > 1:
        with open(sys.argv[1], "w") as f:
            json.dump({"certified": C, "l1_g1": l1(1), "rows": rows, "opnorms": {str(l): {str(k): v for k, v in O[l].items()} for l in LS},
                       "grid": {str(l): {str(k): v for k, v in R[l].items()} for l in LS}, "n2": {str(k): v for k, v in n2.items()},
                       "n3": n3, "n7": n7, "n8": {str(k): v for k, v in n8.items()},
                       "opnorms_prefixed_box": {str(l): {str(k): v for k, v in OO[l].items()} for l in LS}, "plantA": {str(l): pA[l] for l in LS},
                       "plantB_violated": {f"{k[0]}_{k[1]:g}": bool(v) for k, v in pB_viol.items()},
                       "lines": {k: bool(v) for k, v in lines.items()}}, f, indent=1)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
