"""ODD-REMAINDER: uniform-in-l bound on the second-order remainder of the cocycle for ODD-SIGN's family.

D0(q, pi) = (0, -N q), N = 1/(1 + cosh rho).  e_l = (q_l, 0), q_l = exp(-eps h^2) h^2 u  (vector in the q slot).
Pass lines: ..\\logs\\01_prefixed_pass_lines.log (fixed before the first run).
Part C: ball arithmetic (python-flint) for E_lo, N3_hi, K, C0 and (via odd_sign.certified) c.
Part N: float64 grid.  Full cocycle S = exp(tau(X0+gD0)) exp(-tau X0) by exact functional calculus of
        Omega^2 = h^2 + g N  (q'' = -(h^2 + g N) q), cross-checked against a Taylor matrix exponential.
Box: y in [max(-3, -ln 2 lam), 40]; dy = 0.02 (vectors), dy = 0.04 (operator norms on H).
Usage: ..\\..\\even-nonzero\\venv\\Scripts\\python.exe -B odd_remainder.py ..\\out\\odd_remainder.json
"""
import json
import math
import os
import sys
import time

sys.dont_write_bytecode = True          # nothing may be written into odd-sign\
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "odd-sign", "scripts"))

import numpy as np
from flint import arb, ctx

import odd_sign as osn

ctx.prec = 100
TAU, EPS = osn.TAU, osn.EPS
LS = [16, 64, 256, 1024, 16384]
GS = [1.0, 8.0]
SQ2 = math.sqrt(2.0)
YBOT = -3.0
lam_of, lo, hi = osn.lam_of, osn.lo, osn.hi


# ------------------------------------------------------------------ Part C
def sup_mono(power, P=6.0, n=6144):
    """Upper bound for sup_{p>=0} p^power (p coth(pi p))^3 e^(-2 eps p^2), power in {0, 2}.
    p, p coth(pi p) increasing and e^(-2 eps p^2) decreasing: on [a, b] the product is <= b^power (b coth pi b)^3 e^(-2 eps a^2).
    Tail p >= 6: p coth(pi p) <= p + 1/pi, so f <= p^power (p + 1/pi)^3 e^(-3 p^2), decreasing there; value at 6 times 1.01."""
    PI = arb.pi()
    best = 6.0 ** power * (6 + 1 / math.pi) ** 3 * math.exp(-108.0) * 1.01
    arg = 0.0
    for i in range(n):
        a_, b_ = arb(i) * arb(P) / n, arb(i + 1) * arb(P) / n
        cothb = (PI * b_).cosh() / (PI * b_).sinh()
        v = hi(b_ ** power * (b_ * cothb) ** 3 * (-2 * EPS * a_ * a_).exp())
        if v > best:
            best, arg = v, (i + 0.5) * P / n
    return best, arg


def certified():
    PI = arb.pi()
    out, unif, _, _ = osn.certified()
    c = lo(unif["c"])
    A_16, J0_lo = out["A_16"], unif["J0_lo"]
    u2 = osn.integ(lambda y: osn.profile(y)[0] ** 2)
    u2_exact = arb(4) * arb(65536) / arb(109395)                 # B * 2^17 (8!)^2 / 17!
    assert (u2 - u2_exact).abs_upper() < 1e-20
    # E_lo: ||e||^2 = (1/2) int p^3 coth e^(-2 eps p^2) dmu >= (1/2pi) int p^2 e^(-2 eps p^2) dmu >= (1/2pi) J0^2/||u||^2
    E_lo = J0_lo / (2 * PI * u2_exact).sqrt()
    s3, a3 = sup_mono(0)
    s5, a5 = sup_mono(2)
    N3a = (arb(s3) * A_16 / 2).sqrt()                            # measure of hu, mass A_l <= A_16
    N3b = (arb(s5) * u2_exact / 2).sqrt()                        # measure of u, mass ||u||^2
    N3_hi = N3a if hi(N3a) <= hi(N3b) else N3b
    Rb = N3_hi / E_lo
    Kc = (1 + Rb) / 2
    C0 = Kc * PI * PI
    return dict(c=c, E_lo=lo(E_lo), N3_hi=hi(N3_hi), N3a=hi(N3a), N3b=hi(N3b), Rb=hi(Rb), K=hi(Kc), C0=hi(C0),
                s3=s3, s3_argmax=a3, s5=s5, s5_argmax=a5, u2=hi(u2), A_16_hi=hi(A_16), J0_lo=lo(J0_lo),
                E_hi=hi(unif["E_hi"]))


def em1x(x):
    """Upper bound for e^x - 1 - x (x >= 0), float with outward fudge."""
    return (math.expm1(x) - x) * (1 + 1e-12) + 1e-300 if x > 1e-4 else 0.5 * x * x * math.exp(x) * (1 + 1e-12)


# ------------------------------------------------------------------ Part N
def setup(l, dy=0.02, ytop=40.0, ybot=YBOT):
    lam = lam_of(l)
    y = osn.grid(dy, ytop, max(ybot, -math.log(2 * lam)))
    rho = y + math.log(2 * lam)
    pot, lapse = lam ** 2 / np.sinh(rho) ** 2, 1.0 / (1.0 + np.cosh(rho))        # true N (not lam N)
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
    s_ = (y - osn.Y0) / osn.B
    u = np.where(np.abs(s_) < 1, (1 - s_ * s_) ** osn.K, 0.0)
    cu = U.T @ u
    a = p2 * np.exp(-EPS * p2) * cu                               # q_l in the eigenbasis of h
    Nt = (U.T * lapse[None, :]) @ U                               # N in the eigenbasis
    Nt = 0.5 * (Nt + Nt.T)
    c2 = 1.0 / (-np.expm1(-2 * math.pi * p))
    s2 = np.exp(-2 * math.pi * p) * c2
    isq = 1.0 / np.sqrt(p)
    M = isq[:, None] * Nt * isq[None, :]                          # h^(-1/2) N h^(-1/2);  D0: (x, y) -> (0, -M x)

    def F(theta):                                                 # int_0^tau e^(-i s theta) ds
        th = TAU * theta
        return TAU * np.exp(-0.5j * th) * np.sinc(th / (2 * math.pi))

    # D0 z = -(i/2) M (z + zbar);  T z = A z + B zbar,  T = int_0^tau e^(sX0) D0 e^(-sX0) ds,  e^(sX0) z = e^(-ish) z
    A = -0.5j * M * F(p[:, None] - p[None, :])
    Bm = -0.5j * M * F(p[:, None] + p[None, :])
    return dict(l=l, lam=lam, dy=dy, n=n, p=p, a=a, Nt=Nt, M=M, A=A, B=Bm, c2=c2, s2=s2, coth=c2 + s2)


def pert(st, g):
    Om2 = g * st["Nt"]
    Om2[np.diag_indices_from(Om2)] += st["p"] ** 2
    w2, V = np.linalg.eigh(Om2)
    assert w2.min() > 0
    return V, np.sqrt(w2)


def fwd(st, V, w, X, Y):
    """exp(tau(X0 + g D0)) on real coordinates x = h^(1/2) q, y = h^(-1/2) pi (columns): q' = pi, pi' = -(h^2 + gN) q."""
    sp = np.sqrt(st["p"])[:, None]
    xi, eta = V.T @ (X / sp), V.T @ (sp * Y)
    cw, sw, ww = np.cos(TAU * w)[:, None], np.sin(TAU * w)[:, None], w[:, None]
    return sp * (V @ (cw * xi + sw / ww * eta)), (V @ (cw * eta - ww * sw * xi)) / sp


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


def vectors(st, gs, plant=False):
    p, a, coth, dy = st["p"], st["a"], st["coth"], st["dy"]
    xv = np.sqrt(p) * a                                           # psi = h^(1/2) q/sqrt 2 (real):  x = sqrt(p) a, y = 0
    n = st["n"]
    X = np.stack([xv, np.zeros(n)], axis=1)                       # columns: psi, i coth psi
    Y = np.stack([np.zeros(n), coth * xv], axis=1)
    res = {"e": math.sqrt(dy * float(np.sum(coth * xv * xv)) / 2),
           "n3": math.sqrt(dy * float(np.sum(coth ** 3 * xv * xv)) / 2)}
    tx, ty = T_apply(st, X, Y)
    tw = (tx + 1j * ty) / SQ2
    res["beta1"] = beta_norm(st, tw[:, 0], tw[:, 1])              # ||beta1hat e||, unscaled
    for g in gs:
        V, w = pert(st, g)
        sx, sy = S_apply(st, V, w, X, Y)
        sw = ((sx - X) + 1j * (sy - Y)) / SQ2
        rw = sw - g * tw
        res[g] = {"beta": beta_norm(st, sw[:, 0], sw[:, 1]), "rem": beta_norm(st, rw[:, 0], rw[:, 1])}
        if plant:                                                 # PLANT A: cocycle factor exp(-tau X0) dropped
            px, py = fwd(st, V, w, X, Y)
            pw = ((px - X) + 1j * (py - Y)) / SQ2 - g * tw
            res[g]["plant_rem"] = beta_norm(st, pw[:, 0], pw[:, 1])
    return res


def hnorm(st, Mat):
    cd = np.sqrt(np.concatenate([st["coth"], st["coth"]]))
    return float(np.linalg.norm(cd[:, None] * Mat / cd[None, :], 2))


def opnorms(st, gs):
    n = st["n"]
    I, Z = np.eye(n), np.zeros((n, n))
    X, Y = np.hstack([I, Z]), np.hstack([Z, I])
    res = {"D": hnorm(st, np.block([[Z, Z], [-st["M"], Z]]))}
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
    G = np.block([[Z, np.diag(p)], [-np.diag(p) - g * st["M"], Z]])       # x' = h y, y' = -(h + g M) x
    Eex = expm_taylor(TAU * G)
    return float(np.linalg.norm(Efc - Eex) / np.linalg.norm(Eex))


def main():
    t0 = time.time()
    C = certified()
    c, E_lo, N3_hi, Rb, Kc, C0 = C["c"], C["E_lo"], C["N3_hi"], C["Rb"], C["K"], C["C0"]
    print("== Part C: certified constants (ball arithmetic, %d bits; values are outward-rounded floats) ==" % ctx.prec)
    for k, v in C.items():
        print(f"  {k:12s} {v!r}")
    l1 = lambda g: max(16, math.ceil(2 * C0 * 1.001 * abs(g) / c))
    print(f"  K = (1/2)(1 + N3_hi/E_lo) = {Kc:.6f};  C0 = K pi^2 = {C0:.6f}")
    print(f"  positivity threshold |g|/lam < c/(1.001 C0) = {c / (1.001 * C0):.6e};  half-floor threshold c/(2*1.001*C0) = {c / (2.002 * C0):.6e}")
    print(f"  l_1(g) = max(16, ceil({2 * C0 * 1.001 / c:.2f} |g|));  l_1(1) = {l1(1)}, l_1(8) = {l1(8)}   (1.001 covers e^x for x <= 9.9e-4)")
    lines = {}
    lines["C1"] = (all(math.isfinite(v) for v in C.values()) and 0.090 <= E_lo <= 0.097 and 0.19 <= N3_hi <= 0.21
                   and 1.50 <= Kc <= 1.65 and 14.8 <= C0 <= 16.3)
    lines["C2"] = abs(c / 2.038e-3 - 1) < 5e-4

    rem_hi = lambda lam, g: Kc * em1x(TAU * SQ2 * abs(g) / lam)
    low = lambda lam, g: c * abs(g) / lam - rem_hi(lam, g)

    print(f"== Part N: fine grid dy = 0.02, box [{YBOT:g}, 40] ==", flush=True)
    R = {}
    for l in LS:
        st = setup(l)
        R[l] = vectors(st, GS, plant=True)
        print(f"  l={l:6d} n={st['n']} ||e||={R[l]['e']:.6f} ||coth^1.5 psi||={R[l]['n3']:.6f} "
              f"lam*b1/e={st['lam'] * R[l]['beta1'] / R[l]['e']:.8f}", flush=True)
        del st
    O, RC = {}, {}
    print(f"== coarse grid dy = 0.04, box [{YBOT:g}, 40]: operator norms on H ==", flush=True)
    for l in LS:
        st = setup(l, dy=0.04)
        O[l] = opnorms(st, GS)
        RC[l] = vectors(st, GS)
        print(f"  l={l:6d} n={st['n']} lam||D0||={st['lam'] * O[l]['D']:.5f} " + " ".join(
            f"g={g:g}: ||S||={O[l][g]['S']:.8f} lam^2||R2||={st['lam'] ** 2 * O[l][g]['R2']:.5f}" for g in GS), flush=True)
        del st
    info48 = {}
    for l in LS:                                                  # information only: the box EVEN-REMAINDER found unresolved
        st = setup(l, dy=0.04, ybot=-4.8)
        n = st["n"]
        Z = np.zeros((n, n))
        info48[l] = st["lam"] * hnorm(st, np.block([[Z, Z], [-st["M"], Z]]))
        del st

    P = {k: [] for k in ("P1", "P2", "P3", "P4", "P5", "P6")}
    rows = []
    for l in LS:
        lam, r = lam_of(l), R[l]
        P["P1"].append(E_lo <= r["e"] <= C["E_hi"])
        P["P2"].append(r["n3"] <= N3_hi)
        P["P5"].append(O[l]["D"] <= SQ2 / lam)
        for g in GS:
            x = TAU * SQ2 * g / lam
            act, first, rem = r[g]["beta"] / r["e"], g * r["beta1"] / r["e"], r[g]["rem"] / r["e"]
            rh, lw = rem_hi(lam, g), low(lam, g)
            rows.append(dict(l=l, g=g, lam=lam, actual=act, first=first, rem=rem, rem_hi=rh, low=lw, half=c * g / (2 * lam),
                             lam2rem=lam * lam * rem, opR2=O[l][g]["R2"], opS=O[l][g]["S"],
                             plant_rem=r[g]["plant_rem"] / r["e"]))
            P["P3"].append(rem <= rh)
            P["P4"].append(act >= lw)
            P["P5"].append(O[l][g]["S"] <= math.exp(x) * (1 + 1e-12) and O[l][g]["R2"] <= em1x(x))
            P["P6"].append(RC[l][g]["rem"] <= 0.5 * O[l][g]["R2"] * (RC[l]["e"] + RC[l]["n3"]) * (1 + 1e-9))
    for k, v in P.items():
        lines[k] = all(v)
        if not all(v):
            print(f"  {k} detail: {v}")
    row = {(r["l"], r["g"]): r for r in rows}
    g2 = {l: row[(l, 8.0)]["rem"] / row[(l, 1.0)]["rem"] for l in LS}
    lines["N1"] = (all((25 if l == 16 else 40) <= g2[l] <= (110 if l == 16 else 90) for l in LS)
                   and abs(row[(1024, 1.0)]["lam2rem"] / row[(256, 1.0)]["lam2rem"] - 1) <= 0.15)
    ref = {16: 0.297406, 64: 0.303293, 256: 0.304791, 1024: 0.305167}
    n2 = {l: abs(lam_of(l) * R[l]["beta1"] / R[l]["e"] / ref[l] - 1) for l in ref}
    lines["N2"] = max(n2.values()) < 1e-4
    n3 = {f"l{l}_g{g:g}": expm_check(l, g) for l in (16, 64) for g in GS}
    lines["N3"] = max(n3.values()) < 1e-8
    slack = {l: row[(l, 1.0)]["rem_hi"] / row[(l, 1.0)]["rem"] for l in LS}
    lines["N4"] = all(3 <= s <= 2000 for s in slack.values())
    r16k = row[(16384, 1.0)]
    lines["N5"] = (r16k["low"] > r16k["half"] and r16k["actual"] >= r16k["low"]
                   and all(row[(l, 1.0)]["low"] < 0 for l in (16, 64, 256, 1024)))
    lines["N6"] = all(row[(l, 8.0)]["low"] < 0 and row[(l, 8.0)]["actual"] >= row[(l, 8.0)]["half"] for l in LS)
    stf = setup(1024, dy=0.01)
    rf = vectors(stf, [1.0])
    del stf
    n7 = {"rem": abs((rf[1.0]["rem"] / rf["e"]) / row[(1024, 1.0)]["rem"] - 1),
          "ratio": abs((rf[1.0]["beta"] / rf["e"]) / row[(1024, 1.0)]["actual"] - 1)}
    lines["N7"] = n7["rem"] < 0.02 and n7["ratio"] < 1e-3

    # plants
    pA_viol = {l: row[(l, 1.0)]["plant_rem"] > row[(l, 1.0)]["rem_hi"] for l in LS}
    lines["PLANT_A_caught"] = all(pA_viol.values())
    pB = {l: C0 * math.exp(TAU * SQ2 / lam_of(l)) / lam_of(l) ** 3 for l in LS}
    pB_viol = {l: row[(l, 1.0)]["rem"] > pB[l] for l in LS}
    lines["PLANT_B_caught"] = pB_viol[16384]

    print("== table ==  (c = %.4e, K = %.4f, C0 = %.4f, Rb = %.4f, E_lo = %.5f, N3_hi = %.5f)" % (c, Kc, C0, Rb, E_lo, N3_hi))
    print("      l    g   actual ||b e||/||e||   g||b1 e||/||e||   actual rem     proved rem_hi   proved low      c g/(2 lam)   "
          "lam^2*rem   rem_hi/rem  lam^2||R2||_H  bound lam^2(e^x-1-x)  bound")
    for r in rows:
        x = TAU * SQ2 * r["g"] / r["lam"]
        tag = "NON-VACUOUS" if r["low"] > 0 else "vacuous"
        print(f"  {r['l']:5d}  {r['g']:3g}   {r['actual']:.9e}       {r['first']:.9e}   {r['rem']:.6e}   {r['rem_hi']:.6e}   "
              f"{r['low']:+.6e}   {r['half']:.6e}   {r['lam2rem']:.5f}    {r['rem_hi'] / r['rem']:8.2f}   "
              f"{r['lam'] ** 2 * r['opR2']:.5f}        {r['lam'] ** 2 * em1x(x):.5f}             {tag}")
    nv = [(r["l"], r["g"]) for r in rows if r["low"] > 0]
    print(f"  proved lower bound is NON-VACUOUS (positive) at (l, g) = {nv}; vacuous (negative) at every other table row.")
    print(f"  in general: positive for l >= about {1.001 * C0 / c:.0f} |g|, and >= c g/(2 lam) for l >= l_1(g) = {2.002 * C0 / c:.0f} |g|.")
    print("  grid norms: ||e_l|| " + ", ".join(f"{R[l]['e']:.5f}" for l in LS) + f"  (proved range [{E_lo:.5f}, {C['E_hi']:.5f}])")
    print("              ||coth^1.5 psi_l|| " + ", ".join(f"{R[l]['n3']:.5f}" for l in LS) + f"  (proved <= {N3_hi:.5f})")
    print(f"  lam||D0||_H (box [{YBOT:g}, 40], dy 0.04): " + ", ".join(f"{lam_of(l) * O[l]['D']:.5f}" for l in LS) + "  (Lemma A: <= 1.41421)")
    print("  information only, lam||D0||_H on box [-4.8, 40], dy 0.04: " + ", ".join(f"{info48[l]:.5f}" for l in LS))
    print(f"  N1 rem(g=8)/rem(g=1): {g2}")
    print(f"  N2 rel. diff to ODD-SIGN table: {n2}")
    print(f"  N3 expm cross-check: {n3}")
    print(f"  N4 slack rem_hi/rem (g=1): {slack}")
    print(f"  N7 grid stability at l=1024: {n7}")
    print(f"  PLANT A (cocycle factor dropped): planted 'remainder' {[f'{row[(l, 1.0)]['plant_rem']:.3e}' for l in LS]} vs rem_hi "
          f"{[f'{row[(l, 1.0)]['rem_hi']:.3e}' for l in LS]}; violated {list(pA_viol.values())}; caught = {lines['PLANT_A_caught']}")
    print(f"  PLANT B (g^2/lam^3): planted {[f'{pB[l]:.3e}' for l in LS]} vs rem {[f'{row[(l, 1.0)]['rem']:.3e}' for l in LS]}; "
          f"violated {list(pB_viol.values())}; caught = {lines['PLANT_B_caught']}")
    print("== pass lines ==")
    for k, v in lines.items():
        print(f"  {k}: {'PASS' if v else 'FAIL'}")
    ok = all(lines.values())
    print("ALL LINES HELD" if ok else "SOME LINES FAILED", f"({time.time() - t0:.0f} s)")
    if len(sys.argv) > 1:
        with open(sys.argv[1], "w") as f:
            json.dump({"certified": C, "l1_g1": l1(1), "l1_g8": l1(8), "rows": rows,
                       "opnorms": {str(l): {str(k): v for k, v in O[l].items()} for l in LS},
                       "grid": {str(l): {str(k): v for k, v in R[l].items()} for l in LS},
                       "info_D_box48": {str(l): info48[l] for l in LS}, "g2": {str(k): v for k, v in g2.items()},
                       "n2": {str(k): v for k, v in n2.items()}, "n3": n3, "n7": n7,
                       "slack": {str(k): v for k, v in slack.items()},
                       "plantA_violated": {str(k): bool(v) for k, v in pA_viol.items()},
                       "plantB": {str(l): pB[l] for l in LS}, "plantB_violated": {str(k): bool(v) for k, v in pB_viol.items()},
                       "lines": {k: bool(v) for k, v in lines.items()}}, f, indent=1)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
