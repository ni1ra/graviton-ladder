"""ALL-T: first-order floor and remainder at every tau != 0, both parities.

Pass lines: ..\\logs\\01_prefixed_pass_lines.log (fixed before the first run).
Part C: the lanes' own ball arithmetic re-run (even_passage / odd_sign certified(), even_remainder / odd_remainder certified()),
        then c(tau) = [J0^2 t - J0 J1 t^3/3 + J1^2 t^5/20]/(2 Q16 Den E), t = min(|tau|, s*), in ball arithmetic.
Part N: the lanes' float64 grid (even_remainder.setup, odd_remainder.setup), with T(tau) rebuilt for each tau.
Usage: ..\\..\\even-nonzero\\venv\\Scripts\\python.exe -B all_t.py ..\\out\\all_t.json
Reads sibling lane scripts read-only (no bytecode written); writes only the json named on the command line.
"""
import json
import math
import os
import sys
import time

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
for d in ("even-passage", "odd-sign", "even-remainder", "odd-remainder"):
    sys.path.insert(0, os.path.join(HERE, "..", "..", d, "scripts"))

import numpy as np
from flint import arb, ctx

import even_passage as ep
import odd_sign as osn
import even_remainder as er
import odd_remainder as orr

ctx.prec = 100
PI = math.pi
TAUS = [0.05, 0.5, 1.7771, PI, 2 * PI, 20.0]
TAU_NAMES = ["0.05", "0.5", "1.7771", "pi", "2pi", "20"]
LS = [16, 256, 4096]
SQ2 = math.sqrt(2.0)
lo, hi, lam_of = ep.lo, ep.hi, ep.lam_of
MOD = {"even": er, "odd": orr}


# ------------------------------------------------------------------ Part C
def consts(parity):
    if parity == "even":
        out, unif, _, _ = ep.certified()
        C = er.certified()
        K = 0.5 * (1 + C["Rb"])
    else:
        out, unif, _, _ = osn.certified()
        C = orr.certified()
        K = C["K"]
    return dict(J0=unif["J0_lo"], J1=unif["J1_hi"], ss=unif["sstar"], Q=out["Q_16"], Den=unif["Den_hi"], E=unif["E_hi"],
                c_pi=lo(unif["c"]), M_pi=lo(unif["M_lo"]), K=K, E_lo=C["E_lo"], E_hi=C["E_hi"], Den_hi=hi(unif["Den_hi"]),
                c_remlane=C["c"])


def consts_even2():
    """Even family from even-second's 40-digit enclosures of the seven inputs; the chain is re-derived here."""
    with open(os.path.join(HERE, "..", "..", "even-second", "out", "values.json")) as f:
        V = json.load(f)
    up = lambda key: arb(V[key]["hi"]) * (1 + arb(10) ** -38)
    dn = lambda key: arb(V[key]["lo"]) * (1 - arb(10) ** -38)
    A_lo, A_hi, G, Q, Sp, C, s3, u2 = dn("A_inf"), up("A_16"), up("G_hi"), up("Q_16"), up("Sprime"), up("C"), up("s_3"), up("norm_u_sq")
    J0 = A_lo * (-(ep.EPS * G * G / A_lo)).exp()
    J1 = A_hi / (arb(1).exp() * ep.EPS)
    ss = (2 * J0 / J1).sqrt()
    Den, E = (Sp * A_hi).sqrt(), (C * A_hi / 2).sqrt()
    M = 8 * J0 * J0 * ss / (15 * Q)
    E_lo = (A_lo * (-(2 * ep.EPS * G * G / A_lo)).exp() / (2 * arb.pi())).sqrt()
    Kc = (1 + (s3 * u2 / 2).sqrt() / E_lo) / 2
    k = dict(J0=J0, J1=J1, ss=ss, Q=Q, Den=Den, E=E, c_pi=lo(M / (2 * Den * E)), M_pi=lo(M), K=hi(Kc), E_lo=lo(E_lo), E_hi=hi(E),
             Den_hi=hi(Den), c_remlane=lo(M / (2 * Den * E)))
    agree = max(abs(k["c_pi"] / float(V["c"]["lo"]) - 1), abs(hi(Kc * arb.pi() ** 2) / float(V["C0"]["hi"]) - 1))
    return k, agree


def poly_ball(k, tau):
    """int_0^t (J0 - s^2 J1/2)^2 ds at t = min(|tau|, s*_lo), as a ball. Valid lower bound for the true integral of the
    floor because the integrand is increasing in J0 and decreasing in J1 while J0 - s^2 J1/2 >= 0, i.e. s <= s*."""
    t = min(abs(tau), lo(k["ss"]))
    t = arb(t)
    return k["J0"] ** 2 * t - k["J0"] * k["J1"] * t ** 3 / 3 + k["J1"] ** 2 * t ** 5 / 20


def M_lo(k, tau):
    return lo(poly_ball(k, tau) / k["Q"])


def c_of(k, tau):
    return lo(poly_ball(k, tau) / k["Q"] / (2 * k["Den"] * k["E"]))


def c_quad(k, tau, n=200001):
    """Independent float check of the closed form: trapezoid of max(J0 - s^2 J1/2, 0)^2 on [0, |tau|]."""
    J0, J1 = float(k["J0"].mid()), float(k["J1"].mid())
    s = np.linspace(0.0, abs(tau), n)
    f = np.maximum(J0 - s * s * J1 / 2, 0.0) ** 2
    I = float(np.sum((f[1:] + f[:-1]) / 2) * (s[1] - s[0]))
    return I / float(k["Q"].mid()) / (2 * float(k["Den"].mid()) * float(k["E"].mid()))


def eta_of(k):
    return (15 * k["c_pi"] / (8 * lo(k["ss"]))) / (SQ2 * k["K"]) * (1 + 1e-12)


def l1_of(k, g, tau):
    return max(16, math.ceil(2 * math.exp(eta_of(k)) * k["K"] * tau * tau * abs(g) / c_of(k, tau)))


def em1x(x):
    return (math.expm1(x) - x) * (1 + 1e-12) + 1e-300 if x > 1e-4 else 0.5 * x * x * math.exp(x) * (1 + 1e-12)


# ------------------------------------------------------------------ Part N
def retau(st, parity, tau):
    """Rebuild T(tau) z = A z + B zbar, T = int_0^tau e^(sX0) D e^(-sX0) ds (oriented), in the lane's conventions."""
    p, M = st["p"], st["M"]

    def F(theta):                                        # int_0^tau e^(-i s theta) ds, either sign of tau
        th = tau * theta
        return tau * np.exp(-0.5j * th) * np.sinc(th / (2 * PI))

    st["A"] = -0.5j * M * F(p[:, None] - p[None, :])
    st["B"] = (0.5j if parity == "even" else -0.5j) * M * F(p[:, None] + p[None, :])
    MOD[parity].TAU = tau                                # fwd / back0 of the lane read the module global


def tau_free(st, parity):
    """Quantities with no tau in them: z, Den, the quadratic form of the integrand."""
    p, a, dy, lam = st["p"], st["a"], st["dy"], st["lam"]
    if parity == "even":
        z = 1j * a / np.sqrt(2 * p)
        Wq = lam * st["W"]                               # lam W_l in the eigenbasis
    else:
        z = (np.sqrt(p) * a / SQ2).astype(complex)
        Wq = lam * st["Nt"]                              # lam N_l
    # Den^2 = <z, 2 sinh(2 pi h) z>, overflow-free
    z2 = np.abs(z) ** 2
    nz = z2 > 0                                          # z underflows to exactly 0 long before exp(2 pi p) overflows
    den = math.sqrt(dy * float(np.sum(np.exp(2 * PI * p[nz]) * (-np.expm1(-4 * PI * p[nz])) * z2[nz])))
    e = math.sqrt(dy * float(np.sum(st["coth"] * np.abs(z) ** 2)))
    return z, Wq, den, e


def m_of(st, Wq, tau):
    """lam * int_0^tau <a_s, W a_s> ds, a_s = cos(s h) a, oriented; composite Gauss-Legendre."""
    p, a, dy = st["p"], st["a"], st["dy"]
    g, wq = np.polynomial.legendre.leggauss(48)
    npan = max(1, int(math.ceil(abs(tau) / 0.8)))
    hpan = tau / npan
    tot = 0.0
    for j in range(npan):
        for uu, ww in zip(g, wq):
            s = hpan * (j + 0.5 * (uu + 1))
            v = np.cos(s * p) * a
            tot += 0.5 * hpan * ww * dy * float(v @ (Wq @ v))
    return tot


def cell(st, parity, tau, z, Wq, den, e):
    retau(st, parity, tau)
    res = MOD[parity].vectors(st, [])
    Tz = st["A"] @ z + st["B"] @ np.conj(z)
    im = st["lam"] * st["dy"] * float(np.imag(np.vdot(z, Tz)))      # lam * Im<z,Tz>: st["A"], st["B"] carry the true W or N
    m = m_of(st, Wq, tau)
    return dict(r=st["lam"] * res["beta1"] / res["e"], e=res["e"], m=m, im=im, chain=abs(m) / (2 * den * e))


def main():
    t0 = time.time()
    lines = {}
    K = {par: consts(par) for par in ("even", "odd")}
    print("== Part C (ball arithmetic re-run from the lanes' certified(); outward-rounded floats) ==")
    cq = {}
    for par in ("even", "odd"):
        k = K[par]
        print(f"  [{par}] J0_lo={lo(k['J0']):.9f} J1_hi={hi(k['J1']):.9f} s*={lo(k['ss']):.9f} Q16={hi(k['Q']):.6f} Den_hi={k['Den_hi']:.6f} "
              f"E_hi={k['E_hi']:.6f} E_lo={k['E_lo']:.6f} c_pi={k['c_pi']:.6e} K={k['K']:.6f} eta={eta_of(k):.6e} e^eta={math.exp(eta_of(k)):.7f}")
        print(f"         small-tau slope c(tau)/|tau| -> 15 c_pi/(8 s*) = {15 * k['c_pi'] / (8 * lo(k['ss'])):.6e}"
              f"   (= J0^2/(2 Q Den E) = {lo(k['J0'] ** 2 / k['Q'] / (2 * k['Den'] * k['E'])):.6e})")
        for nm, tau in zip(TAU_NAMES, TAUS):
            cq[(par, nm)] = abs(c_quad(k, tau) / c_of(k, tau) - 1)
            print(f"         tau={nm:>6}: c(tau)={c_of(k, tau):.6e}  c/c_pi={c_of(k, tau) / k['c_pi']:.6f}  quad rel.diff={cq[(par, nm)]:.1e}  "
                  f"l_1(1,tau)={l1_of(k, 1, tau)}")
    lines["C1"] = (hi(K["even"]["ss"]) < PI and hi(K["odd"]["ss"]) < PI and abs(K["even"]["c_pi"] / 2.110e-3 - 1) < 5e-4
                   and abs(K["odd"]["c_pi"] / 2.038e-3 - 1) < 5e-4 and K["even"]["K"] <= 1.5715 and K["odd"]["K"] <= 1.5716
                   and max(cq.values()) < 1e-9)
    l1pi = {par: l1_of(K[par], 1, PI) for par in K}
    l1h = {par: l1_of(K[par], 1, 0.5) for par in K}
    l1tp = {par: l1_of(K[par], 1, 2 * PI) for par in K}
    lines["C2"] = (14714 <= l1pi["even"] <= 14716 and 15237 <= l1pi["odd"] <= 15238 and all(600 <= v <= 900 for v in l1h.values())
                   and all(55000 <= v <= 65000 for v in l1tp.values()))

    print("== Part N: grid ==", flush=True)
    ref = {("even", 16): 0.121788, ("even", 256): 0.124501, ("even", 4096): 0.124673,
           ("odd", 16): 0.297406, ("odd", 256): 0.304791, ("odd", 4096): 0.305262}
    T, FREE, SYM, PL, REM = {}, {}, {}, {}, {}
    for par in ("even", "odd"):
        k = K[par]
        for l in LS:
            st = MOD[par].setup(l)
            z, Wq, den, e = tau_free(st, par)
            FREE[(par, l)] = dict(den=den, e=e, n=st["n"])
            for nm, tau in zip(TAU_NAMES, TAUS):
                T[(par, l, nm)] = cell(st, par, tau, z, Wq, den, e)
            PL[(par, l)] = cell(st, par, 0.002, z, Wq, den, e)                      # plant-only tau
            if l == 256:
                for nm, tau in (("0.5", 0.5), ("2pi", 2 * PI)):
                    SYM[(par, nm)] = cell(st, par, -tau, z, Wq, den, e)
            if l == 4096:
                for nm, tau in (("0.5", 0.5), ("2pi", 2 * PI)):
                    retau(st, par, tau)
                    v = MOD[par].vectors(st, [1.0])
                    REM[(par, nm)] = dict(tau=tau, actual=v[1.0]["beta"] / v["e"], first=v["beta1"] / v["e"], rem=v[1.0]["rem"] / v["e"])
            print(f"  [{par}] l={l} n={st['n']} Den={den:.5f} ||e||={e:.5f}  ({time.time() - t0:.0f} s)", flush=True)
            del st, Wq
        # operator norms on H at l = 4096, box [-3, 40], dy = 0.04
        st = er.setup(4096, dy=0.04, ybot=-3.0) if par == "even" else orr.setup(4096, dy=0.04)
        for nm, tau in (("0.5", 0.5), ("2pi", 2 * PI)):
            retau(st, par, tau)
            o = MOD[par].opnorms(st, [1.0])
            REM[(par, nm)].update(opS=o[1.0]["S"], opR2=o[1.0]["R2"], lamD=st["lam"] * o["D"])
        del st

    F = {x: [] for x in ("F1", "F2", "F3", "F4", "F5", "F6", "F7", "R1", "R2", "R3")}
    print("== table: r = lam ||beta1hat(tau) e_l||/||e_l|| (grid) beside the proved c(tau) ==")
    print("  parity      tau      c(tau) proved      r(l=16)       r(l=256)      r(l=4096)    min r/c   lam*m(l=4096)   M_lo(tau)")
    for par in ("even", "odd"):
        k = K[par]
        for nm, tau in zip(TAU_NAMES, TAUS):
            c, Ml = c_of(k, tau), M_lo(k, tau)
            rs = [T[(par, l, nm)]["r"] for l in LS]
            print(f"  {par:5s}  {nm:>8}    {c:.6e}     {rs[0]:.8f}    {rs[1]:.8f}    {rs[2]:.8f}   {min(rs) / c:8.1f}   "
                  f"{T[(par, 4096, nm)]['m']:.6f}       {Ml:.6f}")
            for l in LS:
                t = T[(par, l, nm)]
                F["F1"].append(t["r"] >= c)
                F["F2"].append(t["m"] >= Ml and t["r"] >= t["chain"] * (1 - 1e-9) and abs(t["im"] + t["m"] / 2) <= 1e-6 * abs(t["m"]))
                F["F4"].append(3 <= t["r"] / c <= 3000)
        for l in LS:
            f = FREE[(par, l)]
            F["F3"].append(f["den"] <= k["Den_hi"] and k["E_lo"] <= f["e"] <= k["E_hi"])
            F["F5"].append(abs(T[(par, l, "pi")]["r"] / ref[(par, l)] - 1) < 1e-5)
            F["F6"].append(0.095 <= T[(par, l, "0.05")]["r"] / T[(par, l, "0.5")]["r"] <= 0.20)
    print("  tau-free denominators (grid) vs proved:  " + "; ".join(
        f"[{par} l={l}] Den={FREE[(par, l)]['den']:.5f}<={K[par]['Den_hi']:.5f}, {K[par]['E_lo']:.5f}<=||e||={FREE[(par, l)]['e']:.5f}<={K[par]['E_hi']:.5f}"
        for par in ("even", "odd") for l in LS))
    print("  F5 tau=pi vs lanes (rel.diff): " + ", ".join(
        f"[{par} {l}] {abs(T[(par, l, 'pi')]['r'] / ref[(par, l)] - 1):.1e}" for par in ("even", "odd") for l in LS))
    print("  F6 r(0.05)/r(0.5): " + ", ".join(
        f"[{par} {l}] {T[(par, l, '0.05')]['r'] / T[(par, l, '0.5')]['r']:.5f}" for par in ("even", "odd") for l in LS)
        + f";  c(0.05)/c(0.5) = {c_of(K['even'], 0.05) / c_of(K['even'], 0.5):.5f}")
    print("  F7 sign symmetry at l = 256:")
    for par in ("even", "odd"):
        for nm in ("0.5", "2pi"):
            a, b = SYM[(par, nm)], T[(par, 256, nm)]
            dr, dm = abs(a["r"] / b["r"] - 1), abs(a["m"] / b["m"] + 1)
            F["F7"].append(dr <= 1e-8 and dm <= 1e-9 and abs(a["im"] + a["m"] / 2) <= 1e-6 * abs(a["m"]))
            print(f"    [{par}] tau=-{nm}: r={a['r']:.10f} vs r(+tau)={b['r']:.10f} (rel {dr:.1e});  m(-tau)={a['m']:+.8f}, m(+tau)={b['m']:+.8f}; "
                  f"Im<z,Tz>={a['im']:+.8f}")

    print("== remainder at g = 1, l = 4096 ==")
    lam = lam_of(4096)
    print("  parity  tau     actual ||b e||/||e||   g||b1 e||/||e||     actual rem      K(e^x-1-x)     slack    proved low      c(tau)/(2lam)   "
          "l_1(1,tau)  ||R2||_H       e^x-1-x        ||S||_H")
    for par in ("even", "odd"):
        k = K[par]
        for nm in ("0.5", "2pi"):
            R = REM[(par, nm)]
            x = SQ2 * abs(R["tau"]) / lam
            rh = k["K"] * em1x(x)
            c = c_of(k, R["tau"])
            low, half, l1 = c / lam - rh, c / (2 * lam), l1_of(k, 1, R["tau"])
            R.update(rem_hi=rh, low=low, half=half, l1=l1, x=x)
            F["R1"].append(R["rem"] <= rh and 20 <= rh / R["rem"] <= 1e5)
            F["R2"].append(R["actual"] >= low and ((low >= half > 0 and l1 < 4096) if nm == "0.5" else (low < 0 and l1 > 4096)))
            F["R3"].append(R["opR2"] <= em1x(x) and R["opS"] <= math.exp(x) * (1 + 1e-12))
            print(f"  {par:5s}  {nm:>4}    {R['actual']:.9e}       {R['first']:.9e}    {R['rem']:.6e}    {rh:.6e}  {rh / R['rem']:8.1f}   {low:+.6e}   "
                  f"{half:.6e}    {l1:6d}     {R['opR2']:.6e}   {em1x(x):.6e}   {R['opS']:.9f}")
    for kk, v in F.items():
        lines[kk] = all(v)
        if not all(v):
            print(f"  {kk} detail: {v}")

    # ---- lines added after run 1 (log 03): even family with even-second's enclosures; monotonicity
    k2, agree = consts_even2()
    print("== even family with even-second's enclosures (chain re-derived here) ==")
    print(f"  c_pi={k2['c_pi']:.9e} K={k2['K']:.6f} K pi^2={k2['K'] * PI * PI:.6f} s*={lo(k2['ss']):.9f} Den_hi={k2['Den_hi']:.6f} E_hi={k2['E_hi']:.6f} "
          f"E_lo={k2['E_lo']:.6f} eta={eta_of(k2):.6e}; agreement with that lane's c and C0: {agree:.1e}")
    print(f"  small-tau slope 15 c_pi/(8 s*) = {15 * k2['c_pi'] / (8 * lo(k2['ss'])):.6e}")
    e2f = []
    print("       tau      c(tau) proved    l_1(1,tau)   min_l r/c")
    for nm, tau in zip(TAU_NAMES, TAUS):
        c, Ml = c_of(k2, tau), M_lo(k2, tau)
        rs = [T[("even", l, nm)]["r"] for l in LS]
        print(f"  {nm:>8}    {c:.6e}    {l1_of(k2, 1, tau):8d}    {min(rs) / c:8.1f}")
        e2f += [T[("even", l, nm)]["r"] >= c and T[("even", l, nm)]["m"] >= Ml and 3 <= T[("even", l, nm)]["r"] / c <= 3000 for l in LS]
    e2f += [FREE[("even", l)]["den"] <= k2["Den_hi"] and k2["E_lo"] <= FREE[("even", l)]["e"] <= k2["E_hi"] for l in LS]
    l1p, l1hh = l1_of(k2, 1, PI), l1_of(k2, 1, 0.5)
    lines["E2_C"] = (agree < 1e-12 and 2.1193e-3 <= k2["c_pi"] <= 2.1194e-3 and 1.5692 <= k2["K"] <= 1.5694 and 14625 <= l1p <= 14635
                     and 600 <= l1hh <= 900)
    lines["E2_F"] = all(e2f)
    e2r = []
    for nm in ("0.5", "2pi"):
        R = REM[("even", nm)]
        rh = k2["K"] * em1x(R["x"])
        c = c_of(k2, R["tau"])
        low, half = c / lam - rh, c / (2 * lam)
        e2r.append(R["rem"] <= rh and R["actual"] >= low and ((low >= half) if nm == "0.5" else (low < 0)))
        print(f"  remainder, tau={nm}: actual rem {R['rem']:.6e} <= {rh:.6e};  actual {R['actual']:.9e} >= proved low {low:+.6e};  c/(2lam) = {half:.6e}")
    lines["E2_R"] = all(e2r)
    mono = []
    for par in ("even", "odd"):
        for l in LS:
            ms = [PL[(par, l)]["m"]] + [T[(par, l, nm)]["m"] for nm in TAU_NAMES]
            mono.append(all(b > a > 0 for a, b in zip(ms, ms[1:])))
        cs = [c_of(K[par], t) for t in [0.002] + TAUS]
        mono.append(all(b >= a > 0 for a, b in zip(cs, cs[1:])))
    cs = [c_of(k2, t) for t in [0.002] + TAUS]
    mono.append(all(b >= a > 0 for a, b in zip(cs, cs[1:])))
    lines["MONO"] = all(mono)
    print(f"  MONO: lam*m at l=4096 along tau = 0.002, 0.05, 0.5, 1.7771, pi, 2pi, 20: "
          + "; ".join(f"[{par}] " + ", ".join(f"{v:.6f}" for v in [PL[(par, 4096)]['m']] + [T[(par, 4096, nm)]['m'] for nm in TAU_NAMES])
                      for par in ("even", "odd")))

    # plant: c(tau) continued as the constant c_pi below s*
    print("== PLANT: c(tau) := c_pi for all tau (constant continued below s*) ==")
    pm, pf, p05 = [], [], []
    for par in ("even", "odd"):
        k = K[par]
        for l in LS:
            t05, t002 = T[(par, l, "0.05")], PL[(par, l)]
            pm.append(t05["m"] < k["M_pi"])
            pf.append(t002["r"] < k["c_pi"])
            p05.append(t05["r"] < k["c_pi"])
            print(f"  [{par} l={l}] m-line at 0.05: lam*m={t05['m']:.6f} vs planted M={k['M_pi']:.6f} violated={pm[-1]};  floor at 0.002: r={t002['r']:.6e} "
                  f"vs planted c={k['c_pi']:.6e} violated={pf[-1]};  floor at 0.05: r={t05['r']:.6e} violated={p05[-1]}"
                  f"   [true c(0.002)={c_of(k, 0.002):.6e} holds={t002['r'] >= c_of(k, 0.002)}]")
    lines["PLANT_caught"] = all(pm) and all(pf)
    pred = all(2.2e-3 <= T[(par, l, "0.05")]["r"] <= 3e-2 for par in ("even", "odd") for l in LS)
    print(f"  caught on m-line (all cells) = {all(pm)}; caught on final floor at tau=0.002 (all cells) = {all(pf)}; "
          f"prediction 'not caught on final floor at 0.05, r in [2.2e-3, 3e-2]' held = {pred and not any(p05)}")
    print("== pass lines ==")
    for kk, v in lines.items():
        print(f"  {kk}: {'PASS' if v else 'FAIL'}")
    ok = all(lines.values())
    print("ALL LINES HELD" if ok else "SOME LINES FAILED", f"({time.time() - t0:.0f} s)")
    if len(sys.argv) > 1:
        with open(sys.argv[1], "w") as f:
            json.dump({"consts": {par: {kk: (v if isinstance(v, float) else {"lo": lo(v), "hi": hi(v)}) for kk, v in K[par].items()} for par in K},
                       "c_tau": {par: {nm: c_of(K[par], tau) for nm, tau in zip(TAU_NAMES, TAUS)} for par in K},
                       "l1_g1": {par: {nm: l1_of(K[par], 1, tau) for nm, tau in zip(TAU_NAMES, TAUS)} for par in K},
                       "table": {f"{a}_{b}_{c}": v for (a, b, c), v in T.items()},
                       "free": {f"{a}_{b}": v for (a, b), v in FREE.items()},
                       "sym": {f"{a}_{b}": v for (a, b), v in SYM.items()},
                       "plant_tau_0.002": {f"{a}_{b}": v for (a, b), v in PL.items()},
                       "rem": {f"{a}_{b}": v for (a, b), v in REM.items()},
                       "plant_prediction_held": bool(pred and not any(p05)),
                       "lines": {kk: bool(v) for kk, v in lines.items()}}, f, indent=1)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
