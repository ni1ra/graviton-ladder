# BEYOND-REF: hostile referee's independent check of BEYOND-NORMAL.
# Pass lines: ..\logs\01_pass_lines.log (fixed before the first run).
# Usage: C:\Python314\python.exe -X utf8 beyond_ref.py ..\out\beyond_ref.json
import sys, json, time
import numpy as np

PI = np.pi
LS = [2, 4, 8, 16, 32, 64, 128, 256]
failed, res = [], {}


def line(name, ok, text, expect_fail=False):
    tag = "PASS " if ok else ("FAIL(predicted) " if expect_fail else "FAIL ")
    print(tag + name + ": " + text, flush=True)
    if not ok and not expect_fail:
        failed.append(name + ": " + text)


def grid(hg, rmax, rin=0.0):
    n = int(round(rmax / hg))
    return rin + hg * np.arange(1, n + 1)


def h2_matrix(l, rho, hg):
    lam2 = l * (l + 1)
    pot = lam2 / np.sinh(rho) ** 2
    main = 2.0 / hg**2 + pot
    main[-1] = 1.0 / hg**2 + pot[-1]
    off = -np.ones(len(rho) - 1) / hg**2
    return np.diag(main) + np.diag(off, 1) + np.diag(off, -1), pot


def w_matrix(l, rho, hg):
    lr = np.log(np.tanh(rho / 2))
    lo, hi = np.minimum.outer(lr, lr), np.maximum.outer(lr, lr)
    return np.exp((l + 1) * lo - l * hi) * (1.0 - np.exp((2 * l + 1) * hi)) / (2 * l + 1) * hg


def coth(x):
    return 1.0 / np.tanh(x)


def fun(v, d):
    return (v * d) @ v.T


def inv_sqrt(m):
    w, v = np.linalg.eigh(0.5 * (m + m.T))
    return fun(v, 1.0 / np.sqrt(w))


def sqrt_psd(m):
    w, v = np.linalg.eigh(0.5 * (m + m.T))
    return fun(v, np.sqrt(np.clip(w, 0.0, None))), w[0]


class Sec:
    """All matrices in the vacuum eigenbasis of h^2. forms(beta) -> (Mq, Mp), mu = 1/2 q Mq q + 1/2 p Mp p."""

    def __init__(self, l, hg=0.025, rmax=50.0, cap=None):
        self.l, self.hg = l, hg
        # inner Dirichlet wall where lambda^2/sinh^2 > cap (conditioning; see logs\04_fix_round_pass_lines.log)
        self.rin = 0.0 if cap is None else hg * np.ceil(np.arcsinh(np.sqrt(l * (l + 1) / cap)) / hg)
        self.rho = grid(hg, rmax, self.rin)
        self.h2, self.pot = h2_matrix(l, self.rho, hg)
        e, self.V = np.linalg.eigh(self.h2)
        self.nu = np.sqrt(e)
        self.npot = 1.0 / (1.0 + np.cosh(self.rho))
        self.Nt = self.V.T @ (self.npot[:, None] * self.V)

    def vac(self, beta):
        c = coth(beta * self.nu / 2)
        return np.diag(self.nu * c), np.diag(c / self.nu)

    def odd(self, gcoef=1.0):
        e, v = np.linalg.eigh(np.diag(self.nu**2) + gcoef * self.Nt)
        om = np.sqrt(e)

        def forms(beta):
            c = coth(beta * om / 2)
            return fun(v, om * c), fun(v, c / om)
        return forms, (om, v)

    def even(self):
        w = w_matrix(self.l, self.rho, self.hg)
        self.wnorm = np.linalg.eigvalsh(0.5 * (w + w.T))[-1]
        wt = self.V.T @ w @ self.V
        x = self.nu[:, None] * (np.eye(len(self.nu)) + wt) * self.nu[None, :]
        e, v = np.linalg.eigh(0.5 * (x + x.T))
        om = np.sqrt(e)

        def forms(beta):
            c = coth(beta * om / 2)
            mp = fun(v, om * c) / self.nu[:, None] / self.nu[None, :]
            mq = fun(v, c / om) * self.nu[:, None] * self.nu[None, :]
            return mq, mp
        return forms, (om, v)


def slot_ratios(f0, fg):
    out = []
    for s in (0, 1):
        d = 1.0 / np.sqrt(np.diag(f0[s]))
        x = fg[s] * d[:, None] * d[None, :]
        out.append(np.linalg.eigvalsh(0.5 * (x + x.T)))
    return out


def one_minus_aff(sec, formsg):
    r = slot_ratios(sec.vac(PI), formsg(PI))
    tot = sum(np.sum(np.log1p((np.sqrt(x) - 1.0) ** 2 / (2.0 * np.sqrt(x)))) for x in r)
    return -np.expm1(-0.5 * tot)


def roots(ref, diagonal=False):
    if diagonal:
        return np.diag(1.0 / np.sqrt(np.diag(ref[0]))), np.diag(1.0 / np.sqrt(np.diag(ref[1])))
    return inv_sqrt(ref[0]), inv_sqrt(ref[1])


def shat(rt, f):
    """S-hat = M^(-1/2) (mu + i sigma/2) M^(-1/2), M = 2 mu_ref, made real by the unitary diag(1, i)."""
    rq, rp = rt
    x = rq @ rp
    return 0.5 * np.block([[rq @ f[0] @ rq, x], [x.T, rp @ f[1] @ rp]])


def s0half_exact(nu):
    """Square root of the vacuum S-hat in its own inner product: blocks 1/2 [[1, t], [t, 1]], t = tanh(pi nu)."""
    t = np.tanh(PI * nu)
    a, b = np.sqrt((1 + t) / 2), np.sqrt((1 - t) / 2)
    return np.block([[np.diag((a + b) / 2), np.diag((a - b) / 2)], [np.diag((a - b) / 2), np.diag((a + b) / 2)]])


def ay_quantity(ref, fa, fb, diagonal=False, sa=None):
    rt = roots(ref, diagonal)
    la = 0.0
    if sa is None:
        sa, la = sqrt_psd(shat(rt, fa))
    sb, lb = sqrt_psd(shat(rt, fb))
    return float(np.sum((sa - sb) ** 2)), float(min(la, lb))


def zero_mode(pot_total, hg):
    n = len(pot_total)
    phi = np.zeros(n)
    phi[0] = 1.0
    prev = 0.0
    for i in range(n - 1):
        nxt = (2.0 + hg * hg * pot_total[i]) * phi[i] - prev
        prev = phi[i]
        phi[i + 1] = nxt
    return phi / ((phi[-1] - phi[-2]) / hg)


def sector(l):
    lam2 = l * (l + 1)
    s = Sec(l, cap=5000.0 if l >= 128 else None)
    out = {}
    f0 = s.vac(2 * PI)
    s0h = s0half_exact(s.nu)
    fo, _ = s.odd()
    fe, _ = s.even()
    cl = np.linalg.eigvalsh(s.Nt / s.nu[:, None] / s.nu[None, :])[-1]
    w = s.wnorm
    print(f"  l={l}: rho_in = {s.rin:.3f}", flush=True)
    line(f"R3 l={l}", cl <= 3.52 and abs((l + 1) * cl - 1) <= 1e-3 and abs((2 * l + 3) * w - 1) <= 1e-3,
         f"C_l = {cl:.7f}, (l+1)C_l = {(l+1)*cl:.6f}, (2l+3)||W|| = {(2*l+3)*w:.6f}")
    out.update(C_l=float(cl), W_norm=float(w))
    for par, fg, bq, bp in (("odd", fo(2 * PI), (1.0, 1.0 + cl), (1.0 / (1.0 + cl), 1.0)),
                            ("even", fe(2 * PI), (1.0 / (1.0 + w), 1.0), (1.0, 1.0 + w))):
        rq, rp = slot_ratios(f0, fg)
        ok = rq[0] >= bq[0] - 1e-8 and rq[-1] <= bq[1] + 1e-8 and rp[0] >= bp[0] - 1e-8 and rp[-1] <= bp[1] + 1e-8
        line(f"R2 {par} l={l}", ok, f"q-slot [{rq[0]:.8f}, {rq[-1]:.8f}] in [{bq[0]:.8f}, {bq[1]:.8f}]; "
             f"pi-slot [{rp[0]:.8f}, {rp[-1]:.8f}] in [{bp[0]:.8f}, {bp[1]:.8f}]")
        q0, m0 = ay_quantity(f0, f0, fg, diagonal=True, sa=s0h)
        qg, mg = ay_quantity(fg, f0, fg)
        line(f"R4 {par} l={l}", min(m0, mg) >= -1e-9, f"lambda_min(S-hat) = {min(m0, mg):.3e}")
        a = one_minus_aff(s, fo if par == "odd" else fe)
        line(f"R6b {par} l={l}", 0.5 <= qg / q0 <= 2.0,
             f"lambda^2 Q (ref vac) = {lam2*q0:.6f}, (ref g) = {lam2*qg:.6f}, ratio {qg/q0:.5f}; "
             f"lambda^2 (1-A) = {lam2*a:.6f}; Q/(1-A) = {q0/a:.5f}")
        out[par] = {"Q_ref0": q0, "Q_refg": qg, "one_minus_A": float(a),
                    "q_slot": [float(rq[0]), float(rq[-1])], "p_slot": [float(rp[0]), float(rp[-1])]}
    fpl, _ = s.odd(4.0 / l)
    qp, _ = ay_quantity(f0, f0, fpl(2 * PI), diagonal=True, sa=s0h)
    out["plant_Q"] = qp
    print(f"  plant l={l}: lambda^2 Q = {lam2*qp:.6e}", flush=True)
    if l == 2:
        zero_mode_drift(s)
    return out


def zero_mode_drift(s):
    rho, hg = s.rho, s.hg
    q0 = np.exp(-0.5 * (rho - 5.0) ** 2)
    p0 = 0.7 * np.exp(-0.5 * (rho - 6.0) ** 2)
    # even: pi = nu*y (vacuum basis), y'' = -nu(1+W)nu y, y'(0) = -nu q
    phi = zero_mode(s.pot, hg)
    _, (om, v) = s.even()
    a = v.T @ ((s.V.T @ p0) / s.nu)
    b = v.T @ (-s.nu * (s.V.T @ q0))
    ell0 = phi @ p0
    dr = []
    for t in (1.7, 10.0):
        y = v @ (np.cos(om * t) * a + np.sin(om * t) / om * b)
        dr.append(abs(phi @ (s.V @ (s.nu * y)) / ell0 - 1))
    line("R9a even l=2", max(dr) <= 1e-6, f"ell(0) = {ell0*hg:.5f}; relative drift at tau = 1.7, 10: {dr[0]:.2e}, {dr[1]:.2e}")
    res["R9a_even"] = [float(x) for x in dr]
    # odd: q' = pi, pi' = -(h^2+N) q, zero mode of h^2+N
    phig = zero_mode(s.pot + s.npot, hg)
    _, (om, v) = s.odd()
    a = v.T @ (s.V.T @ q0)
    b = v.T @ (s.V.T @ p0)
    ell0 = phig @ p0
    dr = []
    for t in (1.7, 10.0):
        pt = s.V @ (v @ (-om * np.sin(om * t) * a + np.cos(om * t) * b))
        dr.append(abs(phig @ pt / ell0 - 1))
    line("R9a odd l=2", max(dr) <= 1e-6, f"ell(0) = {ell0*hg:.5f}; relative drift at tau = 1.7, 10: {dr[0]:.2e}, {dr[1]:.2e}")
    res["R9a_odd"] = [float(x) for x in dr]


def zero_mode_unbounded():
    vals = {"even/vacuum": [], "odd": []}
    for rmax in (25.0, 50.0, 100.0):
        hg = 0.05
        rho = grid(hg, rmax)
        h2, pot = h2_matrix(2, rho, hg)
        npot = 1.0 / (1.0 + np.cosh(rho))
        for name, mat, pt in (("even/vacuum", h2, pot), ("odd", h2 + np.diag(npot), pot + npot)):
            e, v = np.linalg.eigh(mat)
            nu = np.sqrt(e)
            phi = zero_mode(pt, hg)
            c = v.T @ phi
            vals[name].append(2.0 * hg * float(np.sum(c * c * nu * np.tanh(PI * nu))))
    for name, x in vals.items():
        line(f"R9b {name} l=2", x[1] / x[0] >= 1.7 and x[2] / x[1] >= 1.7,
             f"sup ell^2/mu at box 25, 50, 100 = {x[0]:.3f}, {x[1]:.3f}, {x[2]:.3f} (ratios {x[1]/x[0]:.3f}, {x[2]/x[1]:.3f})")
    res["R9b"] = vals


def discretisation():
    q = {}
    q["fine"] = (res["l=16"]["odd"]["Q_ref0"], res["l=16"]["even"]["Q_ref0"])
    for key, (hg, rmax) in (("coarse", (0.05, 50.0)), ("wide", (0.05, 100.0))):
        s = Sec(16, hg, rmax)
        f0 = s.vac(2 * PI)
        q[key] = (ay_quantity(f0, f0, s.odd()[0](2 * PI), diagonal=True)[0], ay_quantity(f0, f0, s.even()[0](2 * PI), diagonal=True)[0])
    for i, par in enumerate(("odd", "even")):
        d1 = abs(q["coarse"][i] / q["fine"][i] - 1)
        d2 = abs(q["wide"][i] / q["coarse"][i] - 1)
        line(f"R8 {par} l=16", d1 <= 0.02 and d2 <= 0.02,
             f"Q fine/coarse/wide = {q['fine'][i]:.6e}, {q['coarse'][i]:.6e}, {q['wide'][i]:.6e}; step change {d1:.2e}, box change {d2:.2e}")
    res["R8"] = {k: list(v) for k, v in q.items()}


def slope(ls, vals):
    return float(np.polyfit(np.log(ls), np.log(vals), 1)[0])


if __name__ == "__main__":
    t0 = time.time()
    xs = np.linspace(1e-3, 10.0, 200001)
    r = xs * coth(PI * xs) / (1.0 + xs)
    i = int(np.argmin(r))
    line("R1 lane-claim", r[i] >= 1 / PI - 1e-6,
         f"min x coth(pi x)/(1+x) = {r[i]:.6f} at x = {xs[i]:.4f}; lane's constant 1/pi = {1/PI:.6f}; sup = {r.max():.6f}", expect_fail=True)
    res["R1_min"] = float(r[i])
    for l in LS:
        res[f"l={l}"] = sector(l)
        print(f"  [l={l} done, {time.time()-t0:.0f} s]", flush=True)
    lam2 = {l: l * (l + 1) for l in LS}
    a16 = lam2[16] * res["l=16"]["odd"]["one_minus_A"]
    e2, e16 = res["l=2"]["even"]["one_minus_A"], res["l=16"]["even"]["one_minus_A"]
    line("R0", abs(a16 / 0.04602 - 1) <= 5e-3 and abs(e2 / 7.825418e-04 - 1) <= 1e-3 and abs(e16 / 3.885959e-05 - 1) <= 1e-3,
         f"odd lambda^2(1-A) l=16 = {a16:.5f} (lane 0.04602); even 1-A l=2 = {e2:.6e} (lane 7.825418e-04), l=16 = {e16:.6e} (lane 3.885959e-05)")
    tail = [l for l in LS if l >= 32]
    for par in ("odd", "even"):
        for ref, name in (("Q_ref0", "R5"), ("Q_refg", "R6")):
            v = [lam2[l] * res[f"l={l}"][par][ref] for l in tail]
            s = slope(tail, v)
            line(f"{name} {par}", abs(s) <= 0.15 and v[-1] >= 1e-3,
                 f"slope of log(lambda^2 Q) on l=32..256 = {s:+.4f}; lambda^2 Q = " + ", ".join(f"{x:.6f}" for x in v) + (" -> DIVERGES" if abs(s) <= 0.15 else " -> NOT-DIVERGES"))
            res[f"{name}_{par}_slope"] = s
        rat = [res[f"l={l}"][par]["Q_ref0"] / res[f"l={l}"][par]["one_minus_A"] for l in tail]
        target = 7 * 1.2020569031595942 / (16 * PI**2) if par == "odd" else 0.0135994
        tol = 0.03 if par == "odd" else 0.05
        a256 = lam2[256] * res["l=256"][par]["one_minus_A"]
        line(f"R7 {par}", max(rat) / min(rat) <= 1.10 and abs(a256 / target - 1) <= tol,
             f"Q/(1-A) on l=32..256 = " + ", ".join(f"{x:.5f}" for x in rat) + f"; lambda^2(1-A) at 256 = {a256:.6f} vs record {target:.6f}")
    vp = [lam2[l] * res[f"l={l}"]["plant_Q"] for l in tail]
    sp = slope(tail, vp)
    caught = sp <= -1.5
    print(("CAUGHT " if caught else "MISSED ") + f"PLANT (odd coupling 4/l): slope {sp:+.4f} -> " + ("NOT-DIVERGES" if abs(sp) > 0.15 else "DIVERGES"), flush=True)
    res["plant_slope"], res["plant_caught"] = sp, bool(caught)
    discretisation()
    s = Sec(64, cap=1000.0)
    f0, s0h = s.vac(2 * PI), s0half_exact(s.nu)
    devs = []
    for par, fm in (("odd", s.odd()[0]), ("even", s.even()[0])):
        qc = ay_quantity(f0, f0, fm(2 * PI), diagonal=True, sa=s0h)[0]
        ac = one_minus_aff(s, fm)
        devs += [abs(qc / res["l=64"][par]["Q_ref0"] - 1), abs(ac / res["l=64"][par]["one_minus_A"] - 1)]
    line("R10 wall harmless l=64", max(devs) <= 1e-5, f"rho_in = {s.rin:.3f}; relative change of Q odd, 1-A odd, Q even, 1-A even: " + ", ".join(f"{x:.2e}" for x in devs))
    res["R10"] = [float(x) for x in devs]
    zero_mode_unbounded()
    growth = 0.0
    print("Tower partial sums of (2l+1) Q_l are not formed (8 sample l only); per-decade growth of the sum = 2 ln(10) * lim lambda^2 Q:")
    for par in ("odd", "even"):
        c = lam2[256] * res["l=256"][par]["Q_ref0"]
        print(f"  {par}: lambda^2 Q_256 = {c:.6f}; growth per decade of L ~ {2*np.log(10)*c:.5f}")
    print("FAILED LINES: " + ("none" if not failed else ""))
    for f in failed:
        print("  " + f)
    print(f"plant caught: {caught}; total {time.time()-t0:.0f} s")
    res["failed"] = failed
    if len(sys.argv) > 1:
        with open(sys.argv[1], "w", encoding="utf-8") as fh:
            json.dump(res, fh, indent=1)
    sys.exit(0 if (not failed and caught) else 1)
