# BEYOND-NORMAL (Q11): the quasi-free 2 pi KMS state of the perturbed sector flow.
# Pass lines: ..\logs\01_prefixed_pass_lines.log (fixed before the first run).
# Usage: C:\Python314\python.exe -X utf8 beyond_normal.py ..\out\beyond_normal.json
import sys, json, time
import numpy as np

BETA = 2 * np.pi
HG, RMAX = 0.025, 50.0
N_GRID = int(round(RMAX / HG))
RHO = HG * np.arange(1, N_GRID + 1)
failed, results = [], {}


def line(name, ok, text):
    print(("PASS " if ok else "FAIL ") + name + ": " + text, flush=True)
    if not ok:
        failed.append(name + ": " + text)


def h2_matrix(l):
    lam2 = l * (l + 1)
    main = 2.0 / HG**2 + lam2 / np.sinh(RHO) ** 2
    main[-1] = 1.0 / HG**2 + lam2 / np.sinh(RHO[-1]) ** 2  # Neumann outer wall
    off = -np.ones(N_GRID - 1) / HG**2
    return np.diag(main) + np.diag(off, 1) + np.diag(off, -1)


def n_potential():
    return 1.0 / (1.0 + np.cosh(RHO))


def w_matrix(l):
    lr = np.log1p(-2.0 / (np.exp(RHO) + 1.0))  # log tanh(rho/2)
    lo = np.minimum.outer(lr, lr)
    hi = np.maximum.outer(lr, lr)
    return np.exp((l + 1) * lo - l * hi) * (-np.expm1((2 * l + 1) * hi)) / (2 * l + 1) * HG


def sym_fun(mat, f):
    w, v = np.linalg.eigh(mat)
    return (v * f(w)) @ v.T


class Sector:
    """Stable quadratic flow H = 1/2<pi,(1+Wm)pi> + 1/2<q,h2 q>, normal coordinates q'=A^-1 q, pi'=A pi."""

    def __init__(self, h2, wm=None):
        n = h2.shape[0]
        if wm is None:
            self.A = self.Ai = None
            om2 = h2
        else:
            ww, wv = np.linalg.eigh(np.eye(n) + wm)
            self.A = (wv * np.sqrt(ww)) @ wv.T
            self.Ai = (wv / np.sqrt(ww)) @ wv.T
            om2 = self.A @ h2 @ self.A
            om2 = 0.5 * (om2 + om2.T)
        self.e2, self.V = np.linalg.eigh(om2)
        self.nu = np.sqrt(np.abs(self.e2))

    def _a(self, x):
        return x if self.A is None else self.A @ x

    def _ai(self, x):
        return x if self.Ai is None else self.Ai @ x

    def forms(self, beta):
        """(Mq, Mp): mu(k) = 1/2 q.Mq q + 1/2 pi.Mp pi for the smearing k = (q, pi)."""
        c = 1.0 / np.tanh(beta * self.nu / 2)
        mq = (self.V * (self.nu * c)) @ self.V.T
        mp = (self.V * (c / self.nu)) @ self.V.T
        if self.A is not None:
            mq = self.Ai @ mq @ self.Ai
            mp = self.A @ mp @ self.A
        return 0.5 * (mq + mq.T), 0.5 * (mp + mp.T)

    def flow(self, q, p, tau):
        a = self.V.T @ self._ai(q)
        b = self.V.T @ self._a(p)
        c, s = np.cos(self.nu * tau), np.sin(self.nu * tau)
        return self._a(self.V @ (c * a + s / self.nu * b)), self._ai(self.V @ (-self.nu * s * a + c * b))


def mu(forms, q, p):
    return 0.5 * q @ forms[0] @ q + 0.5 * p @ forms[1] @ p


def gen_eigs(m0, mg):
    w, v = np.linalg.eigh(m0)
    r = (v / np.sqrt(w)) @ v.T
    x = r @ mg @ r
    return np.linalg.eigvalsh(0.5 * (x + x.T))


def one_minus_affinity(vac, pert):
    f0, fg = vac.forms(BETA / 2), pert.forms(BETA / 2)
    tot = 0.0
    for s in (0, 1):
        r = gen_eigs(f0[s], fg[s])
        tot += np.sum(np.log1p((np.sqrt(r) - 1.0) ** 2 / (2.0 * np.sqrt(r))))
    return -np.expm1(-0.5 * tot)


def r_curve(vac_forms, target_forms, pert, q, p, taus):
    m0 = mu(vac_forms, q, p)
    den = mu(target_forms, q, p) - m0
    return [(mu(vac_forms, *pert.flow(q, p, t)) - m0) / den for t in taus], den / m0


HINGE = {2: 0.0251, 16: 0.0460, 64: 0.0510}
TAUS = [5.0, 10.0, 15.0, 20.0, 30.0]


def run_sector(l):
    lam = np.sqrt(l * (l + 1))
    out = {}
    h2 = h2_matrix(l)
    npot = n_potential()
    vac = Sector(h2)
    f0 = vac.forms(BETA)
    lo0 = vac.e2[0]
    rc = np.log(2 * lam) + 1.0
    prof = np.tanh(RHO) ** 3 * (RHO - rc) * np.exp(-0.5 * (RHO - rc) ** 2)
    zero = np.zeros(N_GRID)

    # ---- odd ----
    odd = Sector(h2 + np.diag(npot))
    oddm = Sector(h2 - np.diag(npot))
    line(f"L1 odd l={l}", odd.e2[0] >= lo0 * (1 - 1e-9) > 0,
         f"lambda_min(h^2+N) = {odd.e2[0]:.6e} >= lambda_min(h^2) = {lo0:.6e} > 0")
    line(f"L1b odd g=-1 l={l}", oddm.e2[0] > 0, f"lambda_min(h^2-N) = {oddm.e2[0]:.6e}")
    hinv = (vac.V / vac.nu) @ vac.V.T
    cl = np.linalg.eigvalsh(hinv @ np.diag(npot) @ hinv)[-1]
    line(f"L2o l={l}", cl <= 3.52, f"C_l = lambda_max(h^-1 N h^-1) = {cl:.6f} (lambda*C_l = {lam*cl:.4f})")
    out["odd"] = {"lam_min": float(odd.e2[0]), "lam_min_vac": float(lo0), "C_l": float(cl)}
    if l != 64:
        fg = odd.forms(BETA)
        rq, rp = gen_eigs(f0[0], fg[0]), gen_eigs(f0[1], fg[1])
        line(f"L3o l={l}", rq[0] >= 1 - 1e-9 and rq[-1] <= 4.52 and rp[0] >= 1 / 4.52 and rp[-1] <= 1 + 1e-9,
             f"q-slot ratio in [{rq[0]:.9f}, {rq[-1]:.6f}], pi-slot in [{rp[0]:.6f}, {rp[-1]:.9f}]")
        out["odd"].update(q_slot=[float(rq[0]), float(rq[-1])], p_slot=[float(rp[0]), float(rp[-1])])
        sector_tests(l, "odd", vac, f0, odd, fg, prof, zero, out["odd"])
    ap, am = one_minus_affinity(vac, odd), one_minus_affinity(vac, oddm)
    val = l * l * 0.5 * (ap + am)
    line(f"L6 odd l={l}" + (" (extra)" if l == 64 else ""), abs(val / HINGE[l] - 1) <= 0.04,
         f"l^2 mean(1-A) = {val:.5f} vs HINGE {HINGE[l]} (ratio {val/HINGE[l]:.4f}); "
         f"g=+1: l^2(1-A) = {l*l*ap:.5f}, lambda^2(1-A) = {lam*lam*ap:.5f}; g=-1: {l*l*am:.5f}")
    out["odd"].update(one_minus_A_plus=float(ap), one_minus_A_minus=float(am), l2_mean=float(val))
    if l == 64:
        return out

    # ---- even ----
    wm = w_matrix(l)
    wnorm = np.linalg.eigvalsh(0.5 * (wm + wm.T))[-1]
    line(f"L2e l={l}", abs(wnorm * (2 * l + 3) - 1) <= 0.01, f"||W|| = {wnorm:.8f}, (2l+3)||W|| = {wnorm*(2*l+3):.6f}")
    ev = Sector(h2, wm)
    line(f"L1 even l={l}", ev.e2[0] >= lo0 * (1 - 1e-9) > 0,
         f"lambda_min(A h^2 A) = {ev.e2[0]:.6e} >= lambda_min(h^2) = {lo0:.6e} > 0")
    fg = ev.forms(BETA)
    rq, rp = gen_eigs(f0[0], fg[0]), gen_eigs(f0[1], fg[1])
    w = 1.0 / (2 * l + 3)
    lo, hi = 1 / (np.pi * (1 + w)), np.pi * (1 + w)
    line(f"L3e l={l}", min(rq[0], rp[0]) >= lo and max(rq[-1], rp[-1]) <= hi,
         f"q-slot ratio in [{rq[0]:.6f}, {rq[-1]:.6f}], pi-slot in [{rp[0]:.6f}, {rp[-1]:.6f}]; 1+w = {1+w:.6f}")
    out["even"] = {"lam_min": float(ev.e2[0]), "W_norm": float(wnorm),
                   "q_slot": [float(rq[0]), float(rq[-1])], "p_slot": [float(rp[0]), float(rp[-1])]}
    sector_tests(l, "even", vac, f0, ev, fg, prof, zero, out["even"])
    ae = one_minus_affinity(vac, ev)
    line(f"L6 even l={l}", ae > 0, f"1-A = {ae:.6e}; lambda^2(1-A) = {lam*lam*ae:.5f}, l^2(1-A) = {l*l*ae:.5f} (reported)")
    out["even"]["one_minus_A"] = float(ae)
    return out


def sector_tests(l, par, vac, f0, pert, fg, prof, zero, out):
    q, p = prof, 0.7 * np.roll(prof, 40)
    qt, pt = pert.flow(q, p, 1.7)
    ig = abs(mu(fg, qt, pt) / mu(fg, q, p) - 1)
    i0 = abs(mu(f0, qt, pt) / mu(f0, q, p) - 1)
    line(f"L4 {par} l={l}", ig <= 1e-9 and i0 >= 1e-8, f"mu_g drift {ig:.2e} (<=1e-9); mu_0 drift {i0:.3e} (>=1e-8)")
    out["L4"] = [float(ig), float(i0)]
    for name, (qq, pp) in (("q-data", (prof, zero)), ("pi-data", (zero, prof))):
        rr, rel = r_curve(f0, fg, pert, qq, pp, TAUS)
        line(f"L5 {par} l={l} {name}", abs(rr[3] - 1) <= 0.10,
             "r(5,10,15,20,30) = " + ", ".join(f"{x:.4f}" for x in rr) + f"; (mu_g-mu_0)/mu_0 = {rel:+.4e}")
        out["L5 " + name] = [float(x) for x in rr]
    if l == 2 and par == "odd":
        rr, _ = r_curve(f0, pert.forms(np.pi), pert, prof, zero, TAUS)
        ok = abs(rr[3] - 1) >= 0.5
        print(("CAUGHT " if ok else "MISSED ") + f"P2 plant (target at beta = pi): r(20) = {rr[3]:.5f}", flush=True)
        results["P2_caught"] = bool(ok)
        results["P2_r20"] = float(rr[3])


def toy():
    rng = np.random.default_rng(73)
    n = 6
    qm, _ = np.linalg.qr(rng.standard_normal((n, n)))
    h2 = (qm * np.linspace(0.3, 1.3, n) ** 2) @ qm.T
    b = rng.standard_normal((n, n))
    npos = 0.15 * b @ b.T / n
    c = rng.standard_normal((n, n))
    wpos = 0.2 * c @ c.T / n
    jm = np.block([[np.zeros((n, n)), np.eye(n)], [-np.eye(n), np.zeros((n, n))]])

    def mismatch(h2p, wm, beta_build):
        s = Sector(h2p, wm)
        mq, mp = s.forms(beta_build)
        vop = np.block([[0.5 * mp, np.zeros((n, n))], [np.zeros((n, n)), 0.5 * mq]])
        k = -2j * vop @ np.linalg.inv(jm)
        ew, ev = np.linalg.eig(k)
        xrec = (2.0 / (1j * BETA)) * (ev * (0.5 * np.log((ew + 1) / (ew - 1)))) @ np.linalg.inv(ev)
        kin = np.eye(n) if wm is None else np.eye(n) + wm
        x = jm @ np.block([[h2p, np.zeros((n, n))], [np.zeros((n, n)), kin]])
        return np.linalg.norm(xrec - x) / np.linalg.norm(x)

    res = {}
    for name, (hh, ww) in (("vacuum", (h2, None)), ("odd-type", (h2 + npos, None)), ("even-type", (h2, wpos))):
        m = mismatch(hh, ww, BETA)
        line(f"L7 toy {name}", m <= 1e-8, f"|X_rec - JG|/|JG| = {m:.2e}")
        res[name] = float(m)
    m = mismatch(h2, wpos, np.pi)
    ok = m >= 0.2
    print(("CAUGHT " if ok else "MISSED ") + f"P1 plant (covariance at beta = pi, tested at 2 pi): mismatch {m:.4f}", flush=True)
    results["P1_caught"] = bool(ok)
    results["P1_mismatch"] = float(m)
    results["toy"] = res


if __name__ == "__main__":
    t0 = time.time()
    toy()
    for l in (2, 16, 64):
        results[f"l={l}"] = run_sector(l)
        print(f"  [l={l} done, {time.time()-t0:.0f} s]", flush=True)
    print("FAILED LINES: " + ("none" if not failed else ""))
    for f in failed:
        print("  " + f)
    plants = results.get("P1_caught", False) and results.get("P2_caught", False)
    print(f"plants caught: P1 {results.get('P1_caught')}, P2 {results.get('P2_caught')}")
    results["failed"] = failed
    if len(sys.argv) > 1:
        with open(sys.argv[1], "w", encoding="utf-8") as fh:
            json.dump(results, fh, indent=1)
    sys.exit(0 if (not failed and plants) else 1)
