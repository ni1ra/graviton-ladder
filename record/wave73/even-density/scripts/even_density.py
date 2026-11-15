"""EVEN-DENSITY decisive check: is <pi, h^-2 pi> finite for every C_c^inf pi?
h^2 = -d^2/drho^2 + l(l+1)/sinh^2(rho) on (0, inf).  Pass lines: logs/01_prereg.md.
Usage: python even_density.py <out.json>
"""
import sys, json
import numpy as np
from scipy.linalg import solve_banded
from scipy.integrate import cumulative_simpson, simpson

D_FD = 0.004
BOXES = [8, 16, 32, 64, 128, 256, 512, 1024]   # 1024 only to form D(512)
FAILS = []


def line(name, ok, msg):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {msg}")
    if not ok:
        FAILS.append(f"{name}: {msg}")


def bump(x, a, b):
    y = np.zeros_like(x)
    m = (x > a) & (x < b)
    y[m] = np.exp(-1.0 / ((x[m] - a) * (b - x[m])))
    return y


def legendre_P(l, x):
    p0, p1 = np.ones_like(x), x.copy()
    if l == 0:
        return p0
    for n in range(1, l):
        p0, p1 = p1, ((2 * n + 1) * x * p1 - n * p0) / (n + 1)
    return p1


class Zero:
    """Zero-energy solutions: chi = P_l(coth rho) (-> 1), phi = chi * int_0^rho chi^-2 (regular at 0)."""

    def __init__(self, l, R=30.0, d=5e-4):
        self.l, self.R = l, R
        self.x = np.arange(1, int(round(R / d)) + 1) * d
        self.chi = legendre_P(l, 1.0 / np.tanh(self.x))
        self.J = cumulative_simpson(self.chi ** -2.0, x=self.x, initial=0.0)
        self.phi = self.chi * self.J
        self.JR = self.J[-1]

    def chi_over_phi(self, L):
        if L <= self.R:
            return 1.0 / np.interp(L, self.x, self.J)
        return 1.0 / (self.JR + (L - self.R))      # chi = 1 to double precision beyond R

    def phi_at(self, r):
        return np.interp(r, self.x, self.phi) if r <= self.R else self.JR + (r - self.R)

    def chi_at(self, r):
        return np.interp(r, self.x, self.chi) if r <= self.R else 1.0


def pis(z):
    x = z.x
    b13, b46 = bump(x, 1, 3), bump(x, 4, 6)
    c = simpson(z.phi * b13, x=x) / simpson(z.phi * b46, x=x)
    f = {
        "pi1": lambda r: bump(r, 1, 3),
        "pi2": lambda r: bump(r, 1, 3) - c * bump(r, 4, 6),
        "pi3": lambda r: bump(r, 0.5, 4.5) * np.cos(4 * r),
    }
    return f, c


def green_values(z, f):
    x = z.x
    p = f(x)
    inner = cumulative_simpson(z.phi * p, x=x, initial=0.0)
    I_inf = 2.0 * simpson(p * z.chi * inner, x=x)
    ov = simpson(z.phi * p, x=x)
    return I_inf, ov


def fd_solve(l, L, rhs_funcs, plant=False):
    n = int(round(L / D_FD))
    if plant:   # no potential, Neumann at 0, Dirichlet at L
        r = np.arange(0, n) * D_FD
        V = np.zeros_like(r)
    else:       # Dirichlet at 0 and L
        r = np.arange(1, n) * D_FD
        V = l * (l + 1) / np.sinh(np.minimum(r, 300.0)) ** 2
    m = len(r)
    ab = np.zeros((3, m))
    ab[1] = 2.0 / D_FD ** 2 + V
    ab[0, 1:] = -1.0 / D_FD ** 2
    ab[2, :-1] = -1.0 / D_FD ** 2
    if plant:
        ab[0, 1] = -2.0 / D_FD ** 2      # ghost node: w(-d) = w(d)
    B = np.stack([f(r) for f in rhs_funcs], axis=1)
    Wsol = solve_banded((1, 1), ab, B)
    I = D_FD * np.sum(B * Wsol, axis=0)
    return r, B, Wsol, I


def detector(Ds, scale):
    ok = True
    for a, b in zip(Ds[:-1], Ds[1:]):
        small = abs(a) <= 1e-6 * abs(scale) and abs(b) <= 1e-6 * abs(scale)
        if not (small or (b != 0 and a / b >= 1.5)):
            ok = False
    return "CONVERGES" if ok else "DIVERGES"


def main(out):
    res = {}
    # P1: construction
    print("P1 construction")
    z1 = Zero(1)
    m = (z1.x >= 0.5) & (z1.x <= 6)
    err = np.max(np.abs(z1.phi[m] - (z1.x[m] / np.tanh(z1.x[m]) - 1.0)))
    line("P1 l=1 closed form", err <= 1e-6, f"max|phi - (rho coth rho - 1)| = {err:.2e}")
    for l in (2, 16):
        z = Zero(l)
        dphi, dchi = np.gradient(z.phi, z.x), np.gradient(z.chi, z.x)
        for r0 in (1.0, 2.0, 4.0):
            i = int(round(r0 / 5e-4)) - 1
            wr = dphi[i] * z.chi[i] - z.phi[i] * dchi[i]
            line(f"P1 l={l} Wronskian rho={r0}", abs(wr - 1) <= 1e-5, f"{wr:.8f}")
        cl = 40.0 - z.phi_at(40.0)
        print(f"     l={l}: phi ~ rho - c_l with c_l = {cl:.6f}")
        res[f"c_l={l}"] = cl

    for l in (2, 16):
        lam2 = l * (l + 1)
        z = Zero(l)
        f, c = pis(z)
        names = list(f)
        G = {k: green_values(z, f[k]) for k in names}
        fd = {}
        w_mid = {}
        for L in BOXES:
            r, B, Wsol, I = fd_solve(l, L, [f[k] for k in names])
            fd[L] = dict(zip(names, I))
            if L == 512:
                j = int(round(256 / D_FD)) - 1
                for q, k in enumerate(names):
                    w_mid[k] = (Wsol[j, q], np.max(np.abs(Wsol[:, q])))
        print(f"\n=== l = {l}  (c for pi2 = {c:.6f}) ===")
        for k in names:
            I_inf, ov = G[k]
            print(f"\n {k}: I_inf (Green) = {I_inf:.9e}   <phi,pi> = {ov:.6e}")
            print("     L      I_L (FD)         I_L (Green law)   rel diff    D(L)=I_2L-I_L   D/pred    D(L)/D(2L)")
            Ds, rows = [], []
            for L in BOXES[:-1]:
                pred = I_inf - ov ** 2 * z.chi_over_phi(L)
                rel = abs(fd[L][k] - pred) / abs(I_inf)
                D = fd[2 * L][k] - fd[L][k]
                Dp = ov ** 2 * (z.chi_over_phi(L) - z.chi_over_phi(2 * L))
                Ds.append(D)
                rows.append((L, fd[L][k], pred, rel, D, Dp))
            for i, (L, a, pred, rel, D, Dp) in enumerate(rows):
                rat = Ds[i] / Ds[i + 1] if i + 1 < len(Ds) and Ds[i + 1] != 0 else float("nan")
                dp = D / Dp if abs(Dp) > 1e-300 else float("nan")
                print(f"  {L:5d}  {a:.9e}  {pred:.9e}  {rel:.1e}  {D: .6e}  {dp:8.5f}  {rat:8.4f}")
            relmax = max(rw[3] for rw in rows)
            line(f"P2 l={l} {k}", relmax <= 1e-3, f"max rel diff FD vs Green law = {relmax:.2e}")
            verdict = detector(Ds[1:], fd[1024][k])
            if k in ("pi1", "pi3"):
                rats = [Ds[i] / Ds[i + 1] for i in range(1, len(Ds) - 1)]
                dps = [rw[4] / rw[5] for rw in rows]
                line(f"P3 l={l} {k} ratios", all(1.8 <= x <= 2.6 for x in rats),
                     "D(L)/D(2L), L=16..256: " + ", ".join(f"{x:.4f}" for x in rats))
                line(f"P3 l={l} {k} law", all(0.99 <= x <= 1.01 for x in dps),
                     f"D/pred in [{min(dps):.5f}, {max(dps):.5f}]")
                line(f"P3 l={l} {k} verdict", verdict == "CONVERGES", verdict)
            else:
                mx = max(abs(d) for d in Ds)
                line(f"P4 l={l} {k}", mx <= 1e-6 * abs(I_inf), f"max|D(L)| = {mx:.2e}, I_inf = {I_inf:.4e}")
                line(f"P4 l={l} {k} verdict", verdict == "CONVERGES", verdict)
            p = f[k](z.x)
            b1 = simpson((p * np.sinh(z.x)) ** 2, x=z.x) / lam2
            b2 = 4 * simpson((p * z.x) ** 2, x=z.x)
            line(f"P5 l={l} {k}", 0 < I_inf <= b1 and I_inf <= b2,
                 f"I_inf = {I_inf:.4e} <= {b1:.4e} (sinh bound), <= {b2:.4e} (Hardy bound)")
            res[f"l={l},{k}"] = {"I_inf": I_inf, "overlap": ov, "FD": {str(L): fd[L][k] for L in BOXES},
                                 "verdict": verdict}
        ov1 = G["pi1"][1]
        pw = ov1 * (z.chi_at(256) - z.chi_over_phi(512) * z.phi_at(256))
        line(f"P6 l={l} pi1 plateau", abs(w_mid["pi1"][0] / pw - 1) <= 1e-2,
             f"w_FD(256) = {w_mid['pi1'][0]:.6e}, predicted {pw:.6e} (= <phi,pi>/2 scale: w is not in L^2)")
        line(f"P6 l={l} pi2 compact", abs(w_mid["pi2"][0]) <= 1e-6 * w_mid["pi2"][1],
             f"|w_FD(256)|/max|w| = {abs(w_mid['pi2'][0]) / w_mid['pi2'][1]:.2e}")

    # PLANT
    print("\n=== PLANT: no potential, Neumann at 0 (zero-energy resonance) ===")
    f1 = lambda r: bump(r, 1, 3)
    xs = np.arange(0, 4, 5e-4)
    m0 = simpson(f1(xs), x=xs)
    Ip = {L: fd_solve(0, L, [f1], plant=True)[3][0] for L in BOXES}
    Ds = [Ip[2 * L] - Ip[L] for L in BOXES[:-1]]
    print("     L      I_L (plant)     D(L)          D(L)/(L (int pi)^2)")
    for L, D in zip(BOXES[:-1], Ds):
        print(f"  {L:5d}  {Ip[L]:.6e}  {D:.6e}  {D / (L * m0 ** 2):.5f}")
    verdict = detector(Ds[1:], Ip[1024])
    okp = verdict == "DIVERGES" and all(0.99 <= D / (L * m0 ** 2) <= 1.01 for L, D in zip(BOXES[:-1], Ds))
    line("PLANT caught", okp, f"detector verdict {verdict}; I_L grows linearly")
    res["plant"] = {"I": {str(L): Ip[L] for L in BOXES}, "verdict": verdict}

    res["fails"] = FAILS
    with open(out, "w") as fh:
        json.dump(res, fh, indent=1)
    print("\nALL LINES HELD" if not FAILS else f"\nFAILED LINES ({len(FAILS)}):\n  " + "\n  ".join(FAILS))
    return 0 if not FAILS else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
