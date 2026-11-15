"""EVEN-DICT decisive check. Pass lines are in ..\\logs\\00_passlines.log (fixed before the first run).

r side: exact rationals from the ball definitions (E_1, Dirichlet Green operator of the 3D Laplacian).
rho side: composite Gauss-Legendre quadrature of the hyperbolic forms and of the explicit kernel W.
Usage: python even_dict_check.py ..\\out\\even_dict_check.json
"""
import sys, json, time
from fractions import Fraction as Fr
from mpmath import mp, mpf, cosh, sinh, exp, log1p, expm1

mp.dps = 45
OUT = sys.argv[1] if len(sys.argv) > 1 else None
RES = {"lines": [], "plants": [], "tables": {}}
FAILS = []


def line(name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + name + ("  | " + info if info else ""), flush=True)
    RES["lines"].append({"name": name, "ok": bool(ok), "info": info})
    if not ok:
        FAILS.append(name)


def plant(name, caught, info=""):
    print(("CAUGHT " if caught else "MISSED ") + name + ("  | " + info if info else ""), flush=True)
    RES["plants"].append({"name": name, "caught": bool(caught), "info": info})
    if not caught:
        FAILS.append(name)


# ---------- exact r side: polynomials in r as {power: Fraction} ----------
def padd(a, b):
    c = dict(a)
    for k, v in b.items():
        c[k] = c.get(k, Fr(0)) + v
    return {k: v for k, v in c.items() if v != 0}


def pmul(a, b):
    c = {}
    for i, x in a.items():
        for j, y in b.items():
            c[i + j] = c.get(i + j, Fr(0)) + x * y
    return {k: v for k, v in c.items() if v != 0}


def pscale(a, s):
    return {k: v * s for k, v in a.items() if v * s != 0}


def pder(a):
    return {k - 1: v * k for k, v in a.items() if k != 0}


def pint01(a):
    assert all(k >= 0 for k in a), a
    return sum(v / (k + 1) for k, v in a.items())


NPOL = {0: Fr(1, 2), 2: Fr(-1, 2)}
R2 = {2: Fr(1)}
RM2 = {-2: Fr(1)}
RM1 = {-1: Fr(1)}


def data(l):
    mono = [{l + 2 * k: Fr(1)} for k in range(4)]
    fifth = {l: Fr(3), l + 2: Fr(-7, 2), l + 4: Fr(5, 3), l + 6: Fr(2, 7), l + 8: Fr(-11, 13)}
    return mono + [fifth]


def GD(p, l):
    """U with -Lap_3D(U Y) = p Y, U(1) = 0: Lap(r^(l+2k+2) Y) = (2k+2)(2k+2l+3) r^(l+2k) Y."""
    U = {}
    for pw, c in p.items():
        k2 = pw - l
        assert k2 >= 0 and k2 % 2 == 0
        d = (k2 + 2) * (k2 + 2 * l + 3)
        U = padd(U, {l: c / d, pw + 2: -c / d})
    return U


def rside(l):
    lam2 = l * (l + 1)
    D = data(l)
    n = len(D)
    A = [[Fr(1, 2) * pint01(pmul(NPOL, pmul(R2, pmul(D[i], D[j])))) for j in range(n)] for i in range(n)]
    B = [[Fr(1, 2) * pint01(pmul(NPOL, padd(pmul(R2, pmul(pder(D[i]), pder(D[j]))), pscale(pmul(D[i], D[j]), lam2))))
          + Fr(1, 2) * pint01(pmul(R2, pmul(D[i], D[j]))) for j in range(n)] for i in range(n)]
    # symplectic pairing int psidot_i psi_j: momentum datum i uses D[i], field datum j uses a shifted set so S is not symmetric
    F = [padd(D[j], pscale(D[(j + 1) % n], Fr(j + 1, 3))) for j in range(n)]
    S = [[pint01(pmul(R2, pmul(D[i], F[j]))) for j in range(n)] for i in range(n)]
    U = [GD(p, l) for p in D]
    okdef = True
    for p, u in zip(D, U):
        res = padd(padd(pder(pder(u)), pscale(pmul(RM1, pder(u)), 2)), pscale(pmul(RM2, u), -lam2))
        okdef &= (padd(res, p) == {}) and sum(u.values()) == 0
    Z = [[Fr(1, 2) * pint01(pmul(R2, pmul(D[i], U[j]))) for j in range(n)] for i in range(n)]
    return D, F, A, B, S, Z, okdef


def det(M):
    M = [r[:] for r in M]
    n = len(M)
    d = Fr(1)
    for c in range(n):
        piv = next((r for r in range(c, n) if M[r][c] != 0), None)
        if piv is None:
            return Fr(0)
        if piv != c:
            M[c], M[piv] = M[piv], M[c]
            d = -d
        d *= M[c][c]
        for r in range(c + 1, n):
            f = M[r][c] / M[c][c]
            for k in range(c, n):
                M[r][k] -= f * M[c][k]
    return d


# ---------- rho side ----------
def gl_nodes(n):
    xs, ws = [], []
    for i in range(1, n + 1):
        x = mp.cos(mp.pi * (i - mpf(1) / 4) / (n + mpf(1) / 2))
        for _ in range(100):
            p0, p1 = mpf(1), x
            for k in range(2, n + 1):
                p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
            dp = n * (x * p1 - p0) / (x * x - 1)
            dx = p1 / dp
            x -= dx
            if abs(dx) < mpf(10) ** (-mp.dps + 3):
                break
        p0, p1 = mpf(1), x
        for k in range(2, n + 1):
            p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
        dp = n * (x * p1 - p0) / (x * x - 1)
        xs.append(x)
        ws.append(2 / ((1 - x * x) * dp * dp))
    return xs, ws


PANELS = [mpf(x) for x in ("0 0.5 1 1.5 2 3 4 5 6 8 10 13 16 20 25 30 36 43 50 60 70 80 95 110".split())]


def geom(rho):
    """t = r = tanh(rho/2), 1 - t, N, sinh rho, all without cancellation."""
    x = 2 / (exp(rho) + 1)          # 1 - t
    t = 1 - x
    N = 1 / (1 + cosh(rho))
    return t, x, N, sinh(rho)


def peval(a, r):
    return sum(mpf(v.numerator) / mpf(v.denominator) * r ** k for k, v in a.items())


class Rho:
    def __init__(self, l, D, F, n, pi_pow=1, q_pow=0, neumann=False):
        self.l, self.n = l, n
        xs, ws = gl_nodes(n)
        self.xs, self.ws = xs, ws
        self.nodes = []  # per panel: list of (rho, weight)
        for a, b in zip(PANELS[:-1], PANELS[1:]):
            h, m = (b - a) / 2, (a + b) / 2
            self.nodes.append([(m + h * x, h * w) for x, w in zip(xs, ws)])
        self.D, self.F = D, F
        self.pi_pow, self.q_pow, self.neumann = pi_pow, q_pow, neumann
        self.dD = [pder(p) for p in D]
        self.dF = [pder(p) for p in F]

    # dictionary (with plant exponents): pi = N^pi_pow r p, q = N^q_pow r f
    def pi(self, i, rho):
        t, x, N, sh = geom(rho)
        return N ** self.pi_pow * t * peval(self.D[i], t)

    def q_and_dq(self, P, dP, rho):
        """q = N^q_pow r f and dq/drho, using dr/drho = N, dN/drho = -r N (chain rule only)."""
        t, x, N, sh = geom(rho)
        f, df = peval(P, t), peval(dP, t)
        g = N ** self.q_pow
        q = g * t * f
        dq = g * N * (f + t * df) + (self.q_pow * g * (-t)) * t * f
        return q, dq

    def a(self, rho):
        t = geom(rho)[0]
        return t ** (self.l + 1)

    def b(self, rho):
        t, x, N, sh = geom(rho)
        l = self.l
        if self.neumann:
            return (t ** (-l) + mpf(l) / (l + 1) * t ** (l + 1)) / (2 * l + 1)
        # t^-l - t^(l+1) = t^-l (1 - t^(2l+1)), with 1 - t^(2l+1) = -expm1((2l+1) log1p(-x))
        return t ** (-l) * (-expm1((2 * l + 1) * log1p(-x))) / (2 * l + 1)

    def int1(self, fn):
        return sum(w * fn(r) for pan in self.nodes for r, w in pan)

    def blocks_a(self):
        n = len(self.D)
        lam2 = self.l * (self.l + 1)
        pis = [[[self.pi(i, r) for r, w in pan] for pan in self.nodes] for i in range(n)]
        self.pis = pis
        W = [[w for r, w in pan] for pan in self.nodes]
        qD = [[[self.q_and_dq(self.D[i], self.dD[i], r) for r, w in pan] for pan in self.nodes] for i in range(n)]
        qF = [[[self.q_and_dq(self.F[i], self.dF[i], r) for r, w in pan] for pan in self.nodes] for i in range(n)]
        sh2 = [[sinh(r) ** 2 for r, w in pan] for pan in self.nodes]
        np_, nn = len(self.nodes), self.n
        A = [[sum(W[P][k] * pis[i][P][k] * pis[j][P][k] for P in range(np_) for k in range(nn)) / 2 for j in range(n)] for i in range(n)]
        B = [[sum(W[P][k] * (qD[i][P][k][1] * qD[j][P][k][1] + lam2 * qD[i][P][k][0] * qD[j][P][k][0] / sh2[P][k])
                  for P in range(np_) for k in range(nn)) / 2 for j in range(n)] for i in range(n)]
        S = [[sum(W[P][k] * pis[i][P][k] * qF[j][P][k][0] for P in range(np_) for k in range(nn)) for j in range(n)] for i in range(n)]
        return A, B, S

    def block_b(self):
        """(1/2)<pi_i, W pi_j> as a double integral: W = a(rho_<) b(rho_>); inner integral over rho' < rho by
        full earlier panels plus a remapped Gauss rule on [panel start, rho]; the rho' > rho half by symmetry."""
        n = len(self.D)
        np_, nn = len(self.nodes), self.n
        cum = [mpf(0)] * n               # int_0^{panel start} a pi_j
        Z = [[mpf(0)] * n for _ in range(n)]
        for P in range(np_):
            start = PANELS[P]
            full = [mpf(0)] * n
            for k in range(nn):
                rho, w = self.nodes[P][k]
                # inner: int_start^rho a(rho') pi_j(rho') drho'
                h, m = (rho - start) / 2, (rho + start) / 2
                inner = [mpf(0)] * n
                for x, wx in zip(self.xs, self.ws):
                    rp = m + h * x
                    ap = self.a(rp)
                    for j in range(n):
                        inner[j] += h * wx * ap * self.pi(j, rp)
                Aj = [cum[j] + inner[j] for j in range(n)]
                bk = self.b(rho)
                ak = self.a(rho)
                for i in range(n):
                    pik = self.pis[i][P][k]
                    full[i] += w * ak * pik
                    for j in range(n):
                        # int drho pi_i(rho) b(rho) int_{rho'<rho} a pi_j   (+ transpose added below)
                        Z[i][j] += w * bk * pik * Aj[j]
            cum = [cum[j] + full[j] for j in range(n)]
        return [[(Z[i][j] + Z[j][i]) / 2 for j in range(n)] for i in range(n)]


def relerr(Mr, Mx):
    """Mr numeric (mpf), Mx exact Fractions."""
    ex = [[mpf(v.numerator) / mpf(v.denominator) for v in row] for row in Mx]
    scale = max(abs(v) for row in ex for v in row)
    return max(abs(Mr[i][j] - ex[i][j]) for i in range(len(ex)) for j in range(len(ex[0]))) / scale


def fmt(x):
    return mp.nstr(x, 5)


t0 = time.time()
for l in (2, 3, 4):
    D, F, A, B, S, Z, okdef = rside(l)
    n = len(D)
    line(f"L7 l={l} r-side G_D solves the 3D Dirichlet problem for every datum (exact)", okdef)
    nus = [Fr(1, (j + 1) * (2 * l + 2 * j + 3)) for j in range(4)]
    ok8 = all(det([[Z[i][j] - nu * A[i][j] for j in range(4)] for i in range(4)]) == 0 for nu in nus) \
        and det([[Z[i][j] - Fr(1, 5) * A[i][j] for j in range(4)] for i in range(4)]) != 0
    line(f"L8 l={l} r-side spectrum of Z against (1/2)int N p^2 is 1/((j+1)(2l+2j+3)), j=0..3 (exact; 1/5 is not a root)", ok8)

    R = Rho(l, D, F, 60)
    Ar, Br, Sr = R.blocks_a()
    eA, eB, eS = relerr(Ar, A), relerr(Br, B), relerr(Sr, S)
    line(f"L1 l={l} (a1) momentum block, 5x5 Gram", eA <= mpf("1e-30"), f"rel err {fmt(eA)}; r-side [0][0]={A[0][0]}")
    line(f"L2 l={l} (a2) field block, 5x5 Gram (weight-1 term absorbed)", eB <= mpf("1e-30"), f"rel err {fmt(eB)}; r-side [0][0]={B[0][0]}")
    line(f"L3 l={l} (a3) symplectic pairing, 25 entries, non-symmetric", eS <= mpf("1e-30"), f"rel err {fmt(eS)}; r-side [0][0]={S[0][0]}")
    Zr = R.block_b()
    eZ = relerr(Zr, Z)
    line(f"L4 l={l} (b) (1/2)<psidot,G_D psidot> = (1/2)<pi,W pi>, g=1, 5x5 Gram, double integral", eZ <= mpf("1e-30"),
         f"rel err {fmt(eZ)}; r-side [0][0]={Z[0][0]}")
    gdev = max(abs(Zr[i][j] / (mpf(Z[i][j].numerator) / Z[i][j].denominator) - 1) for i in range(n) for j in range(n))
    line(f"L5 l={l} g read back entrywise: |Z_rho/Z_r - 1|", gdev <= mpf("1e-28"), f"max {fmt(gdev)}")
    R80 = Rho(l, D, F, 80)
    R80.blocks_a()
    Z80 = R80.block_b()
    sc = max(abs(v) for row in Z80 for v in row)
    e68 = max(abs(Zr[i][j] - Z80[i][j]) for i in range(n) for j in range(n)) / sc
    line(f"L6 l={l} quadrature self-check n=60 against n=80 on (b)", e68 <= mpf("1e-32"), f"rel diff {fmt(e68)}")

    P1 = Rho(l, D, F, 60, pi_pow=2)
    A1, _, S1 = P1.blocks_a()
    plant(f"P1 l={l} pi = N^2 r p (psidot = pi/(N^3 sinh rho)) breaks (a1) and (a3)",
          relerr(A1, A) > mpf("1e-3") and relerr(S1, S) > mpf("1e-3"), f"rel err {fmt(relerr(A1, A))}, {fmt(relerr(S1, S))}")
    P2 = Rho(l, D, F, 60, q_pow=1)
    _, B2, _ = P2.blocks_a()
    plant(f"P2 l={l} q = N r f breaks (a2)", relerr(B2, B) > mpf("1e-3"), f"rel err {fmt(relerr(B2, B))}")
    P3 = Rho(l, D, F, 60, neumann=True)
    P3.blocks_a()
    Z3 = P3.block_b()
    plant(f"P3 l={l} Neumann-type kernel breaks (b)", relerr(Z3, Z) > mpf("1e-3"), f"rel err {fmt(relerr(Z3, Z))}")
    RES["tables"][str(l)] = {
        "relerr": {"a1": fmt(eA), "a2": fmt(eB), "a3": fmt(eS), "b": fmt(eZ), "g_dev": fmt(gdev), "n60_vs_n80": fmt(e68)},
        "r_side_exact": {"A": [[str(v) for v in r] for r in A], "B": [[str(v) for v in r] for r in B],
                         "S": [[str(v) for v in r] for r in S], "Z": [[str(v) for v in r] for r in Z]},
        "rho_side_Z": [[mp.nstr(v, 32) for v in r] for r in Zr]}
    print(f"   l={l} done, {time.time()-t0:.0f} s", flush=True)

print("FAILED LINES: " + (", ".join(FAILS) if FAILS else "none"))
RES["failed"] = FAILS
if OUT:
    with open(OUT, "w") as fh:
        json.dump(RES, fh, indent=1)
sys.exit(1 if FAILS else 0)
