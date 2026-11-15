"""Lane RULER-L (wave 73): RULER's computation at l = 2..6, and the one-dimensional inequality.

Part 1 (per l, both parities, exact): Delta = Q_2 - Q_3 with null-cone's forms
Q_p = 2 oint int_0^1 v^2 f |sigma_p|^2, f = v(1 - v), as a bilinear form on Cauchy data at t = 0
(TT gauge), Gram matrices in the chain basis psi = r^(2j) phi, psi_dot = r^(2j) phi (j = 0..3),
divided by N_Y = oint phi^2.  Adapted from ruler\\scripts\\ruler.py; nc_core.py is an unedited copy.

Part 2 (one generator, exact): with h = v^2 (sigma_2 - sigma_3),
    Q_2 - Q_3 = (norm) int_0^1 h^2 / v dv = (norm) int_0^1 v^3 (sigma_2 - sigma_3)^2 dv  >= 0.
Checked on monomials, on polynomial pencils of degree <= d, on random polynomials and on the
family h' = v^2 (1 - v)^N that sends Q_2/Q_3 to infinity.

Usage: python ruler_l.py OUT.json [--l 2,3,4,5,6] [--plant swap] [--no-part2]
"""
import json
import math
import random
import sys
import time
from fractions import Fraction

import sympy as sp

import nc_core as nc
from nc_core import t, v, x, y, z, XS

plant = sys.argv[sys.argv.index("--plant") + 1] if "--plant" in sys.argv else None
LARG = sys.argv[sys.argv.index("--l") + 1] if "--l" in sys.argv else "2,3,4,5,6"
LS = [] if LARG == "none" else [int(a) for a in LARG.split(",")]
NB = 4
P_CHARGE, P_ENTROPY = (3, 2) if plant == "swap" else (2, 3)     # the plant exchanges the two moments
t0 = time.time()
res = {"script": "ruler_l.py", "plant": plant, "checks": [], "per_l": {}}
ok = True


def check(name, cond, detail=""):
    global ok
    cond = bool(cond)
    ok &= cond
    res["checks"].append({"name": name, "pass": cond, "detail": str(detail)})
    print(("PASS " if cond else "FAIL ") + name + (f" | {detail}" if detail != "" else ""), flush=True)


def exact(e):
    """An exact rational number, by cancellation only.  (RULER used sp.nsimplify here; on rationals with
    long numerators it can return a nearby algebraic guess, which run 03 of this lane hit at l = 3, 4.)"""
    e = sp.cancel(e)
    assert e.is_Rational, f"not a rational number: {e}"
    return e


def definiteness(G):
    """Exact verdict from leading principal minors (Sylvester); float eigenvalues only for the record."""
    G = sp.Matrix(G)
    minors = [G[:k, :k].det() for k in range(1, G.shape[0] + 1)]
    try:
        lam_ = sp.Symbol("lam_")
        ev = sorted(float(sp.re(r)) for r in sp.Poly(G.charpoly(lam_).as_expr(), lam_).nroots(n=40, maxsteps=2000))
    except Exception as exc:                                   # the verdict never depends on floats
        ev = [f"float eigenvalues not computed: {type(exc).__name__}"]
    if all(m > 0 for m in minors):
        verdict = "positive definite"
    elif all((-1) ** (k + 1) * m > 0 for k, m in enumerate(minors)):
        verdict = "negative definite"
    else:
        verdict = "not definite"
    return verdict, minors, ev


def rel_spectrum(Q2, S):
    lam_ = sp.Symbol("lam_")
    try:
        pol = sp.Poly((sp.Matrix(Q2) - lam_ * sp.Matrix(S)).det(), lam_)
        return sorted(float(sp.re(r)) for r in pol.nroots(n=40, maxsteps=2000))
    except Exception as exc:
        return [f"not computed: {type(exc).__name__}"]


def mstr(G):
    return [[str(e) for e in G.row(i)] for i in range(G.shape[0])]


def fmt(xs):
    return ["%.6g" % e if isinstance(e, float) else e for e in xs]


# ====================================================================== part 1
def run_l(L):
    tl = time.time()
    print(f"\n==================== l = {L} ====================", flush=True)
    phi = nc.harmonic(L)
    N_Y = nc.surf_int(phi**2)
    s = x**2 + y**2 + z**2
    gphi = nc.grad(phi)
    Hphi = nc.hess(phi)
    W = v**2 * nc.F_CONE
    target = sp.Rational(L + 2, L + 1) ** 2

    def even_tt(k, deg=5):
        """Even-type TT tensor with x.T.x = s^k phi."""
        cs = sp.symbols(f"c0:{4 * (deg + 1)}")
        f = [sum(cs[i * (deg + 1) + m] * s**m for m in range(deg + 1)) for i in range(4)]
        T = sp.Matrix(3, 3, lambda i, j: f[0] * Hphi[i, j] + f[1] * (XS[i] * gphi[j] + XS[j] * gphi[i])
                      + (f[2] * phi if i == j else 0) + f[3] * phi * XS[i] * XS[j])
        eqs = [sp.expand(T.trace())]
        eqs += [sp.expand(sum(sp.diff(T[i, j], XS[i]) for i in range(3))) for j in range(3)]
        eqs += [sp.expand(sum(XS[i] * XS[j] * T[i, j] for i in range(3) for j in range(3)) - s**k * phi)]
        lin = []
        for e in eqs:
            lin += [c for c in sp.Poly(e, x, y, z).coeffs()]
        sol = sp.solve(lin, cs, dict=True)
        assert sol, "no even TT tensor found"
        sol = sol[0]
        full = {c: sol.get(c, c) for c in cs}
        free = set().union(*[sp.sympify(e).free_symbols for e in full.values()]) & set(cs)
        T = T.subs(full).subs({c: 0 for c in free}).applyfunc(sp.expand)
        return T, len(free)

    V = nc.vec_potential(phi)
    K0 = nc.sym_grad(V)

    def odd_tt(k):
        if k == 0:
            return K0
        a = -sp.Rational(2 * (k - 1) + L + 4, (L - 1) * k)
        return sp.Matrix(3, 3, lambda i, j: sp.expand(a * s**k * K0[i, j] + s ** (k - 1) * (XS[i] * V[j] + XS[j] * V[i])))

    def profile(e):
        e = sp.expand(e)
        if e == 0:
            return [sp.Integer(0)] * NB
        p = sp.cancel(e / phi)
        p1 = sp.Poly(sp.expand(p.subs({y: 0, z: 0})), x)
        co = [p1.coeff_monomial(x ** (2 * j)) for j in range(NB)]
        assert sp.expand(e - sum(co[j] * s**j for j in range(NB)) * phi) == 0, "chain field is not a profile times phi"
        return co

    def lowest_power(C):
        return min(min(m[0] for m in sp.Poly(e, v).monoms()) for e in C if e != 0)

    def build(parity):
        data = []
        flags = {"tt": True, "vac": True, "other": True, "unique": True}
        lows = []
        for kind in ("metric", "extrinsic"):
            for i in range(NB):
                if parity == "even":
                    H, nfree = even_tt(i + 1)
                    flags["unique"] &= nfree == 0
                else:
                    H = odd_tt(i + 1 if kind == "metric" else i)
                flags["tt"] &= nc.is_tt(H)
                h = nc.tt_series(H, 0 if kind == "metric" else 1)
                cv = nc.Curv(h)
                flags["vac"] &= cv.ricci_zero()
                Wy = cv.E() if parity == "even" else cv.B()
                Wo = cv.B() if parity == "even" else cv.E()
                chain = sp.expand(sum(XS[a] * XS[b] * Wy[a, b] for a in range(3) for b in range(3)))
                other = sp.expand(sum(XS[a] * XS[b] * Wo[a, b] for a in range(3) for b in range(3)))
                flags["other"] &= other == 0
                row = profile(chain.subs(t, 0)) + profile(sp.diff(chain, t).subs(t, 0))
                C = cv.cone_weyl()
                lows.append(lowest_power(C))
                data.append((kind, C, row, cv))
            print(f"   l={L} {parity} {kind} data built at {time.time() - tl:.0f} s", flush=True)
        if parity == "even":
            check(f"l={L} even: each TT tensor with x.T.x = s^k phi is unique (8 data)", flags["unique"])
        check(f"l={L} {parity}: all 8 data transverse and traceless", flags["tt"])
        check(f"l={L} {parity}: all 8 data solve the linearised vacuum equations exactly", flags["vac"])
        check(f"l={L} {parity}: the other parity's chain field vanishes on all 8 data", flags["other"])
        check(f"l={L} {parity}: the cone Weyl component starts at v^(l-2) or higher on all 8 data",
              min(lows) >= L - 2, f"lowest powers {lows}")
        return data

    def gram(A, B, weight):
        n = len(A)
        G = sp.zeros(n, n)
        for a in range(n):
            for b in range(a, n):
                G[a, b] = G[b, a] = exact(2 * nc.cone_int(nc.norm2(A[a], B[b]), weight) / N_Y)
        return G

    def e1_form():
        cl = 4 * sp.factorial(L - 2) / sp.factorial(L + 2)
        N = (1 - s) / 2
        G = sp.zeros(2 * NB, 2 * NB)
        for a in range(NB):
            for b in range(NB):
                pa, pb = s**a * phi, s**b * phi
                G[a, b] = cl * (nc.ball_int(N * sum(sp.diff(pa, c) * sp.diff(pb, c) for c in XS)) / 2
                                + nc.ball_int(pa * pb) / 2) / N_Y
                G[NB + a, NB + b] = cl * nc.ball_int(N * pa * pb) / 2 / N_Y
        return G.applyfunc(exact)

    S_FORM = e1_form()
    out = {}
    keep = {}
    for parity in ("even", "odd"):
        data = build(parity)
        mq = [nc.moment(d[1], P_CHARGE) for d in data]       # sigma_2 unless planted
        ms = [nc.moment(d[1], P_ENTROPY) for d in data]      # sigma_3 unless planted
        M = sp.Matrix([d[2] for d in data])
        check(f"l={L} {parity}: the eight data span the chain basis (psi, psi_dot; j = 0..3)", M.det() != 0)
        Mi = M.inv()
        G2 = (Mi * gram(mq, mq, W) * Mi.T).applyfunc(exact)
        G3 = (Mi * gram(ms, ms, W) * Mi.T).applyfunc(exact)
        md = [a - b for a, b in zip(mq, ms)]
        GH = (Mi * gram(md, md, v**3) * Mi.T).applyfunc(exact)
        print(f"   l={L} {parity} Gram matrices at {time.time() - tl:.0f} s", flush=True)
        kinds = ["metric" if any(M[a, j] != 0 for a in range(NB)) else "extrinsic" for j in (0, NB)]
        for a in range(2 * NB):
            lo = all(M[a, j] == 0 for j in range(NB))
            hi = all(M[a, j] == 0 for j in range(NB, 2 * NB))
            assert lo != hi, "a datum mixes psi and psi_dot"
        check(f"l={L} {parity}: Q_3 Gram equals c_l E_1 on the chain basis (block diagonal)", (G3 - S_FORM).is_zero_matrix)
        check(f"l={L} {parity}: Q_3 of the first datum equals nc.entropy_slice",
              sp.simplify(2 * nc.cone_int(nc.norm2(ms[0]), W) - nc.entropy_slice(data[0][3], L)) == 0)
        D = G2 - G3
        check(f"l={L} {parity}: identity Delta = 2 oint int v^3 |sigma_2 - sigma_3|^2 holds on the whole 8x8 Gram",
              (D - GH).is_zero_matrix)
        o = {"psi_is": kinds[0], "psi_dot_is": kinds[1], "Q2": mstr(G2), "S": mstr(G3), "Delta": {}}
        for bn, B in {"psi": D[:NB, :NB], "psi_dot": D[NB:, NB:], "full": D}.items():
            verdict, minors, ev = definiteness(B)
            o["Delta"][bn] = {"matrix": mstr(B), "verdict": verdict, "leading_minors": [str(m) for m in minors],
                              "eigenvalues_float": ev}
            exp = "negative definite" if plant else "positive definite"
            check(f"l={L} {parity}: Delta on {bn} block ({kinds[0] if bn == 'psi' else kinds[1] if bn == 'psi_dot' else '8x8'})"
                  f" is {exp} by exact leading minors", verdict == exp, verdict)
            print(f"      minor signs {['+' if m > 0 else '-' if m < 0 else '0' for m in minors]}; eigenvalues {fmt(ev)}")
            if bn != "full":
                sl = slice(0, NB) if bn == "psi" else slice(NB, 2 * NB)
                rs = rel_spectrum(G2[sl, sl], G3[sl, sl])
                o["Delta"][bn]["Q2_over_S_spectrum_float"] = rs
                print(f"      Delta diag {[str(B[i, i]) for i in range(NB)]}; row 0 {[str(e) for e in B.row(0)]}")
                print(f"      spectrum of Q_2 relative to S on this block: {fmt(rs)}")
        X = D[:NB, NB:]
        o["Delta"]["cross"] = {"matrix": mstr(X), "is_zero": bool(X.is_zero_matrix)}
        check(f"l={L} {parity}: cross block of Delta is non-zero", not X.is_zero_matrix,
              f"row 0 {[str(e) for e in X.row(0)]}")
        o["static_ratio_psi00"] = str(G2[0, 0] / G3[0, 0])
        o["first_psi_dot_ratio"] = str(G2[NB, NB] / G3[NB, NB])
        print(f"   [l={L} {parity}] Q_2/S on psi = phi ({kinds[0]}): {o['static_ratio_psi00']};"
              f" on psi_dot = phi ({kinds[1]}): {o['first_psi_dot_ratio']}")
        out[parity] = o
        keep[parity] = (G2, G3)

    # static reproductions with null-cone's own data
    cvE = nc.Curv(nc.metric_electric(phi))
    B1 = nc.Curv(nc.tt_series(K0, 1)).B()
    lam = sp.cancel(B1[0, 0] / Hphi[0, 0])
    check(f"l={L}: B of h = 2tK_0 is a constant multiple of Hess phi", (B1 - lam * Hphi).applyfunc(sp.expand).is_zero_matrix
          and lam.is_number, lam)
    cvB = nc.Curv(nc.tt_series(K0 / lam, 1))
    stat = {}
    for name, cv in (("static electric tide", cvE), ("stationary magnetic tide h = 2tK", cvB)):
        C = cv.cone_weyl()
        q2 = 2 * nc.cone_int(nc.norm2(nc.moment(C, P_CHARGE)), W)
        q3 = 2 * nc.cone_int(nc.norm2(nc.moment(C, P_ENTROPY)), W)
        stat[name] = (q2, q3)
        r = sp.cancel(q2 / q3)
        out.setdefault("static", {})[name] = str(r)
        if plant is None:
            check(f"l={L}: {name} gives Q_2/S = ((l+2)/(l+1))^2 = {target}", r == target, r)
        else:
            check(f"l={L}: PLANT {plant} caught, {name} no longer gives {target}", r != target, r)
    G2e, G3e = keep["even"]
    G2o, G3o = keep["odd"]
    nrm = (L * (L - 1)) ** 2 * N_Y                      # both static chain fields are l(l-1) phi
    q2, q3 = stat["static electric tide"]
    check(f"l={L}: even Gram [psi0, psi0] reproduces the static electric tide (pure metric datum)",
          sp.simplify(q2 - nrm * G2e[0, 0]) == 0 and sp.simplify(q3 - nrm * G3e[0, 0]) == 0 and out["even"]["psi_is"] == "metric")
    q2, q3 = stat["stationary magnetic tide h = 2tK"]
    check(f"l={L}: odd Gram [psi0, psi0] reproduces the stationary magnetic tide (pure extrinsic-curvature datum)",
          sp.simplify(q2 - nrm * G2o[0, 0]) == 0 and sp.simplify(q3 - nrm * G3o[0, 0]) == 0 and out["odd"]["psi_is"] == "extrinsic")
    check(f"l={L}: duality, even and odd chain-basis Gram matrices identical (Q_2 and Q_3)",
          (G2e - G2o).is_zero_matrix and (G3e - G3o).is_zero_matrix)
    if L == 2:
        d00 = G2e[0, 0] - G3e[0, 0]
        if plant is None:
            check("l=2: reproduces RULER's Delta[psi_0, psi_0] = 1/36 and Delta[psi_dot_0, psi_dot_0] = 1/1008",
                  d00 == sp.Rational(1, 36) and G2e[NB, NB] - G3e[NB, NB] == sp.Rational(1, 1008), d00)
    out["seconds"] = round(time.time() - tl)
    print(f"   l={L} done in {out['seconds']} s", flush=True)
    res["per_l"][str(L)] = out


# ====================================================================== part 2
def qform(p, m, n):
    """One generator, C = v^m and v^n:  int_0^1 v^3 (1 - v) sigma_p[v^m] sigma_p[v^n] dv."""
    return (Fraction(1, m + n + 6) - Fraction(1, m + n + 7)) / ((m + p + 1) * (n + p + 1))


def hform(m, n):
    """int_0^1 h_m h_n / v dv with h = v^2 (sigma_2 - sigma_3), h_m = -v^(m+3)/((m+3)(m+4))."""
    return Fraction(1, (m + 3) * (m + 4) * (n + 3) * (n + 4) * (m + n + 6))


def quad(c, p):
    idx = [i for i, a in enumerate(c) if a != 0]
    return sum(c[a] * c[b] * qform(p, a, b) for a in idx for b in idx)


def part2():
    print("\n==================== part 2: one generator ====================", flush=True)
    pq, pe = P_CHARGE, P_ENTROPY
    # (a) the identity, symbolically in the two exponents and exactly on a grid
    m, n = sp.symbols("m n", positive=True)
    q = lambda p: (1 / (m + n + 6) - 1 / (m + n + 7)) / ((m + p + 1) * (n + p + 1))
    ident = sp.simplify(q(pq) - q(pe) - 1 / ((m + 3) * (m + 4) * (n + 3) * (n + 4) * (m + n + 6)))
    check("part 2: Q_2 - Q_3 = int h^2/v on all pairs of monomials (symbolic exponents)", ident == 0, ident)
    check("part 2: the same identity, exact rationals, exponents 0..40",
          all(qform(pq, a, b) - qform(pe, a, b) == hform(a, b) for a in range(41) for b in range(41)))
    # (b) monomials
    mono = {}
    print("   monomials C = v^m: exact Q_2/Q_3 beside the profile law ((m+4)/(m+3))^2 and the bounds 1 < ratio < infinity")
    good = True
    for a in range(9):
        r = qform(pq, a, a) / qform(pe, a, a)
        mono[a] = str(r)
        good &= r == Fraction(a + 4, a + 3) ** 2 and r > 1
        print(f"      m={a}: {r} = {float(r):.6f}")
    check("part 2: monomials m = 0..8 give ((m+4)/(m+3))^2 > 1 exactly", good)
    # (c) exact extremes on polynomial C of degree <= d (generalised eigenvalues of the pencil)
    pen = {}
    lam_ = sp.Symbol("lam_")
    print("   polynomial C of degree <= d: exact pencil extremes of Q_2/Q_3 (floats from exact matrices, 60 digits)")
    for d in range(0, 13):
        A2 = sp.Matrix(d + 1, d + 1, lambda a, b: sp.Rational(qform(pq, a, b).numerator, qform(pq, a, b).denominator))
        A3 = sp.Matrix(d + 1, d + 1, lambda a, b: sp.Rational(qform(pe, a, b).numerator, qform(pe, a, b).denominator))
        roots = sorted(sp.re(r) for r in sp.Poly((A2 - lam_ * A3).det(), lam_).nroots(n=60, maxsteps=5000))
        pen[d] = (float(roots[0]), float(roots[-1]), A2, A3)
        print(f"      d={d:2d}: min {float(roots[0]):.9f}   max {float(roots[-1]):.6f}")
    lo8, hi8, A2, A3 = pen[8]
    pos = lambda G: all(G[:k, :k].det() > 0 for k in range(1, G.shape[0] + 1))
    if plant is None:
        check("part 2: degree <= 8, Q_2 - Q_3 positive definite (exact minors), so every ratio exceeds 1", pos(A2 - A3))
        up, dn = sp.Rational(math.ceil(hi8 * 10**6), 10**6), sp.Rational(math.floor(hi8 * 10**6), 10**6)
        check(f"part 2: degree <= 8, exact bracket of the largest ratio: {dn} < max < {up}",
              pos(up * A3 - A2) and not pos(dn * A3 - A2))
        up, dn = sp.Rational(math.ceil(lo8 * 10**9), 10**9), sp.Rational(math.floor(lo8 * 10**9), 10**9)
        check(f"part 2: degree <= 8, exact bracket of the smallest ratio: {dn} < min < {up}",
              pos(A2 - dn * A3) and not pos(A2 - up * A3))
    else:
        check(f"part 2: PLANT {plant}: degree <= 8, Q_2 - Q_3 positive definite", pos(A2 - A3))
    # (d) random search, exact ratios of random integer polynomials of degree <= 8
    rng = random.Random(7315)
    best = {"monomial basis": [None, None], "Bernstein basis": [None, None]}
    bern = [[Fraction(math.comb(8, a) * math.comb(8 - a, k - a) * (-1) ** (k - a)) if k >= a else Fraction(0)
             for k in range(9)] for a in range(9)]                      # v^a (1 - v)^(8 - a) in monomials
    for trial in range(20000):
        deg = rng.randint(0, 8)
        c = [Fraction(rng.randint(-50, 50)) if k <= deg else Fraction(0) for k in range(9)]
        w = [Fraction(rng.randint(-50, 50)) for _ in range(9)]
        cb = [sum(w[a] * bern[a][k] for a in range(9)) for k in range(9)]
        for key, cc in (("monomial basis", c), ("Bernstein basis", cb)):
            if not any(cc):
                continue
            r = quad(cc, pq) / quad(cc, pe)
            b = best[key]
            b[0] = r if b[0] is None or r < b[0] else b[0]
            b[1] = r if b[1] is None or r > b[1] else b[1]
    smin = min(b[0] for b in best.values())
    smax = max(b[1] for b in best.values())
    for key, b in best.items():
        print(f"   random search, 20000 polynomials, {key}: smallest {float(b[0]):.6f}, largest {float(b[1]):.6f}")
    print(f"   beside: proved infimum 1 (not attained), proved supremum +infinity; degree <= 8 exact range ({lo8:.9f}, {hi8:.6f})")
    check("part 2: every random ratio lies strictly above 1 and inside the exact degree-8 range",
          smin > 1 and float(smin) >= lo8 * (1 - 1e-12) and float(smax) <= hi8 * (1 + 1e-12),
          f"search range [{float(smin):.6f}, {float(smax):.6f}]")
    # (e) the family that sends the ratio to infinity: h' = v^2 (1 - v)^N, i.e. C_N = 4 (1-v)^N - N v (1-v)^(N-1)
    fam = {}
    good = True
    print("   family C_N = 4(1-v)^N - N v (1-v)^(N-1): exact Q_2/Q_3 from the definitions, beside the proved bound 1 + ln(N)/270")
    for N in (1, 2, 4, 8, 16, 32, 64, 128):
        c = [Fraction(0)] * (N + 1)
        for k in range(N + 1):
            c[k] += 4 * math.comb(N, k) * (-1) ** k
            if k >= 1:
                c[k] += -N * math.comb(N - 1, k - 1) * (-1) ** (k - 1)
        r = quad(c, pq) / quad(c, pe)
        bound = 1 + math.log(N) / 270
        fam[N] = float(r)
        good &= float(r) >= bound
        print(f"      N={N:4d}: ratio {float(r):.6f}   bound {bound:.6f}   (ratio - 1)/ln N = {((float(r) - 1) / math.log(N)) if N > 1 else float('nan'):.4f}")
    check("part 2: the family respects the proved lower bound 1 + ln(N)/270 and grows with N",
          good and all(fam[a] < fam[b] for a, b in zip((1, 2, 4, 8, 16, 32, 64), (2, 4, 8, 16, 32, 64, 128))))
    res["part2"] = {"monomials": mono, "pencil": {str(d): [pen[d][0], pen[d][1]] for d in pen},
                    "random": {k: [float(b[0]), float(b[1])] for k, b in best.items()}, "family": fam}


for L_ in LS:
    run_l(L_)
if "--no-part2" not in sys.argv:
    part2()

res["all_pass"] = ok
with open(sys.argv[1], "w", encoding="utf8", newline="\n") as fh:
    json.dump(res, fh, indent=1, sort_keys=True)
nf = sum(1 for c in res["checks"] if not c["pass"])
print(f"\nplant={plant} l={LS} checks={len(res['checks'])} fail={nf} time={time.time() - t0:.0f} s")
