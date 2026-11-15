"""SS-TENSION T1-T4: vacuum two-point functions restricted to the light cone t = r = u,
projected on angular sectors, compared with Satishchandran-Sorce's ultralocal states
(1.4)/(3.4) scalar, (B.13) Maxwell, (B.26) graviton (sign of B.13).

Exact Gaussian-rational arithmetic except one logarithm per integral (mpmath, 60 digits).
Usage: python cone_sector.py            (all tests)
"""
import sys, json, time
import sympy as sp
import mpmath as mp

mp.mp.dps = 60
c, s, Z = sp.symbols('c s Z')
I = sp.I
R = sp.Rational
RESULTS = []


def line(tag, ok, msg, plant=False):
    status = ('CAUGHT' if not ok else 'MISSED') if plant else ('PASS' if ok else 'FAIL')
    RESULTS.append((tag, status, msg))
    print(f"[{status}] {tag}: {msg}", flush=True)


def dot(a, b):
    return -a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3]


def geometry(u1, u2, eps):
    """Points and frames. Point 1 on the generator n1 = z, point 2 at polar angle theta
    in the x-z plane; frames aligned along the connecting great circle."""
    k1 = (1, 0, 0, 1)
    m1 = (0, 1 / sp.sqrt(2), I / sp.sqrt(2), 0)
    k2 = (1, s, 0, c)
    m2 = (0, c / sp.sqrt(2), I / sp.sqrt(2), -s / sp.sqrt(2))
    m2b = (0, c / sp.sqrt(2), -I / sp.sqrt(2), -s / sp.sqrt(2))
    xi = (u1 - u2 - I * eps, -u2 * s, 0, u1 - u2 * c)
    xi2 = sp.expand(-(u1 - u2 - I * eps) ** 2 + u1 ** 2 + u2 ** 2 - 2 * u1 * u2 * c)
    return k1, m1, k2, m2, m2b, xi, xi2


def fder(p, n):
    """n-th derivative of f(Z) = Z^-p."""
    co = sp.Integer(1)
    for j in range(n):
        co *= (-p - j)
    return co * Z ** (-p - n)


def d2(p, xi, a, b):
    sa, sb = 2 * dot(a, xi), 2 * dot(b, xi)
    return fder(p, 2) * sa * sb + fder(p, 1) * 2 * dot(a, b)


def d4(p, xi, a, b, cc, d):
    v = [a, b, cc, d]
    s1 = [2 * dot(x, xi) for x in v]
    s2 = [[2 * dot(x, y) for y in v] for x in v]
    t4 = s1[0] * s1[1] * s1[2] * s1[3]
    t3 = (s2[0][1] * s1[2] * s1[3] + s2[0][2] * s1[1] * s1[3] + s2[0][3] * s1[1] * s1[2]
          + s2[1][2] * s1[0] * s1[3] + s2[1][3] * s1[0] * s1[2] + s2[2][3] * s1[0] * s1[1])
    t2 = s2[0][1] * s2[2][3] + s2[0][2] * s2[1][3] + s2[0][3] * s2[1][2]
    return fder(p, 4) * t4 + fder(p, 3) * t3 + fder(p, 2) * t2


def norm_G(D):
    """Massless scalar Wightman function = norm_G(D) * (xi^2)^-(D-2)/2."""
    return sp.gamma(R(D, 2) - 1) / (4 * sp.pi ** R(D, 2))


def W_scalar(u1s, u2s, eps, D=4, weight=None):
    """<Pi Pi>, Pi = d_u(u^w phi) along the generator, w = (D-2)/2 unless planted."""
    p = R(D - 2, 2)
    w = p if weight is None else weight
    a, b = sp.symbols('a b', positive=True)
    xi2 = -(a - b - I * eps) ** 2 + a ** 2 + b ** 2 - 2 * a * b * c
    e = sp.diff((a * b) ** w * xi2 ** (-p), a, b)
    return sp.simplify(e.subs({a: u1s, b: u2s}))


def W_maxwell(u1, u2, eps, D=4, bar=True):
    """<Phi0(1) conj Phi0(2)>, Phi0 = F(k, m), Feynman gauge <A A> = eta G."""
    p = R(D - 2, 2)
    k1, m1, k2, m2, m2b, xi, xi2 = geometry(u1, u2, eps)
    mm = m2b if bar else m2
    # derivative at point 2 is minus the derivative in xi
    e = -(dot(m1, mm) * d2(p, xi, k1, k2) - dot(m1, k2) * d2(p, xi, k1, mm)
          - dot(k1, mm) * d2(p, xi, m1, k2) + dot(k1, k2) * d2(p, xi, m1, mm))
    return e, xi2


def W_graviton(u1, u2, eps, D=4, bar=True):
    """<Psi0(1) conj Psi0(2)> / G, Psi0 = R(k, m, k, m), de Donder gauge
    <gamma gamma> = 16 pi G (eta eta + eta eta - 2/(D-2) eta eta) G_scalar.
    The trace term drops because k.k = k.m = m.m = 0."""
    p = R(D - 2, 2)
    k1, m1, k2, m2, m2b, xi, xi2 = geometry(u1, u2, eps)
    mm = m2b if bar else m2
    half = R(1, 2)
    T1 = [(1, k1, m1, (m1, k1)), (-half, m1, m1, (k1, k1)), (-half, k1, k1, (m1, m1))]
    T2 = [(1, k2, mm, (mm, k2)), (-half, mm, mm, (k2, k2)), (-half, k2, k2, (mm, mm))]
    e = 0
    for c1, a1, b1, dd1 in T1:
        for c2, a2, b2, dd2 in T2:
            ten = dot(a1, a2) * dot(b1, b2) + dot(a1, b2) * dot(b1, a2)
            if ten == 0:
                continue
            e += c1 * c2 * ten * d4(p, xi, dd1[0], dd1[1], dd2[0], dd2[1])
    return 16 * e, xi2   # times pi * norm_G(D), applied numerically


def reduce_s(expr):
    """Polynomial in c after s^2 -> 1 - c^2; odd powers of s must cancel."""
    expr = sp.expand(expr)
    P = sp.Poly(expr, s)
    out = 0
    for (k,), co in P.terms():
        if k % 2:
            if sp.simplify(co) != 0:
                raise ValueError("odd power of sin(theta) survives")
            continue
        out += co * (1 - c ** 2) ** (k // 2)
    return sp.expand(out)


def dfun(l, sw):
    """Wigner d^l_{ss}(theta) as a polynomial in c = cos(theta)."""
    return sp.expand(((1 + c) / 2) ** sw * sp.jacobi(l - sw, 0, 2 * sw, c))


def project(Wc, xi2, u1, u2, eps, l, sw):
    """g_l = 2 pi int_{-1}^{1} W(c) d^l_{ss}(c) dc, W = sum_m P_m(c) Z^-m, Z = alpha - beta c."""
    alpha = sp.expand(xi2.subs(c, 0))
    beta = -sp.expand(xi2.coeff(c, 1))
    z = sp.symbols('z')
    integrand = sp.expand(Wc * dfun(l, sw))
    # collect powers of Z (negative integers)
    terms = sp.Poly(integrand.subs(Z, 1 / sp.Symbol('Y')), sp.Symbol('Y')).terms()
    laurent = {}
    for (m,), co in terms:
        poly = sp.Poly(sp.expand(co.subs(c, (alpha - z) / beta)), z)
        for (j,), a in poly.terms():
            laurent[j - m] = laurent.get(j - m, 0) + a
    zlo, zhi = sp.expand(alpha - beta), sp.expand(alpha + beta)
    exact = 0
    logco = 0
    for j, a in laurent.items():
        if j == -1:
            logco += a
        else:
            exact += a * (zhi ** (j + 1) - zlo ** (j + 1)) / (j + 1)
    val = to_mpc(exact)
    if logco != 0:
        val += to_mpc(logco) * (mp.log(to_mpc(zhi)) - mp.log(to_mpc(zlo)))
    return 2 * mp.pi * val / to_mpc(beta)


def to_mpc(expr):
    re_, im_ = sp.expand(expr).as_real_imag()
    re_, im_ = sp.Rational(re_), sp.Rational(im_)
    return mp.mpc(mp.mpf(re_.p) / re_.q, mp.mpf(im_.p) / im_.q)


NORM4 = 1 / (4 * mp.pi ** 2)          # norm_G(4)


def run_field(name, sw, ls, builder, target, power, PTS, EPS, tol=1e-6, plant=False, with_u=True, fac=NORM4):
    """target: predicted limit of g_l * (u1 u2 if with_u) * (u1 - u2)^power."""
    worst = 0
    table = []
    for (u1, u2) in PTS:
        for eps in EPS:
            W, xi2 = builder(u1, u2, eps)
            Wc = reduce_s(W)
            row = []
            for l in ls:
                g = project(Wc, xi2, u1, u2, eps, l, sw) * fac
                scaled = g * (mp.mpf(u1 - u2)) ** power * (mp.mpf(u1 * u2) if with_u else 1)
                row.append(scaled)
            table.append((str(u1), str(u2), str(eps), [mp.nstr(x, 12) for x in row]))
            if eps == EPS[-1]:
                for x in row:
                    worst = max(worst, abs(x / target - 1))
    for t in table:
        print("    u1=%s u2=%s eps=%s  scaled g_l (l=%s): %s" % (t[0], t[1], t[2], list(ls), t[3]))
    ok = worst < tol
    line(name, ok, f"target {mp.nstr(target, 12)}, worst relative deviation at smallest eps = {mp.nstr(worst, 3)}",
         plant=plant)
    return worst


def main():
    t0 = time.time()
    PTS = [(R(7, 10), R(3, 10)), (R(2, 5), R(9, 10)), (R(1, 8), R(19, 20)), (R(11, 20), R(1, 2))]
    EPS = [R(1, 10 ** 4), R(1, 10 ** 8), R(1, 10 ** 12)]

    # ---- T1 control (i): scalar, D = 4 --------------------------------------------------
    def b_scalar(u1, u2, eps):
        xi2 = sp.expand(-(u1 - u2 - I * eps) ** 2 + u1 ** 2 + u2 ** 2 - 2 * u1 * u2 * c)
        e = W_scalar(u1, u2, eps)
        # rewrite as polynomial in c over powers of Z
        num, den = sp.fraction(sp.together(e))
        # den = const * (xi2)^3 ; express through Z
        k = sp.Poly(den, c).degree()
        lead = sp.simplify(den / xi2 ** k)
        return sp.expand(num / lead) * Z ** (-k), xi2

    run_field("T1 scalar <Pi Pi>, Pi = d_u(u phi), D=4, l=0..4", 0, range(0, 5), b_scalar,
              -1 / (4 * mp.pi), 2, PTS, EPS, with_u=False)

    def b_scalar_plant(u1, u2, eps):
        xi2 = sp.expand(-(u1 - u2 - I * eps) ** 2 + u1 ** 2 + u2 ** 2 - 2 * u1 * u2 * c)
        e = W_scalar(u1, u2, eps, weight=2)
        num, den = sp.fraction(sp.together(e))
        k = sp.Poly(den, c).degree()
        lead = sp.simplify(den / xi2 ** k)
        return sp.expand(num / lead) * Z ** (-k), xi2

    run_field("P1 plant: scalar with Pi = d_u(u^2 phi)", 0, range(0, 3), b_scalar_plant,
              -1 / (4 * mp.pi), 2, PTS[:1], EPS[1:], plant=True, with_u=False)

    # ---- T2 control (ii): Maxwell D = 4 -------------------------------------------------
    run_field("T2 Maxwell <Phi0 conj Phi0>, D=4, l=1..4 (SS B.13)", 1, range(1, 5),
              lambda a, b, e: W_maxwell(a, b, e), -1 / (4 * mp.pi), 2, PTS, EPS)

    # ---- T3 decisive: graviton D = 4 ----------------------------------------------------
    run_field("T3 graviton <Psi0 conj Psi0>/G, D=4, l=2..6 (SS B.26, kappa=-1/2)", 2, range(2, 7),
              lambda a, b, e: W_graviton(a, b, e), mp.mpf(12), 4, PTS, EPS, fac=mp.pi * NORM4)
    run_field("P3 plant: graviton prediction with kappa=-1 (48)", 2, range(2, 3),
              lambda a, b, e: W_graviton(a, b, e), mp.mpf(48), 4, PTS[:1], EPS[2:], plant=True, fac=mp.pi * NORM4)

    # unbarred correlators must vanish identically (helicity)
    Wg, _ = W_graviton(R(7, 10), R(3, 10), R(1, 10 ** 4), bar=False)
    Wm, _ = W_maxwell(R(7, 10), R(3, 10), R(1, 10 ** 4), bar=False)
    # Z is tied to c; substitute and test at a numeric angle instead of relying on symbol independence
    def at_angle(W, xi2, cv):
        e = W.subs(Z, xi2).subs(s, sp.sqrt(1 - c ** 2)).subs(c, cv)
        return sp.N(e, 30)
    k1, m1, k2, m2, m2b, xi, xi2 = geometry(R(7, 10), R(3, 10), R(1, 10 ** 4))
    vg = at_angle(Wg, xi2, R(1, 3)); vm = at_angle(Wm, xi2, R(1, 3))
    line("T3b unbarred <Psi0 Psi0> and <Phi0 Phi0> vanish (eps = 1e-4, cos = 1/3)",
         abs(vg) < 1e-20 and abs(vm) < 1e-20, f"|graviton| = {sp.N(abs(vg), 3)}, |Maxwell| = {sp.N(abs(vm), 3)}")

    # ---- T4 dimension: unregulated correlator on two different generators ---------------
    print("  T4: correlator at eps = 0, u1 = 7/10, u2 = 3/10, cos(theta) = 1/3 (0 means ultralocal)")
    vals = {}
    for D in (4, 5, 6):
        k1, m1, k2, m2, m2b, xi, xi2 = geometry(R(7, 10), R(3, 10), 0)
        Wm, _ = W_maxwell(R(7, 10), R(3, 10), 0, D=D)
        Wg, _ = W_graviton(R(7, 10), R(3, 10), 0, D=D)
        vm = sp.N(Wm.subs(Z, xi2).subs(s, sp.sqrt(1 - c ** 2)).subs(c, R(1, 3)), 20)
        vg = sp.N(Wg.subs(Z, xi2).subs(s, sp.sqrt(1 - c ** 2)).subs(c, R(1, 3)), 20)
        a_, b_ = sp.symbols('a_ b_', positive=True)
        p = R(D - 2, 2)
        vs = sp.simplify(sp.diff((a_ * b_) ** p * (2 * a_ * b_ * (1 - c)) ** (-p), a_, b_))
        vals[D] = (vs, vm, vg)
        print(f"    D={D}: scalar {vs}, Maxwell {sp.N(abs(vm), 6)}, graviton/G {sp.N(abs(vg), 6)}")
    ok4 = (vals[4][0] == 0 and abs(vals[4][1]) < 1e-15 and abs(vals[4][2]) < 1e-15)
    line("T4a D=4: scalar, Maxwell, graviton all vanish off the diagonal", ok4, "exact zero" if ok4 else "non-zero")
    same = all(vals[D][0] == 0 and abs(vals[D][1]) < 1e-15 and abs(vals[D][2]) < 1e-15 for D in (5, 6))
    # PREREG expected non-zero at D = 5, 6. That expectation was wrong: F(k, m) and R(k, m, k, m) contain
    # only derivatives tangent to the cone, so the off-diagonal correlator vanishes in every D.
    # T4 therefore does not discriminate dimensions; reported as information, not as a pass.
    RESULTS.append(("T4b D=5,6 off-diagonal", "INFO", "also zero: PREREG expectation wrong, test uninformative" if same else "non-zero"))
    print("[INFO] T4b D=5,6: " + RESULTS[-1][2])
    n_bad = sum(1 for r in RESULTS if r[1] in ('FAIL', 'MISSED'))
    print(f"\n{len(RESULTS)} lines, {n_bad} unexpected, {time.time() - t0:.1f} s")
    with open(sys.argv[1] if len(sys.argv) > 1 else 'cone_sector.json', 'w') as fh:
        json.dump([list(r) for r in RESULTS], fh, indent=1)


if __name__ == '__main__':
    main()
