"""SS-TENSION T0, T5, T6, T7 (exact, sympy).

T0  Gaussian null gauge on the cone, D = 4: R_{r A r B} = kappa * r * d_r^2 (s_AB / r).
T5  On polynomial plane-wave solutions: with S the tip-regular solution of T0's equation and
    Q = S - (2/u) int_0^u S, div div Q = a u (x.E.x) and curl div Q = b u (x.B.x) on the cone.
T6  The transform Q -> F = Q - 2u int_u^inf Q/u'^2 preserves the half-line ultralocal form.
T7  Flows on monomials: vacuum flow (weight 0 on Q) versus literal weight-0 flow on S.
"""
import sys, json, time, itertools
import sympy as sp

R = sp.Rational
RESULTS = []


def line(tag, ok, msg, plant=False):
    status = ('CAUGHT' if not ok else 'MISSED') if plant else ('PASS' if ok else 'FAIL')
    RESULTS.append((tag, status, msg))
    print(f"[{status}] {tag}: {msg}", flush=True)


# ------------------------------------------------------------------ T0
def t0():
    w, r, th, ph = X = sp.symbols('w r theta phi', positive=True)
    a = sp.Function('a')(*X)
    b = sp.Function('b')(*X)
    g0 = sp.Matrix([[-1, -1, 0, 0], [-1, 0, 0, 0], [0, 0, r ** 2, 0], [0, 0, 0, r ** 2 * sp.sin(th) ** 2]])
    # s_AB in the frame that is orthonormal for the unit sphere: s_thth = a, s_thph = b, s_phph = -a
    h = sp.zeros(4, 4)
    h[2, 2] = a
    h[2, 3] = h[3, 2] = b * sp.sin(th)
    h[3, 3] = -a * sp.sin(th) ** 2
    gi = g0.inv()

    def gam_low(g):
        return [[[R(1, 2) * (sp.diff(g[l, m], X[n]) + sp.diff(g[l, n], X[m]) - sp.diff(g[m, n], X[l]))
                  for n in range(4)] for m in range(4)] for l in range(4)]

    G0l, G1l = gam_low(g0), gam_low(h)
    gi1 = -gi * h * gi
    G0 = [[[sum(gi[l, k] * G0l[k][m][n] for k in range(4)) for n in range(4)] for m in range(4)] for l in range(4)]
    G1 = [[[sum(gi[l, k] * G1l[k][m][n] + gi1[l, k] * G0l[k][m][n] for k in range(4)) for n in range(4)]
           for m in range(4)] for l in range(4)]

    def riem1(l, s_, m, n):   # first-order R^l_{s m n}
        e = sp.diff(G1[l][n][s_], X[m]) - sp.diff(G1[l][m][s_], X[n])
        for k in range(4):
            e += G0[l][m][k] * G1[k][n][s_] + G1[l][m][k] * G0[k][n][s_]
            e -= G0[l][n][k] * G1[k][m][s_] + G1[l][n][k] * G0[k][m][s_]
        return e

    def low(rho, s_, m, n):
        return sp.simplify(sum(g0[rho, l] * riem1(l, s_, m, n) for l in range(4)))

    Rthth = low(1, 2, 1, 2)
    Rthph = sp.simplify(low(1, 2, 1, 3) / sp.sin(th))
    Rphph = sp.simplify(low(1, 3, 1, 3) / sp.sin(th) ** 2)
    kap = -R(1, 2)
    e1 = sp.simplify(Rthth - kap * r * sp.diff(a / r, r, 2))
    e2 = sp.simplify(Rthph - kap * r * sp.diff(b / r, r, 2))
    e3 = sp.simplify(Rphph + kap * r * sp.diff(a / r, r, 2))
    line("T0 R_{rArB} = -(1/2) r d_r^2(s_AB/r) in Gaussian null gauge (theta-theta, theta-phi, phi-phi)",
         e1 == 0 and e2 == 0 and e3 == 0, f"residuals {e1}, {e2}, {e3}")
    ep = sp.simplify(Rthth + R(1, 2) * sp.diff(a, r, 2))
    line("P0 plant: R_{rArB} = -(1/2) d_r^2 s_AB (no background terms)", ep == 0, f"residual {ep}", plant=True)


# ------------------------------------------------------------------ T5
def t5():
    t, x, y, z = Xc = sp.symbols('t x y z', real=True)
    u, th, ph = sp.symbols('u theta phi', positive=True)
    kap = -R(1, 2)
    consts = []
    for n in range(2, 7):
        for pol in ('plus', 'cross', 'mixed'):
            f = (t - z) ** n
            gam = sp.zeros(4, 4)
            cp, cx = {'plus': (1, 0), 'cross': (0, 1), 'mixed': (2, -3)}[pol]
            gam[1, 1], gam[2, 2] = cp * f, -cp * f
            gam[1, 2] = gam[2, 1] = cx * f

            def Rm(m, al, nu, be):
                return R(1, 2) * (sp.diff(gam[m, be], Xc[al], Xc[nu]) + sp.diff(gam[al, nu], Xc[m], Xc[be])
                                  - sp.diff(gam[al, be], Xc[m], Xc[nu]) - sp.diff(gam[m, nu], Xc[al], Xc[be]))

            Riem = {idx: Rm(*idx) for idx in itertools.product(range(4), repeat=4)}
            cone = {t: u, x: u * sp.sin(th) * sp.cos(ph), y: u * sp.sin(th) * sp.sin(ph), z: u * sp.cos(th)}
            nvec = [sp.sin(th) * sp.cos(ph), sp.sin(th) * sp.sin(ph), sp.cos(th)]
            kv = [1] + nvec
            eth = [0, sp.cos(th) * sp.cos(ph), sp.cos(th) * sp.sin(ph), -sp.sin(th)]
            eph = [0, -sp.sin(ph), sp.cos(ph), 0]

            def con(v1, v2, v3, v4):
                return sum(Riem[(i, j, k, l)].subs(cone) * v1[i] * v2[j] * v3[k] * v4[l]
                           for i, j, k, l in itertools.product(range(4), repeat=4)
                           if Riem[(i, j, k, l)] != 0)

            Cthth, Cthph, Cphph = con(kv, eth, kv, eth), con(kv, eth, kv, eph), con(kv, eph, kv, eph)
            # tip-regular S from S'' = u C / kappa; C is homogeneous of degree n-2 in u
            S1 = Cthth * u ** 3 / (kap * n * (n + 1))
            S2 = Cthph * u ** 3 / (kap * n * (n + 1))
            q1, q2 = S1 * R(n, n + 2), S2 * R(n, n + 2)       # Q = S - (2/u) int_0^u S
            Vth = sp.diff(q1, th) + 2 * sp.cot(th) * q1 + sp.diff(q2, ph) / sp.sin(th)
            Vph = sp.diff(q2, th) + 2 * sp.cot(th) * q2 - sp.diff(q1, ph) / sp.sin(th)
            dd = (sp.diff(sp.sin(th) * Vth, th) + sp.diff(Vph, ph)) / sp.sin(th)
            cd = (sp.diff(sp.sin(th) * Vph, th) - sp.diff(Vth, ph)) / sp.sin(th)
            # E_ij = R_{0i0j}, B_ij = (1/2) eps_ikl R_{kl0j}; psi = x.E.x, x.B.x on the cone
            psiE = sum(Riem[(0, i, 0, j)].subs(cone) * nvec[i - 1] * nvec[j - 1] for i in (1, 2, 3) for j in (1, 2, 3)) * u ** 2
            psiB = sum(R(1, 2) * sp.LeviCivita(i, k, l) * Riem[(k, l, 0, j)].subs(cone) * nvec[i - 1] * nvec[j - 1]
                       for i in (1, 2, 3) for j in (1, 2, 3) for k in (1, 2, 3) for l in (1, 2, 3)) * u ** 2
            tf = sp.simplify(Cthth + Cphph)
            pts = [{u: R(3, 7), th: R(4, 5), ph: R(1, 3)}, {u: R(5, 4), th: R(11, 9), ph: R(-7, 5)},
                   {u: R(2, 3), th: R(2, 1), ph: R(9, 4)}]
            ra, rb = [], []
            for p in pts:
                nE, nB = sp.N((u * psiE).subs(p), 40), sp.N((u * psiB).subs(p), 40)
                ra.append(sp.N(dd.subs(p), 40) / nE)
                rb.append(sp.N(cd.subs(p), 40) / nB)
            consts.append((n, pol, ra, rb, tf))
    a0, b0 = consts[0][2][0], consts[0][3][0]
    dev = max(max(abs(v / a0 - 1) for v in c_[2]) for c_ in consts)
    devb = max(max(abs(v / b0 - 1) for v in c_[3]) for c_ in consts)
    tfok = all(c_[4] == 0 for c_ in consts)
    print(f"    div div Q / (u x.E.x) = {sp.nsimplify(a0, rational=True)}, curl div Q / (u x.B.x) = {sp.nsimplify(b0, rational=True)}")
    line("T5a C_AB trace-free on the cone for all 15 solutions", tfok, "exact")
    line("T5b div div Q = a u (x.E.x), one constant, degrees 2..6, three polarisations, three points",
         dev < R(1, 10 ** 30), f"a = {sp.nsimplify(a0, rational=True)}, max relative spread {sp.N(dev, 3)}")
    line("T5c curl div Q = b u (x.B.x), one constant, same set",
         devb < R(1, 10 ** 30), f"b = {sp.nsimplify(b0, rational=True)}, max relative spread {sp.N(devb, 3)}")
    # plant: use S itself in place of Q (no integral term): constant must move with the degree
    dev_pl = max(abs((c_[2][0] * R(c_[0] + 2, c_[0])) / (a0 * 2) - 1) for c_ in consts)
    line("P5 plant: S in place of Q (drop the integral term)", dev_pl < R(1, 10 ** 30),
         f"max relative spread {sp.N(dev_pl, 3)}", plant=True)


# ------------------------------------------------------------------ T6
def t6():
    k = sp.symbols('k', positive=True)

    def ft(poly_coeffs):   # sum c_m u^m e^{-u}  ->  sum c_m m!/(1+ik)^{m+1}
        return sum(cm * sp.factorial(m) / (1 + sp.I * k) ** (m + 1) for m, cm in poly_coeffs.items())

    def Qc(n):
        return {n: sp.Integer(1)}

    def Fc(n, coef=2):   # Q - coef * u int_u^inf u'^{n-2} e^{-u'} du'
        d = {n: sp.Integer(1)}
        for j in range(0, n - 1):
            d[j + 1] = d.get(j + 1, 0) - coef * sp.factorial(n - 2) / sp.factorial(j)
        return d

    def form(A, B):
        # int_0^inf k conj(A^)(k) B^(k) dk, exactly: every term is k^j/(1+k^2)^M, a Beta function.
        # (sympy's integrate returns unevaluated Meijer G here; see logs/shear_curv_run1.log)
        tot = 0
        for ma, ca in A.items():
            for mb, cb in B.items():
                M = ma + mb + 2
                num = sp.Poly(sp.expand(k * (1 + sp.I * k) ** (ma + 1) * (1 - sp.I * k) ** (mb + 1)), k)
                co = ca * cb * sp.factorial(ma) * sp.factorial(mb)
                for (j,), a_ in num.terms():
                    assert 2 * M > j + 1
                    tot += co * a_ * sp.gamma(R(j + 1, 2)) * sp.gamma(M - R(j + 1, 2)) / sp.gamma(M) / 2
        return sp.simplify(sp.expand(tot))

    ns = range(2, 6)
    worst, worst_pl, wsamp = 0, 0, None
    for n1 in ns:
        for n2 in ns:
            wq = form(Qc(n1), Qc(n2))
            wf = form(Fc(n1), Fc(n2))
            wp = form(Fc(n1, 1), Fc(n2, 1))
            d, dp = sp.simplify(wq - wf), sp.simplify(wq - wp)
            if d != 0:
                worst = 1
            if dp != 0:
                worst_pl = 1
            if (n1, n2) == (2, 3):
                wsamp = (wq, wf, wp)
    print(f"    sample (n1,n2)=(2,3): form on Q = {wsamp[0]}, on F = {wsamp[1]}, planted = {wsamp[2]}")
    line("T6 half-line ultralocal form: W(F_a, F_b) = W(Q_a, Q_b), complex, n = 2..5 (16 entries, exact)",
         worst == 0, "all equal" if worst == 0 else "differs")
    line("P2 plant: coefficient 1 in place of 2 in the transform", worst_pl == 0,
         "equal" if worst_pl == 0 else "differs", plant=True)
    # the integral constraint on the edge-anchored potential
    ok = all(sum(cm * sp.factorial(m) for m, cm in Fc(n).items()) == 0 for n in ns)
    line("T6b int_0^inf F du = 0 for every F = T^-1 Q (the constraint that singles out curvature observables)",
         ok, "exact")


# ------------------------------------------------------------------ T7
def t7():
    u, n = sp.symbols('u n', positive=True)
    T = lambda a: (a - 1) / (a + 1)          # Q = S - (2/u) int_0^u S on u^a
    S = u ** n
    # vacuum flow: N d_u on Q, N = u - u^2, pulled back to S
    dS_vac = T(n) * (n * u ** n / T(n) - n * u ** (n + 1) / T(n + 1))
    dS_lit = (u - u ** 2) * sp.diff(S, u)
    diff = sp.simplify(dS_vac - dS_lit - 2 * sp.integrate(S, (u, 0, u)))
    line("T7a vacuum flow on S = u(1-u) S' + 2 int_0^u S  (all monomials u^n)", diff == 0, f"residual {diff}")
    psi0 = lambda expr: -sp.diff(expr, u, 2) / (2 * u)
    loc = sp.simplify(psi0(dS_vac) - ((u - u ** 2) * sp.diff(psi0(S), u) + (3 - 5 * u) * psi0(S)))
    line("T7b on Psi0 the vacuum flow is local: delta Psi0 = N Psi0' + (3 - 5u) Psi0", loc == 0, f"residual {loc}")
    lit = sp.simplify(psi0(dS_lit) - ((u - u ** 2) * sp.diff(psi0(S), u) + (3 - 5 * u) * psi0(S)))
    line("P7 plant: the literal weight-0 flow on S is the same local flow on Psi0", lit == 0,
         f"residual {sp.factor(lit)}  (= -S'/u, a non-local functional of Psi0)", plant=True)
    # the literal flow does not preserve int F du = 0
    F = sp.Function('F')
    line("T7c literal flow moves the constraint: int N F' du = 2 int u F du (by parts), not 0", True,
         "by hand: int_0^1 u(1-u) F' du = -int (1-2u) F du = 2 int u F du when int F du = 0")


if __name__ == '__main__':
    t_ = time.time()
    t0(); t5(); t6(); t7()
    bad = sum(1 for r_ in RESULTS if r_[1] in ('FAIL', 'MISSED'))
    print(f"\n{len(RESULTS)} lines, {bad} unexpected, {time.time() - t_:.1f} s")
    with open(sys.argv[1] if len(sys.argv) > 1 else 'shear_curv.json', 'w') as fh:
        json.dump([list(r_) for r_ in RESULTS], fh, indent=1)
