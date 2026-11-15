"""SS-LEMMA script B: the moments straight from the definitions (independent of the hand reduction).
Source: point curvature observable, f^{bd} = d_e d_f (Phi^{bedf} delta_p), Phi with Riemann symmetries.
Field:  h_bd = Phi_{b}^{e}{}_{d}^{f} d_e d_f G(x - p)   (trace part drops in R(l,m,l,m): eta(l,l)=eta(l,m)=eta(m,m)=0).
Psi_0 = R(l,m,l,m)[h] = -(1/2) [ (l.d)^2 h(m,m) + (m.d)^2 h(l,l) - CROSS (l.d)(m.d) h(l,m) ],  CROSS = 2.
G(z) = g(sigma), sigma = z.z:  D = 4: delta(sigma);  D = 6: delta'(sigma);  D = 5: (-sigma)_+^(-3/2) (finite part).
Any derivative of g(sigma) is sum_j P_j(z) g^(j)(sigma), d_a g^(j) = 2 eta_aa z_a g^(j+1).
On the generator x = u l: sigma(u) = p.p - 2 u l.p = b (u - u_x), b = -2 l.p > 0, u_x = -p.p/b.
Overall positive/negative constants of G are dropped (only zero / non-zero and ratios are tested).
"""
import itertools, json, os, random, sys
import sympy as sp

R = sp.Rational
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'out', 'lemma_direct.json')
uu = sp.symbols('u')
res = {}


def eta(D):
    return [-1] + [1]*(D - 1)


def dot(a, b):
    e = eta(len(a))
    return sum(e[i]*a[i]*b[i] for i in range(len(a)))


def deriv_tensor(D, order):
    """dict: sorted index tuple of length `order` -> {j: polynomial in z} for d_{a1}..d_{an} g(sigma)."""
    zs = sp.symbols('z0:%d' % D)
    e = eta(D)
    cur = {(): {0: sp.Integer(1)}}
    for _ in range(order):
        nxt = {}
        for idx, fd in cur.items():
            for a in range(D):
                key = tuple(sorted(idx + (a,)))
                if key in nxt:
                    continue
                out = {}
                for j, P in fd.items():
                    dP = sp.diff(P, zs[a])
                    if dP != 0:
                        out[j] = out.get(j, 0) + dP
                    out[j + 1] = out.get(j + 1, 0) + 2*e[a]*zs[a]*P
                nxt[key] = {j: sp.expand(P) for j, P in out.items()}
        cur = nxt
    return zs, cur


def moment_table(D, order, l, p, ks, kind):
    """moments int u^k [d^order g](u l - p) du for every sorted index tuple; exact."""
    zs, T = deriv_tensor(D, order)
    b = -2*dot(l, p)
    ux = -dot(p, p)/b
    sub = {zs[a]: uu*l[a] - p[a] for a in range(D)}
    table = {}
    for key, fd in T.items():
        polys = {j: sp.Poly(sp.expand(P.subs(sub)), uu) for j, P in fd.items()}
        for k in ks:
            tot = 0
            for j, P in polys.items():
                if kind in ('delta', 'ddelta'):
                    jj = j + (1 if kind == 'ddelta' else 0)
                    h = sp.Poly(uu**k, uu)*P
                    tot += (-1)**jj*b**(-jj - 1)*h.diff((uu, jj)).as_expr().subs(uu, ux)
                else:  # D = 5: g = (-sigma)_+^lam, lam = -3/2; result is in units of s^lam, s = b u_x
                    lam = R(-3, 2)
                    cj = (-1)**j*sp.prod([lam - i for i in range(j)])
                    for (n,), c in P.terms():
                        nn = n + k
                        beta = sp.gamma(nn + 1)*sp.gamma(lam - j + 1)/sp.gamma(nn + lam - j + 2)
                        tot += cj*c*b**(-j)*ux**(nn - j + 1)*sp.simplify(beta)
            table[(key, k)] = tot
    return table, ux


def bivector_basis(D):
    return [(i, j) for i in range(D) for j in range(i + 1, D)]


def phi_basis(D):
    """Phi^{bedf} = E_I^{be} E_J^{df} + E_J^{be} E_I^{df}, I <= J, as sparse dicts."""
    bv = bivector_basis(D)
    out = []
    for a, I in enumerate(bv):
        for J in bv[a:]:
            Phi = {}
            for (X, Y) in ((I, J), (J, I)):
                for (b_, e_, s1) in ((X[0], X[1], 1), (X[1], X[0], -1)):
                    for (d_, f_, s2) in ((Y[0], Y[1], 1), (Y[1], Y[0], -1)):
                        Phi[(b_, e_, d_, f_)] = Phi.get((b_, e_, d_, f_), 0) + s1*s2
            out.append(((I, J), Phi))
    return out


def psi0_moments_curv(D, l, m, Phi, table, ks, cross=2):
    """M_k for the curvature point source Phi. Operator: C^{a c e f} d_a d_c d_e d_f G."""
    e = eta(D)
    ll = [e[i]*l[i] for i in range(D)]
    ml = [e[i]*m[i] for i in range(D)]
    # Q_XY^{ef} = X_b Y_d Phi^{bedf}
    def Q(X, Y):
        q = {}
        for (b_, e_, d_, f_), c in Phi.items():
            val = c*X[b_]*Y[d_]
            if val != 0:
                q[(e_, f_)] = q.get((e_, f_), 0) + val
        return q
    Qmm, Qll, Qlm = Q(ml, ml), Q(ll, ll), Q(ll, ml)
    out = {}
    for k in ks:
        tot = 0
        for a in range(D):
            for c in range(D):
                for (pref, q) in ((l[a]*l[c], Qmm), (m[a]*m[c], Qll), (-cross*l[a]*m[c], Qlm)):
                    if pref == 0:
                        continue
                    for (e_, f_), val in q.items():
                        tot += pref*val*table[(tuple(sorted((a, c, e_, f_))), k)]
        out[k] = -sp.sympify(tot)/2   # sympify: an empty sum is the Python int 0, and -0/2 would be the float 0.0 (run01 miscounted these)
    return out


def is_zero(x):
    return sp.simplify(sp.cancel(sp.together(sp.expand(x)))) == 0


ok = True
ks = (0, 1, 2, 3, 4)

# ---------- B1, B2, B5: D = 4, canonical direction, symbolic p ----------
ps = sp.symbols('p0:4', real=True)
l4 = [sp.Integer(1), 0, 0, sp.Integer(1)]
m4 = [0, sp.Integer(1), sp.I, 0]
table, ux = moment_table(4, 4, l4, ps, ks, 'delta')
basis4 = phi_basis(4)
nz = {k: 0 for k in ks}
nz_plant = {k: 0 for k in ks}
for (IJ, Phi) in basis4:
    M = psi0_moments_curv(4, l4, m4, Phi, table, ks)
    Mp = psi0_moments_curv(4, l4, m4, Phi, table, ks, cross=1)
    for k in ks:
        nz[k] += 0 if is_zero(M[k]) else 1
        nz_plant[k] += 0 if is_zero(Mp[k]) else 1
print('B1 D = 4 symbolic p, 21 basis sources: number with M_k not identically zero:', {k: nz[k] for k in (0, 1, 2, 3)})
b1 = all(nz[k] == 0 for k in (0, 1, 2, 3))
print('B2 teeth: sources with M_4 not identically zero: %d of %d' % (nz[4], len(basis4)))
b2 = nz[4] > 0
b5 = any(nz_plant[k] > 0 for k in (0, 1, 2, 3))
print('B5 PLANT (cross term -1 instead of -2): sources with M_k != 0:', {k: nz_plant[k] for k in (0, 1, 2, 3)}, '-> plant caught:', b5)
res.update(B1=b1, B1_counts={str(k): nz[k] for k in ks}, B2=b2, B5_plant_caught=b5, B5_counts={str(k): nz_plant[k] for k in ks})

# ---------- B3: generic rational direction, random rational interior points ----------
lg = [sp.Integer(1), R(2, 3), R(-1, 3), R(2, 3)]
mg = [0, R(1, 3) + sp.I*R(2, 3), R(-2, 3) + sp.I*R(2, 3), R(-2, 3) - sp.I*R(1, 3)]
assert sp.expand(dot(lg, lg)) == 0 and sp.expand(dot(mg, mg)) == 0 and sp.expand(dot(lg, mg)) == 0
random.seed(7312)
pts = []
while len(pts) < 4:
    t = R(random.randint(3, 17), 10)
    xv = [R(random.randint(-9, 9), 20) for _ in range(3)]
    r2 = sum(c**2 for c in xv)
    if r2 < t**2 and r2 < (2 - t)**2 and (t < 1 or True):
        if sp.sqrt(r2) < min(t, 2 - t):
            pts.append([t] + xv)
b3_bad = 0
for p in pts:
    tb, uxg = moment_table(4, 4, lg, p, (0, 1, 2, 3), 'delta')
    assert 0 < uxg < 1
    for (IJ, Phi) in basis4:
        M = psi0_moments_curv(4, lg, mg, Phi, tb, (0, 1, 2, 3))
        b3_bad += sum(0 if sp.simplify(M[k]) == 0 else 1 for k in (0, 1, 2, 3))
print('B3 generic direction, 4 random interior p x 21 sources x 4 moments: non-zero count = %d of %d' % (b3_bad, 4*21*4))
b3 = (b3_bad == 0)
res['B3'] = b3

# ---------- B4: non-conserved point source h_bd = Phi_bd G ----------
p = [R(9, 10), R(1, 5), R(-3, 20), R(1, 4)]
tb2, ux2 = moment_table(4, 2, lg, p, ks, 'delta')
Ph = [[R(3), R(1, 2), R(-2), R(1)], [R(1, 2), R(-1), R(4, 3), R(2)], [R(-2), R(4, 3), R(5), R(-1, 7)], [R(1), R(2), R(-1, 7), R(2, 5)]]
e4 = eta(4)
llg = [e4[i]*lg[i] for i in range(4)]; mlg = [e4[i]*mg[i] for i in range(4)]
hmm = sum(Ph[b_][d_]*mg[b_]*mg[d_] for b_ in range(4) for d_ in range(4))   # Phi with lower indices contracted with m^b m^d
hll = sum(Ph[b_][d_]*lg[b_]*lg[d_] for b_ in range(4) for d_ in range(4))
hlm = sum(Ph[b_][d_]*lg[b_]*mg[d_] for b_ in range(4) for d_ in range(4))
M = {}
for k in ks:
    tot = 0
    for a in range(4):
        for c in range(4):
            tot += (lg[a]*lg[c]*hmm + mg[a]*mg[c]*hll - 2*lg[a]*mg[c]*hlm)*tb2[(tuple(sorted((a, c))), k)]
    M[k] = sp.simplify(-tot/2)
ratio = sp.simplify(M[3]/M[2]) if M[2] != 0 else None
print('B4 non-conserved source: M_0 = %s, M_1 = %s, M_2 = %s, M_3 = %s; M_3/M_2 = %s, 3 u_x = %s' % (M[0], M[1], M[2], M[3], ratio, 3*ux2))
b4 = (M[0] == 0 and M[1] == 0 and M[2] != 0 and M[3] != 0 and sp.simplify(ratio - 3*ux2) == 0)
res['B4'] = bool(b4)

# ---------- B6, B7: D = 6 and D = 5 controls ----------
def control(D, kind):
    l = [sp.Integer(1), 0, 0, sp.Integer(1)] + [0]*(D - 4)
    m = [0, sp.Integer(1), sp.I, 0] + [0]*(D - 4)
    extra = [R(1, 6), R(-1, 9)][:D - 4]
    points = [[R(9, 10), R(1, 5), R(-3, 20), R(1, 4)] + extra, [R(6, 5), R(-1, 4), R(1, 3), R(1, 10)] + extra]
    counts = {k: 0 for k in (0, 1, 2, 3)}
    nsrc = 0
    for p in points:
        tb, _ = moment_table(D, 4, l, p, (0, 1, 2, 3), kind)
        for (IJ, Phi) in phi_basis(D):
            Mk = psi0_moments_curv(D, l, m, Phi, tb, (0, 1, 2, 3))
            nsrc += 1
            for k in (0, 1, 2, 3):
                counts[k] += 0 if sp.simplify(Mk[k]) == 0 else 1
    return counts, nsrc

c6, n6 = control(6, 'ddelta')
print('B6 D = 6 control: curvature sources with M_k != 0 (of %d):' % n6, c6)
b6 = any(v > 0 for v in c6.values())
c5, n5 = control(5, 'frac')
print('B7 D = 5 control: curvature sources with M_k != 0 (of %d):' % n5, c5)
b7 = any(v > 0 for v in c5.values())
res.update(B6=b6, B6_counts={str(k): v for k, v in c6.items()}, B7=b7, B7_counts={str(k): v for k, v in c5.items()})

ok = b1 and b2 and b3 and b4 and b5 and b6 and b7
res['PASS'] = bool(ok)
with open(OUT, 'w', encoding='utf-8') as fh:
    json.dump(res, fh, indent=1)
for name in ('B1', 'B2', 'B3', 'B4', 'B5_plant_caught', 'B6', 'B7'):
    print('  %-16s %s' % (name, 'PASS' if res[name] else 'FAIL'))
print('SCRIPT B', 'PASS' if ok else 'FAIL')
sys.exit(0 if ok else 1)
