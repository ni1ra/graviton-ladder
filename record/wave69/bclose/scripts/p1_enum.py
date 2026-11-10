"""PART 1 of w69-bclose: enumerate every term of the second-order modular side that is linear in h2,
and assert each is the smeared constraint delta(A - kV).  The assertion is the negation of clause (3)'s
failure route, so it can fail.

  python3 p1_enum.py --plants    plants only   -> out/p1_plants.json   (run first)
  python3 p1_enum.py             full run      -> out/p1_enum.json

Method.  h2_{mu nu}(t,x) = a_{mu nu}(x) + t b_{mu nu}(x) + t^2 c_{mu nu}(x)/2 with a, b, c ARBITRARY
functions (sympy Functions).  Every building block of the generator is linearised at eta from its nonlinear
definition and reduced to a list of (function slot, derivative, polynomial coefficient).  Units 16 pi G = 1;
the total is Paper A's "A".  All integrals are exact rationals times pi.

Reused: the P-tensor of wave66/p0test/scripts/s1_transfer.py (S1c), the Green identity of
wave64/btheorem (Lemma 1), the monomial integrals of wave68/q0edge/qcommon.py.
"""
import sys
from fractions import Fraction
import sympy as sp
from sympy.core.function import AppliedUndef
from common import *

PLANTS_ONLY = '--plants' in sys.argv
out = {}
half = sp.Rational(1, 2)
eta = sp.diag(-1, 1, 1, 1)
e_ = [-1, 1, 1, 1]


def fn(kind, m, n):
    return sp.Function(f'{kind}{min(m, n)}{max(m, n)}')(x, y, z)


Am = sp.Matrix(4, 4, lambda m, n: fn('a', m, n))
Bm = sp.Matrix(4, 4, lambda m, n: fn('b', m, n))
Cm = sp.Matrix(4, 4, lambda m, n: fn('c', m, n))
h = Am + t*Bm + t**2/2*Cm


def d4(f, m):
    return sp.diff(f, X4[m])


# ------------------------------------------------------------------ linearisation at eta
hup = sp.Matrix(4, 4, lambda m, n: e_[m]*e_[n]*h[m, n])
trh = sum(e_[m]*h[m, m] for m in range(4))
Gam = [[[half*e_[r]*(d4(h[r, n], m) + d4(h[r, m], n) - d4(h[m, n], r)) for n in range(4)] for m in range(4)]
       for r in range(4)]                                   # Gam[r][m][n] = Gamma^r_{mn}, first order
zeta = [N - t**2/2, -t*x, -t*y, -t*z]
dz = [[d4(zeta[n], a) for n in range(4)] for a in range(4)]             # dz[a][n] = d_a zeta^n
# nabla^m zeta^n = g^{ma}(d_a zeta^n + Gamma^n_{ab} zeta^b): orders 0 and 1
up0 = [[e_[m]*dz[m][n] for n in range(4)] for m in range(4)]
up1 = [[-sum(hup[m, a]*dz[a][n] for a in range(4)) + e_[m]*sum(Gam[n][m][b]*zeta[b] for b in range(4))
        for n in range(4)] for m in range(4)]
W0 = [[up0[m][n] - up0[n][m] for n in range(4)] for m in range(4)]
W1 = [[up1[m][n] - up1[n][m] + trh/2*W0[m][n] for n in range(4)] for m in range(4)]   # delta(sqrt(-g) W)
Qi = [W1[i][0] for i in (1, 2, 3)]                 # delta[ sqrt(-g) (nabla^i zeta^t - nabla^t zeta^i) ]

Ric = [[sum(d4(Gam[r][m][n], r) - d4(Gam[r][m][r], n) for r in range(4)) for n in range(4)] for m in range(4)]
Rs = sum(e_[m]*Ric[m][m] for m in range(4))
G1 = [[Ric[m][n] - half*eta[m, n]*Rs for n in range(4)] for m in range(4)]
theta = [sum(d4(hup[a, b], b) for b in range(4)) - e_[a]*d4(trh, a) for a in range(4)]


def P(a, b, c, d, e, f):
    dl = lambda i, j: e_[i] if i == j else 0
    return (dl(a, e)*dl(f, b)*dl(c, d) - half*dl(a, d)*dl(b, e)*dl(f, c) - half*dl(a, b)*dl(c, d)*dl(e, f)
            - half*dl(b, c)*dl(a, e)*dl(f, d) + half*dl(b, c)*dl(a, d)*dl(e, f))


gam = -2*t*eta                                      # L_zeta eta


def omega_t(g1, g2):
    tot = 0
    for b in range(4):
        for c in range(4):
            for d in range(4):
                for e in range(4):
                    for f in range(4):
                        p = P(0, b, c, d, e, f)
                        if p != 0:
                            tot += p*(g2[b, c]*d4(g1[e, f], d) - g1[b, c]*d4(g2[e, f], d))
    return tot


def at0(e):
    return sp.expand(e.subs(t, 0))


q = [[Am[i + 1, j + 1] for j in range(3)] for i in range(3)]
qkk = q[0][0] + q[1][1] + q[2][2]
qxx = sum(X[i]*X[j]*q[i][j] for i in range(3) for j in range(3))
Rlin_q = sum(sp.diff(q[i][j], X[i], X[j]) for i in range(3) for j in range(3)) - sum(sp.diff(qkk, s, 2) for s in X)

# ------------------------------------------------------------------ the rows, as densities at t = 0
d_R1 = at0(sum(X[i]*Qi[i] for i in range(3)))                       # on S
d_R4 = at0(2*e_[0]*sum(G1[0][b]*zeta[b] for b in range(4)))         # in B: 2 delta(sqrt(-g) G^t_b zeta^b)
d_om = at0(omega_t(h, gam))                                         # in B
d_zth = [at0(zeta[i + 1]*theta[0] - zeta[0]*theta[i + 1]) for i in range(3)]
divQ = at0(sum(sp.diff(Qi[i], X[i]) for i in range(3)))
div_zth = sp.expand(sum(sp.diff(d_zth[i], X[i]) for i in range(3)))

# ------------------------------------------------------------------ gate 1: Iyer-Wald identity, pointwise
gate = {}
for s1 in (1, -1):
    for s2 in (1, -1):
        gate[(s1, s2)] = sp.expand(divQ + d_R4 - s1*d_om - s2*div_zth) == 0
sol = [k for k, v in gate.items() if v]
out['gate_IW_sign_solutions (s_omega, s_theta)'] = str(sol)
out['CHECK gate IW: delta(dQ + C) = s1 omega + s2 d(i_zeta theta) pointwise, unique signs'] = (len(sol) == 1)
s1, s2 = sol[0] if sol else (1, 1)
d_R3 = sp.expand(s1*d_om)
d_R2 = sp.expand(s2*sum(X[i]*d_zth[i] for i in range(3)))            # on S
out['CHECK G_tt(lin) = R_lin[q]/2 pointwise'] = sp.expand(at0(G1[0][0]) - Rlin_q/2) == 0
out['momentum-constraint weight zeta^i at t=0'] = str([at0(zeta[i]) for i in (1, 2, 3)])
out['CHECK R4 momentum weight is zero'] = all(at0(zeta[i]) == 0 for i in (1, 2, 3))

# gate 2: Green identity as a pointwise divergence (Lemma 1 of the B note)
Jv = [N*(sum(sp.diff(q[i][j], X[j]) for j in range(3)) - sp.diff(qkk, X[i]))
      - sum(q[i][j]*sp.diff(N, X[j]) for j in range(3)) + qkk*sp.diff(N, X[i]) for i in range(3)]
out['CHECK gate Green: N R_lin - 2 q_kk = div J pointwise'] = sp.expand(
    sum(sp.diff(Jv[i], X[i]) for i in range(3)) - (N*Rlin_q - 2*qkk)) == 0
fluxJ = sp.expand(sum(X[i]*Jv[i] for i in range(3)))


def on_S_zero(expr):
    """True iff expr (linear in function atoms, polynomial coefficients) vanishes on r = 1."""
    expr = sp.expand(expr)
    ats = sorted(expr.atoms(sp.Derivative) | expr.atoms(AppliedUndef), key=str)
    if not ats:
        return sp.rem(expr, r2 - 1, z) == 0
    Pp = sp.Poly(expr, *ats)
    return all(sp.rem(c, r2 - 1, z) == 0 for c in Pp.coeffs())


def atoms_on_S(expr):
    expr = sp.expand(expr)
    ats = sorted(expr.atoms(sp.Derivative) | expr.atoms(AppliedUndef), key=str)
    Pp = sp.Poly(expr, *ats)
    keep = []
    for mon, c in zip(Pp.monoms(), Pp.coeffs()):
        if sp.rem(c, r2 - 1, z) != 0:
            keep.append(str(ats[mon.index(1)]))
    return keep


out['CHECK gate Green flux: n.J = q_nn - q_kk on S'] = on_S_zero(fluxJ - (qxx - qkk))

# ------------------------------------------------------------------ O1: symbolic residuals
res = {}
res['R1 corner Noether charge'] = on_S_zero(d_R1 + (qkk - qxx))             # = -2 dA density
res['R2 corner potential i_zeta theta'] = on_S_zero(d_R2)
res['R3 bulk non-Killing omega'] = sp.expand(d_R3 + 2*qkk) == 0              # = -4 dV density
res['R4 bulk constraint'] = sp.expand(d_R4 + N*Rlin_q) == 0                  # = -C density
out['R1 atoms surviving on S'] = atoms_on_S(d_R1)
out['R3 density'] = str(d_R3)
out['R1 raw atom count before restriction to S'] = len(sp.expand(d_R1).atoms(sp.Derivative) | sp.expand(d_R1).atoms(AppliedUndef))

# geometric rows (Z2), from the nonlinear area and volume elements
s = sp.Symbol('s')
qs = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f'q{min(i, j)}{max(i, j)}'))
gam3 = sp.eye(3) + s*qs
xv = sp.Matrix(X)
area_el = sp.sqrt(gam3.det()*(xv.T*gam3.inv()*xv)[0])                        # sqrt(gamma) |grad r|_gamma at r = 1
vol_el = sp.sqrt(gam3.det())
dA_lin = sp.expand(sp.diff(area_el, s).subs(s, 0).subs(sp.sqrt(r2), 1))
dV_lin = sp.expand(sp.diff(vol_el, s).subs(s, 0))
qs_kk = qs[0, 0] + qs[1, 1] + qs[2, 2]
qs_xx = sp.expand((xv.T*qs*xv)[0])
tgt = sp.expand(sp.together(dA_lin - half*(qs_kk - qs_xx/r2)*sp.sqrt(r2)))
res['R5 geometric area (Z2)'] = sp.simplify(sp.diff(area_el, s).subs(s, 0) - half*(qs_kk*sp.sqrt(r2) - qs_xx/sp.sqrt(r2))) == 0
res['R6 geometric volume (Z2)'] = sp.expand(dV_lin - half*qs_kk) == 0

# R7: radial relocation.  Generator potentials on a sphere of radius R in eta, units 16 pi G = 1.
Rr = sp.Symbol('R', positive=True)
komar_dens = sum(X[i]*W0[i + 1][0] for i in range(3)).subs(t, 0)                 # x_i (nabla^i zeta^t - nabla^t zeta^i)
komar_R = sp.simplify((komar_dens/Rr).subs({x: Rr, y: 0, z: 0})*4*sp.pi*Rr**2)
G_Z1 = komar_R + 4*sp.Rational(4, 3)*sp.pi*Rr**3
G_Z2 = -2*4*sp.pi*Rr**2 + 4*sp.Rational(4, 3)*sp.pi*Rr**3
dG1 = sp.diff(G_Z1, Rr).subs(Rr, 1)
dG2 = sp.diff(G_Z2, Rr).subs(Rr, 1)
out['R7 Komar charge of fixed zeta on S_R in eta'] = str(komar_R)
out['R7 dG/dR at R=1: Z1 (fixed zeta)'] = str(dG1)
out['R7 dG/dR at R=1: Z2 (A - 2V)'] = str(dG2)
out['CHECK normalisation: G_Z1(1) = G_Z2(1)'] = sp.simplify(G_Z1.subs(Rr, 1) - G_Z2.subs(Rr, 1)) == 0
res['R7 radial relocation (Z2): coefficient'] = (dG2 == 0)
out['O3 coefficient in "8 pi G H - 2V" units (dG1 / -2)'] = str(dG1/(-2))

# R8: time relocation.  (i) slice t = s T(x) in eta has induced metric delta - s^2 dT dT: no O(s) term.
Tf = sp.Function('T')(x, y, z)
ind = sp.Matrix(3, 3, lambda i, j: (1 if i == j else 0) - s**2*sp.diff(Tf, X[i])*sp.diff(Tf, X[j]))
res['R8 time relocation: area and volume (Z2)'] = all(sp.diff(ind[i, j], s).subs(s, 0) == 0 for i in range(3) for j in range(3))
# (ii) Komar flux of the fixed zeta through the cylinder r = 1: J^a = theta^a(eta; L_zeta eta) - zeta^a L, L = 0.
th_bg = [sum(d4(e_[a]*e_[b]*gam[a, b], b) for b in range(4)) - e_[a]*d4(sum(e_[m]*gam[m, m] for m in range(4)), a)
         for a in range(4)]
out['R8 J^a of zeta in eta'] = str(th_bg)
res['R8 time relocation: Komar flux through the cylinder (Z1)'] = all(th_bg[i] == 0 for i in (1, 2, 3))
out['CHECK R7 cross-check: dKomar/dR = J^t * 4 pi'] = sp.simplify(sp.diff(komar_R, Rr).subs(Rr, 1) - th_bg[0]*4*sp.pi) == 0

# linearised mean curvature of r = 1 in delta + q (level function (r^2-1)/2, so grad = x):
#   dH = (1/2) x.grad(q_kk) - x_j d_i q_ij - q_kk + q_nn + (1/2) x.grad(q_nn)   on r = 1
d_dH = (half*sum(X[i]*sp.diff(qkk, X[i]) for i in range(3)) - sum(X[j]*sp.diff(q[i][j], X[i]) for i in range(3) for j in range(3))
        - qkk + qxx + half*sum(X[i]*sp.diff(qxx, X[i]) for i in range(3)))


def check_dH_formula():
    """Independent derivation at a rational point of S, for a polynomial q, from H = div_g(n_g)."""
    qp = [[x*y + z**2, 2*x*z - y, y**2], [0, x**2 - y*z + 1, x*y + z], [0, 0, 3*z**2 + x]]
    qp = [[qp[min(i, j)][max(i, j)] for j in range(3)] for i in range(3)]
    g = sp.eye(3) + s*sp.Matrix(qp)
    gi = sp.eye(3) - s*sp.Matrix(qp)
    w = (xv.T*gi*xv)[0]
    sg = 1 + s*sum(qp[i][i] for i in range(3))/2
    nvec = [sg*sum(gi[i, j]*X[j] for j in range(3))/sp.sqrt(w) for i in range(3)]
    Hs = sum(sp.diff(nvec[i], X[i]) for i in range(3))/sg
    lhs = sp.diff(Hs, s).subs(s, 0).subs({x: sp.Rational(2, 3), y: sp.Rational(1, 3), z: sp.Rational(2, 3)})
    sub = {}
    for i in range(3):
        for j in range(i, 3):
            sub[q[i][j]] = qp[i][j]
    rhs = d_dH.subs(sub).doit().subs({x: sp.Rational(2, 3), y: sp.Rational(1, 3), z: sp.Rational(2, 3)})
    return sp.simplify(lhs - rhs) == 0


out['CHECK linearised mean-curvature formula against div_g(n_g)'] = check_dH_formula()

# ------------------------------------------------------------------ rows as exact linear functionals
def to_row(expr, dom):
    expr = sp.expand(expr)
    ats = sorted(expr.atoms(sp.Derivative) | expr.atoms(AppliedUndef), key=str)
    if not ats:
        assert expr == 0
        return (dom, [])
    Pp = sp.Poly(expr, *ats)
    terms = []
    for mon, c in zip(Pp.monoms(), Pp.coeffs()):
        assert sum(mon) == 1, 'row is not linear in h2'
        at = ats[mon.index(1)]
        dv = [0, 0, 0]
        f = at
        if isinstance(at, sp.Derivative):
            f = at.expr
            for v, n in at.variable_count:
                dv[X.index(v)] += int(n)
        nm = f.func.__name__
        terms.append(((nm[0], int(nm[1]), int(nm[2])), tuple(dv), pd(c)))
    return (dom, terms)


def ev(row, data):
    dom, terms = row
    tot = Fraction(0)
    for slot, dv, cf in terms:
        if slot in data:
            tot += pd_int_prod(cf, pd_diff(data[slot], dv), dom)
    return tot


def subs_q(expr):
    return expr.subs({qs[i, j]: q[i][j] for i in range(3) for j in range(i, 3)})


ROW = {
    'dA': to_row(half*(qkk - qxx), 'S'),
    'dV': to_row(half*qkk, 'B'),
    'C': to_row(N*Rlin_q, 'B'),
    'R1': to_row(d_R1, 'S'),
    'R2': to_row(d_R2, 'S'),
    'R3': to_row(d_R3, 'B'),
    'R4': to_row(d_R4, 'B'),
    'R5': to_row(subs_q(sp.diff(area_el, s).subs(s, 0).subs(sp.sqrt(r2), 1).subs(r2, 1)), 'S'),
    'R6': to_row(subs_q(-2*dV_lin), 'B'),
    'dH': to_row(d_dH, 'S'),
}
f_zon = zonal(2)
f_sec = x**2 - y**2
Rij = [[half*(sum(sp.diff(q[k][j], X[k], X[i]) + sp.diff(q[k][i], X[k], X[j]) for k in range(3))
              - sum(sp.diff(q[i][j], v, 2) for v in X) - sp.diff(qkk, X[i], X[j])) for j in range(3)] for i in range(3)]
Rnn = sum(X[i]*X[j]*Rij[i][j] for i in range(3) for j in range(3))
ROW['KLS_zon'] = to_row(half*f_zon*(qkk - qxx), 'S')
ROW['KLS_sec'] = to_row(half*f_sec*(qkk - qxx), 'S')
ROW['CURV_zon'] = to_row(f_zon*Rnn, 'S')
ROW['CURV_sec'] = to_row(f_sec*Rnn, 'S')
ROW['R4_lapse'] = to_row(-(N + (1 - r2)**2)*Rlin_q, 'B')


def comb(*pairs):
    """A linear combination of rows as a callable on data."""
    return lambda data: sum((Fraction(c)*ev(ROW[k], data) for c, k in pairs), Fraction(0))


dG1_pi = Fraction(int(dG1/sp.pi))                    # dG_Z1/dR / pi
R7_Z1 = lambda data: dG1_pi/8*ev(ROW['dH'], data)    # (dG/dR) * (1/8 pi) oint dH

# ------------------------------------------------------------------ the probe space: all slots, degree <= 3
MONO = [(i, j, k) for i in range(4) for j in range(4) for k in range(4) if i + j + k <= 3]
SLOTS = [(kd, m, n) for kd in 'abc' for m in range(4) for n in range(m, 4)]
BASIS = [(sl, mo) for sl in SLOTS for mo in MONO]


def vec(F):
    return [F({sl: {mo: Fraction(1)}}) for sl, mo in BASIS]


V_A, V_V, V_C = vec(comb((1, 'dA'))), vec(comb((1, 'dV'))), vec(comb((1, 'C')))
out['CHECK Green identity on the probe space: C = -2 dA + 4 dV'] = all(c == -2*a + 4*b for a, b, c in zip(V_A, V_V, V_C))
M2 = sp.Matrix([[sp.Rational(v.numerator, v.denominator) for v in V_A], [sp.Rational(v.numerator, v.denominator) for v in V_V]])


def classify(F):
    """(in_span, a, b, n_bad): is F = a dA + b dV on every probe?  Exact."""
    v = vec(F)
    # two probes that separate dA and dV: q = delta (slot a11+a22+a33 const) is not a single slot, so solve by least squares
    Mv = sp.Matrix([[sp.Rational(u.numerator, u.denominator) for u in v]])
    G = M2*M2.T
    ab = G.LUsolve(M2*Mv.T)
    resid = Mv - ab.T*M2
    bad = [BASIS[i] for i in range(len(BASIS)) if resid[i] != 0]
    return (len(bad) == 0, ab[0], ab[1], bad)


def audit(name, rows_total, r4='R4', extra=None):
    """The pipeline.  rows_total: list of (label, callable) whose sum is the construction's h2-linear part.
    Returns a verdict dict; 'pass' is the negation of the failure route for this construction."""
    rec = {'rows': {}}
    ok = True
    a_tot = b_tot = 0
    for lab, F in rows_total:
        ins, a, b, bad = classify(F)
        rec['rows'][lab] = {'in span{dA,dV}': ins, 'a': str(a), 'b': str(b), 'n_probes_with_residual': len(bad),
                            'residual_slots': sorted({str(sl) for sl, _ in bad})[:12]}
        ok = ok and ins
        a_tot += a
        b_tot += b
    rec['total (a, b)'] = [str(a_tot), str(b_tot)]
    rec['k = -b/a'] = str(-b_tot/a_tot) if a_tot != 0 else 'undefined'
    k_ok = (a_tot != 0 and -b_tot/a_tot == 2)
    # constraint row: weight exactly N, and the integrated Iyer-Wald gate on every probe
    vR4 = vec(comb((1, r4)))
    c_ok = all(u == -c for u, c in zip(vR4, V_C))
    rec['constraint row equals -C on every probe'] = c_ok
    tot = lambda data: sum((F(data) for _, F in rows_total), Fraction(0))
    vT = vec(tot)
    gate_ok = all(u + w == 0 for u, w in zip(vT, vR4))
    rec['integrated identity total = -(constraint row) on every probe'] = gate_ok
    rec['pass'] = bool(ok and k_ok and c_ok and gate_ok)
    rec['_total'] = tot
    return rec


Z1L0 = [('R1', comb((1, 'R1'))), ('R2', comb((-1, 'R2'))), ('R3', comb((-1, 'R3')))]
Z2L1 = [('R5 (x -2)', comb((-2, 'R5'))), ('R6 (x -2, k=2 inside)', comb((-2, 'R6'))),
        ('R7 Z2 (coefficient 0)', lambda data: Fraction(int(dG2/sp.pi))/8*ev(ROW['dH'], data))]
Z1L1 = Z1L0 + [('R7 Z1', R7_Z1)]

# ------------------------------------------------------------------ on-shell data on the scope family
def tt_probe():
    """A source-free TT polynomial tensor of degree <= 3 with nonzero Laplacian (deterministic)."""
    cs, comp = [], {}
    for i in range(3):
        for j in range(i, 3):
            e = 0
            for mo in MONO:
                c = sp.Symbol(f'w{i}{j}_{mo[0]}{mo[1]}{mo[2]}')
                cs.append(c)
                e += c*x**mo[0]*y**mo[1]*z**mo[2]
            comp[(i, j)] = e
    T = [[comp[(min(i, j), max(i, j))] for j in range(3)] for i in range(3)]
    eqs = sp.Poly(T[0][0] + T[1][1] + T[2][2], x, y, z).coeffs()
    for j in range(3):
        eqs += sp.Poly(sum(sp.diff(T[i][j], X[i]) for i in range(3)), x, y, z).coeffs()
    Mx, _ = sp.linear_eq_to_matrix(eqs, cs)
    ns = Mx.nullspace()
    v = sum((((k*k % 7) + 1)*ns[k] for k in range(len(ns))), sp.zeros(len(cs), 1))   # weights changed once at the plant stage: see REPORT, Deviations
    sub = dict(zip(cs, v))
    return [[sp.expand(T[i][j].subs(sub)) for j in range(3)] for i in range(3)]


QTT = tt_probe()
out['CHECK TT probe: traceless, divergence-free, R_lin = 0, Laplacian nonzero'] = (
    sp.expand(QTT[0][0] + QTT[1][1] + QTT[2][2]) == 0
    and all(sp.expand(sum(sp.diff(QTT[i][j], X[i]) for i in range(3))) == 0 for j in range(3))
    and R_lin(QTT) == 0 and any(lap(QTT[i][j]) != 0 for i in range(3) for j in range(3)))
XI = [x/3 + y**2/2 + x*z, z*x - y/5 + x**2*y, 2*x/7 + z/4 - z**3/3 + y*z]
GQ = symgrad(XI)
UH = x**2 - y**2 + 3*x*z + z + sp.Rational(2, 3)
JUNK = {('a', 0, 0): pd(1 + x - 2*y*z + z**3), ('a', 0, 1): pd(x*y - 3), ('a', 0, 2): pd(z**2 + y), ('a', 0, 3): pd(2 - x**2*z),
        ('b', 0, 0): pd(x + 5), ('b', 0, 1): pd(y*z), ('b', 1, 1): pd(1 + x**2), ('b', 1, 2): pd(z - x*y), ('b', 2, 2): pd(3*y),
        ('b', 3, 3): pd(x*z - 2), ('b', 2, 3): pd(1 + y**3), ('c', 1, 1): pd(4 + z), ('c', 0, 2): pd(x)}


def data_of(qm, junk=False):
    dta = {('a', i + 1, j + 1): pd(qm[i][j]) for i in range(3) for j in range(i, 3)}
    if junk:
        dta.update(JUNK)
    return dta


GAUGE_PROBE = data_of(symgrad([x, -y, sp.Integer(0)]))        # l = 2 corner displacement, xi_n = x^2 - y^2
GAUGE_PROBE_Z = data_of(symgrad([-x/2, -y/2, z]))              # l = 2, m = 0 corner displacement
TT_DATA = data_of(QTT)
out['CHECK TT probe overlaps both plant profiles'] = (ev(ROW['CURV_zon'], TT_DATA) != 0 and ev(ROW['CURV_sec'], TT_DATA) != 0)


def scope_gauges(phi):
    K = magnetic_K(phi)
    k2 = K2(K)
    u = inv_lap(-k2/2)
    qa = [[u if i == j else sp.Integer(0) for j in range(3)] for i in range(3)]
    assert R_lin(qa) == k2
    qb = [[sp.expand(qa[i][j] + GQ[i][j]) for j in range(3)] for i in range(3)]
    qc = [[sp.expand(qb[i][j] + (UH if i == j else 0) + QTT[i][j]) for j in range(3)] for i in range(3)]
    assert R_lin(qb) == k2 and R_lin(qc) == k2
    return K, k2, {'G_a': data_of(qa), 'G_b': data_of(qb, True), 'G_c': data_of(qc, True),
                   'G_b_noTT': data_of([[sp.expand(qb[i][j] + (UH if i == j else 0)) for j in range(3)] for i in range(3)], True)}


SCOPE = {}
for name, l, phi in members():
    K, k2, gg = scope_gauges(phi)
    SCOPE[name] = (l, ball(N*k2), sphere(k2), gg)


def fr(v):
    return str(v)


def onshell(tot):
    """Per member: the total in the three gauges, as rationals (times pi)."""
    return {name: {g: tot(SCOPE[name][3][g]) for g in ('G_a', 'G_b', 'G_c', 'G_b_noTT')} for name in SCOPE}


# ------------------------------------------------------------------ PLANTS (run first)
plants = {}


def plant_report(label, rows_total, r4='R4', probe=None, probe_row=None):
    rec = audit(label, rows_total, r4)
    tot = rec.pop('_total')
    os_ = onshell(tot)
    per = {}
    for name in SCOPE:
        A_l = SCOPE[name][1]
        vals = os_[name]
        per[name] = {'A_l': fr(A_l), 'G_a': fr(vals['G_a']), 'G_b': fr(vals['G_b']), 'G_c': fr(vals['G_c']),
                     'gauge_dependent (G_a vs G_b)': vals['G_a'] != vals['G_b'],
                     'moves when the TT term is added (G_b_noTT vs G_c)': vals['G_b_noTT'] != vals['G_c'],
                     'equals A_l in all gauges': all(vals[g] == A_l for g in ('G_a', 'G_b', 'G_c'))}
    rec['on_shell'] = per
    rec['caught_by_audit'] = not rec['pass']
    rec['caught_on_shell_at_every_member'] = all(not p['equals A_l in all gauges'] for p in per.values())
    return rec


for tag in ('zon', 'sec'):
    key = 'KLS_' + tag
    rec = plant_report(key, Z1L0 + [(key, comb((1, key)))])
    rec['value on the pure-gauge corner displacement'] = fr(ev(ROW[key], GAUGE_PROBE_Z if tag == 'zon' else GAUGE_PROBE))
    rec['nonzero on pure gauge'] = ev(ROW[key], GAUGE_PROBE_Z if tag == 'zon' else GAUGE_PROBE) != 0
    rec['gauge dependent on shell at every member'] = all(p['gauge_dependent (G_a vs G_b)'] for p in rec['on_shell'].values())
    plants['PLANT-KLS f=' + tag] = rec
    out[f'CHECK PLANT-KLS f={tag} caught (audit, pure-gauge probe, on shell at every l)'] = (
        rec['caught_by_audit'] and rec['nonzero on pure gauge'] and rec['caught_on_shell_at_every_member']
        and rec['gauge dependent on shell at every member'])
    key = 'CURV_' + tag
    rec = plant_report(key, Z1L0 + [(key, comb((1, key)))])
    gp = [ev(ROW[key], GAUGE_PROBE), ev(ROW[key], GAUGE_PROBE_Z), ev(ROW[key], data_of(GQ))]
    rec['values on three pure-gauge probes'] = [fr(v) for v in gp]
    rec['invisible to the gauge channel'] = all(v == 0 for v in gp)
    rec['value on the source-free TT probe'] = fr(ev(ROW[key], TT_DATA))
    rec['moves with the TT term at every member'] = all(p['moves when the TT term is added (G_b_noTT vs G_c)'] for p in rec['on_shell'].values())
    plants['PLANT-CURV f=' + tag] = rec
    out[f'CHECK PLANT-CURV f={tag} caught (audit; zero on gauge probes; TT term moves the total at every l)'] = (
        rec['caught_by_audit'] and rec['invisible to the gauge channel'] and rec['moves with the TT term at every member']
        and rec['caught_on_shell_at_every_member'])

rec = audit('PLANT-LAPSE', Z1L0, r4='R4_lapse')
rec.pop('_total')
lap_on = {}
for name in SCOPE:
    v = -ev(ROW['R4_lapse'], SCOPE[name][3]['G_a'])
    lap_on[name] = {'-R4_lapse on shell': fr(v), 'A_l': fr(SCOPE[name][1]), 'differs': v != SCOPE[name][1]}
rec['on_shell'] = lap_on
rec['caught_by_audit'] = not rec['pass']
plants['PLANT-LAPSE'] = rec
out['CHECK PLANT-LAPSE caught (constraint row is not -C; gate fails; on-shell value differs at every l)'] = (
    rec['caught_by_audit'] and not rec['constraint row equals -C on every probe']
    and not rec['integrated identity total = -(constraint row) on every probe'] and all(v['differs'] for v in lap_on.values()))
out['plants'] = plants

if PLANTS_ONLY:
    keep = {k: v for k, v in out.items() if k.startswith('CHECK PLANT') or k.startswith('CHECK TT') or k == 'plants'}
    dump('p1_plants', keep)
    sys.exit(0)

# ------------------------------------------------------------------ the unplanted run
out['O1 symbolic residual is zero'] = {k: bool(v) for k, v in res.items()}
out['CHECK O1: every symbolic residual is zero (R1-R8)'] = all(bool(v) for v in res.values())
for cname, rows in (('Z1L0', Z1L0), ('Z2L1', Z2L1)):
    rec = audit(cname, rows)
    tot = rec.pop('_total')
    per = {}
    for name in SCOPE:
        l, A_l, flux, gg = SCOPE[name]
        vals = {g: tot(gg[g]) for g in ('G_a', 'G_b', 'G_c')}
        per[name] = {'A_l': fr(A_l), 'A_l/flux': fr(A_l/flux), 'totals': {g: fr(v) for g, v in vals.items()},
                     'rows': {g: {lab: fr(F(gg[g])) for lab, F in rows} for g in ('G_a', 'G_b', 'G_c')},
                     'equals A_l in all gauges': all(v == A_l for v in vals.values()),
                     'A_l/flux = 1/((2l+1)(2l+3))': A_l/flux == Fraction(1, (2*l + 1)*(2*l + 3))}
    rec['on_shell'] = per
    out[cname] = rec
    out[f'CHECK {cname}: audit passes (every row in span, k = 2, constraint weight N, identity on all probes)'] = rec['pass']
    out[f'CHECK {cname}: total = A_l in every gauge at l = 2, 3, 4, both m'] = all(p['equals A_l in all gauges'] for p in per.values())
    out[f'CHECK {cname}: A_l/flux = 1/((2l+1)(2l+3))'] = all(p['A_l/flux = 1/((2l+1)(2l+3))'] for p in per.values())

# O3: the mixed reading Z1L1
ins, a, b, bad = classify(R7_Z1)
o3 = {'coefficient dG_Z1/dR': str(dG1), 'R7_Z1 = a dA + b dV with (a, b)': [str(a), str(b)],
      'Z1L1 total (a, b) and effective k': [str(-2 + a), str(4 + b), str(-(4 + b)/(-2 + a))], 'in "8 pi G H - 2V" units': str(dG1/(-2)),
      'R7_Z1 is a combination of dA and dV': ins, 'n_probes_with_residual': len(bad),
      'value on pure-gauge probes (l=2 sectoral, l=2 zonal, generic xi)': [fr(R7_Z1(GAUGE_PROBE)), fr(R7_Z1(GAUGE_PROBE_Z)), fr(R7_Z1(data_of(GQ)))],
      'value on a dilation xi = x/10': fr(R7_Z1(data_of(symgrad([x/10, y/10, z/10])))),
      'value on the source-free TT probe': fr(R7_Z1(TT_DATA)), 'on_shell': {}}
for name in SCOPE:
    l, A_l, flux, gg = SCOPE[name]
    o3['on_shell'][name] = {'A_l': fr(A_l), **{g: fr(R7_Z1(gg[g])) for g in ('G_a', 'G_b', 'G_c')},
                           'extra/A_l in G_a': fr(R7_Z1(gg['G_a'])/A_l)}
gauge_dep = R7_Z1(data_of(GQ)) != 0 or R7_Z1(data_of(symgrad([x/10, y/10, z/10]))) != 0
o3['gauge dependent'] = gauge_dep
o3['modifier'] = 'Z1L1-DIRTY' if (gauge_dep or not ins) else 'Z1L1-CLEAN'
out['O3'] = o3
out['CHECK O3 computed (modifier assigned)'] = o3['modifier'] in ('Z1L1-DIRTY', 'Z1L1-CLEAN')
out['probe space'] = f'{len(SLOTS)} function slots x {len(MONO)} monomials (degree <= 3) = {len(BASIS)} probes'
p1 = (out['CHECK O1: every symbolic residual is zero (R1-R8)'] and out['CHECK Z1L0: audit passes (every row in span, k = 2, constraint weight N, identity on all probes)']
      and out['CHECK Z2L1: audit passes (every row in span, k = 2, constraint weight N, identity on all probes)']
      and out['CHECK Z1L0: total = A_l in every gauge at l = 2, 3, 4, both m'] and out['CHECK Z2L1: total = A_l in every gauge at l = 2, 3, 4, both m']
      and all(v for k, v in out.items() if k.startswith('CHECK PLANT')) and all(v for k, v in out.items() if k.startswith('CHECK gate')))
out['PART1_BIN'] = 'P1-PASS' if p1 else 'see failed checks'
dump('p1_enum', out)
print('PART1_BIN', out['PART1_BIN'], '| O3', o3['modifier'])
