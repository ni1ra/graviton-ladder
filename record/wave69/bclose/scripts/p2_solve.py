"""PART 2a of w69-bclose (exact): the O(eps^2) vacuum solution on the scope family, and the canonical energy.

For each member builds g_eps = eta + eps h1 + eps^2 h2 (h1 = 2 t K) with polynomial h2(t, x), and checks that the
FULL nonlinear 4D Ricci tensor vanishes at orders eps and eps^2, all ten components, exactly.
Also computes W = int_B omega^t(eta; h1, L_zeta h1) from the covariant symplectic current (P-tensor of
wave66/p0test s1_transfer.py) and checks |W| = 2 A_l, A_l = int_B N |K|^2  (16 pi G = 1).

Writes out/p2_solve.json and out/p2_metrics.json (the metrics, as strings, for p2_geo.py).
"""
import json, os
import sympy as sp
from common import *

half = sp.Rational(1, 2)
e_ = [-1, 1, 1, 1]
out = {}
metrics = {}


def d4(f, m):
    return sp.diff(f, X4[m])


def gamma_low(H):
    """Gamma_{s m n} = (d_m H_sn + d_n H_sm - d_s H_mn)/2."""
    return [[[sp.expand(half*(d4(H[s][n], m) + d4(H[s][m], n) - d4(H[m][n], s))) for n in range(4)] for m in range(4)]
            for s in range(4)]


def ricci_orders(h1, h2):
    """Exact eps^1 and eps^2 coefficients of the Ricci tensor of eta + eps h1 + eps^2 h2."""
    L1, L2 = gamma_low(h1), gamma_low(h2)
    G1 = [[[e_[r]*L1[r][m][n] for n in range(4)] for m in range(4)] for r in range(4)]
    G2 = [[[sp.expand(e_[r]*L2[r][m][n] - sum(e_[r]*e_[s]*h1[r][s]*L1[s][m][n] for s in range(4))) for n in range(4)]
           for m in range(4)] for r in range(4)]

    def lin(G):
        return [[sp.expand(sum(d4(G[r][m][n], r) - d4(G[r][r][m], n) for r in range(4))) for n in range(4)] for m in range(4)]
    R1, R2 = lin(G1), lin(G2)
    for m in range(4):
        for n in range(4):
            R2[m][n] = sp.expand(R2[m][n] + sum(G1[r][r][s]*G1[s][m][n] - G1[r][n][s]*G1[s][r][m]
                                               for r in range(4) for s in range(4)))
    return R1, R2


def ric3lin(c):
    tr = c[0][0] + c[1][1] + c[2][2]
    return [[sp.expand(half*(sum(sp.diff(c[k][j], X[k], X[i]) + sp.diff(c[k][i], X[k], X[j]) for k in range(3))
                             - lap(c[i][j]) - sp.diff(tr, X[i], X[j]))) for j in range(3)] for i in range(3)]


def embed(sp3):
    return [[sp.Integer(0) if (m == 0 or n == 0) else sp3[m - 1][n - 1] for n in range(4)] for m in range(4)]


def tcoef(expr, n):
    return sp.Poly(expr, t).coeff_monomial(t**n) if expr != 0 else sp.Integer(0)


def solve_h2(K, q2, c1=None):
    """Synchronous-gauge h2_ij = sum_n t^n c_n(x): c_0 = q2, c_1 given (default 0),
    (n+2)(n+1) c_{n+2} = -2 (Quad_n + Ric3lin[c_n]),  Quad = eps^2 Ricci of eta + eps h1 alone."""
    h1 = embed([[2*t*K[i][j] for j in range(3)] for i in range(3)])
    zero = [[sp.Integer(0)]*4 for _ in range(4)]
    _, Quad = ricci_orders(h1, zero)
    cs = {0: q2, 1: c1 if c1 is not None else [[sp.Integer(0)]*3 for _ in range(3)]}
    n = 0
    while True:
        r3 = ric3lin(cs[n])
        nxt = [[sp.expand(-2*(tcoef(Quad[i + 1][j + 1], n) + r3[i][j])/((n + 2)*(n + 1))) for j in range(3)] for i in range(3)]
        cs.setdefault(n + 2, nxt) if n + 2 not in cs else None
        cs[n + 2] = nxt
        n += 1
        if n > 40:
            raise RuntimeError('series did not terminate')
        if all(cs[n][i][j] == 0 for i in range(3) for j in range(3)) and all(cs[n + 1][i][j] == 0 for i in range(3) for j in range(3)) \
                and all(sp.degree(Quad[i + 1][j + 1], t) < n if Quad[i + 1][j + 1] != 0 else True for i in range(3) for j in range(3)):
            break
    h2 = embed([[sp.expand(sum(t**k*c[i][j] for k, c in cs.items())) for j in range(3)] for i in range(3)])
    return h1, h2


def check(h1, h2):
    R1, R2 = ricci_orders(h1, h2)
    return (all(R1[m][n] == 0 for m in range(4) for n in range(4)), all(R2[m][n] == 0 for m in range(4) for n in range(4)))


def lie_gauge(Xi):
    """(L_Xi eta)_{mn} = d_m Xi_n + d_n Xi_m with Xi_m = eta_mm Xi^m."""
    lo = [e_[m]*Xi[m] for m in range(4)]
    return [[sp.expand(d4(lo[n], m) + d4(lo[m], n)) for n in range(4)] for m in range(4)]


def P(a, b, c, d, e, f):
    dl = lambda i, j: e_[i] if i == j else 0
    return (dl(a, e)*dl(f, b)*dl(c, d) - half*dl(a, d)*dl(b, e)*dl(f, c) - half*dl(a, b)*dl(c, d)*dl(e, f)
            - half*dl(b, c)*dl(a, e)*dl(f, d) + half*dl(b, c)*dl(a, d)*dl(e, f))


zeta = [N - t**2/2, -t*x, -t*y, -t*z]


def canonical_W(h1):
    Lz = [[sp.expand(sum(zeta[r]*d4(h1[m][n], r) + h1[r][n]*d4(zeta[r], m) + h1[m][r]*d4(zeta[r], n) for r in range(4)))
           for n in range(4)] for m in range(4)]
    tot = 0
    for b in range(4):
        for c in range(4):
            for d in range(4):
                for e in range(4):
                    for f in range(4):
                        p = P(0, b, c, d, e, f)
                        if p != 0:
                            tot += p*(Lz[b][c]*d4(h1[e][f], d) - h1[b][c]*d4(Lz[e][f], d))
    return ball(sp.expand(tot).subs(t, 0))


def store(name, h1, h2, A_l, flux, note):
    metrics[name] = {'h1': [[str(h1[m][n]) for n in range(4)] for m in range(4)],
                     'h2': [[str(h2[m][n]) for n in range(4)] for m in range(4)],
                     'A_over_pi': str(A_l), 'flux_over_pi': str(flux), 'note': note}


phi_mix = sp.expand(sectoral(2) + zonal(2)/2 + 2*x*z + y*z/3 - x*y)
MEM = [('M1 l=2 m=2', 2, sectoral(2)), ('M2 l=3 m=3', 3, sectoral(3)), ('M3 l=4 m=4', 4, sectoral(4)),
       ('M4 l=2 m=0', 2, zonal(2)), ('M5 l=2 mixed m', 2, phi_mix)]
base = {}
for name, l, phi in MEM:
    K = magnetic_K(phi)
    k2 = K2(K)
    A_l, flux = ball(N*k2), sphere(k2)
    u = inv_lap(-k2/2)
    q2 = [[u if i == j else sp.Integer(0) for j in range(3)] for i in range(3)]
    h1, h2 = solve_h2(K, q2)
    ok1, ok2 = check(h1, h2)
    W = canonical_W(h1)
    rec = {'A_l/pi': str(A_l), 'flux/pi': str(flux), 'A_l/flux': str(A_l/flux),
           'h2 max power of t': max(sp.degree(h2[i][j], t) for i in range(1, 4) for j in range(1, 4) if h2[i][j] != 0),
           'h2 total degree': max(sp.Poly(h2[i][j], *X4).total_degree() for i in range(1, 4) for j in range(1, 4) if h2[i][j] != 0),
           'W/pi = int omega^t(h1, L_zeta h1)': str(W), 'E_can/pi = |W|/2': str(abs(W)/2)}
    out[name] = rec
    out[f'CHECK {name}: Ricci vanishes at O(eps), all 10 components'] = ok1
    out[f'CHECK {name}: Ricci vanishes at O(eps^2), all 10 components (full nonlinear)'] = ok2
    out[f'CHECK {name}: E_can = |W|/2 = A_l exactly'] = (abs(W)/2 == A_l)
    out[f'CHECK {name}: A_l/flux = 1/((2l+1)(2l+3))'] = (A_l/flux == sp.Rational(1, (2*l + 1)*(2*l + 3)))
    store(name, h1, h2, A_l, flux, 'vacuum through O(eps^2), synchronous gauge, time-symmetric')
    base[name] = (K, q2, h1, h2, A_l, flux, u)

# M6: M1 with (a) a t-odd homogeneous second-order term 2 t K', K' the l = 3 zonal magnetic tensor,
#     (b) a second-order gauge term L_Xi eta with Xi not vanishing on S (radial, tangential and time parts).
K, q2, h1, h2, A_l, flux, u = base['M1 l=2 m=2']
Kp = magnetic_K(zonal(3))
c1 = [[sp.expand(2*Kp[i][j]/5) for j in range(3)] for i in range(3)]
h1b, h2b = solve_h2(K, q2, c1)
Xi = [x*y/4 + t*z/5 + sp.Rational(1, 10)*z**2, x/8 + y*z/6 + t*x/7, -y/9 + x**2/7 + t**2/6, z/10 + x*y/5 - t*y/8]
LG = lie_gauge(Xi)
h2c = [[sp.expand(h2b[m][n] + LG[m][n]) for n in range(4)] for m in range(4)]
ok1, ok2 = check(h1b, h2c)
out['M6 M1 dressed'] = {'A_l/pi': str(A_l), 'Xi (t,x,y,z components)': [str(v) for v in Xi],
                        'radial part of Xi on S is nonzero': sphere(sum(X[i]*Xi[i + 1] for i in range(3)).subs(t, 0)**2) != 0,
                        'h2_tt, h2_ti nonzero': any(h2c[0][n] != 0 for n in range(4)),
                        'd_t h2_ij at t=0 nonzero (not time-symmetric)': any(sp.diff(h2c[i][j], t).subs(t, 0) != 0 for i in range(1, 4) for j in range(1, 4))}
out['CHECK M6 M1 dressed: Ricci vanishes at O(eps^2), all 10 components (full nonlinear)'] = ok1 and ok2
store('M6 M1 dressed', h1b, h2c, A_l, flux, 'M1 plus t-odd homogeneous term and second-order gauge term; corner moves at O(eps^2)')

# PLANT-SRC: constraint source scaled by 11/10 (NOT a solution)
q2p = [[sp.Rational(11, 10)*u if i == j else sp.Integer(0) for j in range(3)] for i in range(3)]
h2p = embed([[sp.expand(h2[i + 1][j + 1] + (q2p[i][j] - q2[i][j])) for j in range(3)] for i in range(3)])
ok1, ok2 = check(h1, h2p)
out['PLANT-SRC'] = {'is a solution at O(eps^2)': ok2}
out['CHECK PLANT-SRC is not a solution (exact Ricci check rejects it)'] = (not ok2)
store('PLANT-SRC (M1, source x 11/10)', h1, h2p, A_l, flux, 'plant: violates the O(eps^2) Hamiltonian constraint by 10%')

# PLANT-FLAT: g = psi^* eta, psi: t -> t + eps f(x).  Exactly flat for every eps; the true slice is t = -eps f.
f = sp.expand((x**2 - y**2)*z/2 + x*y/3 + z**2/4 - r2/12)
df = [sp.diff(f, s) for s in X]
h1f = [[sp.Integer(0)]*4 for _ in range(4)]
h2f = [[sp.Integer(0)]*4 for _ in range(4)]
for i in range(3):
    h1f[0][i + 1] = h1f[i + 1][0] = -df[i]
    for j in range(3):
        h2f[i + 1][j + 1] = sp.expand(-df[i]*df[j])
ok1, ok2 = check(h1f, h2f)
out['PLANT-FLAT'] = {'f': str(f), 'exact pullback': 'g = -(dt + eps df)^2 + dx^2'}
out['CHECK PLANT-FLAT: Ricci vanishes at O(eps) and O(eps^2)'] = ok1 and ok2
store('PLANT-FLAT', h1f, h2f, sp.Integer(0), sp.Integer(1), 'exactly flat; true corner at t = -eps f, r = 1')

with open(os.path.join(OUT, 'p2_metrics.json'), 'w') as fh:
    json.dump(metrics, fh)
dump('p2_solve', out)
