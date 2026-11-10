"""S3 (PREREG section 3): the radial sector of the helicity-2 two-point function.

psi = x.E.x.  G(r, r', c; t) = x_i x_j y_k y_l W^{EE}_{ij,kl}(t, x - y), c = cos(angle(x, y)).
Legendre projections g_l = int_{-1}^{1} P_l(c) G dc against d_l = int P_l(c) (1/s) dc, s = t^2 - |x - y|^2
(the common constant -1/(4 pi^2) is dropped from both).
Registered: g_0 = g_1 = 0; g_l/d_l is one constant per l (l = 2, 3, 4), proportional to (l+2)!/(l-2)!.
Also exact: psi transforms under zeta as a conformal scalar of dimension 1.
"""
import itertools
import json
import sys

import sympy as sp
from mpmath import mp, mpf

sys.path.insert(0, __file__.rsplit('/', 1)[0])
from common import T, X1, X2, X3, XI, S, R3, build_W, delta_matrix_part

OUT = __file__.rsplit('/', 2)[0] + '/out/s3_out.json'
out = {}
mp.dps = 40

Wn = build_W()
r, rp, c, t = sp.symbols('r rp c t', real=True)
sn = sp.sqrt(1 - c**2)
x = (0, 0, r)
y = (rp * sn, 0, rp * c)
sub = {T: t, X1: x[0] - y[0], X2: x[1] - y[1], X3: x[2] - y[2]}
G = sum(x[i] * x[j] * y[k] * y[l] * Wn[0][0][i, j, k, l] for i, j, k, l in itertools.product(R3, repeat=4))
G = sp.simplify(sp.expand(G.subs(sub, simultaneous=True)) / (S.subs(sub, simultaneous=True))**5)
s_expr = sp.simplify(S.subs(sub, simultaneous=True))
out['G_closed_form'] = str(sp.factor(G))
out['G_depends_only_on_r_rp_c_t'] = True   # by construction after substitution; sn must have dropped out:
out['G_free_of_sqrt'] = not sp.simplify(G).has(sp.sqrt(1 - c**2))
Gf = sp.lambdify((r, rp, c, t), G, 'mpmath')
sf = sp.lambdify((r, rp, c, t), s_expr, 'mpmath')

pairs = [(mpf('0.3'), mpf('0.7')), (mpf('0.2'), mpf('0.9')), (mpf('0.5'), mpf('0.8')),
         (mpf('1.0'), mpf('2.5')), (mpf('0.4'), mpf('3.0')), (mpf('1.3'), mpf('1.9'))]
rows = []
for (a, b) in pairs:
    for tt in (mpf(0), mpf('0.1')):
        row = {'r': float(a), 'rp': float(b), 't': float(tt)}
        for l in range(0, 5):
            gl = mp.quad(lambda cc: mp.legendre(l, cc) * Gf(a, b, cc, tt), [-1, 1])
            dl = mp.quad(lambda cc: mp.legendre(l, cc) / sf(a, b, cc, tt), [-1, 1])
            row['g%d' % l] = mp.nstr(gl, 20)
            row['d%d' % l] = mp.nstr(dl, 20)
            row['ratio%d' % l] = mp.nstr(gl / dl, 20)
        rows.append(row)
        print(row, flush=True)
out['rows'] = rows
scale = max(abs(mpf(rw['g2'])) for rw in rows)
out['max_abs_g0_over_scale'] = float(max(abs(mpf(rw['g0'])) for rw in rows) / scale)
out['max_abs_g1_over_scale'] = float(max(abs(mpf(rw['g1'])) for rw in rows) / scale)
for l in (2, 3, 4):
    vals = [mpf(rw['ratio%d' % l]) for rw in rows]
    out['ratio%d_mean' % l] = mp.nstr(sum(vals) / len(vals), 20)
    out['ratio%d_spread' % l] = float((max(vals) - min(vals)) / abs(vals[0]))
r2 = mpf(out['ratio2_mean'])
out['ratios_over_ratio2'] = [float(mpf(out['ratio%d_mean' % l]) / r2) for l in (2, 3, 4)]
out['paper_A_prediction_(l+2)!/(l-2)!_over_24'] = [1.0, 5.0, 15.0]
out['pass'] = (out['max_abs_g0_over_scale'] < 1e-10 and out['max_abs_g1_over_scale'] < 1e-10
               and all(out['ratio%d_spread' % l] < 1e-10 for l in (2, 3, 4))
               and all(abs(a - b) < 1e-9 for a, b in zip(out['ratios_over_ratio2'], [1.0, 5.0, 15.0])))

# exact: psi = x.E.x transforms as a dimension-1 scalar under zeta, for any symmetric E, B
tq, a1, a2, a3 = sp.symbols('t a1 a2 a3', real=True)
xs = (tq, a1, a2, a3)
Ef = [[sp.Function('E%d%d' % (min(i, j), max(i, j)))(*xs) for j in R3] for i in R3]
Bf = [[sp.Function('B%d%d' % (min(i, j), max(i, j)))(*xs) for j in R3] for i in R3]
zeta = (sp.Rational(1, 2) * (1 - tq**2 - a1**2 - a2**2 - a3**2), -tq * a1, -tq * a2, -tq * a3)
dE0, _ = delta_matrix_part(zeta, xs, Ef, Bf)
dE = [[sum(zeta[m] * sp.diff(Ef[i][j], xs[m]) for m in range(4)) + dE0[i][j] for j in R3] for i in R3]
sp3 = (a1, a2, a3)
psi = sum(sp3[i] * sp3[j] * Ef[i][j] for i in R3 for j in R3)
lhs = sum(sp3[i] * sp3[j] * dE[i][j] for i in R3 for j in R3)
lam = sum(sp.diff(zeta[m], xs[m]) for m in range(4)) / 4
rhs = sum(zeta[m] * sp.diff(psi, xs[m]) for m in range(4)) + lam * psi
out['psi_is_dimension_1_scalar_under_zeta'] = sp.simplify(sp.expand(lhs - rhs)) == 0
json.dump(out, open(OUT, 'w'), indent=1)
print({k: v for k, v in out.items() if k != 'rows'})
