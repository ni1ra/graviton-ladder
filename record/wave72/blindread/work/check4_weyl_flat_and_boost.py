"""Blind check 4.
(a) Is the Weyl metric with U = eps*z (l=1) flat?  Compute its Kretschmann scalar exactly.
(b) For uniform (l=2) fields: derivative of the paper's 'area term' and 'S' under a small boost of the ball,
    using the paper's own flow rule  dE = 2 sym(beta x B), dB = -2 sym(beta x E)  (eq. 4 with x -> beta).
"""
import sympy as sp

# ---------- (a)
t, rho, z, ph, eps = sp.symbols('t rho z phi epsilon', real=True)
co = [t, rho, z, ph]
U = eps * z
k = -eps**2 * rho**2 / 2          # k_rho = rho(U_rho^2-U_z^2) = -eps^2 rho ; k_z = 0
g = sp.diag(-sp.exp(2 * U), sp.exp(-2 * U + 2 * k), sp.exp(-2 * U + 2 * k), sp.exp(-2 * U) * rho**2)
ginv = g.inv()
n = 4
Gam = [[[sp.simplify(sum(ginv[a, d] * (sp.diff(g[d, b], co[c]) + sp.diff(g[d, c], co[b]) - sp.diff(g[b, c], co[d])) for d in range(n)) / 2)
         for c in range(n)] for b in range(n)] for a in range(n)]
def Riem(a, b, c, d):  # R^a_{bcd}
    e = sp.diff(Gam[a][b][d], co[c]) - sp.diff(Gam[a][b][c], co[d])
    e += sum(Gam[a][c][m] * Gam[m][b][d] - Gam[a][d][m] * Gam[m][b][c] for m in range(n))
    return sp.simplify(e)
Rup = [[[[Riem(a, b, c, d) for d in range(n)] for c in range(n)] for b in range(n)] for a in range(n)]
# Ricci
Ric = sp.Matrix(n, n, lambda b, d: sp.simplify(sum(Rup[a][b][a][d] for a in range(n))))
print('(a) Weyl l=1 (U = eps z, k = -eps^2 rho^2/2): Ricci tensor zero:', Ric == sp.zeros(n, n))
# all-lower Riemann and Kretschmann
Rlow = [[[[sp.simplify(sum(g[a, m] * Rup[m][b][c][d] for m in range(n))) for d in range(n)] for c in range(n)] for b in range(n)] for a in range(n)]
Kr = 0
for a in range(n):
    for b in range(n):
        for c in range(n):
            for d in range(n):
                if Rlow[a][b][c][d] != 0:
                    Kr += Rlow[a][b][c][d]**2 * ginv[a, a] * ginv[b, b] * ginv[c, c] * ginv[d, d]
Kr = sp.simplify(Kr)
print('    Kretschmann scalar =', Kr)
print('    leading order in eps:', sp.series(Kr, eps, 0, 6))
print('    R_{t rho t rho} =', Rlow[0][1][0][1], ' ; series:', sp.series(Rlow[0][1][0][1], eps, 0, 4))

# ---------- (b)
LC = sp.LeviCivita
def symtl(name):
    s = sp.symbols(name + '11 ' + name + '12 ' + name + '13 ' + name + '22 ' + name + '23', real=True)
    return sp.Matrix([[s[0], s[1], s[2]], [s[1], s[3], s[4]], [s[2], s[4], -s[0] - s[3]]])
E = symtl('E'); B = symtl('B')
b = sp.Matrix(sp.symbols('b1 b2 b3', real=True))
def cross(a, T):
    return sp.Matrix(3, 3, lambda i, j: sum(LC(i, k, l) * a[k] * T[l, j] for k in range(3) for l in range(3)))
def sym(M):
    return (M + M.T) / 2
def sq(M):
    return sum(M[i, j]**2 for i in range(3) for j in range(3))
dE = 2 * sym(cross(b, B)); dB = -2 * sym(cross(b, E))
lam = sp.Symbol('lam')
area = sq(E + lam * dE) + sp.Rational(4, 3) * sq(B + lam * dB)
ent = sq(E + lam * dE) + sq(B + lam * dB)
kret = sq(E + lam * dE) - sq(B + lam * dB)
P = sp.Matrix([sum(LC(k, i, l) * (E * B.T)[i, l] for i in range(3) for l in range(3)) for k in range(3)])  # P_k = eps_kil E_ij B_lj
d_area = sp.expand(sp.diff(area, lam).subs(lam, 0))
d_ent = sp.expand(sp.diff(ent, lam).subs(lam, 0))
d_kret = sp.expand(sp.diff(kret, lam).subs(lam, 0))
bP = sp.expand((b.T * P)[0])
print('(b) super-Poynting vector P_k = eps_kil E_ij B_lj')
print('    d/dbeta [E^2 - B^2]        =', sp.simplify(d_kret), ' (invariant, as it must be)')
print('    d/dbeta [E^2 + B^2]        = %s * (beta.P)' % sp.simplify(d_ent / bP))
print('    d/dbeta [E^2 + (4/3) B^2]  = %s * (beta.P)' % sp.simplify(d_area / bP))
comm = E * B - B * E
print('    P = 0 iff [E,B]=0: P components', [sp.factor(p) for p in P], ' ; [E,B] entries (01,02,12):', sp.factor(comm[0, 1]), sp.factor(comm[0, 2]), sp.factor(comm[1, 2]))
