"""C4: Sachs relation on the light cone, by brute-force linearised Riemann.
Minkowski in (u, v, th, ph): t = u+v, r = u-v, ds^2 = -4 du dv + (u-v)^2 dOmega^2. Cone C-: v = 0, u affine.
Perturbation in Gaussian null gauge: h_va = 0; h_uu = v^2 a, h_uA = v b_A (vanish on the cone as B.20 demands);
h_AB generic. Claim to test:  R_{u Ahat u Bhat}|_{v=0} = -(1/2) u^-2 d_u( u^2 d_u shat_AB ),  shat = frame components.
In D = 4 with S = s/u = u shat this is  -(1/(2u)) d_u^2 S_AB.
"""
import sympy as sp

u, v, th, ph = X = sp.symbols('u v theta phi', positive=True)
r = u - v
g = sp.zeros(4)
g[0, 1] = g[1, 0] = -2
g[2, 2] = r**2
g[3, 3] = r**2*sp.sin(th)**2
gi = g.inv()

def F(name):
    return sp.Function(name)(u, v, th, ph)

h = sp.zeros(4)
h[0, 0] = v**2*F('a')
h[0, 2] = h[2, 0] = v*F('b1')
h[0, 3] = h[3, 0] = v*F('b2')
h[2, 2] = F('Hqq'); h[2, 3] = h[3, 2] = F('Hqp'); h[3, 3] = F('Hpp')

Gam = [[[sum(gi[a, d]*(sp.diff(g[d, b], X[c]) + sp.diff(g[d, c], X[b]) - sp.diff(g[b, c], X[d])) for d in range(4))/2
         for c in range(4)] for b in range(4)] for a in range(4)]

def cov_h(d, b, c):  # nabla_d h_bc
    return sp.diff(h[b, c], X[d]) - sum(Gam[e][d][b]*h[e, c] + Gam[e][d][c]*h[b, e] for e in range(4))

covh = [[[cov_h(d, b, c) for c in range(4)] for b in range(4)] for d in range(4)]
# delta Gamma_{a bc} (first index lowered) = (1/2)(nabla_b h_ac + nabla_c h_ab - nabla_a h_bc)
dG = [[[(covh[b][a][c] + covh[c][a][b] - covh[a][b][c])/2 for c in range(4)] for b in range(4)] for a in range(4)]

def cov_dG(e, a, b, c):  # nabla_e dG_{a bc}
    return sp.diff(dG[a][b][c], X[e]) - sum(Gam[f][e][a]*dG[f][b][c] + Gam[f][e][b]*dG[a][f][c] + Gam[f][e][c]*dG[a][b][f] for f in range(4))

def dRiem(a, b, c, d):  # delta R_{abcd} = nabla_c dG_{a db} - nabla_d dG_{a cb}
    return cov_dG(c, a, d, b) - cov_dG(d, a, c, b)

comps = {('th', 'th'): (2, 2, 1), ('th', 'ph'): (2, 3, sp.sin(th)), ('ph', 'ph'): (3, 3, sp.sin(th)**2)}
ok = True
for key, (A, B, norm) in comps.items():
    R = sp.simplify(dRiem(0, A, 0, B).subs(v, 0).doit())
    Rhat = R/(u**2*norm)
    shat = (h[A, B]/(u**2*norm)).subs(v, 0)
    pred = -sp.Rational(1, 2)*sp.diff(u**2*sp.diff(shat, u), u)/u**2
    diff = sp.simplify((Rhat - pred).doit())
    print(key, 'R_{uAuB}(frame) - [-(1/2)u^-2 d_u(u^2 d_u shat)] =', diff)
    ok &= (diff == 0)
    # D = 4 form
    S = sp.Function('S')(u)
    chk = sp.simplify(-sp.Rational(1, 2)*sp.diff(u**2*sp.diff(S/u, u), u)/u**2 + sp.diff(S, u, 2)/(2*u))
    ok &= (chk == 0)
# flat-plane sign convention check: R_{uAuB} = -(1/2) d_u^2 gamma_AB for gamma_AB(u) on a null plane
print('D=4 identity  -(1/2)u^-2 d_u(u^2 d_u(S/u)) = -S\'\'/(2u):', chk == 0)
# D = 6: shat = S/u^2
chk6 = sp.simplify(-sp.Rational(1, 2)*sp.diff(u**2*sp.diff(S/u**2, u), u)/u**2)
print('D=6 form of Psi in terms of S:', sp.expand(chk6))
print('C4 PASS' if ok else 'C4 FAIL')
