"""Route A': memory moments of a diamond observable on the cone, D = 4. Exact rational arithmetic.
Observable: O = R_abcd(p) A^a B^b A'^c B'^d at an interior point p of the diamond (any gamma[f] is a smearing of these).
Its commutator function on the cone is  Psi_0[E O](u l) = const * T^{ab} T'^{cd} d^4 Delta(u l - p),  Delta ~ delta(sigma),
with T from (l, m; A, B), T' from (l, m; A', B')  [mode sum: (1/2) e(v,v) against (1/2) e(w,w'), TT projector, v.v = 0].
sigma(u) = (u l - p)^2 = p^2 - 2 u l.p is LINEAR in u, so
  M_n := int u^n j_F du,  j_F = d_u^2 F = -2 u Psi_0[E O],  is a finite sum of derivatives of polynomials at u_x = p^2/(2 l.p).
Meaning: with F(u) = int_0^u (u-u') j_F the tip-regular shear data of E O,  F(u > u_x) = u M0 - M1  and  int_0^1 F = (M0 - 2 M1 + M2)/2.
Prediction to test (needed for (B.25)/(B.26) to be well defined on gauge-invariant observables):
  M1 = 0, M2 = 0 for every generator; M0(n) = c_ij m^i m^j with a direction-independent c (pure l = 2 "electric" = linear gauge).
"""
import sympy as sp
from common import direction, dot, Twedge, apply_T

R = sp.Rational
uu = sp.symbols('u')
xs = sp.symbols('x0:4')

def moments(par, p, A, B, Ap, Bp, ns=(-1, 0, 1, 2, 3)):
    l, m = direction(par)
    T = Twedge(l, m, A, B)
    Tp = Twedge(l, m, Ap, Bp)
    f4 = apply_T(Tp, apply_T(T, {0: sp.Integer(1)}, xs), xs)
    sub = {xs[a]: uu*l[a] - p[a] for a in range(4)}
    alpha = dot(p, p); b = -2*dot(l, p)
    ux = -alpha/b
    assert 0 < ux < 1 and b > 0
    out = {}
    for n in ns:
        tot = 0
        for k, P in f4.items():
            hfun = uu**(n + 1)*sp.expand(P.subs(sub))
            tot += (-1)**k*b**(-k - 1)*sp.diff(hfun, uu, k).subs(uu, ux)
        out[n] = sp.nsimplify(sp.simplify(tot))
    return out, m, ux

cases = [
    ((R(4, 5), R(1, 5), R(-1, 10), R(3, 10)), (R(1), R(2, 3), R(-1, 2), R(1, 4)), (R(0), R(1), R(1, 3), R(-2)),
     (R(1, 2), R(-1), R(2), R(1, 5)), (R(1, 3), R(1, 7), R(-3), R(1))),
    ((R(6, 5), R(-1, 4), R(1, 3), R(1, 10)), (R(0), R(1), R(0), R(0)), (R(0), R(0), R(1), R(0)),
     (R(0), R(1), R(0), R(0)), (R(0), R(0), R(1), R(0))),            # purely spatial component R_{1212}
    ((R(1, 2), R(1, 10), R(1, 5), R(-1, 10)), (R(1), R(0), R(0), R(0)), (R(0), R(1), R(0), R(0)),
     (R(1), R(0), R(0), R(0)), (R(0), R(0), R(0), R(1))),            # electric component R_{0103}
]
dirs = [(R(1, 3), R(-1, 5)), (R(2, 7), R(3, 4)), (R(-2, 3), R(1, 2)), (R(1, 9), R(-4, 5)), (R(-1, 2), R(-1, 3)),
        (R(3, 2), R(1, 7)), (R(0), R(1, 4)), (R(-5, 4), R(2, 3)), (R(1, 6), R(1, 13)), (R(2), R(-1)), (R(-1, 8), R(5, 3)), (R(3, 5), R(-7, 4))]
ok = True
for ic, (p, A, B, Ap, Bp) in enumerate(cases):
    print('== case %d: p = %s ==' % (ic, p))
    M0s = []
    for par in dirs:
        mom, m, ux = moments(par, p, A, B, Ap, Bp)
        M0s.append((m, mom[0]))
        print('  dir %-14s u_x = %-9s M1 = %s  M2 = %s | M0 = %s  M3 = %s  M_-1 = %s' % (
            str(par), sp.N(ux, 4), mom[1], mom[2], sp.N(mom[0], 6), sp.N(mom[3], 6), sp.N(mom[-1], 6)))
        ok &= (mom[1] == 0 and mom[2] == 0 and mom[-1] == 0 and mom[3] != 0)
    ok &= all(val == 0 for _, val in M0s)
    print('  M0 = 0 on all directions (stronger than the preregistered "M0 = c_ij m^i m^j": c = 0):', all(val == 0 for _, val in M0s))
# random rational cases (polynomial identity testing): interior points, generic vectors, two directions each
import random
random.seed(73)
def rq(lo=-3, hi=3):
    return R(random.randint(lo*6, hi*6), random.randint(2, 9))
nrand = 0
while nrand < 8:
    t = R(random.randint(3, 15), 10)
    xv = [R(random.randint(-9, 9), 20) for _ in range(3)]
    rr = sp.sqrt(sum(c**2 for c in xv))
    if not (rr < t and rr < 2 - t):
        continue
    p = (t, *xv)
    A, B, Ap, Bp = [tuple(rq() for _ in range(4)) for _ in range(4)]
    line = []
    for par in random.sample(dirs, 2):
        try:
            mom, m, ux = moments(par, p, A, B, Ap, Bp)
        except AssertionError:
            continue
        good = all(mom[n] == 0 for n in (-1, 0, 1, 2))
        ok &= good
        line.append('M_-1..M2 all zero: %s, M3 = %s' % (good, sp.N(mom[3], 5)))
    print('== random case %d: p = %s ->' % (nrand, p), ' | '.join(line))
    nrand += 1
print('teeth: M3 is non-zero above, so the vanishing of M_-1, M0, M1, M2 is not an artefact of the machinery.')
print('reading: int u^k Psi_0[E O] du = 0 for k = 0,1,2,3 on every generator, i.e. Psi_0[E O] = d_u^4 chi with chi compactly supported in (0,1).')
print('S4 PASS' if ok else 'S4 FAIL')
