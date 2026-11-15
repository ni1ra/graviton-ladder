"""SS-LEMMA script A: the identities the hand proof rests on, exact symbolic (sympy, no nsimplify).
Signature (-+++). w = -l.y, mu = m.y, s = -y.y, v = mu l + w m.
"""
import json, os, sys
import sympy as sp

R = sp.Rational
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'out', 'lemma_proof.json')
res = {}


def eta(D):
    return [-1] + [1]*(D - 1)


def dot(a, b):
    e = eta(len(a))
    return sum(e[i]*a[i]*b[i] for i in range(len(a)))


def lower(a):
    e = eta(len(a))
    return [e[i]*a[i] for i in range(len(a))]


def zero(expr):
    return sp.simplify(sp.cancel(sp.together(sp.expand(expr)))) == 0


# ---------------- A1: epsilon algebra, symbolic k ----------------
k, eps, mu = sp.symbols('k epsilon mu')
w = sp.symbols('w', positive=True)
fmm, flm, fll = sp.symbols('f_mm f_lm f_ll')
weps = w - eps*mu
f_lm_eps = flm + eps*fmm                 # f(l + eps m, m)
f_ll_eps = fll + 2*eps*flm + eps**2*fmm  # f(l + eps m, l + eps m)
lhs = ((k - 1)*(k - 2)*fmm*w**(1 - k)
       + 2*(k - 2)*sp.diff(f_lm_eps*weps**(1 - k), eps).subs(eps, 0)
       + sp.diff(f_ll_eps*weps**(1 - k), eps, 2).subs(eps, 0))
rhs = k*(k - 1)*w**(-k - 1)*(fmm*w**2 + 2*flm*mu*w + fll*mu**2)
a1 = sp.simplify(sp.expand((lhs - rhs)*w**(k + 1))) == 0
res['A1'] = bool(a1)
print('A1 epsilon algebra, symbolic k: residual zero =', a1)


# ---------------- frames ----------------
def frame(n, e1, e2):
    """l = (1, n), l' = (1, -n), m = (0, e1 + i e2), mb' = 2 conj(m)/(m.conj m)."""
    l = [sp.Integer(1)] + list(n)
    lp = [sp.Integer(1)] + [-c for c in n]
    m = [sp.Integer(0)] + [a + sp.I*b for a, b in zip(e1, e2)]
    mb = [sp.conjugate(c) for c in m]
    N = dot(m, mb)
    mbp = [2*c/N for c in mb]
    X = sp.expand
    assert X(dot(l, l)) == 0 and X(dot(m, m)) == 0 and X(dot(l, m)) == 0 and X(dot(lp, m)) == 0 and X(dot(l, lp)) == -2
    assert sp.simplify(dot(m, mbp) - 2) == 0 and X(dot(l, mbp)) == 0 and X(dot(lp, mbp)) == 0
    return l, lp, m, mbp


canon = frame((0, 0, 1), (1, 0, 0), (0, 1, 0))
# generic rational direction: n = (2/3, -1/3, 2/3); tangent vectors (not orthonormal: only m.m = 0 is needed)
# e1, e2 orthogonal to n, equal length, mutually orthogonal
ngen = (R(2, 3), R(-1, 3), R(2, 3))
e1 = (R(1, 3), R(-2, 3), R(-2, 3))
e2 = (R(2, 3), R(2, 3), R(-1, 3))
generic = frame(ngen, e1, e2)

ys = sp.symbols('y0:4')


def fields(fr, c3=R(1, 3)):
    l, lp, m, mbp = fr
    W = -dot(l, ys); P = -dot(lp, ys); Z = dot(m, ys); Zb = dot(mbp, ys)
    S = -dot(ys, ys)
    assert sp.expand(S - (W*P - Z*Zb)) == 0
    ll, lpl, ml, mbl = lower(l), lower(lp), lower(m), lower(mbp)
    v = [Z*ll[a] + W*ml[a] for a in range(4)]
    xi2 = [Z**2/(2*W**2)*ll[a] + Z/W*ml[a] for a in range(4)]
    xi3 = [-(-P*Z**2/(2*W**2) + c3*Z**3*Zb/W**3)*ll[a] + Z**2/(2*W)*lpl[a]
           + (P*Z/W - Z**2*Zb/(2*W**2))*ml[a] + Z**3/(6*W**2)*mbl[a] for a in range(4)]
    return W, S, v, xi2, xi3


def symgrad_residuals(xi, t):
    bad = 0
    for a in range(4):
        for b in range(a, 4):
            r = (sp.diff(xi[b], ys[a]) + sp.diff(xi[a], ys[b]))/2 - t(a, b)
            if not zero(r):
                bad += 1
    return bad


# ---------------- A2 and A5 ----------------
a2_bad = 0
for name, fr in (('canonical', canon), ('generic', generic)):
    W, S, v, xi2, xi3 = fields(fr)
    b2 = symgrad_residuals(xi2, lambda a, b: v[a]*v[b]/W**3)
    b3 = symgrad_residuals(xi3, lambda a, b: S*v[a]*v[b]/W**4)
    print('A2 %-9s: non-zero residual components  xi2: %d/10   xi3: %d/10' % (name, b2, b3))
    a2_bad += b2 + b3
res['A2'] = (a2_bad == 0)
W, S, v, xi2, xi3p = fields(canon, c3=R(1, 2))
plant_bad = symgrad_residuals(xi3p, lambda a, b: S*v[a]*v[b]/W**4)
res['A5_plant_caught'] = plant_bad > 0
print('A5 PLANT (1/3 -> 1/2 in xi3): non-zero residual components = %d  -> plant caught: %s' % (plant_bad, plant_bad > 0))


# ---------------- A3, A4: linearised Riemann ----------------
def riemann_components(t, yv):
    D = len(yv)
    comps = {}
    pairs = [(a, b) for a in range(D) for b in range(a + 1, D)]
    for i, (a, b) in enumerate(pairs):
        for (c, d) in pairs[i:]:
            comps[(a, b, c, d)] = (sp.diff(t[b][c], yv[a], yv[d]) + sp.diff(t[a][d], yv[b], yv[c])
                                   - sp.diff(t[b][d], yv[a], yv[c]) - sp.diff(t[a][c], yv[b], yv[d]))
    return comps


W, S, v, _, _ = fields(canon)
a3 = {}
for kk in (2, 3, 4):
    t = [[S**(kk - 2)*W**(-kk - 1)*v[a]*v[b] for b in range(4)] for a in range(4)]
    comps = riemann_components(t, ys)
    nz = sum(0 if zero(c) else 1 for c in comps.values())
    a3[kk] = nz
    print('A3 D = 4, k = %d: non-zero linearised Riemann components = %d of %d' % (kk, nz, len(comps)))
res['A3'] = (a3[2] == 0 and a3[3] == 0 and a3[4] > 0)
res['A3_counts'] = {str(a): b for a, b in a3.items()}

y5 = sp.symbols('y0:5')
l5 = [1, 0, 0, 1, 0]; m5 = [0, 1, sp.I, 0, 0]
W5 = -dot(l5, y5); Z5 = dot(m5, y5); S5 = -dot(y5, y5)
v5 = [Z5*lower(l5)[a] + W5*lower(m5)[a] for a in range(5)]
pt = dict(zip(y5, (R(1), R(1, 5), R(-1, 7), R(1, 4), R(1, 3))))
a4 = {}
for kk in (2, 3):
    t = [[S5**(kk - R(5, 2))*W5**(-kk - 1)*v5[a]*v5[b] for b in range(5)] for a in range(5)]
    comps = riemann_components(t, y5)
    nz = sum(0 if sp.simplify(c.subs(pt)) == 0 else 1 for c in comps.values())
    a4[kk] = nz
    print('A4 D = 5 control, k = %d: non-zero components at the rational point = %d of %d' % (kk, nz, len(comps)))
res['A4'] = (a4[2] > 0 and a4[3] > 0)

ok = res['A1'] and res['A2'] and res['A3'] and res['A4'] and res['A5_plant_caught']
res['PASS'] = bool(ok)
with open(OUT, 'w', encoding='utf-8') as fh:
    json.dump(res, fh, indent=1)
print('SCRIPT A', 'PASS' if ok else 'FAIL')
sys.exit(0 if ok else 1)
