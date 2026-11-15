"""Written AFTER reading ss-tension/REPORT.md, to locate the one numerical disagreement (their T4: "off-diagonal
correlator vanishes in every D"; mine: non-zero in D = 5, 6).
Their T4 keeps all vectors in a 4-dimensional subspace (point 2 displaced inside the plane that also carries the
polarisation) and only changes the power of the propagator. Test: eps = 0 kernels for
  A: displacement inside the polarisation plane  (their configuration, embedded in D dimensions)
  B: displacement along a transverse direction orthogonal to the polarisation plane (exists only for D >= 5).
Maxwell closed form at eps = 0 (derived by hand):  <phi_1 phi_2> = -2 p (p-1) (cL - ab) sigma^{-p-1} [x norm],  p = (D-2)/2,
c = m1.m2', L = l1.l2, a = m1.l2, b = l1.m2'.
"""
import sympy as sp
from common import dot, Ttensor, apply_T, power_derivs

R = sp.Rational
I = sp.I

def run(D, cfg, u1=R(7, 10), u2=R(3, 10), cth=R(3, 5), sth=R(4, 5)):
    def vec(**kw):
        v = [sp.Integer(0)]*D
        for k, val in kw.items():
            v[int(k[1:])] = val
        return v
    last = D - 1
    l1 = vec(e0=1, **{'e%d' % last: 1})
    m1 = vec(e1=1, e2=I)
    if cfg == 'A':
        l2 = vec(e0=1, e1=sth, **{'e%d' % last: cth})
        m2b = vec(e1=cth, e2=-I, **{'e%d' % last: -sth})
    else:
        l2 = vec(e0=1, e3=sth, **{'e%d' % last: cth})
        m2b = vec(e1=1, e2=-I)
    assert dot(m2b, l2) == 0 and dot(m2b, m2b) == 0 and dot(l2, l2) == 0
    xs = sp.symbols('x0:%d' % D)
    T = Ttensor(l1, m1, l2, m2b)
    f2 = apply_T(T, {0: sp.Integer(1)}, xs)
    f4 = apply_T(T, f2, xs)
    xv = [u1*l1[a] - u2*l2[a] for a in range(D)]
    sub = dict(zip(xs, xv))
    sigma = dot(xv, xv)
    p = R(D - 2, 2)
    g = power_derivs(p, sigma, 4)
    grav = sp.simplify(sum(P.subs(sub)*g[k] for k, P in f4.items()))
    maxw = sp.simplify(-sum(P.subs(sub)*g[k] for k, P in f2.items()))
    cLab = dot(m1, m2b)*dot(l1, l2) - dot(m1, l2)*dot(l1, m2b)
    pred = sp.simplify(-2*p*(p - 1)*cLab/sigma)
    return grav, maxw, pred, cLab

for D in (4, 5, 6):
    for cfg in ('A', 'B'):
        if cfg == 'B' and D == 4:
            continue
        grav, maxw, pred, cLab = run(D, cfg)
        print('D=%d cfg %s: cL-ab = %s | Maxwell kernel*sigma^p = %s (hand formula %s) | graviton kernel*sigma^p = %s' % (D, cfg, cLab, maxw, pred, grav))
