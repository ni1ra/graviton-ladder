"""Shared helpers of lane SS-SECOND (own code, nothing imported from other lanes).
Exact arithmetic for kernels of the form  T^{ab} T'^{cd} d_a d_b d_c d_d g(sigma),  sigma = eta_ab x^a x^b,
represented as sum_k P_k(x) g^{(k)}(sigma) with d_a g^{(k)} = 2 x_a g^{(k+1)}.
"""
import sympy as sp

def eta(D):
    return [-1] + [1]*(D - 1)

def dot(a, b):
    D = len(a)
    e = eta(D)
    return sum(e[i]*a[i]*b[i] for i in range(D))

def direction(params):
    """Rational point of S^{D-2} from stereographic parameters; returns l = (1, n) and complex null m (unnormalised),
    m = d n/d a_1 + i d n/d a_2, with m.l = 0, m.m = 0."""
    k = len(params)
    a = sp.symbols('a0:%d' % k)
    s = sum(ai**2 for ai in a)
    n = [2*ai/(1 + s) for ai in a] + [(1 - s)/(1 + s)]
    sub = dict(zip(a, params))
    d1 = [sp.diff(c, a[0]) for c in n]
    d2 = [sp.diff(c, a[1]) for c in n]
    nn = [sp.nsimplify(c.subs(sub)) for c in n]
    m = [sp.nsimplify((p + sp.I*q).subs(sub)) for p, q in zip(d1, d2)]
    l = [sp.Integer(1)] + nn
    m = [sp.Integer(0)] + m
    return l, m

def conj(vv):
    return [sp.conjugate(c) for c in vv]

def Ttensor(l1, m1, l2, m2):
    """T^{ab} k_a k_b = v1.v2 with v_i = (l_i.k) m_i - (m_i.k) l_i."""
    D = len(l1)
    c, a_, b_, L = dot(m1, m2), dot(m1, l2), dot(l1, m2), dot(l1, l2)
    return [[l1[a]*l2[b]*c - l1[a]*m2[b]*a_ - m1[a]*l2[b]*b_ + m1[a]*m2[b]*L for b in range(D)] for a in range(D)]

def Twedge(l, m, A, B):
    """T^{ab} k_a k_b = v.w with v = (l.k) m - (m.k) l and w = (A.k) B - (B.k) A (general Riemann component at the other point)."""
    D = len(l)
    return [[l[a]*A[b]*dot(m, B) - l[a]*B[b]*dot(m, A) - m[a]*A[b]*dot(l, B) + m[a]*B[b]*dot(l, A) for b in range(D)] for a in range(D)]

def apply_T(T, f, xs):
    """f: dict k -> polynomial in xs (coefficient of g^{(k)}). Returns dict for T^{ab} d_a d_b f."""
    D = len(xs)
    e = eta(D)
    def d(fd, a):
        out = {}
        for k, P in fd.items():
            dP = sp.diff(P, xs[a])
            if dP != 0:
                out[k] = out.get(k, 0) + dP
            out[k + 1] = out.get(k + 1, 0) + 2*e[a]*xs[a]*P
        return out
    res = {}
    for b in range(D):
        fb = d(f, b)
        for a in range(D):
            if T[a][b] == 0:
                continue
            fab = d(fb, a)
            for k, P in fab.items():
                res[k] = res.get(k, 0) + T[a][b]*P
    return {k: sp.expand(P) for k, P in res.items()}

def power_derivs(p, sigma, kmax):
    """g = sigma^-p: returns [g^{(k)} * sigma^p] for k = 0..kmax (so results stay rational when p is half-integer)."""
    out = []
    c = sp.Integer(1)
    for k in range(kmax + 1):
        out.append(c*sigma**(-k))
        c = c*(-p - k)
    return out
