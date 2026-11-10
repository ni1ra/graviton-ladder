"""Shared exact helpers: polynomial integrals over the ball B_R and sphere S_R in R^3."""
import sympy as sp

x, y, z = X = sp.symbols('x y z', real=True)

def _mono_sphere(i, j, k):
    # \oint_{S^2} x^i y^j z^k dA (unit sphere); zero unless all even
    if i % 2 or j % 2 or k % 2:
        return sp.Integer(0)
    a, b, c = sp.Rational(i+1, 2), sp.Rational(j+1, 2), sp.Rational(k+1, 2)
    return 2*sp.gamma(a)*sp.gamma(b)*sp.gamma(c)/sp.gamma(a+b+c)

def sphere_int(expr, R=1):
    r"""\oint_{S_R} expr dA for polynomial expr."""
    P = sp.Poly(sp.expand(expr), *X)
    tot = 0
    for (i, j, k), c in P.terms():
        tot += c*_mono_sphere(i, j, k)*sp.Integer(R)**(i+j+k+2) if isinstance(R, int) else c*_mono_sphere(i, j, k)*R**(i+j+k+2)
    return sp.simplify(tot)

def ball_int(expr, R=1):
    r"""\int_{B_R} expr dV for polynomial expr."""
    P = sp.Poly(sp.expand(expr), *X)
    tot = 0
    for (i, j, k), c in P.terms():
        n = i+j+k
        tot += c*_mono_sphere(i, j, k)*sp.Integer(1)*R**(n+3)/(n+3)
    return sp.simplify(tot)

def d(f, i):
    return sp.diff(f, X[i])

def R_lin(q):
    """Linearised 3D scalar curvature: d_i d_j q_ij - lap q_kk."""
    tr = sum(q[k][k] for k in range(3))
    s = sum(d(d(q[i][j], i), j) for i in range(3) for j in range(3))
    return sp.expand(s - sum(d(d(tr, k), k) for k in range(3)))

def dA(q, R=1):
    r"""delta A = 1/2 \oint_{S_R} (q_kk - n_i n_j q_ij) dA, n = x/R."""
    tr = sum(q[k][k] for k in range(3))
    qnn = sum(X[i]*X[j]*q[i][j] for i in range(3) for j in range(3))/R**2
    return sp.Rational(1, 2)*sphere_int(tr - qnn, R)

def dV(q, R=1):
    tr = sum(q[k][k] for k in range(3))
    return sp.Rational(1, 2)*ball_int(tr, R)

def symgrad(xi):
    return [[sp.expand(d(xi[j], i) + d(xi[i], j)) for j in range(3)] for i in range(3)]

def rand_poly(rng, deg, scale=5):
    terms = []
    for i in range(deg+1):
        for j in range(deg+1-i):
            for k in range(deg+1-i-j):
                terms.append(sp.Rational(rng.randint(-scale, scale), rng.randint(1, scale))*x**i*y**j*z**k)
    return sp.Add(*terms)
