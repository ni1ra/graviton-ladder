"""Independent exact probe (MODFLOW-REF): Step 6 geometry and the slice form of delta_zeta.
Signature (+,-,-,-). kappa(x) = (x - b x^2)/sigma, b = e1 (spatial), sigma = 1 - 2 b.x + b^2 x^2.
Checks: (1) D kappa . zeta_1 = X_W o kappa, X_W = -(x1-1/2) d_t - t d_1;
(2) (kappa x - kappa y)^2 = (x-y)^2/(sigma(x) sigma(y));
(3) zeta_1 at t=0: lambda=0, A=0, beta=-x, zeta^0=N;
(4) kappa maps sample points of O into the wedge x1 < 1/2 - |t|, and unit sphere to the plane x1=1/2.
"""
import sympy as sp, random
t,x,y,z = sp.symbols('t x y z', real=True)
X = sp.Matrix([t,x,y,z])
def msq(v): return v[0]**2 - v[1]**2 - v[2]**2 - v[3]**2
def kappa(v):
    x2 = msq(v); bx = -v[1]            # b.x with b=(0,1,0,0): -x1
    sig = 1 - 2*bx + (-1)*x2
    b = sp.Matrix([0,1,0,0])
    return (v - b*x2)/sig, sig
K, sig = kappa(X)
zeta = sp.Matrix([(1 - t**2 - x**2 - y**2 - z**2)/2, -t*x, -t*y, -t*z])
J = K.jacobian(X)
lhs = (J*zeta).applyfunc(sp.simplify)
XW = sp.Matrix([-(K[1]-sp.Rational(1,2)), -K[0], 0, 0])
print('1 Dkappa.zeta == X_W o kappa:', (lhs - XW).applyfunc(sp.simplify) == sp.zeros(4,1))
a = sp.symbols('a0:4', real=True); Y = sp.Matrix(a)
KY, sigY = kappa(Y)
print('2 interval identity:', sp.simplify(msq(K-KY) - msq(X-Y)/(sig*sigY)) == 0)
lam = sum(sp.diff(zeta[m], X[m]) for m in range(4))/4
A = sp.Matrix(3,3, lambda i,m: sp.diff(zeta[1+m], X[1+i]) - (lam if i==m else 0))
beta = [sp.diff(zeta[0], X[1+k]) for k in range(3)]
print('3 slice: lam', lam.subs(t,0), 'A', A.subs(t,0), 'beta', [b_.subs(t,0) for b_ in beta], 'zeta0', zeta[0].subs(t,0))
random.seed(1); ok=True; n=0
f = sp.lambdify((t,x,y,z), [K[0],K[1]], 'math')
while n<20000:
    p=[random.uniform(-1,1) for _ in range(4)]
    r=(p[1]**2+p[2]**2+p[3]**2)**.5
    if abs(p[0])+r>=0.999: continue
    n+=1
    q=f(*p)
    if not (q[1] < 0.5-abs(q[0])): ok=False
print('4a O -> wedge on', n, 'samples:', ok)
th,ph = sp.symbols('th ph', real=True)
S = K.subs({t:0, x:sp.cos(th), y:sp.sin(th)*sp.cos(ph), z:sp.sin(th)*sp.sin(ph)})
print('4b sphere -> plane x1=1/2, t=0:', sp.simplify(S[1]-sp.Rational(1,2))==0, S[0]==0)
