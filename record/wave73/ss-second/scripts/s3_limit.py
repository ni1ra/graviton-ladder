"""P2, C1/C2 coefficients, C3(a), sign of (B.26).  G = 1.
What is sampled: the EXACT eps-regulated vacuum kernels of s2 (closed forms verified there), integrated in closed form over
the direction of point 2 against the spin-weight test function (mbar1.m2)^s, then smeared in (u1,u2) on a uniform grid
(trapezoid, N points on [0,1]) with smooth compactly supported test functions. Compared with
  graviton  j = d_u^2 S_m = -2 u Psi_0 :  <j j^bar> -> 48 G delta/(u12 - i0)^4      [ (B.26) with a MINUS sign ]
  Maxwell   d_u calA_m = u phi_0       :  -> -(1/4pi) delta/(u12 - i0)^2           [ (B.13) ]
  scalar    d_u Phi (D = 4, 6)         :  -> -(1/4pi) delta/(u12 - i0)^2           [ (3.4) ]
whose smeared values are computed independently in Fourier space:
  1/(u - i0)^2 = -int_0^inf w e^{iwu} dw,   1/(u - i0)^4 = (1/6) int_0^inf w^3 e^{iwu} dw.
"""
import sys
import numpy as np
import sympy as sp

N = int(sys.argv[1]) if len(sys.argv) > 1 else 4001
EPS = [float(e) for e in sys.argv[2].split(',')] if len(sys.argv) > 2 else [0.08, 0.04, 0.02, 0.01]
ug = np.linspace(0.0, 1.0, N)
h = ug[1] - ug[0]

def bump(u, a, b):
    out = np.zeros_like(u)
    m = (u > a) & (u < b)
    out[m] = np.exp(-(b - a)**2/(4*(u[m] - a)*(b - u[m])) + 1.0)
    return out

K1 = bump(ug, 0.15, 0.75)
K2 = bump(ug, 0.30, 0.95)*(1 + ug)
K3 = np.zeros_like(ug); mm = ug < 0.9
K3[mm] = ug[mm]**2*np.exp(-0.2/(0.9 - ug[mm]))          # O(u^2) at the tip, smooth cut-off before the edge
tests = {'K1,K1': (K1, K1), 'K1,K2': (K1, K2), 'K3,K3 (tip u^2)': (K3, K3), 'K3,K1': (K3, K1)}

# ---------- closed-form angular integrals ----------
def q_grav(w):
    """w^5 * int_0^1 (1-s)^4/(s+w)^5 ds"""
    out = np.empty_like(w)
    big = np.abs(w) >= 3
    wb = w[big]
    winv = 1/wb
    acc = np.zeros_like(wb)
    for n in range(79, -1, -1):          # Horner in -1/w (plain powers overflow to NaN for |w| > 1e4: run05 K3 rows)
        acc = acc*(-winv) + 1.0/(n + 5)
    out[big] = acc
    ws = w[~big]
    Q = np.log((1 + ws)/ws)
    from math import comb
    for k in range(4):
        Q += comb(4, k)*(-1)**k*(1 + ws)**(4 - k)*((1 + ws)**(k - 4) - ws**(k - 4))/(k - 4)
    out[~big] = ws**5*Q
    return out

def J_grav(U1, U2, eps, kappa2=32*np.pi):
    """<j(1) jbar(2)> integrated over direction 2;  j = -2 u Psi_0;  kappa2 = 32 pi G (plant: 16 pi G)."""
    a = 4*U1*U2
    z = 2j*eps*(U1 - U2) + eps**2
    w = z/np.where(a == 0, 1e-300, a)
    # kernel = (kappa2/4) * 16 W'''' eps^4 c^2,  W'''' = 24/(4 pi^2 sigma^5);  angular measure 4 pi ds
    pref = (kappa2/4)*16*24/(4*np.pi**2)*eps**4*4*np.pi
    return 4*U1*U2*pref*q_grav(w)/z**5

s, a_, z_ = sp.symbols('s a z')
_m = sp.integrate((1 - s)**2/(a_*s + z_)**3, (s, 0, 1))
m_fun = sp.lambdify((a_, z_), sp.simplify(_m), 'numpy')
def M_maxw(U1, U2, eps):
    a = 4*U1*U2; z = 2j*eps*(U1 - U2) + eps**2
    return U1*U2*(8*eps**2/np.pi)*m_fun(a + 0j, z)

x1, x2, e_ = sp.symbols('u1 u2 e', positive=True)
zz = 2*sp.I*e_*(x1 - x2) + e_**2
aa = 4*x1*x2
G4 = (sp.log(aa + zz) - sp.log(zz))/(4*sp.pi)
G6 = 4*(x1*x2)**2/sp.pi*sp.integrate(s*(1 - s)/(a_*s + z_)**2, (s, 0, 1)).subs({a_: aa, z_: zz})
S4 = sp.lambdify((x1, x2, e_), sp.diff(G4, x1, x2), 'numpy')
S6 = sp.lambdify((x1, x2, e_), sp.simplify(sp.diff(G6, x1, x2)), 'numpy')

def smear(kfun, Ka, Kb, eps, chunk=400):
    tot = 0j
    idx = np.nonzero(Ka)[0]
    jdx = np.nonzero(Kb)[0]
    u2 = ug[jdx][None, :]
    for i0 in range(0, len(idx), chunk):
        ii = idx[i0:i0 + chunk]
        vals = kfun(ug[ii][:, None] + 0*u2, u2 + 0*ug[ii][:, None], eps)
        tot += np.sum(Ka[ii][:, None]*vals*Kb[jdx][None, :])
    return tot*h*h

_FT = {}
def _ft(K):
    key = id(K)
    if key not in _FT:
        w = np.linspace(0, 1500, 15001)
        ia = np.nonzero(K)[0]
        _FT[key] = np.array([np.sum(K[ia]*np.exp(1j*wi*ug[ia]))*h for wi in w])
    return _FT[key]

def fourier(Ka, Kb, power, damp=0.0):
    """int_0^inf w^power conj(Khat_a(w)) Khat_b(w) exp(-damp w) dw.
    Convention: 1/(u - i e)^2 = -int_0^inf w e^{-iwu - e w} dw, 1/(u - i e)^4 = (1/6) int_0^inf w^3 e^{-iwu - e w} dw,
    so int int Ka(u1) Kb(u2) /(u12 - i e)^n picks conj(Khat_a) Khat_b and damp = e."""
    w = np.linspace(0, 1500, 15001)
    Fa = _ft(Ka); Fb = _ft(Kb)
    integrand = w**power*np.conj(Fa)*Fb*np.exp(-damp*w)
    return np.trapezoid(integrand, w), abs(integrand[-1])/max(abs(integrand).max(), 1e-300)

print('N =', N, ' h =', h)
ok = True
print('\n== pointwise: J_eps / [48/(u12 - i eps/2)^4] ==')
for (p1, p2) in [(0.5, 0.5), (0.5, 0.52), (0.3, 0.7), (0.05, 0.06), (0.9, 0.2)]:
    row = []
    for eps in EPS + [1e-3, 1e-4]:
        Jv = J_grav(np.array([p1]), np.array([p2]), eps)[0]
        row.append(Jv/(48/(p1 - p2 - 0.5j*eps)**4))
    print(' (u1,u2)=(%.2f,%.2f):' % (p1, p2), ' '.join('%.5f%+.5fi' % (r.real, r.imag) for r in row))

def run(label, kfun, target_fun, tsts):
    """err_reg: against the comparator regulated at the same eps (i0 -> i eps/2); err_0: against the eps = 0 target;
    err_R: Richardson extrapolant 2 A(eps_min) - A(2 eps_min) against the eps = 0 target."""
    global ok
    print('\n== %s ==' % label)
    for name, (Ka, Kb) in tsts.items():
        tgt, tail = target_fun(Ka, Kb, 0.0)
        A = [smear(kfun, Ka, Kb, eps) for eps in EPS]
        ereg = [abs(Ai - target_fun(Ka, Kb, eps/2)[0])/abs(tgt) for Ai, eps in zip(A, EPS)]
        e0 = [abs(Ai - tgt)/abs(tgt) for Ai in A]
        eR = abs(2*A[-1] - A[-2] - tgt)/abs(tgt)
        print(' %-16s target(eps=0) = %+.8e%+.3ei  (Fourier tail %.0e)' % (name, tgt.real, tgt.imag, tail))
        print('    eps      :', ' '.join('%9.5f' % e for e in EPS))
        print('    err_reg  :', ' '.join('%9.2e' % e for e in ereg))
        print('    err_0    :', ' '.join('%9.2e' % e for e in e0), '   err_R = %.2e' % eR)
        good = ereg[-1] < 1e-2 and (ereg[-1] < ereg[0] or ereg[0] < 1e-3) and all(e0[i+1] < e0[i] for i in range(len(e0)-1)) and eR < 1e-2
        print('    ', 'ok' if good else 'NOT ok')
        ok &= good

run('P2 graviton D=4: smeared <j jbar> vs 48 G/(u12-i0)^4', J_grav,
    lambda Ka, Kb, d: (lambda r: (48/6*r[0], r[1]))(fourier(Ka, Kb, 3, d)), tests)
F1 = bump(ug, 0.1, 0.8); F2 = bump(ug, 0.25, 0.9)*(1 + ug)
ft = {'F1,F1': (F1, F1), 'F1,F2': (F1, F2)}
tgt2 = lambda Ka, Kb, d: (lambda r: (r[0]/(4*np.pi), r[1]))(fourier(Ka, Kb, 1, d))
run('C2 Maxwell D=4: smeared <dA_m dA_mbar> vs -(1/4pi)/(u12-i0)^2', M_maxw, tgt2, ft)
run('C1 scalar D=4: smeared <dPhi dPhi> vs -(1/4pi)/(u12-i0)^2', lambda a, b, e: S4(a, b, e + 0j), tgt2, ft)
run('C1 scalar D=6: smeared <dPhi dPhi> vs -(1/4pi)/(u12-i0)^2', lambda a, b, e: S6(a, b, e + 0j), tgt2, ft)

print('\n== C3(a) plant: graviton normalisation 16 pi G ==')
tgt, _ = fourier(K1, K1, 3, EPS[-1]/2); tgt = 48/6*tgt
A = smear(lambda a, b, e: J_grav(a, b, e, kappa2=16*np.pi), K1, K1, EPS[-1])
rel = abs(A - tgt)/abs(tgt)
print(' rel.err = %.3f  -> plant caught: %s' % (rel, rel > 0.4))
ok &= rel > 0.4

print('\n== sign of (B.26): int int F F/(u12 - i eps)^2 for real F ==')
for eps in (0.02, 0.01):
    val = smear(lambda a, b, e: 1/(a - b - 1j*e)**2, F1, F1, eps)
    print(' eps = %.2f : %+.6f%+.1ei   (Fourier: %+.6f)' % (eps, val.real, val.imag, -fourier(F1, F1, 1, eps)[0].real))
print(' => printed sign (+1/(32 pi^2 G)) gives omega(gamma[f]^2) < 0; the state is (B.26) with a minus sign.')
print('\nS3 PASS' if ok else '\nS3 FAIL')
