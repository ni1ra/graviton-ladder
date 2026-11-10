"""w71-blind-final: the referee's own checks and plants, one block per proof attacked.

  SQ  Proposition SQ (THEOREM.md Section 4) on a finite-dimensional standard subspace, with exponential vectors
  M   the direct proof of (M): K6 for every l from K4, the conformal weight behind (l+1), the edge vector K 1_B
  C   the support-projection step of (C) on a matrix model whose channel has a non-trivial support
  B   the zero modes of Section 7 (B-d, B-e)
  R   the record lemmas: modflow Lemma 4 and Lemma 5, duality Lemmas A and B
A "plant" is a deliberately wrong theory; it must fail the check that guards it.
Reads wave69/final/scripts/lib.py; writes only wave71/blindfinal/out/p_blind.json.
"""
import json, os, sys
import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', '..', 'wave69', 'final', 'scripts'))
from lib import (x, y, z, X, r2, N, eps, I3, sphere_int, ball_int, lap, mat, mexp, is_zero, sym, curl, xcross,
                 div, contract, ric3, rlin, cotton, harmonics, magnetic, xvec)

checks, plants = [], []


def check(label, ok, detail=None):
    checks.append({'label': label, 'pass': bool(ok), 'detail': None if detail is None else str(detail)})
    if not ok:
        print('FAIL', label, detail)


def plant(label, caught, detail=None):
    plants.append({'label': label, 'caught': bool(caught), 'detail': None if detail is None else str(detail)})
    if not caught:
        print('PLANT NOT CAUGHT', label, detail)


# ---------------------------------------------------------------- SQ
rng = np.random.default_rng(71)
n = 3
lam = rng.uniform(0.2, 5.0, n)
U = np.linalg.qr(rng.normal(size=(2 * n, 2 * n)) + 1j * rng.normal(size=(2 * n, 2 * n)))[0]
D = np.concatenate([lam, 1 / lam])
ip = np.vdot                     # antilinear in the first slot


def inK(a):
    return U @ np.concatenate([a, np.sqrt(lam) * np.conj(a)])


def Jop(v):
    w = np.conj(U).T @ v
    return U @ np.concatenate([np.conj(w[n:]), np.conj(w[:n])])


def Dpow(v, p):
    return U @ (D**p * (np.conj(U).T @ v))


hs = [0.4 * inK(rng.normal(size=n) + 1j * rng.normal(size=n)) for _ in range(5)]
check('SQ toy: sample vectors are fixed by J Delta^(1/2)', max(abs(Jop(Dpow(h, .5)) - h).max() for h in hs) < 1e-12)
# (a) J H is in the symplectic complement
a_ok = max(abs(ip(Jop(h), k).imag) for h in hs for k in hs)
check('SQ (a): Im<J h, k> = 0 on H', a_ok < 1e-12, a_ok)
a_bad = max(abs(ip(1j * Jop(h), k).imag) for h in hs for k in hs)
plant('SQ (a) plant: i J H is not in the symplectic complement', a_bad > 1e-3, a_bad)
# (c) the Weyl two-point function is KMS at beta = -1: f(t - i) = <Delta^(it) k, h>
def f(h, k, zc):
    return ip(h, Dpow(k, 1j * zc))


c_ok = c_bad = 0.0
for t in (-0.8, 0.0, 0.5):
    for h in hs:
        for k in hs:
            target = ip(Dpow(k, 1j * t), h)
            c_ok = max(c_ok, abs(np.exp(-f(h, k, t - 1j)) - np.exp(-target)))
            c_bad = max(c_bad, abs(np.exp(-f(h, k, t + 1j)) - np.exp(-target)))
check('SQ (c): omega(V(h) sigma_{t-i}(V(k))) = omega(sigma_t(V(k)) V(h)), beta = -1', c_ok < 1e-9, c_ok)
plant('SQ (c) plant: beta = +1 (continue into the upper strip) fails', c_bad > 1e-3, c_bad)
# (d) Tomita operator on A Omega, A = sum c_j V(h_j):  ||Gamma(Delta^(1/2)) A Omega|| = ||A* Omega||, J = Gamma(-J_H)
cj = rng.normal(size=5) + 1j * rng.normal(size=5)
nrm = np.array([np.exp(-ip(h, h).real / 2) for h in hs])


def gram(vs_left, vs_right, cl, cr):
    return sum(np.conj(cl[i]) * cr[j] * nrm[i] * nrm[j] * np.exp(ip(vs_left[i], vs_right[j]))
               for i in range(5) for j in range(5))


lhs = gram([Dpow(h, .5) for h in hs], [Dpow(h, .5) for h in hs], cj, cj)
rhs = gram([-h for h in hs], [-h for h in hs], np.conj(cj), np.conj(cj))
check('SQ (d): ||Gamma(Delta^(1/2)) A Omega||^2 = ||A* Omega||^2', abs(lhs - rhs) < 1e-9 * abs(rhs), abs(lhs - rhs))
lhs_bad = gram([Dpow(h, -.5) for h in hs], [Dpow(h, -.5) for h in hs], cj, cj)
plant('SQ (d) plant: Delta^(-1/2) in place of Delta^(1/2)', abs(lhs_bad - rhs) > 1e-3 * abs(rhs), abs(lhs_bad - rhs))
w = hs[0] + 1j * hs[1]
d_ok = max(abs(np.exp(ip(w, -Jop(Dpow(h, .5)))) - np.exp(ip(w, -h))) for h in hs)
d_bad = max(abs(np.exp(ip(w, Jop(Dpow(h, .5)))) - np.exp(ip(w, -h))) for h in hs)
check('SQ (d): <Exp(w), Gamma(-J) Gamma(Delta^(1/2)) Exp(h)> = <Exp(w), Exp(-h)>', d_ok < 1e-9, d_ok)
plant('SQ (d) plant: Gamma(+J) is not the modular conjugation in this Weyl convention', d_bad > 1e-3, d_bad)

# ---------------------------------------------------------------- M
xv = xvec()
# K6 for every l: C f = 2 curl(Ric_3 f - (1/4) delta R_lin f) for every symmetric f, so K6 follows from K4 and Lap K = 0
fpoly = sp.Matrix(3, 3, lambda i, j: 0)
mons = [x**3 * y, y**2 * z**2, x * z**3, x * y * z**2, x**4, y**3 * z]
for idx, (i, j) in enumerate([(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]):
    fpoly[i, j] = fpoly[j, i] = mons[idx] + (idx + 1) * mons[(idx + 2) % 6]
schouten = ric3(fpoly) - sp.Rational(1, 4) * I3 * rlin(fpoly)
check('M K6: C f = 2 curl(Ric_3 f - delta R_lin f / 4) on a generic polynomial f', is_zero(mexp(cotton(fpoly) - 2 * curl(schouten))))
check('M K6: R_lin = trace of Ric_3', sp.expand(rlin(fpoly) - ric3(fpoly).trace()) == 0)
plant('M K6 plant: coefficient 1/2 in place of 1/4', not is_zero(mexp(cotton(fpoly) - 2 * curl(ric3(fpoly) - sp.Rational(1, 2) * I3 * rlin(fpoly)))))
for l in (6, 7):
    v, K = magnetic(harmonics(l)['m=l cos'])
    check(f'M K4 at l={l}: Ric_3(2NK) = (l+2)K', is_zero(mexp(ric3(mexp(2 * N * K))) - (l + 2) * K))
    check(f'M K3 at l={l}: 2 sym(x cross (-curl K)) = (l+1)K', is_zero(mexp(2 * sym(xcross(-curl(K)))) - (l + 1) * K))
v6, K6m = magnetic(harmonics(6)['m=l cos'])
check('M K6 at l=6 directly: C(r^2 K) = -2(l+2) curl K', is_zero(mexp(cotton(mexp(r2 * K6m))) + 2 * 8 * mexp(curl(K6m))))

# the integer l+1: the conformal Lie derivative of the stationary curvature solves the field equations only at weight 3
t = sp.symbols('t', real=True)
for l in (2, 3, 4):
    for name in ('m=l cos', 'm=0'):
        v, K = magnetic(harmonics(l)[name])
        B70 = mexp(curl(K))                      # wave70 sign: B = (1/2) curl pi = curl K; E = 0; static
        def lie(wt, sgn):
            dE = -sgn * 2 * sym(xcross(B70))     # beta = -x, lambda = -t, A = 0
            dB = -t * B70.applyfunc(lambda e: sum(X[i] * sp.diff(e, X[i]) for i in range(3))) - wt * t * B70
            return mexp(dE), mexp(dB)
        def solves(dE, dB):
            return is_zero(mexp(dE.applyfunc(lambda e: sp.diff(e, t)) - curl(dB))) and \
                   is_zero(mexp(dB.applyfunc(lambda e: sp.diff(e, t)) + curl(dE)))
        dE, dB = lie(3, 1)
        check(f'M weight l={l} {name}: delta_zeta(E,B) at weight 3 solves the field equations', solves(dE, dB))
        check(f'M weight l={l} {name}: on the slice delta E = (l+1)K, delta B = 0', is_zero(dE - (l + 1) * K) and is_zero(dB.subs(t, 0)))
        plant(f'M weight l={l} {name} plant: weight 2', not solves(*lie(2, 1)))
        plant(f'M weight l={l} {name} plant: weight 4', not solves(*lie(4, 1)))
        plant(f'M weight l={l} {name} plant: boost term with the opposite sign', not solves(*lie(3, -1)))

# the edge vector a K 1_B: same curvature inside (none), obeys the constraint, would shift S, and has infinite norm
from scipy.special import spherical_jn
from scipy.integrate import quad
for l in (2, 3):
    v, K = magnetic(harmonics(l)['m=l cos'])
    Kx = K * xv
    check(f'M edge l={l}: Ric_3(K) = 0 inside, so K 1_B has no curvature in the open ball', is_zero(mexp(ric3(K))))
    check(f'M edge l={l}: K.x is tangent and divergence-free, so R_lin[K 1_B] = 0 as a distribution',
          sp.expand((xv.T * Kx)[0]) == 0 and sp.expand(sum(sp.diff(Kx[i], X[i]) for i in range(3))) == 0)
    shift = sp.simplify(ball_int(sp.expand(contract(K, K))))
    plant(f'M edge l={l} plant: adding a K 1_B to the cut vector would change S by a int_B |K|^2 / 16G', shift != 0, shift)
    kk = 3.7
    r1 = quad(lambda r: r**(l + 1) * spherical_jn(l - 1, kk * r), 0, 1)[0]
    r2n = quad(lambda r: 0.5 * (1 - r * r) * r**(l + 1) * spherical_jn(l - 1, kk * r), 0, 1)[0]
    check(f'M edge l={l}: radial transform of the jump profile is j_l(k)/k', abs(r1 - spherical_jn(l, kk) / kk) < 1e-10)
    check(f'M edge l={l}: radial transform of the lapse profile N is j_(l+1)(k)/k^2', abs(r2n - spherical_jn(l + 1, kk) / kk**2) < 1e-10)
    def tail(fn, a, b):
        ks = np.linspace(a, b, int((b - a) * 40) + 1)
        return np.trapezoid(fn(ks), ks)
    jump = [tail(lambda k: k * spherical_jn(l, k)**2, 10.0**p, 10.0**(p + 1)) for p in (1, 2, 3)]
    cut = [tail(lambda k: spherical_jn(l + 1, k)**2 / k, 10.0**p, 10.0**(p + 1)) for p in (1, 2, 3)]
    check(f'M edge l={l}: H^(1/2) norm of the cut vector converges (decade tails fall by 100)', cut[2] < 2e-2 * cut[0], cut)
    plant(f'M edge l={l} plant: the edge vector K 1_B has finite norm (each decade adds (1/2) ln 10 instead)',
          all(abs(j - 0.5 * np.log(10)) < 0.05 for j in jump), jump)

# ---------------------------------------------------------------- C
def rand_state(d, seed):
    g = np.random.default_rng(seed)
    a = g.normal(size=(d, d)) + 1j * g.normal(size=(d, d))
    rho = a @ a.conj().T + 0.1 * np.eye(d)
    return rho / np.trace(rho).real


def logm_h(a):
    w_, v_ = np.linalg.eigh(a)
    return v_ @ np.diag(np.log(w_)) @ v_.conj().T


def relent(rho, sig):
    """S(rho || sig) on the common support (both restricted to the range of sig)."""
    w_, v_ = np.linalg.eigh(sig)
    keep = v_[:, w_ > 1e-12]
    r_, s_ = keep.conj().T @ rho @ keep, keep.conj().T @ sig @ keep
    return np.trace(r_ @ (logm_h(r_) - logm_h(s_))).real


d = 3
rho, sig = rand_state(d, 1), rand_state(d, 2)
e11 = np.diag([1.0, 0.0])
# M = M_3 (x) M_2, Phi(a (x) b) = a <e1, b e1>, iota(a) = a (x) 1.  The state omega o Phi has density rho (x) e11.
s_proj = np.kron(np.eye(d), e11)
Phi = lambda m: m.reshape(d, 2, d, 2)[:, 0, :, 0]
m_rand = rng.normal(size=(2 * d, 2 * d)) + 1j * rng.normal(size=(2 * d, 2 * d))
check('C toy: Phi(1 - s) = 0 and Phi(m) = Phi(s m s), s != 1',
      abs(Phi(np.eye(2 * d) - s_proj)).max() < 1e-14 and abs(Phi(m_rand) - Phi(s_proj @ m_rand @ s_proj)).max() < 1e-12)
check('C toy: Phi o iota = id', abs(Phi(np.kron(rho, np.eye(2))) - rho).max() < 1e-14)
S_A = relent(rho, sig)
S_M = relent(np.kron(rho, e11), np.kron(sig, e11))
check('C toy: S_M = S_A when H1-H3 hold and the support is not 1', abs(S_M - S_A) < 1e-10, (S_M, S_A))
pinch = lambda a: np.diag(np.diag(a))
S_noH3 = relent(np.kron(pinch(rho), e11), np.kron(pinch(sig), e11))
plant('C plant: H3 fails (Phi followed by a pinching, no left inverse): S_M drops strictly', S_noH3 < S_A - 1e-3, (S_noH3, S_A))
best = 0.0
for sd in range(40):
    q = np.linalg.qr(np.random.default_rng(100 + sd).normal(size=(d, d)) + 1j * np.random.default_rng(200 + sd).normal(size=(d, d)))[0]
    best = max(best, relent(rho, q @ sig @ q.conj().T))
plant('C plant: H1 fails (a different channel for the vacuum): S_M can exceed S_A', best > S_A + 1e-3, (best, S_A))

# ---------------------------------------------------------------- B (zero modes)
half = sp.Rational(1, 2)
four_pi = 4 * sp.pi
# flat quadratic form on the Poincare modes, from the exact moved ball
# translation by a along x: sigma = a n1 + a^2 (n1^2 - 1)/2; boost by b: tau = b n1, s = b x, sigma = b^2 n1^2 / 2
n1sq = sphere_int(x**2)                  # 4 pi / 3
grad_n1 = sphere_int(1 - x**2)           # |grad_S n1|^2 = 1 - n1^2
Q_trans = half * grad_n1 - n1sq
Q_boost = -half * grad_n1 + ball_int(sp.Integer(1))
check('B zero modes: Q vanishes on a translation of the ball', sp.simplify(Q_trans) == 0)
check('B zero modes: Q vanishes on a boost of the ball (tau = b n1, s = b x)', sp.simplify(Q_boost) == 0)
plant('B plant: tau coefficient -1 instead of -1/2 in Q breaks boost invariance', sp.simplify(-grad_n1 + ball_int(sp.Integer(1))) != 0)
jj = sp.symbols('j', integer=True, nonnegative=True)
check('B zero modes: sigma form (j-1)(j+2)/2 from j(j+1)/2 - 1; tau form -j(j-1)/2 from -j(j+1)/2 + j',
      sp.expand(half * jj * (jj + 1) - 1 - half * (jj - 1) * (jj + 2)) == 0 and sp.expand(-half * jj * (jj + 1) + jj + half * jj * (jj - 1)) == 0)
# B-e, translations: the eps^2 coefficient at the ball moved by a is -(1/2) int_B N |K(x+a)|^2; its derivative is -(1/2) int_B x_a |K|^2
for l in (2, 3):
    for name in ('m=l cos', 'm=0', 'mixed'):
        v, K = magnetic(harmonics(l)[name])
        K2 = sp.expand(contract(K, K))
        grads = [sp.simplify(ball_int(sp.expand(X[a] * K2))) for a in range(3)]
        direct = [sp.simplify(ball_int(sp.expand(N * sp.diff(K2, X[a])))) for a in range(3)]
        check(f'B-e l={l} {name}: translation gradient -(1/2) int_B x_a |K|^2 = 0', all(g == 0 for g in grads))
        check(f'B-e l={l} {name}: it equals int_B N d_a |K|^2 (parts, N = 0 on S)', all(sp.simplify(direct[a] - grads[a]) == 0 for a in range(3)))
K_mix = magnetic(harmonics(2)['mixed'])[1] + magnetic(harmonics(3)['mixed'])[1]
g_mix = [sp.simplify(ball_int(sp.expand(X[a] * contract(K_mix, K_mix)))) for a in range(3)]
plant('B-e plant: a mixture of l = 2 and l = 3 drives the translation mode (single l is needed)', any(g != 0 for g in g_mix), g_mix)
phi = harmonics(2)['m=0']
Ke = mat(lambda i, j: sp.diff(phi, X[i], X[j]))
drive = sphere_int(sp.expand(x * (Ke.trace() - (xv.T * Ke * xv)[0])))
drive0 = sphere_int(sp.expand((Ke.trace() - (xv.T * Ke * xv)[0])**2))
plant('B-d plant: an electric datum K = Hess phi has K_nn != 0, so L_K does not vanish', drive0 != 0, drive0)

# ---------------------------------------------------------------- R (record lemmas)
# modflow Lemma 4, items 3-4, on the toy: e^(2 pi A) = Delta, S = J Delta^(1/2)
us = hs
l4a = max(abs(ip(v_, Dpow(u_, 1.0)) - ip(u_, v_)) for u_ in us for v_ in us)
check('R Lemma 4 item 3: <v, e^(2 pi A) u> = <u, v> on K', l4a < 1e-10, l4a)
xs = [us[0] + 1j * us[1], us[2] - 2j * us[3], us[4] + 0.3j * us[0]]
Sop = lambda xx: Jop(Dpow(xx, .5))
l4b = max(abs(ip(Dpow(Sop(yv), .5), Dpow(Sop(xv_), .5)) - ip(xv_, yv)) for xv_ in xs for yv in xs)
check('R Lemma 4 item 4: J0 = e^(pi A) S is antilinear and isometric on K + iK', l4b < 1e-10, l4b)
l4c = max(abs(ip(Dpow(Sop(yv), .55), Dpow(Sop(xv_), .55)) - ip(xv_, yv)) for xv_ in xs for yv in xs)
plant('R Lemma 4 plant: rate 1.1 (e^(1.1 pi A) S) is not isometric', l4c > 1e-3, l4c)
l4d = max(abs(ip(v_, Dpow(u_, -1.0)) - ip(u_, v_)) for u_ in us for v_ in us)
plant('R Lemma 4 plant: opposite direction (Delta^(-1) at -2 pi i)', l4d > 1e-3, l4d)
graph = [ip(xx, xx).real + ip(Sop(xx), Sop(xx)).real for xx in xs]
parts = [2 * (ip(us[0], us[0]).real + ip(us[1], us[1]).real), 2 * (ip(us[2], us[2]).real + 4 * ip(us[3], us[3]).real),
         2 * (ip(us[4], us[4]).real + 0.09 * ip(us[0], us[0]).real)]
check('R Lemma 4 item 5: graph norm of S on K + iK is sqrt 2 times the norm of K (+) K', max(abs(a - b) for a, b in zip(graph, parts)) < 1e-10)
# modflow Lemma 5: damping on the closed strip, |grad psi| = |q|, Lambda_{i pi} = j
g5 = np.random.default_rng(5)
worst_in, worst_grad, best_out = 0.0, 0.0, 0.0
for _ in range(4000):
    kq = g5.normal(size=3) * 3
    q = np.array([np.linalg.norm(kq), *kq])
    x1 = abs(g5.normal()) + 1e-3
    tt = g5.uniform(-1, 1) * x1
    pt = np.array([tt, x1, g5.normal(), g5.normal()])
    th = g5.uniform(0, np.pi)
    def psi(pnt, th_):
        yv = np.array([pnt[0] * np.cos(th_) + 1j * pnt[1] * np.sin(th_), pnt[1] * np.cos(th_) + 1j * pnt[0] * np.sin(th_), pnt[2], pnt[3]])
        return q[0] * yv[0] - q[1:] @ yv[1:]
    worst_in = max(worst_in, abs(np.exp(1j * psi(pt, th))))
    best_out = max(best_out, abs(np.exp(1j * psi(pt, -th))))
    gr = np.array([q[0] * np.cos(th) - 1j * q[1] * np.sin(th), 1j * q[0] * np.sin(th) - q[1] * np.cos(th), -q[2], -q[3]])
    worst_grad = max(worst_grad, abs(np.vdot(gr, gr).real - q @ q))
check('R Lemma 5: |e^(i p.Lambda_{-z} x)| <= 1 on the closed strip 0 <= theta <= pi, x in the wedge', worst_in <= 1 + 1e-12, worst_in)
check('R Lemma 5: |grad psi|^2 = |q|^2 (Euclidean), uniformly in theta', worst_grad < 1e-9, worst_grad)
plant('R Lemma 5 plant: continuing into the other half-plane (theta < 0) is unbounded', best_out > 10, best_out)
th = sp.symbols('theta', real=True)
Lam = sp.Matrix([[sp.cos(th), sp.I * sp.sin(th), 0, 0], [sp.I * sp.sin(th), sp.cos(th), 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
check('R Lemma 5: Lambda_{i pi} = diag(-1,-1,1,1)', Lam.subs(th, sp.pi) == sp.diag(-1, -1, 1, 1))
plant('R Lemma 5 plant: at theta = pi/2 the point is complex, so no real test field and no reality', any(e.has(sp.I) for e in Lam.subs(th, sp.pi / 2)))
# duality Lemma A, third bullet: a compactly supported TT field has zero mean (so finite H^(-3/2) norm)
bump = (1 - r2)**6
fA = mexp(bump * sp.Matrix([[x * y, z**2, 0], [z**2, -x * y + x, y * z], [0, y * z, -x]]))
TA = mexp(cotton(fA))
check('R Lemma A: T = C f is symmetric, traceless, transverse', is_zero(TA - TA.T) and sp.expand(TA.trace()) == 0 and is_zero(div(TA)))
check('R Lemma A: its mean over the ball is zero (Fourier transform vanishes at k = 0)', all(sp.simplify(ball_int(TA[i, j])) == 0 for i in range(3) for j in range(3)))
plant('R Lemma A plant: a compactly supported field that is not transverse has non-zero mean, hence infinite H^(-3/2) norm',
      any(sp.simplify(ball_int(sp.expand(bump * e))) != 0 for e in sp.Matrix([[1 + x, 0, 0], [0, -1, 0], [0, 0, -x]])))
# duality Lemma B: (P) uses H^2_c = 0.  In a solid torus a toroidal field has a non-zero period against d(angle)/rho.
rho_, zz = sp.symbols('rho zeta', real=True)
s_, al = sp.symbols('s alpha', nonnegative=True)
def disc_int(expr):
    e = sp.expand(expr.subs({rho_: 2 + s_ * sp.cos(al), zz: s_ * sp.sin(al)}) * s_)
    return sp.simplify(sp.integrate(sp.integrate(e, (al, 0, 2 * sp.pi)), (s_, 0, 1)))
gtor = (1 - (rho_ - 2)**2 - zz**2)**3
period_w = 2 * sp.pi * disc_int(gtor)                         # int w . (phi-hat / rho) dV, w = g phi-hat
a_rho, a_z = gtor * (rho_ * zz + 1), gtor * (rho_**2 - zz)      # any compactly supported potential
period_curl = 2 * sp.pi * disc_int(sp.diff(a_rho, zz) - sp.diff(a_z, rho_))
check('R Lemma B: a curl of a potential supported in the torus has zero period', period_curl == 0, period_curl)
plant('R Lemma B plant: in a solid torus a divergence-free toroidal field is not such a curl', period_w != 0, period_w)

# modflow Step 6 (not rerun by run_all.sh): the special conformal map carries zeta to the boost field of the wedge x1 < 1/2 - |t|
tq, xq, yq, zq = sp.symbols('t_q x_q y_q z_q', real=True)
pt4 = sp.Matrix([tq, xq, yq, zq])
sq = tq**2 - xq**2 - yq**2 - zq**2
sig = 1 + 2 * xq - sq                              # b = e_1: b.x = -x1, b^2 = -1
kap = (pt4 - sp.Matrix([0, 1, 0, 0]) * sq) / sig
zeta = sp.Matrix([(1 - tq**2 - xq**2 - yq**2 - zq**2) / 2, -tq * xq, -tq * yq, -tq * zq])
push = kap.jacobian(pt4) * zeta
def boost_field(apex):
    return sp.Matrix([-(kap[1] - apex), -kap[0], 0, 0])
check('R Step 6: D kappa . zeta = X_W o kappa, X_W = -(x1 - 1/2) d_t - t d_1', all(sp.simplify(e) == 0 for e in (push - boost_field(sp.Rational(1, 2)))))
plant('R Step 6 plant: wedge with apex at x1 = 1 instead of 1/2', any(sp.simplify(e) != 0 for e in (push - boost_field(1))))

ok = all(c['pass'] for c in checks) and all(p['caught'] for p in plants)
out = os.path.join(HERE, '..', 'out')
os.makedirs(out, exist_ok=True)
json.dump({'script': 'p_blind', 'ALL_PASS': ok, 'n_checks': len(checks), 'n_plants': len(plants), 'checks': checks, 'plants': plants},
          open(os.path.join(out, 'p_blind.json'), 'w'), indent=1)
print('p_blind', 'ALL_PASS' if ok else 'FAILED', len(checks), 'checks,', len(plants), 'plants')
sys.exit(0 if ok else 1)
