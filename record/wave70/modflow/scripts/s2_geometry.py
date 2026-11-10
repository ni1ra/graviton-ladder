"""S2 (PREREG section 3): the conformal map from the double cone to the wedge, and finite covariance.

kappa_tau = special conformal transformation with parameter tau*b, b = e_1 (R = 1):
  kappa_tau(x) = (x - tau b x^2) / sigma_tau(x),  sigma_tau = 1 - 2 tau b.x + tau^2 b^2 x^2,
the time-tau flow of the conformal Killing field K1 of common.killing_fields.
(The PREREG wrote b = -e_1 in the route but registered the plane x_1 = 1/2 and X_W = -(x_1 - 1/2) d_t - t d_1,
which belong to b = +e_1 in this convention; b = +e_1 is used.)

a. exact identities.   b. numeric finite covariance with the cocycle from the flow ODE.
"""
import itertools
import json
import random
import sys

import sympy as sp
import mpmath
from mpmath import mp, mpf

sys.path.insert(0, __file__.rsplit('/', 1)[0])
from common import XI, S, R3, build_W, killing_fields, delta_matrix_part

OUT = __file__.rsplit('/', 2)[0] + '/out/s2_out_fix.json'
out = {}

# ---------------- a. exact ----------------
t, a, b, c, tau = sp.symbols('t a b c tau', real=True)
xs = (t, a, b, c)


def mink(u, v):
    return u[0] * v[0] - u[1] * v[1] - u[2] * v[2] - u[3] * v[3]


bvec = (0, 1, 0, 0)


def kappa(x, tt=1):
    x2 = mink(x, x)
    sig = 1 - 2 * tt * mink(bvec, x) + tt**2 * mink(bvec, bvec) * x2
    return tuple((x[m] - tt * bvec[m] * x2) / sig for m in range(4)), sig


ys = sp.symbols('y0 y1 y2 y3', real=True)
kx, sx = kappa(xs)
ky, sy = kappa(ys)
dxy = tuple(kx[m] - ky[m] for m in range(4))
out['a1_interval_identity'] = sp.simplify(mink(dxy, dxy) - mink(tuple(xs[m] - ys[m] for m in range(4)),
                                                                tuple(xs[m] - ys[m] for m in range(4))) / (sx * sy)) == 0
out['a2_sigma_is_minus_interval_to_p'] = sp.expand(sx + mink((t, a + 1, b, c), (t, a + 1, b, c))) == 0
# sphere -> plane x1 = 1/2
k1_on_sphere = kx[1].subs(t, 0).subs(c**2, 1 - a**2 - b**2)
out['a3_sphere_to_plane_x1_half'] = sp.simplify(k1_on_sphere - sp.Rational(1, 2)) == 0
out['a3_t0_slice_preserved'] = sp.simplify(kx[0].subs(t, 0)) == 0
# push-forward of zeta is the wedge boost
zeta = (sp.Rational(1, 2) * (1 - t**2 - a**2 - b**2 - c**2), -t * a, -t * b, -t * c)
Jac = sp.Matrix([[sp.diff(kx[m], xs[n]) for n in range(4)] for m in range(4)])
push = Jac * sp.Matrix(zeta)
XW = sp.Matrix([-(kx[1] - sp.Rational(1, 2)), -kx[0], 0, 0])
out['a4_pushforward_zeta_is_wedge_boost'] = sp.simplify(push - XW) == sp.zeros(4, 1)
# one-parameter group: d/dtau kappa_tau = K1(kappa_tau), kappa_0 = id
kt, st = kappa(xs, tau)
K1 = killing_fields(xs)['K1']
K1_at = [sp.sympify(K1[m]).subs(dict(zip(xs, kt)), simultaneous=True) for m in range(4)]
out['a5_kappa_tau_is_flow_of_K1'] = all(sp.simplify(sp.diff(kt[m], tau) - K1_at[m]) == 0 for m in range(4))
# zeta = (P0 - K0)/2
F = killing_fields(xs)
out['a6_zeta_is_(P0-K0)/2'] = all(sp.expand(zeta[m] - (F['P0'][m] - F['K0'][m]) / 2) == 0 for m in range(4))
# sigma_tau > 0 on O for tau in [0,1]: sampled (the report gives the argument)
random.seed(7)
mn = 1e9
for _ in range(20000):
    while True:
        p = [random.uniform(-1, 1) for _ in range(4)]
        r = (p[1]**2 + p[2]**2 + p[3]**2) ** 0.5
        if abs(p[0]) + r < 1:
            break
    tt = random.uniform(0, 1)
    val = 1 + 2 * tt * p[1] - tt * tt * (p[0]**2 - r * r)
    mn = min(mn, val)
out['a7_min_sigma_tau_on_O_sampled'] = mn
# image of O lies in the wedge x1 < 1/2 - |t| (sampled) and the inverse maps the wedge into O
ok = True
kap_num = sp.lambdify(xs, kx, 'math')
for _ in range(5000):
    while True:
        p = [random.uniform(-1, 1) for _ in range(4)]
        r = (p[1]**2 + p[2]**2 + p[3]**2) ** 0.5
        if abs(p[0]) + r < 1:
            break
    q = kap_num(*p)
    if not (q[1] < 0.5 - abs(q[0])):
        ok = False
out['a8_image_of_O_in_wedge_sampled'] = ok
print({k: v for k, v in out.items()}, flush=True)

# ---------------- b. numeric finite covariance ----------------
mp.dps = 30
Wn = build_W()
# orthonormal basis of S^2_0, then (E, B) -> 10 components
s2, s6 = mp.sqrt(2), mp.sqrt(6)
basis = [
    [[1 / s2, 0, 0], [0, -1 / s2, 0], [0, 0, 0]],
    [[1 / s6, 0, 0], [0, 1 / s6, 0], [0, 0, -2 / s6]],
    [[0, 1 / s2, 0], [1 / s2, 0, 0], [0, 0, 0]],
    [[0, 0, 1 / s2], [0, 0, 0], [1 / s2, 0, 0]],
    [[0, 0, 0], [0, 0, 1 / s2], [0, 1 / s2, 0]],
]
idx4 = list(itertools.product(R3, repeat=4))
Wfun = [[sp.lambdify(XI, [Wn[a_][b_][k] for k in idx4], 'mpmath') for b_ in (0, 1)] for a_ in (0, 1)]
Sfun = sp.lambdify(XI, S, 'mpmath')


def W10(xi):
    s = Sfun(*xi)
    M = mp.zeros(10, 10)
    for a_ in (0, 1):
        for b_ in (0, 1):
            vals = dict(zip(idx4, Wfun[a_][b_](*xi)))
            for p in range(5):
                for q in range(5):
                    M[5 * a_ + p, 5 * b_ + q] = mp.fsum(basis[p][i][j] * basis[q][k][l] * vals[i, j, k, l]
                                                        for i, j, k, l in idx4) / s**5
    return M


def Mred_fun(w):
    """10x10 matrix of the zeroth-order part of delta_{K1} in the basis, as a function of the point."""
    Es = sp.symbols('e0:5')
    Bs = sp.symbols('b0:5')
    bs = [sp.Matrix(3, 3, lambda i, j: sp.nsimplify(basis[p][i][j], [sp.sqrt(2), sp.sqrt(6)])) for p in range(5)]
    E = sum((Es[p] * bs[p] for p in range(5)), sp.zeros(3, 3))
    B = sum((Bs[p] * bs[p] for p in range(5)), sp.zeros(3, 3))
    dE, dB = delta_matrix_part(killing_fields(xs)['K1'], xs, E.tolist(), B.tolist(), w=w)
    dE, dB = sp.Matrix(dE), sp.Matrix(dB)
    rows = []
    for D in (dE, dB):
        for p in range(5):
            comp = sp.expand(sum(bs[p][i, j] * D[i, j] for i in R3 for j in R3))
            rows.append([sp.diff(comp, v) for v in list(Es) + list(Bs)])
    return sp.lambdify(xs, sp.Matrix(rows), 'mpmath')


def kap_tau(x, tt):
    x2 = x[0]**2 - x[1]**2 - x[2]**2 - x[3]**2
    sig = 1 + 2 * tt * x[1] - tt * tt * x2
    return [(x[m] - tt * (1 if m == 1 else 0) * x2) / sig for m in range(4)]


def cocycle(x, Mf, nsteps):
    """C with dC/dtau = M(kappa_tau x)^T C, C(0) = 1, tau in [0,1]; RK4."""
    C = mp.eye(10)
    h = mpf(1) / nsteps

    def f(tt, C_):
        return Mf(*kap_tau(x, tt)).T * C_
    for n in range(nsteps):
        tt = n * h
        k1 = f(tt, C)
        k2 = f(tt + h / 2, C + h / 2 * k1)
        k3 = f(tt + h / 2, C + h / 2 * k2)
        k4 = f(tt + h, C + h * k3)
        C = C + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return C


def cocycle_taylor(x, Mf, order=70):
    """Fix round (post hoc): the RK4 of the registered run was not accurate enough near the singular
    cone (out/s2_out_as_registered.json).  Taylor-series integration of dC/dtau = M(kappa_tau x)^T C:
    kappa_tau(x) = (x - tau b x^2)/sigma, sigma = 1 + 2 tau x_1 - tau^2 x^2 is rational in tau, M is
    affine in the point, so the series of M and of C follow from exact recurrences.  Steps are 0.3 of
    the distance to the complex zeros of sigma."""
    x2 = x[0]**2 - x[1]**2 - x[2]**2 - x[3]**2
    roots = mp.polyroots([-x2, 2 * x[1], 1]) if x2 != 0 else [-1 / (2 * x[1])]
    M0 = Mf(mpf(0), mpf(0), mpf(0), mpf(0))
    Mlin = [Mf(*[mpf(1) if n == m else mpf(0) for n in range(4)]) - M0 for m in range(4)]   # M = M0 + sum x^m Mlin[m]
    C = mp.eye(10)
    t0 = mpf(0)
    nsteps = 0
    while t0 < 1:
        rad = min(abs(r - t0) for r in roots)
        h = min(mpf('0.3') * rad, 1 - t0)
        # series in u = tau - t0 of sigma and of the numerators
        sig = [1 + 2 * t0 * x[1] - t0**2 * x2, 2 * x[1] - 2 * t0 * x2, -x2]
        num = [[x[m] - t0 * (x2 if m == 1 else 0), -(x2 if m == 1 else 0)] for m in range(4)]
        # 1/sigma series
        inv = [1 / sig[0]]
        for k in range(1, order + 1):
            acc = sig[1] * inv[k - 1] + (sig[2] * inv[k - 2] if k >= 2 else 0)
            inv.append(-acc / sig[0])
        xser = [[num[m][0] * inv[k] + (num[m][1] * inv[k - 1] if k >= 1 else 0) for k in range(order + 1)] for m in range(4)]
        Aser = []
        for k in range(order + 1):
            A = (M0 if k == 0 else mp.zeros(10, 10))
            for m in range(4):
                A = A + xser[m][k] * Mlin[m]
            Aser.append(A.T)
        Cs = [C]
        for k in range(order):
            acc = mp.zeros(10, 10)
            for j in range(k + 1):
                acc = acc + Aser[j] * Cs[k - j]
            Cs.append(acc / (k + 1))
        Cn = mp.zeros(10, 10)
        hp = mpf(1)
        for k in range(order + 1):
            Cn = Cn + hp * Cs[k]
            hp = hp * h
        tail = max(abs(Cs[order][i, j]) for i in range(10) for j in range(10)) * h**order
        C = Cn
        t0 = t0 + h
        nsteps += 1
    return C, tail


def cocycle_extrap(x, Mf):
    C, tail = cocycle_taylor(x, Mf)
    C2, _ = cocycle_taylor(x, Mf, order=50)
    return C, C2


def rand_pair(kind):
    while True:
        p = [mpf(random.uniform(-1, 1)) for _ in range(4)]
        q = [mpf(random.uniform(-1, 1)) for _ in range(4)]
        ok = all(abs(v[0]) + mp.sqrt(v[1]**2 + v[2]**2 + v[3]**2) < mpf('0.9') for v in (p, q))
        d = [p[m] - q[m] for m in range(4)]
        s = d[0]**2 - d[1]**2 - d[2]**2 - d[3]**2
        if ok and ((kind == 'spacelike' and s < -0.05) or (kind == 'timelike' and s > 0.05)):
            return p, q


random.seed(11)
Mf3, Mf4 = Mred_fun(3), Mred_fun(4)
rows = []
worst = 0
worst_plant = 1e9
for kind in ('spacelike', 'timelike'):
    for n in range(10):
        p, q = rand_pair(kind)
        Cp, Cp_lo = cocycle_extrap(p, Mf3)
        Cq, Cq_lo = cocycle_extrap(q, Mf3)
        kp, kq = kap_tau(p, 1), kap_tau(q, 1)
        lhs = Cp.T * W10([kp[m] - kq[m] for m in range(4)]) * Cq
        rhs = W10([p[m] - q[m] for m in range(4)])
        sc = max(abs(rhs[i, j]) for i in range(10) for j in range(10))
        err = max(abs(lhs[i, j] - rhs[i, j]) for i in range(10) for j in range(10)) / sc
        ode_err = max(abs(Cp[i, j] - Cp_lo[i, j]) for i in range(10) for j in range(10))
        Dp, _ = cocycle_extrap(p, Mf4)
        Dq, _ = cocycle_extrap(q, Mf4)
        lhs4 = Dp.T * W10([kp[m] - kq[m] for m in range(4)]) * Dq
        err4 = max(abs(lhs4[i, j] - rhs[i, j]) for i in range(10) for j in range(10)) / sc
        rows.append({'kind': kind, 'rel_err': float(err), 'ode_order70_vs_50_diff': float(ode_err), 'plant_w4_rel_err': float(err4)})
        worst = max(worst, float(err))
        worst_plant = min(worst_plant, float(err4))
        print(rows[-1], flush=True)
out['b_pairs'] = rows
out['b_worst_rel_err'] = worst
out['b_pass_1e-12'] = worst < 1e-12
out['b_plant_w4_min_rel_err'] = worst_plant
out['b_plant_caught'] = worst_plant > 1e-3
json.dump(out, open(OUT, 'w'), indent=1, default=str)
print('worst', worst, 'plant min', worst_plant)
