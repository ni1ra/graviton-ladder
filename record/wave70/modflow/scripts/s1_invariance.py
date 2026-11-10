"""S1 (PREREG section 3): exact checks on the helicity-2 two-point function.
a. structure of W; field equations on the first slot.
b. Lemma 1: (delta_X (x) 1 + 1 (x) delta_X) W = 0 as a rational identity, 11 generators.
c. [delta_X, delta_Y] = delta_[X,Y], 105 pairs.
Plants: weight 2, weight 4, boost sign +1.
"""
import itertools
import json
import sys
import time
from multiprocessing import Pool

import sympy as sp

sys.path.insert(0, __file__.rsplit('/', 1)[0])
from common import (T, X1, X2, X3, XI, S, R3, eps, build_W, killing_fields, delta_matrix_part)

xs = sp.symbols('x0 x1 x2 x3')
ys = tuple(xs[m] - XI[m] for m in range(4))
OUT = __file__.rsplit('/', 2)[0] + '/out/s1_out.json'

Wn = None


def dnum(N, mu, n=5):
    """numerator of d_mu (N / s^n) over s^(n+1)."""
    return sp.expand(S * sp.diff(N, XI[mu]) - n * N * sp.diff(S, XI[mu]))


def check_structure():
    res = {}
    idx = list(itertools.product(R3, repeat=4))
    res['BB_equals_EE'] = all(sp.expand(Wn[1][1][k] - Wn[0][0][k]) == 0 for k in idx)
    res['BE_equals_minus_EB'] = all(sp.expand(Wn[1][0][k] + Wn[0][1][k]) == 0 for k in idx)
    ok_sym = ok_tr = ok_div = True
    for a in (0, 1):
        for b in (0, 1):
            for i, j, k, l in idx:
                if sp.expand(Wn[a][b][i, j, k, l] - Wn[a][b][j, i, k, l]) != 0:
                    ok_sym = False
                if sp.expand(Wn[a][b][i, j, k, l] - Wn[a][b][i, j, l, k]) != 0:
                    ok_sym = False
            for k, l in itertools.product(R3, repeat=2):
                if sp.expand(sum(Wn[a][b][i, i, k, l] for i in R3)) != 0:
                    ok_tr = False
                if sp.expand(sum(Wn[a][b][k, l, i, i] for i in R3)) != 0:
                    ok_tr = False
                for j in R3:
                    if sp.expand(sum(dnum(Wn[a][b][i, j, k, l], 1 + i) for i in R3)) != 0:
                        ok_div = False
                    if sp.expand(sum(dnum(Wn[a][b][k, l, i, j], 1 + i) for i in R3)) != 0:
                        ok_div = False
    res['symmetric_each_pair'] = ok_sym
    res['traceless_each_pair'] = ok_tr
    res['transverse_each_pair'] = ok_div
    # exchange symmetry W^{ab}_{ij,kl}(xi) vs W^{ba}_{kl,ij}(-xi): the rational function is even/odd so that
    # the commutator is supported on the cone
    flip = {T: -T, X1: -X1, X2: -X2, X3: -X3}
    res['W_ab(xi)_equals_W_ba(-xi)'] = all(
        sp.expand(Wn[a][b][i, j, k, l] - Wn[b][a][k, l, i, j].subs(flip, simultaneous=True)) == 0
        for a in (0, 1) for b in (0, 1) for i, j, k, l in idx)
    # field equations on the first slot: d_T W^{E b} = curl W^{B b}, d_T W^{B b} = -curl W^{E b}
    ok1 = ok2 = True
    for b in (0, 1):
        for i, j, k, l in idx:
            c1 = sum(eps(i, m, n) * dnum(Wn[1][b][n, j, k, l], 1 + m) for m in R3 for n in R3)
            c2 = sum(eps(i, m, n) * dnum(Wn[0][b][n, j, k, l], 1 + m) for m in R3 for n in R3)
            if sp.expand(dnum(Wn[0][b][i, j, k, l], 0) - c1) != 0:
                ok1 = False
            if sp.expand(dnum(Wn[1][b][i, j, k, l], 0) + c2) != 0:
                ok2 = False
    res['Edot_eq_curlB_first_slot'] = ok1
    res['Bdot_eq_minus_curlE_first_slot'] = ok2
    # wave equation
    res['box_W_zero'] = all(
        sp.expand(sp.cancel(sum((1 if m == 0 else -1) * sp.diff(Wn[0][0][k] / S**5, XI[m], 2) for m in range(4)))) == 0
        for k in [(0, 0, 0, 0), (0, 1, 0, 1), (0, 1, 1, 2)])
    return res


def invariance_defect(args):
    """Number of non-zero components of the numerator (over s^6) of
    (delta_X at x on slot 1 + delta_X at y on slot 2) W(x - y)."""
    name, w, s = args
    Xx = killing_fields(xs)[name]
    Xy = tuple(sp.sympify(c).subs(dict(zip(xs, ys)), simultaneous=True) for c in Xx)
    Xy_sym = killing_fields(xs)[name]
    nonzero = 0
    total = 0
    sub_y = dict(zip(xs, ys))
    for b in (0, 1):
        for k, l in itertools.product(R3, repeat=2):
            # first slot matrices for fixed (b, k, l)
            E = [[Wn[0][b][i, j, k, l] for j in R3] for i in R3]
            B = [[Wn[1][b][i, j, k, l] for j in R3] for i in R3]
            dE, dB = delta_matrix_part(Xx, xs, E, B, w=w, s=s)
            first = (dE, dB)
            for a in (0, 1):
                for i, j in itertools.product(R3, repeat=2):
                    # second slot matrices for fixed (a, i, j)
                    E2 = [[Wn[a][0][i, j, m, n] for n in R3] for m in R3]
                    B2 = [[Wn[a][1][i, j, m, n] for n in R3] for m in R3]
                    dE2, dB2 = delta_matrix_part(Xy_sym, xs, E2, B2, w=w, s=s)
                    second = (dE2, dB2)[b][k][l]
                    second = sp.sympify(second).subs(sub_y, simultaneous=True)
                    N = Wn[a][b][i, j, k, l]
                    transport = sum((Xx[m] - Xy[m]) * dnum(N, m) for m in range(4))
                    tot = sp.expand(transport + S * (first[a][i][j] + second))
                    total += 1
                    if tot != 0:
                        nonzero += 1
    return name, w, s, nonzero, total


def lie_algebra_check():
    F = killing_fields(xs)
    names = list(F)
    Es = sp.symbols('e0:6')
    Bs = sp.symbols('b0:6')

    def symm(v):
        return [[v[0], v[1], v[2]], [v[1], v[3], v[4]], [v[2], v[4], v[5]]]

    def flat(M):
        return [M[0][0], M[0][1], M[0][2], M[1][1], M[1][2], M[2][2]]
    vars_ = list(Es) + list(Bs)

    def Mmat(Xf):
        dE, dB = delta_matrix_part(Xf, xs, symm(Es), symm(Bs))
        out = flat(dE) + flat(dB)
        return sp.Matrix([[sp.diff(sp.expand(o), v) for v in vars_] for o in out])
    M = {n: Mmat(F[n]) for n in names}

    def bracket(Xf, Yf):
        return tuple(sp.expand(sum(Xf[m] * sp.diff(Yf[n], xs[m]) - Yf[m] * sp.diff(Xf[n], xs[m]) for m in range(4)))
                     for n in range(4))
    bad = []
    for n1, n2 in itertools.combinations(names, 2):
        Z = bracket(F[n1], F[n2])
        MZ = Mmat(Z)
        lhs = (sum((F[n1][m] * sp.diff(M[n2], xs[m]) - F[n2][m] * sp.diff(M[n1], xs[m]) for m in range(4)),
                   sp.zeros(12, 12)) + M[n1] * M[n2] - M[n2] * M[n1])
        # (delta_X delta_Y - delta_Y delta_X) phi = [X,Y].d phi + (X.dM_Y - Y.dM_X + [M_X, M_Y]) phi
        if sp.simplify(lhs - MZ) != sp.zeros(12, 12):
            bad.append((n1, n2))
    return {'pairs': len(list(itertools.combinations(names, 2))), 'failing_pairs': bad}


def init():
    global Wn
    Wn = build_W()


if __name__ == '__main__':
    t0 = time.time()
    init()
    out = {'a_structure': check_structure()}
    print('a', out['a_structure'], round(time.time() - t0, 1), flush=True)
    gens = ['D', 'N1', 'N2', 'N3', 'J12', 'J13', 'J23', 'K0', 'K1', 'K2', 'K3']
    SGN = int(sys.argv[1]) if len(sys.argv) > 1 else +1   # corrected sign; registered was -1
    jobs = [(g, 3, SGN) for g in gens]
    plants = [(g, 2, SGN) for g in ('D', 'K0', 'K1')] + [(g, 4, SGN) for g in ('D', 'K0', 'K1')] + \
             [(g, 3, -SGN) for g in ('N1', 'K0', 'K1')]
    with Pool(8, initializer=init) as pool:
        res = pool.map(invariance_defect, jobs + plants)
    out['b_invariance'] = [dict(generator=r[0], w=r[1], s=r[2], nonzero_components=r[3], components=r[4])
                           for r in res[:len(jobs)]]
    out['plants'] = [dict(generator=r[0], w=r[1], s=r[2], nonzero_components=r[3], components=r[4])
                     for r in res[len(jobs):]]
    for r in res:
        print(r, flush=True)
    out['c_lie_algebra'] = lie_algebra_check()
    print('c', out['c_lie_algebra'])
    out['b_pass'] = all(r['nonzero_components'] == 0 for r in out['b_invariance'])
    out['plants_caught'] = {
        'w=2': any(r['nonzero_components'] > 0 for r in out['plants'] if r['w'] == 2),
        'w=4': any(r['nonzero_components'] > 0 for r in out['plants'] if r['w'] == 4),
        'opposite_boost_sign': any(r['nonzero_components'] > 0 for r in out['plants'] if r['s'] == -SGN and r['w'] == 3)}
    out['boost_sign_used'] = SGN
    out['seconds'] = round(time.time() - t0, 1)
    json.dump(out, open(OUT, 'w'), indent=1, default=str)
    print('b_pass', out['b_pass'], 'plants', out['plants_caught'], out['seconds'])
