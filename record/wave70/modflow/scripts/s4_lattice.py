"""S4 (PREREG section 3): lattice modular Hamiltonian of the ball in a radial sector,
from the vacuum correlators alone, against the conformal-flow generator.

Chain (a = 1): sites j = 1..ntot, K_jj = 2 + l(l+1)/j^2, K_{j,j+-1} = -1, Dirichlet walls.
Vacuum: X = K^{-1/2}/2, P = K^{1/2}/2.  Ball: sites 1..n, R = n + 1/2.
Modular Hamiltonian H = (pi^T M pi + u^T N u)/2 with
  M = X^{1/2} O (eps/nu) O^T X^{1/2},  N = X^{-1/2} O (eps nu) O^T X^{-1/2},
  X^{1/2} P X^{1/2} = O nu^2 O^T,  eps = log((nu + 1/2)/(nu - 1/2)).
Run with the scratch venv python (mpmath + gmpy2).  usage: s4_lattice.py l n dps
"""
import json
import sys
import time

import mpmath
from mpmath import mp, mpf
from mpmath.matrices.eigen_symmetric import tridiag_eigen


def correlators(l, n, ntot):
    L = l * (l + 1)
    d = [mpf(2) + mpf(L) / (j * j) for j in range(1, ntot + 1)]
    e = [mpf(-1)] * ntot
    z = mp.zeros(n, ntot)
    for i in range(n):
        z[i, i] = 1
    tridiag_eigen(mp, d, e, z)
    X = mp.zeros(n, n)
    P = mp.zeros(n, n)
    sq = [mp.sqrt(x) for x in d]
    for i in range(n):
        for j in range(i, n):
            sx = mp.fsum(z[i, k] * z[j, k] / sq[k] for k in range(ntot)) / 2
            sp_ = mp.fsum(z[i, k] * z[j, k] * sq[k] for k in range(ntot)) / 2
            X[i, j] = X[j, i] = sx
            P[i, j] = P[j, i] = sp_
    return X, P


def modular(X, P):
    n = X.rows
    ex, Qx = mp.eigsy(X)
    Xh = Qx * mp.diag([mp.sqrt(v) for v in ex]) * Qx.T
    Xmh = Qx * mp.diag([1 / mp.sqrt(v) for v in ex]) * Qx.T
    S = Xh * P * Xh
    S = (S + S.T) / 2
    e2, O = mp.eigsy(S)
    nu = [mp.sqrt(v) for v in e2]
    epsv = [mp.log((v + mpf(1) / 2) / (v - mpf(1) / 2)) for v in nu]
    A = Xh * O
    Bm = Xmh * O
    M = A * mp.diag([epsv[k] / nu[k] for k in range(n)]) * A.T
    N = Bm * mp.diag([epsv[k] * nu[k] for k in range(n)]) * Bm.T
    return M, N, nu, epsv


def stats(M, N, l, n, weight='minus', temp=2, Ltarget=None, epsr=0):
    """T1, T2 for the three registered profiles against a target.
    weight: 'minus' -> (R^2 - (1+epsr) r^2)/2R, 'plus' -> (R^2 + r^2)/2R; temp: 2 -> 2 pi, 1 -> pi."""
    R = mpf(n) + mpf(1) / 2
    L = l * (l + 1) if Ltarget is None else Ltarget

    def NR(r):
        r = mpf(r)
        if weight == 'minus':
            return (R * R - (1 + mpf(epsr)) * r * r) / (2 * R)
        return (R * R + r * r) / (2 * R)
    out = []
    for c in (mpf('0.3'), mpf('0.5'), mpf('0.7')):
        g = [mp.exp(-((j - c * R) ** 2) / (2 * (mpf('0.12') * R) ** 2)) for j in range(0, n + 2)]  # g[j] = g(r_j)
        gin = mp.matrix(g[1:n + 1])
        qM = (gin.T * M * gin)[0]
        qN = (gin.T * N * gin)[0]
        tM = temp * mp.pi * mp.fsum(NR(j) * g[j] ** 2 for j in range(1, n + 1))
        gg = [mpf(0)] + g[1:]            # chain has u_0 = 0
        tN = temp * mp.pi * (mp.fsum(NR(j + mpf(1) / 2) * (gg[j + 1] - gg[j]) ** 2 for j in range(0, n + 1))
                             + mp.fsum(NR(j) * L * g[j] ** 2 / (j * j) for j in range(1, n + 1)))
        out.append((float(qM / tM), float(qN / tN)))
    return out


def run(l, n, dps, ratio=6):
    mp.dps = dps
    t0 = time.time()
    X, P = correlators(l, n, ratio * n)
    M, N, nu, epsv = modular(X, P)
    res = {'l': l, 'n': n, 'ntot': ratio * n, 'dps': dps,
           'eps_max': float(max(epsv)), 'eps_min': float(min(epsv)),
           'min_nu_minus_half_log10': float(mp.log10(min(nu) - mpf(1) / 2)),
           'target': stats(M, N, l, n),
           'plant_P1_weight_plus': stats(M, N, l, n, weight='plus'),
           'plant_P2_temp_pi': stats(M, N, l, n, temp=1),
           'plant_P3_L0': stats(M, N, l, n, Ltarget=0),
           'resolution': {str(e): stats(M, N, l, n, epsr=e) for e in (0.02, 0.05, 0.1, 0.2, -0.05, -0.1)},
           # raw profile of the momentum block: diagonal and row sums against 2 pi N_R
           'M_diag_over_2piN': [float(M[j, j] / (2 * mp.pi * ((mpf(n) + mpf(1) / 2) ** 2 - (j + 1) ** 2) / (2 * n + 1)))
                                for j in range(n)],
           'M_rowsum_over_2piN': [float(mp.fsum(M[j, k] for k in range(n)) / (2 * mp.pi * ((mpf(n) + mpf(1) / 2) ** 2 - (j + 1) ** 2) / (2 * n + 1)))
                                  for j in range(n)],
           'M_sample': [[mp.nstr(M[i, j], 15) for j in range(min(n, 6))] for i in range(min(n, 6))],
           'N_sample': [[mp.nstr(N[i, j], 15) for j in range(min(n, 6))] for i in range(min(n, 6))],
           'seconds': None}
    res['seconds'] = round(time.time() - t0, 1)
    return res, M, N


if __name__ == '__main__':
    l, n, dps = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    ratio = int(sys.argv[4]) if len(sys.argv) > 4 else 6
    res, M, N = run(l, n, dps, ratio)
    # stability: second precision
    if len(sys.argv) > 5:
        dps2 = int(sys.argv[5])
        res2, M2, N2 = run(l, n, dps2, ratio)
        mp.dps = 30
        dM = max(abs(M[i, j] - M2[i, j]) for i in range(n) for j in range(n))
        dN = max(abs(N[i, j] - N2[i, j]) for i in range(n) for j in range(n))
        sc = max(abs(M2[i, j]) for i in range(n) for j in range(n))
        scN = max(abs(N2[i, j]) for i in range(n) for j in range(n))
        res['stability'] = {'dps2': dps2, 'max_abs_dM': float(dM), 'max_abs_dN': float(dN),
                            'rel_dM': float(dM / sc), 'rel_dN': float(dN / scN)}
        res['second_precision_target'] = res2['target']
    out = __file__.rsplit('/', 2)[0] + '/out/s4_l%d_n%d%s.json' % (l, n, '' if ratio == 6 else '_ratio%d' % ratio)
    json.dump(res, open(out, 'w'), indent=1)
    print(json.dumps({k: res[k] for k in ('l', 'n', 'dps', 'eps_max', 'min_nu_minus_half_log10', 'target',
                                           'plant_P1_weight_plus', 'plant_P2_temp_pi', 'plant_P3_L0', 'seconds')}))
    if 'stability' in res:
        print(res['stability'])
