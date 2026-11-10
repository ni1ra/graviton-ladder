"""PART 2b of w69-bclose (numerical, float64): the area side at FINITE eps, from the geometry alone.

For a metric g_eps = eta + eps h1 + eps^2 h2 (from out/p2_metrics.json) find the stationary point of
    F[Sigma] = A[boundary of Sigma] - 2 V[Sigma]
over discs  Sigma = { t = T(x), |x| <= 1 + Y(x/|x|) },  Y in real harmonics of degree <= L_Y, T in polynomials of
degree <= d_T (generic bases, no symmetry imposed).  Newton iteration from the coordinate location with the flat
Hessian; its seven Poincare zero modes are projected out and the gradient along them is reported.
Then  F*(eps) - F0 = c1 eps + c2 eps^2 + ...  and the area side is  Delta_2 = -2 c2  (to be compared with A_l).

  python3 p2_geo.py --plants     PLANT-FLAT and PLANT-SRC  -> out/p2_plants.json   (run first)
  python3 p2_geo.py              members M1..M6            -> out/p2_geo.json
"""
import json, os, sys, time
import numpy as np
import sympy as sp
from common import x, y, z, t, X, OUT, lap

PLANTS_ONLY = '--plants' in sys.argv
F0 = 4*np.pi - 8*np.pi/3
EPS_UNIT = [0.10, 0.15, 0.20, 0.25, 0.30]          # in units where the flux oint |K|^2 = 1


def harmonic_basis(L):
    """Real harmonic homogeneous polynomials of degree 0..L (nullspace of the Laplacian)."""
    out = []
    for d in range(L + 1):
        mons = [x**i*y**j*z**(d - i - j) for i in range(d + 1) for j in range(d + 1 - i)]
        cs = sp.symbols(f'c0:{len(mons)}')
        p = sum(c*m for c, m in zip(cs, mons))
        if d < 2:
            out += mons
            continue
        eqs = sp.Poly(lap(p), x, y, z).coeffs()
        Mx, _ = sp.linear_eq_to_matrix(eqs, cs)
        for v in Mx.nullspace():
            out.append(sp.expand(sum(a*m for a, m in zip(v, mons))))
    return out


def lam(expr):
    f = sp.lambdify((t, x, y, z), expr, 'numpy', cse=True)
    return lambda tt, xx, yy, zz: np.broadcast_to(np.asarray(f(tt, xx, yy, zz), dtype=float), xx.shape)


class Metric:
    def __init__(self, rec):
        self.h1 = [[sp.sympify(rec['h1'][m][n]) for n in range(4)] for m in range(4)]
        self.h2 = [[sp.sympify(rec['h2'][m][n]) for n in range(4)] for m in range(4)]
        self.f1 = {(m, n): lam(self.h1[m][n]) for m in range(4) for n in range(m, 4) if self.h1[m][n] != 0}
        self.f2 = {(m, n): lam(self.h2[m][n]) for m in range(4) for n in range(m, 4) if self.h2[m][n] != 0}

    def g(self, eps, tt, xx, yy, zz):
        G = np.zeros((4, 4) + xx.shape)
        G[0, 0] = -1.0
        for i in (1, 2, 3):
            G[i, i] = 1.0
        for src, c in ((self.f1, eps), (self.f2, eps*eps)):
            if c == 0.0:
                continue
            for (m, n), f in src.items():
                v = c*f(tt, xx, yy, zz)
                G[m, n] += v
                if m != n:
                    G[n, m] += v
        return G


class Geo:
    def __init__(self, LY, dT, quad):
        nth, nph, nrho = quad
        mu, wmu = np.polynomial.legendre.leggauss(nth)
        ph = 2*np.pi*np.arange(nph)/nph
        MU, PH = np.meshgrid(mu, ph, indexing='ij')
        st = np.sqrt(1 - MU**2)
        self.n = np.array([st*np.cos(PH), st*np.sin(PH), MU]).reshape(3, -1)
        self.dth = np.array([MU*np.cos(PH), MU*np.sin(PH), -st]).reshape(3, -1)
        self.dph = np.array([-np.sin(PH), np.cos(PH), 0*PH]).reshape(3, -1)
        self.ws = (wmu[:, None]*np.full(nph, 2*np.pi/nph)[None, :]).reshape(-1)
        rho, wr = np.polynomial.legendre.leggauss(nrho)
        self.rho, self.wr = (rho + 1)/2, wr/2
        hb = harmonic_basis(LY)
        nx, ny, nz = self.n
        ev = lambda e: np.broadcast_to(np.asarray(sp.lambdify((x, y, z), e, 'numpy')(nx, ny, nz), dtype=float), nx.shape)
        self.Bv = np.array([ev(p) for p in hb])
        self.Bg = np.array([[ev(sp.diff(p, s)) for s in X] for p in hb])           # (nY, 3, M)
        self.nY = len(hb)
        self.mons = [(i, j, k) for i in range(dT + 1) for j in range(dT + 1 - i) for k in range(dT + 1 - i - j)]
        self.nT = len(self.mons)
        self.npar = self.nY + self.nT
        self.dT = dT

    def Tval(self, tc, P):
        """T and grad T at points P (3, ...)."""
        pw = [[np.ones_like(P[a])] for a in range(3)]
        for a in range(3):
            for _ in range(self.dT):
                pw[a].append(pw[a][-1]*P[a])
        T = np.zeros_like(P[0])
        dTx = np.zeros((3,) + P[0].shape)
        for c, (i, j, k) in zip(tc, self.mons):
            if c == 0.0:
                continue
            T += c*pw[0][i]*pw[1][j]*pw[2][k]
            if i:
                dTx[0] += c*i*pw[0][i - 1]*pw[1][j]*pw[2][k]
            if j:
                dTx[1] += c*j*pw[0][i]*pw[1][j - 1]*pw[2][k]
            if k:
                dTx[2] += c*k*pw[0][i]*pw[1][j]*pw[2][k - 1]
        return T, dTx

    def F(self, p, met, eps, parts=False):
        yc, tc = p[:self.nY], p[self.nY:]
        Y = yc @ self.Bv
        gP = np.einsum('k,kam->am', yc, self.Bg)
        Yth = np.sum(self.dth*gP, axis=0)
        Yph = np.sum(self.dph*gP, axis=0)
        Xs = (1 + Y)*self.n
        e = [Yth*self.n + (1 + Y)*self.dth, Yph*self.n + (1 + Y)*self.dph]
        T, dTx = self.Tval(tc, Xs)
        G = met.g(eps, T, Xs[0], Xs[1], Xs[2])
        E = [np.concatenate([np.sum(dTx*ea, axis=0)[None], ea]) for ea in e]
        s = [[np.einsum('mnp,mp,np->p', G, E[a], E[b]) for b in range(2)] for a in range(2)]
        A = np.sum(self.ws*np.sqrt(s[0][0]*s[1][1] - s[0][1]**2))
        # volume
        P = self.rho[:, None, None]*Xs[None].transpose(1, 0, 2).reshape(1, 3, -1) if False else None
        Pv = np.array([self.rho[:, None]*Xs[a][None, :] for a in range(3)])            # (3, nrho, M)
        T, dTx = self.Tval(tc, Pv)
        G = met.g(eps, T, Pv[0], Pv[1], Pv[2])
        gam = np.empty((3, 3) + T.shape)
        for i in range(3):
            for j in range(3):
                gam[i, j] = G[i + 1, j + 1] + G[0, i + 1]*dTx[j] + G[0, j + 1]*dTx[i] + G[0, 0]*dTx[i]*dTx[j]
        det = (gam[0, 0]*(gam[1, 1]*gam[2, 2] - gam[1, 2]*gam[2, 1]) - gam[0, 1]*(gam[1, 0]*gam[2, 2] - gam[1, 2]*gam[2, 0])
               + gam[0, 2]*(gam[1, 0]*gam[2, 1] - gam[1, 1]*gam[2, 0]))
        V = np.sum(self.wr[:, None]*self.rho[:, None]**2*(self.ws*(1 + Y)**3)[None, :]*np.sqrt(det))
        return (A, V) if parts else A - 2*V

    def grad(self, p, met, eps, d=1e-5):
        g = np.empty(self.npar)
        for i in range(self.npar):
            q = p.copy(); q[i] += d
            fp = self.F(q, met, eps)
            q[i] -= 2*d
            g[i] = (fp - self.F(q, met, eps))/(2*d)
        return g

    def hess0(self, met, d=1e-3):
        n = self.npar
        p0 = np.zeros(n)
        f0 = self.F(p0, met, 0.0)
        f1 = np.empty((n, 2))
        for i in range(n):
            for sgn in (0, 1):
                q = p0.copy(); q[i] = d*(1 - 2*sgn)
                f1[i, sgn] = self.F(q, met, 0.0)
        H = np.empty((n, n))
        for i in range(n):
            H[i, i] = (f1[i, 0] - 2*f0 + f1[i, 1])/d**2
            for j in range(i):
                q = p0.copy(); q[i] = d; q[j] = d
                fpp = self.F(q, met, 0.0)
                q[i] = -d; q[j] = -d
                fmm = self.F(q, met, 0.0)
                H[i, j] = H[j, i] = (fpp + fmm - f1[i, 0] - f1[i, 1] - f1[j, 0] - f1[j, 1] + 2*f0)/(2*d**2)
        w, U = np.linalg.eigh(H)
        zero = np.abs(w) < 1e-4*np.max(np.abs(w))
        self.w, self.U, self.zero = w, U, zero
        self.pinv = (U[:, ~zero]/w[~zero]) @ U[:, ~zero].T
        return int(zero.sum()), float(np.max(np.abs(w[zero])) if zero.any() else 0.0), float(np.min(np.abs(w[~zero])))

    def locate(self, met, eps, tol=2e-9, itmax=40):
        p = np.zeros(self.npar)
        hist = []
        for it in range(itmax):
            g = self.grad(p, met, eps)
            gz = self.U[:, self.zero].T @ g
            gp = g - self.U[:, self.zero] @ gz
            hist.append(float(np.linalg.norm(gp)))
            if hist[-1] < tol:
                break
            p = p - self.pinv @ g
        return p, self.F(p, met, eps), hist, float(np.max(np.abs(gz))) if gz.size else 0.0


def run_member(name, rec, geo, sign_both=True, eps_list=EPS_UNIT):
    met = Metric(rec)
    flux = float(sp.Rational(rec['flux_over_pi']))*np.pi
    A_l = float(sp.Rational(rec['A_over_pi']))*np.pi
    rows = []
    es = []
    for eu in eps_list:
        for sg in ((1, -1) if sign_both else (1,)):
            er = sg*eu/np.sqrt(flux)
            p, Fs, hist, gz = geo.locate(met, er)
            Ffix = geo.F(np.zeros(geo.npar), met, er)
            rows.append({'eps_unit': sg*eu, 'F*-F0': Fs - F0, 'F_fixed-F0': Ffix - F0, 'newton_iterations': len(hist),
                         'final_gradient_norm': hist[-1], 'max|Y coeff|': float(np.max(np.abs(p[:geo.nY]))),
                         'max|T coeff|': float(np.max(np.abs(p[geo.nY:]))),
                         'zero-mode gradient / eps_unit^2': gz/eu**2})
            es.append(sg*eu)
    es = np.array(es)
    dF = np.array([r['F*-F0'] for r in rows])
    deg = min(6, len(es) - 2)
    Vm = np.array([es**k for k in range(1, deg + 1)]).T
    c, res, _, _ = np.linalg.lstsq(Vm, dF, rcond=None)
    fit_resid = float(np.max(np.abs(Vm @ c - dF)))
    out = {'rows': rows, 'fit_degree': deg, 'c1 (O(eps) coefficient, unit-flux eps)': float(c[0]),
           'c2': float(c[1]), 'fit max residual': fit_resid,
           'Delta_2 (raw units) = -2 c2 * flux': float(-2*c[1]*flux), 'A_l': A_l}
    if A_l > 0:
        out['Delta_2/A_l'] = float(-2*c[1]*flux/A_l)
        out['|Delta_2/A_l - 1|'] = abs(out['Delta_2/A_l'] - 1)
        out['c1 / (A_l/flux)'] = float(c[0]*flux/A_l)
    out['all newton runs converged'] = all(r['final_gradient_norm'] < 2e-9 for r in rows)
    return out


def make_geo(LY, dT, quad, flatmet):
    g = Geo(LY, dT, quad)
    nz, zmax, gap = g.hess0(flatmet)
    return g, {'L_Y': LY, 'd_T': dT, 'quadrature (n_theta, n_phi, n_rho)': list(quad), 'parameters': g.npar,
               'flat-Hessian zero modes': nz, 'largest zero-mode |eigenvalue|': zmax, 'smallest nonzero |eigenvalue|': gap}


if __name__ == '__main__':
    t0 = time.time()
    M = json.load(open(os.path.join(OUT, 'p2_metrics.json')))
    flat = Metric({'h1': [['0']*4]*4, 'h2': [['0']*4]*4})
    BASE, QB = (4, 3), (20, 40, 16)
    out = {}
    geo, info = make_geo(*BASE, QB, flat)
    out['geometry base'] = info
    out['CHECK flat Hessian has exactly 7 zero modes (Poincare)'] = (info['flat-Hessian zero modes'] == 7)
    out['F0 quadrature error'] = float(geo.F(np.zeros(geo.npar), flat, 0.0) - F0)
    print(info, 'F0 err', out['F0 quadrature error'], f'{time.time() - t0:.0f}s', flush=True)

    if PLANTS_ONLY:
        # PLANT-FLAT: exactly flat, corner displaced at O(eps)
        met = Metric(M['PLANT-FLAT'])
        rows = []
        for eps in (0.1, 0.2, 0.3):
            p, Fs, hist, gz = geo.locate(met, eps)
            Ffix = geo.F(np.zeros(geo.npar), met, eps)
            rows.append({'eps': eps, 'F*-F0 (located)': Fs - F0, 'F_fixed-F0 (coordinate corner and slice)': Ffix - F0,
                         'newton_iterations': len(hist), 'final_gradient_norm': hist[-1],
                         'max|T coeff|/eps': float(np.max(np.abs(p[geo.nY:]))/eps), 'max|Y coeff|': float(np.max(np.abs(p[:geo.nY])))})
            print('FLAT', rows[-1], flush=True)
        out['PLANT-FLAT'] = rows
        out['CHECK PLANT-FLAT: located |F* - F0| <= 1e-9 at the largest eps'] = abs(rows[-1]['F*-F0 (located)']) <= 1e-9
        out['CHECK PLANT-FLAT: fixed corner is caught (|F_fixed - F0| > 1e-5)'] = abs(rows[-1]['F_fixed-F0 (coordinate corner and slice)']) > 1e-5
        # PLANT-SRC
        r = run_member('PLANT-SRC', M['PLANT-SRC (M1, source x 11/10)'], geo)
        out['PLANT-SRC'] = r
        print('SRC', r['Delta_2/A_l'], flush=True)
        out['CHECK PLANT-SRC: Delta_2/A_l = 1.1 within 1e-4 (the broken constraint is seen)'] = abs(r['Delta_2/A_l'] - 1.1) <= 1e-4
        out['seconds'] = time.time() - t0
        checks = {k: v for k, v in out.items() if k.startswith('CHECK')}
        out['ALL_PASS'] = all(bool(v) for v in checks.values())
        json.dump(out, open(os.path.join(OUT, 'p2_plants.json'), 'w'), indent=1, default=float)
        print('p2_plants ALL_PASS:', out['ALL_PASS'])
        sys.exit(0)

    names = [k for k in M if k.startswith('M')]
    res = {}
    for nm in names:
        res[nm] = {'base': run_member(nm, M[nm], geo)}
        print(nm, 'base', res[nm]['base']['Delta_2/A_l'], f'{time.time() - t0:.0f}s', flush=True)
    geo2, info2 = make_geo(6, 4, QB, flat)
    out['geometry raised basis'] = info2
    for nm in names:
        res[nm]['raised basis'] = run_member(nm, M[nm], geo2, eps_list=EPS_UNIT)
        print(nm, 'basis+', res[nm]['raised basis']['Delta_2/A_l'], f'{time.time() - t0:.0f}s', flush=True)
    geo3, info3 = make_geo(*BASE, (26, 52, 22), flat)
    out['geometry raised quadrature'] = info3
    for nm in names:
        res[nm]['raised quadrature'] = run_member(nm, M[nm], geo3)
        print(nm, 'quad+', res[nm]['raised quadrature']['Delta_2/A_l'], f'{time.time() - t0:.0f}s', flush=True)
    summ = {}
    for nm in names:
        b, rb, rq = (res[nm][k] for k in ('base', 'raised basis', 'raised quadrature'))
        s = {'A_l': b['A_l'], 'Delta_2/A_l': b['Delta_2/A_l'], '|Delta_2/A_l - 1|': b['|Delta_2/A_l - 1|'],
             'change with raised basis (units of A_l)': abs(rb['Delta_2/A_l'] - b['Delta_2/A_l']),
             'change with raised quadrature (units of A_l)': abs(rq['Delta_2/A_l'] - b['Delta_2/A_l']),
             'c1 / (A_l/flux)': b['c1 / (A_l/flux)'],
             'max zero-mode gradient / eps_unit^2': max(abs(r['zero-mode gradient / eps_unit^2']) for r in b['rows']),
             'max |Y coeff| at eps_unit = 0.3': max(r['max|Y coeff|'] for r in b['rows'] if abs(r['eps_unit']) == 0.30),
             'max |T coeff| at eps_unit = 0.3': max(r['max|T coeff|'] for r in b['rows'] if abs(r['eps_unit']) == 0.30),
             'F_fixed - F* at eps_unit = 0.3': max(abs(r['F_fixed-F0'] - r['F*-F0']) for r in b['rows'] if abs(r['eps_unit']) == 0.30),
             'converged': b['all newton runs converged'] and rb['all newton runs converged'] and rq['all newton runs converged']}
        summ[nm] = s
        out[f'CHECK {nm}: |Delta_2/A_l - 1| <= 1e-4'] = s['|Delta_2/A_l - 1|'] <= 1e-4
        out[f'CHECK {nm}: converged; basis and quadrature changes <= 1e-5 A_l'] = (
            s['converged'] and s['change with raised basis (units of A_l)'] <= 1e-5 and s['change with raised quadrature (units of A_l)'] <= 1e-5)
        out[f'CHECK {nm}: O(eps) coefficient zero (|c1| <= 1e-6 in units of A_l/flux)'] = abs(s['c1 / (A_l/flux)']) <= 1e-6
    out['summary'] = summ
    out['detail'] = res
    out['seconds'] = time.time() - t0
    checks = {k: v for k, v in out.items() if k.startswith('CHECK')}
    out['ALL_PASS'] = all(bool(v) for v in checks.values())
    out['PART2_BIN'] = 'P2-PASS' if out['ALL_PASS'] else 'see failed checks'
    json.dump(out, open(os.path.join(OUT, 'p2_geo.json'), 'w'), indent=1, default=float)
    for nm in names:
        print(nm, summ[nm])
    print('p2_geo ALL_PASS:', out['ALL_PASS'], out['PART2_BIN'])
