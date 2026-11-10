"""w52-gravlattice theory side: data from potentials (h = sym curl Y, pi = inc Z, constraints exact),
continuum Weyl fields, chain fields psi = x.C.x per parity, multipole split, and the predictions
  S_T1  = 2 pi sum_l c_l [ int_B N_rho (psidot^2 + |grad psi|^2) + (1/rho) int_B psi^2 ]   (= 4 pi sum c_l Q_rho[u])
in canonical units q = h/sqrt(32 pi G), [q_ij, p_kl] = i Pi_ijkl, H = 1/2 sum_9 (p^2 + |dq|^2)
(S_phys = c_l Q[u_h]/8G = 4 pi c_l Q[u_q]).  Plants: surface sigma = 2 sum_l c_l oint psi_l^2 dA, c_l -> 1/L,
c_l -> 4/L^2, Born coefficient 0.  The chain field is used here only, on the prediction side."""
import sympy as sp, numpy as np, json, sys
from scipy.special import sph_harm_y
from itertools import combinations_with_replacement as cwr

X = sp.symbols('x y z', real=True)
x, y, z = X
r2 = x * x + y * y + z * z
EPS = sp.LeviCivita
TCOMP = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]


def symm(M):
    return sp.Matrix(3, 3, lambda i, j: (M[i, j] + M[j, i]) / 2)


def curl(T):  # (curl T)_ij = eps_ikl d_k T_lj
    return sp.Matrix(3, 3, lambda i, j: sum(EPS(i, k, l) * sp.diff(T[l, j], X[k]) for k in range(3) for l in range(3)))


def symcurl(T):
    return symm(curl(T))


def inc(T):  # eps_ikl eps_jmn d_k d_m T_ln
    return sp.Matrix(3, 3, lambda i, j: sp.expand(sum(EPS(i, k, l) * EPS(j, m, n) * sp.diff(T[l, n], X[k], X[m])
                                                      for k in range(3) for l in range(3) for m in range(3) for n in range(3)
                                                      if EPS(i, k, l) != 0 and EPS(j, m, n) != 0)))


def ric(h):
    tr = h.trace()
    return sp.Matrix(3, 3, lambda i, j: sp.expand(sp.Rational(1, 2) * sum(sp.diff(h[j, k], X[k], X[i]) + sp.diff(h[i, k], X[k], X[j]) - sp.diff(h[i, j], X[k], X[k]) for k in range(3)) - sp.Rational(1, 2) * sp.diff(tr, X[i], X[j])))


def hess(f):
    return sp.Matrix(3, 3, lambda i, j: sp.diff(f, X[i], X[j]))


def monos(deg):
    out = []
    for d in range(deg + 1):
        for c in cwr(range(3), d):
            m = sp.Integer(1)
            for v in c:
                m *= X[v]
            out.append(m)
    return out


def solve_potential(target, op, deg):
    """find symmetric polynomial potential T of degree <= deg with op(T) = target (least-norm particular solution)."""
    ms = monos(deg)
    cs = sp.symbols('a0:%d' % (6 * len(ms)))
    T = sp.zeros(3, 3)
    n = 0
    for (i, j) in TCOMP:
        e = sum(cs[n + k] * ms[k] for k in range(len(ms)))
        n += len(ms)
        T[i, j] = e; T[j, i] = e
    R = op(T) - target
    eqs = []
    for i in range(3):
        for j in range(i, 3):
            eqs += sp.Poly(sp.expand(R[i, j]), *X).coeffs()
    sol = sp.solve(eqs, cs, dict=True)
    assert sol, 'no potential'
    sol = sol[0]
    Ts = T.subs(sol).subs({c: 0 for c in cs})
    assert sp.simplify(op(Ts) - target) == sp.zeros(3, 3)
    return Ts


def fields(Y, Zp):
    q = symcurl(Y)
    p = inc(Zp)
    E = ric(q)
    Kt = p / 2 - p.trace() / 4 * sp.eye(3)
    B = -symcurl(Kt)
    Edot = -symcurl(B)
    Bdot = symcurl(E)
    xv = sp.Matrix(X)
    psi = {'even': sp.expand((xv.T * E * xv)[0]), 'odd': sp.expand((xv.T * B * xv)[0])}
    psid = {'even': sp.expand((xv.T * Edot * xv)[0]), 'odd': sp.expand((xv.T * Bdot * xv)[0])}
    chk = dict(trE=sp.simplify(E.trace()), divE=[sp.simplify(sum(sp.diff(E[i, j], X[i]) for i in range(3))) for j in range(3)],
               trq=sp.simplify(q.trace()), ddq=sp.simplify(sum(sp.diff(q[i, j], X[i], X[j]) for i in range(3) for j in range(3))),
               divp=[sp.simplify(sum(sp.diff(p[i, j], X[i]) for i in range(3))) for j in range(3)])
    return dict(q=q, p=p, E=E, B=B, psi=psi, psid=psid, chk=chk)


# ------------------------------------------------------------ data (unit ball); all but the controls have psi != 0 on S
P2 = (x * x - y * y) / 2 + x * z + sp.Rational(1, 3) * y * z          # harmonic quadratic
P2b = x * y + (y * y - z * z) / 2                                     # harmonic quadratic (B-type)
P3 = x ** 3 - 3 * x * y * y + 3 * x * y * z                          # harmonic cubic
P4 = x ** 4 - 6 * x * x * y * y + y ** 4 + 6 * x * y * z * z - x ** 3 * y - x * y ** 3              # harmonic quartic
P4 = sp.expand(P4)


def harmonic_check(f):
    return sp.simplify(sum(sp.diff(f, v, v) for v in X)) == 0


def data_defs():
    import functools
    zero = sp.zeros(3, 3)
    D = {}
    opE = lambda T: ric(symcurl(T))
    opB = lambda T: -symcurl(inc(T) / 2 - inc(T).trace() / 4 * sp.eye(3))
    YE2 = solve_potential(hess(P2), opE, 3)
    ZB2 = solve_potential(hess(P2b), opB, 3)
    YE3 = solve_potential(hess(P3), opE, 4)
    YE4 = solve_potential(hess(P4), opE, 5)
    D['g1_E2'] = (YE2, zero)                       # uniform tidal E, l = 2
    D['g2_B2'] = (zero, ZB2)                       # uniform B, l = 2
    D['g3_E3'] = (YE3, zero)                       # E linear, l = 3
    D['g4_E23'] = (YE2 + YE3, zero)                # mixed l = 2 + 3 (c_l test)
    D['g5_E4'] = (YE4, zero)                       # l = 4, boundary-heavy
    D['g6_EB'] = (YE2, ZB2)                        # both parities, l = 2
    bump = (1 - r2) ** 4
    D['c1_E'] = (sp.expand(bump * YE2), zero)      # control: vanishes on S
    D['c2_B'] = (zero, sp.expand(bump * ZB2))      # control: vanishes on S
    return D


# ------------------------------------------------------------ multipole projection and predictions
NT, NP, NR = 28, 56, 40
ct, wt = np.polynomial.legendre.leggauss(NT)
TH, PH = np.meshgrid(np.arccos(ct), 2 * np.pi * np.arange(NP) / NP, indexing='ij')
WS = (wt[:, None] * np.full(NP, 2 * np.pi / NP)[None, :])
NV = np.stack([np.sin(TH) * np.cos(PH), np.sin(TH) * np.sin(PH), np.cos(TH)], -1)


def project(fn, rho_nodes, lmax):
    """a[l] = sum_m |a_lm(r)|^2 style blocks: returns complex a_lm(r) array [r, l, m]."""
    out = np.zeros((len(rho_nodes), lmax + 1, 2 * lmax + 1), complex)
    Ys = {(l, m): np.conj(sph_harm_y(l, m, TH, PH)) for l in range(lmax + 1) for m in range(-l, l + 1)}
    for a, r in enumerate(rho_nodes):
        P = r * NV
        v = fn(P[..., 0], P[..., 1], P[..., 2]) * np.ones_like(TH)
        for l in range(lmax + 1):
            for m in range(-l, l + 1):
                out[a, l, m + lmax] = np.sum(WS * v * Ys[(l, m)])
    return out


def cl_true(l):
    L = l * (l + 1)
    return 4.0 / (L * (L - 2)) if l >= 2 else 0.0


CLS = {'true': cl_true, 'oneL': lambda l: 1.0 / (l * (l + 1)) if l >= 2 else 0.0,
       'fourL2': lambda l: 4.0 / (l * (l + 1)) ** 2 if l >= 2 else 0.0}


def predict(F, rhos, lmax):
    """returns, per rho: per-l bulk integrals (weighted energy, Born), surface, and the S variants."""
    xr, wr = np.polynomial.legendre.leggauss(NR)
    res = []
    lam = {}
    for par in ('even', 'odd'):
        ps, pd = F['psi'][par], F['psid'][par]
        dr = sp.expand(sum(X[i] * sp.diff(ps, X[i]) for i in range(3)))  # r d_r psi
        lam[par] = [sp.lambdify(X, e, 'numpy') for e in (ps, pd, dr)]
    for rho in rhos:
        rn = rho * 0.5 * (xr + 1); w = rho * 0.5 * wr
        Nn = (rho * rho - rn * rn) / (2 * rho)
        perl = {}
        for par in ('even', 'odd'):
            a, ad, ar = (project(f, rn, lmax) for f in lam[par])
            ap = ar / rn[:, None, None]  # d_r psi projections
            A2 = (np.abs(a) ** 2).sum(-1); Ad2 = (np.abs(ad) ** 2).sum(-1); Ap2 = (np.abs(ap) ** 2).sum(-1)
            Ls = np.array([l * (l + 1) for l in range(lmax + 1)])
            wN = (w * rn * rn * Nn)[:, None]
            energy = (wN * (Ad2 + Ap2 + Ls[None, :] * A2 / (rn * rn)[:, None])).sum(0)
            born = ((w * rn * rn)[:, None] * A2).sum(0) / rho
            noN = ((w * rn * rn)[:, None] * (Ad2 + Ap2 + Ls[None, :] * A2 / (rn * rn)[:, None])).sum(0)
            asurf = project(lam[par][0], np.array([rho]), lmax)
            surf = (np.abs(asurf[0]) ** 2).sum(-1) * rho * rho
            perl[par] = dict(energy=energy.tolist(), born=born.tolist(), surf=surf.tolist(), energy_flat=noN.tolist())
        S = {}
        for name, cf in CLS.items():
            S[name] = 2 * np.pi * sum(cf(l) * (perl[p]['energy'][l] + perl[p]['born'][l]) for p in perl for l in range(lmax + 1))
        sigma = 2 * sum(cl_true(l) * perl[p]['surf'][l] for p in perl for l in range(lmax + 1))
        S['noBorn'] = 2 * np.pi * sum(cl_true(l) * perl[p]['energy'][l] for p in perl for l in range(lmax + 1))
        S['surf_pi2'] = S['true'] + np.pi / 2 * sigma
        lowl = sum(perl[p]['energy'][l] + perl[p]['born'][l] for p in perl for l in (0, 1))
        res.append(dict(rho=rho, S=S, sigma=sigma, lowl=lowl, perl=perl))
    return res


def lambdify_data(Y, Zp):
    """numpy callables for the potentials, components in TCOMP order (sampled on the lattice)."""
    return ([sp.lambdify(X, Y[i, j], 'numpy') for (i, j) in TCOMP], [sp.lambdify(X, Zp[i, j], 'numpy') for (i, j) in TCOMP])


if __name__ == '__main__':
    out = sys.argv[1]
    for f in (P2, P2b, P3, P4):
        assert harmonic_check(f)
    D = data_defs()
    allres = {}
    for name, (Y, Zp) in D.items():
        F = fields(Y, Zp)
        deg = max(sp.Poly(F['psi'][p] + F['psid'][p] + 1, *X).total_degree() for p in ('even', 'odd'))
        lmax = max(deg, 2) + 1
        rhos = [0.9, 0.95, 1.0, 1.05]
        res = predict(F, rhos, lmax)
        chk = {k: str(v) for k, v in F['chk'].items()}
        allres[name] = dict(lmax=lmax, chk=chk, res=res,
                            Y=[str(Y[i, j]) for (i, j) in TCOMP], Z=[str(Zp[i, j]) for (i, j) in TCOMP])
        r1 = res[2]
        S = r1['S']
        lsh = {p: [round(v, 6) for v in np.array(r1['perl'][p]['energy']) + np.array(r1['perl'][p]['born'])] for p in r1['perl']}
        print(f"{name}: S_T1={S['true']:.6f} sigma={r1['sigma']:.6f} | plant ratios S_plant/S_T1: surf(pi/2)={S['surf_pi2']/S['true']:.4f} "
              f"c_l->1/L={S['oneL']/S['true']:.4f} c_l->4/L^2={S['fourL2']/S['true']:.4f} noBorn={S['noBorn']/S['true']:.4f} | l<=1 part={r1['lowl']:.2e}")
        print('   per-l (energy+Born):', lsh, ' checks:', chk)
    json.dump(allres, open(out, 'w'), indent=1, default=float)
