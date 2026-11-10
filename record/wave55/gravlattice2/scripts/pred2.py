"""w55-gravlattice2 theory side: the first run's eight data (pred.data_defs, unchanged) plus three lower-degree controls
(degree-5/6 potentials from ctrl_solve.py, psi = (1-r^2) P, sigma = 0) and the predictions S_T1 and plants via pred.predict.
usage: python3 pred2.py ../out/controls.json ../out/pred.json"""
import sympy as sp, numpy as np, json, sys
import pred
from pred import X, TCOMP


def data_defs(ctrl_file):
    D = pred.data_defs()
    C = json.load(open(ctrl_file))
    zero = sp.zeros(3, 3)
    for name, comps in C.items():
        T = sp.zeros(3, 3)
        for (i, j), s in zip(TCOMP, comps):
            T[i, j] = T[j, i] = sp.sympify(s, locals={'x': X[0], 'y': X[1], 'z': X[2]})
        D[name] = (T, zero) if '_E' in name else (zero, T)
    return D


if __name__ == '__main__':
    cf, out = sys.argv[1:3]
    D = data_defs(cf)
    allres = {}
    for name, (Y, Zp) in D.items():
        F = pred.fields(Y, Zp)
        deg = max(sp.Poly(F['psi'][p] + F['psid'][p] + 1, *X).total_degree() for p in ('even', 'odd'))
        lmax = max(deg, 2) + 1
        res = pred.predict(F, [0.9, 0.95, 1.0, 1.05], lmax)
        ydeg = max(sp.Poly(Y[i, j] + Zp[i, j] + 1, *X).total_degree() for (i, j) in TCOMP)
        allres[name] = dict(lmax=lmax, potential_degree=ydeg, chk={k: str(v) for k, v in F['chk'].items()}, res=res)
        S = res[2]['S']
        print(f"{name}: potential degree {ydeg} S_T1={S['true']:.6f} sigma={res[2]['sigma']:.3e} | surf(pi/2)={S['surf_pi2']/S['true']:.4f} "
              f"1/L={S['oneL']/S['true']:.4f} 4/L^2={S['fourL2']/S['true']:.4f} noBorn={S['noBorn']/S['true']:.4f} "
              f"c=0.15 plant={(S['true'] + 0.15 * res[2]['sigma'])/S['true']:.4f} | l<=1 {res[2]['lowl']:.1e} | chk {allres[name]['chk']}", flush=True)
    json.dump(allres, open(out, 'w'), indent=1, default=float)
