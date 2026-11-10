"""Dense cross-check of the new controls: wave52's dense path (glat.corr_matrix, glat.commutator, glat.rel_entropy, unchanged)
on the eleven data of pred2, rule B.  usage: python3 run_dense_ctrl.py R s|c outdir"""
import numpy as np, json, sys, os
from glat import *
import pred, pred2

R, ctag, outdir = int(sys.argv[1]), sys.argv[2], sys.argv[3]
HERE = os.path.dirname(os.path.abspath(__file__))
centre2 = np.array([0, 0, 0]) if ctag == 's' else np.array([1, 1, 1])
centre = centre2 / 2.0
Te, Pe = cells(R, centre2, dual=False)
Tb, Pb = cells(R, centre2, dual=True)
inset2 = {(int(t), tuple(int(v) for v in p)) for t, p in zip(Te, Pe)}
ke = stencil_inputs_inside(Te, Pe, INC, inset2)
kb = stencil_inputs_inside(Tb, Pb, SC, inset2)
Te, Pe, Tb, Pb = Te[ke], Pe[ke], Tb[kb], Pb[kb]
Dm = int(max(np.abs(Pe - Pe.min(0)).max(), np.abs(Pb - Pb.min(0)).max()) // 2 + 4)
Garr, tabs = kernel_tables(Dm)
X = corr_matrix(Te, Pe, Garr, tabs)
P = corr_matrix(Tb, Pb, Garr, tabs)
W = commutator(Te, Pe, Tb, Pb)
D = pred2.data_defs(os.path.join(HERE, '..', 'out', 'controls.json'))
names = list(D)
DE = np.zeros((len(Te), len(names))); DB = np.zeros((len(Tb), len(names)))
for n_, name in enumerate(names):
    fy, fz = pred.lambdify_data(*D[name])

    def Yfn(t, p2):
        u = (p2 / 2.0 - centre) / R
        return fy[t](u[:, 0], u[:, 1], u[:, 2]) * np.ones(len(u))

    def Zfn(t, p2):
        u = (p2 / 2.0 - centre) / R
        return fz[t](u[:, 0], u[:, 1], u[:, 2]) * np.ones(len(u))

    def hfn(t, p2):
        return apply_stencil_fn(np.full(len(p2), t), p2, SC, Yfn)

    def pfn(t, p2):
        return apply_stencil_fn(np.full(len(p2), t), p2, INC, Zfn)
    DE[:, n_] = apply_stencil_fn(Te, Pe, INC, hfn)
    DB[:, n_] = apply_stencil_fn(Tb, Pb, SC, pfn)
S, info = rel_entropy(X, P, W, DE, DB)
json.dump(dict(R=R, centre=ctag, names=names, S={str(k): v.tolist() for k, v in S.items()}),
          open(f'{outdir}/dense_R{R}_{ctag}.json', 'w'), indent=1)
print('done', R, ctag)
