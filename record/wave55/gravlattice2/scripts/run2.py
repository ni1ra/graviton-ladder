"""w55-gravlattice2 production: one rule-B row per (R, centre) with all eleven data, via the Z2^3 sector solver.
usage: python3 run2.py R s|c outdir [Garr.npy]"""
import numpy as np, json, sys, time, os
from glat import cells, stencil_inputs_inside, kernel_tables, apply_stencil_fn, INC, SC
from maxlat_core import G_array
import glat2, pred, pred2

R, ctag, outdir = int(sys.argv[1]), sys.argv[2], sys.argv[3]
gcache = sys.argv[4] if len(sys.argv) > 4 else None
HERE = os.path.dirname(os.path.abspath(__file__))
t0 = time.time()
centre2 = np.array([0, 0, 0]) if ctag == 's' else np.array([1, 1, 1])
centre = centre2 / 2.0
Te, Pe = cells(R, centre2, dual=False)
Tb, Pb = cells(R, centre2, dual=True)
nsites = int((Te < 3).sum() / 3)
inset = {(int(t),) + tuple(int(v) for v in p) for t, p in zip(Te, Pe)}
inset2 = {(k[0], (k[1], k[2], k[3])) for k in inset}
ke = stencil_inputs_inside(Te, Pe, INC, inset2)
kb = stencil_inputs_inside(Tb, Pb, SC, inset2)
Te, Pe, Tb, Pb = Te[ke], Pe[ke], Tb[kb], Pb[kb]
row = dict(R=R, centre=ctag, rule='B', solver='Z2^3 sectors', n_e=int(len(Te)), n_b=int(len(Tb)), nsites=nsites,
           R_eff=float((3 * nsites / (4 * np.pi)) ** (1 / 3)))
# sign consistency: the b correlator uses the inc kernel on the dual lattice; its reflection signs must differ from sc's by a t-independent factor
for s in glat2.GROUP:
    f = {glat2.stencil_sign(SC, t, s) * glat2.stencil_sign(INC, t, s) for t in range(6)}
    assert len(f) == 1, ('sign mismatch', s, f)
Dm = int(max(np.abs(Pe - Pe.min(0)).max(), np.abs(Pb - Pb.min(0)).max()) // 2 + 4)
if gcache and os.path.exists(gcache):
    G = np.load(gcache)
    assert G.shape[0] >= Dm + 5, 'G cache too small'
    import glat
    glat.G_array = lambda D: G[:D + 1, :D + 1, :D + 1]
Garr, tabs = kernel_tables(Dm)
perm_e, sig_e, rep_e = glat2.orbits(Te, Pe, centre2, INC)
perm_b, sig_b, rep_b = glat2.orbits(Tb, Pb, centre2, SC)
sec_e, reps_e = glat2.sector_bases(perm_e, sig_e, rep_e)
sec_b, reps_b = glat2.sector_bases(perm_b, sig_b, rep_b)
Xs = glat2.sector_matrix(Te, Pe, reps_e, sec_e, Garr, tabs)
Ps = glat2.sector_matrix(Tb, Pb, reps_b, sec_b, Garr, tabs)
W = glat2.commutator_sparse(Te, Pe, Tb, Pb)
Ws = [(Ve.T @ W @ Vb).toarray() for (_, _, Ve), (_, _, Vb) in zip(sec_e, sec_b)]
# off-sector leakage of W (must vanish): largest |V_k^T W V_k'| for k != k'
leak = max(abs((sec_e[i][2].T @ W @ sec_b[j][2])).max() if sec_e[i][2].shape[1] and sec_b[j][2].shape[1] else 0.0
           for i in range(8) for j in range(8) if i != j)
row['W_offsector_max'] = float(leak)
row['t_build'] = time.time() - t0
D = pred2.data_defs(os.path.join(HERE, '..', 'out', 'controls.json'))
names = list(D)
DE = np.zeros((len(Te), len(names))); DB = np.zeros((len(Tb), len(names)))
for n_, name in enumerate(names):
    Y, Zp = D[name]
    fy, fz = pred.lambdify_data(Y, Zp)

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
dqs = [Ve.T @ DE for (_, _, Ve) in sec_e]
dps = [Vb.T @ DB for (_, _, Vb) in sec_b]
# the sector bases are orthonormal and complete: data norm is preserved
row['data_norm_check'] = float(max(abs(sum((d ** 2).sum(0) for d in dqs) - (DE ** 2).sum(0)).max() / max((DE ** 2).sum(), 1e-300),
                                   abs(sum((d ** 2).sum(0) for d in dps) - (DB ** 2).sum(0)).max() / max((DB ** 2).sum(), 1e-300)))
S, info = glat2.rel_entropy_sectors(Xs, Ps, Ws, dqs, dps)
nq = np.abs(DE).max(0) + 1e-300; npp = np.abs(DB).max(0) + 1e-300
info['dq_null_rel'] = (info.pop('dq_null_abs') / nq * (np.abs(DE).max(0) > 0)).tolist()
info['dp_null_rel'] = (info.pop('dp_null_abs') / npp * (np.abs(DB).max(0) > 0)).tolist()
info['KLq'] = info['KLq'].tolist(); info['KLp'] = info['KLp'].tolist()
row['names'] = names
row['S'] = {str(k): v.tolist() for k, v in S.items()}
row['info'] = info
for n_, name in enumerate(names):
    print(name, R, ctag, 'B', S[30.0][n_], flush=True)
print({k: v for k, v in info.items() if k not in ('dq_null_rel', 'dp_null_rel', 'KLq', 'KLp', 'sectors')},
      'max null', max(info['dq_null_rel'] + info['dp_null_rel']), 'W leak', leak, 'norm', row['data_norm_check'])
print('sectors', [(d['n_e'], d['n_b']) for d in info['sectors']])
row['t_total'] = time.time() - t0
json.dump(row, open(f'{outdir}/row_R{R}_{ctag}_B.json', 'w'), indent=1, default=float)
print('done', R, ctag, 'B', row['t_total'], 'n_e', row['n_e'], 'n_b', row['n_b'])
