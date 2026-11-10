"""w55-gravlattice2 analysis: the PREREG §5 estimator, the §6 verdict, the §7 plants.  Rule B only.
usage: python3 fit2.py ../out/rows ../out/pred.json ../out/fit.json"""
import numpy as np, json, glob, sys

rows_dir, predf, outf = sys.argv[1:4]
pred = json.load(open(predf))
rows = [json.load(open(f)) for f in sorted(glob.glob(rows_dir + '/row_*_B.json'))]
names = rows[0]['names']
CAP = '30.0'
DATA = ['g1_E2', 'g2_B2', 'g3_E3', 'g4_E23', 'g5_E4', 'g6_EB']
CTRL_P = {'g1_E2': 'd1_E2', 'g3_E3': 'd1_E2', 'g4_E23': 'd1_E2', 'g5_E4': 'd1_E2', 'g2_B2': 'd2_B2', 'g6_EB': 'mixEB'}
CTRL_A = {'g1_E2': 'd3_E3', 'g3_E3': 'd3_E3', 'g4_E23': 'd3_E3', 'g5_E4': 'd3_E3', 'g2_B2': 'c2_B', 'g6_EB': 'c1_E'}


def interp(name, key, rho, var=None):
    rs = [r['rho'] for r in pred[name]['res']]
    vs = [(r['S'][var] if var else r[key]) for r in pred[name]['res']]
    return float(np.polyval(np.polyfit(rs, vs, 3), rho))


def S_T1(name, rho, var='true'):
    return interp(name, None, rho, var)


def sigma(name, rho):
    return interp(name, 'sigma', rho)


tab = {}
for r in rows:
    for n, name in enumerate(r['names']):
        caps = [r['S'][c][n] for c in ('25.0', '30.0', '35.0')]
        tab.setdefault(name, {}).setdefault(r['R'], {})[r['centre']] = dict(
            S=r['S'][CAP][n], capspread=(max(caps) - min(caps)) / max(abs(caps[1]), 1e-300), rho=r['R_eff'] / r['R'],
            KL=r['info']['KLq'][n] + r['info']['KLp'][n],
            null=max(r['info']['dq_null_rel'][n], r['info']['dp_null_rel'][n]), nu_min=r['info']['nu_min'])
Rs = sorted(R for R in tab['g1_E2'] if len(tab['g1_E2'][R]) == 2)
Rmax = max(Rs)
wE = S_T1('g1_E2', 1.0); wB = S_T1('g2_B2', 1.0)


def ratio(name, conv='vol', centres=('s', 'c'), var='true', drop_KL=False, inject=0.0, Rlist=None):
    """mean over centrings of S_lat / S_var(rho).  'mixEB' = S_T1-weighted mean of the d1 and d2 control ratios."""
    if name == 'mixEB':
        a = ratio('d1_E2', conv, centres, var, drop_KL, 0.0, Rlist); b = ratio('d2_B2', conv, centres, var, drop_KL, 0.0, Rlist)
        return (wE * a + wB * b) / (wE + wB)
    out = []
    for R in (Rlist or Rs):
        v = []
        for c in centres:
            d = tab[name][R][c]
            rho = 1.0 if conv == 'raw' else d['rho']
            S = d['S'] - (d['KL'] if drop_KL else 0.0) + inject * sigma(name, rho)
            v.append(S / S_T1(name, rho, var))
        out.append(np.mean(v))
    return np.array(out)


def fitwin(ys, lo, deg):
    Rw = np.array([R for R in Rs if lo <= R <= Rmax], float)
    yw = np.array([y for R, y in zip(Rs, ys) if lo <= R <= Rmax])
    A = np.vstack([Rw ** (-k) for k in range(deg + 1)]).T
    return float(np.linalg.lstsq(A, yw, rcond=None)[0][0])


def estimator(ys_vol, ys_raw):
    """PREREG §5: primary = lin on [Rmax-4, Rmax] (vol); bracket = 6 estimates x 2 conventions."""
    prim = fitwin(ys_vol, Rmax - 4, 1)
    br = []
    for ys in (ys_vol, ys_raw):
        br += [fitwin(ys, Rmax - 3, 1), fitwin(ys, Rmax - 4, 1), fitwin(ys, Rmax - 5, 1),
               fitwin(ys, Rmax - 5, 2), fitwin(ys, Rmax - 6, 2), float(ys[-1])]
    return dict(primary=prim, lo=min(br), hi=max(br), all=br)


log = [f'# rule B, R = {Rs}, Rmax = {Rmax}, cap {CAP}; ratio = mean over site and cube of S_lat/S_T1 (vol: S_T1 at R_eff)']
out = dict(R=Rs, Rmax=Rmax, ratios={}, c={}, plants={}, cl={})
for name in names:
    v, rw = ratio(name), ratio(name, 'raw')
    e = estimator(v, rw)
    cs = max(tab[name][R][c]['capspread'] for R in Rs for c in 'sc')
    kl = [float(np.mean([tab[name][R][c]['KL'] / tab[name][R][c]['S'] for c in 'sc'])) for R in Rs]
    out['ratios'][name] = dict(vol=v.tolist(), raw=rw.tolist(), est=e, capspread=cs, KLshare=kl)
    log.append(f"{name:7s} vol {[round(float(x), 3) for x in v]} | lin[{Rmax-4},{Rmax}] {e['primary']:.3f} [{e['lo']:.3f}, {e['hi']:.3f}]"
               f" | capspread {cs:.1e} | KL share R={Rmax}: {kl[-1]:.2f}")

# control convergence (PREREG §6 (i))
ctrl_ok = {}
for cn in ('d1_E2', 'd2_B2'):
    e = out['ratios'][cn]['est']; last = out['ratios'][cn]['vol'][-1]
    ctrl_ok[cn] = bool(abs(e['primary'] - 1) <= 0.10 and abs(last - 1) <= 0.15)
    log.append(f"control {cn}: a_ctrl {e['primary']:.3f}, ratio(Rmax) {last:.3f} -> {'converged' if ctrl_ok[cn] else 'NOT converged'}")
out['ctrl_ok'] = ctrl_ok


def cdiff(name, ctrl, centres=('s', 'c'), drop_KL=False, inject=0.0):
    so = sigma(name, 1.0) / S_T1(name, 1.0)
    cv = (ratio(name, 'vol', centres, drop_KL=drop_KL, inject=inject) - ratio(ctrl, 'vol', centres, drop_KL=drop_KL)) / so
    cr = (ratio(name, 'raw', centres, drop_KL=drop_KL, inject=inject) - ratio(ctrl, 'raw', centres, drop_KL=drop_KL)) / so
    return estimator(cv, cr), cv, so


log.append('\n# surface coefficient c(R) = (ratio_datum - ratio_ctrl)/(sigma/S_T1); c_inf = lin[Rmax-4,Rmax]; bracket 12 estimates')
for name in DATA:
    res = {}
    for tag, ctrl, cen in (('primary', CTRL_P[name], ('s', 'c')), ('alternate', CTRL_A[name], ('s', 'c')),
                           ('site', CTRL_P[name], ('s',)), ('cube', CTRL_P[name], ('c',))):
        e, cv, so = cdiff(name, ctrl, cen)
        res[tag] = dict(ctrl=ctrl, c_R=cv.tolist(), **e)
        log.append(f"{name:7s} {tag:9s} ({ctrl}): c(R) {[round(float(x), 3) for x in cv]} -> c_inf {e['primary']:+.3f} [{e['lo']:+.3f}, {e['hi']:+.3f}]")
    res['sigma_over_S'] = so
    out['c'][name] = res

# verdict (PREREG §6)
P = {n: out['c'][n]['primary'] for n in DATA}
thr_rej = {n: bool(not (P[n]['lo'] <= 0.15 <= P[n]['hi']) and not (P[n]['lo'] <= -0.15 <= P[n]['hi'])) for n in DATA}
pass_ = all(ctrl_ok.values()) and all(abs(P[n]['primary']) < 0.15 for n in DATA) and sum(thr_rej.values()) >= 4


def big(n, tag='primary'):
    e = out['c'][n][tag]
    return abs(e['primary']) >= 0.15 and (e['lo'] > 0 or e['hi'] < 0)


Etype = ['g1_E2', 'g3_E3', 'g4_E23', 'g5_E4']
kill = False
if all(ctrl_ok.values()) and big('g2_B2'):
    sg = np.sign(P['g2_B2']['primary'])
    for n in Etype:
        if big(n) and np.sign(P[n]['primary']) == sg and all(big(m, t) and np.sign(out['c'][m][t]['primary']) == sg
                                                               for m in (n, 'g2_B2') for t in ('site', 'cube', 'alternate')):
            kill = True
verdict = 'PASS' if pass_ else ('KILL' if kill else 'INCONCLUSIVE')
out['verdict'] = dict(verdict=verdict, threshold_rejected=thr_rej, n_threshold_rejected=int(sum(thr_rej.values())),
                      all_c_below=bool(all(abs(P[n]['primary']) < 0.15 for n in DATA)))
log.append(f"\n# VERDICT (PREREG §6): {verdict} | controls converged {ctrl_ok} | all |c_inf| < 0.15: {out['verdict']['all_c_below']}"
           f" | threshold plants c = +-0.15 rejected on {sum(thr_rej.values())}/6: {thr_rej}")

# plants (PREREG §7): a theory plant predicts Delta = S_plant/S_T1(datum) - S_plant/S_T1(ctrl); caught if outside the Delta bracket
log.append('\n# plants: plant value of c (= Delta_plant/(sigma/S)) against the primary c bracket; caught if outside')
for name in DATA:
    e = P[name]; so = out['c'][name]['sigma_over_S']; ctrl = CTRL_P[name]

    def plr(nm, var):
        if nm == 'mixEB':
            return (wE * plr('d1_E2', var) + wB * plr('d2_B2', var)) / (wE + wB)
        return S_T1(nm, 1.0, var) / S_T1(nm, 1.0)
    pl = {}
    for var in ('surf_pi2', 'oneL', 'fourL2', 'noBorn'):
        cp = (plr(name, var) - plr(ctrl, var)) / so
        pl[var] = dict(c_plant=cp, caught=bool(not (e['lo'] <= cp <= e['hi'])), power=bool(abs(cp) > 1e-9))
    for cp in (0.15, -0.15):
        pl[f'thr{cp:+.2f}'] = dict(c_plant=cp, caught=bool(not (e['lo'] <= cp <= e['hi'])), power=True)
    # lattice plant 6: drop the classical centre (KL) from S_lat (datum and control)
    e6, cv6, _ = cdiff(name, ctrl, drop_KL=True)
    pl['noCentre'] = dict(c_inf=e6['primary'], lo=e6['lo'], hi=e6['hi'], caught=bool(e6['lo'] > 0.15 or e6['hi'] < -0.15),
                          power=bool(max(out['ratios'][name]['KLshare']) > 1e-6))
    # lattice plant 7: injection
    inj = {}
    for ci in (0.15, np.pi / 2):
        e7, _, _ = cdiff(name, ctrl, inject=ci)
        inj[f'{ci:.3f}'] = e7['primary'] - e['primary']
    pl['inject'] = inj
    out['plants'][name] = pl
    log.append(f"{name:7s} " + ' '.join(f"{k}: c={v['c_plant']:+.2f} {'CAUGHT' if v['caught'] else ('no power' if not v['power'] else 'missed')}"
                                       for k, v in pl.items() if 'c_plant' in v)
               + f" | noCentre c_inf {pl['noCentre']['c_inf']:+.3f} [{pl['noCentre']['lo']:+.3f}, {pl['noCentre']['hi']:+.3f}] "
               + ('CAUGHT' if pl['noCentre']['caught'] else ('no power' if not pl['noCentre']['power'] else 'missed'))
               + f" | inject recovered {inj}")

# c_l test (secondary): D(R) = ratio(datum)/ratio(g1); estimator on D; plants
log.append('\n# c_l test: D(R) = ratio(datum)/ratio(g1_E2); estimator; plant values')
for name in ('g3_E3', 'g4_E23', 'g5_E4', 'g2_B2'):
    Dv = ratio(name) / ratio('g1_E2'); Dr = ratio(name, 'raw') / ratio('g1_E2', 'raw')
    e = estimator(Dv, Dr)
    plants = {v: (S_T1(name, 1.0, v) / S_T1(name, 1.0)) / (S_T1('g1_E2', 1.0, v) / S_T1('g1_E2', 1.0)) for v in ('oneL', 'fourL2', 'noBorn')}
    caught = {v: bool(not (e['lo'] <= p <= e['hi'])) for v, p in plants.items()}
    out['cl'][name] = dict(D=Dv.tolist(), est=e, plants=plants, caught=caught)
    log.append(f"{name:7s} D(R) {[round(float(x), 3) for x in Dv]} -> {e['primary']:.3f} [{e['lo']:.3f}, {e['hi']:.3f}] | "
               + ' '.join(f"{v}={p:.3f} {'CAUGHT' if caught[v] else 'missed'}" for v, p in plants.items()))
out['max_null'] = max(tab[n][R][c]['null'] for n in names for R in Rs for c in 'sc')
out['nu_min'] = min(tab[n][R][c]['nu_min'] for n in names for R in Rs for c in 'sc')
log.append(f"\nmax null shift {out['max_null']:.1e}; nu_min {out['nu_min']:.12f}")
json.dump(out, open(outf, 'w'), indent=1)
open(outf.replace('.json', '.log'), 'w').write('\n'.join(log) + '\n')
print('\n'.join(log))
