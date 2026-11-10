"""Gate: the sector solver must reproduce every wave52 rule-B row (R = 3-7, s and c) for every datum and cap,
relative 1e-8 (absolute 1e-10 where S < 1e-2).  usage: python3 verify.py ../out/rows ../../../wave52/gravlattice/out/rows ../out/verify.json"""
import json, sys, glob, os
new_dir, old_dir, outf = sys.argv[1:4]
res, worst, ok = [], 0.0, True
for R in range(3, 8):
    for c in 's', 'c':
        fo, fn = f'{old_dir}/row_R{R}_{c}_B.json', f'{new_dir}/row_R{R}_{c}_B.json'
        if not os.path.exists(fn):
            ok = False; res.append(dict(R=R, c=c, missing=True)); continue
        o, n = json.load(open(fo)), json.load(open(fn))
        assert (o['n_e'], o['n_b'], o['nsites']) == (n['n_e'], n['n_b'], n['nsites'])
        for cap in o['S']:
            for k, name in enumerate(o['names']):
                so, sn = o['S'][cap][k], n['S'][cap][n['names'].index(name)]
                d = abs(sn - so); rel = d / abs(so)
                passed = (rel < 1e-8) or (abs(so) < 1e-2 and d < 1e-10)
                ok &= passed
                worst = max(worst, rel)
                res.append(dict(R=R, c=c, cap=cap, name=name, old=so, new=sn, rel=rel, passed=passed))
        for key in ('n_q', 'n_p', 'rank', 'cq', 'cp'):
            ok &= o['info'][key] == n['info'][key]
            if o['info'][key] != n['info'][key]:
                print('MISMATCH', R, c, key, o['info'][key], n['info'][key])
summary = dict(passed=bool(ok), worst_rel=worst, n_compared=len(res))
json.dump(dict(summary=summary, rows=res), open(outf, 'w'), indent=1)
print(summary)
for r in sorted([r for r in res if 'rel' in r], key=lambda r: -r['rel'])[:5]:
    print(r)
