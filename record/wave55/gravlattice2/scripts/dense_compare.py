"""Compare the sector solver with the dense wave52 path on all eleven data (new controls included), R = 5, 6, s and c.
usage: python3 dense_compare.py ../out/dense ../out/rows"""
import json, sys
dd, rd = sys.argv[1:3]
for R in (5, 6):
    for c in 'sc':
        d = json.load(open(f'{dd}/dense_R{R}_{c}.json')); n = json.load(open(f'{rd}/row_R{R}_{c}_B.json'))
        for cap in d['S']:
            print(f'R={R} {c} cap={cap} ' + ' '.join(f"{nm}:{abs(n['S'][cap][k]-d['S'][cap][k])/abs(d['S'][cap][k]):.1e}" for k, nm in enumerate(d['names'])))
