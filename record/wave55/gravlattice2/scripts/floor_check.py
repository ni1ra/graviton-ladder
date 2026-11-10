"""Reproducibility floor of the wave52 reference: the original code (wave52 run.py, unchanged) rerun single-threaded at R = 5, 6 (site),
compared with its shipped rows and with the sector solver.  usage: python3 floor_check.py ../out/floor ../../../wave52/gravlattice/out/rows ../out/rows"""
import json, sys
fl, old, new = sys.argv[1:4]
for R in (5, 6):
    o = json.load(open(f'{old}/row_R{R}_s_B.json')); r = json.load(open(f'{fl}/orig_1thread_row_R{R}_s_B.json')); n = json.load(open(f'{new}/row_R{R}_s_B.json'))
    for cap in o['S']:
        for k, name in enumerate(o['names']):
            a, b, c = o['S'][cap][k], r['S'][cap][k], n['S'][cap][n['names'].index(name)]
            print(f"R={R} cap={cap} {name:7s} rerun-vs-shipped {abs(b-a)/abs(a):.2e}  sector-vs-shipped {abs(c-a)/abs(a):.2e}  sector-vs-rerun {abs(c-b)/abs(b):.2e}")
