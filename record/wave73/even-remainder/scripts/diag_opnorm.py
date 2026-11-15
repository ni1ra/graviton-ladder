"""Diagnostic for run 1's failed line P5: is lam*||D_1||_H > sqrt(2) a sampling artefact at the bottom of the y-box?
Predictions: ..\\logs\\03_note_on_run1_and_diag_prediction.log.
Usage: ..\\..\\even-nonzero\\venv\\Scripts\\python.exe -B diag_opnorm.py
"""
import sys

sys.dont_write_bytecode = True
import numpy as np

import even_remainder as er

l = 1024
vals = {}
for name, ybot, dy in (("D1", -4.8, 0.04), ("D2", -4.8, 0.02), ("D3", -3.0, 0.04), ("D4", -3.0, 0.02)):
    er.ep.YBOT = ybot
    st = er.setup(l, dy=dy)
    n = st["n"]
    Z = np.zeros((n, n))
    cd = np.sqrt(st["coth"])
    v = st["lam"] * float(np.linalg.norm(cd[:, None] * st["M"] / cd[None, :], 2))      # ||D_1||_H = ||coth^1/2 M coth^-1/2||
    vals[name] = v
    print(f"  {name}: ybot={ybot} dy={dy} n={n}  lam*||D_1||_H = {v:.5f}   x_bot*dy/2 = {np.exp(-ybot) * dy / 2:.3f}", flush=True)
ok = {"D1": abs(vals["D1"] - 2.43) <= 0.3, "D2": abs(vals["D2"] - 1.21) <= 0.25, "D3": 0.90 <= vals["D3"] <= 1.05,
      "D4": abs(vals["D4"] / vals["D3"] - 1) <= 0.03}
for k, v in ok.items():
    print(f"  {k}: {'PASS' if v else 'FAIL'}")
