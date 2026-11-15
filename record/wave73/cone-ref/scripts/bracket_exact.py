"""Exact (sympy closed form) version of the bracket check of probe.py, p = 1, 2, 3."""
import sympy as sp
v, u = sp.symbols("v u", positive=True)
def Pi(j, p):
    return sp.expand(v**(1-p) * sp.integrate(sp.expand((u**(p-1) * sp.diff(j, v)).subs(v, u)), (u, 0, v)))
def W(j1, j2):
    P1 = sp.integrate(sp.expand(j1.subs(v, u)), (u, 0, v))
    Pinf = sp.integrate(sp.expand(j1), (v, 0, sp.oo))
    return sp.nsimplify(sp.simplify(sp.integrate(sp.expand((2*P1 - Pinf) * j2), (v, 0, sp.oo))))
js = [(v - v**2/2)*sp.exp(-v), (v**2 - v**3/3)*sp.exp(-v), (v**3 - v**4/4)*sp.exp(-v)]
lines = []
for (i, k) in ((0, 1), (0, 2), (1, 2)):
    for p in (1, 2, 3):
        try:
            r = W(Pi(js[i], p), Pi(js[k], p))
        except Exception as e:
            r = "sympy failed: " + type(e).__name__
        lines.append(f"pair{(i,k)} p={p}: W = {r}"); print(lines[-1], flush=True)
open("logs/bracket_exact.log", "w", encoding="utf-8").write("\n".join(lines) + "\n")
