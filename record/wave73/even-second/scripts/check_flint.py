"""EVEN-SECOND self-check: every integral recomputed by rigorous ball quadrature
(Arb acb.integral) straight from the definitions (no series, no closed forms),
then compared with out/values.json from second.py. Also the plant test.
Run with the python-flint interpreter.
"""
import json, os, sys
from flint import acb, arb, ctx

ctx.prec = 240
HERE = os.path.dirname(os.path.abspath(__file__))
vals = json.load(open(os.path.join(HERE, "..", "out", "values.json")))
LAM2 = 272
lam = arb(LAM2).sqrt()

def t(y): return (y - 2) / 4
def u(y): return (1 - t(y) ** 2) ** 4
def du(y): return -2 * t(y) * (1 - t(y) ** 2) ** 3
def ddu(y): return (-2 * (1 - t(y) ** 2) ** 3 + 12 * t(y) ** 2 * (1 - t(y) ** 2) ** 2) / 4
def Vinf(y): return (-2 * y).exp()
def V16(y): return (-2 * y).exp() * (1 - (-2 * y).exp() / (4 * LAM2)) ** (-2)

def integ(f):
    r = acb.integral(lambda y, _: f(y), -2, 6, rel_tol=arb(10) ** -45, abs_tol=arb(10) ** -45,
                     eval_limit=10 ** 7, depth_limit=10 ** 4)
    return r.real

Q = {
    "norm_u_sq": integ(lambda y: u(y) ** 2),
    "A_inf": integ(lambda y: du(y) ** 2 + Vinf(y) * u(y) ** 2),
    "A_16": integ(lambda y: du(y) ** 2 + V16(y) * u(y) ** 2),
    "G_inf_sq": integ(lambda y: (-ddu(y) + Vinf(y) * u(y)) ** 2),
    "D_sq": integ(lambda y: ((V16(y) - Vinf(y)) * u(y)) ** 2),
    "Q_16": integ(lambda y: y.exp() * (1 + (-y).exp() / (2 * lam)) ** 2
                  * (du(y) ** 2 + V16(y) * u(y) ** 2)),
}

def ball(name):
    lo, hi = arb(vals[name]["lo"]), arb(vals[name]["hi"])
    return lo.union(hi)

ok = True
print(f"{'integral':10s} {'quadrature ball (Arb)':52s} overlap  |mid diff|")
for name, q in Q.items():
    p = ball(name)
    ov = q.overlaps(p)
    ok &= ov
    print(f"{name:10s} {q.str(30):52s} {str(ov):8s} {abs(q.mid() - p.mid()).str(3, radius=False)}")

p = ball("PLANT_A_16")
caught = not Q["A_16"].overlaps(p)
print(f"\nPLANT (bracket exponent -1): closed-form value {p.str(15)}")
print(f"  quadrature A_16 {Q['A_16'].str(15)}; difference {(p - Q['A_16']).str(5)}; "
      f"caught (no overlap): {caught}")
print(f"\nall integrals agree: {ok}; plant caught: {caught}")
sys.exit(0 if (ok and caught) else 1)
