import sys
from fractions import Fraction
import hb
from flint import arb
n, l, prec = 24, 8, 1400
hb.setprec(prec)
XA, PA, dg = hb.vacuum_blocks(n, 4 * n, l, "D2", prec)
mo = hb.modular(XA, PA)
nu, eps = mo["nu"], mo["eps"]
A1 = hb.mid(hb.scale_cols(mo["W"], [1 / x.sqrt() for x in nu]))
A1i = hb.mid(hb.diagmat([x.sqrt() for x in nu]) * mo["Winv"])
a = Fraction(2, 2 * n + 1)
c2 = hb.second_order(eps, A1, a, [1] * n)
print("c2 (second order) =", c2.str(20), " Q2 =", (c2 * l * (l + 1)).str(15))
for den in (1, 4, 16, 64, 256):
    vals = []
    for sgn in (1, -1):
        g = Fraction(sgn, den)
        Mp = hb.perturbed_M(mo["M"], a, [1] * n, g)
        E2, A2, A2i, _ = hb.williamson(Mp, mo["N"])
        vals.append(hb.one_minus_A(hb.affinity_I(eps, A1, A1i, E2, A2, A2i)) * den * den)
    m = (vals[0] + vals[1]) / 2
    print("g=+-1/%d: exact mean (1-A)/g^2 = %s   minus c2: %s   odd part/g: %s" % (den, m.str(18), (m - c2).str(5), ((vals[1] - vals[0]) / 2 * den).str(8)))
