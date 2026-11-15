"""HINGE-BLIND-2: independent numerical evaluation of the planar (Rindler) second-order limit.

Inputs: (i) the second-order formula of hb.second_order (checked against exact affinities on the
lattice, logs/pt2_check_n24_l8.log); (ii) Rindler modes K_{i w}(k rho), eps = 2 pi w, for
K = 2 pi Int rho T_00 with transverse momentum k, and Z = (g/2) Int phi^2 d rho.
Mode overlap G(w,w') = sqrt(sinh(pi w) sinh(pi w')) / (2k cosh(pi(w+w')/2) cosh(pi(w-w')/2)).
k^2 c2 = (pi^2/2) Int Int dw dw' k^2 G^2 [sinh^2((e+e')/4)/(e+e')^2 + sinh^2((e-e')/4)/(e-e')^2]
                                  / (sinh(e/2) sinh(e'/2)),   e = 2 pi w.
Brute-force 2-d quadrature, no analytic simplification, compared with 7 zeta(3)/(16 pi^2).
"""
import numpy as np
from scipy import integrate
from scipy.special import zeta
pi = np.pi
def f(w, v):
    s, d = w + v, w - v
    t1 = np.sinh(pi * s / 2) ** 2 / (2 * pi * s) ** 2
    t2 = 1.0 / 16 if abs(d) < 1e-9 else np.sinh(pi * d / 2) ** 2 / (2 * pi * d) ** 2
    g2 = 1.0 / (4 * np.cosh(pi * s / 2) ** 2 * np.cosh(pi * d / 2) ** 2)   # k^2 G^2 / (sinh sinh)
    return (pi ** 2 / 2) * g2 * (t1 + t2)
import sys
W = float(sys.argv[1]) if len(sys.argv) > 1 else 40.0
val, err = integrate.dblquad(f, 0, W, 0, W, epsabs=1e-11, epsrel=1e-11)
# tail beyond W: only the |d| ~ O(1) strip survives; integrand -> (pi^2/2)(1/(4 pi^2 s^2)) / cosh^2(pi d/2)
tail = 1 / (16 * pi * 2 * W)   # Int_{s>2W} ds/(32 s^2) * (1/2) Int dd sech^2(pi d/2), leading order in 1/W
print("2-d quadrature on [0,%g]^2 : %.9f (quad err %.1e)" % (W, val, err))
print("tail estimate              : %.9f" % tail)
print("sum                        : %.9f" % (val + tail))
print("7 zeta(3) / (16 pi^2)      : %.9f" % (7 * zeta(3) / (16 * pi ** 2)))
