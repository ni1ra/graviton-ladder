"""Post hoc check of the far-region density used to motivate a (ln l)/l term in the large-l fit.

Deep inside the ball the modular Hamiltonian is cold, both states are locally ground states of
fields with mass^2 m^2 = l(l+1)/r^2 and m^2 + g/N(r), N = (1-r^2)/2, and second-order perturbation
theory of the ground-state overlap gives the local density

    d(1-A)/dr = g^2 / (128 N(r)^2 m(r)^3)                       (valid for m N >> 1).

Test: the second-order coefficient with the perturbation restricted to r < rc, on the lattice,
against  l(l+1) Int_0^rc dr r^3 / (128 N^2 (l(l+1))^{3/2}).
Expanding the density near the surface, rho = 1 - r:  rho^-2 (1 - 2 rho + ...)/(128 L^{3/2}),
whose 1/rho term integrates to a contribution  -(ln l)/(64 l)  to Q_l.

  python far_region.py <n> <l>
"""
import sys
from fractions import Fraction
import hb
from flint import arb
from run import default_prec
from mpmath import mp, quad, mpf

n, l = int(sys.argv[1]), int(sys.argv[2])
prec = default_prec(n, l)
hb.setprec(prec)
XA, PA, dg = hb.vacuum_blocks(n, 4 * n, l, "D2", prec)
mo = hb.modular(XA, PA)
nu, eps = mo["nu"], mo["eps"]
A1 = hb.mid(hb.scale_cols(mo["W"], [1 / x.sqrt() for x in nu]))
a = Fraction(2, 2 * n + 1)
L = l * (l + 1)
full = float(hb.second_order(eps, A1, a, [1] * n) * L)
print("n=%d l=%d  Q2 (full ball) = %.8f" % (n, l, full))
print("   rc     sites   Q2[r<rc] lattice   local-density integral   ratio    m*N at rc")
for rc in (0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9):
    jc = int(rc * (n + 0.5))          # sites 1..jc, cell edge at (jc+1/2) a
    redge = (jc + 0.5) / (n + 0.5)
    w = [1] * jc + [0] * (n - jc)
    q = float(hb.second_order(eps, A1, a, w) * L)
    pred = float(quad(lambda r: r ** 3 / (128 * ((1 - r * r) / 2) ** 2), [0, redge]) / mpf(L) ** 0.5)
    mN = (L ** 0.5 / redge) * (1 - redge ** 2) / 2
    print("  %.3f  %4d    %.6e       %.6e          %.4f    %.2f" % (redge, jc, q, pred, q / pred, mN))
