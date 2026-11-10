"""Electric area side: closed form of the second-order expansion, l = 1..12. Written AFTER the blind read.

The blind reader of the short paper reported the closed form
  F_2(a) = 2 pi (l-1)(l+2)/(2l+1) a^2 - 8 pi (l-1)/(2l+1) a + 12 pi (l-1)/((2l+1)(2l+3))
for Area - 2 Volume on r = 1 + eps a P_l in the Weyl solution U = eps r^l P_l. This script checks it against
d1_weyl.F_series (the registered expansion code, unchanged) for l = 1..6, and that its stationary value
equals the entropy formula c_l (l+1)/2 * int psi^2 identically.
Run: python3 d3_weyl_closed_form.py   (exit code 1 on any FAIL)
"""
import sys
import sympy as sp
from sympy import Rational as Q, pi
import d1_weyl as d

a = d.a
l_ = sp.symbols('l', positive=True)
closed = 2 * pi * (l_ - 1) * (l_ + 2) / (2 * l_ + 1) * a**2 - 8 * pi * (l_ - 1) / (2 * l_ + 1) * a + 12 * pi * (l_ - 1) / ((2 * l_ + 1) * (2 * l_ + 3))
d.FAILS.clear()
ok_all = True
for l in range(1, 7):      # the full expansion is slow beyond l = 6; l = 7..12 is covered by ../../blindread/work/check6_weyl_highl.py
    F0, F1, F2 = d.F_series(l)
    ok = sp.simplify(F0 - 4 * pi / 3) == 0 and F1 == 0 and sp.simplify(sp.expand(F2) - sp.expand(closed.subs(l_, l))) == 0
    ok_all &= ok
    print(('PASS ' if ok else 'FAIL ') + f'l={l}: F_2(a) equals the closed form')
d.check('closed form of F_2(a) holds for l = 1..6 (full expansion)', ok_all)

# symbolic in l: stationary point, stationary value, and equality with the entropy formula
astar = sp.solve(sp.diff(closed, a), a)[0]
d.check('stationary shift a* = 2/(l+2), symbolic in l', sp.simplify(astar - 2 / (l_ + 2)) == 0)
F2s = sp.simplify(closed.subs(a, astar))
area_side = sp.simplify(-2 * F2s)                                  # units eps^2/8G
want = 8 * pi * l_ * (l_ - 1) / ((l_ + 2) * (2 * l_ + 1) * (2 * l_ + 3))
d.check('area side = 8 pi l(l-1)/((l+2)(2l+1)(2l+3)), symbolic in l', sp.simplify(area_side - want) == 0)
# entropy: c_l (l+1)/2 * l^2 (l-1)^2 * int_B u^2, with c_l = 4/((l+2)(l+1)l(l-1)) and int_B u^2 = 4 pi/((2l+1)(2l+3))
entropy = 4 / ((l_ + 2) * (l_ + 1) * l_ * (l_ - 1)) * (l_ + 1) / 2 * l_**2 * (l_ - 1)**2 * 4 * pi / ((2 * l_ + 1) * (2 * l_ + 3))
d.check('area side - entropy = 0, symbolic in l (given the closed form)', sp.simplify(area_side - entropy) == 0)
# plant: a wrong constant term must break the l-by-l comparison
bad = closed + pi * (l_ - 1) / (2 * l_ + 1)
F0, F1, F2 = d.F_series(3)
d.check('plant (shifted constant term) is caught at l = 3', sp.simplify(sp.expand(F2) - sp.expand(bad.subs(l_, 3))) != 0)
print('failed:', d.FAILS)
sys.exit(1 if d.FAILS else 0)
