"""s3: the account, and the cross-check of the two routes to the modular side (THEOREM.md, Part W). Exact.

For psi = x.B.x, B = -sym curl K, and E_w[psi] = (1/2) int_B N |grad psi|^2 + (w/2) int_B psi^2 (psi-dot = 0 here):
  W1  transfer route on the scope family: c_l E_1[psi] equals the direct value (l+1)/(l+2) A   [route 1 = route 2]
  W2  A = c_l E_2[psi], so A - S = (c_l/2) int_B psi^2 = A/(l+2)
  W3  the same identity A = c_l E_2[psi] on non-stationary odd momentum data T = f K_l + g P (profile eta = 1 + a r^2)
  W4  uniqueness of the weights: E_w - E_w' = (w - w')/2 int psi^2, and int psi^2 > 0
Plants: c_l with l -> l+1; weight 3/2 on the area side.
"""
import sympy as sp
from lib import *

R = Report('s3_account')
xv = xvec()
half = sp.Rational(1, 2)


def cl(l):
    return sp.Rational(4) * sp.factorial(l - 2) / sp.factorial(l + 2)


def Ew(psi, w):
    g2 = sp.expand(sum(sp.diff(psi, v)**2 for v in X))
    return half * ball_int(sp.expand(N * g2)) + sp.Rational(w) / 2 * ball_int(sp.expand(psi**2))


for l in (2, 3, 4):
    for name in ('m=l cos', 'm=0', 'mixed'):
        phi = harmonics(l)[name]
        tag = f'l={l} {name}: '
        v, K = magnetic(phi)
        A = sp.simplify(ball_int(sp.expand(N * contract(K, K))))
        Bf = -sym(curl(K))
        psi = sp.expand((xv.T * Bf * xv)[0])
        R.check(tag + 'psi harmonic of degree l', sp.expand(lap(psi)) == 0 and psi != 0)
        R.check(tag + 'W0 psi = (l-1) l (l+1)/2 phi', sp.expand(psi - sp.Rational((l - 1) * l * (l + 1), 2) * phi) == 0)
        R.check(tag + 'W0 oint |K|^2 = (l-1)(2l+1) l (l+1)/2 oint phi^2',
                sp.simplify(sphere_int(sp.expand(contract(K, K))) - sp.Rational((l - 1) * (2 * l + 1) * l * (l + 1), 2) * sphere_int(sp.expand(phi**2))) == 0)
        S1 = sp.simplify(cl(l) * Ew(psi, 1))
        R.check(tag + 'W1 c_l E_1[psi] = (l+1)/(l+2) A  (transfer formula = direct value)', sp.simplify(S1 - sp.Rational(l + 1, l + 2) * A) == 0, (S1, A))
        R.check(tag + 'W2 c_l E_2[psi] = A', sp.simplify(cl(l) * Ew(psi, 2) - A) == 0)
        gap = sp.simplify(cl(l) / 2 * ball_int(sp.expand(psi**2)))
        R.check(tag + 'W2 (c_l/2) int psi^2 = A/(l+2) > 0', sp.simplify(gap - A / (l + 2)) == 0 and gap > 0)
        R.plant(tag + 'P c_{l+1} in place of c_l', sp.simplify(cl(l + 1) * Ew(psi, 1) - sp.Rational(l + 1, l + 2) * A) != 0)
        R.plant(tag + 'P weight 3/2 on the area side', sp.simplify(cl(l) * Ew(psi, sp.Rational(3, 2)) - A) != 0)
    # W3: non-stationary momentum datum
    phi = harmonics(l)['mixed']
    v, K = magnetic(phi)
    P = mat(lambda i, j: X[i] * v[j] + v[i] * X[j])
    m = l + 2
    for aa in (sp.Rational(-3, 5), sp.Rational(2)):
        eta = 1 + aa * r2
        g = -2 * aa / m                       # g = -eta'/(m r)
        f = 2 * (eta - r2 * g) / (l - 1)
        T = mexp(f * K + g * P)
        R.check(f'l={l} a={aa}: W3 T is TT', sp.expand(T.trace()) == 0 and is_zero(div(T)))
        A = sp.simplify(ball_int(sp.expand(N * contract(T, T))))
        psi = sp.expand((xv.T * (-sym(curl(T))) * xv)[0])
        R.check(f'l={l} a={aa}: W3 A = int N T.T = c_l E_2[psi]', sp.simplify(cl(l) * Ew(psi, 2) - A) == 0, A)
        R.check(f'l={l} a={aa}: W3 A - c_l E_1 = (c_l/2) int psi^2 > 0', sp.simplify(A - cl(l) * Ew(psi, 1) - cl(l) / 2 * ball_int(sp.expand(psi**2))) == 0)
        R.value(f'l={l} a={aa}: S/A on this datum', sp.simplify(cl(l) * Ew(psi, 1) / A))
R.write()
