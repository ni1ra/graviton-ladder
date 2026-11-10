"""s1: the modular side on the scope family, by the direct route (THEOREM.md, Part M).

Exact. For l = 2..5 and four harmonics each it checks every algebraic identity the direct proof uses:
  M1  K is TT and harmonic, K.x = (l-1)/2 v, x.K.x = 0 (so K_nn = 0 on every sphere)
  M2  the datum is stationary: E = 0, B = -sym curl K, curl K symmetric, d_t E = -sym curl B = 0
  M3  modular flow of the curvature: dE = 2 sym(x cross B) = (l+1) K, dB = 0
  M4  geometric flow of the metric: Ric_3(2 N K) = (l+2) K
  M5  the cut vector: h_v = c 2NK on B, 0 outside, c = (l+1)/(l+2):
        interior linearised scalar curvature 0; surface part of the constraint 2c K_nn = 0;
        its curvature inside B is the modular flow of the datum's curvature
  M6  representative independence: int_B Hess f . N K = 0 for every polynomial f
  M7  S = (1/2) Omega_B(a, delta_C) = c A, A = int_B N |K|^2 = flux/((2l+1)(2l+3))
  M8  the explicit smearing field: C(r^2 K) = alpha_l curl K  (wave69's Cotton-York operator)
Plants: a wrong lapse, a wrong flow coefficient, an electric datum (K_nn != 0).
"""
import sympy as sp
from lib import *

R = Report('s1_modular_side')
xv = xvec()
rng_polys = [x**3 * y - 2 * z**2 + x * y * z + 5, x**4 + y**2 * z**3 - 3 * x * z, z**5 - x**2 * y**2 + y]

for l in range(2, 6):
    c = sp.Rational(l + 1, l + 2)
    for name, phi in harmonics(l).items():
        tag = f'l={l} {name}: '
        R.check(tag + 'phi harmonic, degree l', sp.expand(lap(phi)) == 0 and sp.Poly(phi, x, y, z).is_homogeneous)
        v, K = magnetic(phi)
        # M1
        R.check(tag + 'M1 tr K = 0', sp.expand(K.trace()) == 0)
        R.check(tag + 'M1 div K = 0', is_zero(div(K)))
        R.check(tag + 'M1 Lap K = 0', is_zero(tlap(K)))
        R.check(tag + 'M1 K.x = (l-1)/2 v', is_zero(K * xv - sp.Rational(l - 1, 2) * v))
        R.check(tag + 'M1 x.K.x = 0', sp.expand((xv.T * K * xv)[0]) == 0)
        # M2
        cK = mexp(curl(K))
        R.check(tag + 'M2 curl K symmetric', is_zero(cK - cK.T))
        Bf = -sym(cK)
        R.check(tag + 'M2 B is TT', sp.expand(Bf.trace()) == 0 and is_zero(div(Bf)))
        R.check(tag + 'M2 d_t E = -sym curl B = 0', is_zero(sym(curl(Bf))))
        # M3  (E = 0, d_t B = sym curl E = 0)
        dE = mexp(2 * sym(xcross(Bf)))
        R.check(tag + 'M3 dE = (l+1) K', is_zero(dE - (l + 1) * K))
        # M4
        ric = mexp(ric3(mexp(2 * N * K)))
        R.check(tag + 'M4 Ric_3(2NK) = (l+2) K', is_zero(ric - (l + 2) * K))
        # M5
        H = mexp(c * 2 * N * K)
        R.check(tag + 'M5 interior R_lin[h_v] = 0', sp.expand(rlin(H)) == 0)
        # surface part of d_i d_j h_ij - Lap tr h for h = H 1_B with H|_S = 0 is (-n.div H + d_n tr H) delta_S
        surf = sp.expand(-(xv.T * div(H))[0] + sum(X[i] * sp.diff(H.trace(), X[i]) for i in range(3)))
        R.check(tag + 'M5 surface part of the constraint = 2c x.K.x = 0', surf == 0)
        R.check(tag + 'M5 curvature of h_v inside B = modular flow (E: (l+1)K, B: 0)', is_zero(mexp(ric3(H)) - dE))
        # M6
        ok6 = all(sp.simplify(ball_int(contract(mat(lambda i, j: sp.diff(f, X[i], X[j])), N * K))) == 0 for f in rng_polys)
        R.check(tag + 'M6 int_B Hess f . N K = 0 for three polynomials f', ok6)
        # M7
        K2 = sp.expand(contract(K, K))
        A = sp.simplify(ball_int(N * K2))
        flux = sp.simplify(sphere_int(K2))
        S = sp.simplify(sp.Rational(1, 2) * ball_int(contract(K, H)))
        R.check(tag + 'M7 A = flux/((2l+1)(2l+3))', sp.simplify(A - flux / ((2 * l + 1) * (2 * l + 3))) == 0)
        R.check(tag + 'M7 S = (l+1)/(l+2) A', sp.simplify(S - c * A) == 0)
        R.check(tag + 'M7 A - S = A/(l+2)', sp.simplify(A - S - A / (l + 2)) == 0)
        R.value(tag + 'A', A)
        R.value(tag + 'S', S)
        # M8
        Cf = mexp(cotton(mexp(r2 * K)))
        # solve alpha from one nonzero component, then verify on all
        alpha = None
        for i in range(3):
            for j in range(3):
                if cK[i, j] != 0 and alpha is None:
                    alpha = sp.cancel(Cf[i, j] / cK[i, j])
        R.check(tag + 'M8 C(r^2 K) = -2(l+2) curl K', alpha == -2 * (l + 2) and is_zero(Cf - alpha * cK), alpha)
        R.value(tag + 'alpha', alpha)
        # plants
        Nw = N + sp.Rational(1, 10) * (1 - r2)**2
        ricw = mexp(ric3(mexp(2 * Nw * K)))
        R.plant(tag + 'P1 wrong lapse N + (1-r^2)^2/10: Ric no longer (l+2)K', not is_zero(ricw - (l + 2) * K))
        dEw = mexp(sym(xcross(Bf)))
        R.plant(tag + 'P2 flow coefficient 1 instead of 2: dE != (l+1)K', not is_zero(dEw - (l + 1) * K))
        Sw = sp.simplify(sp.Rational(1, 2) * ball_int(contract(K, mexp(sp.Rational(l, l + 2) * 2 * N * K))))
        R.plant(tag + 'P3 wrong rung l/(l+2): S != c A', sp.simplify(Sw - c * A) != 0)
    # P4: electric datum K = Hess phi has K_nn != 0, so the cut vector violates the constraint and M6 fails
    phi = harmonics(l)['m=l cos']
    Ke = mat(lambda i, j: sp.diff(phi, X[i], X[j]))
    He = mexp(2 * N * Ke)
    surf = sp.expand(-(xv.T * div(He))[0])
    f = rng_polys[0] + phi
    g6 = sp.simplify(ball_int(contract(mat(lambda i, j: sp.diff(f, X[i], X[j])), N * Ke)))
    R.plant(f'l={l} P4 electric datum: surface constraint term nonzero', surf != 0, surf)
    R.plant(f'l={l} P4 electric datum: pairing depends on the representative', g6 != 0, g6)

R.value('alpha closed form', '-2(l+2), so the coherent vector is v_f with f = -(2/(l+2)) chi(r) r^2 K (wave69 Step 1: B[v_f] = C f / 4)')
R.write()
