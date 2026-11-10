"""w55-gravlattice2 theory side: lower-degree controls. Find polynomial potentials of the lowest degree whose chain field
psi = x.C.x vanishes on the unit sphere: psi = (1 - r^2) * P with P a harmonic polynomial (so sigma = 0 exactly).
E-type: Y with psi_even(E(Y)) = (1-r^2) P ; B-type: Z with psi_odd(B(Z)) = (1-r^2) P."""
import sympy as sp, sys, time
from pred import X, r2, monos, TCOMP, ric, symcurl, inc, P2, P2b, P3

xv = sp.Matrix(X)


def solve_scalar(target, op, deg):
    ms = monos(deg)
    cs = sp.symbols('a0:%d' % (6 * len(ms)))
    T = sp.zeros(3, 3); n = 0
    for (i, j) in TCOMP:
        e = sum(cs[n + k] * ms[k] for k in range(len(ms))); n += len(ms)
        T[i, j] = e; T[j, i] = e
    psi = sp.expand((xv.T * op(T) * xv)[0])
    eqs = sp.Poly(sp.expand(psi - target), *X).coeffs()
    sol = sp.solve(eqs, cs, dict=True)
    if not sol:
        return None
    Ts = T.subs(sol[0]).subs({c: 0 for c in cs})
    assert sp.expand((xv.T * op(Ts) * xv)[0] - target) == 0
    return Ts


opE = lambda T: ric(symcurl(T))
opB = lambda T: -symcurl(inc(T) / 2 - inc(T).trace() / 4 * sp.eye(3))
if __name__ == '__main__':
    out = {}
    for name, P, op, degs in (('d1_E2', P2, opE, (4, 5)), ('d2_B2', P2b, opB, (4, 5)), ('d3_E3', P3, opE, (5, 6))):
        for d in degs:
            t0 = time.time()
            T = solve_scalar(sp.expand((1 - r2) * P), op, d)
            print(name, 'deg', d, 'solved' if T is not None else 'no solution', round(time.time() - t0, 1), 's', flush=True)
            if T is not None:
                out[name] = [str(T[i, j]) for (i, j) in TCOMP]
                break
    import json
    json.dump(out, open(sys.argv[1], 'w'), indent=1)
