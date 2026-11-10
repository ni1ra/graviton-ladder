"""w72-decider, second test. The Killing limit of the ladder's two rates.

For TT harmonic K and lapse N:
  geometric rate  G = Ric3(2 N K)
  modular rate    M = 2 sym(grad N x curl K)
Claim: G - M = 2 sym(K . HessN) - (Lap N) K.
Wedge (N = x1): G = M.  Ball (N = (1-r^2)/2): G - M = K.
Run: python3 d2_killing.py   (exit code 1 on any FAIL)
"""
import sys, random
import sympy as sp
from sympy import Rational as Q

X = sp.symbols('x1 x2 x3', real=True)
FAILS = []


def check(name, ok, detail=''):
    print(('PASS ' if ok else 'FAIL ') + name + (('  ' + str(detail)) if detail else ''))
    if not ok:
        FAILS.append(name)


def lap(f):
    return sum(sp.diff(f, v, 2) for v in X)


def random_harmonic(deg, rnd):
    """Random harmonic polynomial of degree deg: real/imag parts of complex powers along random null directions."""
    phi = 0
    for _ in range(3):
        # null vector n = a + i b with a.b = 0, |a| = |b|: use (1, i, 0) rotated by a random integer permutation/sign
        perm = rnd.sample(range(3), 3)
        sg = [rnd.choice((-1, 1)) for _ in range(3)]
        w = sg[0] * X[perm[0]] + sp.I * sg[1] * X[perm[1]]
        term = sp.expand(w**(deg - 1) * (sg[2] * X[perm[2]])) if deg >= 1 and rnd.random() < 0.5 else sp.expand(w**deg)
        phi += rnd.randint(1, 4) * (sp.re(term) if rnd.random() < 0.5 else sp.im(term))
    phi = sp.expand(phi)
    assert sp.expand(lap(phi)) == 0
    return phi


def curl(T):
    return sp.Matrix(3, 3, lambda i, j: sum(sp.LeviCivita(i, k, m) * sp.diff(T[m, j], X[k])
                                            for k in range(3) for m in range(3)))


def sym(T):
    return (T + T.T) / 2


def K_scope(phi):
    xs = sp.Matrix(X)
    v = xs.cross(sp.Matrix([sp.diff(phi, a) for a in X]))
    return sp.Matrix(3, 3, lambda i, j: Q(1, 2) * (sp.diff(v[i], X[j]) + sp.diff(v[j], X[i])))


def hess(phi):
    return sp.Matrix(3, 3, lambda i, j: sp.diff(phi, X[i], X[j]))


def is_tt_harmonic(K):
    z = sp.expand(K.trace()) == 0
    z &= all(sp.expand(sum(sp.diff(K[i, j], X[i]) for i in range(3))) == 0 for j in range(3))
    z &= all(sp.expand(lap(K[i, j])) == 0 for i in range(3) for j in range(3))
    z &= (K - K.T).applyfunc(sp.expand) == sp.zeros(3, 3)
    return z


def ric3(h):
    tr = h.trace()
    div = [sum(sp.diff(h[i, j], X[i]) for i in range(3)) for j in range(3)]
    return sp.Matrix(3, 3, lambda i, j: Q(1, 2) * (-lap(h[i, j]) + sp.diff(div[j], X[i]) + sp.diff(div[i], X[j])
                                                   - sp.diff(tr, X[i], X[j])))


def cross_vec_tensor(g, T):
    return sp.Matrix(3, 3, lambda i, j: sum(sp.LeviCivita(i, k, m) * g[k] * T[m, j] for k in range(3) for m in range(3)))


def rates(N, K, coeff=2):
    g = [sp.diff(N, a) for a in X]
    G = ric3(2 * N * K).applyfunc(sp.expand)
    M = (coeff * sym(cross_vec_tensor(g, curl(K)))).applyfunc(sp.expand)
    H = hess(N)
    rhs = (2 * sym(K * H) - lap(N) * K).applyfunc(sp.expand)
    return G, M, rhs


def zero(M):
    return M.applyfunc(sp.expand) == sp.zeros(3, 3)


if __name__ == '__main__':
    rnd = random.Random(72)
    r2 = sum(v**2 for v in X)
    N_wedge = X[0]
    N_ball = (1 - r2) / 2
    data = []
    for deg in (2, 3, 4, 5):
        for trial in range(2):
            phi = random_harmonic(deg, rnd)
            Ks = K_scope(phi)
            fam = [('scope', Ks), ('hess', hess(phi)), ('curl', curl(Ks))]
            for name, K in fam:
                K = K.applyfunc(sp.expand)
                if K == sp.zeros(3, 3):
                    continue
                assert is_tt_harmonic(K), (name, deg)
                data.append((f'{name} deg(phi)={deg} #{trial}', K))
    print(f'{len(data)} TT harmonic data')
    wedge_ok = ball_ok = gen_ok = plant_caught = True
    for name, K in data:
        G, M, rhs = rates(N_wedge, K)
        wedge_ok &= zero(G - M)
        G, M, rhs = rates(N_ball, K)
        ball_ok &= zero(G - M - K) and zero(G - M - rhs)
        Nq = sum(rnd.randint(-3, 3) * a * b for a in X for b in X) + sum(rnd.randint(-3, 3) * a for a in X) + rnd.randint(-3, 3)
        G, M, rhs = rates(Nq, K)
        gen_ok &= zero(G - M - rhs)
        G, M, rhs = rates(N_wedge, K, coeff=3)
        plant_caught &= not zero(G - M)
    check('wedge (Killing lapse): geometric rate = modular rate on every datum', wedge_ok)
    check('ball: geometric rate - modular rate = K on every datum', ball_ok)
    check('general lapse: G - M = 2 sym(K.HessN) - (Lap N) K', gen_ok)
    check('plant (flow coefficient 2 -> 3) breaks the wedge equality on every datum', plant_caught)
    # the scope family's own integers, as a power check of this code
    for l in (2, 3, 4):
        phi = random_harmonic(l, rnd)
        K = K_scope(phi).applyfunc(sp.expand)
        G, M, _ = rates(N_ball, K)
        check(f'scope family l={l}: M = (l+1)K and G = (l+2)K', zero(M - (l + 1) * K) and zero(G - (l + 2) * K))
    bin_ = 'KILLING-OK' if (wedge_ok and ball_ok and gen_ok) else ('DEFEATED' if not wedge_ok else 'IDENTITY-WRONG')
    print('bin:', bin_, '| failed:', FAILS)
    sys.exit(1 if FAILS else 0)
