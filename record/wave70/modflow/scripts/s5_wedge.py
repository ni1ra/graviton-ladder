"""S5 (added after the PREREG, labelled post hoc; supports Lemma 5 only): the boost generator on (E, B).
M = zeroth-order part of delta_X for X = x_1 d_t + t d_1 on symmetric (E, B) (12 components).
Checks: M is diagonalisable with eigenvalues in {0, +-1, +-2}; exp(i pi M^T) is real and equals
(E, B) -> (R E R, R B R), R = diag(-1, 1, 1): the tensor action of j = diag(-1, -1, 1, 1)."""
import json, sys
import sympy as sp
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from common import killing_fields, delta_matrix_part
xs = sp.symbols('x0 x1 x2 x3')
Es = sp.symbols('e0:6'); Bs = sp.symbols('b0:6')
def symm(v): return [[v[0], v[1], v[2]], [v[1], v[3], v[4]], [v[2], v[4], v[5]]]
def flat(M): return [M[0][0], M[0][1], M[0][2], M[1][1], M[1][2], M[2][2]]
vars_ = list(Es) + list(Bs)
dE, dB = delta_matrix_part(killing_fields(xs)['N1'], xs, symm(Es), symm(Bs))
M = sp.Matrix([[sp.diff(sp.expand(o), v) for v in vars_] for o in flat(dE) + flat(dB)])
out = {'eigenvalues': {str(k): v for k, v in M.eigenvals().items()}, 'diagonalisable': M.is_diagonalizable()}
Dj = sp.simplify((sp.I * sp.pi * M).exp())            # action on the fields; the test-function action is its transpose
R = sp.diag(-1, 1, 1)
RE = (R * sp.Matrix(symm(Es)) * R).tolist(); RB = (R * sp.Matrix(symm(Bs)) * R).tolist()
target = sp.Matrix([[sp.diff(o, v) for v in vars_] for o in flat(RE) + flat(RB)])
out['exp_i_pi_M_is_real'] = all(sp.im(e) == 0 for e in Dj)
out['exp_i_pi_M_equals_R_conjugation'] = sp.simplify(Dj - target) == sp.zeros(12, 12)
out['exp_i_pi_MT_equals_its_transpose_target'] = sp.simplify(Dj.T - target.T) == sp.zeros(12, 12)
json.dump(out, open(__file__.rsplit('/', 2)[0] + '/out/s5_out.json', 'w'), indent=1)
print(out)
