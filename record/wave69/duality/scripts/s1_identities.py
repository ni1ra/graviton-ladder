"""S1: exact symbolic identities behind the helicity-2 duality proof (PREREG section 3, S1).

Constant-coefficient differential operators on R^3 are polynomials in commuting symbols d1,d2,d3
(d_j stands for the partial derivative). Every check is an exact polynomial identity.
"""
import json
import sys
import sympy as sp
from sympy import LeviCivita as eps

d = sp.symbols("d1 d2 d3")
R3 = range(3)
lap = sum(x**2 for x in d)
I3 = sp.eye(3)


def Z():
    return sp.zeros(3, 3)


def curl1(A):
    """curl on the first index: (curl A)_ij = eps_ikl d_k A_lj"""
    M = Z()
    for i in R3:
        for j in R3:
            M[i, j] = sum(eps(i, k, l) * d[k] * A[l, j] for k in R3 for l in R3)
    return M.expand()


def inc(s):
    """inc(s)_ij = eps_ikl eps_jmn d_k d_m s_ln"""
    M = Z()
    for i in R3:
        for j in R3:
            M[i, j] = sum(
                eps(i, k, l) * eps(j, m, n) * d[k] * d[m] * s[l, n]
                for k in R3 for l in R3 for m in R3 for n in R3
            )
    return M.expand()


def cotton(s):
    G = inc(s)
    return curl1(G - sp.Rational(1, 2) * I3 * G.trace())


def cotton_bad(s):  # plant: trace correction dropped
    return curl1(inc(s))


def div1(T):
    """(div T)_j = d_i T_ij"""
    return sp.Matrix([sum(d[i] * T[i, j] for i in R3) for j in R3]).expand()


def killing(v):
    M = Z()
    for i in R3:
        for j in R3:
            M[i, j] = sp.Rational(1, 2) * (d[i] * v[j] + d[j] * v[i])
    return M


def is_zero(M):
    return all(sp.expand(e) == 0 for e in list(M))


# generic symmetric and generic full tensors
s_sym = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f"s{min(i,j)}{max(i,j)}"))
phi_full = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f"p{i}{j}"))
v = sp.symbols("v1 v2 v3")
psi = sp.Symbol("psi")
res = {}


def operator_checks(Cop, tag):
    T = Cop(s_sym)
    out = {}
    out["a_symmetric"] = is_zero(T - T.T)
    out["a_traceless"] = sp.expand(T.trace()) == 0
    out["a_divergence_free"] = is_zero(div1(T))
    out["b_kills_Killing_Kv"] = is_zero(Cop(killing(v)))
    out["b_kills_pure_trace"] = is_zero(Cop(psi * I3))
    # c. formal self-adjointness: for an odd-order real operator, C^dagger = -C(d)^T,
    # so C is formally self-adjoint iff sum_ij g_ij (C f)_ij + f_ij (C g)_ij = 0 as a polynomial
    g_sym = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f"g{min(i,j)}{max(i,j)}"))
    Tg = Cop(g_sym)
    pair = sum(g_sym[i, j] * T[i, j] + s_sym[i, j] * Tg[i, j] for i in R3 for j in R3)
    out["c_formally_self_adjoint"] = sp.expand(pair) == 0
    res[tag] = out
    return out


ok_good = operator_checks(cotton, "C")
ok_bad = operator_checks(cotton_bad, "C_bad_plant")

# d. on TT tensors curl^2 = -lap, and curl lap P_TT f = c * C f.
# P_TT is non-local: multiply through by lap^2. With Pi = lap*I - d d^T (lap times the transverse projector),
# lap^2 P_TT f = Pi f Pi - (1/2) Pi tr(Pi f).
dd = sp.Matrix(3, 3, lambda i, j: d[i] * d[j])
Pi = lap * I3 - dd


def PTT_times_lap2(f):
    return (Pi * f * Pi - sp.Rational(1, 2) * Pi * (Pi * f).trace()).expand()


TT = PTT_times_lap2(s_sym)  # lap^2 * (TT part of s)
dchk = {}
dchk["PTT_output_is_TT"] = is_zero(TT - TT.T) and sp.expand(TT.trace()) == 0 and is_zero(div1(TT))
cc = curl1(curl1(TT))
dchk["curl_curl_equals_minus_lap_on_TT"] = is_zero(cc + lap * TT)
dchk["curl_of_TT_is_symmetric_TT"] = is_zero(curl1(TT) - curl1(TT).T) and is_zero(div1(curl1(TT)))
# curl lap P_TT f  vs  C f :  lap^2 * (curl lap P_TT f) = lap * curl(TT) ; compare with c * lap^2 * C f
lhs = (lap * curl1(TT)).expand()
rhs = (lap**2 * cotton(s_sym)).expand()
# find the constant c from one entry, then test all
c = sp.Symbol("c")
entry = None
for i in R3:
    for j in R3:
        if sp.expand(rhs[i, j]) != 0 and entry is None:
            entry = (i, j)
pl = sp.Poly(lhs[entry], *d, *sorted(s_sym.free_symbols, key=str))
pr = sp.Poly(rhs[entry], *d, *sorted(s_sym.free_symbols, key=str))
cval = sp.Rational(pl.coeffs()[0], 1) / pr.coeff_monomial(pl.monoms()[0])
dchk["c_constant_in_curl_lap_PTT_eq_c_times_C"] = str(cval)
dchk["curl_lap_PTT_f_equals_c_C_f"] = is_zero(lhs - cval * rhs)
# C f depends only on the TT part:  C(P_TT f) = C f
dchk["C_of_TT_part_equals_C"] = is_zero(cotton(TT) - lap**2 * cotton(s_sym))
res["d"] = dchk

# e. Lemma B chain
echk = {}
A = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f"A{i}{j}"))
Tc = curl1(A)
# e1: trace of curl A is the divergence of a_k = eps_kli A_li
a = [sum(eps(k, l, i) * A[l, i] for l in R3 for i in R3) for k in R3]
echk["e1_trace_curlA_is_div_of_dual"] = sp.expand(Tc.trace() - sum(d[k] * a[k] for k in R3)) == 0
# e2: A - grad v (grad on the first index) has the same curl, and its dual vector is a - curl v
gradv = sp.Matrix(3, 3, lambda l, j: d[l] * v[j])
echk["e2_curl_of_gradient_vanishes"] = is_zero(curl1(gradv))
dual_grad = [sum(eps(k, l, i) * gradv[l, i] for l in R3 for i in R3) for k in R3]
curlv = [sum(eps(k, m, n) * d[m] * v[n] for m in R3 for n in R3) for k in R3]
echk["e2_dual_of_gradient_is_curl_v"] = all(sp.expand(dual_grad[k] - curlv[k]) == 0 for k in R3)
# e3: antisymmetric part of curl A (A symmetric) is  d_j A_pj - d_p trA
As = s_sym
Ts = curl1(As)
anti = [sum(eps(p, i, j) * Ts[i, j] for i in R3 for j in R3) for p in R3]
target = [sum(d[j] * As[p, j] for j in R3) - d[p] * As.trace() for p in R3]
echk["e3_symmetry_of_curlA_iff_div_condition"] = all(sp.expand(anti[p] - target[p]) == 0 for p in R3)
# e3b: the same identity for a general (non-symmetric) potential, as used in the Calabi step
Tg_ = curl1(A)
anti_g = [sum(eps(p, i, j) * Tg_[i, j] for i in R3 for j in R3) for p in R3]
target_g = [sum(d[j] * A[p, j] for j in R3) - d[p] * A.trace() for p in R3]
echk["e3b_same_identity_general_potential"] = all(sp.expand(anti_g[p] - target_g[p]) == 0 for p in R3)
# e4: Calabi step. With Mt_pj = eps_jmn d_m phi_pn (phi arbitrary), M = Mt - (1/2) delta tr Mt,
#     curl1(M) = inc(sym part of phi)   exactly (the antisymmetric part of phi drops out)
Mt = sp.Matrix(3, 3, lambda p, j: sum(eps(j, m, n) * d[m] * phi_full[p, n] for m in R3 for n in R3))
M = Mt - sp.Rational(1, 2) * I3 * Mt.trace()
phis = (phi_full + phi_full.T) / 2
echk["e4_curl_M_equals_inc_of_symmetric_part"] = is_zero(curl1(M) - inc(phis))
echk["e4_inc_is_symmetric_and_div_free"] = is_zero(inc(s_sym) - inc(s_sym).T) and is_zero(div1(inc(s_sym)))
# e5: G = A - delta trA  <=>  A = G - (1/2) delta trG
Gm = As - I3 * As.trace()
echk["e5_inverse_trace_shift"] = is_zero(Gm - sp.Rational(1, 2) * I3 * Gm.trace() - As)
res["e"] = echk

# f. conformal Killing currents: d_i (T_ij xi_j) = (d_i T_ij) xi_j + T_ij d_i xi_j ; the second term
# vanishes for symmetric traceless T iff the symmetric traceless part of grad xi vanishes.
x = sp.symbols("x1 x2 x3")
a_, b_ = sp.symbols("a1:4"), sp.symbols("b1:4")
kap = sp.Symbol("kappa")
w = sp.symbols("w1:4")
xsq = sum(xi**2 for xi in x)
bx = sum(b_[i] * x[i] for i in R3)
xi = [
    a_[j]
    + sum(eps(j, m, n) * w[m] * x[n] for m in R3 for n in R3)
    + kap * x[j]
    + b_[j] * xsq - 2 * x[j] * bx
    for j in R3
]
grad = sp.Matrix(3, 3, lambda i, j: sp.diff(xi[j], x[i]))
symtf = (grad + grad.T) / 2 - I3 * grad.trace() / 3
res["f"] = {"ckv_10_parameter_family_has_no_symmetric_traceless_gradient": is_zero(symtf.expand()),
            "parameters": 10}

good_pass = all(ok_good.values()) and all(v_ for k_, v_ in dchk.items() if isinstance(v_, bool)) \
    and all(echk.values()) and res["f"]["ckv_10_parameter_family_has_no_symmetric_traceless_gradient"]
plant_caught = not all(ok_bad.values())
res["S1_ALL_PASS"] = bool(good_pass)
res["S1_PLANT_CAUGHT"] = bool(plant_caught)
res["verdict"] = "PASS" if (good_pass and plant_caught) else ("VOID" if good_pass else "FAIL")
json.dump(res, open(sys.argv[1], "w"), indent=1, default=str)
print(json.dumps(res, indent=1, default=str))
