# CHARGE-REF, unregistered extra: step 5 for the time shift, proved pointwise.
# On the slice t = e X(x) of a first-order solution with TT data (h, K), d_t h = 2K, the induced metric is
#   delta + e h + e^2 (h2 + 2 X K - dX dX) + O(e^3)        (TT gauge, no first-order lapse or shift).
# Claim checked: R_lin[2 X K - dX dX] = 2 K.Hess X + |Hess X|^2 - (Lap X)^2 pointwise, R_lin[m] = d_i d_j m_ij - Lap tr m,
# so  R_lin[h2 + 2XK - dXdX] + R2[h] = |K + Hess X|^2 - (tr(K + Hess X))^2 : the moved slice's area side is the
# fixed-region formula on the data (h, K + Hess X).  Usage: python ref_slice.py
import sympy as sp

x, y, z = sp.symbols("x y z")
XS = (x, y, z)
I3 = range(3)
phi = x * y + z**2 - x**2 + x * y * z + x**3 - 3 * x * y**2          # harmonic, degrees 2 and 3 mixed
g = [sp.diff(phi, c) for c in XS]
v = [XS[1] * g[2] - XS[2] * g[1], XS[2] * g[0] - XS[0] * g[2], XS[0] * g[1] - XS[1] * g[0]]
s = x**2 + y**2 + z**2
# a TT tensor that is neither a Hessian nor a symmetrised gradient: K = curl(s (grad v + grad v^T)) made TT by hand is
# not needed; the identity uses only tr K = 0 and div K = 0, so take K = Hess phi + (grad v + grad v^T) + odd TT piece below.
M = [[sp.diff(v[j], XS[i]) + sp.diff(v[i], XS[j]) for j in I3] for i in I3]
phi2 = x * y
g2 = [sp.diff(phi2, c) for c in XS]
v2 = [XS[1] * g2[2] - XS[2] * g2[1], XS[2] * g2[0] - XS[0] * g2[2], XS[0] * g2[1] - XS[1] * g2[0]]
M2 = [[sp.diff(v2[j], XS[i]) + sp.diff(v2[i], XS[j]) for j in I3] for i in I3]
# odd TT with a non-constant profile at l = 2 (eta = s): S = g1 (x v + v x) + g2 M, g1 = -2/(l+2), g2 = (s + 2 s/(l+2))/(l-1)
S2 = [[sp.Rational(-1, 2) * (XS[i] * v2[j] + v2[i] * XS[j]) + sp.Rational(3, 2) * s * M2[i][j] for j in I3] for i in I3]
K = [[sp.expand(sp.diff(phi, XS[i], XS[j]) + M[i][j] + S2[i][j]) for j in I3] for i in I3]
trK = sp.expand(sum(K[i][i] for i in I3))
divK = [sp.expand(sum(sp.diff(K[i][j], XS[i]) for i in I3)) for j in I3]
okTT = trK == 0 and all(d == 0 for d in divK)
X = x**2 * y + z**3 - 2 * x * y * z + s * x + y**4 + phi * (1 + s) / 3 + x * z
gX = [sp.diff(X, c) for c in XS]
m = [[2 * X * K[i][j] - gX[i] * gX[j] for j in I3] for i in I3]
Rlin = sp.expand(sum(sp.diff(m[i][j], XS[i], XS[j]) for i in I3 for j in I3) - sum(sp.diff(sum(m[i][i] for i in I3), c, 2) for c in XS))
H = [[sp.diff(X, XS[i], XS[j]) for j in I3] for i in I3]
lapX = sum(H[i][i] for i in I3)
rhs = sp.expand(sum(2 * K[i][j] * H[i][j] + H[i][j]**2 for i in I3 for j in I3) - lapX**2)
ok = okTT and sp.expand(Rlin - rhs) == 0 and rhs != 0 and sp.expand(lapX) != 0
print(("PASS" if ok else "FAIL") + f" E01 (unregistered) K TT: {okTT}; R_lin[2XK - dXdX] - (2 K.HessX + |HessX|^2 - (LapX)^2) = {sp.expand(Rlin - rhs)}; S2 non-trivial: {any(sp.expand(S2[i][j]) != 0 for i in I3 for j in I3)}")
# plant: wrong sign of the dX dX term
m = [[2 * X * K[i][j] + gX[i] * gX[j] for j in I3] for i in I3]
Rb = sp.expand(sum(sp.diff(m[i][j], XS[i], XS[j]) for i in I3 for j in I3) - sum(sp.diff(sum(m[i][i] for i in I3), c, 2) for c in XS))
print(("PLANT CAUGHT" if sp.expand(Rb - rhs) != 0 else "PLANT MISSED") + " EP1 induced metric with + dX dX")
