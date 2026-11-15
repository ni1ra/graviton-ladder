"""Lane EVEN-TOWER, Addendum C: the canonical charge (1/2) Omega_B(a, delta^w a) on EVEN TT data, against the
extremised geometric charge.  delta_w and closed4 are copied from charge-id/scripts/charge_id.py.
Usage: python even_canonical.py L [KMAX]"""
import sys, os, runpy
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..", "charge-id", "scripts")))
import polyalg as P
from polyalg import (add, sub, mul, scale, pdiff, diff, sumall, ball, sph, dot, trace, div, tlap, xTx, Tx, ric_lin,
                     rscal_lin, tadd, tscale, tzero, I3, XS, R2, NLAPSE, harmonic, powp)

L = int(sys.argv[1]); KMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 2
src = open(os.path.join(HERE, "even_tower.py")).read()
ns = {"__name__": "lib", "sys": sys, "__file__": os.path.join(HERE, "even_tower.py")}
exec(compile(src.split("# ================================================================== data")[0], "even_tower_head", "exec"), ns)
build_even, unit, xgrad, grad, prof, GD, spoly = (ns[k] for k in ("build_even", "unit", "xgrad", "grad", "prof", "GD", "spoly"))


def delta_w(h, K, lapse=NLAPSE):
    t = P.T_
    Lh, LK = tlap(h), tlap(K)
    H = [[{} for _ in range(4)] for _ in range(4)]
    for i in I3:
        for j in I3:
            H[i + 1][j + 1] = sumall([h[i][j], scale(mul(t, K[i][j]), 2), scale(mul(powp(t, 2), Lh[i][j]), Fr(1, 2)),
                                      scale(mul(powp(t, 3), LK[i][j]), Fr(1, 3))])
    zeta = [sub(lapse, scale(powp(t, 2), Fr(1, 2)))] + [scale(mul(t, XS[i + 1]), -1) for i in I3]
    dz = [[pdiff(zeta[r], m) for m in range(4)] for r in range(4)]
    D = [[{} for _ in range(4)] for _ in range(4)]
    for m in range(4):
        for n in range(4):
            acc = {}
            for r in range(4):
                acc = add(acc, mul(zeta[r], pdiff(H[m][n], r)))
                acc = add(acc, mul(H[r][n], dz[r][m]))
                acc = add(acc, mul(H[m][r], dz[r][n]))
            D[m][n] = add(acc, scale(mul(t, H[m][n]), 2))
    dh = [[P.at_t0(D[i + 1][j + 1]) for j in I3] for i in I3]
    dK = [[sub(scale(P.tcoef(D[i + 1][j + 1], 1), Fr(1, 2)),
               scale(add(P.at_t0(pdiff(D[j + 1][0], i + 1)), P.at_t0(pdiff(D[i + 1][0], j + 1))), Fr(1, 2))) for j in I3] for i in I3]
    return dh, dK, D


def closed4(ha, hb, Ka, Kb, N=NLAPSE):
    gh = sumall([mul(diff(ha[i][j], k + 1), diff(hb[i][j], k + 1)) for i in I3 for j in I3 for k in I3])
    hxa, hxb = Tx(ha), Tx(hb)
    return ball(mul(N, dot(Ka, Kb))) + Fr(1, 4) * ball(mul(N, gh)) - Fr(1, 2) * sph(sumall([mul(x, y) for x, y in zip(hxa, hxb)]))


def line(tag, name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + tag + " " + name + ("  " + str(info) if info else ""), flush=True)


phi = harmonic(L, Fr(3, 5), Fr(4, 5), (1, 2, 3), "re")
nb = KMAX + 1
Tk = [build_even(phi, unit(k)) for k in range(nb)]
Z3 = tzero()
basis = [(T, Z3) for T in Tk] + [(Z3, T) for T in Tk]
n = len(basis)
dl = [delta_w(*b) for b in basis]
om = lambda a, b: ball(sub(dot(a[1], b[0]), dot(b[1], a[0])))
Mom = [[Fr(1, 2) * om(basis[i], dl[j][:2]) for j in range(n)] for i in range(n)]
M4 = [[closed4(basis[i][0], basis[j][0], basis[i][1], basis[j][1]) for j in range(n)] for i in range(n)]
symm = all(Mom[i][j] == Mom[j][i] for i in range(n) for j in range(n))
line("K1a", "Omega_B(a, delta^w b) symmetric on even TT data", symm)
line("K1b", "(1/2) Omega_B(a, delta^w b) = closed form (4) on the even Gram, h-h, K-K, h-K", all(Mom[i][j] == M4[i][j] for i in range(n) for j in range(n)),
     f"A_can(T1 as K)={Mom[nb+1][nb+1]}  A_can(T1 as h)={Mom[1][1]}")
line("K1c", "delta^w H keeps H_00 = 0 and the data stay constraint solutions", all(not d[2][0][0] and not rscal_lin(d[0]) for d in dl))
kap = [xTx(T) for T in Tk]
Gk = [scale(add(xgrad(k), scale(k, 2)), Fr(1, 2)) for k in kap]
okK = all(M4[nb + a][nb + b] - sph(mul(kap[a], kap[b])) / (L * (L - 1)) ==
          ns["BK"](Tk[a], Tk[b]) - sph(mul(kap[a], kap[b])) / (L * (L - 1)) for a in range(nb) for b in range(nb))
line("K2a", "K slot: A_can = int N K.K, so A_can - A_ext^K = + oint kappa kappa/(l(l-1)) (A_ext^K closed form, run01 H3'')", okK,
     f"diff on T1 = {sph(mul(kap[1], kap[1])) / (L * (L - 1))}")
# h slot: A_ext^h = c_l E_1[psi_E] (run01 Q2a)
CL = ns["CL"]
psiE = [xTx(ric_lin(T)) for T in Tk]
gd = lambda a, b: sumall([mul(x, y) for x, y in zip(grad(a), grad(b))])
E1 = [[CL * (Fr(1, 2) * ball(mul(NLAPSE, gd(psiE[a], psiE[b]))) + Fr(1, 2) * ball(mul(psiE[a], psiE[b]))) for b in range(nb)] for a in range(nb)]
okH = all(M4[a][b] - E1[a][b] == -Fr(1, 4) * sph(mul(kap[a], kap[b])) - sph(mul(Gk[a], Gk[b])) / ((L - 1) * (L + 2)) for a in range(nb) for b in range(nb))
line("K2b", "h slot: A_can - A_ext^h = -(1/4) oint (nu.h.nu)^2 - oint G^2/((l-1)(l+2))", okH, f"diff on T1 = {M4[1][1]-E1[1][1]}")
line("K3", "A_can is NOT invariant under the residual TT gauge mode T[1]: its row is non-zero in both slots (A_ext's is zero)",
     any(Mom[0][j] != 0 for j in range(nb)) and any(Mom[nb][nb + j] != 0 for j in range(nb)),
     f"h slot row {[str(Mom[0][j]) for j in range(nb)]} ; K slot row {[str(Mom[nb][nb+j]) for j in range(nb)]}")
# K4 info: transferred flow
for a in range(1, nb):
    dh, dK, _ = dl[nb + a]                     # datum (0, T_a): psi_E = 0, psidot_E = pd
    pd = xTx(ric_lin(tscale(Tk[a], 2)))
    r1 = sub(xTx(ric_lin(dh)), mul(NLAPSE, pd))
    u = mul(spoly(GD(prof(pd, phi))), phi)
    r2 = sub(r1, u)
    print(f"INFO K4 momentum datum T[s^{a}]: psi_E[delta^w h] - N psidot_E is {'zero' if not r1 else 'NON-ZERO'}; minus G_D psidot_E it is "
          f"{'zero' if not r2 else 'non-zero: profile ' + str([str(c) for c in prof(r2, phi)])}")
    dh, dK, _ = dl[a]                          # datum (T_a, 0)
    ps = psiE[a]
    dNg = sumall([diff(mul(NLAPSE, diff(ps, i + 1)), i + 1) for i in I3])
    got = xTx(ric_lin(tscale(dK, 2)))
    ws = [w for w in (-2, -1, 0, 1, 2, 3) if not sub(got, sub(dNg, scale(ps, w)))]
    print(f"INFO K4 field datum T[s^{a}]: psidot_E[delta^w K] = div(N grad psi_E) - w psi_E holds for w in {ws}")
