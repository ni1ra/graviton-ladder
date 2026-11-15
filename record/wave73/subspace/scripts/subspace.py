"""Lane SUBSPACE (wave 73): does the odd-parity transfer (h, K) -> (psi, psi_dot) carry the ball's phase space onto
the scalar one, intertwine the vacuum two-point functions, and turn the conformal flow into weight 1?
Exact rational arithmetic.  Usage: python subspace.py L [KMAX]   (L >= 2; KMAX default 3)
Definitions and pass/fail lines: ../PREREG.md.
polyalg.py is imported read-only from charge-id; KP, fg_from_eta, build_T, psi_of_K, psidot_of_h are copied from
charge-id/scripts/charge_id.py (that script runs at import, so it cannot be imported)."""
import sys, os, json, time
from fractions import Fraction as Fr
from math import factorial

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..", "charge-id", "scripts")))
import polyalg as P
from polyalg import (add, sub, mul, scale, const, pdiff, diff, sumall, ball, gauss2, dot, trace, div, lap,
                     tlap, curl, xTx, ric_lin, tadd, tscale, tmulp, tzero, I3, XS, R2, NLAPSE, eps, harmonic,
                     powp, dfact)

L = int(sys.argv[1])
KMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 3
LL = L * (L + 1)
CL = Fr(4 * factorial(L - 2), factorial(L + 2))
M_ = L + 2
OUT = {"l": L, "kmax": KMAX, "checks": {}}
FAILS = []
PLANTS_MISSED = []
T0 = time.time()


def log(*a):
    print(*a, flush=True)


def check(name, ok, info=""):
    OUT["checks"][name] = {"pass": bool(ok), "info": str(info)}
    log(("PASS " if ok else "FAIL ") + name + ("  " + str(info) if info != "" else ""))
    if not ok:
        FAILS.append(name)


def plant(name, caught, info=""):
    OUT["checks"][name] = {"pass": bool(caught), "plant": True, "info": str(info)}
    log(("CAUGHT " if caught else "MISSED ") + name + ("  " + str(info) if info != "" else ""))
    if not caught:
        PLANTS_MISSED.append(name)


def stop(msg):
    log("STOP: " + msg)
    OUT["fails"] = FAILS
    OUT["stopped"] = msg
    json.dump(OUT, open(os.path.join(HERE, "..", "out", f"subspace_l{L}.json"), "w"), indent=1)
    sys.exit(2)


# ---------------------------------------------------------------- copied helpers (charge-id)
def spoly(coefs):
    return sumall([scale(powp(R2, k), c) for k, c in enumerate(coefs)])


def dcoefs(coefs):
    return [k * c for k, c in enumerate(coefs)][1:] or [Fr(0)]


def padd(a, b):
    n = max(len(a), len(b))
    a = list(a) + [Fr(0)] * (n - len(a))
    b = list(b) + [Fr(0)] * (n - len(b))
    return [x + y for x, y in zip(a, b)]


def KP(phi):
    g = [pdiff(phi, i + 1) for i in I3]
    v = [sumall([scale(mul(XS[j + 1], g[k]), eps(i, j, k)) for j in I3 for k in I3 if eps(i, j, k)])
         for i in I3]
    K = [[scale(add(pdiff(v[j], i + 1), pdiff(v[i], j + 1)), Fr(1, 2)) for j in I3] for i in I3]
    Pt = [[add(mul(XS[i + 1], v[j]), mul(v[i], XS[j + 1])) for j in I3] for i in I3]
    return v, K, Pt


def fg_from_eta(coefs, m=M_):
    coefs = [Fr(c) for c in coefs]
    d = dcoefs(coefs)
    if P.GAUSS[0]:
        d = padd(d, [-c for c in coefs])
    g = [Fr(-2, m) * c for c in d]
    sg = [Fr(0)] + g
    f = [Fr(2, L - 1) * c for c in padd(coefs, [-c for c in sg])]
    return f, g


def build_T(phi, coefs, m=M_):
    v, K, Pt = KP(phi)
    f, g = fg_from_eta(coefs, m)
    return tadd(tmulp(K, spoly(f)), tmulp(Pt, spoly(g)))


def psi_of_K(K):
    return scale(xTx(curl(K)), -1)


def psidot_of_h(h):
    return scale(xTx(curl(tscale(ric_lin(h), -1))), -1)


def grad(p):
    return [diff(p, i + 1) for i in I3]


def is_zero_vec(w):
    return all(not c for c in w)


def is_zero_t(A):
    return all(not A[i][j] for i in I3 for j in I3)


# ---------------------------------------------------------------- new helpers
def sym(A):
    return [[scale(add(A[i][j], A[j][i]), Fr(1, 2)) for j in I3] for i in I3]


def symcurl(A):
    return sym(curl(A))


def symxcross(Mt):
    """sym of (x cross M)_ij = eps_ikl x_k M_lj"""
    C = tzero()
    for i in I3:
        for j in I3:
            acc = {}
            for k in I3:
                for l in I3:
                    e = eps(i, k, l)
                    if e:
                        acc = add(acc, scale(mul(XS[k + 1], Mt[l][j]), e))
            C[i][j] = acc
    return sym(C)


def gdot(pa, pb):
    return sumall([mul(x, y) for x, y in zip(grad(pa), grad(pb))])


def unit(k):
    return [Fr(0)] * k + [Fr(1)]


# Fourier transform of p(x) e^{-r^2}:  pi^{3/2} i^{n0} Q(k) e^{-k^2/4},  x_j -> i d/dk_j.
_H1 = [[Fr(1)]]


def h1(n):
    """1D polynomial h_n with (d/dk)^n e^{-k^2/4} = h_n(k) e^{-k^2/4}: h_{n+1} = h_n' - k h_n / 2."""
    while len(_H1) <= n:
        h = _H1[-1]
        new = [Fr(0)] * (len(h) + 1)
        for i, c in enumerate(h):
            if i:
                new[i - 1] += i * c
            new[i + 1] -= c / 2
        _H1.append(new)
    return _H1[n]


def ft(p):
    """returns (Q, n0); p must be t-free with all monomials of one parity n0."""
    if not p:
        return {}, None
    par = {(a + b + c) % 2 for (t, a, b, c) in p}
    assert len(par) == 1, "mixed parity"
    n0 = par.pop()
    out = {}
    for (t, a, b, c), v in p.items():
        assert t == 0
        sg = (-1) ** ((a + b + c - n0) // 2)
        ha, hb, hc = h1(a), h1(b), h1(c)
        for i, ci in enumerate(ha):
            if not ci:
                continue
            for j, cj in enumerate(hb):
                if not cj:
                    continue
                for k, ck in enumerate(hc):
                    if not ck:
                        continue
                    key = (0, i, j, k)
                    out[key] = out.get(key, 0) + v * sg * ci * cj * ck
    return {k: v for k, v in out.items() if v != 0}, n0


def ftt(A):
    Q = [[None] * 3 for _ in I3]
    pars = set()
    for i in I3:
        for j in I3:
            q, n0 = ft(A[i][j])
            Q[i][j] = q
            if n0 is not None:
                pars.add(n0)
    assert len(pars) <= 1
    return Q


def Iw(Q, s):
    """int d^3k |k|^s Q(k) e^{-k^2/2}; units 4 pi for s = +-1, 4 pi sqrt(pi/2) for s = 0."""
    tot = Fr(0)
    for (t, a, b, c), v in Q.items():
        if a % 2 or b % 2 or c % 2:
            continue
        D = a + b + c
        ang = Fr(dfact(a - 1) * dfact(b - 1) * dfact(c - 1), dfact(D + 1))
        m = D + 2 + s
        rad = dfact(m - 1) if s == 0 else 2 ** ((m - 1) // 2) * factorial((m - 1) // 2)
        tot += v * ang * rad
    return tot


def rank(rows):
    rows = [list(r) for r in rows]
    rk = 0
    ncol = len(rows[0]) if rows else 0
    for c in range(ncol):
        piv = next((i for i in range(rk, len(rows)) if rows[i][c] != 0), None)
        if piv is None:
            continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        pv = rows[rk][c]
        for i in range(len(rows)):
            if i != rk and rows[i][c] != 0:
                fct = rows[i][c] / pv
                rows[i] = [x - fct * y for x, y in zip(rows[i], rows[rk])]
        rk += 1
    return rk


nb = KMAX + 1
phi = harmonic(L, Fr(3, 5), Fr(4, 5))
v1, K1, P1 = KP(phi)
log(f"l={L} c_l={CL} KMAX={KMAX}")

# ================================================================ A: localisation (unit ball, plain polynomials)
okA1 = True
for k in range(nb):
    T = build_T(phi, unit(k))
    okA1 &= is_zero_vec(div(T)) and not trace(T)
    okA1 &= not sub(psi_of_K(T), scale(mul(powp(R2, k), phi), LL))
check("A1 psi[build_T(s^k)] = l(l+1) s^k phi exactly, build_T is TT (K <-> psi local both ways)", okA1)

# plant: wrong constant in the local inverse must break transversality
Tbad = build_T(phi, unit(1), m=L + 1)
plant("PA1 inverse with (l+1) in place of (l+2) is not TT", not is_zero_vec(div(Tbad)))

DEG = 3
basis = [tmulp(K1, powp(R2, i)) for i in range(DEG + 1)] + [tmulp(P1, powp(R2, i)) for i in range(DEG + 1)]
divs = [div(b) for b in basis]
keys = sorted({(c, k) for d in divs for c in I3 for k in d[c]})
mat = [[d[c].get(k, Fr(0)) for d in divs] for (c, k) in keys]      # rows = equations, cols = unknowns
rk = rank(mat)
okA2 = (len(basis) - rk == DEG + 1) and all(not trace(b) for b in basis)
check("A2 TT odd tensors a(s)K_l + b(s)P with deg <= 3 form a space of dimension exactly 4 = span build_T(s^0..s^3)",
      okA2, f"unknowns={len(basis)} rank={rk} nullity={len(basis)-rk}")

kap = []
propto = True
for k in range(nb):
    pd = psidot_of_h(tmulp(P1, powp(R2, k)))
    ref = mul(powp(R2, k), phi)
    key = next(iter(ref))
    c = pd.get(key, Fr(0)) / ref[key]
    propto &= not sub(pd, scale(ref, c))
    kap.append(c)
okA3 = propto and len(set(kap)) == 1 and kap[0] != 0
check("A3 psi_dot[alpha(r) P] = kappa_l alpha phi with one constant (Regge-Wheeler h <-> psi_dot algebraic)", okA3,
      f"proportional={propto} kappa_k={[str(x) for x in kap]}")
OUT["kappa_RW"] = [str(x) for x in kap]

c13 = [Fr(1), Fr(-3), Fr(3), Fr(-1)]


def eta_c(a):
    return [Fr(0)] * a + c13                 # (1-s)^3 s^a


def rw_of(eta):
    f, g = fg_from_eta(eta)
    return tmulp(P1, spoly(padd(g, [-x for x in dcoefs(f)])))


okA4 = True
okA4rw = True
nz = False
for a in range(3):
    Ka = build_T(phi, eta_c(a))
    pa = psi_of_K(Ka)
    for b in range(3):
        hb = build_T(phi, eta_c(b))
        rhs = -CL * ball(mul(pa, psidot_of_h(hb)))
        lhs = ball(dot(Ka, hb))
        lhs_rw = ball(dot(Ka, rw_of(eta_c(b))))
        okA4 &= lhs == rhs
        okA4rw &= lhs_rw == rhs
        nz |= lhs != 0
check("A4a int_B K_a.h_b = -c_l int_B psi_a psidot_b, data vanishing to 2nd order at the boundary, h TT", okA4 and nz)
check("A4b the same with h in Regge-Wheeler gauge", okA4rw and nz)
mism = []
for a in range(3):
    Ka = build_T(phi, unit(a))
    pa = psi_of_K(Ka)
    for b in range(3):
        hb = build_T(phi, unit(b))
        rhs = -CL * ball(mul(pa, psidot_of_h(hb)))
        mism.append((str(ball(dot(Ka, hb)) / rhs) if rhs != 0 else "rhs=0", str(ball(dot(Ka, rw_of(unit(b)))) / rhs) if rhs != 0 else "rhs=0"))
OUT["A4_info_boundary_touching_ratio_TT_and_RW"] = mism
log("INFO A4 boundary-touching data (eta = s^a): Omega_B/(c_l Omega_s,B) in (TT gauge, RW gauge): " + str(mism))

# ---- A5 (post hoc, added after the first run showed ratio 1 for every boundary-touching pair in RW gauge;
# not preregistered): the symplectic densities agree on every sphere, for arbitrary profiles.
def shell(p):
    """angular average of a t-free polynomial as {degree: coefficient of r^degree}, in units of 4 pi."""
    parts = {}
    for k, v in p.items():
        parts.setdefault(k[1] + k[2] + k[3], {})[k] = v
    out = {D: P.sph(q) for D, q in parts.items()}
    return {D: c for D, c in out.items() if c != 0}


okA5 = True
nzA5 = False
for a in range(nb):
    Ka = build_T(phi, unit(a))
    pa = psi_of_K(Ka)
    for b in range(nb):
        hb = tmulp(P1, powp(R2, b))                       # arbitrary Regge-Wheeler profile alpha = s^b
        lhs = shell(dot(Ka, hb))
        rhs = {D: -CL * c for D, c in shell(mul(pa, psidot_of_h(hb))).items()}
        okA5 &= lhs == rhs
        nzA5 |= bool(lhs)
check("A5 (post hoc) sphere by sphere: <K_a . alpha_b P>_S = -c_l <psi_a psidot_b>_S as polynomials in r (RW gauge)",
      okA5 and nzA5)

# ================================================================ C1: Paper A's flow (1) -> weight 1 (unit ball)
eK = [Fr(2), Fr(-1, 3), Fr(5, 7), Fr(1, 2)][:nb]
eh = [Fr(-1), Fr(3, 4), Fr(2, 5), Fr(-7, 3)][:nb]
Kd = build_T(phi, eK)
hd = build_T(phi, eh)
Bt = tscale(symcurl(Kd), -1)
Et = ric_lin(hd)
Edot = tscale(symcurl(Bt), -1)
Bdot = symcurl(Et)
check("C1a dE/dt = -sym curl B equals Ric_lin[2K] on TT data", is_zero_t(tadd(Edot, tscale(ric_lin(tscale(Kd, 2)), -1))))
dE = tadd(tmulp(Edot, NLAPSE), tscale(symxcross(Bt), 2))
dB = tadd(tmulp(Bdot, NLAPSE), tscale(symxcross(Et), -2))
psi = xTx(Bt)
psid = xTx(Bdot)
check("C1a' psi, psi_dot here equal charge-id's psi_of_K, psidot_of_h",
      not sub(psi, psi_of_K(Kd)) and not sub(psid, psidot_of_h(hd)) and bool(psi) and bool(psid))
check("C1b x.dB.x = N psi_dot", not sub(xTx(dB), mul(NLAPSE, psid)))
dpsid = xTx(symcurl(dE))
divNgrad = sumall([diff(mul(NLAPSE, diff(psi, i + 1)), i + 1) for i in I3])


def c1c(w):
    return not sub(dpsid, sub(divNgrad, scale(psi, w)))


check("C1c x.(sym curl dE).x = div(N grad psi) - 1*psi  (weight 1)", c1c(1))
plant("PC1 weight 2 in C1c", not c1c(2))
check("C1d dE, dB are transverse and traceless", is_zero_vec(div(dE)) and is_zero_vec(div(dB)) and not trace(dE) and not trace(dB))

# ================================================================ B: vacuum two-point function on R^3 (Gaussian mode)
P.GAUSS[0] = True
Tg = [build_T(phi, unit(k)) for k in range(nb)]
okTT = all(is_zero_vec(div(t)) and not trace(t) for t in Tg)
check("B-TT Gaussian-profile tensors are TT on R^3", okTT)
ps = [psi_of_K(t) for t in Tg]
pd = [psidot_of_h(t) for t in Tg]
gps = [[gdot(ps[a], ps[b]) for b in range(nb)] for a in range(nb)]
FT_T = [ftt(t) for t in Tg]
FT_ps = [ft(p)[0] for p in ps]
FT_pd = [ft(p)[0] for p in pd]


def kdot(A, B):
    return sumall([mul(A[i][j], B[i][j]) for i in I3 for j in I3])


okB0 = True
for a in range(nb):
    for b in range(nb):
        okB0 &= gauss2(dot(Tg[a], Tg[b])) == Iw(kdot(FT_T[a], FT_T[b]), 0)
        okB0 &= gauss2(mul(ps[a], ps[b])) == Iw(mul(FT_ps[a], FT_ps[b]), 0)
        okB0 &= gauss2(mul(pd[a], pd[b])) == Iw(mul(FT_pd[a], FT_pd[b]), 0)
check("B0 control: Parseval, position space = momentum space, all tensors and scalars", okB0)
if not okB0:
    stop("Parseval control failed: Fourier transform misread")

kset = set()
okE = True
for a in range(nb):
    for b in range(nb):
        og = gauss2(dot(Tg[a], Tg[b]))
        os_ = -gauss2(mul(ps[a], pd[b]))
        kset.add(og / (CL * os_) if os_ != 0 else None)
        okE &= og == CL * Fr(1, 2) * gauss2(gps[a][b])
check("B1a control: Omega_grav = kappa c_l Omega_s with kappa = 1 (charge-id T4)", kset == {Fr(1)}, f"kappa={sorted(str(x) for x in kset)}")
check("B1b control: int K_a.K_b = c_l (1/2) int grad psi_a.grad psi_b (Paper A Prop. 1)", okE)
if kset != {Fr(1)} or not okE:
    stop("known-answer control failed: a definition is misread")

GK = [[2 * Iw(kdot(FT_T[a], FT_T[b]), -1) for b in range(nb)] for a in range(nb)]
SK = [[Iw(mul(FT_ps[a], FT_ps[b]), 1) for b in range(nb)] for a in range(nb)]
Gh = [[Fr(1, 2) * Iw(kdot(FT_T[a], FT_T[b]), 1) for b in range(nb)] for a in range(nb)]
Sh = [[Iw(mul(FT_pd[a], FT_pd[b]), -1) for b in range(nb)] for a in range(nb)]
SKm = [[Iw(mul(FT_ps[a], FT_ps[b]), -1) for b in range(nb)] for a in range(nb)]
Shp = [[Iw(mul(FT_pd[a], FT_pd[b]), 1) for b in range(nb)] for a in range(nb)]


def eqm(A, B, c):
    return all(A[a][b] == c * B[a][b] for a in range(nb) for b in range(nb))


nzB = all(any(x != 0 for r in Mx for x in r) for Mx in (GK, SK, Gh, Sh))
check("B2a K-K block: int 2 K_a.K_b/omega = c_l int omega psi_a psi_b, whole Gram", eqm(GK, SK, CL) and nzB,
      f"GK[0][0]={GK[0][0]} c_l*SK[0][0]={CL*SK[0][0]}")
check("B2b h-h block: int (1/2) omega h_a.h_b = c_l int psidot_a psidot_b/omega, whole Gram", eqm(Gh, Sh, CL) and nzB,
      f"Gh[0][0]={Gh[0][0]} c_l*Sh[0][0]={CL*Sh[0][0]}")
OUT["B2_ratios_diag"] = {"KK": [str(GK[a][a] / SK[a][a]) for a in range(nb)], "hh": [str(Gh[a][a] / Sh[a][a]) for a in range(nb)]}
cbad = CL * Fr(L + 3, L + 2)
plant("PB1 c_l scaled by (l+3)/(l+2)", not eqm(GK, SK, cbad) and not eqm(Gh, Sh, cbad))
# (h, K) -> (2h, K/2): scalar side becomes psi/2, psi_dot*2. Omega identity unchanged, two-point function not.
om_ok = all(gauss2(dot(Tg[a], Tg[b])) == CL * (-gauss2(mul(scale(ps[a], Fr(1, 2)), scale(pd[b], 2))))
            for a in range(nb) for b in range(nb))
plant("PB2 symplectic but not vacuum-preserving map (h,K)->(2h,K/2): Omega identity still holds, B2 fails",
      om_ok and not eqm(GK, SK, CL * Fr(1, 4)) and not eqm(Gh, Sh, CL * 4))
plant("PB3 omega^{+1} and omega^{-1} exchanged on the scalar side", not eqm(GK, SKm, CL) and not eqm(Gh, Shp, CL))

# ================================================================ C3: which weight preserves the scalar vacuum
F = [mul(powp(R2, k), phi) for k in range(nb)]          # prefactors of scalar Gaussian data in sector l
FTF = [ft(p)[0] for p in F]


def dPi(p, w):
    return sub(sumall([diff(mul(NLAPSE, diff(p, i + 1)), i + 1) for i in I3]), scale(p, w))


def xgrad(p):
    return sumall([mul(XS[i + 1], diff(p, i + 1)) for i in I3])


# controls: time translation and dilation are skew for mu_s
okT = all(Iw(mul(ft(lap(F[a]))[0], FTF[b]), -1) + Iw(mul(FTF[a], FTF[b]), 1) == 0 for a in range(nb) for b in range(nb))
check("C3-ctl1 control: time translation (Pi, Lap Phi) is skew for mu_s", okT)
okD = True
for a in range(nb):
    for b in range(nb):
        Da1, Db1 = add(xgrad(F[a]), F[a]), add(xgrad(F[b]), F[b])
        Da2, Db2 = add(xgrad(F[a]), scale(F[a], 2)), add(xgrad(F[b]), scale(F[b], 2))
        okD &= Iw(mul(ft(Da1)[0], FTF[b]), 1) + Iw(mul(FTF[a], ft(Db1)[0]), 1) == 0
        okD &= Iw(mul(ft(Da2)[0], FTF[b]), -1) + Iw(mul(FTF[a], ft(Db2)[0]), -1) == 0
check("C3-ctl2 control: dilation (x.grad Phi + Phi, x.grad Pi + 2 Pi) is skew for mu_s", okD)
okS = all(gauss2(sub(mul(dPi(F[a], w), F[b]), mul(dPi(F[b], w), F[a]))) == 0 for a in range(nb) for b in range(nb) for w in (1, 2))
check("C3-ctl3 control: delta_w is skew for Omega_s at w = 1 and w = 2 (both flows are symplectic)", okS)
if not (okT and okD and okS):
    stop("C3 control failed: a definition is misread")


def skew(w):
    return [[Iw(mul(ft(dPi(F[a], w))[0], FTF[b]), -1) + Iw(mul(FTF[a], ft(mul(NLAPSE, F[b]))[0]), 1)
             for b in range(nb)] for a in range(nb)]


s1, s2 = skew(1), skew(2)
check("C3 delta_1 = (N Pi, div(N grad Phi) - Phi) is skew for mu_s: the weight-1 flow preserves the vacuum",
      all(x == 0 for r in s1 for x in r))
plant("PC3 weight 2 is not skew for mu_s (the area-side flow does not preserve the vacuum two-point function)",
      any(x != 0 for r in s2 for x in r), f"skew_2[0][0]={s2[0][0]}")
P.GAUSS[0] = False

OUT["fails"] = FAILS
OUT["plants_missed"] = PLANTS_MISSED
OUT["seconds"] = time.time() - T0
os.makedirs(os.path.join(HERE, "..", "out"), exist_ok=True)
json.dump(OUT, open(os.path.join(HERE, "..", "out", f"subspace_l{L}.json"), "w"), indent=1)
nchk = len(OUT["checks"])
log(f"DONE l={L}: {nchk-len(FAILS)-len(PLANTS_MISSED)}/{nchk} lines as expected, fails={FAILS}, plants missed={PLANTS_MISSED}, t={time.time()-T0:.1f}s")
