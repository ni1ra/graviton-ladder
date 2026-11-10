"""S3: spanning lemma and the conformal scalar's stress tensor on data (phi, pi) = (0, g).
Usage: python3 s3_span.py OUT.json"""
import json, sys
import sympy as sp

t, x, y, z = sp.symbols('t x y z', real=True)
C = (t, x, y, z)
eta = sp.diag(-1, 1, 1, 1)
xi_c = sp.Rational(1, 6)        # conformal coupling in d = 4
phi = sp.Function('phi')(*C)
out = {}

# improved stress tensor T_mn = d_m phi d_n phi - 1/2 eta_mn (d phi)^2 + xi (eta_mn box - d_m d_n) phi^2  (flat)
dphi = [sp.diff(phi, c) for c in C]
kin = sum(eta[m, m]*dphi[m]**2 for m in range(4))
box = lambda f: sum(eta[m, m]*sp.diff(f, C[m], 2) for m in range(4))
T = [[dphi[m]*dphi[n] - sp.Rational(1, 2)*eta[m, n]*kin + xi_c*(eta[m, n]*box(phi**2) - sp.diff(phi**2, C[m], C[n]))
      for n in range(4)] for m in range(4)]
# trace check (conformal: trace vanishes on shell box phi = 0)
trace = sp.expand(sum(eta[m, m]*T[m][m] for m in range(4)))
onshell = sp.expand(trace.subs(sp.Derivative(phi, (t, 2)), sp.diff(phi, x, 2) + sp.diff(phi, y, 2) + sp.diff(phi, z, 2)))
out["traceless_on_shell"] = bool(sp.simplify(onshell) == 0)

# data on Sigma: phi(0, x) = 0 identically, d_t phi(0, x) = g(x); then all spatial derivatives of phi vanish on Sigma
g = sp.Function('g')(x, y, z)
pt = sp.Symbol('PT')     # stands for d_t phi on Sigma
def on_sigma(e):
    e = sp.expand(e.doit())
    # replace derivatives: spatial-only derivatives of phi -> 0; d_t phi -> g; d_t d_spatial phi -> d_spatial g
    reps = {}
    for D in e.atoms(sp.Derivative):
        vs = D.variables
        nt = sum(1 for v in vs if v == t)
        sp_vars = [v for v in vs if v != t]
        if nt == 0:
            reps[D] = 0
        elif nt == 1:
            reps[D] = sp.diff(g, *sp_vars) if sp_vars else g
        # nt >= 2 terms: keep symbolic (should cancel or multiply phi = 0)
    e = e.subs(reps)
    e = e.subs(phi, 0)
    return sp.simplify(e)
T00 = on_sigma(T[0][0])
T0i = [on_sigma(T[0][i]) for i in (1, 2, 3)]
out["T00_on_data_(0,g)"] = str(T00)
out["T00_equals_half_g2"] = bool(sp.simplify(T00 - g**2/2) == 0)
out["T0i_on_data_(0,g)"] = [str(e) for e in T0i]
out["T0i_zero"] = all(sp.simplify(e) == 0 for e in T0i)
out["Wald_term_phi2_on_S"] = "phi = 0 on Sigma, so the conformal-coupling Wald term -2 pi xi \\oint phi^2 vanishes"

# spanning: for any smooth rho on Bbar and chi = 1 on Bbar,  1/2 (chi + rho)^2 - 1/2 (chi - rho)^2 = 2 chi rho = 2 rho on Bbar
chi, rho = sp.symbols('chi rho')
out["span_identity"] = bool(sp.expand(sp.Rational(1, 2)*(chi + rho)**2 - sp.Rational(1, 2)*(chi - rho)**2 - 2*chi*rho) == 0)

out["ALL_PASS"] = out["traceless_on_shell"] and out["T00_equals_half_g2"] and out["T0i_zero"] and out["span_identity"]
s = json.dumps(out, indent=1); print(s)
open(sys.argv[1], "w").write(s + "\n")
