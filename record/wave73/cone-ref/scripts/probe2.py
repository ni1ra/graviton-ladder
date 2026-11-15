"""CONE-REF probe, part 2 (exact, sympy; a and R rational to keep it fast).
2. finite Moebius covariance: w = v/(a + b v), b = (1-a)/R (fixes tip and edge), Omega = w/v, dw/dv = a Omega^2;
   curvature of dimension 3: C'(w) = a^-2 Omega^-5 C(v).  Is Pi_p(w) dw = Pi_p(v) dv, Pi_p = 2 v sigma_p ?
3. first-law bracket by parts.   4. RULER-L identity with cross terms."""
import sympy as sp
v, u, w = sp.symbols("v u w", positive=True)
out = []
def say(k, val): out.append(f"{k}: {val}"); print(out[-1], flush=True)
def PiC(Cexpr, p, var):
    return -2*var**(1-p) * sp.integrate(sp.expand((u**p * Cexpr).subs(var, u)), (u, 0, var))
for (a, R) in ((sp.Rational(3, 5), sp.Integer(2)), (sp.Rational(7, 3), sp.Integer(1))):
    b = (1 - a)/R
    C = v**2 + 3*v**5
    vw = a*w/(1 - b*w); Om = w/vw
    Cp = sp.cancel((C / a**2).subs(v, vw) * Om**(-5))
    for p in (1, 2, 3):
        lhs = PiC(Cp, p, w)
        rhs = (PiC(C, p, v) / sp.diff(v/(a + b*v), v)).subs(v, vw)
        d = sp.simplify(lhs - rhs)
        say(f"weight-one covariance a={a} R={R} p={p}", d == 0 if d == 0 else f"False (at w=1/3: {sp.nsimplify(d.subs(w, sp.Rational(1,3)))})")
R = sp.Symbol("R", positive=True); n = sp.Symbol("n", positive=True)
dA = -sp.integrate(u*v*(v-u)*u**n, (u, 0, v))
br = dA.subs(v, R) + sp.integrate(sp.expand((2/v - 6/R)*dA), (v, 0, R))
say("IBP: bracket + int v^2 f v^n (symbolic n, R)", sp.simplify(br + sp.integrate(sp.expand(v**3*(1-v/R)*v**n), (v, 0, R))))
say("bracket on dA = v^2 (affine rescaling)", sp.simplify(R**2 + sp.integrate((2/v - 6/R)*v**2, (v, 0, R))))
c1, c2 = sp.symbols("c1 c2", real=True)
Cg = c1*v + c2*v**4
s2 = -v**-2*sp.integrate((u**2*Cg).subs(v, u), (u, 0, v)); s3 = -v**-3*sp.integrate((u**3*Cg).subs(v, u), (u, 0, v))
f = v*(1 - v/R)
say("RULER-L: int v^2 f (s2^2 - s3^2) - int v^3 (s2-s3)^2", sp.simplify(sp.integrate(sp.expand(v**2*f*(s2**2 - s3**2) - v**3*(s2-s3)**2), (v, 0, R))))
open("logs/probe2.log", "w", encoding="utf-8").write("\n".join(out) + "\n")
