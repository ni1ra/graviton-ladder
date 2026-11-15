# REPORT: lane EVEN-CONST (wave 73)

## Lead's check and grade

Sunday 15 November 2026. Small analytic lane (Fable): derive the even-parity per-sector constant in the large-l (Liouville) limit model, or say what kind of number it is; reproduce the odd constant first. The lane returned the text below; the lead saves it.

**Grade: 388 of 420.** It reproduced the odd constant to 17 digits before touching the even one, derived the even constant as a one-dimensional integral and then as a closed form, reported a null PSLQ result in the basis the brief gave before finding the term the basis lacked, and corrected the brief on what the constants are constants OF. Marked down because the agreement with the older ungraded lane hinge-e is not blind (the lane read it before finishing and says so), the second-order affinity formula of its step 2 is its own hand derivation, and that the lattice's finite-l quantity converges to this limit is inherited.

**The result.** In the limit model,

    a_even = 49 zeta(3)/(16 pi^2) - 3 beta(4)/pi^2 + G/8 - (ln 2)/4 = 0.013599405480333791...,

with beta(4) the Dirichlet beta value 0.988944551741... and G Catalan's constant, against a_odd = 7 zeta(3)/(16 pi^2) = 0.0532847998... The lattice value 0.013586 (range 0.013555 to 0.013625) agrees to 0.1%.

**Correction to the lead's brief, accepted.** These constants are not the per-sector constants of the squared Hilbert-Schmidt norm of the antilinear part. They are the constants of the affinity (fidelity) defect, the limit of lambda^2 (1 - A_l), which is what the lattices measured; the sum over sectors of (2l+1)(1 - A_l) grows like 2 a ln L. The first-order Hilbert-Schmidt norm at modular time 1/2 has different limit values (0.29166 odd, 0.05737 even). The record's proof of non-implementability goes through the Hilbert-Schmidt norm and does not use the affinity constants; the two are different measures of the same obstruction. Earlier wording in the record and in the public README that called 0.0533 and 0.0136 "the per-sector constant of the Shale sum" was loose and is corrected in the README.

**Checked by the lead.**
- Arithmetic of the closed form by hand: 0.372993 - 0.300603 + 0.114496 - 0.173287 = 0.013599.
- Independently of the lane's script, in 40-digit mpmath: the one-dimensional integral (1/64) int_0^inf tanh^2(x) Phi(x)/x^2 dx with the lane's Phi equals the closed form to all 20 digits printed (0.013599405480333791098), with beta(4) computed from Hurwitz zeta; and (1/32) int tanh^2(x)/x^2 dx equals 7 zeta(3)/(16 pi^2) to 18 digits. So steps 6 and 7 (integral to closed form) hold numerically at high precision.
- Not checked by the lead: steps 2 to 5 (from the definition of the affinity to the one-dimensional integral), the two cited formulas (DLMF 10.32.17, Gradshteyn-Ryzhik 6.671.14), completeness of the Kontorovich-Lebedev transform.

**Re-run by the lead:** `C:\Python314\python.exe -X utf8 even_const.py <out.json>` from `scripts\`, about 25 s, exit 0, "FAILED LINES: none"; the plant (wrong spectral measure) moves the odd value by 69% and is caught.

**What it changes.** The even constant has a closed form, and the odd constant's "Bessel identity" (an open item since HINGE) is reduced to one cited textbook formula. The paper may use the lane's "may say" sentence. It must not attach either constant to the Hilbert-Schmidt norm.

---

## Lane report (unchanged)

**EVEN-CONST report (2026-11-15)**

The even constant has a closed form, proved here in the limit model: **a_even = 49ζ(3)/(16π²) − 3β(4)/π² + G/8 − (ln 2)/4 = 0.013599405480333791…**, where β(4) is the Dirichlet beta value 0.98894455… and G is Catalan's constant. All pre-fixed lines passed on the first run. hinge-e already holds the same closed form, obtained by the same kind of reduction; I re-derived it by hand and then compared.

**Correction to the brief.** The constant is not ‖β̂₁‖²_HS at τ = π. That norm is 0.2916612 for odd (hinge's J(0.5)/2) and 0.0573715 for even. The constants 7ζ(3)/(16π²) and 0.0136 are the affinity defect, lim λ²(1 − A_l), which is what the lattices measured. So "2a ln L" is the growth of Σ(2l+1)(1 − A_l); the Shale sum grows as J(t) ln L.

## Derivation

Write P = p + p′, q = p − p′, a = πP/2, b = πq/2.

1. **Modes [P, given completeness of the Kontorovich–Lebedev transform].** f_p = √μ(p)·K_ip(e^(−y)) with μ(p) = (2/π²) p sinh πp. Plancherel was checked on x e^(−x) and x² e^(−x).
2. **Affinity at second order [P].** 1 − A = ½ Σ |W_mn|² (√ρ_m − √ρ_n)²/(E_m − E_n)². For W = ∫∫ M (±bb′ ± b†b′† + 2b†b′) in the 2π-thermal state:

   a = (1/2π²) ∫∫ dp dp′ M²/(sinh πp · sinh πp′) · [sinh²(πP/2)/P² + sinh²(πq/2)/q²]

   This is a quarter of the β-def norm with cosh(πν) sin²(πν/2) replaced by sinh²(πν/2).
3. **Odd kernel [P, given DLMF 10.32.17].** M = π v/(2√(pp′)), v = √(μμ′)·∫₀^∞ K_ip K_ip′ dx. Macdonald's product formula, with ∫K_σ dx = π/(2 cos(πσ/2)) and ∫cosh(κt)/cosh t dt = π/(2 cos(πκ/2)), gives ∫K_ip K_ip′ dx = π²/(4 cosh a · cosh b). This is the identity that does the work, and it closes the record's open item at the level of one cited textbook formula.
4. **Even kernel [P, given GR 6.671.14].** Using G = (2/π)∫ sin kx · sin kx′/(1 + k²) dk and ∫K_ip sin(kx) dx/x = π sin(p·arsinh k)/(2p cosh(πp/2)), M_even = M_odd · T with T = tanh(πp/2) tanh(πp′/2).
5. **Double integral [P].**

   a = (1/64) ∫₀^∞ da ∫₋ₐ^ₐ db T^(2e) [tanh²a/(a² cosh²b) + tanh²b/(b² cosh²a)], with e = 0 for odd and 1 for even.

6. **Odd [P].** The inner integrals give (1/32)∫tanh²x/x² dx; the tanh³ terms cancel. With ∫tanh²x/x² dx = 14ζ(3)/π² (Mellin continuation of sech², using ζ′(−2)), a_odd = 7ζ(3)/(16π²).
7. **Even [P].** With T = (cosh a − cosh b)/(cosh a + cosh b), partial fractions give a_even = (1/64)∫tanh²x · Φ(x)/x² dx, where

   Φ(x) = 2 − 4π sech x + 8x/(sinh x cosh x) + 8(x cosh x − sinh x)/sinh³x.

   The Gudermannian and log cosh terms cancel. The four integrals needed are:
   - ∫tanh²x/x² dx = 14ζ(3)/π²
   - ∫tanh x · sech²x/x dx = 7ζ(3)/π² (by parts)
   - ∫(1/(x sinh x) − 1/x²) dx = −ln 2
   - ∫tanh²x · sech x/x² dx = 48β(4)/π³ − 2G/π, from tanh² sech = (sech + sech″)/2, whose Mellin transform is Γ(s)[β(s) + β(s−2)], evaluated at s = −1 with the functional equation.

## Results

- **(1) Odd.** Both quadratures match 7ζ(3)/(16π²) to relative 2.5e-17 and 4e-16 [N]. The identity is now [P] given DLMF 10.32.17, and also checked numerically to 2.6e-26.
- **(2) Even.** The two quadrature settings differ by 3.9e-18 (relative 2.9e-16) [N]. The 1-d reduction equals the closed form to 7.6e-123 [N]. The closed form itself is [P] in the limit model.
- **PSLQ (100 digits, coefficients up to 10⁴).** In the brief's basis {a, 1, ζ3, ζ5, π², ln 2, G, 1/π², ζ3/π², ζ5/π², ln 2/π², G/π²} the result is null. Adding β(4)/π² returns exactly 16a = 49ζ3/π² − 48β4/π² + 2G − 4 ln 2. So the brief's basis could not have found it.
- **(3) Lattice.** a_even is inside [0.013555, 0.013625] and sits 1.34e-5 (0.099%) above 0.013586.
- **(4) Check.** L1–L12, I1 and I2 all passed; no failed lines. The plant (μ = p tanh πp) gives odd = 0.016679, a 68.7% deviation, caught. The finite-difference solve of −w″ + w = K_ip′/x reproduces w(p, p′) to 6e-12.
- **(5) Re-run.** From `record\wave73\even-const\scripts`: `C:\Python314\python.exe -X utf8 even_const.py ..\out\even_const.json` (25 s, exit 0). Logs are `logs\01_prefixed_pass_lines.log` and `logs\02_run1.log`.
- **(6) Weakest step.** The second-order affinity formula of step 2 and the thermofield weights are my hand derivation; they are confirmed only by reproducing the odd constant and by agreeing with hinge-e, which I read before finishing, so the agreement is not blind. That the lattice's finite-l quantity converges to this limit is [I].
- **(7) May say:** "In the Liouville limit model the even-parity per-sector affinity constant is 49ζ(3)/(16π²) − 3β(4)/π² + G/8 − (ln 2)/4 ≈ 0.0135994, against 7ζ(3)/(16π²) for odd parity, and the lattice value 0.013586 agrees to 0.1%."
- **Must not say:** "The even Shale norm ‖β‖²_HS has per-sector constant 0.0136", or that the even divergence is thereby a theorem at finite l.
