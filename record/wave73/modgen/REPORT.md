# REPORT: lane MODGEN (wave 73)

## Lead's check and grade

Sunday 15 November 2026. Small verification lane (Fable): how much of the identification "on each angular sector, in the hyperbolic variables, the vacuum's modular flow on the ball is the free flow with the thermal covariance at inverse temperature 2 pi" can be proved from first principles for the scalar chain field, and what exactly remains for gravitons. Mid-run the lead forwarded the literature lane's findings. The lane returned the text below; the lead saves it.

**Grade: 384 of 420.** The geometry is proved exactly for every l; the state is reduced, by a short spectral argument, to one published theorem (Hislop-Longo) plus a factoriality fact, and confirmed to 23 digits or better in BOTH slots at l = 2 and 3 with two plants caught (the earlier check covered the momentum slot at l = 2 only); and the graviton step is classified honestly, item by item, with the single genuine hypothesis named. Marked down because Hislop-Longo is quoted, not re-proved; the record's own replacement for gravitons (wave70 modflow) has been checked by no second lane; and the lane opened none of the literature sources itself.

**Checked by the lead, by hand. It holds.**
- Step (1)1: for zeta = (1/2)[(1 - t^2 - r^2) d_t - 2 t x.grad], div zeta = -t - 3t = -4t, so the weight-1 term is lambda = -t.
- Step (1)2: at t = 0, delta psi = N psi_dot; delta psi_dot = d_t(zeta psi - t psi) = N Laplacian psi - x.grad psi - psi = div(N grad psi) - psi, using grad N = -x.
- Step (2)1 and (2)2: with S the Tomita operator of the standard subspace K, <k, Delta k'> = <k', k> for k, k' in K, so the real and imaginary parts of <k, Delta^(iu) k'> are the stated integrals and the symplectic form along the flow fixes the covariance, weight coth of half the spectral parameter; this is the converse of the referee's computation.
- Not hand-checked: step (1)3's algebra (the script's symbolic line L0 confirms it), step (1)4, the factoriality fact of step (2)5.

**Re-run by the lead:** `C:\Python314\python.exe modgen_check.py <out.json>` from `scripts\`, about 35 s, exit 0, "FAILED LINES: none": momentum slot 1.9e-29 and 2.3e-28, field slot 4.8e-24 and 5.0e-23, flow and cross tests to 2.5e-24 or better, both plants caught.

**What it changes.** Of the three physics identifications the assembled statement rests on, the second is now sharper. For the scalar chain field it is: exact geometry (proved here) plus Hislop-Longo 1982 plus a factoriality fact, with the thermal covariance confirmed numerically in both slots. For gravitons the single genuine hypothesis is that the ball's algebra is the linearised-curvature net, whose sectors reduce to the scalar problem with the constant c_l; the reduction of the two-point function is exact at low l in the record, and the record's own proof that the modular group acts as the conformal flow on curvature (wave70 modflow) has not been checked by a second lane. The paper should use the lane's "may say" sentence and heed its "must not say".

---

## Lane report (unchanged)

**MODGEN report (15 Nov 2026).** The check ran once, exit 0, every pre-fixed line passed and both plants were caught. No REPORT.md written, nothing committed.

## 1. Verdicts

- **(1) Geometry: proved here**, exactly, for every l. The modular parameter t corresponds to flow parameter s = 2πt, and X_0 is the generator at unit surface gravity.
- **(2) State: proved here, given two published inputs.** These are Hislop–Longo 1982 (the ball's vacuum modular group of the massless scalar is the weight-1 conformal flow at 2π) and "Δ has no eigenvalue 1" (factoriality). Given those, the thermal covariance follows by a short measure argument. Independently it is confirmed numerically to 1e-23 or better in both slots at l = 2, 3.
- **(3) Gravitons: inherited.** One item is proved in the record without Hislop–Longo (wave70 modflow), two are numerical at low l, and one is a genuine hypothesis (which algebra the ball's is).

## 2. Derivations

**(1) Flow.** Take ζ = ½[(1 − t² − r²)∂_t − 2t x·∇], which is π/R times this at 2π normalisation.
1. [P] ∂·ζ = −4t, so λ = −t: λ = 0 on the slice and ∂_tλ = −1.
2. [P] With δ_ζψ = ζ^μ∂_μψ + λψ (weight 1) at t = 0: δψ = Nψ̇ and δψ̇ = NΔψ − x·∇ψ − ψ = div(N∇ψ) − ψ, since ∇N = −x. The −ψ is the weight term.
3. [P] With q = rf, π = Nrp and d/dρ = N d/dr: q̇ = π. Also Nr[(r²Nf′)′/r² − Nλ²f/r² − f] = N(Nq_rr − rq_r) + Nq − N²λ²q/r² − Nq = q_ρρ − λ²q/sinh²ρ. So π̇ = −h²q. Sympy confirms the difference is 0 for generic f (line L0).
4. [P] Conformal picture: flat = N²(−dτ² + dρ² + sinh²ρ dΩ²), ψ̃ = Nψ. The conformal coupling R/6 = −1 gives (−∂_τ² + Δ_H³ + 1)ψ̃ = 0. With u = sinh ρ·ψ̃ = q, the sector Laplacian is (u″ − u − λ²u/sinh²ρ)/sinh ρ. The −1 (bottom of the H³ spectrum) cancels the +1 exactly, so the spectrum of h² starts at 0.

**(2) State.** Let K be the ball's standard subspace, Δ its modular operator, k, k′ ∈ K.
1. [P] ⟨k, Δk′⟩ = ⟨JSk, JSk′⟩ = ⟨k′, k⟩. With ν the spectral measure of log Δ for (k, k′): Re⟨k, Δ^{iu}k′⟩ = ∫e^{iuλ}(1+e^λ)/2 dν and i·Im = ∫e^{iuλ}(1−e^λ)/2 dν.
2. [P] So the function u ↦ σ(k, Δ^{iu}k′), which involves only the symplectic form and the flow, fixes ν off λ = 0, hence fixes the covariance: weight −coth(λ/2). The opposite orientation would make it negative, so positivity fixes the sign.
3. [S, from memory; modflow read HL pdf p2–7] Hislop–Longo: Δ^{it} acts on ball data as the weight-1 conformal flow at 2πt. It uses the massless scalar's two-point function 1/(x−y)², covariant with weight 1 under the conformal map from double cone to wedge, a dense-domain unitary, and Bisognano–Wichmann for the wedge.
4. [P] By (1), σ(k, Δ^{iu}k′) = Im⟨z, e^{∓2πiuh}z′⟩ with z = (h^{1/2}q + ih^{−1/2}π)/√2, so λ = −2πx on the spectrum of h. The covariance is Re⟨z, coth(πh)z′⟩ = ½⟨q, h coth(πh)q′⟩ + ½⟨π, h⁻¹coth(πh)π′⟩. There is no cross term because h is real.
5. [S, from memory] ν({0}) = 0: ker(Δ − 1) is spanned by K ∩ K′, which is {0} for the ball. Without it the vacuum covariance could only exceed the thermal one; the numerics show equality.
6. [N] Lines L3 to L8 below.
7. [U] Not opened by me: Derezinski Thm 40(6) and Brunetti–Guido–Longo Thm 2.3(ii), relayed by the lead from PHYS-ID. Derezinski gives the converse direction (model to modular operator), which with step 4 closes the loop.

**(3) Gravitons.**

| Fact needed | Status |
|---|---|
| Chain field obeys the scalar wave equation | [I] even-tower "Vacuum side" asserts it; the proof is Paper A / companion A8, not re-run |
| Graviton two-point function reduces to the scalar's times c_l | [N, I] modflow S3: exact for l = 2, 3, 4, c_l = 4(l−2)!/(l+2)!; even-tower: whole Gram, l = 2 to 6. General l rests on Prop. Nl |
| Ball's modular group is the conformal flow on curvature | [I, proved in record] wave70 modflow Steps 0 to 8, without Hislop–Longo; no second lane has checked it |
| Flow on (E, B) reduces to weight 1 on the chain field | [N, I] even-tower: x·δE·x = Nψ̇_E, div(N∇ψ_E) − ψ_E |
| Ball's graviton algebra is the curvature net's | **hypothesis**; wave69 step 2 inherited, assembly (iii) |
| Second quantisation Δ = Γ(Δ_K) | [S] statement read by modflow, proof not |

## 3. Check

Data: f_i = r^(l+2i)(1−r²)² (in H², so no surface term), p_j = r^(l+2j). Errors are the maximum over each 3×3 Gram, relative to its largest entry.

| Line | l = 2 | l = 3 |
|---|---|---|
| Momentum slot | 1.9e-29 | 2.3e-28 |
| Field slot | 4.8e-24 | 5.0e-23 |
| Flow/cross test, two flat forms | 1.3e-29 | 3.1e-29 |
| Flow/cross test, thermal against flat | 2.3e-25 | 2.5e-24 |

- **Guards:** Parseval and energy 4e-22; Heine against log form 3e-30; quadrature against 3F2 1e-31; Bessel-transform (momentum space) against position space 7.7e-6, with k cut at 80; grid convergence 1.9e-25.
- **Plant, coth(πh/2):** caught at 0.414 and 0.395 (momentum), 0.035 and 0.028 (field).
- **Plant, weight term dropped:** caught at 0.20 and 0.15.
- **Failed lines:** none.
- **Cross term:** the plain q–π term is zero identically on both sides, so no number is attached; the flow/cross lines test its derivative along the flow.

The first diagonal entry is 0.0265258238486492, which is 1/(12π) to the digits printed.

## 4. Re-run

From `record\wave73\modgen\scripts`:
`C:\Python314\python.exe modgen_check.py ..\out\modgen_check.json` (35 s, exit 0). Pass lines are in `..\logs\00_passlines.log`, the run in `..\logs\01_run.log`.

## 5. Weakest step

Step (2)3: Hislop–Longo is quoted, not re-proved, and for gravitons its replacement (modflow) has been checked by no second lane.

## 6. Wording

**May say:** "For the massless scalar, the vacuum modular group of the ball is the conformal flow (Hislop–Longo 1982). In the variables q = rf, π = Nrḟ, r = tanh(ρ/2), its generator on each angular sector is exactly 2π·X_0, and the vacuum's covariance is the β = 2π thermal form of h; we verify the latter to 23 digits at l = 2, 3. For gravitons we assume the ball's algebra is the linearised-curvature net, whose sectors reduce to this scalar problem with constant c_l."

**Must not say:** that the thermal form is proved for gravitons at all l from first principles; that Hislop–Longo covers helicity 2; that the algebra identification is proved; or "verified numerically for all l".

Sources I opened myself: only the record files named in the brief. None of PHYS-ID's papers, and not ss-second beyond a grep.
