# REPORT: lane EVEN-REMAINDER (wave 73)

## Lead's check and grade

Sunday 15 November 2026. Small proof lane (Fable): bound the second-order remainder of the even-parity flow on EVEN-PASSAGE's explicit family, uniformly in l, so that the first-order floor survives. The lane returned the text below; the lead saves it.

**Grade: 374 of 420.** All three parts (the cocycle remainder, the numerator constants, a lower bound on the vector's norm) are closed with explicit constants, the lane ran a pre-predicted diagnostic on its own failed line, and its list of what is still assumed is complete and unflattering. Marked down for a misprinted factor in its headline inequality and in step 15 (below), and because the proved bound only bites at very large l, so the table tests it non-vacuously at a single l.

**Correction by the lead (the lane's text below is unchanged).** In section 1 and in step 15 the lane writes the factor as 7.755. It is (1/2)(1 + 0.32082/0.14972) = 1.5714. The headline should read: ||beta e_l||/||e_l|| >= c|g|/lambda - 1.5714 (e^x - 1 - x) >= c|g|/lambda - 15.509 e^x g^2/lambda^2. Everything downstream is right as printed: 1.5714 pi^2 = 15.509, C = 15.53, the thresholds 1.36e-4 and 6.80e-5, and l_1(g) = ceil(14714 |g|). The script uses (1/2)(1 + Rb), as its pre-fixed log 01 states, and the table row l = 16384, g = 1 confirms it: 15.509/lambda^2 = 5.78e-8, the printed bound.

**Checked by the lead, by hand.**
- Step 5: d/dr [U(r) U_0(-r)] = g U(r) D_1 U_0(-r), so S - 1 = g int U(r) D_1 U_0(-r) dr and R_2 = g int (U(r) - U_0(r)) D_1 U_0(-r) dr.
- Step 6: with d = ||g D_1|| and x = d tau, ||R_2|| <= int_0^tau d (e^(dr) - 1) dr = e^x - 1 - x <= x^2 e^x/2 = pi^2 e^x g^2/lambda^2 at tau = pi.
- Steps 7 and 8: psi_l = i h^(3/2) e^(-eps h^2) u/sqrt(2); ||u||^2 = 4 * 2^17 (8!)^2/17! = 2.3963, equal to the lane's 4*65536/109395; the bound 0.3208 follows from s_3 <= 0.085899.
- Steps 10 and 13: p coth(pi p) >= 1/pi and Jensen give ||e_l||^2 >= (1/2 pi) A exp(-2 eps G^2/A); with A = 0.938463, G = 0.770250 this is 0.022417, so ||e_l|| >= 0.14972.
- Assembly: c g/lambda - 15.5 g^2/lambda^2 >= c g/(2 lambda) exactly when g/lambda <= c/31.0 = 6.8e-5.
- Not hand-checked: step 3 (bounded perturbation theorem, quoted from memory), step 9's certified sup s_3, step 11 (domain), and the inherited items the lane lists.

**Re-run by the lead:** `..\..\even-nonzero\venv\Scripts\python.exe -B even_remainder.py <out.json>` from `scripts\`, 80 s, exit 1 with exactly one failed line, P5, as the lane reports; P5b and P6b pass; the table equals the lane's. The diagnostic `diag_opnorm.py` was not re-run by the lead.

**On the failed line P5.** On the pre-fixed box the grid operator norm of D_1 exceeded the proved bound at l >= 256 (1.62, 2.19, 2.43 against 1.414). The lane predicted in writing, before running a diagnostic, that this is the unresolved kernel near the bottom of the box (the artefact EVEN-PROOF met), and the predictions held: the excess halves with the grid step and disappears on a box that resolves it. The replacement lines were written after seeing run 1. The proof of the bound (EVEN-PROOF) does not use the grid.

**Where the even divergence stands after this lane.** Every analytic step now has a written proof with a lead-checked core: the uniform bound, the limit operator, non-vanishing, a uniform first-order floor on an explicit family, the definition of the antilinear part, and the remainder. For l >= 14714 |g| the Hilbert-Schmidt norm of the antilinear part is at least c g/(2 lambda), and the sum over sectors of (2l+1)(c g/2 lambda)^2 diverges logarithmically. It is **still not to be called a theorem** until the inherited and quoted items are closed: the dictionary making W_l the even charge (EVEN-PROOF step 1, validated only through the spectrum); density of smooth compactly supported data in H; the two standard theorems (Araki duality and Shale's criterion; the bounded perturbation theorem) read first-hand; the Connes/Shale criterion for weights and the block-diagonal structure of the even tower from HINGE; and an independent re-run of the ball-arithmetic inputs. The paper may say: "for even parity the divergence follows from the steps in the record, each with a written proof; the remaining items are a coordinate dictionary validated numerically and standard theorems not yet re-read."

---

## Lane report (unchanged)

**EVEN-REMAINDER report (2026-11-15)**

The remainder closes for this family: (a), (b) and (c) are each proved, with certified constants. One pre-fixed line failed (P5, first clause); a diagnostic with predictions fixed beforehand showed it to be a grid artefact.

## 1. Final inequality

With c = 2.110·10⁻³, x = √2·π|g|/λ, and for every l ≥ l₀ = 16 and every real g:

‖β^(l)e_l‖/‖e_l‖ ≥ c|g|/λ − 7.755·(eˣ − 1 − x) ≥ c|g|/λ − 15.509·eˣ·g²/λ².

- **Constant C:** for |g|/λ ≤ 2.2·10⁻⁴ one may take C = 15.53.
- **Positive bound:** the right side is positive for |g|/λ < c/C = 1.36·10⁻⁴.
- **Half floor:** ‖β^(l)‖²_HS ≥ (cg/2λ)² for all l ≥ l₁(g) = ⌈14 714·|g|⌉, i.e. |g|/λ ≤ 6.80·10⁻⁵.
- **Where it is useless:** for g = 1 the bound is vacuous below l ≈ 7 400; for g = 8 below l ≈ 59 000.

## 2. Proof

**(a) Second-order remainder of the cocycle.**
1. [P] H is the real Hilbert space with ‖k‖² = ⟨z, coth(πh)z⟩, z = (h^{1/2}q + i·h^{−1/2}π)/√2. The vacuum flow is z ↦ e^{−ish}z, unitary on L² and commuting with coth(πh), so it is isometric and strongly continuous on H in exactly this thermal norm.
2. [I, EVEN-PROOF, lead-checked] d := ‖gD₁‖_H ≤ √2|g|/λ.
3. [S, from memory: bounded perturbation theorem, Engel–Nagel III.1] U(s) = e^{s(X₀+gD₁)} exists and U(s) = U₀(s) + g∫₀^s U(r)D₁U₀(s−r)dr.
4. [P] Gronwall on step 3 gives ‖U(s)‖ ≤ e^{ds}, hence ‖U(r) − U₀(r)‖ ≤ e^{dr} − 1.
5. [P] Multiplying step 3 by U₀(−τ): S − 1 = g∫₀^τ U(r)D₁U₀(−r)dr, so R₂ = g∫₀^τ (U(r) − U₀(r))D₁U₀(−r)dr.
6. [P] ‖R₂‖_H ≤ d∫₀^τ(e^{dr} − 1)dr = eˣ − 1 − x ≤ x²eˣ/2 = π²eˣg²/λ² at τ = π. So K = π²·e^{√2π|g|/λ}, which agrees with HINGE's corollary.

**(b) Numerator constants.**
7. [P] ψ_l = z_l = i·h^{−1/2}a_l/√2 = i·h^{3/2}e^{−εh²}u/√2, and ‖coth^{1/2}ψ_l‖ = ‖e_l‖.
8. [P] ‖coth^{3/2}ψ_l‖² = ½⟨u, h³coth³(πh)e^{−2εh²}u⟩ ≤ ½·s₃·‖u‖², with s₃ = sup_p (p·coth πp)³e^{−3p²}; ‖u‖² = 4·65536/109395 is independent of l.
9. [C] s₃ ≤ 0.085899, from 6 144 cells using that p·coth(πp) increases ([P]: sinh 2x ≥ 2x), plus a tail below 10⁻⁴⁴. Hence ‖coth^{3/2}ψ_l‖ ≤ 0.32082. Ball radii are below 10⁻²⁶.

**(c) Lower bound on ‖e_l‖.**
10. [P] p·coth(πp) ≥ 1/π gives ‖e_l‖² ≥ (1/2π)⟨hu, e^{−2εh²}hu⟩ ≥ (1/2π)·A_l·exp(−2εG_l²/A_l), by Jensen on the spectral measure of hu. This is increasing in A and decreasing in G.
11. [P] u ∈ D(h_l²): u is C³ with compact support in ρ > 0 (needs l ≥ 4), and H² functions of compact support lie in the closure of the minimal operator. This is EVEN-PASSAGE step 1, re-derived here.
12. [I, EVEN-PASSAGE step 2, lead-checked] A_l ≥ A_∞ and G_l ≤ G_hi.
13. [C] A_∞ ≥ 0.938463 and G_hi ≤ 0.770250 give ‖e_l‖ ≥ 0.14972 for all l ≥ 16.

**Assembly.**
14. [P, re-derived from BETA-DEF §6] β = ½((Ŝ−1) + i(Ŝ−1)i) and i·k(ψ) = k(i·coth ψ) + k′(…), so ‖βe − gβ̂₁e‖ ≤ ½‖R₂‖_H(‖e‖ + ‖coth^{3/2}ψ‖). That formula itself is [I, BETA-DEF step 9, not lead-checked].
15. [C] Dividing by ‖e_l‖: the factor is ½(1 + 0.32082/0.14972) = 7.755, so C₀ = 7.755·π² = 15.509.
16. [I, EVEN-PASSAGE, lead-checked except its step 1, step 3's closed form and its ball inputs] The first-order floor c; `certified()` re-run here returns 2.1102·10⁻³.
17. [P] ‖β‖_HS ≥ ‖βe‖/‖e‖.

## 3. Check

Run 2, 82 s. Columns 3 to 7 are on the fine grid (dy = 0.02); λ²‖R₂‖_H is on the box [−3, 40] at dy = 0.04.

| l | g | actual ‖βe‖/‖e‖ | g‖β̂₁e‖/‖e‖ | actual remainder | proved remainder bound | proved lower bound | λ²‖R₂‖_H (bound) |
|---|---|---|---|---|---|---|---|
| 16 | 1 | 7.2017e-3 | 7.3845e-3 | 2.530e-4 | 6.250e-2 | −6.24e-2 | 0.468 (10.82) |
| 64 | 1 | 1.9090e-3 | 1.9218e-3 | 1.776e-5 | 3.815e-3 | −3.78e-3 | 0.589 (10.10) |
| 256 | 1 | 4.8456e-4 | 4.8539e-4 | 1.144e-6 | 2.371e-4 | −2.29e-4 | 0.648 (9.93) |
| 1024 | 1 | 1.21606e-4 | 1.21658e-4 | 7.202e-8 | 1.480e-5 | −1.27e-5 | 0.678 (9.88) |
| 16384 | 1 | 7.60955e-6 | 7.60975e-6 | 2.820e-10 | 5.778e-8 | +7.10e-8 | 0.690 (9.87) |
| 16 | 8 | 4.928e-2 | 5.908e-2 | 1.401e-2 | 8.60 | −8.60 | 28.7 (1489) |
| 1024 | 8 | 9.6995e-4 | 9.7327e-4 | 4.597e-6 | 9.567e-4 | −9.40e-4 | 43.4 (639) |
| 16384 | 8 | 6.0865e-5 | 6.0878e-5 | 1.804e-8 | 3.700e-6 | −2.67e-6 | 44.2 (632) |

- **Inequalities:** every one holds at every l and both g. At the four required l the proved lower bound is negative, so it holds vacuously there; l = 16384 is the only row where it bites (7.10e-8 against c/2λ = 6.44e-8).
- **Slack:** the proved remainder is 205 to 247 times the actual one (λ²·remainder ≈ 0.0757 against 15.5). The true ratio stays at about 0.12/λ, 58 times the floor.
- **Grid norms:** ‖e_l‖ = 0.26768 (proved range 0.1497 to 0.4715); ‖coth^{3/2}ψ_l‖ = 0.29047 (proved ≤ 0.3208).
- **Sanity:** remainder scales as g² (ratio ≈ 64); the first-order column reproduces EVEN-PASSAGE to 2·10⁻⁶; the cocycle agrees with a Taylor matrix exponential to 3·10⁻¹²; halving dy moves the remainder by 6·10⁻⁵.
- **Plant A** (g²/λ³ in place of g²/λ²): violated at l ≥ 256. Caught.
- **Plant B** (2π dropped, τ = 0.5): planted 0.25 against measured 0.47 to 0.69, violated at every l. Caught.
- **FAILED as pre-fixed, P5 (first clause):** on the box [−4.8, 40] at dy = 0.04, λ‖D₁‖_H = 1.62, 2.19, 2.43 at l = 256, 1024, 16384, above √2. The exact kernel is unresolved at the box bottom, the artefact EVEN-PROOF reported for its C2. Predictions were fixed in log 03 before the diagnostic and all four held: the value halves with dy (1.257) and is 1.003 and 0.993 on the box [−3, 40]. Lines P5b, P6b and N8 on that box were written after seeing run 1; they pass, with λ‖D₁‖ from 0.756 to 1.016.

## 4. Re-run

From `record\wave73\even-remainder\scripts`:

`..\..\even-nonzero\venv\Scripts\python.exe -B even_remainder.py ..\out\even_remainder.json` (82 s, exits 1 because of P5)

The diagnostic is `diag_opnorm.py`, run the same way. Logs are `logs\01` to `05`. Nothing was written outside the folder; no commit, no REPORT.md.

## 5. Weakest step

Step 3 with step 2: the perturbed group on H and its Duhamel formula are a textbook theorem quoted from memory, and they need D₁ bounded on all of H, which rests on EVEN-PROOF's assumed density of C_c^∞ in H.

## 6. Still assumed or unverified before "the even Shale sum diverges logarithmically" is a theorem

- **Usefulness only at l ≥ 14 714·|g|.** This suffices for divergence of the tail, but the table tests the bound non-vacuously at one l only.
- **Ball-arithmetic inputs** (this lane's and EVEN-PASSAGE's six): not independently re-run; the cell-sup code is mine.
- **EVEN-PASSAGE:** step 3's closed form for Q_l is not lead-checked.
- **BETA-DEF:** steps 3, 5 and 9 are not hand-checked; step 4 (Araki duality, Shale necessity) is from memory; boundedness of S ⊕ 1 on K is inherited from HINGE.
- **Bounded perturbation theorem** and density of C_c^∞ in H (step 3 above).
- **EVEN-PROOF step 1:** the dictionary making W_l the even charge, with g = 1, is validated only through the spectrum.
- **HINGE items:** the Connes/Shale criterion for weights is not lead-checked; the even tower being block diagonal with a well-defined charge (step 9) and the physics gaps (steps 1, 12) remain open.
- **Sum over sectors:** (2l+1) multiplicity against (c/2λ)² is HINGE step 8, not re-derived here.
