# REPORT: lane EVEN-PROOF (wave 73)

## Lead's check and grade

Sunday 15 November 2026. Small proof lane (Fable): the analogue of HINGE-PROOF's Lemma A for the even-parity charge. The lane returned the text below; the lead saves it.

**Grade: 378 of 420.** One question, answered with a seven-line proof that survives a line-by-line check, a validation of the operator against even-tower's spectrum, and an honest account of two of its own pre-fixed lines failing in the first run. Marked down because the dictionary between the ball variables and the hyperbolic variables (step 1) is the lane's own and is validated only through the spectrum, and because run 02's pass lines were written after run 01 was seen.

**Checked by the lead, by hand, line by line. It holds.**
- Hamilton's equation for Z_E = (1/2)<pi, W pi> gives D_E(q, pi) = (g W pi, 0).
- Step 2: with r = tanh(rho/2), dr = N d rho, sinh rho = r/N, N' = -r. Multiplying -w_rr + lambda^2 w/r^2 = pi/N by N w and integrating by parts gives <N w, pi> = ||h w||^2 + (1/2)||N^(1/2) w||^2; the boundary term vanishes because w(1) = 0.
- Steps 3 and 4: N^2 <= h^2/lambda^2 gives ||h W pi|| <= ||pi||/lambda.
- Step 5: the adjoint gives W^2 <= h^-2/lambda^2; operator monotonicity of the square root gives W <= h^-1/lambda, equivalently W^(1/2) h W^(1/2) <= 1/lambda, hence W h W <= W/lambda <= h^-1/lambda^2.
- Steps 6 and 7: coth(y) <= 1 + 1/y above, x coth(pi x) >= max(x, 1/pi) below; the two terms are each at most 2||k||^2/lambda^2 and the prefactor (1/2)g^2 gives ||D_E k||^2 <= 2 g^2 ||k||^2/lambda^2.

**Re-run by the lead:** `C:\Python314\python.exe even_lemma_a_ext.py <out.json>` from `scripts\`, 2 min 52 s, exit 0, "ALL LINES HELD"; the table for l = 16 to 1024 equals the lane's to the printed digits, the two grids agree to 2.3e-6, and the planted exponent 1/2 is rejected.

**Not checked by the lead:** run 01 (`even_lemma_a.py`, which exits 1 for the two reasons the lane gives), the 30-digit eigenfunction check, the triangularity argument for the spectrum, and the claim that the kernel equals hinge-e's resolvent form beyond the lane's own 3e-5 figure.

**What it changes.** The even charge's perturbation is now a bounded operator with the same bound as the odd one, 2*sqrt(2)*pi*|g|/lambda, so the Duhamel-Gronwall corollary of HINGE-PROOF carries over word for word. The even divergence is still **not proved**: it needs the analogue of Lemma B (a non-zero first-order limit, i.e. the even constant positive as a theorem). The paper may say: "the remainder beyond first order is controlled in both parities; the even first-order limit is computed, not proved."

---

## Lane report (unchanged)

**EVEN-PROOF report (2026-11-15)**

The even analogue of Lemma A holds, with the same constant and exponent as the odd one. Two of my own pre-fixed lines failed in run 01; both are explained in section 4.

## 1. Statement

Write t = tanh(ρ/2), so r = t, N = 1/(1+cosh ρ) = (1−r²)/2 and dr = N dρ. With ψ̇_E = π/(N² sinh ρ):

- Z_E = ½⟨π, Wπ⟩ on L²(dρ), with kernel W(ρ,ρ′) = t_<^{l+1}(t_>^{−l} − t_>^{l+1})/(2l+1).
- Equivalently W = A⁻¹, where A is the operator of the form ∫N⁻¹(w′² + λ²w²/sinh²ρ)dρ; numerically it also equals N^{1/2}(h²+¼)⁻¹N^{1/2}.
- D_E(q, π) = (gWπ, 0), with g = 1 for the area-side charge.

**Lemma A_E.** ‖hW‖ ≤ 1/λ on L²(dρ), hence ‖2πD_E‖ ≤ 2√2·π|g|/λ on H for every l ≥ 1. So C = 2√2·π and a = 1.

## 2. Proof

1. [P] Change of variables. ½∫Nψ̇²r²dr = ½∫π²dρ, and ∫ψ̇G_Dψ̇ d³x = ∫(rψ̇)G_r(rψ̇)dr with rψ̇ = π/N. G_r is the Dirichlet Green function of −d²/dr² + λ²/r² on (0,1), which gives the kernel above.
2. [P] Identity. Take π in C_c^∞ and w = Wπ, so −w_rr + λ²w/r² = π/N, w(1) = 0, w ~ r^{l+1} at 0. Multiply by Nw and integrate by parts, using N′ = −r: ⟨Nw, π⟩ = ‖hw‖² + ½‖N^{1/2}w‖² ≥ ‖hw‖². Here ‖hw‖² = ∫(w_ρ² + λ²w²/sinh²ρ)dρ = ∫N(w_r² + λ²w²/r²)dr.
3. [P] N² ≤ 1/sinh²ρ ≤ h²/λ² (Lemma A step 1, re-derived), so ‖Nw‖ ≤ ‖hw‖/λ.
4. [P] ‖hw‖² ≤ ‖Nw‖‖π‖ ≤ ‖hw‖‖π‖/λ, so ‖hWπ‖ ≤ ‖π‖/λ. This extends to L² because W is bounded and h is closed.
5. [P] Taking adjoints, ‖Wφ‖ ≤ λ⁻¹‖h⁻¹φ‖, i.e. W² ≤ λ⁻²h⁻². Since W ≥ 0, Löwner–Heinz gives W ≤ λ⁻¹h⁻¹, so ‖h^{1/2}Wπ‖² ≤ λ⁻¹⟨π, Wπ⟩ ≤ λ⁻²⟨π, h⁻¹π⟩.
6. [P] x·coth(πx) ≤ x + 1/π, so ⟨Wπ, h coth(πh)Wπ⟩ ≤ λ⁻²(⟨π, h⁻¹π⟩ + ⟨π, h⁻²π⟩/π).
7. [P] coth(πx)/x ≥ max(1/x, 1/(πx²)), so both terms are at most 2‖k‖². Therefore ‖D_E k‖² ≤ 2g²λ⁻²‖k‖². ∎

A weaker route gives exponent ½ only: ‖W‖ = 1/(2l+3) [P] yields ‖D_E‖ ≤ |g|(2l+3)^{−1/2}. That is true but not the scaling.

## 3. Validation against even-tower

- [P/N exact] On w_k = r^{l+1+2k}(1−r²), the problem −w″ + λ²w/r² = μw/N is triangular with diagonal μ_k = (k+1)(2l+2k+3), for l = 2..6. So spec(1+W) = 1 + 1/((j+1)(2l+2j+3)), which is even-tower's spectrum.
- [N] In mpmath at 30 digits, the ρ-kernel applied to the eigenfunctions returns 1/μ_j at l = 2, 3 and j = 0, 1, 2, with error 4·10⁻³⁰.
- [N] Grid eigenvalues are 0.142859, 0.055557, 0.030305 at l = 2 (exact 1/7, 1/18, 1/33).
- [N] The kernel matches N^{1/2}(h²+¼)⁻¹N^{1/2} to 3·10⁻⁵, the form hinge-e's script uses.
- Plant: the Neumann kernel was caught (top eigenvalue 86/7 instead of 1/7).

## 4. Norm table [N]

Run 01 (exact kernel), run 02 (resolvent form, l up to 1024); the two finest grids agree to 2·10⁻⁶ in run 02.

| l | ‖D_E‖ (g = 1) | bound √2/λ | λ‖D_E‖ | λ‖hW‖ | local slope |
|---|---|---|---|---|---|
| 2 | 0.16514 | 0.57735 | 0.4045 | 0.4176 | |
| 4 | 0.11839 | 0.31623 | 0.5294 | 0.5438 | |
| 8 | 0.076733 | 0.16667 | 0.6511 | 0.6653 | |
| 16 | 0.045756 | 0.085749 | 0.7546 | 0.7674 | |
| 32 | 0.025667 | 0.043519 | 0.8341 | 0.8446 | −0.852 |
| 64 | 0.013811 | 0.021926 | 0.8908 | 0.8989 | −0.904 |
| 128 | 0.0072321 | 0.011006 | 0.9293 | 0.9353 | −0.939 |
| 256 | 0.0037223 | 0.0055137 | 0.9548 | 0.9590 | −0.961 |
| 512 | 0.0018951 | 0.0027595 | 0.9712 | 0.9742 | −0.975 |
| 1024 | 0.00095831 | 0.0013804 | 0.9818 | 0.9838 | −0.984 |

- The bound holds everywhere. λ‖D_E‖ rises towards about 1, so the exponent is exactly 1 and √2 is not sharp.
- **Failed as pre-fixed in run 01:**
  - C3: the slope over l = 16..128 was −0.899, outside my window [−1.08, −0.92]. The window was pre-asymptotic; run 02's local slopes rise monotonically to −0.984.
  - C2: λ‖hW‖ = 1.2845 at l = 128 on the coarsest grid (step 0.04). This is a sampling artefact: the exact kernel is unresolved near ρ = 0 once l·step is about 1, predicted value l·step/4 = 1.28, with singular vectors on single grid points. Finer grids give 0.9356.
- Run 02's lines were written after seeing run 01, and all held.
- Plants caught: bound √2/λ² is violated at every l ≥ 3; exponent ½ is rejected (√λ‖D_E‖ falls from 0.111 to 0.031).

## 5. Re-run

From `record\wave73\even-proof\scripts`:

- `C:\Python314\python.exe even_lemma_a.py ..\out\even_lemma_a.json` (about 3.5 min; exits 1 because of C2 and C3 above)
- `C:\Python314\python.exe even_lemma_a_ext.py ..\out\even_lemma_a_ext.json` (about 2.5 min; exits 0)

Logs are `logs\01_even_lemma_a.log` and `logs\02_even_lemma_a_ext.log`. Nothing was committed and no REPORT.md was written.

## 6. Weakest step

Step 1's dictionary ψ̇_E = π/(N² sinh ρ) is my own conformal map, checked only through the spectrum and hinge-e's unverified kernel form. Density of C_c^∞ in H is assumed, as in the odd lemma.

## 7. What Lemma B would need

With y = ρ − ln 2λ, show that λW_l converges strongly to the Liouville operator with kernel sinh(x_<)e^{−x_>}, x = e^{−y} (hinge-e's unverified limit), alongside h_l → (−d²/dy² + e^{−2y})^{1/2} in strong resolvent sense. Then, using ‖λ h_l W_l‖ ≤ 1 from this lemma as the uniform bound, show the first-order matrix element on a fixed wave packet has a non-zero limit, which needs hinge-e's constant a_E > 0 proved, not computed.
