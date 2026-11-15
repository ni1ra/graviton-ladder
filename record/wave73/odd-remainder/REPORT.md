# REPORT: lane ODD-REMAINDER (wave 73)

## Lead's check and grade

Sunday 15 November 2026. Small proof lane (Fable): bound the second-order remainder of the odd-parity flow on ODD-SIGN's explicit family, uniformly in l, so the odd tower's divergence has the same standing as the even tower's. The lane returned the text below; the lead saves it.

**Grade: 386 of 420.** All three parts closed with certified constants, all seventeen pre-fixed lines passing on the first and only run, the box chosen in advance so the grid artefact of the even lane did not appear, and one step where it did NOT copy the even lane because the copy would have been wrong (the Jensen step fails for the q slot; it used Cauchy-Schwarz instead and says why). Marked down because the bound, like the even one, bites only at very large l and is tested non-vacuously at one l, and Lemma A was imported without re-reading.

**Checked by the lead, by hand. It holds.**
- Step 4: the cocycle remainder, identical to EVEN-REMAINDER's steps 5 and 6.
- Step 5: for e = (q, 0), z = h^(1/2) q/sqrt(2) with q = h e^(-eps h^2)(h u), giving the spectral weights p^3 coth(pi p) e^(-2 eps p^2) and p^3 coth^3(pi p) e^(-2 eps p^2) against the measure of h u.
- Step 7: p coth(pi p) >= 1/pi, then J_0 = int (p e^(-eps p^2))(1/p) d mu <= (int p^2 e^(-2 eps p^2) d mu)^(1/2) ||u||, so ||e_l||^2 >= J_0^2/(2 pi ||u||^2). With J_0 >= 0.363568 and ||u||^2 = 2.39631: ||e_l|| >= 0.09370.
- Step 8: N3_hi = sqrt(0.085898 * 0.938807/2) = 0.2008.
- Assembly: K = (1/2)(1 + 0.20080/0.09370) = 1.5715, C_0 = K pi^2 = 15.51, and c g/lambda - 1.001 C_0 g^2/lambda^2 >= c g/(2 lambda) for g/lambda <= 6.563e-5, i.e. l_1(1) = 15237.
- Not hand-checked: step 9 (re-derivation of the formula for beta e), the certified sup s_3 (the same number the blind second code located at 0.085567 for the even family, so this bound is on the safe side).

**Re-run by the lead:** `..\..\even-nonzero\venv\Scripts\python.exe -B odd_remainder.py <out.json>` from `scripts\`, 51 s, exit 0, "ALL LINES HELD"; the non-vacuous row l = 16384, g = 1 gives +6.661e-8 against c/(2 lambda) = 6.220e-8.

**What it changes.** The odd tower now has the same chain as the even one: bound (HINGE-PROOF Lemma A), first-order floor (ODD-SIGN), remainder (here), so for l >= 15237 |g| the Hilbert-Schmidt norm of the sector's antilinear part is at least c g/(2 lambda) with c = 2.038e-3. HINGE-PROOF's Lemma B outline, its Bessel identity and its unchecked corollary are no longer needed for the divergence. What remains is common to both towers and is being assembled by lane ASSEMBLY: boundedness of the extension to the purification, the standard theorems, and the physics identifications.

---

## Lane report (unchanged)

**ODD-REMAINDER: the remainder closes for the odd family.** All 17 pre-fixed lines passed on the first and only run (52 s, exit 0).

## 1. Final inequality

For every l ≥ 16 and real g, with x = √2·π|g|/λ:

‖β^(l)e_l‖/‖e_l‖ ≥ c|g|/λ − K(eˣ − 1 − x) ≥ c|g|/λ − C₀·eˣ·g²/λ²

- **Constants:** c = 2.038·10⁻³, K = ½(1 + 0.20080/0.09370) = 1.5715, C₀ = Kπ² = 15.5105.
- **Positive** for |g|/λ < 1.3127·10⁻⁴, i.e. l above about 7618·|g|.
- **Half floor:** ‖β^(l)‖_HS ≥ c|g|/(2λ) for l ≥ l₁(g) = max(16, ⌈15236.03·|g|⌉), i.e. |g|/λ ≤ 6.563·10⁻⁵; l₁(1) = 15237.
- **K is not copied from the even lane.** It lands near the even 1.5714 because the odd ratio N3_hi/E_lo equals the even one times √(A₁₆/A_∞); the script gives 2.14309.

## 2. Proof

1. [P] The vacuum flow z ↦ e^{−ish}z commutes with coth(πh), so U₀ is isometric on H.
2. [I, HINGE-PROOF Lemma A, lead-checked per brief; I did not read it] d := ‖gD₀‖_H ≤ √2|g|/λ. Grid: λ‖D₀‖_H is 0.78 to 0.99.
3. [S, from memory: bounded perturbation theorem] U(s) = e^{s(X₀+gD₀)} exists, with U(s) = U₀(s) + g∫₀^s U(r)D₀U₀(s−r)dr.
4. [P] Gronwall on step 3 gives ‖U(r)‖ ≤ e^{dr}, so ‖U(r) − U₀(r)‖ ≤ e^{dr} − 1. Multiplying step 3 by U₀(−τ): R₂ = g∫₀^τ(U(r) − U₀(r))D₀U₀(−r)dr, so ‖R₂‖_H ≤ ∫₀^τ d(e^{dr} − 1)dr = eˣ − 1 − x ≤ π²eˣg²/λ².
5. [P] Weights for the q slot: ψ_l = h^{5/2}e^{−εh²}u/√2, real. With μ the spectral measure of hu (mass A_l, ∫p⁻²dμ = ‖u‖² = 4·65536/109395):
   - ‖e_l‖² = ½∫p³coth(πp)e^{−2εp²}dμ
   - ‖coth^{3/2}ψ_l‖² = ½∫p³coth³(πp)e^{−2εp²}dμ
6. [P] ‖coth^{3/2}ψ_l‖² ≤ ½·s₃·A₁₆, with s₃ = sup(p·coth πp)³e^{−2εp²}. A_l ≤ A₁₆ is [I, ODD-SIGN step 6].
7. [P] Lower bound on ‖e_l‖. The even lane's Jensen step fails here: t·e^{−2εt} is concave at the mean t ≈ 0.63 < 1/ε. Instead p·coth(πp) ≥ 1/π and Cauchy–Schwarz give ‖e_l‖² ≥ (1/2π)∫p²e^{−2εp²}dμ ≥ (1/2π)·J₀²/‖u‖². J₀ ≥ J₀_lo is [I, ODD-SIGN step 8, lead-checked].
8. [C] s₃ ≤ 0.085898 (6144 monotone cells plus tail), A₁₆ ≤ 0.938807, J₀_lo = 0.363568. Hence N3_hi = 0.20080 and E_lo = 0.09370.
9. [P, re-derived from BETA-DEF's definitions] i·k(ψ) = k(i·coth ψ) + k′(−i·csch ψ); I checked c·coth − s·csch = c and c·csch − s·coth = s. So βe = ½[k((S−1)ψ) + i·k((S−1)(i·coth ψ))] and ‖βe − gβ̂₁e‖ ≤ ½‖R₂‖_H(‖e‖ + ‖coth^{3/2}ψ‖). The β formula itself is [I, BETA-DEF steps 4 and 5, not all lead-checked].
10. [I, ODD-SIGN, lead-checked except its step 5 and ball inputs] ‖β̂₁e_l‖/‖e_l‖ ≥ c/λ; re-running its `certified()` returns 2.03807·10⁻³.
11. [C] Assembly: K = 1.571546, C₀ = 15.510542. For l ≥ l₁, x ≤ 2.9·10⁻⁴, so eˣ ≤ 1.001 and 1.001·C₀g²/λ² ≤ cg/(2λ).
12. [P] ‖β‖_HS ≥ ‖βe‖/‖e‖.

## 3. Check

Box y ∈ [−3, 40] for every l, chosen before running: decay rate e³ times dy is 0.40 at dy = 0.02 (vectors) and 0.80 at dy = 0.04 (operator norms).

| l | g | actual ‖βe‖/‖e‖ [N] | first order [N] | actual remainder [N] | proved remainder bound [C] | proved lower bound [C] | λ²‖R₂‖_H (bound) |
|---|---|---|---|---|---|---|---|
| 16 | 1 | 1.7523e-2 | 1.8033e-2 | 6.521e-4 | 6.251e-2 | −6.24e-2 | 0.579 (10.82) |
| 64 | 1 | 4.6665e-3 | 4.7024e-3 | 4.579e-5 | 3.816e-3 | −3.78e-3 | 0.646 (10.10) |
| 256 | 1 | 1.18596e-3 | 1.18827e-3 | 2.949e-6 | 2.371e-4 | −2.29e-4 | 0.669 (9.93) |
| 1024 | 1 | 2.97724e-4 | 2.97870e-4 | 1.857e-7 | 1.480e-5 | −1.28e-5 | 0.677 (9.88) |
| 16384 | 1 | 1.863199e-5 | 1.863256e-5 | 7.278e-10 | 5.778e-8 | +6.661e-8 | 0.680 (9.87) |
| 16 | 8 | 1.1654e-1 | 1.4426e-1 | 3.645e-2 | 8.602 | −8.60 | 33.1 (1489) |
| 16384 | 8 | 1.49024e-4 | 1.49060e-4 | 4.653e-8 | 3.700e-6 | −2.71e-6 | 43.5 (632) |

- **Non-vacuous** only at l = 16384, g = 1: +6.661e-8 against c/(2λ) = 6.220e-8. Every other table row holds vacuously, including all of g = 8, where l₁(8) = 121889.
- **Slack:** the proved remainder is 79 to 96 times the actual one at g = 1.
- **Grid norms:** ‖e_l‖ = 0.17456 (proved range 0.0937 to 0.2332); ‖coth^{3/2}ψ_l‖ = 0.18163 (proved ≤ 0.2008).
- **Sanity:** first order matches ODD-SIGN's table to 1.6·10⁻⁶; cocycle matches a Taylor exponential to 2·10⁻¹³; halving dy moves the remainder by 2·10⁻⁵; g² ratio is 55.9 to 63.9.
- **Grid artefact did not appear:** D₀ is a multiplication operator, and on the box [−4.8, 40] λ‖D₀‖_H is still at most 0.997.
- **Plant A** (cocycle factor e^{−τX₀} dropped): "remainder" 0.213 to 0.215 against the bound, violated at every l. Caught.
- **Plant B** (g²/λ³): violated at l = 256, 1024, 16384. Caught.
- **Failed lines:** none.

## 4. Re-run

From `record\wave73\odd-remainder\scripts`:

`..\..\even-nonzero\venv\Scripts\python.exe -B odd_remainder.py ..\out\odd_remainder.json`

Logs are `logs\01_prefixed_pass_lines.log` and `logs\02_odd_remainder_run1.log`. Nothing written outside the folder; no commit, no REPORT.md.

## 5. Weakest step

Steps 2 and 3 together: the perturbed group and its Duhamel formula are a textbook theorem from memory, and they need Lemma A's bound on all of H, which I imported without reading its proof.

## 6. Still assumed before "the odd Shale sum diverges logarithmically" is a theorem

- Lemma A for D₀ (lead-checked, not re-derived here) and density of C_c^∞ in H.
- Bounded perturbation theorem, read first-hand.
- BETA-DEF step 4 (Araki duality, Shale's criterion) and boundedness of S ⊕ 1, from memory or inherited from HINGE.
- ODD-SIGN step 5: essential self-adjointness of h_l², from memory.
- Ball-arithmetic inputs (ODD-SIGN's six, plus s₃ and ‖u‖² here), not independently re-run.
- That D₀(q, π) = (0, −Nq) is the odd charge, block diagonal per sector; the Connes/Shale criterion for weights; the (2l+1) multiplicity sum; HINGE's physics gaps (steps 1, 9, 12).
- The bound bites only for l ≥ 15237·|g| and is tested non-vacuously at one l.
