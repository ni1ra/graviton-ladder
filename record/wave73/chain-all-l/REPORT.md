# REPORT: lane CHAIN-ALL-L (wave 73)

## Lead's check and grade

Sunday 15 November 2026. Small proof lane (Fable): prove for every l >= 2 and both parities that each graviton sector reduces to its scalar chain field: the chain field obeys the scalar wave equation, it and the sector's curvature determine each other, and the helicity-2 one-particle inner product is c_l times the scalar's, c_l = 4 (l-2)!/(l+2)!. Until now this was exact at l = 2 to 6 and cited to an older note for general l. The lane returned the text below; the lead saves it. It sent an interim result to the lead and read the bulletin (B1 to B20), as the protocol asks.

**Grade: 388 of 420.** All three parts are proved for every l by short arguments, the key one a shell-by-shell identity on the momentum sphere that needs no radial derivative; the check is exact at l = 2 to 8 on full Gram matrices, with an unregistered run at l = 9 to 12; it found and stated a correction to the record's wording (the Cauchy datum is determined only up to the residual gauge modes on the ball; the curvature is determined exactly); and it reported that its own registered plant is vacuous at l = 2, with a second plant added and labelled post hoc. Marked down because the span of the tensor families in each sector is quoted from memory, the normalisation was not reconciled with the record's units, and it did not reproduce the record's earlier numbers digit for digit (it could not import the earlier machinery inside its read fence).

**Checked by the lead, by hand. It holds.**
- (i).1: Laplacian(x^i x^j M_ij) = x^i x^j Laplacian M_ij + 4 x^i d^j M_ij + 2 delta^ij M_ij for symmetric M; tracelessness and vanishing divergence remove the last two terms.
- (iii).2 and 3: for Y_AB = (D_A D_B - (1/2) gamma_AB D^2) Y_lm on the unit sphere, D^A D^B Y_AB = (1/2) Lambda Y and the integral of |Y_AB|^2 is (1/2) Lambda times that of Y^2, with Lambda = (l+2)!/(l-2)! (standard). So with psi-hat = -(1/k^2) D^A D^B E-hat_AB on a shell, the integral of |psi-hat|^2 is (Lambda/2)/k^4 times that of |E-hat|^2, and 2 |E-hat|^2/k^4 integrates to (4/Lambda) times |psi-hat|^2: c_l = 4/Lambda.
- Not hand-checked: the explicit radial formulas of (ii).2 to (ii).4 (the script's exact transverse-traceless checks cover them at l = 2 to 8), step (ii).6, the odd-sector curl argument.

**Re-run by the lead:** `C:\Python314\python.exe scripts\chain_all_l.py 2 8`, 3.2 s, exit 0: 113 of 113 lines PASS, 27 plants caught, one registered plant MISSED at l = 2 exactly as the lane reports (the planted factor equals the true one there, since 1! = 0!).

**What it changes.** The general-l reduction, named as open by MODFLOW-REF and MODGEN, is closed at the one-particle level. Together with MODFLOW-REF and MODGEN, the graviton-specific bridge now reads: (1) definition: the ball's graviton algebra is the linearised-curvature net; (2) refereed: its vacuum modular group acts on curvature as the conformal flow; (3) proved here for every l: each sector of that curvature is the scalar chain field, with the inner product scaled by c_l. What remains an assumption is (1) alone, plus second quantisation as a standard step. Wording to fix wherever the record says a "datum" is determined by its chain field: the curvature is; the transverse-traceless datum on the ball is, up to residual gauge.

---

## Lane report (unchanged)

**CHAIN-ALL-L final report (15 Nov 2026)**

## 1. Verdicts

- **(i) Wave equation: proved for all l**, both parities. It needs only E, B in C² (or distributions), symmetric, traceless, divergence-free, and Box E = Box B = 0.
- **(ii) Determination: proved for all l, with one stated correction.** The curvature (E, B) of a sector and the pair (ψ, ψ̇) determine each other, bijectively onto all smooth sector scalars r^l F(r²) Y_lm. The TT Cauchy datum (h, K) is determined only modulo the residual gauge modes on the ball (h = Hess φ or sym∇(x×∇φ), K = Hess φ). On ℝ³ with decay the kernel is zero.
- **(iii) Two-point function: proved for all l**, both parities, c_l = 4(l−2)!/(l+2)!. The proof is at the one-particle level; passing to the algebra uses second quantisation, as before.

## 2. Proofs

Conventions: φ = r^l Y_lm, Λ = (l+2)!/(l−2)!, Dq := q′/r, curl(T)_ij = ε_ikl ∂_k T_lj. The record's signs are E = −½∇²h, B = −curl K, Ė = −curl B, Ḃ = curl E.

**(i)**
1. [P][C] For any symmetric M: ∇²(x·M·x) = x·∇²M·x + 4x·div M + 2 tr M (line L01, symbolic).
2. [P] For M traceless and divergence-free, ∇²(x·M·x) = x·∇²M·x. Since ∂_t commutes with x, Box(x·E·x) = x·Box E·x = 0, and the same for B.
3. [S, from memory] E and B of vacuum linearised gravity have these properties.

**(ii)**
1. [S, from memory] Even-type symmetric tensors of a sector are f₁Hess φ + f₂P_φ + f₃φδ + f₄φ x⊗x. Odd-type are g₁(x⊗v + v⊗x) + g₂(∇v + ∇vᵀ), v = x×∇φ.
2. [I, re-derived, C] Even TT (companion A8 Lemma 7 at n = 3), with x·T·x = qφ:
   - β = −Dq/(l+1), α = (q + r²Dq/(l+1))/l
   - f₂ = (β − 2Dα)/(l+2), f₃ = −Dα − (l+1)f₂
   - f₄ = −(2l f₂ + 3f₃)/r², f₁ = (α − r²f₂)/(l−1)
   
   The linear system for f₂, f₃ has determinant l+2, so T is unique given q. T is smooth iff q is smooth in r². I re-derived the α, β step and the determinant; the rest is the script's TT check.
3. [P][C] Odd TT, new here. S·x = ηv with η = r²g₁ + (l−1)g₂. Using div(x⊗v + v⊗x) = (l+4)v and div(∇v + ∇vᵀ) = 0, div S = 0 gives g₁ = −Dη/(l+2), g₂ = (η + r²Dη/(l+2))/(l−1). Also x·S·x = 0 and x·curl S·x = (x×∇)·(S·x) = −l(l+1)ηφ.
4. [P][C] Maps between the two families, with q̃ = q″ + (2l+2)q′/r:
   - curl S[η] = T[−l(l+1)η], curl T[q] = S[q̃/(l(l+1))]
   - ∇²T[q] = T[q̃], ∇²S[η] = S[η̃]
5. [P] Even sector: E = T[q_E] with q_Eφ = ψ_E, and B = S[η] with ψ̇_E = −x·curl B·x = l(l+1)ηφ. Both are algebraic and local, so (E, B) ↔ (ψ_E, ψ̇_E) is a bijection. The odd sector is the same with E and B exchanged.
6. [P] From Cauchy data: h ↦ E has kernel q̃ = 0, i.e. q constant, i.e. Hess φ (gauge ξ = ∇φ/2); for odd h the kernel is sym∇v. K ↦ −curl K kills only Hess φ. The range is everything, since (r^{2l+2}q′)′ = r^{2l+2}F has a solution smooth in r² for every F smooth in r².

**(iii)**
1. [P] In momentum space, ψ̂ = −∂_{k_i}∂_{k_j}Ê_ij. For Ê tangential and traceless on the k-sphere, ∂_i∂_jÊ_ij = k⁻² D^A D^B Ê_AB (unit-sphere derivatives). There is no radial derivative.
2. [P] For Y_AB = (D_A D_B − ½γ_AB D²)Y_lm: D^A D^B Y_AB = ½ΛY and ∮|Y_AB|² = ½Λ∮Y². The odd harmonic has D^A D^B X_AB = 0, and k̂× maps X_AB to ±Y_AB isometrically.
3. [P] So on each shell and for any radial weight, 2|Ê|²/k⁴ integrated over the sphere equals (4/Λ)∮|ψ̂_E|². Hence ½∫ω|ĥ|² = c_l∫ω|ψ̂_E|², 2∫|K̂|²/ω = c_l∫|ψ̇_E|²/ω and ∫K·h = c_l∫ψ̇_Eψ_E. The odd sector follows from the same computation with the curl.
4. [P] The identity is shell-wise, so the map is isometric up to c_l and onto the scalar sector's one-particle space.
5. [S, from memory] The Fourier transform preserves (l, m, parity), which identifies the momentum-space harmonics with the position-space sectors.

**Normalisation.**
- [C] Graviton form: ½∫ω|ĥ|² + 2∫|K̂|²/ω with Ω = ∫K·h and K = ½ḣ. Scalar form: ∫ω|ψ̂|² + ∫|ψ̇|²/ω with Ω = ∫ψ̇ψ.
- [P] This is equivalent to modflow's ⟨EE⟩ = ½ω⁴P_TT times the scalar kernel.
- [U] By my reckoning this is 16πG = 1 for g = η + h. I did not reconcile it with the record's phrase "1/8G dropped, units of 4π".

## 3. Check

Run 03, l = 2 to 8, exact rationals, φ = Re[(x+iy)^{l−1}]z, profiles 1, s, s² in polynomial and Gaussian mode.

| l | c_l | L lines passed | ½∫ω h·h, entry [0][0] (even) | 2∫K·K/ω, entry [0][0] (odd) |
|---|---|---|---|---|
| 2 | 1/6 | 16/16 | 1/3 | 12/5 |
| 3 | 1/30 | 16/16 | 8/35 | 768/175 |
| 4 | 1/90 | 16/16 | 16/45 | 2560/189 |
| 5 | 1/210 | 16/16 | 1024/1155 | 30720/539 |
| 6 | 1/420 | 16/16 | 3072/1001 | 43008/143 |
| 7 | 1/756 | 16/16 | 40960/3003 | 7340032/3861 |
| 8 | 1/1260 | 16/16 | 16384/221 | 169869312/12155 |

- **Total:** 113 of 113 lines pass, including the symbolic L01; 3.1 s. Box ψ = 0 is verified through t³.
- **Extra run 04**, l = 9 to 12 (not registered): no fails.
- **Failed line:** the registered plant P1 ((l−1)! for (l−2)!) printed MISSED at l = 2 because 1! = 0!, so it is vacuous there. It was caught at l = 3 to 8. P2 and P3 were caught at every l. P1b (c_l scaled by (l+3)/(l+2)), added post hoc in `logs\02_addendum.log`, was caught at every l.
- **Not done as briefed:** I did not reuse the record's machinery, because even-tower's code imports `polyalg` from charge-id, a folder outside my read list. The script is self-contained, so even-tower's 8/3 and 16/15 are not reproduced digit for digit (different Gaussian width and harmonic).

## 4. Re-run

From `record\wave73\chain-all-l`: `C:\Python314\python.exe scripts\chain_all_l.py 2 8`

## 5. Weakest step

Step (ii).1, that the four even and two odd Cartesian tensors span each sector, is quoted from memory; the script's dimension counts (L04, L07) test it only inside a polynomial ansatz of bounded degree.

## 6. Bulletin

I read B1 to B18 at the start and B19 and B20 on notice. No change of course.
- B8 and B17 are consistent with this result.
- B19: the script uses exact QQ arithmetic throughout and never `nsimplify`.
- B20: it names c_l at general l as the open piece, which is (iii) above.

## 7. Sentences

**May use:** "In each sector (l ≥ 2, either parity) the chain field obeys the massless scalar wave equation, determines the sector's linearised curvature and is determined by it, and the helicity-2 one-particle inner product equals c_l = 4(l−2)!/(l+2)! times the scalar's on the chain field; this is proved for every l by a shell-by-shell identity on the momentum sphere and checked exactly for l = 2 to 8."

**Must not use:** "A TT Cauchy datum on the ball is determined by its chain field" (it is determined only up to the residual gauge modes), nor "c_l is verified numerically for l ≤ 4" as the status.

Files are in `record\wave73\chain-all-l\`:
- scripts\chain_all_l.py
- logs\00_passlines.log
- logs\01_run_l2_4.log
- logs\02_addendum.log
- logs\03_run_l2_8.log
- logs\04_run_l9_12.log
- out\chain_all_l_2_8.json

Nothing was committed and no REPORT.md was written.
