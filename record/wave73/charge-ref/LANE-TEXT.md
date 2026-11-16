## Lane report (unchanged)

**CHARGE-ALL-L final report (16 Nov 2026)**

## 1. Verdict per part

- **(a) Odd, all l: proved.** On odd TT data A_ext = c_l E_2[ψ_B, ψ̇_B], with no edge response. ψ_B is built from K alone and ψ̇_B from h alone, so Z = ½c_l∫ψ_B² sits on extrinsic-curvature data.
- **(b) Even momentum block, all l: proved.** A_ext = ∫N|K|² − ∮κ²/(l(l−1)) = c_l[½∫Nψ̇_E² + ½⟨ψ̇_E, G_Dψ̇_E⟩], with the corner coefficient derived for every l.
- **(c) Even field block, all l: proved.** A_ext = c_l[½∫N|∇ψ_E|² + ½∫ψ_E²], with no excess; the stationary shape 2/(l+2) is derived for every l.
- **(d) Cross block: zero, by time reflection.** It has no arithmetic line.

Nothing is "exact at low l only". The stated conditions are l ≥ 2, profiles smooth in r², and the inherited definition of A_ext (steps 1, 2 and 5 below).

## 2. Proofs

Conventions: d = 4, unit ball, N = (1−r²)/2, 1/8G dropped, per N_Y = ∮φ², φ = r^l Y, q̃ = q″ + (2l+2)q′/r, A = −2 × (second-order part of F = Area − 2·Vol).

**Definition chain (l-independent)**
1. [I, companion A8 Prop. 4(i), not refereed; C at l = 2, 3, 4] On TT data A_un = ∫N K·K + ¼∫N|∇h|² − ½∮|h·ν|² + ¼∮(ν·h·ν)². Line X01 recomputes it from Christoffel symbols with no TT simplification, on nine even, odd and mixed data.
2. [P] The second-order metric h₂ is never solved for. Green's identity with ∇∇N = −δ and N = 0 on the edge gives δF[h₂] = −½∫N R_lin[h₂], so only the N-smeared Hamiltonian constraint enters. A solution exists, for example h₂ = 2uδ with u solving a Poisson equation.
3. [P] Edge drive. On the surface r = 1 + σ, the first-order area change of a TT h is −½∮ν_Σ·h·ν_Σ. Its part linear in σ, using div(h·x) = 0, is ∮σG with G = ½(∂_rψ + 2ψ) at r = 1, ψ = x·h·x. So G = ½(q′(1) + (l+2)q(1))Y. [S] The stiffness is ½(l−1)(l+2)∮σ². Hence A_ext = A_un + ∮G²/((l−1)(l+2)). [C] Line X02: A_un(T + 2 sym∇ξ) = A_un(T) − 2∮σG − (l−1)(l+2)∮σ², σ = ξ·ν, at l = 2, 3, 4.
4. [P] Corner and slice. For TT K and any X: ∫N(|K + Hess X|² − (ΔX)²) = ∫N|K|² + 2∮κτ + ∮|∇_tan τ|² − 2∫|∇X|², with τ = X on the edge and κ = x·K·x. This is two integrations by parts and valid for every X [C, line C11, l = 2 to 9]. The same drive −∮τκ follows from the area element at t = ετ, with ∂_tH = 2K. Stationary: X harmonic gives ∫|∇X|² = l∮τ², so the form is 2∮κτ + l(l−1)∮τ², τ = −κ/(l(l−1)), value −∮κ²/(l(l−1)).
5. [I, EVEN-TOWER, by hand, lead-flagged] Moving the region equals a gauge move of the data, and fixed true volume equals Area − 2·Vol at second order.
6. [P] (d). A_un is even under K → −K. By t → −t, h drives only σ and K drives only (τ, X), and the deformation form has no σ–τ term. So the stationary value has no h–K term. E_1 has none either.
7. [P] Odd data have parity (−1)^(l+1) and no scalar of degree l has that parity, so there is no drive: A_ext = A_un.

**Radial reduction (all l)**
8. [I refereed (CHAIN-ALL-L, B22/B27); C at l = 2 to 9] T[q], S[η], ∇²T[q] = T[q̃], ∇²S[η] = S[η̃]. Chain data: ψ_E = −½q̃_h φ, ψ̇_E = −q̃_K φ, ψ_B = l(l+1)η_K φ, ψ̇_B = ½l(l+1)η̃_h φ.
9. [P by hand; C at l = 2 to 9] The ten sphere Gram entries of (Hess φ, P_φ, φδ, φX) and the three odd ones, from ∮∇²F = α(α+1)∮F for F homogeneous of degree α.
10. [P] (b). Pointwise, ⟨T[q], T[q]⟩_r − (c_l/2)q̃² = r^(−2l−2)(r^(2l+2)U)′, with U = ((2l+1)q²/r³ + 2qq′/r² + 2q′²/(r(l+1)(l+2)))/(l(l−1)) [C symbolic in l, S03]. Integrate against N with N′ = −r: ∫N|T|² = (c_l/2)∫Nq̃²r^(2l+2) + q(1)²/(l(l−1)) + (c_l/2)∫r^(2l+2)q′². The corner term cancels q(1)² exactly. Since ∇²((q − q(1))φ) = q̃φ, G_Dψ̇_E = (q − q(1))φ and ⟨ψ̇_E, G_Dψ̇_E⟩ = ∫r^(2l+2)q′².
11. [C symbolic in l, arbitrary q, S05] (c). Density of A_ext minus density of c_l E_1 equals r^(−2l−2)(r^(2l+3)W)′ with W(1) = −(boundary terms) and r^(2l+3)W → 0 at 0. W is companion A8's Theorem 3 potential at n = 3, tested here against independently built densities.
12. [C symbolic, S09a, S09b] (a). K block: the difference is c_l l²(l+1)²(N r^(2l+1)η²)′. h block: companion A8's V.
13. [C symbolic in (l, a, b), S04, S06 to S08] Independently, all four Gram identities hold on r^(2a), r^(2b) as rational-function identities.
14. [P; C symbolic, S12] The static datum −2φδ equals T[q*] − 2 sym∇ξ*, with q* = −l(l−1)r²/(2l+3). Then σ_TT + ξ*·ν = 2/(l+2) and A_ext/N_Y = 2l(l−1)/((l+2)(2l+3)).

**Not re-proved here:** [I] that c_l E_1 is the vacuum modular generator.

## 3. Check

Run 02 (registered, l = 2 to 9): 127 lines pass, 0 fail, 36 plants caught, 0 missed, 21 s.

| l | c_l | static A_ext/N_Y | shape | T[s] momentum A_ext, Z_E | first rung, even / odd |
|---|---|---|---|---|---|
| 2 | 1/6 | 1/7 | 1/2 | 8/27, 1/27 | 8/7, 4/3 |
| 3 | 1/30 | 4/15 | 2/5 | 2/33, 1/165 | 10/9, 5/4 |
| 4 | 1/90 | 4/11 | 1/3 | 4/195, 1/585 | 12/11, 6/5 |
| 5 | 1/210 | 40/91 | 2/7 | 2/225, 1/1575 | 14/13, 7/6 |
| 6 | 1/420 | 1/2 | 1/4 | 8/1785, 1/3570 | 16/15, 8/7 |
| 7 | 1/756 | 28/51 | 2/9 | 1/399, 1/7182 | 18/17, 9/8 |
| 8 | 1/1260 | 56/95 | 1/5 | 2/1323, 1/13230 | 20/19, 10/9 |
| 9 | 1/1980 | 48/77 | 2/11 | 1/1035, 1/22770 | 22/21, 11/10 |

- **EVEN-TOWER reproduced at l = 2:** 1/7, unextremised −6/7, 8/27 × 4/15 = 32/405, Z = 4/405, spectrum 8/7, 19/18, 34/33. Both spectra hold symbolically in l for three rungs.
- **Plant:** shape 2/(l+1) breaks (c) symbolically and at every l; the static datum gives 6/49 against 1/7 at l = 2. Three further plants (corner coefficient, weight 1 on odd K, c_l scaled) are caught everywhere.
- **Failed lines:** run 01 printed FAIL on S04, S06, S07, S08. The cause was my script's monomial integrator, which split coefficient from power wrongly; the pointwise lines had passed. I stopped run 01, replaced the integrator, and those lines pass in run 02. Log 01 is kept and annotated.
- **Departures from the pass lines:** sympy's `.det()` was replaced by an explicit 3×3 expansion, and a non-degeneracy condition was added to S10 and S11.
- **Unregistered extra:** l = 10 to 14, Cartesian only, 60 pass, 0 fail.

## 4. Re-run

From the lane's folder:
`C:\Python314\python.exe -u scripts\charge_all_l.py 2 9 out\charge_all_l_2_9.json`

## 5. Weakest step

Step 5: that deforming the region is a gauge move at second order is inherited and by hand; my Christoffel checks cover its spatial consequences at l ≤ 4 only, and the time shift only through the formula of step 4.

## 6. Bulletin

Read B1 to B27 at the start and B28 on notice. No change of course.
- B27: the role swap is confirmed and used (ψ_B from K, ψ̇_B from h).
- B22: built on.
- B11: its 2/(l+2) now holds for every l.
- B23: no nsimplify anywhere.
- B28: its one unproved input, A_ext = c_l E_2 on odd data, is supplied by (a), subject to steps 1 and 5.

## 7. Sentences

**May use:** "For every l ≥ 2 and either parity, the fixed-volume area-side charge, stationary in the edge shape and the corner's time shift, equals c_l E_1 of the sector's chain field plus a positive form on extrinsic-curvature data, ½c_l∫ψ_B² for odd parity and ½c_l⟨ψ̇_E, G_Dψ̇_E⟩ for even; given the second-order area formula on transverse-traceless data, this is an exact radial identity proved symbolically in l and checked from Cartesian integrals at l = 2 to 9."

**Must not use:** "The sector identity is verified for l ≤ 6 and conjectured beyond", nor "proved from first principles for every l" without naming the inherited area formula and the gauge-move step.
