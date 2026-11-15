# REPORT: lane EVEN-DICT (wave 73)

## Lead's check and grade

Sunday 15 November 2026. Small proof lane (Fable): prove the change of variables that every even-parity proof lane assumed, from the ball's radial sector data to the hyperbolic variables in which the bounds were proved, including that the coupling is g = 1. The lane returned the text below; the lead saves it.

**Grade: 386 of 420.** The dictionary is proved as an exact identity for every l by four lines of integration by parts, checked to 45 digits on full Gram matrices (cross terms included) at l = 2, 3, 4, with three plants caught and the q-map, which EVEN-PROOF never wrote down, stated. Marked down because the lane is clear that it proves the change of variables only: that c_l E_1 with this symplectic form is the vacuum's modular generator on even data, with the thermal norm at inverse temperature 2 pi, remains inherited from EVEN-TOWER and from Hislop-Longo, and the thermal norm was tested here only in the momentum slot at l = 2.

**Checked by the lead, by hand. It holds.**
- Step 1: (1/2) int N p^2 r^2 dr with p = pi/(N r) and dr = N d rho is (1/2) int pi^2 d rho.
- Step 2: int p f r^2 dr with f = q/r is int pi q d rho.
- Step 3: r^2 f'^2 = q_r^2 - (q^2)_r/r + q^2/r^2 and (N/r)' = -1 - N/r^2, so int N r^2 f'^2 dr = int N q_r^2 dr - int q^2 dr; the boundary term vanishes at both ends; the remaining terms are (1/2) int (q_rho^2 + lambda^2 q^2/sinh^2 rho) d rho once the weight-1 term cancels the - int q^2.
- Step 4: the Green function of -d^2/dr^2 + lambda^2/r^2 on (0, 1) with Dirichlet condition at 1, and the measure factors, give <pi, W pi> with the stated kernel and no extra constant.
- Not checked by the lead: that EVEN-TOWER's E_1 contains exactly the weight-1 term used in step 3 (read from the lane's statement), steps 6 and 7.

**Re-run by the lead:** `C:\Python314\python.exe even_dict_check.py <out.json>` from `scripts\`, about 2.5 min, exit 0, "FAILED LINES: none"; relative errors of order 1e-45 at l = 2, 3, 4; plants P1, P2, P3 caught. `thermal_check.py` was not re-run by the lead.

**What it changes.** EVEN-PROOF's step 1 is no longer "validated only through the spectrum": the change of variables is exact and g = 1, so no threshold in EVEN-PASSAGE or EVEN-REMAINDER moves. What remains inherited on the even side is a physics identification, not a computation: that the vacuum modular generator on the sector is c_l E_1 in the thermal norm (EVEN-TOWER, Hislop-Longo).

---

## Lane report (unchanged)

**EVEN-DICT report (15 Nov 2026)**

The dictionary holds as an exact identity for every l, with g = 1; nothing failed. The proof is four lines of integration by parts, and the check agrees with it to 45 digits.

## 1. Theorem

Sector l, ψ = f(r)Y, ψ̇ = p(r)Y, integrals in units of N_Y = ∫Y². Put r = tanh(ρ/2), so N = (1−r²)/2 = 1/(1+cosh ρ), dr = N dρ, r = N sinh ρ, N′ = −r. Define

- **q(ρ) = r·f(r)**, i.e. ψ = q/(N sinh ρ)
- **π(ρ) = N·r·p(r)**, i.e. ψ̇ = π/(N² sinh ρ)

Then, with λ² = l(l+1) and f = O(r^l) at 0, bounded at r = 1:

- (a) c_l·E_1[ψ, ψ̇] = c_l N_Y·[½⟨π, π⟩ + ½⟨q, h²q⟩] on L²(dρ), and the symplectic form c_l∫(ψ̇_a ψ_b − ψ̇_b ψ_a) = c_l N_Y·∫(π_a q_b − π_b q_a)dρ.
- (b) ½c_l⟨ψ̇, G_D ψ̇⟩ = c_l N_Y·½⟨π, Wπ⟩, with the stated kernel and **g = 1**.
- (c) The constant c_l N_Y multiplies both the Hamiltonian and the symplectic form, so it cancels in Hamilton's equations: X_0(q, π) = (π, −h²q) and D_E(q, π) = (Wπ, 0).

## 2. Proof

1. [P] Momentum: ½∫N p² r² dr = ½∫N (π/(N r))² r² N dρ = ½∫π² dρ.
2. [P] Symplectic: ∫p f r² dr = ∫(π/N) q · N dρ = ∫π q dρ.
3. [P] Field: with f = q/r, r²f′² = q_r² − (q²)_r/r + q²/r². Integrating the middle term by parts with (N/r)′ = −1 − N/r² gives ∫N r² f′² dr = ∫N q_r² dr − ∫q² dr. The boundary term N q²/r vanishes at r = 1 because N = 0, and at 0 because q² /r ~ r^(2l+1). The −∫q² cancels the weight-1 term +∫f² r² = ∫q² exactly. What remains is ½∫N(q_r² + λ²q²/r²)dr = ½∫(q_ρ² + λ²q²/sinh²ρ)dρ.
4. [P] Z term: G_D ψ̇ = U(r)Y with w = rU solving −w″ + λ²w/r² = r p = π/N, w(1) = 0, w regular at 0. The solutions r^(l+1), r^(−l) have Wronskian −(2l+1), so the Green function is r_<^(l+1)(r_>^(−l) − r_>^(l+1))/(2l+1). Then w(ρ) = ∫G(r, r′)(π/N)(r′) dr′ = ∫W(ρ, ρ′)π(ρ′)dρ′ and ∫p U r² dr = ∫(π/N) w N dρ = ⟨π, Wπ⟩.
5. [P] Hamilton's equations for ½⟨π, π⟩ + ½⟨q, h²q⟩ + ½⟨π, Wπ⟩ give (c).
6. [I] That c_l E_1 with this symplectic form is the vacuum's modular generator: EVEN-TOWER Q3, lead re-ran at l = 3 only; Hilbert-space steps inherited there.
7. [N]/[I] Same thermal norm: the momentum slot is confirmed numerically (below). The q slot and the general statement rest on Hislop–Longo, inherited, not read first-hand.

The odd lanes justified the same map by HINGE's `s1` exact checks (φ = φ̃/N, u = sinh ρ × radial), which is this q-map. I did not re-run `s1`.

## 3. Check

Main check, `even_dict_check.py`: 5 data (r^(l+2k), k = 0..3, plus a degree-8 mix), full 5×5 Grams. The r side is exact rationals; the ρ side is Gauss–Legendre at 45 digits, with the kernel as a true double integral.

| l | (a1) momentum | (a2) field | (a3) symplectic, 25 entries | (b) Z, g = 1 | max entrywise deviation of Z_ρ/Z_r from 1 |
|---|---|---|---|---|---|
| 2 | 1.7e-45 | 2.7e-45 | 1.9e-45 | 2.4e-45 | 2.5e-45 |
| 3 | 6.2e-46 | 1.5e-45 | 2.7e-45 | 5.0e-46 | 1.6e-45 |
| 4 | 1.0e-45 | 1.2e-45 | 1.0e-45 | 1.1e-45 | 1.9e-45 |

(Relative errors against the largest entry.)

- **Exact side checks:** G_D solves the 3D Dirichlet problem for every datum, and the r-side spectrum is 1/((j+1)(2l+2j+3)), j = 0..3, as in EVEN-TOWER (Z/A = 1/7 at l = 2).
- **Quadrature:** 60 against 80 nodes agree to about 1e-45.
- **Plants, all caught at every l:**
  - π = N²rp (one power of N too many) breaks (a1) by 0.94–0.96 and (a3) by 0.85–0.92.
  - q = N r f breaks (a2) by 0.91–0.96.
  - A Neumann kernel breaks (b) by 1.7–1.8.
- **Thermal norm, secondary** (`thermal_check.py`, l = 2, p ∈ {r², r⁴}): ⟨p, ω⁻¹p⟩ on ℝ³ equals ⟨π, h⁻¹coth(πh)π⟩ to 1.5e-14 on the 2×2 Gram. Parseval holds to 2e-13; the wrong temperature coth(πk/2) is caught at 0.41.

**Failed pass lines: none.** Run 02 crashed on a rounding division in my Legendre-Q formula before printing any line; it was fixed and re-run as log 03. The physics was fine; my logarithm was not.

## 4. Re-run

From `record\wave73\even-dict\scripts`:

- `C:\Python314\python.exe even_dict_check.py ..\out\even_dict_check.json` (about 2.5 min, exit 0)
- `C:\Python314\python.exe thermal_check.py ..\out\thermal_check.json` (about 70 s, exit 0)

Pass lines are in `logs\00_passlines.log`, written before the first run; logs 01 to 03 sit beside it.

## 5. Weakest step

Step 6: the identities (a) to (c) hold for every l, but "c_l E_1 is the vacuum's modular generator on even data, with the β = 2π thermal norm" is still EVEN-TOWER's Q3 plus Hislop–Longo, and the thermal norm is tested here only in the momentum slot at l = 2.

## 6. What must change downstream

Nothing in the mathematics: g = 1, same h, same W, so no threshold moves. I did not read EVEN-LIMIT, EVEN-NONZERO, EVEN-PASSAGE or EVEN-REMAINDER, so this rests on your description of them.

- **Step 1 of EVEN-PROOF:** may be retagged from "own conformal map, validated through the spectrum" to proved for (a) to (c), citing this lane and stating the q-map q = rψ, which EVEN-PROOF never wrote down.
- **Constant:** the physical one-particle norm is c_l N_Y times HINGE-PROOF's (up to the record's factor-2 convention). It is one overall constant and cancels in every operator norm and Hilbert–Schmidt criterion.
- **Still listed as inherited:** the vacuum-side identification (step 6), not the change of variables. The dictionary is now the one part of the bridge nobody needs to take on trust, which makes it the odd one out.
