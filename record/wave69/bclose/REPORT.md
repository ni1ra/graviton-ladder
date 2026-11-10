# REPORT: w69-bclose. Bin BOTH-PASS (modifier Z1L1-DIRTY): no term linear in h⁽²⁾ escapes the smeared constraint in either construction the clause uses, and the finite-ε geometry gives the area side A to ten digits

Lane w69-bclose (Fable 5.2), Tuesday 10 November 2026, 01:31 Oslo time. Claude only. The PREREG was hashed at 01:00:45, before any script existed, and is unchanged (`sha256sum -c` passes). Nothing is priced here, nothing is committed, nothing outside `wave69/bclose` was edited.

## 1. The bin and the numbers

**Bin: BOTH-PASS, with Z1L1-DIRTY.** P1-PASS and P2-PASS as registered. No kill, so the false-kill rule is not triggered.

### Part 1: the route-negating enumeration (`scripts/p1_enum.py`, `out/p1_enum.json`, 24 checks, all pass, exact)

The second-order metric is h⁽²⁾ = a(x) + t·b(x) + t²c(x)/2 with **arbitrary functions** in all 30 slots. Each building block of the generator is linearised from its nonlinear definition by the script. Units 16πG = 1.

| Row | Term linear in h⁽²⁾ | Result (symbolic residual = 0 in every row) |
|---|---|---|
| R1 | corner Noether charge ∮δQ_ζ | −2δA[q₂]. Before restriction to S it contains 13 function atoms (h_tt, shift and derivative terms among them); on S only the six q₂_ij survive |
| R2 | corner potential ∮i_ζθ | 0 on S |
| R3 | bulk non-Killing term ∫ω(η; h⁽²⁾, £_ζη) | density exactly −2q₂_kk, so −4δV[q₂]; no ∂_t h⁽²⁾, lapse or shift |
| R4 | bulk constraint ∫ζ^a δC_a | −∫N R_lin[q₂]: weight exactly N; the momentum constraint's weight ζ^i\|_Σ is 0 |
| R5, R6 | geometric δA and −2δV (field-dependent ξ) | δA and −2δV |
| R7 | radial relocation of the corner | coefficient ∂_R(A − 2V) = 0 |
| R8 | time relocation of corner and slice | coefficient 0 (the induced metric has no first-order term; the Komar flux of ζ through the cylinder r = 1 is 0) |
| R9 | h⁽²⁾ inside the relative entropy | cited (CLR: quadratic form), not computed |

- **Gates.** The Iyer–Wald identity holds pointwise in symbolic h⁽²⁾ for exactly one sign assignment, and the Green identity holds pointwise as a divergence with flux q_nn − q_kk on S.
- **Totals.** Both constructions, Z1L0 (fixed ζ, coordinate corner) and Z2L1 (field-dependent ξ, relational corner), give (a, b) = (−2, 4): **k = 2**, total = ∫N R_lin[q₂]. The same holds on all 600 probes (30 slots × 20 monomials of degree ≤ 3).
- **On the scope family**, in three second-order gauges each (conformal; plus a gauge term not vanishing on S and junk in h_tt, h_ti, ∂_t h; plus a harmonic conformal and a source-free TT term), the total equals A_l every time:

| Member | A_l = ∫N\|K\|² | A_l / ∮\|K\|² | Total in G_a, G_b, G_c |
|---|---|---|---|
| l = 2, m = 2 | 16π/35 | 1/35 | 16π/35 ×3 |
| l = 2, m = 0 | 12π/35 | 1/35 | 12π/35 ×3 |
| l = 3, m = 3 | 128π/105 | 1/63 | 128π/105 ×3 |
| l = 3, m = 0 | 16π/21 | 1/63 | 16π/21 ×3 |
| l = 4, m = 4 | 512π/231 | 1/99 | 512π/231 ×3 |
| l = 4, m = 0 | 40π/33 | 1/99 | 40π/33 ×3 |

The individual rows move with the gauge (l = 2, m = 2: R1 = 4π/5, −56π/45, −296π/45), only the total does not.

**O3, the mixed reading Z1L1 (council 3d, now computed).** With the fixed ζ and a relational corner the generator is not stationary under a radial shift: ∂_R = 4π in "8πG·H − 2V" units, as q0edge derived by hand. The extra term is (1/8π)·(−8π)·∮δH[q₂] and equals **−4δA + 10δV** exactly on all probes. So Z1L1's total is (−6, 14): an **effective k = 7/3**, not a smeared constraint. It is gauge-dependent (4π/5 on a 10% dilation) and its size on the scope family in the conformal gauge is (21/8, 15/4, 77/16)·A_l at l = 2, 3, 4, with other values in other gauges. This is the registered modifier, not the kill: Z1L0 and Z2L1 are clean, and truthcond's scope already says "the unshifted coordinate ball". The consequence is a wording rule: clause (3), P0 and Q0 are stated with Z1L0 or Z2L1, never with a fixed ζ on a moved corner.

### Part 2: whole-clause check by route EIN2 (`p2_solve.py` exact, `p2_geo.py` float64)

- **Exact vacuum solutions.** For five members and one dressed member, a polynomial h⁽²⁾(t, x) is built such that the **full nonlinear 4D Ricci tensor vanishes at O(ε) and O(ε²), all ten components**, exactly (`out/p2_solve.json`, 23 checks). The record had used only the Hamiltonian constraint.
- **Modular side.** The covariant symplectic current gives E_can = ½|W(h⁽¹⁾, £_ζh⁽¹⁾)| = A_l exactly for all five members.
- **Area side at finite ε.** The stationary point of A − 2V is located numerically on g_ε, generic bases, at ten values of ε per member (`out/p2_geo.json`, 19 checks):

| Member | Δ₂ / A_l | Change, raised basis | Change, raised quadrature | O(ε) coefficient |
|---|---|---|---|---|
| M1, l = 2, m = 2 | 1 − 2.2e−11 | 1.4e−12 | 6.2e−11 | 4e−13 |
| M2, l = 3, m = 3 | 1 − 2.2e−10 | 6.1e−12 | 1.2e−10 | −3e−13 |
| M3, l = 4, m = 4 | 1 − 2.4e−10 | 2.1e−12 | 1.6e−10 | −7e−13 |
| M4, l = 2, m = 0 | 1 − 2.5e−11 | 3.7e−12 | 5.6e−11 | −4e−14 |
| M5, l = 2, mixed m | 1 − 2.6e−11 | 3.6e−12 | 6.0e−11 | 7e−13 |
| M6, M1 with a t-odd term and a second-order gauge term | 1 + 5.1e−10 | 2.8e−11 | 6.1e−11 | −7e−10 |

The registered tolerance was 1e−4, with 1e−5 for the convergence changes. The corner does move (|Y| up to 2.7e−4 at the largest ε, 6.6e−4 for M6), and relocating it changes F only at O(ε⁴) (7e−7), as stationarity requires. So area side = E_can = A_l, and against clause (1)'s S_l = ((l+1)/(l+2))·A_l (cited) the miss is A_l/(l+2).

**Zero modes.** The flat Hessian has exactly seven (largest |eigenvalue| 7e−6 against a gap of 0.54). On M1–M5 the gradient along them is zero (≤ 3e−8 per unit ε²), as argued in the PREREG. On M6 it is not: it scales as ε³ (2.2e−4, 4.5e−4, 6.8e−4 per unit ε² at ε = 0.1, 0.2, 0.3), from the t-odd term. The locator projects these directions out, so on M6 the corner's Poincaré position was not relaxed. Its effect on F is O(ε⁴) and does not reach Δ₂.

## 2. The plants (run first)

| Plant | Behaviour |
|---|---|
| **PLANT-KLS**, +½∮f(q_kk − q_nn), f the l = 2 zonal and sectoral harmonics | Caught three ways at both f: the audit finds 15 and 12 probes with residual; nonzero on a pure-gauge corner displacement (−4π/5, −16π/15); the on-shell total is gauge-dependent at every member (l = 3, m = 3, f zonal: 136π/105, 1706π/1575, −274π/1575 against A = 128π/105) |
| **PLANT-CURV**, +∮f n·R^lin[q₂]·n | **Zero on all three pure-gauge probes**, so the gauge channel alone misses it. Caught by the audit (9 and 6 probes) and by the source-free TT term, which moves the total at every member |
| **PLANT-LAPSE**, N → N + (1 − r²)² | Caught: the constraint row is no longer −𝒞, the integrated identity fails, and the on-shell value differs at every member (272π/315 against 16π/35 at l = 2) |
| **PLANT-FLAT**, exact flat metric with the corner displaced at O(ε) | The locator finds T = −εf (leading coefficient 0.5000000) and \|F* − F₀\| = 2e−15 at ε = 0.3. The fixed coordinate corner is off by −4.7e−2 |
| **PLANT-SRC**, constraint source × 11/10 | Rejected by the exact Ricci check; Δ₂/A_l = 1.1 − 3e−11 |

One thing the plants taught: at l = 2, m = 2 in the conformal gauge alone, the sectoral KLS plant and the sectoral curvature plant both return exactly A_l. A single-gauge on-shell check would have passed them. The three-gauge requirement did real work.

## 3. Would a wrong theory also pass?

- **Part 1: yes, for any world whose modular side is the Einstein–Hilbert Iyer–Wald generator.** There the result is a theorem, and the script is a mechanical proof of it with the plants as evidence that the audit is not blind. It does not reach a term outside these constructions.
- **Part 2: yes, for any world in which the area side is not the geometric (A − kV)/4G.** It would fail for a wrong k, a mislocated corner, a K-quadratic term inside A or V, or an h⁽²⁾ that does not solve the equations.

## 4. The weakest step

**Completeness of the enumeration is relative to the construction.** "Every term" means every term of the Iyer–Wald generator in the readings Z1L0 and Z2L1. Three things are outside it and are named, not tested:

1. The corner-ambiguity terms of Harlow–Wu, Chandrasekaran–Speranza and Freidel at h⁽²⁾-linear order. q0edge covered their K² sector only. With a fixed ζ, £_ζh⁽²⁾ does not vanish on S (it is the boost variation of h⁽²⁾), so a term such as β(η; £_ζh⁽²⁾) could be nonzero under Z1; under Z2 it vanishes. Not computed here.
2. The operator-ordered, state-dependent piece of JSS (3.11). No classical computation touches it.
3. R9 is cited.

## 5. What a hostile referee attacks

1. **"It is an identity, not a test."** Correct for Part 1 inside its construction. Its content is the mechanical enumeration, the plants, and the Z1L1 finding.
2. **"Part 2 assumes its conclusion."** See §6.
3. **"k = 7/3 shows the clause depends on a reading."** It does, and the report says which readings are clean. A referee will ask the papers to state Z2 (or the coordinate ball) wherever P0 or Q0 appear.
4. **"The relational corner is not unique."** The seven zero modes mean the corner is located only up to Poincaré motions at this order. On single-l data that is harmless by symmetry (measured: zero). On general data, or M6, a nonzero O(ε³) force exists; it does not reach the O(ε²) area side, but a referee can ask for the next order.
5. **"The audit's probes stop at degree 3."** The symbolic residuals cover arbitrary h⁽²⁾; the probe space is only the common pipeline that plants and real rows share.
6. **Coverage.** Part 2 ran m = l only at l = 3, 4 (plus m = 0 and mixed m at l = 2).

## 6. Is Part 2 genuinely independent of P0, T, U and Q0?

**As inputs, yes. As a hypothesis about what the area side is, no.**

- **T, U, FL1, Q0:** independent. No reflection symmetry is used (M6 is not time-symmetric), there is no matter, the first-order law is an output (c₁ = 0 to 1e−9), and K is kept to all orders at a corner the computation finds itself.
- **P0:** the Taylor transfer is not used; the number is a finite-ε difference. But the route assumes **GEO**: the area side is the geometric (A − 2V)/4G of the extremal corner. GEO is the *conclusion* of the B theorem's Theorem 1, taken from print (JSS §3.3, BMSV) instead of derived, and it is a stronger statement than P0's "some metric functional". EIN2 therefore checks the clause's second sentence ("the area side is A") given GEO and the Einstein equations. It does not check, by a second route, that the gravitating modular side has no other h⁽²⁾ term. A pricer should treat it as a clause-wide computation conditional on GEO, not as an unconditional one.

## 7. Deviations from the PREREG

- **TT probe (plant stage, before any unplanted result was read).** The first source-free TT probe had zero overlap with the sectoral curvature plant, so one of that plant's three required behaviours failed while the audit still caught it. The probe's weights were changed once and a check that it overlaps both profiles was added.
- **One index slip after the unplanted run.** The Komar charge of the fixed ζ was summed over the wrong index and read 0; two bookkeeping checks failed and exposed it. Corrected: −8πR³. It changed O3's coefficient only (to the registered expectation 4π). O1, O2 and the plants were unaffected.
- **Fit.** Ten points (five magnitudes, both signs) with a degree-6 polynomial including odd powers, in place of an even fit over five values, so that the O(ε) coefficient could be measured.
- **Numerical settings not registered:** base L_Y = 4, d_T = 3, quadrature (20, 40, 16); raised to L_Y = 6, d_T = 4 and to (26, 52, 22); Newton stops at a projected gradient of 2e−9, the finite-difference noise floor. One run crashed on a quoting slip in the script and was rerun unchanged.
- **Not enumerated:** the corner-ambiguity terms (§4.1), as the PREREG said.

## 8. Files

All under `<local path>`: `PREREG.md`, `REPORT.md`, `HASHES.sha256`, `scripts/` (`common.py`, `p1_enum.py`, `p2_solve.py`, `p2_geo.py`, `run_all.sh`), `out/` (`p1_plants.json`, `p1_enum.json`, `p2_solve.json`, `p2_metrics.json`, `p2_plants.json`, `p2_geo.json`). Rerun: `scripts/run_all.sh` (plants first; exact parts about 10 s, numerics about 15 minutes, under 1 GB). Python 3.12, sympy 1.14, numpy 2.3.

**Recommended next lane: AXIOM-PRINT on clause C.** Read Uhlmann monotonicity and Haag duality for the helicity-2 curvature net first-hand, quote each with a locator in the clause's scope, and search equally hard for contrary print on A0; in the same lane, do the same for B's premises T, U and FL1 and for GEO. Under RUBRIC12 a conditional derivation is worth exactly 1 until its premises carry outsider-verifiable grades, the top-ranked upgrade is a clause-wide item on C (+8 to +11 points), and this lane has just shown that B's remaining weight sits in premises and constructions that only print can ground, not in any further computation.
