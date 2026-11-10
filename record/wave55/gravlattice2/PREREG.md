# PREREG: w55-gravlattice2. A decisive rerun of the graviton lattice test of Theorem 1 (no surface term), rule B fixed in advance

Lane w55-gravlattice2 (Opus), Mon 9 Nov 2026, written 09:30 Europe/Oslo, before any new lattice row.
Referee change #1 (wave54/rereadA11/REVIEW.md §5): pre-register rule B, reach R = 8–10 with a symmetry-reduced (non-radial) solver,
use lower-degree controls, then let the |c| < 0.15 verdict pass or fail.
First run: wave52/gravlattice (REPORT.md: INCONCLUSIVE; rule B R ≤ 7; controls 0.82–0.84 at R = 7; |c| ≲ 0.2 unresolved).

Already computed before this file (theory side only, no lattice number):
- `scripts/ctrl_solve.py` → `out/controls.json`, `out/ctrl_solve.log`: lower-degree controls (below).
- `scripts/pred2.py` → `out/pred.json`, `out/pred.log`: S_T1, σ and the plant ratios for all eleven data.
- `scripts/glat.py`, `maxlat_core.py`, `pred.py`: copied byte-identical from wave52/gravlattice/scripts (sha256 equal).

## 1. Claim under test (unchanged from wave52)
Theorem 1(c), Paper A: for smooth constraint data near the ball with no l ≤ 1 part, including data whose chain fields do not vanish on the
sphere, the coherent-state relative entropy on the curvature net of the ball is S = Σ c_l Q[u], c_l = 4/(L(L−2)), with no surface term.
Lattice target S_T1 = 2π Σ_l c_l [∫_B N_ρ(ψ̇_l² + |∇ψ_l|²) + ρ⁻¹∫_B ψ_l²] in canonical units (pred.py, unchanged). Prediction: every
ratio S_lat/S_T1 → 1, and the surface coefficient c → 0.

## 2. What is fixed in advance
- **Region rule B is the primary and only verdict rule:** the generators e = inc(h) (primal cells) and b = sym curl π (dual cells) whose whole
  stencil lies on primal cells inside |x − c| < R. Rule A is not run.
- **Lattice, correlators, Gaussian construction: unchanged** (glat.py). ε cap 30 primary; 25/35 for the cap spread.
- **Centres:** site (0,0,0) and cube (½,½,½); every ratio is the mean of the two; each centring is also reported alone.
- **Radii:** R = 3–10, both centrings. R = 11 and 12 are added if memory and time allow; if added, they enter the fits by the window rule below
  (Rmax is the largest R with both centrings).
- **Prediction convention:** "vol" (S_T1 at the site-volume radius R_eff) primary; "raw" (ρ = 1) enters the bracket.

## 3. Memory: the solver (non-radial)
- The ball, the staggered lattice, the vacuum and the stencils are invariant under the three coordinate reflections (about the centre),
  the group Z2³ ⊂ O_h. Each reflection acts on the generators as a signed permutation, computed numerically from the stencils
  (h, π true tensors; derived signs for e and b). Block-diagonalise X = ⟨ee⟩, P = ⟨bb⟩ and W = [e, b]/i into the eight parity sectors.
  This is a cubic-symmetry decomposition of the 3D cubic ball, not a partial-wave chain: every sector mixes all l and both chain parities.
- Only rows X[orbit representatives, all] are built (about n²/8 entries); W stays sparse; each sector runs glat's rel_entropy with the
  null and rank thresholds taken relative to the global (all-sector) maxima, so the sum over sectors is the same computation.
- **Verification gate:** the new solver must reproduce all ten wave52 rule-B rows (R = 3–7, s and c) for every datum and every cap to
  relative 1e-8 (absolute 1e-10 where S < 1e-2) before any row with R ≥ 8 is computed. If it fails after one fix round, the lane stops:
  verdict BLOCKED, no new rows booked.
- Peak memory kept under 9 GB (measured with /usr/bin/time -v per row).

## 4. Data and controls
- Non-vanishing data (ψ ≠ 0 on S), unchanged: g1 uniform tidal E (l = 2), g2 uniform B (l = 2, odd), g3 l = 3, g4 l = 2 + 3, g5 l = 4,
  g6 = g1 + g2. σ/S_T1 = 0.70–0.74.
- **New lower-degree controls (σ = 0 exactly; ψ = (1 − r²)P with P harmonic):** d1_E2 (Y degree 5, P = P2), d2_B2 (Z degree 5, P = P2b),
  d3_E3 (Y degree 6, P = P3). The first run's controls had degree 11. Old controls c1_E, c2_B are kept as secondary.
- **Control assignment (primary):** g1, g3, g4, g5 → d1_E2; g2 → d2_B2; g6 → the S_T1-weighted mean of d1 and d2 ratios
  (weights 3.175 : 3.008, g6's E and B parts). **Alternate:** g3, g4, g5 → d3_E3; g1 → d3_E3; g2 → c2_B; g6 → c1_E.

## 5. Estimator (fixed)
- Per R: ratio_X(R) = mean over centrings of S_lat/S_T1. Differential surface coefficient
  c(R) = (ratio_datum(R) − ratio_ctrl(R)) / (σ/S_T1)_datum.
- **Primary c_∞:** fit c(R) = a + b/R over R ∈ [Rmax − 4, Rmax] (five radii), vol convention; c_∞ = a.
- **Bracket:** min..max of a over {lin a + b/R on [Rmax−3, Rmax], [Rmax−4, Rmax], [Rmax−5, Rmax]; quad a + b/R + d/R² on [Rmax−5, Rmax],
  [Rmax−6, Rmax]; the unextrapolated c(Rmax)}, in both vol and raw conventions (12 numbers).
- Control convergence: the same estimator applied to ratio_ctrl(R) (lin on [Rmax−4, Rmax]) gives a_ctrl.

## 6. Pre-registered outcomes
- **PASS (prediction).** (i) Controls converge: d1 and d2 have |a_ctrl − 1| ≤ 0.10 and |ratio_ctrl(Rmax) − 1| ≤ 0.15.
  (ii) |c_∞| < 0.15 on all six non-vanishing data (primary controls). (iii) The threshold plants c = +0.15 and c = −0.15 both lie outside
  the bracket on at least four of the six data.
- **KILL gate (a continuum surface term |c| ≥ 0.15).** (i) as above; (ii) c_∞ of one sign with |c_∞| ≥ 0.15 on at least one E-type datum
  and on g2 (two families), with the bracket excluding 0; (iii) present in the site-only and cube-only estimates and with the alternate controls.
- **INCONCLUSIVE.** Anything else, including controls that miss (i), or a bracket wider than the gap to be resolved.
- **c_l test (secondary, differential, as wave52):** [ratio](g3, g4, g5)/[ratio](g1) → 1. Plants: c_l → 1/L predicts 2.50, 2.21, 4.50;
  c_l → 4/L² predicts 1.25, 1.20, 1.35; Born 0 predicts 1.125, 1.101, 1.200. A plant is caught when its value is outside the bracket
  (same bracket rule, applied to the differential ratio minus 1).
- **False-kill rule.** Before any negative is booked: centring (site, cube, mean), ε cap, staircase (vol against raw), primary against
  alternate controls, solver verification, data scaling with R, and that the controls pass.

## 7. Plants
| Plant | Side | Predicted ratio S_plant/S_T1 | Must |
|---|---|---|---|
| 1. Surface (π/2)σ | theory | 2.10–2.17 on g1–g6 | be caught on all six |
| 2. c_l → 1/L | theory | 2.50 / 2.21 / 4.50 on g3 / g4 / g5 | be caught on all three |
| 3. c_l → 4/L² | theory | 0.83 / 0.80 / 0.90 (differential 1.25 / 1.20 / 1.35) | low power, stated in advance |
| 4. Born coefficient 0 | theory | 0.67–0.80 | be caught |
| 5. **New: threshold surface c = ±0.15** | theory | 1 ± 0.105–0.111 | the power check of this test: rejected on ≥ 4 of 6 for PASS |
| 6. **New: no classical centre** (drop the KL term of the boundary centre from S_lat) | lattice | g2 and g6 shift by their KL share (0.69 and 0.35 at R = 7 in wave52); E-only data unchanged | be caught on g2 and g6; no power on g1, g3–g5 |
| 7. Injection c_inj σ (c_inj = 0.15, π/2) into S_lat | lattice | recovered by the estimator | linearity check only |

## 8. Legs from structure and price
- Theorem 1 is re-derived analytically (referee), maxlattice found no surface term for helicity 1, and wave52 excluded a natural-size term
  (|c| ≲ 0.2). The remaining alternative is a small spin-2 boundary contribution, 0.15 ≲ |c| ≲ 0.3 (for example an edge or centre effect
  of the curvature net).
- P(PASS | T1 true) ≈ 0.55: the main risk is INCONCLUSIVE from slow convergence (≈ 0.40); a false kill ≈ 0.05.
- P(PASS | |c| ∈ [0.15, 0.3]) ≈ 0.12 (fit bias could hide it). P(KILL | T1 true) ≈ 0.05; P(KILL | |c| ∈ [0.15, 0.3]) ≈ 0.5.
- LR on PASS: about 4 on the narrow clause "no surface term with |c| ≥ 0.15"; on Theorem 1(c) as a paper result (paper value) about 1.3–1.5;
  on the credence rows M and C at most 1.02, since the theorem is already checked analytically. LR on a KILL surviving the false-kill rule:
  about 0.1 on the clause, 0.3 on the paper value. INCONCLUSIVE: LR ≈ 1.
- Eureka outcome: none for a headline claim. A KILL would be a new spin-2 boundary term in modular data and would move Theorem 1 by more
  than 20 points; its honest probability is about 0.05.
