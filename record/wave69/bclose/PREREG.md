# PREREG: w69-bclose. A route-negating enumeration for clause (3) / B, and a whole-clause check that uses none of P0, T, U, Q0

Lane w69-bclose (Fable 5.2), written before any script of this lane exists or runs. Claude only. Spec: `wave68/council2/COUNCIL.md` items 3a and "1. B-CLOSE". No prices or likelihoods are set here; a separate pricer sets them.

Read before writing: GOAL.md; truthcond §2 (clause (3) and §2.3); the council spec; `wave64/btheorem/NOTE.md`; `wave64/blindB/VERDICT.md`; `wave68/q0edge/REPORT.md`; `wave66/p0test/REPORT.md` and its `s1_transfer.py`; `wave39/BR/REPORT.md` (for what "KLS-type" means in the record: the angle-weighted charges δ²A(f), δ²𝒬(f) with f ≠ 1).

## 0. The clause and its failure route

> **(3)** Every term of the second-order modular side that is linear in the second-order metric is the smeared constraint δ(A − kV). So the area side of the second-order first law is A.

Failure route (truthcond §2.3): *a second-order modular term linear in h⁽²⁾ that is not a smeared constraint (for example, KLS-type angle-dependent charges in the ball generator).*

**Setting.** Flat unit ball B in Σ = {t = 0}, S = ∂B, ζ = ½(1 − t² − r²)∂_t − t x·∂, N = (1 − r²)/2. Scope family g_ε = η + ε·2tK + ε²h⁽²⁾, K = sym∇(x × ∇φ), φ a real solid harmonic of degree l. Members: l = 2, 3, 4, each with m = l (φ = Re(x + iy)^l) and m = 0 (zonal). q₂ = h⁽²⁾_ij|_Σ.

**Units.** 16πG = 1 for the enumeration; the totals are quoted in Paper A's units, where the area side reads A_l = ∫_B N|K|². Fixed functionals: δA[q] = ½∮_S(q_kk − q_nn), δV[q] = ½∫_B q_kk, R_lin[q] = ∂_i∂_j q_ij − ∇²q_kk, and the smeared constraint 𝒞[q] := ∫_B N R_lin[q].

**What "is the smeared constraint δ(A − kV)" means operationally (fixed now).** A linear functional T of the second-order data (all ten components of h⁽²⁾ and their derivatives on Σ) passes iff, as an identity in *symbolic* h⁽²⁾ (arbitrary functions, not sampled data):

- (i) T = a·δA[q₂] + b·δV[q₂] + c·𝒞[q₂] with constants a, b, c (no angular weight, no dependence on h_tt, h_ti, ∂_t h, corner jets or curvature), and
- (ii) the construction's total has (a, b) ∝ (1, −2), so that by the Green identity 𝒞 = −2δ(A − 2V) (verified pointwise as a divergence identity in the same script) the total is c_tot·𝒞.

The **residual** of a term is T minus its (a, b, c) part. A nonzero residual is a term linear in h⁽²⁾ that is not a smeared constraint.

## 1. PART 1: the enumeration

The modular side's h⁽²⁾-linear part is the linearisation at η of the generator as it is built. The script expands each building block with symbolic h⁽²⁾ and lists every resulting term. Constructions covered (the ones the record uses):

- **Z1L0**: the Iyer–Wald / JSS (3.11) generator of the fixed coordinate vector ζ at the coordinate corner (w66-p0test S1).
- **Z2L1**: the field-dependent ξ of JSS (constant unit surface gravity in g, attached to the surface), with the corner located by extremising A − kV (JSS §3.3; blindB's repair). Its charge is the geometric A[S] − 2V[Σ_S].
- **Z1L1**: fixed ζ with the relational corner. Registered as an *open observable*, not as part of the kill (see O3).

Rows (each computed from the nonlinear object, then linearised by the script):

| Row | Object | Expected form (to be asserted, can fail) |
|---|---|---|
| R1 | corner Noether charge ∮_S δQ_ζ[h⁽²⁾] | a·δA, a constant; no h_tt, h_tr, h_rr, shift or derivative dependence |
| R2 | corner potential ∮_S i_ζθ[h⁽²⁾] | 0 |
| R3 | bulk non-Killing term ∫_B ω(η; h⁽²⁾, £_ζη) | b·δV; no dependence on ∂_t h⁽²⁾ (the second-order momentum), lapse or shift |
| R4 | bulk constraint term ∫_B ζ^a δC_a[h⁽²⁾] | Hamiltonian part = c·𝒞 with weight exactly N; momentum part has weight 0 |
| R5 | geometric area δA_geom[h⁽²⁾] of the coordinate sphere (Z2) | δA |
| R6 | geometric volume −k δV_geom[h⁽²⁾] of the t = 0 slice (Z2), k = 2 | −2δV |
| R7 | radial relocation of the corner: (∂F₀/∂R)·Ȳ₂[h⁽²⁾], with Y₂ from the linearised extremality condition | coefficient ∂F₀/∂R = 0 for Z2 (F₀ = A − 2V) |
| R8 | time relocation of corner and slice, X₂[h⁽²⁾], T₂[h⁽²⁾] | coefficient 0 (first variation of F₀ in time vanishes) |
| R9 | h⁽²⁾ inside the modular-side relative entropy | absent: the free-net relative entropy is quadratic in the first-order datum (cited: CLR; not computed here) |

Gates (conventions only, not evidence): the off-shell Iyer–Wald identity R3 = R1 − R2 + R4 and the Green identity must hold pointwise in symbolic h⁽²⁾. They fix orientation signs. If a gate fails, the run is INVALID until the convention error is found (one fix round).

**Observables.**

- **O1 (the kill observable).** The residual of each of R1–R8 in Z1L0 and Z2L1, as a symbolic expression. Prediction: all zero.
- **O2.** The totals: (a, b) of Z1L0 (R1, R2, R3) and of Z2L1 (R5–R8); k = −b/a; on the six scope members, with q₂ solving R_lin[q₂] = |K|², the value of each row and of the total in three second-order gauges (G_a: q₂ = uδ; G_b: G_a plus sym∇ξ with ξ polynomial, not vanishing on S and with a monopole radial part; G_c: G_b plus a harmonic conformal term and a source-free TT term). Prediction: k = 2; total = A_l in every gauge at l = 2, 3, 4, both m; A_l/∮|K|² = 1/((2l+1)(2l+3)).
- **O3 (open; council 3d, q0edge §4).** For Z1L1, R7's coefficient is ∂_R(Komar_ζ − 2V) at R = 1. Compute it, compute Ȳ₂[q₂] = (1/8π)∮δH[q₂] from the linearised mean curvature, and report the extra term on the six members in the three gauges. Prediction (q0edge's hand derivation): coefficient 4π in "8πG·H − 2V" units, and the term is gauge-dependent, hence not a smeared constraint.

**Pass and fail criteria, Part 1.**

- **P1-PASS**: O1 all zero in Z1L0 and Z2L1; gates hold; k = 2; totals equal A_l in every gauge at every member; both plants caught.
- **P1-KILL**: a nonzero residual in Z1L0 or Z2L1 that survives the false-kill checklist (§4). This is the registered kill: *a term linear in h⁽²⁾ that is not a smeared constraint.*
- **P1-INVALID**: a gate fails after one fix round, or a plant is not caught.
- **O3 modifier**: **Z1L1-CLEAN** (the extra term vanishes or is a smeared constraint) or **Z1L1-DIRTY** (it is not). DIRTY is *not* the kill, for a reason fixed now: Z1L1 pairs a coordinate-fixed vector with a surface that moves under diffeomorphisms, so it is not a diffeomorphism-invariant definition of a generator, and the record (blindB, q0edge) already restricts Q0 and P0 to Z2. DIRTY is reported as a scope condition on clause (3): "stated with the field-dependent ξ". If O3 comes out CLEAN, the record's hand claim is wrong and the report says so.

**Plants (run first, before any unplanted result is read).**

- **PLANT-KLS** (the spec's): add +1·δA(f)[q₂] = ½∮_S f(n)(q_kk − q_nn) to R1, with f = r²P₂(cos θ) restricted to S (l = 2, m = 0), and again with f the l = 2, m = 2 harmonic. This is the angle-weighted area charge of w39-BR S3 and KLS II (δ²𝒬(f), f ≠ 1). Required: nonzero symbolic residual; nonzero value on a pure-gauge corner displacement; the on-shell total is no longer the same in gauges G_a, G_b, G_c, at each of l = 2, 3, 4.
- **PLANT-CURV** (mine): add +1·∮_S f(n) n_i n_j R^lin_ij[q₂] with the same f, where R^lin_ij is the linearised Ricci tensor of q₂. It is invariant under every smooth gauge transformation, so the gauge channel cannot see it. Required: nonzero symbolic residual; *zero* on the pure-gauge probe (showing the gauge channel alone would have missed it); the on-shell total changes when the source-free TT term of G_c is added, at each l.
- **PLANT-LAPSE** (third, cheap, covers the record's "lapse weight" route): replace N by N + (1 − r²)² in R4. Required: the Iyer–Wald gate fails and the on-shell total differs from A_l.

## 2. PART 2: whole-clause check by route EIN2

**Route.** Compute the second-order area side directly, at finite ε, from the geometry of the scope spacetime, and compare it with the modular side's quadratic energy computed from the symplectic form.

1. **Solve.** For each member, construct an exact polynomial h⁽²⁾(t, x) such that the full nonlinear 4D Ricci tensor of g_ε = η + ε·2tK + ε²h⁽²⁾ vanishes through O(ε²), all ten components, checked exactly (sympy). The record so far used only the Hamiltonian constraint.
2. **Area side, non-perturbatively.** On the metric g_ε at finite ε, find numerically the stationary point of F[Σ] = A[∂Σ] − 2V[Σ] over spacelike discs Σ (interior stationarity is maximality; boundary stationarity is JSS's extremal corner). Parametrise the corner as r = 1 + Y(n), the slice as t = T(x), with Y in real harmonics up to degree L_Y and T in polynomials up to degree d_T, generic (no symmetry imposed). Newton iteration from the coordinate location; the seven Poincaré zero modes of the flat Hessian are projected out and the gradient along them is reported. Define Δ(ε) = −2[F*(ε) − F₀]/ε² and Δ₂ = lim_{ε→0} Δ(ε), from an even-polynomial fit over five values of ε. No formula for δA[q₂], no assumption that the corner stays put, no assumption about K-quadratic corner terms.
3. **Modular side.** E_can = ½|W_Σ(η; h⁽¹⁾, £_ζh⁽¹⁾)| from the covariant Einstein–Hilbert symplectic current on the 4D first-order field, exact (sympy). It involves h⁽¹⁾ only. The relative entropy S_l = ((l+1)/(l+2))·A_l is clause (1), cited, not recomputed.
4. **Compare.** Δ₂ against E_can and against A_l = ∫N|K|² from Part 1.

**Members.** Magnetic: l = 2, 3, 4 with m = l; l = 2 with m = 0; one mixed-m real l = 2 multipole; one m = l, l = 2 member whose h⁽²⁾ carries an extra second-order gauge term with radial, tangential and time parts not vanishing on S (so the corner really moves at O(ε²)) and an extra t-odd homogeneous term (so nothing is time-symmetric).

**Why it is independent of P0, T, U, Q0 (and FL1).**

- **P0** (one Taylor functional; h⁽²⁾ enters only through the first-order functional): not used. Nothing is Taylor-transferred; no h⁽²⁾-linear part is extracted; the number is a finite-ε difference of areas and volumes.
- **T** (time reflection): not used. The bases are generic and one member is deliberately not time-symmetric.
- **U** (matter universality) and **FL1** (first-order first law for matter): not used. There is no matter. The first-order law appears as an output (no O(ε) term), not an input.
- **Q0** (no K-quadratic geometric term): not used. All orders in K are kept at finite ε, at a corner the computation locates itself.
- **What it does assume, stated now so it cannot be discovered later:** (EH2) the vacuum Einstein equations through O(ε²); (GEO) the area side of the first law is the geometric quantity (A − kV)/4G with k = 2, of the corner located by extremising it. GEO is the *conclusion* of the B theorem's Theorem 1, taken here from print (JSS §3.3, BMSV) instead of derived from P0, T, U, FL1. So EIN2 is independent of the four premises as inputs, and is **not** independent of the possibility that the true gravitating area side is not GEO. The report must say this in the same words.

**Pass and fail criteria, Part 2** (fixed tolerances):

- **P2-PASS**: for every member |Δ₂/A_l − 1| ≤ 1e−4, with Δ₂ changing by ≤ 1e−5·A_l when L_Y, d_T and the quadrature order are each raised one step; E_can = A_l exactly; the O(ε) coefficient of F* is zero within the fit; both Part 2 plants behave as required.
- **P2-FAIL**: some member has |Δ₂/A_l − 1| > 1e−2 with the convergence criterion met, surviving §4.
- **P2-INCONCLUSIVE**: anything between, or the numerics do not converge or do not finish within the lane.

**Part 2 plants.**

- **PLANT-FLAT**: g = ψ*η for an exact polynomial diffeomorphism ψ: t → t + εf(x), which displaces the true corner at O(ε). Required: the located stationary value satisfies |F* − F₀| ≤ 1e−9 at the largest ε, while F at the fixed coordinate corner and slice differs from F₀ by more than 1e−5. This shows the locator finds a displaced corner and that a fixed corner is caught.
- **PLANT-SRC**: h⁽²⁾ built with the constraint source scaled by 11/10 (not a solution). Required: Δ₂/A_l = 1.1 within 1e−4, i.e. the comparison is sensitive to the O(ε²) Einstein equations at the 10% level it was broken.

## 3. "Would a wrong theory also pass?" (one guard per outcome)

- **P1-PASS.** Within Einstein–Hilbert covariant phase space the Iyer–Wald identity is a theorem, so a pass is expected in every world in which the modular side is the Z1L0 or Z2L1 generator. The pass has content only as (a) a mechanical enumeration, with symbolic h⁽²⁾, that no building block hides an angle-weighted, lapse-, shift-, momentum- or jet-dependent piece, and (b) the plants. A world where B fails because the true generator contains a term *outside these constructions* (an unprinted edge completion; the operator-ordered, state-dependent piece of JSS (3.11); the corner-ambiguity terms of Harlow–Wu, Chandrasekaran–Speranza and Freidel at h⁽²⁾-linear order, which this lane does not enumerate) **passes too**. The report must list these as outside the enumeration.
- **P1-KILL.** A wrong kill would come from a convention slip (orientation, the factor 2 between the Taylor coefficient and the linear functional, the sign of K, a term dropped because ζ|_S = 0 before differentiating). §4 lists the checks; a kill is provisional until they pass.
- **Z1L1-DIRTY.** Would a correct construction also show it? Yes, by design: it is a property of the mixed reading, not of B.
- **P2-PASS.** A world where the area side is not GEO passes. A world with a K-quadratic term inside A or V themselves, a mislocated corner, a wrong k, or a wrong h⁽²⁾ does not: PLANT-FLAT and PLANT-SRC measure the last two, and k ≠ 2 makes F non-stationary under radial shifts, which the located corner would expose as a Δ₂ that depends on the second-order gauge member.
- **P2-FAIL.** A wrong fail would come from non-convergence, a zero-mode shift of O(1) (see below), an h⁽²⁾ that does not solve the equations, or ε too large for the fit. Each is checked before a fail is believed.
- **Zero modes (fixed now).** The flat Hessian of F has seven zero modes (time translation, three translations, three boosts). If the O(ε²) gradient along one were nonzero, the relational corner would shift at O(1) and Δ₂ would differ from the fixed-corner value. Expectation, argued before running: on single-l scope data the gradient vanishes, by (ε → −ε, t → −t) for the time translation and boosts and by point reflection (with or without t → −t) for translations. The script reports the gradient along the zero modes divided by ε²; a nonzero value is a finding and is reported whichever way the bin falls.

## 4. False-kill checklist (run before any KILL or FAIL is reported)

Orientation of S and Σ; the sign convention of K against h⁽¹⁾ = 2tK; the factor between the O(ε²) Taylor coefficient and the linearisation at h⁽²⁾; whether ζ|_S = 0 was used before or after differentiating; exact integration (closed monomial formulas, no `nsimplify`); completeness of the harmonic and polynomial bases; the units map between 16πG = 1 and Paper A's A; that h⁽²⁾ passes the exact Ricci check; quadrature and basis convergence; ε range.

## 5. Outcome bins

| Bin | Condition |
|---|---|
| **BOTH-PASS** | P1-PASS and P2-PASS |
| **P1-ONLY** | P1-PASS; P2-INCONCLUSIVE |
| **P2-ONLY** | P2-PASS; P1-INVALID |
| **KILL** | P1-KILL, or P2-FAIL, surviving §4 |
| **INVALID** | P1-INVALID and P2 not PASS |

Each bin is reported with the O3 modifier (Z1L1-CLEAN or Z1L1-DIRTY) and with the zero-mode gradients.

## 6. Files and rerun

`scripts/common.py`, `scripts/p1_enum.py`, `scripts/p2_solve.py`, `scripts/p2_geo.py`, `scripts/run_all.sh` (one command: plants first, then the unplanted runs); `out/*.json`; `HASHES.sha256`; `REPORT.md`. Reused: the monomial integration helpers of `wave64/btheorem/scripts/common.py` and `wave68/q0edge/qcommon.py`, and the P-tensor and charge algebra of `wave66/p0test/scripts/s1_transfer.py`.

Stop: REPORT.md exists, or one fix round on a blocker. Fences: nothing outside `wave69/bclose`; no commit; no edit to CALIBRATION.md, LANES.md or the papers.
