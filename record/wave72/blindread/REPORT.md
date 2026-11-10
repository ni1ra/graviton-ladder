# Blind read of the short paper (version 1), with the lead's grading

10 November 2026. One reader (Claude Fable), given only `ladder_short.tex` and its PDF, forbidden to open the scripts, the reports or the record, told to follow the argument as the intended reader and then to break it. It wrote its own scripts from the paper's text. Web access was limited to the reference list and published formulas.

## Verdict

**FOLLOWS WITH FIXES.** Every formula and number the reader could recompute agrees with the paper. Four statements were wrong or broader than the body supports (E1 to E4), and the entropy side rests on fact (b), which the reader could not see.

## Errors found

**E1. Section 4: "returns exactly zero for l = 1 (which is flat space in accelerated coordinates)".** The Weyl metric with U = εz (k = −ε²ρ²/2) is Ricci-flat but not flat: its Kretschmann scalar is 16ε⁴(ε²ρ² + 3)e^{2ε(ερ² + 2z)}. The zero is right for another reason: an O(ε²) vacuum curvature enters F linearly through the Green identity, where R_lin = 0.

**E2. Abstract, eq. (2), "no cross term between the two sectors: it would be odd under parity".** True for the entropy. Not enough for the area side: with E and B both present the stationary ball need not exist at O(ε²). A bilinear E⊗B vector at equal l is odd under parity and under time reversal, which is a boost push. For l = 2, under a boost β of the ball, d/dβ[E² + (4/3)B²] = −(28/3)(β·P), with P_k = ε_kil E_ij B_lj. So the reduced F has a gradient along boosts unless the super-Poynting vector vanishes ([E, B] = 0). The reader obtained this by reduction onto the zero modes, not by computing the force directly.

**E3. Section 1: "For matter in a fixed diamond the two are equal at first order".** Contradicts section 5, where the minimal scalar's entropy exceeds the area side by π∫φ². Should read "for conformal matter".

**E4. Sections 5 and 8 misdescribe references [8] and [2].** Eq. (75) of Jacobson, Senovilla and Speranza (arXiv:1710.07379) is δA|_V = Ω_{d−2}ℓ^{d+2}/[3(d+3)(d²−1)]·[−W + E²(3d(d+3)γ² − 2(d+3)γ + 5/6)], W = E² + B² in d = 4. At the constant-mean-curvature shape γ = 1/(3d) it gives δA|_V = −(πℓ⁶/105)(E² + (4/3)B²): exactly −4G times the paper's area term. So the l = 2 area side is in print. "The prescription of [2]" is a geodesic ball at fixed volume (γ = 0), which gives −(4πℓ⁶/315)(E²/6 + B²); the electric case then fails too. The B² coefficient is shape-independent at fixed volume; the E² coefficient is not.

**E5. Minor.** A missing ε in the cut velocity of 3(c). The symplectic product in 3(a) is half the canonical one and should be labelled so.

## Unclear or unsupported places

1. 3(c), "the profile 2NK vanishes on the sphere" carries the result and was not justified. With only the interior curvature and the constraint across the sphere, the pairing is not unique (35 free gauge parameters, pairing depends on them). Requiring h = 0 on the sphere makes it unique, at 2(l+1)/(l+2) times A (checked l = 2, 3).
2. "The tidal field is a coherent state": a polynomially growing solution has infinite norm; say that only data near the ball matter.
3. "Entanglement (relative) entropy": S is relative entropy; for a coherent state it is the change in modular energy.
4. "Other m follow from rotation invariance": needs Schur's lemma on the irreducible multiplet.
5. "At every multipole order l" for the electric area side: the body showed l = 2 to 6. The reader confirmed l ≤ 12 and the closed form F₂(a) = 2π(l−1)(l+2)/(2l+1)·a² − 8π(l−1)/(2l+1)·a + 12π(l−1)/((2l+1)(2l+3)).
6. "Not a scalar": one line covers all l (every parity-even quadratic Lorentz scalar flips sign under duality). Pre-empt the rewrite: for l = 2, area side = (7/6)S − (πR⁶/2520G)(E² − B²). A scalar gives proportionality at best, never equality.
7. Electric data: state that the dual of the magnetic tide is the electric tide with the same ψ. The reader confirmed the value without duality at l = 2, 3.
8. Cite the scalar-mode reduction of [9] at the second route.
9. Section 8, "what the construction of [10] proposes": a partial read of [10] found a constraint with area alone, no volume term and no stationarity condition. Quote the equation or soften.
10. Say that 2 is the mean curvature of the unit sphere.
11. Notation clashes (Δ, δ, k, U); undefined sym, × and curl; (4) omits the 2π.
12. Lattice row: three digits inside a 13% bracket; the electric numbers are unlabelled.
13. "Casini, Galante and Myers showed when that works": they showed where it fails.

## Verified independently by the reader (all agree unless noted)

K traceless, transverse, harmonic, x·K·x = 0; E = 0 and B = −ε curl K from the 4D linearised Riemann tensor; ψ_B; both rates in 3(c); eq. (4) as the Lie derivative along ζ; ζ conformal Killing with κ = 1; eq. (5); the second-order constraint; A in terms of φ; ψ_E; Q_w; c_l and both tables; the electric table, a = 2/(l+2) and equality with (1) from its own Weyl computation (l = 1 to 6 exact, l ≤ 12 by closed form); the abstract's 8π/105, π/420, 4/3 and R⁶; eq. (6) for a general lapse; the entropy algebra of 3(c) given (a) and (b); the electric entropy without duality (l = 2, 3); the convention of (3) against [3]; the scalar comparison and the Wald term against [6]; ΔF as the area deficit at fixed volume. **Disagree as stated:** "no cross term" on the area side (E2); "l = 1 is flat space" (E1).

## References

All ten arXiv numbers, authors and titles match. Descriptions of [7], [8] and [10] were inexact (items 13, E4, 9).

## Not checked by the reader

Fact (b) and Haag duality for the ball; that P acts as "the unique finite-norm in-ball element"; the second route; the lattice numbers, the 5×10⁻¹⁰ search and the edge-mode statement; the electric area side for l > 12 and the electric direct entropy for l > 3; [8] and [10] beyond the first 100,000 characters of a summarised fetch; the claim that no prior computation of the entropy side exists.

## The reader's strongest objection

Why should −Δ(Area − 2Vol)/4G at the stationary surface equal the graviton modular energy at second order at all? The first law is a first-order identity for a conformal Killing vector, and eq. (6) shows the generator that moves the metric differs from the one that moves the Weyl tensor by exactly the non-Killing term, so a mismatch is expected, as for the minimal scalar. By [8] the second-order vacuum area deficit is shape-ambiguous; the electric agreement holds only at the stationary shape, and the magnetic number is already in [8]. In that language the paper's new content is S = πR⁶W/420G against an area deficit (πR⁶/315G)(W − E²/4). The reader will ask whether "electric holds" is an accident of one shape choice, and whether a graviton analogue of the Wald term would rescue the magnetic case.

---

## Lead's check of the read, and what was done

- **E4 confirmed independently.** The lead fetched arXiv:1710.07379 (ar5iv rendering, first 100,000 characters) and found eq. (75) as quoted, with the authors singling out γ = 1/(3d) from constant extrinsic curvature and noting that no γ makes the deficit proportional to W. Arithmetic at d = 4: prefactor 4πℓ⁶/315, bracket −(3/4)E² − B², so δA = −(πℓ⁶/105)(E² + (4/3)B²). The paper's l = 2 area numbers, including the 420 and the 4/3, are in print. **The claim that the area side was new at l = 2 is withdrawn.**
- **E1, E3, E5 and items 1 to 13:** accepted; version 2 of the paper applies them.
- **E2:** accepted as a correction of scope. Version 2 states the area side per sector and gives the condition for mixed data as an argument, not a computation. The lead did not recompute the boost derivative.
- **Closed form (item 5):** rechecked with the registered expansion code in `decider/scripts/d3_weyl_closed_form.py`.
- **What the read changes in the claim.** "The first law holds on electric data" is true at the stationary surface and not for a geodesic ball; the magnetic coefficient is the shape-independent one. The paper now leads with the entropy being the Bel–Robinson super-energy, S = πR⁶W/420G, which is the new statement, and presents the mismatch as (l+1)/(l+2) on magnetic data.
- **Errors of the lead that the read caught.** The registered control "l = 1 is flat Rindler space" (PREREG, C4) was mis-described: the control's zero is correct, its stated reason was not. The line "the cross term must vanish by parity" in the first report was wrong for the area side.
- **Grade: 405/420.** It recomputed more than twenty items from the text alone, found a prior-art match the programme had queued and never run, found a real scope error in the headline formula, and separated what it verified from what it could not see. Deduction: E2 rests on a reduction it did not verify directly, and two sources were read through a summarising fetch.
