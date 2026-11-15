# REPORT: lane CHAIN-REF (wave 73)

## Lead's check and grade

Sunday 15 November 2026. Adversarial referee lane (Fable), instructed to break CHAIN-ALL-L's proof that each graviton sector reduces to its scalar chain field for every l. The lane returned the text below; the lead saves it. It sent an interim finding to the lead before this report and read the bulletin, as the protocol asks; it discloses that the file read showed it the lead's header before the lane report.

**Grade: 388 of 420.** No break, and the report is worth more than a pass: it found a lemma the proof needs for its use on the ball and never states (locality), supplied its two-line proof, corrected two wordings that would have been wrong in print, did the odd-sector computation the proof waved at, replaced two "from memory" or indirect steps by short proofs, and confirmed the constant with a probe built from different ingredients than the lane's. Marked down because its blindness to the lead's opinion failed again by the same mechanism as the first referee's, and the explicit radial formulas at symbolic l remain unchecked by anyone but their author.

**Checked by the lead.**
- Re-ran `probe_chain_ref.py`: 9 s, exit 0, "DONE fails = 0"; all six ratios equal c_l = 1/6, 1/30, 1/90 at l = 2, 3, 4 on full Gram matrices, for harmonics with m = 0, 2 and 3, and the odd bilinear comes out as -c_l, as the referee says.
- The locality argument by hand: psi_E smeared with g is the curvature E smeared with g x(x)x, so the chain field lies in the curvature subspace of the same region; conversely E = T[q] is a differential operator in q with coefficients smooth in r^2, so integrating by parts moves it onto the test tensor without enlarging the support. It holds, for regions that are origin-centred balls.
- Not checked by the lead: the free-module degree count behind the span claim; the uniqueness argument of defect 5.

**Where the lead was wrong, and accepts it.** The header on CHAIN-ALL-L said the only remaining assumption was the definition of the net "plus second quantisation". It omitted the locality lemma, without which a statement about global one-particle spaces says nothing about the ball. With the lemma (proved here) the conclusion stands; the header should have named it. And the fault that let the referee see the lead's opinion first is the lead's: the header sits in the same file. From now on a referee is pointed at a copy without the header.

**What it changes for the paper.** Use the referee's "may use" sentence. Three corrections to carry: (1) say the chain field and its time derivative, smeared in the ball, generate the same one-particle local subspace as the sector's curvature; (2) c_l multiplies the inner product on data; the chain field's own two-point function is 1/c_l times the scalar's, so never call the chain field "a free scalar field with two-point function c_l times the scalar's"; (3) in the odd sector the roles of the two halves of the data are swapped (the metric datum feeds the chain field's time derivative), which is why the odd bilinear carries a minus sign while the symplectic form does not.

---

## Lane report (unchanged)

**CHAIN-REF verdict: no break found.** The file read showed me the lead's header before the lane report; my verdict was reached on the mathematics anyway.

## 1. Verdicts
- **(i) HOLDS.**
- **(ii) GAP (fixable, write-up only).** Injectivity is proved for every l below; surjectivity onto all smooth q beyond l = 12 rests on the lane's own derivation.
- **(iii) HOLDS**, with two wording defects and one missing lemma for its use on the ball.

## 2. Defects, most serious first
1. **Locality is never stated (bears on (h)).** The one-particle identity is global on R³, but the ball needs: curvature smeared in the origin-centred ball spans the same real subspace as (ψ, ψ̇) smeared there, sector by sector. Proof:
   - ψ_E(g) = E(g x⊗x), so the chain field lies in the curvature subspace.
   - Conversely E = T[q] is a second-order differential operator in q with coefficients smooth in r², so ∫f·E = ∫g ψ with supp g ⊂ supp f after integrating by parts. Sector-finite f are dense.
   - The unscaled map sends the dual vector of a functional to 1/c_l times the scalar's (a positive factor), and the unitary is √c_l times the chain map.
2. **(iii) title, "two-point function".** c_l multiplies the inner product on data. The chain field operator's two-point function is 1/c_l = Λ/4 times the scalar's. Since Λ = L²(L²−2) is a differential operator on the sphere, ψ_E stays local, but it is not a free scalar field across l.
3. **(iii).3, odd sector.** Done rather than accepted: K feeds ψ_B and h feeds ψ̇_B, so ∫K·h = −c_l∫ψ_B ψ̇_B. The antisymmetrised symplectic form is +c_l, so the map is sound; the report's "same computation" hides the swap. The script's L17 already tests −c_l.
4. **(ii).1 "from memory".** Repair: polynomial sector tensors form a free module over R[r²] with generators in degrees l−2, l, l, l+2 (even) and l−1, l+1 (odd). The lane's families have exactly these degrees, so the transition determinant has degree 0. It is non-zero because the families are pointwise independent for r > 0, l ≥ 2. Smooth fields follow by the standard smooth-equivariants theorem. The 6 radial functions split 4 + 2 as claimed.
5. **(ii).2 uniqueness.** Short proof replacing "determinant l+2": if T is TT with x·T·x = 0, then w = T·x is tangential and divergence-free of gradient type, so w = 0. Then T = a(r)Y_AB, whose divergence is non-zero unless a = 0. The odd case is the same with X_AB.
6. **(ii).5** says "algebraic"; it is second-order differential (still local).

## 3. Checked sound
- **(iii).1:** in polar coordinates the first divergence has zero radial component by tracelessness and only sphere derivatives by tangentiality; the second likewise. No radial derivative for any profile; k = 0 has measure zero; the Cauchy-data formulation covers both frequency signs.
- **(iii).2:** D^A D^B Y_AB = ½ΛY, re-derived.
- **(c):** both quadratic forms and the symplectic form scale by the same c_l in both parities; the map commutes with time evolution, so the complex structure is intertwined. ⟨EE⟩ = ½ω⁴P_TT times the scalar kernel is consistent.
- **(b):** ψ_B of even data and ψ_E of odd data vanish; Y- and X-types are orthogonal.
- **(ii).6:** kernels are right. Odd K = sym∇v is correctly not gauge.
- **f₄'s division by r²** is exact for every l (hand check).
- **(f):** Box E = 0 follows from Box h = 0 in TT gauge.
- **(g):** one m suffices by rotation covariance.

## 4. Probe
`chain-ref\scripts\probe_chain_ref.py` builds TT tensors by ∇⁴P_TT on seeds, with its own Fourier transform and none of the lane's radial formulas. It ran l = 3 at m = 0 and m = 3, plus l = 2 and l = 4, with profiles up to s⁴. All six ratios equal c_l on full Gram matrices (1/30, 1/6, 1/90), the odd bilinear is −c_l, exit 0. Log: `chain-ref\logs\01_probe.log`.

## 5. Not checked
- The explicit formulas of (ii).2 to (ii).4 at symbolic l.
- The record's units (c_l is convention-dependent: it scales with the graviton symplectic normalisation).
- Companion A8 and the other lane reports, which I did not open.

## 6. Bulletin and the lead's header
B21 changed my course: because the net is defined by curvature, I added the locality check. B17 means only per-sector structure is needed, so the l-dependent c_l is harmless. B22 is the claim under test.

I disagree with the header in two places:
- "What remains is (1) alone plus second quantisation" omits the locality lemma and the l-dependent rescaling.
- The span claim needs no markdown once defect 4's proof is written in.

## 7. Sentences
**May use:** "In each sector (l ≥ 2, either parity) the chain field and its time derivative, smeared in the ball, generate the same one-particle local subspace as the sector's linearised curvature, and the map from curvature data to chain-field data is symplectic and isometric up to the constant c_l = 4(l−2)!/(l+2)!, intertwining the time evolution."

**Must not use:** "The chain field is a free scalar field with two-point function c_l times the scalar's."
