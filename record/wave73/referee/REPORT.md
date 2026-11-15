# REPORT: lane REFEREE (wave 73)

## Lead's check and grade

Sunday 15 November 2026. Adversarial referee lane (Fable), instructed to break the assembled argument in `assembly\REPORT.md` (the two statements, the chain B.1 to B.9 and section C). The lane returned the text below; the lead saves it. The lane wrote no files and ran no code.

**Grade: 384 of 420.** Verdict HOLDS, and it earned it: every step of section C is addressed with the specific way it could fail, a toy model is used to show what the single-time test can and cannot see, and it found one statement that is false as written and replaced a from-memory hypothesis by a four-line computation. Marked down because its blindness to the lead's opinion failed (the file read returned the lead's header with the lane text, which it discloses), and because its own Tomita computation was by hand with no check.

**Checked by the lead, by hand.**
- The Tomita computation in the first defect: with k(z) = (c z, s z-bar), k(z) + i k(z') = (u, v) has u = c(z + i z'), v = s conj(z - i z'), so k(z) - i k(z') = ((c/s) v-bar, (s/c) u-bar). The adjoint of that antilinear map is (a, b) -> ((s/c) b-bar, (c/s) a-bar), so Delta_H = (s/c)^2 (+) (c/s)^2 = e^(-2 pi h) (+) e^(2 pi h), and Delta_H^(it) k(z) = k(e^(-2 pi i t h) z). It holds.
- C.6: an expression linear in a real s with values in pi Z is zero; the lead's "continuity at s = 0" was unnecessary. Agreed.
- The scope point: Shale and Connes are statements inside the given representation, so the conclusion is about weights and states normal on M. Agreed.

**Where the lead was wrong, and accepts it.**
- The lead's header on ASSEMBLY let the remark "the sign convention for t is immaterial" stand. It is false as written; only the relative orientation of the two flows is immaterial. The orientation must be stated and proved, as the referee does.
- "Implemented by no unitary in the region's algebra" was weaker than what is shown: the map is implemented by no unitary on Fock space at all.
- "KMS implies modular" should not remain a from-memory hypothesis: the one-particle Tomita computation plus the second-quantisation theorem (the modular operator of the algebra of a standard real subspace is the second quantisation of the subspace's) replaces it. That theorem is now the one citation to verify; the literature lane PHYS-ID has been asked for it.

**What it changes for the paper.** Use the referee's "may use" sentence. Say "normal". Define the ambient space K first. State the orientation of the vacuum modular flow and prove it by the Tomita computation. Do not say "no KMS or thermal state exists", do not say "for every t" until lane ALL-T reports, and do not present the identification of M with the ball's graviton algebra as proved. The referee's toy case is worth a sentence: a single-time test can be blind (frequencies differing by integers), so the divergence rests on the proved floors, not on the abstract argument.

---

## Lane report (unchanged)

**1. Verdict: HOLDS.** There is no mathematical break in section C or B.1 to B.9; there are three wording defects and one removable dependency. Disclosure: the file read returned the lead's header together with the lane text, so I saw it before forming my verdict, not after.

**2. Defects, most serious first**

- **B.4, "the sign convention for t is immaterial".** This is false as written. Only the relative orientation of U^g and U⁰ is immaterial; if σ^ω_t were W(U⁰(−2πt)·), the cocycle would be U^g(π)U⁰(π), not S. Repair: state σ^ω_t(W(k)) = W(U⁰(2πt)k) and prove it. In beta-def's model K = L ⊕ L̄, k(z) = (cz, s·z̄):
  - The Tomita operator of H is k(z) + ik(z′) ↦ k(z) − ik(z′), i.e. (u,v) ↦ ((c/s)v̄, (s/c)ū), with c/s = e^{πh}.
  - Hence Δ_H = e^{−2πh} ⊕ e^{2πh}, and Δ_H^{it}k(z) = k(e^{−2πith}z) = k(U⁰(2πt)z), using even-density 4.1.
  - With Δ_M = Γ(Δ_H) (Eckmann–Osterwalder; Leyland–Roberts–Testard) the orientation is correct.
  - This also removes "KMS implies modular" and Kay–Wald from the theorem's hypotheses, so the lane's "weakest step" shrinks to one citation.
- **Scope of the conclusion.** The result excludes weights and states normal on M, i.e. in the thermal folium. It says nothing about states of the Weyl C*-algebra in a disjoint representation, where a KMS state for the perturbed flow may exist. The theorem must say "normal".
- **Setting, "K = closure(H + H′)".** This is circular, since H′ needs an ambient space. Define K = ⊕_{l,m}(L ⊕ L̄) first. The same definition gives the complex orthogonality of blocks that B.9 uses.
- **C.6, cosmetic.** No continuity is needed: s ↦ s·(Im⟨Ŝ̄k,Ŝ̄j⟩ − Im⟨k,j⟩) is linear with values in πℤ, hence zero.

**3. Checked and sound**

- **C.1–C.2:** An escaping subsequence is handled: ‖j_n‖ → ∞ gives ‖W(j_n)Ω − Ω‖² → 2, not 0.
- **C.3:** Correct, because M is defined on Fock space over K; no GNS identification enters.
- **C.4–C.5:** Conjugation by a fixed unitary preserves strong limits; the real-linear boundedness argument and the extension are fine.
- **C.7:** ŜD = D, so the U* argument applies.
- **C.8:** Shale's "implementable" is exactly UW(k)U* = W(Tk); ½(T + iTi) is the antilinear part for the complex structure i and the form Im⟨,⟩.
- **C.9:** D_i = D ∩ iD is a complex subspace; β is well defined there and is the restriction of β̄.
- **B.4:** The composition order is right: σ^ψ_{1/2}∘σ^ω_{−1/2} on W(k) gives W(U^g(π)U⁰(−π)k).
- **B.5:** W(h + h′) = W(h)W(h′) since Im⟨H,H′⟩ = 0, and likewise for Sh, so u implements S ⊕ 1 with no stray phase.
- **B.3:** H ∩ H′ = 0 and H ∩ iH = 0 hold (c, s > 0), so Ŝ is well defined and Ω is separating.
- **B.8:** I re-derived i·k(ψ) = k(i coth ψ) − k′(i csch ψ) from c² − s² = 1; both pieces have norm at most ‖coth^{3/2}ψ‖.
- **B.9:** For antilinear A, Σ‖Af_j‖² = Tr(A†A) is basis-independent, so an orthonormal family gives a lower bound. The tower sum is right.
- **(j):** No automorphism-group hypothesis is needed; only σ^ψ_{1/2} on W(H) is used.

**Toy cases.** No counterexample was found.
- Type I (H complex, H′ = H^⊥): the argument reduces to Shale on H, as it must.
- Countable oscillators with ψ the thermal state of rescaled frequencies h_n(1+ε_n): if every ε_n h_n is an integer then S = 1 at t = ½, although the states are inequivalent unless the ε_n are small. So the t = ½ test is necessary only and can be blind; the divergence rests entirely on the numerical floors c, not on section C.

**4. Not checked**

- The floors c, l₁ and the remainder bounds.
- The theorem numbers in Derezinski and Hiai.
- The Γ(Δ_H) citation.
- The Tomita computation above is by hand; I ran no code and wrote no files.

**5. Disagreements with the lead**

- "Implemented by no unitary in the region's algebra" understates the result: Ŝ is implemented by no unitary on Fock space at all.
- KMS-implies-modular should not stay as a from-memory hypothesis; the one-particle computation above replaces it.
- The lead accepts the sign remark in B.4; I do not.
- The lead's C.6 gloss (continuity at s = 0) is superfluous.
- The grade is fair. A lane refereeing its own section C is like a defendant bringing his own gavel, but here the gavel happens to be sound.

**6. Sentences**

- **May use:** "For every real g ≠ 0, the Bogoliubov map S ⊕ 1 with S = U^g(π)U⁰(−π) is implemented by no unitary on the Fock space over K; consequently no faithful normal semifinite weight on M, in particular no faithful normal state, has modular automorphism at t = ½ equal to W(k) ↦ W(U^g(π)k), under the standard theorems of Shale and Connes and the stated numerical floors."
- **Must not use:** "The perturbed flow admits no KMS or thermal state", "is not a modular flow for any t ≠ 0", or any wording that presents the identification of M with the ball's graviton algebra as proved.
