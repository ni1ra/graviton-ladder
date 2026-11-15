# REPORT: lane BETA-DEF (wave 73)

## Lead's check and grade

Sunday 15 November 2026. Small proof lane (Fable): derive from first principles the antilinear (Bogoliubov) part of the perturbed flow and its first-order term, which HINGE-PROOF, EVEN-NONZERO and EVEN-PASSAGE had used in a form reconstructed from code. The lane returned the text below; the lead saves it.

**Grade: 386 of 420.** The reconstruction is confirmed as exact, the definition side of the check is built without the formula's own ingredients, all 14 pre-fixed lines passed, and the lane corrected the lead's brief (the map must be the cocycle, not the bare exponential). Marked down because the equivalence with implementability inside the region's algebra (step 4) rests on two textbook theorems quoted from memory and on a boundedness statement inherited from HINGE, and a finite-dimensional check cannot test those.

**Checked by the lead, by hand.**
- c^2 + s^2 = coth(pi h) and c^2 - s^2 = 1, so k(z) = (c z, s z-bar) reproduces the thermal norm and the symplectic form (step 2).
- 1/c^2 + 1/s^2 = (1 - e^(-2 pi h)) + (e^(2 pi h) - 1) = 2 sinh(2 pi h); with P + Q = T z, minimising ||c P||^2 + ||s Q||^2 gives the bound of step 7.
- Step 6: Duhamel's formula for the cocycle, and real-linearity of the map from S - 1 to beta.
- Not hand-checked: step 3 (density and trivial intersection), the collection of antilinear terms in step 5 (covered by the check to 1e-15), steps 8 and 9, and the standard facts of steps 1 and 4.

**Re-run by the lead:** `C:\Python314\python.exe beta_def.py <out.json>` from `scripts\`, about 1 s, exit 0, all 14 lines PASS; numbers equal the lane's (exact-formula difference 9.0e-15 to 1.6e-11, decade ratios 10.93, 10.13, 10.01, plants A, B, C caught).

**Correction to the lead's own brief.** The brief wrote S = exp(tau(X_0 + D)). The object is the cocycle S_t = exp(tau(X_0 + g D)) exp(-tau X_0); without the second factor the zeroth order does not cancel. The lane caught this.

**What it changes.** The "weakest step" of EVEN-NONZERO and EVEN-PASSAGE (the form of the first-order term) is closed at the level of a derivation plus two standard theorems. Wording to fix when the paper is written: HINGE-PROOF's beta_1-hat carries a factor lambda that the even lanes' does not; pick one normalisation.

**Where the even divergence stands after this lane.** Still **not a theorem**. Remaining: (1) the remainder for EVEN-PASSAGE's family, sup over l of C(psi_l)/||e_l|| finite, with a proved lower bound on ||e_l||; (2) HINGE-PROOF's Duhamel-Gronwall corollary, used for ||R_2||, not lead-checked; (3) the two standard theorems of step 4 read first-hand, and HINGE's boundedness of the extension; (4) EVEN-PROOF's dictionary making W_l the even charge, validated only through the spectrum.

---

## Lane report (unchanged)

**BETA-DEF report (2026-11-15)**

EVEN-NONZERO's form is **equal** to the first-principles formula, with no factor, conjugation or extension discrepancy. It is in fact exact to all orders: put A, B of S − 1 in place of those of T. Steps 1 and 6 survive. All 14 pre-fixed lines passed.

## 1. Definition and formula

- **One-particle structure.** K = L ⊕ L̄ with i acting as i on both slots, κ(k) = k(z) = (cz, s·z̄). The commutant copy is H′ = {k′(w) = (sw, c·w̄)}.
- **Extension.** Ŝ = S ⊕ 1 on H + H′, where S is the cocycle S_t = e^{τ(X₀+gD)}e^{−τX₀}.
- **Antilinear part.** β = ½(Ŝ + iŜi). Writing (S − 1)z = Az + Bz̄:

  β(u, v) = (cBc·ū − cAs·v̄, s·Ā·c·ū − s·B̄·s·v̄), so β k(z) = (cP, s·Q̄) with P = B(c²z̄) − A(s²z), Q = A(c²z) − B(s²z̄).

- **Squared Hilbert–Schmidt norm.** ‖cBc‖² + ‖cAs‖² + ‖sAc‖² + ‖sBs‖², which is the four terms of `hs_beta`.
- **First order.** β̂₁ is the same expression with A, B taken from T = ∫₀^τ e^{sX₀}De^{−sX₀}ds.

## 2. Derivation

1. [S, from memory: Kay–Wald 1991 Prop. 3.1; Araki–Woods] A quasi-free state with real inner product μ dominating ω has a one-particle structure (K, κ): κ real-linear, Re⟨κk, κk′⟩ = μ, Im = ω/2, κH + iκH dense. It is unique up to a unitary, and the GNS space is Fock space over K.
2. [P] k(z) = (cz, s·z̄) is such a structure: Im⟨k(z), k(z′)⟩ = Im⟨z, (c² − s²)z′⟩ = ω/2, Re = ⟨z, coth(πh)z′⟩, and k(z) + ik(z′) = (cx, s·ȳ) with x, y free.
3. [P] Im⟨k(z), k′(w)⟩ = 0. Any (u, v) with z = cu − s·v̄ and w = c·v̄ − su in the form domain equals k(z) + k′(w), so H + H′ is dense and H ∩ H′ = 0.
4. [S, from memory: Araki duality R(H)′ = R(H′); Shale 1962] A unitary u in M with uV(k)u* = V(Sk) commutes with V(k′), so it implements S ⊕ 1. Conversely, any Fock implementer of S ⊕ 1 commutes with R(H′), hence lies in M. So "S is implemented inside M" is equivalent to "the bounded extension of S ⊕ 1 has a Hilbert–Schmidt antilinear part". Boundedness of the extension is [I] hinge Proposition 1, step 3.
5. [P] Ŝ − 1 sends (u, v) to k((S − 1)z). Substituting z = cu − s·v̄ and collecting the terms antilinear in (u, v) gives the formula in section 1. Setting u = cz, v = s·z̄ gives (cP, s·Q̄).
6. [P] S_t − 1 = ∫₀^τ e^{s(X₀+gD)} gD e^{−sX₀} ds = gT + R₂. The map from S − 1 to β is real-linear, so β = gβ̂₁ + β[R₂].
7. [P] Step 1 of EVEN-NONZERO: P + Q = Tz. Since 1/c² + 1/s² = 2 sinh 2πh, ‖β̂₁e‖ ≥ ‖(2 sinh 2πh)^{−1/2}Tz‖.
8. [P] Step 6: ⟨η, β̂₁e⟩ = ⟨iz, P⟩ + conj⟨iz, Q⟩, whose real part is Im⟨z, Tz⟩.
9. [P] HINGE-PROOF line 63: i·k(ψ) = k(i coth(πh)ψ) + k′(−i csch(πh)ψ).

## 3. Check

The definition side uses no c, s or z: K = C¹² with Gram matrix M + (i/2)Ω, H′ from an SVD null space. n = 6.

| g | exact formula, relative difference | first-order residual r(g) | EVEN-NONZERO form on one vector |
|---|---|---|---|
| 1e-1 | 9.0e-15 | 1.049 | 1.114 |
| 1e-2 | 1.2e-13 | 9.60e-2 | 1.06e-1 |
| 1e-3 | 1.4e-12 | 9.48e-3 | 1.05e-2 |
| 1e-4 | 1.6e-11 | 9.46e-4 | 1.05e-3 |

- **Scaling:** decade ratios are 10.93, 10.13, 10.01, so the absolute difference scales as g².
- **Agreement:** β̂₁ by definition against (cP, s·Q̄): 1.7e-15. Hilbert–Schmidt norm against the four terms: at most 2e-14. Step 1: 3e-16. Step 6: 1e-15. Line 63: 3e-15.
- **Plant A** (c and s swapped): mismatch 1.04. Caught.
- **Plant B** (Q in place of Q̄): mismatch 0.26. Caught.
- **Plant C** (S ⊕ S′ in place of S ⊕ 1): squared norm 0.02886 against 0.02388, so the extension matters.
- **Without the cocycle factor:** e^{τX₀} ⊕ 1 has squared norm 0.717 at g = 0.
- **Failed lines:** none.

## 4. Re-run

From `record\wave73\beta-def\scripts`:

`C:\Python314\python.exe beta_def.py ..\out\beta_def.json`

It takes about 1 s and exits 0. Logs are in `logs\`, numbered 01 to 04; runs 1 and 2 are identical after a dead-code tidy.

## 5. Weakest step

Step 4: Araki duality and Shale's necessity are quoted from memory, and boundedness of the extension of S ⊕ 1 is inherited from hinge; the check is finite-dimensional, so it cannot test either.

## 6. Verdict on point (3)

HINGE-PROOF's "Why this closes the gap" paragraph is correct, with these precisions:

- **Sense of "first order".** For e = k(ψ): βe = ½[k((S−1)ψ) + i·k((S−1)(i coth ψ))]. Hence ‖βe − gβ̂₁e‖ ≤ ½‖R₂‖_H(‖coth^{1/2}ψ‖ + ‖coth^{3/2}ψ‖), with ‖R₂‖_H ≤ 4π²g²t²e^{δt}/λ² from the corollary [I, not lead-checked].
- **Lower bound.** ‖β‖_HS ≥ ‖β‖_op ≥ ‖βe‖/‖e‖ ≥ g‖β̂₁e‖/‖e‖ − C(ψ)g²/(λ²‖e‖).
- **Uniformity [U].** This needs sup_l C(ψ_l)/‖e_l‖ finite, which is still unchecked for EVEN-PASSAGE's family.
- **Extension.** S ⊕ 1 is legitimate, and is the right one for the region's algebra, by step 4.

Wording changes needed; no mathematics changes:

- **Normalisation.** HINGE-PROOF writes "(g/λ)·β̂₁" with liminf‖β̂₁e‖ > 0, so its β̂₁ carries a factor λ. The even lanes use β̂₁ unscaled, with λ‖β̂₁e‖ bounded below. The paper should pick one.
- **Which S.** S must be the cocycle e^{τ(X₀+gD)}e^{−τX₀}, not e^{τ(X₀+D)} as in this lane's brief. Without the factor the zeroth order does not cancel.
- **Even lanes.** EVEN-NONZERO section 5 and EVEN-PASSAGE step 8 can cite steps 2 to 5 above as the definition; their "weakest step" is closed at the level of [P] plus the two [S] facts.
