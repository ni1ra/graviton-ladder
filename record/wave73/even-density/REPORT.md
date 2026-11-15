# REPORT: lane EVEN-DENSITY (wave 73)

## Lead's check and grade

Sunday 15 November 2026. Small proof lane (Fable): prove the functional-analytic facts that EVEN-PROOF, EVEN-REMAINDER and HINGE-PROOF's Lemma A assumed (self-adjointness, what the one-particle space is, density of smooth compactly supported data, extension of the bounds, isometry of the vacuum flow). The lane returned the text below; the lead saves it.

**Grade: 384 of 420.** Four questions, four verdicts, the one place where intuition could fail (the low-frequency weight on the momentum slot) tested directly with a prediction and a plant, and an honest split between what it proved and what it quoted. Marked down because every standard theorem is quoted from memory, and one inherited step (W maps into the domain of h) was not re-proved.

**Checked by the lead, by hand.**
- 1.1: with x = coth rho, dx/d rho = 1 - x^2, so chi = P_l(coth rho) satisfies -chi'' + l(l+1) chi/sinh^2 rho = 0 by Legendre's equation.
- 1.5: <f, h^2 f> = ||f'||^2 + lambda^2 ||f/sinh rho||^2 >= 0, and it vanishes only for f = 0.
- 2.2: |<pi, f>| <= ||pi sinh rho|| ||f/sinh rho|| <= ||pi sinh rho|| ||h f||/lambda, so every pi with pi sinh rho square-integrable lies in the domain of h^-1, which is what the momentum half of the norm needs.
- 3.6: N <= 1/sinh rho gives |<N q, f>| <= ||q|| ||h f||/lambda, so the odd perturbation is defined on all of H with no extension step.
- Not hand-checked: 1.7 (core lemma), 2.3 to 2.9, 3.1, and all quoted theorems.

**Re-run by the lead:** `C:\Python314\python.exe -B even_density.py <out.json>` from `scripts\`, about 3 s, exit 0, "ALL LINES HELD"; the plant (a zero-energy resonance) diverges linearly and is caught.

**What it changes.** "Density assumed" is closed, given the quoted standard theorems: smooth compactly supported data are dense in H, no extra condition on the momentum data is needed, the even and odd perturbations are bounded operators from H into H, and the vacuum flow is a strongly continuous isometry group. The standard theorems it quotes (Weyl's criterion, Reed-Simon X.7 and X.10, the spectral theorem) are being read first-hand by lane STD-THEOREMS.

---

## Lane report (unchanged)

**EVEN-DENSITY report (15 Nov 2026)**

The density that EVEN-PROOF, EVEN-REMAINDER and Lemma A assume is true, and no extra condition on smooth compactly supported π is needed. The check ran once, exit 0, every pre-fixed line held and the plant was caught.

## 1. Verdicts

1. **Q1:** proved for essential self-adjointness, h² ≥ 0 and no eigenvalue at 0, given two quoted theorems. "Purely absolutely continuous, spectrum [0, ∞)" is quoted from memory, not proved, and nothing downstream uses it.
2. **Q2:** proved. Every π in C_c^∞ lies in H, with ‖h⁻¹π‖ ≤ ‖π sinh ρ‖/λ. C_c^∞ ⊕ C_c^∞ is dense in H.
3. **Q3:** proved, given EVEN-PROOF steps 2 to 5 and Lemma A steps 1 to 2 as inherited. D_E and D map H into H itself.
4. **Q4:** proved, given the spectral theorem.

## 2. Proofs

Tags: [P] proved here, [S] standard theorem quoted from memory, [N] numeric, [I] inherited. Notation: a(x) = x coth(πx), b(x) = coth(πx)/x, E(f) = ‖f′‖² + λ²‖f/sinh ρ‖².

**Q1.**
- 1.1 [P] Zero-energy solutions are χ = P_l(coth ρ) and φ = χ∫₀^ρ χ⁻², with Wronskian φ′χ − φχ′ = 1. Near 0, χ ~ cρ⁻ˡ and φ ~ ρ^{l+1}/(2l+1); at infinity χ → 1 and φ ~ ρ − c_l. [N] c_l = 1.500000 at l = 2 and 3.380729 at l = 16.
- 1.2 [P] χ is not L² near 0 because 2l ≥ 2 > 1; φ is not L² near infinity.
- 1.3 [S, memory: Weyl alternative] One non-L² solution at one energy gives limit-point at that end. The Reed–Simon X.10 criterion at 0 is V ≥ 3/(4ρ²) near 0; here ρ²V → λ² = l(l+1) ≥ 2 > 3/4. At infinity V is bounded below (X.8).
- 1.4 [S, memory: Reed–Simon X.7] Limit-point at both ends gives essential self-adjointness on C_c^∞.
- 1.5 [P] ⟨f, h²f⟩ = E(f) ≥ 0 on the core, so h² ≥ 0. If h²f = 0 then E(f) = 0 on the form domain, so f = 0. Hence h⁻¹ is a densely defined positive operator, a(h) ≥ 1/π, b(h) is injective, and the norm is a norm.
- 1.6 [S, memory: Weidmann; Weyl] V ∈ L¹(1, ∞) gives purely absolutely continuous spectrum on (0, ∞); V → 0 gives essential spectrum [0, ∞).
- 1.7 [P] Core lemma: a core for A ≥ 0 is a core for A^θ, 0 < θ ≤ 1. Proof: ‖A^θf‖² ≤ ‖f‖² + ‖Af‖², and D(A) is a core for A^θ by spectral cut-offs. So C_c^∞ is a core for h and for h^{1/2}.

**Q2.**
- 2.1 [P] Lemma: π ∈ D(h⁻¹) if and only if |⟨π, f⟩| ≤ C‖hf‖ on a core of h, and then ‖h⁻¹π‖ is the best C. Proof: Riesz on the dense set h(core), then closure.
- 2.2 [P] For f in C_c^∞: |⟨π, f⟩| ≤ ‖π sinh ρ‖·‖f/sinh ρ‖ ≤ ‖π sinh ρ‖·‖hf‖/λ. So every π with π sinh ρ ∈ L², in particular every C_c^∞ π, is in D(h⁻¹).
- 2.3 [P] ⟨π, h⁻²π⟩ = ⟨π, Gπ⟩ with G = φ(ρ_<)χ(ρ_>). Proof: w = Gπ has E(w) = ⟨π, w⟩ finite and ⟨π, f⟩ = E(w, f), which gives "≤". Cut-offs wη_R have E(w − wη_R) = O(1/R), which gives equality.
- 2.4 [P] w tends to the constant ⟨φ, π⟩, so h⁻²π is not in L² unless ⟨φ, π⟩ = 0. That orthogonality is the condition for D(h⁻²), not for H.
- 2.5 [P] H = H_q ⊕ H_π. Since max(x, 1/π) ≤ a(x) ≤ x + 1/π, H_q = D(h^{1/2}) with an equivalent graph norm, already complete inside L². H_π is the completion of D(h⁻¹) in ‖b(h)^{1/2}π‖; in a spectral representation it is {f : ∫b(x)|f|²dμ < ∞}.
- 2.6 [P] H_π is neither inside L² nor contains it: the weight is ~1/x at high frequency (distributions) and ~1/(πx²) at low frequency.
- 2.7 [P] Equivalently, k ↦ z = (h^{1/2}q + i·h^{−1/2}π)/√2 is a real-orthogonal map of H onto D(h^{−1/2}) ⊂ L²(ℂ) with norm ‖coth^{1/2}(πh)z‖. Onto because a^{1/2}(h)H_q and the closure of b^{1/2}(h)D(h⁻¹) are both all of real L²; h is a real operator.
- 2.8 [P] Density in H_q: 1.7 with θ = 1/4.
- 2.9 [P] Density in H_π: B = b(h)^{1/2}h² satisfies ‖Bf‖ ≤ C(‖f‖ + ‖h²f‖). So B·C_c^∞ is dense in B·D(h²), which contains every spectral subspace E[a, b]L² with 0 < a < b, and their union is dense because E{0} = 0. Hence h²C_c^∞ ⊂ C_c^∞ is dense in H_π; even the subfamily with ⟨φ, π⟩ = 0 is dense.

**Q3.**
- 3.1 [P] W is Hilbert–Schmidt on L²(dρ): its kernel is at most C·e^{−ρ_>}.
- 3.2 [I] W maps L² into D(h) with ‖hW‖ ≤ 1/λ (EVEN-PROOF steps 2 to 4).
- 3.3 [P] EVEN-PROOF steps 5 to 7 need π ∈ D(h⁻¹). Step 2.2 supplies this; it was the unstated hypothesis.
- 3.4 [P] For π in C_c^∞, Wπ ∈ D(h) ⊂ H_q, so D_E maps the dense set into H.
- 3.5 [S, memory: bounded linear extension] A bounded map from a dense subspace into a complete space extends uniquely with the same norm. So D_E is bounded on H with ‖D_E‖ ≤ √2|g|/λ and range in H_q ⊕ 0.
- 3.6 [P] Odd case: for any q in H_q ⊂ L², |⟨Nq, f⟩| ≤ ‖q‖·‖hf‖/λ, so Nq ∈ L² ∩ D(h⁻¹) ⊂ H_π. D is defined on all of H with no extension step. Lemma A's form inequalities pass from C_c^∞ to D(h^{1/2}) by 1.7.

**Q4.**
- 4.1 [P] In the z picture the flow is z ↦ e^{−ish}z, checked by direct substitution.
- 4.2 [S, memory: spectral theorem, Stone] e^{−ish} is a strongly continuous unitary group on L² commuting with coth(πh).
- 4.3 [P] So it preserves D(h^{−1/2}) and its norm, and is strongly continuous there (apply 4.2 to coth^{1/2}z). Transport by 2.7.

## 3. Check

Finite differences (step 0.004, Dirichlet box) against the Green-function law I_L = I_∞ − ⟨φ, π⟩²χ(L)/φ(L). Here I_L is ⟨π, w⟩ on the box of size L.

| l | π | ⟨φ, π⟩ | I_8 | I_64 | I_512 | I_∞ | (I_2L − I_L) ratios, L = 16..256 | verdict |
|---|---|---|---|---|---|---|---|---|
| 2 | π1 = bump(1,3) | 0.30922 | 0.0972083 | 0.1103891 | 0.1117317 | 0.1119189 | 2.155 → 2.009 | converges like 1/L |
| 2 | π2, orthogonal | 0 | 0.0828171 | 0.0828171 | 0.0828171 | 0.0828169 | differences ≤ 3e-16 | constant |
| 2 | π3 = bump(0.5,4.5)·cos 4ρ | −0.30979 | 0.0728711 | 0.0860998 | 0.0874473 | 0.0876342 | 2.155 → 2.009 | converges |
| 16 | π1 | 7.9505e-3 | 6.706836e-3 | 6.719477e-3 | 6.720396e-3 | 6.720513e-3 | 2.402 → 2.020 | converges |
| 16 | π2 | 0 | 6.676479e-3 | 6.676479e-3 | 6.676479e-3 | 6.676473e-3 | differences ≤ 2e-17 | constant |
| 16 | π3 | −0.12266 | 0.0457122 | 0.0487210 | 0.0489396 | 0.0489686 | 2.402 → 2.020 | converges |

- **Agreement:** finite differences match the law to 1.2e-5 relative or better; successive differences match the predicted ones to five digits.
- **Bounds of 2.2:** I_∞ ≤ ‖π sinh ρ‖²/λ² holds in all six cases, tightest at l = 16, π1 (6.72e-3 against 8.11e-3).
- **Plateau (2.4):** on the box L = 512, w(256) = 0.15507 at l = 2 against ⟨φ, π1⟩/2-scale prediction 0.15507; for π2 it is below 2e-7 of max|w|.
- **Plant** (no potential, Neumann at 0, i.e. a zero-energy resonance): I_L = 1.14, 12.2, 100.5 at L = 8, 64, 512, growing as L(∫π)² to five digits. The same detector returns DIVERGES. Caught.
- **Failed lines:** none.
- **Departures from the first prereg text, both made before the first run:** the ratio window was widened from [1.8, 2.3] to [1.8, 2.6] (logged with the arithmetic), and a box L = 1024 was added only to form the last difference.

## 4. Re-run

From `record\wave73\even-density\scripts`:

`C:\Python314\python.exe -B even_density.py ..\out\even_density.json`

It takes about 3 s and exits 0. Logs are `logs\01_prereg.md` and `logs\02_run1.log`. Nothing was committed and no REPORT.md was written.

## 5. Weakest step

Step 3.2: I inherited that W maps L² into D(h) with ‖hW‖ ≤ 1/λ from EVEN-PROOF and proved only that W is Hilbert–Schmidt; every [S] item is quoted from memory, not read.

## 6. Corrections to other lanes

- **EVEN-PROOF section 6, EVEN-REMAINDER section 5 and list item "density":** replace "assumed" by a citation of 2.8, 2.9 and 3.5. EVEN-PROOF steps 5 to 7 should cite 2.2 for ⟨π, h⁻²π⟩ < ∞.
- **HINGE-PROOF Lemma A:** no extension is needed (3.6).
- **EVEN-REMAINDER step 1:** correct as written. Add that z ranges over D(h^{−1/2}) ⊂ L², while the π slot alone may be a distribution.
- **No mathematics changes; no extra condition on π.** The zero-energy solution's only job is to decide whether h⁻²π is square-integrable, and nobody asked it.
- **Family vectors:** EVEN-REMAINDER step 7 gives EVEN-PASSAGE's vector as q = 0, π = a_l = h²e^{−εh²}u. It is not compactly supported, but h⁻¹a_l = h·e^{−εh²}u ∈ L², so it is in H directly. I did not read EVEN-NONZERO or EVEN-PASSAGE (outside my read fence), so their vectors are otherwise unverified; any π with π sinh ρ ∈ L², or π ∈ Ran(h), qualifies.

The plant diverged exactly on schedule, which makes it the most punctual thing in the folder.
