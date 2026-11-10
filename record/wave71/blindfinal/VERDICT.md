# Blind referee verdict on the ladder theorem's new proofs

Lane w71-blind-final (Fable), Tuesday 10 November 2026, Oslo. Blind: read only `wave69/final/THEOREM.md` Sections 1–8, `wave69/duality/REPORT.md` §§1–2, `wave70/modflow/REPORT.md` §§1–2, the scripts of `wave69/final/scripts/` with their outputs, and open texts in `sources/`. Nothing committed.

**Overall verdict: no step is broken.** Every step I was asked to attack holds. Nine places need a sentence added or corrected; none changes a statement, a number or a sign of the theorem. The value S = ((l+1)/(l+2))·A/8G survives an independent check of its integer (below, M3).

**What was run.**
- `wave69/final/scripts/run_all.sh`, rerun here: exit 0, 77 s. 304 + 81 + 72 + 4 = **461 checks**, 56 + 32 + 18 + 3 = **109 plants**, all behaving; ledger audit 62 nodes pass; the eight record scripts rerun. Log: `run_all.log`. Twelve of the thirteen output files are byte-identical to the stored ones; `upstream/modflow_s1.json` differs only in its wall-clock field `"seconds"` (17.5 → 20.3), which is why `SHA256SUMS` changed.
- My own script `scripts/p_blind.py`: **62 checks, 39 plants**, all behaving (44 s). Output `out/p_blind.json`. At least one plant per proof; they are listed with each row.

---

## 1. Verdict by step

| # | Step | Verdict | The fix, or what was checked |
|---|---|---|---|
| 1a | SQ (a): H′ = J_H H, H′ standard | **holds** | Checked line by line, including the antilinear adjoint S_H* = J_HΔ_H^{−1/2} and the two orthogonality claims. |
| 1b | SQ (b): Ω cyclic and separating | **holds** | The identity-theorem step needs Exp to be analytic (not antianalytic) in its argument; it is. |
| 1c | SQ (c): Δ_{R(H)}^{it} = Γ(Δ_H^{it}) | **holds with a fix** | Fix F1 below: the KMS condition is verified on the span 𝔄 of Weyl operators, but the theorem cited (T7) needs it on all of R(H) or on entire elements. Four lines close it. The sign β = −1, the strip −1 ≤ Im z ≤ 0 and the boundary value f(t − i) = ⟨Δ^{it}k, h⟩ are right. |
| 1d | SQ (d): J_{R(H)} = Γ(−J_H), R(H)′ = R(H′) | **holds with a fix** | Fix F2: "by the functional calculus" should be the analytic-continuation criterion (now on file). Γ(−J_H), not Γ(+J_H), is correct in the convention V(u)Ω ∝ Exp(u). |
| 1e | SQ (e): additivity | **holds** | — |
| 2a | M-a: the state k, and independence of the representative | **holds with a fix** | Fix F3: (K6) is used for every l but supported only by a script at l ≤ 5. It follows from (K4) in one line. |
| 2b | Lemma E (flow acts geometrically where smooth) | **holds** | Checked: Lemma 2(i) is applied on U = O ∪ (a neighbourhood of supp G), where the flow is smooth for |s| ≤ s₀; Lemma 1 holds for all test fields; density closes it. |
| 2c | Lemma G (generator) | **holds** | Checked the bound: the terms of δ_ζ†G_n are a(x)t^mδ_n^{(q)}(t), m ≤ 2, q ≤ 1; their t-transforms are ≤ C(1 + ω) uniformly in n; the wave function gains ω^{3/2} and a ∈ C_c^∞, so dominated convergence applies. Closedness of log Δ and Stone are used correctly. |
| 2d | Equation (2): (E, B)[i log Δ k] = 2π((l+1)K, 0) on B | **holds** | Independent check M3: weight 3 is forced, and gives l + 1. |
| 2e | M-c: the cut vector c·2NK, c = (l+1)/(l+2) | **holds** | Checked: finite Ḣ^{1/2} norm; R_lin = 0 as a distribution (no surface term because N = 0 and K_nn = 0 on S); Ric₃ gauge invariant; (K4). Remark R1: the number rests on factoriality, which excludes a competing vector only by a logarithm. |
| 2f | H(O) factorial (no boost-invariant vector) | **holds** | Checked line by line. |
| 2g | M-d: the value | **holds** | Normalisation Im⟨u, u′⟩ = (1/32πG)∫(h·p′ − h′·p) re-derived from the Einstein–Hilbert action; I1 read in Ciolli–Longo–Ruzzi (eq. (10), Cor. 2.3, Thm 4.5): S_k = −Im(k, P_H i log Δ k), inner product linear in the second slot, as used. |
| 3 | (C), step 2: the support projection | **holds with a fix** | Fix F4: "S_M = S_{sMs}" is asserted without a locator or proof; give either. The rest (Φ(1 − s) = 0, Φ(m) = Φ(sms), same s for both states, the two applications of I3) is right. |
| 4a | B-d: zero modes of Q, and the value at every first-order stationary configuration | **holds with a fix** | Fixes F5, F6, F7 (wording). The form Q, its seven-dimensional kernel and the conclusion are right. |
| 4b | B-e: the Poincaré modes are not driven | **holds with a fix** | Fix F8: 𝒬 must be defined for first-order data that do not vanish on the slice. The parity argument is right. Single l is essential (plant). |
| 5a | modflow Lemma 4, item 5 (cores) | **holds** | Checked; the core claim for e^{πA} is true and a two-line proof is given below (N1). |
| 5b | modflow Lemma 5, continuity at the strip edges | **holds** | Checked: the bound C_N(F)(1 + |q|)^{−N} is uniform in θ ∈ [0, π] including both edges, so χ^{(z)} is L²-continuous on the closed strip by dominated convergence; at θ = π the vector is e_{F′} with F′ real and supported in the opposite wedge. |
| 5c | duality Lemma A | **holds** | Checked: dilation scales the Ḣ^{−3/2} norm by λ³; mollification is a contraction converging strongly; T̂(0) = 0 for a compactly supported TT field. |
| 5d | duality Lemma B | **holds** | Checked each of the five steps of the chain with indices, and (P). |
| — | Reproducibility | **holds with a fix** | Fix F9: `SHA256SUMS` is not stable across reruns. |

---

## 2. The fixes

**F1 (SQ (c): from the Weyl span to the whole algebra).** Sorce's Theorem 3.2 (the proof of T7 now on file) and Bratteli–Robinson's F-form both ask for the strip condition for every pair A, B in the von Neumann algebra; the definition through a dense subalgebra asks for entire elements, and V(h) is not entire for σ. Add:

> *Lemma (density).* Let σ_t = Ad U(t) with U(t)Ω = Ω, and let 𝔄 ⊂ R be a σ-weakly dense *-subalgebra such that for A, B ∈ 𝔄 there is F_{A,B}, bounded and continuous on −1 ≤ Im z ≤ 0, analytic inside, with F(t) = ω(Aσ_t(B)) and F(t − i) = ω(σ_t(B)A). Then the same holds for all A, B ∈ R.
> *Proof.* Both boundary values are bounded by ‖A‖‖B‖, so |F_{A,B}| ≤ ‖A‖‖B‖ on the strip (three lines). By Kaplansky's density theorem take nets A_α → A, B_α → B in 𝔄, strongly-* and bounded in norm. Then ω(A_ασ_t(B_α)) = ⟨A_α*Ω, U(t)B_αΩ⟩ and ω(σ_t(B_α)A_α) = ⟨U(t)B_α*Ω, A_αΩ⟩ converge uniformly in t. By the three-lines theorem applied to differences, F_{A_α,B_α} converges uniformly on the closed strip; the limit is F_{A,B}. ∎

Also say in (c) that Δ_H^{it}H = H and J_HΔ_H^{it} = Δ_H^{it}J_H are one-particle facts (from S_H² = 1), not consequences of T6.

**F2 (SQ (d)).** Replace "by the functional calculus" with: z ↦ Exp(Δ_H^{iz}h) is continuous on −½ ≤ Im z ≤ 0 and analytic inside (h ∈ Dom Δ_H^{1/2}, and Exp is entire), and equals Δ_{R(H)}^{it}Exp(h) on the real line; so Exp(h) ∈ Dom Δ_{R(H)}^{1/2} with image Exp(Δ_H^{1/2}h) (Sorce Thm 2.13, pdf p15).

**F3 (K6 for every l).** For every symmetric f, Cf = 2 curl(Ric₃(f) − ¼δ R_lin(f)), because inc f = 2Ric₃(f) − δR_lin(f) (checked exactly on a generic polynomial). With r²K = K − 2NK, Ric₃(K) = −½ΔK = 0, (K4) and tr K = 0: Ric₃(r²K) = −(l+2)K and R_lin = 0, so C(r²K) = −2(l+2) curl K for every l. Put this line in Section 1 next to (K6). (Checked directly at l = 6 as well.)

**F4 (C, step 2).** Either cite or prove "S_M(ω‖φ) = S_{sMs}(ω‖φ) when both supports equal s". A proof: in a standard form (M, ℋ, J, P) put q = s·JsJ; then (sMs, qℋ, qJq, qP) is a standard form of sMs, the vector representatives of both states lie in qℋ, and the relative modular operator Δ_{φ,ω} has support s(φ)·Js(ω)J = q and restricts to the relative modular operator of the reduced algebra. Note also that Reible's Theorem IV.3.7 (read here, pdf p68) asks for a cyclic vector only for the first state ω_i, and for ω₁(1) = ω₂(1); both hold. So the support step is needed for ω_k∘Φ only, as written.

**F5 (H5's last sentence is false as written).** "Only the Hamiltonian constraint on t = 0 at O(ε²) is used" holds for B-b and B-c. B-d and B-e apply B-c to P_a*g_ε for time translations and boosts, that is, they use the Hamiltonian constraint of g_ε on boosted and shifted slices. That is G_{μν}n^μn^ν = 0 for a neighbourhood of normals, which contains the momentum constraint. Replace the sentence by: "B-b and B-c use only the Hamiltonian constraint on t = 0 at O(ε²); B-d and B-e use it on every Poincaré image of that slice."

**F6 (B-d, the Poincaré part).** The sentence "P_a*g_ε is a family whose first-order datum differs from 2tK by O(ε), which changes F at O(ε³)" should read: with a = εa₁, P_a*g_ε = η + 2εtK + ε²(h⁽²⁾ + a₁£_X(2tK)) + O(ε³). The first-order part is unchanged; the change is a second-order term, and B-b, B-c hold whatever h⁽²⁾ is. This is exact and simpler than the argument printed.

**F7 (B-d, "nondegenerate" in infinite dimensions).** State the sense: on the complement of the seven modes the eigenvalues of Q are ½(j−1)(j+2) ≥ 2 on σ, −½j(j−1) ≤ −1 on τ and the Dirichlet form on s, all growing like j²; so DQ is an isomorphism H¹ → H⁻¹ there, and "gradient o(ε) implies d^⊥ = o(ε)" is meant in H¹. Also: H4 says "the stationary configuration nearest the coordinate ball", but what is proved is stationarity through O(ε²) and the value at every configuration stationary to first order. Say in H4 that the configuration is determined up to the Poincaré family and to the order shown; the value does not depend on the choice.

**F8 (B-e, the functional 𝒬).** For time translations and boosts, P_a*h⁽¹⁾ does not vanish on t = 0 (already at first order in a its spatial part is 2aK for a time translation), so "B-b and B-c applied to P_a*g_ε" does not define 𝒬 there. Add: for a general first-order solution the ε-coefficient of F at the coordinate ball is −½∫_B N R_lin[q₁] = 0 by the linearised constraint, and the ε²-coefficient is (δ²A − 2δ²V)[q₁, q₁] − ½∫_B N R_lin[q₂], where the second-order Hamiltonian constraint gives R_lin[q₂] as a quadratic expression in the first-order data. So 𝒬 is a quadratic functional of h⁽¹⁾ alone, invariant under the symmetries of (η, B), and the parity argument applies. For spatial translations the printed argument is complete and I checked it in closed form: the ε²-coefficient at the ball moved by a is −½∫_B N|K(x + a)|², with gradient −½∫_B x_a|K|² = 0 for a single l.

**F9 (reproducibility).** `wave70/modflow/scripts/s1_invariance.py` writes its run time into `s1_out.json`, so `out/SHA256SUMS` changes on every rerun. Drop the field from the hashed copy.

---

## 3. Notes on the steps that hold

**N1 (Lemma 4, item 5: why D_c is a core for e^{πA}).** Let x ∈ Dom(e^{πA}) be orthogonal to D_c in the graph inner product. For d ∈ D_c and g ∈ C_c^∞(ℝ), g(A)d ∈ D_c, so 0 = ⟨x, g(A)d⟩ + ⟨e^{πA}x, e^{πA}g(A)d⟩ = ⟨ḡ(A)(1 + e^{2πA})x, d⟩, where ḡ(A)(1 + e^{2πA}) is bounded. D_c is dense, so ḡ(A)(1 + e^{2πA})x = 0 for every g, hence x = 0. The invariance of D_c under g(A) is true: g·φ̌ = ψ̌ with ψ = ψ₁ + iψ₂, ψ_j real and ψ̌_j ∈ C_c^∞, so g(A)u_φ = u_{ψ₁} + iu_{ψ₂}. Items 1–4 were also checked: the reflection, the bilinear extension, the identity ⟨v, e^{2πA}u⟩ = ⟨u, v⟩ and the isometry of J₀ = e^{πA}S (expanded in real and imaginary parts, four terms each side).

**R1 (the number c rests on factoriality, by a logarithm).** The field aK·1_B (a constant, K cut sharply at the sphere) has no curvature in the open ball and obeys the linearised scalar constraint as a distribution, because K·x is tangent to S and surface-divergence-free. If it were a vector of ℋ it would lie in H(O) ∩ H(O)′, the decomposition of i log Δ k would not be unique, and S would move by a∫_B|K|²/16G. It is not a vector: the radial transform of the jump profile is j_l(k)/k, so its Ḣ^{1/2} norm is ∫k j_l(k)² dk, which diverges like ½ ln k (measured: each decade adds 1.15 ≈ ½ ln 10), while the cut vector's profile N gives j_{l+1}(k)/k² and ∫j_{l+1}(k)²/k dk converges. This agrees with the abstract proof of factoriality and is the one place where the result is sensitive to the norm. One sentence in M-c would pre-empt the question. It is also the edge mode that an algebra violating H2 or H3 could add.

**M3 (independent check of the integer l + 1).** For the stationary solution E = 0, B = −curl K, the conformal Lie derivative along ζ with weight w is, at all times, δE = (boost term) and δB = −t(l − 2 + w)B. It solves the field equations if and only if w = 3, and then δE = (l+1)K on the slice. Weights 2 and 4, and the boost term with the opposite sign, fail (exact, l = 2, 3, 4). So l + 1 = (l − 2) + 3: the degree of the curvature plus its conformal weight, where the metric has l + 2 from (K4).

**N2 (what I could not make fail in (M)).** The sign: S ≥ 0 fixes it, as the text says. The 2π: Lemma 4's strip is −π ≤ Im z ≤ 0 and its toy model rejects the rates 1.1/2π, −1/2π and 1/π. The direction of the flow in Lemma 5: the damping |e^{ip·Λ_{−z}x}| ≤ 1 holds on the lower strip and fails by 10¹¹ on the upper one. The geometry: Dκ·ζ = X_𝒲∘κ with the wedge apex at x¹ = ½, exact (this check is not in `run_all.sh`; it is in my script).

**N3 (single l is essential for (B)).** A mixture of l = 2 and l = 3 drives the translation mode: −½∫_B x_a|K|² ≠ 0 (exact values in `out/p_blind.json`). The theorem's scope line says "any mixture of m at fixed l", which is right. For mixed l, B-d still holds as printed (tr K = 0 and K_nn = 0 are linear in K, so L_K ≡ 0), but B-e fails: the coordinate ball is no longer stationary through O(ε²). Worth one sentence in the Scope paragraph saying that fixed l is used in B-e and nowhere else in (B).

---

## 4. My plants (39, all caught)

| Proof | Plant | Result |
|---|---|---|
| SQ (a) | iJ_HH taken for the symplectic complement | Im⟨·,·⟩ = 3.2 ≠ 0 |
| SQ (c) | KMS continued into the upper strip (β = +1) | mismatch 229 |
| SQ (d) | Δ^{−1/2} in the Tomita operator; Γ(+J_H) as the conjugation | mismatch 245; 23.9 |
| M | Schouten coefficient ½ in Cf; weight 2; weight 4; boost term of opposite sign (six harmonics each) | all fail the field equations |
| M | the edge vector aK·1_B added to the cut vector | shifts S by a·16π/5 (l = 2); norm diverges, ½ ln 10 per decade |
| C | H3 fails (pinching after Φ) | S drops from 0.473 to 0.093 |
| C | H1 fails (another channel for the vacuum) | S rises to 2.08 |
| B | τ-coefficient −1 in Q | boost no longer a zero mode |
| B | mixed l = 2, 3 | translation gradient (−256π/35, 0, 768π/7) |
| B | electric datum K = Hess φ | K_nn ≠ 0, L_K ≠ 0 |
| Lemma 4 | rate 1.1; opposite direction | not isometric (0.31); mismatch 8.6 |
| Lemma 5 | upper half-plane; θ = π/2 | unbounded (4.9e11); complex point |
| Step 6 | wedge apex at x¹ = 1 | Dκ·ζ ≠ boost field |
| Lemma A | non-transverse compactly supported field | non-zero mean, so infinite Ḣ^{−3/2} norm |
| Lemma B | solid torus | toroidal field has period π²/2; curls of supported potentials have period 0 |

The matrix model for (C) has M = M₃ ⊗ M₂, Φ(a ⊗ b) = a⟨e₁, be₁⟩, ι(a) = a ⊗ 1: H1–H3 hold, the support projection is 1 ⊗ e₁₁ ≠ 1, and S_M = S_𝒜 to 1e−10. So the support step is not empty.

---

## 5. Sources fetched

| File in `sources/` | What it proves | Locator |
|---|---|---|
| `sorce-2023-intuitive-construction-of-modular-flow-2309.16766.pdf` (arXiv 2309.16766v2) | T6 and T7, with full proofs | T7: Def. 3.1 pdf p19, Thm 3.2 pdf p20 (proof pp20–23). T6: Thm 4.15 pdf p37, Cor. 4.18 pdf p39, §4.4 pdf pp39–40 (JAJ = A′). Thm 2.13 pdf p15 for F2. |
| `summers-2005-tomita-takesaki-modular-theory-math-ph-0511034.pdf` (arXiv math-ph/0511034) | statements only | Thm 1.1 pdf p2 (T6); §2.1 pdf p2, Thm 2.1 and Thm 2.2 pdf p3 (T7, and the dense-subalgebra form of the KMS condition) |

Both are in `sources/INDEX.md`, new section 6, ids checked against arXiv's API. THEOREM.md's line "these locators are from memory; the books are not on file" can now point at Sorce. Sorce states the KMS condition on the whole algebra, which is why F1 is needed.

---

## 6. What this verdict does not cover

- Sections 9–12 of THEOREM.md, by the blindness rule.
- The physical hypotheses H1–H4, and Proposition R.
- The record's Lemma 1 beyond its statement: I read its proof and found no fault, and its exact script reran, but I did not re-derive the 324-component identity by hand.
- I1b (Araki's lower semicontinuity) and Reible's §III.4: statements read, proofs not.
- The second route to (M) (Paper A's transfer and Longo–Morsella): not read; only `s3`'s arithmetic reran.

## 7. Recommended next lane

**w72-lattice (Fable, with one Sonnet runner): an independent numerical measurement of (l+1)/(l+2).** Both routes to (M) share I1 and Proposition SQ, and every check so far is algebra on the same formulas. Put the odd-parity sector of degree l on a radial lattice in curvature variables, compute the vacuum modular Hamiltonian of the ball from the Gaussian state, and evaluate the relative entropy of the stationary multipole at l = 2, 3, 4 with a continuum extrapolation. Pre-register: the target ratio S·8G/A = 3/4, 4/5, 5/6, a tolerance, and two plants (the canonical-energy value 1, and the edge-vector shift of R1). It also tests R1 directly, since a lattice algebra that keeps the edge mode should show the logarithm.

Before that, the lead applies F1–F9 to THEOREM.md; they are text only and need no lane.

## Files

- `VERDICT.md` (this file)
- `scripts/p_blind.py`, `out/p_blind.json`: my checks and plants
- `run_all.log`: the rerun of `wave69/final/scripts/run_all.sh`
- `sources/INDEX.md` (section 6 added) and the two PDFs above
