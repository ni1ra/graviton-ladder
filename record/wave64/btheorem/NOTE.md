# The area side of the second-order first law of a flat ball is fixed by the first-order first law

Lane w64-btheorem, 9 November 2026. A self-contained note for a hostile referee. Scripts and outputs are in `scripts/` and `out/` (exact sympy; every script runs in under 3 s).

## 1. Setting and the object

Linearised gravity on Minkowski space at leading order in G (the G⁰ counting of JSS, arXiv 2306.01837). Σ = {t = 0}; B the unit ball; S = ∂B; D(B) its diamond; ζ = ½(1 − t² − r²)∂_t − t x·∂ its conformal Killing vector, with ζ|_Σ = N∂_t, N = (1 − r²)/2. Cauchy data of a metric perturbation are (q, p) = (δh_ij, δK_ij), smooth on a neighbourhood of B̄. Write

- R_lin[q] = ∂_i∂_j q_ij − ∇²q_kk (the linearised scalar curvature of the slice),
- δA[q] = ½∮_S (q_kk − n_i n_j q_ij) dA, δV[q] = ½∫_B q_kk dV, with n = x on S,
- gauge: (q, p) → (q + ∂_iξ_j + ∂_jξ_i, p + Hess f).

Expand a one-parameter family of configurations as g_ε = η + εh⁽¹⁾ + ε²h⁽²⁾ + O(ε³). The **second-order first law** compares an *area side* (geometric) with a *modular side* (the coherent-state relative entropy of the field net) at O(ε²). On the scope family of the ladder, h⁽¹⁾ = 2tK with K transverse, traceless and harmonic (a magnetic multipole, l ≥ 2), so h⁽¹⁾|_Σ = 0 and the first-order data are (q₁, p₁) = (0, K).

The free-net relative entropy is an exact quadratic form in the first-order datum, so h⁽²⁾ can enter the comparison only through the gravitating, geometric part. The claim to prove, in **total form**:

> **(BR₄)** The total part of the second-order first law (area side minus modular side) that is linear in h⁽²⁾ equals −(1/4G)·δ(A − 2V)[q₂] = (1/8G)∫_B N R_lin[q₂], where q₂ is the Cauchy datum of h⁽²⁾. On the scope family this equals A/8G, A = ∫_B N|K|², i.e. "A" in the units of Paper A (which drop 1/8G).

Clause (3) of the LADDER-WHY truth condition says that every term of the second-order modular side linear in h⁽²⁾ is δ(A − kV), so that the area side is A. LADDER-PHYS uses (3) only through "the area side is A", which (BR₄) gives. Given clause (2) (the physical modular side is the free relative entropy, so it has no h⁽²⁾), (BR₄) is clause (3) verbatim. Without (2), a gravitational term τΨ[h⁽²⁾] on the modular side is allowed by (BR₄) only together with −τΨ[h⁽²⁾] on the area side. That is bookkeeping between B and C, not a change in the ladder.

## 2. Premises

**(P0) One Taylor functional.** At G⁰ the first-law defect is 𝔉(ε) = Φ[g_ε] − 𝒮(ε). Here Φ is a twice-differentiable functional of the metric's Cauchy data on the (relationally fixed) region, independent of the state, and 𝒮 is the modular-side functional of the field state. Write 𝒟 := DΦ|_η, a linear functional on all smooth (q, p) (no continuity or locality is assumed). Lemma 7 shows that h⁽²⁾ then enters 𝔉 at O(ε²) only as 𝒟(q₂, p₂).

**(T) Time reflection.** t → −t is an isometry of η that preserves D(B) and the vacuum Ω (antiunitarily). The first-law functional is invariant under it, as entropies and relative entropies are invariant under antiunitary symmetries. Since q → q and p → −p, this gives 𝒟(q, p) = 𝒟(q, −p).

**(U) Universality.** Adding a free conformally coupled scalar does not change 𝒟 on configurations where the scalar's field value vanishes on Σ. Its only correction to the gravitational entropy functional is the Wald term −2πξ∮_S φ², which then vanishes (Lemma 3).

**(FL1) First-order first law with respect to the vacuum, on the physical slice 𝒮.**
- *The states.* Let ω_g be the coherent state of the conformal scalar with Cauchy data (φ, φ̇) = (0, g), g smooth, and ρ_λ = (1 − λ)Ω + λω_g.
- *The law.* For every first-order metric datum (q, p) ∈ 𝒮 solving the linearised constraints with source λ⟨T_μν⟩_{ω_g}, the first-order area side equals the first-order vacuum modular energy:
  𝒟(q, p) = λ·S(ω_g‖Ω) = 2πλ∫_B N·½g².
- *The slice.* 𝒮 is either the space of all data (a relationally fixed region, "corner diffeomorphisms are gauge"), or the kernel 𝒮_q × 𝒮_p of the linearisation of a physical edge rule, such as JSS's fixed-volume condition (H_ξ^Σ)^(−1) = 0, §3.3. A slice is admitted if it meets every gauge orbit and contains p = 0.

FL1 is the first-order content of JSS's first law (5.35) (their assumption A5), of the Bueno–Min–Speranza–Visser first law of causal diamonds (1612.04374), and of Jacobson's entanglement equilibrium (1505.04753). Those sources show that the Iyer–Wald/JSS functional satisfies FL1 (existence; Lemma 6 re-checks it). FL1 is used here as a property of the **unknown** gravitating functional.

**(Q0) No K-quadratic geometric term.** Φ's quadratic part vanishes on first-order data with q₁ = 0. It holds for any functional of the spatial metric alone (JSS's A/4G and V). The Einstein–Hilbert Wald entropy has no K dependence. At Σ, ∇·ζ = 0 removes the JKM and Y ambiguities (w52-b5 S3). The Noether charge of ζ at S is algebraic in g there (w39-BR S1). Q0 is used only by Corollary 2. *Added after the PREREG*: Corollary 2 needs it, and the PREREG did not register it.

**What the theorem is, stated plainly.** FL1 quantifies over every datum in 𝒮 with a given source, and the linearised constraint makes that source R_lin[q]/16πG. So FL1 *is* first-order (BR), stated on a spanning set of data, and Lemma 6 is its linear extension to all of 𝒮. This extension covers all angular components, the corner jets and the corner-displaced data, and needs no rotation invariance, corner gauge, locality or jet bound. The result is therefore:

> first-order (BR) [printed physics: BMSV; JSS §3.3 and (5.35); Wald] + one Taylor functional (P0) + (T, U) ⇒ second-order (BR), and with (Q0), "the area side is A".

The open part of B was always the second-order transfer: whether the coefficient of h⁽²⁾ in the gravitating first law is the first-order functional. That transfer is Lemma 7. The computations below check identities and the algebra of the extension. They pass identically in every world, so they carry no power on B. The evidence is the reduction itself.

## 3. Theorem

**Theorem 1.** Under (P0), (T), (U) and (FL1):

    𝒟(q, p) = −(1/4G)·δ(A − 2V)[q] = (1/8G)∫_B N R_lin[q]    for all (q, p) with q ∈ 𝒮_q.

Hence, by Lemma 7 and Lemma 8, (BR₄) holds: the h⁽²⁾-linear part of the second-order first law is (1/8G)∫_B N R_lin[q₂], up to terms quadratic in first-order data that come from the edge rule's second-order part.

**Corollary 2 (scope family).** Assume also (Q0). Then the second-order area side on h⁽¹⁾ = 2tK is exactly A/8G, with A = ∫_B N|K|². [COMPUTED: 16π/35 and 128π/105 for l = 2, 3, with ratio 1/((2l+1)(2l+3)) to the flux ∮|K|², matching Paper A's 1/35 at unit flux; `out/s9_scope.json`.]

## 4. Lemmas and proofs

**Lemma 1 (Green identity).** For smooth w and symmetric q on B̄,

    ∫_B w R_lin[q] = ∫_B q_ij(∂_i∂_j w − δ_ij∇²w) + ∮_S [w n_i∂_jq_ij − n_j q_ij ∂_i w − w ∂_n q_kk + q_kk ∂_n w].

If w|_S = 0 and ∂_i w|_S = −κ n_i, the boundary term is −2κ δA[q]. For w = N (κ = 1, ∂_i∂_jN = −δ_ij, ∇²N = −3):

    ∫_B N R_lin[q] = −2δA[q] + 4δV[q] = −2δ(A − 2V)[q]     for every smooth q.

*Proof.* Integrate by parts twice in each term of R_lin. With w = 0 on S only the terms −n_j q_ij ∂_i w + q_kk ∂_n w survive, giving κ(q_nn − q_kk) = −κ(q_kk − q_nn). ∎ [COMPUTED on three random degree-4 tensors, all six components, plus the general-w form with w = N + (3/7)(1 − r²)². The plant k = 2.2 is caught. `out/s1s2_green_gauge.json`]

**Lemma 2 (corner invariance).** For every smooth ξ on B̄, with no condition at S, δ(A − 2V)[∂_iξ_j + ∂_jξ_i] = 0, while δA and δV separately are not invariant.

*Proof.* δV = ∫_B ∂_kξ_k = ∮ ξ·n. δA = ∮(∇·ξ − n·∂_nξ) = ∮ div_S ξ = ∮(div_S ξ_∥ + H ξ_n), with mean curvature H = 2. The first term integrates to zero on a closed surface. So δA = 2∮ξ_n. ∎ This holds because k = 2 is the mean curvature of S: A − kV is stationary under normal displacements exactly at the CMC radius. [COMPUTED on four random cubic ξ and on the dilation (δA = 8π, δV = 4π). The plants are caught: k = 2.2 fails, and at radius 2 the invariant is A − V, so A − 2V fails.]

**Lemma 3 (the matter data).** For the conformally coupled scalar (ξ = 1/6) with data (φ, φ̇)|_Σ = (0, g), the improved stress tensor on Σ is T₀₀ = ½g² and T₀ᵢ = 0, and the Wald term ∝ ∮φ² vanishes. [COMPUTED from the covariant T_μν with φ a sympy `Function`; the trace vanishes on shell; `out/s3_span.json`]

**Lemma 4 (spanning).** Every smooth ρ on B̄ is a real linear combination of profiles ½g² with g smooth on R³. *Proof.* Extend ρ smoothly and take χ ∈ C_c^∞ with χ ≡ 1 near B̄. Then ½(χ + ρ)² − ½(χ − ρ)² = 2ρ on B̄. ∎

**Lemma 5 (first law for mixtures) [proof from CITED inputs; the von Neumann form of Donald's identity was not re-read this session].** Let ω, σ be normal states on a von Neumann algebra with S(ω‖σ) < ∞, and ρ_λ = (1 − λ)σ + λω. Then S(ρ_λ‖σ) = o(λ).

*Proof.*
- Donald's compensation identity (M. J. Donald, CMP 105 (1986) 13 and Math. Proc. Camb. Phil. Soc. 101 (1987) 363; cf. Ohya–Petz, *Quantum Entropy and Its Use*, Ch. 5) gives λS(ω‖σ) = S(ρ_λ‖σ) + λS(ω‖ρ_λ) + (1 − λ)S(σ‖ρ_λ). All terms are finite: ρ_λ ≥ (1 − λ)σ gives S(σ‖ρ_λ) ≤ −log(1 − λ), and joint convexity gives S(ω‖ρ_λ) ≤ (1 − λ)S(ω‖σ).
- Hence S(ρ_λ‖σ)/λ ≤ S(ω‖σ) − S(ω‖ρ_λ).
- As λ → 0, ρ_λ → σ in norm. Lower semicontinuity of Araki's relative entropy gives lim inf S(ω‖ρ_λ) ≥ S(ω‖σ).
- So lim sup S(ρ_λ‖σ)/λ ≤ 0. ∎

The modular side, a relative entropy, therefore has no first-order term, and the first-order modular energy of ρ_λ is λ·S(ω_g‖Ω). For coherent states the reduced state is a local displacement of the vacuum's, so ΔS = 0 and Δ⟨K⟩ = S. Coherent families cannot test a first law, because their relative entropy starts at O(amplitude²); the mixtures ρ_λ can. For data (0, g), Longo–Morsella's formula (arXiv 2012.00565; the form used in Paper A §1 and Cor. 2(iii)) gives S(ω_g‖Ω) = π∫_B N g² = 2π∫_B N T₀₀. Any boundary or "Born" term is ∝ ∮φ² and vanishes. [CITED, not re-read this session, including its validity for data not localised in B]

**Lemma 6 (linear extension of FL1).** Assume (T) and (FL1). Then 𝒟(q, p) = (1/8G)∫_B N R_lin[q] for q ∈ 𝒮_q.

*Proof.*
1. By (T) and linearity, 𝒟(0, p) = 0, so 𝒟(q, p) = 𝒟(q, 0).
2. Fix q ∈ 𝒮_q and put ρ = R_lin[q]/16πG, which is smooth on B̄. By Lemma 4, ρ = ½·(½g₊²) − ½·(½g₋²) with g_± = χ ± ρ.
3. For each sign choose q_± ∈ 𝒮_q with R_lin[q_±] = 16πG·½g_±². This is possible: q = uδ with −2∇²u = source solves it, and a gauge term, which leaves R_lin unchanged, moves it into 𝒮_q.
4. The data (q_±, 0) satisfy the momentum constraint with source T₀ᵢ = 0 (Lemma 3), so FL1's O(λ) coefficient gives 𝒟(q_±, 0) = 2π∫N·½g_±².
5. q₀ := q − ½q₊ + ½q₋ ∈ 𝒮_q is source-free. Since q₊ + q₀ is another datum in 𝒮 with the source of g₊, FL1 gives 𝒟(q₀, 0) = 0.
6. Linearity then gives 𝒟(q, 0) = 2π∫_B N ρ = (1/8G)∫_B N R_lin[q]. By Lemma 1 this equals −(1/4G)δ(A − 2V)[q]. ∎

Nothing in the proof uses rotation invariance, corner gauge, locality, a jet bound, the graviton's modular flow or any KMS state for the graviton. [COMPUTED in `out/s10_power.json` (units 1/4G = 1). This illustrates the extension; given FL1 ⇔ 𝒟 = −δ(A − 2V) on 𝒮 it is tautological, and it has no power on B.
- −δ(A − 2V) satisfies FL1 on every configuration tried: two matter profiles, one with ρ ≠ 0 on S and one angle-dependent; an even-parity TT source-free datum; the dilation and an l = 2 corner displacement; a matter datum plus TT plus gauge; and the JSS fixed-volume slice.
- Each rival fails FL1 on at least one of them:
  - an angle-dependent area charge ∮Y₂₀q_AA;
  - the corner monopole δA;
  - the weight N + (1 − r²)²;
  - a corner jet ∮R_lin|_S;
  - the normalisation 11/10.
- JSS's own functional −δA passes on its slice and fails off it, as it should.
- A T-odd term ∫p_kk is invisible to FL1 and is excluded by (T) alone.]

**Lemma 7 (Taylor transfer; P0's content).** Expand 𝔉(ε). Then h⁽²⁾ appears at O(ε²) only in 𝒟(q₂, p₂). The three ways it could appear otherwise all fail:
- (i) **Through the modular side.** The coherent-state relative entropy of a free net is exactly a quadratic form Q in the datum (Ciolli–Longo–Ruzzi; Paper A Thm 1). Including a second-order displacement gives Q(εd₁ + ε²d₂) = ε²Q(d₁) + O(ε³).
- (ii) **Through a data-dependent region, frame or flow** R_ε = R₀ + εR₁ + …. In Φ[g_ε; R_ε], h⁽²⁾ multiplies only the background derivative at O(ε²). Cross terms such as ∂_R DΦ[h⁽¹⁾]·R₁ are quadratic in first-order data. This is the order-counting spine of w51-b4.
- (iii) **Through matter-loop dependence of Φ on the metric** (vacuum entanglement terms). This is O(G⁰)·h⁽²⁾ with h⁽²⁾ = O(G) data, which is O(G) relative to the 1/G area term and drops out at the ladder's order. ∎ (i) is cited, (ii) is proved, (iii) is the G counting.

**Lemma 8 (the edge rule's second-order part).** Let the physical edge rule be E[g] = 0, with E twice differentiable, E[η] = 0 and L = DE[η] surjective. Then 𝒮 = ker L, and at O(ε²), L(q₂, p₂) = −½D²E[(q₁, p₁), (q₁, p₁)]. With a right inverse L⁺, (q₂, p₂) = (element of 𝒮) − L⁺(½D²E[·,·]). So 𝒟(q₂, p₂) = (1/8G)∫N R_lin[q₂] plus a term quadratic in first-order data. For JSS's fixed-volume rule E = V − V₀, the quadratic part ½D²V[q₁, q₁] vanishes on the scope family (q₁ = 0), so q₂ ∈ 𝒮 exactly. ∎

*Proof of Theorem 1.* Lemma 6 gives 𝒟 on 𝒮. Lemmas 7 and 8 carry it to h⁽²⁾. Corollary 2 follows from Lemma 1 and the O(ε²) Hamiltonian constraint ³R + (tr K)² − K_ijK^ij = 0. With γ = δ + ε²q₂ and K_ij = εK + O(ε²), it reads R_lin[q₂] = |K|². So (1/8G)∫N R_lin[q₂] = A/8G. [COMPUTED: the exact Ricci scalar of δ + ε²q has no O(ε) term and its O(ε²) term is R_lin[q]. An explicit q₂ = uδ gives −2δ(A − 2V)[q₂] = A at l = 2, 3. Adding harmonic and gauge terms leaves it unchanged, while δA alone moves.] ∎

## 5. What this settles about (c3), and what it does not use

- **Corner gauge (c3) is not assumed.**
  - *Relational region.* FL1 must hold on corner-displaced configurations too, since they are first-order data with the same source. Lemma 6 gives a functional that does annihilate every corner displacement (Lemma 2). This is automatic for the round unit sphere, because k = 2 is its mean curvature, so it is not a test. It is why the relational and sliced readings agree.
  - *Physical edge rule.* 𝒟 is pinned on the physical slice, and h⁽²⁾ lies on it up to quadratic terms (Lemma 8).
  - Either way, the five corner monopoles that w39-BR S2 found invariant without corner gauge (δA among them) are excluded from the linear clause by FL1 (planted δA caught in `s10_power`), or, on a slice where one of them coincides with the right functional (δA on JSS's fixed-volume slice), are identically harmless. This is an argument, not new evidence.
  - What is left of the edge question is quadratic in first-order data, which is clause (2)/C.
- **Not used:**
  - SO(3) invariance of the generator (record condition (a)): FL1 with non-spherical matter pins every angular component.
  - The lapse weight from graviton Bisognano–Wichmann, or Takesaki on the curvature net (condition (c)/(c′)): the weight N comes from matter, whose vacuum modular flow is geometric (Hislop–Longo).
  - JSS's A4 for the graviton, the state Ψ_ζ and the converse-cocycle lemma (w60-cprime, w61-cocycle).
  - Any truncation in degree or corner jets.
- **If one insists on JSS's construction anyway,** w61's edge-regularised cure replaces N by w_f only in a shell of width δ at S, with |w_f − N| ≤ C(1 − r) there. So |∫(w_f − N)R_lin[q₂]| ≤ C‖R_lin[q₂]‖_∞·4π∫_{1−δ}^1 (1 − r)dr = O(δ²), which vanishes as the regulator is removed.

## 6. Cross-check by the record's route (exact)

The record's route assumes (a) SO(3) invariance and (c3) corner gauge, and classifies invariants. That classification is now exact with no truncation [COMPUTED with sympy `Function` profiles; `out/s11_l0.json`].
- *Reduction.* For l = 0 data h = aδ + b nn and K = cδ + d nn, the smooth gauge maps ξ = −½a(r)x and f = −∫₀^r s c(s) ds remove a and c for arbitrary profiles. The invariants are β = b − ra′ and γ = d − rc′.
- *Constraints.* R_lin = (2/r²)(rβ)′ and n·div π = 2γ/r. Regularity inverts both, β = (2r)⁻¹∫₀^r s²R_lin, because the homogeneous solution C/r is singular.
- *Conclusion.* Every linear functional F (of any jet order, local or not) that is SO(3)-invariant and invariant under all smooth gauge transformations satisfies F(q, p) = F(0, β nn, 0, γ nn). It is therefore a functional of (R̄_lin, n·div π̄) alone, i.e. a smeared constraint.
- FL1 restricted to radial matter then pins its weight to N, and (T) removes the momentum part.
- This reproduces Theorem 1 on the record's premises and closes the truncation of w39-BR S2 and w60-br6 T2.

## 7. Where a referee should push

1. **FL1 is a premise.**
   - It is the first-order first law, read as a property of whatever the gravitating area side is.
   - Theorem 1 says that only (A − 2V)/4G has it.
   - If FL1 failed (no first-order equilibrium relative to the vacuum; for example the ladder's reference state is not Ω, or first-order edge terms survive), the theorem says nothing.
2. **P0, "one Taylor functional".**
   - Lemma 7 closes the three routes by which h⁽²⁾ could enter otherwise at the ladder's order.
   - A referee may still ask for a definition of the G⁰ gravitating second-order first law that is independent of JSS's crossed product. The ladder's comparison (geometric functional against CLR relative entropy) is such a definition.
3. **(Q0), Corollary 2's extra hypothesis.** A geometric term quadratic in K at the corner would be an edge term at O(ε²). That is clause (2)'s territory, not (3)'s.
4. **The weakest place is P0 with Q0, not FL1.** FL1 is printed first-order physics. A world where B fails while FL1 holds must have a second-order gravitating area side that is not the Taylor continuation of one fixed metric functional, or that has a K-quadratic corner term.
