# REPORT: w70-modflow. The modular group of the ball in the free helicity-2 curvature net is the conformal flow, with a proof on file

Lane w70-modflow (Fable), Tue 10 Nov 2026, Oslo. PREREG.md hashed 10:13 (`PREREG.sha256`), before any script existed.

**Bin: PROOF** for the one-particle statement (the modular group of the standard subspace K(O), which is what Paper A uses). On the strictest reading of the PREREG wording, PROOF-WITH-CONDITION, for two reasons that I state rather than argue away:
- the step from K(O) to the von Neumann algebra A(O) imports the second-quantisation theorem (statement read, proof not read), and the setting is wave69's, whose own weakest step 2 is inherited;
- three registered items missed as registered (§4): the boost sign I registered was wrong by a convention (corrected, and the registered sign is now a plant that fails); a numerical tolerance was missed by my integrator (fixed post hoc, as-registered output kept); one monotonicity clause of the lattice test fails at the 1e−4 level.

Nobody else has checked the proof. It uses neither Hislop 1988, nor Bisognano–Wichmann's papers, nor a representation of the conformal group, nor Paper A's transfer to scalar waves.

## 1. The theorem

**Setting** (wave69/duality/REPORT.md §1, with one sign made explicit).
- H = TT Cauchy data u = (h, π); σ(u, u′) = ∫(h·π′ − h′·π); E[u] = −½Δh, B[u] = ½ curl π; Ė = curl B, Ḃ = −curl E; T_t the time evolution of data.
- **Complex structure.** I(h, π) = (−ω⁻¹π, ωh), which is −J in wave69's notation. In the wave function χ_u(k) = ω^{1/2}ĥ(k) + iω^{−1/2}π̂(k) it is multiplication by i, and ⟨u, u′⟩ = ∫ conj(χ_u)·χ_{u′} d³k/(2π)³ has Im⟨u, u′⟩ = σ(u, u′) and wave69's norm. The sign is forced by positive energy: the observable "E at time t smeared with f" is the functional u ↦ σ(u, T_{−t}v_f), so moving an observable forward by a acts on its vector by T_{−a}, which is e^{+iωa} on χ. With +J the energy is negative. wave69's theorem (a statement about symplectic complements) does not see this sign; the direction of a modular flow does.
- **Vectors of smeared curvature.** For a real test field F = (f, g) on spacetime with values in S²₀ ⊕ S²₀, e_F := ∫dt T_{−t}(v_{f_t} + w_{g_t}), so that σ(u, e_F) = ∫d⁴x (f·E + g·B)[T_t u] =: C_u(F), the curvature of u smeared with F. K(U) := closed real span of {e_F : supp F ⊂ U}. This is wave69's K(O).
- **Conformal action.** For a conformal Killing field X (∂_μX_ν + ∂_νX_μ = 2λη_μν, λ = ∂·X/4), with A_im = ∂_iX^m − λδ_im and β_k = ∂_kX⁰,

  (δ_X E)_ij = X^μ∂_μE_ij + A_im E_mj + A_jm E_im + 3λE_ij + (ε_ikl β_k B_lj + ε_jkl β_k B_li),
  (δ_X B)_ij = X^μ∂_μB_ij + A_im B_mj + A_jm B_im + 3λB_ij − (ε_ikl β_k E_lj + ε_jkl β_k E_li).

 Write δ_X = X·∂ + M_X(x) and δ_X†F := −∂_μ(X^μF) + M_XᵀF. The push-forward of a test field along the flow κ_τ of X is the solution F_τ = κ_{τ*}F of ∂_τF_τ = δ_X†F_τ, F₀ = F; explicitly (κ_{τ*}F)(κ_τx) = J_τ(x)⁻¹ C_τ(x) F(x), J the Jacobian and dC/dτ = M_X(κ_τx)ᵀC.

**Theorem.** Let O_R = D(B_R) be the double cone over the ball of radius R centred at the origin of the slice t = 0, and

  ζ_R = (1/2R)[(R² − t² − r²)∂_t − 2t x·∇], equal to N_R ∂_t at t = 0, N_R = (R² − r²)/2R,

with flow Λ_s (it maps O_R onto itself for every s).
- **(a)** K(O_R) is a standard subspace and Δ_{K(O_R)}^{it} e_F = e_{Λ_{−2πt*}F} for every real F ∈ C_c^∞(O_R).
- **(b)** For every u ∈ H and F ∈ C_c^∞(O_R): C_{Δ^{it}u}(F) = C_u(Λ_{2πt*}F). The curvature of Δ^{it}u inside O_R is the curvature of u inside O_R transported along Λ_{2πt} with weight 3, and d/dt C_{Δ^{it}u} = 2π δ_ζ C_{Δ^{it}u} there. On the slice t = 0:

  δ_ζE = N_R curl B − (2/R) sym(x × B), δ_ζB = −N_R curl E + (2/R) sym(x × E), (x × T)_ij = ε_ikl x_k T_lj.

 With Paper A's magnetic field B_A = −B this is Paper A's eq. (flow).
- **(c)** With wave69's C1, K(O_R) = N(B̄_R): the group acts on the finite-energy data whose E and B vanish outside the closed ball, and it is the flow with speed (R² − r²)/2R at the slice.
- **(d)** Given the second-quantisation theorem, the vacuum modular group of A(O_R) is Γ(Δ^{it}); on observables σ_t(Φ(e_F)) = Φ(e_{Λ_{−2πt*}F}): Hislop–Longo's geometric flow, run backwards at 2π per unit of modular time, exactly as for the scalar.

The form of (a) and (d) does not depend on conventions. The sign in (b) does: it flips if the vectors are identified with solutions through −σ, or the Hilbert space is conjugated, which is how Paper A (after Longo–Morsella) writes Δ^{it} = U(Λ_B(−2πt)) on waves.

## 2. The proof

All numbered script checks are exact computer algebra unless marked numeric.

**Step 0. Wave functions and the two-point function.** From wave69 Step 1, v_f = (0, −½ΔP_TT f) and w_g = (−½ curl P_TT g, 0), so χ_{v_f} = (i/2)ω^{3/2}P f̂ and χ_{w_g} = −(i/2)ω^{3/2} k̂ × P ĝ, with P = P_TT(k̂) and (k̂ × T)_ij = ε_iab k̂_a T_bj. Since T_{−t} is e^{iωt},

  χ_{e_F}(k) = 𝒩(p)F̃(p), 𝒩(p)(f, g) = (i/2)ω^{3/2}P(k̂)[f − k̂ × g], F̃(p) = ∫d⁴x e^{ip·x}F(x), p = (ω, k), p·x = ωt − k·x.

Hence W(F, G) := ⟨e_F, e_G⟩ = ∫∫F(x)ᵀ𝒲(x − y)G(y) with 𝒲(ξ) = ∫d³k/(2π)³ 𝒩†𝒩 e^{−ip·ξ}. The entries of ω⁻¹·2𝒩†𝒩 are polynomials Q(ω, k) of degree 4 on the cone (S1), so

  𝒲(ξ) = lim_{ε↓0} 𝒲_rat(ξ − iεe₀), 𝒲_rat(z) = −(1/4π²) Q(i∂_t, −i∇)(1/z²), z² = t² − |x|²,

a matrix of rational functions with denominator (z²)⁵. S1a checks on 𝒲_rat: 𝒲^{BB} = 𝒲^{EE}, 𝒲^{BE} = −𝒲^{EB}; each index pair symmetric, traceless and transverse; Ė = curl B and Ḃ = −curl E in the first slot; and 𝒲_rat,ab(ξ) = 𝒲_rat,ba(−ξ).

*Locality and Huygens.* Im W(G, F) = (1/2i)∫∫Gᵀ[𝒲(x − y) − 𝒲ᵀ(y − x)]F, and by the last identity the bracket is the difference of the two boundary values of one rational function. It vanishes where that function is regular: for supports that are spacelike or timelike separated.

**Step 1. Lemma 1 (infinitesimal conformal invariance).** For every conformal Killing field X and all F, G ∈ C_c^∞(ℝ⁴): W(δ_X†F, G) + W(F, δ_X†G) = 0.

Proof. S1b: (δ_X at x on the first slot + δ_X at y on the second slot) 𝒲_rat(x − y) = 0 identically, as a polynomial identity in (x, x − y) after clearing (ξ²)⁶, for dilation, three boosts, three rotations and four special conformal fields (translations are trivial): 324 components each, all zero. The identity therefore holds for complex x. Put x′ = x − iεe₀ and 𝒲_ε(x, y) = 𝒲_rat(x′ − y). Then (δ_x ⊗ 1 + 1 ⊗ δ_y)𝒲_ε = [(X^μ(x) − X^μ(x′))∂_μ + M_X(x) − M_X(x′)]𝒲_ε = ε R_ε 𝒲_ε, with R_ε a first-order operator whose coefficients are polynomials in (x, ε). Integrating against F ⊗ G and by parts, the left side tends to W(δ†F, G) + W(F, δ†G), and the right side is ε times a quantity that converges (to W(R₀†F, G)), so it tends to 0. ∎

**Lemma 1′.** [δ_X, δ_Y] = δ_{[X,Y]} for all pairs of the 15 generators (S1c, 105 pairs). Taking adjoints, [δ_Y†, δ_X†] = δ†_{[X,Y]}.

**Step 2. Lemma 2 (finite covariance along a flow).** Let κ_τ be the flow of X, and let U be open with κ_τ smooth on U for τ ∈ [0, τ₁]. Then for F, G ∈ C_c^∞(U):
- (i) W(κ_{τ*}F, κ_{τ*}G) = W(F, G) for τ ∈ [0, τ₁];
- (ii) for every conformal Killing Y, κ_{τ*}δ_Y†F = δ†_{Y_τ}κ_{τ*}F, where Y_τ is the push-forward of Y by κ_τ (again a conformal Killing field, Ẏ_τ = −[X, Y_τ]).

Proof. (i) τ ↦ F_τ is C¹ into C_c^∞ of a fixed compact set and W is continuous, so d/dτ W(F_τ, G_τ) = W(δ_X†F_τ, G_τ) + W(F_τ, δ_X†G_τ) = 0 by Lemma 1. (ii) G(τ) := δ†_{Y_τ}κ_{τ*}F − κ_{τ*}δ_Y†F has G(0) = 0 and G′ = δ†_{Ẏ_τ}F_τ + [δ†_{Y_τ}, δ_X†]F_τ + δ_X†G = δ_X†G by Lemma 1′; a transport equation has one solution. ∎

For a flow defined on all of ℝ⁴ (Poincaré transformations, dilations) (i) holds for all test fields, so e_F ↦ e_{κ_*F} is a unitary of H mapping K(U) onto K(κU), and by (ii) it carries the flow of Y to the flow of κ_*Y. S2b checks the pointwise form of (i), C(x)ᵀ𝒲_rat(κx − κy)C(y) = 𝒲_rat(x − y), numerically for the map κ of Step 6.

**Step 3. Lemma 3 (one-particle Reeh–Schlieder).** For every non-empty open U ⊂ ℝ⁴ the complex span of {e_F : supp F ⊂ U} is dense in H.

Proof. Let x be orthogonal to it, and let F₀ be supported in a ball whose small translates stay in U. g(a) := ⟨x, e_{F₀(·−a)}⟩ = ∫ conj(χ_x)·𝒩F̃₀ e^{ip·a} d³k/(2π)³ is continuous on Im a ∈ V̄₊, analytic for Im a ∈ V₊ (the integrand is in L¹ and |e^{ip·a}| ≤ 1), and vanishes for small real a. For fixed small real spatial a, the function of a⁰ is analytic in the upper half-plane, continuous on its closure and zero on a real interval, hence zero. At a⁰ = iη, η > 0, the function of spatial a is real-analytic and zero near the origin, hence zero on ℝ³. Then for every real spatial a the function of a⁰ vanishes on the imaginary axis, hence on the real axis. So x ⊥ e_F for every translate, and by a partition of unity for every F ∈ C_c^∞(ℝ⁴). The closed span of all e_F contains every v_f and w_g (wave69 Step 7), whose real span L(ℝ³) has L(ℝ³)′ = N(∅) = {0} (wave69 Step 2), so it is H. Hence x = 0. ∎

**Step 4. Lemma 4 (a KMS criterion for standard subspaces; proved here).** Let K be a standard subspace, V(s) = e^{isA} a unitary group with V(s)K = K, and K₀ ⊂ K a dense real-linear V-invariant subspace. Suppose that for all u, v ∈ K₀ the function s ↦ ⟨v, V(s)u⟩ extends to a bounded continuous function f on −π ≤ Im z ≤ 0, analytic inside, with f(s − iπ) real for all real s. Then V(s) = Δ_K^{is/2π}.

Proof.
1. *Reflection.* f(z) := conj f(z̄ − 2πi) extends f analytically to −2π < Im z < 0, with f(s − 2πi) = ⟨V(s)u, v⟩. By the three-lines theorem |f| ≤ ‖u‖‖v‖ on the closed strip.
2. *All of K.* f depends real-bilinearly on (v, u) and is bounded by ‖u‖‖v‖, so it extends to u, v ∈ K by uniform limits, with the same properties.
3. *Entire vectors.* For real φ with φ̌ ∈ C_c^∞, u_φ := ∫φ(s)V(s)u ds = φ̌(A)u lies in K and is entire for A. ⟨v, e^{izA}u_φ⟩ is entire and equals f on the real line, hence on the strip. At z = −2πi: ⟨v, e^{2πA}u_φ⟩ = ⟨u_φ, v⟩ for all v ∈ K. Let K_c be the real span of the u_φ and D_c = K_c + iK_c. K_c is dense in K; D_c is invariant under g(A) for g ∈ C_c^∞(ℝ) and under e^{±πA}.
4. *An antiunitary.* With S = S_K (S(u₁ + iu₂) = u₁ − iu₂), the identity of item 3 gives, for x, y ∈ D_c, ⟨e^{πA}Sy, e^{πA}Sx⟩ = ⟨Sy, e^{2πA}Sx⟩ = ⟨x, y⟩ (expand in real and imaginary parts). So J₀ := e^{πA}S on D_c is antilinear and isometric with range e^{πA}D_c = D_c, hence extends to an antiunitary. It commutes with V(s), so J₀e^{πA} = e^{−πA}J₀, and S = J₀e^{πA} on D_c.
5. *Cores.* The graph norm of S on K + iK is √2 times the norm of K ⊕ K, so D_c is a core for S. D_c is dense and invariant under g(A), g ∈ C_c^∞, so it is a core for e^{πA}. Two closed operators that agree on a common core are equal: S = J₀e^{πA}. Uniqueness of the polar decomposition gives Δ_K^{1/2} = e^{πA}. ∎

(The statement is the standard one-particle KMS uniqueness; I found no printed locator in the files on hand, so the proof is given.)

**Step 5. Lemma 5 (one-particle Bisognano–Wichmann for this net).** Let W_R = {x¹ > |t|}, X = x¹∂_t + t∂₁ with flow Λ_τ (future-directed in W_R). Then K(W_R) is standard and Δ_{K(W_R)}^{it} e_F = e_{Λ_{−2πt*}F} for real F ∈ C_c^∞(W_R).

Proof. X is a Killing field with constant M = M_X, so Λ_{τ*}F(y) = e^{τMᵀ}F(Λ_{−τ}y), and by Lemma 2 V(s)e_F := e_{Λ_{−s*}F} defines a unitary group with V(s)K(W_R) = K(W_R). K(W_R) is standard: cyclic by Lemma 3; separating because K(W_L) ⊂ K(W_R)′ by locality and K(W_L) is cyclic.

For z = s − iθ, 0 ≤ θ ≤ π, put F̃^{(z)}(p) := e^{−zMᵀ}∫e^{ip·Λ_{−z}x}F(x)d⁴x and χ^{(z)} := 𝒩F̃^{(z)}.
- Λ_{−z} = Λ_{−s}Λ_{iθ} and Λ_{iθ}x = (t cos θ + i x¹ sin θ, x¹ cos θ + i t sin θ, x²,x³). With q = Λ_s p (real, on the cone), |e^{ip·Λ_{−z}x}| = e^{−sin θ (q⁰x¹ − q¹t)} ≤ 1 for x ∈ W_R.
- The phase ψ(x) = q·Λ_{iθ}x is linear in x with |∇ψ|² = |q|²_Euclid, so N integrations by parts give |∫e^{iψ}F| ≤ C_N(F)(1 + |q|)^{−N}, uniformly in θ ∈ [0, π].
- M is diagonalisable with eigenvalues 0, ±1, ±2 (S5), so e^{−zMᵀ} is entire and bounded for z in compacts.

So z ↦ χ^{(z)} ∈ H is continuous on the closed strip and analytic inside, and χ^{(s)} is the wave function of V(s)e_F. For fixed s, V(s)χ^{(z)} and χ^{(z+s)} are analytic and agree on the real line, hence in the strip, so ‖χ^{(s−iθ)}‖ = ‖χ^{(−iθ)}‖ is bounded. Therefore f(z) := ⟨e_G, χ^{(z)}⟩, G real in W_R, is bounded, continuous on the closed strip, analytic inside, and f(s) = ⟨e_G, V(s)e_F⟩.

At θ = π, Λ_{iπ} = j := diag(−1, −1, 1, 1) and e^{iπMᵀ} is the real matrix (E, B) ↦ (RER, RBR), R = diag(−1, 1, 1) (S5). So χ^{(s−iπ)} is the wave function of e_{F′} with F′ = Λ_{−s*}(j_*F), j_*F(y) := e^{iπMᵀ}F(jy), a real test field supported in jW_R = W_L. By locality f(s − iπ) = ⟨e_G, e_{F′}⟩ is real. Lemma 4 with K₀ = span{e_F} gives V(s) = Δ^{is/2π}. ∎

By Lemma 2 the same holds for every wedge gW_R, g a Poincaré transformation, with X replaced by g_*X.

**Step 6. Geometry (R = 1).** Let κ_τ(x) = (x − τb x²)/σ_τ(x), σ_τ = 1 − 2τ b·x + τ²b²x², b = e₁: the special conformal flow of the Killing field K₁ = 2(b·x)x − x²b. S2a (exact):
- σ_τ(x) = −τ²(x − p_τ)² with p_τ = −e₁/τ on the slice t = 0. Every point of O = O₁ is spacelike to p_τ for 0 < τ ≤ 1 (|x⃗ − p⃗_τ| ≥ 1/τ − r ≥ 1 − r > |t|), so σ_τ > 0 on O and κ_τ is smooth on a neighbourhood of every compact subset of O for all τ ∈ [0, 1].
- (κx − κy)² = (x − y)²/(σ(x)σ(y)); κ := κ₁ preserves t = 0 and maps the unit sphere (minus p₁) onto the plane x¹ = ½.
- κ maps O onto the wedge 𝒲 = {x¹ < ½ − |t|}: O is connected, contains 0 = κ(0), and its points are spacelike to the sphere, so κ(O) lies in the component of the set spacelike to the plane that contains 0; conversely κ₋₁ is smooth on 𝒲 (its singular cone is that of +e₁/τ, spacelike to 𝒲), and its image is open, connected, contains 0 and is spacelike to the sphere, so it lies in O.
- Dκ·ζ₁ = X_𝒲∘κ with X_𝒲 = −(x¹ − ½)∂_t − t∂₁, the future-directed boost field of 𝒲. So the flow Λ_s of ζ₁ on O is κ⁻¹∘(boost of 𝒲)∘κ; in particular it is complete on O.

**Step 7. The theorem for R = 1.** On the complex span D_O of {e_F : F ∈ C_c^∞(O)} define U_κ Σc_j e_{F_j} := Σc_j e_{κ_*F_j}. By Lemma 2(i), applied on a neighbourhood of the supports with τ₁ = 1, it preserves inner products, so it is well defined and isometric. D_O is dense (Lemma 3); κ_* is a bijection of C_c^∞(O) onto C_c^∞(𝒲) (inverse κ₋₁*), so the range is the complex span of {e_G : G ∈ C_c^∞(𝒲)}, dense by Lemma 3. Hence U_κ is a unitary of H with U_κK(O) = K(𝒲). No conformal transformation is implemented outside O, and nothing is said about the cover of compactified Minkowski space: the density of D_O does the work a global representation would do.

So K(O) is standard and Δ_{K(O)}^{it} = U_κ*Δ_{K(𝒲)}^{it}U_κ. By Lemma 5, Δ_{K(𝒲)}^{it}e_G = e_{Λ^𝒲_{−2πt*}G}. By Lemma 2(ii) with Y = ζ₁ and Step 6, κ_*δ†_ζ = δ†_{X_𝒲}κ_*, so s ↦ κ_*Λ_{s*}F and s ↦ Λ^𝒲_{s*}κ_*F solve the same transport equation with the same initial value and are equal. Therefore

  Δ_{K(O)}^{it}e_F = U_κ* e_{Λ^𝒲_{−2πt*}κ_*F} = U_κ* e_{κ_*Λ_{−2πt*}F} = e_{Λ_{−2πt*}F}. This is (a).

**Step 8. The rest.**
- *(b).* Δ^{it} is unitary, so σ(Δ^{it}u, e_F) = σ(u, Δ^{−it}e_F) = σ(u, e_{Λ_{2πt*}F}). Differentiating, d/ds C(Λ_{s*}F) = C(δ_ζ†F_s) = (δ_ζC)(F_s). At the slice t = 0 the field ζ_R has λ = 0, A = 0, ζ⁰ = N_R and β = −x/R, which gives the displayed generator.
- *Radius R and other centres.* Dilations and translations are global flows; by Lemma 2 they give unitaries mapping K(O₁) onto K(O_R + a) and carrying ζ₁ to (D_R)_*ζ₁ = ζ_R.
- *(c)* is wave69's C1. *(d)*: K standard ⇔ Ω cyclic and separating for R(K), and Δ_{R(K)}^{it} = Γ(Δ_K^{it}).

### Premise ledger

| Premise | Status | Locator |
|---|---|---|
| Step 0 formulas for χ_{e_F} | derived here from wave69 Step 1; the signs are cross-checked by the field equations on 𝒲 | `scripts/common.py`, S1a |
| Lemma 1 rational identity, Lemma 1′ | proved here, exact | `scripts/s1_invariance.py`, `out/s1_out.json` |
| Lemmas 2, 3, 4, 5 and Steps 6–8 | proved here | this file; S2a, S5 exact |
| Boundary values: e^{−εω} regularisation gives 𝒲 = b.v. of 𝒲_rat; the two-point function of the massless scalar is −1/(4π²(ξ² − iεξ⁰)) | textbook, not sourced | — |
| Three-lines theorem, Schwarz reflection, uniqueness of polar decomposition, invariant-domain core argument | textbook, not sourced | — |
| wave69 Steps 1, 2, 7 (vectors of E(f), B(g); L(ℝ³)′ = {0}; time-zero vectors are limits of e_F) and C1 | proved in wave69; not re-audited here | wave69/duality/REPORT.md §2 |
| Second quantisation: Δ_{R(K)} = Γ(Δ_K), standard ⇔ cyclic and separating | statement read; proof not read | Longo–Morsella 2012.00565 pdf p2 ("the vacuum modular group … is the second quantisation … of the modular group of the local standard subspace"); Brunetti–Guido–Longo math-ph/0203021 pdf p9 (standard ⇔ Reeh–Schlieder, citing Araki, Eckmann–Osterwalder, Leyland–Roberts–Testard) |
| The TT Fock space is the vacuum representation of the curvature field | inherited from wave69 (its weakest step 2) | — |

## 3. The other routes, and where the superworld went

- **Hislop–Longo 1982 (read, pdf p2–7).** Same architecture for the scalar: a unitary U_ρ for the ray inversion, defined on a dense domain of test-function vectors (their Thm 1), mapping double cone, wedge and forward cone onto each other; the wedge from Bisognano–Wichmann and the cone from Buchholz, both cited. Step 7 above is their Thm 1 in a weaker form (one map, one region), which is all the modular group needs; Lemma 5 replaces the citation.
- **Lledó + Brunetti–Guido–Longo 1993.** Lledó proves SU(2,2) covariance of his helicity nets (Thm 5.7, pdf p24) and lists "the modular groups of M(K₁), M(W_r), M(V₊) act geometrically" (item (ii), pdf p26) as a consequence of BGL, without proof. That route is sound in outline but needs his net identified with the curvature net; I did not do that and did not need it.
- **Where an obstruction could have sat.** The modular *group* of a double cone only needs a unitary from K(O) to K(wedge); the map κ is regular on O, and Reeh–Schlieder makes the locally defined isometry a unitary. The cover of compactified Minkowski space matters for the modular *conjugation*: J_O = U_κ*J_𝒲U_κ sends K(O) to U_κ*K(𝒲′), and κ⁻¹(𝒲′) is the spacelike complement of O together with regions timelike to O. That J_O K(O) = K(O′) is wave69's duality theorem (equivalently Huygens' principle, Step 0); it is not used here.
- **Integer helicity** enters in one place: e^{iπM} must be real (Lemma 5), which holds because the boost weights are integers (S5).
- **Paper A's own derivation (Corollary mod).** Transfer each (l, parity) sector unitarily onto the scalar l-waves (Thm transfer (a), (b), resting on Prop. Nl), then Hislop–Longo for the scalar. That is a second, independent proof of the same statement; the two agree, including the generator 2πδ_ζ (after B_A = −B and the conjugate complex structure). S3 below re-derives its constant c_l from the tensor two-point function.

## 4. Numerical and computer-algebra checks

### 4.1 S1: two-point function and conformal invariance (`scripts/s1_invariance.py`): PASS after one correction

- **a.** All nine structure checks hold (`a_structure` in `out/s1_out.json`).
- **b. As registered: FAIL, and the cause is mine.** I registered weight w = 3 and boost sign s = −1, the latter transcribed from Paper A's eq. (flow). With s = −1 the three boosts and four special conformal fields give non-zero numerators (276 to 304 of 324 components), while the registered "plant" s = +1 passes on every generator tested (`out/s1_out_as_registered.json`, `out/s1_as_registered.log`). Paper A's conventions are ∂_tE = −curl B, ∂_tB = curl E; wave69's are Ė = curl B, Ḃ = −curl E (confirmed on 𝒲 by S1a). So B_A = −B and the sign of the boost term flips. Nothing physical differs.
- **b. Corrected run (post hoc, `out/s1_out.json`):** w = 3, s = +1: zero numerator in all 324 components for each of D, N₁, N₂, N₃, J₁₂, J₁₃, J₂₃, K₀, K₁, K₂, K₃.
- **c.** [δ_X, δ_Y] = δ_{[X,Y]}: 105 of 105 pairs.
- **Plants caught.** w = 2 and w = 4: 294 of 324 components non-zero for each of D, K₀, K₁. Opposite boost sign: 276, 294, 304 non-zero for N₁, K₀, K₁. Exactly one (w, s) makes the two-point function invariant.

### 4.2 S2: geometry and finite covariance (`scripts/s2_geometry.py`): PASS after one integrator fix

- **a.** All exact identities of Step 6 hold; σ_τ > 0 on 20,000 sampled (x, τ) (minimum 0.019); 5,000 sampled points of O land in the wedge.
- A registration slip: the PREREG's route says b = −e₁ while the registered plane x¹ = ½ and X_𝒲 belong to b = +e₁ in the convention used. I used b = +e₁.
- **b. As registered: MISS.** With my RK4 cocycle the relative error was below 1e−12 on 11 of 20 pairs and as large as 7.5e−8 on the worst (`out/s2_out_as_registered.json`). The error tracked the integrator's own extrapolation error (up to 2e−5 for points near the singular cone), not the identity.
- **b. Fix (post hoc, `out/s2_out_fix.json`):** Taylor-series integration of the cocycle (order 70, steps 0.3 of the distance to the complex zeros of σ_τ). Worst relative error 2.2e−29 over ten spacelike and ten timelike pairs.
- **Plant caught.** Weight 4 in the cocycle: relative error between 0.063 and 74 on every pair (threshold 1e−3).

### 4.3 S3: the radial sector of the tensor two-point function (`scripts/s3_sector.py`): PASS

ψ = x·E·x. Exact closed form of G = x_i x_j y_k y_l 𝒲^{EE}_{ij,kl} in `out/s3_out.json`. Legendre projections at 6 radius pairs and t ∈ {0, 0.1}, 40 digits:
- g₀ = g₁ = 0 (below 1e−41 of the l = 2 scale): the curvature has no l ≤ 1 part.
- g_l/d_l = 6, 30, 90 for l = 2, 3, 4, with zero spread over the 12 points. That is (l+2)!/(4(l−2)!) = 1/c_l with Paper A's c_l = 4(l−2)!/(l+2)!: an independent confirmation of Prop. Nl's constant for l = 2, 3, 4, from the two-point function instead of Gaussian integrals.
- ψ transforms under ζ as a conformal scalar of dimension 1 for any symmetric E, B (exact; Paper A's Lemma stepA is the t-derivative of this).

So in each (l, parity) sector the chain field u = rψ has 1/c_l times the scalar l-wave vacuum correlators, and S4's chain is the helicity-2 two-point function in that sector.

### 4.4 S4: lattice modular Hamiltonian of the ball, l = 2 and 3 (`scripts/s4_lattice.py`): PASS on the thresholds, one clause missed

Chain, correlators, modular Hamiltonian and statistics exactly as registered. T1 tests the momentum block against 2πN_R, T2 the field block against 2π[−∂N_R∂ + N_R l(l+1)/r²]; three Gaussian profiles centred at 0.3R, 0.5R, 0.7R. Entries are (T1, T2); 1 is agreement.

| l | n | digits | largest ε_k | c = 0.3 | c = 0.5 | c = 0.7 |
|---|---|---|---|---|---|---|
| 2 | 10 | 80 | 83 | 0.9910, 0.9868 | 0.9939, 0.9959 | 0.9954, 0.9973 |
| 2 | 20 | 200 | 165 | 0.9992, 0.9975 | 0.9987, 0.9994 | 0.9989, 0.9995 |
| 2 | 40 | 400 | 330 | 1.0009, 1.0001 | 0.9999, 1.0002 | 0.9998, 1.0000 |
| 3 | 10 | 80 | 93 | 0.9718, 0.9505 | 0.9850, 0.9859 | 0.9893, 0.9909 |
| 3 | 20 | 200 | 184 | 0.9948, 0.9841 | 0.9962, 0.9965 | 0.9972, 0.9977 |
| 3 | 40 | 400 | 366 | 0.9993, 0.9912 | 0.9991, 0.9992 | 0.9993, 0.9994 |
N80ROWS

- **Thresholds.** At the largest n every |T − 1| is at most 0.0088 against the registered 0.05. M and N are stable between the two precisions to relative 1e−133 or better (registered 1e−8).
- **The clause "errors do not grow from n = 20 to the largest n" is missed in a few cells at the 1e−4 level:** MONOTONE. The errors there pass through zero and settle at a few parts in 10⁴, the size of the finite-volume effect of the wall at 6R. I registered the clause without allowing for a sign change; as written it is a miss.
- **Plants caught** (n = 40; values for the three profiles):
 - P1, weight (R² + r²)/2R: T1 = 0.82, 0.59, 0.34 and T2 = 0.87, 0.60, 0.33 (l = 2); the same to two digits for l = 3. Two of three profiles miss 1 by more than 0.25, as registered; the innermost profile, where the two weights differ least, misses by 0.18.
 - P2, π instead of 2π: 2.00 everywhere.
 - P3, l(l+1) → 0 in the T2 target: T2 = 4.04, 1.81, 1.40 (l = 2) and 7.02, 2.62, 1.81 (l = 3).
- **Resolution.** With N_R → (R² − (1+ε)r²)/2R the outer profile reads 1.020 at ε = 0.02 and 1.052 at ε = 0.05: the test resolves the coefficient of r² to about 5% at the registered threshold, and to about 0.5% if the observed scatter (≤ 0.2% for the two outer profiles) is taken as the yardstick.
- **What S4 is evidence for.** The weight, the temperature 2π, the sector, and the sign (the lattice modular Hamiltonian is +2π times a positive energy, which fixes the direction in (d)). It is the scalar l-wave lattice; S3 is what ties it to helicity 2. It does not test Lemma 4 or Step 7.

### 4.5 S5: the boost generator (`scripts/s5_wedge.py`, added after the PREREG)

Supports Lemma 5 only: M is diagonalisable with eigenvalues 0 (×4), ±1 (×2), ±2 (×2) on symmetric (E, B); e^{iπM} is real and equals (E, B) ↦ (RER, RBR).

## 5. Consequence for Paper A's Theorem 1 (clause M)

- **What Theorem 1 actually uses.** In reviseA18 the proof of Theorem transfer uses (i) that rotations and the inversion commute with Δ_{H(B)} (general), (ii) the *scalar* modular group and entropy (Hislop–Longo 1982, Longo–Morsella), and (iii) the transfer V_σ. The helicity-2 modular group is then *derived* (Corollary mod). Hislop 1988 is cited for it in the introduction ("by conformal covariance so does that of a diamond") and in the remark after Corollary mod, and for duality in Prop. alg. So M does not logically rest on Hislop's modular result; the brief's premise is stricter than the paper needs.
- **What this lane adds.** A direct proof of the statement of Corollary mod, with no transfer and no Hislop: two independent proofs now give the same group and the same generator. Together with wave69, both of Paper A's citations to Hislop 1988 (modular group, duality) can be backed by proofs on file.
- **S3 backs Prop. Nl.** c_l = 4(l−2)!/(l+2)! is reproduced for l = 2, 3, 4 from the two-point function.
- **Suggested wording.** In the introduction: "… and by conformal covariance so does that of a diamond [HL, Hislop, BGL]; for helicity 2 a direct proof is given in Appendix Y." Appendix Y is §2 of this file. After Corollary mod: "A proof that does not pass through the scalar is in Appendix Y."
- **Convention warning for whoever merges appendices.** wave69 and this file use Ė = curl B; Paper A uses ∂_tE = −curl B (B_A = −B). wave69 calls its complex structure J; the positive-energy one is −J. Paper A writes the modular group on waves with the conjugate convention (Δ^{it} = U(Λ_B(−2πt)), A = 2πδ_ζ); its entropy formula is fixed by S ≥ 0 and is not affected. I did not audit Paper A's signs beyond this.
- **Not touched.** Theorem 1(c) (the entropy for data not localised in the ball, the cutting projection, the Born-term cancellation) is not re-proved here.

## 6. Weakest step

1. **Lemma 4 and Lemma 5 are mine and unchecked.** They are the two places where a sign or a domain could be wrong: the core argument in Lemma 4 item 5, and the claim in Lemma 5 that χ^{(z)} is continuous in H up to both edges of the strip (it rests on the uniform integration-by-parts bound). The lattice confirms the outcome (weight, 2π, sign), not these arguments.
2. **Step 0's boundary-value statement** and the scalar two-point function are quoted from textbooks without a locator.
3. **The link to the algebra** is the imported second-quantisation theorem.
4. **Inherited:** wave69's Steps 1, 2, 7 and its model identification.
5. **S4 is a radial chain.** It would not see an error that S3 and the chain share; it is tied to the tensor field only through S3's projection, which is exact for l ≤ 4 and not run beyond.

## 7. Recommended next lane

**w71-entropy-direct (Fable): Paper A's Theorem 1(c) for ball-localised data, straight from the generator.** With A = 2πδ_ζ now proved without the transfer, the coherent-state entropy of data supported in the ball is S = 2π·½σ(δ_ζu, u) (up to the convention sign), a quadratic form in (h, π) that can be evaluated in closed form: it should be the ζ-canonical energy with no boundary term, and sector by sector Σ c_l Q[u]. That is an independent proof of the localised half of clause M that never touches V_σ or Longo–Morsella, and it isolates what remains specific to the transfer: the cutting projection for data that do not vanish on the sphere. Deliverable: the identity S = Σ c_l Q[u] for u ∈ K(O) derived from δ_ζ, an exact check on polynomial-times-bump TT data for l = 2, 3, and a planted boundary term that must be caught.

Cheaper companion, not a physics lane: a blind check of §2 by a fresh model given only the theorem, Steps 0–8 and the ledger, asked to break Lemma 4 item 5, Lemma 5's continuity at the strip edges, and Step 6's claim κ(O) = 𝒲.

## Files
- `PREREG.md`, `PREREG.sha256`
- `scripts/common.py` (two-point function, conformal action)
- `scripts/s1_invariance.py`; `out/s1_out.json` (corrected sign), `out/s1_out_as_registered.json`, `out/s1_as_registered.log`
- `scripts/s2_geometry.py`; `out/s2_out_fix.json`, `out/s2_fix_run.log` (Taylor integrator), `out/s2_out_as_registered.json`, `out/s2_as_registered.log` (RK4)
- `scripts/s3_sector.py`; `out/s3_out.json`, `out/s3_run.log`
- `scripts/s4_lattice.py`; `out/s4_l{2,3}_n{10,20,40}.json`, `out/s4_run_*.log`N80FILES
- `scripts/s5_wedge.py`; `out/s5_out.json`
- S2–S4 were run with a scratch virtualenv (mpmath with gmpy2, sympy); S1 and S5 with the system python.
- `REPORT.md` (this file)
