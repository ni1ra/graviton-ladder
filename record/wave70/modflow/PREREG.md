# PREREG: w70-modflow. The modular group of the ball in the free helicity-2 curvature net

Lane w70-modflow (Fable). Written before any script in this folder exists or runs.

Read before writing: wave69/duality/REPORT.md (all); wave36/reviseA18/NOTE.tex lines 44–165 (Setup, Prop. Nl, Lemma stepA, Theorem transfer, Corollary mod); Hislop–Longo 1982 §1–3 (pdf p1–7); Longo–Morsella 2012.00565 §4.1 (Thm 4.1); Brunetti–Guido–Longo math-ph/0203021 §3 and §4 remark (2); Lledó math-ph/0006018 Thm 5.7, Prop. 5.8 and the list after it. Two web searches (arXiv) for modular theory of free fields of helicity ≥ 1 on a double cone: nothing newer than La Piana–Morsella 2405.07888 (helicity ½, quoting Hislop 1988 second hand). Hislop 1988 and Bisognano–Wichmann 1975/76 are not readable here.

## 1. Claim under test

Setting as in wave69/duality/REPORT.md §1: H = TT Cauchy data u = (h, π), σ(u, u′) = ∫(h·π′ − h′·π), E[u] = −½Δh, B[u] = ½ curl π; for a real test field F = (f, g) on spacetime with values in S²₀ ⊕ S²₀, e_F ∈ H is the vector with σ(u, e_F) = ∫dt∫d³x (f·E + g·B)[T_t u]; K(U) = closed real span of {e_F : supp F ⊂ U}. The Hilbert-space complex structure is the one in which the time translation of observables has positive energy (I record below that this is −J in wave69's notation; wave69's theorem does not depend on the sign).

**Theorem (to prove or refute).** Let O = D(B_R) be the double cone over the ball of radius R centred at the origin of the t = 0 slice, and

  ζ_R = (1/2R)[(R² − t² − r²)∂_t − 2t x·∇]  (equal to ((R² − r²)/2R)∂_t at t = 0).

Then K(O) is standard and its modular group is the conformal flow of ζ_R:

  Δ_{K(O)}^{it} e_F = e_{F_{−2πt}} for every real F ∈ C_c^∞(O),

where F_s is F pushed forward along the flow Λ_s of ζ_R by the weighted action dual to δ_X(E, B) = X·∂(E, B) + (local Lorentz part of ∂X acting on the tensor indices) + 3λ(E, B), λ = ∂·X/4. Equivalently: for every u ∈ H the curvature of Δ^{it}u inside O is the curvature of u inside O transported by Λ_{2πt}; on smooth data in B, d/dt Δ^{it}u|₀ = 2π δ_ζ u with, at t = 0, δ_ζE = N curl B + 2 sym(x × B)/R, δ_ζB = −N curl E − 2 sym(x × E)/R, N = (R² − r²)/2R (Paper A eq. (flow) at R = 1). By second quantisation the vacuum modular group of A(O) is Γ(Δ^{it}).

**Planned route.**
1. Explicit one-particle wave functions of e_F and the two-point distribution W(F, G) = ⟨e_F, e_G⟩ = boundary value of Q(i∂) applied to −1/(4π²ξ²), Q a degree-4 polynomial matrix.
2. Lemma 1 (infinitesimal conformal invariance): (δ_X ⊗ 1 + 1 ⊗ δ_X)W = 0 for all 15 conformal Killing fields, as a rational identity (computer algebra) plus a boundary-value argument.
3. Lemma 2 (finite covariance along flows): W(κ_*F, κ_*G) = W(F, G) whenever the flow of X is non-singular near the supports. Needs [δ_X, δ_Y] = δ_{[X,Y]}.
4. Lemma 3 (one-particle Reeh–Schlieder) and Lemma 4 (a KMS criterion for standard subspaces, proved here).
5. Lemma 5 (one-particle Bisognano–Wichmann for the wedge of this net), by Lemma 4 with the explicit wave functions.
6. Geometry: the special conformal transformation κ with b = −e₁/R maps O onto a wedge and is the time-1 map of a flow that is non-singular on compacts of O; it carries ζ_R to the wedge's boost field. Then U_κ e_F := e_{κ_*F} is unitary from the complex span of K(O) onto that of K(W), and the theorem follows. No global representation of the conformal group, and no statement about the superworld, is used.

The alternative (Paper A's own Corollary mod: transfer to scalar l-waves, then Hislop–Longo) is compared, not re-proved.

## 2. Bins

- **PROOF.** Every step is proved here or read in a printed argument with a locator; S1–S3 pass on every registered identity; S4 passes its registered thresholds; every plant is caught.
- **PROOF-WITH-CONDITION.** The proof stands given named premises that are neither proved here nor read in print. The condition is stated.
- **REFUTED-WITH-MECHANISM.** A step fails and the failure exhibits how the modular group differs from the flow (for example a weight for which the two-point function is not invariant, a non-geometric piece from the superworld, or a lattice modular Hamiltonian that converges to something else).
- **OPEN.** None of the above.

Expectation before computing: PROOF, or PROOF-WITH-CONDITION on the imported second-quantisation theorem (Araki; Eckmann–Osterwalder) and on wave69's setting.

## 3. Computations and plants

### S1: two-point function and conformal invariance (`scripts/s1_invariance.py`, sympy, exact)
- a. W^{EE}, W^{EB}, W^{BE}, W^{BB} from the TT mode sum as Q(i∂)(1/ξ²); check W^{BB} = W^{EE}, W^{BE} = −W^{EB}, symmetry, tracelessness and transversality in each index pair, and the field equations Ė = curl B, Ḃ = −curl E on the first slot.
- b. Lemma 1: the numerator of (δ_X ⊗ 1 + 1 ⊗ δ_X)W vanishes identically for dilation, the three boosts, the three rotations and the four special conformal fields, with weight w = 3 and boost sign s = −1 (δE ∋ s(ε_ikl β_k B_lj + i↔j), β = ∇X⁰).
- c. [δ_X, δ_Y] = δ_{[X,Y]} for all 105 pairs of the 15 generators.
- **Plants S1-P1, P2, P3:** w = 2, w = 4, and s = +1. Each must give a non-zero numerator for at least one generator. If any plant passes, S1 has no power and is VOID.

### S2: geometry and finite covariance (`scripts/s2_geometry.py`)
- a. Exact: κ = SCT with b = −e₁ (R = 1) has (κx − κy)² = (x − y)²/(σ(x)σ(y)); it maps the unit sphere at t = 0 to the plane x₁ = ½; Dκ·ζ = X_W∘κ with X_W = −(x₁ − ½)∂_t − t∂₁; σ_τ > 0 on O for all τ ∈ [0, 1] (shown analytically in the report, sampled here).
- b. Numeric (mpmath, 30 digits): integrate the flow and the 10 × 10 cocycle A_τ(x) for 20 random pairs x, y ∈ O, ten spacelike and ten timelike separated, and compare A(x) W(x − y) A(y)ᵀ with W(κx − κy) (in the convention the script fixes). Pass: relative error < 1e−12 on every pair.
- **Plant S2-P:** the same with weight 4 in the cocycle. Must fail by more than 1e−3.

### S3: the radial sector of the helicity-2 two-point function (`scripts/s3_sector.py`)
ψ = x·E·x. From S1's W^{EE}: G(r, r′, cos γ; t) = x_i x_j y_k y_l W^{EE}_{ij,kl}(t, x − y). Its Legendre projections g_l = ∫P_l(c)G dc are compared with the scalar's d_l = ∫P_l(c)D⁺ dc at 6 pairs (r, r′) and t ∈ {0, 0.1} (spacelike).
- Registered: g_0 = g_1 = 0; for l = 2, 3, 4 the ratio g_l/d_l is one constant per l, independent of (r, r′, t), to 1e−10. (Paper A's Prop. Nl predicts the ratio ∝ (l+2)!/(l−2)!; I register the proportionality to that with one overall constant.)
- This licenses S4: in each (l, parity) sector the chain field u = rψ has the scalar l-wave vacuum correlators, so its lattice is the helicity-2 two-point function in that sector.
- Also here (exact): Paper A's Lemma stepA from S1's δ_ζ: x·δE·x = Nψ̇ and x·δĖ·x = NΔψ − x·∇ψ − ψ on TT solutions.

### S4: lattice modular Hamiltonian of the ball, radial sectors l = 2 and 3 (`scripts/s4_lattice.py`, mpmath high precision)
- Chain: sites j = 1..n_tot, a = 1, K_jj = 2 + l(l+1)/j², K_{j,j±1} = −1, Dirichlet walls; vacuum X = ½K^{−1/2}, P = ½K^{1/2}; ball = sites 1..n, R = n + ½; n_tot = 6n; n ∈ {10, 20, 40} (80 if time allows).
- Modular Hamiltonian ½(πᵀMπ + uᵀNu) of the Gaussian state from X_A, P_A alone (Peschel; Casini–Huerta 2009), with precision raised until M, N are stable to 1e−8 between two precisions.
- Target (conformal generator): M ↔ 2π N_R(r) δ(r − r′), N ↔ 2π[−∂_r N_R ∂_r + N_R l(l+1)/r²], N_R = (R² − r²)/2R.
- Statistics, three smooth profiles g_a(r) = exp(−(r − c_a R)²/(2(0.12R)²)), c = 0.3, 0.5, 0.7:
 T1_a = gᵀMg / [2π Σ_j N_R(j) g_j²], T2_a = gᵀNg / [2π Σ_j N_R(j+½)(g_{j+1} − g_j)² + 2π Σ_j N_R(j) l(l+1) g_j²/j²].
- **Pass:** at the largest n, |T1_a − 1| ≤ 0.05 and |T2_a − 1| ≤ 0.05 for all a and both l, and the errors do not grow from n = 20 to the largest n.
- **Plants (must be caught, i.e. the same statistic computed against the planted target misses 1 by more than 0.25 for at least two of the three profiles):**
 - S4-P1 wrong weight: N_R → (R² + r²)/2R.
 - S4-P2 wrong temperature: 2π → π.
 - S4-P3 wrong sector in T2: target built with l(l+1) → 0.
- Resolution probe (reported, not a pass condition): the smallest ε for which the target with N_R → (R² − (1+ε)r²)/2R is rejected at the 0.05 threshold.
- If S4 fails with the plants caught, that is evidence against the theorem or against the lattice statistic; one fix round is allowed and will be labelled post hoc.

## 4. What would change my mind

- S1b failing for every weight: the two-point function is not conformally invariant in the form used, and route 6 is dead.
- S2a failing: κ does not carry ζ_R to a boost; the flow is misidentified.
- S4 converging to a different weight with S1–S3 passing: an error in the proof (most likely in the sign or normalisation of Lemma 5), to be found before any bin is given.
