# PRINT_READ — the five imported statements, read in print

Date: 10 Nov 2026. Sub-lane of w69-final. Read-only pass over `R/sources`; page numbers are PDF pages.
Method: `pdftotext -layout` for every file; page images viewed where the text layer drops formulas
(Araki 1976 pp1–2, Araki 1977 pp5–9, Hislop–Longo pp5–6). One numerical check was run (I4).

Status vocabulary
- READ-COMPLETE: statement and proof read line by line; every step is in the file or is elementary.
- READ-COMPLETE-MODULO X: proof read line by line; it rests on X, which the file cites and does not prove.
- READ-STATEMENT-ONLY: statement located; no proof in the file (citation only), or proof not read.
- DEFERS-TO X: the file's proof is a pointer to X.

Summary table

| Import | Verdict |
|---|---|
| I1 CLR Thm 2.2, Cor 2.3, Prop 2.4, Prop 2.5, Lemma 4.3, Thm 4.5 | READ-COMPLETE-MODULO (a) Araki's formula as the definition, (b) lower semicontinuity of relative entropy, (c) `Δ_R(H) = Γ(Δ_H)` and `R(H′) ⊂ R(H)′`, (d) tensor additivity (replaceable by I3) |
| I1 second-quantisation facts on CLR p19 | READ-STATEMENT-ONLY ("we recall", no proof, no specific citation) |
| I2 `R(K)′ = R(K′)` | Garbarz–Palau: READ, proof is for the pair form `R_W(K1,K2)` only, DEFERS-TO Eckmann–Osterwalder Lemma 4 for one step, one acknowledged soft spot. Camassa, BGL: READ-STATEMENT-ONLY |
| I2 lattice facts ∨, ∧; "Ω cyclic and separating iff K standard"; `Δ_R(K) = Γ(Δ_K)` | READ-STATEMENT-ONLY everywhere (cited to Araki 1963, Eckmann–Osterwalder 1973, Leyland–Roberts–Testard 1978). No file read contains a proof |
| I3 Reible Thm IV.3.7 | READ-COMPLETE-MODULO Löwner–Heinz, bounded Hansen–Pedersen (both cited to textbooks) and the spatial-derivative facts of Reible §III.4 (partly cited to Hiai), which I skimmed and did not read line by line |
| I3 Araki 1976 (1.6), Araki 1977 Thm 3.8(2) | READ-STATEMENT-ONLY; restriction monotonicity only for three special kinds of subalgebra |
| I4 Longo–Morsella Thm 4.1 | DEFERS-TO Hislop–Longo (sketch) |
| I4 Thm 4.7, Thm 4.12, Thm 5.4, Cor 5.5 | READ; the algebra was redone here. Thm 5.4 / Cor 5.5 coefficient `πD/R` is right; Thm 4.12's `2π(d−1)/R` is wrong by a factor 4. Printed proofs of Lemma 4.6 and Thm 5.4 are sketches with gaps, closed below by a direct argument |
| I5 Hislop–Longo Thm 2, Cor 2 | READ-COMPLETE-MODULO Buchholz 1978 (forward-cone modular group; cited, conference proceedings, odd d) and Araki 1964 (time-zero algebras). Bisognano–Wichmann is cited, not proved, and is not used in the proof for odd d |

---

## I1. Ciolli–Longo–Ruzzi, "The information in a wave", arXiv 1906.01707

File: `sources/ciolli-longo-ruzzi-2019-information-in-a-wave-1906.01707.pdf` (24 pages, all read except §3.1 pp12–14, the wedge computation, which is outside the scope).

### Conventions (as printed)

- Inner product `(·,·)` is antilinear in the FIRST slot, linear in the second. Not stated as a sentence; it is forced by eq (11), p10: `ℑ(h, Ah) = −(h, log Δ h)` with `A = −i log Δ`, and by eq (31), p20: `−(h, log Δ_H h) = i(h, i log Δ_H h) = −ℑ(h, i log Δ_H h)`.
- Symplectic complement (p6): `H′ = (iH)^{⊥_R}`, orthogonal with respect to `ℜ(·,·)`; equivalently `H′ = {k : ℑ(k,h) = 0 ∀h ∈ H}`.
- Standard: closed, `H + iH` dense, `H ∩ iH = {0}`. Factorial: `H ∩ H′ = {0}`, equivalently `H + H′` dense, equivalently 1 is not an eigenvalue of `Δ_H` (cited to Figliolini–Guido). `S = J Δ^{1/2}`, `Δ^{is}H = H`, `JH = H′`.
- Weyl unitaries (p18): `V(h) e^k/‖e^k‖ = e^{h+k}/‖e^{h+k}‖`; eq (25) `V(h+k) = e^{iℑ(h,k)} V(h)V(k)`; eq (26) `V(h)ξ = e^{−(h,h)/2} e^h`; eq (27) `ϕ(V(h)) = e^{−‖h‖²/2}`.
- `R(H) = {V(h) : h ∈ H}″` (p19).
- Relative entropy, eq (23), p17: `S(ϕ1‖ϕ2) = −(ξ1, log Δ_{ξ2,ξ1} ξ1)` for normal faithful states, understood through the spectral integral (24). First argument is "the state", second is the reference: `S(ϕ1‖ϕ2) = tr ρ1(log ρ1 − log ρ2)` in finite dimension.
- Coherent state: `ϕ_k = (V(k)ξ, · V(k)ξ)` restricted to `R(H)`; `ϕ` is the vacuum.

### Printed statements

- Cutting projection (p6): for `H` standard factorial and `k ∈ H + H′`, `P_H k = h` where `k = h + h′` is the unique decomposition. Domain `H + H′`, real linear, densely defined.
- Lemma 2.1 (pp6–7): `P_H` closed; `P_H² = P_H`; `P_H + P_{H′} = 1` on `H + H′`; `P_H^* = P_{iH} = −i P_H i` (real adjoint); `Δ^{is} P_H = P_H Δ^{is}`.
- Thm 2.2 (p7): with `a(λ) = λ^{−1/2}(λ^{−1/2} − λ^{1/2})^{−1} = 1/(1−λ)`, `b(λ) = (λ^{−1/2} − λ^{1/2})^{−1} = λ^{1/2}/(1−λ)`, `D = Dom(a(Δ)) ∩ Dom(b(Δ))`: "We have `D ⊂ H + H′` and `P_H|_D = Δ^{−1/2}(Δ^{−1/2} − Δ^{1/2})^{−1} + J(Δ^{−1/2} − Δ^{1/2})^{−1}`. Moreover, `D` is a core for `P_H` … Indeed, already `D0` is a core", `D0 = ∪_ε E_ε𝓗`, `E_ε` the spectral projection of `Δ` for `(ε, 1−ε) ∪ ((1−ε)^{−1}, ε^{−1})`.
- Cor 2.3 (p8): for continuous `f` with a zero of order ≥ 1 at `λ = 1` and growth of order ≤ 1/2 at `∞`, `Ran f(Δ) ⊂ H + H′`. "In particular, for any `k ∈ 𝓗`, we have `k − Δ^{is}k ∈ H + H′`, and `log Δ k ∈ H + H′`, `k ∈ Dom(log Δ)`."
- Definition (10), p9: for ANY `k ∈ 𝓗`, `S_k = ℑ(k, P_H A k)`, `A = −i log Δ_H`, where the right side is DEFINED by
  `−(k, P_H A k) ≡ i∫ a(λ) log λ d(k, E(λ)k) − i∫ b(λ) log λ d(k, J E(λ)k)`.
  Unwound: `S_k = ∫ (−log λ)/(1−λ) d(k,E(λ)k) + ℜ∫ λ^{1/2} log λ/(1−λ) d(k, J E(λ)k)`. The first integrand is ≥ 0, the second is bounded, so `S_k ∈ (−∞, +∞]` is always defined.
- Prop 2.4 (p9): `S_k < +∞` iff `−∫_0^1 log λ d(k,E(λ)k) < +∞`, i.e. iff `k ∈ Dom(√|log Δ| E_−)`, `E_−` the negative spectral projection of `log Δ`. All vectors in `Dom(log Δ)` have finite entropy.
- Prop 2.5 (p10): (a) if `k = h + h′`, then `S_k = −(h, log Δ h)`; (b) `S_k = ℜ(k, i P_H i log Δ k) = −ℜ(k, P_H^* log Δ k)`; (c) `S_k ≥ 0` and `S_k = 0` iff `k ∈ H′`; (d) additivity under direct sums; (e) `S_k = lim_{ε→0+} S_{k_ε}`, non-decreasing, `k_ε = E_ε k`; (f) `S_k = −lim_{s→0} ℑ(k, P_H(Δ^{is}k − k))/s`; (g) continuity in the graph norm of `√|log Δ| E_−`.
- Prop 4.2 (p20): "Let `h ∈ H`. If `h ∈ Dom(log Δ_H)`, the relative entropy on `R(H)` between `ϕ` and `ϕ_h` is given by `S(ϕ_h‖ϕ) = S_h = −(h, log Δ_H h)`."
- Lemma 4.3 (p20): for `k ∈ 𝓗_ε`, `S(ϕ_k‖ϕ) = S_k`.
- **Thm 4.5 (p21): "Let `k ∈ 𝓗`. We have `S(ϕ_k‖ϕ) = S_k = −ℑ(k, P_H i log Δ k)`."** Stated for every `k` in the Hilbert space; both sides may be `+∞`.
- Second-quantisation facts (p19): "we recall that, if `H` a standard subspace, the following hold: `K` dense in `H` ⇒ `R(K) = R(H)`; `ξ` is cyclic and separating for `R(H)`; `Δ^{is}_{R(H)} = Γ(Δ^{is}_H)`, `J_{R(H)} = Γ(J_H)`; `R(H′) = R(H)′`." No proof and no specific citation at that place; the only pointer is the general "We refer to [14, 15]" (Longo's lecture notes) at the start of §2.

### Hypotheses

`H` standard and factorial ("In this paper, we shall always assume the standard subspaces to be factorial", p6; the non-factorial and non-standard cases are reduced by the direct-sum decompositions (3), (4), stated without proof). `𝓗` any complex Hilbert space. No separability assumption is printed.

### Proof status, line by line

- Lemma 2.1, Thm 2.2, Cor 2.3, Prop 2.4: complete in the file. Two one-line gaps, both fine:
  - eq (5) `E_ε H ⊂ H` is attributed to "Lemma 2.1"; the actual reason is that `E_ε` is a real, even function of `log Δ`, hence a strong limit of real combinations of `Δ^{is}`, and `Δ^{is}H = H`.
  - p9 says "`a(λ) log λ` is … positive on (0,1)"; it is negative there. What is meant and used is `−a(λ) log λ ≥ 0`.
- Prop 2.5: complete. The general case of (a) is obtained from (e) by monotone convergence; `−(h, log Δ h)` is then a form value (it may be `+∞`; the part of the integral over `λ > 1` is finite because `h ∈ Dom(Δ^{1/2})`). The proof of (e) cites "Theorem 6", a dead reference (it means Thm 2.2).
- **Independence of `H′`** is printed twice: abstractly in Prop 2.5(a),(c), and at the algebra level in the proof of Lemma 4.3: `S(ϕ·Ad V(k)^*‖ϕ) = S(ϕ·Ad V(h)^*V(h′)^*‖ϕ) = S(ϕ·Ad V(h)^*‖ϕ)`. The second step needs only `V(h′) ∈ R(H)′`, the EASY inclusion `R(H′) ⊂ R(H)′`, which follows from (25). Full duality is not used anywhere in §4.
- Prop 4.2: complete once two things are granted. (i) `Δ_{ξ, V(h)ξ} = Δ_ξ` for `V(h) ∈ R(H)`: elementary (`S_{ξ,uξ} = u S_ξ` for a unitary `u` in the algebra), though the printed line `Δ^{is}_{V(h)ξ,ξ} = V(h)^* Δ^{is}_ξ V(h) Δ^{−is}_ξ` is garbled. (ii) `Δ_{R(H)} = Γ(Δ_H)`: imported, see I2. The differentiation `i d/ds exp(−‖h‖² + (h, Δ^{is}h))|_{s=0} = −(h, log Δ h)` is right. The pointer to Ohya–Petz in Prop 4.1 is not needed in the case used: `h ∈ Dom(log Δ_H)` makes `s ↦ (V(h)ξ, Γ(Δ^{is})V(h)ξ)` twice differentiable, so `V(h)ξ ∈ Dom(log Δ_{R(H)})` and the derivative is the expectation value by the spectral theorem. Lemma 4.3 applies Prop 4.2 only to `h ∈ E_ε H`, where `log Δ` is bounded.
- Lemma 4.3 contains a sign slip: it prints `S_k = ℑ(h+h′, P_H i log Δ (h+h′)) = … = ℑ(h, i log Δ h) = −S_h`. With (10) the first expression carries a minus sign and the chain ends in `S_k = S_h`. The conclusion of the lemma is right.
- Lemma 4.4 (tensor additivity of relative entropy): asserted, not proved, not cited.
- Thm 4.5: `S(ϕ_k‖ϕ) = S(ϕ_{k_ε}‖ϕ) + S(ϕ_{(1−E_ε)k}‖ϕ) ≥ S(ϕ_{k_ε}‖ϕ)` (Lemma 4.4 and positivity), then "by the semicontinuity of the relative entropy" `S(ϕ_k‖ϕ) ≤ lim inf S(ϕ_{k_ε}‖ϕ)`, hence `S(ϕ_k‖ϕ) = lim S(ϕ_{k_ε}‖ϕ) = lim S_{k_ε} = S_k`. Lower semicontinuity is imported without citation (it is Araki 1976 property (2), eq (1.4), for faithful functionals and norm convergence; `V(k_ε)ξ → V(k)ξ` in norm, and all states involved are faithful on `R(H)`). I read the statement of (1.4) and the first lines of its proof only.
- Remark (mine, not printed): the inequality `S(ϕ_k‖ϕ) ≥ S(ϕ_{k_ε}‖ϕ)` also follows from I3 alone, without Lemma 4.4: `R(H_ε) ⊂ R(H)` because `E_ε H ⊂ H`, and `ϕ_k` restricted to `R(H_ε)` equals `ϕ_{k_ε}` because `ℑ(k,h) = ℑ(k_ε,h)` for `h ∈ H_ε`. So the imports of Thm 4.5 reduce to: Araki's formula, `Δ_{R(H)} = Γ(Δ_H)`, restriction monotonicity, lower semicontinuity.

### How `k` outside the domain is handled

There is no domain condition in Thm 4.5. `S_k` is defined for every `k ∈ 𝓗` by the spectral integrals (10); it equals `+∞` exactly when `k ∉ Dom(√|log Δ| E_−)` (Prop 2.4), and the theorem then says the relative entropy is `+∞`. For `k ∈ H + H′`, `S_k = −(h, log Δ h)` with `h = P_H k` (Prop 2.5(a)). For `k ∈ Dom(log Δ)`, `log Δ k ∈ H + H′` (Cor 2.3) and the operator expression `−ℑ(k, P_H i log Δ k)` is literal: `|a(λ) log λ|` and `|b(λ) log λ|` are bounded by `C(1 + |log λ|)`, so `log Δ k ∈ D` and Thm 2.2 applies. For general `k`, `S_k = lim_ε S_{E_ε k}`.

### Fit to our setting

- `K = K(D(B))` standard and factorial: both hypotheses met as stated in the brief.
- `k` arbitrary in the one-particle space: covered.
- Needs from elsewhere: `Δ_{R(K)} = Γ(Δ_K)` (I2), lower semicontinuity (Araki 1976), monotonicity or tensor additivity.
- The formula depends on `k` only through `k` mod `K′` and needs only the easy inclusion `R(K′) ⊂ R(K)′`.

### Traps

1. **Weyl normalisation.** CLR: `V(h)V(k) = e^{−iℑ(h,k)} V(h+k)`, vacuum `e^{−‖h‖²/2}`. Brunetti–Guido–Longo (Remark 3.6, pdf p9) print `V(h)V(k) = e^{−(i/2)ℑ(h,k)} V(h+k)`, vacuum `e^{−‖h‖²/4}`. If our `W(k)` has the BGL normalisation, then `W(k) = V(k/√2)` and the entropy is HALF of `S_k`. Fix our convention against `ω(W(k)) = e^{−‖k‖²/2}` before quoting the formula.
2. **Sign.** With the inner product antilinear in the first slot the correct formula is `S_k = −ℑ(k, P_H i log Δ_H k)` (CLR (1), Thm 4.5). Longo–Morsella eq (7) and (17) print `+ℑ(k, P_H i log Δ_H k)` with the same inner-product convention; that is a misprint there (it is negative on `H`). See I4.
3. **Argument order.** CLR's `S(ϕ1‖ϕ2)` has the reference state second. Araki 1976 writes the same number as `S(ϕ2/ϕ1)` (reference first); Araki 1977 writes `S(ψ/φ)` with the reference `ψ` first as well. CLR eq (30) and the Introduction print `−(ξ, log Δ_{V(h)ξ,ξ} ξ)`, which by their own (23) is `S(ϕ‖ϕ_h)`; this is harmless only because `S(ϕ_h‖ϕ) = S(ϕ‖ϕ_{−h})` and `S_k` is even in `k`.
4. `P_H` is unbounded; `−ℑ(k, P_H i log Δ k)` is a quadratic form, not an operator expression, unless `k ∈ Dom(log Δ)`.
5. Eq (29) prints `S(ϕ_h‖ϕ_k) = S(ϕ_{k−h}‖ϕ)`; the automorphism algebra gives `ϕ_{h−k}`. Same value, by evenness.

---

## I2. Araki's duality and the second-quantisation lattice facts

### Garbarz–Palau, "A note on Haag duality", arXiv 2108.01257, §5.3 (pdf pp26–34, read line by line)

What is actually proved. Setting: a REAL Hilbert space `K`, two closed subspaces `K1, K2 ⊂ K` in generic position (`K1 ∩ K2 = K1 ∩ K2^⊥ = K1^⊥ ∩ K2 = K1^⊥ ∩ K2^⊥ = {0}`), fields `ϕ(f)`, `π(g)` on the Fock space, `R_W(K1,K2) = {e^{iϕ(f)}, e^{iπ(g)} : f ∈ K1, g ∈ K2}″`. Claim (eq (49), p24 and the discussion after (54), p27): `R_W(K1,K2)′ = R_W(K2^⊥, K1^⊥)`. By eq (47), `R_W(K1,K2) = R_S(K1 ⊕ βK2)`, so this is the duality `R(K)′ = R(K′)` for real subspaces of the special "position ⊕ momentum" form with respect to a fixed conjugation. It is the Eckmann–Osterwalder theorem, not the statement for an arbitrary standard subspace.

Chain of the proof:
- Thm 5.2 (Tomita–Takesaki): cited to Bratteli–Robinson and Takesaki. Fine as an import.
- `Ω` cyclic for `R_W(S_R(B), S_I(B))`: from Reeh–Schlieder, cited; separating by footnote 23, which needs cyclicity for the complementary pair as well. Cyclic-and-separating is an assumption of the abstract argument ("from now on we always assume this").
- Lemma 5.3 (`J`, `S`, `S^*` commute with the particle-number projections `E^{(k)}`; `S` preserves the real one-particle space): "The part of this Lemma referring to `S` and `J` is proven in [7]" = Eckmann–Osterwalder. **Not proved in the file.** Only the `S^*` part is proved.
- Lemma 5.4 (`K1 + K2 ⊂ Dom(S)`, eq (57)): proved, by differentiating Weyl unitaries.
- Thm 5.5 (Halmos two projections): cited.
- Lemma 5.6 (`Ĵ K1 = K2^⊥`, `Ĵ K2 = K1^⊥`): computed from the matrix form of `S` on `K1 + K2`. **Soft spot, acknowledged in print (p29):** "In [7] the authors claim that `K1 + K2 = D ∩ F_r^{(1)}`, however we do not see how to prove `D ∩ F_r^{(1)} ⊆ K1 + K2`. In any case it will be sufficient with what we have claimed above." It is not sufficient as written: the polar decomposition of `S` restricted to the one-particle space is computed from its restriction to `K1 + K2`, which requires `K1 + K2` to be the full one-particle domain (a core in the graph norm is what is needed, density in the Hilbert norm is what is shown). The step can be closed (the one-particle Tomita operator of the commutant contains the analogous operator for `(K2^⊥, K1^⊥)`, and `S ⊂ (S′)^*`), but that argument is not in the file.
- Lemma 5.7 and Cor 5.7.1 (`J ϕ(f) J = ϕ(Ĵf)`, hence `J e^{iϕ(f)} J = e^{iϕ(Ĵf)}`): proved in detail (affiliation, Nelson analytic vectors). Read; no gap found beyond the reliance on Lemma 5.6.

Verdict: READ; a proof of the pair-form duality that DEFERS-TO Eckmann–Osterwalder (Lemma 4 there) for one lemma and has one self-flagged gap. Additivity `R_S(∨Y_λ) = (∪ R_S(Y_λ))″` is cited ("see theorem 1 in [14]", Araki 1963, p24). The intersection property, `Δ = Γ(Δ_K)` and "cyclic and separating iff standard" are not proved; footnote 25 (p29) mentions standardness without proof.

Fit: our `K ⊂ H` is a standard subspace of a complex Hilbert space given by smeared curvature of helicity ±2 data, with no evident splitting `K = K1 ⊕ βK2` in generic position with respect to a conjugation. Garbarz–Palau does not cover it without an extra argument.

### Camassa, math-ph/0607057, Thm 2.1 (pdf pp2–3)

Printed statement: "Theorem 2.1. (Araki) Let `H0` be a complex Hilbert space, `V, W, V_α ∈ subsp(H0)` (closed real subspaces of `H0`), then 1. `R(V) ⊃ R(W) ⟺ V ⊃ W` and `R(V) = R(W) ⟺ V = W` (isotony), 2. `∨_α R(V_α) = R(∨_α V_α)` (additivity), 3. `∧_α R(V_α) = R(∧_α V_α)`, 4. `R(V)′ = R(V′)` (abstract duality). It was originally proven in [1], Th. 1; see also [22, 3]." ([1] Araki 1963, [22] Leyland–Roberts–Testard, [3] Araki 1964.) This is exactly the general statement we need, for arbitrary closed real subspaces of a complex Hilbert space, symplectic form `ℑ(·,·)`. **No proof in the file.** READ-STATEMENT-ONLY.

### Brunetti–Guido–Longo, math-ph/0203021 (pdf pp8–9)

Cor 3.5 is one-particle duality `K_{O′} = K_O′` for the BGL net, proved from wedge duality (their Thm 2.5, not read). Remark 3.6(1): definition of `R(O)`, the Weyl relations in the BGL normalisation (trap 1 of I1), and "The standard property of `K_O` is equivalent to the Reeh–Schlieder property for `R(O)` (cf. [1, 15, 33])" (Araki, Eckmann–Osterwalder, Leyland–Roberts–Testard). **Citation only.** `Δ_{R} = Γ(Δ_K)` is not stated on these pages.

### Which facts have a complete proof in a file I read

| Fact | Proved in a file read? |
|---|---|
| `R(K′) ⊂ R(K)′` (locality) | Elementary from the Weyl relations; used as such in CLR and Garbarz–Palau |
| `R(K)′ = R(K′)`, general closed real `K` | No. Camassa Thm 2.1(4), CLR p19: statement only. Garbarz–Palau: pair form only, with the two reservations above |
| `R(∨K_a) = ∨R(K_a)` | No. Camassa says "immediate consequences of the definition" for isotony and additivity (Prop 2.4); Garbarz–Palau cites Araki 1963 Thm 1 |
| `R(∧K_a) = ∧R(K_a)` | No. Camassa Thm 2.1(3): statement only. (Given duality for all closed real subspaces and additivity, it follows by taking commutants; that derivation is not printed in the files read.) |
| `Ω` cyclic and separating for `R(K)` iff `K` standard | No. CLR p19 (one direction), BGL Remark 3.6, Garbarz–Palau fn 25: statements only |
| `Δ_{R(K)} = Γ(Δ_K)`, `J_{R(K)} = Γ(J_K)` | No. CLR p19: "we recall", no proof, no citation. Not in Camassa, not in the BGL pages read, not in Garbarz–Palau |

Consequence for the final theorem: the second-quantisation facts are either a named hypothesis ("Araki 1963 / Eckmann–Osterwalder 1973 / Leyland–Roberts–Testard 1978, not on file") or need a proof in the repo. What CLR Thm 4.5 actually consumes is small: `Ω` cyclic and separating for `R(K)`, `Δ_{R(K)}^{is} = Γ(Δ_K^{is})`, and the easy inclusion. Remark (mine, not printed, unchecked beyond this sketch): `Δ^{is}_{R(K)} = Γ(Δ^{is}_K)` has a short direct proof by the KMS condition on Weyl operators, since for `h, k ∈ K` the function `s ↦ (h, Δ_K^{s}k) = (Δ^{s/2}h, Δ^{s/2}k)` is defined on `[0,1]` with value `(k,h)` at `s = 1`, and Takesaki's theorem identifies the modular group once `Ω` is separating. That route does not need full duality.

---

## I3. Uhlmann monotonicity

### Reible, arXiv 2501.04564, Thm IV.3.7 (pdf pp68–71; read: Def IV.1.1, Lemma IV.1.3, Def IV.2.1, Prop IV.2.3, Thm IV.2.4, Lemmas IV.3.1–IV.3.6, Thm IV.3.7 with proof, Remarks IV.3.8, Cor IV.3.9–IV.3.12, Remark after Prop IV.4.2)

Printed statement: "For `i ∈ {1,2}`, let `M_i ⊆ B(H_i)` be a von Neumann algebra acting on a Hilbert space `H_i` and `ϕ_i, ω_i ∈ (M_i)_*^+`. Assume that `ω_i = ω_{Ω_i}` is given by a cyclic vector `Ω_i ∈ H_i` for the algebra `M_i`, and that `ω_1(Id) = ω_2(Id)`. Let `α : M_1 → M_2` be a Schwarz mapping such that `ω_2∘α ≤ ω_1` and `ϕ_2∘α ≤ ϕ_1` (IV.14). Finally, assume that at least one of: (U1) the vector `Ω_i` is separating for `M_i`; (U2) the mapping `α` is 2-positive. Then `S^spa_{M_1}(ω_1, ϕ_1) ≤ S^spa_{M_2}(ω_2, ϕ_2)` (IV.15)."

Conventions and hypotheses:
- `S^spa_M(ω, ϕ) = −⟨Ω, log Δ(ϕ/ω′) Ω⟩` if `Ω ∈ supp(ϕ)`, `+∞` otherwise (Def IV.2.1); `Δ(ϕ/ω′)` is Connes' spatial derivative. `ω` (first argument) is "the state" and must be a vector state of a CYCLIC vector; `ϕ` (second argument) is the reference and is an arbitrary normal positive functional, faithful or not. Same order as CLR's `S(ω‖ϕ)`. In standard form, `S^std(ψ,ϕ) = −⟨Ψ, log Δ_{Φ,Ψ} Ψ⟩` (Def IV.1.1), which is CLR (23); `S^std = S^spa` is printed for `ψ` faithful (Prop IV.2.3, resting on Prop III.4.13, whose proof is cited to Hiai Prop 11.6).
- Direction: the map goes `M_1 → M_2`, the states live on both, and entropy on the SOURCE is ≤ entropy on the TARGET. For an inclusion `M_1 ⊂ M_2` (Cor IV.3.9): `S_{M_1}(ω|M_1, ϕ|M_1) ≤ S_{M_2}(ω, ϕ)`, with `Ω` cyclic for the smaller algebra `M_1`.
- "Schwarz mapping" is by Def II.2.10 UNITAL and satisfies `α(A)^*α(A) ≤ α(A^*A)`. Unitality is used twice (Lemma IV.3.6 and `TΩ_1 = α(1)Ω_2 = Ω_2`); the remark after Prop IV.4.2 says that under (U1) it can be weakened to `α(1)Ω_2 = Ω_2`.
- Normality of `α` is not assumed (Reible's introduction describes Uhlmann's original as for unital normal Schwarz maps); the hypotheses `ω_2∘α ≤ ω_1`, `ϕ_2∘α ≤ ϕ_1` with normal functionals do the work.
- No faithfulness of `ϕ_i`; no σ-finiteness or separability stated. Faithfulness of `ω_i` enters only through (U1).

Proof status. Read line by line: Statement 1 (`T(AΩ_1) = α(A)Ω_2` is a well-defined contraction), Statement 2 (`T Dom(Δ_1^{1/2}) ⊂ Dom(Δ_2^{1/2})`, `‖Δ_2^{1/2}Tξ‖ ≤ ‖Δ_1^{1/2}ξ‖`, in both cases (U1) and (U2)), interpolation (Lemma IV.3.4) to exponent `t/2`, then Uhlmann's representation (Thm IV.2.4) `S = −lim_{t↘0} t^{−1}(‖Δ^{t/2}Ω‖² − ‖Ω‖²)`. I found no error. It is self-contained in the file EXCEPT:
- Lemma IV.3.1 (Löwner–Heinz, unbounded form): cited to Schmüdgen Prop 10.14, not proved.
- Lemma IV.3.2 (bounded Hansen–Jensen–Pedersen): cited to Bhatia Thm V.2.3 and Ohya–Petz Lem 1.2, not proved.
- The spatial-derivative facts used (`M_iΩ_i` is a core for `Δ_i^{1/2}`; `‖Δ_i^{1/2}AΩ_i‖² = ϕ_i(Θ(AΩ_i))`; `Θ(AΩ) = A[M′Ω]A^*`, Prop III.4.10; support of `Δ(ϕ/ω′)` is `s(ϕ)`): §III.4, partly proved there, partly cited to Hiai's lecture notes. I skimmed §III.4 and did not read it line by line.
- Minor: Lemma IV.3.3 is stated for an arbitrary contraction `K`, for which `Dom(T^{1/2}K)` need not be dense; in its only use (Lemma IV.3.4) the domain is dense by hypothesis (α). I did not read Para. C.2.1, which fixes the meaning of `≤` for unbounded forms.
- Remark IV.3.8(3) says the proof in Ohya–Petz Thm 5.3 "seems to contain some inaccuracies" (unstated separating/standard-form assumptions). Uhlmann's original (CMP 54, 1977, Prop. 18) is not on file.

Fit to our setting. For `K_1 ⊂ K_2` standard and `α` the inclusion `R(K_1) ⊂ R(K_2)` (completely positive, unital): `ω = ω_k`, vector `W(k)Ω`. This vector is cyclic AND separating for every `R(K)` with `K` standard and every `k ∈ H`, because `W(k) R(K) W(k)^* = R(K)` (Weyl relations) and `Ω` is cyclic and separating. So both (U1) and (U2) hold, and `S_{R(K_1)}(ω_k‖ω) ≤ S_{R(K_2)}(ω_k‖ω)`. The coherent state must sit in the FIRST slot (it carries the cyclic vector); the vacuum is the reference. All hypotheses are met, given "Ω cyclic and separating for `R(K)`" (I2).

### Araki 1976 (page images pp1–2 read)

(1.1) `S(φ/ψ) ≡ −(Ψ, (log Δ_{Φ,Ψ})Ψ)` for normal FAITHFUL positive functionals; (1.2) `S(ρ_φ/ρ_ψ) = tr(ρ_ψ log ρ_ψ) − tr(ρ_ψ log ρ_φ)`. So Araki's `S(φ/ψ)` is the modern `S(ψ‖φ)`: reference state FIRST. Property (2), eq (1.4): lower semicontinuity under norm convergence, `lim inf S(φ_n/ψ_n) ≥ S(φ/ψ)`. Property (4), eq (1.6): `S(E_𝔑φ/E_𝔑ψ) ≤ S(φ/ψ)` for the restriction to a subalgebra `𝔑`, **only** when (α) `𝔑 = 𝔄′ ∩ 𝔐` with `𝔄` finite-dimensional abelian, (β) `𝔐 = 𝔑 ⊗ 𝔑_1`, or (γ) `𝔑` approximately finite. Proofs (pp3 on) not read.

### Araki 1977 (page images pp5–9 read)

Def 3.1: `S(ψ/φ) = −∫ log λ d(Φ, E_λ^{Φ,Ψ}Φ)` if `s(ψ) ≥ s(φ)`, `+∞` otherwise; (3.3) `S(ψ/φ) = −(Φ, log Δ_{Ψ,Φ} Φ)`; eq (1.1) `S(ψ/φ) = φ(log ρ_φ − log ρ_ψ)`. Again the reference `ψ` is the first argument and the support condition is on the reference. Standing assumption: `M` has a faithful normal state (σ-finite).
- Thm 3.6(4) is monotonicity in the FUNCTIONAL: `ψ_1 ≥ ψ_2 ⇒ S(ψ_1/φ) ≤ S(ψ_2/φ)`. Proof read (p179–180), complete. It is not the subalgebra statement. `sources/INDEX.md` line 15 describes item (4) as "the monotonicity step `S(E_Nφ‖E_Nψ) ≤ S(φ‖ψ)`"; that description is wrong.
- Thm 3.7 (p180): lower semicontinuity for non-faithful functionals; "We shall give proof of this Theorem in the next section" (not read).
- Thm 3.8(2), eq (3.9) (p181): restriction monotonicity `S(E_Nψ/E_Nφ) ≤ S(ψ/φ)`, again only for the three types (α), (β), (γ); the proof reduces to the faithful case of Araki 1976 by Thm 3.7.

Verdict: Araki's two papers do not give restriction monotonicity for a general subalgebra. Case (γ) would need "`R(K_1)` is approximately finite", a separate fact. Use Reible Thm IV.3.7 / Cor IV.3.9 for the import, and Araki 1976 (1.4) for lower semicontinuity.

Traps: (i) the argument order flips between Araki and CLR/Reible; (ii) Reible's theorem is stated for `S^spa` with a cyclic vector in the first slot, CLR uses the relative-modular formula, and the identification is Prop IV.2.3 (first state faithful, satisfied here); (iii) the direction: restriction to the smaller algebra lowers the entropy.

---

## I4. Longo–Morsella, "The massless modular Hamiltonian", arXiv 2012.00565

File: `sources/longo-morsella-2020-massless-modular-hamiltonian-2012.00565.pdf`; pp2–21 read.

### Setting and conventions

Free massless real scalar, `d ≥ 2` space dimensions. `S` = real smooth compactly supported functions on `R^d`; `T` = waves `Φ = w(f,g)` with Cauchy data `f = Φ|_{x0=0}`, `g = ∂_0Φ|_{x0=0}` in `S`. `β(Φ,Ψ) = ½∫(Ψ∂_0Φ − Φ∂_0Ψ)dx = ½∫(f_2 g_1 − f_1 g_2)dx` (24),(38). Complex structure `ı_0⟨f,g⟩ = ⟨μ^{−1}g, −μf⟩`, `μ = |p|`; `(Φ,Ψ) = β(Φ, ı_0Ψ) + iβ(Φ,Ψ)`, which is antilinear in the first slot (check: `(Φ, ı_0Ψ) = i(Φ,Ψ)`), same as CLR. `H = H^{1/2} ⊕ H^{−1/2}` (homogeneous Sobolev spaces). `H(B)` = closure of `T(B)` (data supported in `B`) `= H^{1/2}(B) ⊕ H^{−1/2}(B)`. `H(B)′ = H(B′)` (29) "by duality [1], see also [12, 20]" (Araki, Eckmann–Osterwalder, Leyland–Roberts–Testard): cited. `D = (d−1)/2`. `⟨T_00⟩_Φ = ½((∂_0Φ)² + |∇Φ|²)`.

### Printed statements

- Thm 4.1 (p12): "The modular group of `H(B)` is given by `Δ_B^{−is} = V(2πs)`", with `(V(s)Φ)(u,v) = γ(u,v;s) Φ(Z(u,s), Z(v,s))`, `u = x0 + r`, `v = x0 − r`, `Z = g/f`, `f(z,s) = ((1+z) + e^{−s}(1−z))/2`, `g(z,s) = ((1+z) − e^{−s}(1−z))/2`, `γ = F(u,s)F(−v,−s)`, `F = f^{−D}`.
- Prop 4.4 (p14): `K = d/ds V(s)|_0` acts on `S²` as `K⟨f,g⟩ = ⟨½(1−r²)g, ½(1−r²)∇²f − r∂_r f − Df⟩` (37).
- Prop 4.5 (p15): with `A = −ı_0K`, `(Φ, AΦ) = β(Φ, KΦ) = ½∫(1−r²)⟨T_00⟩_Φ dx + (D/2)∫Φ² dx` (41), integrals over all of `R^d`, `Φ ∈ T`.
- Thm 4.7 (p16): "`log Δ_B` is given by `−2πA = log Δ_B` with `K = [[0, M],[L, 0]] = [[0, M],[−M^*, 0]]` and `A = −ı_0K`", `M` = closure of multiplication by `½(1−r²)`, `L` = closure of `½(1−r²)∇² − r∂_r − D`.
- Thm 4.12 (p17): "`−ℜ(Φ, log Δ_R Ψ) = 2π∫_{B_R} (R²−r²)/(2R) ⟨T_00⟩_{Φ,Ψ} dx + 2π (d−1)/R ∫_{B_R} ΦΨ dx`, where `Φ, Ψ ∈ T(B_R)`." Proof: "Immediate by the above discussion."
- Prop 5.1–5.3 (pp17–18): `P_B` is multiplication of both Cauchy data by `χ_B`; for `Φ ∈ T`, `Ψ = KΦ ∈ H(B) + H(B′)` and `P_BΨ = w(χ_B f, χ_B g)`; `S_Φ = π∫_B(ΨΦ′ − ΦΨ′)dx` (47), `Ψ = KΦ`.
- Thm 5.4 (p18): "The entropy `S_Φ` of the wave `Φ ∈ T` in the unit ball `B`, i.e. with respect to `H(B)`, is given by `S_Φ = 2π∫_B (1−r²)/2 ⟨T_00⟩_Φ dx + πD∫_B Φ² dx` (48). In particular, `S_Φ` is finite for all `Φ ∈ T`."
- Cor 5.5 (p19): `S_Φ(R) = π∫_{B_R(x̄)} (R²−r²)/R ⟨T_00⟩_Φ dx + (πD/R)∫_{B_R(x̄)} Φ² dx` at time `t`.
- Cor 5.6 (p20): `S(ϕ_Φ‖ϕ)` on `A(O_R)` equals `S_Φ(R)` for the one-particle wave `Φ`; proof: "As shown in [11]" (CLR).

### What Thm 5.4 assumes about the data

`Φ ∈ T`: real, smooth, COMPACTLY SUPPORTED Cauchy data anywhere in `R^d`, not required to be localised in `B` and not required to vanish on the sphere. The integrals in (48) are over `B` only.

### The coefficient: `πD/R` is right, Thm 4.12 is wrong by a factor 4

Redone from (37) and (38):

`2β(Φ, KΦ) = ∫[(KΦ)∂_0Φ − Φ∂_0(KΦ)] = ∫½(1−r²)g² − ∫½(1−r²) f∇²f + ∫ r f∂_r f + D∫f²` (this is (39)).

One integration by parts: `−∫½(1−r²) f∇²f = ∫∇(½(1−r²)f)·∇f = ∫½(1−r²)|∇f|² − ∫ r f∂_r f`. The `r f∂_r f` terms cancel, giving

`β(Φ,KΦ) = ½∫(1−r²)·½(g² + |∇f|²) + (D/2)∫f²`, i.e. (41),

and `−(Φ, log Δ_B Φ) = 2π(Φ,AΦ) = 2π∫ ½(1−r²)⟨T_00⟩ + πD∫Φ²`: coefficient **`πD = π(d−1)/2`** at `R = 1`, as in (5), Thm 5.4.

Radius `R`: with `δ_R Φ = R^D Φ(R·)` (data `R^D f(R·)`, `R^{D+1} g(R·)`), `∫_{B_1}(δ_RΦ)² = R^{2D−d}∫_{B_R}Φ² = R^{−1}∫_{B_R}Φ²` and `∫_{B_1}½(1−ρ²)⟨T_00⟩_{δ_RΦ} = ∫_{B_R} (R²−r²)/(2R) ⟨T_00⟩_Φ`. Hence the second term is `(πD/R)∫_{B_R}Φ² = π(d−1)/(2R)∫Φ²`, as in Cor 5.5 and the Introduction (p4). Thm 4.12's `2π(d−1)/R` should read `2π(d−1)/(4R) = πD/R`. **Thm 4.12 as printed is wrong by a factor 4 in the `ΦΨ` term; Thm 5.4 and Cor 5.5 are right.**

Numerical check run here (pure Python, Simpson, `d = 3`, radial data `f = e^{−r²}(1+r²/3)`, `g = e^{−2r²}(2−r²)`, neither vanishing on the sphere): `π∫_{B_R}((K_RΦ) g − f ∂_0(K_RΦ))` against the two candidate right sides.

| R | left side | with `πD/R` | with `2π(d−1)/R` |
|---|---|---|---|
| 1.0 | 9.36422314 | 9.36422315 | 27.67327646 |
| 1.5 | 16.04612232 | 16.04612231 | 34.52678984 |
| 0.4 | 2.05516906 | 2.05516906 | 7.61518682 |

For `R ≠ 1` the operator `K_R` was taken from the scaling (45), so those two rows test the scaling step; the `R = 1` row tests the printed (37) directly.

### How data not vanishing on the sphere are handled

Printed route (proof of Thm 5.4, p18): define `q_B` by the right side of (48); assert `q_B(Φ,Ψ) = ℜ(Φ, P_B^* A P_B Ψ)` "similarly as in [11, Thm. 3.5]"; assert closability; "(49) holds if `Φ ∈ T(B)` by Theorem 4.12. Then it holds if `Φ ∈ D(log Δ_B)` by the same argument, using Prop. 5.3. As `D(log Δ_B)` is a form core for `S` … we conclude that `S = q_B` on `T`." This is a sketch: it leans on the (misprinted) Thm 4.12, on Prop 5.3 whose own proof is "in analogy with the proof given in [11]", and on Prop 5.1 whose proof is three lines ("`P_B` acts as in (46) on vectors in `H(B) + H(B′)` given by smooth functions, hence on all `D(P_B)` because `P_B` is a closed operator"). There is also a typo `P_0^*AP_B`.

Direct route (my check; every step elementary, given Thm 4.7 on `T`): for `Φ = w(f,g) ∈ T`, `KΦ = w(½(1−r²)g, Lf)`.
- The first component `½(1−r²)g` VANISHES on the unit sphere, so `χ_B·½(1−r²)g` is Lipschitz with compact support, lies in `H^{1/2}(B)`, and so does the exterior part in `H^{1/2}(B′)`. (A cut `χ_B f` of a smooth `f` that does not vanish on the sphere is NOT in `H^{1/2}`; the weight `1−r²` is what saves the cut.) The second component `Lf` is bounded with compact support, and `χ_B Lf ∈ L²(B) ⊂ H^{−1/2}(B)` for `d ≥ 2`.
- So `KΦ = w(χ_B·, χ_B·) + w(χ_{B′}·, χ_{B′}·) ∈ H(B) + H(B′)` explicitly, with `H(B′) ⊂ H(B)′` by locality (the easy inclusion), and `P_B KΦ = w(χ_B½(1−r²)g, χ_B Lf)` by the definition of the cutting projection.
- `S_Φ = −ℑ(Φ, P_B ı_0 log Δ_B Φ) = 2π β(Φ, P_B KΦ) = π∫_B[½(1−r²)g² − f(½(1−r²)∇²f − r∂_r f − Df)]dx`, using `ı_0 log Δ_B = −2πK`.
- Integration by parts ON THE BALL: `−∫_B ½(1−r²) f∇²f = ∫_B ½(1−r²)|∇f|² − ∫_B r f∂_r f − ∮_{∂B} ½(1−r²) f ∂_n f`, and the surface term is zero because `1−r² = 0` on the sphere, whatever `f` and `∂_n f` are there. Result: (48), with no boundary term.

So Thm 5.4 holds for all `Φ ∈ T` by a direct computation; the mechanism for non-vanishing boundary data is the zero of `1−r²` on the sphere, twice (membership of the cut vector in `H^{1/2}`, and the vanishing surface term).

### What is imported

- Thm 4.1: DEFERS-TO Hislop–Longo [18]. The printed proof: "That `Δ_B^{−is/2π}` is associated with the flow `Z` (30) and the cocycle `γ` follows rather directly from the test function formula [18], but for the determination of the value of the constant `D`"; `D` is then fixed by complex linearity (Cor 4.8). Not a self-contained proof.
- `H(B)′ = H(B′)` (29): cited (Araki / Eckmann–Osterwalder / Leyland–Roberts–Testard). Needed for `Dom(P_B) = H(B) + H(B′)` (Prop 5.1); the direct route above needs only `H(B′) ⊂ H(B)′`.
- CLR for `P_H`, the entropy form and Cor 5.6.
- Lemma 4.6 (`S²` is a core for `K`, hence `T ⊂ D(K)` with `K` acting as (37)): the printed proof says `C_0^∞(B)² + C_0^∞(B′)²` is "globally `Δ_B^{is}`-invariant". For the exterior part this is not literally true: the conformal flow moves spatial infinity, so compact support in `B′` is not preserved. The interior part `T(B)` is invariant. The conclusion can be recovered without the exterior: `T(B) + ı_0T(B)` is dense, invariant and in the domain, hence a core for `K`; and the differential expression (37) is skew for the full complex inner product on `S²`, which is the identity `μMμ = −L`. I verified that identity in Fourier space: `½|p|(1+Δ_p)|p| = ½p² + p·∇_p + ½p²Δ_p + (d−1)/2` and `−L̂ = ½p² + p·∇_p + ½p²Δ_p + D`, equal iff `D = (d−1)/2`. Then every `Φ ∈ T` is in `Dom(K^*) = Dom(K)` and `KΦ` is given by (37). This repair is mine, not printed.

### Fit to our setting

Scalar field only. For the helicity-2 net these theorems are the template, not the statement; the corresponding facts must come from the repo's own proofs (wave70 modflow for the modular group; the analogue of (37) and of the vanishing-weight mechanism for the entropy formula). What transfers verbatim is the abstract part: CLR's `S_k = −ℑ(k, P_K i log Δ_K k)` and the fact that a generator whose "position" component carries a weight vanishing on the sphere makes the cut vector lie in `K + K′`.

### Traps

1. Thm 4.12's coefficient (factor 4), above. Do not cite Thm 4.12 for the `Φ²` term.
2. Sign of the entropy formula: eq (7) and (17) print `S_k = +ℑ(k, P_H i log Δ_H k)`; with their (and CLR's) inner product that is `≤ 0` on `H`. The correct sign is CLR's minus. Eq (47) with `Ψ = KΦ` is consistent with CLR's sign; eq (8) in the Introduction says `Ψ = i log Δ_B Φ`, which is `−2πKΦ`, inconsistent with (47). Likewise p8 says "`a(λ) log λ` is … positive on (0,1)" (it is negative).
3. Sign in eq (1)–(2): they print `log Δ_B = −2πi [[0,M],[L,0]]`, i.e. `−2π ı_0 K`. Thm 4.7 gives `log Δ_B = −2πA = +2π ı_0 K`. Positivity of `−(Φ, log Δ_B Φ) = 2πβ(Φ,KΦ)` on `T(B)` confirms Thm 4.7; eq (1) has the opposite sign.
4. `Δ_B^{−is} = V(2πs)`: the minus sign in the exponent matters for (3).
5. `β` is `½∫(…)`, so `ℑ(Φ,Ψ)` carries a factor ½ relative to the usual symplectic form `∫(f_1g_2 − f_2g_1)` and has the opposite sign.

---

## I5. Hislop–Longo 1982, CMP 84, 71–85

File: `sources/hislop-longo-1982-modular-structure-free-massless-scalar-cmp84.pdf` (16 pages, all read; the text layer drops most displayed formulas, so pp5–6 were read as images; the two-point-function formulas of Appendix 1 were not legible in the text layer and were not viewed as images).

### Printed statements

- Thm 1 (p4): the ray-inversion operator `U_ρ` (`ρ(x) = −x/x²`, `f_ρ(x) = −(x²)^{−3} f(ρ(x))` in `d = 3`) "extends to a unitary operator on `𝓗`". Proved in §5 (Lemmas 1, 2) from the two-point functions (A1)–(A5), citing Swieca–Völkel for the time-zero unitarity.
- Cor 1 (p5): "`R(O_1)`, `R(V_+)` and `R(W_R)` are spatially isomorphic. The equivalences are given by unitary operators on `Γ(𝓗)` induced by unitary operators on `𝓗`."
- Quoted input (p5): "For `O = V_+`, the modular operators `{J_{V_+}, Δ^{iλ}_{V_+}}` in the vacuum representation were shown by Buchholz [2] to be given by `Δ^{−iλ}_{V_+} = Γ[D(2πλ)]`, `J_{V_+} = Γ(−Θ)`", `D(λ)` the dilatations.
- Thm 2 (p6): "`V(λ)`, `λ ∈ R`, is a one-parameter strongly continuous unitary group on `𝓗`. The modular operators for `R(O_1)` are: `Δ^{−iλ}_{O_1} = Γ[V(2πλ)]`, `λ ∈ R`, `J_{O_1} = Γ[I_t U_ρ]`", where `V(λ)φ[f]Ω = φ[f_λ]Ω`, `f_λ(x) = γ(x;λ) f(v_{−λ}(x))`, `v_λ` the flow `z_i ↦ ((1+z_i) − e^{λ}(1−z_i))/((1+z_i) + e^{λ}(1−z_i))` in `z_{1,2} = x_0 ± x_3`, and `γ(z_1,z_2;λ) = 2^{6}[(1+z_1) + e^{−λ}(1−z_1)]^{−3}[(1−z_2) + e^{λ}(1+z_2)]^{−3}` for `d = 3` test functions.
- Cor 2 (p6): "The local algebra `R(O_1)` satisfies spacelike duality: `R(O_1)′ = R(O_1^c)`."

### Proof status and imports

- Thm 2's proof: conjugate Buchholz's forward-cone result by `Γ[T(1)U_ρT(½)]` and check two identities between geometric operators, "the theorem follows by elementary computation" (the computation is not written out). Imports:
  - **Buchholz 1978** (ref [2], Leipzig conference proceedings): modular group of `R(V_+)` is the dilatations, `J_{V_+} = Γ(−Θ)`, and timelike duality `R(V_+)′ = R(V_−)`, for odd spatial dimension `d > 1`. CITED, not proved in the file.
  - Thm 1 (proved in the file; depends on Appendix 1 formulas and Swieca–Völkel [4] for Lemma 1).
- **Bisognano–Wichmann** (refs [1], [3]) is CITED in the Introduction and in the remark that `v_λ` is conformally conjugate to boosts; it is not proved in the file and is not what the proof of Thm 2 uses in odd `d`. The wedge enters through Cor 1 only as an equivalence.
- Even `d` (end of §4, pp9–10): timelike commutativity fails and `U_ρ` does not preserve the real one-particle space; the modular objects of `R(V_+)` are said to be obtainable from the wedge through `U_ρ` ("the modular operators for `R(V_+)` and `R(V_−)` can be obtained"), with a twist `Γ(Z)`, `Z = i`. This is a sketch; Thm 4 states the twisted timelike duality and `J_{O_1} = Γ[Z I_t U_ρ]`. For even `d` the argument therefore rests on Bisognano–Wichmann, cited.
- Cor 2's proof: `R(O_1)′ = J R(O_1) J = J R(B) J = R[ρ(B)] = R(O_1^c)` with `B = O_1 ∩ S`, using Lemma 1 of Appendix 2 (`R(O_1) = R(O_1 ∩ S)`), which itself imports Araki 1964 [13] (`R(O_1) = R_0(B_1)`, the time-zero algebra). Complete given Thm 2 and that import.
- The paper works directly on Fock space with fields `φ[f]`; it does not go through `Δ_{R(K)} = Γ(Δ_K)`, because Buchholz's input is already second-quantised.

Verdict: READ-COMPLETE-MODULO Buchholz 1978 (odd `d`), Araki 1964, Swieca–Völkel; the geometric identity check is left to the reader. The wedge theorem is cited, not proved.

### Fit and traps

- Scalar field only. Remark 5 (p7) records the formula for `Γ(U_ρ)` on the free electromagnetic field (attributed to Buchholz, private communication, and Gross 1964); nothing for helicity 2.
- `d = 3` (our case, odd): the chain is Buchholz → ray inversion → double cone. A theorem that wants no unread premise cannot rest on this file alone, since Buchholz 1978 is not on file.
- Sign convention: `Δ^{−iλ}_{O_1} = Γ[V(2πλ)]`, same as Longo–Morsella Thm 4.1; note that Hislop–Longo's `v_λ` uses `e^{λ}` where Longo–Morsella's `Z(z,s)` uses `e^{−s}`, with `f_λ = γ·f∘v_{−λ}` acting on test functions against `Φ∘Z(·,s)` acting on waves. The two are consistent only after that transposition; do not copy a flow direction from one paper into the other's formula.

---

## What I did not read

- CLR §3.1 (pp12–14), the wedge entropy computation.
- Garbarz–Palau outside pp24–35; its Appendix D (generic position of `S_R(B)`, `S_I(B)`).
- Camassa beyond pp2–4; BGL beyond pp8–10 (in particular its Thm 2.5, wedge duality).
- Reible: chapters II, III and the appendices, except the statements of Def II.2.10 and the skim of §III.4; nothing after §IV.4's first remark.
- Araki 1976 after p3 (all proofs); Araki 1977 pp1–4 as text only (garbled), §4 (proof of Thm 3.7) and later.
- Hislop–Longo Appendix 1 formulas (illegible in the text layer, not viewed as images).
- Not on file at all: Araki 1963, Araki 1964, Eckmann–Osterwalder 1973, Leyland–Roberts–Testard 1978, Buchholz 1978, Bisognano–Wichmann 1975/76, Uhlmann 1977, Ohya–Petz, Longo's "Entropy of coherent excitations" (CLR ref [19]).
