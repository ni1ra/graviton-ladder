# PREREG: w64-btheorem. Pillar B at d = 4 as a theorem

Lane w64-btheorem (Opus), Mon 9 Nov 2026, written 18:30–18:50 Oslo, before any derivation script ran (only `scripts/price.py` has run). Claude only.

## Read before writing (not blind)
- `wave45/truthcond/REPORT.md` §2.2–2.3 (clause (3) verbatim); `wave39/BR/` (PREREG, REPORT); `wave60/br6/REPORT.md`; `wave60/cprime/REPORT.md`; `wave61/cocycle/` (PREREG, REPORT); `wave61/edgelog/REPORT.md`; `wave51/b4/REPORT.md`; `wave52/b5/REPORT.md` and its `price_out.json`; `wave31/backreact/REPORT.md`; Paper A r8 (`wave36/reviseA8/NOTE.md`) §1–3, §8–9, App. D; CALIBRATION bookings w39-BR, w47-hforced, w49-hwforce, w51-b4, w52-b5, w60-br6, w60-cprime, w61-cocycle, w61-edgelog.
- JSS 2306.01837 (library PDF, extracted): §2.1 (A1–A6), §3.2, §3.3 (region fixing: (H_ξ^Σ)^(−1) = 0 is a fixed-volume gauge condition for the ball; A − kV proposal), §5.3 (eq. (5.35) first law, (5.36) S_gen), App. A eqs. (A.20)–(A.26).
- **The route below was conceived before this file was written**, while reading. The derivation scripts have not run. Its likelihoods are priced from structure.

## 1. The object (fixed exactly)

Clause (3) of the truth condition: "Every term of the second-order modular side that is linear in the second-order metric is the smeared constraint δ(A − kV). So the area side of the second-order first law is A."

- The free-net CLR relative entropy is an exact quadratic form in the first-order datum. On that object (3) is vacuous. So (3) is a claim about the **G⁰ gravitating** first law: the geometric (area) side, and any gravitational term on the modular side (JSS's area operator, observer energy, dressing, edge charge).
- **Theorem form ("total form").** Let 𝔄 be the total h⁽²⁾-linear part of the second-order first law, area side minus modular side, at O(ε²). Claim: 𝔄 = −(1/4G)·δ(A − 2V)[h⁽²⁾] = (1/8G)∫_B N R_lin[h⁽²⁾]. On the second-order constraint surface of the scope family (h⁽¹⁾ = 2tK) this is A/8G, i.e. A in Paper A's units.
- **Relation to (3) verbatim.** LADDER-PHYS uses (3) only through "the area side is A". The total form gives exactly that. Given clause (2) (C: the physical modular side equals the free relative entropy, hence has no h⁽²⁾), the total form is (3) verbatim. Without C, a term τΨ[h⁽²⁾] could sit on the modular side only balanced by −τΨ on the area side; the NOTE states this as a restatement, not a move.

## 2. The theorem (to be proved in NOTE.md)

**Setting.** d = 4, Minkowski, Σ = {t = 0}, unit ball B, S = ∂B, ζ the diamond's conformal Killing vector, N = (1 − r²)/2. First-order Cauchy data (q, p) = (δh_ij, δK_ij), smooth on a neighbourhood of B̄ (topology: C^∞(B̄), any linear functional; no continuity or locality is needed). δA[q] = ½∮_S(q_kk − q_nn), δV[q] = ½∫_B q_kk.

**Premises.**
- **(P0) Single Taylor functional at G⁰.** At O(ε²), h⁽²⁾ enters the first law only through 𝒟(h⁽²⁾), where 𝒟 is the linear functional by which the metric perturbation enters the same first law at O(ε). Sub-lemmas proved in the NOTE: (i) the free-field relative entropy is exactly quadratic in coherent data (no h⁽²⁾); (ii) a data-dependent region, frame or flow multiplies h⁽²⁾ only at O(ε³); (iii) matter-loop metric dependence is O(G) relative to the 1/G area term.
- **(T) Time reflection.** t → −t is an isometry of η preserving D(B) and Ω; the first-law functional is invariant (entropies are invariant under antiunitary symmetries). So 𝒟(q, p) = 𝒟(q, −p).
- **(U) Universality.** Adding a free conformally coupled scalar does not change 𝒟 on configurations where the scalar vanishes on Σ (its Wald term ∝ ∮φ² is then zero).
- **(FL1) First-order first law w.r.t. the vacuum, on the physical slice 𝒮.** For the mixtures ρ_λ = (1 − λ)Ω + λω_g, ω_g the coherent state of the scalar with data (φ, π) = (0, g), and every first-order metric datum (q, p) ∈ 𝒮 solving the linearised constraints with source λ⟨T⟩_{ω_g}: 𝒟(q, p) = λ·S(ω_g‖Ω) = 2πλ∫_B N·½g². (This is JSS's A5 first law (5.35) at first order and the BMSV first law of causal diamonds, read as a property of the unknown gravitating functional.) 𝒮 is either all data (relational region; c3) or any linear slice 𝒮_q × 𝒮_p that meets every gauge orbit and contains p = 0 (a physical edge rule, e.g. JSS's fixed-volume gauge).

**Claim.** 𝒟(q, p) = −(1/4G)δ(A − 2V)[q] on 𝒮, uniquely. Hence the total form of (3), and the area side on the scope family is A. The non-circularity point: BMSV's identity shows the JSS/Iyer–Wald functional satisfies FL1 (existence); the spanning lemma shows nothing else does (uniqueness).

## 3. Step list, each marked

| # | Step | Status before running |
|---|---|---|
| S1 | Green identity ∫_B N R_lin[q] = −2δ(A − 2V)[q], every smooth q (and the general-w form ∫wR_lin = ∫q_ij(∂_i∂_jw − δ_ij∇²w) − 2κδA when w\|_S = 0, ∂_rw\|_S = −κ) | proved in record (w31-backreact, w52-b5); re-verify here exactly, plant k = 2.2 |
| S2 | δ(A − 2V)[sym∇ξ] = 0 for every smooth ξ, including ξ ≠ 0 on S; δA, δV separately not invariant | to be proved here (exact, symbolic), with plants k ≠ 2 and a deformed region |
| S3 | Spanning lemma: the profiles ½g² (g smooth on a neighbourhood of B̄) span all smooth functions on B̄; for data (0, g) the improved T₀₀ is ½g², T₀ᵢ = 0, Wald term 0 | to be proved here (exact) |
| S4 | First law for mixtures: S(ρ_λ‖Ω) = o(λ) whenever S(ω‖Ω) < ∞ (Donald's compensation identity + lower semicontinuity) | to be proved here; inputs cited (Araki, Donald 1987, Ohya–Petz) |
| S5 | Coherent-state entropy S(ω_g‖Ω) = 2π∫N·½g² for data (0, g) (no Born/boundary term since φ = 0) | cited (Longo–Morsella; Paper A §1, Cor. 2(iii) structure) |
| S6 | Uniqueness: FL1 + S3 + T ⇒ 𝒟 = −(1/4G)δ(A − 2V) on 𝒮 | to be proved here |
| S7 | Taylor transfer to h⁽²⁾ (P0 sub-lemmas i–iii) | (i) cited (CLR / Paper A Thm 1); (ii) record (w51-b4 spine), re-argued; (iii) argued |
| S8 | Physical h⁽²⁾ lies in 𝒮 up to O(h⁽¹⁾²) corrections from the slice condition's quadratic part; these vanish on the scope family when the condition depends on q only (e.g. fixed volume) | to be proved here |
| S9 | Scope family: R_lin[q₂] = \|K\|² from the exact ADM Hamiltonian constraint at O(ε²) with h⁽¹⁾ = 2tK; ∫N\|K\|² = A = 1/35 at l = 2, unit flux | to be computed here (exact), gate vs Paper A |
| S10 | Power controls: each rival functional violates FL1 on an explicit configuration — angle-dependent area charge ∮Y₂₀q_AA, corner monopole αδA, weight N + (1 − r²)², corner jet ∮R_lin\|_S, T-odd term ∫K_nn, wrong normalisation | to be computed here (exact rationals) |
| S11 | Exact l = 0 classification closure (Route II cross-check): the gauge maps ξ = ra/2, f′ = rc remove a, c for **arbitrary** profiles; β = b − ra′, γ = d − rc′; R_lin = (2/r²)(rβ)′, n·div π = 2γ/r; regularity inverts both | to be proved here with sympy `Function` profiles (no polynomial ansatz) |

**The four open items the lead named:**
1. Converse-cocycle sufficiency at d = 4: **bypassed, prediction.** No step uses Ψ_ζ or JSS's A4 for gravitons. The weight comes from matter, where A4 holds exactly (HL: Ψ = Ω), and is carried to h⁽²⁾ by P0. w61's cocycle bears on the algebraic interpretation (C), not on (3).
2. (c′) at d = 4: **bypassed, prediction.** FL1 needs only the matter vacuum's HL generator, not graviton BW or Takesaki on the curvature net. No equivalent premise enters: FL1 is a first-order statement about matter perturbations.
3. Corner gauge (c3): **to be shown harmless for the linear clause** (S2, S6, S8): FL1 on any slice meeting every gauge orbit pins 𝒟 there; residual dependence is quadratic in first-order data (C territory). Failable test: S2 (consistency of FL1 across gauge orbits needs δ(A − 2V) to annihilate every corner displacement) and the planted αδA in S10.
4. Truncation closure: **bypassed** (S6 pins 𝒟 on all data, any jet order, non-local included) **and closed independently** (S11).

## 4. What counts as failure
- **A step false:** S1 or S2 fails for some smooth q or ξ; S3's span misses a profile; T₀₀ for data (0, g) is not ½g² or a Wald/boundary term survives; S9 gives an area side ≠ A; S11's gauge map fails for some profile. Apply the false-kill rule first (check conventions, normalisations, which object), then report plainly.
- **A step that cannot be closed:** S4 or S5 cannot be supported by a citable theorem; or P0, T, U or FL1 turns out to be non-standard in a way a referee would reject (e.g. FL1 contradicted by a printed first law). Name it as open.
- **The route is circular:** if FL1 can only be stated by presupposing 𝒟 = −(1/4G)δ(A − 2V), the theorem is a restatement (LR 1).

## 5. Pricing (`scripts/price.py`, `out/price_out.json`, run before any derivation)
- Slices of ¬B = 0.128 (w52-b5's table through b5's PASS rates): c3 0.0711, b4_q 0.0254, b4_sc 0.0167, b5_ref 0.0081, b5_bw 0.0061, b5_geo 0.0006.
- A blind referee confirms the implication, not the premises. P(confirm | slice s) = q_c(1 − p_prem[s]) + f·p_prem[s], with p_prem = P(P0 ∧ T ∧ U ∧ FL1 | s): c3 0.60, b4_sc 0.70, b4_q 0.40, b5_ref 0.10, b5_bw 0.80, b5_geo 0.80 (reasons in the script). P(prem | ¬B) = 0.552.

| Outcome | P(·∣B) | P(·∣¬B) | LR | B → | LADDER-WHY → |
|---|---|---|---|---|---|
| complete proof, pending blind check (f = 0.30) | 0.90 | 0.569 | **1.58** | 0.915 | 39.26% |
| complete proof, blind-confirmed (f = 0.08) | 0.70 | 0.358 | **1.96** | 0.930 | 39.91% |
| a step open (named) | — | — | 1.00 | — | — |
| a step false (after false-kill rule) | — | — | ≤ 0.80 | — | — |

- Caps: removing all of c3 gives 2.25; without c3 the cap is 1.22. **c3 is the lever.** b4_q and b5_ref are definitional and keep weight in ¬B.
- Double-count guard: c3 is also inside C (w39-BR flag; the lead booked no restatement). Only the part of c3 that sits in B's linear clause moves; the quadratic edge remnant stays in C. b5_bw is not removed inside B (it is the separately held graviton-BW premise); its p_prem is high only because the route does not use it.
- **Guard: would a wrong (BR) also yield a "complete proof"?** Yes, with probability about 0.36 after a blind check: whenever ¬B lives in ¬FL1 or ¬P0, the implication still holds.

## 6. Eureka
None at 20 points. The value is a referee-checkable reduction of B to the first-order first law.
