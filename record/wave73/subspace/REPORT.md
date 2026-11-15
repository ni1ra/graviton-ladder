# REPORT: lane SUBSPACE (wave 73)

## Lead's check and grade (2026-11-15)

- Re-run by the lead: `python subspace.py 4 3` gave 27 of 27 lines as expected in 4.9 s (21 PASS, 6 plants CAUGHT, none missed), including A3 (kappa = -180 = -2/c_4) and both blocks of B2 on the whole Gram matrix. This is the lane's own code re-executed, not an independent implementation. `hl_gamma.py` and the other l values were not re-run by the lead.
- Correction to the lead's brief, raised by the lane: this input was not "never rechecked"; `wave70/modflow/REPORT.md` holds an independent proof. Neither that proof nor this one has had a second reader.
- Standing after this lane: inside the theorem, what remains is general l for c_l (record's symbolic computation), the standard second-quantisation step (cited, unread), and the choice of the curvature net as the ball's algebra. Outside it, step 12 (that this is the charge JSS mean) is the largest open gap; the even tower is not covered.
- Lane's own disclosures worth keeping in view: TT gauge is the wrong representative for data reaching the edge (gauge-dependent boundary term); A5 was added post hoc.
- Grade: 372 of 420. Exact arithmetic with a plant that isolates the complex structure, first-hand reading of Hislop-Longo and Longo-Morsella, repaired a thin step in the record, disclosed its own script errors; marked down because general l is a pattern at l = 2 to 6.

The lane's report follows unchanged.

**SUBSPACE report (2026-11-15)**

Tags: [P] proved this session, [N] numeric or exact-arithmetic this session, [S1] checked against a source read first-hand, [S2] source quoted second-hand or from memory, [I] inherited unchecked, [B] belief.

## 1. Verdict

**The input holds; I found no gap inside it.** On each odd sector the transfer (h, K) → (ψ, ψ̇) is a symplectic bijection with constant c_l that is local in r in both directions, it carries the graviton vacuum's two-point function onto c_l times the massless scalar's, and for that scalar Hislop–Longo's theorem gives the ball's modular generator as 2π·E₁. Hence K_Ω = 2π·c_l·E₁[ψ]. All of this is exact at l = 2 to 6 on polynomial and Gaussian profiles, with controls passing and every plant caught.

What stays conditional sits outside the three parts:
- the value of c_l at general l rests on the record's symbolic computation;
- the step from the one-particle subspace to the von Neumann algebra is the standard second-quantisation theorem, cited and not read;
- "the ball's algebra" means the curvature net.

One correction to the brief: this input was not "never rechecked". `wave70/modflow/REPORT.md` holds an independent proof that the ball's modular group in the helicity-2 curvature net is the conformal flow, using neither Hislop 1988 nor the transfer. Its S3 measured the sector two-point function at 1/c_l = 6, 30, 90 for l = 2, 3, 4. I read only its sections 1 and 4, so I am reporting that it exists, not vouching for it. Nobody has checked that proof either.

## 2. The record's Theorem 1

"Paper A" is `wave36/reviseA18/NOTE.md` (the companions cite revision 14). Theorem 1 is at line 57, Corollary 1 at line 109.

**Statement.** For a real sector σ = (l, parity), l ≥ 2, with V_σ = Φ/N_l and N_l² = 4πG(l+2)!/(l−2)!:
- (a) V_σ extends to a unitary from the helicity-2 sector onto the scalar l-sector, intertwining time translations, rotations and dilations.
- (b) V_σ(H(B) ∩ H_σ) = H_s(B) ∩ H_{s,l}, and H(B) is the direct sum of its sector parts. H(B) is the closure of the vectors of Weyl curvature smeared in D(B).
- (c) The coherent-state relative entropy is Σ c_l·Q[u].
- Corollary 1: H(B) is standard and factorial, its modular group is the unitary conformal flow along ζ, and A = 2π·δ_ζ on smooth data compactly supported in B.

**What its proof covers** [re-read line by line this session, not machine-checked]:
- Isometry comes from Proposition 1 (Ω = c_l·Φ*Ω_chain). Onto comes from Schur's lemma on the momentum fibres.
- The inclusion V(H(B)) ⊂ H_s(B) holds because curvature data of e_F are compactly supported in B.
- The reverse inclusion uses symplectic complements: V(H(B)′) ⊂ H_s(B)′ because ψ vanishes where the curvature does, then double complement. No Haag duality is used.
- I found each step complete. One step is thin: Step 1 says V commutes with the complex structures because it commutes with the time-translation generators, on a domain where V is not yet known bounded. The repair is in section 3(b).

**What it does not cover:**
- the constant c_l at general l by anything other than its own symbolic computation;
- the passage from H(B) to the algebra (Lemma 2 cites Araki);
- whether the curvature net is the right ball algebra (Proposition 2 cites Hislop 1988 for duality, which nobody here has read).

## 3. Verdicts for (a), (b), (c)

Script: `scripts/subspace.py L 3`, exact rational arithmetic. All 27 lines came out as expected at each of l = 2, 3, 4, 5, and at l = 6 with KMAX = 2. Conventions are charge-id's.

### (a) Symplectic isomorphism onto, localisation both ways: PROVED at l = 2 to 6 [P/N]; general l rests on the record [I]

| Line (pre-registered unless marked) | Result |
|---|---|
| A1: ψ[build_T(s^k)] = l(l+1)·s^k·φ, and build_T is TT | PASS |
| A2: TT odd tensors a(s)K_l + b(s)P of degree ≤ 3 form a space of dimension exactly 4 | PASS (8 unknowns, rank 4) |
| A3: ψ̇[α(r)·P] = κ_l·α·φ with one constant | PASS, κ_l = −12, −60, −180, −420, −840 = −2/c_l |
| A4a/b: ∫_B K_a·h_b = −c_l∫_B ψ_a ψ̇_b for data vanishing to second order at the edge, h in TT and in Regge–Wheeler gauge | PASS |
| A5 (post hoc): sphere-averaged symplectic densities agree as polynomials in r, Regge–Wheeler gauge, arbitrary profiles | PASS |
| Plant PA1: (l+1) in place of (l+2) in the inverse | CAUGHT (not transverse) |

What this proves:
- **K ↔ ψ is a local bijection.** ψ is algebraic in K's profile η, and K is recovered from η and η′ with no integration: g = −η′/((l+2)r), f = 2(η − r²g)/(l−1). So supp K = supp ψ on every shell.
- **h ↔ ψ̇ is algebraic in Regge–Wheeler gauge.** With h = α(r)(x⊗v + v⊗x), ψ̇ = −(2/c_l)·α·r^l Y. This matches Benedetti–Casini's eq. (5.24), where the scalar momentum is a constant times h₁ [S1].
- **Onto, as a statement about the ball.** Every smooth (ψ, ψ̇) supported in a closed shell has a preimage supported in the same shell. At Hilbert-space level, onto is Paper A's Steps 2 to 5.
- **Data touching the boundary.** H(B) is a closure, so nothing special happens to vectors. In Regge–Wheeler gauge the symplectic identity holds sphere by sphere, with no boundary term (A5). In TT gauge it fails for data reaching the edge: the ratio Ω_B/(c_l·Ω_s,B) came out −2, −20/7, −13/9 and similar at l = 2. That is a gauge-dependent boundary term; the TT representative is the wrong one for edge questions. It does not affect charge-id's value identities.
- Polynomial data cut off at r = 1 are not vectors of H(B) (a jump is not in H^{1/2}). They enter only through E₁, E₂ as forms and through Paper A's compact completion.

### (b) Vacuum two-point function: PROVED given (a) [P]; confirmed exactly at l = 2 to 6 [N]

**Proof.** ψ(t, x) = x·B(t, x)·x is pointwise in t with real, time-independent coefficients. So the positive-frequency part of ψ[k] is ψ of the positive-frequency part of k, and Φ∘J = J_s∘Φ on solutions with Schwartz data. With Ω = c_l·Φ*Ω_s this gives μ(a, b) = Ω(a, Jb) = c_l·Ω_s(Φa, J_sΦb) = c_l·μ_s(Φa, Φb). The vacuum is the quasi-free state fixed by μ, so the two-point functions correspond. ∎ This also repairs Step 1 of Paper A.

**Test B2** (pre-registered line: equality on the whole Gram matrix, both blocks, non-zero), on R³ with profiles q(s)e^{−s}:
- K–K block: ∫2K_a·K_b/ω = c_l∫ω ψ_aψ_b. PASS; for example 24/5 = 24/5 at l = 2.
- h–h block: ∫½ω h_a·h_b = c_l∫ψ̇_aψ̇_b/ω. PASS; 48/5 = 48/5 at l = 2.

**Controls, all pass:**
- B0: Parseval, position space against momentum space.
- B1a: Ω_grav = κ·c_l·Ω_s with κ = 1 (charge-id's T4).
- B1b: ∫K_a·K_b = c_l·½∫∇ψ_a·∇ψ_b (Paper A, Proposition 1).

**Plants, all caught:**
- PB1: c_l scaled by (l+3)/(l+2).
- PB2: the map (h, K) → (2h, K/2), which keeps the Ω identity and breaks B2. This shows the test sees the complex structure, not only the symplectic form.
- PB3: ω and 1/ω exchanged.

### (c) Modular generator: STANDARD [S1], hypotheses checked; weight confirmed three ways [P/N]

- **The field.** ψ obeys □ψ = 0 in 3+1 flat space (charge-id T5) with the scalar vacuum two-point function by (b). A free massless scalar in flat 3+1 is the dimension-1 conformal field; there is no separate "conformally coupled" case because R = 0.
- **Hislop–Longo 1982, Theorem 2** (CMP 84, p. 76; `sources/` PDF, pages 71–76 read). Free massless scalar, three space dimensions treated explicitly, unit double cone at the origin: Δ^{−iλ} = Γ[V(2πλ)], with the multiplier γ quoted in `hl_gamma.py`.
- **Longo–Morsella 2012.00565** [S1, text extraction of pp. 12–16]. Theorem 4.1: Δ_B^{−is} = V(2πs). Proposition 4.4: K(f, g) = (½(1−r²)g, ½(1−r²)∇²f − r∂_r f − Df), D = (d−1)/2 = 1. Proposition 4.5: (Φ, AΦ) = ½∫(1−r²)T₀₀ + (D/2)∫Φ². Theorem 4.7: log Δ_B = −2πA. That is 2π·E₁ with w = 1.
- **C2** (`hl_gamma.py`, sympy): dγ/dλ at 0 is −3x₀; the generator of v_λ is −ζ; γ = (conformal factor)^{4−Δ} with Δ = 1; the induced field generator is ζ·∂φ + ¼(∂·ζ)φ. PASS. Plants Δ = 2, 3: CAUGHT.
- **C1** (exact, ball polynomials): Paper A's flow (1) on (E, B) gives x·δB·x = Nψ̇ and x·(sym curl δE)·x = ∇·(N∇ψ) − ψ, and stays TT. PASS. Plant w = 2: CAUGHT.
- **C3** (exact, scalar Gaussian data): δ₁ is skew for the vacuum form μ_s. PASS. Controls pass: time translation and dilation are skew for μ_s; δ_w is skew for Ω_s at w = 1 and 2. Plant w = 2: CAUGHT, with skew₂[0][0] = −2/15 at l = 2. So the weight-2 flow is symplectic but does not preserve the vacuum.
- **Sector restriction.** Rotations and inversion preserve D(B), so the sector projections commute with the modular data.

**What LATTICE-L already tests.** It tests (c) only, on the reduced chain, and only the field block on one static profile per l (to 10⁻⁷). wave70's S4 tests both blocks on three bump profiles at l = 2, 3, to under 1%. Neither touches (a). (b) on the tensor side is wave70's S3 and this lane's B2.

**2π bookkeeping** [P, by hand]: A/8G generates 2π·δʷ with respect to Ω/16πG, and c_l/16πG cancels in the flows. So the two generators are 2πE₁ and 2πE₂ in the same variables, and Z = ½∫ψ².

## 4. The flagship sentence

**May say** (charge-id section 4, with the conditional replaced):

"On every odd-parity graviton sector of a flat ball, the second-order area-side charge (fixed-volume, and equally the canonical charge of the conformal Killing flow with the background rescaling removed) is exactly c_l·E₂ of the scalar chain field, and its Hamiltonian flow is the weight-2 scalar flow. The vacuum modular generator of the ball's curvature algebra is 2π·c_l·E₁ in the same variables: the transfer is a symplectic bijection, local shell by shell, that carries the graviton vacuum onto the scalar vacuum, and Hislop–Longo's theorem then applies. Both statements are proved in the record for all l and checked independently in exact arithmetic at l = 2 to 6. Hence the area-side flow is not the modular flow of any normal state or weight in the vacuum sector: the Shale sum diverges logarithmically with per-sector constant 7ζ(3)/(16π²), confirmed on an independent lattice. It is the spacelike-slice counterpart of Satishchandran–Sorce's null-slice theorem, which left gravitons open."

If asked what is assumed, add: "the ball algebra is the one generated by the linearised curvature, and the odd tower only."

**Must not say** (unchanged): "We have proved that the Jensen–Sorce–Speranza assumption fails for gravitons, or that no equilibrium state exists for a ball." Also not: "this is the charge JSS define."

## 5. Unchecked, ranked, and the next experiment

1. Step 12, that JSS's charge is this one. Untouched here, and still the largest gap in the flagship.
2. General l for c_l, and for A1, A3, A5. The pattern κ_l = −2/c_l is observed at l = 2 to 6 only [B for all l].
3. Second quantisation (one-particle modular group to the algebra) and the choice of the curvature net as the ball's algebra [S2, from memory; the record's wave69].
4. The Hilbert-space "onto" argument was re-read by me, not checked by a second reader. wave70's independent proof is likewise unchecked by anyone else.
5. Hislop 1988 remains unread [S2]. Nothing in this chain needs it.
6. Even parity and d ≥ 5: not examined.

**Next experiment for this input:** prove A1, A3 and A5 with l as a symbol in Regge–Wheeler gauge. That turns the transfer into a short derivation that can be cited against Benedetti–Casini (5.24)–(5.25) and closes item 2.

## Disclosures

- `hl_gamma.py` failed two lines on its first run through my own script errors: sympy did not reduce a true zero under a fractional power, and I had the sign of dγ wrong in one scripted formula. The pre-registered statement did not change. First run kept in `logs/hl_gamma_first_run.log`.
- A5 is post hoc, added after the first run's INFO line showed ratio 1 in Regge–Wheeler gauge. First-run logs are in `logs/run1/`.
- One exploratory shell command of mine (a disk-wide `find`) went to the background; I stopped it with the task tool. No process was killed by pid.
- Native Windows Python 3.14 only; no WSL, no subagents, no web requests (all sources were on file). Nothing committed; no report file written.

## Reproduce

From `record\wave73\subspace\scripts`, each run 2 to 7 s:
- `python subspace.py 2 3`, then `3 3`, `4 3`, `5 3`, `6 2`
- `python hl_gamma.py`

Files added, all in `record\wave73\subspace\`:
- `PREREG.md`
- `scripts\subspace.py`
- `scripts\hl_gamma.py`
- `logs\subspace_l2.log`, `subspace_l3.log`, `subspace_l4.log`, `subspace_l5.log`, `subspace_l6_k2.log`
- `logs\run1\subspace_l2.log` to `subspace_l5.log`
- `logs\hl_gamma.log`, `logs\hl_gamma_first_run.log`
- `out\subspace_l2.json` to `subspace_l6.json`
