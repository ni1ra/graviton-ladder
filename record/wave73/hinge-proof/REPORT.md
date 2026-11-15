# REPORT: lane HINGE-PROOF (wave 73)

## Lead's check and grade (2026-11-15)

- Lemma A re-checked by the lead by hand, line by line: coth(y) <= 1/y + 1; x coth(pi x) >= max(x, 1/pi); N <= h/lambda is equivalent to N h^-1 N <= h/lambda^2; the constant 2 g^2 / lambda^2 follows. It holds.
- Not re-checked by the lead: the Duhamel-Gronwall corollary, Lemma B (an outline in the lane's own words), the numerical test (not re-run), and every cited source (quoted from memory or a tool summary, per the lane).
- Standing of the flagship claim after this lane: the mathematical core holds given the reduction to scalar sectors; three physics gaps remain (steps 1, 9, 12). The honest public sentence is the lane's "may say" sentence in Section 5.
- Grade: 360 of 420. Closed the named gap with a proof, reported its own failed check, searched for counterexamples; marked down for sources not read first-hand and Lemma B left as an outline.

The lane's report follows unchanged.

**HINGE-PROOF report (2026-11-15)**

Tags: [P] proved this session, [S] checked against a cited source, [N] numeric this session, [I] inherited unchecked, [B] belief.

## 1. Verdict

The non-existence claim **stands with named gaps, and the gap HINGE named is not one of them any more**. The mathematical core holds: given the reduction to scalar sectors, the flow of K_Ω + 2πZ is the modular flow of no faithful normal semifinite weight, hence of no normal state. The uniform remainder beyond first order is now controlled by an operator-norm lemma I proved, and my own exact code confirms it. What remains open is on the physics side:

- the reduction itself is a record input I did not re-prove;
- the charge is only defined on the odd-parity tower;
- nobody has shown this charge is the one Jensen–Sorce–Speranza (JSS) mean.

I found no solvable case where the reasoning gives a wrong answer.

## 2. Step table

| # | Step | Verdict | Reason |
|---|---|---|---|
| 1 | Reduction: odd (l, m) graviton sector is a scalar l-wave sector; charge = 2πE₂ | GAP [I] | `wave69/final/THEOREM.md` proves it only for stationary magnetic data. The extension to all odd data is "from the record, not re-proved here" (line 113), and the area side rests on hypothesis H4 (line 111). |
| 2 | Frame identification: ball vacuum = β = 2π thermal state on ℝ × H³, Z = potential gN(ρ) | STANDARD + [P] + [N] | Hislop–Longo 1982 and Casini–Huerta–Myers 2011, quoted from memory. I re-derived the radial operator h_l² = −d²/dρ² + λ²/sinh²ρ. The blind lattice in the flat frame gives the same constant. |
| 3 | The flow maps H to H boundedly | PROVED [P] | Not addressed by HINGE. Lemma A below. |
| 4 | Criterion: weight exists ⇒ S_t ⊕ 1 bounded with Hilbert–Schmidt antilinear part | PROVED [P], given Connes 1973 and Shale 1962 (memory) | Weights are covered because Connes's theorem is for weights. One line is missing in HINGE: Shale needs S⁻¹ bounded, which follows from u_t* by the same argument. Satishchandran–Sorce Section 5 uses the same route for weights [S, via tool summary]. |
| 5 | Sector formula with thermofield weights | NUMERIC [N] | My independent exact code reproduces HINGE's table at g = +1: 0.2443 against 0.244 (l = 16), 0.2777 against 0.278 (l = 64). |
| 6 | Limit constant positive; a₁ = 7ζ(3)/(16π²) | NUMERIC [I + N] | hinge-blind for a₁. My λ²·E1_l rises 0.230, 0.255, 0.270, 0.280, 0.285, 0.288 towards J(0.5)/2 = 0.2917. |
| 7 | Uniformity in l beyond first order | PROVED for the remainder [P]; first-order limit by standard theorems, outline only [P-outline] | Lemmas A and B below, plus the numeric test. |
| 8 | Sum over sectors diverges | PROVED [P] | Blocks are orthogonal and invariant; multiplicity 2l+1 against c/λ². A log divergence suffices because the criterion is "finite or not". |
| 9 | Even/electric tower "can only add" | GAP | True if the charge exists there and is block diagonal. The record says the area side is gauge dependent on electric data. The theorem is about the odd tower. |
| 10 | Not an artefact of reference or gauge | PROVED for reference and edge terms [P]; GAP for the choice of charge | See Section 3d. |
| 11 | State corollary via product states | PROVED in outline [P], redundant with 4 | Tail-norm convergence plus multiplicativity of fidelity; no semifiniteness needed in this direction. |
| 12 | This charge is JSS's H^g_ξ, so their assumption A4 fails | GAP [B] | Not written anywhere in print; see Section 4e. |

Not examined: HINGE Sections 6 and 7 (lapse tuning, relative entropy, the replacement functional) and the exponents.

## 3. Proofs and the numerical test

**Setting.** Sector l, u = sinh ρ × radial function, h² = −d²/dρ² + λ²/sinh²ρ on L²(0, ∞), N = 1/(1 + cosh ρ). Data k = (q, π) in H have norm ‖k‖² = ½⟨q, h coth(πh) q⟩ + ½⟨π, h⁻¹coth(πh) π⟩. The perturbed generator is 2π(X₀ + D), with D(q, π) = (0, −gNq).

**Lemma A [P].** ‖2πD‖ ≤ 2√2·π|g|/λ on H, for every l ≥ 1.

1. N²sinh²ρ = tanh²(ρ/2) < 1, so N² ≤ 1/sinh²ρ ≤ λ⁻²h² as forms, since −d²/dρ² ≥ 0.
2. Hence ‖h⁻¹N‖ ≤ 1/λ, and by Löwner–Heinz N ≤ h/λ, so ⟨Nq, h⁻¹Nq⟩ ≤ λ⁻²⟨q, hq⟩.
3. coth(πx)/x ≤ 1/(πx²) + 1/x, so ⟨Nq, h⁻¹coth(πh)Nq⟩ ≤ λ⁻²(‖q‖²/π + ⟨q, hq⟩).
4. x·coth(πx) ≥ max(x, 1/π), so both terms are at most 2‖k‖². Therefore ‖Dk‖² ≤ 2g²λ⁻²‖k‖². ∎

**Corollary [P].** The vacuum flow is isometric on H, so by Duhamel and Gronwall:

- T_t is a bounded symplectic group on H, with ‖S_t‖ ≤ e^{δt}, δ = 2√2·π|g|/λ;
- ‖S_t − 1 − B₁‖ ≤ e^{δt} − 1 − δt ≤ 4π²g²t²e^{δt}/λ², where B₁ is the first-order term.

**Why this closes the gap [P].** A Hilbert–Schmidt norm is bounded below by its value on one unit vector, so no Hilbert–Schmidt control of the remainder is needed.

- Take e = k(ψ) in H with ψ in the domain of coth^{3/2}(πh). Then i·e = k(coth(πh)·iψ) + k′(−csch(πh)·iψ) lies in H + H′.
- So βe = ½(Ŝe + iŜ(ie)) is explicit: the zeroth order cancels, the first order is (g/λ)·β̂₁e, and the rest is at most C(ψ)·g²/λ².
- If liminf ‖β̂₁^(l) e_l‖ = c > 0, then ‖β^(l)‖²_HS ≥ (cg/2λ)² for large l, and Σ(2l+1)/λ² diverges.

**Lemma B [P-outline, not written line by line].** The liminf is positive.

- With y = ρ − ln 2λ, the operators h_l² converge to −d²/dy² + e^{−2y} on the common core C_c^∞(ℝ), hence in strong resolvent sense (Reed–Simon VIII.25 in Weidmann's form).
- λN_l → e^{−y} uniformly on compacts, and ‖λN_l h_l⁻¹‖ ≤ 1 from Lemma A turns this into weak convergence of λN_l G(h_l)χ.
- So first-order matrix elements between fixed wave packets converge to those of the Liouville model. Its β̂₁ is non-zero because its squared norm is J(t)/2 > 0. That last number rests on HINGE's Bessel identity [I, numeric there].

**Numerical test [N].** `scripts/u1_uniformity.py`, own code, exact in g, n = 1250, t = 0.5, l = 4 to 128, g = ±1. Lines were fixed in the PREREG before the run.

| Check | Pre-registered line | Result |
|---|---|---|
| U1 | λ‖N h⁻¹‖ ≤ 1 | PASS: 0.633 rising to 0.948 |
| U2 | ‖S_t − 1‖ ≤ e^{δt} − 1 | PASS: 0.0087 against 0.0352 at l = 128 |
| U3 | λ²‖S_t − 1 − gB₁‖ ≤ 4π²g²t²e^{δt}·1.02 | PASS: 0.42 to 0.66 against 10.4 or more |
| U4 | deviation of exact from first order at l = 128 at most 0.25 × that at l = 16; λ²E1 not falling | PASS: ratios 0.136 and 0.128; deviation 0.54% at l = 128 |
| U5 | l² × mean over g = ±1 within 3% of HINGE's table | **FAIL as written**: 1.044 and 1.010. HINGE's β rows are g = +1, not the mean; I misread the header. Against g = +1 the ratios are 1.001 and 0.999. |
| C1 | lapse-weighted control falls by at least half from l = 16 to 64; sharp does not fall | PASS: 0.273 against 1.102 |

**3d. Reference and gauge.**

- Independent of the reference state [P]: all weights' modular groups differ by cocycles in M, so the vacuum is only a convenience.
- Edge terms do not help [P]: ∫_B F = ⅓(∮F − ∫N∇²F), so Z is an edge term plus a lapse-weighted dimension-4 term. The edge term does not change the flow on H, and the bulk term lies outside JSS's finite window.
- The choice of charge does matter: the divergence is proportional to g², and g measures how far the charge sits from the vacuum's own generator by an improvement-type term. Whether gravitons allow such a redefinition is physics outside the theorem.

## 4. Counterexample search

The argument would be broken if it gave non-existence in (c) or (d), or existence in (b). It did neither.

- **(a) Killing case:** Z ≡ 0, nothing fires [I, HINGE T5].
- **(b) Rescaled temperature:** non-existence is known exactly (modular automorphisms of a type III₁ factor are outer). The reasoning gives a power-law divergence, correct [P, sketch].
- **(c) Mass-type, lapse-weighted term in d = 4:** existence is the published expectation (Eckmann–Fröhlich 1974, from memory; JSS Section 4.3). My control gives sector terms falling as l⁻³, a convergent sum [N]. This is a consistency check, not a strict test, because that flow is not the massive vacuum's modular flow.
- **(d) 2d thermal interval:** a normal state whose lapse differs from the vacuum's at second order at the edge. The reasoning has one sector and a finite constant, so no false non-existence [B, by analogy].
- **(e) Satishchandran–Sorce, arXiv:2607.20678, "Uniqueness of null-local modular flow"** [S, tool summary of the HTML, not read first-hand]:
  - Theorem 5.1: for the massless scalar, D ≥ 3, the only weight in the vacuum sector with null-local modular flow is the vacuum. Same Connes-cocycle and Hilbert–Schmidt logic.
  - They build states with other null flows outside the vacuum sector.
  - Appendix B.2 builds a graviton state with geometric null flow, says it is not the vacuum's restriction, and does not assert uniqueness for gravitons.
  - The spacelike-slice case is called open.
  - Consequence: our claim agrees in direction and sits in the case they left open. But for a non-Killing ξ the null charge and the spacelike charge differ, so "no state in the vacuum folium" does not mean "no state".

Limit of the search: no exactly solvable normal state with an unweighted edge term is known to me, so the search is not decisive.

## 5. What may and may not be claimed

**May say:** "Granting our reduction of the odd-parity graviton sectors on a flat ball to scalar sectors, the flow of the area-side charge is not the modular flow of any normal state or weight in the vacuum sector: the Shale sum diverges logarithmically with per-sector constant 7ζ(3)/(16π²), confirmed on an independent lattice, and it is the spacelike-slice counterpart of Satishchandran–Sorce's null-slice theorem, which left gravitons open."

**Must not say:** "We have proved that the Jensen–Sorce–Speranza assumption fails for gravitons, or that no equilibrium state exists for a ball."

## 6. Unchecked, ranked, and the next experiment

1. That the area-side charge is JSS's H^g_ξ on the slice (step 12).
2. The record's reduction for all odd data, including momenta and the constant c_l (step 1).
3. The even tower (step 9).
4. Lemma B written line by line, and HINGE's Bessel identity proved rather than computed.
5. Satishchandran–Sorce read first-hand; Connes, Shale, Hislop–Longo and Eckmann–Fröhlich locators are from memory.
6. HINGE Sections 6 and 7 and the exponents: untouched.

**Next experiment:** compute the second-order canonical charge of linearised gravity for ξ = the conformal Killing field on the t = 0 slice, symbolically, on general odd-parity data (field and momentum) at l = 2, 3, 4, and test whether it equals c_l·E₂ in the transfer variables. It decides items 1 and 2 together.

Files are in `record\wave73\hinge-proof\`:
- `PREREG.md`
- `scripts\u1_uniformity.py`
- `logs\u1_run1.log`
- `out\u1_uniformity.json`

Nothing committed, no report file written.

Sources: [arXiv:2607.20678](https://arxiv.org/abs/2607.20678), [arXiv:2306.01837](https://arxiv.org/pdf/2306.01837)
