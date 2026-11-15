# ODD-DICT: the odd change of variables and coupling, every l

## Lead's check and grade

**Grade: 388 of 420.** Accepted.

**Checked by the lead.** (1) Re-ran `C:\Python314\python.exe odd_dict_check.py <out.json>` from `scripts\`: exit 0, "FAILED LINES: none", all five plants caught at l = 4 as the lane's table states. (2) By hand: r = tanh(ρ/2) gives dr/dρ = (1 − r²)/2 = N, so ½∫f²r²dr = ½∫q²dr = ½∫N q² dρ with q = r f: step 2 is a one-line substitution with no boundary term, coupling 1. (3) Step 3's reading g = w − 1 agrees with the even dictionary's structure (the free generator is weight 1; weight 2 adds one copy of N q).

**Not checked by the lead.** Step 5's bilinear identity beyond the lane's own exact rows; it rests on CHAIN-ALL-L and CHAIN-REF for general l, as the lane says. Step 6 is inherited from EVEN-DICT.

**What it changes.** The odd tower's "g = 1, D₀(q, π) = (0, −Nq), block diagonal" is no longer inherited from HINGE-PROOF: it is proved for every l, given the sector identity A_ext = c_l E₂ (card Q02, still running). No constant in `odd-sign` or `odd-remainder` moves. The extra term sits on extrinsic-curvature data in both parities.

**Process note.** Pass lines were fixed before the run; one run; no fix round. The lane's first plant is the one the card asked for.

---

## Lane report (unchanged)

**ODD-DICT report (Q08, 16 Nov 2026)**

**1. Verdict.** The odd change of variables holds as an exact identity for every l ≥ 2, with g = 1 and D₀(q, π) = (0, −Nq). No constant moves, so nothing changes in `odd-sign` or `odd-remainder`. All 39 pre-fixed lines passed on the first and only run, and all 15 plant instances were caught.

**2. Theorem and proof.** Odd sector l, ψ_B = f(r)Y, ψ̇_B = p(r)Y, integrals in units of N_Y. Put r = tanh(ρ/2), N = (1−r²)/2 = 1/(1+cosh ρ), q = r·f, π = N·r·p, λ² = l(l+1). Input, not proved here: A_ext = c_l·E₂ on odd data (Q02).

1. [P] **(a)** E₁ is the same functional of (ψ, ψ̇) as in the even sector, so EVEN-DICT steps 1 to 3 apply verbatim [I, lead-checked]: c_l·E₁ = c_l·N_Y·[½⟨π, π⟩ + ½⟨q, h²q⟩], and c_l∫(ψ̇_aψ_b − ψ̇_bψ_a) = c_l·N_Y·∫(π_a q_b − π_b q_a)dρ. Re-verified here from both sides' own definitions [C].
2. [P] **(b)** Z = E₂ − E₁ = ½∫_B ψ². Then ½∫f²r²dr = ½∫q²dr = ½∫N q²dρ, using dr = N dρ. So Z_odd = c_l·N_Y·½⟨q, Nq⟩, g = 1, with no integration by parts and no boundary term. It converges because N ~ 2e^(−ρ) while q → f(1).
3. [P] **(c)** On the ball, X = (Nψ̇, ∇·(N∇ψ) − 2ψ) satisfies Ω(b, Xa) = dE₂(a)(b); the boundary term vanishes because N = 0 at the edge. With q_ρρ = N(N q_r)_r, one line gives N·r·[(r²N f′)′/r² − Nλ²f/r² − w·f] = q_ρρ − λ²q/sinh²ρ − (w−1)N·q. At w = 2 the field is (π, −h²q − Nq) = X₀ + D₀. The constant c_l·N_Y multiplies Hamiltonian and symplectic form alike and cancels. g = w − 1 = 1 holds whichever sign convention is used for ω, since X₀ and D₀ come from one Hamiltonian.
4. [P] Block diagonal: Z = ½∫ψ² and the Y_lm are orthogonal.
5. [P]/[C] **(d)** ψ_B = −x·curl K·x is linear in K alone; ψ̇_B = x·curl(−½∇²h)·x in h alone. So q carries extrinsic curvature and π carries the metric, which is the swap of B27. On TT tensors S[η] (CHAIN-ALL-L's family, re-verified TT): ψ_B = l(l+1)η·φ, ψ̇_B = ½l(l+1)(η″ + (2l+2)η′/r)·φ, and ∫_B K_a·h_b = −c_l∫_B ψ_B[K_a]ψ̇_B[h_b], hence Ω_grav = +c_l·Ω_scalar. Z is a form on K alone, and its flow D₀ moves only h.
6. [I, not checked here] By EVEN-DICT, the even Z sits on π, which there is also the K half. So Z sits on extrinsic-curvature data in both parities.

**3. Check** (`odd_dict_check.py`, 13 s, exit 0). Five f data and five p data (r^(l+2k), k = 0..3, plus a mix). The r side is exact rationals computed twice (3D Cartesian integrals with a harmonic polynomial, and the radial reduction; they agree exactly). The ρ side is Gauss–Legendre at 45 digits.

| l | (a1) momentum | (a2) field | (a3) symplectic, 25 entries | (b) Z, g = 1 | 10×10 E₂ Gram | g from pointwise flow, max \|g − 1\| | ∫K·h / ∫ψ_Bψ̇_B |
|---|---|---|---|---|---|---|---|
| 2 | 1.4e-46 | 1.2e-46 | 1.8e-46 | 2.3e-46 | 1.9e-46 | 4.5e-43 | −1/6 |
| 3 | 1.1e-46 | 1.4e-46 | 1.2e-46 | 1.8e-46 | 1.8e-46 | 5.5e-43 | −1/30 |
| 4 | 1.5e-46 | 1.6e-46 | 1.5e-46 | 2.4e-46 | 1.4e-46 | 7.4e-43 | −1/90 |

- **Exact lines:** L0, L8, L9, L11, L12 and L13 hold with difference 0.
- **Quadrature:** 60 against 80 nodes agree to 3e-46.
- **Plants, caught at every l:**
  - P1 (the card's): Z with N² deviates 0.83 to 0.92.
  - P2: π = N²rp deviates 0.86 to 0.98.
  - P3: q = Nrf deviates 0.91 to 0.99.
  - P4: weight 3 reads back g = 2.0.
  - P5: +c_l in the bilinear gives a non-zero difference.
- **Failed lines:** none. No fix round was used.

**4. Re-run.** From the lane's `scripts` folder:
`C:\Python314\python.exe odd_dict_check.py ..\out\odd_dict_check.json`

**5. Weakest step.** Step 5's bilinear identity is exact only at l = 2, 3, 4 on three profiles vanishing to third order at the edge; for all l it rests on CHAIN-ALL-L and CHAIN-REF [I, refereed]. The thermal norm on (q, π) is untested here (MODGEN, Hislop–Longo).

**6. Bulletin.** Read B1 to B27; no steer. B27's sign was taken as the pre-fixed line L13 and came out as stated.

**7. Sentences.**
- **May use:** "In the variables q = rψ_B, π = N r ψ̇_B, r = tanh(ρ/2), the odd charge c_l E₂ is, for every l, the free generator plus ½⟨q, Nq⟩ with coupling exactly 1, block diagonal per sector; q is built from the extrinsic curvature alone, so the extra term sits on extrinsic-curvature data."
- **Must not use:** "The odd dictionary proves the charge is c_l E₂" (that is Q02's input), nor "the metric half carries the extra term".

**Files written**, all in the lane's folder `odd-dict\`:
- `logs\00_passlines.log`
- `logs\01_odd_dict_run1.log`
- `scripts\odd_dict_check.py`
- `out\odd_dict_check.json`
