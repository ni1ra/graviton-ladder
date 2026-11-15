# REPORT: lane CHARGE-ID (wave 73)

## Lead's check and grade (2026-11-15)

- Re-run by the lead: `python charge_id.py 3 3` gave 32 of 32 checks passing in 3.4 s, including T1a/T1b (difference 0 on all Gram entries), the flow identities T2a/T2b, the stationary control C1 (E1/E2 = 4/5) and all planted errors caught. This is the lane's own code re-executed, not an independent implementation. Other l values and the record re-runs were not re-run by the lead.
- Standing after this lane: referee step 1 (the reduction, at the level of the charge and its flow) is closed for l = 2 to 6 in exact arithmetic and rests on the record's symbolic proof for general l. Step 12 is sharper but still open: the literal JSS formula (3.28) gives A + (1/2) int |h|^2, which is not gauge invariant; the fixed-volume reading gives c_l E_2. Step 9 (even tower) is now known to be unsupported.
- New remaining grant: the standard-subspace transfer with K_Omega = 2 pi E_1 (Paper A Theorem 1), inherited and unchecked.
- Grade: 375 of 420. Exact arithmetic, controls and plants, own derivation of the Lie derivative, first-hand reading of JSS, caught its own vacuous plant; marked down because general l leans on the record's scripts.

The lane's report follows unchanged.

**CHARGE-ID report (2026-11-15)**

Tags: [P] proved this session (exact rational or polynomial identity), [N] numeric this session, [S1] checked against a source read first-hand, [I] inherited unchecked, [B] belief.

## 1. Verdict

**Outcome (i): the reduction holds identically, for l = 2, 3, 4, on general odd-parity data with field and momentum independent.** The canonical charge, the fixed-volume geometric area side and c_l·E₂[ψ] are the same quadratic form, with zero difference in exact rational arithmetic on every entry of an 8×8 Gram matrix (h-h, K-K and h-K blocks; profiles 1, r², r⁴, r⁶ on each slot), and on a random mixed-m datum [P]. The same holds at l = 5 and 6, and at l = 2 up to r¹² [P]. Beyond the value, the flow of the charge is pointwise the weight-2 scalar flow in the transfer variables, and the transfer is symplectic with constant exactly c_l [P]. Certainty is high for these l on polynomial profiles; for all l and all smooth profiles it rests on the record's symbolic-l radial identity, whose scripts I re-ran but did not re-derive. Step 1 is closed at the level of the charge and its flow; it was never really open, since the record already held a proof (companion A8, Theorem 2) that the referee had not opened.

## 2. Definitions used, and their sources

All in d = 4, unit ball, t = 0, N = (1−r²)/2, factor 1/8G dropped.

- **Data.** a = (h, K), with K the extrinsic curvature. Odd tensors are T = f(r)K_l + g(r)P, with v = x×∇φ, K_l = sym∇v, P = x⊗v + v⊗x (`wave36/companionA8/NOTE.md` line 83). The TT profile parametrisation is Lemma 4 (line 89); I checked div T = 0 and tr T = 0 exactly rather than trusting it.
- **Symplectic form.** Ω(1,2) = ∫_B(K₁·h₂ − K₂·h₁) (NOTE.md line 19).
- **Q1, canonical charge.** A_can = ½Ω(a, δʷa), δʷH = £_ζH + 2tH on the TT-gauge solution, read back as Cauchy data (NOTE.md line 69; `wave69/final/THEOREM.md` line 105). I computed the Lie derivative from the four-dimensional component formula, not from the record's Lemma 3.
- **Q2, geometric fixed-volume area side.** A_geo = −2 × (ε² coefficient of Area(∂B) − 2Vol(B)) at the coordinate ball, using the N-smeared second-order Hamiltonian constraint (NOTE.md Proposition 4, line 133; THEOREM.md H4 line 69 and Section 7). In my run this is ∫N(K·K − (tr K)²) − ∫N·R⁽²⁾[h] − 2δ²F[h], with R⁽²⁾ computed from Christoffel symbols and no TT simplification.
- **Q3, scalar side.** ψ = x·B·x, B = −sym curl K; ψ̇ = −x·curl(K̇)·x with K̇ = −Ric_lin[h]; E_w = ½∫_B N(ψ̇² + |∇ψ|²) + (w/2)∫_B ψ²; c_l = 4(l−2)!/(l+2)! (THEOREM.md line 105; NOTE.md eq. (2) and line 67; `wave73/hinge/REPORT.md` line 24).

## 3. The computation

**Identity [P].** A_can = A_geo = c_l·E₂[ψ] on every odd TT datum tested; there is no h-K cross term; the record's closed form (4) equals A_can.

| Check | Result (l = 2, 3, 4; also 5, 6) |
|---|---|
| T1 value identity, all Gram entries | exact, difference 0 |
| T2 flow: ψ[δK] = N·ψ̇ and ψ̇[δh] = ∇·(N∇ψ) − 2ψ | exact polynomial identities |
| T3 Ω_B(a, δʷb) symmetric, so δʷ is the Hamiltonian field of A on the ball | exact |
| T4 on ℝ³ with Gaussian profiles: Ω_grav = κ·c_l·Ω_scalar | κ = 1 on all 16 pairs |
| T5 (added after the first run): ψ obeys the scalar wave equation | exact |
| S1 random mixed-harmonic datum | exact |
| S2 floating-point Gauss quadrature of the same datum | relative error about 1e−14 |

**Controls.**
- C1, stationary magnetic datum: A = ∮|K|²/((2l+1)(2l+3)) = c_l·E₂ and c_l·E₁/A = (l+1)/(l+2), exact. A/4π = 16/35, 32/105, 512/2079.
- C2, the record's value E₁/E₂ = 204/265 at l = 2 on the profile 1 − 3r²/5: reproduced exactly.

**Plants, all caught.** Wrong lapse (1 − 1.1r²)/2; weight 1 instead of 2 (in the value and in the flow); c_l scaled by (l+3)/(l+2); the +2tH term dropped.

One preregistered plant was badly chosen: c_l → 4(l−1)!/(l+2)! equals c₂ at l = 2 because 0! = 1!, so it is vacuous there. It is caught at l ≥ 3; I added the (l+3)/(l+2) plant after seeing this.

**Gauge.**
- Computations are in TT gauge. In the odd sector K is gauge invariant and h moves by 2·sym∇(a(r)v).
- ψ and ψ̇ are exactly invariant [P].
- A_geo is exactly invariant under a random odd gauge move and in Regge–Wheeler gauge (h_AB = 0, verified) [P]. The reason is that odd gauge vectors are tangent to the sphere, so the coordinate ball maps to itself.
- The closed form (4) is not invariant off TT gauge (for example 20.16 becomes −1003 at l = 4). "A read on the TT representative" is therefore a real restriction for the canonical formula; the geometric definition does not need it.
- The residual TT mode h → h + c·K_l has an all-zero Gram row [P].

**Pure Lie derivative [P].** Without the +2tH term the charge is A + ½∫_B|h|², exactly. That extra term changes under h → h + c·K_l, so it is not gauge invariant.

**General l and profile.** I re-ran the record's `p1_momentum_general_d.py` and `p5_field_half_all_l.py` (symbolic l and n, arbitrary profile η): all residuals are 0. That is a re-run of the record's proof, not an independent derivation. A proof independent of the record would need T2 at symbolic l, which is a short Euler-relation computation.

## 4. What this does to the flagship claim

**May say:** "On every odd-parity graviton sector of a flat ball, the second-order area-side charge (fixed-volume, and equally the canonical charge of the conformal Killing flow with the background rescaling removed) is exactly c_l·E₂ of the scalar chain field, and its Hamiltonian flow is the weight-2 scalar flow. This is proved in the record for all l and checked independently in exact arithmetic at l = 2 to 6. Given the record's identification of the vacuum modular generator with weight 1, this flow is not the modular flow of any normal state or weight in the vacuum sector: the Shale sum diverges logarithmically with per-sector constant 7ζ(3)/(16π²), confirmed on an independent lattice. It is the spacelike-slice counterpart of Satishchandran–Sorce's null-slice theorem, which left gravitons open."

The change from the referee's sentence: "granting our reduction" is replaced by a proved statement for the charge. The grant that remains is smaller: that the transfer carries the ball's standard subspace onto the scalar one with K_Ω = 2π·E₁ (Paper A Theorem 1, [I]). My T4 and T5 support it but do not prove it.

**Must not say:** "We have proved that the Jensen–Sorce–Speranza assumption fails for gravitons, or that no equilibrium state exists for a ball." Also not: "this is the charge JSS define"; see Section 5.

## 5. Secondary findings

**Step 12, JSS (arXiv:2306.01837, pages 15–17 and 21–28 read first-hand as PDF pages) [S1].**
- A4 (p. 16) reads: "The flow generated by H^g_ξ coincides with the modular flow for some state".
- Eq. (3.28): the graviton part obeys {h, (H^g)⁽¹⁾} = £_ξh, the pure Lie derivative, with no rescaling term.
- Section 4 (p. 28) calls H^g_ξ a local integral over the complete Cauchy surface of matter and graviton stress tensors weighted by ξ; H^g_ξ = H^Σ − H^Σ̄ is two-sided.
- Section 3.3 does not fix the region. It proposes (H^Σ_ξ)⁽⁻¹⁾ = 0, which for a diamond is first-order volume fixing, and extremising κA/8πG + 𝒱[ξ], modelled on A − kV, "not treated in detail".

Consequence [P + B]:
- Taken literally, (3.28) on our data gives A + ½∫_B|h|², which is not c_l·E₂ and not gauge invariant. So (3.28) does not define a charge on the gauge-invariant graviton algebra when ξ is not Killing.
- Their own Section 3.3 proposal is the record's fixed-volume construction, which gives c_l·E₂.
- So the record's charge is the natural gauge-invariant reading of JSS, but JSS do not write it, and the literal formula differs. Step 12 stays a gap, now a sharper one. Whether the ball-only charge and the two-sided H^g agree on the ball algebra is also unchecked.

**Step 9, even tower [P by hand, two lines; not machine-checked].** δʷ(£_χη) = £_[ζ,χ]η − 2χ⁰η. Odd gauge vectors have χ⁰ = 0, so δʷ respects gauge there. Even ones do not, so the canonical construction does not carry over as it stands. The record's A_ext on even data rests on the extremisation, and on the even momentum half it reports a spectrum that is not one weight (NOTE.md Section 5, [I]). "The even tower can only add" is unsupported.

## 6. Unchecked, ranked, and the next experiment

1. Step 12: that the gauge-invariant fixed-volume charge is what JSS's A4 is about, including the two-sided H^Σ − H^Σ̄ against the ball-only charge.
2. The standard-subspace half of the reduction: transfer of H(B) onto the scalar ball subspace and K_Ω = 2π·E₁ (Paper A Theorem 1) [I].
3. Step 9, the even tower: no gauge-invariant Hamiltonian charge exhibited.
4. Extremisation (Proposition 4(ii), Schur argument) not computed here [I, B: sound]. Existence of the second-order metric h₂ is assumed.
5. General l by an independent proof; the 2π and 1/16πG bookkeeping between A/8G and 2π·E₂ was checked by hand only.
6. Referee items 4 to 6 (Lemma B, source locators, HINGE Sections 6 and 7): untouched.

**Next experiment:** on the same odd data, compute the charge of JSS's literal prescription in two TT representatives and in Regge–Wheeler gauge, then repeat with their Section 3.3 gauge fixing imposed, to show in exact arithmetic that the gauge fixing removes exactly ½∫_B|h|². That turns step 12 from "not written anywhere" into a stated, checkable identification.

## Reproduce

From `record\wave73\charge-id\scripts`, native Windows Python 3.14, each run 3 to 30 s:
- `python charge_id.py 2 3`, then `3 3`, `4 3`, `5 3`, `6 2`, and `2 6` for the high-degree run.
- Record re-run, from `wave36\companionA8\proof`: `python p1_momentum_general_d.py` and `python p5_field_half_all_l.py`.

Files added, all in `record\wave73\charge-id\`:
- `PREREG.md`
- `scripts\polyalg.py`
- `scripts\charge_id.py`
- `logs\charge_id_l2.log`, `charge_id_l3.log`, `charge_id_l4.log`, `charge_id_l5.log`, `charge_id_l6.log`, `charge_id_l2_k6.log`
- `logs\record_p1_rerun.log`, `record_p5_rerun.log`
- `out\charge_id_l2.json`, `charge_id_l3.json`, `charge_id_l4.json`, `charge_id_l5.json`, `charge_id_l6.json`, `charge_id_l2_k6.json`

`out\charge_id_l2_k6.json` predates the T5 addition; its log is complete. Nothing committed, no report file written. The JSS PDF was saved by the fetch tool outside the repo, in the session's tool-results folder.
