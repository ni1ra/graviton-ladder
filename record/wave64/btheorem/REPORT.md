# REPORT: w64-btheorem (saved by the lead; the harness refused the write)

**w64-btheorem: report** (the harness blocked REPORT.md, so this is the full report)

**Path and sha.** `<local path>`, final commit `b389c3bea`, pushed 18:55 Oslo. The PREREG is `6005d1de7`, pushed 18:43:05 before any derivation script ran; its hash is in `out/prereg.sha256`. The folder holds NOTE.md (theorem and proof for a hostile referee), PREREG.md, HASHES.md (every hash verifies), and `scripts/` and `out/`.

**Verdict**
1. B at d = 4 now follows from four premises: the first-order first law against the vacuum (FL1, already in print: BMSV and JSS (5.35)), one Taylor functional at G⁰ (P0), time reflection (T) and matter universality (U). With a fifth premise, Q0 (no corner term quadratic in K, added after the PREREG), the area side on the scope family is exactly A.
2. The proof needs none of: SO(3) invariance, corner gauge, the graviton's modular flow, Takesaki or BW on the curvature net, JSS's A4 for the graviton, Ψ_ζ, the converse-cocycle lemma, or any truncation. The weight N comes from matter (Hislop–Longo) at first order.
3. No step is false. Every exact computation passes and every plant is caught. The computations are identities, though, and give no evidence on B by themselves; the evidence is the reduction.
4. A referee will spot this at once, and I concede it: FL1 already is first-order (BR) on a set of data that spans everything, and Lemma 6 extends it linearly. The open part of B was always the step to second order (Lemma 7).

**Steps and their status**
| Step | Status |
|---|---|
| S1 Green identity, ∫N R_lin = −2δ(A−2V) for every q | Proved; computed exactly; the k = 2.2 plant is caught |
| S2 δ(A−2V) is unchanged by every corner displacement | Proved; computed. Not a real test: it holds automatically for a round unit sphere (k = mean curvature) |
| S3 Spanning lemma; for scalar data (0, g), T₀₀ = ½g², T₀ᵢ = 0, Wald term 0 | Proved; computed |
| S4 First law for mixtures, S(ρ_λ‖Ω) = o(λ) | Proved from cited inputs (Donald's identity, convexity, lower semicontinuity). The von Neumann form of Donald's identity was not re-read this session |
| S5 S(ω_g‖Ω) = 2π∫N·½g² | Cited (Longo–Morsella), not re-read this session |
| S6 FL1 + T fix the metric-linear part on the physical slice | Proved. The rival functionals fail as they must, which is automatic given FL1 |
| S7 Taylor transfer to second order (P0) | Quadratic free-net entropy cited (CLR, Paper A Thm 1); order counting proved; matter-loop terms argued from G counting |
| S8 The second-order metric lies on the physical slice up to quadratic terms | Proved; exact for JSS's fixed-volume slice on the scope family |
| S9 Scope family area side = A | Computed exactly: 16π/35 at l = 2 and 128π/105 at l = 3; ratio to flux 1/35 and 1/63, matching Paper A's gate |
| S11 Record route: l = 0 classification | Exact with arbitrary profiles, no polynomial ansatz. This closes the truncation in w39 S2 and br6 T2 |

**The four open items you named**
- Converse-cocycle lemma: not needed. If one insists on JSS's construction, w61's edge cure changes the result by only O(δ²) on smooth data.
- (c′) at d = 4: bypassed; no equivalent premise comes in.
- Corner gauge (c3): harmless for the linear clause, by argument rather than new evidence. The five monopoles from w39 S2 are either excluded by FL1 or harmless on the slice (δA on JSS's slice). What remains is quadratic, and that belongs to C.
- Truncation: bypassed by the proof, and closed separately by S11.

**LR to book after a blind check: 1.87 on B.** B goes 0.872 → 0.927 and LADDER-WHY 37.41% → 39.78%.
- The PREREG registered 1.96 before Q0 existed. Pending the blind check the figure is 1.54; don't book that one.
- The range is 1.48–2.26 as p_prem(c3) runs from 0.3 to 0.8. The whole LR rests on these judgements, and I set them after I had the route.
- If you reject the c3 argument, the floor is 1.23 (→ 38.31%).
- Your bookkeeping ruling on w39's B/C double count (B's factor read as P(clause 3 | clause 2)) would be worth up to 2.43 by restatement alone. It overlaps the c3 part of the 1.87, so never multiply the two.
- A simpler floor you can book from your own beliefs: P(B | C) ≥ P(P0 ∧ T ∧ U ∧ FL1 ∧ Q0 | C).

**Weakest step: P0 together with Q0.** Is the second-order area side at G⁰ really the Taylor continuation of one fixed metric functional, with no K-quadratic corner term? FL1 is published first-order physics, so any world where B fails has to live in a failure of P0 or Q0.

**Disclosures**
- Two tooling errors, both fixed before any reported number. The first s1s2 run reported a failure caused by an `nsimplify` bug in the integration helper. The first s11 run crashed on a sympy derivative with respect to √(x²+y²+z²).
- Changes after the PREREG:
  - Q0 added as a premise.
  - `price.py` amended. Only a block was appended, and a rerun-diff confirmed the registered output is unchanged.
  - NOTE reworded after the advisor review.
- Not done: no blind re-derivation, no outside model.

**Suggested wording for Paper A App. D:** "(BR) follows from the first-order first law of the diamond (BMSV; JSS (5.35)) by Taylor transfer, with no rotation-invariance, corner-gauge or graviton-modular-flow premise."
