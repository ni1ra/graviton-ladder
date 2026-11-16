# BEYOND-NORMAL: is there any state at all for the area-side flow?

## Lead's check and grade

**Grade: 372 of 420.** Accepted as "exists, quasi-free, disjoint given the theorem"; not yet refereed.

**Regraded 352 of 420 on 16 Nov 2026** (was 372) after its referee (`beyond-ref\REPORT.md`): the uniqueness claims in 2(a) and 2(c) are false as stated (correct: unique among states continuous in the thermal norm), and one inequality in 2(b) is false (the referee's proof replaces the step). Existence, the modular flow and disjointness stand.

**Checked by the lead.** (1) Re-ran `C:\Python314\python.exe -X utf8 beyond_normal.py <out.json>` from `scripts\`: exit 0, "FAILED LINES: none", both plants caught, the l = 64 affinity row equals the lane's (0.05102 against HINGE's 0.0510). (2) By hand: stability is immediate from the dictionaries (odd: h² + N with N > 0; even: (1 + W) with W ≥ 0), so each sector's generator is a positive quadratic Hamiltonian and a thermal quasi-free state at inverse temperature 2π exists per sector. (3) The logic of disjointness: if the state were quasi-equivalent to the vacuum it would be a faithful normal state with the excluded modular flow, which Theorems E/O forbid.

**Not checked by the lead.** The uniform equivalence of the two norms (the step that makes the product a state on the same Weyl algebra): Hardy's inequality constant 3.52 and the interpolation step; the two standard facts quoted from memory (second quantisation of the modular operator, primarity of quasi-free states); uniqueness among all KMS states; anything about the state's physical acceptability. The direct Araki–Yamagami criterion was not computed.

**What it changes.** The question "which state, if any, does the area-side flow belong to" has an answer: one exists, per sector and on the tower, and it lives outside the vacuum representation of the ball. The record's affinity constants (EVEN-CONST, HINGE) are distances from the vacuum to exactly this state. The sentence "no equilibrium state exists" stays banned, now because it is false as stated for the Weyl algebra. A hostile referee (card BEYOND-REF) is assigned before the paper leans on this beyond one hedged sentence.

**Bearing on the wider question.** The gravitational-algebra literature conjectures states whose modular flow is geometric. For free gravitons on a ball this says: such a state exists for the area-side flow, and it is not in the vacuum's folium. Whether it is the conjectured state, and whether it is physically acceptable, is open.

---

## Lane report (unchanged)

**Q11 BEYOND-NORMAL report (16 Nov 2026)**

**1. Answer.** For g = 1, both parities, every l, the area-side sector flow is stable, so a quasi-free 2π-KMS state ψ_g exists on the tower's Weyl algebra, and its modular flow on its own von Neumann algebra is the area-side flow. Given the record's theorem, ψ_g is disjoint from the vacuum representation; it is the natural candidate, and it lives outside the vacuum sector.

**2. Argument.**

(a) Stability and covariance.
- [P] Odd: H = ½⟨π,π⟩ + ½⟨q,(h²+N)q⟩ with N > 0, so h_g² = h²+N ≥ h² > 0.
- [P] Even: H = ½⟨π,(1+W)π⟩ + ½⟨q,h²q⟩ with W ≥ 0. Its norm is ‖W_l‖ = 1/(2l+3), with positive eigenfunction π = r^(l+1)(1−r²).
- [P] Normal coordinates q′ = A⁻¹q, π′ = Aπ, A = (1+W)^½ (A = 1 for odd), Ω² = A h_g² A. Then μ_g(k) = ½⟨q′, Ω coth(πΩ) q′⟩ + ½⟨π′, Ω⁻¹coth(πΩ) π′⟩. It is the only quasi-free 2π-KMS covariance, since Ker Ω = 0.

(b) Product state and disjointness.
- [P] μ_g and the thermal norm μ_0 are equivalent uniformly in l, so ψ_g = ⊗ψ_g^(l,m) is a state on the same Weyl algebra over H.
  - The energy forms satisfy E_0 ≤ E_g ≤ (1+C)E_0, and E(X⁻¹·) likewise, with C = 1/(2l+3) for even.
  - For odd, C ≤ 3.52, from N ≤ 0.88/ρ² and Hardy's inequality.
  - coth(πx)/x lies between 1/π and 1 times (x⁻² + x⁻¹).
  - [S, from memory] Hilbert-space interpolation (Heinz) covers the middle power.
- [P given S] The referee's Tomita computation with Ω in place of h gives σ^ψ_t(W(k)) = W(U^g(2πt)k) on M_g = π_ψ(𝔚)″. The second-quantisation citation (Eckmann–Osterwalder; Leyland–Roberts–Testard) is from memory.
- [P given I: Theorems E/O, refereed] If ψ_g were quasi-equivalent to ω, transporting it would give a faithful normal state on M with the excluded modular automorphism at t = ½. So it is not quasi-equivalent.
- [S, from memory] ψ_g and ω are primary (Ker Ω = 0; Araki, Manuceau–Verbeure), so not quasi-equivalent means disjoint.
- Which theorem: none of the three is needed, and none applies to the Shale sum directly. The cocycle S is not the map carrying ω to ψ_g; the record's link is Connes plus Shale.
  - [S, from memory, not opened] Araki–Yamagami 1982 is the direct criterion for two quasi-free CCR states: same topology and √S_g − √S_0 Hilbert–Schmidt.
  - van Daele 1971 is the earlier CCR version.
  - Powers–Størmer 1970 is for fermions (CAR) and must not be cited here.
- [U] The Araki–Yamagami condition itself was not computed.

(c) Candidate and what is unphysical.
- [S, from memory: Rocca–Sirugue–Testard 1970, hypotheses unchecked] Uniqueness among all KMS states, not only quasi-free ones.
- [N] Return to equilibrium: the vacuum evolved by the area-side flow converges to ψ_g on localised data.
- [U, HINGE only] Infinite relative entropy to the vacuum, with the defect logarithmic at the edge; Hadamard behaviour in the interior is not examined.
- [P] In chain-field variables the even generator contains the ball's Dirichlet Green function, so that flow is non-local.

**3. Check.** Grid ρ ≤ 50, step 0.025, HINGE's walls; first run, exit 0, 39 s.

| Line | l = 2 | l = 16 |
|---|---|---|
| λ_min: vacuum / odd / even (×10⁻³) | 1.0484 / 1.0795 / 1.0485 | 1.1347 / 1.1415 / 1.1347 |
| Odd sup N/h² (bound 3.52) | 0.33334 | 0.05882 |
| (2l+3)‖W‖ | 1.00007 | 1.00009 |
| μ_g/μ_0 extremes, odd | 0.817 to 1.089 | 0.961 to 1.018 |
| μ_g/μ_0 extremes, even | 0.932 to 1.061 | 0.986 to 1.013 |
| μ_g drift under U^g(1.7), worst | 4.0e-15 | 3.1e-15 |
| μ_0 drift under U^g(1.7), odd / even | 9.3e-3 / 5.3e-3 | 2.4e-3 / 2.9e-3 |
| r(τ) at 5, worst case; at 20 | 1.036; 1.0000 | 1.011; 1.0000 |
| Odd l²·mean(1−A) against HINGE | 0.02514 (0.0251) | 0.04600 (0.0460) |
| Even λ²(1−A), reported only | 0.00470 | 0.01057 |

- r(τ) = [μ_0(U^g(τ)k) − μ_0(k)]/[μ_g(k) − μ_0(k)]; all 8 cases reach 1.0000 by τ = 20. This is what I verified in place of the card's cocycle wording.
- The odd affinity matches HINGE's table at l = 64 too (0.05102 against 0.0510), so the record's affinity constants are distances to exactly this state.
- Toy, n = 6: the generator read off the covariance alone equals JG to 5e-13 in all three types.
- Odd sup N/h² equals 1/(l+1) to five digits at l = 2, 16, 64; this is an observation, not proved.
- Plants (wrong temperature π): toy mismatch 0.5000, caught; return target r(20) = 0.354, caught.
- Failed lines: none.

**4. Re-run.** From the lane's `scripts` folder:

`C:\Python314\python.exe -X utf8 beyond_normal.py ..\out\beyond_normal.json`

**5. Weakest step.** Disjointness is inherited from Theorems E/O plus two from-memory standard facts (second quantisation of the modular operator; primarity); the direct Araki–Yamagami condition is not computed.

**6. Bulletin.** Read (B1 to B27); no steer.

**7. Sentences.**
- **May use:** "The area-side flow is stable in every sector and has a quasi-free KMS state at inverse temperature 2π on the Weyl algebra, whose modular flow is that flow; by the theorem this state is disjoint from the vacuum representation of the ball, and whether it is physically acceptable (its relative entropy to the vacuum appears to diverge at the edge) is open."
- **Must not use:** "No equilibrium state exists", "the Shale sum shows the two states are inequivalent by Powers–Størmer", or "the state is Hadamard / unphysical" as a proved statement.
