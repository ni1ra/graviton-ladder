# CONE-REF: hostile referee of the cone forms

## Lead's check and grade

**Grade: 384 of 420.** Accepted. No break; one real gap found and closed by the referee.

**Checked by the lead.** (1) Re-ran both probes (`python scripts\bracket_exact.py`, `python scripts\probe2.py`, from the lane's folder): exit 0; logs equal the report (bracket form 1/8, 15/32, 15/32 for p = 1 and p = 3 against 1/4, 15/16, 25/32 for p = 2; Möbius covariance true only for p = 3; first-law integration by parts and the RULER-L identity with zero residual). (2) The multiplier lemma by hand: (s + p − 1)(s + 3 − p) = s² + 2s + (p − 1)(3 − p), so the product s(s+2)/((s+p−1)(s+3−p)) equals 1 exactly when p = 1 or p = 3.

**Not checked by the lead.** The cone symplectic form and its normalisation (taken as standard by the referee, sources not opened); the extension of the weight-3 flow from curvature smeared inside the double cone to its boundary cone (assumed by continuity); the lemma as a theorem with function spaces (exact on powers and three exact pairs only).

**What it changes.** "S = Q_3" on the cone is no longer a match at l ≤ 4: with the lemma it is derived for every l from the refereed weight-3 flow (B21) and the canonical bracket, with no free constant, and the weight is forced (p = 2, the shear, is not a canonical current). "Q_2 is the cone first-law charge" and duality hold. Two wording corrections for the ruler reports: "diagonal in parity" is parity orthogonality, not duality; "the only linear functional" means on area profiles.

**Still hypotheses.** The ball's algebra is the linearised-curvature net; the standard cone symplectic form; continuity to the boundary cone.

---

## Lane report (unchanged)

**CONE-REF final report.** No break found. "S = Q_3" had one real gap, which I closed with a short lemma.

**1. Verdicts**
- **S = Q_3: GAP (fixable, repair below).** The weight p = 3 is derived, not fitted. With the repair it is proved for every l, given the weight-3 flow of B21 and the standard cone symplectic form. The report itself only matched it on 18 data at l = 2, 3, 4.
- **Q_2 is the cone first-law charge: HOLDS, proved.** It is an identity in the source profile.
- **Duality: HOLDS, proved.** One wording slip, listed below.

**2. Defects, most serious first**
- **Missing step in "So S = Q_3".** The report shows only that Π = 2vσ_3 transforms as a weight-one current. That Q_3 is the generator also needs Π to have the canonical bracket, which is not automatic: Π is non-local in the canonical variable j = ∂_v(vγ) = 2vσ_1.
  - Repair: on powers j = v^s, vσ_p = [s/(s+p−1)]·vσ_1. The bracket ∫∫ j_1 sgn(v−v′) j_2 pairs s with −2−s. The product of multipliers is s(s+2)/((s+p−1)(s+3−p)), which equals 1 iff (p−1)(3−p) = 0.
  - So exactly p = 1 and p = 3 are canonical currents, and the shear p = 2 is not. Covariance under the diamond's Möbius maps with curvature dimension 3 selects p = 3; p = 1 is the cone-local state.
  - A canonical current of weight one has the unique quadratic generator ∫fΠ². With Wall's normalisation (bracket 1/32πG in vγ) this gives S = 2π·(1/32πG)∮∫f|Π|² = (1/4G)∮∫v²f|σ_3|², with no free constant.
  - This answers "could another kernel pass": no.
- **The tag "UNVERIFIED locator" can go.** C′(w) = a⁻²Ω⁻⁵C(v) follows from φ*C′_abcd = Ω·C_abcd (dimension 3, four lower indices) with dw = aΩ²dv.
- **Chain route proved only at l ≤ 4.** The pointwise identity between ∂_v(vψ_E,B) and Π is checked at l = 2, 3, 4 only. The route above does not need it; if the paper cites the chain route for the cone, that identity needs a general-l proof.
- **Wording, section 7.** "Every cone form is diagonal in parity by duality" is wrong: that is parity orthogonality of the even and odd harmonics. Duality gives equality of the electric and magnetic values.
- **Wording, section 3.3.** "The bracket is the only linear functional" is unique as a functional on area profiles; say that.

**3. Checked by hand and found sound**
- f = v(1−v/R) from the diamond's conformal Killing vector, and the scalar check that ∂_v(vφ) is the weight-one current.
- First-order δθ = 0, so the dropped (δθ)² term is fourth order and the cone loses nothing there.
- δA(v) = −∮∫_0^v u·v(v−u)S du, and the bracket equals −∮∫u²f·S.
- σ_2 is gauge invariant (zero background shear) and the linearised Weyl component is invariant, so Q_2 is gauge invariant at second order.
- Duality: C_AB → ε_A^C C_CB commutes with every σ_p and preserves |·|² pointwise.
- Fairness: same C, same [0,R], same affine v, same G, same contraction σ_ABσ^AB. The bracket is also invariant under affine rescaling.
- RULER-L's identity, by hand with J = ∫u²σ_2.
- The report's flow δψ = fψ′ + (2/R)∫ψ is equivalent to δΠ = (fΠ)′.

**4. Probe (exact sympy)**
- Bracket form on three data pairs: p = 1 gives 1/8, 15/32, 15/32; p = 3 gives the same; p = 2 gives 1/4, 15/16, 25/32.
- Finite Möbius covariance: true for p = 3, false for p = 1 and 2, at two (a, R) values.
- First-law integration by parts: zero residual for symbolic n and R.
- RULER-L identity with cross terms: zero residual.
- Files: `scripts\bracket_exact.py`, `scripts\probe2.py`, and `logs\` in the lane's folder.

**5. Not checked**
- The cited TeX sources (arXiv:1703.10656, Wall, Wang) were not opened.
- The graviton cone symplectic form was taken as standard, not derived.
- B21 is stated for curvature smeared inside the open double cone; its extension to the boundary cone is assumed by continuity.
- The lane's scripts were not re-run; I read `nc_core.py` only, and none of the sibling reports.
- Section 3.4 (surfaces, Wang's table) and 3.5 were not refereed.
- The multiplier lemma is exact on powers and three exact pairs, not written as a theorem with function spaces.

**6. Bulletin**
B21 changed my route: I prove S = Q_3 directly from the refereed weight-3 curvature flow, not by transporting the chain scalar (B17, B22, B27 not needed). B23's identity is confirmed independently. Nothing else changed.

**7. Sentences**
- **May use:** "On the cone, the vacuum modular energy of the graviton is S = (1/4G)∮∫v²f|σ_3|², because 2vσ_3 is the unique canonical current that the weight-3 conformal flow transports locally; the first-law area functional equals the same expression with the shear σ_2, and exceeds it by (1/4G)∮∫v³|σ_2−σ_3|² for every non-zero cone curvature, identically for electric and magnetic data."
- **Must not use:** "The entropy on the cone was verified to be Q_3" (it was matched at l ≤ 4), or any sentence saying a "cone proof" shows this without naming the hypothesis that the ball's algebra is the linearised-curvature net.
