# REPORT: lane SS-TENSION (wave 73)

## Lead's check and grade (2026-11-15)

- Re-run by the lead: `python cone_sector.py` gave 8 lines, 0 unexpected, in 4.4 s. The scaled vacuum kernel reads 12.0 for l = 2 to 6 (worst relative deviation 4.4e-11 at the smallest regulator), the plant "48" is caught, the unbarred correlators vanish. This is the lane's own code re-executed, not an independent implementation. `shear_curv.py` was not re-run by the lead.
- What this lane claims is strong: that one sentence of a published 2026 appendix is too broad at D = 4. The evidence is a sampled high-precision check plus a hand argument, by one lane, using its own reading of (B.26) including a sign it corrected. Until a second, independent derivation agrees, this is to be presented as a question to the authors, not as a correction.
- Standing: the apparent contradiction with our input is resolved in our favour on current evidence; the input (vacuum flow geometric on the curvature, weight 1) now has three confirmations (wave70 proof, subspace exact arithmetic, this cone computation).
- Open items the lane itself lists and the lead keeps: T3 is sampled, not symbolic; tip and edge behaviour is by hand; D other than 4 untested; T4 was uninformative.
- Grade: 368 of 420. First-hand reading with page images, vacuum computed from the covariant propagator with nothing imported, controls on scalar and Maxwell, honest about a wrong pre-registered expectation and an overwritten log; marked down for the overwritten log and the sampled (not symbolic) decisive test.

The lane's report follows unchanged.

**SS-TENSION report (lane, 2026-11-15)**

Tags: [P] proved by hand this session, [N] exact or high-precision arithmetic this session, [S1] source read first-hand, [S2] second-hand, [I] inherited unchecked, [B] belief.

## 1. Verdict

**(R-d) at D = 4, with one (R-b) qualifier. Our input stands, and it now has a third independent confirmation.**

- **The state.** On the gauge-invariant graviton algebra of the diamond in D = 4, Satishchandran–Sorce's state (B.26) is the Minkowski vacuum, normalisation included. Their p. 32 sentence that it "is not the restriction of the vacuum state" is wrong at D = 4.
- **Why.** The D = 4 linearised Weyl tensor is conformally covariant, as the Maxwell field strength is; they list Maxwell D = 4 as the exception and missed this one.
- **The qualifier.** "Geometric modular flow on C-" is true for the curvature (weight 1 on our ψ), but not as a diffeomorphism acting on their shear potential S.
- **Certainty.** High for the D = 4 computation: no free parameter, vacuum kernel constant 12.0 in every sector l = 2..6 at four point pairs. This is a sampled check, not a symbolic identity. D ≠ 4 was not tested, and I do not contest their sentence there.

## 2. Both claims, exactly

**Satishchandran–Sorce, arXiv:2607.20678v1** [S1: PDF fetched, text extraction, page images of printed pp. 24 and 32; v1 is the only version listed]

- **Theorems.** Scalar only: massless Klein–Gordon, D ≥ 3, unit diamond, past null cone C-, variable Π = ∂_u(u^{(D−2)/2}φ). Vacuum: ⟨ΠΠ⟩ = −(1/4π)δ/(u₁−u₂−iε)² (1.4, 3.4), flow 2πu(1−u)∂_u. Theorem 5.1 (p. 24) is uniqueness within H_vac.
- **p. 3.** The spacelike case is "an open question".
- **p. 24.** "massless gauge fields like the free graviton or the free Maxwell field (which is non-conformal in D ≠ 4)". The text extraction drops the slash of ≠; the page image shows it. The parenthesis is attached to Maxwell only.
- **p. 26.** "except in conformally invariant settings (like the Maxwell field in D = 4), the state … will not be the vacuum state"; whether these states lie in the vacuum sector needs "additional work".
- **B.2 (pp. 30–32), any D.**
  - Algebra: γ[f] with f symmetric and divergence free (B.16). These are gauge-invariant observables, the same content as the linearised curvature. So (R-a) is not the answer.
  - Variables: Gaussian null gauge (B.20); free data s_AB, the trace-free angular part in a unit-sphere frame (B.22); Ω = (1/32πG)∫u^{D−6}(s₁∂_u s₂ − s₂∂_u s₁) (B.23); S = u^{(D−6)/2}s (B.24); γ[f] = −(1/16πG)∫F^{AB}∂_u S_AB (B.25).
  - State (B.26): (1/32π²G)∫F₁F₂ δ/(u₁−u₂−i0⁺)².
- **p. 32, verbatim.** "It is straightforward to show that ω₀^GR defines a state on A^GR(D) which has geometric modular flow on C⁻. Since the linearized gravitational field is not conformally invariant, ω₀^GR is not the restriction of the vacuum state to D."
- **Sign.** (B.26) is printed with + where (B.13) and (3.4) carry −. With + the form is not positive, so I read it with the minus.

**Ours**

- `wave70\modflow\REPORT.md` lines 24–38: D = 4, K(O_R) is the closed span of smeared (E, B); Δ^{it}e_F = e_{Λ_{−2πt*}F}; the curvature is transported with weight 3. Line 56: the curvature two-point function is invariant under all 15 conformal Killing fields.
- `wave73\subspace\REPORT.md` line 19 and section 3: ψ = x·B·x, vacuum ↦ c_l × scalar vacuum, K_Ω = 2π·c_l·E₁[ψ].

Both sides speak about the same gauge-invariant algebra of the same diamond.

## 3. Computation and results

All on the cone t = r = u, with k = ∂_t + ∂_r and Ψ₀ = R(k, m, k, m). The vacuum comes from the textbook de Donder propagator only, so nothing is imported from earlier lanes.

**T0 [N, exact].** In their gauge at D = 4, R_{uAuB} = −½·u·∂_u²(s_AB/u), so Ψ₀ = −S″/(2u). Then (B.26) predicts ⟨Ψ₀Ψ̄₀⟩ = 12G·δ_{S²}/(u₁u₂(u₁−u₂−i0)⁴) in every sector.

**T3, decisive [N].** The regulated vacuum correlator is projected on d^l_{22}, with exact Gaussian rationals plus one logarithm. The scaled value g_l·u₁u₂(u₁−u₂)⁴/G should be 12:

| ε | result |
|---|---|
| 10⁻⁴ | 11.99999… |
| 10⁻⁸ | 12.0 + O(10⁻⁷)i |
| 10⁻¹² | 12.0, worst relative deviation 4.4·10⁻¹¹ |

- This holds for l = 2, 3, 4, 5, 6 at four (u₁, u₂) pairs.
- The unbarred ⟨Ψ₀Ψ₀⟩ vanishes (10⁻¹³⁷).
- So the two states agree outright: no weight difference, no soft-sector difference, no normalisation difference.

**Controls [N].**
- (i) Scalar: g_l(u₁−u₂)² → −1/(4π) for l = 0..4 (6.8·10⁻¹¹).
- (ii) Maxwell D = 4: g_l·u₁u₂(u₁−u₂)² → −1/(4π) for l = 1..4. This is (B.13), where they and we agree.

**Plants, all caught.**
- Scalar with weight 2 (l-dependent, off by a factor of about 5).
- Graviton constant 48.
- T0 without background terms.
- S used in place of Q in T5.
- Coefficient 1 for 2 in the transform T.
- The literal S-flow passed off as local on Ψ₀.

**T5 [N, exact].** With Q_AB = S_AB − (2/u)∫₀ᵘS_AB:
- div div Q = −2u(x·E·x) and curl div Q = +2u(x·B·x) on the cone.
- The constants are the same on 15 polynomial plane-wave solutions: degrees 2..6, three polarisations, hence l ≤ 6 and both parities.
- So their potential and our chain field are related by one integration along the generator; S itself is not u·ψ.

**T6 [N, exact].** That transform preserves the half-line ultralocal form: 16 complex entries are equal, so both the symmetric and the symplectic parts agree. Its Mellin multiplier (s−1)/(s+1) has modulus 1 [P]. Every curvature observable has ∫F du = 0.

**T7 [N on monomials, P].**
- The vacuum flow is weight 0 on Q = uψ, i.e. weight 1 on ψ.
- On S it reads δS = u(1−u)S′ + 2∫₀ᵘS.
- On Ψ₀ it is local: δΨ₀ = NΨ₀′ + (3−5u)Ψ₀, so uΨ₀ is a weight-2 primary with the 1/u⁴ kernel.
- The literal weight-0 flow on S differs by −S′/u on Ψ₀, which is non-local in the curvature, and it does not preserve ∫F du = 0 [P].
- For Maxwell the field strength is a first derivative of a weight-0 potential, so no such term arises. That is why their Maxwell statement is fine and the graviton one is not.

**Which weight their flow carries in our variables.** Weight 1, because their state is the vacuum. The literal S-diffeomorphism is not a weight-w flow on ψ for any w.

**By-hand cross-check [P, not machine-checked].** (B.26) with T5 and the standard harmonic eigenvalue ½(l+2)!/(l−2)! gives 1/c_l = 4πG(l+2)!/(l−2)!. At 16πG = 1 that is 6, 30, 90, the record's values.

**Disclosures.**
- T4's pre-registered expectation was wrong. Ψ₀ contains only derivatives tangent to the cone, so the off-diagonal correlator vanishes in every D. The test is uninformative and is reported as INFO.
- T6 failed its first run because sympy returned unevaluated Meijer G functions. A second attempt with `nsimplify` produced spurious closed forms; that log was overwritten. The final version uses exact Gamma values. The registered statement did not change.
- T3 was widened after the first run to l = 2..6 and four point pairs.

## 4. Flagship

**Amended "may say"** (jss-id Section 5 paragraph unchanged up to "…the vacuum's modular generator", then add):

"This input is confirmed a third way: the vacuum two-point function of Ψ₀ on the diamond's null cone, computed from the covariant propagator, is 12G·δ/(u₁u₂(u₁−u₂)⁴) in every sector l = 2..6, the two-point function of a weight-2 primary, so the vacuum's flow is geometric on the curvature."

Replace the last sentence of subspace Section 4 ("It is the spacelike-slice counterpart … which left gravitons open") with:

**The one sentence on Satishchandran–Sorce.** "In D = 4 the state they construct for gravitons on the null cone (their B.26) is the Minkowski vacuum restricted to the diamond's gauge-invariant algebra, since the linearised Weyl tensor is conformally covariant there just as the Maxwell field strength is; its modular flow is the conformal flow on the curvature (weight 1 on ψ), which on their shear potential is u(1−u)∂_u plus a term 2∫₀ᵘS and not a pure diffeomorphism."

**Must not say:**
- "Satishchandran–Sorce contradict us", or "they showed the graviton null-local state is not the vacuum".
- "Their paper is wrong": one sentence of an appendix is too broad at D = 4, and the theorems are untouched.
- "Their uniqueness theorem covers gravitons": Theorem 5.1 is scalar only.
- Anything about D ≠ 4.
- "Proved": l = 2..6 at sampled points to 10⁻¹¹ is evidence, backed by wave70's exact conformal-invariance lemma.

**Open question.** "In D = 4, isn't B.26 just the vacuum? uΨ₀ = −S″/2 has the 1/u⁴ two-point function of a weight-2 primary, which is what the de Donder propagator gives on the cone."

## 5. Unchecked, ranked

1. T3 is at sampled points, not a symbolic limit [N].
2. That Ψ₀ on C- generates the whole diamond algebra is their own premise plus the characteristic initial-value problem [B, standard].
3. Tip and edge behaviour [P, by hand]: in the tip-regular gauge F = Ef carries a tail b(x)·u to the edge, so literal (B.26) diverges logarithmically; the tip-anchored potential's b-mode has infinite vacuum variance. Not machine-checked.
4. D ≠ 4 untested. Whether their "geometric flow" claim survives the same derivative-of-a-potential issue there is open [B: probably not literally].
5. The harmonic eigenvalue used in the c_l cross-check is from memory.
6. Literature: Longo–Morinelli–Preta–Rehren's abstract says finite-helicity representations "extend to the conformal group" [S1, abstract only]. Benedetti–Casini take the curvature-generated algebra [S2, search summary]. Hislop 1988 remains unread.

**Next experiment.** Take the ε → 0 limit of T3 symbolically in (u₁, u₂) for l = 2, then repeat T3 at D = 6 with tensor harmonics on S⁴ to show the kernel becomes l-dependent. That would turn "D = 4 is the exception" into a two-sided statement.

## Reproduce

From `record\wave73\ss-tension\scripts`, native Windows Python 3.14:
- `python cone_sector.py ..\out\cone_sector.json` (about 5 s; 7 lines as expected plus one INFO)
- `python shear_curv.py ..\out\shear_curv.json` (about 20 s; 13 lines as expected)

Files added, all in `record\wave73\ss-tension\`:
- `PREREG.md` (with a post-run errata section)
- `scripts\cone_sector.py`
- `scripts\shear_curv.py`
- `logs\cone_sector.log`
- `logs\cone_sector_run1.log`
- `logs\shear_curv.log`
- `logs\shear_curv_run1.log`
- `out\cone_sector.json`
- `out\shear_curv.json`

Nothing committed, no report file written, no subagents, no processes killed. The fetched PDF sits in the session scratchpad, outside the repo.

Sources: [arXiv:2607.20678](https://arxiv.org/abs/2607.20678), [arXiv:1806.04967](https://arxiv.org/abs/1806.04967), [arXiv:1908.01800](https://arxiv.org/pdf/1908.01800)
