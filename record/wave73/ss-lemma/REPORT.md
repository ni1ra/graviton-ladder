# SS-LEMMA: the D = 4 four-moment lemma, proved

## Lead's check and grade

**Grade: 386 of 420.** Accepted.

**Checked by the lead.** Re-ran both scripts from the lane's folder (`C:\Python314\python.exe scripts\lemma_proof.py`, then `scripts\lemma_direct.py`): "SCRIPT A PASS", "SCRIPT B PASS", exit 0 each; the D = 5 control lines read 55 of 55 non-zero and the k = 4 line 21 of 21 non-zero, as the table states. The structure of the proof by reading: steps 2 to 4 reduce each moment to k(k−1) times one integral, so the first two moments vanish for any symmetric source; steps 5 exhibit the k = 2 and k = 3 integrands as symmetrised gradients, which integrate to zero against a divergence-free source.

**Not checked by the lead.** The two explicit potentials ξ in step 5 by hand (the script verifies them component by component in two frames, and its plant catches a changed coefficient); step 1, which the lane recalls from memory.

**What it changes.** The question the paper puts to the authors of arXiv:2607.20678 no longer rests on eleven configurations: in D = 4 the four moments vanish for every divergence-free source supported in the diamond, and the lemma is sharp (the fifth moment does not vanish; D = 5 and 6 fail). The question stays a question about one sentence; their theorems are scalar-only and nothing here says any is wrong.

**Process note.** Two code faults were fixed between runs (an assert, a float-zero miscount); neither changed a pass line, and both first-run logs are kept.

---

## Lane report (unchanged)

**SS-LEMMA (Q12) report**

**1. Verdict.** Proved for all sources in the class, in D = 4: for every smooth, symmetric, divergence-free f compactly supported in the open diamond, ∫u^k Ψ₀[Ef](u·l) du = 0 for k = 0, 1, 2, 3 on every cone generator (ss-second's M₋₁..M₂). No failing source exists. The proof never uses "curvature smearings exhaust γ[f]", so ss-second's unchecked item 2 is not needed for the lemma.

**2. Proof.** Signature (−+++), l = (1, n), m complex null tangent at n, w = −l·y, μ = m·y, s = −y·y, v = μl + wm.
1. [S, from memory] On the cone, Ef is (up to sign) the advanced solution h = G∗f̄, with G ∝ δ(σ). The trace part drops because η(l,l) = η(l,m) = η(m,m) = 0, and Ψ₀ is gauge-invariant.
2. [P] On the support, δ((ul−y)²) fires at u_y = s/(2w), which lies in a compact subset of (0,1). So h has compact u-support on each generator and ∫u^j h_ab(ul) du = c·2^(−j−1) ∫f_ab s^j w^(−j−1) d⁴y.
3. [P] Ψ₀ = −½[(l·∂)²h_mm + (m·∂)²h_ll − 2(l·∂)(m·∂)h_lm] contains only cone-tangential derivatives: l·∂ = ∂_u and m·∂ = u⁻¹ times the derivative along n → n + εm, which stays on the complexified sphere. Integrating by parts in u gives M_k ∝ ∫u^(k−2)[(k−1)(k−2)H_mm + 2(k−2)∂_εH_lm + ∂_ε²H_ll].
4. [P, C: A1] Inserting step 2 with l_ε = l + εm and w_ε = w − εμ gives M_k = const·k(k−1) ∫ s^(k−2) w^(−k−1) f(v,v) d⁴y. Hence M₀ = M₁ = 0 for any symmetric f, conserved or not.
5. [P, C: A2] For k = 2, v_a v_b/w³ = ∂_(a ξ_b) with ξ = (μ²/2w²)l + (μ/w)m. For k = 3, s·v_a v_b/w⁴ = ∂_(a ξ_b) with, for p = −l′·y, ζ̄ = m̄·y, l·l′ = −2, m·m̄ = 2: ξ = (pμ²/2w² − μ³ζ̄/3w³)l + (μ²/2w)l′ + (pμ/w − μ²ζ̄/2w²)m + (μ³/6w²)m̄. Both are smooth on w > 0, so ∫f^ab ∂_a ξ_b = 0 by ∂_a f^ab = 0. ∎
6. [C] For k = 4 the tensor has non-zero linearised Riemann, so the lemma is sharp at four moments.

**3. Checks** (pass lines fixed in `logs\00_passlines.log` before any run):

| Line | Content | Result |
|---|---|---|
| A1 | ε-algebra, symbolic k | residual 0 |
| A2 | sym-grad ξ = target, 2 frames, 40 components | all 0 |
| A3 | linearised Riemann of the k-tensor | 0 at k = 2, 3; 21/21 non-zero at k = 4 |
| A4 | D = 5 control, reduction | 55/55 non-zero at k = 2 and k = 3 |
| A5 | plant, 1/3 → 1/2 in ξ₃ | caught, 7 residuals |
| B1 | direct from definitions, symbolic point p, all 21 curvature components | M₀..M₃ ≡ 0, 84 identities |
| B2 | teeth | M₄ ≢ 0 for 21/21 |
| B3 | generic direction, 4 random p | 0 of 336 non-zero |
| B4 | non-conserved source | M₀ = M₁ = 0, M₂, M₃ ≠ 0, M₃/M₂ = 3u_x = 411/220 exactly |
| B5 | plant, cross term −1 | caught |
| B6 | D = 6 direct | M₃ ≠ 0 for 20/240 sources; M₀..M₂ = 0 |
| B7 | D = 5 direct, finite part | all four moments non-zero, 110/110 |

B1 is identically zero in p, so it covers every point multipole and every smooth curvature smearing. That is a second, independent proof for that subclass.

Failed lines: none in the record runs. Two code faults were fixed, neither changing a pass line: `01_proof_run01` aborted on an unexpanded assert, and `02_direct_run01` miscounted float 0.0 as non-zero in B6 (30/30/30/50 printed, true counts 0/0/0/20).

**4. Re-run.** From the lane's folder: `C:\Python314\python.exe scripts\lemma_proof.py` (under a minute), then `C:\Python314\python.exe scripts\lemma_direct.py` (about 3 minutes).

**5. Weakest step.** Step 1, that the commutator function on the cone is the scalar advanced Green function acting on conserved f in de Donder gauge, is standard but recalled from memory; B7 additionally rests on the finite-part form of the D = 5 Green function, checked only by A4's agreement.

**6. Bulletin.** Read, B1 to B27; no steer.

**7. Sentences.**
- May use, as a question: "In D = 4, for every divergence-free f supported in the diamond, the first four u-moments of R(l,m,l,m)[Ef] vanish on each generator of the cone, so (B.26) converges and equals the vacuum two-point function on the γ[f] algebra; is the sentence that ω₀^GR is not the restriction of the vacuum meant for D ≠ 4 only?"
- Must not use: "the appendix, or any theorem of arXiv:2607.20678, is wrong", or "the lemma holds in every dimension" (it fails in D = 5 and 6).
