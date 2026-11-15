# REPORT: lane SS-SECOND (wave 73)

## Lead's check and grade (2026-11-15)

- Re-run by the lead: `s2_kernel.py` ends "S2 PASS" (exact rationals: D = 4 kernel zero off coincidence, closed form at three regulators, D = 5/6 flagged non-zero) and `s5_offplane.py` runs clean. Lane's own code re-executed; `s3_limit.py` (the smeared limit, about 5 minutes) and `s4_moments.py` were not re-run by the lead.
- Independence: this lane was barred from the first lane's report and scripts until its own result was fixed, and used a different route (TT mode sum, exact kernel, smeared limit) from the first lane's (de Donder propagator, harmonic projection, sampled points). Both give 12G delta/(u1 u2 (u1-u2-i0)^4) in D = 4.
- Result on the first lane (ss-tension): its headline stands; two of its side statements are wrong and are withdrawn here: (1) its erratum that off-diagonal correlators vanish in every D (they do not in D = 5, 6; its configuration was degenerate); (2) its unchecked item that (B.26) diverges logarithmically at the edge (it converges for f supported inside the diamond). Also overstated there: that the authors "missed" D = 4.
- Standing: in D = 4 the published graviton state (B.26, sign as in B.13) is the Minkowski vacuum on the gamma[f] algebra, by two independent computations; in D = 5 and 6 the paper's sentence is right. Still a question for the authors, not a correction: the four-moment lemma is exact on 11 sample sources, not proved for all f, and some supporting literature is recalled, not read.
- The fetched PDF of the paper was moved out of the repo by the lead and is not committed.
- Grade: 385 of 420. Genuinely independent route, a method that discriminates dimensions, caught two errors in the lane it was checking, numbered logs, declared amendments; marked down for the sampled lemma and unread references.

The lane's report follows unchanged.

**SS-SECOND report (2026-11-15)**

Tags: [P] proved by hand this session, [N] exact or numeric this session, [S1] source read first-hand, [S2] second-hand or from memory, [B] belief.

## 1. Verdict

**D = 4: CONFIRMED, high certainty.** On the algebra the paper itself defines (γ[f], f symmetric, divergence-free, supported in the diamond), the state (B.26) with its sign corrected is the Minkowski vacuum, coefficient included. The one-line reason: the vacuum two-point function of the Weyl component Ψ₀ = R(k,m,k,m) on the cone is exactly ultralocal, 12G·δ_{S²}/(u₁u₂(u₁−u₂−i0)⁴), which is ∂₁∂₂ of their kernel because Ψ₀ = −S″/(2u).

**D = 5 and 6: the paper's sentence is right.** The same exact kernel is non-zero between different generators, so the vacuum is not ω₀^GR there. The method therefore discriminates dimensions, which the first lane's could not.

Two things in the first lane's report are wrong, though neither touches its headline (section 4). The first lane's conclusion survived an independent attack; its erratum did not, which is the only time today a correction needed correcting.

## 2. Their statements, exactly [S1, page images of printed pp. 2–3, 7–11, 23–33]

- **Version.** arXiv:2607.20678v1 is the only version (submitted 22 Jul 2026). INSPIRE lists 3 citing papers (2610.37904, 2610.20398, 2609.01706); none is a comment or erratum by title, and I did not read them.
- **Scalar (3.4).** ⟨∂_uΦ ∂_uΦ⟩ = −(1/4π)δ/(u₁₂−iε)², with Φ = u^{(D−2)/2}φ.
- **p. 24.** "the free graviton or the free Maxwell field (which is non-conformal in D ≠ 4)". The parenthesis attaches to Maxwell only.
- **p. 26.** "except in conformally invariant settings (like the Maxwell field in D = 4), the state 'ω₀' … will not be the vacuum state."
- **Algebra (B.16).** "A_GR is the algebra of polynomials in γ[f]", f symmetric and divergence-free. So it is gauge-invariant by construction, with no soft, edge or tip modes in it. The diamond is assumed "invariantly defined" (p. 31).
- **Variables.** Ω = (1/32πG)∫u^{D−6}(s₁∂_us₂ − s₂∂_us₁) (B.23); S = u^{(D−6)/2}s (B.24); γ[f] = −(1/16πG)∫F^{AB}∂_uS_AB (B.25), with F from Ef "in Gaussian null gauge". Footnote 12: s = O(u²) at the tip.
- **(B.26).** Printed as +(1/(32π²G))∫F₁F₂ δ/(u₁−u₂−i0⁺)².
- **p. 32, verbatim.** "It is straightforward to show that ω₀^GR defines a state on A^GR(D) which has geometric modular flow on C⁻. Since the linearized gravitational field is not conformally invariant, ω₀^GR is not the restriction of the vacuum state to D." The D is calligraphic (the diamond); no dimension restriction is stated.

**The sign.** The printed plus is a misprint, and nothing else depends on it.
- [P] 1/(u−i0)² = −∫₀^∞ ω e^{−iωu}dω, so ∫∫F F/(u₁₂−i0)² = −∫ω|F̂|² ≤ 0. The printed form gives ω(γ[f]²) < 0.
- [N] Numerically −3.904 against −3.904 from the Fourier side.
- [P] The magnitude 1/(32π²G) is forced by (B.23) and (B.25): (1/16πG)²·8G. There is no free normalisation.

**Can the sentence be true at D = 4 under their definitions?** No. The only reading I found that rescues it is a larger algebra of all ∫F ∂_uS with arbitrary F. That contains tip-anchored shear observables, which are not gauge-invariant and have divergent vacuum variance (1/ε² at the tip) [P]. It is not the algebra they define.

## 3. Derivation, computations, controls

**Route A: TT mode sum [P, then N].**
- The action (1/16πG)√g R gives ⟨γγ⟩_TT = 32πG·Π_ijkl × scalar mode sum.
- For a plane wave, R(l,m,l,m) = ½e(v,v) with v = (l·k)m − (m·k)l. With m complex null this gives, exactly and regulated: ⟨Ψ(1)Ψ(2)⟩ = 8πG·T^{ab}T^{cd}∂_a∂_b∂_c∂_dW_D, where T^{ab}k_ak_b = v₁·v₂.
- On the cone, T·x·x = (cL−ab)(−u₁u₂L + iεu₁₂) − ε²c, with c = m₁·m₂′, L = l₁·l₂, a = m₁·l₂, b = l₁·m₂′.
- In D = 4 with opposite helicities, cL−ab ≡ 0. So the kernel is (768G/π)·ε⁴(m₁·m̄₂)²/σ⁵, with σ = 4u₁u₂sin²(θ/2) + 2iεu₁₂ + ε².
- This is a nascent delta on the sphere. Its angular integral is π/(64u₁u₂(u₁₂−iε/2)⁴), so the coefficient is 768/64 = 12G.

**Checks [N].**

| Script | What it shows |
|---|---|
| `s1` | R_{uÂuB̂} = −½u⁻²∂_u(u²∂_uŝ) on the cone, from a brute-force linearised Riemann tensor. |
| `s2` | Exact rationals: the D = 4 kernel is identically 0 off coincidence; the closed form above holds at three ε; the same-helicity kernel is 0. |
| `s3` run06 | Smeared in (u₁,u₂) on an 8001-point grid, four test pairs including one that is O(u²) at the tip, ε = 0.04 down to 0.005. |

In `s3` the error against the same-ε comparator is 1.9e-3 to 5.9e-3 at the smallest ε, the error against the ε = 0 target halves with ε, and the Richardson extrapolant is within 3.8e-4 to 6.0e-3.

**Controls, all pass [N].**
- Scalar, D = 4 and 6: −1/(4π).
- Maxwell D = 4: −1/(4π), ultralocal.
- Maxwell D = 5, 6: non-ultralocal, with the hand formula −2p(p−1)(cL−ab)σ^{−p−1}, which is visibly proportional to D−4.
- Plants caught: 16πG normalisation (relative error 0.50), and the D = 5 kernel fed to the D = 4 test.

**D = 5, 6 graviton [N, exact].** The kernel is non-zero at ε = 0 between distinct generators (for example 12500/49 and 160000/147 in units of σ^{−p}).

**Route A′: commutator function [N exact on samples, P for consequences].**
- For a curvature observable at an interior point, ∫u^kΨ₀[EO]du = 0 for k = 0, 1, 2, 3 on every generator, and k = 4 is non-zero.
- Tested on 3 hand-chosen cases × 12 generators plus 8 random cases × 2.
- So F is compactly supported inside (0,1) with ∫F du = 0, (B.25)/(B.26) converge as written, and the integration by parts to the Ψ₀ form has no boundary terms.

**Modular flow [P, using standard KMS and Möbius facts].**
- The weight-1 flow F → F∘ψ⁻¹ on ∂_uS does not preserve ∫F du = 0, so it is not an automorphism of A^GR(D).
- The vacuum flow there is the weight-2 Möbius flow on ∂_u²S, i.e. conformal on the curvature.
- On the potential it reads δS = 2π[u(1−u)S′ + 2∫₀ᵘS], up to a + bu.
- For Maxwell the invariant is a first derivative, so no such issue arises.

**Declared amendments.**
- My first Fourier comparator had the wrong complex conjugate. It was caught because the controls failed identically.
- The P2 pass line was made precise and ε was pushed lower after the first run.
- `s3` run05 has NaN in two rows (series overflow); run06 is the record.
- `s4` came out stronger than preregistered: M₀ = 0, not merely l = 2.

## 4. Comparison with the first lane's REPORT.md (read after my result was fixed)

**Agreements.**
- Ψ₀ = −S″/(2u), the constant 12G, and the vanishing unbarred correlator.
- Controls at −1/(4π).
- The sign of (B.26) is a misprint.
- ∫F du = 0 for curvature observables.
- The flow on S (theirs omits my 2π), and δΨ₀ = NΨ₀′ + (3−5u)Ψ₀, which follows from my weight-2 statement.

**Wrong 1: their T4 erratum.** It says the off-diagonal correlator "vanishes in every D" and that D ≠ 4 was therefore untestable.
- Their configuration keeps every vector in a 4-dimensional subspace, with the displacement inside the polarisation plane. There cL−ab = 0 for any propagator power.
- `s5` reproduces their zero in that configuration and gets non-zero for a displacement orthogonal to the polarisation plane, which exists only for D ≥ 5.
- Their original preregistered expectation was right.

**Wrong 2: their unchecked item 3.** It says F carries a tail b(x)·u to the edge, so literal (B.26) diverges logarithmically. For f supported inside the diamond there is no tail (M₀ = M₁ = 0 exactly), and (B.26) converges.

**Overstated.**
- "Missed this one" as a description of the authors: p. 24 shows they class the graviton as non-conformal generically; whether D = 4 was considered is unknown.
- The lead's line "D other than 4 untested" is now closed.

## 5. What may and may not be claimed

**May say.** "In D = 4 (B.26) of arXiv:2607.20678, with the sign of (B.13), is the vacuum on the γ[f] algebra: the vacuum two-point function of R(k,m,k,m̄) on the cone is exactly 12G δ_{S²}/(u₁u₂(u₁₂−i0)⁴), which is ∂₁∂₂ of its kernel since Ψ₀ = −S″/(2u); in D = 5 and 6 it is not ultralocal, so the sentence holds there."

**Must not say.**
- That the paper or its theorems are wrong.
- That the state has no geometric flow.
- "Off-diagonal correlators vanish in every D", or "(B.26) diverges at the edge".
- "Proved" for the moment lemma: it is exact on 11 source configurations, not a symbolic proof for all f.

**Question for the authors.** "For f supported inside the diamond in D = 4, isn't ω₀^GR(γ[f₁]γ[f₂]) the vacuum value, and which flow is meant by 'geometric': the diffeomorphism on ∂_uS does not preserve ∫F du = 0, which every such f satisfies, whereas the weight-2 action on ∂_u²S does?"

## 6. Unchecked, ranked

1. The four-moment lemma for all f: exact on samples, with no symbolic proof.
2. That curvature smearings exhaust γ[f] (compact-support Poincaré lemma for the Calabi complex), and that cone data generate the diamond von Neumann algebra [S2/B].
3. The flow argument rests on continuity of F ↦ ∫F du on the local subspace and on results for derivative-of-current subnets (Guido–Longo–Wiesbrock), which I recalled from memory and did not read.
4. In D ≠ 4, whether ω₀^GR has any geometric flow on A^GR(D), i.e. whether their "straightforward" claim holds.
5. The ε → 0 limit near the tip is a scaling argument plus numerics, not a theorem.
6. Literature: Longo–Morinelli–Preta–Rehren read as abstract only ("Finite helicity Poincaré representations extend to the conformal group") [S1]; Hislop 1988 and Mack 1977 unread.

**Next experiment.** Compute the admissible F-space in D = 6 by the `s4` method and test whether the paper's flow preserves it. That decides whether their graviton flow claim holds anywhere.

## Reproduce

From `record\wave73\ss-second\`, native Windows Python 3.14:
- `python scripts\s1_sachs.py`
- `python scripts\s2_kernel.py`
- `python scripts\s3_limit.py 8001 0.04,0.02,0.01,0.005` (about 5 minutes)
- `python scripts\s4_moments.py`
- `python scripts\s5_offplane.py`

Nothing committed, no report file written, no subagents, no processes killed. The fetched PDF sits inside the lane folder (0.9 MB, uncommitted); decide whether it stays before any commit.

Files are in `record\wave73\ss-second\`:
- PREREG.md (result in section 6, comparison in section 7)
- scripts\common.py
- scripts\s1_sachs.py
- scripts\s2_kernel.py
- scripts\s3_limit.py
- scripts\s4_moments.py
- scripts\s5_offplane.py
- logs\s1_run01.log
- logs\s2_run01.log
- logs\s3_run01_N2001.log to s3_run04_N2001.log
- logs\s3_run05_N8001.log
- logs\s3_run06_N8001.log
- logs\s4_run01.log
- logs\s4_run02.log
- logs\s5_run01.log
- paper\2607.20678v1.pdf
