# REPORT: lane ODD-SIGN (wave 73)

## Lead's check and grade

Sunday 15 November 2026. Small proof lane (Fable): carry the sign-definite argument of the even lanes to the odd-parity tower, so that the first-order lower bound no longer rests on HINGE-PROOF's Lemma B outline and a numerically evaluated Bessel identity. The lane returned the text below; the lead saves it.

**Grade: 388 of 420.** One question, answered with an explicit family and constant, all fifteen pre-fixed lines passing on the first and only run, and an exact list of what it replaces in HINGE-PROOF and what it leaves. Marked down because the six certified inputs are the lane's own ball arithmetic, and the remainder for this family (a proved lower bound on the vector's norm, and the Duhamel-Gronwall corollary) is still open for the odd tower.

**Checked by the lead, by hand. It holds.**
- Step 1: omega(k, D_0 k) = <q, -N q> - <pi, 0> = -<q, N q>, definite, on the q slot.
- Step 2: for e = (q, 0), q_s = cos(s h) q under the vacuum flow.
- Step 7: (1 + cosh rho)/lambda = e^y (1 + e^(-y)/(2 lambda))^2, decreasing in l; Cauchy-Schwarz in the weight lambda N gives <b, lambda N b> >= <b, u>^2/<u, (lambda N)^-1 u>.
- Steps 8 and 9: identical to EVEN-PASSAGE's steps 5 to 7, hand-checked there.
- Step 10: with z = h^(1/2) q/sqrt(2) and q = h e^(-eps h^2)(h u), the spectral weights are p^3 sinh(2 pi p) e^(-2 eps p^2) and p^3 coth(pi p) e^(-2 eps p^2) <= (p^3 + p^2/pi) e^(-2 eps p^2).
- Arithmetic of step 11 recomputed from the six inputs: M_lo = 0.004656, Den <= 4.8985, ||e|| <= 0.2332, c = 2.038e-3. Agrees.
- Not hand-checked: step 5's domain statement, the six ball-arithmetic inputs.

**Re-run by the lead:** `..\..\even-nonzero\venv\Scripts\python.exe odd_sign.py <out.json>` from `scripts\`, 33 s, exit 0, "ALL LINES HELD"; the table equals the lane's (0.2974 rising to 0.3053 against the proved 2.038e-3; sign identity 1.000000000 at every l; both plants caught).

**What it changes.** For the odd tower the first-order lower bound is now a written proof on an explicit family, uniform in l >= 16, with no special-function identity. The Bessel identity is needed only for the value 7 zeta(3)/(16 pi^2) of the per-sector constant, not for the divergence. Still open for the odd tower: the second-order remainder on this family (the even tower's was closed by EVEN-REMAINDER; the same method should apply but has not been run), and the inherited items common to both towers.

---

## Lane report (unchanged)

**ODD-SIGN report (2026-11-15)**

The sign argument carries over to the odd tower with an explicit family and constant. All 15 pre-fixed lines passed on the first and only run.

## 1. Statement

Let u(y) = (1 − ((y−2)/4)²)⁴ on −2 ≤ y ≤ 6 and 0 elsewhere, with y = ρ − ln 2λ, ε = 1.5, τ = π. Take

e_l = (q_l, 0), q_l = e^{−εh_l²} h_l² u.

Then for every l ≥ l₀ = 16: λ‖β̂₁^(l) e_l‖/‖e_l‖ ≥ c = 2.038·10⁻³. Here β̂₁ is BETA-DEF's (unscaled), built from T = ∫₀^τ e^{sX₀}D₀e^{−sX₀}ds with D₀(q, π) = (0, −Nq), i.e. g factored out. On the grid the value is 0.297 rising to 0.305.

## 2. Proof

Notation: z = (h^{1/2}q + i·h^{−1/2}π)/√2, ω = ⟨q, π′⟩ − ⟨π, q′⟩ = 2 Im⟨z, z′⟩, A_l = ‖h_l u‖², G_l = ‖h_l² u‖, V_l = λ²/sinh²ρ, Q_l = ⟨u, (λN_l)⁻¹u⟩.

1. [P] **Sign and slot.** ω(k, D₀k) = ⟨q, −Nq⟩ − ⟨π, 0⟩ = −⟨q, Nq⟩ ≤ 0, so the definite form sits on the q slot. With the brief's D it is −g⟨q, Nq⟩.
2. [P] The vacuum flow is z ↦ e^{−ish}z, which is symplectic. So ω(e, Te) = −m_l(τ), with m_l = ∫₀^τ⟨q_s, N q_s⟩ds, and Im⟨z, Tz⟩ = −m_l/2. For e = (q, 0), q_s = cos(sh)q.
3. [I, BETA-DEF steps 5 and 7, lead-checked; the minimisation re-derived here] β̂₁k(z) = (cP, s·Q̄) with P + Q = Tz, so ‖β̂₁e‖ ≥ ‖(2 sinh 2πh)^{−1/2}Tz‖. T is bounded by Lemma A [I, lead-checked].
4. [P] Cauchy–Schwarz: m_l/2 ≤ |⟨z, Tz⟩| ≤ Den_l·‖β̂₁e‖, where Den_l² = ⟨z, 2 sinh(2πh)z⟩.
5. [P] u is C³ with compact support in ρ > 0 for l ≥ 16, so u ∈ D(h_l²) with h_l²u = −u″ + V_l u. [S, from memory] h_l² is essentially self-adjoint on C_c^∞ because λ² ≥ 3/4.
6. [P] sinh ρ = λe^y(1 − e^{−2y}/4λ²), so V_l decreases to e^{−2y}. Hence A_∞ ≤ A_l ≤ A_16 and G_l ≤ G_∞ + ‖(V_16 − V_∞)u‖ =: G_hi.
7. [P] (1 + cosh ρ)/λ = e^y(1 + e^{−y}/2λ)², decreasing in l, so Q_l ≤ Q_16. Cauchy–Schwarz in the weight λN: ⟨b, λN_l b⟩ ≥ ⟨b, u⟩²/Q_16 for every b.
8. [P] ⟨q_s, u⟩ = ⟨hu, cos(sh)e^{−εh²}hu⟩ ≥ J₀ − s²J₁/2, using cos x ≥ 1 − x²/2. Jensen on the spectral measure of hu (mass A_l, second moment G_l²) gives J₀ ≥ A_∞·exp(−εG_hi²/A_∞) =: J₀_lo. Also J₁ ≤ A_16/(eε) =: J₁_hi.
9. [P] The integrand is non-negative, so λm_l ≥ (8/15)·J₀_lo²·s*/Q_16 =: M_lo, with s* = (2J₀_lo/J₁_hi)^{1/2} < π.
10. [P] Den_l² = ⟨hu, h³sinh(2πh)e^{−2εh²}hu⟩ ≤ S′·A_16, with S′ = sup p³e^{−2εp²}sinh 2πp. And ‖e_l‖² = ⟨z, coth(πh)z⟩ ≤ ½C·A_16, with C = sup(p³ + p²/π)e^{−2εp²}.
11. [C] Ball arithmetic at 100 bits, radii below 1e-26:

| Input | Value | Derived | Value |
|---|---|---|---|
| A_∞ | 0.938464 | J₀_lo | 0.363568 |
| A_16 | 0.938807 | J₁_hi | 0.230245 |
| G_hi | 0.770249 | s* | 1.77710 |
| Q_16 | 26.90624 | M_lo | 0.0046562 |
| S′ | ≤ 25.55923 | Den | ≤ 4.89849 |
| C | ≤ 0.115848 | ‖e‖ | ≤ 0.23319 |

   These give c = M_lo/(2·Den·‖e‖) = 2.038·10⁻³. With the l = ∞ constants the bound is 2.094·10⁻³.

Domain: z_l = h^{5/2}e^{−εh²}u/√2 lies in the domain of coth^{3/2}(πh) and of e^{πh}, which is HINGE-PROOF's hypothesis.

## 3. Check

| l | λ‖β̂₁e_l‖/‖e_l‖ [N] | uniform floor [C] | floor at this l [C] | λm_l | Den_l | ‖e_l‖ |
|---|---|---|---|---|---|---|
| 16 | 0.297406 | 2.038e-3 | 2.085e-3 | 0.022745 | 2.32772 | 0.17460 |
| 64 | 0.303293 | 2.038e-3 | 2.091e-3 | 0.023250 | 2.32695 | 0.17457 |
| 256 | 0.304791 | 2.038e-3 | 2.093e-3 | 0.023383 | 2.32690 | 0.17456 |
| 1024 | 0.305167 | 2.038e-3 | 2.093e-3 | 0.023416 | 2.32690 | 0.17456 |
| 4096 | 0.305262 | 2.038e-3 | 2.094e-3 | 0.023425 | 2.32690 | 0.17456 |
| limit | 0.305293 | | 2.094e-3 | 0.023427 | 2.32690 | 0.17456 |

- **Floor:** holds at every l; the grid value is about 146 to 150 times the floor.
- **Sign identity:** −2 Im⟨z, Tz⟩/m = 1.000000000 at every l.
- **Intermediate inequalities:** all held; the pointwise integrand is at least 3.86 times its floor.
- **Approach to the limit:** λ·|r_l/r_∞ − 1| is 0.426, 0.423, 0.422, 0.421, 0.421.
- **Grid stability:** halving dy moves r by 9.3e-6; doubling the box by 7e-12.
- **Plant A** (factor λ dropped): values fall from 1.8e-2 to 7.5e-5, below the floor from l = 256 on. Caught.
- **Plant B** (wrong slot: even-parity sign on the z̄ part of T, l = 64): the sign identity is off by 18.3% of m, and r shifts 1.2%. Caught.
- **Failed lines:** none.

## 4. Re-run

From `record\wave73\odd-sign\scripts`:

`..\..\even-nonzero\venv\Scripts\python.exe odd_sign.py ..\out\odd_sign.json`

It takes 25 s and exits 0. Logs are `logs\01_prefixed_pass_lines.log` and `logs\02_odd_sign_run1.log`. Nothing committed, no REPORT.md.

## 5. Weakest step

Step 3: the bound is on BETA-DEF's formula for β̂₁, whose link to implementability rests on two theorems quoted there from memory; within this lane the only non-[P] inputs are the six ball-arithmetic constants, not independently re-run.

## 6. What this replaces in HINGE-PROOF, and what it leaves

**Replaces:**
- Lemma B entirely: the strong-resolvent outline, the weak-convergence bullet and the appeal to J(t)/2 > 0.
- Item 4 of its Section 6 ("Lemma B written line by line, and HINGE's Bessel identity proved") as a prerequisite for the divergence. The Bessel identity is now needed only for the value of the constant 7ζ(3)/(16π²), not for positivity.
- The "[P-outline]" half of step 7 of its table, for the first-order lower bound.

**Leaves open:**
- The remainder for this family [U]: sup_l C(ψ_l)/‖e_l‖ finite needs a proved lower bound on ‖e_l‖ (0.1746 on the grid, [N] only) and the Duhamel–Gronwall corollary, which is not lead-checked.
- The sharp per-sector constant (table step 6) stays numeric.
- Connes, Shale and Araki duality are from memory; physics gaps of steps 1, 9 and 12 are untouched.
- Normalisation: HINGE-PROOF's β̂₁ carries a factor λ; this lane uses BETA-DEF's unscaled one.
