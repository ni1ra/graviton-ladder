# REPORT: lane EVEN-NONZERO (wave 73)

## Lead's check and grade

Sunday 15 November 2026. Small proof lane (Fable): in the large-l limit model, is the first-order term of the even-parity flow non-zero on some vector? The lane returned the text below; the lead saves it.

**Grade: 388 of 420.** The question is answered by a short analytic argument that is stronger than asked (non-zero on every non-zero vector, every t > 0) and needs no special-function identity. The lane reported its own failed pass line and separated the theorem from the certified numbers. Marked down because the explicit form of the first-order term is the lane's reconstruction from HINGE-PROOF's code rather than a written definition, and because the certified enclosure rests on a spectral density quoted from memory (step 9).

**Checked by the lead, by hand: steps 1, 4, 5 and 6, which give the theorem. They hold, given the lane's form of the first-order term.**
- Step 1: with T z = A z + B z-bar, P = B(c^2 z-bar) - A(s^2 z), Q = A(c^2 z) - B(s^2 z-bar) and c^2 - s^2 = 1, P + Q = T z. Since c and s are positive functions of h, T e != 0 forces the first-order term to be non-zero.
- Step 4: omega(k, D k) = <q, 0> - <pi, W pi> = -<pi, W pi>; the vacuum flow preserves omega, so omega(e, T e) = -integral over [0, tau] of <pi_s, W pi_s>.
- Step 5: W_inf is positive and injective, so the integral vanishes only if pi_s = 0 on the whole interval, which forces e = 0.
- Not hand-checked: steps 2, 3, 7 and 8 (positivity and the bounds on W_inf, used for the enclosure), step 9 (quoted), the ball arithmetic of step 10.

**Re-run by the lead:** `..\venv\Scripts\python.exe even_nonzero.py <out.json>` from `scripts\`, 20 s, exit 1 with exactly one failed line, N7, as the lane reports. Numbers equal the lane's: R_inf = 0.01354748, m = 0.06080047, Im<z, T z> = -0.03040023 = -m/2, finite-l values at l = 64, 256, 1024 approaching with lambda times the gap = 1.002, 1.005, 1.005; plant A shifts R by 35.7%, plant B by 1.3%. The venv (python-flint 0.9.0) is in the lane folder and is not committed.

**On the failed line N7.** The lane required the "coth dropped" plant to move the number by more than 20%; it moved it by 1.3%. The plant is detected (about 200 times the grid tolerance) but the pre-fixed threshold was not met. This weakens the check as a test of the coth factor; it does not touch the theorem, which uses no numerics.

**Two things worth acting on.**
- The lane's remark (tagged unverified by it): steps 1, 4 and 5 use only that the extra term has a definite sign, so the same argument should replace the numerically evaluated Bessel identity in the odd Lemma B. That would be a separate small lane.
- R_inf = 0.013547 is close to the even per-sector constant (about 0.0136, measured 0.013586). The lane says these are different objects and that it has no explanation. Do not present the closeness as a confirmation of anything.

**Where the even half stands after this lane.** Proved and lead-checked: the uniform bound (EVEN-PROOF), the limit operator (EVEN-LIMIT, statement a), non-vanishing of the first-order term in the limit model (here). Still missing for the even divergence as a theorem: the passage from finite l to the limit (strong resolvent convergence of h_l, so that the liminf is positive), the remainder constant checked for the chosen vector, and a written definition of the first-order term to replace the reconstruction. The paper must still say "not proved" for the even divergence.

---

## Lane report (unchanged)

**EVEN-NONZERO report (2026-11-15)**

The statement holds in the limit model, by an analytic argument that needs no Bessel identity: β̂₁e ≠ 0 for every non-zero e in H and every t > 0. One pre-fixed line failed (N7, plant B).

## 1. Definition and statement

From HINGE-PROOF Section 3, lines 63–65: β = ½(Ŝ + iŜi) with Ŝ = S ⊕ 1 on H + H′, and i·k(ψ) = k(coth(πh)·iψ) + k′(−csch(πh)·iψ). The first order is β̂₁ = ½(T̂ + iT̂i), T̂ = T ⊕ 0, T = B₁(τ) = ∫₀^τ e^{sX₀} D e^{−sX₀} ds, τ = 2πt.

I made this concrete from `hs_beta` in `u1_uniformity.py`: z = (h^{1/2}q + i·h^{−1/2}π)/√2, k(z) = (cz, s·z̄), c² = 1/(1−e^{−2πh}), s² = c² − 1. Writing Tz = Az + Bz̄:

β̂₁k(z) = (c·P, s·Q̄), with P = B(c²z̄) − A(s²z), Q = A(c²z) − B(s²z̄).

This reproduces line 63 and the four terms of `hs_beta`.

**Theorem.** In the limit model with D(q, π) = (W_inf π, 0): β̂₁e ≠ 0 for every e ≠ 0 in H and every τ > 0.

**Explicit instance.** π₀ = x·e^{−x}, e = (0, e^{−h²/4}π₀), t = 0.5, η = (iz/c, conj(iz)/s). Then Re⟨η, β̂₁e⟩ lies in [−0.098012, −0.014876], and ‖β̂₁e‖/‖e‖ ≥ 5.11·10⁻⁴.

## 2. Argument

1. [P] P + Q = Tz, since c² − s² = 1. So Te ≠ 0 implies β̂₁e ≠ 0, and ‖β̂₁e‖ ≥ ‖(2 sinh 2πh)^{−1/2}Tz‖.
2. [P] W_inf is positive and injective: ⟨f, Wf⟩ = ∫(w′² + w²)dx, where −w″ + w = f/x, w(0) = 0.
3. [P] ‖h_inf W_inf‖ ≤ 1. From h²w = xf − xw_x and ∫ww_x dx = 0: ⟨xw, f⟩ = ‖hw‖², and x² ≤ h². So D and T are bounded on H, as in EVEN-PROOF.
4. [P] With ω = ⟨q, π′⟩ − ⟨π, q′⟩ = 2 Im⟨z, z′⟩: ω(k, Dk) = −⟨π, Wπ⟩. The vacuum flow is symplectic, so ω(e, Te) = −m(τ), where m(τ) = ∫₀^τ⟨π_s, Wπ_s⟩ds and (q_s, π_s) = e^{−sX₀}e.
5. [P] m(τ) > 0 unless π_s ≡ 0 on [0, τ]. That forces sin(sh)·z = 0 for all s, so e = 0. With step 1 this gives the theorem.
6. [P] Matrix element: Re⟨η, β̂₁e⟩ = Im⟨z, Tz⟩ = −m(τ)/2. Here η is in K because e^{πh}z is in L².
7. [P] ‖W_inf‖ ≤ ½, from ∫w²/x ≤ ½∫(w′² + w²) (complete the square with x·e^{−x}). Also Wπ₀ = π₀/2 exactly.
8. [P] For a = φ(h)π₀: ⟨a, Wa⟩ ≥ 2⟨a, π₀⟩² and |⟨a_s, Wa_s⟩ − ⟨a, Wa⟩| ≤ ½s²‖h²a‖‖a‖. m is non-decreasing, so m(π) ≥ ⅔·w_lo·s*.
9. [I, from memory] The Kontorovich–Lebedev spectral density of π₀ is 2p³/sinh(πp) dp (GR 6.611.3 plus Plancherel). The moments 1/4, 1/2, 17/8 match my x-space integrals, but the density itself is not proved here.
10. [C, given 9] Ball arithmetic at 120 bits, with head and tail bounded by hand: m(π) ∈ [0.029752, 0.196024]; ‖β̂₁e‖/‖e‖ ≥ 5.11·10⁻⁴.
11. [N] m(π) = 0.060800 and R = ‖β̂₁e‖²/‖e‖² = 0.013547.

Remark [U, not written out]: steps 1, 4 and 5 use only that the charge is sign-definite, so they should also replace the Bessel identity in the odd Lemma B.

## 3. Check

| | R | m | \|R_l/R_inf − 1\| | λ × that |
|---|---|---|---|---|
| limit | 0.01354748 | 0.06080047 | | |
| l = 64 | 0.01333702 | 0.06024470 | 1.554e-2 | 1.002 |
| l = 256 | 0.01349442 | 0.06066235 | 3.917e-3 | 1.005 |
| l = 1024 | 0.01353419 | 0.06076600 | 9.813e-4 | 1.005 |

- **Grid stability:** halving dy changes R by 6.0e-5; doubling the box changes it by 9e-12.
- **Identity check:** Im⟨z, Tz⟩ = −0.03040023, equal to −m/2 as step 6 predicts.
- **Plant A** (e^{−y} in place of W_inf): R = 0.018386, a 35.7% shift. Caught.
- **Plant B** (coth dropped): R = 0.013370, a 1.3% shift. **N7 FAILED as pre-fixed**: I asked for more than 20%. The shift is still about 200 times the grid tolerance.
- C1–C3 and N1–N6 passed.
- R_inf is close to hinge-e's 0.0136, but they are different objects and I have no explanation for the closeness.

## 4. Re-run

From `record\wave73\even-nonzero\scripts`:

`..\venv\Scripts\python.exe even_nonzero.py ..\out\even_nonzero.json`

It takes 21 s and exits 1 because of N7. The venv is in the lane folder (python-flint 0.9.0). Logs are in `logs\`: `01_prefixed_pass_lines.log`, `02_even_nonzero_run1.log`, `03_even_nonzero_run2.log`, `04_note_on_runs.log`. Run 1 crashed only when writing the JSON; run 2 is the same numerics with that cast fixed.

## 5. Weakest step

The form of β̂₁ on H + H′ is my reconstruction from `hs_beta` and line 63, not a written definition. The theorem needs only that; the certified numbers also rest on step 9.

## 6. Still missing for the even divergence as a theorem

- The passage from finite l to the limit: liminf ‖β̂₁^(l)e_l‖ > 0. By step 8 it is enough that ⟨a, λW_l a⟩, ‖a‖, ‖h_l²a‖, ‖e_l‖ and ‖(2 sinh 2πh_l)^{1/2}z_l‖ converge, which needs strong resolvent convergence of h_l. The table shows a clean 1/λ approach [N].
- HINGE-PROOF's remainder constant C(ψ) checked for this e.
- A proof of step 9, if the numbers are quoted.
- EVEN-PROOF's dictionary (its step 1) and the physics gaps of HINGE-PROOF steps 1, 9 and 12.
