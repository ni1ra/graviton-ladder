# REPORT: lane ALL-T (wave 73)

## Lead's check and grade

Sunday 15 November 2026. Small proof lane (Fable): extend the first-order floor and the remainder from the single modular time t = 1/2 to every non-zero time, both signs, both parities. The lane returned the text below; the lead saves it. Mid-run the lead forwarded the referee's findings (state the orientation of the vacuum flow; a single-time test can be blind in a toy model, so say why the floor cannot vanish at isolated times); the lane added lines for them after its first run and says so.

**Grade: 384 of 420.** The extension is clean (the bound runs through a monotone quantity, so it cannot return to zero the way the referee's toy cocycle can), the constants are closed-form in tau, it corrected the lead's brief (the threshold l_1 does not blow up as tau goes to 0; it falls to 16), and it reports a failed pre-fixed line and two script bugs in its first run with the predictions it logged before fixing them. Marked down because the bound is tested non-vacuously at only one (tau, l), the inherited formula for beta and the perturbed group were not re-derived, and four pass lines were added after run 1.

**Checked by the lead, by hand. It holds.**
- Step 1: with Theta(q, pi) = (q, -pi), Theta X_0 Theta = -X_0 and Theta D Theta = -D for both D(q, pi) = (W pi, 0) and (0, -N q), and Theta is an isometry of H, so S(-tau) = Theta S(tau) Theta.
- Step 3: W >= 0 and N >= 0 make m(tau) non-decreasing for tau > 0.
- Step 5: integral from 0 to t of (J_0 - s^2 J_1/2)^2 ds = J_0^2 t - J_0 J_1 t^3/3 + J_1^2 t^5/20; with s*^2 = 2 J_0/J_1 this is (8/15) J_0^2 s* Phi(t/s*), Phi(sigma) = (15 sigma - 10 sigma^3 + 3 sigma^5)/8, Phi(1) = 1.
- Slope at 0: c_pi * 15/(8 s*) = 2.1193e-3 * 1.875/1.7771 = 2.236e-3.
- Step 10: K x^2 e^x/2 with x^2 = 2 g^2 tau^2/lambda^2 is at most c(tau)|g|/(2 lambda) exactly when lambda >= 2 e^x K tau^2 |g|/c(tau).
- Not hand-checked: step 6's re-derivation, step 7's constants, the bound eta on x.

**Re-run by the lead:** `..\..\even-nonzero\venv\Scripts\python.exe -B all_t.py <out.json>` from `scripts\`, 50 s, exit 1 with exactly one failed line, C2, as the lane reports; all floor, remainder, symmetry and monotonicity lines pass and the plant is caught.

**On the failed line C2.** The lane predicted the odd threshold l_1(1, pi) would be 15237 or 15238 and computed 15236, because it used a slightly sharper factor than the record's rounded 1.001. A wrong prediction of an integer by one, not a violated inequality. ODD-REMAINDER's published 15237 (with the rounded factor) remains valid as a threshold.

**What it changes.** "Only modular time 1/2 is covered" is closed: for every t != 0 and g != 0 the cocycle has a non-Hilbert-Schmidt antilinear part over the tower, with an explicit floor that is positive for t != 0 and linear in |t| near 0. The paper may use the lane's sentence in its section 6. This also answers the referee's toy case for these operators: the bound cannot vanish at isolated times.

---

## Lane report (unchanged)

**ALL-T report (15 Nov 2026)**

The floor and the remainder both extend to every τ ≠ 0, either sign, both parities; no range of τ fails. One pre-fixed line failed (C2, my prediction of a ceiling).

## 1. Statement

Orientation: σ^ω_t(W(k)) = W(U⁰(2πt)k), so S(τ) = U^g(τ)U⁰(−τ) with τ = 2πt, the same definition for τ < 0. T(τ) = ∫₀^τ e^{sX₀}De^{−sX₀}ds is the oriented integral.

For every l ≥ 16, real g and τ ≠ 0, with x = √2|g||τ|/λ:

- **Floor:** λ‖β̂₁(τ)e_l‖/‖e_l‖ ≥ c(τ) = c_π·Φ(min(|τ|, s*)/s*), where Φ(σ) = (15σ − 10σ³ + 3σ⁵)/8 and s* = 1.777101.
- **Constants:** even c_π = 2.1193·10⁻³, K = 1.5693 (even-second's enclosures); odd c_π = 2.0381·10⁻³, K = 1.5716 (odd-sign's own).
- **Shape of c(τ):** equal to c_π for |τ| ≥ s*, strictly increasing below, with no zero except τ = 0. It vanishes linearly: c(τ) ≈ 2.2361·10⁻³|τ| (even), 2.1503·10⁻³|τ| (odd).
- **Remainder:** ‖βe_l‖/‖e_l‖ ≥ c(τ)|g|/λ − K(eˣ − 1 − x).
- **Half floor:** the right side is at least c(τ)|g|/(2λ) for l ≥ l₁(g, τ) = max(16, ⌈2e^η K τ²|g|/c(τ)⌉), η = 15c_π/(8√2·s*·K) ≈ 1.0·10⁻³.
- **Conclusion:** Σ_{l≥l₁}(2l+1)(c(τ)g/2λ)² = ∞ for every τ ≠ 0 and g ≠ 0.

**l₁ does not blow up as τ → 0**, contrary to the brief: the remainder is O(τ²) and the floor O(τ), so l₁ ≈ 1405|τ||g| (even) and falls to 16. It grows like 1482τ²|g| (even) and 1544τ²|g| (odd) as τ → ∞. What degenerates at τ → 0 is the coefficient c(τ)² of the divergent sum. Neither matters, because for fixed τ the sum is a tail of the harmonic series.

## 2. Proof

1. [P] Time reversal Θ(q, π) = (q, −π) is isometric on H, with ΘX₀Θ = −X₀ and ΘDΘ = −D for both D(q,π) = (Wπ, 0) and (0, −Nq). So S(−τ) = ΘS(τ)Θ and T(−τ) = ΘT(τ)Θ, and Θe_l = ∓e_l. Every bound at τ transfers to −τ.
2. [P] Independently of step 1: ω(e, T(τ)e) = ∫₀^τ ω(U⁰(−s)e, DU⁰(−s)e)ds = −m(τ), with m(τ) = ∫₀^τ⟨a_s, Wa_s⟩ds (odd: ⟨q_s, Nq_s⟩) and a_s = cos(sh)a even in s. So m(−τ) = −m(τ) and |Im⟨z, Tz⟩| = m(|τ|)/2.
3. [P] W ≥ 0 and N ≥ 0, so m is non-decreasing on τ > 0. This is why the floor cannot vanish at isolated τ, unlike the referee's toy case where S itself returns to 1: the bound runs through a monotone quantity.
4. [I, even-passage steps 4–6 and odd-sign steps 7–8, lead-checked] The integrand is at least (J₀_lo − s²J₁_hi/2)²/Q₁₆ for s ≤ s*.
5. [P] Integrating to t = min(|τ|, s*): λm(|τ|) ≥ [J₀²t − J₀J₁t³/3 + J₁²t⁵/20]/Q₁₆ = M_lo·Φ(t/s*). For |τ| > s* step 3 keeps the value at s*.
6. [P, re-derived] ‖β̂₁e‖ ≥ ‖(2 sinh 2πh)^{−1/2}Tz‖, since c²|P|² + s²|Q|² ≥ |P+Q|²/(c⁻² + s⁻²) and c⁻² + s⁻² = 2 sinh 2πh. Then m/2 ≤ Den·‖β̂₁e‖ with Den² = ⟨z, 2 sinh(2πh)z⟩. Den and ‖e‖ contain no τ.
7. [C] Constants: lanes' ball arithmetic re-run; even chain re-derived from even-second's enclosures, agreeing with that lane's c and C₀ to 1·10⁻¹⁵.
8. [P] Remainder: R₂ = g∫₀^τ(U(r) − U₀(r))DU₀(−r)dr and ‖U(r)‖ ≤ e^{d|r|} for both signs (Gronwall on ±generator), so ‖R₂‖_H ≤ eˣ − 1 − x. [I] ‖D‖ ≤ √2/λ and the bounded perturbation theorem.
9. [I, beta-def/assembly] ‖βe − gβ̂₁e‖ ≤ ½‖R₂‖(‖e‖ + ‖coth^{3/2}ψ‖); no τ enters, so K is unchanged.
10. [P] l₁: eˣ − 1 − x ≤ x²eˣ/2 = g²τ²eˣ/λ². If λ ≥ 2e^ηKτ²|g|/c(τ) then x ≤ c(τ)/(√2K|τ|e^η) ≤ η, because c(τ)/|τ| ≤ 15c_π/(8s*). Hence K(eˣ − 1 − x) ≤ c(τ)|g|/(2λ).
11. [I, assembly B.9] The (2l+1) multiplicity sum.

## 3. Check (run 2, 51 s)

Grid r = λ‖β̂₁e_l‖/‖e_l‖ [N] beside c(τ) [C]. The even c(τ) column is the first code's; with even-second's it is 0.43% higher and every cell still holds.

| parity | τ | c(τ) | l=16 | l=256 | l=4096 | min r/c |
|---|---|---|---|---|---|---|
| even | 0.05 | 1.1127e-4 | 0.006129 | 0.006349 | 0.006364 | 55.1 |
| even | 0.5 | 1.0559e-3 | 0.059278 | 0.061365 | 0.061503 | 56.1 |
| even | 1.7771 | 2.1102e-3 | 0.143955 | 0.147943 | 0.148199 | 68.2 |
| even | π | 2.1102e-3 | 0.121788 | 0.124501 | 0.124673 | 57.7 |
| even | 2π | 2.1102e-3 | 0.104541 | 0.107084 | 0.107246 | 49.5 |
| even | 20 | 2.1102e-3 | 0.104922 | 0.107346 | 0.107500 | 49.7 |
| odd | 0.05 | 1.0746e-4 | 0.010720 | 0.011105 | 0.011130 | 99.8 |
| odd | 0.5 | 1.0198e-3 | 0.104420 | 0.108095 | 0.108336 | 102.4 |
| odd | 1.7771 | 2.0381e-3 | 0.277999 | 0.286145 | 0.286670 | 136.4 |
| odd | π | 2.0381e-3 | 0.297406 | 0.304791 | 0.305262 | 145.9 |
| odd | 2π | 2.0381e-3 | 0.256158 | 0.262363 | 0.262758 | 125.7 |
| odd | 20 | 2.0381e-3 | 0.243317 | 0.249271 | 0.249650 | 119.4 |

Remainder at g = 1, l = 4096:

| parity | τ | actual remainder | proved K(eˣ−1−x) | actual ‖βe‖/‖e‖ | proved lower bound | l₁(1, τ) |
|---|---|---|---|---|---|---|
| even | 0.5 | 1.137e-10 | 2.338e-8 | 1.5013e-5 | +2.355e-7 (half floor 1.294e-7) | 741 |
| even | 2π | 3.654e-9 | 3.694e-6 | 2.6178e-5 | −3.18e-6 (vacuous) | 58524 |
| odd | 0.5 | 1.918e-10 | 2.341e-8 | 2.6446e-5 | +2.255e-7 (half floor 1.245e-7) | 772 |
| odd | 2π | 1.239e-8 | 3.700e-6 | 6.4136e-5 | −3.20e-6 (vacuous) | 60943 |

- **l₁(1, π):** 14631 (even, even-second's constants), 14714 (even, first code), 15236 (odd).
- **Operator norm:** ‖R₂‖_H is 2.9e-9 against 1.49e-8 at τ = 0.5, and 6.9e-8 against 2.35e-6 at 2π.
- **Sign symmetry:** r(−τ) = r(τ) to every printed digit and m(−τ) = −m(τ), at l = 256, τ = 0.5 and 2π.
- **Monotonicity:** λm is strictly increasing along all seven τ. r itself is not monotone; the floor only needs m.
- **Plant** (c continued as the constant c_π below s*): caught on the m-line at τ = 0.05 in all six cells (0.00125 against the planted 0.00750), and on the final floor at the plant-only τ = 0.002 (r = 2.5e-4 against 2.11e-3). As predicted, it is not caught on the final floor at τ = 0.05.
- **FAILED as pre-fixed, C2:** I required odd l₁(1, π) in [15237, 15238] and got 15236, because my per-parity e^η = 1.000968 is below the record's blanket 1.001. A wrong prediction, not a violated inequality.
- **Run 1** failed F2, F3 and F7 from two bugs in my script (an overflow in the even Den, a missing λ on Im⟨z,Tz⟩). Both were fixed with predictions logged first, and no r digit moved. Lines E2_C, E2_F, E2_R and MONO were added after run 1 on the lead's update; they pass.

## 4. Re-run

From `record\wave73\all-t\scripts`:

`..\..\even-nonzero\venv\Scripts\python.exe -B all_t.py ..\out\all_t.json` (51 s, exits 1 because of C2)

Logs are `logs\01` to `05`. Nothing was written outside `all-t\`; no commit, no REPORT.md.

## 5. Weakest step

Step 9 with step 8's inherited half: the β formula and the perturbed group are imported, and I checked only that no τ enters them, not their proofs; the bound is also tested non-vacuously only at τ = 0.5, l = 4096.

## 6. Sentence for the paper

"For every t ≠ 0 and every g ≠ 0 the cocycle U^g(2πt)U⁰(−2πt) has a non-Hilbert–Schmidt antilinear part over the tower, with per-sector floor c(2πt)|g|/(2λ) for l ≥ l₁(g, 2πt), where c is explicit, positive for t ≠ 0 and linear in |t| near 0; hence the flow agrees with the modular group of no faithful normal semifinite weight on M at any single nonzero time."
