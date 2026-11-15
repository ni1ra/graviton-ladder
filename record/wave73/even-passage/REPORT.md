# REPORT: lane EVEN-PASSAGE (wave 73)

## Lead's check and grade

Sunday 15 November 2026. Small proof lane (Fable): carry the non-vanishing of the even-parity first-order term from the limit model to finite l, with a lower bound uniform in l. The lane returned the text below; the lead saves it.

**Grade: 390 of 420.** It found a shorter route than the one briefed (no resolvent convergence, no use of EVEN-LIMIT), gave an explicit family and an explicit constant valid for every l >= 16, reported its own failed pass line, and noticed that EVEN-NONZERO's example vector probably fails HINGE-PROOF's domain hypothesis. Marked down because the bound is still about the reconstructed form of the first-order term, and the constant's inputs are certified by the lane's own ball arithmetic only.

**Checked by the lead, by hand: steps 2 and 4 to 10. They hold.**
- Step 2: sinh(rho) = lambda e^y (1 - e^(-2y)/(4 lambda^2)), so V_l = e^(-2y)(1 - e^(-2y)/(4 lambda^2))^(-2), decreasing in l; hence the bounds on A_l and G_l.
- Step 4: <b, u> = <b, W v> <= <b, W b>^(1/2) <v, W v>^(1/2) with <v, W v> = lambda Q_l.
- Steps 5 and 6: pi_s = cos(s h) a; cos x >= 1 - x^2/2; Jensen on the spectral measure of h u (total mass A, second moment G^2); sup of p^2 e^(-eps p^2) is 1/(e eps).
- Step 7: the integral of J0^2 (1 - s^2/s*^2)^2 over [0, s*] is (8/15) J0^2 s*, and s* = 1.777 < pi.
- Steps 9 and 10: Cauchy-Schwarz with the weight (2 sinh 2 pi h)^(1/2), and the two sup bounds; both sups are finite because 2 eps p^2 = 3 p^2 beats 2 pi p.
- Arithmetic of step 11 recomputed by hand from the six certified inputs: J0_lo = 0.36357, J1_hi = 0.23025, s* = 1.7771, M_lo = 0.007500, Den <= 3.7685, ||e|| <= 0.47154, c = 2.110e-3. Agrees.
- Not hand-checked: step 1 (domain), step 3's closed form for Q_l, and the six ball-arithmetic inputs themselves.

**Re-run by the lead:** `..\..\even-nonzero\venv\Scripts\python.exe even_passage.py <out.json>` from `scripts\`, 27 s, exit 1 with exactly one failed line, C2, as the lane reports. The table equals the lane's: 0.12179, 0.12395, 0.12450, 0.12464, 0.12467 at l = 16 to 4096 against the proved 2.110e-3; both plants caught; grid values of A and G lie inside the certified ranges.

**On the failed line C2.** The lane predicted the constant would land in [2.2e-3, 2.7e-3]; it came out 2.110e-3 because the certified version uses a slightly cruder majorant. A wrong prediction of a constant's size, not a violated inequality.

**Where the even half stands after this lane.** Lead-checked: uniform bound (EVEN-PROOF), limit operator (EVEN-LIMIT), non-vanishing in the limit (EVEN-NONZERO), and now a uniform positive lower bound at finite l on an explicit family. The measured value is about 58 times the proved floor, so the floor is safe but far from sharp. **Still not a theorem.** Missing, in order of weight: (1) a written definition of the first-order term on H + H', including that the first-order part of beta is (g/lambda) times it, to replace the reconstruction from code; (2) the remainder constant C(e_l)/||e_l|| bounded uniformly in l for this family, with a lower bound on ||e_l||; (3) EVEN-PROOF's dictionary making W_l the even charge, validated only through the spectrum; (4) HINGE-PROOF's Duhamel-Gronwall corollary and the Shale/Connes criterion, not lead-checked. The paper must still say "not proved" for the even divergence; it may say that the uniform bound and a uniform first-order lower bound are established.

---

## Lane report (unchanged)

**EVEN-PASSAGE report (2026-11-15)**

The liminf is positive with an explicit family and an explicit constant, by the shorter route: no strong resolvent convergence and no use of EVEN-LIMIT is needed. One pre-fixed line failed (C2, my prediction of the size of the constant).

## 1. Statement

Let u(y) = (1 − ((y−2)/4)²)⁴ on −2 ≤ y ≤ 6 and 0 elsewhere, placed at y = ρ − ln 2λ, so its support lies in ρ > 0 for every l ≥ 16. Let ε = 1.5, τ = π (t = 0.5), and

e_l = (q, π) = (0, a_l), a_l = e^{−ε h_l²} h_l² u.

Then for every l ≥ 16: λ‖β̂₁^(l) e_l‖ / ‖e_l‖ ≥ 2.110·10⁻³, and the liminf is at least 2.163·10⁻³. On the grid the actual value is 0.1218 rising to 0.1247.

The family also satisfies HINGE-PROOF's domain hypothesis: h_l⁻² a_l = e^{−ε h_l²} u, so z_l lies in the domain of coth^{3/2}(πh_l) and of e^{πh_l}.

## 2. Proof

Notation: A_l = ‖h_l u‖², G_l = ‖h_l² u‖, V_l = λ²/sinh²ρ = e^{−2y}(1 − e^{−2y}/4λ²)⁻², Q_l = ⟨u, W_l⁻¹u⟩/λ.

1. [P] u is C³ with compact support in ρ > 0, so u is in D(h_l²) with h_l²u = −u″ + V_l u (h_l² is essentially self-adjoint on C_c^∞ since λ² ≥ 3/4).
2. [P] V_l decreases pointwise to e^{−2y} as l grows. Hence A_∞ ≤ A_l ≤ A_16 and G_l ≤ G_∞ + ‖(V_16 − V_∞)u‖ =: G_hi.
3. [P] W_l is the Dirichlet Green operator of −d²/dr² + λ²/r² on (0,1), with dρ = dr/N. So u = W_l v with v = N(−u_rr + λ²u/r²), and ⟨v, W_l v⟩ = ∫(u_r² + λ²u²/r²)dr = λQ_l, where Q_l = ∫e^y(1 + e^{−y}/2λ)²(u_y² + V_l u²)dy ≤ Q_16.
4. [P, given W_l ≥ 0 bounded: I, EVEN-PROOF and EVEN-LIMIT lead-checked parts] Cauchy–Schwarz in the W-form: λ⟨b, W_l b⟩ ≥ ⟨b, u⟩²/Q_l for every b in L².
5. [P] With a_s = cos(s h)a: ⟨a_s, u⟩ = ⟨hu, cos(sh)e^{−εh²}hu⟩ ≥ J₀ − s²J₁/2, where J₀ = ⟨hu, e^{−εh²}hu⟩ and J₁ = ⟨hu, h²e^{−εh²}hu⟩ (cos x ≥ 1 − x²/2).
6. [P] Jensen on the spectral measure of hu: J₀ ≥ A_l·exp(−εG_l²/A_l) ≥ A_∞·exp(−εG_hi²/A_∞) =: J₀_lo. Also J₁ ≤ A_l/(eε) ≤ A_16/(eε) =: J₁_hi.
7. [P] The integrand of m_l is non-negative, so λm_l(τ) ≥ ∫₀^{s*}(J₀_lo − s²J₁_hi/2)²ds/Q_16 = (8/15)J₀_lo² s*/Q_16 =: M_lo, with s* = (2J₀_lo/J₁_hi)^{1/2} < τ.
8. [I, EVEN-NONZERO steps 1, 4, 6, lead-checked; re-derived here at finite l] ‖β̂₁e‖ ≥ ‖(2 sinh 2πh)^{−1/2}Tz‖ and Im⟨z, Tz⟩ = −m_l/2. T is bounded on z-space because ‖h^{1/2}W_l h^{1/2}‖ ≤ 1/λ (EVEN-PROOF step 5, lead-checked).
9. [P] Cauchy–Schwarz: ‖β̂₁e_l‖ ≥ m_l/(2·Den_l), Den_l² = ⟨a, h⁻¹sinh(2πh)a⟩ ≤ S′·A_16, with S′ = sup p·e^{−2εp²}sinh(2πp).
10. [P] ‖e_l‖² = ½⟨a, h⁻¹coth(πh)a⟩ ≤ ½C·A_16, with C = sup (p + 1/π)e^{−2εp²}.
11. [C] Ball arithmetic at 100 bits (radii below 1e-26): A_∞ = 0.938464, A_16 = 0.938807, G_hi = 0.770249, Q_16 = 16.70472, S′ ≤ 15.12692, C ≤ 0.473687. These give J₀_lo = 0.363568, J₁_hi = 0.230245, s* = 1.7771, M_lo = 0.0074997, Den ≤ 3.76846, ‖e_l‖ ≤ 0.47154, and c = M_lo/(2·Den·‖e‖) = 2.110·10⁻³.

The unbounded functions never need to converge: e^{πh} and coth(πh)/h appear only multiplied by p²e^{−2εp²}, as bounded functions whose sup is l-independent.

## 3. Check

| l | λ‖β̂₁e_l‖/‖e_l‖ [N] | uniform bound [C] | bound at this l [C] | λm_l | Den_l | ‖e_l‖ | λ·\|r_l/r_∞ − 1\| |
|---|---|---|---|---|---|---|---|
| 16 | 0.121788 | 2.110e-3 | 2.158e-3 | 0.029496 | 2.64584 | 0.26769 | 0.383 |
| 64 | 0.123952 | 2.110e-3 | 2.162e-3 | 0.030009 | 2.64478 | 0.26768 | 0.379 |
| 256 | 0.124501 | 2.110e-3 | 2.163e-3 | 0.030140 | 2.64471 | 0.26768 | 0.378 |
| 1024 | 0.124639 | 2.110e-3 | 2.163e-3 | 0.030173 | 2.64470 | 0.26768 | 0.378 |
| 4096 | 0.124673 | 2.110e-3 | 2.163e-3 | 0.030182 | 2.64470 | 0.26768 | 0.378 |
| limit | 0.124685 | | 2.164e-3 | 0.030184 | 2.64470 | 0.26768 | |

- **Intermediate inequalities:** every one held at every l (P1–P6). The pointwise integrand is at least 3.1 times its proved floor.
- **Grid stability:** halving dy moves r by 2.4e-5; doubling the box by 3e-14.
- **Plant A** (factor λ dropped): values fall from 7.4e-3 to 3.0e-5, below the bound from l = 64 on. Caught.
- **Plant B** (lapse weight dropped in Q, i.e. W⁻¹ confused with h²): the planted M_lo = 0.1334 is violated at every l. Caught.
- **FAILED as pre-fixed, C2:** I required 2.2e-3 ≤ c ≤ 2.7e-3 and got 2.110e-3. My scan used the true sup of p·coth(πp)e^{−2εp²}; the script certifies the majorant (p + 1/π)e^{−2εp²}, which costs 18%. The bound is valid; my size prediction was wrong. Nothing was re-run.

## 4. Re-run

From `record\wave73\even-passage\scripts`:

`..\..\even-nonzero\venv\Scripts\python.exe even_passage.py ..\out\even_passage.json` (28 s, exits 1 because of C2)

Logs are in `logs\`: `01_prefixed_pass_lines.log`, `02_even_passage_run1.log`, `03_note_on_run1.log`. `scripts\scan_constants.py` is the exploration that chose (b, y₀, ε); it builds no operator. Nothing committed, no REPORT.md.

## 5. Weakest step

Step 8: the bound is about β̂₁ in EVEN-NONZERO's reconstructed form (cP, sQ̄), which is still not a written definition.

## 6. Still missing before "the even Shale sum diverges" is a theorem

- **Remainder constant for this family:** sup_l C(e_l)/‖e_l‖ < ∞ is not checked. It looks reachable, since ⟨z, coth³(πh)z⟩ ≤ ½·sup[p³coth³(πp)e^{−2εp²}]·‖u‖² [U]. It also needs a lower bound on ‖e_l‖ (0.2677 on the grid, [N] only).
- **EVEN-NONZERO's own instance** e = (0, e^{−h²/4}π₀) appears not to lie in D(coth^{3/2}), if its quoted spectral density is right [U]. Use this lane's family in the paper.
- **A written definition of β̂₁ on H + H′,** including that the first-order part of β is (g/λ)·β̂₁.
- **EVEN-PROOF step 1** (the dictionary making W_l the even charge): validated only through the spectrum.
- **HINGE-PROOF items:** the Duhamel–Gronwall corollary and the Shale/Connes criterion are not lead-checked; the even tower being block diagonal with a well-defined charge (step 9) and the physics gaps (steps 1, 12) remain.
- The ball-arithmetic step [C] has not been independently re-run by the lead.
