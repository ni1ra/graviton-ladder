# REPORT: lane EVEN-TOWER (wave 73)

## Lead's check and grade (2026-11-15)

- Re-run by the lead: `python even_tower.py 3 3` exits clean with no FAIL and no missed plant; the momentum block equals c_l[(1/2) int N psidot^2 + (1/2) int psidot G_D psidot] on the whole Gram with difference 0, the spectrum is 10/9, 23/22, 40/39, and the gauge-moved random datum returns the same charge. Lane's own code re-executed; `shale_check.py`, `even_canonical.py` and the HINGE-E re-runs were not re-run by the lead.
- This lane corrects the lead twice. (1) The brief called the even tower untouched; GEN-LAW and `hinge-e/` had already worked on it. (2) The lead had told the author that stretching gravitons balance and that the mismatch is twist-only. That is true only on stationary data. On general data both parities balance on metric data and fail on extrinsic-curvature data. The ladder's contrast is metric against extrinsic curvature, not electric against magnetic.
- Standing: the area-side charge is now defined and gauge invariant on every graviton sector of the ball, and differs from the vacuum's generator on every sector. The odd-tower non-existence result is checked to the standard of the earlier lanes; the even-tower divergence (constant about 0.0136) rests on HINGE-E's unreported scripts plus this lane's re-run on the same grid, with no proof of the remainder and no blind measurement yet.
- Kept in view: the lane's own Dirichlet-wall control failed as pre-registered and was diagnosed (box size, Neumann wall) with lines fixed in advance; the pull-back argument for "moving the region is a gauge move" is by hand.
- Grade: 378 of 420. Computed the geometric charge instead of taking it on trust, gauge tests on non-TT data, a failed control reported and diagnosed properly, numbered logs; marked down because the even Shale constant is inherited and general l leans on the record.

The lane's report follows unchanged.

**EVEN-TOWER report (2026-11-15)**

Tags: [P] proved this session by hand, [N] exact or numeric this session, [S1] source read first-hand, [S2] second-hand, [I] inherited unchecked, [B] belief.

## 1. Verdict

The even tower neither behaves like the odd one nor balances: it balances on metric data and fails on extrinsic-curvature data, with a non-local defect.

- **Q1: yes.** A gauge-invariant second-order area-side charge exists on even data. It is the stationary value of Area − 2·Vol over the edge's normal shift, the corner's time shift and the slice, and it is unchanged under explicit even gauge moves [N, l = 2 to 6]. No residual freedom is left for l ≥ 2. The canonical charge read in TT gauge is gauge dependent on even data [N].
- **Q2: not one weight.** A_ext = c_l·E_1[ψ_E, ψ̇_E] + Z_E, with Z_E = ½c_l⟨ψ̇_E, G_D ψ̇_E⟩ and G_D the Dirichlet Green operator of the ball [N, every Gram entry]. The field block is exactly weight 1. The record's "spectrum not one weight" is reproduced: 1 + 1/((j+1)(2l+2j+3)).
- **Q3: yes.** The vacuum side is 2π·c_l·E_1[ψ_E] by the subspace lane's method on even data [N for the identities; Hilbert-space steps I].
- **Q4: a different, non-geometric flow.** It is well defined per sector but is not the vacuum's. Its Shale sum diverges logarithmically with per-sector constant a_E ≈ 0.0136 against the odd 7ζ(3)/(16π²) ≈ 0.0533. The constant is HINGE-E's [I]; I re-ran its scripts and reproduced its exact sector numbers with my own code at l = 16 to 128 [N]. Its remainder is not proved to HINGE-PROOF's standard.

**Theorem scope, one sentence:** on every graviton sector of the flat ball, both parities, the gauge-invariant fixed-volume area-side charge equals the vacuum's modular generator on metric data and exceeds it on extrinsic-curvature data by a positive form (local ½c_l∫ψ_B² for odd parity, non-local ½c_l⟨ψ̇_E, G_Dψ̇_E⟩ for even), so its flow is the vacuum's on no sector and each tower's Shale sum diverges logarithmically.

Certainty: high for the charge identities at l = 2 to 6 on polynomial profiles. General l rests on the record's symbolic proofs [I]. Medium for the even Shale constant as a theorem.

One correction to the brief: this half was not untouched. GEN-LAW held the closed form for the even momentum half, and `hinge-e/` holds a Shale test of it (PREREG, scripts and outputs, no report on file). Both took the geometric charge on trust; this lane computes it.

## 2. Definitions as used

Conventions are charge-id's: d = 4, unit ball, t = 0, N = (1−r²)/2, 1/8G dropped, integrals in units of 4π.

- **Even TT tensors.** T[q] = f₁Hessφ + f₂P_φ + f₃φδ + f₄φX, with x·T·x = qφ (`wave36/companionA8/NOTE.md` line 163, Lemma 7). The script checks div, trace and x·T·x exactly. q = 1 is the gauge mode Hessφ/(l(l−1)).
- **Even gauge moves** [P]. h → h + 2 sym∇ξ, K → K + Hess X⁰, with ξ = a(s)∇φ + b(s)φx and X⁰ = c(s)φ. σ = ξ·ν at r = 1 moves the edge, τ = X⁰ at r = 1 moves the corner in time, and X⁰ inside moves the slice.
- **Geometric charge.** A_un(h, K) = ∫N(K·K − (tr K)²) − ∫N·R⁽²⁾[h] − 2δ²F[h], F = Area − 2·Vol (NOTE.md line 133, Proposition 4; charge-id's Q2). No TT simplification is used.
- **Extremisation** [P]. By diffeomorphism invariance, F of the region displaced by εX equals F of the coordinate ball in the pulled-back metric. Its first-order data are the gauge-moved (h, K), and the second-order metric again drops out through the N-smeared constraint. So A_ext(a) is the stationary value of A_un over the gauge orbit of a.
- **Fixed true volume** [P]. With no l = 0 part at first order, Area at fixed volume and Area − 2·Vol have the same second-order coefficient.
- **Chain variables.** ψ_E = x·Ric[h]·x, ψ̇_E = x·Ric[2K]·x (NOTE.md lines 161 and 201). E_w as in the brief.

## 3. Computations

Script `even_tower.py L 3`: 55 of 55 checks and 10 of 10 plants at each of l = 2, 3, 4, 5, and at l = 6 with KMAX 2, on the first run (3 to 18 s).

**How the ball moves and what the extremal surface is.**
- A_un is unchanged by every ξ tangent to the edge or supported inside; it depends on ξ only through σ [N].
- On pure gauge data the formula returns flat-space geometry: −(l−1)(l+2)∮σ² for the edge, ∮|∇_tan τ|² − 2∫|∇X⁰|² for the time shift [N].
- For TT K: A_un(K + Hess X⁰) = ∫N|K|² + 2∮κτ + ∮|∇_tan τ|² − 2∫|∇X⁰|², κ = x·K·x [P, then N].
- Stationary point: X⁰ harmonic (maximal slice, tr K = 0), τ = −κ/(l(l−1)) (the moved K has no normal-normal part on the edge), and σ set by companion A8's Lemma 8 [N].
- It is a saddle: minimum in σ and the slice, maximum in τ.

**The charge.**

| Block | Result (exact, whole Gram) |
|---|---|
| Field (h) | A_ext = c_l·[½∫N\|∇ψ_E\|² + ½∫ψ_E²], weight 1 |
| Momentum (K) | A_ext = ∫N\|K\|² − ∮κ²/(l(l−1)) = c_l·[½∫Nψ̇_E² + ½∫ψ̇_E·G_Dψ̇_E] |
| h–K | zero by construction (σ couples to h, X⁰ to K) |

- At l = 2 on T[s] as momentum: A_ext = 32/405, Z_E = 4/405, ratio to the vacuum part 1/7.
- Spectrum of the momentum block against the vacuum's: 8/7, 19/18, 34/33 at l = 2; 10/9, 23/22, 40/39 at l = 3; 12/11, 27/26, 46/45 at l = 4.
- Flow on the curvature: x·Ric[2N·K̃]·x = Nψ̇_E + G_Dψ̇_E, with K̃ the stationary representative. It is non-local in ψ̇_E.
- Both towers share one ladder: 1 + 2/(n(2l+n+1)), n odd for odd parity and n even for even parity [N, first three rungs]. The headline (l+2)/(l+1) is rung n = 1.

**Gauge tests.**
- Random datum with h and K independent, moved by random ξ (a, b of degree 2) and X⁰ (degree 2), pipeline run on the non-TT data: A_ext unchanged, for example 35.595804 and 3.5974376 at l = 2.
- The unextremised values move: 33.93 → −30.97 and 7.39 → 71.27.
- The stationary point holds against directions outside the solve family. A two-harmonic datum also passes.

**Canonical charge** (`even_canonical.py`, exploratory, l = 2, 3, 4; makes charge-id's Step 9 concrete).
- ½Ω_B(a, δʷa) on even TT data is symmetric and equals the record's closed form (4).
- It is not invariant under the residual TT gauge mode in either slot, and is negative on a field datum (−1/5 against A_ext = 7/15 at l = 2).
- It differs from A_ext by +∮κ²/(l(l−1)) on the K slot and −¼∮(ν·h·ν)² − ∮G²/((l−1)(l+2)) on the h slot.

**Vacuum side.**
- On ℝ³ with Gaussian profiles: ∫K_a·h_b = c_l∫ψ̇_a ψ_b, one constant, κ = +1.
- Two-point function, both blocks, whole Gram: ∫½ω h·h = c_l∫ω ψψ and ∫2K·K/ω = c_l∫ψ̇ψ̇/ω (8/3 and 16/15 at l = 2).
- Paper A's flow on (E, B) gives x·δE·x = Nψ̇_E and div(N∇ψ_E) − 1·ψ_E: weight 1.
- Local both ways: E = T[q_E] with q_Eφ = ψ_E, and B is the odd-type tensor with profile ψ̇_E/(l(l+1)). ψ_E obeys the wave equation.

**Controls and plants.**
- Odd data through the same code, at degree l and at l+1 (same polynomial parity as the gauge family): zero drive, A = c_l·E_2 on the Gram, and (l+1)/(l+2) on the stationary datum.
- Stationary electric datum h = −2φδ (not TT) through the full pipeline: A_ext/N_Y = 1/7, 4/15, 4/11, 40/91, 1/2 at l = 2 to 6 = 2l(l−1)/((l+2)(2l+3)), ratio 1. Unextremised it is −6/7 at l = 2.
- Plants caught: corner time shift frozen; slice term dropped; wrong lapse; weight 2 on the field block; Neumann inverse; coefficient 2 on Z_E; three on the vacuum side.

**Shale.**
- HINGE-E's scripts re-run here: e1 10/10, e2 20/20, e3 15/16. The failing line is its g = 0.25 test at l = 512, where 1 − A is 3·10⁻⁹ against a g = 0 noise of 10⁻¹⁰; the g = 1 numbers agree to 10⁻⁵ up to l = 256.
- Own code (`shale_check.py`, covariance-matrix affinity): its odd-tower control failed as pre-registered with a Dirichlet wall (1.535·10⁻⁴ against 1.692·10⁻⁴), so nothing was claimed from that run.
- Diagnosis tested with lines fixed in advance: doubling the box halved the deficit (9.3% to 4.5%), and a Neumann wall reproduces the control and HINGE-E's even numbers to four digits.
- l²(1 − A_l) = 0.0102, 0.0115, 0.0124, 0.0129 at l = 16 to 128, rising towards 0.0136.

## 4. The flagship paragraph, amended

**May say:** "On every graviton sector of a flat ball, of either parity, the second-order area deficit at fixed volume, which is what Jensen–Sorce–Speranza's Section 3.3 proposal gives once the edge is extremal and the slice maximal, is gauge invariant and equals c_l·E_1 of the sector's scalar chain field on metric data, and exceeds it on extrinsic-curvature data by a positive form Z. For odd parity Z = ½c_l∫ψ_B², so the charge is c_l·E_2 and is also the unique gauge-invariant canonical charge of the conformal Killing flow. For even parity Z = ½c_l⟨ψ̇_E, G_Dψ̇_E⟩ with G_D the ball's Dirichlet Green operator, which is non-local and no weight. The conformal Killing diffeomorphism acting on the Weyl curvature is the weight-1 flow in both parities, generated by c_l·E_1, the vacuum's modular generator. No gauge or region choice in their text removes Z. Hence the area-side flow is not the modular flow of any normal state or weight in the vacuum sector: the Shale sum diverges logarithmically in each tower, with per-sector constants 7ζ(3)/(16π²) for odd parity and about 0.0136 for even parity." Then continue with ss-second's Section 5 sentence unchanged.

If asked what is assumed: the ball algebra is the curvature net; exact checks are at l = 2 to 6; the even constant's remainder is not proved.

**Must not say:**
- "A4 fails for gravitons", "no equilibrium state exists", "this is the charge JSS define" (all unchanged).
- "The electric sectors balance" or "the mismatch is a twist-only phenomenon": that holds only on data with K = 0.
- "The even charge is c_l·E_w" for any w, or "the even tower has weight 1".
- "The canonical charge" for even parity: there it is gauge dependent, and only the extremised one is defined.
- "Odd-parity sectors only" is no longer needed; "proved" for the even Shale divergence is not yet earned.

## 5. Stretch versus twist, and the known weakness

- The ladder's contrast was never parity. The stationary electric datum is pure h and the stationary magnetic datum is pure K. On general data both parities balance on h and fail on K.
- The defect is a functional of the slice's magnetic curvature in both towers: ψ_B = x·B·x for odd, and ψ̇_E, which is algebraic in the even sector's B, for even. GEN-LAW's item 4 said this [I]; it is now checked from the geometric definition.
- The weakness stands and is sharper. The vacuum's generator is duality symmetric (E_1 in both towers). The area rule is not: it sees K as bare ∫N·K̃·K̃ on a maximal slice. So the asymmetry belongs to the fixed-volume maximal-slice rule, now as "metric against extrinsic curvature" rather than "electric against magnetic". A duality-blind rule such as the cone cannot reproduce it.

## 6. Unchecked, ranked, and the next experiment

1. The even Shale divergence as a theorem: no analogue of HINGE-PROOF's Lemmas A and B for W. HINGE-E has no report or lead grade, and its flat-lattice run (e4) did not close under its pre-registered convention. My check covers the same grid only.
2. That JSS's Section 3.3 means extremising over the corner's time shift as well. For even data this is decisive: without it the value is gauge dependent [I from jss-id's reading; JSS not re-read here].
3. General l: l = 2 to 6 here. The record's symbolic proofs (companion A8 Theorem 3; GEN-LAW's radial proof) were not re-run by me.
4. The pull-back argument for "moving the region is a gauge move" is by hand; the machine checks are its consequences.
5. Hilbert-space steps of the vacuum side (onto, second quantisation, curvature net) [I].
6. Cross terms between parities and between different l: by symmetry, not computed. d ≥ 5 untouched.
7. Sources: only the record was read first-hand. Benedetti–Casini (arXiv:1908.01800) was seen as a tool summary of the abstract [S2]: two modes per l with scalar dynamics, parities not stated there.

**Next experiment:** a blind flat-frame lattice measurement of the even tower's per-sector constant, as HINGE-BLIND did for the odd one.

## Reproduce

From `record\wave73\even-tower\scripts`, native Windows Python 3.14:
- `python even_tower.py 2 3`, then `3 3`, `4 3`, `5 3`, `6 2`
- `python even_canonical.py 2 2`, then `3 2`, `4 2`
- `python shale_check.py neumann` (about 30 s); `dirichlet` and `dirichlet2` for the failed control and its diagnosis
- HINGE-E re-run, from `..\..\hinge-e\scripts`: `python e1_structure.py ..\..\even-tower\out\rerun_e1_structure.json`, likewise `e2_kernel`, `e3_hyperbolic`

Files added, all in `record\wave73\even-tower\`:
- `PREREG.md` (with Addenda A, B, C, each written before its script ran)
- `scripts\even_tower.py`, `scripts\even_canonical.py`, `scripts\shale_check.py`
- `logs\run01_l2.log` to `run01_l6.log`; `run02_hinge_e_e1.log`, `_e2`, `_e3`; `run03_shale_check.log` (failed control); `run04_shale_dirichlet2.log`; `run05_shale_neumann.log` (all lines printed, crashed writing JSON); `run06_shale_neumann.log`; `run07_canonical_l2.log`, `_l3` (script error); `run08_canonical_l2.log`, `_l3`, `_l4`
- `out\even_tower_l2.json` to `_l6.json`; `rerun_e1_structure.json`, `rerun_e2_kernel.json`, `rerun_e3_hyperbolic.json`; `shale_check.json`, `shale_check_dirichlet2.json`, `shale_check_neumann.json`

Nothing committed, no report file written, no subagents, nothing outside this folder modified.
