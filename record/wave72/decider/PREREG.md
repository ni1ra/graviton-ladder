# w72-decider: pre-registration

Tuesday 10 November 2026, Oslo. Main thread (Fable), no lanes. Written and pushed before any script in this folder ran.

## Question

The ladder theorem (wave69/final/THEOREM.md) says the magnetic sector misses the second-order first law of the unit ball by A/((l+2)·8G). Whether that is new physics or a known effect renamed depends on the electric sector, which the theorem leaves out:

- If the electric area side misses by the opposite amount, the total miss is proportional to B² − E², a Lorentz scalar, and Jacobson's local scalar term can absorb it (the Casini–Galante–Myers class of effect).
- If the electric side balances, or misses with the same sign, the miss is not a scalar and cannot be absorbed that way.

## Method (new to the record)

The exact static axisymmetric vacuum solution (Weyl class):

  ds² = −e^{2U}dt² + e^{−2U}[e^{2k}(dρ² + dz²) + ρ²dφ²],

with U = ε r^l P_l(cos θ) harmonic and k fixed by quadrature (k_ρ = ρ(U_ρ² − U_z²), k_z = 2ρU_ρU_z). Its linear part is the static electric multipole E_ij = ∂_i∂_jU, B = 0, the duality partner of the scope family. The slice t = 0 is time symmetric, so stationary surfaces of Area − 2·Volume lie in it. The script expands F = Area − 2·Volume on r = 1 + ε a P_l to O(ε²) in exact rational arithmetic and takes the stationary value in a. A stationary value of an invariant functional does not depend on coordinates, so no gauge choice enters.

Modular side: electric-magnetic duality is a real rotation of (E, B). It preserves each local subspace H(U) and commutes with the conformal group, so by Theorem F it commutes with the ball modular group, and S(electric datum) = S(dual magnetic datum) = (c_l/8G)·E₁[ψ], with ψ = x·E·x = l(l−1)U.

Reported number: X_l = [−(1/4G)·Δ₂F] ÷ [(c_l/8G)·½∫_B ψ²], the area side in units where the modular side is l + 1.

## Bins (exact arithmetic, so equality is the test), at l = 2, 3, 4, 5, 6

| Bin | X_l | Area side in units 1/8G, U = r^l P_l | Meaning |
|---|---|---|---|
| MOD | l + 1 | 4π/35, 16π/105, 16π/99 at l = 2, 3, 4 | electric side balances; the miss is B-only, not a scalar |
| CAN | l + 2 | | area side is duality symmetric; the miss is of Bel–Robinson type, not a scalar |
| SCAL | l | | total miss ∝ B² − E², absorbable; the ladder is a known class of effect |
| OTHER | anything else | | reported as found |

The record expects MOD (CALIBRATION: w39-BR, w49, w57 give 4π/35 and 16π/105 by other routes; the w57 comparison was not pre-registered). My own lens from this morning (area weights 0, 1, 2) put real weight on SCAL. My odds before running: MOD 0.6, SCAL 0.2, CAN 0.1, OTHER 0.1.

## Controls that must pass for the result to count

- C1. First-law identity δ(Area − 2Vol) = −½∫_B N δR on random polynomial metric perturbations of the flat ball.
- C2. Normalisation: x·(curl K)·x = ±½(l−1)l(l+1)φ and ∫_B N|K|² = c_l E₂[ψ] on the scope family, so that E and B are normalised alike.
- C3. The Weyl 3-metric has R = 0 through O(ε²).
- C4. l = 1 (U = εz, flat Rindler space): Δ₂F must be exactly 0.
- C5. Plants: drop k; drop the surface shift; change e^{−2U} to e^{−U}. Each must change X_l.

## Forced items (no evidential weight, stated for the paper)

- Scalar control: by C1, a minimally coupled scalar with a static harmonic datum has area side 2π·E₀, while Longo–Morsella gives S = 2π·E₁. Miss +π∫φ², a local scalar density of positive sign.
- Lattice rows of wave55 against the rival ratio (l+2)/(l+1): old data read again, no new run.

## What I will call it (answer to The author, 10 Nov: "if this goes thru, will u agree on eureka?")

Yes, if all of these hold:

1. C1–C5 pass.
2. X_l lands in MOD or CAN at every l computed.
3. The scalar control shows the opposite sign and a scalar density, so the graviton miss is not the scalar effect under a new name.

Then I agree on eureka for this statement: in free linearised gravity on a flat ball, the second-order first law holds on static electric multipoles and misses on magnetic ones by exactly A/((l+2)·8G), and the miss is not a local scalar. One hypothesis stays open either way and is said in the paper: H0, that the free net is the G → 0 limit of the gravitational ball algebra.

If X_l lands in SCAL: no eureka. The theorem still stands, and the paper is rewritten as "the graviton case of a known correction".

Nothing here books credence; the rubric owns that.
