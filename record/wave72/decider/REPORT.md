# REPORT: w72-decider. The electric side balances exactly; the ladder's miss is magnetic only and not a scalar

Main thread (Fable), Tuesday 10 November 2026, Oslo. Claude only. Pre-registered in `PREREG.md` (commit 90ed148ae, pushed before the script ran). Nothing is booked here; the rubric owns credence.

## Result

Bin **MOD** at every l computed. On the exact static Weyl solution with U = ε r^l P_l, the stationary value of Area − 2·Volume drops at O(ε²) by exactly the amount the modular side predicts:

| l | Area side, units ε²/8G | Surface shift a* | X_l | Modular side l + 1 |
|---|---|---|---|---|
| 2 | 4π/35 | 1/2 | 3 | 3 |
| 3 | 16π/105 | 2/5 | 4 | 4 |
| 4 | 16π/99 | 1/3 | 5 | 5 |
| 5 | 160π/1001 | 2/7 | 6 | 6 |
| 6 | 2π/13 | 1/4 | 7 | 7 |

- The l = 2, 3 values equal the record's 4π/35 and 16π/105 (w39-BR, w57), reached there by other routes. l = 4 was registered as 16π/99 before the run. l = 5, 6 are new.
- The stationary surface is r = 1 + ε·(2/(l+2))·P_l.
- All controls passed: C1 (first-law identity, three random data), C2 (normalisation of E against B, l = 2, 3, 4, two harmonics each), C3 (R = 0 through ε², l = 2, 3), C4 (l = 1 Rindler gives exactly 0), C5 (three plants at l = 2, 3 each move X_l: to −9, −18, 39/4 and −6, −20, 17/2).
- Output: `out/d1_weyl.txt`. Rerun: `python3 scripts/d1_weyl.py` (14 s, exact rational arithmetic, exit code 0).

## What it means

In units where the entropy of a degree-l multipole is l + 1:

```
                 entropy     area side     first law at O(eps^2)
electric (E)      l + 1        l + 1        holds
magnetic (B)      l + 1        l + 2        misses by 1 part in l + 2
```

- The entropy cannot tell E from B: duality is a symmetry of the free net, its local subspaces and the ball's modular group. The area side can. So the two cannot agree on both sectors, and they agree on E.
- The miss is ∝ B² alone at l = 2. A local Lorentz scalar would have to be ∝ B² − E². So Jacobson's scalar term cannot absorb it, and it is not the Casini–Galante–Myers correction under another name.
- SCAL, the outcome I feared this morning, is excluded exactly. My "area weights 0, 1, 2" picture was wrong for the electric row: its weight is 1.

## Forced items (stated for the paper, no evidential weight)

- **Scalar control.** By C1, a minimally coupled scalar with a static harmonic datum φ has area side 2π·½∫N|∇φ|², and Longo–Morsella gives S = 2π(½∫N|∇φ|² + ½∫φ²). The mismatch is +π∫φ²: entropy above area side, a local scalar density. The graviton miss has the other sign (entropy below area side) and is not a scalar density. A conformally coupled scalar balances once the Wald term −(π/3)∮φ² is counted.
- **Lattice, read again.** wave55's 3D graviton lattice has one purely magnetic datum, g2 (uniform B, l = 2): S_lat/S_T1 → 0.999, bracket [0.891, 1.021]. If the entropy equalled the area side the ratio would be 4/3. The other single-l rows (g1, g3, g5) are electric, where area side and entropy agree, so they test the entropy formula (1.019, 1.008, 0.991) but not the miss.

## Correction to my own earlier statement

On 10 Nov I told the author the lattice had measured "the l = 2, 3, 4 rungs" against rivals 1.33, 1.25, 1.20. Only l = 2 is a magnetic row. The l = 3 and l = 4 rows are electric and carry no rival.

## Limits

- The modular side on electric data is taken from duality plus Theorem F (wave70), written in PREREG; it agrees with the record's direct electric chain value. It was not re-derived here.
- Axisymmetric data only. Other m at fixed l follow because the second-order stationary value is an SO(3)-invariant quadratic form on an irreducible representation, so it is a multiple of ∮U². That argument is one line and was not scripted.
- Static electric data only, the duality partner of the theorem's stationary magnetic data.
- No blind check. The controls are in the script.
- Not done: the comparison with Jacobson–Senovilla–Speranza (1710.07379) eq. (27), a different prescription (fixed-radius geodesic ball) that is also not duality symmetric. Their formula was read from ar5iv only (wave18).
- H0 and H4 are untouched. This run settles what the miss is not; it does not settle whether the free net is the G → 0 limit, or whether extremal Area − 2·Volume is the physical area side.

## Second test: the Killing limit (PREREG2.md, `scripts/d2_killing.py`, `out/d2_killing.txt`)

the author asked for H4 and H0 to be tested. H4 can be attacked where the answer is known: on a Killing horizon the second-order first law must balance, so the same machinery must give no miss on the Rindler wedge.

**Result: bin KILLING-OK, with one registered control failing as written.**

- For any lapse N and any TT harmonic K: Ric₃(2NK) − 2 sym(∇N × curl K) = 2 sym(K·∇∇N) − (ΔN)·K. Verified on 24 random data with a random quadratic lapse.
- Wedge, N = x₁: geometric rate = modular rate on every datum. No miss.
- Ball: the difference is exactly K on every datum. So the ladder's "+1" between l + 1 and l + 2 is 2 sym(K·∇∇N) − (ΔN)K with ∇∇N = −δ: the failure of ζ to be Killing, nothing else.
- **The plant control failed as registered.** I wrote "the plant must break the wedge equality on every datum". It was uncaught on 16 of 24. Post hoc diagnostic (`scripts/d2_diag.py`, `out/d2_diag.txt`, written after the failure): those 16 are the Hessian and double-curl families, which have curl K = 0 or curl curl K = 0, so both rates vanish on the wedge and no coefficient can show. The plant is caught on all 8 data with a non-zero modular rate (the scope family, degrees 2 to 5). The control was badly designed, by me: two of my three "families" are pure gauge as momenta. The test has power on 8 data, not 24.
- **Order of events, stated exactly.** PREREG2 was committed before the run (local commit 43e8d37e3), but its push was rejected because the remote had moved. It reached the remote after the run, rebased, as 3f79a6be4. The file's own line "pushed before the script ran" is therefore false; the commit time is before the run, the push time is not.
- Weight: small. I had derived the identity by hand before running, and said so. Its value is structural: it shows the machinery balances where balance is known, and it locates the miss.

**Also (algebra, no script).** At this order the stationary value of Area − 2·Volume equals the area deficit at fixed volume of the isoperimetric surface (Jacobson 2015): fixing the volume sets the second-order radial shift b = −V₂(a)/4π, and then A₂ + 8πb = F₂(a) for every a.

**H0 was not tested.** No script can. It needs a construction of the gravitational ball algebra at second order in which the vacuum modular generator is compared with the constraint (gate G1c). The identity above sharpens what H0 has to face: on the ball the diffeomorphism flow of ζ on the metric and the vacuum modular flow move the curvature at different rates, so the vacuum is not a state whose modular flow is that geometric flow.

## Eureka terms (PREREG)

All three registered conditions hold: controls pass, X_l is in MOD at every l, and the scalar control has the opposite sign and is a scalar density. So I agree on eureka for the registered statement: in free linearised gravity on a flat ball, the second-order first law holds on static electric multipoles and misses on magnetic ones by exactly A/((l+2)·8G), and the miss is not a local scalar. H0 stays open and is said in the paper.

## Recommended next

1. A blind read of the short paper (`wave72/paper/`), per GOAL.md, before it goes to an outside reader. Not run: the author asked for main thread only.
2. The JSS 2017 comparison above: one script on the same Weyl solution with their prescription.
3. Mixed E + B at fixed l: the cross term must vanish by parity; one line, unscripted.
