# REPORT: lane RULER-WHY (wave 73)

## Lead's check and grade

Sunday 15 November 2026. Small exact-computation lane (Fable): which ingredient of the fixed-volume rule (maximal slice, volume constraint, extremal edge) makes the metric-versus-extrinsic-curvature asymmetry, on the two static tides, and is RULER's "weight v^3 f" coincidence an identity? The lane returned the text below; the lead saves it.

**Grade: 366 of 420.** A clear answer with exact numbers at l = 2, 3, 4, the inherited values reproduced by an independent brute-force route (the exact Weyl metric, not the earlier lane's constraint identity), a pre-registration that discloses what the lane had already seen, a plant caught, and a correct "must not say". Marked down for three stated limits: the "size switched alone" row rests on one chosen definition of size by radius, with a hand-derived term that only one printed number tests; the maximal-slice ingredient was handled by an argument, not computed; and everything is on the two static tides, not on general data.

**Checked by the lead.**
- The pre-registration predates the run logs and discloses that the lane had opened RULER's plant output first.
- The magnetic column's logic holds by parity: the stationary magnetic datum has no first-order scalar on the slice, so no rotation-invariant rule can drive an edge shape or time shift with it; only the size rule can move its value.
- Consistency: case (a) reproduces 1 and (l+2)/(l+1); case (b) reproduces the earlier 2/9, 5/32, 3/25; the stationary shape a = 2/(l+2) is recovered (1/2, 2/5, 1/3).
- Not checked by the lead: the geodesic-bending and deflection terms behind rows (c1) and (c2).

**Re-run by the lead:** `C:\Python314\python.exe why.py <out.json>` from `scripts\`, 18 s, exit 0, 46 checks, 0 failures; the table at l = 3 and 4 and the v^3 f rows at l = 2 to 5 equal the lane's. The plant run was not repeated by the lead.

**What it changes.** On the static tides the asymmetry is a combination, not one ingredient: "equals the entropy on metric data" needs both the volume constraint and the extremal edge (either switch alone destroys it), while the excess on extrinsic-curvature data is set by the size rule alone, because no edge can couple to that datum. The v^3 f coincidence is an identity on the two lowest profiles only and matches the odd tower only; it is not an identity of forms and must not be quoted as one.

---

## Lane report (unchanged)

**RULER-WHY**

**(1) Answer.** The asymmetry is a combination, not one ingredient: "equals S on metric data" dies under either single switch (edge not extremised, or size not volume), while the excess on extrinsic-curvature data moves only when the size rule changes. The maximal-slice ingredient does nothing at first order in tr K, by argument only.

**(2) Table: area side / S (electric = metric, magnetic = K)**

| Setting | l=2 | l=3 | l=4 |
|---|---|---|---|
| (a) all three on | 1, 4/3 [I] | 1, 5/4 [I] | 1, 6/5 [I] |
| (b) volume, geodesic ball | 2/9, 4/3 [I] | 5/32, 5/4 [I] | 3/25, 6/5 [I] |
| (b′) volume, diamond edge | 2/9, 4/3 [I] | 5/8, 5/4 [I] | 39/50, 6/5 [I] |
| (c1) geodesic ball, fixed geodesic radius | 14/45, 28/15 [I] | 45/224, 45/28 [C] | 11/75, 22/15 [C] |
| (c2) fixed mean geodesic radius, edge extremised | 28/15, 28/15 [C] | 15/8, 45/28 [C] | 53/30, 22/15 [C] |
| (d) constant-mean-curvature slice, first order | no change [U] | same | same |
| Proper time between tips, bare edge (Wang) | −56/9, 0 [I*] | – | – |

Every [I] entry was reproduced by the script; the (c1) l=2 pair is Wang's printed pair.

**(3) Switches.**
- **(a)** Area − 2·Volume on the t = 0 slice, stationary in the edge shape a (found at a = 2/(l+2)).
- **(b)** Same functional on a fixed shape: geodesic ball a = 1/(l+1), diamond edge a = 2/(l+1).
- **(c1)** Geodesic ball, bare area at fixed geodesic radius; this switches size and edge together.
- **(c2)** Size is the geodesic radius averaged over directions at the centre, with the edge extremised under it; this switches size alone.
- **(d)** Not computed: the maximal slice maximises volume, so the first-order change vanishes; second order was not attempted.
- **Magnetic column:** the datum has no scalar on the slice (checked), so no rule can drive an edge shape or time shift, and only the size rule can move it.

**(4) v³f verdict.** It is an identity on the two lowest profiles, not of forms.
- Static tides give (l+2)/(l+1) at l = 2 to 5, and for symbolic n by the profile law. [C]
- The ψ̇ = φ rung gives exactly 1 at l = 2 to 5, in both parities. [C; unexplained]
- It matches the odd tower only: the even K rung's fixed-volume value is 1 + 1/(2l+3), not 1. [I, not recomputed]
- At l = 2 the full form agrees with the fixed-volume odd form in only two entries. [C on RULER's plant matrices, I]

**(5) Plant and gates.** Replacing "− 2 Volume" by "− 1 Volume" moves (a) to 101/6 and 2/3 at l = 2; caught at l = 2, 3, 4. The final run has 46 checks and 0 failures; no gate failed. The first run covered l = 2, 3 only; ranges were widened afterwards with no formula change. Before the pre-registration I had seen RULER's plant matrices, which is disclosed there.

**(6) Re-run.** From `record\wave73\ruler-why\scripts`:
`C:\Python314\python.exe why.py ..\out\why.json` (19 s); add `--plant k1` for the plant. This uses the global sympy 1.14, not a venv.

**(7) Weakest step.** Row (c2) rests on one chosen definition of "size by radius" (mean over directions at the centre, with a hand-derived bending and deflection term); the bending term is tested by Wang's 14/45 (98/45 without it), the deflection term is not, and the 28/15 = 28/15 coincidence at l = 2 is unexplained.

**(8) Paper.**
- May say: "On the static tides, balance on metric data requires both the volume constraint and the extremal edge, whereas the excess on extrinsic-curvature data is fixed by the volume constraint alone, since no edge can couple to that datum."
- Must not say: "The weight v³f reproduces the fixed-volume charge", nor that the maximal slice was tested beyond a first-order stationarity argument.

Files are in `record\wave73\ruler-why\`:
- scripts\why.py
- scripts\nc_core.py (unedited copy)
- out\why.json
- out\why_plant_k1.json
- logs\00_prereg.md
- logs\01_run.log
- logs\02_run.log
- logs\03_plant_k1.log
- logs\04_outcome.md
