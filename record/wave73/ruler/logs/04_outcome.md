# 04 Outcome against 00_prereg (lane RULER)

Runs: 01 (script error, see note at its end), 02 (clean, 68 checks, 0 fail, 75 s), 03 (plant w3, caught).

Landed: outcome B1 + D of 00_prereg.
- Delta = Q_2 - Q_3 is nonzero and positive definite on the metric block and on the extrinsic-curvature block,
  in both parities (exact leading minors all positive), and positive definite on the full 8-dim space.
- Even and odd chain-basis Gram matrices are identical (duality), so in h/K terms the two blocks swap between
  parities: psi block = metric (even) = extrinsic (odd); psi_dot block = extrinsic (even) = metric (odd).
- Cross block (psi x psi_dot) is NOT zero.
- Static reproductions: 16/9 for the static electric tide (Newtonian gauge and TT-gauge Gram entry agree) and
  for the stationary magnetic tide. psi_dot = phi gives 7/4 (null-cone's scalar momentum rung and vector field rung).
- Plant w3 (v^3 f in Q_2): static value 16/9 -> 4/3, psi_dot = phi value 7/4 -> 1. Caught.
  Unexplained and not pursued: those two planted numbers equal the odd tower's fixed-volume values.
- Spectrum of Q_2 relative to S on the 4-dim blocks: psi 1.0500, 1.1394, 1.4280, 3.1354; psi_dot 1.0383,
  1.0950, 1.2900, 3.2173 (floats from exact characteristic polynomials). 16/9 is neither a lower nor an upper bound.

Not done: l > 2; any basis beyond r^(2j) phi, j <= 3; no proof of positivity beyond this 8-dim space.
