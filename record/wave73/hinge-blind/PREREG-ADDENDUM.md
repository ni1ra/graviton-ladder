# HINGE-BLIND addendum: rule for the final lattice extrapolation, fixed mid-run

Written 2026-11-15 while the n = 192 rows were still computing. PREREG.md is unchanged and its
rules are applied and reported as written. This file fixes, before the remaining rows are
analysed, the second analysis that will be reported next to it, and records what was known when
it was fixed. The sealed folder is still unopened and the predicted constant is still unknown to
this lane.

## What was known when this was written

- Lattice rows for n = 24, 32, 48, 64, 96, 128 at l = 2..32 (primary and boundary-site sets),
  Srednicki-form rows to n = 96, planted-error rows, second-order coefficient checked against the
  exact affinities at small g.
- The PREREG extrapolation in n (three-parameter fits) wastes the precision of the rows: exact-fit
  Richardson extrapolation in the basis {1, x, x^2 ln x, x^2, x^3 ln x, ...} is stable to about
  1e-6 at l = 2 and 1e-4 at l = 16, and the two boundary conventions ('mid', 'site') reach the same
  continuum values that way (l = 2: 0.0376997 against 0.0376972).
- The PREREG forms in l (powers of 1/l or 1/(l(l+1))) drift with fit order and window (0.0521 to
  0.0532 so far). A far-region argument, checked on the lattice (far_region.py), gives a term
  -(ln l)/(64 l) in Q_l; with that term the fits on l <= 16 gave 0.05328 to 0.05345, and a fitted
  log coefficient of -0.0155 on one window against -1/64 = -0.0156.
- After seeing those fits, this lane derived a closed form of its own for the large-l (planar,
  second-order) limit: with Rindler modes K_{i omega}(k rho),
  Q_inf = (1/32) Int_0^inf tanh^2(y)/y^2 dy = 7 zeta(3)/(16 pi^2) = 0.053285.
  That value is this lane's own derivation, not read from anywhere, and it was obtained before the
  rule below was fixed. It is therefore reported separately from the lattice number, and the rule
  below contains no reference to it.

## Rule for the second (post hoc) lattice analysis

Step 1, n -> infinity at fixed l. For each l use the lattices with n >= max(24, 2l), the eight
largest of them if there are more. Central value: exact fit in the basis
{1, x, x^2 ln x, x^2, x^3 ln x, x^3, x^4, x^4 ln x} (truncated to the number of lattices), x = a.
Bracket: largest deviation among (smallest lattice dropped), (largest lattice dropped), (pure
polynomial basis), (basis {1, x ln x, x, x^2 ln x, ...}), same lattices.

Step 2, l -> infinity. Use every l whose step-1 bracket is below 5e-5. Four series: second-order
coefficient Q2, mean of g = +-1, g = -1, g = +1. With lam = l + 1/2, u = 1/lam, least-squares
fits of
  logfix-3, logfix-4 :  Q_inf + sum_{j<=k} c_j u^j - (1/64) u ln(lam)
  logfix2-3, logfix2-4: the same plus c' u^2 ln(lam)
  logfree-3          :  Q_inf + sum_{j<=3} c_j u^j + b' u ln(lam), b' fitted
on windows l >= 3, l >= 4, l >= 6 (a fit is kept when it has at least one more point than
parameters). Central value: the median of all kept fits over the four series. Bracket: their full
range. The pure-polynomial fits (PREREG forms) are reported beside it, not folded in.

The first line of the lane's report quotes this lattice number and bracket. The closed form above
is quoted after it as a separate, theory-side statement by this lane.
