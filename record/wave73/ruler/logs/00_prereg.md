# 00 PREREG, lane RULER (wave 73). Written before the first run of ruler.py.

15 Nov 2026.

## Question
l = 2, both parities. Delta = Q_2 - Q_3 as a quadratic form on Cauchy data at t = 0:
(i) metric block (h TT, K = 0), (ii) extrinsic-curvature block (h = 0, K TT), (iii) cross block.

## Set-up fixed in advance
- Forms: null-cone's, unchanged (nc_core.py copied, not edited). Q_p = 2 oint int_0^1 v^2 f |sigma_p|^2,
  f = v(1-v), R = 1, units 1/(8G). S = Q_3.
- Data: TT gauge. Even parity: T_k = f1 Hess(phi) + f2 (x d phi + d phi x) + f3 phi delta + f4 phi xx with
  x.T.x = r^(2k) phi, k = 1..4, as h(0) (metric) and as K (extrinsic). Odd parity: O_0 = sym grad(x cross grad phi),
  O_(k+1) = a_k r^(2k+2) O_0 + r^(2k) (xV + Vx); O_1..O_4 as h(0), O_0..O_3 as K. Time evolution by nc.tt_series.
- Gram matrices are reported in the chain basis: chain field psi = x.E.x (even), x.B.x (odd); basis
  psi = r^(2j) phi and psi_dot = r^(2j) phi, j = 0..3, divided by N_Y = oint phi^2.
  Even parity: psi is metric data, psi_dot is K data. Odd parity: psi_B is K data, psi_dot_B is metric data.
- Gates that must pass before any Delta is read: every datum solves the vacuum equations; Q_3 Gram equals
  c_l E_1 (the slice entropy, block diagonal, same matrix in both parities); static electric and stationary
  magnetic give Q_2/Q_3 = 16/9; plant (v^2 f -> v^3 f in Q_2) moves 16/9.

## Expectation written before running (by hand)
Duality (E -> B, B -> -E) maps even-metric data onto odd-K data and even-K onto odd-metric, and the cone forms
are duality blind [inherited from NULL-CONE]. So in the chain basis the even and odd Gram matrices must be
IDENTICAL, and the cone rule cannot have a metric-versus-K asymmetry common to both parities: whatever it does
on the psi block it does on metric data in one parity and on K data in the other. Since the static value is
16/9, not 1, I expect: balance on NEITHER block. Prior 90%.
Unknown before the run: sign of Delta beyond the monomial profiles; whether the cross block vanishes
(the cone is not time-symmetric, so no symmetry forces it; prior 60% nonzero).

## What each outcome would mean
A. Delta = 0 on metric block, both parities, > 0 on K block: the cone ruler has the same asymmetry; the
   asymmetry is not a property of the fixed-volume ruler. (Contradicts duality + 16/9; would signal a bug first.)
B. Delta != 0 on both blocks, even and odd chain Grams identical: the cone ruler has NO metric-versus-K
   asymmetry; it has at most a psi-versus-psi_dot (profile) dependence, which is parity-swapped in h/K terms.
   The asymmetry belongs to the fixed-volume maximal-slice ruler.
   B1. Delta positive definite on each block and on the whole 8-dim space: the cone rule always over-counts.
   B2. Delta indefinite somewhere: the cone excess is profile dependent in sign; 16/9 is not a bound.
C. Delta = 0 on one block in one parity only: a parity asymmetry on the cone; contradicts duality; bug first.
D. Cross block nonzero: on the cone, metric and K data do not decouple, so "block" verdicts are basis
   statements about t = 0 only; the full-form signature is then the invariant statement.
E. Any gate fails: nothing is claimed; report the failure.

Plant fixed in advance: `--plant w3` replaces the weight v^2 f by v^3 f in Q_2 only. Must move the static ratio
off 16/9 (hand value: Q_2 picks up the factor (int v^9 f)/(int v^8 f) relative... computed by the script).
