# EVEN-BLIND pre-registration

Lane: EVEN-BLIND, 2026-11-15. Written before any number for the even-tower quantity
had been computed (no pilot, no test run of the even perturbation). Read before writing: only
`..\hinge-blind\` (PREREG.md, both addenda, REPORT.md, hb.py, run.py, continuum.py, validate.py,
addendum_rule.py, two timing logs). Not opened, listed or searched: `hinge-e\`, `even-tower\`,
`gen-law\`, `hinge\`, or any other folder. The claimed even-tower value is unknown to this lane.
The odd-tower value 7 zeta(3)/(16 pi^2) = 0.0532848 is known (it is in hinge-blind's REPORT.md).
No closed form for the even case has been derived or attempted.

## Quantity

Free massless scalar, 3+1 flat dimensions, sector l, radial field phi_l = r x (radial part),
momentum pi_l conjugate to it, H_l = 1/2 Int dr [pi^2 + (d_r phi)^2 + l(l+1) phi^2/r^2].
rho_l = vacuum reduced to the unit ball, rho_l ∝ exp(-K_l). Perturbation

    Z_even = (g/2) <Pi, G_D Pi>_ball = (g/2) Int_0^1 dr Int_0^1 dr' pi_l(r) G_l(r,r') pi_l(r'),

G_l = inverse of D_l = -d^2/dr^2 + l(l+1)/r^2 on (0,1) with U(0) = 0, U(1) = 0.
(Derivation: 3-d u = U(r) Y/r, Pi = pi_l(r) Y/r; -Lap u = Pi becomes D_l U = pi_l, and
Int d^3x Pi u = Int dr pi_l U. D_l is the same operator whose lattice form is K below.)
rho'_l ∝ exp(-K_l - 2 pi Z_even). Measured:

    A_l = tr(rho_l^{1/2} rho'_l^{1/2}),   Q^E_l = l(l+1)(1 - A_l)/g^2,   also l^2 (1 - A_l) at g = 1.

## How Z_even enters

Lattice variables as in hinge-blind: r_j = j a, phi_j = phi(r_j), pi_j = a pi(r_j),
H = (1/2a)(pi.pi + phi.K.phi), K_mod = 1/2 (pi.N.pi + phi.M.phi) from the vacuum covariances.
Discretising D_l U = pi as (K_B/a^2) U = pi_j/a gives U = a K_B^{-1} pi, and
Int dr pi U = Sum_j pi_j U_j, so

    Z_even = (g/2) a pi^T G pi,   G = K_B^{-1}  (ball sites only, zero elsewhere),
    N' = N + 2 pi g a G,   M' = M.

(Odd case for comparison: M' = M + 2 pi g a diag(w), N' = N.) The perturbed Gaussian state is
obtained from the Williamson form of (M, N') with hb.williamson; the affinity with hb.affinity_I,
unchanged. If N' is not positive definite (possible for g < 0) the state does not exist on that
lattice and that is reported.

Second order in g: the canonical map Q -> Pi, Pi -> -Q leaves K = Sum eps (Q^2 + Pi^2)/2 invariant,
so hb.second_order's formula holds with its matrix replaced by Gt = A^{-1} G A^{-T}
(pi = A^{-T} Pi). Q2 = l(l+1) lim (1-A)/g^2 is computed that way and checked against exact
affinities at small g.

## Discretisation of G_D and boundary conventions

K is hinge-blind's primary "D2": K_jj = 2 + l(l+1)/j^2, K_{j,j+1} = -1. Region A = sites 1..n.
- PRIMARY "ghost": a = 1/(n + 1/2) (surface midway between n and n+1, as hinge-blind's primary).
  K_B = leading n x n block of K with K_B[n,n] increased by 1 (antisymmetric ghost
  U_{n+1} = -U_n: the zero of U sits on the surface r = 1 to second order). All n sites carry pi.
- Alternative "trunc": same a; K_B = leading n x n block unchanged (U_{n+1} = 0, zero at
  r = 1 + a/2, a first-order boundary error).
- Alternative "site": a = 1/n, surface on site n, U_n = 0: K_B = leading (n-1) x (n-1) block,
  G has a zero row and column at site n.
Chosen as primary before seeing any number, on the a-priori ground of second-order accuracy with
hinge-blind's primary geometry. All three must reach the same continuum value; if they do not
within brackets, that is the first caveat of the report.
N_tot = 4n. Arithmetic: python-flint arb, precision from hinge-blind's default_prec (>= 220
digits beyond the smallest nu^2 - 1/4), rows rejected if the margin is below 100 digits.

## Lattices

n in {24, 32, 40, 48, 56, 64, 80, 96, 112, 128, 144, 160, 176, 192};
l in {2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24, 32, 40, 48}; only rows with n >= 2l.
Per row: Q2 for ghost, trunc, site; exact affinity for ghost at g = +1 and g = -1.
Larger n or l only if cost allows; any such row is used by the same rules.

## Series

q2 (primary, the g -> 0 coefficient), p1 (g = +1), m1 (g = -1) and mean of +-1 where N' > 0,
p1L = l^2 (1 - A) at g = 1 (= p1 x l/(l+1)).

## Step 1: n -> infinity at fixed l (hinge-blind's addendum rule, unchanged)

Lattices n >= max(24, 2l), the eight largest. Central: exact fit in
{1, x, x^2 ln x, x^2, x^3 ln x, x^3, x^4, x^4 ln x} truncated to the number of lattices, x = a.
Bracket: largest deviation among drop-smallest, drop-largest, pure polynomial basis, basis led by
x ln x. An l is usable when bracket <= 1e-3 x |central| (hinge-blind's 5e-5 on 0.05, made relative
because the scale here is unknown). If the series does not converge in n (bracket > 10 percent at
every l) no continuum limit is claimed and d Q/d ln n is reported.

## Step 2: l -> infinity. Competing forms, lam = l + 1/2, u = 1/lam

Constant family (limit q):
  C1 q + c1 u;  C2 + c2 u^2;  C3 + c3 u^3;
  CL1 q + c1 u + b u ln lam;  CL2 + c2 u^2;  CL3 + c3 u^3.
Zero family (limit 0):
  ZP c lam^s (s fitted, required s < 0);  ZP2 lam^s (c + d u);
  ZL1 c1 u + b u ln lam;  ZL2 + c2 u^2;  ZL3 + c3 u^3;
  ZG c/ln lam + d/ln^2 lam   (decays, but too slowly for the sum to converge).
Growth family:
  GL q + b ln lam + c1 u;  GP = ZP/ZP2 landing on s > 0.
Linear forms by least squares; ZP/ZP2 by a scan in s with linear least squares in the rest.

Deciding criterion: out-of-sample prediction. Splits (usable l only):
  S1 train 3 <= l <= 12, test l >= 16;   S2 train 3 <= l <= 20, test l >= 24;
  S3 (secondary) train 4 <= l <= 16, test l >= 20.
A form is fitted when train points >= parameters + 1 (CL3 on S1 has 7 points for 5 parameters).
A form PASSES a split when, at every test l, |prediction - Q_l| <= 2 x (step-1 bracket at that l).
Its score is the largest miss/bracket. Done for series q2 (deciding) and p1 (reported).

Gates:
- "finite non-zero constant" (sum over sectors diverges logarithmically): (a) at least one
  constant-family form passes S1 and S2; (b) no zero-family or growth-family form passes both;
  (c) the limit, defined as the median over [passing constant forms] x [windows l >= 3, 4, 6, all
  usable l] x [series q2, p1, and mean if it exists] with bracket = full range, excludes zero and
  has half-range <= 20 percent of the median.
- "tends to zero": a zero-family form passes S1 and S2 and either no constant-family form passes
  or the passing constant forms' range of q includes zero. Convergence of the sum then follows
  from the fitted decay (power s < 0 or u ln lam: converges; ZG: still diverges, as ln ln).
- "grows": a growth-family form passes S1 and S2 and no constant-family form passes.
- anything else: undetermined, with the full score table. The gate is not edited after data.
If fewer than two test points exist in a split, that split is reported as not evaluable and the
gate is undetermined.

Consequence: Q -> c finite non-zero gives Sum_{l<=L}(2l+1)(1-A_l) = 2 c g^2 ln L + const.

## Controls (all required)

(i) Odd pipeline intact: Z_odd through this lane's copy at (n,l) = (24,8), (32,4), (48,16),
    g = +-1, primary 'mid' set, compared digit by digit with hinge-blind's stored raw JSON.
    If any differs beyond the 28th digit: stop and report.
(ii) g -> 0: exact 1 - A at g = +-1/4, +-1/16 at n in {24, 48, 96}, l in {2, 8, 32}; A -> 1 and
    (1-A)/g^2 -> Q2/(l(l+1)).
(iii) Local momentum contrast: Z = (g/2) Int_0^1 N(r) pi^2 dr, N = (1 - r^2)/2, i.e.
    N' = N + 2 pi g diag(N(r_j))/a, a = 1/(n+1/2), g = +-1, +-1/4, same n, l as (ii) plus
    n = 192 if cheap. Reported: 1 - A against n and l. Expectation stated before running
    [reasoning only]: this rescales the momentum term of K_l itself, a temperature-type change
    felt by every Rindler mode, so -ln A should grow like ln n at fixed l (no continuum limit),
    unlike the non-local term.
(iv) Outside size: N_tot = 8n against 4n (and 2n) at n in {24, 48, 96}, l in {2, 4, 8, 32}.
(v) Planted errors at n in {24, 48, 96}, l in {2, 8, 32}: Neumann instead of Dirichlet at the
    surface (K_B[n,n] decreased by 1); factor 2 (Z without the 1/2); centrifugal term dropped from
    G_D. Each must move Q well outside the step-1 bracket.
Also: hinge-blind's validate.py re-run on the copy; a brute-force Fock-space check of a
momentum-block perturbation (two modes) for affinity_I and for the second-order formula; precision
reruns at 1.4x bits on (24,8), (48,16), (96,32).

## Prior (weakly held, reasoning only)

Planar scaling: near the surface N ~ x, while G_D acts as ~1/k^2 on functions of width 1/k, so the
relative size of the perturbation is ~ g/l and 1 - A ~ g^2/l^2: a finite non-zero Q^E_inf is
expected, of order 0.01 to 1. For g = -1 positivity of N' is not guaranteed at small l.
