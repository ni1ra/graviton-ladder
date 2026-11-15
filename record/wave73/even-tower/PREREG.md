# PREREG: wave73 / even-tower (the even-parity half of the graviton on the ball)

Lane EVEN-TOWER, Sunday 15 November 2026. Written before any script of this lane was written or run.
Folder: `record/wave73/even-tower/`. No report file is written by this lane; the report text goes back to the lead.

## Disclosure

- Read before writing this: `charge-id/REPORT.md`, `jss-id/REPORT.md`, `subspace/REPORT.md`, `gen-law/REPORT.md`,
  `hinge-proof/REPORT.md`, `ss-second/REPORT.md` Section 5, companion A8 (`wave36/companionA8/NOTE.md`) Sections 2, 3,
  Proposition 4 and Section 5, `hinge-e/PREREG.md` and the tails of its three logs, and the scripts `polyalg.py`,
  `charge_id.py`, `subspace.py`. So this lane is not blind: it knows GEN-LAW's closed form for the even momentum half,
  companion A8's Theorem 3 for the even field half and HINGE-E's Shale constant.
- The predictions below were derived by hand in this session (integration by parts, written out in the lane's final
  message) before any code. The scripts check hand algebra.
- A folder `hinge-e/` already holds a Shale test of the even tower that takes GEN-LAW's form on trust. It has no report
  on file. What it took on trust is exactly what this lane computes from the geometric definition.

## Conventions (charge-id's, unchanged)

d = 4, unit ball B, t = 0, N = (1 - r^2)/2, factor 1/8G dropped, integrals in units of 4 pi. Data a = (h, K), K the
extrinsic curvature, K_ij = (1/2)(d_t H_ij - d_i H_0j - d_j H_0i). Omega(1,2) = int_B (K_1.h_2 - K_2.h_1).
c_l = 4 (l-2)!/(l+2)!. phi = r^l Y a harmonic polynomial of degree l, s = r^2.

## Definitions

**D1. Even TT tensors (companion A8 Lemma 7, n = 3).** T[q] = f1 Hess(phi) + f2 P_phi + f3 phi delta + f4 phi X,
P_phi = x (x) grad phi + grad phi (x) x, X = x (x) x, with (primes are d/ds)
beta = -2 q'/(l+1), alpha = (q + 2 s q'/(l+1))/l, f2 = (beta - 4 alpha')/(l+2), f3 = -2 alpha' - (l+1) f2,
f4 = -(2 l f2 + 3 f3)/s, f1 = (alpha - s f2)/(l-1). Then x.T.x = q phi. Basis: q = s^k, k = 0..KMAX; q = 1 is the
gauge mode Hess(phi)/(l(l-1)). The script checks div T = 0, tr T = 0, x.T.x = q phi exactly and does not trust the lemma.

**D2. Even gauge moves.** A diffeomorphism X = (X^0, xi) acts on first-order data by
h -> h + 2 sym grad xi, K -> K + Hess(X^0) (derived from D_mu X_nu + D_nu X_mu and the definition of K above).
Sector vectors: xi = a(s) grad phi + b(s) phi x, X^0 = c(s) phi. Normal shift of the edge sigma = xi.nu at r = 1,
corner time shift tau = X^0 at r = 1, slice deformation X^0 in the bulk.

**D3. The geometric charge (Q1).** F(S, Sigma) = Area(S) - 2 Vol(Sigma), Sigma a slice with edge S, in the metric
eta + eps H + eps^2 H_2. At the coordinate ball,
  A_un(h, K) = -2 [eps^2 coefficient of F] = int_B N (K.K - (tr K)^2) - int_B N R^(2)[h] - 2 delta^2 F[h],
the second-order metric eliminated by the N-smeared second-order Hamiltonian constraint (companion A8 Proposition 4;
charge-id Q2). No TT simplification: R^(2) from Christoffel symbols, tr K kept.
*Moving the region is a gauge move of the data.* By diffeomorphism invariance, F of the region displaced by eps X in
the metric g equals F of the coordinate ball in the pulled-back metric, whose first-order data are D2 applied to (h, K)
and whose second-order data again drop out through the smeared constraint. So the second-order value at every
displaced region is A_un on the gauge orbit, and
  **A_ext(a) = stationary value over all sector gauge moves g of A_un(a + g).**
Stationary = saddle (minimum in sigma and in the bulk slice, maximum in tau): the extremal edge on the maximal slice.
Fixed true volume: with no l = 0 part at first order, Area at Vol = Vol_0 and Area - 2 Vol have the same second-order
coefficient (the l = 0 shift is O(eps^2) and F_0 is stationary under it).

**D4. Chain variables.** psi_E = x.Ric_lin[h].x, psidot_E = x.Ric_lin[2K].x (companion A8 Section 5, Corollary 3).
E_w = (1/2) int_B N (psidot^2 + |grad psi|^2) + (w/2) int_B psi^2. G_D = inverse Dirichlet Laplacian of the unit ball,
-Lap(G_D f) = f, G_D f = 0 at r = 1.

## Hand predictions (frozen)

- H1. A_un is unchanged by every gauge move with sigma = 0 and X^0 = 0 (xi tangent to the edge, anything inside).
  It depends on xi only through sigma.
- H2. On pure gauge data the smeared-constraint formula must return elementary flat-space geometry:
  A_un(2 sym grad xi, 0) = -(l-1)(l+2) oint sigma^2, and
  A_un(0, Hess X^0) = oint |grad_tan tau|^2 - 2 int_B |grad X^0|^2.
- H3. Drives. For TT K: A_un(K + Hess X^0) = int N|K|^2 + 2 oint kappa tau + oint|grad_tan tau|^2 - 2 int|grad X^0|^2,
  kappa = x.K.x. Stationary point: X^0 harmonic, tau = -kappa/(l(l-1)) on the edge (so the moved K has zero
  normal-normal component there; the bulk condition is tr K = 0, the maximal slice). Value:
  **A_ext^K = int_B N|K|^2 - oint kappa^2/(l(l-1))** on TT K.
  For TT h (companion A8 Lemma 8): A_ext^h = A_un^h + oint G^2/((l-1)(l+2)), G = (1/2)(d_r + 2)(x.h.x) at r = 1.
- H4. No sigma-tau cross term and no h-K cross term, so A_ext(h, K) = A_ext^h(h) + A_ext^K(K).

## The four questions as testable identities

**Q1 (exists, gauge invariant).**
- Q1a. The stationarity system is solvable and its value unique (rank test: the gauge family's form has rank exactly
  1 in the xi directions; non-degenerate in the X^0 directions).
- Q1b. A_ext(a + g) = A_ext(a) exactly for random even gauge moves g (xi with a, b of degree 2 in s; X^0 with c of
  degree 2), on a random datum with h and K both non-zero, computed by the same pipeline from the moved (non-TT) data.
- Q1c. The stationary point found in a finite polynomial family is stationary against extra directions outside it.
- Q1d. Non-vacuity: A_un itself is not invariant under the same moves.
- FAIL if Q1b fails: then report what the residual dependence is.

**Q2 (what it is).** On the Gram matrix of the basis D1 (k = 0..KMAX on the h slot and on the K slot), exact:
- Q2a. Field block: A_ext^h = c_l E_1[psi_E] = c_l [(1/2) int N |grad psi_E|^2 + (1/2) int psi_E^2]. (Companion A8 Thm 3.)
- Q2b. Momentum block: A_ext^K = c_l [(1/2) int N psidot_E^2 + (1/2) int psidot_E G_D psidot_E]. (GEN-LAW item 1.)
- Q2c. h-K block zero.
- Q2d. Spectrum of the momentum block against c_l (1/2) int N psidot^2: 1 + 1/((j+1)(2l+2j+3)), j = 0, 1, 2
  (the record's "not one weight"). Field block against c_l E_1: all 1.
- Q2e. Flow level: with Kt = K + Hess X^0* the stationary representative, x.Ric_lin[2 N Kt].x = N psidot_E + G_D psidot_E.
- Outcome bins: ONE-WEIGHT (momentum block equals c_l E_w for some w: impossible unless the G_D term vanishes, since E_w
  has no w on that block); NON-LOCAL (Q2b holds); OTHER (neither: report the form found).

**Q3 (vacuum side, subspace-lane method on even data).**
- Q3a. On R^3, Gaussian profiles: int K_a.h_b = kappa c_l int psidot_a psi_b with one constant kappa = +1.
- Q3b. Two-point function, whole Gram, both blocks: int (1/2) omega h_a.h_b = c_l int omega psi_a psi_b and
  int 2 K_a.K_b/omega = c_l int psidot_a psidot_b/omega (momentum space, exact). Parseval control first.
- Q3c. Paper A's flow (1), dE = N E_dot + 2 sym(x cross B), dB = N B_dot - 2 sym(x cross E), on even data:
  x.dE.x = N psidot_E and -x.(sym curl dB).x = div(N grad psi_E) - 1 psi_E (weight 1); dE, dB stay TT.
- Q3d. Locality both ways: E = T[q_E] with q_E phi = psi_E (algebraic), and B = -sym curl K is the odd-type TT tensor
  with profile eta = F/(l(l+1)), psidot_E = F phi (algebraic). psi_E obeys the wave equation (psi_E[Lap h] = Lap psi_E[h]).
- If Q3a-c hold, the vacuum generator is 2 pi c_l E_1[psi_E, psidot_E] by the subspace lane's argument (Hislop-Longo
  for the scalar), unchanged.

**Q4 (consequence).** If Q2a-b and Q3 hold: A_ext - K_Omega/2pi = Z_E = (c_l/2) <psidot_E, G_D psidot_E>, on the
momentum half only, positive, non-local, in canonical variables Z_E = (1/2)<pi, G_D pi> with coefficient exactly 1.
That is the operator HINGE-E assumed. Its Shale test is then inherited (and its structure script re-run with output
in this folder); this lane does not re-derive the constant.

## Controls (a failure stops the run: a definition is misread)

- C1. ODD data through the same code path: odd TT tensors (charge-id's build_T) at degree l and at degree l+1 (same
  polynomial parity as the even gauge family, so the vanishing is not by inversion parity) give zero drive against
  every even gauge direction, and A_ext = A_un = c_l E_2[psi_B] on the odd Gram (h and K slots).
- C2. Stationary electric datum h = -2 phi delta, K = 0 (not TT), through the full pipeline:
  A_ext = 2 l (l-1)/((l+2)(2l+3)) oint Y^2 = c_l E_1[psi_E], psi_E = l(l-1) phi: ratio information/area = 1
  (companion A8 Proposition 3(b)).
- C3. Green identity int N R_lin[q] = -2 delta F[q] on a random symmetric q (charge-id D1).
- C4. Companion A8 Corollary 2(ii): unextremised A_un(T_0, T[q]) = -(l/4)(q'(1) + (l+2) q(1)) oint Y^2 with T_0 = Hess phi
  (r-derivative), non-zero; extremised pairing zero.

## Plants (each must be CAUGHT)

- PL1 (extremisation). Corner time shift left out (tau frozen at 0, slice still free): Q1b must fail under X^0 moves.
- PL2 (extremisation). Slice not extremised (X^0 restricted to harmonic c phi only is the true answer; the plant
  instead drops the bulk term -2 int|grad X^0|^2 from the stiffness): Q2b must fail.
- PL3. Lapse (1 - 1.1 r^2)/2 in A_un: Q2a must fail.
- PL4. Weight 2 on the even field block: Q2a must fail.
- PL5. Neumann-type inverse (constant fixed by d_r u = 0 at r = 1) in place of G_D, and coefficient 2 in place of 1: Q2b must fail.
- PL6. Vacuum side: (h, K) -> (2h, K/2) keeps Q3a and breaks Q3b; weight 2 in Q3c fails.

## Scope

l = 2, 3, 4, KMAX = 3, exact rationals, native Windows Python. polyalg imported read-only from charge-id.
Not attempted: symbolic l; d >= 5; a fresh derivation of the Shale constant; existence of the second-order metric.

## Addendum A (written after `logs/run01_*` had passed, before `scripts/shale_check.py` was written or run)

The main run confirmed the operator HINGE-E assumed. Its Shale constant is still inherited. One independent number is
added, own code, sharing nothing with HINGE-E's scripts except the model:

- Sector l in the hyperbolic frame, grid 0 < rho <= 50, step 0.025, Dirichlet at both ends (HINGE-E used a Neumann wall).
  M0 = -d^2/drho^2 + l(l+1)/sinh^2 rho, N = 1/(1 + cosh rho), W = N^(1/2) (M0 + 1/4)^-1 N^(1/2).
  H_0 = (1/2) pi.pi + (1/2) phi.M0.phi, H_1 = (1/2) pi.(1 + g W).pi + (1/2) phi.M0.phi, each thermal at beta = 2 pi.
- Affinity A_l = tr(rho_0^(1/2) rho_1^(1/2)) from covariance matrices:
  A = prod coth(pi w/2)^(1/2) over both spectra / sqrt(det(V_0 + V_1)), V_i the beta = pi thermal covariance.
- S1 control: g = 0 gives |log A| < 1e-8.
- S2 control (known answer): the same routine with the odd tower's potential (H_1 = H_0 + (g/2) phi.N.phi) reproduces
  HINGE-E's stored odd-control number 1 - A = 1.6919e-4 at l = 16, g = 1, within 2 %.
- S3: 1 - A_l at g = 1 within 2 % of HINGE-E's stored values at l = 16, 32, 64, 128.
- S4: l^2 (1 - A_l) (mean of g = +1, -1) is increasing in l over 16..128 and lies in [0.6, 1.05] a_E, a_E = 0.0135994;
  so sum (2l+1)(1 - A_l) diverges logarithmically.
- Plant: W -> N (the odd tower's operator put in the kinetic slot) must miss S3.
If S2 fails the routine is wrong and nothing is claimed from it.

## Addendum B (after `logs/run03_shale_check.log`, where control S2 FAILED as written: 1.535e-4 against 1.692e-4)

Per Addendum A nothing is claimed from run03. Diagnosis offered before any re-run: the Dirichlet wall puts the box
frequencies at p = k pi/L_box (an endpoint rule that drops the half cell at p = 0, where the thermal weights are
largest), HINGE-E's Neumann wall puts them at (k - 1/2) pi/L_box (a midpoint rule); with L_box about 46 the expected
deficit is of order 5 to 10 %. Two tests of the diagnosis, lines fixed here:
- B1. Same routine, Neumann wall at rho = 50 (last diagonal entry 1/h^2): S2, S3, S4 and the plant as in Addendum A, unchanged thresholds.
- B2. Dirichlet wall with the box doubled (rho = 100, step 0.05 so the matrix size is unchanged), odd control at l = 16:
  the deficit against HINGE-E's number must fall to between 0.3 and 0.7 of run03's deficit (9.3 %).
If B1 fails, the routine is wrong and the Shale constant stays inherited only. The script takes the wall as an argument;
run03's log is kept.

## Addendum C (exploratory, written after runs 01 to 06, before `scripts/even_canonical.py` was written or run)

Charge-id's Step 9 says the canonical construction (1/2) Omega(a, delta^w a), delta^w H = Lie_zeta H + 2tH on the
TT-gauge solution, does not carry over to even data. Make that concrete on the even TT basis (delta_w copied from
charge_id.py). Expectations, by hand from companion A8 Lemma 3 ("for every TT datum") and this lane's H3:
- K1. Omega_B(a, delta^w b) is symmetric on even TT data and equals the closed form (4):
  int N K.K + (1/4) int N |grad h|^2 - (1/2) oint |h.nu|^2.
- K2. A_can - A_ext = + oint kappa_a kappa_b/(l(l-1)) on the K slot and
  - (1/4) oint (nu.h.nu)_a (nu.h.nu)_b - oint G_a G_b/((l-1)(l+2)) on the h slot.
- K3. A_can is not invariant under the residual TT gauge (the row of T[1] = Hess phi/(l(l-1)) is non-zero in both slots),
  so on even data it is not a function of the curvature; A_ext is.
- K4 (info only, no pass line): what delta^w does to (psi_E, psidot_E).
