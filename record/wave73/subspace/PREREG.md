# PREREG: lane SUBSPACE (wave 73)

2026-11-15. Written before any script of this lane exists or runs. Things done before this file: reading
(charge-id, hinge, hinge-proof, lattice-l head, wave70/modflow sections 1 and 4, Paper A = `wave36/reviseA18/NOTE.md`
sections 2-3, Hislop-Longo 1982 pages 71-76 from `sources/`), and hand algebra. No number of this lane exists yet.

## The input under test

Paper A Theorem 1 (a)-(b) with Corollary 1 (`wave36/reviseA18/NOTE.md` lines 57, 109): on each real sector
sigma = (l, parity), l >= 2, V = Phi/N_l is unitary from the helicity-2 one-particle sector onto the scalar l-wave
sector, V(H(B) ∩ H_sigma) = H_s(B) ∩ H_{s,l}, and the modular group of H(B) is the conformal flow along zeta, so
that on the odd tower K_Omega = 2 pi c_l E_1[psi] in the units of charge-id (factor 1/8G dropped), psi = x.B.x,
B = -sym curl K, psi_dot = x.curl(Ric_lin[h]).x, c_l = 4(l-2)!/(l+2)!.

## Exact statements

Conventions are charge-id's: TT data a = (h, K), Omega(1,2) = int (K1.h2 - K2.h1), h_dot = 2K, K_dot = (1/2) Lap h;
scalar Omega_s(1,2) = int (psidot1 psi2 - psidot2 psi1), psi_ddot = Lap psi.

**(a) Symplectic isomorphism onto, localisation both ways.**
- (a1) Omega = c_l Phi^* Omega_s on odd-sector data (compact support, or Schwartz on R^3).
- (a2) Phi is a bijection between regular odd TT K in sector (l, m) and scalar fields psi = G(r) r^l Y, given by a
  formula local in r in both directions: K = build_T(eta) with psi = l(l+1) eta r^l Y. Hence supp K = supp psi.
- (a3) The same for the h half modulo gauge: in Regge-Wheeler gauge h = alpha(r) (x⊗v + v⊗x), psi_dot is a local
  expression in alpha with a local inverse. Prediction (hand argument: a radial ODE operator with trivial kernel must
  have order 0): psi_dot[alpha P] = kappa_l alpha r^l Y with kappa_l a constant independent of the profile.
- Failure of (a): a kappa in a1 that is not one constant; a regular TT K with psi = 0 on a shell and K != 0 there
  (a1/a2 broken); a compactly supported psi or psi_dot in B whose preimage needs support outside B (onto fails as a
  statement about the ball); a differential operator of positive order with a regular non-gauge kernel in a3.

**(b) Vacuum two-point function.** With J = T/|T| the vacuum complex structure, mu(a,b) = Omega(a, J b):
 mu_grav(a,b) = int [2 K_a.K_b / omega + (1/2) omega h_a.h_b],  mu_s = int [psidot_a psidot_b/omega + omega psi_a psi_b],
 and the claim is mu_grav = c_l mu_s on Phi-related data, block by block (K-K and h-h; no cross block on either side).
- Failure of (b): any Gram entry with mu_grav != c_l mu_s; or agreement that depends on the profile (a kappa(omega)).

**(c) Modular generator.** For the free massless scalar in 3+1 dimensions the vacuum modular group of the unit
double cone is the conformal flow of a dimension-1 field (Hislop-Longo 1982, Theorem 2), whose generator on Cauchy
data is 2 pi (N Pi, div(N grad Phi) - w Phi) with w = 1, the Hamiltonian flow of 2 pi E_1. psi obeys the massless
wave equation (charge-id T5), so it is that field.
- Failure of (c): the multiplier in Hislop-Longo's Theorem 2 does not give dimension 1 (w != 1); or the dimension-1
  flow fails to preserve the vacuum two-point function; or the conformal flow (1) of Paper A on (E, B) does not
  induce w = 1 on (psi, psi_dot).

## Tests, with pass/fail lines (all exact rational arithmetic unless marked)

Script `scripts/subspace.py L [KMAX]`, L = 2, 3, 4 (5 as an extra), KMAX = 3 (profiles s^0..s^3, s = r^2).
It imports `charge-id/scripts/polyalg.py` read-only and re-implements charge-id's `KP`, `fg_from_eta`, `build_T`,
`psi_of_K`, `psidot_of_h` by copy (charge_id.py runs at import).

A (localisation, unit ball, polynomial profiles):
- A1 PASS iff psi[build_T(s^k)] - l(l+1) s^k phi == 0 for k = 0..KMAX, and build_T is TT.
- A2 PASS iff the null space of {a(s) K_l + b(s) P : div = 0}, deg a, b <= 3, has dimension exactly 4 and equals the
  span of build_T(s^k), k = 0..3 (no TT odd tensor outside the image of the local inverse).
- A3 PASS iff psi_dot[s^k P] == kappa_l s^k phi with one constant kappa_l != 0 for k = 0..KMAX. If it fails, report
  the operator and its kernel; the lane does not stop.
- A4 PASS iff for eta_a = (1-s)^3 s^a (data vanishing to second order at the boundary) int_B K_a.h_b ==
  -c_l int_B psi_a psidot_b for all a, b, with h TT and again with h in Regge-Wheeler gauge. INFO only: the same with
  eta = s^a (boundary-touching), where a mismatch is expected and its size is reported.

B (two-point function, R^3, profiles q(s) e^{-s}, Fourier transforms exact as polynomial times e^{-k^2/4}):
- B0 control (Parseval, known answer): position-space int f_a f_b == momentum-space value, all tensors and scalars
  used. If this fails the lane stops: the transform is misread.
- B1 control (known answers from the record): Omega_grav = kappa c_l Omega_s with kappa = 1 on all pairs (charge-id
  T4), and int K_a.K_b = c_l (1/2) int grad psi_a . grad psi_b (Paper A Prop. 1 energy identity, momentum half).
- B2 PASS iff GK_ab == c_l SK_ab and Gh_ab == c_l Sh_ab for all a, b = 0..KMAX, all four matrices non-zero, where
  GK = int 2 K_a.K_b/omega, SK = int omega psi_a psi_b, Gh = int (1/2) omega h_a.h_b, Sh = int psidot_a psidot_b/omega.
- Plants (each must FAIL the B2 line): PB1 c_l -> c_l (l+3)/(l+2); PB2 the symplectic but not vacuum-preserving map
  (h, K) -> (2h, K/2), for which the script first confirms B1's Omega identity still holds; PB3 omega^{+1} and
  omega^{-1} exchanged on the scalar side.

C (weight):
- C1 PASS iff, on polynomial odd TT data in the ball, with B = -sym curl K, E = Ric_lin[h], dE/dt = -sym curl B,
  dB/dt = sym curl E and Paper A's flow (1): (C1a) dE/dt, dB/dt agree with the data's own evolution; (C1b)
  x.dB.x == N psi_dot; (C1c) x.(sym curl dE).x == div(N grad psi) - w psi with w = 1; (C1d) dE, dB are TT.
  Plant PC1: w = 2 in C1c must FAIL.
- C2 (sympy, `scripts/hl_gamma.py`) PASS iff d gamma/d lambda at 0 equals -3 x_0 for Hislop-Longo's multiplier and
  the induced field transformation is zeta.d phi + (1/4)(d.zeta) phi (dimension 1), zeta the diamond's conformal
  Killing field up to the sign of the flow parameter. Plant: dimension 2 must not match.
- C3 PASS iff on scalar Gaussian data in sector l the generator delta_w(Phi, Pi) = (N Pi, div(N grad Phi) - w Phi)
  is skew for mu_s (mu(delta a, b) + mu(a, delta b) == 0 on the whole Gram) at w = 1. Controls (known answers):
  time translation (Pi, Lap Phi) and dilation (x.grad Phi + Phi, x.grad Pi + 2 Pi) are skew for mu_s; delta_w is
  skew for Omega_s at w = 1 and w = 2. Plant PC3: w = 2 must FAIL skewness for mu_s. If a control fails, stop.

## What is not tested here, by design
- General l beyond 5 and non-polynomial profiles: rest on the fibre argument (Schur plus dilation) in Paper A.
- The lattice modular Hamiltonian: already measured (wave70/modflow S4, both blocks, l = 2, 3; LATTICE-L, field
  block on the static profile, l = 2, 3, 4). Not repeated; the report says what those do and do not test.
- Second quantisation (one-particle modular group to the von Neumann algebra), even parity, d >= 5.

## Verdict rule
Each of (a), (b), (c) gets PROVED / STANDARD / NUMERIC ONLY / GAP / WRONG. A verdict of PROVED needs a written
argument that does not lean on the tests; the tests are then confirmation at l = 2..5. Any failed non-plant line is
reported with its output, and the verdict for that part cannot be better than GAP until explained.
