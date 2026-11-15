# PREREG: lane CHARGE-ID (wave 73, 2026-11-15)

Written before any computation of this lane. Question: is the reduction "area-side charge on odd-parity graviton data = c_l E_2 in the transfer variables" true on general odd data (field and momentum), l = 2, 3, 4?

## 1. Definitions, and where each comes from

Setting: d = 4, flat slice t = 0, unit ball B (radius R by dilation), N = (1 - r^2)/2, zeta = (1/2)(1 - t^2 - r^2) d_t - t x.d (wave69/final/THEOREM.md line 35). Common factor 1/8G dropped.

Data. a = (h, K): h the spatial metric perturbation, K the extrinsic curvature (K = (1/2) d_t h in TT gauge), both symmetric tensors on a neighbourhood of the closed ball. Odd sector of degree l: with phi a harmonic polynomial homogeneous of degree l, v = x cross grad(phi), K_l = sym grad v ((sym grad v)_ij = (1/2)(d_i v_j + d_j v_i)), P = x(x)v + v(x)x, every odd symmetric tensor is T = f(r) K_l + g(r) P (wave36/companionA8/NOTE.md line 83). I do not trust the record's TT parametrisation: I build T from a profile eta polynomial in r^2 by g = -eta'/((l+2) r), f = 2(eta - r^2 g)/(l-1) (NOTE.md Lemma 4) and then CHECK div T = 0, tr T = 0 exactly in Cartesian components. If that check fails the run stops.

Symplectic form. Omega(1,2) = int_B (K_1.h_2 - K_2.h_1) (NOTE.md line 19 with pi = K on traceless data).

Charge, three independent definitions, all to be computed from primitives:

- (Q1) canonical: A_can(a) = (1/2) Omega(a, delta^w a), where delta^w H = Lie_zeta H + 2 t H acts on the TT-gauge solution H_ij = h + 2tK + (1/2) t^2 Lap h + (1/3) t^3 Lap K, H_0mu = 0, and is read back as Cauchy data at t = 0 by delta h_ij = (delta^w H)_ij, delta K_ij = (1/2) d_t (delta^w H)_ij - d_(i (delta^w H)_j)0 (NOTE.md line 69). I compute the Lie derivative component by component in (t, x, y, z); I do not use the record's Lemma 3 formula. Lemma 3's closed form (4), A = int N K.K + (1/4) int N |grad h|^2 - (1/2) oint |h.nu|^2, is then a derived comparison (Q1b).
- (Q2) geometric, fixed volume: A_geo(a) = -2 (dF[h_2] + d2F[h]), F = Area(dB) - 2 Vol(B) at the fixed coordinate ball, h_2 any second-order metric obeying the N-smeared second-order Hamiltonian constraint (NOTE.md Proposition 4; THEOREM.md H4 and Section 7). Using Green's identity int_B N R_lin[q] = -2 dF[q] (verified in the run on random q) this is
  A_geo = int_B N (K.K - (tr K)^2) - int_B N R2[h] - 2 d2F[h],
  with R2[h] the O(eps^2) coefficient of the scalar curvature of delta + eps h, computed from Christoffel symbols with no TT simplification, and d2F[h] the O(eps^2) coefficient of Area - 2 Vol: volume density (1/8)(tr h)^2 - (1/4)|h|^2, area density (1/8)(tr_tan h)^2 - (1/4)|h_tan|^2 on the unit sphere.
- (Q3) scalar side: psi = x.B.x, B = -sym curl K, (curl T)_ij = eps_ikl d_k T_lj; psi_dot = x.(-sym curl K_dot).x with K_dot = -Ric_lin[h] (= (1/2) Lap h on TT h); E_w = (1/2) int_B N (psi_dot^2 + |grad psi|^2) + (w/2) int_B psi^2; c_l = 4 (l-2)!/(l+2)! (THEOREM.md line 105; NOTE.md eq. (2), line 67).

Transfer variables: (psi, psi_dot) as above (HINGE's "canonical chain field" phi with E_w, wave73/hinge/REPORT.md line 24).

## 2. Identities to be tested

- T1 (value): A_can = A_geo = c_l E_2[psi] on every odd TT datum (h, K) of degree l.
- T2 (flow, pointwise): the transferred flow is the weight-2 scalar flow: psi[delta^w K part] = N psi_dot and psi_dot[delta^w h part] = div(N grad psi) - w psi with w = 2, as polynomial identities. (psi_dot of the non-TT tensor 2NK is defined through Ric_lin, which is gauge invariant.)
- T3 (Hamiltonian): Omega_B(a, delta^w b) is symmetric in (a, b) on TT data, so delta^w is the Hamiltonian vector field of A_can on the ball.
- T4 (symplectic normalisation, only if T1 to T3 are done): on data on all of R^3 with Gaussian-weighted polynomial profiles, Omega_grav(a, b) = kappa c_l int (psi_a psi_dot_b - psi_b psi_dot_a) with one constant kappa independent of data and of l. Report kappa.

## 3. Data

- l = 2, 3, 4. Harmonics: real and imaginary parts of (n.x)^l for rational complex null vectors n, and random rational mixtures of them.
- Profiles: eta = r^(2k), k = 0..3, for the field h and for the momentum K independently: an 8-dimensional data space per harmonic. The three quadratic forms are computed as full 8 x 8 Gram matrices by polarisation in exact rational arithmetic (integrals of polynomials over the ball and the sphere by the Gamma-function formula).
- Spot checks: (a) random rational data with different harmonics (mixed m) on the h and K slots and random profiles of degree up to r^6; (b) the same integrals by Gauss quadrature in floating point, independent of the monomial formula.

## 4. Pass / fail lines

- "Equal" in exact arithmetic means the difference of the two rational numbers is exactly 0, for every Gram entry (h-h, K-K and h-K blocks). Outcome (i) needs T1 exact on all entries at l = 2, 3, 4. Any non-zero entry is reported with its block and value: non-zero only off the stationary element is outcome (ii); the control failing is a misread definition and stops the lane.
- Floating-point spot check: relative difference at most 1e-9.
- Control C1 (known answer): stationary magnetic datum (0, K_l): A = int_B N|K_l|^2 = oint|K_l|^2/((2l+1)(2l+3)) = c_l E_2[psi], and c_l E_1/A = (l+1)/(l+2) (THEOREM.md K5, line 107). Must pass exactly.
- Control C2 (record value): momentum profile eta = 1 - 3r^2/5 at l = 2: E_1/E_2 = 204/265 (THEOREM.md line 250). Informative; a miss is reported, not fatal, because the record's profile convention may differ.
- Plants, each must be caught (non-zero difference): PL1 lapse (1 - (11/10) r^2)/2 on the gravity side only; PL2 weight w = 1 instead of 2; PL3 c_l replaced by 4 (l-1)!/(l+2)!; PL4 the +2tH term dropped from delta^w.
- Gauge. Computations are in TT gauge (the record's). Odd gauge: K is invariant (a gauge move changes K by Hess xi_0, scalar type); h moves by sym grad w, w = a(r) v. Tests: G1 psi, psi_dot invariant exactly; G2 A_geo invariant exactly under h -> h + 2 sym grad w, including the move to Regge-Wheeler gauge a = -f/2... (the a that kills the K_l component, h_AB = 0); expected because w is tangent to the sphere and divergence free, so the diffeomorphism maps the coordinate ball to itself; G3 the closed form (4) applied naively off TT gauge: report whether it moves (expected: it does; then "A is read on the TT representative" is a real restriction and is stated).

## 5. What would change the flagship sentence

- T1, T2, T3 exact: step 1 closed at l = 2, 3, 4 on polynomial profiles up to r^6; say what general l and general profile need.
- T1 exact but T2 or T3 failing: the value identity holds but the flow in transfer variables is not the weight-2 flow; HINGE's operator statement is then unsupported.
- T1 failing: counter-data given explicitly.

## 6. Fences

Writes only in this folder. Native Windows Python, no WSL, no subagents, nothing committed.
