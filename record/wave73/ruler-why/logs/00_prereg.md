# 00 PREREG, lane RULER-WHY (wave 73). Written before the first run of why.py.

15 Nov 2026.

## Disclosure
Before writing this I read the sibling reports and opened `ruler\out\ruler_plant_w3.json` (RULER's own plant
run). It already shows that at l = 2 the form Q_2[v^3 f] - S is NOT zero off the [psi_dot_0, psi_dot_0] entry
and has a nonzero cross block. So the "form identity" half of the v^3 f question was seen before this file;
the l = 3, 4 values were not.

## Set-up fixed in advance
Flat unit ball, units 1/(8G): area side = -2 x (eps^2 coefficient of the area at fixed size).
Two data: static electric tide (exact static Weyl metric, U = eps r^l P_l, K = 0: pure metric data) and
stationary magnetic tide h = 2tK (h = 0 on the slice: pure extrinsic-curvature data).
Because both data are l >= 2 at first order, "area at fixed size Z" = A_2 - (dA/dR)/(dZ/dR) Z_2.

- Electric: brute force in (r, cos theta) from the exact Weyl metric
  e^{-2U}[e^{2k}(dr^2 + r^2 dtheta^2) + r^2 sin^2 dphi^2], k from the Weyl quadrature. Edge r = 1 + eps a P_l + eps^2 c.
  This does not use D-DEP's constraint identity, so agreement with D-DEP's F2(a) for every a is a real test.
- Magnetic: no first-order scalar on the slice (odd parity, K_rr = 0, tr K = 0), so no edge shape or time shift
  can be driven under ANY rotation-invariant rule. Only the l = 0 part of the second-order metric enters, and the
  Hamiltonian constraint fixes it: areal gauge, dr^2 (1 + 2m/r), 4m'/r^2 = <|K|^2>.
  fixed volume: -2 dA = 16 pi int_0^1 r m dr;  fixed geodesic radius: -2 dA = 16 pi int_0^1 (m/r) dr.

Switches:
 (a) maximal slice, volume, edge extremised.
 (b) edge not extremised, volume kept: geodesic ball a = 1/(l+1); also diamond edge a = 2/(l+1) and the
     isotropic-coordinate sphere a = 0 for the record.
 (c1) geodesic ball, bare area at fixed geodesic radius (size AND edge switched; Wang's "GCD bare area").
 (c2) size = mean geodesic radius over directions at the centre, edge extremised over a under that constraint
      (size switched alone). Geodesic distance includes the second-order bending of the geodesic, and the mean is
      over initial directions (Jacobian of the first-order deflection included).
 (d) slice: constant-mean-curvature slice through the same edge, first order in tr K. By hand before the run:
     the maximal slice maximises volume, so d Vol/d(tr K) = 0 at tr K = 0 in any spacetime; nothing moves at
     first order. This is an argument, not a script output. Second order is not attempted.

## Priors and what each outcome means
 (a) 1 and (l+2)/(l+1) (99%). If not: I am using the machinery wrongly; nothing else is claimed.
 (b) electric 2/9 (l=2), 5/32 (l=3) (95%, inherited); magnetic unchanged, 4/3, 5/4 (exact by parity).
     Meaning: switch (3) alone destroys "equals S on metric data" and leaves the K excess alone.
 (c1) 14/45 and 28/15 at l = 2 (Wang in print) (85%). A miss means my geodesic-radius bookkeeping is wrong.
 (c2) electric: not 1 (prior 90%), value unknown. Magnetic: (l+2)(2l+3)/((l+1)(2l+1)) = 28/15, 45/28.
     If electric c2 = 1: balance on metric data would be due to edge extremisation alone, not to volume.
     If electric c2 != 1: balance on metric data needs volume AND the extremal edge.
 Expected verdict (70%): metric balance is destroyed by either switch (2) or (3); the K excess is changed by
 switch (2) only; switch (1) does nothing at first order. So the asymmetry is a combination: the K excess is
 set by the volume rule alone, the metric balance by volume plus extremal edge.

## v^3 f
 By hand: for single-power cone data gamma ~ v^n, Q_2[v^3 f]/Q_3[v^2 f] = ((n+2)/(n+1))^2 (n+1)/(n+2)
 = (n+2)/(n+1); static tide n = l. So "4/3" is (l+2)/(l+1) for every l on the static tides (97%).
 "7/4 -> 1" on psi_dot = phi (odd metric rung): unknown at l = 3, 4; prior 50% that it stays 1.
 Identity outcome: static (l+2)/(l+1) AND rung = 1 at l = 2, 3, 4. Even then it is an identity on the two
 lowest profiles only, not of forms (see disclosure), and it can match one parity only (cone forms are
 duality blind; the fixed-volume even momentum rung is 8/7, not 1).
 Accident outcome: rung != 1 at l = 3 or 4.

## Plant
 `--plant k1`: "Area - 2 Volume" becomes "Area - 1 Volume" (coordinate size c = 0). Must move both entries of (a).
