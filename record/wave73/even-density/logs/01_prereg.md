# 01 PREREG, lane EVEN-DENSITY (written before the first run, 2026-11-15)

Claim under test (proved on paper before the run): for h^2 = -d^2/drho^2 + l(l+1)/sinh^2(rho) on (0, inf),
EVERY pi in C_c^inf(0, inf) has <pi, h^-2 pi> finite. No orthogonality condition is needed.
Reason: the zero-energy solution regular at 0 (phi ~ rho^(l+1)/(2l+1)) grows linearly at infinity, so 0 is not a
resonance; the Green function is G = phi(rho_<) chi(rho_>) with chi = P_l(coth rho) -> 1, Wronskian 1, and
<pi, h^-2 pi> = <pi, G pi>. The solution w = G pi tends to the constant <phi, pi> (so w is NOT in L^2 unless
<phi, pi> = 0), but <pi, w> is finite all the same.

Predicted box law (Dirichlet at rho = L): I_L = I_inf - <phi, pi>^2 chi(L)/phi(L), so I_inf - I_L ~ <phi,pi>^2 / L.

Test functions, bump(a,b)(rho) = exp(-1/((rho-a)(b-rho))) on (a,b), else 0:
- pi1 = bump(1,3)                      (positive, <phi,pi1> != 0)
- pi2 = bump(1,3) - c bump(4,6), c chosen so that <phi, pi2> = 0   (orthogonal to the regular zero-energy solution)
- pi3 = bump(0.5,4.5) cos(4 rho)       (generic, sign changing)
l = 2 and l = 16. Boxes L = 8, 16, 32, 64, 128, 256, 512. Finite differences, step 0.004 (independent of G).
D(L) := I_2L - I_L (finite differences on both boxes).

Pass lines
- P1 (construction): Wronskian phi' chi - phi chi' = 1 within 1e-5 at rho = 1, 2, 4; for l = 1 the same code
  returns phi = rho coth rho - 1 within 1e-6 on [0.5, 6].
- P2 (finite differences against Green function): |I_L(FD) - (I_inf - <phi,pi>^2 chi(L)/phi(L))| <= 1e-3 * I_inf
  for every l, pi, L.
- P3 (convergence, pi1 and pi3): D(L)/D(2L) in [1.8, 2.3] for every L from 16 up, and
  D(L) / (<phi,pi>^2 (chi(L)/phi(L) - chi(2L)/phi(2L))) in [0.99, 1.01] for every L. Verdict CONVERGES.
- P4 (orthogonal case pi2): |D(L)| <= 1e-6 * I_inf for every L: the box value does not move at all. Verdict CONVERGES.
- P5 (proved bounds): 0 < I_inf <= ||pi sinh rho||^2 / lambda^2 and I_inf <= 4 ||rho pi||^2, every l and pi.
- P6 (w not in L^2): for pi1, w(L/2) on the box L = 512 is within 2% of <phi,pi1> (1 - (L/2)/L) scale, i.e.
  w_L(rho) = <phi,pi1> (chi(rho) - chi(L) phi(rho)/phi(L)); checked at rho = 256 to 1%. For pi2, |w(256)| <= 1e-6 max|w|.
- PLANT (must be caught): remove the potential and replace Dirichlet by Neumann at rho = 0 (a zero-energy resonance,
  the whole-line situation). Predicted: D(L) = L (int pi1)^2, so D(L)/D(2L) = 0.5 and the same convergence detector
  must return DIVERGES for pi1. Caught = verdict DIVERGES and D(L)/(L (int pi)^2) in [0.99, 1.01].

Convergence detector (fixed): CONVERGES if for all consecutive doublings from L = 16 up either D(L)/D(2L) >= 1.5
or |D| <= 1e-6 * |I at the largest box|; otherwise DIVERGES.

Amendment, still before the first run (no number seen yet): phi ~ rho - c_l at infinity with c_1 = 1 and c_l growing
like a harmonic number (about 3.4 at l = 16 if c_l = H_l, a guess), so the predicted D(16)/D(32) at l = 16 is about
0.5 * (64 - c)/(16 - c) = 2.4. The P3 window for D(L)/D(2L) is therefore [1.8, 2.6], not [1.8, 2.3].
P6 reads precisely: on the box L = 512, w_FD(256) = <phi,pi1> (chi(256) - chi(512) phi(256)/phi(512)) within 1%
(about <phi,pi1>/2: a non-decaying plateau), and for pi2 |w_FD(256)| <= 1e-6 max|w_FD|.

Exit code 0 only if every line holds and the plant is caught.
