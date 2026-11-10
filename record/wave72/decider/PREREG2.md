# w72-decider, second registration: an attack on H4 through the Killing limit

Tuesday 10 November 2026, Oslo. Main thread (Fable). Written and pushed before `scripts/d2_killing.py` ran. The author, after the first result: "ok test them then!" (H4 and H0).

## What can be tested and what cannot

- **H4** (the area side is the stationary value of Area − 2·Volume) can be attacked where the answer is known. On a true Killing horizon the second-order first law must balance: the modular flow is a symmetry of the metric theory, and relative entropy equals canonical energy. The same machinery that gives the ladder on the ball must give no miss on the Rindler wedge. If it gives a miss there, the machinery is wrong and the ladder falls with it.
- **H0** (the free net is the G → 0 limit of the gravitational ball algebra) cannot be tested by a script. It needs a construction of the gravitational algebra at second order in which the vacuum modular generator is compared with the constraint (the record's open gate G1c). Nothing here touches it.

## The test

In the theorem's proof the two integers are two rates on a transverse, traceless, harmonic K:

- modular rate: the curvature velocity under the flow, 2 sym(∇N × curl K) with the sign conventions of THEOREM.md (K3 reads (l+1)K for the ball, N = (1 − r²)/2);
- geometric rate: Ric₃(2NK), the curvature of the lapse-weighted metric profile (K4 reads (l+2)K for the ball).

Claim to test, derived by hand before the run: for any lapse N and any TT harmonic K,

  Ric₃(2NK) − 2 sym(∇N × curl K) = 2 sym(K·∇∇N) − (ΔN)·K.

Consequences:
- Rindler wedge, N = x₁ (a Killing lapse, ∇∇N = 0): the two rates are equal identically. No miss.
- Ball, ∇∇N = −δ: the difference is exactly K, for every TT harmonic K of every degree, not only the scope family.

## Bins

| Bin | Outcome | Meaning |
|---|---|---|
| KILLING-OK | wedge difference ≡ 0 on all random data; ball difference ≡ K | the machinery balances where balance is known; the miss is exactly the failure of ζ to be Killing |
| DEFEATED | wedge difference ≠ 0 on any datum | the machinery contradicts the Killing-horizon first law; the ladder is not credible |
| IDENTITY-WRONG | wedge 0 but the general formula fails | my hand derivation is wrong; report what holds |

My odds before running: KILLING-OK 0.9 (I derived it by hand, so the script mostly checks my algebra), IDENTITY-WRONG 0.08, DEFEATED 0.02. Because the hand derivation exists, this is close to forced and carries little evidential weight. Its value is structural: it locates the miss in ∇∇N.

## Data and controls

- Random TT harmonic polynomial tensors of degree 1 to 4 from three families: sym∇(x × ∇φ) (magnetic scope type), Hess φ (electric type), and curl of the first.
- Lapses: x₁ (wedge), (1 − r²)/2 (ball), and a random quadratic polynomial (general identity).
- Each K is asserted TT and harmonic before use.
- Plant: change the flow coefficient 2 → 3 in the modular rate. The wedge equality must then fail.

## Also stated, no script (algebraic)

At this order the stationary value of Area − 2·Volume equals the area deficit at fixed volume of the isoperimetric surface, Jacobson's 2015 prescription: with the surface r = 1 + εaP_l + ε²b, fixing the volume sets b = −V₂(a)/4π, and then A₂ + 8πb = A₂(a) − 2V₂(a) = F₂(a) for every a. So H4 and the fixed-volume prescription give the same number.
