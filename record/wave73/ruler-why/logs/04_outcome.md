# 04 Outcome against 00_prereg (lane RULER-WHY)

Runs: 01 (l = 2, 3; v^3 f at l = 2..4; 32 checks, 0 fail, 6 s), 02 (final: l = 2, 3, 4; v^3 f at l = 2..5 plus the
even K rung; 46 checks, 0 fail, 19 s), 03 (plant k1, caught at l = 2, 3, 4: (a) becomes 101/6 and 2/3 at l = 2).
Between 01 and 02 only the l ranges and one extra datum were added; no formula changed.

Area side / S (electric = metric data, magnetic = K data):

| setting | l=2 | l=3 | l=4 |
|---|---|---|---|
| (a) maximal slice, volume, extremal edge | 1, 4/3 | 1, 5/4 | 1, 6/5 |
| (b) volume, geodesic ball | 2/9, 4/3 | 5/32, 5/4 | 3/25, 6/5 |
| (b') volume, diamond edge | 2/9, 4/3 | 5/8, 5/4 | 39/50, 6/5 |
| (c1) geodesic ball at fixed geodesic radius | 14/45, 28/15 | 45/224, 45/28 | 11/75, 22/15 |
| (c2) fixed mean geodesic radius, extremal edge | 28/15, 28/15 | 15/8, 45/28 | 53/30, 22/15 |
| (d) CMC slice, first order in tr K | no change (argument: volume is stationary at the maximal slice); not scripted |

Landed: the expected verdict. Metric balance is destroyed by switch (3) alone and by switch (2) alone; the K
excess moves only under switch (2). (c1) reproduces Wang's printed 14/45 and 28/15, which tests the
geodesic-radius bookkeeping (without the bending term the electric value would be 98/45).
Not expected: at l = 2 the (c2) electric and magnetic values coincide (28/15); they differ at l = 3, 4.

v^3 f: static tides give (l+2)/(l+1) at l = 2..5 (and for symbolic n by the profile law); the psi_dot = phi rung
gives exactly 1 at l = 2..5 in both parities. The fixed-volume values there are 1 (odd metric) and 1 + 1/(2l+3)
(even K), so the match is with the odd tower only, and at l = 2 it holds on the two lowest profile entries only
([psi_0, psi_0] and [psi_dot_0, psi_dot_0]); the rest of the form differs (row 0 of the psi block is off by
1/270, 1/165, 1/130). Why the rung gives exactly 1 for every l tested is not explained.

Not done: (d) beyond first order; size = proper time between tips (Wang's -56/9, 0); time-dependent data.
