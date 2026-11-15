# EVEN-BLIND addendum 2: two exploratory items, fixed before they were run (2026-11-15, about 20:55)

PREREG.md and PREREG-ADDENDUM-1.md are unchanged. Nothing here enters any gate or the headline
range; both items are reported under the label "exploratory, post hoc".

## What was known when this was written

- The complete PREREG lattice set (n = 24..192, 186 rows) and its PREREG analysis
  (results\24_prereg_analysis.txt): usable l (cut 1e-3) = 2..12, so S1, S2, S3 are NOT EVALUABLE and
  the PREREG gate reads "undetermined". Reference fits of the constant family on l <= 12: C1 0.0117,
  C2 0.0125, C3 0.0129, CL1 0.0134, CL2 0.01357, CL3 0.01359; fitted log coefficient b = -0.0117
  (CL1), -0.0147 (CL2), -0.0153 (CL3).
- An interim run of the addendum-1 secondary analysis on n <= 160 (results\19_*): CL2 and CL3 passed
  S1 and S2 (misses 0.1 to 0.2 brackets), every zero-family and growth-family form failed by 7 to
  240 brackets, C3 and CL1 failed or nearly failed.
- Six of nine n = 224 rows (step 1 only looked at).
- The odd lane's far-region term -(ln l)/(64 l) (hinge-blind addendum; -1/64 = -0.015625). I have
  no derivation that the even perturbation has the same coefficient; the free fits above happen to
  land near it. No closed form for the even limit has been derived or attempted, and the claimed
  value is still unknown to this lane.

## Exploratory items

X1. Fixed-log forms XF2 = q + c1 u + c2 u^2 - (1/64) u ln lam and XF3 = XF2 + c3 u^3, run through
    the same splits and windows as the other forms, listed separately (family "explore"), never
    counted among the constant-family forms of any gate, never folded into the headline range.
X2. The limit of the series p1L = l^2 (1 - A_l) at g = 1 (the normalisation a previous lane quoted)
    and of p1 = l(l+1)(1 - A_l) at g = 1, fitted with the forms that pass the secondary gate,
    windows l >= 3, 4, 6, listed per series.
X3. Secondary analysis repeated unweighted (A3 off) and with the q2-only rows left out, to show how
    much the weighting and the extra lattices move the answer.
