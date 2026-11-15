# HINGE-BLIND addendum 2: continuation by a second lane (HINGE-BLIND-2)

Written 2026-11-15 by the lane that took over the interrupted run, while the n = 192 rows were
computing on this machine and before any n = 192 number, any outside-size row, any precision
rerun or any second-order (pt2) row had been looked at. PREREG.md and PREREG-ADDENDUM.md are
unchanged and are applied as written. The sealed folders are unopened; the HINGE lane's predicted
constant is unknown to this lane.

## What this lane knew when this was written

- PREREG.md, PREREG-ADDENDUM.md, the scripts and the inherited logs, read in full. That includes
  the interrupted lane's own closed form 7 zeta(3)/(16 pi^2) = 0.053285 stated in the first
  addendum. This lane re-derived that planar second-order integral by hand and by 2-d quadrature
  (planar_check.py) before the n = 192 rows existed. Every fit below is therefore made by someone
  who knows a candidate value; the rules are fixed here so that the choice of fit cannot follow it.
- Inherited rows n = 24..128 (primary, boundary-site), Srednicki rows to n = 96, three planted
  rows per n. One test run of addendum_rule.py restricted to n <= 128, in which only l <= 12 passed
  the 5e-5 bracket and three fits were kept (0.053240, 0.053287, 0.053335).

## Deviations and additions, with reasons

1. Environment. Rows are computed natively on Windows (Python 3.14, python-flint 0.9.0). One
   inherited row (n = 24, l = 8) was recomputed first and agreed to all 20 printed digits.
2. Second-order series. results/pt2 was not in the checkout and the inherited grid rows carry no
   Q2 field, so the Q2 series of the addendum rule is recomputed here (plans pt2_small, pt2_128;
   the new n = 192 rows carry Q2 themselves).
3. Lattice set for the addendum rule. Reported first: the rule on the PREREG lattice set
   n in {24, 32, 48, 64, 96, 128, 192}. If the extra lattices already named in run.py
   (n = 40, 56, 80, 112, 160, and l = 48, 64) are computed, the same rule on the enlarged set is
   reported beside it, not instead of it. Neither is chosen after seeing which is closer to
   anything.
4. Theory-free estimate, reported beside the addendum number. The addendum rule imposes a
   -(1/64) u ln(lam) term taken from an argument, in four of its five forms. A sceptic can reject
   that. Fixed here: the "free-log" estimate is the median and full range of the fits logfree-2,
   logfree-3 (log coefficient fitted) over the same windows (l >= 3, 4, 6), the same four series
   and the same 5e-5 inclusion cut. The pure-polynomial fits (poly-2, 3, 4) are listed as well.
5. Headline. The PREREG rule decides the pre-registered gates and is reported as written, pass or
   fail. The number handed to the lead for comparison is the addendum-rule median with its full
   range on the PREREG lattice set; the free-log range and the PREREG-rule value are given next to
   it, and if they disagree with the addendum range that disagreement is the first caveat.
6. Unfinished items are listed as unfinished, not extrapolated around.
