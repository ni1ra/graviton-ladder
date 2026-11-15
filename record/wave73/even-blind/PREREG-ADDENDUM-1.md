# EVEN-BLIND addendum 1: secondary analysis fixed mid-run (2026-11-15, about 19:45 local)

PREREG.md is unchanged and is applied and reported as written. This file fixes a second analysis
that will be reported beside it, and records what was known when it was fixed.

## What was known when this was written

- All rows with n <= 144 (130 rows), controls (i)-(v), precision and outside-size rows.
  n = 160, 176, 192 rows were computing; none had been looked at.
- Interim PREREG analyses on n <= 96 and n <= 144 (results\13_*, results\16_*). Seen in them:
  step-1 values Q2_l = 0.005325 (l=2), 0.007410 (4), 0.009312 (8), 0.010811 (16), 0.01148 (24),
  0.01185 (32), rising and decelerating; the three boundary conventions agree in the continuum at
  small l; the relative step-1 bracket exceeds the PREREG cut 1e-3 from l = 10 upward on n <= 144
  (2.3e-3 at l = 16, 4.4e-3 at l = 24, 8.8e-3 at l = 32), so no out-of-sample split was evaluable;
  the reference line "all constant forms", median 0.0125, range 0.0108 to 0.0136 (forms C1 0.0111,
  C2 0.0122, C3 0.0127, CL1 0.0133, CL2 0.0136).
- No closed form for the even case has been derived. The sealed folders are unopened and the
  claimed value is unknown. No fit below is chosen with reference to any target value.

## Why an addendum

The PREREG usability cut (bracket <= 1e-3 relative) was copied from the odd lane without a pilot.
Lattice artefacts here scale with l a and are several times larger than in the odd case (the lattice
value at n = 144 is 6 percent above the continuum at l = 16 and 15 percent at l = 32), so with
n <= 192 the cut will probably leave fewer than two test points in S1 and S2 and the PREREG gate
will then read "not evaluable / undetermined". That outcome will be reported as such.

## Secondary analysis (rules fixed here)

A1. Extra lattices, same rules. Rows with only the second-order coefficients (tag q2only: Q2 for
    ghost, trunc, site; no exact affinities) at n = 224 for l in {8, 10, 12, 16, 20, 24, 32, 40, 48},
    and at n = 256 for as many of l in {16, 20, 24, 32, 40, 48} as cost allows. They enter the q2
    series only, through the unchanged step-1 rule (eight largest lattices). The PREREG analysis is
    reported both without and with them.
A2. Usable cut relaxed to bracket <= 1e-2 x |central|. Step-1 rule otherwise unchanged.
A3. Fits weighted: every least-squares fit in step 2 uses weights 1/sigma_l,
    sigma_l = max(bracket_l, 1e-4 |Q_l|).
A4. Out-of-sample test: same forms, same splits S1, S2, S3, same pass criterion (every test l
    within 2 x bracket). Added, because wide brackets make passing easy: for every split the
    report gives, per test l, the spread of all forms' predictions in units of the bracket; two
    forms are "separated" by a split only if their predictions differ by more than 4 brackets at
    some test l. A family is called excluded only if all its forms fail both S1 and S2.
A5. Gate, secondary: same wording as PREREG gates (a)-(c) with A2-A4 in force, limit = median and
    full range over [constant forms passing S1 and S2] x [windows l >= 3, 4, 6] x [series q2, p1,
    mean]. If a zero-family or growth-family form also passes both splits, the verdict is
    "undetermined", and the report says which forms could not be separated.
A6. Descriptive, no gate: doubling differences D_l = Q_2l - Q_l and their ratios D_2l/D_l for
    l = 2, 4, 8, 16 (pure logarithmic growth gives ratio 1; an approach to a constant like (ln l)/l
    gives ratios falling below 1 towards 1/2); and the convention check (ghost, trunc, site
    continuum values of Q2 with brackets).

The headline of the report is whichever of the PREREG gate and this secondary gate is evaluable;
if both are, the PREREG gate is the headline and this one is given beside it. The number quoted
for the limit carries the full range over forms, not the best fit.
