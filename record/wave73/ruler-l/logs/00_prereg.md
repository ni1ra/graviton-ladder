# 00 PREREG, lane RULER-L (wave 73). Written before the first run of ruler_l.py.

15 Nov 2026. Units and forms as in RULER: R = 1, Q_p = 2 oint int_0^1 v^2 f |sigma_p|^2 dv, f = v(1 - v),
sigma_p = -v^(-p) int_0^v u^p C du, S = Q_3, Delta = Q_2 - Q_3. nc_core.py is an unedited copy.

## Hand derivation done before any run (part 2)
Put h(v) = v^2 (sigma_2 - sigma_3). Then h(0) = 0, h' = v sigma_3, sigma_2 = h'/v + h/v^2, so
Q_3 ~ int v(1-v) h'^2, Q_2 ~ int v(1-v) (h' + h/v)^2, and after one integration by parts
    Delta = 2 oint int_0^1 h^2 / v dv = 2 oint int_0^1 v^3 |sigma_2 - sigma_3|^2 dv   (R drops out).
So I expect, for every l and every datum:
1. Delta > 0 unless C = 0. Every Gram block of Delta is positive definite; the full 8 by 8 form too.
2. inf Q_2/Q_3 = 1, not attained (monomials C = v^m give ((m+4)/(m+3))^2 -> 1).
3. sup Q_2/Q_3 = +infinity, even on polynomial C: family C_N with h' = v^2 (1-v)^N should give
   Q_2/Q_3 - 1 >= ln(N)/270. No upper bound; growth only logarithmic in the degree.
4. Per l: gates pass as at l = 2; even and odd Grams identical; static ratio ((l+2)/(l+1))^2 =
   25/16, 36/25, 49/36, 64/49 at l = 3, 4, 5, 6; cross block non-zero.
5. New exact gate: Delta Gram from Q_2 - Q_3 equals the Gram of 2 oint int v^3 (sigma_2 - sigma_3)^2.
6. l = 2 reproduction first: Delta[psi_0, psi_0] = 1/36, static 16/9.
Prior that 1 to 6 all hold: 90%. Least sure: run time at l = 6 (may stop at 5), and whether the chain
profile of the data stays inside j = 0..3 for l > 2 (the script asserts it).

## What each outcome would mean
A. All as expected: "the cone rule over-counts on every non-trivial datum, at every l, in both halves" is a
   theorem given the cone forms; the Gram tables are illustrations of it.
B. Gate 5 fails (Delta Gram differs from the h-form): my identity is wrong; part 2 is withdrawn, part 1 stands.
C. A Delta block not positive definite while gate 5 passes: impossible unless the data are linearly dependent
   on the cone (a zero mode); check the data count first.
D. Q_3 Gram differs from c_l E_1 at some l: the entropy identification fails at that l; stop and report.
E. Even and odd Grams differ: duality fails at that l; bug first.
F. Random search finds a ratio below 1 or above the degree-8 pencil maximum: the proof or the pencil is wrong.

## Plant
--plant swap: p = 2 and p = 3 exchanged (entropy form built from sigma_2, charge from sigma_3).
Expected: gate "Q_3 Gram equals c_l E_1" fails, static ratio becomes ((l+1)/(l+2))^2, Delta negative definite.
