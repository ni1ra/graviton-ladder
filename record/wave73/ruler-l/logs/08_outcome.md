# 08 OUTCOME, lane RULER-L

Against 00_prereg.md: outcome A, all six expectations held.
- l = 2 (run 07, fixed script, with part 2): 40 checks, 0 failures; Delta[psi_0, psi_0] = 1/36, static 16/9.
- l = 3, 4, 5, 6 (runs 05, fixed script, four processes in parallel): 31 checks each, 0 failures.
  Times under that load: 307, 297, 804, 700 s. l = 5 and 6 exceed the 10 minute guide when run side by side;
  not timed alone.
- Part 2 (runs 01 and 07): identity exact; monomials ((m+4)/(m+3))^2; degree <= 8 range 1.011980867 to 5.724642;
  random search 1.056 to 3.217 (weak: it reaches neither end); family ratio 8.79 at N = 128.
- Plant (run 06, --plant swap, l = 2): 12 failures, including "Q_3 Gram equals c_l E_1"; static 16/9 -> 9/16.
- Void: runs 03 (see 04_note_nsimplify.md).
Not expected beforehand: the degree-d maximum grows faster than the proved ln(N)/270 family suggests
(5.72 at d = 8, 7.31 at d = 12); no growth law is claimed.
