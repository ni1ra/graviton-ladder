# 04 NOTE: run 03 (l = 3..6, first attempt) is void; cause found and fixed

Run 03 used RULER's `sp.nsimplify` to turn exact Gram entries into rationals. On rationals with long
numerators nsimplify returned a nearby algebraic guess (for example an entry containing 2**(27/55) at l = 4),
which made two gates fail at l = 3, 4: "Q_3 Gram equals c_l E_1" (l = 4 even) and the new identity gate
(l = 3 even, l = 4 even). This is an arithmetic-tool fault, not a physics result. At l = 2 the entries are
short and the fault does not occur (run 02 reproduces RULER exactly).
Fix: `exact()` = sp.cancel plus an assertion that the result is a rational number; no nsimplify anywhere.
The 03 logs are kept as they are. Runs 05 onward use the fixed script.

The void 03 run at l = 6 hung in the float eigenvalue step on a corrupted matrix and was stopped by this lane (its own process, found by its command line uler_l.py ../out/l6.json). The void 03 run at l = 5 finished with 0 failures (823 s); 03 at l = 3 had 1 failure, at l = 4 had 4.
