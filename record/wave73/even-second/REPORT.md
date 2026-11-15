# REPORT: lane EVEN-SECOND (wave 73)

## Lead's check and grade

Sunday 15 November 2026. Blind second code (Fable): recompute, from definitions written in its brief, the seven constants and every derived quantity that EVEN-PASSAGE and EVEN-REMAINDER certified, without opening those folders or any report. The lane returned the text below; the lead saves it.

**Grade: 392 of 420.** Blindness was by construction (the brief held only definitions, no values, and fenced off the folders that hold them), the method differs from the first code (exact rational closed forms and series with explicit tails, against the first code's ball quadrature), each number was also computed a second way inside the lane, and the plant was caught. Marked down only because the lead cannot verify that the lane obeyed the read fence beyond its own statement and its file list.

**Comparison made by the lead (first code's certified bound, then this lane's enclosure, rounded):**

| Quantity | Needed as | First code | Second code | Consistent |
|---|---|---|---|---|
| A_inf | lower | 0.938463 / 0.938464 | 0.9384639805 | yes |
| A_16 | upper | 0.938807 | 0.9388067264 | yes |
| G_hi | upper | 0.770249 / 0.770250 | 0.7702493170 | yes |
| Q_16 | upper | 16.70472 | 16.7047178144 | yes |
| S' | upper | <= 15.12692 | 15.0220114366 | yes, first is looser |
| C | upper | <= 0.473687 | 0.4729137631 | yes, first is looser |
| s_3 | upper | <= 0.085899 | 0.0855673037 | yes, first is looser |
| c (first-order floor) | lower | 2.110e-3 | 2.1193e-3 | yes, first is lower |
| E_lo | lower | 0.14972 | 0.1497222621 | yes |
| N3_hi | upper | 0.32082 | 0.3201917593 | yes, first is looser |
| C_0 | upper | 15.509 | 15.4882294337 | yes, first is looser |

Every first-code bound lies on the safe side of the second code's enclosure: the three sups were bounded by cell covers in the first code and located exactly here, which accounts for all the differences. The published floor c = 2.110e-3 and remainder constant 15.509 stand; the sharper values the definitions allow are 2.1193e-3 and 15.488.

**Re-run by the lead:** `powershell -File scripts\run_all.ps1`, about 70 s, exit 0: "all sup checks agree: True", "all integrals agree: True; plant caught: True". The run appended `logs\004-run.log`, which is the lead's.

**Not checked by the lead:** the series tails and the monotonicity arguments for the sups were read in the report, not re-derived.

**What it changes.** "Certified constants not re-derived independently" is closed for the even family. The odd family's constants (ODD-SIGN: Q_16, S', C differ) have no second code yet.

---

## Lane report (unchanged)

All seven constants and every derived quantity are enclosed to about 40 digits, and the in-lane second route agrees on every one. The planted error is caught.

**Enclosures** (15 digits, lower column rounded down, upper column rounded up; true widths are 1e-40 or smaller)

| Quantity | Lower | Upper | Use as | Method |
|---|---|---|---|---|
| ‖u‖² | 2.39630696101284 | 2.39630696101285 | – | exact rational polynomial integral |
| A_inf | 0.938463980485440 | 0.938463980485441 | lower: 0.938463980485440 | exact rational closed form times interval e^4, e^-12 |
| A_16 | 0.938806726367441 | 0.938806726367442 | upper: 0.938806726367442 | series in q = e^(-2y)/1088, 40 terms, explicit positive tail |
| G_inf | 0.766842297451488 | 0.766842297451489 | – | closed form of ∫(u'')² − 2u''u e^(-2y) + u² e^(-4y), then sqrt |
| ‖(V_16−V_inf)u‖ | 0.00340701952158627 | 0.00340701952158628 | – | series for ((1−q)^(-2) − 1)², coefficients C(n+3,3) − 2(n+1), plus tail |
| G_hi | 0.770249316973074 | 0.770249316973075 | upper: 0.770249316973075 | sum of the two rows above |
| Q_16 | 16.7047178143806 | 16.7047178143807 | upper: 16.7047178143807 | weight expanded as e^y + 1/λ + e^(-y)/1088, same series and tail |
| S' | 15.0220114365668 | 15.0220114365669 | upper: 15.0220114365669 | log-derivative strictly decreasing, verified sign change at p = 1.18754392853915… |
| C | 0.472913763057772 | 0.472913763057773 | upper: 0.472913763057773 | exact maximiser, root of 6p² + 6p/π − 1 = 0, p = 0.279019636929870… |
| s_3 | 0.0855673036727704 | 0.0855673036727705 | upper: 0.0855673036727705 | cube of sup of p coth(πp) e^(-p²) = 0.440658970707057…, at p = 0.658813285776856… |
| J0_lo | 0.363567554126733 | 0.363567554126734 | lower | interval propagation |
| J1_hi | 0.230245129242697 | 0.230245129242698 | upper | interval propagation |
| s* | 1.77710190033248 | 1.77710190033249 | lower | interval propagation |
| M_lo | 0.00749966997515246 | 0.00749966997515247 | lower | interval propagation |
| Den_hi | 3.75536487977101 | 3.75536487977102 | upper | interval propagation |
| E_hi | 0.471155293799388 | 0.471155293799389 | upper | interval propagation |
| c | 0.00211931709539644 | 0.00211931709539645 | lower: 0.00211931709539644 | interval propagation |
| E_lo | 0.149722262105171 | 0.149722262105172 | lower: 0.149722262105171 | interval propagation |
| N3_hi | 0.320191759283526 | 0.320191759283527 | upper: 0.320191759283527 | interval propagation |
| C0 | 15.4882294336885 | 15.4882294336886 | upper: 15.4882294336886 | interval propagation |
| c/(2 C0) | 6.84170228905152e-5 | 6.84170228905153e-5 | lower | interval propagation |

For s_3 the log-derivative is not globally monotone, so the proof has three parts:
- **Window:** on J = [0.6, 0.75] the log-derivative is strictly decreasing with a verified sign change.
- **Outside the window:** 32 cells covering [0, 0.6] and [0.75, 3] each bound the function strictly below the value at the maximiser, using that x coth x is increasing.
- **Tail:** for p ≥ 3 the bound (p + 1/π) e^(-p²) is decreasing and already far below.

**Self-check**
- **Integrals:** Arb rigorous quadrature (`acb.integral`) straight from the definitions, with no series or closed forms, overlaps the primary value for all six integrals (‖u‖², A_inf, A_16, G_inf², the squared difference norm, Q_16). Quadrature radii are 1e-29 to 1e-35; midpoints differ by 1e-44 or less.
- **Sups:** an independent cell cover of the whole half-line, tolerance 1e-9 relative, contains the primary value for S', C and the s_3 base sup (547,907, 114,168 and 145,968 cells).
- **Plant:** with bracket exponent −1 instead of −2, the closed form gives A_16 = 0.938635114759952, which is 1.7161e-4 below the quadrature value with no overlap, so it is caught.

**Re-run** (about 70 s, writes a numbered log):
`powershell -File record\wave73\even-second\scripts\run_all.ps1`

**Ambiguities resolved**
- G_hi is G_inf plus the unsquared norm of (V_16 − V_inf)u; the derived formulas then square G_hi.
- "16" means λ = √272 exactly, so 4λ² = 1088 and 1/λ is an interval square root in Q_16.
- The e in J1_hi is Euler's number.
- Derived quantities are true enclosures of the stated formulas evaluated at the exact constants, not re-rounded at each step. If the other lane rounded the inputs outward first, its derived values will sit slightly on the safe side of mine.

**Weakest enclosure:** S' (width 4e-40) and its dependants Den_hi and C0, limited only by the 1e-42 bracket I chose for the numerically located maximiser. The weakest second-route check is the sup cell cover, at 9 digits by design.

Files are in record\wave73\even-second\:
- scripts\second.py
- scripts\check_flint.py
- scripts\run_all.ps1
- out\values.json
- logs\001-run.log to 003-run.log (003 is the clean run; 001 and 002 are the two runs before output-formatting fixes, with the same numbers)
