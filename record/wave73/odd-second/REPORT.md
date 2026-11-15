# REPORT: lane ODD-SECOND (wave 73)

## Lead's check and grade

Sunday 15 November 2026. Blind second code (Fable): recompute, from definitions written in its brief, the constants and derived quantities that ODD-SIGN and ODD-REMAINDER certified for the odd family, without opening any other folder of the record. The lane returned the text below; the lead saves it.

**Grade: 390 of 420.** Blind by construction (definitions only, every other folder fenced off), a different method from the first code (exact rationals and closed forms against ball quadrature; sups located by concavity and exclusion covers), each number computed a second way inside the lane, two plants caught. Marked down only because the lead cannot verify the read fence beyond the lane's statement and file list, and the second-way check of the sups reaches 8 to 9 digits, not 30.

**Comparison made by the lead (first code's certified bound, then this lane's enclosure, rounded):**

| Quantity | Needed as | First code | Second code | Consistent |
|---|---|---|---|---|
| A_inf | lower | 0.938464 | 0.9384639805 | yes |
| A_16 | upper | 0.938807 | 0.9388067264 | yes |
| G_hi | upper | 0.770249 | 0.7702493170 | yes |
| Q_16 (odd weight) | upper | 26.90624 | 26.9062410592 | yes |
| S' (odd) | upper | <= 25.55923 | 25.3498796541 | yes, first is looser |
| C (odd) | upper | <= 0.115848 | 0.1153953872 | yes, first is looser |
| s_3 | upper | <= 0.085898 | 0.0855673037 | yes, first is looser |
| c (first-order floor) | lower | 2.038e-3 | 2.0505e-3 | yes, first is lower |
| E_lo | lower | 0.09370 | 0.0936966276 | yes (first rounded up in the 5th digit; its script value 0.093697 is what the assembly used) |
| N3_hi | upper | 0.20080 | 0.2004135228 | yes, first is looser |
| K | upper | 1.5715 | 1.5694809832 | yes, first is looser |
| C_0 | upper | 15.5105 | 15.4901564195 | yes, first is looser |
| l_1(1) | upper | 15237 | 15108.83 | yes, first is larger |

Every first-code bound lies on the safe side. The published odd floor c = 2.038e-3 and threshold l_1 = 15237 |g| stand; the sharper values the definitions allow are 2.0505e-3 and 15109 |g|. One note on E_lo: ODD-REMAINDER's report prints 0.09370, which is the lower bound 0.0936966 rounded to four figures in the unsafe direction; its script uses the unrounded value, so no inequality is affected, but the paper should print 0.09369.

**Re-run by the lead:** `scripts\odd_second.py` (under a second, exit 0) and `scripts\compare.py` (exit 0: "all sup checks agree: True", "planted sup mismatch caught: True"). The 3.5-minute cover script `sup_cover.py` was not re-run by the lead; its stored output was used by `compare.py`.

**What it changes.** "The odd family's constants have no second code" is closed. Both families' certified constants are now independently reproduced.

---

## Lane report (unchanged)

All 23 quantities are enclosed rigorously to about 30 digits, both self-checks agree, and both planted errors are caught. Every script exited 0 on the final run.

## Enclosures

Each value lies in [lo, hi], outward-rounded; the digits shown are common to both ends, with the last digit's lo/hi in brackets. The full 30-digit bounds are in `out\primary.json`.

| Quantity | Enclosure | Side used | Method |
|---|---|---|---|
| ‖u‖² | 2.39630696101284336578454225512[9,32] | exact | rational 262144/109395 |
| ‖u‖ | 1.5480009563992017656173786644[14,17] | exact | square root of the rational |
| ∫u′² | 0.72745032745032745032… | exact | rational 32768/45045 |
| A_inf | 0.9384639804854404918096878613[69,72] | LOWER | closed form (integration by parts, exact rationals times e⁴ and e⁻¹²) |
| A_16 | 0.9388067263674416664370385326[00,03] | UPPER | series in e^(-2y)/1088, 90 closed-form terms, tail ≤ 2.3e-116 |
| G_inf | 0.7668422974514880631530023079[96,99] | UPPER | square root of closed form |
| ‖(V_16−V_inf)u‖ | 0.0034070195215862755346129494[04,07] | UPPER | series for ((1−x)⁻²−1)², closed-form terms plus tail |
| G_hi | 0.7702493169730743386876152574[01,04] | UPPER | sum of the two above |
| Qodd_inf | 26.7604931153857749537124973959[90,93] | info | closed form |
| Qodd_16 | 26.9062410592414062903468757220[52,55] | UPPER | closed form, λ = √272 |
| Sodd | 25.3498796540859469957531711634[85,88] | UPPER | decreasing log-derivative, maximiser bracketed to 1e-60 at p ≈ 1.40345998181517648 |
| Codd | 0.1153953872029401706471180291[83,86] | UPPER | decreasing log-derivative, maximiser at p ≈ 0.66799062651859448 |
| s_3 | 0.0855673036727704287377864887[25,28] | UPPER | see below; maximiser at p ≈ 0.65881328577685656 |
| J0_lo | 0.3635675541267330432134447779[00,03] | LOWER | derived |
| J1_hi | 0.2302451292426970703645960108[88,91] | UPPER | derived |
| s* | 1.7771019003324851909595125929[47,50] | LOWER | derived |
| M_lo | 0.0046561639866403963015372193[86,89] | LOWER | derived |
| Den_hi | 4.8783847256915933449969946195[88,91] | UPPER | derived |
| E_hi | 0.2327380133303277474636953525[04,07] | UPPER | derived |
| c_odd | 0.0020504769976226179464733484[81,84] | LOWER | derived |
| E_lo | 0.0936966275916733522411926203[29,32] | LOWER | derived |
| N3_hi | 0.2004135228036301623841300291[44,47] | UPPER | derived |
| K | 1.5694809832271943224575656443[88,91] | UPPER | derived |
| C0 | 15.4901564194851709207575840392[61,64] | UPPER | derived |
| c_odd/(2 C0) | 0.0000661864522892522013362000[65,68] | LOWER | derived |
| l1 | 15108.8321765568735020229419151583[90,93] | UPPER | derived |

For s_3 the log-derivative of h = p coth(πp) e^(−p²) is not monotone on (0, ∞), so the proof has four parts:
- **Concavity:** (log h)″ < 0 on [0.5, 0.9], verified on 400 cells, so there is one critical point there, bracketed to 1e-60.
- **Near zero:** on [0, 0.05], coth x ≤ 1/x + x/3 gives h below the maximum.
- **Exclusion cover:** 8,700 cells on [0.05, 0.5] and [0.9, 3] each prove h below the maximum.
- **Tail:** for p ≥ 3, h ≤ 3 coth(3π) e⁻⁹.

## Self-checks

- **Integrals:** closed form against rigorous quadrature (`acb.integral` at 200 bits) overlaps for all seven: A_inf, A_16, G_inf², ‖(V_16−V_inf)u‖², Qodd_16, Qodd_inf and ‖u‖². The quadrature radii are 1e-24 to 1e-30.
- **Sups:** an independent branch-and-bound cell cover in mpmath interval arithmetic, with no derivatives, contains the primary enclosure for all three.

| Sup | Cover enclosure | Relative width |
|---|---|---|
| Sodd | [25.3498796540859, 25.3498797037872] | 2.0e-9 |
| Codd | [0.11539538720294, 0.11539538741831] | 1.9e-9 |
| s_3 | [0.085567303672770, 0.085567304133137] | 5.4e-9 |

## Planted errors

- **Qodd weight e^y replaced by e^(2y):** the planted closed form gives 630.6476… against quadrature 26.9062…, a 96% relative difference, flagged MISMATCH.
- **Wrong sup:** the primary Codd compared against the Sodd cover is flagged MISMATCH.

## Re-run

From `record\wave73\odd-second`:

```
& record\wave73\even-nonzero\venv\Scripts\python.exe scripts\odd_second.py
& C:\Python314\python.exe scripts\sup_cover.py
& C:\Python314\python.exe scripts\compare.py
```

The first takes under a second, the second about 3.5 minutes, the third is instant. Logs are `logs\01_primary.log`, `02_sup_cover.log` and `03_compare.log`; results are `out\primary.json` and `out\sup_cover.json`.

## Ambiguities resolved

- **‖u‖ in E_lo:** read as the L² norm, the square root of ∫u².
- **λ in Qodd_16:** λ = √272, the same λ² = l(l+1) = 272 used in V_16.
- **G norms:** taken over the whole line. u vanishes to order 4 at −2 and 6, so −u″ + Vu has no boundary singular part and the norm reduces to [−2, 6].
- **eps:** exactly 3/2.
- **s_3:** computed as (sup h)³, which equals the sup of the cube since h > 0.
- **Derived quantities:** each is a two-sided ball of the exact expression in the brief; the "side" column says which end to quote.

## Weakest enclosure

The primary enclosures are all far beyond 8 digits; the loosest are the sups and everything downstream of them (relative width about 3e-59, set by the 1e-60 bracket). The weakest evidence is the second-way cover for s_3 at 5.4e-9 relative. Its naive interval evaluation has first-order overestimation, so each extra digit costs about three times the cells (2.8 million evaluations already). It confirms 8 digits independently, not 30.
