# REPORT: lane EVEN-BLIND (wave 73)

## Lead's check and grade

Sunday 15 November 2026. The measuring lane was stopped mid-analysis and wrote no report. A separate read-back lane (Sonnet, read-only, confined to this folder) wrote the close-out below from the lane's own files. The lead then re-ran the analysis and spot-checked the close-out.

**Grade: 310 of 420 for the result as evidence** (below the 336 pass line as a blind confirmation; usable as supporting evidence with the caveats stated). Read-back lane: 365.

**What the result is, in one paragraph.** The pre-registered gate, applied as written, reads *undetermined*: with the pre-registered accuracy cut only l = 2..12 are usable and no out-of-sample split can be evaluated. A secondary analysis, fixed in `PREREG-ADDENDUM-1.md` while the run was in progress, gives the even constant as **0.013586, range 0.013555 to 0.013625** (forms CL2 and CL3, all eight n = 224 rows included). On the deciding series (q2, splits S1 and S2) every form that tends to zero or grows fails by 15 to 479 brackets; on the two boundary-convention series, which do not decide, the closest such form still fails, by 2.3 brackets. The value the lane was blind to is about 0.0136 and lies inside the range.

**Why it is not a full blind confirmation.**
- The rules of the secondary analysis were fixed after 130 of the 172 main rows (all n <= 144) had been seen, including an interim CL2 value of 0.0136. Only the n >= 160 and n = 224 rows are out of sample for those rules.
- The addendum's title says "about 19:45"; its file time is 19:14:32 and its contents match the 19:14 state (130 rows, numbers from file 16, no mention of file 19). The 31-minute conflict is unexplained. Timestamps are the only evidence of order; git has one commit for everything.
- CL2 and CL3 carry the log structure of the odd lane's addendum, which was read before `PREREG.md` was written.
- The outside-box effect (5.0e-6 at l = 2, n = 96) is larger than the step-1 bracket (3.5e-7) and is not in the quoted spread.
- n = 256 never ran, n = 224 at l = 10 is missing, control (iii) at n = 192 was not run.

**Lead's checks (run 15 Nov 2026, C:\Python314, mpmath 1.3.0, numpy 2.3.4).**
- `python eanalyze.py 192` reproduces `results\24_prereg_analysis.txt` line for line; the only difference is an appended exploratory block (addendum 2) that the on-disk script now prints.
- `python eanalyze.py 0 1e-2 1 1` reproduces `results\26_interim_secondary_partial224.txt` except at l = 40 and l = 48, whose n = 224 rows landed after file 26 was written. With them: CL2 and CL3 still pass S1 and S2; gate (c) limit median 0.013586, range [0.013555, 0.013625], half-range 0.3%; q2 series alone 0.013587 [0.013566, 0.013590]. File 26 itself says 0.013580 with the same range.
- Spot-checks of the close-out that held: 172 main rows; 130 of them dated at or before 19:14:32; eight n = 224 rows with l = 10 absent; file times of results 13, 16, 19, 24, 26 and of `eanalyze.py` (20:51:36, later than both result files); the quoted lines 16:146 and 19:322; the gate wording at `PREREG.md` 103-114.
- Not checked by the lead: the per-line citations into file 14 and the logs, and the outside-box and control (iii) remarks.

**How to cite it.** "A pre-registered lattice measurement was inconclusive under its own gate (too few usable sectors). A secondary analysis, fixed mid-run after about three quarters of the data, gives 0.01359 with a spread of 0.3% and rejects decay to zero and growth." Do not write "blind-reproduced" for the even constant.

---

## Lane report (read-back lane, unchanged)

**Even-blind close-out (read-back lane).** Paths are under `record\wave73\even-blind\`. Citations are `file:line`: `24` is `results\24_prereg_analysis.txt`, `26` is `results\26_interim_secondary_partial224.txt`, `19` is `results\19_interim_n160_secondary.txt`, `16` is `results\16_interim_n144_analysis.txt`, `14` is `results\14_controls.txt`. I ran nothing and wrote nothing.

**Value.** The measured quantity is Q2_l = l(l+1)·lim(1-A)/g² as g→0, on ball lattices with n=24..192 and N=4n, for l=2..48 (PREREG.md:24). The even constant is its l→∞ limit.
- **Secondary headline:** median 0.013580, full range [0.013555, 0.013625], half-range 0.3% (26:382). It comes from CL2 and CL3 only, over windows l≥3, 4, 6 and series q2, p1 and mean (26:383-387). The q2 series alone gives 0.013576 [0.013555, 0.013584] (26:383). Usable l is 2..48 (26:144), with n=224 present for l=8, 12, 16, 20, 24, 32 only (26:9-15).
- **PREREG run:** n≤192, usable l=2..12. CL2 gives 0.013572 [0.013558, 0.013633] (24:175). CL3 gives 0.013593 [0.013586, 0.013621] (24:176). All six constant forms pooled give 0.012910 [0.011259, 0.013633], labelled "spread, not gate" (24:166-167).
- **Drift as data were added:** 0.013551 on n≤144 (16:146), 0.013563 on n≤160 (19:322), 0.013580 on all rows (26:382).
- **Fitted log coefficient b:** CL2 -0.014774 and CL3 -0.015223 (26:401-402), against -1/64 = -0.015625 (ADDENDUM-2:17).
- **Not inside the quoted range:** CL1, which passes only S3 (26:181), has median 0.013407 (26:396). C3 has median 0.013031 (26:395).

**Gates as written (PREREG.md:103-114).**
- **Finite non-zero constant:**
  - (a) Undetermined. No split was evaluable, because there were fewer than two test points (24:145-147, 164).
  - (b) Vacuously met, since no zero-family or growth-family form passed (24:164).
  - (c) Undefined. It is defined over passing forms and none passed. The pooled reference, half-range 9.2% (24:167), is explicitly not the gate.
- **Tends to zero:** not met (24:164, zero = []).
- **Grows:** not met (24:164, grow = []).
- **Overall:** undetermined by PREREG.md:113-114. The word "undetermined" appears in ADDENDUM-2:9-10, not in file 24.
- **Boundary-convention agreement (PREREG.md:57-58):** the trunc and site continuum values differ from ghost by less than their brackets (24:116-129). This is step 1 only, because the trunc and site series have too few usable l for step 2 (24:159-161).
- **Controls:**
  - (i) Odd pipeline identical to 33 of 33 digits (`logs\03_odd_compare.log`).
  - (ii) g→0 converges (14:2-29).
  - (iv) Outside size is reported (14:60-70).
  - (v) Planted errors move Q by factors of 3.6 to 195 (14:51-58).
  - (iii) The local-momentum control was run, but n=192 was not (14:40-49).

**Secondary analysis (ADDENDUM-1, run in file 26).**
- **A2/A3:** relative cut 1e-2 and weighted fits were used (26:1).
- **A4/A5:** CL2 and CL3 pass both S1 and S2 (misses 0.0 to 0.3 brackets, 26:150-151, 166-167). Every zero-family and growth-family form fails by 12 to 479 brackets (26:152-158, 168-174). By the A5 wording this reads as the constant outcome, but file 26 does not print a verdict word.
- **Forms not separated:** C3|CL1, CL1|CL2, CL2|CL3 and ZP2|GL on S1 (26:160).
- **A6:** doubling-difference ratios fall from 0.912 to 0.706 for l=2, 4, 8, and from 0.834 to 0.661 for l=3, 6, 12 (26:132-140). That is descriptive only.
- **File 26 is explicitly partial224.**

**Timing of ADDENDUM-1.**
- Its file creation and last-write times are both 19:14:32, with no later edit. File 16 was written at 19:13:39.
- The addendum quotes numbers that match file 16 exactly: brackets 2.3e-3, 4.4e-3 and 8.8e-3 (16:12-15), and the CL2 value 0.0136 (16:146). It does not mention file 19.
- It says 130 rows were known. At 19:14:32 exactly 130 of the 172 main rows existed, all n≤144, and none with n≥160 (inferred from raw-file times). The first n=160 row is 19:18:46 (`logs\15_big160_early.log`), the first n=192 row 19:34:22 (`logs\10_big192.log`), and the first n=224 row 19:59:52 (`logs\18_q2only_224.log`).
- The addendum's own title says "about 19:45", which contradicts the file time by 31 minutes. By 19:45 there were 148 rows (inferred).
- **Reading:** the evidence points to writing at about 19:14, after n≤144 had been seen and after the PREREG cut was seen to leave nothing testable, but before any n≥160 or n=224 row. Its rules were therefore fixed with knowledge of 130 of the 172 rows. Timestamps are the only evidence, and the title conflict is unexplained.
- PREREG.md (18:46:14) predates the first even-perturbation number (18:48:34).

**Not done.**
- n=224 rows were planned for 9 values of l. Only 8 exist, and l=10 is missing. Rows l=40 and l=48 landed at 20:54 and 20:56, after file 26 was written (20:50:44), so they are unanalysed.
- n=256: no rows. The attempts in logs 22 and 23 wrote nothing, and `logs\00_pids.txt:3` records one stopped for memory.
- A1's "PREREG analysis reported with the extra rows" has no result file.
- ADDENDUM-2 items X1 to X3 have no output. The X code exists in `eanalyze.py` (created 20:51:36), but file 26 contains no X output.
- Control (iii) at n=192 was not run.

**Concerns.**
- ADDENDUM-2:8 says 186 rows. The lattice has 172 (inferred from the raw-file count).
- The outside-box effect is not in the quoted spread. For N=8n against 4n, l=2 shifts by 5.0e-6 at n=96 (14:67), against a step-1 bracket of 3.5e-7 (24:4).
- File 14 prints (1-A) in the (ii) block (lines 4-5, 8-9 and others) without exponents. The values only make sense as e-5.
- The control (iii) expectation of ln n growth (PREREG.md:129-131) is not what was seen. -lnA roughly doubles per doubling of n (14:47: 0.355, 0.691, 1.359).
- The on-disk `eanalyze.py` postdates results 24 and 26, so those outputs came from an earlier revision. I did not rerun it.
- The code uses windows l≥3, 4, 6 only, while PREREG.md:107 also lists "all usable l".
- ZP is fitted with s free and then classified by sign (eanalyze.py:181-182). A constrained s<0 power law is never tested as such.
- GL has 3 parameters against 4 to 5 for CL2 and CL3, so the exclusion is not like-for-like.
- The series q2, p1 and mean in the median are not independent.
- CL2 and CL3 carry the log structure of the odd lane's addendum, which was read before PREREG.md was written.
- File 26 existed 33 seconds before ADDENDUM-2. ADDENDUM-2:11 says step 1 only was looked at for the n=224 rows. This is inferred.
