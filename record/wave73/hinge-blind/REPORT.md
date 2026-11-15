# REPORT: lane HINGE-BLIND-2 (wave 73)

## Lead's comparison and grade (2026-11-15)

- The claim under test (`../hinge/REPORT.md`, lines 75 and 112): limit a1 = 7 zeta(3)/(16 pi^2) = 0.0532848, partial-sum slope 0.10657 g^2 per e-fold. The quantity HINGE asked to have measured is the one this lane measured: lambda^2 (1 - A_l) on the radial lattice.
- Measured: 0.05333, full range 0.05325 to 0.05349 (pre-registered lattice set); 0.05329, range 0.05325 to 0.05336 (enlarged set). The claimed value lies inside both ranges; the medians sit 0.09 and 0.015 percent above it.
- Not a clean pass: the pre-registered convergence gate failed as written (ratios 0.844, 0.738, 0.819 against a 0.6 threshold), so the pre-registered verdict is "undetermined". The headline comes from a rule fixed after data up to n = 128 had been seen.
- Not fully blind: PREREG-ADDENDUM.md, which the lead told the lane to read, contains the interrupted lane's own closed form 0.053285. The lane reported this itself.
- What this does and does not check: it confirms the numerical constant in the scaling limit of one scalar sector. It does not check HINGE's step from that constant to "no normal state or weight exists for free gravitons", nor the proof gap HINGE names.
- Grade: 340 of 420. Careful controls, honest about the failed gate and the weakened blind; marked down for no independent second code and a limit that is extrapolated from l up to 24.

The lane's report follows unchanged.

**HINGE-BLIND-2 report (2026-11-15)**

Tags: [M] measured this session, [I] inherited from the interrupted lane (re-run status stated), [R] reasoning only.

## 1. Verdict

The measured constant is **Q_inf = 0.05333, full range 0.05325 to 0.05349**, from the addendum rule on the pre-registered lattice set [M]. The same rule on the enlarged lattice set gives 0.05329, range 0.05325 to 0.05336 [M].

- **Pre-registered gate "tends to a finite nonzero constant": FAILED as written.** The doubling ratios are 0.844, 0.738 and 0.819, and the gate needs the last two at or below 0.6 [M]. No other branch is met either, so the PREREG verdict is "undetermined".
- **PREREG-form value:** 0.0509, with alternates 0.0532 and 0.0516, so roughly 0.0509 ± 0.0028. It excludes zero and its half-width is 5.5 percent, so only the ratio criterion fails [M].
- **Why the gate fails:** the approach to the limit behaves like (ln l)/l, for which that ratio sits near 0.75 by construction [R]. The gate was poorly designed for this approach; that is my judgement and I did not change the gate.
- **Sign and prior:** g = -1 lies above g = +1 at every l, as the prior said. The value is inside the prior's plausible range (0.03 to 3) but below its "order 0.1 to 1" [M].
- **Consequence for the sum:** logarithmic divergence with slope 2·c·g² ≈ 0.1066·g² per e-fold of the l cutoff [R from M].
- **Blind status:** not broken with respect to the HINGE folders. I never opened or read anything under `hinge\` or `hinge-e\`. I did list the `wave73` directory once, which showed folder names only.
- **Blind weakened in one way:** PREREG-ADDENDUM.md, which you told me to read in full, states the interrupted lane's own closed form 7ζ(3)/(16π²) = 0.053285. I knew that candidate before analysing n = 192. If it equals HINGE's prediction, this run was not blind to it.

## 2. How the number was obtained

**Quantity.** PREREG.md defines a free massless scalar, angular sector l, on the unit ball, with Q_l = l(l+1)(1 − A_l)/g². It is not a graviton computation; your brief said gravitons.

**Lattices.**
- [I, not re-run] n = 24, 32, 48, 64, 96, and n = 128 for l ≤ 24. I recomputed n = 24, l = 8 and it agreed to 20 digits [M]. Nine of these rows were recomputed at about 1.4 times the bits and agreed to all 30 stored digits [M].
- [M] n = 128 l = 32; n = 192 for l = 2, 3, 4, 6, 8, 12, 16, 24, 32; extra n = 40, 56, 80, 112, 160; l = 48 and 64 at n = 64 to 192.
- [M] The second-order series Q2 was recomputed, because `results\pt2` was not in the checkout.
- N_tot = 4n throughout; 193 rows completed with no failures and no precision-margin rejections; precision margin was 240 to 360 digits on the rows whose logs I read.

**Step 1, n → ∞ at fixed l.** Exact-fit Richardson in {1, x, x² ln x, x², …}. Brackets on the mean series [M]:

| l | PREREG lattice set | Enlarged set |
|---|---|---|
| 2 | 1.2e-6 | 5e-7 |
| 16 | 4.0e-5 | 1.6e-5 |
| 24 | 1.1e-4 (excluded) | 3.2e-5 |
| 32 | 5.6e-4 (excluded) | 6.0e-5 (excluded) |

The two sets agree within brackets at every l. Continuum values on the enlarged set, mean of g = ±1: 0.037699 (l=2), 0.042270 (4), 0.046090 (8), 0.048875 (16), 0.050043 (24), 0.050703 (32).

**Step 2, l → ∞.** Spread by choice of form and lattice set [M]:

| Analysis | PREREG lattice set | Enlarged set |
|---|---|---|
| Addendum rule (median, range) | 0.053331, 0.053251–0.053489 (20 fits) | 0.053293, 0.053254–0.053361 (40 fits) |
| Free-log, log coefficient fitted | 0.053408, 0.053262–0.053713 | 0.053305, 0.053244–0.053560 |
| Polynomial order 2 | 0.05196–0.05257 | 0.05218–0.05268 |
| Polynomial order 3 | 0.05255–0.05275 | 0.05267–0.05289 |
| Polynomial order 4 | 0.05283–0.05285 | 0.05288–0.05294 |

- Per-series medians on the enlarged set agree: Q2 0.053293, mean 0.053296, g = -1 0.053295, g = +1 0.053296.
- The fitted log coefficient on the enlarged set is -0.0151 to -0.0185 (logfree-3 fits); the imposed value is -1/64 = -0.0156.
- The polynomial fits sit 1 to 2 percent lower and climb with order.
- PREREG rule in n (three-parameter fits): effective exponent p = 1.2 to 1.6, brackets 6e-6 (l = 2) to 4.5e-4 (l = 32).

**Out-of-sample test** (written after the main fits, not pre-registered) [M]. Fits on 3 ≤ l ≤ 16, predicting l = 24 (bracket 3e-5) and l = 32 (bracket 6e-5):

| Form | Miss at l = 24 | Miss at l = 32 |
|---|---|---|
| Log forms (fixed or free) | 2e-7 to 1e-5; 2.4e-5 for one fit | 2e-6 to 2e-5; 4.5e-5 for one fit |
| Polynomial order 3 | about 5.5e-5 (outside bracket) | about 1.1e-4 (outside bracket) |
| Polynomial order 4 | 1.9e-5 | 5e-5 |

At l = 48 (fit to l ≤ 24, bracket 2e-4) the log forms miss by about 3e-5 and polynomial order 3 by 1.5e-4.

**Theory side, separate from the lattice number** [R, plus quadrature M]. I re-derived the planar second-order limit by hand from the code's second-order formula and Rindler modes. It reduces to (1/32)·∫tanh²(y)/y² dy = 7ζ(3)/(16π²) = 0.0532848. Brute-force 2-d quadrature gives 0.0532846. The lattice median on the enlarged set is 1.5e-4 relative above this.

## 3. Controls

- **Known-answer validation** [I, re-run M; identical apart from float roundoff at 1e-15]:
  - thermal chain: 2e-15;
  - squeezed vacuum: 2e-16;
  - brute-force Fock space, one and two modes and a small lattice: ≤ 3e-15;
  - formula I against formula II on the ball: agreement to 437 to 1071 digits;
  - Srednicki large-l entropy ratio: 0.993 to 0.998.
- **Precision** [M]: 9 cases at 1.4 times the bits, |ΔQ| = 0 at 30 digits.
- **Outside size** [M], change in Q from N = 4n to 8n: l = 2 +1.5e-5; l = 4 +8e-7; l = 8 +2e-7; l = 32 ≤ 1.4e-8. From 8n to 16n: 3e-7 at l = 2.
- **Radial discretisation** [I, n ≤ 96 only], Srednicki minus direct second difference at n = 96: -3.3e-7 (l = 2) to +2.7e-8 (l = 32), shrinking roughly as n⁻².
- **Boundary-site variant** [M], continuum values: l = 2 0.0376986 against 0.0376995 (bracket 1.1e-5); l = 4 0.0422673 against 0.0422702 (bracket 3.2e-5). It converges much more slowly and fails the 5e-5 cut for l ≥ 6.
- **Planted errors** [I for 7 rows, M for 2], ratio to the correct value at n = 96:
  - boundary weight 0: 0.826 (l = 2), 0.340 (l = 32);
  - boundary weight 2: 1.221 (l = 2), 2.375 (l = 32);
  - factor 1/2 dropped from Z: 5.80 (l = 2), 4.01 (l = 32).
- **Float64 pipeline** [I, re-run M] at n = 24, l = 8: 10 of 24 modes are unphysical. With the usual clip the error is +1.9 percent.
- **Scaling in g** [M]:
  - at n = 96 the mean at |g| = 1/16 differs from Q2 by 1.4e-5 (l = 2) and 1.5e-7 (l = 32);
  - continuum Q(mean ±1/4)/Q(mean ±1) = 0.905 (l = 2), 0.990 (l = 8), 0.999 (l = 32);
  - odd part at |g| = 1 = 0.0128 (l = 2), 0.0052 (l = 8), 0.0015 (l = 32). The departure from quadratic at |g| = 1, l = 2 is 10 percent and real.
- **Far-region density check** [I, re-run M] at n = 64, l = 16: lattice over local-density ratio 0.87 to 1.05. That is only rough support for the -1/64 coefficient.

## 4. What a sceptic should doubt, ranked

1. **Is this the HINGE constant at all?** The lattice measures a scalar sector with a φ² perturbation, normalised as l(l+1)(1 − A)/g² with Z = (g/2)∫φ². A different normalisation moves the number by factors like 4. I cannot check this blind.
2. **The limit is extrapolated, not reached.** The largest usable l is 24, where Q is still 6 percent below the limit. The answer depends on allowing a (ln l)/l term: with it 0.0532 to 0.0536, without it 0.0520 to 0.0529. The evidence for the log is the fitted coefficient and the out-of-sample test.
3. **The pre-registered gate failed.** The headline comes from a post hoc rule the interrupted lane fixed after seeing data up to n = 128 and after deriving its own closed form.
4. **Candidate value known.** I knew 0.053285 before analysing n = 192. I fixed my rules in PREREG-ADDENDUM-2.md before looking, but had already seen one test fit on n ≤ 128.
5. **Enlarged set is tighter but not the pre-registered one.** Addendum 2 names the PREREG-set number as the headline; the enlarged set is reported beside it.
6. **Most rows with n ≤ 128 are inherited.** Only the spot check and the nine precision reruns were recomputed.
7. **Outside size N = 4n** biases l = 2 by 1.5e-5. That is small against the range, and l = 2 is outside the fit windows.

**Unfinished:**
- Srednicki rows at n = 128 and 192.
- n = 256 rows.
- l = 64 has too few lattices to use.
- No independent second code. The same `hb.py` produced everything, checked only by the validation cases.

## 5. Reproduce

Python 3.14 with `pip install python-flint numpy scipy mpmath` (python-flint 0.9.0), run from `record\wave73\hinge-blind`:

```
python run.py one D2 128 4 32
python run.py plan grid_192 6
python run.py plan pt2_small 5
python run.py plan pt2_128 5
python run.py plan planted 4
python run.py plan precision 4
python run.py plan outside 4
python run.py plan extra_n 3
python run.py plan extra_160 4
python run.py plan hi_l 4
python addendum_rule.py mid D2 prereg     (headline)
python addendum_rule.py mid D2 all        (enlarged set)
python analyze.py                         (PREREG gates)
python controls.py
python predict_check.py 16
python planar_check.py 100
```

Wall time was about 75 minutes, with peak memory about 6 GB across workers.

**Files added (nothing existing modified, nothing committed):**
- `PREREG-ADDENDUM-2.md`
- `addendum_rule.py`, `predict_check.py`, `planar_check.py`
- `logs\b2_*.log`, `logs\b2_controls.txt`, `logs\b2_pids.txt`
- `results\b2_tables_prereg.txt`
- `results\b2_addendum_prereg_set.txt`, `results\b2_addendum_all_lattices.txt`, `results\b2_addendum_all_lattices_lmax64.txt`
- `results\b2_predict_check_lfit16.txt`, `results\b2_predict_check_lfit24.txt`
- `results\pt2\*.json` (69 files)
- `results\raw\*.json` (120 new, 195 total)

No REPORT.md was written. The virtualenv lives in the session scratchpad, outside the repo.
