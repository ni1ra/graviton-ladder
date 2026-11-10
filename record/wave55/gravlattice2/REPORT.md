# REPORT: w55-gravlattice2. Graviton lattice rerun, rule B fixed in advance, sector solver to R = 10, degree-5 controls: INCONCLUSIVE by the pre-registered power criterion (Paper A, Theorem 1)

Lane w55-gravlattice2 (Opus), Mon 9 Nov 2026, 09:20–10:05 Europe/Oslo. Claude only. The lead books credence; this report only proposes.

## Verdict: INCONCLUSIVE under PREREG §6 as written (Rmax = 10)

| §6 condition | Required | Result |
|---|---|---|
| (i) Controls converge | \|a_ctrl − 1\| ≤ 0.10 and \|ratio(Rmax) − 1\| ≤ 0.15 | **Met.** d1_E2: a_ctrl 1.078, ratio(10) 0.954. d2_B2: 0.958, 1.035 |
| (ii) All \|c_∞\| < 0.15 | on all six data | **Met.** −0.123 to +0.056 |
| (iii) Threshold plants c = ±0.15 outside the bracket | on ≥ 4 of 6 data | **Missed: 2 of 6** (g1, g6) |
| KILL | \|c_∞\| ≥ 0.15 with the bracket excluding 0 | **Not met.** No datum reaches it; every bracket contains 0 |

The test lacks power; the prediction did not fail. On g2–g5 the brackets are too wide to exclude |c| = 0.15. Their wide ends come from the quadratic fits and the raw convention. No negative is booked, so the false-kill rule is not triggered.

## 1. Pre-registration and freezing
- `PREREG.md` was pushed at 09:30 as `c14a5a975`, before any new row (sha256 in `out/prereg_time.sha256`). Before it, only the theory-side `out/controls.json` and `out/pred.json` existed.
- Fixed in advance:
  - rule B as the only verdict rule;
  - R = 3–10, with 11–12 optional;
  - the degree-5 controls;
  - the estimator;
  - the bins;
  - the plants, including two new ones (threshold c = ±0.15, no classical centre).
- `glat2.py` and `run2.py` were frozen at 09:31, before the first test row (`out/frozen_code.sha256`).
- `fit2.py` was written at 09:34. Its only later edit, at 09:36, changed number printing after a partial-data dry run; the estimator, bins and plants are untouched.

## 2. The solver: cubic symmetry, not the chain
- **Unchanged physics.** `glat.py`, `maxlat_core.py` and `pred.py` are byte-identical copies of the wave52 files.
- **New: Z2³ parity sectors (`glat2.py`).**
  - The three reflections about the centre (a subgroup of O_h) act on e = inc(h) and b = sym curl π as signed permutations. The signs are read off the stencils.
  - X, P and W are projected onto the eight characters.
  - Only the rows X[orbit reps, :] are built, and W stays sparse.
  - glat's Gaussian construction then runs sector by sector, with global null and rank thresholds.
  - Run-time checks: W off-sector ≤ 3e-16, data norm preserved to ≤ 2e-15, and dual-kernel sign consistency asserted.
- **Why it is not the chain.** Every sector mixes all l and both chain parities, and nothing is projected onto spherical harmonics.
- **Cost.**
  - At R = 7 a row takes 17 s and 0.4 GB, against 143 s and 3.0 GB for the dense code.
  - Peak memory per row: R = 8 0.78 GB, R = 9 1.53 GB, R = 10 2.57 GB (7.2 min, two rows in parallel).
- **R = 11.** The site row was stopped unfinished on the lead's instruction, so no R = 11 number exists. R = 12 would need about 7.7 GB per row.

## 3. Verification gate (PREREG §3)
- **Against the ten shipped wave52 rule-B rows (R = 3–7, s and c):**
  - g1–g6 match to ≤ 8.6e-12 at every cap: **pass.**
  - The old degree-11 controls c1 and c2 miss 1e-8: 6e-8 at cap 25, 2e-6 at cap 30 and 3e-5 at cap 35 (`out/verify.*`).
- **The deviation is the old code's own conditioning floor.** The unchanged wave52 code, run again single-threaded at 09:34 (before any R ≥ 8 row), misses its own shipped c1 and c2 values by 1e-8 to 4e-6 (`out/floor_check.log`, `out/floor/`). The cause is modes with 2ν − 1 ≈ 1e-11 to 1e-15, whose capped logarithm is set by rounding.
- **No primary-control row is affected.** The dense wave52 path, run on the new controls at R = 5 and 6 for both centrings (`out/dense*`), agrees as follows:
  - d1_E2 and d2_B2: ≤ 2.7e-9 at every cap;
  - d3_E3: ≤ 3.2e-9 at cap 30, with one entry at 1.4e-8 at cap 35;
  - g1–g6: ≤ 1.2e-11.
- c1 and c2 enter only the alternate estimates for g2 and g6.

## 4. Results (`out/fit.log`, `out/fit.json`)

### 4a. Ratio S_lat/S_T1
Rule B, mean of site and cube, vol convention, cap 30.

| Datum | R = 3 → 10 | lin[6,10] [bracket] |
|---|---|---|
| g1 uniform tidal E | 0.796, 0.888, 0.939, 0.944, 0.976, 0.972, 0.962, 0.983 | 1.019 [0.952, 1.030] |
| g2 uniform B | 1.082, 0.994, 1.008, 0.995, 1.038, 1.012, 0.990, 1.009 | 0.999 [0.891, 1.021] |
| g3 l = 3 | 0.686, 0.815, 0.901, 0.929, 0.960, 0.958, 0.941, 0.973 | 1.008 [0.886, 1.042] |
| g4 l = 2 + 3 | 0.708, 0.829, 0.908, 0.932, 0.963, 0.961, 0.945, 0.975 | 1.010 [0.899, 1.040] |
| g5 l = 4 | 0.707, 0.807, 0.827, 0.929, 0.936, 0.958, 0.932, 0.962 | 0.991 [0.661, 1.102] |
| g6 E + B | 0.935, 0.939, 0.973, 0.969, 1.006, 0.992, 0.976, 0.995 | 1.010 [0.937, 1.026] |
| d1_E2 control (deg 5) | 0.483, 0.681, 0.817, 0.859, 0.915, 0.915, 0.930, 0.954 | 1.078 [0.939, 1.086] |
| d2_B2 control (deg 5) | 1.315, 1.139, 1.074, 1.059, 1.090, 1.061, 1.011, 1.035 | 0.958 [0.721, 1.047] |
| d3_E3 control (deg 6) | 0.512, 0.679, 0.797, 0.821, 0.907, 0.896, 0.908, 0.944 | 1.093 [0.939, 1.093] |
| c1, c2 old (deg 11) | 0.338 … 0.900; 0.946 … 0.912 | 1.115; 1.105 |

### 4b. Surface coefficient
c(R) = (ratio_datum − ratio_ctrl)/(σ/S_T1). c_∞ is the lin fit on [6, 10], with the 12-estimate bracket.

| Datum (control) | σ/S | c(10) | c_∞ | Bracket | ±0.15 both outside? |
|---|---|---|---|---|---|
| g1 (d1) | 0.743 | +0.039 | −0.079 | [−0.105, +0.039] | yes |
| g2 (d2) | 0.743 | −0.035 | +0.056 | [−0.092, +0.229] | no |
| g3 (d1) | 0.716 | +0.026 | −0.097 | [−0.212, +0.026] | no |
| g4 (d1) | 0.721 | +0.029 | −0.094 | [−0.191, +0.029] | no |
| g5 (d1) | 0.700 | +0.012 | −0.123 | [−0.568, +0.127] | no |
| g6 (d1/d2 mix) | 0.743 | +0.003 | −0.013 | [−0.036, +0.064] | yes |

- **Site-only, cube-only and alternate-control estimates:** every c_∞ lies in [−0.146, +0.074], and every bracket contains 0.
- **Post hoc, does not move the bin.**
  - c(R) falls monotonically towards 0 on every datum (g1 0.42 → 0.04, g2 −0.31 → −0.04).
  - The negative c_∞ values are the lin fit overshooting past 0.
  - Unextrapolated, c(10) lies in [−0.035, +0.039].

## 5. Plants
A plant is caught when its value lies outside the primary c bracket.

| Plant | Result |
|---|---|
| 1. Surface (π/2)σ | Caught on all six |
| 2. c_l → 1/L | Caught on all six. In the c_l differential test (D = ratio/ratio(g1)) it is caught on g3, g4 and g5: plant 2.50, 2.21 and 4.50 against D = 0.989, 0.991 and 0.972 |
| 3. c_l → 4/L² | Caught on all six. In the differential test it is caught on g3, g4 and g5 (better than the pre-registered "low power") |
| 4. Born coefficient 0 | Caught on five of six (missed on g5). In the differential test it is caught on g3, g4 and g5 |
| 5. **New:** threshold c = ±0.15 | Both signs caught only on g1 and g6. This is the power criterion that sets INCONCLUSIVE |
| 6. **New:** drop the classical centre | Caught on g2 (c_∞ −0.62 [−0.86, −0.21]) and g6 (−0.34 [−0.40, −0.15]). No power on E-only data, as pre-registered. The centre holds 62% of S at R = 10 for uniform B |
| 7. Injection of c = 0.15 and π/2 | Recovered to ≤ 2e-6 (a linearity check) |

## 6. Deviations from PREREG
- **Verification gate booked as passed, not BLOCKED.** The evidence is in §3; the 1e-8 criterion on c1 and c2 is below the reference code's own reproducibility. This was decided before any R ≥ 8 row existed.
- **R = 11 started one centring at a time.** A parallel pair projected to about 9 GB. It was stopped unfinished on the lead's instruction, so the verdict uses PREREG's primary set, Rmax = 10.
- **Nothing else changed.** All 16 rows are shipped.

## 7. Proposed likelihood ratio
- **Narrow clause.** The pre-registered bin is INCONCLUSIVE, so LR ≈ 1 on "no surface term with |c| ≥ 0.15" (PREREG §8).
- **Paper value.** LR ≈ 1.15 on Theorem 1(c)'s "no natural-size surface term, correct c_l", slightly above wave52's 1.1.
  - What lifts it: converged controls, a rule fixed in advance, the 4/L² and no-centre plants now caught, and every datum within 4% of S_T1 at R = 10.
  - It stays below the PASS value (1.3–1.5) because the threshold-power criterion failed.
- **C and M:** unchanged (≤ 1.01).
- **Both legs.** A surface term of |c| ≥ 0.5, or a wrong c_l, would have given a KILL (LR ≈ 0.1); neither appeared. 0.15 ≤ |c| ≤ 0.25 stays formally open on g2–g5.
- **Step most likely to need correction.** The extrapolation form: the lin fit overshoots while c(R) falls monotonically. Reaching R ≈ 14 needs a lighter X-row build, since R = 12 already projects to about 8 GB.
- **For the paper.** It may say: "An independent 3D lattice test, with the region rule fixed in advance and converged controls, finds |c_∞| < 0.15 on all six data. It excludes a natural-size surface term and confirms the l-dependence of c_l. A 0.15 resolution is not reached on four of the six data."

## Files
- `PREREG.md`, `REPORT.md`, `HASHES.md`.
- `scripts/`:
  - wave52 copies: `glat.py`, `maxlat_core.py`, `pred.py`;
  - theory side: `ctrl_solve.py`, `pred2.py`;
  - solver: `glat2.py`, `run2.py`, `run_set.sh`;
  - gate: `verify.py`, `floor_check.py`, `run_dense_ctrl.py`, `dense_compare.py`;
  - estimator: `fit2.py`;
  - `hashes.sh`.
- `out/rows/row_R{3..10}_{s,c}_B.json`, with logs that carry the `/usr/bin/time -v` time and memory.
- `out/verify.*`, `out/floor*`, `out/dense*`, `out/fit.*`, `out/pred.*`, `out/controls.json`, `out/ctrl_solve.log`, `out/frozen_code.sha256`, `out/prereg_time.sha256`.

To rerun, from `scripts/`:
1. `python3 ctrl_solve.py ../out/controls.json`
2. `python3 pred2.py ../out/controls.json ../out/pred.json`
3. `./run_set.sh s 3 4 5 6 7 8 9 10`
4. `./run_set.sh c 3 4 5 6 7 8 9 10`
5. `python3 verify.py ../out/rows ../../../wave52/gravlattice/out/rows ../out/verify.json`
6. `python3 fit2.py ../out/rows ../out/pred.json ../out/fit.json`

`$GCACHE` (a `.npy` of `G_array(36)`) is optional. Expect about 20 minutes on 8 cores and 2.6 GB peak per row.
