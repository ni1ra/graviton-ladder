# REPORT: lane RULER (wave 73)

## Lead's check and grade

Sunday 15 November 2026. Small exact-computation lane (Fable): does the light-cone size rule (the one the printed first-law proof uses) share the metric-versus-extrinsic-curvature asymmetry of the fixed-volume rule, at l = 2? The lane returned the text below; the lead saves it.

**Grade: 380 of 420.** One question, an exact answer, the expectation and its reason written down before the first run (`logs\00_prereg.md`: duality swaps metric data of one parity with extrinsic-curvature data of the other, and the cone forms are duality blind, so the cone rule cannot carry an asymmetry common to both parities), gates passed before any result was read, and a plant that moves the static value. Marked down because positivity is shown on four polynomial profiles per block rather than proved, and because the split into "metric" and "extrinsic-curvature" halves is a t = 0 split the cone does not respect (the cross block is non-zero), which the lane states itself.

**Checked by the lead.**
- The pre-registration file predates the run log and states the expectation and the 90% prior.
- Consistency by hand: Delta[psi_0, psi_0] = 1/36 with the static ratio 16/9 means S[psi_0, psi_0] = (1/36)/(7/9) = 1/28; Delta[psi_dot_0, psi_dot_0] = 1/1008 with ratio 7/4 means S = (1/1008)/(3/4) = 1/756. Both are what the lane's gate "Q_3 Gram equals c_l E_1" must produce for it to pass; the lead did not recompute c_l E_1 independently.
- The logical core (duality forbids a parity-independent asymmetry on the cone) is a two-line argument and holds given NULL-CONE's duality statement, which was proved and checked there.

**Re-run by the lead:** `C:\Python314\python.exe ruler.py <out.json>` from `scripts\`, 76 s, exit 0, 68 checks, 0 failures; static electric and stationary magnetic both 16/9; duality gate passes; the full 8 by 8 form is positive definite in both parities (smallest eigenvalue 9.6e-10, exact leading minors positive). The plant run was not repeated by the lead.

**What it changes.** The ruler question is answered at l = 2: the cone rule balances on neither half of the data and treats both halves alike, so the metric-versus-extrinsic-curvature asymmetry is a property of the fixed-volume, maximal-slice rule, not of the gravitons. This sharpens the paper's statement of where the ladder comes from. Open: l = 3 to 6 (expected to follow from the same duality argument), and which single ingredient of the fixed-volume rule removes the excess on metric data. The lane's remark that the planted weight reproduces the fixed-volume odd numbers (4/3 and 1) is tagged unverified by it and must not be quoted as a result.

---

## Lane report (unchanged)

**RULER, l = 2, both parities**

**(1) Answer.** The cone ruler has no metric-versus-extrinsic-curvature asymmetry: it over-counts on both halves of the Cauchy data in both parities, so the asymmetry belongs to the fixed-volume ruler.

- Metric data: no balance; Δ = Q₂ − S is positive definite, both parities. [C]
- Extrinsic-curvature data: no balance; Δ is positive definite, both parities. [C]
- Cross block: not zero; Δ is still positive definite on the full 8-dimensional space. [C]

**(2) Evidence.**
- **Data.** TT gauge, evolved with `nc.tt_series`, cone forms from an unedited copy of `nc_core.py`. Even: the unique TT tensor with x·T·x = r^(2k)φ, k = 1..4, used as h and as K. Odd: O₀ = sym∇(x×∇φ) and three rungs above it, O₁..O₄ as h and O₀..O₃ as K. [C: TT and vacuum equations exact on all 16]
- **Basis.** Chain field ψ = r^(2j)φ and ψ̇ = r^(2j)φ, j = 0..3, divided by ∮φ². For even parity ψ is metric data and ψ̇ is K data; for odd parity it is the reverse. [C]
- **Gates.** The Q₃ Gram equals c_l·E₁ exactly and is block diagonal. [C] The even and odd Gram matrices are identical, as duality requires. [C; the duality argument itself is I, NULL-CONE]
- **Δ on the ψ block** (even metric, odd K), upper triangle:
  `[1/36, 1/36, 1/36, 1/36; 4/135, 61/1980, 37/1170; 657/20020, 187/5460; 7552/208845]`
- **Δ on the ψ̇ block** (even K, odd metric):
  `[1/1008, 1/1296, 1/1584, 1/1872; 131/213840, 157/308880, 61/140400; 307/720720, 227/618800; 504337/1587222000]`
- **Cross block** (rows ψ, columns ψ̇):
  `[-1/252, -1/324, -1/396, -1/468; -1/216, -131/35640, -157/51480, -61/23400; -1/198, -157/38610, -307/90090, -227/77350; -5/936, -61/14040, -227/61880, -504337/158722200]`
- **Signs.** All leading minors are positive, exactly. [C] The ratio Q₂/S on each block runs from about 1.05 to 3.14 (ψ) and 1.04 to 3.22 (ψ̇), so 16/9 is not a bound in either direction. [C, floats from exact polynomials]

**(3) Static reproductions and plant.**
- Static electric tide (Newtonian gauge, pure metric): 16/9, equal to the even Gram entry [ψ₀,ψ₀]. [C]
- Stationary magnetic tide h = 2tK (pure K): 16/9, equal to the odd entry. [C]
- ψ̇ = φ gives 7/4, NULL-CONE's two 7/4 rows. [C]
- Profile law ((n+2)/(n+1))² holds for symbolic n. [C] The static tide is the single-power n = 2 case; general data mix powers along the cone, hence the spread.
- Plant (v³f in Q₂): 16/9 becomes 4/3 and 7/4 becomes 1. Caught. Those two planted values equal the odd tower's fixed-volume numbers. [U, not pursued]

**(4) Re-run.** From `record\wave73\ruler\scripts`:
`C:\Python314\python.exe ruler.py ..\out\ruler.json` (75 s, 68 checks, 0 failures); add `--plant w3` for the plant.
Run 01 crashed in a float eigenvalue routine and carried one mis-stated check; the note is at the end of `logs\01_run.log`. The outcomes were fixed beforehand in `logs\00_prereg.md`.

**(5) Weakest step.** Positivity is shown only on four polynomial profiles per block, and because the cross block is nonzero the "metric" and "K" verdicts are statements about the t = 0 split, which the cone does not respect.

**(6) Next question.** Which single ingredient of the fixed-volume rule (maximal slice, volume constraint, or extremal edge) kills the excess on metric data — and is the planted v³f coincidence telling us it is the weight?
