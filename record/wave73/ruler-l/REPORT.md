# REPORT: lane RULER-L (wave 73)

## Lead's check and grade

Sunday 15 November 2026. Small exact-computation lane (Fable): repeat RULER's light-cone computation at l = 3 to 6, and turn "positive on four profiles" into a statement for all profiles. The lane returned the text below; the lead saves it. It sent two interim results to the lead before this report and read the bulletin through B22, as the protocol asks.

**Grade: 386 of 420.** The four-profile observation is replaced by a one-line identity valid for every cone curvature and every l, with the infimum and supremum of the ratio settled and a "must not say" that withdraws a number the earlier lane's profiles had suggested; the higher-l runs pass every gate; and it found, diagnosed and documented a real fault in machinery inherited from RULER (a simplification call that turned long exact rationals into algebraic guesses), voiding its own first runs rather than reporting them. Marked down because the l = 5 and 6 runs exceeded the brief's ten-minute guide (run four in parallel, not timed alone), the combined re-run command was never run as one job, and the link from "non-zero cone curvature" to "non-trivial perturbation in the diamond" is checked only on the eight-dimensional test spaces, as the lane states.

**Checked by the lead, by hand. It holds.**
- Steps 1 and 2: with sigma_p = -v^(-p) int_0^v u^p C du and h = v^2 (sigma_2 - sigma_3) = -int u^2 C + v^-1 int u^3 C, differentiation gives h' = -v^-2 int u^3 C = v sigma_3, hence sigma_3 = h'/v and sigma_2 = h'/v + h/v^2.
- Steps 3 to 5 (R = 1): Q_p is proportional to int v^3 (1 - v) sigma_p^2, so Q_3 ~ int v(1-v) h'^2 and Q_2 ~ int v(1-v)(h' + h/v)^2 = Q_3 + int (1-v)(h^2)' + int (1-v) h^2/v; the middle term is int h^2 after integration by parts; the sum is int h^2/v = int v^3 (sigma_2 - sigma_3)^2.
- On monomials C = v^(n-2) both sides equal 1/((n+1)^2 (n+2)^2 (2n+2)) (checked earlier, on the interim result).
- The static ratios 25/16, 36/25, 49/36, 64/49 are ((l+2)/(l+1))^2 at l = 3 to 6.
- Not hand-checked: step 8's estimates for the unbounded family; the l = 3 to 6 Gram entries.

**Re-run by the lead:** `C:\Python314\python.exe ruler_l.py --l none <out.json>` (part 2), 39 s, exit 0, 8 checks, 0 failures; the unbounded family's ratio grows from 4.0 to 8.8 as N goes from 8 to 128. The l = 3 to 6 runs (five to thirteen minutes each) were NOT re-run by the lead; their stored outputs are in `out\`. The plant run was not repeated.

**Note for the record.** RULER's published l = 2 numbers are unaffected by the simplification fault (reproduced exactly with the fixed script), but `ruler\scripts\ruler.py` still contains the call; anyone extending it to higher l must use this lane's script.

**What it changes.** "Which ruler" is now a theorem about the cone forms: the cone first-law charge strictly exceeds the entropy for every non-zero cone curvature at every l, identically in both parities, with the ratio unbounded above and with infimum 1 not attained. Nothing about the cone rule is special to four profiles or to l = 2. The earlier impression that the excess lies "between a few percent and a factor of three" is withdrawn.

---

## Lane report (unchanged)

**RULER-L final report** (bulletin read through B22, no change of course, except that I adopt B22's wording and say "non-zero cone curvature", never "non-zero datum").

**(1) Answer.** At l = 3, 4, 5, 6 the result is the same as at l = 2: every gate passes in both parities, and Δ = Q₂ − S is positive definite on the ψ block, the ψ̇ block and the full 8×8 space, with a non-zero cross block. Positivity is now an identity for all profiles: Δ = (1/4G)∮∫₀ᴿ v³|σ₂ − σ₃|² dv, so Δ > 0 whenever C is not identically zero; inf Q₂/Q₃ = 1 (not attained) and sup = +∞.

**(2) Per-l results [C].** Stopped at l = 6. Both parities give identical numbers.

| l | checks / fail | static ratio | Δ[ψ₀,ψ₀] | Δ[ψ̇₀,ψ̇₀] | cross [0,0] | Q₂/S on ψ̇ = φ | Q₂/S range, ψ block | Q₂/S range, ψ̇ block |
|---|---|---|---|---|---|---|---|---|
| 2 | 40 / 0 (incl. part 2) | 16/9 | 1/36 | 1/1008 | −1/252 | 7/4 | 1.050–3.135 | 1.038–3.217 |
| 3 | 31 / 0 | 25/16 | 1/240 | 1/10800 | −1/2160 | 31/20 | 1.042–2.585 | 1.032–2.666 |
| 4 | 31 / 0 | 36/25 | 1/900 | 1/59400 | −1/9900 | 43/30 | 1.036–2.257 | 1.028–2.331 |
| 5 | 31 / 0 | 49/36 | 1/2520 | 1/229320 | −1/32760 | 19/14 | 1.032–2.040 | 1.025–2.106 |
| 6 | 31 / 0 | 64/49 | 1/5880 | 1/705600 | −1/88200 | 73/56 | 1.028–1.885 | 1.022–1.945 |

- **Gates at every l:** all 16 data are TT and solve the vacuum equations; the Q₃ Gram equals c_l·E₁ and is block diagonal; even and odd Grams are identical; all leading minors are positive (4, 4 and 8).
- **New gate:** Δ equals the Gram of 2∮∫v³(σ₂−σ₃)² exactly on each 8×8.
- **Observed:** C starts at v^(l−2) on every datum.

**(3) The inequality.** Per generator, R = 1, with C satisfying ∫₀u²|C|du < ∞:

1. [P] Define h = v²(σ₂ − σ₃) = −∫₀ᵛu²C du + v⁻¹∫₀ᵛu³C du. Then |h(v)| ≤ 2∫₀ᵛu²|C| → 0, so h(0) = 0, and h(1) is finite.
2. [P] h′ = −v⁻²∫₀ᵛu³C = vσ₃, so σ₃ = h′/v and σ₂ = h′/v + h/v².
3. [P] Q₃ ∝ ∫v(1−v)h′², and Q₂ ∝ ∫v(1−v)(h′ + h/v)² = Q₃ + ∫(1−v)(h²)′ + ∫(1−v)h²/v.
4. [P] ∫₀¹(1−v)(h²)′ = [(1−v)h²]₀¹ + ∫h² = ∫h². Work on [ε,1] and let ε → 0; the boundary term (1−ε)h(ε)² → 0.
5. [P] So Δ ∝ ∫h² + ∫(1−v)h²/v = ∫₀¹h²/v dv ≥ 0. For general R the two 1/R terms cancel in the same way. [C: symbolic exponents; exponents 0..40; each 8×8 Gram, l = 2..6]
6. [P] Δ = 0 iff h ≡ 0 iff σ₃ ≡ 0 iff C = 0 almost everywhere.
7. [P][C] Infimum: C = vᵐ gives ((m+4)/(m+3))² → 1, and by step 6 the value 1 is not attained. Exact for m = 0..8: 16/9, 25/16, 36/25, 49/36, 64/49, 81/64, 100/81, 121/100, 144/121.
8. [P] Supremum: take h′ = v²(1−v)ᴺ, that is C_N = 4(1−v)ᴺ − Nv(1−v)ᴺ⁻¹. Then Q₃ ∝ B(6, 2N+2) ≤ (15/8)N⁻⁶. Also h(1/N) ≥ 1/(12N³) and h is increasing, so Δ ≥ ln N/(144N⁶). Hence Q₂/Q₃ − 1 ≥ ln N/270 → ∞. With h′ = vˡ(1−v)ᴺ the same holds inside v^(l−2)·polynomials, so it holds at every l.
9. [C] Exact family values: 2.27, 3.22, 4.99, 7.43, 8.79 at N = 1, 4, 16, 64, 128.
10. [C] For polynomial C of degree ≤ 8 the exact range is 1.011980867 to 5.724642 (bracketed by exact minors); at degree ≤ 12 the maximum is 7.31.
11. [N] Random search over 20,000 integer polynomials in each of two bases: 1.056 to 3.217. It stays inside the exact range but reaches neither end.
12. [I] The cone forms themselves, S = Q₃, and the duality argument come from NULL-CONE.
13. [U] That C ≡ 0 on the cone forces zero curvature in the diamond for general data. It is checked only on the 8-dimensional spaces, through the positive-definite 8×8.

**(4) Plant and failed gates.**
- **Plant:** `--plant swap` at l = 2 gives 12 failures, including "Q₃ Gram equals c_l E₁"; the static ratio 16/9 becomes 9/16. Caught.
- **Void first attempt:** the first l = 3..6 runs (logs `03_*`) are void. `sp.nsimplify`, inherited from `ruler.py`, turned long exact rationals into algebraic guesses and failed gates at l = 3 and 4. I replaced it with exact cancellation plus an is-rational assertion; see `logs\04_note_nsimplify.md`.
- **RULER itself:** its l = 2 numbers are unaffected; they are reproduced exactly with the fixed script.
- **Stopped process:** I stopped one hung process, this lane's own void l = 6 run.

**(5) Re-run.** From `record\wave73\ruler-l\scripts`:

`C:\Python314\python.exe ruler_l.py ..\out\all.json`

This runs l = 2..6 and part 2 in sequence. That combined command was not run as one job. I ran `--l N --no-part2` per l, four in parallel: 307, 297, 804 and 700 s for l = 3..6 under that load, so l = 5 and 6 exceeded the 10-minute guide and were not timed alone. Part 2 alone is `--l none` (40 s). The plant is `--l 2 --plant swap`.

Files are in `record\wave73\ruler-l\`:
- `scripts\ruler_l.py`
- `logs\00_prereg.md`
- `logs\08_outcome.md`
- `out\ruler_l_3.json` (and `_4`, `_5`, `_6`, `_2_part2`)

**(6) Weakest step.** The theorem is about the cone forms as given: S = Q₃ and the normalisations are inherited, and "non-zero cone curvature" is linked to "non-trivial vacuum perturbation" only on the 8-dimensional test spaces.

**(7) For the paper.**
- **May say:** "Given the cone forms, Q₂ − S = (1/4G)∮∫₀ᴿ v³|σ₂ − σ₃|² dv, so the cone charge strictly exceeds the entropy for every non-zero cone curvature at every l, identically in both parities; the ratio has infimum 1, not attained, and no upper bound."
- **Must not say:** that ((l+2)/(l+1))² or any other number bounds Q₂/S, or that the excess is "between 4% and a factor 3" — those were properties of four profiles.

**Part (3) of the brief, in one sentence.** Yes: given the cone forms, "balances on no data" is a theorem for every l (Δ = 0 only for C ≡ 0), and "prefers neither half" follows from it — positivity holds on any subspace, metric or extrinsic-curvature — while the stronger "treats both halves identically" rests on NULL-CONE's duality [I], checked exactly here at l = 2..6 [C].

The inequality turned out to be one integration by parts, which is a humbling return on four profiles and 75 seconds per run. And nsimplify, asked to simplify a fraction, invented 2^(27/55) — the first time this programme has been lied to by its own calculator rather than by a referee.
