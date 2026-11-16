# CHARGE-REF: hostile referee of the sector identity

## Lead's grade

**Grade: 408 of 420.**

**Defects priced (start 420).**
- −6: the Dirichlet kernel used in the independent probe is quoted from memory (its jump condition was checked by hand).
- −6: the extension beyond smooth data (closure, the Hardy-type bound) is a sketch.

**Verified by the lead (no points).** Re-ran `C:\Python314\python.exe -u scripts\ref_probe.py 16 exp <out.json>` from the referee's folder: "pass 15 fail 0 plant CAUGHT", 7.7 s, exit 0. The referee read a header-free copy of the lane's text.

**What it changes.** The sector identity survives a hostile referee in every part, and one of its two inherited inputs ("moving the region is a gauge move; fixed volume") is no longer inherited: the referee proves it. CHARGE-ALL-L is regraded 398 (was 386): what remains against it is three wording gaps and the inherited completeness of the data family. Wording to adopt everywhere: the charge is the second-order part of Area − 2·Volume taken at its SADDLE over three freedoms: the edge's shape, the corner's time shift, and the slice (maximal). Never "the ball's area is extremal (maximal or minimal) at fixed volume"; never the four-term area formula on data that are not transverse-traceless.

**Still inherited.** That the two families T[q], S[η] exhaust the smooth transverse-traceless data of a sector (companion A8 Lemma 7; CHAIN-ALL-L); that c_l E_1 is the modular generator.

---

## Lane report (unchanged)

**CHARGE-REF referee report on CHARGE-ALL-L**

I could not break it: every part and both inherited items hold, with three wording gaps. All registered lines passed (70 of 70, plants caught); pass lines were fixed in `logs\00_passlines.log` before the first run.

## Verdicts

| Item | Verdict |
|---|---|
| (a) odd, all l | HOLDS |
| (b) even momentum block | HOLDS |
| (c) even field block | HOLDS |
| (d) cross block | HOLDS, GAP in wording |
| Inherited (i): area formula on TT data | HOLDS |
| Inherited (ii): region move is a gauge move; fixed volume | HOLDS, now [P]; GAP in wording |

## Defects, most serious first

1. **The slice is a third freedom, and the stationary point is a saddle (GAP).** "Stationary in the edge shape and the corner's time shift" omits the bulk slice X. Stationarity in X means ΔX = 0, the maximal slice [P]. In A the point is a maximum in σ, a maximum in the slice and a minimum in τ [P]. Repair: add "and the slice (maximal, tr K = 0); a saddle". Z ≥ 0 comes from the formula, not from a variational principle.
2. **Step 5 need not be inherited (GAP, closed here).**
   - Space: the scalar curvature vanishes at orders 0 and 1, so its second-order coefficient is unchanged pointwise by a first-order diffeomorphism [P].
   - Time: the slice t = εX has induced second-order metric h₂ + 2XK − ∂X∂X, and R_lin[2XK − ∂X∂X] = 2K·Hess X + |Hess X|² − (ΔX)² pointwise [P; C, E01].
   - The gauge route equals the direct route for the edge [C, A05, l = 2, 3] and for the corner [P; C, P05].
   - The drive −∮τκ and the stiffness ½(l−1)(l+2) are confirmed by hand [P].
   - Fixed volume: the rescaling needed is second order and F is stationary under it [P].
3. **(d) is stated too narrowly (GAP).** Time reflection does kill hK, hτ, hX, Kσ, στ and σX with the edge moving in time [P]. The even–odd and m–m′ cross terms need one more sentence (parity and Schur).
4. **The area formula is right term by term** [P by hand; C, A01], and the coarea form gives the same second-order area [P]. Off TT gauge it is false: h = −2φδ at l = 2 gives −1/7 against the true −6/7 [C, A03]. Repair: quote the gauge-free line A = ∫N(K·K − (tr K)²) − ∫N R⁽²⁾ − 2δ²F whenever the data are not TT.
5. **"Profiles smooth in r²" undersells the result.** The potentials are identities for arbitrary q, so C³ profiles bounded at 0 suffice [P]. Beyond smooth data A_ext exists only by closure, and Z_E ≤ C·c_l E_1 by a Hardy-type bound [P, sketch].

Step 10 holds: r^(2l+3)U·l(l−1) = (r^(2l+1)q²)′ + 2r^(2l+2)q′²/((l+1)(l+2)), the corner cancels q(1)² exactly, and r^(2l+2)U → 0 at 0 [P].

## Probe

Own code, a harmonic of order m = l−2, TT data from an exact null-space solve (not the lane's f, g formulas), and the textbook Dirichlet kernel [S, from memory, jump condition checked by hand]. Values are per N_Y, entry AA.

| l, weight | odd K (a) | even K (b) | Z_E / vacuum | lhs − rhs |
|---|---|---|---|---|
| 15, polynomial | ≈ 514.84 | ≈ 4.9507 | ≈ 0.03019 (< 1/33) | exactly 0 [C] |
| 16, polynomial | ≈ 605.92 | ≈ 4.9308 | ≈ 0.02848 (< 1/35) | exactly 0 [C] |
| 15, exp(−r²) | 5.0414e11 | 34.4927 | 0.0303019 | ≤ 1e-49 [N] |
| 16, exp(−r²) | 9.1692e11 | 33.9867 | 0.0285705 | ≤ 1e-49 [N] |

- The odd h and even h blocks pass on all three entries (AA, AB, BB) in every run, as does the corner identity with a generic non-harmonic X.
- **Plant:** the free-space Green function in place of the Dirichlet one was caught in all four runs (for example −7676100608/151737987629655 at l = 15).
- **Failed lines:** none.
- **Unregistered extras:** A00 (my curvature code returns R = 0 on an exactly flat pulled-back metric) and E01 with its own plant, caught.

## Re-run

From the referee's folder:

```
C:\Python314\python.exe -u scripts\ref_probe.py 16 exp out\probe_l16_exp.json
```

The same script takes `15` or `16` and `poly` or `exp`. The other two are `scripts\ref_area.py out\area_formula.json` (about 4 minutes) and `scripts\ref_slice.py`.

## Weakest remaining step

Completeness: that T[q] and S[η] exhaust the smooth TT data of a sector is still inherited [I, A8 Lemma 7, CHAIN-ALL-L]. So is the claim that c_l E_1 is the modular generator [I].

## Sentences

- **May use:** "For every l ≥ 2 and either parity, the second-order part of Area − 2·Volume, taken at its saddle point over the edge's shape, the corner's time shift and the slice, equals c_l E_1 of the chain field plus a non-negative form on extrinsic-curvature data."
- **Must not use:** "The ball's area is extremal (maximal, or minimal) at fixed volume", nor the four-term area formula on data that are not transverse-traceless.
