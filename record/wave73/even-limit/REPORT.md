# REPORT: lane EVEN-LIMIT (wave 73)

## Lead's check and grade

Sunday 15 November 2026. Small proof lane (Fable): the first step of the even-parity analogue of HINGE-PROOF's Lemma B, the large-l limit of the even charge's operator. The lane returned the text below; the lead saves it.

**Grade: 382 of 420.** One question, answered more strongly than asked (Hilbert-Schmidt convergence, not only strong), with the inherited limit kernel confirmed rather than assumed, pass lines fixed before the run, and a planted error caught. Marked down because the 1/lambda rate (parts b and c) rests on constants the lane says it did not write out to the last inequality, and because its domination check L6 only shows "no violation" (both sides underflow to 0 at the worst point).

**Checked by the lead, by hand: steps 1 to 6, which give statement (a). They hold.**
- Step 1: with t = e^(-tau), mu = l + 1/2, xi = mu*tau, the kernel t_<^(l+1)(t_>^(-l) - t_>^(l+1))/(2l+1) equals (t t')^(1/2) sinh(xi_<) e^(-xi_>)/mu, and -ln tanh(rho/2) = 2 artanh(e^(-rho)) = 2 artanh(x/2 lambda).
- Step 2: d xi/dx = (mu/lambda)/(1 - (x/2 lambda)^2) >= mu/lambda > 1.
- Step 3: (t t')^(1/2) = 1 - (x + x')/(2 lambda) + O(lambda^-2).
- Step 4: sinh(a) e^(-b) <= (1/2) e^(-(b-a)) and <= a for a <= b, so the kernel is dominated by min(2 x_<, (1/2) e^(-|x-x'|)).
- Step 5, first region: the integral of 4 e^(-2 max(y,y')) over y, y' > 0 is 2. The other two regions were read, not recomputed.
- Step 6: dominated convergence in L^2(dy dy').

**Re-run by the lead:** `C:\Python314\python.exe even_limit.py <out.json>` from `scripts\`, 65 s, exit 0, "ALL LINES HELD"; the table equals the lane's to the printed digits (slopes -0.9982 and -0.9977, plant ratios 197.5 and 1275.5).

**Not checked by the lead:** steps 7 and 8 (the rate and the leading correction), beyond their agreement with the table.

**What it changes.** Together with EVEN-PROOF the even half now has: the uniform bound (proved), and the limit operator (proved, statement a). Still missing for the even divergence as a theorem: passing the first-order matrix element to the limit (needs strong resolvent convergence of h_l), and that the limit is non-zero (the even constant positive as a theorem). The paper must still say "computed, not proved" for the even divergence.

---

## Lane report (unchanged)

**EVEN-LIMIT report (2026-11-15)**

hinge-e's limit kernel is correct as stated, and the convergence is stronger than asked: Hilbert–Schmidt norm, not merely strong. All seven pre-fixed lines held on the first run; none failed.

## 1. Statement

Let λ = √(l(l+1)), μ = l+½, y = ρ − ln 2λ, x = e^(−y). Extend W_l by zero to L²(ℝ, dy), since ρ > 0 means y > −ln 2λ. Let W_inf have kernel sinh(x_<)e^(−x_>) with respect to dy dy′, so w = W_inf f solves −w_xx + w = f/x, w(0) = 0. No factor or measure correction is needed.

- **(a)** W_inf is Hilbert–Schmidt and ‖λW_l − W_inf‖_HS → 0, hence also in operator norm and strongly on all of L²(ℝ, dy).
- **(b)** For f in L² with support in a fixed compact y-interval, ‖λW_l f − W_inf f‖ ≤ C(supp f)‖f‖/λ.
- **(c)** For such f, λ(λW_l − W_inf)f → −½(XW_inf + W_inf X)f in L², with X = multiplication by e^(−y). So the rate 1/l is sharp whenever that vector is non-zero.
- A uniform rate is not established: I expect O(λ^(−1/2)) in HS norm and O(1/λ) in operator norm, both [U].

## 2. Proof

1. [P] Exact rewriting. With τ = −ln t = 2 artanh(x/2λ) and ξ = μτ: λW_l(y,y′) = (λ/μ)(tt′)^(1/2) G(ξ,ξ′), where G(a,b) = sinh(a_<)e^(−a_>) = ½e^(−|a−b|)(1 − e^(−2a_<)). The smaller t corresponds to the larger x.
2. [P] ξ(x) = 2μ artanh(x/2λ) has dξ/dx = (μ/λ)/(1−(x/2λ)²) ≥ 1, so δ = ξ − x ≥ 0 is increasing and |ξ−ξ′| ≥ |x−x′|. For x ≤ λ, δ ≤ x/(8λ²) + c·x³/λ² and ξ ≤ 1.17x.
3. [P] Pointwise limit. For fixed x, x′: ξ → x and (λ/μ)(tt′)^(1/2) = 1 − (x+x′)/(2λ) + O(λ⁻²), so λW_l → G(x,x′) and λ(λW_l − G) → −½(x+x′)G.
4. [P] Domination. 0 ≤ λW_l ≤ D := min(2x_<, ½e^(−|x−x′|)), from step 2, λ/μ < 1 and t ≤ 1; for x_< > λ use ½ ≤ 2x_<. The same bound holds for G.
5. [P] D is in L²(dy dy′) = L²(dx dx′/(xx′)), in three regions:
   - x_> ≤ 1: D² ≤ 4x_<², and ∬4e^(−2max(y,y′)) = 2.
   - x_< < ½ ≤ x_>: D² ≤ x_< e^(½−x_>), integrable.
   - both ≥ ½: ∫e^(−2|x−x′|)dx′/x′ ≤ C/x, and ∫dx/x² < ∞.
6. [P] Steps 3 to 5 with dominated convergence give (a). Below the cut the difference is −G ≤ D.
7. [P] (b). For x′ in [a,b] ⊂ (0,λ], split ΔK = [(λ/μ)(tt′)^(1/2) − 1]G(ξ,ξ′) + [G(ξ,ξ′) − G(x,x′)].
   - First term: at most (1.1(x+b)/2λ + 1/8λ²)·D for x ≤ λ.
   - Second term: at most (δ + sinh(x)δ′) ≤ C_b·x/λ² for x < x′, and e^(−(x−x′))(δ+δ′) for x′ < x ≤ λ.
   - For x > λ, including below the cut, everything is at most e^(−(x−b)).
   - Each bound is in L²(dx/x) uniformly in x′ with norm C_{a,b}/λ; Cauchy–Schwarz in y′ finishes.
8. [P] (c). Step 3's pointwise limit with step 7's bounds as dominant.
9. [N] Sharpness and constants: table below.

## 3. Check

err = ‖λW_l f − W_inf f‖ on L²(dy). f1 is a bump on |y| < 1; f2 is a bump times cos 3y on (−0.5, 2.5). F = −½(XW_inf + W_inf X)f.

| l | err f1 | λ·err/‖F‖ | plant f1 | err f2 | λ·err/‖F‖ | plant f2 |
|---|---|---|---|---|---|---|
| 16 | 1.0257e-2 | 0.9787 | 4.59e-2 | 8.676e-4 | 0.9642 | 1.67e-2 |
| 64 | 2.6652e-3 | 0.9945 | 3.63e-2 | 2.279e-4 | 0.9907 | 1.80e-2 |
| 256 | 6.7294e-4 | 0.9986 | 3.39e-2 | 5.772e-5 | 0.9976 | 1.84e-2 |
| 1024 | 1.6866e-4 | 0.9997 | 3.33e-2 | 1.448e-5 | 0.9994 | 1.85e-2 |

- Slopes from l = 256 to 1024 are −0.998 and −0.998.
- ‖λ·diff − F‖/‖F‖ at l = 1024 is 4.0e-4 and 6.5e-4.
- The kernel code matches the literal formula in 60-digit mpmath to 1.4e-13.
- Refining the quadrature changes err by at most 2e-14 relative.
- **Plant caught:** with the shift ln λ instead of ln 2λ, the error stalls at 3.3e-2 and 1.8e-2, which is 198 and 1276 times the true error at l = 1024.
- L6, the sampled domination check, passed, but its reported maximum is exactly 0 (both sides underflow at some sample points), so it only shows no violation.

## 4. Re-run

From `record\wave73\even-limit\scripts`:

`C:\Python314\python.exe even_limit.py ..\out\even_limit.json` (65 s, exit 0)

The pass lines were written to `logs\01_prefixed_pass_lines.log` before the run; the run log is `logs\02_even_limit_run1.log`. Nothing was committed and no REPORT.md was written.

## 5. Weakest step

Step 7's constants were bounded by hand and not written out to the last inequality; the 1/λ rate rests on that sketch plus the table, whereas (a) needs only steps 1 to 6.

## 6. What remains for Lemma B

Combine (a) with strong resolvent convergence of h_l and the uniform bound ‖λh_lW_l‖ ≤ 1 from EVEN-PROOF to pass the first-order matrix element on a fixed wave packet to the Liouville limit. Then show that limit is non-zero, which means hinge-e's constant a_E > 0 as a theorem; that part is untouched here.
