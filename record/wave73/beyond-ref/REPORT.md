# BEYOND-REF: hostile referee of the state outside the vacuum's representation

## Lead's grade

**Grade: 394 of 420.**

**Defects priced (start 420).**
- −10: the direct criterion (Araki–Yamagami 1982) is stated from memory, not opened.
- −8: the divergence of the direct sum is an eight-point trend to l = 256, not a proved floor (the referee says so).
- −8: both registered runs exit 1: one of the referee's own thresholds was wrong (set at 1e-5 where the quantity is good to about 0.5%), and run 2 took 795 s against the ten-minute fence.

**Verified by the lead (no points).** Nothing re-run. The referee was not blind to the lead's header: the header-free copy did not exist yet when it started, which is the lead's fault.

**Verdicts (referee's).** (1) Uniform norm equivalence: HOLDS, with a better proof than the lane's (operator monotonicity; N ≤ h²/(l+1) for every l from the explicit function tanh^(l+1)(ρ/2)). (2) The state's modular flow is the area-side flow: HOLDS, with one citation gap (the GNS space is Fock space over the new one-particle space; three lines from memory). (3) Disjointness from the theorem: HOLDS given Theorems E/O. (4) Direct criterion, computed independently: the Hilbert–Schmidt sum over sectors grows logarithmically (about 0.64 per decade odd, 0.22 even), supporting disjointness; the per-sector affinities head to the record's constants. (5) Uniqueness: BROKEN as stated.

**What was wrong in BEYOND-NORMAL.** (a) "The only quasi-free KMS covariance" and uniqueness among all KMS states are false: a flow-invariant linear functional unbounded in the thermal norm (built from the zero-energy solution Q_l(coth ρ)) gives a one-parameter family of further KMS states and Gaussian mixtures of them, the Bose-condensate mechanism. Correct statement: unique among states continuous in the thermal norm. (b) One inequality in its norm-equivalence step is false (the minimum of x·coth(πx)/(1+x) is 0.297, below 1/π); the referee's proof replaces that step. BEYOND-NORMAL is regraded 352 (was 372) for these two defects.

**What stands, and the sentence to use.** "The area-side flow has a quasi-free 2π-KMS state on the Weyl algebra, unique among states continuous in the thermal norm, whose modular flow is that flow; it is disjoint from the vacuum representation by the theorem, and numerically by the Araki–Yamagami sum, which grows logarithmically." Never "the unique KMS state" or any wording that makes it THE state.

**Detail.** The referee's table (per-sector Hilbert–Schmidt quantity and affinity at l = 2 to 256, both parities), plant, failed lines and the re-run command are in `logs\`, `out\beyond_ref.json` and `scripts\beyond_ref.py`; re-run from `scripts\`: `C:\Python314\python.exe -X utf8 beyond_ref.py ..\out\beyond_ref.json`.
