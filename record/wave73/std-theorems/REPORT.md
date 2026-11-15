# REPORT: lane STD-THEOREMS (wave 73)

## Lead's check and grade

Sunday 15 November 2026. Literature-verification lane (Opus, web access): find the seven textbook theorems the proof lanes quoted from memory in sources that can be opened, record the statements and hypotheses, and check them against how the programme uses them. The lane returned the text below; the lead saves it. The lane wrote no files.

**Grade: 352 of 420.** Useful and honest: every item has a source, a locator and a verdict, and the list of what was NOT verified first-hand is long and explicit. Marked down for exactly that list: four of the seven items rest on secondary sources (lecture notes or a paper that reproduces the statement), two were seen only through a fetch summary or a search excerpt, and the originals (Shale 1962, Araki 1963, Kay-Wald 1991, Connes 1973, Engel-Nagel, Reed-Simon) were not opened. This lane reduces "quoted from memory" to "quoted from a reputable secondary source"; it does not make the standard theorems first-hand.

**Checked by the lead.**
- The four arXiv identifiers resolve to the papers the lane names: math-ph/0511030 (Derezinski, representations of the CCR and CAR), 2004.02383 (Hiai, lectures on von Neumann algebras), math-ph/0504012 (continuity of KMS states), 1804.00899 (absolute value of unbounded operators). The theorem numbers inside them were not checked by the lead.
- Item 7's elementary check is right: for -d^2/dx^2 + c/x^2 the solutions are x^a with a = 1/2 +- sqrt(1/4 + c), and exactly one is square-integrable at 0 when c >= 3/4.
- One correction to the lane: under item 6 it says "in hinge-proof, N is not shown to be positive". N = 1/(1 + cosh rho) is multiplication by a positive bounded function, so it is a positive bounded self-adjoint operator and Loewner-Heinz applies as HINGE-PROOF used it. No change is needed there.

**The finding that matters.** Item 1's condition: Shale's theorem is about a BOUNDED INVERTIBLE symplectic map of the whole purification space K. The proof lanes showed the flow is bounded on the region's one-particle space H, and inherited boundedness of its extension S (+) 1 to K from an older lane nobody re-checked. Because H and its complement H' are not at a positive angle, boundedness on H does not obviously give boundedness on K. This is a real gap in the chain. It has a likely clean bypass (if a unitary implements the map, the map is automatically continuous, so non-Hilbert-Schmidt or unbounded both exclude implementability); lane ASSEMBLY is working on it.

**What the paper must state as hypotheses** is the lane's own list near the end of its report; the lead adopts it, with the correction above.

---

## Lane report (unchanged)

**STD-THEOREMS report.** Five of the seven are covered as used once a stated condition is added, and two (items 2 and 4) are covered outright. None fails. Every source below is one I opened unless I say otherwise.

**1. Shale.** Source: Dereziński, "Introduction to Representations of the CCR and CAR", Lect. Notes Phys. 695 (2006), [arXiv:math-ph/0511030](https://arxiv.org/abs/math-ph/0511030), Thms 23–24, §6.3–6.4. The 1962 original at ams.org returned 403.
- Setting: a complex Hilbert space Y treated as a real symplectic space. Sp(Y) is the set of *bounded invertible* linear maps preserving the symplectic form.
- Statement: for r in Sp(Y), a Fock unitary U with UW(y)U* = W(ry) exists iff any one of these holds: j − rjr⁻¹ is Hilbert–Schmidt, [r,j] is HS, the antilinear part q is HS, or p*p − 1 is trace class. U is then unique up to a phase.
- The inverse is bounded automatically (r⁻¹ = −j r^# j), but r must be onto. "S*S − 1 is HS" is not in the list. It is equivalent, because r^#r − 1 = j⁻¹r⁻¹[j,r] when r^# is the adjoint for Re⟨,⟩.
- **Verdict: covered, with a condition.** S₀⊕1 must be a bounded bijection of all of K, symplectic for K's form, with j the Fock complex structure of K.

**2. Araki duality.** Same notes, Thm 54(6), §12.3. The proof given is Eckmann–Osterwalder's; Araki 1963 itself was not opened.
- Statement: for every closed real subspace V, M(V)′ = M(iV^⊥), and iV^⊥ is the symplectic complement. Neither standardness nor general position is needed; Thm 52 handles the general case by decomposition.
- Applied to our case: U W(k′)U* = W(k′) exactly, so U ∈ M(H′)′ = M(H), using H″ = H for closed H.
- Parts (1), (3), (4), (5) and (7) of the same theorem give: closure(H+H′) = K ⇔ H∩H′ = {0} ⇔ R(H) is a factor.
- **Verdict: covered.** It still needs H closed and S₀H = H.

**3. One-particle structure.** Kay–Wald itself is paywalled. It is reproduced with that citation in [arXiv:math-ph/0504012](https://arxiv.org/abs/math-ph/0504012), whose text I read through the fetch tool's summary rather than the PDF itself.
- Statement there: σ² ≤ 4μμ; Re⟨Kφ,Kφ′⟩ = μ; σ = 2Im (so Im = σ/2); KS + iKS is dense; ω(W(φ)) = e^{−μ(φ,φ)/2}; KS dense ⇔ the state is pure.
- Uniqueness up to unitary is *not* stated there. First-hand, Dereziński's Thm 41 gives the Araki–Woods purification, assuming σ is nondegenerate on the completion.
- **Verdict: covered, with a condition.** Use the Kay–Wald normalisation. Dereziński's Thm 39 uses e^{−¼yηy}, so do not mix factors between the two.

**4. Connes cocycle.** NUMDAM hosts Connes 1973 openly, but the PDF exceeded the fetch limit. First-hand source: Hiai, [arXiv:2004.02383](https://arxiv.org/abs/2004.02383), Thm 7.6 (citing Connes, Bull. Sci. Math. 97 (1973)) and Thm 12.20 (Connes' converse, ASENS Thm 1.2.4).
- Statement: for faithful normal semifinite weights φ, ψ, the cocycle u_t = (Dψ:Dφ)_t is a strongly* continuous family of unitaries in M, with σ^ψ_t = Ad u_t ∘ σ^φ_t and u_{s+t} = u_s σ^φ_s(u_t).
- What follows: if α = σ^ψ, then α_t∘σ^φ_{−t} is inner **for every t**. So failure at one single t is enough. Continuity and the cocycle identity are conclusions of the theorem, not inputs we need.
- What does not follow: the converse. Being inner at each t yields a weight only with a strongly continuous choice of unitaries obeying the cocycle identity (Thm 12.20).
- "Normal state" here means a faithful one.
- **Verdict: covered.**

**5. Bounded perturbation.** The Engel–Nagel host refused the connection. First-hand source: Schnaubelt, *Evolution Equations* lecture notes (KIT, 2023), Thm 3.4, which cites Engel–Nagel §III.1.
- Statement: if A generates T(t) with ‖T(t)‖ ≤ Me^{ωt} and B is bounded, then A+B on D(A) generates S(t) with ‖S(t)‖ ≤ Me^{(ω+M‖B‖)t}. Both Duhamel formulas and the Dyson–Phillips series hold, and the solution is unique.
- These notes assume *complex* Banach spaces.
- Transfer to our case: get the group by applying the theorem to ±(A+B). Then complexify. On a real Hilbert space ‖B_C‖ = ‖B‖, and each Dyson–Phillips term maps X into X, so the bound and Duhamel restrict to X.
- **Verdict: covered, with a condition.** The paper must state the complexification remark. With A generating isometries, the bound is e^{‖B‖|t|}.

**6. Löwner–Heinz.** Secondary only: Boucif–Dehimi–Mortad, [arXiv:1804.00899](https://arxiv.org/abs/1804.00899) §2, citing Schmüdgen GTM 265 without a number.
- Statement there: S ≥ T ≥ 0 in the form sense (D(S^½) ⊂ D(T^½) and ‖T^½x‖ ≤ ‖S^½x‖) implies S^½ ≥ T^½.
- Applied with S = B², T = A²: the hypothesis means D(B) ⊂ D(A) and ‖Ax‖ ≤ ‖Bx‖.
- **Verdict: covered, with a condition.** A and B must be self-adjoint, not merely symmetric. In hinge-proof, N is not shown to be positive. Apply the theorem to |N| to get ±N ≤ h/λ, which is what the next line really uses.

**7. Limit point / limit circle.** Reed–Simon was not opened; Thm X.10 is quoted in [arXiv:gr-qc/0401009](https://arxiv.org/abs/gr-qc/0401009), Thm 3, which I saw only as a search excerpt.
- Statement: for V continuous and positive near 0, V ≥ ¾x⁻² gives limit point, and V ≤ (¾−ε)x⁻² gives limit circle.
- Elementary check: the solutions are x^a with a = ½ ± √(¼+c), so the endpoint is limit point iff c ≥ ¾.
- Our case: λ²/sinh²ρ = λ²/ρ² − λ²/3 + O(ρ²), and λ² ≥ 2, so the bound holds near 0.
- **Verdict: covered.** The operator must first be put in Liouville normal form on L²(dρ).

**What the paper must state as hypotheses**
- S₀⊕1 is a bounded bijection of K that preserves its symplectic form; S₀ maps H onto H.
- H is closed and σ is nondegenerate on it. This is equivalent to K = closure(H+H′) and to R(H) being a factor.
- The Kay–Wald normalisation (Im = σ/2, e^{−μ/2}) is used throughout.
- The states considered are faithful; the argument uses only Connes' forward direction.
- The perturbation theorem is transferred to the real space by complexification.
- In Löwner–Heinz, the operators are positive self-adjoint and ordered in the form sense; |N| is used.
- The radial operator is in Liouville normal form, and l ≥ 1.

**Not verified first-hand**
- Shale 1962, Araki 1963, Leyland–Roberts–Testard.
- Kay–Wald Prop. 3.1, including its uniqueness clause.
- Connes' theorem numbers 1.2.1 and 1.2.4, which I have from Hiai and other papers' citations.
- Takesaki VIII.3.3, Bratteli–Robinson, Ruijsenaars, and the theorem numbers in the Dereziński–Gérard book.
- Engel–Nagel III.1.3, Pazy and Kato numbering, Schmüdgen's numbering, and the text of Reed–Simon X.10.

I wrote no notes file and touched nothing in git.