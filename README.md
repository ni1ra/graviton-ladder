# The graviton ladder

**Entropy cannot tell electric tides from magnetic ones. Area can.**

Take a ball in flat space and put a weak tidal gravitational field through it. Two things respond at second order in the field:

- the **entropy side**: the relative entropy of the gravitons inside the ball (how well an observer confined to the ball's causal diamond can tell the tidal field from empty space);
- the **area side**: the drop of the stationary value of Area − 2·Volume, divided by 4G.

The first law of causal diamonds, carried to second order, would make them equal. This repository contains an exact computation of both in free linearised gravity, with every check rerunnable in one command.

| Tidal field of multipole order l | Entropy | Area side | First law at second order |
|---|---|---|---|
| static electric | l + 1 | l + 1 | holds |
| stationary magnetic | l + 1 | l + 2 | misses by one part in l + 2 |
| (Rindler wedge, any data) | equal rates | equal rates | no miss |

For a uniform tidal field in a ball of radius R (l = 2, ħ = c = 1) the entropy is the Bel–Robinson super-energy of the ball:

$$S=\frac{\pi R^{6}}{420\,G}\left(\mathcal E_{ij}\mathcal E_{ij}+\mathcal B_{ij}\mathcal B_{ij}\right)$$

while the area side is $\frac{\pi R^{6}}{420\,G}\,\mathcal E_{ij}\mathcal E_{ij}$ for a purely electric field and $\frac{\pi R^{6}}{420\,G}\cdot\tfrac43\,\mathcal B_{ij}\mathcal B_{ij}$ for a purely magnetic one.

The entropy is blind to electric–magnetic duality. The area is not. No local Lorentz scalar can remove the mismatch (a scalar can move it between the two sectors, never cancel it). It is exactly the failure of the ball's conformal Killing vector to be a Killing vector: on the Rindler wedge, where the vector is Killing, the same computation gives no mismatch.

**What is new and what is not.** The area numbers at l = 2 are already in print: Jacobson, Senovilla and Speranza, [arXiv:1710.07379](https://arxiv.org/abs/1710.07379), eq. (75) at the constant-mean-curvature shape. They also showed the area deficit is not the Bel–Robinson energy for any shape. What we did not find in print is the entropy side, and so the comparison. "Holds on electric data" is true at the stationary surface; for a geodesic ball the electric case fails too, while the magnetic coefficient is the same for every shape.

**Start here:** [`paper/ladder_short.pdf`](paper/ladder_short.pdf) (8 pages, version 3).

## Status

Read this before quoting the result.

- **Proved here, with scripts:** the entropy formula (two routes), the area side on magnetic data, the area side on electric data from the exact Weyl solution, the wedge identity.
- **Assumed, and stated as assumptions:** the theory is free linearised gravity at order G⁰ on a flat ball in four dimensions; the gravitational term of the first law at this order is the stationary value of Area − 2·Volume; three conditions on what counts as "the algebra of the ball"; Einstein's equations to second order. The full list is in the paper, section 8, and in [`record/wave69/final/THEOREM.md`](record/wave69/final/THEOREM.md), section 2.
- **Open:** whether the ball algebra of full quantum gravity reduces to the free one as G → 0. Nothing here tests that.
- **Not peer reviewed.** Every check so far was made inside this programme. No outside referee has read the proofs.
- **Corrected once already.** Version 1 of the paper was public for about an hour before a blind read came back. It confirmed every formula and found four statements that were wrong or too broad, including that part of the area result is in print. Version 2 fixed them; version 3 only shortens the text. The read is in [`record/wave72/blindread/`](record/wave72/blindread/REPORT.md).

## Since version 3: newer work, not yet in the paper

The folders under [`record/wave73/`](record/wave73/) hold work done after version 3. A version 4 of the paper is in preparation. Until it is here, the paper and the summary above do not reflect this work; where they differ, the newer record is the more careful statement. None of it is run by `./run_all.sh` yet, and each report lists what its lead did and did not re-check.

What the newer record supports, with its stated gaps:

- **The contrast is metric data against extrinsic-curvature data, not electric against magnetic.** On every graviton sector of the flat ball, of either parity, the gauge-invariant fixed-volume area-side charge equals the vacuum's modular generator on metric data and exceeds it on extrinsic-curvature data by a positive form (local for odd parity, non-local for even). The ladder above is that statement restricted to stationary magnetic data. Exact checks are at l = 2 to 6. See `even-tower`, `charge-id`, `subspace`, `jss-id`.
- **Odd parity.** Granting the reduction of the odd sectors to scalar sectors, the flow of that charge is not the modular flow of any normal state or weight in the vacuum sector: the Shale sum over sectors diverges logarithmically, with per-sector constant 7ζ(3)/(16π²) ≈ 0.0533. See `hinge-proof`. The constant was reproduced by a lattice run that was meant to be blind and was weakened in one way its report states (`hinge-blind`). That proof left its first-order lower bound as an outline resting on a numerically evaluated identity; `odd-sign` replaces the outline by a written argument on an explicit family, uniform in l. The second-order remainder on that family is bounded in `odd-remainder`, so the odd chain no longer uses the earlier outline or its unchecked corollary.
- **Even parity.** The per-sector constant is about 0.0136. A pre-registered lattice measurement was inconclusive under its own gate; a secondary analysis, fixed mid-run after about three quarters of the data, gives 0.01359 with a spread of 0.3% (`even-blind`). The even divergence is **not proved**. Established so far: the operator bound, the large-l limit of the operator, a non-zero first-order term, a lower bound on it uniform in l, a first-principles definition of the antilinear part, and a bound on the second-order remainder (`even-proof`, `even-limit`, `even-nonzero`, `even-passage`, `beta-def`, `even-remainder`). Each step has a written proof whose core the lead checked by hand. Three gaps named in earlier versions of this paragraph are closed: a blind second code reproduces every certified constant, each on the safe side (`even-second`); the density and domain facts the steps assumed are proved (`even-density`); and the change of variables that identifies the operator with the even charge is an exact identity, with coupling 1 (`even-dict`). The assembled statement and what it rests on are under "Both parities, assembled" below.
- **Both parities, assembled.** `assembly` states the result once per parity with its hypotheses sorted: for every non-zero coupling, the flow generated by the vacuum's generator plus the extra term is, at modular time 1/2, implemented by no unitary in the region's algebra, so it is the modular group of no faithful normal semifinite weight and of no faithful normal state. The link that was only inherited (Shale's criterion needs a bounded map of the whole purification space) is closed there: a unitary implementer forces the map to be continuous, so no boundedness input is needed. What the statement rests on: (1) steps proved in this record, each with a written proof whose core the lead checked by hand; (2) standard theorems (Shale's necessity direction, Connes' cocycle theorem, the bounded perturbation theorem) located in secondary sources with their hypotheses (`std-theorems`), and a few still quoted from memory, among them that a KMS group of a faithful normal state is its modular group; (3) three identifications that are physics, not proved here: what the ball's graviton algebra is, that the vacuum's modular generator on each sector is the form used, and that the area-side charge is the perturbed generator with coupling 1. Still open: only modular time 1/2 is covered; the odd family's constants have no second code; the assembly's own closing argument has been checked by the lead and by no independent reader; the lower bounds bite only for l above about 15,000 times the coupling, which suffices for the divergence of the tail and is tested non-vacuously at a single l.
- **Which ruler.** On the ball's light cone, with the size rule the printed first-law proof uses, the area-side form exceeds the entropy on metric data and on extrinsic-curvature data alike at l = 2, in both parities, and electric-magnetic duality forbids it any asymmetry common to both parities. So the metric-against-extrinsic-curvature contrast above belongs to the fixed-volume, maximal-slice rule. Exact on four polynomial data per block; l = 3 and up not run. See `ruler`; it builds on a folder of the working record that is not published here. On the two static tides, `ruler-why` switches the fixed-volume rule's ingredients one at a time (l = 2, 3, 4): balance on metric data needs both the volume constraint and the extremal edge, while the excess on extrinsic-curvature data is set by the size rule alone, since no edge can couple to that datum. The maximal-slice ingredient was handled there by an argument, not computed.
- **Not claimed.** That an assumption of Jensen, Sorce and Speranza fails for gravitons, or that no equilibrium state exists for a ball. That the charge here is the one they define: `jss-id` examines the relation and leaves a gap. That any theorem of Satishchandran and Sorce is wrong: `ss-second` checks one statement of their graviton appendix in D = 4, 5 and 6 and withdraws two side statements made in `ss-tension`.

These reports cite sibling folders of the working record that are not published here (`hinge`, `hinge-e`, `gen-law` and some earlier waves). The scripts were run on Windows with Python 3.14; a few need `python-flint` in addition to `requirements.txt`, and say so in their report.

## Reproduce

```bash
pip install -r requirements.txt
./run_all.sh
```

About three minutes on a laptop. Exit code 0 means every check passed and every planted error was caught: 461 checks and 109 planted errors behind the theorem, the electric-sector and Killing-limit scripts, and the blind reader's own scripts. `FULL=1 ./run_all.sh` also reruns the slow numerical search for the extremal surface (about 30 minutes). The same command runs on every push (see the Actions tab).

All scripts use exact rational arithmetic (sympy) except four numerical checks and the lattice.

## What is in here

| Path | What it is |
|---|---|
| [`paper/`](paper/) | The short paper, source and PDF. |
| [`record/wave69/final/`](record/wave69/final/) | `THEOREM.md`: the full statement and proof, with six explicit hypotheses and a ledger of every premise (`out/LEDGER.md`). `PRINT_READ.md`: what was and was not read in the cited literature. |
| [`record/wave69/duality/`](record/wave69/duality/) | Haag duality for the helicity-2 field on a ball: proof and checks. |
| [`record/wave70/modflow/`](record/wave70/modflow/) | The modular group of the ball acts on curvature as the conformal flow: proof and checks. |
| [`record/wave69/bclose/`](record/wave69/bclose/) | The area side on exact second-order vacuum solutions, including a numerical search for the extremal surface. |
| [`record/wave64/btheorem/`](record/wave64/btheorem/) | What supports the choice of area side (Proposition R). |
| [`record/wave71/blindfinal/`](record/wave71/blindfinal/) | An adversarial read of the proofs by a checker instructed to break them: its brief and its verdict ("no step is broken", nine text fixes, and remarks the authors had missed). |
| [`record/wave72/decider/`](record/wave72/decider/) | The electric sector from the exact Weyl solution and the Killing-limit test, each with its pre-registration and a report that lists what went wrong as well as what passed. |
| [`record/wave72/blindread/`](record/wave72/blindread/) | A blind read of the short paper by a checker that saw only the paper: its findings, the scripts it wrote from the paper's text (`work/`), and how each finding was handled. |
| [`record/wave73/`](record/wave73/) | Work since version 3, one folder per question, each with a report that opens with the lead's check and grade: the area-side charge on every sector (`even-tower`, `charge-id`, `subspace`, `jss-id`), the odd-parity divergence and its lattice reproduction (`hinge-proof`, `hinge-blind`, `odd-sign`, `odd-remainder`), the even-parity measurement and proof steps (`even-blind`, `even-proof`, `even-limit`, `even-nonzero`, `even-passage`, `beta-def`, `even-remainder`, `even-second`, `even-density`, `even-dict`), the light-cone size rule at l = 2 and which ingredient of the fixed-volume rule matters (`ruler`, `ruler-why`), the standard theorems and the conditions of their use (`std-theorems`), the assembled statement with its hypothesis list (`assembly`), and a check against a published null-cone construction (`ss-tension`, `ss-second`). Not yet wired into `run_all.sh`. |
| [`record/wave55/gravlattice2/`](record/wave55/gravlattice2/) | A three-dimensional lattice computation of the graviton entropy with no multipole reduction. Stored outputs only: a rerun needs hours and several GB of memory. Its own verdict was INCONCLUSIVE for the question it was registered to answer; the paper quotes its rows for a different comparison and says so. |
| [`docs/MANIFESTO.md`](docs/MANIFESTO.md) | The standing rules of this repository. |

The folder names under `record/` are those of the working record the files were copied from, kept so that the scripts and the ledger run unchanged. Some documents cite other files of that record which are not published here yet (see below).

## How to break it

The result would be refuted by any of these, and we would like to hear about it:

1. An error in the proof that the modular flow of the helicity-2 field on a ball is the conformal flow (`record/wave70/modflow/REPORT.md`), or in Haag duality for the ball (`record/wave69/duality/REPORT.md`).
2. A second-order first law in print that balances on magnetic data with the same area prescription.
3. A demonstration that the stationary surface for magnetic data is not the one used.
4. A script here that passes when it should fail. Each one contains planted errors it must catch; if you find a wrong variant it does not catch, that is a finding.

Open an issue with three headings: **Why** (what you think is wrong), **What** (the step, equation or script), **Where** (file and line). Corrections are recorded, not erased: see `record/wave72/decider/REPORT.md` for examples in our own work.

## How this was made

The derivations, computations and checks were carried out by AI agents (Anthropic's Claude models) under the direction of the author, who set the questions, the scope and the rules. The record says which agent did what. The working rules were: outcome bins registered before each run; every script must catch a planted error; a result is credited only through content that was read and rerun; blind checkers are told to break the result; a second method counts for more than a second reader.

This is not a substitute for peer review, and one limit should be stated plainly: all checkers belong to one family of models and can share blind spots. That is why the scripts, the hypotheses and the failure conditions are published together.

## Not here yet

- The second proof of the entropy formula (a unitary transfer of each multipole sector onto scalar partial waves), cited in the paper as a draft.
- The earlier lattice code that `wave55` builds on.
- Results in other dimensions and for other spins.

## Citing

See [`CITATION.cff`](CITATION.cff).

## Licence

Code (everything under `record/**/scripts/` and `run_all.sh`): MIT, see [`LICENSE`](LICENSE). Text, papers and figures: CC BY 4.0, see [`LICENSE-docs.md`](LICENSE-docs.md).
