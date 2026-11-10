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

**Start here:** [`paper/ladder_short.pdf`](paper/ladder_short.pdf) (8 pages, version 2).

## Status

Read this before quoting the result.

- **Proved here, with scripts:** the entropy formula (two routes), the area side on magnetic data, the area side on electric data from the exact Weyl solution, the wedge identity.
- **Assumed, and stated as assumptions:** the theory is free linearised gravity at order G⁰ on a flat ball in four dimensions; the gravitational term of the first law at this order is the stationary value of Area − 2·Volume; three conditions on what counts as "the algebra of the ball"; Einstein's equations to second order. The full list is in the paper, section 8, and in [`record/wave69/final/THEOREM.md`](record/wave69/final/THEOREM.md), section 2.
- **Open:** whether the ball algebra of full quantum gravity reduces to the free one as G → 0. Nothing here tests that.
- **Not peer reviewed.** Every check so far was made inside this programme. No outside referee has read the proofs.
- **Corrected once already.** Version 1 of the paper was public for about an hour before a blind read came back. It confirmed every formula and found four statements that were wrong or too broad, including that part of the area result is in print. Version 2 fixes them and lists them in its last section. The read is in [`record/wave72/blindread/`](record/wave72/blindread/REPORT.md).

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
