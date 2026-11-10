# Manifesto

Standing rules of this repository. Temporary state belongs in issues, not here.

## What this repository is

One result and everything needed to check it: the exact comparison of graviton relative entropy with the second-order area term for a flat ball. Anything that can stand alone as another result gets its own repository or paper.

## Rules for claims

1. Every claim carries its scope. "Proved" means proved under the hypotheses listed in `record/wave69/final/THEOREM.md`, section 2. A sentence that drops the scope is a bug.
2. Every number traces to a script and a stored output. A number that exists only in prose does not exist.
3. Every script contains planted errors that it must catch. A check that cannot fail is removed.
4. A test's possible outcomes and their meaning are written down and committed before the test runs. The file stays in the repository whatever the outcome.
5. Corrections are added, never substituted. A wrong statement, a failed control or a slipped procedure is recorded where it happened, with what was done about it.
6. The record states what was not read and what was not checked.
7. Nobody's opinion moves a claim, including ours and including a referee's. An argument or a computation does.

## Rules for changes

- `./run_all.sh` passes on every commit to `main`. Continuous integration runs it on every push and pull request.
- Files under `record/` are historical. They change only to fix a rerun on a new machine or to remove private wording, and the commit says which. New work goes in a new folder.
- Issues and pull requests use three headings: **Why**, **What**, **Where**.
- A pull request that weakens a check must say so in its title.

## Architecture

- `paper/`: the short paper.
- `record/<wave>/<name>/`: one piece of work each, with its pre-registration, report, `scripts/` and `out/`. Folder names are those of the working record, kept so scripts and the ledger run unchanged.
- `run_all.sh`: the single entry point. It calls `record/wave69/final/scripts/run_all.sh` (the theorem's bundle and ledger audit) and the `wave72` scripts.
- `record/wave69/final/scripts/s5_ledger.py`: the audit. It holds every premise of the theorem as a typed graph and fails on a cycle, a dangling dependency, a missing or failing output, or a script without a caught plant.
