#!/usr/bin/env bash
# One command reruns every exact check behind the paper, and the quick deterministic checks of wave 73
# (what that section leaves out, and why, is in record/wave73/wave73.manifest).
#   ./run_all.sh          about eight minutes in a Linux container with eight CPUs; the wave 73 section is four of them
#   FULL=1 ./run_all.sh   also reruns the slow numerical search for the extremal surface (about 30 minutes)
# Exit code 0 means every check passed and every planted error was caught.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"

echo "### 1. The theorem: modular side, area side, the account, the ledger audit"
bash "$ROOT/record/wave69/final/scripts/run_all.sh"

echo "### 2. Electric sector from the exact Weyl solution"
( cd "$ROOT/record/wave72/decider/scripts" && python3 d1_weyl.py )

echo "### 3. Killing limit (Rindler wedge against the ball)"
( cd "$ROOT/record/wave72/decider/scripts" && python3 d2b_killing_fixed_control.py )

echo "### 4. Electric expansion in closed form (l = 1..6, and symbolic in l)"
( cd "$ROOT/record/wave72/decider/scripts" && python3 d3_weyl_closed_form.py )

echo "### 5. The blind reader's own scripts, written from the paper's text alone"
( cd "$ROOT/record/wave72/blindread/work" && for f in check1_magnetic check1b_green check2_weyl check3_electric_direct check4_weyl_flat_and_boost check5_magnetic_cut_unique check6_weyl_highl; do
    python3 "$f.py" > "../out/$f.txt"; echo "   $f ok"; done )

echo "### 6. Wave 73: the quick deterministic checks, and the stored lattice analyses rederived from their stored rows"
echo "###    (what is rerun and what is not: record/wave73/wave73.manifest)"
# FULL=1 also runs the checks the manifest marks slow (between fifteen and forty more minutes depending on the load)
SLOW_CHECKS="$([ "${FULL:-0}" = "1" ] && echo 1 || true)" bash "$ROOT/record/wave73/run_wave73.sh"

echo "### ALL CHECKS PASSED"
