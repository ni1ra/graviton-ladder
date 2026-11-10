#!/usr/bin/env bash
# One command reruns every exact check behind the paper.
#   ./run_all.sh          about three minutes
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

echo "### ALL CHECKS PASSED"
