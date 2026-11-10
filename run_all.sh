#!/usr/bin/env bash
# One command reruns every exact check behind the paper.
#   ./run_all.sh          about two minutes
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

echo "### ALL CHECKS PASSED"
