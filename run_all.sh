#!/usr/bin/env bash
# One command reruns every exact check behind the paper, and the quick deterministic checks of wave 73
# (what that section leaves out, and why, is in record/wave73/wave73.manifest).
#   ./run_all.sh          about twelve minutes in a Linux container limited to four CPUs; the wave 73 section takes about four minutes of that
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
# FULL=1 also runs the checks the manifest marks slow (many more minutes: all 80 checks took 39 minutes in one loaded container run)
# The default run has a time budget: past it the section fails and names what was still running, instead of
# the CI job being cancelled at its 20-minute limit with nothing printed.
if [ "${FULL:-0}" = "1" ]; then
  SLOW_CHECKS=1 SECTION_BUDGET= bash "$ROOT/record/wave73/run_wave73.sh"
else
  SECTION_BUDGET="${SECTION_BUDGET:-600}" bash "$ROOT/record/wave73/run_wave73.sh"
fi

echo "### ALL CHECKS PASSED"
