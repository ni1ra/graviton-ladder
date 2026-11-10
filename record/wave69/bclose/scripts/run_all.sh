#!/usr/bin/env bash
# w69-bclose: one command.  Plants first, then the unplanted runs.  Exact parts ~10 s; numerics tens of minutes.
set -euo pipefail
cd "$(dirname "$0")"
python3 p1_enum.py --plants
python3 p1_enum.py
python3 p2_solve.py
python3 p2_geo.py --plants
python3 p2_geo.py
