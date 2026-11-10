#!/usr/bin/env bash
# w69-final: one command reruns everything the ladder theorem's ledger points at.
#   ./run_all.sh          exact scripts of this lane, exact scripts of the record (rerun in a scratch copy), ledger audit
#   FULL=1 ./run_all.sh   also reruns the slow numerics of the record in the scratch copy (about 30 minutes)
# Nothing outside wave69/final is written: record scripts are copied to a temporary directory and run there.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
FINAL="$(dirname "$HERE")"
ROOT="$(cd "$FINAL/../.." && pwd)"
OUT="$FINAL/out"; UP="$OUT/upstream"
mkdir -p "$OUT" "$UP"
cd "$HERE"
echo "== this lane (exact; planted errors inside each script)"
python3 s1_modular_side.py
python3 s2_area_side.py
python3 s3_account.py
python3 s4_kms_toy.py

echo "== record scripts, rerun in a scratch copy"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
for d in wave69/duality wave70/modflow wave69/bclose wave64/btheorem wave45/truthcond; do
  mkdir -p "$TMP/$d/out"; cp -r "$ROOT/$d/scripts" "$TMP/$d/"
done
( cd "$TMP"
  python3 wave69/duality/scripts/s1_identities.py wave69/duality/out/s1_out.json > /dev/null
  python3 wave70/modflow/scripts/s1_invariance.py > /dev/null
  python3 wave70/modflow/scripts/s5_wedge.py > /dev/null
  python3 wave69/bclose/scripts/p1_enum.py --plants > /dev/null
  python3 wave69/bclose/scripts/p1_enum.py > /dev/null
  python3 wave69/bclose/scripts/p2_solve.py > /dev/null
  python3 wave64/btheorem/scripts/s1s2_green_gauge.py wave64/btheorem/out/s1s2.json > /dev/null
  python3 wave45/truthcond/scripts/check_E2.py > wave45/truthcond/out/check_E2.json
  if [ "${FULL:-0}" = "1" ]; then
    python3 wave69/bclose/scripts/p2_geo.py --plants > /dev/null
    python3 wave69/bclose/scripts/p2_geo.py > /dev/null
    cp wave69/bclose/out/p2_geo.json "$UP/bclose_p2_geo_rerun.json"
  fi )
cp "$TMP/wave69/duality/out/s1_out.json"      "$UP/duality_s1.json"
cp "$TMP/wave70/modflow/out/s1_out.json"      "$UP/modflow_s1.json"
cp "$TMP/wave70/modflow/out/s5_out.json"      "$UP/modflow_s5.json"
cp "$TMP/wave69/bclose/out/p1_plants.json"    "$UP/bclose_p1_plants.json"
cp "$TMP/wave69/bclose/out/p1_enum.json"      "$UP/bclose_p1_enum.json"
cp "$TMP/wave69/bclose/out/p2_solve.json"     "$UP/bclose_p2_solve.json"
cp "$TMP/wave64/btheorem/out/s1s2.json"       "$UP/btheorem_s1s2.json"
cp "$TMP/wave45/truthcond/out/check_E2.json"  "$UP/truthcond_check_E2.json"
echo "   record reruns done"

echo "== ledger audit"
python3 s5_ledger.py
( cd "$OUT" && sha256sum s1_modular_side.json s2_area_side.json s3_account.json s4_kms_toy.json s5_ledger.json upstream/*.json > SHA256SUMS )
echo "ALL DONE: outputs in $OUT"
