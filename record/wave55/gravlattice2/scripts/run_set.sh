#!/bin/bash
# usage: run_set.sh centre R1 R2 ... ; rule-B rows via the sector solver; G cache path in $GCACHE (optional, derivable by G_array)
cd "$(dirname "$0")"
c=$1; shift
for R in "$@"; do
  OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 /usr/bin/time -v python3 run2.py $R $c ../out/rows $GCACHE > ../out/rows/log_R${R}_${c}_B.txt 2>&1
done
