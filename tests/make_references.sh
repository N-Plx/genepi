#!/bin/bash
# (Re)create the regression references after an INTENDED physics change, then commit them
# together with the change:
#   tests/ref/xsec_points.ref    cross-section at fixed kinematics
#   tests/ref/distributions.json statistics of large samples (REF_EVENTS per sample, other seed)
# usage: tests/make_references.sh [REF_EVENTS]      (default 20000; takes several minutes)
set -eu
TOP=$(cd "$(dirname "$0")/.." && pwd)
export GENEPI=$TOP
N=${1:-20000}
W=$TOP/tests/out/references; rm -rf "$W"; mkdir -p "$W" "$TOP/tests/ref"
"$TOP/tests/xsec_points" "$TOP/gpd_table.dat" "$TOP/tests/ref/xsec_points.ref" --write
rm -f "$TOP/tests/ref/distributions.json"
while IFS='|' read -r name n opts; do
  [ -z "$name" ] && continue
  "$TOP/genepi" --trig $N --seed 99991 --output "$W/$name.lund" $opts > "$W/$name.log" 2>&1 &
done < "$TOP/tests/regression_samples.txt"
wait
while IFS='|' read -r name n opts; do
  [ -z "$name" ] && continue
  python3 "$TOP/tests/check_distributions.py" "$name" "$W/$name.lund" "$TOP/tests/ref/distributions.json" --write
done < "$TOP/tests/regression_samples.txt"
