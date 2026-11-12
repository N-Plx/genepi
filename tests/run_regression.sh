#!/bin/bash
# Physics regression tests:
#  1) DVCS/BH cross-section at fixed kinematics vs tests/ref/xsec_points.ref
#  2) distributions of short fixed-seed samples vs tests/ref/distributions.json
# usage: tests/run_regression.sh      (from the top directory, after `make` and `make -C tests`)
set -u
TOP=$(cd "$(dirname "$0")/.." && pwd)
export GENEPI=$TOP
W=${REGRESSION_DIR:-$TOP/tests/out/regression}
rm -rf "$W"; mkdir -p "$W"
nfail=0

"$TOP/tests/xsec_points" "$TOP/gpd_table.dat" "$TOP/tests/ref/xsec_points.ref" || nfail=$((nfail+1))

# name | events | genepi options  (keep in sync with tests/make_references.sh)
while IFS='|' read -r name n opts; do
  [ -z "$name" ] && continue
  "$TOP/genepi" --trig $n --seed 20251 --output "$W/$name.lund" $opts > "$W/$name.log" 2>&1 \
    || { echo "FAIL  $name: genepi exit $?"; nfail=$((nfail+1)); continue; }
  python3 "$TOP/tests/check_distributions.py" "$name" "$W/$name.lund" "$TOP/tests/ref/distributions.json" || nfail=$((nfail+1))
done < "$TOP/tests/regression_samples.txt"

echo "regression tests: $nfail failing group(s)"
exit $((nfail > 0))
