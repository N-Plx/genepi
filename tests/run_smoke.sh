#!/bin/bash
# Smoke tests: short fixed-seed runs for every target/process, LUND consistency checks,
# command-line validation and exit codes, seed reproducibility.
# usage: tests/run_smoke.sh            (from the repository top directory, after `make`)
set -u
TOP=$(cd "$(dirname "$0")/.." && pwd)
GENEPI_EXE=$TOP/genepi
export GENEPI=$TOP
W=${SMOKE_DIR:-$TOP/tests/out/smoke}
rm -rf "$W"; mkdir -p "$W"; cd "$W"
nfail=0
ok()   { echo "PASS  $*"; }
bad()  { echo "FAIL  $*"; nfail=$((nfail+1)); }

# expect NAME EXPECTED_RC -- genepi args...
expect() {
  local name=$1 want=$2; shift 3
  "$GENEPI_EXE" "$@" > "$name.log" 2>&1; local rc=$?
  if [ $rc -eq $want ]; then ok "$name (exit $rc)"; else bad "$name: exit $rc, expected $want"; tail -5 "$name.log"; fi
}

# ---- 1) physics runs: every target and process, checked with check_lund.py
run_and_check() {
  local name=$1 target=$2; shift 2
  local extra=()
  "$GENEPI_EXE" --trig 50 --seed 1234 --output "$name.lund" "$@" > "$name.log" 2>&1; local rc=$?
  if [ $rc -ne 0 ]; then bad "$name: genepi exit $rc"; tail -5 "$name.log"; return; fi
  local n; n=$(awk 'NF==10' "$name.lund" | wc -l | tr -d ' ')
  [ "$n" -eq 50 ] || { bad "$name: $n events instead of 50"; return; }
  # pass raster settings to the checker when given
  local args=("$@"); local rx=0 ry=0
  for ((k=0; k<${#args[@]}; k++)); do
    [ "${args[$k]}" = "--raster_x" ] && rx=${args[$((k+1))]}
    [ "${args[$k]}" = "--raster_y" ] && ry=${args[$((k+1))]}
  done
  if python3 "$TOP/tests/check_lund.py" "$name.lund" --target "$target" --raster "$rx" "$ry" > "$name.check" 2>&1
  then ok "$name ($(head -1 "$name.check" | sed 's/.*: //'))"; else bad "$name"; cat "$name.check"; fi
}
run_and_check p_dvcs     p   --targ_A 1  --targ_Z 1
run_and_check n_dvcs     n   --targ_A 1  --targ_Z 0
run_and_check d_dvcs     d   --targ_A 2  --targ_Z 1
run_and_check he4_dvcs   he4 --targ_A 4  --targ_Z 4
run_and_check n14_dvcs   n14 --targ_A 14 --targ_Z 7
run_and_check nh3_dvcs   nh3 --targ_A 17 --targ_Z 10
run_and_check nd3_dvcs   nd3 --targ_A 20 --targ_Z 10
run_and_check p_pi0      p   --targ_A 1  --targ_Z 1  --process 1 --meson 0
run_and_check d_pi0      d   --targ_A 2  --targ_Z 1  --process 1 --meson 0
run_and_check nh3_pi0    nh3 --targ_A 17 --targ_Z 10 --process 1 --meson 0
run_and_check p_eta      p   --targ_A 1  --targ_Z 1  --process 1 --meson 1
run_and_check p_raster   p   --targ_A 1  --targ_Z 1  --raster_x 0.6 --raster_y 0.4
# known issue, not tested yet: free-neutron meson output (pi0 written with only one decay photon)

# ---- 2) reproducibility
"$GENEPI_EXE" --trig 20 --seed 77 --process 1 --output r1.lund > r1.log 2>&1
"$GENEPI_EXE" --trig 20 --seed 77 --process 1 --output r2.lund > r2.log 2>&1
cmp -s r1.lund r2.lund && ok "same --seed gives identical output" || bad "same --seed gives different output"
"$GENEPI_EXE" --trig 20 --output a1.lund > a1.log 2>&1 & "$GENEPI_EXE" --trig 20 --output a2.lund > a2.log 2>&1 & wait
cmp -s a1.lund a2.lund && bad "two runs without --seed are identical" || ok "two simultaneous runs without --seed differ"

# ---- 3) command-line validation and exit codes
expect bad_option      1 -- --trig 5 --Ebem 10.6
expect missing_value   1 -- --trig 5 --Ebeam
expect garbage_value   1 -- --trig 5 --Ebeam abc
expect stray_argument  1 -- --trig 5 junk
expect inverted_range  1 -- --trig 5 --x_min 0.5 --x_max 0.1
expect bad_seed        1 -- --trig 5 --seed 2147483648
expect bad_target      1 -- --trig 5 --targ_A 1 --targ_Z 2
expect bad_photon_cut  1 -- --trig 5 --gamma_beam_min 200
expect empty_phase     2 -- --trig 5 --seed 1 --Ebeam 2 --max_trials 20000
expect max_trials      3 -- --trig 1000 --seed 1 --max_trials 5000
GENEPI=/nonexistent "$GENEPI_EXE" --trig 5 > nogpd.log 2>&1
[ $? -eq 1 ] && ok "missing GPD table (exit 1)" || bad "missing GPD table: wrong exit code"

echo "smoke tests: $nfail failure(s)"
exit $((nfail > 0))
