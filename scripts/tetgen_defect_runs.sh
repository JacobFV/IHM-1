#!/usr/bin/env bash
# One TetGen run, logged under a name the report parser understands.
#
#   scripts/tetgen_defect_runs.sh <variant> <label> <flags> <condition> <i> [extra env...]
#
# variant   : suffix of data/runtime/tetgen-asan/build/tetgen_driver_<variant>
# label     : PLC label under data/derived/tetgen-defect-v1/plc/<label>.plcbin
# flags     : TetGen switch string, e.g. pY
# condition : free text recorded in the log name, e.g. base, norandom, perturb
#
# TG_PREFIX wraps the command (e.g. TG_PREFIX="setarch -R").
# TG_TIMEOUT caps the run in seconds (default 7200); a timeout is recorded as exit 124.
#
# Runs in a scratch CWD so TetGen's tetgen-tmpfile_skipped.* droppings never land in the
# repository, and cleans them up afterwards.
set -uo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="$ROOT/data/derived/tetgen-defect-v1"
variant="$1"; label="$2"; flags="$3"; cond="$4"; i="$5"; shift 5
mkdir -p "$OUT/logs"
scratch="$(mktemp -d)"
log="$OUT/logs/${variant}__${label}__${flags}__${cond}__${i}.log"
plc="$OUT/plc/${label}.plcbin"
[ -f "$plc" ] || { echo "no plc $plc" >&2; exit 2; }
if [ "$variant" = wheel ]; then
  # the shipped libigl wheel, driven through the probe script so it prints the same
  # RESULT line the native driver prints
  cmd=("$ROOT/data/runtime/geometry/libigl-2.6.2/venv/bin/python"
       "$ROOT/scripts/probe_tetgen_boundary_recovery_defect.py"
       --wheel "$label" --flags "$flags" --out "$OUT")
else
  drv="$ROOT/data/runtime/tetgen-asan/build/tetgen_driver_${variant}"
  [ -x "$drv" ] || { echo "no driver $drv" >&2; exit 2; }
  cmd=("$drv" "$plc" "$flags" 1)
fi
( cd "$scratch" && env "$@" ${TG_PREFIX:-} nice -n 15 timeout "${TG_TIMEOUT:-7200}" "${cmd[@]}" ) > "$log" 2>&1
rc=$?
echo "# exit_code $rc" >> "$log"
rm -rf "$scratch"
grep -hE '^RESULT|corrupted|not recovered' "$log" | head -4
echo "-> $log (exit $rc)"
