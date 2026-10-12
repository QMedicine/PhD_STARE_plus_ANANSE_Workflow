#!/usr/bin/env bash
set -euo pipefail

# This file is executed *inside* the container.
# It overrides STARE.sh via a shell function defined in BASH_ENV.

# Resolve the real STARE.sh (recorded by beforeScript); fall back to PATH if needed
REAL="${PWD}/.shim/.real_path"
if [[ -s "$REAL" ]]; then
  REAL_BIN="$(cat "$REAL")"
else
  REAL_BIN="$(command -v STARE.sh || true)"
fi
if [[ -z "${REAL_BIN:-}" || ! -x "$REAL_BIN" ]]; then
  echo "[wrapper] ERROR: cannot resolve real STARE.sh" >&2
  exit 127
fi

# Parse args to capture -o (output/sample id) and -b (peaks)
outprefix=""
bedfile=""
declare -a passthru
while (( "$#" )); do
  case "$1" in
    -o|--output)
      outprefix="$2"; passthru+=("$1" "$2"); shift 2;;
    -b|--bed_file)
      bedfile="$2"; shift 2;;
    *)
      passthru+=("$1"); shift;;
  esac
done

# If a BED was provided, trim to 3 columns (STARE’s preproc expects BED3)
if [[ -n "$bedfile" ]]; then
  bed3=".shim/$(basename "$bedfile").bed3"
  awk 'BEGIN{OFS="\t"} !/^#/ {print $1,$2,$3}' "$bedfile" > "$bed3"
  passthru+=("-b" "$bed3")
fi

# Look up an ABC file for this sample
abc=""
if [[ -n "$outprefix" && -s ".shim/abc_map.csv" ]]; then
  # CSV format: sample,abc_path
  abc="$(awk -F',' -v s="$outprefix" 'NR>1 && $1==s {print $2; exit}' .shim/abc_map.csv || true)"
fi

if [[ -n "$abc" && -s "$abc" ]]; then
  echo "[wrapper] Using ABC for sample '${outprefix}': $abc" >&2
  exec "$REAL_BIN" "${passthru[@]}" -r "$abc"
else
  echo "[wrapper] No ABC mapping for '${outprefix}' — running STARE normally." >&2
  exec "$REAL_BIN" "${passthru[@]}"
fi
