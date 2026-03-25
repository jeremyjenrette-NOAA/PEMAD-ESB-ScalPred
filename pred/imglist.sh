#!/usr/bin/env bash
set -euo pipefail

BASE="/z/data/habcam/proc/Images"
OUTDIR="./img_inventory_out"
LOG="$OUTDIR/img_inventory.log"

mkdir -p "$OUTDIR"

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <year>"
  echo "Example: $0 2022"
  exit 1
fi

yr="$1"
year_dir="$BASE/$yr"
out_txt="$OUTDIR/img_inventory_${yr}.txt"
out_count="$OUTDIR/img_inventory_${yr}_count.txt"

log() {
  echo "$(date '+%F %T') | $*" | tee -a "$LOG"
}

log "Starting image inventory build for $yr"

if [[ -f "$out_txt" ]]; then
  log "Checkpoint exists for $yr, skipping"
  exit 0
fi

if [[ ! -d "$year_dir" ]]; then
  log "Missing directory: $year_dir"
  exit 1
fi

log "Scanning $year_dir"

find -L "$year_dir" -type f -iname '*.png' 2>/dev/null | \
  grep -v $'\r' | \
  sort > "$out_txt"

wc -l < "$out_txt" > "$out_count"
n=$(cat "$out_count")

log "Finished $yr | files found: $n"
log "Saved file list to $out_txt"
log "Saved count to $out_count"
log "Inventory complete for $yr"
