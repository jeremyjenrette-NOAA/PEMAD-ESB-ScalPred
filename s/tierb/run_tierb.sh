#!/bin/bash
# Tier B (independent 2026 survey test): split -> two-stage calibration.
# Run from the repo root or anywhere:  bash s/tierb/run_tierb.sh [--mode whole]
# First copy from ARC (train_arc/output/2224scallop_yolo12n_slurm_598580/eval_tierb/tierb2026/)
# into data/raw_model/2224scallop_yolo12n_slurm_598580/tierb2026/:
#   autotest_tierb2026_yolo12n.csv  tierb_manifest.csv  (others optional)
set -eo pipefail
cd "$(dirname "$0")/.." || exit 1          # -> s/
SRC=../data/raw_model/2224scallop_yolo12n_slurm_598580/tierb2026
for f in tierb_manifest.csv autotest_tierb2026_yolo12n.csv; do
  [[ -f "$SRC/$f" ]] || { echo "Missing $SRC/$f (copy it from ARC first)"; exit 1; }
done
mkdir -p logs
STAMP=$(date +%Y%m%d_%H%M%S)
caffeinate -w $$ 2>/dev/null &
python3 tierb/split_tierb2026.py --manifest "$SRC/tierb_manifest.csv" \
  --out ../data/raw/dataset_split_tierb2026.csv "$@" | tee "logs/${STAMP}_tierb_split.log"
Rscript tierb/tierb_calibration.R 2>&1 | tee "logs/${STAMP}_tierb_calibration.log"
