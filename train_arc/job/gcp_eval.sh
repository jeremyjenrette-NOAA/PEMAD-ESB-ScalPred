#!/bin/bash
# ==============================================================================
# Re-run Stage 1 eval + two-stage cascade eval on an EXISTING run (no retraining).
#
#   ./job/gcp_eval.sh <system> <run_name>
#   e.g. ./job/gcp_eval.sh crab 2426crabdata_yolo12n_gcp_20260831_143915
#
# <run_name> is a folder under output/. Uses the same job/systems/<system>.env
# profile as gcp_yolo_job.sh. Run from train_arc/.
# ==============================================================================
set -euo pipefail

SYSTEM="${1:?usage: $0 <system> <run_name>}"
RUN_NAME="${2:?usage: $0 <system> <run_name>}"
PROFILE="$(dirname "$0")/systems/${SYSTEM}.env"
[ -f "${PROFILE}" ] || { echo "ERROR: no system profile ${PROFILE}"; exit 1; }
# shellcheck disable=SC1090
source "${PROFILE}"

export PYTORCH_ALLOC_CONF=expandable_segments:True
export CUDA_VISIBLE_DEVICES=0

DATA_DIR="${LABEL}${YEAR}"
ORIGINAL_YOLO_ROOT="${DATA_DIR}/yolo"
BROAD_YOLO_ROOT="${DATA_DIR}/yolo_broad"
TAXONOMY_JSON="config/${LABEL}_taxonomy.json"
RUN_DIR="$(pwd)/output/${RUN_NAME}"
MODEL_TAG=$(basename "$MODEL" .pt | tr '.-' '__')
BEST_WEIGHTS="${RUN_DIR}/weights/best.pt"
STAGE2_WEIGHTS="${RUN_DIR}/weights/${LABEL}_tax.pt"

for f in "${BEST_WEIGHTS}" "${STAGE2_WEIGHTS}" "${TAXONOMY_JSON}"; do
    [ -f "$f" ] || { echo "ERROR: missing $f"; exit 1; }
done

python -c "
from config.taxonomy_utils import load_taxonomy, validate_taxonomy_against_yolo
t = load_taxonomy('${TAXONOMY_JSON}')
validate_taxonomy_against_yolo(t, '${ORIGINAL_YOLO_ROOT}/data.yaml', taxonomy_json_path='${TAXONOMY_JSON}')
print('Taxonomy OK')
"

source config/write_val.sh
write_val_image_csv "${BROAD_YOLO_ROOT}" "${RUN_DIR}/eval/val_images${YEAR}_${MODEL_TAG}.csv"

python config/eval_yolo_detections_hung_multi.py \
    --year "${YEAR}" --model_name "${MODEL_TAG}" --run_dir "${RUN_DIR}" \
    --weights "${BEST_WEIGHTS}" \
    --data_root "${BROAD_YOLO_ROOT}" --gt_data_root "${ORIGINAL_YOLO_ROOT}" \
    --taxonomy_json "${TAXONOMY_JSON}" \
    --out_csv "${RUN_DIR}/eval/autotest.csv" \
    --gt_out_csv "${RUN_DIR}/eval/mantest.csv" \
    --out_fn_csv "${RUN_DIR}/eval/fn.csv" \
    --imgsz 1024 --conf 0.01 --nms_iou 0.65 --match_iou 0.1 \
    --max_det "${MAX_DET}" --device 0 --batch 4 --spname ${BROAD_CLASS_NAME}

python config/eval_two_stage_predictions.py \
    --autotest_csv "${RUN_DIR}/eval/autotest.csv" \
    --val_img_dir "${ORIGINAL_YOLO_ROOT}/images/val" \
    --stage2_weights "${STAGE2_WEIGHTS}" \
    --taxonomy_json "${TAXONOMY_JSON}" \
    --out_csv "${RUN_DIR}/eval/autotest_two_stage_cascade.csv"

echo "=== Eval Completed Successfully (${SYSTEM}) ==="
