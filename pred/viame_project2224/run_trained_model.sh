#!/bin/bash

# -----------------------------
# User-configurable settings
# -----------------------------

# Root VIAME install
export VIAME_INSTALL="/projects/sharkpulse/archived/viame" # wherever your VIAME is installed
# Project folder = directory this script lives in
export VIAME_PROJECT_DIR="/projects/sharkpulse/archived/PEMAD-ESB-ScalPred/pred/viame_project2224" # wherever you are running this script
# Core processing options
export INPUT_LIST="input_list_test.txt"
JOB_TAG="${SLURM_JOB_ID:-manual}"
# export OUTPUT_DIRECTORY="pred_out_2024_${JOB_TAG}"
export OUTPUT_DIRECTORY="pred_out"
export INPUT_FRAME_RATE=1
export PROCESS_FRAME_RATE=1

# GPU / parallel settings
export TOTAL_GPU_COUNT=1
export PIPES_PER_GPU=1

# Chunk / resume options
export CHUNK_SIZE=5000
export RESUME_FILE="${OUTPUT_DIRECTORY}/completed.txt"
export PROGRESS_FILE="${OUTPUT_DIRECTORY}/progress.log"
export MASTER_CSV="${OUTPUT_DIRECTORY}/all_detections.csv"

# Stereo split behavior
export IMAGE_SPLIT="right"

# -----------------------------
# Environment setup
# -----------------------------
module reset
source ~/.bashrc
conda deactivate

source "${VIAME_INSTALL}/setup_viame.sh"

mkdir -p "${OUTPUT_DIRECTORY}"

echo "======================================"
echo "Starting VIAME inference"
echo "VIAME_INSTALL      = ${VIAME_INSTALL}"
echo "VIAME_PROJECT_DIR  = ${VIAME_PROJECT_DIR}"
echo "INPUT_LIST         = ${INPUT_LIST}"
echo "OUTPUT_DIRECTORY   = ${OUTPUT_DIRECTORY}"
echo "TOTAL_GPU_COUNT    = ${TOTAL_GPU_COUNT}"
echo "PIPES_PER_GPU      = ${PIPES_PER_GPU}"
echo "CHUNK_SIZE         = ${CHUNK_SIZE}"
echo "IMAGE_SPLIT        = ${IMAGE_SPLIT}"
echo "======================================"

cd "${VIAME_PROJECT_DIR}"

python "${VIAME_PROJECT_DIR}/predict_habcam_viame.py" \
  -l "${INPUT_LIST}" \
  -ifrate "${INPUT_FRAME_RATE}" \
  -frate "${PROCESS_FRAME_RATE}" \
  -p "${VIAME_PROJECT_DIR}/pipelines/detector_project_folder.pipe" \
  -o "${OUTPUT_DIRECTORY}" \
  --no-reset-prompt \
  -gpus "${TOTAL_GPU_COUNT}" \
  -pipes-per-gpu "${PIPES_PER_GPU}" \
  --chunk-size "${CHUNK_SIZE}" \
  --resume-file "${RESUME_FILE}" \
  --progress-file "${PROGRESS_FILE}" \
  --master-csv "${MASTER_CSV}" \
  --split "${IMAGE_SPLIT}"

echo "Inference complete."